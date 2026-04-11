
! =========================================================================
! curve_fill_solver — Euler solver on a Hilbert AMR mesh
!
! State arrays are cell-centred, indexed 1..mesh%length.
! Solid cells (is_solid==1) carry no physical state and are skipped in the
! flux loop; ghost cells have their state overwritten after each flux step.
! =========================================================================

module cf_solver_module
  use types
  use euler_flux
  use iso_c_binding
  implicit none

  interface
    subroutine emit_cf_state(s) bind(C, name="emit_cf_state")
      use iso_c_binding
      use types
      implicit none
      type(cf_state_c), intent(in) :: s
    end subroutine emit_cf_state
  end interface

  ! -----------------------------------------------------------------------
  ! Flow state — all arrays length = mesh%length
  ! -----------------------------------------------------------------------
  type :: cf_state
    integer :: ncells
    ! conserved variables
    real(8), allocatable :: ro(:), rovx(:), rovy(:), roe(:)
    ! increments accumulated during one flux pass
    real(8), allocatable :: dro(:), drovx(:), drovy(:), droe(:)
    ! secondary (derived each step)
    real(8), allocatable :: p(:), vx(:), vy(:), hstag(:)
    ! per-cell geometry (computed once from mesh)
    real(8), allocatable :: area(:)
    real(8), allocatable :: dt(:)
  end type cf_state

  contains

  ! -----------------------------------------------------------------------
  ! Allocate and initialise state arrays.
  ! Freestream: uniform flow at Mach Ma, angle alpha_deg, stagnation
  ! conditions pstag/tstag (same as block solver bconds convention).
  ! -----------------------------------------------------------------------
  subroutine init_state(mesh, av, bcs, st)
    type(lod_mesh),    intent(in)  :: mesh
    type(t_appvars),   intent(in)  :: av
    type(t_bconds),    intent(in)  :: bcs
    type(cf_state),    intent(out) :: st

    integer :: i
    real(8) :: gam, rgas, pstag, tstag, alpha
    real(8) :: ro0, p0, T0, c0, Ma, u0, v0, roe0

    st%ncells = mesh%length
    allocate(st%ro   (st%ncells), st%rovx (st%ncells), &
             st%rovy (st%ncells), st%roe  (st%ncells))
    allocate(st%dro  (st%ncells), st%drovx(st%ncells), &
             st%drovy(st%ncells), st%droe (st%ncells))
    allocate(st%p    (st%ncells), st%vx   (st%ncells), &
             st%vy   (st%ncells), st%hstag(st%ncells))
    allocate(st%area (st%ncells))
    allocate(st%dt   (st%ncells))

    ! --- cell areas (trivial for axis-aligned cells) ----------------------
    do i = 1, st%ncells
      st%area(i) = (mesh%cells(i)%xmax - mesh%cells(i)%xmin) * &
                   (mesh%cells(i)%ymax - mesh%cells(i)%ymin)
    end do

    ! --- freestream from stagnation conditions ----------------------------
    gam   = real(av%gam,  8)
    rgas  = real(av%rgas, 8)
    pstag = real(bcs%pstag, 8)
    tstag = real(bcs%tstag, 8)
    alpha = real(bcs%alpha, 8) * acos(-1.0D0) / 180.0D0   ! deg → rad

    ! Derive freestream Mach from p_out/pstag isentropically,
    ! matching the block solver convention (rfin is a relaxation factor, not Ma)
    p0  = real(bcs%p_out, 8)
    Ma  = sqrt(2.0D0/(gam-1.0D0) * ((pstag/p0)**((gam-1.0D0)/gam) - 1.0D0))
    T0  = tstag / (1.0D0 + 0.5D0*(gam-1.0D0)*Ma**2)
    ro0 = p0 / (rgas * T0)
    c0  = sqrt(gam * p0 / ro0)
    u0  = Ma * c0 * cos(alpha)
    v0  = Ma * c0 * sin(alpha)
    roe0 = p0/(gam-1.0D0) + 0.5D0*ro0*(u0**2 + v0**2)

    print '(A,F6.4,A,F8.1,A,F8.1)', ' CF freestream: Ma=', Ma, '  p0=', p0, '  T0=', T0

    st%ro   = ro0
    st%rovx = ro0 * u0
    st%rovy = ro0 * v0
    st%roe  = roe0
    ! zero solid cells (they carry no physical meaning)
    do i = 1, st%ncells
      if (mesh%cells(i)%is_solid == 1) then
        st%ro(i) = 0.0D0;  st%rovx(i) = 0.0D0
        st%rovy(i) = 0.0D0; st%roe(i) = 0.0D0
      end if
    end do

  end subroutine init_state

  ! -----------------------------------------------------------------------
  ! Secondary variables from conserved.  Skips solid cells.
  ! -----------------------------------------------------------------------
  subroutine set_secondary(mesh, av, st)
    type(lod_mesh),  intent(in)    :: mesh
    type(t_appvars), intent(in)    :: av
    type(cf_state),  intent(inout) :: st

    integer :: i
    real(8) :: gam, pi

    real(8), parameter :: RO_MIN = 1.0D-6   ! positivity floor
    real(8), parameter :: P_MIN  = 1.0D0    ! 1 Pa floor

    gam = real(av%gam, 8)
    do i = 1, st%ncells
      if (mesh%cells(i)%is_solid == 1) cycle
      ! Positivity: clamp density and recompute if needed
      if (st%ro(i) < RO_MIN) st%ro(i) = RO_MIN
      st%p(i) = pressure(st%ro(i), st%rovx(i), st%rovy(i), st%roe(i), gam)
      if (st%p(i) < P_MIN) then
        st%p(i)   = P_MIN
        st%roe(i) = P_MIN/(gam-1.0D0) + 0.5D0*(st%rovx(i)**2 + st%rovy(i)**2)/st%ro(i)
      end if
      st%vx(i)    = st%rovx(i) / st%ro(i)
      st%vy(i)    = st%rovy(i) / st%ro(i)
      st%hstag(i) = (st%roe(i) + st%p(i)) / st%ro(i)
    end do
  end subroutine set_secondary

  ! -----------------------------------------------------------------------
  ! Ghost cell boundary condition.
  ! For each ghost cell, interpolate the flow state at the mirror point
  ! (bilinear from the four nearest cell centres), then reflect the normal
  ! velocity to enforce the slip wall condition.
  !
  ! Simple fallback when the mirror point hits a solid or boundary region:
  ! use the ghost cell's own current state with just the normal velocity
  ! negated (zeroth-order ghost).
  ! -----------------------------------------------------------------------
  subroutine apply_ghost_bc(mesh, st)
    type(lod_mesh), intent(in)    :: mesh
    type(cf_state), intent(inout) :: st

    integer  :: g, ci, j, best
    real(8)  :: nx, ny, mx, my, d2, best_d2
    real(8)  :: cx, cy, un
    real(8)  :: ro_m, rovx_m, rovy_m, roe_m

    do g = 1, mesh%ghost_count
      ci = mesh%ghost_indices(g)
      nx = mesh%ghost_normals(g,1)
      ny = mesh%ghost_normals(g,2)
      mx = mesh%ghost_mirror (g,1)
      my = mesh%ghost_mirror (g,2)

      ! Find the fluid cell whose centre is closest to the mirror point.
      ! (Later: proper bilinear interpolation using neigh_indices.)
      best_d2 = 1.0D300
      best    = ci   ! fallback: self
      do j = 1, mesh%length
        if (mesh%cells(j)%is_solid == 1) cycle
        if (j == ci) cycle
        cx = 0.5D0*(mesh%cells(j)%xmin + mesh%cells(j)%xmax)
        cy = 0.5D0*(mesh%cells(j)%ymin + mesh%cells(j)%ymax)
        d2 = (cx-mx)**2 + (cy-my)**2
        if (d2 < best_d2) then
          best_d2 = d2
          best    = j
        end if
      end do

      ro_m   = st%ro  (best)
      rovx_m = st%rovx(best)
      rovy_m = st%rovy(best)
      roe_m  = st%roe (best)

      ! Reflect normal velocity component: u_ghost = u_mirror - 2*(u·n)*n
      un = (rovx_m*nx + rovy_m*ny) / ro_m
      st%ro  (ci) = ro_m
      st%rovx(ci) = rovx_m - 2.0D0*ro_m*un*nx
      st%rovy(ci) = rovy_m - 2.0D0*ro_m*un*ny
      st%roe (ci) = roe_m

    end do
  end subroutine apply_ghost_bc

  ! -----------------------------------------------------------------------
  ! Per-cell CFL timestep.  Global minimum is used for explicit update.
  ! -----------------------------------------------------------------------
  subroutine compute_dt(mesh, av, st, dt_global)
    type(lod_mesh),  intent(in)    :: mesh
    type(t_appvars), intent(in)    :: av
    type(cf_state),  intent(inout) :: st
    real(8),         intent(out)   :: dt_global

    integer :: i
    real(8) :: gam, cfl, dx, dy, lam_x, lam_y, p_i, c_i

    gam = real(av%gam, 8)
    cfl = min(real(av%cfl, 8), 0.4D0)   ! Lax-Friedrichs 2D stability limit
    dt_global = 1.0D300

    do i = 1, st%ncells
      if (mesh%cells(i)%is_solid == 1) cycle
      p_i  = st%p(i)
      c_i  = sqrt(max(gam*p_i/st%ro(i), 0.0D0))
      dx   = mesh%cells(i)%xmax - mesh%cells(i)%xmin
      dy   = mesh%cells(i)%ymax - mesh%cells(i)%ymin
      lam_x = abs(st%vx(i)) + c_i
      lam_y = abs(st%vy(i)) + c_i
      st%dt(i) = cfl / max(lam_x/dx + lam_y/dy, 1.0D-30)
      dt_global = min(dt_global, st%dt(i))
    end do
  end subroutine compute_dt

  ! -----------------------------------------------------------------------
  ! One explicit Euler flux step.
  !
  ! For each fluid cell i, loop over its neighbours j via neigh_indices.
  ! The face flux (Lax-Friedrichs) is scaled by face_len and added to the
  ! increment of i; by Newton's 3rd law the same flux with opposite sign
  ! goes to j.  We skip j-side accumulation for solid cells.
  !
  ! After accumulating, conserved vars are updated:
  !   Q_new = Q_old - dt/area * sum(F * face_len)
  ! -----------------------------------------------------------------------
  subroutine euler_step(mesh, av, st)
    type(lod_mesh),  intent(in)    :: mesh
    type(t_appvars), intent(in)    :: av
    type(cf_state),  intent(inout) :: st

    integer :: i, k, j, offset
    real(8) :: gam, nx, ny, face_len
    real(8) :: f_ro, f_rovx, f_rovy, f_roe
    real(8) :: dtA_i

    gam = real(av%gam, 8)

    st%dro   = 0.0D0
    st%drovx = 0.0D0
    st%drovy = 0.0D0
    st%droe  = 0.0D0

    do i = 1, st%ncells
      if (mesh%cells(i)%is_solid == 1) cycle

      offset = mesh%cells(i)%neigh_offset

      do k = 0, mesh%cells(i)%neigh_count - 1
        j = mesh%neigh_indices(offset + k)

        ! Skip solid neighbours entirely — ghost BC handles wall condition
        if (mesh%cells(j)%is_solid == 1) cycle
        ! Only process each pair once: skip if j already handled its side
        if (j < i) cycle

        call face_geometry(mesh%cells(i), mesh%cells(j), face_len, nx, ny)

        call lax_friedrichs(                                    &
          st%ro(i), st%rovx(i), st%rovy(i), st%roe(i),         &
          st%ro(j), st%rovx(j), st%rovy(j), st%roe(j),         &
          nx, ny, gam,                                          &
          f_ro, f_rovx, f_rovy, f_roe)

        ! Accumulate signed flux (flux leaves i, enters j)
        st%dro  (i) = st%dro  (i) - f_ro   * face_len
        st%drovx(i) = st%drovx(i) - f_rovx * face_len
        st%drovy(i) = st%drovy(i) - f_rovy * face_len
        st%droe (i) = st%droe (i) - f_roe  * face_len

        st%dro  (j) = st%dro  (j) + f_ro   * face_len
        st%drovx(j) = st%drovx(j) + f_rovx * face_len
        st%drovy(j) = st%drovy(j) + f_rovy * face_len
        st%droe (j) = st%droe (j) + f_roe  * face_len
      end do
    end do

    ! Apply increments with per-cell local timestep
    do i = 1, st%ncells
      if (mesh%cells(i)%is_solid == 1) cycle
      dtA_i = st%dt(i) / st%area(i)
      st%ro  (i) = st%ro  (i) + dtA_i * st%dro  (i)
      st%rovx(i) = st%rovx(i) + dtA_i * st%drovx(i)
      st%rovy(i) = st%rovy(i) + dtA_i * st%drovy(i)
      st%roe (i) = st%roe (i) + dtA_i * st%droe (i)
    end do

  end subroutine euler_step

  ! -----------------------------------------------------------------------
  ! Farfield BC — matches the block solver's validated approach:
  !
  ! INLET (open sides on upstream half of domain):
  !   Relax boundary density toward rostag via rfin, then derive all other
  !   quantities isentropically from stagnation conditions — exactly as
  !   apply_bconds does in the block solver.
  !   ro is capped at 0.9999*rostag to prevent transient crashes.
  !
  ! OUTLET (open sides on downstream half):
  !   Fix static pressure to p_out only; density and momentum extrapolated
  !   from interior (1 characteristic enters from outside for subsonic flow).
  ! -----------------------------------------------------------------------
  subroutine apply_farfield_bc(mesh, av, bcs, st)
    type(lod_mesh),  intent(in)    :: mesh
    type(t_appvars), intent(in)    :: av
    type(t_bconds),  intent(in)    :: bcs
    type(cf_state),  intent(inout) :: st

    integer :: i, total_sides
    real(8) :: gam, cp, pstag, tstag, alpha, rfin, rostag, p_out
    real(8) :: ro_i, Tstatic, Vinlet, cx

    gam    = real(av%gam,    8)
    cp     = real(av%cp,     8)
    pstag  = real(bcs%pstag, 8)
    tstag  = real(bcs%tstag, 8)
    alpha  = real(bcs%alpha, 8) * acos(-1.0D0) / 180.0D0
    rfin   = real(bcs%rfin,  8)
    rostag = real(bcs%rostag,8)
    p_out  = real(bcs%p_out, 8)

    do i = 1, st%ncells
      if (mesh%cells(i)%is_solid == 1) cycle
      total_sides = int(mesh%cells(i)%side_count(1)) + &
                    int(mesh%cells(i)%side_count(2)) + &
                    int(mesh%cells(i)%side_count(3)) + &
                    int(mesh%cells(i)%side_count(4))
      if (total_sides >= 4) cycle

      cx = 0.5D0*(mesh%cells(i)%xmin + mesh%cells(i)%xmax)

      if (cx <= 1.5D0) then
        ! --- Inlet: stagnation-condition driven, same as block solver ---
        ! Relax boundary density toward interior value (filtered by rfin)
        ro_i = rfin * st%ro(i) + (1.0D0 - rfin) * st%ro(i)  ! = st%ro(i) first step
        ro_i = min(ro_i, 0.9999D0 * rostag)

        Tstatic = tstag * (ro_i / rostag)**(gam - 1.0D0)
        Vinlet  = sqrt(max(2.0D0 * cp * (tstag - Tstatic), 0.0D0))

        st%ro  (i) = ro_i
        st%rovx(i) = ro_i * Vinlet * cos(alpha)
        st%rovy(i) = ro_i * Vinlet * sin(alpha)
        st%roe (i) = ro_i * (cp/gam * Tstatic + 0.5D0 * Vinlet**2)
        ! p = pstag*(ro/rostag)^gam — consistent with isentropic inlet
      else
        ! --- Outlet: fix static pressure only ---
        st%roe(i) = p_out/(gam-1.0D0) + &
                    0.5D0*(st%rovx(i)**2 + st%rovy(i)**2) / max(st%ro(i), 1.0D-6)
      end if
    end do
  end subroutine apply_farfield_bc

  ! -----------------------------------------------------------------------
  ! Send flow state to the Qt GUI over the C bridge.
  ! Packs pointers into cf_state_c and calls emit_cf_state.
  ! -----------------------------------------------------------------------
  subroutine cf_state_to_qt(st)
    type(cf_state), intent(in), target :: st
    type(cf_state_c) :: sc

    sc%length = int(st%ncells, c_int)
    sc%ro     = c_loc(st%ro   (1))
    sc%rovx   = c_loc(st%rovx (1))
    sc%rovy   = c_loc(st%rovy (1))
    sc%roe    = c_loc(st%roe  (1))
    sc%p      = c_loc(st%p    (1))
    sc%vx     = c_loc(st%vx   (1))
    sc%vy     = c_loc(st%vy   (1))
    sc%hstag  = c_loc(st%hstag(1))

    call emit_cf_state(sc)
  end subroutine cf_state_to_qt

end module cf_solver_module


! =========================================================================
! C-callable entry point
! =========================================================================
subroutine curve_fill_solver(av_c, bcs_c, g_c) bind(C, name="curve_fill_solver")
  use iso_c_binding
  use types
  use conversion
  use mesh_gen
  use io_module
  use cf_solver_module
  implicit none

  type(t_appvars_c), intent(in)    :: av_c
  type(t_bconds_c),  intent(in)    :: bcs_c
  type(t_grid_c),    intent(inout) :: g_c   ! unused — block mesh grid passed by C++ harness

  type(t_appvars) :: av
  type(t_bconds)  :: bcs
  type(t_grid)    :: g_dummy   ! bconds_from_c needs a grid for inlet array sizing
  type(lod_mesh)  :: mesh
  type(cf_state)  :: st
  real(8)         :: dt, d_max, d_avg
  integer         :: step

  print *, 'entered curve_fill_solver'

  call appvars_from_c(av_c, av)
  call grid_from_c   (g_c, g_dummy)
  call bconds_from_c (bcs_c, bcs, g_dummy)

  call generate_cmesh(1000, 1000, "2412", 1.0D0, real(6.0,8), 3.0D0, 3.0D0, mesh)
  call lod_mesh_to_qt(mesh)

  call init_state   (mesh, av, bcs, st)
  call set_secondary(mesh, av, st)

  print *, 'Starting Euler iterations, nsteps=', av%nsteps

  do step = 1, av%nsteps

    ! BCs must be applied before flux step so boundary cells are
    ! correct when euler_step reads them as neighbours
    call apply_farfield_bc(mesh, av, bcs, st)
    call apply_ghost_bc   (mesh, st)
    call set_secondary    (mesh, av, st)
    call compute_dt       (mesh, av, st, dt)
    call euler_step       (mesh, av, st)

    if (step == 1) &
      print '(A,ES10.3)', ' dt step 1 = ', dt

    if (mod(step, 100) == 0) then
      d_max = maxval(abs(st%dro))
      d_avg = sum   (abs(st%dro)) / real(st%ncells, 8)
      print '(A,I6,A,ES10.3,A,ES10.3,A,ES10.3)', &
        ' step ', step, '  dro_max=', d_max, '  dro_avg=', d_avg, '  dt=', dt
      call cf_state_to_qt(st)
    end if

  end do

  print *, 'curve_fill_solver done'

end subroutine curve_fill_solver

module mesh_traverse
    use types
    use general_utils
    use mesh_utils
    implicit none

    ! ---------------------------------------------------------------------
    ! Curvature anchors: extra refinement sources at high-curvature points
    ! along the polygon.  Each anchor contributes its own log-distance term
    ! to the stop_level via min(); since min of Lipschitz functions is
    ! Lipschitz, the 2:1 balance invariant is preserved.
    ! ---------------------------------------------------------------------
    real(8), allocatable, save :: anchor_x(:), anchor_y(:)

    ! D_MIN for anchors.  Each factor of 2 over WALL_DMIN extends the
    ! finest-cell zone by one refinement level outward from the anchor.
    !   anchor_dmin = WALL_DMIN  → no extra refinement
    !   anchor_dmin = 2·WALL_DMIN → +1 level (cells half-size out to 2·WALL_DMIN)
    !   anchor_dmin = 4·WALL_DMIN → +2 levels, etc.
    real(8),              save :: anchor_dmin = 2.0D0

    ! Stop-level constants (kept module-scope so they're tunable from one place)
    real(8), parameter :: WALL_DMIN = 1.0D0             ! grid units
    real(8), parameter :: ANCHOR_KAPPA_FRAC = 0.10D0    ! select vertices with κ > frac · κ_max

    contains

    ! Build curvature anchors from the polygon.  Vertices with curvature above
    ! ANCHOR_KAPPA_FRAC · kappa_max are kept; each one becomes a fine-mesh
    ! source with its own anchor_dmin.  For NACA-style airfoils this picks
    ! out the leading edge cleanly.
    subroutine setup_curvature_anchors(poly, kappa_max)
        real(8), intent(in) :: poly(:,:)
        real(8), intent(in) :: kappa_max

        integer :: i, np, te_idx, count
        real(8) :: kappa_i, kappa_thresh, dist_dummy

        np = size(poly, 1)
        te_idx = maxloc(poly(:,1), 1)
        kappa_thresh = ANCHOR_KAPPA_FRAC * kappa_max

        if (allocated(anchor_x)) deallocate(anchor_x)
        if (allocated(anchor_y)) deallocate(anchor_y)

        ! Two-pass: count first, then allocate exact size and fill
        count = 0
        do i = 2, np - 1
            if (abs(i - te_idx) <= CURVATURE_STENCIL_W) cycle
            call curvature_at_idx(poly(i,1), poly(i,2), poly, i, dist_dummy, kappa_i)
            if (kappa_i > kappa_thresh) count = count + 1
        end do

        allocate(anchor_x(count), anchor_y(count))
        count = 0
        do i = 2, np - 1
            if (abs(i - te_idx) <= CURVATURE_STENCIL_W) cycle
            call curvature_at_idx(poly(i,1), poly(i,2), poly, i, dist_dummy, kappa_i)
            if (kappa_i > kappa_thresh) then
                count = count + 1
                anchor_x(count) = poly(i, 1)
                anchor_y(count) = poly(i, 2)
            end if
        end do

        print '(A,I5,A,F8.5,A,F6.3,A)', '  curvature anchors: ', count, &
            '  (κ_threshold=', kappa_thresh, ', D_MIN_anchor=', anchor_dmin, ')'
    end subroutine setup_curvature_anchors

    ! -----------------------------------------------------------------------
    ! Decide the coarsest level at which a cell may be a leaf.
    !
    ! Logarithmic-distance refinement:
    !
    !   stop_level = floor( log2(d / D_MIN) )   (clamped to [0, max_level])
    !
    ! Why this, instead of the older exponential blend?  Because log2 is
    ! Lipschitz with gradient 1 / (d · ln 2), the change in stop_level over
    ! one cell-width is bounded above by 1.  That bound is exactly what is
    ! needed for 2:1 AMR balance: adjacent cells can differ by at most one
    ! refinement level by construction — no post-pass split needed.
    !
    ! Geometrically: a level-L cell sits at d ≈ 2^L · D_MIN; level (L+2)
    ! sits at ≥ 2^(L+2) · D_MIN.  Their combined half-widths are only
    ! O(2^L · finest_size), which is much less than the 3 · 2^L · D_MIN
    ! gap, so a 4:1 jump is geometrically impossible.
    !
    ! Curvature: handled by separate anchor sources (see anchor_x/y above).
    ! Each anchor contributes its own log-distance term and the final
    ! stop_level is the min over wall + all anchors.  min of Lipschitz
    ! functions is Lipschitz, so 2:1 balance is preserved.
    !
    ! D_MIN is in grid units (n,m space).  finest cell size is base/2^max_level
    ! = max(n,m)/2^max_level ≈ 0.5 grid units, so D_MIN = 1.0 forces the
    ! finest level for any cell whose centre is within ~half the cell of the
    ! wall — i.e. for every cell that actually straddles the surface.
    ! -----------------------------------------------------------------------
    function calc_stop_level(dist, px, py, n, m) result(stop_level)
        real(8), intent(in) :: dist, px, py
        integer, intent(in) :: n, m
        integer :: stop_level

        integer :: max_level, i, sl_a
        real(8) :: d_eff, d_a, ln2

        ln2 = log(2.0_8)
        max_level = ceiling(log(real(max(n,m),8)) / ln2) + 1

        ! Wall log-distance term (Lipschitz → 2:1 balanced for the wall).
        d_eff = max(dist, WALL_DMIN)
        stop_level = floor(log(d_eff / WALL_DMIN) / ln2)

        ! Curvature anchors.  min over Lipschitz contributions is Lipschitz,
        ! so the resulting field still satisfies the 2:1 invariant.
        if (allocated(anchor_x)) then
            do i = 1, size(anchor_x)
                d_a = sqrt((px - anchor_x(i))**2 + (py - anchor_y(i))**2)
                d_a = max(d_a, anchor_dmin)
                sl_a = floor(log(d_a / anchor_dmin) / ln2)
                if (sl_a < stop_level) stop_level = sl_a
            end do
        end if

        stop_level = max(0, min(max_level, stop_level))
    end function calc_stop_level

    recursive subroutine traverse_ncells(x0, y0, xi, xj, yi, yj, level, n, m, poly, &
                                         dist_ref, kappa_min, kappa_max, ncells)
        implicit none
        real(8), intent(in) :: x0, y0, xi, xj, yi, yj
        integer, intent(in) :: level, n, m
        real(8), intent(in) :: poly(:,:)
        real(8), intent(in) :: dist_ref, kappa_min, kappa_max
        integer, intent(inout) :: ncells

        real(8) :: rx0, rx1, ry0, ry1
        real(8) :: px, py, dist, curvature
        integer :: pidx, stop_level
        logical :: inside

        rx0 = min(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        rx1 = max(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        ry0 = min(y0, y0 + xj, y0 + yj, y0 + xj + yj)
        ry1 = max(y0, y0 + xj, y0 + yj, y0 + xj + yj)

        if (rx1 <= 0.0_8 .or. real(m,8) <= rx0 .or. &
            ry1 <= 0.0_8 .or. real(n,8) <= ry0)  return

        px = x0 + 0.5_8*(xi + yi)
        py = y0 + 0.5_8*(xj + yj)

        inside = (px >= 0.0_8 .and. px < real(m,8) .and. py >= 0.0_8 .and. py < real(n,8))
        if (.not. inside) return

        call nearest_idx(px, py, poly, pidx)
        call dist_xy_to_xy(px, py, poly(pidx,1), poly(pidx,2), dist)
        call curvature_at_idx(px, py, poly, pidx, dist, curvature)

        stop_level = calc_stop_level(dist, px, py, n, m)

        if (level <= stop_level) then
            ncells = ncells + 1
            return
        end if

        call traverse_ncells(x0,                      y0,                      &
            yi/2.0_8, yj/2.0_8, xi/2.0_8, xj/2.0_8, level-1, n, m, poly,     &
            dist_ref, kappa_min, kappa_max, ncells)

        call traverse_ncells(x0 + xi/2.0_8,           y0 + xj/2.0_8,          &
            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8,  level-1, n, m, poly,    &
            dist_ref, kappa_min, kappa_max, ncells)

        call traverse_ncells(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8,  level-1, n, m, poly,    &
            dist_ref, kappa_min, kappa_max, ncells)

        call traverse_ncells(x0 + xi/2.0_8 + yi,      y0 + xj/2.0_8 + yj,     &
            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, poly, &
            dist_ref, kappa_min, kappa_max, ncells)

    end subroutine traverse_ncells

    recursive subroutine traverse_cells(x0, y0, xi, xj, yi, yj, level, n, m, poly, &
                                        dist_ref, kappa_min, kappa_max, cidx, hmesh)
        implicit none
        real(8), intent(in) :: x0, y0, xi, xj, yi, yj
        integer, intent(in) :: level, n, m
        real(8), intent(in) :: poly(:,:)
        real(8), intent(in) :: dist_ref, kappa_min, kappa_max
        integer, intent(inout) :: cidx
        type(helper_lod_mesh), intent(inout) :: hmesh

        real(8) :: rx0, rx1, ry0, ry1
        real(8) :: px, py, dist, curvature
        integer :: pidx, stop_level
        logical :: inside

        rx0 = min(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        rx1 = max(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        ry0 = min(y0, y0 + xj, y0 + yj, y0 + xj + yj)
        ry1 = max(y0, y0 + xj, y0 + yj, y0 + xj + yj)

        if (rx1 <= 0.0_8 .or. real(m,8) <= rx0 .or. &
            ry1 <= 0.0_8 .or. real(n,8) <= ry0)  return

        px = x0 + 0.5_8*(xi + yi)
        py = y0 + 0.5_8*(xj + yj)

        inside = (px >= 0.0_8 .and. px < real(m,8) .and. py >= 0.0_8 .and. py < real(n,8))
        if (.not. inside) return

        call nearest_idx(px, py, poly, pidx)
        call dist_xy_to_xy(px, py, poly(pidx,1), poly(pidx,2), dist)
        call curvature_at_idx(px, py, poly, pidx, dist, curvature)

        stop_level = calc_stop_level(dist, px, py, n, m)

        if (level <= stop_level) then
            hmesh%cells(cidx)%xmin = x0
            hmesh%cells(cidx)%xmax = x0 + xi + yi
            hmesh%cells(cidx)%ymin = y0
            hmesh%cells(cidx)%ymax = y0 + xj + yj
            hmesh%cells(cidx)%level = level

            hmesh%nearest_cell_poly_curvature(cidx) = curvature
            hmesh%nearest_cell_poly_distance(cidx)  = dist
            hmesh%nearest_cell_poly_idx(cidx)        = pidx
            cidx = cidx + 1
            return
        end if

        call traverse_cells(x0,                      y0,                      &
            yi/2.0_8, yj/2.0_8, xi/2.0_8, xj/2.0_8, level-1, n, m, poly,    &
            dist_ref, kappa_min, kappa_max, cidx, hmesh)

        call traverse_cells(x0 + xi/2.0_8,           y0 + xj/2.0_8,          &
            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8,  level-1, n, m, poly,   &
            dist_ref, kappa_min, kappa_max, cidx, hmesh)

        call traverse_cells(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8,  level-1, n, m, poly,   &
            dist_ref, kappa_min, kappa_max, cidx, hmesh)

        call traverse_cells(x0 + xi/2.0_8 + yi,      y0 + xj/2.0_8 + yj,     &
            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, poly, &
            dist_ref, kappa_min, kappa_max, cidx, hmesh)

    end subroutine traverse_cells

end module mesh_traverse

module mesh_alloc
    use types
    use mesh_utils
    use mesh_traverse
    implicit none

    contains

    subroutine alloc_ncells(n, m, poly, extra_global_levels, hmesh)
        implicit none
        integer, intent(in) :: n, m
        real(8), intent(in) :: poly(:,:)
        integer, intent(in), optional :: extra_global_levels
        type(helper_lod_mesh), intent(out) :: hmesh

        integer :: fine_bits, base_bits, extra_levels
        real(8) :: Nb
        real(8) :: dist_max, kappa_min, kappa_max, dist_ref

        extra_levels = 1
        if (present(extra_global_levels)) extra_levels = extra_global_levels

        base_bits = ceiling(log(real(max(n, m),8)) / log(2.0_8))
        fine_bits  = base_bits + extra_levels
        Nb         = real(max(n,m), 8)

        call poly_stats(poly, n, m, dist_max, kappa_min, kappa_max)
        dist_ref = dist_max * 0.25D0

        print *, 'poly_stats: dist_max=', dist_max, ' kappa=[', kappa_min, ',', kappa_max, ']'

        hmesh%length = 0
        call traverse_ncells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, poly, &
                              dist_ref, kappa_min, kappa_max, hmesh%length)

        allocate(hmesh%cells(hmesh%length))
        allocate(hmesh%nearest_cell_poly_idx(hmesh%length))
        allocate(hmesh%nearest_cell_poly_distance(hmesh%length))
        allocate(hmesh%nearest_cell_poly_curvature(hmesh%length))

    end subroutine alloc_ncells

    subroutine alloc_nindex(n, m, cells, nindex)
        implicit none
        integer, intent(in) :: n, m
        type(cell2d), intent(inout) :: cells(:)
        integer, intent(out) :: nindex
        nindex = 0
    end subroutine alloc_nindex

end module mesh_alloc

module mesh_build
    use types
    use mesh_utils
    use mesh_traverse
    implicit none

    contains

    subroutine build_cells(n, m, poly, extra_global_levels, hmesh)
        implicit none
        integer, intent(in) :: n, m
        real(8), intent(in) :: poly(:,:)
        integer, intent(in), optional :: extra_global_levels
        type(helper_lod_mesh), intent(inout) :: hmesh

        integer :: fine_bits, base_bits, extra_levels
        real(8) :: Nb
        real(8) :: dist_max, kappa_min, kappa_max, dist_ref
        real(8) :: tempxmin, tempxmax, tempymin, tempymax
        integer :: cidx, i

        extra_levels = 1
        if (present(extra_global_levels)) extra_levels = extra_global_levels

        base_bits = ceiling(log(real(max(n, m),8)) / log(2.0_8))
        fine_bits  = base_bits + extra_levels
        Nb         = real(max(n,m), 8)

        call poly_stats(poly, n, m, dist_max, kappa_min, kappa_max)
        dist_ref = dist_max * 0.25D0

        cidx = 1
        call traverse_cells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, poly, &
                             dist_ref, kappa_min, kappa_max, cidx, hmesh)

        ! fix min/max signs (Hilbert rotations can invert axes)
        do i = 1, hmesh%length
            tempxmin = min(hmesh%cells(i)%xmin, hmesh%cells(i)%xmax)
            tempxmax = max(hmesh%cells(i)%xmin, hmesh%cells(i)%xmax)
            tempymin = min(hmesh%cells(i)%ymin, hmesh%cells(i)%ymax)
            tempymax = max(hmesh%cells(i)%ymin, hmesh%cells(i)%ymax)

            hmesh%cells(i)%xmin = tempxmin
            hmesh%cells(i)%xmax = tempxmax
            hmesh%cells(i)%ymin = tempymin
            hmesh%cells(i)%ymax = tempymax
        end do

    end subroutine build_cells

end module mesh_build

! =============================================================================
! Single-pass mesh builder: dynamic buffer, one traversal, MOVE_ALLOC on trim.
! Equivalent to alloc_ncells + build_cells but visits the recursion tree once.
! =============================================================================
module mesh_build_v2
    use types
    use mesh_utils
    use mesh_traverse
    implicit none

    contains

    recursive subroutine traverse_build(x0, y0, xi, xj, yi, yj, level, n, m, poly, &
                                        dist_ref, kappa_min, kappa_max, hmesh, capacity)
        implicit none
        real(8), intent(in)    :: x0, y0, xi, xj, yi, yj
        integer, intent(in)    :: level, n, m
        real(8), intent(in)    :: poly(:,:)
        real(8), intent(in)    :: dist_ref, kappa_min, kappa_max
        type(helper_lod_mesh), intent(inout) :: hmesh
        integer,               intent(inout) :: capacity

        real(8) :: rx0, rx1, ry0, ry1, px, py, dist, curvature
        integer :: pidx, stop_level, new_cap, k
        logical :: inside
        type(cell2d), allocatable         :: tmp_cells(:)
        integer,      allocatable         :: tmp_idx(:)
        real(8),      allocatable         :: tmp_dist(:), tmp_curv(:)

        rx0 = min(x0, x0+xi, x0+yi, x0+xi+yi)
        rx1 = max(x0, x0+xi, x0+yi, x0+xi+yi)
        ry0 = min(y0, y0+xj, y0+yj, y0+xj+yj)
        ry1 = max(y0, y0+xj, y0+yj, y0+xj+yj)

        if (rx1 <= 0.0_8 .or. real(m,8) <= rx0 .or. &
            ry1 <= 0.0_8 .or. real(n,8) <= ry0) return

        px = x0 + 0.5_8*(xi + yi)
        py = y0 + 0.5_8*(xj + yj)

        inside = (px >= 0.0_8 .and. px < real(m,8) .and. py >= 0.0_8 .and. py < real(n,8))
        if (.not. inside) return

        call nearest_idx(px, py, poly, pidx)
        call dist_xy_to_xy(px, py, poly(pidx,1), poly(pidx,2), dist)
        call curvature_at_idx(px, py, poly, pidx, dist, curvature)

        stop_level = calc_stop_level(dist, px, py, n, m)

        if (level <= stop_level) then
            ! grow buffer if needed (double capacity)
            if (hmesh%length == capacity) then
                new_cap = capacity * 2
                allocate(tmp_cells(new_cap), tmp_idx(new_cap), tmp_dist(new_cap), tmp_curv(new_cap))
                tmp_cells(1:capacity) = hmesh%cells(1:capacity)
                tmp_idx  (1:capacity) = hmesh%nearest_cell_poly_idx(1:capacity)
                tmp_dist (1:capacity) = hmesh%nearest_cell_poly_distance(1:capacity)
                tmp_curv (1:capacity) = hmesh%nearest_cell_poly_curvature(1:capacity)
                call move_alloc(tmp_cells, hmesh%cells)
                call move_alloc(tmp_idx,   hmesh%nearest_cell_poly_idx)
                call move_alloc(tmp_dist,  hmesh%nearest_cell_poly_distance)
                call move_alloc(tmp_curv,  hmesh%nearest_cell_poly_curvature)
                capacity = new_cap
            end if

            hmesh%length = hmesh%length + 1
            k = hmesh%length
            hmesh%cells(k)%xmin  = x0
            hmesh%cells(k)%xmax  = x0 + xi + yi
            hmesh%cells(k)%ymin  = y0
            hmesh%cells(k)%ymax  = y0 + xj + yj
            hmesh%cells(k)%level = level
            hmesh%nearest_cell_poly_curvature(k) = curvature
            hmesh%nearest_cell_poly_distance(k)  = dist
            hmesh%nearest_cell_poly_idx(k)       = pidx
            return
        end if

        call traverse_build(x0,                       y0,                       &
            yi/2.0_8, yj/2.0_8, xi/2.0_8, xj/2.0_8,  level-1, n, m, poly,     &
            dist_ref, kappa_min, kappa_max, hmesh, capacity)

        call traverse_build(x0 + xi/2.0_8,            y0 + xj/2.0_8,           &
            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8,   level-1, n, m, poly,    &
            dist_ref, kappa_min, kappa_max, hmesh, capacity)

        call traverse_build(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8,   level-1, n, m, poly,    &
            dist_ref, kappa_min, kappa_max, hmesh, capacity)

        call traverse_build(x0 + xi/2.0_8 + yi,       y0 + xj/2.0_8 + yj,      &
            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, poly,  &
            dist_ref, kappa_min, kappa_max, hmesh, capacity)

    end subroutine traverse_build

end module mesh_build_v2

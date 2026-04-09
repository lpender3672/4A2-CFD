module mesh_traverse
    use types
    use general_utils
    use mesh_utils
    implicit none
    contains

    ! -----------------------------------------------------------------------
    ! Decide the coarsest level at which a cell may be a leaf.
    !
    !   w_dist  = exp(-dist / dist_ref)            1 at surface, decays outward
    !   w_kappa = (kappa - kappa_min) / range      0 = flat, 1 = sharpest
    !
    !   refinement = alpha*w_dist + (1-alpha)*w_kappa   (alpha ~ 0.7)
    !   stop_level = round( (1 - refinement) * max_level )
    !
    ! A lower stop_level means finer cells (we recurse until level <= stop_level).
    ! max_level is derived from n,m identically to fine_bits in the callers so it
    ! does not need to be threaded through the recursion as a parameter.
    ! -----------------------------------------------------------------------
    pure function calc_stop_level(dist, kappa, dist_ref, kappa_min, kappa_max, n, m) &
        result(stop_level)

        real(8), intent(in) :: dist, kappa
        real(8), intent(in) :: dist_ref, kappa_min, kappa_max
        integer, intent(in) :: n, m
        integer :: stop_level

        real(8), parameter :: ALPHA = 0.7D0, EPS = 1.0D-12
        real(8) :: w_dist, w_kappa, refinement
        integer :: max_level

        max_level = ceiling(log(real(max(n,m),8)) / log(2.0_8)) + 1

        w_dist = exp(-dist / max(dist_ref, EPS))

        if (kappa_max > kappa_min + EPS) then
            w_kappa = (kappa - kappa_min) / (kappa_max - kappa_min)
        else
            w_kappa = 0.0D0
        end if
        w_kappa = max(0.0D0, min(1.0D0, w_kappa))

        refinement = ALPHA * w_dist + (1.0D0 - ALPHA) * w_kappa
        stop_level = nint((1.0D0 - refinement) * real(max_level, 8))
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

        stop_level = calc_stop_level(dist, curvature, dist_ref, kappa_min, kappa_max, n, m)

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

        stop_level = calc_stop_level(dist, curvature, dist_ref, kappa_min, kappa_max, n, m)

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

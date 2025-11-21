module mesh_traverse
    use types
    use general_utils
    use mesh_utils
    implicit none
    contains

    pure function calc_stop_level(dist, curvature) result(level)
        implicit none
        real(8), intent(in) :: dist, curvature
        integer :: level
        real(8) :: d_thresh, c_thresh

        ! thresholds (tunable)
        d_thresh = 5.0_8
        c_thresh = 0.1_8

        ! simple heuristic: higher level (finer) for smaller distance and higher curvature
        level = int( max(0.0_8, min(10.0_8, &
                    10.0_8 - (dist / d_thresh)*5.0_8 - (curvature / c_thresh)*5.0_8)) )

    end function calc_stop_level

    recursive subroutine traverse_ncells(x0, y0, xi, xj, yi, yj, level, n, m, poly, ncells)
        implicit none
        real(8), intent(in) :: x0, y0, xi, xj, yi, yj
        integer, intent(in) :: level, n, m
        real(8), intent(in) :: poly(:,:) ! (np,2)
        integer, intent(inout) :: ncells

        real(8) :: rx0, rx1, ry0, ry1
        real(8) :: px, py, dist, curvature
        integer :: pidx, stop_level
        logical :: inside

        ! Compute region bounding box (fine grid coords)
        rx0 = min(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        rx1 = max(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        ry0 = min(y0, y0 + xj, y0 + yj, y0 + xj + yj)
        ry1 = max(y0, y0 + xj, y0 + yj, y0 + xj + yj)

        ! Skip if outside domain
        if (rx1 <= 0.0_8 .or. real(m,8) <= rx0 .or. &
            ry1 <= 0.0_8 .or. real(n,8) <= ry0)  return

        ! Region center
        px = x0 + 0.5_8*(xi + yi)
        py = y0 + 0.5_8*(xj + yj)

        ! Domain bounds check
        inside = (px >= 0.0_8 .and. px < real(m,8) .and. py >= 0.0_8 .and. py < real(n,8))
        if (.not. inside) return

        ! Lookup LOD index
        call nearest_idx(px, py, poly, pidx)
        call dist_xy_to_xy(px, py, poly(pidx,1), poly(pidx,2), dist)
        call curvature_at_idx(px, py, poly, pidx, dist, curvature)

        stop_level = calc_stop_level(dist, curvature)

        if (level <= stop_level) then
            ncells = ncells + 1
            return
        end if

        ! Recurse in Hilbert order (same pattern as Python)
        call traverse_ncells(x0,                  y0,                  yi/2.0_8,  yj/2.0_8,  xi/2.0_8,  xj/2.0_8,  level-1, n, m, poly, ncells)

        call traverse_ncells(x0 + xi/2.0_8,       y0 + xj/2.0_8,       xi/2.0_8,  xj/2.0_8,  yi/2.0_8,  yj/2.0_8,  level-1, n, m, poly, ncells)

        call traverse_ncells(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
                            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8, level-1, n, m, poly, ncells)

        call traverse_ncells(x0 + xi/2.0_8 + yi,  y0 + xj/2.0_8 + yj, &
                            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, poly, ncells)

    end subroutine traverse_ncells

    recursive subroutine traverse_cells(x0, y0, xi, xj, yi, yj, level, n, m, poly, cidx, hmesh)
        implicit none
        real(8), intent(in) :: x0, y0, xi, xj, yi, yj
        integer, intent(in) :: level, n, m
        real(8), intent(in) :: poly(:,:) ! (np,2)
        integer, intent(inout) :: cidx
        type(helper_lod_mesh), intent(inout) :: hmesh

        real(8) :: rx0, rx1, ry0, ry1
        real(8) :: px, py, dist, curvature
        integer :: pidx, stop_level
        logical :: inside

        ! Compute region bounding box (fine grid coords)
        rx0 = min(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        rx1 = max(x0, x0 + xi, x0 + yi, x0 + xi + yi)
        ry0 = min(y0, y0 + xj, y0 + yj, y0 + xj + yj)
        ry1 = max(y0, y0 + xj, y0 + yj, y0 + xj + yj)

        ! Skip if outside domain
        if (rx1 <= 0.0_8 .or. real(m,8) <= rx0 .or. &
            ry1 <= 0.0_8 .or. real(n,8) <= ry0)  return

        ! Region center
        px = x0 + 0.5_8*(xi + yi)
        py = y0 + 0.5_8*(xj + yj)

        ! Domain bounds check
        inside = (px >= 0.0_8 .and. px < real(m,8) .and. py >= 0.0_8 .and. py < real(n,8))
        if (.not. inside) return

        ! Lookup LOD index
        call nearest_idx(px, py, poly, pidx)
        call dist_xy_to_xy(px, py, poly(pidx,1), poly(pidx,2), dist)
        call curvature_at_idx(px, py, poly, pidx, dist, curvature)

        stop_level = calc_stop_level(dist, curvature)

        if (level <= stop_level) then
            ! fill cells array her
            hmesh%cells(cidx)%xmin = x0
            hmesh%cells(cidx)%xmax = x0 + xi + yi
            hmesh%cells(cidx)%ymin = y0
            hmesh%cells(cidx)%ymax = y0 + xj + yj
            hmesh%cells(cidx)%level = level

            hmesh%nearest_cell_poly_curvature(cidx) = curvature
            hmesh%nearest_cell_poly_distance(cidx) = dist
            hmesh%nearest_cell_poly_idx(cidx) = pidx
            ! etc
            cidx = cidx + 1
            return
        end if

        ! Recurse in Hilbert order (same pattern as Python)
        call traverse_cells(x0,                  y0,                  yi/2.0_8,  yj/2.0_8,  xi/2.0_8,  xj/2.0_8,  level-1, n, m, poly, cidx, hmesh)

        call traverse_cells(x0 + xi/2.0_8,       y0 + xj/2.0_8,       xi/2.0_8,  xj/2.0_8,  yi/2.0_8,  yj/2.0_8,  level-1, n, m, poly, cidx, hmesh)

        call traverse_cells(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
                            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8, level-1, n, m, poly, cidx, hmesh)

        call traverse_cells(x0 + xi/2.0_8 + yi,  y0 + xj/2.0_8 + yj, &
                            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, poly, cidx, hmesh)

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
        real(8), intent(in) :: poly(:,:) ! (np,2)
        integer, intent(in), optional :: extra_global_levels
        type(helper_lod_mesh), intent(out) :: hmesh

        integer :: fine_bits, base_bits
        real(8) :: Nb
        integer :: extra_levels
        integer, parameter :: unit_out = 33

        ! For now, ignore airfoil_poly (no geometry checks)
        extra_levels = 1
        if (present(extra_global_levels)) extra_levels = extra_global_levels

        base_bits = ceiling(log(real(max(n, m),8)) / log(2.0_8))
        fine_bits = base_bits + extra_levels
        ! Nb = 2.0_8**fine_bits
        Nb = real(max(n,m), 8)

        hmesh%length = 0
        
        call traverse_ncells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, poly, hmesh%length)

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

        ! traverse neighbouring cells incrementing nindex

        ! this will fill up neigh_index of all neighbours
        ! then theres two arrays of length ncells neigh_start and neigh_count

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
        real(8), intent(in) :: poly(:,:) ! (np,2)
        integer, intent(in), optional :: extra_global_levels
        type(helper_lod_mesh), intent(inout) :: hmesh

        integer :: fine_bits, base_bits
        real(8) :: Nb
        real(8) :: tempxmin, tempxmax, tempymin, tempymax
        integer :: extra_levels
        integer :: cidx, i
        
        extra_levels = 1
        if (present(extra_global_levels)) extra_levels = extra_global_levels

        base_bits = ceiling(log(real(max(n, m),8)) / log(2.0_8))
        fine_bits = base_bits + extra_levels
        ! Nb = 2.0_8**fine_bits
        Nb = real(max(n,m), 8)

        cidx = 1
        
        call traverse_cells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, poly, cidx, hmesh)

        ! fix min and max
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
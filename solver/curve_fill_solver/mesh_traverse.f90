module mesh_traverse
    use types
    use mesh_utils
    implicit none
    contains

    recursive subroutine traverse_ncells(x0, y0, xi, xj, yi, yj, level, n, m, ILOD, ncells)
        implicit none
        real(8), intent(in) :: x0, y0, xi, xj, yi, yj
        integer, intent(in) :: level, n, m
        integer, intent(in) :: ILOD(n,m)
        integer, intent(inout) :: ncells

        real(8) :: rx0, rx1, ry0, ry1
        real(8) :: px, py
        integer :: stop_level
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
        stop_level = lod_from_xy(px, py, n, m, ILOD)

        if (level <= stop_level) then
            ncells = ncells + 1
            return
        end if

        ! Recurse in Hilbert order (same pattern as Python)
        call traverse_ncells(x0,                  y0,                  yi/2.0_8,  yj/2.0_8,  xi/2.0_8,  xj/2.0_8,  level-1, n, m, ILOD, ncells)

        call traverse_ncells(x0 + xi/2.0_8,       y0 + xj/2.0_8,       xi/2.0_8,  xj/2.0_8,  yi/2.0_8,  yj/2.0_8,  level-1, n, m, ILOD, ncells)

        call traverse_ncells(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
                            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8, level-1, n, m, ILOD, ncells)

        call traverse_ncells(x0 + xi/2.0_8 + yi,  y0 + xj/2.0_8 + yj, &
                            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, ILOD, ncells)

    end subroutine traverse_ncells

    recursive subroutine traverse_cells(x0, y0, xi, xj, yi, yj, level, n, m, ILOD, cidx, cells)
        implicit none
        real(8), intent(in) :: x0, y0, xi, xj, yi, yj
        integer, intent(in) :: level, n, m
        integer, intent(in) :: ILOD(n,m)
        integer, intent(inout) :: cidx
        type(cell2d), intent(inout) :: cells(:)

        real(8) :: rx0, rx1, ry0, ry1
        real(8) :: px, py
        integer :: stop_level
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
        stop_level = lod_from_xy(px, py, n, m, ILOD)

        if (level <= stop_level) then
            ! fill cells array her
            cells(cidx)%xmin = x0
            cells(cidx)%xmax = x0 + xi + yi
            cells(cidx)%ymin = y0
            cells(cidx)%ymax = y0 + xj + yj
            cells(cidx)%level = level

            ! etc
            cidx = cidx + 1
            return
        end if

        ! Recurse in Hilbert order (same pattern as Python)
        call traverse_cells(x0,                  y0,                  yi/2.0_8,  yj/2.0_8,  xi/2.0_8,  xj/2.0_8,  level-1, n, m, ILOD, cidx, cells)

        call traverse_cells(x0 + xi/2.0_8,       y0 + xj/2.0_8,       xi/2.0_8,  xj/2.0_8,  yi/2.0_8,  yj/2.0_8,  level-1, n, m, ILOD, cidx, cells)

        call traverse_cells(x0 + xi/2.0_8 + yi/2.0_8, y0 + xj/2.0_8 + yj/2.0_8, &
                            xi/2.0_8, xj/2.0_8, yi/2.0_8, yj/2.0_8, level-1, n, m, ILOD, cidx, cells)

        call traverse_cells(x0 + xi/2.0_8 + yi,  y0 + xj/2.0_8 + yj, &
                            -yi/2.0_8, -yj/2.0_8, -xi/2.0_8, -xj/2.0_8, level-1, n, m, ILOD, cidx, cells)

    end subroutine traverse_cells

end module mesh_traverse

module mesh_alloc
    use types
    use mesh_utils
    use mesh_traverse
    implicit none

    contains

    subroutine alloc_ncells(n, m, ILOD, extra_global_levels, ncells, cells)
        implicit none
        integer, intent(in) :: n, m
        integer, intent(in) :: ILOD(n,m)             ! precomputed LOD map
        integer, intent(in), optional :: extra_global_levels
        integer, intent(out) :: ncells
        type(cell2d), allocatable, intent(out) :: cells(:)

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

        ncells = 0
        
        call traverse_ncells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, ILOD, ncells)

        allocate(cells(ncells))

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

    subroutine build_cells(n, m, ILOD, extra_global_levels, cells)
        implicit none
        integer, intent(in) :: n, m
        integer, intent(in) :: ILOD(n,m)             ! precomputed LOD map
        integer, intent(in), optional :: extra_global_levels
        type(cell2d), intent(out) :: cells(:)

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
        
        call traverse_cells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, ILOD, cidx, cells)

        ! fix min and max
        do i = 1, size(cells)
            tempxmin = min(cells(i)%xmin, cells(i)%xmax)
            tempxmax = max(cells(i)%xmin, cells(i)%xmax)
            tempymin = min(cells(i)%ymin, cells(i)%ymax)
            tempymax = max(cells(i)%ymin, cells(i)%ymax)

            cells(i)%xmin = tempxmin
            cells(i)%xmax = tempxmax
            cells(i)%ymin = tempymin
            cells(i)%ymax = tempymax
        end do

    end subroutine build_cells

end module mesh_build
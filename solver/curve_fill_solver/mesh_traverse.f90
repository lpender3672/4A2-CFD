module mesh_traverse

    implicit none
    contains
    
end module mesh_traverse

module mesh_alloc
    use mesh_utils
    use types
    implicit none

    contains

    subroutine calc_ncells(n, m, ILOD, extra_global_levels, ncells)
        implicit none
        integer, intent(in) :: n, m
        integer, intent(in) :: ILOD(n,m)             ! precomputed LOD map
        integer, intent(in), optional :: extra_global_levels
        integer, intent(out) :: ncells

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
        
        open(unit_out, file="hilbert_cells.dat", status="replace", action="write", form="formatted")

        call traverse_ncells(0.0_8, 0.0_8, Nb, 0.0_8, 0.0_8, Nb, fine_bits, n, m, ILOD, ncells)

        close(unit_out)

    end subroutine calc_ncells


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
        integer, parameter :: unit_out = 33

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
            write(unit_out,'(I6,1X,F12.6,1X,F12.6,1X,I3)') level, px, py, stop_level
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

    subroutine calc_nindex(n, m, cells, nindex)

        implicit none
        integer, intent(in) :: n, m
        type(cell2d), intent(inout) :: cells(:)
        integer, intent(out) :: nindex

        ! traverse neighbouring cells incrementing nindex

        ! this will fill up neigh_index of all neighbours
        ! then theres two arrays of length ncells neigh_start and neigh_count

    end subroutine calc_nindex

    integer function lod_from_xy(x, y, n, m, ILOD) result(lod)
        implicit none
        real(8), intent(in) :: x, y
        integer, intent(in) :: n, m
        integer, intent(in) :: ILOD(n,m)
        integer :: ix, iy

        ix = int(max(0.0_8, min(real(m-1,8), real(floor(x * real(m-1,8) / real(m,8)), 8))))
        iy = int(max(0.0_8, min(real(n-1,8), real(floor(y * real(n-1,8) / real(n,8)), 8))))

        lod = ILOD(iy+1, ix+1)  ! fortran arrays are 1-based
    end function lod_from_xy

end module mesh_alloc
module mesh_gen
    use types
    use mesh_utils
    use mesh_alloc
    ! use iso_fortran_env, only: real64

    implicit none
    contains

    subroutine generate_cmesh(n, m, naca_code, chord, aoa, cells)
        implicit none
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: naca_code
        real(8), intent(in) :: chord, aoa
        type(cell2d), allocatable, intent(out) :: cells(:)

        integer :: ILOD(n,m)
        integer :: ncells, max_level, i
        real(8), allocatable :: foil(:,:), foil_t(:,:)


        ! Build an airfoil in domain units as you like:
        call naca4_airfoil("2412", 400, .true., foil)
        allocate(foil_t(size(foil,1),2))
        call transform_airfoil(foil, chord, aoa, origin=[real(200, 8), real(260, 8)], poly_out=foil_t)

        ! Generate airfoil geometry

        max_level = 5
        ! Build LOD map based on airfoil geometry and NACA code
        call calc_lod(n, m, foil_t, max_level, ILOD)

        ! print ILOD

        do i=1, n
            write(*,'(100i4)') ILOD(i,1:m)
        end do

        ! Calculate number of cells needed
        call calc_ncells(n, m, ILOD, 1, ncells)

        allocate(cells(ncells))

        print *, 'Number of cells to allocate: ', ncells

    end subroutine generate_cmesh
    
end module mesh_gen
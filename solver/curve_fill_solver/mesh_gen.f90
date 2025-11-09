module mesh_gen
    use types
    use mesh_utils
    use mesh_alloc
    use mesh_build
    ! use iso_fortran_env, only: real64

    implicit none
    contains

    subroutine generate_cmesh(n, m, naca_code, chord, aoa, mesh)
        implicit none
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: naca_code
        real(8), intent(in) :: chord, aoa
        !type(cell2d), allocatable, intent(out) :: cells(:)
        type(lod_mesh), intent(out) :: mesh

        integer :: ILOD(n,m)
        integer :: ncells, max_level, i, j, unit
        real(8), allocatable :: foil(:,:), foil_t(:,:)


        ! Build an airfoil in domain units as you like:
        call naca4_airfoil("2412", 400, .true., foil)
        allocate(foil_t(size(foil,1),2))
        call transform_airfoil(foil, chord, aoa, origin=[real(n*0.5, 8), real(m*0.5, 8)], poly_out=foil_t)

        ! Generate airfoil geometry

        max_level = 5
        ! Build LOD map based on airfoil geometry and NACA code
        call calc_lod(n, m, foil_t, max_level, ILOD)

        ! Calculate number of cells needed
        call alloc_ncells(n, m, ILOD, 1, mesh%length, mesh%cells)

        call build_cells(n, m, ILOD, 1, mesh%cells)

        print *, 'Number of cells to allocate: ', mesh%length

    end subroutine generate_cmesh
    
end module mesh_gen
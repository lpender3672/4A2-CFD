module mesh_gen
    use types
    use mesh_utils
    use mesh_alloc
    use mesh_build
    use neighbouring
    ! use iso_fortran_env, only: real64

    implicit none
    contains

    subroutine write_mesh_csv(mesh, filename)
        type(lod_mesh), intent(in) :: mesh
        character(*),   intent(in) :: filename

        integer :: i
        integer :: unit

        open(newunit=unit, file=filename, status='replace', action='write', form='formatted')
        write(unit,'(A)') 'i,xmin,xmax,ymin,ymax'
        do i = 1, mesh%length
            write(unit,'(I6,1X,4ES20.10)') i, mesh%cells(i)%xmin, mesh%cells(i)%xmax, &
                                            mesh%cells(i)%ymin, mesh%cells(i)%ymax
        end do
        close(unit)
    end subroutine write_mesh_csv

    subroutine generate_cmesh(n, m, naca_code, chord, aoa, mesh)
        implicit none
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: naca_code
        real(8), intent(in) :: chord, aoa

        type(lod_mesh) :: full_mesh
        type(lod_mesh), intent(out) :: mesh

        integer :: ILOD(n,m)
        integer :: ncells, max_level, i, j, unit
        real(8), allocatable :: foil(:,:), foil_t(:,:), DIST(:,:), KAPPA(:,:)


        ! Build an airfoil in domain units as you like:
        call naca4_airfoil("2412", 400, .true., foil)
        allocate(foil_t(size(foil,1),2))
        call transform_airfoil(foil, chord, aoa, origin=[real(n*0.5, 8), real(m*0.5, 8)], poly_out=foil_t)

        ! Generate airfoil geometry

        max_level = 8
        ! Build LOD map based on airfoil geometry and NACA code
        ! also gives distance and curvature maps for walls
        call calc_lod(n, m, foil_t, max_level, DIST, KAPPA, ILOD)

        ! Calculate full mesh
        call alloc_ncells(n, m, ILOD, 1, full_mesh%length, full_mesh%cells)
        call build_cells(n, m, ILOD, 1, full_mesh%cells)

        print *, 'Full mesh cells allocated: ', full_mesh%length
        !call write_mesh_csv(full_mesh, 'mesh.csv')

        ! now function to reduce mesh
        call build_walls(full_mesh, foil_t, DIST, KAPPA, mesh)

        print *, 'Reduced mesh cells allocated: ', mesh%length

        ! now build indicies
        call build_indicies(mesh)

        print *, size(mesh%neigh_indices)
        
        call print_first_five_neighbors(mesh)

    end subroutine generate_cmesh

    subroutine build_walls(full_mesh, poly, distance, curvature, reduced_mesh)
        type(lod_mesh), intent(in) :: full_mesh

        ! precomputed distance and curvature fields used for lod decision
        real(8), intent(in) :: poly(:,:), distance(:,:), curvature(:,:)

        type(lod_mesh), intent(out) :: reduced_mesh
        real(8) :: phis(full_mesh%length)

        ! thresholds to decide whether to use linear or curved intersection
        real(8) :: distance_threshold, curvature_threshold
        integer :: i

        curvature_threshold = 1d-1 ! idk

        phis = 0.0d0

        do i=1, full_mesh%length
            ! if we're not within a cell size, skip
            distance_threshold = max(full_mesh%cells(i)%xmax - full_mesh%cells(i)%xmin, &
                                     full_mesh%cells(i)%ymax - full_mesh%cells(i)%ymin)

            if (interp_from_cell(full_mesh%cells(i), distance) > distance_threshold) cycle

            ! if curvature below arbitrary threshold
            if (interp_from_cell(full_mesh%cells(i), curvature) < curvature_threshold) then
                ! simple straight line intersection
                call linear_cell_phi(full_mesh%cells(i), poly, phis(i))
            else
                ! complex curve intersection
                call poly_cell_phi(full_mesh%cells(i), poly, phis(i))
            end if
        end do

        ! loop through phis and get length of reduced cells
        ! allocate cells
        ! fill reduced_mesh
        ! calculate normals as well, maybe with phi or maybe after

        reduced_mesh = full_mesh
    end subroutine

    subroutine linear_cell_phi(cell, poly, phi)
        type(cell2d), intent(in) :: cell
        real(8), intent(in) :: poly(:,:)
        real(8), intent(out) :: phi

        !TODO

        phi = 0.0d0
    end subroutine linear_cell_phi

    subroutine poly_cell_phi(cell, poly, phi)
        type(cell2d), intent(in) :: cell
        real(8), intent(in) :: poly(:,:)
        real(8), intent(out) :: phi

        !TODO

        phi = 0.0d0
    end subroutine poly_cell_phi
    
end module mesh_gen
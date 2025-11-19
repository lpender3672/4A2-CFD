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
        type(helper_lod_mesh), intent(in) :: mesh
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

    subroutine generate_cmesh(n, m, naca_code, chord, aoa, fmesh)
        implicit none
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: naca_code
        real(8), intent(in) :: chord, aoa

        type(helper_lod_mesh) :: hmesh ! helper mesh
        type(lod_mesh), intent(out) :: fmesh ! final mesh

        integer :: ILOD(n,m)
        integer :: ncells, max_level, i, j, unit
        real(8), allocatable :: foil(:,:), foil_t(:,:), DIST(:,:), KAPPA(:,:)


        ! Build an airfoil in domain units as you like:
        call naca4_airfoil("2412", 400, .true., foil)
        allocate(foil_t(size(foil,1),2))
        call transform_airfoil(foil, chord, aoa, origin=[real(n*0.5, 8), real(m*0.5, 8)], poly_out=foil_t)

        print *, 'Airfoil generated with ', size(foil_t,1), ' points.'

        ! Generate airfoil geometry

        max_level = 4
        ! Build LOD map based on airfoil geometry and NACA code
        ! also gives distance and curvature maps for walls
        ! call calc_lod(n, m, foil_t, max_level, DIST, KAPPA, ILOD)

        ! Calculate full mesh
        call alloc_ncells(n, m, foil_t, 1, hmesh)
        call build_cells(n, m, foil_t, 1, hmesh)

        print *, 'Full mesh cells allocated: ', hmesh%length
        !call write_mesh_csv(hmesh, 'mesh.csv')

        ! now function to reduce mesh
        call build_walls(hmesh, foil_t, fmesh)

        print *, 'Reduced mesh cells allocated: ', fmesh%length

        ! now build indicies
        call build_indicies(fmesh)

        print *, size(fmesh%neigh_indices)
        
        call print_first_five_neighbors(fmesh)

    end subroutine generate_cmesh

    subroutine build_walls(hmesh, poly, fmesh)
        type(helper_lod_mesh), intent(inout) :: hmesh

        ! precomputed distance and curvature fields used for lod decision
        real(8), intent(in) :: poly(:,:)

        type(lod_mesh), intent(out) :: fmesh
        real(8) :: phis(hmesh%length)

        ! thresholds to decide whether to use linear or curved intersection
        real(8) :: distance_threshold, curvature_threshold
        real(8) :: line_point(2), line_dir(2)
        integer :: i, cidx, widx

        real(8) :: EPS
        EPS = 1.0d-16

        curvature_threshold = 1d-1 ! idk

        phis = 0.0d0

        do i=1, hmesh%length
            ! if we're not within a cell size, skip
            distance_threshold = max(hmesh%cells(i)%xmax - hmesh%cells(i)%xmin, &
                                     hmesh%cells(i)%ymax - hmesh%cells(i)%ymin)

            if (hmesh%nearest_cell_poly_distance(i) > distance_threshold) cycle

            ! if curvature below arbitrary threshold
            if (hmesh%nearest_cell_poly_curvature(i) < curvature_threshold) then
                ! simple straight line intersection
                call line_cell_phi(hmesh%cells(i), line_point, line_dir, phis(i))
            else
                ! complex curve intersection
                call poly_cell_phi(hmesh%cells(i), poly, phis(i))
            end if

            if (phis(i) + EPS < 1.0D0) then
                fmesh%length = fmesh%length + 1
                if (phis(i) > EPS) then
                    fmesh%wall_count = fmesh%wall_count + 1
                end if
            end if
        end do

        allocate(fmesh%cells(fmesh%length))
        allocate(fmesh%wall_indices(fmesh%length))
        allocate(fmesh%solid_fractions(fmesh%wall_count))
        allocate(fmesh%wall_normals(fmesh%wall_count, 2))

        fmesh%wall_indices = 0
        cidx = 0
        widx = 0

        ! loop through phis and get length of reduced cells
        do i = 1, hmesh%length
            if (phis(i) + EPS < 1.0D0) then
                cidx = cidx + 1
                fmesh%cells(cidx) = hmesh%cells(i)
                if (phis(i) > EPS) then
                    widx = widx + 1
                    fmesh%wall_indices(cidx) = widx
                    fmesh%solid_fractions(widx) = phis(i)
                    fmesh%wall_normals(widx, :) = 0 ! TODO
                end if
            end if
        end do
        ! allocate cells
        ! fill reduced_mesh
        ! calculate normals as well, maybe with phi or maybe after

    end subroutine

    subroutine line_cell_phi(cel, line_point, line_dir, phi)
        type(cell2d), intent(in) :: cel
        real(8), intent(in) :: line_point(2), line_dir(2)
        real(8), intent(out) :: phi

        !TODO

        phi = 0.0d0
    end subroutine line_cell_phi

    subroutine poly_cell_phi(cel, poly, phi)
        type(cell2d), intent(in) :: cel
        real(8), intent(in) :: poly(:,:)
        real(8), intent(out) :: phi

        !TODO

        phi = 0.0d0
    end subroutine poly_cell_phi
    
end module mesh_gen
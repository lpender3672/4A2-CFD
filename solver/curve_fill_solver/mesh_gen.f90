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
        write(unit,'(A)') 'i,xmin,xmax,ymin,ymax,curvature,distance'
        do i = 1, mesh%length
            write(unit,'(I6,1X,6ES20.10)') i, mesh%cells(i)%xmin, mesh%cells(i)%xmax, &
                                            mesh%cells(i)%ymin, mesh%cells(i)%ymax, &
                                            mesh%nearest_cell_poly_curvature(i), &
                                            mesh%nearest_cell_poly_distance(i)
        end do
        close(unit)
    end subroutine write_mesh_csv

    subroutine generate_cmesh(n, m, naca_code, chord, aoa, fmesh)
        implicit none
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: naca_code
        real(8), intent(in) :: chord, aoa

        type(helper_lod_mesh) :: hmesh
        type(lod_mesh), intent(out) :: fmesh

        real(8), allocatable :: foil(:,:), foil_t(:,:)

        call naca4_airfoil("2412", 400, .true., foil)
        allocate(foil_t(size(foil,1),2))
        call transform_airfoil(foil, chord, aoa, origin=[real(n*0.5, 8), real(m*0.5, 8)], poly_out=foil_t)

        print *, 'Airfoil generated with ', size(foil_t,1), ' points.'

        call alloc_ncells(n, m, foil_t, 1, hmesh)
        call build_cells(n, m, foil_t, 1, hmesh)

        print *, 'Full mesh cells allocated: ', hmesh%length
        call write_mesh_csv(hmesh, 'mesh.csv')

        call build_ghost_cells(hmesh, foil_t, fmesh)

        print *, 'Mesh cells: ', fmesh%length, '  ghost cells: ', fmesh%ghost_count

        call build_indicies(fmesh)

        print *, 'Total neighbour count ', size(fmesh%neigh_indices)

        call print_first_five_neighbors(fmesh)

    end subroutine generate_cmesh

    ! ------------------------------------------------------------------
    ! Populate fmesh from hmesh:
    !   - copy all cells, mark each as solid (centre inside poly) or fluid
    !   - identify ghost cells: fluid cells whose centre is within one
    !     cell-diagonal of the airfoil surface
    !   - for each ghost cell store the outward wall normal and the
    !     mirror-point coordinates used by the solver to set the BC state
    ! ------------------------------------------------------------------
    subroutine build_ghost_cells(hmesh, poly, fmesh)
        type(helper_lod_mesh), intent(in)  :: hmesh
        real(8),               intent(in)  :: poly(:,:)
        type(lod_mesh),        intent(out) :: fmesh

        integer :: i, gcnt
        real(8) :: cx, cy, cell_diag, nx, ny, qx, qy, dn

        fmesh%length = hmesh%length
        allocate(fmesh%cells(fmesh%length))
        fmesh%cells = hmesh%cells

        ! --- mark solid / fluid -------------------------------------------
        do i = 1, fmesh%length
            cx = 0.5D0*(fmesh%cells(i)%xmin + fmesh%cells(i)%xmax)
            cy = 0.5D0*(fmesh%cells(i)%ymin + fmesh%cells(i)%ymax)
            if (point_in_polygon(cx, cy, poly)) then
                fmesh%cells(i)%is_solid = 1
            else
                fmesh%cells(i)%is_solid = 0
            end if
        end do

        ! --- count ghost cells --------------------------------------------
        ! A fluid cell is a ghost cell when its centre is closer to the
        ! airfoil than the cell diagonal (i.e. the surface cuts through
        ! the cell neighbourhood).
        gcnt = 0
        do i = 1, fmesh%length
            if (fmesh%cells(i)%is_solid == 1) cycle
            cell_diag = sqrt((fmesh%cells(i)%xmax - fmesh%cells(i)%xmin)**2 + &
                             (fmesh%cells(i)%ymax - fmesh%cells(i)%ymin)**2)
            if (hmesh%nearest_cell_poly_distance(i) < cell_diag) gcnt = gcnt + 1
        end do

        fmesh%ghost_count = gcnt
        allocate(fmesh%ghost_indices(gcnt))
        allocate(fmesh%ghost_normals(gcnt, 2))
        allocate(fmesh%ghost_mirror (gcnt, 2))

        ! --- fill ghost cell data -----------------------------------------
        gcnt = 0
        do i = 1, fmesh%length
            if (fmesh%cells(i)%is_solid == 1) cycle
            cell_diag = sqrt((fmesh%cells(i)%xmax - fmesh%cells(i)%xmin)**2 + &
                             (fmesh%cells(i)%ymax - fmesh%cells(i)%ymin)**2)
            if (hmesh%nearest_cell_poly_distance(i) >= cell_diag) cycle

            gcnt = gcnt + 1
            cx = 0.5D0*(fmesh%cells(i)%xmin + fmesh%cells(i)%xmax)
            cy = 0.5D0*(fmesh%cells(i)%ymin + fmesh%cells(i)%ymax)

            call nearest_edge_normal(cx, cy, poly, qx, qy, nx, ny)

            ! signed normal distance from wall (positive into fluid)
            dn = (cx - qx)*nx + (cy - qy)*ny

            fmesh%ghost_indices(gcnt)   = i
            fmesh%ghost_normals(gcnt,1) = nx
            fmesh%ghost_normals(gcnt,2) = ny
            fmesh%ghost_mirror (gcnt,1) = cx - 2.0D0*dn*nx
            fmesh%ghost_mirror (gcnt,2) = cy - 2.0D0*dn*ny
        end do

    end subroutine build_ghost_cells

    ! ------------------------------------------------------------------
    ! Ray-casting point-in-polygon test.
    ! ------------------------------------------------------------------
    logical function point_in_polygon(px, py, poly)
        real(8), intent(in) :: px, py, poly(:,:)
        integer :: i, j, n
        n = size(poly, 1)
        point_in_polygon = .false.
        j = n
        do i = 1, n
            if (((poly(i,2) > py) .neqv. (poly(j,2) > py)) .and. &
                (px < (poly(j,1)-poly(i,1))*(py-poly(i,2)) / &
                      (poly(j,2)-poly(i,2)) + poly(i,1))) then
                point_in_polygon = .not. point_in_polygon
            end if
            j = i
        end do
    end function point_in_polygon

    ! ------------------------------------------------------------------
    ! Find the closest point (qx,qy) on the polygon edges to (px,py),
    ! and return the outward unit normal (nx,ny) at that edge.
    ! "Outward" = pointing toward (px,py), i.e. into the fluid.
    ! ------------------------------------------------------------------
    subroutine nearest_edge_normal(px, py, poly, qx, qy, nx, ny)
        real(8), intent(in)  :: px, py, poly(:,:)
        real(8), intent(out) :: qx, qy, nx, ny

        integer :: i, j, n
        real(8) :: dx, dy, len2, t, ex, ey, d2, best_d2, inv_len

        n = size(poly, 1)
        best_d2 = 1.0D300
        qx = poly(1,1);  qy = poly(1,2)
        nx = 0.0D0;      ny = 1.0D0

        do i = 1, n
            j = mod(i, n) + 1
            dx = poly(j,1) - poly(i,1)
            dy = poly(j,2) - poly(i,2)
            len2 = dx*dx + dy*dy
            if (len2 < 1.0D-24) cycle

            t  = ((px - poly(i,1))*dx + (py - poly(i,2))*dy) / len2
            t  = max(0.0D0, min(1.0D0, t))
            ex = poly(i,1) + t*dx
            ey = poly(i,2) + t*dy
            d2 = (px-ex)**2 + (py-ey)**2

            if (d2 < best_d2) then
                best_d2 = d2
                qx = ex;  qy = ey
                inv_len = 1.0D0 / sqrt(len2)
                ! Two candidate outward normals perpendicular to edge tangent.
                ! Pick the one with a positive dot product toward (px,py).
                if ((-dy*inv_len)*(px-qx) + (dx*inv_len)*(py-qy) > 0.0D0) then
                    nx = -dy*inv_len;  ny =  dx*inv_len
                else
                    nx =  dy*inv_len;  ny = -dx*inv_len
                end if
            end if
        end do

    end subroutine nearest_edge_normal

end module mesh_gen
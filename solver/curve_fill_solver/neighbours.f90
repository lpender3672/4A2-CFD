
module neighbouring
  use iso_fortran_env, only: real64
  use types
  implicit none

  contains

  integer pure function cells_share_edge(a, b, eps)

    type(cell2d), intent(in) :: a, b
    real(c_double), intent(in) :: eps

    real(c_double) :: x1a, x1b, y1a, y1b
    real(c_double) :: x2a, x2b, y2a, y2b
    logical :: horiz_touch, horiz_overlap, vert_touch, vert_overlap

    cells_share_edge = 0

    ! Bounds from your cell2d structure
    x1a = a%xmin;  x1b = a%xmax
    y1a = a%ymin;  y1b = a%ymax
    x2a = b%xmin;  x2b = b%xmax
    y2a = b%ymin;  y2b = b%ymax

    ! Horizontal edge contact (cells stacked vertically)
    horiz_touch   = (abs(y1a - y2b) < eps) .or. (abs(y1b - y2a) < eps)
    horiz_overlap = min(x1b, x2b) > max(x1a, x2a) + eps

    ! Vertical edge contact (cells side by side)
    vert_touch    = (abs(x1a - x2b) < eps) .or. (abs(x1b - x2a) < eps)
    vert_overlap  = min(y1b, y2b) > max(y1a, y2a) + eps

    if (vert_touch .and. vert_overlap .and. x1a > x2a) cells_share_edge = 1  ! left
    if (vert_touch .and. vert_overlap .and. x1a < x2a) cells_share_edge = 2  ! right
    if (horiz_touch .and. horiz_overlap .and. y1a > y2a) cells_share_edge = 3 ! bottom
    if (horiz_touch .and. horiz_overlap .and. y1a < y2a) cells_share_edge = 4 ! top
    
  end function cells_share_edge

  subroutine build_indicies(mesh)

    type(lod_mesh), intent(inout) :: mesh

    real(c_double) :: tol
    integer :: n_cells, i, j, s, total_neigh, pos
    integer :: side, tmp_count
    ! tempory buffers to store neighbour indicies before writing them in side order
    ! necessary because looping through a curve could give neighbours at any side
    integer, parameter :: MAX_SIDE_NEIGH = 32
    integer :: tmp_sides(4, MAX_SIDE_NEIGH)
    integer :: cached_side(mesh%length, mesh%length)

    tol = 1d-9
    n_cells = mesh%length
    tmp_sides = 0
    cached_side = 0

    ! count neighbours to determine size to allocate
    do i = 1, n_cells
        tmp_count = 0

        do j = 1, n_cells
            if (j == i) cycle

            side = cells_share_edge(mesh%cells(i), mesh%cells(j), tol)
            if (side > 0) tmp_count = tmp_count + 1
            cached_side(i, j) = side
        end do

        mesh%cells(i)%neigh_count = tmp_count
    end do

    mesh%cells(1)%neigh_offset = 1
    do i = 2, n_cells
        mesh%cells(i)%neigh_offset = mesh%cells(i-1)%neigh_offset + mesh%cells(i-1)%neigh_count
    end do
    total_neigh = mesh%cells(n_cells)%neigh_offset + mesh%cells(n_cells)%neigh_count - 1
    
    allocate(mesh%neigh_indices(total_neigh))
    mesh%neigh_indices = 0

    ! fill 
    do i = 1, n_cells
        ! zero side counts
        mesh%cells(i)%side_count = 0

        ! right now we need to sort the neighbours
        do j = 1, n_cells
            if (j == i) cycle

            side = cached_side(i, j)
            if (side > 0) then
                mesh%cells(i)%side_count(side) = mesh%cells(i)%side_count(side) + 1
                tmp_sides(side, mesh%cells(i)%side_count(side)) = j
            end if
        end do

        tmp_count   = 0
        pos         = mesh%cells(i)%neigh_offset

        do s = 1, 4
            do j = 1, mesh%cells(i)%side_count(s)
                mesh%neigh_indices(pos + tmp_count) = tmp_sides(s, j);
                tmp_count = tmp_count + 1;
            end do
        end do
    end do

  end subroutine build_indicies

  subroutine print_first_five_neighbors(mesh)
    use iso_c_binding
    implicit none
    type(lod_mesh), intent(in) :: mesh

    integer :: i, start, count, last, n_show

    if (.not.allocated(mesh%cells)) then
        print *, "Cells not allocated"
        return
    endif
    if (.not.allocated(mesh%neigh_indices)) then
        print *, "Neighbor indices not allocated"
        return
    endif

    n_show = min(mesh%length, 5)

    print *, "---- First ", n_show, " cells ----"
    do i = 1, n_show
        count = mesh%cells(i)%neigh_count
        start = mesh%cells(i)%neigh_offset
        if (count > 0) then
            last = start + count - 1
            print "(A,I3,A,*(I5))", "Cell ", i, " neighbors:", mesh%neigh_indices(start:last)
        else
            print "(A,I3,A)", "Cell ", i, " has no neighbors."
        endif
    end do
    print *, "---------------------------"
  end subroutine print_first_five_neighbors

end module neighbouring
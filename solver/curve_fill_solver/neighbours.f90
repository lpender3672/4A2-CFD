
module neighbouring
  use iso_fortran_env, only: real64
  use types
  implicit none

  contains

  logical function cells_share_edge(a, b, eps, side)

    type(cell2d), intent(in) :: a, b
    real(c_double), intent(in) :: eps
    integer, intent(out) :: side

    real(c_double) :: x1a, x1b, y1a, y1b
    real(c_double) :: x2a, x2b, y2a, y2b
    logical :: horiz_touch, horiz_overlap, vert_touch, vert_overlap

    ! Bounds from your cell2d structure
    x1a = a%xmin;  x1b = a%xmax
    y1a = a%ymin;  y1b = a%ymax
    x2a = b%xmin;  x2b = b%xmax
    y2a = b%ymin;  y2b = b%ymax

    ! Horizontal edge contact (cells stacked vertically)
    horiz_touch   = (abs(y1a - y2b) < eps) .or. (abs(y1b - y2a) < eps)
    horiz_overlap = (x1b > x2a + eps) .and. (x2b > x1a + eps)

    ! Vertical edge contact (cells side by side)
    vert_touch    = (abs(x1a - x2b) < eps) .or. (abs(x1b - x2a) < eps)
    vert_overlap  = (y1b > y2a + eps) .and. (y2b > y1a + eps)

    cells_share_edge = (horiz_touch .and. horiz_overlap) .or. &
                       (vert_touch  .and. vert_overlap)

    if (vert_touch .and. x1a > x2a) side = 1  ! left
    if (vert_touch .and. x1a < x2a) side = 2  ! right
    if (horiz_touch .and. y1a > y2a) side = 3 ! bottom
    if (horiz_touch .and. y1a < y2a) side = 4 ! top
    
  end function cells_share_edge

  subroutine build_indicies(mesh)

    type(lod_mesh), intent(inout) :: mesh

    real(c_double) :: tol
    integer :: n_cells, i, j, k, total_neigh, pos
    integer, parameter :: MAX_TEMP_NEIGH = 64
    integer :: tmp_neigh(MAX_TEMP_NEIGH)
    integer :: tmp_count
    real(c_double) :: x0a, x1a, y0a, y1a
    real(c_double) :: x0b, x1b, y0b, y1b
    real(c_double) :: dx, dy

    integer, allocatable :: neigh_count(:)
    integer, allocatable :: neigh_offset(:)
    integer :: side

    tol = 1d-6
    n_cells = mesh%length

    allocate(neigh_count(n_cells))
    neigh_count = 0

    ! count neighbours to determine size to allocate
    do i = 1, n_cells
        tmp_count = 0

        do j = 1, n_cells
            if (j == i) cycle

            if (cells_share_edge(mesh%cells(i), mesh%cells(j), tol, side)) then
                tmp_count = tmp_count + 1
            end if
        end do

        neigh_count(i) = tmp_count
    end do

    allocate(neigh_offset(n_cells))
    neigh_offset(1) = 1
    do i = 2, n_cells
        neigh_offset(i) = neigh_offset(i-1) + neigh_count(i-1)
    end do
    total_neigh = neigh_offset(n_cells) + neigh_count(n_cells) - 1
    allocate(mesh%neigh_indices(total_neigh))

    ! fill 
    do i = 1, n_cells
        tmp_count = 0

        ! zero side counts
        do j = 1, 4
            mesh%cells(i)%side_count(j) = 0
        end do

        do j = 1, n_cells
            if (j == i) cycle

            if (cells_share_edge(mesh%cells(i), mesh%cells(j), tol, side)) then
                tmp_count = tmp_count + 1
                pos = neigh_offset(i) + tmp_count - 1
                mesh%neigh_indices(pos) = j

                mesh%cells(i)%side_count(side) = mesh%cells(i)%side_count(side) + 1
            end if
        end do

        mesh%cells(i)%neigh_offset = neigh_offset(i)
        mesh%cells(i)%neigh_count  = int(tmp_count, kind=c_int8_t)

    end do

    deallocate(neigh_offset, neigh_count)

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

    n_show = min(mesh%length, 10)

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
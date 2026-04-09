
module neighbouring
  use iso_fortran_env, only: real64
  use types
  implicit none

  contains

  ! -----------------------------------------------------------------------
  ! Quicksort argsort helpers — sort index array by cell coordinate field
  ! -----------------------------------------------------------------------

  recursive subroutine qsort_xmin(cells, idx, lo, hi)
    type(cell2d), intent(in)    :: cells(:)
    integer,      intent(inout) :: idx(:)
    integer,      intent(in)    :: lo, hi
    integer       :: i, j, tmp
    real(c_double):: piv
    if (lo >= hi) return
    piv = cells(idx((lo+hi)/2))%xmin
    i = lo;  j = hi
    do
      do while (cells(idx(i))%xmin < piv);  i = i + 1;  end do
      do while (cells(idx(j))%xmin > piv);  j = j - 1;  end do
      if (i >= j) exit
      tmp = idx(i);  idx(i) = idx(j);  idx(j) = tmp
      i = i + 1;     j = j - 1
    end do
    call qsort_xmin(cells, idx, lo, j)
    call qsort_xmin(cells, idx, j+1, hi)
  end subroutine

  recursive subroutine qsort_ymin(cells, idx, lo, hi)
    type(cell2d), intent(in)    :: cells(:)
    integer,      intent(inout) :: idx(:)
    integer,      intent(in)    :: lo, hi
    integer       :: i, j, tmp
    real(c_double):: piv
    if (lo >= hi) return
    piv = cells(idx((lo+hi)/2))%ymin
    i = lo;  j = hi
    do
      do while (cells(idx(i))%ymin < piv);  i = i + 1;  end do
      do while (cells(idx(j))%ymin > piv);  j = j - 1;  end do
      if (i >= j) exit
      tmp = idx(i);  idx(i) = idx(j);  idx(j) = tmp
      i = i + 1;     j = j - 1
    end do
    call qsort_ymin(cells, idx, lo, j)
    call qsort_ymin(cells, idx, j+1, hi)
  end subroutine

  ! -----------------------------------------------------------------------
  ! Binary search: first position in sorted idx where xmin >= target - tol
  ! -----------------------------------------------------------------------

  integer function bs_xmin(cells, idx, n, target, tol) result(lo)
    type(cell2d),  intent(in) :: cells(:)
    integer,       intent(in) :: idx(:), n
    real(c_double),intent(in) :: target, tol
    integer :: hi, mid
    lo = 1;  hi = n + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (cells(idx(mid))%xmin < target - tol) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
  end function

  integer function bs_ymin(cells, idx, n, target, tol) result(lo)
    type(cell2d),  intent(in) :: cells(:)
    integer,       intent(in) :: idx(:), n
    real(c_double),intent(in) :: target, tol
    integer :: hi, mid
    lo = 1;  hi = n + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (cells(idx(mid))%ymin < target - tol) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
  end function

  ! -----------------------------------------------------------------------
  ! Main routine
  ! -----------------------------------------------------------------------

  subroutine build_indicies(mesh)

    type(lod_mesh), intent(inout) :: mesh

    ! MAX_NEIGH: max neighbours any single cell can accumulate across both passes.
    ! In a Hilbert AMR mesh each cell touches at most O(2^level_diff) fine cells
    ! per face; 128 is generous even for aggressive refinement ratios.
    integer, parameter :: MAX_NEIGH = 128

    real(c_double) :: tol
    integer :: n_cells, i, j, k, s, lo, pos, total_neigh

    integer, allocatable :: srt_x(:), srt_y(:)   ! argsort indices
    ! Per-cell temporary neighbour lists (neighbour cell index + which side)
    integer, allocatable :: tmp_j(:,:), tmp_s(:,:), tmp_cnt(:)

    tol    = 1d-9
    n_cells = mesh%length

    allocate(srt_x(n_cells), srt_y(n_cells))
    allocate(tmp_j  (MAX_NEIGH, n_cells), source=0)
    allocate(tmp_s  (MAX_NEIGH, n_cells), source=0)
    allocate(tmp_cnt(n_cells),            source=0)

    ! Build sorted index arrays
    do i = 1, n_cells
      srt_x(i) = i
      srt_y(i) = i
    end do
    call qsort_xmin(mesh%cells, srt_x, 1, n_cells)
    call qsort_ymin(mesh%cells, srt_y, 1, n_cells)

    ! -----------------------------------------------------------------------
    ! Right / left contacts: cell j is to the right of cell i when
    !   xmin[j] == xmax[i]   (within tol)   AND   y-intervals overlap.
    ! We find the contacts for all i in one sorted sweep; each found (i,j)
    ! pair is registered from both sides immediately.
    ! -----------------------------------------------------------------------
    do i = 1, n_cells
      lo = bs_xmin(mesh%cells, srt_x, n_cells, mesh%cells(i)%xmax, tol)
      do k = lo, n_cells
        j = srt_x(k)
        if (mesh%cells(j)%xmin > mesh%cells(i)%xmax + tol) exit  ! past window
        if (j == i) cycle
        if (min(mesh%cells(i)%ymax, mesh%cells(j)%ymax) > &
            max(mesh%cells(i)%ymin, mesh%cells(j)%ymin) + tol) then
          ! j is right neighbour of i  (side 2)
          tmp_cnt(i) = tmp_cnt(i) + 1
          tmp_j(tmp_cnt(i), i) = j;  tmp_s(tmp_cnt(i), i) = 2
          ! i is left  neighbour of j  (side 1)
          tmp_cnt(j) = tmp_cnt(j) + 1
          tmp_j(tmp_cnt(j), j) = i;  tmp_s(tmp_cnt(j), j) = 1
        end if
      end do
    end do

    ! -----------------------------------------------------------------------
    ! Top / bottom contacts: cell j is above cell i when
    !   ymin[j] == ymax[i]   AND   x-intervals overlap.
    ! -----------------------------------------------------------------------
    do i = 1, n_cells
      lo = bs_ymin(mesh%cells, srt_y, n_cells, mesh%cells(i)%ymax, tol)
      do k = lo, n_cells
        j = srt_y(k)
        if (mesh%cells(j)%ymin > mesh%cells(i)%ymax + tol) exit
        if (j == i) cycle
        if (min(mesh%cells(i)%xmax, mesh%cells(j)%xmax) > &
            max(mesh%cells(i)%xmin, mesh%cells(j)%xmin) + tol) then
          ! j is top    neighbour of i  (side 4)
          tmp_cnt(i) = tmp_cnt(i) + 1
          tmp_j(tmp_cnt(i), i) = j;  tmp_s(tmp_cnt(i), i) = 4
          ! i is bottom neighbour of j  (side 3)
          tmp_cnt(j) = tmp_cnt(j) + 1
          tmp_j(tmp_cnt(j), j) = i;  tmp_s(tmp_cnt(j), j) = 3
        end if
      end do
    end do

    ! -----------------------------------------------------------------------
    ! Populate cell metadata from temp buffers
    ! -----------------------------------------------------------------------
    do i = 1, n_cells
      mesh%cells(i)%neigh_count = tmp_cnt(i)
      mesh%cells(i)%side_count  = 0
      do k = 1, tmp_cnt(i)
        s = tmp_s(k, i)
        mesh%cells(i)%side_count(s) = mesh%cells(i)%side_count(s) + 1
      end do
    end do

    ! CSR offsets
    mesh%cells(1)%neigh_offset = 1
    do i = 2, n_cells
      mesh%cells(i)%neigh_offset = mesh%cells(i-1)%neigh_offset + mesh%cells(i-1)%neigh_count
    end do
    total_neigh = mesh%cells(n_cells)%neigh_offset + mesh%cells(n_cells)%neigh_count - 1

    allocate(mesh%neigh_indices(total_neigh))
    mesh%neigh_indices = 0

    ! Write neigh_indices in side order 1→2→3→4 so the solver can slice
    ! the CSR row by side using side_count.
    do i = 1, n_cells
      pos = mesh%cells(i)%neigh_offset
      k   = 0
      do s = 1, 4
        do j = 1, tmp_cnt(i)
          if (tmp_s(j, i) == s) then
            mesh%neigh_indices(pos + k) = tmp_j(j, i)
            k = k + 1
          end if
        end do
      end do
    end do

    deallocate(srt_x, srt_y, tmp_j, tmp_s, tmp_cnt)

  end subroutine build_indicies

  ! -----------------------------------------------------------------------

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


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

  ! -----------------------------------------------------------------------
  ! Mesh quality for axis-aligned AMR meshes.  Computes per-cell:
  !   ar       — aspect ratio max(dx/dy, dy/dx); should be 1 on Cartesian AMR
  !   ratio    — max size ratio with any neighbour (= 2:1-balance metric;
  !              well-balanced AMR keeps this ≤ 2)
  !   nsolid   — count of solid neighbours.  ≥ 2 means a ghost cell sitting
  !              at a concave stair-step corner (where the wall turns 90°);
  !              these dominate pressure ringing along curved surfaces.
  !   nfluid   — count of fluid neighbours.  0 → orphan, unsolvable.
  ! Reports summary stats to stdout and writes a CSV for postprocessing.
  ! -----------------------------------------------------------------------
  subroutine mesh_quality_report(mesh, csv_path)
    type(lod_mesh),   intent(in) :: mesh
    character(len=*), intent(in) :: csv_path

    integer :: i, k, j, offset
    real(8) :: dx_i, dy_i, dx_j, dy_j, ar_i, r, ratio_i
    real(8) :: ar_max, ar_max_loc(2), area_min, area_max
    real(8) :: ratio_max, ratio_max_loc(2)
    integer :: n_fluid, n_solid, n_ghost, n_orphan
    integer :: n_jump_gt2, n_jump_gt4
    integer :: n_concave        ! ghost cells with ≥ 2 solid sides
    integer :: nf_i, ns_i, ar_max_idx, ratio_max_idx, unit
    real(8) :: cx, cy

    real(8), allocatable :: cell_ratio(:), cell_ar(:)
    integer, allocatable :: cell_nf(:), cell_ns(:)

    allocate(cell_ratio(mesh%length), cell_ar(mesh%length))
    allocate(cell_nf(mesh%length), cell_ns(mesh%length))

    ar_max     = 0.0D0;  ar_max_idx     = 1
    ratio_max  = 1.0D0;  ratio_max_idx  = 1
    area_min   = huge(1.0D0)
    area_max   = 0.0D0
    n_fluid    = 0;  n_solid = 0;  n_ghost = 0;  n_orphan = 0
    n_jump_gt2 = 0;  n_jump_gt4 = 0;  n_concave = 0

    do i = 1, mesh%length
      dx_i = mesh%cells(i)%xmax - mesh%cells(i)%xmin
      dy_i = mesh%cells(i)%ymax - mesh%cells(i)%ymin
      ar_i = max(dx_i/dy_i, dy_i/dx_i)
      cell_ar(i) = ar_i
      area_min   = min(area_min, dx_i*dy_i)
      area_max   = max(area_max, dx_i*dy_i)
      if (ar_i > ar_max) then
        ar_max = ar_i;  ar_max_idx = i
      end if

      if (mesh%cells(i)%is_solid == 1) then
        n_solid = n_solid + 1
        cell_ratio(i) = 1.0D0;  cell_nf(i) = 0;  cell_ns(i) = 0
        cycle
      end if
      n_fluid = n_fluid + 1

      offset  = mesh%cells(i)%neigh_offset
      ratio_i = 1.0D0
      nf_i    = 0;  ns_i = 0
      do k = 0, mesh%cells(i)%neigh_count - 1
        j = mesh%neigh_indices(offset + k)
        if (mesh%cells(j)%is_solid == 1) then
          ns_i = ns_i + 1
          cycle
        end if
        nf_i = nf_i + 1
        dx_j = mesh%cells(j)%xmax - mesh%cells(j)%xmin
        dy_j = mesh%cells(j)%ymax - mesh%cells(j)%ymin
        r = max(dx_i/dx_j, dx_j/dx_i, dy_i/dy_j, dy_j/dy_i)
        if (r > ratio_i) ratio_i = r
      end do
      cell_ratio(i) = ratio_i
      cell_nf   (i) = nf_i
      cell_ns   (i) = ns_i

      if (ns_i >= 1)            n_ghost  = n_ghost  + 1
      if (ns_i >= 2)            n_concave = n_concave + 1
      if (nf_i == 0)            n_orphan = n_orphan + 1
      if (ratio_i > 2.0D0 + 1.0D-9) n_jump_gt2 = n_jump_gt2 + 1
      if (ratio_i > 4.0D0 + 1.0D-9) n_jump_gt4 = n_jump_gt4 + 1
      if (ratio_i > ratio_max) then
        ratio_max = ratio_i;  ratio_max_idx = i
      end if
    end do

    ar_max_loc = [0.5D0*(mesh%cells(ar_max_idx   )%xmin + mesh%cells(ar_max_idx   )%xmax), &
                  0.5D0*(mesh%cells(ar_max_idx   )%ymin + mesh%cells(ar_max_idx   )%ymax)]
    ratio_max_loc = [0.5D0*(mesh%cells(ratio_max_idx)%xmin + mesh%cells(ratio_max_idx)%xmax), &
                     0.5D0*(mesh%cells(ratio_max_idx)%ymin + mesh%cells(ratio_max_idx)%ymax)]

    print *, '---- mesh quality ----'
    print '(A,I7,A,I7,A,I7,A,I7)', '  cells: total=', mesh%length, &
        '  fluid=', n_fluid, '  solid=', n_solid, '  ghost(wall-adj)=', n_ghost
    print '(A,ES10.3,A,ES10.3,A,ES10.3)', '  area: min=', area_min, &
        '  max=', area_max, '  max/min=', area_max/max(area_min, 1.0D-30)
    print '(A,F6.3,A,I6,A,F6.3,A,F6.3,A)', '  aspect ratio max=', ar_max, &
        '  at cell ', ar_max_idx, ' @(', ar_max_loc(1), ',', ar_max_loc(2), ')'
    print '(A,F6.3,A,I6,A,F6.3,A,F6.3,A)', '  AMR size-jump max=', ratio_max, &
        '  at cell ', ratio_max_idx, ' @(', ratio_max_loc(1), ',', ratio_max_loc(2), ')'
    print '(A,I6,A,I6)', '  cells with jump>2: ', n_jump_gt2, &
        '   cells with jump>4: ', n_jump_gt4
    print '(A,I6,A,I6)', '  concave stair corners (ghost ns>=2): ', n_concave, &
        '   orphan fluid cells (nf=0): ', n_orphan
    print *, '----------------------'

    ! Per-cell CSV — colour-plot any column over the mesh to localise issues
    open(newunit=unit, file=csv_path, status='replace', action='write', form='formatted')
    write(unit,'(A)') 'i,xmin,xmax,ymin,ymax,is_solid,nfluid,nsolid,aspect_ratio,size_ratio'
    do i = 1, mesh%length
      cx = mesh%cells(i)%xmin;  cy = mesh%cells(i)%ymin
      write(unit,'(I7,1X,4ES20.10,3(1X,I3),2(1X,ES20.10))') &
        i, mesh%cells(i)%xmin, mesh%cells(i)%xmax, &
        mesh%cells(i)%ymin, mesh%cells(i)%ymax, &
        mesh%cells(i)%is_solid, cell_nf(i), cell_ns(i), &
        cell_ar(i), cell_ratio(i)
    end do
    close(unit)

    deallocate(cell_ratio, cell_ar, cell_nf, cell_ns)
  end subroutine mesh_quality_report

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

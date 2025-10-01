
module sfc_quadtree_airfoil
  use types
  implicit none

contains

  !======================== morton/z-order basics ========================

  integer(i4) function quantize(v, vmin, vmax, bits) result(q)
    real(rk), intent(in) :: v, vmin, vmax
    integer,  intent(in) :: bits
    real(rk) :: t
    if (vmax <= vmin) then
       q = 0; return
    end if
    t = (v - vmin) / (vmax - vmin)
    if (t < 0.0_rk) t = 0.0_rk
    if (t > 1.0_rk) t = 1.0_rk
    q = int( t * real(2**bits - 1, rk), i4)
  end function quantize

  integer(i8) function morton2(ix, iy, bits) result(k)
    integer(i4), intent(in) :: ix, iy
    integer,    intent(in)  :: bits
    integer :: b
    integer(i8) :: bit
    k = 0_i8
    do b = 0, bits-1
       bit = merge(1_i8, 0_i8, btest(ix,b))
       k = ior(k, ishft(bit, 2*b+1))
       bit = merge(1_i8, 0_i8, btest(iy,b))
       k = ior(k, ishft(bit, 2*b  ))
    end do
  end function morton2

  !======================== small geometry helpers ========================

  logical function segs_intersect(ax,ay,bx,by, cx,cy,dx,dy) result(hit)
    real(rk), intent(in) :: ax,ay,bx,by,cx,cy,dx,dy
    real(rk) :: d, s, t
    real(rk) :: r_px, r_py, s_px, s_py, q_px, q_py
    r_px = bx - ax; r_py = by - ay
    s_px = dx - cx; s_py = dy - cy
    q_px = cx - ax; q_py = cy - ay
    d = r_px*s_py - r_py*s_px
    if (abs(d) < 1.0e-15_rk) then
       hit = .false.; return   ! parallel or colinear (treat as no-cross)
    end if
    s = (q_px*s_py - q_py*s_px) / d
    t = (q_px*r_py - q_py*r_px) / d
    hit = (s >= 0.0_rk .and. s <= 1.0_rk .and. t >= 0.0_rk .and. t <= 1.0_rk)
  end function segs_intersect

  logical function point_in_box(x,y, xmin,xmax,ymin,ymax) result(inside)
    real(rk), intent(in) :: x,y,xmin,xmax,ymin,ymax
    inside = (x >= xmin .and. x <= xmax .and. y >= ymin .and. y <= ymax)
  end function point_in_box

  logical function segment_hits_box(x1,y1,x2,y2, xmin,xmax,ymin,ymax) result(hit)
    real(rk), intent(in) :: x1,y1,x2,y2, xmin,xmax,ymin,ymax
    ! trivial accept: endpoint inside
    if (point_in_box(x1,y1,xmin,xmax,ymin,ymax)) then; hit = .true.; return; end if
    if (point_in_box(x2,y2,xmin,xmax,ymin,ymax)) then; hit = .true.; return; end if
    ! test against 4 box edges
    hit = .false.
    hit = hit .or. segs_intersect(x1,y1,x2,y2, xmin,ymin, xmin,ymax)
    hit = hit .or. segs_intersect(x1,y1,x2,y2, xmin,ymax, xmax,ymax)
    hit = hit .or. segs_intersect(x1,y1,x2,y2, xmax,ymax, xmax,ymin)
    hit = hit .or. segs_intersect(x1,y1,x2,y2, xmax,ymin, xmin,ymin)
  end function segment_hits_box

  logical function polyline_hits_box(xp,yp,n, xmin,xmax,ymin,ymax) result(hit)
    real(rk), intent(in) :: xp(n), yp(n)
    integer,  intent(in) :: n
    real(rk), intent(in) :: xmin,xmax,ymin,ymax
    integer :: i
    hit = .false.
    do i=1,n-1
       if (segment_hits_box(xp(i),yp(i), xp(i+1),yp(i+1), xmin,xmax,ymin,ymax)) then
          hit = .true.; return
       end if
    end do
  end function polyline_hits_box

  real(rk) function point_segment_distance(px,py, ax,ay, bx,by) result(d)
    real(rk), intent(in) :: px,py, ax,ay, bx,by
    real(rk) :: vx,vy, wx,wy, c1, c2, t, dx,dy, denom
    vx = bx - ax; vy = by - ay
    wx = px - ax; wy = py - ay
    denom = vx*vx + vy*vy
    if (denom <= 0.0_rk) then
       dx = px - ax; dy = py - ay; d = sqrt(dx*dx + dy*dy); return
    end if
    c1 = wx*vx + wy*vy
    if (c1 <= 0.0_rk) then
       dx = px - ax; dy = py - ay; d = sqrt(dx*dx + dy*dy); return
    end if
    c2 = denom
    if (c2 <= c1) then
      dx = px - bx; dy = py - by; d = sqrt(dx*dx + dy*dy); return
    end if
    t = c1 / c2
    dx = px - (ax + t*vx)
    dy = py - (ay + t*vy)
    d = sqrt(dx*dx + dy*dy)
  end function point_segment_distance

  real(rk) function distance_to_polyline(px,py, xp,yp,n) result(dmin)
    real(rk), intent(in) :: px,py, xp(n), yp(n)
    integer,  intent(in) :: n
    integer :: i
    real(rk) :: d
    dmin = huge(dmin)
    do i=1,n-1
      d = point_segment_distance(px,py, xp(i),yp(i), xp(i+1),yp(i+1))
      if (d < dmin) dmin = d
    end do
  end function distance_to_polyline

  !======================== sizing law & refinement ========================

  logical function need_refine(c, xp,yp,n, hmin, beta, max_level) result(refine)
    type(cell2d), intent(in) :: c
    real(rk),     intent(in) :: xp(n), yp(n), hmin, beta
    integer,      intent(in) :: n, max_level
    real(rk) :: xc, yc, h, d
    logical :: touches

    ! centroid and nominal size
    xc = 0.5_rk*(c%xmin + c%xmax)
    yc = 0.5_rk*(c%ymin + c%ymax)
    h  = max(c%xmax - c%xmin, c%ymax - c%ymin)

    ! fast refine if box touches the polyline (captures the boundary)
    touches = polyline_hits_box(xp,yp,n, c%xmin,c%xmax,c%ymin,c%ymax)

    ! distance-based target size
    d = distance_to_polyline(xc,yc, xp,yp,n)
    ! local target: linear growth away from airfoil, clamped by hmin
    ! beta ~ growth factor (e.g. 0.5 → fairly tight near wall)
    ! you can tune this to your solver requirements
    if (beta <= 0.0_rk) then
       refine = (touches .and. c%level < max_level)
       return
    end if

    ! target size at centroid
    ! h_target = max(hmin, beta * d)
    refine = ( (touches .or. (h > max(hmin, beta*d))) .and. (c%level < max_level) )
  end function need_refine

  !======================== quadtree builder ========================

  subroutine build_quadtree(xp,yp,n, domain_bb, max_level, hmin, beta, leaves, nleaf)
    real(rk), intent(in) :: xp(n), yp(n)
    integer,  intent(in) :: n
    real(rk), intent(in) :: domain_bb(4)
    integer,  intent(in) :: max_level
    real(rk), intent(in) :: hmin, beta
    type(cell2d), allocatable, intent(out) :: leaves(:)
    integer, intent(out) :: nleaf

    type(cell2d), allocatable :: stack(:)
    integer :: top, cap

    ! start clean
    if (allocated(leaves)) deallocate(leaves)
    nleaf = 0

    call init_stack(stack, cap, top)
    call push_box(stack, cap, top, domain_bb(1),domain_bb(2),domain_bb(3),domain_bb(4), 0)

    do while (top > 0)
      call process_one(xp, yp, n, max_level, hmin, beta, stack, cap, top, leaves, nleaf)
    end do

    ! ensure we return a valid but empty array if nothing was created (shouldn’t happen)
    if (.not. allocated(leaves)) then
      allocate(leaves(0))
      nleaf = 0
    end if
  end subroutine build_quadtree

  subroutine init_stack(stack, cap, top)
    type(cell2d), allocatable, intent(out) :: stack(:)
    integer, intent(out) :: cap, top
    cap = 1024
    allocate(stack(cap))
    top = 0
  end subroutine init_stack

  subroutine push_box(stack, cap, top, xmin,xmax,ymin,ymax, level)
    type(cell2d), allocatable, intent(inout) :: stack(:)
    integer,      intent(inout) :: cap, top
    real(rk), intent(in) :: xmin,xmax,ymin,ymax
    integer,  intent(in) :: level
    if (top >= cap) call grow(stack, cap)
    top = top + 1
    stack(top)%xmin = xmin; stack(top)%xmax = xmax
    stack(top)%ymin = ymin; stack(top)%ymax = ymax
    stack(top)%level = level
    stack(top)%key = 0_i8; stack(top)%id = 0
  end subroutine push_box

  subroutine grow(stack, cap)
    type(cell2d), allocatable, intent(inout) :: stack(:)
    integer, intent(inout) :: cap
    type(cell2d), allocatable :: tmp(:)
    integer :: newcap
    newcap = 2*cap
    allocate(tmp(newcap))
    tmp(1:cap) = stack(1:cap)
    deallocate(stack)
    call move_alloc(tmp, stack)
    cap = newcap
  end subroutine grow

  subroutine process_one(xp,yp,n, max_level, hmin, beta, stack, cap, top, leaves, nleaf)
    real(rk), intent(in) :: xp(n), yp(n)
    integer,  intent(in) :: n, max_level
    real(rk), intent(in) :: hmin, beta
    type(cell2d), allocatable, intent(inout) :: stack(:)
    type(cell2d), allocatable, intent(inout) :: leaves(:)
    integer, intent(inout) :: cap, top, nleaf
    type(cell2d) :: c
    real(rk) :: xmid, ymid

    c = stack(top); top = top - 1

    if (.not. need_refine(c, xp,yp,n, hmin, beta, max_level)) then
      call append_leaf(leaves, nleaf, c)
      return
    end if

    xmid = 0.5_rk*(c%xmin + c%xmax)
    ymid = 0.5_rk*(c%ymin + c%ymax)

    call push_box(stack, cap, top, c%xmin, xmid, c%ymin, ymid, c%level+1)
    call push_box(stack, cap, top, xmid,  c%xmax, c%ymin, ymid, c%level+1)
    call push_box(stack, cap, top, c%xmin, xmid, ymid,  c%ymax, c%level+1)
    call push_box(stack, cap, top, xmid,  c%xmax, ymid,  c%ymax, c%level+1)
  end subroutine process_one

  subroutine append_leaf(leaves, nleaf, c)
  type(cell2d), allocatable, intent(inout) :: leaves(:)
  integer, intent(inout) :: nleaf
  type(cell2d), intent(in) :: c
  type(cell2d), allocatable :: tmp(:)

  if (.not. allocated(leaves)) then
     allocate(leaves(1))
     leaves(1) = c
     nleaf = 1
  else
     allocate(tmp(nleaf+1))
     tmp(1:nleaf) = leaves
     tmp(nleaf+1) = c
     call move_alloc(tmp, leaves)
     nleaf = nleaf + 1
  end if
end subroutine append_leaf

  !======================== key building & sorting ========================

  subroutine build_keys_and_sort(leaves, nleaf, bits, bb4)
    type(cell2d), intent(inout) :: leaves(:)
    integer,      intent(in)    :: nleaf, bits
    real(rk),     intent(in)    :: bb4(4)
    integer :: i
    integer(i4) :: ix, iy
    real(rk) :: xc, yc

    do i=1,nleaf
      xc = 0.5_rk*(leaves(i)%xmin + leaves(i)%xmax)
      yc = 0.5_rk*(leaves(i)%ymin + leaves(i)%ymax)
      ix = quantize(xc, bb4(1), bb4(2), bits)
      iy = quantize(yc, bb4(3), bb4(4), bits)
      leaves(i)%key = morton2(ix,iy,bits)
      leaves(i)%id  = i
    end do

    call mergesort_cells(leaves, 1, nleaf)
  end subroutine build_keys_and_sort

  recursive subroutine mergesort_cells(a, l, r)
    type(cell2d), intent(inout) :: a(:)
    integer, intent(in) :: l, r
    integer :: m
    if (l >= r) return
    m = (l+r)/2
    call mergesort_cells(a, l, m)
    call mergesort_cells(a, m+1, r)
    call merge_cells(a, l, m, r)
  end subroutine mergesort_cells

  subroutine merge_cells(a, l, m, r)
    type(cell2d), intent(inout) :: a(:)
    integer, intent(in) :: l, m, r
    integer :: n1, n2, i, j, p
    type(cell2d), allocatable :: left(:), right(:)

    n1 = m - l + 1
    n2 = r - m
    allocate(left(n1), right(n2))
    left  = a(l:m)
    right = a(m+1:r)
    i=1; j=1; p=l
    do while (i<=n1 .and. j<=n2)
      if (left(i)%key <= right(j)%key) then
        a(p) = left(i); i=i+1
      else
        a(p) = right(j); j=j+1
      end if
      p=p+1
    end do
    do while (i<=n1); a(p)=left(i);  p=p+1; i=i+1; end do
    do while (j<=n2); a(p)=right(j); p=p+1; j=j+1; end do
    deallocate(left,right)
  end subroutine merge_cells

end module sfc_quadtree_airfoil

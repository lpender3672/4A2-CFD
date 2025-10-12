module sfc_reorder
   use types
   implicit none
   intrinsic :: merge, ieor

contains

   
  !==============================
  ! Coordinate normalization
  !==============================
  real(rk) function normalize(v, vmin, vmax) result(t)
    real(rk), intent(in) :: v, vmin, vmax
    if (vmax <= vmin) then
       t = 0.0_rk
       return
    end if
    t = (v - vmin) / (vmax - vmin)
    if (t < 0.0_rk) t = 0.0_rk
    if (t > 1.0_rk) t = 1.0_rk
  end function normalize

  integer(i4) function quantize_lod(t, bits) result(q)
    real(rk), intent(in) :: t
    integer,  intent(in) :: bits
    q = int(t * real(2**bits - 1, rk), i4)
  end function quantize_lod

  !==============================
  ! Morton (Z-order)
  !==============================
  integer(i8) function morton2(ix, iy, bits) result(key)
    integer(i4), intent(in) :: ix, iy
    integer,    intent(in)  :: bits
    integer :: i
    integer(i8) :: k
    key = 0_i8
    do i = 0, bits-1
       k = merge(1_i8, 0_i8, btest(ix,i))
       key = ior(key, ishft(k,2*i))
       k = merge(1_i8, 0_i8, btest(iy,i))
       key = ior(key, ishft(k,2*i+1))
    end do
  end function morton2

  integer(i8) function morton3(ix, iy, iz, bits) result(key)
    integer(i4), intent(in) :: ix, iy, iz
    integer,    intent(in)  :: bits
    integer :: i
    integer(i8) :: k
    key = 0_i8
    do i = 0, bits-1
       k = merge(1_i8, 0_i8, btest(ix,i)); key = ior(key, ishft(k,3*i    ))
       k = merge(1_i8, 0_i8, btest(iy,i)); key = ior(key, ishft(k,3*i+1))
       k = merge(1_i8, 0_i8, btest(iz,i)); key = ior(key, ishft(k,3*i+2))
    end do
  end function morton3

  !==============================
  ! Hilbert (2D)
  !==============================
  integer(i8) function hilbert2(ix, iy, bits) result(h)
    integer(i4), intent(in) :: ix, iy
    integer,    intent(in)  :: bits
    integer(i8) :: X(2), t, Q, P
    integer :: i, n

    X(1) = int(ix,i8)
    X(2) = int(iy,i8)
    n = 2
    Q = ishft(1_i8, bits-1)

    do while (Q > 1)
       P = Q - 1
       do i = 1, n
          if (iand(X(i), Q) /= 0) then
             X(1) = ieor(X(1), P)
          else
             t = iand(ieor(X(1), X(i)), P)
             X(1) = ieor(X(1), t)
             X(i) = ieor(X(i), t)
          end if
       end do
       Q = ishft(Q, -1)
    end do

    do i = n, 2, -1
       X(i) = ieor(X(i), X(i-1))
    end do
    t = 0_i8
    do Q = 2, ishft(1_i8, bits), Q*2
       if (iand(X(n), Q) /= 0) t = ieor(t, Q-1)
    end do
    do i = 1, n
       X(i) = ieor(X(i), t)
    end do

    h = 0_i8
    do i = bits-1, 0, -1
       do n = 1, 2
          h = ishft(h, 1)
          if (btest(X(n), i)) h = ior(h, 1_i8)
       end do
    end do
  end function hilbert2

  !==============================
  ! Hilbert (3D)
  !==============================
  integer(i8) function hilbert3(ix, iy, iz, bits) result(h)
    integer(i4), intent(in) :: ix, iy, iz
    integer,    intent(in)  :: bits
    integer(i8) :: X(3), t, Q, P
    integer :: i, n, j

    X(1) = int(ix,i8)
    X(2) = int(iy,i8)
    X(3) = int(iz,i8)
    n = 3
    Q = ishft(1_i8, bits-1)

    do while (Q > 1)
       P = Q - 1
       do i = 1, n
          if (iand(X(i), Q) /= 0) then
             X(1) = ieor(X(1), P)
          else
             t = iand(ieor(X(1), X(i)), P)
             X(1) = ieor(X(1), t)
             X(i) = ieor(X(i), t)
          end if
       end do
       Q = ishft(Q, -1)
    end do

    do i = n, 2, -1
       X(i) = ieor(X(i), X(i-1))
    end do
    t = 0_i8
    do Q = 2, ishft(1_i8, bits), Q*2
       if (iand(X(n), Q) /= 0) t = ieor(t, Q-1)
    end do
    do i = 1, n
       X(i) = ieor(X(i), t)
    end do

    h = 0_i8
    do i = bits-1, 0, -1
       do j = 1, n
          h = ishft(h, 1)
          if (btest(X(j), i)) h = ior(h, 1_i8)
       end do
    end do
  end function hilbert3

  !==============================
  ! Build keys at arbitrary LOD
  !==============================
  subroutine build_keys(mesh, bits, bb, use_hilbert, dim)
    type(cell), intent(inout) :: mesh(:)
    integer, intent(in)       :: bits
    real(rk), intent(in)      :: bb(:)      ! bounding box
    logical, intent(in), optional :: use_hilbert
    integer, intent(in), optional :: dim    ! 2 or 3D

    integer(i4) :: ix, iy, iz
    integer :: i
    logical :: hilb
    integer :: nd

    hilb = .false.; if (present(use_hilbert)) hilb = use_hilbert
    nd   = 3;       if (present(dim)) nd = dim

    do i=1,size(mesh)
       ix = quantize_lod(normalize(mesh(i)%x, bb(1), bb(2)), bits)
       iy = quantize_lod(normalize(mesh(i)%y, bb(3), bb(4)), bits)

       if (nd == 2) then
          if (hilb) then
             mesh(i)%key = hilbert2(ix, iy, bits)
          else
             mesh(i)%key = morton2(ix, iy, bits)
          end if
       else
          iz = quantize_lod(normalize(mesh(i)%z, bb(5), bb(6)), bits)
          if (hilb) then
             mesh(i)%key = hilbert3(ix, iy, iz, bits)
          else
             mesh(i)%key = morton3(ix, iy, iz, bits)
          end if
       end if
    end do
  end subroutine build_keys


  recursive subroutine mergesort(mesh,l,r)
    type(cell), intent(inout) :: mesh(:)
    integer, intent(in) :: l,r
    integer :: m
    if (l >= r) return
    m = (l+r)/2
    call mergesort(mesh,l,m)
    call mergesort(mesh,m+1,r)
    call merge_cells(mesh,l,m,r)
  end subroutine mergesort

  subroutine merge_cells(mesh,l,m,r)
    type(cell), intent(inout) :: mesh(:)
    integer, intent(in) :: l,m,r
    integer :: n1,n2,i,j,pos,t
    type(cell), allocatable :: left(:), right(:)
    n1 = m-l+1
    n2 = r-m
    allocate(left(n1), right(n2))
    left = mesh(l:m)
    right = mesh(m+1:r)
    i=1; j=1; pos=l
    do while (i<=n1 .and. j<=n2)
       if (left(i)%key <= right(j)%key) then
          mesh(pos)=left(i); i=i+1
       else
          mesh(pos)=right(j); j=j+1
       end if
       pos=pos+1
    end do
    do t=i,n1
       mesh(pos)=left(t); pos=pos+1
    end do
    do t=j,n2
       mesh(pos)=right(t); pos=pos+1
    end do
    deallocate(left,right)
  end subroutine merge_cells

end module sfc_reorder

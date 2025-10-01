module sfc_reorder
   use types
   implicit none
   intrinsic :: merge

contains

  integer(i4) function quantize(v,vmin,vmax,bits) result(q)
    real(rk), intent(in) :: v, vmin, vmax
    integer,  intent(in) :: bits
    real(rk) :: t
    if (vmax <= vmin) then
       q = 0
       return
    end if
    t = (v - vmin) / (vmax - vmin)
    if (t < 0.0_rk) t = 0.0_rk
    if (t > 1.0_rk) t = 1.0_rk
    q = int(t * real(2**bits - 1, rk), i4)
  end function quantize

  integer(i8) function morton3(ix,iy,iz,bits) result(key)
    integer(i4), intent(in) :: ix,iy,iz
    integer,    intent(in)  :: bits
    integer :: i
    integer(i8) :: k
    key = 0_i8
    do i = 0, bits-1
       k = merge(1_i8, 0_i8, btest(ix,i)); key = ior(key, ishft(k,3*i+2))
       k = merge(1_i8, 0_i8, btest(iy,i)); key = ior(key, ishft(k,3*i+1))
       k = merge(1_i8, 0_i8, btest(iz,i)); key = ior(key, ishft(k,3*i  ))
    end do
  end function morton3

  subroutine build_keys(mesh, max_bits, bb)
    type(cell), intent(inout) :: mesh(:)
    integer, intent(in)       :: max_bits
    real(rk), intent(in)      :: bb(6)
    integer(i4) :: ix,iy,iz
    integer :: i
    do i=1,size(mesh)
       ix = quantize(mesh(i)%x, bb(1), bb(2), max_bits)
       iy = quantize(mesh(i)%y, bb(3), bb(4), max_bits)
       iz = quantize(mesh(i)%z, bb(5), bb(6), max_bits)
       mesh(i)%key = morton3(ix,iy,iz,max_bits)
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

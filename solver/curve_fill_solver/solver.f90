
subroutine curve_fill_solver() bind(C, name="curve_fill_solver")
  
  use iso_c_binding, only: c_loc
  use types
  use sfc_quadtree_airfoil
  use mesh_gen
  use waffle
  use io_module
  implicit none

  ! example: simple symmetric airfoil polyline (replace with your data)
  integer, parameter :: n = 101
  real(rk) :: xa(n), ya(n)

  ! domain and meshing params
  real(rk) :: bb(4)
  integer  :: max_level, bits
  real(rk) :: hmin, beta

  type(lod_mesh) :: m
  type(cell2d), allocatable :: leaves(:)
  integer :: i

  type(cell2d), allocatable :: cells(:)

  print *, 'entered curve_fill_solver'

  call generate_cmesh(300, 300, "2412", real(2, 8), real(6.0, 8), cells)

  call build_example_airfoil(xa, ya)

  ! choose a farfield box (e.g., chord in [0,1], farfield [-5,10] x [-5,5])
  bb = (/ -5.0_rk, 10.0_rk, -5.0_rk, 5.0_rk /)

  ! meshing parameters
  max_level = 10        ! maximum quadtree depth
  hmin      = 0.01_rk   ! minimum target cell size near the airfoil
  beta      = 0.5_rk    ! growth factor in h_target = max(hmin, beta * dist)
  bits      = 20        ! quantization bits for morton keys (>= max_level is fine)

  call build_quadtree(xa,ya,n, bb, max_level, hmin, beta, m%cells, m%length)
  print *, 'built mesh starting sort'
  call build_keys_and_sort(m%cells, m%length, bits, bb)

  print *, 'nleaf = ', m%length
  print *, 'first 10 morton keys (key, level, centroid):'
  do i=1, min(m%length,10)
     write(*,'(i12,1x,i3,1x,3f12.6)') m%cells(i)%key, m%cells(i)%level, &
          0.5_rk*(m%cells(i)%xmin+m%cells(i)%xmax), 0.5_rk*(m%cells(i)%ymin+m%cells(i)%ymax), 0.0_rk
  end do

  ! emit mesh to C++ side (for visualization, etc)
  call lod_mesh_to_qt(m)

end subroutine curve_fill_solver

module waffle
  use sfc_quadtree_airfoil, only: rk
  implicit none
  contains
  subroutine build_example_airfoil(xa, ya)
    
    implicit none
    real(rk), intent(out) :: xa(:), ya(:)
    integer :: i, n
    real(rk) :: x, t, yt, c
    ! quick naca 0012-ish polyline (upper then lower) over x in [0,1]
    n = size(xa)
    c = 1.0_rk
    do i=1, n/2
      x = real(i-1, rk)/real(n/2-1, rk) * c
      t = 0.12_rk
      yt = 5.0_rk*t*(0.2969_rk*sqrt(x) - 0.1260_rk*x - 0.3516_rk*x*x + 0.2843_rk*x*x*x - 0.1015_rk*x*x*x*x)
      xa(i)   = x
      ya(i)   = yt
    end do
    do i=1, n/2
      x = real(n/2 - i, rk)/real(n/2-1, rk) * c
      t = 0.12_rk
      yt = -5.0_rk*t*(0.2969_rk*sqrt(x) - 0.1260_rk*x - 0.3516_rk*x*x + 0.2843_rk*x*x*x - 0.1015_rk*x*x*x*x)
      xa(n/2 + i) = x
      ya(n/2 + i) = yt
    end do
  end subroutine build_example_airfoil
end module waffle


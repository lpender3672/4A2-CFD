
subroutine curve_fill_solver() bind(C, name="curve_fill_solver")
  
  use iso_c_binding, only: c_loc
  use types
  use mesh_gen
  use io_module
  implicit none

  ! example: simple symmetric airfoil polyline (replace with your data)
  integer, parameter :: n = 101
  real(rk) :: xa(n), ya(n)

  ! domain and meshing params
  real(rk) :: bb(4)
  integer  :: max_level, bits
  real(rk) :: hmin, beta

  type(lod_mesh) :: mesh
  type(cell2d), allocatable :: leaves(:)
  integer :: i

  type(cell2d), allocatable :: cells(:)

  print *, 'entered curve_fill_solver'

  call generate_cmesh(512, 512, "2412", real(2, 8), real(6.0, 8), mesh)

  ! emit mesh to C++ side (for visualization, etc)
  call lod_mesh_to_qt(mesh)


end subroutine curve_fill_solver

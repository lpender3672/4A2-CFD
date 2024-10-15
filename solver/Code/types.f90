
      module types
      
      use iso_c_binding, only: C_INT, C_FLOAT, C_PTR


!     Define the derived data types used in the main program and subroutines,
!     you many need to create further variables for specific extensions

!     Appvars type contains all application variables and gas constants
      type t_appvars

!         Case name
          character(len=:), allocatable :: casename

!         Gas Properties
          real :: rgas, gam, cp, cv, fgam

!         Timestepping, smoothing and other run options
          real ::  cfl, sfac, dt, d_max, d_avg
          integer :: nsteps, nstep

!         Reference values of the primary flow variables
          real :: ro_ref, roe_ref, rov_ref

!         Size of the mesh for single block test cases
          integer :: ni, nj

!         Number of blocks and matching patches for multi-block extension
          integer :: nn, nm

      end type t_appvars

!     Boundary condition type contains inlet and outlet data
      type t_bconds

!         Single value floats at the inlet
          real :: pstag, tstag, alpha, rfin, rostag

!         Vectors of floats at the inlet
          real, dimension(:), allocatable :: p, ro

!         Single value float at the outlet
          real :: p_out

!         Block numbers of the inlet and outlet in multi-block extension
          integer :: n_in, n_out

      end type t_bconds

!     Matching patch type contains block to block interface data
      type t_match

!         Length of the matching patch and adjoining block numbers
          integer :: nk, n_1, n_2

!         Lists of indices that are coincident on both sides
          integer, dimension(:), allocatable :: i_1, j_1, i_2, j_2

      end type t_match

!     Geometry type contains the boundaries of the case geometry
      type t_geometry

!         Curve length
          integer :: ni_a, ni_b

!         Coordinate data for upper and lower domain boundaries
          real, dimension(:), allocatable :: x_a, y_a, x_b, y_b

      end type t_geometry

!     Grid type contains coordinate data and all flow variables
      type t_grid

!         Mesh size for all cases
          integer(C_INT) :: ni, nj

!         Mesh coordinate data in 2D matrices
          real(C_FLOAT), dimension(:,:), pointer :: x, y, area, lx_i, ly_i, &
              lx_j, ly_j
          real(C_FLOAT)  ::  l_min

!         Primary variables at nodes
          real(C_FLOAT), dimension(:,:), pointer :: ro, roe, rovx, rovy
!         Variables to hold cell increments
          real(C_FLOAT), dimension(:,:), pointer :: dro, droe, drovx, drovy
!         Secondary variables at nodes
          real(C_FLOAT), dimension(:,:), pointer :: p, hstag, vx, vy 
!         Logical array to store wall locations for the nodes
          logical(C_INT), dimension(:,:), pointer :: wall

      end type t_grid

      type, bind(C) :: t_grid_c
            integer(C_INT) :: ni, nj
            type(C_PTR) :: x, y, area, lx_i, ly_i, lx_j, ly_j
            real(C_FLOAT) :: l_min
            type(C_PTR) :: ro, roe, rovx, rovy
            type(C_PTR) :: dro, droe, drovx, drovy
            type(C_PTR) :: p, hstag, vx, vy
            type(C_PTR) :: wall
      end type t_grid_c

      end module types

      module conversion
      
      use iso_c_binding, only: C_INT, C_FLOAT, c_loc, c_f_pointer
      use types
      implicit none

      contains

      subroutine grid_from_c(grid_c, g)
        
            type(t_grid_c), intent(in) :: grid_c  ! Incoming C-compatible struct
            type(t_grid), intent(out) :: g     ! Fortran-specific struct
            integer :: ni, nj
        
            ! Assign scalar components directly
            ni = grid_c%ni
            nj = grid_c%nj
            g%ni = ni
            g%nj = nj

            g%l_min = grid_c%l_min
        
            !     Node coordinates in the mesh
            allocate(g%x(ni,nj),g%y(ni,nj))

            !     Cell areas and projected side lengths at the centre of each respectively
            allocate(g%area(ni-1,nj-1),g%lx_i(ni,nj-1),g%ly_i(ni,nj-1), &
            g%lx_j(ni-1,nj),g%ly_j(ni-1,nj))
      
            !     Primary flow variables in the mesh
            allocate(g%ro(ni,nj),g%rovx(ni,nj),g%rovy(ni,nj),g%roe(ni,nj))
      
            !     Cell centred primary increments
            allocate(g%dro(ni-1,nj-1),g%drovx(ni-1,nj-1), &
            g%drovy(ni-1,nj-1),g%droe(ni-1,nj-1))
      
            !     Secondary variables stored at the nodes for convenience
            allocate(g%p(ni,nj),g%hstag(ni,nj),g%vx(ni,nj),g%vy(ni,nj))

            ! allocate wall now
            allocate(g%wall(ni,nj))
            
      end subroutine grid_from_c

      subroutine grid_to_c(grid, grid_c)
        
            type(t_grid), intent(in) :: grid     ! Fortran-specific struct
            type(t_grid_c), intent(out) :: grid_c  ! Outgoing C-compatible struct
        
            ! Assign scalar components directly
            grid_c%ni = grid%ni
            grid_c%nj = grid%nj

            grid_c%l_min = grid%l_min
            
            !     Conversion to C-compatible pointers
            grid_c%x = c_loc(grid%x)
            grid_c%y = c_loc(grid%y)
            grid_c%area = c_loc(grid%area)
            grid_c%lx_i = c_loc(grid%lx_i)
            grid_c%ly_i = c_loc(grid%ly_i)
            grid_c%lx_j = c_loc(grid%lx_j)
            grid_c%ly_j = c_loc(grid%ly_j)
            grid_c%ro = c_loc(grid%ro)
            grid_c%rovx = c_loc(grid%rovx)
            grid_c%rovy = c_loc(grid%rovy)
            grid_c%roe = c_loc(grid%roe)
            grid_c%dro = c_loc(grid%dro)
            grid_c%drovx = c_loc(grid%drovx)
            grid_c%drovy = c_loc(grid%drovy)
            grid_c%droe = c_loc(grid%droe)
            grid_c%p = c_loc(grid%p)
            grid_c%hstag = c_loc(grid%hstag)
            grid_c%vx = c_loc(grid%vx)
            grid_c%vy = c_loc(grid%vy)
            grid_c%wall = c_loc(grid%wall)
        
      end subroutine grid_to_c
      
      end module conversion





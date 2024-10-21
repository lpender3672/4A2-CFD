
!     advection

!     Change to the source code directory and compile this program with 
!     "gfortran -g -o advection.x advection.f90 routines.f90". It uses two
!     Fortran source files to create the single executable program "advection.x"

!     Change to the data directory you want to run this program in and execute 
!     with "path_to_solver/advection.x"

!     You can program Fortran in free format with the .f90 extension, this means
!     you can place your statements anywhere you like. However it is useful to
!     stick to a style to improve readability. In this project I:
!         - Start the text on column 7 and don't go beyond column 80 so that you
!           can view two files side by side easily.
!         - Indent all loops, if statements and types by 4 spaces.
!         - Start a comment line with an exclamation mark in column 1.
!         - End all programs, modules, types, subroutines, do and if statements
!           explicitly.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module types
            use iso_c_binding, only: C_INT, C_FLOAT, C_PTR
            
            implicit none

      !     In Fortran you can define user derived data types. These are similar to
      !     Python dictionaries or Matlab structures and can be used to store any mix
      !     of variables you wish. They must be defined in a module that is included
      !     within every main program and subroutine that they are used.
      !     
      !     This single type stores everything useful for this program.
            type :: t_grid
                  real(C_FLOAT) :: cfl, dt_total, a, phi_inlet, phi_start
                  integer(C_INT) :: ni
                  real(C_FLOAT), dimension(:), pointer :: x, phi
            end type t_grid

            type, bind(C) :: t_grid_c
                  real(C_FLOAT) :: cfl, dt_total, a, phi_inlet, phi_start
                  integer(C_INT) :: ni
                  type(C_PTR) :: x, phi
            end type t_grid_c

      end module types

      subroutine grid_from_c(grid_c, grid)
            use iso_c_binding, only: C_INT, C_FLOAT
            use types
        
            implicit none
        
            type(t_grid_c), intent(in) :: grid_c  ! Incoming C-compatible struct
            type(t_grid), intent(out) :: grid     ! Fortran-specific struct
        
            ! Assign scalar components directly
            grid%cfl = grid_c%cfl
            grid%dt_total = grid_c%dt_total
            grid%a = grid_c%a
            grid%phi_inlet = grid_c%phi_inlet
            grid%phi_start = grid_c%phi_start
            grid%ni = grid_c%ni
        
            allocate(grid%x(grid%ni), grid%phi(grid%ni))
                
      end subroutine grid_from_c

      subroutine grid_to_c(grid, grid_c)
            use iso_c_binding, only: C_INT, C_FLOAT, c_loc, c_f_pointer
            use types
        
            implicit none
        
            type(t_grid), intent(in) :: grid     ! Fortran-specific struct
            type(t_grid_c), intent(out) :: grid_c  ! Outgoing C-compatible struct
        
            ! Assign scalar components directly
            grid_c%cfl = grid%cfl
            grid_c%dt_total = grid%dt_total
            grid_c%a = grid%a
            grid_c%phi_inlet = grid%phi_inlet
            grid_c%phi_start = grid%phi_start
            grid_c%ni = grid%ni

            grid_c%x = c_loc(grid%x)
            grid_c%phi = c_loc(grid%phi)
            
        
      end subroutine grid_to_c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     The main program begins

      subroutine advection(grid_in, grid_out) bind(C, name="advection")
            use iso_c_binding, only: C_INT, C_FLOAT

      !     Use the derived data type by using the module that contains the definition
            use types

      !     Use a library of predefined useful subroutines for this course
            use routines

      !     Never use historical implicit variable naming
            implicit none

      !     Explicitly declare the required variables, including the grid type named
      !     "g". You can also assign values to a variable at the same time as
      !     declaring them if you wanted to with "=".
            type(t_grid_c), intent(in) :: grid_in
            type(t_grid_c), intent(out) :: grid_out
            type(t_grid) :: g


            integer :: ni, n, nsteps
            real :: x_min, x_max, dt_total, dt, dx

      !     You can print to "standard out" with this command, unless you specify a
      !     file in the command you execute your program with it will just print text 
      !     in the terminal.

            call grid_from_c(grid_in, g)

            write(6,*) 'Started the advection example program '

      !     Choose the grid resolution and coordinate limits and assign them into
      !     local variables. You can write multiple statements on one line if they 
      !     are separated by a semicolon.
            x_min = 0; x_max = 1; ni = 51;

      !     Set the total time to run the simulation for
            dt_total = 0.4

      !     Fill in the run settings in the grid data type, you can access the
      !     different fields with "%". 
            g%cfl = 0.4; g%ni = ni; g%a = 1;  ! These are now set in the C struct
            g%phi_inlet = 1; g%phi_start = 0;

      !     Allocate the size of the grid and the scalar variable phi in the memory of
      !     the program now that you know their size
            ! allocate(g%x(ni),g%phi(ni))

      !     Create the grid with linearly spaced coordinates between two limits, here
      !     a separate file with useful subroutines is used called "routines.f90". One
      !     of these is linspace, which is similar to what you are used to with Numpy.
      !     The subroutine must be called with the inputs and outputs on the same
      !     line. In this case "x_min" and "x_max" are inputs, and "g%x" is both input
      !     and output, the size of the array is taken into the subroutine, the values
      !     are calculated and then returned into the same variable "g%x".
            call linspace(x_min,x_max,g%x)

      !     Calculate the grid spacing to be used for spatial derivatives
            dx = g%x(2) - g%x(1)
            write(6,*) '  Grid spacing =', dx

      !     Calculate the timestep from the CFL number "g%cfl", the convection speed
      !     "g%a" and the grid spacing "dx"
            dt = g%cfl * dx / g%a
            write(6,*) '  Timestep =', dt

      !     Calculate the number of timesteps required to cover the time period by
      !     rounding up using the Fortran intrinsic function "ceiling". If you need to
      !     continue your code on to another line place an "&" before the end.
            nsteps = ceiling(dt_total / dt)
            write(6,*) '  Total number of timesteps required =', nsteps, &
            ' for runtime of dt_total =', dt_total

      !     Set the initial value of "phi" everywhere, although "phi" is an array you 
      !     can set every element in it to the scalar value very easily
            g%phi = g%phi_start

      !     Start the time stepping do loop for "nsteps". This is the heart of the
      !     program at marches through time modifying the value of "phi"
            do n = 0, nsteps

      !         Apply the boundary condition by setting the first value to "phi_inlet"
            g%phi(1) = g%phi_inlet

      !         Complete an iteration by marching forward in time with a first order
      !         upwind scheme
            call upwind(g%phi,dt,dx,g%a,ni)

      !         Print the progress of the solution every 10 steps by using an if 
      !         statement and the mod command to check the remainder from a division.
      !         Formatting is used in the write command to help visualise the data.
            if(mod(n,10) == 0) then
                  write(6,'(A,I4)') '  Step', n
                  write(6,'(*(F4.2,1X))') g%phi
            end if

            end do

      !     At the end of the program the result can be written to file. The first
      !     step is to open a "unit". Some have special significance, like 5 & 6 which
      !     are the standard in and out files. In this case we want a new file so
      !     select unit 1. Without any further options this will give us an easy to
      !     read ASCII text file
            open(unit=1,file='cases/advection_output.txt', status='old')

      !     Write both the coordinate and scalar data to the file, each "write" 
      !     operation prints a single line in the file where each element of the array 
      !     is separated by whitespace 
            write(1,*) g%x; write(1,*) g%phi;

            call grid_to_c(g, grid_out)

            close(1)

      end subroutine advection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Subroutines are used in place of functions as they can take multiple
!     inputs and return multiple outputs. Both inputs and outputs are defined on 
!     the first line where the subroutine begins, they can be in any order.

      subroutine upwind(phi,dt,dx,a,ni)

!     A single sided upwind difference is calculated in the variable phi, the
!     value of phi is then updated by marching forward in time by one step

!     All of the variables used in this subroutine are explicitly declared. The 
!     "in", "out" or "inout" intent statements are used to specify whether a 
!     variable can be changed by the subroutine, the compiler will check that
!     this is true and can help detect bugs. 
      implicit none
      real, dimension(ni), intent(inout) :: phi
      real, intent(in) :: dt, dx, a
      integer, intent(in) :: ni 
      real, dimension(ni-1) :: dphi

!     Calculate the change in phi by indexing the array twice to create two new 
!     arrays offset by one position and then taking the difference with a single 
!     elementwise subtraction
      dphi = phi(2:ni) - phi(1:ni-1)

!     Update the value of phi by marching forward in time by dt
      phi(2:ni) = phi(2:ni) - dt * a * dphi / dx 

!     These two lines are repeated many times and as such it is important that
!     they are computed as quickly as possible. Fortran is the fastest language 
!     for these simple stenciling operations and that is why we are using it.
!     Using the vectorised expressions also makes the code extremely concise and 
!     easily parellisable if required.

      end subroutine upwind



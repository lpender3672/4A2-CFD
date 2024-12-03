
      subroutine solver(av_c, bcs_c, g_c) bind(C, name="solver")

!     The main body of the CFD solver, calls all subroutines to march towards a
!     flow solution

!     Change to the directory you want to run your case within and execute with
!     "path_to_solver/solver.x < input_casename.txt" to run with screen output
!     or "path_to_solver/solver.x < input_casename.txt > log_casename.txt &" to
!     run in the background

!     Use derived data types defined in a separate module
      use iso_c_binding, only: c_char, c_int, c_float, c_ptr
      use types
      use routines
      use io_module
      use conversion
      use mesh_io
      use patches
      use bconds
      use guess
      use gen_mesh
      use check_conv_mod

!     Don't use historical implicit variable naming
      implicit none

!     Explicitly declare the required variables
      type(t_appvars_c), intent(in) :: av_c
      type(t_bconds_c), intent(in) :: bcs_c
      type(t_grid_c), intent(inout) :: g_c

      integer :: i, len_path
      character(len=128) :: msg_bfr
      
      type(t_appvars) :: av
      type(t_bconds) :: bcs
      type(t_match), allocatable, target :: p(:)
      type(t_geometry) :: geom
      type(t_grid), allocatable, target :: g(:)
      type(t_conv_point) :: conv_point
      real :: d_max = 1, d_avg = 1, avg_of_hist = 1
      integer :: nstep, nconv = 50, ncheck = 10
      integer :: nrkuts = 4
      integer :: nrkut, n
      integer :: ni, nj, m

      integer :: ng, np

      real, allocatable :: d_avg_hist(:)
      allocate(d_avg_hist(100))
      d_avg_hist(1:100) = 1
      d_avg_hist(100) = 0
      
      call appvars_from_c(av_c, av)

      ! call write_to_qt('Hello World!')
      write(msg_bfr,*) 'Solver preprocessing started'
      call write_to_qt(msg_bfr)

!     Read in the data on the run settings
      call read_settings(av,bcs)

!     Determine whether to generate the mesh within this Fortran program or read
!     it directly from a binary file written in Python
      if(av%ni /= -1) then

          allocate(g(1))
          av%nn = 1
          av%nm = 0
          bcs%n_in = 1
          bcs%n_out = 1

          write(msg_bfr,*) 'Generating Mesh'
          call write_to_qt(msg_bfr)
!         Now the size of the grid is known, the space in memory can be 
!         allocated within the grid type
          call allocate_arrays(av,g(1),bcs)

!         Read in the case geometry
          call read_geom(av,geom)

!         Set up the mesh coordinates, interpolated between the geometry curves
          call generate_mesh(geom,g)

      else 
          write(msg_bfr,*) 'Reading Mesh from file'
          call write_to_qt(msg_bfr)
!         Read the mesh coordinates directly from file - used for extension
          call read_mesh(g, p, av, bcs)

      end if

!     Calculate cell areas and facet lengths
      do n = 1, av%nn
            call calc_areas(g(n))
      end do

      write(msg_bfr,*) 'Minimum mesh size found to be ', g%l_min
      call write_to_qt(msg_bfr)

!     Optional output call to inspect the mesh you have generated
      call write_output(av,g(1),1)

!     Check that the areas and projected lengths are correct
      call check_mesh(g(1), av)

!     Calculate the initial guess of the flowfield in the domain. There are two
!     options that can be chosen with the input argument "guesstype":
!         1. Uniform flow properties when "guesstype = 1", this is completed
!            for you already, it will allow you to get the solver started but
!            convergence will take more iterations.
!         2. A 1D varying flowfield when "guesstype = 2", assuming isentropic
!            flow in the i-direction allows a calculation of a better
!            approximation to the converged flowfield and so the time to
!            solution will be reduced. You will need to complete this option.
      do ng = 1, av%nn
          call flow_guess(av,g(ng),bcs,2)
          call set_secondary(av,g(ng))
      end do

!     Optional output call to inspect the initial guess of the flowfield
      call write_output(av,g(1),2)

      ! print grid y values here
      !write(msg_bfr,*) 'Grid y values:'
      !call write_to_qt(msg_bfr)
      !do i = 1, g%ni
      !    write(msg_bfr,*) g%y(i,:)
      !    call write_to_qt(msg_bfr)
      !end do

!     Set the length of the timestep, initially this is a constant based on a 
!     conservative guess of the mach number
      call set_timestep(av,g(1),bcs)

!     Open file to store the convergence history. This is human readable during
!     a long run by using "tail -f conv_example.csv" in a terminal window
      open(unit=3,file= trim(av%casefolder) // '/conv_' // trim(av%casename) // '.csv')

!     Initialise the "stopit" file, during long runs you can request an output
!     is written by setting the value to 1, or to terminate the calculation if
!     set to 2
      open(unit=11,file='stopit')
      write(11,*) 0; close(11);


      write(msg_bfr,*) 'Calculation started'
      call write_to_qt(msg_bfr)
!     Start the time stepping do loop for "nsteps". This is now the heart of the
!     program, you should aim to program anything inside this loop to operate as
!     efficiently as you can.

      call grids_to_qt(g, av%nn)

      do nstep = 1, av%nsteps

!         Update record of nstep to use in subroutines
          av%nstep = nstep

          do ng = 1, av%nn
            g(ng)%ro_start = g(ng)%ro
            g(ng)%roe_start = g(ng)%roe
            g(ng)%rovx_start = g(ng)%rovx
            g(ng)%rovy_start = g(ng)%rovy
          end do          

          do nrkut = 1,nrkuts
              av%dt = av%dt_total / (1 + nrkuts - nrkut)

              do ng = 1, av%nn
                  call set_secondary(av,g(ng))
              end do

              call apply_bconds(av,g,bcs)

              do ng = 1, av%nn
                  call euler_iteration(av,g(ng))
              end do

              do np = 1, av%nm
                  call apply_patch(g,p(np))
              end do
          end do

!         Write out summary every "nconv" steps and update "davg" and "dmax" 
          if(mod(av%nstep,nconv) == 0) then
              call check_conv(av,g,d_avg,d_max)
              conv_point%iters = av%nstep
              conv_point%d_max = d_max
              conv_point%d_avg = d_avg

              d_avg_hist(1:99) = d_avg_hist(2:100)
              d_avg_hist(100) = d_avg

              call conv_point_to_qt(conv_point)
                  
          end if

!         Check the solution hasn't diverged or a stop has been requested 
!         every "ncheck" steps
          if(mod(av%nstep,ncheck) == 0) then
              call check_stop(av,g(1))
          end if

!         Stop marching if converged to the desired tolerance "conlim"
          !if(d_max < av%d_max .and. d_avg < av%d_avg) then
          avg_of_hist = sum(d_avg_hist)/100
          if (all(abs(d_avg_hist / avg_of_hist - 1) < av%d_var)) then
            ! this code is modified. The original code is commented out above
            ! the calculation stops when the variation of the average residual is less than 1%
              if(d_max < av%d_max .and. d_avg < av%d_avg) then
                  write(msg_bfr,*) 'Calculation converged within bounds in', nstep,'iterations'
                  call write_to_qt(msg_bfr)
              else
                  write(msg_bfr,*) 'Calculation converged outside bounds in', nstep,'iterations'
                  call write_to_qt(msg_bfr)
              end if
              
              exit
          end if
          

          if (av%crashed) then
              exit
          end if

      end do

!     Calculation finished, call "write_output" to write the final, not 
!     necessarily converged flowfield
      write(msg_bfr,*) 'Calculation completed after', av%nstep,'iterations'
      call write_to_qt(msg_bfr)
      call write_output(av,g(1),3)

      !call grids_to_qt(g, av%nn)

!
!     Close open convergence history file
      close(3)

      end subroutine solver



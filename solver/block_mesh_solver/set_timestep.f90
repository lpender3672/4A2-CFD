
module timestep
      implicit none
      contains

      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      use io_module
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
      real :: astag
      integer :: ni, nj

      real, allocatable, dimension(:,:) :: vx, vy

      character(len=256) :: msg_bfr

      ni = g%ni; nj = g%nj
      allocate(vx(ni-1,nj-1), vy(ni-1,nj-1))

!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
      astag = sqrt(av%gam * av%rgas * bcs%tstag)

!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number

      if (av%tstep_method == 1) then
            !     Calculate the timestep using the CFL number and store it in "av%dt"
            g%dt_total = av%cfl * av%l_min / (2 * astag)
      else
            vx = 0.25 * (g%vx(1:ni-1,1:nj) + g%vx(2:ni,1:nj) + g%vx(1:ni-1,2:nj) + g%vx(2:ni,2:nj))
            vy = 0.25 * (g%vy(1:ni-1,1:nj) + g%vy(2:ni,1:nj) + g%vy(1:ni-1,2:nj) + g%vy(2:ni,2:nj))

      !     Calculate the timestep using the CFL number and store it in "av%dt"
            g%dt_total = av%cfl * g%l_min / (astag + sqrt(vx**2 + vy**2))
      
      end if

!     Print the calculated timestep and some intermediate values
      ! Dont print anymore this is run in a loop
      !write(msg_bfr,*) 'Time step set to ', av%dt_total, ' seconds'
      !call write_to_qt(msg_bfr)

      end subroutine set_timestep

end module timestep

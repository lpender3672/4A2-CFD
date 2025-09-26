
module check_stop_mod

      use types
      use io_module
      use write_output_mod
      implicit none

      contains

      subroutine check_stop(av,g)

!     This subroutine checks for divergence in the calculation by the presence
!     of NaNs in the density field

!     Explicitly declare the required variables
      use types
      use io_module
      use write_output_mod
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), allocatable, intent(in) :: g(:)
      integer :: ifstop, ng
      character(len=64) :: msg_bfr

!     Check the stop file
      !open(unit=11,file='stopit')
      !read(11,*) ifstop; close(11);

!     Check for NaNs in the density
      do ng = 1, av%nn
            if(isnan(sum(g(ng)%ro))) then
            av%crashed = 1
            write(msg_bfr,*) 'NaN detected, stopping the solver'
            call write_to_qt(msg_bfr)
            exit
            end if
      end do
     
!     Write output file if stop file is not zero
      !if(ifstop > 1) then
!      if (av%crashed) then
!          write(msg_bfr,*) '"av%crashed" modified, writing an output'
!          call write_to_qt(msg_bfr)
!          call write_output(av,g,3)
!      end if
 
!     Finish the calculation if stop file equals 2
        if(av%crashed /= 0) then
              write(msg_bfr,*) 'Solver stopped prematurely'
              call write_to_qt(msg_bfr)
              !stop
              av%crashed = 1
        end if

!     Reset the stop file      
      !ifstop = 0
      !open(unit=11,file='stopit')      
      !write(11,*) ifstop; close(11);

      end subroutine check_stop

end module check_stop_mod

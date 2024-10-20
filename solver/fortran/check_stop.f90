
      subroutine check_stop(av,g)

!     This subroutine checks for divergence in the calculation by the presence
!     of NaNs in the density field

!     Explicitly declare the required variables
      use types
      use io_module
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      integer :: ifstop
      character(len=64) :: msg_bfr

!     Check the stop file
      open(unit=11,file='stopit')
      read(11,*) ifstop; close(11);

!     Check for NaNs in the density
      if(isnan(sum(g%ro))) then
          ifstop = 2
          write(msg_bfr,*) 'NaN detected, stopping the solver'
          call write_to_qt(msg_bfr)
      end if
     
!     Write output file if stop file is not zero
      if(ifstop > 1) then
          write(msg_bfr,*) '"ifstop" modified, writing an output'
          call write_to_qt(msg_bfr)
          call write_output(av,g,3)
      end if
 
!     Finish the calculation if stop file equals 2
      if(ifstop == 2) then
          write(msg_bfr,*) 'Solver stopped prematurely'
          call write_to_qt(msg_bfr)
          !stop
          av%crashed = .true.
      end if

!     Reset the stop file      
      ifstop = 0
      open(unit=11,file='stopit')      
      write(11,*) ifstop; close(11);

      end subroutine check_stop



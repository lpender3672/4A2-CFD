      
module write_output_mod

      use types
      implicit none

      contains

      subroutine write_output(av,g,outtype)
      
!     Explicitly declare the required variables
      
      type(t_appvars), intent(in) :: av
      type(t_grid), allocatable, intent(in) :: g(:)
      integer, intent(in) :: outtype
      character(len=5) :: outname
      integer :: ng

!     Check what data to write to file depending on the value contained within
!     "outtype", options are to output grid coordinates only, grid + initial
!     guess, or the complete works 
      if(outtype == 1) then
          outname = 'coord'
      elseif(outtype == 2) then
          outname = 'guess'
      elseif(outtype == 3) then
          outname = 'final'
      end if
      
!     Open a new file to write the data into, it is an unformatted binary file
!     that takes up minimal space but contains all of the grid, flow and
!     residual variables contained within the "g" variable. It is relatively
!     straightforward to read into other programs as long as you know the
!     structure.
      open(unit=7,file= trim(av%casefolder) // '/out_' // outname // '_' // trim(av%casename) // '.bin', &
          form='unformatted',access='stream',status='replace')

      do ng = 1, av%nn
!     Write the size of the mesh
            write(7) [g(ng)%ni, g(ng)%nj]

!     Always write mesh coordinates
            write(7) g(ng)%x; write(7) g(ng)%y; 

!     Always write cell areas and projected facet lengths 
            write(7) g(ng)%area; write(7) g(ng)%lx_i; write(7) g(ng)%ly_i; 
            write(7) g(ng)%lx_j; write(7) g(ng)%ly_j;

!     Always write the wall array
            write(7) g(ng)%wall

!     Write flow solution if it has been initialised with an initial guess or 
!     has been completely solved
            if(outtype > 1) then

      !         Write primary flow variables only
            write(7) g(ng)%ro; write(7) g(ng)%roe;           
            write(7) g(ng)%rovx; write(7) g(ng)%rovy;           

            end if

!     Write cell increment data only if the flow has been solved, if you want to 
!     include any other detailed data about the solution or time marching 
!     calculations then you should write the variables here
            if(outtype > 2) then
            
      !         Write cell increments, note these arrays are smaller than the node  
      !         data that has been written above
            write(7) g(ng)%dro; write(7) g(ng)%droe;           
            write(7) g(ng)%drovx; write(7) g(ng)%drovy;           

            end if
      end do

!     Close the unit
      close(7)

      end subroutine write_output

end module write_output_mod

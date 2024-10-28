      
      subroutine apply_bconds(av,g,bcs)

!     This subroutine applies both the inlet and outlet boundary conditions, as
!     it modifies both the primary and secondary flow variables they must be
!     calculated first

!     Explicitly declare the required variables
      use types
      use debug
      use io_module
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs
      character(len=64) :: msg_bfr

!     Declare the other variables you need here
      
      real :: Tstatic(g%nj), Vinlet(g%nj)

!     At the inlet boundary the change in density is driven towards "rostag",
!     which is then used to obtain the other flow properties to match the
!     specified stagnation pressure, temperature and flow angle. 

!     To help prevent instabilities forming at the inlet boundary condition the 
!     changes in inlet density are relaxed by a factor "rfin" normally set to 
!     0.25 but it can be reduced further.

!     It is also worth checking if "ro" is greater than "rostag" and limiting 
!     the values to be slightly less than "rostag". This can prevent the solver 
!     crashing during severe transients.

      !write(msg_bfr,*) 'before applying bcs'
      !call write_to_qt(msg_bfr)
      !call find_all_NaN(g,msg_bfr)

      if(av%nstep == 1) then
          bcs%ro = g%ro(1,:)
      else
          bcs%ro = bcs%rfin * g%ro(1,:) + (1 - bcs%rfin) * bcs%ro
      endif
      bcs%ro = min(bcs%ro,0.9999 * bcs%rostag)

!     Calculate "p(1,:)", "rovx(1,:)", "rovy(1,:)" and "roe(1,:)" from the inlet 
!     "ro(:)", "pstag", "tstag" and "alpha". Also set "vx(1,:)", "vy(1,:)" and 
!     "hstag(1,:)"

      Tstatic = bcs%tstag * (bcs%ro / bcs%rostag)**(av%gam - 1)
      Vinlet = sqrt(2 * av%cp * (bcs%tstag - Tstatic))

      g%vx(1,:) = Vinlet * cos(bcs%alpha)
      g%vy(1,:) = Vinlet * sin(bcs%alpha)
      g%rovx(1,:) = bcs%ro * g%vx(1,:)
      g%rovy(1,:) = bcs%ro * g%vy(1,:)
      g%p(1,:) = bcs%pstag * (bcs%ro / bcs%rostag ) ** av%gam
      g%roe(1,:) = bcs%ro * ( av%cv * Tstatic + 0.5 * Vinlet**2 )
      g%hstag(1,:) = (g%roe(1,:) + g%p(1,:)) / bcs%ro
      ! CHECK AGAIN

      !write(msg_bfr,*) 'after applying bcs'
      !call write_to_qt(msg_bfr)
      !call find_all_NaN(g,msg_bfr)
      
!     For the outlet boundary condition set the value of "p(ni,:)" to the
!     specified value of static pressure "p_out" in "bcs"
      g%p(g%ni,1:g%nj) = bcs%p_out

      end subroutine apply_bconds



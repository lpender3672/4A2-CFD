      
      subroutine flow_guess(av,g,bcs,guesstype)

!     This calculates an initial guess of the primary flowfield variables and
!     stores them at the nodes within the mesh dataytype

!     Explicitly declare the required variables
      use types
      use routines
      use io_module
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid
      character(len=128) :: msg_bfr
      
!     Variables required for the crude guess
      real :: t_out, v_out, ro_out, lx, ly, l

!     Variables required for the improved guess, you will need to add to these
      real :: l_temp(g%nj), l_i(g%ni), v_guess(g%ni), ro_guess(g%ni), t_guess(g%ni)
      real :: mdot_exit, mach_lim, t_lim
      
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Assuming isentropic flow to the the exit plane calculate the static
!     temperature and the exit velocity
      t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
      v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
      ro_out = bcs%p_out / (av%rgas * t_out)

!     Determine which guess calcation method to use by the value of "guesstype"
      if(guesstype == 1) then

!         Store the exit density and internal energy as if they were uniform 
          g%ro = ro_out 
          g%roe  = g%ro * (av%cv * t_out + 0.5 * v_out**2)

!         Calculate the gradient of the mesh lines in the centre of the domain
!         to determine the assumed direction of the flow
          j_mid = nj / 2
          do i = 1,ni-1
              lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
              l = hypot(lx,ly)
              g%rovx(i,:) = g%ro(i,:) * v_out * ly / l
              g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
          end do

!         Copy the values to the "i = ni" nodes as an approximation
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)

!         Print the guess that has been calculated
          write(6,*) 'Crude flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
              g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)

      else if(guesstype == 2) then 

!         Calculate the length of each "i = const" line between the "j = 1" and 
!         "j = nj" boundaries of the domain and store it in the local variable
!         "l_i". You could calculate the length along each i-facet from the x 
!         and y projected lengths with "hypot" and then sum them up in the
!         second dimension with "sum". 
          l_i = sum(hypot(g%lx_i,g%ly_i),2)

!         Use the exit temperature, density and velocity calculated for the 
!         crude guess with "l_i" to estimate the mass flow rate at the exit
          mdot_exit = ro_out * v_out * l_i(ni)

!         Set a limit to the maximum allowable mach number in the initial
!         guess, call this "mach_lim", calculate the corresponding temperature,
!         called "t_lim"
          mach_lim = 0.9
          t_lim = bcs%tstag / (1 + 0.5 * (av%gam - 1) * mach_lim ** 2)

!         Now estimate the velocity and density at every "i = const" line, call 
!         the velocity "v_guess(i)" and the density "ro_guess(i)":
!             1. Assume density is constant at the exit value
!             2. Use continuity and "l_i(i)" to estimate the velocity
!             3. Assume stagnation temperature is constant for static temp
!             4. Limit the static temperature, lookup intrinsic "max"
!             5. Calculate the density throughout "ro_guess(i)"
!             6. Update the estimate of the velocity "v_guess(i)" 
          ! av%fgam = (av%gam - 1.0) / av%gam
          ro_guess(ni) = ro_out
          v_guess(ni) = v_out
          t_guess(ni) = max(t_lim, bcs%tstag - 0.5 * v_guess(ni)**2 / av%cp)
          do i = 1, ni-1
              v_guess(i) = mdot_exit / (ro_out * l_i(i))
              t_guess(i) = max(t_lim, bcs%tstag - 0.5 * v_guess(i)**2 / av%cp)
              ! this is pstatic / (rgas * tstatic)
              ro_guess(i) = ( bcs%pstag * (t_guess(i) / bcs%tstag)**(av%fgam) ) / (av%rgas * t_guess(i))
              ! update the velocity guess
              v_guess(i) = mdot_exit / (ro_guess(i) * l_i(i))
          end do

!         Direct the calculated velocity to be parallel to the "j = const"
!         gridlines for all values of i and j. This can be achieved with a 
!         similar calculation to the "j = nj/2" one that was performed in the 
!         crude guess. Then set all of ro, roe, rovx and rovy, note that roe 
!         includes the kinetic energy component of the internal energy.

          ! direction of flow is now taken at each point j instead of just the middle
          ! This is then also updated at every point
          do j = 1, nj
            do i = 1,ni-1
                  lx = g%lx_j(i,j);
                  ly = g%ly_j(i,j); 
                  l = hypot(lx,ly)
                  g%ro(i,j) = ro_guess(i)
                  g%rovx(i, j) = g%ro(i, j) * v_guess(i) * ly / l
                  g%rovy(i, j) = -g%ro(i, j) * v_guess(i) * lx / l
                  g%roe(i, j) = g%ro(i, j) * (av%cv * t_guess(i) + 0.5 * v_guess(i)**2)
              end do
          end do
              
!         Make sure the guess has been copied for the "i = ni" values too
          g%ro(ni,:) = g%ro(ni-1,:)
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)
          g%roe(ni,:) = g%roe(ni-1,:)


!         Print the first elements of the guess like for the crude guess
          write(6,*) 'Flow guess calaculated using improved method'
          write(6,*) 'At g(1,1): ro =', g%ro(1,1), 'roe =', g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*) 'At g(ni,1): ro =', g%ro(ni,1), 'roe =', g%roe(ni,1), 'rovx =', g%rovx(ni,1), 'rovy =', g%rovy(ni,1)
          write(6,*)

      end if

!     The initial guess values derived from the boundary conditions are also
!     useful as a reference to non-dimensionalise the convergence
      av%ro_ref = sum(g%ro(1,:)) / nj
      av%roe_ref = sum(g%roe(1,:)) / nj
      av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj

      end subroutine flow_guess



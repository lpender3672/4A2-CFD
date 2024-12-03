
module check_conv_mod
      
      use types
      use routines
      use io_module

      implicit none
      
      contains

      subroutine check_conv(av,g,d_avg,d_max)

!     This subroutine checks the residuals in all primary flow variables and
!     prints their values, you should not need to change this subroutine

!     Explicitly declare the required variables

      type(t_appvars), intent(in) :: av
      type(t_grid), allocatable, intent(in) :: g(:)
      real, intent(out) :: d_avg, d_max
      real, allocatable, dimension(:,:) :: dro, droe, drovx, drovy
      integer :: ijx_max(2), ijy_max(2), ij_max(2), ncells
      integer :: ijx_max_grid(2), ijy_max_grid(2)
      integer :: ng, ncells_sum
      real :: dro_max, drovx_max, drovy_max, droe_max, dro_avg, drovx_avg, &
          drovy_avg, droe_avg
      real :: dro_sum, droe_sum, drovx_sum, drovy_sum, drovx_max_grid, drovy_max_grid
      character(len=100) :: fmt_step
      character(len=128) :: msg_bfr

!     Get the number of cells from the size of the residual arrays

      dro_max = 0
      droe_max = 0
      drovx_max = 0
      drovy_max = 0

      dro_sum = 0
      drovx_sum = 0
      drovy_sum = 0
      droe_sum = 0

      ncells_sum = 0

      do ng = 1, av%nn

            ncells_sum = ncells_sum + size(g(ng)%dro)

!     Use "abs" to make all residual values positive and store locally
            dro = abs(g(ng)%dro);
            droe = abs(g(ng)%droe);
            drovx = abs(g(ng)%drovx); 
            drovy = abs(g(ng)%drovy);

!     Calculate the mean changes for each variable
            dro_sum = dro_sum + sum(abs(dro))
            droe_sum = droe_sum + sum(abs(droe))
            drovx_sum = drovx_sum + sum(abs(drovx))
            drovy_sum = drovy_sum + sum(abs(drovy))

!     Find the maximum value of change for the momenta and the positions
            dro_max = max(dro_max, maxval(dro) / av%ro_ref)
            droe_max = max(droe_max, maxval(droe) / av%roe_ref)
            ijx_max_grid = maxloc(drovx);
            drovx_max_grid = drovx(ijx_max_grid(1),ijx_max_grid(2)) / av%rov_ref
            if (drovx_max_grid > drovx_max) then
                drovx_max = drovx_max_grid
                ijx_max = ijx_max_grid
            end if
            ijy_max_grid = maxloc(drovy);
            drovy_max_grid = drovy(ijy_max_grid(1),ijy_max_grid(2)) / av%rov_ref
            if (drovy_max_grid > drovy_max) then
                drovy_max = drovy_max_grid
                ijy_max = ijy_max_grid
            end if

      end do 

      dro_avg = dro_sum / (ncells_sum * av%ro_ref)
      droe_avg = droe_sum / (ncells_sum * av%roe_ref)
      drovx_avg = drovx_sum / (ncells_sum * av%rov_ref)
      drovy_avg = drovy_sum / (ncells_sum * av%rov_ref)

!     Store single values as the maximum of either the x or y-momentum
      if(drovx_avg > drovy_avg) then
          d_max = drovx_max; d_avg = drovx_avg; ij_max = ijx_max;
      else
          d_max = drovy_max; d_avg = drovy_avg; ij_max = ijy_max;
      end if

!     Write the average and maximum changes in the primary variables to unit 3
!     for convergenge plotting
      write(3,'(i13,8e15.6)') av%nstep, dro_avg, droe_avg, drovx_avg, &
          drovy_avg, dro_max, droe_max, drovx_max, drovy_max

!     Write a short human readable output summary to the screen.
      !write(6,*) 'Time step number ', av%nstep

      fmt_step = '(a,i6,a,e10.3,a,i4,a,i4,a,e10.3)'
      write(msg_bfr,fmt_step) 'iteration ', av%nstep, ', d_max =', d_max, ' at i =', ij_max(1), ', j =', ij_max(2), ', d_avg =', d_avg
      call write_to_qt(msg_bfr)

      end subroutine check_conv


end module check_conv_mod
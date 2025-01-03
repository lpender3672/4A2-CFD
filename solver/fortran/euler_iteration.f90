
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      use io_module
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i, mom_xi, mom_yi, av_p_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j, mom_xj, mom_yj, av_p_j
      real, dimension(g%ni-1,g%nj-1) :: dt_over_area
      integer :: i, j, ni, nj
      character(len=64) :: msg_bfr

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

      !write(msg_bfr,*) 'before fluxing'
      !call write_to_qt(msg_bfr)
      !call find_all_NaN(g,msg_bfr)

      av_p_i = (g%p(1:ni,1:nj-1) + g%p(1:ni,2:nj)) / 2
      av_p_j = (g%p(1:ni-1,1:nj) + g%p(2:ni,1:nj)) / 2

      dt_over_area = g%dt / g%area
      
      ! g%lx_i(ni,nj-1),g%ly_i(ni,nj-1), g%lx_j(ni-1,nj),g%ly_j(ni-1,nj)

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
      mass_i = 0.5 * (g%rovx(1:ni,1:nj-1) + g%rovx(1:ni,2:nj)) * g%lx_i + 0.5 * (g%rovy(1:ni,1:nj-1) + g%rovy(1:ni,2:nj)) * g%ly_i
      mass_j = 0.5 * (g%rovx(1:ni-1,1:nj) + g%rovx(2:ni,1:nj)) * g%lx_j + 0.5 * (g%rovy(1:ni-1,1:nj) + g%rovy(2:ni,1:nj)) * g%ly_j
     
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
      call sum_fluxes(av,mass_i,mass_j, dt_over_area, g%ro, g%ro_start, g%dro_1, g%dro)

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
      flux_i = mass_i * 0.5 * (g%hstag(1:ni,1:nj-1) + g%hstag(1:ni,2:nj))
      flux_j = mass_j * 0.5 * (g%hstag(1:ni-1,1:nj) + g%hstag(2:ni,1:nj))

!     Update the internal energy with enthalpy fluxes
      call sum_fluxes(av,flux_i,flux_j, dt_over_area, g%roe, g%roe_start, g%droe_1, g%droe)

!     Setup the x-momentum equation including momentum flux and pressure forces
      mom_xi = mass_i * 0.5 * (g%vx(1:ni,1:nj-1) + g%vx(1:ni,2:nj)) + av_p_i * g%lx_i
      mom_xj = mass_j * 0.5 * (g%vx(1:ni-1,1:nj) + g%vx(2:ni,1:nj)) + av_p_j * g%lx_j

!     Update the x-momentum with momentum flux
      call sum_fluxes(av,mom_xi,mom_xj, dt_over_area, g%rovx, g%rovx_start, g%drovx_1, g%drovx)

!     Setup the y-momentum equation including momentum flux and pressure forces
      mom_yi = mass_i * 0.5 * (g%vy(1:ni,1:nj-1) + g%vy(1:ni,2:nj)) + av_p_i * g%ly_i
      mom_yj = mass_j * 0.5 * (g%vy(1:ni-1,1:nj) + g%vy(2:ni,1:nj)) + av_p_j * g%ly_j

!     Update the y-momentum with momentum flux
      call sum_fluxes(av,mom_yi,mom_yj, dt_over_area, g%rovy, g%rovy_start, g%drovy_1, g%drovy)

      !write(msg_bfr,*) 'after fluxing'
      !call write_to_qt(msg_bfr)
      !call find_all_NaN(g,msg_bfr)

!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro, g%corr_ro)
      call smooth_array(av,g%roe, g%corr_roe)
      call smooth_array(av,g%rovx, g%corr_rovx)
      call smooth_array(av,g%rovy, g%corr_rovy)

      end subroutine euler_iteration



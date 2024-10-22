
      subroutine euler_iteration(av,g)

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use stencils
      use io_module
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i, mom_xi, mom_yi, av_vx_i, av_vy_i, av_hstag_i, av_p_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j, mom_xj, mom_yj, av_vx_j, av_vy_j, av_hstag_j, av_p_j
      integer :: i, j, ni, nj
      character(len=64) :: msg_bfr

!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

      av_vx_i = (g%vx(1:ni,1:nj-1) + g%vx(1:ni,2:nj)) / 2
      av_vx_j = (g%vx(1:ni-1,1:nj) + g%vx(2:ni,1:nj)) / 2
      av_vy_i = (g%vy(1:ni,1:nj-1) + g%vy(1:ni,2:nj)) / 2
      av_vy_j = (g%vy(1:ni-1,1:nj) + g%vy(2:ni,1:nj)) / 2
      av_p_i = (g%p(1:ni,1:nj-1) + g%p(1:ni,2:nj)) / 2
      av_p_j = (g%p(1:ni-1,1:nj) + g%p(2:ni,1:nj)) / 2
      av_hstag_i = (g%hstag(1:ni,1:nj-1) + g%hstag(1:ni,2:nj)) / 2
      av_hstag_j = (g%hstag(1:ni-1,1:nj) + g%hstag(2:ni,1:nj)) / 2
      
      ! g%lx_i(ni,nj-1),g%ly_i(ni,nj-1), g%lx_j(ni-1,nj),g%ly_j(ni-1,nj)

!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
      mass_i = g%ro(1:ni,1:nj-1) * (av_vx_i * g%lx_i + av_vy_i * g%ly_i)
      mass_j = g%ro(1:ni-1,1:nj) * (av_vx_j * g%lx_j + av_vy_j * g%ly_j)
     
!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 

!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT
      call sum_fluxes(av,mass_i,mass_j, g%area, g%ro, g%dro)

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT
      flux_i = mass_i * av_hstag_i
      flux_j = mass_j * av_hstag_j

!     Update the internal energy with enthalpy fluxes
!     INSERT
      call sum_fluxes(av,flux_i,flux_j, g%area, g%roe, g%droe)

!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT
      mom_xi = mass_i * av_vx_i + av_p_i * g%lx_i
      mom_xj = mass_j * av_vx_j + av_p_j * g%lx_j

!     Update the x-momentum with momentum flux
!     INSERT
      call sum_fluxes(av,mom_xi,mom_xj, g%area, g%rovx, g%drovx)

!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT
      mom_yi = mass_i * av_vy_i + av_p_i * g%ly_i
      mom_yj = mass_j * av_vy_j + av_p_j * g%ly_j

!     Update the y-momentum with momentum flux
!     INSERT
      call sum_fluxes(av,mom_yi,mom_yj, g%area, g%rovy, g%drovy)

!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro)
      call smooth_array(av,g%roe)
      call smooth_array(av,g%rovx)
      call smooth_array(av,g%rovy)
      

      end subroutine euler_iteration



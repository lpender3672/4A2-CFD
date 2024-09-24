      
      subroutine allocate_arrays(av,g,bcs)

!     Allocate memory for all arrays in the grid and bcs datatype, this has been
!     completed for the basic solver. If you add further arrays to the code in
!     the extensions you will need to allocate them here.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs
      integer :: ni, nj

!     Get the size of the mesh and store locally for convenience
      ni = av%ni; nj = av%nj;

!     Copy the mesh size to the grid datatype
      g%ni = ni; g%nj = nj;

!     Wall flag array
      allocate(g%wall(ni,nj))

!     Arrays to store static conditions at the inlet plane
      allocate(bcs%ro(nj),bcs%p(nj))

!     Node coordinates in the mesh
      allocate(g%x(ni,nj),g%y(ni,nj))

!     Cell areas and projected side lengths at the centre of each respectively
      allocate(g%area(ni-1,nj-1),g%lx_i(ni,nj-1),g%ly_i(ni,nj-1), &
          g%lx_j(ni-1,nj),g%ly_j(ni-1,nj))

!     Primary flow variables in the mesh
      allocate(g%ro(ni,nj),g%rovx(ni,nj),g%rovy(ni,nj),g%roe(ni,nj))

!     Cell centred primary increments
      allocate(g%dro(ni-1,nj-1),g%drovx(ni-1,nj-1), &
          g%drovy(ni-1,nj-1),g%droe(ni-1,nj-1))

!     Secondary variables stored at the nodes for convenience
      allocate(g%p(ni,nj),g%hstag(ni,nj),g%vx(ni,nj),g%vy(ni,nj))


      end subroutine allocate_arrays



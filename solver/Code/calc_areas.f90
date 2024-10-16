      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      real, dimension(:, :), allocatable :: ax, ay, bx, by
      integer :: ni, nj

!     Declare integers or any extra variables you need here
      allocate(ax(g%ni-1,g%nj-1), ay(g%ni-1,g%nj-1))
      allocate(bx(g%ni-1,g%nj-1), by(g%ni-1,g%nj-1))

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
!     Cross product of diagonal vectors 0.5 * abs( ax*by - ay*bx )
!     with a = p_i+1,j+1 - p_i,j and b = p_i,j+1 - p_i,j
      ax = g%x(2:,2:) - g%x(1:-1,1:-1)
      ay = g%y(2:,2:) - g%y(1:-1,1:-1)
      bx = g%x(1:-1,2:) - g%x(1:-1,1:-1)
      by = g%y(1:-1,2:) - g%y(1:-1,1:-1)

      g%area = 0.5 * abs( ax*by - ay*bx )


!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
      g%lx_i = g%y(:,2:) - g%y(:,1:-1)
      g%ly_i = g%x(:,1:-1) - g%x(:,2:)

!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
      g%lx_j = g%y(1:-1,:) - g%y(2:,:)
      g%ly_j = g%x(2:,:) - g%x(1:-1,:)

!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.

      g%l_min = min(minval( hypot(g%lx_i, g%ly_i) ), minval( hypot(g%lx_j, g%ly_j)))
!
!     Print the overall minimum length size that has been calculated
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', g%l_min
      write(6,*)

      end subroutine calc_areas

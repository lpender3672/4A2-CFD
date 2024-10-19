      
      subroutine check_mesh(g, av)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g
      type(t_appvars), intent(inout) :: av
      real :: area_min, dx_error, dy_error, tol
      integer :: ni, nj

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Exact checking of floating point numbers never goes well, define a
!     small tolerance to use for a comparative operation instead
      tol = 1e-4 * g%l_min

!     Check that all of the cell areas are positive, either with the intrinsic
!     "minval" function or with nested do loops. Print the output to the screen
!     and flag negative numbers as an error with an if statement to "stop" the
!     program
      area_min = minval(g%area)
      if (area_min < tol) then
          write(6,*) 'Negative cell area found: ', area_min, ' at i,j = ', &
              merge(ni, ni-1, area_min == g%area(ni,nj)), &
              merge(nj, nj-1, area_min == g%area(ni,nj))
          !stop
          av%crashed = .true.
      end if

!     Next check that the sum of the edge vectors around every quadrilateral is 
!     very nearly zero in both the x and y-coordinate directions. You can
!     complete this with some elementwise addition of the arrays and use of the
!     "maxval" and "abs" intrinsic functions.

!      g%area(ni-1,nj-1) g%lx_i(ni,nj-1) g%ly_i(ni,nj-1) g%lx_j(ni-1,nj) g%ly_j(ni-1,nj))

      dx_error = maxval(g%lx_i(:-1,:) - g%lx_i(2:, :) + g%lx_j(:,:-1) - g%lx_j(:, 2:))
      dy_error = maxval(g%ly_i(:-1,:) - g%ly_i(2:, :) + g%ly_j(:,:-1) - g%ly_j(:, 2:))
      if (dx_error > tol .or. dy_error > tol) then
          write(6,*) 'Edge vector sum error: ', dx_error, dy_error
          !stop
          av%crashed = .true.
      end if

!     It may be worthwhile to complete some other checks, the prevous call to
!     the "write_output" subroutine has written a file that you can read and
!     postprocess using the Python script plot_mesh.py. This program also has
!     access to all of the mesh parameters used within the solver that you could
!     inspect graphically.

!     Print a blank 
      write(6,*) 'Mesh check passed'
      write(6,*)

      end subroutine check_mesh

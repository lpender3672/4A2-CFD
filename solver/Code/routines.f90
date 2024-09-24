
      module routines

!     Some basic routines to manipulate the data, including them in a module
!     allows them to use assumed shape arrays to recieve the data from other
!     subroutines 

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dist(x,y,scale_flag,s)

!     Calculate relative distance along a curve in 2D
     
!     Explicitly declare the required variables, the size is assumed
      implicit none
      real, intent(in) :: x(:), y(:)
      integer, intent(in) :: scale_flag
      real, intent(out) :: s(:)
      real, allocatable :: d(:)
      integer :: n, nn

!     Length of the vectors
      nn = size(x)
      allocate(d(nn-1))

!     Length of each segment
      d = hypot(x(2:nn) - x(1:nn-1),y(2:nn) - y(1:nn-1))

!     Calculate cumilative distance by adding each length in turn
      s(1) = 0
      do n = 2,nn
          s(n) = s(n-1) + d(n-1) 
      end do

!     Scale the cumilative distance to end at 1
      if(scale_flag == 1) s = s / s(nn)

      end subroutine dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine linspace(x1,x2,x)

!     Create a vector of linearly space values between two limits
     
!     Explicitly declare the required variables, the size is assumed
      implicit none
      real, intent(in) :: x1, x2
      real, intent(out) :: x(:)
      integer :: n, nn
      real :: dx

!     Length of the output vector
      nn = size(x)

!     Step size between points     
      dx = (x2 - x1) / float(nn - 1)

!     Loop over all indices of x and calculate its value
      do n = 1,nn
          x(n) = x1 + (n - 1) * dx 
      end do

      end subroutine linspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      subroutine interp(x,y,xp,yp)

!     Piecewise linear 1D interpolation of a curve to new values, input data
!     must be monotonically increasing, can be non-uniformly spaced
     
!     Explicitly declare the required variables, the sizes are assumed 
      implicit none
      real, intent(in) :: x(:), y(:), xp(:)
      real, intent(out) :: yp(:)
      integer :: n, p, nn, np
      real :: s

!     Length of the input vectors
      nn = size(x)
      np = size(xp)
     
!     Loop over all new values of x 
      do p = 1,np
      
!         Find the index of the interval that the new x value lies within
          n = 1 
          do while(x(n+1) < xp(p))
              n = n + 1
          end do
     
!         Calculate weight between values at the ends of the interval 
          s = (xp(p) - x(n)) / (x(n+1) - x(n))

!         Interpolated value using weight value
          yp(p) = (1.0 - s) * y(n) + s * y(n+1)
      
      end do

      end subroutine interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module routines



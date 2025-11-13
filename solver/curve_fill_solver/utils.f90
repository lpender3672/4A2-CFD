
module general_utils

  use iso_fortran_env, only: real64
  use types
  implicit none

  real(8), parameter :: PI = 3.1415926535897932384626433832795D0
  real(8), parameter :: EPS = 1.0D-16

  contains

  subroutine gaussian_kernel_2d(size_k, sigma, kernel)
    integer, intent(in) :: size_k
    real(8), intent(in) :: sigma
    real(8), allocatable, intent(out) :: kernel(:,:)  ! (size_k,size_k)

    integer :: i, j, half, idx, jdx
    real(8) :: s2, val, sumv

    half = size_k/2
    allocate(kernel(size_k, size_k))
    s2 = 2.0D0 * sigma*sigma
    sumv = 0.0D0
    do i = 1, size_k
      do j = 1, size_k
        idx = i - half - 1
        jdx = j - half - 1
        val = exp(-(real(idx,8)**2 + real(jdx,8)**2) / s2)
        kernel(i,j) = val
        sumv = sumv + val
      end do
    end do
    if (sumv <= 0.0D0) sumv = 1.0D0
    kernel = kernel / sumv
  end subroutine gaussian_kernel_2d

  subroutine smooth2d(data, size_k, sigma, out)
    real(8), intent(in)  :: data(:, :)
    integer, intent(in)  :: size_k
    real(8), intent(in)  :: sigma
    real(8), intent(out) :: out(size(data,1), size(data,2))

    integer :: n, m, i, j, ii, jj, di, dj, half
    real(8), allocatable :: ker(:,:)
    real(8) :: acc

    n = size(data,1); m = size(data,2)
    call gaussian_kernel_2d(size_k, sigma, ker)
    half = size_k/2

    do i = 1, n
      do j = 1, m
        acc = 0.0D0
        do di = -half, +half
          do dj = -half, +half
            ii = reflect_index(i + di, n)
            jj = reflect_index(j + dj, m)
            acc = acc + data(ii, jj) * ker(di+half+1, dj+half+1)
          end do
        end do
        out(i,j) = acc
      end do
    end do
  end subroutine smooth2d


  subroutine rotate_points(points_in, angle_rad, points_out)
    ! points_in/out: (N,2)
    real(8), intent(in)  :: points_in(:, :)
    real(8), intent(in)  :: angle_rad
    real(8), intent(out) :: points_out(size(points_in,1), 2)
    integer :: n
    real(8) :: c, s
    n = size(points_in,1)
    c = cos(angle_rad)
    s = sin(angle_rad)
    ! R^T = [[c, s],[-s, c]]
    points_out(:,1) =  points_in(:,1)*c + points_in(:,2)*s
    points_out(:,2) = -points_in(:,1)*s + points_in(:,2)*c
  end subroutine rotate_points

    pure subroutine linspace(a, b, n, x)
    real(8), intent(in) :: a, b
    integer, intent(in) :: n
    real(8), intent(out) :: x(n)
    integer :: i
    if (n == 1) then
      x(1) = a
    else
      do i = 1, n
        x(i) = a + (b - a) * real(i - 1,8) / real(n - 1,8)
      end do
    end if
  end subroutine linspace

  pure subroutine meshgrid_xy(xv, yv, X, Y)

    real(8), intent(in)  :: xv(:), yv(:)
    real(8), intent(out) :: X(size(yv), size(xv)), Y(size(yv), size(xv))
    integer :: i, j, n, m

    n = size(yv)
    m = size(xv)

    do i = 1, n
        do j = 1, m
        X(i,j) = xv(j)
        Y(i,j) = yv(i)
        end do
    end do
    end subroutine meshgrid_xy

      pure function clamp(x, a, b) result(y)
    real(8), intent(in) :: x, a, b
    real(8) :: y
    y = max(a, min(b, x))
  end function clamp

  pure function reflect_index(i, n) result(r)

    integer, intent(in) :: i, n
    integer :: r, t
    if (n <= 1) then
      r = 1
      return
    end if
    t = i
    do
      if (t < 1) then
        t = 2 - t
      else if (t > n) then
        t = 2*n - t
      else
        exit
      end if
    end do
    r = t
  end function reflect_index

end module
module mesh_utils
  implicit none
  private
  public :: rotate_points, naca4_airfoil, transform_airfoil
  public :: curvature_at_nearest_point
  public :: gaussian_kernel_2d, smooth2d
  public :: calc_lod

  real(8), parameter :: PI = 3.1415926535897932384626433832795D0
  real(8), parameter :: EPS = 1.0D-16

contains

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

  subroutine naca4_airfoil(code, n, closed_te, poly)
    character(len=*), intent(in) :: code
    integer, intent(in) :: n
    logical, intent(in) :: closed_te
    real(8), allocatable, intent(out) :: poly(:,:)  ! (2*n, 2)

    integer :: i
    real(8) :: m, p, t, k5
    real(8), allocatable :: beta(:), x(:), yc(:), dyc_dx(:), yt(:)
    real(8), allocatable :: theta(:), xu(:), yu(:), xl(:), yl(:)

    if (len_trim(code) /= 4) stop "naca4_airfoil: code must have length 4"
    if (.not. all([(lge(code(i:i),'0') .and. lle(code(i:i),'9'), i=1,4)])) &
         stop "naca4_airfoil: code must be numeric"

    m = real(iachar(code(1:1)) - iachar('0'),8) / 100.0D0
    p = real(iachar(code(2:2)) - iachar('0'),8) / 10.0D0
    t = real( (iachar(code(3:3)) - iachar('0'))*10 + (iachar(code(4:4)) - iachar('0')), 8 ) / 100.0D0

    allocate(beta(n), x(n), yc(n), dyc_dx(n), yt(n), theta(n), xu(n), yu(n), xl(n), yl(n))
    call linspace(0.0D0, PI, n, beta)
    do i = 1, n
      x(i) = 0.5D0 * (1.0D0 - cos(beta(i)))
    end do

    if (closed_te) then
      k5 = -0.1015D0
    else
      k5 = -0.1036D0
    end if

    do i = 1, n
      ! clip x to [0,1] implicitly since construction yields that range
      yt(i) = 5.0D0 * t * ( 0.2969D0*sqrt(max(0.0D0, x(i))) - 0.1260D0*x(i) &
                 - 0.3516D0*x(i)**2 + 0.2843D0*x(i)**3 + k5*x(i)**4 )
    end do

    yc     = 0.0D0
    dyc_dx = 0.0D0
    if (m > 0.0D0 .and. p > 0.0D0) then
      do i = 1, n
        if (x(i) < p) then
          yc(i)     = m/(p**2)       * (2.0D0*p*x(i) - x(i)**2)
          dyc_dx(i) = 2.0D0*m/(p**2) * (p - x(i))
        else
          yc(i)     = m/((1.0D0-p)**2) * ((1.0D0 - 2.0D0*p) + 2.0D0*p*x(i) - x(i)**2)
          dyc_dx(i) = 2.0D0*m/((1.0D0-p)**2) * (p - x(i))
        end if
      end do
    end if

    do i = 1, n
      theta(i) = atan(dyc_dx(i))
      xu(i) = x(i) - yt(i) * sin(theta(i))
      yu(i) = yc(i) + yt(i) * cos(theta(i))
      xl(i) = x(i) + yt(i) * sin(theta(i))
      yl(i) = yc(i) - yt(i) * cos(theta(i))
    end do

    allocate(poly(2*n, 2))
    ! upper: 1..n in order
    do i = 1, n
      poly(i,1) = xu(i)
      poly(i,2) = yu(i)
    end do
    ! lower reversed: n..1
    do i = 1, n
      poly(n+i,1) = xl(n - i + 1)
      poly(n+i,2) = yl(n - i + 1)
    end do
  end subroutine naca4_airfoil

  subroutine transform_airfoil(poly_in, chord, alpha_deg, origin, poly_out)
    real(8), intent(in)  :: poly_in(:, :)
    real(8), intent(in)  :: chord, alpha_deg
    real(8), intent(in)  :: origin(2)
    real(8), intent(out) :: poly_out(size(poly_in,1), 2)
    real(8), allocatable :: tmp(:,:)
    real(8) :: angle_rad
    integer :: n

    n = size(poly_in,1)
    allocate(tmp(n,2))
    tmp(:,1) = poly_in(:,1) * chord
    tmp(:,2) = poly_in(:,2) * chord
    angle_rad = alpha_deg * PI / 180.0D0
    call rotate_points(tmp, angle_rad, poly_out)
    poly_out(:,1) = poly_out(:,1) + origin(1)
    poly_out(:,2) = poly_out(:,2) + origin(2)
  end subroutine transform_airfoil

  subroutine curvature_at_nearest_point(px, py, poly, dist, curvature)
    real(8), intent(in) :: px, py
    real(8), intent(in) :: poly(:, :)        ! (N,2)
    real(8), intent(out):: dist, curvature

    integer :: n, i, idx, i_min
    real(8) :: vx, vy, wx, wy, vv, t, projx, projy, dx, dy, d2, d2_min
    real(8) :: x_prev, x_curr, x_next, y_prev, y_curr, y_next
    real(8) :: dx1, dy1, ddx, ddy, denom

    n = size(poly,1)
    if (n < 3) then
      dist = 0.0D0; curvature = 0.0D0
      return
    end if

    d2_min = 1.0D300
    i_min = 1
    t = 0.0D0

    do i = 1, n-1
      vx = poly(i+1,1) - poly(i,1)
      vy = poly(i+1,2) - poly(i,2)
      wx = px - poly(i,1)
      wy = py - poly(i,2)
      vv = vx*vx + vy*vy
      if (vv <= 0.0D0) vv = EPS
      ! projection parameter clamped
      t = (wx*vx + wy*vy) / vv
      if (t < 0.0D0) t = 0.0D0
      if (t > 1.0D0) t = 1.0D0

      projx = poly(i,1) + t*vx
      projy = poly(i,2) + t*vy
      dx = px - projx
      dy = py - projy
      d2 = dx*dx + dy*dy
      if (d2 < d2_min) then
        d2_min = d2
        i_min = i
      end if
    end do

    dist = sqrt(d2_min)

    ! choose closer endpoint of segment i_min
    if (t < 0.5D0) then
      idx = i_min
    else
      idx = i_min + 1
    end if
    idx = max(2, min(n-1, idx))

    x_prev = poly(idx-1,1); x_curr = poly(idx,1); x_next = poly(idx+1,1)
    y_prev = poly(idx-1,2); y_curr = poly(idx,2); y_next = poly(idx+1,2)

    dx1 = x_next - x_prev
    dy1 = y_next - y_prev
    ddx = x_next - 2.0D0*x_curr + x_prev
    ddy = y_next - 2.0D0*y_curr + y_prev
    denom = (dx1*dx1 + dy1*dy1)**1.5D0
    if (denom <= 0.0D0) denom = EPS
    curvature = abs(dx1*ddy - dy1*ddx) / denom
  end subroutine curvature_at_nearest_point

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

  subroutine calc_lod(n, m, foil, max_level, X, Y, ILOD)
    integer, intent(in) :: n, m, max_level
    real(8), intent(in) :: foil(:, :)              ! (Nf,2) polyline
    real(8), intent(out):: X(n, m), Y(n, m)
    integer, intent(out):: ILOD(n, m)

    real(8), allocatable :: xv(:), yv(:), KAPPA(:,:), DIST(:,:), LOD(:,:)
    real(8), allocatable :: NLOD(:,:), tmp(:,:)
    real(8) :: min_k_pos, max_dist, min_nlod, max_nlod
    integer :: i, j

    allocate(xv(m), yv(n))
    if (m > 1) then
      call linspace(0.0D0, real(m,8), m, xv)
    else
      xv(1) = 0.0D0
    end if
    if (n > 1) then
      call linspace(0.0D0, real(n,8), n, yv)
    else
      yv(1) = 0.0D0
    end if

    call meshgrid_xy(xv, yv, X, Y)

    allocate(KAPPA(n,m), DIST(n,m), LOD(n,m))
    KAPPA = 0.0D0; DIST = 0.0D0; LOD = 0.0D0

    do i = 1, n
      do j = 1, m
        call curvature_at_nearest_point(X(i,j), Y(i,j), foil, DIST(i,j), KAPPA(i,j))
      end do
    end do

    ! LKAPPA = log(KAPPA / min(KAPPA[KAPPA>0]))
    min_k_pos = 1.0D300
    do i = 1, n
      do j = 1, m
        if (KAPPA(i,j) > 0.0D0) min_k_pos = min(min_k_pos, KAPPA(i,j))
      end do
    end do
    if (min_k_pos == 1.0D300) min_k_pos = max(1.0D0, maxval(KAPPA))  ! fallback
    if (min_k_pos <= 0.0D0) min_k_pos = 1.0D0

    ! IDIST = 1 - DIST/max(DIST)
    max_dist = max(maxval(DIST), EPS)

    do i = 1, n
      do j = 1, m
        LOD(i,j) = ( 1.0D0 - DIST(i,j)/max_dist ) * log( max(KAPPA(i,j)/min_k_pos, EPS) )
      end do
    end do

    ! Smooth: size=101, sigma=10
    allocate(NLOD(n,m))
    call smooth2d(LOD, 101, 10.0D0, NLOD)

    min_nlod = minval(NLOD)
    max_nlod = maxval(NLOD)
    if (max_nlod - min_nlod <= EPS) then
      NLOD = 0.0D0
    else
      NLOD = (NLOD - min_nlod) / (max_nlod - min_nlod) * real(max_level,8)
    end if

    ! ILOD = max_level - round( clip(NLOD, 0, max_level) )
    do i = 1, n
      do j = 1, m
        NLOD(i,j) = clamp(NLOD(i,j), 0.0D0, real(max_level,8))
        ILOD(i,j) = max_level - nint(NLOD(i,j))
      end do
    end do
  end subroutine calc_lod

end module mesh_utils

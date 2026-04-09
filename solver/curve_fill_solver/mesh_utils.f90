module mesh_utils
  use iso_fortran_env, only: real64
  use types
  use general_utils
  implicit none
  private
  public :: naca4_airfoil, transform_airfoil
  public :: nearest_idx, curvature_at_idx, dist_xy_to_xy
  public :: poly_stats

contains

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

  subroutine nearest_idx(px, py, poly, idx_near)
    real(8), intent(in) :: px, py
    real(8), intent(in) :: poly(:, :)        ! (N,2)
    integer, intent(out):: idx_near

    integer :: n, i
    real(8) :: dx, dy, d2, d2_min

    n = size(poly,1)
    d2_min = 1.0D300
    idx_near = 1

    do i = 1, n
      dx = px - poly(i,1)
      dy = py - poly(i,2)
      d2 = dx*dx + dy*dy
      if (d2 < d2_min) then
        d2_min = d2
        idx_near = i
      end if
    end do
  end subroutine nearest_idx

  subroutine dist_xy_to_xy(ax, ay, bx, by, dist)
    real(8), intent(in) :: ax, ay, bx, by
    real(8), intent(out) :: dist

    real(8) :: dx, dy, d2

    dx = ax - bx
    dy = ay - by
    d2 = dx*dx + dy*dy
    dist = sqrt(d2)
  end subroutine dist_xy_to_xy
  
  subroutine curvature_at_idx(px, py, poly, idx, dist, curvature)
    real(8), intent(in)  :: px, py
    real(8), intent(in)  :: poly(:, :)    ! (N,2)
    integer, intent(in) :: idx
    real(8), intent(out) :: dist, curvature

    integer :: n, im1, ip1, cidx
    real(8) :: dx, dy, d2
    real(8) :: dx1, dy1, ddx, ddy, denom
    real(8), parameter :: EPS = 1.0D-12

    n = size(poly,1)
    if (n < 3) then
      dist = 0.0D0
      curvature = 0.0D0
      return
    end if

    cidx = max(2, min(n-1, idx))

    im1 = cidx-1
    ip1 = cidx+1

    dx1 = poly(ip1,1) - poly(im1,1)
    dy1 = poly(ip1,2) - poly(im1,2)
    ddx = poly(ip1,1) - 2.0D0*poly(cidx,1) + poly(im1,1)
    ddy = poly(ip1,2) - 2.0D0*poly(cidx,2) + poly(im1,2)

    denom = (dx1*dx1 + dy1*dy1)**1.5D0
    if (denom <= 0.0D0) denom = EPS

    curvature = abs(dx1*ddy - dy1*ddx) / denom
  end subroutine curvature_at_idx

  ! Compute curvature and distance extremes from the poly and domain bounds.
  ! Called once before traversal so calc_stop_level can normalise properly.
  !
  !   kappa_min / kappa_max  — curvature range over all interior poly vertices
  !   dist_max               — distance from the farthest domain corner to the poly
  !
  subroutine poly_stats(poly, n, m, dist_max, kappa_min, kappa_max)
    real(8), intent(in)  :: poly(:,:)   ! (np, 2)
    integer, intent(in)  :: n, m
    real(8), intent(out) :: dist_max, kappa_min, kappa_max

    integer :: i, np, pidx
    real(8) :: d, kappa, kappa_dummy, dist_dummy
    real(8) :: corners(4, 2)

    np = size(poly, 1)

    ! --- curvature range from interior poly vertices ----------------------
    kappa_min =  1.0D300
    kappa_max = -1.0D300
    do i = 2, np - 1
      call curvature_at_idx(poly(i,1), poly(i,2), poly, i, dist_dummy, kappa)
      kappa_min = min(kappa_min, kappa)
      kappa_max = max(kappa_max, kappa)
    end do
    if (kappa_min ==  1.0D300) kappa_min = 0.0D0
    if (kappa_max == -1.0D300) kappa_max = 0.0D0

    ! --- dist_max from the four domain corners ----------------------------
    corners(1,:) = [0.0D0,    0.0D0   ]
    corners(2,:) = [real(m,8),0.0D0   ]
    corners(3,:) = [0.0D0,    real(n,8)]
    corners(4,:) = [real(m,8),real(n,8)]

    dist_max = 0.0D0
    do i = 1, 4
      call nearest_idx(corners(i,1), corners(i,2), poly, pidx)
      call dist_xy_to_xy(corners(i,1), corners(i,2), poly(pidx,1), poly(pidx,2), d)
      dist_max = max(dist_max, d)
    end do

  end subroutine poly_stats

end module mesh_utils

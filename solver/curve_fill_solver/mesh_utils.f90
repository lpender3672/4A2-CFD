module mesh_utils
  use iso_fortran_env, only: real64
  use types
  use general_utils
  implicit none
  private
  public :: lod_from_xy, interp_from_cell
  public :: naca4_airfoil, transform_airfoil
  public :: nearest_idx, dist_curvature_at_idx
  public :: calc_lod

contains

  integer function lod_from_xy(x, y, n, m, ILOD) result(lod)
      implicit none
      real(8), intent(in) :: x, y        ! coordinates in full-resolution domain [0, m), [0, n)
      integer, intent(in) :: n, m        ! fine domain size (for normalization)
      integer, intent(in) :: ILOD(:,:)   ! coarse integer LOD map (nL × mL)
      integer :: nL, mL
      real(8) :: rx, ry, gx, gy
      integer :: ix, iy
      real(8) :: f00, f10, f01, f11, wx, wy, interp

      ! size of coarse grid
      nL = size(ILOD, 1)
      mL = size(ILOD, 2)

      ! fractional position in coarse grid
      gx = x * real(mL - 1,8) / real(m,8)
      gy = y * real(nL - 1,8) / real(n,8)

      ! integer base indices (0-based)
      ix = int(floor(gx))
      iy = int(floor(gy))

      ! fractional parts
      wx = gx - real(ix,8)
      wy = gy - real(iy,8)

      ! clamp to valid interior range
      ix = max(0, min(mL-2, ix))
      iy = max(0, min(nL-2, iy))

      ! get corner values (convert to real for interpolation)
      f00 = real(ILOD(iy+1, ix+1), 8)
      f10 = real(ILOD(iy+1, ix+2), 8)
      f01 = real(ILOD(iy+2, ix+1), 8)
      f11 = real(ILOD(iy+2, ix+2), 8)

      ! bilinear interpolation
      interp = f00*(1.0-wx)*(1.0-wy) + f10*wx*(1.0-wy) + f01*(1.0-wx)*wy + f11*wx*wy

      ! round back to nearest integer
      lod = int(nint(interp))
    end function lod_from_xy


    real(8) function interp_from_cell(cell, map) result(phi)
    implicit none
    type(cell2d), intent(in) :: cell
    real(8), intent(in) :: map(:, :)
    integer :: n, m

    real(8) :: gx, gy
    real(8) :: wx, wy
    integer :: ix, iy
    real(8) :: d00, d10, d01, d11

    n = size(map, 1)
    m = size(map, 2)

    gx = 0.5d0 * (cell%xmin + cell%xmax)
    gy = 0.5d0 * (cell%ymin + cell%ymax)

    ix = int(floor(gx))
    iy = int(floor(gy))

    ix = max(0, min(m-2, ix))
    iy = max(0, min(n-2, iy))

    wx = gx - real(ix, 8)
    wy = gy - real(iy, 8)

    d00 = map(iy+1, ix+1)
    d10 = map(iy+1, ix+2)
    d01 = map(iy+2, ix+1)
    d11 = map(iy+2, ix+2)

    phi = d00*(1d0-wx)*(1d0-wy) + &
          d10*(wx)*(1d0-wy)     + &
          d01*(1d0-wx)*(wy)     + &
          d11*(wx)*(wy)

  end function interp_from_cell
  

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

  
  subroutine dist_curvature_at_idx(px, py, poly, idx, dist, curvature)
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

    dx = px - poly(idx,1)
    dy = py - poly(idx,2)
    d2 = dx*dx + dy*dy
    dist = sqrt(d2)

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
  end subroutine dist_curvature_at_idx

  subroutine calc_lod(n, m, foil, max_level, DIST, KAPPA, ILOD)
    integer, intent(in) :: n, m, max_level
    real(8), intent(in) :: foil(:, :)              ! (Nf,2) polyline
    integer, intent(out):: ILOD(n, m)
    real(8), allocatable, intent(out):: DIST(:,:), KAPPA(:, :)

    real(8) :: X(n, m), Y(n, m)
    real(8), allocatable :: xv(:), yv(:), LOD(:,:)
    real(8), allocatable :: NLOD(:,:)
    real(8) :: min_k_pos, max_dist, min_nlod, max_nlod
    integer :: i, j, idx

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
        call nearest_idx(X(i,j), Y(i,j),foil, idx)
        call dist_curvature_at_idx(X(i,j), Y(i,j), foil, idx, DIST(i,j), KAPPA(i,j))
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

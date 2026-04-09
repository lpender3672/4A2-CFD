module euler_flux
  use types, only: cell2d
  implicit none

  contains

  ! -----------------------------------------------------------------------
  ! Pressure from conserved variables: p = (gam-1)*(roe - 0.5*|rv|^2/ro)
  ! -----------------------------------------------------------------------
  pure real(8) function pressure(ro, rovx, rovy, roe, gam)
    real(8), intent(in) :: ro, rovx, rovy, roe, gam
    pressure = (gam - 1.0D0) * (roe - 0.5D0*(rovx**2 + rovy**2)/ro)
  end function

  ! -----------------------------------------------------------------------
  ! Physical Euler flux dotted with outward unit normal (nx, ny).
  ! Flux convention: positive = leaving the left state.
  ! -----------------------------------------------------------------------
  pure subroutine euler_flux_n(ro, rovx, rovy, roe, nx, ny, gam, &
                                f_ro, f_rovx, f_rovy, f_roe)
    real(8), intent(in)  :: ro, rovx, rovy, roe, nx, ny, gam
    real(8), intent(out) :: f_ro, f_rovx, f_rovy, f_roe
    real(8) :: p, un

    p  = pressure(ro, rovx, rovy, roe, gam)
    un = (rovx*nx + rovy*ny) / ro   ! normal velocity

    f_ro   = ro   * un
    f_rovx = rovx * un + p*nx
    f_rovy = rovy * un + p*ny
    f_roe  = (roe + p) * un
  end subroutine euler_flux_n

  ! -----------------------------------------------------------------------
  ! Maximum wave speed |u·n| + c at a state.  Used for CFL and dissipation.
  ! -----------------------------------------------------------------------
  pure real(8) function max_wave_speed(ro, rovx, rovy, roe, nx, ny, gam)
    real(8), intent(in) :: ro, rovx, rovy, roe, nx, ny, gam
    real(8) :: p, un, c

    p  = pressure(ro, rovx, rovy, roe, gam)
    un = (rovx*nx + rovy*ny) / ro
    c  = sqrt(max(gam*p/ro, 0.0D0))
    max_wave_speed = abs(un) + c
  end function

  ! -----------------------------------------------------------------------
  ! Lax-Friedrichs numerical flux between left and right states.
  ! (nx,ny) is the outward normal from L toward R.
  !
  !   F_LF = 0.5*(F_L + F_R) - 0.5*lambda*(Q_R - Q_L)
  !
  ! lambda = max wave speed across the face.  This is equivalent to the
  ! scalar dissipation the block-mesh solver adds via its smoothing pass,
  ! but applied locally per face.
  ! -----------------------------------------------------------------------
  pure subroutine lax_friedrichs(ro_L, rovx_L, rovy_L, roe_L, &
                                  ro_R, rovx_R, rovy_R, roe_R, &
                                  nx, ny, gam,                  &
                                  f_ro, f_rovx, f_rovy, f_roe)
    real(8), intent(in)  :: ro_L, rovx_L, rovy_L, roe_L
    real(8), intent(in)  :: ro_R, rovx_R, rovy_R, roe_R
    real(8), intent(in)  :: nx, ny, gam
    real(8), intent(out) :: f_ro, f_rovx, f_rovy, f_roe

    real(8) :: fL_ro, fL_rovx, fL_rovy, fL_roe
    real(8) :: fR_ro, fR_rovx, fR_rovy, fR_roe
    real(8) :: lam

    call euler_flux_n(ro_L, rovx_L, rovy_L, roe_L, nx, ny, gam, &
                      fL_ro, fL_rovx, fL_rovy, fL_roe)
    call euler_flux_n(ro_R, rovx_R, rovy_R, roe_R, nx, ny, gam, &
                      fR_ro, fR_rovx, fR_rovy, fR_roe)

    lam = max(max_wave_speed(ro_L, rovx_L, rovy_L, roe_L, nx, ny, gam), &
              max_wave_speed(ro_R, rovx_R, rovy_R, roe_R, nx, ny, gam))

    f_ro   = 0.5D0*(fL_ro   + fR_ro  ) - 0.5D0*lam*(ro_R   - ro_L  )
    f_rovx = 0.5D0*(fL_rovx + fR_rovx) - 0.5D0*lam*(rovx_R - rovx_L)
    f_rovy = 0.5D0*(fL_rovy + fR_rovy) - 0.5D0*lam*(rovy_R - rovy_L)
    f_roe  = 0.5D0*(fL_roe  + fR_roe ) - 0.5D0*lam*(roe_R  - roe_L )
  end subroutine lax_friedrichs

  ! -----------------------------------------------------------------------
  ! Shared face geometry for two axis-aligned cells.
  ! (nx,ny) is the outward unit normal from a toward b.
  ! face_len is the length of the shared edge (overlap interval).
  ! -----------------------------------------------------------------------
  pure subroutine face_geometry(a, b, face_len, nx, ny)
    type(cell2d), intent(in)  :: a, b
    real(8),      intent(out) :: face_len, nx, ny
    real(8), parameter :: tol = 1.0D-9

    if (abs(a%xmax - b%xmin) < tol) then           ! b right of a
      nx = 1.0D0;  ny = 0.0D0
      face_len = min(a%ymax, b%ymax) - max(a%ymin, b%ymin)
    else if (abs(a%xmin - b%xmax) < tol) then       ! b left of a
      nx = -1.0D0; ny = 0.0D0
      face_len = min(a%ymax, b%ymax) - max(a%ymin, b%ymin)
    else if (abs(a%ymax - b%ymin) < tol) then       ! b above a
      nx = 0.0D0;  ny = 1.0D0
      face_len = min(a%xmax, b%xmax) - max(a%xmin, b%xmin)
    else                                             ! b below a
      nx = 0.0D0;  ny = -1.0D0
      face_len = min(a%xmax, b%xmax) - max(a%xmin, b%xmin)
    end if
  end subroutine face_geometry

end module euler_flux

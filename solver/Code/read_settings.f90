      
      subroutine read_settings(fpath, av,bcs)

!     Read in the application variables and gas constants, they are in the
!     standard input file which is already assigned to unit 5

!     Explicitly declare the required variables
      use types
      implicit none
      character(len=*), intent(in) :: fpath
      character(len=:), allocatable :: folderpath
      type(t_appvars), intent(out) :: av
      type(t_bconds), intent(out) :: bcs
      character(len=80) :: tempname
      integer :: last_slash_index

!     Read the case name and trim to the required length
      open(unit=5, file=fpath, status='old')
      read(5,*) tempname
      av%casename = trim(tempname)

      last_slash_index = scan(fpath, '/', back=.true.)
      if (last_slash_index > 0) then
          folderpath = fpath(1:last_slash_index)
      else
          folderpath = ''
      end if

      av%casefolder = folderpath

!     You should read in the following variables sequentially and store them in
!     the dervived "av" datatype with the % syntax:
!         rgas, gam
!         cfl, sfac, d_max
!         nsteps
!         ni, nj
      read(5,*) av%rgas, av%gam
      read(5,*) av%cfl, av%sfac, av%d_max
      read(5,*) av%nsteps
      read(5,*) av%ni, av%nj

!     Calculate other gas constants used throughout the calculation
      av%cp = av%rgas * av%gam / (av%gam - 1.0)
      av%cv = av%cp / av%gam
      av%fgam = (av%gam - 1.0) / av%gam

!     Scale the smoothing factor and convergence limits by cfl number, the 
!     changes over a timestep should be proportional to dt, which is 
!     proportional to cfl
      av%sfac = av%sfac * av%cfl
      av%d_max = av%d_max * av%cfl

!     Average convergence limit on residuals is set to half of the maximum
      av%d_avg = 0.5 * av%d_max

!     Read the inlet boundary condition data and store into the "bcs" datatype
!         pstag, tstag, alpha, rfin
      read(5,*) bcs%pstag, bcs%tstag, bcs%alpha, bcs%rfin

!     Convert the inlet angle to radians
      bcs%alpha = bcs%alpha * 3.14159 / 180.0

!     Calculate the inlet stagnation density "rostag"
      bcs%rostag = bcs%pstag / (av%rgas * bcs%tstag)

!     Read the outlet static pressure and store into the "bcs" datatype
      read(5,*) bcs%p_out

      close(5)

!     Print the settings to check they have been read, you can use this syntax
!     anywhere else you want in the program to debug your code
      write(6,*)
      write(6,*) 'Solver begins on ', av%casename, ' case in folder ', av%casefolder
      write(6,*)
      write(6,*) 'Read application variables from file'
      write(6,*) '  rgas =', av%rgas, 'cp =', av%cp, 'cv =', av%cv
      write(6,*) '  CFL =', av%cfl, 'sfac =', av%sfac
      write(6,*) '  Convergence  d_max =', av%d_max
      write(6,*) '  Mesh size  ni =', av%ni, 'nj =', av%nj
      write(6,*) '  Inlet  pstag =', bcs%pstag, 'tstag =', bcs%tstag, &
          'alpha = ', bcs%alpha
      write(6,*) '  Outlet  p_out =', bcs%p_out
      write(6,*)

      end subroutine read_settings



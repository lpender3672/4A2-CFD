      
      subroutine read_settings(av,bcs)

!     Read in the application variables and gas constants, they are in the
!     standard input file which is already assigned to unit 5

!     Explicitly declare the required variables
      use types
      use io_module
      use routines
      implicit none
      character(len=:), allocatable :: folderpath, temp_casefolder, temp_casename, fpath
      type(t_appvars), intent(inout) :: av
      type(t_bconds), intent(inout) :: bcs
      character(len=80) :: tempname
      character(len=256) :: msg_bfr 
      integer :: last_slash_index

!     Read the case name and trim to the required length

      open(5, file = av%casefolder // '/input_' // av%casename // '.txt', status='old')
      ! open(5, file = 'cases/bump/input_bump.txt', status='old')
      read(5,*) tempname

      av%crashed = 0
      av%l_min = 1.0e10

!     You should read in the following variables sequentially and store them in
!     the dervived "av" datatype with the % syntax:
!         gam, rgas
!         cfl, sfac, d_max
!         nsteps
!         ni, nj
      read(5,*) av%rgas, av%gam
      read(5,*) av%cfl, av%sfac, av%sfac_res, av%d_max, av%d_var, av%facsec, av%fcorr
      read(5,*) av%nsteps, av%nrkuts, av%guess_method, av%tstep_method
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
      write(msg_bfr,*) 'Solver begins on ', trim(av%casename), ' case in folder ', trim(av%casefolder)
      call write_to_qt(msg_bfr)
      write(6,*)
      write(msg_bfr,*) 'Read application variables from file'
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  rgas =', av%rgas, 'cp =', av%cp, 'cv =', av%cv, 'gam =', av%gam
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  CFL =', av%cfl, 'sfac =', av%sfac, 'rfin =', bcs%rfin
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  sfac_res =', av%sfac_res, 'fcorr =', av%fcorr, 'facsec =', av%facsec
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  flow_guess_method =', av%guess_method, 'tstep method =', av%tstep_method, 'nrkuts = ', av%nrkuts
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  Convergence  d_max =', av%d_max
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  Mesh size  ni =', av%ni, 'nj =', av%nj
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  Inlet  pstag =', bcs%pstag, 'tstag =', bcs%tstag, 'alpha = ', bcs%alpha
      call write_to_qt(msg_bfr)
      write(msg_bfr,*) '  Outlet  p_out =', bcs%p_out
      call write_to_qt(msg_bfr)
      write(6,*)

      end subroutine read_settings



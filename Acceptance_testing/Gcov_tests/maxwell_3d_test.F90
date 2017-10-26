! ______________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! MAXWELL_3D_FDTD_TEST.F90
!
! Test Maxwell solver of PICSAR in 3D geometry by considering propagation of
! a quasi-monochromatic (Hanning temporal profile) gaussian laser in vacuum. 
! The test compares spatial amplitude profile to the gaussian theoretical solution. 
! It is performed for FDTD solvers parallelized with the standard cartesian domain 
! decomposition. 
!
! @author: Henri VINCENTI, 2017.10
! ______________________________________________________________________________

PROGRAM maxwell_3d_test
  USE constants
  USE fields
  USE particles
  USE params
  USE shared_data
  USE mpi_routines
  USE control_file
  USE time_stat
  USE diagnostics
#if defined(FFTW)
  USE mpi_fftw3
  USE fourier
  USE fastfft
  USE fftw3_fortran
#endif
  IMPLICIT NONE
  REAL(num) :: ttime, Emax, laser_w0, x0, y0, z0, lambda_laser, zslice1, zslice2, &
  tpeak, epsilon
  REAL(num), DIMENSION(3) :: err_vec1, err_vec2
  CHARACTER(len=64)            :: title 
  TYPE(particle_antenna), POINTER :: curr_laser
  LOGICAL(isp) :: passed

  ! --- default init
  CALL default_init

  ! --- reads input_file
  CALL read_input_file
  
  ! --- reads from command line 
  CALL read_from_cl

  ! --- mpi init communicator
#if defined(FFTW)
  IF (fftw_with_mpi) THEN 
    CALL mpi_minimal_init_fftw
  ELSE
#endif 
    CALL mpi_minimal_init
#if defined(FFTW)
  ENDIF 
#endif 
  IF (rank .EQ. 0) THEN
    write(0,*) "_________________________________________________________________"
    write(0,*) ""
    write(0,*) " TEST MAXWELL 3D"
    write(0,*) "_________________________________________________________________"
  ENDIF

  ! --- Check domain decomposition / Create Cartesian communicator / Allocate grid arrays
  CALL mpi_initialise

  ! --- allocates and inits particle distributions (on each subdomain)
  CALL initall

  ! --- Diagnostics
  CALL init_diags

  !----------------------------------------------
  ! ADVANCE MAXWELL EQUATIONS IN TIME 
  !----------------------------------------------
  IF (rank .EQ. 0) startsim=MPI_WTIME()
  CALL step(nsteps)

  IF (rank .EQ. 0) endsim=MPI_WTIME()

  !--- Compare computed values to theoretical value (Gaussian beam)
  err_vec1=0._num; err_vec2=0._num
  curr_laser=>species_parray(1)%antenna_params
  Emax=curr_laser%Emax
  IF (rank .EQ. 0) PRINT *, "LASER EMAX = ", Emax
  laser_w0=curr_laser%laser_w0
  x0=curr_laser%spot_x
  y0=curr_laser%spot_y
  z0=curr_laser%spot_z
  lambda_laser=curr_laser%lambda_laser
  tpeak=curr_laser%t_peak
  zslice1=dz/2.0_num
  zslice2=dz
  ttime=it*dt
  CALL check_amplitude_slice_error_loc(ttime,Emax,laser_w0,x0,y0,z0,lambda_laser, tpeak, &
  zslice1, err_vec1)
  CALL check_amplitude_slice_error_loc(ttime,Emax,laser_w0,x0,y0,z0,lambda_laser, tpeak, &
  zslice2, err_vec2)
  !--- Display test result on standard output 
  IF (l_spectral) THEN 
    IF (fftw_with_mpi) THEN 
      title='Maxwell solver PSAOTD - FFTW WITH MPI'
    ELSE
      title='Maxwell solver PSAOTD - FFTW WITH standard DD'
    ENDIF 
  ELSE
  	title='Maxwell solver Yee standard DD'
  ENDIF 
  CALL output_check(title,'Slice 1',err_vec1)
  CALL output_check(title,'Slice 2',err_vec2)
  
#if defined(FFTW)
  ! --- Destroy fftw plans
  IF(l_spectral) THEN
    IF(fftw_with_mpi) THEN
      CALL DFFTW_DESTROY_PLAN(plan_r2c_mpi)
      CALL DFFTW_DESTROY_PLAN(plan_c2r_mpi)
    ELSE
      CALL fast_fftw_destroy_plan_dft(plan_r2c)
      CALL fast_fftw_destroy_plan_dft(plan_c2r)
    ENDIF
  ENDIF
#endif

  ! --- MPI CLOSE 
  CALL mpi_close
  
  ! --- Decide test issue (success or failure)
  ! --- If failure, abort program 
  epsilon=1e-7
  passed=.TRUE.
  IF (err_vec1(1) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec1(2) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec1(3) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec2(1) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec2(2) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec2(3) .gt. epsilon) passed = (passed .AND. .FALSE.)

  write(0,*)
  PRINT *, "rank, error", rank, err_vec1, err_vec2
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST MAXWELL SOLVER 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST MAXWELL SOLVER 3D FAILED **********  \e[0m \n"')
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')  

END PROGRAM maxwell_3d_test

! ________________________________________________________________________________________
! External functions

! ________________________________________________________________
!
! Check the results by comparison with the original cases
!
! ________________________________________________________________
SUBROUTINE check_amplitude_slice_error_loc(ttime,Emax,w0,x0,y0,z0,lambda_laser,tpeak, &
zslice, err_vec)
  USE constants
  USE fields
  USE shared_data
  IMPLICIT NONE
  REAL(num), INTENT(IN) :: zslice,w0,x0,y0,z0,lambda_laser, ttime , Emax, tpeak
  REAL(num), DIMENSION(3), INTENT(OUT) :: err_vec
  REAL(num) :: wslice, Rslice, zetas, omega_laser, k_laser, zslice_g, zr0
  REAL(num), ALLOCATABLE, DIMENSION(:) :: x_ax_loc, y_ax_loc
  INTEGER(idp) :: islice, ix, iy
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: eth, eslice, diff

  
  ! -- Test if slice is on current CPU domain 
  IF ((zslice .LT. z_max_local) .AND. (zslice .GE. z_min_local)) THEN 
    ! -- Get field in the slice 
    ALLOCATE(eslice(nx,ny))
    islice=NINT((zslice-z_min_local)/dz,idp)
    zslice_g=islice*dz+z_min_local+dz/2.0_num
    eslice(1:nx,1:ny)=ex(0:nx-1,0:ny-1,islice)
    ! -- x and y axis local 
    ALLOCATE(x_ax_loc(nx), y_ax_loc(ny))
	DO ix = 1, nx
	  x_ax_loc(ix) = x_min_local-x0 + (ix-1) * dx
	ENDDO
	DO iy = 1, ny
	  y_ax_loc(iy) = y_min_local-y0 + (iy-1) * dy
	ENDDO
    ! -- Get theoretical field in the slice 
    ALLOCATE(eth(nx,ny))
    k_laser=2.0_num*pi/lambda_laser
    omega_laser=clight*k_laser
    zr0=pi*w0**2/lambda_laser
    wslice=w0*SQRT(1+((zslice_g-z0)/zr0)**2)
    Rslice=zslice_g*(1+(zr0/(zslice_g-z0))**2)
    zetas =ATAN((zslice_g-z0)/zr0)
    PRINT *, "zetas, x0, y0, z0, ATAN(0.0_num), zslice, zslice_g", zetas, x0, y0, z0, ATAN(0.0_num), zslice, zslice_g
	DO ix = 1, nx
	  DO iy = 1, ny
	      eth(ix,iy)= Emax*w0/wslice*EXP(-(x_ax_loc(ix)**2+y_ax_loc(iy)**2)/wslice**2)   &
	      *COS(zetas+omega_laser*(ttime-tpeak)-k_laser*zslice_g-k_laser*(x_ax_loc(ix)**2+&
	      y_ax_loc(iy)**2)/(2*Rslice))
	  END DO 
	END DO 
	PRINT *, "MAXVAL ETH, zslice", MAXVAL(eth), zslice
	PRINT *, "MAXVAL Eslice, zslice", MAXVAL(eslice), zslice
    ! -- Error point by point
    ALLOCATE(diff(nx,ny))
    diff = eslice - eth
    ! -- Average error on the slice 
    err_vec(1) = SUM(ABS(diff/eth),MASK=((eth .NE. 0._num) .AND. (eslice .NE. 0._num)))/ &
    MAX(1,COUNT(((eth .NE. 0._num) .AND. (eslice .NE. 0._num))))
    ! -- Min error on the slice 
    err_vec(2) = MAXVAL(ABS(diff/eth),MASK=((eth .NE. 0._num) .AND. (eslice .NE. 0._num)))
    ! -- Max error on the slice 
    err_vec(3) = MINVAL(ABS(diff/eth),MASK=((eth .NE. 0._num) .AND. (eslice .NE. 0._num)))
    DEALLOCATE(x_ax_loc, y_ax_loc, eth, eslice, diff)
  ELSE 
    err_vec(1) = 0._num 
    err_vec(2) = 0._num 
    err_vec(3) = 0._num 
  ENDIF   	
END SUBROUTINE
! ________________________________________________________________

! ________________________________________________________________
!
! This subroutine output the results of the checks
!
SUBROUTINE output_check(title,position,err_vec)
! ________________________________________________________________

  USE constants
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)            :: title, position
  REAL(num), DIMENSION (3)                  :: err_vec

  WRITE(0,*)
  WRITE(0,'(A40)') title
  WRITE(0,'(A40, 6(X,A13))') "Slice", "ave err", "min err", "max err"
  WRITE(0,'(" _____________________________________________________")')
  WRITE(0,'(A40,6(X,E13.5))') position, err_vec(1), err_vec(2), err_vec(3)

END SUBROUTINE
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
! MAXWELL_2D_FDTD_TEST.F90
!
! Test Maxwell solver of PICSAR in 3D geometry by considering propagation of
! a plane wave laser in vacuum. 
! The test compares spatial amplitude profile to the gaussian theoretical solution. 
! It is performed for FDTD solvers parallelized with the standard cartesian domain 
! decomposition. 
!
! @author: Haithem Kallala
! ______________________________________________________________________________

PROGRAM maxwell_2d_test
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
  REAL(num) :: ttime, Emax, laser_w0, x0, z0, lambda_laser, zslice1, zslice2, &
  tpeak, epsilon
  REAL(num), DIMENSION(3) :: err_vec1, err_vec2
  REAL(num) :: err_min, err_max, err_avg
  CHARACTER(len=64)            :: title 
  TYPE(particle_antenna), POINTER :: curr_laser
  LOGICAL(isp) :: passed, has_passed

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
    write(0,*) " TEST MAXWELL 2D PLANE WAVE NPROCS =", nproc
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

  !----------------------------------------------
  ! COMPARE THEORETICAL VALUES TO COMPUTED VALUES 
  !----------------------------------------------
  err_vec1=0._num; err_vec2=0._num
  curr_laser=>species_parray(1)%antenna_params
  Emax=curr_laser%Emax
  laser_w0=curr_laser%laser_w0
  x0=curr_laser%spot_x
  z0=curr_laser%spot_z
  lambda_laser=curr_laser%lambda_laser
  tpeak=curr_laser%t_peak
  zslice1=30._num*dz
  zslice2=50._num*dz
  ttime=it*dt
  CALL check_amplitude_slice_error_loc_2d(ttime,Emax,laser_w0,x0,z0,lambda_laser, tpeak, &
  zslice1, err_vec1)
  CALL check_amplitude_slice_error_loc_2d(ttime,Emax,laser_w0,x0,z0,lambda_laser, tpeak, &
  zslice2, err_vec2)
  !--- Display test result on standard output 
  IF (l_spectral) THEN 
    IF (fftw_with_mpi) THEN 
      title='Test Maxwell solver PSAOTD - FFTW WITH MPI'
    ELSE
      title='Test Maxwell solver PSAOTD - FFTW WITH standard DD'
    ENDIF 
  ELSE
  	title='Test Maxwell solver Yee standard DD'
  ENDIF 
  
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
  
  !----------------------------------------------
  ! DECIDE TEST ISSUE (if failure abort program)
  !----------------------------------------------
  epsilon=5e-2
  passed=.TRUE.
  IF (err_vec1(1) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec1(2) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec1(3) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec2(1) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec2(2) .gt. epsilon) passed = (passed .AND. .FALSE.)
  IF (err_vec2(3) .gt. epsilon) passed = (passed .AND. .FALSE.)

  ! REDUCE RESULTS TO PROC 0
  IF (rank .EQ. 0) CALL output_test_title(title)
  CALL MPI_REDUCE(err_vec1(1),err_min,1_isp,mpidbl,MPI_MIN,0_isp,comm,errcode)
  CALL MPI_REDUCE(err_vec1(2),err_max,1_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)
  CALL MPI_REDUCE(err_vec1(3)/nproc,err_avg,1_isp,mpidbl,MPI_SUM,0_isp,comm,errcode)
  IF (rank .EQ. 0) CALL output_check('Slice 1',(/err_min, err_max, err_avg/))
  CALL MPI_REDUCE(err_vec2(1),err_min,1_isp,mpidbl,MPI_MIN,0_isp,comm,errcode)
  CALL MPI_REDUCE(err_vec2(2),err_max,1_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)
  CALL MPI_REDUCE(err_vec2(3)/nproc,err_avg,1_isp,mpidbl,MPI_SUM,0_isp,comm,errcode)
  IF (rank .EQ. 0) CALL output_check('Slice 2',(/err_min, err_max, err_avg/))

  ! REDUCED PASSED RESULTS 
  CALL MPI_REDUCE(passed,has_passed,1_isp,MPI_LOGICAL,MPI_LAND,0_isp,comm,errcode)

  write(0,*)
  IF (rank .EQ. 0) THEN 
    IF (passed) THEN
      CALL system('printf "\e[32m ********** TEST MAXWELL SOLVER 3D PASSED **********  & 
      \e[0m \n"')
    ELSE
      CALL system('printf "\e[31m ********** TEST MAXWELL SOLVER 3D FAILED **********  & 
      \e[0m \n"')
      CALL EXIT(9)
    ENDIF
    write(0,'(" ______________________________________________________________________")')  
  ENDIF 

  ! --- MPI CLOSE 
  CALL mpi_close

END PROGRAM maxwell_2d_test

! ________________________________________________________________________________________
! External functions

! _______________________________________________________________________
! Routine that check the results by comparison with theoretical solution 
! _______________________________________________________________________
SUBROUTINE check_amplitude_slice_error_loc_2d(ttime,Emax,w0,x0,z0,lambda_laser,tpeak, &
zslice, err_vec)
  USE constants
  USE fields
  USE shared_data
  IMPLICIT NONE
  REAL(num), INTENT(IN) :: zslice,w0,x0,z0,lambda_laser, ttime , Emax, tpeak
  REAL(num), DIMENSION(3), INTENT(OUT) :: err_vec
  REAL(num) :: wslice, Rslice, zetas, omega_laser, k_laser, zslice_g, zr0
  REAL(num), ALLOCATABLE, DIMENSION(:) :: x_ax_loc
  INTEGER(idp) :: islice, ix, iy
  REAL(num), ALLOCATABLE, DIMENSION(:) :: eth, eslice, diff

  ! -- Test if slice is on current CPU domain 
  IF ((zslice .LT. z_max_local) .AND. (zslice .GE. z_min_local)) THEN 
    ! -- Get field in the slice 
    ALLOCATE(eslice(nx))
    islice=NINT((zslice-z_min_local)/dz,idp)
    zslice_g=islice*dz+z_min_local+dz/2.0_num
    eslice(1:nx)=ex(0:nx-1,0,islice)
    ! -- x and y axis local 
    ALLOCATE(x_ax_loc(nx))
    DO ix = 1, nx
      x_ax_loc(ix) = x_min_local-x0 + (ix-1) * dx
    ENDDO
    ! -- Get theoretical field in the slice (here w0>>lambda_laser to simulate plane wave)
    ALLOCATE(eth(nx))
    k_laser=2.0_num*pi/lambda_laser
    omega_laser=clight*k_laser
    zr0=pi*w0**2/lambda_laser
    wslice=w0*SQRT(1+((zslice_g-z0)/zr0)**2)
    Rslice=(zslice_g-z0)*(1+(zr0/(zslice_g-z0))**2)
    zetas =ATAN((zslice_g-z0)/zr0)
	DO ix = 1, nx
	      eth(ix)= Emax*w0/wslice*EXP(-(x_ax_loc(ix)**2)/wslice**2)                      &
	      *(-1.0_num)*SIN(zetas+omega_laser*(ttime-tpeak)-k_laser*(zslice_g-z0)      -   &
	      k_laser*(x_ax_loc(ix)**2/(2*Rslice)))
	END DO 
    ! -- Error point by point
    ALLOCATE(diff(nx))
    diff = eslice - eth
    ! -- Average error on the slice 
    err_vec(1) = SUM(ABS(diff/eth),MASK=((eth .NE. 0._num) .AND. (eslice .NE. 0._num)))/ &
    MAX(1,COUNT(((eth .NE. 0._num) .AND. (eslice .NE. 0._num))))
    ! -- Min error on the slice 
    err_vec(2) = MINVAL(ABS(diff/eth),MASK=((eth .NE. 0._num) .AND. (eslice .NE. 0._num)))
    ! -- Max error on the slice 
    err_vec(3) = MAXVAL(ABS(diff/eth),MASK=((eth .NE. 0._num) .AND. (eslice .NE. 0._num)))
    DEALLOCATE(x_ax_loc, eth, eslice, diff)
  ELSE 
    err_vec(1) = 0._num 
    err_vec(2) = 0._num 
    err_vec(3) = 0._num 
  ENDIF   	
END SUBROUTINE check_amplitude_slice_error_loc_2d
! ________________________________________________________________

! ________________________________________________________________
!
! Write test title 
! ________________________________________________________________
SUBROUTINE output_test_title(title)
  USE constants
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)            :: title
  
  WRITE(0,*)
  WRITE(0,'(A40)') title
  WRITE(0,'(A40, 6(X,A13))') "Slice", "ave err", "min err", "max err"
  WRITE(0,'(" _____________________________________________________")')
END SUBROUTINE

! ________________________________________________________________
!
! This subroutine output the results of the checks
! ________________________________________________________________
SUBROUTINE output_check(position,err_vec)
  USE constants
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)            :: position
  REAL(num), DIMENSION (3)                  :: err_vec

  WRITE(0,*)
  WRITE(0,'(A40,6(X,E13.5))') position, err_vec(1), err_vec(2), err_vec(3)
END SUBROUTINE


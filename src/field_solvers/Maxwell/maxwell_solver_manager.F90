! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
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
! MAXWELL_SOLVER_MANAGER.F90
!
! Purpose:
! This file contains subroutines for the Maxwell solvers.
!
! Authors:
! Henri Vincenti
! Mathieu Lobet
!
! Date:
! Creation 2015
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> PUSH B field half a time step
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE push_bfield
  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2).AND.(norderz.eq.2)) then
    CALL pxrpush_em3d_bvec(ex, ey, ez, bx, by, bz, 0.5_num*dt/dx, 0.5_num*dt/dy,      &
    0.5_num*dt/dz, nx, ny, nz, nxguards, nyguards, nzguards, nxs, nys, nzs,           &
    l_nodalgrid)
    ! Yee scheme arbitrary order
  ELSE
    CALL pxrpush_em3d_bvec_norder(ex, ey, ez, bx, by, bz, 0.5_num*dt/dx*xcoeffs,      &
    0.5_num*dt/dy*ycoeffs, 0.5_num*dt/dz*zcoeffs, nx, ny, nz, norderx, nordery,       &
    norderz, nxguards, nyguards, nzguards, nxs, nys, nzs, l_nodalgrid)
  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(5) = localtimes(5) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE push_bfield
! ________________________________________________________________________________________
!> @brief
!> Computes EM energy
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________
SUBROUTINE compute_em_energy
  USE shared_data
  USE constants
  USE fields
  USE params
  USE mpi
  IMPLICIT NONE

  electro_energy_mpi = 0.0_num
  magnetic_energy_mpi = 0.0_num
  electromagn_energy_mpi = 0.0_num

  electro_energy_mpi = SUM(ABS(ex(0:nx-1, 0:ny-1, 0:nz-1)**2))*dx*dy*dz
  electro_energy_mpi = electro_energy_mpi + SUM(ABS(ey(0:nx-1, 0:ny-1,                &
  0:nz-1)**2))*dx*dy*dz
  electro_energy_mpi = electro_energy_mpi + SUM(ABS(ez(0:nx-1, 0:ny-1,                &
  0:nz-1)**2))*dx*dy*dz
  electro_energy_mpi = electro_energy_mpi*0.5_num*eps0

  magnetic_energy_mpi = SUM(ABS(bx(0:nx-1, 0:ny-1, 0:nz-1)**2))*dx*dy*dz
  magnetic_energy_mpi = magnetic_energy_mpi+ SUM(ABS(by(0:nx-1, 0:ny-1,               &
  0:nz-1)**2))*dx*dy*dz
  magnetic_energy_mpi = magnetic_energy_mpi+ SUM(ABS(bz(0:nx-1, 0:ny-1,               &
  0:nz-1)**2))*dx*dy*dz
  magnetic_energy_mpi = magnetic_energy_mpi*0.5_num/mu0

  electromagn_energy_mpi =  magnetic_energy_mpi + electro_energy_mpi

  electro_energy_total = 0.0_num
  magneto_energy_total = 0.0_num
  electromagn_energy_total = 0.0_num

  CALL MPI_ALLREDUCE(electro_energy_mpi, electro_energy_total, 1_isp, MPI_DOUBLE,     &
  MPI_SUM, comm, errcode)
  CALL MPI_ALLREDUCE(magnetic_energy_mpi, magneto_energy_total, 1_isp, MPI_DOUBLE,    &
  MPI_SUM, comm, errcode)
  electromagn_energy_total= electro_energy_total+magneto_energy_total
END SUBROUTINE


! ________________________________________________________________________________________
!> @brief
!> Field damping in pml region
!> Damps fields in pml region
!> if vaccum then does nothing
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________

SUBROUTINE field_damping_bcs
  USE fields
  USE shared_data
  USE constants
  USE omp_lib
  USE time_stat
  INTEGER(idp)  :: ix,iy,iz
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
  DO ix = -nxguards,nx+nxguards-1
    DO iy = -nyguards,ny+nyguards-1
      DO iz = -nzguards,nz+nzguards-1
        exy(ix,iy,iz) = sigma_y_e(iy) *exy(ix,iy,iz)
        exz(ix,iy,iz) = sigma_z_e(iz) *exz(ix,iy,iz)
        eyx(ix,iy,iz) = sigma_x_e(ix) *eyx(ix,iy,iz)
        eyz(ix,iy,iz) = sigma_z_e(iz) *eyz(ix,iy,iz)
        ezx(ix,iy,iz) = sigma_x_e(ix) *ezx(ix,iy,iz)
        ezy(ix,iy,iz) = sigma_y_e(iy) *ezy(ix,iy,iz)
        bxy(ix,iy,iz) = sigma_y_b(iy) *bxy(ix,iy,iz)
        bxz(ix,iy,iz) = sigma_z_b(iz) *bxz(ix,iy,iz)
        byx(ix,iy,iz) = sigma_x_b(ix) *byx(ix,iy,iz)
        byz(ix,iy,iz) = sigma_z_b(iz) *byz(ix,iy,iz)
        bzx(ix,iy,iz) = sigma_x_b(ix) *bzx(ix,iy,iz)
        bzy(ix,iy,iz) = sigma_y_b(iy) *bzy(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END subroutine field_damping_bcs


! ________________________________________________________________________________________
!> @brief
!> Mege splitted fields of when using absorbing_bcs to compute real EM field
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________

SUBROUTINE merge_fields()
  USE fields
  USE shared_data
  USE omp_lib
  USE time_stat
 
  INTEGER(idp)  :: ix,iy,iz
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF


  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
  DO ix = -nxguards,nx+nxguards
    DO iy = -nyguards,ny+nyguards
      DO iz = -nzguards,nz+nzguards
        ex(ix,iy,iz) = exy(ix,iy,iz) + exz(ix,iy,iz)
        ey(ix,iy,iz) = eyx(ix,iy,iz) + eyz(ix,iy,iz)
        ez(ix,iy,iz) = ezx(ix,iy,iz) + ezy(ix,iy,iz)
        bx(ix,iy,iz) = bxy(ix,iy,iz) + bxz(ix,iy,iz)
        by(ix,iy,iz) = byx(ix,iy,iz) + byz(ix,iy,iz)
        bz(ix,iy,iz) = bzx(ix,iy,iz) + bzy(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END subroutine merge_fields
   
! ________________________________________________________________________________________
!> @brief
!> PUSH E field a full  time step
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE push_efield
  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2).AND.(norderz.eq.2)) then
    CALL pxrpush_em3d_evec(ex, ey, ez, bx, by, bz, jx, jy, jz, clight**2*mu0*dt,      &
    clight**2*dt/dx, clight**2*dt/dy, clight**2*dt/dz, nx, ny, nz, nxguards,          &
    nyguards, nzguards, nxs, nys, nzs, l_nodalgrid)

  ELSE
    ! Yee scheme arbitrary order
    CALL pxrpush_em3d_evec_norder(ex, ey, ez, bx, by, bz, jx, jy, jz,                 &
    clight**2*mu0*dt, clight**2*dt/dx*xcoeffs, clight**2*dt/dy*ycoeffs,               &
    clight**2*dt/dz*zcoeffs, nx, ny, nz, norderx, nordery, norderz, nxguards,         &
    nyguards, nzguards, nxs, nys, nzs, l_nodalgrid)
  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(7) = localtimes(7) + (MPI_WTIME() - tmptime)
  ENDIF
END SUBROUTINE push_efield


! ________________________________________________________________________________________
!> @brief
!> PUSH B field half a time step in 2D
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE push_bfield_2d
  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(norderz.eq.2)) then

    CALL pxrpush_em2d_bvec(ex, ey, ez, bx, by, bz, 0.5_num*dt/dx, 0._num,             &
    0.5_num*dt/dz, nx, ny, nz, nxguards, 0_idp, nzguards, nxs, 0_idp, nzs,         &
    l_nodalgrid)

    ! Yee scheme arbitrary order
  ELSE

    CALL pxrpush_em2d_bvec_norder(ex, ey, ez, bx, by, bz, 0.5_num*dt/dx*xcoeffs,      &
    0.5_num*dt/dy*ycoeffs, 0.5_num*dt/dz*zcoeffs, nx, ny, nz, norderx, nordery,       &
    norderz, nxguards, nyguards, nzguards, nxs, nys, nzs, l_nodalgrid)

  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(5) = localtimes(5) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE push_bfield_2d


! ________________________________________________________________________________________
!> @brief
!> PUSH B field half a time step in 2D
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE push_efield_2d
  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime,mdt
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF
  mdt = mu0*clight**2*dt
  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(norderz.eq.2)) then
    CALL pxrpush_em2d_evec(ex, ey, ez, bx, by, bz,jx,jy,jz,mdt, clight**2*dt/dx,clight**2*dt/dy, &
    clight**2*dt/dz, nx,ny,nz, nxguards, nyguards, nzguards,nxs,0_idp,nzs,&
    l_nodalgrid)

    ! Yee scheme arbitrary order
  ELSE

    CALL pxrpush_em2d_evec_norder(ex, ey, ez, bx, by, bz,jx,jy,jz,mdt, clight**2*dt/dx*xcoeffs,&
    clight**2*dt/dy*ycoeffs, clight**2*dt/dz*zcoeffs, nx, ny, nz, norderx, nordery,&
    norderz, nxguards, nyguards, nzguards, nxs, 0_idp, nzs, l_nodalgrid)

  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(7) = localtimes(7) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE push_efield_2d


! ________________________________________________________________________________________
!> @brief
!> PUSH E, B PSAOTD a full time step
!> This subroutine pushes the electric and the magnetic fields using
!> the PSATD solver.
!
!> @details
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation March 29 2017
! ________________________________________________________________________________________
SUBROUTINE push_psatd_ebfield_3d() bind(C, name='push_psatd_ebfield_3d_')
  USE constants
  USE time_stat
  USE params
  USE shared_data
  USE fields
#if defined(FFTW)
  USE fourier_psaotd
  USE matrix_coefficients
#endif
  IMPLICIT NONE

  REAL(num) :: tmptime, tmptime_m

#if defined(DEBUG)
  WRITE(0, *) "push psatd ebfield 3d: start"
#endif

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

#if defined(FFTW)
  ! - Fourier Transform R2C
  IF (fftw_with_mpi) THEN
    IF(fftw_hybrid) THEN
      CALL get_Ffields_mpi_lb!  -global-hybrid balanced FFT
    ELSE
      CALL get_Ffields_mpi! - global-hybrid FFT
    ENDIF
  ELSE
    CALL get_Ffields! - local FFT
  ENDIF
  IF(g_spectral) THEN
    IF (it.ge.timestat_itstart) THEN
      tmptime_m = MPI_WTIME()
    ENDIF
    CALL multiply_mat_vector(nmatrixes)
    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime_m)
    ENDIF
  ELSE 
    CALL push_psaotd_ebfielfs! - PUSH PSATD
  ENDIF
  ! - Inverse Fourier Transform C2R
  IF (fftw_with_mpi) THEN
    IF(fftw_hybrid) THEN
      CALL get_fields_mpi_lb! -global-hybrid balanced IFFT
    ELSE
      CALL get_fields_mpi! global-hybrid IFFT
    ENDIF
  ELSE
    CALL get_fields! local IFFT
  ENDIF
#endif
  IF (it.ge.timestat_itstart) THEN
    localtimes(24) = localtimes(24) + (MPI_WTIME() - tmptime)
  ENDIF

#if defined(DEBUG)
  WRITE(0, *) "push psatd ebfield 3d: end"
#endif

END SUBROUTINE

!> @brief
!> PUSH E, B PSAOTD a full time step
!> This subroutine pushes the electric and the magnetic fields using
!> the PSATD solver.
!
!> @details
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation March 29 2017
! ________________________________________________________________________________________
SUBROUTINE push_psatd_ebfield_2d() bind(C, name='push_psatd_ebfield_2d_')
  USE constants
  USE time_stat
  USE params
  USE shared_data
  USE fields
#if defined(FFTW)
  USE fourier_psaotd
  USE matrix_coefficients
#endif
  IMPLICIT NONE

  REAL(num) :: tmptime, tmptime_m

#if defined(DEBUG)
  WRITE(0, *) "push psatd ebfield 2d: start"
#endif

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF
#if defined(FFTW)
  ! - Fourier Transform R2C
  IF (fftw_with_mpi) THEN
    IF(fftw_hybrid) THEN
      CALL get_Ffields_mpi_lb ! -global-hybrid balanced FFT
    ELSE
      CALL get_Ffields_mpi! - global-hybrid FFT
    ENDIF

  ELSE
    CALL get_Ffields! - local FFT
  ENDIF

  IF(g_spectral) THEN
    IF (it.ge.timestat_itstart) THEN
      tmptime_m = MPI_WTIME()
    ENDIF

    CALL multiply_mat_vector(nmatrixes)

    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime_m)
    ENDIF

  ELSE
    CALL push_psaotd_ebfielfs_2d! - PUSH PSATD
  ENDIF
  ! - Inverse Fourier Transform C2R
  IF (fftw_with_mpi) THEN
    IF(fftw_hybrid) THEN
      CALL get_fields_mpi_lb! global-hybrid balanced IFFT
    ELSE
      CALL get_fields_mpi! global-hybrid IFFT
    ENDIF
  ELSE
    CALL get_fields! local IFFT
  ENDIF
#endif
  IF (it.ge.timestat_itstart) THEN
    localtimes(24) = localtimes(24) + (MPI_WTIME() - tmptime)
  ENDIF

#if defined(DEBUG)
  WRITE(0, *) "push psatd ebfield 2d: end"
#endif

END SUBROUTINE



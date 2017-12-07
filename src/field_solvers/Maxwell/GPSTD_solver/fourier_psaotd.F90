! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016,
! The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a paid-up,
! nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
! copies to the public, prepare derivative works, and perform publicly and display
! publicly, and to permit other to do so.
!
!
! fourier_psaotd.F90
!
! Purpose:
! This file contains subroutines for:
! (i)  Fourier init (allocation and init of Fourier k-vectors,
! creation of planes for FFTW etc.),
! (ii) Forward/Backward Fourier transform of EM fields using,
! FFTW (Distributed/Shared version),
! (iii) Maxwell push in the spectral space using the Pseudo-Spectral Arbitrary Order
! Analytical Time Domain (PSAOTD) Maxwell solver.
!
! Authors:
! Henri Vincenti
! Jean-Luc Vay
!
! Date:
! Creation March, 2017
!
! ________________________________________________________________________________________

! ---

MODULE fourier_psaotd
  USE gpstd_solver
  USE matrix_coefficients

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE init_plans_fourier_mpi(nopenmp)
    USE PICSAR_precision
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE group_parameters
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nopenmp
    INTEGER(C_INT) :: nopenmp_cint, iret
    INTEGER(C_INTPTR_T) :: nx_cint, ny_cint, nz_cint
    INTEGER(idp)        :: i
    INTEGER(isp)        :: planner_flag_1, planner_flag_2
    nopenmp_cint=nopenmp

    IF  (fftw_threads_ok) THEN
      CALL  DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    ENDIF
    planner_flag_1 = FFTW_MEASURE
    planner_flag_2 = FFTW_MEASURE
    IF(fftw_mpi_transpose) THEN
      planner_flag_1 = FFTW_MPI_TRANSPOSED_OUT
      planner_flag_2 = FFTW_MPI_TRANSPOSED_IN
    ENDIF
    IF(.NOT. fftw_hybrid) THEN
      nz_cint=nz_global
      ny_cint=ny_global
      nx_cint=nx_global
      IF(c_dim == 3) THEN 
        plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint, ny_cint, nx_cint, ex_r, exf,   &
        comm, planner_flag_1)
        plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(nz_cint, ny_cint, nx_cint, exf, ex_r,   &
        comm, planner_flag_2)
      ELSE IF(c_dim == 2) THEN
        plan_r2c_mpi = fftw_mpi_plan_dft_r2c_2d(nz_cint, nx_cint,ex_r, exf,   &
        comm, planner_flag_1)
        plan_c2r_mpi = fftw_mpi_plan_dft_c2r_2d(nz_cint,nx_cint, exf,ex_r,   &
        comm, planner_flag_2)
      ENDIF
    ELSE
      nz_cint = nz_group
      ny_cint = ny_group
      nx_cint = nx_group
      DO i=1, nb_group
        IF(MPI_COMM_GROUP_ID(i) .NE. MPI_COMM_NULL) THEN
          IF(c_dim == 3) THEN
            plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint, ny_cint, nx_cint, ex_r,    &
            exf, MPI_COMM_GROUP_ID(i), planner_flag_1)
            plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(nz_cint, ny_cint, nx_cint, exf,     &
            ex_r, MPI_COMM_GROUP_ID(i), planner_flag_2)
          ELSE IF(c_dim == 2) THEN
            plan_r2c_mpi = fftw_mpi_plan_dft_r2c_2d(nz_cint, nx_cint,ex_r,    &
            exf, MPI_COMM_GROUP_ID(i), planner_flag_1)
            plan_c2r_mpi = fftw_mpi_plan_dft_c2r_2d(nz_cint, nx_cint,exf,     &
            ex_r, MPI_COMM_GROUP_ID(i), planner_flag_2)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE init_plans_fourier_mpi

  SUBROUTINE get_Ffields()
    USE shared_data
    USE fields
    USE fourier
    USE fastfft
    USE time_stat
    USE params

    IMPLICIT NONE
    INTEGER(idp) :: nfftx, nffty, nfftz, nxx, nyy, nzz
    REAL(num)    :: tmptime
    nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
    nfftz=nz+2*nzguards+1
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nz+2*nzguards
#endif
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Init fourier fields fields
#if !defined(LIBRARY)
    CALL copy_field(ex_r, nfftx, nffty, nfftz, ex, nxx, nyy, nzz)
    CALL copy_field(ey_r, nfftx, nffty, nfftz, ey, nxx, nyy, nzz)
    CALL copy_field(ez_r, nfftx, nffty, nfftz, ez, nxx, nyy, nzz)
    CALL copy_field(bx_r, nfftx, nffty, nfftz, bx, nxx, nyy, nzz)
    CALL copy_field(by_r, nfftx, nffty, nfftz, by, nxx, nyy, nzz)
    CALL copy_field(bz_r, nfftx, nffty, nfftz, bz, nxx, nyy, nzz)
    CALL copy_field(jx_r, nfftx, nffty, nfftz, jx, nxx, nyy, nzz)
    CALL copy_field(jy_r, nfftx, nffty, nfftz, jy, nxx, nyy, nzz)
    CALL copy_field(jz_r, nfftx, nffty, nfftz, jz, nxx, nyy, nzz)
    CALL copy_field(rho_r, nfftx, nffty, nfftz, rho, nxx, nyy, nzz)
    CALL copy_field(rhoold_r, nfftx, nffty, nfftz, rhoold, nxx, nyy, nzz)
#endif
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    ! Do Fourier Transform
    CALL fft_forward_r2c_local(nfftx,nffty,nfftz)

  END SUBROUTINE get_Ffields
 
  SUBROUTINE get_Ffields_mpi_lb()
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    USE field_boundary

    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz
    REAL(num)    :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
jy_r=0.
call mpi_barrier(comm,errcode)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=1,sizes_to_exchange_f_to_recv(z_coords+1)
      DO iy=iy_min_r,iy_max_r
        DO ix =ix_min_r,ix_max_r
      !     ex_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
      !           ex(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     ey_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
      !           ey(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     ez_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
      !           ez(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     bx_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
      !           bx(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     by_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
      !           by(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     bz_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =& 
      !          bz(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     jx_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =& 
      !          jx(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
           jy_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
                 jy(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     jz_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =& 
      !          jz(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     rho_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) = &
      !           rho(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
      !     rhoold_r(ix,iy,iz-1+f_first_cell_to_recv(z_coords+1)) =&
      !           rhoold(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_send(z_coords+1))
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    CALL generalized_comms_group_r2f()
call mpi_barrier(comm,errcode)
    ! Get global Fourier transform of all fields components and currents
!    CALL fft_forward_r2c_mpi() 
    CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r, jyf)


  END SUBROUTINE get_Ffields_mpi_lb 

  SUBROUTINE get_Ffields_mpi
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    USE field_boundary
    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz
    REAL(num)    :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    ! Copy array values before FFT
    IF(.NOT. fftw_hybrid) THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=iz_min_r, iz_max_r
        DO iy=iy_min_r, iy_max_r
          DO ix=ix_min_r, ix_max_r
            ex_r(ix, iy, iz)=ex(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            ey_r(ix, iy, iz)=ey(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            ez_r(ix, iy, iz)=ez(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            bx_r(ix, iy, iz)=bx(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            by_r(ix, iy, iz)=by(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            bz_r(ix, iy, iz)=bz(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            jx_r(ix, iy, iz)=jx(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            jy_r(ix, iy, iz)=jy(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            jz_r(ix, iy, iz)=jz(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            rho_r(ix, iy, iz)=rho(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
            rhoold_r(ix, iy, iz)=rhoold(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ELSE
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=iz_min_r, iz_max_r
        DO iy=iy_min_r, iy_max_r
          DO ix=ix_min_r, ix_max_r
            ex_r(ix, iy, iz)=ex(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            ey_r(ix, iy, iz)=ey(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            ez_r(ix, iy, iz)=ez(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            bx_r(ix, iy, iz)=bx(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            by_r(ix, iy, iz)=by(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            bz_r(ix, iy, iz)=bz(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            jx_r(ix, iy, iz)=jx(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            jy_r(ix, iy, iz)=jy(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            jz_r(ix, iy, iz)=jz(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,           &
            iz-iz_min_r)
            rho_r(ix, iy, iz)=rho(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,         &
            iz-iz_min_r)
            rhoold_r(ix, iy, iz)=rhoold(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards,   &
            iz-iz_min_r)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    IF(fftw_hybrid) THEN
      CALL ebj_field_bcs_groups()
    ENDIF
    ! Get global Fourier transform of all fields components and currents
    CALL fft_forward_r2c_mpi

  END SUBROUTINE get_Ffields_mpi


  SUBROUTINE get_fields()
    USE params
    USE shared_data
    USE fields
    USE fourier
    USE fastfft
    USE time_stat
    IMPLICIT NONE
    REAL(num) :: coeff_norm, tmptime
    INTEGER(idp) :: ix, iy, iz, nxx, nyy, nzz, nfftx, nffty, nfftz
    nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
    nfftz=nz+2*nzguards+1
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nz+2*nzguards
#endif
    coeff_norm=1.0_num/(nfftx*nffty*nfftz)

    ! Get Inverse Fourier transform of all fields components 
    CALL fft_backward_c2r_local(nfftx,nffty,nfftz)

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
#if defined (LIBRARY)

    CALL normalize_Fourier(ex_r, nfftx, nffty, nfftz, ex_r, nfftx, nffty, nfftz,      &
    coeff_norm)
    CALL normalize_Fourier(ey_r, nfftx, nffty, nfftz, ey_r, nfftx, nffty, nfftz,      &
    coeff_norm)
    CALL normalize_Fourier(ez_r, nfftx, nffty, nfftz, ez_r, nfftx, nffty, nfftz,      &
    coeff_norm)
    CALL normalize_Fourier(bx_r, nfftx, nffty, nfftz, bx_r, nfftx, nffty, nfftz,      &
    coeff_norm)
    CALL normalize_Fourier(by_r, nfftx, nffty, nfftz, by_r, nfftx, nffty, nfftz,      &
    coeff_norm)
    CALL normalize_Fourier(bz_r, nfftx, nffty, nfftz, bz_r, nfftx, nffty, nfftz,      &
    coeff_norm)
#else
    CALL normalize_Fourier(ex, nxx, nyy, nzz, ex_r, nfftx, nffty, nfftz, coeff_norm)
    CALL normalize_Fourier(ey, nxx, nyy, nzz, ey_r, nfftx, nffty, nfftz, coeff_norm)
    CALL normalize_Fourier(ez, nxx, nyy, nzz, ez_r, nfftx, nffty, nfftz, coeff_norm)
    CALL normalize_Fourier(bx, nxx, nyy, nzz, bx_r, nfftx, nffty, nfftz, coeff_norm)
    CALL normalize_Fourier(by, nxx, nyy, nzz, by_r, nfftx, nffty, nfftz, coeff_norm)
    CALL normalize_Fourier(bz, nxx, nyy, nzz, bz_r, nfftx, nffty, nfftz, coeff_norm)
#endif
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE get_fields


  SUBROUTINE get_fields_mpi
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    IMPLICIT NONE
    REAL(num) :: coeff_norm, tmptime
    INTEGER(idp) :: ix, iy, iz

    ! Get global Fourier transform of all fields components 
    CALL fft_forward_c2r_mpi

    IF(fftw_hybrid) THEN 
      coeff_norm  = 1.0_num/(nx_group*ny_group*nz_group)
    ELSE 
       coeff_norm = 1.0_num/(nx_global*ny_global*nz_global)
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    IF(.NOT. fftw_hybrid) THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=iz_min_r, iz_max_r
        DO iy=iy_min_r, iy_max_r
          DO ix=ix_min_r, ix_max_r
            ex(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=ex_r(ix, iy, iz)*coeff_norm
            ey(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=ey_r(ix, iy, iz)*coeff_norm
            ez(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=ez_r(ix, iy, iz)*coeff_norm
            bx(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=bx_r(ix, iy, iz)*coeff_norm
            by(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=by_r(ix, iy, iz)*coeff_norm
            bz(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=bz_r(ix, iy, iz)*coeff_norm
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ELSE
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=iz_min_r, iz_max_r
        DO iy=iy_min_r, iy_max_r
          DO ix=ix_min_r, ix_max_r
            ex(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards, iz-iz_min_r) = ex_r(ix,    &
            iy, iz)*coeff_norm
            ey(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards, iz-iz_min_r) = ey_r(ix,    &
            iy, iz)*coeff_norm
            ez(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards, iz-iz_min_r) = ez_r(ix,    &
            iy, iz)*coeff_norm
            bx(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards, iz-iz_min_r) = bx_r(ix,    &
            iy, iz)*coeff_norm
            by(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards, iz-iz_min_r) = by_r(ix,    &
            iy, iz)*coeff_norm
            bz(ix-ix_min_r-nxguards, iy-iy_min_r-nyguards, iz-iz_min_r) = bz_r(ix,    &
            iy, iz)*coeff_norm
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE get_fields_mpi


  SUBROUTINE get_fields_mpi_lb
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    USE field_boundary
    IMPLICIT NONE
    REAL(num) :: coeff_norm, tmptime
    INTEGER(idp) :: ix, iy, iz
real(num)       :: pprime(4)
    ! Get global Fourier transform of all fields components 
    !CALL fft_forward_c2r_mpi
    CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jyf, jy_r)
    coeff_norm = 1.0_num/(nx_group*ny_group*nz_group)
coeff_norm=1.0_num
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
jy_r=0.0_num
call mpi_barrier(comm,errcode)
pprime(1)=5.0_num
pprime(2)=7.0_num
pprime(3)=13.0_num
pprime(4)=19.0_num
print*,"array_of_ranks_to_send_to",array_of_ranks_to_send_to
if(rank==1)print*,r_first_cell_to_recv
jy_r(103,1,:)=pprime(rank+1)!1._num
call mpi_barrier(comm,errcode)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=1,sizes_to_exchange_r_to_recv(z_coords+1)
      DO iy=iy_min_r, iy_max_r
        DO ix=ix_min_r, ix_max_r
!          ex(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1)) =&
!                 ex_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!          ey(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1)) =& 
!                ey_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!          ez(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1)) =&
!                 ez_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!          bx(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1)) =& 
!                bx_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!          by(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1)) =&
!                 by_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!          bz(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1)) =&
!                 bz_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!
!
!          jx(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1))=&
!                jx_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
          jy(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1))=jy_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!          jz(ix-ix_min_r-nxguards,iy-iy_min_r-nyguards,iz-1+r_first_cell_to_recv(z_coords+1))=&
!                 jz_r(ix,iy,iz-1+f_first_cell_to_send(z_coords+1))*coeff_norm
!

        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
  CALL generalized_comms_group_f2r(coeff_norm)
call mpi_barrier(comm,errcode)
  END SUBROUTINE get_fields_mpi_lb
  
  SUBROUTINE fft_forward_r2c_local(nfftx,nffty,nfftz) 
    USE fields
    USE fastfft
    USE fourier
    USE time_stat
    USE shared_data
    USE params
    REAL(num)   :: tmptime
    INTEGER(idp), INTENT(IN)    ::   nfftx,nffty,nfftz
    
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ex_r, exf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ey_r, eyf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ez_r, ezf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bx_r, bxf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, by_r, byf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bz_r, bzf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jx_r, jxf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jy_r, jyf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jz_r, jzf, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rhoold_r, rhooldf,plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rho_r, rhof, plan_r2c) 
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_forward_r2c_local

  SUBROUTINE fft_forward_r2c_mpi()
    USE fields
    USE fastfft
    USE mpi_fftw3
    USE time_stat
    USE shared_data
    USE params
    REAL(num)   :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ex_r, exf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ey_r, eyf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ez_r, ezf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bx_r, bxf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, by_r, byf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bz_r, bzf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r, jxf)
    CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r, jyf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r, jzf)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r, rhof)
  !  CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r, rhooldf)
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_forward_r2c_mpi

  SUBROUTINE fft_backward_c2r_local(nfftx,nffty,nfftz)
    USE fields
    USE fourier
    USE fastfft
    USE time_stat
    USE shared_data
    USE params
    REAL(num)   :: tmptime
    INTEGER(idp), INTENT(IN)     :: nfftx,nffty,nfftz

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, exf, ex_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, eyf, ey_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, ezf, ez_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, bxf, bx_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, byf, by_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, bzf, bz_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, jxf, jx_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, jyf, jy_r, plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, jzf, jz_r, plan_c2r)
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_backward_c2r_local

  SUBROUTINE fft_forward_c2r_mpi()
    USE fields
    USE fastfft
    USE mpi_fftw3
    USE time_stat
    USE shared_data
    USE params
    REAL(num)   :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, exf, ex_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, eyf, ey_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, ezf, ez_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, bxf, bx_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, byf, by_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, bzf, bz_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jxf, jx_r)
    CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jyf, jy_r)
  !  CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jzf, jz_r)

    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_forward_c2r_mpi

  SUBROUTINE push_psaotd_ebfielfs_2d() bind(C, name='push_psaotd_ebfields_2d')
    USE shared_data
    USE fields
    USE fourier
    USE time_stat
    USE params
    USE mpi_fftw3
    IMPLICIT NONE
    INTEGER(idp) ::  ix, iy, iz, nxx, nzz
    REAL(num) :: tmptime
    COMPLEX(cpx) :: bxfold, byfold, bzfold, exfold, eyfold, ezfold

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    nxx=size(exf(:, 1, 1))
    nzz=size(exf(1, 1, :))
    iy=1_idp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iz, exfold, eyfold, ezfold,     &
    !$OMP bxfold, byfold, bzfold) COLLAPSE(2)
    DO iz=1, nzz
        DO ix=1, nxx
          ! - Bx
          exfold=exf(ix, iy, iz)
          eyfold=eyf(ix, iy, iz)
          ezfold=ezf(ix, iy, iz)
          bxfold=bxf(ix, iy, iz)
          byfold=byf(ix, iy, iz)
          bzfold=bzf(ix, iy, iz)

          bxf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(4, 4)%block3dc(ix, iy,   &
          iz)*bxfold + cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc(ix, iy,        &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc(ix, iy,        &
          iz)*jyf(ix, iy, iz)  
          ! - By
          byf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc(ix, iy,   &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc(ix, iy,        &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(5, 3)%block3dc(ix, iy,        &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc(ix, iy,        &
          iz)*jxf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc(ix,   &
          iy, iz)*jzf(ix, iy, iz)


          ! - Bz
          bzf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc(ix, iy,   &
          iz)*bzfold + cc_mat(nmatrixes)%block_matrix2d(6, 2)%block3dc(ix, iy,        &
          iz)*eyfold +         cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc(ix,    &
          iy, iz)*jyf(ix, iy, iz)

          ! Push E a full time step
          ! - Ex
          exf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(1, 1)%block3dc(ix, iy,   &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(1, 7)%block3dc(ix, iy,        &
          iz)*jxf(ix, iy, iz)     + cc_mat(nmatrixes)%block_matrix2d(1,               &
          11)%block3dc(ix, iy, iz)*rhof(ix, iy, iz) +                                 &
          cc_mat(nmatrixes)%block_matrix2d(1, 10)%block3dc(ix, iy, iz)*rhooldf(ix,    &
          iy, iz)

          ! - Ey
          eyf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(2, 2)%block3dc(ix, iy,   &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc(ix, iy,        &
          iz)*bxfold  + cc_mat(nmatrixes)%block_matrix2d(2, 6)%block3dc(ix, iy,       &
          iz)*bzfold  + cc_mat(nmatrixes)%block_matrix2d(2, 8)%block3dc(ix, iy,       &
          iz)*jyf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(2, 11)%block3dc(ix,  &
          iy, iz)*rhof(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(2,              &
          10)%block3dc(ix, iy, iz)*rhooldf(ix, iy, iz)


          ! - Ez
          ezf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(3, 3)%block3dc(ix, iy,   &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(3, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(3, 9)%block3dc(ix, iy,        &
          iz)*jzf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(3, 11)%block3dc(ix,  &
          iy, iz)*rhof(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(3,              &
          10)%block3dc(ix, iy, iz)*rhooldf(ix, iy, iz)


        END DO
    END DO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE push_psaotd_ebfielfs_2d




  SUBROUTINE push_psaotd_ebfielfs() bind(C, name='push_psaotd_ebfields')
    USE shared_data
    USE fields
    USE fourier
    USE time_stat
    USE params
    USE mpi_fftw3
    IMPLICIT NONE
    INTEGER(idp) ::  ix, iy, iz, nxx, nyy, nzz
    REAL(num) :: tmptime
    COMPLEX(cpx) :: bxfold, byfold, bzfold, exfold, eyfold, ezfold

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    nxx=size(exf(:, 1, 1))
    nyy=size(exf(1, :, 1))
    nzz=size(exf(1, 1, :))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, exfold, eyfold, ezfold,     &
    !$OMP bxfold, byfold, bzfold) COLLAPSE(3)
    DO iz=1, nzz
      DO iy=1, nyy
        DO ix=1, nxx
          ! - Bx
          exfold=exf(ix, iy, iz)
          eyfold=eyf(ix, iy, iz)
          ezfold=ezf(ix, iy, iz)
          bxfold=bxf(ix, iy, iz)
          byfold=byf(ix, iy, iz)
          bzfold=bzf(ix, iy, iz)

          bxf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(4, 4)%block3dc(ix, iy,   &
          iz)*bxfold + cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc(ix, iy,        &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc(ix, iy,        &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc(ix, iy,        &
          iz)*jyf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc(ix,   &
          iy, iz)*jzf(ix, iy, iz)

          ! - By
          byf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc(ix, iy,   &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc(ix, iy,        &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(5, 3)%block3dc(ix, iy,        &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc(ix, iy,        &
          iz)*jxf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc(ix,   &
          iy, iz)*jzf(ix, iy, iz)


          ! - Bz
          bzf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc(ix, iy,   &
          iz)*bzfold + cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc(ix, iy,        &
          iz)*exfold+ cc_mat(nmatrixes)%block_matrix2d(6, 2)%block3dc(ix, iy,         &
          iz)*eyfold+ cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc(ix, iy,         &
          iz)*jxf(ix, iy, iz)+ cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc(ix,    &
          iy, iz)*jyf(ix, iy, iz)

          ! Push E a full time step
          ! - Ex
          exf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(1, 1)%block3dc(ix, iy,   &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc(ix, iy,        &
          iz)*bzfold + cc_mat(nmatrixes)%block_matrix2d(1, 7)%block3dc(ix, iy,        &
          iz)*jxf(ix, iy, iz)     + cc_mat(nmatrixes)%block_matrix2d(1,               &
          11)%block3dc(ix, iy, iz)*rhof(ix, iy, iz) +                                 &
          cc_mat(nmatrixes)%block_matrix2d(1, 10)%block3dc(ix, iy, iz)*rhooldf(ix,    &
          iy, iz)

          ! - Ey
          eyf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(2, 2)%block3dc(ix, iy,   &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc(ix, iy,        &
          iz)*bxfold  + cc_mat(nmatrixes)%block_matrix2d(2, 6)%block3dc(ix, iy,       &
          iz)*bzfold  + cc_mat(nmatrixes)%block_matrix2d(2, 8)%block3dc(ix, iy,       &
          iz)*jyf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(2, 11)%block3dc(ix,  &
          iy, iz)*rhof(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(2,              &
          10)%block3dc(ix, iy, iz)*rhooldf(ix, iy, iz)


          ! - Ez
          ezf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(3, 3)%block3dc(ix, iy,   &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc(ix, iy,        &
          iz)*bxfold + cc_mat(nmatrixes)%block_matrix2d(3, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(3, 9)%block3dc(ix, iy,        &
          iz)*jzf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(3, 11)%block3dc(ix,  &
          iy, iz)*rhof(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(3,              &
          10)%block3dc(ix, iy, iz)*rhooldf(ix, iy, iz)


        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE push_psaotd_ebfielfs



  SUBROUTINE init_plans_blocks() bind(C, name='init_plans_blocks_pxr')
    USE shared_data
    USE fastfft
    Use fourier
    USE fftw3_fortran
    USE fields
    USE omp_lib
    USE params

    INTEGER(idp) :: nfftx, nffty, nfftz, nopenmp
#ifdef _OPENMP
    nopenmp=OMP_GET_MAX_THREADS()
#else
    nopenmp=1
#endif

    IF (fftw_with_mpi) THEN
      nfftx=nx_global
      nffty=ny_global
      nfftz=nz_global
    ELSE
#if defined(LIBRARY)
      nfftx=nx+2*nxguards+1
      nffty=ny+2*nyguards+1
      nfftz=nz+2*nzguards+1
#else
      nfftx=nx+2*nxguards
      nffty=ny+2*nyguards
      nfftz=nz+2*nzguards
#endif
    ENDIF
    CALL init_gpstd()
    IF(rank==0) WRITE(0, *) 'INIT GPSTD MATRIX DONE'
    IF (fftw_with_mpi) THEN
      CALL init_plans_fourier_mpi(nopenmp)
    ELSE
      IF(c_dim ==3) THEN
        CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp, nfftx, nffty, nfftz,ex_r, exf,  &
        plan_r2c, INT(FFTW_MEASURE, idp), INT(FFTW_FORWARD, idp))
        CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp, nfftx, nffty, nfftz, exf,ex_r,  &
        plan_c2r, INT(FFTW_MEASURE, idp), INT(FFTW_BACKWARD, idp))
      ELSE IF(c_dim == 2) THEN
        CALL fast_fftw_create_plan_r2c_2d_dft(nopenmp, nfftx, nfftz, ex_r, exf,&
        plan_r2c, INT(FFTW_MEASURE, idp), INT(FFTW_FORWARD, idp))
        CALL fast_fftw_create_plan_c2r_2d_dft(nopenmp, nfftx, nfftz, exf, ex_r,&
        plan_c2r, INT(FFTW_MEASURE, idp), INT(FFTW_BACKWARD, idp))
      ENDIF
    ENDIF
    IF(rank==0) WRITE(0, *) 'INIT GPSTD PLANS DONE'
  END SUBROUTINE init_plans_blocks
END MODULE fourier_psaotd

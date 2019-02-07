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
  USE fields, ONLY: bx, by, bz, ex, ey, ez, l_nodalgrid, norderx, nordery, norderz,  &
    nxguards, nxs, nyguards, nys, nzguards, nzs, xcoeffs, ycoeffs, zcoeffs
  USE mpi
  USE params, ONLY: dt, it
  USE picsar_precision, ONLY: num
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE

  REAL(num) :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2).AND.(norderz.eq.2)) then
    CALL pxrpush_em3d_bvec( &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         ex, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ey, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ez, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bx, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         by, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bz, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         0.5_num*dt/dx, 0.5_num*dt/dy, 0.5_num*dt/dz)

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
  USE constants, ONLY: eps0, mu0
  USE fields, ONLY: bx, by, bz, electro_energy_mpi, electro_energy_total,            &
    electromagn_energy_mpi, electromagn_energy_total, ex, ey, ez,                    &
    magnetic_energy_mpi, magneto_energy_total
  USE mpi
  USE picsar_precision, ONLY: isp, num
  USE shared_data, ONLY: comm, dx, dy, dz, errcode, nx, ny, nz
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
  USE fields, ONLY: bxy, bxz, byx, byz, bzx, bzy, exy, exz, eyx, eyz, ezx, ezy,      &
    nxguards, nyguards, nzguards, sigma_x_b, sigma_x_e, sigma_y_b, sigma_y_e,        &
    sigma_z_b, sigma_z_e
  USE mpi
  USE omp_lib
  USE params, ONLY: it
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: c_dim, nx, ny, nz
  USE time_stat, ONLY: localtimes, timestat_itstart

  IMPLICIT NONE
  INTEGER(idp)  :: ix,iy,iz
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF
  IF(c_dim == 3) THEN 
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz = -nzguards,nz+nzguards-1
      DO iy = -nyguards,ny+nyguards-1
        DO ix = -nxguards,nx+nxguards-1
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
  ELSE IF(c_dim==2) THEN
    iy=0_idp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,  iz) COLLAPSE(2)
    DO iz = -nzguards,nz+nzguards-1
        DO ix = -nxguards,nx+nxguards-1
          exz(ix,iy,iz) = sigma_z_e(iz) *exz(ix,iy,iz)
          eyx(ix,iy,iz) = sigma_x_e(ix) *eyx(ix,iy,iz)
          eyz(ix,iy,iz) = sigma_z_e(iz) *eyz(ix,iy,iz)
          ezx(ix,iy,iz) = sigma_x_e(ix) *ezx(ix,iy,iz)
          bxz(ix,iy,iz) = sigma_z_b(iz) *bxz(ix,iy,iz)
          byx(ix,iy,iz) = sigma_x_b(ix) *byx(ix,iy,iz)
          byz(ix,iy,iz) = sigma_z_b(iz) *byz(ix,iy,iz)
          bzx(ix,iy,iz) = sigma_x_b(ix) *bzx(ix,iy,iz)
        ENDDO
      ENDDO
    !$OMP END PARALLEL DO
  ENDIF

   
  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE field_damping_bcs


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
  USE fields, ONLY: bx, bxy, bxz, by, byx, byz, bz, bzx, bzy, ex, exy, exz, ey, eyx, &
    eyz, ez, ezx, ezy
  USE mpi
  USE omp_lib
  USE params, ONLY: it
  USE picsar_precision, ONLY: idp, num
  USE time_stat, ONLY: localtimes, timestat_itstart

  IMPLICIT NONE 
  INTEGER(idp)  :: ix,iy,iz,ixx,iyy,izz,&
  ubound_s(3), ubound_f(3), lbound_s(3), lbound_f(3)
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ubound_s = UBOUND(exy)
  lbound_s = LBOUND(exy)
  ubound_f = UBOUND(ex)
  lbound_f = LBOUND(ex)

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz,ixx ,iyy, izz) COLLAPSE(3)
  DO iz = lbound_f(3),ubound_f(3)
    DO iy = lbound_f(2),ubound_f(2)
      DO ix = lbound_f(1),ubound_f(1)
        ixx = ix - lbound_f(1) + lbound_s(1)
        iyy = iy - lbound_f(2) + lbound_s(2)
        izz = iz - lbound_f(3) + lbound_s(3)

        ex(ix,iy,iz) = exy(ixx,iyy,izz) + exz(ixx,iyy,izz)
        ey(ix,iy,iz) = eyx(ixx,iyy,izz) + eyz(ixx,iyy,izz)
        ez(ix,iy,iz) = ezx(ixx,iyy,izz) + ezy(ixx,iyy,izz)
        bx(ix,iy,iz) = bxy(ixx,iyy,izz) + bxz(ixx,iyy,izz)
        by(ix,iy,iz) = byx(ixx,iyy,izz) + byz(ixx,iyy,izz)
        bz(ix,iy,iz) = bzx(ixx,iyy,izz) + bzy(ixx,iyy,izz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE merge_fields

! ________________________________________________________________________________________
!> @brief
!> Mege Electric splitted fields when using absorbing_bcs to compute real electric field
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________

SUBROUTINE merge_e_fields()
  USE fields, ONLY: ex, exy, exz, ey, eyx, eyz, ez, ezx, ezy
  USE mpi
  USE omp_lib
  USE params, ONLY: it
  USE picsar_precision, ONLY: idp, num
  USE time_stat, ONLY: localtimes, timestat_itstart

  IMPLICIT NONE
  INTEGER(idp)  :: ix,iy,iz,ixx,iyy,izz,&
  ubound_s(3), ubound_f(3), lbound_s(3), lbound_f(3)
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ubound_s = UBOUND(exy)
  lbound_s = LBOUND(exy)
  ubound_f = UBOUND(ex)
  lbound_f = LBOUND(ex)

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz,ixx ,iyy, izz) COLLAPSE(3)
  DO iz = lbound_f(3),ubound_f(3)
    DO iy = lbound_f(2),ubound_f(2)
      DO ix = lbound_f(1),ubound_f(1)
        ixx = ix - lbound_f(1) + lbound_s(1)
        iyy = iy - lbound_f(2) + lbound_s(2)
        izz = iz - lbound_f(3) + lbound_s(3)

        ex(ix,iy,iz) = exy(ixx,iyy,izz) + exz(ixx,iyy,izz)
        ey(ix,iy,iz) = eyx(ixx,iyy,izz) + eyz(ixx,iyy,izz)
        ez(ix,iy,iz) = ezx(ixx,iyy,izz) + ezy(ixx,iyy,izz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE merge_e_fields

! ________________________________________________________________________________________
!> @brief
!> Mege Electric splitted fields when using absorbing_bcs to compute real magnetic field
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________

SUBROUTINE merge_b_fields()
  USE fields, ONLY: bx, bxy, bxz, by, byx, byz, bz, bzx, bzy
  USE mpi
  USE omp_lib
  USE params, ONLY: it
  USE picsar_precision, ONLY: idp, num
  USE time_stat, ONLY: localtimes, timestat_itstart

  IMPLICIT NONE
  INTEGER(idp)  :: ix,iy,iz,ixx,iyy,izz,&
  ubound_s(3), ubound_f(3), lbound_s(3), lbound_f(3)
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ubound_s = UBOUND(bxy)
  lbound_s = LBOUND(bxy)
  ubound_f = UBOUND(bx)
  lbound_f = LBOUND(bx)

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz,ixx ,iyy, izz) COLLAPSE(3)
  DO iz = lbound_f(3),ubound_f(3)
    DO iy = lbound_f(2),ubound_f(2)
      DO ix = lbound_f(1),ubound_f(1)
        ixx = ix - lbound_f(1) + lbound_s(1)
        iyy = iy - lbound_f(2) + lbound_s(2)
        izz = iz - lbound_f(3) + lbound_s(3)

        bx(ix,iy,iz) = bxy(ixx,iyy,izz) + bxz(ixx,iyy,izz)
        by(ix,iy,iz) = byx(ixx,iyy,izz) + byz(ixx,iyy,izz)
        bz(ix,iy,iz) = bzx(ixx,iyy,izz) + bzy(ixx,iyy,izz)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE merge_b_fields


   
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
  USE constants, ONLY: clight, mu0
  USE fields, ONLY: bx, by, bz, ex, ey, ez, jx, jy, jz, l_nodalgrid, norderx,        &
    nordery, norderz, nxguards, nxs, nyguards, nys, nzguards, nzs, xcoeffs, ycoeffs, &
    zcoeffs
  USE mpi
  USE params, ONLY: dt, it
  USE picsar_precision, ONLY: num
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE

  REAL(num) :: tmptime
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2).AND.(norderz.eq.2)) then

    CALL pxrpush_em3d_evec( &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         ex, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ey, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ez, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bx, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         by, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bz, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         jx, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         jy, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         jz, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         clight**2*mu0*dt, clight**2*dt/dx , clight**2*dt/dy, clight**2*dt/dz)

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
  USE fields, ONLY: bx, by, bz, ex, ey, ez, l_nodalgrid, norderx, nordery, norderz,  &
    nxguards, nxs, nyguards, nys, nzguards, nzs, xcoeffs, ycoeffs, zcoeffs
  USE mpi
  USE params, ONLY: dt, it
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE

  REAL(num) :: tmptime
  INTEGER(idp) :: iy = 0

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(norderz.eq.2)) then

    CALL pxrpush_em2d_bvec( (/-nxs, -nzs/), (/nx+nxs, nz+nzs/), (/-nxs, -nzs/), &
         (/nx+nxs, nz+nzs/), (/-nxs, -nzs/), (/nx+nxs, nz+nzs/), &
         ex(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         ey(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         ez(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         bx(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         by(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         bz(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         0.5_num*dt/dx ,0._num, 0.5_num*dt/dz)

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
  USE constants, ONLY: clight, mu0
  USE fields, ONLY: bx, by, bz, ex, ey, ez, jx, jy, jz, l_nodalgrid, norderx,        &
    nordery, norderz, nxguards, nxs, nyguards, nzguards, nzs, xcoeffs, ycoeffs,      &
    zcoeffs
  USE mpi
  USE params, ONLY: dt, it
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE

  REAL(num) :: tmptime,mdt
  INTEGER(idp) :: iy = 0

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF
  mdt = mu0*clight**2*dt
  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(norderz.eq.2)) then

    CALL pxrpush_em2d_evec( (/-nxs, -nzs/), (/nx+nxs, nz+nzs/), (/-nxs, -nzs/), &
         (/nx+nxs, nz+nzs/), (/-nxs, -nzs/), (/nx+nxs, nz+nzs/), &
         ex(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         ey(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         ez(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         bx(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         by(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         bz(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         jx(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         jy(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         jz(:,iy,:), (/-nxguards, -nzguards/), (/nx+nxguards, nz+nzguards/), &
         mdt, clight**2*dt/dx ,0., clight**2*dt/dz)

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

  SUBROUTINE push_psatd_ebfield
  USE fields, ONLY: g_spectral
#if defined(FFTW)
  USE fourier_psaotd
  USE matrix_data, ONLY: nmatrixes
#endif
  USE mpi
  USE params, ONLY: it
  USE picsar_precision, ONLY: num
  USE shared_data, ONLY: c_dim, fftw_hybrid, fftw_with_mpi
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE

  REAL(num) :: tmptime, tmptime_m

#if defined(DEBUG)
  WRITE(0, *) "push psatd ebfield : start"
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
    IF(c_dim == 3) THEN
      CALL push_psaotd_ebfielfs_3d! - PUSH PSATD
    ELSE IF(c_dim == 2) THEN
      CALL push_psaotd_ebfielfs_2d
    ENDIF
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
  WRITE(0, *) "push psatd ebfield : end"
#endif

  END SUBROUTINE push_psatd_ebfield

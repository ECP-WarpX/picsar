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
  USE params
  USE constants

  IMPLICIT NONE
  INTEGER(idp)  :: ix,iy,iz
  REAL(num)     :: tmptime, ieps0
  INTEGER(idp) :: lbound_e(3), ubound_e(3)
  
  ieps0 = 1._num/eps0
  IF(.NOT. u_pml) THEN 
    lbound_e = LBOUND(exy)
    ubound_e = UBOUND(exy) 

  ELSE IF(u_pml) THEN
    lbound_e = LBOUND(dex)
    ubound_e = UBOUND(dex)
  ENDIF   
  ubound_e = ubound_e - 1_idp 
  IF(c_dim == 2) THEN
    ubound_e(2) = lbound_e(2)
  ENDIF
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF


  IF(.NOT. u_pml) THEN
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz = lbound_e(3),ubound_e(3)
      DO iy = lbound_e(2),ubound_e(2)
        DO ix = lbound_e(1),ubound_e(1)
          exy(ix,iy,iz) = a_y_e(iy) *exy(ix,iy,iz)
          exz(ix,iy,iz) = a_z_e(iz) *exz(ix,iy,iz)
          eyx(ix,iy,iz) = a_x_e(ix)*eyx(ix,iy,iz)
          eyz(ix,iy,iz) = a_z_e(iz) *eyz(ix,iy,iz)
          ezx(ix,iy,iz) = a_x_e(ix) *ezx(ix,iy,iz)
          ezy(ix,iy,iz) = a_y_e(iy) *ezy(ix,iy,iz)
          bxy(ix,iy,iz) = a_y_b(iy) *bxy(ix,iy,iz)
          bxz(ix,iy,iz) = a_z_b(iz) *bxz(ix,iy,iz)
          byx(ix,iy,iz) = a_x_b(ix) *byx(ix,iy,iz)
          byz(ix,iy,iz) = a_z_b(iz) *byz(ix,iy,iz)
          bzx(ix,iy,iz) = a_x_b(ix) *bzx(ix,iy,iz)
          bzy(ix,iy,iz) = a_y_b(iy) *bzy(ix,iy,iz)

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ELSE IF(u_pml) THEN
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz = lbound_e(3),ubound_e(3)
      DO iy = lbound_e(2),ubound_e(2)
        DO ix = lbound_e(1),ubound_e(1)
         ! - DAMP bx ,by ,bz, dex, dey, dez along y z x respectively

          dex(ix,iy,iz) = a_y_e(iy) *dex(ix,iy,iz)
          dey(ix,iy,iz) = a_z_e(iz) *dey(ix,iy,iz)
          dez(ix,iy,iz) = a_x_e(ix) *dez(ix,iy,iz)
          bx(ix,iy,iz) = a_y_b(iy) *bx(ix,iy,iz)
          by(ix,iy,iz) = a_z_b(iz) *by(ix,iy,iz)
          bz(ix,iy,iz) = a_x_b(ix) *bz(ix,iy,iz)
          
          ! Compute E
          ex(ix,iy,iz) = a_z_e(iz)*ex(ix,iy,iz) + b_z_e(iz)/b_x_b(ix)*(dex(ix,iy,iz) - a_x_b(ix)*dexold(ix,iy,iz))*ieps0
          ey(ix,iy,iz) = a_x_e(ix)*ey(ix,iy,iz) + b_x_e(ix)/b_y_b(iy)*(dey(ix,iy,iz) - a_y_b(iy)*deyold(ix,iy,iz))*ieps0
          ez(ix,iy,iz) = a_y_e(iy)*ez(ix,iy,iz) + b_y_e(iy)/b_z_b(iz)*(dez(ix,iy,iz) - a_z_b(iz)*dezold(ix,iy,iz))*ieps0
          ! Compute H
          hx(ix,iy,iz) = a_z_b(iz)*hx(ix,iy,iz) + b_z_b(iz)/b_x_e(ix)*(bx(ix,iy,iz) - a_x_e(ix)*bxold(ix,iy,iz))*imu0
          hy(ix,iy,iz) = a_x_b(ix)*hy(ix,iy,iz) + b_x_b(ix)/b_y_e(iy)*(by(ix,iy,iz) - a_y_e(iy)*byold(ix,iy,iz))*imu0
          hz(ix,iy,iz) = a_y_b(iy)*hz(ix,iy,iz) + b_y_b(iy)/b_z_e(iz)*(bz(ix,iy,iz) - a_z_e(iz)*bzold(ix,iy,iz))*imu0

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  
  ENDIF 
  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END subroutine field_damping_bcs


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

SUBROUTINE damp_e_field
  USE fields
  USE shared_data
  USE constants
  USE omp_lib
  USE time_stat
  USE params

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
          exy(ix,iy,iz) = EXP(-dt*sigma_y_e(iy)) *exy(ix,iy,iz)
          exz(ix,iy,iz) = EXP(-dt*sigma_z_e(iz)) *exz(ix,iy,iz)
          eyx(ix,iy,iz) = EXP(-dt*sigma_x_e(ix)) *eyx(ix,iy,iz)
          eyz(ix,iy,iz) = EXP(-dt*sigma_z_e(iz)) *eyz(ix,iy,iz)
          ezx(ix,iy,iz) = EXP(-dt*sigma_x_e(ix)) *ezx(ix,iy,iz)
          ezy(ix,iy,iz) = EXP(-dt*sigma_y_e(iy)) *ezy(ix,iy,iz)

        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ELSE IF(c_dim==2) THEN
    iy=0_idp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,  iz) COLLAPSE(2)
    DO iz = -nzguards,nz+nzguards-1
        DO ix = -nzguards,nx+nxguards-1
          ex(ix,iy,iz)  = EXP(-dt*sigma_z_e(iz)) *ex(ix,iy,iz)
          eyx(ix,iy,iz) = EXP(-dt*sigma_x_e(ix)) *eyx(ix,iy,iz)
          eyz(ix,iy,iz) = EXP(-dt*sigma_z_e(iz)) *eyz(ix,iy,iz)
          ez(ix,iy,iz)  = EXP(-dt*sigma_x_e(ix)) *ez(ix,iy,iz)

        ENDDO
      ENDDO
    !$OMP END PARALLEL DO
  ENDIF

   
  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END subroutine damp_e_field

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

SUBROUTINE damp_b_field
  USE fields
  USE shared_data
  USE constants
  USE omp_lib
  USE time_stat
  USE params

  IMPLICIT NONE
  INTEGER(idp)  :: ix,iy,iz
  REAL(num)     :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF
  IF(c_dim == 3) THEN 
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz = -nzguards,nz-nzguards-1
      DO iy = -nyguards,ny+nyguards-1
        DO ix = -nxguards,nx+nxguards-1
          bxy(ix,iy,iz) = EXP(-dt*sigma_y_b(iy)) *bxy(ix,iy,iz)
          bxz(ix,iy,iz) = EXP(-dt*sigma_z_b(iz)) *bxz(ix,iy,iz)
          byx(ix,iy,iz) = EXP(-dt*sigma_x_b(ix)) *byx(ix,iy,iz)
          byz(ix,iy,iz) = EXP(-dt*sigma_z_b(iz)) *byz(ix,iy,iz)
          bzx(ix,iy,iz) = EXP(-dt*sigma_x_b(ix)) *bzx(ix,iy,iz)
          bzy(ix,iy,iz) = EXP(-dt*sigma_y_b(iy)) *bzy(ix,iy,iz)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ELSE IF(c_dim==2) THEN
    iy=0_idp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,  iz) COLLAPSE(2)
    DO iz = -nzguards,nz+nzguards-1
        DO ix = -nxguards,nx+nxguards-1
          bx(ix,iy,iz)  = EXP(-dt*sigma_z_b(iz)) *bx(ix,iy,iz)
          byx(ix,iy,iz) = EXP(-dt*sigma_x_b(ix)) *byx(ix,iy,iz)
          byz(ix,iy,iz) = EXP(-dt*sigma_z_b(iz)) *byz(ix,iy,iz)
          bz(ix,iy,iz)  = EXP(-dt*sigma_x_b(ix)) *bz(ix,iy,iz)
        ENDDO
      ENDDO
    !$OMP END PARALLEL DO
  ENDIF

   
  IF (it.ge.timestat_itstart) THEN
    localtimes(26) = localtimes(26) + (MPI_WTIME() - tmptime)
  ENDIF

END subroutine damp_b_field

SUBROUTINE set_b_to_0

USE fields
use shared_data

bx(-nxguards:-1,:,:) = 0.0_num
bx(nx+1:nx+nxguards,:,:) = 0.0_num

by(:,-nyguards:-1,:) = 0.0_num
by(:,ny+1:ny+nyguards,:) = 0.0_num

bz(:,:,-nzguards:-1) = 0.0_num
bz(:,:,nz+1:nz+nzguards) = 0.0_num
if( u_pml) THEN
hx(-nxguards:-1,:,:) = 0.0_num
hx(nx+1:nx+nxguards,:,:) = 0.0_num
hy(:,-nyguards:-1,:) = 0.0_num
hy(:,ny+1:ny+nyguards,:) = 0.0_num
hz(:,:,-nzguards:-1) = 0.0_num
hz(:,:,nz+1:nz+nzguards) = 0.0_num
else

byx(:,-nyguards:-1,:) = 0.0_num
byx(:,ny+1:ny+nyguards,:) = 0.0_num
byz(:,-nyguards:-1,:) = 0.0_num
byz(:,ny+1:ny+nyguards,:) = 0.0_num


endif
END SUBROUTINE set_b_to_0

SUBROUTINE set_e_to_0

USE fields
USE shared_data

ey(-nxguards:-1,:,:) = 0.0_num
ey(nx+1:nx+nxguards,:,:) = 0.0_num

ez(-nxguards:-1,:,:) = 0.0_num
ez(nx+1:nx+nxguards,:,:) = 0.0_num

ez(:,-nyguards:-1,:) = 0.0_num
ez(:,ny+1:ny+nyguards,:) = 0.0_num

ex(:,-nyguards:-1,:) = 0.0_num
ex(:,ny+1:ny+nyguards,:) = 0.0_num

ex(:,:,-nzguards:-1) = 0.0_num
ex(:,:,nz+1:nz+nzguards) = 0.0_num

ey(:,:,-nzguards:-1) = 0.0_num
ey(:,:,nz+1:nz+nzguards) = 0.0_num

if(u_pml) then
dey(-nxguards:-1,:,:) = 0.0_num
dey(nx+1:nx+nxguards,:,:) = 0.0_num
dez(-nxguards:-1,:,:) = 0.0_num
dez(nx+1:nx+nxguards,:,:) = 0.0_num
dez(:,-nyguards:-1,:) = 0.0_num
dez(:,ny+1:ny+nyguards,:) = 0.0_num
dex(:,-nyguards:-1,:) = 0.0_num
dex(:,ny+1:ny+nyguards,:) = 0.0_num
dex(:,:,-nzguards:-1) = 0.0_num
dex(:,:,nz+1:nz+nzguards) = 0.0_num
dey(:,:,-nzguards:-1) = 0.0_num
dey(:,:,nz+1:nz+nzguards) = 0.0_num
else
eyz(-nxguards:-1,:,:) = 0.0_num
eyz(nx+1:nx+nxguards,:,:) = 0.0_num

eyx(-nxguards:-1,:,:) = 0.0_num
eyx(nx+1:nx+nxguards,:,:) = 0.0_num

eyz(:,:,-nzguards:-1) = 0.0_num
eyz(:,:,nz+1:nz+nzguards) = 0.0_num

eyx(:,:,-nzguards:-1) = 0.0_num
eyx(:,:,nz+1:nz+nzguards) = 0.0_num

endif


END SUBROUTINE set_e_to_0



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
  USE params

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
  IF(absorbing_bcs) THEN
    IF(.NOT. u_pml) THEN
        CALL pxrpush_em_standard_pml_2d_bvec(ex,ey,ez,bx,by,bz,exy,exz,eyx,eyz,bzx,bzy,bxy,bxz,byx,byz,bzx,bzy,           &
        1._num/dx, 0._num, 1._num/dz, nx, ny, nz, nxguards, 0_idp, nzguards, nxs, & 
        0_idp, nzs,   l_nodalgrid)
    ELSE
      CALL pxrpush_em_upml_2d_bvec(hx,hy,hz,ex,ey,ez,bx,by,bz,dt,   &
      dt/dx, 0._num, dt/dz, nx, ny, nz, nxguards, 0_idp, nzguards, nxs, &
      0_idp, nzs,   l_nodalgrid) 
    ENDIF                       
  ELSE
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


  IF(absorbing_bcs) THEN
    IF(.NOT. u_pml) THEN
        CALL pxrpush_em_standard_pml_2d_evec(ex,ey,ez,bx,by,bz,exy,exz,eyx, eyz,ezx, ezy, bxy, bxz,byx, byz, bzx, &
        bzy ,jx,jy,jz,mdt, clight**2*1._num/dx,clight**2*1._num/dy, &
        clight**2*1._num/dz, nx,ny,nz, nxguards, nyguards, nzguards,nxs,0_idp,nzs,&
        l_nodalgrid)
    ELSE
      CALL pxrpush_em_upml_2d_evec(dex,dey,dez,hx,hy,hz,&
      ex,ey,ez,bx,by,bz,jx,jy,jz,dt,mdt, dt/dx,dt/dy, &
      dt/dz, nx,ny,nz, nxguards, nyguards, nzguards,nxs,0_idp,nzs,&
      l_nodalgrid)
    ENDIF
  ELSE
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
SUBROUTINE push_psatd_ebfield() 
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
integer :: i
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
  WRITE(0, *) "push psatd ebfield: end"
#endif

END SUBROUTINE push_psatd_ebfield



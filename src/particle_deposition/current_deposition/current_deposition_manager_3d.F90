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
! CURRENT_DEPOSITION_MANAGER_3D.F90
!
! Developers
! Henri Vincenti, ! Mathieu Lobet
!
! Description:
! This file contains subroutines for managing the current deposition in 3D.
!
! List of subroutines:
! - pxrdepose_currents_on_grid_jxjyjz
!
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq
! - pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp
! - pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Generic subroutine for current deposition on one tile
!>
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox, noy,    &
  noz, current_depo_algo)
  USE picsar_precision, ONLY: idp, num
  USE fields, ONLY: l_nodalgrid
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard,           &
  current_depo_algo
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(3), nvalid(3)
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)

  call depose_jxjyjz_generic( jx, nguard, nvalid, jy, nguard, nvalid, jz, nguard,     &
  nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy,  &
  dz, nox, noy, noz, l_nodalgrid, current_depo_algo)
END SUBROUTINE depose_jxjyjz

! ________________________________________________________________________________________
!> @brief
!> Esirkepov subroutine for current deposition on one tile
!>
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, &
  w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,   &
  noy, noz)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin

  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(3), nvalid(3)
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)

  CALL pxr_depose_jxjyjz_esirkepov_n( jx, nguard, nvalid, jy, nguard, nvalid, jz,     &
  nguard, nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt,  &
  dx, dy, dz, nox, noy, noz, .TRUE._idp, .FALSE._idp)
END SUBROUTINE depose_jxjyjz_esirkepov

! ________________________________________________________________________________________
!> @brief
!> Generic subroutines for current deposition, adapted for field
!> arrays having different sizes depending on their nodal/cell-centered nature
!>
!> @details
!> This routine calls the relevant current deposition routine depending
!> on the order of the particle shape and the selected algorithm.
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_generic( jx, jx_nguard, jx_nvalid, jy, jy_nguard, jy_nvalid, &
  jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin,    &
  zmin, dt, dx, dy, dz, nox, noy, noz, l_nodal, current_depo_algo)     !#do not wrap
  USE picsar_precision, ONLY: idp, num, lp
  IMPLICIT NONE
  INTEGER(idp) :: np, nox, noy, noz, current_depo_algo
  INTEGER(idp), intent(in)                :: jx_nguard(3), jx_nvalid(3),              &
  jy_nguard(3), jy_nvalid(3), jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1,                                          &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1,                                          &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,                                          &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  LOGICAL(lp), intent(in)  :: l_nodal
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin

  ! Maintain variables nx, ny, nz, nxguard, nyguard, nzguard for compilation
  ! and for compatibility with automated tests, although they will not be used
  ! in the future
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard
  nx = jx_nvalid(1)-1
  ny = jx_nvalid(2)-1
  nz = jx_nvalid(3)-1
  nxguard = jx_nguard(1)
  nyguard = jx_nguard(2)
  nzguard = jx_nguard(3)

  Select CASE(current_depo_algo)
    ! Scalar classical current deposition subroutines
  CASE(3)

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      CALL depose_jxjyjz_scalar_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)
    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      CALL depose_jxjyjz_scalar_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)
    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      CALL depose_jxjyjz_scalar_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)
    ELSE
      CALL pxr_depose_jxjyjz_esirkepov_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, nox, noy, noz, .TRUE._idp, .FALSE._idp)
    ENDIF
    ! Optimized classical current deposition
  CASE(2)

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      CALL depose_jxjyjz_vecHVv2_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)
    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      CALL depose_jxjyjz_vecHVv2_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)
    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      CALL depose_jxjyjz_vecHVv3_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)
    ELSE
      CALL pxr_depose_jxjyjz_esirkepov_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, nox, noy, noz, .TRUE._idp, .FALSE._idp)
    ENDIF
    ! Esirkepov non optimized
  CASE(1)

    CALL pxr_depose_jxjyjz_esirkepov_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
    jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
    xmin, ymin, zmin, dt, dx, dy, dz, nox, noy, noz, .TRUE._idp, .FALSE._idp)
    ! Optimized Esirkepov
  CASE DEFAULT

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      CALL depose_jxjyjz_esirkepov_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz)
    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      CALL depose_jxjyjz_esirkepov_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz)
    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      CALL depose_jxjyjz_esirkepov_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz)
    ELSE
      CALL pxr_depose_jxjyjz_esirkepov_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz, nox, noy, noz, .TRUE._idp, .FALSE._idp)
    ENDIF

  END SELECT

END SUBROUTINE depose_jxjyjz_generic

! ________________________________________________________________________________________
!> @brief
!> Main subroutine for managing the current deposition across tiles
!
!> @details
!> This subroutine is called in submain.F90 in step().
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2015-2016
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz
  USE fields, ONLY: jx, jy, jz, nox, noy, noz, nxjguards, nyjguards, nzjguards
  USE mpi
  USE params, ONLY: currdepo, dt, it, lvec_curr_depo
  USE particle_properties, ONLY: nspecies
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz, xmin, ymin, zmin
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE
  REAL(num) :: tdeb, tend

  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE
    ! ____________________________________________________________________________________
    ! Generic current deposition routine
    SUBROUTINE depose_jxjyjz(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
      xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,     &
      noy, noz, current_depo_algo)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard,       &
      current_depo_algo
      REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                  &
      -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
    END SUBROUTINE

    ! Interface for subroutine with no reduction - classical deposition order 1
    SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1(jxcells, jycells, jzcells, np, ncells,   &
      xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx,    &
      ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
      REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
      REAL(num), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      REAL(num), INTENT(IN) :: q, dt, dx, dy, dz, xmin, ymin, zmin
      LOGICAL(lp), INTENT(IN) :: l_nodal
      INTEGER(idp) :: ncx, ncy, ncz
      INTEGER(idp) :: lvect
    END SUBROUTINE

    ! Interface for subroutine with no reduction - classical deposition order 2
    SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2(jxcells, jycells, jzcells, np, ncells,   &
      xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx,    &
      ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
      REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
      REAL(num), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      REAL(num), INTENT(IN) :: q, dt, dx, dy, dz, xmin, ymin, zmin
      LOGICAL(lp), INTENT(IN) :: l_nodal
      INTEGER(idp) :: ncx, ncy, ncz
      INTEGER(idp) :: lvect
    END SUBROUTINE

    ! Interface for subroutine with no reduction - classical deposition order 3
    SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3(jxcells, jycells, jzcells, np, ncells,   &
      xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx,    &
      ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
      REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
      REAL(num), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      REAL(num), INTENT(IN) :: q, dt, dx, dy, dz, xmin, ymin, zmin
      LOGICAL(lp), INTENT(IN) :: l_nodal
      INTEGER(idp) :: ncx, ncy, ncz
      INTEGER(idp) :: lvect
    END SUBROUTINE

    SUBROUTINE current_reduction_1_1_1(jx, jy, jz, jxcells, jycells, jzcells, ncells, &
      nx, ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                 :: nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
      INTEGER(idp), INTENT(IN)                 :: nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN), DIMENSION(8, ncells):: jxcells, jycells, jzcells
    END SUBROUTINE

    SUBROUTINE current_reduction_2_2_2(jx, jy, jz, jxcells, jycells, jzcells, ncells, &
      nx, ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                 :: nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
      INTEGER(idp), INTENT(IN)                 :: nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN), DIMENSION(8, ncells):: jxcells, jycells, jzcells
    END SUBROUTINE

    SUBROUTINE current_reduction_3_3_3(jx, jy, jz, jxcells, jycells, jzcells, ncells, &
      nx, ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                 :: nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
      INTEGER(idp), INTENT(IN)                 :: nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN), DIMENSION(8, ncells):: jxcells, jycells, jzcells
    END SUBROUTINE

  END INTERFACE
  ! ___________________________________________________________________________

#if defined(DEBUG)
  WRITE(0, *) "Depose_currents_on_grid: start"
#endif

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==2
  CALL start_vtune_collection()
#endif
#if SDE==2
  CALL start_vtune_collection()
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  jx = 0.0_num
  jy = 0.0_num
  jz = 0.0_num

  ! Current deposition branches
  ! _______________________________________________________
  ! Classical current deposition, non-optimized/no tiling
  IF (currdepo.EQ.5) THEN

    IF ((nox.eq.noy).AND.(noy.eq.noz)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(depose_jxjyjz, jx, jy, &
      jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, &
      3_idp)
      ! The last argument is 3:
      ! this means the scalar routines will be used inside `depose_jxjyjz`
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, jx,  &
      jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, &
      dt, 1_idp)
      ! The last argument is 1:
      ! this means the generic esirkepov routine will be used inside `depose_jxjyjz`
    ENDIF


    ! _______________________________________________________
    ! Classical current deposition, non-optimized/tiling
  ELSE IF (currdepo.EQ.4) THEN

    IF ((nox.eq.noy).AND.(noy.eq.noz)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz, jx,  &
      jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, &
      dt, 3_idp)
      ! The last argument is 3:
      ! this means the scalar routines will be used inside `depose_jxjyjz`
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, jx,  &
      jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, &
      dt, 1_idp)
      ! The last argument is 1:
      ! this means the generic esirkepov routine will be used inside `depose_jxjyjz`
    ENDIF

    ! _______________________________________________________
    ! Classical current deposition, parallel, vectorized
  ELSE IF (currdepo.EQ.3) THEN

    IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
      ! Version with reduction for each species
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2(                 &
      depose_jxjyjz_vecHV_vnr_3_3_3, current_reduction_3_3_3, jx, jy, jz, nx, ny, nz, &
      nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, lvec_curr_depo)
    ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN
      ! Version with reduction for each species
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2(                 &
      depose_jxjyjz_vecHV_vnr_2_2_2, current_reduction_2_2_2, jx, jy, jz, nx, ny, nz, &
      nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, lvec_curr_depo)
    ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN
      ! Version with reduction for each species
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2(                 &
      depose_jxjyjz_vecHV_vnr_1_1_1, current_reduction_1_1_1, jx, jy, jz, nx, ny, nz, &
      nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, lvec_curr_depo)
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, jx,  &
      jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, &
      dt, 1_idp)
    ENDIF
    ! _______________________________________________________
    ! Esirkepov sequential version
  ELSE IF (currdepo.EQ.2) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq(jx, jy, jz, nx, ny, nz,  &
    nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt)

    ! _______________________________________________________
    ! Esirkepov tiling version
  ELSE IF (currdepo.EQ.1) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, jx,    &
    jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz,   &
    dt, 0_idp)
    ! The last argument is 0:
    ! this means the optimized esirkepov routine will be used inside `depose_jxjyjz`

    ! _______________________________________________________
    ! Default - Esirkepov parallel version with OPENMP/tiling and optimizations
  ELSE IF (currdepo .EQ. 0) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, jx,    &
    jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz,   &
    dt, 0_idp)
    ! The last argument is 1:
    ! this means the optimized esirkepov routine will be used inside `depose_jxjyjz`

    ! Arbitrary order
  ELSE
    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, jx,    &
    jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz,   &
    dt, 1_idp)
    ! The last argument is 1:
    ! this means the generic esirkepov routine will be used inside `depose_jxjyjz`
  ENDIF

  !!! --- Stop Vtune analysis
#if VTUNE==2
  CALL stop_vtune_collection()
#endif
#if SDE==2
  CALL stop_sde_collection()
#endif

  IF (it.ge.timestat_itstart) THEN
    tend = MPI_WTIME()
    localtimes(3)=localtimes(3)+(tend-tdeb)
  ENDIF

#if defined(DEBUG)
  WRITE(0, *) "Depose_current_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz

! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile with the classical method using an external given function
!
!> @details
!> OpenMP version. Avoids conflict while reducing tile currents in the global
!> current array.
!> This subroutine uses an external function represented by the argument func_order
!> for the current deposition method.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(func_order, jxg,    &
  jyg, jzg, nxx, nyy, nzz, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy,    &
  dzz, dtt, current_depo_algo)
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  IMPLICIT NONE
  ! Interfaces for func_order
  INTERFACE
    SUBROUTINE func_order(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,    &
      xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,     &
      noy, noz, current_depo_algo)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard,       &
      current_depo_algo
      REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                  &
      -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
    END SUBROUTINE
  END INTERFACE

  ! Parameters
  INTEGER(idp), INTENT(IN)  :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN)  :: noxx, noyy, nozz
  REAL(num), INTENT(IN)     :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)              :: ispecies, ix, iy, iz, count
  INTEGER(idp)              :: current_depo_algo
  INTEGER(idp)              :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)              :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                     :: isdeposited=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, current_depo_algo, nzjguard, dxx, dyy,    &
  !$OMP dzz, dtt, jxg, jyg, jzg, noxx, noyy, nozz, aofgrid_tiles) PRIVATE(ix, iy, iz, &
  !$OMP ispecies, curr, currg, curr_tile, count, jmin, jmax, kmin, kmax, lmin, lmax,  &
  !$OMP jminc, jmaxc, kminc, kmaxc, lminc, lmaxc, nxc, nyc, nzc, nxjg, nyjg, nzjg,    &
  !$OMP isdeposited)
  !! Current deposition
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        nxjg=curr_tile%nxg_tile
        nyjg=curr_tile%nyg_tile
        nzjg=curr_tile%nzg_tile
        jmin=curr_tile%nx_tile_min
        jmax=curr_tile%nx_tile_max
        kmin=curr_tile%ny_tile_min
        kmax=curr_tile%ny_tile_max
        lmin=curr_tile%nz_tile_min
        lmax=curr_tile%nz_tile_max
        nxc=curr_tile%nx_cells_tile
        nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0.
        currg%arr2=0.
        currg%arr3=0.
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF
          ! Depose current in jtile
          CALL func_order(currg%arr1, currg%arr2, currg%arr3, count,            &
          curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,    &
          curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv,                &
          curr_tile%pid(1, wpid), curr%charge, curr_tile%x_grid_tile_min,             &
          curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz,   &
          nxc, nyc, nzc, nxjg, nyjg, nzjg, noxx, noyy, nozz, current_depo_algo )
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
          jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr2(0:nxc, 0:nyc, 0:nzc)
          jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr3(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- X
          jxg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jxg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jxg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jxg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- X
          jyg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jyg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jyg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jyg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- X
          jzg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jzg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jzg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jzg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Y
          jxg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jxg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jxg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jxg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Y
          jyg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jyg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jyg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jyg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Y
          jzg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jzg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jzg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jzg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Z
          jxg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jxg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          jxg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jxg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Z
          jyg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jyg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr2(0:nxc, 0:nyc, -nzjg:-1)
          jyg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jyg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr2(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Z
          jzg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jzg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr3(0:nxc, 0:nyc, -nzjg:-1)
          jzg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jzg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr3(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp


! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile with the classical method with an external given function.
!> In this second version, the transient current arrays are reduced
!> after the current deposition for all species and not for each species.
!
!> @details
!> OpenMP version. Avoids conflict while reducing tile currents in the global
!> current array.
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!
!> @date
!> 2016
!
!> @param[in] func_order represent the subroutine to be used for current
!> deposition depending on the selected order
!> @param[in] jxg, jyg, jzg current arrays
!> @param[in] nxx, nyy, nzz cell number in each direction
!> @param[in] nxjguard, nyjguard, nzjguard guard cells
!> @param[in] noxx, noyy, nozz orders for current deposition
!> @param[in] dxx, dyy, dzz, dtt space and time steps
!> @param[in] lvect vector size
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( curr_depo_sub,  &
  curr_reduc_sub, jxg, jyg, jzg, nxx, nyy, nzz, nxjguard, nyjguard, nzjguard, noxx,     &
  noyy, nozz, dxx, dyy, dzz, dtt, lvect)
  USE fields, ONLY: l_nodalgrid
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  IMPLICIT NONE
  ! _______________________________________________________________________
  ! Interfaces for curr_depo_sub and curr_reduc_sub
  INTERFACE
    SUBROUTINE curr_depo_sub(jxcells, jycells, jzcells, np, ncells, xp, yp, zp, uxp,  &
      uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard,    &
      nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
      REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
      REAL(num), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      REAL(num), INTENT(IN)                         :: q, dt, dx, dy, dz, xmin, ymin, &
      zmin
      INTEGER(idp)                                  :: ncx, ncy, ncz
      INTEGER(idp)                                  :: lvect
      LOGICAL(lp), INTENT(IN)                       :: l_nodal
    END SUBROUTINE

    SUBROUTINE curr_reduc_sub(jx, jy, jz, jxcells, jycells, jzcells, ncells, nx, ny,  &
      nz, nxguard, nyguard, nzguard, ncx, ncy, ncz) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                 :: nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
      INTEGER(idp), INTENT(IN)                 :: nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN), DIMENSION(8, ncells):: jxcells, jycells, jzcells
    END SUBROUTINE

  END INTERFACE

  ! _______________________________________________________________________
  ! Parameters
  INTEGER(idp), INTENT(IN)               :: nxx, nyy, nzz, nxjguard, nyjguard,        &
  nzjguard
  INTEGER(idp), INTENT(IN)               :: noxx, noyy, nozz
  INTEGER(idp), INTENT(IN)               :: lvect
  REAL(num), INTENT(IN)                  :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: jxcells, jycells, jzcells
  INTEGER(idp)                           :: ispecies, ix, iy, iz, np, ncells
  INTEGER(idp)                           :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)                           :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  INTEGER(idp)                           :: ncz, ncy, ncx
  TYPE(particle_species), POINTER        :: curr
  TYPE(particle_tile), POINTER           :: curr_tile
  TYPE(grid_tile), POINTER               :: currg
  INTEGER(idp)                           :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                            :: isdeposited=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN
  ! _______________________________________________________________________
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, jxg, jyg,   &
  !$OMP jzg, noxx, noyy, nozz, aofgrid_tiles) FIRSTPRIVATE(lvect, l_nodalgrid) PRIVATE(ix, iy, iz, &
  !$OMP ispecies, ncells, curr, currg, curr_tile, np, jmin, jmax, kmin, kmax, lmin,   &
  !$OMP lmax, jminc, jmaxc, kminc, kmaxc, lminc, lmaxc, nxc, nyc, nzc, nxjg, nyjg,    &
  !$OMP nzjg, isdeposited, jxcells, jycells, jzcells, ncx, ncy, ncz)
  !! Current deposition
  ! Loop on the tiles
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        nxjg=curr_tile%nxg_tile
        nyjg=curr_tile%nyg_tile
        nzjg=curr_tile%nzg_tile
        jmin=curr_tile%nx_tile_min
        jmax=curr_tile%nx_tile_max
        kmin=curr_tile%ny_tile_min
        kmax=curr_tile%ny_tile_max
        lmin=curr_tile%nz_tile_min
        lmax=curr_tile%nz_tile_max
        nxc=curr_tile%nx_cells_tile
        nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0.
        currg%arr2=0.
        currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
        isdeposited=.FALSE.

        IF ((noxx.eq.3).and.(noyy.eq.3).and.(nozz.eq.3))  THEN
          ncx=nxc+5; ncy=nyc+4; ncz=nzc+3
        ELSE IF ((noxx.eq.2).and.(noyy.eq.2).and.(nozz.eq.2)) THEN
          ! Originally was ncx=nxc+4; ncy=nyc+4; ncz=nzc+4
          ! But we need one more for the algorithm at order 2
          ncx=nxc+5; ncy=nyc+5; ncz=nzc+5
        ELSE IF ((noxx.eq.1).and.(noyy.eq.1).and.(nozz.eq.1)) THEN
          ncx=nxc+3; ncy=nyc+3; ncz=nzc+3
        ENDIF

        ncells = ncx*ncy*ncz
        ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
        jxcells=0.0_num
        jycells=0.0_num
        jzcells=0.0_num

        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF

          ! Depose current in jtile
          CALL curr_depo_sub(jxcells, jycells, jzcells, np, ncells, curr_tile%part_x, &
          curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy,   &
          curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%pid(1, wpid),           &
          curr%charge, curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,          &
          curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,   &
          nzjg, ncx, ncy, ncz, lvect, l_nodalgrid)

          !print*, sum(jxcells), sum(jycells), sum(jzcells)

        END DO! END LOOP ON SPECIES

        CALL curr_reduc_sub(currg%arr1, currg%arr2, currg%arr3, jxcells,        &
        jycells, jzcells, ncells, nxc, nyc, nzc, nxjg, nyjg, nzjg, ncx, ncy, ncz)

        DEALLOCATE(jxcells, jzcells, jycells)

        IF (isdeposited) THEN
          jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
          jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr2(0:nxc, 0:nyc, 0:nzc)
          jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr3(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- X
          jxg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jxg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jxg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jxg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- X
          jyg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jyg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jyg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jyg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- X
          jzg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jzg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jzg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jzg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Y
          jxg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jxg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jxg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jxg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Y
          jyg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jyg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jyg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jyg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Y
          jzg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jzg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jzg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jzg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Z
          jxg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jxg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          jxg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jxg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Z
          jyg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jyg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr2(0:nxc, 0:nyc, -nzjg:-1)
          jyg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jyg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr2(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Z
          jzg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jzg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr3(0:nxc, 0:nyc, -nzjg:-1)
          jzg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jzg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr3(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2

! ________________________________________________________________________________________
!
!> @brief
!> Deposit current in each tile with the classical method version 3
!
!> @details
!> OpenMP version. Avoids conflict while reducing tile currents in the global
!> current array.
!> In this second version, the transient current arrays are reduced
!> after the current deposition for all species and not for each species.
!> The loop over the species is also firt loop and not inside the tile loops.
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!
!> @param[in] func_order represent the subroutine to be used for current
!> deposition depending on the selected order
!> @param[in] curr_reduc_sub subroutine to be used for the reduction
!> @param[in] jxg, jyg, jzg current arrays
!> @param[in] nxx, nyy, nzz cell number in each direction
!> @param[in] nxjguard, nyjguard, nzjguard guard cells
!> @param[in] noxx, noyy, nozz orders for current deposition
!> @param[in] dxx, dyy, dzz, dtt space and time steps
!> @param[in] lvect vector size
!
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v3( curr_depo_sub,  &
  curr_reduc_sub, jxg, jyg, jzg, nxx, nyy, nzz, nxjguard, nyjguard, nzjguard, noxx,     &
  noyy, nozz, dxx, dyy, dzz, dtt, lvect)
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  ! ______________________________________________________________________________
  IMPLICIT NONE

  ! _______________________________________________________________________
  ! Interfaces for curr_depo_sub and curr_reduc_sub
  INTERFACE
    SUBROUTINE curr_depo_sub(jxcells, jycells, jzcells, np, ncells, xp, yp, zp, uxp,  &
      uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard,    &
      nyguard, nzguard, ncx, ncy, ncz, lvect)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
      REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
      REAL(num), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      REAL(num), INTENT(IN)                         :: q, dt, dx, dy, dz, xmin, ymin, &
      zmin
      INTEGER(idp)                                  :: ncx, ncy, ncz
      INTEGER(idp)                                  :: lvect
    END SUBROUTINE

    SUBROUTINE curr_reduc_sub(jx, jy, jz, jxcells, jycells, jzcells, ncells, nx, ny,  &
      nz, nxguard, nyguard, nzguard, ncx, ncy, ncz) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT(IN)                 :: nx, ny, nz, ncells
      INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
      INTEGER(idp), INTENT(IN)                 :: nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN OUT) ::                                                    &
      jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT(IN), DIMENSION(8, ncells):: jxcells, jycells, jzcells
    END SUBROUTINE

  END INTERFACE

  ! _______________________________________________________________________
  ! Parameters
  INTEGER(idp), INTENT(IN)               :: nxx, nyy, nzz, nxjguard, nyjguard,        &
  nzjguard
  INTEGER(idp), INTENT(IN)               :: noxx, noyy, nozz
  INTEGER(idp), INTENT(IN)               :: lvect
  REAL(num), INTENT(IN)                  :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: jxcells, jycells, jzcells
  INTEGER(idp)                           :: ispecies, ix, iy, iz, np, ncells
  INTEGER(idp)                           :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)                           :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  INTEGER(idp)                           :: ncz, ncy, ncx
  TYPE(particle_species), POINTER        :: curr
  TYPE(particle_tile), POINTER           :: curr_tile
  TYPE(grid_tile), POINTER               :: currg
  INTEGER(idp)                           :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                            :: isdeposited=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN
  ! _______________________________________________________________________
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, jxg, jyg,   &
  !$OMP jzg, noxx, noyy, nozz, aofgrid_tiles) FIRSTPRIVATE(lvect) PRIVATE(ix, iy, iz, &
  !$OMP ispecies, ncells, curr, currg, curr_tile, np, jmin, jmax, kmin, kmax, lmin,   &
  !$OMP lmax, jminc, jmaxc, kminc, kmaxc, lminc, lmaxc, nxc, nyc, nzc, nxjg, nyjg,    &
  !$OMP nzjg, isdeposited, jxcells, jycells, jzcells, ncx, ncy, ncz)
  DO ispecies=1, nspecies! LOOP ON SPECIES

    !! Current deposition
    ! Loop on the tiles
    !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
    DO iz=1, ntilez
      DO iy=1, ntiley
        DO ix=1, ntilex
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jmin=curr_tile%nx_tile_min
          jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min
          kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min
          lmax=curr_tile%nz_tile_max
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          currg=>aofgrid_tiles(ix, iy, iz)
          currg%arr1=0.
          currg%arr2=0.
          currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
          isdeposited=.FALSE.

          IF ((noxx.eq.3).and.(noyy.eq.3).and.(nozz.eq.3))  THEN
            ncx=nxc+5; ncy=nyc+4; ncz=nzc+3
          ELSE IF ((noxx.eq.2).and.(noyy.eq.2).and.(nozz.eq.2)) THEN
            ! Originally was ncx=nxc+4; ncy=nyc+4; ncz=nzc+4
            ! But we need one more for the algorithm at order 2
            ncx=nxc+5; ncy=nyc+5; ncz=nzc+5
          ELSE IF ((noxx.eq.1).and.(noyy.eq.1).and.(nozz.eq.1)) THEN
            ncx=nxc+3; ncy=nyc+3; ncz=nzc+3
          ENDIF

          ncells = ncx*ncy*ncz
          ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
          jxcells=0.0_num
          jycells=0.0_num
          jzcells=0.0_num

          IF (np .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF

          ! Depose current in jtile
          CALL curr_depo_sub(jxcells, jycells, jzcells, np, ncells, curr_tile%part_x, &
          curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy,   &
          curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%pid(1, wpid),           &
          curr%charge, curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,          &
          curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,   &
          nzjg, ncx, ncy, ncz, lvect)

          !print*, sum(jxcells), sum(jycells), sum(jzcells)


          CALL curr_reduc_sub(currg%arr1, currg%arr2, currg%arr3, jxcells,      &
          jycells, jzcells, ncells, nxc, nyc, nzc, nxjg, nyjg, nzjg, ncx, ncy, ncz)

          DEALLOCATE(jxcells, jzcells, jycells)

          IF (isdeposited) THEN
            jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,            &
            lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
            jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,            &
            lmin:lmax)+currg%arr2(0:nxc, 0:nyc, 0:nzc)
            jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,            &
            lmin:lmax)+currg%arr3(0:nxc, 0:nyc, 0:nzc)
          ENDIF
        END DO
      END DO
    END DO!END LOOP ON TILES
    !$OMP END DO

  END DO! END LOOP ON SPECIES

  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- X
          jxg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jxg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jxg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jxg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- X
          jyg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jyg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jyg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jyg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- X
          jzg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jzg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jzg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jzg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Y
          jxg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jxg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jxg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jxg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Y
          jyg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jyg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jyg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jyg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Y
          jzg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jzg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jzg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jzg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np=curr_tile%np_tile(1)
          IF (np .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Z
          jxg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jxg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          jxg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jxg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Z
          jyg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jyg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr2(0:nxc, 0:nyc, -nzjg:-1)
          jyg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jyg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr2(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Z
          jzg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jzg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr3(0:nxc, 0:nyc, -nzjg:-1)
          jzg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jzg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr3(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v3

! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile with Esirkepov method
!
!> @details
!> This subroutine is called from Fortram main program and contains an interface argument
!> OpenMP version. Avoids conflict while reducing tile currents in the global
!> current array.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2016
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(func_order, jxg,    &
  jyg, jzg, nxx, nyy, nzz, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy,    &
  dzz, dtt, current_depo_algo)
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  REAL(num), INTENT(IN) :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp) :: ispecies, ix, iy, iz, count, current_depo_algo
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER :: curr_tile
  TYPE(grid_tile), POINTER :: currg
  INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)  :: isdeposited=.FALSE.

  ! Interfaces for func_order
  INTERFACE
    SUBROUTINE func_order(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,    &
      xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,     &
      noy, noz, current_depo_algo)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard,       &
      current_depo_algo
      REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                  &
      -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
    END SUBROUTINE
  END INTERFACE

  IF (nspecies .EQ. 0_idp) RETURN

  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, current_depo_algo, nzjguard, dxx, dyy,    &
  !$OMP dzz, dtt, jxg, jyg, jzg, noxx, noyy, nozz, aofgrid_tiles, c_dim) PRIVATE(ix,  &
  !$OMP iy, iz, ispecies, curr, currg, curr_tile, count, jmin, jmax, kmin, kmax,      &
  !$OMP lmin, lmax, jminc, jmaxc, kminc, kmaxc, lminc, lmaxc, nxc, nyc, nzc, nxjg,    &
  !$OMP nyjg, nzjg, isdeposited)
  !! Current deposition
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        nxjg=curr_tile%nxg_tile
        nyjg=curr_tile%nyg_tile
        nzjg=curr_tile%nzg_tile
        jmin=curr_tile%nx_tile_min
        jmax=curr_tile%nx_tile_max
        kmin=curr_tile%ny_tile_min
        kmax=curr_tile%ny_tile_max
        lmin=curr_tile%nz_tile_min
        lmax=curr_tile%nz_tile_max
        nxc=curr_tile%nx_cells_tile;
        nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0.
        currg%arr2=0.
        currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF

          ! Depose current in jtile
          CALL func_order(currg%arr1, currg%arr2, currg%arr3, count,            &
          curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,    &
          curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv,                &
          curr_tile%pid(1, wpid), curr%charge, curr_tile%x_grid_tile_min,             &
          curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz,   &
          nxc, nyc, nzc, nxjg, nyjg, nzjg, noxx, noyy, nozz, current_depo_algo)

        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
          jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr2(0:nxc, 0:nyc, 0:nzc)
          jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr3(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- X
          jxg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jxg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jxg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jxg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- X
          jyg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jyg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jyg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jyg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- X
          jzg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jzg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jzg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jzg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Y
          jxg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jxg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jxg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jxg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Y
          jyg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jyg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jyg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jyg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Y
          jzg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jzg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jzg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jzg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Z
          jxg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jxg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          jxg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jxg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Z
          jyg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jyg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr2(0:nxc, 0:nyc, -nzjg:-1)
          jyg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jyg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr2(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Z
          jzg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jzg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr3(0:nxc, 0:nyc, -nzjg:-1)
          jzg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jzg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr3(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp

! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile with Esirkepov method
!
!> @details
!> This subroutine is called from Python and does not have interface arguments
!> OpenMP version. Avoids conflict while reducing tile currents in the global
!> current array.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jxg, jyg, jzg, nxx, nyy, nzz, &
  nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt)
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  IMPLICIT NONE

  INTEGER(idp), INTENT(IN)        :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN)        :: noxx, noyy, nozz
  REAL(num), INTENT(IN)           :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT)       :: jxg(-nxjguard:nxx+nxjguard,                      &
  -nyjguard:nyy+nyjguard, -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT)       :: jyg(-nxjguard:nxx+nxjguard,                      &
  -nyjguard:nyy+nyjguard, -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT)       :: jzg(-nxjguard:nxx+nxjguard,                      &
  -nyjguard:nyy+nyjguard, -nzjguard:nzz+nzjguard)
  INTEGER(idp)                    :: ispecies, ix, iy, iz, count
  INTEGER(idp)                    :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)                    :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg

  INTEGER(idp)                    :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                     :: isdeposited=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN

  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, jxg, jyg,   &
  !$OMP jzg, noxx, noyy, nozz, aofgrid_tiles, c_dim) PRIVATE(ix, iy, iz, ispecies,    &
  !$OMP curr, currg, curr_tile, count, jmin, jmax, kmin, kmax, lmin, lmax, jminc,     &
  !$OMP jmaxc, kminc, kmaxc, lminc, lmaxc, nxc, nyc, nzc, nxjg, nyjg, nzjg,           &
  !$OMP isdeposited)
  !! Current deposition
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        nxjg=curr_tile%nxg_tile
        nyjg=curr_tile%nyg_tile
        nzjg=curr_tile%nzg_tile
        jmin=curr_tile%nx_tile_min
        jmax=curr_tile%nx_tile_max
        kmin=curr_tile%ny_tile_min
        kmax=curr_tile%ny_tile_max
        lmin=curr_tile%nz_tile_min
        lmax=curr_tile%nz_tile_max
        nxc=curr_tile%nx_cells_tile;
        nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0.
        currg%arr2=0.
        currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF
          ! Depose current in jtile
          SELECT CASE (c_dim)
          CASE (2)
            CALL depose_jxjyjz_esirkepov_2d( currg%arr1(:, 0, :), currg%arr2(:,   &
            0, :), currg%arr3(:, 0, :), count, curr_tile%part_x, curr_tile%part_y,  &
            curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy,                   &
            curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%pid(1, wpid),         &
            curr%charge, curr_tile%x_grid_tile_min, curr_tile%z_grid_tile_min, dtt,   &
            dxx, dzz, nxc, nzc, nxjg, nzjg, noxx, nozz)
          CASE DEFAULT
            CALL depose_jxjyjz_esirkepov( currg%arr1, currg%arr2, currg%arr3,   &
            count, curr_tile%part_x, curr_tile%part_y, curr_tile%part_z,              &
            curr_tile%part_ux, curr_tile%part_uy, curr_tile%part_uz,                  &
            curr_tile%part_gaminv, curr_tile%pid(1, wpid), curr%charge,               &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg, &
            nzjg, noxx, noyy, nozz)
          END SELECT
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
          jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr2(0:nxc, 0:nyc, 0:nzc)
          jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr3(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- X
          jxg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jxg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jxg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jxg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- X
          jyg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jyg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jyg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jyg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr2(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- X
          jzg(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = jzg(jminc:jmin-1,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(-nxjg:-1, -nyjg:nyc+nyjg,           &
          -nzjg:nzc+nzjg)
          jzg(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = jzg(jmax+1:jmaxc,             &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr3(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,     &
          -nzjg:nzc+nzjg)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Y
          jxg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jxg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jxg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jxg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Y
          jyg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jyg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jyg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jyg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr2(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Y
          jzg(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = jzg(jmin:jmax, kminc:kmin-1,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          jzg(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = jzg(jmin:jmax, kmax+1:kmaxc,    &
          lminc:lmaxc)+ currg%arr3(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE.
        END DO
        IF (isdeposited) THEN
          currg=>aofgrid_tiles(ix, iy, iz)
          curr => species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
          kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
          lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
          nxjg=curr_tile%nxg_tile
          nyjg=curr_tile%nyg_tile
          nzjg=curr_tile%nzg_tile
          jminc=jmin-nxjg; jmaxc=jmax+nxjg
          kminc=kmin-nyjg; kmaxc=kmax+nyjg
          lminc=lmin-nzjg; lmaxc=lmax+nzjg
          nxc=curr_tile%nx_cells_tile
          nyc=curr_tile%ny_cells_tile
          nzc=curr_tile%nz_cells_tile
          ! ----- Add guardcells in adjacent tiles
          ! --- JX
          ! - FACES +/- Z
          jxg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jxg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          jxg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jxg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JY
          ! - FACES +/- Z
          jyg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jyg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr2(0:nxc, 0:nyc, -nzjg:-1)
          jyg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jyg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr2(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
          ! --- JZ
          ! - FACES +/- Z
          jzg(jmin:jmax, kmin:kmax, lminc:lmin-1) = jzg(jmin:jmax, kmin:kmax,         &
          lminc:lmin-1)+ currg%arr3(0:nxc, 0:nyc, -nzjg:-1)
          jzg(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = jzg(jmin:jmax, kmin:kmax,         &
          lmax+1:lmaxc)+ currg%arr3(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_sub_openmp

! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile sequentially
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(func_order, jxg, jyg,  &
  jzg, nxx, nyy, nzz, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz,    &
  dtt, currrent_depo_algo )
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling

  IMPLICIT NONE
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  REAL(num), INTENT(IN) :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)                    :: ispecies, ix, iy, iz, count
  INTEGER(idp)                    :: currrent_depo_algo
  INTEGER(idp)                    :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                     :: isdeposited=.FALSE.

  ! Interfaces for func_order
  INTERFACE
    SUBROUTINE func_order(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,    &
      xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,     &
      noy, noz, current_depo_algo)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard,       &
      current_depo_algo
      REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                  &
      -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
    END SUBROUTINE
  END INTERFACE

  IF (nspecies .EQ. 0_idp) RETURN
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        nxjg=curr_tile%nxg_tile
        nyjg=curr_tile%nyg_tile
        nzjg=curr_tile%nzg_tile
        jmin=curr_tile%nx_tile_min-nxjg
        jmax=curr_tile%nx_tile_max+nxjg
        kmin=curr_tile%ny_tile_min-nyjg
        kmax=curr_tile%ny_tile_max+nyjg
        lmin=curr_tile%nz_tile_min-nzjg
        lmax=curr_tile%nz_tile_max+nzjg
        nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0.
        currg%arr2=0.
        currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF
          ! Depose current in jtile
          CALL func_order(currg%arr1, currg%arr2, currg%arr3, count,            &
          curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,    &
          curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv,                &
          curr_tile%pid(1, wpid), curr%charge, curr_tile%x_grid_tile_min,             &
          curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz,   &
          nxc, nyc, nzc, nxjg, nyjg, nzjg, noxx, noyy, nozz, currrent_depo_algo)
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr1
          jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr2
          jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr3
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES


END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq


! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile
!> Sequential version
!
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq(jxg, jyg, jzg, nxx,    &
  nyy, nzz, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt)
USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
USE particle_properties, ONLY: nspecies, wpid
USE particle_speciesmodule, ONLY: particle_species
USE particle_tilemodule, ONLY: particle_tile
USE particles, ONLY: species_parray
USE picsar_precision, ONLY: idp, lp, num
USE tile_params, ONLY: ntilex, ntiley, ntilez
USE tiling


  IMPLICIT NONE
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  REAL(num), INTENT(IN) :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER :: curr_tile
  TYPE(grid_tile), POINTER :: currg
  INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)  :: isdeposited=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        nxjg=curr_tile%nxg_tile
        nyjg=curr_tile%nyg_tile
        nzjg=curr_tile%nzg_tile
        jmin=curr_tile%nx_tile_min-nxjg
        jmax=curr_tile%nx_tile_max+nxjg
        kmin=curr_tile%ny_tile_min-nyjg
        kmax=curr_tile%ny_tile_max+nyjg
        lmin=curr_tile%nz_tile_min-nzjg
        lmax=curr_tile%nz_tile_max+nzjg
        nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0.
        currg%arr2=0.
        currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF
          ! Depose current in jtile
          CALL depose_jxjyjz_esirkepov( currg%arr1, currg%arr2, currg%arr3,     &
          count, curr_tile%part_x, curr_tile%part_y, curr_tile%part_z,                &
          curr_tile%part_ux, curr_tile%part_uy, curr_tile%part_uz,                    &
          curr_tile%part_gaminv, curr_tile%pid(1, wpid), curr%charge,                 &
          curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                       &
          curr_tile%z_grid_tile_min, dtt, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,   &
          nzjg, noxx, noyy, nozz)
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          jxg(jmin:jmax, kmin:kmax, lmin:lmax)=jxg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr1
          jyg(jmin:jmax, kmin:kmax, lmin:lmax)=jyg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr2
          jzg(jmin:jmax, kmin:kmax, lmin:lmax)=jzg(jmin:jmax, kmin:kmax,              &
          lmin:lmax)+currg%arr3
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq

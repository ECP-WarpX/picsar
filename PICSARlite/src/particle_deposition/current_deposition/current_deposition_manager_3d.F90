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
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq
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
  noz)
  USE constants
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

  call depose_jxjyjz_generic( jx, nguard, nvalid, jy, nguard, nvalid, jz, nguard,     &
  nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy,  &
  dz, nox, noy, noz)
END SUBROUTINE

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
  zmin, dt, dx, dy, dz, nox, noy, noz)     !#do not wrap
  USE constants
  IMPLICIT NONE
  INTEGER(idp) :: np, nox, noy, noz
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

    ! Scalar classical current deposition subroutines
!       CALL depose_jxjyjz_scalar_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
!       jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
!       q, xmin, ymin, zmin, dt, dx, dy, dz)

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      CALL depose_jxjyjz_scalar_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz)
    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      CALL depose_jxjyjz_scalar_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz)
    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      CALL depose_jxjyjz_scalar_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
      jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
      q, xmin, ymin, zmin, dt, dx, dy, dz)
!     ELSE
!       CALL pxr_depose_jxjyjz_esirkepov_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
!       jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
!       q, xmin, ymin, zmin, dt, dx, dy, dz, nox, noy, noz, .TRUE._idp, .FALSE._idp)
    ENDIF


END SUBROUTINE

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
  USE fields
  USE particles
  USE shared_data
  USE params
  USE time_stat
#if defined(VTUNE) && VTUNE==2
  USE ITT_FORTRAN
#endif
#if defined(SDE) && SDE==2
  USE SDE_FORTRAN
#endif
  IMPLICIT NONE
  REAL(num) :: tdeb, tend

  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE
    ! ____________________________________________________________________________________
    ! Generic current deposition routine
    SUBROUTINE depose_jxjyjz(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
      xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,     &
      noy, noz)  !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard
      REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                  &
      -nzguard:nz+nzguard), intent(in out) :: jx, jy, jz
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
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
  IF ((nox.eq.noy).AND.(noy.eq.noz).AND.(nox>0.and.nox<4)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(depose_jxjyjz, jx, jy, &
      jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt)
  ELSE
       PRINT *, 'CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp'
       STOP
  ENDIF
    ! _______________________________________________________
    ! Classical current deposition, non-optimized/tiling

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
  dtt)
  USE particles
  USE constants
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
      noy, noz)  !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard
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
          nxc, nyc, nzc, nxjg, nyjg, nzjg, noxx, noyy, nozz)
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

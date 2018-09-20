! ________________________________________________________________________________________
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
! CURRENT_DEPOSITION_MANAGER_2D.F90
!
! Developers
! Henri Vincenti, ! Mathieu Lobet
!
! Description:
! This file contains subroutines for managing the current deposition in 2D.
!
! List of subroutines
! - pxrdepose_currents_on_grid_jxjyjz_2d
!
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Generic subroutine for current deposition on one tile
!>
!> @details
!> This routine calls the relevant current deposition routine depending
!> on the order of the particle shape and the selected algorithm.
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_2d(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,  &
  xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz, lvect)
  USE picsar_precision, ONLY: idp, num
  implicit none
  integer(idp)                          :: np, nx, nz, nox, noz, nxguard, nzguard
  integer(idp)                          :: lvect
  real(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard), intent(inout) ::    &
  jx, jy, jz
  real(num), dimension(np)              :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin

  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(2), nvalid(2)
  nguard = (/ nxguard, nzguard /)
  nvalid = (/ nx+1, nz+1 /)

  call depose_jxjyjz_generic_2d( jx, nguard, nvalid, jy, nguard, nvalid, jz, nguard,  &
  nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, zmin, dt, dx, dz, nox,   &
  noz, lvect)

END SUBROUTINE depose_jxjyjz_2d

! ______________________________________________________________________________
!> @brief
!> Generic subroutines for current deposition, adapted for field
!> arrays having different sizes depending on their nodal/cell-centered nature
!>
!> @details
!> This routine calls the relevant current deposition routine depending
!> on the order of the particle shape and the selected algorithm.
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_generic_2d( jx, jx_nguard, jx_nvalid, jy, jy_nguard,         &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, nox, noz, lvect)     !#do not wrap
  USE picsar_precision, ONLY: idp, num, lp
  implicit none
  integer(idp)                          :: np, nox, noz
  INTEGER(idp), intent(in)              :: jx_nguard(2), jx_nvalid(2), jy_nguard(2),  &
  jy_nvalid(2), jz_nguard(2), jz_nvalid(2)
  integer(idp)                          :: lvect
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1 )
  real(num), dimension(np)              :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin

  IF ((nox.eq.1).and.(noz.eq.1)) THEN
    CALL pxr_depose_jxjyjz_esirkepov2d_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,  &
    jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, q,     &
    xmin, zmin, dt, dx, dz, lvect, .TRUE._lp, .FALSE._lp, .FALSE._lp, 0_idp)
  ELSE IF ((nox.eq.2).and.(noz.eq.2)) THEN
    CALL pxr_depose_jxjyjz_esirkepov2d_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,  &
    jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, q,     &
    xmin, zmin, dt, dx, dz, lvect, .TRUE._lp, .FALSE._lp, .FALSE._lp, 0_idp)
  ELSE IF ((nox.eq.3).and.(noz.eq.3)) THEN
    CALL pxr_depose_jxjyjz_esirkepov2d_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,  &
    jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, q,     &
    xmin, zmin, dt, dx, dz, lvect, .TRUE._lp, .FALSE._lp, .FALSE._lp, 0_idp)
  ELSE
    CALL pxr_depose_jxjyjz_esirkepov2d_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
    jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
    xmin, zmin, dt, dx, dz, nox, noz, .TRUE._lp, .FALSE._lp, .FALSE._lp, 0_idp)
  ENDIF

END SUBROUTINE depose_jxjyjz_generic_2d

! ________________________________________________________________________________________
!> @brief
!> Esirkepov subroutine for current deposition on one tile
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_2d(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp,      &
  gaminv, w, q, xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz)
  USE picsar_precision, ONLY: idp, num, lp
  IMPLICIT NONE
  integer(idp)                          :: np, nx, nz, nox, noz, nxguard, nzguard
  real(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard), intent(inout) ::    &
  jx, jy, jz
  real(num), dimension(np)              :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin

  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(2), nvalid(2)
  nguard = (/ nxguard, nzguard /)
  nvalid = (/ nx+1, nz+1 /)

  CALL pxr_depose_jxjyjz_esirkepov2d_n( jx, nguard, nvalid, jy, nguard, nvalid, jz,   &
  nguard, nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, zmin, dt, dx,    &
  dz, nox, noz, .TRUE._lp, .FALSE._lp, .FALSE._lp, 0_idp)

END SUBROUTINE depose_jxjyjz_esirkepov_2d

! ________________________________________________________________________________________
!> @brief
!> Main subroutine for the current deposition called in submain in 2d.

!> @details
!> This subroutine determines which model and algorithm to use depending on the parameter
!> currdepo and the interpolation order.

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d
  USE fields, ONLY: nyjguards, nox, noy, noz, jz, nxjguards, nzjguards, jy, jx,      &
    l4symtry
  USE mpi
  USE params, ONLY: currdepo, lvec_curr_depo, dt
  USE picsar_precision, ONLY: idp, num, lp
  USE shared_data, ONLY: nz, ny, nx, xmin, zmin, dx, dy, dz
  USE time_stat, ONLY: localtimes
  IMPLICIT NONE

  ! __ Parameter declaration __________________________________________________
  REAL(num) :: tdeb, tend

  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE

    SUBROUTINE depose_jxjyjz_2d(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, &
      q, xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz, lvect)!#do not parse
      USE PICSAR_precision
      USE constants
      implicit none
      integer(idp)                          :: np, nx, nz, nox, noz, nxguard, nzguard
      integer(idp)                          :: lvect
      real(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard), intent(inout)   &
      :: jx, jy, jz
      real(num), dimension(np)              :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      real(num)                             :: q, dt, dx, dz, xmin, zmin
    END SUBROUTINE

#if defined(DEV)
    subroutine pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3(jx, jy, jz, np, xp, zp, uxp,   &
      uyp, uzp, gaminv, w, q, xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox,    &
      noz, lvect, l_particles_weight, l4symtry, l_2drz, type_rz_depose)    !#do not parse
      USE PICSAR_precision
      USE constants
      implicit none
      integer(idp)                          :: np, nx, nz, nox, noz, nxguard,         &
      nzguard, type_rz_depose
      integer(idp)                          :: lvect
      real(num), dimension((1+nx+2*nxguard)*(1+nz+2*nzguard)), intent(in out) :: jx,  &
      jy, jz
      real(num), dimension(np)              :: xp, zp, uxp, uyp, uzp, gaminv, w
      real(num)                             :: q, dt, dx, dz, xmin, zmin
      LOGICAL(lp)                           :: l_particles_weight, l4symtry, l_2drz
    END SUBROUTINE
#endif

  END INTERFACE
  ! ___________________________________________________________________________

  ! For debugging
#if defined(DEBUG)
  WRITE(0, *) "Depose_currents_on_grid: start"
#endif

  ! For time statistics
  tdeb=MPI_WTIME()

  ! For profiling with Vtune/SDE
#if PROFILING==2
  CALL start_collection()
#endif

  jx = 0.0_num
  jy = 0.0_num
  jz = 0.0_num

  ! __ Current deposition ________________________________________________________________

  ! _______________________________________________________
  ! Esirkepov general order subroutine
  IF (currdepo.EQ.2) THEN
    CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx, jy, jz, nx, ny, nz,         &
    nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt)
    ! _________________________________________________________________
    ! Esirkepov OpenMP/tiling version non-vectorized but more optimized
    ! than the general order subroutine
  ELSE IF (currdepo.EQ.1) THEN

    IF ((nox.eq.noz)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(depose_jxjyjz_2d, &
      jx, jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, &
      dz, dt, lvec_curr_depo)
      ! Order n
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx, jy, jz, nx, ny, nz,       &
      nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt)
    ENDIF

    ! _________________________________________________________________________
    ! Default - Esirkepov parallel version with OPENMP/tiling and optimizations
  ELSE

    IF ((nox.eq.noz)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(depose_jxjyjz_2d, &
      jx, jy, jz, nx, ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, &
      dz, dt, lvec_curr_depo)
      ! Order n
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx, jy, jz, nx, ny, nz,       &
      nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt)
    ENDIF

  ENDIF

  ! Stop Vtune/SDE analysis
#if PROFILING==2
  CALL stop_collection()
#endif

  ! For time statistics
  tend = MPI_WTIME()
  localtimes(3)=localtimes(3)+(tend-tdeb)

  ! For debugging
#if defined(DEBUG)
  WRITE(0, *) "Depose_current_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d



! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile with Esirkepov method in 2D
!
!> @details
!> This subroutine is called from Fortran main program and contains an interface argument
!> OpenMP version. Avoids conflict while reducing tile currents in the global
!> current array.

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(curr_depo_sub,    &
  jxg, jyg, jzg, nxx, nyy, nzz, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx,    &
  dyy, dzz, dtt, lvect)
  USE grid_tilemodule, ONLY: grid_tile, aofgrid_tiles
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, num, lp
  USE tile_params, ONLY: ntilez, ntilex, ntiley
  USE tiling
  IMPLICIT NONE
  ! __ Parameter declaration _____________________________________________________________
  INTEGER(idp), INTENT(IN)  :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN)  :: noxx, noyy, nozz
  INTEGER(idp), INTENT(IN)  :: lvect
  REAL(num), INTENT(IN)     :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)              :: ispecies, ix, iy, iz, count
  INTEGER(idp)              :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)              :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER  :: currg
  INTEGER(idp)              :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)               :: isdeposited=.FALSE.

  ! ___ Interface ________________________________________________________________________
  ! For the func_order input function
  INTERFACE
    SUBROUTINE curr_depo_sub(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
      xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz, lvect)!#do not parse !#do not parse !#do not parse
      USE PICSAR_precision
      USE constants
      implicit none
      integer(idp)                          :: np, nx, nz, nox, noz, nxguard, nzguard
      integer(idp)                          :: lvect
      real(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard), intent(inout)   &
      :: jx, jy, jz
      real(num), dimension(np)              :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
      real(num)                             :: q, dt, dx, dz, xmin, zmin
    END SUBROUTINE
  END INTERFACE

  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, jxg, jyg,   &
  !$OMP jzg, noxx, noyy, nozz, aofgrid_tiles, c_dim, lvect) PRIVATE(ix, iy, iz,       &
  !$OMP ispecies, curr, currg, curr_tile, count, jmin, jmax, kmin, kmax, lmin, lmax,  &
  !$OMP jminc, jmaxc, kminc, kmaxc, lminc, lmaxc, nxc, nyc, nzc, nxjg, nyjg, nzjg,    &
  !$OMP isdeposited)
  !! Current deposition
  !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO ix=1, ntilex
      curr => species_parray(1)
      curr_tile=>curr%array_of_tiles(ix, 1, iz)
      nxjg=curr_tile%nxg_tile
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
      currg=>aofgrid_tiles(ix, 1, iz)
      currg%arr1=0.
      currg%arr2=0.
      currg%arr3=0.!jzg(jmin:jmax, kmin:kmax, lmin:lmax)
      isdeposited=.FALSE.
      DO ispecies=1, nspecies! LOOP ON SPECIES
        curr => species_parray(ispecies)
        IF (.NOT. curr%ldodepos) CYCLE
        curr_tile=>curr%array_of_tiles(ix, 1, iz)
        count=curr_tile%np_tile(1)
        IF (count .EQ. 0) THEN
          CYCLE
        ELSE
          isdeposited=.TRUE.
        ENDIF

        ! Depose current in jtile
        CALL curr_depo_sub( currg%arr1(:,0,:), currg%arr2(:,0,:), currg%arr3(:,0,:), count,          &
        curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,      &
        curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%pid(1, &
        wpid), curr%charge, curr_tile%x_grid_tile_min, curr_tile%z_grid_tile_min,     &
        dtt, dxx, dzz, nxc, nzc, nxjg, nzjg, noxx, nozz, lvect)

      END DO! END LOOP ON SPECIES
      IF (isdeposited) THEN
        jxg(jmin:jmax, 0, lmin:lmax)=jxg(jmin:jmax, 0, lmin:lmax)+currg%arr1(0:nxc, &
        0, 0:nzc)
        jyg(jmin:jmax, 0, lmin:lmax)=jyg(jmin:jmax, 0, lmin:lmax)+currg%arr2(0:nxc, &
        0, 0:nzc)
        jzg(jmin:jmax, 0, lmin:lmax)=jzg(jmin:jmax, 0, lmin:lmax)+currg%arr3(0:nxc, &
        0, 0:nzc)
      ENDIF
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

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp

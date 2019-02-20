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
!
! CHARGE_DEPOSITION_MANAGER.F90
!
! Author
! Henri Vincenti, ! Mathieu Lobet
!
! Brief description:
! File containing subroutines to manage the charge deposition.
! These subroutines use the tiling and call the charge deposition subroutines in
! charge_deposition_2d.F90 and charge_deposition_3d.F90.
!
! List of suboutines:
!
! Main manager subroutine:
! - pxrdepose_rho_on_grid
!
! tile deposition:
! - pxrdepose_rho_on_grid_sub_openmp_2d
! - pxrdepose_rho_on_grid_sub_openmp_3d
! - pxrdepose_rho_on_grid_sub_openmp_3d_n
!
! - pxrdepose_rho_on_grid_sub_openmp_3d_scalar
! - pxrdepose_rho_on_grid_sub_openmp_3d_vecto
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Main subroutine for the charge deposition
!
!> @details
!> This subroutine is called in main.F90 and controle all the algorithm.
!> The parameter rhodepo enable to select a specific algorithm.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
!> last update 09/13/2016
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_rho_on_grid
  USE fields, ONLY: nox, noy, noz, nxjguards, nyjguards, nzjguards
  USE mpi
  USE params, ONLY: dt, it, lvec_charge_depo, rhodepo
  USE particle_properties, ONLY: nspecies
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: c_dim, dx, dy, dz, nx, ny, nz, rho, xmin, ymin, zmin
  USE time_stat, ONLY: localtimes, timestat_itstart
  IMPLICIT NONE

  INTEGER(idp) :: c_rho_old
  REAL(num)    :: tmptime

  INTERFACE

    SUBROUTINE depose_rho_scalar_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin,   &
      dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      INTEGER(idp), INTENT (IN) :: lvect
      REAL(num), INTENT (IN) :: q, dx, dy, dz, xmin, ymin, zmin
      REAL(num), INTENT (IN) :: xp(np), yp(np), zp(np), w(np)

    END SUBROUTINE

    SUBROUTINE depose_rho_scalar_2_2_2(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin,   &
      dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      INTEGER(idp), INTENT (IN) :: lvect
      REAL(num), INTENT (IN) :: q, dx, dy, dz, xmin, ymin, zmin
      REAL(num), INTENT (IN) :: xp(np), yp(np), zp(np), w(np)

    END SUBROUTINE

    SUBROUTINE depose_rho_scalar_3_3_3(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin,   &
      dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      INTEGER(idp), INTENT (IN) :: lvect
      REAL(num), INTENT (IN) :: q, dx, dy, dz, xmin, ymin, zmin
      REAL(num), INTENT (IN) :: xp(np), yp(np), zp(np), w(np)

    END SUBROUTINE

    SUBROUTINE depose_rho_vecHVv2_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin,  &
      dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      INTEGER(idp), INTENT (IN) :: lvect
      REAL(num), INTENT (IN) :: q, dx, dy, dz, xmin, ymin, zmin
      REAL(num), INTENT (IN) :: xp(np), yp(np), zp(np), w(np)

    END SUBROUTINE

    SUBROUTINE depose_rho_vecHVv2_2_2_2(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin,  &
      dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
      INTEGER(idp), INTENT (IN) :: lvect
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT (IN) :: xp(np), yp(np), zp(np), w(np)
      REAL(num), INTENT (IN) :: q, dx, dy, dz, xmin, ymin, zmin
    END SUBROUTINE

    SUBROUTINE depose_rho_vecHVv4_3_3_3(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin,  &
      dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
      INTEGER(idp), INTENT (IN) :: lvect
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), INTENT (IN) :: xp(np), yp(np), zp(np), w(np)
      REAL(num), INTENT (IN) :: q, dx, dy, dz, xmin, ymin, zmin

    END SUBROUTINE

  END INTERFACE

  IF (nspecies .EQ. 0_idp) RETURN
  ! ______________________________________
  ! Parameters

  ! For debugging
#if defined(DEBUG)
  WRITE(0, *) "pxrdepose_rho_on_grid: start"
#endif

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  c_rho_old = 0
  rho = 0.0_num

  ! ______________________________________
  ! DEPOSIT Charge on the grid

  SELECT CASE (c_dim)
    ! ___ In 2D _________________________________________
  CASE (2)

    CALL pxrdepose_rho_on_grid_sub_openmp_2d(rho, nx, ny, nz, nxjguards, nyjguards,   &
    nzjguards, nox, noy, noz, dx, dy, dz, dt, c_rho_old)

    ! ___ In 3D _________________________________________
  CASE DEFAULT

    ! ___ Optimized functions ______________________
    IF (rhodepo.EQ.0) THEN

      IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
        CALL pxrdepose_rho_on_grid_sub_openmp_3d(depose_rho_vecHVv4_3_3_3, rho, nx,   &
        ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,       &
        LVEC_charge_depo, c_rho_old)
      ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN
        CALL pxrdepose_rho_on_grid_sub_openmp_3d(depose_rho_vecHVv2_2_2_2, rho, nx,   &
        ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,       &
        LVEC_charge_depo, c_rho_old)
      ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN
        CALL pxrdepose_rho_on_grid_sub_openmp_3d(depose_rho_vecHVv2_1_1_1, rho, nx,   &
        ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,       &
        LVEC_charge_depo, c_rho_old)
      ELSE
        CALL pxrdepose_rho_on_grid_sub_openmp_3d_n(rho, nx, ny, nz, nxjguards,        &
        nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, c_rho_old)
      ENDIF

      ! ___ Scalar subroutines _______________________
    ELSE IF (rhodepo.EQ.1) THEN

      IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
        CALL pxrdepose_rho_on_grid_sub_openmp_3d(depose_rho_scalar_3_3_3, rho, nx,    &
        ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,       &
        LVEC_charge_depo, c_rho_old)
      ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN
        CALL pxrdepose_rho_on_grid_sub_openmp_3d(depose_rho_scalar_2_2_2, rho, nx,    &
        ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,       &
        LVEC_charge_depo, c_rho_old)
      ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN
        CALL pxrdepose_rho_on_grid_sub_openmp_3d(depose_rho_scalar_1_1_1, rho, nx,    &
        ny, nz, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,       &
        LVEC_charge_depo, c_rho_old)
      ELSE
        CALL pxrdepose_rho_on_grid_sub_openmp_3d_n(rho, nx, ny, nz, nxjguards,        &
        nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, c_rho_old)
      ENDIF

      ! ___ Non-optimized general function ____________________
    ELSE

      CALL pxrdepose_rho_on_grid_sub_openmp_3d_n(rho, nx, ny, nz, nxjguards,          &
      nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt, c_rho_old)

    ENDIF

  END SELECT
  IF (it.ge.timestat_itstart) THEN
    localtimes(12) = localtimes(12) + (MPI_WTIME() - tmptime)
  ENDIF

  ! For debugging
#if defined(DEBUG)
  WRITE(0, *) "pxrdepose_rho_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_rho_on_grid


! ________________________________________________________________________________________
!> @brief
!> Deposit rho in each tile in 3D with the subroutine pxr_depose_rho_n()
!
!> @details
!> This subroutine perform the charge deposition among the tiles using OpenMP version
!> in 3D.
!> It avoids conflict while reducing tile charge in the global charge array.
!> This subroutine uses only the general order function pxr_depose_rho_n().
!>
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
!
!> @param[inout] rhog global array for the charge
!> @param[in] nxx, nyy, nzz number of cells
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells
!> @param[in] noxx, noyy, nozz interpolation order
!> @param[in] dxx, dyy, dzz space discretization steps
!> @param[in] dtt time step
!> @param[in] c_rho_old
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d_n(rhog, nxx, nyy, nzz, nxjguard,       &
  nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, c_rho_old)
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
  ! Declarations
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz, c_rho_old
  REAL(num), INTENT(IN)    :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT):: rhog(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)             :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                     :: isdeposited=.FALSE._lp

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, rhog, noxx, &
  !$OMP noyy, nozz, aofgrid_tiles, c_dim, c_rho_old) PRIVATE(ix, iy, iz, ispecies,    &
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
        currg%arr1=0._num
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE._lp
          ENDIF
          ! Depose charge in rhotile
          CALL pxr_depose_rho_n(currg%arr1, count, curr_tile%part_x,               &
          curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,    &
          curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                       &
          curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg, nzjg,  &
          noxx, noyy, nozz, .TRUE._lp, .FALSE._lp)
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          rhog(jmin:jmax, kmin:kmax, lmin:lmax)=rhog(jmin:jmax, kmin:kmax,            &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding charge from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- X
          rhog(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = rhog(jminc:jmin-1,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,          &
          -nzjg:nzc+nzjg)
          rhog(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = rhog(jmax+1:jmaxc,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,    &
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Y
          rhog(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = rhog(jmin:jmax, kminc:kmin-1,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          rhog(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = rhog(jmin:jmax, kmax+1:kmaxc,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Z
          rhog(jmin:jmax, kmin:kmax, lminc:lmin-1) = rhog(jmin:jmax, kmin:kmax,       &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          rhog(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = rhog(jmin:jmax, kmin:kmax,       &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d_n


! ________________________________________________________________________________________
!> @brief
!> Deposit rho in each tile in 3D with the subroutine given in parameter
!
!> @details
!> This subroutine perform the charge deposition among the tiles using OpenMP version.
!> It avoids conflict while reducing tile charge in the global charge array.
!>
!> This version uses arbitrary charge deposition subroutines specified as a parameter:
!> func_order.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] func_order subroutine for the charge deposition
!> @param[in] rhog global array for the charge
!> @param[in] nxx, nyy, nzz number of cells
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells
!> @param[in] noxx, noyy, nozz interpolation order
!> @param[in] dxx, dyy, dzz space discretization steps
!> @param[in] dtt time step
!> @param[in] lvectt vector length
!> @param[in] c_rho_old
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d(func_order, rhog, nxx, nyy, nzz,       &
  nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, lvectt,           &
  c_rho_old)
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
  ! Interfaces for func_order
  INTERFACE
    SUBROUTINE func_order(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx, dy, dz,    &
      nx, ny, nz, nxguard, nyguard, nzguard, lvect) !#do not parse

      USE PICSAR_precision
      USE constants
      IMPLICIT NONE

      INTEGER(idp), INTENT (IN)    :: np, nx, ny, nz, nxguard, nyguard, nzguard
      REAL(num), INTENT(IN OUT) ::                                                    &
      rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      INTEGER(idp), INTENT (IN)    :: lvect
      REAL(num), INTENT(IN)     :: xp(np), yp(np), zp(np), w(np)
      REAL(num), INTENT(IN)     :: q, dx, dy, dz, xmin, ymin, zmin
    END SUBROUTINE
  END INTERFACE
  ! ______________________________________________________________________
  ! Declarations
  INTEGER(idp), INTENT(IN)  :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN)  :: noxx, noyy, nozz, c_rho_old, lvectt
  REAL(num), INTENT(IN)     :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: rhog(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,   &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)              :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)  :: isdeposited=.FALSE._lp

  IF (nspecies .EQ. 0_idp) RETURN

  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, rhog, noxx, &
  !$OMP noyy, nozz, aofgrid_tiles, c_dim, c_rho_old) FIRSTPRIVATE(lvectt) PRIVATE(ix, &
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
        nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0._num
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE._lp
          ENDIF
          ! Depose charge in rhotile
          CALL func_order(currg%arr1, count, curr_tile%part_x, curr_tile%part_y,   &
          curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,                      &
          curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                       &
          curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg, nzjg,  &
          lvectt)
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          rhog(jmin:jmax, kmin:kmax, lmin:lmax)= rhog(jmin:jmax, kmin:kmax,           &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding charge from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- X
          rhog(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = rhog(jminc:jmin-1,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,          &
          -nzjg:nzc+nzjg)
          rhog(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = rhog(jmax+1:jmaxc,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,    &
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Y
          rhog(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = rhog(jmin:jmax, kminc:kmin-1,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          rhog(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = rhog(jmin:jmax, kmax+1:kmaxc,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Z
          rhog(jmin:jmax, kmin:kmax, lminc:lmin-1) = rhog(jmin:jmax, kmin:kmax,       &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          rhog(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = rhog(jmin:jmax, kmin:kmax,       &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d


! ________________________________________________________________________________________
!> @brief
!> Deposit rho in each tile in 2D
!
!> @details
!> This subroutine perform the charge deposition among the tiles
!> using OpenMP version in 2D.
!> It avoids conflict while reducing tile charge in the global charge array.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
!
!> @param[inout] rhog global array for the charge
!> @param[in] nxx, nyy, nzz number of cells
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells
!> @param[in] noxx, noyy, nozz interpolation order
!> @param[in] dxx, dyy, dzz space discretization steps
!> @param[in] dtt time step
!> @param[in] c_rho_old
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_2d(rhog, nxx, nyy, nzz, nxjguard,         &
  nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, c_rho_old)
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
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz, c_rho_old
  REAL(num), INTENT(IN) :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT) :: rhog(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,   &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER :: curr_tile
  TYPE(grid_tile), POINTER :: currg
  INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)  :: isdeposited=.FALSE._lp

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, rhog, noxx, &
  !$OMP noyy, nozz, aofgrid_tiles, c_dim, c_rho_old) PRIVATE(ix, iy, iz, ispecies,    &
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
        nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0._num
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE._lp
          ENDIF
          ! Depose charge in rhotile
          SELECT CASE (c_rho_old)
          CASE(1)! Rho at older time
            CALL pxr_depose_rhoold_n_2dxz(currg%arr1(:, 0, :), count,              &
            curr_tile%part_x, curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy, &
            curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%pid(1, wpid),         &
            curr%charge, curr_tile%x_grid_tile_min, curr_tile%z_grid_tile_min, dtt,   &
            dxx, dzz, nxc, nzc, nxjg, nzjg, noxx, nozz, .TRUE._lp, .FALSE._lp)
          CASE DEFAULT! Rho at current time
            CALL pxr_depose_rho_n_2dxz(currg%arr1(:, 0, :), count,                 &
            curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1,    &
            wpid), curr%charge, curr_tile%x_grid_tile_min, curr_tile%z_grid_tile_min, &
            dxx, dzz, nxc, nzc, nxjg, nzjg, noxx, nozz, .TRUE._lp, .FALSE._lp,        &
            .FALSE._lp, 0_idp)
          END SELECT
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          rhog(jmin:jmax, kmin:kmax, lmin:lmax)= rhog(jmin:jmax, kmin:kmax,           &
          lmin:lmax)+ currg%arr1(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding charge from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- X
          rhog(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = rhog(jminc:jmin-1,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,          &
          -nzjg:nzc+nzjg)
          rhog(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = rhog(jmax+1:jmaxc,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,    &
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Y
          rhog(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = rhog(jmin:jmax, kminc:kmin-1,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          rhog(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = rhog(jmin:jmax, kmax+1:kmaxc,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Z
          rhog(jmin:jmax, kmin:kmax, lminc:lmin-1) = rhog(jmin:jmax, kmin:kmax,       &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          rhog(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = rhog(jmin:jmax, kmin:kmax,       &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_2d


! ________________________________________________________________________________________
!> @brief
!> Deposit rho in each tile in 3D with the scalar subroutine
!
!> @details
!> This subroutine perform the charge deposition among the tiles using OpenMP version in 3D.
!> It avoids conflict while reducing tile charge in the global charge array.
!> This subroutine uses only the scalar subroutines
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[inout] rhog global array for the charge
!> @param[in] nxx, nyy, nzz number of cells
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells
!> @param[in] noxx, noyy, nozz interpolation order
!> @param[in] dxx, dyy, dzz space discretization steps
!> @param[in] dtt time step
!> @param[in] c_rho_old
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d_scalar(rhog, nxx, nyy, nzz, nxjguard,  &
  nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, c_rho_old)
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
  ! Declarations
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz, c_rho_old
  REAL(num), INTENT(IN)    :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT):: rhog(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)             :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                     :: isdeposited=.FALSE._lp
  
  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, rhog, noxx, &
  !$OMP noyy, nozz, aofgrid_tiles, c_dim, c_rho_old) PRIVATE(ix, iy, iz, ispecies,    &
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
        nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0._num
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE._lp
          ENDIF
          ! Depose charge in rhotile
          IF ((noxx.eq.3).AND.(noyy.eq.3).AND.(nozz.eq.3)) THEN
            CALL depose_rho_scalar_3_3_3(currg%arr1, count, curr_tile%part_x,      &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, 0_idp)
          ELSE IF ((noxx.eq.2).AND.(noyy.eq.2).AND.(nozz.eq.2)) THEN
            CALL depose_rho_scalar_2_2_2(currg%arr1, count, curr_tile%part_x,      &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, 0_idp)
          ELSE IF ((noxx.eq.1).AND.(noyy.eq.1).AND.(nozz.eq.1)) THEN
            CALL depose_rho_scalar_1_1_1(currg%arr1, count, curr_tile%part_x,      &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, 0_idp)
          ELSE
            CALL pxr_depose_rho_n(currg%arr1, count, curr_tile%part_x,             &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, noxx, noyy, nozz, .TRUE._lp, .FALSE._lp)
          ENDIF
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          rhog(jmin:jmax, kmin:kmax, lmin:lmax)=rhog(jmin:jmax, kmin:kmax,            &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO

  !! Adding charge from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- X
          rhog(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = rhog(jminc:jmin-1,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,          &
          -nzjg:nzc+nzjg)
          rhog(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = rhog(jmax+1:jmaxc,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,    &
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Y
          rhog(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = rhog(jmin:jmax, kminc:kmin-1,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          rhog(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = rhog(jmin:jmax, kmax+1:kmaxc,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Z
          rhog(jmin:jmax, kmin:kmax, lminc:lmin-1) = rhog(jmin:jmax, kmin:kmax,       &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          rhog(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = rhog(jmin:jmax, kmin:kmax,       &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d_scalar


! ________________________________________________________________________________________
!> @brief
!> Deposit rho in each tile in 3D with the vectorized subroutine
!
!> @details
!> This subroutine perform the charge deposition among the tiles using OpenMP version in 3D.
!> It avoids conflict while reducing tile charge in the global charge array.
!> This subroutine uses only the vectorized subroutines
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[inout] rhog global array for the charge
!> @param[in] nxx, nyy, nzz number of cells
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells
!> @param[in] noxx, noyy, nozz interpolation order
!> @param[in] dxx, dyy, dzz space discretization steps
!> @param[in] dtt time step
!> @param[in] c_rho_old
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d_vecto(rhog, nxx, nyy, nzz, nxjguard,   &
  nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, c_rho_old, lvect)
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
  ! Declarations
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxjguard, nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz, c_rho_old, lvect
  REAL(num), INTENT(IN)    :: dxx, dyy, dzz, dtt
  REAL(num), INTENT(IN OUT):: rhog(-nxjguard:nxx+nxjguard, -nyjguard:nyy+nyjguard,    &
  -nzjguard:nzz+nzjguard)
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)             :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)                     :: isdeposited=.FALSE._lp

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex, ntiley, ntilez, nspecies,               &
  !$OMP species_parray, nxjguard, nyjguard, nzjguard, dxx, dyy, dzz, dtt, rhog, noxx, &
  !$OMP noyy, nozz, aofgrid_tiles, c_dim, c_rho_old, lvect) PRIVATE(ix, iy, iz,       &
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
        nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        currg=>aofgrid_tiles(ix, iy, iz)
        currg%arr1=0._num
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE._lp
          ENDIF
          ! Depose charge in rhotile
          IF ((noxx.eq.3).AND.(noyy.eq.3).AND.(nozz.eq.3)) THEN
            CALL depose_rho_vecHVv4_3_3_3(currg%arr1, count, curr_tile%part_x,     &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, LVECT)
          ELSE IF ((noxx.eq.2).AND.(noyy.eq.2).AND.(nozz.eq.2)) THEN
            CALL depose_rho_vecHVv2_2_2_2(currg%arr1, count, curr_tile%part_x,     &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, lvect)
          ELSE IF ((noxx.eq.1).AND.(noyy.eq.1).AND.(nozz.eq.1)) THEN
            CALL depose_rho_vecHVv2_1_1_1(currg%arr1, count, curr_tile%part_x,     &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, lvect)
          ELSE
            CALL pxr_depose_rho_n(currg%arr1, count, curr_tile%part_x,             &
            curr_tile%part_y, curr_tile%part_z, curr_tile%pid(1, wpid), curr%charge,  &
            curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,                     &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, nxc, nyc, nzc, nxjg, nyjg,      &
            nzjg, noxx, noyy, nozz, .TRUE._lp, .FALSE._lp)
          ENDIF
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          rhog(jmin:jmax, kmin:kmax, lmin:lmax)=rhog(jmin:jmax, kmin:kmax,            &
          lmin:lmax)+currg%arr1(0:nxc, 0:nyc, 0:nzc)
        ENDIF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO

  !! Adding charge from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- X
          rhog(jminc:jmin-1, kminc:kmaxc, lminc:lmaxc) = rhog(jminc:jmin-1,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(-nxjg:-1, -nyjg:nyc+nyjg,          &
          -nzjg:nzc+nzjg)
          rhog(jmax+1:jmaxc, kminc:kmaxc, lminc:lmaxc) = rhog(jmax+1:jmaxc,           &
          kminc:kmaxc, lminc:lmaxc)+ currg%arr1(nxc+1:nxc+nxjg, -nyjg:nyc+nyjg,    &
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Y
          rhog(jmin:jmax, kminc:kmin-1, lminc:lmaxc) = rhog(jmin:jmax, kminc:kmin-1,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, -nyjg:-1, -nzjg:nzc+nzjg)
          rhog(jmin:jmax, kmax+1:kmaxc, lminc:lmaxc) = rhog(jmin:jmax, kmax+1:kmaxc,  &
          lminc:lmaxc)+ currg%arr1(0:nxc, nyc+1:nyc+nyjg, -nzjg:nzc+nzjg)
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
        isdeposited=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isdeposited=.TRUE._lp
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
          ! --- RHO
          ! - FACES +/- Z
          rhog(jmin:jmax, kmin:kmax, lminc:lmin-1) = rhog(jmin:jmax, kmin:kmax,       &
          lminc:lmin-1)+ currg%arr1(0:nxc, 0:nyc, -nzjg:-1)
          rhog(jmin:jmax, kmin:kmax, lmax+1:lmaxc) = rhog(jmin:jmax, kmin:kmax,       &
          lmax+1:lmaxc)+ currg%arr1(0:nxc, 0:nyc, nzc+1:nzc+nzjg)
        END IF
      END DO
    END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrdepose_rho_on_grid_sub_openmp_3d_vecto

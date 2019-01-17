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
! FIELD_GATHERING_MANAGER.F90
!
! This file contains subroutines to manage the field gathering in 3D.
!
! Developers:
! - Henri vincenti
! - Mathieu Lobet
!
! options:
! - DEV: activates developer's secret subroutines
! - DEBUG: activates DEBUG prints and outputs
!
! List of subroutines:
!
! - field_gathering
! - field_gathering_sub
! - geteb3d_energy_conserving
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Field gathering main subroutine in 3D called in the main loop when not coupled
!> with the particle pusher.
! ________________________________________________________________________________________
SUBROUTINE field_gathering
  USE fields, ONLY: bx, by, bz, ex, ey, ez, l_lower_order_in_v, nox, noy, noz,       &
    nxguards, nxjguards, nyguards, nyjguards, nzguards, nzjguards
  USE params, ONLY: dt
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz
  IMPLICIT NONE

#if defined(DEBUG)
  WRITE(0, *) "Field gathering: start"
#endif

  CALL field_gathering_sub(ex, ey, ez, bx, by, bz, nx, ny, nz, nxguards, nyguards,    &
  nzguards, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,           &
  l_lower_order_in_v)


#if defined(DEBUG)
  WRITE(0, *) "Field gathering: stop"
#endif

END SUBROUTINE field_gathering


! ________________________________________________________________________________________
!> @brief
!> This subroutine performs the field gathering in 3D only
! ________________________________________________________________________________________
SUBROUTINE field_gathering_sub(exg, eyg, ezg, bxg, byg, bzg, nxx, nyy, nzz, nxguard,  &
  nyguard, nzguard, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, &
  l_lower_order_in_v_in)
  USE grid_tilemodule, ONLY: aofgrid_tiles
  USE mpi
  USE particle_properties, ONLY: nspecies
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  USE time_stat, ONLY: localtimes, timestat_itstart
  ! Vtune/SDE profiling
  IMPLICIT NONE

  ! ___ Parameter declaration ________________________________________
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard,     &
  nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  LOGICAL(lp)              :: l_lower_order_in_v_in
  REAL(num), INTENT(IN)    :: exg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,         &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: eyg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,         &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: ezg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,         &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bxg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,         &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: byg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,         &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bzg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,         &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: dxx, dyy, dzz, dtt
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                :: tdeb, tend
  INTEGER(idp)             :: nxc, nyc, nzc
  INTEGER(idp)             :: nxjg, nyjg, nzjg
  INTEGER(idp)             :: nxt, nyt, nzt
  INTEGER(idp)             :: nxt_o, nyt_o, nzt_o
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile, bytile, bztile
  LOGICAL(lp)                   :: isgathered=.FALSE._lp

  IF (nspecies .EQ. 0_idp) RETURN

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==3
  CALL start_vtune_collection()
#endif
#if SDE==3
  CALL start_sde_collection()
#endif

  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex,                                         &
  !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, nxjguard, nyjguard,  &
  !$OMP nzjguard, nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, dxx, dyy,  &
  !$OMP dzz, dtt, noxx, noyy, nozz, c_dim, l_lower_order_in_v_in, fieldgathe,         &
  !$OMP LVEC_fieldgathe) PRIVATE(ix, iy, iz, ispecies, curr, curr_tile, count,        &
  !$OMP extile, eytile, eztile, bxtile, bytile, bztile, nxt, nyt, nzt,                &
  !$OMP jmin, jmax, kmin, kmax, lmin, lmax, nxc, nyc, nzc, nxjg, nyjg, nzjg,          &
  !$OMP isgathered, nxt_o, nyt_o, nzt_o)      
  nxt_o=0_idp
  nyt_o=0_idp
  nzt_o=0_idp   
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)               
  DO iz=1, ntilez! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr=>species_parray(1)
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
        nxc=curr_tile%nx_cells_tile
        nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        isgathered=.FALSE._lp
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr=>species_parray(ispecies)
          IF (curr%is_antenna) CYCLE
          IF (curr%lfreeze) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isgathered=.TRUE.
        END DO
        IF (isgathered) THEN
          nxt=jmax-jmin+1_idp
          nyt=kmax-kmin+1_idp
          nzt=lmax-lmin+1_idp
          ! - Only resize temporary private grid tile arrays if tile size has changed
          ! - i.e (nxt!=nxt_o or nyt!= nyt_o or nzt!=nzt_o)
          ! - If tile array not allocated yet, allocate tile array with sizes nxt, nyt,
          ! - nzt
          IF (.NOT. ALLOCATED(extile)) THEN
            ALLOCATE(extile(nxt,nyt,nzt))
            ALLOCATE(eytile(nxt,nyt,nzt))
            ALLOCATE(eztile(nxt,nyt,nzt))
            ALLOCATE(bxtile(nxt,nyt,nzt))
            ALLOCATE(bytile(nxt,nyt,nzt))
            ALLOCATE(bztile(nxt,nyt,nzt))
          ELSE
            IF ((nxt .NE. nxt_o) .OR. (nyt .NE. nyt_o) .OR. (nzt .NE. nzt_o)) THEN
              DEALLOCATE(extile,eytile,eztile,bxtile,bytile,bztile)
              ALLOCATE(extile(nxt,nyt,nzt))
              ALLOCATE(eytile(nxt,nyt,nzt))
              ALLOCATE(eztile(nxt,nyt,nzt))
              ALLOCATE(bxtile(nxt,nyt,nzt))
              ALLOCATE(bytile(nxt,nyt,nzt))
              ALLOCATE(bztile(nxt,nyt,nzt))
            ENDIF
          ENDIF
          nxt_o=nxt
          nyt_o=nyt
          nzt_o=nzt
          ! - Copy values of field arrays in temporary grid tile arrays 
          extile=exg(jmin:jmax, kmin:kmax, lmin:lmax)
          eytile=eyg(jmin:jmax, kmin:kmax, lmin:lmax)
          eztile=ezg(jmin:jmax, kmin:kmax, lmin:lmax)
          bxtile=bxg(jmin:jmax, kmin:kmax, lmin:lmax)
          bytile=byg(jmin:jmax, kmin:kmax, lmin:lmax)
          bztile=bzg(jmin:jmax, kmin:kmax, lmin:lmax)
          DO ispecies=1, nspecies! LOOP ON SPECIES
            ! - Get current tile properties
            ! - Init current tile variables

            curr=>species_parray(ispecies)
            IF (curr%is_antenna) CYCLE
            IF (curr%lfreeze) CYCLE
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            IF (count .EQ. 0) CYCLE
            curr_tile%part_ex(1:count) = 0.0_num
            curr_tile%part_ey(1:count) = 0.0_num
            curr_tile%part_ez(1:count) = 0.0_num
            curr_tile%part_bx(1:count)=0.0_num
            curr_tile%part_by(1:count)=0.0_num
            curr_tile%part_bz(1:count)=0.0_num
            !!! ---- Loop by blocks over particles in a tile (blocking)
            !!! --- Gather electric field on particles

            !!! --- Gather electric and magnetic fields on particles
            CALL geteb3d_energy_conserving(count, curr_tile%part_x, curr_tile%part_y, &
            curr_tile%part_z, curr_tile%part_ex, curr_tile%part_ey,                   &
            curr_tile%part_ez, curr_tile%part_bx, curr_tile%part_by,                  &
            curr_tile%part_bz, curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,  &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, curr_tile%nx_cells_tile,        &
            curr_tile%ny_cells_tile, curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, noxx, &
            noyy, nozz, extile, eytile, eztile, bxtile, bytile, bztile, .FALSE._lp,   &
            l_lower_order_in_v_in, LVEC_fieldgathe, fieldgathe)
          END DO! END LOOP ON SPECIES
        ENDIF
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END DO 
  IF (ALLOCATED(extile)) THEN ! Deallocation of tile arrays 
    DEALLOCATE(extile,eytile,eztile,bxtile,bytile,bztile)
  ENDIF
  !$OMP END PARALLEL
  
#if VTUNE==3
  CALL stop_vtune_collection()
#endif
#if SDE==3
  CALL stop_sde_collection()
#endif

  IF (it.ge.timestat_itstart) THEN
    tend=MPI_WTIME()
    localtimes(14) = localtimes(14) + (tend-tdeb)
  ENDIF

END SUBROUTINE field_gathering_sub


! ________________________________________________________________________________________
!> @brief
!> General subroutines for the 3D field gathering
!
!> @details
!> This subroutine controls the different algorithms for the field gathering
!> Choice of an algorithm is done using the argument field_gathe_algo.
!> This subroutine is called in the subroutine field_gathering_sub().
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2015-2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle positions
!> @param[out] ex, ey, ez particle electric field
!> @param[out] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile origin
!> @param[in] dx, dy, dz space discretization
!> @param[in] nx, ny, nz number of cells in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] nox, noy, noz shape factor order
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] ll4symtry
!> @param[in] l_lower_order_in_v
!> @param[in] field_gathe_algo gathering algorithm
!> @param[in] lvect vector length
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving(np, xp, yp, zp, ex, ey, ez, bx, by, bz, xmin,    &
  ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox, noy, noz, exg,    &
  eyg, ezg, bxg, byg, bzg, ll4symtry, l_lower_order_in_v, lvect, field_gathe_algo)
  USE picsar_precision, ONLY: idp, lp, num
  USE fields, ONLY: l_nodalgrid
  implicit none

  integer(idp)                  :: field_gathe_algo
  integer(idp)                  :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard,   &
  nzguard
  LOGICAL(lp), intent(in)      :: ll4symtry, l_lower_order_in_v
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  real(num)                     :: xmin, ymin, zmin, dx, dy, dz
  integer(idp)                  :: lvect
  real(num), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  real(num), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg

  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(3), nvalid(3)
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)

  call geteb3d_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by, bz,      &
  xmin, ymin, zmin, dx, dy, dz, nox, noy, noz, exg, nguard, nvalid, eyg, nguard,      &
  nvalid, ezg, nguard, nvalid, bxg, nguard, nvalid, byg, nguard, nvalid, bzg, nguard, &
  nvalid, ll4symtry, l_lower_order_in_v, l_nodalgrid, lvect, field_gathe_algo)
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> General subroutines for the 3D field gathering, adapted for field
!> arrays having different sizes depending on their nodal/cell-centered nature
!>
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle positions
!> @param[out] ex, ey, ez particle electric field
!> @param[out] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile origin
!> @param[in] dx, dy, dz space discretization
!> @param[in] nox, noy, noz shape factor order
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!> exg, eyg, ezg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 3 integers)
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of
!> the bxg, byg, bzg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 3 integers)
!> @param[in] ll4symtry
!> @param[in] l_lower_order_in_v
!> @param[in] field_gathe_algo gathering algorithm
!> @param[in] lvect vector length
!>
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by, bz,  &
  xmin, ymin, zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg,        &
  eyg_nguard, eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid,     &
  byg, byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry,                  &
  l_lower_order_in_v, l_nodal, lvect, field_gathe_algo)            !#do not wrap
  USE params, ONLY: fieldgathe, lvec_fieldgathe
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp)                  :: field_gathe_algo
  integer(idp)                  :: np, nox, noy, noz
  integer(idp), intent(IN)      :: exg_nguard(3), exg_nvalid(3), eyg_nguard(3),       &
  eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3), bxg_nguard(3), bxg_nvalid(3),          &
  byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  LOGICAL(lp), intent(in)      :: ll4symtry, l_lower_order_in_v, l_nodal
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  real(num)                     :: xmin, ymin, zmin, dx, dy, dz
  integer(idp)                  :: lvect
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1,                                       &
  -exg_nguard(3):exg_nvalid(3)+exg_nguard(3)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1,                                       &
  -eyg_nguard(3):eyg_nvalid(3)+eyg_nguard(3)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1,                                       &
  -ezg_nguard(3):ezg_nvalid(3)+ezg_nguard(3)-1)
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1,                                       &
  -bxg_nguard(3):bxg_nvalid(3)+bxg_nguard(3)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1,                                       &
  -byg_nguard(3):byg_nvalid(3)+byg_nguard(3)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1,                                       &
  -bzg_nguard(3):bzg_nvalid(3)+bzg_nguard(3)-1)

  IF (np .EQ. 0_idp) RETURN

#if defined(DEV)
  ! Maintain variables nx, ny, nz, nxguard, nyguard, nzguard
  ! for compatibility with the DEV routines
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard
  nx = exg_nvalid(1)-1
  ny = exg_nvalid(2)-1
  nz = exg_nvalid(3)-1
  nxguard = exg_nguard(1)
  nyguard = exg_nguard(2)
  nzguard = exg_nguard(3)
#endif

  SELECT CASE(field_gathe_algo)

    ! ______________________________________________________________________________
    ! Developer's functions (experimental or under development)
#if defined(DEV)

    ! ______________________________________
    ! Experimental vectorized field gathering subroutines
  CASE(7)

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

      CALL geteb3d_energy_conserving_vecV3_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,   &
      eyg, ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

      CALL geteb3d_energy_conserving_vecV3_2_2_2(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,   &
      eyg, ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      CALL geteb3d_energy_conserving_vecV2_3_3_3(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,   &
      eyg, ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)

      ! Arbitrary order
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF


    ! ______________________________________
    ! Vectorized non efficient field gathering subroutines
  CASE(6)

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      !!! --- Gather electric fields on particles
      CALL geteb3d_energy_conserving_vecV1_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,   &
      eyg, ezg, bxg, byg, bzg, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      !!! --- Gather electric fields on particles
      CALL geteb3d_energy_conserving_vecV1_2_2_2(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,   &
      eyg, ezg, bxg, byg, bzg, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      !!! --- Gather electric fields on particles
      CALL geteb3d_energy_conserving_vec_3_3_3(np, xp, yp, zp, ex, ey, ez, bx, by,    &
      bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,   &
      eyg, ezg, bxg, byg, bzg, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF

    ! ______________________________________
    ! Vectorized field gathering subroutines, separated E and B functions
  CASE(5)

    IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      !!! --- Gather electric fields on particles
      CALL gete3d_energy_conserving_vec2_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin,      &
      ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg, ezg,   &
      LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL getb3d_energy_conserving_vec2_3_3_3(np, xp, yp, zp, bx, by, bz, xmin,      &
      ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, bxg, byg, bzg,   &
      LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF

    ! ______________________________________
    ! Vectorized field gathering subroutine by block splited into smaller loops
  CASE(4)

    IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      !!! --- Gather electric and magnetic fields on particles
      CALL geteb3d_energy_conserving_blockvec2_3_3_3(np, xp, yp, zp, ex, ey, ez, bx,  &
      by, bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard,    &
      exg, eyg, ezg, bxg, byg, bzg, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF

    ! ______________________________________
    ! Linearized field gathering subroutines
  CASE(3)

    IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      !!! --- Gather electric field on particles
      CALL gete3d_energy_conserving_linear_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin,    &
      ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg, ezg,   &
      l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL getb3d_energy_conserving_linear_3_3_3(np, xp, yp, zp, bx, by, bz, xmin,    &
      ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, bxg, byg, bzg,   &
      l_lower_order_in_v, l_nodal)
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF

#endif
    ! ______________________________________________________________________________

    ! ______________________________________________
    ! Arbitrary order, non-optimized subroutines
  CASE(2)

    !!! --- Gather electric field on particles
    CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin, zmin, &
    dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,          &
    eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin, zmin, &
    dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,          &
    byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)

    ! ______________________________________________
    ! Non-optimized scalar subroutines

  CASE(1)

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      !!! --- Gather electric field on particles
      CALL gete3d_energy_conserving_scalar_1_1_1(np, xp, yp, zp, ex, ey, ez, xmin,    &
      ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,           &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL getb3d_energy_conserving_scalar_1_1_1(np, xp, yp, zp, bx, by, bz, xmin,    &
      ymin, zmin, dx, dy, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,           &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      !!! --- Gather electric field on particles
      CALL gete3d_energy_conserving_scalar_2_2_2(np, xp, yp, zp, ex, ey, ez, xmin,    &
      ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,           &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL getb3d_energy_conserving_scalar_2_2_2(np, xp, yp, zp, bx, by, bz, xmin,    &
      ymin, zmin, dx, dy, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,           &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      !!! --- Gather electric field on particles
      CALL gete3d_energy_conserving_scalar_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin,    &
      ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,           &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL getb3d_energy_conserving_scalar_3_3_3(np, xp, yp, zp, bx, by, bz, xmin,    &
      ymin, zmin, dx, dy, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,           &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF


    ! ________________________________________
    ! Optimized subroutines, E and B in the same vectorized loop, default
  CASE DEFAULT

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      !!! --- Gather electric and magnetic fields on particles
      CALL geteb3d_energy_conserving_vecV4_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,      &
      byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      !!! --- Gather electric and magnetic fields on particles
      CALL geteb3d_energy_conserving_vecV4_2_2_2(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,      &
      byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN
      !!! --- Gather electric and magnetic fields on particles
      CALL geteb3d_energy_conserving_vecV3_3_3_3(np, xp, yp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,      &
      byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v, l_nodal)
    ELSE
      !!! --- Gather electric field on particles
      CALL pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,  &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
      zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,  &
      byg_nvalid, bzg, bzg_nguard, bzg_nvalid, ll4symtry, l_lower_order_in_v, l_nodal)
    ENDIF
  END SELECT
END SUBROUTINE

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
! PARTICLE_PUSHER_MANAGER_2D.F90
!
! Subroutines for managing the particle pushers in 2d.
!
! Developers:
! Henri Vincenti
! Mathieu Lobet
!
! Date:
! Creation 2015
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> This subroutine is the 2D version of the particle pusher + field gathering
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] nxx, nyy, nzz number of cells in each direction for the grids,
!> nyy should be equal to 1.
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> for the grids
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells for the current grids
!> @param[in] noxx, noyy, nozz interpolation orders
!> @param[in] dxx, dyy, dzz space steps
!> @param[in] dtt time step
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher_sub_2d(exg, eyg, ezg, bxg, byg, bzg,  &
  nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard, nyjguard, nzjguard, noxx, noyy, &
  nozz, dxx, dyy, dzz, dtt)
  USE grid_tilemodule, ONLY: aofgrid_tiles
  USE mpi
  USE output_data, ONLY: pushtime
  USE particle_properties, ONLY: bxoldpid, byoldpid, bzoldpid, exoldpid, eyoldpid,   &
    ezoldpid, nspecies, particle_pusher
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  USE time_stat, ONLY: localtimes
  ! Vtune/SDE profiling

  IMPLICIT NONE
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard,     &
  nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  REAL(num), INTENT(IN) :: exg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,            &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: eyg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,            &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: ezg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,            &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: bxg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,            &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: byg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,            &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: bzg(-nxguard:nxx+nxguard, -nyguard:nyy+nyguard,            &
  -nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: dxx, dyy, dzz, dtt
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER :: curr_tile
  REAL(num) :: tdeb, tend
  INTEGER(idp) :: nxc, nyc, nzc, ipmin, ipmax, ip
  INTEGER(idp) :: nxjg, nzjg
  INTEGER(idp)             :: nxt, nyt, nzt
  INTEGER(idp)             :: nxt_o, nyt_o, nzt_o
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile, bytile, bztile
  LOGICAL(lp)  :: isgathered=.FALSE.

  tdeb=MPI_WTIME()
#if PROFILING==3
  CALL start_collection()
#endif
  !$OMP PARALLEL DEFAULT(NONE) SHARED(particle_pusher, ntilex,                        &
  !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, nxjguard, nyjguard,  &
  !$OMP nzjguard, nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, dxx, dyy,  &
  !$OMP dzz, dtt, noxx, noyy, nozz, c_dim, fieldgathe, LVEC_fieldgathe, exoldpid,     &
  !$OMP eyoldpid, ezoldpid, bxoldpid, byoldpid, bzoldpid) PRIVATE(ix,                 &
  !$OMP iy, iz, ispecies, curr, curr_tile, count, jmin, jmax, kmin, kmax,             &
  !$OMP lmin, lmax, nxc, nyc, nzc, ipmin, ipmax, ip, nxjg, nzjg, isgathered,          &
  !$OMP extile, eytile, eztile, bxtile, bytile, bztile, nxt, nyt, nzt, nxt_o, nyt_o,  &
  !$OMP nzt_o)
  nxt_o=0_idp
  nyt_o=0_idp
  nzt_o=0_idp
  !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
  DO iz=1, ntilez! LOOP ON TILES
    DO ix=1, ntilex
      curr=>species_parray(1)
      curr_tile=>curr%array_of_tiles(ix, 1, iz)
      nxjg=curr_tile%nxg_tile
      nzjg=curr_tile%nzg_tile
      jmin=curr_tile%nx_tile_min-nxjg
      jmax=curr_tile%nx_tile_max+nxjg
      kmin=0! 2D Case
      kmax=0! 2D case
      lmin=curr_tile%nz_tile_min-nzjg
      lmax=curr_tile%nz_tile_max+nzjg
      nxc=curr_tile%nx_cells_tile
      nzc=curr_tile%nz_cells_tile
      isgathered=.FALSE.
      DO ispecies=1, nspecies! LOOP ON SPECIES
        curr=>species_parray(ispecies)
        IF (curr%is_antenna) CYCLE
        IF (curr%lfreeze) CYCLE
        curr_tile=>curr%array_of_tiles(ix, 1, iz)
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
          curr_tile=>curr%array_of_tiles(ix, 1, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) CYCLE
          IF (fieldgathe.gt.-1) then
            curr_tile%part_ex(1:count) = 0.0_num
            curr_tile%part_ey(1:count) = 0.0_num
            curr_tile%part_ez(1:count) = 0.0_num
            curr_tile%part_bx(1:count)=0.0_num
            curr_tile%part_by(1:count)=0.0_num
            curr_tile%part_bz(1:count)=0.0_num
            !!! ---- Loop by blocks over particles in a tile (blocking)
            !!! --- Gather electric field on particles
            !!! --- Gather electric and magnetic fields on particles
            CALL geteb2dxz_energy_conserving(count, curr_tile%part_x,                 &
            curr_tile%part_y, curr_tile%part_z, curr_tile%part_ex, curr_tile%part_ey, &
            curr_tile%part_ez, curr_tile%part_bx, curr_tile%part_by,                  &
            curr_tile%part_bz, curr_tile%x_grid_tile_min, curr_tile%y_grid_tile_min,  &
            curr_tile%z_grid_tile_min, dxx, dyy, dzz, curr_tile%nx_cells_tile,        &
            curr_tile%ny_cells_tile, curr_tile%nz_cells_tile, nxjg, 0_idp, nzjg, noxx,&
            noyy, nozz, extile, eytile, eztile, bxtile,                               &
            bytile, bztile , .FALSE., .TRUE., LVEC_fieldgathe,                        &
            fieldgathe)
          END IF

	  SELECT CASE (particle_pusher)
	  !! Vay pusher -- Full push
	  CASE (1_idp)
	    CALL pxr_ebcancelpush3d(count, curr_tile%part_ux, curr_tile%part_uy,     &
	    curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,             &
	    curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                 &
	    curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt,       &
	    0_idp)
	  !! Boris pusher with RR (S09 model, according to VRANIC2016, https://doi.org/10.1016/j.cpc.2016.04.002) -- Full push
	  CASE (2_idp)
	    CALL pxr_boris_push_rr_S09_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,&
	    curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,                &
	    curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                    &
	    curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)
	  !! Boris pusher with RR (B08 model, according to VRANIC2016, https://doi.org/10.1016/j.cpc.2016.04.002) -- Full push
	  CASE (3_idp)
	    CALL pxr_boris_push_rr_B08_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,&
	    curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,                &
	    curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                    &
	    curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)
	  !! Boris pusher with RR (B08 model, according to VRANIC2016, https://doi.org/10.1016/j.cpc.2016.04.002) -- Full push
	  CASE (4_idp)

	    CALL pxr_boris_push_rr_LL_u_3d(count, curr_tile%part_ux, 		    &
	    curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv, 	    &
	    curr_tile%pid(1:count,exoldpid), curr_tile%pid(1:count,eyoldpid),	    &
	    curr_tile%pid(1:count,ezoldpid), curr_tile%pid(1:count,bxoldpid),	    &
	    curr_tile%pid(1:count,byoldpid), curr_tile%pid(1:count,bzoldpid), 	    &
	    curr_tile%part_ex, curr_tile%part_ey, curr_tile%part_ez, 		    &
	    curr_tile%part_bx, curr_tile%part_by, curr_tile%part_bz, 		    &
	    curr%charge, curr%mass, dtt)

	    curr_tile%pid(1:count,exoldpid) = curr_tile%part_ex
	    curr_tile%pid(1:count,eyoldpid) = curr_tile%part_ey
	    curr_tile%pid(1:count,ezoldpid) = curr_tile%part_ez
	    curr_tile%pid(1:count,bxoldpid) = curr_tile%part_bx
	    curr_tile%pid(1:count,byoldpid) = curr_tile%part_by
	    curr_tile%pid(1:count,bzoldpid) = curr_tile%part_bz
	  !! Boris pusher -- Full push

	  CASE DEFAULT
	    !! Push momentum using the Boris method in a single subroutine
	    CALL pxr_boris_push_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,   &
	    curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
	    curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
	    curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)
	  END SELECT

          CALL pxr_push2dxz(count, curr_tile%part_x, curr_tile%part_z,                &
          curr_tile%part_ux, curr_tile%part_uy, curr_tile%part_uz,                    &
          curr_tile%part_gaminv, dtt)
        END DO! END LOOP ON SPECIES
      ENDIF
    END DO
  END DO! END LOOP ON TILES
  !$OMP END DO
  IF (ALLOCATED(extile)) THEN ! Deallocation of tile arrays
    DEALLOCATE(extile,eytile,eztile,bxtile,bytile,bztile)
  ENDIF
  !$OMP END PARALLEL

#if PROFILING==3
  CALL stop_collection()
#endif

  tend=MPI_WTIME()
  pushtime=pushtime+(tend-tdeb)
  localtimes(1) = localtimes(1) + (tend-tdeb)
END SUBROUTINE field_gathering_plus_particle_pusher_sub_2d

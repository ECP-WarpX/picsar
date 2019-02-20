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
! PARTICLE_PUSHER_MANAGER.F90
!
! Subroutines for managing the particle pushers.
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Main subroutine for the field subroutine + particle pusher called
!> in the main loop (in submain.F90)
!
!> @details
!> This routine calls the subroutines for the different
!> particle pusher + field gathering algorithms:
!> * field_gathering_plus_particle_pusher_sub_2d()
!> * field_gathering_plus_particle_pusher_sub()
!> * field_gathering_plus_particle_pusher_cacheblock_sub()
!> * field_gathering_sub()
!>
!> Choice of the method is done via the parameter fg_p_pp_separated
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2015
!> Revision 10.06.2015
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher
  USE fields, ONLY: bx_p, by_p, bz_p, ex_p, ey_p, ez_p, l_lower_order_in_v, nox,     &
    noy, noz, nxguards, nxjguards, nyguards, nyjguards, nzguards, nzjguards
  USE mpi
  USE params, ONLY: dt, fg_p_pp_separated
  USE particle_properties, ONLY: nspecies, particle_pusher
  USE picsar_precision, ONLY: idp
  USE shared_data, ONLY: c_dim, dx, dy, dz, nx, ny, nz
  IMPLICIT NONE

#if defined(DEBUG)
  WRITE(0, *) "Field gathering + Push_particles: start"
#endif
  IF (nspecies .EQ. 0_idp) RETURN
  SELECT CASE (c_dim)

    ! ___________________________________________________________
    ! 2D CASE X Z
  CASE (2)

    ! Particle advance (one time step)
    CALL field_gathering_plus_particle_pusher_sub_2d(ex_p, ey_p, ez_p, bx_p, by_p, bz_p, &
    nx, ny, nz, nxguards, nyguards, nzguards, nxjguards, nyjguards, nzjguards,           &
    nox, noy, noz, dx, dy, dz, dt)

    ! ___________________________________________________________
    ! 3D CASE
  CASE DEFAULT

    ! The field gathering and the particle pusher are performed together
    IF (fg_p_pp_separated.eq.0) THEN

      CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex_p, ey_p, ez_p, &
      bx_p, by_p, bz_p, nx, ny, nz, nxguards, nyguards, nzguards,                &
      nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz, dt,            &
      l_lower_order_in_v)

    ELSE IF (fg_p_pp_separated.eq.1) THEN

      CALL field_gathering_plus_particle_pusher_sub(ex_p, ey_p, ez_p, bx_p, by_p, bz_p, &
      nx, ny,   nz, nxguards, nyguards, nzguards, nxjguards, nyjguards, nzjguards,      &
      nox, noy,    noz, dx, dy, dz, dt, l_lower_order_in_v)

      ! The field gathering and the particle pusher are performed separately
    ELSE

      CALL field_gathering_sub(ex_p, ey_p, ez_p, bx_p, by_p, bz_p, nx, ny, nz, nxguards, &
      nyguards, nzguards, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz,    &
      dt, l_lower_order_in_v)

      CALL particle_pusher_sub(ex_p, ey_p, ez_p, bx_p, by_p, bz_p, nx, ny, nz, nxguards, &
      nyguards, nzguards, nxjguards, nyjguards, nzjguards, nox, noy, noz, dx, dy, dz,    &
      dt, l_lower_order_in_v)

    ENDIF

  END SELECT

#if defined(DEBUG)
  WRITE(0, *) "Field gathering + Push_particles: stop"
#endif

END SUBROUTINE field_gathering_plus_particle_pusher


! ________________________________________________________________________________________
!> @brief
!> Particle pusher in 3D called by the main function push_particle
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
!> Revision 10.06.2015
!
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] nxx, nyy, nzz number of cells in each direction for the grids
!> @param[in] nxguard, nyguard, nzguard number of guard cells in
!> each direction for the grids
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells for the current grids
!> @param[in] noxx, noyy, nozz interpolation orders
!> @param[in] dxx, dyy, dzz space steps
!> @param[in] dtt time step
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!>
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher_sub(exg, eyg, ezg, bxg, byg, bzg,     &
  nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard, nyjguard, nzjguard, noxx, noyy, &
  nozz, dxx, dyy, dzz, dtt, l_lower_order_in_v_in)
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
  USE time_stat, ONLY: localtimes, timestat_itstart
  ! Vtune/SDE profiling
  IMPLICIT NONE

  ! ___ Parameter declaration __________________________________________

  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard,     &
  nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  LOGICAL(lp), INTENT(IN)  :: l_lower_order_in_v_in
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

  ! Local parameters
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                :: tdeb, tend
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin, ipmax, ip
  INTEGER(idp)             :: nxjg, nyjg, nzjg
  INTEGER(idp)             :: nxt, nyt, nzt
  INTEGER(idp)             :: nxt_o, nyt_o, nzt_o
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile, bytile, bztile
  LOGICAL(lp)              :: isgathered=.FALSE.

  tdeb=MPI_WTIME()

#if PROFILING==3
  CALL start_collection()
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex,                                         &
  !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, nxjguard, nyjguard,  &
  !$OMP nzjguard, nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, dxx, dyy,  &
  !$OMP dzz, dtt, noxx, noyy, nozz, c_dim, l_lower_order_in_v_in, particle_pusher,    &
  !$OMP fieldgathe, LVEC_fieldgathe, exoldpid, eyoldpid, ezoldpid, bxoldpid, 	      &
  !$OMP byoldpid, bzoldpid) PRIVATE(ix, iy, iz, ispecies, curr, curr_tile,   	      &
  !$OMP count, jmin, jmax, kmin, kmax, lmin, lmax, nxc, nyc, nzc, ipmin,              &
  !$OMP extile, eytile, eztile, bxtile, bytile, bztile,                               &
  !$OMP nxt, nyt, nzt, ipmax, ip, nxjg, nyjg, nzjg, isgathered, nxt_o, nyt_o, nzt_o)
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
        isgathered=.FALSE.
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
            SELECT CASE (c_dim)
            CASE (2)! 2D CASE X Z
              !!! --- Gather electric and magnetic fields on particles
              CALL geteb2dxz_energy_conserving(count, curr_tile%part_x,               &
              curr_tile%part_y, curr_tile%part_z, curr_tile%part_ex,                  &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr_tile%x_grid_tile_min,        &
              curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,    &
              curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                       &
              curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, noxx, noyy, nozz,            &
              extile, eytile, eztile, bxtile, bytile,                                 &
              bztile, .FALSE., l_lower_order_in_v_in, LVEC_fieldgathe,                &
              fieldgathe)

            CASE DEFAULT! 3D CASE

              !!! --- Gather electric and magnetic fields on particles
              CALL geteb3d_energy_conserving(count, curr_tile%part_x,                 &
              curr_tile%part_y, curr_tile%part_z, curr_tile%part_ex,                  &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr_tile%x_grid_tile_min,        &
              curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,    &
              curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                       &
              curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, noxx, noyy, nozz,            &
              extile, eytile, eztile, bxtile, bytile,                                 &
              bztile , .FALSE., l_lower_order_in_v_in, lvec_fieldgathe,               &
              fieldgathe)
            END SELECT

            SELECT CASE (particle_pusher)
            !! Vay pusher -- Full push
            CASE (1_idp)
              CALL pxr_ebcancelpush3d(count, curr_tile%part_ux, curr_tile%part_uy,    &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt,      &
              0_idp)

            !! Boris pusher with RR (S09 model, according to VRANIC2016,
            !! https://doi.org/10.1016/j.cpc.2016.04.002)-- Full push
	        CASE (2_idp)
              CALL pxr_boris_push_rr_S09_u_3d(count, curr_tile%part_ux,               &
              curr_tile%part_uy,                                                      &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)

	        !! Boris pusher with RR (B08 model, according to VRANIC2016,
	        !! https://doi.org/10.1016/j.cpc.2016.04.002)-- Full push
	        CASE (3_idp)
              CALL pxr_boris_push_rr_B08_u_3d(count, curr_tile%part_ux,               &
              curr_tile%part_uy,                                                      &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)

            !! Boris pusher with RR (LL model, according to VRANIC2016,
            !! https://doi.org/10.1016/j.cpc.2016.04.002)-- Full push
	        CASE (4_idp)
              CALL pxr_boris_push_rr_LL_u_3d(count, curr_tile%part_ux, 		      &
              curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv, 	      &
              curr_tile%pid(1:count,exoldpid), curr_tile%pid(1:count,eyoldpid),   &
  	          curr_tile%pid(1:count,ezoldpid), curr_tile%pid(1:count,bxoldpid),	  &
	          curr_tile%pid(1:count,byoldpid), curr_tile%pid(1:count,bzoldpid),   &
              curr_tile%part_ex, curr_tile%part_ey, curr_tile%part_ez, 		      &
              curr_tile%part_bx, curr_tile%part_by, curr_tile%part_bz, 		      &
              curr%charge, curr%mass, dtt)

			  ! Store fields on particle position for next step
			  ! to compute dE/dt and dB/dt terms in LL RR force
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
            !!!! --- push particle species positions a time step
            CALL pxr_pushxyz(count, curr_tile%part_x, curr_tile%part_y,               &
            curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy,                   &
            curr_tile%part_uz, curr_tile%part_gaminv, dtt)
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

#if PROFILING==3
  CALL stop_collection()
#endif

  IF (it.ge.timestat_itstart) THEN
    tend=MPI_WTIME()
    localtimes(1) = localtimes(1) + (tend-tdeb)
  ENDIF
  pushtime=pushtime+(MPI_WTIME()-tdeb)

END SUBROUTINE field_gathering_plus_particle_pusher_sub


! ________________________________________________________________________________________
!> @brief
!> Particle pusher in 3D called by the main function push_particle for the subroutines
!> with cache blocking
!
!> @author
!> Mathieu Lobet

!> @date
!> Creation 2015
!> Revision 10.06.2015
!
!>
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] bxg, byg, bzg electric field grids
!> @param[in] nxx, nyy, nzz number of cells in each direction for the grids
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> for the grids
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells for the current grids
!> @param[in] noxx, noyy, nozz interpolation orders
!> @param[in] dxx, dyy, dzz space steps
!> @param[in] dtt time step
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!>
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher_cacheblock_sub(exg, eyg, ezg, bxg,    &
  byg, bzg, nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard, nyjguard, nzjguard,   &
  noxx, noyy, nozz, dxx, dyy, dzz, dtt, l_lower_order_in_v_in)
  USE grid_tilemodule, ONLY: aofgrid_tiles
  USE mpi
  USE output_data, ONLY: pushtime
  USE params, ONLY: fieldgathe, it, lvec_fieldgathe
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

  ! ___ Parameter declaration __________________________________________
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard,     &
  nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  LOGICAL(lp)                   :: l_lower_order_in_v_in
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
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin, ipmax, ip
  INTEGER(idp)             :: nxjg, nyjg, nzjg
  INTEGER(idp)             :: nxt, nyt, nzt
  INTEGER(idp)             :: nxt_o, nyt_o, nzt_o
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile, bytile, bztile
  LOGICAL(lp)              :: isgathered=.FALSE.

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==3
  CALL start_vtune_collection()
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL  DEFAULT(NONE) SHARED(ntilex,                                        &
  !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, nxjguard, nyjguard,  &
  !$OMP nzjguard, nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, dxx, dyy,  &
  !$OMP dzz, dtt, noxx, noyy, nozz, c_dim, lvec_fieldgathe, l_lower_order_in_v)       &
  !$OMP PRIVATE(ix, iy, iz, ispecies, curr, curr_tile, count, jmin, jmax,             &
  !$OMP extile, eytile, eztile, bxtile, bytile, bztile, nxt, nyt, nzt,                &
  !$OMP kmin, kmax, lmin, lmax, nxc, nyc, nzc, ipmin, ipmax, ip, nxjg, nyjg, nzjg,    &
  !$OMP isgathered) FIRSTPRIVATE(nxt_o, nyt_o, nzt_o)
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
        isgathered=.FALSE.

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
            curr_tile%part_ex(1:count)=0.0_num
            curr_tile%part_ey(1:count)=0.0_num
            curr_tile%part_ez(1:count)=0.0_num
            curr_tile%part_bx(1:count)=0.0_num
            curr_tile%part_by(1:count)=0.0_num
            curr_tile%part_bz(1:count)=0.0_num

            IF ((noxx.eq.1).and.(noyy.eq.1).and.(nozz.eq.1)) THEN

              !!! ---- Loop by blocks over particles in a tile (blocking)
              CALL field_gathering_plus_particle_pusher_1_1_1(count,                  &
              curr_tile%part_x, curr_tile%part_y, curr_tile%part_z,                   &
              curr_tile%part_ux, curr_tile%part_uy, curr_tile%part_uz,                &
              curr_tile%part_gaminv, curr_tile%part_ex, curr_tile%part_ey,            &
              curr_tile%part_ez, curr_tile%part_bx, curr_tile%part_by,                &
              curr_tile%part_bz, curr_tile%x_grid_tile_min,                           &
              curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,    &
              dtt, curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                  &
              curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, extile, eytile,              &
              eztile, bxtile, bytile, bztile, curr%charge,                            &
              curr%mass, lvec_fieldgathe, l_lower_order_in_v)

            ELSE IF ((noxx.eq.2).and.(noyy.eq.2).and.(nozz.eq.2)) THEN

              !!! ---- Loop by blocks over particles in a tile (blocking)
              CALL field_gathering_plus_particle_pusher_2_2_2(count,                  &
              curr_tile%part_x, curr_tile%part_y, curr_tile%part_z,                   &
              curr_tile%part_ux, curr_tile%part_uy, curr_tile%part_uz,                &
              curr_tile%part_gaminv, curr_tile%part_ex, curr_tile%part_ey,            &
              curr_tile%part_ez, curr_tile%part_bx, curr_tile%part_by,                &
              curr_tile%part_bz, curr_tile%x_grid_tile_min,                           &
              curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,    &
              dtt, curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                  &
              curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, extile, eytile,              &
              eztile, bxtile, bytile, bztile, curr%charge,                            &
              curr%mass, lvec_fieldgathe, l_lower_order_in_v)

            ELSE IF ((noxx.eq.3).and.(noyy.eq.3).and.(nozz.eq.3)) THEN

              !!! ---- Loop by blocks over particles in a tile (blocking)
              CALL field_gathering_plus_particle_pusher_3_3_3(count,                  &
              curr_tile%part_x, curr_tile%part_y, curr_tile%part_z,                   &
              curr_tile%part_ux, curr_tile%part_uy, curr_tile%part_uz,                &
              curr_tile%part_gaminv, curr_tile%part_ex, curr_tile%part_ey,            &
              curr_tile%part_ez, curr_tile%part_bx, curr_tile%part_by,                &
              curr_tile%part_bz, curr_tile%x_grid_tile_min,                           &
              curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,    &
              dtt, curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                  &
              curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, extile, eytile,              &
              eztile, bxtile, bytile, bztile, curr%charge,                            &
              curr%mass, lvec_fieldgathe, l_lower_order_in_v)

            ENDIF

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

  IF (it.ge.timestat_itstart) THEN
    tend=MPI_WTIME()
    localtimes(1) = localtimes(1) + (tend-tdeb)
  ENDIF
  pushtime=pushtime+(tend-tdeb)
END SUBROUTINE field_gathering_plus_particle_pusher_cacheblock_sub

! ________________________________________________________________________________________
!> @brief
!> Particle pusher in 3D called by the main function push_particle
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 10.06.2015
!>
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] bxg, byg, bzg electric field grids
!> @param[in] nxx, nyy, nzz number of cells in each direction for the grids
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> for the grids
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells for the current grids
!> @param[in] noxx, noyy, nozz interpolation orders
!> @param[in] dxx, dyy, dzz space steps
!> @param[in] dtt time step
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!
! ________________________________________________________________________________________
SUBROUTINE particle_pusher_sub(exg, eyg, ezg, bxg, byg, bzg, nxx, nyy, nzz, nxguard,  &
  nyguard, nzguard, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy, dzz, dtt, &
  l_lower_order_in_v_in)
  USE grid_tilemodule, ONLY: aofgrid_tiles
  USE mpi
  USE output_data, ONLY: pushtime
  USE particle_properties, ONLY: nspecies, particle_pusher
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  USE time_stat, ONLY: localtimes, timestat_itstart
  ! Vtune/SDE profiling
  IMPLICIT NONE
  ! ___ Parameter declaration __________________________________________
  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard,     &
  nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  LOGICAL(lp)                   :: l_lower_order_in_v_in
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
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin, ipmax, ip
  INTEGER(idp)             :: nxjg, nyjg, nzjg
  INTEGER(idp)             :: nxt, nyt, nzt
  INTEGER(idp)             :: nxt_o, nyt_o, nzt_o
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile, bytile, bztile
  LOGICAL(lp)              :: isgathered=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN

  tdeb=MPI_WTIME()

#if VTUNE==3
  CALL start_vtune_collection()
#endif


  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex,                                         &
  !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, nxjguard, nyjguard,  &
  !$OMP nzjguard, nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, dxx, dyy,  &
  !$OMP dzz, dtt, noxx, noyy, nozz, c_dim, particle_pusher) PRIVATE(ix, iy, iz,       &
  !$OMP ispecies, curr, curr_tile, count, jmin, jmax, kmin, kmax, lmin, lmax,         &
  !$OMP extile, eytile, eztile, bxtile, bytile, bztile, nxt, nyt, nzt,                &
  !$OMP nxc, nyc, nzc, ipmin, ipmax, ip, nxjg, nyjg, nzjg, isgathered,                &
  !$OMP nxt_o, nyt_o, nzt_o)
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
        isgathered=.FALSE.
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
            SELECT CASE (particle_pusher)
              !! Vay pusher -- Full push
            CASE (1_idp)
              CALL pxr_ebcancelpush3d(count, curr_tile%part_ux, curr_tile%part_uy,    &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt,      &
              0_idp)

              !! Boris pusher with RR (S09 model, according to VRANIC2016, https://doi.org/10.1016/j.cpc.2016.04.002)-- Full push
	    CASE (2_idp)
              CALL pxr_boris_push_rr_S09_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,&
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)

              !! Boris pusher with RR (B08 model, according to VRANIC2016, https://doi.org/10.1016/j.cpc.2016.04.002)model-- Full push
	    CASE (3_idp)
              CALL pxr_boris_push_rr_B08_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,&
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)
	      !! Boris pusher with RR (LL model, according to VRANIC2016, https://doi.org/10.1016/j.cpc.2016.04.002)-- Full push
	    CASE (4_idp)
              CALL pxr_boris_push_rr_LL_u_3d(count, curr_tile%part_ux, 		      &
              curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv, 	      &
              curr_tile%pid(1:count,4), curr_tile%pid(1:count,5),		      &
  	      curr_tile%pid(1:count,6), curr_tile%pid(1:count,7),		      &
	      curr_tile%pid(1:count,8), curr_tile%pid(1:count,9),                     &
              curr_tile%part_ex, curr_tile%part_ey, curr_tile%part_ez, 		      &
              curr_tile%part_bx, curr_tile%part_by, curr_tile%part_bz, 		      &
              curr%charge, curr%mass, dtt)

              !! Boris pusher -- Full push
            CASE DEFAULT

              !! Push momentum using the Boris method in a single subroutine
              CALL pxr_boris_push_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,   &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)
            END SELECT
            !!!! --- push particle species positions a time step
            CALL pxr_pushxyz(count, curr_tile%part_x, curr_tile%part_y,               &
            curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy,                   &
            curr_tile%part_uz, curr_tile%part_gaminv, dtt)
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

  tend=MPI_WTIME()
  pushtime=pushtime+(tend-tdeb)
  IF (it.ge.timestat_itstart) THEN
    localtimes(1) = localtimes(1) + (tend-tdeb)
  ENDIF

#if defined(DEBUG)
  WRITE(0, *) "Push_particles: stop"
#endif

END SUBROUTINE particle_pusher_sub

! ________________________________________________________________________________________
!> @brief
!>  Field gathering+ (E & B) Push half a time step
!
!> @details
!> This subroutine is wrapper for pxrpush_particles_part1_sub()
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_particles_part1
  USE fields, ONLY: bx_p, by_p, bz_p, ex_p, ey_p, ez_p, l4symtry,                    &
    l_lower_order_in_v, nox, noy, noz, nxguards, nxjguards, nyguards, nyjguards,     &
    nzguards, nzjguards
  USE params, ONLY: dt, fieldgathe, lvec_fieldgathe
  USE shared_data, ONLY: dx, dy, dz, nx, ny, nz
  IMPLICIT NONE

#if defined(DEBUG)
  WRITE(0, *) "pxrpush_particles_part1: start"
#endif

  CALL pxrpush_particles_part1_sub(ex_p, ey_p, ez_p, bx_p, by_p, bz_p, nx, ny, nz, &
  nxguards,  nyguards, nzguards, nxjguards, nyjguards, nzjguards, nox, noy, noz,   &
  dx, dy, dz, dt, l4symtry, l_lower_order_in_v, LVEC_fieldgathe, fieldgathe)

#if defined(DEBUG)
  WRITE(0, *) "pxrpush_particles_part1: stop"
#endif

END SUBROUTINE pxrpush_particles_part1

! ________________________________________________________________________________________
!> @brief
!> Perform the field gathering + (E & B) Push half a time step
!
!> @details
!> This subroutine is called in pxrpush_particles_part1()
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!>
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] bxg, byg, bzg electric field grids
!> @param[in] nxx, nyy, nzz number of cells in each direction for the grids
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction for the grids
!> @param[in] nxjguard, nyjguard, nzjguard number of guard cells for the current grids
!> @param[in] noxx, noyy, nozz interpolation orders
!> @param[in] dxx, dyy, dzz space steps
!> @param[in] dtt time step
!> @param[in] l4symtry_in
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!>
! ________________________________________________________________________________________
SUBROUTINE pxrpush_particles_part1_sub(exg, eyg, ezg, bxg, byg, bzg, nxx, nyy, nzz,   &
  nxguard, nyguard, nzguard, nxjguard, nyjguard, nzjguard, noxx, noyy, nozz, dxx, dyy,&
  dzz, dtt, l4symtry_in, l_lower_order_in_v_in, lvect, field_gathe_algo)
  USE grid_tilemodule, ONLY: aofgrid_tiles
  USE particle_properties, ONLY: nspecies, particle_pusher
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, lp, num
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  IMPLICIT NONE

  INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxguard, nyguard, nzguard, nxjguard,     &
  nyjguard, nzjguard
  INTEGER(idp), INTENT(IN) :: noxx, noyy, nozz
  INTEGER(idp), INTENT(IN) :: lvect, field_gathe_algo
  LOGICAL(lp), INTENT(IN) :: l4symtry_in, l_lower_order_in_v_in
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
  INTEGER(idp)          :: ispecies, ix, iy, iz, count
  INTEGER(idp)          :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  INTEGER(idp)                    :: nxc, nyc, nzc
  INTEGER(idp)                    :: nxjg, nyjg, nzjg
  INTEGER(idp)             :: nxt, nyt, nzt
  INTEGER(idp)             :: nxt_o, nyt_o, nzt_o
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile, bytile, bztile
  LOGICAL(lp)                     :: isgathered=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DEFAULT(NONE) SHARED(ntilex,                                         &
  !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, nxjguard, nyjguard,  &
  !$OMP nzjguard, exg, eyg, ezg, bxg, byg, bzg, dxx, dyy, dzz, dtt, noxx, noyy,       &
  !$OMP l4symtry_in, l_lower_order_in_v_in, nozz, c_dim, fieldgathe, particle_pusher, &
  !$OMP field_gathe_algo, lvect) PRIVATE(ix, iy, iz, ispecies, curr, curr_tile,       &
  !$OMP extile, eytile, eztile, bxtile, bytile, bztile, nxt, nyt, nzt,                &
  !$OMP count, jmin, jmax, kmin, kmax, lmin, lmax, nxc, nyc, nzc, nxjg, nyjg,         &
  !$OMP nzjg, isgathered, nxt_o, nyt_o, nzt_o)
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
        isgathered=.FALSE.
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

            IF (field_gathe_algo.gt.-1) then

              curr_tile%part_ex(1:count) = 0.0_num
              curr_tile%part_ey(1:count) = 0.0_num
              curr_tile%part_ez(1:count) = 0.0_num
              curr_tile%part_bx(1:count) = 0.0_num
              curr_tile%part_by(1:count) = 0.0_num
              curr_tile%part_bz(1:count) = 0.0_num

              !!! ---- Loop by blocks over particles in a tile (blocking)
              SELECT CASE (c_dim)
              CASE (2)! 2D CASE
                !!! --- Gather electric and magnetic fields on particles
                CALL geteb2dxz_energy_conserving(count, curr_tile%part_x,             &
                curr_tile%part_y, curr_tile%part_z, curr_tile%part_ex,                &
                curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,              &
                curr_tile%part_by, curr_tile%part_bz, curr_tile%x_grid_tile_min,      &
                curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,  &
                curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                     &
                curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, noxx, noyy, nozz,          &
                extile, eytile, eztile, bxtile, bytile,                               &
                bztile , l4symtry_in, l_lower_order_in_v_in, lvect,                   &
                field_gathe_algo)
              CASE DEFAULT! 3D CASE
                !!! --- Gather electric and magnetic fields on particles
                CALL geteb3d_energy_conserving(count, curr_tile%part_x,               &
                curr_tile%part_y, curr_tile%part_z, curr_tile%part_ex,                &
                curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,              &
                curr_tile%part_by, curr_tile%part_bz, curr_tile%x_grid_tile_min,      &
                curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min, dxx, dyy, dzz,  &
                curr_tile%nx_cells_tile, curr_tile%ny_cells_tile,                     &
                curr_tile%nz_cells_tile, nxjg, nyjg, nzjg, noxx, noyy, nozz,          &
                extile, eytile, eztile, bxtile, bytile,                               &
                bztile , l4symtry_in, l_lower_order_in_v_in, lvect,                   &
                field_gathe_algo)
              END SELECT

            END IF

            SELECT CASE (particle_pusher)
              !! Vay pusher -- half push part 1

            CASE (1_idp)
              CALL pxr_ebcancelpush3d(count, curr_tile%part_ux, curr_tile%part_uy,    &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt,      &
              1_idp)
              !! Boris pusher -- half push part 1
            CASE DEFAULT
              !! --- Push velocity with E half step
              CALL pxr_epush_v(count, curr_tile%part_ux, curr_tile%part_uy,           &
              curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,                &
              curr_tile%part_ez, curr%charge, curr%mass, dtt*0.5_num)
              !! --- Set gamma of particles
              CALL pxr_set_gamma(count, curr_tile%part_ux, curr_tile%part_uy,         &
              curr_tile%part_uz, curr_tile%part_gaminv)
              !! --- Push velocity with B half step
              CALL pxr_bpush_v(count, curr_tile%part_ux, curr_tile%part_uy,           &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_bx,            &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass,           &
              dtt*0.5_num)
            END SELECT
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

END SUBROUTINE pxrpush_particles_part1_sub

! ________________________________________________________________________________________
!> @brief
!> (B & E) Push half a time step + XYZ push half a time step
!>
!> @author
!> Henri Vincenti
!>
!> @date
!> Creation 2015
!> Revision 06.10.2016
! ________________________________________________________________________________________
SUBROUTINE pxrpush_particles_part2
  USE fields, ONLY: bx, by, bz, ex, ey, ez, nxjguards, nyjguards, nzjguards
  USE mpi
  USE output_data, ONLY: pushtime
  USE params, ONLY: dt
  USE particle_properties, ONLY: nspecies, particle_pusher
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: c_dim, dx, dy, dz, x, y, z
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  USE tiling
  IMPLICIT NONE
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                       :: tdeb, tend
  INTEGER(idp)                    :: nxc, nyc, nzc

#if defined(DEBUG)
  WRITE(0, *) "pxrpush_particles_part2: start"
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  tdeb=MPI_WTIME()
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) SHARED(ntilex,        &
  !$OMP ntiley, ntilez, nspecies, species_parray, nxjguards, nyjguards, nzjguards,    &
  !$OMP ex, ey, ez, bx, by, bz, dx, dy, dz, dt, c_dim, particle_pusher) PRIVATE(ix,   &
  !$OMP iy, iz, ispecies, curr, curr_tile, count, jmin, jmax, kmin, kmax, lmin, lmax, &
  !$OMP nxc, nyc, nzc)
  DO iz=1, ntilez! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        DO ispecies=1, nspecies! LOOP ON SPECIES
          ! - Get current tile properties
          ! - Init current tile variables
          curr=>species_parray(ispecies)
          IF (curr%is_antenna) CYCLE
          IF (curr%lfreeze) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) CYCLE
          SELECT CASE (particle_pusher)
            !! Vay pusher -- half push part 2
          CASE (1_idp)
            CALL pxr_ebcancelpush3d(count, curr_tile%part_ux, curr_tile%part_uy,      &
            curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,              &
            curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                  &
            curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dt, 2_idp)
          CASE DEFAULT
            !! Boris pusher -- half push part 2
            !!! --- Push velocity with B half step
            CALL pxr_bpush_v(count, curr_tile%part_ux, curr_tile%part_uy,             &
            curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_bx,              &
            curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dt*0.5_num)
            !! --- Push velocity with E half step
            CALL pxr_epush_v(count, curr_tile%part_ux, curr_tile%part_uy,             &
            curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,                  &
            curr_tile%part_ez, curr%charge, curr%mass, dt*0.5_num)
            !! --- Sets gamma of particles
            CALL pxr_set_gamma(count, curr_tile%part_ux, curr_tile%part_uy,           &
            curr_tile%part_uz, curr_tile%part_gaminv)
          END SELECT

          SELECT CASE (c_dim)
          CASE (2)! 2D CASE
            !! --- Advance particle position of one time step
            CALL pxr_pushxz(count, curr_tile%part_x, curr_tile%part_z,                &
            curr_tile%part_ux, curr_tile%part_uz, curr_tile%part_gaminv, dt)
          CASE DEFAULT! 3D CASE
            !! --- Advance particle position of one time step
            CALL pxr_pushxyz(count, curr_tile%part_x, curr_tile%part_y,               &
            curr_tile%part_z, curr_tile%part_ux, curr_tile%part_uy,                   &
            curr_tile%part_uz, curr_tile%part_gaminv, dt)
          END SELECT
        END DO! END LOOP ON SPECIES
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO
  tend=MPI_WTIME()
  pushtime=pushtime+(tend-tdeb)

#if defined(DEBUG)
  WRITE(0, *) "pxrpush_particles_part2: stop"
#endif

END SUBROUTINE pxrpush_particles_part2



! ________________________________________________________________________________________
!> @brief
!> This function combined the field gathering and the particle pusher
!> in 3D for CIC (order 1) shape factor.
!> The field gathering and the particle pusher are done in the same particle loop.
!
!> @details
!> This function is vectorized.

!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revision 10.06.2015
!
!> @warning
!> Only l_lower_order_in_v=True is implemented
!
! Input parameters:
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[in] uxp, uyp, uzp particle momentum
!> @param[in] gaminv inverse of the particle Lorentz factor
!> @param[in] ex, ey, ez particle electric field
!> @param[in] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dtt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v performe the field interpolation at a lower order
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher_1_1_1(np, xp, yp, zp, uxp, uyp, uzp,  &
  gaminv, ex, ey, ez, bx, by, bz, xmin, ymin, zmin, dx, dy, dz, dtt, nx, ny, nz,        &
  nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, q, m, lvect,                 &
  l_lower_order_in_v)
  USE constants, ONLY: clight
  USE omp_lib
  USE params, ONLY: dt
  USE particle_properties, ONLY: particle_pusher
  USE picsar_precision, ONLY: idp, isp, lp, num

  IMPLICIT NONE

  ! ___ Parameter declaration ____________________________________

  ! Input/Output parameters
  INTEGER(idp), INTENT(IN)                :: np, nx, ny, nz, nxguard, nyguard,        &
  nzguard
  INTEGER(idp), INTENT(IN)                :: lvect
  REAL(num), INTENT(IN)                   :: q, m
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp, yp, zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: ex, ey, ez
  REAL(num), DIMENSION(np), INTENT(INOUT) :: bx, by, bz
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp, uyp, uzp, gaminv
  LOGICAL(lp), INTENT(IN)                 :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN)                              :: exg, eyg, ezg,     &
  bxg, byg, bzg
  REAL(num), INTENT(IN)                   :: xmin, ymin, zmin, dx, dy, dz, dtt

  ! Local parameters
  INTEGER(isp)                         :: j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn
  INTEGER(idp)                         :: blocksize
  REAL(num)                            :: dxi, dyi, dzi
  REAL(num)                            :: x, y, z
  REAL(num)                            :: a
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: clghtisq, const1
  REAL(num), DIMENSION(0:1)            :: sx, sy, sz
  REAL(num), DIMENSION(0:1)            :: sx0, sy0, sz0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: sx, sy, sz, sx0, sy0, sz0
#endif
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  ! ___ Parameter initialization ____________________________________
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  const1 = 0.5_num*q*dtt/m
  clghtisq = 1.0_num/clight**2

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  sx0(0) = 1.0_num
  sy0(0) = 1.0_num
  sz0(0) = 1.0_num

  ! ____________________________________________________________________________
  ! Loop on block of particles of size lvect
  DO ip=1, np, lvect

    blocksize = MIN(lvect, np-ip+1)

    ! __________________________________________________________________________
    ! Field gathering

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD private(sx, sy, sz, sx0, sy0, sz0)
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO nn=ip, blocksize+ip-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x)
      k=floor(y)
      k0=floor(y)
      l=floor(z)
      l0=floor(z)

      xint=x-j
      yint=y-k
      zint=z-l

      ! Compute shape factors
      sx(0) = 1.0_num-xint
      sx(1) = xint
      sy(0) = 1.0_num-yint
      sy(1) = yint
      sz(0) = 1.0_num-zint
      sz(1) = zint

      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0

      ! Compute Ex on particle
      a = (sy(0)*exg(j0, k, l)   + sy(1)*exg(j0, k+1, l))*sz(0) + (sy(0)*exg(j0, k,   &
      l+1) + sy(1)*exg(j0, k+1, l+1))*sz(1)
      ex(nn) = ex(nn) + a

      ! Compute Ey on particle
      a = (sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0, l))*sz(0) + (sx(0)*eyg(j, k0,     &
      l+1) + sx(1)*eyg(j+1, k0, l+1))*sz(1)
      ey(nn) = ey(nn) + a

      ! Compute Ez on particle
      a = (sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k, l0))*sy(0) + (sx(0)*ezg(j, k+1,    &
      l0) + sx(1)*ezg(j+1, k+1, l0))*sy(1)
      ez(nn) = ez(nn) + a

      ! Compute Bx on particle
      a = (sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0, l0))
      bx(nn) = bx(nn) + a

      ! Compute By on particle
      a = (sy(0)*byg(j0, k, l0) + sy(1)*byg(j0, k+1, l0))
      by(nn) = by(nn) + a

      ! Compute Bz on particle
      a = (sz(0)*bzg(j0, k0, l) + sz(1)*bzg(j0, k0, l+1))
      bz(nn) = bz(nn) + a

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! __________________________________________________________________________
    ! Particle pusher

    SELECT CASE (particle_pusher)
      !! Vay pusher -- Full push
    CASE (1_idp)
      CALL pxr_ebcancelpush3d(blocksize, uxp(ip:ip+blocksize-1),                      &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt,  &
      0_idp)

      !! Boris pusher -- Full push
    CASE DEFAULT

      !! Push momentum using the Boris method in a single subroutine

      CALL pxr_boris_push_u_3d(blocksize, uxp(ip:ip+blocksize-1),                     &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize, uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),   &
      uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1))

    END SELECT
    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize, xp(ip:ip+blocksize-1), yp(ip:ip+blocksize-1),         &
    zp(ip:ip+blocksize-1), uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),            &
    uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1), dt)
  ENDDO

  RETURN
END SUBROUTINE field_gathering_plus_particle_pusher_1_1_1


! ________________________________________________________________________________________
!
!> @brief
!> This function combined the field gathering and the particle pusher
!> in 3D for order 2 shape factor.
!
!> @details
!> The field gathering and the particle pusher are done in the same particle loop.
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revision 10.06.2015
!
!> @warning
!> Only l_lower_order_in_v=True is implemented
!
! Input parameters:
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[in] uxp, uyp, uzp particle momentum
!> @param[in] gaminv inverse of the particle Lorentz factor
!> @param[in] ex, ey, ez particle electric field
!> @param[in] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dtt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v performe the field interpolation at a lower order
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher_2_2_2(np, xp, yp, zp, uxp, uyp, uzp,  &
  gaminv, ex, ey, ez, bx, by, bz, xmin, ymin, zmin, dx, dy, dz, dtt, nx, ny, nz,        &
  nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, q, m, lvect,                 &
  l_lower_order_in_v)
  USE constants, ONLY: clight
  USE omp_lib
  USE params, ONLY: dt
  USE particle_properties, ONLY: particle_pusher
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN)                :: np, nx, ny, nz, nxguard, nyguard,        &
  nzguard
  INTEGER(idp), INTENT(IN)                :: lvect
  REAL(num), INTENT(IN)                   :: q, m
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp, yp, zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: ex, ey, ez
  REAL(num), DIMENSION(np), INTENT(INOUT) :: bx, by, bz
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp, uyp, uzp, gaminv
  LOGICAL(lp), INTENT(IN)                 :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN)                              :: exg, eyg, ezg,     &
  bxg, byg, bzg
  REAL(num), INTENT(IN)                   :: xmin, ymin, zmin, dx, dy, dz, dtt

  ! Local parameters
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn
  INTEGER(idp)                         :: blocksize
  INTEGER(isp)                         :: j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq, yintsq, zintsq
  REAL(num)                            :: a
  REAL(num)                            :: clghtisq, const1
  REAL(num), DIMENSION(-1:1)           :: sx, sy, sz
  REAL(num), DIMENSION(-1:1)           :: sx0, sy0, sz0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: sx, sy, sz, sx0, sy0, sz0
#endif
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  const1 = 0.5_num*q*dtt/m

  clghtisq = 1.0_num/clight**2

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  ! ___ Loop on particles _______________________
  DO ip=1, np, lvect

    blocksize = MIN(lvect, np-ip+1)

    ! __________________________________________________________________________
    ! Field gathering

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD private(sx, sy, sz, sx0, sy0, sz0)
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO nn=ip, blocksize+ip-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! Compute index of particle
      j=nint(x)
      j0=floor(x-0.5_num)
      k=nint(y)
      k0=floor(y-0.5_num)
      l=nint(z)
      l0=floor(z-0.5_num)

      xint=x-j
      yint=y-k
      zint=z-l

      ! Compute shape factors
      xintsq = xint*xint
      sx(-1) = 0.5_num*(0.5_num-xint)**2
      sx( 0) = 0.75_num-xintsq
      sx( 1) = 0.5_num*(0.5_num+xint)**2

      yintsq = yint*yint
      sy(-1) = 0.5_num*(0.5_num-yint)**2
      sy( 0) = 0.75_num-yintsq
      sy( 1) = 0.5_num*(0.5_num+yint)**2

      zintsq = zint*zint
      sz(-1) = 0.5_num*(0.5_num-zint)**2
      sz( 0) = 0.75_num-zintsq
      sz( 1) = 0.5_num*(0.5_num+zint)**2

      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0

      sx0( 0) = 1.0_num-xint
      sx0( 1) = xint

      sy0( 0) = 1.0_num-yint
      sy0( 1) = yint

      sz0( 0) = 1.0_num-zint
      sz0( 1) = zint

      ! Compute Ex on particle
      a = (sx0(0)*exg(j0, k-1, l-1) + sx0(1)*exg(j0+1, k-1, l-1))*sy(-1)
      a = a + (sx0(0)*exg(j0, k, l-1) + sx0(1)*exg(j0+1, k, l-1))*sy(0)
      a = a + (sx0(0)*exg(j0, k+1, l-1) + sx0(1)*exg(j0+1, k+1, l-1))*sy(1)
      ex(nn) = ex(nn) + a*sz(-1)
      a = (sx0(0)*exg(j0, k-1, l) + sx0(1)*exg(j0+1, k-1, l))*sy(-1)
      a = a + (sx0(0)*exg(j0, k, l) + sx0(1)*exg(j0+1, k, l))*sy(0)
      a = a + (sx0(0)*exg(j0, k+1, l) + sx0(1)*exg(j0+1, k+1, l))*sy(1)
      ex(nn) = ex(nn) + a*sz(0)
      a = (sx0(0)*exg(j0, k-1, l+1) + sx0(1)*exg(j0+1, k-1, l+1))*sy(-1)
      a = a + (sx0(0)*exg(j0, k, l+1) + sx0(1)*exg(j0+1, k, l+1))*sy(0)
      a = a + (sx0(0)*exg(j0, k+1, l+1) + sx0(1)*exg(j0+1, k+1, l+1))*sy(1)
      ex(nn) = ex(nn) + a*sz(1)

      ! Compute Ey on particle
      a = (sx(-1)*eyg(j-1, k0, l-1) + sx(0)*eyg(j, k0, l-1) + sx(1)*eyg(j+1, k0,      &
      l-1))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l-1) + sx(0)*eyg(j, k0+1, l-1) + sx(1)*eyg(j+1,  &
      k0+1, l-1))*sy0(1)
      ey(nn) = ey(nn) + a*sz(-1)
      a = (sx(-1)*eyg(j-1, k0, l) + sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0,          &
      l))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l) + sx(0)*eyg(j, k0+1, l) + sx(1)*eyg(j+1,      &
      k0+1, l))*sy0(1)
      ey(nn) = ey(nn) + a*sz(0)
      a = (sx(-1)*eyg(j-1, k0, l+1) + sx(0)*eyg(j, k0, l+1) + sx(1)*eyg(j+1, k0,      &
      l+1))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l+1) + sx(0)*eyg(j, k0+1, l+1) + sx(1)*eyg(j+1,  &
      k0+1, l+1))*sy0(1)
      ey(nn) = ey(nn) + a*sz(1)

      ! Compute Ez on particle
      a = (sx(-1)*ezg(j-1, k-1, l0) + sx(0)*ezg(j, k-1, l0) + sx(1)*ezg(j+1, k-1,     &
      l0))*sy(-1)
      a = a + (sx(-1)*ezg(j-1, k, l0) + sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k,       &
      l0))*sy(0)
      a = a + (sx(-1)*ezg(j-1, k+1, l0) + sx(0)*ezg(j, k+1, l0) + sx(1)*ezg(j+1, k+1, &
      l0))*sy(1)
      ez(nn) = ez(nn) + a*sz0(0)
      a = (sx(-1)*ezg(j-1, k-1, l0+1) + sx(0)*ezg(j, k-1, l0+1) + sx(1)*ezg(j+1, k-1, &
      l0+1))*sy(-1)
      a = a + (sx(-1)*ezg(j-1, k, l0+1) + sx(0)*ezg(j, k, l0+1) + sx(1)*ezg(j+1, k,   &
      l0+1))*sy(0)
      a = a + (sx(-1)*ezg(j-1, k+1, l0+1) + sx(0)*ezg(j, k+1, l0+1) + sx(1)*ezg(j+1,  &
      k+1, l0+1))*sy(1)
      ez(nn) = ez(nn) + a*sz0(1)

      ! Compute Bx on particle
      a = (sx(-1)*bxg(j-1, k0, l0) + sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0,        &
      l0))*sy0(0)
      a = a + (sx(-1)*bxg(j-1, k0+1, l0) + sx(0)*bxg(j, k0+1, l0) + sx(1)*bxg(j+1,    &
      k0+1, l0))*sy0(1)
      bx(nn) = bx(nn) + a*sz0(0)
      a = (sx(-1)*bxg(j-1, k0, l0+1) + sx(0)*bxg(j, k0, l0+1) + sx(1)*bxg(j+1, k0,    &
      l0+1))*sy0(0)
      a = a + (sx(-1)*bxg(j-1, k0+1, l0+1) + sx(0)*bxg(j, k0+1, l0+1) +               &
      sx(1)*bxg(j+1, k0+1, l0+1))*sy0(1)
      bx(nn) = bx(nn) + a*sz0(1)

      ! Compute By on particle
      a = (sx0(0)*byg(j0, k-1, l0) + sx0(1)*byg(j0+1, k-1, l0))*sy(-1)
      a = a + (sx0(0)*byg(j0, k, l0) + sx0(1)*byg(j0+1, k, l0))*sy(0)
      a = a + (sx0(0)*byg(j0, k+1, l0) + sx0(1)*byg(j0+1, k+1, l0))*sy(1)
      by(nn) = by(nn) + a*sz0(0)
      a = (sx0(0)*byg(j0, k-1, l0+1) + sx0(1)*byg(j0+1, k-1, l0+1))*sy(-1)
      a = a + (sx0(0)*byg(j0, k, l0+1) + sx0(1)*byg(j0+1, k, l0+1))*sy(0)
      a = a + (sx0(0)*byg(j0, k+1, l0+1) + sx0(1)*byg(j0+1, k+1, l0+1))*sy(1)
      by(nn) = by(nn) + a*sz0(1)

      ! Compute Bz on particle
      a = (sx0(0)*bzg(j0, k0, l-1) + sx0(1)*bzg(j0+1, k0, l-1))*sy0(0)
      a = a + (sx0(0)*bzg(j0, k0+1, l-1) + sx0(1)*bzg(j0+1, k0+1, l-1))*sy0(1)
      bz(nn) = bz(nn) + a*sz(-1)
      a = (sx0(0)*bzg(j0, k0, l) + sx0(1)*bzg(j0+1, k0, l))*sy0(0)
      a = a + (sx0(0)*bzg(j0, k0+1, l) + sx0(1)*bzg(j0+1, k0+1, l))*sy0(1)
      bz(nn) = bz(nn) + a*sz(0)
      a = (sx0(0)*bzg(j0, k0, l+1) + sx0(1)*bzg(j0+1, k0, l+1))*sy0(0)
      a = a + (sx0(0)*bzg(j0, k0+1, l+1) + sx0(1)*bzg(j0+1, k0+1, l+1))*sy0(1)
      bz(nn) = bz(nn) + a*sz(1)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! __________________________________________________________________________
    ! Particle pusher

    SELECT CASE (particle_pusher)
      !! Vay pusher -- Full push
    CASE (1_idp)
      CALL pxr_ebcancelpush3d(blocksize, uxp(ip:ip+blocksize-1),                      &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt,  &
      0_idp)

      !! Boris pusher -- Full push
    CASE DEFAULT

      !! Push momentum using the Boris method in a single subroutine

      CALL pxr_boris_push_u_3d(blocksize, uxp(ip:ip+blocksize-1),                     &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize, uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),   &
      uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1))

    END SELECT
    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize, xp(ip:ip+blocksize-1), yp(ip:ip+blocksize-1),         &
    zp(ip:ip+blocksize-1), uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),            &
    uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1), dt)
  ENDDO

  RETURN

END SUBROUTINE field_gathering_plus_particle_pusher_2_2_2


! ________________________________________________________________________________________
!> @brief
!> This function combined the field gathering and the particle pusher
!> in 3D for CIC (order 1) shape factor.
!> The field gathering and the particle pusher are done in the same particle loop.
!
!> @detail
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revision 12/03/2016 - more efficient gather algorithm for vectorization
!> Revision 10/06/2016 - new pusher
!
!> @warning
!> Only l_lower_order_in_v=True is implemented
!
! Input parameters:
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[in] uxp, uyp, uzp particle momentum
!> @param[in] gaminv inverse of the particle Lorentz factor
!> @param[in] ex, ey, ez particle electric field
!> @param[in] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dtt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v performe the field interpolation at a lower order
!
! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher_3_3_3(np, xp, yp, zp, uxp, uyp, uzp,  &
  gaminv, ex, ey, ez, bx, by, bz, xmin, ymin, zmin, dx, dy, dz, dtt, nx, ny, nz,        &
  nxguard, nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg, q, m, lvect,                 &
  l_lower_order_in_v)
  USE constants, ONLY: clight
  USE omp_lib
  USE params, ONLY: dt
  USE particle_properties, ONLY: particle_pusher
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN)                :: np, nx, ny, nz, nxguard, nyguard,        &
  nzguard
  INTEGER(idp), INTENT(IN)                :: lvect
  REAL(num), INTENT(IN)                   :: q, m
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp, yp, zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: ex, ey, ez
  REAL(num), DIMENSION(np), INTENT(INOUT) :: bx, by, bz
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp, uyp, uzp, gaminv
  LOGICAL(lp), INTENT(IN)                 :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN)                              :: exg, eyg, ezg,     &
  bxg, byg, bzg
  REAL(num), INTENT(IN)                   :: xmin, ymin, zmin, dx, dy, dz, dtt
  ! Local parameters
  INTEGER(isp)                         :: ip
  INTEGER(idp)                         :: blocksize
  INTEGER(isp)                         :: nn
  INTEGER(isp)                         :: j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq, oxint, yintsq, oyint, zintsq,       &
  ozint, oxintsq, oyintsq, ozintsq
  REAL(num)                            :: clghtisq, const1
  REAL(num)                            :: a
  REAL(num), DIMENSION(-1:2)           :: sx
  REAL(num), DIMENSION(-1:2)           :: sy
  REAL(num), DIMENSION(-1:2)           :: sz
  REAL(num), DIMENSION(-1:1)           :: sx0
  REAL(num), DIMENSION(-1:1)           :: sy0
  REAL(num), DIMENSION(-1:1)           :: sz0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: sx, sy, sz, sx0, sy0, sz0
#endif
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  const1 = 0.5_num*q*dtt/m

  clghtisq = 1.0_num/clight**2

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  ! ___ Loop on partciles _______________________
  DO ip=1, np, lvect

    blocksize = MIN(lvect, np-ip+1)

    ! __________________________________________________________________________
    ! Field gathering

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD private(sx, sy, sz, sx0, sy0, sz0)
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO nn=ip, blocksize+ip-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x)
      k=floor(y)
      k0=floor(y)
      l=floor(z)
      l0=floor(z)

      xint=x-j
      yint=y-k
      zint=z-l

      ! Compute shape factors
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0

      xintsq = xint*xint
      sx0(-1) = 0.5_num*(0.5_num-xint)**2
      sx0( 0) = 0.75_num-xintsq
      sx0( 1) = 0.5_num*(0.5_num+xint)**2

      yintsq = yint*yint
      sy0(-1) = 0.5_num*(0.5_num-yint)**2
      sy0( 0) = 0.75_num-yintsq
      sy0( 1) = 0.5_num*(0.5_num+yint)**2

      zintsq = zint*zint
      sz0(-1) = 0.5_num*(0.5_num-zint)**2
      sz0( 0) = 0.75_num-zintsq
      sz0( 1) = 0.5_num*(0.5_num+zint)**2

      ! Compute Ex on particle
      a = (sx0(-1)*exg(j0-1, k-1, l-1) + sx0(0)*exg(j0, k-1, l-1) + sx0(1)*exg(j0+1,  &
      k-1, l-1))*sy(-1)
      a = a + (sx0(-1)*exg(j0-1, k, l-1) + sx0(0)*exg(j0, k, l-1) + sx0(1)*exg(j0+1,  &
      k, l-1))*sy(0)
      a = a + (sx0(-1)*exg(j0-1, k+1, l-1) + sx0(0)*exg(j0, k+1, l-1) +               &
      sx0(1)*exg(j0+1, k+1, l-1))*sy(1)
      a = a + (sx0(-1)*exg(j0-1, k+2, l-1) + sx0(0)*exg(j0, k+2, l-1) +               &
      sx0(1)*exg(j0+1, k+2, l-1))*sy(2)
      ex(nn) = ex(nn) + a*sz(-1)
      a = (sx0(-1)*exg(j0-1, k-1, l) + sx0(0)*exg(j0, k-1, l) + sx0(1)*exg(j0+1, k-1, &
      l))*sy(-1)
      a = a + (sx0(-1)*exg(j0-1, k, l) + sx0(0)*exg(j0, k, l) + sx0(1)*exg(j0+1, k,   &
      l))*sy(0)
      a = a + (sx0(-1)*exg(j0-1, k+1, l) + sx0(0)*exg(j0, k+1, l) + sx0(1)*exg(j0+1,  &
      k+1, l))*sy(1)
      a = a + (sx0(-1)*exg(j0-1, k+2, l) + sx0(0)*exg(j0, k+2, l) + sx0(1)*exg(j0+1,  &
      k+2, l))*sy(2)
      ex(nn) = ex(nn) + a*sz(0)
      a = (sx0(-1)*exg(j0-1, k-1, l+1) + sx0(0)*exg(j0, k-1, l+1) + sx0(1)*exg(j0+1,  &
      k-1, l+1))*sy(-1)
      a = a + (sx0(-1)*exg(j0-1, k, l+1) + sx0(0)*exg(j0, k, l+1) + sx0(1)*exg(j0+1,  &
      k, l+1))*sy(0)
      a = a + (sx0(-1)*exg(j0-1, k+1, l+1) + sx0(0)*exg(j0, k+1, l+1) +               &
      sx0(1)*exg(j0+1, k+1, l+1))*sy(1)
      a = a + (sx0(-1)*exg(j0-1, k+2, l+1) + sx0(0)*exg(j0, k+2, l+1) +               &
      sx0(1)*exg(j0+1, k+2, l+1))*sy(2)
      ex(nn) = ex(nn) + a*sz(1)
      a = (sx0(-1)*exg(j0-1, k-1, l+2) + sx0(0)*exg(j0, k-1, l+2) + sx0(1)*exg(j0+1,  &
      k-1, l+2))*sy(-1)
      a = a + (sx0(-1)*exg(j0-1, k, l+2) + sx0(0)*exg(j0, k, l+2) + sx0(1)*exg(j0+1,  &
      k, l+2))*sy(0)
      a = a + (sx0(-1)*exg(j0-1, k+1, l+2) + sx0(0)*exg(j0, k+1, l+2) +               &
      sx0(1)*exg(j0+1, k+1, l+2))*sy(1)
      a = a + (sx0(-1)*exg(j0-1, k+2, l+2) + sx0(0)*exg(j0, k+2, l+2) +               &
      sx0(1)*exg(j0+1, k+2, l+2))*sy(2)
      ex(nn) = ex(nn) + a*sz(2)

      ! Compute Ey on particle
      a = (sx(-1)*eyg(j-1, k0-1, l-1) + sx(0)*eyg(j, k0-1, l-1) + sx(1)*eyg(j+1,      &
      k0-1, l-1) + sx(2)*eyg(j+2, k0-1, l-1))*sy0(-1)
      a = a + (sx(-1)*eyg(j-1, k0, l-1) + sx(0)*eyg(j, k0, l-1) + sx(1)*eyg(j+1, k0,  &
      l-1) + sx(2)*eyg(j+2, k0, l-1))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l-1) + sx(0)*eyg(j, k0+1, l-1) + sx(1)*eyg(j+1,  &
      k0+1, l-1) + sx(2)*eyg(j+2, k0+1, l-1))*sy0(1)
      ey(nn) = ey(nn) + a*sz(-1)
      a = (sx(-1)*eyg(j-1, k0-1, l) + sx(0)*eyg(j, k0-1, l) + sx(1)*eyg(j+1, k0-1, l) &
      + sx(2)*eyg(j+2, k0-1, l))*sy0(-1)
      a = a + (sx(-1)*eyg(j-1, k0, l) + sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0, l) + &
      sx(2)*eyg(j+2, k0, l))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l) + sx(0)*eyg(j, k0+1, l) + sx(1)*eyg(j+1,      &
      k0+1, l) + sx(2)*eyg(j+2, k0+1, l))*sy0(1)
      ey(nn) = ey(nn) + a*sz(0)
      a = (sx(-1)*eyg(j-1, k0-1, l+1) + sx(0)*eyg(j, k0-1, l+1) + sx(1)*eyg(j+1,      &
      k0-1, l+1) + sx(2)*eyg(j+2, k0-1, l+1))*sy0(-1)
      a = a + (sx(-1)*eyg(j-1, k0, l+1) + sx(0)*eyg(j, k0, l+1) + sx(1)*eyg(j+1, k0,  &
      l+1) + sx(2)*eyg(j+2, k0, l+1))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l+1) + sx(0)*eyg(j, k0+1, l+1) + sx(1)*eyg(j+1,  &
      k0+1, l+1) + sx(2)*eyg(j+2, k0+1, l+1))*sy0(1)
      ey(nn) = ey(nn) + a*sz(1)
      a = (sx(-1)*eyg(j-1, k0-1, l+2) + sx(0)*eyg(j, k0-1, l+2) + sx(1)*eyg(j+1,      &
      k0-1, l+2) + sx(2)*eyg(j+2, k0-1, l+2))*sy0(-1)
      a = a + (sx(-1)*eyg(j-1, k0, l+2) + sx(0)*eyg(j, k0, l+2) + sx(1)*eyg(j+1, k0,  &
      l+2) + sx(2)*eyg(j+2, k0, l+2))*sy0(0)
      a = a + (sx(-1)*eyg(j-1, k0+1, l+2) + sx(0)*eyg(j, k0+1, l+2) + sx(1)*eyg(j+1,  &
      k0+1, l+2) + sx(2)*eyg(j+2, k0+1, l+2))*sy0(1)
      ey(nn) = ey(nn) + a*sz(2)

      ! Compute Ez on particle
      a = (sx(-1)*ezg(j-1, k-1, l0-1) + sx(0)*ezg(j, k-1, l0-1) + sx(1)*ezg(j+1, k-1, &
      l0-1) + sx(2)*ezg(j+2, k-1, l0-1))*sy(-1)
      a = a + (sx(-1)*ezg(j-1, k, l0-1) + sx(0)*ezg(j, k, l0-1) + sx(1)*ezg(j+1, k,   &
      l0-1) + sx(2)*ezg(j+2, k, l0-1))*sy(0)
      a = a + (sx(-1)*ezg(j-1, k+1, l0-1) + sx(0)*ezg(j, k+1, l0-1) + sx(1)*ezg(j+1,  &
      k+1, l0-1) + sx(2)*ezg(j+2, k+1, l0-1))*sy(1)
      a = a + (sx(-1)*ezg(j-1, k+2, l0-1) + sx(0)*ezg(j, k+2, l0-1) + sx(1)*ezg(j+1,  &
      k+2, l0-1) + sx(2)*ezg(j+2, k+2, l0-1))*sy(2)
      ez(nn) = ez(nn) + a*sz0(-1)
      a = (sx(-1)*ezg(j-1, k-1, l0) + sx(0)*ezg(j, k-1, l0) + sx(1)*ezg(j+1, k-1, l0) &
      + sx(2)*ezg(j+2, k-1, l0))*sy(-1)
      a = a + (sx(-1)*ezg(j-1, k, l0) + sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k, l0) + &
      sx(2)*ezg(j+2, k, l0))*sy(0)
      a = a + (sx(-1)*ezg(j-1, k+1, l0) + sx(0)*ezg(j, k+1, l0) + sx(1)*ezg(j+1, k+1, &
      l0) + sx(2)*ezg(j+2, k+1, l0))*sy(1)
      a = a + (sx(-1)*ezg(j-1, k+2, l0) + sx(0)*ezg(j, k+2, l0) + sx(1)*ezg(j+1, k+2, &
      l0) + sx(2)*ezg(j+2, k+2, l0))*sy(2)
      ez(nn) = ez(nn) + a*sz0(0)
      a = (sx(-1)*ezg(j-1, k-1, l0+1) + sx(0)*ezg(j, k-1, l0+1) + sx(1)*ezg(j+1, k-1, &
      l0+1) + sx(2)*ezg(j+2, k-1, l0+1))*sy(-1)
      a = a + (sx(-1)*ezg(j-1, k, l0+1) + sx(0)*ezg(j, k, l0+1) + sx(1)*ezg(j+1, k,   &
      l0+1) + sx(2)*ezg(j+2, k, l0+1))*sy(0)
      a = a + (sx(-1)*ezg(j-1, k+1, l0+1) + sx(0)*ezg(j, k+1, l0+1) + sx(1)*ezg(j+1,  &
      k+1, l0+1) + sx(2)*ezg(j+2, k+1, l0+1))*sy(1)
      a = a + (sx(-1)*ezg(j-1, k+2, l0+1) + sx(0)*ezg(j, k+2, l0+1) + sx(1)*ezg(j+1,  &
      k+2, l0+1) + sx(2)*ezg(j+2, k+2, l0+1))*sy(2)
      ez(nn) = ez(nn) + a*sz0(1)

      ! Compute Bx on particle
      a = (sx(-1)*bxg(j-1, k0-1, l0-1) + sx(0)*bxg(j, k0-1, l0-1) + sx(1)*bxg(j+1,    &
      k0-1, l0-1) + sx(2)*bxg(j+2, k0-1, l0-1))*sy0(-1)
      a = a + (sx(-1)*bxg(j-1, k0, l0-1) + sx(0)*bxg(j, k0, l0-1) + sx(1)*bxg(j+1,    &
      k0, l0-1) + sx(2)*bxg(j+2, k0, l0-1))*sy0(0)
      a = a + (sx(-1)*bxg(j-1, k0+1, l0-1) + sx(0)*bxg(j, k0+1, l0-1) +               &
      sx(1)*bxg(j+1, k0+1, l0-1) + sx(2)*bxg(j+2, k0+1, l0-1))*sy0(1)
      bx(nn) = bx(nn) + a*sz0(-1)
      a = (sx(-1)*bxg(j-1, k0-1, l0) + sx(0)*bxg(j, k0-1, l0) + sx(1)*bxg(j+1, k0-1,  &
      l0) + sx(2)*bxg(j+2, k0-1, l0))*sy0(-1)
      a = a + (sx(-1)*bxg(j-1, k0, l0) + sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0,    &
      l0) + sx(2)*bxg(j+2, k0, l0))*sy0(0)
      a = a + (sx(-1)*bxg(j-1, k0+1, l0) + sx(0)*bxg(j, k0+1, l0) + sx(1)*bxg(j+1,    &
      k0+1, l0) + sx(2)*bxg(j+2, k0+1, l0))*sy0(1)
      bx(nn) = bx(nn) + a*sz0(0)
      a = (sx(-1)*bxg(j-1, k0-1, l0+1) + sx(0)*bxg(j, k0-1, l0+1) + sx(1)*bxg(j+1,    &
      k0-1, l0+1) + sx(2)*bxg(j+2, k0-1, l0+1))*sy0(-1)
      a = a + (sx(-1)*bxg(j-1, k0, l0+1) + sx(0)*bxg(j, k0, l0+1) + sx(1)*bxg(j+1,    &
      k0, l0+1) + sx(2)*bxg(j+2, k0, l0+1))*sy0(0)
      a = a + (sx(-1)*bxg(j-1, k0+1, l0+1) + sx(0)*bxg(j, k0+1, l0+1) +               &
      sx(1)*bxg(j+1, k0+1, l0+1) + sx(2)*bxg(j+2, k0+1, l0+1))*sy0(1)
      bx(nn) = bx(nn) + a*sz0(1)

      ! Compute By on particle
      a = (sx0(-1)*byg(j0-1, k-1, l0-1) + sx0(0)*byg(j0, k-1, l0-1) +                 &
      sx0(1)*byg(j0+1, k-1, l0-1))*sy(-1)
      a = a + (sx0(-1)*byg(j0-1, k, l0-1) + sx0(0)*byg(j0, k, l0-1) +                 &
      sx0(1)*byg(j0+1, k, l0-1))*sy(0)
      a = a + (sx0(-1)*byg(j0-1, k+1, l0-1) + sx0(0)*byg(j0, k+1, l0-1) +             &
      sx0(1)*byg(j0+1, k+1, l0-1))*sy(1)
      a = a + (sx0(-1)*byg(j0-1, k+2, l0-1) + sx0(0)*byg(j0, k+2, l0-1) +             &
      sx0(1)*byg(j0+1, k+2, l0-1))*sy(2)
      by(nn) = by(nn) + a*sz0(-1)
      a = (sx0(-1)*byg(j0-1, k-1, l0) + sx0(0)*byg(j0, k-1, l0) + sx0(1)*byg(j0+1,    &
      k-1, l0))*sy(-1)
      a = a + (sx0(-1)*byg(j0-1, k, l0) + sx0(0)*byg(j0, k, l0) + sx0(1)*byg(j0+1, k, &
      l0))*sy(0)
      a = a + (sx0(-1)*byg(j0-1, k+1, l0) + sx0(0)*byg(j0, k+1, l0) +                 &
      sx0(1)*byg(j0+1, k+1, l0))*sy(1)
      a = a + (sx0(-1)*byg(j0-1, k+2, l0) + sx0(0)*byg(j0, k+2, l0) +                 &
      sx0(1)*byg(j0+1, k+2, l0))*sy(2)
      by(nn) = by(nn) + a*sz0(0)
      a = (sx0(-1)*byg(j0-1, k-1, l0+1) + sx0(0)*byg(j0, k-1, l0+1) +                 &
      sx0(1)*byg(j0+1, k-1, l0+1))*sy(-1)
      a = a + (sx0(-1)*byg(j0-1, k, l0+1) + sx0(0)*byg(j0, k, l0+1) +                 &
      sx0(1)*byg(j0+1, k, l0+1))*sy(0)
      a = a + (sx0(-1)*byg(j0-1, k+1, l0+1) + sx0(0)*byg(j0, k+1, l0+1) +             &
      sx0(1)*byg(j0+1, k+1, l0+1))*sy(1)
      a = a + (sx0(-1)*byg(j0-1, k+2, l0+1) + sx0(0)*byg(j0, k+2, l0+1) +             &
      sx0(1)*byg(j0+1, k+2, l0+1))*sy(2)
      by(nn) = by(nn) + a*sz0(1)

      ! Compute Bz on particle
      a = (sx0(-1)*bzg(j0-1, k0-1, l-1) + sx0(0)*bzg(j0, k0-1, l-1) +                 &
      sx0(1)*bzg(j0+1, k0-1, l-1))*sy0(-1)
      a = a + (sx0(-1)*bzg(j0-1, k0, l-1) + sx0(0)*bzg(j0, k0, l-1) +                 &
      sx0(1)*bzg(j0+1, k0, l-1))*sy0(0)
      a = a + (sx0(-1)*bzg(j0-1, k0+1, l-1) + sx0(0)*bzg(j0, k0+1, l-1) +             &
      sx0(1)*bzg(j0+1, k0+1, l-1))*sy0(1)
      bz(nn) = bz(nn) + a*sz(-1)
      a = (sx0(-1)*bzg(j0-1, k0-1, l) + sx0(0)*bzg(j0, k0-1, l) + sx0(1)*bzg(j0+1,    &
      k0-1, l))*sy0(-1)
      a = a + (sx0(-1)*bzg(j0-1, k0, l) + sx0(0)*bzg(j0, k0, l) + sx0(1)*bzg(j0+1,    &
      k0, l))*sy0(0)
      a = a + (sx0(-1)*bzg(j0-1, k0+1, l) + sx0(0)*bzg(j0, k0+1, l) +                 &
      sx0(1)*bzg(j0+1, k0+1, l))*sy0(1)
      bz(nn) = bz(nn) + a*sz(0)
      a = (sx0(-1)*bzg(j0-1, k0-1, l+1) + sx0(0)*bzg(j0, k0-1, l+1) +                 &
      sx0(1)*bzg(j0+1, k0-1, l+1))*sy0(-1)
      a = a + (sx0(-1)*bzg(j0-1, k0, l+1) + sx0(0)*bzg(j0, k0, l+1) +                 &
      sx0(1)*bzg(j0+1, k0, l+1))*sy0(0)
      a = a + (sx0(-1)*bzg(j0-1, k0+1, l+1) + sx0(0)*bzg(j0, k0+1, l+1) +             &
      sx0(1)*bzg(j0+1, k0+1, l+1))*sy0(1)
      bz(nn) = bz(nn) + a*sz(1)
      a = (sx0(-1)*bzg(j0-1, k0-1, l+2) + sx0(0)*bzg(j0, k0-1, l+2) +                 &
      sx0(1)*bzg(j0+1, k0-1, l+2))*sy0(-1)
      a = a + (sx0(-1)*bzg(j0-1, k0, l+2) + sx0(0)*bzg(j0, k0, l+2) +                 &
      sx0(1)*bzg(j0+1, k0, l+2))*sy0(0)
      a = a + (sx0(-1)*bzg(j0-1, k0+1, l+2) + sx0(0)*bzg(j0, k0+1, l+2) +             &
      sx0(1)*bzg(j0+1, k0+1, l+2))*sy0(1)
      bz(nn) = bz(nn) + a*sz(2)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! __________________________________________________________________________
    ! Particle pusher

    SELECT CASE (particle_pusher)
      !! Vay pusher -- Full push
    CASE (1_idp)
      CALL pxr_ebcancelpush3d(blocksize, uxp(ip:ip+blocksize-1),                      &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt,  &
      0_idp)

      !! Boris pusher -- Full push
    CASE DEFAULT

      !! Push momentum using the Boris method in a single subroutine

      CALL pxr_boris_push_u_3d(blocksize, uxp(ip:ip+blocksize-1),                     &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize, uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),   &
      uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1))
    END SELECT

    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize, xp(ip:ip+blocksize-1), yp(ip:ip+blocksize-1),         &
    zp(ip:ip+blocksize-1), uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),            &
    uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1), dt)
  ENDDO

  RETURN
END SUBROUTINE field_gathering_plus_particle_pusher_3_3_3

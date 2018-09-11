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
  USE fields
  USE shared_data
  USE params
  USE particles
  USE time_stat
  IMPLICIT NONE
  
#if defined(DEBUG)
  WRITE(0, *) "Field gathering + Push_particles: start"
#endif
  IF (nspecies .EQ. 0_idp) RETURN
    ! ___________________________________________________________
    ! 3D CASE
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
  USE particles
  USE constants
  USE tiling
  USE time_stat
  ! Vtune/SDE profiling
#if defined(PROFILING) && PROFILING==3
  USE ITT_SDE_FORTRAN
#endif
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

            !! Boris pusher -- Full push
              !! Push momentum using the Boris method in a single subroutine
              CALL pxr_boris_push_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,   &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)

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
  USE particles
  USE constants
  USE tiling
  USE time_stat
  USE params
  ! Vtune/SDE profiling
#if defined(PROFILING) && PROFILING==3
  USE ITT_SDE_FORTRAN
#endif
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
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            IF (count .EQ. 0) CYCLE
            curr_tile%part_ex(1:count)=0.0_num
            curr_tile%part_ey(1:count)=0.0_num
            curr_tile%part_ez(1:count)=0.0_num
            curr_tile%part_bx(1:count)=0.0_num
            curr_tile%part_by(1:count)=0.0_num
            curr_tile%part_bz(1:count)=0.0_num

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
  USE particles
  USE constants
  USE tiling
  USE time_stat
  ! Vtune/SDE profiling
#if defined(PROFILING) && PROFILING==3
  USE ITT_SDE_FORTRAN
#endif
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

#if defined(DEBUG)
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
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            IF (count .EQ. 0) CYCLE
              !! Push momentum using the Boris method in a single subroutine
              CALL pxr_boris_push_u_3d(count, curr_tile%part_ux, curr_tile%part_uy,   &
              curr_tile%part_uz, curr_tile%part_gaminv, curr_tile%part_ex,            &
              curr_tile%part_ey, curr_tile%part_ez, curr_tile%part_bx,                &
              curr_tile%part_by, curr_tile%part_bz, curr%charge, curr%mass, dtt)

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
  USE omp_lib
  USE constants
  USE params
  USE particles

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

      !! Push momentum using the Boris method in a single subroutine

      CALL pxr_boris_push_u_3d(blocksize, uxp(ip:ip+blocksize-1),                     &
      uyp(ip:ip+blocksize-1), uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1),      &
      ex(ip:ip+blocksize-1), ey(ip:ip+blocksize-1), ez(ip:ip+blocksize-1),            &
      bx(ip:ip+blocksize-1), by(ip:ip+blocksize-1), bz(ip:ip+blocksize-1), q, m, dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize, uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),   &
      uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1))

    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize, xp(ip:ip+blocksize-1), yp(ip:ip+blocksize-1),         &
    zp(ip:ip+blocksize-1), uxp(ip:ip+blocksize-1), uyp(ip:ip+blocksize-1),            &
    uzp(ip:ip+blocksize-1), gaminv(ip:ip+blocksize-1), dt)
  ENDDO

  RETURN
END SUBROUTINE field_gathering_plus_particle_pusher_1_1_1

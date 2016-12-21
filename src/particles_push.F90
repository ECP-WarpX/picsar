! ________________________________________________________________________________________
! PARTICLES_PUSH.F90
!
! Subroutines for the particle pusher in 3D
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

SUBROUTINE field_gathering_plus_particle_pusher
! ________________________________________________________________________________________
  USE fields
  USE shared_data
  USE params
  USE particles
  USE time_stat
  IMPLICIT NONE

#if defined(DEBUG)
  WRITE(0,*) "Field gathering + Push_particles: start"
#endif
  IF (nspecies .EQ. 0_idp) RETURN
  SELECT CASE (c_dim)
  
    ! ___________________________________________________________
    ! 2D CASE X Z
    CASE (2)

      ! Particle advance (one time step)
      CALL field_gathering_plus_particle_pusher_sub_2d(ex,ey,ez,bx,by,bz,nx,ny,nz, &
                                                       nxguards,nyguards,nzguards, &
                                                       nxjguards,nyjguards,nzjguards,&
                                                       nox,noy,noz,dx,dy,dz,dt)

    ! ___________________________________________________________
    ! 3D CASE
    CASE DEFAULT

      ! The field gathering and the particle pusher are performed together
      IF (fg_p_pp_separated.eq.0) THEN

        CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,&
        nxguards,nyguards,nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

      ELSE IF (fg_p_pp_separated.eq.1) THEN

        CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
         nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

      ! The field gathering and the particle pusher are performed separately
      ELSE

        CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
         nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

        CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
         nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

      ENDIF

  END SELECT

#if defined(DEBUG)
  WRITE(0,*) "Field gathering + Push_particles: stop"
#endif

END SUBROUTINE field_gathering_plus_particle_pusher

! ________________________________________________________________________________________
!> @brief
!> Particle pusher in 3D called by the main function push_particle

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2015
!> Revision 10.06.2015
!
!> @param[in] exg,eyg,ezg electric field grids
!> @param[in] bxg,byg,bzg electric field grids
!> @param[in] nxx,nyy,nzz number of cells in each direction for the grids
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction for the grids
!> @param[in] nxjguard,nyjguard,nzjguard number of guard cells for the current grids
!> @param[in] noxx,noyy,nozz interpolation orders
!> @param[in] dxx,dyy,dzz space steps
!> @param[in] dtt time step
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!>
SUBROUTINE field_gathering_plus_particle_pusher_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
           nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,  &
           nozz,dxx,dyy,dzz,dtt,l_lower_order_in_v_in)
! ________________________________________________________________________________________

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
  INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
  LOGICAL(lp)                   :: l_lower_order_in_v_in
  REAL(num), INTENT(IN)    :: exg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: eyg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: ezg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bxg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: byg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bzg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: dxx,dyy,dzz, dtt
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(grid_tile), POINTER        :: currg
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                :: tdeb, tend
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin,ipmax, ip
  INTEGER(idp)             :: nxjg,nyjg,nzjg
  LOGICAL(lp)              :: isgathered=.FALSE.

  tdeb=MPI_WTIME()

#if PROFILING==3
  CALL start_collection()
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
  !$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,&
  !$OMP bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz,c_dim,l_lower_order_in_v_in,&
  !$OMP particle_pusher,fieldgathe,LVEC_fieldgathe) &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile, currg, count,jmin,jmax,kmin,kmax,lmin, &
  !$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,nxjg,nyjg,nzjg, isgathered)
  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr=>species_parray(1)
        curr_tile=>curr%array_of_tiles(ix,iy,iz)
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

        DO ispecies=1, nspecies ! LOOP ON SPECIES
          curr=>species_parray(ispecies)
          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isgathered=.TRUE.
        END DO

        IF (isgathered) THEN
          currg=>aofgrid_tiles(ix,iy,iz)
          currg%extile=exg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%eytile=eyg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%eztile=ezg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bxtile=bxg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bytile=byg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bztile=bzg(jmin:jmax,kmin:kmax,lmin:lmax)
          DO ispecies=1, nspecies ! LOOP ON SPECIES
            ! - Get current tile properties
            ! - Init current tile variables
            curr=>species_parray(ispecies)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
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
            CASE (2) ! 2D CASE X Z
              !!! --- Gather electric and magnetic fields on particles
              CALL geteb2dxz_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,            &
                          curr_tile%part_z, curr_tile%part_ex,                                     &
                          curr_tile%part_ey,curr_tile%part_ez,                                     &
                          curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz,                  &
                          curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                     &
                          curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,          &
                          curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,               &
                          nzjg,noxx,noyy,nozz,currg%extile,currg%eytile,                           &
                          currg%eztile,                                                            &
                          currg%bxtile,currg%bytile,currg%bztile,                                  &
                          .FALSE.,                                                                 &
                          l_lower_order_in_v_in,                                                   &
                          LVEC_fieldgathe,                                                         &
                          fieldgathe)
                          
            CASE DEFAULT ! 3D CASE

              !!! --- Gather electric and magnetic fields on particles
              CALL geteb3d_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,            &
                          curr_tile%part_z, curr_tile%part_ex,                                   &
                          curr_tile%part_ey,curr_tile%part_ez,                                   &
                          curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz,                &
                          curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                   &
                          curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,        &
                          curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,             &
                          nzjg,noxx,noyy,nozz,currg%extile,currg%eytile,                         &
                          currg%eztile,                                                          &
                          currg%bxtile,currg%bytile,currg%bztile                                 &
                          ,.FALSE.,l_lower_order_in_v_in,                                        &
                          LVEC_fieldgathe,                                                       &
                          fieldgathe)
            END SELECT

            SELECT CASE (particle_pusher)
            !! Vay pusher -- Full push
            CASE (1_idp)
              CALL pxr_ebcancelpush3d(count,curr_tile%part_ux, curr_tile%part_uy,&
              curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_ex,        &
              curr_tile%part_ey,                                                  &
              curr_tile%part_ez,curr_tile%part_bx, curr_tile%part_by,            &
              curr_tile%part_bz,curr%charge,curr%mass,dtt,0_idp)
            !! Boris pusher -- Full push
            CASE DEFAULT
              !! --- Push velocity with E half step
              CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,      &
              curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,           &
              curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
              !! --- Set gamma of particles
              CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,    &
              curr_tile%part_uz, curr_tile%part_gaminv)
              !! --- Push velocity with B half step
              CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,       &
              curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx,        &
              curr_tile%part_by,                                                 &
              curr_tile%part_bz, curr%charge,curr%mass,dtt)
              !!! --- Push velocity with E half step
              CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,       &
              curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,           &
              curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
              !! --- Set gamma of particles
              CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,     &
              curr_tile%part_uz, curr_tile%part_gaminv)
            END SELECT
            !!!! --- push particle species positions a time step
            CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y,          &
            curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy,             &
            curr_tile%part_uz,curr_tile%part_gaminv,dtt)
          END DO! END LOOP ON SPECIES
        ENDIF
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO


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
!> Particle pusher in 3D called by the main function push_particle for the subroutines with cache blocking
!
!> @author
!> Mathieu Lobet

!> @date
!> Creation 2015
!> Revision 10.06.2015
!
!>
!> @param[in] exg,eyg,ezg electric field grids
!> @param[in] bxg,byg,bzg electric field grids
!> @param[in] nxx,nyy,nzz number of cells in each direction for the grids
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction for the grids
!> @param[in] nxjguard,nyjguard,nzjguard number of guard cells for the current grids
!> @param[in] noxx,noyy,nozz interpolation orders
!> @param[in] dxx,dyy,dzz space steps
!> @param[in] dtt time step
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!>
SUBROUTINE field_gathering_plus_particle_pusher_cacheblock_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
      nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt,l_lower_order_in_v_in)
! ________________________________________________________________________________________

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
  INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
  LOGICAL(lp)                   :: l_lower_order_in_v_in
  REAL(num), INTENT(IN)    :: exg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: eyg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: ezg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bxg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: byg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bzg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: dxx,dyy,dzz, dtt
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(grid_tile), POINTER        :: currg
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                :: tdeb, tend
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin,ipmax, ip
  INTEGER(idp)             :: nxjg,nyjg,nzjg
  LOGICAL(lp)              :: isgathered=.FALSE.

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==3
  CALL start_vtune_collection()
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
  !$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,&
  !$OMP bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz,c_dim,lvec_fieldgathe,l_lower_order_in_v) &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile, currg, count,jmin,jmax,kmin,kmax,lmin, &
  !$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,nxjg,nyjg,nzjg, isgathered)
  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr=>species_parray(1)
        curr_tile=>curr%array_of_tiles(ix,iy,iz)
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

        DO ispecies=1, nspecies ! LOOP ON SPECIES
          curr=>species_parray(ispecies)
          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isgathered=.TRUE.
        END DO

        IF (isgathered) THEN
            currg=>aofgrid_tiles(ix,iy,iz)
          currg%extile=exg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%eytile=eyg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%eztile=ezg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bxtile=bxg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bytile=byg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bztile=bzg(jmin:jmax,kmin:kmax,lmin:lmax)
          DO ispecies=1, nspecies ! LOOP ON SPECIES
            ! - Get current tile properties
            ! - Init current tile variables
            curr=>species_parray(ispecies)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
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
               CALL field_gathering_plus_particle_pusher_1_1_1(count,   &
                    curr_tile%part_x,curr_tile%part_y,curr_tile%part_z, &
                    curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz, &
                    curr_tile%part_gaminv, &
                    curr_tile%part_ex,curr_tile%part_ey,curr_tile%part_ez, &
                    curr_tile%part_bx,curr_tile%part_by,curr_tile%part_bz, &
                    curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,curr_tile%z_grid_tile_min, &
                    dxx,dyy,dzz,dtt,&
                    curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,&
                    nxjg,nyjg,nzjg, &
                    currg%extile,currg%eytile,currg%eztile,&
                    currg%bxtile,currg%bytile,currg%bztile,&
                    curr%charge,curr%mass,lvec_fieldgathe,l_lower_order_in_v)

            ELSE IF ((noxx.eq.2).and.(noyy.eq.2).and.(nozz.eq.2)) THEN

              !!! ---- Loop by blocks over particles in a tile (blocking)
               CALL field_gathering_plus_particle_pusher_2_2_2(count,   &
                    curr_tile%part_x,curr_tile%part_y,curr_tile%part_z, &
                    curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz, &
                    curr_tile%part_gaminv, &
                    curr_tile%part_ex,curr_tile%part_ey,curr_tile%part_ez, &
                    curr_tile%part_bx,curr_tile%part_by,curr_tile%part_bz, &
                    curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,curr_tile%z_grid_tile_min, &
                    dxx,dyy,dzz,dtt,&
                    curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,&
                    nxjg,nyjg,nzjg, &
                    currg%extile,currg%eytile,currg%eztile,&
                    currg%bxtile,currg%bytile,currg%bztile,&
                    curr%charge,curr%mass,lvec_fieldgathe,l_lower_order_in_v)

            ELSE IF ((noxx.eq.3).and.(noyy.eq.3).and.(nozz.eq.3)) THEN

              !!! ---- Loop by blocks over particles in a tile (blocking)
               CALL field_gathering_plus_particle_pusher_3_3_3(count,   &
                    curr_tile%part_x,curr_tile%part_y,curr_tile%part_z, &
                    curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz, &
                    curr_tile%part_gaminv, &
                    curr_tile%part_ex,curr_tile%part_ey,curr_tile%part_ez, &
                    curr_tile%part_bx,curr_tile%part_by,curr_tile%part_bz, &
                    curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,curr_tile%z_grid_tile_min, &
                    dxx,dyy,dzz,dtt,&
                    curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,&
                    nxjg,nyjg,nzjg, &
                    currg%extile,currg%eytile,currg%eztile,&
                    currg%bxtile,currg%bytile,currg%bztile,&
                    curr%charge,curr%mass,lvec_fieldgathe,l_lower_order_in_v)

            ENDIF

          END DO! END LOOP ON SPECIES
        ENDIF
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO

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
!> @param[in] exg,eyg,ezg electric field grids
!> @param[in] bxg,byg,bzg electric field grids
!> @param[in] nxx,nyy,nzz number of cells in each direction for the grids
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction for the grids
!> @param[in] nxjguard,nyjguard,nzjguard number of guard cells for the current grids
!> @param[in] noxx,noyy,nozz interpolation orders
!> @param[in] dxx,dyy,dzz space steps
!> @param[in] dtt time step
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!
SUBROUTINE particle_pusher_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
      nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,&
      noxx,noyy,nozz,dxx,dyy,dzz,dtt,l_lower_order_in_v_in)
! ________________________________________________________________________________________
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
  INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
  LOGICAL(lp)                   :: l_lower_order_in_v_in
  REAL(num), INTENT(IN)    :: exg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: eyg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: ezg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bxg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: byg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: bzg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN)    :: dxx,dyy,dzz, dtt
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(grid_tile), POINTER        :: currg
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                :: tdeb, tend
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin,ipmax, ip
  INTEGER(idp)             :: nxjg,nyjg,nzjg
  LOGICAL(lp)              :: isgathered=.FALSE.

  IF (nspecies .EQ. 0_idp) RETURN

  tdeb=MPI_WTIME()

#if VTUNE==3
  CALL start_vtune_collection()
#endif

#if defined(DEBUG)
#endif


  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
  !$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,bxg,byg, &
  !$OMP bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz,c_dim, particle_pusher) &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile, currg, count,jmin,jmax,kmin,kmax,lmin, &
  !$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,nxjg,nyjg,nzjg, isgathered)
  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr=>species_parray(1)
        curr_tile=>curr%array_of_tiles(ix,iy,iz)
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
        DO ispecies=1, nspecies ! LOOP ON SPECIES
          curr=>species_parray(ispecies)
          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isgathered=.TRUE.
        END DO
        IF (isgathered) THEN
            currg=>aofgrid_tiles(ix,iy,iz)
            currg%extile=exg(jmin:jmax,kmin:kmax,lmin:lmax)
            currg%eytile=eyg(jmin:jmax,kmin:kmax,lmin:lmax)
            currg%eztile=ezg(jmin:jmax,kmin:kmax,lmin:lmax)
            currg%bxtile=bxg(jmin:jmax,kmin:kmax,lmin:lmax)
            currg%bytile=byg(jmin:jmax,kmin:kmax,lmin:lmax)
            currg%bztile=bzg(jmin:jmax,kmin:kmax,lmin:lmax)
            DO ispecies=1, nspecies ! LOOP ON SPECIES
              ! - Get current tile properties
              ! - Init current tile variables
              curr=>species_parray(ispecies)
              curr_tile=>curr%array_of_tiles(ix,iy,iz)
              count=curr_tile%np_tile(1)
              IF (count .EQ. 0) CYCLE
              SELECT CASE (particle_pusher)
              !! Vay pusher -- Full push
              CASE (1_idp)
                  CALL pxr_ebcancelpush3d(count,curr_tile%part_ux, curr_tile%part_uy,&
                  curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_ex,        &
                  curr_tile%part_ey,                                                  &
                  curr_tile%part_ez,curr_tile%part_bx, curr_tile%part_by,            &
                  curr_tile%part_bz,curr%charge,curr%mass,dtt,0_idp)
              !! Boris pusher -- Full push
              CASE DEFAULT
                !! --- Push velocity with E half step
                CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,          &
                curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,               &
                curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
                !! --- Set gamma of particles
                CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,        &
                curr_tile%part_uz, curr_tile%part_gaminv)
                !! --- Push velocity with B half step
                CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,          &
                curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx,           &
                curr_tile%part_by,                                                    &
                curr_tile%part_bz, curr%charge,curr%mass,dtt)
                !!! --- Push velocity with E half step
                CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,          &
                curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,              &
                curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
                !! --- Set gamma of particles
                CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,        &
                curr_tile%part_uz, curr_tile%part_gaminv)
              END SELECT
              !!!! --- push particle species positions a time step
              CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y,             &
              curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy,                &
              curr_tile%part_uz,curr_tile%part_gaminv,dtt)
            END DO! END LOOP ON SPECIES
          ENDIF
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO


#if VTUNE==3
  CALL stop_vtune_collection()
#endif

  tend=MPI_WTIME()
  pushtime=pushtime+(tend-tdeb)
  IF (it.ge.timestat_itstart) THEN
    localtimes(1) = localtimes(1) + (tend-tdeb)
  ENDIF

#if defined(DEBUG)
  WRITE(0,*) "Push_particles: stop"
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
SUBROUTINE pxrpush_particles_part1
! ________________________________________________________________________________________
  USE fields
  USE shared_data
  USE params
  IMPLICIT NONE

#if defined(DEBUG)
  WRITE(0,*) "pxrpush_particles_part1: start"
#endif

  CALL pxrpush_particles_part1_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l4symtry,l_lower_order_in_v)

#if defined(DEBUG)
  WRITE(0,*) "pxrpush_particles_part1: stop"
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
!> @param[in] exg,eyg,ezg electric field grids
!> @param[in] bxg,byg,bzg electric field grids
!> @param[in] nxx,nyy,nzz number of cells in each direction for the grids
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction for the grids
!> @param[in] nxjguard,nyjguard,nzjguard number of guard cells for the current grids
!> @param[in] noxx,noyy,nozz interpolation orders
!> @param[in] dxx,dyy,dzz space steps
!> @param[in] dtt time step
!> @param[in] l4symtry_in 
!> @param[in] l_lower_order_in_v_in flag to activate interpolation at a lower order
!>
SUBROUTINE pxrpush_particles_part1_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
      nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy, &
      dzz,dtt,l4symtry_in, l_lower_order_in_v_in)
! ________________________________________________________________________________________

  USE particles
  USE constants
  USE tiling
  USE timing
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
  LOGICAL(lp) , INTENT(IN) :: l4symtry_in, l_lower_order_in_v_in
  REAL(num), INTENT(IN) :: exg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: eyg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: ezg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: bxg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: byg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: bzg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
  REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
  INTEGER(idp)          :: ispecies, ix, iy, iz, count
  INTEGER(idp)          :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nxc, nyc, nzc
  INTEGER(idp)                    :: nxjg,nyjg,nzjg
  LOGICAL(lp)                     :: isgathered=.FALSE.


  IF (nspecies .EQ. 0_idp) RETURN
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
  !$OMP nxjguard,nyjguard,nzjguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy, &
  !$OMP l4symtry_in, l_lower_order_in_v_in, nozz,c_dim,fieldgathe,particle_pusher) &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,currg,count,jmin,jmax,kmin,kmax,lmin, &
  !$OMP lmax,nxc,nyc,nzc,nxjg,nyjg,nzjg,isgathered)
  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr=>species_parray(1)
        curr_tile=>curr%array_of_tiles(ix,iy,iz)
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
        DO ispecies=1, nspecies ! LOOP ON SPECIES
          curr=>species_parray(ispecies)
          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          count=curr_tile%np_tile(1)
          IF (count .GT. 0) isgathered=.TRUE.
        END DO
        IF (isgathered) THEN
          currg=>aofgrid_tiles(ix,iy,iz)
          currg%extile=exg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%eytile=eyg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%eztile=ezg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bxtile=bxg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bytile=byg(jmin:jmax,kmin:kmax,lmin:lmax)
          currg%bztile=bzg(jmin:jmax,kmin:kmax,lmin:lmax)
          DO ispecies=1, nspecies ! LOOP ON SPECIES
            ! - Get current tile properties
            ! - Init current tile variables
            curr=>species_parray(ispecies)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
            count=curr_tile%np_tile(1)
            IF (count .EQ. 0) CYCLE

            curr_tile%part_ex(1:count) = 0.0_num
            curr_tile%part_ey(1:count) = 0.0_num
            curr_tile%part_ez(1:count) = 0.0_num
            curr_tile%part_bx(1:count)=0.0_num
            curr_tile%part_by(1:count)=0.0_num
            curr_tile%part_bz(1:count)=0.0_num

            !!! ---- Loop by blocks over particles in a tile (blocking)
            SELECT CASE (c_dim)
            CASE (2) ! 2D CASE
              !!! --- Gather electric and magnetic fields on particles
              CALL geteb2dxz_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,           &
                    curr_tile%part_z, curr_tile%part_ex,                              &
                    curr_tile%part_ey,curr_tile%part_ez,                         &
                    curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz,       &
                    curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,              &
                    curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,   &
                    curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,        &
                    nzjg,noxx,noyy,nozz,currg%extile,currg%eytile,           &
                    currg%eztile,                                                &
                    currg%bxtile,currg%bytile,currg%bztile                       &
                    ,l4symtry_in,l_lower_order_in_v_in)
            CASE DEFAULT ! 3D CASE
              !!! --- Gather electric and magnetic fields on particles
              CALL geteb3d_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,                &
                    curr_tile%part_z, curr_tile%part_ex,                               &
                    curr_tile%part_ey,curr_tile%part_ez,                              &
                    curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz,              &
                    curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                 &
                    curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile, &
                    curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,           &
                    nzjg,noxx,noyy,nozz,currg%extile,currg%eytile,                &
                    currg%eztile,                                                  &
                    currg%bxtile,currg%bytile,currg%bztile                         &
                    ,l4symtry_in,l_lower_order_in_v_in,fieldgathe)
            END SELECT
            SELECT CASE (particle_pusher)
            !! Vay pusher -- half push part 1

          CASE (1_idp)
              CALL pxr_ebcancelpush3d(count,curr_tile%part_ux, curr_tile%part_uy,&
              curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_ex,        &
              curr_tile%part_ey,                                                  &
              curr_tile%part_ez,curr_tile%part_bx, curr_tile%part_by,            &
              curr_tile%part_bz,curr%charge,curr%mass,dtt,1_idp)
              !! Boris pusher -- half push part 1
            CASE DEFAULT
              !! --- Push velocity with E half step
              CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                    &
              curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,               &
              curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
              !! --- Set gamma of particles
              CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
              curr_tile%part_uz, curr_tile%part_gaminv)
              !! --- Push velocity with B half step
              CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                   &
              curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by,  &
              curr_tile%part_bz, curr%charge,curr%mass,dtt*0.5_num)
            END SELECT
          END DO! END LOOP ON SPECIES
        ENDIF
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO

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
SUBROUTINE pxrpush_particles_part2
! ________________________________________________________________________________________
USE particles
USE constants
USE fields
USE params
USE shared_data
USE tiling
IMPLICIT NONE
INTEGER(idp) :: ispecies, ix, iy, iz, count
INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER    :: curr_tile
REAL(num)                       :: tdeb, tend
INTEGER(idp)                    :: nxc, nyc, nzc

#if defined(DEBUG)
  WRITE(0,*) "pxrpush_particles_part2: start"
#endif

IF (nspecies .EQ. 0_idp) RETURN
tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$OMP nxjguards,nyjguards,nzjguards,ex,ey,ez,bx,by,bz,dx,dy,dz,dt,c_dim, particle_pusher) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,nxc,nyc,nzc)
DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
        DO ix=1, ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                ! - Get current tile properties
                ! - Init current tile variables
                curr=>species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .EQ. 0) CYCLE
                SELECT CASE (particle_pusher)
                !! Vay pusher -- half push part 2
              CASE (1_idp)
                  CALL pxr_ebcancelpush3d(count,curr_tile%part_ux, curr_tile%part_uy,&
                  curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_ex,        &
                  curr_tile%part_ey,                                                  &
                  curr_tile%part_ez,curr_tile%part_bx, curr_tile%part_by,            &
                  curr_tile%part_bz,curr%charge,curr%mass,dt,2_idp)
                CASE DEFAULT
                  !! Boris pusher -- half push part 2
                  !!! --- Push velocity with B half step
                  CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                      &
                  curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by,     &
                  curr_tile%part_bz, curr%charge,curr%mass,dt*0.5_num)
                  !! --- Push velocity with E half step
                  CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                      &
                  curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,                  &
                  curr_tile%part_ez, curr%charge,curr%mass,dt*0.5_num)
                  !! --- Sets gamma of particles
                  CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                    &
                  curr_tile%part_uz, curr_tile%part_gaminv)
                END SELECT

                SELECT CASE (c_dim)
                CASE (2) ! 2D CASE
                  !! --- Advance particle position of one time step
                  CALL pxr_pushxz(count,curr_tile%part_x,                                  &
                  curr_tile%part_z, curr_tile%part_ux,                               &
                  curr_tile%part_uz,curr_tile%part_gaminv,dt)
                CASE DEFAULT ! 3D CASE
                  !! --- Advance particle position of one time step
                  CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y,                     &
                  curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy,                 &
                  curr_tile%part_uz,curr_tile%part_gaminv,dt)
                END SELECT
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)

#if defined(DEBUG)
  WRITE(0,*) "pxrpush_particles_part2: stop"
#endif

END SUBROUTINE pxrpush_particles_part2



! ________________________________________________________________________________________
!> @brief
!> Advance particle positions
!
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] ux,uy,uz normalized momentum in each direction
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] dt time step
SUBROUTINE pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
! ________________________________________________________________________________________
  USE constants
  USE omp_lib
  
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num)      :: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np), gaminv(np)
  REAL(num)      :: dt
  INTEGER(idp)   :: ip

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
  DO ip=1,np
    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

  RETURN
END SUBROUTINE pxr_pushxyz


! ________________________________________________________________________________________
!> @brief
!> Push the particle velocity with E field
!
!> @details
!>fast b-field rotation algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex,ey,ez particle electric fields in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
SUBROUTINE pxr_epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
! ________________________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp) :: np
  REAL(num)    :: uxp(np),uyp(np),uzp(np)
  REAL(num)    :: ex(np),ey(np),ez(np)
  REAL(num)    :: q,m,dt
  INTEGER(idp) :: ip
  REAL(num)    :: const

  const = q*dt/m

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,ex,ey,ez)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
  !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
  !DIR$ SIMD
#endif
  DO ip=1,np
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif
  RETURN
END SUBROUTINE pxr_epush_v

! ________________________________________________________________________________________
!> @brief
!> Push the particle velocity with B field (Boris algorithm)
!
!> @details
!>fast b-field rotation algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @gparam[in] gaminv particle Lorentz factors
!> @param[in] bx,by,bz particle magnetic fields in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
SUBROUTINE pxr_bpush_v(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,dt)
! ________________________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num)      :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num)      :: bx(np), by(np), bz(np)
  REAL(num)      :: q,m,dt
  INTEGER(idp)   :: ip
  REAL(num)      :: const,sx,sy,sz,tx,ty,tz,tsqi,uxppr,uyppr,uzppr

  const = q*dt*0.5_num/m

#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
    !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
    !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64,uxp,uyp,uzp)
    !IBM* ALIGN(64,bx,by,bz)
    !IBM* ALIGN(64,gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
  DO ip=1,np
    tx = gaminv(ip)*bx(ip)*const
    ty = gaminv(ip)*by(ip)*const
    tz = gaminv(ip)*bz(ip)*const
    tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
    sx = tx*tsqi
    sy = ty*tsqi
    sz = tz*tsqi
    uxppr = uxp(ip) + uyp(ip)*tz - uzp(ip)*ty
    uyppr = uyp(ip) + uzp(ip)*tx - uxp(ip)*tz
    uzppr = uzp(ip) + uxp(ip)*ty - uyp(ip)*tx
    uxp(ip) = uxp(ip) + uyppr*sz - uzppr*sy
    uyp(ip) = uyp(ip) + uzppr*sx - uxppr*sz
    uzp(ip) = uzp(ip) + uxppr*sy - uyppr*sx
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

  RETURN

END SUBROUTINE pxr_bpush_v


! ________________________________________________________________________________________
!> @brief
!>  Push the particle velocity with B field (Boris algorithm)
!
!> @details
!> fast b-field rotation algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!
SUBROUTINE pxr_set_gamma(np,uxp,uyp,uzp,gaminv)
! ________________________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp)   :: ip, np
  REAL(num) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num) :: clghtisq, usq

  clghtisq = 1.0_num/clight**2

#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
    !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64,uxp,uyp,uzp)
    !IBM* ALIGN(64,gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif

  DO ip=1,np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
  END DO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

RETURN

END SUBROUTINE pxr_set_gamma


! ________________________________________________________________________________________
!> @brief
!> Push the particle velocity with E and B fields, assuming Vmid = 0.5*(Vold+Vnew),
!> solving directly for the new gamma.
!>
!> @details
!> This offers better cancellation of E+VxB than the Boris velocity push.
!> Question: should we recompute gamma from the new u, in order to prevent roundoff errors
!> to create mismatched values of u and gamma?
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
!
!>
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gi
!> @param[in] exp,eyp,ezp particle electric field values in each direction
!> @param[in] bxp,byp,bzp particle electric field values in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
!> @param[in] which algorithm
!
SUBROUTINE pxr_ebcancelpush3d(np,uxp,uyp,uzp,gi,exp,eyp,ezp,bxp,byp,bzp,q,m,dt,which)
! ________________________________________________________________________________________

  USE constants
  
  INTEGER(idp) :: np,which
  REAL(num)    :: uxp(np),uyp(np),uzp(np),gi(np)
  REAL(num)    :: exp(np),eyp(np),ezp(np),bxp(np),byp(np),bzp(np)
  REAL(num)    :: q,m,dt
  INTEGER(idp) :: ip
  REAL(num)    :: const,bconst,s,gisq,invclight,invclightsq,gprsq
  REAL(num)    :: tx,ty,tz,tu,uxpr,uypr,uzpr,bg,vx,vy,vz
  REAL(num)    :: taux,tauy,tauz,tausq,ust,sigma

  invclight   = 1./clight
  invclightsq = 1./(clight*clight)

  IF (which==0) THEN
    !     --- full push
    const = q*dt/m
    bconst = 0.5_num*const


#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
    DO ip=1,np
      ! --- get tau
      taux = bconst*bxp(ip)
      tauy = bconst*byp(ip)
      tauz = bconst*bzp(ip)
      tausq = taux*taux+tauy*tauy+tauz*tauz
      ! --- get U',gamma'^2
      uxpr = uxp(ip) + const*exp(ip) + (uyp(ip)*tauz-uzp(ip)*tauy)*gi(ip)
      uypr = uyp(ip) + const*eyp(ip) + (uzp(ip)*taux-uxp(ip)*tauz)*gi(ip)
      uzpr = uzp(ip) + const*ezp(ip) + (uxp(ip)*tauy-uyp(ip)*taux)*gi(ip)
      gprsq = (1._num+(uxpr*uxpr+uypr*uypr+uzpr*uzpr)*invclightsq)
      !       --- get u*
      ust = (uxpr*taux+uypr*tauy+uzpr*tauz)*invclight
      ! --- get new gamma
      sigma = gprsq-tausq
      gisq = 2._num/(sigma+sqrt(sigma*sigma+4._num*(tausq+ust*ust)))
      gi(ip) = sqrt(gisq)
      !       --- get t,s
      bg = bconst*gi(ip)
      tx = bg*bxp(ip)
      ty = bg*byp(ip)
      tz = bg*bzp(ip)
      s = 1._num/(1._num+tausq*gisq)
      !  --- get t.u'
      tu = tx*uxpr+ty*uypr+tz*uzpr
      ! --- get new U
      uxp(ip) = s*(uxpr+tx*tu+uypr*tz-uzpr*ty)
      uyp(ip) = s*(uypr+ty*tu+uzpr*tx-uxpr*tz)
      uzp(ip) = s*(uzpr+tz*tu+uxpr*ty-uypr*tx)
    END DO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

  ELSE IF(which==1) THEN
  
    !     --- first half push
    const = 0.5_num*q*dt/m
    
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
    DO ip=1,np
      ! --- get new U
      vx = uxp(ip)*gi(ip)
      vy = uyp(ip)*gi(ip)
      vz = uzp(ip)*gi(ip)
      uxp(ip) = uxp(ip) + const*( exp(ip) + vy*bzp(ip)-vz*byp(ip) )
      uyp(ip) = uyp(ip) + const*( eyp(ip) + vz*bxp(ip)-vx*bzp(ip) )
      uzp(ip) = uzp(ip) + const*( ezp(ip) + vx*byp(ip)-vy*bxp(ip) )
      gi(ip) = 1./sqrt(1.+(uxp(ip)*uxp(ip)+uyp(ip)*uyp(ip)+uzp(ip)*uzp(ip))*invclightsq)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif


  ELSE IF(which==2) THEN
  !     --- second half push
    const = 0.5_num*q*dt/m
    bconst = const

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif

    DO ip=1,np
      !     --- get U'
      uxpr = uxp(ip) + const*exp(ip)
      uypr = uyp(ip) + const*eyp(ip)
      uzpr = uzp(ip) + const*ezp(ip)
      gprsq = (1_num+(uxpr*uxpr+uypr*uypr+uzpr*uzpr)*invclightsq)
      !       --- get tau
      taux = bconst*bxp(ip)
      tauy = bconst*byp(ip)
      tauz = bconst*bzp(ip)
      tausq = taux*taux+tauy*tauy+tauz*tauz
      !       --- get u*
      ust = (uxpr*taux+uypr*tauy+uzpr*tauz)*invclight
      !       --- get new gamma
      sigma = gprsq-tausq
      gisq = 2._num/(sigma+sqrt(sigma*sigma+4._num*(tausq+ust*ust)))
      gi(ip) = sqrt(gisq)
      !       --- get t,s
      bg = bconst*gi(ip)
      tx = bg*bxp(ip)
      ty = bg*byp(ip)
      tz = bg*bzp(ip)
      s = 1._num/(1._num+tausq*gisq)
      !       --- get t.u'
      tu = tx*uxpr+ty*uypr+tz*uzpr
      !       --- get new U
      uxp(ip) = s*(uxpr+tx*tu+uypr*tz-uzpr*ty)
      uyp(ip) = s*(uypr+ty*tu+uzpr*tx-uxpr*tz)
      uzp(ip) = s*(uzpr+tz*tu+uxpr*ty-uypr*tx)
      ENDDO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif
  ENDIF
  RETURN
END SUBROUTINE pxr_ebcancelpush3d


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
!> @param[in] xp,yp,zp particle position
!> @param[in] uxp,uyp,uzp particle momentum
!> @param[in] gaminv inverse of the particle Lorentz factor
!> @param[in] ex,ey,ez particle electric field
!> @param[in] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dtt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v performe the field interpolation at a lower order
!
SUBROUTINE field_gathering_plus_particle_pusher_1_1_1(np,xp,yp,zp,uxp,uyp,uzp,gaminv, &
                                      ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,dtt,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,q,m,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  USE particles

  ! ___ Parameter declaration ____________________________________
  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num)                            :: q,m
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv
  LOGICAL(lp)                          :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz,dtt

  INTEGER(isp)                         :: j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn,n
  INTEGER(idp)                         :: blocksize
  REAL(num)                            :: dxi, dyi, dzi
  REAL(num)                            :: x,y,z
  REAL(num)                            :: a
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: clghtisq,const1
  REAL(num), DIMENSION(lvect,0:1)      :: sx,sy,sz
  REAL(num), DIMENSION(lvect,0:1)      :: sx0,sy0,sz0
#if defined __INTEL_COMPILER
    !dir$ attributes align:64 :: sx,sy,sz,sx0,sy0,sz0
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

  sx0(:,0) = 1.0_num
  sy0(:,0) = 1.0_num
  sz0(:,0) = 1.0_num

  ! ______________________________________________________________________________________
  ! Loop on block of particles of size lvect
  DO ip=1,np,lvect

    blocksize = MIN(lvect,np-ip+1)

    ! ____________________________________________________________________________________
    ! Field gathering

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,blocksize
      nn=ip+n-1

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
      sx(n,0) = 1.0_num-xint
      sx(n,1) = xint
      sy(n,0) = 1.0_num-yint
      sy(n,1) = yint
      sz(n,0) = 1.0_num-zint
      sz(n,1) = zint

      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0

      ! Compute Ex on particle
      a = (sy(n,0)*exg(j0,k,l)   + sy(n,1)*exg(j0,k+1,l))*sz(n,0)  &
       + (sy(n,0)*exg(j0,k,l+1) + sy(n,1)*exg(j0,k+1,l+1))*sz(n,1)
      ex(nn) = ex(nn) + a*sx0(n,0)

      ! Compute Ey on particle
      a = (sx(n,0)*eyg(j,k0,l) + sx(n,1)*eyg(j+1,k0,l))*sz(n,0) &
       + (sx(n,0)*eyg(j,k0,l+1) + sx(n,1)*eyg(j+1,k0,l+1))*sz(n,1)
      ey(nn) = ey(nn) + a*sy0(n,0)

      ! Compute Ez on particle
      a = (sx(n,0)*ezg(j,k,l0) + sx(n,1)*ezg(j+1,k,l0))*sy(n,0) &
        + (sx(n,0)*ezg(j,k+1,l0) + sx(n,1)*ezg(j+1,k+1,l0))*sy(n,1)
      ez(nn) = ez(nn) + a*sz0(n,0)

      ! Compute Bx on particle
      a = (sx(n,0)*bxg(j,k0,l0) + sx(n,1)*bxg(j+1,k0,l0))*sy0(n,0)
      bx(nn) = bx(nn) + a*sz0(n,0)
    
      ! Compute By on particle
      a = (sy(n,0)*byg(j0,k,l0) + sy(n,1)*byg(j0,k+1,l0))*sx0(n,0)
      by(nn) = by(nn) + a*sz0(n,0)
    
      ! Compute Bz on particle
      a = (sz(n,0)*bzg(j0,k0,l) + sz(n,1)*bzg(j0,k0,l+1))*sx0(n,0)
      bz(nn) = bz(nn) + a*sy0(n,0)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

    ! ____________________________________________________________________________________
    ! Particle pusher
    
    SELECT CASE (particle_pusher)
    !! Vay pusher -- Full push
    CASE (1_idp)
      CALL pxr_ebcancelpush3d(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 gaminv(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),  &
                                 bx(ip:ip+blocksize-1),  &
                                 by(ip:ip+blocksize-1),  &
                                 bz(ip:ip+blocksize-1),q,m,dt,0_idp)
        
    !! Boris pusher -- Full push
    CASE DEFAULT

      ! ___ Push with E ___
      CALL pxr_epush_v(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),q,m,0.5_num*dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1))    

      ! ___ Push with B ___
      CALL pxr_bpush_v(blocksize,uxp(ip:ip+blocksize-1),   &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1),&
                                 bx(ip:ip+blocksize-1),    &
                                 by(ip:ip+blocksize-1),    &
                                 bz(ip:ip+blocksize-1),q,m,dt)

      ! ___ Push with E ___
      CALL pxr_epush_v(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),q,m,0.5_num*dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1))    
    END SELECT
    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize,xp(ip:ip+blocksize-1),  &
                               yp(ip:ip+blocksize-1),  &
                               zp(ip:ip+blocksize-1),  &
                               uxp(ip:ip+blocksize-1), &
                               uyp(ip:ip+blocksize-1),  &
                               uzp(ip:ip+blocksize-1),  &
                               gaminv(ip:ip+blocksize-1),dt)

    ! ___ Push with E + gamma ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,ex,ey,ez)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       uxp(nn) = uxp(nn) + ex(nn)*const1
!       uyp(nn) = uyp(nn) + ey(nn)*const1
!       uzp(nn) = uzp(nn) + ez(nn)*const1
! 
!       usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
!       gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
! 
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Push with B ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       const2 = gaminv(nn)*const1
!       tx = bx(nn)*const2
!       ty = by(nn)*const2
!       tz = bz(nn)*const2
!       tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
!       wx = tx*tsqi
!       wy = ty*tsqi
!       wz = tz*tsqi
!       uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
!       uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
!       uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
!       uxp(nn) = uxp(nn) + uyppr*wz - uzppr*wy
!       uyp(nn) = uyp(nn) + uzppr*wx - uxppr*wz
!       uzp(nn) = uzp(nn) + uxppr*wy - uyppr*wx
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Push with E + gamma ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,ex,ey,ez)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       uxp(nn) = uxp(nn) + ex(nn)*const1
!       uyp(nn) = uyp(nn) + ey(nn)*const1
!       uzp(nn) = uzp(nn) + ez(nn)*const1
! 
!       usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
!       gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Update position ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       const2 = gaminv(nn)*dtt
!       xp(nn) = xp(nn) + uxp(nn)*const2
!       yp(nn) = yp(nn) + uyp(nn)*const2
!       zp(nn) = zp(nn) + uzp(nn)*const2
! 
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

  ! End loop on particles
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
!> @param[in] xp,yp,zp particle position
!> @param[in] uxp,uyp,uzp particle momentum
!> @param[in] gaminv inverse of the particle Lorentz factor
!> @param[in] ex,ey,ez particle electric field
!> @param[in] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dtt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v performe the field interpolation at a lower order
SUBROUTINE field_gathering_plus_particle_pusher_2_2_2(np,xp,yp,zp,uxp,uyp,uzp,gaminv, &
                                      ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,dtt,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,q,m,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  USE particles

  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num)                            :: q,m
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv
  LOGICAL(lp)                          :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz,dtt
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn,n
  INTEGER(idp)                         :: blocksize
  INTEGER(isp)                         :: j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq,yintsq,zintsq
  REAL(num)                            :: a
  REAL(num)                            :: clghtisq,const1
  REAL(num), DIMENSION(lvect,-1:1)     :: sx,sy,sz
  REAL(num), DIMENSION(lvect,-1:1)     :: sx0,sy0,sz0
#if defined __INTEL_COMPILER
    !dir$ attributes align:64 :: sx,sy,sz,sx0,sy0,sz0
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
  DO ip=1,np,lvect

    blocksize = MIN(lvect,np-ip+1)

    ! ____________________________________________________________________________________
    ! Field gathering

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,blocksize
      nn=ip+n-1

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
      sx(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx(n, 0) = 0.75_num-xintsq
      sx(n, 1) = 0.5_num*(0.5_num+xint)**2

      yintsq = yint*yint
      sy(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy(n, 0) = 0.75_num-yintsq
      sy(n, 1) = 0.5_num*(0.5_num+yint)**2

      zintsq = zint*zint
      sz(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz(n, 0) = 0.75_num-zintsq
      sz(n, 1) = 0.5_num*(0.5_num+zint)**2

      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0

      sx0(n, 0) = 1.0_num-xint
      sx0(n, 1) = xint

      sy0(n, 0) = 1.0_num-yint
      sy0(n, 1) = yint

      sz0(n, 0) = 1.0_num-zint
      sz0(n, 1) = zint

      ! Compute Ex on particle
      a = (sx0(n,0)*exg(j0,k-1,l-1) &
          + sx0(n,1)*exg(j0+1,k-1,l-1))*sy(n,-1)
      a = a + (sx0(n,0)*exg(j0,k,l-1) &
          + sx0(n,1)*exg(j0+1,k,l-1))*sy(n,0)
      a = a + (sx0(n,0)*exg(j0,k+1,l-1) &
          + sx0(n,1)*exg(j0+1,k+1,l-1))*sy(n,1)
      ex(nn) = ex(nn) + a*sz(n,-1)
      a = (sx0(n,0)*exg(j0,k-1,l) &
          + sx0(n,1)*exg(j0+1,k-1,l))*sy(n,-1)
      a = a + (sx0(n,0)*exg(j0,k,l) &
          + sx0(n,1)*exg(j0+1,k,l))*sy(n,0)
      a = a + (sx0(n,0)*exg(j0,k+1,l) &
          + sx0(n,1)*exg(j0+1,k+1,l))*sy(n,1)
      ex(nn) = ex(nn) + a*sz(n,0)
      a = (sx0(n,0)*exg(j0,k-1,l+1) &
          + sx0(n,1)*exg(j0+1,k-1,l+1))*sy(n,-1)
      a = a + (sx0(n,0)*exg(j0,k,l+1) &
          + sx0(n,1)*exg(j0+1,k,l+1))*sy(n,0)
      a = a + (sx0(n,0)*exg(j0,k+1,l+1) &
          + sx0(n,1)*exg(j0+1,k+1,l+1))*sy(n,1)
      ex(nn) = ex(nn) + a*sz(n,1)
    
      ! Compute Ey on particle
      a = (sx(n,-1)*eyg(j-1,k0,l-1) &
          + sx(n,0)*eyg(j,k0,l-1) &
          + sx(n,1)*eyg(j+1,k0,l-1))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l-1) &
          + sx(n,0)*eyg(j,k0+1,l-1) &
          + sx(n,1)*eyg(j+1,k0+1,l-1))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,-1)
      a = (sx(n,-1)*eyg(j-1,k0,l) &
          + sx(n,0)*eyg(j,k0,l) &
          + sx(n,1)*eyg(j+1,k0,l))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l) &
          + sx(n,0)*eyg(j,k0+1,l) &
          + sx(n,1)*eyg(j+1,k0+1,l))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,0)
      a = (sx(n,-1)*eyg(j-1,k0,l+1) &
          + sx(n,0)*eyg(j,k0,l+1) &
          + sx(n,1)*eyg(j+1,k0,l+1))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l+1) &
          + sx(n,0)*eyg(j,k0+1,l+1) &
          + sx(n,1)*eyg(j+1,k0+1,l+1))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,1)
    
      ! Compute Ez on particle
      a = (sx(n,-1)*ezg(j-1,k-1,l0) &
          + sx(n,0)*ezg(j,k-1,l0) &
          + sx(n,1)*ezg(j+1,k-1,l0))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0) &
          + sx(n,0)*ezg(j,k,l0) &
          + sx(n,1)*ezg(j+1,k,l0))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0) &
          + sx(n,0)*ezg(j,k+1,l0) &
          + sx(n,1)*ezg(j+1,k+1,l0))*sy(n,1)
      ez(nn) = ez(nn) + a*sz0(n,0)
      a = (sx(n,-1)*ezg(j-1,k-1,l0+1) &
          + sx(n,0)*ezg(j,k-1,l0+1) &
          + sx(n,1)*ezg(j+1,k-1,l0+1))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0+1) &
          + sx(n,0)*ezg(j,k,l0+1) &
          + sx(n,1)*ezg(j+1,k,l0+1))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0+1) &
          + sx(n,0)*ezg(j,k+1,l0+1) &
          + sx(n,1)*ezg(j+1,k+1,l0+1))*sy(n,1)
      ez(nn) = ez(nn) + a*sz0(n,1)
    
      ! Compute Bx on particle
      a = (sx(n,-1)*bxg(j-1,k0,l0) &
          + sx(n,0)*bxg(j,k0,l0) &
          + sx(n,1)*bxg(j+1,k0,l0))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0) &
          + sx(n,0)*bxg(j,k0+1,l0) &
          + sx(n,1)*bxg(j+1,k0+1,l0))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,0)
      a = (sx(n,-1)*bxg(j-1,k0,l0+1) &
          + sx(n,0)*bxg(j,k0,l0+1) &
          + sx(n,1)*bxg(j+1,k0,l0+1))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0+1) &
          + sx(n,0)*bxg(j,k0+1,l0+1) &
          + sx(n,1)*bxg(j+1,k0+1,l0+1))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,1)
    
      ! Compute By on particle
      a = (sx0(n,0)*byg(j0,k-1,l0) &
          + sx0(n,1)*byg(j0+1,k-1,l0))*sy(n,-1)
      a = a + (sx0(n,0)*byg(j0,k,l0) &
          + sx0(n,1)*byg(j0+1,k,l0))*sy(n,0)
      a = a + (sx0(n,0)*byg(j0,k+1,l0) &
          + sx0(n,1)*byg(j0+1,k+1,l0))*sy(n,1)
      by(nn) = by(nn) + a*sz0(n,0)
      a = (sx0(n,0)*byg(j0,k-1,l0+1) &
          + sx0(n,1)*byg(j0+1,k-1,l0+1))*sy(n,-1)
      a = a + (sx0(n,0)*byg(j0,k,l0+1) &
          + sx0(n,1)*byg(j0+1,k,l0+1))*sy(n,0)
      a = a + (sx0(n,0)*byg(j0,k+1,l0+1) &
          + sx0(n,1)*byg(j0+1,k+1,l0+1))*sy(n,1)
      by(nn) = by(nn) + a*sz0(n,1)
    
      ! Compute Bz on particle
      a = (sx0(n,0)*bzg(j0,k0,l-1) &
          + sx0(n,1)*bzg(j0+1,k0,l-1))*sy0(n,0)
      a = a + (sx0(n,0)*bzg(j0,k0+1,l-1) &
          + sx0(n,1)*bzg(j0+1,k0+1,l-1))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,-1)
      a = (sx0(n,0)*bzg(j0,k0,l) &
          + sx0(n,1)*bzg(j0+1,k0,l))*sy0(n,0)
      a = a + (sx0(n,0)*bzg(j0,k0+1,l) &
          + sx0(n,1)*bzg(j0+1,k0+1,l))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,0)
      a = (sx0(n,0)*bzg(j0,k0,l+1) &
          + sx0(n,1)*bzg(j0+1,k0,l+1))*sy0(n,0)
      a = a + (sx0(n,0)*bzg(j0,k0+1,l+1) &
          + sx0(n,1)*bzg(j0+1,k0+1,l+1))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,1)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

    ! ____________________________________________________________________________________
    ! Particle pusher

    SELECT CASE (particle_pusher)
    !! Vay pusher -- Full push
    CASE (1_idp)
      CALL pxr_ebcancelpush3d(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 gaminv(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),  &
                                 bx(ip:ip+blocksize-1),  &
                                 by(ip:ip+blocksize-1),  &
                                 bz(ip:ip+blocksize-1),q,m,dt,0_idp)
        
    !! Boris pusher -- Full push
    CASE DEFAULT

      ! ___ Push with E ___
      CALL pxr_epush_v(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),q,m,0.5_num*dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1))    

      ! ___ Push with B ___
      CALL pxr_bpush_v(blocksize,uxp(ip:ip+blocksize-1),   &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1),&
                                 bx(ip:ip+blocksize-1),    &
                                 by(ip:ip+blocksize-1),    &
                                 bz(ip:ip+blocksize-1),q,m,dt)

      ! ___ Push with E ___
      CALL pxr_epush_v(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),q,m,0.5_num*dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1))    
    END SELECT
    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize,xp(ip:ip+blocksize-1),  &
                               yp(ip:ip+blocksize-1),  &
                               zp(ip:ip+blocksize-1),  &
                               uxp(ip:ip+blocksize-1), &
                               uyp(ip:ip+blocksize-1),  &
                               uzp(ip:ip+blocksize-1),  &
                               gaminv(ip:ip+blocksize-1),dt)

    ! ___ Push with E + gamma ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,ex,ey,ez)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       uxp(nn) = uxp(nn) + ex(nn)*const1
!       uyp(nn) = uyp(nn) + ey(nn)*const1
!       uzp(nn) = uzp(nn) + ez(nn)*const1
! 
!       usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
!       gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
! 
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Push with B ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       const2 = gaminv(nn)*const1
!       tx = bx(nn)*const2
!       ty = by(nn)*const2
!       tz = bz(nn)*const2
!       tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
!       wx = tx*tsqi
!       wy = ty*tsqi
!       wz = tz*tsqi
!       uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
!       uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
!       uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
!       uxp(nn) = uxp(nn) + uyppr*wz - uzppr*wy
!       uyp(nn) = uyp(nn) + uzppr*wx - uxppr*wz
!       uzp(nn) = uzp(nn) + uxppr*wy - uyppr*wx
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Push with E + gamma ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,ex,ey,ez)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       uxp(nn) = uxp(nn) + ex(nn)*const1
!       uyp(nn) = uyp(nn) + ey(nn)*const1
!       uzp(nn) = uzp(nn) + ez(nn)*const1
! 
!       usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
!       gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Update position ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       const2 = gaminv(nn)*dtt
!       xp(nn) = xp(nn) + uxp(nn)*const2
!       yp(nn) = yp(nn) + uyp(nn)*const2
!       zp(nn) = zp(nn) + uzp(nn)*const2
! 
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

  ! End loop on particles
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
!> @param[in] xp,yp,zp particle position
!> @param[in] uxp,uyp,uzp particle momentum
!> @param[in] gaminv inverse of the particle Lorentz factor
!> @param[in] ex,ey,ez particle electric field
!> @param[in] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dtt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v performe the field interpolation at a lower order
!
SUBROUTINE field_gathering_plus_particle_pusher_3_3_3(np,xp,yp,zp,uxp,uyp,uzp,gaminv, &
                                      ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,dtt,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,q,m,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  USE particles

  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num)                            :: q,m
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv
  LOGICAL(lp)                          :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz,dtt
  INTEGER(isp)                         :: ip
  INTEGER(idp)                         :: blocksize
  INTEGER(isp)                         :: nn,n
  INTEGER(isp)                         :: j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  REAL(num)                            :: clghtisq,const1
  REAL(num)                            :: a
  REAL(num), DIMENSION(lvect,-1:2)     :: sx
  REAL(num), DIMENSION(lvect,-1:2)     :: sy
  REAL(num), DIMENSION(lvect,-1:2)     :: sz
  REAL(num), DIMENSION(lvect,-1:1)     :: sx0
  REAL(num), DIMENSION(lvect,-1:1)     :: sy0
  REAL(num), DIMENSION(lvect,-1:1)     :: sz0
#if defined __INTEL_COMPILER
    !dir$ attributes align:64 :: sx,sy,sz,sx0,sy0,sz0
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
  DO ip=1,np,lvect

    blocksize = MIN(lvect,np-ip+1)

    ! ____________________________________________________________________________________
    ! Field gathering

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,blocksize
      nn=ip+n-1

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
      sx(n,-1) = onesixth*oxintsq*oxint
      sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx(n, 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(n,-1) = onesixth*oyintsq*oyint
      sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy(n, 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(n,-1) = onesixth*ozintsq*ozint
      sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz(n, 2) = onesixth*zintsq*zint

      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0

      xintsq = xint*xint
      sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx0(n, 0) = 0.75_num-xintsq
      sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

      yintsq = yint*yint
      sy0(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy0(n, 0) = 0.75_num-yintsq
      sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

      zintsq = zint*zint
      sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz0(n, 0) = 0.75_num-zintsq
      sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

      ! Compute Ex on particle
      a = (sx0(n,-1)*exg(j0-1,k-1,l-1) &
          + sx0(n,0)*exg(j0,k-1,l-1) &
          + sx0(n,1)*exg(j0+1,k-1,l-1))*sy(n,-1)
      a = a + (sx0(n,-1)*exg(j0-1,k,l-1) &
          + sx0(n,0)*exg(j0,k,l-1) &
          + sx0(n,1)*exg(j0+1,k,l-1))*sy(n,0)
      a = a + (sx0(n,-1)*exg(j0-1,k+1,l-1) &
          + sx0(n,0)*exg(j0,k+1,l-1) &
          + sx0(n,1)*exg(j0+1,k+1,l-1))*sy(n,1)
      a = a + (sx0(n,-1)*exg(j0-1,k+2,l-1) &
          + sx0(n,0)*exg(j0,k+2,l-1) &
          + sx0(n,1)*exg(j0+1,k+2,l-1))*sy(n,2)
      ex(nn) = ex(nn) + a*sz(n,-1)
      a = (sx0(n,-1)*exg(j0-1,k-1,l) &
          + sx0(n,0)*exg(j0,k-1,l) &
          + sx0(n,1)*exg(j0+1,k-1,l))*sy(n,-1)
      a = a + (sx0(n,-1)*exg(j0-1,k,l) &
          + sx0(n,0)*exg(j0,k,l) &
          + sx0(n,1)*exg(j0+1,k,l))*sy(n,0)
      a = a + (sx0(n,-1)*exg(j0-1,k+1,l) &
          + sx0(n,0)*exg(j0,k+1,l) &
          + sx0(n,1)*exg(j0+1,k+1,l))*sy(n,1)
      a = a + (sx0(n,-1)*exg(j0-1,k+2,l) &
          + sx0(n,0)*exg(j0,k+2,l) &
          + sx0(n,1)*exg(j0+1,k+2,l))*sy(n,2)
      ex(nn) = ex(nn) + a*sz(n,0)
      a = (sx0(n,-1)*exg(j0-1,k-1,l+1) &
          + sx0(n,0)*exg(j0,k-1,l+1) &
          + sx0(n,1)*exg(j0+1,k-1,l+1))*sy(n,-1)
      a = a + (sx0(n,-1)*exg(j0-1,k,l+1) &
          + sx0(n,0)*exg(j0,k,l+1) &
          + sx0(n,1)*exg(j0+1,k,l+1))*sy(n,0)
      a = a + (sx0(n,-1)*exg(j0-1,k+1,l+1) &
          + sx0(n,0)*exg(j0,k+1,l+1) &
          + sx0(n,1)*exg(j0+1,k+1,l+1))*sy(n,1)
      a = a + (sx0(n,-1)*exg(j0-1,k+2,l+1) &
          + sx0(n,0)*exg(j0,k+2,l+1) &
          + sx0(n,1)*exg(j0+1,k+2,l+1))*sy(n,2)
      ex(nn) = ex(nn) + a*sz(n,1)
      a = (sx0(n,-1)*exg(j0-1,k-1,l+2) &
          + sx0(n,0)*exg(j0,k-1,l+2) &
          + sx0(n,1)*exg(j0+1,k-1,l+2))*sy(n,-1)
      a = a + (sx0(n,-1)*exg(j0-1,k,l+2) &
          + sx0(n,0)*exg(j0,k,l+2) &
          + sx0(n,1)*exg(j0+1,k,l+2))*sy(n,0)
      a = a + (sx0(n,-1)*exg(j0-1,k+1,l+2) &
          + sx0(n,0)*exg(j0,k+1,l+2) &
          + sx0(n,1)*exg(j0+1,k+1,l+2))*sy(n,1)
      a = a + (sx0(n,-1)*exg(j0-1,k+2,l+2) &
          + sx0(n,0)*exg(j0,k+2,l+2) &
          + sx0(n,1)*exg(j0+1,k+2,l+2))*sy(n,2)
      ex(nn) = ex(nn) + a*sz(n,2)

      ! Compute Ey on particle
      a = (sx(n,-1)*eyg(j-1,k0-1,l-1) &
          + sx(n,0)*eyg(j,k0-1,l-1) &
          + sx(n,1)*eyg(j+1,k0-1,l-1) &
          + sx(n,2)*eyg(j+2,k0-1,l-1))*sy0(n,-1)
      a = a + (sx(n,-1)*eyg(j-1,k0,l-1) &
          + sx(n,0)*eyg(j,k0,l-1) &
          + sx(n,1)*eyg(j+1,k0,l-1) &
          + sx(n,2)*eyg(j+2,k0,l-1))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l-1) &
          + sx(n,0)*eyg(j,k0+1,l-1) &
          + sx(n,1)*eyg(j+1,k0+1,l-1) &
          + sx(n,2)*eyg(j+2,k0+1,l-1))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,-1)
      a = (sx(n,-1)*eyg(j-1,k0-1,l) &
          + sx(n,0)*eyg(j,k0-1,l) &
          + sx(n,1)*eyg(j+1,k0-1,l) &
          + sx(n,2)*eyg(j+2,k0-1,l))*sy0(n,-1)
      a = a + (sx(n,-1)*eyg(j-1,k0,l) &
          + sx(n,0)*eyg(j,k0,l) &
          + sx(n,1)*eyg(j+1,k0,l) &
          + sx(n,2)*eyg(j+2,k0,l))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l) &
          + sx(n,0)*eyg(j,k0+1,l) &
          + sx(n,1)*eyg(j+1,k0+1,l) &
          + sx(n,2)*eyg(j+2,k0+1,l))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,0)
      a = (sx(n,-1)*eyg(j-1,k0-1,l+1) &
          + sx(n,0)*eyg(j,k0-1,l+1) &
          + sx(n,1)*eyg(j+1,k0-1,l+1) &
          + sx(n,2)*eyg(j+2,k0-1,l+1))*sy0(n,-1)
      a = a + (sx(n,-1)*eyg(j-1,k0,l+1) &
          + sx(n,0)*eyg(j,k0,l+1) &
          + sx(n,1)*eyg(j+1,k0,l+1) &
          + sx(n,2)*eyg(j+2,k0,l+1))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l+1) &
          + sx(n,0)*eyg(j,k0+1,l+1) &
          + sx(n,1)*eyg(j+1,k0+1,l+1) &
          + sx(n,2)*eyg(j+2,k0+1,l+1))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,1)
      a = (sx(n,-1)*eyg(j-1,k0-1,l+2) &
          + sx(n,0)*eyg(j,k0-1,l+2) &
          + sx(n,1)*eyg(j+1,k0-1,l+2) &
          + sx(n,2)*eyg(j+2,k0-1,l+2))*sy0(n,-1)
      a = a + (sx(n,-1)*eyg(j-1,k0,l+2) &
          + sx(n,0)*eyg(j,k0,l+2) &
          + sx(n,1)*eyg(j+1,k0,l+2) &
          + sx(n,2)*eyg(j+2,k0,l+2))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l+2) &
          + sx(n,0)*eyg(j,k0+1,l+2) &
          + sx(n,1)*eyg(j+1,k0+1,l+2) &
          + sx(n,2)*eyg(j+2,k0+1,l+2))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,2)

      ! Compute Ez on particle
      a = (sx(n,-1)*ezg(j-1,k-1,l0-1) &
          + sx(n,0)*ezg(j,k-1,l0-1) &
          + sx(n,1)*ezg(j+1,k-1,l0-1) &
          + sx(n,2)*ezg(j+2,k-1,l0-1))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0-1) &
          + sx(n,0)*ezg(j,k,l0-1) &
          + sx(n,1)*ezg(j+1,k,l0-1) &
          + sx(n,2)*ezg(j+2,k,l0-1))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0-1) &
          + sx(n,0)*ezg(j,k+1,l0-1) &
          + sx(n,1)*ezg(j+1,k+1,l0-1) &
          + sx(n,2)*ezg(j+2,k+1,l0-1))*sy(n,1)
      a = a + (sx(n,-1)*ezg(j-1,k+2,l0-1) &
          + sx(n,0)*ezg(j,k+2,l0-1) &
          + sx(n,1)*ezg(j+1,k+2,l0-1) &
          + sx(n,2)*ezg(j+2,k+2,l0-1))*sy(n,2)
      ez(nn) = ez(nn) + a*sz0(n,-1)
      a = (sx(n,-1)*ezg(j-1,k-1,l0) &
          + sx(n,0)*ezg(j,k-1,l0) &
          + sx(n,1)*ezg(j+1,k-1,l0) &
          + sx(n,2)*ezg(j+2,k-1,l0))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0) &
          + sx(n,0)*ezg(j,k,l0) &
          + sx(n,1)*ezg(j+1,k,l0) &
          + sx(n,2)*ezg(j+2,k,l0))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0) &
          + sx(n,0)*ezg(j,k+1,l0) &
          + sx(n,1)*ezg(j+1,k+1,l0) &
          + sx(n,2)*ezg(j+2,k+1,l0))*sy(n,1)
      a = a + (sx(n,-1)*ezg(j-1,k+2,l0) &
          + sx(n,0)*ezg(j,k+2,l0) &
          + sx(n,1)*ezg(j+1,k+2,l0) &
          + sx(n,2)*ezg(j+2,k+2,l0))*sy(n,2)
      ez(nn) = ez(nn) + a*sz0(n,0)
      a = (sx(n,-1)*ezg(j-1,k-1,l0+1) &
          + sx(n,0)*ezg(j,k-1,l0+1) &
          + sx(n,1)*ezg(j+1,k-1,l0+1) &
          + sx(n,2)*ezg(j+2,k-1,l0+1))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0+1) &
          + sx(n,0)*ezg(j,k,l0+1) &
          + sx(n,1)*ezg(j+1,k,l0+1) &
          + sx(n,2)*ezg(j+2,k,l0+1))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0+1) &
          + sx(n,0)*ezg(j,k+1,l0+1) &
          + sx(n,1)*ezg(j+1,k+1,l0+1) &
          + sx(n,2)*ezg(j+2,k+1,l0+1))*sy(n,1)
      a = a + (sx(n,-1)*ezg(j-1,k+2,l0+1) &
          + sx(n,0)*ezg(j,k+2,l0+1) &
          + sx(n,1)*ezg(j+1,k+2,l0+1) &
          + sx(n,2)*ezg(j+2,k+2,l0+1))*sy(n,2)
      ez(nn) = ez(nn) + a*sz0(n,1)

      ! Compute Bx on particle
      a = (sx(n,-1)*bxg(j-1,k0-1,l0-1) &
          + sx(n,0)*bxg(j,k0-1,l0-1) &
          + sx(n,1)*bxg(j+1,k0-1,l0-1) &
          + sx(n,2)*bxg(j+2,k0-1,l0-1))*sy0(n,-1)
      a = a + (sx(n,-1)*bxg(j-1,k0,l0-1) &
          + sx(n,0)*bxg(j,k0,l0-1) &
          + sx(n,1)*bxg(j+1,k0,l0-1) &
          + sx(n,2)*bxg(j+2,k0,l0-1))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0-1) &
          + sx(n,0)*bxg(j,k0+1,l0-1) &
          + sx(n,1)*bxg(j+1,k0+1,l0-1) &
          + sx(n,2)*bxg(j+2,k0+1,l0-1))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,-1)
      a = (sx(n,-1)*bxg(j-1,k0-1,l0) &
          + sx(n,0)*bxg(j,k0-1,l0) &
          + sx(n,1)*bxg(j+1,k0-1,l0) &
          + sx(n,2)*bxg(j+2,k0-1,l0))*sy0(n,-1)
      a = a + (sx(n,-1)*bxg(j-1,k0,l0) &
          + sx(n,0)*bxg(j,k0,l0) &
          + sx(n,1)*bxg(j+1,k0,l0) &
          + sx(n,2)*bxg(j+2,k0,l0))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0) &
          + sx(n,0)*bxg(j,k0+1,l0) &
          + sx(n,1)*bxg(j+1,k0+1,l0) &
          + sx(n,2)*bxg(j+2,k0+1,l0))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,0)
      a = (sx(n,-1)*bxg(j-1,k0-1,l0+1) &
          + sx(n,0)*bxg(j,k0-1,l0+1) &
          + sx(n,1)*bxg(j+1,k0-1,l0+1) &
          + sx(n,2)*bxg(j+2,k0-1,l0+1))*sy0(n,-1)
      a = a + (sx(n,-1)*bxg(j-1,k0,l0+1) &
          + sx(n,0)*bxg(j,k0,l0+1) &
          + sx(n,1)*bxg(j+1,k0,l0+1) &
          + sx(n,2)*bxg(j+2,k0,l0+1))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0+1) &
          + sx(n,0)*bxg(j,k0+1,l0+1) &
          + sx(n,1)*bxg(j+1,k0+1,l0+1) &
          + sx(n,2)*bxg(j+2,k0+1,l0+1))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,1)

      ! Compute By on particle
      a = (sx0(n,-1)*byg(j0-1,k-1,l0-1) &
          + sx0(n,0)*byg(j0,k-1,l0-1) &
          + sx0(n,1)*byg(j0+1,k-1,l0-1))*sy(n,-1)
      a = a + (sx0(n,-1)*byg(j0-1,k,l0-1) &
          + sx0(n,0)*byg(j0,k,l0-1) &
          + sx0(n,1)*byg(j0+1,k,l0-1))*sy(n,0)
      a = a + (sx0(n,-1)*byg(j0-1,k+1,l0-1) &
          + sx0(n,0)*byg(j0,k+1,l0-1) &
          + sx0(n,1)*byg(j0+1,k+1,l0-1))*sy(n,1)
      a = a + (sx0(n,-1)*byg(j0-1,k+2,l0-1) &
          + sx0(n,0)*byg(j0,k+2,l0-1) &
          + sx0(n,1)*byg(j0+1,k+2,l0-1))*sy(n,2)
      by(nn) = by(nn) + a*sz0(n,-1)
      a = (sx0(n,-1)*byg(j0-1,k-1,l0) &
          + sx0(n,0)*byg(j0,k-1,l0) &
          + sx0(n,1)*byg(j0+1,k-1,l0))*sy(n,-1)
      a = a + (sx0(n,-1)*byg(j0-1,k,l0) &
          + sx0(n,0)*byg(j0,k,l0) &
          + sx0(n,1)*byg(j0+1,k,l0))*sy(n,0)
      a = a + (sx0(n,-1)*byg(j0-1,k+1,l0) &
          + sx0(n,0)*byg(j0,k+1,l0) &
          + sx0(n,1)*byg(j0+1,k+1,l0))*sy(n,1)
      a = a + (sx0(n,-1)*byg(j0-1,k+2,l0) &
          + sx0(n,0)*byg(j0,k+2,l0) &
          + sx0(n,1)*byg(j0+1,k+2,l0))*sy(n,2)
      by(nn) = by(nn) + a*sz0(n,0)
      a = (sx0(n,-1)*byg(j0-1,k-1,l0+1) &
          + sx0(n,0)*byg(j0,k-1,l0+1) &
          + sx0(n,1)*byg(j0+1,k-1,l0+1))*sy(n,-1)
      a = a + (sx0(n,-1)*byg(j0-1,k,l0+1) &
          + sx0(n,0)*byg(j0,k,l0+1) &
          + sx0(n,1)*byg(j0+1,k,l0+1))*sy(n,0)
      a = a + (sx0(n,-1)*byg(j0-1,k+1,l0+1) &
          + sx0(n,0)*byg(j0,k+1,l0+1) &
          + sx0(n,1)*byg(j0+1,k+1,l0+1))*sy(n,1)
      a = a + (sx0(n,-1)*byg(j0-1,k+2,l0+1) &
          + sx0(n,0)*byg(j0,k+2,l0+1) &
          + sx0(n,1)*byg(j0+1,k+2,l0+1))*sy(n,2)
      by(nn) = by(nn) + a*sz0(n,1)

      ! Compute Bz on particle
      a = (sx0(n,-1)*bzg(j0-1,k0-1,l-1) &
          + sx0(n,0)*bzg(j0,k0-1,l-1) &
          + sx0(n,1)*bzg(j0+1,k0-1,l-1))*sy0(n,-1)
      a = a + (sx0(n,-1)*bzg(j0-1,k0,l-1) &
          + sx0(n,0)*bzg(j0,k0,l-1) &
          + sx0(n,1)*bzg(j0+1,k0,l-1))*sy0(n,0)
      a = a + (sx0(n,-1)*bzg(j0-1,k0+1,l-1) &
          + sx0(n,0)*bzg(j0,k0+1,l-1) &
          + sx0(n,1)*bzg(j0+1,k0+1,l-1))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,-1)
      a = (sx0(n,-1)*bzg(j0-1,k0-1,l) &
          + sx0(n,0)*bzg(j0,k0-1,l) &
          + sx0(n,1)*bzg(j0+1,k0-1,l))*sy0(n,-1)
      a = a + (sx0(n,-1)*bzg(j0-1,k0,l) &
          + sx0(n,0)*bzg(j0,k0,l) &
          + sx0(n,1)*bzg(j0+1,k0,l))*sy0(n,0)
      a = a + (sx0(n,-1)*bzg(j0-1,k0+1,l) &
          + sx0(n,0)*bzg(j0,k0+1,l) &
          + sx0(n,1)*bzg(j0+1,k0+1,l))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,0)
      a = (sx0(n,-1)*bzg(j0-1,k0-1,l+1) &
          + sx0(n,0)*bzg(j0,k0-1,l+1) &
          + sx0(n,1)*bzg(j0+1,k0-1,l+1))*sy0(n,-1)
      a = a + (sx0(n,-1)*bzg(j0-1,k0,l+1) &
          + sx0(n,0)*bzg(j0,k0,l+1) &
          + sx0(n,1)*bzg(j0+1,k0,l+1))*sy0(n,0)
      a = a + (sx0(n,-1)*bzg(j0-1,k0+1,l+1) &
          + sx0(n,0)*bzg(j0,k0+1,l+1) &
          + sx0(n,1)*bzg(j0+1,k0+1,l+1))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,1)
      a = (sx0(n,-1)*bzg(j0-1,k0-1,l+2) &
          + sx0(n,0)*bzg(j0,k0-1,l+2) &
          + sx0(n,1)*bzg(j0+1,k0-1,l+2))*sy0(n,-1)
      a = a + (sx0(n,-1)*bzg(j0-1,k0,l+2) &
          + sx0(n,0)*bzg(j0,k0,l+2) &
          + sx0(n,1)*bzg(j0+1,k0,l+2))*sy0(n,0)
      a = a + (sx0(n,-1)*bzg(j0-1,k0+1,l+2) &
          + sx0(n,0)*bzg(j0,k0+1,l+2) &
          + sx0(n,1)*bzg(j0+1,k0+1,l+2))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,2)


    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

    ! ____________________________________________________________________________________
    ! Particle pusher

    SELECT CASE (particle_pusher)
    !! Vay pusher -- Full push
    CASE (1_idp)
      CALL pxr_ebcancelpush3d(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 gaminv(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),  &
                                 bx(ip:ip+blocksize-1),  &
                                 by(ip:ip+blocksize-1),  &
                                 bz(ip:ip+blocksize-1),q,m,dt,0_idp)
        
    !! Boris pusher -- Full push
    CASE DEFAULT

      ! ___ Push with E ___
      CALL pxr_epush_v(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),q,m,0.5_num*dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1))    

      ! ___ Push with B ___
      CALL pxr_bpush_v(blocksize,uxp(ip:ip+blocksize-1),   &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1),&
                                 bx(ip:ip+blocksize-1),    &
                                 by(ip:ip+blocksize-1),    &
                                 bz(ip:ip+blocksize-1),q,m,dt)

      ! ___ Push with E ___
      CALL pxr_epush_v(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1), &
                                 uzp(ip:ip+blocksize-1), &
                                 ex(ip:ip+blocksize-1),  &
                                 ey(ip:ip+blocksize-1),  &
                                 ez(ip:ip+blocksize-1),q,m,0.5_num*dt)

      ! ___ compute Gamma ___
      CALL pxr_set_gamma(blocksize,uxp(ip:ip+blocksize-1), &
                                 uyp(ip:ip+blocksize-1),   &
                                 uzp(ip:ip+blocksize-1),   &
                                 gaminv(ip:ip+blocksize-1))    
    END SELECT
    ! ___ Update position ___
    CALL pxr_pushxyz(blocksize,xp(ip:ip+blocksize-1),  &
                               yp(ip:ip+blocksize-1),  &
                               zp(ip:ip+blocksize-1),  &
                               uxp(ip:ip+blocksize-1), &
                               uyp(ip:ip+blocksize-1),  &
                               uzp(ip:ip+blocksize-1),  &
                               gaminv(ip:ip+blocksize-1),dt)
    ! ___ Push with E + gamma ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,ex,ey,ez)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       uxp(nn) = uxp(nn) + ex(nn)*const1
!       uyp(nn) = uyp(nn) + ey(nn)*const1
!       uzp(nn) = uzp(nn) + ez(nn)*const1
! 
!       usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
!       gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
! 
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Push with B ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       const2 = gaminv(nn)*const1
!       tx = bx(nn)*const2
!       ty = by(nn)*const2
!       tz = bz(nn)*const2
!       tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
!       wx = tx*tsqi
!       wy = ty*tsqi
!       wz = tz*tsqi
!       uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
!       uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
!       uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
!       uxp(nn) = uxp(nn) + uyppr*wz - uzppr*wy
!       uyp(nn) = uyp(nn) + uzppr*wx - uxppr*wz
!       uzp(nn) = uzp(nn) + uxppr*wy - uyppr*wx
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Push with E + gamma ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,ex,ey,ez)
!       !IBM* ALIGN(64,bx,by,bz)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       uxp(nn) = uxp(nn) + ex(nn)*const1
!       uyp(nn) = uyp(nn) + ey(nn)*const1
!       uzp(nn) = uzp(nn) + ez(nn)*const1
! 
!       usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
!       gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

    ! ___ Update position ___
! #if defined __INTEL_COMPILER
!       !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
!       !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
!       !DIR$ ASSUME_ALIGNED gaminv:64
!       !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)
!       !DIR$ SIMD VECREMAINDER
! #elif  defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP SIMD
! #endif
! #elif defined __IBMBGQ__
!       !IBM* ALIGN(64,xp,yp,zp)
!       !IBM* ALIGN(64,uxp,uyp,uzp)
!       !IBM* ALIGN(64,gaminv)
!       !IBM* SIMD_LEVEL
! #endif
!     DO nn=ip,ip+MIN(lvect,np-ip+1)-1
!       const2 = gaminv(nn)*dtt
!       xp(nn) = xp(nn) + uxp(nn)*const2
!       yp(nn) = yp(nn) + uyp(nn)*const2
!       zp(nn) = zp(nn) + uzp(nn)*const2
! 
!     END DO
! #if defined __INTEL_COMPILER
! #elif defined _OPENMP && _OPENMP>=201307
! #ifndef NOVEC
!   !$OMP END SIMD
! #endif
! #endif

  ! End loop on particles
  ENDDO

RETURN
END SUBROUTINE field_gathering_plus_particle_pusher_3_3_3

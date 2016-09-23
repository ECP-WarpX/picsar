! ________________________________________________________________________________________
! PARTICLES_PUSH.F90
!
! Subroutines for the particle pusher in 3D
! ________________________________________________________________________________________


! ________________________________________________________________________________________
SUBROUTINE field_gathering_plus_particle_pusher
!
! Main subroutine for the field subroutine + particle pusher called 
! in the main loop (in submain.F90)
!
! ________________________________________________________________________________________
  USE fields
  USE shared_data
  USE params
  USE time_stat
  IMPLICIT NONE

#if defined(DEBUG)
  WRITE(0,*) "Field gathering + Push_particles: start"
#endif

  SELECT CASE (c_dim)
    CASE (2) ! 2D CASE X Z 

    ! Particle advance (one time step)
    CALL field_gathering_plus_particle_pusher_sub_2d(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
     nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

    CASE DEFAULT ! 3D CASE 
  
      ! The field gathering and the particle pusher are performed together
      IF (fg_p_pp_seperated.eq.0) THEN

        CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,&
        nxguards,nyguards,nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

      ELSE IF (fg_p_pp_seperated.eq.1) THEN

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
SUBROUTINE field_gathering_plus_particle_pusher_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
			     nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,  &
			     nozz,dxx,dyy,dzz,dtt,l_lower_order_in_v_in)
! Particle pusher in 3D called by the main function push_particle
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
  LOGICAL                  :: l_lower_order_in_v_in
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
  LOGICAL(idp)             :: isgathered=.FALSE.

    tdeb=MPI_WTIME()

#if PROFILING==3               
  CALL start_collection()      
#endif                         

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
!$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,&
!$OMP bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz,c_dim,l_lower_order_in_v_in) &
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
            CALL geteb2dxz_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,     			&
											  curr_tile%part_z, curr_tile%part_ex,                            	&
											  curr_tile%part_ey,curr_tile%part_ez,                   			&
											  curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz, 			&
											  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,              &
											  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,   &
											  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,        &
											  nzjg,noxx,noyy,nozz,currg%extile,currg%eytile, 					&
											  currg%eztile,                                          			&
											  currg%bxtile,currg%bytile,currg%bztile                 			&
											  ,.FALSE.,l_lower_order_in_v_in)										  
										  										  		
					CASE DEFAULT ! 3D CASE
					 
						!!! --- Gather electric and magnetic fields on particles						
						CALL geteb3d_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,     			&
											  curr_tile%part_z, curr_tile%part_ex,                            	&
											  curr_tile%part_ey,curr_tile%part_ez,                   			&
											  curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz, 			&
											  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,              &
											  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,   &
											  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,        &
											  nzjg,noxx,noyy,nozz,currg%extile,currg%eytile, 					&
											  currg%eztile,                                          			&
											  currg%bxtile,currg%bytile,currg%bztile                 			&
											  ,.FALSE.,l_lower_order_in_v_in)
					END SELECT
					
					!! --- Push velocity with E half step
					CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                    &
					curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey, 					    &
					curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
					!! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
					curr_tile%part_uz, curr_tile%part_gaminv)
					!! --- Push velocity with B half step
					CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                    &
					curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by,  &
					curr_tile%part_bz, curr%charge,curr%mass,dtt)
          !!! --- Push velocity with E half step
          CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                       &
          curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,                        &
          curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
          !! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
					curr_tile%part_uz, curr_tile%part_gaminv)
          !!!! --- push particle species positions a time step
          CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y, &
          curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy, &
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
SUBROUTINE field_gathering_plus_particle_pusher_cacheblock_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
			nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt,l_lower_order_in_v_in)
! Particle pusher in 3D called by the main function push_particle
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
  LOGICAL                  :: l_lower_order_in_v_in
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
  LOGICAL(idp)             :: isgathered=.FALSE.

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==3               
  CALL start_vtune_collection()      
#endif                         

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
                curr%charge,curr%mass,lvec_fieldgathe,LOGICAL(l_lower_order_in_v,isp))

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
                curr%charge,curr%mass,lvec_fieldgathe,LOGICAL(l_lower_order_in_v,isp))
                
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
                curr%charge,curr%mass,lvec_fieldgathe,LOGICAL(l_lower_order_in_v,isp))
          
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
SUBROUTINE particle_pusher_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
			nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,&
			noxx,noyy,nozz,dxx,dyy,dzz,dtt,l_lower_order_in_v_in)
!			
! Particle pusher in 3D called by the main function push_particle
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
  LOGICAL                  :: l_lower_order_in_v_in
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
  LOGICAL(idp)             :: isgathered=.FALSE.


  tdeb=MPI_WTIME()

#if VTUNE==3               
  CALL start_vtune_collection()      
#endif                         

#if defined(DEBUG)
  WRITE(0,*) "Push_particles: start"
#endif

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
!$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz,c_dim) &
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
					
					!! --- Push velocity with E half step
					CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                    &
					curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey, 					    &
					curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
					!! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
					curr_tile%part_uz, curr_tile%part_gaminv)
					!! --- Push velocity with B half step
					CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                    &
					curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by,  &
					curr_tile%part_bz, curr%charge,curr%mass,dtt)
                    !!! --- Push velocity with E half step
                    CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                       &
                    curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,                        &
                    curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
                    !! --- Set gamma of particles 
					          CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
					curr_tile%part_uz, curr_tile%part_gaminv)
                  !!!! --- push particle species positions a time step
                    CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y, &
                    curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy, &
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



!===============================================================================
!>  Field gathering+ (E & B) Push half a time step 
SUBROUTINE pxrpush_particles_part1
! ________________________________________________________________________________________
USE fields
USE shared_data
USE params
IMPLICIT NONE

CALL pxrpush_particles_part1_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
END SUBROUTINE pxrpush_particles_part1

! ________________________________________________________________________________________
!> Perform the field gathering + (E & B) Push half a time step 
!> @brief
!
!> This subroutine is called in pxrpush_particles_part1()
!> @details
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
SUBROUTINE pxrpush_particles_part1_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
			nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt)
! ________________________________________________________________________________________
	USE particles
	USE constants
	USE tiling
	USE timing
	IMPLICIT NONE
	INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard
	INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
	REAL(num), INTENT(IN) :: exg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
	REAL(num), INTENT(IN) :: eyg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
	REAL(num), INTENT(IN) :: ezg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
	REAL(num), INTENT(IN) :: bxg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
	REAL(num), INTENT(IN) :: byg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
	REAL(num), INTENT(IN) :: bzg(-nxguard:nxx+nxguard,-nyguard:nyy+nyguard,-nzguard:nzz+nzguard)
	REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
	INTEGER(idp) :: ispecies, ix, iy, iz, count
	INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
	TYPE(particle_species), POINTER :: curr
	TYPE(particle_tile), POINTER :: curr_tile
	TYPE(grid_tile), POINTER :: currg
	REAL(num) :: tdeb, tend
	INTEGER(idp) :: nxc, nyc, nzc, ipmin,ipmax, np,ip
	INTEGER(idp) :: nxjg,nyjg,nzjg
	LOGICAL(idp) :: isgathered=.FALSE.


!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
!$OMP nxjguard,nyjguard,nzjguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz,c_dim) &
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
            CALL geteb2dxz_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,     			&
											  curr_tile%part_z, curr_tile%part_ex,                            	&
											  curr_tile%part_ey,curr_tile%part_ez,                   			&
											  curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz, 			&
											  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,              &
											  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,   &
											  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,        &
											  nzjg,noxx,noyy,nozz,currg%extile,currg%eytile, 					&
											  currg%eztile,                                          			&
											  currg%bxtile,currg%bytile,currg%bztile                 			&
											  ,.FALSE.,.TRUE.)		
					CASE DEFAULT ! 3D CASE 
						!!! --- Gather electric and magnetic fields on particles
						CALL geteb3d_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,     		   	&
											  curr_tile%part_z, curr_tile%part_ex,                      	   	&
											  curr_tile%part_ey,curr_tile%part_ez,                   		   	&
											  curr_tile%part_bx, curr_tile%part_by,curr_tile%part_bz,          	&
											  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,             	&
											  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,  	&
											  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,       	&
											  nzjg,noxx,noyy,nozz,currg%extile,currg%eytile, 				   	&
											  currg%eztile,                                          			&
											  currg%bxtile,currg%bytile,currg%bztile                 			&
											  ,.FALSE.,.TRUE.)	
					END SELECT 				

					!! --- Push velocity with E half step
					CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                    &
					curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey, 					    &
					curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
					!! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
					curr_tile%part_uz, curr_tile%part_gaminv)
					!! --- Push velocity with B half step
					CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                 	&
					curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by,  &
					curr_tile%part_bz, curr%charge,curr%mass,dtt*0.5_num)
					!END DO
				END DO! END LOOP ON SPECIES
			ENDIF
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO

END SUBROUTINE pxrpush_particles_part1_sub


!===============================================================================
!  (B & E) Push half a time step + XYZ push half a time step 
!===============================================================================
SUBROUTINE pxrpush_particles_part2
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
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, ipmin,ipmax, np,ip

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$OMP nxjguards,nyjguards,nzjguards,ex,ey,ez,bx,by,bz,dx,dy,dz,dt,c_dim) &
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
                !!! ---- Loop by blocks over particles in a tile (blocking)
                !!! --- Push velocity with B half step
				CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                      &
				curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by, 	  &
				curr_tile%part_bz, curr%charge,curr%mass,dt*0.5_num)
				!! --- Push velocity with E half step
				CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                      &
				curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,				          &
				curr_tile%part_ez, curr%charge,curr%mass,dt*0.5_num)
				!! --- Sets gamma of particles
				CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                    &
				curr_tile%part_uz, curr_tile%part_gaminv)
				SELECT CASE (c_dim)
				CASE (2) ! 2D CASE 
					!! --- Advance particle position of one time step 
					CALL pxr_pushxz(count,curr_tile%part_x,                        				  &
					curr_tile%part_z, curr_tile%part_ux,   				          				  & 
					curr_tile%part_uz,curr_tile%part_gaminv,dt)				
				CASE DEFAULT ! 3D CASE 
					!! --- Advance particle position of one time step 
					CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y,                     &
					curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy,   				      & 
					curr_tile%part_uz,curr_tile%part_gaminv,dt)
				END SELECT 
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)
END SUBROUTINE pxrpush_particles_part2



!===============================================================================
!  Advance particle positions
SUBROUTINE pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
!===============================================================================
  USE constants
  USE omp_lib
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num) :: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np), gaminv(np)
  REAL(num) :: dt
  INTEGER(idp)  :: ip

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


!===============================================================================
!  Push the particle velocity with E field
SUBROUTINE pxr_epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
!===============================================================================

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

!===============================================================================
!  Push the particle velocity with B field (Boris algorithm)
! --- fast b-field rotation algorithm
SUBROUTINE pxr_bpush_v(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,dt)
!===============================================================================

USE constants
IMPLICIT NONE
INTEGER(idp)   :: np
REAL(num) :: uxp(np), uyp(np), uzp(np), gaminv(np)
REAL(num) :: bx(np), by(np), bz(np)
REAL(num) :: q,m,dt
INTEGER(idp)   :: ip
REAL(num) :: const,sx,sy,sz,tx,ty,tz,tsqi,uxppr,uyppr,uzppr

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


!===============================================================================
!  Push the particle velocity with B field (Boris algorithm)
! --- fast b-field rotation algorithm
SUBROUTINE pxr_set_gamma(np,uxp,uyp,uzp,gaminv)
!===============================================================================

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


! Push the particle velocity with E and B fields, assuming Vmid = 0.5*(Vold+Vnew),
! solving directly for the new gamma.
! This offers better cancellation of E+VxB than the Boris velocity push.
! Question: should we recompute gamma from the new u, in order to prevent roundoff errors
! to create mismatched values of u and gamma?
SUBROUTINE pxr_ebcancelpush3d(np,uxp,uyp,uzp,exp,eyp,ezp,bxp,byp,bzp,q,m,dt,which)
USE constants
INTEGER(idp) :: np,which
REAL(num)    :: uxp(np),uyp(np),uzp(np)
REAL(num)    :: exp(np),eyp(np),ezp(np),bxp(np),byp(np),bzp(np)
REAL(num)    :: q,m,dt
INTEGER(idp) :: ip
REAL(num)    :: const,bconst,s,gi,gisq,invclight,invclightsq,gprsq
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
		! 		--- get gi
		gi = 1_num/sqrt(1_num+(uxp(ip)*uxp(ip)+uyp(ip)*uyp(ip)+uzp(ip)*uzp(ip))*invclightsq)

		!       --- get tau
		taux = bconst*bxp(ip)
		tauy = bconst*byp(ip)
		tauz = bconst*bzp(ip)
		tausq = taux*taux+tauy*tauy+tauz*tauz
		!       --- get U',gamma'^2
		uxpr = uxp(ip) + const*exp(ip) + (uyp(ip)*tauz-uzp(ip)*tauy)*gi
		uypr = uyp(ip) + const*eyp(ip) + (uzp(ip)*taux-uxp(ip)*tauz)*gi
		uzpr = uzp(ip) + const*ezp(ip) + (uxp(ip)*tauy-uyp(ip)*taux)*gi
		gprsq = (1_num+(uxpr*uxpr+uypr*uypr+uzpr*uzpr)*invclightsq)
		!       --- get u*
		ust = (uxpr*taux+uypr*tauy+uzpr*tauz)*invclight
		!       --- get new gamma
		sigma = gprsq-tausq
		gisq = 2_num/(sigma+sqrt(sigma*sigma+4_num*(tausq+ust*ust)))
		gi = sqrt(gisq)
		!       --- get t,s
		bg = bconst*gi
		tx = bg*bxp(ip)
		ty = bg*byp(ip)
		tz = bg*bzp(ip)
		s = 1_num/(1_num+tausq*gisq)
		!       --- get t.u'
		tu = tx*uxpr+ty*uypr+tz*uzpr
		!       --- get new U
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
		! 		--- get gi
		gi = 1_num/sqrt(1_num+(uxp(ip)*uxp(ip)+uyp(ip)*uyp(ip)+uzp(ip)*uzp(ip))*invclightsq)
		!     --- get new U
		vx = uxp(ip)*gi
		vy = uyp(ip)*gi
		vz = uzp(ip)*gi
		uxp(ip) = uxp(ip) + const*( exp(ip) + vy*bzp(ip)-vz*byp(ip) )
		uyp(ip) = uyp(ip) + const*( eyp(ip) + vz*bxp(ip)-vx*bzp(ip) )
		uzp(ip) = uzp(ip) + const*( ezp(ip) + vx*byp(ip)-vy*bxp(ip) )
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
		gisq = 2_num/(sigma+sqrt(sigma*sigma+4_num*(tausq+ust*ust)))
		gi = sqrt(gisq)
		!       --- get t,s
		bg = bconst*gi
		tx = bg*bxp(ip)
		ty = bg*byp(ip)
		tz = bg*bzp(ip)
		s = 1_num/(1_num+tausq*gisq)
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
!> This function combined the field gathering and the particle pusher 
!> in 3D for CIC (order 1) shape factor.
!> The field gathering and the particle pusher are done in the same particle loop.
!> @brief
!
!> This function is vectorized.
!> @detail
! 
! Input parameters:
! - np: number of particles
! - xp,yp,zp: particle position
! - uxp,uyp,uzp: particle momentum
! - gaminv: inverse of the particle Lorentz factor
! - ex,ey,ez: particle electric field
! - bx,by,bz: particle magnetic field
! - xmin,ymin,zmin: tile minimum grid position
! - dx,dy,dz: space step
! - dtt: time step
! - nx,ny,nz: number of grid points in each direction
! - nxguard, nyguard, nzguard: number of guard cells in each direction 
! - exg,eyg,ezg: electric field grid
! - bxg,byg,bzg: magnetic field grid
! - lvect: vector size for cache blocking
! - l_lower_order_in_v: 
SUBROUTINE field_gathering_plus_particle_pusher_1_1_1(np,xp,yp,zp,uxp,uyp,uzp,gaminv, &
                                      ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,dtt,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,q,m,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  
  ! ___ Parameter declaration ____________________________________
  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num)                            :: q,m
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv
  LOGICAL(isp)                         :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz,dtt

  INTEGER(isp), DIMENSION(lvect)       :: j, k, l
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j,k,l
#endif  
  INTEGER(isp), DIMENSION(lvect)       :: j0, k0, l0
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j0,k0,l0
#endif    
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn,n
  INTEGER(isp)                         :: jj, kk, ll
  REAL(num)                            :: dxi, dyi, dzi
  REAL(num)                            :: x,y,z  
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: clghtisq,const1,const2,usq
  REAL(num)                            :: tx,ty,tz,tsqi
  REAL(num)                            :: wx,wy,wz    
  REAL(num)                            :: uxppr,uyppr,uzppr
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

  ! ______________________________________________________________________________________
  ! Loop on block of particles of size lvect
  DO ip=1,np,lvect
  
    ! ____________________________________________________________________________________
    ! Field gathering  
  
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
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
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1
      
			x = (xp(nn)-xmin)*dxi
			y = (yp(nn)-ymin)*dyi
			z = (zp(nn)-zmin)*dzi
	
			! Compute index of particle
			j(n)=floor(x)
			j0(n)=floor(x)
			k(n)=floor(y)
			k0(n)=floor(y)
			l(n)=floor(z)
			l0(n)=floor(z)
	
			xint=x-j(n)
			yint=y-k(n)
			zint=z-l(n)
	
			! Compute shape factors
			sx(n,0) = 1.0_num-xint
			sx(n,1) = xint
			sy(n,0) = 1.0_num-yint
			sy(n,1) = yint
			sz(n,0) = 1.0_num-zint
			sz(n,1) = zint
	
			xint=x-0.5_num-j0(n)
			yint=y-0.5_num-k0(n)
			zint=z-0.5_num-l0(n)
	
			sx0(n,0) = 1.0_num
			sy0(n,0) = 1.0_num
			sz0(n,0) = 1.0_num
			
		ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif


#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif 
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1 
      
			! Compute Ex on particle
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0(n),k(n),l(n))
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0(n),k(n)+1,l(n))
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0(n),k(n),l(n)+1)
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0(n),k(n)+1,l(n)+1)

			! Compute Ey on particle
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j(n),k0(n),l(n))
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j(n)+1,k0(n),l(n))
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j(n),k0(n),l(n)+1)
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j(n)+1,k0(n),l(n)+1)

			! Compute Ez on particle
			ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j(n),k(n),l0(n))
			ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j(n)+1,k(n),l0(n))
			ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j(n),k(n)+1,l0(n))
			ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j(n)+1,k(n)+1,l0(n))

    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1    

      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j(n),k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j(n)+1,k0(n),l0(n))
    
      ! Compute By on particle
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0(n),k(n),l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0(n),k(n)+1,l0(n))
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0(n),k0(n),l(n))
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0(n),k0(n),l(n)+1) 
      
    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ____________________________________________________________________________________
    ! Particle pusher
    
    ! ___ Push with E + gamma ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64      
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      uxp(nn) = uxp(nn) + ex(nn)*const1
      uyp(nn) = uyp(nn) + ey(nn)*const1
      uzp(nn) = uzp(nn) + ez(nn)*const1 
      
      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)      

    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Push with B ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64 
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64 
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      const2 = gaminv(nn)*const1
      tx = bx(nn)*const2
      ty = by(nn)*const2
      tz = bz(nn)*const2
      tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
      wx = tx*tsqi
      wy = ty*tsqi
      wz = tz*tsqi
      uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
      uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
      uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
      uxp(nn) = uxp(nn) + uyppr*wz - uzppr*wy
      uyp(nn) = uyp(nn) + uzppr*wx - uxppr*wz
      uzp(nn) = uzp(nn) + uxppr*wy - uyppr*wx
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Push with E + gamma ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64      
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      uxp(nn) = uxp(nn) + ex(nn)*const1
      uyp(nn) = uyp(nn) + ey(nn)*const1
      uzp(nn) = uzp(nn) + ez(nn)*const1 

      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq) 
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Update position ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      const2 = gaminv(nn)*dtt
      xp(nn) = xp(nn) + uxp(nn)*const2
      yp(nn) = yp(nn) + uyp(nn)*const2
      zp(nn) = zp(nn) + uzp(nn)*const2
        
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

  ! End loop on particles    
  ENDDO

RETURN
END SUBROUTINE field_gathering_plus_particle_pusher_1_1_1


! ________________________________________________________________________________________
!
! This function combined the field gathering and the particle pusher 
! in 3D for order 2 shape factor.
! The field gathering and the particle pusher are done in the same particle loop.
! This function is vectorized.
! 
! Input parameters:
! - np: number of particles
! - xp,yp,zp: particle position
! - uxp,uyp,uzp: particle momentum
! - gaminv: inverse of the particle Lorentz factor
! - ex,ey,ez: particle electric field
! - bx,by,bz: particle magnetic field
! - xmin,ymin,zmin: tile minimum grid position
! - dx,dy,dz: space step
! - dtt: time step
! - nx,ny,nz: number of grid points in each direction
! - nxguard, nyguard, nzguard: number of guard cells in each direction 
! - exg,eyg,ezg: electric field grid
! - bxg,byg,bzg: magnetic field grid
! - lvect: vector size for cache blocking
! - l_lower_order_in_v: 
SUBROUTINE field_gathering_plus_particle_pusher_2_2_2(np,xp,yp,zp,uxp,uyp,uzp,gaminv, &
                                      ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,dtt,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,q,m,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  
  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num)                            :: q,m
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv
  LOGICAL(isp)                         :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz,dtt
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn,n
  INTEGER(isp)                         :: jj, kk, ll
  INTEGER(isp), DIMENSION(lvect)       :: j, k, l
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j,k,l
#endif  
  INTEGER(isp), DIMENSION(lvect)       :: j0, k0, l0
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j0,k0,l0
#endif
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  REAL(num)                            :: clghtisq,const1,const2,usq
  REAL(num)                            :: tx,ty,tz,tsqi
  REAL(num)                            :: wx,wy,wz    
  REAL(num)                            :: uxppr,uyppr,uzppr
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

    ! ____________________________________________________________________________________
    ! Field gathering  
  
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
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
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1
      
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
    
			! Compute index of particle
			j(n)=nint(x)
			j0(n)=floor(x-0.5_num)
			k(n)=nint(y)
			k0(n)=floor(y-0.5_num)
			l(n)=nint(z)
			l0(n)=floor(z-0.5_num)
    
      xint=x-j(n)
      yint=y-k(n)
      zint=z-l(n)
      
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
		
			xint=x-0.5_num-j0(n)
			yint=y-0.5_num-k0(n)
			zint=z-0.5_num-l0(n)
		
			sx0(n, 0) = 1.0_num-xint
			sx0(n, 1) = xint
		
			sy0(n, 0) = 1.0_num-yint
			sy0(n, 1) = yint
		
			sz0(n, 0) = 1.0_num-zint
			sz0(n, 1) = zint
			
		ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif


#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ex:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1 
      
			! Compute Ex on particle
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0(n),k(n)-1,l(n)-1)
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0(n)+1,k(n)-1,l(n)-1)
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0(n),k(n),l(n)-1)
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0(n)+1,k(n),l(n)-1)
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0(n),k(n)+1,l(n)-1)
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0(n)+1,k(n)+1,l(n)-1)
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0(n),k(n)-1,l(n))
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0(n)+1,k(n)-1,l(n))
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0(n),k(n),l(n))
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0(n)+1,k(n),l(n))
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0(n),k(n)+1,l(n))
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0(n)+1,k(n)+1,l(n))
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0(n),k(n)-1,l(n)+1)
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0(n)+1,k(n)-1,l(n)+1)
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0(n),k(n),l(n)+1)
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0(n)+1,k(n),l(n)+1)
			ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0(n),k(n)+1,l(n)+1)
			ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0(n)+1,k(n)+1,l(n)+1)
	ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif 
#endif
	
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ey:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1 
    
			! Compute Ey on particle
			ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j(n)-1,k0(n),l(n)-1)
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j(n),k0(n),l(n)-1)
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j(n)+1,k0(n),l(n)-1)
			ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j(n)-1,k0(n)+1,l(n)-1)
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j(n),k0(n)+1,l(n)-1)
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j(n)+1,k0(n)+1,l(n)-1)
			ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j(n)-1,k0(n),l(n))
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j(n),k0(n),l(n))
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j(n)+1,k0(n),l(n))
			ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j(n)-1,k0(n)+1,l(n))
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j(n),k0(n)+1,l(n))
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j(n)+1,k0(n)+1,l(n))
			ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j(n)-1,k0(n),l(n)+1)
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j(n),k0(n),l(n)+1)
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j(n)+1,k0(n),l(n)+1)
			ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j(n)-1,k0(n)+1,l(n)+1)
			ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j(n),k0(n)+1,l(n)+1)
			ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j(n)+1,k0(n)+1,l(n)+1)
	ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif 
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ez:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1     
			! Compute Ez on particle
			ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j(n)-1,k(n)-1,l0(n))
			ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j(n),k(n)-1,l0(n))
			ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j(n)+1,k(n)-1,l0(n))
			ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j(n)-1,k(n),l0(n))
			ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j(n),k(n),l0(n))
			ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j(n)+1,k(n),l0(n))
			ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j(n)-1,k(n)+1,l0(n))
			ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j(n),k(n)+1,l0(n))
			ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j(n)+1,k(n)+1,l0(n))
			ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j(n)-1,k(n)-1,l0(n)+1)
			ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j(n),k(n)-1,l0(n)+1)
			ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j(n)+1,k(n)-1,l0(n)+1)
			ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j(n)-1,k(n),l0(n)+1)
			ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j(n),k(n),l0(n)+1)
			ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j(n)+1,k(n),l0(n)+1)
			ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j(n)-1,k(n)+1,l0(n)+1)
			ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j(n),k(n)+1,l0(n)+1)
			ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j(n)+1,k(n)+1,l0(n)+1)

    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED bx:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1    

			! Compute Bx on particle
			bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j(n)-1,k0(n),l0(n))
			bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j(n),k0(n),l0(n))
			bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j(n)+1,k0(n),l0(n))
			bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j(n)-1,k0(n)+1,l0(n))
			bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j(n),k0(n)+1,l0(n))
			bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j(n)+1,k0(n)+1,l0(n))
			bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j(n)-1,k0(n),l0(n)+1)
			bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j(n),k0(n),l0(n)+1)
			bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j(n)+1,k0(n),l0(n)+1)
			bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j(n)-1,k0(n)+1,l0(n)+1)
			bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j(n),k0(n)+1,l0(n)+1)
			bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j(n)+1,k0(n)+1,l0(n)+1)

    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED by:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1  
    
			! Compute By on particle
			by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0(n),k(n)-1,l0(n))
			by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0(n)+1,k(n)-1,l0(n))
			by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0(n),k(n),l0(n))
			by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0(n)+1,k(n),l0(n))
			by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0(n),k(n)+1,l0(n))
			by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0(n)+1,k(n)+1,l0(n))
			by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0(n),k(n)-1,l0(n)+1)
			by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0(n)+1,k(n)-1,l0(n)+1)
			by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0(n),k(n),l0(n)+1)
			by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0(n)+1,k(n),l0(n)+1)
			by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0(n),k(n)+1,l0(n)+1)
			by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0(n)+1,k(n)+1,l0(n)+1)
	ENDDO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1  
    
			! Compute Bz on particle
			bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0(n),k0(n),l(n)-1)
			bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0(n)+1,k0(n),l(n)-1)
			bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0(n),k0(n)+1,l(n)-1)
			bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0(n)+1,k0(n)+1,l(n)-1)
			bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0(n),k0(n),l(n))
			bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0(n)+1,k0(n),l(n))
			bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0(n),k0(n)+1,l(n))
			bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0(n)+1,k0(n)+1,l(n))
			bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0(n),k0(n),l(n)+1)
			bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0(n)+1,k0(n),l(n)+1)
			bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0(n),k0(n)+1,l(n)+1)
			bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0(n)+1,k0(n)+1,l(n)+1)
      
    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ____________________________________________________________________________________
    ! Particle pusher
    
    ! ___ Push with E + gamma ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64      
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      uxp(nn) = uxp(nn) + ex(nn)*const1
      uyp(nn) = uyp(nn) + ey(nn)*const1
      uzp(nn) = uzp(nn) + ez(nn)*const1 
      
      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)      

    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Push with B ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64 
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64 
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      const2 = gaminv(nn)*const1
      tx = bx(nn)*const2
      ty = by(nn)*const2
      tz = bz(nn)*const2
      tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
      wx = tx*tsqi
      wy = ty*tsqi
      wz = tz*tsqi
      uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
      uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
      uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
      uxp(nn) = uxp(nn) + uyppr*wz - uzppr*wy
      uyp(nn) = uyp(nn) + uzppr*wx - uxppr*wz
      uzp(nn) = uzp(nn) + uxppr*wy - uyppr*wx
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Push with E + gamma ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64      
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      uxp(nn) = uxp(nn) + ex(nn)*const1
      uyp(nn) = uyp(nn) + ey(nn)*const1
      uzp(nn) = uzp(nn) + ez(nn)*const1 

      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq) 
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Update position ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      const2 = gaminv(nn)*dtt
      xp(nn) = xp(nn) + uxp(nn)*const2
      yp(nn) = yp(nn) + uyp(nn)*const2
      zp(nn) = zp(nn) + uzp(nn)*const2
        
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

  ! End loop on particles    
  ENDDO

RETURN

END SUBROUTINE field_gathering_plus_particle_pusher_2_2_2


! ________________________________________________________________________________________
!> This function combined the field gathering and the particle pusher 
!> in 3D for CIC (order 1) shape factor.
!> The field gathering and the particle pusher are done in the same particle loop.
!> @brief
! 
!> This function is vectorized.
!> @detail
! 
! Input parameters:
! - np: number of particles
! - xp,yp,zp: particle position
! - uxp,uyp,uzp: particle momentum
! - gaminv: inverse of the particle Lorentz factor
! - ex,ey,ez: particle electric field
! - bx,by,bz: particle magnetic field
! - xmin,ymin,zmin: tile minimum grid position
! - dx,dy,dz: space step
! - dtt: time step
! - nx,ny,nz: number of grid points in each direction
! - nxguard, nyguard, nzguard: number of guard cells in each direction 
! - exg,eyg,ezg: electric field grid
! - bxg,byg,bzg: magnetic field grid
! - lvect: vector size for cache blocking
! - l_lower_order_in_v: 
SUBROUTINE field_gathering_plus_particle_pusher_3_3_3(np,xp,yp,zp,uxp,uyp,uzp,gaminv, &
                                      ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,dtt,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,q,m,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  
  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num)                            :: q,m
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv
  LOGICAL(isp)                         :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz,dtt
  INTEGER(isp)                         :: ip
  INTEGER(isp)                         :: nn,n
  INTEGER(isp)                         :: jj, kk, ll
  INTEGER(isp), DIMENSION(lvect)       :: j, k, l
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j,k,l
#endif  
  INTEGER(isp), DIMENSION(lvect)       :: j0, k0, l0
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j0,k0,l0
#endif    
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  REAL(num)                            :: clghtisq,const1,const2,usq
  REAL(num)                            :: tx,ty,tz,tsqi
  REAL(num)                            :: wx,wy,wz    
  REAL(num)                            :: uxppr,uyppr,uzppr
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

    ! ____________________________________________________________________________________
    ! Field gathering  
  
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
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
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1
      
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
    
      ! Compute index of particle
      j(n)=floor(x)
      j0(n)=floor(x)
      k(n)=floor(y)
      k0(n)=floor(y)
      l(n)=floor(z)
      l0(n)=floor(z)
    
      xint=x-j(n)
      yint=y-k(n)
      zint=z-l(n)
    
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
    
      xint=x-0.5_num-j0(n)
      yint=y-0.5_num-k0(n)
      zint=z-0.5_num-l0(n)
    
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
			
		ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif


#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ex:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1 
      
      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,-1)*exg(j0(n)-1,k(n)-1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0(n),k(n)-1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0(n)+1,k(n)-1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,-1)*exg(j0(n)-1,k(n),l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0(n),k(n),l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0(n)+1,k(n),l(n)-1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,-1)*exg(j0(n)-1,k(n)+1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0(n),k(n)+1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0(n)+1,k(n)+1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,2)*sz(n,-1)*exg(j0(n)-1,k(n)+2,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,2)*sz(n,-1)*exg(j0(n),k(n)+2,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,2)*sz(n,-1)*exg(j0(n)+1,k(n)+2,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,0)*exg(j0(n)-1,k(n)-1,l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0(n),k(n)-1,l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0(n)+1,k(n)-1,l(n))
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,0)*exg(j0(n)-1,k(n),l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0(n),k(n),l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0(n)+1,k(n),l(n))
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,0)*exg(j0(n)-1,k(n)+1,l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0(n),k(n)+1,l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0(n)+1,k(n)+1,l(n))
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,2)*sz(n,0)*exg(j0(n)-1,k(n)+2,l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,2)*sz(n,0)*exg(j0(n),k(n)+2,l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,2)*sz(n,0)*exg(j0(n)+1,k(n)+2,l(n))
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,1)*exg(j0(n)-1,k(n)-1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0(n),k(n)-1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0(n)+1,k(n)-1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,1)*exg(j0(n)-1,k(n),l(n)+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0(n),k(n),l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0(n)+1,k(n),l(n)+1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,1)*exg(j0(n)-1,k(n)+1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0(n),k(n)+1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0(n)+1,k(n)+1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,2)*sz(n,1)*exg(j0(n)-1,k(n)+2,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,2)*sz(n,1)*exg(j0(n),k(n)+2,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,2)*sz(n,1)*exg(j0(n)+1,k(n)+2,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,2)*exg(j0(n)-1,k(n)-1,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,2)*exg(j0(n),k(n)-1,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,2)*exg(j0(n)+1,k(n)-1,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,2)*exg(j0(n)-1,k(n),l(n)+2)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,2)*exg(j0(n),k(n),l(n)+2)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,2)*exg(j0(n)+1,k(n),l(n)+2)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,2)*exg(j0(n)-1,k(n)+1,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,2)*exg(j0(n),k(n)+1,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,2)*exg(j0(n)+1,k(n)+1,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,-1)*sy(n,2)*sz(n,2)*exg(j0(n)-1,k(n)+2,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,2)*sz(n,2)*exg(j0(n),k(n)+2,l(n)+2)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,2)*sz(n,2)*exg(j0(n)+1,k(n)+2,l(n)+2)
      
	ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif 
#endif
	
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ey:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1 
    
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,-1)*eyg(j(n)-1,k0(n)-1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,-1)*eyg(j(n),k0(n)-1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,-1)*eyg(j(n)+1,k0(n)-1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,-1)*sz(n,-1)*eyg(j(n)+2,k0(n)-1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j(n)-1,k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j(n),k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j(n)+1,k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,0)*sz(n,-1)*eyg(j(n)+2,k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j(n)-1,k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j(n),k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j(n)+1,k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,1)*sz(n,-1)*eyg(j(n)+2,k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,0)*eyg(j(n)-1,k0(n)-1,l(n))
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,0)*eyg(j(n),k0(n)-1,l(n))
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,0)*eyg(j(n)+1,k0(n)-1,l(n))
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,-1)*sz(n,0)*eyg(j(n)+2,k0(n)-1,l(n))
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j(n)-1,k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j(n),k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j(n)+1,k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,0)*sz(n,0)*eyg(j(n)+2,k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j(n)-1,k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j(n),k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j(n)+1,k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,1)*sz(n,0)*eyg(j(n)+2,k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,1)*eyg(j(n)-1,k0(n)-1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,1)*eyg(j(n),k0(n)-1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,1)*eyg(j(n)+1,k0(n)-1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,-1)*sz(n,1)*eyg(j(n)+2,k0(n)-1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j(n)-1,k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j(n),k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j(n)+1,k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,0)*sz(n,1)*eyg(j(n)+2,k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j(n)-1,k0(n)+1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j(n),k0(n)+1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j(n)+1,k0(n)+1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,1)*sz(n,1)*eyg(j(n)+2,k0(n)+1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,2)*eyg(j(n)-1,k0(n)-1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,2)*eyg(j(n),k0(n)-1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,2)*eyg(j(n)+1,k0(n)-1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,-1)*sz(n,2)*eyg(j(n)+2,k0(n)-1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,2)*eyg(j(n)-1,k0(n),l(n)+2)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,2)*eyg(j(n),k0(n),l(n)+2)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,2)*eyg(j(n)+1,k0(n),l(n)+2)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,0)*sz(n,2)*eyg(j(n)+2,k0(n),l(n)+2)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,2)*eyg(j(n)-1,k0(n)+1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,2)*eyg(j(n),k0(n)+1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,2)*eyg(j(n)+1,k0(n)+1,l(n)+2)
      ey(nn) = ey(nn) + sx(n,2)*sy0(n,1)*sz(n,2)*eyg(j(n)+2,k0(n)+1,l(n)+2)
      
	ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif 
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ez:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1     
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,-1)*ezg(j(n)-1,k(n)-1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,-1)*ezg(j(n),k(n)-1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,-1)*ezg(j(n)+1,k(n)-1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,-1)*sz0(n,-1)*ezg(j(n)+2,k(n)-1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,-1)*ezg(j(n)-1,k(n),l0(n)-1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,-1)*ezg(j(n),k(n),l0(n)-1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,-1)*ezg(j(n)+1,k(n),l0(n)-1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,0)*sz0(n,-1)*ezg(j(n)+2,k(n),l0(n)-1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,-1)*ezg(j(n)-1,k(n)+1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,-1)*ezg(j(n),k(n)+1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,-1)*ezg(j(n)+1,k(n)+1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,1)*sz0(n,-1)*ezg(j(n)+2,k(n)+1,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,2)*sz0(n,-1)*ezg(j(n)-1,k(n)+2,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,2)*sz0(n,-1)*ezg(j(n),k(n)+2,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,2)*sz0(n,-1)*ezg(j(n)+1,k(n)+2,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,2)*sz0(n,-1)*ezg(j(n)+2,k(n)+2,l0(n)-1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j(n)-1,k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j(n),k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j(n)+1,k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,2)*sy(n,-1)*sz0(n,0)*ezg(j(n)+2,k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j(n)-1,k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j(n),k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j(n)+1,k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,2)*sy(n,0)*sz0(n,0)*ezg(j(n)+2,k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j(n)-1,k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j(n),k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j(n)+1,k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,2)*sy(n,1)*sz0(n,0)*ezg(j(n)+2,k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,2)*sz0(n,0)*ezg(j(n)-1,k(n)+2,l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,2)*sz0(n,0)*ezg(j(n),k(n)+2,l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,2)*sz0(n,0)*ezg(j(n)+1,k(n)+2,l0(n))
      ez(nn) = ez(nn) + sx(n,2)*sy(n,2)*sz0(n,0)*ezg(j(n)+2,k(n)+2,l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j(n)-1,k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j(n),k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j(n)+1,k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,-1)*sz0(n,1)*ezg(j(n)+2,k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j(n)-1,k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j(n),k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j(n)+1,k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,0)*sz0(n,1)*ezg(j(n)+2,k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j(n)-1,k(n)+1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j(n),k(n)+1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j(n)+1,k(n)+1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,1)*sz0(n,1)*ezg(j(n)+2,k(n)+1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,2)*sz0(n,1)*ezg(j(n)-1,k(n)+2,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,2)*sz0(n,1)*ezg(j(n),k(n)+2,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,2)*sz0(n,1)*ezg(j(n)+1,k(n)+2,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,2)*sy(n,2)*sz0(n,1)*ezg(j(n)+2,k(n)+2,l0(n)+1)

    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED bx:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1    

      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,-1)*bxg(j(n)-1,k0(n)-1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,-1)*bxg(j(n),k0(n)-1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,-1)*bxg(j(n)+1,k0(n)-1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,-1)*sz0(n,-1)*bxg(j(n)+2,k0(n)-1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,-1)*bxg(j(n)-1,k0(n),l0(n)-1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,-1)*bxg(j(n),k0(n),l0(n)-1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,-1)*bxg(j(n)+1,k0(n),l0(n)-1)
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,0)*sz0(n,-1)*bxg(j(n)+2,k0(n),l0(n)-1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,-1)*bxg(j(n)-1,k0(n)+1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,-1)*bxg(j(n),k0(n)+1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,-1)*bxg(j(n)+1,k0(n)+1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,1)*sz0(n,-1)*bxg(j(n)+2,k0(n)+1,l0(n)-1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,0)*bxg(j(n)-1,k0(n)-1,l0(n))
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,0)*bxg(j(n),k0(n)-1,l0(n))
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,0)*bxg(j(n)+1,k0(n)-1,l0(n))
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,-1)*sz0(n,0)*bxg(j(n)+2,k0(n)-1,l0(n))
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j(n)-1,k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j(n),k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j(n)+1,k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,0)*sz0(n,0)*bxg(j(n)+2,k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j(n)-1,k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j(n),k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j(n)+1,k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,1)*sz0(n,0)*bxg(j(n)+2,k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,1)*bxg(j(n)-1,k0(n)-1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,1)*bxg(j(n),k0(n)-1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,1)*bxg(j(n)+1,k0(n)-1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,-1)*sz0(n,1)*bxg(j(n)+2,k0(n)-1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j(n)-1,k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j(n),k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j(n)+1,k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,0)*sz0(n,1)*bxg(j(n)+2,k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j(n)-1,k0(n)+1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j(n),k0(n)+1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j(n)+1,k0(n)+1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,2)*sy0(n,1)*sz0(n,1)*bxg(j(n)+2,k0(n)+1,l0(n)+1)

    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED by:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1  
      ! Compute By on particle
      by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,-1)*byg(j0(n)-1,k(n)-1,l0(n)-1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,-1)*byg(j0(n),k(n)-1,l0(n)-1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,-1)*byg(j0(n)+1,k(n)-1,l0(n)-1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,-1)*byg(j0(n)-1,k(n),l0(n)-1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,-1)*byg(j0(n),k(n),l0(n)-1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,-1)*byg(j0(n)+1,k(n),l0(n)-1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,-1)*byg(j0(n)-1,k(n)+1,l0(n)-1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,-1)*byg(j0(n),k(n)+1,l0(n)-1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,-1)*byg(j0(n)+1,k(n)+1,l0(n)-1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,2)*sz0(n,-1)*byg(j0(n)-1,k(n)+2,l0(n)-1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,2)*sz0(n,-1)*byg(j0(n),k(n)+2,l0(n)-1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,2)*sz0(n,-1)*byg(j0(n)+1,k(n)+2,l0(n)-1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,0)*byg(j0(n)-1,k(n)-1,l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0(n),k(n)-1,l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0(n)+1,k(n)-1,l0(n))
      by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,0)*byg(j0(n)-1,k(n),l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0(n),k(n),l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0(n)+1,k(n),l0(n))
      by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,0)*byg(j0(n)-1,k(n)+1,l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0(n),k(n)+1,l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0(n)+1,k(n)+1,l0(n))
      by(nn) = by(nn) + sx0(n,-1)*sy(n,2)*sz0(n,0)*byg(j0(n)-1,k(n)+2,l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,2)*sz0(n,0)*byg(j0(n),k(n)+2,l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,2)*sz0(n,0)*byg(j0(n)+1,k(n)+2,l0(n))
      by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,1)*byg(j0(n)-1,k(n)-1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0(n),k(n)-1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0(n)+1,k(n)-1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,1)*byg(j0(n)-1,k(n),l0(n)+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0(n),k(n),l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0(n)+1,k(n),l0(n)+1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,1)*byg(j0(n)-1,k(n)+1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0(n),k(n)+1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0(n)+1,k(n)+1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,-1)*sy(n,2)*sz0(n,1)*byg(j0(n)-1,k(n)+2,l0(n)+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,2)*sz0(n,1)*byg(j0(n),k(n)+2,l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,2)*sz0(n,1)*byg(j0(n)+1,k(n)+2,l0(n)+1)
	ENDDO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED bz:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1  
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,-1)*bzg(j0(n)-1,k0(n)-1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,-1)*bzg(j0(n),k0(n)-1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,-1)*bzg(j0(n)+1,k0(n)-1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,-1)*bzg(j0(n)-1,k0(n),l(n)-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0(n),k0(n),l(n)-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0(n)+1,k0(n),l(n)-1)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,-1)*bzg(j0(n)-1,k0(n)+1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0(n),k0(n)+1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0(n)+1,k0(n)+1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,0)*bzg(j0(n)-1,k0(n)-1,l(n))
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,0)*bzg(j0(n),k0(n)-1,l(n))
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,0)*bzg(j0(n)+1,k0(n)-1,l(n))
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,0)*bzg(j0(n)-1,k0(n),l(n))
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0(n),k0(n),l(n))
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0(n)+1,k0(n),l(n))
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,0)*bzg(j0(n)-1,k0(n)+1,l(n))
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0(n),k0(n)+1,l(n))
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0(n)+1,k0(n)+1,l(n))
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,1)*bzg(j0(n)-1,k0(n)-1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,1)*bzg(j0(n),k0(n)-1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,1)*bzg(j0(n)+1,k0(n)-1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,1)*bzg(j0(n)-1,k0(n),l(n)+1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0(n),k0(n),l(n)+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0(n)+1,k0(n),l(n)+1)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,1)*bzg(j0(n)-1,k0(n)+1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0(n),k0(n)+1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0(n)+1,k0(n)+1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,2)*bzg(j0(n)-1,k0(n)-1,l(n)+2)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,2)*bzg(j0(n),k0(n)-1,l(n)+2)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,2)*bzg(j0(n)+1,k0(n)-1,l(n)+2)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,2)*bzg(j0(n)-1,k0(n),l(n)+2)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,2)*bzg(j0(n),k0(n),l(n)+2)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,2)*bzg(j0(n)+1,k0(n),l(n)+2)
      bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,2)*bzg(j0(n)-1,k0(n)+1,l(n)+2)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,2)*bzg(j0(n),k0(n)+1,l(n)+2)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,2)*bzg(j0(n)+1,k0(n)+1,l(n)+2)
      
    END DO
#if  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ____________________________________________________________________________________
    ! Particle pusher
    
    ! ___ Push with E + gamma ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64      
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      uxp(nn) = uxp(nn) + ex(nn)*const1
      uyp(nn) = uyp(nn) + ey(nn)*const1
      uzp(nn) = uzp(nn) + ez(nn)*const1 
      
      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)      

    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Push with B ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64 
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64 
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      const2 = gaminv(nn)*const1
      tx = bx(nn)*const2
      ty = by(nn)*const2
      tz = bz(nn)*const2
      tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
      wx = tx*tsqi
      wy = ty*tsqi
      wz = tz*tsqi
      uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
      uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
      uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
      uxp(nn) = uxp(nn) + uyppr*wz - uzppr*wy
      uyp(nn) = uyp(nn) + uzppr*wx - uxppr*wz
      uzp(nn) = uzp(nn) + uxppr*wy - uyppr*wx
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Push with E + gamma ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64      
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)      
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      uxp(nn) = uxp(nn) + ex(nn)*const1
      uyp(nn) = uyp(nn) + ey(nn)*const1
      uzp(nn) = uzp(nn) + ez(nn)*const1 

      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq) 
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

    ! ___ Update position ___
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED gaminv:64
      !DIR VECTOR NONTEMPORAL(xp,yp,zp,ex,ey,ez,bx,by,bz,uxp,uyp,uzp,gaminv)       
      !DIR$ SIMD VECREMAINDER 
#elif  defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP SIMD 
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,gaminv)
      !IBM* SIMD_LEVEL      
#endif 
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
      const2 = gaminv(nn)*dtt
      xp(nn) = xp(nn) + uxp(nn)*const2
      yp(nn) = yp(nn) + uyp(nn)*const2
      zp(nn) = zp(nn) + uzp(nn)*const2
        
    END DO
#if defined __INTEL_COMPILER 
#elif defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
	!$OMP END SIMD 
#endif
#endif

  ! End loop on particles    
  ENDDO

RETURN
END SUBROUTINE field_gathering_plus_particle_pusher_3_3_3

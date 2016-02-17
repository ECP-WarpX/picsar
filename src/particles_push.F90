
!===============================================================================
!  Advance particles a full time step
!===============================================================================

SUBROUTINE push_particles
USE fields
USE shared_data
USE params
IMPLICIT NONE

! Particle advance (one time step)
CALL push_particles_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
	 nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
	 
END SUBROUTINE push_particles


SUBROUTINE push_particles_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
			nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt)
USE particles
USE constants
USE tiling
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
TYPE(grid_tile), POINTER :: currg
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, ipmin,ipmax, np,ip
INTEGER(idp) :: nxjg,nyjg,nzjg
LOGICAL(idp) :: isgathered=.FALSE.

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles, &
!$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile, currg, count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,np,nxjg,nyjg,nzjg, isgathered)
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
					CALL pxr_gete3d_n_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,                                     &
										  curr_tile%part_z, curr_tile%part_ex,                                                       &
										  curr_tile%part_ey,curr_tile%part_ez,                   									 &
										  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                            			 &
										  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                  			 &
										  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,             				     &
										  nzjg,noxx,noyy,nozz,currg%extile,currg%eytile, 											 &
										  currg%eztile,.FALSE.,.TRUE.)
					!!! --- Gather magnetic fields on particles
					CALL pxr_getb3d_n_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,                                     &
									  curr_tile%part_z, curr_tile%part_bx,                    									     &
									  curr_tile%part_by,curr_tile%part_bz,                    										 &
									  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                              			 &
									  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                   			 &
									  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,                					 &
									  nzjg,noxx,noyy,nozz,currg%bxtile,currg%bytile,  	 											 &
									  currg%bztile,.FALSE.,.TRUE.)
					!! --- Push velocity with E half step
					CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                   &
					curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey, 					   &
					curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
					!! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                 &
					curr_tile%part_uz, curr_tile%part_gaminv)
					!! --- Push velocity with B half step
					CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                   &
					curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by, &
					curr_tile%part_bz, curr%charge,curr%mass,dtt)
                    !!! --- Push velocity with E half step
                    CALL pxr_epush_v(np,curr_tile%part_ux, curr_tile%part_uy,                       &
                    curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,                        &
                    curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
                    !! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                  &
					curr_tile%part_uz, curr_tile%part_gaminv)
                    !!!! --- push particle species positions a time step
                    CALL pxr_pushxyz(np,curr_tile%part_x,curr_tile%part_y,                          &
                    curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy,   				        &
                    curr_tile%part_uz,curr_tile%part_gaminv,dtt)
                END DO! END LOOP ON SPECIES
            ENDIF
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)
END SUBROUTINE push_particles_sub

!===============================================================================
!  Field gathering+ (E & B) Push half a time step 
!===============================================================================
SUBROUTINE pxrpush_particles_part1
USE fields
USE shared_data
USE params
IMPLICIT NONE

CALL pxrpush_particles_part1_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
END SUBROUTINE pxrpush_particles_part1


SUBROUTINE pxrpush_particles_part1_sub(exg,eyg,ezg,bxg,byg,bzg,nxx,nyy,nzz, &
			nxguard,nyguard,nzguard,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt)
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
!$OMP nxjguard,nyjguard,nzjguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,noxx,noyy,nozz) &
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
					!!! --- Gather electric field on particles
					CALL pxr_gete3d_n_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,                                     &
										  curr_tile%part_z, curr_tile%part_ex,                                                       &
										  curr_tile%part_ey,curr_tile%part_ez,                   									 &
										  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                            			 &
										  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                  			 &
										  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,             				     &
										  nzjg,noxx,noyy,nozz,currg%extile,currg%eytile, &
										  currg%eztile,.FALSE.,.TRUE.)
					!!! --- Gather magnetic fields on particles
					CALL pxr_getb3d_n_energy_conserving(count,curr_tile%part_x,curr_tile%part_y,                                     &
									  curr_tile%part_z, curr_tile%part_bx,                    									     &
									  curr_tile%part_by,curr_tile%part_bz,                    										 &
									  curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                              			 &
									  curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                   			 &
									  curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjg,nyjg,                					 &
									  nzjg,noxx,noyy,nozz,currg%bxtile,currg%bytile,  	 &
									  currg%bztile,.FALSE.,.TRUE.)
					!! --- Push velocity with E half step
					CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                 &
					curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey, 					 &
					curr_tile%part_ez, curr%charge,curr%mass,dtt*0.5_num)
					!! --- Set gamma of particles 
					CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                 &
					curr_tile%part_uz, curr_tile%part_gaminv)
					!! --- Push velocity with B half step
					CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                 &
					curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by,				     &
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
!$OMP nxjguards,nyjguards,nzjguards,ex,ey,ez,bx,by,bz,dx,dy,dz,dt) &
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
				CALL pxr_bpush_v(count,curr_tile%part_ux, curr_tile%part_uy,                 &
				curr_tile%part_uz,curr_tile%part_gaminv, curr_tile%part_bx, curr_tile%part_by, 				     &
				curr_tile%part_bz, curr%charge,curr%mass,dt*0.5_num)
				!! --- Push velocity with E half step
				CALL pxr_epush_v(count,curr_tile%part_ux, curr_tile%part_uy,                 &
				curr_tile%part_uz, curr_tile%part_ex, curr_tile%part_ey,				     &
				curr_tile%part_ez, curr%charge,curr%mass,dt*0.5_num)
				!! --- Sets gamma of particles
				CALL pxr_set_gamma(count,curr_tile%part_ux, curr_tile%part_uy,                 &
				curr_tile%part_uz, curr_tile%part_gaminv)
				!! --- Advance particle position of one time step 
				CALL pxr_pushxyz(count,curr_tile%part_x,curr_tile%part_y,                    &
				curr_tile%part_z, curr_tile%part_ux,curr_tile%part_uy,   				     &
				curr_tile%part_uz,curr_tile%part_gaminv,dt)
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

!!$OMP PARALLEL DO PRIVATE(ip)
DO ip=1,np
    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt
ENDDO
!!$OMP END PARALLEL DO

RETURN
END SUBROUTINE pxr_pushxyz

!===============================================================================
!  Push the particle velocity with E field
SUBROUTINE pxr_epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
!===============================================================================

USE constants
IMPLICIT NONE
INTEGER(idp) :: np
REAL(num):: uxp(np),uyp(np),uzp(np)
REAL(num):: ex(np),ey(np),ez(np)
REAL(num):: q,m,dt

INTEGER(idp) :: ip
REAL(num):: const

const = q*dt/m
!!$OMP PARALLEL DO PRIVATE(ip)
DO ip=1,np
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const
ENDDO
!!$OMP END PARALLEL DO

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

!!$OMP PARALLEL DO PRIVATE(ip, tx, ty, tz, tsqi, sx, sy, sz, uxppr, uyppr, uzppr)
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
!!$OMP END PARALLEL DO

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

!!$OMP PARALLEL DO PRIVATE(ip)
DO ip=1,np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
END DO
!!$OMP END PARALLEL DO

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
ELSE IF(which==1) THEN
	!     --- first half push
    const = 0.5_num*q*dt/m
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
ELSE IF(which==2) THEN
	!     --- second half push
    const = 0.5_num*q*dt/m
    bconst = const
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
ENDIF
RETURN
END SUBROUTINE pxr_ebcancelpush3d










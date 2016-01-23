SUBROUTINE push_particles
USE fields
USE shared_data
USE params
IMPLICIT NONE

CALL push_particles_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
	 nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
END SUBROUTINE push_particles


!===============================================================================
!  Advance particles a full time step
!===============================================================================
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
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, ipmin,ipmax, np,ip
INTEGER(idp) :: nblk=64

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$OMP nxjguard,nyjguard,nzjguard,nxguard,nyguard,nzguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,nblk,noxx,noyy,nozz) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,np)
DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
        DO ix=1, ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                ! - Get current tile properties
                ! - Init current tile variables
                curr=>species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min-nxjguard
                jmax=curr_tile%nx_tile_max+nxjguard
                kmin=curr_tile%ny_tile_min-nyjguard
                kmax=curr_tile%ny_tile_max+nyjguard
                lmin=curr_tile%nz_tile_min-nzjguard
                lmax=curr_tile%nz_tile_max+nzjguard
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                curr_tile%part_ex = 0.0_num
                curr_tile%part_ey = 0.0_num
                curr_tile%part_ez = 0.0_num
                curr_tile%part_bx=0.0_num
                curr_tile%part_by=0.0_num
                curr_tile%part_bz=0.0_num
                !!! ---- Loop by blocks over particles in a tile (blocking)
                DO ip=1,count,nblk
                    np=MIN(count-ip+1,nblk)
                    ipmin=ip
                    ipmax=ip+np-1
                    !!! --- Gather electric field on particles
                     CALL pxrgete3d_n_energy_conserving(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),                &
                                       curr_tile%part_z(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax),                                  &
                                       curr_tile%part_ey(ipmin:ipmax),curr_tile%part_ez(ipmin:ipmax),                                  &
                                       curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                                            &
                                       curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                                 &
				                       curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguard,nyjguard,              				   &
				    				   nzjguard,noxx,noyy,nozz,exg(jmin:jmax,kmin:kmax,lmin:lmax),eyg(jmin:jmax,kmin:kmax,lmin:lmax),  &
                                      ezg(jmin:jmax,kmin:kmax,lmin:lmax),.TRUE.)
                    !!! --- Gather magnetic fields on particles
                    CALL pxrgetb3d_n_energy_conserving(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),     			  &
                                      curr_tile%part_z(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax),                    			  &
				                      curr_tile%part_by(ipmin:ipmax),curr_tile%part_bz(ipmin:ipmax),                    			  &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                             				  &
                                      curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                   			  &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguard,nyjguard,               				  &
                                      nzjguard,noxx,noyy,nozz,bxg(jmin:jmax,kmin:kmax,lmin:lmax),byg(jmin:jmax,kmin:kmax,lmin:lmax),  &
                                      bzg(jmin:jmax,kmin:kmax,lmin:lmax), .TRUE.)
                    !! --- Push velocity with E half step
                    CALL pxr_epush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax), curr_tile%part_ey(ipmin:ipmax), &
                    curr_tile%part_ez(ipmin:ipmax), curr%charge,curr%mass,dtt*0.5_num)
                    !!! --- Push velocity with B half step
                    CALL pxr_bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax), &
                    curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dtt)
                    !!! --- Push velocity with E half step
                    CALL pxr_epush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax), curr_tile%part_ey(ipmin:ipmax), &
                    curr_tile%part_ez(ipmin:ipmax), curr%charge,curr%mass,dtt*0.5_num)
                    !!!! --- push particle species positions a time step
                    CALL pxr_pushxyz(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),                    &
                    curr_tile%part_z(ipmin:ipmax), curr_tile%part_ux(ipmin:ipmax),curr_tile%part_uy(ipmin:ipmax),   &
                    curr_tile%part_uz(ipmin:ipmax),dtt)
                END DO
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)
END SUBROUTINE push_particles_sub

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
INTEGER(idp) :: nblk=64

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$OMP nxjguards,nyjguards,nzjguards,ex,ey,ez,bx,by,bz,dx,dy,dz,dt,nblk) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,np)
DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
        DO ix=1, ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                ! - Get current tile properties
                ! - Init current tile variables
                curr=>species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min-nxjguards
                jmax=curr_tile%nx_tile_max+nxjguards
                kmin=curr_tile%ny_tile_min-nyjguards
                kmax=curr_tile%ny_tile_max+nyjguards
                lmin=curr_tile%nz_tile_min-nzjguards
                lmax=curr_tile%nz_tile_max+nzjguards
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                !!! ---- Loop by blocks over particles in a tile (blocking)
                DO ip=1,count,nblk
                    np=MIN(count-ip+1,nblk)
                    ipmin=ip
                    ipmax=ip+np-1
                    !!! --- Push velocity with B half step
                    CALL pxr_bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax), &
                    curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dt*0.5_num)
                    !! --- Push velocity with E half step
                    CALL pxr_epush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax), curr_tile%part_ey(ipmin:ipmax), &
                    curr_tile%part_ez(ipmin:ipmax), curr%charge,curr%mass,dt*0.5_num)
                    CALL pxr_pushxyz(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),                    &
                    curr_tile%part_z(ipmin:ipmax), curr_tile%part_ux(ipmin:ipmax),curr_tile%part_uy(ipmin:ipmax),   &
                    curr_tile%part_uz(ipmin:ipmax),dt)
                END DO
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)
END SUBROUTINE pxrpush_particles_part2


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
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, ipmin,ipmax, np,ip
INTEGER(idp) :: nblk=64


tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$OMP nxjguard,nyjguard,nzjguard,exg,eyg,ezg,bxg,byg,bzg,dxx,dyy,dzz,dtt,nblk,noxx,noyy,nozz) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,nxc,nyc,nzc, ipmin,ipmax,ip,np)
DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
        DO ix=1, ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                ! - Get current tile properties
                ! - Init current tile variables
                curr=>species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min-nxjguard
                jmax=curr_tile%nx_tile_max+nxjguard
                kmin=curr_tile%ny_tile_min-nyjguard
                kmax=curr_tile%ny_tile_max+nyjguard
                lmin=curr_tile%nz_tile_min-nzjguard
                lmax=curr_tile%nz_tile_max+nzjguard
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                curr_tile%part_ex = 0.0_num
                curr_tile%part_ey = 0.0_num
                curr_tile%part_ez = 0.0_num
                curr_tile%part_bx=0.0_num
                curr_tile%part_by=0.0_num
                curr_tile%part_bz=0.0_num
                !!! ---- Loop by blocks over particles in a tile (blocking)
                DO ip=1,count,nblk
                    np=MIN(count-ip+1,nblk)
                    ipmin=ip
                    ipmax=ip+np-1
                    !!! --- Gather electric field on particles
                    CALL pxrgete3d_n_energy_conserving(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),    &
                                      curr_tile%part_z(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax),                   &
                                      curr_tile%part_ey(ipmin:ipmax),curr_tile%part_ez(ipmin:ipmax),                   &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                             &
                                      curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                  &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguard,nyjguard,               &
									  nzjguard,noxx,noyy,nozz,exg(jmin:jmax,kmin:kmax,lmin:lmax),eyg(jmin:jmax,kmin:kmax,lmin:lmax),  &
							          ezg(jmin:jmax,kmin:kmax,lmin:lmax),.TRUE.)
                    !!! --- Gather magnetic fields on particles
                    CALL pxrgetb3d_n_energy_conserving(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),     &
                                      curr_tile%part_z(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax),                    &
                                      curr_tile%part_by(ipmin:ipmax),curr_tile%part_bz(ipmin:ipmax),                    &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                              &
                                      curr_tile%z_grid_tile_min, dxx,dyy,dzz,curr_tile%nx_cells_tile,                   &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguard,nyjguard,                &
									  nzjguard,noxx,noyy,nozz,bxg(jmin:jmax,kmin:kmax,lmin:lmax),byg(jmin:jmax,kmin:kmax,lmin:lmax),  &
									  bzg(jmin:jmax,kmin:kmax,lmin:lmax),.TRUE.)
                    CALL pxr_epush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax), curr_tile%part_ey(ipmin:ipmax), &
                    curr_tile%part_ez(ipmin:ipmax), curr%charge,curr%mass,dtt*0.5_num)
                    !!! --- Push velocity with B half step
                    CALL pxr_bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax), &
                    curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dtt*0.5_num)
                END DO
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)
END SUBROUTINE pxrpush_particles_part1_sub
!===============================================================================
!  Advance particle positions
SUBROUTINE pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,dt)
!===============================================================================
USE constants
USE omp_lib
IMPLICIT NONE
INTEGER(idp)   :: np
REAL(num) :: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np)
REAL(num) :: dt,gaminv,clghtisq,usq
INTEGER(idp)  :: ip

clghtisq = 1.0_num/clight**2
!!$OMP PARALLEL DO PRIVATE(ip, usq, gaminv)
DO ip=1,np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv = 1.0_num/sqrt(1.0_num + usq)
    xp(ip) = xp(ip) + uxp(ip)*gaminv*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv*dt
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
SUBROUTINE pxr_bpush_v(np,uxp,uyp,uzp,bx,by,bz,q,m,dt)
!===============================================================================

USE constants
IMPLICIT NONE
INTEGER(idp)   :: np
REAL(num) :: uxp(np), uyp(np), uzp(np)
REAL(num) :: bx(np), by(np), bz(np)
REAL(num) :: q,m,dt,gaminv
INTEGER(idp)   :: ip
REAL(num) :: const,clghtisq,sx,sy,sz,tx,ty,tz,tsqi,uxppr,uyppr,uzppr,usq

const = q*dt*0.5_num/m
clghtisq = 1.0_num/clight**2

!!$OMP PARALLEL DO PRIVATE(ip, tx, ty, tz, tsqi, sx, sy, sz, uxppr, uyppr, uzppr, usq, gaminv)
DO ip=1,np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv = 1.0_num/sqrt(1.0_num + usq)
    tx = gaminv*bx(ip)*const
    ty = gaminv*by(ip)*const
    tz = gaminv*bz(ip)*const
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












!===============================================================================
!  Advance particles a full time step
!===============================================================================
SUBROUTINE push_particles
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
INTEGER(idp) :: nblk=1000000

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
                    CALL gete3d_energy_conserving_1_1_1(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),&
                                      curr_tile%part_z(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax),                   &
                                      curr_tile%part_ey(ipmin:ipmax),curr_tile%part_ez(ipmin:ipmax),                   &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                             &
                                      curr_tile%z_grid_tile_min, dx,dy,dz,curr_tile%nx_cells_tile,                     &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,             &
                                      nzjguards,ex(jmin:jmax,kmin:kmax,lmin:lmax),ey(jmin:jmax,kmin:kmax,lmin:lmax),   &
                                      ez(jmin:jmax,kmin:kmax,lmin:lmax))
                    !!! --- Gather magnetic fields on particles
                    CALL getb3d_energy_conserving_1_1_1(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),&
                                      curr_tile%part_z(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax),                      &
                                      curr_tile%part_by(ipmin:ipmax),curr_tile%part_bz(ipmin:ipmax),                      &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                                &
                                      curr_tile%z_grid_tile_min, dx,dy,dz,curr_tile%nx_cells_tile,                        &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,                &
                                      nzjguards,bx(jmin:jmax,kmin:kmax,lmin:lmax),by(jmin:jmax,kmin:kmax,lmin:lmax),      &
                                      bz(jmin:jmax,kmin:kmax,lmin:lmax))
                    !!! --- Push velocity with B half step
                    CALL bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax), &
                    curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dt*0.5_num)
                    !!! --- Push velocity with E half step
                    CALL epush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax), curr_tile%part_ey(ipmin:ipmax), &
                    curr_tile%part_ez(ipmin:ipmax), curr%charge,curr%mass,dt)
                    !!! --- Push velocity with B half step
                    CALL bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),                 &
                    curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax), &
                    curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dt*0.5_num)
                    !!!! --- push particle species positions a time step
                    CALL pushxyz(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),                    &
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
END SUBROUTINE push_particles

!===============================================================================
!  Advance particle positions
SUBROUTINE pushxyz(np,xp,yp,zp,uxp,uyp,uzp,dt)
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
END SUBROUTINE pushxyz

!===============================================================================
!  Push the particle velocity with E field
SUBROUTINE epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
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
END SUBROUTINE epush_v

!===============================================================================
!  Push the particle velocity with B field (Boris algorithm)
! --- fast b-field rotation algorithm
SUBROUTINE bpush_v(np,uxp,uyp,uzp,bx,by,bz,q,m,dt)
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

END SUBROUTINE bpush_v

!===============================================================================
!  Advance velocity half a time step
!===============================================================================
SUBROUTINE push_particles_v
USE particles
USE constants
USE params
USE tiling
IMPLICIT NONE
INTEGER :: nblock=32
INTEGER ispecies, ix, iy, iz, count,subcount,ip,ipmax
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez, nspecies,species_parray, dt,nblock) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,ip,subcount,ipmax)
DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
        DO ix=1, ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=> curr%array_of_tiles(ix,iy,iz)
                count= curr_tile%np_tile
                DO ip=1,count,nblock
                    subcount=MIN(nblock,count-ip+1)
                    ipmax=subcount+ip
		    !! -- Merged routine
                    CALL merged_push_v(subcount,curr_tile%part_ux(ip:ipmax), curr_tile%part_uy(ip:ipmax), &
                    curr_tile%part_uz(ip:ipmax), curr_tile%part_bx(ip:ipmax), curr_tile%part_by(ip:ipmax), &
                    curr_tile%part_bz(ip:ipmax), curr_tile%part_ex(ip:ipmax), curr_tile%part_ey(ip:ipmax), &
		    curr_tile%part_ez(ip:ipmax), curr%charge,curr%mass,dt*0.5_num)
                    !! -- Push velocity with B half step
		    ! CALL bpush_v(subcount,curr_tile%part_ux(ip:ipmax), curr_tile%part_uy(ip:ipmax), &
                    ! curr_tile%part_uz(ip:ipmax), curr_tile%part_bx(ip:ipmax), curr_tile%part_by(ip:ipmax), &
                    ! curr_tile%part_bz(ip:ipmax), curr%charge,curr%mass,dt*0.5_num)
                    !! -- Push velocity with E half step
                    ! CALL epush_v(subcount,curr_tile%part_ux(ip:ipmax), curr_tile%part_uy(ip:ipmax), &
                    ! curr_tile%part_uz(ip:ipmax), curr_tile%part_ex(ip:ipmax), curr_tile%part_ey(ip:ipmax), &
                    ! curr_tile%part_ez(ip:ipmax), curr%charge,curr%mass,dt*0.5_num)
                END DO
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)

END SUBROUTINE push_particles_v

!===============================================================================
!  Advance particles a full time step
!===============================================================================
SUBROUTINE push_particles_xyz
USE particles
USE constants
USE params
USE tiling
IMPLICIT NONE
INTEGER ispecies, ix, iy, iz, count
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,dt) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count)
DO iz=1,ntilez !LOOP ON TILES
    DO iy=1,ntiley
        DO ix=1,ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=> curr%array_of_tiles(ix,iy,iz)
                count= curr_tile%np_tile
                !!!! --- push particle species positions a time step
                CALL pushxyz(count,curr_tile%part_x(1),curr_tile%part_y(1), &
                curr_tile%part_z(1), curr_tile%part_ux(1),curr_tile%part_uy(1),&
                curr_tile%part_uz(1),dt)
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)

END SUBROUTINE push_particles_xyz





!===============================================================================
!  Advance particle positions
SUBROUTINE pushxyz(np,xp,yp,zp,uxp,uyp,uzp,dt)
!===============================================================================
USE constants
USE omp_lib
IMPLICIT NONE
INTEGER   :: np
REAL(num) :: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np)
REAL(num) :: dt,gaminv,clghtisq,usq
INTEGER   :: ip

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
END

!===============================================================================
!  Push the particle velocity with E field
SUBROUTINE epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
!===============================================================================

USE constants
IMPLICIT NONE
INTEGER :: np
REAL(num):: uxp(np),uyp(np),uzp(np)
REAL(num):: ex(np),ey(np),ez(np)
REAL(num):: q,m,dt

INTEGER :: ip
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
INTEGER   :: np
REAL(num) :: uxp(np), uyp(np), uzp(np)
REAL(num) :: bx(np), by(np), bz(np)
REAL(num) :: q,m,dt,gaminv
INTEGER   :: ip
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

!===============================================================================
!  Merged push routine
SUBROUTINE merged_push_v(np,uxp,uyp,uzp,bx,by,bz,ex,ey,ez,q,m,dt)
!===============================================================================

USE constants
IMPLICIT NONE
INTEGER   :: np
REAL(num) :: uxp(np), uyp(np), uzp(np)
REAL(num) :: bx(np), by(np), bz(np), ex(np), ey(np), ez(np)
REAL(num) :: q,m,dt,gaminv
INTEGER   :: ip
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

    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const
ENDDO
!!$OMP END PARALLEL DO

RETURN

END SUBROUTINE merged_push_v

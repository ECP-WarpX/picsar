!===============================================================================
!  Advance particle positions
subroutine pushxyz(np,xp,yp,zp,uxp,uyp,uzp,dt)
!===============================================================================
use constants
use omp_lib
implicit none
integer:: np
real(num):: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np)
real(num):: dt,gaminv,clghtisq,usq
integer:: ip

clghtisq = 1.0_num/clight**2
!$OMP PARALLEL DO PRIVATE(ip, usq, gaminv)
do ip=1,np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv = 1.0_num/sqrt(1.0_num + usq)
    xp(ip) = xp(ip) + uxp(ip)*gaminv*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv*dt
enddo
!$OMP END PARALLEL DO

return
end

!===============================================================================
!  Push the particle velocity with E field
subroutine epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
!===============================================================================

use constants
implicit none
integer :: np
real(num):: uxp(np),uyp(np),uzp(np)
real(num):: ex(np),ey(np),ez(np)
real(num):: q,m,dt

integer :: ip
real(num):: const

const = q*dt/m
!$OMP PARALLEL DO PRIVATE(ip)
do ip=1,np
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const
enddo
!$OMP END PARALLEL DO

return
end subroutine epush_v

!===============================================================================
!  Push the particle velocity with B field (Boris algorithm)
! --- fast b-field rotation algorithm
subroutine bpush_v(np,uxp,uyp,uzp,bx,by,bz,q,m,dt)
!===============================================================================

use constants
implicit none
integer :: np
real(num):: uxp(np), uyp(np), uzp(np)
real(num):: bx(np), by(np), bz(np)
real(num):: q,m,dt,gaminv
integer :: ip
real(num):: const,clghtisq,sx,sy,sz,tx,ty,tz,tsqi,uxppr,uyppr,uzppr,usq

const = q*dt*0.5_num/m
clghtisq = 1.0_num/clight**2

!$OMP PARALLEL DO PRIVATE(ip, tx, ty, tz, tsqi, sx, sy, sz, uxppr, uyppr, uzppr, usq, gaminv)
do ip=1,np
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
enddo
!$OMP END PARALLEL DO

return

end subroutine bpush_v

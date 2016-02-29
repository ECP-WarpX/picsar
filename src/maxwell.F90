!===============================================================================
! PUSH B field half a time step
!==============================================================================
SUBROUTINE push_bfield
USE constants
USE params
USE fields
USE shared_data
USE time_stat
IMPLICIT NONE

REAL(num) :: tmptime
tmptime = MPI_WTIME()

CALL pxrpush_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                       &
    0.5_num*dt/dx*xcoeffs,0.5_num*dt/dy*ycoeffs,0.5_num*dt/dz*zcoeffs,  &
    nx,ny,nz, norderx,nordery,norderz,                                  &
    nxguards,nyguards,nzguards,nxs,nys,nzs,                             &
    l_nodalgrid)

localtimes(5) = localtimes(5) + (MPI_WTIME() - tmptime)

END SUBROUTINE push_bfield


!===============================================================================
! PUSH E field a full  time step
!==============================================================================
SUBROUTINE push_efield
USE constants
USE params
USE fields
USE shared_data
USE time_stat
IMPLICIT NONE

REAL(num) :: tmptime
tmptime = MPI_WTIME()

CALL pxrpush_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,clight**2*mu0*dt,        &
    clight**2*dt/dx*xcoeffs,clight**2*dt/dy*ycoeffs,                           &
    clight**2*dt/dz*zcoeffs,nx,ny,nz,                                          &
    norderx,nordery,norderz,                                                   &
    nxguards,nyguards,nzguards,nxs,nys,nzs,                                    &
    l_nodalgrid)

localtimes(7) = localtimes(7) + (MPI_WTIME() - tmptime)

END SUBROUTINE push_efield


!===============================================================================
! PUSH ELECTRIC FIELD YEE 3D ARBITRARY ORDER
subroutine pxrpush_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,    &
                                 dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                 norderx,nordery,norderz,             &
                                 nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                 l_nodalgrid)
!===============================================================================
use constants
use omp_lib
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: Jx, Jy, Jz
real(num), intent(IN) :: mudt,dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
    ist = 0
else
    ist = 1
end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j,i)
!$OMP DO COLLAPSE(3)
! advance Ex
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            Ex(j,k,l) = Ex(j,k,l) - mudt  * Jx(j,k,l)
            do i = 1, nordery/2
                Ex(j,k,l) = Ex(j,k,l) + dtsdy(i) * (Bz(j,k+i-ist,l)   - Bz(j,k-i,l  ))
            end do
            do i = 1, norderz/2
                Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i-ist)   - By(j,k  ,l-i))
            end do
        end do
    end do
end do
!$OMP END DO
!$OMP DO COLLAPSE(3)
! advance Ey
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            Ey(j,k,l) = Ey(j,k,l) - mudt  * Jy(j,k,l)
            do i = 1, norderx/2
                Ey(j,k,l) = Ey(j,k,l) - dtsdx(i) * (Bz(j+i-ist,k,l)   - Bz(j-i,k,l))
            end do
            do i = 1, norderz/2
                Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist)   - Bx(j,k,l-i))
            end do
        end do
    end do
end do
!$OMP END DO
!$OMP DO COLLAPSE(3)
! advance Ez
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            Ez(j,k,l) = Ez(j,k,l) - mudt  * Jz(j,k,l)
            do i = 1, norderx/2
                Ez(j,k,l) = Ez(j,k,l) + dtsdx(i) * (By(j+i-ist,k,l) - By(j-i,k  ,l))
            end do
            do i = 1, nordery/2
                Ez(j,k,l) = Ez(j,k,l) - dtsdy(i) * (Bx(j,k+i-ist,l) - Bx(j  ,k-i,l))
            end do
        end do
    end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pxrpush_em3d_evec_norder

subroutine pxrpush_em2d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,    &
                                 dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                 norderx,nordery,norderz,             &
                                 nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                 l_nodalgrid)
!===============================================================================
use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: Jx, Jy, Jz
real(num), intent(IN) :: mudt,dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
    ist = 0
else
    ist = 1
end if

k = 0

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j,i)
!!$OMP DO COLLAPSE(3)
! advance Ex

do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            Ex(j,k,l) = Ex(j,k,l) - mudt  * Jx(j,k,l)
            do i = 1, norderz/2
                Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i-ist)   - By(j,k  ,l-i))
            end do
        end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance Ey
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            Ey(j,k,l) = Ey(j,k,l) - mudt  * Jy(j,k,l)
            do i = 1, norderx/2
                Ey(j,k,l) = Ey(j,k,l) - dtsdx(i) * (Bz(j+i-ist,k,l)   - Bz(j-i,k,l))
            end do
            do i = 1, norderz/2
                Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist)   - Bx(j,k,l-i))
            end do
    end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance Ez
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            Ez(j,k,l) = Ez(j,k,l) - mudt  * Jz(j,k,l)
            do i = 1, norderx/2
                Ez(j,k,l) = Ez(j,k,l) + dtsdx(i) * (By(j+i-ist,k,l) - By(j-i,k  ,l))
            end do
        end do
end do
!!$OMP END DO
!!$OMP END PARALLEL
return
end subroutine pxrpush_em2d_evec_norder



subroutine pxrpush_em2d_evec(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,    &
                                 dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                 nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                 l_nodalgrid)
!===============================================================================
use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: Jx, Jy, Jz
real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
    ist = 0
else
    ist = 1
end if

k = 0

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j,i)
!!$OMP DO COLLAPSE(3)
! advance Ex

do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            Ex(j,k,l) = Ex(j,k,l) - mudt  * Jx(j,k,l)
			Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l+1-ist)   - By(j,k  ,l-1))
        end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance Ey
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            Ey(j,k,l) = Ey(j,k,l) - mudt  * Jy(j,k,l)
			Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j+1-ist,k,l)   - Bz(j-1,k,l))
			Ey(j,k,l) = Ey(j,k,l) + dtsdz * (Bx(j,k,l+1-ist)   - Bx(j,k,l-1))
    end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance Ez
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            Ez(j,k,l) = Ez(j,k,l) - mudt  * Jz(j,k,l)
			Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j+1-ist,k,l) - By(j-1,k  ,l))
        end do
end do
!!$OMP END DO
!!$OMP END PARALLEL
return
end subroutine pxrpush_em2d_evec

!===============================================================================
! PUSH ELECTRIC FIELD YEE 3D ORDER 2
!===============================================================================
subroutine pxrpush_em3d_evec(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,    &
                                 dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                 nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                 l_nodalgrid)
use constants

integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: jx,jy,jz
real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz
integer(idp):: j,k,l
logical :: l_nodalgrid

! advance Ex
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l)   - Bz(j,k-1,l  )) &
            - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
            - mudt  * jx(j,k,l)
        end do
    end do
end do

  ! advance Ey
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
            + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
            - mudt  * jy(j,k,l)
        end do
    end do
end do

  ! advance Ez 
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
            - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
            - mudt  * jz(j,k,l)
        end do
    end do
end do

return
end subroutine pxrpush_em3d_evec


!===============================================================================
! PUSH MAGNETIC FIELD YEE 3D ARBITRARY ORDER
subroutine pxrpush_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                  &
                                dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                norderx,nordery,norderz,             &
                                nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                l_nodalgrid)
!===============================================================================
use constants
use omp_lib
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN) :: dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
ist = 0
else
ist = 1
end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j,i)
!$OMP DO COLLAPSE(3)
! advance Bx
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            do i = 1, nordery/2
                Bx(j,k,l) = Bx(j,k,l) - dtsdy(i) * (Ez(j,k+i,l  ) - Ez(j,k-i+ist,l))
            end do
            do i = 1, norderz/2
                Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,  l+i) - Ey(j,k,l-i+ist))
            end do
        end do
    end do
end do
!$OMP END DO
!$OMP DO COLLAPSE(3)
! advance By
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            do i = 1, norderx/2
                By(j,k,l) = By(j,k,l) + dtsdx(i) * (Ez(j+i,k,l  ) - Ez(j-i+ist,k,l))
            end do
            do i = 1, norderz/2
                By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j  ,k,l+i) - Ex(j,k,l-i+ist))
            end do
        end do
    end do
end do
!$OMP END DO
!$OMP DO COLLAPSE(3)
! advance Bz
do l = -nzs, nz+nzs
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs
            do i = 1, norderx/2
                Bz(j,k,l) = Bz(j,k,l) - dtsdx(i) * (Ey(j+i,k,l) - Ey(j-i+ist,k,l))
            end do
            do i = 1, nordery/2
                Bz(j,k,l) = Bz(j,k,l) + dtsdy(i) * (Ex(j,k+i,l) - Ex(j,k-i+ist,l))
            end do
        end do
    end do
end do
!$OMP END DO
!$OMP END PARALLEL
return

end subroutine pxrpush_em3d_bvec_norder

subroutine pxrpush_em2d_bvec_norder(ex,ey,ez,bx,by,bz,                  &
                                dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                norderx,nordery,norderz,             &
                                nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                l_nodalgrid)
!===============================================================================
use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN) :: dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
ist = 0
else
ist = 1
end if

k = 0

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j,i)
!!$OMP DO COLLAPSE(3)
! advance Bx
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            do i = 1, norderz/2
                Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,  l+i) - Ey(j,k,l-i+ist))
            end do
        end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance By
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            do i = 1, norderx/2
                By(j,k,l) = By(j,k,l) + dtsdx(i) * (Ez(j+i,k,l  ) - Ez(j-i+ist,k,l))
            end do
            do i = 1, norderz/2
                By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j  ,k,l+i) - Ex(j,k,l-i+ist))
            end do
        end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance Bz
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
            do i = 1, norderx/2
                Bz(j,k,l) = Bz(j,k,l) - dtsdx(i) * (Ey(j+i,k,l) - Ey(j-i+ist,k,l))
            end do
        end do
end do
!!$OMP END DO
!!$OMP END PARALLEL
return

end subroutine pxrpush_em2d_bvec_norder


subroutine pxrpush_em2d_bvec(ex,ey,ez,bx,by,bz,                  &
                                dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                                nxguard,nyguard,nzguard,nxs,nys,nzs, &
                                l_nodalgrid)
!===============================================================================
use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
ist = 0
else
ist = 1
end if

k = 0

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j,i)
!!$OMP DO COLLAPSE(3)
! advance Bx
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
                Bx(j,k,l) = Bx(j,k,l) + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l-1+ist))
        end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance By
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
                By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j-1+ist,k,l))
                By(j,k,l) = By(j,k,l) - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l-1+ist))
        end do
end do
!!$OMP END DO
!!$OMP DO COLLAPSE(3)
! advance Bz
do l = -nzs, nz+nzs
        do j = -nxs, nx+nxs
                Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j-1+ist,k,l))
        end do
end do
!!$OMP END DO
!!$OMP END PARALLEL
return

end subroutine pxrpush_em2d_bvec

!===============================================================================
! PUSH MAGNETIC FIELD YEE 3D ORDER 2
!===============================================================================
subroutine pxrpush_em3d_bvec(ex,ey,ez,bx,by,bz,                   &
                          dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
                          nxguard,nyguard,nzguard,nxs,nys,nzs, &
                          l_nodalgrid)
use constants
integer(idp):: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs
real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(num), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(idp) :: j,k,l
logical :: l_nodalgrid

! advance Bx
do l = -nzs, nz+nzs-1
    do k = -nys, ny+nys-1
        do j = -nxs, nx+nxs
            Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j,k+1,l  ) - Ez(j,k,l)) &
            + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
        end do
    end do
end do

! advance By
do l = -nzs, nz+nzs-1
    do k = -nys, ny+nys
        do j = -nxs, nx+nxs-1
            By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &
            - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l))
        end do
    end do
end do

! advance Bz
do l = -nzs, nz+nzs
    do k = -nys, ny+nys-1
        do j = -nxs, nx+nxs-1
            Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j,k,l)) &
            + dtsdy * (Ex(j,k+1,l) - Ex(j,k,l))
        end do
    end do
end do
end subroutine pxrpush_em3d_bvec



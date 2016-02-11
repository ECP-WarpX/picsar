!=================================================================================
! Gathering of electric field from Yee grid ("energy conserving") on particles
! at order 1
SUBROUTINE gete3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg)
!=================================================================================

USE omp_lib
USE constants
IMPLICIT NONE
INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,ex,ey,ez
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1) :: sx
REAL(num), DIMENSION(0:1) :: sy
REAL(num), DIMENSION(0:1) :: sz
REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
ALLOCATE(sx0(0:1),sy0(0:1),sz0(0:1))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np
    
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=floor(x)
    j0=floor(x-0.5_num)
    k=floor(y)
    k0=floor(y-0.5_num)
    l=floor(z)
    l0=floor(z-0.5_num)
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
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
    ex(ip) = ex(ip) + sx(0)*sy0(0)*sz0(0)*exg(j0,k,l)
    ex(ip) = ex(ip) + sx(1)*sy0(0)*sz0(0)*exg(j0+1,k,l)
    ex(ip) = ex(ip) + sx(0)*sy0(1)*sz0(0)*exg(j0,k+1,l)
    ex(ip) = ex(ip) + sx(1)*sy0(1)*sz0(0)*exg(j0+1,k+1,l)
    ex(ip) = ex(ip) + sx(0)*sy0(0)*sz0(1)*exg(j0,k,l+1)
    ex(ip) = ex(ip) + sx(1)*sy0(0)*sz0(1)*exg(j0+1,k,l+1)
    ex(ip) = ex(ip) + sx(0)*sy0(1)*sz0(1)*exg(j0,k+1,l+1)
    ex(ip) = ex(ip) + sx(1)*sy0(1)*sz0(1)*exg(j0+1,k+1,l+1)
    
    ! Compute Ey on particle
    ey(ip) = ey(ip) + sx0(0)*sy(0)*sz0(0)*eyg(j,k0,l)
    ey(ip) = ey(ip) + sx0(1)*sy(0)*sz0(0)*eyg(j+1,k0,l)
    ey(ip) = ey(ip) + sx0(0)*sy(1)*sz0(0)*eyg(j,k0+1,l)
    ey(ip) = ey(ip) + sx0(1)*sy(1)*sz0(0)*eyg(j+1,k0+1,l)
    ey(ip) = ey(ip) + sx0(0)*sy(0)*sz0(1)*eyg(j,k0,l+1)
    ey(ip) = ey(ip) + sx0(1)*sy(0)*sz0(1)*eyg(j+1,k0,l+1)
    ey(ip) = ey(ip) + sx0(0)*sy(1)*sz0(1)*eyg(j,k0+1,l+1)
    ey(ip) = ey(ip) + sx0(1)*sy(1)*sz0(1)*eyg(j+1,k0+1,l+1)
    
    ! Compute Ez on particle
    ez(ip) = ez(ip) + sx0(0)*sy0(0)*sz(0)*ezg(j,k,l0)
    ez(ip) = ez(ip) + sx0(1)*sy0(0)*sz(0)*ezg(j+1,k,l0)
    ez(ip) = ez(ip) + sx0(0)*sy0(1)*sz(0)*ezg(j,k+1,l0)
    ez(ip) = ez(ip) + sx0(1)*sy0(1)*sz(0)*ezg(j+1,k+1,l0)
    ez(ip) = ez(ip) + sx0(0)*sy0(0)*sz(1)*ezg(j,k,l0+1)
    ez(ip) = ez(ip) + sx0(1)*sy0(0)*sz(1)*ezg(j+1,k,l0+1)
    ez(ip) = ez(ip) + sx0(0)*sy0(1)*sz(1)*ezg(j,k+1,l0+1)
    ez(ip) = ez(ip) + sx0(1)*sy0(1)*sz(1)*ezg(j+1,k+1,l0+1)
END DO
!!$OMP END PARALLEL DO
DEALLOCATE(sx0,sz0)
RETURN
END SUBROUTINE gete3d_energy_conserving_1_1_1



!=================================================================================
! Gathering of Magnetic field from Yee grid ("energy conserving") on particles
! at order 1
SUBROUTINE getb3d_energy_conserving_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg)
!=================================================================================

USE omp_lib
USE constants
IMPLICIT NONE
INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,bx,by,bz
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1) :: sx
REAL(num), DIMENSION(0:1) :: sy
REAL(num), DIMENSION(0:1) :: sz
REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
ALLOCATE(sx0(0:1),sy0(0:1),sz0(0:1))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np
    
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=floor(x)
    j0=floor(x-0.5_num)
    k=floor(y)
    k0=floor(y-0.5_num)
    l=floor(z)
    l0=floor(z-0.5_num)
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
    sy0( 0) = 1.0_num-yint
    sy0( 1) = yint
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint
    
    ! Compute Bx on particle
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)
    
    ! Compute By on particle
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)
    
    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)
END DO
!!$OMP END PARALLEL DO
DEALLOCATE(sx0,sz0)
RETURN
END SUBROUTINE getb3d_energy_conserving_1_1_1


!=================================================================================
! Gathering of electric field from Yee grid ("energy conserving") on particles
! At arbitrary order. WARNING: Highly unoptimized routine
SUBROUTINE pxrgete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,       &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,exg,eyg,ezg,l_lower_order_in_v)
!=================================================================================
USE omp_lib
USE constants
IMPLICIT NONE
INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), dimension(np) :: xp,yp,zp,ex,ey,ez
LOGICAL :: l4symtry,l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
REAL(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
REAL(num), dimension(:), allocatable :: sx0,sy0,sz0
REAL(num), parameter :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num


dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
iymin = -int(noy/2)
iymax =  int((noy+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1

IF (l_lower_order_in_v) THEN
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
ELSE
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
END IF

ALLOCATE(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

signx = 1.0_num
signy = 1.0_num
!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint, &
!!$OMP   sx,sy,sz,sx0,sy0,sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq,ozint,zintsq,ozintsq)
DO ip=1,np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    IF (l_lower_order_in_v) THEN
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x-0.5_num)
        ELSE
            j=floor(x)
            j0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y-0.5_num)
        ELSE
            k=floor(y)
            k0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z-0.5_num)
        ELSE
            l=floor(z)
            l0=floor(z)
        END IF
    ELSE
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x)
        ELSE
            j=floor(x)
            j0=floor(x-0.5_num)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y)
        ELSE
            k=floor(y)
            k0=floor(y-0.5_num)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z)
        ELSE
            l=floor(z)
            l0=floor(z-0.5_num)
        END IF
    END IF

    xint=x-j
    yint=y-k
    zint=z-l

    IF (nox==1) THEN
        sx( 0) = 1.0_num-xint
        sx( 1) = xint
    ELSEIF (nox==2) THEN
        xintsq = xint*xint
        sx(-1) = 0.5_num*(0.5_num-xint)**2
        sx( 0) = 0.75_num-xintsq
        sx( 1) = 0.5_num*(0.5_num+xint)**2
    ELSEIF (nox==3) THEN
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx( 2) = onesixth*xintsq*xint
    END IF

    IF (noy==1) THEN
        sy( 0) = 1.0_num-yint
        sy( 1) = yint
    ELSEIF (noy==2) THEN
        yintsq = yint*yint
        sy(-1) = 0.5_num*(0.5_num-yint)**2
        sy( 0) = 0.75_num-yintsq
        sy( 1) = 0.5_num*(0.5_num+yint)**2
    ELSEIF (noy==3) THEN
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy( 2) = onesixth*yintsq*yint
    END IF

    IF (noz==1) THEN
        sz( 0) = 1.0_num-zint
        sz( 1) = zint
    ELSEIF (noz==2) THEN
        zintsq = zint*zint
        sz(-1) = 0.5_num*(0.5_num-zint)**2
        sz( 0) = 0.75_num-zintsq
        sz( 1) = 0.5_num*(0.5_num+zint)**2
    ELSEIF (noz==3) THEN
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz( 2) = onesixth*zintsq*zint
    END IF

    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0

    IF (l_lower_order_in_v) THEN

        IF (nox==1) THEN
            sx0( 0) = 1.0_num
        ELSEIF (nox==2) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==3) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num
        ELSEIF (noy==2) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==3) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        END IF

        IF (noz==1) THEN
            sz0( 0) = 1.0_num
        ELSEIF (noz==2) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==3) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        END IF

    ELSE

        IF (nox==1) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==2) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        ELSEIF (nox==3) THEN
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==2) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        ELSEIF (noy==3) THEN
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END IF

        IF (noz==1) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==2) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        ELSEIF (noz==3) THEN
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END IF

    END IF

    DO ll = izmin, izmax+1
        DO kk = iymin, iymax+1
            DO jj = ixmin0, ixmax0
                ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)*signx
            END DO
        END DO
    END DO

    DO ll = izmin, izmax+1
        DO kk = iymin0, iymax0
            DO jj = ixmin, ixmax+1
                ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)*signy
            END DO
        END DO
    END DO

    DO ll = izmin0, izmax0
        DO kk = iymin, iymax+1
            DO jj = ixmin, ixmax+1
                ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj,k+kk,l0+ll)
            END DO
        END DO
    END DO

END DO
!!$OMP END PARALLEL DO
DEALLOCATE(sx0,sy0,sz0)

RETURN
END SUBROUTINE pxrgete3d_n_energy_conserving

!=================================================================================
! Gathering of Magnetic field from Yee grid ("energy conserving") on particles
! At arbitrary order. WARNING: Highly unoptimized routine
SUBROUTINE pxrgetb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,    &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,bxg,byg,bzg,l_lower_order_in_v)
!=================================================================================
USE omp_lib
USE constants
IMPLICIT NONE
INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,bx,by,bz
LOGICAL :: l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
REAL(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
iymin = -int(noy/2)
iymax =  int((noy+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1


IF (l_lower_order_in_v) THEN
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
ELSE
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
END IF
ALLOCATE(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

signx = 1.0_num
signy = 1.0_num

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    IF (l_lower_order_in_v) THEN
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x-0.5_num)
        ELSE
            j=floor(x)
            j0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y-0.5_num)
        ELSE
            k=floor(y)
            k0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z-0.5_num)
        ELSE
            l=floor(z)
            l0=floor(z)
        END IF
    ELSE
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x)
        ELSE
            j=floor(x)
            j0=floor(x-0.5_num)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y)
        ELSE
            k=floor(y)
            k0=floor(y-0.5_num)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z)
        ELSE
            l=floor(z)
            l0=floor(z-0.5_num)
        END IF
    END IF

    xint=x-j
    yint=y-k
    zint=z-l

    IF (nox==1) THEN
        sx( 0) = 1.0_num-xint
        sx( 1) = xint
    ELSEIF (nox==2) THEN
        xintsq = xint*xint
        sx(-1) = 0.5_num*(0.5_num-xint)**2
        sx( 0) = 0.75_num-xintsq
        sx( 1) = 0.5_num*(0.5_num+xint)**2
    ELSEIF (nox==3) THEN
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx( 2) = onesixth*xintsq*xint
    END IF

    IF (noy==1) THEN
        sy( 0) = 1.0_num-yint
        sy( 1) = yint
    ELSEIF (noy==2) THEN
        yintsq = yint*yint
        sy(-1) = 0.5_num*(0.5_num-yint)**2
        sy( 0) = 0.75_num-yintsq
        sy( 1) = 0.5_num*(0.5_num+yint)**2
    ELSEIF (noy==3) THEN
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy( 2) = onesixth*yintsq*yint
    END IF

    IF (noz==1) THEN
        sz( 0) = 1.0_num-zint
        sz( 1) = zint
    ELSEIF (noz==2) THEN
        zintsq = zint*zint
        sz(-1) = 0.5_num*(0.5_num-zint)**2
        sz( 0) = 0.75_num-zintsq
        sz( 1) = 0.5_num*(0.5_num+zint)**2
    ELSEIF (noz==3) THEN
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz( 2) = onesixth*zintsq*zint
    END IF

    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0

    IF (l_lower_order_in_v) THEN
        IF (nox==1) THEN
            sx0( 0) = 1.0_num
        ELSEIF (nox==2) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==3) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num
        ELSEIF (noy==2) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==3) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        END IF

        IF (noz==1) THEN
            sz0( 0) = 1.0_num
        ELSEIF (noz==2) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==3) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        END IF
    ELSE

        IF (nox==1) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==2) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        ELSEIF (nox==3) THEN
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==2) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        ELSEIF (noy==3) THEN
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END IF
 
        IF (noz==1) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==2) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        ELSEIF (noz==3) THEN
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END IF
    END IF

    DO ll = izmin0, izmax0
        DO kk = iymin0, iymax0
            DO jj = ixmin, ixmax+1
                bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)*signx
            END DO
        END DO
    END DO

    DO ll = izmin0, izmax0
        DO kk = iymin, iymax+1
            DO jj = ixmin0, ixmax0
                by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)*signy
            END DO
        END DO
    END DO

    DO ll = izmin, izmax+1
        DO kk = iymin0, iymax0
            DO jj = ixmin0, ixmax0
                bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj,k0+kk,l+ll)
            END DO
        END DO
    END DO
END DO
!!OMP END PARALLEL DO
DEALLOCATE(sx0,sz0)

RETURN
END SUBROUTINE pxrgetb3d_n_energy_conserving

subroutine pxr_getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,bxg,byg,bzg,l4symtry,l_lower_order_in_v)
      use constants
      implicit none
      integer(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      real(num), dimension(np) :: xp,yp,zp,bx,by,bz
      logical(idp) :: l4symtry,l_lower_order_in_v
      real(num), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(num) :: xmin,ymin,zmin,dx,dy,dz
      integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
      real(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
                      xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
      real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
      real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(num), dimension(:), allocatable :: sx0,sy0,sz0
      real(num), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox+1)/2)-1
      iymin = -int(noy/2)
      iymax =  int((noy+1)/2)-1
      izmin = -int(noz/2)
      izmax =  int((noz+1)/2)-1


      if (l_lower_order_in_v) then
        ixmin0 = -int((nox-1)/2)
        ixmax0 =  int((nox)/2)
        iymin0 = -int((noy-1)/2)
        iymax0 =  int((noy)/2)
        izmin0 = -int((noz-1)/2)
        izmax0 =  int((noz)/2)
      else
        ixmin0 = -int((nox)/2)
        ixmax0 =  int((nox+1)/2)
        iymin0 = -int((noy)/2)
        iymax0 =  int((noy+1)/2)
        izmin0 = -int((noz)/2)
        izmax0 =  int((noz+1)/2)
      end if
      allocate(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

      signx = 1.
      signy = 1.

      sx=0
      sy=0.
      sz=0.
      sx0=0.
      sy0=0.
      sz0=0.

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
          if (y<0.) then
            y = -y
            signy = -1.
          else
            signy = 1.
          end if
        end if

        if (l_lower_order_in_v) then
         if (nox==2*(nox/2)) then
          j=nint(x)
          j0=floor(x-0.5)
         else
          j=floor(x)
          j0=floor(x)
         end if
         if (noy==2*(noy/2)) then
          k=nint(y)
          k0=floor(y-0.5)
         else
          k=floor(y)
          k0=floor(y)
         end if
         if (noz==2*(noz/2)) then
          l=nint(z)
          l0=floor(z-0.5)
         else
          l=floor(z)
          l0=floor(z)
         end if
        else
          if (nox==2*(nox/2)) then
            j=nint(x)
            j0=floor(x)
          else
            j=floor(x)
            j0=floor(x-0.5)
          end if
          if (noy==2*(noy/2)) then
            k=nint(y)
            k0=floor(y)
          else
            k=floor(y)
            k0=floor(y-0.5)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z)
          else
            l=floor(z)
            l0=floor(z-0.5)
          end if
        end if

        xint=x-j
        yint=y-k
        zint=z-l
        
        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy( 0) = 1.-yint
          sy( 1) = yint
        elseif (noy==2) then
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

        xint=x-0.5-j0
        yint=y-0.5-k0
        zint=z-0.5-l0

        if (l_lower_order_in_v) then

         if (nox==1) then
          sx0( 0) = 1.
         elseif (nox==2) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         end if

         if (noy==1) then
          sy0( 0) = 1.
         elseif (noy==2) then
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         elseif (noy==3) then
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         end if

         if (noz==1) then
          sz0( 0) = 1.
         elseif (noz==2) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         end if

        else

         if (nox==1) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==2) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
         end if

         if (noy==1) then
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         elseif (noy==2) then
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
         end if
         
         if (noz==1) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==2) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
         end if

        end if
        
        do ll = izmin0, izmax0
          do kk = iymin0, iymax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)*signx
            end do
          end do
        end do

        do ll = izmin0, izmax0
          do kk = iymin, iymax+1
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)*signy
            end do
          end do
        end do

        do ll = izmin, izmax+1
          do kk = iymin0, iymax0
            do jj = ixmin0, ixmax0
              bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj,k0+kk,l+ll)
            end do
          end do
        end do
                
     end do
     deallocate(sx0,sz0)

   return
 end subroutine pxr_getb3d_n_energy_conserving


  subroutine pxr_gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,exg,eyg,ezg,l4symtry,l_lower_order_in_v)
      use constants
      implicit none
      integer(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      real(num), dimension(np) :: xp,yp,zp,ex,ey,ez
      logical(idp) :: l4symtry,l_lower_order_in_v
      real(num), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(num) :: xmin,ymin,zmin,dx,dy,dz
      integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
      real(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
                      xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
      real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
      real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(num), dimension(:), allocatable :: sx0,sy0,sz0
      real(num), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox+1)/2)-1
      iymin = -int(noy/2)
      iymax =  int((noy+1)/2)-1
      izmin = -int(noz/2)
      izmax =  int((noz+1)/2)-1

      if (l_lower_order_in_v) then
        ixmin0 = -int((nox-1)/2)
        ixmax0 =  int((nox)/2)
        iymin0 = -int((noy-1)/2)
        iymax0 =  int((noy)/2)
        izmin0 = -int((noz-1)/2)
        izmax0 =  int((noz)/2)
      else
        ixmin0 = -int((nox)/2)
        ixmax0 =  int((nox+1)/2)
        iymin0 = -int((noy)/2)
        iymax0 =  int((noy+1)/2)
        izmin0 = -int((noz)/2)
        izmax0 =  int((noz+1)/2)
      end if
      allocate(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

      signx = 1.
      signy = 1.

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
          if (y<0.) then
            y = -y
            signy = -1.
          else
            signy = 1.
          end if
        end if
         
        if (l_lower_order_in_v) then
          if (nox==2*(nox/2)) then
            j=nint(x)
            j0=floor(x-0.5)
          else
            j=floor(x)
            j0=floor(x)
          end if
          if (noy==2*(noy/2)) then
            k=nint(y)
            k0=floor(y-0.5)
          else
            k=floor(y)
            k0=floor(y)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z-0.5)
          else
            l=floor(z)
            l0=floor(z)
          end if
        else
          if (nox==2*(nox/2)) then
            j=nint(x)
            j0=floor(x)
          else
            j=floor(x)
            j0=floor(x-0.5)
          end if
          if (noy==2*(noy/2)) then
            k=nint(y)
            k0=floor(y)
          else
            k=floor(y)
            k0=floor(y-0.5)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z)
          else
            l=floor(z)
            l0=floor(z-0.5)
          end if
        end if

        xint=x-j
        yint=y-k
        zint=z-l

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy( 0) = 1.-yint
          sy( 1) = yint
        elseif (noy==2) then
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

        xint=x-0.5-j0
        yint=y-0.5-k0
        zint=z-0.5-l0

        if (l_lower_order_in_v) then

         if (nox==1) then
          sx0( 0) = 1.
         elseif (nox==2) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         end if

         if (noy==1) then
          sy0( 0) = 1.
         elseif (noy==2) then
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         elseif (noy==3) then
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         end if

         if (noz==1) then
          sz0( 0) = 1.
         elseif (noz==2) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         end if

        else

         if (nox==1) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==2) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
         end if

         if (noy==1) then
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         elseif (noy==2) then
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
         end if

        
         if (noz==1) then
            sz0( 0) = 1.-zint
            sz0( 1) = zint
         elseif (noz==2) then
            zintsq = zint*zint
            sz0(-1) = 0.5*(0.5-zint)**2
            sz0( 0) = 0.75-zintsq
            sz0( 1) = 0.5*(0.5+zint)**2
         elseif (noz==3) then
            ozint = 1.-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.-zint/2)
            sz0( 1) = twothird-ozintsq*(1.-ozint/2)
            sz0( 2) = onesixth*zintsq*zint
         end if

        end if
        
        do ll = izmin, izmax+1
          do kk = iymin, iymax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)*signx
            end do
          end do
        end do

        do ll = izmin, izmax+1
          do kk = iymin0, iymax0
            do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)*signy
            end do
          end do
        end do

        do ll = izmin0, izmax0
          do kk = iymin, iymax+1
            do jj = ixmin, ixmax+1
              ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj,k+kk,l0+ll)
            end do
          end do
        end do
                     
     end do
     deallocate(sx0,sy0,sz0)

   return
 end subroutine pxr_gete3d_n_energy_conserving






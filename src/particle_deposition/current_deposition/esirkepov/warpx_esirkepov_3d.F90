! ______________________________________________________________________________
!> @brief
!> Esirkepov current deposition algorithm for linear, quadratic or cubic splines
!
!> @details
!> This subroutine can be used for several orders
!> WARNING: Highly unoptimized routine ---> USE INLINED ROUTINE
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
SUBROUTINE warpx_pxr_depose_jxjyjz_esirkepov_n( &
  jx, jx_nguard, jx_nvalid, &
  jy, jy_nguard, jy_nvalid, &
  jz, jz_nguard, jz_nvalid, &
  np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
  xmin,ymin,zmin,dt,dx,dy,dz,nox,noy,noz,l_particles_weight,l4symtry)
! ______________________________________________________________________________

  USE constants
  IMPLICIT NONE

  INTEGER, INTENT(IN)      :: jx_nguard(3),jy_nguard(3),jz_nguard(3), &
  	   		      jx_nvalid(3),jy_nvalid(3),jz_nvalid(3)
  INTEGER(idp)             :: np,nox,noy,noz
  REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
  REAL(num)                :: q,dt,dx,dy,dz,xmin,ymin,zmin
  REAL(num), INTENT(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1, &
  	     	       	         -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1, &
  	     	       	         -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1)
  REAL(num), INTENT(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1, &
  	     	       	         -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1, &
  	     	       	         -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1)
  REAL(num), INTENT(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1, &
  	     	       	         -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1, &
  	     	       	         -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1)
  REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
  REAL(num) :: xold,yold,zold,x,y,z,wq,wqx,wqy,wqz,vx,vy,vz,dts2dx,dts2dy,dts2dz
  REAL(num) :: invvol,invdtdx,invdtdy,invdtdz
  REAL(num) :: oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
  REAL(num) :: dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(:), ALLOCATABLE :: sx, sx0, dsx
  REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
  REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
  INTEGER(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,i,j,k,ic,jc,kc, &
  ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx, ndtody, ndtodz, &
  xl,xu,yl,yu,zl,zu
  LOGICAL(lp)  :: l_particles_weight,l4symtry

! PARAMETER INIT
  ndtodx = int(clight*dt/dx)
  ndtody = int(clight*dt/dy)
  ndtodz = int(clight*dt/dz)
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  dts2dx0 = 0.5_num*dtsdx0
  dts2dy0 = 0.5_num*dtsdy0
  dts2dz0 = 0.5_num*dtsdz0
  invvol = 1.0_num/(dx*dy*dz)
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)

  xl = -int(nox/2)-1-ndtodx
  xu = int((nox+1)/2)+1+ndtodx
  yl = -int(noy/2)-1-ndtody
  yu = int((noy+1)/2)+1+ndtody
  zl = -int(noz/2)-1-ndtodz
  zu = int((noz+1)/2)+1+ndtodz

ALLOCATE(sdx(xl:xu,yl:yu,zl:zu),sdy(xl:xu,yl:yu,zl:zu),sdz(xl:xu,yl:yu,zl:zu))
ALLOCATE(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
ALLOCATE(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
ALLOCATE(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num


DO ip=1,np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- applies 4-fold symmetry
    IF (l4symtry) THEN
        x=abs(x)
        y=abs(y)
        xold=abs(xold)
        yold=abs(yold)
        vx = (x-xold)/dtsdx0
        vy = (y-yold)/dtsdy0
    END IF
    !computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdy = dtsdy0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dy = dts2dy0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    y=yold
    z=zold
    DO icell = 1,ncells
        xold = x
        yold = y
        zold = z
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz

        ! --- computes particles "weights"
        IF (l_particles_weight) THEN
            wq=q*w(ip)
        ELSE
            wq=q*w(1)
        END IF
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

! --- finds node of cell containing particles for current positions
! --- (different for odd/even spline orders)
        IF (nox==2*(nox/2)) THEN
            iixp0=nint(x)
        ELSE
            iixp0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            ijxp0=nint(y)
        ELSE
            ijxp0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            ikxp0=nint(z)
        ELSE
            ikxp0=floor(z)
        END IF
! --- computes distance between particle and node for current positions
        xint=x-iixp0
        yint=y-ijxp0
        zint=z-ikxp0

! --- computes coefficients for node centered quantities
        SELECT CASE(nox)
        CASE(0)
            sx0( 0) = 1.0_num
        CASE(1)
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        CASE(2)
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        CASE(3)
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END SELECT

        SELECT CASE(noy)
        CASE(0)
            sy0( 0) = 1.0_num
        CASE(1)
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        CASE(2)
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        CASE(3)
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END SELECT

        SELECT CASE(noz)
        CASE(0)
            sz0( 0) = 1.0_num
        CASE(1)
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        CASE(2)
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        CASE(3)
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END SELECT

! --- finds node of cell containing particles for old positions
! --- (different for odd/even spline orders)
        IF (nox==2*(nox/2)) THEN
            iixp=nint(xold)
        ELSE
            iixp=floor(xold)
        END IF
        IF (noy==2*(noy/2)) THEN
            ijxp=nint(yold)
        ELSE
            ijxp=floor(yold)
        END IF
        IF (noz==2*(noz/2)) THEN
            ikxp=nint(zold)
        ELSE
            ikxp=floor(zold)
        END IF

! --- computes distance between particle and node for old positions
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.0_num;sy=0.0_num;sz=0.0_num

! --- computes coefficients for quantities centered between nodes
        SELECT CASE(nox)
        CASE(0)
            sx( 0+dix) = 1.0_num
        CASE(1)
            sx( 0+dix) = 1.0_num-xint
            sx( 1+dix) = xint
        CASE(2)
            xintsq = xint*xint
            sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
            sx( 0+dix) = 0.75_num-xintsq
            sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
        CASE(3)
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx(-1+dix) = onesixth*oxintsq*oxint
            sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx( 2+dix) = onesixth*xintsq*xint
        END SELECT

        SELECT CASE(noy)
        CASE(0)
            sy( 0+diy) = 1.0_num
        CASE(1)
            sy( 0+diy) = 1.0_num-yint
            sy( 1+diy) = yint
        CASE(2)
            yintsq = yint*yint
            sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
            sy( 0+diy) = 0.75_num-yintsq
            sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
        CASE(3)
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy(-1+diy) = onesixth*oyintsq*oyint
            sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy( 2+diy) = onesixth*yintsq*yint
        END SELECT

        SELECT CASE(noz)
        CASE(0)
            sz( 0+diz) = 1.0_num
        CASE(1)
            sz( 0+diz) = 1.0_num-zint
            sz( 1+diz) = zint
        CASE(2)
            zintsq = zint*zint
            sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
            sz( 0+diz) = 0.75_num-zintsq
            sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
        CASE(3)
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz(-1+diz) = onesixth*ozintsq*ozint
            sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz( 2+diz) = onesixth*zintsq*zint
        END SELECT

! --- computes coefficients difference
        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0

! --- computes min/max positions of current contributions
        ixmin = min(0_idp,dix)-int(nox/2)
        ixmax = max(0_idp,dix)+int((nox+1)/2)
        iymin = min(0_idp,diy)-int(noy/2)
        iymax = max(0_idp,diy)+int((noy+1)/2)
        izmin = min(0_idp,diz)-int(noz/2)
        izmax = max(0_idp,diz)+int((noz+1)/2)
! --- add current contributions
        DO k=izmin, izmax
            DO j=iymin, iymax
                DO i=ixmin, ixmax
                    ic = iixp0+i
                    jc = ijxp0+j
                    kc = ikxp0+k
                    IF(i<ixmax) THEN
                        sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
                        (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
                        IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                        jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)
                    END IF
                    IF(j<iymax) THEN
                        sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
                        (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
                        IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                        jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)
                    END IF
                    IF(k<izmax) THEN
                        sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
                        (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
                        IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                        jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)
                    END IF
                END DO
            END DO
        END DO

    END DO
END DO

DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)

RETURN
END SUBROUTINE warpx_pxr_depose_jxjyjz_esirkepov_n

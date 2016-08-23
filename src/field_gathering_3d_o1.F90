! ________________________________________________________________________________________
! 
! FIELD_GATHERING_3D_O1.F90
! 
! Field gathering subroutines in 3D at order 1
!
! List of subroutines:
!
! 
! - gete3d_energy_conserving_scalar_1_1_1
! - getb3d_energy_conserving_scalar_1_1_1
! - gete3d_energy_conserving_linear_1_1_1
! - getb3d_energy_conserving_linear_1_1_1
! - gete3d_energy_conserving_1_1_1
! - getb3d_energy_conserving_1_1_1
! 
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> Scalar version: gathering of electric field from Yee grid ("energy conserving") on particles
!> at order 1
!> @brief
SUBROUTINE gete3d_energy_conserving_scalar_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
!
! ________________________________________________________________________________________

USE omp_lib
USE constants
USE params

IMPLICIT NONE
INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez
logical                              :: l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
INTEGER(isp)                         :: ip, j, k, l
INTEGER(isp)                         :: ixmin, ixmax, iymin, iymax, izmin, izmax
INTEGER(isp)                         :: ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0
INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
REAL(num)                            :: dxi, dyi, dzi, x, y, z
REAL(num)                            :: xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1)            :: sx, sx0
REAL(num), DIMENSION(0:1)            :: sy, sy0
REAL(num), DIMENSION(0:1)            :: sz, sz0
REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = 0
ixmax = 0
iymin = 0
iymax = 0
izmin = 0
izmax = 0

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num

IF (l_lower_order_in_v) THEN

  ixmin0 = 0
  ixmax0 = 0
  iymin0 = 0
  iymax0 = 0
  izmin0 = 0
  izmax0 = 0

  DO ip=1,np
  
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
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
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
    
    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0
    
    sx0( 0) = 1.0_num
    sy0( 0) = 1.0_num
    sz0( 0) = 1.0_num
    
    do ll = izmin, izmax+1
      do kk = iymin, iymax+1
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)
        end do
      end do
    end do

    do ll = izmin, izmax+1
      do kk = iymin0, iymax0
        do jj = ixmin, ixmax+1
          ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)
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

    ! Debugging
!     IF (it.gt.0) THEN
!       print*,'ex,ey,ez',ex(ip),ey(ip),ez(ip)
!       print*,'j',j,k,l
!       print*,'j0',j0,k0,l0  
!       print*,'sx',sx(:)
!       print*,'sy',sy(:)   
!       print*,'sz',sz(:)   
!       print*,'sx0',sx0(:)
!       print*,'sy0',sy0(:)   
!       print*,'sz0',sz0(:)            
!       read*
!     ENDIF

  END DO

! __ l_lower_order_in_v false  _____________________________
ELSE

  ixmin0 = 0
  ixmax0 = 1
  iymin0 = 0
  iymax0 = 1
  izmin0 = 0
  izmax0 = 1

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
    
    do ll = izmin, izmax+1
      do kk = iymin, iymax+1
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)
        end do
      end do
    end do

    do ll = izmin, izmax+1
      do kk = iymin0, iymax0
        do jj = ixmin, ixmax+1
          ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)
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
      
  END DO
ENDIF

RETURN
END SUBROUTINE gete3d_energy_conserving_scalar_1_1_1



! ________________________________________________________________________________________
!> Scalar version: Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 1
!> @brief
!
!> This function is vectorized 
!> @details
SUBROUTINE getb3d_energy_conserving_scalar_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
! ________________________________________________________________________________________

USE omp_lib
USE constants
IMPLICIT NONE

INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np)             :: xp,yp,zp,bx,by,bz
logical                              :: l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
INTEGER(idp)                         :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num)                            :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1)            :: sx, sx0
REAL(num), DIMENSION(0:1)            :: sy, sy0
REAL(num), DIMENSION(0:1)            :: sz, sz0
REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num

ixmin = 0
ixmax = 0
iymin = 0
iymax = 0
izmin = 0
izmax = 0


IF (l_lower_order_in_v) THEN

  ixmin0 = 0
  ixmax0 = 0
  iymin0 = 0
  iymax0 = 0
  izmin0 = 0
  izmax0 = 0

  DO ip=1,np
  
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
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
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
    
    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0
    
    sx0( 0) = 1.0_num
    sy0( 0) = 1.0_num
    sz0( 0) = 1.0_num
    
    do ll = izmin0, izmax0
      do kk = iymin0, iymax0
        do jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)
        end do
      end do
    end do

    do ll = izmin0, izmax0
      do kk = iymin, iymax+1
        do jj = ixmin0, ixmax0
          by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)
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
  
  END DO

ELSE

  ixmin0 = 0
  ixmax0 = 1
  iymin0 = 0
  iymax0 = 1
  izmin0 = 0
  izmax0 = 1

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
    
      ! Compute shape factors
      xint=x-j
      yint=y-k
      zint=z-l    
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
    
      do ll = izmin0, izmax0
        do kk = iymin0, iymax0
          do jj = ixmin, ixmax+1
            bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)
          end do
        end do
      end do

      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          do jj = ixmin0, ixmax0
            by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)
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
  END DO
ENDIF
RETURN
END SUBROUTINE getb3d_energy_conserving_scalar_1_1_1

!=================================================================================
!> Gathering of electric field from Yee grid ("energy conserving") on particles
!> at order 1.
!> @brief
SUBROUTINE gete3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
!=================================================================================

USE omp_lib
USE constants
USE params
IMPLICIT NONE

INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez
LOGICAL                              :: l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
INTEGER(isp)                         :: ip, j, k, l
INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
REAL(num)                            :: dxi, dyi, dzi, x, y, z
REAL(num)                            :: xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1)            :: sx, sx0
REAL(num), DIMENSION(0:1)            :: sy, sy0
REAL(num), DIMENSION(0:1)            :: sz, sz0
REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num

IF (l_lower_order_in_v) THEN

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,ex,ey,ez)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR IVDEP
!DIR DISTRIBUTE POINT
#endif
  DO ip=1,np
  
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
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
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
    
    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0
    
    sx0( 0) = 1.0_num
    sy0( 0) = 1.0_num
    sz0( 0) = 1.0_num
    
    ! Compute Ex on particle
    ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
    ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
    ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
    ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
    
    ! Compute Ey on particle
    ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
    ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
    ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
    ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
    
    ! Compute Ez on particle
    ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
    ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
    ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
    ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)

    ! Debugging
!     IF (it.gt.0) THEN
!       print*,'ex,ey,ez',ex(ip),ey(ip),ez(ip)
!       print*,'j',j,k,l
!       print*,'j0',j0,k0,l0  
!       print*,'sx',sx(:)
!       print*,'sy',sy(:)   
!       print*,'sz',sz(:)   
!       print*,'sx0',sx0(:)
!       print*,'sy0',sy0(:)   
!       print*,'sz0',sz0(:)            
!       read*
!     ENDIF

  END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif

! l_lower_order_in_v false    
ELSE

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,ex,ey,ez)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR IVDEP
!DIR DISTRIBUTE POINT
#endif
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
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(0)*exg(j0+1,k,l)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(0)*exg(j0+1,k+1,l)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(1)*exg(j0+1,k,l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(1)*exg(j0+1,k+1,l+1)
    
      ! Compute Ey on particle
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(0)*eyg(j,k0+1,l)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(0)*eyg(j+1,k0+1,l)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(1)*eyg(j,k0+1,l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(1)*eyg(j+1,k0+1,l+1)
    
      ! Compute Ez on particle
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(1)*ezg(j,k,l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(1)*ezg(j+1,k,l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(1)*ezg(j,k+1,l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(1)*ezg(j+1,k+1,l0+1)
      
  END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif
ENDIF

RETURN
END SUBROUTINE gete3d_energy_conserving_1_1_1


!=================================================================================
!> Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 1.
!> @brief
!
! This function is vectorized                                      
!> @details
SUBROUTINE getb3d_energy_conserving_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
!=================================================================================

USE omp_lib
USE constants
IMPLICIT NONE

INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np)             :: xp,yp,zp,bx,by,bz
LOGICAL                              :: l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
INTEGER(isp)                         :: ip, j, k, l
INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
REAL(num)                            :: dxi, dyi, dzi
REAL(num)                            :: x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1)            :: sx, sx0
REAL(num), DIMENSION(0:1)            :: sy, sy0
REAL(num), DIMENSION(0:1)            :: sz, sz0
REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num

IF (l_lower_order_in_v) THEN

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,bx,by,bz)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR IVDEP
!DIR DISTRIBUTE POINT
#endif
  DO ip=1,np
  
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
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
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
    
    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0
    
    sx0( 0) = 1.0_num
    sy0( 0) = 1.0_num
    sz0( 0) = 1.0_num
    
    ! Compute Bx on particle
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
    
    ! Compute By on particle
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
    
    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1) 
  
  END DO
#if defined _OPENMP && _OPENMP>=201307
        !$OMP END SIMD 
#endif

ELSE

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,bx,by,bz)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR IVDEP
!DIR DISTRIBUTE POINT
#endif
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
    
      ! Compute shape factors
      xint=x-j
      yint=y-k
      zint=z-l    
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
#if defined _OPENMP && _OPENMP>=201307
        !$OMP END SIMD 
#endif
ENDIF
RETURN
END SUBROUTINE getb3d_energy_conserving_1_1_1

! ________________________________________________________________________________________
!> Field gathering CIC (order 1) with gathering of E and B merged in a single loop
!> @brief
!
!> This function is vectorized
!@details
SUBROUTINE geteb3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
! Input parameters:
! - np: number of particles
! - xp,yp,zp: particle position
! - ex,ey,ez: particle electric field
! - bx,by,bz: particle magnetic field
! - xmin,ymin,zmin: tile minimum grid position
! - dx,dy,dz: space step
! - dt: time step
! - nx,ny,nz: number of grid points in each direction
! - nxguard, nyguard, nzguard: number of guard cells in each direction 
! - exg,eyg,ezg: electric field grid
! - bxg,byg,bzg: magnetic field grid
! - lvect: vector size for cache blocking
! - l_lower_order_in_v: 
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  USE params
  
  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez,bx,by,bz
  LOGICAL                              :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: nn,n
  INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint, &
                xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  REAL(num), DIMENSION(0:1)            :: sx, sx0
  REAL(num), DIMENSION(0:1)            :: sy, sy0
  REAL(num), DIMENSION(0:1)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

IF (l_lower_order_in_v) THEN

  ! ___ Loop on partciles _______________________
  DO ip=1,np,lvect
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64      
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR IVDEP
!DIR DISTRIBUTE POINT
#endif
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1
  
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
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sy( 0) = 1.0_num-yint
      sy( 1) = yint
      sz( 0) = 1.0_num-zint
      sz( 1) = zint
    
      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0
    
      sx0( 0) = 1.0_num
      sy0( 0) = 1.0_num
      sz0( 0) = 1.0_num
    
      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
      ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
      ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
      ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
    
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
      ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
      ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
      ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
    
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
      ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
      ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
      ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
    
      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
      bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
    
      ! Compute By on particle
      by(nn) = by(nn) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
      by(nn) = by(nn) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
      bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1) 

    END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif
  ENDDO

ELSE

  ! ___ Loop on partciles _______________________
  DO ip=1,np,lvect
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64      
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* ALIGN(64,bx,by,bz)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR IVDEP
!DIR DISTRIBUTE POINT
#endif
    DO nn=ip,ip+MIN(lvect,np-ip+1)-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
    
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
      ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
      ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(0)*exg(j0+1,k,l)
      ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
      ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(0)*exg(j0+1,k+1,l)
      ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
      ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(1)*exg(j0+1,k,l+1)
      ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
      ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(1)*exg(j0+1,k+1,l+1)
    
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
      ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
      ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(0)*eyg(j,k0+1,l)
      ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(0)*eyg(j+1,k0+1,l)
      ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
      ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
      ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(1)*eyg(j,k0+1,l+1)
      ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(1)*eyg(j+1,k0+1,l+1)
    
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
      ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
      ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
      ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
      ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(1)*ezg(j,k,l0+1)
      ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(1)*ezg(j+1,k,l0+1)
      ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(1)*ezg(j,k+1,l0+1)
      ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(1)*ezg(j+1,k+1,l0+1)

      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
      bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
      bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
      bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
      bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
      bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
      bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
      bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)
    
      ! Compute By on particle
      by(nn) = by(nn) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
      by(nn) = by(nn) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
      by(nn) = by(nn) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
      by(nn) = by(nn) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
      by(nn) = by(nn) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
      by(nn) = by(nn) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
      by(nn) = by(nn) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
      by(nn) = by(nn) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
      bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
      bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
      bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
      bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
      bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
      bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
      bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)

    END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif
  ENDDO

ENDIF

RETURN
END SUBROUTINE geteb3d_energy_conserving_1_1_1
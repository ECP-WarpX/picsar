! ________________________________________________________________________________________
! 
! FIELD_GATHERING_3D_O1.F90
! 
! Field gathering subroutines in 3D at order 1
!
! List of subroutines:
!
! 
! - gete3d_energy_conserving_scalar_2_2_2
! - getb3d_energy_conserving_scalar_2_2_2
! - gete3d_energy_conserving_2_2_2
! - getb3d_energy_conserving_2_2_2
! - geteb3d_energy_conserving_2_2_2
! 
! ________________________________________________________________________________________

!=================================================================================
SUBROUTINE gete3d_energy_conserving_2_2_2(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,       &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
! Gathering of electric field from Yee grid ("energy conserving") on particles
! at order 2
! This function is vectorized                                      
!=================================================================================                                      
  USE omp_lib
  USE constants
  IMPLICIT NONE
  
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez
  LOGICAL                              :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
  REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  REAL(num), DIMENSION(-1:1) :: sx,sx0
  REAL(num), DIMENSION(-1:1) :: sy,sy0
  REAL(num), DIMENSION(-1:1) :: sz,sz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

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

#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
  DO ip=1,np
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    
    j=nint(x)
    j0=floor(x-0.5_num)
    k=nint(y)
    k0=floor(y-0.5_num)
    l=nint(z)
    l0=floor(z-0.5_num)
    
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
    xintsq = xint*xint
    sx(-1) = 0.5_num*(0.5_num-xint)**2
    sx( 0) = 0.75_num-xintsq
    sx( 1) = 0.5_num*(0.5_num+xint)**2
    
    yintsq = yint*yint
    sy(-1) = 0.5_num*(0.5_num-yint)**2
    sy( 0) = 0.75_num-yintsq
    sy( 1) = 0.5_num*(0.5_num+yint)**2
    
    zintsq = zint*zint
    sz(-1) = 0.5_num*(0.5_num-zint)**2
    sz( 0) = 0.75_num-zintsq
    sz( 1) = 0.5_num*(0.5_num+zint)**2
    
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
    ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(-1)*exg(j0,k-1,l-1)
    ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(-1)*exg(j0+1,k-1,l-1)
    ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(-1)*exg(j0,k,l-1)
    ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(-1)*exg(j0+1,k,l-1)
    ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(-1)*exg(j0,k+1,l-1)
    ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(-1)*exg(j0+1,k+1,l-1)
    ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(0)*exg(j0,k-1,l)
    ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(0)*exg(j0+1,k-1,l)
    ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
    ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(0)*exg(j0+1,k,l)
    ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
    ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(0)*exg(j0+1,k+1,l)
    ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(1)*exg(j0,k-1,l+1)
    ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(1)*exg(j0+1,k-1,l+1)
    ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
    ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(1)*exg(j0+1,k,l+1)
    ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
    ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(1)*exg(j0+1,k+1,l+1)
    
    ! Compute Ey on particle
    ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(-1)*eyg(j-1,k0,l-1)
    ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(-1)*eyg(j,k0,l-1)
    ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(-1)*eyg(j+1,k0,l-1)
    ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(-1)*eyg(j-1,k0+1,l-1)
    ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(-1)*eyg(j,k0+1,l-1)
    ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(-1)*eyg(j+1,k0+1,l-1)
    ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(0)*eyg(j-1,k0,l)
    ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
    ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
    ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(0)*eyg(j-1,k0+1,l)
    ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(0)*eyg(j,k0+1,l)
    ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(0)*eyg(j+1,k0+1,l)
    ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(1)*eyg(j-1,k0,l+1)
    ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
    ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
    ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(1)*eyg(j-1,k0+1,l+1)
    ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(1)*eyg(j,k0+1,l+1)
    ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(1)*eyg(j+1,k0+1,l+1)
    
    ! Compute Ez on particle
    ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(0)*ezg(j-1,k-1,l0)
    ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(0)*ezg(j,k-1,l0)
    ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(0)*ezg(j+1,k-1,l0)
    ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(0)*ezg(j-1,k,l0)
    ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
    ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
    ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(0)*ezg(j-1,k+1,l0)
    ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
    ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
    ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(1)*ezg(j-1,k-1,l0+1)
    ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(1)*ezg(j,k-1,l0+1)
    ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(1)*ezg(j+1,k-1,l0+1)
    ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(1)*ezg(j-1,k,l0+1)
    ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(1)*ezg(j,k,l0+1)
    ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(1)*ezg(j+1,k,l0+1)
    ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(1)*ezg(j-1,k+1,l0+1)
    ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(1)*ezg(j,k+1,l0+1)
    ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(1)*ezg(j+1,k+1,l0+1)

  END DO
#if defined _OPENMP && _OPENMP>=201307
			!$OMP END SIMD 
#endif

ELSE

#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
DO ip=1,np
  
		x = (xp(ip)-xmin)*dxi
		y = (yp(ip)-ymin)*dyi
		z = (zp(ip)-zmin)*dzi
	
		! Compute index of particle
		j=nint(x)
		j0=floor(x)
		k=nint(y)
		k0=floor(y)
		l=nint(z)
		l0=floor(z)
	
		xint=x-j
		yint=y-k
		zint=z-l
  
		! Compute shape factors
		xintsq = xint*xint
		sx(-1) = 0.5_num*(0.5_num-xint)**2
		sx( 0) = 0.75_num-xintsq
		sx( 1) = 0.5_num*(0.5_num+xint)**2
	
		yintsq = yint*yint
		sy(-1) = 0.5_num*(0.5_num-yint)**2
		sy( 0) = 0.75_num-yintsq
		sy( 1) = 0.5_num*(0.5_num+yint)**2
	
		zintsq = zint*zint
		sz(-1) = 0.5_num*(0.5_num-zint)**2
		sz( 0) = 0.75_num-zintsq
		sz( 1) = 0.5_num*(0.5_num+zint)**2
	
		xint=x-0.5_num-j0
		yint=y-0.5_num-k0
		zint=z-0.5_num-l0
  
		xintsq = xint*xint
		sx0(-1) = 0.5_num*(0.5_num-xint)**2
		sx0( 0) = 0.75_num-xintsq
		sx0( 1) = 0.5_num*(0.5_num+xint)**2
	
		yintsq = yint*yint
		sy0(-1) = 0.5_num*(0.5_num-yint)**2
		sy0( 0) = 0.75_num-yintsq
		sy0( 1) = 0.5_num*(0.5_num+yint)**2
	
		zintsq = zint*zint
		sz0(-1) = 0.5_num*(0.5_num-zint)**2
		sz0( 0) = 0.75_num-zintsq
		sz0( 1) = 0.5_num*(0.5_num+zint)**2
  
		! Compute Ex on particle
		ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(-1)*exg(j0-1,k-1,l-1)
		ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(-1)*exg(j0,k-1,l-1)
		ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(-1)*exg(j0+1,k-1,l-1)
		ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(-1)*exg(j0-1,k,l-1)
		ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(-1)*exg(j0,k,l-1)
		ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(-1)*exg(j0+1,k,l-1)
		ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(-1)*exg(j0-1,k+1,l-1)
		ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(-1)*exg(j0,k+1,l-1)
		ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(-1)*exg(j0+1,k+1,l-1)
		ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(0)*exg(j0-1,k-1,l)
		ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(0)*exg(j0,k-1,l)
		ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(0)*exg(j0+1,k-1,l)
		ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(0)*exg(j0-1,k,l)
		ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
		ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(0)*exg(j0+1,k,l)
		ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(0)*exg(j0-1,k+1,l)
		ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
		ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(0)*exg(j0+1,k+1,l)
		ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(1)*exg(j0-1,k-1,l+1)
		ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(1)*exg(j0,k-1,l+1)
		ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(1)*exg(j0+1,k-1,l+1)
		ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(1)*exg(j0-1,k,l+1)
		ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
		ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(1)*exg(j0+1,k,l+1)
		ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(1)*exg(j0-1,k+1,l+1)
		ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
		ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(1)*exg(j0+1,k+1,l+1)
  
		! Compute Ey on particle
		ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(-1)*eyg(j-1,k0-1,l-1)
		ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(-1)*eyg(j,k0-1,l-1)
		ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(-1)*eyg(j+1,k0-1,l-1)
		ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(-1)*eyg(j-1,k0,l-1)
		ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(-1)*eyg(j,k0,l-1)
		ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(-1)*eyg(j+1,k0,l-1)
		ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(-1)*eyg(j-1,k0+1,l-1)
		ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(-1)*eyg(j,k0+1,l-1)
		ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(-1)*eyg(j+1,k0+1,l-1)
		ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(0)*eyg(j-1,k0-1,l)
		ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(0)*eyg(j,k0-1,l)
		ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(0)*eyg(j+1,k0-1,l)
		ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(0)*eyg(j-1,k0,l)
		ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
		ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
		ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(0)*eyg(j-1,k0+1,l)
		ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(0)*eyg(j,k0+1,l)
		ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(0)*eyg(j+1,k0+1,l)
		ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(1)*eyg(j-1,k0-1,l+1)
		ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(1)*eyg(j,k0-1,l+1)
		ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(1)*eyg(j+1,k0-1,l+1)
		ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(1)*eyg(j-1,k0,l+1)
		ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
		ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
		ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(1)*eyg(j-1,k0+1,l+1)
		ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(1)*eyg(j,k0+1,l+1)
		ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(1)*eyg(j+1,k0+1,l+1)
  
		! Compute Ez on particle
		ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(-1)*ezg(j-1,k-1,l0-1)
		ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(-1)*ezg(j,k-1,l0-1)
		ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(-1)*ezg(j+1,k-1,l0-1)
		ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(-1)*ezg(j-1,k,l0-1)
		ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(-1)*ezg(j,k,l0-1)
		ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(-1)*ezg(j+1,k,l0-1)
		ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(-1)*ezg(j-1,k+1,l0-1)
		ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(-1)*ezg(j,k+1,l0-1)
		ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(-1)*ezg(j+1,k+1,l0-1)
		ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(0)*ezg(j-1,k-1,l0)
		ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(0)*ezg(j,k-1,l0)
		ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(0)*ezg(j+1,k-1,l0)
		ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(0)*ezg(j-1,k,l0)
		ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
		ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
		ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(0)*ezg(j-1,k+1,l0)
		ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
		ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
		ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(1)*ezg(j-1,k-1,l0+1)
		ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(1)*ezg(j,k-1,l0+1)
		ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(1)*ezg(j+1,k-1,l0+1)
		ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(1)*ezg(j-1,k,l0+1)
		ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(1)*ezg(j,k,l0+1)
		ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(1)*ezg(j+1,k,l0+1)
		ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(1)*ezg(j-1,k+1,l0+1)
		ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(1)*ezg(j,k+1,l0+1)
		ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(1)*ezg(j+1,k+1,l0+1)
  
  END DO
#if defined _OPENMP && _OPENMP>=201307
			!$OMP END SIMD 
#endif
ENDIF


RETURN
END SUBROUTINE gete3d_energy_conserving_2_2_2





!=================================================================================
!> Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 2
!> @brief
!  
!> This function is vectorized
!> @details
SUBROUTINE getb3d_energy_conserving_2_2_2(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,       &
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
  REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  REAL(num), DIMENSION(-1:1)           :: sx,sx0
  REAL(num), DIMENSION(-1:1)           :: sy,sy0
  REAL(num), DIMENSION(-1:1)           :: sz,sz0
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

#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif  
  DO ip=1,np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=nint(x)
    j0=floor(x-0.5_num)
    k=nint(y)
    k0=floor(y-0.5_num)
    l=nint(z)
    l0=floor(z-0.5_num)
    
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
    xintsq = xint*xint
    sx(-1) = 0.5_num*(0.5_num-xint)**2
    sx( 0) = 0.75_num-xintsq
    sx( 1) = 0.5_num*(0.5_num+xint)**2
    
    yintsq = yint*yint
    sy(-1) = 0.5_num*(0.5_num-yint)**2
    sy( 0) = 0.75_num-yintsq
    sy( 1) = 0.5_num*(0.5_num+yint)**2
    
    zintsq = zint*zint
    sz(-1) = 0.5_num*(0.5_num-zint)**2
    sz( 0) = 0.75_num-zintsq
    sz( 1) = 0.5_num*(0.5_num+zint)**2
    
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
    bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(0)*bxg(j-1,k0,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
    bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(0)*bxg(j-1,k0+1,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
    bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(1)*bxg(j-1,k0,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
    bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(1)*bxg(j-1,k0+1,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)
    
    ! Compute By on particle
    by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(0)*byg(j0,k-1,l0)
    by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(0)*byg(j0+1,k-1,l0)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
    by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(1)*byg(j0,k-1,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(1)*byg(j0+1,k-1,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)
    
    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(-1)*bzg(j0,k0,l-1)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(-1)*bzg(j0+1,k0,l-1)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(-1)*bzg(j0,k0+1,l-1)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(-1)*bzg(j0+1,k0+1,l-1)
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
  
  ELSE

#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif  
  DO ip=1,np
    
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=nint(x)
    j0=floor(x)
    k=nint(y)
    k0=floor(y)
    l=nint(z)
    l0=floor(z)
    
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
    xintsq = xint*xint
    sx(-1) = 0.5_num*(0.5_num-xint)**2
    sx( 0) = 0.75_num-xintsq
    sx( 1) = 0.5_num*(0.5_num+xint)**2
    
    yintsq = yint*yint
    sy(-1) = 0.5_num*(0.5_num-yint)**2
    sy( 0) = 0.75_num-yintsq
    sy( 1) = 0.5_num*(0.5_num+yint)**2
    
    zintsq = zint*zint
    sz(-1) = 0.5_num*(0.5_num-zint)**2
    sz( 0) = 0.75_num-zintsq
    sz( 1) = 0.5_num*(0.5_num+zint)**2
    
    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0
    
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2
    
    yintsq = yint*yint
    sy0(-1) = 0.5_num*(0.5_num-yint)**2
    sy0( 0) = 0.75_num-yintsq
    sy0( 1) = 0.5_num*(0.5_num+yint)**2
    
    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = 0.75_num-zintsq
    sz0( 1) = 0.5_num*(0.5_num+zint)**2
  
    ! Compute Bx on particle
    bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(-1)*bxg(j-1,k0-1,l0-1)
    bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(-1)*bxg(j,k0-1,l0-1)
    bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(-1)*bxg(j+1,k0-1,l0-1)
    bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(-1)*bxg(j-1,k0,l0-1)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(-1)*bxg(j,k0,l0-1)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(-1)*bxg(j+1,k0,l0-1)
    bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(-1)*bxg(j-1,k0+1,l0-1)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(-1)*bxg(j,k0+1,l0-1)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(-1)*bxg(j+1,k0+1,l0-1)
    bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(0)*bxg(j-1,k0-1,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(0)*bxg(j,k0-1,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(0)*bxg(j+1,k0-1,l0)
    bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(0)*bxg(j-1,k0,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
    bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(0)*bxg(j-1,k0+1,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
    bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(1)*bxg(j-1,k0-1,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(1)*bxg(j,k0-1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(1)*bxg(j+1,k0-1,l0+1)
    bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(1)*bxg(j-1,k0,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
    bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(1)*bxg(j-1,k0+1,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)
    
    ! Compute By on particle
    by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(-1)*byg(j0-1,k-1,l0-1)
    by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(-1)*byg(j0,k-1,l0-1)
    by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(-1)*byg(j0+1,k-1,l0-1)
    by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(-1)*byg(j0-1,k,l0-1)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(-1)*byg(j0,k,l0-1)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(-1)*byg(j0+1,k,l0-1)
    by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(-1)*byg(j0-1,k+1,l0-1)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(-1)*byg(j0,k+1,l0-1)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(-1)*byg(j0+1,k+1,l0-1)
    by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(0)*byg(j0-1,k-1,l0)
    by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(0)*byg(j0,k-1,l0)
    by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(0)*byg(j0+1,k-1,l0)
    by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(0)*byg(j0-1,k,l0)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
    by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(0)*byg(j0-1,k+1,l0)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
    by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(1)*byg(j0-1,k-1,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(1)*byg(j0,k-1,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(1)*byg(j0+1,k-1,l0+1)
    by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(1)*byg(j0-1,k,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
    by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(1)*byg(j0-1,k+1,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)
    
    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(-1)*bzg(j0-1,k0-1,l-1)
    bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(-1)*bzg(j0,k0-1,l-1)
    bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(-1)*bzg(j0+1,k0-1,l-1)
    bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(-1)*bzg(j0-1,k0,l-1)
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(-1)*bzg(j0,k0,l-1)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(-1)*bzg(j0+1,k0,l-1)
    bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(-1)*bzg(j0-1,k0+1,l-1)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(-1)*bzg(j0,k0+1,l-1)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(-1)*bzg(j0+1,k0+1,l-1)
    bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(0)*bzg(j0-1,k0-1,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(0)*bzg(j0,k0-1,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(0)*bzg(j0+1,k0-1,l)
    bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(0)*bzg(j0-1,k0,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
    bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(0)*bzg(j0-1,k0+1,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
    bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(1)*bzg(j0-1,k0-1,l+1)
    bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(1)*bzg(j0,k0-1,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(1)*bzg(j0+1,k0-1,l+1)
    bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(1)*bzg(j0-1,k0,l+1)
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
    bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(1)*bzg(j0-1,k0+1,l+1)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)
    
  END DO
#if defined _OPENMP && _OPENMP>=201307
			!$OMP END SIMD 
#endif  
  ENDIF
  
  RETURN
END SUBROUTINE getb3d_energy_conserving_2_2_2

! ________________________________________________________________________________________
!> Field gathering TSC (order 2) with gathering of E and B merged in a single loop
!> @brief
!
!> This function is vectorized
!> @details
SUBROUTINE geteb3d_energy_conserving_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
!
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

  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np,nx,ny,nz,nxguard,nyguard,nzguard
  REAL(num), DIMENSION(np)               :: xp,yp,zp,ex,ey,ez,bx,by,bz
  INTEGER(idp)                           :: lvect
  LOGICAL                                :: l_lower_order_in_v 
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg  
  REAL(num)                              :: xmin,ymin,zmin,dx,dy,dz
  INTEGER(isp)                           :: ip, j, k, l 
  INTEGER(isp)                           :: jj, kk, ll, j0, k0, l0
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  INTEGER(isp)                           :: nn              
  REAL(num), DIMENSION(-1:1)             :: sx,sx0
  REAL(num), DIMENSION(-1:1)             :: sy,sy0
  REAL(num), DIMENSION(-1:1)             :: sz,sz0
  REAL(num), PARAMETER                   :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                   :: twothird=2.0_num/3.0_num

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
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif

    DO nn=ip,ip+MIN(lvect,np-ip+1)-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

			! Compute index of particle
			j=nint(x)
			j0=floor(x-0.5_num)
			k=nint(y)
			k0=floor(y-0.5_num)
			l=nint(z)
			l0=floor(z-0.5_num)
		
			xint=x-j
			yint=y-k
			zint=z-l
		
			! Compute shape factors
			xintsq = xint*xint
			sx(-1) = 0.5_num*(0.5_num-xint)**2
			sx( 0) = 0.75_num-xintsq
			sx( 1) = 0.5_num*(0.5_num+xint)**2
		
			yintsq = yint*yint
			sy(-1) = 0.5_num*(0.5_num-yint)**2
			sy( 0) = 0.75_num-yintsq
			sy( 1) = 0.5_num*(0.5_num+yint)**2
		
			zintsq = zint*zint
			sz(-1) = 0.5_num*(0.5_num-zint)**2
			sz( 0) = 0.75_num-zintsq
			sz( 1) = 0.5_num*(0.5_num+zint)**2
		
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
			ex(nn) = ex(nn) + sx0(0)*sy(-1)*sz(-1)*exg(j0,k-1,l-1)
			ex(nn) = ex(nn) + sx0(1)*sy(-1)*sz(-1)*exg(j0+1,k-1,l-1)
			ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(-1)*exg(j0,k,l-1)
			ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(-1)*exg(j0+1,k,l-1)
			ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(-1)*exg(j0,k+1,l-1)
			ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(-1)*exg(j0+1,k+1,l-1)
			ex(nn) = ex(nn) + sx0(0)*sy(-1)*sz(0)*exg(j0,k-1,l)
			ex(nn) = ex(nn) + sx0(1)*sy(-1)*sz(0)*exg(j0+1,k-1,l)
			ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
			ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(0)*exg(j0+1,k,l)
			ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
			ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(0)*exg(j0+1,k+1,l)
			ex(nn) = ex(nn) + sx0(0)*sy(-1)*sz(1)*exg(j0,k-1,l+1)
			ex(nn) = ex(nn) + sx0(1)*sy(-1)*sz(1)*exg(j0+1,k-1,l+1)
			ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
			ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(1)*exg(j0+1,k,l+1)
			ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
			ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(1)*exg(j0+1,k+1,l+1)
    
			! Compute Ey on particle
			ey(nn) = ey(nn) + sx(-1)*sy0(0)*sz(-1)*eyg(j-1,k0,l-1)
			ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(-1)*eyg(j,k0,l-1)
			ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(-1)*eyg(j+1,k0,l-1)
			ey(nn) = ey(nn) + sx(-1)*sy0(1)*sz(-1)*eyg(j-1,k0+1,l-1)
			ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(-1)*eyg(j,k0+1,l-1)
			ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(-1)*eyg(j+1,k0+1,l-1)
			ey(nn) = ey(nn) + sx(-1)*sy0(0)*sz(0)*eyg(j-1,k0,l)
			ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
			ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
			ey(nn) = ey(nn) + sx(-1)*sy0(1)*sz(0)*eyg(j-1,k0+1,l)
			ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(0)*eyg(j,k0+1,l)
			ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(0)*eyg(j+1,k0+1,l)
			ey(nn) = ey(nn) + sx(-1)*sy0(0)*sz(1)*eyg(j-1,k0,l+1)
			ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
			ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
			ey(nn) = ey(nn) + sx(-1)*sy0(1)*sz(1)*eyg(j-1,k0+1,l+1)
			ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(1)*eyg(j,k0+1,l+1)
			ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(1)*eyg(j+1,k0+1,l+1)
    
			! Compute Ez on particle
			ez(nn) = ez(nn) + sx(-1)*sy(-1)*sz0(0)*ezg(j-1,k-1,l0)
			ez(nn) = ez(nn) + sx(0)*sy(-1)*sz0(0)*ezg(j,k-1,l0)
			ez(nn) = ez(nn) + sx(1)*sy(-1)*sz0(0)*ezg(j+1,k-1,l0)
			ez(nn) = ez(nn) + sx(-1)*sy(0)*sz0(0)*ezg(j-1,k,l0)
			ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
			ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
			ez(nn) = ez(nn) + sx(-1)*sy(1)*sz0(0)*ezg(j-1,k+1,l0)
			ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
			ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
			ez(nn) = ez(nn) + sx(-1)*sy(-1)*sz0(1)*ezg(j-1,k-1,l0+1)
			ez(nn) = ez(nn) + sx(0)*sy(-1)*sz0(1)*ezg(j,k-1,l0+1)
			ez(nn) = ez(nn) + sx(1)*sy(-1)*sz0(1)*ezg(j+1,k-1,l0+1)
			ez(nn) = ez(nn) + sx(-1)*sy(0)*sz0(1)*ezg(j-1,k,l0+1)
			ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(1)*ezg(j,k,l0+1)
			ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(1)*ezg(j+1,k,l0+1)
			ez(nn) = ez(nn) + sx(-1)*sy(1)*sz0(1)*ezg(j-1,k+1,l0+1)
			ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(1)*ezg(j,k+1,l0+1)
			ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(1)*ezg(j+1,k+1,l0+1)

			! Compute Bx on particle
			bx(nn) = bx(nn) + sx(-1)*sy0(0)*sz0(0)*bxg(j-1,k0,l0)
			bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
			bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
			bx(nn) = bx(nn) + sx(-1)*sy0(1)*sz0(0)*bxg(j-1,k0+1,l0)
			bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
			bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
			bx(nn) = bx(nn) + sx(-1)*sy0(0)*sz0(1)*bxg(j-1,k0,l0+1)
			bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
			bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
			bx(nn) = bx(nn) + sx(-1)*sy0(1)*sz0(1)*bxg(j-1,k0+1,l0+1)
			bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
			bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)
    
			! Compute By on particle
			by(nn) = by(nn) + sx0(0)*sy(-1)*sz0(0)*byg(j0,k-1,l0)
			by(nn) = by(nn) + sx0(1)*sy(-1)*sz0(0)*byg(j0+1,k-1,l0)
			by(nn) = by(nn) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
			by(nn) = by(nn) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
			by(nn) = by(nn) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
			by(nn) = by(nn) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
			by(nn) = by(nn) + sx0(0)*sy(-1)*sz0(1)*byg(j0,k-1,l0+1)
			by(nn) = by(nn) + sx0(1)*sy(-1)*sz0(1)*byg(j0+1,k-1,l0+1)
			by(nn) = by(nn) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
			by(nn) = by(nn) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
			by(nn) = by(nn) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
			by(nn) = by(nn) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)
    
			! Compute Bz on particle
			bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(-1)*bzg(j0,k0,l-1)
			bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(-1)*bzg(j0+1,k0,l-1)
			bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(-1)*bzg(j0,k0+1,l-1)
			bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(-1)*bzg(j0+1,k0+1,l-1)
			bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
			bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
			bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
			bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
			bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
			bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
			bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
			bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)
      
    ENDDO
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
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif

			DO nn=ip,ip+MIN(lvect,np-ip+1)-1

				x = (xp(nn)-xmin)*dxi
				y = (yp(nn)-ymin)*dyi
				z = (zp(nn)-zmin)*dzi
	
				! Compute index of particle
				j=nint(x)
				j0=floor(x)
				k=nint(y)
				k0=floor(y)
				l=nint(z)
				l0=floor(z)
	
				xint=x-j
				yint=y-k
				zint=z-l

				! Compute shape factors
				xintsq = xint*xint
				sx(-1) = 0.5_num*(0.5_num-xint)**2
				sx( 0) = 0.75_num-xintsq
				sx( 1) = 0.5_num*(0.5_num+xint)**2
	
				yintsq = yint*yint
				sy(-1) = 0.5_num*(0.5_num-yint)**2
				sy( 0) = 0.75_num-yintsq
				sy( 1) = 0.5_num*(0.5_num+yint)**2
	
				zintsq = zint*zint
				sz(-1) = 0.5_num*(0.5_num-zint)**2
				sz( 0) = 0.75_num-zintsq
				sz( 1) = 0.5_num*(0.5_num+zint)**2
	
				xint=x-0.5_num-j0
				yint=y-0.5_num-k0
				zint=z-0.5_num-l0

				xintsq = xint*xint
				sx0(-1) = 0.5_num*(0.5_num-xint)**2
				sx0( 0) = 0.75_num-xintsq
				sx0( 1) = 0.5_num*(0.5_num+xint)**2
	
				yintsq = yint*yint
				sy0(-1) = 0.5_num*(0.5_num-yint)**2
				sy0( 0) = 0.75_num-yintsq
				sy0( 1) = 0.5_num*(0.5_num+yint)**2
	
				zintsq = zint*zint
				sz0(-1) = 0.5_num*(0.5_num-zint)**2
				sz0( 0) = 0.75_num-zintsq
				sz0( 1) = 0.5_num*(0.5_num+zint)**2

				! Compute Ex on particle
				ex(nn) = ex(nn) + sx0(-1)*sy(-1)*sz(-1)*exg(j0-1,k-1,l-1)
				ex(nn) = ex(nn) + sx0(0)*sy(-1)*sz(-1)*exg(j0,k-1,l-1)
				ex(nn) = ex(nn) + sx0(1)*sy(-1)*sz(-1)*exg(j0+1,k-1,l-1)
				ex(nn) = ex(nn) + sx0(-1)*sy(0)*sz(-1)*exg(j0-1,k,l-1)
				ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(-1)*exg(j0,k,l-1)
				ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(-1)*exg(j0+1,k,l-1)
				ex(nn) = ex(nn) + sx0(-1)*sy(1)*sz(-1)*exg(j0-1,k+1,l-1)
				ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(-1)*exg(j0,k+1,l-1)
				ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(-1)*exg(j0+1,k+1,l-1)
				ex(nn) = ex(nn) + sx0(-1)*sy(-1)*sz(0)*exg(j0-1,k-1,l)
				ex(nn) = ex(nn) + sx0(0)*sy(-1)*sz(0)*exg(j0,k-1,l)
				ex(nn) = ex(nn) + sx0(1)*sy(-1)*sz(0)*exg(j0+1,k-1,l)
				ex(nn) = ex(nn) + sx0(-1)*sy(0)*sz(0)*exg(j0-1,k,l)
				ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(0)*exg(j0,k,l)
				ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(0)*exg(j0+1,k,l)
				ex(nn) = ex(nn) + sx0(-1)*sy(1)*sz(0)*exg(j0-1,k+1,l)
				ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(0)*exg(j0,k+1,l)
				ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(0)*exg(j0+1,k+1,l)
				ex(nn) = ex(nn) + sx0(-1)*sy(-1)*sz(1)*exg(j0-1,k-1,l+1)
				ex(nn) = ex(nn) + sx0(0)*sy(-1)*sz(1)*exg(j0,k-1,l+1)
				ex(nn) = ex(nn) + sx0(1)*sy(-1)*sz(1)*exg(j0+1,k-1,l+1)
				ex(nn) = ex(nn) + sx0(-1)*sy(0)*sz(1)*exg(j0-1,k,l+1)
				ex(nn) = ex(nn) + sx0(0)*sy(0)*sz(1)*exg(j0,k,l+1)
				ex(nn) = ex(nn) + sx0(1)*sy(0)*sz(1)*exg(j0+1,k,l+1)
				ex(nn) = ex(nn) + sx0(-1)*sy(1)*sz(1)*exg(j0-1,k+1,l+1)
				ex(nn) = ex(nn) + sx0(0)*sy(1)*sz(1)*exg(j0,k+1,l+1)
				ex(nn) = ex(nn) + sx0(1)*sy(1)*sz(1)*exg(j0+1,k+1,l+1)
	
				! Compute Ey on particle
				ey(nn) = ey(nn) + sx(-1)*sy0(-1)*sz(-1)*eyg(j-1,k0-1,l-1)
				ey(nn) = ey(nn) + sx(0)*sy0(-1)*sz(-1)*eyg(j,k0-1,l-1)
				ey(nn) = ey(nn) + sx(1)*sy0(-1)*sz(-1)*eyg(j+1,k0-1,l-1)
				ey(nn) = ey(nn) + sx(-1)*sy0(0)*sz(-1)*eyg(j-1,k0,l-1)
				ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(-1)*eyg(j,k0,l-1)
				ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(-1)*eyg(j+1,k0,l-1)
				ey(nn) = ey(nn) + sx(-1)*sy0(1)*sz(-1)*eyg(j-1,k0+1,l-1)
				ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(-1)*eyg(j,k0+1,l-1)
				ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(-1)*eyg(j+1,k0+1,l-1)
				ey(nn) = ey(nn) + sx(-1)*sy0(-1)*sz(0)*eyg(j-1,k0-1,l)
				ey(nn) = ey(nn) + sx(0)*sy0(-1)*sz(0)*eyg(j,k0-1,l)
				ey(nn) = ey(nn) + sx(1)*sy0(-1)*sz(0)*eyg(j+1,k0-1,l)
				ey(nn) = ey(nn) + sx(-1)*sy0(0)*sz(0)*eyg(j-1,k0,l)
				ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(0)*eyg(j,k0,l)
				ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(0)*eyg(j+1,k0,l)
				ey(nn) = ey(nn) + sx(-1)*sy0(1)*sz(0)*eyg(j-1,k0+1,l)
				ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(0)*eyg(j,k0+1,l)
				ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(0)*eyg(j+1,k0+1,l)
				ey(nn) = ey(nn) + sx(-1)*sy0(-1)*sz(1)*eyg(j-1,k0-1,l+1)
				ey(nn) = ey(nn) + sx(0)*sy0(-1)*sz(1)*eyg(j,k0-1,l+1)
				ey(nn) = ey(nn) + sx(1)*sy0(-1)*sz(1)*eyg(j+1,k0-1,l+1)
				ey(nn) = ey(nn) + sx(-1)*sy0(0)*sz(1)*eyg(j-1,k0,l+1)
				ey(nn) = ey(nn) + sx(0)*sy0(0)*sz(1)*eyg(j,k0,l+1)
				ey(nn) = ey(nn) + sx(1)*sy0(0)*sz(1)*eyg(j+1,k0,l+1)
				ey(nn) = ey(nn) + sx(-1)*sy0(1)*sz(1)*eyg(j-1,k0+1,l+1)
				ey(nn) = ey(nn) + sx(0)*sy0(1)*sz(1)*eyg(j,k0+1,l+1)
				ey(nn) = ey(nn) + sx(1)*sy0(1)*sz(1)*eyg(j+1,k0+1,l+1)

				! Compute Ez on particle
				ez(nn) = ez(nn) + sx(-1)*sy(-1)*sz0(-1)*ezg(j-1,k-1,l0-1)
				ez(nn) = ez(nn) + sx(0)*sy(-1)*sz0(-1)*ezg(j,k-1,l0-1)
				ez(nn) = ez(nn) + sx(1)*sy(-1)*sz0(-1)*ezg(j+1,k-1,l0-1)
				ez(nn) = ez(nn) + sx(-1)*sy(0)*sz0(-1)*ezg(j-1,k,l0-1)
				ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(-1)*ezg(j,k,l0-1)
				ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(-1)*ezg(j+1,k,l0-1)
				ez(nn) = ez(nn) + sx(-1)*sy(1)*sz0(-1)*ezg(j-1,k+1,l0-1)
				ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(-1)*ezg(j,k+1,l0-1)
				ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(-1)*ezg(j+1,k+1,l0-1)
				ez(nn) = ez(nn) + sx(-1)*sy(-1)*sz0(0)*ezg(j-1,k-1,l0)
				ez(nn) = ez(nn) + sx(0)*sy(-1)*sz0(0)*ezg(j,k-1,l0)
				ez(nn) = ez(nn) + sx(1)*sy(-1)*sz0(0)*ezg(j+1,k-1,l0)
				ez(nn) = ez(nn) + sx(-1)*sy(0)*sz0(0)*ezg(j-1,k,l0)
				ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(0)*ezg(j,k,l0)
				ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(0)*ezg(j+1,k,l0)
				ez(nn) = ez(nn) + sx(-1)*sy(1)*sz0(0)*ezg(j-1,k+1,l0)
				ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(0)*ezg(j,k+1,l0)
				ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(0)*ezg(j+1,k+1,l0)
				ez(nn) = ez(nn) + sx(-1)*sy(-1)*sz0(1)*ezg(j-1,k-1,l0+1)
				ez(nn) = ez(nn) + sx(0)*sy(-1)*sz0(1)*ezg(j,k-1,l0+1)
				ez(nn) = ez(nn) + sx(1)*sy(-1)*sz0(1)*ezg(j+1,k-1,l0+1)
				ez(nn) = ez(nn) + sx(-1)*sy(0)*sz0(1)*ezg(j-1,k,l0+1)
				ez(nn) = ez(nn) + sx(0)*sy(0)*sz0(1)*ezg(j,k,l0+1)
				ez(nn) = ez(nn) + sx(1)*sy(0)*sz0(1)*ezg(j+1,k,l0+1)
				ez(nn) = ez(nn) + sx(-1)*sy(1)*sz0(1)*ezg(j-1,k+1,l0+1)
				ez(nn) = ez(nn) + sx(0)*sy(1)*sz0(1)*ezg(j,k+1,l0+1)
				ez(nn) = ez(nn) + sx(1)*sy(1)*sz0(1)*ezg(j+1,k+1,l0+1)

				! Compute Bx on particle
				bx(nn) = bx(nn) + sx(-1)*sy0(-1)*sz0(-1)*bxg(j-1,k0-1,l0-1)
				bx(nn) = bx(nn) + sx(0)*sy0(-1)*sz0(-1)*bxg(j,k0-1,l0-1)
				bx(nn) = bx(nn) + sx(1)*sy0(-1)*sz0(-1)*bxg(j+1,k0-1,l0-1)
				bx(nn) = bx(nn) + sx(-1)*sy0(0)*sz0(-1)*bxg(j-1,k0,l0-1)
				bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(-1)*bxg(j,k0,l0-1)
				bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(-1)*bxg(j+1,k0,l0-1)
				bx(nn) = bx(nn) + sx(-1)*sy0(1)*sz0(-1)*bxg(j-1,k0+1,l0-1)
				bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(-1)*bxg(j,k0+1,l0-1)
				bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(-1)*bxg(j+1,k0+1,l0-1)
				bx(nn) = bx(nn) + sx(-1)*sy0(-1)*sz0(0)*bxg(j-1,k0-1,l0)
				bx(nn) = bx(nn) + sx(0)*sy0(-1)*sz0(0)*bxg(j,k0-1,l0)
				bx(nn) = bx(nn) + sx(1)*sy0(-1)*sz0(0)*bxg(j+1,k0-1,l0)
				bx(nn) = bx(nn) + sx(-1)*sy0(0)*sz0(0)*bxg(j-1,k0,l0)
				bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
				bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
				bx(nn) = bx(nn) + sx(-1)*sy0(1)*sz0(0)*bxg(j-1,k0+1,l0)
				bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
				bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
				bx(nn) = bx(nn) + sx(-1)*sy0(-1)*sz0(1)*bxg(j-1,k0-1,l0+1)
				bx(nn) = bx(nn) + sx(0)*sy0(-1)*sz0(1)*bxg(j,k0-1,l0+1)
				bx(nn) = bx(nn) + sx(1)*sy0(-1)*sz0(1)*bxg(j+1,k0-1,l0+1)
				bx(nn) = bx(nn) + sx(-1)*sy0(0)*sz0(1)*bxg(j-1,k0,l0+1)
				bx(nn) = bx(nn) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
				bx(nn) = bx(nn) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
				bx(nn) = bx(nn) + sx(-1)*sy0(1)*sz0(1)*bxg(j-1,k0+1,l0+1)
				bx(nn) = bx(nn) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
				bx(nn) = bx(nn) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)

				! Compute By on particle
				by(nn) = by(nn) + sx0(-1)*sy(-1)*sz0(-1)*byg(j0-1,k-1,l0-1)
				by(nn) = by(nn) + sx0(0)*sy(-1)*sz0(-1)*byg(j0,k-1,l0-1)
				by(nn) = by(nn) + sx0(1)*sy(-1)*sz0(-1)*byg(j0+1,k-1,l0-1)
				by(nn) = by(nn) + sx0(-1)*sy(0)*sz0(-1)*byg(j0-1,k,l0-1)
				by(nn) = by(nn) + sx0(0)*sy(0)*sz0(-1)*byg(j0,k,l0-1)
				by(nn) = by(nn) + sx0(1)*sy(0)*sz0(-1)*byg(j0+1,k,l0-1)
				by(nn) = by(nn) + sx0(-1)*sy(1)*sz0(-1)*byg(j0-1,k+1,l0-1)
				by(nn) = by(nn) + sx0(0)*sy(1)*sz0(-1)*byg(j0,k+1,l0-1)
				by(nn) = by(nn) + sx0(1)*sy(1)*sz0(-1)*byg(j0+1,k+1,l0-1)
				by(nn) = by(nn) + sx0(-1)*sy(-1)*sz0(0)*byg(j0-1,k-1,l0)
				by(nn) = by(nn) + sx0(0)*sy(-1)*sz0(0)*byg(j0,k-1,l0)
				by(nn) = by(nn) + sx0(1)*sy(-1)*sz0(0)*byg(j0+1,k-1,l0)
				by(nn) = by(nn) + sx0(-1)*sy(0)*sz0(0)*byg(j0-1,k,l0)
				by(nn) = by(nn) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
				by(nn) = by(nn) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
				by(nn) = by(nn) + sx0(-1)*sy(1)*sz0(0)*byg(j0-1,k+1,l0)
				by(nn) = by(nn) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
				by(nn) = by(nn) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
				by(nn) = by(nn) + sx0(-1)*sy(-1)*sz0(1)*byg(j0-1,k-1,l0+1)
				by(nn) = by(nn) + sx0(0)*sy(-1)*sz0(1)*byg(j0,k-1,l0+1)
				by(nn) = by(nn) + sx0(1)*sy(-1)*sz0(1)*byg(j0+1,k-1,l0+1)
				by(nn) = by(nn) + sx0(-1)*sy(0)*sz0(1)*byg(j0-1,k,l0+1)
				by(nn) = by(nn) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
				by(nn) = by(nn) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
				by(nn) = by(nn) + sx0(-1)*sy(1)*sz0(1)*byg(j0-1,k+1,l0+1)
				by(nn) = by(nn) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
				by(nn) = by(nn) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)

				! Compute Bz on particle
				bz(nn) = bz(nn) + sx0(-1)*sy0(-1)*sz(-1)*bzg(j0-1,k0-1,l-1)
				bz(nn) = bz(nn) + sx0(0)*sy0(-1)*sz(-1)*bzg(j0,k0-1,l-1)
				bz(nn) = bz(nn) + sx0(1)*sy0(-1)*sz(-1)*bzg(j0+1,k0-1,l-1)
				bz(nn) = bz(nn) + sx0(-1)*sy0(0)*sz(-1)*bzg(j0-1,k0,l-1)
				bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(-1)*bzg(j0,k0,l-1)
				bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(-1)*bzg(j0+1,k0,l-1)
				bz(nn) = bz(nn) + sx0(-1)*sy0(1)*sz(-1)*bzg(j0-1,k0+1,l-1)
				bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(-1)*bzg(j0,k0+1,l-1)
				bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(-1)*bzg(j0+1,k0+1,l-1)
				bz(nn) = bz(nn) + sx0(-1)*sy0(-1)*sz(0)*bzg(j0-1,k0-1,l)
				bz(nn) = bz(nn) + sx0(0)*sy0(-1)*sz(0)*bzg(j0,k0-1,l)
				bz(nn) = bz(nn) + sx0(1)*sy0(-1)*sz(0)*bzg(j0+1,k0-1,l)
				bz(nn) = bz(nn) + sx0(-1)*sy0(0)*sz(0)*bzg(j0-1,k0,l)
				bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
				bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
				bz(nn) = bz(nn) + sx0(-1)*sy0(1)*sz(0)*bzg(j0-1,k0+1,l)
				bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
				bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
				bz(nn) = bz(nn) + sx0(-1)*sy0(-1)*sz(1)*bzg(j0-1,k0-1,l+1)
				bz(nn) = bz(nn) + sx0(0)*sy0(-1)*sz(1)*bzg(j0,k0-1,l+1)
				bz(nn) = bz(nn) + sx0(1)*sy0(-1)*sz(1)*bzg(j0+1,k0-1,l+1)
				bz(nn) = bz(nn) + sx0(-1)*sy0(0)*sz(1)*bzg(j0-1,k0,l+1)
				bz(nn) = bz(nn) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
				bz(nn) = bz(nn) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
				bz(nn) = bz(nn) + sx0(-1)*sy0(1)*sz(1)*bzg(j0-1,k0+1,l+1)
				bz(nn) = bz(nn) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
				bz(nn) = bz(nn) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)

			ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif  
		ENDDO  
  
  ENDIF
  
  RETURN     
END SUBROUTINE geteb3d_energy_conserving_2_2_2
! ________________________________________________________________________________________
! 
! FIELD_GATHERING_3D_O2.F90
! 
! Field gathering subroutines in 3D at order 2
!
! List of subroutines:
!
! 
! - gete3d_energy_conserving_scalar_2_2_2
! - getb3d_energy_conserving_scalar_2_2_2
!
! - gete3d_energy_conserving_2_2_2
! - getb3d_energy_conserving_2_2_2
!
! - geteb3d_energy_conserving_2_2_2
! 
! - geteb3d_energy_conserving_vec_2_2_2
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Scalar version: gathering of electric field from Yee grid ("energy conserving")
!> on particles at order 2.
!
!> @details
!> This subroutine is NOT vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position
!> @param[inout] ex,ey,ez particle electric field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction 
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
SUBROUTINE gete3d_energy_conserving_scalar_2_2_2(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
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
  REAL(num), DIMENSION(-1:1)           :: sx, sx0
  REAL(num), DIMENSION(-1:1)           :: sy, sy0
  REAL(num), DIMENSION(-1:1)           :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  sx0 = 0
  sy0 = 0
  sz0 = 0
  sx = 0
  sy = 0
  sz = 0

  ixmin = -1
  ixmax = 0
  iymin = -1
  iymax = 0
  izmin = -1
  izmax = 0

  IF (l_lower_order_in_v) THEN

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

      j=nint(x)
      j0=floor(x-0.5_num)

      k=nint(y)
      k0=floor(y-0.5_num)

      l=nint(z)
      l0=floor(z-0.5_num)

      xint=x-j
      yint=y-k
      zint=z-l

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

      DO ll = izmin, izmax+1
          DO kk = iymin, iymax+1
              DO jj = ixmin0, ixmax0
                  ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)
              END DO
          END DO
      END DO

      DO ll = izmin, izmax+1
          DO kk = iymin0, iymax0
              DO jj = ixmin, ixmax+1
                  ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)
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

  ! __ l_lower_order_in_v false  _____________________________
  ELSE

    ixmin0 = -1
    ixmax0 = 1
    iymin0 = -1
    iymax0 = 1
    izmin0 = -1
    izmax0 = 1

    DO ip=1,np
    
      x = (xp(ip)-xmin)*dxi
      y = (yp(ip)-ymin)*dyi
      z = (zp(ip)-zmin)*dzi

      j=nint(x)
      j0=floor(x)

      k=nint(y)
      k0=floor(y)

      l=nint(z)
      l0=floor(z)

      xint=x-j
      yint=y-k
      zint=z-l

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

      do ll = izmin, izmax+1
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)
          end do
        end do
      end do

      do ll = izmin, izmax+1
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)
          end do
        end do
      end do

      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj,k+kk,l0+ll)
          end do
        end do
      end do
      
    END DO
  ENDIF

RETURN
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> Scalar version: Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 2
!
!> @details
!> This function is vectorized 
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/05/2016
!
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position arrays
!> @param[inout] bx,by,bz particle magnetic field arrays
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space steps in every directions
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction 
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
SUBROUTINE getb3d_energy_conserving_scalar_2_2_2(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
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
  REAL(num), DIMENSION(-1:1)            :: sx, sx0
  REAL(num), DIMENSION(-1:1)            :: sy, sy0
  REAL(num), DIMENSION(-1:1)            :: sz, sz0
  REAL(num), PARAMETER                  :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  ixmin = -1
  ixmax = 0
  iymin = -1
  iymax = 0
  izmin = -1
  izmax = 0

  IF (l_lower_order_in_v) THEN

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

      j=nint(x)
      j0=floor(x-0.5_num)

      k=nint(y)
      k0=floor(y-0.5_num)

      l=nint(z)
      l0=floor(z-0.5_num)

      xint=x-j
      yint=y-k
      zint=z-l

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

      do ll = izmin0, izmax0
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)
          end do
        end do
      end do

      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)
          end do
        end do
      end do

      do ll = izmin, izmax+1
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj,k0+kk,l+ll)
          end do
        end do
      end do
  
    END DO

  ELSE

    ixmin0 = -1
    ixmax0 = 1
    iymin0 = -1
    iymax0 = 1
    izmin0 = -1
    izmax0 = 1

    DO ip=1,np
    
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        j=nint(x)
        j0=floor(x)

        k=nint(y)
        k0=floor(y)

        l=nint(z)
        l0=floor(z)

        xint=x-j
        yint=y-k
        zint=z-l

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

    
        do ll = izmin0, izmax0
          do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)
            end do
          end do
        end do

        do ll = izmin0, izmax0
          do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)
            end do
          end do
        end do

        do ll = izmin, izmax+1
          do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
            do jj = ixmin0, ixmax0
              bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj,k0+kk,l+ll)
            end do
          end do
        end do
    END DO
  ENDIF
  RETURN
END SUBROUTINE


#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Gathering of electric field from Yee grid ("energy conserving") on particles
!> at order 2.
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/05/2016
!
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position
!> @param[inout] ex,ey,ez particle electric field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction 
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
SUBROUTINE gete3d_energy_conserving_vec_2_2_2(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,&
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________
  USE omp_lib
  USE constants
  IMPLICIT NONE
  
  INTEGER(idp)                         :: np,nx,ny,nz
  INTEGER(idp)                         :: nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp,yp,zp,ex,ey,ez
  LOGICAL                              :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: jj, kk, ll
  INTEGER(isp)                         :: j0, k0, l0
  INTEGER(isp)                         :: n,nn
  REAL(num)                            :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                            :: xintsq,oxint,yintsq,oyint,zintsq
  REAL(num)                            :: ozint,oxintsq,oyintsq,ozintsq
  REAL(num), DIMENSION(lvect,-1:1)     :: sx,sx0
  REAL(num), DIMENSION(lvect,-1:1)     :: sy,sy0
  REAL(num), DIMENSION(lvect,-1:1)     :: sz,sz0
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

  ! ___ Loop on partciles _______________________
  DO ip=1,np,lvect

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
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif

    ! Loop over the particles inside a block
    DO n=1,MIN(lvect,np-ip+1)

      nn=ip+n-1

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
      sx(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx(n, 0) = 0.75_num-xintsq
      sx(n, 1) = 0.5_num*(0.5_num+xint)**2
    
      yintsq = yint*yint
      sy(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy(n, 0) = 0.75_num-yintsq
      sy(n, 1) = 0.5_num*(0.5_num+yint)**2
    
      zintsq = zint*zint
      sz(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz(n, 0) = 0.75_num-zintsq
      sz(n, 1) = 0.5_num*(0.5_num+zint)**2
    
      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0
    
      sx0(n, 0) = 1.0_num-xint
      sx0(n, 1) = xint
    
      sy0(n, 0) = 1.0_num-yint
      sy0(n, 1) = yint
    
      sz0(n, 0) = 1.0_num-zint
      sz0(n, 1) = zint

      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0,k-1,l-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0+1,k-1,l-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0,k,l-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0+1,k,l-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0,k+1,l-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0+1,k+1,l-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0,k-1,l)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0+1,k-1,l)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0,k,l)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0+1,k,l)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0,k+1,l)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0+1,k+1,l)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0,k-1,l+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0+1,k-1,l+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0,k,l+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0+1,k,l+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0,k+1,l+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0+1,k+1,l+1)
    
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j-1,k0,l-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j,k0,l-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j+1,k0,l-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j-1,k0+1,l-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j,k0+1,l-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j+1,k0+1,l-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j-1,k0,l)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j,k0,l)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j+1,k0,l)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j-1,k0+1,l)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j,k0+1,l)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j+1,k0+1,l)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j-1,k0,l+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j,k0,l+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j+1,k0,l+1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j-1,k0+1,l+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j,k0+1,l+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j+1,k0+1,l+1)
    
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j-1,k-1,l0)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j,k-1,l0)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j+1,k-1,l0)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j-1,k,l0)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j,k,l0)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j+1,k,l0)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j-1,k+1,l0)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j,k+1,l0)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j+1,k+1,l0)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j-1,k-1,l0+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j,k-1,l0+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j+1,k-1,l0+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j-1,k,l0+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j,k,l0+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j+1,k,l0+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j-1,k+1,l0+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j,k+1,l0+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j+1,k+1,l0+1)

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
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1
  
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
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy(n, 0) = 0.75_num-yintsq
        sy(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2
  
        xint=x-0.5_num-j0
        yint=y-0.5_num-k0
        zint=z-0.5_num-l0

        xintsq = xint*xint
        sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy0(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,-1)*exg(j0-1,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0+1,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,-1)*exg(j0-1,k,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0,k,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0+1,k,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,-1)*exg(j0-1,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0+1,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,0)*exg(j0-1,k-1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0,k-1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0+1,k-1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,0)*exg(j0-1,k,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0,k,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0+1,k,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,0)*exg(j0-1,k+1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0,k+1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0+1,k+1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,1)*exg(j0-1,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0+1,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,1)*exg(j0-1,k,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0,k,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0+1,k,l+1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,1)*exg(j0-1,k+1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0,k+1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0+1,k+1,l+1)
    
        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,-1)*eyg(j-1,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,-1)*eyg(j,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,-1)*eyg(j+1,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j-1,k0,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j,k0,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j+1,k0,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j-1,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j+1,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,0)*eyg(j-1,k0-1,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,0)*eyg(j,k0-1,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,0)*eyg(j+1,k0-1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j-1,k0,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j,k0,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j+1,k0,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j-1,k0+1,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j,k0+1,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j+1,k0+1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,1)*eyg(j-1,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,1)*eyg(j,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,1)*eyg(j+1,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j-1,k0,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j,k0,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j+1,k0,l+1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j-1,k0+1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j,k0+1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j+1,k0+1,l+1)
    
        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,-1)*ezg(j-1,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,-1)*ezg(j,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,-1)*ezg(j+1,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,-1)*ezg(j-1,k,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,-1)*ezg(j,k,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,-1)*ezg(j+1,k,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,-1)*ezg(j-1,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,-1)*ezg(j,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,-1)*ezg(j+1,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j-1,k-1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j,k-1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j+1,k-1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j-1,k,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j,k,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j+1,k,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j-1,k+1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j,k+1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j+1,k+1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j-1,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j+1,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j-1,k,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j,k,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j+1,k,l0+1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j-1,k+1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j,k+1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j+1,k+1,l0+1)
  
      END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif
    ENDDO
  ENDIF

RETURN
END SUBROUTINE
#endif



#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 2
!  
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/05/2016
!
! Input parameters:
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position
!> @param[inout] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
SUBROUTINE getb3d_energy_conserving_vec_2_2_2(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,       &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,lvect,l_lower_order_in_v)
! ________________________________________________________________________________________
                                    
  USE omp_lib
  USE constants
  IMPLICIT NONE
  INTEGER(idp)                         :: np,nx,ny,nz,nxguard,nyguard,nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp,yp,zp,bx,by,bz
  LOGICAL                              :: l_lower_order_in_v
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
  REAL(num)                            :: xmin,ymin,zmin,dx,dy,dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
  INTEGER(isp)                         :: n,nn
  REAL(num)                            :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                            :: xintsq,oxint,yintsq,oyint,zintsq
  REAL(num)                            :: ozint,oxintsq,oyintsq,ozintsq
  REAL(num), DIMENSION(lvect,-1:1)     :: sx,sx0
  REAL(num), DIMENSION(lvect,-1:1)     :: sy,sy0
  REAL(num), DIMENSION(lvect,-1:1)     :: sz,sz0
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

  ! ___ Loop on partciles _______________________
  DO ip=1,np,lvect

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
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif  

    ! Loop over the particles inside a block
    DO n=1,MIN(lvect,np-ip+1)

      nn=ip+n-1

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
      sx(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx(n, 0) = 0.75_num-xintsq
      sx(n, 1) = 0.5_num*(0.5_num+xint)**2
    
      yintsq = yint*yint
      sy(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy(n, 0) = 0.75_num-yintsq
      sy(n, 1) = 0.5_num*(0.5_num+yint)**2
    
      zintsq = zint*zint
      sz(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz(n, 0) = 0.75_num-zintsq
      sz(n, 1) = 0.5_num*(0.5_num+zint)**2
    
      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0
    
      sx0(n, 0) = 1.0_num-xint
      sx0(n, 1) = xint
    
      sy0(n, 0) = 1.0_num-yint
      sy0(n, 1) = yint
    
      sz0(n, 0) = 1.0_num-zint
      sz0(n, 1) = zint
    
      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j-1,k0,l0)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j,k0,l0)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j+1,k0,l0)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j-1,k0+1,l0)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j,k0+1,l0)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j+1,k0+1,l0)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j-1,k0,l0+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j,k0,l0+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j+1,k0,l0+1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j-1,k0+1,l0+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j,k0+1,l0+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j+1,k0+1,l0+1)
    
      ! Compute By on particle
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0,k-1,l0)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0+1,k-1,l0)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0,k,l0)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0+1,k,l0)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0,k+1,l0)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0+1,k+1,l0)
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0,k-1,l0+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0+1,k-1,l0+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0,k,l0+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0+1,k,l0+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0,k+1,l0+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0+1,k+1,l0+1)
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0,k0,l-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0+1,k0,l-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0,k0+1,l-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0+1,k0+1,l-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0,k0,l)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0+1,k0,l)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0,k0+1,l)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0+1,k0+1,l)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0,k0,l+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0+1,k0,l+1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0,k0+1,l+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0+1,k0+1,l+1)


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
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif  
      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1
    
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
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy(n, 0) = 0.75_num-yintsq
        sy(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2
  
        xint=x-0.5_num-j0
        yint=y-0.5_num-k0
        zint=z-0.5_num-l0

        xintsq = xint*xint
        sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy0(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2
  
        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,-1)*bxg(j-1,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,-1)*bxg(j,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,-1)*bxg(j+1,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,-1)*bxg(j-1,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,-1)*bxg(j,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,-1)*bxg(j+1,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,-1)*bxg(j-1,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,-1)*bxg(j,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,-1)*bxg(j+1,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,0)*bxg(j-1,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,0)*bxg(j,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,0)*bxg(j+1,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j-1,k0,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j,k0,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j+1,k0,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j-1,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j+1,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,1)*bxg(j-1,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,1)*bxg(j,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,1)*bxg(j+1,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j-1,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j+1,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j-1,k0+1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j,k0+1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j+1,k0+1,l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,-1)*byg(j0-1,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,-1)*byg(j0,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,-1)*byg(j0+1,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,-1)*byg(j0-1,k,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,-1)*byg(j0,k,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,-1)*byg(j0+1,k,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,-1)*byg(j0-1,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,-1)*byg(j0,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,-1)*byg(j0+1,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,0)*byg(j0-1,k-1,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0,k-1,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0+1,k-1,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,0)*byg(j0-1,k,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0,k,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0+1,k,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,0)*byg(j0-1,k+1,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0,k+1,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0+1,k+1,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,1)*byg(j0-1,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0+1,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,1)*byg(j0-1,k,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0,k,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0+1,k,l0+1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,1)*byg(j0-1,k+1,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0,k+1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0+1,k+1,l0+1)
    
        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,-1)*bzg(j0-1,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,-1)*bzg(j0,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,-1)*bzg(j0+1,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,-1)*bzg(j0-1,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0+1,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,-1)*bzg(j0-1,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0+1,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,0)*bzg(j0-1,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,0)*bzg(j0,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,0)*bzg(j0+1,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,0)*bzg(j0-1,k0,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0,k0,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0+1,k0,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,0)*bzg(j0-1,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0+1,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,1)*bzg(j0-1,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,1)*bzg(j0,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,1)*bzg(j0+1,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,1)*bzg(j0-1,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0+1,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,1)*bzg(j0-1,k0+1,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0,k0+1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0+1,k0+1,l+1)
    
      END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif  
    ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE 
#endif

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Field gathering TSC (order 2) with gathering of E and B merged in a single loop
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/05/2016
!
! Input parameters:
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position
!> @param[inout] ex,ey,ez particle electric field
!> @param[inout] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction 
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
SUBROUTINE geteb3d_energy_conserving_vecV1_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
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
  INTEGER(isp)                           :: nn,n 
  REAL(num), DIMENSION(lvect,-1:1)       :: sx,sx0
  REAL(num), DIMENSION(lvect,-1:1)       :: sy,sy0
  REAL(num), DIMENSION(lvect,-1:1)       :: sz,sz0
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

    ! Loop over the particles inside a block
    DO n=1,MIN(lvect,np-ip+1)

      nn=ip+n-1

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
      sx(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx(n, 0) = 0.75_num-xintsq
      sx(n, 1) = 0.5_num*(0.5_num+xint)**2
    
      yintsq = yint*yint
      sy(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy(n, 0) = 0.75_num-yintsq
      sy(n, 1) = 0.5_num*(0.5_num+yint)**2
    
      zintsq = zint*zint
      sz(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz(n, 0) = 0.75_num-zintsq
      sz(n, 1) = 0.5_num*(0.5_num+zint)**2
    
      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0
    
      sx0(n, 0) = 1.0_num-xint
      sx0(n, 1) = xint
    
      sy0(n, 0) = 1.0_num-yint
      sy0(n, 1) = yint
    
      sz0(n, 0) = 1.0_num-zint
      sz0(n, 1) = zint

      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0,k-1,l-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0+1,k-1,l-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0,k,l-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0+1,k,l-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0,k+1,l-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0+1,k+1,l-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0,k-1,l)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0+1,k-1,l)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0,k,l)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0+1,k,l)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0,k+1,l)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0+1,k+1,l)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0,k-1,l+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0+1,k-1,l+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0,k,l+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0+1,k,l+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0,k+1,l+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0+1,k+1,l+1)
    
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j-1,k0,l-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j,k0,l-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j+1,k0,l-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j-1,k0+1,l-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j,k0+1,l-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j+1,k0+1,l-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j-1,k0,l)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j,k0,l)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j+1,k0,l)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j-1,k0+1,l)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j,k0+1,l)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j+1,k0+1,l)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j-1,k0,l+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j,k0,l+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j+1,k0,l+1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j-1,k0+1,l+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j,k0+1,l+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j+1,k0+1,l+1)
    
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j-1,k-1,l0)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j,k-1,l0)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j+1,k-1,l0)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j-1,k,l0)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j,k,l0)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j+1,k,l0)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j-1,k+1,l0)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j,k+1,l0)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j+1,k+1,l0)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j-1,k-1,l0+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j,k-1,l0+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j+1,k-1,l0+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j-1,k,l0+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j,k,l0+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j+1,k,l0+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j-1,k+1,l0+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j,k+1,l0+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j+1,k+1,l0+1)
    
      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j-1,k0,l0)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j,k0,l0)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j+1,k0,l0)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j-1,k0+1,l0)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j,k0+1,l0)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j+1,k0+1,l0)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j-1,k0,l0+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j,k0,l0+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j+1,k0,l0+1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j-1,k0+1,l0+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j,k0+1,l0+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j+1,k0+1,l0+1)
    
      ! Compute By on particle
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0,k-1,l0)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0+1,k-1,l0)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0,k,l0)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0+1,k,l0)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0,k+1,l0)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0+1,k+1,l0)
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0,k-1,l0+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0+1,k-1,l0+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0,k,l0+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0+1,k,l0+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0,k+1,l0+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0+1,k+1,l0+1)
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0,k0,l-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0+1,k0,l-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0,k0+1,l-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0+1,k0+1,l-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0,k0,l)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0+1,k0,l)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0,k0+1,l)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0+1,k0+1,l)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0,k0,l+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0+1,k0,l+1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0,k0+1,l+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0+1,k0+1,l+1)
      
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

      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1

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
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy(n, 0) = 0.75_num-yintsq
        sy(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2
  
        xint=x-0.5_num-j0
        yint=y-0.5_num-k0
        zint=z-0.5_num-l0

        xintsq = xint*xint
        sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy0(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,-1)*exg(j0-1,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0+1,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,-1)*exg(j0-1,k,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0,k,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0+1,k,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,-1)*exg(j0-1,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0+1,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,0)*exg(j0-1,k-1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0,k-1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0+1,k-1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,0)*exg(j0-1,k,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0,k,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0+1,k,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,0)*exg(j0-1,k+1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0,k+1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0+1,k+1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,1)*exg(j0-1,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0+1,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,1)*exg(j0-1,k,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0,k,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0+1,k,l+1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,1)*exg(j0-1,k+1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0,k+1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0+1,k+1,l+1)
    
        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,-1)*eyg(j-1,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,-1)*eyg(j,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,-1)*eyg(j+1,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j-1,k0,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j,k0,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j+1,k0,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j-1,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j+1,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,0)*eyg(j-1,k0-1,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,0)*eyg(j,k0-1,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,0)*eyg(j+1,k0-1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j-1,k0,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j,k0,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j+1,k0,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j-1,k0+1,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j,k0+1,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j+1,k0+1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,1)*eyg(j-1,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,1)*eyg(j,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,1)*eyg(j+1,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j-1,k0,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j,k0,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j+1,k0,l+1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j-1,k0+1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j,k0+1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j+1,k0+1,l+1)
    
        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,-1)*ezg(j-1,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,-1)*ezg(j,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,-1)*ezg(j+1,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,-1)*ezg(j-1,k,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,-1)*ezg(j,k,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,-1)*ezg(j+1,k,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,-1)*ezg(j-1,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,-1)*ezg(j,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,-1)*ezg(j+1,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j-1,k-1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j,k-1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j+1,k-1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j-1,k,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j,k,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j+1,k,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j-1,k+1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j,k+1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j+1,k+1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j-1,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j+1,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j-1,k,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j,k,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j+1,k,l0+1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j-1,k+1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j,k+1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j+1,k+1,l0+1)

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,-1)*bxg(j-1,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,-1)*bxg(j,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,-1)*bxg(j+1,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,-1)*bxg(j-1,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,-1)*bxg(j,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,-1)*bxg(j+1,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,-1)*bxg(j-1,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,-1)*bxg(j,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,-1)*bxg(j+1,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,0)*bxg(j-1,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,0)*bxg(j,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,0)*bxg(j+1,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j-1,k0,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j,k0,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j+1,k0,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j-1,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j+1,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,1)*bxg(j-1,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,1)*bxg(j,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,1)*bxg(j+1,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j-1,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j+1,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j-1,k0+1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j,k0+1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j+1,k0+1,l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,-1)*byg(j0-1,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,-1)*byg(j0,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,-1)*byg(j0+1,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,-1)*byg(j0-1,k,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,-1)*byg(j0,k,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,-1)*byg(j0+1,k,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,-1)*byg(j0-1,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,-1)*byg(j0,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,-1)*byg(j0+1,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,0)*byg(j0-1,k-1,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0,k-1,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0+1,k-1,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,0)*byg(j0-1,k,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0,k,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0+1,k,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,0)*byg(j0-1,k+1,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0,k+1,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0+1,k+1,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,1)*byg(j0-1,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0+1,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,1)*byg(j0-1,k,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0,k,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0+1,k,l0+1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,1)*byg(j0-1,k+1,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0,k+1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0+1,k+1,l0+1)
    
        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,-1)*bzg(j0-1,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,-1)*bzg(j0,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,-1)*bzg(j0+1,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,-1)*bzg(j0-1,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0+1,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,-1)*bzg(j0-1,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0+1,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,0)*bzg(j0-1,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,0)*bzg(j0,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,0)*bzg(j0+1,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,0)*bzg(j0-1,k0,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0,k0,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0+1,k0,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,0)*bzg(j0-1,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0+1,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,1)*bzg(j0-1,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,1)*bzg(j0,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,1)*bzg(j0+1,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,1)*bzg(j0-1,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0+1,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,1)*bzg(j0-1,k0+1,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0,k0+1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0+1,k0+1,l+1)

      ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif  
    ENDDO  
  
  ENDIF
  
  RETURN     
END SUBROUTINE
#endif

! ________________________________________________________________________________________
!
!> @brief
!> Field gathering (order 2) with gathering of E and B merged in a single loop
!
!> @detail
!> This function is vectorized but not very efficient.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/05/2016
!
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position
!> @param[inout] ex,ey,ez particle electric field
!> @param[inout] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction 
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
SUBROUTINE geteb3d_energy_conserving_vecV2_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz, &
                                           xmin,ymin,zmin,       &
                                           dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                           exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v )
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  IMPLICIT NONE

  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np,nx,ny,nz,nxguard,nyguard,nzguard
  REAL(num), DIMENSION(np)               :: xp,yp,zp,ex,ey,ez,bx,by,bz
  INTEGER(idp)                           :: lvect
  LOGICAL                                :: l_lower_order_in_v 
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg  
  REAL(num)                              :: xmin,ymin,zmin,dx,dy,dz
  INTEGER(isp)                           :: ip
  INTEGER(isp), DIMENSION(lvect)         :: j, k, l
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j,k,l
#endif  
  INTEGER(isp), DIMENSION(lvect)         :: j0, k0, l0
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: j0,k0,l0
#endif    
  INTEGER(idp)                           :: jj, kk, ll
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
  INTEGER(isp)                           :: nn,n          
  REAL(num), DIMENSION(lvect,-1:1)       :: sx,sy,sz
  REAL(num), DIMENSION(lvect,-1:1)       :: sx0,sy0,sz0
#if defined __INTEL_COMPILER 
    !dir$ attributes align:64 :: sx,sy,sz,sx0,sy0,sz0
#endif    
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

  ! ___ Loop on partciles _______________________
  DO ip=1,np,lvect
  
#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!!DIR DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
    
      ! Compute index of particle
      j(n)=nint(x)
      j0(n)=floor(x-0.5_num)
      k(n)=nint(y)
      k0(n)=floor(y-0.5_num)
      l(n)=nint(z)
      l0(n)=floor(z-0.5_num)
    
      xint=x-j(n)
      yint=y-k(n)
      zint=z-l(n)
      
      ! Compute shape factors
      xintsq = xint*xint
      sx(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx(n, 0) = 0.75_num-xintsq
      sx(n, 1) = 0.5_num*(0.5_num+xint)**2
    
      yintsq = yint*yint
      sy(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy(n, 0) = 0.75_num-yintsq
      sy(n, 1) = 0.5_num*(0.5_num+yint)**2
    
      zintsq = zint*zint
      sz(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz(n, 0) = 0.75_num-zintsq
      sz(n, 1) = 0.5_num*(0.5_num+zint)**2
    
      xint=x-0.5_num-j0(n)
      yint=y-0.5_num-k0(n)
      zint=z-0.5_num-l0(n)
    
      sx0(n, 0) = 1.0_num-xint
      sx0(n, 1) = xint
    
      sy0(n, 0) = 1.0_num-yint
      sy0(n, 1) = yint
    
      sz0(n, 0) = 1.0_num-zint
      sz0(n, 1) = zint
      
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif  

#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64  
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,ex,ey,ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1
      
      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0(n),k(n)-1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0(n)+1,k(n)-1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0(n),k(n),l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0(n)+1,k(n),l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0(n),k(n)+1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0(n)+1,k(n)+1,l(n)-1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0(n),k(n)-1,l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0(n)+1,k(n)-1,l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0(n),k(n),l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0(n)+1,k(n),l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0(n),k(n)+1,l(n))
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0(n)+1,k(n)+1,l(n))
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0(n),k(n)-1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0(n)+1,k(n)-1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0(n),k(n),l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0(n)+1,k(n),l(n)+1)
      ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0(n),k(n)+1,l(n)+1)
      ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0(n)+1,k(n)+1,l(n)+1)
    
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j(n)-1,k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j(n),k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j(n)+1,k0(n),l(n)-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j(n)-1,k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j(n),k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j(n)+1,k0(n)+1,l(n)-1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j(n)-1,k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j(n),k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j(n)+1,k0(n),l(n))
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j(n)-1,k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j(n),k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j(n)+1,k0(n)+1,l(n))
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j(n)-1,k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j(n),k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j(n)+1,k0(n),l(n)+1)
      ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j(n)-1,k0(n)+1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j(n),k0(n)+1,l(n)+1)
      ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j(n)+1,k0(n)+1,l(n)+1)
    
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j(n)-1,k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j(n),k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j(n)+1,k(n)-1,l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j(n)-1,k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j(n),k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j(n)+1,k(n),l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j(n)-1,k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j(n),k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j(n)+1,k(n)+1,l0(n))
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j(n)-1,k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j(n),k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j(n)+1,k(n)-1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j(n)-1,k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j(n),k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j(n)+1,k(n),l0(n)+1)
      ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j(n)-1,k(n)+1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j(n),k(n)+1,l0(n)+1)
      ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j(n)+1,k(n)+1,l0(n)+1)
      
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif  


#if defined __INTEL_COMPILER 
      !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64   
      !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64     
      !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
      !DIR$ ASSUME_ALIGNED j:64,k:64,l:64
      !DIR$ ASSUME_ALIGNED j0:64,k0:64,l0:64         
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,bx,by,bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
    DO n=1,MIN(lvect,np-ip+1)
      nn=ip+n-1

      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j(n)-1,k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j(n),k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j(n)+1,k0(n),l0(n))
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j(n)-1,k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j(n),k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j(n)+1,k0(n)+1,l0(n))
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j(n)-1,k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j(n),k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j(n)+1,k0(n),l0(n)+1)
      bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j(n)-1,k0(n)+1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j(n),k0(n)+1,l0(n)+1)
      bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j(n)+1,k0(n)+1,l0(n)+1)
    
      ! Compute By on particle
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0(n),k(n)-1,l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0(n)+1,k(n)-1,l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0(n),k(n),l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0(n)+1,k(n),l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0(n),k(n)+1,l0(n))
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0(n)+1,k(n)+1,l0(n))
      by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0(n),k(n)-1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0(n)+1,k(n)-1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0(n),k(n),l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0(n)+1,k(n),l0(n)+1)
      by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0(n),k(n)+1,l0(n)+1)
      by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0(n)+1,k(n)+1,l0(n)+1)
    
      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0(n),k0(n),l(n)-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0(n)+1,k0(n),l(n)-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0(n),k0(n)+1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0(n)+1,k0(n)+1,l(n)-1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0(n),k0(n),l(n))
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0(n)+1,k0(n),l(n))
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0(n),k0(n)+1,l(n))
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0(n)+1,k0(n)+1,l(n))
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0(n),k0(n),l(n)+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0(n)+1,k0(n),l(n)+1)
      bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0(n),k0(n)+1,l(n)+1)
      bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0(n)+1,k0(n)+1,l(n)+1)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif 

  ENDDO
  
  RETURN 
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> Field gathering TSC (order 2) version 3with gathering of E and B merged 
!> in a single loop with an efficient vectorization.
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/05/2016
!
! Input parameters:
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position
!> @param[inout] ex,ey,ez particle electric field
!> @param[inout] bx,by,bz particle magnetic field
!> @param[in] xmin,ymin,zmin tile minimum grid position
!> @param[in] dx,dy,dz space step
!> @param[in] dt time step
!> @param[in] nx,ny,nz number of grid points in each direction
!> @param[in] nxguard,nyguard,nzguard number of guard cells in each direction 
!> @param[in] exg,eyg,ezg electric field grid
!> @param[in] bxg,byg,bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
SUBROUTINE geteb3d_energy_conserving_vecV3_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
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
  REAL(num)                              :: xintsq,oxint,yintsq,oyint,zintsq
  REAL(num)                              :: ozint,oxintsq,oyintsq,ozintsq
  REAL(num)                              :: a
  INTEGER(isp)                           :: nn,n 
  REAL(num), DIMENSION(lvect,-1:1)       :: sx,sx0
  REAL(num), DIMENSION(lvect,-1:1)       :: sy,sy0
  REAL(num), DIMENSION(lvect,-1:1)       :: sz,sz0
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

    ! Loop over the particles inside a block
    DO n=1,MIN(lvect,np-ip+1)

      nn=ip+n-1

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
      sx(n,-1) = 0.5_num*(0.5_num-xint)**2
      sx(n, 0) = 0.75_num-xintsq
      sx(n, 1) = 0.5_num*(0.5_num+xint)**2
    
      yintsq = yint*yint
      sy(n,-1) = 0.5_num*(0.5_num-yint)**2
      sy(n, 0) = 0.75_num-yintsq
      sy(n, 1) = 0.5_num*(0.5_num+yint)**2
    
      zintsq = zint*zint
      sz(n,-1) = 0.5_num*(0.5_num-zint)**2
      sz(n, 0) = 0.75_num-zintsq
      sz(n, 1) = 0.5_num*(0.5_num+zint)**2
    
      xint=x-0.5_num-j0
      yint=y-0.5_num-k0
      zint=z-0.5_num-l0
    
      sx0(n, 0) = 1.0_num-xint
      sx0(n, 1) = xint
    
      sy0(n, 0) = 1.0_num-yint
      sy0(n, 1) = yint
    
      sz0(n, 0) = 1.0_num-zint
      sz0(n, 1) = zint

      ! Compute Ex on particle
      a = (sx0(n,0)*exg(j0,k-1,l-1) &
          + sx0(n,1)*exg(j0+1,k-1,l-1))*sy(n,-1)
      a = a + (sx0(n,0)*exg(j0,k,l-1) &
          + sx0(n,1)*exg(j0+1,k,l-1))*sy(n,0)
      a = a + (sx0(n,0)*exg(j0,k+1,l-1) &
          + sx0(n,1)*exg(j0+1,k+1,l-1))*sy(n,1)
      ex(nn) = ex(nn) + a*sz(n,-1)
      a = (sx0(n,0)*exg(j0,k-1,l) &
          + sx0(n,1)*exg(j0+1,k-1,l))*sy(n,-1)
      a = a + (sx0(n,0)*exg(j0,k,l) &
          + sx0(n,1)*exg(j0+1,k,l))*sy(n,0)
      a = a + (sx0(n,0)*exg(j0,k+1,l) &
          + sx0(n,1)*exg(j0+1,k+1,l))*sy(n,1)
      ex(nn) = ex(nn) + a*sz(n,0)
      a = (sx0(n,0)*exg(j0,k-1,l+1) &
          + sx0(n,1)*exg(j0+1,k-1,l+1))*sy(n,-1)
      a = a + (sx0(n,0)*exg(j0,k,l+1) &
          + sx0(n,1)*exg(j0+1,k,l+1))*sy(n,0)
      a = a + (sx0(n,0)*exg(j0,k+1,l+1) &
          + sx0(n,1)*exg(j0+1,k+1,l+1))*sy(n,1)
      ex(nn) = ex(nn) + a*sz(n,1)
    
      ! Compute Ey on particle
      a = (sx(n,-1)*eyg(j-1,k0,l-1) &
          + sx(n,0)*eyg(j,k0,l-1) &
          + sx(n,1)*eyg(j+1,k0,l-1))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l-1) &
          + sx(n,0)*eyg(j,k0+1,l-1) &
          + sx(n,1)*eyg(j+1,k0+1,l-1))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,-1)
      a = (sx(n,-1)*eyg(j-1,k0,l) &
          + sx(n,0)*eyg(j,k0,l) &
          + sx(n,1)*eyg(j+1,k0,l))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l) &
          + sx(n,0)*eyg(j,k0+1,l) &
          + sx(n,1)*eyg(j+1,k0+1,l))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,0)
      a = (sx(n,-1)*eyg(j-1,k0,l+1) &
          + sx(n,0)*eyg(j,k0,l+1) &
          + sx(n,1)*eyg(j+1,k0,l+1))*sy0(n,0)
      a = a + (sx(n,-1)*eyg(j-1,k0+1,l+1) &
          + sx(n,0)*eyg(j,k0+1,l+1) &
          + sx(n,1)*eyg(j+1,k0+1,l+1))*sy0(n,1)
      ey(nn) = ey(nn) + a*sz(n,1)
    
      ! Compute Ez on particle
      a = (sx(n,-1)*ezg(j-1,k-1,l0) &
          + sx(n,0)*ezg(j,k-1,l0) &
          + sx(n,1)*ezg(j+1,k-1,l0))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0) &
          + sx(n,0)*ezg(j,k,l0) &
          + sx(n,1)*ezg(j+1,k,l0))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0) &
          + sx(n,0)*ezg(j,k+1,l0) &
          + sx(n,1)*ezg(j+1,k+1,l0))*sy(n,1)
      ez(nn) = ez(nn) + a*sz0(n,0)
      a = (sx(n,-1)*ezg(j-1,k-1,l0+1) &
          + sx(n,0)*ezg(j,k-1,l0+1) &
          + sx(n,1)*ezg(j+1,k-1,l0+1))*sy(n,-1)
      a = a + (sx(n,-1)*ezg(j-1,k,l0+1) &
          + sx(n,0)*ezg(j,k,l0+1) &
          + sx(n,1)*ezg(j+1,k,l0+1))*sy(n,0)
      a = a + (sx(n,-1)*ezg(j-1,k+1,l0+1) &
          + sx(n,0)*ezg(j,k+1,l0+1) &
          + sx(n,1)*ezg(j+1,k+1,l0+1))*sy(n,1)
      ez(nn) = ez(nn) + a*sz0(n,1)
    
      ! Compute Bx on particle
      a = (sx(n,-1)*bxg(j-1,k0,l0) &
          + sx(n,0)*bxg(j,k0,l0) &
          + sx(n,1)*bxg(j+1,k0,l0))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0) &
          + sx(n,0)*bxg(j,k0+1,l0) &
          + sx(n,1)*bxg(j+1,k0+1,l0))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,0)
      a = (sx(n,-1)*bxg(j-1,k0,l0+1) &
          + sx(n,0)*bxg(j,k0,l0+1) &
          + sx(n,1)*bxg(j+1,k0,l0+1))*sy0(n,0)
      a = a + (sx(n,-1)*bxg(j-1,k0+1,l0+1) &
          + sx(n,0)*bxg(j,k0+1,l0+1) &
          + sx(n,1)*bxg(j+1,k0+1,l0+1))*sy0(n,1)
      bx(nn) = bx(nn) + a*sz0(n,1)
    
      ! Compute By on particle
      a = (sx0(n,0)*byg(j0,k-1,l0) &
          + sx0(n,1)*byg(j0+1,k-1,l0))*sy(n,-1)
      a = a + (sx0(n,0)*byg(j0,k,l0) &
          + sx0(n,1)*byg(j0+1,k,l0))*sy(n,0)
      a = a + (sx0(n,0)*byg(j0,k+1,l0) &
          + sx0(n,1)*byg(j0+1,k+1,l0))*sy(n,1)
      by(nn) = by(nn) + a*sz0(n,0)
      a = (sx0(n,0)*byg(j0,k-1,l0+1) &
          + sx0(n,1)*byg(j0+1,k-1,l0+1))*sy(n,-1)
      a = a + (sx0(n,0)*byg(j0,k,l0+1) &
          + sx0(n,1)*byg(j0+1,k,l0+1))*sy(n,0)
      a = a + (sx0(n,0)*byg(j0,k+1,l0+1) &
          + sx0(n,1)*byg(j0+1,k+1,l0+1))*sy(n,1)
      by(nn) = by(nn) + a*sz0(n,1)
    
      ! Compute Bz on particle
      a = (sx0(n,0)*bzg(j0,k0,l-1) &
          + sx0(n,1)*bzg(j0+1,k0,l-1))*sy0(n,0)
      a = a + (sx0(n,0)*bzg(j0,k0+1,l-1) &
          + sx0(n,1)*bzg(j0+1,k0+1,l-1))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,-1)
      a = (sx0(n,0)*bzg(j0,k0,l) &
          + sx0(n,1)*bzg(j0+1,k0,l))*sy0(n,0)
      a = a + (sx0(n,0)*bzg(j0,k0+1,l) &
          + sx0(n,1)*bzg(j0+1,k0+1,l))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,0)
      a = (sx0(n,0)*bzg(j0,k0,l+1) &
          + sx0(n,1)*bzg(j0+1,k0,l+1))*sy0(n,0)
      a = a + (sx0(n,0)*bzg(j0,k0+1,l+1) &
          + sx0(n,1)*bzg(j0+1,k0+1,l+1))*sy0(n,1)
      bz(nn) = bz(nn) + a*sz(n,1)
      
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

      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1

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
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy(n, 0) = 0.75_num-yintsq
        sy(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2
  
        xint=x-0.5_num-j0
        yint=y-0.5_num-k0
        zint=z-0.5_num-l0

        xintsq = xint*xint
        sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2
  
        yintsq = yint*yint
        sy0(n,-1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2
  
        zintsq = zint*zint
        sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,-1)*exg(j0-1,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,-1)*exg(j0,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,-1)*exg(j0+1,k-1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,-1)*exg(j0-1,k,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,-1)*exg(j0,k,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,-1)*exg(j0+1,k,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,-1)*exg(j0-1,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,-1)*exg(j0,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,-1)*exg(j0+1,k+1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,0)*exg(j0-1,k-1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,0)*exg(j0,k-1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,0)*exg(j0+1,k-1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,0)*exg(j0-1,k,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,0)*exg(j0,k,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,0)*exg(j0+1,k,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,0)*exg(j0-1,k+1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,0)*exg(j0,k+1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,0)*exg(j0+1,k+1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,-1)*sz(n,1)*exg(j0-1,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,-1)*sz(n,1)*exg(j0,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,-1)*sz(n,1)*exg(j0+1,k-1,l+1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,0)*sz(n,1)*exg(j0-1,k,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,0)*sz(n,1)*exg(j0,k,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,0)*sz(n,1)*exg(j0+1,k,l+1)
        ex(nn) = ex(nn) + sx0(n,-1)*sy(n,1)*sz(n,1)*exg(j0-1,k+1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sy(n,1)*sz(n,1)*exg(j0,k+1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sy(n,1)*sz(n,1)*exg(j0+1,k+1,l+1)
    
        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,-1)*eyg(j-1,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,-1)*eyg(j,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,-1)*eyg(j+1,k0-1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,-1)*eyg(j-1,k0,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,-1)*eyg(j,k0,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,-1)*eyg(j+1,k0,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,-1)*eyg(j-1,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,-1)*eyg(j,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,-1)*eyg(j+1,k0+1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,0)*eyg(j-1,k0-1,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,0)*eyg(j,k0-1,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,0)*eyg(j+1,k0-1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,0)*eyg(j-1,k0,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,0)*eyg(j,k0,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,0)*eyg(j+1,k0,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,0)*eyg(j-1,k0+1,l)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,0)*eyg(j,k0+1,l)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,0)*eyg(j+1,k0+1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,-1)*sz(n,1)*eyg(j-1,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,-1)*sz(n,1)*eyg(j,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,-1)*sz(n,1)*eyg(j+1,k0-1,l+1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,0)*sz(n,1)*eyg(j-1,k0,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,0)*sz(n,1)*eyg(j,k0,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,0)*sz(n,1)*eyg(j+1,k0,l+1)
        ey(nn) = ey(nn) + sx(n,-1)*sy0(n,1)*sz(n,1)*eyg(j-1,k0+1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sy0(n,1)*sz(n,1)*eyg(j,k0+1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sy0(n,1)*sz(n,1)*eyg(j+1,k0+1,l+1)
    
        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,-1)*ezg(j-1,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,-1)*ezg(j,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,-1)*ezg(j+1,k-1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,-1)*ezg(j-1,k,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,-1)*ezg(j,k,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,-1)*ezg(j+1,k,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,-1)*ezg(j-1,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,-1)*ezg(j,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,-1)*ezg(j+1,k+1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,0)*ezg(j-1,k-1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,0)*ezg(j,k-1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,0)*ezg(j+1,k-1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,0)*ezg(j-1,k,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,0)*ezg(j,k,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,0)*ezg(j+1,k,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,0)*ezg(j-1,k+1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,0)*ezg(j,k+1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,0)*ezg(j+1,k+1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,-1)*sz0(n,1)*ezg(j-1,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,-1)*sz0(n,1)*ezg(j,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,-1)*sz0(n,1)*ezg(j+1,k-1,l0+1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,0)*sz0(n,1)*ezg(j-1,k,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,0)*sz0(n,1)*ezg(j,k,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,0)*sz0(n,1)*ezg(j+1,k,l0+1)
        ez(nn) = ez(nn) + sx(n,-1)*sy(n,1)*sz0(n,1)*ezg(j-1,k+1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sy(n,1)*sz0(n,1)*ezg(j,k+1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sy(n,1)*sz0(n,1)*ezg(j+1,k+1,l0+1)

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,-1)*bxg(j-1,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,-1)*bxg(j,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,-1)*bxg(j+1,k0-1,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,-1)*bxg(j-1,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,-1)*bxg(j,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,-1)*bxg(j+1,k0,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,-1)*bxg(j-1,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,-1)*bxg(j,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,-1)*bxg(j+1,k0+1,l0-1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,0)*bxg(j-1,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,0)*bxg(j,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,0)*bxg(j+1,k0-1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,0)*bxg(j-1,k0,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,0)*bxg(j,k0,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,0)*bxg(j+1,k0,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,0)*bxg(j-1,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,0)*bxg(j,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,0)*bxg(j+1,k0+1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,-1)*sz0(n,1)*bxg(j-1,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,-1)*sz0(n,1)*bxg(j,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,-1)*sz0(n,1)*bxg(j+1,k0-1,l0+1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,0)*sz0(n,1)*bxg(j-1,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,0)*sz0(n,1)*bxg(j,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,0)*sz0(n,1)*bxg(j+1,k0,l0+1)
        bx(nn) = bx(nn) + sx(n,-1)*sy0(n,1)*sz0(n,1)*bxg(j-1,k0+1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sy0(n,1)*sz0(n,1)*bxg(j,k0+1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sy0(n,1)*sz0(n,1)*bxg(j+1,k0+1,l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,-1)*byg(j0-1,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,-1)*byg(j0,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,-1)*byg(j0+1,k-1,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,-1)*byg(j0-1,k,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,-1)*byg(j0,k,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,-1)*byg(j0+1,k,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,-1)*byg(j0-1,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,-1)*byg(j0,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,-1)*byg(j0+1,k+1,l0-1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,0)*byg(j0-1,k-1,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,0)*byg(j0,k-1,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,0)*byg(j0+1,k-1,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,0)*byg(j0-1,k,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,0)*byg(j0,k,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,0)*byg(j0+1,k,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,0)*byg(j0-1,k+1,l0)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,0)*byg(j0,k+1,l0)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,0)*byg(j0+1,k+1,l0)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,-1)*sz0(n,1)*byg(j0-1,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,-1)*sz0(n,1)*byg(j0,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,-1)*sz0(n,1)*byg(j0+1,k-1,l0+1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,0)*sz0(n,1)*byg(j0-1,k,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,0)*sz0(n,1)*byg(j0,k,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,0)*sz0(n,1)*byg(j0+1,k,l0+1)
        by(nn) = by(nn) + sx0(n,-1)*sy(n,1)*sz0(n,1)*byg(j0-1,k+1,l0+1)
        by(nn) = by(nn) + sx0(n,0)*sy(n,1)*sz0(n,1)*byg(j0,k+1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sy(n,1)*sz0(n,1)*byg(j0+1,k+1,l0+1)
    
        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,-1)*bzg(j0-1,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,-1)*bzg(j0,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,-1)*bzg(j0+1,k0-1,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,-1)*bzg(j0-1,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,-1)*bzg(j0,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,-1)*bzg(j0+1,k0,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,-1)*bzg(j0-1,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,-1)*bzg(j0,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,-1)*bzg(j0+1,k0+1,l-1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,0)*bzg(j0-1,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,0)*bzg(j0,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,0)*bzg(j0+1,k0-1,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,0)*bzg(j0-1,k0,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,0)*bzg(j0,k0,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,0)*bzg(j0+1,k0,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,0)*bzg(j0-1,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,0)*bzg(j0,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,0)*bzg(j0+1,k0+1,l)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,-1)*sz(n,1)*bzg(j0-1,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,-1)*sz(n,1)*bzg(j0,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,-1)*sz(n,1)*bzg(j0+1,k0-1,l+1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,0)*sz(n,1)*bzg(j0-1,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,0)*sz(n,1)*bzg(j0,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,0)*sz(n,1)*bzg(j0+1,k0,l+1)
        bz(nn) = bz(nn) + sx0(n,-1)*sy0(n,1)*sz(n,1)*bzg(j0-1,k0+1,l+1)
        bz(nn) = bz(nn) + sx0(n,0)*sy0(n,1)*sz(n,1)*bzg(j0,k0+1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sy0(n,1)*sz(n,1)*bzg(j0+1,k0+1,l+1)

      ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD 
#endif  
    ENDDO  
  
  ENDIF
  
  RETURN     
END SUBROUTINE
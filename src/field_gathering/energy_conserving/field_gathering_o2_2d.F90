! ______________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! FIELD_GATHERING_O2_3D.F90
!
! Field gathering subroutines in 2D at order 2.
!
! List of subroutines:
!
! - pxr_gete2dxz_energy_conserving_vect_2_2
! - pxr_getb2dxz_energy_conserving_vect_2_2
!
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Field gathering cartesian in 2D for the electric field at order 2
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 12/01/2016
!
!> @param[in] np Number of particles
!> @param[in] xp,zp particle position arrays
!> @param[inout] ex,ey,ez electric field particle arrays
!> @param[in] xmin,zmin tile boundaries
!> @param[in] dx,dz space steps
!> @param[in] nx,nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] exg, eyg,ezg field arrays
!> @param[in] lvect vector size for the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
subroutine pxr_gete2dxz_energy_conserving_vect_2_2(     &
    np,xp,zp,ex,ey,ez,xmin,zmin,dx,dz,                  &
    exg,exg_nguard,exg_nvalid,                          &
    eyg,eyg_nguard,eyg_nvalid,                          &
    ezg,ezg_nguard,ezg_nvalid,                          &
    lvect,l_lower_order_in_v)  !#do not wrap
! ______________________________________________________________________________
  use constants
  implicit none

  integer(idp)                  :: np
  integer(idp), intent(IN)      :: exg_nguard(2),exg_nvalid(2),&
                                   eyg_nguard(2),eyg_nvalid(2),&
                                   ezg_nguard(2),ezg_nvalid(2)
  integer(idp)                  :: lvect
  real(num), dimension(np)      :: xp,zp,ex,ey,ez
  logical(idp)                  :: l_lower_order_in_v
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,1, &
                              -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,1, &
                              -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,1, &
                              -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  real(num)                     :: xmin,zmin,dx,dz
  integer(isp)                  :: ip, j, l, j0, l0
  integer(isp)                  :: nn,n
  real(num)                     :: dxi, dzi, x, z, xint, zint
  real(num)                     :: xintsq,zintsq
  real(num), DIMENSION(lvect,-1:1)    :: sx,sx0
  real(num), DIMENSION(lvect,-1:1)    :: sz,sz0
  real(num), parameter          :: onesixth=1./6.,twothird=2./3.

  dxi = 1./dx
  dzi = 1./dz

  ! ___________________________
  IF (l_lower_order_in_v) THEN

    ! Loop over the particles by block
    DO ip=1,np,lvect

#if defined __INTEL_COMPILER
!!DIR$ IVDEP
!!DIR$ DISTRIBUTE POINT
!DIR$ ASSUME_ALIGNED xp:64,zp:64
!DIR$ ASSUME_ALIGNED sx:64,sz:64
!DIR$ ASSUME_ALIGNED sx0:64,sz0:64
!DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=nint(x)
        j0=floor(x-0.5_num)
        l=nint(z)
        l0=floor(z-0.5_num)
        xint=x-j
        zint=z-l

        ! Compute shape factors
        xintsq = xint*xint
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint

        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n,0)*sz(n,-1)*exg(j0,1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sz(n,-1)*exg(j0+1,1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sz(n,0)*exg(j0,1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sz(n,0)*exg(j0+1,1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sz(n,1)*exg(j0,1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sz(n,1)*exg(j0+1,1,l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n,-1)*sz(n,-1)*eyg(j-1,1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sz(n,-1)*eyg(j,1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sz(n,-1)*eyg(j+1,1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sz(n,0)*eyg(j-1,1,l)
        ey(nn) = ey(nn) + sx(n,0)*sz(n,0)*eyg(j,1,l)
        ey(nn) = ey(nn) + sx(n,1)*sz(n,0)*eyg(j+1,1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sz(n,1)*eyg(j-1,1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sz(n,1)*eyg(j,1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sz(n,1)*eyg(j+1,1,l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n,-1)*sz0(n,0)*ezg(j-1,1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sz0(n,0)*ezg(j,1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sz0(n,0)*ezg(j+1,1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sz0(n,1)*ezg(j-1,1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sz0(n,1)*ezg(j,1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sz0(n,1)*ezg(j+1,1,l0+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    ENDDO

  ! ___________________________
  ! l_lower_order_in_v is false
  ELSE

    ! Loop over the particles by block
    DO ip=1,np,lvect

#if defined __INTEL_COMPILER
!!DIR$ IVDEP
!!DIR$ DISTRIBUTE POINT
!DIR$ ASSUME_ALIGNED xp:64,zp:64
!DIR$ ASSUME_ALIGNED sx:64,sz:64
!DIR$ ASSUME_ALIGNED sx0:64,sz0:64
!DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=nint(x)
        j0=floor(x)
        l=nint(z)
        l0=floor(z)
        xint=x-j
        zint=z-l

        ! Compute shape factors
        xintsq = xint*xint
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

        xintsq = xint*xint
        sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        zintsq = zint*zint
        sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n,-1)*sz(n,-1)*exg(j0-1,1,l-1)
        ex(nn) = ex(nn) + sx0(n,0)*sz(n,-1)*exg(j0,1,l-1)
        ex(nn) = ex(nn) + sx0(n,1)*sz(n,-1)*exg(j0+1,1,l-1)
        ex(nn) = ex(nn) + sx0(n,-1)*sz(n,0)*exg(j0-1,1,l)
        ex(nn) = ex(nn) + sx0(n,0)*sz(n,0)*exg(j0,1,l)
        ex(nn) = ex(nn) + sx0(n,1)*sz(n,0)*exg(j0+1,1,l)
        ex(nn) = ex(nn) + sx0(n,-1)*sz(n,1)*exg(j0-1,1,l+1)
        ex(nn) = ex(nn) + sx0(n,0)*sz(n,1)*exg(j0,1,l+1)
        ex(nn) = ex(nn) + sx0(n,1)*sz(n,1)*exg(j0+1,1,l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n,-1)*sz(n,-1)*eyg(j-1,1,l-1)
        ey(nn) = ey(nn) + sx(n,0)*sz(n,-1)*eyg(j,1,l-1)
        ey(nn) = ey(nn) + sx(n,1)*sz(n,-1)*eyg(j+1,1,l-1)
        ey(nn) = ey(nn) + sx(n,-1)*sz(n,0)*eyg(j-1,1,l)
        ey(nn) = ey(nn) + sx(n,0)*sz(n,0)*eyg(j,1,l)
        ey(nn) = ey(nn) + sx(n,1)*sz(n,0)*eyg(j+1,1,l)
        ey(nn) = ey(nn) + sx(n,-1)*sz(n,1)*eyg(j-1,1,l+1)
        ey(nn) = ey(nn) + sx(n,0)*sz(n,1)*eyg(j,1,l+1)
        ey(nn) = ey(nn) + sx(n,1)*sz(n,1)*eyg(j+1,1,l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n,-1)*sz0(n,-1)*ezg(j-1,1,l0-1)
        ez(nn) = ez(nn) + sx(n,0)*sz0(n,-1)*ezg(j,1,l0-1)
        ez(nn) = ez(nn) + sx(n,1)*sz0(n,-1)*ezg(j+1,1,l0-1)
        ez(nn) = ez(nn) + sx(n,-1)*sz0(n,0)*ezg(j-1,1,l0)
        ez(nn) = ez(nn) + sx(n,0)*sz0(n,0)*ezg(j,1,l0)
        ez(nn) = ez(nn) + sx(n,1)*sz0(n,0)*ezg(j+1,1,l0)
        ez(nn) = ez(nn) + sx(n,-1)*sz0(n,1)*ezg(j-1,1,l0+1)
        ez(nn) = ez(nn) + sx(n,0)*sz0(n,1)*ezg(j,1,l0+1)
        ez(nn) = ez(nn) + sx(n,1)*sz0(n,1)*ezg(j+1,1,l0+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    ENDDO

  ENDIF

end subroutine


! ______________________________________________________________________________
!> @brief
!> Field gathering cartesian in 2D for the magnetic field at order 2
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 12/01/2016
!
!> @param[in] np Number of particles
!> @param[in] xp,zp particle position arrays
!> @param[inout] bx,by,bz magnetic field particle arrays
!> @param[in] xmin,zmin tile boundaries
!> @param[in] dx,dz space steps
!> @param[in] nx,nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] bxg, byg,bzg field arrays
!> @param[in] lvect the vector length of the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
subroutine pxr_getb2dxz_energy_conserving_vect_2_2(     &
    np,xp,zp,bx,by,bz,xmin,zmin,dx,dz,                  &
    bxg,bxg_nguard,bxg_nvalid,                          &
    byg,byg_nguard,byg_nvalid,                          &
    bzg,bzg_nguard,bzg_nvalid,                          &
    lvect,l_lower_order_in_v) !#do not wrap
! ______________________________________________________________________________

  use constants
  implicit none

  ! __ Parameter declaration ___________________________________________
  integer(idp)                       :: np
  integer(idp), intent(IN)                :: bxg_nguard(2),bxg_nvalid(2),&
                                             byg_nguard(2),byg_nvalid(2),&
                                             bzg_nguard(2),bzg_nvalid(2)
  integer(idp)                       :: lvect
  real(num), dimension(np)           :: xp,zp,bx,by,bz
  logical(idp)                       :: l_lower_order_in_v
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,1, &
                              -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,1, &
                              -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,1, &
                              -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)
  real(num)                          :: xmin,zmin,dx,dz
  integer(idp)                       :: ip, j, l
  integer(idp)                       :: j0, l0
  integer(idp)                       :: n,nn
  real(num)                          :: dxi, dzi, x, z, xint, zint
  real(num)                          :: xintsq,zintsq
  real(num), DIMENSION(lvect,-1:1)   :: sx, sx0
  real(num), DIMENSION(lvect,-1:1)   :: sz, sz0
  real(num), parameter               :: onesixth=1./6.,twothird=2./3.

  ! ___________________________
  ! Compute parameters

  dxi = 1./dx
  dzi = 1./dz

  sx=0
  sz=0.
  sx0=0.
  sz0=0.

  ! ___________________________
  IF (l_lower_order_in_v) THEN

    ! Loop over the particles by block
    DO ip=1,np,lvect

#if defined __INTEL_COMPILER
!!DIR$ IVDEP
!!DIR$ DISTRIBUTE POINT
!DIR$ ASSUME_ALIGNED xp:64,zp:64
!DIR$ ASSUME_ALIGNED sx:64,sz:64
!DIR$ ASSUME_ALIGNED sx0:64,sz0:64
!DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        j=nint(x)
        j0=floor(x-0.5)

        l=nint(z)
        l0=floor(z-0.5)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        xintsq = xint*xint
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2
        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint
        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n,-1)*sz0(n,0)*bxg(j-1,1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sz0(n,0)*bxg(j,1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sz0(n,0)*bxg(j+1,1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sz0(n,1)*bxg(j-1,1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sz0(n,1)*bxg(j,1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sz0(n,1)*bxg(j+1,1,l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n,0)*sz0(n,0)*byg(j0,1,l0)
        by(nn) = by(nn) + sx0(n,1)*sz0(n,0)*byg(j0+1,1,l0)
        by(nn) = by(nn) + sx0(n,0)*sz0(n,1)*byg(j0,1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sz0(n,1)*byg(j0+1,1,l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n,0)*sz(n,-1)*bzg(j0,1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sz(n,-1)*bzg(j0+1,1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sz(n,0)*bzg(j0,1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sz(n,0)*bzg(j0+1,1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sz(n,1)*bzg(j0,1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sz(n,1)*bzg(j0+1,1,l+1)

      end do
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    end do

  ! ___________________________
  ! l_lower_order_in_v is false
  else

    ! Loop over the particles by block
    DO ip=1,np,lvect

#if defined __INTEL_COMPILER
!!DIR$ IVDEP
!!DIR$ DISTRIBUTE POINT
!DIR$ ASSUME_ALIGNED xp:64,zp:64
!DIR$ ASSUME_ALIGNED sx:64,sz:64
!DIR$ ASSUME_ALIGNED sx0:64,sz0:64
!DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
      ! Loop over the particles inside a block
      DO n=1,MIN(lvect,np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        j=nint(x)
        j0=floor(x)

        l=nint(z)
        l0=floor(z)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        xintsq = xint*xint
        sx(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx(n, 0) = 0.75_num-xintsq
        sx(n, 1) = 0.5_num*(0.5_num+xint)**2

        zintsq = zint*zint
        sz(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz(n, 0) = 0.75_num-zintsq
        sz(n, 1) = 0.5_num*(0.5_num+zint)**2

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

        xintsq = xint*xint
        sx0(n,-1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        zintsq = zint*zint
        sz0(n,-1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n,-1)*sz0(n,0)*bxg(j-1,1,l0)
        bx(nn) = bx(nn) + sx(n,0)*sz0(n,0)*bxg(j,1,l0)
        bx(nn) = bx(nn) + sx(n,1)*sz0(n,0)*bxg(j+1,1,l0)
        bx(nn) = bx(nn) + sx(n,-1)*sz0(n,1)*bxg(j-1,1,l0+1)
        bx(nn) = bx(nn) + sx(n,0)*sz0(n,1)*bxg(j,1,l0+1)
        bx(nn) = bx(nn) + sx(n,1)*sz0(n,1)*bxg(j+1,1,l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n,0)*sz0(n,0)*byg(j0,1,l0)
        by(nn) = by(nn) + sx0(n,1)*sz0(n,0)*byg(j0+1,1,l0)
        by(nn) = by(nn) + sx0(n,0)*sz0(n,1)*byg(j0,1,l0+1)
        by(nn) = by(nn) + sx0(n,1)*sz0(n,1)*byg(j0+1,1,l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n,0)*sz(n,-1)*bzg(j0,1,l-1)
        bz(nn) = bz(nn) + sx0(n,1)*sz(n,-1)*bzg(j0+1,1,l-1)
        bz(nn) = bz(nn) + sx0(n,0)*sz(n,0)*bzg(j0,1,l)
        bz(nn) = bz(nn) + sx0(n,1)*sz(n,0)*bzg(j0+1,1,l)
        bz(nn) = bz(nn) + sx0(n,0)*sz(n,1)*bzg(j0,1,l+1)
        bz(nn) = bz(nn) + sx0(n,1)*sz(n,1)*bzg(j0+1,1,l+1)

      enddo
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    end do

  ENDIF
  return

end subroutine

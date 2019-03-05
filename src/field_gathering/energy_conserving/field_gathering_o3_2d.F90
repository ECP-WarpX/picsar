! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software, ! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! FIELD_GATHERING_O3_2D.F90
!
! Field gathering subroutines in 2D at order 3.
!
! List of subroutines:
!
! - pxr_gete2dxz_energy_conserving_scalar_3_3
! - pxr_getb2dxz_energy_conserving_scalar_3_3
! - pxr_gete2dxz_energy_conserving_vect_3_3
! - pxr_getb2dxz_energy_conserving_vect_3_3
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Scalar subroutine for the cartesian electric field gathering in 2D
!
!> @details
!> This function is not vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 12/01/2016
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!> exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_gete2dxz_energy_conserving_scalar_3_3( np, xp, zp, ex, ey, ez, xmin,   &
  zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,          &
  ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp)             :: np
  integer(idp), intent(IN) :: exg_nguard(2), exg_nvalid(2), eyg_nguard(2),            &
  eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  real(num), dimension(np) :: xp, zp, ex, ey, ez
  logical(lp)              :: l_lower_order_in_v, l_nodal
  real(num)                :: stagger_shift
  real(num)                :: xmin, zmin, dx, dz
  integer(idp)             :: ip, j, l, ixmin, ixmax, izmin, izmax
  integer(idp)             :: ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0, l0
  real(num)                :: dxi, dzi, x, z, xint, zint
  real(num)                :: xintsq, oxint, zintsq, ozint, oxintsq, ozintsq
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1, 1,        &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1, 1,        &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, 1,        &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  REAL(num), DIMENSION(-1:2)            :: sx, sx0
  REAL(num), DIMENSION(-1:2)            :: sz, sz0
  real(num), parameter                  :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz

  ixmin = -1
  ixmax =  1
  izmin = -1
  izmax =  1

  if (l_lower_order_in_v) then
    ixmin0 = -1
    ixmax0 =  1
    izmin0 = -1
    izmax0 =  1

  else
    ixmin0 = -1
    ixmax0 =  2
    izmin0 = -1
    izmax0 =  2
  end if

  if (l_lower_order_in_v) then

    !$acc parallel deviceptr(exg, eyg, ezg, xp, zp, ex, ey, ez)
    !$acc loop gang vector private(sx(-1:2), sz(-1:2), sx0(-1:2), sz0(-1:2))
    do ip=1, np

      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      j=floor(x)
      j0=floor(x+0.5_num-stagger_shift)

      l=floor(z)
      l0=floor(z+0.5_num-stagger_shift)

      xint=x-j
      zint=z-l

      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.-xint/2)
      sx( 1) = twothird-oxintsq*(1.-oxint/2)
      sx( 2) = onesixth*xintsq*xint

      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      xintsq = xint*xint
      sx0(-1) = 0.5*(0.5-xint)**2
      sx0( 0) = 0.75-xintsq
      sx0( 1) = 0.5*(0.5+xint)**2

      zintsq = zint*zint
      sz0(-1) = 0.5*(0.5-zint)**2
      sz0( 0) = 0.75-zintsq
      sz0( 1) = 0.5*(0.5+zint)**2

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j0+jj, 1, l+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin, ixmax+1
          ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj, 1, l+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
        do jj = ixmin, ixmax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj, 1, l0+ll)
        end do
      end do
      !$acc end loop

    enddo
    !$acc end loop
    !$acc end parallel

 else

    !$acc parallel deviceptr(exg, eyg, ezg, xp, zp, ex, ey, ez)
    !$acc loop gang vector private(sx(-1:2), sz(-1:2), sx0(-1:2), sz0(-1:2))
    do ip=1, np

      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      j=floor(x)
      j0=floor(x-stagger_shift)

      l=floor(z)
      l0=floor(z-stagger_shift)

      xint=x-j
      zint=z-l

      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.-xint/2)
      sx( 1) = twothird-oxintsq*(1.-oxint/2)
      sx( 2) = onesixth*xintsq*xint

      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.-xint/2)
      sx0( 1) = twothird-oxintsq*(1.-oxint/2)
      sx0( 2) = onesixth*xintsq*xint

      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.-zint/2)
      sz0( 1) = twothird-ozintsq*(1.-ozint/2)
      sz0( 2) = onesixth*zintsq*zint

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j0+jj, 1, l+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin, ixmax+1

          ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj, 1, l+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin, ixmax+1
          ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj, 1, l0+ll)
        end do
      end do
      !$acc end loop
    enddo
    !$acc end loop
    !$acc end parallel
  endif

end subroutine

! ________________________________________________________________________________________
!> @brief
!> Field gathering cartesian in 2D for the electric field
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
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!>  exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] lvect vector size for the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_gete2dxz_energy_conserving_vect_3_3( np, xp, zp, ex, ey, ez, xmin,     &
  zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,          &
  ezg_nguard, ezg_nvalid, lvect, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  ! ______________________________________________________________________________
  implicit none
  integer(idp)                  :: np
  integer(idp), intent(IN)      :: exg_nguard(2), exg_nvalid(2), eyg_nguard(2),       &
  eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  integer(idp)                  :: lvect
  real(num), dimension(np)      :: xp, zp, ex, ey, ez
  logical(lp)                   :: l_lower_order_in_v, l_nodal
  real(num)                     :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1, 1,        &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1, 1,        &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, 1,        &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  real(num)                     :: xmin, zmin, dx, dz
  integer(isp)                  :: ip, j, l, j0, l0
  integer(isp)                  :: nn, n
  real(num)                     :: dxi, dzi, x, z, xint, zint
  real(num)                     :: xintsq, oxint, zintsq, ozint, oxintsq, ozintsq
  real(num), DIMENSION(lvect, -1:2)    :: sx, sx0
  real(num), DIMENSION(lvect, -1:2)    :: sz, sz0
  real(num), parameter          :: onesixth=1./6., twothird=2./3.

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dzi = 1./dz

  ! ___________________________
  IF (l_lower_order_in_v) THEN

    ! Loop over the particles by block
    DO ip=1, np, lvect

#if defined __INTEL_COMPILER
      !!DIR$ IVDEP
      !!DIR$ DISTRIBUTE POINT
      !DIR$ ASSUME_ALIGNED xp:64, zp:64
      !DIR$ ASSUME_ALIGNED sx:64, sz:64
      !DIR$ ASSUME_ALIGNED sx0:64, sz0:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
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
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)

        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, -1)*exg(j0-1, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, -1)*exg(j0, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, -1)*exg(j0+1, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, 0)*exg(j0-1, 1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 0)*exg(j0, 1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 0)*exg(j0+1, 1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, 1)*exg(j0-1, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 1)*exg(j0, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 1)*exg(j0+1, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, 2)*exg(j0-1, 1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 2)*exg(j0, 1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 2)*exg(j0+1, 1, l+2)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, -1)*eyg(j-1, 1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, -1)*eyg(j, 1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, -1)*eyg(j+1, 1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, -1)*eyg(j+2, 1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, 0)*eyg(j-1, 1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 0)*eyg(j, 1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 0)*eyg(j+1, 1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, 0)*eyg(j+2, 1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, 1)*eyg(j-1, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 1)*eyg(j, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 1)*eyg(j+1, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, 1)*eyg(j+2, 1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, 2)*eyg(j-1, 1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 2)*eyg(j, 1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 2)*eyg(j+1, 1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, 2)*eyg(j+2, 1, l+2)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, -1)*ezg(j-1, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, -1)*ezg(j, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, -1)*ezg(j+1, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, -1)*ezg(j+2, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, 0)*ezg(j-1, 1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 0)*ezg(j, 1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 0)*ezg(j+1, 1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, 0)*ezg(j+2, 1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, 1)*ezg(j-1, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 1)*ezg(j, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 1)*ezg(j+1, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, 1)*ezg(j+2, 1, l0+1)
      END DO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    ENDDO
    ! ___________________________
    ! l_lower_order_in_v is false
  ELSE

    DO ip=1, np

#if defined __INTEL_COMPILER
      !!DIR$ IVDEP
      !!DIR$ DISTRIBUTE POINT
      !DIR$ ASSUME_ALIGNED xp:64, zp:64
      !DIR$ ASSUME_ALIGNED sx:64, sz:64
      !DIR$ ASSUME_ALIGNED sx0:64, sz0:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
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
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)
        l=floor(z)
        l0=floor(z-stagger_shift)

        xint=x-j
        zint=z-l

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, -1)*exg(j0-1, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, -1)*exg(j0, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, -1)*exg(j0+1, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sz(n, -1)*exg(j0+2, 1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, 0)*exg(j0-1, 1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 0)*exg(j0, 1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 0)*exg(j0+1, 1, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sz(n, 0)*exg(j0+2, 1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, 1)*exg(j0-1, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 1)*exg(j0, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 1)*exg(j0+1, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sz(n, 1)*exg(j0+2, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sz(n, 2)*exg(j0-1, 1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 2)*exg(j0, 1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 2)*exg(j0+1, 1, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sz(n, 2)*exg(j0+2, 1, l+2)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, -1)*eyg(j-1, 1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, -1)*eyg(j, 1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, -1)*eyg(j+1, 1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, -1)*eyg(j+2, 1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, 0)*eyg(j-1, 1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 0)*eyg(j, 1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 0)*eyg(j+1, 1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, 0)*eyg(j+2, 1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, 1)*eyg(j-1, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 1)*eyg(j, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 1)*eyg(j+1, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, 1)*eyg(j+2, 1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sz(n, 2)*eyg(j-1, 1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 2)*eyg(j, 1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 2)*eyg(j+1, 1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sz(n, 2)*eyg(j+2, 1, l+2)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, -1)*ezg(j-1, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, -1)*ezg(j, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, -1)*ezg(j+1, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, -1)*ezg(j+2, 1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, 0)*ezg(j-1, 1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 0)*ezg(j, 1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 0)*ezg(j+1, 1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, 0)*ezg(j+2, 1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, 1)*ezg(j-1, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 1)*ezg(j, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 1)*ezg(j+1, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, 1)*ezg(j+2, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sz0(n, 2)*ezg(j-1, 1, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 2)*ezg(j, 1, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 2)*ezg(j+1, 1, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sz0(n, 2)*ezg(j+2, 1, l0+2)

      end do
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    ENDDO

  ENDIF
  return
end subroutine pxr_gete2dxz_energy_conserving_vect_3_3


! ________________________________________________________________________________________
!> @brief
!> Scalar Cartesian subroutine for the magnetic field gathering in 2D at order 3.
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
!> @param[in] xp, zp particle position arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the
!> bxg, byg, bzg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 2 integers)
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_getb2dxz_energy_conserving_scalar_3_3( np, xp, zp, bx, by, bz, xmin,   &
  zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,          &
  bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  ! __ Parameter declaration ___________________________________________
  integer(idp), intent(in)                :: np
  integer(idp), intent(IN)                :: bxg_nguard(2), bxg_nvalid(2),            &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  real(num), dimension(np), intent(in)    :: xp, zp
  real(num), dimension(np), intent(inout) :: bx, by, bz
  logical(lp) , intent(in)                :: l_lower_order_in_v, l_nodal
  real(num)          :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1, 1,        &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1, 1,        &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1, 1,        &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)

  real(num)                               :: xmin, zmin, dx, dz
  integer(idp)             :: ip, j, l, ixmin, ixmax, izmin, izmax
  integer(idp)             :: ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0, l0
  real(num)                :: dxi, dzi, x, z, xint, zint
  real(num)                :: xintsq, oxint, zintsq, ozint, oxintsq, ozintsq
  REAL(num), DIMENSION(-1:2)            :: sx, sx0
  REAL(num), DIMENSION(-1:2)            :: sz, sz0
  real(num), parameter                  :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz

  ixmin = -1
  ixmax =  1
  izmin = -1
  izmax =  1

  if (l_lower_order_in_v) then
    ixmin0 = -1
    ixmax0 =  1
    izmin0 = -1
    izmax0 =  1
  else
    ixmin0 = -1
    ixmax0 =  2
    izmin0 = -1
    izmax0 =  2
  end if

  if (l_lower_order_in_v) then

    !$acc parallel deviceptr(bxg, byg, bzg, xp, zp, bx, by, bz)
    !$acc loop gang vector private(sx(-1:2), sz(-1:2), sx0(-1:2), sz0(-1:2))
    do ip=1, np
      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      j=floor(x)
      j0=floor(x+0.5_num-stagger_shift)

      l=floor(z)
      l0=floor(z+0.5_num-stagger_shift)

      xint=x-j
      zint=z-l

      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.-xint/2)
      sx( 1) = twothird-oxintsq*(1.-oxint/2)
      sx( 2) = onesixth*xintsq*xint

      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      xintsq = xint*xint
      sx0(-1) = 0.5*(0.5-xint)**2
      sx0( 0) = 0.75-xintsq
      sx0( 1) = 0.5*(0.5+xint)**2

      zintsq = zint*zint
      sz0(-1) = 0.5*(0.5-zint)**2
      sz0( 0) = 0.75-zintsq
      sz0( 1) = 0.5*(0.5+zint)**2

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj, 1, l0+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin0, ixmax0
          by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj, 1, l0+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin0, ixmax0
          bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j0+jj, 1, l+ll)
        end do
      end do
      !$acc end loop

    enddo
    !$acc end loop
    !$acc end parallel

  else

    !$acc parallel deviceptr(bxg, byg, bzg, xp, zp, bx, by, bz)
    !$acc loop gang vector private(sx(-1:2), sz(-1:2), sx0(-1:2), sz0(-1:2))
    do ip=1, np
      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      j=floor(x)
      j0=floor(x-stagger_shift)

      l=floor(z)
      l0=floor(z-stagger_shift)

      xint=x-j
      zint=z-l

      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.-xint/2)
      sx( 1) = twothird-oxintsq*(1.-oxint/2)
      sx( 2) = onesixth*xintsq*xint

      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.-xint/2)
      sx0( 1) = twothird-oxintsq*(1.-oxint/2)
      sx0( 2) = onesixth*xintsq*xint

      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.-zint/2)
      sz0( 1) = twothird-ozintsq*(1.-ozint/2)
      sz0( 2) = onesixth*zintsq*zint

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj, 1, l0+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin0, ixmax0
          by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj, 1, l0+ll)
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
        ! Prevent wrong vectorization from the compiler
        !DIR$ NOVECTOR
        do jj = ixmin0, ixmax0
          bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j0+jj, 1, l+ll)
        end do
      end do
      !$acc end loop

    enddo
    !$acc end loop
    !$acc end parallel

  end if
end subroutine

! ________________________________________________________________________________________
!> @brief
!> Field gathering cartesian in 2D for the magnetic field at order 3
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 12/01/2016
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the
!> bxg, byg, bzg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 2 integers)
!> @param[in] lvect the vector length of the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_getb2dxz_energy_conserving_vect_3_3( np, xp, zp, bx, by, bz, xmin,     &
  zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,          &
  bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  ! ______________________________________________________________________________

  implicit none

  ! __ Parameter declaration ___________________________________________
  integer(idp)                       :: np
  integer(idp), intent(IN)                :: bxg_nguard(2), bxg_nvalid(2),            &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  integer(idp)                       :: lvect
  real(num), dimension(np)           :: xp, zp, bx, by, bz
  logical(lp)                        :: l_lower_order_in_v, l_nodal
  real(num)                          :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1, 1,        &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1, 1,        &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1, 1,        &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)
  real(num)                          :: xmin, zmin, dx, dz
  integer(idp)                       :: ip, j, l
  integer(idp)                       :: j0, l0
  integer(idp)                       :: n, nn
  real(num)                          :: dxi, dzi, x, z, xint, zint
  real(num)                          :: xintsq, oxint, zintsq, ozint, oxintsq,        &
  ozintsq
  real(num), DIMENSION(lvect, -1:2)   :: sx, sx0
  real(num), DIMENSION(lvect, -1:2)   :: sz, sz0
  real(num), parameter               :: onesixth=1./6., twothird=2./3.

  ! ___________________________
  ! Compute parameters

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dzi = 1./dz

  sx=0
  sz=0.
  sx0=0.
  sz0=0.

  ! ___________________________
  IF (l_lower_order_in_v) THEN

    ! Loop over the particles by block
    DO ip=1, np, lvect

#if defined __INTEL_COMPILER
      !!DIR$ IVDEP
      !!DIR$ DISTRIBUTE POINT
      !DIR$ ASSUME_ALIGNED xp:64, zp:64
      !DIR$ ASSUME_ALIGNED sx:64, sz:64
      !DIR$ ASSUME_ALIGNED sx0:64, sz0:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
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
      DO n=1, MIN(lvect, np-ip+1)


        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)

        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, -1)*bxg(j-1, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, -1)*bxg(j, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, -1)*bxg(j+1, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, -1)*bxg(j+2, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, 0)*bxg(j-1, 1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 0)*bxg(j, 1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 0)*bxg(j+1, 1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, 0)*bxg(j+2, 1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, 1)*bxg(j-1, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 1)*bxg(j, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 1)*bxg(j+1, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, 1)*bxg(j+2, 1, l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, -1)*byg(j0-1, 1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, -1)*byg(j0, 1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, -1)*byg(j0+1, 1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, 0)*byg(j0-1, 1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 0)*byg(j0, 1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, 0)*byg(j0+1, 1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, 1)*byg(j0-1, 1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 1)*byg(j0, 1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, 1)*byg(j0+1, 1, l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, -1)*bzg(j0-1, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, -1)*bzg(j0, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, -1)*bzg(j0+1, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, 0)*bzg(j0-1, 1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 0)*bzg(j0, 1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, 0)*bzg(j0+1, 1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, 1)*bzg(j0-1, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 1)*bzg(j0, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, 1)*bzg(j0+1, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, 2)*bzg(j0-1, 1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 2)*bzg(j0, 1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, 2)*bzg(j0+1, 1, l+2)
      end do
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    end do

    ! ___________________________
    ! l_lower_order_in_v is false
  ELSE


    ! Loop over the particles by block
    DO ip=1, np, lvect

#if defined __INTEL_COMPILER
      !!DIR$ IVDEP
      !!DIR$ DISTRIBUTE POINT
      !DIR$ ASSUME_ALIGNED xp:64, zp:64
      !DIR$ ASSUME_ALIGNED sx:64, sz:64
      !DIR$ ASSUME_ALIGNED sx0:64, sz0:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
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
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)

        l=floor(z)
        l0=floor(z-stagger_shift)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, -1)*bxg(j-1, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, -1)*bxg(j, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, -1)*bxg(j+1, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, -1)*bxg(j+2, 1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, 0)*bxg(j-1, 1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 0)*bxg(j, 1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 0)*bxg(j+1, 1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, 0)*bxg(j+2, 1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, 1)*bxg(j-1, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 1)*bxg(j, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 1)*bxg(j+1, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, 1)*bxg(j+2, 1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sz0(n, 2)*bxg(j-1, 1, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 2)*bxg(j, 1, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 2)*bxg(j+1, 1, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sz0(n, 2)*bxg(j+2, 1, l0+2)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, -1)*byg(j0-1, 1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, -1)*byg(j0, 1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, -1)*byg(j0+1, 1, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sz0(n, -1)*byg(j0+2, 1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, 0)*byg(j0-1, 1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 0)*byg(j0, 1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, 0)*byg(j0+1, 1, l0)
        by(nn) = by(nn) + sx0(n, 2)*sz0(n, 0)*byg(j0+2, 1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, 1)*byg(j0-1, 1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 1)*byg(j0, 1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, 1)*byg(j0+1, 1, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sz0(n, 1)*byg(j0+2, 1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sz0(n, 2)*byg(j0-1, 1, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 2)*byg(j0, 1, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sz0(n, 2)*byg(j0+1, 1, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sz0(n, 2)*byg(j0+2, 1, l0+2)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, -1)*bzg(j0-1, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, -1)*bzg(j0, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, -1)*bzg(j0+1, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sz(n, -1)*bzg(j0+2, 1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, 0)*bzg(j0-1, 1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 0)*bzg(j0, 1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, 0)*bzg(j0+1, 1, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sz(n, 0)*bzg(j0+2, 1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, 1)*bzg(j0-1, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 1)*bzg(j0, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, 1)*bzg(j0+1, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sz(n, 1)*bzg(j0+2, 1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sz(n, 2)*bzg(j0-1, 1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 2)*bzg(j0, 1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sz(n, 2)*bzg(j0+1, 1, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sz(n, 2)*bzg(j0+2, 1, l+2)

      enddo
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    end do

  ENDIF
  return
end subroutine pxr_getb2dxz_energy_conserving_vect_3_3

! ________________________________________________________________________________________
!> @brief
!> Cartesian vectorized field gathering in 2D for the magnetic and the electric field
!> at order 3 in the same loop.
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 12/01/2016
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nx, nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of
!> the exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the
!> bxg, byg, bzg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 2 integers)
!> @param[in] lvect the vector length of the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_geteb2dxz_energy_conserving_vect_3_3( np, xp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
  ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,            &
  byg_nvalid, bzg, bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v, l_nodal)        !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  ! __ Parameter declaration ___________________________________________
  integer(idp), intent(in)                :: np
  integer(idp), intent(IN)      :: exg_nguard(2), exg_nvalid(2), eyg_nguard(2),       &
  eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2), bxg_nguard(2), bxg_nvalid(2),          &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  integer(idp), intent(in)                :: lvect
  real(num), dimension(np), intent(in)    :: xp, zp
  real(num), dimension(np), intent(inout) :: ex, ey, ez
  real(num), dimension(np), intent(inout) :: bx, by, bz
  logical(lp) , intent(in)                :: l_lower_order_in_v, l_nodal
  real(num)          :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1, 1,        &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1, 1,        &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, 1,        &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1, 1,        &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1, 1,        &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1, 1,        &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)
  real(num)                          :: xmin, zmin, dx, dz
  integer(idp)                       :: ip, j, l
  integer(idp)                       :: j0, l0
  integer(idp)                       :: n, nn
  real(num)                          :: dxi, dzi, x, z, xint, zint
  real(num)                          :: a
  real(num)                          :: xintsq, oxint, zintsq, ozint, oxintsq,        &
  ozintsq
  real(num), DIMENSION(lvect, -1:2)   :: sx, sx0
  real(num), DIMENSION(lvect, -1:2)   :: sz, sz0
  real(num), parameter               :: onesixth=1./6., twothird=2./3.

  ! ___________________________
  ! Compute parameters

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dzi = 1./dz

  sx=0
  sz=0.
  sx0=0.
  sz0=0.

  ! ___________________________
  IF (l_lower_order_in_v) THEN

    ! Loop over the particles by block
    DO ip=1, np, lvect

#if defined __INTEL_COMPILER
      !!DIR$ IVDEP
      !!DIR$ DISTRIBUTE POINT
      !DIR$ ASSUME_ALIGNED xp:64, zp:64
      !DIR$ ASSUME_ALIGNED sx:64, sz:64
      !DIR$ ASSUME_ALIGNED sx0:64, sz0:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
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
      DO n=1, MIN(lvect, np-ip+1)


        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)

        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2


        ! Compute Ex on particle
        a = (sx0(n, -1)*exg(j0-1, 1, l-1) + sx0(n, 0)*exg(j0, 1, l-1) + sx0(n,        &
        1)*exg(j0+1, 1, l-1))
        ex(nn) = ex(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*exg(j0-1, 1, l) + sx0(n, 0)*exg(j0, 1, l) + sx0(n,            &
        1)*exg(j0+1, 1, l))
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*exg(j0-1, 1, l+1) + sx0(n, 0)*exg(j0, 1, l+1) + sx0(n,        &
        1)*exg(j0+1, 1, l+1))
        ex(nn) = ex(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*exg(j0-1, 1, l+2) + sx0(n, 0)*exg(j0, 1, l+2) + sx0(n,        &
        1)*exg(j0+1, 1, l+2))
        ex(nn) = ex(nn) + a*sz(n, 2)

        ! Compute Ey on particle
        a = (sx(n, -1)*eyg(j-1, 1, l-1) + sx(n, 0)*eyg(j, 1, l-1) + sx(n, 1)*eyg(j+1, &
        1, l-1) + sx(n, 2)*eyg(j+2, 1, l-1))
        ey(nn) = ey(nn) + a*sz(n, -1)
        a = (sx(n, -1)*eyg(j-1, 1, l) + sx(n, 0)*eyg(j, 1, l) + sx(n, 1)*eyg(j+1, 1,  &
        l) + sx(n, 2)*eyg(j+2, 1, l))
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = (sx(n, -1)*eyg(j-1, 1, l+1) + sx(n, 0)*eyg(j, 1, l+1) + sx(n, 1)*eyg(j+1, &
        1, l+1) + sx(n, 2)*eyg(j+2, 1, l+1))
        ey(nn) = ey(nn) + a*sz(n, 1)
        a = (sx(n, -1)*eyg(j-1, 1, l+2) + sx(n, 0)*eyg(j, 1, l+2) + sx(n, 1)*eyg(j+1, &
        1, l+2) + sx(n, 2)*eyg(j+2, 1, l+2))
        ey(nn) = ey(nn) + a*sz(n, 2)

        ! Compute Ez on particle
        a = (sx(n, -1)*ezg(j-1, 1, l0-1) + sx(n, 0)*ezg(j, 1, l0-1) + sx(n,           &
        1)*ezg(j+1, 1, l0-1) + sx(n, 2)*ezg(j+2, 1, l0-1))
        ez(nn) = ez(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*ezg(j-1, 1, l0) + sx(n, 0)*ezg(j, 1, l0) + sx(n, 1)*ezg(j+1,   &
        1, l0) + sx(n, 2)*ezg(j+2, 1, l0))
        ez(nn) = ez(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*ezg(j-1, 1, l0+1) + sx(n, 0)*ezg(j, 1, l0+1) + sx(n,           &
        1)*ezg(j+1, 1, l0+1) + sx(n, 2)*ezg(j+2, 1, l0+1))
        ez(nn) = ez(nn) + a*sz0(n, 1)

        ! Compute Bx on particle
        a = (sx(n, -1)*bxg(j-1, 1, l0-1) + sx(n, 0)*bxg(j, 1, l0-1) + sx(n,           &
        1)*bxg(j+1, 1, l0-1) + sx(n, 2)*bxg(j+2, 1, l0-1))
        bx(nn) = bx(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*bxg(j-1, 1, l0) + sx(n, 0)*bxg(j, 1, l0) + sx(n, 1)*bxg(j+1,   &
        1, l0) + sx(n, 2)*bxg(j+2, 1, l0))
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*bxg(j-1, 1, l0+1) + sx(n, 0)*bxg(j, 1, l0+1) + sx(n,           &
        1)*bxg(j+1, 1, l0+1) + sx(n, 2)*bxg(j+2, 1, l0+1))
        bx(nn) = bx(nn) + a*sz0(n, 1)

        ! Compute By on particle
        a = (sx0(n, -1)*byg(j0-1, 1, l0-1) + sx0(n, 0)*byg(j0, 1, l0-1) + sx0(n,      &
        1)*byg(j0+1, 1, l0-1))
        by(nn) = by(nn) + a*sz0(n, -1)
        a = (sx0(n, -1)*byg(j0-1, 1, l0) + sx0(n, 0)*byg(j0, 1, l0) + sx0(n,          &
        1)*byg(j0+1, 1, l0))
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, -1)*byg(j0-1, 1, l0+1) + sx0(n, 0)*byg(j0, 1, l0+1) + sx0(n,      &
        1)*byg(j0+1, 1, l0+1))
        by(nn) = by(nn) + a*sz0(n, 1)

        ! Compute Bz on particle
        a = (sz(n, -1)*bzg(j0-1, 1, l-1) + sz(n, 0)*bzg(j0-1, 1, l) + sz(n,           &
        1)*bzg(j0-1, 1, l+1) + sz(n, 2)*bzg(j0-1, 1, l+2))
        bz(nn) = bz(nn) + a*sx0(n, -1)
        a = (sz(n, -1)*bzg(j0, 1, l-1) + sz(n, 0)*bzg(j0, 1, l) + sz(n, 1)*bzg(j0, 1, &
        l+1) + sz(n, 2)*bzg(j0, 1, l+2))
        bz(nn) = bz(nn) + a*sx0(n, 0)
        a = (sz(n, -1)*bzg(j0+1, 1, l-1) + sz(n, 0)*bzg(j0+1, 1, l) + sz(n,           &
        1)*bzg(j0+1, 1, l+1) + sz(n, 2)*bzg(j0+1, 1, l+2))
        bz(nn) = bz(nn) + a*sx0(n, 1)

      end do
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    end do

    ! ___________________________
    ! l_lower_order_in_v is false
  ELSE


    ! Loop over the particles by block
    DO ip=1, np, lvect

#if defined __INTEL_COMPILER
      !!DIR$ IVDEP
      !!DIR$ DISTRIBUTE POINT
      !DIR$ ASSUME_ALIGNED xp:64, zp:64
      !DIR$ ASSUME_ALIGNED sx:64, sz:64
      !DIR$ ASSUME_ALIGNED sx0:64, sz0:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
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
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)

        l=floor(z)
        l0=floor(z-stagger_shift)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        a = (sx0(n, -1)*exg(j0-1, 1, l-1) + sx0(n, 0)*exg(j0, 1, l-1) + sx0(n,        &
        1)*exg(j0+1, 1, l-1) + sx0(n, 2)*exg(j0+2, 1, l-1))
        ex(nn) = ex(nn) + a*sz(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, 1, l) + sx0(n, 0)*exg(j0, 1, l) + sx0(n,        &
        1)*exg(j0+1, 1, l) + sx0(n, 2)*exg(j0+2, 1, l))
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, 1, l+1) + sx0(n, 0)*exg(j0, 1, l+1) + sx0(n,    &
        1)*exg(j0+1, 1, l+1) + sx0(n, 2)*exg(j0+2, 1, l+1))
        ex(nn) = ex(nn) + a*sz(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, 1, l+2) + sx0(n, 0)*exg(j0, 1, l+2) + sx0(n,    &
        1)*exg(j0+1, 1, l+2) + sx0(n, 2)*exg(j0+2, 1, l+2))
        ex(nn) = ex(nn) + a*sz(n, 2)

        ! Compute Ey on particle
        a = (sx(n, -1)*eyg(j-1, 1, l-1) + sx(n, 0)*eyg(j, 1, l-1) + sx(n, 1)*eyg(j+1, &
        1, l-1) + sx(n, 2)*eyg(j+2, 1, l-1))
        ey(nn) = ey(nn) + a*sz(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, 1, l) + sx(n, 0)*eyg(j, 1, l) + sx(n, 1)*eyg(j+1, &
        1, l) + sx(n, 2)*eyg(j+2, 1, l))
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, 1, l+1) + sx(n, 0)*eyg(j, 1, l+1) + sx(n,         &
        1)*eyg(j+1, 1, l+1) + sx(n, 2)*eyg(j+2, 1, l+1))
        ey(nn) = ey(nn) + a*sz(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, 1, l+2) + sx(n, 0)*eyg(j, 1, l+2) + sx(n,         &
        1)*eyg(j+1, 1, l+2) + sx(n, 2)*eyg(j+2, 1, l+2))
        ey(nn) = ey(nn) + a*sz(n, 2)

        ! Compute Ez on particle
        a = (sx(n, -1)*ezg(j-1, 1, l0-1) + sx(n, 0)*ezg(j, 1, l0-1) + sx(n,           &
        1)*ezg(j+1, 1, l0-1) + sx(n, 2)*ezg(j+2, 1, l0-1))
        ey(nn) = ey(nn) + a*sz0(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, 1, l0) + sx(n, 0)*ezg(j, 1, l0) + sx(n,           &
        1)*ezg(j+1, 1, l0) + sx(n, 2)*ezg(j+2, 1, l0))
        ey(nn) = ey(nn) + a*sz0(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, 1, l0+1) + sx(n, 0)*ezg(j, 1, l0+1) + sx(n,       &
        1)*ezg(j+1, 1, l0+1) + sx(n, 2)*ezg(j+2, 1, l0+1))
        ey(nn) = ey(nn) + a*sz0(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, 1, l0+2) + sx(n, 0)*ezg(j, 1, l0+2) + sx(n,       &
        1)*ezg(j+1, 1, l0+2) + sx(n, 2)*ezg(j+2, 1, l0+2))
        ey(nn) = ey(nn) + a*sz0(n, 2)


        ! Compute Bx on particle
        a = (sx(n, -1)*bxg(j-1, 1, l0-1) + sx(n, 0)*bxg(j, 1, l0-1) + sx(n,           &
        1)*bxg(j+1, 1, l0-1) + sx(n, 2)*bxg(j+2, 1, l0-1))
        bx(nn) = bx(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*bxg(j-1, 1, l0) + sx(n, 0)*bxg(j, 1, l0) + sx(n, 1)*bxg(j+1,   &
        1, l0) + sx(n, 2)*bxg(j+2, 1, l0))
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*bxg(j-1, 1, l0+1) + sx(n, 0)*bxg(j, 1, l0+1) + sx(n,           &
        1)*bxg(j+1, 1, l0+1) + sx(n, 2)*bxg(j+2, 1, l0+1))
        bx(nn) = bx(nn) + a*sz0(n, 1)
        a = (sx(n, -1)*bxg(j-1, 1, l0+2) + sx(n, 0)*bxg(j, 1, l0+2) + sx(n,           &
        1)*bxg(j+1, 1, l0+2) + sx(n, 2)*bxg(j+2, 1, l0+2))
        bx(nn) = bx(nn) + a*sz0(n, 2)

        ! Compute By on particle
        a = (sx0(n, -1)*byg(j0-1, 1, l0-1) + sx0(n, 0)*byg(j0, 1, l0-1) + sx0(n,      &
        1)*byg(j0+1, 1, l0-1) + sx0(n, 2)*byg(j0+2, 1, l0-1))
        by(nn) = by(nn) + a*sz0(n, -1)
        a = (sx0(n, -1)*byg(j0-1, 1, l0) + sx0(n, 0)*byg(j0, 1, l0) + sx0(n,          &
        1)*byg(j0+1, 1, l0) + sx0(n, 2)*byg(j0+2, 1, l0))
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, -1)*byg(j0-1, 1, l0+1) + sx0(n, 0)*byg(j0, 1, l0+1) + sx0(n,      &
        1)*byg(j0+1, 1, l0+1) + sx0(n, 2)*byg(j0+2, 1, l0+1))
        by(nn) = by(nn) + a*sz0(n, 1)
        a = (sx0(n, -1)*byg(j0-1, 1, l0+2) + sx0(n, 0)*byg(j0, 1, l0+2) + sx0(n,      &
        1)*byg(j0+1, 1, l0+2) + sx0(n, 2)*byg(j0+2, 1, l0+2))
        by(nn) = by(nn) + a*sz0(n, 2)

        ! Compute Bz on particle
        a = (sx0(n, -1)*bzg(j0-1, 1, l-1) + sx0(n, 0)*bzg(j0, 1, l-1) + sx0(n,        &
        1)*bzg(j0+1, 1, l-1) + sx0(n, 2)*bzg(j0+2, 1, l-1))
        bx(nn) = bx(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*bzg(j0-1, 1, l) + sx0(n, 0)*bzg(j0, 1, l) + sx0(n,            &
        1)*bzg(j0+1, 1, l) + sx0(n, 2)*bzg(j0+2, 1, l))
        bx(nn) = bx(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*bzg(j0-1, 1, l+1) + sx0(n, 0)*bzg(j0, 1, l+1) + sx0(n,        &
        1)*bzg(j0+1, 1, l+1) + sx0(n, 2)*bzg(j0+2, 1, l+1))
        bx(nn) = bx(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*bzg(j0-1, 1, l+2) + sx0(n, 0)*bzg(j0, 1, l+2) + sx0(n,        &
        1)*bzg(j0+1, 1, l+2) + sx0(n, 2)*bzg(j0+2, 1, l+2))
        bx(nn) = bx(nn) + a*sz(n, 2)
      enddo
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    end do
  ENDIF
  return
end subroutine pxr_geteb2dxz_energy_conserving_vect_3_3

! ________________________________________________________________________________________
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
! FIELD_GATHERING_O1_2D.F90
!
! Field gathering subroutines in 2D at order 1
!
! List of subroutines:
!
! - pxr_gete2dxz_energy_conserving_scalar_1_1
! - pxr_getb2dxz_energy_conserving_scalar_1_1
! - pxr_gete2dxz_energy_conserving_vect_1_1
! - pxr_getb2dxz_energy_conserving_vect_1_1
!
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Scalar version: gathering of electric field from Yee grid ("energy conserving")
!> on particles at order 1.
!
!> @details
!> This subroutine is NOT vectorized.
!
!> @param[in] np number of particles
!> @param[in] xp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[in] xmin, zmin tile minimum grid position
!> @param[in] dx, dz space step
!> @param[in] dt time step
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!>  exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE pxr_gete2dxz_energy_conserving_scalar_1_1(np, xp, zp, ex, ey, ez, xmin,    &
  zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
  ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np
  INTEGER(idp), intent(in)             :: exg_nguard(2), exg_nvalid(2),               &
  eyg_nguard(2), eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  REAL(num), DIMENSION(np)             :: xp, zp, ex, ey, ez
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  REAL(num)                                :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  REAL(num)                            :: xmin, zmin, dx, dz
  INTEGER(idp)                         :: ip, j, l
  INTEGER(idp)                         :: ixmin, ixmax, izmin, izmax
  INTEGER(idp)                         :: ixmin0, ixmax0, izmin0,     &
  izmax0
  INTEGER(idp)                         :: jj, ll, j0, l0
  REAL(num)                            :: dxi, dzi, x, z
  REAL(num)                            :: xint, zint
  REAL(num), DIMENSION(0:1)            :: sx, sx0
  REAL(num), DIMENSION(0:1)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz

  ixmin = 0
  ixmax = 0
  izmin = 0
  izmax = 0

  IF (l_lower_order_in_v) THEN

    ixmin0 = 0
    ixmax0 = 0
    izmin0 = 0
    izmax0 = 0

!$acc parallel deviceptr(exg, eyg, ezg, xp, zp, ex, ey, ez)
!$acc loop gang vector private(sx(0:1), sz(0:1), sx0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle

      j=floor(x)
      j0=floor(x+0.5_num-stagger_shift)
      l=floor(z)
      l0=floor(z+0.5_num-stagger_shift)

      xint=x-j
      zint=z-l

      ! Compute shape factors
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sz( 0) = 1.0_num-zint
      sz( 1) = zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      sx0( 0) = 1.0_num
      sz0( 0) = 1.0_num
      sx0( 1) = 0.0_num
      sz0( 1) = 0.0_num

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j0+jj, l+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj, l+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj, l0+ll)
          end do
      end do
      !$acc end loop

    END DO
    !$acc end loop
    !$acc end parallel

    ! __ l_lower_order_in_v false  _____________________________
  ELSE

    ixmin0 = 0
    ixmax0 = 1
    izmin0 = 0
    izmax0 = 1

    !$acc parallel deviceptr(exg, eyg, ezg, xp, zp, ex, ey, ez)
    !$acc loop gang vector private(sx(0:1), sz(0:1), sx0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x-stagger_shift)
      l=floor(z)
      l0=floor(z-stagger_shift)
      xint=x-j
      zint=z-l

      ! Compute shape factors
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sz( 0) = 1.0_num-zint
      sz( 1) = zint
      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0
      sx0( 0) = 1.0_num-xint
      sx0( 1) = xint
      sz0( 0) = 1.0_num-zint
      sz0( 1) = zint

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j0+jj, l+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj, l+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj, l0+ll)
          end do
      end do
      !$acc end loop

    END DO
    !$acc end loop
    !$acc end parallel
  ENDIF
  RETURN
END SUBROUTINE pxr_gete2dxz_energy_conserving_scalar_1_1

! ______________________________________________________________________________
!> @brief
!> Scalar version: Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 1
!
!> @details
!> This function is NOT vectorized
!
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] bx, by, bz particle magnetic field arrays
!> @param[in] xmin, zmin tile minimum grid position
!> @param[in] dx, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] bxg, byg, bzg electric field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the bxg, byg, bzg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 2 integers)
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
! ________________________________________________________________________________________
SUBROUTINE pxr_getb2dxz_energy_conserving_scalar_1_1(np, xp, zp, bx, by, bz, xmin,    &
  zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
  bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np
  INTEGER(idp), intent(in)                :: bxg_nguard(2), bxg_nvalid(2),            &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  REAL(num), DIMENSION(np)             :: xp, zp, bx, by, bz
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  REAL(num)                                :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)
  REAL(num)                            :: xmin, zmin, dx, dz
  INTEGER(idp)                         :: ip, j, l, ixmin, ixmax,    &
  izmin, izmax, ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0,   &
  l0
  REAL(num)                            :: dxi, dzi, x, z, xint, zint
  REAL(num), DIMENSION(0:1)            :: sx, sx0
  REAL(num), DIMENSION(0:1)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num,                   &
  twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz

  ixmin = 0
  ixmax = 0
  izmin = 0
  izmax = 0

  IF (l_lower_order_in_v) THEN

    ixmin0 = 0
    ixmax0 = 0
    izmin0 = 0
    izmax0 = 0

!$acc parallel deviceptr(bxg, byg, bzg, xp, zp, bx, by, bz)
!$acc loop gang vector private(sx(0:1), sz(0:1), sx0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x+0.5_num-stagger_shift)
      l=floor(z)
      l0=floor(z+0.5_num-stagger_shift)

      xint=x-j
      zint=z-l

      ! Compute shape factors
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sz( 0) = 1.0_num-zint
      sz( 1) = zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      sx0( 0) = 1.0_num
      sz0( 0) = 1.0_num
      sx0( 1) = 0.0_num
      sz0( 1) = 0.0_num

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj, l0+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj, l0+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j0+jj, l+ll)
          end do
      end do
      !$acc end loop

    END DO
    !$acc end loop
    !$acc end parallel

  ELSE

    ixmin0 = 0
    ixmax0 = 1
    izmin0 = 0
    izmax0 = 1

    !$acc parallel deviceptr(bxg, byg, bzg, xp, zp, bx, by, bz)
    !$acc loop gang vector private(sx(0:1), sz(0:1), sx0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x-stagger_shift)
      l=floor(z)
      l0=floor(z-stagger_shift)

      ! Compute shape factors
      xint=x-j
      zint=z-l
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sz( 0) = 1.0_num-zint
      sz( 1) = zint

      xint=x-stagger_shift-j0
      zint=z-stagger_shift-l0

      sx0( 0) = 1.0_num-xint
      sx0( 1) = xint
      sz0( 0) = 1.0_num-zint
      sz0( 1) = zint

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj, l0+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin0, izmax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj, l0+ll)
          end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(2)
      do ll = izmin, izmax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j0+jj, l+ll)
          end do
      end do
      !$acc end loop
    END DO
    !$acc end loop
    !$acc end parallel
  ENDIF
  RETURN
END SUBROUTINE pxr_getb2dxz_energy_conserving_scalar_1_1

! ________________________________________________________________________________________
!> @brief
!> Field gathering cartesian in 2D for the electric field at order 1
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
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] lvect vector size for the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_gete2dxz_energy_conserving_vect_1_1( np, xp, zp, ex, ey, ez, xmin,     &
  zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,          &
  ezg_nguard, ezg_nvalid, lvect, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp), intent(in)                :: np
  integer(idp), intent(IN)                :: exg_nguard(2), exg_nvalid(2),            &
  eyg_nguard(2), eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  integer(idp), intent(in)                :: lvect
  real(num), dimension(np), intent(in)    :: xp, zp
  real(num), dimension(np), intent(inout) :: ex, ey, ez
  logical(lp)                             :: l_lower_order_in_v, l_nodal
  real(num)                               :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1, 1,        &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1, 1,        &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, 1,        &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  real(num)                          :: xmin, zmin, dx, dz
  integer(idp)                       :: ip, j, l, j0, l0
  integer(idp)                       :: nn, n
  real(num)                          :: dxi, dzi, x, z, xint, zint
  real(num), DIMENSION(lvect, 0:1)    :: sx, sx0
  real(num), DIMENSION(lvect, 0:1)    :: sz, sz0
  real(num), parameter               :: onesixth=1./6., twothird=2./3.

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
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num

        sz0(n, 0) = 1.0_num

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 0)*exg(j0, 1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 1)*exg(j0, 1, l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 0)*eyg(j, 1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 0)*eyg(j+1, 1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 1)*eyg(j, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 1)*eyg(j+1, 1, l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 0)*ezg(j, 1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 0)*ezg(j+1, 1, l0)

      ENDDO
    ENDDO

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

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint

        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 0)*exg(j0, 1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 0)*exg(j0+1, 1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sz(n, 1)*exg(j0, 1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sz(n, 1)*exg(j0+1, 1, l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 0)*eyg(j, 1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 0)*eyg(j+1, 1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sz(n, 1)*eyg(j, 1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sz(n, 1)*eyg(j+1, 1, l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 0)*ezg(j, 1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 0)*ezg(j+1, 1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sz0(n, 1)*ezg(j, 1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sz0(n, 1)*ezg(j+1, 1, l0+1)

      ENDDO
    ENDDO

  end if

end subroutine


! ________________________________________________________________________________________
!> @brief
!> Field gathering cartesian in 2D for the magnetic field at order 1
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
!> @param[in] lvect the vector length of the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_getb2dxz_energy_conserving_vect_1_1( np, xp, zp, bx, by, bz, xmin,     &
  zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,          &
  bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  ! __ Parameter declaration ___________________________________________
  integer(idp), intent(in)                :: np
  integer(idp), intent(IN)                :: bxg_nguard(2), bxg_nvalid(2),            &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  integer(idp), intent(in)                :: lvect
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
  real(num)                          :: xmin, zmin, dx, dz
  integer(idp)                       :: ip, j, l
  integer(idp)                       :: j0, l0
  integer(idp)                       :: n, nn
  real(num)                          :: dxi, dzi, x, z, xint, zint
  real(num), DIMENSION(lvect, 0:1)    :: sx, sx0
  real(num), DIMENSION(lvect, 0:1)    :: sz, sz0
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
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num
        sz0(n, 0) = 1.0_num

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 0)*bxg(j, 1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 0)*bxg(j+1, 1, l0)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 0)*byg(j0, 1, l0)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 0)*bzg(j0, 1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 1)*bzg(j0, 1, l+1)

      ENDDO
    ENDDO

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
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint

        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, 0)*sz0(n, 0)*bxg(j, 1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sz0(n, 0)*bxg(j+1, 1, l0)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, 0)*sz0(n, 0)*byg(j0, 1, l0)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 0)*bzg(j0, 1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sz(n, 1)*bzg(j0, 1, l+1)

      ENDDO
    ENDDO
  ENDIF
  return
end subroutine

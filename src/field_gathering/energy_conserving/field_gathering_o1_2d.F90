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
! FIELD_GATHERING_O1_3D.F90
!
! Field gathering subroutines in 2D at order 1
!
! List of subroutines:
!
! - pxr_gete2dxz_energy_conserving_vect_1_1
! - pxr_getb2dxz_energy_conserving_vect_1_1
!
! ________________________________________________________________________________________


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
  ezg_nguard, ezg_nvalid, lvect, l_lower_order_in_v)     !#do not wrap
  USE PICSAR_precision
  USE constants
  implicit none

  integer(idp), intent(in)                :: np
  integer(idp), intent(IN)                :: exg_nguard(2), exg_nvalid(2),            &
  eyg_nguard(2), eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  integer(idp), intent(in)                :: lvect
  real(num), dimension(np), intent(in)    :: xp, zp
  real(num), dimension(np), intent(inout) :: ex, ey, ez
  logical(idp)                            :: l_lower_order_in_v
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1, 1,        &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1, 1,        &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, 1,        &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  real(num)                          :: xmin, zmin, dx, dz
  integer(isp)                       :: ip, j, l, j0, l0
  integer(isp)                       :: nn, n
  real(num)                          :: dxi, dzi, x, z, xint, zint
  real(num), DIMENSION(lvect, 0:1)    :: sx, sx0
  real(num), DIMENSION(lvect, 0:1)    :: sz, sz0
  real(num), parameter               :: onesixth=1./6., twothird=2./3.

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
        j0=floor(x)
        l=floor(z)
        l0=floor(z)
        xint=x-j
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

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
        j0=floor(x-0.5_num)
        l=floor(z)
        l0=floor(z-0.5_num)
        xint=x-j
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

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
  bzg_nguard, bzg_nvalid, lvect, l_lower_order_in_v)     !#do not wrap
  USE PICSAR_precision
  USE constants
  implicit none
  ! __ Parameter declaration ___________________________________________
  integer(idp), intent(in)                :: np
  integer(idp), intent(IN)                :: bxg_nguard(2), bxg_nvalid(2),            &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  integer(idp), intent(in)                :: lvect
  real(num), dimension(np), intent(in)    :: xp, zp
  real(num), dimension(np), intent(inout) :: bx, by, bz
  logical(idp), intent(in)                :: l_lower_order_in_v
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
        j0=floor(x)

        l=floor(z)
        l0=floor(z)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

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
        j0=floor(x-0.5_num)

        l=floor(z)
        l0=floor(z-0.5_num)

        xint=x-j
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-0.5_num-j0
        zint=z-0.5_num-l0

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

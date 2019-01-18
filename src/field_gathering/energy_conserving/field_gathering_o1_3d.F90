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
! Developers:
! - Henri vincenti
! - Mathieu Lobet
!
! Field gathering subroutines in 3D at order 1
!
! List of subroutines:
!
! - gete3d_energy_conserving_scalar_1_1_1
! - getb3d_energy_conserving_scalar_1_1_1
!
! - gete3d_energy_conserving_vec_1_1_1
! - getb3d_energy_conserving_vec_1_1_1
!
! - geteb3d_energy_conserving_vecV1_1_1_1
! - geteb3d_energy_conserving_vecV2_1_1_1
! - geteb3d_energy_conserving_vecV3_1_1_1
! - geteb3d_energy_conserving_vec_1_1_1_v2
! - geteb3d_energy_conserving_vec_1_1_1_sub
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
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!>  exg, eyg, ezg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 3 integers)
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE gete3d_energy_conserving_scalar_1_1_1(np, xp, yp, zp, ex, ey, ez, xmin,    &
  ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
  ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np
  INTEGER(idp), intent(in)             :: exg_nguard(3), exg_nvalid(3),               &
  eyg_nguard(3), eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3)
  REAL(num), DIMENSION(np)             :: xp, yp, zp, ex, ey, ez
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  REAL(num)                                :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1,                                       &
  -exg_nguard(3):exg_nvalid(3)+exg_nguard(3)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1,                                       &
  -eyg_nguard(3):eyg_nvalid(3)+eyg_nguard(3)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1,                                       &
  -ezg_nguard(3):ezg_nvalid(3)+ezg_nguard(3)-1)
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: ixmin, ixmax, iymin, iymax, izmin, izmax
  INTEGER(isp)                         :: ixmin0, ixmax0, iymin0, iymax0, izmin0,     &
  izmax0
  INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num), DIMENSION(0:1)            :: sx, sx0
  REAL(num), DIMENSION(0:1)            :: sy, sy0
  REAL(num), DIMENSION(0:1)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

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

!$acc parallel deviceptr(exg, eyg, ezg, xp, yp, zp, ex, ey, ez)
!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      y = (yp(ip)-ymin)*dyi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle

      j=floor(x)
      j0=floor(x+0.5_num-stagger_shift)
      k=floor(y)
      k0=floor(y+0.5_num-stagger_shift)
      l=floor(z)
      l0=floor(z+0.5_num-stagger_shift)

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

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      sx0( 0) = 1.0_num
      sy0( 0) = 1.0_num
      sz0( 0) = 1.0_num
      sx0( 1) = 0.0_num
      sy0( 1) = 0.0_num
      sz0( 1) = 0.0_num

      !$acc loop seq independent collapse(3)
      do ll = izmin, izmax+1
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj, k+kk, l+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin, izmax+1
        do kk = iymin0, iymax0
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj, k0+kk, l+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj, k+kk, l0+ll)
          end do
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
    iymin0 = 0
    iymax0 = 1
    izmin0 = 0
    izmax0 = 1

    !$acc parallel deviceptr(exg, eyg, ezg, xp, yp, zp, ex, ey, ez)
    !$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      y = (yp(ip)-ymin)*dyi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x-stagger_shift)
      k=floor(y)
      k0=floor(y-stagger_shift)
      l=floor(z)
      l0=floor(z-stagger_shift)
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
      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0
      sx0( 0) = 1.0_num-xint
      sx0( 1) = xint
      sy0( 0) = 1.0_num-yint
      sy0( 1) = yint
      sz0( 0) = 1.0_num-zint
      sz0( 1) = zint

      !$acc loop seq independent collapse(3)
      do ll = izmin, izmax+1
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj, k+kk, l+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin, izmax+1
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj, k0+kk, l+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj, k+kk, l0+ll)
          end do
        end do
      end do
      !$acc end loop

    END DO
    !$acc end loop
    !$acc end parallel
  ENDIF
  RETURN
END SUBROUTINE gete3d_energy_conserving_scalar_1_1_1

! ______________________________________________________________________________
!> @brief
!> Scalar version: Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> at order 1
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] bx, by, bz particle magnetic field arrays
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] bxg, byg, bzg electric field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the bxg, byg, bzg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 3 integers)
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
! ________________________________________________________________________________________
SUBROUTINE getb3d_energy_conserving_scalar_1_1_1(np, xp, yp, zp, bx, by, bz, xmin,    &
  ymin, zmin, dx, dy, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
  bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np
  INTEGER(idp), intent(in)                :: bxg_nguard(3), bxg_nvalid(3),            &
  byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  REAL(num), DIMENSION(np)             :: xp, yp, zp, bx, by, bz
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  REAL(num)                                :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1,                                       &
  -bxg_nguard(3):bxg_nvalid(3)+bxg_nguard(3)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1,                                       &
  -byg_nguard(3):byg_nvalid(3)+byg_nguard(3)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1,                                       &
  -bzg_nguard(3):bzg_nvalid(3)+bzg_nguard(3)-1)
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(idp)                         :: ip, j, k, l, ixmin, ixmax, iymin, iymax,    &
  izmin, izmax, ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0,   &
  l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num), DIMENSION(0:1)            :: sx, sx0
  REAL(num), DIMENSION(0:1)            :: sy, sy0
  REAL(num), DIMENSION(0:1)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num,                   &
  twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

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

!$acc parallel deviceptr(bxg, byg, bzg, xp, yp, zp, bx, by, bz)
!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      y = (yp(ip)-ymin)*dyi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x+0.5_num-stagger_shift)
      k=floor(y)
      k0=floor(y+0.5_num-stagger_shift)
      l=floor(z)
      l0=floor(z+0.5_num-stagger_shift)

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

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      sx0( 0) = 1.0_num
      sy0( 0) = 1.0_num
      sz0( 0) = 1.0_num
      sx0( 1) = 0.0_num
      sy0( 1) = 0.0_num
      sz0( 1) = 0.0_num

      !$acc loop seq independent collapse(3)
      do ll = izmin0, izmax0
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj, k0+kk, l0+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj, k+kk, l0+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin, izmax+1
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj, k0+kk, l+ll)
          end do
        end do
      end do
      !$acc end loop

    END DO
    !$acc end loop
    !$acc end parallel

  ELSE

    ixmin0 = 0
    ixmax0 = 1
    iymin0 = 0
    iymax0 = 1
    izmin0 = 0
    izmax0 = 1

    !$acc parallel deviceptr(bxg, byg, bzg, xp, yp, zp, bx, by, bz)
    !$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))
    DO ip=1, np

      x = (xp(ip)-xmin)*dxi
      y = (yp(ip)-ymin)*dyi
      z = (zp(ip)-zmin)*dzi

      ! Compute index of particle
      j=floor(x)
      j0=floor(x-stagger_shift)
      k=floor(y)
      k0=floor(y-stagger_shift)
      l=floor(z)
      l0=floor(z-stagger_shift)

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

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      sx0( 0) = 1.0_num-xint
      sx0( 1) = xint
      sy0( 0) = 1.0_num-yint
      sy0( 1) = yint
      sz0( 0) = 1.0_num-zint
      sz0( 1) = zint

      !$acc loop seq independent collapse(3)
      do ll = izmin0, izmax0
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin, ixmax+1
            bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj, k0+kk, l0+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin0, izmax0
        do kk = iymin, iymax+1
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj, k+kk, l0+ll)
          end do
        end do
      end do
      !$acc end loop

      !$acc loop seq independent collapse(3)
      do ll = izmin, izmax+1
        do kk = iymin0, iymax0
          ! Prevent wrong vectorization from the compiler
          !DIR$ NOVECTOR
          do jj = ixmin0, ixmax0
            bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj, k0+kk, l+ll)
          end do
        end do
      end do
      !$acc end loop
    END DO
    !$acc end loop
    !$acc end parallel
  ENDIF
  RETURN
END SUBROUTINE getb3d_energy_conserving_scalar_1_1_1

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Gathering of electric field from Yee grid ("energy conserving") on particles
!> at order 1.
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] ex, ey, ez particle electric field arrays
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v lower order for the interpolation
! ________________________________________________________________________________________
SUBROUTINE gete3d_energy_conserving_vec_1_1_1(np, xp, yp, zp, ex, ey, ez, xmin, ymin, &
  zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg, ezg, lvect,        &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                           :: lvect
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez
  LOGICAL(lp)                                :: l_lower_order_in_v, l_nodal
  REAL(num)                                  :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, j, k, l
  INTEGER(isp)                           :: jj, kk, ll, j0, k0, l0
  INTEGER(isp)                           :: n, nn
  REAL(num)                              :: dxi, dyi, dzi, x, y, z
  REAL(num)                              :: xint, yint, zint, xintsq, oxint, yintsq,  &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(lvect, 0:1)        :: sx, sx0
  REAL(num), DIMENSION(lvect, 0:1)        :: sy, sy0
  REAL(num), DIMENSION(lvect, 0:1)        :: sz, sz0
  REAL(num), PARAMETER                   :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                   :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

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

    ! Loop over the particles by block
    DO ip=1, np, lvect

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num
        sy0(n, 0) = 1.0_num
        sz0(n, 0) = 1.0_num

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

    ! l_lower_order_in_v false
  ELSE

    ! Loop over the particles by block
    DO ip=1, np, lvect

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)
        k=floor(y)
        k0=floor(y-stagger_shift)
        l=floor(z)
        l0=floor(z-stagger_shift)
        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint
        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0
        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint
        sy0(n, 0) = 1.0_num-yint
        sy0(n, 1) = yint
        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0+1, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0+1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0+1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0+1, k+1, l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j+1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j+1, k0+1, l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j+1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j+1, k+1, l0+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
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
!> at order 1.
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] bx, by, bz particle magnetic field arrays
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
! ________________________________________________________________________________________
SUBROUTINE getb3d_energy_conserving_vec_1_1_1(np, xp, yp, zp, bx, by, bz, xmin, ymin, &
  zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, bxg, byg, bzg, lvect,        &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp, yp, zp, bx, by, bz
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  REAL(num)                                :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: jj, kk, ll
  INTEGER(isp)                         :: j0, k0, l0
  INTEGER(isp)                         :: n, nn
  REAL(num)                            :: dxi, dyi, dzi
  REAL(num)                            :: x, y, z, xint, yint, zint, xintsq, oxint,   &
  yintsq, oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(lvect, 0:1)      :: sx, sx0
  REAL(num), DIMENSION(lvect, 0:1)      :: sy, sy0
  REAL(num), DIMENSION(lvect, 0:1)      :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num,                   &
  twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

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

    sx0(:, 0) = 1.0_num
    sy0(:, 0) = 1.0_num
    sz0(:, 0) = 1.0_num

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
      ! #if defined __INTEL_COMPILER
      ! !DIR IVDEP
      ! !DIR DISTRIBUTE POINT
      ! #endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

  ELSE

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
      ! #if defined __INTEL_COMPILER
      ! !DIR IVDEP
      ! !DIR DISTRIBUTE POINT
      ! #endif
      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)
        k=floor(y)
        k0=floor(y-stagger_shift)
        l=floor(z)
        l0=floor(z-stagger_shift)

        ! Compute shape factors
        xint=x-j
        yint=y-k
        zint=z-l
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint
        sy0(n, 0) = 1.0_num-yint
        sy0(n, 1) = yint
        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j+1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j+1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j+1, k0+1, l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0+1, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0+1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0, k, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0+1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0+1, k+1, l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0+1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0+1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0+1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0+1, k0+1, l+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE
#endif

! ______________________________________________________________________________
!> @brief
!> Field gathering CIC (order 1) with gathering of E and B merged in a single loop
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vecV1_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg,    &
  ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! __ Parameter declaration _____________________________________________________
  INTEGER(idp)                         :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp, yp, zp, ex, ey, ez, bx, by, bz
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  REAL(num)                                :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg, bxg, byg, bzg
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: nn, n
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num), DIMENSION(lvect, 0:1)      :: sx, sx0
  REAL(num), DIMENSION(lvect, 0:1)      :: sy, sy0
  REAL(num), DIMENSION(lvect, 0:1)      :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

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

    sx0(:, 0) = 1.0_num
    sy0(:, 0) = 1.0_num
    sz0(:, 0) = 1.0_num

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle

        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

  ELSE

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)
        k=floor(y)
        k0=floor(y-stagger_shift)
        l=floor(z)
        l0=floor(z-stagger_shift)
        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint
        sy0(n, 0) = 1.0_num-yint
        sy0(n, 1) = yint
        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0+1, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0+1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0+1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0+1, k+1, l+1)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j+1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j+1, k0+1, l+1)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j+1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j+1, k+1, l0+1)

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j+1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j+1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j+1, k0+1, l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0+1, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0+1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0, k, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0+1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0+1, k+1, l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0+1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0+1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0+1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0+1, k0+1, l+1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE

! ________________________________________________________________________________________
!
!> @brief
!> Field gathering (order 1) with gathering of E and B merged in a single loop.
!
!> @detail
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vecV2_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg,    &
  ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np
  INTEGER(idp)                           :: nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                                :: l_lower_order_in_v, l_nodal
  REAL(num)                                  :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp), DIMENSION(lvect)         :: j, k, l
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j, k, l
#endif
  INTEGER(isp), DIMENSION(lvect)         :: j0, k0, l0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j0, k0, l0
#endif
  REAL(num)                              :: x, y, z, xint, yint, zint
  REAL(num)                              :: dxi, dyi, dzi
  INTEGER(isp)                           :: n, ip, nn
  REAL(num), DIMENSION(lvect, 0:1)        :: sx, sy, sz
  REAL(num), DIMENSION(lvect, 0:1)        :: sx0, sy0, sz0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: sx, sy, sz, sx0, sy0, sz0
#endif
  REAL(num), PARAMETER                   :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                   :: twothird=2.0_num/3.0_num


  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  ! ___ Loop on all the partciles by block of size lvect _______________________

  DO ip=1, np, lvect

    ! ____________________________________________________________________________________
    ! Field gathering

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! Compute index of particle
      j(n)=floor(x)
      j0(n)=floor(x+0.5_num-stagger_shift)
      k(n)=floor(y)
      k0(n)=floor(y+0.5_num-stagger_shift)
      l(n)=floor(z)
      l0(n)=floor(z+0.5_num-stagger_shift)

      xint=x-j(n)
      yint=y-k(n)
      zint=z-l(n)

      ! Compute shape factors
      sx(n, 0) = 1.0_num-xint
      sx(n, 1) = xint
      sy(n, 0) = 1.0_num-yint
      sy(n, 1) = yint
      sz(n, 0) = 1.0_num-zint
      sz(n, 1) = zint

      xint=x-stagger_shift-j0(n)
      yint=y-stagger_shift-k0(n)
      zint=z-stagger_shift-l0(n)

      sx0(n, 0) = 1.0_num
      sy0(n, 0) = 1.0_num
      sz0(n, 0) = 1.0_num
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0(n), k(n), l(n))
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0(n), k(n)+1, l(n))
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0(n), k(n), l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0(n), k(n)+1, l(n)+1)

      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j(n), k0(n), l(n))
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j(n)+1, k0(n), l(n))
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j(n), k0(n), l(n)+1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j(n)+1, k0(n), l(n)+1)

      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j(n), k(n), l0(n))
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j(n)+1, k(n), l0(n))
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j(n), k(n)+1, l0(n))
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j(n)+1, k(n)+1, l0(n))

    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, bx, by, bz)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j(n), k0(n), l0(n))
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j(n)+1, k0(n), l0(n))

      ! Compute By on particle
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0(n), k(n), l0(n))
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0(n), k(n)+1, l0(n))

      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0(n), k0(n), l(n))
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0(n), k0(n), l(n)+1)

    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

  ENDDO

  RETURN
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> Field gathering CIC (order 1) with gathering of E and B merged in a single loop
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vecV3_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg,    &
  ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE

  ! __ Parameter declaration _____________________________________________________
  INTEGER(idp), intent(in)                         :: np, nx, ny, nz, nxguard,        &
  nyguard, nzguard
  INTEGER(idp), intent(in)                         :: lvect
  REAL(num), DIMENSION(np), intent(in)             :: xp, yp, zp
  REAL(num), DIMENSION(np), intent(inout)          :: ex, ey, ez, bx, by, bz
  LOGICAL(lp), intent(in)                          :: l_lower_order_in_v, l_nodal
  real(num)          :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), intent(in) :: exg, eyg, ezg, bxg, byg, bzg

  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: nn, n
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: a
  REAL(num), DIMENSION(lvect, 0:1)      :: sx, sx0
  REAL(num), DIMENSION(lvect, 0:1)      :: sy, sy0
  REAL(num), DIMENSION(lvect, 0:1)      :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

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

    sx0(:, 0) = 1.0_num
    sy0(:, 0) = 1.0_num
    sz0(:, 0) = 1.0_num

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle

        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        ! Compute Ex on particle
        a = (sy(n, 0)*exg(j0, k, l)   + sy(n, 1)*exg(j0, k+1, l))*sz(n, 0) + (sy(n,   &
        0)*exg(j0, k, l+1) + sy(n, 1)*exg(j0, k+1, l+1))*sz(n, 1)
        ex(nn) = ex(nn) + a

        ! Compute Ey on particle
        a = (sx(n, 0)*eyg(j, k0, l) + sx(n, 1)*eyg(j+1, k0, l))*sz(n, 0) + (sx(n,     &
        0)*eyg(j, k0, l+1) + sx(n, 1)*eyg(j+1, k0, l+1))*sz(n, 1)
        ey(nn) = ey(nn) + a

        ! Compute Ez on particle
        a = (sx(n, 0)*ezg(j, k, l0) + sx(n, 1)*ezg(j+1, k, l0))*sy(n, 0) + (sx(n,     &
        0)*ezg(j, k+1, l0) + sx(n, 1)*ezg(j+1, k+1, l0))*sy(n, 1)
        ez(nn) = ez(nn) + a

        ! Compute Bx on particle
        a = (sx(n, 0)*bxg(j, k0, l0) + sx(n, 1)*bxg(j+1, k0, l0))
        bx(nn) = bx(nn) + a

        ! Compute By on particle
        a = (sy(n, 0)*byg(j0, k, l0) + sy(n, 1)*byg(j0, k+1, l0))
        by(nn) = by(nn) + a

        ! Compute Bz on particle
        a = (sz(n, 0)*bzg(j0, k0, l) + sz(n, 1)*bzg(j0, k0, l+1))
        bz(nn) = bz(nn) + a

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

  ELSE


    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO n=1, MIN(lvect, np-ip+1)

        nn=ip+n-1

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)
        k=floor(y)
        k0=floor(y-stagger_shift)
        l=floor(z)
        l0=floor(z-stagger_shift)
        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(n, 0) = 1.0_num-xint
        sx(n, 1) = xint
        sy(n, 0) = 1.0_num-yint
        sy(n, 1) = yint
        sz(n, 0) = 1.0_num-zint
        sz(n, 1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        sx0(n, 0) = 1.0_num-xint
        sx0(n, 1) = xint
        sy0(n, 0) = 1.0_num-yint
        sy0(n, 1) = yint
        sz0(n, 0) = 1.0_num-zint
        sz0(n, 1) = zint

        ! Compute Ex on particle
        a = (sx0(n, 0)*exg(j0, k, l) + sx0(n, 1)*exg(j0+1, k, l))*sy(n, 0)
        a = a + (sx0(n, 0)*exg(j0, k+1, l) + sx0(n, 1)*exg(j0+1, k+1, l))*sy(n, 1)
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = (sx0(n, 0)*exg(j0, k, l+1) + sx0(n, 1)*exg(j0+1, k, l+1))*sy(n, 0)
        a = a + (sx0(n, 0)*exg(j0, k+1, l+1) + sx0(n, 1)*exg(j0+1, k+1, l+1))*sy(n,   &
        1)
        ex(nn) = ex(nn) + a*sz(n, 1)

        ! Compute Ey on particle
        a = (sx(n, 0)*eyg(j, k0, l) + sx(n, 1)*eyg(j+1, k0, l))*sy0(n, 0)
        a = a + (sx(n, 0)*eyg(j, k0+1, l) + sx(n, 1)*eyg(j+1, k0+1, l))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = (sx(n, 0)*eyg(j, k0, l+1) + sx(n, 1)*eyg(j+1, k0, l+1))*sy0(n, 0)
        a = a + (sx(n, 0)*eyg(j, k0+1, l+1) + sx(n, 1)*eyg(j+1, k0+1, l+1))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 1)

        ! Compute Ez on particle
        a = (sx(n, 0)*ezg(j, k, l0) + sx(n, 1)*ezg(j+1, k, l0))*sy(n, 0)
        a = a + (sx(n, 0)*ezg(j, k+1, l0) + sx(n, 1)*ezg(j+1, k+1, l0))*sy(n, 1)
        ez(nn) = ez(nn) + a*sz0(n, 0)
        a = (sx(n, 0)*ezg(j, k, l0+1) + sx(n, 1)*ezg(j+1, k, l0+1))*sy(n, 0)
        a = a + (sx(n, 0)*ezg(j, k+1, l0+1) + sx(n, 1)*ezg(j+1, k+1, l0+1))*sy(n, 1)
        ez(nn) = ez(nn) + a*sz0(n, 1)

        ! Compute Bx on particle
        a = (sx(n, 0)*bxg(j, k0, l0) + sx(n, 1)*bxg(j+1, k0, l0))*sy0(n, 0)
        a = a + (sx(n, 0)*bxg(j, k0+1, l0) + sx(n, 1)*bxg(j+1, k0+1, l0))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, 0)*bxg(j, k0, l0+1) + sx(n, 1)*bxg(j+1, k0, l0+1))*sy0(n, 0)
        a = a + (sx(n, 0)*bxg(j, k0+1, l0+1) + sx(n, 1)*bxg(j+1, k0+1, l0+1))*sy0(n,  &
        1)
        bx(nn) = bx(nn) + a*sz0(n, 1)

        ! Compute By on particle
        a = (sx0(n, 0)*byg(j0, k, l0) + sx0(n, 1)*byg(j0+1, k, l0))*sy(n, 0)
        a = a + (sx0(n, 0)*byg(j0, k+1, l0) + sx0(n, 1)*byg(j0+1, k+1, l0))*sy(n, 1)
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, 0)*byg(j0, k, l0+1) + sx0(n, 1)*byg(j0+1, k, l0+1))*sy(n, 0)
        a = a + (sx0(n, 0)*byg(j0, k+1, l0+1) + sx0(n, 1)*byg(j0+1, k+1, l0+1))*sy(n, &
        1)
        by(nn) = by(nn) + a*sz0(n, 1)

        ! Compute Bz on particle
        a = (sx0(n, 0)*bzg(j0, k0, l) + sx0(n, 1)*bzg(j0+1, k0, l))*sy0(n, 0)
        a = a + (sx0(n, 0)*bzg(j0, k0+1, l) + sx0(n, 1)*bzg(j0+1, k0+1, l))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 0)
        a = (sx0(n, 0)*bzg(j0, k0, l+1) + sx0(n, 1)*bzg(j0+1, k0, l+1))*sy0(n, 0)
        a = a + (sx0(n, 0)*bzg(j0, k0+1, l+1) + sx0(n, 1)*bzg(j0+1, k0+1,             &
        l+1))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> Field gathering CIC (order 1) with gathering of E and B merged
!> in a single loop.
!> this 4th version uses !$OMP SIMD private for sx, sy, sz, sx0, sy0, sz0
!> instead of having a new dimension
!
!> @details
!> This function is vectorized.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Revison 12/04/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid grid points
!> @param[in] bxg, byg, bzg electric field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid grid points
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vecV4_1_1_1(np, xp, yp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,       &
  eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,            &
  byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, lvect,                           &
  l_lower_order_in_v, l_nodal)        !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! __ Parameter declaration _____________________________________________________
  INTEGER(idp), intent(in)                :: np
  INTEGER(idp), intent(in)                :: exg_nguard(3), exg_nvalid(3),            &
  eyg_nguard(3), eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3), bxg_nguard(3),          &
  bxg_nvalid(3), byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  INTEGER(idp), intent(in)                :: lvect
  REAL(num), DIMENSION(np), intent(in)    :: xp, yp, zp
  REAL(num), DIMENSION(np), intent(inout) :: ex, ey, ez, bx, by, bz
  LOGICAL(lp), intent(in)                 :: l_lower_order_in_v, l_nodal
  real(num)          :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1,                                       &
  -exg_nguard(3):exg_nvalid(3)+exg_nguard(3)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1,                                       &
  -eyg_nguard(3):eyg_nvalid(3)+eyg_nguard(3)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1,                                       &
  -ezg_nguard(3):ezg_nvalid(3)+ezg_nguard(3)-1)
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1,                                       &
  -bxg_nguard(3):bxg_nvalid(3)+bxg_nguard(3)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1,                                       &
  -byg_nguard(3):byg_nvalid(3)+byg_nguard(3)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1,                                       &
  -bzg_nguard(3):bzg_nvalid(3)+bzg_nguard(3)-1)

  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: nn
  INTEGER(isp)                         :: j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: a
  REAL(num), DIMENSION(0:1)            :: sx, sx0
  REAL(num), DIMENSION(0:1)            :: sy, sy0
  REAL(num), DIMENSION(0:1)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

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

    sx0(0) = 1.0_num
    sy0(0) = 1.0_num
    sz0(0) = 1.0_num

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD private(sx, sy, sz, sx0, sy0, sz0)
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO nn=ip, MIN(ip+lvect-1, np)

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle

        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)

        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(0) = 1.0_num-xint
        sx(1) = xint
        sy(0) = 1.0_num-yint
        sy(1) = yint
        sz(0) = 1.0_num-zint
        sz(1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        ! Compute Ex on particle
        a = (sy(0)*exg(j0, k, l)   + sy(1)*exg(j0, k+1, l))*sz(0) + (sy(0)*exg(j0, k, &
        l+1) + sy(1)*exg(j0, k+1, l+1))*sz(1)
        ex(nn) = ex(nn) + a

        ! Compute Ey on particle
        a = (sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0, l))*sz(0) + (sx(0)*eyg(j, k0,   &
        l+1) + sx(1)*eyg(j+1, k0, l+1))*sz(1)
        ey(nn) = ey(nn) + a

        ! Compute Ez on particle
        a = (sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k, l0))*sy(0) + (sx(0)*ezg(j, k+1,  &
        l0) + sx(1)*ezg(j+1, k+1, l0))*sy(1)
        ez(nn) = ez(nn) + a

        ! Compute Bx on particle
        a = (sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0, l0))
        bx(nn) = bx(nn) + a

        ! Compute By on particle
        a = (sy(0)*byg(j0, k, l0) + sy(1)*byg(j0, k+1, l0))
        by(nn) = by(nn) + a

        ! Compute Bz on particle
        a = (sz(0)*bzg(j0, k0, l) + sz(1)*bzg(j0, k0, l+1))
        bz(nn) = bz(nn) + a

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

  ELSE

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD private(sx, sy, sz, sx0, sy0, sz0)
#endif
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, xp, yp, zp)
      !IBM* ALIGN(64, ex, ey, ez)
      !IBM* ALIGN(64, bx, by, bz)
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      ! !DIR DISTRIBUTE POINT
#endif

      ! Loop over the particles inside a block
      DO nn=ip, MIN(ip-1+lvect, np)

        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi

        ! Compute index of particle
        j=floor(x)
        j0=floor(x-stagger_shift)
        k=floor(y)
        k0=floor(y-stagger_shift)
        l=floor(z)
        l0=floor(z-stagger_shift)
        xint=x-j
        yint=y-k
        zint=z-l

        ! Compute shape factors
        sx(0) = 1.0_num-xint
        sx(1) = xint
        sy(0) = 1.0_num-yint
        sy(1) = yint
        sz(0) = 1.0_num-zint
        sz(1) = zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        sx0(0) = 1.0_num-xint
        sx0(1) = xint
        sy0(0) = 1.0_num-yint
        sy0(1) = yint
        sz0(0) = 1.0_num-zint
        sz0(1) = zint

        ! Compute Ex on particle
        a = (sx0(0)*exg(j0, k, l) + sx0(1)*exg(j0+1, k, l))*sy(0)
        a = a + (sx0(0)*exg(j0, k+1, l) + sx0(1)*exg(j0+1, k+1, l))*sy(1)
        ex(nn) = ex(nn) + a*sz(0)
        a = (sx0(0)*exg(j0, k, l+1) + sx0(1)*exg(j0+1, k, l+1))*sy(0)
        a = a + (sx0(0)*exg(j0, k+1, l+1) + sx0(1)*exg(j0+1, k+1, l+1))*sy(1)
        ex(nn) = ex(nn) + a*sz(1)

        ! Compute Ey on particle
        a = (sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0, l))*sy0(0)
        a = a + (sx(0)*eyg(j, k0+1, l) + sx(1)*eyg(j+1, k0+1, l))*sy0(1)
        ey(nn) = ey(nn) + a*sz(0)
        a = (sx(0)*eyg(j, k0, l+1) + sx(1)*eyg(j+1, k0, l+1))*sy0(0)
        a = a + (sx(0)*eyg(j, k0+1, l+1) + sx(1)*eyg(j+1, k0+1, l+1))*sy0(1)
        ey(nn) = ey(nn) + a*sz(1)

        ! Compute Ez on particle
        a = (sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k, l0))*sy(0)
        a = a + (sx(0)*ezg(j, k+1, l0) + sx(1)*ezg(j+1, k+1, l0))*sy(1)
        ez(nn) = ez(nn) + a*sz0(0)
        a = (sx(0)*ezg(j, k, l0+1) + sx(1)*ezg(j+1, k, l0+1))*sy(0)
        a = a + (sx(0)*ezg(j, k+1, l0+1) + sx(1)*ezg(j+1, k+1, l0+1))*sy(1)
        ez(nn) = ez(nn) + a*sz0(1)

        ! Compute Bx on particle
        a = (sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0, l0))*sy0(0)
        a = a + (sx(0)*bxg(j, k0+1, l0) + sx(1)*bxg(j+1, k0+1, l0))*sy0(1)
        bx(nn) = bx(nn) + a*sz0(0)
        a = (sx(0)*bxg(j, k0, l0+1) + sx(1)*bxg(j+1, k0, l0+1))*sy0(0)
        a = a + (sx(0)*bxg(j, k0+1, l0+1) + sx(1)*bxg(j+1, k0+1, l0+1))*sy0(1)
        bx(nn) = bx(nn) + a*sz0(1)

        ! Compute By on particle
        a = (sx0(0)*byg(j0, k, l0) + sx0(1)*byg(j0+1, k, l0))*sy(0)
        a = a + (sx0(0)*byg(j0, k+1, l0) + sx0(1)*byg(j0+1, k+1, l0))*sy(1)
        by(nn) = by(nn) + a*sz0(0)
        a = (sx0(0)*byg(j0, k, l0+1) + sx0(1)*byg(j0+1, k, l0+1))*sy(0)
        a = a + (sx0(0)*byg(j0, k+1, l0+1) + sx0(1)*byg(j0+1, k+1, l0+1))*sy(1)
        by(nn) = by(nn) + a*sz0(1)

        ! Compute Bz on particle
        a = (sx0(0)*bzg(j0, k0, l) + sx0(1)*bzg(j0+1, k0, l))*sy0(0)
        a = a + (sx0(0)*bzg(j0, k0+1, l) + sx0(1)*bzg(j0+1, k0+1, l))*sy0(1)
        bz(nn) = bz(nn) + a*sz(0)
        a = (sx0(0)*bzg(j0, k0, l+1) + sx0(1)*bzg(j0+1, k0, l+1))*sy0(0)
        a = a + (sx0(0)*bzg(j0, k0+1, l+1) + sx0(1)*bzg(j0+1, k0+1, l+1))*sy0(1)
        bz(nn) = bz(nn) + a*sz(1)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE

#if defined (DEV)
! ________________________________________________________________________________________
!
!> @brief
!> Field gathering (order 1) with gathering of E and B merged in a single loop
!
!> @detail
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vec_1_1_1_v2(np, xp, yp, zp, ex, ey, ez, bx, by, &
  bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg,    &
  ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                                :: l_lower_order_in_v, l_nodal
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, ipmax
  INTEGER(idp)                           :: size
  REAL(num)                              :: dxi, dyi, dzi

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  ! ___ Loop on all the partciles by block of size lvect _______________________
  DO ip=1, np, lvect
    size = MIN(lvect, np-ip+1)
    ipmax = ip+size-1
    ! Field gathering on the block lvect
    CALL geteb3d_energy_conserving_vec_1_1_1_sub(size, xp(ip:ipmax), yp(ip:ipmax),    &
    zp(ip:ipmax), ex(ip:ipmax), ey(ip:ipmax), ez(ip:ipmax), bx(ip:ipmax),             &
    by(ip:ipmax), bz(ip:ipmax), xmin, ymin, zmin, dxi, dyi, dzi, nx, ny, nz, nxguard, &
    nyguard, nzguard, exg, eyg, ezg, bxg, byg, bzg)

  ENDDO

  RETURN
END SUBROUTINE
#endif

#if defined (DEV)
! ________________________________________________________________________________________
!
!> @brief
!> This vectorized subroutine perform the field gathering on a subblock of particles.
!
!> @detail
!> This subroutine is called in geteb3d_energy_conserving_vec_1_1_1().
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] size size of the particle block
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dxi, dyi, dzi inverse of the space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vec_1_1_1_sub(size, xp, yp, zp, ex, ey, ez, bx,  &
  by, bz, xmin, ymin, zmin, dxi, dyi, dzi, nx, ny, nz, nxguard, nyguard, nzguard, exg,  &
  eyg, ezg, bxg, byg, bzg)
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: size
  INTEGER(idp)                           :: nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(size)             :: xp, yp, zp, ex, ey, ez, bx, by, bz
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dxi, dyi, dzi
  INTEGER(isp), DIMENSION(size)          :: j, k, l
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j, k, l
#endif
  INTEGER(isp), DIMENSION(size)         :: j0, k0, l0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j0, k0, l0
#endif
  INTEGER(idp)                           :: jj, kk, ll
  REAL(num)                              :: x, y, z, xint, yint, zint
  INTEGER(isp)                           :: n
  REAL(num), DIMENSION(size, 0:1)         :: sx, sy, sz
  REAL(num), DIMENSION(size, 0:1)         :: sx0, sy0, sz0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: sx, sy, sz, sx0, sy0, sz0
#endif
  REAL(num), PARAMETER                   :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                   :: twothird=2.0_num/3.0_num

  ! ___ Computation _________________________________________________

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
  !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
  !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
  !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
  !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, xp, yp, zp)
  !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
  !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
  !DIR$ IVDEP
  !!DIR DISTRIBUTE POINT
#endif
  DO n=1, size

    x = (xp(n)-xmin)*dxi
    y = (yp(n)-ymin)*dyi
    z = (zp(n)-zmin)*dzi

    ! Compute index of particle
    j(n)=floor(x)
    j0(n)=floor(x+0.5_num)
    k(n)=floor(y)
    k0(n)=floor(y+0.5_num)
    l(n)=floor(z)
    l0(n)=floor(z+0.5_num)

    xint=x-j(n)
    yint=y-k(n)
    zint=z-l(n)

    ! Compute shape factors
    sx(n, 0) = 1.0_num-xint
    sx(n, 1) = xint
    sy(n, 0) = 1.0_num-yint
    sy(n, 1) = yint
    sz(n, 0) = 1.0_num-zint
    sz(n, 1) = zint

    xint=x-j0(n)
    yint=y-k0(n)
    zint=z-l0(n)

    sx0(n, 0) = 1.0_num
    sy0(n, 0) = 1.0_num
    sz0(n, 0) = 1.0_num

  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
  !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
  !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
  !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
  !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, ex, ey, ez)
  !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
  !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
  !DIR$ IVDEP
  !!DIR DISTRIBUTE POINT
#endif
  DO n=1, size

    ! Compute Ex on particle
    ex(n) = ex(n) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0(n), k(n), l(n))
    ex(n) = ex(n) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0(n), k(n)+1, l(n))
    ex(n) = ex(n) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0(n), k(n), l(n)+1)
    ex(n) = ex(n) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0(n), k(n)+1, l(n)+1)

    ! Compute Ey on particle
    ey(n) = ey(n) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j(n), k0(n), l(n))
    ey(n) = ey(n) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j(n)+1, k0(n), l(n))
    ey(n) = ey(n) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j(n), k0(n), l(n)+1)
    ey(n) = ey(n) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j(n)+1, k0(n), l(n)+1)

    ! Compute Ez on particle
    ez(n) = ez(n) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j(n), k(n), l0(n))
    ez(n) = ez(n) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j(n)+1, k(n), l0(n))
    ez(n) = ez(n) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j(n), k(n)+1, l0(n))
    ez(n) = ez(n) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j(n)+1, k(n)+1, l0(n))

  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
  !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
  !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
  !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
  !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, bx, by, bz)
  !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
  !DIR$ SIMD
#endif
#if defined __INTEL_COMPILER
  !DIR$ IVDEP
  !!DIR DISTRIBUTE POINT
#endif
  DO n=1, size

    ! Compute Bx on particle
    bx(n) = bx(n) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j(n), k0(n), l0(n))
    bx(n) = bx(n) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j(n)+1, k0(n), l0(n))

    ! Compute By on particle
    by(n) = by(n) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0(n), k(n), l0(n))
    by(n) = by(n) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0(n), k(n)+1, l0(n))

    ! Compute Bz on particle
    bz(n) = bz(n) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0(n), k0(n), l(n))
    bz(n) = bz(n) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0(n), k0(n), l(n)+1)

  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE
#endif

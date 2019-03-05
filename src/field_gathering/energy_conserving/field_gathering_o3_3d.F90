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
! FIELD_GATHERING_O3_3D.F90
!
! Field gathering subroutines in 3D at order 3
!
! Developers:
! - Henri vincenti
! - Mathieu Lobet
!
! List of subroutines:
!
! - gete3d_energy_conserving_scalar_3_3_3
! - getb3d_energy_conserving_scalar_3_3_3

! - gete3d_energy_conserving_linear_3_3_3
! - getb3d_energy_conserving_linear_3_3_3

! - gete3d_energy_conserving_vec_3_3_3
! - getb3d_energy_conserving_vec_3_3_3

! - gete3d_energy_conserving_vec2_3_3_3
! - getb3d_energy_conserving_vec2_3_3_3

! - geteb3d_energy_conserving_vec_3_3_3

! - geteb3d_energy_conserving_vec2_3_3_3

! - gete3d_energy_conserving_vecblock_3_3_3
! - getb3d_energy_conserving_vecblock_3_3_3

! - gete3d_energy_conserving_vecblock2_3_3_3
! - getb3d_energy_conserving_vecblock2_3_3_3
!
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Scalar version: gathering of electric field from Yee grid ("energy conserving")
!> on particles at order 3
!
!> @details
!> This function is not vectorized but performs better than general order subroutines.
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
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!> exg, eyg, ezg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 3 integers)
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE gete3d_energy_conserving_scalar_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin,    &
  ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
  ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np
  INTEGER(idp), intent(in)             :: exg_nguard(3), exg_nvalid(3),               &
  eyg_nguard(3), eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3)
  REAL(num), DIMENSION(np)             :: xp, yp, zp, ex, ey, ez
  logical(lp)                          :: l_lower_order_in_v, l_nodal
  real(num)                            :: stagger_shift
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
  REAL(num)                            :: xint, yint, zint, xintsq, oxint, yintsq,    &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(-1:2)            :: sx, sx0
  REAL(num), DIMENSION(-1:2)            :: sy, sy0
  REAL(num), DIMENSION(-1:2)            :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1._num/dx
  dyi = 1._num/dy
  dzi = 1._num/dz

  ixmin = -1
  ixmax =  1
  iymin = -1
  iymax =  1
  izmin = -1
  izmax =  1

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  IF (l_lower_order_in_v) THEN

    ixmin0 = -1
    ixmax0 =  1
    iymin0 = -1
    iymax0 =  1
    izmin0 = -1
    izmax0 =  1

    !$acc parallel deviceptr(exg, eyg, ezg, xp, yp, zp, ex, ey, ez)
    !$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2), sz0(-1:2))
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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

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

    ENDDO
    !$acc end loop
    !$acc end parallel

  ELSE

    ixmin0 = -1
    ixmax0 =  2
    iymin0 = -1
    iymax0 =  2
    izmin0 = -1
    izmax0 =  2

    !$acc parallel deviceptr(exg, eyg, ezg, xp, yp, zp, ex, ey, ez)
    !$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2), sz0(-1:2))
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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx0( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy0(-1) = onesixth*oyintsq*oyint
      sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy0( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz0( 2) = onesixth*zintsq*zint

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
    ENDDO! end loop on particles
    !$acc end loop
    !$acc end parallel
  ENDIF
  RETURN
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> Scalar version: Gathering of Magnetic field from Yee grid ("energy conserving")
!>  on particles at order 3
!
!> @details
!> This function is NOT vectorized but performs better than general order subroutines.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of
!> guard cells of the bxg, byg, bzg arrays in each direction
!>(1d arrays containing 3 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 3 integers)
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE getb3d_energy_conserving_scalar_3_3_3(np, xp, yp, zp, bx, by, bz, xmin,    &
  ymin, zmin, dx, dy, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
  bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np
  INTEGER(idp), intent(in)             :: bxg_nguard(3), bxg_nvalid(3),               &
  byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  REAL(num), DIMENSION(np)             :: xp, yp, zp, bx, by, bz
  logical(lp)                          :: l_lower_order_in_v, l_nodal
  real(num)                            :: stagger_shift
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
  REAL(num)                            :: dxi, dyi, dzi, x, y, z, xint, yint, zint,   &
  xintsq, oxint, yintsq, oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(-1:2)           :: sx, sx0
  REAL(num), DIMENSION(-1:2)           :: sy, sy0
  REAL(num), DIMENSION(-1:2)           :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dyi = 1./dy
  dzi = 1./dz

  ixmin = -1
  ixmax =  1
  iymin = -1
  iymax =  1
  izmin = -1
  izmax =  1

  if (l_lower_order_in_v) then

    ixmin0 = -1
    ixmax0 =  1
    iymin0 = -1
    iymax0 =  1
    izmin0 = -1
    izmax0 =  1

    !$acc parallel deviceptr(bxg, byg, bzg, xp, yp, zp, bx, by, bz)
    !$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2), sz0(-1:2))
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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0
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

    ENDDO
    !$acc end loop
    !$acc end parallel

  ELSE

    ixmin0 = -1
    ixmax0 =  2
    iymin0 = -1
    iymax0 =  2
    izmin0 = -1
    izmax0 =  2

    !$acc parallel deviceptr(bxg, byg, bzg, xp, yp, zp, bx, by, bz)
    !$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2), sz0(-1:2))
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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint
      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx0( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy0(-1) = onesixth*oyintsq*oyint
      sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy0( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz0( 2) = onesixth*zintsq*zint

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

    ENDDO
    !$acc end loop
    !$acc end parallel

  ENDIF

  RETURN
END SUBROUTINE

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Scalar version: gathering of electric field from Yee grid ("energy conserving").
!> on particles at order 3
!
!> @details
!> This function is NOT vectorized
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
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE gete3d_energy_conserving_linear_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin,    &
  ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg, ezg,         &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)             :: xp, yp, zp, ex, ey, ez
  logical(lp)                          :: l_lower_order_in_v, l_nodal
  real(num)                            :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: jj, kk, ll, j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint, xintsq, oxint, yintsq,    &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(-1:2)           :: sx, sx0
  REAL(num), DIMENSION(-1:2)           :: sy, sy0
  REAL(num), DIMENSION(-1:2)           :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1._num/dx
  dyi = 1._num/dy
  dzi = 1._num/dz

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  !   write(0, *) 'l_lower_order_in_v ', l_lower_order_in_v
  !   write(0, *) 'sum(xp)', sum(xp), sum(yp), sum(zp)
  !   write(0, *) 'sum(exg)', sum(exg), sum(eyg), sum(ezg)

  IF (l_lower_order_in_v) THEN

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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

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
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(-1)*exg(j0-1, k-1, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(-1)*exg(j0, k-1, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(-1)*exg(j0+1, k-1, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(-1)*exg(j0-1, k, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(-1)*exg(j0, k, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(-1)*exg(j0+1, k, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(-1)*exg(j0-1, k+1, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(-1)*exg(j0, k+1, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(-1)*exg(j0+1, k+1, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(-1)*exg(j0-1, k+2, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(-1)*exg(j0, k+2, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(-1)*exg(j0+1, k+2, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(0)*exg(j0-1, k-1, l)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(0)*exg(j0, k-1, l)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(0)*exg(j0+1, k-1, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(0)*exg(j0-1, k, l)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(0)*exg(j0, k, l)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(0)*exg(j0+1, k, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(0)*exg(j0-1, k+1, l)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(0)*exg(j0, k+1, l)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(0)*exg(j0+1, k+1, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(0)*exg(j0-1, k+2, l)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(0)*exg(j0, k+2, l)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(0)*exg(j0+1, k+2, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(1)*exg(j0-1, k-1, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(1)*exg(j0, k-1, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(1)*exg(j0+1, k-1, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(1)*exg(j0-1, k, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(1)*exg(j0, k, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(1)*exg(j0+1, k, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(1)*exg(j0-1, k+1, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(1)*exg(j0, k+1, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(1)*exg(j0+1, k+1, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(1)*exg(j0-1, k+2, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(1)*exg(j0, k+2, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(1)*exg(j0+1, k+2, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(2)*exg(j0-1, k-1, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(2)*exg(j0, k-1, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(2)*exg(j0+1, k-1, l+2)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(2)*exg(j0-1, k, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(2)*exg(j0, k, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(2)*exg(j0+1, k, l+2)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(2)*exg(j0-1, k+1, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(2)*exg(j0, k+1, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(2)*exg(j0+1, k+1, l+2)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(2)*exg(j0-1, k+2, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(2)*exg(j0, k+2, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(2)*exg(j0+1, k+2, l+2)

      ! Compute Ey on particle
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(-1)*eyg(j-1, k0-1, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(-1)*eyg(j, k0-1, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(-1)*eyg(j+1, k0-1, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(-1)*eyg(j+2, k0-1, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(-1)*eyg(j-1, k0, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(-1)*eyg(j, k0, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(-1)*eyg(j+1, k0, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(-1)*eyg(j+2, k0, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(-1)*eyg(j-1, k0+1, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(-1)*eyg(j, k0+1, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(-1)*eyg(j+1, k0+1, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(-1)*eyg(j+2, k0+1, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(0)*eyg(j-1, k0-1, l)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(0)*eyg(j, k0-1, l)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(0)*eyg(j+1, k0-1, l)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(0)*eyg(j+2, k0-1, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(0)*eyg(j-1, k0, l)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(0)*eyg(j, k0, l)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(0)*eyg(j+1, k0, l)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(0)*eyg(j+2, k0, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(0)*eyg(j-1, k0+1, l)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(0)*eyg(j, k0+1, l)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(0)*eyg(j+1, k0+1, l)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(0)*eyg(j+2, k0+1, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(1)*eyg(j-1, k0-1, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(1)*eyg(j, k0-1, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(1)*eyg(j+1, k0-1, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(1)*eyg(j+2, k0-1, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(1)*eyg(j-1, k0, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(1)*eyg(j, k0, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(1)*eyg(j+1, k0, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(1)*eyg(j+2, k0, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(1)*eyg(j-1, k0+1, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(1)*eyg(j, k0+1, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(1)*eyg(j+1, k0+1, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(1)*eyg(j+2, k0+1, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(2)*eyg(j-1, k0-1, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(2)*eyg(j, k0-1, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(2)*eyg(j+1, k0-1, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(2)*eyg(j+2, k0-1, l+2)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(2)*eyg(j-1, k0, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(2)*eyg(j, k0, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(2)*eyg(j+1, k0, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(2)*eyg(j+2, k0, l+2)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(2)*eyg(j-1, k0+1, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(2)*eyg(j, k0+1, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(2)*eyg(j+1, k0+1, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(2)*eyg(j+2, k0+1, l+2)

      ! Compute Ez on particle
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(-1)*ezg(j-1, k-1, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(-1)*ezg(j, k-1, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(-1)*ezg(j+1, k-1, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(-1)*ezg(j+2, k-1, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(-1)*ezg(j-1, k, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(-1)*ezg(j, k, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(-1)*ezg(j+1, k, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(-1)*ezg(j+2, k, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(-1)*ezg(j-1, k+1, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(-1)*ezg(j, k+1, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(-1)*ezg(j+1, k+1, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(-1)*ezg(j+2, k+1, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(-1)*ezg(j-1, k+2, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(-1)*ezg(j, k+2, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(-1)*ezg(j+1, k+2, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(-1)*ezg(j+2, k+2, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(0)*ezg(j-1, k-1, l0)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(0)*ezg(j, k-1, l0)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(0)*ezg(j+1, k-1, l0)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(0)*ezg(j+2, k-1, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(0)*ezg(j-1, k, l0)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(0)*ezg(j, k, l0)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(0)*ezg(j+1, k, l0)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(0)*ezg(j+2, k, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(0)*ezg(j-1, k+1, l0)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(0)*ezg(j, k+1, l0)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(0)*ezg(j+1, k+1, l0)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(0)*ezg(j+2, k+1, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(0)*ezg(j-1, k+2, l0)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(0)*ezg(j, k+2, l0)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(0)*ezg(j+1, k+2, l0)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(0)*ezg(j+2, k+2, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(1)*ezg(j-1, k-1, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(1)*ezg(j, k-1, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(1)*ezg(j+1, k-1, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(1)*ezg(j+2, k-1, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(1)*ezg(j-1, k, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(1)*ezg(j, k, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(1)*ezg(j+1, k, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(1)*ezg(j+2, k, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(1)*ezg(j-1, k+1, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(1)*ezg(j, k+1, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(1)*ezg(j+1, k+1, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(1)*ezg(j+2, k+1, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(1)*ezg(j-1, k+2, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(1)*ezg(j, k+2, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(1)*ezg(j+1, k+2, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(1)*ezg(j+2, k+2, l0+1)
    ENDDO

  ELSE

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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx0( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy0(-1) = onesixth*oyintsq*oyint
      sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy0( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz0( 2) = onesixth*zintsq*zint


      ! Compute Ex on particle
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(-1)*exg(j0-1, k-1, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(-1)*exg(j0, k-1, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(-1)*exg(j0+1, k-1, l-1)
      ex(ip) = ex(ip) + sx0(2)*sy(-1)*sz(-1)*exg(j0+2, k-1, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(-1)*exg(j0-1, k, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(-1)*exg(j0, k, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(-1)*exg(j0+1, k, l-1)
      ex(ip) = ex(ip) + sx0(2)*sy(0)*sz(-1)*exg(j0+2, k, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(-1)*exg(j0-1, k+1, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(-1)*exg(j0, k+1, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(-1)*exg(j0+1, k+1, l-1)
      ex(ip) = ex(ip) + sx0(2)*sy(1)*sz(-1)*exg(j0+2, k+1, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(-1)*exg(j0-1, k+2, l-1)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(-1)*exg(j0, k+2, l-1)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(-1)*exg(j0+1, k+2, l-1)
      ex(ip) = ex(ip) + sx0(2)*sy(2)*sz(-1)*exg(j0+2, k+2, l-1)
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(0)*exg(j0-1, k-1, l)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(0)*exg(j0, k-1, l)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(0)*exg(j0+1, k-1, l)
      ex(ip) = ex(ip) + sx0(2)*sy(-1)*sz(0)*exg(j0+2, k-1, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(0)*exg(j0-1, k, l)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(0)*exg(j0, k, l)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(0)*exg(j0+1, k, l)
      ex(ip) = ex(ip) + sx0(2)*sy(0)*sz(0)*exg(j0+2, k, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(0)*exg(j0-1, k+1, l)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(0)*exg(j0, k+1, l)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(0)*exg(j0+1, k+1, l)
      ex(ip) = ex(ip) + sx0(2)*sy(1)*sz(0)*exg(j0+2, k+1, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(0)*exg(j0-1, k+2, l)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(0)*exg(j0, k+2, l)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(0)*exg(j0+1, k+2, l)
      ex(ip) = ex(ip) + sx0(2)*sy(2)*sz(0)*exg(j0+2, k+2, l)
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(1)*exg(j0-1, k-1, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(1)*exg(j0, k-1, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(1)*exg(j0+1, k-1, l+1)
      ex(ip) = ex(ip) + sx0(2)*sy(-1)*sz(1)*exg(j0+2, k-1, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(1)*exg(j0-1, k, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(1)*exg(j0, k, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(1)*exg(j0+1, k, l+1)
      ex(ip) = ex(ip) + sx0(2)*sy(0)*sz(1)*exg(j0+2, k, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(1)*exg(j0-1, k+1, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(1)*exg(j0, k+1, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(1)*exg(j0+1, k+1, l+1)
      ex(ip) = ex(ip) + sx0(2)*sy(1)*sz(1)*exg(j0+2, k+1, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(1)*exg(j0-1, k+2, l+1)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(1)*exg(j0, k+2, l+1)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(1)*exg(j0+1, k+2, l+1)
      ex(ip) = ex(ip) + sx0(2)*sy(2)*sz(1)*exg(j0+2, k+2, l+1)
      ex(ip) = ex(ip) + sx0(-1)*sy(-1)*sz(2)*exg(j0-1, k-1, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(-1)*sz(2)*exg(j0, k-1, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(-1)*sz(2)*exg(j0+1, k-1, l+2)
      ex(ip) = ex(ip) + sx0(2)*sy(-1)*sz(2)*exg(j0+2, k-1, l+2)
      ex(ip) = ex(ip) + sx0(-1)*sy(0)*sz(2)*exg(j0-1, k, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(0)*sz(2)*exg(j0, k, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(0)*sz(2)*exg(j0+1, k, l+2)
      ex(ip) = ex(ip) + sx0(2)*sy(0)*sz(2)*exg(j0+2, k, l+2)
      ex(ip) = ex(ip) + sx0(-1)*sy(1)*sz(2)*exg(j0-1, k+1, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(1)*sz(2)*exg(j0, k+1, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(1)*sz(2)*exg(j0+1, k+1, l+2)
      ex(ip) = ex(ip) + sx0(2)*sy(1)*sz(2)*exg(j0+2, k+1, l+2)
      ex(ip) = ex(ip) + sx0(-1)*sy(2)*sz(2)*exg(j0-1, k+2, l+2)
      ex(ip) = ex(ip) + sx0(0)*sy(2)*sz(2)*exg(j0, k+2, l+2)
      ex(ip) = ex(ip) + sx0(1)*sy(2)*sz(2)*exg(j0+1, k+2, l+2)
      ex(ip) = ex(ip) + sx0(2)*sy(2)*sz(2)*exg(j0+2, k+2, l+2)

      ! Compute Ey on particle
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(-1)*eyg(j-1, k0-1, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(-1)*eyg(j, k0-1, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(-1)*eyg(j+1, k0-1, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(-1)*eyg(j+2, k0-1, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(-1)*eyg(j-1, k0, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(-1)*eyg(j, k0, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(-1)*eyg(j+1, k0, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(-1)*eyg(j+2, k0, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(-1)*eyg(j-1, k0+1, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(-1)*eyg(j, k0+1, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(-1)*eyg(j+1, k0+1, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(-1)*eyg(j+2, k0+1, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(2)*sz(-1)*eyg(j-1, k0+2, l-1)
      ey(ip) = ey(ip) + sx(0)*sy0(2)*sz(-1)*eyg(j, k0+2, l-1)
      ey(ip) = ey(ip) + sx(1)*sy0(2)*sz(-1)*eyg(j+1, k0+2, l-1)
      ey(ip) = ey(ip) + sx(2)*sy0(2)*sz(-1)*eyg(j+2, k0+2, l-1)
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(0)*eyg(j-1, k0-1, l)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(0)*eyg(j, k0-1, l)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(0)*eyg(j+1, k0-1, l)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(0)*eyg(j+2, k0-1, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(0)*eyg(j-1, k0, l)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(0)*eyg(j, k0, l)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(0)*eyg(j+1, k0, l)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(0)*eyg(j+2, k0, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(0)*eyg(j-1, k0+1, l)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(0)*eyg(j, k0+1, l)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(0)*eyg(j+1, k0+1, l)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(0)*eyg(j+2, k0+1, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(2)*sz(0)*eyg(j-1, k0+2, l)
      ey(ip) = ey(ip) + sx(0)*sy0(2)*sz(0)*eyg(j, k0+2, l)
      ey(ip) = ey(ip) + sx(1)*sy0(2)*sz(0)*eyg(j+1, k0+2, l)
      ey(ip) = ey(ip) + sx(2)*sy0(2)*sz(0)*eyg(j+2, k0+2, l)
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(1)*eyg(j-1, k0-1, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(1)*eyg(j, k0-1, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(1)*eyg(j+1, k0-1, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(1)*eyg(j+2, k0-1, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(1)*eyg(j-1, k0, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(1)*eyg(j, k0, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(1)*eyg(j+1, k0, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(1)*eyg(j+2, k0, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(1)*eyg(j-1, k0+1, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(1)*eyg(j, k0+1, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(1)*eyg(j+1, k0+1, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(1)*eyg(j+2, k0+1, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(2)*sz(1)*eyg(j-1, k0+2, l+1)
      ey(ip) = ey(ip) + sx(0)*sy0(2)*sz(1)*eyg(j, k0+2, l+1)
      ey(ip) = ey(ip) + sx(1)*sy0(2)*sz(1)*eyg(j+1, k0+2, l+1)
      ey(ip) = ey(ip) + sx(2)*sy0(2)*sz(1)*eyg(j+2, k0+2, l+1)
      ey(ip) = ey(ip) + sx(-1)*sy0(-1)*sz(2)*eyg(j-1, k0-1, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(-1)*sz(2)*eyg(j, k0-1, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(-1)*sz(2)*eyg(j+1, k0-1, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(-1)*sz(2)*eyg(j+2, k0-1, l+2)
      ey(ip) = ey(ip) + sx(-1)*sy0(0)*sz(2)*eyg(j-1, k0, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(0)*sz(2)*eyg(j, k0, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(0)*sz(2)*eyg(j+1, k0, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(0)*sz(2)*eyg(j+2, k0, l+2)
      ey(ip) = ey(ip) + sx(-1)*sy0(1)*sz(2)*eyg(j-1, k0+1, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(1)*sz(2)*eyg(j, k0+1, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(1)*sz(2)*eyg(j+1, k0+1, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(1)*sz(2)*eyg(j+2, k0+1, l+2)
      ey(ip) = ey(ip) + sx(-1)*sy0(2)*sz(2)*eyg(j-1, k0+2, l+2)
      ey(ip) = ey(ip) + sx(0)*sy0(2)*sz(2)*eyg(j, k0+2, l+2)
      ey(ip) = ey(ip) + sx(1)*sy0(2)*sz(2)*eyg(j+1, k0+2, l+2)
      ey(ip) = ey(ip) + sx(2)*sy0(2)*sz(2)*eyg(j+2, k0+2, l+2)

      ! Compute Ez on particle
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(-1)*ezg(j-1, k-1, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(-1)*ezg(j, k-1, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(-1)*ezg(j+1, k-1, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(-1)*ezg(j+2, k-1, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(-1)*ezg(j-1, k, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(-1)*ezg(j, k, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(-1)*ezg(j+1, k, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(-1)*ezg(j+2, k, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(-1)*ezg(j-1, k+1, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(-1)*ezg(j, k+1, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(-1)*ezg(j+1, k+1, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(-1)*ezg(j+2, k+1, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(-1)*ezg(j-1, k+2, l0-1)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(-1)*ezg(j, k+2, l0-1)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(-1)*ezg(j+1, k+2, l0-1)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(-1)*ezg(j+2, k+2, l0-1)
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(0)*ezg(j-1, k-1, l0)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(0)*ezg(j, k-1, l0)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(0)*ezg(j+1, k-1, l0)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(0)*ezg(j+2, k-1, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(0)*ezg(j-1, k, l0)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(0)*ezg(j, k, l0)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(0)*ezg(j+1, k, l0)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(0)*ezg(j+2, k, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(0)*ezg(j-1, k+1, l0)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(0)*ezg(j, k+1, l0)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(0)*ezg(j+1, k+1, l0)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(0)*ezg(j+2, k+1, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(0)*ezg(j-1, k+2, l0)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(0)*ezg(j, k+2, l0)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(0)*ezg(j+1, k+2, l0)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(0)*ezg(j+2, k+2, l0)
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(1)*ezg(j-1, k-1, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(1)*ezg(j, k-1, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(1)*ezg(j+1, k-1, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(1)*ezg(j+2, k-1, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(1)*ezg(j-1, k, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(1)*ezg(j, k, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(1)*ezg(j+1, k, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(1)*ezg(j+2, k, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(1)*ezg(j-1, k+1, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(1)*ezg(j, k+1, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(1)*ezg(j+1, k+1, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(1)*ezg(j+2, k+1, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(1)*ezg(j-1, k+2, l0+1)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(1)*ezg(j, k+2, l0+1)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(1)*ezg(j+1, k+2, l0+1)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(1)*ezg(j+2, k+2, l0+1)
      ez(ip) = ez(ip) + sx(-1)*sy(-1)*sz0(2)*ezg(j-1, k-1, l0+2)
      ez(ip) = ez(ip) + sx(0)*sy(-1)*sz0(2)*ezg(j, k-1, l0+2)
      ez(ip) = ez(ip) + sx(1)*sy(-1)*sz0(2)*ezg(j+1, k-1, l0+2)
      ez(ip) = ez(ip) + sx(2)*sy(-1)*sz0(2)*ezg(j+2, k-1, l0+2)
      ez(ip) = ez(ip) + sx(-1)*sy(0)*sz0(2)*ezg(j-1, k, l0+2)
      ez(ip) = ez(ip) + sx(0)*sy(0)*sz0(2)*ezg(j, k, l0+2)
      ez(ip) = ez(ip) + sx(1)*sy(0)*sz0(2)*ezg(j+1, k, l0+2)
      ez(ip) = ez(ip) + sx(2)*sy(0)*sz0(2)*ezg(j+2, k, l0+2)
      ez(ip) = ez(ip) + sx(-1)*sy(1)*sz0(2)*ezg(j-1, k+1, l0+2)
      ez(ip) = ez(ip) + sx(0)*sy(1)*sz0(2)*ezg(j, k+1, l0+2)
      ez(ip) = ez(ip) + sx(1)*sy(1)*sz0(2)*ezg(j+1, k+1, l0+2)
      ez(ip) = ez(ip) + sx(2)*sy(1)*sz0(2)*ezg(j+2, k+1, l0+2)
      ez(ip) = ez(ip) + sx(-1)*sy(2)*sz0(2)*ezg(j-1, k+2, l0+2)
      ez(ip) = ez(ip) + sx(0)*sy(2)*sz0(2)*ezg(j, k+2, l0+2)
      ez(ip) = ez(ip) + sx(1)*sy(2)*sz0(2)*ezg(j+1, k+2, l0+2)
      ez(ip) = ez(ip) + sx(2)*sy(2)*sz0(2)*ezg(j+2, k+2, l0+2)

    ENDDO! end loop on particles
  ENDIF

  RETURN
END SUBROUTINE
#endif

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
! Scalar version: Gathering of Magnetic field from Yee grid ("energy conserving")
! on particles
! at order 3
!
!> @details
!> This function is not vectorized
!
!> @date
!> Creation 2016
!
!> @author
!> Mathieu Lobet
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE getb3d_energy_conserving_linear_3_3_3(np, xp, yp, zp, bx, by, bz, xmin,    &
  ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, bxg, byg, bzg,         &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)             :: xp, yp, zp, bx, by, bz
  logical(lp)                          :: l_lower_order_in_v, l_nodal
  real(num)                            :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(idp)                         :: ip, j, k, l, jj, kk, ll, j0, k0, l0
  REAL(num)                            :: dxi, dyi, dzi, x, y, z, xint, yint, zint,   &
  xintsq, oxint, yintsq, oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(-1:2)           :: sx, sx0
  REAL(num), DIMENSION(-1:2)           :: sy, sy0
  REAL(num), DIMENSION(-1:2)           :: sz, sz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                 :: twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dyi = 1./dy
  dzi = 1./dz

  if (l_lower_order_in_v) then

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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0
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
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(-1)*bxg(j-1, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(-1)*bxg(j, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(-1)*bxg(j+1, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(-1)*bxg(j+2, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(-1)*bxg(j-1, k0, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(-1)*bxg(j, k0, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(-1)*bxg(j+1, k0, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(-1)*bxg(j+2, k0, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(-1)*bxg(j-1, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(-1)*bxg(j, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(-1)*bxg(j+1, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(-1)*bxg(j+2, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(0)*bxg(j-1, k0-1, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(0)*bxg(j, k0-1, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(0)*bxg(j+1, k0-1, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(0)*bxg(j+2, k0-1, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(0)*bxg(j-1, k0, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j, k0, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1, k0, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(0)*bxg(j+2, k0, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(0)*bxg(j-1, k0+1, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j, k0+1, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1, k0+1, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(0)*bxg(j+2, k0+1, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(1)*bxg(j-1, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(1)*bxg(j, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(1)*bxg(j+1, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(1)*bxg(j+2, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(1)*bxg(j-1, k0, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j, k0, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1, k0, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(1)*bxg(j+2, k0, l0+1)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(1)*bxg(j-1, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(1)*bxg(j+2, k0+1, l0+1)

      ! Compute By on particle
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(-1)*byg(j0-1, k-1, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(-1)*byg(j0, k-1, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(-1)*byg(j0+1, k-1, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(-1)*byg(j0-1, k, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(-1)*byg(j0, k, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(-1)*byg(j0+1, k, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(-1)*byg(j0-1, k+1, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(-1)*byg(j0, k+1, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(-1)*byg(j0+1, k+1, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(-1)*byg(j0-1, k+2, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(-1)*byg(j0, k+2, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(-1)*byg(j0+1, k+2, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(0)*byg(j0-1, k-1, l0)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(0)*byg(j0, k-1, l0)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(0)*byg(j0+1, k-1, l0)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(0)*byg(j0-1, k, l0)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0, k, l0)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1, k, l0)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(0)*byg(j0-1, k+1, l0)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0, k+1, l0)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1, k+1, l0)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(0)*byg(j0-1, k+2, l0)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(0)*byg(j0, k+2, l0)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(0)*byg(j0+1, k+2, l0)
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(1)*byg(j0-1, k-1, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(1)*byg(j0, k-1, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(1)*byg(j0+1, k-1, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(1)*byg(j0-1, k, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0, k, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1, k, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(1)*byg(j0-1, k+1, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0, k+1, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1, k+1, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(1)*byg(j0-1, k+2, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(1)*byg(j0, k+2, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(1)*byg(j0+1, k+2, l0+1)

      ! Compute Bz on particle
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(-1)*bzg(j0-1, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(-1)*bzg(j0, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(-1)*bzg(j0+1, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(-1)*bzg(j0-1, k0, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(-1)*bzg(j0, k0, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(-1)*bzg(j0+1, k0, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(-1)*bzg(j0-1, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(-1)*bzg(j0, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(-1)*bzg(j0+1, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(0)*bzg(j0-1, k0-1, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(0)*bzg(j0, k0-1, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(0)*bzg(j0+1, k0-1, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(0)*bzg(j0-1, k0, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0, k0, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1, k0, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(0)*bzg(j0-1, k0+1, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(0)*bzg(j0, k0+1, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1, k0+1, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(1)*bzg(j0-1, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(1)*bzg(j0, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(1)*bzg(j0+1, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(1)*bzg(j0-1, k0, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0, k0, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1, k0, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(1)*bzg(j0-1, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(1)*bzg(j0, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(2)*bzg(j0-1, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(2)*bzg(j0, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(2)*bzg(j0+1, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(2)*bzg(j0-1, k0, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(2)*bzg(j0, k0, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(2)*bzg(j0+1, k0, l+2)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(2)*bzg(j0-1, k0+1, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(2)*bzg(j0, k0+1, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(2)*bzg(j0+1, k0+1, l+2)

    ENDDO

  ELSE

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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint
      xint=x-stagger_shift-j0
      yint=y-stagger_shift-k0
      zint=z-stagger_shift-l0

      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx0( 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy0(-1) = onesixth*oyintsq*oyint
      sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy0( 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz0( 2) = onesixth*zintsq*zint

      ! Compute Bx on particle
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(-1)*bxg(j-1, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(-1)*bxg(j, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(-1)*bxg(j+1, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(-1)*bxg(j+2, k0-1, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(-1)*bxg(j-1, k0, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(-1)*bxg(j, k0, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(-1)*bxg(j+1, k0, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(-1)*bxg(j+2, k0, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(-1)*bxg(j-1, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(-1)*bxg(j, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(-1)*bxg(j+1, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(-1)*bxg(j+2, k0+1, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(2)*sz0(-1)*bxg(j-1, k0+2, l0-1)
      bx(ip) = bx(ip) + sx(0)*sy0(2)*sz0(-1)*bxg(j, k0+2, l0-1)
      bx(ip) = bx(ip) + sx(1)*sy0(2)*sz0(-1)*bxg(j+1, k0+2, l0-1)
      bx(ip) = bx(ip) + sx(2)*sy0(2)*sz0(-1)*bxg(j+2, k0+2, l0-1)
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(0)*bxg(j-1, k0-1, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(0)*bxg(j, k0-1, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(0)*bxg(j+1, k0-1, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(0)*bxg(j+2, k0-1, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(0)*bxg(j-1, k0, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j, k0, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1, k0, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(0)*bxg(j+2, k0, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(0)*bxg(j-1, k0+1, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j, k0+1, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1, k0+1, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(0)*bxg(j+2, k0+1, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(2)*sz0(0)*bxg(j-1, k0+2, l0)
      bx(ip) = bx(ip) + sx(0)*sy0(2)*sz0(0)*bxg(j, k0+2, l0)
      bx(ip) = bx(ip) + sx(1)*sy0(2)*sz0(0)*bxg(j+1, k0+2, l0)
      bx(ip) = bx(ip) + sx(2)*sy0(2)*sz0(0)*bxg(j+2, k0+2, l0)
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(1)*bxg(j-1, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(1)*bxg(j, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(1)*bxg(j+1, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(1)*bxg(j+2, k0-1, l0+1)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(1)*bxg(j-1, k0, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j, k0, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1, k0, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(1)*bxg(j+2, k0, l0+1)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(1)*bxg(j-1, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(1)*bxg(j+2, k0+1, l0+1)
      bx(ip) = bx(ip) + sx(-1)*sy0(2)*sz0(1)*bxg(j-1, k0+2, l0+1)
      bx(ip) = bx(ip) + sx(0)*sy0(2)*sz0(1)*bxg(j, k0+2, l0+1)
      bx(ip) = bx(ip) + sx(1)*sy0(2)*sz0(1)*bxg(j+1, k0+2, l0+1)
      bx(ip) = bx(ip) + sx(2)*sy0(2)*sz0(1)*bxg(j+2, k0+2, l0+1)
      bx(ip) = bx(ip) + sx(-1)*sy0(-1)*sz0(2)*bxg(j-1, k0-1, l0+2)
      bx(ip) = bx(ip) + sx(0)*sy0(-1)*sz0(2)*bxg(j, k0-1, l0+2)
      bx(ip) = bx(ip) + sx(1)*sy0(-1)*sz0(2)*bxg(j+1, k0-1, l0+2)
      bx(ip) = bx(ip) + sx(2)*sy0(-1)*sz0(2)*bxg(j+2, k0-1, l0+2)
      bx(ip) = bx(ip) + sx(-1)*sy0(0)*sz0(2)*bxg(j-1, k0, l0+2)
      bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(2)*bxg(j, k0, l0+2)
      bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(2)*bxg(j+1, k0, l0+2)
      bx(ip) = bx(ip) + sx(2)*sy0(0)*sz0(2)*bxg(j+2, k0, l0+2)
      bx(ip) = bx(ip) + sx(-1)*sy0(1)*sz0(2)*bxg(j-1, k0+1, l0+2)
      bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(2)*bxg(j, k0+1, l0+2)
      bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(2)*bxg(j+1, k0+1, l0+2)
      bx(ip) = bx(ip) + sx(2)*sy0(1)*sz0(2)*bxg(j+2, k0+1, l0+2)
      bx(ip) = bx(ip) + sx(-1)*sy0(2)*sz0(2)*bxg(j-1, k0+2, l0+2)
      bx(ip) = bx(ip) + sx(0)*sy0(2)*sz0(2)*bxg(j, k0+2, l0+2)
      bx(ip) = bx(ip) + sx(1)*sy0(2)*sz0(2)*bxg(j+1, k0+2, l0+2)
      bx(ip) = bx(ip) + sx(2)*sy0(2)*sz0(2)*bxg(j+2, k0+2, l0+2)

      ! Compute By on particle
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(-1)*byg(j0-1, k-1, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(-1)*byg(j0, k-1, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(-1)*byg(j0+1, k-1, l0-1)
      by(ip) = by(ip) + sx0(2)*sy(-1)*sz0(-1)*byg(j0+2, k-1, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(-1)*byg(j0-1, k, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(-1)*byg(j0, k, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(-1)*byg(j0+1, k, l0-1)
      by(ip) = by(ip) + sx0(2)*sy(0)*sz0(-1)*byg(j0+2, k, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(-1)*byg(j0-1, k+1, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(-1)*byg(j0, k+1, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(-1)*byg(j0+1, k+1, l0-1)
      by(ip) = by(ip) + sx0(2)*sy(1)*sz0(-1)*byg(j0+2, k+1, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(-1)*byg(j0-1, k+2, l0-1)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(-1)*byg(j0, k+2, l0-1)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(-1)*byg(j0+1, k+2, l0-1)
      by(ip) = by(ip) + sx0(2)*sy(2)*sz0(-1)*byg(j0+2, k+2, l0-1)
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(0)*byg(j0-1, k-1, l0)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(0)*byg(j0, k-1, l0)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(0)*byg(j0+1, k-1, l0)
      by(ip) = by(ip) + sx0(2)*sy(-1)*sz0(0)*byg(j0+2, k-1, l0)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(0)*byg(j0-1, k, l0)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0, k, l0)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1, k, l0)
      by(ip) = by(ip) + sx0(2)*sy(0)*sz0(0)*byg(j0+2, k, l0)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(0)*byg(j0-1, k+1, l0)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0, k+1, l0)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1, k+1, l0)
      by(ip) = by(ip) + sx0(2)*sy(1)*sz0(0)*byg(j0+2, k+1, l0)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(0)*byg(j0-1, k+2, l0)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(0)*byg(j0, k+2, l0)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(0)*byg(j0+1, k+2, l0)
      by(ip) = by(ip) + sx0(2)*sy(2)*sz0(0)*byg(j0+2, k+2, l0)
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(1)*byg(j0-1, k-1, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(1)*byg(j0, k-1, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(1)*byg(j0+1, k-1, l0+1)
      by(ip) = by(ip) + sx0(2)*sy(-1)*sz0(1)*byg(j0+2, k-1, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(1)*byg(j0-1, k, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0, k, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1, k, l0+1)
      by(ip) = by(ip) + sx0(2)*sy(0)*sz0(1)*byg(j0+2, k, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(1)*byg(j0-1, k+1, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0, k+1, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1, k+1, l0+1)
      by(ip) = by(ip) + sx0(2)*sy(1)*sz0(1)*byg(j0+2, k+1, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(1)*byg(j0-1, k+2, l0+1)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(1)*byg(j0, k+2, l0+1)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(1)*byg(j0+1, k+2, l0+1)
      by(ip) = by(ip) + sx0(2)*sy(2)*sz0(1)*byg(j0+2, k+2, l0+1)
      by(ip) = by(ip) + sx0(-1)*sy(-1)*sz0(2)*byg(j0-1, k-1, l0+2)
      by(ip) = by(ip) + sx0(0)*sy(-1)*sz0(2)*byg(j0, k-1, l0+2)
      by(ip) = by(ip) + sx0(1)*sy(-1)*sz0(2)*byg(j0+1, k-1, l0+2)
      by(ip) = by(ip) + sx0(2)*sy(-1)*sz0(2)*byg(j0+2, k-1, l0+2)
      by(ip) = by(ip) + sx0(-1)*sy(0)*sz0(2)*byg(j0-1, k, l0+2)
      by(ip) = by(ip) + sx0(0)*sy(0)*sz0(2)*byg(j0, k, l0+2)
      by(ip) = by(ip) + sx0(1)*sy(0)*sz0(2)*byg(j0+1, k, l0+2)
      by(ip) = by(ip) + sx0(2)*sy(0)*sz0(2)*byg(j0+2, k, l0+2)
      by(ip) = by(ip) + sx0(-1)*sy(1)*sz0(2)*byg(j0-1, k+1, l0+2)
      by(ip) = by(ip) + sx0(0)*sy(1)*sz0(2)*byg(j0, k+1, l0+2)
      by(ip) = by(ip) + sx0(1)*sy(1)*sz0(2)*byg(j0+1, k+1, l0+2)
      by(ip) = by(ip) + sx0(2)*sy(1)*sz0(2)*byg(j0+2, k+1, l0+2)
      by(ip) = by(ip) + sx0(-1)*sy(2)*sz0(2)*byg(j0-1, k+2, l0+2)
      by(ip) = by(ip) + sx0(0)*sy(2)*sz0(2)*byg(j0, k+2, l0+2)
      by(ip) = by(ip) + sx0(1)*sy(2)*sz0(2)*byg(j0+1, k+2, l0+2)
      by(ip) = by(ip) + sx0(2)*sy(2)*sz0(2)*byg(j0+2, k+2, l0+2)

      ! Compute Bz on particle
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(-1)*bzg(j0-1, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(-1)*bzg(j0, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(-1)*bzg(j0+1, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(2)*sy0(-1)*sz(-1)*bzg(j0+2, k0-1, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(-1)*bzg(j0-1, k0, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(-1)*bzg(j0, k0, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(-1)*bzg(j0+1, k0, l-1)
      bz(ip) = bz(ip) + sx0(2)*sy0(0)*sz(-1)*bzg(j0+2, k0, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(-1)*bzg(j0-1, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(-1)*bzg(j0, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(-1)*bzg(j0+1, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(2)*sy0(1)*sz(-1)*bzg(j0+2, k0+1, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(2)*sz(-1)*bzg(j0-1, k0+2, l-1)
      bz(ip) = bz(ip) + sx0(0)*sy0(2)*sz(-1)*bzg(j0, k0+2, l-1)
      bz(ip) = bz(ip) + sx0(1)*sy0(2)*sz(-1)*bzg(j0+1, k0+2, l-1)
      bz(ip) = bz(ip) + sx0(2)*sy0(2)*sz(-1)*bzg(j0+2, k0+2, l-1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(0)*bzg(j0-1, k0-1, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(0)*bzg(j0, k0-1, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(0)*bzg(j0+1, k0-1, l)
      bz(ip) = bz(ip) + sx0(2)*sy0(-1)*sz(0)*bzg(j0+2, k0-1, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(0)*bzg(j0-1, k0, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0, k0, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1, k0, l)
      bz(ip) = bz(ip) + sx0(2)*sy0(0)*sz(0)*bzg(j0+2, k0, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(0)*bzg(j0-1, k0+1, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(0)*bzg(j0, k0+1, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1, k0+1, l)
      bz(ip) = bz(ip) + sx0(2)*sy0(1)*sz(0)*bzg(j0+2, k0+1, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(2)*sz(0)*bzg(j0-1, k0+2, l)
      bz(ip) = bz(ip) + sx0(0)*sy0(2)*sz(0)*bzg(j0, k0+2, l)
      bz(ip) = bz(ip) + sx0(1)*sy0(2)*sz(0)*bzg(j0+1, k0+2, l)
      bz(ip) = bz(ip) + sx0(2)*sy0(2)*sz(0)*bzg(j0+2, k0+2, l)
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(1)*bzg(j0-1, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(1)*bzg(j0, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(1)*bzg(j0+1, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(2)*sy0(-1)*sz(1)*bzg(j0+2, k0-1, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(1)*bzg(j0-1, k0, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0, k0, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1, k0, l+1)
      bz(ip) = bz(ip) + sx0(2)*sy0(0)*sz(1)*bzg(j0+2, k0, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(1)*bzg(j0-1, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(1)*bzg(j0, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(2)*sy0(1)*sz(1)*bzg(j0+2, k0+1, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(2)*sz(1)*bzg(j0-1, k0+2, l+1)
      bz(ip) = bz(ip) + sx0(0)*sy0(2)*sz(1)*bzg(j0, k0+2, l+1)
      bz(ip) = bz(ip) + sx0(1)*sy0(2)*sz(1)*bzg(j0+1, k0+2, l+1)
      bz(ip) = bz(ip) + sx0(2)*sy0(2)*sz(1)*bzg(j0+2, k0+2, l+1)
      bz(ip) = bz(ip) + sx0(-1)*sy0(-1)*sz(2)*bzg(j0-1, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(-1)*sz(2)*bzg(j0, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(-1)*sz(2)*bzg(j0+1, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(2)*sy0(-1)*sz(2)*bzg(j0+2, k0-1, l+2)
      bz(ip) = bz(ip) + sx0(-1)*sy0(0)*sz(2)*bzg(j0-1, k0, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(2)*bzg(j0, k0, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(2)*bzg(j0+1, k0, l+2)
      bz(ip) = bz(ip) + sx0(2)*sy0(0)*sz(2)*bzg(j0+2, k0, l+2)
      bz(ip) = bz(ip) + sx0(-1)*sy0(1)*sz(2)*bzg(j0-1, k0+1, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(2)*bzg(j0, k0+1, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(2)*bzg(j0+1, k0+1, l+2)
      bz(ip) = bz(ip) + sx0(2)*sy0(1)*sz(2)*bzg(j0+2, k0+1, l+2)
      bz(ip) = bz(ip) + sx0(-1)*sy0(2)*sz(2)*bzg(j0-1, k0+2, l+2)
      bz(ip) = bz(ip) + sx0(0)*sy0(2)*sz(2)*bzg(j0, k0+2, l+2)
      bz(ip) = bz(ip) + sx0(1)*sy0(2)*sz(2)*bzg(j0+1, k0+2, l+2)
      bz(ip) = bz(ip) + sx0(2)*sy0(2)*sz(2)*bzg(j0+2, k0+2, l+2)
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE
#endif

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Vectorized Gathering of electric field from Yee grid ("energy conserving") on particles
!> at order 3.
!
!> @detail
!> This function is vectorized
!
!> @date
!> Creation 2016
!
!> @author
!> Mathieu Lobet
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE gete3d_energy_conserving_vec_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin, ymin, &
  zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg, ezg, lvect,        &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE

  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                           :: lvect
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, j, k, l
  INTEGER(isp)                           :: jj, kk, ll
  INTEGER(isp)                           :: j0, k0, l0
  INTEGER(isp)                           :: nn, n
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq,     &
  ozint
  REAL(num)                              :: oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(lvect, -1:2)       :: sx, sx0
  REAL(num), DIMENSION(lvect, -1:2)       :: sy, sy0
  REAL(num), DIMENSION(lvect, -1:2)       :: sz, sz0
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

  IF (l_lower_order_in_v ) THEN

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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        yintsq = yint*yint
        sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, -1)*exg(j0-1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, -1)*exg(j0, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, -1)*exg(j0+1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, -1)*exg(j0-1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, -1)*exg(j0, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, -1)*exg(j0+1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, -1)*exg(j0-1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, -1)*exg(j0, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, -1)*exg(j0+1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, -1)*exg(j0-1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, -1)*exg(j0, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, -1)*exg(j0+1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 0)*exg(j0-1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 0)*exg(j0, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 0)*exg(j0+1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 0)*exg(j0-1, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0+1, k, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 0)*exg(j0-1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0+1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 0)*exg(j0-1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 0)*exg(j0, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 0)*exg(j0+1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 1)*exg(j0-1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 1)*exg(j0, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 1)*exg(j0+1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 1)*exg(j0-1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0+1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 1)*exg(j0-1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0+1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 1)*exg(j0-1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 1)*exg(j0, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 1)*exg(j0+1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 2)*exg(j0-1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 2)*exg(j0, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 2)*exg(j0+1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 2)*exg(j0-1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 2)*exg(j0, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 2)*exg(j0+1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 2)*exg(j0-1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 2)*exg(j0, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 2)*exg(j0+1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 2)*exg(j0-1, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 2)*exg(j0, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 2)*exg(j0+1, k+2, l+2)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, -1)*eyg(j-1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, -1)*eyg(j, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, -1)*eyg(j+1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, -1)*eyg(j+2, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, -1)*eyg(j-1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, -1)*eyg(j, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, -1)*eyg(j+1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, -1)*eyg(j+2, k0, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, -1)*eyg(j-1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, -1)*eyg(j, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, -1)*eyg(j+1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, -1)*eyg(j+2, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 0)*eyg(j-1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 0)*eyg(j, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 0)*eyg(j+1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 0)*eyg(j+2, k0-1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 0)*eyg(j-1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 0)*eyg(j+2, k0, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 0)*eyg(j-1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j+1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 0)*eyg(j+2, k0+1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 1)*eyg(j-1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 1)*eyg(j, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 1)*eyg(j+1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 1)*eyg(j+2, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 1)*eyg(j-1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 1)*eyg(j+2, k0, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 1)*eyg(j-1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j+1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 1)*eyg(j+2, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 2)*eyg(j-1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 2)*eyg(j, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 2)*eyg(j+1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 2)*eyg(j+2, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 2)*eyg(j-1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 2)*eyg(j, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 2)*eyg(j+1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 2)*eyg(j+2, k0, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 2)*eyg(j-1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 2)*eyg(j, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 2)*eyg(j+1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 2)*eyg(j+2, k0+1, l+2)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, -1)*ezg(j-1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, -1)*ezg(j, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, -1)*ezg(j+1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, -1)*ezg(j+2, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, -1)*ezg(j-1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, -1)*ezg(j, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, -1)*ezg(j+1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, -1)*ezg(j+2, k, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, -1)*ezg(j-1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, -1)*ezg(j, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, -1)*ezg(j+1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, -1)*ezg(j+2, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, -1)*ezg(j-1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, -1)*ezg(j, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, -1)*ezg(j+1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, -1)*ezg(j+2, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 0)*ezg(j-1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 0)*ezg(j, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 0)*ezg(j+1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 0)*ezg(j+2, k-1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 0)*ezg(j-1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 0)*ezg(j+2, k, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 0)*ezg(j-1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 0)*ezg(j+2, k+1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 0)*ezg(j-1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 0)*ezg(j, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 0)*ezg(j+1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 0)*ezg(j+2, k+2, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 1)*ezg(j-1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 1)*ezg(j, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 1)*ezg(j+1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 1)*ezg(j+2, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 1)*ezg(j-1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j+1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 1)*ezg(j+2, k, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 1)*ezg(j-1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j+1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 1)*ezg(j+2, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 1)*ezg(j-1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 1)*ezg(j, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 1)*ezg(j+1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 1)*ezg(j+2, k+2, l0+1)
      ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif

    ENDDO

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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(n, -1) = onesixth*oyintsq*oyint
        sy0(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, -1)*exg(j0-1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, -1)*exg(j0, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, -1)*exg(j0+1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, -1)*exg(j0+2, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, -1)*exg(j0-1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, -1)*exg(j0, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, -1)*exg(j0+1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, -1)*exg(j0+2, k, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, -1)*exg(j0-1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, -1)*exg(j0, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, -1)*exg(j0+1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, -1)*exg(j0+2, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, -1)*exg(j0-1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, -1)*exg(j0, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, -1)*exg(j0+1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, -1)*exg(j0+2, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 0)*exg(j0-1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 0)*exg(j0, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 0)*exg(j0+1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, 0)*exg(j0+2, k-1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 0)*exg(j0-1, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0+1, k, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, 0)*exg(j0+2, k, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 0)*exg(j0-1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0+1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, 0)*exg(j0+2, k+1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 0)*exg(j0-1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 0)*exg(j0, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 0)*exg(j0+1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, 0)*exg(j0+2, k+2, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 1)*exg(j0-1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 1)*exg(j0, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 1)*exg(j0+1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, 1)*exg(j0+2, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 1)*exg(j0-1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0+1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, 1)*exg(j0+2, k, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 1)*exg(j0-1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0+1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, 1)*exg(j0+2, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 1)*exg(j0-1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 1)*exg(j0, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 1)*exg(j0+1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, 1)*exg(j0+2, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 2)*exg(j0-1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 2)*exg(j0, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 2)*exg(j0+1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, 2)*exg(j0+2, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 2)*exg(j0-1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 2)*exg(j0, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 2)*exg(j0+1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, 2)*exg(j0+2, k, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 2)*exg(j0-1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 2)*exg(j0, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 2)*exg(j0+1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, 2)*exg(j0+2, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 2)*exg(j0-1, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 2)*exg(j0, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 2)*exg(j0+1, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, 2)*exg(j0+2, k+2, l+2)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, -1)*eyg(j-1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, -1)*eyg(j, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, -1)*eyg(j+1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, -1)*eyg(j+2, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, -1)*eyg(j-1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, -1)*eyg(j, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, -1)*eyg(j+1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, -1)*eyg(j+2, k0, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, -1)*eyg(j-1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, -1)*eyg(j, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, -1)*eyg(j+1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, -1)*eyg(j+2, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, -1)*eyg(j-1, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, -1)*eyg(j, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, -1)*eyg(j+1, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, -1)*eyg(j+2, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 0)*eyg(j-1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 0)*eyg(j, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 0)*eyg(j+1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 0)*eyg(j+2, k0-1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 0)*eyg(j-1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 0)*eyg(j+2, k0, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 0)*eyg(j-1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j+1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 0)*eyg(j+2, k0+1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, 0)*eyg(j-1, k0+2, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, 0)*eyg(j, k0+2, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, 0)*eyg(j+1, k0+2, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, 0)*eyg(j+2, k0+2, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 1)*eyg(j-1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 1)*eyg(j, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 1)*eyg(j+1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 1)*eyg(j+2, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 1)*eyg(j-1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 1)*eyg(j+2, k0, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 1)*eyg(j-1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j+1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 1)*eyg(j+2, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, 1)*eyg(j-1, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, 1)*eyg(j, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, 1)*eyg(j+1, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, 1)*eyg(j+2, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 2)*eyg(j-1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 2)*eyg(j, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 2)*eyg(j+1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 2)*eyg(j+2, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 2)*eyg(j-1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 2)*eyg(j, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 2)*eyg(j+1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 2)*eyg(j+2, k0, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 2)*eyg(j-1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 2)*eyg(j, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 2)*eyg(j+1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 2)*eyg(j+2, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, 2)*eyg(j-1, k0+2, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, 2)*eyg(j, k0+2, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, 2)*eyg(j+1, k0+2, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, 2)*eyg(j+2, k0+2, l+2)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, -1)*ezg(j-1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, -1)*ezg(j, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, -1)*ezg(j+1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, -1)*ezg(j+2, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, -1)*ezg(j-1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, -1)*ezg(j, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, -1)*ezg(j+1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, -1)*ezg(j+2, k, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, -1)*ezg(j-1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, -1)*ezg(j, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, -1)*ezg(j+1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, -1)*ezg(j+2, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, -1)*ezg(j-1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, -1)*ezg(j, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, -1)*ezg(j+1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, -1)*ezg(j+2, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 0)*ezg(j-1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 0)*ezg(j, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 0)*ezg(j+1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 0)*ezg(j+2, k-1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 0)*ezg(j-1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 0)*ezg(j+2, k, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 0)*ezg(j-1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 0)*ezg(j+2, k+1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 0)*ezg(j-1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 0)*ezg(j, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 0)*ezg(j+1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 0)*ezg(j+2, k+2, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 1)*ezg(j-1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 1)*ezg(j, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 1)*ezg(j+1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 1)*ezg(j+2, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 1)*ezg(j-1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j+1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 1)*ezg(j+2, k, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 1)*ezg(j-1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j+1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 1)*ezg(j+2, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 1)*ezg(j-1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 1)*ezg(j, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 1)*ezg(j+1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 1)*ezg(j+2, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 2)*ezg(j-1, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 2)*ezg(j, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 2)*ezg(j+1, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 2)*ezg(j+2, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 2)*ezg(j-1, k, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 2)*ezg(j, k, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 2)*ezg(j+1, k, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 2)*ezg(j+2, k, l0+2)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 2)*ezg(j-1, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 2)*ezg(j, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 2)*ezg(j+1, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 2)*ezg(j+2, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 2)*ezg(j-1, k+2, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 2)*ezg(j, k+2, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 2)*ezg(j+1, k+2, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 2)*ezg(j+2, k+2, l0+2)
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
!> Vectorized gathering of magnetic field from Yee grid ("energy conserving") on particles
!> at order 3
!
!> @date
!> Creation 2016
!
!> @author
!> Mathieu Lobet
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE getb3d_energy_conserving_vec_3_3_3(np, xp, yp, zp, bx, by, bz, xmin, ymin, &
  zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, bxg, byg, bzg, lvect,        &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp, yp, zp, bx, by, bz
  LOGICAL(lp)                          :: l_lower_order_in_v, l_nodal
  REAL(num)                            :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  INTEGER(isp)                         :: jj, kk, ll
  INTEGER(isp)                         :: n, nn
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint, xintsq, oxint, yintsq,    &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq
  REAL(num), DIMENSION(lvect, -1:2)     :: sx, sx0
  REAL(num), DIMENSION(lvect, -1:2)     :: sy, sy0
  REAL(num), DIMENSION(lvect, -1:2)     :: sz, sz0
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

    ! Loop over the particles by block
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0
        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2
        yintsq = yint*yint
        sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2
        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, -1)*bxg(j-1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, -1)*bxg(j, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, -1)*bxg(j+1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, -1)*bxg(j+2, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, -1)*bxg(j-1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, -1)*bxg(j, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, -1)*bxg(j+1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, -1)*bxg(j+2, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, -1)*bxg(j-1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, -1)*bxg(j, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, -1)*bxg(j+1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, -1)*bxg(j+2, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 0)*bxg(j-1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 0)*bxg(j, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 0)*bxg(j+1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 0)*bxg(j+2, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 0)*bxg(j-1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 0)*bxg(j+2, k0, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 0)*bxg(j-1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j+1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 0)*bxg(j+2, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 1)*bxg(j-1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 1)*bxg(j, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 1)*bxg(j+1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 1)*bxg(j+2, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 1)*bxg(j-1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j+1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 1)*bxg(j+2, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 1)*bxg(j-1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j+1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 1)*bxg(j+2, k0+1, l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, -1)*byg(j0-1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, -1)*byg(j0, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, -1)*byg(j0+1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, -1)*byg(j0-1, k, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, -1)*byg(j0, k, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, -1)*byg(j0+1, k, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, -1)*byg(j0-1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, -1)*byg(j0, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, -1)*byg(j0+1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, -1)*byg(j0-1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, -1)*byg(j0, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, -1)*byg(j0+1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 0)*byg(j0-1, k-1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 0)*byg(j0, k-1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 0)*byg(j0+1, k-1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 0)*byg(j0-1, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0+1, k, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 0)*byg(j0-1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0+1, k+1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 0)*byg(j0-1, k+2, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 0)*byg(j0, k+2, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 0)*byg(j0+1, k+2, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 1)*byg(j0-1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 1)*byg(j0, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 1)*byg(j0+1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 1)*byg(j0-1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0, k, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0+1, k, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 1)*byg(j0-1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0+1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 1)*byg(j0-1, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 1)*byg(j0, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 1)*byg(j0+1, k+2, l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, -1)*bzg(j0-1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, -1)*bzg(j0, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, -1)*bzg(j0+1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, -1)*bzg(j0-1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, -1)*bzg(j0, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, -1)*bzg(j0+1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, -1)*bzg(j0-1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, -1)*bzg(j0, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, -1)*bzg(j0+1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 0)*bzg(j0-1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 0)*bzg(j0, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 0)*bzg(j0+1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 0)*bzg(j0-1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0+1, k0, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 0)*bzg(j0-1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0+1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 1)*bzg(j0-1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 1)*bzg(j0, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 1)*bzg(j0+1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 1)*bzg(j0-1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0+1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 1)*bzg(j0-1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0+1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 2)*bzg(j0-1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 2)*bzg(j0, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 2)*bzg(j0+1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 2)*bzg(j0-1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 2)*bzg(j0, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 2)*bzg(j0+1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 2)*bzg(j0-1, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 2)*bzg(j0, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 2)*bzg(j0+1, k0+1, l+2)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

  ELSE
    ! Loop over the particles by block
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
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint
        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(n, -1) = onesixth*oyintsq*oyint
        sy0(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, -1)*bxg(j-1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, -1)*bxg(j, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, -1)*bxg(j+1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, -1)*bxg(j+2, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, -1)*bxg(j-1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, -1)*bxg(j, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, -1)*bxg(j+1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, -1)*bxg(j+2, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, -1)*bxg(j-1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, -1)*bxg(j, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, -1)*bxg(j+1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, -1)*bxg(j+2, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, -1)*bxg(j-1, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, -1)*bxg(j, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, -1)*bxg(j+1, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, -1)*bxg(j+2, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 0)*bxg(j-1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 0)*bxg(j, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 0)*bxg(j+1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 0)*bxg(j+2, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 0)*bxg(j-1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 0)*bxg(j+2, k0, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 0)*bxg(j-1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j+1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 0)*bxg(j+2, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, 0)*bxg(j-1, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, 0)*bxg(j, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, 0)*bxg(j+1, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, 0)*bxg(j+2, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 1)*bxg(j-1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 1)*bxg(j, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 1)*bxg(j+1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 1)*bxg(j+2, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 1)*bxg(j-1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j+1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 1)*bxg(j+2, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 1)*bxg(j-1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j+1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 1)*bxg(j+2, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, 1)*bxg(j-1, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, 1)*bxg(j, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, 1)*bxg(j+1, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, 1)*bxg(j+2, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 2)*bxg(j-1, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 2)*bxg(j, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 2)*bxg(j+1, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 2)*bxg(j+2, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 2)*bxg(j-1, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 2)*bxg(j, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 2)*bxg(j+1, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 2)*bxg(j+2, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 2)*bxg(j-1, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 2)*bxg(j, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 2)*bxg(j+1, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 2)*bxg(j+2, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, 2)*bxg(j-1, k0+2, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, 2)*bxg(j, k0+2, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, 2)*bxg(j+1, k0+2, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, 2)*bxg(j+2, k0+2, l0+2)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, -1)*byg(j0-1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, -1)*byg(j0, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, -1)*byg(j0+1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, -1)*byg(j0+2, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, -1)*byg(j0-1, k, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, -1)*byg(j0, k, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, -1)*byg(j0+1, k, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, -1)*byg(j0+2, k, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, -1)*byg(j0-1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, -1)*byg(j0, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, -1)*byg(j0+1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, -1)*byg(j0+2, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, -1)*byg(j0-1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, -1)*byg(j0, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, -1)*byg(j0+1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, -1)*byg(j0+2, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 0)*byg(j0-1, k-1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 0)*byg(j0, k-1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 0)*byg(j0+1, k-1, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, 0)*byg(j0+2, k-1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 0)*byg(j0-1, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0+1, k, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, 0)*byg(j0+2, k, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 0)*byg(j0-1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0+1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, 0)*byg(j0+2, k+1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 0)*byg(j0-1, k+2, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 0)*byg(j0, k+2, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 0)*byg(j0+1, k+2, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, 0)*byg(j0+2, k+2, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 1)*byg(j0-1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 1)*byg(j0, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 1)*byg(j0+1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, 1)*byg(j0+2, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 1)*byg(j0-1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0, k, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0+1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, 1)*byg(j0+2, k, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 1)*byg(j0-1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0+1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, 1)*byg(j0+2, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 1)*byg(j0-1, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 1)*byg(j0, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 1)*byg(j0+1, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, 1)*byg(j0+2, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 2)*byg(j0-1, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 2)*byg(j0, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 2)*byg(j0+1, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, 2)*byg(j0+2, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 2)*byg(j0-1, k, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 2)*byg(j0, k, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 2)*byg(j0+1, k, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, 2)*byg(j0+2, k, l0+2)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 2)*byg(j0-1, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 2)*byg(j0, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 2)*byg(j0+1, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, 2)*byg(j0+2, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 2)*byg(j0-1, k+2, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 2)*byg(j0, k+2, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 2)*byg(j0+1, k+2, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, 2)*byg(j0+2, k+2, l0+2)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, -1)*bzg(j0-1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, -1)*bzg(j0, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, -1)*bzg(j0+1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, -1)*bzg(j0+2, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, -1)*bzg(j0-1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, -1)*bzg(j0, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, -1)*bzg(j0+1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, -1)*bzg(j0+2, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, -1)*bzg(j0-1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, -1)*bzg(j0, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, -1)*bzg(j0+1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, -1)*bzg(j0+2, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, -1)*bzg(j0-1, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, -1)*bzg(j0, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, -1)*bzg(j0+1, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, -1)*bzg(j0+2, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 0)*bzg(j0-1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 0)*bzg(j0, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 0)*bzg(j0+1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, 0)*bzg(j0+2, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 0)*bzg(j0-1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0+1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, 0)*bzg(j0+2, k0, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 0)*bzg(j0-1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0+1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, 0)*bzg(j0+2, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, 0)*bzg(j0-1, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, 0)*bzg(j0, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, 0)*bzg(j0+1, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, 0)*bzg(j0+2, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 1)*bzg(j0-1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 1)*bzg(j0, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 1)*bzg(j0+1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, 1)*bzg(j0+2, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 1)*bzg(j0-1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0+1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, 1)*bzg(j0+2, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 1)*bzg(j0-1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0+1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, 1)*bzg(j0+2, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, 1)*bzg(j0-1, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, 1)*bzg(j0, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, 1)*bzg(j0+1, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, 1)*bzg(j0+2, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 2)*bzg(j0-1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 2)*bzg(j0, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 2)*bzg(j0+1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, 2)*bzg(j0+2, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 2)*bzg(j0-1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 2)*bzg(j0, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 2)*bzg(j0+1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, 2)*bzg(j0+2, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 2)*bzg(j0-1, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 2)*bzg(j0, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 2)*bzg(j0+1, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, 2)*bzg(j0+2, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, 2)*bzg(j0-1, k0+2, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, 2)*bzg(j0, k0+2, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, 2)*bzg(j0+1, k0+2, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, 2)*bzg(j0+2, k0+2, l+2)
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
!> Vectorized gathering of the electric fields (version 2) from Yee grid
!> ("energy conserving")
!> on particles at order 3.
!
!> @details
!> This subroutine is the second vectorized version of
!> the electric field gathering gete3d_energy_conserving_3_3_3.
!>
!> The gathering has been optimized to perform less operations, less memory access
!> and more FMA so that the vectorization is at the end more efficient.
!
!> @date
!> Creation 12/03/2016
!
!> @author
!> Mathieu Lobet
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] ex, ey, ez particle electric field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE gete3d_energy_conserving_vec2_3_3_3(np, xp, yp, zp, ex, ey, ez, xmin,      &
  ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg, ezg, lvect,  &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                           :: lvect
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, j, k, l
  INTEGER(isp)                           :: jj, kk, ll
  INTEGER(isp)                           :: j0, k0, l0
  INTEGER(isp)                           :: nn, n
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq,     &
  ozint
  REAL(num)                              :: oxintsq, oyintsq, ozintsq
  REAL(num)                              :: a
  REAL(num), DIMENSION(lvect, -1:2)       :: sx, sx0
  REAL(num), DIMENSION(lvect, -1:2)       :: sy, sy0
  REAL(num), DIMENSION(lvect, -1:2)       :: sz, sz0
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

  !write(0, *) 'l_lower_order_in_v ', l_lower_order_in_v
  !write(0, *) 'sum(xp)', sum(xp), sum(yp), sum(zp)
  !write(0, *) 'sum(exg)', sum(exg), sum(eyg), sum(ezg)

  IF (l_lower_order_in_v ) THEN

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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        yintsq = yint*yint
        sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        a = (sx0(n, -1)*exg(j0-1, k-1, l-1) + sx0(n, 0)*exg(j0, k-1, l-1) + sx0(n,    &
        1)*exg(j0+1, k-1, l-1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l-1) + sx0(n, 0)*exg(j0, k, l-1) + sx0(n,    &
        1)*exg(j0+1, k, l-1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l-1) + sx0(n, 0)*exg(j0, k+1, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l-1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l-1) + sx0(n, 0)*exg(j0, k+2, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l-1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l) + sx0(n, 0)*exg(j0, k-1, l) + sx0(n,        &
        1)*exg(j0+1, k-1, l))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l) + sx0(n, 0)*exg(j0, k, l) + sx0(n,        &
        1)*exg(j0+1, k, l))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l) + sx0(n, 0)*exg(j0, k+1, l) + sx0(n,    &
        1)*exg(j0+1, k+1, l))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l) + sx0(n, 0)*exg(j0, k+2, l) + sx0(n,    &
        1)*exg(j0+1, k+2, l))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+1) + sx0(n, 0)*exg(j0, k-1, l+1) + sx0(n,    &
        1)*exg(j0+1, k-1, l+1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+1) + sx0(n, 0)*exg(j0, k, l+1) + sx0(n,    &
        1)*exg(j0+1, k, l+1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+1) + sx0(n, 0)*exg(j0, k+1, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+1) + sx0(n, 0)*exg(j0, k+2, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+2) + sx0(n, 0)*exg(j0, k-1, l+2) + sx0(n,    &
        1)*exg(j0+1, k-1, l+2))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+2) + sx0(n, 0)*exg(j0, k, l+2) + sx0(n,    &
        1)*exg(j0+1, k, l+2))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+2) + sx0(n, 0)*exg(j0, k+1, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+2))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+2) + sx0(n, 0)*exg(j0, k+2, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+2))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 2)

        ! Compute Ey on particle
        a = (sx(n, -1)*eyg(j-1, k0-1, l-1) + sx(n, 0)*eyg(j, k0-1, l-1) + sx(n,       &
        1)*eyg(j+1, k0-1, l-1) + sx(n, 2)*eyg(j+2, k0-1, l-1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l-1) + sx(n, 0)*eyg(j, k0, l-1) + sx(n,       &
        1)*eyg(j+1, k0, l-1) + sx(n, 2)*eyg(j+2, k0, l-1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l-1) + sx(n, 0)*eyg(j, k0+1, l-1) + sx(n,   &
        1)*eyg(j+1, k0+1, l-1) + sx(n, 2)*eyg(j+2, k0+1, l-1))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, -1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l) + sx(n, 0)*eyg(j, k0-1, l) + sx(n,           &
        1)*eyg(j+1, k0-1, l) + sx(n, 2)*eyg(j+2, k0-1, l))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l) + sx(n, 0)*eyg(j, k0, l) + sx(n,           &
        1)*eyg(j+1, k0, l) + sx(n, 2)*eyg(j+2, k0, l))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l) + sx(n, 0)*eyg(j, k0+1, l) + sx(n,       &
        1)*eyg(j+1, k0+1, l) + sx(n, 2)*eyg(j+2, k0+1, l))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+1) + sx(n, 0)*eyg(j, k0-1, l+1) + sx(n,       &
        1)*eyg(j+1, k0-1, l+1) + sx(n, 2)*eyg(j+2, k0-1, l+1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+1) + sx(n, 0)*eyg(j, k0, l+1) + sx(n,       &
        1)*eyg(j+1, k0, l+1) + sx(n, 2)*eyg(j+2, k0, l+1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+1) + sx(n, 0)*eyg(j, k0+1, l+1) + sx(n,   &
        1)*eyg(j+1, k0+1, l+1) + sx(n, 2)*eyg(j+2, k0+1, l+1))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+2) + sx(n, 0)*eyg(j, k0-1, l+2) + sx(n,       &
        1)*eyg(j+1, k0-1, l+2) + sx(n, 2)*eyg(j+2, k0-1, l+2))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+2) + sx(n, 0)*eyg(j, k0, l+2) + sx(n,       &
        1)*eyg(j+1, k0, l+2) + sx(n, 2)*eyg(j+2, k0, l+2))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+2) + sx(n, 0)*eyg(j, k0+1, l+2) + sx(n,   &
        1)*eyg(j+1, k0+1, l+2) + sx(n, 2)*eyg(j+2, k0+1, l+2))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 2)

        ! Compute Ez on particle
        a = (sx(n, -1)*ezg(j-1, k-1, l0-1) + sx(n, 0)*ezg(j, k-1, l0-1) + sx(n,       &
        1)*ezg(j+1, k-1, l0-1) + sx(n, 2)*ezg(j+2, k-1, l0-1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0-1) + sx(n, 0)*ezg(j, k, l0-1) + sx(n,       &
        1)*ezg(j+1, k, l0-1) + sx(n, 2)*ezg(j+2, k, l0-1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0-1) + sx(n, 0)*ezg(j, k+1, l0-1) + sx(n,   &
        1)*ezg(j+1, k+1, l0-1) + sx(n, 2)*ezg(j+2, k+1, l0-1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0-1) + sx(n, 0)*ezg(j, k+2, l0-1) + sx(n,   &
        1)*ezg(j+1, k+2, l0-1) + sx(n, 2)*ezg(j+2, k+2, l0-1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*ezg(j-1, k-1, l0) + sx(n, 0)*ezg(j, k-1, l0) + sx(n,           &
        1)*ezg(j+1, k-1, l0) + sx(n, 2)*ezg(j+2, k-1, l0))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0) + sx(n, 0)*ezg(j, k, l0) + sx(n,           &
        1)*ezg(j+1, k, l0) + sx(n, 2)*ezg(j+2, k, l0))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0) + sx(n, 0)*ezg(j, k+1, l0) + sx(n,       &
        1)*ezg(j+1, k+1, l0) + sx(n, 2)*ezg(j+2, k+1, l0))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0) + sx(n, 0)*ezg(j, k+2, l0) + sx(n,       &
        1)*ezg(j+1, k+2, l0) + sx(n, 2)*ezg(j+2, k+2, l0))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*ezg(j-1, k-1, l0+1) + sx(n, 0)*ezg(j, k-1, l0+1) + sx(n,       &
        1)*ezg(j+1, k-1, l0+1) + sx(n, 2)*ezg(j+2, k-1, l0+1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0+1) + sx(n, 0)*ezg(j, k, l0+1) + sx(n,       &
        1)*ezg(j+1, k, l0+1) + sx(n, 2)*ezg(j+2, k, l0+1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0+1) + sx(n, 0)*ezg(j, k+1, l0+1) + sx(n,   &
        1)*ezg(j+1, k+1, l0+1) + sx(n, 2)*ezg(j+2, k+1, l0+1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0+1) + sx(n, 0)*ezg(j, k+2, l0+1) + sx(n,   &
        1)*ezg(j+1, k+2, l0+1) + sx(n, 2)*ezg(j+2, k+2, l0+1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 1)
      ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif

    ENDDO

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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(n, -1) = onesixth*oyintsq*oyint
        sy0(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        a = (sx0(n, -1)*exg(j0-1, k-1, l-1) + sx0(n, 0)*exg(j0, k-1, l-1) + sx0(n,    &
        1)*exg(j0+1, k-1, l-1) + sx0(n, 2)*exg(j0+2, k-1, l-1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l-1) + sx0(n, 0)*exg(j0, k, l-1) + sx0(n,    &
        1)*exg(j0+1, k, l-1) + sx0(n, 2)*exg(j0+2, k, l-1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l-1) + sx0(n, 0)*exg(j0, k+1, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l-1) + sx0(n, 2)*exg(j0+2, k+1, l-1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l-1) + sx0(n, 0)*exg(j0, k+2, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l-1) + sx0(n, 2)*exg(j0+2, k+2, l-1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l) + sx0(n, 0)*exg(j0, k-1, l) + sx0(n,        &
        1)*exg(j0+1, k-1, l) + sx0(n, 2)*exg(j0+2, k-1, l))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l) + sx0(n, 0)*exg(j0, k, l) + sx0(n,        &
        1)*exg(j0+1, k, l) + sx0(n, 2)*exg(j0+2, k, l))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l) + sx0(n, 0)*exg(j0, k+1, l) + sx0(n,    &
        1)*exg(j0+1, k+1, l) + sx0(n, 2)*exg(j0+2, k+1, l))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l) + sx0(n, 0)*exg(j0, k+2, l) + sx0(n,    &
        1)*exg(j0+1, k+2, l) + sx0(n, 2)*exg(j0+2, k+2, l))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+1) + sx0(n, 0)*exg(j0, k-1, l+1) + sx0(n,    &
        1)*exg(j0+1, k-1, l+1) + sx0(n, 2)*exg(j0+2, k-1, l+1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+1) + sx0(n, 0)*exg(j0, k, l+1) + sx0(n,    &
        1)*exg(j0+1, k, l+1) + sx0(n, 2)*exg(j0+2, k, l+1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+1) + sx0(n, 0)*exg(j0, k+1, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+1) + sx0(n, 2)*exg(j0+2, k+1, l+1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+1) + sx0(n, 0)*exg(j0, k+2, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+1) + sx0(n, 2)*exg(j0+2, k+2, l+1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+2) + sx0(n, 0)*exg(j0, k-1, l+2) + sx0(n,    &
        1)*exg(j0+1, k-1, l+2) + sx0(n, 2)*exg(j0+2, k-1, l+2))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+2) + sx0(n, 0)*exg(j0, k, l+2) + sx0(n,    &
        1)*exg(j0+1, k, l+2) + sx0(n, 2)*exg(j0+2, k, l+2))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+2) + sx0(n, 0)*exg(j0, k+1, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+2) + sx0(n, 2)*exg(j0+2, k+1, l+2))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+2) + sx0(n, 0)*exg(j0, k+2, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+2) + sx0(n, 2)*exg(j0+2, k+2, l+2))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 2)

        ! Compute Ey on particle
        a = (sx(n, -1)*eyg(j-1, k0-1, l-1) + sx(n, 0)*eyg(j, k0-1, l-1) + sx(n,       &
        1)*eyg(j+1, k0-1, l-1) + sx(n, 2)*eyg(j+2, k0-1, l-1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l-1) + sx(n, 0)*eyg(j, k0, l-1) + sx(n,       &
        1)*eyg(j+1, k0, l-1) + sx(n, 2)*eyg(j+2, k0, l-1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l-1) + sx(n, 0)*eyg(j, k0+1, l-1) + sx(n,   &
        1)*eyg(j+1, k0+1, l-1) + sx(n, 2)*eyg(j+2, k0+1, l-1))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l-1) + sx(n, 0)*eyg(j, k0+2, l-1) + sx(n,   &
        1)*eyg(j+1, k0+2, l-1) + sx(n, 2)*eyg(j+2, k0+2, l-1))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, -1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l) + sx(n, 0)*eyg(j, k0-1, l) + sx(n,           &
        1)*eyg(j+1, k0-1, l) + sx(n, 2)*eyg(j+2, k0-1, l))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l) + sx(n, 0)*eyg(j, k0, l) + sx(n,           &
        1)*eyg(j+1, k0, l) + sx(n, 2)*eyg(j+2, k0, l))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l) + sx(n, 0)*eyg(j, k0+1, l) + sx(n,       &
        1)*eyg(j+1, k0+1, l) + sx(n, 2)*eyg(j+2, k0+1, l))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l) + sx(n, 0)*eyg(j, k0+2, l) + sx(n,       &
        1)*eyg(j+1, k0+2, l) + sx(n, 2)*eyg(j+2, k0+2, l))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+1) + sx(n, 0)*eyg(j, k0-1, l+1) + sx(n,       &
        1)*eyg(j+1, k0-1, l+1) + sx(n, 2)*eyg(j+2, k0-1, l+1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+1) + sx(n, 0)*eyg(j, k0, l+1) + sx(n,       &
        1)*eyg(j+1, k0, l+1) + sx(n, 2)*eyg(j+2, k0, l+1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+1) + sx(n, 0)*eyg(j, k0+1, l+1) + sx(n,   &
        1)*eyg(j+1, k0+1, l+1) + sx(n, 2)*eyg(j+2, k0+1, l+1))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l+1) + sx(n, 0)*eyg(j, k0+2, l+1) + sx(n,   &
        1)*eyg(j+1, k0+2, l+1) + sx(n, 2)*eyg(j+2, k0+2, l+1))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, 1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+2) + sx(n, 0)*eyg(j, k0-1, l+2) + sx(n,       &
        1)*eyg(j+1, k0-1, l+2) + sx(n, 2)*eyg(j+2, k0-1, l+2))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+2) + sx(n, 0)*eyg(j, k0, l+2) + sx(n,       &
        1)*eyg(j+1, k0, l+2) + sx(n, 2)*eyg(j+2, k0, l+2))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+2) + sx(n, 0)*eyg(j, k0+1, l+2) + sx(n,   &
        1)*eyg(j+1, k0+1, l+2) + sx(n, 2)*eyg(j+2, k0+1, l+2))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l+2) + sx(n, 0)*eyg(j, k0+2, l+2) + sx(n,   &
        1)*eyg(j+1, k0+2, l+2) + sx(n, 2)*eyg(j+2, k0+2, l+2))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, 2)

        ! Compute Ez on particle
        a = (sx(n, -1)*ezg(j-1, k-1, l0-1) + sx(n, 0)*ezg(j, k-1, l0-1) + sx(n,       &
        1)*ezg(j+1, k-1, l0-1) + sx(n, 2)*ezg(j+2, k-1, l0-1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0-1) + sx(n, 0)*ezg(j, k, l0-1) + sx(n,       &
        1)*ezg(j+1, k, l0-1) + sx(n, 2)*ezg(j+2, k, l0-1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0-1) + sx(n, 0)*ezg(j, k+1, l0-1) + sx(n,   &
        1)*ezg(j+1, k+1, l0-1) + sx(n, 2)*ezg(j+2, k+1, l0-1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0-1) + sx(n, 0)*ezg(j, k+2, l0-1) + sx(n,   &
        1)*ezg(j+1, k+2, l0-1) + sx(n, 2)*ezg(j+2, k+2, l0-1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*ezg(j-1, k-1, l0) + sx(n, 0)*ezg(j, k-1, l0) + sx(n,           &
        1)*ezg(j+1, k-1, l0) + sx(n, 2)*ezg(j+2, k-1, l0))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0) + sx(n, 0)*ezg(j, k, l0) + sx(n,           &
        1)*ezg(j+1, k, l0) + sx(n, 2)*ezg(j+2, k, l0))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0) + sx(n, 0)*ezg(j, k+1, l0) + sx(n,       &
        1)*ezg(j+1, k+1, l0) + sx(n, 2)*ezg(j+2, k+1, l0))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0) + sx(n, 0)*ezg(j, k+2, l0) + sx(n,       &
        1)*ezg(j+1, k+2, l0) + sx(n, 2)*ezg(j+2, k+2, l0))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*ezg(j-1, k-1, l0+1) + sx(n, 0)*ezg(j, k-1, l0+1) + sx(n,       &
        1)*ezg(j+1, k-1, l0+1) + sx(n, 2)*ezg(j+2, k-1, l0+1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0+1) + sx(n, 0)*ezg(j, k, l0+1) + sx(n,       &
        1)*ezg(j+1, k, l0+1) + sx(n, 2)*ezg(j+2, k, l0+1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0+1) + sx(n, 0)*ezg(j, k+1, l0+1) + sx(n,   &
        1)*ezg(j+1, k+1, l0+1) + sx(n, 2)*ezg(j+2, k+1, l0+1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0+1) + sx(n, 0)*ezg(j, k+2, l0+1) + sx(n,   &
        1)*ezg(j+1, k+2, l0+1) + sx(n, 2)*ezg(j+2, k+2, l0+1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 1)
        a = (sx(n, -1)*ezg(j-1, k-1, l0+2) + sx(n, 0)*ezg(j, k-1, l0+2) + sx(n,       &
        1)*ezg(j+1, k-1, l0+2) + sx(n, 2)*ezg(j+2, k-1, l0+2))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0+2) + sx(n, 0)*ezg(j, k, l0+2) + sx(n,       &
        1)*ezg(j+1, k, l0+2) + sx(n, 2)*ezg(j+2, k, l0+2))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0+2) + sx(n, 0)*ezg(j, k+1, l0+2) + sx(n,   &
        1)*ezg(j+1, k+1, l0+2) + sx(n, 2)*ezg(j+2, k+1, l0+2))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0+2) + sx(n, 0)*ezg(j, k+2, l0+2) + sx(n,   &
        1)*ezg(j+1, k+2, l0+2) + sx(n, 2)*ezg(j+2, k+2, l0+2))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 2)
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
!> Vectorized gathering of magnetic field (version 2) from Yee grid ("energy conserving")
!> on particles at order 3.
!
!> @details
!> This subroutine is the second vectorized version of
!> the electric field gathering gete3d_energy_conserving_3_3_3.
!>
!> The gathering has been optimized to perform less operations, less memory access
!> and more FMA so that the vectorization is at the end more efficient.
!
!> @date
!> Creation 12/03/2016
!
!> @author
!> Mathieu Lobet
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position
!> @param[inout] bx, by, bz particle magnetic field
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space step
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v decrease the interpolation order if True
!
! ________________________________________________________________________________________
SUBROUTINE getb3d_energy_conserving_vec2_3_3_3(np, xp, yp, zp, bx, by, bz, xmin,      &
  ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, bxg, byg, bzg, lvect,  &
  l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                         :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                         :: lvect
  REAL(num), DIMENSION(np)             :: xp, yp, zp, bx, by, bz
  LOGICAL(lp)                          :: l_lower_order_in_v, l_nodal
  REAL(num)                            :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                            :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                         :: ip, j, k, l
  INTEGER(isp)                         :: j0, k0, l0
  INTEGER(isp)                         :: jj, kk, ll
  INTEGER(isp)                         :: n, nn
  REAL(num)                            :: dxi, dyi, dzi, x, y, z
  REAL(num)                            :: xint, yint, zint
  REAL(num)                            :: xintsq, oxint, yintsq, oyint, zintsq
  REAL(num)                            :: ozint, oxintsq, oyintsq, ozintsq
  REAL(num)                            :: a
  REAL(num), DIMENSION(lvect, -1:2)     :: sx, sx0
  REAL(num), DIMENSION(lvect, -1:2)     :: sy, sy0
  REAL(num), DIMENSION(lvect, -1:2)     :: sz, sz0
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

    ! Loop over the particles by block
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0
        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2
        yintsq = yint*yint
        sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2
        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Bx on particle
        a = (sx(n, -1)*bxg(j-1, k0-1, l0-1) + sx(n, 0)*bxg(j, k0-1, l0-1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0-1) + sx(n, 2)*bxg(j+2, k0-1, l0-1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0-1) + sx(n, 0)*bxg(j, k0, l0-1) + sx(n,     &
        1)*bxg(j+1, k0, l0-1) + sx(n, 2)*bxg(j+2, k0, l0-1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0-1) + sx(n, 0)*bxg(j, k0+1, l0-1) + sx(n, &
        1)*bxg(j+1, k0+1, l0-1) + sx(n, 2)*bxg(j+2, k0+1, l0-1))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0) + sx(n, 0)*bxg(j, k0-1, l0) + sx(n,         &
        1)*bxg(j+1, k0-1, l0) + sx(n, 2)*bxg(j+2, k0-1, l0))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0) + sx(n, 0)*bxg(j, k0, l0) + sx(n,         &
        1)*bxg(j+1, k0, l0) + sx(n, 2)*bxg(j+2, k0, l0))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0) + sx(n, 0)*bxg(j, k0+1, l0) + sx(n,     &
        1)*bxg(j+1, k0+1, l0) + sx(n, 2)*bxg(j+2, k0+1, l0))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0+1) + sx(n, 0)*bxg(j, k0-1, l0+1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0+1) + sx(n, 2)*bxg(j+2, k0-1, l0+1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0+1) + sx(n, 0)*bxg(j, k0, l0+1) + sx(n,     &
        1)*bxg(j+1, k0, l0+1) + sx(n, 2)*bxg(j+2, k0, l0+1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0+1) + sx(n, 0)*bxg(j, k0+1, l0+1) + sx(n, &
        1)*bxg(j+1, k0+1, l0+1) + sx(n, 2)*bxg(j+2, k0+1, l0+1))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, 1)

        ! Compute By on particle
        a = (sx0(n, -1)*byg(j0-1, k-1, l0-1) + sx0(n, 0)*byg(j0, k-1, l0-1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0-1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0-1) + sx0(n, 0)*byg(j0, k, l0-1) + sx0(n,  &
        1)*byg(j0+1, k, l0-1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0-1) + sx0(n, 0)*byg(j0, k+1, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0-1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0-1) + sx0(n, 0)*byg(j0, k+2, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0-1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, -1)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0) + sx0(n, 0)*byg(j0, k-1, l0) + sx0(n,      &
        1)*byg(j0+1, k-1, l0))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0) + sx0(n, 0)*byg(j0, k, l0) + sx0(n,      &
        1)*byg(j0+1, k, l0))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0) + sx0(n, 0)*byg(j0, k+1, l0) + sx0(n,  &
        1)*byg(j0+1, k+1, l0))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0) + sx0(n, 0)*byg(j0, k+2, l0) + sx0(n,  &
        1)*byg(j0+1, k+2, l0))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0+1) + sx0(n, 0)*byg(j0, k-1, l0+1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0+1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0+1) + sx0(n, 0)*byg(j0, k, l0+1) + sx0(n,  &
        1)*byg(j0+1, k, l0+1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0+1) + sx0(n, 0)*byg(j0, k+1, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0+1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0+1) + sx0(n, 0)*byg(j0, k+2, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0+1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 1)

        ! Compute Bz on particle
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l-1) + sx0(n, 0)*bzg(j0, k0-1, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l-1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l-1) + sx0(n, 0)*bzg(j0, k0, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0, l-1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l-1) + sx0(n, 0)*bzg(j0, k0+1, l-1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l-1))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l) + sx0(n, 0)*bzg(j0, k0-1, l) + sx0(n,      &
        1)*bzg(j0+1, k0-1, l))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l) + sx0(n, 0)*bzg(j0, k0, l) + sx0(n,      &
        1)*bzg(j0+1, k0, l))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l) + sx0(n, 0)*bzg(j0, k0+1, l) + sx0(n,  &
        1)*bzg(j0+1, k0+1, l))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+1) + sx0(n, 0)*bzg(j0, k0-1, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+1) + sx0(n, 0)*bzg(j0, k0, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0, l+1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+1) + sx0(n, 0)*bzg(j0, k0+1, l+1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+1))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+2) + sx0(n, 0)*bzg(j0, k0-1, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+2))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+2) + sx0(n, 0)*bzg(j0, k0, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0, l+2))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+2) + sx0(n, 0)*bzg(j0, k0+1, l+2) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+2))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 2)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO

  ELSE

    ! Loop over the particles by block
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
#if defined __INTEL_COMPILER
      !DIR$ IVDEP
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint
        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(n, -1) = onesixth*oyintsq*oyint
        sy0(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Bx on particle
        a = (sx(n, -1)*bxg(j-1, k0-1, l0-1) + sx(n, 0)*bxg(j, k0-1, l0-1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0-1) + sx(n, 2)*bxg(j+2, k0-1, l0-1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0-1) + sx(n, 0)*bxg(j, k0, l0-1) + sx(n,     &
        1)*bxg(j+1, k0, l0-1) + sx(n, 2)*bxg(j+2, k0, l0-1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0-1) + sx(n, 0)*bxg(j, k0+1, l0-1) + sx(n, &
        1)*bxg(j+1, k0+1, l0-1) + sx(n, 2)*bxg(j+2, k0+1, l0-1))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0-1) + sx(n, 0)*bxg(j, k0+2, l0-1) + sx(n, &
        1)*bxg(j+1, k0+2, l0-1) + sx(n, 2)*bxg(j+2, k0+2, l0-1))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0) + sx(n, 0)*bxg(j, k0-1, l0) + sx(n,         &
        1)*bxg(j+1, k0-1, l0) + sx(n, 2)*bxg(j+2, k0-1, l0))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0) + sx(n, 0)*bxg(j, k0, l0) + sx(n,         &
        1)*bxg(j+1, k0, l0) + sx(n, 2)*bxg(j+2, k0, l0))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0) + sx(n, 0)*bxg(j, k0+1, l0) + sx(n,     &
        1)*bxg(j+1, k0+1, l0) + sx(n, 2)*bxg(j+2, k0+1, l0))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0) + sx(n, 0)*bxg(j, k0+2, l0) + sx(n,     &
        1)*bxg(j+1, k0+2, l0) + sx(n, 2)*bxg(j+2, k0+2, l0))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0+1) + sx(n, 0)*bxg(j, k0-1, l0+1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0+1) + sx(n, 2)*bxg(j+2, k0-1, l0+1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0+1) + sx(n, 0)*bxg(j, k0, l0+1) + sx(n,     &
        1)*bxg(j+1, k0, l0+1) + sx(n, 2)*bxg(j+2, k0, l0+1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0+1) + sx(n, 0)*bxg(j, k0+1, l0+1) + sx(n, &
        1)*bxg(j+1, k0+1, l0+1) + sx(n, 2)*bxg(j+2, k0+1, l0+1))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0+1) + sx(n, 0)*bxg(j, k0+2, l0+1) + sx(n, &
        1)*bxg(j+1, k0+2, l0+1) + sx(n, 2)*bxg(j+2, k0+2, l0+1))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, 1)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0+2) + sx(n, 0)*bxg(j, k0-1, l0+2) + sx(n,     &
        1)*bxg(j+1, k0-1, l0+2) + sx(n, 2)*bxg(j+2, k0-1, l0+2))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0+2) + sx(n, 0)*bxg(j, k0, l0+2) + sx(n,     &
        1)*bxg(j+1, k0, l0+2) + sx(n, 2)*bxg(j+2, k0, l0+2))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0+2) + sx(n, 0)*bxg(j, k0+1, l0+2) + sx(n, &
        1)*bxg(j+1, k0+1, l0+2) + sx(n, 2)*bxg(j+2, k0+1, l0+2))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0+2) + sx(n, 0)*bxg(j, k0+2, l0+2) + sx(n, &
        1)*bxg(j+1, k0+2, l0+2) + sx(n, 2)*bxg(j+2, k0+2, l0+2))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, 2)

        ! Compute By on particle
        a = (sx0(n, -1)*byg(j0-1, k-1, l0-1) + sx0(n, 0)*byg(j0, k-1, l0-1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0-1) + sx0(n, 2)*byg(j0+2, k-1, l0-1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0-1) + sx0(n, 0)*byg(j0, k, l0-1) + sx0(n,  &
        1)*byg(j0+1, k, l0-1) + sx0(n, 2)*byg(j0+2, k, l0-1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0-1) + sx0(n, 0)*byg(j0, k+1, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0-1) + sx0(n, 2)*byg(j0+2, k+1, l0-1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0-1) + sx0(n, 0)*byg(j0, k+2, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0-1) + sx0(n, 2)*byg(j0+2, k+2, l0-1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, -1)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0) + sx0(n, 0)*byg(j0, k-1, l0) + sx0(n,      &
        1)*byg(j0+1, k-1, l0) + sx0(n, 2)*byg(j0+2, k-1, l0))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0) + sx0(n, 0)*byg(j0, k, l0) + sx0(n,      &
        1)*byg(j0+1, k, l0) + sx0(n, 2)*byg(j0+2, k, l0))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0) + sx0(n, 0)*byg(j0, k+1, l0) + sx0(n,  &
        1)*byg(j0+1, k+1, l0) + sx0(n, 2)*byg(j0+2, k+1, l0))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0) + sx0(n, 0)*byg(j0, k+2, l0) + sx0(n,  &
        1)*byg(j0+1, k+2, l0) + sx0(n, 2)*byg(j0+2, k+2, l0))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0+1) + sx0(n, 0)*byg(j0, k-1, l0+1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0+1) + sx0(n, 2)*byg(j0+2, k-1, l0+1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0+1) + sx0(n, 0)*byg(j0, k, l0+1) + sx0(n,  &
        1)*byg(j0+1, k, l0+1) + sx0(n, 2)*byg(j0+2, k, l0+1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0+1) + sx0(n, 0)*byg(j0, k+1, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0+1) + sx0(n, 2)*byg(j0+2, k+1, l0+1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0+1) + sx0(n, 0)*byg(j0, k+2, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0+1) + sx0(n, 2)*byg(j0+2, k+2, l0+1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 1)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0+2) + sx0(n, 0)*byg(j0, k-1, l0+2) + sx0(n,  &
        1)*byg(j0+1, k-1, l0+2) + sx0(n, 2)*byg(j0+2, k-1, l0+2))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0+2) + sx0(n, 0)*byg(j0, k, l0+2) + sx0(n,  &
        1)*byg(j0+1, k, l0+2) + sx0(n, 2)*byg(j0+2, k, l0+2))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0+2) + sx0(n, 0)*byg(j0, k+1, l0+2) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0+2) + sx0(n, 2)*byg(j0+2, k+1, l0+2))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0+2) + sx0(n, 0)*byg(j0, k+2, l0+2) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0+2) + sx0(n, 2)*byg(j0+2, k+2, l0+2))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 2)

        ! Compute Bz on particle
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l-1) + sx0(n, 0)*bzg(j0, k0-1, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l-1) + sx0(n, 2)*bzg(j0+2, k0-1, l-1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l-1) + sx0(n, 0)*bzg(j0, k0, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0, l-1) + sx0(n, 2)*bzg(j0+2, k0, l-1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l-1) + sx0(n, 0)*bzg(j0, k0+1, l-1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l-1) + sx0(n, 2)*bzg(j0+2, k0+1, l-1))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l-1) + sx0(n, 0)*bzg(j0, k0+2, l-1) +     &
        sx0(n, 1)*bzg(j0+1, k0+2, l-1) + sx0(n, 2)*bzg(j0+2, k0+2, l-1))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l) + sx0(n, 0)*bzg(j0, k0-1, l) + sx0(n,      &
        1)*bzg(j0+1, k0-1, l) + sx0(n, 2)*bzg(j0+2, k0-1, l))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l) + sx0(n, 0)*bzg(j0, k0, l) + sx0(n,      &
        1)*bzg(j0+1, k0, l) + sx0(n, 2)*bzg(j0+2, k0, l))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l) + sx0(n, 0)*bzg(j0, k0+1, l) + sx0(n,  &
        1)*bzg(j0+1, k0+1, l) + sx0(n, 2)*bzg(j0+2, k0+1, l))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l) + sx0(n, 0)*bzg(j0, k0+2, l) + sx0(n,  &
        1)*bzg(j0+1, k0+2, l) + sx0(n, 2)*bzg(j0+2, k0+2, l))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+1) + sx0(n, 0)*bzg(j0, k0-1, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+1) + sx0(n, 2)*bzg(j0+2, k0-1, l+1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+1) + sx0(n, 0)*bzg(j0, k0, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0, l+1) + sx0(n, 2)*bzg(j0+2, k0, l+1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+1) + sx0(n, 0)*bzg(j0, k0+1, l+1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+1) + sx0(n, 2)*bzg(j0+2, k0+1, l+1))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l+1) + sx0(n, 0)*bzg(j0, k0+2, l+1) +     &
        sx0(n, 1)*bzg(j0+1, k0+2, l+1) + sx0(n, 2)*bzg(j0+2, k0+2, l+1))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+2) + sx0(n, 0)*bzg(j0, k0-1, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+2) + sx0(n, 2)*bzg(j0+2, k0-1, l+2))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+2) + sx0(n, 0)*bzg(j0, k0, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0, l+2) + sx0(n, 2)*bzg(j0+2, k0, l+2))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+2) + sx0(n, 0)*bzg(j0, k0+1, l+2) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+2) + sx0(n, 2)*bzg(j0+2, k0+1, l+2))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l+2) + sx0(n, 0)*bzg(j0, k0+2, l+2) +     &
        sx0(n, 1)*bzg(j0+1, k0+2, l+2) + sx0(n, 2)*bzg(j0+2, k0+2, l+2))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, 2)

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
!
!> @brief
!> Vectorized field gathering at order 3 with gathering of E and B merged in a single loop
!
!> @details
!> This function is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!> Last modified 12/02/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] ex, ey, ez particle electric field arrays
!> @param[inout] bx, by, bz particle magnetic field arrays
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vec_3_3_3(np, xp, yp, zp, ex, ey, ez, bx, by,    &
  bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg,    &
  ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, j, k, l
  INTEGER(idp)                           :: jj, kk, ll
  INTEGER(idp)                           :: j0, k0, l0
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq,     &
  ozint, oxintsq, oyintsq, ozintsq
  INTEGER(isp)                           :: nn, n
  REAL(num), DIMENSION(lvect, -1:2)       :: sx, sx0
  REAL(num), DIMENSION(lvect, -1:2)       :: sy, sy0
  REAL(num), DIMENSION(lvect, -1:2)       :: sz, sz0
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

  IF (l_lower_order_in_v ) THEN

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
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        yintsq = yint*yint
        sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, -1)*exg(j0-1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, -1)*exg(j0, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, -1)*exg(j0+1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, -1)*exg(j0-1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, -1)*exg(j0, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, -1)*exg(j0+1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, -1)*exg(j0-1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, -1)*exg(j0, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, -1)*exg(j0+1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, -1)*exg(j0-1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, -1)*exg(j0, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, -1)*exg(j0+1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 0)*exg(j0-1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 0)*exg(j0, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 0)*exg(j0+1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 0)*exg(j0-1, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0+1, k, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 0)*exg(j0-1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0+1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 0)*exg(j0-1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 0)*exg(j0, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 0)*exg(j0+1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 1)*exg(j0-1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 1)*exg(j0, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 1)*exg(j0+1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 1)*exg(j0-1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0+1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 1)*exg(j0-1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0+1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 1)*exg(j0-1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 1)*exg(j0, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 1)*exg(j0+1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 2)*exg(j0-1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 2)*exg(j0, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 2)*exg(j0+1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 2)*exg(j0-1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 2)*exg(j0, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 2)*exg(j0+1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 2)*exg(j0-1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 2)*exg(j0, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 2)*exg(j0+1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 2)*exg(j0-1, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 2)*exg(j0, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 2)*exg(j0+1, k+2, l+2)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, -1)*eyg(j-1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, -1)*eyg(j, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, -1)*eyg(j+1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, -1)*eyg(j+2, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, -1)*eyg(j-1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, -1)*eyg(j, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, -1)*eyg(j+1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, -1)*eyg(j+2, k0, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, -1)*eyg(j-1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, -1)*eyg(j, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, -1)*eyg(j+1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, -1)*eyg(j+2, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 0)*eyg(j-1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 0)*eyg(j, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 0)*eyg(j+1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 0)*eyg(j+2, k0-1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 0)*eyg(j-1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 0)*eyg(j+2, k0, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 0)*eyg(j-1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j+1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 0)*eyg(j+2, k0+1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 1)*eyg(j-1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 1)*eyg(j, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 1)*eyg(j+1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 1)*eyg(j+2, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 1)*eyg(j-1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 1)*eyg(j+2, k0, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 1)*eyg(j-1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j+1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 1)*eyg(j+2, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 2)*eyg(j-1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 2)*eyg(j, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 2)*eyg(j+1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 2)*eyg(j+2, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 2)*eyg(j-1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 2)*eyg(j, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 2)*eyg(j+1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 2)*eyg(j+2, k0, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 2)*eyg(j-1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 2)*eyg(j, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 2)*eyg(j+1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 2)*eyg(j+2, k0+1, l+2)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, -1)*ezg(j-1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, -1)*ezg(j, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, -1)*ezg(j+1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, -1)*ezg(j+2, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, -1)*ezg(j-1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, -1)*ezg(j, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, -1)*ezg(j+1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, -1)*ezg(j+2, k, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, -1)*ezg(j-1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, -1)*ezg(j, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, -1)*ezg(j+1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, -1)*ezg(j+2, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, -1)*ezg(j-1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, -1)*ezg(j, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, -1)*ezg(j+1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, -1)*ezg(j+2, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 0)*ezg(j-1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 0)*ezg(j, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 0)*ezg(j+1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 0)*ezg(j+2, k-1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 0)*ezg(j-1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 0)*ezg(j+2, k, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 0)*ezg(j-1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 0)*ezg(j+2, k+1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 0)*ezg(j-1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 0)*ezg(j, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 0)*ezg(j+1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 0)*ezg(j+2, k+2, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 1)*ezg(j-1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 1)*ezg(j, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 1)*ezg(j+1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 1)*ezg(j+2, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 1)*ezg(j-1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j+1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 1)*ezg(j+2, k, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 1)*ezg(j-1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j+1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 1)*ezg(j+2, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 1)*ezg(j-1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 1)*ezg(j, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 1)*ezg(j+1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 1)*ezg(j+2, k+2, l0+1)

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, -1)*bxg(j-1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, -1)*bxg(j, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, -1)*bxg(j+1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, -1)*bxg(j+2, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, -1)*bxg(j-1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, -1)*bxg(j, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, -1)*bxg(j+1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, -1)*bxg(j+2, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, -1)*bxg(j-1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, -1)*bxg(j, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, -1)*bxg(j+1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, -1)*bxg(j+2, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 0)*bxg(j-1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 0)*bxg(j, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 0)*bxg(j+1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 0)*bxg(j+2, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 0)*bxg(j-1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 0)*bxg(j+2, k0, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 0)*bxg(j-1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j+1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 0)*bxg(j+2, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 1)*bxg(j-1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 1)*bxg(j, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 1)*bxg(j+1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 1)*bxg(j+2, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 1)*bxg(j-1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j+1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 1)*bxg(j+2, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 1)*bxg(j-1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j+1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 1)*bxg(j+2, k0+1, l0+1)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, -1)*byg(j0-1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, -1)*byg(j0, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, -1)*byg(j0+1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, -1)*byg(j0-1, k, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, -1)*byg(j0, k, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, -1)*byg(j0+1, k, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, -1)*byg(j0-1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, -1)*byg(j0, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, -1)*byg(j0+1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, -1)*byg(j0-1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, -1)*byg(j0, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, -1)*byg(j0+1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 0)*byg(j0-1, k-1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 0)*byg(j0, k-1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 0)*byg(j0+1, k-1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 0)*byg(j0-1, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0+1, k, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 0)*byg(j0-1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0+1, k+1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 0)*byg(j0-1, k+2, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 0)*byg(j0, k+2, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 0)*byg(j0+1, k+2, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 1)*byg(j0-1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 1)*byg(j0, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 1)*byg(j0+1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 1)*byg(j0-1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0, k, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0+1, k, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 1)*byg(j0-1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0+1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 1)*byg(j0-1, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 1)*byg(j0, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 1)*byg(j0+1, k+2, l0+1)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, -1)*bzg(j0-1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, -1)*bzg(j0, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, -1)*bzg(j0+1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, -1)*bzg(j0-1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, -1)*bzg(j0, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, -1)*bzg(j0+1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, -1)*bzg(j0-1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, -1)*bzg(j0, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, -1)*bzg(j0+1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 0)*bzg(j0-1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 0)*bzg(j0, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 0)*bzg(j0+1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 0)*bzg(j0-1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0+1, k0, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 0)*bzg(j0-1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0+1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 1)*bzg(j0-1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 1)*bzg(j0, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 1)*bzg(j0+1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 1)*bzg(j0-1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0+1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 1)*bzg(j0-1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0+1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 2)*bzg(j0-1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 2)*bzg(j0, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 2)*bzg(j0+1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 2)*bzg(j0-1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 2)*bzg(j0, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 2)*bzg(j0+1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 2)*bzg(j0-1, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 2)*bzg(j0, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 2)*bzg(j0+1, k0+1, l+2)

      ENDDO
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
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(n, -1) = onesixth*oyintsq*oyint
        sy0(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, -1)*exg(j0-1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, -1)*exg(j0, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, -1)*exg(j0+1, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, -1)*exg(j0+2, k-1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, -1)*exg(j0-1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, -1)*exg(j0, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, -1)*exg(j0+1, k, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, -1)*exg(j0+2, k, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, -1)*exg(j0-1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, -1)*exg(j0, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, -1)*exg(j0+1, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, -1)*exg(j0+2, k+1, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, -1)*exg(j0-1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, -1)*exg(j0, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, -1)*exg(j0+1, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, -1)*exg(j0+2, k+2, l-1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 0)*exg(j0-1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 0)*exg(j0, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 0)*exg(j0+1, k-1, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, 0)*exg(j0+2, k-1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 0)*exg(j0-1, k, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0, k, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0+1, k, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, 0)*exg(j0+2, k, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 0)*exg(j0-1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0+1, k+1, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, 0)*exg(j0+2, k+1, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 0)*exg(j0-1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 0)*exg(j0, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 0)*exg(j0+1, k+2, l)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, 0)*exg(j0+2, k+2, l)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 1)*exg(j0-1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 1)*exg(j0, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 1)*exg(j0+1, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, 1)*exg(j0+2, k-1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 1)*exg(j0-1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0+1, k, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, 1)*exg(j0+2, k, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 1)*exg(j0-1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0+1, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, 1)*exg(j0+2, k+1, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 1)*exg(j0-1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 1)*exg(j0, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 1)*exg(j0+1, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, 1)*exg(j0+2, k+2, l+1)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 2)*exg(j0-1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 2)*exg(j0, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 2)*exg(j0+1, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, -1)*sz(n, 2)*exg(j0+2, k-1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 2)*exg(j0-1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 2)*exg(j0, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 2)*exg(j0+1, k, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 0)*sz(n, 2)*exg(j0+2, k, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 2)*exg(j0-1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 2)*exg(j0, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 2)*exg(j0+1, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 1)*sz(n, 2)*exg(j0+2, k+1, l+2)
        ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 2)*exg(j0-1, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 2)*exg(j0, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 2)*exg(j0+1, k+2, l+2)
        ex(nn) = ex(nn) + sx0(n, 2)*sy(n, 2)*sz(n, 2)*exg(j0+2, k+2, l+2)

        ! Compute Ey on particle
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, -1)*eyg(j-1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, -1)*eyg(j, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, -1)*eyg(j+1, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, -1)*eyg(j+2, k0-1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, -1)*eyg(j-1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, -1)*eyg(j, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, -1)*eyg(j+1, k0, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, -1)*eyg(j+2, k0, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, -1)*eyg(j-1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, -1)*eyg(j, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, -1)*eyg(j+1, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, -1)*eyg(j+2, k0+1, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, -1)*eyg(j-1, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, -1)*eyg(j, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, -1)*eyg(j+1, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, -1)*eyg(j+2, k0+2, l-1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 0)*eyg(j-1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 0)*eyg(j, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 0)*eyg(j+1, k0-1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 0)*eyg(j+2, k0-1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 0)*eyg(j-1, k0, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j, k0, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j+1, k0, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 0)*eyg(j+2, k0, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 0)*eyg(j-1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j+1, k0+1, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 0)*eyg(j+2, k0+1, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, 0)*eyg(j-1, k0+2, l)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, 0)*eyg(j, k0+2, l)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, 0)*eyg(j+1, k0+2, l)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, 0)*eyg(j+2, k0+2, l)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 1)*eyg(j-1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 1)*eyg(j, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 1)*eyg(j+1, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 1)*eyg(j+2, k0-1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 1)*eyg(j-1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j+1, k0, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 1)*eyg(j+2, k0, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 1)*eyg(j-1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j+1, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 1)*eyg(j+2, k0+1, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, 1)*eyg(j-1, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, 1)*eyg(j, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, 1)*eyg(j+1, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, 1)*eyg(j+2, k0+2, l+1)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 2)*eyg(j-1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 2)*eyg(j, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 2)*eyg(j+1, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 2)*eyg(j+2, k0-1, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 2)*eyg(j-1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 2)*eyg(j, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 2)*eyg(j+1, k0, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 2)*eyg(j+2, k0, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 2)*eyg(j-1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 2)*eyg(j, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 2)*eyg(j+1, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 2)*eyg(j+2, k0+1, l+2)
        ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 2)*sz(n, 2)*eyg(j-1, k0+2, l+2)
        ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 2)*sz(n, 2)*eyg(j, k0+2, l+2)
        ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 2)*sz(n, 2)*eyg(j+1, k0+2, l+2)
        ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 2)*sz(n, 2)*eyg(j+2, k0+2, l+2)

        ! Compute Ez on particle
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, -1)*ezg(j-1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, -1)*ezg(j, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, -1)*ezg(j+1, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, -1)*ezg(j+2, k-1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, -1)*ezg(j-1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, -1)*ezg(j, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, -1)*ezg(j+1, k, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, -1)*ezg(j+2, k, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, -1)*ezg(j-1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, -1)*ezg(j, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, -1)*ezg(j+1, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, -1)*ezg(j+2, k+1, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, -1)*ezg(j-1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, -1)*ezg(j, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, -1)*ezg(j+1, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, -1)*ezg(j+2, k+2, l0-1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 0)*ezg(j-1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 0)*ezg(j, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 0)*ezg(j+1, k-1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 0)*ezg(j+2, k-1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 0)*ezg(j-1, k, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j, k, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j+1, k, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 0)*ezg(j+2, k, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 0)*ezg(j-1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j+1, k+1, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 0)*ezg(j+2, k+1, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 0)*ezg(j-1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 0)*ezg(j, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 0)*ezg(j+1, k+2, l0)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 0)*ezg(j+2, k+2, l0)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 1)*ezg(j-1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 1)*ezg(j, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 1)*ezg(j+1, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 1)*ezg(j+2, k-1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 1)*ezg(j-1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j+1, k, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 1)*ezg(j+2, k, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 1)*ezg(j-1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j+1, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 1)*ezg(j+2, k+1, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 1)*ezg(j-1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 1)*ezg(j, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 1)*ezg(j+1, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 1)*ezg(j+2, k+2, l0+1)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 2)*ezg(j-1, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 2)*ezg(j, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 2)*ezg(j+1, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 2)*ezg(j+2, k-1, l0+2)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 2)*ezg(j-1, k, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 2)*ezg(j, k, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 2)*ezg(j+1, k, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 2)*ezg(j+2, k, l0+2)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 2)*ezg(j-1, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 2)*ezg(j, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 2)*ezg(j+1, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 2)*ezg(j+2, k+1, l0+2)
        ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 2)*ezg(j-1, k+2, l0+2)
        ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 2)*ezg(j, k+2, l0+2)
        ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 2)*ezg(j+1, k+2, l0+2)
        ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 2)*ezg(j+2, k+2, l0+2)

        ! Compute Bx on particle
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, -1)*bxg(j-1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, -1)*bxg(j, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, -1)*bxg(j+1, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, -1)*bxg(j+2, k0-1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, -1)*bxg(j-1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, -1)*bxg(j, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, -1)*bxg(j+1, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, -1)*bxg(j+2, k0, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, -1)*bxg(j-1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, -1)*bxg(j, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, -1)*bxg(j+1, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, -1)*bxg(j+2, k0+1, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, -1)*bxg(j-1, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, -1)*bxg(j, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, -1)*bxg(j+1, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, -1)*bxg(j+2, k0+2, l0-1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 0)*bxg(j-1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 0)*bxg(j, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 0)*bxg(j+1, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 0)*bxg(j+2, k0-1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 0)*bxg(j-1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j, k0, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j+1, k0, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 0)*bxg(j+2, k0, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 0)*bxg(j-1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j+1, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 0)*bxg(j+2, k0+1, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, 0)*bxg(j-1, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, 0)*bxg(j, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, 0)*bxg(j+1, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, 0)*bxg(j+2, k0+2, l0)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 1)*bxg(j-1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 1)*bxg(j, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 1)*bxg(j+1, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 1)*bxg(j+2, k0-1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 1)*bxg(j-1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j+1, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 1)*bxg(j+2, k0, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 1)*bxg(j-1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j+1, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 1)*bxg(j+2, k0+1, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, 1)*bxg(j-1, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, 1)*bxg(j, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, 1)*bxg(j+1, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, 1)*bxg(j+2, k0+2, l0+1)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 2)*bxg(j-1, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 2)*bxg(j, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 2)*bxg(j+1, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 2)*bxg(j+2, k0-1, l0+2)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 2)*bxg(j-1, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 2)*bxg(j, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 2)*bxg(j+1, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 2)*bxg(j+2, k0, l0+2)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 2)*bxg(j-1, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 2)*bxg(j, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 2)*bxg(j+1, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 2)*bxg(j+2, k0+1, l0+2)
        bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 2)*sz0(n, 2)*bxg(j-1, k0+2, l0+2)
        bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 2)*sz0(n, 2)*bxg(j, k0+2, l0+2)
        bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 2)*sz0(n, 2)*bxg(j+1, k0+2, l0+2)
        bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 2)*sz0(n, 2)*bxg(j+2, k0+2, l0+2)

        ! Compute By on particle
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, -1)*byg(j0-1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, -1)*byg(j0, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, -1)*byg(j0+1, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, -1)*byg(j0+2, k-1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, -1)*byg(j0-1, k, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, -1)*byg(j0, k, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, -1)*byg(j0+1, k, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, -1)*byg(j0+2, k, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, -1)*byg(j0-1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, -1)*byg(j0, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, -1)*byg(j0+1, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, -1)*byg(j0+2, k+1, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, -1)*byg(j0-1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, -1)*byg(j0, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, -1)*byg(j0+1, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, -1)*byg(j0+2, k+2, l0-1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 0)*byg(j0-1, k-1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 0)*byg(j0, k-1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 0)*byg(j0+1, k-1, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, 0)*byg(j0+2, k-1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 0)*byg(j0-1, k, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0, k, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0+1, k, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, 0)*byg(j0+2, k, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 0)*byg(j0-1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0, k+1, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0+1, k+1, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, 0)*byg(j0+2, k+1, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 0)*byg(j0-1, k+2, l0)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 0)*byg(j0, k+2, l0)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 0)*byg(j0+1, k+2, l0)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, 0)*byg(j0+2, k+2, l0)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 1)*byg(j0-1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 1)*byg(j0, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 1)*byg(j0+1, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, 1)*byg(j0+2, k-1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 1)*byg(j0-1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0, k, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0+1, k, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, 1)*byg(j0+2, k, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 1)*byg(j0-1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0+1, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, 1)*byg(j0+2, k+1, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 1)*byg(j0-1, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 1)*byg(j0, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 1)*byg(j0+1, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, 1)*byg(j0+2, k+2, l0+1)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 2)*byg(j0-1, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 2)*byg(j0, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 2)*byg(j0+1, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, -1)*sz0(n, 2)*byg(j0+2, k-1, l0+2)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 2)*byg(j0-1, k, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 2)*byg(j0, k, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 2)*byg(j0+1, k, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 0)*sz0(n, 2)*byg(j0+2, k, l0+2)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 2)*byg(j0-1, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 2)*byg(j0, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 2)*byg(j0+1, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 1)*sz0(n, 2)*byg(j0+2, k+1, l0+2)
        by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 2)*byg(j0-1, k+2, l0+2)
        by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 2)*byg(j0, k+2, l0+2)
        by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 2)*byg(j0+1, k+2, l0+2)
        by(nn) = by(nn) + sx0(n, 2)*sy(n, 2)*sz0(n, 2)*byg(j0+2, k+2, l0+2)

        ! Compute Bz on particle
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, -1)*bzg(j0-1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, -1)*bzg(j0, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, -1)*bzg(j0+1, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, -1)*bzg(j0+2, k0-1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, -1)*bzg(j0-1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, -1)*bzg(j0, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, -1)*bzg(j0+1, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, -1)*bzg(j0+2, k0, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, -1)*bzg(j0-1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, -1)*bzg(j0, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, -1)*bzg(j0+1, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, -1)*bzg(j0+2, k0+1, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, -1)*bzg(j0-1, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, -1)*bzg(j0, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, -1)*bzg(j0+1, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, -1)*bzg(j0+2, k0+2, l-1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 0)*bzg(j0-1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 0)*bzg(j0, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 0)*bzg(j0+1, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, 0)*bzg(j0+2, k0-1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 0)*bzg(j0-1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0, k0, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0+1, k0, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, 0)*bzg(j0+2, k0, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 0)*bzg(j0-1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0+1, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, 0)*bzg(j0+2, k0+1, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, 0)*bzg(j0-1, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, 0)*bzg(j0, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, 0)*bzg(j0+1, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, 0)*bzg(j0+2, k0+2, l)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 1)*bzg(j0-1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 1)*bzg(j0, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 1)*bzg(j0+1, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, 1)*bzg(j0+2, k0-1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 1)*bzg(j0-1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0+1, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, 1)*bzg(j0+2, k0, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 1)*bzg(j0-1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0+1, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, 1)*bzg(j0+2, k0+1, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, 1)*bzg(j0-1, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, 1)*bzg(j0, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, 1)*bzg(j0+1, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, 1)*bzg(j0+2, k0+2, l+1)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 2)*bzg(j0-1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 2)*bzg(j0, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 2)*bzg(j0+1, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, -1)*sz(n, 2)*bzg(j0+2, k0-1, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 2)*bzg(j0-1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 2)*bzg(j0, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 2)*bzg(j0+1, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 0)*sz(n, 2)*bzg(j0+2, k0, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 2)*bzg(j0-1, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 2)*bzg(j0, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 2)*bzg(j0+1, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 1)*sz(n, 2)*bzg(j0+2, k0+1, l+2)
        bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 2)*sz(n, 2)*bzg(j0-1, k0+2, l+2)
        bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 2)*sz(n, 2)*bzg(j0, k0+2, l+2)
        bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 2)*sz(n, 2)*bzg(j0+1, k0+2, l+2)
        bz(nn) = bz(nn) + sx0(n, 2)*sy0(n, 2)*sz(n, 2)*bzg(j0+2, k0+2, l+2)
      ENDDO
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

! ________________________________________________________________________________________
!
!> @brief
!> Vectorized field gathering at order 3 with gathering of E and B merged in a single loop
!
!> @details
!> This function is vectorized and is the version 2 using a more optimized
!> gathering method more efficient for the vectorization.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 12/03/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] ex, ey, ez particle electric field arrays
!> @param[inout] bx, by, bz particle magnetic field arrays
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] nx, ny, nz number of grid points in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells in each direction
!> @param[in] exg, eyg, ezg electric field grid
!> @param[in] bxg, byg, bzg magnetic field grid
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
! ________________________________________________________________________________________
SUBROUTINE geteb3d_energy_conserving_vecV2_3_3_3(np, xp, yp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg, eyg,    &
  ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, j, k, l
  INTEGER(idp)                           :: j0, k0, l0
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq
  REAL(num)                              :: ozint, oxintsq, oyintsq, ozintsq
  REAL(num)                              :: a
  INTEGER(isp)                           :: nn, n
  REAL(num), DIMENSION(lvect, -1:2)       :: sx, sx0
  REAL(num), DIMENSION(lvect, -1:2)       :: sy, sy0
  REAL(num), DIMENSION(lvect, -1:2)       :: sz, sz0
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

  IF (l_lower_order_in_v ) THEN

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
      !!DIR PREFETCH ex:1:1
      !!DIR PREFETCH ey:1:1
      !!DIR PREFETCH ez:1:1
      !!DIR PREFETCH bx:1:1
      !!DIR PREFETCH by:1:1
      !!DIR PREFETCH bz:1:1
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
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
        sx0(n, 0) = 0.75_num-xintsq
        sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

        yintsq = yint*yint
        sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
        sy0(n, 0) = 0.75_num-yintsq
        sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

        zintsq = zint*zint
        sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
        sz0(n, 0) = 0.75_num-zintsq
        sz0(n, 1) = 0.5_num*(0.5_num+zint)**2

        ! Compute Ex on particle
        a = (sx0(n, -1)*exg(j0-1, k-1, l-1) + sx0(n, 0)*exg(j0, k-1, l-1) + sx0(n,    &
        1)*exg(j0+1, k-1, l-1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l-1) + sx0(n, 0)*exg(j0, k, l-1) + sx0(n,    &
        1)*exg(j0+1, k, l-1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l-1) + sx0(n, 0)*exg(j0, k+1, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l-1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l-1) + sx0(n, 0)*exg(j0, k+2, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l-1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l) + sx0(n, 0)*exg(j0, k-1, l) + sx0(n,        &
        1)*exg(j0+1, k-1, l))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l) + sx0(n, 0)*exg(j0, k, l) + sx0(n,        &
        1)*exg(j0+1, k, l))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l) + sx0(n, 0)*exg(j0, k+1, l) + sx0(n,    &
        1)*exg(j0+1, k+1, l))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l) + sx0(n, 0)*exg(j0, k+2, l) + sx0(n,    &
        1)*exg(j0+1, k+2, l))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+1) + sx0(n, 0)*exg(j0, k-1, l+1) + sx0(n,    &
        1)*exg(j0+1, k-1, l+1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+1) + sx0(n, 0)*exg(j0, k, l+1) + sx0(n,    &
        1)*exg(j0+1, k, l+1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+1) + sx0(n, 0)*exg(j0, k+1, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+1) + sx0(n, 0)*exg(j0, k+2, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+2) + sx0(n, 0)*exg(j0, k-1, l+2) + sx0(n,    &
        1)*exg(j0+1, k-1, l+2))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+2) + sx0(n, 0)*exg(j0, k, l+2) + sx0(n,    &
        1)*exg(j0+1, k, l+2))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+2) + sx0(n, 0)*exg(j0, k+1, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+2))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+2) + sx0(n, 0)*exg(j0, k+2, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+2))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 2)

        ! Compute Ey on particle
        a = (sx(n, -1)*eyg(j-1, k0-1, l-1) + sx(n, 0)*eyg(j, k0-1, l-1) + sx(n,       &
        1)*eyg(j+1, k0-1, l-1) + sx(n, 2)*eyg(j+2, k0-1, l-1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l-1) + sx(n, 0)*eyg(j, k0, l-1) + sx(n,       &
        1)*eyg(j+1, k0, l-1) + sx(n, 2)*eyg(j+2, k0, l-1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l-1) + sx(n, 0)*eyg(j, k0+1, l-1) + sx(n,   &
        1)*eyg(j+1, k0+1, l-1) + sx(n, 2)*eyg(j+2, k0+1, l-1))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, -1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l) + sx(n, 0)*eyg(j, k0-1, l) + sx(n,           &
        1)*eyg(j+1, k0-1, l) + sx(n, 2)*eyg(j+2, k0-1, l))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l) + sx(n, 0)*eyg(j, k0, l) + sx(n,           &
        1)*eyg(j+1, k0, l) + sx(n, 2)*eyg(j+2, k0, l))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l) + sx(n, 0)*eyg(j, k0+1, l) + sx(n,       &
        1)*eyg(j+1, k0+1, l) + sx(n, 2)*eyg(j+2, k0+1, l))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+1) + sx(n, 0)*eyg(j, k0-1, l+1) + sx(n,       &
        1)*eyg(j+1, k0-1, l+1) + sx(n, 2)*eyg(j+2, k0-1, l+1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+1) + sx(n, 0)*eyg(j, k0, l+1) + sx(n,       &
        1)*eyg(j+1, k0, l+1) + sx(n, 2)*eyg(j+2, k0, l+1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+1) + sx(n, 0)*eyg(j, k0+1, l+1) + sx(n,   &
        1)*eyg(j+1, k0+1, l+1) + sx(n, 2)*eyg(j+2, k0+1, l+1))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+2) + sx(n, 0)*eyg(j, k0-1, l+2) + sx(n,       &
        1)*eyg(j+1, k0-1, l+2) + sx(n, 2)*eyg(j+2, k0-1, l+2))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+2) + sx(n, 0)*eyg(j, k0, l+2) + sx(n,       &
        1)*eyg(j+1, k0, l+2) + sx(n, 2)*eyg(j+2, k0, l+2))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+2) + sx(n, 0)*eyg(j, k0+1, l+2) + sx(n,   &
        1)*eyg(j+1, k0+1, l+2) + sx(n, 2)*eyg(j+2, k0+1, l+2))*sy0(n, 1)
        ey(nn) = ey(nn) + a*sz(n, 2)

        ! Compute Ez on particle
        a = (sx(n, -1)*ezg(j-1, k-1, l0-1) + sx(n, 0)*ezg(j, k-1, l0-1) + sx(n,       &
        1)*ezg(j+1, k-1, l0-1) + sx(n, 2)*ezg(j+2, k-1, l0-1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0-1) + sx(n, 0)*ezg(j, k, l0-1) + sx(n,       &
        1)*ezg(j+1, k, l0-1) + sx(n, 2)*ezg(j+2, k, l0-1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0-1) + sx(n, 0)*ezg(j, k+1, l0-1) + sx(n,   &
        1)*ezg(j+1, k+1, l0-1) + sx(n, 2)*ezg(j+2, k+1, l0-1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0-1) + sx(n, 0)*ezg(j, k+2, l0-1) + sx(n,   &
        1)*ezg(j+1, k+2, l0-1) + sx(n, 2)*ezg(j+2, k+2, l0-1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*ezg(j-1, k-1, l0) + sx(n, 0)*ezg(j, k-1, l0) + sx(n,           &
        1)*ezg(j+1, k-1, l0) + sx(n, 2)*ezg(j+2, k-1, l0))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0) + sx(n, 0)*ezg(j, k, l0) + sx(n,           &
        1)*ezg(j+1, k, l0) + sx(n, 2)*ezg(j+2, k, l0))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0) + sx(n, 0)*ezg(j, k+1, l0) + sx(n,       &
        1)*ezg(j+1, k+1, l0) + sx(n, 2)*ezg(j+2, k+1, l0))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0) + sx(n, 0)*ezg(j, k+2, l0) + sx(n,       &
        1)*ezg(j+1, k+2, l0) + sx(n, 2)*ezg(j+2, k+2, l0))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*ezg(j-1, k-1, l0+1) + sx(n, 0)*ezg(j, k-1, l0+1) + sx(n,       &
        1)*ezg(j+1, k-1, l0+1) + sx(n, 2)*ezg(j+2, k-1, l0+1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0+1) + sx(n, 0)*ezg(j, k, l0+1) + sx(n,       &
        1)*ezg(j+1, k, l0+1) + sx(n, 2)*ezg(j+2, k, l0+1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0+1) + sx(n, 0)*ezg(j, k+1, l0+1) + sx(n,   &
        1)*ezg(j+1, k+1, l0+1) + sx(n, 2)*ezg(j+2, k+1, l0+1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0+1) + sx(n, 0)*ezg(j, k+2, l0+1) + sx(n,   &
        1)*ezg(j+1, k+2, l0+1) + sx(n, 2)*ezg(j+2, k+2, l0+1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 1)

        ! Compute Bx on particle
        a = (sx(n, -1)*bxg(j-1, k0-1, l0-1) + sx(n, 0)*bxg(j, k0-1, l0-1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0-1) + sx(n, 2)*bxg(j+2, k0-1, l0-1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0-1) + sx(n, 0)*bxg(j, k0, l0-1) + sx(n,     &
        1)*bxg(j+1, k0, l0-1) + sx(n, 2)*bxg(j+2, k0, l0-1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0-1) + sx(n, 0)*bxg(j, k0+1, l0-1) + sx(n, &
        1)*bxg(j+1, k0+1, l0-1) + sx(n, 2)*bxg(j+2, k0+1, l0-1))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0) + sx(n, 0)*bxg(j, k0-1, l0) + sx(n,         &
        1)*bxg(j+1, k0-1, l0) + sx(n, 2)*bxg(j+2, k0-1, l0))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0) + sx(n, 0)*bxg(j, k0, l0) + sx(n,         &
        1)*bxg(j+1, k0, l0) + sx(n, 2)*bxg(j+2, k0, l0))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0) + sx(n, 0)*bxg(j, k0+1, l0) + sx(n,     &
        1)*bxg(j+1, k0+1, l0) + sx(n, 2)*bxg(j+2, k0+1, l0))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0+1) + sx(n, 0)*bxg(j, k0-1, l0+1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0+1) + sx(n, 2)*bxg(j+2, k0-1, l0+1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0+1) + sx(n, 0)*bxg(j, k0, l0+1) + sx(n,     &
        1)*bxg(j+1, k0, l0+1) + sx(n, 2)*bxg(j+2, k0, l0+1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0+1) + sx(n, 0)*bxg(j, k0+1, l0+1) + sx(n, &
        1)*bxg(j+1, k0+1, l0+1) + sx(n, 2)*bxg(j+2, k0+1, l0+1))*sy0(n, 1)
        bx(nn) = bx(nn) + a*sz0(n, 1)

        ! Compute By on particle
        a = (sx0(n, -1)*byg(j0-1, k-1, l0-1) + sx0(n, 0)*byg(j0, k-1, l0-1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0-1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0-1) + sx0(n, 0)*byg(j0, k, l0-1) + sx0(n,  &
        1)*byg(j0+1, k, l0-1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0-1) + sx0(n, 0)*byg(j0, k+1, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0-1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0-1) + sx0(n, 0)*byg(j0, k+2, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0-1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, -1)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0) + sx0(n, 0)*byg(j0, k-1, l0) + sx0(n,      &
        1)*byg(j0+1, k-1, l0))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0) + sx0(n, 0)*byg(j0, k, l0) + sx0(n,      &
        1)*byg(j0+1, k, l0))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0) + sx0(n, 0)*byg(j0, k+1, l0) + sx0(n,  &
        1)*byg(j0+1, k+1, l0))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0) + sx0(n, 0)*byg(j0, k+2, l0) + sx0(n,  &
        1)*byg(j0+1, k+2, l0))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0+1) + sx0(n, 0)*byg(j0, k-1, l0+1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0+1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0+1) + sx0(n, 0)*byg(j0, k, l0+1) + sx0(n,  &
        1)*byg(j0+1, k, l0+1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0+1) + sx0(n, 0)*byg(j0, k+1, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0+1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0+1) + sx0(n, 0)*byg(j0, k+2, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0+1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 1)

        ! Compute Bz on particle
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l-1) + sx0(n, 0)*bzg(j0, k0-1, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l-1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l-1) + sx0(n, 0)*bzg(j0, k0, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0, l-1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l-1) + sx0(n, 0)*bzg(j0, k0+1, l-1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l-1))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l) + sx0(n, 0)*bzg(j0, k0-1, l) + sx0(n,      &
        1)*bzg(j0+1, k0-1, l))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l) + sx0(n, 0)*bzg(j0, k0, l) + sx0(n,      &
        1)*bzg(j0+1, k0, l))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l) + sx0(n, 0)*bzg(j0, k0+1, l) + sx0(n,  &
        1)*bzg(j0+1, k0+1, l))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+1) + sx0(n, 0)*bzg(j0, k0-1, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+1) + sx0(n, 0)*bzg(j0, k0, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0, l+1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+1) + sx0(n, 0)*bzg(j0, k0+1, l+1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+1))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+2) + sx0(n, 0)*bzg(j0, k0-1, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+2))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+2) + sx0(n, 0)*bzg(j0, k0, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0, l+2))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+2) + sx0(n, 0)*bzg(j0, k0+1, l+2) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+2))*sy0(n, 1)
        bz(nn) = bz(nn) + a*sz(n, 2)

      ENDDO
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
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(n, -1) = onesixth*oxintsq*oxint
        sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(n, -1) = onesixth*oyintsq*oyint
        sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(n, -1) = onesixth*ozintsq*ozint
        sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz(n, 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(n, -1) = onesixth*oxintsq*oxint
        sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0(n, 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(n, -1) = onesixth*oyintsq*oyint
        sy0(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0(n, 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(n, -1) = onesixth*ozintsq*ozint
        sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0(n, 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        a = (sx0(n, -1)*exg(j0-1, k-1, l-1) + sx0(n, 0)*exg(j0, k-1, l-1) + sx0(n,    &
        1)*exg(j0+1, k-1, l-1) + sx0(n, 2)*exg(j0+2, k-1, l-1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l-1) + sx0(n, 0)*exg(j0, k, l-1) + sx0(n,    &
        1)*exg(j0+1, k, l-1) + sx0(n, 2)*exg(j0+2, k, l-1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l-1) + sx0(n, 0)*exg(j0, k+1, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l-1) + sx0(n, 2)*exg(j0+2, k+1, l-1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l-1) + sx0(n, 0)*exg(j0, k+2, l-1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l-1) + sx0(n, 2)*exg(j0+2, k+2, l-1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l) + sx0(n, 0)*exg(j0, k-1, l) + sx0(n,        &
        1)*exg(j0+1, k-1, l) + sx0(n, 2)*exg(j0+2, k-1, l))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l) + sx0(n, 0)*exg(j0, k, l) + sx0(n,        &
        1)*exg(j0+1, k, l) + sx0(n, 2)*exg(j0+2, k, l))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l) + sx0(n, 0)*exg(j0, k+1, l) + sx0(n,    &
        1)*exg(j0+1, k+1, l) + sx0(n, 2)*exg(j0+2, k+1, l))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l) + sx0(n, 0)*exg(j0, k+2, l) + sx0(n,    &
        1)*exg(j0+1, k+2, l) + sx0(n, 2)*exg(j0+2, k+2, l))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+1) + sx0(n, 0)*exg(j0, k-1, l+1) + sx0(n,    &
        1)*exg(j0+1, k-1, l+1) + sx0(n, 2)*exg(j0+2, k-1, l+1))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+1) + sx0(n, 0)*exg(j0, k, l+1) + sx0(n,    &
        1)*exg(j0+1, k, l+1) + sx0(n, 2)*exg(j0+2, k, l+1))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+1) + sx0(n, 0)*exg(j0, k+1, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+1) + sx0(n, 2)*exg(j0+2, k+1, l+1))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+1) + sx0(n, 0)*exg(j0, k+2, l+1) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+1) + sx0(n, 2)*exg(j0+2, k+2, l+1))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*exg(j0-1, k-1, l+2) + sx0(n, 0)*exg(j0, k-1, l+2) + sx0(n,    &
        1)*exg(j0+1, k-1, l+2) + sx0(n, 2)*exg(j0+2, k-1, l+2))*sy(n, -1)
        a = a + (sx0(n, -1)*exg(j0-1, k, l+2) + sx0(n, 0)*exg(j0, k, l+2) + sx0(n,    &
        1)*exg(j0+1, k, l+2) + sx0(n, 2)*exg(j0+2, k, l+2))*sy(n, 0)
        a = a + (sx0(n, -1)*exg(j0-1, k+1, l+2) + sx0(n, 0)*exg(j0, k+1, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+1, l+2) + sx0(n, 2)*exg(j0+2, k+1, l+2))*sy(n, 1)
        a = a + (sx0(n, -1)*exg(j0-1, k+2, l+2) + sx0(n, 0)*exg(j0, k+2, l+2) +       &
        sx0(n, 1)*exg(j0+1, k+2, l+2) + sx0(n, 2)*exg(j0+2, k+2, l+2))*sy(n, 2)
        ex(nn) = ex(nn) + a*sz(n, 2)

        ! Compute Ey on particle
        a = (sx(n, -1)*eyg(j-1, k0-1, l-1) + sx(n, 0)*eyg(j, k0-1, l-1) + sx(n,       &
        1)*eyg(j+1, k0-1, l-1) + sx(n, 2)*eyg(j+2, k0-1, l-1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l-1) + sx(n, 0)*eyg(j, k0, l-1) + sx(n,       &
        1)*eyg(j+1, k0, l-1) + sx(n, 2)*eyg(j+2, k0, l-1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l-1) + sx(n, 0)*eyg(j, k0+1, l-1) + sx(n,   &
        1)*eyg(j+1, k0+1, l-1) + sx(n, 2)*eyg(j+2, k0+1, l-1))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l-1) + sx(n, 0)*eyg(j, k0+2, l-1) + sx(n,   &
        1)*eyg(j+1, k0+2, l-1) + sx(n, 2)*eyg(j+2, k0+2, l-1))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, -1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l) + sx(n, 0)*eyg(j, k0-1, l) + sx(n,           &
        1)*eyg(j+1, k0-1, l) + sx(n, 2)*eyg(j+2, k0-1, l))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l) + sx(n, 0)*eyg(j, k0, l) + sx(n,           &
        1)*eyg(j+1, k0, l) + sx(n, 2)*eyg(j+2, k0, l))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l) + sx(n, 0)*eyg(j, k0+1, l) + sx(n,       &
        1)*eyg(j+1, k0+1, l) + sx(n, 2)*eyg(j+2, k0+1, l))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l) + sx(n, 0)*eyg(j, k0+2, l) + sx(n,       &
        1)*eyg(j+1, k0+2, l) + sx(n, 2)*eyg(j+2, k0+2, l))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, 0)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+1) + sx(n, 0)*eyg(j, k0-1, l+1) + sx(n,       &
        1)*eyg(j+1, k0-1, l+1) + sx(n, 2)*eyg(j+2, k0-1, l+1))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+1) + sx(n, 0)*eyg(j, k0, l+1) + sx(n,       &
        1)*eyg(j+1, k0, l+1) + sx(n, 2)*eyg(j+2, k0, l+1))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+1) + sx(n, 0)*eyg(j, k0+1, l+1) + sx(n,   &
        1)*eyg(j+1, k0+1, l+1) + sx(n, 2)*eyg(j+2, k0+1, l+1))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l+1) + sx(n, 0)*eyg(j, k0+2, l+1) + sx(n,   &
        1)*eyg(j+1, k0+2, l+1) + sx(n, 2)*eyg(j+2, k0+2, l+1))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, 1)
        a = (sx(n, -1)*eyg(j-1, k0-1, l+2) + sx(n, 0)*eyg(j, k0-1, l+2) + sx(n,       &
        1)*eyg(j+1, k0-1, l+2) + sx(n, 2)*eyg(j+2, k0-1, l+2))*sy0(n, -1)
        a = a + (sx(n, -1)*eyg(j-1, k0, l+2) + sx(n, 0)*eyg(j, k0, l+2) + sx(n,       &
        1)*eyg(j+1, k0, l+2) + sx(n, 2)*eyg(j+2, k0, l+2))*sy0(n, 0)
        a = a + (sx(n, -1)*eyg(j-1, k0+1, l+2) + sx(n, 0)*eyg(j, k0+1, l+2) + sx(n,   &
        1)*eyg(j+1, k0+1, l+2) + sx(n, 2)*eyg(j+2, k0+1, l+2))*sy0(n, 1)
        a = a + (sx(n, -1)*eyg(j-1, k0+2, l+2) + sx(n, 0)*eyg(j, k0+2, l+2) + sx(n,   &
        1)*eyg(j+1, k0+2, l+2) + sx(n, 2)*eyg(j+2, k0+2, l+2))*sy0(n, 2)
        ey(nn) = ey(nn) + a*sz(n, 2)

        ! Compute Ez on particle
        a = (sx(n, -1)*ezg(j-1, k-1, l0-1) + sx(n, 0)*ezg(j, k-1, l0-1) + sx(n,       &
        1)*ezg(j+1, k-1, l0-1) + sx(n, 2)*ezg(j+2, k-1, l0-1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0-1) + sx(n, 0)*ezg(j, k, l0-1) + sx(n,       &
        1)*ezg(j+1, k, l0-1) + sx(n, 2)*ezg(j+2, k, l0-1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0-1) + sx(n, 0)*ezg(j, k+1, l0-1) + sx(n,   &
        1)*ezg(j+1, k+1, l0-1) + sx(n, 2)*ezg(j+2, k+1, l0-1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0-1) + sx(n, 0)*ezg(j, k+2, l0-1) + sx(n,   &
        1)*ezg(j+1, k+2, l0-1) + sx(n, 2)*ezg(j+2, k+2, l0-1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*ezg(j-1, k-1, l0) + sx(n, 0)*ezg(j, k-1, l0) + sx(n,           &
        1)*ezg(j+1, k-1, l0) + sx(n, 2)*ezg(j+2, k-1, l0))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0) + sx(n, 0)*ezg(j, k, l0) + sx(n,           &
        1)*ezg(j+1, k, l0) + sx(n, 2)*ezg(j+2, k, l0))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0) + sx(n, 0)*ezg(j, k+1, l0) + sx(n,       &
        1)*ezg(j+1, k+1, l0) + sx(n, 2)*ezg(j+2, k+1, l0))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0) + sx(n, 0)*ezg(j, k+2, l0) + sx(n,       &
        1)*ezg(j+1, k+2, l0) + sx(n, 2)*ezg(j+2, k+2, l0))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*ezg(j-1, k-1, l0+1) + sx(n, 0)*ezg(j, k-1, l0+1) + sx(n,       &
        1)*ezg(j+1, k-1, l0+1) + sx(n, 2)*ezg(j+2, k-1, l0+1))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0+1) + sx(n, 0)*ezg(j, k, l0+1) + sx(n,       &
        1)*ezg(j+1, k, l0+1) + sx(n, 2)*ezg(j+2, k, l0+1))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0+1) + sx(n, 0)*ezg(j, k+1, l0+1) + sx(n,   &
        1)*ezg(j+1, k+1, l0+1) + sx(n, 2)*ezg(j+2, k+1, l0+1))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0+1) + sx(n, 0)*ezg(j, k+2, l0+1) + sx(n,   &
        1)*ezg(j+1, k+2, l0+1) + sx(n, 2)*ezg(j+2, k+2, l0+1))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 1)
        a = (sx(n, -1)*ezg(j-1, k-1, l0+2) + sx(n, 0)*ezg(j, k-1, l0+2) + sx(n,       &
        1)*ezg(j+1, k-1, l0+2) + sx(n, 2)*ezg(j+2, k-1, l0+2))*sy(n, -1)
        a = a + (sx(n, -1)*ezg(j-1, k, l0+2) + sx(n, 0)*ezg(j, k, l0+2) + sx(n,       &
        1)*ezg(j+1, k, l0+2) + sx(n, 2)*ezg(j+2, k, l0+2))*sy(n, 0)
        a = a + (sx(n, -1)*ezg(j-1, k+1, l0+2) + sx(n, 0)*ezg(j, k+1, l0+2) + sx(n,   &
        1)*ezg(j+1, k+1, l0+2) + sx(n, 2)*ezg(j+2, k+1, l0+2))*sy(n, 1)
        a = a + (sx(n, -1)*ezg(j-1, k+2, l0+2) + sx(n, 0)*ezg(j, k+2, l0+2) + sx(n,   &
        1)*ezg(j+1, k+2, l0+2) + sx(n, 2)*ezg(j+2, k+2, l0+2))*sy(n, 2)
        ez(nn) = ez(nn) + a*sz0(n, 2)

        ! Compute Bx on particle
        a = (sx(n, -1)*bxg(j-1, k0-1, l0-1) + sx(n, 0)*bxg(j, k0-1, l0-1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0-1) + sx(n, 2)*bxg(j+2, k0-1, l0-1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0-1) + sx(n, 0)*bxg(j, k0, l0-1) + sx(n,     &
        1)*bxg(j+1, k0, l0-1) + sx(n, 2)*bxg(j+2, k0, l0-1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0-1) + sx(n, 0)*bxg(j, k0+1, l0-1) + sx(n, &
        1)*bxg(j+1, k0+1, l0-1) + sx(n, 2)*bxg(j+2, k0+1, l0-1))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0-1) + sx(n, 0)*bxg(j, k0+2, l0-1) + sx(n, &
        1)*bxg(j+1, k0+2, l0-1) + sx(n, 2)*bxg(j+2, k0+2, l0-1))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, -1)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0) + sx(n, 0)*bxg(j, k0-1, l0) + sx(n,         &
        1)*bxg(j+1, k0-1, l0) + sx(n, 2)*bxg(j+2, k0-1, l0))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0) + sx(n, 0)*bxg(j, k0, l0) + sx(n,         &
        1)*bxg(j+1, k0, l0) + sx(n, 2)*bxg(j+2, k0, l0))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0) + sx(n, 0)*bxg(j, k0+1, l0) + sx(n,     &
        1)*bxg(j+1, k0+1, l0) + sx(n, 2)*bxg(j+2, k0+1, l0))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0) + sx(n, 0)*bxg(j, k0+2, l0) + sx(n,     &
        1)*bxg(j+1, k0+2, l0) + sx(n, 2)*bxg(j+2, k0+2, l0))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, 0)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0+1) + sx(n, 0)*bxg(j, k0-1, l0+1) + sx(n,     &
        1)*bxg(j+1, k0-1, l0+1) + sx(n, 2)*bxg(j+2, k0-1, l0+1))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0+1) + sx(n, 0)*bxg(j, k0, l0+1) + sx(n,     &
        1)*bxg(j+1, k0, l0+1) + sx(n, 2)*bxg(j+2, k0, l0+1))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0+1) + sx(n, 0)*bxg(j, k0+1, l0+1) + sx(n, &
        1)*bxg(j+1, k0+1, l0+1) + sx(n, 2)*bxg(j+2, k0+1, l0+1))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0+1) + sx(n, 0)*bxg(j, k0+2, l0+1) + sx(n, &
        1)*bxg(j+1, k0+2, l0+1) + sx(n, 2)*bxg(j+2, k0+2, l0+1))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, 1)
        a = (sx(n, -1)*bxg(j-1, k0-1, l0+2) + sx(n, 0)*bxg(j, k0-1, l0+2) + sx(n,     &
        1)*bxg(j+1, k0-1, l0+2) + sx(n, 2)*bxg(j+2, k0-1, l0+2))*sy0(n, -1)
        a = a + (sx(n, -1)*bxg(j-1, k0, l0+2) + sx(n, 0)*bxg(j, k0, l0+2) + sx(n,     &
        1)*bxg(j+1, k0, l0+2) + sx(n, 2)*bxg(j+2, k0, l0+2))*sy0(n, 0)
        a = a + (sx(n, -1)*bxg(j-1, k0+1, l0+2) + sx(n, 0)*bxg(j, k0+1, l0+2) + sx(n, &
        1)*bxg(j+1, k0+1, l0+2) + sx(n, 2)*bxg(j+2, k0+1, l0+2))*sy0(n, 1)
        a = a + (sx(n, -1)*bxg(j-1, k0+2, l0+2) + sx(n, 0)*bxg(j, k0+2, l0+2) + sx(n, &
        1)*bxg(j+1, k0+2, l0+2) + sx(n, 2)*bxg(j+2, k0+2, l0+2))*sy0(n, 2)
        bx(nn) = bx(nn) + a*sz0(n, 2)

        ! Compute By on particle
        a = (sx0(n, -1)*byg(j0-1, k-1, l0-1) + sx0(n, 0)*byg(j0, k-1, l0-1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0-1) + sx0(n, 2)*byg(j0+2, k-1, l0-1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0-1) + sx0(n, 0)*byg(j0, k, l0-1) + sx0(n,  &
        1)*byg(j0+1, k, l0-1) + sx0(n, 2)*byg(j0+2, k, l0-1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0-1) + sx0(n, 0)*byg(j0, k+1, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0-1) + sx0(n, 2)*byg(j0+2, k+1, l0-1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0-1) + sx0(n, 0)*byg(j0, k+2, l0-1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0-1) + sx0(n, 2)*byg(j0+2, k+2, l0-1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, -1)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0) + sx0(n, 0)*byg(j0, k-1, l0) + sx0(n,      &
        1)*byg(j0+1, k-1, l0) + sx0(n, 2)*byg(j0+2, k-1, l0))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0) + sx0(n, 0)*byg(j0, k, l0) + sx0(n,      &
        1)*byg(j0+1, k, l0) + sx0(n, 2)*byg(j0+2, k, l0))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0) + sx0(n, 0)*byg(j0, k+1, l0) + sx0(n,  &
        1)*byg(j0+1, k+1, l0) + sx0(n, 2)*byg(j0+2, k+1, l0))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0) + sx0(n, 0)*byg(j0, k+2, l0) + sx0(n,  &
        1)*byg(j0+1, k+2, l0) + sx0(n, 2)*byg(j0+2, k+2, l0))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 0)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0+1) + sx0(n, 0)*byg(j0, k-1, l0+1) + sx0(n,  &
        1)*byg(j0+1, k-1, l0+1) + sx0(n, 2)*byg(j0+2, k-1, l0+1))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0+1) + sx0(n, 0)*byg(j0, k, l0+1) + sx0(n,  &
        1)*byg(j0+1, k, l0+1) + sx0(n, 2)*byg(j0+2, k, l0+1))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0+1) + sx0(n, 0)*byg(j0, k+1, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0+1) + sx0(n, 2)*byg(j0+2, k+1, l0+1))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0+1) + sx0(n, 0)*byg(j0, k+2, l0+1) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0+1) + sx0(n, 2)*byg(j0+2, k+2, l0+1))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 1)
        a = (sx0(n, -1)*byg(j0-1, k-1, l0+2) + sx0(n, 0)*byg(j0, k-1, l0+2) + sx0(n,  &
        1)*byg(j0+1, k-1, l0+2) + sx0(n, 2)*byg(j0+2, k-1, l0+2))*sy(n, -1)
        a = a + (sx0(n, -1)*byg(j0-1, k, l0+2) + sx0(n, 0)*byg(j0, k, l0+2) + sx0(n,  &
        1)*byg(j0+1, k, l0+2) + sx0(n, 2)*byg(j0+2, k, l0+2))*sy(n, 0)
        a = a + (sx0(n, -1)*byg(j0-1, k+1, l0+2) + sx0(n, 0)*byg(j0, k+1, l0+2) +     &
        sx0(n, 1)*byg(j0+1, k+1, l0+2) + sx0(n, 2)*byg(j0+2, k+1, l0+2))*sy(n, 1)
        a = a + (sx0(n, -1)*byg(j0-1, k+2, l0+2) + sx0(n, 0)*byg(j0, k+2, l0+2) +     &
        sx0(n, 1)*byg(j0+1, k+2, l0+2) + sx0(n, 2)*byg(j0+2, k+2, l0+2))*sy(n, 2)
        by(nn) = by(nn) + a*sz0(n, 2)

        ! Compute Bz on particle
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l-1) + sx0(n, 0)*bzg(j0, k0-1, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l-1) + sx0(n, 2)*bzg(j0+2, k0-1, l-1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l-1) + sx0(n, 0)*bzg(j0, k0, l-1) + sx0(n,  &
        1)*bzg(j0+1, k0, l-1) + sx0(n, 2)*bzg(j0+2, k0, l-1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l-1) + sx0(n, 0)*bzg(j0, k0+1, l-1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l-1) + sx0(n, 2)*bzg(j0+2, k0+1, l-1))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l-1) + sx0(n, 0)*bzg(j0, k0+2, l-1) +     &
        sx0(n, 1)*bzg(j0+1, k0+2, l-1) + sx0(n, 2)*bzg(j0+2, k0+2, l-1))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, -1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l) + sx0(n, 0)*bzg(j0, k0-1, l) + sx0(n,      &
        1)*bzg(j0+1, k0-1, l) + sx0(n, 2)*bzg(j0+2, k0-1, l))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l) + sx0(n, 0)*bzg(j0, k0, l) + sx0(n,      &
        1)*bzg(j0+1, k0, l) + sx0(n, 2)*bzg(j0+2, k0, l))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l) + sx0(n, 0)*bzg(j0, k0+1, l) + sx0(n,  &
        1)*bzg(j0+1, k0+1, l) + sx0(n, 2)*bzg(j0+2, k0+1, l))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l) + sx0(n, 0)*bzg(j0, k0+2, l) + sx0(n,  &
        1)*bzg(j0+1, k0+2, l) + sx0(n, 2)*bzg(j0+2, k0+2, l))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, 0)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+1) + sx0(n, 0)*bzg(j0, k0-1, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+1) + sx0(n, 2)*bzg(j0+2, k0-1, l+1))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+1) + sx0(n, 0)*bzg(j0, k0, l+1) + sx0(n,  &
        1)*bzg(j0+1, k0, l+1) + sx0(n, 2)*bzg(j0+2, k0, l+1))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+1) + sx0(n, 0)*bzg(j0, k0+1, l+1) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+1) + sx0(n, 2)*bzg(j0+2, k0+1, l+1))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l+1) + sx0(n, 0)*bzg(j0, k0+2, l+1) +     &
        sx0(n, 1)*bzg(j0+1, k0+2, l+1) + sx0(n, 2)*bzg(j0+2, k0+2, l+1))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, 1)
        a = (sx0(n, -1)*bzg(j0-1, k0-1, l+2) + sx0(n, 0)*bzg(j0, k0-1, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0-1, l+2) + sx0(n, 2)*bzg(j0+2, k0-1, l+2))*sy0(n, -1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0, l+2) + sx0(n, 0)*bzg(j0, k0, l+2) + sx0(n,  &
        1)*bzg(j0+1, k0, l+2) + sx0(n, 2)*bzg(j0+2, k0, l+2))*sy0(n, 0)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+1, l+2) + sx0(n, 0)*bzg(j0, k0+1, l+2) +     &
        sx0(n, 1)*bzg(j0+1, k0+1, l+2) + sx0(n, 2)*bzg(j0+2, k0+1, l+2))*sy0(n, 1)
        a = a + (sx0(n, -1)*bzg(j0-1, k0+2, l+2) + sx0(n, 0)*bzg(j0, k0+2, l+2) +     &
        sx0(n, 1)*bzg(j0+1, k0+2, l+2) + sx0(n, 2)*bzg(j0+2, k0+2, l+2))*sy0(n, 2)
        bz(nn) = bz(nn) + a*sz(n, 2)

      ENDDO
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
!> Vectorized field gathering at order 3 with gathering of E and B
!> merged in a single loop.
!> This version use the pragma !$OMP SIMD private for sx, sy, sy, sx0, sy0, sz0.
!
!> @details
!> This function is vectorized and is the version 2 using a more optimized
!> gathering method more efficient for the vectorization.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 12/03/2016
!
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] ex, ey, ez particle electric field arrays
!> @param[inout] bx, by, bz particle magnetic field arrays
!> @param[in] xmin, ymin, zmin tile minimum grid position
!> @param[in] dx, dy, dz space steps in every directions
!> @param[in] dt time step
!> @param[in] exg, eyg, ezg electric field grids
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!> exg, eyg, ezg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 3 integers)
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the
!> bxg, byg, bzg arrays in each direction (1d arrays containing 3 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 3 integers)
!> @param[in] lvect vector size for cache blocking
!> @param[in] l_lower_order_in_v lower order for the interpolation
!
SUBROUTINE geteb3d_energy_conserving_vecV3_3_3_3(np, xp, yp, zp, ex, ey, ez, bx, by,  &
  bz, xmin, ymin, zmin, dx, dy, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,       &
  eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,            &
  byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, lvect,                           &
  l_lower_order_in_v, l_nodal)        !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp), intent(in)                :: np
  INTEGER(idp), intent(in)                :: exg_nguard(3), exg_nvalid(3),            &
  eyg_nguard(3), eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3), bxg_nguard(3),          &
  bxg_nvalid(3), byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
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
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip, j, k, l
  INTEGER(idp)                           :: j0, k0, l0
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq
  REAL(num)                              :: ozint, oxintsq, oyintsq, ozintsq
  REAL(num)                              :: a
  INTEGER(isp)                           :: nn
  REAL(num), DIMENSION(-1:2)             :: sx, sx0
  REAL(num), DIMENSION(-1:2)             :: sy, sy0
  REAL(num), DIMENSION(-1:2)             :: sz, sz0
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

  IF (l_lower_order_in_v ) THEN

    ! ___ Loop on partciles _______________________
    DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
      !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
      !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
      !!DIR PREFETCH ex:1:1
      !!DIR PREFETCH ey:1:1
      !!DIR PREFETCH ez:1:1
      !!DIR PREFETCH bx:1:1
      !!DIR PREFETCH by:1:1
      !!DIR PREFETCH bz:1:1
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
      !!DIR DISTRIBUTE POINT
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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx( 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy( 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz( 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        xintsq = xint*xint
        sx0( -1) = 0.5_num*(0.5_num-xint)**2
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
        a = (sx0(-1)*exg(j0-1, k-1, l-1) + sx0(0)*exg(j0, k-1, l-1) +                 &
        sx0(1)*exg(j0+1, k-1, l-1))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l-1) + sx0(0)*exg(j0, k, l-1) +                 &
        sx0(1)*exg(j0+1, k, l-1))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l-1) + sx0(0)*exg(j0, k+1, l-1) +             &
        sx0(1)*exg(j0+1, k+1, l-1))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l-1) + sx0(0)*exg(j0, k+2, l-1) +             &
        sx0(1)*exg(j0+1, k+2, l-1))*sy(2)
        ex(nn) = ex(nn) + a*sz(-1)
        a = (sx0(-1)*exg(j0-1, k-1, l) + sx0(0)*exg(j0, k-1, l) + sx0(1)*exg(j0+1,    &
        k-1, l))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l) + sx0(0)*exg(j0, k, l) + sx0(1)*exg(j0+1, k, &
        l))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l) + sx0(0)*exg(j0, k+1, l) +                 &
        sx0(1)*exg(j0+1, k+1, l))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l) + sx0(0)*exg(j0, k+2, l) +                 &
        sx0(1)*exg(j0+1, k+2, l))*sy(2)
        ex(nn) = ex(nn) + a*sz(0)
        a = (sx0(-1)*exg(j0-1, k-1, l+1) + sx0(0)*exg(j0, k-1, l+1) +                 &
        sx0(1)*exg(j0+1, k-1, l+1))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l+1) + sx0(0)*exg(j0, k, l+1) +                 &
        sx0(1)*exg(j0+1, k, l+1))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l+1) + sx0(0)*exg(j0, k+1, l+1) +             &
        sx0(1)*exg(j0+1, k+1, l+1))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l+1) + sx0(0)*exg(j0, k+2, l+1) +             &
        sx0(1)*exg(j0+1, k+2, l+1))*sy(2)
        ex(nn) = ex(nn) + a*sz(1)
        a = (sx0(-1)*exg(j0-1, k-1, l+2) + sx0(0)*exg(j0, k-1, l+2) +                 &
        sx0(1)*exg(j0+1, k-1, l+2))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l+2) + sx0(0)*exg(j0, k, l+2) +                 &
        sx0(1)*exg(j0+1, k, l+2))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l+2) + sx0(0)*exg(j0, k+1, l+2) +             &
        sx0(1)*exg(j0+1, k+1, l+2))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l+2) + sx0(0)*exg(j0, k+2, l+2) +             &
        sx0(1)*exg(j0+1, k+2, l+2))*sy(2)
        ex(nn) = ex(nn) + a*sz(2)

        ! Compute Ey on particle
        a = (sx(-1)*eyg(j-1, k0-1, l-1) + sx(0)*eyg(j, k0-1, l-1) + sx(1)*eyg(j+1,    &
        k0-1, l-1) + sx(2)*eyg(j+2, k0-1, l-1))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l-1) + sx(0)*eyg(j, k0, l-1) + sx(1)*eyg(j+1,    &
        k0, l-1) + sx(2)*eyg(j+2, k0, l-1))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l-1) + sx(0)*eyg(j, k0+1, l-1) +               &
        sx(1)*eyg(j+1, k0+1, l-1) + sx(2)*eyg(j+2, k0+1, l-1))*sy0(1)
        ey(nn) = ey(nn) + a*sz(-1)
        a = (sx(-1)*eyg(j-1, k0-1, l) + sx(0)*eyg(j, k0-1, l) + sx(1)*eyg(j+1, k0-1,  &
        l) + sx(2)*eyg(j+2, k0-1, l))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l) + sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0, l) &
        + sx(2)*eyg(j+2, k0, l))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l) + sx(0)*eyg(j, k0+1, l) + sx(1)*eyg(j+1,    &
        k0+1, l) + sx(2)*eyg(j+2, k0+1, l))*sy0(1)
        ey(nn) = ey(nn) + a*sz(0)
        a = (sx(-1)*eyg(j-1, k0-1, l+1) + sx(0)*eyg(j, k0-1, l+1) + sx(1)*eyg(j+1,    &
        k0-1, l+1) + sx(2)*eyg(j+2, k0-1, l+1))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l+1) + sx(0)*eyg(j, k0, l+1) + sx(1)*eyg(j+1,    &
        k0, l+1) + sx(2)*eyg(j+2, k0, l+1))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l+1) + sx(0)*eyg(j, k0+1, l+1) +               &
        sx(1)*eyg(j+1, k0+1, l+1) + sx(2)*eyg(j+2, k0+1, l+1))*sy0(1)
        ey(nn) = ey(nn) + a*sz(1)
        a = (sx(-1)*eyg(j-1, k0-1, l+2) + sx(0)*eyg(j, k0-1, l+2) + sx(1)*eyg(j+1,    &
        k0-1, l+2) + sx(2)*eyg(j+2, k0-1, l+2))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l+2) + sx(0)*eyg(j, k0, l+2) + sx(1)*eyg(j+1,    &
        k0, l+2) + sx(2)*eyg(j+2, k0, l+2))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l+2) + sx(0)*eyg(j, k0+1, l+2) +               &
        sx(1)*eyg(j+1, k0+1, l+2) + sx(2)*eyg(j+2, k0+1, l+2))*sy0(1)
        ey(nn) = ey(nn) + a*sz(2)

        ! Compute Ez on particle
        a = (sx(-1)*ezg(j-1, k-1, l0-1) + sx(0)*ezg(j, k-1, l0-1) + sx(1)*ezg(j+1,    &
        k-1, l0-1) + sx(2)*ezg(j+2, k-1, l0-1))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0-1) + sx(0)*ezg(j, k, l0-1) + sx(1)*ezg(j+1, k, &
        l0-1) + sx(2)*ezg(j+2, k, l0-1))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0-1) + sx(0)*ezg(j, k+1, l0-1) +               &
        sx(1)*ezg(j+1, k+1, l0-1) + sx(2)*ezg(j+2, k+1, l0-1))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0-1) + sx(0)*ezg(j, k+2, l0-1) +               &
        sx(1)*ezg(j+1, k+2, l0-1) + sx(2)*ezg(j+2, k+2, l0-1))*sy(2)
        ez(nn) = ez(nn) + a*sz0(-1)
        a = (sx(-1)*ezg(j-1, k-1, l0) + sx(0)*ezg(j, k-1, l0) + sx(1)*ezg(j+1, k-1,   &
        l0) + sx(2)*ezg(j+2, k-1, l0))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0) + sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k, l0) &
        + sx(2)*ezg(j+2, k, l0))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0) + sx(0)*ezg(j, k+1, l0) + sx(1)*ezg(j+1,    &
        k+1, l0) + sx(2)*ezg(j+2, k+1, l0))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0) + sx(0)*ezg(j, k+2, l0) + sx(1)*ezg(j+1,    &
        k+2, l0) + sx(2)*ezg(j+2, k+2, l0))*sy(2)
        ez(nn) = ez(nn) + a*sz0(0)
        a = (sx(-1)*ezg(j-1, k-1, l0+1) + sx(0)*ezg(j, k-1, l0+1) + sx(1)*ezg(j+1,    &
        k-1, l0+1) + sx(2)*ezg(j+2, k-1, l0+1))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0+1) + sx(0)*ezg(j, k, l0+1) + sx(1)*ezg(j+1, k, &
        l0+1) + sx(2)*ezg(j+2, k, l0+1))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0+1) + sx(0)*ezg(j, k+1, l0+1) +               &
        sx(1)*ezg(j+1, k+1, l0+1) + sx(2)*ezg(j+2, k+1, l0+1))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0+1) + sx(0)*ezg(j, k+2, l0+1) +               &
        sx(1)*ezg(j+1, k+2, l0+1) + sx(2)*ezg(j+2, k+2, l0+1))*sy(2)
        ez(nn) = ez(nn) + a*sz0(1)

        ! Compute Bx on particle
        a = (sx(-1)*bxg(j-1, k0-1, l0-1) + sx(0)*bxg(j, k0-1, l0-1) + sx(1)*bxg(j+1,  &
        k0-1, l0-1) + sx(2)*bxg(j+2, k0-1, l0-1))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0-1) + sx(0)*bxg(j, k0, l0-1) + sx(1)*bxg(j+1,  &
        k0, l0-1) + sx(2)*bxg(j+2, k0, l0-1))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0-1) + sx(0)*bxg(j, k0+1, l0-1) +             &
        sx(1)*bxg(j+1, k0+1, l0-1) + sx(2)*bxg(j+2, k0+1, l0-1))*sy0(1)
        bx(nn) = bx(nn) + a*sz0(-1)
        a = (sx(-1)*bxg(j-1, k0-1, l0) + sx(0)*bxg(j, k0-1, l0) + sx(1)*bxg(j+1,      &
        k0-1, l0) + sx(2)*bxg(j+2, k0-1, l0))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0) + sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0,  &
        l0) + sx(2)*bxg(j+2, k0, l0))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0) + sx(0)*bxg(j, k0+1, l0) + sx(1)*bxg(j+1,  &
        k0+1, l0) + sx(2)*bxg(j+2, k0+1, l0))*sy0(1)
        bx(nn) = bx(nn) + a*sz0(0)
        a = (sx(-1)*bxg(j-1, k0-1, l0+1) + sx(0)*bxg(j, k0-1, l0+1) + sx(1)*bxg(j+1,  &
        k0-1, l0+1) + sx(2)*bxg(j+2, k0-1, l0+1))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0+1) + sx(0)*bxg(j, k0, l0+1) + sx(1)*bxg(j+1,  &
        k0, l0+1) + sx(2)*bxg(j+2, k0, l0+1))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0+1) + sx(0)*bxg(j, k0+1, l0+1) +             &
        sx(1)*bxg(j+1, k0+1, l0+1) + sx(2)*bxg(j+2, k0+1, l0+1))*sy0(1)
        bx(nn) = bx(nn) + a*sz0(1)

        ! Compute By on particle
        a = (sx0(-1)*byg(j0-1, k-1, l0-1) + sx0(0)*byg(j0, k-1, l0-1) +               &
        sx0(1)*byg(j0+1, k-1, l0-1))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0-1) + sx0(0)*byg(j0, k, l0-1) +               &
        sx0(1)*byg(j0+1, k, l0-1))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0-1) + sx0(0)*byg(j0, k+1, l0-1) +           &
        sx0(1)*byg(j0+1, k+1, l0-1))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0-1) + sx0(0)*byg(j0, k+2, l0-1) +           &
        sx0(1)*byg(j0+1, k+2, l0-1))*sy(2)
        by(nn) = by(nn) + a*sz0(-1)
        a = (sx0(-1)*byg(j0-1, k-1, l0) + sx0(0)*byg(j0, k-1, l0) + sx0(1)*byg(j0+1,  &
        k-1, l0))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0) + sx0(0)*byg(j0, k, l0) + sx0(1)*byg(j0+1,  &
        k, l0))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0) + sx0(0)*byg(j0, k+1, l0) +               &
        sx0(1)*byg(j0+1, k+1, l0))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0) + sx0(0)*byg(j0, k+2, l0) +               &
        sx0(1)*byg(j0+1, k+2, l0))*sy(2)
        by(nn) = by(nn) + a*sz0(0)
        a = (sx0(-1)*byg(j0-1, k-1, l0+1) + sx0(0)*byg(j0, k-1, l0+1) +               &
        sx0(1)*byg(j0+1, k-1, l0+1))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0+1) + sx0(0)*byg(j0, k, l0+1) +               &
        sx0(1)*byg(j0+1, k, l0+1))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0+1) + sx0(0)*byg(j0, k+1, l0+1) +           &
        sx0(1)*byg(j0+1, k+1, l0+1))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0+1) + sx0(0)*byg(j0, k+2, l0+1) +           &
        sx0(1)*byg(j0+1, k+2, l0+1))*sy(2)
        by(nn) = by(nn) + a*sz0(1)

        ! Compute Bz on particle
        a = (sx0(-1)*bzg(j0-1, k0-1, l-1) + sx0(0)*bzg(j0, k0-1, l-1) +               &
        sx0(1)*bzg(j0+1, k0-1, l-1))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l-1) + sx0(0)*bzg(j0, k0, l-1) +               &
        sx0(1)*bzg(j0+1, k0, l-1))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l-1) + sx0(0)*bzg(j0, k0+1, l-1) +           &
        sx0(1)*bzg(j0+1, k0+1, l-1))*sy0(1)
        bz(nn) = bz(nn) + a*sz(-1)
        a = (sx0(-1)*bzg(j0-1, k0-1, l) + sx0(0)*bzg(j0, k0-1, l) + sx0(1)*bzg(j0+1,  &
        k0-1, l))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l) + sx0(0)*bzg(j0, k0, l) + sx0(1)*bzg(j0+1,  &
        k0, l))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l) + sx0(0)*bzg(j0, k0+1, l) +               &
        sx0(1)*bzg(j0+1, k0+1, l))*sy0(1)
        bz(nn) = bz(nn) + a*sz(0)
        a = (sx0(-1)*bzg(j0-1, k0-1, l+1) + sx0(0)*bzg(j0, k0-1, l+1) +               &
        sx0(1)*bzg(j0+1, k0-1, l+1))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l+1) + sx0(0)*bzg(j0, k0, l+1) +               &
        sx0(1)*bzg(j0+1, k0, l+1))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l+1) + sx0(0)*bzg(j0, k0+1, l+1) +           &
        sx0(1)*bzg(j0+1, k0+1, l+1))*sy0(1)
        bz(nn) = bz(nn) + a*sz(1)
        a = (sx0(-1)*bzg(j0-1, k0-1, l+2) + sx0(0)*bzg(j0, k0-1, l+2) +               &
        sx0(1)*bzg(j0+1, k0-1, l+2))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l+2) + sx0(0)*bzg(j0, k0, l+2) +               &
        sx0(1)*bzg(j0+1, k0, l+2))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l+2) + sx0(0)*bzg(j0, k0+1, l+2) +           &
        sx0(1)*bzg(j0+1, k0+1, l+2))*sy0(1)
        bz(nn) = bz(nn) + a*sz(2)

      ENDDO
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
      !!DIR DISTRIBUTE POINT
#endif
      ! Loop over the particles inside a block
      DO nn=ip, MIN(ip+lvect-1, np)

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
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx( 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy( 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz( -1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz( 2) = onesixth*zintsq*zint

        xint=x-stagger_shift-j0
        yint=y-stagger_shift-k0
        zint=z-stagger_shift-l0

        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0( -1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0( 2) = onesixth*xintsq*xint

        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0( -1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0( 2) = onesixth*yintsq*yint

        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0( 2) = onesixth*zintsq*zint

        ! Compute Ex on particle
        a = (sx0(-1)*exg(j0-1, k-1, l-1) + sx0(0)*exg(j0, k-1, l-1) +                 &
        sx0(1)*exg(j0+1, k-1, l-1) + sx0(2)*exg(j0+2, k-1, l-1))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l-1) + sx0(0)*exg(j0, k, l-1) +                 &
        sx0(1)*exg(j0+1, k, l-1) + sx0(2)*exg(j0+2, k, l-1))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l-1) + sx0(0)*exg(j0, k+1, l-1) +             &
        sx0(1)*exg(j0+1, k+1, l-1) + sx0(2)*exg(j0+2, k+1, l-1))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l-1) + sx0(0)*exg(j0, k+2, l-1) +             &
        sx0(1)*exg(j0+1, k+2, l-1) + sx0(2)*exg(j0+2, k+2, l-1))*sy(2)
        ex(nn) = ex(nn) + a*sz(-1)
        a = (sx0(-1)*exg(j0-1, k-1, l) + sx0(0)*exg(j0, k-1, l) + sx0(1)*exg(j0+1,    &
        k-1, l) + sx0(2)*exg(j0+2, k-1, l))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l) + sx0(0)*exg(j0, k, l) + sx0(1)*exg(j0+1, k, &
        l) + sx0(2)*exg(j0+2, k, l))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l) + sx0(0)*exg(j0, k+1, l) +                 &
        sx0(1)*exg(j0+1, k+1, l) + sx0(2)*exg(j0+2, k+1, l))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l) + sx0(0)*exg(j0, k+2, l) +                 &
        sx0(1)*exg(j0+1, k+2, l) + sx0(2)*exg(j0+2, k+2, l))*sy(2)
        ex(nn) = ex(nn) + a*sz(0)
        a = (sx0(-1)*exg(j0-1, k-1, l+1) + sx0(0)*exg(j0, k-1, l+1) +                 &
        sx0(1)*exg(j0+1, k-1, l+1) + sx0(2)*exg(j0+2, k-1, l+1))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l+1) + sx0(0)*exg(j0, k, l+1) +                 &
        sx0(1)*exg(j0+1, k, l+1) + sx0(2)*exg(j0+2, k, l+1))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l+1) + sx0(0)*exg(j0, k+1, l+1) +             &
        sx0(1)*exg(j0+1, k+1, l+1) + sx0(2)*exg(j0+2, k+1, l+1))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l+1) + sx0(0)*exg(j0, k+2, l+1) +             &
        sx0(1)*exg(j0+1, k+2, l+1) + sx0(2)*exg(j0+2, k+2, l+1))*sy(2)
        ex(nn) = ex(nn) + a*sz(1)
        a = (sx0(-1)*exg(j0-1, k-1, l+2) + sx0(0)*exg(j0, k-1, l+2) +                 &
        sx0(1)*exg(j0+1, k-1, l+2) + sx0(2)*exg(j0+2, k-1, l+2))*sy(-1)
        a = a + (sx0(-1)*exg(j0-1, k, l+2) + sx0(0)*exg(j0, k, l+2) +                 &
        sx0(1)*exg(j0+1, k, l+2) + sx0(2)*exg(j0+2, k, l+2))*sy(0)
        a = a + (sx0(-1)*exg(j0-1, k+1, l+2) + sx0(0)*exg(j0, k+1, l+2) +             &
        sx0(1)*exg(j0+1, k+1, l+2) + sx0(2)*exg(j0+2, k+1, l+2))*sy(1)
        a = a + (sx0(-1)*exg(j0-1, k+2, l+2) + sx0(0)*exg(j0, k+2, l+2) +             &
        sx0(1)*exg(j0+1, k+2, l+2) + sx0(2)*exg(j0+2, k+2, l+2))*sy(2)
        ex(nn) = ex(nn) + a*sz(2)

        ! Compute Ey on particle
        a = (sx(-1)*eyg(j-1, k0-1, l-1) + sx(0)*eyg(j, k0-1, l-1) + sx(1)*eyg(j+1,    &
        k0-1, l-1) + sx(2)*eyg(j+2, k0-1, l-1))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l-1) + sx(0)*eyg(j, k0, l-1) + sx(1)*eyg(j+1,    &
        k0, l-1) + sx(2)*eyg(j+2, k0, l-1))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l-1) + sx(0)*eyg(j, k0+1, l-1) +               &
        sx(1)*eyg(j+1, k0+1, l-1) + sx(2)*eyg(j+2, k0+1, l-1))*sy0(1)
        a = a + (sx(-1)*eyg(j-1, k0+2, l-1) + sx(0)*eyg(j, k0+2, l-1) +               &
        sx(1)*eyg(j+1, k0+2, l-1) + sx(2)*eyg(j+2, k0+2, l-1))*sy0(2)
        ey(nn) = ey(nn) + a*sz(-1)
        a = (sx(-1)*eyg(j-1, k0-1, l) + sx(0)*eyg(j, k0-1, l) + sx(1)*eyg(j+1, k0-1,  &
        l) + sx(2)*eyg(j+2, k0-1, l))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l) + sx(0)*eyg(j, k0, l) + sx(1)*eyg(j+1, k0, l) &
        + sx(2)*eyg(j+2, k0, l))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l) + sx(0)*eyg(j, k0+1, l) + sx(1)*eyg(j+1,    &
        k0+1, l) + sx(2)*eyg(j+2, k0+1, l))*sy0(1)
        a = a + (sx(-1)*eyg(j-1, k0+2, l) + sx(0)*eyg(j, k0+2, l) + sx(1)*eyg(j+1,    &
        k0+2, l) + sx(2)*eyg(j+2, k0+2, l))*sy0(2)
        ey(nn) = ey(nn) + a*sz(0)
        a = (sx(-1)*eyg(j-1, k0-1, l+1) + sx(0)*eyg(j, k0-1, l+1) + sx(1)*eyg(j+1,    &
        k0-1, l+1) + sx(2)*eyg(j+2, k0-1, l+1))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l+1) + sx(0)*eyg(j, k0, l+1) + sx(1)*eyg(j+1,    &
        k0, l+1) + sx(2)*eyg(j+2, k0, l+1))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l+1) + sx(0)*eyg(j, k0+1, l+1) +               &
        sx(1)*eyg(j+1, k0+1, l+1) + sx(2)*eyg(j+2, k0+1, l+1))*sy0(1)
        a = a + (sx(-1)*eyg(j-1, k0+2, l+1) + sx(0)*eyg(j, k0+2, l+1) +               &
        sx(1)*eyg(j+1, k0+2, l+1) + sx(2)*eyg(j+2, k0+2, l+1))*sy0(2)
        ey(nn) = ey(nn) + a*sz(1)
        a = (sx(-1)*eyg(j-1, k0-1, l+2) + sx(0)*eyg(j, k0-1, l+2) + sx(1)*eyg(j+1,    &
        k0-1, l+2) + sx(2)*eyg(j+2, k0-1, l+2))*sy0(-1)
        a = a + (sx(-1)*eyg(j-1, k0, l+2) + sx(0)*eyg(j, k0, l+2) + sx(1)*eyg(j+1,    &
        k0, l+2) + sx(2)*eyg(j+2, k0, l+2))*sy0(0)
        a = a + (sx(-1)*eyg(j-1, k0+1, l+2) + sx(0)*eyg(j, k0+1, l+2) +               &
        sx(1)*eyg(j+1, k0+1, l+2) + sx(2)*eyg(j+2, k0+1, l+2))*sy0(1)
        a = a + (sx(-1)*eyg(j-1, k0+2, l+2) + sx(0)*eyg(j, k0+2, l+2) +               &
        sx(1)*eyg(j+1, k0+2, l+2) + sx(2)*eyg(j+2, k0+2, l+2))*sy0(2)
        ey(nn) = ey(nn) + a*sz(2)

        ! Compute Ez on particle
        a = (sx(-1)*ezg(j-1, k-1, l0-1) + sx(0)*ezg(j, k-1, l0-1) + sx(1)*ezg(j+1,    &
        k-1, l0-1) + sx(2)*ezg(j+2, k-1, l0-1))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0-1) + sx(0)*ezg(j, k, l0-1) + sx(1)*ezg(j+1, k, &
        l0-1) + sx(2)*ezg(j+2, k, l0-1))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0-1) + sx(0)*ezg(j, k+1, l0-1) +               &
        sx(1)*ezg(j+1, k+1, l0-1) + sx(2)*ezg(j+2, k+1, l0-1))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0-1) + sx(0)*ezg(j, k+2, l0-1) +               &
        sx(1)*ezg(j+1, k+2, l0-1) + sx(2)*ezg(j+2, k+2, l0-1))*sy(2)
        ez(nn) = ez(nn) + a*sz0(-1)
        a = (sx(-1)*ezg(j-1, k-1, l0) + sx(0)*ezg(j, k-1, l0) + sx(1)*ezg(j+1, k-1,   &
        l0) + sx(2)*ezg(j+2, k-1, l0))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0) + sx(0)*ezg(j, k, l0) + sx(1)*ezg(j+1, k, l0) &
        + sx(2)*ezg(j+2, k, l0))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0) + sx(0)*ezg(j, k+1, l0) + sx(1)*ezg(j+1,    &
        k+1, l0) + sx(2)*ezg(j+2, k+1, l0))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0) + sx(0)*ezg(j, k+2, l0) + sx(1)*ezg(j+1,    &
        k+2, l0) + sx(2)*ezg(j+2, k+2, l0))*sy(2)
        ez(nn) = ez(nn) + a*sz0(0)
        a = (sx(-1)*ezg(j-1, k-1, l0+1) + sx(0)*ezg(j, k-1, l0+1) + sx(1)*ezg(j+1,    &
        k-1, l0+1) + sx(2)*ezg(j+2, k-1, l0+1))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0+1) + sx(0)*ezg(j, k, l0+1) + sx(1)*ezg(j+1, k, &
        l0+1) + sx(2)*ezg(j+2, k, l0+1))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0+1) + sx(0)*ezg(j, k+1, l0+1) +               &
        sx(1)*ezg(j+1, k+1, l0+1) + sx(2)*ezg(j+2, k+1, l0+1))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0+1) + sx(0)*ezg(j, k+2, l0+1) +               &
        sx(1)*ezg(j+1, k+2, l0+1) + sx(2)*ezg(j+2, k+2, l0+1))*sy(2)
        ez(nn) = ez(nn) + a*sz0(1)
        a = (sx(-1)*ezg(j-1, k-1, l0+2) + sx(0)*ezg(j, k-1, l0+2) + sx(1)*ezg(j+1,    &
        k-1, l0+2) + sx(2)*ezg(j+2, k-1, l0+2))*sy(-1)
        a = a + (sx(-1)*ezg(j-1, k, l0+2) + sx(0)*ezg(j, k, l0+2) + sx(1)*ezg(j+1, k, &
        l0+2) + sx(2)*ezg(j+2, k, l0+2))*sy(0)
        a = a + (sx(-1)*ezg(j-1, k+1, l0+2) + sx(0)*ezg(j, k+1, l0+2) +               &
        sx(1)*ezg(j+1, k+1, l0+2) + sx(2)*ezg(j+2, k+1, l0+2))*sy(1)
        a = a + (sx(-1)*ezg(j-1, k+2, l0+2) + sx(0)*ezg(j, k+2, l0+2) +               &
        sx(1)*ezg(j+1, k+2, l0+2) + sx(2)*ezg(j+2, k+2, l0+2))*sy(2)
        ez(nn) = ez(nn) + a*sz0(2)

        ! Compute Bx on particle
        a = (sx(-1)*bxg(j-1, k0-1, l0-1) + sx(0)*bxg(j, k0-1, l0-1) + sx(1)*bxg(j+1,  &
        k0-1, l0-1) + sx(2)*bxg(j+2, k0-1, l0-1))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0-1) + sx(0)*bxg(j, k0, l0-1) + sx(1)*bxg(j+1,  &
        k0, l0-1) + sx(2)*bxg(j+2, k0, l0-1))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0-1) + sx(0)*bxg(j, k0+1, l0-1) +             &
        sx(1)*bxg(j+1, k0+1, l0-1) + sx(2)*bxg(j+2, k0+1, l0-1))*sy0(1)
        a = a + (sx(-1)*bxg(j-1, k0+2, l0-1) + sx(0)*bxg(j, k0+2, l0-1) +             &
        sx(1)*bxg(j+1, k0+2, l0-1) + sx(2)*bxg(j+2, k0+2, l0-1))*sy0(2)
        bx(nn) = bx(nn) + a*sz0(-1)
        a = (sx(-1)*bxg(j-1, k0-1, l0) + sx(0)*bxg(j, k0-1, l0) + sx(1)*bxg(j+1,      &
        k0-1, l0) + sx(2)*bxg(j+2, k0-1, l0))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0) + sx(0)*bxg(j, k0, l0) + sx(1)*bxg(j+1, k0,  &
        l0) + sx(2)*bxg(j+2, k0, l0))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0) + sx(0)*bxg(j, k0+1, l0) + sx(1)*bxg(j+1,  &
        k0+1, l0) + sx(2)*bxg(j+2, k0+1, l0))*sy0(1)
        a = a + (sx(-1)*bxg(j-1, k0+2, l0) + sx(0)*bxg(j, k0+2, l0) + sx(1)*bxg(j+1,  &
        k0+2, l0) + sx(2)*bxg(j+2, k0+2, l0))*sy0(2)
        bx(nn) = bx(nn) + a*sz0(0)
        a = (sx(-1)*bxg(j-1, k0-1, l0+1) + sx(0)*bxg(j, k0-1, l0+1) + sx(1)*bxg(j+1,  &
        k0-1, l0+1) + sx(2)*bxg(j+2, k0-1, l0+1))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0+1) + sx(0)*bxg(j, k0, l0+1) + sx(1)*bxg(j+1,  &
        k0, l0+1) + sx(2)*bxg(j+2, k0, l0+1))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0+1) + sx(0)*bxg(j, k0+1, l0+1) +             &
        sx(1)*bxg(j+1, k0+1, l0+1) + sx(2)*bxg(j+2, k0+1, l0+1))*sy0(1)
        a = a + (sx(-1)*bxg(j-1, k0+2, l0+1) + sx(0)*bxg(j, k0+2, l0+1) +             &
        sx(1)*bxg(j+1, k0+2, l0+1) + sx(2)*bxg(j+2, k0+2, l0+1))*sy0(2)
        bx(nn) = bx(nn) + a*sz0(1)
        a = (sx(-1)*bxg(j-1, k0-1, l0+2) + sx(0)*bxg(j, k0-1, l0+2) + sx(1)*bxg(j+1,  &
        k0-1, l0+2) + sx(2)*bxg(j+2, k0-1, l0+2))*sy0(-1)
        a = a + (sx(-1)*bxg(j-1, k0, l0+2) + sx(0)*bxg(j, k0, l0+2) + sx(1)*bxg(j+1,  &
        k0, l0+2) + sx(2)*bxg(j+2, k0, l0+2))*sy0(0)
        a = a + (sx(-1)*bxg(j-1, k0+1, l0+2) + sx(0)*bxg(j, k0+1, l0+2) +             &
        sx(1)*bxg(j+1, k0+1, l0+2) + sx(2)*bxg(j+2, k0+1, l0+2))*sy0(1)
        a = a + (sx(-1)*bxg(j-1, k0+2, l0+2) + sx(0)*bxg(j, k0+2, l0+2) +             &
        sx(1)*bxg(j+1, k0+2, l0+2) + sx(2)*bxg(j+2, k0+2, l0+2))*sy0(2)
        bx(nn) = bx(nn) + a*sz0(2)

        ! Compute By on particle
        a = (sx0(-1)*byg(j0-1, k-1, l0-1) + sx0(0)*byg(j0, k-1, l0-1) +               &
        sx0(1)*byg(j0+1, k-1, l0-1) + sx0(2)*byg(j0+2, k-1, l0-1))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0-1) + sx0(0)*byg(j0, k, l0-1) +               &
        sx0(1)*byg(j0+1, k, l0-1) + sx0(2)*byg(j0+2, k, l0-1))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0-1) + sx0(0)*byg(j0, k+1, l0-1) +           &
        sx0(1)*byg(j0+1, k+1, l0-1) + sx0(2)*byg(j0+2, k+1, l0-1))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0-1) + sx0(0)*byg(j0, k+2, l0-1) +           &
        sx0(1)*byg(j0+1, k+2, l0-1) + sx0(2)*byg(j0+2, k+2, l0-1))*sy(2)
        by(nn) = by(nn) + a*sz0(-1)
        a = (sx0(-1)*byg(j0-1, k-1, l0) + sx0(0)*byg(j0, k-1, l0) + sx0(1)*byg(j0+1,  &
        k-1, l0) + sx0(2)*byg(j0+2, k-1, l0))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0) + sx0(0)*byg(j0, k, l0) + sx0(1)*byg(j0+1,  &
        k, l0) + sx0(2)*byg(j0+2, k, l0))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0) + sx0(0)*byg(j0, k+1, l0) +               &
        sx0(1)*byg(j0+1, k+1, l0) + sx0(2)*byg(j0+2, k+1, l0))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0) + sx0(0)*byg(j0, k+2, l0) +               &
        sx0(1)*byg(j0+1, k+2, l0) + sx0(2)*byg(j0+2, k+2, l0))*sy(2)
        by(nn) = by(nn) + a*sz0(0)
        a = (sx0(-1)*byg(j0-1, k-1, l0+1) + sx0(0)*byg(j0, k-1, l0+1) +               &
        sx0(1)*byg(j0+1, k-1, l0+1) + sx0(2)*byg(j0+2, k-1, l0+1))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0+1) + sx0(0)*byg(j0, k, l0+1) +               &
        sx0(1)*byg(j0+1, k, l0+1) + sx0(2)*byg(j0+2, k, l0+1))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0+1) + sx0(0)*byg(j0, k+1, l0+1) +           &
        sx0(1)*byg(j0+1, k+1, l0+1) + sx0(2)*byg(j0+2, k+1, l0+1))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0+1) + sx0(0)*byg(j0, k+2, l0+1) +           &
        sx0(1)*byg(j0+1, k+2, l0+1) + sx0(2)*byg(j0+2, k+2, l0+1))*sy(2)
        by(nn) = by(nn) + a*sz0(1)
        a = (sx0(-1)*byg(j0-1, k-1, l0+2) + sx0(0)*byg(j0, k-1, l0+2) +               &
        sx0(1)*byg(j0+1, k-1, l0+2) + sx0(2)*byg(j0+2, k-1, l0+2))*sy(-1)
        a = a + (sx0(-1)*byg(j0-1, k, l0+2) + sx0(0)*byg(j0, k, l0+2) +               &
        sx0(1)*byg(j0+1, k, l0+2) + sx0(2)*byg(j0+2, k, l0+2))*sy(0)
        a = a + (sx0(-1)*byg(j0-1, k+1, l0+2) + sx0(0)*byg(j0, k+1, l0+2) +           &
        sx0(1)*byg(j0+1, k+1, l0+2) + sx0(2)*byg(j0+2, k+1, l0+2))*sy(1)
        a = a + (sx0(-1)*byg(j0-1, k+2, l0+2) + sx0(0)*byg(j0, k+2, l0+2) +           &
        sx0(1)*byg(j0+1, k+2, l0+2) + sx0(2)*byg(j0+2, k+2, l0+2))*sy(2)
        by(nn) = by(nn) + a*sz0(2)

        ! Compute Bz on particle
        a = (sx0(-1)*bzg(j0-1, k0-1, l-1) + sx0(0)*bzg(j0, k0-1, l-1) +               &
        sx0(1)*bzg(j0+1, k0-1, l-1) + sx0(2)*bzg(j0+2, k0-1, l-1))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l-1) + sx0(0)*bzg(j0, k0, l-1) +               &
        sx0(1)*bzg(j0+1, k0, l-1) + sx0(2)*bzg(j0+2, k0, l-1))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l-1) + sx0(0)*bzg(j0, k0+1, l-1) +           &
        sx0(1)*bzg(j0+1, k0+1, l-1) + sx0(2)*bzg(j0+2, k0+1, l-1))*sy0(1)
        a = a + (sx0(-1)*bzg(j0-1, k0+2, l-1) + sx0(0)*bzg(j0, k0+2, l-1) +           &
        sx0(1)*bzg(j0+1, k0+2, l-1) + sx0(2)*bzg(j0+2, k0+2, l-1))*sy0(2)
        bz(nn) = bz(nn) + a*sz(-1)
        a = (sx0(-1)*bzg(j0-1, k0-1, l) + sx0(0)*bzg(j0, k0-1, l) + sx0(1)*bzg(j0+1,  &
        k0-1, l) + sx0(2)*bzg(j0+2, k0-1, l))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l) + sx0(0)*bzg(j0, k0, l) + sx0(1)*bzg(j0+1,  &
        k0, l) + sx0(2)*bzg(j0+2, k0, l))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l) + sx0(0)*bzg(j0, k0+1, l) +               &
        sx0(1)*bzg(j0+1, k0+1, l) + sx0(2)*bzg(j0+2, k0+1, l))*sy0(1)
        a = a + (sx0(-1)*bzg(j0-1, k0+2, l) + sx0(0)*bzg(j0, k0+2, l) +               &
        sx0(1)*bzg(j0+1, k0+2, l) + sx0(2)*bzg(j0+2, k0+2, l))*sy0(2)
        bz(nn) = bz(nn) + a*sz(0)
        a = (sx0(-1)*bzg(j0-1, k0-1, l+1) + sx0(0)*bzg(j0, k0-1, l+1) +               &
        sx0(1)*bzg(j0+1, k0-1, l+1) + sx0(2)*bzg(j0+2, k0-1, l+1))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l+1) + sx0(0)*bzg(j0, k0, l+1) +               &
        sx0(1)*bzg(j0+1, k0, l+1) + sx0(2)*bzg(j0+2, k0, l+1))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l+1) + sx0(0)*bzg(j0, k0+1, l+1) +           &
        sx0(1)*bzg(j0+1, k0+1, l+1) + sx0(2)*bzg(j0+2, k0+1, l+1))*sy0(1)
        a = a + (sx0(-1)*bzg(j0-1, k0+2, l+1) + sx0(0)*bzg(j0, k0+2, l+1) +           &
        sx0(1)*bzg(j0+1, k0+2, l+1) + sx0(2)*bzg(j0+2, k0+2, l+1))*sy0(2)
        bz(nn) = bz(nn) + a*sz(1)
        a = (sx0(-1)*bzg(j0-1, k0-1, l+2) + sx0(0)*bzg(j0, k0-1, l+2) +               &
        sx0(1)*bzg(j0+1, k0-1, l+2) + sx0(2)*bzg(j0+2, k0-1, l+2))*sy0(-1)
        a = a + (sx0(-1)*bzg(j0-1, k0, l+2) + sx0(0)*bzg(j0, k0, l+2) +               &
        sx0(1)*bzg(j0+1, k0, l+2) + sx0(2)*bzg(j0+2, k0, l+2))*sy0(0)
        a = a + (sx0(-1)*bzg(j0-1, k0+1, l+2) + sx0(0)*bzg(j0, k0+1, l+2) +           &
        sx0(1)*bzg(j0+1, k0+1, l+2) + sx0(2)*bzg(j0+2, k0+1, l+2))*sy0(1)
        a = a + (sx0(-1)*bzg(j0-1, k0+2, l+2) + sx0(0)*bzg(j0, k0+2, l+2) +           &
        sx0(1)*bzg(j0+1, k0+2, l+2) + sx0(2)*bzg(j0+2, k0+2, l+2))*sy0(2)
        bz(nn) = bz(nn) + a*sz(2)
      ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Field gathering by block at order 3 with gathering of E and B
!> splited in several small subloops.
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
SUBROUTINE geteb3d_energy_conserving_blockvec_3_3_3(np, xp, yp, zp, ex, ey, ez, bx,   &
  by, bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,     &
  eyg, ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip
  INTEGER(isp), DIMENSION(lvect)         :: j, k, l
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j, k, l
#endif
  INTEGER(isp), DIMENSION(lvect)         :: j0, k0, l0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j0, k0, l0
#endif
  INTEGER(idp)                           :: jj, kk, ll
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq,     &
  ozint, oxintsq, oyintsq, ozintsq
  INTEGER(isp)                           :: nn, n
  REAL(num), DIMENSION(lvect, -1:2)       :: sx
  REAL(num), DIMENSION(lvect, -1:2)       :: sy
  REAL(num), DIMENSION(lvect, -1:2)       :: sz
  REAL(num), DIMENSION(lvect, -1:1)       :: sx0
  REAL(num), DIMENSION(lvect, -1:1)       :: sy0
  REAL(num), DIMENSION(lvect, -1:1)       :: sz0
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

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  ! ___ Loop on partciles _______________________
  DO ip=1, np, lvect

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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(n, -1) = onesixth*oxintsq*oxint
      sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx(n, 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(n, -1) = onesixth*oyintsq*oyint
      sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy(n, 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(n, -1) = onesixth*ozintsq*ozint
      sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz(n, 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0(n)
      yint=y-stagger_shift-k0(n)
      zint=z-stagger_shift-l0(n)

      xintsq = xint*xint
      sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
      sx0(n, 0) = 0.75_num-xintsq
      sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

      yintsq = yint*yint
      sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
      sy0(n, 0) = 0.75_num-yintsq
      sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

      zintsq = zint*zint
      sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
      sz0(n, 0) = 0.75_num-zintsq
      sz0(n, 1) = 0.5_num*(0.5_num+zint)**2
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ex:64
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

      ! Compute Ex on particle
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, -1)*exg(j0(n)-1, k(n)-1, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, -1)*exg(j0(n), k(n)-1, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, -1)*exg(j0(n)+1, k(n)-1, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, -1)*exg(j0(n)-1, k(n), l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, -1)*exg(j0(n), k(n), l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, -1)*exg(j0(n)+1, k(n), l(n)-1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, -1)*exg(j0(n)-1, k(n)+1, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, -1)*exg(j0(n), k(n)+1, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, -1)*exg(j0(n)+1, k(n)+1, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, -1)*exg(j0(n)-1, k(n)+2, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, -1)*exg(j0(n), k(n)+2, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, -1)*exg(j0(n)+1, k(n)+2, l(n)-1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 0)*exg(j0(n)-1, k(n)-1, l(n))
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 0)*exg(j0(n), k(n)-1, l(n))
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 0)*exg(j0(n)+1, k(n)-1, l(n))
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 0)*exg(j0(n)-1, k(n), l(n))
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 0)*exg(j0(n), k(n), l(n))
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 0)*exg(j0(n)+1, k(n), l(n))
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 0)*exg(j0(n)-1, k(n)+1, l(n))
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 0)*exg(j0(n), k(n)+1, l(n))
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 0)*exg(j0(n)+1, k(n)+1, l(n))
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 0)*exg(j0(n)-1, k(n)+2, l(n))
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 0)*exg(j0(n), k(n)+2, l(n))
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 0)*exg(j0(n)+1, k(n)+2, l(n))
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 1)*exg(j0(n)-1, k(n)-1, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 1)*exg(j0(n), k(n)-1, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 1)*exg(j0(n)+1, k(n)-1, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 1)*exg(j0(n)-1, k(n), l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 1)*exg(j0(n), k(n), l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 1)*exg(j0(n)+1, k(n), l(n)+1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 1)*exg(j0(n)-1, k(n)+1, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 1)*exg(j0(n), k(n)+1, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 1)*exg(j0(n)+1, k(n)+1, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 1)*exg(j0(n)-1, k(n)+2, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 1)*exg(j0(n), k(n)+2, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 1)*exg(j0(n)+1, k(n)+2, l(n)+1)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, -1)*sz(n, 2)*exg(j0(n)-1, k(n)-1, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, -1)*sz(n, 2)*exg(j0(n), k(n)-1, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, -1)*sz(n, 2)*exg(j0(n)+1, k(n)-1, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 0)*sz(n, 2)*exg(j0(n)-1, k(n), l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 0)*sz(n, 2)*exg(j0(n), k(n), l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 0)*sz(n, 2)*exg(j0(n)+1, k(n), l(n)+2)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 1)*sz(n, 2)*exg(j0(n)-1, k(n)+1, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 1)*sz(n, 2)*exg(j0(n), k(n)+1, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 1)*sz(n, 2)*exg(j0(n)+1, k(n)+1, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, -1)*sy(n, 2)*sz(n, 2)*exg(j0(n)-1, k(n)+2, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 0)*sy(n, 2)*sz(n, 2)*exg(j0(n), k(n)+2, l(n)+2)
      ex(nn) = ex(nn) + sx0(n, 1)*sy(n, 2)*sz(n, 2)*exg(j0(n)+1, k(n)+2, l(n)+2)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ey:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1
      ! Compute Ey on particle
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, -1)*eyg(j(n)-1, k0(n)-1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, -1)*eyg(j(n), k0(n)-1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, -1)*eyg(j(n)+1, k0(n)-1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, -1)*eyg(j(n)+2, k0(n)-1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, -1)*eyg(j(n)-1, k0(n), l(n)-1)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, -1)*eyg(j(n), k0(n), l(n)-1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, -1)*eyg(j(n)+1, k0(n), l(n)-1)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, -1)*eyg(j(n)+2, k0(n), l(n)-1)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, -1)*eyg(j(n)-1, k0(n)+1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, -1)*eyg(j(n), k0(n)+1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, -1)*eyg(j(n)+1, k0(n)+1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, -1)*eyg(j(n)+2, k0(n)+1, l(n)-1)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 0)*eyg(j(n)-1, k0(n)-1, l(n))
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 0)*eyg(j(n), k0(n)-1, l(n))
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 0)*eyg(j(n)+1, k0(n)-1, l(n))
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 0)*eyg(j(n)+2, k0(n)-1, l(n))
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 0)*eyg(j(n)-1, k0(n), l(n))
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 0)*eyg(j(n), k0(n), l(n))
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 0)*eyg(j(n)+1, k0(n), l(n))
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 0)*eyg(j(n)+2, k0(n), l(n))
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 0)*eyg(j(n)-1, k0(n)+1, l(n))
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 0)*eyg(j(n), k0(n)+1, l(n))
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 0)*eyg(j(n)+1, k0(n)+1, l(n))
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 0)*eyg(j(n)+2, k0(n)+1, l(n))
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 1)*eyg(j(n)-1, k0(n)-1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 1)*eyg(j(n), k0(n)-1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 1)*eyg(j(n)+1, k0(n)-1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 1)*eyg(j(n)+2, k0(n)-1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 1)*eyg(j(n)-1, k0(n), l(n)+1)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 1)*eyg(j(n), k0(n), l(n)+1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 1)*eyg(j(n)+1, k0(n), l(n)+1)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 1)*eyg(j(n)+2, k0(n), l(n)+1)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 1)*eyg(j(n)-1, k0(n)+1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 1)*eyg(j(n), k0(n)+1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 1)*eyg(j(n)+1, k0(n)+1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 1)*eyg(j(n)+2, k0(n)+1, l(n)+1)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, -1)*sz(n, 2)*eyg(j(n)-1, k0(n)-1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, -1)*sz(n, 2)*eyg(j(n), k0(n)-1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, -1)*sz(n, 2)*eyg(j(n)+1, k0(n)-1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, -1)*sz(n, 2)*eyg(j(n)+2, k0(n)-1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 0)*sz(n, 2)*eyg(j(n)-1, k0(n), l(n)+2)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 0)*sz(n, 2)*eyg(j(n), k0(n), l(n)+2)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 0)*sz(n, 2)*eyg(j(n)+1, k0(n), l(n)+2)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 0)*sz(n, 2)*eyg(j(n)+2, k0(n), l(n)+2)
      ey(nn) = ey(nn) + sx(n, -1)*sy0(n, 1)*sz(n, 2)*eyg(j(n)-1, k0(n)+1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, 0)*sy0(n, 1)*sz(n, 2)*eyg(j(n), k0(n)+1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, 1)*sy0(n, 1)*sz(n, 2)*eyg(j(n)+1, k0(n)+1, l(n)+2)
      ey(nn) = ey(nn) + sx(n, 2)*sy0(n, 1)*sz(n, 2)*eyg(j(n)+2, k0(n)+1, l(n)+2)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ez:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1
      ! Compute Ez on particle
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, -1)*ezg(j(n)-1, k(n)-1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, -1)*ezg(j(n), k(n)-1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, -1)*ezg(j(n)+1, k(n)-1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, -1)*ezg(j(n)+2, k(n)-1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, -1)*ezg(j(n)-1, k(n), l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, -1)*ezg(j(n), k(n), l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, -1)*ezg(j(n)+1, k(n), l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, -1)*ezg(j(n)+2, k(n), l0(n)-1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, -1)*ezg(j(n)-1, k(n)+1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, -1)*ezg(j(n), k(n)+1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, -1)*ezg(j(n)+1, k(n)+1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, -1)*ezg(j(n)+2, k(n)+1, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, -1)*ezg(j(n)-1, k(n)+2, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, -1)*ezg(j(n), k(n)+2, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, -1)*ezg(j(n)+1, k(n)+2, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, -1)*ezg(j(n)+2, k(n)+2, l0(n)-1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 0)*ezg(j(n)-1, k(n)-1, l0(n))
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 0)*ezg(j(n), k(n)-1, l0(n))
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 0)*ezg(j(n)+1, k(n)-1, l0(n))
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 0)*ezg(j(n)+2, k(n)-1, l0(n))
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 0)*ezg(j(n)-1, k(n), l0(n))
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 0)*ezg(j(n), k(n), l0(n))
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 0)*ezg(j(n)+1, k(n), l0(n))
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 0)*ezg(j(n)+2, k(n), l0(n))
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 0)*ezg(j(n)-1, k(n)+1, l0(n))
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 0)*ezg(j(n), k(n)+1, l0(n))
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 0)*ezg(j(n)+1, k(n)+1, l0(n))
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 0)*ezg(j(n)+2, k(n)+1, l0(n))
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 0)*ezg(j(n)-1, k(n)+2, l0(n))
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 0)*ezg(j(n), k(n)+2, l0(n))
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 0)*ezg(j(n)+1, k(n)+2, l0(n))
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 0)*ezg(j(n)+2, k(n)+2, l0(n))
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, -1)*sz0(n, 1)*ezg(j(n)-1, k(n)-1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, -1)*sz0(n, 1)*ezg(j(n), k(n)-1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, -1)*sz0(n, 1)*ezg(j(n)+1, k(n)-1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, -1)*sz0(n, 1)*ezg(j(n)+2, k(n)-1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 0)*sz0(n, 1)*ezg(j(n)-1, k(n), l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 0)*sz0(n, 1)*ezg(j(n), k(n), l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 0)*sz0(n, 1)*ezg(j(n)+1, k(n), l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 0)*sz0(n, 1)*ezg(j(n)+2, k(n), l0(n)+1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 1)*sz0(n, 1)*ezg(j(n)-1, k(n)+1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 1)*sz0(n, 1)*ezg(j(n), k(n)+1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 1)*sz0(n, 1)*ezg(j(n)+1, k(n)+1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 1)*sz0(n, 1)*ezg(j(n)+2, k(n)+1, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, -1)*sy(n, 2)*sz0(n, 1)*ezg(j(n)-1, k(n)+2, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 0)*sy(n, 2)*sz0(n, 1)*ezg(j(n), k(n)+2, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 1)*sy(n, 2)*sz0(n, 1)*ezg(j(n)+1, k(n)+2, l0(n)+1)
      ez(nn) = ez(nn) + sx(n, 2)*sy(n, 2)*sz0(n, 1)*ezg(j(n)+2, k(n)+2, l0(n)+1)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED bx:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Bx on particle
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, -1)*bxg(j(n)-1, k0(n)-1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, -1)*bxg(j(n), k0(n)-1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, -1)*bxg(j(n)+1, k0(n)-1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, -1)*bxg(j(n)+2, k0(n)-1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, -1)*bxg(j(n)-1, k0(n), l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, -1)*bxg(j(n), k0(n), l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, -1)*bxg(j(n)+1, k0(n), l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, -1)*bxg(j(n)+2, k0(n), l0(n)-1)
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, -1)*bxg(j(n)-1, k0(n)+1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, -1)*bxg(j(n), k0(n)+1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, -1)*bxg(j(n)+1, k0(n)+1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, -1)*bxg(j(n)+2, k0(n)+1, l0(n)-1)
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 0)*bxg(j(n)-1, k0(n)-1, l0(n))
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 0)*bxg(j(n), k0(n)-1, l0(n))
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 0)*bxg(j(n)+1, k0(n)-1, l0(n))
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 0)*bxg(j(n)+2, k0(n)-1, l0(n))
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 0)*bxg(j(n)-1, k0(n), l0(n))
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 0)*bxg(j(n), k0(n), l0(n))
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 0)*bxg(j(n)+1, k0(n), l0(n))
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 0)*bxg(j(n)+2, k0(n), l0(n))
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 0)*bxg(j(n)-1, k0(n)+1, l0(n))
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 0)*bxg(j(n), k0(n)+1, l0(n))
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 0)*bxg(j(n)+1, k0(n)+1, l0(n))
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 0)*bxg(j(n)+2, k0(n)+1, l0(n))
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, -1)*sz0(n, 1)*bxg(j(n)-1, k0(n)-1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, -1)*sz0(n, 1)*bxg(j(n), k0(n)-1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, -1)*sz0(n, 1)*bxg(j(n)+1, k0(n)-1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, -1)*sz0(n, 1)*bxg(j(n)+2, k0(n)-1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 0)*sz0(n, 1)*bxg(j(n)-1, k0(n), l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 0)*sz0(n, 1)*bxg(j(n), k0(n), l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 0)*sz0(n, 1)*bxg(j(n)+1, k0(n), l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 0)*sz0(n, 1)*bxg(j(n)+2, k0(n), l0(n)+1)
      bx(nn) = bx(nn) + sx(n, -1)*sy0(n, 1)*sz0(n, 1)*bxg(j(n)-1, k0(n)+1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 0)*sy0(n, 1)*sz0(n, 1)*bxg(j(n), k0(n)+1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 1)*sy0(n, 1)*sz0(n, 1)*bxg(j(n)+1, k0(n)+1, l0(n)+1)
      bx(nn) = bx(nn) + sx(n, 2)*sy0(n, 1)*sz0(n, 1)*bxg(j(n)+2, k0(n)+1, l0(n)+1)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED by:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1
      ! Compute By on particle
      by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, -1)*byg(j0(n)-1, k(n)-1, l0(n)-1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, -1)*byg(j0(n), k(n)-1, l0(n)-1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, -1)*byg(j0(n)+1, k(n)-1, l0(n)-1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, -1)*byg(j0(n)-1, k(n), l0(n)-1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, -1)*byg(j0(n), k(n), l0(n)-1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, -1)*byg(j0(n)+1, k(n), l0(n)-1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, -1)*byg(j0(n)-1, k(n)+1, l0(n)-1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, -1)*byg(j0(n), k(n)+1, l0(n)-1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, -1)*byg(j0(n)+1, k(n)+1, l0(n)-1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, -1)*byg(j0(n)-1, k(n)+2, l0(n)-1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, -1)*byg(j0(n), k(n)+2, l0(n)-1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, -1)*byg(j0(n)+1, k(n)+2, l0(n)-1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 0)*byg(j0(n)-1, k(n)-1, l0(n))
      by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 0)*byg(j0(n), k(n)-1, l0(n))
      by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 0)*byg(j0(n)+1, k(n)-1, l0(n))
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 0)*byg(j0(n)-1, k(n), l0(n))
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 0)*byg(j0(n), k(n), l0(n))
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 0)*byg(j0(n)+1, k(n), l0(n))
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 0)*byg(j0(n)-1, k(n)+1, l0(n))
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 0)*byg(j0(n), k(n)+1, l0(n))
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 0)*byg(j0(n)+1, k(n)+1, l0(n))
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 0)*byg(j0(n)-1, k(n)+2, l0(n))
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 0)*byg(j0(n), k(n)+2, l0(n))
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 0)*byg(j0(n)+1, k(n)+2, l0(n))
      by(nn) = by(nn) + sx0(n, -1)*sy(n, -1)*sz0(n, 1)*byg(j0(n)-1, k(n)-1, l0(n)+1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, -1)*sz0(n, 1)*byg(j0(n), k(n)-1, l0(n)+1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, -1)*sz0(n, 1)*byg(j0(n)+1, k(n)-1, l0(n)+1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 0)*sz0(n, 1)*byg(j0(n)-1, k(n), l0(n)+1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 0)*sz0(n, 1)*byg(j0(n), k(n), l0(n)+1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 0)*sz0(n, 1)*byg(j0(n)+1, k(n), l0(n)+1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 1)*sz0(n, 1)*byg(j0(n)-1, k(n)+1, l0(n)+1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 1)*sz0(n, 1)*byg(j0(n), k(n)+1, l0(n)+1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 1)*sz0(n, 1)*byg(j0(n)+1, k(n)+1, l0(n)+1)
      by(nn) = by(nn) + sx0(n, -1)*sy(n, 2)*sz0(n, 1)*byg(j0(n)-1, k(n)+2, l0(n)+1)
      by(nn) = by(nn) + sx0(n, 0)*sy(n, 2)*sz0(n, 1)*byg(j0(n), k(n)+2, l0(n)+1)
      by(nn) = by(nn) + sx0(n, 1)*sy(n, 2)*sz0(n, 1)*byg(j0(n)+1, k(n)+2, l0(n)+1)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED bz:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Bz on particle
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, -1)*bzg(j0(n)-1, k0(n)-1, l(n)-1)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, -1)*bzg(j0(n), k0(n)-1, l(n)-1)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, -1)*bzg(j0(n)+1, k0(n)-1, l(n)-1)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, -1)*bzg(j0(n)-1, k0(n), l(n)-1)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, -1)*bzg(j0(n), k0(n), l(n)-1)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, -1)*bzg(j0(n)+1, k0(n), l(n)-1)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, -1)*bzg(j0(n)-1, k0(n)+1, l(n)-1)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, -1)*bzg(j0(n), k0(n)+1, l(n)-1)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, -1)*bzg(j0(n)+1, k0(n)+1, l(n)-1)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 0)*bzg(j0(n)-1, k0(n)-1, l(n))
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 0)*bzg(j0(n), k0(n)-1, l(n))
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 0)*bzg(j0(n)+1, k0(n)-1, l(n))
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 0)*bzg(j0(n)-1, k0(n), l(n))
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 0)*bzg(j0(n), k0(n), l(n))
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 0)*bzg(j0(n)+1, k0(n), l(n))
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 0)*bzg(j0(n)-1, k0(n)+1, l(n))
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 0)*bzg(j0(n), k0(n)+1, l(n))
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 0)*bzg(j0(n)+1, k0(n)+1, l(n))
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 1)*bzg(j0(n)-1, k0(n)-1, l(n)+1)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 1)*bzg(j0(n), k0(n)-1, l(n)+1)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 1)*bzg(j0(n)+1, k0(n)-1, l(n)+1)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 1)*bzg(j0(n)-1, k0(n), l(n)+1)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 1)*bzg(j0(n), k0(n), l(n)+1)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 1)*bzg(j0(n)+1, k0(n), l(n)+1)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 1)*bzg(j0(n)-1, k0(n)+1, l(n)+1)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 1)*bzg(j0(n), k0(n)+1, l(n)+1)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 1)*bzg(j0(n)+1, k0(n)+1, l(n)+1)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, -1)*sz(n, 2)*bzg(j0(n)-1, k0(n)-1, l(n)+2)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, -1)*sz(n, 2)*bzg(j0(n), k0(n)-1, l(n)+2)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, -1)*sz(n, 2)*bzg(j0(n)+1, k0(n)-1, l(n)+2)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 0)*sz(n, 2)*bzg(j0(n)-1, k0(n), l(n)+2)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 0)*sz(n, 2)*bzg(j0(n), k0(n), l(n)+2)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 0)*sz(n, 2)*bzg(j0(n)+1, k0(n), l(n)+2)
      bz(nn) = bz(nn) + sx0(n, -1)*sy0(n, 1)*sz(n, 2)*bzg(j0(n)-1, k0(n)+1, l(n)+2)
      bz(nn) = bz(nn) + sx0(n, 0)*sy0(n, 1)*sz(n, 2)*bzg(j0(n), k0(n)+1, l(n)+2)
      bz(nn) = bz(nn) + sx0(n, 1)*sy0(n, 1)*sz(n, 2)*bzg(j0(n)+1, k0(n)+1, l(n)+2)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

  ENDDO

  RETURN
END SUBROUTINE
#endif

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Field gathering by block at order 3 with gathering of E and B
!> splited in several small subloops, version 2.
!
!> @detail
!> This function is vectorized and is the version 2 using a more optimized
!> gathering method more efficient for the vectorization.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 12/03/2016
!
!> @warning
!> Only l_lower_order_in_v=True is implemented
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
SUBROUTINE geteb3d_energy_conserving_blockvec2_3_3_3(np, xp, yp, zp, ex, ey, ez, bx,  &
  by, bz, xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, exg,     &
  eyg, ezg, bxg, byg, bzg, lvect, l_lower_order_in_v, l_nodal)
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _________________________________________________
  INTEGER(idp)                           :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(np)               :: xp, yp, zp, ex, ey, ez, bx, by, bz
  INTEGER(idp)                           :: lvect
  LOGICAL(lp)                            :: l_lower_order_in_v, l_nodal
  REAL(num)                              :: stagger_shift
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  REAL(num)                              :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(isp)                           :: ip
  INTEGER(isp), DIMENSION(lvect)         :: j, k, l
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j, k, l
#endif
  INTEGER(isp), DIMENSION(lvect)         :: j0, k0, l0
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: j0, k0, l0
#endif
  INTEGER(idp)                           :: jj, kk, ll
  REAL(num)                              :: dxi, dyi, dzi, x, y, z, xint, yint, zint
  REAL(num)                              :: xintsq, oxint, yintsq, oyint, zintsq
  REAL(num)                              :: ozint, oxintsq, oyintsq, ozintsq
  REAL(num)                              :: a
  INTEGER(isp)                           :: nn, n
  REAL(num), DIMENSION(lvect, -1:2)       :: sx
  REAL(num), DIMENSION(lvect, -1:2)       :: sy
  REAL(num), DIMENSION(lvect, -1:2)       :: sz
  REAL(num), DIMENSION(lvect, -1:1)       :: sx0
  REAL(num), DIMENSION(lvect, -1:1)       :: sy0
  REAL(num), DIMENSION(lvect, -1:1)       :: sz0
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

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num

  ! ___ Loop on partciles _______________________
  DO ip=1, np, lvect

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
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(n, -1) = onesixth*oxintsq*oxint
      sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx(n, 2) = onesixth*xintsq*xint

      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(n, -1) = onesixth*oyintsq*oyint
      sy(n, 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy(n, 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy(n, 2) = onesixth*yintsq*yint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(n, -1) = onesixth*ozintsq*ozint
      sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz(n, 2) = onesixth*zintsq*zint

      xint=x-stagger_shift-j0(n)
      yint=y-stagger_shift-k0(n)
      zint=z-stagger_shift-l0(n)

      xintsq = xint*xint
      sx0(n, -1) = 0.5_num*(0.5_num-xint)**2
      sx0(n, 0) = 0.75_num-xintsq
      sx0(n, 1) = 0.5_num*(0.5_num+xint)**2

      yintsq = yint*yint
      sy0(n, -1) = 0.5_num*(0.5_num-yint)**2
      sy0(n, 0) = 0.75_num-yintsq
      sy0(n, 1) = 0.5_num*(0.5_num+yint)**2

      zintsq = zint*zint
      sz0(n, -1) = 0.5_num*(0.5_num-zint)**2
      sz0(n, 0) = 0.75_num-zintsq
      sz0(n, 1) = 0.5_num*(0.5_num+zint)**2
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ex:64
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

      ! Compute Ex on particle
      a = (sx0(n, -1)*exg(j0(n)-1, k(n)-1, l(n)-1) + sx0(n, 0)*exg(j0(n), k(n)-1,     &
      l(n)-1) + sx0(n, 1)*exg(j0(n)+1, k(n)-1, l(n)-1))*sy(n, -1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n), l(n)-1) + sx0(n, 0)*exg(j0(n), k(n),     &
      l(n)-1) + sx0(n, 1)*exg(j0(n)+1, k(n), l(n)-1))*sy(n, 0)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+1, l(n)-1) + sx0(n, 0)*exg(j0(n), k(n)+1, &
      l(n)-1) + sx0(n, 1)*exg(j0(n)+1, k(n)+1, l(n)-1))*sy(n, 1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+2, l(n)-1) + sx0(n, 0)*exg(j0(n), k(n)+2, &
      l(n)-1) + sx0(n, 1)*exg(j0(n)+1, k(n)+2, l(n)-1))*sy(n, 2)
      ex(nn) = ex(nn) + a*sz(n, -1)
      a = (sx0(n, -1)*exg(j0(n)-1, k(n)-1, l(n)) + sx0(n, 0)*exg(j0(n), k(n)-1, l(n)) &
      + sx0(n, 1)*exg(j0(n)+1, k(n)-1, l(n)))*sy(n, -1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n), l(n)) + sx0(n, 0)*exg(j0(n), k(n), l(n)) &
      + sx0(n, 1)*exg(j0(n)+1, k(n), l(n)))*sy(n, 0)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+1, l(n)) + sx0(n, 0)*exg(j0(n), k(n)+1,   &
      l(n)) + sx0(n, 1)*exg(j0(n)+1, k(n)+1, l(n)))*sy(n, 1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+2, l(n)) + sx0(n, 0)*exg(j0(n), k(n)+2,   &
      l(n)) + sx0(n, 1)*exg(j0(n)+1, k(n)+2, l(n)))*sy(n, 2)
      ex(nn) = ex(nn) + a*sz(n, 0)
      a = (sx0(n, -1)*exg(j0(n)-1, k(n)-1, l(n)+1) + sx0(n, 0)*exg(j0(n), k(n)-1,     &
      l(n)+1) + sx0(n, 1)*exg(j0(n)+1, k(n)-1, l(n)+1))*sy(n, -1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n), l(n)+1) + sx0(n, 0)*exg(j0(n), k(n),     &
      l(n)+1) + sx0(n, 1)*exg(j0(n)+1, k(n), l(n)+1))*sy(n, 0)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+1, l(n)+1) + sx0(n, 0)*exg(j0(n), k(n)+1, &
      l(n)+1) + sx0(n, 1)*exg(j0(n)+1, k(n)+1, l(n)+1))*sy(n, 1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+2, l(n)+1) + sx0(n, 0)*exg(j0(n), k(n)+2, &
      l(n)+1) + sx0(n, 1)*exg(j0(n)+1, k(n)+2, l(n)+1))*sy(n, 2)
      ex(nn) = ex(nn) + a*sz(n, 1)
      a = (sx0(n, -1)*exg(j0(n)-1, k(n)-1, l(n)+2) + sx0(n, 0)*exg(j0(n), k(n)-1,     &
      l(n)+2) + sx0(n, 1)*exg(j0(n)+1, k(n)-1, l(n)+2))*sy(n, -1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n), l(n)+2) + sx0(n, 0)*exg(j0(n), k(n),     &
      l(n)+2) + sx0(n, 1)*exg(j0(n)+1, k(n), l(n)+2))*sy(n, 0)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+1, l(n)+2) + sx0(n, 0)*exg(j0(n), k(n)+1, &
      l(n)+2) + sx0(n, 1)*exg(j0(n)+1, k(n)+1, l(n)+2))*sy(n, 1)
      a = a + (sx0(n, -1)*exg(j0(n)-1, k(n)+2, l(n)+2) + sx0(n, 0)*exg(j0(n), k(n)+2, &
      l(n)+2) + sx0(n, 1)*exg(j0(n)+1, k(n)+2, l(n)+2))*sy(n, 2)
      ex(nn) = ex(nn) + a*sz(n, 2)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ey:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Ey on particle
      a = (sx(n, -1)*eyg(j(n)-1, k0(n)-1, l(n)-1) + sx(n, 0)*eyg(j(n), k0(n)-1,       &
      l(n)-1) + sx(n, 1)*eyg(j(n)+1, k0(n)-1, l(n)-1) + sx(n, 2)*eyg(j(n)+2, k0(n)-1, &
      l(n)-1))*sy0(n, -1)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n), l(n)-1) + sx(n, 0)*eyg(j(n), k0(n),       &
      l(n)-1) + sx(n, 1)*eyg(j(n)+1, k0(n), l(n)-1) + sx(n, 2)*eyg(j(n)+2, k0(n),     &
      l(n)-1))*sy0(n, 0)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n)+1, l(n)-1) + sx(n, 0)*eyg(j(n), k0(n)+1,   &
      l(n)-1) + sx(n, 1)*eyg(j(n)+1, k0(n)+1, l(n)-1) + sx(n, 2)*eyg(j(n)+2, k0(n)+1, &
      l(n)-1))*sy0(n, 1)
      ey(nn) = ey(nn) + a*sz(n, -1)
      a = (sx(n, -1)*eyg(j(n)-1, k0(n)-1, l(n)) + sx(n, 0)*eyg(j(n), k0(n)-1, l(n)) + &
      sx(n, 1)*eyg(j(n)+1, k0(n)-1, l(n)) + sx(n, 2)*eyg(j(n)+2, k0(n)-1,             &
      l(n)))*sy0(n, -1)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n), l(n)) + sx(n, 0)*eyg(j(n), k0(n), l(n)) + &
      sx(n, 1)*eyg(j(n)+1, k0(n), l(n)) + sx(n, 2)*eyg(j(n)+2, k0(n), l(n)))*sy0(n,   &
      0)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n)+1, l(n)) + sx(n, 0)*eyg(j(n), k0(n)+1,     &
      l(n)) + sx(n, 1)*eyg(j(n)+1, k0(n)+1, l(n)) + sx(n, 2)*eyg(j(n)+2, k0(n)+1,     &
      l(n)))*sy0(n, 1)
      ey(nn) = ey(nn) + a*sz(n, 0)
      a = (sx(n, -1)*eyg(j(n)-1, k0(n)-1, l(n)+1) + sx(n, 0)*eyg(j(n), k0(n)-1,       &
      l(n)+1) + sx(n, 1)*eyg(j(n)+1, k0(n)-1, l(n)+1) + sx(n, 2)*eyg(j(n)+2, k0(n)-1, &
      l(n)+1))*sy0(n, -1)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n), l(n)+1) + sx(n, 0)*eyg(j(n), k0(n),       &
      l(n)+1) + sx(n, 1)*eyg(j(n)+1, k0(n), l(n)+1) + sx(n, 2)*eyg(j(n)+2, k0(n),     &
      l(n)+1))*sy0(n, 0)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n)+1, l(n)+1) + sx(n, 0)*eyg(j(n), k0(n)+1,   &
      l(n)+1) + sx(n, 1)*eyg(j(n)+1, k0(n)+1, l(n)+1) + sx(n, 2)*eyg(j(n)+2, k0(n)+1, &
      l(n)+1))*sy0(n, 1)
      ey(nn) = ey(nn) + a*sz(n, 1)
      a = (sx(n, -1)*eyg(j(n)-1, k0(n)-1, l(n)+2) + sx(n, 0)*eyg(j(n), k0(n)-1,       &
      l(n)+2) + sx(n, 1)*eyg(j(n)+1, k0(n)-1, l(n)+2) + sx(n, 2)*eyg(j(n)+2, k0(n)-1, &
      l(n)+2))*sy0(n, -1)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n), l(n)+2) + sx(n, 0)*eyg(j(n), k0(n),       &
      l(n)+2) + sx(n, 1)*eyg(j(n)+1, k0(n), l(n)+2) + sx(n, 2)*eyg(j(n)+2, k0(n),     &
      l(n)+2))*sy0(n, 0)
      a = a + (sx(n, -1)*eyg(j(n)-1, k0(n)+1, l(n)+2) + sx(n, 0)*eyg(j(n), k0(n)+1,   &
      l(n)+2) + sx(n, 1)*eyg(j(n)+1, k0(n)+1, l(n)+2) + sx(n, 2)*eyg(j(n)+2, k0(n)+1, &
      l(n)+2))*sy0(n, 1)
      ey(nn) = ey(nn) + a*sz(n, 2)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED ez:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Ez on particle
      a = (sx(n, -1)*ezg(j(n)-1, k(n)-1, l0(n)-1) + sx(n, 0)*ezg(j(n), k(n)-1,        &
      l0(n)-1) + sx(n, 1)*ezg(j(n)+1, k(n)-1, l0(n)-1) + sx(n, 2)*ezg(j(n)+2, k(n)-1, &
      l0(n)-1))*sy(n, -1)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n), l0(n)-1) + sx(n, 0)*ezg(j(n), k(n),        &
      l0(n)-1) + sx(n, 1)*ezg(j(n)+1, k(n), l0(n)-1) + sx(n, 2)*ezg(j(n)+2, k(n),     &
      l0(n)-1))*sy(n, 0)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n)+1, l0(n)-1) + sx(n, 0)*ezg(j(n), k(n)+1,    &
      l0(n)-1) + sx(n, 1)*ezg(j(n)+1, k(n)+1, l0(n)-1) + sx(n, 2)*ezg(j(n)+2, k(n)+1, &
      l0(n)-1))*sy(n, 1)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n)+2, l0(n)-1) + sx(n, 0)*ezg(j(n), k(n)+2,    &
      l0(n)-1) + sx(n, 1)*ezg(j(n)+1, k(n)+2, l0(n)-1) + sx(n, 2)*ezg(j(n)+2, k(n)+2, &
      l0(n)-1))*sy(n, 2)
      ez(nn) = ez(nn) + a*sz0(n, -1)
      a = (sx(n, -1)*ezg(j(n)-1, k(n)-1, l0(n)) + sx(n, 0)*ezg(j(n), k(n)-1, l0(n)) + &
      sx(n, 1)*ezg(j(n)+1, k(n)-1, l0(n)) + sx(n, 2)*ezg(j(n)+2, k(n)-1,              &
      l0(n)))*sy(n, -1)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n), l0(n)) + sx(n, 0)*ezg(j(n), k(n), l0(n)) + &
      sx(n, 1)*ezg(j(n)+1, k(n), l0(n)) + sx(n, 2)*ezg(j(n)+2, k(n), l0(n)))*sy(n, 0)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n)+1, l0(n)) + sx(n, 0)*ezg(j(n), k(n)+1,      &
      l0(n)) + sx(n, 1)*ezg(j(n)+1, k(n)+1, l0(n)) + sx(n, 2)*ezg(j(n)+2, k(n)+1,     &
      l0(n)))*sy(n, 1)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n)+2, l0(n)) + sx(n, 0)*ezg(j(n), k(n)+2,      &
      l0(n)) + sx(n, 1)*ezg(j(n)+1, k(n)+2, l0(n)) + sx(n, 2)*ezg(j(n)+2, k(n)+2,     &
      l0(n)))*sy(n, 2)
      ez(nn) = ez(nn) + a*sz0(n, 0)
      a = (sx(n, -1)*ezg(j(n)-1, k(n)-1, l0(n)+1) + sx(n, 0)*ezg(j(n), k(n)-1,        &
      l0(n)+1) + sx(n, 1)*ezg(j(n)+1, k(n)-1, l0(n)+1) + sx(n, 2)*ezg(j(n)+2, k(n)-1, &
      l0(n)+1))*sy(n, -1)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n), l0(n)+1) + sx(n, 0)*ezg(j(n), k(n),        &
      l0(n)+1) + sx(n, 1)*ezg(j(n)+1, k(n), l0(n)+1) + sx(n, 2)*ezg(j(n)+2, k(n),     &
      l0(n)+1))*sy(n, 0)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n)+1, l0(n)+1) + sx(n, 0)*ezg(j(n), k(n)+1,    &
      l0(n)+1) + sx(n, 1)*ezg(j(n)+1, k(n)+1, l0(n)+1) + sx(n, 2)*ezg(j(n)+2, k(n)+1, &
      l0(n)+1))*sy(n, 1)
      a = a + (sx(n, -1)*ezg(j(n)-1, k(n)+2, l0(n)+1) + sx(n, 0)*ezg(j(n), k(n)+2,    &
      l0(n)+1) + sx(n, 1)*ezg(j(n)+1, k(n)+2, l0(n)+1) + sx(n, 2)*ezg(j(n)+2, k(n)+2, &
      l0(n)+1))*sy(n, 2)
      ez(nn) = ez(nn) + a*sz0(n, 1)


    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED bx:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Bx on particle
      a = (sx(n, -1)*bxg(j(n)-1, k0(n)-1, l0(n)-1) + sx(n, 0)*bxg(j(n), k0(n)-1,      &
      l0(n)-1) + sx(n, 1)*bxg(j(n)+1, k0(n)-1, l0(n)-1) + sx(n, 2)*bxg(j(n)+2,        &
      k0(n)-1, l0(n)-1))*sy0(n, -1)
      a = a + (sx(n, -1)*bxg(j(n)-1, k0(n), l0(n)-1) + sx(n, 0)*bxg(j(n), k0(n),      &
      l0(n)-1) + sx(n, 1)*bxg(j(n)+1, k0(n), l0(n)-1) + sx(n, 2)*bxg(j(n)+2, k0(n),   &
      l0(n)-1))*sy0(n, 0)
      a = a + (sx(n, -1)*bxg(j(n)-1, k0(n)+1, l0(n)-1) + sx(n, 0)*bxg(j(n), k0(n)+1,  &
      l0(n)-1) + sx(n, 1)*bxg(j(n)+1, k0(n)+1, l0(n)-1) + sx(n, 2)*bxg(j(n)+2,        &
      k0(n)+1, l0(n)-1))*sy0(n, 1)
      bx(nn) = bx(nn) + a*sz0(n, -1)
      a = (sx(n, -1)*bxg(j(n)-1, k0(n)-1, l0(n)) + sx(n, 0)*bxg(j(n), k0(n)-1, l0(n)) &
      + sx(n, 1)*bxg(j(n)+1, k0(n)-1, l0(n)) + sx(n, 2)*bxg(j(n)+2, k0(n)-1,          &
      l0(n)))*sy0(n, -1)
      a = a + (sx(n, -1)*bxg(j(n)-1, k0(n), l0(n)) + sx(n, 0)*bxg(j(n), k0(n), l0(n)) &
      + sx(n, 1)*bxg(j(n)+1, k0(n), l0(n)) + sx(n, 2)*bxg(j(n)+2, k0(n),              &
      l0(n)))*sy0(n, 0)
      a = a + (sx(n, -1)*bxg(j(n)-1, k0(n)+1, l0(n)) + sx(n, 0)*bxg(j(n), k0(n)+1,    &
      l0(n)) + sx(n, 1)*bxg(j(n)+1, k0(n)+1, l0(n)) + sx(n, 2)*bxg(j(n)+2, k0(n)+1,   &
      l0(n)))*sy0(n, 1)
      bx(nn) = bx(nn) + a*sz0(n, 0)
      a = (sx(n, -1)*bxg(j(n)-1, k0(n)-1, l0(n)+1) + sx(n, 0)*bxg(j(n), k0(n)-1,      &
      l0(n)+1) + sx(n, 1)*bxg(j(n)+1, k0(n)-1, l0(n)+1) + sx(n, 2)*bxg(j(n)+2,        &
      k0(n)-1, l0(n)+1))*sy0(n, -1)
      a = a + (sx(n, -1)*bxg(j(n)-1, k0(n), l0(n)+1) + sx(n, 0)*bxg(j(n), k0(n),      &
      l0(n)+1) + sx(n, 1)*bxg(j(n)+1, k0(n), l0(n)+1) + sx(n, 2)*bxg(j(n)+2, k0(n),   &
      l0(n)+1))*sy0(n, 0)
      a = a + (sx(n, -1)*bxg(j(n)-1, k0(n)+1, l0(n)+1) + sx(n, 0)*bxg(j(n), k0(n)+1,  &
      l0(n)+1) + sx(n, 1)*bxg(j(n)+1, k0(n)+1, l0(n)+1) + sx(n, 2)*bxg(j(n)+2,        &
      k0(n)+1, l0(n)+1))*sy0(n, 1)
      bx(nn) = bx(nn) + a*sz0(n, 1)

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED by:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute By on particle
      a = (sx0(n, -1)*byg(j0(n)-1, k(n)-1, l0(n)-1) + sx0(n, 0)*byg(j0(n), k(n)-1,    &
      l0(n)-1) + sx0(n, 1)*byg(j0(n)+1, k(n)-1, l0(n)-1))*sy(n, -1)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n), l0(n)-1) + sx0(n, 0)*byg(j0(n), k(n),    &
      l0(n)-1) + sx0(n, 1)*byg(j0(n)+1, k(n), l0(n)-1))*sy(n, 0)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n)+1, l0(n)-1) + sx0(n, 0)*byg(j0(n),        &
      k(n)+1, l0(n)-1) + sx0(n, 1)*byg(j0(n)+1, k(n)+1, l0(n)-1))*sy(n, 1)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n)+2, l0(n)-1) + sx0(n, 0)*byg(j0(n),        &
      k(n)+2, l0(n)-1) + sx0(n, 1)*byg(j0(n)+1, k(n)+2, l0(n)-1))*sy(n, 2)
      by(nn) = by(nn) + a*sz0(n, -1)
      a = (sx0(n, -1)*byg(j0(n)-1, k(n)-1, l0(n)) + sx0(n, 0)*byg(j0(n), k(n)-1,      &
      l0(n)) + sx0(n, 1)*byg(j0(n)+1, k(n)-1, l0(n)))*sy(n, -1)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n), l0(n)) + sx0(n, 0)*byg(j0(n), k(n),      &
      l0(n)) + sx0(n, 1)*byg(j0(n)+1, k(n), l0(n)))*sy(n, 0)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n)+1, l0(n)) + sx0(n, 0)*byg(j0(n), k(n)+1,  &
      l0(n)) + sx0(n, 1)*byg(j0(n)+1, k(n)+1, l0(n)))*sy(n, 1)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n)+2, l0(n)) + sx0(n, 0)*byg(j0(n), k(n)+2,  &
      l0(n)) + sx0(n, 1)*byg(j0(n)+1, k(n)+2, l0(n)))*sy(n, 2)
      by(nn) = by(nn) + a*sz0(n, 0)
      a = (sx0(n, -1)*byg(j0(n)-1, k(n)-1, l0(n)+1) + sx0(n, 0)*byg(j0(n), k(n)-1,    &
      l0(n)+1) + sx0(n, 1)*byg(j0(n)+1, k(n)-1, l0(n)+1))*sy(n, -1)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n), l0(n)+1) + sx0(n, 0)*byg(j0(n), k(n),    &
      l0(n)+1) + sx0(n, 1)*byg(j0(n)+1, k(n), l0(n)+1))*sy(n, 0)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n)+1, l0(n)+1) + sx0(n, 0)*byg(j0(n),        &
      k(n)+1, l0(n)+1) + sx0(n, 1)*byg(j0(n)+1, k(n)+1, l0(n)+1))*sy(n, 1)
      a = a + (sx0(n, -1)*byg(j0(n)-1, k(n)+2, l0(n)+1) + sx0(n, 0)*byg(j0(n),        &
      k(n)+2, l0(n)+1) + sx0(n, 1)*byg(j0(n)+1, k(n)+2, l0(n)+1))*sy(n, 2)
      by(nn) = by(nn) + a*sz0(n, 1)


    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED bz:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED j:64, k:64, l:64
    !DIR$ ASSUME_ALIGNED j0:64, k0:64, l0:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, bz)
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
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1

      ! Compute Bz on particle
      a = (sx0(n, -1)*bzg(j0(n)-1, k0(n)-1, l(n)-1) + sx0(n, 0)*bzg(j0(n), k0(n)-1,   &
      l(n)-1) + sx0(n, 1)*bzg(j0(n)+1, k0(n)-1, l(n)-1))*sy0(n, -1)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n), l(n)-1) + sx0(n, 0)*bzg(j0(n), k0(n),   &
      l(n)-1) + sx0(n, 1)*bzg(j0(n)+1, k0(n), l(n)-1))*sy0(n, 0)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n)+1, l(n)-1) + sx0(n, 0)*bzg(j0(n),        &
      k0(n)+1, l(n)-1) + sx0(n, 1)*bzg(j0(n)+1, k0(n)+1, l(n)-1))*sy0(n, 1)
      bz(nn) = bz(nn) + a*sz(n, -1)
      a = (sx0(n, -1)*bzg(j0(n)-1, k0(n)-1, l(n)) + sx0(n, 0)*bzg(j0(n), k0(n)-1,     &
      l(n)) + sx0(n, 1)*bzg(j0(n)+1, k0(n)-1, l(n)))*sy0(n, -1)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n), l(n)) + sx0(n, 0)*bzg(j0(n), k0(n),     &
      l(n)) + sx0(n, 1)*bzg(j0(n)+1, k0(n), l(n)))*sy0(n, 0)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n)+1, l(n)) + sx0(n, 0)*bzg(j0(n), k0(n)+1, &
      l(n)) + sx0(n, 1)*bzg(j0(n)+1, k0(n)+1, l(n)))*sy0(n, 1)
      bz(nn) = bz(nn) + a*sz(n, 0)
      a = (sx0(n, -1)*bzg(j0(n)-1, k0(n)-1, l(n)+1) + sx0(n, 0)*bzg(j0(n), k0(n)-1,   &
      l(n)+1) + sx0(n, 1)*bzg(j0(n)+1, k0(n)-1, l(n)+1))*sy0(n, -1)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n), l(n)+1) + sx0(n, 0)*bzg(j0(n), k0(n),   &
      l(n)+1) + sx0(n, 1)*bzg(j0(n)+1, k0(n), l(n)+1))*sy0(n, 0)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n)+1, l(n)+1) + sx0(n, 0)*bzg(j0(n),        &
      k0(n)+1, l(n)+1) + sx0(n, 1)*bzg(j0(n)+1, k0(n)+1, l(n)+1))*sy0(n, 1)
      bz(nn) = bz(nn) + a*sz(n, 1)
      a = (sx0(n, -1)*bzg(j0(n)-1, k0(n)-1, l(n)+2) + sx0(n, 0)*bzg(j0(n), k0(n)-1,   &
      l(n)+2) + sx0(n, 1)*bzg(j0(n)+1, k0(n)-1, l(n)+2))*sy0(n, -1)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n), l(n)+2) + sx0(n, 0)*bzg(j0(n), k0(n),   &
      l(n)+2) + sx0(n, 1)*bzg(j0(n)+1, k0(n), l(n)+2))*sy0(n, 0)
      a = a + (sx0(n, -1)*bzg(j0(n)-1, k0(n)+1, l(n)+2) + sx0(n, 0)*bzg(j0(n),        &
      k0(n)+1, l(n)+2) + sx0(n, 1)*bzg(j0(n)+1, k0(n)+1, l(n)+2))*sy0(n, 1)
      bz(nn) = bz(nn) + a*sz(n, 2)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
  ENDDO
  RETURN
END SUBROUTINE
#endif

! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
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
! FIELD_GATHERING_MANAGER_3D.F90
!
! This file contains subroutines for the field gathering in 3d
! at arbitrary order.
!
! Developers:
! - Henri vincenti
! - Mathieu Lobet
!
! options:
! - DEV: activates developer's secret subroutines
! - DEBUG: activates DEBUG prints and outputs
!
! List of subroutines:
!
! - pxr_getb3d_n_energy_conserving
! - pxr_gete3d_n_energy_conserving
! ________________________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Gathering of electric field from Yee grid ("energy conserving") on particles
!> at arbitrary order. WARNING: Highly unoptimized routine
!
!> @details
!> This subroutine is inherited from Warp
!
!> @author
!> From Warp
!
!> @date
!> 2015
! ________________________________________________________________________________________
SUBROUTINE pxrgete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,      &
  zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,        &
  eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE omp_lib
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp)             :: np
  INTEGER(idp), intent(in) :: exg_nguard(3), exg_nvalid(3), eyg_nguard(3),            &
  eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3)
  INTEGER(idp)             :: nox, noy, noz
  REAL(num), dimension(np) :: xp, yp, zp, ex, ey, ez
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1,                                       &
  -exg_nguard(3):exg_nvalid(3)+exg_nguard(3)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1,                                       &
  -eyg_nguard(3):eyg_nvalid(3)+eyg_nguard(3)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1,                                       &
  -ezg_nguard(3):ezg_nvalid(3)+ezg_nguard(3)-1)
  LOGICAL(lp)              :: l_lower_order_in_v, l_nodal
  REAL(num)                :: stagger_shift
  REAL(num) :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, ixmin0,      &
  ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
  REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, xintsq, oxint, yintsq,       &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq, signx, signy
  REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  REAL(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
  REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  REAL(num), dimension(:), allocatable :: sx0, sy0, sz0
  REAL(num), parameter :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num


  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)-1
  iymin = -int(noy/2)
  iymax =  int((noy+1)/2)-1
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)-1

  IF (l_lower_order_in_v) THEN
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
  ELSE
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
  END IF

  ALLOCATE(sx0(ixmin0:ixmax0), sy0(iymin0:iymax0), sz0(izmin0:izmax0))

  signx = 1.0_num
  signy = 1.0_num
  DO ip=1, np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    IF (l_lower_order_in_v) THEN
      IF (nox==2*(nox/2)) THEN
        j=nint(x)
        j0=floor(x-stagger_shift)
      ELSE
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
      END IF
      IF (noy==2*(noy/2)) THEN
        k=nint(y)
        k0=floor(y-stagger_shift)
      ELSE
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
      END IF
      IF (noz==2*(noz/2)) THEN
        l=nint(z)
        l0=floor(z-stagger_shift)
      ELSE
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)
      END IF
    ELSE
      IF (nox==2*(nox/2)) THEN
        j=nint(x)
        j0=floor(x+0.5_num-stagger_shift)
      ELSE
        j=floor(x)
        j0=floor(x-stagger_shift)
      END IF
      IF (noy==2*(noy/2)) THEN
        k=nint(y)
        k0=floor(y+0.5_num-stagger_shift)
      ELSE
        k=floor(y)
        k0=floor(y-stagger_shift)
      END IF
      IF (noz==2*(noz/2)) THEN
        l=nint(z)
        l0=floor(z+0.5_num-stagger_shift)
      ELSE
        l=floor(z)
        l0=floor(z-stagger_shift)
      END IF
    END IF

    xint=x-j
    yint=y-k
    zint=z-l

    IF (nox==1) THEN
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
    ELSEIF (nox==2) THEN
      xintsq = xint*xint
      sx(-1) = 0.5_num*(0.5_num-xint)**2
      sx( 0) = 0.75_num-xintsq
      sx( 1) = 0.5_num*(0.5_num+xint)**2
    ELSEIF (nox==3) THEN
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint
    END IF

    IF (noy==1) THEN
      sy( 0) = 1.0_num-yint
      sy( 1) = yint
    ELSEIF (noy==2) THEN
      yintsq = yint*yint
      sy(-1) = 0.5_num*(0.5_num-yint)**2
      sy( 0) = 0.75_num-yintsq
      sy( 1) = 0.5_num*(0.5_num+yint)**2
    ELSEIF (noy==3) THEN
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint
    END IF

    IF (noz==1) THEN
      sz( 0) = 1.0_num-zint
      sz( 1) = zint
    ELSEIF (noz==2) THEN
      zintsq = zint*zint
      sz(-1) = 0.5_num*(0.5_num-zint)**2
      sz( 0) = 0.75_num-zintsq
      sz( 1) = 0.5_num*(0.5_num+zint)**2
    ELSEIF (noz==3) THEN
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint
    END IF

    xint=x-stagger_shift-j0
    yint=y-stagger_shift-k0
    zint=z-stagger_shift-l0

    IF (l_lower_order_in_v) THEN

      IF (nox==1) THEN
        sx0( 0) = 1.0_num
      ELSEIF (nox==2) THEN
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
      ELSEIF (nox==3) THEN
        xintsq = xint*xint
        sx0(-1) = 0.5_num*(0.5_num-xint)**2
        sx0( 0) = 0.75_num-xintsq
        sx0( 1) = 0.5_num*(0.5_num+xint)**2
      END IF

      IF (noy==1) THEN
        sy0( 0) = 1.0_num
      ELSEIF (noy==2) THEN
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
      ELSEIF (noy==3) THEN
        yintsq = yint*yint
        sy0(-1) = 0.5_num*(0.5_num-yint)**2
        sy0( 0) = 0.75_num-yintsq
        sy0( 1) = 0.5_num*(0.5_num+yint)**2
      END IF

      IF (noz==1) THEN
        sz0( 0) = 1.0_num
      ELSEIF (noz==2) THEN
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint
      ELSEIF (noz==3) THEN
        zintsq = zint*zint
        sz0(-1) = 0.5_num*(0.5_num-zint)**2
        sz0( 0) = 0.75_num-zintsq
        sz0( 1) = 0.5_num*(0.5_num+zint)**2
      END IF

    ELSE

      IF (nox==1) THEN
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
      ELSEIF (nox==2) THEN
        xintsq = xint*xint
        sx0(-1) = 0.5_num*(0.5_num-xint)**2
        sx0( 0) = 0.75_num-xintsq
        sx0( 1) = 0.5_num*(0.5_num+xint)**2
      ELSEIF (nox==3) THEN
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(-1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0( 2) = onesixth*xintsq*xint
      END IF

      IF (noy==1) THEN
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
      ELSEIF (noy==2) THEN
        yintsq = yint*yint
        sy0(-1) = 0.5_num*(0.5_num-yint)**2
        sy0( 0) = 0.75_num-yintsq
        sy0( 1) = 0.5_num*(0.5_num+yint)**2
      ELSEIF (noy==3) THEN
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0( 2) = onesixth*yintsq*yint
      END IF

      IF (noz==1) THEN
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint
      ELSEIF (noz==2) THEN
        zintsq = zint*zint
        sz0(-1) = 0.5_num*(0.5_num-zint)**2
        sz0( 0) = 0.75_num-zintsq
        sz0( 1) = 0.5_num*(0.5_num+zint)**2
      ELSEIF (noz==3) THEN
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0( 2) = onesixth*zintsq*zint
      END IF

    END IF

    DO ll = izmin, izmax+1
      DO kk = iymin, iymax+1
        DO jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj, k+kk, l+ll)*signx
        END DO
      END DO
    END DO

    DO ll = izmin, izmax+1
      DO kk = iymin0, iymax0
        DO jj = ixmin, ixmax+1
          ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj, k0+kk, l+ll)*signy
        END DO
      END DO
    END DO

    DO ll = izmin0, izmax0
      DO kk = iymin, iymax+1
        DO jj = ixmin, ixmax+1
          ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj, k+kk, l0+ll)
        END DO
      END DO
    END DO

  END DO
  DEALLOCATE(sx0, sy0, sz0)

  RETURN
END SUBROUTINE pxrgete3d_n_energy_conserving

! ________________________________________________________________________________________
!> @brief
!> Gathering of Magnetic field from Yee grid ("energy conserving") on particles
!> At arbitrary order. WARNING: Highly unoptimized routine
! ________________________________________________________________________________________
SUBROUTINE pxrgetb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,      &
  zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,        &
  byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE omp_lib
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nox, noy, noz
  INTEGER(idp), intent(in)             :: bxg_nguard(3), bxg_nvalid(3),               &
  byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  REAL(num), DIMENSION(np) :: xp, yp, zp, bx, by, bz
  LOGICAL(lp)  :: l_lower_order_in_v, l_nodal
  REAL(num)    :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1,                                       &
  -bxg_nguard(3):bxg_nvalid(3)+bxg_nguard(3)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1,                                       &
  -byg_nguard(3):byg_nvalid(3)+byg_nguard(3)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1,                                       &
  -bzg_nguard(3):bzg_nvalid(3)+bzg_nguard(3)-1)
  REAL(num) :: xmin, ymin, zmin, dx, dy, dz
  INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, ixmin0,      &
  ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
  REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, xintsq, oxint, yintsq,       &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq, signx, signy
  REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  REAL(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
  REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  REAL(num), DIMENSION(:), ALLOCATABLE :: sx0, sy0, sz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)-1
  iymin = -int(noy/2)
  iymax =  int((noy+1)/2)-1
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)-1


  IF (l_lower_order_in_v) THEN
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
  ELSE
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
  END IF
  ALLOCATE(sx0(ixmin0:ixmax0), sy0(iymin0:iymax0), sz0(izmin0:izmax0))

  signx = 1.0_num
  signy = 1.0_num

  sx=0.0_num
  sy=0.0_num
  sz=0.0_num
  sx0=0.0_num
  sy0=0.0_num
  sz0=0.0_num
  DO ip=1, np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    IF (l_lower_order_in_v) THEN
      IF (nox==2*(nox/2)) THEN
        j=nint(x)
        j0=floor(x-stagger_shift)
      ELSE
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
      END IF
      IF (noy==2*(noy/2)) THEN
        k=nint(y)
        k0=floor(y-stagger_shift)
      ELSE
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
      END IF
      IF (noz==2*(noz/2)) THEN
        l=nint(z)
        l0=floor(z-stagger_shift)
      ELSE
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)
      END IF
    ELSE
      IF (nox==2*(nox/2)) THEN
        j=nint(x)
        j0=floor(x+0.5_num-stagger_shift)
      ELSE
        j=floor(x)
        j0=floor(x-stagger_shift)
      END IF
      IF (noy==2*(noy/2)) THEN
        k=nint(y)
        k0=floor(y+0.5_num-stagger_shift)
      ELSE
        k=floor(y)
        k0=floor(y-stagger_shift)
      END IF
      IF (noz==2*(noz/2)) THEN
        l=nint(z)
        l0=floor(z+0.5_num-stagger_shift)
      ELSE
        l=floor(z)
        l0=floor(z-stagger_shift)
      END IF
    END IF

    xint=x-j
    yint=y-k
    zint=z-l

    IF (nox==1) THEN
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
    ELSEIF (nox==2) THEN
      xintsq = xint*xint
      sx(-1) = 0.5_num*(0.5_num-xint)**2
      sx( 0) = 0.75_num-xintsq
      sx( 1) = 0.5_num*(0.5_num+xint)**2
    ELSEIF (nox==3) THEN
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2) = onesixth*xintsq*xint
    END IF

    IF (noy==1) THEN
      sy( 0) = 1.0_num-yint
      sy( 1) = yint
    ELSEIF (noy==2) THEN
      yintsq = yint*yint
      sy(-1) = 0.5_num*(0.5_num-yint)**2
      sy( 0) = 0.75_num-yintsq
      sy( 1) = 0.5_num*(0.5_num+yint)**2
    ELSEIF (noy==3) THEN
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
      sy( 2) = onesixth*yintsq*yint
    END IF

    IF (noz==1) THEN
      sz( 0) = 1.0_num-zint
      sz( 1) = zint
    ELSEIF (noz==2) THEN
      zintsq = zint*zint
      sz(-1) = 0.5_num*(0.5_num-zint)**2
      sz( 0) = 0.75_num-zintsq
      sz( 1) = 0.5_num*(0.5_num+zint)**2
    ELSEIF (noz==3) THEN
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2) = onesixth*zintsq*zint
    END IF

    xint=x-stagger_shift-j0
    yint=y-stagger_shift-k0
    zint=z-stagger_shift-l0

    IF (l_lower_order_in_v) THEN
      IF (nox==1) THEN
        sx0( 0) = 1.0_num
      ELSEIF (nox==2) THEN
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
      ELSEIF (nox==3) THEN
        xintsq = xint*xint
        sx0(-1) = 0.5_num*(0.5_num-xint)**2
        sx0( 0) = 0.75_num-xintsq
        sx0( 1) = 0.5_num*(0.5_num+xint)**2
      END IF

      IF (noy==1) THEN
        sy0( 0) = 1.0_num
      ELSEIF (noy==2) THEN
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
      ELSEIF (noy==3) THEN
        yintsq = yint*yint
        sy0(-1) = 0.5_num*(0.5_num-yint)**2
        sy0( 0) = 0.75_num-yintsq
        sy0( 1) = 0.5_num*(0.5_num+yint)**2
      END IF

      IF (noz==1) THEN
        sz0( 0) = 1.0_num
      ELSEIF (noz==2) THEN
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint
      ELSEIF (noz==3) THEN
        zintsq = zint*zint
        sz0(-1) = 0.5_num*(0.5_num-zint)**2
        sz0( 0) = 0.75_num-zintsq
        sz0( 1) = 0.5_num*(0.5_num+zint)**2
      END IF
    ELSE

      IF (nox==1) THEN
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
      ELSEIF (nox==2) THEN
        xintsq = xint*xint
        sx0(-1) = 0.5_num*(0.5_num-xint)**2
        sx0( 0) = 0.75_num-xintsq
        sx0( 1) = 0.5_num*(0.5_num+xint)**2
      ELSEIF (nox==3) THEN
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(-1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx0( 2) = onesixth*xintsq*xint
      END IF

      IF (noy==1) THEN
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
      ELSEIF (noy==2) THEN
        yintsq = yint*yint
        sy0(-1) = 0.5_num*(0.5_num-yint)**2
        sy0( 0) = 0.75_num-yintsq
        sy0( 1) = 0.5_num*(0.5_num+yint)**2
      ELSEIF (noy==3) THEN
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
        sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
        sy0( 2) = onesixth*yintsq*yint
      END IF

      IF (noz==1) THEN
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint
      ELSEIF (noz==2) THEN
        zintsq = zint*zint
        sz0(-1) = 0.5_num*(0.5_num-zint)**2
        sz0( 0) = 0.75_num-zintsq
        sz0( 1) = 0.5_num*(0.5_num+zint)**2
      ELSEIF (noz==3) THEN
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz0( 2) = onesixth*zintsq*zint
      END IF
    END IF

    DO ll = izmin0, izmax0
      DO kk = iymin0, iymax0
        DO jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj, k0+kk, l0+ll)*signx
        END DO
      END DO
    END DO

    DO ll = izmin0, izmax0
      DO kk = iymin, iymax+1
        DO jj = ixmin0, ixmax0
          by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj, k+kk, l0+ll)*signy
        END DO
      END DO
    END DO

    DO ll = izmin, izmax+1
      DO kk = iymin0, iymax0
        DO jj = ixmin0, ixmax0
          bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj, k0+kk, l+ll)
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE(sx0, sz0)

  RETURN
END SUBROUTINE pxrgetb3d_n_energy_conserving

! ________________________________________________________________________________________
!> @brief
!> Gathering of magnetic field from Yee grid ("energy conserving") on particles
!> at arbitrary order. WARNING: Highly unoptimized routine
!
!> @details
!> This subroutine is inherited from Warp.
!
!> @Author
!> From Warp
!
!> @date
!> 2015
! ________________________________________________________________________________________
subroutine pxr_getb3d_n_energy_conserving(np, xp, yp, zp, bx, by, bz, xmin, ymin,     &
  zmin, dx, dy, dz, nox, noy, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard,        &
  byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l4symtry, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp)                     :: np, nox, noy, noz
  INTEGER(idp), intent(in)             :: bxg_nguard(3), bxg_nvalid(3),               &
  byg_nguard(3), byg_nvalid(3), bzg_nguard(3), bzg_nvalid(3)
  real(num), dimension(np)         :: xp, yp, zp, bx, by, bz
  LOGICAL(lp)       :: l4symtry, l_lower_order_in_v, l_nodal
  REAL(num)         :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1,                                       &
  -bxg_nguard(3):bxg_nvalid(3)+bxg_nguard(3)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1,                                       &
  -byg_nguard(3):byg_nvalid(3)+byg_nguard(3)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1,                                       &
  -bzg_nguard(3):bzg_nvalid(3)+bzg_nguard(3)-1)
  real(num) :: xmin, ymin, zmin, dx, dy, dz
  integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, ixmin0,      &
  ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
  real(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, xintsq, oxint, yintsq,       &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq, signx, signy
  real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  real(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
  real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  real(num), dimension(:), allocatable :: sx0, sy0, sz0
  real(num), parameter :: onesixth=1./6., twothird=2./3.

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dyi = 1./dy
  dzi = 1./dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)-1
  iymin = -int(noy/2)
  iymax =  int((noy+1)/2)-1
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)-1


  if (l_lower_order_in_v) then
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
  else
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
  end if
  allocate(sx0(ixmin0:ixmax0), sy0(iymin0:iymax0), sz0(izmin0:izmax0))

  signx = 1.
  signy = 1.

  sx=0
  sy=0.
  sz=0.
  sx0=0.
  sy0=0.
  sz0=0.

  do ip=1, np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    if (l4symtry) then
      if (x<0.) then
        x = -x
        signx = -1.
      else
        signx = 1.
      end if
      if (y<0.) then
        y = -y
        signy = -1.
      else
        signy = 1.
      end if
    end if

    if (l_lower_order_in_v) then
      if (nox==2*(nox/2)) then
        j=nint(x)
        j0=floor(x-stagger_shift)
      else
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
      end if
      if (noy==2*(noy/2)) then
        k=nint(y)
        k0=floor(y-stagger_shift)
      else
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
      end if
      if (noz==2*(noz/2)) then
        l=nint(z)
        l0=floor(z-stagger_shift)
      else
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)
      end if
    else
      if (nox==2*(nox/2)) then
        j=nint(x)
        j0=floor(x+0.5_num-stagger_shift)
      else
        j=floor(x)
        j0=floor(x-stagger_shift)
      end if
      if (noy==2*(noy/2)) then
        k=nint(y)
        k0=floor(y+0.5_num-stagger_shift)
      else
        k=floor(y)
        k0=floor(y-stagger_shift)
      end if
      if (noz==2*(noz/2)) then
        l=nint(z)
        l0=floor(z+0.5_num-stagger_shift)
      else
        l=floor(z)
        l0=floor(z-stagger_shift)
      end if
    end if

    xint=x-j
    yint=y-k
    zint=z-l

    if (nox==1) then
      sx( 0) = 1.-xint
      sx( 1) = xint
    elseif (nox==2) then
      xintsq = xint*xint
      sx(-1) = 0.5*(0.5-xint)**2
      sx( 0) = 0.75-xintsq
      sx( 1) = 0.5*(0.5+xint)**2
    elseif (nox==3) then
      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.-xint/2)
      sx( 1) = twothird-oxintsq*(1.-oxint/2)
      sx( 2) = onesixth*xintsq*xint
    end if

    if (noy==1) then
      sy( 0) = 1.-yint
      sy( 1) = yint
    elseif (noy==2) then
      yintsq = yint*yint
      sy(-1) = 0.5*(0.5-yint)**2
      sy( 0) = 0.75-yintsq
      sy( 1) = 0.5*(0.5+yint)**2
    elseif (noy==3) then
      oyint = 1.-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.-yint/2)
      sy( 1) = twothird-oyintsq*(1.-oyint/2)
      sy( 2) = onesixth*yintsq*yint
    end if

    if (noz==1) then
      sz( 0) = 1.-zint
      sz( 1) = zint
    elseif (noz==2) then
      zintsq = zint*zint
      sz(-1) = 0.5*(0.5-zint)**2
      sz( 0) = 0.75-zintsq
      sz( 1) = 0.5*(0.5+zint)**2
    elseif (noz==3) then
      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint
    end if

    xint=x-stagger_shift-j0
    yint=y-stagger_shift-k0
    zint=z-stagger_shift-l0

    if (l_lower_order_in_v) then

      if (nox==1) then
        sx0( 0) = 1.
      elseif (nox==2) then
        sx0( 0) = 1.-xint
        sx0( 1) = xint
      elseif (nox==3) then
        xintsq = xint*xint
        sx0(-1) = 0.5*(0.5-xint)**2
        sx0( 0) = 0.75-xintsq
        sx0( 1) = 0.5*(0.5+xint)**2
      end if

      if (noy==1) then
        sy0( 0) = 1.
      elseif (noy==2) then
        sy0( 0) = 1.-yint
        sy0( 1) = yint
      elseif (noy==3) then
        yintsq = yint*yint
        sy0(-1) = 0.5*(0.5-yint)**2
        sy0( 0) = 0.75-yintsq
        sy0( 1) = 0.5*(0.5+yint)**2
      end if

      if (noz==1) then
        sz0( 0) = 1.
      elseif (noz==2) then
        sz0( 0) = 1.-zint
        sz0( 1) = zint
      elseif (noz==3) then
        zintsq = zint*zint
        sz0(-1) = 0.5*(0.5-zint)**2
        sz0( 0) = 0.75-zintsq
        sz0( 1) = 0.5*(0.5+zint)**2
      end if

    else

      if (nox==1) then
        sx0( 0) = 1.-xint
        sx0( 1) = xint
      elseif (nox==2) then
        xintsq = xint*xint
        sx0(-1) = 0.5*(0.5-xint)**2
        sx0( 0) = 0.75-xintsq
        sx0( 1) = 0.5*(0.5+xint)**2
      elseif (nox==3) then
        oxint = 1.-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(-1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.-xint/2)
        sx0( 1) = twothird-oxintsq*(1.-oxint/2)
        sx0( 2) = onesixth*xintsq*xint
      end if

      if (noy==1) then
        sy0( 0) = 1.-yint
        sy0( 1) = yint
      elseif (noy==2) then
        yintsq = yint*yint
        sy0(-1) = 0.5*(0.5-yint)**2
        sy0( 0) = 0.75-yintsq
        sy0( 1) = 0.5*(0.5+yint)**2
      elseif (noy==3) then
        oyint = 1.-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.-yint/2)
        sy0( 1) = twothird-oyintsq*(1.-oyint/2)
        sy0( 2) = onesixth*yintsq*yint
      end if

      if (noz==1) then
        sz0( 0) = 1.-zint
        sz0( 1) = zint
      elseif (noz==2) then
        zintsq = zint*zint
        sz0(-1) = 0.5*(0.5-zint)**2
        sz0( 0) = 0.75-zintsq
        sz0( 1) = 0.5*(0.5+zint)**2
      elseif (noz==3) then
        ozint = 1.-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.-zint/2)
        sz0( 1) = twothird-ozintsq*(1.-ozint/2)
        sz0( 2) = onesixth*zintsq*zint
      end if

    end if

    do ll = izmin0, izmax0
      do kk = iymin0, iymax0
        do jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj, k0+kk, l0+ll)*signx
        end do
      end do
    end do

    do ll = izmin0, izmax0
      do kk = iymin, iymax+1
        do jj = ixmin0, ixmax0
          by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj, k+kk, l0+ll)*signy
        end do
      end do
    end do

    do ll = izmin, izmax+1
      do kk = iymin0, iymax0
        do jj = ixmin0, ixmax0
          bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj, k0+kk, l+ll)
        end do
      end do
    end do

  end do
  deallocate(sx0, sz0)

  return
end subroutine pxr_getb3d_n_energy_conserving

! ________________________________________________________________________________________
!> Gathering of electric field from Yee grid ("energy conserving") on particles
!> at arbitrary order. WARNING: Highly unoptimized routine
!> @brief
!
!> This subroutine is inherited from Warp
!> @details
!
!> @Author
!> From Warp
!
!> @date
!> 2015
! ________________________________________________________________________________________
subroutine pxr_gete3d_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, ymin,     &
  zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,        &
  eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, l4symtry, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp) :: np, nox, noy, noz
  INTEGER(idp), intent(in) :: exg_nguard(3), exg_nvalid(3), eyg_nguard(3),            &
  eyg_nvalid(3), ezg_nguard(3), ezg_nvalid(3)
  real(num), dimension(np) :: xp, yp, zp, ex, ey, ez
  LOGICAL(lp)       :: l4symtry, l_lower_order_in_v, l_nodal
  REAL(num)         :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1,                                       &
  -exg_nguard(3):exg_nvalid(3)+exg_nguard(3)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1,                                       &
  -eyg_nguard(3):eyg_nvalid(3)+eyg_nguard(3)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1,                                       &
  -ezg_nguard(3):ezg_nvalid(3)+ezg_nguard(3)-1)
  real(num) :: xmin, ymin, zmin, dx, dy, dz
  integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, ixmin0,      &
  ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
  real(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, xintsq, oxint, yintsq,       &
  oyint, zintsq, ozint, oxintsq, oyintsq, ozintsq, signx, signy
  real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  real(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
  real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  real(num), dimension(:), allocatable :: sx0, sy0, sz0
  real(num), parameter :: onesixth=1./6., twothird=2./3.

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dyi = 1./dy
  dzi = 1./dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)-1
  iymin = -int(noy/2)
  iymax =  int((noy+1)/2)-1
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)-1

  if (l_lower_order_in_v) then
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
  else
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
  end if

  allocate(sx0(ixmin0:ixmax0), sy0(iymin0:iymax0), sz0(izmin0:izmax0))

  signx = 1.
  signy = 1.

  do ip=1, np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    if (l4symtry) then
      if (x<0.) then
        x = -x
        signx = -1.
      else
        signx = 1.
      end if
      if (y<0.) then
        y = -y
        signy = -1.
      else
        signy = 1.
      end if
    end if

    if (l_lower_order_in_v) then
      if (nox==2*(nox/2)) then
        j=nint(x)
        j0=floor(x-stagger_shift)
      else
        j=floor(x)
        j0=floor(x+0.5_num-stagger_shift)
      end if
      if (noy==2*(noy/2)) then
        k=nint(y)
        k0=floor(y-stagger_shift)
      else
        k=floor(y)
        k0=floor(y+0.5_num-stagger_shift)
      end if
      if (noz==2*(noz/2)) then
        l=nint(z)
        l0=floor(z-stagger_shift)
      else
        l=floor(z)
        l0=floor(z+0.5_num-stagger_shift)
      end if
    else
      if (nox==2*(nox/2)) then
        j=nint(x)
        j0=floor(x+0.5_num-stagger_shift)
      else
        j=floor(x)
        j0=floor(x-stagger_shift)
      end if
      if (noy==2*(noy/2)) then
        k=nint(y)
        k0=floor(y+0.5_num-stagger_shift)
      else
        k=floor(y)
        k0=floor(y-stagger_shift)
      end if
      if (noz==2*(noz/2)) then
        l=nint(z)
        l0=floor(z+0.5_num-stagger_shift)
      else
        l=floor(z)
        l0=floor(z-stagger_shift)
      end if
    end if

    xint=x-j
    yint=y-k
    zint=z-l

    if (nox==1) then
      sx( 0) = 1.-xint
      sx( 1) = xint
    elseif (nox==2) then
      xintsq = xint*xint
      sx(-1) = 0.5*(0.5-xint)**2
      sx( 0) = 0.75-xintsq
      sx( 1) = 0.5*(0.5+xint)**2
    elseif (nox==3) then
      oxint = 1.-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.-xint/2)
      sx( 1) = twothird-oxintsq*(1.-oxint/2)
      sx( 2) = onesixth*xintsq*xint
    end if

    if (noy==1) then
      sy( 0) = 1.-yint
      sy( 1) = yint
    elseif (noy==2) then
      yintsq = yint*yint
      sy(-1) = 0.5*(0.5-yint)**2
      sy( 0) = 0.75-yintsq
      sy( 1) = 0.5*(0.5+yint)**2
    elseif (noy==3) then
      oyint = 1.-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.-yint/2)
      sy( 1) = twothird-oyintsq*(1.-oyint/2)
      sy( 2) = onesixth*yintsq*yint
    end if

    if (noz==1) then
      sz( 0) = 1.-zint
      sz( 1) = zint
    elseif (noz==2) then
      zintsq = zint*zint
      sz(-1) = 0.5*(0.5-zint)**2
      sz( 0) = 0.75-zintsq
      sz( 1) = 0.5*(0.5+zint)**2
    elseif (noz==3) then
      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint
    end if

    xint=x-stagger_shift-j0
    yint=y-stagger_shift-k0
    zint=z-stagger_shift-l0

    if (l_lower_order_in_v) then

      if (nox==1) then
        sx0( 0) = 1.
      elseif (nox==2) then
        sx0( 0) = 1.-xint
        sx0( 1) = xint
      elseif (nox==3) then
        xintsq = xint*xint
        sx0(-1) = 0.5*(0.5-xint)**2
        sx0( 0) = 0.75-xintsq
        sx0( 1) = 0.5*(0.5+xint)**2
      end if

      if (noy==1) then
        sy0( 0) = 1.
      elseif (noy==2) then
        sy0( 0) = 1.-yint
        sy0( 1) = yint
      elseif (noy==3) then
        yintsq = yint*yint
        sy0(-1) = 0.5*(0.5-yint)**2
        sy0( 0) = 0.75-yintsq
        sy0( 1) = 0.5*(0.5+yint)**2
      end if

      if (noz==1) then
        sz0( 0) = 1.
      elseif (noz==2) then
        sz0( 0) = 1.-zint
        sz0( 1) = zint
      elseif (noz==3) then
        zintsq = zint*zint
        sz0(-1) = 0.5*(0.5-zint)**2
        sz0( 0) = 0.75-zintsq
        sz0( 1) = 0.5*(0.5+zint)**2
      end if

    else

      if (nox==1) then
        sx0( 0) = 1.-xint
        sx0( 1) = xint
      elseif (nox==2) then
        xintsq = xint*xint
        sx0(-1) = 0.5*(0.5-xint)**2
        sx0( 0) = 0.75-xintsq
        sx0( 1) = 0.5*(0.5+xint)**2
      elseif (nox==3) then
        oxint = 1.-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(-1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.-xint/2)
        sx0( 1) = twothird-oxintsq*(1.-oxint/2)
        sx0( 2) = onesixth*xintsq*xint
      end if

      if (noy==1) then
        sy0( 0) = 1.-yint
        sy0( 1) = yint
      elseif (noy==2) then
        yintsq = yint*yint
        sy0(-1) = 0.5*(0.5-yint)**2
        sy0( 0) = 0.75-yintsq
        sy0( 1) = 0.5*(0.5+yint)**2
      elseif (noy==3) then
        oyint = 1.-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.-yint/2)
        sy0( 1) = twothird-oyintsq*(1.-oyint/2)
        sy0( 2) = onesixth*yintsq*yint
      end if


      if (noz==1) then
        sz0( 0) = 1.-zint
        sz0( 1) = zint
      elseif (noz==2) then
        zintsq = zint*zint
        sz0(-1) = 0.5*(0.5-zint)**2
        sz0( 0) = 0.75-zintsq
        sz0( 1) = 0.5*(0.5+zint)**2
      elseif (noz==3) then
        ozint = 1.-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.-zint/2)
        sz0( 1) = twothird-ozintsq*(1.-ozint/2)
        sz0( 2) = onesixth*zintsq*zint
      end if

    end if

    do ll = izmin, izmax+1
      do kk = iymin, iymax+1
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj, k+kk, l+ll)*signx
        end do
      end do
    end do

    do ll = izmin, izmax+1
      do kk = iymin0, iymax0
        do jj = ixmin, ixmax+1
          ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj, k0+kk, l+ll)*signy
        end do
      end do
    end do

    do ll = izmin0, izmax0
      do kk = iymin, iymax+1
        do jj = ixmin, ixmax+1
          ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj, k+kk, l0+ll)
        end do
      end do
    end do

  end do
  deallocate(sx0, sy0, sz0)

  return
end subroutine pxr_gete3d_n_energy_conserving

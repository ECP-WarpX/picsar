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
!
! CHARGE_DEPOSITION_3D.F90
!
! Developers:
! Henri Vincenti, ! Mathieu Lobet
!
! Brief description:
! File containing subroutines for the charge deposition itself in 3D.
!
! List of suboutines:
!
! Scalar subroutines:
! - depose_rho_scalar_1_1_1
! - depose_rho_scalar_2_2_2
! - depose_rho_scalar_3_3_3
!
! Vectorized subroutines:
! - depose_rho_vecSH_1_1_1
! - depose_rho_vecNOY_1_1_1
! - depose_rho_vecHV_1_1_1
! - depose_rho_vecHVv2_1_1_1
! - depose_rho_vecHVv2_2_2_2
! - depose_rho_vecHVv3_3_3_3
! - depose_rho_vecHVv4_3_3_3
!
! General order:
! - pxr_depose_rho_n
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Order 1 3D scalar charge deposition routine
!
!> @details
!> This subroutine computes the charge density sequentially on grid at order 1.
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] lvect: vector length (useless here, just for interface compatibility)
! ________________________________________________________________________________________
SUBROUTINE depose_rho_scalar_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,   &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN OUT) :: rho
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dx, dy, dz, xmin, ymin, zmin
  INTEGER(idp), INTENT (IN) :: lvect
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint
  REAL(num) :: x, y, z, wq, invvol
  REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, ip

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi

  ! Prevent the compiler to vectorize (dependencies)
  !DIR$ NOVECTOR
  !$acc parallel deviceptr(rho, xp, yp, zp, w)
  !$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1))
  DO ip=1, np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- finds node of cell containing particles for current positions
    j=floor(x)
    k=floor(y)
    l=floor(z)
    ! --- computes distance between particle and node for current positions
    xint = x-j
    yint = y-k
    zint = z-l
    ! --- computes particles weights
    wq=q*w(ip)*invvol
    ! --- computes coefficients for node centered quantities
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint
    ! --- add charge density contributions
    !$acc atomic update
    rho(j, k, l)      = rho(j, k, l)+sx(0)*sy(0)*sz(0)*wq
    !$acc atomic update
    rho(j+1, k, l)    = rho(j+1, k, l)+sx(1)*sy(0)*sz(0)*wq
    !$acc atomic update
    rho(j, k+1, l)    = rho(j, k+1, l)+sx(0)*sy(1)*sz(0)*wq
    !$acc atomic update
    rho(j+1, k+1, l)  = rho(j+1, k+1, l)+sx(1)*sy(1)*sz(0)*wq
    !$acc atomic update
    rho(j, k, l+1)    = rho(j, k, l+1)+sx(0)*sy(0)*sz(1)*wq
    !$acc atomic update
    rho(j+1, k, l+1)  = rho(j+1, k, l+1)+sx(1)*sy(0)*sz(1)*wq
    !$acc atomic update
    rho(j, k+1, l+1)  = rho(j, k+1, l+1)+sx(0)*sy(1)*sz(1)*wq
    !$acc atomic update
    rho(j+1, k+1, l+1)= rho(j+1, k+1, l+1)+sx(1)*sy(1)*sz(1)*wq
  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_rho_scalar_1_1_1


! ________________________________________________________________________________________
!> @brief
!> Order 2 3D scalar charge deposition routine
!
!> @details
!> This subroutine computes the charge density sequentially on grid at order 1.
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] lvect: vector length (useless here, just for interface compatibility)
! ________________________________________________________________________________________
SUBROUTINE depose_rho_scalar_2_2_2(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,   &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN OUT) :: rho
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dx, dy, dz, xmin, ymin, zmin
  INTEGER(idp), INTENT (IN) :: lvect
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, xintsq, yintsq, zintsq
  REAL(num) :: x, y, z, wq, invvol, sx1, sx2, sx3, sx4, sx5, sx6, sx7, sx8, sx9
  REAL(num) :: sx(-1:1), sy(-1:1), sz(-1:1)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, ip

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  !DIR$ NOVECTOR
  !$acc parallel deviceptr(rho, xp, yp, zp, w)
  !$acc loop gang vector private(sx(-1:1), sy(-1:1), sz(-1:1))
  DO ip=1, np

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- finds node of cell containing particles for current positions
    j=nint(x)
    k=nint(y)
    l=nint(z)
    ! --- computes distance between particle and node for current positions
    xint = x-j
    yint = y-k
    zint = z-l
    ! --- computes particles weights
    wq=q*w(ip)*invvol
    ! --- computes coefficients for node centered quantities
    xintsq = xint*xint
    sx(-1) = 0.5_num*(0.5_num-xint)**2
    sx( 0) = 0.75_num-xintsq
    sx( 1) = 0.5_num*(0.5_num+xint)**2
    yintsq = yint*yint
    sy(-1) = 0.5_num*(0.5_num-yint)**2
    sy( 0) = 0.75_num-yintsq
    sy( 1) = 0.5_num*(0.5_num+yint)**2
    zintsq = zint*zint
    sz(-1) = 0.5_num*(0.5_num-zint)**2*wq
    sz( 0) = (0.75_num-zintsq)*wq
    sz( 1) = 0.5_num*(0.5_num+zint)**2*wq
    sx1=sx(-1)*sy(-1)
    sx2=sx(0)*sy(-1)
    sx3=sx(1)*sy(-1)
    sx4=sx(-1)*sy(0)
    sx5=sx(0)*sy(0)
    sx6=sx(1)*sy(0)
    sx7=sx(-1)*sy(1)
    sx8=sx(0)*sy(1)
    sx9=sx(1)*sy(1)
    ! --- add charge density contributions to the 27
    ! --- nearest vertices
    !$acc atomic update
    rho(j-1, k-1, l-1) = rho(j-1, k-1, l-1)   + sx1*sz(-1)
    !$acc atomic update
    rho(j, k-1, l-1)   = rho(j, k-1, l-1)     + sx2*sz(-1)
    !$acc atomic update
    rho(j+1, k-1, l-1) = rho(j+1, k-1, l-1)   + sx3*sz(-1)
    !$acc atomic update
    rho(j-1, k, l-1)   = rho(j-1, k, l-1)     + sx4*sz(-1)
    !$acc atomic update
    rho(j, k, l-1)     = rho(j, k, l-1)       + sx5*sz(-1)
    !$acc atomic update
    rho(j+1, k, l-1)   = rho(j+1, k, l-1)     + sx6*sz(-1)
    !$acc atomic update
    rho(j-1, k+1, l-1) = rho(j-1, k+1, l-1)   + sx7*sz(-1)
    !$acc atomic update
    rho(j, k+1, l-1)   = rho(j, k+1, l-1)     + sx8*sz(-1)
    !$acc atomic update
    rho(j+1, k+1, l-1) = rho(j+1, k+1, l-1)   + sx9*sz(-1)
    !$acc atomic update
    rho(j-1, k-1, l)   = rho(j-1, k-1, l)     + sx1*sz(0)
    !$acc atomic update
    rho(j, k-1, l)     = rho(j, k-1, l)       + sx2*sz(0)
    !$acc atomic update
    rho(j+1, k-1, l)   = rho(j+1, k-1, l)     + sx3*sz(0)
    !$acc atomic update
    rho(j-1, k, l)     = rho(j-1, k, l)       + sx4*sz(0)
    !$acc atomic update
    rho(j, k, l)       = rho(j, k, l)         + sx5*sz(0)
    !$acc atomic update
    rho(j+1, k, l)     = rho(j+1, k, l)       + sx6*sz(0)
    !$acc atomic update
    rho(j-1, k+1, l)   = rho(j-1, k+1, l)     + sx7*sz(0)
    !$acc atomic update
    rho(j, k+1, l)     = rho(j, k+1, l)       + sx8*sz(0)
    !$acc atomic update
    rho(j+1, k+1, l)   = rho(j+1, k+1, l)     + sx9*sz(0)
    !$acc atomic update
    rho(j-1, k-1, l+1) = rho(j-1, k-1, l+1)   + sx1*sz(1)
    !$acc atomic update
    rho(j, k-1, l+1)   = rho(j, k-1, l+1)     + sx2*sz(1)
    !$acc atomic update
    rho(j+1, k-1, l+1) = rho(j+1, k-1, l+1)   + sx3*sz(1)
    !$acc atomic update
    rho(j-1, k, l+1)   = rho(j-1, k, l+1)     + sx4*sz(1)
    !$acc atomic update
    rho(j, k, l+1)     = rho(j, k, l+1)       + sx5*sz(1)
    !$acc atomic update
    rho(j+1, k, l+1)   = rho(j+1, k, l+1)     + sx6*sz(1)
    !$acc atomic update
    rho(j-1, k+1, l+1) = rho(j-1, k+1, l+1)   + sx7*sz(1)
    !$acc atomic update
    rho(j, k+1, l+1)   = rho(j, k+1, l+1)     + sx8*sz(1)
    !$acc atomic update
    rho(j+1, k+1, l+1) = rho(j+1, k+1, l+1)   + sx9*sz(1)
  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_rho_scalar_2_2_2

! ________________________________________________________________________________________
!> @brief
!> Order 3 3D scalar charge deposition routine
!
!> @details
!> This subroutine computes the charge density sequentially on grid at order 3.
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] lvect: vector length (useless here, just for interface compatibility)
!>
! ________________________________________________________________________________________
SUBROUTINE depose_rho_scalar_3_3_3(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,   &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN OUT) :: rho
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dx, dy, dz, xmin, ymin, zmin
  INTEGER(idp), INTENT (IN) :: lvect
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq,  &
  zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, wq, invvol
  REAL(num) :: sx(-1:2), sy(-1:2), sz(-1:2)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, ip

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  !DIR$ NOVECTOR
  !$acc parallel deviceptr(rho, xp, yp, zp, w)
  !$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2))
  DO ip=1, np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- finds node of cell containing particles for current positions
    j=floor(x)
    k=floor(y)
    l=floor(z)
    ! --- computes distance between particle and node for current positions
    xint = x-j
    yint = y-k
    zint = z-l
    ! --- computes particles weights
    wq=q*w(ip)*invvol
    ! --- computes coefficients for node centered quantities
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1) = onesixth*oxintsq*oxint
    sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
    sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
    sx( 2) = onesixth*xintsq*xint
    oyint = 1.0_num-yint
    yintsq = yint*yint
    oyintsq = oyint*oyint
    sy(-1) = onesixth*oyintsq*oyint
    sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
    sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
    sy( 2) = onesixth*yintsq*yint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1) = onesixth*ozintsq*ozint
    sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
    sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
    sz( 2) = onesixth*zintsq*zint
    ! --- add charge density contributions to the 64
    ! --- nearest vertices
    ! --- add charge density contributions
    !$acc atomic update
    rho(j-1, k-1, l-1)=rho(j-1, k-1, l-1)+sx(-1)*sy(-1)*sz(-1)*wq
    !$acc atomic update
    rho(j, k-1, l-1)=rho(j, k-1, l-1)+sx(0)*sy(-1)*sz(-1)*wq
    !$acc atomic update
    rho(j+1, k-1, l-1)=rho(j+1, k-1, l-1)+sx(1)*sy(-1)*sz(-1)*wq
    !$acc atomic update
    rho(j+2, k-1, l-1)=rho(j+2, k-1, l-1)+sx(2)*sy(-1)*sz(-1)*wq
    !$acc atomic update
    rho(j-1, k, l-1)=rho(j-1, k, l-1)+sx(-1)*sy(0)*sz(-1)*wq
    !$acc atomic update
    rho(j, k, l-1)=rho(j, k, l-1)+sx(0)*sy(0)*sz(-1)*wq
    !$acc atomic update
    rho(j+1, k, l-1)=rho(j+1, k, l-1)+sx(1)*sy(0)*sz(-1)*wq
    !$acc atomic update
    rho(j+2, k, l-1)=rho(j+2, k, l-1)+sx(2)*sy(0)*sz(-1)*wq
    !$acc atomic update
    rho(j-1, k+1, l-1)=rho(j-1, k+1, l-1)+sx(-1)*sy(1)*sz(-1)*wq
    !$acc atomic update
    rho(j, k+1, l-1)=rho(j, k+1, l-1)+sx(0)*sy(1)*sz(-1)*wq
    !$acc atomic update
    rho(j+1, k+1, l-1)=rho(j+1, k+1, l-1)+sx(1)*sy(1)*sz(-1)*wq
    !$acc atomic update
    rho(j+2, k+1, l-1)=rho(j+2, k+1, l-1)+sx(2)*sy(1)*sz(-1)*wq
    !$acc atomic update
    rho(j-1, k+2, l-1)=rho(j-1, k+2, l-1)+sx(-1)*sy(2)*sz(-1)*wq
    !$acc atomic update
    rho(j, k+2, l-1)=rho(j, k+2, l-1)+sx(0)*sy(2)*sz(-1)*wq
    !$acc atomic update
    rho(j+1, k+2, l-1)=rho(j+1, k+2, l-1)+sx(1)*sy(2)*sz(-1)*wq
    !$acc atomic update
    rho(j+2, k+2, l-1)=rho(j+2, k+2, l-1)+sx(2)*sy(2)*sz(-1)*wq
    !$acc atomic update
    rho(j-1, k-1, l)=rho(j-1, k-1, l)+sx(-1)*sy(-1)*sz(0)*wq
    !$acc atomic update
    rho(j, k-1, l)=rho(j, k-1, l)+sx(0)*sy(-1)*sz(0)*wq
    !$acc atomic update
    rho(j+1, k-1, l)=rho(j+1, k-1, l)+sx(1)*sy(-1)*sz(0)*wq
    !$acc atomic update
    rho(j+2, k-1, l)=rho(j+2, k-1, l)+sx(2)*sy(-1)*sz(0)*wq
    !$acc atomic update
    rho(j-1, k, l)=rho(j-1, k, l)+sx(-1)*sy(0)*sz(0)*wq
    !$acc atomic update
    rho(j, k, l)=rho(j, k, l)+sx(0)*sy(0)*sz(0)*wq
    !$acc atomic update
    rho(j+1, k, l)=rho(j+1, k, l)+sx(1)*sy(0)*sz(0)*wq
    !$acc atomic update
    rho(j+2, k, l)=rho(j+2, k, l)+sx(2)*sy(0)*sz(0)*wq
    !$acc atomic update
    rho(j-1, k+1, l)=rho(j-1, k+1, l)+sx(-1)*sy(1)*sz(0)*wq
    !$acc atomic update
    rho(j, k+1, l)=rho(j, k+1, l)+sx(0)*sy(1)*sz(0)*wq
    !$acc atomic update
    rho(j+1, k+1, l)=rho(j+1, k+1, l)+sx(1)*sy(1)*sz(0)*wq
    !$acc atomic update
    rho(j+2, k+1, l)=rho(j+2, k+1, l)+sx(2)*sy(1)*sz(0)*wq
    !$acc atomic update
    rho(j-1, k+2, l)=rho(j-1, k+2, l)+sx(-1)*sy(2)*sz(0)*wq
    !$acc atomic update
    rho(j, k+2, l)=rho(j, k+2, l)+sx(0)*sy(2)*sz(0)*wq
    !$acc atomic update
    rho(j+1, k+2, l)=rho(j+1, k+2, l)+sx(1)*sy(2)*sz(0)*wq
    !$acc atomic update
    rho(j+2, k+2, l)=rho(j+2, k+2, l)+sx(2)*sy(2)*sz(0)*wq
    !$acc atomic update
    rho(j-1, k-1, l+1)=rho(j-1, k-1, l+1)+sx(-1)*sy(-1)*sz(1)*wq
    !$acc atomic update
    rho(j, k-1, l+1)=rho(j, k-1, l+1)+sx(0)*sy(-1)*sz(1)*wq
    !$acc atomic update
    rho(j+1, k-1, l+1)=rho(j+1, k-1, l+1)+sx(1)*sy(-1)*sz(1)*wq
    !$acc atomic update
    rho(j+2, k-1, l+1)=rho(j+2, k-1, l+1)+sx(2)*sy(-1)*sz(1)*wq
    !$acc atomic update
    rho(j-1, k, l+1)=rho(j-1, k, l+1)+sx(-1)*sy(0)*sz(1)*wq
    !$acc atomic update
    rho(j, k, l+1)=rho(j, k, l+1)+sx(0)*sy(0)*sz(1)*wq
    !$acc atomic update
    rho(j+1, k, l+1)=rho(j+1, k, l+1)+sx(1)*sy(0)*sz(1)*wq
    !$acc atomic update
    rho(j+2, k, l+1)=rho(j+2, k, l+1)+sx(2)*sy(0)*sz(1)*wq
    !$acc atomic update
    rho(j-1, k+1, l+1)=rho(j-1, k+1, l+1)+sx(-1)*sy(1)*sz(1)*wq
    !$acc atomic update
    rho(j, k+1, l+1)=rho(j, k+1, l+1)+sx(0)*sy(1)*sz(1)*wq
    !$acc atomic update
    rho(j+1, k+1, l+1)=rho(j+1, k+1, l+1)+sx(1)*sy(1)*sz(1)*wq
    !$acc atomic update
    rho(j+2, k+1, l+1)=rho(j+2, k+1, l+1)+sx(2)*sy(1)*sz(1)*wq
    !$acc atomic update
    rho(j-1, k+2, l+1)=rho(j-1, k+2, l+1)+sx(-1)*sy(2)*sz(1)*wq
    !$acc atomic update
    rho(j, k+2, l+1)=rho(j, k+2, l+1)+sx(0)*sy(2)*sz(1)*wq
    !$acc atomic update
    rho(j+1, k+2, l+1)=rho(j+1, k+2, l+1)+sx(1)*sy(2)*sz(1)*wq
    !$acc atomic update
    rho(j+2, k+2, l+1)=rho(j+2, k+2, l+1)+sx(2)*sy(2)*sz(1)*wq
    !$acc atomic update
    rho(j-1, k-1, l+2)=rho(j-1, k-1, l+2)+sx(-1)*sy(-1)*sz(2)*wq
    !$acc atomic update
    rho(j, k-1, l+2)=rho(j, k-1, l+2)+sx(0)*sy(-1)*sz(2)*wq
    !$acc atomic update
    rho(j+1, k-1, l+2)=rho(j+1, k-1, l+2)+sx(1)*sy(-1)*sz(2)*wq
    !$acc atomic update
    rho(j+2, k-1, l+2)=rho(j+2, k-1, l+2)+sx(2)*sy(-1)*sz(2)*wq
    !$acc atomic update
    rho(j-1, k, l+2)=rho(j-1, k, l+2)+sx(-1)*sy(0)*sz(2)*wq
    !$acc atomic update
    rho(j, k, l+2)=rho(j, k, l+2)+sx(0)*sy(0)*sz(2)*wq
    !$acc atomic update
    rho(j+1, k, l+2)=rho(j+1, k, l+2)+sx(1)*sy(0)*sz(2)*wq
    !$acc atomic update
    rho(j+2, k, l+2)=rho(j+2, k, l+2)+sx(2)*sy(0)*sz(2)*wq
    !$acc atomic update
    rho(j-1, k+1, l+2)=rho(j-1, k+1, l+2)+sx(-1)*sy(1)*sz(2)*wq
    !$acc atomic update
    rho(j, k+1, l+2)=rho(j, k+1, l+2)+sx(0)*sy(1)*sz(2)*wq
    !$acc atomic update
    rho(j+1, k+1, l+2)=rho(j+1, k+1, l+2)+sx(1)*sy(1)*sz(2)*wq
    !$acc atomic update
    rho(j+2, k+1, l+2)=rho(j+2, k+1, l+2)+sx(2)*sy(1)*sz(2)*wq
    !$acc atomic update
    rho(j-1, k+2, l+2)=rho(j-1, k+2, l+2)+sx(-1)*sy(2)*sz(2)*wq
    !$acc atomic update
    rho(j, k+2, l+2)=rho(j, k+2, l+2)+sx(0)*sy(2)*sz(2)*wq
    !$acc atomic update
    rho(j+1, k+2, l+2)=rho(j+1, k+2, l+2)+sx(1)*sy(2)*sz(2)*wq
    !$acc atomic update
    rho(j+2, k+2, l+2)=rho(j+2, k+2, l+2)+sx(2)*sy(2)*sz(2)*wq
  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_rho_scalar_3_3_3


#if defined (DEV)
! ________________________________________________________________________________________
!> @brief
!> Order 1 3D vector charge deposition routine
!
!> @details
!> Computes charge density on grid vectorized with Schwarzmeier and Hewitt scheme
!> This routine does vectorize on SIMD architecture but poor performances
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecSH_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,    &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT) ::                                                        &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq,  &
  zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, wq, invvol, sx0, sy0, sz0, sx1, sy1, sz1
  REAL(num), ALLOCATABLE :: ww(:, :)
  INTEGER(idp), ALLOCATABLE :: ll(:, :)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, nn, ip, n, m, ixmin, ixmax, iymin, iymax, izmin, izmax
  INTEGER(idp) :: nblk
  INTEGER(idp) :: nnx, nnxy, ind0
  INTEGER(idp) :: moff(1:8)

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  ! Vectorization parameters init
  nblk=4! multiple of vector instruction length (32 bytes on Avx)
  ALLOCATE(ww(1:8, 1:nblk), ll(1:8, 1:nblk))

  nnx = nx + 1 + 2*nxguard
  nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
  moff(1) = 0_idp
  moff(2) = 1_idp
  moff(3) = nnx
  moff(4) = nnx+1_idp
  moff(5) = nnxy
  moff(6) = nnxy+1_idp
  moff(7) = nnxy+nnx
  moff(8) = nnxy+nnx+1_idp

  DO ip=1, np, nblk

    ! Alignement
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED w:64, ww:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, ww, xp, yp, zp, w)
#endif
#endif

    ! Vectorization
#ifndef NOVEC
#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
#endif
    DO n=ip, MIN(ip+nblk-1, np)!!!! vector
      nn=n-ip+1
      !- Computations relative to particle n
      ! --- computes current position in grid units
      x = (xp(n)-xmin)*dxi
      y = (yp(n)-ymin)*dyi
      z = (zp(n)-zmin)*dzi
      ! --- finds node of cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint = y-k
      zint = z-l
      ! --- computes particles weights
      wq=q*w(n)*invvol
      ! --- computes coefficients for node centered quantities
      sx0 = 1.0_num-xint
      sx1 = xint
      sy0 = 1.0_num-yint
      sy1 = yint
      sz0 = 1.0_num-zint
      sz1 = zint
      ! --- computes weight for each of the 8-vertices surrounding particle n
      ind0 = (j+nxguard+1) + (k+nyguard+1)*nnx + (l+nzguard+1)*nnxy
      ww(1, nn) = sx0*sy0*sz0*wq
      ll(1, nn) = ind0+moff(1)
      ww(2, nn) = sx1*sy0*sz0*wq
      ll(2, nn) = ind0+moff(2)
      ww(3, nn) = sx0*sy1*sz0*wq
      ll(3, nn) = ind0+moff(3)
      ww(4, nn) = sx1*sy1*sz0*wq
      ll(4, nn) = ind0+moff(4)
      ww(5, nn) = sx0*sy0*sz1*wq
      ll(5, nn) = ind0+moff(5)
      ww(6, nn) = sx1*sy0*sz1*wq
      ll(6, nn) = ind0+moff(6)
      ww(7, nn) = sx0*sy1*sz1*wq
      ll(7, nn) = ind0+moff(7)
      ww(8, nn) = sx1*sy1*sz1*wq
      ll(8, nn) = ind0+moff(8)
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! --- add charge density contributions
    DO m= 1, MIN(nblk, np-ip+1)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED ww:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, ww)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO l=1, 8!!!! Vector
        rho(ll(l, m)) = rho(ll(l, m))+ww(l, m)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(ww, ll)
  RETURN
END SUBROUTINE depose_rho_vecSH_1_1_1
#endif

! ________________________________________________________________________________________
!> @brief
!> Order 1 3D vector charge deposition routine
!
!> @details
!> Computes charge density on grid vectorized with Nishiguchi, Orii and Yabe scheme (NOY)
!> This routine does vectorize on SIMD architecture but poor performances.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecNOY_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,   &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN OUT) :: rho
  REAL(num), DIMENSION(:, :, :, :), ALLOCATABLE :: rho1
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint
  REAL(num) :: x, y, z, wq, invvol
  REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, vv, n, ip, jj
  INTEGER(idp), INTENT(IN) :: lvect
  REAL(num), DIMENSION(lvect, 8) :: ww

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  ALLOCATE(rho1(1:lvect, -nxguard:nx+nxguard, -nyguard:ny+nyguard,                    &
  -nzguard:nz+nzguard))
  rho1=0.0_num
  DO ip=1, np, lvect
    !DIR$ ASSUME_ALIGNED xp:32
    !DIR$ ASSUME_ALIGNED yp:32
    !DIR$ ASSUME_ALIGNED zp:32
    !DIR$ ASSUME_ALIGNED w:32
    !DIR$ ASSUME_ALIGNED ww:32
    !DIR$ IVDEP
    DO vv=1, MIN(lvect, np-ip+1)!!! Vector
      n=vv+ip-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x = (xp(n)-xmin)*dxi
      y = (yp(n)-ymin)*dyi
      z = (zp(n)-zmin)*dzi
      ! --- finds node of cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint = y-k
      zint = z-l
      ! --- computes particles weights
      wq=q*w(n)*invvol
      ! --- computes coefficients for node centered quantities
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sy( 0) = 1.0_num-yint
      sy( 1) = yint
      sz( 0) = (1.0_num-zint)*wq
      sz( 1) = zint*wq
      ww(vv, 1) = sx(0)*sy(0)*sz(0)
      ww(vv, 2) = sx(1)*sy(0)*sz(0)
      ww(vv, 3) = sx(0)*sy(1)*sz(0)
      ww(vv, 4) = sx(1)*sy(1)*sz(0)
      ww(vv, 5) = sx(0)*sy(0)*sz(1)
      ww(vv, 6) = sx(1)*sy(0)*sz(1)
      ww(vv, 7) = sx(0)*sy(1)*sz(1)
      ww(vv, 8) = sx(1)*sy(1)*sz(1)
    END DO
  END DO

  DO jj=-nxguard, nxguard+nx!!! Vector
    DO vv=1, lvect
      rho(jj, :, :)=rho(jj, :, :)+rho1(vv, jj, :, :)
    END DO
  END DO
  DEALLOCATE(rho1)
  RETURN
END SUBROUTINE depose_rho_vecNOY_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 1 3D vector charge deposition routine
!>
!> @details
!> Computes charge density on grid vectorized (HV-SCHEME v2)
!> This routine does vectorize on SIMD architecture but poor performances
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecHV_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,    &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT) ::                                                        &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: rhocells
  INTEGER(idp), PARAMETER :: LVEC2=8
  INTEGER(idp), DIMENSION(LVEC2) :: ICELL
  REAL(num) :: wq
  INTEGER(idp) :: NCELLS
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi
  REAL(num) :: xint, yint, zint
  REAL(num) :: x, y, z, invvol
  REAL(num) :: sx(0:1), sy(0:1), sz(0:1)
  REAL(num) :: ww(1:LVEC2, 8)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: ic, j, k, l, n, ip, nv, nn
  INTEGER(idp) :: nnx, nnxy
  INTEGER(idp) :: moff(1:8)

  ! Init parameters
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  NCELLS=(2*nxguard+nx)*(2*nyguard+ny)*(2*nzguard+nz)
  ALLOCATE(rhocells(8, NCELLS))
  rhocells=0.0_num
  nnx = nx + 1 + 2*nxguard
  nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
  moff(1) = 0_idp
  moff(2) = 1_idp
  moff(3) = nnx
  moff(4) = nnx+1_idp
  moff(5) = nnxy
  moff(6) = nnxy+1_idp
  moff(7) = nnxy+nnx
  moff(8) = nnxy+nnx+1_idp

  ! FIRST LOOP: computes cell index of particle and their weight on vertices
  DO ip=1, np, LVEC2
    !DIR$ ASSUME_ALIGNED xp:64
    !DIR$ ASSUME_ALIGNED yp:64
    !DIR$ ASSUME_ALIGNED zp:64
    !DIR$ ASSUME_ALIGNED w:64
    !DIR$ ASSUME_ALIGNED ww:64
    !DIR$ ASSUME_ALIGNED ICELL:64
    !DIR$ IVDEP
    DO n=1, MIN(LVEC2, np-ip+1)
      nn=ip+n-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x= (xp(nn)-xmin)*dxi
      z = (zp(nn)-zmin)*dzi
      y = (yp(nn)-ymin)*dyi
      ! --- finds cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ICELL(n)=1+j+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint = y-k
      zint = z-l
      ! --- computes particles weights
      wq=q*w(nn)*invvol
      ! --- Compute weight
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
      sy( 0) = 1.0_num-yint
      sy( 1) = yint
      sz( 0) = (1.0_num-zint)*wq
      sz( 1) = zint*wq
      ww(n, 1) = sx(0)*sy(0)*sz(0)
      ww(n, 2) = sx(1)*sy(0)*sz(0)
      ww(n, 3) = sx(0)*sy(1)*sz(0)
      ww(n, 4) = sx(1)*sy(1)*sz(0)
      ww(n, 5) = sx(0)*sy(0)*sz(1)
      ww(n, 6) = sx(1)*sy(0)*sz(1)
      ww(n, 7) = sx(0)*sy(1)*sz(1)
      ww(n, 8) = sx(1)*sy(1)*sz(1)
    END DO
    ! Current deposition on vertices
    DO n=1, MIN(LVEC2, np-ip+1)
      ! --- add charge density contributions to vertices of the current cell
      ic=ICELL(n)
      !DIR$ ASSUME_ALIGNED rhocells:64
      !DIR$ ASSUME_ALIGNED ww:64
      !DIR NOUNROLL
      DO nv=1, 8!!! - VECTOR
        rhocells(nv, ic)=rhocells(nv, ic)+ww(n, nv)
      END DO
    END DO
  END DO
  DO nv=1, 8
    !DIR$ ASSUME_ALIGNED rhocells:64
    !DIR$ ASSUME_ALIGNED rho:64
    !DIR$ IVDEP
    DO ic=1, NCELLS!!! VECTOR
      rho(ic+moff(nv))=rho(ic+moff(nv))+rhocells(nv, ic)
    END DO
  END DO
  DEALLOCATE(rhocells)
  RETURN
END SUBROUTINE depose_rho_vecHV_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 1 3D vector charge deposition routine
!
!> @details
!> Computes charge density on grid vectorized at order 1 (HV-SCHEME v2).
!> This routine does vectorize on SIMD architecture with good performances:
!> Speedup>2 on AVX 256 bits.
!> The parameter lvect is the vector length and is originally at 64 for order 1.
!>
!> @image html charge_deposition_grid.jpg "Description of the rho structure
!> (2d vision) and internal parameters"
!>
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] lvect vector length
!
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecHVv2_1_1_1(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,  &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  !bind(C, name="depose_rho_vecHVv2_1_1_1")
  IMPLICIT NONE

  INTEGER(idp), INTENT (IN) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT)  ::                                                       &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  INTEGER(idp), INTENT (IN) :: lvect
  REAL(num), INTENT (IN)    :: xp(np), yp(np), zp(np), w(np)
  REAL(num), INTENT (IN)    :: q, dx, dy, dz, xmin, ymin, zmin

  INTEGER(idp), DIMENSION(lvect) :: ICELL
  REAL(num)                 :: ww
  INTEGER(idp) :: NCELLS
  REAL(num) :: dxi, dyi, dzi
  REAL(num) :: x, y, z, invvol
  REAL(num) :: sx(lvect), sy(lvect), sz(lvect), wq(lvect)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: ic, igrid, j, k, l, n, ip, nv, nn
  INTEGER(idp) :: nnx, nnxy
  INTEGER(idp) :: moff(1:8)
  REAL(num):: mx(1:8), my(1:8), mz(1:8), sgn(1:8)
  INTEGER(idp) :: orig, jorig, korig, lorig
  INTEGER(idp) :: ncx, ncy, ncxy, ncz, ix, iy, iz, ngridx, ngridy, ngx, ngxy
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: rhocells
  !DIR$ ATTRIBUTES ALIGN:64 :: rhocells

  ! Init parameters
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
  ncx=nx+2;ncy=ny+2;ncz=nz+2
  NCELLS=ncx*ncy*ncz
  ALLOCATE(rhocells(8, NCELLS))
  rhocells=0.0_num
  nnx = ngridx
  nnxy = nnx*ngridy
  moff = (/0_idp, 1_idp, nnx, nnx+1_idp, nnxy, nnxy+1_idp, nnxy+nnx, nnxy+nnx+1_idp/)
  mx=(/1_num, 0_num, 1_num, 0_num, 1_num, 0_num, 1_num, 0_num/)
  my=(/1_num, 1_num, 0_num, 0_num, 1_num, 1_num, 0_num, 0_num/)
  mz=(/1_num, 1_num, 1_num, 1_num, 0_num, 0_num, 0_num, 0_num/)
  sgn=(/-1_num, 1_num, 1_num, -1_num, 1_num, -1_num, -1_num, 1_num/)
  jorig=-1; korig=-1;lorig=-1
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! ________________________________________________________________________
  ! FIRST LOOP: computes cell index of particle and their weight on vertices
  DO ip=1, np, lvect
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64
    !DIR$ ASSUME_ALIGNED yp:64
    !DIR$ ASSUME_ALIGNED zp:64
    !DIR$ ASSUME_ALIGNED w:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp, w, ICELL)
#endif
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
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x= (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
      ! --- finds cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
      ! --- computes distance between particle and node for current positions
      sx(n) = x-j
      sy(n) = y-k
      sz(n) = z-l
      ! --- computes particles weights
      wq(n)=q*w(nn)*invvol
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! Current deposition on vertices
    DO n=1, MIN(lvect, np-ip+1)
      ! --- add charge density contributions to vertices of the current cell
      ic=ICELL(n)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED rhocells:64
      !DIR$ ASSUME_ALIGNED sx:64
      !DIR$ ASSUME_ALIGNED sy:64
      !DIR$ ASSUME_ALIGNED sz:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, rhocells, sx, sy, sz)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 8!!! - VECTOR
        ww=(-mx(nv)+sx(n))*(-my(nv)+sy(n))* (-mz(nv)+sz(n))*wq(n)*sgn(nv)
        rhocells(nv, ic)=rhocells(nv, ic)+ww
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! - reduction of rhocells in rho
  DO iz=1, ncz
    DO iy=1, ncy
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix=1, ncx!! VECTOR (take ncx multiple of vector length)
        ic=ix+(iy-1)*ncx+(iz-1)*ncxy
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy+orig
        rho(igrid+moff(1))=rho(igrid+moff(1))+rhocells(1, ic)
        rho(igrid+moff(2))=rho(igrid+moff(2))+rhocells(2, ic)
        rho(igrid+moff(3))=rho(igrid+moff(3))+rhocells(3, ic)
        rho(igrid+moff(4))=rho(igrid+moff(4))+rhocells(4, ic)
        rho(igrid+moff(5))=rho(igrid+moff(5))+rhocells(5, ic)
        rho(igrid+moff(6))=rho(igrid+moff(6))+rhocells(6, ic)
        rho(igrid+moff(7))=rho(igrid+moff(7))+rhocells(7, ic)
        rho(igrid+moff(8))=rho(igrid+moff(8))+rhocells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(rhocells)
  RETURN
END SUBROUTINE depose_rho_vecHVv2_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 2 3D vector charge deposition routine
!
!> @details
!> Computes charge density on grid vectorized at order 2 (HV-SCHEME)
!> This routine does vectorize on SIMD architecture with good performances
!> Speedup>2 on AVX 256 bits
!> lvect, the vector length was originally at 16
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] lvect vector length
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecHVv2_2_2_2(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,  &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE

  INTEGER(idp), INTENT (IN)    :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp), INTENT (IN)    :: lvect
  REAL(num), INTENT(IN OUT)     ::                                                    &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT (IN)       :: xp(np), yp(np), zp(np), w(np)
  REAL(num), INTENT (IN)       :: q, dx, dy, dz, xmin, ymin, zmin

  REAL(num), DIMENSION(:, :), ALLOCATABLE :: rhocells
  INTEGER(idp), DIMENSION(lvect) :: ICELL, IG
  REAL(num)                      :: ww
  INTEGER(idp) :: NCELLS
  REAL(num) :: dxi, dyi, dzi
  REAL(num) :: xint, yint, zint, xintsq, yintsq, zintsq
  REAL(num) :: x, y, z, invvol, wq0, wq, szy, syy0, syy1, syy2, szz0, szz1, szz2
  REAL(num) :: sx0(lvect), sx1(lvect), sx2(lvect)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: ic, igrid, j, k, l, n, ip, nv, nn
  INTEGER(idp) :: nnx, nnxy
  INTEGER(idp) :: moff(1:8)
  REAL(num):: ww0(1:lvect, 1:8), www(1:lvect, 1:8)
  INTEGER(idp) :: orig, jorig, korig, lorig
  INTEGER(idp) :: ncx, ncy, ncxy, ncz, ix, iy, iz, ngridx, ngridy, ngx, ngxy

  ! Init parameters
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  wq0=q*invvol
  ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
  ncx=nx+3;ncy=ny+3;ncz=nz+3
  NCELLS=ncx*ncy*ncz
  ALLOCATE(rhocells(8, NCELLS))
  rhocells=0.0_num
  nnx = nx + 1 + 2*nxguard
  nnxy = nnx*(ny+1+2*nyguard)
  moff = (/-nnx-nnxy, -nnxy, nnx-nnxy, -nnx, nnx, -nnx+nnxy, nnxy, nnx+nnxy/)
  ww0=0.0_num
  jorig=-1; korig=-1;lorig=-1
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! FIRST LOOP: computes cell index of particle and their weight on vertices
  DO ip=1, np, lvect

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED w:64, sx0:64, sx1:64, sx2:64
    !DIR$ ASSUME_ALIGNED ICELL:64, IG:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp, w, sx0, sx1, sx2, ICELL, IG)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x= (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
      ! --- finds cell containing particles for current positions
      j=nint(x)
      k=nint(y)
      l=nint(z)
      ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
      IG(n)=ICELL(n)+(k-korig)*ngx+(l-lorig)*ngxy
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint = y-k
      zint = z-l
      xintsq=xint**2
      yintsq=yint**2
      zintsq=zint**2
      ! --- computes particles weights
      wq=w(nn)*wq0
      sx0(n)=0.5_num*(0.5_num-xint)**2
      sx1(n)=(0.75_num-xintsq)
      sx2(n)=0.5_num*(0.5_num+xint)**2
      syy0=0.5_num*(0.5_num-yint)**2
      syy1=(0.75_num-yintsq)
      syy2=0.5_num*(0.5_num+yint)**2
      szz0=0.5_num*(0.5_num-zint)**2*wq
      szz1=(0.75_num-zintsq)*wq
      szz2=0.5_num*(0.5_num+zint)**2*wq
      www(n, 1) = syy0*szz0
      www(n, 2) = syy1*szz0
      www(n, 3) = syy2*szz0
      www(n, 4) = syy0*szz1
      www(n, 5) = syy2*szz1
      www(n, 6) = syy0*szz2
      www(n, 7) = syy1*szz2
      www(n, 8) = syy2*szz2
      szy=syy1*szz1! central point
      ww0(n, 1)=szy*sx0(n)
      ww0(n, 2)=szy*sx1(n)
      ww0(n, 3)=szy*sx2(n)
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! Current deposition on vertices
    DO n=1, MIN(lvect, np-ip+1)
      ! --- add charge density contributions to vertices of the current cell

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED rhocells:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, rhocells)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 8!!! - VECTOR
        ww=www(n, nv)
        ! Loop on (i=-1, j, k)
        rhocells(nv, ICELL(n)-1)=rhocells(nv, ICELL(n)-1)+ww*sx0(n)
        ! Loop on (i=0, j, k)
        rhocells(nv, ICELL(n))=rhocells(nv, ICELL(n))+ww*sx1(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ICELL(n)+1)=rhocells(nv, ICELL(n)+1)+ww*sx2(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 4
        rho(orig+IG(n)+nv-2)=rho(orig+IG(n)+nv-2)+ww0(n, nv)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  ! - reduction of rhocells in rho
  DO iz=1, ncz
    DO iy=1, ncy
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix=1, ncx!! VECTOR (take ncx multiple of vector length)
        ic=ix+(iy-1)*ncx+(iz-1)*ncxy
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
        rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1, ic)
        rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2, ic)
        rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3, ic)
        rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4, ic)
        rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5, ic)
        rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6, ic)
        rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7, ic)
        rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(rhocells)
  RETURN
END SUBROUTINE depose_rho_vecHVv2_2_2_2

#if defined (DEV)
! ________________________________________________________________________________________
!> @brief
!> Order 3 3D vector charge deposition routine version 2
!
!> @details
!> Computes charge density on grid vectorized (HV-SCHEME)
!> This routine does vectorize on SIMD architecture with good performances
!> Speedup>2 on AVX 256 bits
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecHVv2_3_3_3(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,  &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT) ::                                                        &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: rhocells
  INTEGER(idp), PARAMETER :: LVEC2=8
  INTEGER(idp), DIMENSION(LVEC2) :: ICELL
  REAL(num) :: ww, wwx, wwy, wwz
  INTEGER(idp) :: NCELLS
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq,  &
  zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, invvol, wq0, wq
  REAL(num) :: sx1(LVEC2), sx2(LVEC2), sx3(LVEC2), sx4(LVEC2)
  REAL(num) :: sy(-1:2), sz(-1:2)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: ic, j, k, l, vv, n, ip, jj, kk, ll, nv, nn
  INTEGER(idp) :: nnx, nnxy, off0, ind0
  INTEGER(idp) :: moff(1:16)
  REAL(num):: www(1:LVEC2, 1:16)

  ! Init parameters
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  wq0=q*invvol
  NCELLS=(2*nxguard+nx)*(2*nyguard+ny)*(2*nzguard+nz)
  ALLOCATE(rhocells(16, NCELLS))
  rhocells=0.0_num
  nnx = nx + 1 + 2*nxguard
  nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
  moff = (/-1_idp-nnx-nnxy, -nnx-nnxy, 1_idp-nnx-nnxy, -1_idp-nnxy, -nnxy,            &
  1_idp-nnxy, -1_idp+nnx-nnxy, nnx-nnxy, -1_idp-nnx-nnxy, -nnx-nnxy, 1_idp-nnx-nnxy,  &
  -1_idp-nnxy, -nnxy, 1_idp-nnxy, -1_idp+nnx-nnxy, nnx-nnxy/)
  off0=1+nnx+nnxy

  ! FIRST LOOP: computes cell index of particle and their weight on vertices
  DO ip=1, np, LVEC2

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED w:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp, w, ICELL)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
    DO n=1, MIN(LVEC2, np-ip+1)
      nn=ip+n-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x= (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
      ! --- finds cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ICELL(n)=1+j+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)
      wq=w(nn)*wq0
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint = y-k
      zint = z-l
      ! --- computes coefficients for node centered quantities
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx1(n) = onesixth*oxintsq*oxint
      sx2(n) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx3(n) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx4(n) = onesixth*xintsq*xint
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
      sz(-1) = onesixth*ozintsq*ozint*wq
      sz( 0) = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
      sz( 1) = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
      sz( 2) = onesixth*zintsq*zint*wq
      www(n, 1)=sy(-1)*sz(-1)
      www(n, 2)=sy(0)*sz(-1)
      www(n, 3)=sy(1)*sz(-1)
      www(n, 4)=sy(2)*sz(-1)
      www(n, 5)=sy(-1)*sz(0)
      www(n, 6)=sy(0)*sz(0)
      www(n, 7)=sy(1)*sz(0)
      www(n, 8)=sy(2)*sz(0)
      www(n, 9)=sy(-1)*sz(1)
      www(n, 10)=sy(0)*sz(1)
      www(n, 11)=sy(1)*sz(1)
      www(n, 12)=sy(2)*sz(1)
      www(n, 13)=sy(-1)*sz(2)
      www(n, 14)=sy(0)*sz(2)
      www(n, 15)=sy(1)*sz(2)
      www(n, 16)=sy(2)*sz(2)
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! Current deposition on vertices
    DO n=1, MIN(LVEC2, np-ip+1)
      ! --- add charge density contributions to vertices of the current cell
      ic=ICELL(n)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED rhocells:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, rhocells)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 16!!! - VECTOR
        ww=www(n, nv)
        !ww=0.1_num
        ! Loop on (i=-1, j, k)
        rhocells(nv, ic-1)=rhocells(nv, ic-1)+ww*sx1(n)
        ! Loop on (i=0, j, k)
        rhocells(nv, ic)=rhocells(nv, ic)+ww*sx2(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic+1)=rhocells(nv, ic+1)+ww*sx3(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic+2)=rhocells(nv, ic+2)+ww*sx4(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DO nv=1, 16
    ind0=off0+moff(nv)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED rhocells:64
    !DIR$ ASSUME_ALIGNED rho:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, rhocells, rho)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
    DO ic=1, NCELLS!!! VECTOR
      rho(ic+ind0)=rho(ic+ind0)+rhocells(nv, ic)
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
  END DO
  DEALLOCATE(rhocells)
  RETURN
END SUBROUTINE depose_rho_vecHVv2_3_3_3
#endif

#if defined (DEV)
! ________________________________________________________________________________________
!> @brief
!> Order 3 3D vector charge deposition routine
!
!> @details
!> Computes charge density on grid vectorized (HV-SCHEME)
!> This routine does vectorize on SIMD architecture with good performances
!> Speedup>2 on AVX 256 bits
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecHVv3_3_3_3(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,  &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE

  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT) ::                                                        &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: rhocells
  INTEGER(idp), PARAMETER :: LVEC2=16
  INTEGER(idp), DIMENSION(LVEC2) :: ICELL
  REAL(num) :: ww, wwx, wwy, wwz
  INTEGER(idp) :: NCELLS
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint(1:LVEC2), oxint, oyint, ozint, xintsq, &
  yintsq, zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, invvol, wq0, wq
  REAL(num) :: sx1(LVEC2), sx2(LVEC2), sx3(LVEC2), sx4(LVEC2), sy1(LVEC2),            &
  sy2(LVEC2), sy3(LVEC2), sy4(LVEC2), sz1, sz2, sz3, sz4
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: ic, igrid, ic0, j, k, l, vv, n, ip, jj, kk, ll, nv, nn
  INTEGER(idp) :: nnx, nnxy, off0, ind0
  INTEGER(idp) :: moff(1:8)
  REAL(num):: www(1:16, 1:LVEC2), zdec(1:8), h1(1:8), h11(1:8), h12(1:8), sgn(1:8),   &
  szz(1:8)
  INTEGER(idp) :: orig, jorig, korig, lorig
  INTEGER(idp) :: ncx, ncy, ncxy, ncz, ix, iy, iz, ngridx, ngridy, ngx, ngxy

  ! Init parameters
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  wq0=q*invvol
  ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
  ncx=nx+5; ncy=ny+4; ncz=nz+2
  NCELLS=ncx*ncy*ncz
  ALLOCATE(rhocells(8, NCELLS))
  rhocells=0_num; www=0.0_num
  nnx = ngridx
  nnxy = ngridx*ngridy
  moff = (/-nnxy, 0_idp, nnxy, 2_idp*nnxy, nnx-nnxy, nnx, nnx+nnxy, nnx+2_idp*nnxy/)
  jorig=-2_idp; korig=-2_idp;lorig=-1_idp
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  h1(1:4)=(/1_num, 0_num, 1_num, 0_num/); sgn(1:4)=(/1_num, -1_num, 1_num, -1_num/)
  h11(1:4)=(/0_num, 1_num, 1_num, 0_num/); h12(1:4)=(/1_num, 0_num, 0_num, 1_num/)

  ! FIRST LOOP: computes cell index of particle and their weight on vertices
  DO ip=1, np, LVEC2

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED w:64, sx1:64, sx2:64, sx3:64, sx4:64
    !DIR$ ASSUME_ALIGNED w:64, sy1:64, sy2:64, sy3:64, sy4:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp, w, sx1, sx2, sx3, sx4, sy1, sy2, sy3, sy4, ICELL)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
    DO n=1, MIN(LVEC2, np-ip+1)
      nn=ip+n-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x= (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
      ! --- finds cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
      wq=w(nn)*wq0
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint= y-k
      zint(n) = z-l
      ! --- computes coefficients for node centered quantities
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx1(n) = onesixth*oxintsq*oxint
      sx2(n) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx3(n) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx4(n) = onesixth*xintsq*xint
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy1(n) = onesixth*oyintsq*oyint*wq
      sy2(n) = (twothird-yintsq*(1.0_num-yint*0.5_num))*wq
      sy3(n) = (twothird-oyintsq*(1.0_num-oyint*0.5_num))*wq
      sy4(n) = onesixth*yintsq*yint*wq
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    DO n=1, MIN(LVEC2, np-ip+1)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !$DIR ASSUME_ALIGNED www:64, h1:64, h11:64, h12:64, zdec:64, sgn:64, szz:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, www, h1, h11, h12, zdec, sgn, szz)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 4!!! Vector
        zdec(nv)     = (h1(nv)-zint(n))*sgn(nv)
        szz(nv)      = (twothird-zdec(nv)**2*(1.0_num-zdec(nv)*0.5_num))*h11(nv)+     &
        onesixth*zdec(nv)**3*h12(nv)
        www(nv, n)    = szz(nv)*sy1(n)
        www(nv+4, n)  = szz(nv)*sy2(n)
        www(nv+8, n)  = szz(nv)*sy3(n)
        www(nv+12, n) = szz(nv)*sy4(n)
      ENDDO
    END DO
    ! Current deposition on vertices
    DO n=1, MIN(LVEC2, np-ip+1)
      ! --- add charge density contributions to vertices of the current cell
      ic=ICELL(n)
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED rhocells:64, www:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, rhocells, www)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 8!!! - VECTOR
        ! Loop on (i=-1, j, k)
        rhocells(nv, ic-ncx-1) = rhocells(nv, ic-ncx-1) + www(nv, n)*sx1(n)
        ! Loop on (i=0, j, k)
        rhocells(nv, ic-ncx)   = rhocells(nv, ic-ncx)   + www(nv, n)*sx2(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic-ncx+1) = rhocells(nv, ic-ncx+1) + www(nv, n)*sx3(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic-ncx+2) = rhocells(nv, ic-ncx+2) + www(nv, n)*sx4(n)
        ! Loop on (i=-1, j, k)
        rhocells(nv, ic+ncx-1) = rhocells(nv, ic+ncx-1) + www(nv+8, n)*sx1(n)
        ! Loop on (i=0, j, k)
        rhocells(nv, ic+ncx)   = rhocells(nv, ic+ncx)   + www(nv+8, n)*sx2(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic+ncx+1) = rhocells(nv, ic+ncx+1) + www(nv+8, n)*sx3(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic+ncx+2) = rhocells(nv, ic+ncx+2) + www(nv+8, n)*sx4(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  ! - reduction of rhocells in rho
  DO iz=1, ncz
    DO iy=1, ncy
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix=1, ncx!! VECTOR (take ncx multiple of vector length)
        ic=ix+(iy-1)*ncx+(iz-1)*ncxy
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
        rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1, ic)
        rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2, ic)
        rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3, ic)
        rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4, ic)
        rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5, ic)
        rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6, ic)
        rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7, ic)
        rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(rhocells)
  RETURN
END SUBROUTINE depose_rho_vecHVv3_3_3_3
#endif

! ________________________________________________________________________________________
!> @brief
!> Order 3 3D vector charge deposition routine
!
!> @details
!> Computes charge density on grid vectorized (HV-SCHEME)
!> This routine does vectorize on SIMD architecture with good performances
!> Speedup>2 on AVX 256 bits
!> lvect, the vector length was originally at 16
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] lvect vector length
! ________________________________________________________________________________________
SUBROUTINE depose_rho_vecHVv4_3_3_3(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx,  &
  dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, lvect)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp), INTENT (IN)    :: np, nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp), INTENT (IN)    :: lvect
  REAL(num), INTENT(IN OUT)     ::                                                    &
  rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT (IN)       :: xp(np), yp(np), zp(np), w(np)
  REAL(num), INTENT (IN)       :: q, dx, dy, dz, xmin, ymin, zmin

  REAL(num), DIMENSION(:, :), ALLOCATABLE:: rhocells
  INTEGER(idp), DIMENSION(lvect) :: ICELL
  INTEGER(idp) :: NCELLS
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint(1:lvect), oxint, oyint, ozint, xintsq, &
  yintsq, zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, invvol, wq0, wq
  REAL(num) :: sx1(lvect), sx2(lvect), sx3(lvect), sx4(lvect), sy1, sy2, sy3, sy4,    &
  sz1, sz2, sz3, sz4, w1, w2
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp):: ic, igrid, j, k, l, n, ip, nv, nn
  INTEGER(idp) :: nnx, nnxy
  INTEGER(idp) :: moff(1:8)
  REAL(num):: www1(lvect, 8), www2(lvect, 8)
  INTEGER(idp) :: orig, jorig, korig, lorig
  INTEGER(idp) :: ncx, ncy, ncxy, ncz, ix, iy, iz, ngridx, ngridy, ngx, ngxy

  ! Init parameters
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  wq0=q*invvol
  ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
  ncx=nx+5; ncy=ny+4; ncz=nz+2
  NCELLS=ncx*ncy*ncz
  ALLOCATE(rhocells(8, NCELLS))
  rhocells=0_num
  nnx = ngridx
  nnxy = ngridx*ngridy
  moff = (/-nnxy, 0_idp, nnxy, 2_idp*nnxy, nnx-nnxy, nnx, nnx+nnxy, nnx+2_idp*nnxy/)
  jorig=-2_idp; korig=-2_idp;lorig=-1_idp
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! FIRST LOOP: computes cell index of particle and their weight on vertices
  DO ip=1, np, lvect

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED w:64, sx1:64, sx2:64, sx3:64, sx4:64
    !DIR$ ASSUME_ALIGNED ICELL:64, www1:64, www2:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp, w, sx1, sx2, sx3, sx4, ICELL, www1, www2)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
    DO n=1, MIN(lvect, np-ip+1)
      nn=ip+n-1
      ! Calculation relative to particle n
      ! --- computes current position in grid units
      x= (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
      ! --- finds cell containing particles for current positions
      j=floor(x)
      k=floor(y)
      l=floor(z)
      ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
      wq=w(nn)*wq0
      ! --- computes distance between particle and node for current positions
      xint = x-j
      yint= y-k
      zint(n) = z-l
      ! --- computes coefficients for node centered quantities
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx1(n) = onesixth*oxintsq*oxint
      sx2(n) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx3(n) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx4(n) = onesixth*xintsq*xint
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy1 = onesixth*oyintsq*oyint
      sy2 = (twothird-yintsq*(1.0_num-yint*0.5_num))
      sy3 = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
      sy4 = onesixth*yintsq*yint
      ozint = 1.0_num-zint(n)
      zintsq = zint(n)*zint(n)
      ozintsq = ozint*ozint
      sz1 = onesixth*ozintsq*ozint*wq
      sz2 = (twothird-zintsq*(1.0_num-zint(n)*0.5_num))*wq
      sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
      sz4 = onesixth*zintsq*zint(n)*wq
      www1(n, 1)=sz1*sy1
      www1(n, 2)=sz2*sy1
      www1(n, 3)=sz3*sy1
      www1(n, 4)=sz4*sy1
      www1(n, 5)=sz1*sy2
      www1(n, 6)=sz2*sy2
      www1(n, 7)=sz3*sy2
      www1(n, 8)=sz4*sy2
      www2(n, 1)=sz1*sy3
      www2(n, 2)=sz2*sy3
      www2(n, 3)=sz3*sy3
      www2(n, 4)=sz4*sy3
      www2(n, 5)=sz1*sy4
      www2(n, 6)=sz2*sy4
      www2(n, 7)=sz3*sy4
      www2(n, 8)=sz4*sy4
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! Current deposition on vertices
    DO n=1, MIN(lvect, np-ip+1)
      ! --- add charge density contributions to vertices of the current cell
      ic=ICELL(n)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED rhocells:64, www1:64, www2:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, rhocells, www1, www2)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 8!!! - VECTOR
        w1=www1(n, nv)
        ! Loop on (i=-1, j, k)
        rhocells(nv, ic-ncx-1) = rhocells(nv, ic-ncx-1) + w1*sx1(n)
        ! Loop on (i=0, j, k)
        rhocells(nv, ic-ncx)   = rhocells(nv, ic-ncx)   + w1*sx2(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic-ncx+1) = rhocells(nv, ic-ncx+1) + w1*sx3(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic-ncx+2) = rhocells(nv, ic-ncx+2) + w1*sx4(n)

        w2=www2(n, nv)
        ! Loop on (i=-1, j, k)
        rhocells(nv, ic+ncx-1) = rhocells(nv, ic+ncx-1) + w2*sx1(n)
        ! Loop on (i=0, j, k)
        rhocells(nv, ic+ncx)   = rhocells(nv, ic+ncx)   + w2*sx2(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic+ncx+1) = rhocells(nv, ic+ncx+1) + w2*sx3(n)
        !Loop on (i=1, j, k)
        rhocells(nv, ic+ncx+2) = rhocells(nv, ic+ncx+2) + w2*sx4(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  ! - reduction of rhocells in rho
  DO iz=1, ncz
    DO iy=1, ncy
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix=1, ncx!! VECTOR (take ncx multiple of vector length)
        ic=ix+(iy-1)*ncx+(iz-1)*ncxy
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
        rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1, ic)
        rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2, ic)
        rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3, ic)
        rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4, ic)
        rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5, ic)
        rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6, ic)
        rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7, ic)
        rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(rhocells)
  RETURN
END SUBROUTINE depose_rho_vecHVv4_3_3_3


! ______________________________________________________________________________
!> @brief
!> General charge deposition routine (Warning: Highly unoptimized routine)
!
!> @details
!> Computes charge density on grid at arbitrary orders nox, noy and noz
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nx, ny, nz number of cells
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] nox, noy, noz interpolation order
!> @param[in] l_particle_weight flag to activate the use of the particle weight
!> @param[in] l4symtry
! ________________________________________________________________________________________
SUBROUTINE pxr_depose_rho_n(rho, np, xp, yp, zp, w, q, xmin, ymin, zmin, dx, dy, dz,  &
  nx, ny, nz, nxguard, nyguard, nzguard, nox, noy, noz, l_particles_weight, l4symtry)
  USE picsar_precision, ONLY: idp, lp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), intent(in out) :: rho
  REAL(num) :: xp(np), yp(np), zp(np), w(np)
  REAL(num) :: q, dx, dy, dz, xmin, ymin, zmin
  LOGICAL(lp) :: l_particles_weight, l4symtry
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq,  &
  zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, wq, invvol
  REAL(num) :: sx(-int(nox/2):int((nox+1)/2)), sy(-int(noy/2):int((noy+1)/2)),        &
  sz(-int(noz/2):int((noz+1)/2))
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, ip, jj, kk, ll, ixmin, ixmax, iymin, iymax, izmin, izmax

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi

  ixmin = -int(nox/2)
  ixmax = int((nox+1)/2)
  iymin = -int(noy/2)
  iymax = int((noy+1)/2)
  izmin = -int(noz/2)
  izmax = int((noz+1)/2)

  DO ip=1, np

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    ! --- applies 4-fold symmetry
    IF (l4symtry) THEN
      x=abs(x)
      y=abs(y)
    END IF

    ! --- finds node of cell containing particles for current positions
    ! --- (different for odd/even spline orders)
    IF (nox==2*(nox/2)) THEN
      j=nint(x)
    ELSE
      j=floor(x)
    END IF
    IF (noy==2*(noy/2)) THEN
      k=nint(y)
    ELSE
      k=floor(y)
    END IF
    IF(noz==2*(noz/2)) THEN
      l=nint(z)
    ELSE
      l=floor(z)
    END IF

    ! --- computes distance between particle and node for current positions
    xint = x-j
    yint = y-k
    zint = z-l

    ! --- computes particles "weights"
    IF (l_particles_weight) THEN
      wq=q*w(ip)*invvol
    ELSE
      wq=q*invvol*w(1)
    ENDIF

    ! --- computes coefficients for node centered quantities
    SELECT CASE(nox)
    CASE(0)
      sx( 0) = 1.0_num
    CASE(1)
      sx( 0) = 1.0_num-xint
      sx( 1) = xint
    CASE(2)
      xintsq = xint*xint
      sx(-1) = 0.5_num*(0.5_num-xint)**2
      sx( 0) = 0.75_num-xintsq
      sx( 1) = 0.5_num*(0.5_num+xint)**2
    CASE(3)
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
      sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
      sx( 2) = onesixth*xintsq*xint
    END SELECT

    SELECT CASE(noy)
    CASE(0)
      sy( 0) = 1.0_num
    CASE(1)
      sy( 0) = 1.0_num-yint
      sy( 1) = yint
    CASE(2)
      yintsq = yint*yint
      sy(-1) = 0.5_num*(0.5_num-yint)**2
      sy( 0) = 0.75_num-yintsq
      sy( 1) = 0.5_num*(0.5_num+yint)**2
    CASE(3)
      oyint = 1.0_num-yint
      yintsq = yint*yint
      oyintsq = oyint*oyint
      sy(-1) = onesixth*oyintsq*oyint
      sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
      sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
      sy( 2) = onesixth*yintsq*yint
    END SELECT

    SELECT CASE(noz)
    CASE(0)
      sz( 0) = 1.0_num
    CASE(1)
      sz( 0) = 1.0_num-zint
      sz( 1) = zint
    CASE(2)
      zintsq = zint*zint
      sz(-1) = 0.5_num*(0.5_num-zint)**2
      sz( 0) = 0.75_num-zintsq
      sz( 1) = 0.5_num*(0.5_num+zint)**2
    CASE(3)
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
      sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
      sz( 2) = onesixth*zintsq*zint
    END SELECT

    ! --- add charge density contributions
    DO ll = izmin, izmax
      DO kk = iymin, iymax
        DO jj = ixmin, ixmax
          rho(j+jj, k+kk, l+ll)=rho(j+jj, k+kk, l+ll)+sx(jj)*sy(kk)*sz(ll)*wq
        END DO
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE pxr_depose_rho_n

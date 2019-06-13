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
! DIRECT_CURRENT_DEPOSITION_2D.F90
!
! Developers
! Remi Lehe
!
! Description:
! This file contains subroutines for the direct current deposition in 2D.
!
! List of subroutines:
!
! Classical non-optimized:
! - depose_jxjyjz_scalar2d_1_1_1
! - depose_jxjyjz_scalar2d_2_2_2
! - depose_jxjyjz_scalar2d_3_3_3
!
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Order 1 2D scalar direct current deposition routine (rho*v)
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!
!> @param[in] np Number of particles
!> @param[in] xp 1D array of x-coordinates of particles
!> @param[in] zp 1D array of x-coordinates of particles
!> @param[in] uxp 1D array of ux-velocity components of particles
!> @param[in] uyp 1D array of ux-velocity components of particles
!> @param[in] uzp 1D array of ux-velocity components of particles
!> @param[in] gaminv 1D array of the inverse 1/gamma-factor of particles
!> @param[in] w 1D array of the weghts of particles
!> @param[in] q charge of current species (scalar)
!> @param[in] xmin x-minimum boundary of current tile
!> @param[in] zmin z-minimum boundary of current tile
!> @param[in] dt time step (scalar)
!> @param[in] dx mesh size along x (scalar)
!> @param[in] dz mesh size along z (scalar)
!> @param[inout] jx x-current component (2D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 2 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jx array (1d array containing 2 integers)
!> @param[inout] jy y-current component (2D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!>  (1d array containing 2 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 2 integers)
!> @param[inout] jz z-current component (2D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 2 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 2 integers)
!> @warning arrays jx, jy, jz should be set to 0 before entering this subroutine.
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_scalar2d_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp)             :: np
  LOGICAL(idp)             :: l_nodal
  REAL(num)                :: stagger_shift
  INTEGER(idp), intent(in) :: jx_nguard(2), jx_nvalid(2), jy_nguard(2), jy_nvalid(2), &
  jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1)
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1)
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1)
  REAL(num), DIMENSION(np) :: xp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                :: q, dt, dx, dz, xmin, zmin
  REAL(num)                :: dxi, dzi, xint, zint
  REAL(num)                :: x, z, xmid, zmid, vx, vy, vz, invvol, dts2dx, dts2dz
  REAL(num)                :: wq, wqx, wqy, wqz, clightsq
  REAL(num), DIMENSION(2)  :: sx(0:1), sz(0:1), sx0(0:1), sz0(0:1)
  REAL(num), PARAMETER     :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp)             :: j, l, j0, l0, ip

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  invvol = dxi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
!  sx=0.0_num;sz=0.0_num;
!  sx0=0.0_num;sz0=0.0_num;

  ! LOOP ON PARTICLES
  ! Prevent loop to vectorize (dependencies)
  !DIR$ NOVECTOR
!$acc parallel deviceptr(jx, jy, jz, xp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(0:1), sz(0:1), sx0(0:1), sz0(0:1) )
  DO ip=1, np

    ! --- computes position in  grid units at (n+1)
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi

    ! Computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes particles weights
    wq=q*w(ip)
    wqx=wq*invvol*vx
    wqy=wq*invvol*vy
    wqz=wq*invvol*vz

    ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
    xmid=x-dts2dx*vx
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j=floor(xmid)
    l=floor(zmid)
    j0=floor(xmid-stagger_shift)
    l0=floor(zmid-stagger_shift)

    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    zint = zmid-l
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-stagger_shift-j0
    zint = zmid-stagger_shift-l0
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
    ! - JX
    !$acc atomic update
    jx(j0, l  )    = jx(j0, l  )  +   sx0(0)*sz(0)*wqx
    !$acc atomic update
    jx(j0+1, l  )    = jx(j0+1, l  )  +   sx0(1)*sz(0)*wqx
    !$acc atomic update
    jx(j0, l+1)    = jx(j0, l+1)  +   sx0(0)*sz(1)*wqx
    !$acc atomic update
    jx(j0+1, l+1)    = jx(j0+1, l+1)  +   sx0(1)*sz(1)*wqx

    ! - JY
    !$acc atomic update
    jy(j, l  )    = jy(j, l  )  +   sx(0)*sz(0)*wqy
    !$acc atomic update
    jy(j+1, l  )    = jy(j+1, l  )  +   sx(1)*sz(0)*wqy
    !$acc atomic update
    jy(j, l+1)    = jy(j, l+1)  +   sx(0)*sz(1)*wqy
    !$acc atomic update
    jy(j+1, l+1)    = jy(j+1, l+1)  +   sx(1)*sz(1)*wqy

    ! - JZ
    !$acc atomic update
    jz(j, l0  )    = jz(j, l0  )  +   sx(0)*sz0(0)*wqz
    !$acc atomic update
    jz(j+1, l0  )    = jz(j+1, l0  )  +   sx(1)*sz0(0)*wqz
    !$acc atomic update
    jz(j, l0+1)    = jz(j, l0+1)  +   sx(0)*sz0(1)*wqz
    !$acc atomic update
    jz(j+1, l0+1)    = jz(j+1, l0+1)  +   sx(1)*sz0(1)*wqz
  END DO
!$acc end loop
!$acc end parallel

  RETURN
END SUBROUTINE depose_jxjyjz_scalar2d_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 2 2D scalar current deposition routine (jx*v)
!
!> @details
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2015
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_scalar2d_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np
  LOGICAL(idp) :: l_nodal
  REAL(num)    :: stagger_shift
  INTEGER(idp), intent(in) :: jx_nguard(2), jx_nvalid(2), jy_nguard(2), jy_nvalid(2), &
  jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1)
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1)
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1)
  REAL(num), DIMENSION(np) :: xp, zp, uxp, uyp, uzp, gaminv, w
  REAL(num) :: q, dt, dx, dz, xmin, zmin
  REAL(num) :: dxi, dzi, xint, zint
  REAL(num) :: xintsq, zintsq
  REAL(num) :: x, z, xmid, zmid, vx, vy, vz, invvol, dts2dx, dts2dz
  REAL(num) :: wq, wqx, wqy, wqz, clightsq
  REAL(num), DIMENSION(3) :: sx(-1:1), sz(-1:1), sx0(-1:1),      &
  sz0(-1:1)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, l, j0, l0, ip

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  invvol = dxi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
  sx=0.0_num;sz=0.0_num;
  sx0=0.0_num;sz0=0.0_num;

  ! LOOP ON PARTICLES
  !$acc parallel deviceptr(jx, jy, jz, xp, zp, uxp, uyp, uzp, w, gaminv)
  !$acc loop gang vector private(sx(-1:1), sz(-1:1), sx0(-1:1), sz0(-1:1))
  DO ip=1, np
    ! --- computes position in  grid units at (n+1)
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi

    ! Computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes particles weights
    wq=q*w(ip)*invvol
    wqx=wq*vx
    wqy=wq*vy
    wqz=wq*vz

    ! Gets position in grid units at (n+1/2) for computing jx(n+1/2)
    xmid=x-dts2dx*vx
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j=nint(xmid)
    l=nint(zmid)
    j0=nint(xmid-stagger_shift)
    l0=nint(zmid-stagger_shift)
    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    zint = zmid-l
    xintsq = xint*xint
    sx(-1) = 0.5_num*(0.5_num-xint)**2
    sx( 0) = 0.75_num-xintsq
    sx( 1) = 0.5_num*(0.5_num+xint)**2
    zintsq = zint*zint
    sz(-1) = 0.5_num*(0.5_num-zint)**2
    sz( 0) = (0.75_num-zintsq)
    sz( 1) = 0.5_num*(0.5_num+zint)**2

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-stagger_shift-j0
    zint = zmid-stagger_shift-l0
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2
    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = (0.75_num-zintsq)
    sz0( 1) = 0.5_num*(0.5_num+zint)**2

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
    ! - JX
    !$acc atomic update
    jx(j0-1, l-1)  = jx(j0-1, l-1)  +   sx0(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0, l-1)  = jx(j0, l-1)  +   sx0(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, l-1)  = jx(j0+1, l-1)  +   sx0(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, l  )  = jx(j0-1, l  )  +   sx0(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0, l  )  = jx(j0, l  )  +   sx0(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, l  )  = jx(j0+1, l  )  +   sx0(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, l+1)  = jx(j0-1, l+1)  +   sx0(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0, l+1)  = jx(j0, l+1)  +   sx0(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, l+1)  = jx(j0+1, l+1)  +   sx0(1 )*sz(1 )*wqx

    !        ! - JY
    !$acc atomic update
    jy(j-1, l-1)  = jy(j-1, l-1)  +   sx(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j, l-1)  = jy(j, l-1)  +   sx(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, l-1)  = jy(j+1, l-1)  +   sx(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, l  )  = jy(j-1, l  )  +   sx(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j, l  )  = jy(j, l  )  +   sx(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, l  )  = jy(j+1, l  )  +   sx(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, l+1)  = jy(j-1, l+1)  +   sx(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j, l+1)  = jy(j, l+1)  +   sx(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, l+1)  = jy(j+1, l+1)  +   sx(1 )*sz(1 )*wqy

    ! - JZ
    !$acc atomic update
    jz(j-1, l0-1)  = jz(j-1, l0-1)  +   sx(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j, l0-1)  = jz(j, l0-1)  +   sx(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, l0-1)  = jz(j+1, l0-1)  +   sx(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, l0  )  = jz(j-1, l0  )  +   sx(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j, l0  )  = jz(j, l0  )  +   sx(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, l0  )  = jz(j+1, l0  )  +   sx(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, l0+1)  = jz(j-1, l0+1)  +   sx(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j, l0+1)  = jz(j, l0+1)  +   sx(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, l0+1)  = jz(j+1, l0+1)  +   sx(1 )*sz0(1 )*wqz
  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_scalar2d_2_2_2

! ______________________________________________________________________________
!> @brief
!> Order 3 2D scalar current deposition routine (rho*v)
!
!> @details
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_scalar2d_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np
  LOGICAL(idp) :: l_nodal
  REAL(num)    :: stagger_shift
  INTEGER(idp), intent(in) :: jx_nguard(2), jx_nvalid(2), jy_nguard(2), jy_nvalid(2), &
  jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1)
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1)
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1)
  REAL(num), DIMENSION(np) :: xp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dz, xmin, zmin
  REAL(num) :: dxi, dzi, xint, zint, oxint, ozint, xintsq,  &
  zintsq, oxintsq, ozintsq
  REAL(num) :: x, z, xmid, zmid, vx, vy, vz, invvol, dts2dx, dts2dz
  REAL(num) :: wq, wqx, wqy, wqz, clightsq
  REAL(num), DIMENSION(4) :: sx(-1:2), sz(-1:2), sx0(-1:2),      &
  sz0(-1:2)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, l, j0, l0, ip


  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  invvol = dxi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
  sx=0.0_num;sz=0.0_num;
  sx0=0.0_num;sz0=0.0_num;

  ! LOOP ON PARTICLES
  !$acc parallel deviceptr(jx, jy, jz, xp, zp, uxp, uyp, uzp, w, gaminv)
  !$acc loop gang vector private(sx(-1:2), sz(-1:2), sx0(-1:2), sz0(-1:2))
  DO ip=1, np
    ! --- computes position in  grid units at (n+1)
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi

    ! Computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes particles weights
    wq=q*w(ip)*invvol
    wqx=wq*vx
    wqy=wq*vy
    wqz=wq*vz

    ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
    xmid=x-dts2dx*vx
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j=floor(xmid)
    l=floor(zmid)
    j0=floor(xmid-stagger_shift)
    l0=floor(zmid-stagger_shift)

    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    zint = zmid-l
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1) = onesixth*oxintsq*oxint
    sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx( 2) = onesixth*xintsq*xint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1) = onesixth*ozintsq*ozint
    sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz( 2) = onesixth*zintsq*zint

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-stagger_shift-j0
    zint = zmid-stagger_shift-l0
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx0(-1) = onesixth*oxintsq*oxint
    sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx0( 2) = onesixth*xintsq*xint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz0(-1) = onesixth*ozintsq*ozint
    sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz0( 2) = onesixth*zintsq*zint

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
    ! - JX
    !$acc atomic update
    jx(j0-1, l-1)  = jx(j0-1, l-1)  +   sx0(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0, l-1)  = jx(j0, l-1)  +   sx0(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, l-1)  = jx(j0+1, l-1)  +   sx0(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+2, l-1)  = jx(j0+2, l-1)  +   sx0(2 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, l  )  = jx(j0-1, l  )  +   sx0(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0, l  )  = jx(j0, l  )  +   sx0(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, l  )  = jx(j0+1, l  )  +   sx0(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+2, l  )  = jx(j0+2, l  )  +   sx0(2 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, l+1)  = jx(j0-1, l+1)  +   sx0(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0, l+1)  = jx(j0, l+1)  +   sx0(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, l+1)  = jx(j0+1, l+1)  +   sx0(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+2, l+1)  = jx(j0+2, l+1)  +   sx0(2 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, l+2)  = jx(j0-1, l+2)  +   sx0(-1)*sz(2 )*wqx
    !$acc atomic update
    jx(j0, l+2)  = jx(j0, l+2)  +   sx0(0 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+1, l+2)  = jx(j0+1, l+2)  +   sx0(1 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+2, l+2)  = jx(j0+2, l+2)  +   sx0(2 )*sz(2 )*wqx

    ! - JY
    !$acc atomic update
    jy(j-1, l-1)  = jy(j-1, l-1)  +   sx(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j, l-1)  = jy(j, l-1)  +   sx(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, l-1)  = jy(j+1, l-1)  +   sx(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+2, l-1)  = jy(j+2, l-1)  +   sx(2 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, l  )  = jy(j-1, l  )  +   sx(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j, l  )  = jy(j, l  )  +   sx(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, l  )  = jy(j+1, l  )  +   sx(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+2, l  )  = jy(j+2, l  )  +   sx(2 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, l+1)  = jy(j-1, l+1)  +   sx(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j, l+1)  = jy(j, l+1)  +   sx(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, l+1)  = jy(j+1, l+1)  +   sx(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+2, l+1)  = jy(j+2, l+1)  +   sx(2 )*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, l+2)  = jy(j-1, l+2)  +   sx(-1)*sz(2 )*wqy
    !$acc atomic update
    jy(j, l+2)  = jy(j, l+2)  +   sx(0 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+1, l+2)  = jy(j+1, l+2)  +   sx(1 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+2, l+2)  = jy(j+2, l+2)  +   sx(2 )*sz(2 )*wqy

    ! - JZ
    !$acc atomic update
    jz(j-1, l0-1)  = jz(j-1, l0-1)  +   sx(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j, l0-1)  = jz(j, l0-1)  +   sx(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, l0-1)  = jz(j+1, l0-1)  +   sx(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+2, l0-1)  = jz(j+2, l0-1)  +   sx(2 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, l0  )  = jz(j-1, l0  )  +   sx(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j, l0  )  = jz(j, l0  )  +   sx(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, l0  )  = jz(j+1, l0  )  +   sx(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+2, l0  )  = jz(j+2, l0  )  +   sx(2 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, l0+1)  = jz(j-1, l0+1)  +   sx(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j, l0+1)  = jz(j, l0+1)  +   sx(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, l0+1)  = jz(j+1, l0+1)  +   sx(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+2, l0+1)  = jz(j+2, l0+1)  +   sx(2 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, l0+2)  = jz(j-1, l0+2)  +   sx(-1)*sz0(2 )*wqz
    !$acc atomic update
    jz(j, l0+2)  = jz(j, l0+2)  +   sx(0 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+1, l0+2)  = jz(j+1, l0+2)  +   sx(1 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+2, l0+2)  = jz(j+2, l0+2)  +   sx(2 )*sz0(2 )*wqz

  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_scalar2d_3_3_3

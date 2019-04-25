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
! ESIRKEPOV_CURRENT_DEPOSITION_2D.F90
!
! Developers
! Henri Vincenti, ! Mathieu Lobet
!
! Description:
! This file contains subroutines for Esirkepov current deposition in 2D.
!
! List of subroutines:
!
! - pxr_depose_jxjyjz_esirkepov2d_n
! - pxr_depose_jxjyjz_esirkepov2d_1_1
! - pxr_depose_jxjyjz_esirkepov2d_2_2
! - pxr_depose_jxjyjz_esirkepov2d_3_3
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition esirkepov n order (from 0 to 3)
!
!> @details
!> This subroutine is adapted from the version of WARP.
!> This subroutine is called in pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp().

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2016

! Input parameters:
!> @param[inout] jx x-current component (2D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 2 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the jx array
!> (1d array containing 2 integers)
!> @param[inout] jy y-current component (2D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!>  (1d array containing 2 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the jy array
!> (1d array containing 2 integers)
!> @param[inout] jz z-current component (2D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 2 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the jz array
!> (1d array containing 2 integers)
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________
subroutine pxr_depose_jxjyjz_esirkepov2d_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,  &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, nox, noz, l_particles_weight, l4symtry, l_2drz,               &
  type_rz_depose)      !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp) :: np, nox, noz, type_rz_depose
  INTEGER(idp), intent(in)              :: jx_nguard(2), jx_nvalid(2), jy_nguard(2),  &
  jy_nvalid(2), jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1 )
  real(num), dimension(np) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  real(num)    :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp)  :: l_particles_weight, l4symtry, l_2drz

  real(num) :: dxi, dzi, dtsdx, dtsdz, xint, zint
  real(num), dimension(:, :), allocatable :: sdx, sdz
  real(num) :: xold, yold, zold, rold, x, y, z, r, c, s, wq, wqx, wqz
  real(num) :: vx, vy, vz, dts2dx, dts2dz
  real(num) :: invvol, invdtdx, invdtdz
  real(num) :: oxint, ozint, xintsq, zintsq, oxintsq, ozintsq
  real(num) :: dtsdx0, dtsdz0, dts2dx0, dts2dz0
  real(num), parameter :: onesixth=1./6.
  real(num), parameter :: twothird=2./3.
  real(num), dimension(:), allocatable :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp) :: iixp0, ikxp0, iixp, ikxp, ip, dix, diz, i, k, ic, kc
  integer(idp) :: ixmin, ixmax, izmin, izmax, icell, ncells, ndtodx, ndtodz
  integer(idp) :: xl, xu, zl, zu

  ndtodx = int(clight*dt/dx)
  ndtodz = int(clight*dt/dz)
  xl = -int(nox/2)-1-ndtodx
  xu = int((nox+1)/2)+1+ndtodx
  zl = -int(noz/2)-1-ndtodz
  zu = int((noz+1)/2)+1+ndtodz
  allocate(sdx(xl:xu, zl:zu), sdz(xl:xu, zl:zu))
  allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
  allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

  sx0=0.;sz0=0.
  sdx=0.;sdz=0.

  ! Davoine method : limited to order 1 in r
  if (type_rz_depose==2) then
    nox = 1
  endif

  dxi = 1./dx
  dzi = 1./dz
  invvol = 1./(dx*dz)
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  dts2dx0 = 0.5*dtsdx0
  dts2dz0 = 0.5*dtsdz0
  invdtdx = 1./(dt*dz)
  invdtdz = 1./(dt*dx)

  do ip=1, np

    ! --- computes current position in grid units
    x = xp(ip)
    if (l_2drz) then
      y = yp(ip)
      r=sqrt(x*x+y*y)
      if (r*dxi>1.e-10) then
        c = x/r
        s = y/r
      else
        c = 1.
        s = 0.
      end if
      x = r
    end if
    x=x*dxi
    z = zp(ip)*dzi

    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes old position in grid units
    if (l_2drz) then
      xold = xp(ip)-dt*vx
      yold = yp(ip)-dt*vy
      rold = sqrt(xold*xold+yold*yold)
      xold=rold*dxi
      vy = -vx*s+vy*c
      vx = (x-xold)/dtsdx0
    else
      xold=x-dtsdx0*vx
    end if
    zold=z-dtsdz0*vz

    ! --- applies 4-fold symmetry
    if (l4symtry) then
      x=abs(x)
      xold=abs(xold)
      vx = (x-xold)/dtsdx0
    end if

    ! --- sets positions relative to grid  start
    x = x-xmin*dxi
    z = z-zmin*dzi
    xold = xold-xmin*dxi
    zold = zold-zmin*dzi

    ! computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    z=zold

    do icell = 1, ncells

      xold = x
      zold = z

      x = x+dtsdx*vx
      z = z+dtsdz*vz

      ! --- computes particles "weights"
      if (l_particles_weight) then
        wq=q*w(ip)
      else
        wq=q*w(1)
      end if
      wqx = wq*invdtdx
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      ! --- (different for odd/even spline orders)
      if (nox==2*(nox/2)) then
        iixp0=nint(x)
      else
        iixp0=floor(x)
      end if
      if (noz==2*(noz/2)) then
        ikxp0=nint(z)
      else
        ikxp0=floor(z)
      end if

      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
      if (type_rz_depose == 2) then! Davoine method, modified particle shapes in r
        sx0(0) = 1. - xint  + 1./(4*iixp0+2)*( -xint + xint**2 )
        sx0(1) = 1. - sx0(0)
      else! Standard method, canonical shapes in r
        select case(nox)
        case(0)
          sx0( 0) = 1.
        case(1)
          sx0( 0) = 1.-xint
          sx0( 1) = xint
        case(2)
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end select
      endif

      select case(noz)
      case(0)
        sz0( 0) = 1.
      case(1)
        sz0( 0) = 1.-zint
        sz0( 1) = zint
      case(2)
        zintsq = zint*zint
        sz0(-1) = 0.5*(0.5-zint)**2
        sz0( 0) = 0.75-zintsq
        sz0( 1) = 0.5*(0.5+zint)**2
      case(3)
        ozint = 1.-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.-zint/2)
        sz0( 1) = twothird-ozintsq*(1.-ozint/2)
        sz0( 2) = onesixth*zintsq*zint
      end select

      ! --- finds node of cell containing particles for old positions
      ! --- (different for odd/even spline orders)
      if (nox==2*(nox/2)) then
        iixp=nint(xold)
      else
        iixp=floor(xold)
      end if
      if (noz==2*(noz/2)) then
        ikxp=nint(zold)
      else
        ikxp=floor(zold)
      end if

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients
      ! --- (needed because of different dix and diz for each particle)
      sx=0.;sz=0.

      ! --- computes coefficients for quantities centered between nodes
      if (type_rz_depose == 2) then! Davoine method, modified particle shapes in r
        sx(0) = 1. - xint  + 1./(4*iixp+2)*( -xint + xint**2 )
        sx(1) = 1. - sx(0)
      else! Standard method, canonical shapes in r
        select case(nox)
        case(0)
          sx( 0+dix) = 1.
        case(1)
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
        case(2)
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
        case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end select
      endif

      select case(noz)
      case(0)
        sz( 0+diz) = 1.
      case(1)
        sz( 0+diz) = 1.-zint
        sz( 1+diz) = zint
      case(2)
        zintsq = zint*zint
        sz(-1+diz) = 0.5*(0.5-zint)**2
        sz( 0+diz) = 0.75-zintsq
        sz( 1+diz) = 0.5*(0.5+zint)**2
      case(3)
        ozint = 1.-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1+diz) = onesixth*ozintsq*ozint
        sz( 0+diz) = twothird-zintsq*(1.-zint/2)
        sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
        sz( 2+diz) = onesixth*zintsq*zint
      end select

      ! --- computes coefficients difference
      dsx = sx - sx0
      dsz = sz - sz0

      ! --- computes min/max positions of current contributions
      ixmin = min(0_idp, dix)-int(nox/2)
      ixmax = max(0_idp, dix)+int((nox+1)/2)
      izmin = min(0_idp, diz)-int(noz/2)
      izmax = max(0_idp, diz)+int((noz+1)/2)

      ! --- add current contributions
      ! --- NB : the current is later divided by the cylindrical cell volume in applybc_j
      do k=izmin, izmax
        do i=ixmin, ixmax
          ic = iixp0+i
          kc = ikxp0+k

          ! -- Jx
          if(i<ixmax) then
            sdx(i, k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )! Wx coefficient from esirkepov
            if (i>ixmin) sdx(i, k)=sdx(i, k)+sdx(i-1, k)! Integration of Wx along x
            jx(ic, kc) = jx(ic, kc) + sdx(i, k)! Deposition on the current
          end if

          ! -- Jy (2D Esirkepov scheme)
          jy(ic, kc) = jy(ic, kc) + wq*vy*invvol/ncells* ( (sz0(k)+0.5*dsz(k))*sx0(i) &
          + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )

          ! -- Jz
          if(k<izmax) then
            sdz(i, k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))! Wz coefficient from esirkepov
            if (k>izmin) sdz(i, k)=sdz(i, k)+sdz(i, k-1)! Integration of Wz along z
            jz(ic, kc) = jz(ic, kc) + sdz(i, k)! Deposition on the current
          end if
        end do
      end do
    end do
  end do

  deallocate(sdx, sdz, sx, sx0, dsx, sz, sz0, dsz)

  return
END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_n

! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 1
!
!> @details
!> This function is not optimized but provides better performances than
!> using the abitrary order function

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2016

!
! Input parameters:
!> @param[inout] jx x-current component (2D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 2 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the jx array
!> (1d array containing 2 integers)
!> @param[inout] jy y-current component (2D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!> (1d array containing 2 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the jy array
!>  (1d array containing 2 integers)
!> @param[inout] jz z-current component (2D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 2 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the jz array
!>  (1d array containing 2 integers)
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] lvect: vector length (useless since not vectorized)
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here but kept for common interface)
!> @param[in] l_2drz  (useless here but kept for common interface)
!> @param[in] type_rz_depose (useless here but kept for common interface)
! ________________________________________________________________________________________
SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1( jx, jx_nguard, jx_nvalid, jy,           &
  jy_nguard, jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, &
  q, xmin, zmin, dt, dx, dz, lvect, l_particles_weight, l4symtry, l_2drz,               &
  type_rz_depose)      !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  implicit none
  ! __ Parameter declaration ________________________________________________________
  integer(idp)                          :: np, type_rz_depose
  integer(idp)                          :: lvect
  INTEGER(idp), intent(in)              :: jx_nguard(2), jx_nvalid(2), jy_nguard(2),  &
  jy_nvalid(2), jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1 )
  real(num), dimension(np)              :: xp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp)                           :: l_particles_weight, l4symtry, l_2drz
  real(num)                             :: dxi, dzi, xint, zint
  real(num)                             :: xold, zold, x, z, wq, wqx, wqz
  real(num)                             :: vx, vy, vz
  real(num)                             :: invvol, invdtdx, invdtdz
  real(num)                             :: dtsdx0, dtsdz0
  real(num), parameter                  :: onesixth=1./6., twothird=2./3.
  real(num), parameter                  :: onethird=1./3.
  integer(idp)                          :: iixp0, ikxp0, iixp, ikxp, ip, i, k, ic, kc
  integer(isp)                          :: dix, diz
  integer(isp)                          :: ixmin, ixmax, izmin, izmax

  real(num), dimension(4) :: sx(-1:2), sx0(-1:2), dsx(-1:2)
  real(num), dimension(4) :: sz(-1:2), sz0(-1:2), dsz(-1:2)
  real(num)               :: sdxi, sdxim1
  real(num), dimension(4) :: sdzk(-1:2), sdzkm1(-1:2)

  ! __ Parameter initialization ______________________________
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  sdxi   = 0.0_num
  sdxim1 = 0.0_num
  sdzk   = 0.0_num
  sdzkm1 = 0.0_num

!$acc parallel deviceptr(jx, jy, jz, xp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(-1:2), sz(-1:2), &
!$acc&                         sx0(-1:2), sz0(-1:2), &
!$acc&                         dsx(-1:2), dsz(-1:2), &
!$acc&                         sdxi, sdxim1, &
!$acc&                         sdzk(-1:2), sdzkm1(-1:2) )
  DO ip=1, np
     sx  = 0.0_num
     sz  = 0.0_num
     sx0 = 0.0_num
     sz0 = 0.0_num

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi

    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz

    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz

    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ikxp0=floor(z)

    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint

    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint

    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ikxp=floor(zold)

    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp

    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0

    ! --- zero out coefficients
    ! --- (needed because of different dix and diz for each particle)
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num

    ! --- computes coefficients for quantities centered between nodes
    sx( 0+dix) = 1.0_num-xint
    sx( 1+dix) = xint

    sz( 0+diz) = 1.0_num-zint
    sz( 1+diz) = zint

    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0

    ! --- computes min/max positions of current contributions
    ixmin = min(0, dix)-0
    ixmax = max(0, dix)+1
    izmin = min(0, diz)-0
    izmax = max(0, diz)+1

    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k

        ! --- Jx
        IF(i<ixmax) THEN
          sdxi = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) ) ! Wx coefficient from esirkepov
          if (i>ixmin) sdxi = sdxi + sdxim1 ! Integration of Wx along x
           !$acc atomic update
           jx(ic, kc) = jx(ic, kc) + sdxi ! Deposition on the current
        END IF

        ! -- Jy (2D Esirkepov scheme)
        !$acc atomic update
        jy(ic, kc) = jy(ic, kc) + wq*vy*invvol* ( (sz0(k)+0.5*dsz(k))*sx0(i) +        &
        (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )

        ! --- Jz
        IF(k<izmax) THEN
          sdzk(i) = wqz*dsz(k)*(sx0(i)+0.5*dsx(i)) ! Wz coefficient from esirkepov&
          if (k>izmin) sdzk(i) = sdzk(i) + sdzkm1(i) ! Integration of Wz along z
            !$acc atomic update
          jz(ic, kc) = jz(ic, kc) + sdzk(i) ! Deposition on the current
        END IF
        sdxim1 = sdxi
      END DO
      sdzkm1 = sdzk
    END DO

  END DO
!$acc end loop
!$acc end parallel
  RETURN

END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1

! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 2
!> This function is not optimized but provides better performances than
!> using the abitrary order function
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx x-current component (2D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 2 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the jx array
!> (1d array containing 2 integers)
!> @param[inout] jy y-current component (2D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!> (1d array containing 2 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the jy array
!>  (1d array containing 2 integers)
!> @param[inout] jz z-current component (2D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 2 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the jz array
!>  (1d array containing 2 integers)
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________
SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_2_2( jx, jx_nguard, jx_nvalid, jy,           &
  jy_nguard, jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, &
  q, xmin, zmin, dt, dx, dz, lvect, l_particles_weight, l4symtry, l_2drz,               &
  type_rz_depose)      !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  implicit none
  ! __ Parameter declaration ________________________________________________________
  integer(idp)                          :: np, type_rz_depose
  integer(idp)                          :: lvect
  INTEGER(idp), intent(in)              :: jx_nguard(2), jx_nvalid(2), jy_nguard(2),  &
  jy_nvalid(2), jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1 )
  real(num), dimension(np)              :: xp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp)                           :: l_particles_weight, l4symtry, l_2drz
  real(num)                             :: dxi, dzi, xint, zint
  real(num)                             :: xold, zold, x, z, wq, wqx, wqz
  real(num)                             :: vx, vy, vz
  real(num)                             :: invvol, invdtdx, invdtdz
  real(num)                             :: xintsq, zintsq
  real(num)                             :: dtsdx0, dtsdz0
  real(num), parameter                  :: onesixth=1./6., twothird=2./3.
  real(num), parameter                  :: onethird=1./3.
  integer(idp)                          :: iixp0, ikxp0, iixp, ikxp, ip, i, k, ic, kc
  integer(isp)                          :: dix, diz
  integer(isp)                          :: ixmin, ixmax, izmin, izmax

  real(num), dimension(5) :: sx(-2:2), sx0(-2:2), dsx(-2:2)
  real(num), dimension(5) :: sz(-2:2), sz0(-2:2), dsz(-2:2)
  real(num)               :: sdxi, sdxim1
  real(num), dimension(5) :: sdzk(-2:2), sdzkm1(-2:2)

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  sdxi   = 0.0_num
  sdxim1 = 0.0_num
  sdzk   = 0.0_num
  sdzkm1 = 0.0_num
!$acc parallel deviceptr(jx, jy, jz, xp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(-2:2), sz(-2:2), &
!$acc&                         sx0(-2:2), sz0(-2:2), &
!$acc&                         dsx(-2:2), dsz(-2:2), &
!$acc&                         sdxi, sdxim1, &
!$acc&                         sdzk(-2:2), sdzkm1(-2:2) )
  DO ip=1, np

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi

    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz

    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz

    ! --- finds node of cell containing particles for current positions
    iixp0=nint(x)
    ikxp0=nint(z)

    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2

    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = 0.75_num-zintsq
    sz0( 1) = 0.5_num*(0.5_num+zint)**2

    ! --- finds node of cell containing particles for old positions
    iixp=nint(xold)
    ikxp=nint(zold)

    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp

    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0

    ! --- zero out coefficients
    ! --- (needed because of different dix and diz for each particle)
    sx(-2)=0.0_num
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sz(-2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num

    ! --- computes coefficients for quantities centered between nodes
    xintsq = xint*xint
    sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
    sx( 0+dix) = 0.75_num-xintsq
    sx( 1+dix) = 0.5_num*(0.5_num+xint)**2

    zintsq = zint*zint
    sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
    sz( 0+diz) = 0.75_num-zintsq
    sz( 1+diz) = 0.5_num*(0.5_num+zint)**2

    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0

    ! --- computes min/max positions of current contributions
    ixmin = min(0, dix)-1
    ixmax = max(0, dix)+1
    izmin = min(0, diz)-1
    izmax = max(0, diz)+1

    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k

        ! --- Jx
        IF(i<ixmax) THEN
          sdxi = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) ) ! Wx coefficient from esirkepov
          if (i>ixmin) sdxi = sdxi + sdxim1 ! Integration of Wx along x
          !$acc atomic update
          jx(ic, kc) = jx(ic, kc) + sdxi ! Deposition on the current
        END IF

        ! -- Jy (2D Esirkepov scheme)
        !$acc atomic update
        jy(ic, kc) = jy(ic, kc) + wq*vy*invvol* ( (sz0(k)+0.5*dsz(k))*sx0(i) +        &
        (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )

        ! --- Jz
        IF(k<izmax) THEN
          sdzk(i) = wqz*dsz(k)*(sx0(i)+0.5*dsx(i)) ! Wz coefficient from esirkepov&
          if (k>izmin) sdzk(i) = sdzk(i) + sdzkm1(i) ! Integration of Wz along z
          !$acc atomic update
          jz(ic, kc) = jz(ic, kc) + sdzk(i) ! Deposition on the current
        END IF
        sdxim1 = sdxi
      END DO
      sdzkm1 = sdzk
    END DO
  END DO
!$acc end loop
!$acc end parallel
  RETURN

END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_2_2

! ______________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 3
!> This function is not optimized but provides better performances than
!> using the abitrary order function
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx x-current component (2D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 2 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the jx array
!> (1d array containing 2 integers)
!> @param[inout] jy y-current component (2D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!> (1d array containing 2 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 2 integers)
!> @param[inout] jz z-current component (2D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 2 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells)
!> of the jz array (1d array containing 2 integers)
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
!
! ________________________________________________________________________________________
subroutine pxr_depose_jxjyjz_esirkepov2d_3_3( jx, jx_nguard, jx_nvalid, jy,           &
  jy_nguard, jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, zp, uxp, uyp, uzp, gaminv, w, &
  q, xmin, zmin, dt, dx, dz, lvect, l_particles_weight, l4symtry, l_2drz,               &
  type_rz_depose)      !#do not wrap
  USE picsar_precision, ONLY: idp, isp, lp, num
  implicit none
  ! __ Parameter declaration _______________________________________________
  integer(idp)                          :: np, type_rz_depose
  integer(idp)                          :: lvect
  INTEGER(idp), intent(in)              :: jx_nguard(2), jx_nvalid(2), jy_nguard(2),  &
  jy_nvalid(2), jz_nguard(2), jz_nvalid(2)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1 )
  real(num), dimension(np)              :: xp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp)                           :: l_particles_weight, l4symtry, l_2drz
  real(num)                             :: dxi, dzi, xint, zint
  real(num)                             :: xold, zold, x, z, wq, wqx, wqz
  real(num)                             :: vx, vy, vz
  real(num)                             :: invvol, invdtdx, invdtdz
  real(num)                             :: oxint, ozint, xintsq, zintsq, oxintsq,     &
  ozintsq
  real(num)                             :: dtsdx0, dtsdz0
  real(num), parameter                  :: onesixth=1./6., twothird=2./3.
  integer(idp)                          :: iixp0, ikxp0, iixp, ikxp, ip, i, k, ic, kc
  integer(isp)                          :: dix, diz
  integer(isp)                          :: ixmin, ixmax, izmin, izmax

  real(num), dimension(6) :: sx(-2:3), sx0(-2:3), dsx(-2:3)
  real(num), dimension(6) :: sz(-2:3), sz0(-2:3), dsz(-2:3)
  real(num)               :: sdxi, sdxim1
  real(num), dimension(6) :: sdzk(-2:3), sdzkm1(-2:3)

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi

  sdxi   = 0.0_num
  sdxim1 = 0.0_num
  sdzk   = 0.0_num
  sdzkm1 = 0.0_num
!$acc parallel deviceptr(jx, jy, jz, xp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(-2:3), sz(-2:3), &
!$acc&                         sx0(-2:3), sz0(-2:3), &
!$acc&                         dsx(-2:3), dsz(-2:3), &
!$acc&                         sdxi, sdxim1, &
!$acc&                         sdzk(-2:3), sdzkm1(-2:3) )
  DO ip=1, np
     sx  = 0.0_num
     sz  = 0.0_num
     sx0 = 0.0_num
     sz0 = 0.0_num

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi

    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz

    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz

    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ikxp0=floor(z)

    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
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

    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ikxp=floor(zold)

    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp

    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0

    ! --- zero out coefficients
    ! --- (needed because of different dix and diz for each particle)
    sx(-2)=0.0_num
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sx(3)=0.0_num
    sz(-2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
    sz(3)=0.0_num

    ! --- computes coefficients for quantities centered between nodes
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1+dix) = onesixth*oxintsq*oxint
    sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx( 2+dix) = onesixth*xintsq*xint

    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1+diz) = onesixth*ozintsq*ozint
    sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz( 2+diz) = onesixth*zintsq*zint

    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0

    ! --- computes min/max positions of current contributions
    ixmin = min(0, dix)-1
    ixmax = max(0, dix)+2
    izmin = min(0, diz)-1
    izmax = max(0, diz)+2

    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k

        ! --- Jx
        IF(i<ixmax) THEN
          sdxi = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) ) ! Wx coefficient from esirkepov
          if (i>ixmin) sdxi = sdxi + sdxim1 ! Integration of Wx along x
          !$acc atomic update
          jx(ic, kc) = jx(ic, kc) + sdxi ! Deposition on the current
        END IF

        ! -- Jy (2D Esirkepov scheme)
        !$acc atomic update
        jy(ic, kc) = jy(ic, kc) + wq*vy*invvol* ( (sz0(k)+0.5*dsz(k))*sx0(i) +        &
        (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )

        ! --- Jz
        IF(k<izmax) THEN
          sdzk(i) = wqz*dsz(k)*(sx0(i)+0.5*dsx(i)) ! Wz coefficient from esirkepov
          if (k>izmin) sdzk(i) = sdzk(i) + sdzkm1(i) ! Integration of Wz along z
          !$acc atomic update
          jz(ic, kc) = jz(ic, kc) + sdzk(i) ! Deposition on the current
        END IF
        sdxim1 = sdxi
      END DO
      sdzkm1 = sdzk
    END DO

  END DO
!$acc end loop
!$acc end parallel
  RETURN

END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_3_3


#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 3
!> This function is semi-vectorized: only the first part
!> with the computation of indexes is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx, jy, jz current arrays
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx, nz tile grid size
!> @param[in] nxguard, nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] lvect vector length for vectorization
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________
subroutine pxr_depose_jxjyjz_esirkepov2d_svec_3_3(jx, jy, jz, np, xp, zp, uxp, uyp,   &
  uzp, gaminv, w, q, xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz, lvect, &
  l_particles_weight, l4symtry, l_2drz, type_rz_depose)
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, lp, num
  implicit none
  ! __ Parameter declaration _______________________________________________
  integer(idp)                          :: np, nx, nz, nox, noz, nxguard, nzguard,    &
  type_rz_depose
  integer(idp)                          :: lvect
  real(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard), intent(in out) ::   &
  jx, jy, jz
  real(num), dimension(np)              :: xp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp)                           :: l_particles_weight, l4symtry, l_2drz
  real(num)                             :: dxi, dzi, dtsdx, dtsdz, xint, zint
  real(num), dimension(:, :, :), allocatable :: sdx, sdz, sdy
  real(num)                             :: xold, zold, rold, xmid, zmid, x, z, c, s,  &
  wq, wqx, wqz
  real(num)                             :: tmp, vx, vy, vz, dts2dx, dts2dz
  real(num)                             :: invvol, invdtdx, invdtdz
  real(num)                             :: oxint, ozint, xintsq, zintsq, oxintsq,     &
  ozintsq
  real(num)                             :: dtsdx0, dtsdz0, dts2dx0, dts2dz0
  real(num), parameter                  :: onesixth=1./6., twothird=2./3.
  real(num), parameter                  :: onethird=1./3.
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(isp)                          :: iixp, ikxp, ip, dix, diz, idx, idz, i, k,  &
  ic, kc
  integer(isp)                          :: icell, ndtodx, ndtodz
  integer(isp), dimension(lvect)        :: ixmin, ixmax, izmin, izmax, iixp0, ikxp0
  integer(isp)                          :: xl, xu, zl, zu, n, nn

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(lvect, -2:3, -2:3), sdz(lvect, -2:3, -2:3))
  ALLOCATE(sx(-2:3), sx0(-2:3), dsx(-2:3))
  ALLOCATE(sz(-2:3), sz0(-2:3), dsz(-2:3))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num

  ! Outer loop on particles with period LVEC
  DO ip=1, np, LVEC

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, zp:64, gaminv:64
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, zp)
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, ICELL)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif

    ! Inner loop on particle
    DO n=1, MIN(LVECT, np-ip+1)
      nn=ip+n-1

      ! --- computes current position in grid units
      x = (xp(nn)-xmin)*dxi
      z = (zp(nn)-zmin)*dzi

      ! --- computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)

      ! --- computes old position in grid units
      xold=x-dtsdx0*vx
      zold=z-dtsdz0*vz

      ! --- computes particles weights
      wq=q*w(nn)
      wqx = wq*invdtdx
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      iixp0(n)=floor(x)
      ikxp0(n)=floor(z)

      ! --- computes distance between particle and node for current positions
      xint=x-iixp0(n)
      zint=z-ikxp0(n)

      ! --- computes coefficients for node centered quantities
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

      ! --- finds node of cell containing particles for old positions
      iixp=floor(xold)
      ikxp=floor(zold)

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0(n)
      diz = ikxp-ikxp0(n)

      ! --- zero out coefficients
      ! --- (needed because of different dix and diz for each particle)
      sx(-2)=0.0_num
      sx(-1)=0.0_num
      sx(0)=0.0_num
      sx(1)=0.0_num
      sx(2)=0.0_num
      sx(3)=0.0_num
      sz(-2)=0.0_num
      sz(-1)=0.0_num
      sz(0)=0.0_num
      sz(1)=0.0_num
      sz(2)=0.0_num
      sz(3)=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1+dix) = onesixth*oxintsq*oxint
      sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2+dix) = onesixth*xintsq*xint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1+diz) = onesixth*ozintsq*ozint
      sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2+diz) = onesixth*zintsq*zint

      ! --- computes coefficients difference
      dsx(-2) = sx(-2) - sx0(-2)
      dsz(-2) = sz(-2) - sz0(-2)
      dsx(-1) = sx(-1) - sx0(-1)
      dsz(-1) = sz(-1) - sz0(-1)
      dsx(0) = sx(0) - sx0(0)
      dsz(0) = sz(0) - sz0(0)
      dsx(1) = sx(1) - sx0(1)
      dsz(1) = sz(1) - sz0(1)
      dsx(2) = sx(2) - sx0(2)
      dsz(2) = sz(2) - sz0(2)
      dsx(3) = sx(3) - sx0(3)
      dsz(3) = sz(3) - sz0(3)

      ! --- computes min/max positions of current contributions
      ixmin(n) = min(0, dix)-1
      ixmax(n) = max(0, dix)+2
      izmin(n) = min(0, diz)-1
      izmax(n) = max(0, diz)+2

      sdx(n, -2, -2)  = wqx*dsx(-2)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, -1, -2)  = wqx*dsx(-1)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, -1, -2)=sdx(n, -1, -2)+sdx(n, -1-1, -2)
      sdx(n, 0, -2)  = wqx*dsx(0)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 0, -2)=sdx(n, 0, -2)+sdx(n, 0-1, -2)
      sdx(n, 1, -2)  = wqx*dsx(1)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 1, -2)=sdx(n, 1, -2)+sdx(n, 1-1, -2)
      sdx(n, 2, -2)  = wqx*dsx(2)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 2, -2)=sdx(n, 2, -2)+sdx(n, 2-1, -2)
      sdx(n, -2, -1)  = wqx*dsx(-2)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, -1, -1)  = wqx*dsx(-1)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, -1, -1)=sdx(n, -1, -1)+sdx(n, -1-1, -1)
      sdx(n, 0, -1)  = wqx*dsx(0)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 0, -1)=sdx(n, 0, -1)+sdx(n, 0-1, -1)
      sdx(n, 1, -1)  = wqx*dsx(1)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 1, -1)=sdx(n, 1, -1)+sdx(n, 1-1, -1)
      sdx(n, 2, -1)  = wqx*dsx(2)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 2, -1)=sdx(n, 2, -1)+sdx(n, 2-1, -1)
      sdx(n, -2, 0)  = wqx*dsx(-2)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, -1, 0)  = wqx*dsx(-1)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, -1, 0)=sdx(n, -1, 0)+sdx(n, -1-1, 0)
      sdx(n, 0, 0)  = wqx*dsx(0)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 0, 0)=sdx(n, 0, 0)+sdx(n, 0-1, 0)
      sdx(n, 1, 0)  = wqx*dsx(1)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 1, 0)=sdx(n, 1, 0)+sdx(n, 1-1, 0)
      sdx(n, 2, 0)  = wqx*dsx(2)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 2, 0)=sdx(n, 2, 0)+sdx(n, 2-1, 0)
      sdx(n, -2, 1)  = wqx*dsx(-2)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, -1, 1)  = wqx*dsx(-1)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, -1, 1)=sdx(n, -1, 1)+sdx(n, -1-1, 1)
      sdx(n, 0, 1)  = wqx*dsx(0)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 0, 1)=sdx(n, 0, 1)+sdx(n, 0-1, 1)
      sdx(n, 1, 1)  = wqx*dsx(1)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 1, 1)=sdx(n, 1, 1)+sdx(n, 1-1, 1)
      sdx(n, 2, 1)  = wqx*dsx(2)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 2, 1)=sdx(n, 2, 1)+sdx(n, 2-1, 1)
      sdx(n, -2, 2)  = wqx*dsx(-2)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, -1, 2)  = wqx*dsx(-1)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, -1, 2)=sdx(n, -1, 2)+sdx(n, -1-1, 2)
      sdx(n, 0, 2)  = wqx*dsx(0)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 0, 2)=sdx(n, 0, 2)+sdx(n, 0-1, 2)
      sdx(n, 1, 2)  = wqx*dsx(1)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 1, 2)=sdx(n, 1, 2)+sdx(n, 1-1, 2)
      sdx(n, 2, 2)  = wqx*dsx(2)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 2, 2)=sdx(n, 2, 2)+sdx(n, 2-1, 2)
      sdx(n, -2, 3)  = wqx*dsx(-2)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, -1, 3)  = wqx*dsx(-1)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, -1, 3)=sdx(n, -1, 3)+sdx(n, -1-1, 3)
      sdx(n, 0, 3)  = wqx*dsx(0)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 0, 3)=sdx(n, 0, 3)+sdx(n, 0-1, 3)
      sdx(n, 1, 3)  = wqx*dsx(1)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 1, 3)=sdx(n, 1, 3)+sdx(n, 1-1, 3)
      sdx(n, 2, 3)  = wqx*dsx(2)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 2, 3)=sdx(n, 2, 3)+sdx(n, 2-1, 3)

      sdy(n, -2, -2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(-2) +                &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-2) )
      sdy(n, -1, -2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(-1) +                &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-1) )
      sdy(n, 0, -2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(0) +                  &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(0) )
      sdy(n, 1, -2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(1) +                  &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(1) )
      sdy(n, 2, -2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(2) +                  &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(2) )
      sdy(n, 3, -2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(3) +                  &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(3) )
      sdy(n, -2, -1) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(-2) +                &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-2) )
      sdy(n, -1, -1) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(-1) +                &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-1) )
      sdy(n, 0, -1) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(0) +                  &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(0) )
      sdy(n, 1, -1) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(1) +                  &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(1) )
      sdy(n, 2, -1) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(2) +                  &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(2) )
      sdy(n, 3, -1) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(3) +                  &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(3) )
      sdy(n, -2, 0) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(-2) +                   &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(-2) )
      sdy(n, -1, 0) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(-1) +                   &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(-1) )
      sdy(n, 0, 0) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(0) +                     &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(0) )
      sdy(n, 1, 0) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(1) +                     &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(1) )
      sdy(n, 2, 0) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(2) +                     &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(2) )
      sdy(n, 3, 0) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(3) +                     &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(3) )
      sdy(n, -2, 1) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(-2) +                   &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(-2) )
      sdy(n, -1, 1) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(-1) +                   &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(-1) )
      sdy(n, 0, 1) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(0) +                     &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(0) )
      sdy(n, 1, 1) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(1) +                     &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(1) )
      sdy(n, 2, 1) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(2) +                     &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(2) )
      sdy(n, 3, 1) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(3) +                     &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(3) )
      sdy(n, -2, 2) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(-2) +                   &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(-2) )
      sdy(n, -1, 2) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(-1) +                   &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(-1) )
      sdy(n, 0, 2) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(0) +                     &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(0) )
      sdy(n, 1, 2) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(1) +                     &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(1) )
      sdy(n, 2, 2) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(2) +                     &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(2) )
      sdy(n, 3, 2) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(3) +                     &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(3) )
      sdy(n, -2, 3) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(-2) +                   &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(-2) )
      sdy(n, -1, 3) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(-1) +                   &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(-1) )
      sdy(n, 0, 3) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(0) +                     &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(0) )
      sdy(n, 1, 3) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(1) +                     &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(1) )
      sdy(n, 2, 3) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(2) +                     &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(2) )
      sdy(n, 3, 3) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(3) +                     &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(3) )

      sdz(n, -2, -2)=wqz*dsz(-2)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, -1, -2)=wqz*dsz(-2)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, 0, -2)=wqz*dsz(-2)*(sx0(0)+0.5*dsx(0))
      sdz(n, 1, -2)=wqz*dsz(-2)*(sx0(1)+0.5*dsx(1))
      sdz(n, 2, -2)=wqz*dsz(-2)*(sx0(2)+0.5*dsx(2))
      sdz(n, 3, -2)=wqz*dsz(-2)*(sx0(3)+0.5*dsx(3))
      sdz(n, -2, -1)=wqz*dsz(-1)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, -2, -1)=sdz(n, -2, -1)+sdz(n, -2, -1-1)
      sdz(n, -1, -1)=wqz*dsz(-1)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, -1, -1)=sdz(n, -1, -1)+sdz(n, -1, -1-1)
      sdz(n, 0, -1)=wqz*dsz(-1)*(sx0(0)+0.5*dsx(0))
      sdz(n, 0, -1)=sdz(n, 0, -1)+sdz(n, 0, -1-1)
      sdz(n, 1, -1)=wqz*dsz(-1)*(sx0(1)+0.5*dsx(1))
      sdz(n, 1, -1)=sdz(n, 1, -1)+sdz(n, 1, -1-1)
      sdz(n, 2, -1)=wqz*dsz(-1)*(sx0(2)+0.5*dsx(2))
      sdz(n, 2, -1)=sdz(n, 2, -1)+sdz(n, 2, -1-1)
      sdz(n, 3, -1)=wqz*dsz(-1)*(sx0(3)+0.5*dsx(3))
      sdz(n, 3, -1)=sdz(n, 3, -1)+sdz(n, 3, -1-1)
      sdz(n, -2, 0)=wqz*dsz(0)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, -2, 0)=sdz(n, -2, 0)+sdz(n, -2, 0-1)
      sdz(n, -1, 0)=wqz*dsz(0)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, -1, 0)=sdz(n, -1, 0)+sdz(n, -1, 0-1)
      sdz(n, 0, 0)=wqz*dsz(0)*(sx0(0)+0.5*dsx(0))
      sdz(n, 0, 0)=sdz(n, 0, 0)+sdz(n, 0, 0-1)
      sdz(n, 1, 0)=wqz*dsz(0)*(sx0(1)+0.5*dsx(1))
      sdz(n, 1, 0)=sdz(n, 1, 0)+sdz(n, 1, 0-1)
      sdz(n, 2, 0)=wqz*dsz(0)*(sx0(2)+0.5*dsx(2))
      sdz(n, 2, 0)=sdz(n, 2, 0)+sdz(n, 2, 0-1)
      sdz(n, 3, 0)=wqz*dsz(0)*(sx0(3)+0.5*dsx(3))
      sdz(n, 3, 0)=sdz(n, 3, 0)+sdz(n, 3, 0-1)
      sdz(n, -2, 1)=wqz*dsz(1)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, -2, 1)=sdz(n, -2, 1)+sdz(n, -2, 1-1)
      sdz(n, -1, 1)=wqz*dsz(1)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, -1, 1)=sdz(n, -1, 1)+sdz(n, -1, 1-1)
      sdz(n, 0, 1)=wqz*dsz(1)*(sx0(0)+0.5*dsx(0))
      sdz(n, 0, 1)=sdz(n, 0, 1)+sdz(n, 0, 1-1)
      sdz(n, 1, 1)=wqz*dsz(1)*(sx0(1)+0.5*dsx(1))
      sdz(n, 1, 1)=sdz(n, 1, 1)+sdz(n, 1, 1-1)
      sdz(n, 2, 1)=wqz*dsz(1)*(sx0(2)+0.5*dsx(2))
      sdz(n, 2, 1)=sdz(n, 2, 1)+sdz(n, 2, 1-1)
      sdz(n, 3, 1)=wqz*dsz(1)*(sx0(3)+0.5*dsx(3))
      sdz(n, 3, 1)=sdz(n, 3, 1)+sdz(n, 3, 1-1)
      sdz(n, -2, 2)=wqz*dsz(2)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, -2, 2)=sdz(n, -2, 2)+sdz(n, -2, 2-1)
      sdz(n, -1, 2)=wqz*dsz(2)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, -1, 2)=sdz(n, -1, 2)+sdz(n, -1, 2-1)
      sdz(n, 0, 2)=wqz*dsz(2)*(sx0(0)+0.5*dsx(0))
      sdz(n, 0, 2)=sdz(n, 0, 2)+sdz(n, 0, 2-1)
      sdz(n, 1, 2)=wqz*dsz(2)*(sx0(1)+0.5*dsx(1))
      sdz(n, 1, 2)=sdz(n, 1, 2)+sdz(n, 1, 2-1)
      sdz(n, 2, 2)=wqz*dsz(2)*(sx0(2)+0.5*dsx(2))
      sdz(n, 2, 2)=sdz(n, 2, 2)+sdz(n, 2, 2-1)
      sdz(n, 3, 2)=wqz*dsz(2)*(sx0(3)+0.5*dsx(3))
      sdz(n, 3, 2)=sdz(n, 3, 2)+sdz(n, 3, 2-1)
    ENDDO

    ! Inner loop on particle
    DO n=1, MIN(LVECT, np-ip+1)
      ! --- add current contributions
      DO k=izmin(n), izmax(n)
        DO i=ixmin(n), ixmax(n)
          ic = iixp0(n)+i
          kc = ikxp0(n)+k
          ! --- Jx
          jx(ic, kc) = jx(ic, kc) + sdx(n, i, k)! Deposition on the current
          ! -- Jy (2D Esirkepov scheme)
          jy(ic, kc) = jy(ic, kc) + sdy(n, i, k)
          ! --- Jz
          jz(ic, kc) = jz(ic, kc) + sdz(n, i, k)! Deposition on the current
        END DO
      END DO
    ENDDO
  END DO
  DEALLOCATE(sdx, sdz, sx, sx0, dsx, sz, sz0, dsz)
END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_svec_3_3
#endif

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Vectorized 2D Current deposition with the method of Esirkepov at order 3.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx, jy, jz current arrays
!> @param[in] np number of particles
!> @param[in] xp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx, nz tile grid size
!> @param[in] nxguard, nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] lvect vector length for vectorization
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________
subroutine pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3(jx, jy, jz, np, xp, zp, uxp, uyp,  &
  uzp, gaminv, w, q, xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz, lvect, &
  l_particles_weight, l4symtry, l_2drz, type_rz_depose)   !#do not parse
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, lp, num
  implicit none
  ! __ Parameter declaration ____________________________________________________
  ! In/out parameters
  integer(idp)                          :: np, nx, nz, nox, noz, nxguard, nzguard,    &
  type_rz_depose
  integer(idp)                          :: lvect
  real(num), dimension((1+nx+2*nxguard)*(1+nz+2*nzguard)), intent(in out) :: jx, jy,  &
  jz
  real(num), dimension(np)              :: xp, zp, uxp, uyp, uzp, gaminv, w
  real(num)                             :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp)                           :: l_particles_weight, l4symtry, l_2drz
  ! Local parameters
  real(num)                             :: dxi, dzi, dtsdx, dtsdz, xint, zint
  real(num)                             :: xold, zold, rold, xmid, zmid, x, z, c, s,  &
  wq, wqx, wqz
  real(num)                             :: tmp, vx, vy, vz, dts2dx, dts2dz
  real(num)                             :: invvol, invdtdx, invdtdz
  real(num)                             :: oxint, ozint, xintsq, zintsq, oxintsq,     &
  ozintsq
  real(num)                             :: dtsdx0, dtsdz0, dts2dx0, dts2dz0
  real(num), parameter                  :: onesixth=1./6., twothird=2./3.
  real(num), parameter                  :: onethird=1./3.
  real(num), dimension(-2:3)            :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0, ikxp0, iixp, ikxp, ip, dix, diz,    &
  idx, idz, i, k, ic, kc
  integer(idp)                          :: ixmin, ixmax, izmin, izmax, ndtodx, ndtodz
  integer(idp)                          :: xl, xu, zl, zu
  integer(isp)                          :: ngridx, ncx, ncz, ncells
  integer(isp)                          :: n, nn, nv
  integer(isp)                          :: iixporig, ikxporig, igrid, ix, iz, orig,   &
  nnx
  INTEGER(isp), DIMENSION(LVEC, 3)       :: ICELL
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: jxcells, jycells, jzcells
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: sdx, sdy, sdz

  ! __ Parameter initialization __________________________________________________

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi

  sx0(-2:3)=0._num
  sz0(-2:3)=0._num

  ngridx=nx+1+2*nxguard
  ncx=nx+1+2*nxguard
  ncz=nz+1+2*nzguard
  NCELLS=ncx*ncz
  iixporig=-nxguard
  ikxporig=-nzguard
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))


  jxcells = 0._num
  jycells = 0._num
  jzcells = 0._num

  ALLOCATE(sdx(lvect, 48), sdy(lvect, 48), sdz(lvect, 48))

  sdx = 0._num
  sdz = 0._num

  nnx = ngridx
  orig=(nxguard+iixporig) + (nzguard+ikxporig)*nnx

  ! Outer loop on particles with period LVEC
  DO ip=1, np, LVEC

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, zp:64, gaminv:64
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, ICELL)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif

    ! Inner loop on particle
    DO n=1, MIN(LVECT, np-ip+1)
      nn=ip+n-1

      ! --- computes current position in grid units
      x = (xp(nn)-xmin)*dxi
      z = (zp(nn)-zmin)*dzi

      ! --- computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)

      ! --- computes old position in grid units
      xold=x-dtsdx0*vx
      zold=z-dtsdz0*vz

      ! --- computes particles weights
      wq=q*w(ip)
      wqx = wq*invdtdx
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      iixp0=floor(x)
      ikxp0=floor(z)

      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
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

      ! --- finds node of cell containing particles for old positions
      iixp=floor(xold)
      ikxp=floor(zold)

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients (needed because of different dix and diz
      ! --- for each particle)
      sx(-2)=0.0_num
      sx(-1)=0.0_num
      sx(0)=0.0_num
      sx(1)=0.0_num
      sx(2)=0.0_num
      sx(3)=0.0_num
      sz(-2)=0.0_num
      sz(-1)=0.0_num
      sz(0)=0.0_num
      sz(1)=0.0_num
      sz(2)=0.0_num
      sz(3)=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1+dix) = onesixth*oxintsq*oxint
      sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2+dix) = onesixth*xintsq*xint

      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1+diz) = onesixth*ozintsq*ozint
      sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2+diz) = onesixth*zintsq*zint

      ! --- computes coefficients difference
      dsx(-2) = sx(-2) - sx0(-2)
      dsz(-2) = sz(-2) - sz0(-2)
      dsx(-1) = sx(-1) - sx0(-1)
      dsz(-1) = sz(-1) - sz0(-1)
      dsx(0) = sx(0) - sx0(0)
      dsz(0) = sz(0) - sz0(0)
      dsx(1) = sx(1) - sx0(1)
      dsz(1) = sz(1) - sz0(1)
      dsx(2) = sx(2) - sx0(2)
      dsz(2) = sz(2) - sz0(2)
      dsx(3) = sx(3) - sx0(3)
      dsz(3) = sz(3) - sz0(3)

      ! --- Position of the first cell (-2, -2, -2)
      icell(n, 1) = (iixp0-iixporig-1)+(ikxp0-ikxporig-2)*ncx

      ! --- Weight
      sdx(n, 1)  = wqx*dsx(-2)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 2)  = wqx*dsx(-1)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 2) = sdx(n, 2)+sdx(n, 1)
      sdx(n, 3)  = wqx*dsx(0)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 3) = sdx(n, 3)+sdx(n, 2)
      sdx(n, 4)  = wqx*dsx(1)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 4) = sdx(n, 4)+sdx(n, 3)
      sdx(n, 5)  = wqx*dsx(2)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n, 5) = sdx(n, 5)+sdx(n, 4)
      sdx(n, 9)  = wqx*dsx(-2)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 10)  = wqx*dsx(-1)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 10) = sdx(n, 10)+sdx(n, 9)
      sdx(n, 11)  = wqx*dsx(0)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 11) = sdx(n, 11)+sdx(n, 10)
      sdx(n, 12)  = wqx*dsx(1)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 12) = sdx(n, 12)+sdx(n, 11)
      sdx(n, 13)  = wqx*dsx(2)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n, 13) = sdx(n, 13)+sdx(n, 12)
      sdx(n, 17)  = wqx*dsx(-2)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 18)  = wqx*dsx(-1)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 18) = sdx(n, 18)+sdx(n, 17)
      sdx(n, 19)  = wqx*dsx(0)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 19) = sdx(n, 19)+sdx(n, 18)
      sdx(n, 20)  = wqx*dsx(1)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 20) = sdx(n, 20)+sdx(n, 19)
      sdx(n, 21)  = wqx*dsx(2)*( sz0(0) + 0.5*dsz(0) )
      sdx(n, 21) = sdx(n, 21)+sdx(n, 20)
      sdx(n, 25)  = wqx*dsx(-2)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 26)  = wqx*dsx(-1)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 26) = sdx(n, 26)+sdx(n, 25)
      sdx(n, 27)  = wqx*dsx(0)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 27) = sdx(n, 27)+sdx(n, 26)
      sdx(n, 28)  = wqx*dsx(1)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 28) = sdx(n, 28)+sdx(n, 27)
      sdx(n, 29)  = wqx*dsx(2)*( sz0(1) + 0.5*dsz(1) )
      sdx(n, 29) = sdx(n, 29)+sdx(n, 28)
      sdx(n, 33)  = wqx*dsx(-2)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 34)  = wqx*dsx(-1)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 34) = sdx(n, 34)+sdx(n, 33)
      sdx(n, 35)  = wqx*dsx(0)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 35) = sdx(n, 35)+sdx(n, 34)
      sdx(n, 36)  = wqx*dsx(1)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 36) = sdx(n, 36)+sdx(n, 35)
      sdx(n, 37)  = wqx*dsx(2)*( sz0(2) + 0.5*dsz(2) )
      sdx(n, 37) = sdx(n, 37)+sdx(n, 36)
      sdx(n, 41)  = wqx*dsx(-2)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 42)  = wqx*dsx(-1)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 42) = sdx(n, 42)+sdx(n, 41)
      sdx(n, 43)  = wqx*dsx(0)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 43) = sdx(n, 43)+sdx(n, 42)
      sdx(n, 44)  = wqx*dsx(1)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 44) = sdx(n, 44)+sdx(n, 43)
      sdx(n, 45)  = wqx*dsx(2)*( sz0(3) + 0.5*dsz(3) )
      sdx(n, 45) = sdx(n, 45)+sdx(n, 44)

      sdy(n, 1) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(-2) +                     &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-2))
      sdy(n, 2) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(-1) +                     &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-1))
      sdy(n, 3) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(0) +                      &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(0))
      sdy(n, 4) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(1) +                      &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(1))
      sdy(n, 5) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(2) +                      &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(2))
      sdy(n, 6) = wq*vy*invvol* ( (sz0(-2)+0.5*dsz(-2))*sx0(3) +                      &
      (0.5*sz0(-2)+onethird*dsz(-2))*dsx(3))
      sdy(n, 9) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(-2) +                     &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-2))
      sdy(n, 10) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(-1) +                    &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 11) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(0) +                     &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 12) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(1) +                     &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 13) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(2) +                     &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 14) = wq*vy*invvol* ( (sz0(-1)+0.5*dsz(-1))*sx0(3) +                     &
      (0.5*sz0(-1)+onethird*dsz(-1))*dsx(3))
      sdy(n, 17) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(-2) +                      &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(-2))
      sdy(n, 18) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(-1) +                      &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 19) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(0) +                       &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 20) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(1) +                       &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 21) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(2) +                       &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 22) = wq*vy*invvol* ( (sz0(0)+0.5*dsz(0))*sx0(3) +                       &
      (0.5*sz0(0)+onethird*dsz(0))*dsx(3))
      sdy(n, 25) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(-2) +                      &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(-2))
      sdy(n, 26) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(-1) +                      &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 27) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(0) +                       &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 28) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(1) +                       &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 29) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(2) +                       &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 30) = wq*vy*invvol* ( (sz0(1)+0.5*dsz(1))*sx0(3) +                       &
      (0.5*sz0(1)+onethird*dsz(1))*dsx(3))
      sdy(n, 33) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(-2) +                      &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(-2))
      sdy(n, 34) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(-1) +                      &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 35) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(0) +                       &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 36) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(1) +                       &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 37) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(2) +                       &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 38) = wq*vy*invvol* ( (sz0(2)+0.5*dsz(2))*sx0(3) +                       &
      (0.5*sz0(2)+onethird*dsz(2))*dsx(3))
      sdy(n, 41) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(-2) +                      &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(-2))
      sdy(n, 42) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(-1) +                      &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(-1))
      sdy(n, 43) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(0) +                       &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(0))
      sdy(n, 44) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(1) +                       &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(1))
      sdy(n, 45) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(2) +                       &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(2))
      sdy(n, 46) = wq*vy*invvol* ( (sz0(3)+0.5*dsz(3))*sx0(3) +                       &
      (0.5*sz0(3)+onethird*dsz(3))*dsx(3))

      sdz(n, 1)  = wqz*dsz(-2)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, 2)  = wqz*dsz(-2)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, 3)  = wqz*dsz(-2)*(sx0(0)+0.5*dsx(0))
      sdz(n, 4)  = wqz*dsz(-2)*(sx0(1)+0.5*dsx(1))
      sdz(n, 5)  = wqz*dsz(-2)*(sx0(2)+0.5*dsx(2))
      sdz(n, 6)  = wqz*dsz(-2)*(sx0(3)+0.5*dsx(3))
      sdz(n, 9)  = wqz*dsz(-1)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, 9) = sdz(n, 9)+sdz(n, 8)
      sdz(n, 10)  = wqz*dsz(-1)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, 10) = sdz(n, 10)+sdz(n, 9)
      sdz(n, 11)  = wqz*dsz(-1)*(sx0(0)+0.5*dsx(0))
      sdz(n, 11) = sdz(n, 11)+sdz(n, 10)
      sdz(n, 12)  = wqz*dsz(-1)*(sx0(1)+0.5*dsx(1))
      sdz(n, 12) = sdz(n, 12)+sdz(n, 11)
      sdz(n, 13)  = wqz*dsz(-1)*(sx0(2)+0.5*dsx(2))
      sdz(n, 13) = sdz(n, 13)+sdz(n, 12)
      sdz(n, 14)  = wqz*dsz(-1)*(sx0(3)+0.5*dsx(3))
      sdz(n, 14) = sdz(n, 14)+sdz(n, 13)
      sdz(n, 17)  = wqz*dsz(0)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, 17) = sdz(n, 17)+sdz(n, 16)
      sdz(n, 18)  = wqz*dsz(0)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, 18) = sdz(n, 18)+sdz(n, 17)
      sdz(n, 19)  = wqz*dsz(0)*(sx0(0)+0.5*dsx(0))
      sdz(n, 19) = sdz(n, 19)+sdz(n, 18)
      sdz(n, 20)  = wqz*dsz(0)*(sx0(1)+0.5*dsx(1))
      sdz(n, 20) = sdz(n, 20)+sdz(n, 19)
      sdz(n, 21)  = wqz*dsz(0)*(sx0(2)+0.5*dsx(2))
      sdz(n, 21) = sdz(n, 21)+sdz(n, 20)
      sdz(n, 22)  = wqz*dsz(0)*(sx0(3)+0.5*dsx(3))
      sdz(n, 22) = sdz(n, 22)+sdz(n, 21)
      sdz(n, 25)  = wqz*dsz(1)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, 25) = sdz(n, 25)+sdz(n, 24)
      sdz(n, 26)  = wqz*dsz(1)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, 26) = sdz(n, 26)+sdz(n, 25)
      sdz(n, 27)  = wqz*dsz(1)*(sx0(0)+0.5*dsx(0))
      sdz(n, 27) = sdz(n, 27)+sdz(n, 26)
      sdz(n, 28)  = wqz*dsz(1)*(sx0(1)+0.5*dsx(1))
      sdz(n, 28) = sdz(n, 28)+sdz(n, 27)
      sdz(n, 29)  = wqz*dsz(1)*(sx0(2)+0.5*dsx(2))
      sdz(n, 29) = sdz(n, 29)+sdz(n, 28)
      sdz(n, 30)  = wqz*dsz(1)*(sx0(3)+0.5*dsx(3))
      sdz(n, 30) = sdz(n, 30)+sdz(n, 29)
      sdz(n, 33)  = wqz*dsz(2)*(sx0(-2)+0.5*dsx(-2))
      sdz(n, 33) = sdz(n, 33)+sdz(n, 32)
      sdz(n, 34)  = wqz*dsz(2)*(sx0(-1)+0.5*dsx(-1))
      sdz(n, 34) = sdz(n, 34)+sdz(n, 33)
      sdz(n, 35)  = wqz*dsz(2)*(sx0(0)+0.5*dsx(0))
      sdz(n, 35) = sdz(n, 35)+sdz(n, 34)
      sdz(n, 36)  = wqz*dsz(2)*(sx0(1)+0.5*dsx(1))
      sdz(n, 36) = sdz(n, 36)+sdz(n, 35)
      sdz(n, 37)  = wqz*dsz(2)*(sx0(2)+0.5*dsx(2))
      sdz(n, 37) = sdz(n, 37)+sdz(n, 36)
      sdz(n, 38)  = wqz*dsz(2)*(sx0(3)+0.5*dsx(3))
      sdz(n, 38) = sdz(n, 38)+sdz(n, 37)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
    !$OMP END SIMD
#endif

    ! Add weights to nearest vertices
    DO n=1, MIN(LVECT, np-ip+1)
      ! Alignment
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, jxcells, jycells, jzcells)
#endif
#endif

#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO nv=1, 8

        !print*, nv, ICELL(n, 1)

        jxcells(nv, ICELL(n, 1)) = jxcells(nv, ICELL(n, 1)) + sdx(n, nv)
        jxcells(nv, ICELL(n, 1)+ncx) = jxcells(nv, ICELL(n, 1)+ncx) + sdx(n, nv+8)
        jxcells(nv, ICELL(n, 1)+2*ncx) = jxcells(nv, ICELL(n, 1)+2*ncx) + sdx(n,      &
        nv+16)
        jxcells(nv, ICELL(n, 1)+3*ncx) = jxcells(nv, ICELL(n, 1)+3*ncx) + sdx(n,      &
        nv+24)
        jxcells(nv, ICELL(n, 1)+4*ncx) = jxcells(nv, ICELL(n, 1)+4*ncx) + sdx(n,      &
        nv+32)
        jxcells(nv, ICELL(n, 1)+5*ncx) = jxcells(nv, ICELL(n, 1)+5*ncx) + sdx(n,      &
        nv+40)

        jycells(nv, ICELL(n, 1)) = jycells(nv, ICELL(n, 1)) + sdy(n, nv)
        jycells(nv, ICELL(n, 1)+ncx) = jycells(nv, ICELL(n, 1)+ncx) + sdy(n, nv+8)
        jycells(nv, ICELL(n, 1)+2*ncx) = jycells(nv, ICELL(n, 1)+2*ncx) + sdy(n,      &
        nv+16)
        jycells(nv, ICELL(n, 1)+3*ncx) = jycells(nv, ICELL(n, 1)+3*ncx) + sdy(n,      &
        nv+24)
        jycells(nv, ICELL(n, 1)+4*ncx) = jycells(nv, ICELL(n, 1)+4*ncx) + sdy(n,      &
        nv+32)
        jycells(nv, ICELL(n, 1)+5*ncx) = jycells(nv, ICELL(n, 1)+5*ncx) + sdy(n,      &
        nv+40)

        jzcells(nv, ICELL(n, 1)) = jzcells(nv, ICELL(n, 1)) + sdz(n, nv)
        jzcells(nv, ICELL(n, 1)+ncx) = jzcells(nv, ICELL(n, 1)+ncx) + sdz(n, nv+8)
        jzcells(nv, ICELL(n, 1)+2*ncx) = jzcells(nv, ICELL(n, 1)+2*ncx) + sdz(n,      &
        nv+16)
        jzcells(nv, ICELL(n, 1)+3*ncx) = jzcells(nv, ICELL(n, 1)+3*ncx) + sdz(n,      &
        nv+24)
        jzcells(nv, ICELL(n, 1)+4*ncx) = jzcells(nv, ICELL(n, 1)+4*ncx) + sdz(n,      &
        nv+32)
      ENDDO
#if defined _OPENMP && _OPENMP>=201307
      !$OMP END SIMD
#endif
    ENDDO
  ENDDO

  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
  DO iz=1, ncz-2
#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !$DIR SIMD
#endif
    DO ix=1, ncx-8!! VECTOR (take ncx multiple of vector length)
      ic=ix+(iz-1)*ncx
      igrid =orig+ix+(iz-1)*nnx
      ! jx
      jx(igrid)=jx(igrid)+jxcells(1, ic)
      jx(igrid+1)=jx(igrid+1)+jxcells(2, ic)
      jx(igrid+2)=jx(igrid+2)+jxcells(3, ic)
      jx(igrid+3)=jx(igrid+3)+jxcells(4, ic)
      jx(igrid+4)=jx(igrid+4)+jxcells(5, ic)
      jx(igrid+5)=jx(igrid+5)+jxcells(6, ic)
      jx(igrid+6)=jx(igrid+6)+jxcells(7, ic)
      jx(igrid+7)=jx(igrid+7)+jxcells(8, ic)
      ! jy
      jy(igrid)  =jy(igrid)  +jycells(1, ic)
      jy(igrid+1)=jy(igrid+1)+jycells(2, ic)
      jy(igrid+2)=jy(igrid+2)+jycells(3, ic)
      jy(igrid+3)=jy(igrid+3)+jycells(4, ic)
      jy(igrid+4)=jy(igrid+4)+jycells(5, ic)
      jy(igrid+5)=jy(igrid+5)+jycells(6, ic)
      jy(igrid+6)=jy(igrid+6)+jycells(7, ic)
      jy(igrid+7)=jy(igrid+7)+jycells(8, ic)
      ! jz
      jz(igrid)  =jz(igrid)  +jzcells(1, ic)
      jz(igrid+1)=jz(igrid+1)+jzcells(2, ic)
      jz(igrid+2)=jz(igrid+2)+jzcells(3, ic)
      jz(igrid+3)=jz(igrid+3)+jzcells(4, ic)
      jz(igrid+4)=jz(igrid+4)+jzcells(5, ic)
      jz(igrid+5)=jz(igrid+5)+jzcells(6, ic)
      jz(igrid+6)=jz(igrid+6)+jzcells(7, ic)
      jz(igrid+7)=jz(igrid+7)+jzcells(8, ic)
    END DO
#if defined _OPENMP && _OPENMP>=201307
    !$OMP END SIMD
#endif
  ENDDO

  DEALLOCATE(jxcells, jycells, jzcells)
  DEALLOCATE(sdx, sdy, sdz)

END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3
#endif

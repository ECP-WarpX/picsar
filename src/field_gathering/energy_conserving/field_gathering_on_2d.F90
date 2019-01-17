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
! FIELD_GATHERING_MANAGER_2D.F90
!
! This file contains subroutines for the field gathering in 2D.
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
! - pxr_gete2dxz_n_energy_conserving
! - pxr_getb2dxz_n_energy_conserving
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> 2D electric field non-optimized gathering routine
!
!> @details
!> This function is similar to what is implemented in WARP
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nx, nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] exg, eyg, ezg field arrays
!> @param[in] lvect vector size for the block of particles
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,  &
  dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,      &
  ezg_nguard, ezg_nvalid, l4symtry, l_2drz, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp)             :: np, nox, noz
  integer(idp), intent(IN) :: exg_nguard(2), exg_nvalid(2), eyg_nguard(2),            &
  eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  real(num), dimension(np) :: xp, yp, zp, ex, ey, ez
  logical(lp)              :: l4symtry, l_2drz, l_lower_order_in_v, l_nodal
  real(num)                :: stagger_shift
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1, 1,        &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1, 1,        &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, 1,        &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  real(num)                :: xmin, zmin, dx, dz, costheta, sintheta
  integer(idp)             :: ip, j, l, ixmin, ixmax, izmin, izmax, ixmin0, ixmax0,   &
  izmin0, izmax0, jj, ll, j0, l0
  real(num) :: dxi, dzi, x, y, z, r, xint, zint, xintsq, oxint, zintsq, ozint,        &
  oxintsq, ozintsq, signx
  real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  real(num), dimension(:), allocatable :: sx0, sz0
  real(num), parameter :: onesixth=1./6., twothird=2./3.

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dzi = 1./dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)-1
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)-1

  if (l_lower_order_in_v) then
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
  else
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
  end if
  allocate(sx0(ixmin0:ixmax0), sz0(izmin0:izmax0))

  signx = 1.

  do ip=1, np

    if (l_2drz) then
      x = xp(ip)
      y = yp(ip)
      r=sqrt(x*x+y*y)
      if (r*dxi>1.e-20) then
        costheta=x/r
        sintheta=y/r
      else
        costheta=1.
        sintheta=0.
      end if
      x = (r-xmin)*dxi
    else
      x = (xp(ip)-xmin)*dxi
    end if

    z = (zp(ip)-zmin)*dzi

    if (l4symtry) then
      if (x<0.) then
        x = -x
        signx = -1.
      else
        signx = 1.
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
      if (noz==2*(noz/2)) then
        l=nint(z)
        l0=floor(z+0.5_num-stagger_shift)
      else
        l=floor(z)
        l0=floor(z-stagger_shift)
      end if
    end if

    xint=x-j
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

    if (l_2drz) then

      !          write(0, *) 'field gathering needs to be done for fstype=4 in EM-RZ'
      !          stop
      do ll = izmin, izmax+1
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sz(ll)*sx0(jj)*(exg(j0+jj, 1, l+ll)*costheta-eyg(j0+jj,   &
          1, l+ll)*sintheta)
          ey(ip) = ey(ip) + sz(ll)*sx0(jj)*(exg(j0+jj, 1, l+ll)*sintheta+eyg(j0+jj,   &
          1, l+ll)*costheta)
        end do
      end do

    else

      do ll = izmin, izmax+1
        do jj = ixmin0, ixmax0
          ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j0+jj, 1, l+ll)*signx
        end do
      end do

      do ll = izmin, izmax+1
        do jj = ixmin, ixmax+1
          ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj, 1, l+ll)
        end do
      end do

    end if

    do ll = izmin0, izmax0
      do jj = ixmin, ixmax+1
        ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj, 1, l0+ll)
      end do
    end do

  end do
  deallocate(sx0, sz0)

  return
end subroutine pxr_gete2dxz_n_energy_conserving

! ________________________________________________________________________________________
!> @brief
!> 2D magnetic field gathering routine for arbitrary orders
!
!> @details
!> This function is not vectorized and not optimized.
!
!> @author
!> Henri Vincenti
!
!> @date
!> 12/01/2016
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nx, nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] exg, eyg, ezg field arrays
!> @param[in] nox, noz interpolation order
!> @param[in] l4symtry
!> @param[in] l_2drz use the 2d cylindrical geometry system
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!
! ________________________________________________________________________________________
subroutine pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,  &
  dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,      &
  bzg_nguard, bzg_nvalid, l4symtry, l_2drz, l_lower_order_in_v, l_nodal)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp) :: np, nox, noz
  integer(idp), intent(IN)                :: bxg_nguard(2), bxg_nvalid(2),            &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  real(num), dimension(np) :: xp, yp, zp, bx, by, bz
  logical(lp)  :: l4symtry, l_2drz, l_lower_order_in_v, l_nodal
  real(num)    :: stagger_shift
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1, 1,        &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1, 1,        &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1, 1,        &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)
  real(num) :: xmin, zmin, dx, dz
  integer(idp) :: ip, j, l, ixmin, ixmax, izmin, izmax, ixmin0, ixmax0, izmin0,       &
  izmax0, jj, ll, j0, l0
  real(num) :: dxi, dzi, x, y, z, xint, zint, xintsq, oxint, zintsq, ozint, oxintsq,  &
  ozintsq, signx, r, costheta, sintheta
  real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  real(num), dimension(:), allocatable :: sx0, sz0
  real(num), parameter :: onesixth=1./6., twothird=2./3.

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1./dx
  dzi = 1./dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)-1
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)-1

  if (l_lower_order_in_v) then
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
  else
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
  end if
  allocate(sx0(ixmin0:ixmax0), sz0(izmin0:izmax0))

  signx = 1.

  sx=0
  sz=0.
  sx0=0.
  sz0=0.

  do ip=1, np

    if (l_2drz) then
      x = xp(ip)
      y = yp(ip)
      r=sqrt(x*x+y*y)
      if (r*dxi>1.e-20) then
        costheta=x/r
        sintheta=y/r
      else
        costheta=1.
        sintheta=0.
      end if
      x = (r-xmin)*dxi
    else
      x = (xp(ip)-xmin)*dxi
    end if

    z = (zp(ip)-zmin)*dzi

    if (l4symtry) then
      if (x<0.) then
        x = -x
        signx = -1.
      else
        signx = 1.
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
      if (noz==2*(noz/2)) then
        l=nint(z)
        l0=floor(z+0.5_num-stagger_shift)
      else
        l=floor(z)
        l0=floor(z-stagger_shift)
      end if
    end if

    xint=x-j
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

    if (l_2drz) then

      do ll = izmin0, izmax0
        do jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sz0(ll)*(bxg(j+jj, 1, l0+ll)*costheta-byg(j+jj, 1, &
          l0+ll)*sintheta)
          by(ip) = by(ip) + sx(jj)*sz0(ll)*(bxg(j+jj, 1, l0+ll)*sintheta+byg(j+jj, 1, &
          l0+ll)*costheta)
        end do
      end do

    else

      do ll = izmin0, izmax0
        do jj = ixmin, ixmax+1
          bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj, 1, l0+ll)*signx
        end do
      end do

      do ll = izmin0, izmax0
        do jj = ixmin0, ixmax0
          by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj, 1, l0+ll)
        end do
      end do

    end if

    do ll = izmin, izmax+1
      do jj = ixmin0, ixmax0
        bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j0+jj, 1, l+ll)
      end do
    end do

  end do
  deallocate(sx0, sz0)
  return
end subroutine pxr_getb2dxz_n_energy_conserving

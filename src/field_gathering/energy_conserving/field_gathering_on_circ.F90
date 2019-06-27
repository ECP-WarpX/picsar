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
! FIELD_GATHERING_MANAGER_circ.F90
!
! This file contains subroutines for the field gathering in RZ.
! at arbitrary order.
!
! Developers:
! - Jean-Luc Vay
! - Henri vincenti
! - Mathieu Lobet
! - David Grote
!
! List of subroutines:
!
! - pxr_gete2drz_n_energy_conserving
! - pxr_getb2drz_n_energy_conserving
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
!> David Grote
!
!> @date
!> Creation 2019
!
!> @param[in] np Number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nmodes number of modes (including mode 0)
!> @param[in] nox, noz order of interpolation
!> @param[in] *_nguard number of guard cells
!> @param[in] *_nvalid number of valid cells
!> @param[in] erg, etg, ezg field arrays
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order in r
!> @param[in] l_nodal whether is is nodal or Yee staggered
!
! ________________________________________________________________________________________
subroutine pxr_gete2drz_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, &
                                            xmin, zmin, dx, dz, nmodes, nox, noz, &
                                            erg, erg_nguard, erg_nvalid, &
                                            etg, etg_nguard, etg_nvalid, &
                                            ezg, ezg_nguard, ezg_nvalid, &
                                            l_lower_order_in_v, l_nodal) !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num, cpx
  implicit none

  integer(idp), intent(IN) :: np, nmodes, nox, noz
  integer(idp), intent(IN) :: erg_nguard(2), erg_nvalid(2), etg_nguard(2)
  integer(idp), intent(IN) :: etg_nvalid(2), ezg_nguard(2), ezg_nvalid(2)
  REAL(num), intent(IN), dimension(np) :: xp, yp, zp
  REAL(num), intent(IN OUT), dimension(np) :: ex, ey, ez
  logical(lp), intent(IN)  :: l_lower_order_in_v, l_nodal
  complex(num), intent(IN) :: erg(-erg_nguard(1):erg_nvalid(1)+erg_nguard(1)-1, &
                                  -erg_nguard(2):erg_nvalid(2)+erg_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN) :: etg(-etg_nguard(1):etg_nvalid(1)+etg_nguard(1)-1, &
                                  -etg_nguard(2):etg_nvalid(2)+etg_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN) :: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1, &
                                  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1,0:nmodes-1)
  real(num), intent(IN)    :: xmin, zmin, dx, dz

  real(num) :: stagger_shift
  integer(idp) :: ip, j, l, m, ixmin, ixmax, izmin, izmax, ixmin0, ixmax0
  integer(idp) :: izmin0, izmax0, jj, ll, j0, l0
  real(num) :: dxi, dzi, x, y, z, r, xint, zint
  REAL(num) :: xintsq, oxint, zintsq, ozint, oxintsq, ozintsq, exc, eyc, ezc
  real(num) :: costheta, sintheta
  real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  real(num), dimension(:), allocatable :: sx0, sz0
  real(num), parameter :: onesixth=1./6., twothird=2./3.
  COMPLEX(cpx) :: xy, xy0

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

  do ip=1, np

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
    xy0 = cmplx(costheta, -sintheta)

    z = (zp(ip)-zmin)*dzi

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

    ! Mode m = 0
    do ll = izmin, izmax+1
      do jj = ixmin0, ixmax0
        ex(ip) = ex(ip) + sz(ll)*sx0(jj)*real(erg(j0+jj,l+ll,0)*costheta-etg(j0+jj,l+ll,0)*sintheta)
        ey(ip) = ey(ip) + sz(ll)*sx0(jj)*real(erg(j0+jj,l+ll,0)*sintheta+etg(j0+jj,l+ll,0)*costheta)
      end do
    end do

    do ll = izmin0, izmax0
      do jj = ixmin, ixmax+1
        ez(ip) = ez(ip) + sx(jj)*sz0(ll)*real(ezg(j+jj,l0+ll,0))
      end do
    end do

    ! Modes m>0
    xy = 1.
    do m = 1, nmodes-1
      xy = xy*xy0
      do ll = izmin, izmax+1
        do jj = ixmin0, ixmax0
          exc = real(erg(j+jj,l+ll,m)*xy, num)
          eyc = real(etg(j+jj,l+ll,m)*xy, num)
          ex(ip) = ex(ip) + sx0(jj)*sz(ll)*(exc*costheta - eyc*sintheta)
          ey(ip) = ey(ip) + sx0(jj)*sz(ll)*(exc*sintheta + eyc*costheta)
        end do
      end do

      do ll = izmin0, izmax0
        do jj = ixmin, ixmax+1
          ezc = real(ezg(j+jj,l+ll,m)*xy, num)
          ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezc
        end do
      end do
    end do

  end do
  deallocate(sx0, sz0)

  return
end subroutine pxr_gete2drz_n_energy_conserving

! ________________________________________________________________________________________
!> @brief
!> 2D magnetic field gathering routine for arbitrary orders
!
!> @details
!> This function is not vectorized and not optimized.
!
!> @author
!> Henri Vincenti
!> David Grote
!
!> @date
!> 12/01/2019
!
!> @param[in] np Number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nmodes number of modes (including mode 0)
!> @param[in] nox, noz order of interpolation
!> @param[in] *_nguard number of guard cells
!> @param[in] *_nvalid number of valid cells
!> @param[in] brg, btg, bzg field arrays
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order in r
!> @param[in] l_nodal whether is is nodal or Yee staggered
!
! ________________________________________________________________________________________
subroutine pxr_getb2drz_n_energy_conserving(np, xp, yp, zp, bx, by, bz, &
                                            xmin, zmin, dx, dz, nmodes, nox, noz, &
                                            brg, brg_nguard, brg_nvalid, &
                                            btg, btg_nguard, btg_nvalid, &
                                            bzg, bzg_nguard, bzg_nvalid, &
                                            l_lower_order_in_v, l_nodal) !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num, cpx
  implicit none

  integer(idp), intent(IN) :: np, nmodes, nox, noz
  integer(idp), intent(IN) :: brg_nguard(2), brg_nvalid(2), btg_nguard(2)
  integer(idp), intent(IN) :: btg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  REAL(num), intent(IN), dimension(np) :: xp, yp, zp
  REAL(num), intent(IN OUT), dimension(np) :: bx, by, bz
  logical(lp), intent(IN)  :: l_lower_order_in_v, l_nodal
  complex(num), intent(IN) :: brg(-brg_nguard(1):brg_nvalid(1)+brg_nguard(1)-1, &
                                  -brg_nguard(2):brg_nvalid(2)+brg_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN) :: btg(-btg_nguard(1):btg_nvalid(1)+btg_nguard(1)-1, &
                                  -btg_nguard(2):btg_nvalid(2)+btg_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN) :: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1, &
                                  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1,0:nmodes-1)
  real(num), intent(IN)    :: xmin, zmin, dx, dz

  real(num)    :: stagger_shift
  integer(idp) :: ip, j, l, m, ixmin, ixmax, izmin, izmax
  integer(idp) :: ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0, l0
  real(num) :: dxi, dzi, x, y, z, xint, zint
  REAL(num) :: xintsq, oxint, zintsq, ozint, oxintsq, ozintsq, bxc, byc, bzc
  real(num) :: r, costheta, sintheta
  real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  real(num), dimension(:), allocatable :: sx0, sz0
  real(num), parameter :: onesixth=1./6., twothird=2./3.
  COMPLEX(cpx) :: xy, xy0

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

  sx=0
  sz=0.
  sx0=0.
  sz0=0.

  do ip=1, np

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
    xy0 = cmplx(costheta, -sintheta)

    z = (zp(ip)-zmin)*dzi

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


    ! Mode m = 0
    do ll = izmin0, izmax0
      do jj = ixmin, ixmax+1
        bx(ip) = bx(ip) + sx(jj)*sz0(ll)*real(brg(j+jj,l0+ll,0)*costheta-btg(j+jj,l0+ll,0)*sintheta)
        by(ip) = by(ip) + sx(jj)*sz0(ll)*real(brg(j+jj,l0+ll,0)*sintheta+btg(j+jj,l0+ll,0)*costheta)
      end do
    end do

    do ll = izmin, izmax+1
      do jj = ixmin0, ixmax0
        bz(ip) = bz(ip) + sx0(jj)*sz(ll)*real(bzg(j0+jj,l+ll,0))
      end do
    end do

    ! Modes m>0
    xy = 1.
    do m = 1, nmodes-1
      xy = xy*xy0
      do ll = izmin0, izmax0
        do jj = ixmin, ixmax+1
          bxc = real(brg(j+jj,l+ll,m)*xy, num)
          byc = real(btg(j+jj,l+ll,m)*xy, num)
          bx(ip) = bx(ip) + sx(jj)*sz0(ll)*(bxc*costheta - byc*sintheta)
          by(ip) = by(ip) + sx(jj)*sz0(ll)*(bxc*sintheta + byc*costheta)
        end do
      end do
      do ll = izmin, izmax+1
        do jj = ixmin0, ixmax0
          bzc = real(bzg(j+jj,l+ll,m)*xy, num)
          bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzc
        end do
      end do
    end do

  end do
  deallocate(sx0, sz0)
  return
end subroutine pxr_getb2drz_n_energy_conserving

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
! CHARGE_DEPOSITION_2D.F90
!
! Developers:
! Henri Vincenti, ! Mathieu Lobet
!
! Brief description:
! File containing subroutines for the charge deposition itself in 2D.
!
! List of suboutines:
!
! General order:
! - pxr_depose_rho_n_2dxy
! - pxr_depose_rhoold_n_2dxy
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> 2D Charge deposition at nth order
!>
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE pxr_depose_rho_n_2dxz(rho, np, xp, yp, zp, w, q, xmin, zmin, dx, dz, nx,   &
  nz, nxguard, nzguard, nox, noz, l_particles_weight, l4symtry, l_2drz, type_rz_depose)
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp) :: np, nx, nz, nox, noz, nxguard, nzguard, type_rz_depose
  real(num), dimension(-nxguard:nx+nxguard, 0:0, -nzguard:nz+nzguard), intent(in out) &
  :: rho
  real(num), dimension(np) :: xp, yp, zp, w
  real(num) :: q, dx, dz, xmin, zmin
  LOGICAL(lp)  :: l_particles_weight, l4symtry, l_2drz

  real(num) :: dxi, dzi, xint, zint, oxint, ozint, xintsq, zintsq, oxintsq, ozintsq
  real(num) :: x, z, r, wq, invvol
  real(num) :: sx(-int(nox/2):int((nox+1)/2)), sz(-int(noz/2):int((noz+1)/2))
  real(num), parameter :: onesixth=1./6., twothird=2./3.
  integer(idp) :: j, l, ip, jj, ll, ixmin, ixmax, izmin, izmax

  dxi = 1./dx
  dzi = 1./dz
  invvol = dxi*dzi

  ! Davoine method : limited to order 1 in r
  if (type_rz_depose==2) then
    nox = 1
  endif

  ixmin = -int(nox/2)
  ixmax = int((nox+1)/2)
  izmin = -int(noz/2)
  izmax = int((noz+1)/2)

  do ip=1, np

    ! --- computes current position in grid units
    if (l_2drz) then
      r = sqrt(xp(ip)*xp(ip)+yp(ip)*yp(ip))
      x = (r-xmin)*dxi
      z = (zp(ip)-zmin)*dzi
    else
      x = (xp(ip)-xmin)*dxi
      z = (zp(ip)-zmin)*dzi
    end if

    ! --- applies 4-fold symmetry
    if (l4symtry) then
      x=abs(x)
    end if

    ! --- finds node of cell containing particles for current positions
    ! --- (different for odd/even spline orders)
    if (nox==2*(nox/2)) then
      j=nint(x)
    else
      j=floor(x)
    end if
    if (noz==2*(noz/2)) then
      l=nint(z)
    else
      l=floor(z)
    end if

    ! --- computes distance between particle and node for current positions
    xint = x-j
    zint = z-l

    ! --- computes particles "weights"
    if (l_particles_weight) then
      wq=q*w(ip)*invvol
    else
      wq=q*invvol
    end if

    ! --- computes coefficients for node centered quantities
    if (type_rz_depose == 2) then! Davoine method, modified particle shapes in r
      sx(0) = 1. - xint  + 1./(4*j+2)*( -xint + xint**2 )
      sx(1) = 1. - sx(0)
    else! Standard method, canonical shapes in r
      select case(nox)
      case(0)
        sx( 0) = 1.
      case(1)
        sx( 0) = 1.-xint
        sx( 1) = xint
      case(2)
        xintsq = xint*xint
        sx(-1) = 0.5*(0.5-xint)**2
        sx( 0) = 0.75-xintsq
        sx( 1) = 0.5*(0.5+xint)**2
      case(3)
        oxint = 1.-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.-xint/2)
        sx( 1) = twothird-oxintsq*(1.-oxint/2)
        sx( 2) = onesixth*xintsq*xint
      end select
    endif

    select case(noz)
    case(0)
      sz( 0) = 1.
    case(1)
      sz( 0) = 1.-zint
      sz( 1) = zint
    case(2)
      zintsq = zint*zint
      sz(-1) = 0.5*(0.5-zint)**2
      sz( 0) = 0.75-zintsq
      sz( 1) = 0.5*(0.5+zint)**2
    case(3)
      ozint = 1.-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird-zintsq*(1.-zint/2)
      sz( 1) = twothird-ozintsq*(1.-ozint/2)
      sz( 2) = onesixth*zintsq*zint
    end select

    ! --- add charge density contributions
    do ll = izmin, izmax
      do jj = ixmin, ixmax
        rho(j+jj, 0, l+ll)=rho(j+jj, 0, l+ll)+sx(jj)*sz(ll)*wq
      end do
    end do

  end do

  return
END SUBROUTINE pxr_depose_rho_n_2dxz

! ________________________________________________________________________________________
!> @brief
!> Old 2D Charge deposition subroutine at nth order
!>
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE pxr_depose_rhoold_n_2dxz(rhoold, np, xp, zp, ux, uy, uz, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, nx, nz, nxguard, nzguard, nox, noz, l_particles_weight,       &
  l4symtry)
  USE picsar_precision, ONLY: idp, lp, num
  implicit none
  integer(idp) :: np, nx, nz, nox, noz, nxguard, nzguard
  real(num), dimension(-nxguard:nx+nxguard, 0:0, -nzguard:nz+nzguard), intent(in out) &
  :: rhoold
  real(num), dimension(np) :: xp, zp, w, ux, uy, uz, gaminv
  real(num) :: q, dt, dx, dz, xmin, zmin
  LOGICAL(lp) :: l_particles_weight, l4symtry

  real(num) :: dxi, dzi, xint, zint, oxint, ozint, xintsq, zintsq, oxintsq, ozintsq
  real(num) :: xintold, zintold, oxintold, ozintold
  real(num) :: x, z, xold, zold, wq, invvol, vx, vy, vz
  real(num) :: sx(-int(nox/2):int((nox+1)/2)), sz(-int(noz/2):int((noz+1)/2))
  real(num) :: sxold(-int(nox/2):int((nox+1)/2)), szold(-int(noz/2):int((noz+1)/2))
  real(num), parameter :: onesixth=1./6., twothird=2./3.
  integer(idp) :: j, l, ip, jj, ll, jold, lold, ixmin, ixmax, izmin, izmax, ndt, idt
  real(num) :: dxp, dzp, x0, z0, x1, z1

  dxi = 1./dx
  dzi = 1./dz
  invvol = dxi*dzi

  ixmin = -int(nox/2)
  ixmax = int((nox+1)/2)
  izmin = -int(noz/2)
  izmax = int((noz+1)/2)
  ndt = 1

  do ip=1, np

    vx = ux(ip)*gaminv(ip)
    vy = uy(ip)*gaminv(ip)
    vz = uz(ip)*gaminv(ip)

    x1 = (xp(ip)-xmin)*dxi
    z1 = (zp(ip)-zmin)*dzi
    x0 = x1 - vx*dt*dxi
    z0 = z1 - vz*dt*dzi

    dxp=(x1-x0)/ndt
    dzp=(z1-z0)/ndt

    xold=x0
    zold=z0

    do idt=1, ndt

      if (idt>1) then
        xold=x
        zold=z
      end if
      x=xold+dxp
      z=zold+dzp

      if (l4symtry) then
        x=abs(x)
        xold=abs(xold)
      end if

      ! --- finds node of cell containing particles for current positions
      ! --- (different for odd/even spline orders)
      if (nox==2*(nox/2)) then
        j=nint(x)
      else
        j=floor(x)
      end if
      if (noz==2*(noz/2)) then
        l=nint(z)
      else
        l=floor(z)
      end if

      if (nox==2*(nox/2)) then
        jold=nint(xold)
      else
        jold=floor(xold)
      end if
      if (noz==2*(noz/2)) then
        lold=nint(zold)
      else
        lold=floor(zold)
      end if

      xint = x-j
      zint = z-l
      xintold = xold-jold
      zintold = zold-lold

      if (l_particles_weight) then
        wq=q*w(ip)*invvol
      else
        wq=q*w(1)*invvol
      end if

      select case(nox)
      case(0)
        sxold( 0) = 1.
      case(1)
        sxold( 0) = 1.-xintold
        sxold( 1) = xintold
      case(2)
        xintsq = xintold*xintold
        sxold(-1) = 0.5*(0.5-xintold)**2
        sxold( 0) = 0.75-xintsq
        sxold( 1) = 0.5*(0.5+xintold)**2
      case(3)
        oxintold = 1.-xintold
        xintsq = xintold*xintold
        oxintsq = oxintold*oxintold
        sxold(-1) = onesixth*oxintsq*oxintold
        sxold( 0) = twothird-xintsq*(1.-xintold/2)
        sxold( 1) = twothird-oxintsq*(1.-oxintold/2)
        sxold( 2) = onesixth*xintsq*xintold
      end select

      select case(noz)
      case(0)
        szold( 0) = 1.
      case(1)
        szold( 0) = 1.-zintold
        szold( 1) = zintold
      case(2)
        zintsq = zintold*zintold
        szold(-1) = 0.5*(0.5-zintold)**2
        szold( 0) = 0.75-zintsq
        szold( 1) = 0.5*(0.5+zintold)**2
      case(3)
        ozintold = 1.-zintold
        zintsq = zintold*zintold
        ozintsq = ozintold*ozintold
        szold(-1) = onesixth*ozintsq*ozintold
        szold( 0) = twothird-zintsq*(1.-zintold/2)
        szold( 1) = twothird-ozintsq*(1.-ozintold/2)
        szold( 2) = onesixth*zintsq*zintold
      end select

      select case(nox)
      case(0)
        sx( 0) = 1.
      case(1)
        sx( 0) = 1.-xint
        sx( 1) = xint
      case(2)
        xintsq = xint*xint
        sx(-1) = 0.5*(0.5-xint)**2
        sx( 0) = 0.75-xintsq
        sx( 1) = 0.5*(0.5+xint)**2
      case(3)
        oxint = 1.-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.-xint/2)
        sx( 1) = twothird-oxintsq*(1.-oxint/2)
        sx( 2) = onesixth*xintsq*xint
      end select

      select case(noz)
      case(0)
        sz( 0) = 1.
      case(1)
        sz( 0) = 1.-zint
        sz( 1) = zint
      case(2)
        zintsq = zint*zint
        sz(-1) = 0.5*(0.5-zint)**2
        sz( 0) = 0.75-zintsq
        sz( 1) = 0.5*(0.5+zint)**2
      case(3)
        ozint = 1.-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.-zint/2)
        sz( 1) = twothird-ozintsq*(1.-ozint/2)
        sz( 2) = onesixth*zintsq*zint
      end select

      do ll = izmin, izmax
        do jj = ixmin, ixmax

          rhoold(jold+jj, 0, lold+ll) = rhoold(jold+jj, 0, lold+ll) +                 &
          sxold(jj)*szold(ll)*wq

        end do
      end do
    end do
  end do

  return
END SUBROUTINE pxr_depose_rhoold_n_2dxz


! ______________________________________________________________________________
!> @brief
!> Applies the inverse cell volume scaling to charge density.
!>
!> @details
!> Applies the inverse cell volume scaling. It is more efficient to apply
!> the scaling afterward rather than with the particles.
! ________________________________________________________________________________________
SUBROUTINE apply_rz_volume_scaling_rho( rho, rho_nguard, rho_nvalid, rmin, dr, type_rz_depose)     !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  USE constants, ONLY: pi
  implicit none
  INTEGER(idp), intent(in) :: rho_nguard(2), rho_nvalid(2)
  REAL(num), intent(IN OUT):: rho(-rho_nguard(1):rho_nvalid(1)+rho_nguard(1)-1,  &
                                  -rho_nguard(2):rho_nvalid(2)+rho_nguard(2)-1 )
  real(num), intent(in)    :: dr, rmin
  INTEGER(idp), intent(in) :: type_rz_depose

  INTEGER(idp) :: j
  real(num):: r

  if (rmin == 0.) then
     rho(1:rho_nguard(1),:) = rho(1:rho_nguard(1),:) + rho(-1:-rho_nguard(1):-1,:)
  end if

  ! In rz geometry, divide the current by the cell volume

  ! In the lower guard cells in x
  do j=-rho_nguard(1),-1
     r = abs(rmin + j*dr)
     rho(j,:) = rho(j,:)/(2.*pi*r)
  end do

  ! On the lower boundary
  j = 0
  if (rmin == 0.) then
     ! On axis
     if (type_rz_depose == 1) then ! Verboncoeur JCP 164, 421-427 (2001) : corrected volumes
        rho(j,:) = rho(j,:)/(pi*dr/3.)
     else                          ! Standard volume
        rho(j,:) = rho(j,:)/(pi*dr/4.)
     endif
  else
     ! Not the axis
     r = abs(rmin + j*dr)
     rho(j,:) = rho(j,:)/(2.*pi*r)
  end if

  ! In the rest of the grid
  do j=1,rho_nvalid(1) + rho_nguard(1)-1
     r = abs(rmin + j*dr)
     rho(j,:) = rho(j,:)/(2.*pi*r)
  end do

  return
END SUBROUTINE apply_rz_volume_scaling_rho

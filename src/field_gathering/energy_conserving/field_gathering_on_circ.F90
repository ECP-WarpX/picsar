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
! - David Grote
!
! List of subroutines:
!
! - pxr_getf2drz_n_energy_conserving
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> RZ field non-optimized gathering routine (for both E and B)
!
!> @details
!> This function is similar to what is implemented in WARP
!
!> @author
!> Jean-Luc Vay
!> David Grote
!
!> @date
!> Creation 2019
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nmodes number of modes (including mode 0)
!> @param[in] nox, noz order of interpolation
!> @param[in] *_nguard number of guard cells
!> @param[in] *_nvalid number of valid cells
!> @param[in] frg, ftg, flg field arrays
!
! ________________________________________________________________________________________
subroutine pxr_getf2drz_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, xmin, zmin,  &
  dx, dz, nmodes, nox, noz, frg, frg_nguard, frg_nvalid, ftg, ftg_nguard, ftg_nvalid, flg,      &
  flg_nguard, flg_nvalid)     !#do not wrap
  USE picsar_precision, ONLY: idp, num, cpx
  implicit none

  INTEGER(idp)             :: np, nmodes, nox, noz
  INTEGER(idp), intent(IN) :: frg_nguard(2), frg_nvalid(2), ftg_nguard(2),            &
  ftg_nvalid(2), flg_nguard(2), flg_nvalid(2)
  REAL(num), dimension(np) :: xp, yp, zp, ex, ey, ez
  REAL(num)                :: stagger_shift
  COMPLEX(cpx), intent(IN):: frg(-frg_nguard(1):frg_nvalid(1)+frg_nguard(1)-1,        &
  -frg_nguard(2):frg_nvalid(2)+frg_nguard(2)-1, 0:nmodes-1)
  COMPLEX(cpx), intent(IN):: ftg(-ftg_nguard(1):ftg_nvalid(1)+ftg_nguard(1)-1,        &
  -ftg_nguard(2):ftg_nvalid(2)+ftg_nguard(2)-1, 0:nmodes-1)
  COMPLEX(cpx), intent(IN):: flg(-flg_nguard(1):flg_nvalid(1)+flg_nguard(1)-1,        &
  -flg_nguard(2):flg_nvalid(2)+flg_nguard(2)-1, 0:nmodes-1)
  REAL(num)                :: xmin, zmin, dx, dz
   
  INTEGER(idp) :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                  ixmin0, ixmax0, izmin0, izmax0, jj, ll, m
  REAL(num) :: dxi, dzi, x, y, z, xint, zint, r, costheta, sintheta, invr, stot
  REAL(num) :: xintsq, oxint, zintsq, ozint, oxintsq, ozintsq, signx, exc, eyc, ezc
  REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
  REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
  REAL(num), parameter :: onesixth = 1./6., twothird = 2./3.
  COMPLEX(cpx) :: xy, xy0

  dxi = 1./dx
  dzi = 1./dz

  ixmin = -int(nox/2)
  ixmax =  int((nox+1)/2)
  izmin = -int(noz/2)
  izmax =  int((noz+1)/2)

  signx = 1.
  
  do ip = 1, np

    x = xp(ip)
    y = yp(ip)
    r = sqrt(x*x + y*y)
    if (r*dxi > 1.e-20) then
      invr = 1./r
      costheta = x*invr
      sintheta = y*invr
    else  
      costheta = 1.
      sintheta = 0.
    end if
    r = (r - xmin)*dxi
    xy0 = cmplx(costheta, -sintheta)

    z = (zp(ip) - zmin)*dzi

    ! --- finds node of cell containing particles for current positions 
    ! --- (different for odd/even spline orders)
    if (nox == 2*(nox/2)) then
      j = nint(r)
    else
      j = floor(r)
    end if
    if (noz == 2*(noz/2)) then
      l = nint(z)
    else
      l = floor(z)
    end if

    xint = r - j
    zint = z - l

    select case(nox)
     case(0)
      sx( 0) = 1.
     case(1)
      sx( 0) = 1. - xint
      sx( 1) = xint
     case(2)
      xintsq = xint*xint
      sx(-1) = 0.5*(0.5 - xint)**2
      sx( 0) = 0.75 - xintsq
      sx( 1) = 0.5*(0.5 + xint)**2
     case(3)
      oxint = 1. - xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1) = onesixth*oxintsq*oxint
      sx( 0) = twothird - xintsq*(1. - xint/2)
      sx( 1) = twothird - oxintsq*(1. - oxint/2)
      sx( 2) = onesixth*xintsq*xint
    end select        

    select case(noz)
     case(0)
      sz( 0) = 1.
     case(1)
      sz( 0) = 1. - zint
      sz( 1) = zint
     case(2)
      zintsq = zint*zint
      sz(-1) = 0.5*(0.5 - zint)**2
      sz( 0) = 0.75 - zintsq
      sz( 1) = 0.5*(0.5 + zint)**2
     case(3)
      ozint = 1. - zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1) = onesixth*ozintsq*ozint
      sz( 0) = twothird - zintsq*(1. - zint/2)
      sz( 1) = twothird - ozintsq*(1. - ozint/2)
      sz( 2) = onesixth*zintsq*zint
    end select        

    ! Mode m = 0
    do ll = izmin, izmax
      do jj = ixmin, ixmax
        ex(ip) = ex(ip) + sx(jj)*sz(ll)*(frg(j+jj,l+ll,0)*costheta - ftg(j+jj,l+ll,0)*sintheta)
        ey(ip) = ey(ip) + sx(jj)*sz(ll)*(frg(j+jj,l+ll,0)*sintheta + ftg(j+jj,l+ll,0)*costheta)
        ez(ip) = ez(ip) + sx(jj)*sz(ll)*flg(j+jj,l+ll,0)
      end do
    end do
    
    ! Modes m>0
    xy = 1.
    do m = 1, nmodes-1
      xy = xy*xy0
      do ll = izmin, izmax
        do jj = ixmin, ixmax
          exc = real(frg(j+jj,l+ll,m)*xy, num)
          eyc = real(ftg(j+jj,l+ll,m)*xy, num) 
          ezc = real(flg(j+jj,l+ll,m)*xy, num)
          ex(ip) = ex(ip) + sx(jj)*sz(ll)*(exc*costheta - eyc*sintheta)
          ey(ip) = ey(ip) + sx(jj)*sz(ll)*(exc*sintheta + eyc*costheta)
          ez(ip) = ez(ip) + sx(jj)*sz(ll)*ezc
        end do
      end do
    end do

  end do

  return
end subroutine pxr_getf2drz_n_energy_conserving

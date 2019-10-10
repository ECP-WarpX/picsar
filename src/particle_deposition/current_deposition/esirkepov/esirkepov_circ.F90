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
! ESIRKEPOV_CURRENT_DEPOSITION_CIRC.F90
!
! Developers
! Jean-Luc Vay
! David Grote
! Henri Vincenti
!
! Description:
! This file contains subroutines for Esirkepov current deposition in RZ multimode.
!
! List of subroutines:
!
! - pxr_depose_jrjtjz_esirkepov_n_2d_circ
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> RZ multimode Current deposition esirkepov n order (from 0 to 3)
!
!> @details
!> This subroutine is adapted from the version of WARP.

!> @author
!> Jean-Luc Vay
!> David Grote
!> Henri Vincenti

!> @date
!> Creation 2019

! Input parameters:
!> @param[inout] jr r-current component (3D array)
!> @param[in] jr_nguard number of guard cells of the jr array in each direction
!> (1d array containing 2 integers)
!> @param[in] jr_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jr array (1d array containing 2 integers)
!> @param[inout] jt theta-current component (3D array)
!> @param[in] jt_nguard number of guard cells of the jt array in each direction
!> (1d array containing 2 integers)
!> @param[in] jt_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jt array (1d array containing 2 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 2 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 2 integers)
!> @param[in] nmodes number of theta modes
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] rmin, zmin minimal boundaries of the tile
!> @param[in] dt, dr, dz time and space discretization
!> @param[in] nox, noz shape factor order
!> @param[in] l_particles_weight flags whether each particle has its own weight
!> @param[in] type_rz_depose Flags what order to use for radial deposition
! ________________________________________________________________________________________
subroutine pxr_depose_jrjtjz_esirkepov_n_2d_circ( &
              jr, jr_nguard, jr_nvalid, jt, jt_nguard, jt_nvalid, jz, jz_nguard, jz_nvalid, &
              nmodes, &
              np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
              rmin, zmin, dt, dr, dz, &
              nox, noz, l_particles_weight, type_rz_depose) !#do not wrap
  use constants, only: clight
  use picsar_precision, only: idp, num
  implicit none
  integer(idp), intent(in) :: jr_nguard(2), jr_nvalid(2)
  integer(idp), intent(in) :: jt_nguard(2), jt_nvalid(2)
  integer(idp), intent(in) :: jz_nguard(2), jz_nvalid(2)
  complex(num), intent(IN OUT):: jr(-jr_nguard(1):jr_nvalid(1)+jr_nguard(1)-1, &
                                    -jr_nguard(2):jr_nvalid(2)+jr_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN OUT):: jt(-jt_nguard(1):jt_nvalid(1)+jt_nguard(1)-1, &
                                    -jt_nguard(2):jt_nvalid(2)+jt_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1, &
                                    -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,0:nmodes-1)
  integer(idp) :: np, nox, noz, nmodes, type_rz_depose
  real(num), dimension(np) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  real(num) :: q, dt, dr, dz, rmin, zmin
  logical(idp) :: l_particles_weight

  real(num) :: dri, dzi, dtsdr, dtsdz, rint, zint, invr, dti
  real(num), dimension(:,:), allocatable :: sdr, sdz
  real(num), dimension(1:nmodes-1) :: wqt, invdtm
  real(num) :: xold, yold, zold, rold, xmid, ymid, zmid, x, y, z, r, c, s, wq, wqx, wqz, &
       tmp, vx, vy, vz, dts2dr, dts2dz, &
       s1x, s2x, s1z, s2z, invvol, invdtdr, invdtdz, &
       orint, ozint, rintsq, zintsq, orintsq, ozintsq, &
       dtsdr0, dtsdz0, dts2dr0, dts2dz0, rmid, cold, cmid, sold, smid
  real(num), parameter :: onesixth = 1._num/6._num
  real(num), parameter :: twothird = 2._num/3._num
  real(num), dimension(:), allocatable :: sr, sr0, dsr, sz, sz0, dsz
  integer(idp) :: iixp0, ikxp0, iixp, ikxp, ip, dir, diz, idr, idz, i, k, ic, kc, &
       irmin, irmax, izmin, izmax, icell, ncells, m, ndtodr, ndtodz, &
                   xl, xu, zl, zu
  complex(num) :: xymid, xymid0, xy, xy0, xyold, xyold0, im
  integer:: alloc_status

  im = cmplx(0._num, 1._num)

  ndtodr = int(clight*dt/dr)
  ndtodz = int(clight*dt/dz)
  xl = -int(nox/2) - 1 - ndtodr
  xu = int((nox + 1)/2) + 1 + ndtodr
  zl = -int(noz/2) - 1 - ndtodz
  zu = int((noz + 1)/2) + 1 + ndtodz
  allocate(sdr(xl:xu,zl:zu),sdz(xl:xu,zl:zu), &
           sr(xl:xu), sr0(xl:xu), dsr(xl:xu), &
           sz(zl:zu), sz0(zl:zu), dsz(zl:zu), stat=alloc_status)
  if (alloc_status /= 0) then
    print*,"Error:pxr_depose_jrjtjz_esirkepov_n_2d_circ: sdr et al could not be allocated"
    stop
  endif

  sr0 = 0._num
  sz0 = 0._num
  sdr = 0._num
  sdz = 0._num

  ! Davoine method : limited to order 1 in r
  if (type_rz_depose == 2) then
     nox = 1
  endif

  dri = 1._num/dr
  dzi = 1._num/dz
  dti = 1._num/dt
  invvol = 1._num/(dr*dz)
  dtsdr0 = dt*dri
  dtsdz0 = dt*dzi
  dts2dr0 = 0.5_num*dtsdr0
  dts2dz0 = 0.5_num*dtsdz0
  invdtdr = 1._num/(dt*dz)
  invdtdz = 1._num/(dt*dr)
  do m = 1, nmodes - 1
     invdtm(m) = invdtdr/m
  enddo

  do ip=1, np

     ! --- computes current position in grid units
     x = xp(ip)
     y = yp(ip)
     xmid = 0.5_num*x
     ymid = 0.5_num*y
     r = sqrt(x*x + y*y)
     if (r*dri > 1.e-10_num) then
        invr = 1._num/r
        c = x*invr 
        s = y*invr
     else
        c = 1._num
        s = 0._num
     end if
     xy0 = cmplx(c, s)
     x = r
     x = x*dri
     z = zp(ip)*dzi

     ! --- computes velocity
     vx = uxp(ip)*gaminv(ip)
     vy = uyp(ip)*gaminv(ip)
     vz = uzp(ip)*gaminv(ip)

     ! --- computes old position in grid units
     xold = xp(ip) - dt*vx
     yold = yp(ip) - dt*vy
     rold = sqrt(xold*xold + yold*yold)
     if (rold*dri > 1.e-10_num) then
        invr = 1._num/rold
        cold = xold*invr 
        sold = yold*invr
     else
        cold = 1._num
        sold = 0._num
     end if
     xyold0 = cmplx(cold, sold)
     xmid = xmid + 0.5_num*xold
     ymid = ymid + 0.5_num*yold
     rmid = sqrt(xmid*xmid + ymid*ymid)
     if (rmid*dri > 1.e-10_num) then
        invr = 1._num/rmid
        cmid = xmid*invr 
        smid = ymid*invr
     else
        cmid = 1._num
        smid = 0._num
     end if
     xymid0 = cmplx(cmid, smid)
     xold = rold*dri
     vy = -vx*smid + vy*cmid
     vx = (x - xold)*dr*dti
     zold = z - dtsdz0*vz

     ! --- sets positions relative to grid start
     x = x - rmin*dri
     z = z - zmin*dzi
     xold = xold - rmin*dri
     zold = zold - zmin*dzi

     ! computes maximum number of cells traversed by particle in a given dimension
     ncells = 1!+max( int(abs(x-xold)), int(abs(z-zold)))
     dtsdr = dtsdr0/ncells
     dtsdz = dtsdz0/ncells
     dts2dr = dts2dr0/ncells
     dts2dz = dts2dz0/ncells

     x = xold
     z = zold

     do icell = 1, ncells

        xold = x
        zold = z

        x = x + dtsdr*vx
        z = z + dtsdz*vz

        ! --- computes particles "weights"
        if (l_particles_weight) then
           wq = q*w(ip)
        else
           wq = q*w(1)
        end if
        wqx = wq*invdtdr
        wqz = wq*invdtdz
        wqt(:) = wq*invdtm(:)

        ! --- finds node of cell containing particles for current positions 
        ! --- (different for odd/even spline orders)
        if (nox == 2*(nox/2)) then
           iixp0 = nint(x)
        else
           iixp0 = floor(x)
        end if
        if (noz == 2*(noz/2)) then
           ikxp0 = nint(z)
        else
           ikxp0 = floor(z)
        end if

        ! --- computes distance between particle and node for current positions
        rint = x - iixp0
        zint = z - ikxp0

        ! --- computes coefficients for node centered quantities
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sr0(0) = 1._num - rint + 1._num/(4_num*iixp0 + 2_num)*( -rint + rint**2 )
           sr0(1) = 1._num - sr0(0)
        else                          ! Standard method, canonical shapes in r
           select case(nox)
           case(0)
              sr0( 0) = 1._num
           case(1)
              sr0( 0) = 1._num - rint
              sr0( 1) = rint
           case(2)
              rintsq = rint*rint
              sr0(-1) = 0.5_num*(0.5_num - rint)**2
              sr0( 0) = 0.75_num - rintsq
              sr0( 1) = 0.5_num*(0.5_num + rint)**2
           case(3)
              orint = 1._num - rint
              rintsq = rint*rint
              orintsq = orint*orint
              sr0(-1) = onesixth*orintsq*orint
              sr0( 0) = twothird - rintsq*(1._num - rint/2_num)
              sr0( 1) = twothird - orintsq*(1._num - orint/2_num)
              sr0( 2) = onesixth*rintsq*rint
           end select
        endif

        select case(noz)
        case(0)
           sz0( 0) = 1._num
        case(1)
           sz0( 0) = 1._num - zint
           sz0( 1) = zint
        case(2)
           zintsq = zint*zint
           sz0(-1) = 0.5_num*(0.5_num - zint)**2
           sz0( 0) = 0.75_num-zintsq
           sz0( 1) = 0.5_num*(0.5_num + zint)**2
        case(3)
           ozint = 1._num - zint
           zintsq = zint*zint
           ozintsq = ozint*ozint
           sz0(-1) = onesixth*ozintsq*ozint
           sz0( 0) = twothird - zintsq*(1._num - zint/2_num)
           sz0( 1) = twothird - ozintsq*(1._num - ozint/2_num)
           sz0( 2) = onesixth*zintsq*zint
        end select

        ! --- finds node of cell containing particles for old positions 
        ! --- (different for odd/even spline orders)
        if (nox == 2*(nox/2)) then
           iixp = nint(xold)
        else
           iixp = floor(xold)
        end if
        if (noz == 2*(noz/2)) then
           ikxp = nint(zold)
        else
           ikxp = floor(zold)
        end if

        ! --- computes distance between particle and node for old positions
        rint = xold - iixp
        zint = zold - ikxp

        ! --- computes node separation between old and current positions
        dir = iixp - iixp0
        diz = ikxp - ikxp0

        ! --- zero out coefficients (needed because of different dir and diz for each particle)
        sr = 0._num
        sz = 0._num

        ! --- computes coefficients for quantities centered between nodes
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sr(0+dir) = 1._num - rint + 1._num/(4_num*iixp + 2_num)*( -rint + rint**2 )
           sr(1+dir) = 1._num - sr(0+dir)
        else! Standard method, canonical shapes in r 
           select case(nox)
           case(0)
              sr( 0+dir) = 1._num
           case(1)
              sr( 0+dir) = 1._num - rint
              sr( 1+dir) = rint
           case(2)
              rintsq = rint*rint
              sr(-1+dir) = 0.5_num*(0.5_num - rint)**2
              sr( 0+dir) = 0.75_num-rintsq
              sr( 1+dir) = 0.5_num*(0.5_num + rint)**2
           case(3)
              orint = 1._num - rint
              rintsq = rint*rint
              orintsq = orint*orint
              sr(-1+dir) = onesixth*orintsq*orint
              sr( 0+dir) = twothird - rintsq*(1._num - rint/2_num)
              sr( 1+dir) = twothird - orintsq*(1._num - orint/2_num)
              sr( 2+dir) = onesixth*rintsq*rint
           end select
        endif

        select case(noz)
        case(0)
           sz( 0+diz) = 1._num
        case(1)
           sz( 0+diz) = 1._num - zint
           sz( 1+diz) = zint
        case(2)
           zintsq = zint*zint
           sz(-1+diz) = 0.5_num*(0.5_num - zint)**2
           sz( 0+diz) = 0.75_num - zintsq
           sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
        case(3)
           ozint = 1._num - zint
           zintsq = zint*zint
           ozintsq = ozint*ozint
           sz(-1+diz) = onesixth*ozintsq*ozint
           sz( 0+diz) = twothird - zintsq*(1._num - zint/2_num)
           sz( 1+diz) = twothird - ozintsq*(1._num - ozint/2_num)
           sz( 2+diz) = onesixth*zintsq*zint
        end select

        ! --- computes coefficients difference
        dsr = sr - sr0
        dsz = sz - sz0

        ! --- computes min/max positions of current contributions
        irmin = min(0, dir) - int(nox/2)
        irmax = max(0, dir) + int((nox + 1)/2)
        izmin = min(0, diz) - int(noz/2)
        izmax = max(0, diz) + int((noz + 1)/2)

        ! --- add current contributions
        ! -- NB : the current is later divided by the cylindrical cell volume in applybc_j
        do k=izmin, izmax
           do i=irmin, irmax
              ic = iixp0 + i
              kc = ikxp0 + k

              ! -- Jr
              if(i < irmax) then
                 sdr(i,k)  = wqx*dsr(i)*( sz0(k) + 0.5_num*dsz(k) )    ! Wr coefficient from esirkepov
                 if (i > irmin) sdr(i,k) = sdr(i,k) + sdr(i-1,k)         ! Integration of Wr along r
                 jr(ic,kc,0) = jr(ic,kc,0) + sdr(i,k)              ! Deposition on the mode m = 0
                 xymid = xymid0 ! Throughout the following loop, xymid takes the value e^{i m theta}
                 do m = 1, nmodes - 1                                ! Deposition on the modes m>0
                    jr(ic,kc,m) = jr(ic,kc,m) + 2._num*sdr(i,k)*xymid
                    ! The factor 2 comes from the normalization of the modes
                    xymid = xymid*xymid0
                 enddo
              end if

              ! -- Jtheta
              ! Mode m = 0 : similar to the 2D Esirkepov scheme
              jt(ic,kc,0) = jt(ic,kc,0) + wq*vy*invvol/ncells* &
                   ( (sz0(k) + 0.5_num*dsz(k))*sr0(i) + (0.5_num*sz0(k) + 1._num/3._num*dsz(k))*dsr(i) )
              ! Mode m > 0 : see Davidson et al. JCP 281 (2014)
              xy = xy0 
              xymid = xymid0 
              xyold = xyold0
              ! Throughout the following loop, xy_ takes the value e^{i m theta_}
              do m = 1, nmodes - 1
                 jt(ic,kc,m) = jt(ic,kc,m) - 2_num*im*(ic + rmin*dri)*wqt(m) * &
                      ( sr0(i)*sz0(k)*(xy - xymid) + sr(i)*sz(k)*(xymid - xyold) )
                 ! The factor 2 comes from the normalization of the modes
                 ! The minus sign comes from the different convention with respect to Davidson et al.
                 xy = xy*xy0 
                 xymid = xymid*xymid0 
                 xyold = xyold*xyold0
              enddo

              ! -- Jz
              if(k < izmax) then
                 sdz(i,k)  = wqz*dsz(k)*(sr0(i) + 0.5_num*dsr(i))        ! Wz coefficient from esirkepov
                 if (k > izmin) sdz(i,k) = sdz(i,k) + sdz(i,k-1)         ! Integration of Wz along z
                 jz(ic,kc,0) = jz(ic,kc,0) + sdz(i,k)              ! Deposition on the mode m=0
                 xymid = xymid0 ! Throughout the following loop, xymid takes the value e^{i m theta}
                 do m = 1, nmodes - 1                                ! Deposition on the modes m>0
                    jz(ic,kc,m) = jz(ic,kc,m) + 2._num*sdz(i,k)*xymid
                    ! The factor 2 comes from the normalization of the modes
                    xymid = xymid*xymid0
                 enddo
              end if
           end do
        end do

     end do

  end do

  deallocate(sdr, sdz, sr, sr0, dsr, sz, sz0, dsz)

  return
end subroutine pxr_depose_jrjtjz_esirkepov_n_2d_circ

! ______________________________________________________________________________
!> @brief
!> Applies the inverse cell volume scaling to current density.
!>
!> @details
!> Applies the inverse cell volume scaling. It is more efficient to apply
!> the scaling afterward rather than with the particles.
! ________________________________________________________________________________________
subroutine apply_2dcirc_volume_scaling_j( &
              jr, jr_nguard, jr_nvalid, jt, jt_nguard, jt_nvalid, jz, jz_nguard, jz_nvalid, &
              nmodes, &
              rmin, dr, &
              type_rz_depose) !#do not wrap
  use constants, only: pi
  use picsar_precision, only: idp, num
  implicit none

  integer(idp), intent(in) :: jr_nguard(2), jr_nvalid(2)
  integer(idp), intent(in) :: jt_nguard(2), jt_nvalid(2)
  integer(idp), intent(in) :: jz_nguard(2), jz_nvalid(2)
  complex(num), intent(IN OUT):: jr(-jr_nguard(1):jr_nvalid(1)+jr_nguard(1)-1, &
                                    -jr_nguard(2):jr_nvalid(2)+jr_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN OUT):: jt(-jt_nguard(1):jt_nvalid(1)+jt_nguard(1)-1, &
                                    -jt_nguard(2):jt_nvalid(2)+jt_nguard(2)-1,0:nmodes-1)
  complex(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1, &
                                    -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,0:nmodes-1)
  integer(idp) :: nmodes, type_rz_depose
  real(num) :: dr, rmin

  integer(idp) :: j, m, ifact
  real(num) :: r
  complex(num) :: I = (0._num, 1._num)

  ! MODE 0 : Fetch the current deposited in the guards cells and add it to the grid (fold back)

  ! In rz geometry, for the guards cells below the axis
  if (rmin == 0._num) then
    do m = 0, nmodes-1
      if (mod(m, 2) == 0) then
        ifact = 1
      else
        ifact = -1
      end if
      ! Fields that are located on the boundary
      jt(1:jt_nguard(1),:,m) = jt(1:jt_nguard(1),:,m) + ifact*jt(-1:-jt_nguard(1):-1,:,m)
      jz(1:jz_nguard(1),:,m) = jz(1:jz_nguard(1),:,m) + ifact*jz(-1:-jz_nguard(1):-1,:,m)
      ! Fields that are located off the boundary
      jr(0:jr_nguard(1)-1,:,m) = jr(0:jr_nguard(1)-1,:,m) - ifact*jr(-1:-jr_nguard(1):-1,:,m)
    end do
  endif

  ! Divide the current by the cell volume

  ! -- Jr

  ! Since Jr is not cell centered in r, no need for distinction
  ! between on axis and off-axis factors
  do j=-jr_nguard(1), jr_nvalid(1)+jr_nguard(1)-1
     r = abs(rmin + (j + 0.5_num)*dr)
     jr(j,:,:) = jr(j,:,:)/(2._num*pi*r)
  end do

  ! -- Jtheta and Jz

  ! In the lower guard cells in radius
  do j=-jt_nguard(1), -1
     r = abs(rmin + j*dr)
     jt(j,:,:) = jt(j,:,:)/(2._num*pi*r)
  end do
  do j=-jz_nguard(1), -1
     r = abs(rmin + j*dr)
     jz(j,:,:) = jz(j,:,:)/(2._num*pi*r)
  end do

  ! On the lower boundary
  if (rmin == 0._num) then
     ! On axis
     ! Jz, mode 0
     if (type_rz_depose == 1) then ! Verboncoeur JCP 164, 421-427 (2001) : corrected volume
        jz(0,:,0) = jz(0,:,0)/(pi*dr/3._num)
     else                          ! Standard volume
        jz(0,:,0) = jz(0,:,0)/(pi*dr/4._num)
     endif
     ! Jz, modes > 0
     jz(0,:,1:) = 0._num ! Mode > 0 : Jz is zero on axis.
     ! Jt, mode 0 and modes > 1
     jt(0,:,0) = 0._num ! Mode 0 : Jt is zero on axis.
     ! Jt, mode 1
     if (nmodes > 1) jt(0,:,1) = -I*jr(0,:,1)
     ! Because the previous line uses Jr, it is important that Jr be properly calculated first
     if (nmodes > 2) jt(0,:,2:) = 0._num ! Modes > 1 : Jt = 0.
  else
     ! Not the axis
     r = abs(rmin + 0*dr)
     jt(0,:,:) = jt(0,:,:)/(2._num*pi*r)
     jz(0,:,:) = jz(0,:,:)/(2._num*pi*r)
  end if

  ! In the rest of the grid
  do j=1, jt_nvalid(1) + jt_nguard(1) - 1
     r = abs(rmin + j*dr)
     jt(j,:,:) = jt(j,:,:)/(2._num*pi*r)
  end do
  do j=1, jz_nvalid(1) + jz_nguard(1) - 1
     r = abs(rmin + j*dr)
     jz(j,:,:) = jz(j,:,:)/(2._num*pi*r)
  end do

  return
end subroutine apply_2dcirc_volume_scaling_j

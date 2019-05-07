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
!> @param[in] xmin, zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nox, noz shape factor order
!> @param[in] l_particles_weight flags whether each particle has its own weight
!> @param[in] type_rz_depose Flags what order to use for radial deposition
! ________________________________________________________________________________________
subroutine pxr_depose_jrjtjz_esirkepov_n_2d_circ( &
              jr, jr_nguard, jr_nvalid, jt, jt_nguard, jt_nvalid, jz, jz_nguard, jz_nvalid, &
              nmodes, &
              np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
              xmin,zmin,dt,dx,dz, &
              nox,noz,l_particles_weight,type_rz_depose) !#do not wrap
  use constants, only: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp), intent(in) :: jr_nguard(2), jr_nvalid(2)
  INTEGER(idp), intent(in) :: jt_nguard(2), jt_nvalid(2)
  INTEGER(idp), intent(in) :: jz_nguard(2), jz_nvalid(2)
  COMPLEX(num), intent(IN OUT):: jr(-jr_nguard(1):jr_nvalid(1)+jr_nguard(1)-1, &
                                    -jr_nguard(2):jr_nvalid(2)+jr_nguard(2)-1,0:nmodes-1)
  COMPLEX(num), intent(IN OUT):: jt(-jt_nguard(1):jt_nvalid(1)+jt_nguard(1)-1, &
                                    -jt_nguard(2):jt_nvalid(2)+jt_nguard(2)-1,0:nmodes-1)
  COMPLEX(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1, &
                                    -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,0:nmodes-1)
  integer(idp) :: np,nox,noz,nmodes,type_rz_depose
  real(num), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
  real(num) :: q,dt,dx,dz,xmin,zmin
  logical(idp) :: l_particles_weight

  real(num) :: dxi,dzi,dtsdx,dtsdz,xint,yint,zint, invr, dti
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num), dimension(1:nmodes-1) :: wqt, invdtm
  real(num) :: xold,yold,zold,rold,xmid,ymid,zmid,x,y,z,r,c,s,wq,wqx,wqz, &
       tmp,vx,vy,vz,dts2dx,dts2dz, &
       s1x,s2x,s1z,s2z,invvol,invdtdx,invdtdz, &
       oxint,ozint,xintsq,zintsq,oxintsq,ozintsq, &
       dtsdx0,dtsdz0,dts2dx0,dts2dz0,rmid,cold,cmid,sold,smid
  real(num), parameter :: onesixth=1./6.,twothird=2./3.
  real(num), dimension(:), allocatable :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp) :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc, &
       ixmin, ixmax, izmin, izmax, icell, ncells, m, ndtodx, ndtodz, &
                   xl,xu,zl,zu
  complex(num) :: xymid,xymid0,xy,xy0,xyold,xyold0, im
  integer:: alloc_status

  im = cmplx(0.,1.)

  ndtodx = int(clight*dt/dx)
  ndtodz = int(clight*dt/dz)
  xl = -int(nox/2)-1-ndtodx
  xu = int((nox+1)/2)+1+ndtodx
  zl = -int(noz/2)-1-ndtodz
  zu = int((noz+1)/2)+1+ndtodz
  allocate(sdx(xl:xu,zl:zu),sdz(xl:xu,zl:zu), &
           sx(xl:xu), sx0(xl:xu), dsx(xl:xu), &
           sz(zl:zu), sz0(zl:zu), dsz(zl:zu), stat=alloc_status)
  if (alloc_status /= 0) then
    print*,"Error:pxr_depose_jrjtjz_esirkepov_n_2d_circ: sdx et al could not be allocated"
    stop
  endif

  sx0=0.;sz0=0.
  sdx=0.;sdz=0.

  ! Davoine method : limited to order 1 in r
  if (type_rz_depose==2) then
     nox = 1
  endif

  dxi = 1./dx
  dzi = 1./dz
  dti = 1./dt
  invvol = 1./(dx*dz)
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  dts2dx0 = 0.5*dtsdx0
  dts2dz0 = 0.5*dtsdz0
  invdtdx = 1./(dt*dz)
  invdtdz = 1./(dt*dx)
  do m = 1, nmodes-1
     invdtm(m) = invdtdx/m
  enddo

  do ip=1,np

     ! --- computes current position in grid units
     x = xp(ip)
     y = yp(ip)
     xmid = 0.5*x
     ymid = 0.5*y
     r=sqrt(x*x+y*y)
     if (r*dxi>1.e-10) then
        invr = 1./r
        c = x*invr 
        s = y*invr
     else
        c = 1.
        s = 0.
     end if
     xy0 = cmplx(c,s)
     x = r
     x = x*dxi
     z = zp(ip)*dzi

     ! --- computes velocity
     vx = uxp(ip)*gaminv(ip)
     vy = uyp(ip)*gaminv(ip)
     vz = uzp(ip)*gaminv(ip)

     ! --- computes old position in grid units
     xold = xp(ip)-dt*vx
     yold = yp(ip)-dt*vy
     rold = sqrt(xold*xold+yold*yold)
     if (rold*dxi>1.e-10) then
        invr = 1./rold
        cold = xold*invr 
        sold = yold*invr
     else
        cold = 1.
        sold = 0.
     end if
     xyold0 = cmplx(cold, sold)
     xmid = xmid + 0.5*xold
     ymid = ymid + 0.5*yold
     rmid=sqrt(xmid*xmid+ymid*ymid)
     if (rmid*dxi>1.e-10) then
        invr = 1./rmid
        cmid = xmid*invr 
        smid = ymid*invr
     else
        cmid = 1.
        smid = 0.
     end if
     xymid0 = cmplx(cmid,smid)
     xold=rold*dxi
     vy = -vx*smid+vy*cmid
     vx = (x-xold)*dx*dti
     zold=z-dtsdz0*vz

     ! --- sets positions relative to grid start
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

     do icell = 1,ncells

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
        wqt(:) = wq*invdtm(:)

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
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sx0(0) = 1. - xint + 1./(4*iixp0+2)*( -xint + xint**2 )
           sx0(1) = 1. - sx0(0)
        else                          ! Standard method, canonical shapes in r
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

        ! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.;sz=0.

        ! --- computes coefficients for quantities centered between nodes
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sx(0+dix) = 1. - xint + 1./(4*iixp+2)*( -xint + xint**2 )
           sx(1+dix) = 1. - sx(0+dix)
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
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)

        ! --- add current contributions
        ! -- NB : the current is later divided by the cylindrical cell volume in applybc_j
        do k=izmin, izmax
           do i=ixmin, ixmax
              ic = iixp0+i
              kc = ikxp0+k

              ! -- Jr
              if(i<ixmax) then
                 sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wr coefficient from esirkepov
                 if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wr along r
                 jr(ic,kc,0) = jr(ic,kc,0) + sdx(i,k)              ! Deposition on the mode m = 0
                 xymid = xymid0 ! Throughout the following loop, xymid takes the value e^{i m theta}
                 do m = 1, nmodes-1                                ! Deposition on the modes m>0
                    jr(ic,kc,m) = jr(ic,kc,m) + 2.*sdx(i,k)*xymid
                    ! The factor 2 comes from the normalization of the modes
                    xymid = xymid*xymid0
                 enddo
              end if

              ! -- Jtheta
              ! Mode m = 0 : similar to the 2D Esirkepov scheme
              jt(ic,kc,0) = jt(ic,kc,0) + wq*vy*invvol/ncells* &
                   ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )
              ! Mode m > 0 : see Davidson et al. JCP 281 (2014)
              xy = xy0 ; xymid = xymid0 ; xyold = xyold0
              ! Throughout the following loop, xy_ takes the value e^{i m theta_}
              do m = 1, nmodes-1
                 jt(ic,kc,m) = jt(ic,kc,m) - 2*im*(ic+xmin*dxi)*wqt(m) * &
                      ( sx0(i)*sz0(k)*(xy-xymid) + sx(i)*sz(k)*(xymid-xyold) )
                 ! The factor 2 comes from the normalization of the modes
                 ! The minus sign comes from the different convention with respect to Davidson et al.
                 xy = xy*xy0 ; xymid = xymid*xymid0 ; xyold = xyold*xyold0
              enddo

              ! -- Jz
              if(k<izmax) then
                 sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov
                 if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
                 jz(ic,kc,0) = jz(ic,kc,0) + sdz(i,k)              ! Deposition on the mode m=0
                 xymid = xymid0 ! Throughout the following loop, xymid takes the value e^{i m theta}
                 do m = 1, nmodes-1                                ! Deposition on the modes m>0
                    jz(ic,kc,m) = jz(ic,kc,m) + 2.*sdz(i,k)*xymid
                    ! The factor 2 comes from the normalization of the modes
                    xymid = xymid*xymid0
                 enddo
              end if
           end do
        end do

     end do

  end do

  deallocate(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)

  return
end subroutine pxr_depose_jrjtjz_esirkepov_n_2d_circ


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
! ESIRKEPOV_CURRENT_DEPOSITION_3D.F90
!
! Developers
! Henri Vincenti, ! Mathieu Lobet
!
! Description:
! This file contains subroutines for Esirkepov current deposition in 3D.
!
! List of subroutines:
!
! Esirkepov:
! - depose_jxjyjz_esirkepov_1_1_1
! - depose_jxjyjz_esirkepov_2_2_2
! - depose_jxjyjz_esirkepov_3_3_3
! - warp_depose_jxjyjz_esirkepov_n
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Esirkepov scalar current deposition algorithm at order 1 in x, y, z (nox=noy=noz=1)
!>
!> @detail
!> This function is not vectorized
!>
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!>
!> @date
!> Revision 10/09/2016
!>
!> @param[inout] jx x-current component (3D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 3 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jx array (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!> (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 3 integers)
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv particle Lorentz factor arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nox, noy, noz interpolation order
!> @param[in] l_particles_weight use the particle weigth
!> @param[in] l4symtry
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz)    !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE

  ! Input/output parameters
  INTEGER(idp)             :: np
  INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), jy_nguard(3), jy_nvalid(3), &
  jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1,                                          &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1,                                          &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,                                          &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                :: q, dt, dx, dy, dz, xmin, ymin, zmin

  ! Internal parameters
  REAL(num)                                :: dxi, dyi, dzi
  REAL(num)                                :: xold, yold, zold, x, y, z, wq
  REAL(num), PARAMETER :: clghtisq=1.0_num/clight**2
  REAL(num)                                :: invvol, invdtdx, invdtdy, invdtdz
  REAL(num)                                :: dtsdx0, dtsdy0, dtsdz0
  REAL(num), PARAMETER                 :: onesixth=1.0_num/6.0_num,                   &
  twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(4) :: sx(-1:2), sx0(-1:2)
  REAL(num), DIMENSION(4) :: sy(-1:2), sy0(-1:2)
  REAL(num), DIMENSION(4) :: sz(-1:2), sz0(-1:2)
  INTEGER(idp)                         :: iixp0, ijxp0, ikxp0
  INTEGER(idp)                         :: iixp, ijxp, ikxp
  INTEGER(idp)                         :: ip, i, j, k, ic, jc, kc,     &
  ixmin, ixmax, iymin, iymax, izmin, izmax
  REAL(num), PARAMETER :: onethird = 1.0_num/3.0_num
  REAL(num)                 :: sdxim1,            sdxi
  REAL(num), DIMENSION(4)   :: sdyjm1(-1:2),      sdyj(-1:2)
  REAL(num), DIMENSION(4,4) :: sdzkm1(-1:2,-1:2), sdzk(-1:2,-1:2)

  ! PARAMETER INIT
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dy*dz)
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)
  dtsdz0 = dt*dzi

  sdxi=0.0_num
  sdxim1=0.0_num
  sdyj=0.0_num
  sdyjm1=0.0_num
  sdzk = 0.0_num
  sdzkm1 = 0.0_num
!$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2), sdxi, sdxim1, &
!$acc&                         sx0(-1:2), sy0(-1:2), sz0(-1:2), &
!$acc&                         sdyj(-1:2), sdyjm1(-1:2), &
!$acc&                         sdzk(-1:2,-1:2), sdzkm1(-1:2,-1:2) )
  DO ip=1, np
    sx = 0.0_num
    sy = 0.0_num
    sz = 0.0_num
    sx0 = 0.0_num
    sy0 = 0.0_num
    sz0 = 0.0_num
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes old position in grid units
    xold=x-dtsdx0*uxp(ip)*gaminv(ip)
    yold=y-dtsdy0*uyp(ip)*gaminv(ip)
    zold=z-dtsdz0*uzp(ip)*gaminv(ip)
    ! --- computes particles weights
    wq=q*w(ip)
    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ijxp0=floor(y)
    ikxp0=floor(z)
    ! --- computes coefficients for node centered quantities
    ! --- x-iixp0: distance between particle and node for current positions
    sx0( 0) = 1.0_num-(x-iixp0)
    sx0( 1) = x-iixp0
    sy0( 0) = 1.0_num-(y-ijxp0)
    sy0( 1) = y-ijxp0
    sz0( 0) = 1.0_num-(z-ikxp0)
    sz0( 1) = z-ikxp0
    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ijxp=floor(yold)
    ikxp=floor(zold)
    ! --- computes coefficients for quantities centered between nodes
    ! --- iixp-iixp0: node separation between old and current positions 
    sx( 0+iixp-iixp0) = 1.0_num-(xold-iixp)
    sx( 1+iixp-iixp0) = xold-iixp
    sy( 0+ijxp-ijxp0) = 1.0_num-(yold-ijxp)
    sy( 1+ijxp-ijxp0) = yold-ijxp
    sz( 0+ikxp-ikxp0) = 1.0_num-(zold-ikxp)
    sz( 1+ikxp-ikxp0) = zold-ikxp
    ! --- computes min/max positions of current contributions
    ixmin = min(0_idp, iixp-iixp0)
    ixmax = max(0_idp, iixp-iixp0)+1
    iymin = min(0_idp, ijxp-ijxp0)
    iymax = max(0_idp, ijxp-ijxp0)+1
    izmin = min(0_idp, ikxp-ikxp0)
    izmax = max(0_idp, ikxp-ikxp0)+1
    ! --- add current contributions
    DO k=izmin, izmax
      DO j=iymin, iymax
        DO i=ixmin, ixmax
          ic = iixp0+i
          jc = ijxp0+j
          kc = ikxp0+k
          IF(i<ixmax) THEN
            sdxi  = wq*invdtdx*(sx(i)-sx0(i))*((sy0(j)+0.5_num*(sy(j)-sy0(j)))*sz0(k) +              &
            (0.5_num*sy0(j)+onethird*(sy(j)-sy0(j)))*(sz(k)-sz0(k)))
            IF (i>ixmin) sdxi = sdxi + sdxim1
            !$acc atomic update
            jx(ic, jc, kc) = jx(ic, jc, kc) + sdxi
          END IF
          IF(j<iymax) THEN
            sdyj(i) = wq*invdtdy*(sy(j)-sy0(j))*((sz0(k)+0.5_num*(sz(k)-sz0(k)))*sx0(i) +              &
            (0.5_num*sz0(k)+onethird*(sz(k)-sz0(k)))*(sx(i)-sx0(i)))
            IF (j>iymin) sdyj(i) = sdyj(i)+sdyjm1(i)
            !$acc atomic update
            jy(ic, jc, kc) = jy(ic, jc, kc) + sdyj(i)
          END IF
          IF(k<izmax) THEN
            sdzk(i, j)  = wq*invdtdz*(sz(k)-sz0(k))*((sx0(i)+0.5_num*(sx(i)-sx0(i)))*sy0(j) +              &
            (0.5_num*sx0(i)+onethird*(sx(i)-sx0(i)))*(sy(j)-sy0(j)))
            IF (k>izmin) sdzk(i, j)=sdzk(i, j)+sdzkm1(i, j)
            !$acc atomic update
            jz(ic, jc, kc) = jz(ic, jc, kc) + sdzk(i, j)
          END IF
          ! REAL
          sdxim1 = sdxi
        END DO
        ! array of 4 REAL
        sdyjm1 = sdyj
      END DO
      ! array of 4 x 4 REAL
      sdzkm1 = sdzk
    END DO

  ENDDO
!$acc end loop
!$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_1_1_1

#if defined (DEV)
! ________________________________________________________________________________________
!> @brief
!> Esirkepov current deposition optimized at order 1
!
!> @details
!> This function gives slightly better performances with AVX512 vector registers.
!> We can expect 30% speedup on KNL however performances are bad with small vector
!> registers.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @warning
!> PROBLEM: CORRECTION REQUIRED
!> DO NOT USE
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_vecHV_1_1_1(jx, jy, jz, np, xp, yp, zp, uxp, uyp,  &
  uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard,    &
  nzguard, l_particles_weight, l4symtry)
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, lp, num
  !USE precomputed
  IMPLICIT NONE
  INTEGER(idp)             :: np, nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT) ::                                                        &
  jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                :: q, dt, dx, dy, dz, xmin, ymin, zmin
  ! Useless here but need to be passed in argument to match func_order arguments
  LOGICAL(lp)              :: l_particles_weight, l4symtry
  REAL(num)                :: xint, yint, zint
  REAL(num)                :: oxint, oyint, ozint, xintsq, yintsq, zintsq, oxintsq,   &
  oyintsq, ozintsq
  REAL(num)                :: x, y, z, xmid, ymid, zmid
  REAL(num)                :: ww, wwx, wwy, wwz
  REAL(num), PARAMETER     :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER     :: twothird=2.0_num/3.0_num
  REAL(num), PARAMETER     :: onethird=1.0_num/3.0_num
  INTEGER(isp)             :: j, k, l, j0, k0, l0, ip, NCELLS, ic, ix, iy, iz
  INTEGER(isp)             :: nnx, nnxy, ngridx, ngridy, n, nn, nv
  INTEGER(isp)             :: moffjx(1:8), moffjy(1:8), moffjz(1:8)

  INTEGER(isp), DIMENSION(LVEC, 3) :: ICELL
  REAL(num), DIMENSION(LVEC, 48)   :: sdx, sdy, sdz
  REAL(num)                :: vx, vy, vz
  REAL(num)                :: wwwx(LVEC, 16), wwwy(LVEC, 16), wwwz(LVEC, 16), wq
  REAL(num)                :: sx1(LVEC), sx2(LVEC), sx3(LVEC), sx4(LVEC)
  REAL(num)                :: sx01(LVEC), sx02(LVEC), sx03(LVEC), sx04(LVEC)
  REAL(num)                :: sy01, sy02, sy03, sy04, sz01, sz02, sz03, sz04
  REAL(num), DIMENSION(4)  :: szz, zdec, h1, h11, h12, sgn
  REAL(num)                :: wx1, wx2, wy1, wy2, wz1, wz2
  INTEGER(isp)             :: orig, ncxy, ncx, ncy, ncz, ngx, ngxy, igrid, jorig,     &
  korig, lorig
  REAL(num), DIMENSION(:), ALLOCATABLE :: sx, sx0, dsx
  REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
  REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz

  INTEGER(isp)             :: dix, diy, diz
  INTEGER(isp)             :: iixp0, ijxp0, ikxp0, iixp, ijxp, ikxp
  INTEGER(isp)             :: iixporig, ijxporig, ikxporig
  INTEGER(isp)             :: errcode
  REAL(num)                :: invdtdx, invdtdy, invdtdz
  REAL(num)                :: wqx, wqy, wqz
  REAL(num)                :: xold, yold, zold, oldsum
  REAL(num)                :: dxi, dyi, dzi
  REAL(num)                :: dtsdx0, dtsdy0, dtsdz0

  ! ___________________________________________
  ! Computation of the parameters

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)

  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard
  ncx=nx+1+2*nxguard
  ncy=ny+1+2*nyguard
  ncz=nz+1+2*nzguard
  NCELLS=ncx*ncy*ncz
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
  ALLOCATE(sx(-1:2), sx0(-1:2), dsx(-1:2))
  ALLOCATE(sy(-1:2), sy0(-1:2), dsy(-1:2))
  ALLOCATE(sz(-1:2), sz0(-1:2), dsz(-1:2))
  jxcells=0.0_num
  jycells=0.0_num
  jzcells=0.0_num
  nnx = ngridx
  nnxy = ngridx*ngridy
  iixporig=-nxguard
  ijxporig=-nyguard
  ikxporig=-nzguard
  !orig=iixporig+nxguard+nnx*(ijxporig+nyguard)+(ikxporig+nzguard)*nnxy
  orig=(nxguard+iixporig) + (nyguard+ijxporig)*nnx + (nzguard+ikxporig)*nnxy

  ngx=(ngridx-ncx)
  ncxy=ncx*ncy
  ngxy=(ngridx*ngridy-ncxy)

  moffjx = (/0_isp, nnx, 2_isp*nnx, 3_isp*nnx, nnxy, nnx+nnxy, 2_isp*nnx+nnxy,        &
  3_isp*nnx+nnxy/)
  moffjy = (/0_isp, 1_isp, 2_isp, 3_isp, nnxy, 1_isp+nnxy, 2_isp+nnxy, 3_isp+nnxy/)
  moffjz = (/0_isp, 1_isp, 2_isp, 3_isp, nnx, 1_isp+nnx, 2_isp+nnx, 3_isp+nnx/)

  h1=(/1_num, 0_num, 1_num, 0_num/); sgn=(/1_num, -1_num, 1_num, -1_num/)
  h11=(/0_num, 1_num, 1_num, 0_num/); h12=(/1_num, 0_num, 0_num, 1_num/)

  sx0 = 0._num
  sy0 = 0._num
  sz0 = 0._num

  !#if DEBUG==1
  !        print*, 'Compute weights'
  !#endif

  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC

#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, ICELL)
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
    DO n=1, MIN(LVEC, np-ip+1)
      nn=ip+n-1
      ! --- computes position in  grid units at (n+1)
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! --- Computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)

      ! --- computes old position in grid units
      xold=x-dtsdx0*vx
      yold=y-dtsdy0*vy
      zold=z-dtsdz0*vz

      ! --- computes particles weights
      wq=q*w(ip)
      wqx = wq*invdtdx
      wqy = wq*invdtdy
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      iixp0=floor(x)
      ijxp0=floor(y)
      ikxp0=floor(z)
      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      yint=y-ijxp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
      sx0(0) = 1.0_num-xint
      sx0(1) = xint
      sy0(0) = 1.0_num-yint
      sy0(1) = yint
      sz0(0) = 1.0_num-zint
      sz0(1) = zint

      ! --- finds node of cell containing particles for old positions
      iixp=floor(xold)
      ijxp=floor(yold)
      ikxp=floor(zold)

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      yint = yold-ijxp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diy = ijxp-ijxp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients
      ! --- (needed because of different dix and diz for each particle)
      sx(-1)=0.0_num
      sx(0)=0.0_num
      sx(1)=0.0_num
      sx(2)=0.0_num

      sy(-1)=0.0_num
      sy(0)=0.0_num
      sy(1)=0.0_num
      sy(2)=0.0_num

      sz(-1)=0.0_num
      sz(0)=0.0_num
      sz(1)=0.0_num
      sz(2)=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      !#if DEBUG==1
      !  if ((dix>1).or.(diy>1).or.(diz>1)) then
      !    print*, '', dix, diy, diz
      !  end if
      !#endif
      sx( 0+dix) = 1.0_num-xint
      sx( 1+dix) = xint
      sy( 0+diy) = 1.0_num-yint
      sy( 1+diy) = yint
      sz( 0+diz) = 1.0_num-zint
      sz( 1+diz) = zint

      ! --- computes coefficients difference
      dsx(-1) = sx(-1) - sx0(-1)
      dsx(0) = sx(0) - sx0(0)
      dsx(1) = sx(1) - sx0(1)
      dsx(2) = sx(2) - sx0(2)

      dsy(-1) = sy(-1) - sy0(-1)
      dsy(0) = sy(0) - sy0(0)
      dsy(1) = sy(1) - sy0(1)
      dsy(2) = sy(2) - sy0(2)

      dsz(-1) = sz(-1) - sz0(-1)
      dsz(0) = sz(0) - sz0(0)
      dsz(1) = sz(1) - sz0(1)
      dsz(2) = sz(2) - sz0(2)

      ! Icell like in the previous function
      !ICELL(n, 1)=1+(iixp0-iixporig)+(ijxp0-ijxporig)*ncx+(ikxp0-ikxporig)*ncxy
      ! With the shift
      ICELL(n, 1)=(iixp0-iixporig)+(ijxp0-ijxporig-1)*ncx+(ikxp0-ikxporig-1)*ncxy
      !ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      !ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! Weight
      sdx(n, 1)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 2)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 3)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 4)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 5)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 6)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 7)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 8)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 9)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 10)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 11)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 12)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 13)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 14)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 15)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 16)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 17)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 17)=sdx(n, 17)+sdx(n, 1)
      sdx(n, 18)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 18)=sdx(n, 18)+sdx(n, 2)
      sdx(n, 19)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 19)=sdx(n, 19)+sdx(n, 3)
      sdx(n, 20)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 20)=sdx(n, 20)+sdx(n, 4)
      sdx(n, 21)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 21)=sdx(n, 21)+sdx(n, 5)
      sdx(n, 22)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 22)=sdx(n, 22)+sdx(n, 6)
      sdx(n, 23)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 23)=sdx(n, 23)+sdx(n, 7)
      sdx(n, 24)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 24)=sdx(n, 24)+sdx(n, 8)
      sdx(n, 25)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 25)=sdx(n, 25)+sdx(n, 9)
      sdx(n, 26)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 26)=sdx(n, 26)+sdx(n, 10)
      sdx(n, 27)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 27)=sdx(n, 27)+sdx(n, 11)
      sdx(n, 28)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 28)=sdx(n, 28)+sdx(n, 12)
      sdx(n, 29)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 29)=sdx(n, 29)+sdx(n, 13)
      sdx(n, 30)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 30)=sdx(n, 30)+sdx(n, 14)
      sdx(n, 31)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 31)=sdx(n, 31)+sdx(n, 15)
      sdx(n, 32)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 32)=sdx(n, 32)+sdx(n, 16)
      sdx(n, 33)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 33)=sdx(n, 33)+sdx(n, 17)
      sdx(n, 34)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 34)=sdx(n, 34)+sdx(n, 18)
      sdx(n, 35)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 35)=sdx(n, 35)+sdx(n, 19)
      sdx(n, 36)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 36)=sdx(n, 36)+sdx(n, 20)
      sdx(n, 37)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 37)=sdx(n, 37)+sdx(n, 21)
      sdx(n, 38)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 38)=sdx(n, 38)+sdx(n, 22)
      sdx(n, 39)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 39)=sdx(n, 39)+sdx(n, 23)
      sdx(n, 40)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 40)=sdx(n, 40)+sdx(n, 24)
      sdx(n, 41)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 41)=sdx(n, 41)+sdx(n, 25)
      sdx(n, 42)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 42)=sdx(n, 42)+sdx(n, 26)
      sdx(n, 43)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 43)=sdx(n, 43)+sdx(n, 27)
      sdx(n, 44)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 44)=sdx(n, 44)+sdx(n, 28)
      sdx(n, 45)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 45)=sdx(n, 45)+sdx(n, 29)
      sdx(n, 46)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 46)=sdx(n, 46)+sdx(n, 30)
      sdx(n, 47)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 47)=sdx(n, 47)+sdx(n, 31)
      sdx(n, 48)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 48)=sdx(n, 48)+sdx(n, 32)

      sdy(n, 1)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 2)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 3)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 4)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 5)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 6)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 7)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 8)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 9)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 10)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 11)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 12)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 13)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 14)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 15)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 16)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 17)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 17)=sdy(n, 17)+sdy(n, 1)
      sdy(n, 18)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 18)=sdy(n, 18)+sdy(n, 2)
      sdy(n, 19)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 19)=sdy(n, 19)+sdy(n, 3)
      sdy(n, 20)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 20)=sdy(n, 20)+sdy(n, 4)
      sdy(n, 21)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 21)=sdy(n, 21)+sdy(n, 5)
      sdy(n, 22)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 22)=sdy(n, 22)+sdy(n, 6)
      sdy(n, 23)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 23)=sdy(n, 23)+sdy(n, 7)
      sdy(n, 24)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 24)=sdy(n, 24)+sdy(n, 8)
      sdy(n, 25)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 25)=sdy(n, 25)+sdy(n, 9)
      sdy(n, 26)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 26)=sdy(n, 26)+sdy(n, 10)
      sdy(n, 27)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 27)=sdy(n, 27)+sdy(n, 11)
      sdy(n, 28)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 28)=sdy(n, 28)+sdy(n, 12)
      sdy(n, 29)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 29)=sdy(n, 29)+sdy(n, 13)
      sdy(n, 30)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 30)=sdy(n, 30)+sdy(n, 14)
      sdy(n, 31)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 31)=sdy(n, 31)+sdy(n, 15)
      sdy(n, 32)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 32)=sdy(n, 32)+sdy(n, 16)
      sdy(n, 33)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 33)=sdy(n, 33)+sdy(n, 17)
      sdy(n, 34)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 34)=sdy(n, 34)+sdy(n, 18)
      sdy(n, 35)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 35)=sdy(n, 35)+sdy(n, 19)
      sdy(n, 36)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 36)=sdy(n, 36)+sdy(n, 20)
      sdy(n, 37)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 37)=sdy(n, 37)+sdy(n, 21)
      sdy(n, 38)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 38)=sdy(n, 38)+sdy(n, 22)
      sdy(n, 39)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 39)=sdy(n, 39)+sdy(n, 23)
      sdy(n, 40)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 40)=sdy(n, 40)+sdy(n, 24)
      sdy(n, 41)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 41)=sdy(n, 41)+sdy(n, 25)
      sdy(n, 42)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 42)=sdy(n, 42)+sdy(n, 26)
      sdy(n, 43)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 43)=sdy(n, 43)+sdy(n, 27)
      sdy(n, 44)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 44)=sdy(n, 44)+sdy(n, 28)
      sdy(n, 45)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 45)=sdy(n, 45)+sdy(n, 29)
      sdy(n, 46)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 46)=sdy(n, 46)+sdy(n, 30)
      sdy(n, 47)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 47)=sdy(n, 47)+sdy(n, 31)
      sdy(n, 48)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 48)=sdy(n, 48)+sdy(n, 32)

      sdz(n, 1)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 2)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 3)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 4)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 5)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 6)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 7)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 8)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 9)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 10)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 11)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 12)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 13)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 14)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 15)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 16)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 17)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 17)=sdz(n, 17)+sdz(n, 1)
      sdz(n, 18)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 18)=sdz(n, 18)+sdz(n, 2)
      sdz(n, 19)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 19)=sdz(n, 19)+sdz(n, 3)
      sdz(n, 20)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 20)=sdz(n, 20)+sdz(n, 4)
      sdz(n, 21)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 21)=sdz(n, 21)+sdz(n, 5)
      sdz(n, 22)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 22)=sdz(n, 22)+sdz(n, 6)
      sdz(n, 23)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 23)=sdz(n, 23)+sdz(n, 7)
      sdz(n, 24)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 24)=sdz(n, 24)+sdz(n, 8)
      sdz(n, 25)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 25)=sdz(n, 25)+sdz(n, 9)
      sdz(n, 26)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 26)=sdz(n, 26)+sdz(n, 10)
      sdz(n, 27)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 27)=sdz(n, 27)+sdz(n, 11)
      sdz(n, 28)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 28)=sdz(n, 28)+sdz(n, 12)
      sdz(n, 29)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 29)=sdz(n, 29)+sdz(n, 13)
      sdz(n, 30)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 30)=sdz(n, 30)+sdz(n, 14)
      sdz(n, 31)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 31)=sdz(n, 31)+sdz(n, 15)
      sdz(n, 32)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 32)=sdz(n, 32)+sdz(n, 16)
      sdz(n, 33)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 33)=sdz(n, 33)+sdz(n, 17)
      sdz(n, 34)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 34)=sdz(n, 34)+sdz(n, 18)
      sdz(n, 35)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 35)=sdz(n, 35)+sdz(n, 19)
      sdz(n, 36)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 36)=sdz(n, 36)+sdz(n, 20)
      sdz(n, 37)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 37)=sdz(n, 37)+sdz(n, 21)
      sdz(n, 38)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 38)=sdz(n, 38)+sdz(n, 22)
      sdz(n, 39)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 39)=sdz(n, 39)+sdz(n, 23)
      sdz(n, 40)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 40)=sdz(n, 40)+sdz(n, 24)
      sdz(n, 41)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 41)=sdz(n, 41)+sdz(n, 25)
      sdz(n, 42)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 42)=sdz(n, 42)+sdz(n, 26)
      sdz(n, 43)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 43)=sdz(n, 43)+sdz(n, 27)
      sdz(n, 44)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 44)=sdz(n, 44)+sdz(n, 28)
      sdz(n, 45)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 45)=sdz(n, 45)+sdz(n, 29)
      sdz(n, 46)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 46)=sdz(n, 46)+sdz(n, 30)
      sdz(n, 47)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 47)=sdz(n, 47)+sdz(n, 31)
      sdz(n, 48)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 48)=sdz(n, 48)+sdz(n, 32)
    END DO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! Add weights to nearest vertices
    DO n=1, MIN(LVEC, np-ip+1)
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, jxcells, jycells, jzcells)
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
      DO nv=1, 8
        ! --- JX
        ! ICELL = (-1, -1, -1)
        ! Loop on (i=-1, j=-1, k=-1)
        jxcells(nv, ICELL(n, 1)) = jxcells(nv, ICELL(n, 1)) + sdx(n, nv)
        ! Loop on (i=-1, j=-1, k=1)
        jxcells(nv, ICELL(n, 1)+2*ncxy)   = jxcells(nv, ICELL(n, 1)+2*ncxy) + sdx(n,  &
        8+nv)
        !Loop on (i=0, j, k=-1)
        jxcells(nv, ICELL(n, 1)+1) = jxcells(nv, ICELL(n, 1)+1) + sdx(n, 16+nv)
        !Loop on (i=0, j, k=1)
        jxcells(nv, ICELL(n, 1)+1+2*ncxy) = jxcells(nv, ICELL(n, 1)+1+2*ncxy) +       &
        sdx(n, 24+nv)
        ! Loop on (i=1, j, k=-1)
        jxcells(nv, ICELL(n, 1)+2) = jxcells(nv, ICELL(n, 1)+2) + sdx(n, 32+nv)
        ! Loop on (i=1, j, k=1)
        jxcells(nv, ICELL(n, 1)+2+2*ncxy) = jxcells(nv, ICELL(n, 1)+2+2*ncxy) +       &
        sdx(n, 40+nv)

        ! --- JY
        ! Loop on (i=-1, j=-1, k=-1)
        jycells(nv, ICELL(n, 1))           = jycells(nv, ICELL(n, 1)) + sdy(n, nv)
        ! Loop on (i=-1, j=-1, k=1)
        jycells(nv, ICELL(n, 1)+2*ncxy)    = jycells(nv, ICELL(n, 1)+2*ncxy) + sdy(n, &
        nv+8)
        !Loop on (i=-1, j=0, k=-1)
        jycells(nv, ICELL(n, 1)+ncx)       = jycells(nv, ICELL(n, 1)+ncx) + sdy(n,    &
        nv+16)
        !Loop on (i=-1, j=0, k=1)
        jycells(nv, ICELL(n, 1)+ncx+2*ncxy) = jycells(nv, ICELL(n, 1)+ncx+2*ncxy) +   &
        sdy(n, nv+24)
        ! Loop on (i=-1, j=1, k=-1)
        jycells(nv, ICELL(n, 1)+2*ncx)      = jycells(nv, ICELL(n, 1)+2*ncx)  +       &
        sdy(n, nv+32)
        ! Loop on (i=-1, j=1, k=1)
        jycells(nv, ICELL(n, 1)+2*ncx+2*ncxy)=jycells(nv, ICELL(n, 1)+2*ncx+2*ncxy)+  &
        sdy(n, nv+40)
        ! --- JZ
        ! Loop on (i=-1, j=-1, k=-1)
        jzcells(nv, ICELL(n, 1))     = jzcells(nv, ICELL(n, 1)) + sdz(n, nv)
        ! Loop on (i=-1, j=1, k=-1)
        jzcells(nv, ICELL(n, 1)+2*ncx) = jzcells(nv, ICELL(n, 1)+2*ncx) + sdz(n,      &
        nv+8)
        !Loop on (i=-1, j=-1, k=0)
        jzcells(nv, ICELL(n, 1)+ncxy) = jzcells(nv, ICELL(n, 1)+ncxy)+ sdz(n, nv+16)
        !Loop on (i=-1, j=1, k=0)
        jzcells(nv, ICELL(n, 1)+2*ncx+ncxy) = jzcells(nv, ICELL(n, 1)+2*ncx+ncxy) +   &
        sdz(n, nv+24)
        ! Loop on (i=-1, j, k=1)
        jzcells(nv, ICELL(n, 1)+2*ncxy) = jzcells(nv, ICELL(n, 1)+2*ncxy) + sdz(n,    &
        nv+32)
        ! Loop on (i=-1, j, k=1)
        jzcells(nv, ICELL(n, 1)+2*ncx+2*ncxy)=jzcells(nv, ICELL(n,                    &
        1)+2*ncx+2*ncxy)+sdz(n, nv+40)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
  DO iz=1, ncz-2
    DO iy=1, ncy-2
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix=1, ncx-2!! VECTOR (take ncx multiple of vector length)
        ic=ix+(iy-1)*ncx+(iz-1)*ncxy
        !igrid=orig+ic+(iy-1)*ngx+(iz-1)*ngxy
        igrid =orig+ix+(iy-1)*nnx+(iz-1)*nnxy

        ! jx
        jx(igrid+moffjx(1))=jx(igrid+moffjx(1))+jxcells(1, ic)
        jx(igrid+moffjx(2))=jx(igrid+moffjx(2))+jxcells(2, ic)
        jx(igrid+moffjx(3))=jx(igrid+moffjx(3))+jxcells(3, ic)
        jx(igrid+moffjx(4))=jx(igrid+moffjx(4))+jxcells(4, ic)
        jx(igrid+moffjx(5))=jx(igrid+moffjx(5))+jxcells(5, ic)
        jx(igrid+moffjx(6))=jx(igrid+moffjx(6))+jxcells(6, ic)
        jx(igrid+moffjx(7))=jx(igrid+moffjx(7))+jxcells(7, ic)
        jx(igrid+moffjx(8))=jx(igrid+moffjx(8))+jxcells(8, ic)
        ! jy
        jy(igrid+moffjy(1))=jy(igrid+moffjy(1))+jycells(1, ic)
        jy(igrid+moffjy(2))=jy(igrid+moffjy(2))+jycells(2, ic)
        jy(igrid+moffjy(3))=jy(igrid+moffjy(3))+jycells(3, ic)
        jy(igrid+moffjy(4))=jy(igrid+moffjy(4))+jycells(4, ic)
        jy(igrid+moffjy(5))=jy(igrid+moffjy(5))+jycells(5, ic)
        jy(igrid+moffjy(6))=jy(igrid+moffjy(6))+jycells(6, ic)
        jy(igrid+moffjy(7))=jy(igrid+moffjy(7))+jycells(7, ic)
        jy(igrid+moffjy(8))=jy(igrid+moffjy(8))+jycells(8, ic)
        ! jz
        jz(igrid+moffjz(1))=jz(igrid+moffjz(1))+jzcells(1, ic)
        jz(igrid+moffjz(2))=jz(igrid+moffjz(2))+jzcells(2, ic)
        jz(igrid+moffjz(3))=jz(igrid+moffjz(3))+jzcells(3, ic)
        jz(igrid+moffjz(4))=jz(igrid+moffjz(4))+jzcells(4, ic)
        jz(igrid+moffjz(5))=jz(igrid+moffjz(5))+jzcells(5, ic)
        jz(igrid+moffjz(6))=jz(igrid+moffjz(6))+jzcells(6, ic)
        jz(igrid+moffjz(7))=jz(igrid+moffjz(7))+jzcells(7, ic)
        jz(igrid+moffjz(8))=jz(igrid+moffjz(8))+jzcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(jxcells, jycells, jzcells)
  RETURN

END SUBROUTINE depose_jxjyjz_esirkepov_vecHV_1_1_1
#endif


! __ Developer zone _________
#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Esirkepov current deposition at order 1
!
!> @details
!> This function is a test of optimization but does not give good performance
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @warning
!> DO NOT USE
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_vecHVv2_1_1_1(jx, jy, jz, np, xp, yp, zp, uxp,     &
  uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard,        &
  nyguard, nzguard, nox, noy, noz, l_particles_weight, l4symtry)
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, lp, num
  USE precomputed, ONLY: dtsdx0, dtsdy0, dtsdz0, dxi, dyi, dzi
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard, nox, noy, noz
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                      &
  -nzguard:nz+nzguard), INTENT(IN OUT) :: jx, jy, jz

  REAL(num), DIMENSION(:, :), ALLOCATABLE:: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq, zintsq,         &
  oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, xmid, ymid, zmid
  REAL(num) ::   ww, wwx, wwy, wwz
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: i, j, k, l, j0, k0, l0, ip, NCELLS, ic, ix, iy, iz
  INTEGER(isp) :: i2, j2, k2, jc, kc
  INTEGER(idp) :: nnx, nnxy, ngridx, ngridy, n, nn, nv
  INTEGER(idp) :: moffjx(1:8), moffjy(1:8), moffjz(1:8)

  INTEGER(idp), DIMENSION(LVEC, 3) :: ICELL
  REAL(num), DIMENSION(LVEC) :: vx, vy, vz
  REAL(num) ::  wwwx(LVEC, 16), wwwy(LVEC, 16), wwwz(LVEC, 16), wq
  REAL(num) :: sx1(LVEC), sx2(LVEC), sx3(LVEC), sx4(LVEC)
  REAL(num) :: sx01(LVEC), sx02(LVEC), sx03(LVEC), sx04(LVEC)
  REAL(num) :: sy1, sy2, sy3, sy4, sz1, sz2, sz3, sz4
  REAL(num) :: sy01, sy02, sy03, sy04, sz01, sz02, sz03, sz04
  REAL(num), DIMENSION(4) :: szz, zdec, h1, h11, h12, sgn
  REAL(num):: wx1, wx2, wy1, wy2, wz1, wz2
  INTEGER(idp) :: orig, ncxy, ncx, ncy, ncz, ngx, ngxy, igrid, jorig, korig, lorig
  ! Useless here but need to be passed in argument to match func_order arguments
  LOGICAL(lp)  :: l_particles_weight, l4symtry
  REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx
  REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
  REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
  REAL(num), PARAMETER :: onethird = 1.0_num/3.0_num

  REAL(num), DIMENSION(LVEC, 48) :: sdx, sdy, sdz

  INTEGER(idp) :: dix, diy, diz
  INTEGER(idp) :: iixp0(LVEC), ijxp0(LVEC), ikxp0(LVEC), iixp, ijxp, ikxp
  INTEGER(idp) :: iixporig, ijxporig, ikxporig
  REAL(num) :: invdtdx, invdtdy, invdtdz
  REAL(num) :: wqx, wqy, wqz
  REAL(num) :: xold, yold, zold


  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard
  ncx=nx+1+2*nxguard
  ncy=ny+1+2*nyguard
  ncz=nz+1+2*nzguard
  NCELLS=ncx*ncy*ncz
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
  ALLOCATE(sx(-1:2), sx0(-1:2), dsx(-1:2))
  ALLOCATE(sy(-1:2), sy0(-1:2), dsy(-1:2))
  ALLOCATE(sz(-1:2), sz0(-1:2), dsz(-1:2))
  jxcells=0.0_num
  jycells=0.0_num
  jzcells=0.0_num
  nnx = ngridx
  nnxy = ngridx*ngridy
  iixporig=-nxguard
  ijxporig=-nyguard
  ikxporig=-nzguard
  orig=(nxguard+iixporig) + (nyguard+ijxporig)*nnx + (nzguard+ikxporig)*nnxy

  ngx=(ngridx-ncx)
  ncxy=ncx*ncy
  ngxy=(ngridx*ngridy-ncxy)

  moffjx = (/0_idp, nnx, 2_idp*nnx, 3_idp*nnx, nnxy, nnx+nnxy, 2_idp*nnx+nnxy,        &
  3_idp*nnx+nnxy/)
  moffjy = (/0_idp, 1_idp, 2_idp, 3_idp, nnxy, 1_idp+nnxy, 2_idp+nnxy, 3_idp+nnxy/)
  moffjz = (/0_idp, 1_idp, 2_idp, 3_idp, nnx, 1_idp+nnx, 2_idp+nnx, 3_idp+nnx/)

  h1=(/1_num, 0_num, 1_num, 0_num/); sgn=(/1_num, -1_num, 1_num, -1_num/)
  h11=(/0_num, 1_num, 1_num, 0_num/); h12=(/1_num, 0_num, 0_num, 1_num/)

  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED vx:64, vy:64, vz:64
    !DIR$ ASSUME_ALIGNED sx1:64, sx2:64, sx3:64, sx4:64
    !DIR$ ASSUME_ALIGNED sx01:64, sx02:64, sx03:64, sx04:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, vx, vy, vz)
    !IBM* ALIGN(64, sx1, sx2, sx3, sx4)
    !IBM* ALIGN(64, sx01, sx02, sx03, sx04)
    !IBM* ALIGN(64, ICELL)
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
    DO n=1, MIN(LVEC, np-ip+1)
      nn=ip+n-1
      ! --- computes position in  grid units at (n+1)
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! --- Computes velocity
      vx(n) = uxp(nn)*gaminv(nn)
      vy(n) = uyp(nn)*gaminv(nn)
      vz(n) = uzp(nn)*gaminv(nn)

      ! --- computes old position in grid units
      xold=x-dtsdx0*vx(n)
      yold=y-dtsdy0*vy(n)
      zold=z-dtsdz0*vz(n)

      ! --- computes particles weights
      wq=q*w(ip)
      wqx = wq*invdtdx
      wqy = wq*invdtdy
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      iixp0(n)=floor(x)
      ijxp0(n)=floor(y)
      ikxp0(n)=floor(z)
      ! --- computes distance between particle and node for current positions
      xint=x-iixp0(n)
      yint=y-ijxp0(n)
      zint=z-ikxp0(n)

      ! --- computes coefficients for node centered quantities
      sx0(0) = 1.0_num-xint
      sx0(1) = xint
      sy0(0) = 1.0_num-yint
      sy0(1) = yint
      sz0(0) = 1.0_num-zint
      sz0(1) = zint

      ! --- finds node of cell containing particles for old positions
      iixp=floor(xold)
      ijxp=floor(yold)
      ikxp=floor(zold)

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      yint = yold-ijxp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0(n)
      diy = ijxp-ijxp0(n)
      diz = ikxp-ikxp0(n)

      ! --- zero out coefficients
      ! ---  (needed because of different dix and diz for each particle)
      sx(-1)=0.0_num
      sx(0)=0.0_num
      sx(1)=0.0_num
      sx(2)=0.0_num

      sy(-1)=0.0_num
      sy(0)=0.0_num
      sy(1)=0.0_num
      sy(2)=0.0_num

      sz(-1)=0.0_num
      sz(0)=0.0_num
      sz(1)=0.0_num
      sz(2)=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      sx( 0+dix) = 1.0_num-xint
      sx( 1+dix) = xint
      sy( 0+diy) = 1.0_num-yint
      sy( 1+diy) = yint
      sz( 0+diz) = 1.0_num-zint
      sz( 1+diz) = zint

      ! --- computes coefficients difference
      dsx(-1) = sx(-1) - sx0(-1)
      dsx(0) = sx(0) - sx0(0)
      dsx(1) = sx(1) - sx0(1)
      dsx(2) = sx(2) - sx0(2)

      dsy(-1) = sy(-1) - sy0(-1)
      dsy(0) = sy(0) - sy0(0)
      dsy(1) = sy(1) - sy0(1)
      dsy(2) = sy(2) - sy0(2)

      dsz(-1) = sz(-1) - sz0(-1)
      dsz(0) = sz(0) - sz0(0)
      dsz(1) = sz(1) - sz0(1)
      dsz(2) = sz(2) - sz0(2)

      ! Icell like in the previous function
      ! With the shift
      ICELL(n,                                                                        &
      1)=(iixp0(n)-iixporig)+(ijxp0(n)-ijxporig-1)*ncx+(ikxp0(n)-ikxporig-1)*ncxy

      ! Weight
      sdx(n, 1)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 2)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 2)=sdx(n, 2)+sdx(n, 1)
      sdx(n, 3)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 3)=sdx(n, 3)+sdx(n, 2)
      sdx(n, 4)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 5)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 5)=sdx(n, 5)+sdx(n, 4)
      sdx(n, 6)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 6)=sdx(n, 6)+sdx(n, 5)
      sdx(n, 7)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 8)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 8)=sdx(n, 8)+sdx(n, 7)
      sdx(n, 9)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 9)=sdx(n, 9)+sdx(n, 8)
      sdx(n, 10)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 11)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 11)=sdx(n, 11)+sdx(n, 10)
      sdx(n, 12)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 12)=sdx(n, 12)+sdx(n, 11)
      sdx(n, 13)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 14)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 14)=sdx(n, 14)+sdx(n, 13)
      sdx(n, 15)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 15)=sdx(n, 15)+sdx(n, 14)
      sdx(n, 16)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 17)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 17)=sdx(n, 17)+sdx(n, 16)
      sdx(n, 18)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 18)=sdx(n, 18)+sdx(n, 17)
      sdx(n, 19)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 20)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 20)=sdx(n, 20)+sdx(n, 19)
      sdx(n, 21)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 21)=sdx(n, 21)+sdx(n, 20)
      sdx(n, 22)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 23)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 23)=sdx(n, 23)+sdx(n, 22)
      sdx(n, 24)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 24)=sdx(n, 24)+sdx(n, 23)
      sdx(n, 25)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 26)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 26)=sdx(n, 26)+sdx(n, 25)
      sdx(n, 27)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 27)=sdx(n, 27)+sdx(n, 26)
      sdx(n, 28)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 29)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 29)=sdx(n, 29)+sdx(n, 28)
      sdx(n, 30)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 30)=sdx(n, 30)+sdx(n, 29)
      sdx(n, 31)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 32)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 32)=sdx(n, 32)+sdx(n, 31)
      sdx(n, 33)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 33)=sdx(n, 33)+sdx(n, 32)
      sdx(n, 34)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 35)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 35)=sdx(n, 35)+sdx(n, 34)
      sdx(n, 36)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 36)=sdx(n, 36)+sdx(n, 35)
      sdx(n, 37)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 38)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 38)=sdx(n, 38)+sdx(n, 37)
      sdx(n, 39)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 39)=sdx(n, 39)+sdx(n, 38)
      sdx(n, 40)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 41)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 41)=sdx(n, 41)+sdx(n, 40)
      sdx(n, 42)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 42)=sdx(n, 42)+sdx(n, 41)
      sdx(n, 43)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 44)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 44)=sdx(n, 44)+sdx(n, 43)
      sdx(n, 45)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 45)=sdx(n, 45)+sdx(n, 44)
      sdx(n, 46)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 47)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 47)=sdx(n, 47)+sdx(n, 46)
      sdx(n, 48)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 48)=sdx(n, 48)+sdx(n, 47)

      ! Weight for y
      sdy(n, 1)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 2)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 3)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 4)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 5)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 5)=sdy(n, 5)+sdy(n, 1)
      sdy(n, 6)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                     &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 6)=sdy(n, 6)+sdy(n, 2)
      sdy(n, 7)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                     &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 7)=sdy(n, 7)+sdy(n, 3)
      sdy(n, 8)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                     &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 8)=sdy(n, 8)+sdy(n, 4)
      sdy(n, 9)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 9)=sdy(n, 9)+sdy(n, 5)
      sdy(n, 10)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 10)=sdy(n, 10)+sdy(n, 6)
      sdy(n, 11)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 11)=sdy(n, 11)+sdy(n, 7)
      sdy(n, 12)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 12)=sdy(n, 12)+sdy(n, 8)
      sdy(n, 13)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                    &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 14)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 15)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 16)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 17)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 17)=sdy(n, 17)+sdy(n, 13)
      sdy(n, 18)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 18)=sdy(n, 18)+sdy(n, 14)
      sdy(n, 19)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 19)=sdy(n, 19)+sdy(n, 15)
      sdy(n, 20)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 20)=sdy(n, 20)+sdy(n, 16)
      sdy(n, 21)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 21)=sdy(n, 21)+sdy(n, 17)
      sdy(n, 22)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 22)=sdy(n, 22)+sdy(n, 18)
      sdy(n, 23)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 23)=sdy(n, 23)+sdy(n, 19)
      sdy(n, 24)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 24)=sdy(n, 24)+sdy(n, 20)
      sdy(n, 25)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                    &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 26)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 27)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 28)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 29)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 29)=sdy(n, 29)+sdy(n, 25)
      sdy(n, 30)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 30)=sdy(n, 30)+sdy(n, 26)
      sdy(n, 31)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 31)=sdy(n, 31)+sdy(n, 27)
      sdy(n, 32)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 32)=sdy(n, 32)+sdy(n, 28)
      sdy(n, 33)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 33)=sdy(n, 33)+sdy(n, 29)
      sdy(n, 34)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 34)=sdy(n, 34)+sdy(n, 30)
      sdy(n, 35)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 35)=sdy(n, 35)+sdy(n, 31)
      sdy(n, 36)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 36)=sdy(n, 36)+sdy(n, 32)
      sdy(n, 37)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 38)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 39)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 40)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 41)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 41)=sdy(n, 41)+sdy(n, 37)
      sdy(n, 42)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 42)=sdy(n, 42)+sdy(n, 38)
      sdy(n, 43)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 43)=sdy(n, 43)+sdy(n, 39)
      sdy(n, 44)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 44)=sdy(n, 44)+sdy(n, 40)
      sdy(n, 45)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 45)=sdy(n, 45)+sdy(n, 41)
      sdy(n, 46)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 46)=sdy(n, 46)+sdy(n, 42)
      sdy(n, 47)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 47)=sdy(n, 47)+sdy(n, 43)
      sdy(n, 48)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                      &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 48)=sdy(n, 48)+sdy(n, 44)

      ! Weight for z
      sdz(n, 1)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 2)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 3)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 4)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 5)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 6)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 7)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 8)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 9)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 10)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 11)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 12)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 13)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 14)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 15)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 16)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 17)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 17)=sdz(n, 17)+sdz(n, 1)
      sdz(n, 18)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 18)=sdz(n, 18)+sdz(n, 2)
      sdz(n, 19)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 19)=sdz(n, 19)+sdz(n, 3)
      sdz(n, 20)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 20)=sdz(n, 20)+sdz(n, 4)
      sdz(n, 21)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 21)=sdz(n, 21)+sdz(n, 5)
      sdz(n, 22)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 22)=sdz(n, 22)+sdz(n, 6)
      sdz(n, 23)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 23)=sdz(n, 23)+sdz(n, 7)
      sdz(n, 24)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 24)=sdz(n, 24)+sdz(n, 8)
      sdz(n, 25)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 25)=sdz(n, 25)+sdz(n, 9)
      sdz(n, 26)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 26)=sdz(n, 26)+sdz(n, 10)
      sdz(n, 27)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 27)=sdz(n, 27)+sdz(n, 11)
      sdz(n, 28)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 28)=sdz(n, 28)+sdz(n, 12)
      sdz(n, 29)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 29)=sdz(n, 29)+sdz(n, 13)
      sdz(n, 30)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 30)=sdz(n, 30)+sdz(n, 14)
      sdz(n, 31)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 31)=sdz(n, 31)+sdz(n, 15)
      sdz(n, 32)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 32)=sdz(n, 32)+sdz(n, 16)
      sdz(n, 33)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 33)=sdz(n, 33)+sdz(n, 17)
      sdz(n, 34)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 34)=sdz(n, 34)+sdz(n, 18)
      sdz(n, 35)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 35)=sdz(n, 35)+sdz(n, 19)
      sdz(n, 36)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 36)=sdz(n, 36)+sdz(n, 20)
      sdz(n, 37)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 37)=sdz(n, 37)+sdz(n, 21)
      sdz(n, 38)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 38)=sdz(n, 38)+sdz(n, 22)
      sdz(n, 39)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 39)=sdz(n, 39)+sdz(n, 23)
      sdz(n, 40)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 40)=sdz(n, 40)+sdz(n, 24)
      sdz(n, 41)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 41)=sdz(n, 41)+sdz(n, 25)
      sdz(n, 42)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 42)=sdz(n, 42)+sdz(n, 26)
      sdz(n, 43)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 43)=sdz(n, 43)+sdz(n, 27)
      sdz(n, 44)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 44)=sdz(n, 44)+sdz(n, 28)
      sdz(n, 45)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 45)=sdz(n, 45)+sdz(n, 29)
      sdz(n, 46)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 46)=sdz(n, 46)+sdz(n, 30)
      sdz(n, 47)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 47)=sdz(n, 47)+sdz(n, 31)
      sdz(n, 48)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 48)=sdz(n, 48)+sdz(n, 32)

    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! Add weights to nearest vertices
    DO n=1, MIN(LVEC, np-ip+1)

      ! --- add current contributions
      i2 = 0
      j2 = 0
      k2 = 0
      DO k=-1, 2
        DO j=-1, 2
          DO i=-1, 2
            ic = iixp0(n)+i
            jc = ijxp0(n)+j
            kc = ikxp0(n)+k
            IF(i<2) THEN
              i2 = i2+1
              jx(ic, jc, kc) = jx(ic, jc, kc) + sdx(n, i2)
            END IF
            IF(j<2) THEN
              j2 = j2+1
              jy(ic, jc, kc) = jy(ic, jc, kc) + sdy(n, j2)
            END IF
            IF(k<2) THEN
              k2 = k2 + 1
              jz(ic, jc, kc) = jz(ic, jc, kc) + sdz(n, k2)
            END IF
          END DO
        END DO
      END DO

    END DO
  END DO

  RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_vecHVv2_1_1_1
#endif

! ________________________________________________________________________________________
!> @brief
!> Esirkepov scalar current deposition algorithm at order 2 in x, y, z (nox=noy=noz=2)
!>
!> @detail
!> This function is not vectorized
!>
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!>
!> @date
!> Revision 10/09/2016
!>
!> @param[inout] jx x-current component (3D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 3 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the jx array
!> (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!> (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the jy array
!> (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!>  (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the jz array
!>  (1d array containing 3 integers)s
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv particle Lorentz factor arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nox, noy, noz interpolation order
!> @param[in] l_particles_weight use the particle weigth
!> @param[in] l4symtry
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz)    !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER :: np
  INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), jy_nguard(3), jy_nvalid(3), &
  jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1,                                          &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1,                                          &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,                                          &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint
  REAL(num) :: clghtisq, xold, yold, zold, x, y, z, wq, wqx, wqy, wqz, vx, vy, vz
  REAL(num) :: invvol, invdtdx, invdtdy, invdtdz
  REAL(num) :: xintsq, yintsq, zintsq
  REAL(num) :: dtsdx0, dtsdy0, dtsdz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER :: twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(5)   ::  sx(-2:2),  sy(-2:2),  sz(-2:2)
  REAL(num), DIMENSION(5)   :: sx0(-2:2), sy0(-2:2), sz0(-2:2)
  REAL(num), DIMENSION(5)   :: dsx(-2:2), dsy(-2:2), dsz(-2:2)
  REAL(num)                 :: sdxim1,            sdxi
  REAL(num), DIMENSION(5)   :: sdyjm1(-2:2),      sdyj(-2:2)
  REAL(num), DIMENSION(5,5) :: sdzkm1(-2:2,-2:2), sdzk(-2:2,-2:2)
  INTEGER :: iixp0, ijxp0, ikxp0, iixp, ijxp, ikxp, ip, dix, diy, diz, i, j, k, ic,   &
  jc, kc
  INTEGER :: ixmin, ixmax, iymin, iymax, izmin, izmax

  ! PARAMETER INIT
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dy*dz)
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)
  clghtisq = 1.0_num/clight**2
  dtsdz0 = dt*dzi
  sdxi=0.0_num
  sdxim1=0.0_num
  sdyj=0.0_num
  sdyjm1=0.0_num
  sdzk = 0.0_num
  sdzkm1 = 0.0_num
!$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(-2:2), sy(-2:2), sz(-2:2), sdxi, sdxim1, &
!$acc&                         sx0(-2:2), sy0(-2:2), sz0(-2:2), &
!$acc&                         dsx(-2:2), dsy(-2:2), dsz(-2:2), &
!$acc&                         sdyj(-2:2), sdyjm1(-2:2), &
!$acc&                         sdzk(-2:2,-2:2), sdzkm1(-2:2,-2:2) )
  DO ip=1, np
    sx = 0.0_num
    sy = 0.0_num
    sz = 0.0_num
    sx0 = 0.0_num
    sy0 = 0.0_num
    sz0 = 0.0_num
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqy = wq*invdtdy
    wqz = wq*invdtdz
    ! --- finds node of cell containing particles for current positions
    iixp0=nint(x)
    ijxp0=nint(y)
    ikxp0=nint(z)
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    yint=y-ijxp0
    zint=z-ikxp0
    ! --- computes coefficients for node centered quantities
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2
    yintsq = yint*yint
    sy0(-1) = 0.5_num*(0.5_num-yint)**2
    sy0( 0) = 0.75_num-yintsq
    sy0( 1) = 0.5_num*(0.5_num+yint)**2
    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = 0.75_num-zintsq
    sz0( 1) = 0.5_num*(0.5_num+zint)**2
    ! --- finds node of cell containing particles for old positions
    iixp=nint(xold)
    ijxp=nint(yold)
    ikxp=nint(zold)
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    yint = yold-ijxp
    zint = zold-ikxp
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diy = ijxp-ijxp0
    diz = ikxp-ikxp0
    ! --- zero out coefficients
    ! --- computes coefficients for quantities centered between nodes
    xintsq = xint*xint
    sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
    sx( 0+dix) = 0.75_num-xintsq
    sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
    yintsq = yint*yint
    sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
    sy( 0+diy) = 0.75_num-yintsq
    sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
    zintsq = zint*zint
    sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
    sz( 0+diz) = 0.75_num-zintsq
    sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsy = sy - sy0
    dsz = sz - sz0

    ! --- computes min/max positions of current contributions
    ixmin = min(0, dix)-1
    ixmax = max(0, dix)+1
    iymin = min(0, diy)-1
    iymax = max(0, diy)+1
    izmin = min(0, diz)-1
    izmax = max(0, diz)+1

    ! --- add current contributions
    DO k=izmin, izmax
      DO j=iymin, iymax
        DO i=ixmin, ixmax
          ic = iixp0+i
          jc = ijxp0+j
          kc = ikxp0+k
          IF(i<ixmax) THEN
            sdxi  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) +              &
            (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
            IF (i>ixmin) sdxi = sdxi + sdxim1
            !$acc atomic update
            jx(ic, jc, kc) = jx(ic, jc, kc) + sdxi
          END IF
          IF(j<iymax) THEN
            sdyj(i)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) +              &
            (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
            IF (j>iymin) sdyj(i) = sdyj(i) + sdyjm1(i)
            !$acc atomic update
            jy(ic, jc, kc) = jy(ic, jc, kc) + sdyj(i)
          END IF
          IF(k<izmax) THEN
            sdzk(i,j)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) +              &
            (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
            IF (k>izmin) sdzk(i,j) = sdzk(i,j) + sdzkm1(i,j)
            !$acc atomic update
            jz(ic, jc, kc) = jz(ic, jc, kc) + sdzk(i,j)
          END IF
          sdxim1 = sdxi
        END DO
        sdyjm1 = sdyj
      END DO
      sdzkm1 = sdzk
    END DO
  END DO
!$acc end loop
!$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_2_2_2

#if defined (DEV)
! ______________________________________________________________________________
!> @brief
!> Current deposition with the Esirkepov method
!
!> @details
!> Implementation is based on Vincenti's method used for the classical current deposition.
!> Despite vectorization, this subroutine does not exhibit gain in performance
!> with AVX512 architecture.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2015
!
!> @warning
!> PROBLEM: CORRECTION REQUIRED
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_vecHV_2_2_2(jx, jy, jz, np, xp, yp, zp, uxp, uyp,  &
  uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard,    &
  nzguard, nox, noy, noz, l_particles_weight, l4symtry)
  USE constants, ONLY: clight, lvec
  USE picsar_precision, ONLY: idp, isp, lp, num
  IMPLICIT NONE
  INTEGER(idp)                             :: np, nx, ny, nz, nox, noy, noz, nxguard, &
  nyguard, nzguard
  REAL(num), INTENT(IN OUT)         ::                                                &
  jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT)         ::                                                &
  jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT)         ::                                                &
  jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), DIMENSION(np)                 :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                                :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num)                                :: dxi, dyi, dzi, dtsdx, dtsdy, dtsdz,     &
  xint, yint, zint
  REAL(num)                                :: clghtisq, usq, xold, yold, zold
  REAL(num)                                :: xmid, ymid, zmid
  REAL(num)                                :: x, y, z, wq, wqx, wqy, wqz, tmp, vx,    &
  vy, vz
  REAL(num)                                :: invvol, invdtdx, invdtdy, invdtdz,      &
  oxint, oyint, ozint, xintsq, yintsq, zintsq, oxintsq, oyintsq, ozintsq, dtsdx0,     &
  dtsdy0, dtsdz0
  REAL(num), DIMENSION(:, :), ALLOCATABLE   :: jxcells, jycells, jzcells
  REAL(num), PARAMETER                     :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                     :: onethird=1.0_num/3.0_num
  REAL(num), PARAMETER                     :: twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(:), ALLOCATABLE     :: sx, sx0, dsx
  REAL(num), DIMENSION(:), ALLOCATABLE     :: sy, sy0, dsy
  REAL(num), DIMENSION(:), ALLOCATABLE     :: sz, sz0, dsz
  REAL(num), DIMENSION(LVEC, 120)           :: sdx, sdy, sdz
  INTEGER(isp), DIMENSION(LVEC, 3)          :: ICELL
  INTEGER(isp), DIMENSION(14)              :: moffjxc, moffjyc, moffjzc
  INTEGER(isp)                             :: iixporig, ijxporig, ikxporig
  INTEGER(isp)                             :: NCELLS
  INTEGER(isp)                             :: iixp0, ijxp0, ikxp0
  INTEGER(isp)                             :: iixp, ijxp, ikxp
  INTEGER(isp)                             :: ip, dix, diy, diz, idx, idy, idz, i, j, &
  k, ic, jc, kc, ixmin, ixmax, iymin, iymax, izmin, izmax
  INTEGER(isp)                             :: ncx, ncxy, ncy, ncz
  INTEGER(isp)                             :: ngridx, ngridy, nnx, nnxy
  INTEGER(isp)                             :: n, nn, nv
  INTEGER(isp)                             :: igrid, ix, iy, iz, orig
  INTEGER(isp)                             :: moffjx(1:8), moffjy(1:8), moffjz(1:8)
  LOGICAL(lp)                              :: l_particles_weight, l4symtry

  ! __________________________________________________________
  ! Computation of the parameters

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dy*dz)
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)
  clghtisq = 1.0_num/clight**2

  ALLOCATE(sx(-2:2), sx0(-2:2), dsx(-2:2))
  ALLOCATE(sy(-2:2), sy0(-2:2), dsy(-2:2))
  ALLOCATE(sz(-2:2), sz0(-2:2), dsz(-2:2))
  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdy=0.0_num;sdz=0.0_num

  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard
  ncx=nx+1+2*nxguard
  ncy=ny+1+2*nyguard
  ncz=nz+1+2*nzguard
  ncxy=ncx*ncy
  NCELLS=ncx*ncy*ncz
  nnx = ngridx
  nnxy = ngridx*ngridy
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))

  iixporig=-nxguard
  ijxporig=-nyguard
  ikxporig=-nzguard

  moffjx = (/0_isp, 1_isp, 2_isp, 3_isp, nnx, 1_isp+nnx, 2_isp+nnx, 3_isp+nnx/)
  moffjy = (/0_isp, nnx, 2_isp*nnx, 3_isp*nnx, 1_isp, 1_isp+nnx, 1_isp+2_isp*nnx,     &
  1_isp+3_isp*nnx/)
  moffjz = (/0_isp, nnxy, 2_isp*nnxy, 3_isp*nnxy, 1_isp, 1_isp+nnxy,                  &
  1_isp+2_isp*nnxy, 1_isp+3_isp*nnxy/)

  moffjxc = (/2_isp*ncx, 4_isp*ncx, ncxy, ncxy+2_isp*ncx, ncxy+4_isp*ncx, 2_isp*ncxy, &
  2_isp*ncxy+2_isp*ncx, 2_isp*ncxy+4_isp*ncx, 3_isp*ncxy, 3_isp*ncxy+2_isp*ncx,       &
  3_isp*ncxy+4_isp*ncx, 4_isp*ncxy, 4_isp*ncxy+2_isp*ncx, 4_isp*ncxy+4_isp*ncx/)

  orig=(nxguard+iixporig) + (nyguard+ijxporig)*nnx + (nzguard+ikxporig)*nnxy

  ! ______________________________________________________
  ! Loop ober the particles
  DO ip=1, np, LVEC
#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, ICELL)
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
    DO n=1, MIN(LVEC, np-ip+1)
      nn=ip+n-1
      ! --- computes current position in grid units
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi
      ! --- computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)
      ! --- computes old position in grid units
      xold=x-dtsdx0*vx
      yold=y-dtsdy0*vy
      zold=z-dtsdz0*vz
      ! --- computes particles weights
      wq=q*w(ip)
      wqx = wq*invdtdx
      wqy = wq*invdtdy
      wqz = wq*invdtdz
      ! --- finds node of cell containing particles for current positions
      iixp0=nint(x)
      ijxp0=nint(y)
      ikxp0=nint(z)

      ! Cell position with shift (-2, -2, -2)
      ICELL(n, 1)=(iixp0-iixporig-1)+(ijxp0-ijxporig-2)*ncx+(ikxp0-ikxporig-2)*ncxy

      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      yint=y-ijxp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
      xintsq = xint*xint
      sx0(-1) = 0.5_num*(0.5_num-xint)**2
      sx0( 0) = 0.75_num-xintsq
      sx0( 1) = 0.5_num*(0.5_num+xint)**2
      yintsq = yint*yint
      sy0(-1) = 0.5_num*(0.5_num-yint)**2
      sy0( 0) = 0.75_num-yintsq
      sy0( 1) = 0.5_num*(0.5_num+yint)**2
      zintsq = zint*zint
      sz0(-1) = 0.5_num*(0.5_num-zint)**2
      sz0( 0) = 0.75_num-zintsq
      sz0( 1) = 0.5_num*(0.5_num+zint)**2
      ! --- finds node of cell containing particles for old positions
      iixp=nint(xold)
      ijxp=nint(yold)
      ikxp=nint(zold)
      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      yint = yold-ijxp
      zint = zold-ikxp
      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diy = ijxp-ijxp0
      diz = ikxp-ikxp0
      ! --- zero out coefficients
      ! --- (needed because of different dix and diz for each particle)
      sx(-2)=0.0_num;sy(-2)=0.0_num;sz(-2)=0.0_num
      sx(-1)=0.0_num;sy(-1)=0.0_num;sz(-1)=0.0_num
      sx(0)=0.0_num;sy(0)=0.0_num;sz(0)=0.0_num
      sx(1)=0.0_num;sy(1)=0.0_num;sz(1)=0.0_num
      sx(2)=0.0_num;sy(2)=0.0_num;sz(2)=0.0_num
      ! --- computes coefficients for quantities centered between nodes
      xintsq = xint*xint
      sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
      sx( 0+dix) = 0.75_num-xintsq
      sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
      yintsq = yint*yint
      sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
      sy( 0+diy) = 0.75_num-yintsq
      sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
      zintsq = zint*zint
      sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
      sz( 0+diz) = 0.75_num-zintsq
      sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
      ! --- computes coefficients difference
      dsx = sx - sx0
      dsy = sy - sy0
      dsz = sz - sz0

      ! --- Debugging
      IF (ICELL(n, 1) .gt. NCELLS) THEN
        print*, 'Particle', n, nn
        print*, 'ICELL', ICELL(n, 1), NCELLS, ncx, ncy, ncz, ncxy
        print*, 'iixporig', iixporig, ijxporig, ikxporig
        print*, 'iixp0', iixp0, ijxp0, ikxp0
        stop
      ENDIF

      ! --- Compute weights for x
      sdx(n, 1)  = wqx*dsx(-2)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-2) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-2))
      sdx(n, 2)  = wqx*dsx(-1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-2) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-2))
      sdx(n, 2)=sdx(n, 2)+sdx(n, 1)
      sdx(n, 3)  = wqx*dsx(0)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-2) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-2))
      sdx(n, 3)=sdx(n, 3)+sdx(n, 2)
      sdx(n, 4)  = wqx*dsx(1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-2) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-2))
      sdx(n, 4)=sdx(n, 4)+sdx(n, 3)
      sdx(n, 5)  = wqx*dsx(-2)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-2) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-2))
      sdx(n, 6)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-2) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-2))
      sdx(n, 6)=sdx(n, 6)+sdx(n, 5)
      sdx(n, 7)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-2) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-2))
      sdx(n, 7)=sdx(n, 7)+sdx(n, 6)
      sdx(n, 8)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-2) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-2))
      sdx(n, 8)=sdx(n, 8)+sdx(n, 7)
      sdx(n, 9)  = wqx*dsx(-2)*((sy0(0)+0.5_num*dsy(0))*sz0(-2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-2))
      sdx(n, 10)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(-2) +                    &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-2))
      sdx(n, 10)=sdx(n, 10)+sdx(n, 9)
      sdx(n, 11)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(-2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-2))
      sdx(n, 11)=sdx(n, 11)+sdx(n, 10)
      sdx(n, 12)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(-2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-2))
      sdx(n, 12)=sdx(n, 12)+sdx(n, 11)
      sdx(n, 13)  = wqx*dsx(-2)*((sy0(1)+0.5_num*dsy(1))*sz0(-2) +                    &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-2))
      sdx(n, 14)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(-2) +                    &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-2))
      sdx(n, 14)=sdx(n, 14)+sdx(n, 13)
      sdx(n, 15)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(-2) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-2))
      sdx(n, 15)=sdx(n, 15)+sdx(n, 14)
      sdx(n, 16)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(-2) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-2))
      sdx(n, 16)=sdx(n, 16)+sdx(n, 15)
      sdx(n, 17)  = wqx*dsx(-2)*((sy0(2)+0.5_num*dsy(2))*sz0(-2) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-2))
      sdx(n, 18)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(-2) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-2))
      sdx(n, 18)=sdx(n, 18)+sdx(n, 17)
      sdx(n, 19)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(-2) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-2))
      sdx(n, 19)=sdx(n, 19)+sdx(n, 18)
      sdx(n, 20)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(-2) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-2))
      sdx(n, 20)=sdx(n, 20)+sdx(n, 19)
      sdx(n, 21)  = 0.
      sdx(n, 22)  = 0.
      sdx(n, 23)  = 0.
      sdx(n, 24)  = 0.
      sdx(n, 25)  = wqx*dsx(-2)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-1) +                  &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-1))
      sdx(n, 26)  = wqx*dsx(-1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-1) +                  &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-1))
      sdx(n, 26)=sdx(n, 26)+sdx(n, 25)
      sdx(n, 27)  = wqx*dsx(0)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-1) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-1))
      sdx(n, 27)=sdx(n, 27)+sdx(n, 26)
      sdx(n, 28)  = wqx*dsx(1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(-1) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(-1))
      sdx(n, 28)=sdx(n, 28)+sdx(n, 27)
      sdx(n, 29)  = wqx*dsx(-2)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                  &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 30)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                  &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 30)=sdx(n, 30)+sdx(n, 29)
      sdx(n, 31)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 31)=sdx(n, 31)+sdx(n, 30)
      sdx(n, 32)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(-1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(-1))
      sdx(n, 32)=sdx(n, 32)+sdx(n, 31)
      sdx(n, 33)  = wqx*dsx(-2)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                    &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 34)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                    &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 34)=sdx(n, 34)+sdx(n, 33)
      sdx(n, 35)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 35)=sdx(n, 35)+sdx(n, 34)
      sdx(n, 36)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(-1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(-1))
      sdx(n, 36)=sdx(n, 36)+sdx(n, 35)
      sdx(n, 37)  = wqx*dsx(-2)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                    &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 38)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                    &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 38)=sdx(n, 38)+sdx(n, 37)
      sdx(n, 39)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 39)=sdx(n, 39)+sdx(n, 38)
      sdx(n, 40)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(-1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(-1))
      sdx(n, 40)=sdx(n, 40)+sdx(n, 39)
      sdx(n, 41)  = wqx*dsx(-2)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 42)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 42)=sdx(n, 42)+sdx(n, 41)
      sdx(n, 43)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 43)=sdx(n, 43)+sdx(n, 42)
      sdx(n, 44)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(-1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(-1))
      sdx(n, 44)=sdx(n, 44)+sdx(n, 43)
      sdx(n, 45)  = 0.
      sdx(n, 46)  = 0.
      sdx(n, 47)  = 0.
      sdx(n, 48)  = 0.
      sdx(n, 49)  = wqx*dsx(-2)*((sy0(-2)+0.5_num*dsy(-2))*sz0(0) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(0))
      sdx(n, 50)  = wqx*dsx(-1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(0) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(0))
      sdx(n, 50)=sdx(n, 50)+sdx(n, 49)
      sdx(n, 51)  = wqx*dsx(0)*((sy0(-2)+0.5_num*dsy(-2))*sz0(0) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(0))
      sdx(n, 51)=sdx(n, 51)+sdx(n, 50)
      sdx(n, 52)  = wqx*dsx(1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(0) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(0))
      sdx(n, 52)=sdx(n, 52)+sdx(n, 51)
      sdx(n, 53)  = wqx*dsx(-2)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 54)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 54)=sdx(n, 54)+sdx(n, 53)
      sdx(n, 55)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 55)=sdx(n, 55)+sdx(n, 54)
      sdx(n, 56)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(0) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(0))
      sdx(n, 56)=sdx(n, 56)+sdx(n, 55)
      sdx(n, 57)  = wqx*dsx(-2)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 58)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 58)=sdx(n, 58)+sdx(n, 57)
      sdx(n, 59)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 59)=sdx(n, 59)+sdx(n, 58)
      sdx(n, 60)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(0) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(0))
      sdx(n, 60)=sdx(n, 60)+sdx(n, 59)
      sdx(n, 61)  = wqx*dsx(-2)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 62)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 62)=sdx(n, 62)+sdx(n, 61)
      sdx(n, 63)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 63)=sdx(n, 63)+sdx(n, 62)
      sdx(n, 64)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(0) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(0))
      sdx(n, 64)=sdx(n, 64)+sdx(n, 63)
      sdx(n, 65)  = wqx*dsx(-2)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 66)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 66)=sdx(n, 66)+sdx(n, 65)
      sdx(n, 67)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 67)=sdx(n, 67)+sdx(n, 66)
      sdx(n, 68)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(0) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(0))
      sdx(n, 68)=sdx(n, 68)+sdx(n, 67)
      sdx(n, 69)  = 0.
      sdx(n, 70)  = 0.
      sdx(n, 71)  = 0.
      sdx(n, 72)  = 0.
      sdx(n, 73)  = wqx*dsx(-2)*((sy0(-2)+0.5_num*dsy(-2))*sz0(1) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(1))
      sdx(n, 74)  = wqx*dsx(-1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(1) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(1))
      sdx(n, 74)=sdx(n, 74)+sdx(n, 73)
      sdx(n, 75)  = wqx*dsx(0)*((sy0(-2)+0.5_num*dsy(-2))*sz0(1) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(1))
      sdx(n, 75)=sdx(n, 75)+sdx(n, 74)
      sdx(n, 76)  = wqx*dsx(1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(1) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(1))
      sdx(n, 76)=sdx(n, 76)+sdx(n, 75)
      sdx(n, 77)  = wqx*dsx(-2)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 78)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 78)=sdx(n, 78)+sdx(n, 77)
      sdx(n, 79)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 79)=sdx(n, 79)+sdx(n, 78)
      sdx(n, 80)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(1) +                    &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(1))
      sdx(n, 80)=sdx(n, 80)+sdx(n, 79)
      sdx(n, 81)  = wqx*dsx(-2)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 82)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 82)=sdx(n, 82)+sdx(n, 81)
      sdx(n, 83)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 83)=sdx(n, 83)+sdx(n, 82)
      sdx(n, 84)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(1) +                      &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(1))
      sdx(n, 84)=sdx(n, 84)+sdx(n, 83)
      sdx(n, 85)  = wqx*dsx(-2)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 86)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 86)=sdx(n, 86)+sdx(n, 85)
      sdx(n, 87)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 87)=sdx(n, 87)+sdx(n, 86)
      sdx(n, 88)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(1) +                      &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(1))
      sdx(n, 88)=sdx(n, 88)+sdx(n, 87)
      sdx(n, 89)  = wqx*dsx(-2)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 90)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 90)=sdx(n, 90)+sdx(n, 89)
      sdx(n, 91)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 91)=sdx(n, 91)+sdx(n, 90)
      sdx(n, 92)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(1) +                      &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(1))
      sdx(n, 92)=sdx(n, 92)+sdx(n, 91)
      sdx(n, 93)  = 0.
      sdx(n, 94)  = 0.
      sdx(n, 95)  = 0.
      sdx(n, 96)  = 0.
      sdx(n, 97)  = wqx*dsx(-2)*((sy0(-2)+0.5_num*dsy(-2))*sz0(2) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(2))
      sdx(n, 98)  = wqx*dsx(-1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(2) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(2))
      sdx(n, 98)=sdx(n, 98)+sdx(n, 97)
      sdx(n, 99)  = wqx*dsx(0)*((sy0(-2)+0.5_num*dsy(-2))*sz0(2) +                    &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(2))
      sdx(n, 99)=sdx(n, 99)+sdx(n, 98)
      sdx(n, 100)  = wqx*dsx(1)*((sy0(-2)+0.5_num*dsy(-2))*sz0(2) +                   &
      (0.5_num*sy0(-2)+onethird*dsy(-2))*dsz(2))
      sdx(n, 100)=sdx(n, 100)+sdx(n, 99)
      sdx(n, 101)  = wqx*dsx(-2)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                  &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 102)  = wqx*dsx(-1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                  &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 102)=sdx(n, 102)+sdx(n, 101)
      sdx(n, 103)  = wqx*dsx(0)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 103)=sdx(n, 103)+sdx(n, 102)
      sdx(n, 104)  = wqx*dsx(1)*((sy0(-1)+0.5_num*dsy(-1))*sz0(2) +                   &
      (0.5_num*sy0(-1)+onethird*dsy(-1))*dsz(2))
      sdx(n, 104)=sdx(n, 104)+sdx(n, 103)
      sdx(n, 105)  = wqx*dsx(-2)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                    &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 106)  = wqx*dsx(-1)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                    &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 106)=sdx(n, 106)+sdx(n, 105)
      sdx(n, 107)  = wqx*dsx(0)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 107)=sdx(n, 107)+sdx(n, 106)
      sdx(n, 108)  = wqx*dsx(1)*((sy0(0)+0.5_num*dsy(0))*sz0(2) +                     &
      (0.5_num*sy0(0)+onethird*dsy(0))*dsz(2))
      sdx(n, 108)=sdx(n, 108)+sdx(n, 107)
      sdx(n, 109)  = wqx*dsx(-2)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                    &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 110)  = wqx*dsx(-1)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                    &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 110)=sdx(n, 110)+sdx(n, 109)
      sdx(n, 111)  = wqx*dsx(0)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 111)=sdx(n, 111)+sdx(n, 110)
      sdx(n, 112)  = wqx*dsx(1)*((sy0(1)+0.5_num*dsy(1))*sz0(2) +                     &
      (0.5_num*sy0(1)+onethird*dsy(1))*dsz(2))
      sdx(n, 112)=sdx(n, 112)+sdx(n, 111)
      sdx(n, 113)  = wqx*dsx(-2)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 114)  = wqx*dsx(-1)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                    &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 114)=sdx(n, 114)+sdx(n, 113)
      sdx(n, 115)  = wqx*dsx(0)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 115)=sdx(n, 115)+sdx(n, 114)
      sdx(n, 116)  = wqx*dsx(1)*((sy0(2)+0.5_num*dsy(2))*sz0(2) +                     &
      (0.5_num*sy0(2)+onethird*dsy(2))*dsz(2))
      sdx(n, 116)=sdx(n, 116)+sdx(n, 115)
      sdx(n, 117)  = 0.
      sdx(n, 118)  = 0.
      sdx(n, 119)  = 0.
      sdx(n, 120)  = 0.

      ! --- Compute weights for y
      sdy(n, 1)  = wqy*dsy(-2)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-2) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-2))
      sdy(n, 2)  = wqy*dsy(-1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-2) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-2))
      sdy(n, 2)=sdy(n, 2)+sdy(n, 1)
      sdy(n, 3)  = wqy*dsy(0)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-2) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-2))
      sdy(n, 3)=sdy(n, 3)+sdy(n, 2)
      sdy(n, 4)  = wqy*dsy(1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-2) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-2))
      sdy(n, 4)=sdy(n, 4)+sdy(n, 3)
      sdy(n, 5)  = wqy*dsy(-2)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-1) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-1))
      sdy(n, 6)  = wqy*dsy(-1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-1) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-1))
      sdy(n, 6)=sdy(n, 6)+sdy(n, 5)
      sdy(n, 7)  = wqy*dsy(0)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-1) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-1))
      sdy(n, 7)=sdy(n, 7)+sdy(n, 6)
      sdy(n, 8)  = wqy*dsy(1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(-1) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(-1))
      sdy(n, 8)=sdy(n, 8)+sdy(n, 7)
      sdy(n, 9)  = wqy*dsy(-2)*((sz0(-2)+0.5_num*dsz(-2))*sx0(0) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(0))
      sdy(n, 10)  = wqy*dsy(-1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(0) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(0))
      sdy(n, 10)=sdy(n, 10)+sdy(n, 9)
      sdy(n, 11)  = wqy*dsy(0)*((sz0(-2)+0.5_num*dsz(-2))*sx0(0) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(0))
      sdy(n, 11)=sdy(n, 11)+sdy(n, 10)
      sdy(n, 12)  = wqy*dsy(1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(0) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(0))
      sdy(n, 12)=sdy(n, 12)+sdy(n, 11)
      sdy(n, 13)  = wqy*dsy(-2)*((sz0(-2)+0.5_num*dsz(-2))*sx0(1) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(1))
      sdy(n, 14)  = wqy*dsy(-1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(1) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(1))
      sdy(n, 14)=sdy(n, 14)+sdy(n, 13)
      sdy(n, 15)  = wqy*dsy(0)*((sz0(-2)+0.5_num*dsz(-2))*sx0(1) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(1))
      sdy(n, 15)=sdy(n, 15)+sdy(n, 14)
      sdy(n, 16)  = wqy*dsy(1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(1) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(1))
      sdy(n, 16)=sdy(n, 16)+sdy(n, 15)
      sdy(n, 17)  = wqy*dsy(-2)*((sz0(-2)+0.5_num*dsz(-2))*sx0(2) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(2))
      sdy(n, 18)  = wqy*dsy(-1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(2) +                   &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(2))
      sdy(n, 18)=sdy(n, 18)+sdy(n, 17)
      sdy(n, 19)  = wqy*dsy(0)*((sz0(-2)+0.5_num*dsz(-2))*sx0(2) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(2))
      sdy(n, 19)=sdy(n, 19)+sdy(n, 18)
      sdy(n, 20)  = wqy*dsy(1)*((sz0(-2)+0.5_num*dsz(-2))*sx0(2) +                    &
      (0.5_num*sz0(-2)+onethird*dsz(-2))*dsx(2))
      sdy(n, 20)=sdy(n, 20)+sdy(n, 19)
      sdy(n, 21)  = 0.
      sdy(n, 22)  = 0.
      sdy(n, 23)  = 0.
      sdy(n, 24)  = 0.
      sdy(n, 25)  = wqy*dsy(-2)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-2) +                  &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-2))
      sdy(n, 26)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-2) +                  &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-2))
      sdy(n, 26)=sdy(n, 26)+sdy(n, 25)
      sdy(n, 27)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-2) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-2))
      sdy(n, 27)=sdy(n, 27)+sdy(n, 26)
      sdy(n, 28)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-2) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-2))
      sdy(n, 28)=sdy(n, 28)+sdy(n, 27)
      sdy(n, 29)  = wqy*dsy(-2)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                  &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 30)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                  &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 30)=sdy(n, 30)+sdy(n, 29)
      sdy(n, 31)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 31)=sdy(n, 31)+sdy(n, 30)
      sdy(n, 32)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(-1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(-1))
      sdy(n, 32)=sdy(n, 32)+sdy(n, 31)
      sdy(n, 33)  = wqy*dsy(-2)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 34)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 34)=sdy(n, 34)+sdy(n, 33)
      sdy(n, 35)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 35)=sdy(n, 35)+sdy(n, 34)
      sdy(n, 36)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(0) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(0))
      sdy(n, 36)=sdy(n, 36)+sdy(n, 35)
      sdy(n, 37)  = wqy*dsy(-2)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 38)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 38)=sdy(n, 38)+sdy(n, 37)
      sdy(n, 39)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 39)=sdy(n, 39)+sdy(n, 38)
      sdy(n, 40)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(1) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(1))
      sdy(n, 40)=sdy(n, 40)+sdy(n, 39)
      sdy(n, 41)  = wqy*dsy(-2)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 42)  = wqy*dsy(-1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                   &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 42)=sdy(n, 42)+sdy(n, 41)
      sdy(n, 43)  = wqy*dsy(0)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 43)=sdy(n, 43)+sdy(n, 42)
      sdy(n, 44)  = wqy*dsy(1)*((sz0(-1)+0.5_num*dsz(-1))*sx0(2) +                    &
      (0.5_num*sz0(-1)+onethird*dsz(-1))*dsx(2))
      sdy(n, 44)=sdy(n, 44)+sdy(n, 43)
      sdy(n, 45)  = 0.
      sdy(n, 46)  = 0.
      sdy(n, 47)  = 0.
      sdy(n, 48)  = 0.
      sdy(n, 49)  = wqy*dsy(-2)*((sz0(0)+0.5_num*dsz(0))*sx0(-2) +                    &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-2))
      sdy(n, 50)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(-2) +                    &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-2))
      sdy(n, 50)=sdy(n, 50)+sdy(n, 49)
      sdy(n, 51)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(-2) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-2))
      sdy(n, 51)=sdy(n, 51)+sdy(n, 50)
      sdy(n, 52)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(-2) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-2))
      sdy(n, 52)=sdy(n, 52)+sdy(n, 51)
      sdy(n, 53)  = wqy*dsy(-2)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                    &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 54)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                    &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 54)=sdy(n, 54)+sdy(n, 53)
      sdy(n, 55)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 55)=sdy(n, 55)+sdy(n, 54)
      sdy(n, 56)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(-1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(-1))
      sdy(n, 56)=sdy(n, 56)+sdy(n, 55)
      sdy(n, 57)  = wqy*dsy(-2)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 58)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 58)=sdy(n, 58)+sdy(n, 57)
      sdy(n, 59)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 59)=sdy(n, 59)+sdy(n, 58)
      sdy(n, 60)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(0) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(0))
      sdy(n, 60)=sdy(n, 60)+sdy(n, 59)
      sdy(n, 61)  = wqy*dsy(-2)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 62)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 62)=sdy(n, 62)+sdy(n, 61)
      sdy(n, 63)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 63)=sdy(n, 63)+sdy(n, 62)
      sdy(n, 64)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(1) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(1))
      sdy(n, 64)=sdy(n, 64)+sdy(n, 63)
      sdy(n, 65)  = wqy*dsy(-2)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 66)  = wqy*dsy(-1)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                     &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 66)=sdy(n, 66)+sdy(n, 65)
      sdy(n, 67)  = wqy*dsy(0)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 67)=sdy(n, 67)+sdy(n, 66)
      sdy(n, 68)  = wqy*dsy(1)*((sz0(0)+0.5_num*dsz(0))*sx0(2) +                      &
      (0.5_num*sz0(0)+onethird*dsz(0))*dsx(2))
      sdy(n, 68)=sdy(n, 68)+sdy(n, 67)
      sdy(n, 69)  = 0.
      sdy(n, 70)  = 0.
      sdy(n, 71)  = 0.
      sdy(n, 72)  = 0.
      sdy(n, 73)  = wqy*dsy(-2)*((sz0(1)+0.5_num*dsz(1))*sx0(-2) +                    &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-2))
      sdy(n, 74)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(-2) +                    &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-2))
      sdy(n, 74)=sdy(n, 74)+sdy(n, 73)
      sdy(n, 75)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(-2) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-2))
      sdy(n, 75)=sdy(n, 75)+sdy(n, 74)
      sdy(n, 76)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(-2) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-2))
      sdy(n, 76)=sdy(n, 76)+sdy(n, 75)
      sdy(n, 77)  = wqy*dsy(-2)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                    &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 78)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                    &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 78)=sdy(n, 78)+sdy(n, 77)
      sdy(n, 79)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 79)=sdy(n, 79)+sdy(n, 78)
      sdy(n, 80)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(-1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(-1))
      sdy(n, 80)=sdy(n, 80)+sdy(n, 79)
      sdy(n, 81)  = wqy*dsy(-2)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 82)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 82)=sdy(n, 82)+sdy(n, 81)
      sdy(n, 83)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 83)=sdy(n, 83)+sdy(n, 82)
      sdy(n, 84)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(0) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(0))
      sdy(n, 84)=sdy(n, 84)+sdy(n, 83)
      sdy(n, 85)  = wqy*dsy(-2)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 86)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 86)=sdy(n, 86)+sdy(n, 85)
      sdy(n, 87)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 87)=sdy(n, 87)+sdy(n, 86)
      sdy(n, 88)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(1) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(1))
      sdy(n, 88)=sdy(n, 88)+sdy(n, 87)
      sdy(n, 89)  = wqy*dsy(-2)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 90)  = wqy*dsy(-1)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                     &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 90)=sdy(n, 90)+sdy(n, 89)
      sdy(n, 91)  = wqy*dsy(0)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 91)=sdy(n, 91)+sdy(n, 90)
      sdy(n, 92)  = wqy*dsy(1)*((sz0(1)+0.5_num*dsz(1))*sx0(2) +                      &
      (0.5_num*sz0(1)+onethird*dsz(1))*dsx(2))
      sdy(n, 92)=sdy(n, 92)+sdy(n, 91)
      sdy(n, 93)  = 0.
      sdy(n, 94)  = 0.
      sdy(n, 95)  = 0.
      sdy(n, 96)  = 0.
      sdy(n, 97)  = wqy*dsy(-2)*((sz0(2)+0.5_num*dsz(2))*sx0(-2) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-2))
      sdy(n, 98)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(-2) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-2))
      sdy(n, 98)=sdy(n, 98)+sdy(n, 97)
      sdy(n, 99)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(-2) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-2))
      sdy(n, 99)=sdy(n, 99)+sdy(n, 98)
      sdy(n, 100)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(-2) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-2))
      sdy(n, 100)=sdy(n, 100)+sdy(n, 99)
      sdy(n, 101)  = wqy*dsy(-2)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                   &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 102)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                   &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 102)=sdy(n, 102)+sdy(n, 101)
      sdy(n, 103)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 103)=sdy(n, 103)+sdy(n, 102)
      sdy(n, 104)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(-1) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(-1))
      sdy(n, 104)=sdy(n, 104)+sdy(n, 103)
      sdy(n, 105)  = wqy*dsy(-2)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 106)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 106)=sdy(n, 106)+sdy(n, 105)
      sdy(n, 107)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 107)=sdy(n, 107)+sdy(n, 106)
      sdy(n, 108)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(0) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(0))
      sdy(n, 108)=sdy(n, 108)+sdy(n, 107)
      sdy(n, 109)  = wqy*dsy(-2)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 110)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 110)=sdy(n, 110)+sdy(n, 109)
      sdy(n, 111)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 111)=sdy(n, 111)+sdy(n, 110)
      sdy(n, 112)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(1) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(1))
      sdy(n, 112)=sdy(n, 112)+sdy(n, 111)
      sdy(n, 113)  = wqy*dsy(-2)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 114)  = wqy*dsy(-1)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                    &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 114)=sdy(n, 114)+sdy(n, 113)
      sdy(n, 115)  = wqy*dsy(0)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 115)=sdy(n, 115)+sdy(n, 114)
      sdy(n, 116)  = wqy*dsy(1)*((sz0(2)+0.5_num*dsz(2))*sx0(2) +                     &
      (0.5_num*sz0(2)+onethird*dsz(2))*dsx(2))
      sdy(n, 116)=sdy(n, 116)+sdy(n, 115)
      sdy(n, 117)  = 0.
      sdy(n, 118)  = 0.
      sdy(n, 119)  = 0.
      sdy(n, 120)  = 0.

      ! --- Compute weights for z
      sdz(n, 1)  = wqz*dsz(-2)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-2) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-2))
      sdz(n, 2)  = wqz*dsz(-1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-2) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-2))
      sdz(n, 2)=sdz(n, 2)+sdz(n, 1)
      sdz(n, 3)  = wqz*dsz(0)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-2) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-2))
      sdz(n, 3)=sdz(n, 3)+sdz(n, 2)
      sdz(n, 4)  = wqz*dsz(1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-2) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-2))
      sdz(n, 4)=sdz(n, 4)+sdz(n, 3)
      sdz(n, 5)  = wqz*dsz(-2)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-2) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-2))
      sdz(n, 6)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-2) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-2))
      sdz(n, 6)=sdz(n, 6)+sdz(n, 5)
      sdz(n, 7)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-2) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-2))
      sdz(n, 7)=sdz(n, 7)+sdz(n, 6)
      sdz(n, 8)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-2) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-2))
      sdz(n, 8)=sdz(n, 8)+sdz(n, 7)
      sdz(n, 9)  = wqz*dsz(-2)*((sx0(0)+0.5_num*dsx(0))*sy0(-2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-2))
      sdz(n, 10)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(-2) +                    &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-2))
      sdz(n, 10)=sdz(n, 10)+sdz(n, 9)
      sdz(n, 11)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(-2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-2))
      sdz(n, 11)=sdz(n, 11)+sdz(n, 10)
      sdz(n, 12)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(-2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-2))
      sdz(n, 12)=sdz(n, 12)+sdz(n, 11)
      sdz(n, 13)  = wqz*dsz(-2)*((sx0(1)+0.5_num*dsx(1))*sy0(-2) +                    &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-2))
      sdz(n, 14)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(-2) +                    &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-2))
      sdz(n, 14)=sdz(n, 14)+sdz(n, 13)
      sdz(n, 15)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(-2) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-2))
      sdz(n, 15)=sdz(n, 15)+sdz(n, 14)
      sdz(n, 16)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(-2) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-2))
      sdz(n, 16)=sdz(n, 16)+sdz(n, 15)
      sdz(n, 17)  = wqz*dsz(-2)*((sx0(2)+0.5_num*dsx(2))*sy0(-2) +                    &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-2))
      sdz(n, 18)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(-2) +                    &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-2))
      sdz(n, 18)=sdz(n, 18)+sdz(n, 17)
      sdz(n, 19)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(-2) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-2))
      sdz(n, 19)=sdz(n, 19)+sdz(n, 18)
      sdz(n, 20)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(-2) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-2))
      sdz(n, 20)=sdz(n, 20)+sdz(n, 19)
      sdz(n, 21)  = 0.
      sdz(n, 22)  = 0.
      sdz(n, 23)  = 0.
      sdz(n, 24)  = 0.
      sdz(n, 25)  = wqz*dsz(-2)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-1) +                  &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-1))
      sdz(n, 26)  = wqz*dsz(-1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-1) +                  &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-1))
      sdz(n, 26)=sdz(n, 26)+sdz(n, 25)
      sdz(n, 27)  = wqz*dsz(0)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-1) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-1))
      sdz(n, 27)=sdz(n, 27)+sdz(n, 26)
      sdz(n, 28)  = wqz*dsz(1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(-1) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(-1))
      sdz(n, 28)=sdz(n, 28)+sdz(n, 27)
      sdz(n, 29)  = wqz*dsz(-2)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                  &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 30)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                  &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 30)=sdz(n, 30)+sdz(n, 29)
      sdz(n, 31)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 31)=sdz(n, 31)+sdz(n, 30)
      sdz(n, 32)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(-1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(-1))
      sdz(n, 32)=sdz(n, 32)+sdz(n, 31)
      sdz(n, 33)  = wqz*dsz(-2)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                    &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 34)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                    &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 34)=sdz(n, 34)+sdz(n, 33)
      sdz(n, 35)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 35)=sdz(n, 35)+sdz(n, 34)
      sdz(n, 36)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(-1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(-1))
      sdz(n, 36)=sdz(n, 36)+sdz(n, 35)
      sdz(n, 37)  = wqz*dsz(-2)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                    &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 38)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                    &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 38)=sdz(n, 38)+sdz(n, 37)
      sdz(n, 39)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 39)=sdz(n, 39)+sdz(n, 38)
      sdz(n, 40)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(-1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(-1))
      sdz(n, 40)=sdz(n, 40)+sdz(n, 39)
      sdz(n, 41)  = wqz*dsz(-2)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                    &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 42)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                    &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 42)=sdz(n, 42)+sdz(n, 41)
      sdz(n, 43)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 43)=sdz(n, 43)+sdz(n, 42)
      sdz(n, 44)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(-1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(-1))
      sdz(n, 44)=sdz(n, 44)+sdz(n, 43)
      sdz(n, 45)  = 0.
      sdz(n, 46)  = 0.
      sdz(n, 47)  = 0.
      sdz(n, 48)  = 0.
      sdz(n, 49)  = wqz*dsz(-2)*((sx0(-2)+0.5_num*dsx(-2))*sy0(0) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(0))
      sdz(n, 50)  = wqz*dsz(-1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(0) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(0))
      sdz(n, 50)=sdz(n, 50)+sdz(n, 49)
      sdz(n, 51)  = wqz*dsz(0)*((sx0(-2)+0.5_num*dsx(-2))*sy0(0) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(0))
      sdz(n, 51)=sdz(n, 51)+sdz(n, 50)
      sdz(n, 52)  = wqz*dsz(1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(0) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(0))
      sdz(n, 52)=sdz(n, 52)+sdz(n, 51)
      sdz(n, 53)  = wqz*dsz(-2)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 54)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 54)=sdz(n, 54)+sdz(n, 53)
      sdz(n, 55)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 55)=sdz(n, 55)+sdz(n, 54)
      sdz(n, 56)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(0) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(0))
      sdz(n, 56)=sdz(n, 56)+sdz(n, 55)
      sdz(n, 57)  = wqz*dsz(-2)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 58)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 58)=sdz(n, 58)+sdz(n, 57)
      sdz(n, 59)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 59)=sdz(n, 59)+sdz(n, 58)
      sdz(n, 60)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(0) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(0))
      sdz(n, 60)=sdz(n, 60)+sdz(n, 59)
      sdz(n, 61)  = wqz*dsz(-2)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 62)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 62)=sdz(n, 62)+sdz(n, 61)
      sdz(n, 63)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 63)=sdz(n, 63)+sdz(n, 62)
      sdz(n, 64)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(0) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(0))
      sdz(n, 64)=sdz(n, 64)+sdz(n, 63)
      sdz(n, 65)  = wqz*dsz(-2)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 66)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 66)=sdz(n, 66)+sdz(n, 65)
      sdz(n, 67)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 67)=sdz(n, 67)+sdz(n, 66)
      sdz(n, 68)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(0) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(0))
      sdz(n, 68)=sdz(n, 68)+sdz(n, 67)
      sdz(n, 69)  = 0.
      sdz(n, 70)  = 0.
      sdz(n, 71)  = 0.
      sdz(n, 72)  = 0.
      sdz(n, 73)  = wqz*dsz(-2)*((sx0(-2)+0.5_num*dsx(-2))*sy0(1) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(1))
      sdz(n, 74)  = wqz*dsz(-1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(1) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(1))
      sdz(n, 74)=sdz(n, 74)+sdz(n, 73)
      sdz(n, 75)  = wqz*dsz(0)*((sx0(-2)+0.5_num*dsx(-2))*sy0(1) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(1))
      sdz(n, 75)=sdz(n, 75)+sdz(n, 74)
      sdz(n, 76)  = wqz*dsz(1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(1) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(1))
      sdz(n, 76)=sdz(n, 76)+sdz(n, 75)
      sdz(n, 77)  = wqz*dsz(-2)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 78)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 78)=sdz(n, 78)+sdz(n, 77)
      sdz(n, 79)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 79)=sdz(n, 79)+sdz(n, 78)
      sdz(n, 80)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(1) +                    &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(1))
      sdz(n, 80)=sdz(n, 80)+sdz(n, 79)
      sdz(n, 81)  = wqz*dsz(-2)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 82)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 82)=sdz(n, 82)+sdz(n, 81)
      sdz(n, 83)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 83)=sdz(n, 83)+sdz(n, 82)
      sdz(n, 84)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(1) +                      &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(1))
      sdz(n, 84)=sdz(n, 84)+sdz(n, 83)
      sdz(n, 85)  = wqz*dsz(-2)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 86)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 86)=sdz(n, 86)+sdz(n, 85)
      sdz(n, 87)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 87)=sdz(n, 87)+sdz(n, 86)
      sdz(n, 88)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(1) +                      &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(1))
      sdz(n, 88)=sdz(n, 88)+sdz(n, 87)
      sdz(n, 89)  = wqz*dsz(-2)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 90)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 90)=sdz(n, 90)+sdz(n, 89)
      sdz(n, 91)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 91)=sdz(n, 91)+sdz(n, 90)
      sdz(n, 92)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(1) +                      &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(1))
      sdz(n, 92)=sdz(n, 92)+sdz(n, 91)
      sdz(n, 93)  = 0.
      sdz(n, 94)  = 0.
      sdz(n, 95)  = 0.
      sdz(n, 96)  = 0.
      sdz(n, 97)  = wqz*dsz(-2)*((sx0(-2)+0.5_num*dsx(-2))*sy0(2) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(2))
      sdz(n, 98)  = wqz*dsz(-1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(2) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(2))
      sdz(n, 98)=sdz(n, 98)+sdz(n, 97)
      sdz(n, 99)  = wqz*dsz(0)*((sx0(-2)+0.5_num*dsx(-2))*sy0(2) +                    &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(2))
      sdz(n, 99)=sdz(n, 99)+sdz(n, 98)
      sdz(n, 100)  = wqz*dsz(1)*((sx0(-2)+0.5_num*dsx(-2))*sy0(2) +                   &
      (0.5_num*sx0(-2)+onethird*dsx(-2))*dsy(2))
      sdz(n, 100)=sdz(n, 100)+sdz(n, 99)
      sdz(n, 101)  = wqz*dsz(-2)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                  &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 102)  = wqz*dsz(-1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                  &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 102)=sdz(n, 102)+sdz(n, 101)
      sdz(n, 103)  = wqz*dsz(0)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 103)=sdz(n, 103)+sdz(n, 102)
      sdz(n, 104)  = wqz*dsz(1)*((sx0(-1)+0.5_num*dsx(-1))*sy0(2) +                   &
      (0.5_num*sx0(-1)+onethird*dsx(-1))*dsy(2))
      sdz(n, 104)=sdz(n, 104)+sdz(n, 103)
      sdz(n, 105)  = wqz*dsz(-2)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                    &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 106)  = wqz*dsz(-1)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                    &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 106)=sdz(n, 106)+sdz(n, 105)
      sdz(n, 107)  = wqz*dsz(0)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 107)=sdz(n, 107)+sdz(n, 106)
      sdz(n, 108)  = wqz*dsz(1)*((sx0(0)+0.5_num*dsx(0))*sy0(2) +                     &
      (0.5_num*sx0(0)+onethird*dsx(0))*dsy(2))
      sdz(n, 108)=sdz(n, 108)+sdz(n, 107)
      sdz(n, 109)  = wqz*dsz(-2)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                    &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 110)  = wqz*dsz(-1)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                    &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 110)=sdz(n, 110)+sdz(n, 109)
      sdz(n, 111)  = wqz*dsz(0)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 111)=sdz(n, 111)+sdz(n, 110)
      sdz(n, 112)  = wqz*dsz(1)*((sx0(1)+0.5_num*dsx(1))*sy0(2) +                     &
      (0.5_num*sx0(1)+onethird*dsx(1))*dsy(2))
      sdz(n, 112)=sdz(n, 112)+sdz(n, 111)
      sdz(n, 113)  = wqz*dsz(-2)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                    &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 114)  = wqz*dsz(-1)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                    &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 114)=sdz(n, 114)+sdz(n, 113)
      sdz(n, 115)  = wqz*dsz(0)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 115)=sdz(n, 115)+sdz(n, 114)
      sdz(n, 116)  = wqz*dsz(1)*((sx0(2)+0.5_num*dsx(2))*sy0(2) +                     &
      (0.5_num*sx0(2)+onethird*dsx(2))*dsy(2))
      sdz(n, 116)=sdz(n, 116)+sdz(n, 115)
      sdz(n, 117)  = 0.
      sdz(n, 118)  = 0.
      sdz(n, 119)  = 0.
      sdz(n, 120)  = 0.
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! Add weights to nearest vertices
    DO n=1, MIN(LVEC, np-ip+1)
#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, jxcells, jycells, jzcells)
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
      DO nv=1, 8

        ! --- JX
        jxcells(nv, ICELL(n, 1)) = jxcells(nv, ICELL(n, 1)) + sdx(n, nv)
        jxcells(nv, ICELL(n, 1)+2*ncx) = jxcells(nv, ICELL(n, 1)+2*ncx) + sdx(n,      &
        8+nv)
        jxcells(nv, ICELL(n, 1)+4*ncx) = jxcells(nv, ICELL(n, 1)+4*ncx) + sdx(n,      &
        16+nv)

        jxcells(nv, ICELL(n, 1)+ncxy) = jxcells(nv, ICELL(n, 1)+ncxy) + sdx(n, 24+nv)
        jxcells(nv, ICELL(n, 1)+ncxy+2*ncx) = jxcells(nv, ICELL(n, 1)+ncxy+2*ncx) +   &
        sdx(n, 32+nv)
        jxcells(nv, ICELL(n, 1)+ncxy+4*ncx) = jxcells(nv, ICELL(n, 1)+ncxy+4*ncx) +   &
        sdx(n, 40+nv)

        jxcells(nv, ICELL(n, 1)+2*ncxy) = jxcells(nv, ICELL(n, 1)+2*ncxy) + sdx(n,    &
        48+nv)
        jxcells(nv, ICELL(n, 1)+2*ncxy+2*ncx) = jxcells(nv, ICELL(n, 1)+2*ncxy+2*ncx) &
        + sdx(n, 56+nv)
        jxcells(nv, ICELL(n, 1)+2*ncxy+4*ncx) = jxcells(nv, ICELL(n, 1)+2*ncxy+4*ncx) &
        + sdx(n, 64+nv)

        jxcells(nv, ICELL(n, 1)+3*ncxy) = jxcells(nv, ICELL(n, 1)+3*ncxy) + sdx(n,    &
        72+nv)
        jxcells(nv, ICELL(n, 1)+3*ncxy+2*ncx) = jxcells(nv, ICELL(n, 1)+3*ncxy+2*ncx) &
        + sdx(n, 80+nv)
        jxcells(nv, ICELL(n, 1)+3*ncxy+4*ncx) = jxcells(nv, ICELL(n, 1)+3*ncxy+4*ncx) &
        + sdx(n, 88+nv)

        jxcells(nv, ICELL(n, 1)+4*ncxy) = jxcells(nv, ICELL(n, 1)+4*ncxy) + sdx(n,    &
        96+nv)
        jxcells(nv, ICELL(n, 1)+4*ncxy+2*ncx) = jxcells(nv, ICELL(n, 1)+4*ncxy+2*ncx) &
        + sdx(n, 104+nv)
        jxcells(nv, ICELL(n, 1)+4*ncxy+4*ncx) = jxcells(nv, ICELL(n, 1)+4*ncxy+4*ncx) &
        + sdx(n, 112+nv)

        ! --- JY
        jycells(nv, ICELL(n, 1)) = jycells(nv, ICELL(n, 1)) + sdy(n, nv)
        jycells(nv, ICELL(n, 1)+2) = jycells(nv, ICELL(n, 1)+2) + sdy(n, nv+8)
        jycells(nv, ICELL(n, 1)+4) = jycells(nv, ICELL(n, 1)+4) + sdy(n, nv+16)

        jycells(nv, ICELL(n, 1)+ncxy) = jycells(nv, ICELL(n, 1)+ncxy) + sdy(n, nv+24)
        jycells(nv, ICELL(n, 1)+ncxy+2) = jycells(nv, ICELL(n, 1)+ncxy+2) + sdy(n,    &
        nv+32)
        jycells(nv, ICELL(n, 1)+ncxy+4) = jycells(nv, ICELL(n, 1)+ncxy+4) + sdy(n,    &
        nv+40)

        jycells(nv, ICELL(n, 1)+2*ncxy) = jycells(nv, ICELL(n, 1)+2*ncxy) + sdy(n,    &
        nv+48)
        jycells(nv, ICELL(n, 1)+2*ncxy+2) = jycells(nv, ICELL(n, 1)+2*ncxy+2) +       &
        sdy(n, nv+56)
        jycells(nv, ICELL(n, 1)+2*ncxy+4) = jycells(nv, ICELL(n, 1)+2*ncxy+4) +       &
        sdy(n, nv+64)

        jycells(nv, ICELL(n, 1)+3*ncxy) = jycells(nv, ICELL(n, 1)+3*ncxy) + sdy(n,    &
        nv+72)
        jycells(nv, ICELL(n, 1)+3*ncxy+2) = jycells(nv, ICELL(n, 1)+3*ncxy+2) +       &
        sdy(n, nv+80)
        jycells(nv, ICELL(n, 1)+3*ncxy+4) = jycells(nv, ICELL(n, 1)+3*ncxy+4) +       &
        sdy(n, nv+88)

        jycells(nv, ICELL(n, 1)+4*ncxy) = jycells(nv, ICELL(n, 1)+4*ncxy) + sdy(n,    &
        nv+96)
        jycells(nv, ICELL(n, 1)+4*ncxy+2) = jycells(nv, ICELL(n, 1)+4*ncxy+2) +       &
        sdy(n, nv+104)
        jycells(nv, ICELL(n, 1)+4*ncxy+4) = jycells(nv, ICELL(n, 1)+4*ncxy+4) +       &
        sdy(n, nv+112)

        ! --- JZ
        jzcells(nv, ICELL(n, 1)) = jzcells(nv, ICELL(n, 1)) + sdz(n, nv)
        jzcells(nv, ICELL(n, 1)+2) = jzcells(nv, ICELL(n, 1)+2) + sdz(n, nv+8)
        jzcells(nv, ICELL(n, 1)+4) = jzcells(nv, ICELL(n, 1)+4) + sdz(n, nv+16)

        jzcells(nv, ICELL(n, 1)+ncx) = jzcells(nv, ICELL(n, 1)+ncx) + sdz(n, nv+24)
        jzcells(nv, ICELL(n, 1)+ncx+2) = jzcells(nv, ICELL(n, 1)+ncx+2) + sdz(n,      &
        nv+32)
        jzcells(nv, ICELL(n, 1)+ncx+4) = jzcells(nv, ICELL(n, 1)+ncx+4) + sdz(n,      &
        nv+40)

        jzcells(nv, ICELL(n, 1)+2*ncx) = jzcells(nv, ICELL(n, 1)+2*ncx) + sdz(n,      &
        nv+48)
        jzcells(nv, ICELL(n, 1)+2*ncx+2) = jzcells(nv, ICELL(n, 1)+2*ncx+2) + sdz(n,  &
        nv+56)
        jzcells(nv, ICELL(n, 1)+2*ncx+4) = jzcells(nv, ICELL(n, 1)+2*ncx+4) + sdz(n,  &
        nv+64)

        jzcells(nv, ICELL(n, 1)+3*ncx) = jzcells(nv, ICELL(n, 1)+3*ncx) + sdz(n,      &
        nv+72)
        jzcells(nv, ICELL(n, 1)+3*ncx+2) = jzcells(nv, ICELL(n, 1)+3*ncx+2) + sdz(n,  &
        nv+80)
        jzcells(nv, ICELL(n, 1)+3*ncx+4) = jzcells(nv, ICELL(n, 1)+3*ncx+4) + sdz(n,  &
        nv+88)

        jzcells(nv, ICELL(n, 1)+4*ncx) = jzcells(nv, ICELL(n, 1)+4*ncx) + sdz(n,      &
        nv+96)
        jzcells(nv, ICELL(n, 1)+4*ncx+2) = jzcells(nv, ICELL(n, 1)+4*ncx+2) + sdz(n,  &
        nv+104)
        jzcells(nv, ICELL(n, 1)+4*ncx+4) = jzcells(nv, ICELL(n, 1)+4*ncx+4) + sdz(n,  &
        nv+112)
      ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    ENDDO
  END DO
  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
  DO iz=1, ncz-3
    DO iy=1, ncy-3
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix=1, ncx-3!! VECTOR (take ncx multiple of vector length)
        ic=ix+(iy-1)*ncx+(iz-1)*ncxy
        igrid =orig+ix+(iy-1)*nnx+(iz-1)*nnxy

        ! jx
        jx(igrid+moffjx(1))=jx(igrid+moffjx(1))+jxcells(1, ic)
        jx(igrid+moffjx(2))=jx(igrid+moffjx(2))+jxcells(2, ic)
        jx(igrid+moffjx(3))=jx(igrid+moffjx(3))+jxcells(3, ic)
        jx(igrid+moffjx(4))=jx(igrid+moffjx(4))+jxcells(4, ic)
        jx(igrid+moffjx(5))=jx(igrid+moffjx(5))+jxcells(5, ic)
        jx(igrid+moffjx(6))=jx(igrid+moffjx(6))+jxcells(6, ic)
        jx(igrid+moffjx(7))=jx(igrid+moffjx(7))+jxcells(7, ic)
        jx(igrid+moffjx(8))=jx(igrid+moffjx(8))+jxcells(8, ic)
        ! jy
        jy(igrid+moffjy(1))=jy(igrid+moffjy(1))+jycells(1, ic)
        jy(igrid+moffjy(2))=jy(igrid+moffjy(2))+jycells(2, ic)
        jy(igrid+moffjy(3))=jy(igrid+moffjy(3))+jycells(3, ic)
        jy(igrid+moffjy(4))=jy(igrid+moffjy(4))+jycells(4, ic)
        jy(igrid+moffjy(5))=jy(igrid+moffjy(5))+jycells(5, ic)
        jy(igrid+moffjy(6))=jy(igrid+moffjy(6))+jycells(6, ic)
        jy(igrid+moffjy(7))=jy(igrid+moffjy(7))+jycells(7, ic)
        jy(igrid+moffjy(8))=jy(igrid+moffjy(8))+jycells(8, ic)
        ! jz
        jz(igrid+moffjz(1))=jz(igrid+moffjz(1))+jzcells(1, ic)
        jz(igrid+moffjz(2))=jz(igrid+moffjz(2))+jzcells(2, ic)
        jz(igrid+moffjz(3))=jz(igrid+moffjz(3))+jzcells(3, ic)
        jz(igrid+moffjz(4))=jz(igrid+moffjz(4))+jzcells(4, ic)
        jz(igrid+moffjz(5))=jz(igrid+moffjz(5))+jzcells(5, ic)
        jz(igrid+moffjz(6))=jz(igrid+moffjz(6))+jzcells(6, ic)
        jz(igrid+moffjz(7))=jz(igrid+moffjz(7))+jzcells(7, ic)
        jz(igrid+moffjz(8))=jz(igrid+moffjz(8))+jzcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  DEALLOCATE(jxcells, jycells, jzcells)
  RETURN

END SUBROUTINE depose_jxjyjz_esirkepov_vecHV_2_2_2
#endif

! ________________________________________________________________________________________
!> @brief
!> Esirkepov scalar current deposition at order 3 in x, y, z (nox=noy=noz=3)
!>
!> @detail
!> This function is not vectorized
!>
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!>
!> @date
!> Revision 10/09/2016
!>
!> @param[inout] jx x-current component (3D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 3 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the jx array
!> (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!> (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the jy array
!> (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the jz array
!> (1d array containing 3 integers)
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv particle Lorentz factor arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] nox, noy, noz interpolation order
!> @param[in] l_particles_weight use the particle weigth
!> @param[in] l4symtry
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz)    !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER :: np
  INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), jy_nguard(3), jy_nvalid(3), &
  jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1,                                          &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1,                                          &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,                                          &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: sdx, sdy, sdz
  REAL(num) :: clghtisq, xold, yold, zold
  REAL(num) :: x, y, z, wq, wqx, wqy, wqz, vx, vy, vz
  REAL(num) :: invvol, invdtdx, invdtdy, invdtdz
  REAL(num) :: oxint, oyint, ozint, xintsq, yintsq, zintsq, oxintsq, oyintsq, ozintsq
  REAL(num)            :: dtsdx0, dtsdy0, dtsdz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(6)   ::  sx(-2:3),  sy(-2:3),  sz(-2:3)
  REAL(num), DIMENSION(6)   :: sx0(-2:3), sy0(-2:3), sz0(-2:3)
  REAL(num), DIMENSION(6)   :: dsx(-2:3), dsy(-2:3), dsz(-2:3)
  REAL(num)                 :: sdxim1,            sdxi
  REAL(num), DIMENSION(6)   :: sdyjm1(-2:3),      sdyj(-2:3) 
  REAL(num), DIMENSION(6,6) :: sdzkm1(-2:3,-2:3), sdzk(-2:3,-2:3)
  INTEGER :: iixp0, ijxp0, ikxp0, iixp, ijxp, ikxp, ip, dix, diy, diz, i, j, k, ic,   &
  jc, kc, ixmin, ixmax, iymin, iymax, izmin, izmax

  ! PARAMETER INIT
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dy*dz)
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)
  clghtisq = 1.0_num/clight**2
  dtsdz0 = dt*dzi
  sdxi=0.0_num
  sdxim1=0.0_num
  sdyj=0.0_num
  sdyjm1=0.0_num
  sdzk = 0.0_num
  sdzkm1 = 0.0_num
!$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(-2:3), sy(-2:3), sz(-2:3), sdxi, sdxim1, &
!$acc&                         sx0(-2:3), sy0(-2:3), sz0(-2:3), &
!$acc&                         dsx(-2:3), dsy(-2:3), dsz(-2:3), &
!$acc&                         sdyj(-2:3), sdyjm1(-2:3), &
!$acc&                         sdzk(-2:3,-2:3), sdzkm1(-2:3,-2:3) )
  DO ip=1, np
    sx = 0.0_num
    sy = 0.0_num
    sz = 0.0_num
    sx0 = 0.0_num
    sy0 = 0.0_num
    sz0 = 0.0_num
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqy = wq*invdtdy
    wqz = wq*invdtdz
    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ijxp0=floor(y)
    ikxp0=floor(z)
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    yint=y-ijxp0
    zint=z-ikxp0
    ! --- computes coefficients for node centered quantities
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx0(-1) = onesixth*oxintsq*oxint
    sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
    sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
    sx0( 2) = onesixth*xintsq*xint
    oyint = 1.0_num-yint
    yintsq = yint*yint
    oyintsq = oyint*oyint
    sy0(-1) = onesixth*oyintsq*oyint
    sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
    sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
    sy0( 2) = onesixth*yintsq*yint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz0(-1) = onesixth*ozintsq*ozint
    sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
    sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
    sz0( 2) = onesixth*zintsq*zint
    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ijxp=floor(yold)
    ikxp=floor(zold)
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    yint = yold-ijxp
    zint = zold-ikxp
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diy = ijxp-ijxp0
    diz = ikxp-ikxp0
    ! --- computes coefficients for quantities centered between nodes
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1+dix) = onesixth*oxintsq*oxint
    sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
    sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
    sx( 2+dix) = onesixth*xintsq*xint
    oyint = 1.0_num-yint
    yintsq = yint*yint
    oyintsq = oyint*oyint
    sy(-1+diy) = onesixth*oyintsq*oyint
    sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
    sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
    sy( 2+diy) = onesixth*yintsq*yint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1+diz) = onesixth*ozintsq*ozint
    sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
    sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
    sz( 2+diz) = onesixth*zintsq*zint
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsy = sy - sy0
    dsz = sz - sz0

    ! --- computes min/max positions of current contributions
    ixmin = min(0, dix)-1
    ixmax = max(0, dix)+2
    iymin = min(0, diy)-1
    iymax = max(0, diy)+2
    izmin = min(0, diz)-1
    izmax = max(0, diz)+2

    ! --- add current contributions
    DO k=izmin, izmax
      DO j=iymin, iymax
        DO i=ixmin, ixmax
          ic = iixp0+i
          jc = ijxp0+j
          kc = ikxp0+k
          IF(i<ixmax) THEN
            sdxi  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) +              &
            (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
            IF (i>ixmin) sdxi = sdxi + sdxim1
            !$acc atomic update
            jx(ic, jc, kc) = jx(ic, jc, kc) + sdxi
          END IF
          IF(j<iymax) THEN
            sdyj(i)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) +              &
            (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
            IF (j>iymin) sdyj(i) = sdyj(i) + sdyjm1(i)
            !$acc atomic update
            jy(ic, jc, kc) = jy(ic, jc, kc) + sdyj(i)
          END IF
          IF(k<izmax) THEN
            sdzk(i,j)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) +              &
            (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
            IF (k>izmin) sdzk(i,j) = sdzk(i,j) + sdzkm1(i,j)
            !$acc atomic update
            jz(ic, jc, kc) = jz(ic, jc, kc) + sdzk(i,j)
          END IF
          sdxim1 = sdxi
        END DO
        sdyjm1 = sdyj
      END DO
      sdzkm1 = sdzk
    END DO
  END DO
!$acc end loop
!$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_3_3_3

! ________________________________________________________________________________________
!> @brief
!> Esirkepov current deposition algorithm for linear, quadratic or cubic splines
!
!> @details
!> This subroutine can be used for several orders
!> WARNING: Highly unoptimized routine ---> USE INLINED ROUTINE
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
! ________________________________________________________________________________________
SUBROUTINE pxr_depose_jxjyjz_esirkepov_n( jx, jx_nguard, jx_nvalid, jy, jy_nguard,    &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, nox, noy, noz, l_particles_weight,                  &
  l4symtry)!#do not wrap
USE constants, ONLY: clight
USE picsar_precision, ONLY: idp, lp, num

  IMPLICIT NONE

  INTEGER(idp)             :: np, nox, noy, noz
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                :: q, dt, dx, dy, dz, xmin, ymin, zmin
  INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), jy_nguard(3), jy_nvalid(3), &
  jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,           &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1,                                          &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,           &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1,                                          &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,           &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,                                          &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  REAL(num) :: dxi, dyi, dzi, dtsdx, dtsdy, dtsdz, xint, yint, zint
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: sdx, sdy, sdz
  REAL(num) :: xold, yold, zold, x, y, z, wq, wqx, wqy, wqz, vx, vy, vz, dts2dx,      &
  dts2dy, dts2dz
  REAL(num) :: invvol, invdtdx, invdtdy, invdtdz
  REAL(num) :: oxint, oyint, ozint, xintsq, yintsq, zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: dtsdx0, dtsdy0, dtsdz0, dts2dx0, dts2dy0, dts2dz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(:), ALLOCATABLE :: sx, sx0, dsx
  REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
  REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
  INTEGER(idp) :: iixp0, ijxp0, ikxp0, iixp, ijxp, ikxp, ip, dix, diy, diz, i, j, k,  &
  ic, jc, kc, ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx,        &
  ndtody, ndtodz, xl, xu, yl, yu, zl, zu
  LOGICAL(lp)  :: l_particles_weight, l4symtry

  ! PARAMETER INIT
  ndtodx = int(clight*dt/dx)
  ndtody = int(clight*dt/dy)
  ndtodz = int(clight*dt/dz)
  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  dts2dx0 = 0.5_num*dtsdx0
  dts2dy0 = 0.5_num*dtsdy0
  dts2dz0 = 0.5_num*dtsdz0
  invvol = 1.0_num/(dx*dy*dz)
  invdtdx = 1.0_num/(dt*dy*dz)
  invdtdy = 1.0_num/(dt*dx*dz)
  invdtdz = 1.0_num/(dt*dx*dy)

  xl = -int(nox/2)-1-ndtodx
  xu = int((nox+1)/2)+1+ndtodx
  yl = -int(noy/2)-1-ndtody
  yu = int((noy+1)/2)+1+ndtody
  zl = -int(noz/2)-1-ndtodz
  zu = int((noz+1)/2)+1+ndtodz

  ALLOCATE(sdx(xl:xu, yl:yu, zl:zu), sdy(xl:xu, yl:yu, zl:zu), sdz(xl:xu, yl:yu,      &
  zl:zu))
  ALLOCATE(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
  ALLOCATE(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
  ALLOCATE(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdy=0.0_num;sdz=0.0_num


  DO ip=1, np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- applies 4-fold symmetry
    IF (l4symtry) THEN
      x=abs(x)
      y=abs(y)
      xold=abs(xold)
      yold=abs(yold)
      vx = (x-xold)/dtsdx0
      vy = (y-yold)/dtsdy0
    END IF
    !computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdy = dtsdy0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dy = dts2dy0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    y=yold
    z=zold
    DO icell = 1, ncells
      xold = x
      yold = y
      zold = z
      x = x+dtsdx*vx
      y = y+dtsdy*vy
      z = z+dtsdz*vz

      ! --- computes particles "weights"
      IF (l_particles_weight) THEN
        wq=q*w(ip)
      ELSE
        wq=q*w(1)
      END IF
      wqx = wq*invdtdx
      wqy = wq*invdtdy
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      ! --- (different for odd/even spline orders)
      IF (nox==2*(nox/2)) THEN
        iixp0=nint(x)
      ELSE
        iixp0=floor(x)
      END IF
      IF (noy==2*(noy/2)) THEN
        ijxp0=nint(y)
      ELSE
        ijxp0=floor(y)
      END IF
      IF (noz==2*(noz/2)) THEN
        ikxp0=nint(z)
      ELSE
        ikxp0=floor(z)
      END IF
      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      yint=y-ijxp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
      SELECT CASE(nox)
      CASE(0)
        sx0( 0) = 1.0_num
      CASE(1)
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
      CASE(2)
        xintsq = xint*xint
        sx0(-1) = 0.5_num*(0.5_num-xint)**2
        sx0( 0) = 0.75_num-xintsq
        sx0( 1) = 0.5_num*(0.5_num+xint)**2
      CASE(3)
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(-1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx0( 2) = onesixth*xintsq*xint
      END SELECT

      SELECT CASE(noy)
      CASE(0)
        sy0( 0) = 1.0_num
      CASE(1)
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
      CASE(2)
        yintsq = yint*yint
        sy0(-1) = 0.5_num*(0.5_num-yint)**2
        sy0( 0) = 0.75_num-yintsq
        sy0( 1) = 0.5_num*(0.5_num+yint)**2
      CASE(3)
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy0( 2) = onesixth*yintsq*yint
      END SELECT

      SELECT CASE(noz)
      CASE(0)
        sz0( 0) = 1.0_num
      CASE(1)
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint
      CASE(2)
        zintsq = zint*zint
        sz0(-1) = 0.5_num*(0.5_num-zint)**2
        sz0( 0) = 0.75_num-zintsq
        sz0( 1) = 0.5_num*(0.5_num+zint)**2
      CASE(3)
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz0( 2) = onesixth*zintsq*zint
      END SELECT

      ! --- finds node of cell containing particles for old positions
      ! --- (different for odd/even spline orders)
      IF (nox==2*(nox/2)) THEN
        iixp=nint(xold)
      ELSE
        iixp=floor(xold)
      END IF
      IF (noy==2*(noy/2)) THEN
        ijxp=nint(yold)
      ELSE
        ijxp=floor(yold)
      END IF
      IF (noz==2*(noz/2)) THEN
        ikxp=nint(zold)
      ELSE
        ikxp=floor(zold)
      END IF

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      yint = yold-ijxp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diy = ijxp-ijxp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients
      ! --- (needed because of different dix and diz for each particle)
      sx=0.0_num;sy=0.0_num;sz=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      SELECT CASE(nox)
      CASE(0)
        sx( 0+dix) = 1.0_num
      CASE(1)
        sx( 0+dix) = 1.0_num-xint
        sx( 1+dix) = xint
      CASE(2)
        xintsq = xint*xint
        sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
        sx( 0+dix) = 0.75_num-xintsq
        sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
      CASE(3)
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1+dix) = onesixth*oxintsq*oxint
        sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx( 2+dix) = onesixth*xintsq*xint
      END SELECT

      SELECT CASE(noy)
      CASE(0)
        sy( 0+diy) = 1.0_num
      CASE(1)
        sy( 0+diy) = 1.0_num-yint
        sy( 1+diy) = yint
      CASE(2)
        yintsq = yint*yint
        sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
        sy( 0+diy) = 0.75_num-yintsq
        sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
      CASE(3)
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1+diy) = onesixth*oyintsq*oyint
        sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy( 2+diy) = onesixth*yintsq*yint
      END SELECT

      SELECT CASE(noz)
      CASE(0)
        sz( 0+diz) = 1.0_num
      CASE(1)
        sz( 0+diz) = 1.0_num-zint
        sz( 1+diz) = zint
      CASE(2)
        zintsq = zint*zint
        sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
        sz( 0+diz) = 0.75_num-zintsq
        sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
      CASE(3)
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1+diz) = onesixth*ozintsq*ozint
        sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz( 2+diz) = onesixth*zintsq*zint
      END SELECT

      ! --- computes coefficients difference
      dsx = sx - sx0
      dsy = sy - sy0
      dsz = sz - sz0

      ! --- computes min/max positions of current contributions
      ixmin = min(0_idp, dix)-int(nox/2)
      ixmax = max(0_idp, dix)+int((nox+1)/2)
      iymin = min(0_idp, diy)-int(noy/2)
      iymax = max(0_idp, diy)+int((noy+1)/2)
      izmin = min(0_idp, diz)-int(noz/2)
      izmax = max(0_idp, diz)+int((noz+1)/2)
      ! --- add current contributions
      DO k=izmin, izmax
        DO j=iymin, iymax
          DO i=ixmin, ixmax
            ic = iixp0+i
            jc = ijxp0+j
            kc = ikxp0+k
            IF(i<ixmax) THEN
              sdx(i, j, k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) +            &
              (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
              IF (i>ixmin) sdx(i, j, k)=sdx(i, j, k)+sdx(i-1, j, k)
              jx(ic, jc, kc) = jx(ic, jc, kc) + sdx(i, j, k)
            END IF
            IF(j<iymax) THEN
              sdy(i, j, k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) +            &
              (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
              IF (j>iymin) sdy(i, j, k)=sdy(i, j, k)+sdy(i, j-1, k)
              jy(ic, jc, kc) = jy(ic, jc, kc) + sdy(i, j, k)
            END IF
            IF(k<izmax) THEN
              sdz(i, j, k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) +            &
              (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
              IF (k>izmin) sdz(i, j, k)=sdz(i, j, k)+sdz(i, j, k-1)
              jz(ic, jc, kc) = jz(ic, jc, kc) + sdz(i, j, k)
            END IF
          END DO
        END DO
      END DO

    END DO
  END DO

  DEALLOCATE(sdx, sdy, sdz, sx, sx0, dsx, sy, sy0, dsy, sz, sz0, dsz)

  RETURN
END SUBROUTINE pxr_depose_jxjyjz_esirkepov_n

#if defined (DEV)
! ________________________________________________________________________________________
!> @brief
!> Warp fonction for esirkepov
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine warp_depose_jxjyjz_esirkepov_n(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp,  &
  w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox,   &
  noy, noz, l_particles_weight, l4symtry)
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  implicit none
  integer(8) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard
  real(kind=8), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                   &
  -nzguard:nz+nzguard), intent(in out) :: jx
  real(kind=8), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                   &
  -nzguard:nz+nzguard), intent(in out) :: jy
  real(kind=8), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                   &
  -nzguard:nz+nzguard), intent(in out) :: jz
  real(kind=8), dimension(np) :: xp, yp, zp, uxp, uyp, uzp, v, w
  real(kind=8) :: q, dt, dx, dy, dz, xmin, ymin, zmin, clghtisq, usq, gaminv
  logical(8) :: l_particles_weight, l4symtry

  real(kind=8) :: dxi, dyi, dzi, dtsdx, dtsdy, dtsdz, xint, yint, zint
  real(kind=8), dimension(:, :, :), allocatable :: sdx, sdy, sdz
  real(kind=8) :: xold, yold, zold, xmid, ymid, zmid, x, y, z, wq, wqx, wqy, wqz,     &
  tmp, vx, vy, vz, dts2dx, dts2dy, dts2dz, s1x, s2x, s1y, s2y, s1z, s2z, invvol,      &
  invdtdx, invdtdy, invdtdz, oxint, oyint, ozint, xintsq, yintsq, zintsq, oxintsq,    &
  oyintsq, ozintsq, dtsdx0, dtsdy0, dtsdz0, dts2dx0, dts2dy0, dts2dz0
  real(kind=8), parameter :: onesixth=1./6., twothird=2./3.
  real(kind=8), DIMENSION(:), allocatable :: sx, sx0, dsx
  real(kind=8), DIMENSION(:), allocatable :: sy, sy0, dsy
  real(kind=8), DIMENSION(:), allocatable :: sz, sz0, dsz
  integer(8) :: iixp0, ijxp0, ikxp0, iixp, ijxp, ikxp, ip, dix, diy, diz, idx, idy,   &
  idz, i, j, k, ic, jc, kc, ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells,  &
  ndtodx, ndtody, ndtodz, xl, xu, yl, yu, zl, zu
  ndtodx = int(clight*dt/dx)
  ndtody = int(clight*dt/dy)
  ndtodz = int(clight*dt/dz)
  xl = -int(nox/2)-1-ndtodx
  xu = int((nox+1)/2)+1+ndtodx
  yl = -int(noy/2)-1-ndtody
  yu = int((noy+1)/2)+1+ndtody
  zl = -int(noz/2)-1-ndtodz
  zu = int((noz+1)/2)+1+ndtodz
  allocate(sdx(xl:xu, yl:yu, zl:zu), sdy(xl:xu, yl:yu, zl:zu), sdz(xl:xu, yl:yu,      &
  zl:zu))
  allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
  allocate(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
  allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))
  clghtisq = 1.0_num/clight**2
  sx0=0.;sy0=0.;sz0=0.
  sdx=0.;sdy=0.;sdz=0.

  dxi = 1./dx
  dyi = 1./dy
  dzi = 1./dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  dts2dx0 = 0.5*dtsdx0
  dts2dy0 = 0.5*dtsdy0
  dts2dz0 = 0.5*dtsdz0
  invvol = 1./(dx*dy*dz)
  invdtdx = 1./(dt*dy*dz)
  invdtdy = 1./(dt*dx*dz)
  invdtdz = 1./(dt*dx*dy)

  do ip=1, np

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq
    gaminv = 1.0_num/sqrt(1.0_num + usq)
    ! --- computes velocity
    vx = uxp(ip)*gaminv
    vy = uyp(ip)*gaminv
    vz = uzp(ip)*gaminv

    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz

    ! --- applies 4-fold symmetry
    if (l4symtry) then
      x=abs(x)
      y=abs(y)
      xold=abs(xold)
      yold=abs(yold)
      vx = (x-xold)/dtsdx0
      vy = (y-yold)/dtsdy0
    end if

    ! computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdy = dtsdy0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dy = dts2dy0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    y=yold
    z=zold

    do icell = 1, ncells

      xold = x
      yold = y
      zold = z

      x = x+dtsdx*vx
      y = y+dtsdy*vy
      z = z+dtsdz*vz

      ! --- computes particles "weights"
      if (l_particles_weight) then
        wq=q*w(ip)
      else
        wq=q*w(1)
      end if
      wqx = wq*invdtdx
      wqy = wq*invdtdy
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      ! --- (different for odd/even spline orders)
      if (nox==2*(nox/2)) then
        iixp0=nint(x)
      else
        iixp0=floor(x)
      end if
      if (noy==2*(noy/2)) then
        ijxp0=nint(y)
      else
        ijxp0=floor(y)
      end if
      if (noz==2*(noz/2)) then
        ikxp0=nint(z)
      else
        ikxp0=floor(z)
      end if

      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      yint=y-ijxp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
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

      select case(noy)
      case(0)
        sy0( 0) = 1.
      case(1)
        sy0( 0) = 1.-yint
        sy0( 1) = yint
      case(2)
        yintsq = yint*yint
        sy0(-1) = 0.5*(0.5-yint)**2
        sy0( 0) = 0.75-yintsq
        sy0( 1) = 0.5*(0.5+yint)**2
      case(3)
        oyint = 1.-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.-yint/2)
        sy0( 1) = twothird-oyintsq*(1.-oyint/2)
        sy0( 2) = onesixth*yintsq*yint
      end select

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
      if (noy==2*(noy/2)) then
        ijxp=nint(yold)
      else
        ijxp=floor(yold)
      end if
      if (noz==2*(noz/2)) then
        ikxp=nint(zold)
      else
        ikxp=floor(zold)
      end if

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      yint = yold-ijxp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diy = ijxp-ijxp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients
      ! (needed because of different dix and diz for each particle)
      sx=0.;sy=0.;sz=0.

      ! --- computes coefficients for quantities centered between nodes
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

      select case(noy)
      case(0)
        sy( 0+diy) = 1.
      case(1)
        sy( 0+diy) = 1.-yint
        sy( 1+diy) = yint
      case(2)
        yintsq = yint*yint
        sy(-1+diy) = 0.5*(0.5-yint)**2
        sy( 0+diy) = 0.75-yintsq
        sy( 1+diy) = 0.5*(0.5+yint)**2
      case(3)
        oyint = 1.-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1+diy) = onesixth*oyintsq*oyint
        sy( 0+diy) = twothird-yintsq*(1.-yint/2)
        sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
        sy( 2+diy) = onesixth*yintsq*yint
      end select

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
      dsy = sy - sy0
      dsz = sz - sz0

      ! --- computes min/max positions of current contributions
      ixmin = min(0_idp, dix)-int(nox/2)
      ixmax = max(0_idp, dix)+int((nox+1)/2)
      iymin = min(0_idp, diy)-int(noy/2)
      iymax = max(0_idp, diy)+int((noy+1)/2)
      izmin = min(0_idp, diz)-int(noz/2)
      izmax = max(0_idp, diz)+int((noz+1)/2)

      ! --- add current contributions
      do k=izmin, izmax
        do j=iymin, iymax
          do i=ixmin, ixmax
            ic = iixp0+i
            jc = ijxp0+j
            kc = ikxp0+k
            if(i<ixmax) then
              sdx(i, j, k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) +               &
              (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
              if (i>ixmin) sdx(i, j, k)=sdx(i, j, k)+sdx(i-1, j, k)
              jx(ic, jc, kc) = jx(ic, jc, kc) + sdx(i, j, k)
            end if
            if(j<iymax) then
              sdy(i, j, k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) +               &
              (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              if (j>iymin) sdy(i, j, k)=sdy(i, j, k)+sdy(i, j-1, k)
              jy(ic, jc, kc) = jy(ic, jc, kc) + sdy(i, j, k)
            end if
            if(k<izmax) then
              sdz(i, j, k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) +               &
              (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
              if (k>izmin) sdz(i, j, k)=sdz(i, j, k)+sdz(i, j, k-1)
              jz(ic, jc, kc) = jz(ic, jc, kc) + sdz(i, j, k)
            end if
          end do
        end do
      end do

    end do

  end do

  deallocate(sdx, sdy, sdz, sx, sx0, dsx, sy, sy0, dsy, sz, sz0, dsz)

  return
END SUBROUTINE warp_depose_jxjyjz_esirkepov_n
#endif

#if defined (DEV)
! ________________________________________________________________________________________
!> @brief
!> Warp fonction for esirkepov
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine picsar_depose_jxjyjz_esirkepov_n(cj, np, xp, yp, zp, uxp, uyp, uzp,        &
  gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard,         &
  nzguard, nox, noy, noz, l_particles_weight, l4symtry)
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp
  implicit none
  integer(8) :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard, nzguard
  real(kind=8), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                   &
  -nzguard:nz+nzguard, 3), intent(in out) :: cj
  real(kind=8), dimension(np) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  real(kind=8) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  logical(8) :: l_particles_weight, l4symtry

  real(kind=8) :: dxi, dyi, dzi, dtsdx, dtsdy, dtsdz, xint, yint, zint
  real(kind=8), dimension(:, :, :), allocatable :: sdx, sdy, sdz
  real(kind=8) :: xold, yold, zold, xmid, ymid, zmid, x, y, z, wq, wqx, wqy, wqz,     &
  tmp, vx, vy, vz, dts2dx, dts2dy, dts2dz, s1x, s2x, s1y, s2y, s1z, s2z, invvol,      &
  invdtdx, invdtdy, invdtdz, oxint, oyint, ozint, xintsq, yintsq, zintsq, oxintsq,    &
  oyintsq, ozintsq, dtsdx0, dtsdy0, dtsdz0, dts2dx0, dts2dy0, dts2dz0
  real(kind=8), parameter :: onesixth=1./6., twothird=2./3.
  real(kind=8), DIMENSION(:), allocatable :: sx, sx0, dsx
  real(kind=8), DIMENSION(:), allocatable :: sy, sy0, dsy
  real(kind=8), DIMENSION(:), allocatable :: sz, sz0, dsz
  integer(8) :: iixp0, ijxp0, ikxp0, iixp, ijxp, ikxp, ip, dix, diy, diz, idx, idy,   &
  idz, i, j, k, ic, jc, kc, ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells,  &
  ndtodx, ndtody, ndtodz, xl, xu, yl, yu, zl, zu
  ndtodx = int(clight*dt/dx)
  ndtody = int(clight*dt/dy)
  ndtodz = int(clight*dt/dz)
  xl = -int(nox/2)-1-ndtodx
  xu = int((nox+1)/2)+1+ndtodx
  yl = -int(noy/2)-1-ndtody
  yu = int((noy+1)/2)+1+ndtody
  zl = -int(noz/2)-1-ndtodz
  zu = int((noz+1)/2)+1+ndtodz
  allocate(sdx(xl:xu, yl:yu, zl:zu), sdy(xl:xu, yl:yu, zl:zu), sdz(xl:xu, yl:yu,      &
  zl:zu))
  allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
  allocate(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
  allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

  sx0=0.;sy0=0.;sz0=0.
  sdx=0.;sdy=0.;sdz=0.

  dxi = 1./dx
  dyi = 1./dy
  dzi = 1./dz
  dtsdx0 = dt*dxi
  dtsdy0 = dt*dyi
  dtsdz0 = dt*dzi
  dts2dx0 = 0.5*dtsdx0
  dts2dy0 = 0.5*dtsdy0
  dts2dz0 = 0.5*dtsdz0
  invvol = 1./(dx*dy*dz)
  invdtdx = 1./(dt*dy*dz)
  invdtdy = 1./(dt*dx*dz)
  invdtdz = 1./(dt*dx*dy)

  do ip=1, np

    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)

    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz

    ! --- applies 4-fold symmetry
    if (l4symtry) then
      x=abs(x)
      y=abs(y)
      xold=abs(xold)
      yold=abs(yold)
      vx = (x-xold)/dtsdx0
      vy = (y-yold)/dtsdy0
    end if

    ! computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdy = dtsdy0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dy = dts2dy0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    y=yold
    z=zold

    do icell = 1, ncells

      xold = x
      yold = y
      zold = z

      x = x+dtsdx*vx
      y = y+dtsdy*vy
      z = z+dtsdz*vz

      ! --- computes particles "weights"
      if (l_particles_weight) then
        wq=q*w(ip)
      else
        wq=q*w(1)
      end if
      wqx = wq*invdtdx
      wqy = wq*invdtdy
      wqz = wq*invdtdz

      ! --- finds node of cell containing particles for current positions
      ! --- (different for odd/even spline orders)
      if (nox==2*(nox/2)) then
        iixp0=nint(x)
      else
        iixp0=floor(x)
      end if
      if (noy==2*(noy/2)) then
        ijxp0=nint(y)
      else
        ijxp0=floor(y)
      end if
      if (noz==2*(noz/2)) then
        ikxp0=nint(z)
      else
        ikxp0=floor(z)
      end if

      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      yint=y-ijxp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
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

      select case(noy)
      case(0)
        sy0( 0) = 1.
      case(1)
        sy0( 0) = 1.-yint
        sy0( 1) = yint
      case(2)
        yintsq = yint*yint
        sy0(-1) = 0.5*(0.5-yint)**2
        sy0( 0) = 0.75-yintsq
        sy0( 1) = 0.5*(0.5+yint)**2
      case(3)
        oyint = 1.-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.-yint/2)
        sy0( 1) = twothird-oyintsq*(1.-oyint/2)
        sy0( 2) = onesixth*yintsq*yint
      end select

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
      if (noy==2*(noy/2)) then
        ijxp=nint(yold)
      else
        ijxp=floor(yold)
      end if
      if (noz==2*(noz/2)) then
        ikxp=nint(zold)
      else
        ikxp=floor(zold)
      end if

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      yint = yold-ijxp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diy = ijxp-ijxp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients
      ! --- (needed because of different dix and diz for each particle)
      sx=0.;sy=0.;sz=0.

      ! --- computes coefficients for quantities centered between nodes
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

      select case(noy)
      case(0)
        sy( 0+diy) = 1.
      case(1)
        sy( 0+diy) = 1.-yint
        sy( 1+diy) = yint
      case(2)
        yintsq = yint*yint
        sy(-1+diy) = 0.5*(0.5-yint)**2
        sy( 0+diy) = 0.75-yintsq
        sy( 1+diy) = 0.5*(0.5+yint)**2
      case(3)
        oyint = 1.-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1+diy) = onesixth*oyintsq*oyint
        sy( 0+diy) = twothird-yintsq*(1.-yint/2)
        sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
        sy( 2+diy) = onesixth*yintsq*yint
      end select

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
      dsy = sy - sy0
      dsz = sz - sz0

      ! --- computes min/max positions of current contributions
      ixmin = min(0_idp, dix)-int(nox/2)
      ixmax = max(0_idp, dix)+int((nox+1)/2)
      iymin = min(0_idp, diy)-int(noy/2)
      iymax = max(0_idp, diy)+int((noy+1)/2)
      izmin = min(0_idp, diz)-int(noz/2)
      izmax = max(0_idp, diz)+int((noz+1)/2)

      ! --- add current contributions
      do k=izmin, izmax
        do j=iymin, iymax
          do i=ixmin, ixmax
            ic = iixp0+i
            jc = ijxp0+j
            kc = ikxp0+k
            if(i<ixmax) then
              sdx(i, j, k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) +               &
              (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
              if (i>ixmin) sdx(i, j, k)=sdx(i, j, k)+sdx(i-1, j, k)
              cj(ic, jc, kc, 1) = cj(ic, jc, kc, 1) + sdx(i, j, k)
            end if
            if(j<iymax) then
              sdy(i, j, k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) +               &
              (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              if (j>iymin) sdy(i, j, k)=sdy(i, j, k)+sdy(i, j-1, k)
              cj(ic, jc, kc, 2) = cj(ic, jc, kc, 2) + sdy(i, j, k)
            end if
            if(k<izmax) then
              sdz(i, j, k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) +               &
              (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
              if (k>izmin) sdz(i, j, k)=sdz(i, j, k)+sdz(i, j, k-1)
              cj(ic, jc, kc, 3) = cj(ic, jc, kc, 3) + sdz(i, j, k)
            end if
          end do
        end do
      end do

    end do

  end do

  deallocate(sdx, sdy, sdz, sx, sx0, dsx, sy, sy0, dsy, sz, sz0, dsz)

  return
END SUBROUTINE picsar_depose_jxjyjz_esirkepov_n
#endif

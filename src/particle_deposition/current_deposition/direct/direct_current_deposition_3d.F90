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
! DIRECT_CURRENT_DEPOSITION_3D.F90
!
! Developers
! Henri Vincenti, ! Mathieu Lobet
!
! Description:
! This file contains subroutines for the direct current deposition in 3D.
!
! List of subroutines:
!
! Classical non-optimized:
! - depose_jxjyjz_scalar_1_1_1
! - depose_jxjyjz_scalar_2_2_2
! - depose_jxjyjz_scalar_3_3_3
!
! Classical parallel/vectorized
! - depose_jxjyjz_vecHVv2_1_1_1
! - depose_jxjyjz_vecHVv2_2_2_2
! - depose_jxjyjz_vecHVv2_3_3_3
! - depose_jxjyjz_vecHVv3_3_3_3
!
! Subroutines with no reduction inside
! - depose_jxjyjz_vecHV_vnr_1_1_1
! - depose_jxjyjz_vecHV_vnr_2_2_2
! - depose_jxjyjz_vecHV_vnr_3_3_3
!
! - current_reduction_1_1_1
! - current_reduction_2_2_2
! - current_reduction_3_3_3
!
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Order 1 3D scalar direct current deposition routine (rho*v)
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
!> @param[in] yp 1D array of x-coordinates of particles
!> @param[in] zp 1D array of x-coordinates of particles
!> @param[in] uxp 1D array of ux-velocity components of particles
!> @param[in] uyp 1D array of ux-velocity components of particles
!> @param[in] uzp 1D array of ux-velocity components of particles
!> @param[in] gaminv 1D array of the inverse 1/gamma-factor of particles
!> @param[in] w 1D array of the weghts of particles
!> @param[in] q charge of current species (scalar)
!> @param[in] xmin x-minimum boundary of current tile
!> @param[in] ymin y-minimum boundary of current tile
!> @param[in] zmin z-minimum boundary of current tile
!> @param[in] dt time step (scalar)
!> @param[in] dx mesh size along x (scalar)
!> @param[in] dy mesh size along y (scalar)
!> @param[in] dz mesh size along z (scalar)
!> @param[inout] jx x-current component (3D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 3 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jx array (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!>  (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 3 integers)
!> @warning arrays jx, jy, jz should be set to 0 before entering this subroutine.
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_scalar_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp)             :: np
  LOGICAL(idp)             :: l_nodal
  REAL(num)                :: stagger_shift
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
  REAL(num)                :: dxi, dyi, dzi, xint, yint, zint
  REAL(num)                :: x, y, z, xmid, ymid, zmid, vx, vy, vz, invvol, dts2dx,  &
  dts2dy, dts2dz
  REAL(num)                :: wq, wqx, wqy, wqz, clightsq
  REAL(num), DIMENSION(2)  :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
  REAL(num), PARAMETER     :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp)             :: j, k, l, j0, k0, l0, ip

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
!  sx=0.0_num;sy=0.0_num;sz=0.0_num;
!  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;

  ! LOOP ON PARTICLES
  ! Prevent loop to vectorize (dependencies)
  !DIR$ NOVECTOR
!$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)
!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1) )
  DO ip=1, np

    ! --- computes position in  grid units at (n+1)
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
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
    ymid=y-dts2dy*vy
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j=floor(xmid)
    k=floor(ymid)
    l=floor(zmid)
    j0=floor(xmid-stagger_shift)
    k0=floor(ymid-stagger_shift)
    l0=floor(zmid-stagger_shift)

    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    yint = ymid-k
    zint = zmid-l
    sx( 0) = 1.0_num-xint
    sx( 1) = xint
    sy( 0) = 1.0_num-yint
    sy( 1) = yint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-stagger_shift-j0
    yint = ymid-stagger_shift-k0
    zint = zmid-stagger_shift-l0
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
    sy0( 0) = 1.0_num-yint
    sy0( 1) = yint
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
    ! - JX
    !$acc atomic update
    jx(j0, k, l  )    = jx(j0, k, l  )  +   sx0(0)*sy(0)*sz(0)*wqx
    !$acc atomic update
    jx(j0+1, k, l  )    = jx(j0+1, k, l  )  +   sx0(1)*sy(0)*sz(0)*wqx
    !$acc atomic update
    jx(j0, k+1, l  )    = jx(j0, k+1, l  )  +   sx0(0)*sy(1)*sz(0)*wqx
    !$acc atomic update
    jx(j0+1, k+1, l  )    = jx(j0+1, k+1, l  )  +   sx0(1)*sy(1)*sz(0)*wqx
    !$acc atomic update
    jx(j0, k, l+1)    = jx(j0, k, l+1)  +   sx0(0)*sy(0)*sz(1)*wqx
    !$acc atomic update
    jx(j0+1, k, l+1)    = jx(j0+1, k, l+1)  +   sx0(1)*sy(0)*sz(1)*wqx
    !$acc atomic update
    jx(j0, k+1, l+1)    = jx(j0, k+1, l+1)  +   sx0(0)*sy(1)*sz(1)*wqx
    !$acc atomic update
    jx(j0+1, k+1, l+1)    = jx(j0+1, k+1, l+1)  +   sx0(1)*sy(1)*sz(1)*wqx

    ! - JY
    !$acc atomic update
    jy(j, k0, l  )    = jy(j, k0, l  )  +   sx(0)*sy0(0)*sz(0)*wqy
    !$acc atomic update
    jy(j+1, k0, l  )    = jy(j+1, k0, l  )  +   sx(1)*sy0(0)*sz(0)*wqy
    !$acc atomic update
    jy(j, k0+1, l  )    = jy(j, k0+1, l  )  +   sx(0)*sy0(1)*sz(0)*wqy
    !$acc atomic update
    jy(j+1, k0+1, l  )    = jy(j+1, k0+1, l  )  +   sx(1)*sy0(1)*sz(0)*wqy
    !$acc atomic update
    jy(j, k0, l+1)    = jy(j, k0, l+1)  +   sx(0)*sy0(0)*sz(1)*wqy
    !$acc atomic update
    jy(j+1, k0, l+1)    = jy(j+1, k0, l+1)  +   sx(1)*sy0(0)*sz(1)*wqy
    !$acc atomic update
    jy(j, k0+1, l+1)    = jy(j, k0+1, l+1)  +   sx(0)*sy0(1)*sz(1)*wqy
    !$acc atomic update
    jy(j+1, k0+1, l+1)    = jy(j+1, k0+1, l+1)  +   sx(1)*sy0(1)*sz(1)*wqy

    ! - JZ
    !$acc atomic update
    jz(j, k, l0  )    = jz(j, k, l0  )  +   sx(0)*sy(0)*sz0(0)*wqz
    !$acc atomic update
    jz(j+1, k, l0  )    = jz(j+1, k, l0  )  +   sx(1)*sy(0)*sz0(0)*wqz
    !$acc atomic update
    jz(j, k+1, l0  )    = jz(j, k+1, l0  )  +   sx(0)*sy(1)*sz0(0)*wqz
    !$acc atomic update
    jz(j+1, k+1, l0  )    = jz(j+1, k+1, l0  )  +   sx(1)*sy(1)*sz0(0)*wqz
    !$acc atomic update
    jz(j, k, l0+1)    = jz(j, k, l0+1)  +   sx(0)*sy(0)*sz0(1)*wqz
    !$acc atomic update
    jz(j+1, k, l0+1)    = jz(j+1, k, l0+1)  +   sx(1)*sy(0)*sz0(1)*wqz
    !$acc atomic update
    jz(j, k+1, l0+1)    = jz(j, k+1, l0+1)  +   sx(0)*sy(1)*sz0(1)*wqz
    !$acc atomic update
    jz(j+1, k+1, l0+1)    = jz(j+1, k+1, l0+1)  +   sx(1)*sy(1)*sz0(1)*wqz
  END DO
!$acc end loop
!$acc end parallel

  RETURN
END SUBROUTINE depose_jxjyjz_scalar_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 1 3D vector direct current deposition routine (rho*v)
!> This versions have good performances on SIMD architectures
!> Providing that OpenMP 4.0 is available (Directive SIMD)
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!
!> @param[in] np Number of particles
!> @param[in] xp 1D array of x-coordinates of particles
!> @param[in] yp 1D array of x-coordinates of particles
!> @param[in] zp 1D array of x-coordinates of particles
!> @param[in] uxp 1D array of ux-velocity components of particles
!> @param[in] uyp 1D array of ux-velocity components of particles
!> @param[in] uzp 1D array of ux-velocity components of particles
!> @param[in] gaminv 1D array of the inverse 1/gamma-factor of particles
!> @param[in] w 1D array of the weghts of particles
!> @param[in] q charge of current species (scalar)
!> @param[in] xmin x-minimum boundary of current tile
!> @param[in] ymin y-minimum boundary of current tile
!> @param[in] zmin z-minimum boundary of current tile
!> @param[in] dt time step (scalar)
!> @param[in] dx mesh size along x (scalar)
!> @param[in] dy mesh size along y (scalar)
!> @param[in] dz mesh size along z (scalar)
!> @param[inout] jx x-current component (3D array)
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 3 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jx array (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!>  (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 3 integers)
!> @warning arrays jx, jy, jz should be set to 0 before entering this subroutine.
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHVv2_1_1_1( jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  ! ___ Parameter declaration ______________________________________
  INTEGER(idp)             :: np
  LOGICAL(idp)             :: l_nodal
  REAL(num)                :: stagger_shift
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
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  REAL(num)                :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num)                :: dxi, dyi, dzi
  REAL(num)                :: x, y, z, xmid, ymid, zmid, invvol, dts2dx, dts2dy,      &
  dts2dz

  INTEGER(idp)                    :: j, k, l, j0, k0, l0, ip, NCELLS, ic
  INTEGER(idp)                    :: n, nn, nv
  INTEGER(idp), DIMENSION(LVEC, 3) :: ICELL
  REAL(num), DIMENSION(LVEC)      :: sx, sy, sz, sx0, sy0, sz0, wqx, wqy, wqz
  REAL(num)                       :: wwx, wwy, wwz, wq, vx, vy, vz, wx, wx0, wy, wy0, &
  wz, wz0
  INTEGER(idp)                    :: jorig, korig, lorig
  INTEGER(idp)                    :: ncx, ncy, ncxy, ncz, ix, iy, iz
  ! Coefficients for the computation of the weights
  REAL(num), DIMENSION(8), PARAMETER   :: mx=(/1_num, 0_num, 1_num, 0_num, 1_num,     &
  0_num, 1_num, 0_num/)
  REAL(num), DIMENSION(8), PARAMETER   :: my=(/1_num, 1_num, 0_num, 0_num, 1_num,     &
  1_num, 0_num, 0_num/)
  REAL(num), DIMENSION(8), PARAMETER   :: mz=(/1_num, 1_num, 1_num, 1_num, 0_num,     &
  0_num, 0_num, 0_num/)
  REAL(num), DIMENSION(8), PARAMETER   :: sgn=(/-1_num, 1_num, 1_num, -1_num, 1_num,  &
  -1_num, -1_num, 1_num/)
  ! Positions of the nodes on the grids
  INTEGER(idp), DIMENSION(8), PARAMETER  :: mxoff=(/0_idp, 1_idp, 0_idp, 1_idp,       &
  0_idp, 1_idp, 0_idp, 1_idp/)
  INTEGER(idp), DIMENSION(8), PARAMETER  :: myoff=(/0_idp, 0_idp, 1_idp, 1_idp,       &
  0_idp, 0_idp, 1_idp, 1_idp/)
  INTEGER(idp), DIMENSION(8), PARAMETER  :: mzoff=(/0_idp, 0_idp, 0_idp, 0_idp,       &
  1_idp, 1_idp, 1_idp, 1_idp/)

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi
  sx=0.0_num;sy=0.0_num;sz=0.0_num
  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
  ! Find the maximal number of cells in each direction, so as to
  ! allocate the linearized 1D arrays
  jorig=-2; korig=-2;lorig=-2
  ncx = MAX( jx_nvalid(1), jy_nvalid(1), jz_nvalid(1) ) + 3
  ncy = MAX( jx_nvalid(2), jy_nvalid(2), jz_nvalid(2) ) + 3
  ncz = MAX( jx_nvalid(3), jy_nvalid(3), jz_nvalid(3) ) + 3
  ncxy = ncx*ncy
  NCELLS = ncx*ncy*ncz
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
  jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num;

  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED w:64, gaminv:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, sx, sy, sz)
    !IBM* ALIGN(64, sy0, sz0)
    !IBM* ALIGN(64, w, gaminv)
    !IBM* ALIGN(64, ICELL)
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

      ! Computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)

      ! --- computes particles weights
      wq=q*w(nn)*invvol
      wqx(n)=wq*vx
      wqy(n)=wq*vy
      wqz(n)=wq*vz

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx
      ymid=y-dts2dy*vy
      zmid=z-dts2dz*vz

      ! --- finds node of cell containing particles for current positions
      j=floor(xmid)
      k=floor(ymid)
      l=floor(zmid)
      j0=floor(xmid-stagger_shift)
      k0=floor(ymid-stagger_shift)
      l0=floor(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! --- computes set of coefficients for node centered quantities
      sx(n) = xmid-j
      sy(n) = ymid-k
      sz(n) = zmid-l

      ! --- computes set of coefficients for staggered quantities
      sx0(n) = xmid-stagger_shift-j0
      sy0(n) = ymid-stagger_shift-k0
      sz0(n) = zmid-stagger_shift-l0
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    DO n=1, MIN(LVEC, np-ip+1)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64, jxcells, jycells, jzcells, mx, my, mz)
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
        wx=-mx(nv)+sx(n)
        wx0=-mx(nv)+sx0(n)
        wy=-my(nv)+sy(n)
        wy0=-my(nv)+sy0(n)
        wz=-mz(nv)+sz(n)
        wz0=-mz(nv)+sz0(n)
        wwx=wx0*wy*wz*wqx(n)*sgn(nv)
        wwy=wx*wy0*wz*wqy(n)*sgn(nv)
        wwz=wx*wy*wz0*wqz(n)*sgn(nv)
        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! - JX
        jxcells(nv, ICELL(n, 1))=jxcells(nv, ICELL(n, 1))+wwx
        ! - JY
        jycells(nv, ICELL(n, 2))=jycells(nv, ICELL(n, 2))+wwy
        ! - JZ
        jzcells(nv, ICELL(n, 3))=jzcells(nv, ICELL(n, 3))+wwz
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! ----------------------------------------------------
  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
  ! ----------------------------------------------------

  ! Reduction of jxcells in jx
  ! Note: the bounds below make sure that we never go out-of-bound for jxcells,
  ! even when taking into account the additional 0/1 offsets in mxoff, myoff, mzoff
  DO iz = MAX(lorig, -jx_nguard(3)), MIN(lorig+ncz-1, jx_nvalid(3)+jx_nguard(3)-2)
    DO iy = MAX(korig, -jx_nguard(2)), MIN(korig+ncy-1, jx_nvalid(2)+jx_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix = MAX(jorig, -jx_nguard(1)), MIN(jorig+ncx-1, jx_nvalid(1)+jx_nguard(1)-2)
        ! Compute linearized index
        ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
        jx(ix+mxoff(1), iy+myoff(1), iz+mzoff(1)) = jx(ix+mxoff(1), iy+myoff(1),      &
        iz+mzoff(1)) + jxcells(1, ic)
        jx(ix+mxoff(2), iy+myoff(2), iz+mzoff(2)) = jx(ix+mxoff(2), iy+myoff(2),      &
        iz+mzoff(2)) + jxcells(2, ic)
        jx(ix+mxoff(3), iy+myoff(3), iz+mzoff(3)) = jx(ix+mxoff(3), iy+myoff(3),      &
        iz+mzoff(3)) + jxcells(3, ic)
        jx(ix+mxoff(4), iy+myoff(4), iz+mzoff(4)) = jx(ix+mxoff(4), iy+myoff(4),      &
        iz+mzoff(4)) + jxcells(4, ic)
        jx(ix+mxoff(5), iy+myoff(5), iz+mzoff(5)) = jx(ix+mxoff(5), iy+myoff(5),      &
        iz+mzoff(5)) + jxcells(5, ic)
        jx(ix+mxoff(6), iy+myoff(6), iz+mzoff(6)) = jx(ix+mxoff(6), iy+myoff(6),      &
        iz+mzoff(6)) + jxcells(6, ic)
        jx(ix+mxoff(7), iy+myoff(7), iz+mzoff(7)) = jx(ix+mxoff(7), iy+myoff(7),      &
        iz+mzoff(7)) + jxcells(7, ic)
        jx(ix+mxoff(8), iy+myoff(8), iz+mzoff(8)) = jx(ix+mxoff(8), iy+myoff(8),      &
        iz+mzoff(8)) + jxcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jycells in jy
  ! Note: the bounds below make sure that we never go out-of-bound for jxcells,
  ! even when taking into account the additional 0/1 offsets in mxoff, myoff, mzoff
  DO iz = MAX(lorig, -jy_nguard(3)), MIN(lorig+ncz-1, jy_nvalid(3)+jy_nguard(3)-2)
    DO iy = MAX(korig, -jy_nguard(2)), MIN(korig+ncy-1, jy_nvalid(2)+jy_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix = MAX(jorig, -jy_nguard(1)), MIN(jorig+ncx-1, jy_nvalid(1)+jy_nguard(1)-2)
        ! Compute linearized index
        ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
        jy(ix+mxoff(1), iy+myoff(1), iz+mzoff(1)) = jy(ix+mxoff(1), iy+myoff(1),      &
        iz+mzoff(1)) + jycells(1, ic)
        jy(ix+mxoff(2), iy+myoff(2), iz+mzoff(2)) = jy(ix+mxoff(2), iy+myoff(2),      &
        iz+mzoff(2)) + jycells(2, ic)
        jy(ix+mxoff(3), iy+myoff(3), iz+mzoff(3)) = jy(ix+mxoff(3), iy+myoff(3),      &
        iz+mzoff(3)) + jycells(3, ic)
        jy(ix+mxoff(4), iy+myoff(4), iz+mzoff(4)) = jy(ix+mxoff(4), iy+myoff(4),      &
        iz+mzoff(4)) + jycells(4, ic)
        jy(ix+mxoff(5), iy+myoff(5), iz+mzoff(5)) = jy(ix+mxoff(5), iy+myoff(5),      &
        iz+mzoff(5)) + jycells(5, ic)
        jy(ix+mxoff(6), iy+myoff(6), iz+mzoff(6)) = jy(ix+mxoff(6), iy+myoff(6),      &
        iz+mzoff(6)) + jycells(6, ic)
        jy(ix+mxoff(7), iy+myoff(7), iz+mzoff(7)) = jy(ix+mxoff(7), iy+myoff(7),      &
        iz+mzoff(7)) + jycells(7, ic)
        jy(ix+mxoff(8), iy+myoff(8), iz+mzoff(8)) = jy(ix+mxoff(8), iy+myoff(8),      &
        iz+mzoff(8)) + jycells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jzcells in jz
  ! Note: the bounds below make sure that we never go out-of-bound for jxcells,
  ! even when taking into account the additional 0/1 offsets in mxoff, myoff, mzoff
  DO iz = MAX(lorig, -jz_nguard(3)), MIN(lorig+ncz-1, jz_nvalid(3)+jz_nguard(3)-2)
    DO iy = MAX(korig, -jz_nguard(2)), MIN(korig+ncy-1, jz_nvalid(2)+jz_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix = MAX(jorig, -jz_nguard(1)), MIN(jorig+ncx-1, jz_nvalid(1)+jz_nguard(1)-2)
        ! Compute linearized index
        ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
        jz(ix+mxoff(1), iy+myoff(1), iz+mzoff(1)) = jz(ix+mxoff(1), iy+myoff(1),      &
        iz+mzoff(1)) + jzcells(1, ic)
        jz(ix+mxoff(2), iy+myoff(2), iz+mzoff(2)) = jz(ix+mxoff(2), iy+myoff(2),      &
        iz+mzoff(2)) + jzcells(2, ic)
        jz(ix+mxoff(3), iy+myoff(3), iz+mzoff(3)) = jz(ix+mxoff(3), iy+myoff(3),      &
        iz+mzoff(3)) + jzcells(3, ic)
        jz(ix+mxoff(4), iy+myoff(4), iz+mzoff(4)) = jz(ix+mxoff(4), iy+myoff(4),      &
        iz+mzoff(4)) + jzcells(4, ic)
        jz(ix+mxoff(5), iy+myoff(5), iz+mzoff(5)) = jz(ix+mxoff(5), iy+myoff(5),      &
        iz+mzoff(5)) + jzcells(5, ic)
        jz(ix+mxoff(6), iy+myoff(6), iz+mzoff(6)) = jz(ix+mxoff(6), iy+myoff(6),      &
        iz+mzoff(6)) + jzcells(6, ic)
        jz(ix+mxoff(7), iy+myoff(7), iz+mzoff(7)) = jz(ix+mxoff(7), iy+myoff(7),      &
        iz+mzoff(7)) + jzcells(7, ic)
        jz(ix+mxoff(8), iy+myoff(8), iz+mzoff(8)) = jz(ix+mxoff(8), iy+myoff(8),      &
        iz+mzoff(8)) + jzcells(8, ic)
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
END SUBROUTINE depose_jxjyjz_vecHVv2_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 1 3D vector current deposition routine (rho*v) with no reduction
!
!> @details
!> This versions have good performances on SIMD architectures
!> Providing that OpenMP 4.0 is available (Directive SIMD)
!> This subroutine is similar to depose_jxjyjz_vecHVv2_1_1_1
!> without the reduction process at the end
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] jxcells, jycells, jzcells temporary current arrays
!> @param[in] np Number of particles
!> @param[in] ncells number of cells in the tile
!> @param[in] xp 1D array of x-coordinates of particles
!> @param[in] yp 1D array of x-coordinates of particles
!> @param[in] zp 1D array of x-coordinates of particles
!> @param[in] uxp 1D array of ux-velocity components of particles
!> @param[in] uyp 1D array of ux-velocity components of particles
!> @param[in] uzp 1D array of ux-velocity components of particles
!> @param[in] gaminv 1D array of the inverse 1/gamma-factor of particles
!> @param[in] w 1D array of the weghts of particles
!> @param[in] q charge of current species (scalar)
!> @param[in] xmin x-minimum boundary of current tile
!> @param[in] ymin y-minimum boundary of current tile
!> @param[in] zmin z-minimum boundary of current tile
!> @param[in] dt time step (scalar)
!> @param[in] dx mesh size along x (scalar)
!> @param[in] dy mesh size along y (scalar)
!> @param[in] dz mesh size along z (scalar)
!> @param[in] nx number of cells along x (scalar)
!> @param[in] ny number of cells along y (scalar)
!> @param[in] nz number of cells along z (scalar)
!> @param[in] nxguard number of guard cells along x (scalar)
!> @param[in] nyguard number of guard cells along y (scalar)
!> @param[in] nzguard number of guard cells along z (scalar)
!> @param[in] ncx, ncy, ncz tile cell extended number (depends on the order)
!> @param[in] lvect vector length
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1(jxcells, jycells, jzcells, np, ncells, xp,   &
  yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz,    &
  nxguard, nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not wrap
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
 LOGICAL(idp)                      :: l_nodal
 REAL(num)                         :: stagger_shift
  INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
  INTEGER(idp), INTENT(IN)                      :: lvect
  REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np), INTENT(IN)          :: xp, yp, zp
  REAL(num), DIMENSION(np), INTENT(IN)          :: uxp, uyp, uzp, gaminv, w
  REAL(num), INTENT(IN)                         :: q, dt, dx, dy, dz, xmin, ymin,     &
  zmin
  REAL(num)                                     :: x, y, z, xmid, ymid, zmid
  INTEGER(isp)                                  :: j, k, l, j0, k0, l0, ip
  INTEGER(isp)                                  :: n, nn, nv
  REAL(num)                                     :: mx(1:8), my(1:8), mz(1:8),         &
  sgn(1:8)
  REAL(num)                                     :: invvol, dxi, dyi, dzi
  REAL(num)                                     :: dts2dx, dts2dy, dts2dz
  INTEGER(isp), DIMENSION(LVECT, 3)              :: ICELL
  REAL(num), DIMENSION(LVECT)                   :: sx, sy, sz
  REAL(num), DIMENSION(LVECT)                   :: sx0, sy0, sz0, wqx, wqy, wqz
  REAL(num) :: wwx, wwy, wwz, wq, vx, vy, vz, wx, wx0, wy, wy0, wz, wz0
  INTEGER(isp)                                  :: jorig, korig, lorig
  INTEGER(isp)                                  :: ncx, ncy, ncxy, ncz

  ! _____________________________________________
  ! Computation of the parameters
  ncxy=ncx*ncy

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi

  sx=0.0_num;sy=0.0_num;sz=0.0_num
  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num

  mx=(/1_num, 0_num, 1_num, 0_num, 1_num, 0_num, 1_num, 0_num/)
  my=(/1_num, 1_num, 0_num, 0_num, 1_num, 1_num, 0_num, 0_num/)
  mz=(/1_num, 1_num, 1_num, 1_num, 0_num, 0_num, 0_num, 0_num/)
  sgn=(/-1_num, 1_num, 1_num, -1_num, 1_num, -1_num, -1_num, 1_num/)

  jorig=-2
  korig=-2
  lorig=-2

  ! ____________________________________________
  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED w:64, gaminv:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, sx, sy, sz)
    !IBM* ALIGN(64, sy0, sz0)
    !IBM* ALIGN(64, w, gaminv)
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

      ! --- computes particles weights
      wq=q*w(nn)*invvol
      wqx(n)=wq*vx
      wqy(n)=wq*vy
      wqz(n)=wq*vz

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx
      ymid=y-dts2dy*vy
      zmid=z-dts2dz*vz

      ! --- finds node of cell containing particles for current positions
      j=floor(xmid)
      k=floor(ymid)
      l=floor(zmid)
      j0=floor(xmid-stagger_shift)
      k0=floor(ymid-stagger_shift)
      l0=floor(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! --- computes set of coefficients for node centered quantities
      sx(n) = xmid-j
      sy(n) = ymid-k
      sz(n) = zmid-l

      ! --- computes set of coefficients for staggered quantities
      sx0(n) = xmid-stagger_shift-j0
      sy0(n) = ymid-stagger_shift-k0
      sz0(n) = zmid-stagger_shift-l0
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    DO n=1, MIN(LVEC, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
      !DIR$ ASSUME_ALIGNED mx:64, my:64, mz:64, sgn:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(32, jxcells, jycells, jzcells)
      !IBM* ALIGN(32, mx, my, mz, sg)
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
        wx=-mx(nv)+sx(n)
        wx0=-mx(nv)+sx0(n)
        wy=-my(nv)+sy(n)
        wy0=-my(nv)+sy0(n)
        wz=-mz(nv)+sz(n)
        wz0=-mz(nv)+sz0(n)
        wwx=wx0*wy*wz*wqx(n)*sgn(nv)
        wwy=wx*wy0*wz*wqy(n)*sgn(nv)
        wwz=wx*wy*wz0*wqz(n)*sgn(nv)
        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! - JX
        jxcells(nv, ICELL(n, 1))=jxcells(nv, ICELL(n, 1))+wwx
        ! - JY
        jycells(nv, ICELL(n, 2))=jycells(nv, ICELL(n, 2))+wwy
        ! - JZ
        jzcells(nv, ICELL(n, 3))=jzcells(nv, ICELL(n, 3))+wwz
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  RETURN
END SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1

! ________________________________________________________________________________________
!> @brief
!> Order 2 3D scalar current deposition routine (jx*v)
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
SUBROUTINE depose_jxjyjz_scalar_2_2_2( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np
  LOGICAL(idp) :: l_nodal
  REAL(num)    :: stagger_shift
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
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, gaminv, w
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint
  REAL(num) :: xintsq, yintsq, zintsq
  REAL(num) :: x, y, z, xmid, ymid, zmid, vx, vy, vz, invvol, dts2dx, dts2dy, dts2dz
  REAL(num) :: wq, wqx, wqy, wqz, clightsq
  REAL(num), DIMENSION(3) :: sx(-1:1), sy(-1:1), sz(-1:1), sx0(-1:1), sy0(-1:1),      &
  sz0(-1:1)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, j0, k0, l0, ip

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
  sx=0.0_num;sy=0.0_num;sz=0.0_num;
  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;

  ! LOOP ON PARTICLES
  !$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)
  !$acc loop gang vector private(sx(-1:1), sy(-1:1), sz(-1:1), sx0(-1:1), sy0(-1:1), sz0(-1:1))
  DO ip=1, np
    ! --- computes position in  grid units at (n+1)
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
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
    ymid=y-dts2dy*vy
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j=nint(xmid)
    k=nint(ymid)
    l=nint(zmid)
    j0=nint(xmid-stagger_shift)
    k0=nint(ymid-stagger_shift)
    l0=nint(zmid-stagger_shift)
    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    yint = ymid-k
    zint = zmid-l
    xintsq = xint*xint
    sx(-1) = 0.5_num*(0.5_num-xint)**2
    sx( 0) = 0.75_num-xintsq
    sx( 1) = 0.5_num*(0.5_num+xint)**2
    yintsq = yint*yint
    sy(-1) = 0.5_num*(0.5_num-yint)**2
    sy( 0) = 0.75_num-yintsq
    sy( 1) = 0.5_num*(0.5_num+yint)**2
    zintsq = zint*zint
    sz(-1) = 0.5_num*(0.5_num-zint)**2
    sz( 0) = (0.75_num-zintsq)
    sz( 1) = 0.5_num*(0.5_num+zint)**2

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-stagger_shift-j0
    yint = ymid-stagger_shift-k0
    zint = zmid-stagger_shift-l0
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
    sz0( 0) = (0.75_num-zintsq)
    sz0( 1) = 0.5_num*(0.5_num+zint)**2

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
    ! --- to the 27 nearest vertices
    ! - JX
    !$acc atomic update
    jx(j0-1, k-1, l-1)  = jx(j0-1, k-1, l-1)  +   sx0(-1)*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k-1, l-1)  = jx(j0, k-1, l-1)  +   sx0(0 )*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k-1, l-1)  = jx(j0+1, k-1, l-1)  +   sx0(1 )*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k, l-1)  = jx(j0-1, k, l-1)  +   sx0(-1)*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k, l-1)  = jx(j0, k, l-1)  +   sx0(0 )*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k, l-1)  = jx(j0+1, k, l-1)  +   sx0(1 )*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k+1, l-1)  = jx(j0-1, k+1, l-1)  +   sx0(-1)*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k+1, l-1)  = jx(j0, k+1, l-1)  +   sx0(0 )*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k+1, l-1)  = jx(j0+1, k+1, l-1)  +   sx0(1 )*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k-1, l  )  = jx(j0-1, k-1, l  )  +   sx0(-1)*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k-1, l  )  = jx(j0, k-1, l  )  +   sx0(0 )*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k-1, l  )  = jx(j0+1, k-1, l  )  +   sx0(1 )*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k, l  )  = jx(j0-1, k, l  )  +   sx0(-1)*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k, l  )  = jx(j0, k, l  )  +   sx0(0 )*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k, l  )  = jx(j0+1, k, l  )  +   sx0(1 )*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k+1, l  )  = jx(j0-1, k+1, l  )  +   sx0(-1)*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k+1, l  )  = jx(j0, k+1, l  )  +   sx0(0 )*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k+1, l  )  = jx(j0+1, k+1, l  )  +   sx0(1 )*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k-1, l+1)  = jx(j0-1, k-1, l+1)  +   sx0(-1)*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k-1, l+1)  = jx(j0, k-1, l+1)  +   sx0(0 )*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k-1, l+1)  = jx(j0+1, k-1, l+1)  +   sx0(1 )*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, k, l+1)  = jx(j0-1, k, l+1)  +   sx0(-1)*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k, l+1)  = jx(j0, k, l+1)  +   sx0(0 )*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k, l+1)  = jx(j0+1, k, l+1)  +   sx0(1 )*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, k+1, l+1)  = jx(j0-1, k+1, l+1)  +   sx0(-1)*sy(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k+1, l+1)  = jx(j0, k+1, l+1)  +   sx0(0 )*sy(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k+1, l+1)  = jx(j0+1, k+1, l+1)  +   sx0(1 )*sy(1 )*sz(1 )*wqx

    !        ! - JY
    !$acc atomic update
    jy(j-1, k0-1, l-1)  = jy(j-1, k0-1, l-1)  +   sx(-1)*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0-1, l-1)  = jy(j, k0-1, l-1)  +   sx(0 )*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0-1, l-1)  = jy(j+1, k0-1, l-1)  +   sx(1 )*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0, l-1)  = jy(j-1, k0, l-1)  +   sx(-1)*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0, l-1)  = jy(j, k0, l-1)  +   sx(0 )*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0, l-1)  = jy(j+1, k0, l-1)  +   sx(1 )*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0+1, l-1)  = jy(j-1, k0+1, l-1)  +   sx(-1)*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0+1, l-1)  = jy(j, k0+1, l-1)  +   sx(0 )*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0+1, l-1)  = jy(j+1, k0+1, l-1)  +   sx(1 )*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0-1, l  )  = jy(j-1, k0-1, l  )  +   sx(-1)*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0-1, l  )  = jy(j, k0-1, l  )  +   sx(0 )*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0-1, l  )  = jy(j+1, k0-1, l  )  +   sx(1 )*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0, l  )  = jy(j-1, k0, l  )  +   sx(-1)*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0, l  )  = jy(j, k0, l  )  +   sx(0 )*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0, l  )  = jy(j+1, k0, l  )  +   sx(1 )*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0+1, l  )  = jy(j-1, k0+1, l  )  +   sx(-1)*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0+1, l  )  = jy(j, k0+1, l  )  +   sx(0 )*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0+1, l  )  = jy(j+1, k0+1, l  )  +   sx(1 )*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0-1, l+1)  = jy(j-1, k0-1, l+1)  +   sx(-1)*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0-1, l+1)  = jy(j, k0-1, l+1)  +   sx(0 )*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0-1, l+1)  = jy(j+1, k0-1, l+1)  +   sx(1 )*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, k0, l+1)  = jy(j-1, k0, l+1)  +   sx(-1)*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0, l+1)  = jy(j, k0, l+1)  +   sx(0 )*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0, l+1)  = jy(j+1, k0, l+1)  +   sx(1 )*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, k0+1, l+1)  = jy(j-1, k0+1, l+1)  +   sx(-1)*sy0(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0+1, l+1)  = jy(j, k0+1, l+1)  +   sx(0 )*sy0(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0+1, l+1)  = jy(j+1, k0+1, l+1)  +   sx(1 )*sy0(1 )*sz(1 )*wqy

    ! - JZ
    !$acc atomic update
    jz(j-1, k-1, l0-1)  = jz(j-1, k-1, l0-1)  +   sx(-1)*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k-1, l0-1)  = jz(j, k-1, l0-1)  +   sx(0 )*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k-1, l0-1)  = jz(j+1, k-1, l0-1)  +   sx(1 )*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k, l0-1)  = jz(j-1, k, l0-1)  +   sx(-1)*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k, l0-1)  = jz(j, k, l0-1)  +   sx(0 )*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k, l0-1)  = jz(j+1, k, l0-1)  +   sx(1 )*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k+1, l0-1)  = jz(j-1, k+1, l0-1)  +   sx(-1)*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k+1, l0-1)  = jz(j, k+1, l0-1)  +   sx(0 )*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k+1, l0-1)  = jz(j+1, k+1, l0-1)  +   sx(1 )*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k-1, l0  )  = jz(j-1, k-1, l0  )  +   sx(-1)*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k-1, l0  )  = jz(j, k-1, l0  )  +   sx(0 )*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k-1, l0  )  = jz(j+1, k-1, l0  )  +   sx(1 )*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k, l0  )  = jz(j-1, k, l0  )  +   sx(-1)*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k, l0  )  = jz(j, k, l0  )  +   sx(0 )*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k, l0  )  = jz(j+1, k, l0  )  +   sx(1 )*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k+1, l0  )  = jz(j-1, k+1, l0  )  +   sx(-1)*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k+1, l0  )  = jz(j, k+1, l0  )  +   sx(0 )*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k+1, l0  )  = jz(j+1, k+1, l0  )  +   sx(1 )*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k-1, l0+1)  = jz(j-1, k-1, l0+1)  +   sx(-1)*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k-1, l0+1)  = jz(j, k-1, l0+1)  +   sx(0 )*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k-1, l0+1)  = jz(j+1, k-1, l0+1)  +   sx(1 )*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, k, l0+1)  = jz(j-1, k, l0+1)  +   sx(-1)*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k, l0+1)  = jz(j, k, l0+1)  +   sx(0 )*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k, l0+1)  = jz(j+1, k, l0+1)  +   sx(1 )*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, k+1, l0+1)  = jz(j-1, k+1, l0+1)  +   sx(-1)*sy(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k+1, l0+1)  = jz(j, k+1, l0+1)  +   sx(0 )*sy(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k+1, l0+1)  = jz(j+1, k+1, l0+1)  +   sx(1 )*sy(1 )*sz0(1 )*wqz
  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_scalar_2_2_2


! ________________________________________________________________________________________
!> @brief
!> Order 2 3D vector direct current deposition routine (rho*v)
!
!> @details
!> This versions have good performances on SIMD architectures
!> Providing that OpenMP 4.0 is available (Directive SIMD)
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!
!> @param[in] np Number of particles
!> @param[in] xp 1D array of x-coordinates of particles
!> @param[in] yp 1D array of x-coordinates of particles
!> @param[in] zp 1D array of x-coordinates of particles
!> @param[in] uxp 1D array of ux-velocity components of particles
!> @param[in] uyp 1D array of ux-velocity components of particles
!> @param[in] uzp 1D array of ux-velocity components of particles
!> @param[in] gaminv 1D array of the inverse 1/gamma-factor of particles
!> @param[in] w 1D array of the weghts of particles
!> @param[in] q charge of current species (scalar)
!> @param[in] xmin x-minimum boundary of current tile
!> @param[in] ymin y-minimum boundary of current tile
!> @param[in] zmin z-minimum boundary of current tile
!> @param[in] dt time step (scalar)
!> @param[in] dx mesh size along x (scalar)
!> @param[in] dy mesh size along y (scalar)
!> @param[in] dz mesh size along z (scalar)
!> jx array (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!>  (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 3 integers)
!> @warning arrays jx, jy, jz should be set to 0 before entering this subroutine.
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHVv2_2_2_2(jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,  &
  xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE

  INTEGER(idp)             :: np
  LOGICAL(idp)             :: l_nodal
  REAL(num)                :: stagger_shift
  INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), &
  jy_nguard(3), jy_nvalid(3), &
  jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1, &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1, &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1, &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1, &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1, &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1, &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  REAL(num), DIMENSION(:, :), ALLOCATABLE :: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num)                :: xint, yint, zint
  REAL(num)                :: xintsq, yintsq, zintsq
  REAL(num)                :: x, y, z, xmid, ymid, zmid
  REAL(num)                :: wqx, wqy, wqz, wwx, wwy, wwz
  REAL(num)                :: invvol, dxi, dyi, dzi
  REAL(num)                :: dts2dx, dts2dy, dts2dz
  REAL(num), PARAMETER     :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp)             :: j, k, l, j0, k0, l0, ip, NCELLS, ic
  INTEGER(idp)             :: n, nn, nv
  INTEGER(idp), DIMENSION(LVEC, 3) :: ICELL, pos_jx, pos_jy, pos_jz
  REAL(num)                :: vx, vy, vz
  REAL(num)                :: ww0x(LVEC, 4), ww0y(LVEC, 4), ww0z(LVEC, 4), wwwx(LVEC, &
  8), wwwy(LVEC, 8), wwwz(LVEC, 8), wq
  REAL(num)                :: sx0(LVEC), sx1(LVEC), sx2(LVEC)
  REAL(num)                :: sx00(LVEC), sx01(LVEC), sx02(LVEC)
  REAL(num)                :: sy0, sy1, sy2, sy00, sy01, sy02
  REAL(num)                :: sz0, sz1, sz2, sz00, sz01, sz02, syz
  INTEGER(idp)             :: jorig, korig, lorig
  INTEGER(idp)             :: ncx, ncy, ncxy, ncz, ix, iy, iz

  ! Relative position of the 8 nodes (over which the algorithm vectorizes)
  ! in respect to the particle computed node (the nodes are only shifted in
  ! y and z, not in x ; see arXiv:1601.02056, Fig. 3, for order 2)
  INTEGER(idp), DIMENSION(8), PARAMETER  :: myoff=(/-1_idp, 0_idp, 1_idp, &
  -1_idp, 1_idp, -1_idp, 0_idp, 1_idp/)
  INTEGER(idp), DIMENSION(8), PARAMETER  :: mzoff=(/-1_idp, -1_idp, -1_idp, &
  0_idp, 0_idp, 1_idp, 1_idp, 1_idp/)

  ! ___ Parameter initialization _________________

  ww0x=0._num; ww0y=0._num; ww0z=0._num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi

  ! Lower integer position where the particles can deposit
  jorig=-2; korig=-2;lorig=-2
  ! Find the maximal number of cells in each direction, so as to
  ! allocate the linearized 1D arrays
  ncx = MAX( jx_nvalid(1), jy_nvalid(1), jz_nvalid(1) ) + 4
  ncy = MAX( jx_nvalid(2), jy_nvalid(2), jz_nvalid(2) ) + 4
  ncz = MAX( jx_nvalid(3), jy_nvalid(3), jz_nvalid(3) ) + 4
  ncxy = ncx*ncy
  NCELLS = ncx*ncy*ncz
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
  jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num

  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED gaminv:64
    !DIR$ ASSUME_ALIGNED sx0:64, sx1:64, sx2:64
    !DIR$ ASSUME_ALIGNED sx00:64, sx01:64, sx02:64
    !DIR$ ASSUME_ALIGNED w:64, wwwx:64, wwwy:64, wwwz:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, sx0, sx1, sx2)
    !IBM* ALIGN(64, sx00, sx01, sx02)
    !IBM* ALIGN(64, w, wwwx, wwwy, wwwz)
    !IBM* ALIGN(64, ICELL)
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
    DO n=1, MIN(LVEC, np-ip+1)
      nn=ip+n-1
      ! --- computes position in  grid units at (n+1)
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! Computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)

      ! --- computes particles weights
      wq=q*w(nn)*invvol
      wqx=wq*vx
      wqy=wq*vy
      wqz=wq*vz

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx
      ymid=y-dts2dy*vy
      zmid=z-dts2dz*vz

      ! --- finds node of cell containing particles for current positions
      j=nint(xmid)
      k=nint(ymid)
      l=nint(zmid)
      j0=nint(xmid-stagger_shift)
      k0=nint(ymid-stagger_shift)
      l0=nint(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! Store integer positions of the LVEC particles with respect to
      ! the staggered grids of jx, jy and jz
      pos_jx(n,:) = (/ j0, k, l /)
      pos_jy(n,:) = (/ j, k0, l /)
      pos_jz(n,:) = (/ j, k, l0 /)

      ! --- computes set of coefficients for node centered quantities
      xint = xmid-j
      yint = ymid-k
      zint = zmid-l
      xintsq= xint**2
      yintsq= yint**2
      zintsq= zint**2
      sx0(n)=0.5_num*(0.5_num-xint)**2
      sx1(n)=(0.75_num-xintsq)
      sx2(n)=0.5_num*(0.5_num+xint)**2
      sy0=0.5_num*(0.5_num-yint)**2
      sy1=(0.75_num-yintsq)
      sy2=0.5_num*(0.5_num+yint)**2
      sz0=0.5_num*(0.5_num-zint)**2
      sz1=(0.75_num-zintsq)
      sz2=0.5_num*(0.5_num+zint)**2

      ! --- computes set of coefficients for staggered quantities
      xint = xmid-stagger_shift-j0
      yint = ymid-stagger_shift-k0
      zint = zmid-stagger_shift-l0
      xintsq= xint**2
      yintsq= yint**2
      zintsq= zint**2
      sx00(n)=0.5_num*(0.5_num-xint)**2
      sx01(n)=(0.75_num-xintsq)
      sx02(n)=0.5_num*(0.5_num+xint)**2
      sy00=0.5_num*(0.5_num-yint)**2
      sy01=(0.75_num-yintsq)
      sy02=0.5_num*(0.5_num+yint)**2
      sz00=0.5_num*(0.5_num-zint)**2
      sz01=(0.75_num-zintsq)
      sz02=0.5_num*(0.5_num+zint)**2

      ! -- Weights for planes of 8  vertices
      ! Weights - X
      wwwx(n, 1) = sy0*sz0*wqx
      wwwx(n, 2) = sy1*sz0*wqx
      wwwx(n, 3) = sy2*sz0*wqx
      wwwx(n, 4) = sy0*sz1*wqx
      wwwx(n, 5) = sy2*sz1*wqx
      wwwx(n, 6) = sy0*sz2*wqx
      wwwx(n, 7) = sy1*sz2*wqx
      wwwx(n, 8) = sy2*sz2*wqx

      ! Weights - Y
      wwwy(n, 1) = sy00*sz0*wqy
      wwwy(n, 2) = sy01*sz0*wqy
      wwwy(n, 3) = sy02*sz0*wqy
      wwwy(n, 4) = sy00*sz1*wqy
      wwwy(n, 5) = sy02*sz1*wqy
      wwwy(n, 6) = sy00*sz2*wqy
      wwwy(n, 7) = sy01*sz2*wqy
      wwwy(n, 8) = sy02*sz2*wqy

      ! Weights - Z
      wwwz(n, 1) = sy0*sz00*wqz
      wwwz(n, 2) = sy1*sz00*wqz
      wwwz(n, 3) = sy2*sz00*wqz
      wwwz(n, 4) = sy0*sz01*wqz
      wwwz(n, 5) = sy2*sz01*wqz
      wwwz(n, 6) = sy0*sz02*wqz
      wwwz(n, 7) = sy1*sz02*wqz
      wwwz(n, 8) = sy2*sz02*wqz

      ! -- 3 remaining central points
      syz=sz1*sy1*wqx
      ww0x(n, 1)=syz*sx00(n)
      ww0x(n, 2)=syz*sx01(n)
      ww0x(n, 3)=syz*sx02(n)
      syz=sz1*sy01*wqy
      ww0y(n, 1)=syz*sx0(n)
      ww0y(n, 2)=syz*sx1(n)
      ww0y(n, 3)=syz*sx2(n)
      syz=sz01*sy1*wqz
      ww0z(n, 1)=syz*sx0(n)
      ww0z(n, 2)=syz*sx1(n)
      ww0z(n, 3)=syz*sx2(n)

    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! LOOP OVER LVEC CONSECUTIVE PARTICLES
    DO n=1, MIN(LVEC, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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
      !DIR$ SIMD
#endif
      DO nv=1, 8
        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! - JX
        wwx=wwwx(n, nv)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)-1) = jxcells(nv, ICELL(n, 1)-1) +wwx*sx00(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1))   = jxcells(nv, ICELL(n, 1))   +wwx*sx01(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+1) = jxcells(nv, ICELL(n, 1)+1) +wwx*sx02(n)
        ! - JY
        wwy=wwwy(n, nv)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)-1) = jycells(nv, ICELL(n, 2)-1) +wwy*sx0(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2))   = jycells(nv, ICELL(n, 2))   +wwy*sx1(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+1) = jycells(nv, ICELL(n, 2)+1) +wwy*sx2(n)
        ! - JZ
        wwz=wwwz(n, nv)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)-1) = jzcells(nv, ICELL(n, 3)-1) +wwz*sx0(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3))   = jzcells(nv, ICELL(n, 3))   +wwz*sx1(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+1) = jzcells(nv, ICELL(n, 3)+1) +wwz*sx2(n)

      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif

    END DO
    DO n=1, MIN(LVEC, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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
      ! For each of the LVEC particles, add the contribution for the 3
      ! central points which are not vectorized (see arXiv:1601.02056,
      ! Fig. 3, for order 2). These 3 points are 3 consecutive gridpoints in x
      ! (and correspond to the same position in y and z)
      DO nv=1, 4
        jx( pos_jx(n,1)+nv-2, pos_jx(n,2), pos_jx(n,3) ) = &
        jx( pos_jx(n,1)+nv-2, pos_jx(n,2), pos_jx(n,3) ) + ww0x(n, nv)
        jy( pos_jy(n,1)+nv-2, pos_jy(n,2), pos_jy(n,3) ) = &
        jy( pos_jy(n,1)+nv-2, pos_jy(n,2), pos_jy(n,3) ) + ww0y(n, nv)
        jz( pos_jz(n,1)+nv-2, pos_jz(n,2), pos_jz(n,3) ) = &
        jz( pos_jz(n,1)+nv-2, pos_jz(n,2), pos_jz(n,3) ) + ww0z(n, nv)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! ---------------------------------------------------------------
  ! Vectorized reduction of jxcells, jycells, jzcells in jx, jy, jz
  ! ---------------------------------------------------------------

  ! Reduction of jxcells in jx
  ! Note: the bounds below make sure that we never go out-of-bound for jxcells,
  ! even when taking into account the additional +/- 1 offsets in myoff and mzoff
  DO iz = MAX(lorig, -jx_nguard(3)+1), MIN(lorig+ncz-1, jx_nvalid(3)+jx_nguard(3)-2)
    DO iy = MAX(korig, -jx_nguard(2)+1), MIN(korig+ncy-1, jx_nvalid(2)+jx_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
      DO ix = MAX(jorig, -jx_nguard(1)), MIN(jorig+ncx-1, jx_nvalid(1)+jx_nguard(1)-2)
        ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
        ! Compute linearized index
        jx(ix,iy+myoff(1),iz+mzoff(1))=jx(ix,iy+myoff(1),iz+mzoff(1))+jxcells(1, ic)
        jx(ix,iy+myoff(2),iz+mzoff(2))=jx(ix,iy+myoff(2),iz+mzoff(2))+jxcells(2, ic)
        jx(ix,iy+myoff(3),iz+mzoff(3))=jx(ix,iy+myoff(3),iz+mzoff(3))+jxcells(3, ic)
        jx(ix,iy+myoff(4),iz+mzoff(4))=jx(ix,iy+myoff(4),iz+mzoff(4))+jxcells(4, ic)
        jx(ix,iy+myoff(5),iz+mzoff(5))=jx(ix,iy+myoff(5),iz+mzoff(5))+jxcells(5, ic)
        jx(ix,iy+myoff(6),iz+mzoff(6))=jx(ix,iy+myoff(6),iz+mzoff(6))+jxcells(6, ic)
        jx(ix,iy+myoff(7),iz+mzoff(7))=jx(ix,iy+myoff(7),iz+mzoff(7))+jxcells(7, ic)
        jx(ix,iy+myoff(8),iz+mzoff(8))=jx(ix,iy+myoff(8),iz+mzoff(8))+jxcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jycells in jy
  ! Note: the bounds below make sure that we never go out-of-bound for jycells,
  ! even when taking into account the additional +/- 1 offsets in myoff and mzoff
  DO iz = MAX(lorig, -jy_nguard(3)+1), MIN(lorig+ncz-1, jy_nvalid(3)+jy_nguard(3)-2)
    DO iy = MAX(korig, -jy_nguard(2)+1), MIN(korig+ncy-1, jy_nvalid(2)+jy_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
    DO ix = MAX(jorig, -jy_nguard(1)), MIN(jorig+ncx-1, jy_nvalid(1)+jy_nguard(1)-2)
        ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
        ! Compute linearized index
        jy(ix,iy+myoff(1),iz+mzoff(1))=jy(ix,iy+myoff(1),iz+mzoff(1))+jycells(1, ic)
        jy(ix,iy+myoff(2),iz+mzoff(2))=jy(ix,iy+myoff(2),iz+mzoff(2))+jycells(2, ic)
        jy(ix,iy+myoff(3),iz+mzoff(3))=jy(ix,iy+myoff(3),iz+mzoff(3))+jycells(3, ic)
        jy(ix,iy+myoff(4),iz+mzoff(4))=jy(ix,iy+myoff(4),iz+mzoff(4))+jycells(4, ic)
        jy(ix,iy+myoff(5),iz+mzoff(5))=jy(ix,iy+myoff(5),iz+mzoff(5))+jycells(5, ic)
        jy(ix,iy+myoff(6),iz+mzoff(6))=jy(ix,iy+myoff(6),iz+mzoff(6))+jycells(6, ic)
        jy(ix,iy+myoff(7),iz+mzoff(7))=jy(ix,iy+myoff(7),iz+mzoff(7))+jycells(7, ic)
        jy(ix,iy+myoff(8),iz+mzoff(8))=jy(ix,iy+myoff(8),iz+mzoff(8))+jycells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jzcells in jz
  ! Note: the bounds below make sure that we never go out-of-bound for jzcells,
  ! even when taking into account the additional +/- 1 offsets in myoff and mzoff
  DO iz = MAX(lorig, -jz_nguard(3)+1), MIN(lorig+ncz-1, jz_nvalid(3)+jz_nguard(3)-2)
    DO iy = MAX(korig, -jz_nguard(2)+1), MIN(korig+ncy-1, jz_nvalid(2)+jz_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
    DO ix = MAX(jorig, -jz_nguard(1)), MIN(jorig+ncx-1, jz_nvalid(1)+jz_nguard(1)-2)
        ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
        ! Compute linearized index
        jz(ix,iy+myoff(1),iz+mzoff(1))=jz(ix,iy+myoff(1),iz+mzoff(1))+jzcells(1, ic)
        jz(ix,iy+myoff(2),iz+mzoff(2))=jz(ix,iy+myoff(2),iz+mzoff(2))+jzcells(2, ic)
        jz(ix,iy+myoff(3),iz+mzoff(3))=jz(ix,iy+myoff(3),iz+mzoff(3))+jzcells(3, ic)
        jz(ix,iy+myoff(4),iz+mzoff(4))=jz(ix,iy+myoff(4),iz+mzoff(4))+jzcells(4, ic)
        jz(ix,iy+myoff(5),iz+mzoff(5))=jz(ix,iy+myoff(5),iz+mzoff(5))+jzcells(5, ic)
        jz(ix,iy+myoff(6),iz+mzoff(6))=jz(ix,iy+myoff(6),iz+mzoff(6))+jzcells(6, ic)
        jz(ix,iy+myoff(7),iz+mzoff(7))=jz(ix,iy+myoff(7),iz+mzoff(7))+jzcells(7, ic)
        jz(ix,iy+myoff(8),iz+mzoff(8))=jz(ix,iy+myoff(8),iz+mzoff(8))+jzcells(8, ic)
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
END SUBROUTINE depose_jxjyjz_vecHVv2_2_2_2


! ________________________________________________________________________________________
! depose_jxjyjz_vecHV_vnr_2_2_2
!
!> @brief
!> Order 2 3D vector current deposition routine (rho*v) with no reduction
!
!> @details
!> This versions have good performances on SIMD architectures
!> providing that OpenMP 4.0 is available (Directive SIMD).
!> This subroutine is similar to depose_jxjyjz_vecHVv2_1_1_1
!> without the reduction process at the end.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2016
!>
!> @param[inout] jxcells, jycells, jzcells transient current arrays
!> @param[in] np particle number
!> @param[in] ncells number of cells in the tile
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse Lorentz factor arrays
!> @param[in] w particle wight arrays
!> @param[in] q charge
!> @param[in] xmin, ymin, zmin tile minimum positions
!> @param[in] dt, dx, dy, dz time and space steps
!> @param[in] nx, ny, nz tile cell numbers in each direction
!> @param[in] nxguard, nyguard, nzguard guard cells
!> @param[in] ncx, ncy, ncz tile cell extended number (depends on the order)
!> @param[in] lvect vector length
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2(jxcells, jycells, jzcells, np, ncells, xp,   &
  yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz,    &
  nxguard, nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not wrap
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE
  ! ____ Parameter initialization _____________________________________
  INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
 LOGICAL(idp)                      :: l_nodal
 REAL(num)                         :: stagger_shift
  INTEGER(idp), INTENT(IN)                      :: ncx, ncy, ncz
  INTEGER(idp)                                  :: nxguard, nyguard, nzguard
  INTEGER(idp), INTENT(IN)                      :: lvect
  REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np)                      :: xp, yp, zp, uxp, uyp, uzp, gaminv, &
  w
  REAL(num)                                     :: q, dt, dx, dy, dz, xmin, ymin,     &
  zmin
  REAL(num)                                     :: xint, yint, zint
  REAL(num)                                     :: xintsq, yintsq, zintsq
  REAL(num)                                     :: x, y, z, xmid, ymid, zmid
  REAL(num)                                     ::   wqx, wqy, wqz, wwx, wwy, wwz
  REAL(num)                                     :: invvol, dxi, dyi, dzi
  REAL(num)                                     :: dts2dx, dts2dy, dts2dz
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(isp)                                  :: j, k, l, j0, k0, l0, ip
  INTEGER(isp)                                  :: nnx, nnxy, n, nn, nv
  INTEGER(isp), DIMENSION(LVECT, 3)              :: ICELL, IG
  REAL(num)                                     :: vx, vy, vz
  REAL(num)                                     :: ww0x(LVECT, 4), ww0y(LVECT, 4),    &
  ww0z(LVECT, 4)
  REAL(num)                                     :: wwwx(LVECT, 8), wwwy(LVECT, 8),    &
  wwwz(LVECT, 8)
  REAL(num)                                     :: wq
  REAL(num)                                     :: sx0(LVECT), sx1(LVECT), sx2(LVECT)
  REAL(num)                                     :: sx00(LVECT), sx01(LVECT),          &
  sx02(LVECT)
  REAL(num)                                     :: sy0, sy1, sy2, sy00, sy01, sy02
  REAL(num)                                     :: sz0, sz1, sz2, sz00, sz01, sz02,   &
  syz
  INTEGER(isp)                                  :: orig, jorig, korig, lorig
  INTEGER(isp)                                  :: ncxy
  INTEGER(isp)                                  :: ngridx, ngridy, ngx, ngxy

  ! __________________________________________________________
  ! Parameters

  ww0x=0._num
  ww0y=0._num
  ww0z=0._num

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi

  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard

  nnx = nx + 1 + 2*nxguard
  nnxy = nnx*(ny+1+2*nyguard)

  jorig=-2
  korig=-2
  lorig=-2
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! ___ Computation ______________________________________________

  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED sx0:64, sx1:64, sx2:64
    !DIR$ ASSUME_ALIGNED sx00:64, sx01:64, sx02:64
    !DIR$ ASSUME_ALIGNED w:64, wwwx:64, wwwy:64, wwwz:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, sx0, sx1, sx2)
    !IBM* ALIGN(64, sx00, sx01, sx02)
    !IBM* ALIGN(64, w, wwwx, wwwy, wwwz)
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

      ! Computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)

      ! --- computes particles weights
      wq=q*w(nn)*invvol
      wqx=wq*vx
      wqy=wq*vy
      wqz=wq*vz

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx
      ymid=y-dts2dy*vy
      zmid=z-dts2dz*vz

      ! --- finds node of cell containing particles for current positions
      j=nint(xmid)
      k=nint(ymid)
      l=nint(zmid)
      j0=nint(xmid-stagger_shift)
      k0=nint(ymid-stagger_shift)
      l0=nint(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy
      IG(n, 1) = ICELL(n, 1) + ncx + ncxy
      IG(n, 2) = ICELL(n, 2) + ncx + ncxy
      IG(n, 3) = ICELL(n, 3) + ncx + ncxy

      ! --- computes set of coefficients for node centered quantities
      xint = xmid-j
      yint = ymid-k
      zint = zmid-l
      xintsq= xint**2
      yintsq= yint**2
      zintsq= zint**2
      sx0(n)=0.5_num*(0.5_num-xint)**2
      sx1(n)=(0.75_num-xintsq)
      sx2(n)=0.5_num*(0.5_num+xint)**2
      sy0=0.5_num*(0.5_num-yint)**2
      sy1=(0.75_num-yintsq)
      sy2=0.5_num*(0.5_num+yint)**2
      sz0=0.5_num*(0.5_num-zint)**2
      sz1=(0.75_num-zintsq)
      sz2=0.5_num*(0.5_num+zint)**2

      ! --- computes set of coefficients for staggered quantities
      xint = xmid-stagger_shift-j0
      yint = ymid-stagger_shift-k0
      zint = zmid-stagger_shift-l0
      xintsq= xint**2
      yintsq= yint**2
      zintsq= zint**2
      sx00(n)=0.5_num*(0.5_num-xint)**2
      sx01(n)=(0.75_num-xintsq)
      sx02(n)=0.5_num*(0.5_num+xint)**2
      sy00=0.5_num*(0.5_num-yint)**2
      sy01=(0.75_num-yintsq)
      sy02=0.5_num*(0.5_num+yint)**2
      sz00=0.5_num*(0.5_num-zint)**2
      sz01=(0.75_num-zintsq)
      sz02=0.5_num*(0.5_num+zint)**2

      ! -- Weights for planes of 8  vertices
      ! Weights - X
      wwwx(n, 1) = sy0*sz0*wqx
      wwwx(n, 2) = sy1*sz0*wqx
      wwwx(n, 3) = sy2*sz0*wqx
      wwwx(n, 4) = sy0*sz1*wqx
      wwwx(n, 5) = sy2*sz1*wqx
      wwwx(n, 6) = sy0*sz2*wqx
      wwwx(n, 7) = sy1*sz2*wqx
      wwwx(n, 8) = sy2*sz2*wqx

      ! Weights - Y
      wwwy(n, 1) = sy00*sz0*wqy
      wwwy(n, 2) = sy01*sz0*wqy
      wwwy(n, 3) = sy02*sz0*wqy
      wwwy(n, 4) = sy00*sz1*wqy
      wwwy(n, 5) = sy02*sz1*wqy
      wwwy(n, 6) = sy00*sz2*wqy
      wwwy(n, 7) = sy01*sz2*wqy
      wwwy(n, 8) = sy02*sz2*wqy

      ! Weights - Z
      wwwz(n, 1) = sy0*sz00*wqz
      wwwz(n, 2) = sy1*sz00*wqz
      wwwz(n, 3) = sy2*sz00*wqz
      wwwz(n, 4) = sy0*sz01*wqz
      wwwz(n, 5) = sy2*sz01*wqz
      wwwz(n, 6) = sy0*sz02*wqz
      wwwz(n, 7) = sy1*sz02*wqz
      wwwz(n, 8) = sy2*sz02*wqz

      ! -- 3 remaining central points
      syz=sz1*sy1*wqx
      ww0x(n, 1)=syz*sx00(n)
      ww0x(n, 2)=syz*sx01(n)
      ww0x(n, 3)=syz*sx02(n)
      syz=sz1*sy01*wqy
      ww0y(n, 1)=syz*sx0(n)
      ww0y(n, 2)=syz*sx1(n)
      ww0y(n, 3)=syz*sx2(n)
      syz=sz01*sy1*wqz
      ww0z(n, 1)=syz*sx0(n)
      ww0z(n, 2)=syz*sx1(n)
      ww0z(n, 3)=syz*sx2(n)
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    DO n=1, MIN(LVEC, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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
        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! - JX
        wwx=wwwx(n, nv)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)-1) = jxcells(nv, ICELL(n, 1)-1) +wwx*sx00(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1))   = jxcells(nv, ICELL(n, 1))   +wwx*sx01(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+1) = jxcells(nv, ICELL(n, 1)+1) +wwx*sx02(n)
        ! - JY
        wwy=wwwy(n, nv)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)-1) = jycells(nv, ICELL(n, 2)-1) +wwy*sx0(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2))   = jycells(nv, ICELL(n, 2))   +wwy*sx1(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+1) = jycells(nv, ICELL(n, 2)+1) +wwy*sx2(n)
        ! - JZ
        wwz=wwwz(n, nv)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)-1) = jzcells(nv, ICELL(n, 3)-1) +wwz*sx0(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3))   = jzcells(nv, ICELL(n, 3))   +wwz*sx1(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+1) = jzcells(nv, ICELL(n, 3)+1) +wwz*sx2(n)
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
        !IF (IG(n, 1)+nv-2 .gt.ncells) THEN
        !  print*, nv, n
        !endif
        jxcells(1, IG(n, 1)+nv-2) = jxcells(1, IG(n, 1)+nv-2) + ww0x(n, nv)
        jycells(1, IG(n, 2)+nv-2) = jycells(1, IG(n, 2)+nv-2) + ww0y(n, nv)
        jzcells(1, IG(n, 3)+nv-2) = jzcells(1, IG(n, 3)+nv-2) + ww0z(n, nv)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  RETURN
END SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2


! ______________________________________________________________________________
!> @brief
!> Order 3 3D scalar current deposition routine (rho*v)
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
SUBROUTINE depose_jxjyjz_scalar_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,       &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np
  LOGICAL(idp) :: l_nodal
  REAL(num)    :: stagger_shift
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
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq,  &
  zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, xmid, ymid, zmid, vx, vy, vz, invvol, dts2dx, dts2dy, dts2dz
  REAL(num) :: wq, wqx, wqy, wqz, clightsq
  REAL(num), DIMENSION(4) :: sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2),      &
  sz0(-1:2)
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, j0, k0, l0, ip


  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
  sx=0.0_num;sy=0.0_num;sz=0.0_num;
  sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;

  ! LOOP ON PARTICLES
  !$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)
  !$acc loop gang vector private(sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2), sz0(-1:2))
  DO ip=1, np
    ! --- computes position in  grid units at (n+1)
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
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
    ymid=y-dts2dy*vy
    zmid=z-dts2dz*vz

    ! --- finds node of cell containing particles for current positions
    j=floor(xmid)
    k=floor(ymid)
    l=floor(zmid)
    j0=floor(xmid-stagger_shift)
    k0=floor(ymid-stagger_shift)
    l0=floor(zmid-stagger_shift)

    ! --- computes set of coefficients for node centered quantities
    xint = xmid-j
    yint = ymid-k
    zint = zmid-l
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1) = onesixth*oxintsq*oxint
    sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx( 2) = onesixth*xintsq*xint
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
    sz(-1) = onesixth*ozintsq*ozint
    sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz( 2) = onesixth*zintsq*zint

    ! --- computes set of coefficients for staggered quantities
    xint = xmid-stagger_shift-j0
    yint = ymid-stagger_shift-k0
    zint = zmid-stagger_shift-l0
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx0(-1) = onesixth*oxintsq*oxint
    sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx0( 2) = onesixth*xintsq*xint
    oyint = 1.0_num-yint
    yintsq = yint*yint
    oyintsq = oyint*oyint
    sy0(-1) = onesixth*oyintsq*oyint
    sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
    sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
    sy0( 2) = onesixth*yintsq*yint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz0(-1) = onesixth*ozintsq*ozint
    sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz0( 2) = onesixth*zintsq*zint

    ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
    ! --- to the 64 nearest vertices
    ! - JX
    !$acc atomic update
    jx(j0-1, k-1, l-1)  = jx(j0-1, k-1, l-1)  +   sx0(-1)*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k-1, l-1)  = jx(j0, k-1, l-1)  +   sx0(0 )*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k-1, l-1)  = jx(j0+1, k-1, l-1)  +   sx0(1 )*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0+2, k-1, l-1)  = jx(j0+2, k-1, l-1)  +   sx0(2 )*sy(-1)*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k, l-1)  = jx(j0-1, k, l-1)  +   sx0(-1)*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k, l-1)  = jx(j0, k, l-1)  +   sx0(0 )*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k, l-1)  = jx(j0+1, k, l-1)  +   sx0(1 )*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+2, k, l-1)  = jx(j0+2, k, l-1)  +   sx0(2 )*sy(0 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k+1, l-1)  = jx(j0-1, k+1, l-1)  +   sx0(-1)*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k+1, l-1)  = jx(j0, k+1, l-1)  +   sx0(0 )*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k+1, l-1)  = jx(j0+1, k+1, l-1)  +   sx0(1 )*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+2, k+1, l-1)  = jx(j0+2, k+1, l-1)  +   sx0(2 )*sy(1 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k+2, l-1)  = jx(j0-1, k+2, l-1)  +   sx0(-1)*sy(2 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0, k+2, l-1)  = jx(j0, k+2, l-1)  +   sx0(0 )*sy(2 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+1, k+2, l-1)  = jx(j0+1, k+2, l-1)  +   sx0(1 )*sy(2 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0+2, k+2, l-1)  = jx(j0+2, k+2, l-1)  +   sx0(2 )*sy(2 )*sz(-1)*wqx
    !$acc atomic update
    jx(j0-1, k-1, l  )  = jx(j0-1, k-1, l  )  +   sx0(-1)*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k-1, l  )  = jx(j0, k-1, l  )  +   sx0(0 )*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k-1, l  )  = jx(j0+1, k-1, l  )  +   sx0(1 )*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0+2, k-1, l  )  = jx(j0+2, k-1, l  )  +   sx0(2 )*sy(-1)*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k, l  )  = jx(j0-1, k, l  )  +   sx0(-1)*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k, l  )  = jx(j0, k, l  )  +   sx0(0 )*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k, l  )  = jx(j0+1, k, l  )  +   sx0(1 )*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+2, k, l  )  = jx(j0+2, k, l  )  +   sx0(2 )*sy(0 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k+1, l  )  = jx(j0-1, k+1, l  )  +   sx0(-1)*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k+1, l  )  = jx(j0, k+1, l  )  +   sx0(0 )*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k+1, l  )  = jx(j0+1, k+1, l  )  +   sx0(1 )*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+2, k+1, l  )  = jx(j0+2, k+1, l  )  +   sx0(2 )*sy(1 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k+2, l  )  = jx(j0-1, k+2, l  )  +   sx0(-1)*sy(2 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0, k+2, l  )  = jx(j0, k+2, l  )  +   sx0(0 )*sy(2 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+1, k+2, l  )  = jx(j0+1, k+2, l  )  +   sx0(1 )*sy(2 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0+2, k+2, l  )  = jx(j0+2, k+2, l  )  +   sx0(2 )*sy(2 )*sz(0 )*wqx
    !$acc atomic update
    jx(j0-1, k-1, l+1)  = jx(j0-1, k-1, l+1)  +   sx0(-1)*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k-1, l+1)  = jx(j0, k-1, l+1)  +   sx0(0 )*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k-1, l+1)  = jx(j0+1, k-1, l+1)  +   sx0(1 )*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0+2, k-1, l+1)  = jx(j0+2, k-1, l+1)  +   sx0(2 )*sy(-1)*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, k, l+1)  = jx(j0-1, k, l+1)  +   sx0(-1)*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k, l+1)  = jx(j0, k, l+1)  +   sx0(0 )*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k, l+1)  = jx(j0+1, k, l+1)  +   sx0(1 )*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+2, k, l+1)  = jx(j0+2, k, l+1)  +   sx0(2 )*sy(0 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, k+1, l+1)  = jx(j0-1, k+1, l+1)  +   sx0(-1)*sy(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k+1, l+1)  = jx(j0, k+1, l+1)  +   sx0(0 )*sy(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k+1, l+1)  = jx(j0+1, k+1, l+1)  +   sx0(1 )*sy(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+2, k+1, l+1)  = jx(j0+2, k+1, l+1)  +   sx0(2 )*sy(1 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, k+2, l+1)  = jx(j0-1, k+2, l+1)  +   sx0(-1)*sy(2 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0, k+2, l+1)  = jx(j0, k+2, l+1)  +   sx0(0 )*sy(2 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+1, k+2, l+1)  = jx(j0+1, k+2, l+1)  +   sx0(1 )*sy(2 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0+2, k+2, l+1)  = jx(j0+2, k+2, l+1)  +   sx0(2 )*sy(2 )*sz(1 )*wqx
    !$acc atomic update
    jx(j0-1, k-1, l+2)  = jx(j0-1, k-1, l+2)  +   sx0(-1)*sy(-1)*sz(2 )*wqx
    !$acc atomic update
    jx(j0, k-1, l+2)  = jx(j0, k-1, l+2)  +   sx0(0 )*sy(-1)*sz(2 )*wqx
    !$acc atomic update
    jx(j0+1, k-1, l+2)  = jx(j0+1, k-1, l+2)  +   sx0(1 )*sy(-1)*sz(2 )*wqx
    !$acc atomic update
    jx(j0+2, k-1, l+2)  = jx(j0+2, k-1, l+2)  +   sx0(2 )*sy(-1)*sz(2 )*wqx
    !$acc atomic update
    jx(j0-1, k, l+2)  = jx(j0-1, k, l+2)  +   sx0(-1)*sy(0 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0, k, l+2)  = jx(j0, k, l+2)  +   sx0(0 )*sy(0 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+1, k, l+2)  = jx(j0+1, k, l+2)  +   sx0(1 )*sy(0 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+2, k, l+2)  = jx(j0+2, k, l+2)  +   sx0(2 )*sy(0 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0-1, k+1, l+2)  = jx(j0-1, k+1, l+2)  +   sx0(-1)*sy(1 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0, k+1, l+2)  = jx(j0, k+1, l+2)  +   sx0(0 )*sy(1 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+1, k+1, l+2)  = jx(j0+1, k+1, l+2)  +   sx0(1 )*sy(1 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+2, k+1, l+2)  = jx(j0+2, k+1, l+2)  +   sx0(2 )*sy(1 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0-1, k+2, l+2)  = jx(j0-1, k+2, l+2)  +   sx0(-1)*sy(2 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0, k+2, l+2)  = jx(j0, k+2, l+2)  +   sx0(0 )*sy(2 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+1, k+2, l+2)  = jx(j0+1, k+2, l+2)  +   sx0(1 )*sy(2 )*sz(2 )*wqx
    !$acc atomic update
    jx(j0+2, k+2, l+2)  = jx(j0+2, k+2, l+2)  +   sx0(2 )*sy(2 )*sz(2 )*wqx

    ! - JY
    !$acc atomic update
    jy(j-1, k0-1, l-1)  = jy(j-1, k0-1, l-1)  +   sx(-1)*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0-1, l-1)  = jy(j, k0-1, l-1)  +   sx(0 )*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0-1, l-1)  = jy(j+1, k0-1, l-1)  +   sx(1 )*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j+2, k0-1, l-1)  = jy(j+2, k0-1, l-1)  +   sx(2 )*sy0(-1)*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0, l-1)  = jy(j-1, k0, l-1)  +   sx(-1)*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0, l-1)  = jy(j, k0, l-1)  +   sx(0 )*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0, l-1)  = jy(j+1, k0, l-1)  +   sx(1 )*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+2, k0, l-1)  = jy(j+2, k0, l-1)  +   sx(2 )*sy0(0 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0+1, l-1)  = jy(j-1, k0+1, l-1)  +   sx(-1)*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0+1, l-1)  = jy(j, k0+1, l-1)  +   sx(0 )*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0+1, l-1)  = jy(j+1, k0+1, l-1)  +   sx(1 )*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+2, k0+1, l-1)  = jy(j+2, k0+1, l-1)  +   sx(2 )*sy0(1 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0+2, l-1)  = jy(j-1, k0+2, l-1)  +   sx(-1)*sy0(2 )*sz(-1)*wqy
    !$acc atomic update
    jy(j, k0+2, l-1)  = jy(j, k0+2, l-1)  +   sx(0 )*sy0(2 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+1, k0+2, l-1)  = jy(j+1, k0+2, l-1)  +   sx(1 )*sy0(2 )*sz(-1)*wqy
    !$acc atomic update
    jy(j+2, k0+2, l-1)  = jy(j+2, k0+2, l-1)  +   sx(2 )*sy0(2 )*sz(-1)*wqy
    !$acc atomic update
    jy(j-1, k0-1, l  )  = jy(j-1, k0-1, l  )  +   sx(-1)*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0-1, l  )  = jy(j, k0-1, l  )  +   sx(0 )*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0-1, l  )  = jy(j+1, k0-1, l  )  +   sx(1 )*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j+2, k0-1, l  )  = jy(j+2, k0-1, l  )  +   sx(2 )*sy0(-1)*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0, l  )  = jy(j-1, k0, l  )  +   sx(-1)*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0, l  )  = jy(j, k0, l  )  +   sx(0 )*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0, l  )  = jy(j+1, k0, l  )  +   sx(1 )*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+2, k0, l  )  = jy(j+2, k0, l  )  +   sx(2 )*sy0(0 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0+1, l  )  = jy(j-1, k0+1, l  )  +   sx(-1)*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0+1, l  )  = jy(j, k0+1, l  )  +   sx(0 )*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0+1, l  )  = jy(j+1, k0+1, l  )  +   sx(1 )*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+2, k0+1, l  )  = jy(j+2, k0+1, l  )  +   sx(2 )*sy0(1 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0+2, l  )  = jy(j-1, k0+2, l  )  +   sx(-1)*sy0(2 )*sz(0 )*wqy
    !$acc atomic update
    jy(j, k0+2, l  )  = jy(j, k0+2, l  )  +   sx(0 )*sy0(2 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+1, k0+2, l  )  = jy(j+1, k0+2, l  )  +   sx(1 )*sy0(2 )*sz(0 )*wqy
    !$acc atomic update
    jy(j+2, k0+2, l  )  = jy(j+2, k0+2, l  )  +   sx(2 )*sy0(2 )*sz(0 )*wqy
    !$acc atomic update
    jy(j-1, k0-1, l+1)  = jy(j-1, k0-1, l+1)  +   sx(-1)*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0-1, l+1)  = jy(j, k0-1, l+1)  +   sx(0 )*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0-1, l+1)  = jy(j+1, k0-1, l+1)  +   sx(1 )*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j+2, k0-1, l+1)  = jy(j+2, k0-1, l+1)  +   sx(2 )*sy0(-1)*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, k0, l+1)  = jy(j-1, k0, l+1)  +   sx(-1)*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0, l+1)  = jy(j, k0, l+1)  +   sx(0 )*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0, l+1)  = jy(j+1, k0, l+1)  +   sx(1 )*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+2, k0, l+1)  = jy(j+2, k0, l+1)  +   sx(2 )*sy0(0 )*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, k0+1, l+1)  = jy(j-1, k0+1, l+1)  +   sx(-1)*sy0(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0+1, l+1)  = jy(j, k0+1, l+1)  +   sx(0 )*sy0(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0+1, l+1)  = jy(j+1, k0+1, l+1)  +   sx(1 )*sy0(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+2, k0+1, l+1)  = jy(j+2, k0+1, l+1)  +   sx(2 )*sy0(1 )*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, k0+2, l+1)  = jy(j-1, k0+2, l+1)  +   sx(-1)*sy0(2 )*sz(1 )*wqy
    !$acc atomic update
    jy(j, k0+2, l+1)  = jy(j, k0+2, l+1)  +   sx(0 )*sy0(2 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+1, k0+2, l+1)  = jy(j+1, k0+2, l+1)  +   sx(1 )*sy0(2 )*sz(1 )*wqy
    !$acc atomic update
    jy(j+2, k0+2, l+1)  = jy(j+2, k0+2, l+1)  +   sx(2 )*sy0(2 )*sz(1 )*wqy
    !$acc atomic update
    jy(j-1, k0-1, l+2)  = jy(j-1, k0-1, l+2)  +   sx(-1)*sy0(-1)*sz(2 )*wqy
    !$acc atomic update
    jy(j, k0-1, l+2)  = jy(j, k0-1, l+2)  +   sx(0 )*sy0(-1)*sz(2 )*wqy
    !$acc atomic update
    jy(j+1, k0-1, l+2)  = jy(j+1, k0-1, l+2)  +   sx(1 )*sy0(-1)*sz(2 )*wqy
    !$acc atomic update
    jy(j+2, k0-1, l+2)  = jy(j+2, k0-1, l+2)  +   sx(2 )*sy0(-1)*sz(2 )*wqy
    !$acc atomic update
    jy(j-1, k0, l+2)  = jy(j-1, k0, l+2)  +   sx(-1)*sy0(0 )*sz(2 )*wqy
    !$acc atomic update
    jy(j, k0, l+2)  = jy(j, k0, l+2)  +   sx(0 )*sy0(0 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+1, k0, l+2)  = jy(j+1, k0, l+2)  +   sx(1 )*sy0(0 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+2, k0, l+2)  = jy(j+2, k0, l+2)  +   sx(2 )*sy0(0 )*sz(2 )*wqy
    !$acc atomic update
    jy(j-1, k0+1, l+2)  = jy(j-1, k0+1, l+2)  +   sx(-1)*sy0(1 )*sz(2 )*wqy
    !$acc atomic update
    jy(j, k0+1, l+2)  = jy(j, k0+1, l+2)  +   sx(0 )*sy0(1 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+1, k0+1, l+2)  = jy(j+1, k0+1, l+2)  +   sx(1 )*sy0(1 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+2, k0+1, l+2)  = jy(j+2, k0+1, l+2)  +   sx(2 )*sy0(1 )*sz(2 )*wqy
    !$acc atomic update
    jy(j-1, k0+2, l+2)  = jy(j-1, k0+2, l+2)  +   sx(-1)*sy0(2 )*sz(2 )*wqy
    !$acc atomic update
    jy(j, k0+2, l+2)  = jy(j, k0+2, l+2)  +   sx(0 )*sy0(2 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+1, k0+2, l+2)  = jy(j+1, k0+2, l+2)  +   sx(1 )*sy0(2 )*sz(2 )*wqy
    !$acc atomic update
    jy(j+2, k0+2, l+2)  = jy(j+2, k0+2, l+2)  +   sx(2 )*sy0(2 )*sz(2 )*wqy

    ! - JZ
    !$acc atomic update
    jz(j-1, k-1, l0-1)  = jz(j-1, k-1, l0-1)  +   sx(-1)*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k-1, l0-1)  = jz(j, k-1, l0-1)  +   sx(0 )*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k-1, l0-1)  = jz(j+1, k-1, l0-1)  +   sx(1 )*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j+2, k-1, l0-1)  = jz(j+2, k-1, l0-1)  +   sx(2 )*sy(-1)*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k, l0-1)  = jz(j-1, k, l0-1)  +   sx(-1)*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k, l0-1)  = jz(j, k, l0-1)  +   sx(0 )*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k, l0-1)  = jz(j+1, k, l0-1)  +   sx(1 )*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+2, k, l0-1)  = jz(j+2, k, l0-1)  +   sx(2 )*sy(0 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k+1, l0-1)  = jz(j-1, k+1, l0-1)  +   sx(-1)*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k+1, l0-1)  = jz(j, k+1, l0-1)  +   sx(0 )*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k+1, l0-1)  = jz(j+1, k+1, l0-1)  +   sx(1 )*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+2, k+1, l0-1)  = jz(j+2, k+1, l0-1)  +   sx(2 )*sy(1 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k+2, l0-1)  = jz(j-1, k+2, l0-1)  +   sx(-1)*sy(2 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j, k+2, l0-1)  = jz(j, k+2, l0-1)  +   sx(0 )*sy(2 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+1, k+2, l0-1)  = jz(j+1, k+2, l0-1)  +   sx(1 )*sy(2 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j+2, k+2, l0-1)  = jz(j+2, k+2, l0-1)  +   sx(2 )*sy(2 )*sz0(-1)*wqz
    !$acc atomic update
    jz(j-1, k-1, l0  )  = jz(j-1, k-1, l0  )  +   sx(-1)*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k-1, l0  )  = jz(j, k-1, l0  )  +   sx(0 )*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k-1, l0  )  = jz(j+1, k-1, l0  )  +   sx(1 )*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j+2, k-1, l0  )  = jz(j+2, k-1, l0  )  +   sx(2 )*sy(-1)*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k, l0  )  = jz(j-1, k, l0  )  +   sx(-1)*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k, l0  )  = jz(j, k, l0  )  +   sx(0 )*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k, l0  )  = jz(j+1, k, l0  )  +   sx(1 )*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+2, k, l0  )  = jz(j+2, k, l0  )  +   sx(2 )*sy(0 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k+1, l0  )  = jz(j-1, k+1, l0  )  +   sx(-1)*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k+1, l0  )  = jz(j, k+1, l0  )  +   sx(0 )*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k+1, l0  )  = jz(j+1, k+1, l0  )  +   sx(1 )*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+2, k+1, l0  )  = jz(j+2, k+1, l0  )  +   sx(2 )*sy(1 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k+2, l0  )  = jz(j-1, k+2, l0  )  +   sx(-1)*sy(2 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j, k+2, l0  )  = jz(j, k+2, l0  )  +   sx(0 )*sy(2 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+1, k+2, l0  )  = jz(j+1, k+2, l0  )  +   sx(1 )*sy(2 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j+2, k+2, l0  )  = jz(j+2, k+2, l0  )  +   sx(2 )*sy(2 )*sz0(0 )*wqz
    !$acc atomic update
    jz(j-1, k-1, l0+1)  = jz(j-1, k-1, l0+1)  +   sx(-1)*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k-1, l0+1)  = jz(j, k-1, l0+1)  +   sx(0 )*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k-1, l0+1)  = jz(j+1, k-1, l0+1)  +   sx(1 )*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j+2, k-1, l0+1)  = jz(j+2, k-1, l0+1)  +   sx(2 )*sy(-1)*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, k, l0+1)  = jz(j-1, k, l0+1)  +   sx(-1)*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k, l0+1)  = jz(j, k, l0+1)  +   sx(0 )*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k, l0+1)  = jz(j+1, k, l0+1)  +   sx(1 )*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+2, k, l0+1)  = jz(j+2, k, l0+1)  +   sx(2 )*sy(0 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, k+1, l0+1)  = jz(j-1, k+1, l0+1)  +   sx(-1)*sy(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k+1, l0+1)  = jz(j, k+1, l0+1)  +   sx(0 )*sy(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k+1, l0+1)  = jz(j+1, k+1, l0+1)  +   sx(1 )*sy(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+2, k+1, l0+1)  = jz(j+2, k+1, l0+1)  +   sx(2 )*sy(1 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, k+2, l0+1)  = jz(j-1, k+2, l0+1)  +   sx(-1)*sy(2 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j, k+2, l0+1)  = jz(j, k+2, l0+1)  +   sx(0 )*sy(2 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+1, k+2, l0+1)  = jz(j+1, k+2, l0+1)  +   sx(1 )*sy(2 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j+2, k+2, l0+1)  = jz(j+2, k+2, l0+1)  +   sx(2 )*sy(2 )*sz0(1 )*wqz
    !$acc atomic update
    jz(j-1, k-1, l0+2)  = jz(j-1, k-1, l0+2)  +   sx(-1)*sy(-1)*sz0(2 )*wqz
    !$acc atomic update
    jz(j, k-1, l0+2)  = jz(j, k-1, l0+2)  +   sx(0 )*sy(-1)*sz0(2 )*wqz
    !$acc atomic update
    jz(j+1, k-1, l0+2)  = jz(j+1, k-1, l0+2)  +   sx(1 )*sy(-1)*sz0(2 )*wqz
    !$acc atomic update
    jz(j+2, k-1, l0+2)  = jz(j+2, k-1, l0+2)  +   sx(2 )*sy(-1)*sz0(2 )*wqz
    !$acc atomic update
    jz(j-1, k, l0+2)  = jz(j-1, k, l0+2)  +   sx(-1)*sy(0 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j, k, l0+2)  = jz(j, k, l0+2)  +   sx(0 )*sy(0 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+1, k, l0+2)  = jz(j+1, k, l0+2)  +   sx(1 )*sy(0 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+2, k, l0+2)  = jz(j+2, k, l0+2)  +   sx(2 )*sy(0 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j-1, k+1, l0+2)  = jz(j-1, k+1, l0+2)  +   sx(-1)*sy(1 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j, k+1, l0+2)  = jz(j, k+1, l0+2)  +   sx(0 )*sy(1 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+1, k+1, l0+2)  = jz(j+1, k+1, l0+2)  +   sx(1 )*sy(1 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+2, k+1, l0+2)  = jz(j+2, k+1, l0+2)  +   sx(2 )*sy(1 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j-1, k+2, l0+2)  = jz(j-1, k+2, l0+2)  +   sx(-1)*sy(2 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j, k+2, l0+2)  = jz(j, k+2, l0+2)  +   sx(0 )*sy(2 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+1, k+2, l0+2)  = jz(j+1, k+2, l0+2)  +   sx(1 )*sy(2 )*sz0(2 )*wqz
    !$acc atomic update
    jz(j+2, k+2, l0+2)  = jz(j+2, k+2, l0+2)  +   sx(2 )*sy(2 )*sz0(2 )*wqz

  END DO
  !$acc end loop
  !$acc end parallel
  RETURN
END SUBROUTINE depose_jxjyjz_scalar_3_3_3

#if defined (DEV)
! _______________________________________________________________________________________
!> @brief
!> Order 3 3D vector current deposition routine (rho*v)
!
!> @details
!> This versions have good performances on SIMD architectures
!> Providing that OpenMP 4.0 is available (Directive SIMD)
!> Use with nox=4
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHVv2_3_3_3(jx, jy, jz, np, xp, yp, zp, uxp, uyp, uzp,     &
  gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz, nxguard, nyguard,         &
  nzguard, l_nodal)  !#do not wrap
  USE constants, ONLY: clight, lvec
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, ny, nz, nxguard, nyguard, nzguard
  LOGICAL(idp) :: l_nodal
  REAL(num)    :: stagger_shift
  REAL(num), INTENT(IN OUT) ::                                                        &
  jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, oxint, oyint, xintsq, yintsq, oxintsq,      &
  oyintsq
  REAL(num) :: x, y, z, xmid, ymid, zmid, invvol, dts2dx, dts2dy, dts2dz
  REAL(num) ::   ww, wwx, wwy, wwz, usq, clightsq
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, j0, k0, l0, ip, NCELLS, ic, ix, iy, iz
  INTEGER(idp) :: nnx, nnxy, ngridx, ngridy, n, nn, nv
  INTEGER(idp) :: moff(1:8)
  INTEGER(idp), PARAMETER :: LVEC2=32
  REAL(num) :: zint(LVEC), zint0(LVEC2)
  INTEGER(idp), DIMENSION(LVEC2, 3) :: ICELL
  REAL(num), DIMENSION(LVEC2) :: vx, vy, vz
  REAL(num) ::  wwwx(LVEC2, 16), wwwy(LVEC2, 16), wwwz(LVEC2, 16), wq
  REAL(num) :: sx1(LVEC2), sx2(LVEC2), sx3(LVEC2), sx4(LVEC2)
  REAL(num) :: sx01(LVEC2), sx02(LVEC2), sx03(LVEC2), sx04(LVEC2)
  REAL(num) :: sy1(LVEC2), sy2(LVEC2), sy3(LVEC2), sy4(LVEC2)
  REAL(num) :: sy01(LVEC2), sy02(LVEC2), sy03(LVEC2), sy04(LVEC2)
  REAL(num), DIMENSION(4) :: szz, zdec, h1, h11, h12, sgn
  INTEGER(idp) :: orig, ncxy, ncx, ncy, ncz, ngx, ngxy, igrid, jorig, korig, lorig

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
  ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
  ncx=nx+5; ncy=ny+5; ncz=nz+5
  NCELLS=ncx*ncy*ncz
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
  jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num;
  nnx = ngridx
  nnxy = ngridx*ngridy
  moff = (/-nnxy, 0_idp, nnxy, 2_idp*nnxy, nnx-nnxy, nnx, nnx+nnxy, nnx+2_idp*nnxy/)
  jorig=-3_idp; korig=-3_idp;lorig=-3_idp
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  h1=(/1_num, 0_num, 1_num, 0_num/); sgn=(/1_num, -1_num, 1_num, -1_num/)
  h11=(/0_num, 1_num, 1_num, 0_num/); h12=(/1_num, 0_num, 0_num, 1_num/)
  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC2
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED vx:64, vy:64, vz:64
    !DIR$ ASSUME_ALIGNED sx1:64, sx2:64, sx3:64, sx4:64
    !DIR$ ASSUME_ALIGNED sy1:64, sy2:64, sy3:64, sy4:64
    !DIR$ ASSUME_ALIGNED sx01:64, sx02:64, sx03:64, sx04:64
    !DIR$ ASSUME_ALIGNED sy01:64, sy02:64, sy03:64, sy04:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, vx, vy, vz)
    !IBM* ALIGN(64, sx1, sx2, sx3, sx4)
    !IBM* ALIGN(64, sy1, sy2, sy3, sy4)
    !IBM* ALIGN(64, sx01, sx02, sx03, sx04)
    !IBM* ALIGN(64, sy01, sy02, sy03, sy04)
    !IBM* ALIGN(64, ICELL:64)
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
      ! --- computes position in  grid units at (n+1)
      x = (xp(nn)-xmin)*dxi
      y = (yp(nn)-ymin)*dyi
      z = (zp(nn)-zmin)*dzi

      ! Computes velocity
      vx(n) = uxp(nn)*gaminv(nn)
      vy(n) = uyp(nn)*gaminv(nn)
      vz(n) = uzp(nn)*gaminv(nn)

      ! --- computes particles weights
      wq=q*w(nn)*invvol

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx(n)
      ymid=y-dts2dy*vy(n)
      zmid=z-dts2dz*vz(n)

      ! --- finds node of cell containing particles for current positions
      j=floor(xmid)
      k=floor(ymid)
      l=floor(zmid)
      j0=floor(xmid-stagger_shift)
      k0=floor(ymid-stagger_shift)
      l0=floor(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! --- computes set of coefficients for node centered quantities
      xint    = xmid-j
      yint    = ymid-k
      zint(n) = zmid-l
      oxint   = 1.0_num-xint
      xintsq  = xint*xint
      oxintsq = oxint*oxint
      sx1(n)  = onesixth*oxintsq*oxint
      sx2(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx3(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx4(n)  = onesixth*xintsq*xint
      oyint   = 1.0_num-yint
      yintsq  = yint*yint
      oyintsq = oyint*oyint
      sy1(n)  = onesixth*oyintsq*oyint*wq
      sy2(n)  = (twothird-yintsq*(1.0_num-yint*0.5_num))*wq
      sy3(n)  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))*wq
      sy4(n)  = onesixth*yintsq*yint*wq

      ! --- computes set of coefficients for staggered quantities
      xint     = xmid-stagger_shift-j0
      yint     = ymid-stagger_shift-k0
      zint0(n) = zmid-stagger_shift-l0
      oxint    = 1.0_num-xint
      xintsq   = xint*xint
      oxintsq  = oxint*oxint
      sx01(n)  = onesixth*oxintsq*oxint
      sx02(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx03(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx04(n)  = onesixth*xintsq*xint
      oyint    = 1.0_num-yint
      yintsq   = yint*yint
      oyintsq  = oyint*oyint
      sy01(n)  = onesixth*oyintsq*oyint*wq
      sy02(n)  = (twothird-yintsq*(1.0_num-yint*0.5_num))*wq
      sy03(n)  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))*wq
      sy04(n)  = onesixth*yintsq*yint*wq
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
    ! Compute weights
    DO n=1, MIN(LVEC2, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED w:64, wwwx:64, wwwy:64, wwwz:64
#elif defined  __IBMBGQ__
      !IBM* ALIGN(64, w, wwwx, wwwy, wwwz)
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
        ! - Weiths for jx
        zdec(nv)      = (h1(nv)-zint(n))*sgn(nv)
        szz(nv)       = (twothird-zdec(nv)**2*(1.0_num-zdec(nv)*0.5_num))*h11(nv)     &
        +onesixth*zdec(nv)**3*h12(nv)
        wwwx(nv, n)    = szz(nv)*sy1(n)*vx(n)
        wwwx(nv+4, n)  = szz(nv)*sy2(n)*vx(n)
        wwwx(nv+8, n)  = szz(nv)*sy3(n)*vx(n)
        wwwx(nv+12, n) = szz(nv)*sy4(n)*vx(n)
        ! - Weiths for jy
        wwwy(nv, n)    = szz(nv)*sy01(n)*vy(n)
        wwwy(nv+4, n)  = szz(nv)*sy02(n)*vy(n)
        wwwy(nv+8, n)  = szz(nv)*sy03(n)*vy(n)
        wwwy(nv+12, n) = szz(nv)*sy04(n)*vy(n)
        ! - Weiths for jz
        zdec(nv)      = (h1(nv)-zint0(n))*sgn(nv)
        szz(nv)       = (twothird-zdec(nv)**2*(1.0_num-zdec(nv)*0.5_num))*h11(nv)     &
        +onesixth*zdec(nv)**3*h12(nv)
        wwwz(nv, n)    = szz(nv)*sy1(n)*vz(n)
        wwwz(nv+4, n)  = szz(nv)*sy2(n)*vz(n)
        wwwz(nv+8, n)  = szz(nv)*sy3(n)*vz(n)
        wwwz(nv+12, n) = szz(nv)*sy4(n)*vz(n)

      ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO

    ! Add weights to nearest vertices
    DO n=1, MIN(LVEC2, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx-1) = jxcells(nv, ICELL(n, 1)-ncx-1) + wwwx(nv,    &
        n)*sx01(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1)-ncx)   = jxcells(nv, ICELL(n, 1)-ncx)   + wwwx(nv,    &
        n)*sx02(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx+1) = jxcells(nv, ICELL(n, 1)-ncx+1) + wwwx(nv,    &
        n)*sx03(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx+2) = jxcells(nv, ICELL(n, 1)-ncx+2) + wwwx(nv,    &
        n)*sx04(n)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx-1) = jxcells(nv, ICELL(n, 1)+ncx-1) + wwwx(nv+8,  &
        n)*sx01(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1)+ncx)   = jxcells(nv, ICELL(n, 1)+ncx)   + wwwx(nv+8,  &
        n)*sx02(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx+1) = jxcells(nv, ICELL(n, 1)+ncx+1) + wwwx(nv+8,  &
        n)*sx03(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx+2) = jxcells(nv, ICELL(n, 1)+ncx+2) + wwwx(nv+8,  &
        n)*sx04(n)

        ! --- JY
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)-ncx-1) = jycells(nv, ICELL(n, 2)-ncx-1) + wwwy(nv,    &
        n)*sx1(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2)-ncx)   = jycells(nv, ICELL(n, 2)-ncx)   + wwwy(nv,    &
        n)*sx2(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)-ncx+1) = jycells(nv, ICELL(n, 2)-ncx+1) + wwwy(nv,    &
        n)*sx3(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)-ncx+2) = jycells(nv, ICELL(n, 2)-ncx+2) + wwwy(nv,    &
        n)*sx4(n)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)+ncx-1) = jycells(nv, ICELL(n, 2)+ncx-1) + wwwy(nv+8,  &
        n)*sx1(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2)+ncx)   = jycells(nv, ICELL(n, 2)+ncx)   + wwwy(nv+8,  &
        n)*sx2(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+ncx+1) = jycells(nv, ICELL(n, 2)+ncx+1) + wwwy(nv+8,  &
        n)*sx3(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+ncx+2) = jycells(nv, ICELL(n, 2)+ncx+2) + wwwy(nv+8,  &
        n)*sx4(n)

        ! --- JZ
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx-1) = jzcells(nv, ICELL(n, 3)-ncx-1) + wwwz(nv,    &
        n)*sx1(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3)-ncx)   = jzcells(nv, ICELL(n, 3)-ncx)   + wwwz(nv,    &
        n)*sx2(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx+1) = jzcells(nv, ICELL(n, 3)-ncx+1) + wwwz(nv,    &
        n)*sx3(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx+2) = jzcells(nv, ICELL(n, 3)-ncx+2) + wwwz(nv,    &
        n)*sx4(n)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx-1) = jzcells(nv, ICELL(n, 3)+ncx-1) + wwwz(nv+8,  &
        n)*sx1(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3)+ncx)   = jzcells(nv, ICELL(n, 3)+ncx)   + wwwz(nv+8,  &
        n)*sx2(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx+1) = jzcells(nv, ICELL(n, 3)+ncx+1) + wwwz(nv+8,  &
        n)*sx3(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx+2) = jzcells(nv, ICELL(n, 3)+ncx+2) + wwwz(nv+8,  &
        n)*sx4(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
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
        ! jx
        jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1, ic)
        jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2, ic)
        jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3, ic)
        jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4, ic)
        jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5, ic)
        jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6, ic)
        jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7, ic)
        jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8, ic)
        ! jy
        jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1, ic)
        jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2, ic)
        jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3, ic)
        jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4, ic)
        jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5, ic)
        jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6, ic)
        jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7, ic)
        jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8, ic)
        ! jz
        jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1, ic)
        jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2, ic)
        jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3, ic)
        jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4, ic)
        jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5, ic)
        jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6, ic)
        jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7, ic)
        jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8, ic)
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
END SUBROUTINE depose_jxjyjz_vecHVv2_3_3_3
#endif

! ______________________________________________________________________________
!> @brief
!> Order 3 3D vector current deposition routine (rho*v)
!>
!> @details
!> This versions have good performances on SIMD architectures
!> Providing that OpenMP 4.0 is available (Directive SIMD)
!> Use with nox=3
!>
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revison 10/09/2016
!>
!> @param[inout] jx, jy, jz current arrays
!> @param[in] np number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv particle Lorentz factor arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin, ymin, zmin tile grid minimum position
!> @param[in] dx, dy, dz space discretization steps
!> @param[in] jx_nguard number of guard cells of the jx array in each direction
!> (1d array containing 3 integers)
!> @param[in] jx_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jx array (1d array containing 3 integers)
!> @param[inout] jy y-current component (3D array)
!> @param[in] jy_nguard number of guard cells of the jy array in each direction
!>  (1d array containing 3 integers)
!> @param[in] jy_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jy array (1d array containing 3 integers)
!> @param[inout] jz z-current component (3D array)
!> @param[in] jz_nguard number of guard cells of the jz array in each direction
!> (1d array containing 3 integers)
!> @param[in] jz_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jz array (1d array containing 3 integers)
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHVv3_3_3_3( jx, jx_nguard, jx_nvalid, jy, jy_nguard,      &
  jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, ymin, zmin, dt, dx, dy, dz, l_nodal)     !#do not wrap
  USE constants, ONLY: clight, lvec
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  ! ___ Parameter declaration _______________________________________
  INTEGER(idp)             :: np
  LOGICAL(idp)             :: l_nodal
  REAL(num)                :: stagger_shift
  INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), &
  jy_nguard(3), jy_nvalid(3), jz_nguard(3), jz_nvalid(3)
  REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1, &
  -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1, &
  -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1, &
  -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1, &
  -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )
  REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1, &
  -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1, &
  -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )
  REAL(num), DIMENSION(:, :), ALLOCATABLE:: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dy, dz, xmin, ymin, zmin
  REAL(num) :: dxi, dyi, dzi, xint, yint, zint, oxint, oyint, ozint, xintsq, yintsq,  &
  zintsq, oxintsq, oyintsq, ozintsq
  REAL(num) :: x, y, z, xmid, ymid, zmid, invvol, dts2dx, dts2dy, dts2dz
  REAL(num) ::   clightsq
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num, twothird=2.0_num/3.0_num
  INTEGER(idp) :: j, k, l, j0, k0, l0, ip, NCELLS, ic
  INTEGER(idp) :: n, nn, nv
  INTEGER(idp), DIMENSION(LVEC, 3) :: ICELL
  REAL(num), DIMENSION(LVEC) :: vx, vy, vz
  REAL(num) :: wq
  REAL(num) :: sx1(LVEC), sx2(LVEC), sx3(LVEC), sx4(LVEC)
  REAL(num) :: sx01(LVEC), sx02(LVEC), sx03(LVEC), sx04(LVEC)
  REAL(num) :: sy1, sy2, sy3, sy4, sz1, sz2, sz3, sz4
  REAL(num) :: sy01, sy02, sy03, sy04, sz01, sz02, sz03, sz04
  REAL(num), DIMENSION(4) :: h1, h11, h12, sgn
  REAL(num) :: wwwx1(LVEC, 8), wwwx2(LVEC, 8), wwwy1(LVEC, 8), wwwy2(LVEC, 8),         &
  wwwz1(LVEC, 8), wwwz2(LVEC, 8)
  REAL(num) :: wx1, wx2, wy1, wy2, wz1, wz2
  INTEGER(idp) :: jorig, korig, lorig
  INTEGER(idp) :: ncx, ncy, ncxy, ncz, ix, iy, iz
  ! Relative position of the 8 nodes (over which the algorithm vectorizes)
  ! in respect to the particle computed node (the nodes are only shifted in
  ! y and z, not in x ; see arXiv:1601.02056, Fig. 3, for order 3)
  INTEGER(idp), DIMENSION(8), PARAMETER  :: mzoff = &
  (/ -1_idp, 0_idp, 1_idp, 2_idp, -1_idp, 0_idp, 1_idp, 2_idp /)
  INTEGER(idp), DIMENSION(8), PARAMETER  :: myoff = &
  (/  0_idp, 0_idp, 0_idp, 0_idp,  1_idp, 1_idp, 1_idp, 1_idp /)

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi
  clightsq = 1.0_num/clight**2
  ! Lower integer position where the particles can deposit
  jorig=-2_idp; korig=-2_idp; lorig=-2_idp
  ! Find the maximal number of cells in each direction, so as to
  ! allocate the linearized 1D arrays
  ncx = MAX( jx_nvalid(1), jy_nvalid(1), jz_nvalid(1) ) + 5
  ncy = MAX( jx_nvalid(2), jy_nvalid(2), jz_nvalid(2) ) + 4
  ncz = MAX( jx_nvalid(3), jy_nvalid(3), jz_nvalid(3) ) + 3
  ncxy = ncx*ncy
  NCELLS = ncx*ncy*ncz
  ALLOCATE(jxcells(8, NCELLS), jycells(8, NCELLS), jzcells(8, NCELLS))
  jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num;

  h1=(/1_num, 0_num, 1_num, 0_num/); sgn=(/1_num, -1_num, 1_num, -1_num/)
  h11=(/0_num, 1_num, 1_num, 0_num/); h12=(/1_num, 0_num, 0_num, 1_num/)
  ! LOOP ON PARTICLES
  DO ip=1, np, LVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED vx:64, vy:64, vz:64
    !DIR$ ASSUME_ALIGNED gaminv:64
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

      ! Computes velocity
      vx(n) = uxp(nn)*gaminv(nn)
      vy(n) = uyp(nn)*gaminv(nn)
      vz(n) = uzp(nn)*gaminv(nn)

      ! --- computes particles weights
      wq=q*w(nn)*invvol

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx(n)
      ymid=y-dts2dy*vy(n)
      zmid=z-dts2dz*vz(n)

      ! --- finds node of cell containing particles for current positions
      j=floor(xmid)
      k=floor(ymid)
      l=floor(zmid)
      j0=floor(xmid-stagger_shift)
      k0=floor(ymid-stagger_shift)
      l0=floor(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! --- computes set of coefficients for node centered quantities
      xint    = xmid-j
      yint    = ymid-k
      zint    = zmid-l
      oxint   = 1.0_num-xint
      xintsq  = xint*xint
      oxintsq = oxint*oxint
      sx1(n)  = onesixth*oxintsq*oxint
      sx2(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx3(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx4(n)  = onesixth*xintsq*xint
      oyint   = 1.0_num-yint
      yintsq  = yint*yint
      oyintsq = oyint*oyint
      sy1     = onesixth*oyintsq*oyint
      sy2  = (twothird-yintsq*(1.0_num-yint*0.5_num))
      sy3  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
      sy4  = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz1 = onesixth*ozintsq*ozint*wq
      sz2 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
      sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
      sz4 = onesixth*zintsq*zint*wq

      ! --- computes set of coefficients for staggered quantities
      xint     = xmid-stagger_shift-j0
      yint     = ymid-stagger_shift-k0
      zint     = zmid-stagger_shift-l0
      oxint    = 1.0_num-xint
      xintsq   = xint*xint
      oxintsq  = oxint*oxint
      sx01(n)  = onesixth*oxintsq*oxint
      sx02(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx03(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx04(n)  = onesixth*xintsq*xint
      oyint    = 1.0_num-yint
      yintsq   = yint*yint
      oyintsq  = oyint*oyint
      sy01     = onesixth*oyintsq*oyint
      sy02     = (twothird-yintsq*(1.0_num-yint*0.5_num))
      sy03  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
      sy04  = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz01 = onesixth*ozintsq*ozint*wq
      sz02 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
      sz03 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
      sz04 = onesixth*zintsq*zint*wq
      ! --- computes weights
      ! - X
      wwwx1(n, 1)=sz1*sy1
      wwwx1(n, 2)=sz2*sy1
      wwwx1(n, 3)=sz3*sy1
      wwwx1(n, 4)=sz4*sy1
      wwwx1(n, 5)=sz1*sy2
      wwwx1(n, 6)=sz2*sy2
      wwwx1(n, 7)=sz3*sy2
      wwwx1(n, 8)=sz4*sy2
      wwwx2(n, 1)=sz1*sy3
      wwwx2(n, 2)=sz2*sy3
      wwwx2(n, 3)=sz3*sy3
      wwwx2(n, 4)=sz4*sy3
      wwwx2(n, 5)=sz1*sy4
      wwwx2(n, 6)=sz2*sy4
      wwwx2(n, 7)=sz3*sy4
      wwwx2(n, 8)=sz4*sy4
      ! - Y
      wwwy1(n, 1)=sz1*sy01
      wwwy1(n, 2)=sz2*sy01
      wwwy1(n, 3)=sz3*sy01
      wwwy1(n, 4)=sz4*sy01
      wwwy1(n, 5)=sz1*sy02
      wwwy1(n, 6)=sz2*sy02
      wwwy1(n, 7)=sz3*sy02
      wwwy1(n, 8)=sz4*sy02
      wwwy2(n, 1)=sz1*sy03
      wwwy2(n, 2)=sz2*sy03
      wwwy2(n, 3)=sz3*sy03
      wwwy2(n, 4)=sz4*sy03
      wwwy2(n, 5)=sz1*sy04
      wwwy2(n, 6)=sz2*sy04
      wwwy2(n, 7)=sz3*sy04
      wwwy2(n, 8)=sz4*sy04
      ! - Z
      wwwz1(n, 1)=sz01*sy1
      wwwz1(n, 2)=sz02*sy1
      wwwz1(n, 3)=sz03*sy1
      wwwz1(n, 4)=sz04*sy1
      wwwz1(n, 5)=sz01*sy2
      wwwz1(n, 6)=sz02*sy2
      wwwz1(n, 7)=sz03*sy2
      wwwz1(n, 8)=sz04*sy2
      wwwz2(n, 1)=sz01*sy3
      wwwz2(n, 2)=sz02*sy3
      wwwz2(n, 3)=sz03*sy3
      wwwz2(n, 4)=sz04*sy3
      wwwz2(n, 5)=sz01*sy4
      wwwz2(n, 6)=sz02*sy4
      wwwz2(n, 7)=sz03*sy4
      wwwz2(n, 8)=sz04*sy4
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! Add weights to nearest vertices
    DO n=1, MIN(LVEC, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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
        wx1=wwwx1(n, nv); wx2=wwwx2(n, nv)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx-1) = jxcells(nv, ICELL(n, 1)-ncx-1) +             &
        wx1*sx01(n)*vx(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1)-ncx)   = jxcells(nv, ICELL(n, 1)-ncx)   +             &
        wx1*sx02(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx+1) = jxcells(nv, ICELL(n, 1)-ncx+1) +             &
        wx1*sx03(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx+2) = jxcells(nv, ICELL(n, 1)-ncx+2) +             &
        wx1*sx04(n)*vx(n)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx-1) = jxcells(nv, ICELL(n, 1)+ncx-1) +             &
        wx2*sx01(n)*vx(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1)+ncx)   = jxcells(nv, ICELL(n, 1)+ncx)   +             &
        wx2*sx02(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx+1) = jxcells(nv, ICELL(n, 1)+ncx+1) +             &
        wx2*sx03(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx+2) = jxcells(nv, ICELL(n, 1)+ncx+2) +             &
        wx2*sx04(n)*vx(n)

        ! --- JY
        wy1=wwwy1(n, nv); wy2=wwwy2(n, nv)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)-ncx-1) = jycells(nv, ICELL(n, 2)-ncx-1) +             &
        wy1*sx1(n)*vy(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2)-ncx)   = jycells(nv, ICELL(n, 2)-ncx)   +             &
        wy1*sx2(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)-ncx+1) = jycells(nv, ICELL(n, 2)-ncx+1) +             &
        wy1*sx3(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)-ncx+2) = jycells(nv, ICELL(n, 2)-ncx+2) +             &
        wy1*sx4(n)*vy(n)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)+ncx-1) = jycells(nv, ICELL(n, 2)+ncx-1) +             &
        wy2*sx1(n)*vy(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2)+ncx)   = jycells(nv, ICELL(n, 2)+ncx)   +             &
        wy2*sx2(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+ncx+1) = jycells(nv, ICELL(n, 2)+ncx+1) +             &
        wy2*sx3(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+ncx+2) = jycells(nv, ICELL(n, 2)+ncx+2) +             &
        wy2*sx4(n)*vy(n)

        ! --- JZ
        wz1=wwwz1(n, nv); wz2=wwwz2(n, nv)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx-1) = jzcells(nv, ICELL(n, 3)-ncx-1) +             &
        wz1*sx1(n)*vz(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3)-ncx)   = jzcells(nv, ICELL(n, 3)-ncx)   +             &
        wz1*sx2(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx+1) = jzcells(nv, ICELL(n, 3)-ncx+1) +             &
        wz1*sx3(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx+2) = jzcells(nv, ICELL(n, 3)-ncx+2) +             &
        wz1*sx4(n)*vz(n)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx-1) = jzcells(nv, ICELL(n, 3)+ncx-1) +             &
        wz2*sx1(n)*vz(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3)+ncx)   = jzcells(nv, ICELL(n, 3)+ncx)   +             &
        wz2*sx2(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx+1) = jzcells(nv, ICELL(n, 3)+ncx+1) +             &
        wz2*sx3(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx+2) = jzcells(nv, ICELL(n, 3)+ncx+2) +             &
        wz2*sx4(n)*vz(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  ! ---------------------------------------------------------------
  ! Vectorized reduction of jxcells, jycells, jzcells in jx, jy, jz
  ! ---------------------------------------------------------------

  ! Reduction of jxcells in jx
  ! Note: the bounds below make sure that we never go out-of-bound for jxcells,
  ! even when taking into account the additional offsets in myoff and mzoff
  DO iz = MAX(lorig, -jx_nguard(3)+1), MIN(lorig+ncz-1, jx_nvalid(3)+jx_nguard(3)-2)
    DO iy = MAX(korig, -jx_nguard(2)), MIN(korig+ncy-1, jx_nvalid(2)+jx_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
       DO ix = MAX(jorig, -jx_nguard(1)), MIN(jorig+ncx-1, jx_nvalid(1)+jx_nguard(1)-2)
           ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
           ! Compute linearized index
           jx(ix,iy+myoff(1),iz+mzoff(1))=jx(ix,iy+myoff(1),iz+mzoff(1))+jxcells(1, ic)
           jx(ix,iy+myoff(2),iz+mzoff(2))=jx(ix,iy+myoff(2),iz+mzoff(2))+jxcells(2, ic)
           jx(ix,iy+myoff(3),iz+mzoff(3))=jx(ix,iy+myoff(3),iz+mzoff(3))+jxcells(3, ic)
           jx(ix,iy+myoff(4),iz+mzoff(4))=jx(ix,iy+myoff(4),iz+mzoff(4))+jxcells(4, ic)
           jx(ix,iy+myoff(5),iz+mzoff(5))=jx(ix,iy+myoff(5),iz+mzoff(5))+jxcells(5, ic)
           jx(ix,iy+myoff(6),iz+mzoff(6))=jx(ix,iy+myoff(6),iz+mzoff(6))+jxcells(6, ic)
           jx(ix,iy+myoff(7),iz+mzoff(7))=jx(ix,iy+myoff(7),iz+mzoff(7))+jxcells(7, ic)
           jx(ix,iy+myoff(8),iz+mzoff(8))=jx(ix,iy+myoff(8),iz+mzoff(8))+jxcells(8, ic)
       END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jycells in jy
  ! Note: the bounds below make sure that we never go out-of-bound for jycells,
  ! even when taking into account the additional offsets in myoff and mzoff
  DO iz = MAX(lorig, -jy_nguard(3)+1), MIN(lorig+ncz-1, jy_nvalid(3)+jy_nguard(3)-2)
    DO iy = MAX(korig, -jy_nguard(2)), MIN(korig+ncy-1, jy_nvalid(2)+jy_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
       DO ix = MAX(jorig, -jy_nguard(1)), MIN(jorig+ncx-1, jy_nvalid(1)+jy_nguard(1)-2)
           ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
           ! Compute linearized index
           jy(ix,iy+myoff(1),iz+mzoff(1))=jy(ix,iy+myoff(1),iz+mzoff(1))+jycells(1, ic)
           jy(ix,iy+myoff(2),iz+mzoff(2))=jy(ix,iy+myoff(2),iz+mzoff(2))+jycells(2, ic)
           jy(ix,iy+myoff(3),iz+mzoff(3))=jy(ix,iy+myoff(3),iz+mzoff(3))+jycells(3, ic)
           jy(ix,iy+myoff(4),iz+mzoff(4))=jy(ix,iy+myoff(4),iz+mzoff(4))+jycells(4, ic)
           jy(ix,iy+myoff(5),iz+mzoff(5))=jy(ix,iy+myoff(5),iz+mzoff(5))+jycells(5, ic)
           jy(ix,iy+myoff(6),iz+mzoff(6))=jy(ix,iy+myoff(6),iz+mzoff(6))+jycells(6, ic)
           jy(ix,iy+myoff(7),iz+mzoff(7))=jy(ix,iy+myoff(7),iz+mzoff(7))+jycells(7, ic)
           jy(ix,iy+myoff(8),iz+mzoff(8))=jy(ix,iy+myoff(8),iz+mzoff(8))+jycells(8, ic)
       END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  ! Reduction of jzcells in jz
  ! Note: the bounds below make sure that we never go out-of-bound for jzcells,
  ! even when taking into account the additional offsets in myoff and mzoff
  DO iz = MAX(lorig, -jz_nguard(3)+1), MIN(lorig+ncz-1, jz_nvalid(3)+jz_nguard(3)-2)
    DO iy = MAX(korig, -jz_nguard(2)), MIN(korig+ncy-1, jz_nvalid(2)+jz_nguard(2)-2)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP SIMD
#endif
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
      !$DIR SIMD
#endif
       DO ix = MAX(jorig, -jz_nguard(1)), MIN(jorig+ncx-1, jz_nvalid(1)+jz_nguard(1)-2)
           ic = 1 + (ix-jorig) + (iy-korig)*ncx + (iz-lorig)*ncxy
           ! Compute linearized index
           jz(ix,iy+myoff(1),iz+mzoff(1))=jz(ix,iy+myoff(1),iz+mzoff(1))+jzcells(1, ic)
           jz(ix,iy+myoff(2),iz+mzoff(2))=jz(ix,iy+myoff(2),iz+mzoff(2))+jzcells(2, ic)
           jz(ix,iy+myoff(3),iz+mzoff(3))=jz(ix,iy+myoff(3),iz+mzoff(3))+jzcells(3, ic)
           jz(ix,iy+myoff(4),iz+mzoff(4))=jz(ix,iy+myoff(4),iz+mzoff(4))+jzcells(4, ic)
           jz(ix,iy+myoff(5),iz+mzoff(5))=jz(ix,iy+myoff(5),iz+mzoff(5))+jzcells(5, ic)
           jz(ix,iy+myoff(6),iz+mzoff(6))=jz(ix,iy+myoff(6),iz+mzoff(6))+jzcells(6, ic)
           jz(ix,iy+myoff(7),iz+mzoff(7))=jz(ix,iy+myoff(7),iz+mzoff(7))+jzcells(7, ic)
           jz(ix,iy+myoff(8),iz+mzoff(8))=jz(ix,iy+myoff(8),iz+mzoff(8))+jzcells(8, ic)
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
END SUBROUTINE depose_jxjyjz_vecHVv3_3_3_3

! ________________________________________________________________________________________
!> @brief
!> Order 3 3D vector current deposition routine (rho*v).
!
!> @details
!> This versions have good performances on SIMD architectures
!> Providing that OpenMP 4.0 is available (Directive SIMD)
!> This subroutine is similar to depose_jxjyjz_vecHVv2_1_1_1
!> without the reduction process at the end
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!
!> @date
!> Creation 2016
!
! Input / output arguments:
!> @param[inout] jxcells, jycells, jzcells temporary current arrays
!> @param[in] np particle number
!> @param[in] ncells number of cells in the tile
!> @param[in] xp, yp, zp particle position arrays
!> @param[in] uxp, uyp, uzp particle momentum arrays
!> @param[in] gaminv inverse Lorentz factor arrays
!> @param[in] w particle wight arrays
!> @param[in] q charge
!> @param[in] xmin, ymin, zmin tile minimum positions
!> @param[in] dt, dx, dy, dz time and space steps
!> @param[in] nx, ny, nz tile cell numbers in each direction
!> @param[in] nxguard, nyguard, nzguard guard cells
!> @param[in] ncx, ncy, ncz tile cell extended number (depends on the order)
!> @param[in] lvect vector size
! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3(jxcells, jycells, jzcells, np, ncells, xp,   &
  yp, zp, uxp, uyp, uzp, gaminv, w, q, xmin, ymin, zmin, dt, dx, dy, dz, nx, ny, nz,    &
  nxguard, nyguard, nzguard, ncx, ncy, ncz, lvect, l_nodal)  !#do not wrap
  USE constants, ONLY: lvec
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE

  INTEGER(idp), INTENT(IN)                      :: np, nx, ny, nz, ncells
 LOGICAL(idp)                      :: l_nodal
 REAL(num)                         :: stagger_shift
  INTEGER(idp), INTENT(IN)                      :: lvect
  INTEGER(idp), INTENT(IN)                      :: nxguard, nyguard, nzguard
  REAL(num), DIMENSION(8, ncells), INTENT(INOUT) :: jxcells, jycells, jzcells
  REAL(num), DIMENSION(np), INTENT(IN)          :: xp, yp, zp, uxp, uyp, uzp, gaminv, &
  w
  REAL(num), INTENT(IN)                         :: q, dt, dx, dy, dz, xmin, ymin,     &
  zmin
  INTEGER(idp), INTENT(IN)                      :: ncx, ncy, ncz

  REAL(num)                                     :: xint, yint, zint, oxint, oyint,    &
  ozint, xintsq, yintsq, zintsq, oxintsq, oyintsq, ozintsq
  REAL(num)                                     :: x, y, z, xmid, ymid, zmid
  REAL(num), PARAMETER                          :: onesixth=1.0_num/6.0_num
  REAL(num), PARAMETER                          :: twothird=2.0_num/3.0_num
  REAL(num)                                     :: invvol, dxi, dyi, dzi
  REAL(num)                                     :: dts2dx, dts2dy, dts2dz
  INTEGER(isp)                                  :: j, k, l, j0, k0, l0, ip
  INTEGER(isp)                                  :: nnx, nnxy, ngridx, ngridy
  INTEGER(isp)                                  :: n, nn, nv
  INTEGER(isp), DIMENSION(lvect, 3)              :: ICELL
  REAL(num), DIMENSION(lvect)                   :: vx, vy, vz
  REAL(num)                                     :: wq
  REAL(num)                                     :: sx1(lvect), sx2(lvect),            &
  sx3(lvect), sx4(lvect)
  REAL(num) :: sx01(lvect), sx02(lvect), sx03(lvect), sx04(lvect)
  REAL(num) :: sy1, sy2, sy3, sy4, sz1, sz2, sz3, sz4
  REAL(num) :: sy01, sy02, sy03, sy04, sz01, sz02, sz03, sz04
  REAL(num), DIMENSION(4)                       :: h1, h11, h12, sgn
  REAL(num):: wwwx1(lvect, 8), wwwx2(lvect, 8), wwwy1(lvect, 8), wwwy2(lvect, 8),     &
  wwwz1(lvect, 8), wwwz2(lvect, 8)
  REAL(num)                                     :: wx1, wx2, wy1, wy2, wz1, wz2
  INTEGER(isp)                                  :: orig, ncxy, ngx, ngxy
  INTEGER(isp)                                  :: jorig, korig, lorig

  ! ___________________________________________________________
  ! Parameters
  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard

  IF (l_nodal) THEN
    stagger_shift = 0_num
  ELSE
    stagger_shift = 0.5_num
  ENDIF

  dxi = 1.0_num/dx
  dyi = 1.0_num/dy
  dzi = 1.0_num/dz
  invvol = dxi*dyi*dzi
  dts2dx = 0.5_num*dt*dxi
  dts2dy = 0.5_num*dt*dyi
  dts2dz = 0.5_num*dt*dzi

  nnx = ngridx
  nnxy = ngridx*ngridy
  jorig=-2_idp
  korig=-2_idp
  lorig=-2_idp
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  h1=(/1_num, 0_num, 1_num, 0_num/)
  sgn=(/1_num, -1_num, 1_num, -1_num/)
  h11=(/0_num, 1_num, 1_num, 0_num/)
  h12=(/1_num, 0_num, 0_num, 1_num/)

  ! LOOP ON PARTICLES

  DO ip=1, np, LVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED vx:64, vy:64, vz:64
    !DIR$ ASSUME_ALIGNED gaminv:64
    !DIR$ ASSUME_ALIGNED sx1:64, sx2:64, sx3:64, sx4:64
    !DIR$ ASSUME_ALIGNED sx01:64, sx02:64, sx03:64, sx04:64
    !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, xp, yp, zp)
    !IBM* ALIGN(64, vx, vy, vz)
    !IBM* ALIGN(64, sx1, sx2, sx3, sx4)
    !DIR$ ALIGN(64, sx01, sx02, sx03, sx04)
    !DIR$ ALIGN(64, ICELL)
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

      ! Computes velocity
      vx(n) = uxp(nn)*gaminv(nn)
      vy(n) = uyp(nn)*gaminv(nn)
      vz(n) = uzp(nn)*gaminv(nn)

      ! --- computes particles weights
      wq=q*w(nn)*invvol

      ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
      xmid=x-dts2dx*vx(n)
      ymid=y-dts2dy*vy(n)
      zmid=z-dts2dz*vz(n)

      ! --- finds node of cell containing particles for current positions
      j=floor(xmid)
      k=floor(ymid)
      l=floor(zmid)
      j0=floor(xmid-stagger_shift)
      k0=floor(ymid-stagger_shift)
      l0=floor(zmid-stagger_shift)
      ICELL(n, 1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
      ICELL(n, 3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

      ! --- computes set of coefficients for node centered quantities
      xint    = xmid-j
      yint    = ymid-k
      zint    = zmid-l
      oxint   = 1.0_num-xint
      xintsq  = xint*xint
      oxintsq = oxint*oxint
      sx1(n)  = onesixth*oxintsq*oxint
      sx2(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx3(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx4(n)  = onesixth*xintsq*xint
      oyint   = 1.0_num-yint
      yintsq  = yint*yint
      oyintsq = oyint*oyint
      sy1     = onesixth*oyintsq*oyint
      sy2     = (twothird-yintsq*(1.0_num-yint*0.5_num))
      sy3     = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
      sy4     = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz1 = onesixth*ozintsq*ozint*wq
      sz2 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
      sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
      sz4 = onesixth*zintsq*zint*wq

      ! --- computes set of coefficients for staggered quantities
      xint     = xmid-stagger_shift-j0
      yint     = ymid-stagger_shift-k0
      zint     = zmid-stagger_shift-l0
      oxint    = 1.0_num-xint
      xintsq   = xint*xint
      oxintsq  = oxint*oxint
      sx01(n)  = onesixth*oxintsq*oxint
      sx02(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx03(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx04(n)  = onesixth*xintsq*xint
      oyint    = 1.0_num-yint
      yintsq   = yint*yint
      oyintsq  = oyint*oyint
      sy01  = onesixth*oyintsq*oyint
      sy02  = (twothird-yintsq*(1.0_num-yint*0.5_num))
      sy03  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
      sy04  = onesixth*yintsq*yint
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz01 = onesixth*ozintsq*ozint*wq
      sz02 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
      sz03 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
      sz04 = onesixth*zintsq*zint*wq
      ! --- computes weights
      ! - X
      wwwx1(n, 1)=sz1*sy1
      wwwx1(n, 2)=sz2*sy1
      wwwx1(n, 3)=sz3*sy1
      wwwx1(n, 4)=sz4*sy1
      wwwx1(n, 5)=sz1*sy2
      wwwx1(n, 6)=sz2*sy2
      wwwx1(n, 7)=sz3*sy2
      wwwx1(n, 8)=sz4*sy2
      wwwx2(n, 1)=sz1*sy3
      wwwx2(n, 2)=sz2*sy3
      wwwx2(n, 3)=sz3*sy3
      wwwx2(n, 4)=sz4*sy3
      wwwx2(n, 5)=sz1*sy4
      wwwx2(n, 6)=sz2*sy4
      wwwx2(n, 7)=sz3*sy4
      wwwx2(n, 8)=sz4*sy4
      ! - Y
      wwwy1(n, 1)=sz1*sy01
      wwwy1(n, 2)=sz2*sy01
      wwwy1(n, 3)=sz3*sy01
      wwwy1(n, 4)=sz4*sy01
      wwwy1(n, 5)=sz1*sy02
      wwwy1(n, 6)=sz2*sy02
      wwwy1(n, 7)=sz3*sy02
      wwwy1(n, 8)=sz4*sy02
      wwwy2(n, 1)=sz1*sy03
      wwwy2(n, 2)=sz2*sy03
      wwwy2(n, 3)=sz3*sy03
      wwwy2(n, 4)=sz4*sy03
      wwwy2(n, 5)=sz1*sy04
      wwwy2(n, 6)=sz2*sy04
      wwwy2(n, 7)=sz3*sy04
      wwwy2(n, 8)=sz4*sy04
      ! - Z
      wwwz1(n, 1)=sz01*sy1
      wwwz1(n, 2)=sz02*sy1
      wwwz1(n, 3)=sz03*sy1
      wwwz1(n, 4)=sz04*sy1
      wwwz1(n, 5)=sz01*sy2
      wwwz1(n, 6)=sz02*sy2
      wwwz1(n, 7)=sz03*sy2
      wwwz1(n, 8)=sz04*sy2
      wwwz2(n, 1)=sz01*sy3
      wwwz2(n, 2)=sz02*sy3
      wwwz2(n, 3)=sz03*sy3
      wwwz2(n, 4)=sz04*sy3
      wwwz2(n, 5)=sz01*sy4
      wwwz2(n, 6)=sz02*sy4
      wwwz2(n, 7)=sz03*sy4
      wwwz2(n, 8)=sz04*sy4
    END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

    ! Add weights to nearest vertices

    DO n=1, MIN(LVEC, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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
        wx1=wwwx1(n, nv); wx2=wwwx2(n, nv)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx-1) = jxcells(nv, ICELL(n, 1)-ncx-1) +             &
        wx1*sx01(n)*vx(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1)-ncx)   = jxcells(nv, ICELL(n, 1)-ncx)   +             &
        wx1*sx02(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx+1) = jxcells(nv, ICELL(n, 1)-ncx+1) +             &
        wx1*sx03(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)-ncx+2) = jxcells(nv, ICELL(n, 1)-ncx+2) +             &
        wx1*sx04(n)*vx(n)
        ! Loop on (i=-1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx-1) = jxcells(nv, ICELL(n, 1)+ncx-1) +             &
        wx2*sx01(n)*vx(n)
        ! Loop on (i=0, j, k)
        jxcells(nv, ICELL(n, 1)+ncx)   = jxcells(nv, ICELL(n, 1)+ncx)   +             &
        wx2*sx02(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx+1) = jxcells(nv, ICELL(n, 1)+ncx+1) +             &
        wx2*sx03(n)*vx(n)
        !Loop on (i=1, j, k)
        jxcells(nv, ICELL(n, 1)+ncx+2) = jxcells(nv, ICELL(n, 1)+ncx+2) +             &
        wx2*sx04(n)*vx(n)

        ! --- JY
        wy1=wwwy1(n, nv); wy2=wwwy2(n, nv)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)-ncx-1) = jycells(nv, ICELL(n, 2)-ncx-1) +             &
        wy1*sx1(n)*vy(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2)-ncx)   = jycells(nv, ICELL(n, 2)-ncx)   +             &
        wy1*sx2(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)-ncx+1) = jycells(nv, ICELL(n, 2)-ncx+1) +             &
        wy1*sx3(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)-ncx+2) = jycells(nv, ICELL(n, 2)-ncx+2) +             &
        wy1*sx4(n)*vy(n)
        ! Loop on (i=-1, j, k)
        jycells(nv, ICELL(n, 2)+ncx-1) = jycells(nv, ICELL(n, 2)+ncx-1) +             &
        wy2*sx1(n)*vy(n)
        ! Loop on (i=0, j, k)
        jycells(nv, ICELL(n, 2)+ncx)   = jycells(nv, ICELL(n, 2)+ncx)   +             &
        wy2*sx2(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+ncx+1) = jycells(nv, ICELL(n, 2)+ncx+1) +             &
        wy2*sx3(n)*vy(n)
        !Loop on (i=1, j, k)
        jycells(nv, ICELL(n, 2)+ncx+2) = jycells(nv, ICELL(n, 2)+ncx+2) +             &
        wy2*sx4(n)*vy(n)

        ! --- JZ
        wz1=wwwz1(n, nv); wz2=wwwz2(n, nv)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx-1) = jzcells(nv, ICELL(n, 3)-ncx-1) +             &
        wz1*sx1(n)*vz(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3)-ncx)   = jzcells(nv, ICELL(n, 3)-ncx)   +             &
        wz1*sx2(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx+1) = jzcells(nv, ICELL(n, 3)-ncx+1) +             &
        wz1*sx3(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)-ncx+2) = jzcells(nv, ICELL(n, 3)-ncx+2) +             &
        wz1*sx4(n)*vz(n)
        ! Loop on (i=-1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx-1) = jzcells(nv, ICELL(n, 3)+ncx-1) +             &
        wz2*sx1(n)*vz(n)
        ! Loop on (i=0, j, k)
        jzcells(nv, ICELL(n, 3)+ncx)   = jzcells(nv, ICELL(n, 3)+ncx)   +             &
        wz2*sx2(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx+1) = jzcells(nv, ICELL(n, 3)+ncx+1) +             &
        wz2*sx3(n)*vz(n)
        !Loop on (i=1, j, k)
        jzcells(nv, ICELL(n, 3)+ncx+2) = jzcells(nv, ICELL(n, 3)+ncx+2) +             &
        wz2*sx4(n)*vz(n)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  RETURN
END SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3


! ________________________________________________________________________________________
!> @brief
!> This subroutine performs the reduction of jxcellx, jycells and jzcells
!> into jx, jy and jz.
!
!> @details
!> This subroutine is called after the loop on particles where
!> depose_jxjyjz_vecHV_vnr_1_1_1() is performed for each species
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2016
!
!> @param[inout] jx, jy, jz global current grids
!> @param[in] jxcells, jycells, jzcells temporary current arrays
!> @param[in] ncells tile cell numbers
!> @param[in] nx, ny, nz: tile cell numbers in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] ncx, ncy, ncz
! ________________________________________________________________________________________
SUBROUTINE current_reduction_1_1_1(jx, jy, jz, jxcells, jycells, jzcells, ncells, nx, &
  ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz)
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN)                 :: nx, ny, nz, ncells
  INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
  INTEGER(idp), INTENT(IN)                 :: nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT) ::                                                        &
  jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN), DIMENSION(8, ncells):: jxcells, jycells, jzcells
  INTEGER(isp)                             :: nnx, nnxy
  INTEGER(isp)                             :: moff(1:8)
  INTEGER(isp)                             :: orig, jorig, korig, lorig
  INTEGER(isp)                             :: igrid, ic
  INTEGER(isp) :: ncxy, ix, iy, iz, ngridx, ngridy, ngx, ngxy
  ! _____________________________________________________________
  ! Parameters

  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard

  nnx = ngridx
  nnxy = nnx*ngridy
  moff = (/0_isp, 1_isp, nnx, nnx+1_isp, nnxy, nnxy+1_isp, nnxy+nnx, nnxy+nnx+1_isp/)

  jorig=-2
  korig=-2
  lorig=-2
  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! ____________________________________________________________
  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
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
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy + orig
        ! jx
        jx(igrid+moff(1))=jx(igrid+moff(1))+jxcells(1, ic)
        jx(igrid+moff(2))=jx(igrid+moff(2))+jxcells(2, ic)
        jx(igrid+moff(3))=jx(igrid+moff(3))+jxcells(3, ic)
        jx(igrid+moff(4))=jx(igrid+moff(4))+jxcells(4, ic)
        jx(igrid+moff(5))=jx(igrid+moff(5))+jxcells(5, ic)
        jx(igrid+moff(6))=jx(igrid+moff(6))+jxcells(6, ic)
        jx(igrid+moff(7))=jx(igrid+moff(7))+jxcells(7, ic)
        jx(igrid+moff(8))=jx(igrid+moff(8))+jxcells(8, ic)
        ! jy
        jy(igrid+moff(1))=jy(igrid+moff(1))+jycells(1, ic)
        jy(igrid+moff(2))=jy(igrid+moff(2))+jycells(2, ic)
        jy(igrid+moff(3))=jy(igrid+moff(3))+jycells(3, ic)
        jy(igrid+moff(4))=jy(igrid+moff(4))+jycells(4, ic)
        jy(igrid+moff(5))=jy(igrid+moff(5))+jycells(5, ic)
        jy(igrid+moff(6))=jy(igrid+moff(6))+jycells(6, ic)
        jy(igrid+moff(7))=jy(igrid+moff(7))+jycells(7, ic)
        jy(igrid+moff(8))=jy(igrid+moff(8))+jycells(8, ic)
        ! jz
        jz(igrid+moff(1))=jz(igrid+moff(1))+jzcells(1, ic)
        jz(igrid+moff(2))=jz(igrid+moff(2))+jzcells(2, ic)
        jz(igrid+moff(3))=jz(igrid+moff(3))+jzcells(3, ic)
        jz(igrid+moff(4))=jz(igrid+moff(4))+jzcells(4, ic)
        jz(igrid+moff(5))=jz(igrid+moff(5))+jzcells(5, ic)
        jz(igrid+moff(6))=jz(igrid+moff(6))+jzcells(6, ic)
        jz(igrid+moff(7))=jz(igrid+moff(7))+jzcells(7, ic)
        jz(igrid+moff(8))=jz(igrid+moff(8))+jzcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO
  RETURN
END SUBROUTINE current_reduction_1_1_1

! ______________________________________________________________________________
!> @brief
!> This subroutine performs the reduction of jxcellx, jycells
!> and jzcells into jx, jy and jz.
!
!> @details
!> This subroutine is called after the loop on particles where
!> depose_jxjyjz_vecHV_vnr_2_2_2 is performed for each species
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[inout] jx, jy, jz global current grids
!> @param[in] jxcells, jycells, jzcells temporary current arrays
!> @param[in] ncells tile cell numbers
!> @param[in] nx, ny, nz: tile cell numbers in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] ncx, ncy, ncz
! ________________________________________________________________________________________
SUBROUTINE current_reduction_2_2_2(jx, jy, jz, jxcells, jycells, jzcells, ncells, nx, &
  ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz)
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN)               :: nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp), INTENT(IN)               :: ncx, ncy, ncz
  REAL(num), INTENT(IN OUT) ::                                                        &
  jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  INTEGER(idp), INTENT(IN)               :: ncells
  REAL(num), INTENT(IN), DIMENSION(8, ncells)  :: jxcells, jycells, jzcells
  INTEGER(isp)                           :: ic
  INTEGER(isp)                           :: moff(1:8)
  INTEGER(isp)                           :: igrid, orig, jorig, korig, lorig
  INTEGER(isp)                           :: ncxy, nnx, nnxy
  INTEGER(isp)                           :: ix, iy, iz, ngridx, ngridy, ngx, ngxy

  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard

  nnx = nx + 1 + 2*nxguard
  nnxy = nnx*(ny+1+2*nyguard)
  moff = (/-nnx-nnxy, -nnxy, nnx-nnxy, -nnx, nnx, -nnx+nnxy, nnxy, nnx+nnxy/)

  jorig=-2
  korig=-2
  lorig=-2

  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy

  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
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
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy + orig
        ! jx
        jx(igrid+moff(1))=jx(igrid+moff(1))+jxcells(1, ic)
        jx(igrid+moff(2))=jx(igrid+moff(2))+jxcells(2, ic)
        jx(igrid+moff(3))=jx(igrid+moff(3))+jxcells(3, ic)
        jx(igrid+moff(4))=jx(igrid+moff(4))+jxcells(4, ic)
        jx(igrid+moff(5))=jx(igrid+moff(5))+jxcells(5, ic)
        jx(igrid+moff(6))=jx(igrid+moff(6))+jxcells(6, ic)
        jx(igrid+moff(7))=jx(igrid+moff(7))+jxcells(7, ic)
        jx(igrid+moff(8))=jx(igrid+moff(8))+jxcells(8, ic)
        ! jy
        jy(igrid+moff(1))=jy(igrid+moff(1))+jycells(1, ic)
        jy(igrid+moff(2))=jy(igrid+moff(2))+jycells(2, ic)
        jy(igrid+moff(3))=jy(igrid+moff(3))+jycells(3, ic)
        jy(igrid+moff(4))=jy(igrid+moff(4))+jycells(4, ic)
        jy(igrid+moff(5))=jy(igrid+moff(5))+jycells(5, ic)
        jy(igrid+moff(6))=jy(igrid+moff(6))+jycells(6, ic)
        jy(igrid+moff(7))=jy(igrid+moff(7))+jycells(7, ic)
        jy(igrid+moff(8))=jy(igrid+moff(8))+jycells(8, ic)
        ! jz
        jz(igrid+moff(1))=jz(igrid+moff(1))+jzcells(1, ic)
        jz(igrid+moff(2))=jz(igrid+moff(2))+jzcells(2, ic)
        jz(igrid+moff(3))=jz(igrid+moff(3))+jzcells(3, ic)
        jz(igrid+moff(4))=jz(igrid+moff(4))+jzcells(4, ic)
        jz(igrid+moff(5))=jz(igrid+moff(5))+jzcells(5, ic)
        jz(igrid+moff(6))=jz(igrid+moff(6))+jzcells(6, ic)
        jz(igrid+moff(7))=jz(igrid+moff(7))+jzcells(7, ic)
        jz(igrid+moff(8))=jz(igrid+moff(8))+jzcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  RETURN
END SUBROUTINE current_reduction_2_2_2

! ________________________________________________________________________________________
!> @brief
!> This subroutine performs the reduction of jxcellx, jycells and jzcells
!> into jx, jy and jz.
!
!> @details
!> This subroutine is called after the loop on particles where
!> depose_jxjyjz_vecHV_vnr_3_3_3 is performed for each species
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[inout] jx, jy, jz global current grids
!> @param[in] jxcells, jycells, jzcells temporary current arrays
!> @param[in] ncells tile cell numbers
!> @param[in] nx, ny, nz: tile cell numbers in each direction
!> @param[in] nxguard, nyguard, nzguard number of guard cells
!> @param[in] ncx, ncy, ncz
! ________________________________________________________________________________________
SUBROUTINE current_reduction_3_3_3(jx, jy, jz, jxcells, jycells, jzcells, ncells, nx, &
  ny, nz, nxguard, nyguard, nzguard, ncx, ncy, ncz)
  USE picsar_precision, ONLY: idp, isp, num
  IMPLICIT NONE

  INTEGER(idp)                              :: nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp)                              :: ncx, ncy, ncz
  REAL(num), INTENT(IN OUT) ::                                                        &
  jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN OUT) ::                                                        &
  jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
  REAL(num), INTENT(IN), DIMENSION(8, ncells) :: jxcells, jycells, jzcells
  INTEGER(idp), INTENT(IN)                  :: ncells

  INTEGER(isp)                              :: moff(1:8)
  INTEGER(isp)                              :: ic
  INTEGER(isp)                              :: igrid, orig, jorig, korig, lorig
  INTEGER(isp)                              :: ncxy, nnx, nnxy
  INTEGER(isp)                              :: ix, iy, iz, ngridx, ngridy, ngx, ngxy

  ngridx=nx+1+2*nxguard
  ngridy=ny+1+2*nyguard

  nnx = ngridx
  nnxy = ngridx*ngridy

  moff = (/-nnxy, 0_isp, nnxy, 2_isp*nnxy, nnx-nnxy, nnx, nnx+nnxy, nnx+2_isp*nnxy/)

  jorig=-2_isp
  korig=-2_isp
  lorig=-2_isp

  orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy

  ngx=(ngridx-ncx)
  ngxy=(ngridx*ngridy-ncx*ncy)
  ncxy=ncx*ncy

  ! Reduction of jxcells, jycells, jzcells in jx, jy, jz
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
        igrid=ic+(iy-1)*ngx+(iz-1)*ngxy + orig
        ! jx
        jx(igrid+moff(1))=jx(igrid+moff(1))+jxcells(1, ic)
        jx(igrid+moff(2))=jx(igrid+moff(2))+jxcells(2, ic)
        jx(igrid+moff(3))=jx(igrid+moff(3))+jxcells(3, ic)
        jx(igrid+moff(4))=jx(igrid+moff(4))+jxcells(4, ic)
        jx(igrid+moff(5))=jx(igrid+moff(5))+jxcells(5, ic)
        jx(igrid+moff(6))=jx(igrid+moff(6))+jxcells(6, ic)
        jx(igrid+moff(7))=jx(igrid+moff(7))+jxcells(7, ic)
        jx(igrid+moff(8))=jx(igrid+moff(8))+jxcells(8, ic)
        ! jy
        jy(igrid+moff(1))=jy(igrid+moff(1))+jycells(1, ic)
        jy(igrid+moff(2))=jy(igrid+moff(2))+jycells(2, ic)
        jy(igrid+moff(3))=jy(igrid+moff(3))+jycells(3, ic)
        jy(igrid+moff(4))=jy(igrid+moff(4))+jycells(4, ic)
        jy(igrid+moff(5))=jy(igrid+moff(5))+jycells(5, ic)
        jy(igrid+moff(6))=jy(igrid+moff(6))+jycells(6, ic)
        jy(igrid+moff(7))=jy(igrid+moff(7))+jycells(7, ic)
        jy(igrid+moff(8))=jy(igrid+moff(8))+jycells(8, ic)
        ! jz
        jz(igrid+moff(1))=jz(igrid+moff(1))+jzcells(1, ic)
        jz(igrid+moff(2))=jz(igrid+moff(2))+jzcells(2, ic)
        jz(igrid+moff(3))=jz(igrid+moff(3))+jzcells(3, ic)
        jz(igrid+moff(4))=jz(igrid+moff(4))+jzcells(4, ic)
        jz(igrid+moff(5))=jz(igrid+moff(5))+jzcells(5, ic)
        jz(igrid+moff(6))=jz(igrid+moff(6))+jzcells(6, ic)
        jz(igrid+moff(7))=jz(igrid+moff(7))+jzcells(7, ic)
        jz(igrid+moff(8))=jz(igrid+moff(8))+jzcells(8, ic)
      END DO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
      !$OMP END SIMD
#endif
#endif
    END DO
  END DO

  RETURN
END SUBROUTINE current_reduction_3_3_3

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
! DIRECT_CURRENT_DEPOSITION_CIRC.F90
!
! Developers
! Jean-Luc Vay
!
! Description:
! This file contains subroutines for the direct current deposition in 3D.
!
! List of subroutines:
!
! Classical non-optimized:
! - depose_jrjtjl_scalar_1_1_1
! 
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Order 1 3D scalar direct current deposition routine (rho*v)
!> This version does not vectorize on SIMD architectures
!
!> @author
!> Jean-Luc Vay
!
!> @date
!> Creation 2019
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
!> @param[in] dr mesh size along x (scalar)
!> @param[in] dz mesh size along z (scalar)
!> @param[inout] jx x-current component (3D array)
!> @param[in] jr_nguard number of guard cells of the jr array in each direction
!> (1d array containing 3 integers)
!> @param[in] jr_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jr array (1d array containing 3 integers)
!> @param[inout] jt t-current component (3D array)
!> @param[in] jt_nguard number of guard cells of the jt array in each direction
!>  (1d array containing 3 integers)
!> @param[in] jt_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jt array (1d array containing 3 integers)
!> @param[inout] jl z-current component (3D array)
!> @param[in] jl_nguard number of guard cells of the jl array in each direction
!> (1d array containing 3 integers)
!> @param[in] jl_nvalid number of valid gridpoints (i.e. not guard cells) of the
!> jl array (1d array containing 3 integers)
!> @warning arrays jx, jy, jl should be set to 0 before entering this subroutine.
! ________________________________________________________________________________________
SUBROUTINE depose_jrjtjl_scalar_1_1_1( jr, jr_nguard, jr_nvalid, &
                                            jt, jt_nguard, jt_nvalid, &
                                            jl, jl_nguard, jl_nvalid, nmodes, &
                                            np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
                                            rmin, zmin, dt, dr, dz)     !#do not wrap
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num, cpx
  IMPLICIT NONE
  INTEGER(idp)             :: np, nmodes
  INTEGER(idp), intent(in) :: jr_nguard(2), jr_nvalid(2), &
                              jt_nguard(2), jt_nvalid(2), &
                              jl_nguard(2), jl_nvalid(2)
  COMPLEX(cpx), intent(IN OUT):: jr(-jr_nguard(1):jr_nvalid(1)+jr_nguard(1)-1, &
                                    -jr_nguard(2):jr_nvalid(2)+jr_nguard(2)-1, &
                                    0:nmodes-1 )
  COMPLEX(cpx), intent(IN OUT):: jt(-jt_nguard(1):jt_nvalid(1)+jt_nguard(1)-1, &
                                    -jt_nguard(2):jt_nvalid(2)+jt_nguard(2)-1, &
                                    0:nmodes-1 )
  COMPLEX(cpx), intent(IN OUT):: jl(-jl_nguard(1):jl_nvalid(1)+jl_nguard(1)-1, &
                                    -jl_nguard(2):jl_nvalid(2)+jl_nguard(2)-1, &
                                    0:nmodes-1 )
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num)                :: q, dt, dr, dz, rmin, zmin, c, s
  REAL(num)                :: dri, dzi, rint, zint
  REAL(num)                :: x, y, r, z, vx, vy, vr, vt, vz, invvol
  REAL(num)                :: wq, clightsq
  REAL(num), DIMENSION(2)  :: sr(0:1), sz(0:1)
  INTEGER(idp)             :: j, l, m, ip
  COMPLEX(cpx)             :: wqr, wqt, wql, exptheta_m

  dri = 1.0_num/dr
  dzi = 1.0_num/dz
  invvol = dri*dzi
  clightsq = 1.0_num/clight**2
  sr=0.0_num;sz=0.0_num;

  ! LOOP ON PARTICLES
  ! Prevent loop to vectorize (dependencies)
  !DIR$ NOVECTOR
  DO ip=1, np

    ! Computes velocity (Cartesian)
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    
    ! --- computes position in  grid units at (n+1/2)
    x = xp(ip)-0.5*vx*dt
    y = yp(ip)-0.5*vy*dt
    r=sqrt(x*x+y*y)
    if (r*dri>1.e-10) then
      c = x/r 
      s = y/r
    else
      c = 1.
      s = 0.
    end if
    exptheta_m = cmplx(c,s)
    r = (r-rmin)*dri-0.5 ! grid starts at 0.5*dx in x

    z = (zp(ip)-0.5*vz*dt-zmin)*dzi

    ! Computes velocity (radial and azimuthal)
    vr = c*vx + s*vz
    vt = c*vy - s*vx

    ! --- computes particles weights
    wq=q*w(ip)
    wqr=wq*invvol*vr
    wqt=wq*invvol*vt
    wql=wq*invvol*vz

    ! --- finds node of cell containing particles for current positions
    j=floor(r)
    l=floor(z)

    ! --- computes set of coefficients for node centered quantities
    rint = r-j
    zint = z-l
    sr( 0) = 1.0_num-rint
    sr( 1) = rint
    sz( 0) = 1.0_num-zint
    sz( 1) = zint

    DO m = 0, nmodes-1
        ! --- add current contributions in the form rho(n+1/2)v(n+1/2) to modes m=0...circ_m
        ! - jr
        jr(j,   l,   m) = jr(j,   l,   m  )  +   sr(0)*sz(0)*wqr
        jr(j+1, l,   m) = jr(j+1, l,   m  )  +   sr(1)*sz(0)*wqr
        jr(j,   l+1, m) = jr(j,   l+1, m  )  +   sr(0)*sz(1)*wqr
        jr(j+1, l+1, m) = jr(j+1, l+1, m  )  +   sr(1)*sz(1)*wqr

        ! - jt
        jt(j,   l,   m) = jt(j,   l,   m  )  +   sr(0)*sz(0)*wqt
        jt(j+1, l,   m) = jt(j+1, l,   m  )  +   sr(1)*sz(0)*wqt
        jt(j,   l+1, m) = jt(j,   l+1, m  )  +   sr(0)*sz(1)*wqt
        jt(j+1, l+1, m) = jt(j+1, l+1, m  )  +   sr(1)*sz(1)*wqt
        
        ! - jl
        jl(j,   l,   m) = jl(j,   l,   m  )  +   sr(0)*sz(0)*wql
        jl(j+1, l,   m) = jl(j+1, l,   m  )  +   sr(1)*sz(0)*wql
        jl(j,   l+1, m) = jl(j,   l+1, m  )  +   sr(0)*sz(1)*wql
        jl(j+1, l+1, m) = jl(j+1, l+1, m  )  +   sr(1)*sz(1)*wql

        wqr = wqr*exptheta_m
        wqt = wqt*exptheta_m
        wql = wql*exptheta_m
        if (m==0) then 
            wqr = 2*wqr
            wqt = 2*wqt
            wql = 2*wql
        END IF
    END DO

  END DO
  RETURN
END SUBROUTINE depose_jrjtjl_scalar_1_1_1


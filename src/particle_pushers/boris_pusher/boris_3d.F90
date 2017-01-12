! ______________________________________________________________________________
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
! BORIS_3D.F90
!
! Subroutines for the Boris particle pushers in 3d.
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Advance particle positions
!
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] ux,uy,uz normalized momentum in each direction
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] dt time step
SUBROUTINE pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
! ______________________________________________________________________________
  USE constants
  USE omp_lib
  
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num)      :: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np), gaminv(np)
  REAL(num)      :: dt
  INTEGER(idp)   :: ip

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,xp,yp,zp)
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,gaminv)
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
  DO ip=1,np
    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

  RETURN
END SUBROUTINE pxr_pushxyz


! ______________________________________________________________________________
!> @brief
!> Push the particle velocity with E field
!
!> @details
!>fast b-field rotation algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex,ey,ez particle electric fields in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
SUBROUTINE pxr_epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)
! ______________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp) :: np
  REAL(num)    :: uxp(np),uyp(np),uzp(np)
  REAL(num)    :: ex(np),ey(np),ez(np)
  REAL(num)    :: q,m,dt
  INTEGER(idp) :: ip
  REAL(num)    :: const

  const = q*dt/m

#if defined __INTEL_COMPILER
      !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
      !DIR$ ASSUME_ALIGNED ex:64,ey:64,ez:64
#elif defined __IBMBGQ__
      !IBM* ALIGN(64,uxp,uyp,uzp)
      !IBM* ALIGN(64,ex,ey,ez)
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
  DO ip=1,np
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif
  RETURN
END SUBROUTINE pxr_epush_v


! ______________________________________________________________________________
!> @brief
!> Push the particle velocity with B field (Boris algorithm)
!
!> @details
!>fast b-field rotation algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @gparam[in] gaminv particle Lorentz factors
!> @param[in] bx,by,bz particle magnetic fields in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
SUBROUTINE pxr_bpush_v(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,dt)
! ______________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num)      :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num)      :: bx(np), by(np), bz(np)
  REAL(num)      :: q,m,dt
  INTEGER(idp)   :: ip
  REAL(num)      :: const,sx,sy,sz,tx,ty,tz,tsqi,uxppr,uyppr,uzppr

  const = q*dt*0.5_num/m

#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
    !DIR$ ASSUME_ALIGNED bx:64,by:64,bz:64
    !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64,uxp,uyp,uzp)
    !IBM* ALIGN(64,bx,by,bz)
    !IBM* ALIGN(64,gaminv)
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
  DO ip=1,np
    tx = gaminv(ip)*bx(ip)*const
    ty = gaminv(ip)*by(ip)*const
    tz = gaminv(ip)*bz(ip)*const
    tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
    sx = tx*tsqi
    sy = ty*tsqi
    sz = tz*tsqi
    uxppr = uxp(ip) + uyp(ip)*tz - uzp(ip)*ty
    uyppr = uyp(ip) + uzp(ip)*tx - uxp(ip)*tz
    uzppr = uzp(ip) + uxp(ip)*ty - uyp(ip)*tx
    uxp(ip) = uxp(ip) + uyppr*sz - uzppr*sy
    uyp(ip) = uyp(ip) + uzppr*sx - uxppr*sz
    uzp(ip) = uzp(ip) + uxppr*sy - uyppr*sx
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

  RETURN

END SUBROUTINE pxr_bpush_v


! ______________________________________________________________________________
!> @brief
!>  Push the particle velocity with B field (Boris algorithm)
!
!> @details
!> fast b-field rotation algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 06.10.2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp,uyp,uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!
SUBROUTINE pxr_set_gamma(np,uxp,uyp,uzp,gaminv)
! ______________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp)   :: ip, np
  REAL(num) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num) :: clghtisq, usq

  clghtisq = 1.0_num/clight**2

#if defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
    !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64,uxp,uyp,uzp)
    !IBM* ALIGN(64,gaminv)
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

  DO ip=1,np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
  END DO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

RETURN

END SUBROUTINE pxr_set_gamma



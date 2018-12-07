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
! BORIS_2D.F90
!
! Subroutines for the Boris particle pushers in 2d.
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Advance particle positions 2D Case, serial version.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
!
!> @param[in] np number of particles
!> @param[inout] xp, zp particle position arrays
!> @param[in] uxp, uzp particle momentum arrays
!> @param[in] gaminv particle inverse Lorentz factor
!> @param[in] dt time step
! ________________________________________________________________________________________
SUBROUTINE pxr_pushxz(np, xp, zp, uxp, uzp, gaminv, dt)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num) :: xp(np), zp(np), uxp(np), uzp(np), gaminv(np)
  REAL(num) :: dt
  INTEGER(idp)  :: ip

  !$acc parallel deviceptr(xp, zp, uxp, uzp, gaminv)
  !$acc loop gang vector 
  DO ip=1, np
    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt
  ENDDO
  !$acc end loop
  !$acc end parallel

  RETURN
END SUBROUTINE pxr_pushxz


! ________________________________________________________________________________________
!> @brief
!> Advance particle positions 2D Case, vectorized version.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
!
!> @param[in] np number of particles
!> @param[inout] xp, zp particle position arrays
!> @param[in] uxp, uzp particle momentum arrays
!> @param[in] gaminv particle inverse Lorentz factor
!> @param[in] dt time step
! ________________________________________________________________________________________
SUBROUTINE pxr_push2dxz(np, xp, zp, uxp, uyp, uzp, gaminv, dt)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp)   :: np
  REAL(num) :: xp(np), zp(np), uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num) :: dt
  INTEGER(idp)  :: ip


#if defined _OPENMP && _OPENMP>=201307
  !$OMP SIMD
#elif defined __IBMBGQ__
  !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
  !$DIR SIMD
#endif
  !$acc parallel deviceptr(xp, zp, uxp, uzp, gaminv)
  !$acc loop gang vector 
  DO ip=1, np
    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt
 ENDDO
 !$acc end loop
 !$acc end parallel
#if defined _OPENMP && _OPENMP>=201307
  !$OMP END SIMD
#endif

  RETURN
END SUBROUTINE pxr_push2dxz

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
! BORIS_3D.F90
!
! Subroutines for the Boris particle pushers in 3d.
!
! Developers:
! Henri Vincenti
! Mathieu Lobet
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the momentum to the next step using the Boris pusher.
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!
!> @date
!> Creation January 16 2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex, ey, ez particle electric field arrays
!> @param[in] bx, by, bz particle magnetic field arrays
!> @param[in] dt time step
!
! ________________________________________________________________________________________
SUBROUTINE pxr_boris_push_u_3d(np, uxp, uyp, uzp, gaminv, ex, ey, ez, bx, by, bz, q,  &
  m, dt)
  USE constants
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: np
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np)
  REAL(num), INTENT(IN)    :: bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)    :: q, m, dt
  ! Local variables
  INTEGER(idp)             :: ip
  REAL(num)                :: const
  REAL(num)                :: clghtisq, usq, tsqi
  REAL(num)                :: tx, ty, tz
  REAL(num)                :: sx, sy, sz
  REAL(num)                :: uxppr, uyppr, uzppr
  REAL(num)                :: gaminvtmp
  ! Initialization
  const = q*dt*0.5_num/m
  clghtisq = 1.0_num/clight**2

  ! Loop over the particles
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
  !DIR$ ASSUME_ALIGNED gaminv:64
  !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
  !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, uxp, uyp, uzp)
  !IBM* ALIGN(64, gaminv)
  !IBM* ALIGN(64, ex, ey, ez)
  !IBM* ALIGN(64, bx, by, bz)
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
  DO ip=1, np
    ! Push using the electric field
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const

    ! Compute temporary Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminvtmp = 1.0_num/sqrt(1.0_num + usq)

    ! Magnetic rotation
    tx = gaminvtmp*bx(ip)*const
    ty = gaminvtmp*by(ip)*const
    tz = gaminvtmp*bz(ip)*const
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

    ! Push using the electric field
    uxp(ip) = uxp(ip) + ex(ip)*const
    uyp(ip) = uyp(ip) + ey(ip)*const
    uzp(ip) = uzp(ip) + ez(ip)*const

    ! Compute final Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)

  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE


! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the momentum to the next step using the Boris pusher.
!> The particle loop is divided by block of particles for
!> cache-blocking purposes.
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!
!> @date
!> Creation January 16 2016
!
!> @param[in] np number of super-particles
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex, ey, ez particle electric field arrays
!> @param[in] bx, by, bz particle magnetic field arrays
!> @param[in] dt time step
!
! ________________________________________________________________________________________
SUBROUTINE pxr_boris_push_u_3d_block(np, uxp, uyp, uzp, gaminv, ex, ey, ez, bx, by,   &
  bz, q, m, dt, lvect)
  USE constants
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: np
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np)
  REAL(num), INTENT(IN)    :: bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)    :: q, m, dt
  INTEGER(idp), INTENT(IN) :: lvect
  ! Local variables
  INTEGER(idp)             :: ip, nn, n, blocksize
  REAL(num)                :: const
  REAL(num)                :: clghtisq, usq, tsqi
  REAL(num)                :: tx, ty, tz
  REAL(num)                :: sx, sy, sz
  REAL(num)                :: uxppr, uyppr, uzppr
  REAL(num)                :: gaminvtmp

  ! Initialization
  const = q*dt*0.5_num/m
  clghtisq = 1.0_num/clight**2

  ! Loop over the particles
  DO ip=1, np, lvect

    ! Size of the block
    blocksize = MIN(lvect, np-ip+1)

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED gaminv:64
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, gaminv)
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
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

    ! Loop on the particles of blocks
    DO n=1, blocksize
      nn=ip+n-1

      ! Push using the electric field
      uxp(nn) = uxp(nn) + ex(nn)*const
      uyp(nn) = uyp(nn) + ey(nn)*const
      uzp(nn) = uzp(nn) + ez(nn)*const

      ! Compute temporary Gamma
      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminvtmp = 1.0_num/sqrt(1.0_num + usq)

      ! Magnetic rotation
      tx = gaminvtmp*bx(nn)*const
      ty = gaminvtmp*by(nn)*const
      tz = gaminvtmp*bz(nn)*const
      tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
      sx = tx*tsqi
      sy = ty*tsqi
      sz = tz*tsqi
      uxppr = uxp(nn) + uyp(nn)*tz - uzp(nn)*ty
      uyppr = uyp(nn) + uzp(nn)*tx - uxp(nn)*tz
      uzppr = uzp(nn) + uxp(nn)*ty - uyp(nn)*tx
      uxp(nn) = uxp(nn) + uyppr*sz - uzppr*sy
      uyp(nn) = uyp(nn) + uzppr*sx - uxppr*sz
      uzp(nn) = uzp(nn) + uxppr*sy - uyppr*sx

      ! Push using the electric field
      uxp(nn) = uxp(nn) + ex(nn)*const
      uyp(nn) = uyp(nn) + ey(nn)*const
      uzp(nn) = uzp(nn) + ez(nn)*const

      ! Compute final Gamma
      usq = (uxp(nn)**2 + uyp(nn)**2+ uzp(nn)**2)*clghtisq
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
  ENDDO
END SUBROUTINE

SUBROUTINE pxr_boris_push_u_3d_block_rr


END SUBROUTINE 
! ________________________________________________________________________________________
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
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] dt time step
! ________________________________________________________________________________________
SUBROUTINE pxr_pushxyz(np, xp, yp, zp, uxp, uyp, uzp, gaminv, dt)
  USE constants
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN)   :: np
  REAL(num), INTENT(INOUT)   :: xp(np), yp(np), zp(np)
  REAL(num), INTENT(IN)      :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)      :: dt
  ! Local parameters
  INTEGER(idp)               :: ip

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
  !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
  !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, xp, yp, zp)
  !IBM* ALIGN(64, uxp, uyp, uzp)
  !IBM* ALIGN(64, gaminv)
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
  DO ip=1, np
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

! ________________________________________________________________________________________
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
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex, ey, ez particle electric fields in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
! ________________________________________________________________________________________
SUBROUTINE pxr_epush_v(np, uxp, uyp, uzp, ex, ey, ez, q, m, dt)
  USE constants
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: np
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np)
  REAL(num), INTENT(IN)    :: q, m, dt
  ! Local parameters
  INTEGER(idp) :: ip
  REAL(num)    :: const

  const = q*dt/m

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
  !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, uxp, uyp, uzp)
  !IBM* ALIGN(64, ex, ey, ez)
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
  DO ip=1, np
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


! ________________________________________________________________________________________
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
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @gparam[in] gaminv particle Lorentz factors
!> @param[in] bx, by, bz particle magnetic fields in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
! ________________________________________________________________________________________
SUBROUTINE pxr_bpush_v(np, uxp, uyp, uzp, gaminv, bx, by, bz, q, m, dt)
  USE constants
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN)   :: np
  REAL(num), INTENT(INOUT)   :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)      :: bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)      :: q, m, dt
  ! Local parameters
  INTEGER(idp)   :: ip
  REAL(num)      :: const, sx, sy, sz, tx, ty, tz, tsqi, uxppr, uyppr, uzppr

  const = q*dt*0.5_num/m

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
  !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
  !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, uxp, uyp, uzp)
  !IBM* ALIGN(64, bx, by, bz)
  !IBM* ALIGN(64, gaminv)
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
  DO ip=1, np
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


! ________________________________________________________________________________________
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
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!
! ________________________________________________________________________________________
SUBROUTINE pxr_set_gamma(np, uxp, uyp, uzp, gaminv)
  USE constants
  IMPLICIT NONE

  ! Input/output parameters
  INTEGER(idp), INTENT(IN)   :: np
  REAL(num), INTENT(IN)      :: uxp(np), uyp(np), uzp(np)
  REAL(num), INTENT(INOUT)   :: gaminv(np)

  ! Local parameters
  INTEGER(idp)   :: ip
  REAL(num)      :: clghtisq, usq

  clghtisq = 1.0_num/clight**2

#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
  !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, uxp, uyp, uzp)
  !IBM* ALIGN(64, gaminv)
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

  DO ip=1, np
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

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
! VAY_3D.F90
!
! Subroutines for the J.L. Vay particle pusher in 3D
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Push the particle velocity with E and B fields, assuming Vmid = 0.5*(Vold+Vnew),
!> solving directly for the new gamma.
!>
!> @details
!> This offers better cancellation of E+VxB than the Boris velocity push.
!> Question: should we recompute gamma from the new u, in order to prevent
!> roundoff errors
!> to create mismatched values of u and gamma?
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
!
!>
!> @param[in] np number of super-particles
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gi
!> @param[in] ex, ey, ez particle electric field values in each direction
!> @param[in] bx, by, bz particle electric field values in each direction
!> @param[in] q charge
!> @param[in] m masse
!> @param[in] dt time step
!> @param[in] which algorithm
!
! ________________________________________________________________________________________
SUBROUTINE pxr_ebcancelpush3d(np, uxp, uyp, uzp, gi, ex, ey, ez, bx, by, bz, q, &
  m, dt, which)
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  ! Input/Ooutput parameters
  INTEGER(idp), INTENT(IN) :: np, which
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np), gi(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np), bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)    :: q, m, dt

  ! Local parameters
  INTEGER(idp) :: ip
  REAL(num)    :: const, bconst, s, gisq, invclight, invclightsq, gprsq
  REAL(num)    :: tx, ty, tz, tu, uxpr, uypr, uzpr, bg, vx, vy, vz
  REAL(num)    :: taux, tauy, tauz, tausq, ust, sigma

  invclight   = 1./clight
  invclightsq = 1./(clight*clight)

  IF (which==0) THEN
    !     --- full push
    const = q*dt/m
    bconst = 0.5_num*const

#ifndef NOVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
    !DIR$ ASSUME_ALIGNED gi:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
    !IBM* ALIGN(64, gi)
#endif
#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#endif
!$acc parallel deviceptr(uxp, uyp, uzp, gi, ex, ey, ez, bx, by, bz)
!$acc loop gang vector
    DO ip=1, np
      ! --- get tau
      taux = bconst*bx(ip)
      tauy = bconst*by(ip)
      tauz = bconst*bz(ip)
      tausq = taux*taux+tauy*tauy+tauz*tauz
      ! --- get U', gamma'^2
      uxpr = uxp(ip) + const*ex(ip) + (uyp(ip)*tauz-uzp(ip)*tauy)*gi(ip)
      uypr = uyp(ip) + const*ey(ip) + (uzp(ip)*taux-uxp(ip)*tauz)*gi(ip)
      uzpr = uzp(ip) + const*ez(ip) + (uxp(ip)*tauy-uyp(ip)*taux)*gi(ip)
      gprsq = (1._num+(uxpr*uxpr+uypr*uypr+uzpr*uzpr)*invclightsq)
      !       --- get u*
      ust = (uxpr*taux+uypr*tauy+uzpr*tauz)*invclight
      ! --- get new gamma
      sigma = gprsq-tausq
      gisq = 2._num/(sigma+sqrt(sigma*sigma+4._num*(tausq+ust*ust)))
      gi(ip) = sqrt(gisq)
      !       --- get t, s
      bg = bconst*gi(ip)
      tx = bg*bx(ip)
      ty = bg*by(ip)
      tz = bg*bz(ip)
      s = 1._num/(1._num+tausq*gisq)
      !  --- get t.u'
      tu = tx*uxpr+ty*uypr+tz*uzpr
      ! --- get new U
      uxp(ip) = s*(uxpr+tx*tu+uypr*tz-uzpr*ty)
      uyp(ip) = s*(uypr+ty*tu+uzpr*tx-uxpr*tz)
      uzp(ip) = s*(uzpr+tz*tu+uxpr*ty-uypr*tx)
    END DO
!$acc end loop
!$acc end parallel

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif

  ELSE IF(which==1) THEN

    !     --- first half push
    const = 0.5_num*q*dt/m

#ifndef NOVEC
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
    !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
    !DIR$ ASSUME_ALIGNED gi:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, ex, ey, ez)
    !IBM* ALIGN(64, bx, by, bz)
    !IBM* ALIGN(64, gi)
#endif
#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#endif
!$acc parallel deviceptr(uxp, uyp, uzp, gi, ex, ey, ez, bx, by, bz)
!$acc loop gang vector
    DO ip=1, np
      ! --- get new U
      vx = uxp(ip)*gi(ip)
      vy = uyp(ip)*gi(ip)
      vz = uzp(ip)*gi(ip)
      uxp(ip) = uxp(ip) + const*( ex(ip) + vy*bz(ip)-vz*by(ip) )
      uyp(ip) = uyp(ip) + const*( ey(ip) + vz*bx(ip)-vx*bz(ip) )
      uzp(ip) = uzp(ip) + const*( ez(ip) + vx*by(ip)-vy*bx(ip) )
      gi(ip) =                                                                        &
      1./sqrt(1.+(uxp(ip)*uxp(ip)+uyp(ip)*uyp(ip)+uzp(ip)*uzp(ip))*invclightsq)
    ENDDO
!$acc end loop
!$acc end parallel
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif


  ELSE IF(which==2) THEN
    !     --- second half push
    const = 0.5_num*q*dt/m
    bconst = const

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
      !     --- get U'
      uxpr = uxp(ip) + const*ex(ip)
      uypr = uyp(ip) + const*ey(ip)
      uzpr = uzp(ip) + const*ez(ip)
      gprsq = (1_num+(uxpr*uxpr+uypr*uypr+uzpr*uzpr)*invclightsq)
      !       --- get tau
      taux = bconst*bx(ip)
      tauy = bconst*by(ip)
      tauz = bconst*bz(ip)
      tausq = taux*taux+tauy*tauy+tauz*tauz
      !       --- get u*
      ust = (uxpr*taux+uypr*tauy+uzpr*tauz)*invclight
      !       --- get new gamma
      sigma = gprsq-tausq
      gisq = 2._num/(sigma+sqrt(sigma*sigma+4._num*(tausq+ust*ust)))
      gi(ip) = sqrt(gisq)
      !       --- get t, s
      bg = bconst*gi(ip)
      tx = bg*bx(ip)
      ty = bg*by(ip)
      tz = bg*bz(ip)
      s = 1._num/(1._num+tausq*gisq)
      !       --- get t.u'
      tu = tx*uxpr+ty*uypr+tz*uzpr
      !       --- get new U
      uxp(ip) = s*(uxpr+tx*tu+uypr*tz-uzpr*ty)
      uyp(ip) = s*(uypr+ty*tu+uzpr*tx-uxpr*tz)
      uzp(ip) = s*(uzpr+tz*tu+uxpr*ty-uypr*tx)
    ENDDO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !$OMP END SIMD
#endif
#endif
  ENDIF
  RETURN
END SUBROUTINE pxr_ebcancelpush3d

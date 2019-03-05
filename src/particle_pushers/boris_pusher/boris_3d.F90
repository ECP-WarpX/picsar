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
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
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
  
!$acc parallel deviceptr(uxp, uyp, uzp, gaminv, ex, ey, ez, bx, by, bz)
!$acc loop gang vector
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
!$acc end loop
!$acc end parallel
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE pxr_boris_push_u_3d


! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the momentum to the next step using the Boris pusher taking into account the classical Radiation Reaction (RR) force (S09 model).
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!> Murad Abuzarli
!
!> @date
!> Creation April 16 2018
!
!> @param[in] np number of super-particles
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex, ey, ez particle electric field arrays
!> @param[in] bx, by, bz particle magnetic field arrays
!> @param[in] dt time step
!
! ________________________________________________________________________________________
SUBROUTINE pxr_boris_push_rr_S09_u_3d(np, uxp, uyp, uzp, gaminv, ex, ey, ez, bx, by, bz, q,  &
  m, dt)
  USE constants, ONLY: clight, mu0, pi
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: np
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np)
  REAL(num), INTENT(IN)    :: bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)    :: q, m, dt
  ! Local variables
  INTEGER(idp)             :: ip
  REAL(num)                :: const, dtinv
  REAL(num)                :: clghtisq, usq, tsqi
  REAL(num)                :: tx, ty, tz
  REAL(num)                :: sx, sy, sz
  REAL(num)                :: uxppr, uyppr, uzppr
  REAL(num)                :: gaminvtmp, gamtmp
  REAL(num)                :: uxold, uyold, uzold
  REAL(num)                :: alx, aly, alz
  REAL(num)                :: urx, ury, urz
  REAL(num)                :: dvx, dvy, dvz
  REAL(num)                :: ddv
  REAL(num)                :: dval
  REAL(num)                :: tau0inv
  REAL(num)                :: ual

  ! Initialization
  const = q*dt*0.5_num/m
  dtinv=1_num/dt
  clghtisq = 1.0_num/clight**2
  ! tau zero
  tau0inv=6_num*pi*m*clight/(mu0*q**2)

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

    ! --- Pushing with the Lorentz force

    ! Save previous moments
    uxold=uxp(ip)
    uyold=uyp(ip)
    uzold=uzp(ip)

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

    ! --- Pushing with the RR force
    
    ! Lorentz acceleration
    alx=(uxp(ip)-uxold)*dtinv
    aly=(uyp(ip)-uyold)*dtinv
    alz=(uzp(ip)-uzold)*dtinv

    ! Momentum/mass at half step
    uxp(ip) = (uxp(ip) + uxold)*0.5_num
    uyp(ip) = (uyp(ip) + uyold)*0.5_num
    uzp(ip) = (uzp(ip) + uzold)*0.5_num

    ! Compute temporary Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gamtmp=sqrt(1.0_num + usq)
    gaminvtmp = 1.0_num/gamtmp

    ! u dot al
    ual=uxp(ip)*alx+uyp(ip)*aly+uzp(ip)*alz
    
    ! denominator of delta u
    ddv=1.0_num/(tau0inv+clghtisq*gaminvtmp*ual)

    ! delta v
    dvx=ddv*(alx-clghtisq*gaminvtmp**2*ual*uxp(ip))
    dvy=ddv*(aly-clghtisq*gaminvtmp**2*ual*uyp(ip))
    dvz=ddv*(alz-clghtisq*gaminvtmp**2*ual*uzp(ip))

    ! delta v dot al
    dval=dvx*alx+dvy*aly+dvz*alz

    ! RR correction
    urx=2.0_num*const*(dvy*bz(ip)-dvz*by(ip))-clghtisq*gamtmp*dval*uxp(ip)*dt
    ury=2.0_num*const*(dvz*bx(ip)-dvx*bz(ip))-clghtisq*gamtmp*dval*uyp(ip)*dt
    urz=2.0_num*const*(dvx*by(ip)-dvy*bx(ip))-clghtisq*gamtmp*dval*uzp(ip)*dt

    ! Push using the RR force
    uxp(ip) = uxold+urx+dt*alx
    uyp(ip) = uyold+ury+dt*aly
    uzp(ip) = uzold+urz+dt*alz

    ! Compute final Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE pxr_boris_push_rr_S09_u_3d


! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the momentum to the next step using the Boris pusher taking into account the classical Radiation Reaction (RR) force (B08 model).
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!> Murad Abuzarli
!
!> @date
!> Creation April 19 2018
!
!> @param[in] np number of super-particles
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] ex, ey, ez particle electric field arrays
!> @param[in] bx, by, bz particle magnetic field arrays
!> @param[in] dt time step
!
! ________________________________________________________________________________________
SUBROUTINE pxr_boris_push_rr_B08_u_3d(np, uxp, uyp, uzp, gaminv, ex, ey, ez, bx, by, bz, q,  &
  m, dt)
  USE constants, ONLY: clight, eps0, pi
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: np
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np)
  REAL(num), INTENT(IN)    :: bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)    :: q, m, dt
  ! Local variables
  INTEGER(idp)             :: ip
  REAL(num)                :: const, dtinv, crr
  REAL(num)                :: clghtisq, usq, tsqi
  REAL(num)                :: tx, ty, tz
  REAL(num)                :: sx, sy, sz
  REAL(num)                :: uxppr, uyppr, uzppr
  REAL(num)                :: gaminvtmp, gamtmp
  REAL(num)                :: uxold, uyold, uzold
  REAL(num)                :: etx, ety, etz
  REAL(num)                :: urx, ury, urz
  REAL(num)                :: ulx, uly, ulz
  REAL(num)                :: epar
  REAL(num)                :: ddv
  REAL(num)                :: dval
  REAL(num)                :: minv, upinvsq
  REAL(num)                :: fsq

  ! Initialization
  const = q*dt*0.5_num/m
  dtinv=1_num/dt
  clghtisq = 1.0_num/clight**2
  
  ! inverse of mass
  minv=1_num/m

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

    ! --- Pushing with the Lorentz force

    ! Save previous moments
    uxold=uxp(ip)
    uyold=uyp(ip)
    uzold=uzp(ip)

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
!write(*,*)"check 1",uxp(ip),uyp(ip),uzp(ip)
    ! Compute final Gamma
   usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
   gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)

    ! --- Pushing with the RR force

    ! Lorentz momentum
    ulx = uxp(ip)
    uly = uyp(ip)
    ulz = uzp(ip)

    ! Momentum/mass at half step
    uxp(ip) = (uxp(ip) + uxold)*0.5_num
    uyp(ip) = (uyp(ip) + uyold)*0.5_num
    uzp(ip) = (uzp(ip) + uzold)*0.5_num

    ! Compute temporary Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gamtmp=sqrt(1.0_num + usq)
    gaminvtmp = 1.0_num/gamtmp

    ! RR force constant
    crr=-(4_num*const*q**3*gamtmp)/(3_num*m**2*clight**5*4*pi*eps0)
    
    ! 'up' inverse square norm 
    upinvsq=1_num/(uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)

    ! E dot u / u square
    epar=upinvsq*(ex(ip)*uxp(ip)+ey(ip)*uyp(ip)+ez(ip)*uzp(ip))

    ! Transverse E field 
    etx=(ex(ip)-epar*uxp(ip))
    ety=(ey(ip)-epar*uyp(ip))
    etz=(ez(ip)-epar*uzp(ip))

    ! Square of fields
    fsq=crr*((etx+gaminvtmp*(uyp(ip)*bz(ip)-uzp(ip)*by(ip)))**2+ & 
    (ety+gaminvtmp*(uzp(ip)*bx(ip)-uxp(ip)*bz(ip)))**2+          & 
    (etz+gaminvtmp*(uxp(ip)*by(ip)-uyp(ip)*bx(ip)))**2)


    ! RR momentum
    urx=fsq*uxp(ip)
    ury=fsq*uyp(ip)
    urz=fsq*uzp(ip)

    ! Push using the RR force
    uxp(ip) = urx+ulx
    uyp(ip) = ury+uly
    uzp(ip) = urz+ulz

    ! Compute final Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE pxr_boris_push_rr_B08_u_3d


! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the momentum to the next step using the Boris pusher taking 
!> into account the classical Radiation Reaction (RR) force (LL model).
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!> Murad Abuzarli
!
!> @date
!> Creation April 27 2018
!
!> @param[in] np number of super-particles
!> @param[in] uxp, uyp, uzp normalized momentum in each direction
!> @param[in] gaminv particle Lorentz factors
!> @param[in] exold, eyold, ezold particle electric field arrays at previous iteration
!> @param[in] bxold, byold, bzold particle magnetic field arrays at previous iteration
!> @param[in] ex, ey, ez particle electric field arrays
!> @param[in] bx, by, bz particle magnetic field arrays
!> @param[in] dt time step
!
! ________________________________________________________________________________________
SUBROUTINE pxr_boris_push_rr_LL_u_3d(np, uxp, uyp, uzp, gaminv, exold, eyold, ezold, & 
           bxold, byold, bzold, ex, ey, ez, bx, by, bz, q, m, dt)
  USE constants, ONLY: clight, eps0, pi
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  ! Input/Output parameters
  INTEGER(idp), INTENT(IN) :: np
  REAL(num), INTENT(INOUT) :: uxp(np), uyp(np), uzp(np), gaminv(np)
  REAL(num), INTENT(IN)    :: exold(np), eyold(np), ezold(np)
  REAL(num), INTENT(IN)    :: bxold(np), byold(np), bzold(np)
  REAL(num), INTENT(IN)    :: ex(np), ey(np), ez(np)
  REAL(num), INTENT(IN)    :: bx(np), by(np), bz(np)
  REAL(num), INTENT(IN)    :: q, m, dt
  ! Local variables
  INTEGER(idp)             :: ip
  REAL(num)                :: const, qminv, crr, dtinv
  REAL(num)                :: clghtisq, usq, tsqi
  REAL(num)                :: tx, ty, tz
  REAL(num)                :: sx, sy, sz
  REAL(num)                :: uxppr, uyppr, uzppr
  REAL(num)                :: gaminvtmp, gamtmp
  REAL(num)                :: uxold, uyold, uzold
  REAL(num)                :: urx, ury, urz
  REAL(num)                :: ulx, uly, ulz
  REAL(num)                :: eup
  REAL(num)                :: minv, upinvsq
  REAL(num)                :: fsq
  REAL(num)                :: dex, dey, dez
  REAL(num)                :: dbx, dby, dbz
  REAL(num)                :: ebx, eby, ebz
  REAL(num)                :: bubx, buby, bubz
  REAL(num)                :: ubx, uby, ubz

  ! Initialization
  const = q*dt*0.5_num/m
  dtinv=1.0_num/dt
  clghtisq = 1.0_num/clight**2
  crr=q**3*dt/(6.0_num*m**2*clight**3*pi*eps0)
  ! inverse of mass
  qminv=q/m

  ! Loop over the particles
#if defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
  !DIR$ ASSUME_ALIGNED gaminv:64
  !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
  !DIR$ ASSUME_ALIGNED bx:64, by:64, bz:64
  !DIR$ ASSUME_ALIGNED exold:64, eyold:64, ezold:64
  !DIR$ ASSUME_ALIGNED bxold:64, byold:64, bzold:64

#elif defined __IBMBGQ__
  !IBM* ALIGN(64, uxp, uyp, uzp)
  !IBM* ALIGN(64, gaminv)
  !IBM* ALIGN(64, exold, eyold, ezold)
  !IBM* ALIGN(64, bxold, byold, bzold)
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

    ! --- Pushing with the Lorentz force

    ! Save previous moments
    uxold=uxp(ip)
    uyold=uyp(ip)
    uzold=uzp(ip)

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
    !write(*,*)"check 1",uxp(ip),uyp(ip),uzp(ip)
    ! Compute final Gamma
   usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
   gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)

    ! --- Pushing with the RR force
    
    ! Lorentz push
    ulx=(uxp(ip)-uxold)
    uly=(uyp(ip)-uyold)
    ulz=(uzp(ip)-uzold)
    
    ! Field derivatives
    dex = dtinv*(ex(ip)-exold(ip))
    dey = dtinv*(ey(ip)-eyold(ip))
    dez = dtinv*(ez(ip)-ezold(ip))
    dbx = dtinv*(bx(ip)-bxold(ip))
    dby = dtinv*(by(ip)-byold(ip))
    dbz = dtinv*(bz(ip)-bzold(ip))

    ! Momentum/mass at half step
    uxp(ip) = (uxp(ip) + uxold)*0.5_num
    uyp(ip) = (uyp(ip) + uyold)*0.5_num
    uzp(ip) = (uzp(ip) + uzold)*0.5_num

    ! Compute temporary Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gamtmp=sqrt(1.0_num + usq)
    gaminvtmp = 1.0_num/gamtmp

    ! 'up' cross B
    ubx=uyp(ip)*bz(ip)-uzp(ip)*by(ip)
    uby=uzp(ip)*bx(ip)-uxp(ip)*bz(ip)
    ubz=uxp(ip)*by(ip)-uyp(ip)*bx(ip)
    ! E cross B
    ebx=ey(ip)*bz(ip)-ez(ip)*by(ip)
    eby=ez(ip)*bx(ip)-ex(ip)*bz(ip)
    ebz=ex(ip)*by(ip)-ey(ip)*bx(ip)

    ! 'up' dot e
    eup=ex(ip)*uxp(ip)+ey(ip)*uyp(ip)+ez(ip)*uzp(ip)

    ! B cross up cross B
    bubx=by(ip)*ubz-bz(ip)*uby
    buby=bz(ip)*ubx-bx(ip)*ubz
    bubz=bx(ip)*uby-by(ip)*ubx

    ! Square of fields
    fsq=(ex(ip)+gaminvtmp*ubx)**2 +   		 & 
    (ey(ip)+gaminvtmp*uby)**2 +        		 &
    (ez(ip)+gaminvtmp*ubz)**2
    !write(*,*)"check 1",by(ip),ex(ip)
    
    ! RR force 
    urx= crr*(gamtmp*(dex+gaminvtmp*(uyp(ip)*dbz-uzp(ip)*dby))+		&
    qminv*((ebx-gaminvtmp*bubx+clghtisq*gaminvtmp*eup*ex(ip))-			&
    gamtmp*clghtisq*(fsq-gaminvtmp**2*clghtisq*(eup)**2)*uxp(ip)))       				
    
    ury= crr*(gamtmp*(dey+gaminvtmp*(uzp(ip)*dbx-uxp(ip)*dbz))+		&
    qminv*((eby-gaminvtmp*buby+clghtisq*gaminvtmp*eup*ey(ip))-			&
    gamtmp*clghtisq*(fsq-gaminvtmp**2*clghtisq*(eup)**2)*uyp(ip)))       							
    
    urz= crr*(gamtmp*(dez+gaminvtmp*(uxp(ip)*dby-uyp(ip)*dbx))+		&
    qminv*((ebz-gaminvtmp*bubz+clghtisq*gaminvtmp*eup*ez(ip))-			&
    gamtmp*clghtisq*(fsq-gaminvtmp**2*clghtisq*(eup)**2)*uzp(ip)))       
    
    ! Push using the RR force
    uxp(ip) = uxold+ulx+urx
    uyp(ip) = uyold+uly+ury
    uzp(ip) = uzold+ulz+urz

    ! Compute final Gamma
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
  ENDDO

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE pxr_boris_push_rr_LL_u_3d



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
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
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
  REAL(num)                :: uxold, uyold, uzold

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
END SUBROUTINE pxr_boris_push_u_3d_block

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
  USE picsar_precision, ONLY: idp, num
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
!$acc parallel deviceptr(xp, yp, zp, uxp, uyp, uzp, gaminv)
!$acc loop gang vector
  DO ip=1, np
    xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dt
    yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dt
    zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dt
  ENDDO
!$acc end loop
!$acc end parallel
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
  USE picsar_precision, ONLY: idp, num
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
  USE picsar_precision, ONLY: idp, num
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
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
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

!$acc parallel deviceptr(uxp, uyp, uzp, gaminv)
!$acc loop gang vector
  DO ip=1, np
    usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
    gaminv(ip) = 1.0_num/sqrt(1.0_num + usq)
  END DO
!$acc end loop
!$acc end parallel

#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif
  RETURN
END SUBROUTINE pxr_set_gamma

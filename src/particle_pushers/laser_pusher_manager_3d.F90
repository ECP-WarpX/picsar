! ________________________________________________________________________________________
!> @brief
!> Main subroutine for  laser and antenna
!> called in the main loop (see. submain.F90)
!
!> @details
!> This subroutine calls subroutines for the laser pusher
!
!> @author
!> H. Vincenti
!> @date
!> Creation 2017
! ________________________________________________________________________________________

SUBROUTINE push_laser_particles
  USE particles
  USE PICSAR_precision
  USE constants
  USE fields
  USE params
  USE shared_data
  USE tiling
  USE time_stat
  IMPLICIT NONE
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                       :: tdeb, tend, real_time

#if defined(DEBUG)
  WRITE(0, *) "push_laser_particles: start"
#endif
  real_time=it*dt   -dt/2.0_num

  IF (nspecies .EQ. 0_idp) RETURN
  tdeb=MPI_WTIME()
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) SHARED(ntilex,        &
  !$OMP ntiley, ntilez, nspecies, species_parray, npid, dt, real_time) PRIVATE(ix,    &
  !$OMP iy, iz, ispecies, curr, curr_tile, count)
  DO iz=1, ntilez! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        DO ispecies=1, nspecies! LOOP ON SPECIES
          ! - Get current tile properties
          ! - Init current tile variables
          curr=>species_parray(ispecies)
          IF (.NOT. curr%is_antenna) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) CYCLE
          IF(curr%antenna_params%time_window == 0) THEN
            CALL laserp_pusher_gaussian(count, npid, curr_tile%pid(1:count, 1:npid),  &
            curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,  &
            curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv, dt, 8_idp,   &
            curr%antenna_params%Emax, curr%antenna_params%Emax_laser_1,               &
            curr%antenna_params%Emax_laser_2, curr%antenna_params%polvector1,         &
            curr%antenna_params%polvector2, curr%antenna_params%k0_laser,             &
            curr%antenna_params%q_z, curr%antenna_params%laser_tau, real_time,        &
            curr%antenna_params%t_peak, curr%antenna_params%temporal_order,           &
            curr%antenna_params%polangle)
          ELSE
            CALL laserp_pusher_hanning(count, npid, curr_tile%pid(1:count, 1:npid),   &
            curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,  &
            curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv, dt, 8_idp,   &
            curr%antenna_params%Emax, curr%antenna_params%Emax_laser_1,               &
            curr%antenna_params%Emax_laser_2, curr%antenna_params%polvector1,         &
            curr%antenna_params%polvector2, curr%antenna_params%k0_laser,             &
            curr%antenna_params%q_z, real_time, curr%antenna_params%t_peak,           &
            curr%antenna_params%temporal_order, curr%antenna_params%polangle)
          ENDIF
        END DO! END LOOP ON SPECIES
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO
  tend=MPI_WTIME()
  pushtime=pushtime+(tend-tdeb)
  IF (it.ge.timestat_itstart) THEN
    tend=MPI_WTIME()
    localtimes(1) = localtimes(1) + (tend-tdeb)
  ENDIF

#if defined(DEBUG)
  WRITE(0, *) "push_laser_particles: stop"
#endif
END SUBROUTINE push_laser_particles

! ________________________________________________________________________________________
!> @brief
!> Subroutine for pushing particles of type antenna
!
!
!> @author
!> Haithem Kallala
!> @date
!> Creation 2017
! ________________________________________________________________________________________

SUBROUTINE laserp_pusher_gaussian(np, npidd, pid, xp, yp, zp, uxp, uyp, uzp, gaminv,  &
  dtt, lvect, emax, emax1, emax2, polvector1, polvector2, k0_laser, q_z, laser_tau,     &
  real_time, t_peak, temporal_order, polangle)
  USE shared_data
  USE omp_lib
  USE PICSAR_precision
  USE constants
  USE params
  USE particles
  USE particle_speciesmodule
  USE particle_properties
  USE antenna
  USE particle_tilemodule
  use fields
  INTEGER(idp), INTENT(IN)                :: np
  INTEGER(idp), INTENT(IN)                :: npidd
  INTEGER(idp), INTENT(IN)                :: lvect
  REAL(num), DIMENSION(1:np, 1:npidd), INTENT(IN)  :: pid
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp, yp, zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp, uyp, uzp, gaminv
  REAL(num), INTENT(IN)                   :: dtt
  REAL(num), DIMENSION(3), INTENT(IN)     :: polvector1, polvector2
  REAL(num), INTENT(IN)                   :: emax, emax1, emax2, k0_laser, laser_tau,  &
  real_time, t_peak, polangle
  COMPLEX(cpx), INTENT(IN)                :: q_z
  INTEGER(idp), INTENT(IN)                :: temporal_order
  INTEGER(idp)                            :: n, nn, ip, blocksize
  REAL(num)                               :: amp1, amp2, amp3
  REAL(num)                               :: xx, yy, clightsq, coeff_ampli,      &
  disp_max
  disp_max   = 0.01_num*clight
  coeff_ampli = disp_max / emax
  clightsq = 1._num/clight**2

  !____________________________________________________________________________
  ! Loop on block of particles of size lvect
  DO ip=1, np, lvect
    blocksize = MIN(lvect, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !!$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
    DO n=1, blocksize
      nn=ip+n-1
      xx = pid(nn, 2)
      yy = pid(nn, 3)
      CALL gaussian_profile(xx, yy, amp1, amp2, amp3, emax, emax1, emax2, polvector1, &
      polvector2, k0_laser, q_z, laser_tau, real_time, t_peak, temporal_order,        &
      polangle)
      ! --- Update particle momenta based on laser electric field
      uxp(nn) = amp1*coeff_ampli
      uyp(nn) = amp2*coeff_ampli
      uzp(nn) = amp3*coeff_ampli
      ! --- Update gaminv
      gaminv(nn) = 1.0_num
      ! --- Push x, y, z
      xp(nn)  = xp(nn) + dt*uxp(nn)
      IF(c_dim == 3) THEN
        yp(nn)  = yp(nn) + dt*uyp(nn)
      ENDIF
      zp(nn)  = zp(nn) + dt*uzp(nn)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !!$OMP END SIMD
#endif
#endif
  ENDDO
END SUBROUTINE laserp_pusher_gaussian
! ________________________________________________________________________________________
!> @brief
!> Subroutine for computing Hamming laser profile in time and space
!
!> @author
!> Haithem Kallala
!> @date
!> Creation 2017
! ________________________________________________________________________________________
SUBROUTINE laserp_pusher_hanning(np, npidd, pid, xp, yp, zp, uxp, uyp, uzp, gaminv,   &
  dtt, lvect, emax, emax1, emax2, polvector1, polvector2, k0_laser, q_z, real_time,     &
  t_peak, temporal_order, polangle)
  USE shared_data
  USE omp_lib
  USE PICSAR_precision
  USE constants
  USE params
  USE particles
  USE particle_speciesmodule
  USE particle_properties
  USE antenna
  USE particle_tilemodule
  use fields
  INTEGER(idp), INTENT(IN)                :: np
  INTEGER(idp), INTENT(IN)                :: npidd
  INTEGER(idp), INTENT(IN)                :: lvect
  REAL(num), DIMENSION(1:np, 1:npidd), INTENT(IN)  :: pid
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp, yp, zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp, uyp, uzp, gaminv
  REAL(num), INTENT(IN)                   :: dtt
  REAL(num), DIMENSION(3), INTENT(IN)    :: polvector1, polvector2
  REAL(num), INTENT(IN)                  :: emax, emax1, emax2, k0_laser, real_time,  &
  t_peak, polangle
  COMPLEX(cpx), INTENT(IN)                :: q_z
  INTEGER(idp), INTENT(IN)                :: temporal_order
  INTEGER(idp)                            :: n, nn, ip, blocksize
  REAL(num)                               :: amp1, amp2, amp3
  REAL(num)                               :: xx, yy, clightsq, coeff_ampli,      &
  disp_max
  disp_max   = 0.01_num*clight
  coeff_ampli = disp_max / emax
  clightsq = 1._num/clight**2

  !____________________________________________________________________________
  ! Loop on block of particles of size lvect
  DO ip=1, np, lvect
    blocksize = MIN(lvect, np-ip+1)
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED uxp:64, uyp:64, uzp:64
    !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
    !IBM* ALIGN(64, uxp, uyp, uzp)
    !IBM* ALIGN(64, gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !!$OMP SIMD
#endif
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
    DO n=1, blocksize
      nn=ip+n-1
      xx = pid(nn, 2)
      yy = pid(nn, 3)
      CALL hanning_profile(xx, yy, amp1, amp2, amp3, emax, emax1, emax2, polvector1,  &
      polvector2, k0_laser, q_z, real_time, t_peak, temporal_order, polangle)
      ! --- Update particle momenta based on laser electric field
      uxp(nn) = amp1*coeff_ampli
      uyp(nn) = amp2*coeff_ampli
      uzp(nn) = amp3*coeff_ampli
      ! --- Update gaminv
      gaminv(nn) = 1.0_num
      ! --- Push x, y, z
      xp(nn)  = xp(nn) + dt*uxp(nn)
      IF(c_dim ==3) THEN
        yp(nn)  = yp(nn) + dt*uyp(nn)
      ENDIF
      zp(nn)  = zp(nn) + dt*uzp(nn)
    ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
    !!$OMP END SIMD
#endif
#endif
  ENDDO
END SUBROUTINE laserp_pusher_hanning

! ________________________________________________________________________________________
!> @brief
!> Subroutine for computing gaussian laser profile in time and space
!
!> @author
!> Haithem Kallala
!> @date
!> Creation 2017
! ________________________________________________________________________________________

SUBROUTINE gaussian_profile(xx, yy, amp1, amp2, amp3, emax, emax1, emax2, polvector1, &
  polvector2, k0_laser, q_z, laser_tau, real_time, t_peak, temporal_order, polangle)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !!$OMP DECLARE SIMD(gaussian_profile) UNIFORM(emax, emax1, emax2, polvector1,        &
  !!$OMP polvector2, k0_laser, q_z, laser_tau, real_time, t_peak, temporal_order,      &
  !!$OMP polangle)
#endif
#elif defined __INTEL_COMPILER
  !DIR$ ATTRIBUTES VECTOR :                                  &
  !DIR$ UNIFORM(emax, emax1, emax2, emax, polvector1, polvector2, &
  !DIR$ k0_laser, q_z, laser_tau, real_time, t_peak, &
  !DIR$ temporal_order, polangle)  :: gaussian_profile
#endif
  USE PICSAR_precision
  USE constants
  USE params
  USE shared_data
  USE omp_lib
  REAL(num), INTENT(INOUT)   :: amp1, amp2, amp3
  REAL(num), DIMENSION(3), INTENT(IN)       :: polvector1, polvector2
  REAL(num), INTENT(IN)                      :: emax, emax1, emax2, k0_laser,         &
  laser_tau, real_time, t_peak, polangle
  COMPLEX(cpx), INTENT(IN)                  :: q_z
  REAL(num), INTENT(IN)                     :: xx, yy
  INTEGER(idp), INTENT(IN)                  :: temporal_order
  COMPLEX(cpx), DIMENSION(3)                :: arg
  COMPLEX(cpx)                               :: j, u1, u2

  j=(0.0_num, 1.0_num)
  IF (temporal_order .EQ. 0_idp) THEN 
    u1 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser*(xx**2+yy**2)/(2*q_z)
    u2 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser*(xx**2+yy**2)/(2*q_z)
  ELSE
    u1 = j*k0_laser*clight*(real_time-t_peak)- j*k0_laser*(xx**2+yy**2)/(2*q_z) -       &
    ((real_time - t_peak )/laser_tau)**temporal_order
    u2 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser*(xx**2+yy**2)/(2*q_z) -      &
    ((real_time - t_peak )/laser_tau)**temporal_order+polangle*2.0_num*pi*j
  ENDIF 
  u1 = EXP(u1)*emax1
  u2 = EXP(u2)*emax2
  arg(1) = (u1*polvector1(1) + u2*polvector2(1))
  amp1 = AIMAG(arg(1))
  arg(2) = (u1*polvector1(2) + u2*polvector2(2))
  amp2 = AIMAG(arg(2))
  arg(3) = (u1*polvector1(3) + u2*polvector2(3))
  amp3 = AIMAG(arg(3))

END SUBROUTINE gaussian_profile
! ________________________________________________________________________________________
!> @brief
!> Subroutine for computing Hamming laser profile in time and space
!
!> @author
!> Haithem Kallala
!> @date
!> Creation 2017
! ________________________________________________________________________________________

SUBROUTINE hanning_profile(xx, yy, amp1, amp2, amp3, emax, emax1, emax2, polvector1,  &
  polvector2, k0_laser, q_z, real_time, t_peak, temporal_order, polangle)
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !!$OMP DECLARE SIMD(hanning_profile) UNIFORM(emax, emax1, emax2, polvector1,         &
  !!$OMP polvector2, k0_laser, q_z, real_time, t_peak, temporal_order, polangle)
#endif
#elif defined __INTEL_COMPILER
  !DIR$ ATTRIBUTES VECTOR :                       &
  !DIR$ UNIFORM(emax, emax1, emax2, emax, &
  !DIR$ polvector1, polvector2, k0_laser, q_z, &
  !DIR$ real_time, t_peak, temporal_order, polangle) &
  !DIR$ :: hanning_profile
#endif

  USE shared_data
  USE PICSAR_precision
  USE constants
  USE params
  USE omp_lib
  REAL(num), INTENT(INOUT) :: amp1, amp2, amp3
  REAL(num), DIMENSION(3), INTENT(IN)    :: polvector1, polvector2
  REAL(num), INTENT(IN)                  :: emax, emax1, emax2, k0_laser, real_time,  &
  t_peak, polangle
  COMPLEX(cpx), INTENT(IN)                :: q_z
  REAL(num), INTENT(IN)                   :: xx, yy
  INTEGER(idp), INTENT(IN)                 :: temporal_order
  COMPLEX(cpx), DIMENSION(3)              :: arg
  COMPLEX(cpx)                            :: j, u1, u2
  REAL(num)                               :: Idd

  j=(0.0_num, 1.0_num)
  IF(real_time .LE. 2.0_num*t_peak) THEN
    Idd=1.0_num
  ELSE
    Idd=0.0_num
    amp1 = 0.0_num
    amp2 = 0.0_num
    amp3 = 0.0_num
    RETURN
  ENDIF
  u1 = -j*k0_laser*((xx**2+yy**2)/(2*q_z) - clight*(real_time-t_peak))
  u2 = -j*k0_laser*((xx**2+yy**2)/(2*q_z)                                             &
  -clight*(real_time-t_peak))+j*polangle*2.0_num*pi
  u1=Idd*EXP(u1)*emax1*(0.5_num -0.5_num*COS(real_time*2.0_num*pi/(2.0_num*t_peak)))
  u2=Idd*EXP(u2)*emax2*(0.5_num -0.5_num*COS(real_time*2.0_num*pi/(2.0_num*t_peak)))
  arg(1) = (u1*polvector1(1) + u2*polvector2(1))
  amp1 = REAL(arg(1), num)
  arg(2) = (u1*polvector1(2) + u2*polvector2(2))
  amp2 = REAL(arg(2), num)
  arg(3) = (u1*polvector1(3) + u2*polvector2(3))
  amp3 = REAL(arg(3), num)

END SUBROUTINE hanning_profile

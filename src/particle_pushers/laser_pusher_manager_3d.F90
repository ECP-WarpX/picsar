! ______________________________________________________________________________
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
SUBROUTINE push_laser_particles
  ! ______________________________________________________________________________
  USE particles
  USE constants
  USE fields
  USE params
  USE shared_data
  USE tiling
  IMPLICIT NONE
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                       :: tdeb, tend, real_time
  INTEGER(idp)                    :: nxc, nyc, nzc
  
#if defined(DEBUG)
  WRITE(0,*) "push_laser_particles: start"
#endif
  real_time=it*dt
  
  IF (nspecies .EQ. 0_idp) RETURN
  tdeb=MPI_WTIME()
  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,npid,dt,real_time) &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count)
  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        DO ispecies=1, nspecies ! LOOP ON SPECIES
          ! - Get current tile properties
          ! - Init current tile variables
          curr=>species_parray(ispecies)
          IF (.NOT. curr%is_antenna) CYCLE 
          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) CYCLE
          CALL laserp_pusher(count,npid,curr_tile%pid(1:count,1:npid),      &
          curr_tile%part_x,                                                 &
          curr_tile%part_y,curr_tile%part_z, curr_tile%part_ux,             &
          curr_tile%part_uy,                                                &
          curr_tile%part_uz,curr_tile%part_gaminv,dt,8_idp,                 &
          curr%antenna_params%Emax,                                         &
          curr%antenna_params%Emax_laser_1,                                 &
          curr%antenna_params%Emax_laser_2,                                 & 
          curr%antenna_params%polvector1,                                   & 
          curr%antenna_params%polvector2,                                   &
          curr%antenna_params%k0_laser,                                     &
          curr%antenna_params%q_z,                                          &
          curr%antenna_params%laser_tau,                                    &   
          real_time,                                                        &    
          curr%antenna_params%t_peak,                                       &
          curr%antenna_params%temporal_order,                               &
          curr%antenna_params%polangle)                  
        END DO! END LOOP ON SPECIES
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO
  tend=MPI_WTIME()
  pushtime=pushtime+(tend-tdeb)
#if defined(DEBUG)
  WRITE(0,*) "push_laser_particles: stop"
#endif
END SUBROUTINE push_laser_particles 

! ______________________________________________________________________________
!> @brief
!> Subroutine for pushing particles of type antenna
!
!> @details
!> This routine calls the subroutines for the different
!
!> @author
!> Haithem Kallala
!> @date
!> Creation 2017
SUBROUTINE laserp_pusher(np,npidd,pid,xp,yp,zp,uxp,uyp,uzp,gaminv,&
  dtt,lvect,emax,emax1,emax2,polvector1,polvector2,              &
  k0_laser,q_z,laser_tau,real_time,t_peak,temporal_order,polangle)
  USE shared_data
  USE omp_lib
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
  REAL(num), DIMENSION(1:np,1:npidd), INTENT(IN)  :: pid
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp,yp,zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp,uyp,uzp,gaminv
  REAL(num), INTENT(IN)                   :: dtt 
  REAL(num) , DIMENSION(3), INTENT(IN)    :: polvector1,polvector2
  REAL(num) , INTENT(IN)                  :: emax,emax1,emax2,k0_laser, laser_tau, &
  real_time,t_peak, polangle 
  COMPLEX(cpx), INTENT(IN)                :: q_z
  INTEGER(idp), INTENT(IN)                :: temporal_order 
  INTEGER(idp)                            :: n,nn,ip,i,j,k,blocksize
  REAL(num), DIMENSION(3)                 :: amp
  REAL(num)                               :: xx,yy,clightsq,usq, coeff_ampli, disp_max
  
  disp_max   = 0.01_num*clight
  coeff_ampli = disp_max / emax
  clightsq = 1._num/clight**2
  
  !____________________________________________________________________________
  ! Loop on block of particles of size lvect
  DO ip=1,np,lvect
    blocksize = MIN(lvect,np-ip+1)
    
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
#endif
    
#if defined _OPENMP && _OPENMP>=201307
#ifndef defined __IBMBGQ__
    !IBM* ALIGN(64,xp,yp,zp)
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
    !DIR$ SIMD
#endif
#endif
#if defined __INTEL_COMPILER
    !DIR$ IVDEP
    !!DIR DISTRIBUTE POINT
#endif
    DO n=1,blocksize
      
      nn=ip+n-1
      xx = pid(nn,2) 
      yy = pid(nn,3)
      CALL gaussian_profile(xx,yy,amp,emax,emax1,emax2,polvector1,polvector2,  &
      k0_laser,q_z,laser_tau,real_time,t_peak,temporal_order,polangle)
      ! --- Update particle momenta based on laser electric field 
      uxp(nn) = amp(1)*coeff_ampli
      uyp(nn) = amp(2)*coeff_ampli
      uzp(nn) = amp(3)*coeff_ampli
      ! --- Update gaminv 
      gaminv(nn) = 1.0_num
      ! --- Push x,y,z
      xp(nn)  = xp(nn) + dt*uxp(nn)! + dt*source_v(1) 
      yp(nn)  = yp(nn) + dt*uyp(nn)! + dt*source_v(2)
      zp(nn)  = zp(nn) + dt*uzp(nn)! + dt*source_v(3)
    ENDDO 
  ENDDO
 !print*,it,"max v ",maxval(abs(jy))/maxval(abs(uyp)),maxval(abs(ey))/maxval(abs(jy)),maxval(abs(ey))/emax1
END SUBROUTINE laserp_pusher 


SUBROUTINE gaussian_profile(xx,yy,amp,emax,emax1,emax2,polvector1,polvector2,  &
  k0_laser,q_z,laser_tau,real_time,t_peak,temporal_order,polangle)
  USE shared_data
  USE omp_lib
  USE constants
  USE params
  USE particles
  USE particle_speciesmodule
  REAL(num) , DIMENSION(3), INTENT(INOUT) :: amp
  REAL(num) , DIMENSION(3), INTENT(IN)    :: polvector1,polvector2
  REAL(num) , INTENT(IN)                  :: emax,emax1,emax2, k0_laser, laser_tau, &
  real_time,t_peak, polangle 
  COMPLEX(cpx), INTENT(IN)                :: q_z
  REAL(num) ,INTENT(IN)                   :: xx,yy
  INTEGER(idp), INTENT(IN)                :: temporal_order
  COMPLEX(cpx) , DIMENSION(3)             :: arg
  COMPLEX(cpx)                            :: j,u1,u2
  INTEGER(idp)                            :: i 
  j=(0.,1.)
  u1 = j*k0_laser*clight*(real_time-t_peak)- j*k0_laser*(xx**2+yy**2)/(2*q_z) &
  - ((real_time - t_peak )/laser_tau)**temporal_order
  
  u2 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser*(xx**2+yy**2)/(2*q_z) &
  - ((real_time - t_peak )/laser_tau)**temporal_order+polangle*2_num*pi
  u1 = EXP(u1)*emax1
  u2 = EXP(u2)*emax2
  DO i=1,3
    arg(i) = (u1*polvector1(i) + u2*polvector2(i))
  END DO
  amp = REAL(arg)
END SUBROUTINE 

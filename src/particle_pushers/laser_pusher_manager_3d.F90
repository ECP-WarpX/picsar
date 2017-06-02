! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Main subroutine for  laser and antenna
!> in the main loop (in submain.F90)
!
!> @details
!> This routine calls the subroutines for the different
!
!> @author
!> Haithem Kallala
!> @date
!> Creation 2017
SUBROUTINE laserp_pusher(np,npidd,pid,xp,yp,zp,uxp,uyp,uzp,gaminv,&
                        dtt,lvect)
  USE shared_data
  USE omp_lib
  USE constants
  USE params
  USE particles       
  USE particle_speciesmodule
  USE particle_properties
  USE antenna
  USE particle_tilemodule
  INTEGER(idp), INTENT(IN)                :: np
  INTEGER(idp), INTENT(IN)                :: npidd
  INTEGER(idp), INTENT(IN)                :: lvect
  REAL(num), DIMENSION(1:np,1:npidd), INTENT(IN)  :: pid
  REAL(num), DIMENSION(np), INTENT(INOUT) :: xp,yp,zp
  REAL(num), DIMENSION(np), INTENT(INOUT) :: uxp,uyp,uzp,gaminv
  REAL(num), INTENT(IN)                   :: dtt
  INTEGER(idp)                            :: n,nn,ip,i,j,k,blocksize
  REAL(num), DIMENSION(3)                 :: amp
  REAL(num)                               :: xx,yy,clightsq,usq

!    print*,(real_time-t_peak)/laser_tau,temporal_order,"temporal env",exp(- ((real_time - t_peak )/laser_tau)**temporal_order)
!    print*,"time phase ",cos(k0_laser*clight*(real_time-t_peak))
!    print*,"alltim cont",exp(- ((real_time - t_peak )/laser_tau)**temporal_order)*cos(k0_laser*clight*(real_time-t_peak))
!    print*,"waist and diff",laser_w0,diffract_factor
  clightsq = 1._num/clight**2

       print*,"exp = ",temporal_order

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
      CALL gaussian_profile(xx,yy,amp)
      uxp(nn) = amp(1) 
!print*,max(usq,abs(amp(2)),"a"
      uyp(nn) = amp(2)
      uzp(nn) = amp(3)
      usq = (uxp(nn)**2+uyp(nn)**2+uzp(nn)**2)*clightsq
      !gaminv(nn) = 1._num
      gaminv(nn) = 1.0_num/sqrt(1.0_num + usq)
      gaminv(nn) = 1.0_num
! print*,uyp(nn),(q_z/k0_laser),(real_time - t_peak ),laser_tau
!      xp(nn)  = xp(nn) + dt*uxp(nn)! + dt*source_v(1) 
!      yp(nn)  = yp(nn) + dt*uyp(nn)! + dt*source_v(2)
!      zp(nn)  = zp(nn) + dt*uzp(nn)! + dt*source_v(3)
       
    ENDDO 
  ENDDO
END SUBROUTINE laserp_pusher 
  

SUBROUTINE gaussian_profile(xx,yy,amp)

  USE shared_data
  USE omp_lib
  USE constants
  USE params
  USE particles
  USE particle_speciesmodule
  USE antenna

  REAL(num) , DIMENSION(3) ,INTENT(INOUT) :: amp
  REAL(num) ,INTENT(IN)                   :: xx,yy
  COMPLEX(cpx) , DIMENSION(3)             :: arg
  COMPLEX(cpx)                            :: j,u1,u2
  INTEGER(idp)                            :: i 
  

  j=(0.,1.)
  u1 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser*(xx**2+yy**2)/(2*q_z) &
        - ((real_time - t_peak )/laser_tau)**temporal_order

  u2 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser*(xx**2+yy**2)/(2*q_z) &
        - ((real_time - t_peak )/laser_tau)**temporal_order 


 ! u2 = j*k0_laser*clight*(real_time-t_peak) - j*k0_laser(xx**2+yy**2)/(2*q_z) &
!        - ((real_time - t_peak )/laser_tau)**temporal_order
 
!  print*,xx,yy,abs(exp(-j*k0_laser*(xx**2+yy**2)/(2*q_z))),1./sqrt(real(j*k0_laser/(2*q_z)))!,1./(log(abs(exp(-k0_laser*(xx**2+yy**2)/(2*q_z)))))*(x**2+y**2)
  u1 = EXP(u1)*laser_a_1*clight!Emax_laser_1
  u2 = EXP(u2)*laser_a_2*clight!/Emax_laser_2
  DO i=1,3
    arg(i) = (u1*polvector1(i) + u2*polvector2(i))/diffract_factor
  END DO 
  amp = REAL(arg)

END SUBROUTINE 
   
  













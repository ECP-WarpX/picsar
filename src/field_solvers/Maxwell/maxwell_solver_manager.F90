! ______________________________________________________________________________
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
! MAXWELL_SOLVER_MANAGER.F90
!
! Purpose:
! This file contains subroutines for the Maxwell solvers.
!
! Authors:
! Henri Vincenti
! Mathieu Lobet
!
! Date:
! Creation 2015
! ______________________________________________________________________________

! ______________________________________________________________________________
!> @brief
!> PUSH B field half a time step
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
SUBROUTINE push_bfield
! ______________________________________________________________________________
  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime

  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2).AND.(norderz.eq.2)) then
    CALL pxrpush_em3d_bvec(ex,ey,ez,bx,by,bz,                   &
              0.5_num*dt/dx,0.5_num*dt/dy,0.5_num*dt/dz,&
              nx,ny,nz,nxguards,nyguards,nzguards,nxs,nys,nzs, &
              l_nodalgrid)
  ! Yee scheme arbitrary order
  ELSE
    CALL pxrpush_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                       &
        0.5_num*dt/dx*xcoeffs,0.5_num*dt/dy*ycoeffs,0.5_num*dt/dz*zcoeffs,  &
        nx,ny,nz, norderx,nordery,norderz,                                  &
        nxguards,nyguards,nzguards,nxs,nys,nzs,                             &
        l_nodalgrid)
  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(5) = localtimes(5) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE push_bfield


! ______________________________________________________________________________
!> @brief
!> PUSH E field a full  time step
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
SUBROUTINE push_efield
! ______________________________________________________________________________
  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2).AND.(norderz.eq.2)) then
  CALL pxrpush_em3d_evec(ex,ey,ez,bx,by,bz,jx,jy,jz,clight**2*mu0*dt,        &
      clight**2*dt/dx,clight**2*dt/dy,                           &
      clight**2*dt/dz,nx,ny,nz,                                          &
      nxguards,nyguards,nzguards,nxs,nys,nzs,                                    &
      l_nodalgrid)

  ELSE
  ! Yee scheme arbitrary order
  CALL pxrpush_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,clight**2*mu0*dt,        &
      clight**2*dt/dx*xcoeffs,clight**2*dt/dy*ycoeffs,                           &
      clight**2*dt/dz*zcoeffs,nx,ny,nz,                                          &
      norderx,nordery,norderz,                                                   &
      nxguards,nyguards,nzguards,nxs,nys,nzs,                                    &
      l_nodalgrid)
  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(7) = localtimes(7) + (MPI_WTIME() - tmptime)
  ENDIF
END SUBROUTINE push_efield


! ______________________________________________________________________________
!> @brief
!> PUSH B field half a time step in 2D
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
SUBROUTINE push_bfield_2d
! ______________________________________________________________________________

  USE constants
  USE params
  USE fields
  USE shared_data
  USE time_stat
  IMPLICIT NONE

  REAL(num) :: tmptime
  IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
  ENDIF

  ! Yee scheme at order 2
  IF ((norderx.eq.2).AND.(nordery.eq.2)) then

    CALL pxrpush_em2d_bvec(ex,ey,ez,bx,by,bz,                  &
                           0.5_num*dt/dx,0._num,0.5_num*dt/dz,nx,0_idp,nz,&
                           nxguards,0_idp,nzguards,nxs,0_idp,nzs, &
                           l_nodalgrid)

  ! Yee scheme arbitrary order
  ELSE

    CALL pxrpush_em2d_bvec_norder(ex,ey,ez,bx,by,bz,                       &
      0.5_num*dt/dx*xcoeffs,0.5_num*dt/dy*ycoeffs,0.5_num*dt/dz*zcoeffs,  &
      nx,ny,nz, norderx,nordery,norderz,                                  &
      nxguards,nyguards,nzguards,nxs,nys,nzs,                             &
      l_nodalgrid)

  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(5) = localtimes(5) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE push_bfield_2d

! ______________________________________________________________________________
!> @brief
!> PUSH E,B PSAOTD a full time step 
!> This subroutine pushes the electric and the magnetic fields using
!> the PSATD solver.
!
!> @details
!> Time spent in this subroutine is stored in the electric field timer.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation March 29 2017
SUBROUTINE push_psatd_ebfield_3d
! ______________________________________________________________________________

  USE constants
  USE time_stat
  USE params
  USE shared_data
#if defined(FFTW)
  USE fourier_psaotd
#endif
  IMPLICIT NONE

  REAL(num) :: tmptime
  IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
  ENDIF
#if defined(FFTW)
   ! - Fourier Transform R2C
   IF (fftw_with_mpi) THEN 
		CALL get_Ffields_mpi ! - global FFT 
   ELSE
		CALL get_Ffields ! - local FFT  
   ENDIF 

   CALL push_psaotd_ebfielfs ! - PUSH PSATD 

   ! - Inverse Fourier Transform C2R
   IF (fftw_with_mpi) THEN 
		CALL get_fields_mpi  ! global IFFT
   ELSE
		CALL get_fields  ! local IFFT
   ENDIF 
#endif 
  IF (it.ge.timestat_itstart) THEN
    localtimes(7) = localtimes(7) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE


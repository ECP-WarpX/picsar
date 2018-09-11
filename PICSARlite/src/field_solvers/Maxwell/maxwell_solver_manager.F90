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
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> PUSH B field half a time step
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE push_bfield
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
    CALL pxrpush_em3d_bvec( &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         ex, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ey, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ez, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bx, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         by, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bz, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         0.5_num*dt/dx, 0.5_num*dt/dy, 0.5_num*dt/dz)
  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(5) = localtimes(5) + (MPI_WTIME() - tmptime)
  ENDIF

END SUBROUTINE push_bfield

! ________________________________________________________________________________________
!> @brief
!> PUSH E field a full  time step
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE push_efield
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

    CALL pxrpush_em3d_evec( &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         (/-nxs, -nys, -nzs/), (/nx+nxs, ny+nys, nz+nzs/), &
         ex, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ey, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         ez, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bx, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         by, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         bz, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         jx, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         jy, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         jz, (/-nxguards, -nyguards, -nzguards/), (/nx+nxguards, ny+nyguards, nz+nzguards/), &
         clight**2*mu0*dt, clight**2*dt/dx , clight**2*dt/dy, clight**2*dt/dz)
  ENDIF

  IF (it.ge.timestat_itstart) THEN
    localtimes(7) = localtimes(7) + (MPI_WTIME() - tmptime)
  ENDIF
END SUBROUTINE push_efield

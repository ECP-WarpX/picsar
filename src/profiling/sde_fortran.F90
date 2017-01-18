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
! SDE_FORTRAN.F90
!
! Tools to use the profiling tool Intel SDE with Picsar.
! This module also contains function for the Intel Design Forward Project.
!
! Author:
! Mathieu Lobet
!
! Date
! Creation 2016
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Module that contains interface with C functions used by the profiling
!> tool SDE.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
MODULE SDE_FORTRAN
! ______________________________________________________________________________
USE, INTRINSIC :: ISO_C_BINDING

INTERFACE

   SUBROUTINE FORTRAN_SDE_START() &
      BIND(C, NAME='fortran_sde_start')
   END SUBROUTINE FORTRAN_SDE_START

   SUBROUTINE FORTRAN_SDE_STOP() &
      BIND(C, NAME='fortran_sde_stop')
   END SUBROUTINE FORTRAN_SDE_STOP

   SUBROUTINE FORTRAN_DFP_INIT_START() &
      BIND(C, NAME='fortran_init_start')
   END SUBROUTINE FORTRAN_DFP_INIT_START

   SUBROUTINE FORTRAN_DFP_INIT_STOP() &
      BIND(C, NAME='fortran_init_stop')
   END SUBROUTINE FORTRAN_DFP_INIT_STOP

   SUBROUTINE FORTRAN_DFP_MAIN_START() &
      BIND(C, NAME='fortran_main_start')
   END SUBROUTINE FORTRAN_DFP_MAIN_START

   SUBROUTINE FORTRAN_DFP_MAIN_STOP() &
      BIND(C, NAME='fortran_main_stop')
   END SUBROUTINE FORTRAN_DFP_MAIN_STOP

   SUBROUTINE FORTRAN_DFP_FINAL_START() &
      BIND(C, NAME='fortran_finalize_start')
   END SUBROUTINE FORTRAN_DFP_FINAL_START

   SUBROUTINE FORTRAN_DFP_FINAL_STOP() &
      BIND(C, NAME='fortran_finalize_stop')
   END SUBROUTINE FORTRAN_DFP_FINAL_STOP

END INTERFACE

contains

  ! ____________________________________________________________________________
  !> @brief
  !> This function tells SDE to start the profiling collection.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  subroutine start_sde_collection()
  ! ____________________________________________________________________________
    write(0,*) "SDE profiling: start collecting data"
    call fortran_sde_start()
  end subroutine start_sde_collection

  ! ____________________________________________________________________________
  !> @brief
  !> This function tells SDE to stop the profiling collection.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  subroutine stop_sde_collection()
  ! ____________________________________________________________________________
    call fortran_sde_stop()
    write(0,*) "SDE profiling: stop collecting data"
  end subroutine stop_sde_collection

  ! ____________________________________________________________________________
  !> @brief
  !> This function used for the Design Forward Project should be called
  !> before the initialization.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
   Subroutine dfp_init_start()
  ! ____________________________________________________________________________
     CALL FORTRAN_DFP_INIT_START()
   end subroutine

  ! ____________________________________________________________________________
  !> @brief
  !> This function used for the Design Forward Project should be called
  !> after the initialization to stop the profiling of this specific part.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  Subroutine DFP_INIT_STOP()
  ! ____________________________________________________________________________
    CALL FORTRAN_DFP_INIT_STOP()
  end subroutine

  ! ____________________________________________________________________________
  !> @brief
  !> This function used for the Design Forward Project should be called
  !> before the main loop to start the profiling of this specific part.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  Subroutine DFP_MAIN_START()
  ! ____________________________________________________________________________
   CALL FORTRAN_DFP_MAIN_START()
  end subroutine

  ! ____________________________________________________________________________
  !> @brief
  !> This function used for the Design Forward Project should be called
  !> after the main loop to stop the profiling of this specific part.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  Subroutine DFP_MAIN_STOP()
  ! ____________________________________________________________________________
    CALL FORTRAN_DFP_MAIN_STOP()
  end subroutine

  ! ____________________________________________________________________________
  !> @brief
  !> This function used for the Design Forward Project should be called
  !> before the finalization to start the profiling of this specific part.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  Subroutine DFP_FINAL_START()
  ! ____________________________________________________________________________
    CALL FORTRAN_DFP_FINAL_START()
  end subroutine

  ! ____________________________________________________________________________
  !> @brief
  !> This function used for the Design Forward Project should be called
  !> after the finalization to stop the profiling of this specific part.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  Subroutine DFP_FINAL_STOP()
  ! ____________________________________________________________________________
    CALL FORTRAN_DFP_FINAL_STOP()
  end subroutine

END MODULE

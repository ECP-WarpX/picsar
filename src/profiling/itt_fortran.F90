! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016,
! The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a paid-up,
! nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
! copies to the public, prepare derivative works, and perform publicly and display
! publicly, and to permit other to do so.
!
! ITT_SDE_FORTRAN.F90
!
! This file contains interface to use Vtune functions with Fortran.
!
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> This module contains functions to use C Vtune routines with Fortran.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
MODULE ITT_FORTRAN
! ________________________________________________________________________________________
USE, INTRINSIC :: ISO_C_BINDING

INTERFACE

   SUBROUTINE FORTRAN_ITT_RESUME() &
      BIND(C, NAME='fortran_itt_resume')
   END SUBROUTINE FORTRAN_ITT_RESUME

   SUBROUTINE FORTRAN_ITT_PAUSE() &
      BIND(C, NAME='fortran_itt_pause')
   END SUBROUTINE FORTRAN_ITT_PAUSE
END INTERFACE

contains

  ! _____________________________________________________________________________________
  !> @brief
  !> Vtune starts collecting data when this subroutine is called.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  subroutine start_vtune_collection()
  ! _____________________________________________________________________________________
    write(0,*) "Vtune profiling: start collecting data"
    call fortran_itt_resume()
  end subroutine start_vtune_collection

  ! _____________________________________________________________________________________
  !> @brief
  !> Vtune starts collecting data when this subroutine is called.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  subroutine stop_vtune_collection()
  ! _____________________________________________________________________________________
    call fortran_itt_pause()
    write(0,*) "Vtune profiling: stop collecting data"
  end subroutine stop_vtune_collection

END MODULE

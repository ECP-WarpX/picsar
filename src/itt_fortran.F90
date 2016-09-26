! __________________________________________________________________
! ITT_SDE_FORTRAN.F90
! 
! Tools to use Vtune functions with fortran
! 
! __________________________________________________________________

MODULE ITT_FORTRAN
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

   subroutine start_vtune_collection()
     write(0,*) "Vtune profiling: start collecting data"
     call fortran_itt_resume()
   end subroutine start_vtune_collection

   subroutine stop_vtune_collection() 
    call fortran_itt_pause()
    write(0,*) "Vtune profiling: stop collecting data"
   end subroutine stop_vtune_collection

END MODULE

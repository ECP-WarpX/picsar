! __________________________________________________________________
! SDE_FORTRAN.F90
! 
! Tools to use SDE with Picsar
! 
! __________________________________________________________________

MODULE SDE_FORTRAN
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

   subroutine start_sde_collection()
     write(0,*) "SDE profiling: start collecting data"
     call fortran_sde_start()
   end subroutine start_sde_collection

   subroutine stop_sde_collection() 
    call fortran_sde_stop()
    write(0,*) "SDE profiling: stop collecting data"
   end subroutine stop_sde_collection

   Subroutine dfp_init_start()
     CALL FORTRAN_DFP_INIT_START()
   end subroutine

   Subroutine DFP_INIT_STOP()
     CALL FORTRAN_DFP_INIT_STOP()
   end subroutine

   Subroutine DFP_MAIN_START()
     CALL FORTRAN_DFP_MAIN_START()
   end subroutine

   Subroutine DFP_MAIN_STOP()
     CALL FORTRAN_DFP_MAIN_STOP()
   end subroutine

   Subroutine DFP_FINAL_START()
     CALL FORTRAN_DFP_FINAL_START()
   end subroutine

   Subroutine DFP_FINAL_STOP()
     CALL FORTRAN_DFP_FINAL_STOP()
   end subroutine

END MODULE

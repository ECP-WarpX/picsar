#include "ittnotify.h"

/* ____________________________ 

   VTUNE
  ____________________________ */

// Start Vtune analysis
void fortran_itt_resume()
{
    __itt_resume();
}

// Stop Vtune analysis
void fortran_itt_pause()
{
    __itt_pause();
}

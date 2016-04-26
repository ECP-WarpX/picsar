#include "ittnotify.h"

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

// Start sde analysis
void fortran_sde_start()
{
  __SSC_MARK(0x111);
}  

// Stop sde analysis
void fortran_sde_stop()
{
  __SSC_MARK(0x222);
} 

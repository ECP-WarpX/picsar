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

/* ____________________________ 

   SDE
  ____________________________ */

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

/* ____________________________ 

   OTHER INTEL TOOLS
  ____________________________ */
  
void fortran_init_start()
{
  __SSC_MARK(0x220);
}

void fortran_init_stop()
{
  __SSC_MARK(0x221);
}

void fortran_main_start()
{
  __SSC_MARK(0x44332211);
} 

void fortran_main_stop()
{
  __SSC_MARK(0x55332211);
} 

void fortran_finalize_start()
{
  __SSC_MARK(0x330);
}

void fortran_finalize_stop()
{
  __SSC_MARK(0x331);
}

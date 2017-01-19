/* _______________________________________________________________________________________

 *** Copyright Notice ***

 “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016,
 The Regents of the University of California, through Lawrence Berkeley National
 Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
 All rights reserved.

 If you have questions about your rights to use or distribute this software,
 please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

 NOTICE.
 This Software was developed under funding from the U.S. Department of Energy
 and the U.S. Government consequently retains certain rights. As such, the U.S.
 Government has been granted for itself and others acting on its behalf a paid-up,
 nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
 copies to the public, prepare derivative works, and perform publicly and display
 publicly, and to permit other to do so.

_______________________________________________________________________________________ */

#include "ittnotify.h"

/* _______________________________________________________________________________________

   VTUNE functions
_______________________________________________________________________________________ */

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

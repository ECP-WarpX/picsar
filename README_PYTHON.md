**INSTALLATION OF PYTHON MODULE PICSAR**
============================================================

**Compiling**
-------------
First edit the file Makefile_Forthon and indicate the following fields: 
FCOMP: fortran compiler 
FCOMPEXEC: openmpi Fortran wrapper 

Make -f Makefile_Forthon parse
Make -f Makefile_Forthon 


**Test***
-------------

To test: 
Make -f Makefile_forthon test
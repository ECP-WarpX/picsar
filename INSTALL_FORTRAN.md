**INSTALLING PICSAR IN FULL FORTRAN 90**
========================================


**Makefile config**
-------------------------

First edit the file Makefile and indicate the following environment variables:

* FC: your MPI Fortran compiler wrapper (e.g mpif90, mpiifort, ftn etc.),

* FARGS: your compiler arguments (optimization flags etc.). To get OpenMP version of PICSAR use the flag -fopenmp (with gfortran) and -openmp (Cray, Intel). 

**Compiling and Testing**
-------------------------

To compile and test, invoke the rule "all": 

* Make -f Makefile all


**INSTALLATION OF PYTHON MODULE PICSAR**
============================================================


**Makefile_Forthon config**
-------------------------

First edit the file Makefile_Forthon and indicate the following environment variables:

- FCOMP: your fortran compiler (e.g gfortran),

- FCOMPEXEC: your MPI Fortran wrapper (e.g mpif90),

- LIBDIR: your library folder containing MPI libraries (e.g /usr/local/Cellar/open-mpi/1.8.6/lib/ for an Homebrew install of open-mpi on MACOSX),

- LIBS: required libraries for the install. With Open-MPI, the compilation of picsar requires the following libraries: -lmpi, -lmpi_usempi, -lmpi_mpifh, -lgomp. Depending on 
your version of open-mpi, you should use -lmpi_usempif08 instead of -lmpi_usempi.  


**Compiling and Testing**
-------------------------

To compile and test, invoke the rule "all": 

- Make -f Makefile_Forthon all


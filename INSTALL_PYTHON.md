**INSTALLATION OF PYTHON MODULE PICSAR**
============================================================


**2. Installing python and packages **
--------------------------------------

* Install python 2 or later. We recommend python from anaconda (`http://docs.continuum.io/anaconda/install`)
* Install numpy (pip install numpy)
* Install mpi4py (pip install mpi4py)


**2. Installing Forthon **
-------------------------

Before creating the python module picsarpy for picsar, you must install the Forthon compiler. To do so: 

* Copy the last stable version of Forthon by typing: `git clone https://github.com/dpgrote/Forthon.git`

* Follow installation steps detailed in README 



**3. Makefile_Forthon config**
------------------------------

First edit the file Makefile_Forthon and indicate the following environment variables:

- FCOMP: your fortran compiler (e.g gfortran),

- FCOMPEXEC: your MPI Fortran wrapper (e.g mpif90),

- LIBDIR: your library folder containing MPI libraries (e.g /usr/local/Cellar/open-mpi/1.8.6/lib/ for an Homebrew install of open-mpi on MACOSX),

- LIBS: required libraries for the install. With Open-MPI, the compilation of picsar requires the following libraries: -lmpi, -lmpi_usempi, -lmpi_mpifh, -lgomp. Depending on 
your version of open-mpi, you should use -lmpi_usempif08 instead of -lmpi_usempi.  


**4. Compiling and Testing**
----------------------------

To compile and test, invoke the rule "all": 

- Make -f Makefile_Forthon all


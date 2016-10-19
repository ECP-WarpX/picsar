**INSTALLATION OF PYTHON MODULE PICSAR**
============================================================


**1. Installing python and packages **
--------------------------------------

* Install python 2 or later. We recommend python from anaconda (`http://docs.continuum.io/anaconda/install`)
* Install numpy (pip install numpy)
* Install mpi4py (pip install mpi4py)

NB: **On the cluster Edison at NERSC**:  
These packages cannot be installed, but have to be loaded instead.
To do so, please enter the following lines in your `.bashrc.ext`:
```
if [ "$NERSC_HOST" == "edison" ]
then
  module load python/2.7-anaconda
  module swap PrgEnv-intel PrgEnv-gnu
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
  export PYTHONPATH=$HOME/picsar/python_libs:$HOME/picsar/python_bin:$PYTHONPATH
fi
```
then, from your `$HOME`, type `source .bashrc`. The first line loads
`python`, `numpy`, `mpi4py`, etc. while the next lines prepare the
environment for the next installation steps.

**2. Installing Forthon **
-------------------------

Before creating the python module picsarpy for picsar, you must install the Forthon compiler. To do so: 

* Copy the last stable version of Forthon by typing: `git clone https://github.com/dpgrote/Forthon.git`

* Follow installation steps detailed in README

NB: **On the cluster Edison at NERSC**:  
Simply type `pip install Forthon --user`

**3. Makefile_Forthon config**
------------------------------

In order to automatically configure the installation, type
```
./configure
```

Then type
```
make -f Makefile_Forthon
```

**If the configure step failed:**   
You will need to edit the file `Makefile_Forthon` and indicate the following environment variables:

- FCOMP: your fortran compiler (e.g gfortran),

- FCOMPEXEC: your MPI Fortran wrapper (e.g mpif90),

- FARGS: arguments of the $(FCOMP) compiler. To get OpenMP version of PICSAR use the flag -fopenmp (with gfortran) and -openmp (Cray, Intel). NB: this version of PICSAR requires at least **OpenMP 4.0**.  

- LIBDIR: your library folder containing MPI libraries (e.g /usr/local/Cellar/open-mpi/1.8.6/lib/ for an Homebrew install of open-mpi on MACOSX, /opt/local/lib/mpich-mp/ for a Macports install of mpich),

- LIBS: required libraries for the install. With Open-MPI, the compilation of picsar requires the following libraries: -lmpi, -lmpi_usempi, -lmpi_mpifh, -lgomp. For open-mpi>1.8.x, you should use -lmpi_usempif08 instead of -lmpi_usempi. For a Macports install of mpich, you should use -lmpifort -lmpi -lpmpi.   

**On the cluster Edison at NERSC**:  
Clone `picsar` in your `$HOME` folder. Then modify the `Makefile_Forthon` to use the following configuration
```
SRCDIR= src
BINDIR = python_bin
APPNAME=picsar
PYTHON_NAME=picsarpy
UTIL=utils
FC=Forthon
FCARGS=-v --no2underscores  --nowritemodules
FCOMP=gfortran
FCOMPEXEC=ftn
FARGS="-O3 -fopenmp -ffree-line-length-none -ftree-vectorize -ftree-vectorizer-verbose=0"
LIBDIR=
LIBS= -lgomp
TESTDIR=example_scripts_python
```

**4. Compiling and installing**
----------------------------

To compile and test, invoke the rule "all": 
```
make -f Makefile_Forthon all
```
Then make sure that the folders `python_libs` and `python_bin` are in
your `$PYTHONPATH`. (On Edison, this is ensured by the code added in
your `.bashrc.ext`.)

**5. Testing**
----------------------------

Testing the code after compilation is highly recommended.

To test the compilation/execution, you can use the makefile (py.test is required):

```
make -f Makefile_Forthon test2
```

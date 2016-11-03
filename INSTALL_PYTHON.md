# **INSTALLATION OF PYTHON MODULE PICSAR**

## **1. Installing python and packages**

### **On MacOS**

If you have already installed Warp, this step is already done.
If not, we simply recommend to first install Warp.

Install python 2.7.
We recommend python from Anaconda (`http://docs.continuum.io/anaconda/install`)

Then, install numpy

``` 
pip install numpy
```

Finally, install mpi4py to be able to use MPI with Python.

```
pip install mpi4py
```

### **On Edison**
 
These packages are available via the modules.
To do so, please enter the following lines in your `.bashrc.ext`.

Using the GNU compiler:

```
if [ "$NERSC_HOST" == "edison" ]
then
  module load python/2.7-anaconda
  module swap PrgEnv-intel PrgEnv-gnu
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
  export PYTHONPATH=$HOME/picsar/python_libs:$HOME/picsar/python_bin:$PYTHONPATH
fi
```

Using the Intel compiler:

```
if [ "$NERSC_HOST" == "edison" ]
then
  module load python/2.7-anaconda
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
  export PYTHONPATH=$HOME/picsar/python_libs:$HOME/picsar/python_bin:$PYTHONPATH
fi
```

then, from your `$HOME`, type `source .bashrc`. The first line loads
`python`, `numpy`, `mpi4py`, etc. while the next lines prepare the
environment for the next installation steps.

### **On Cori**

These packages are available via the modules.
To do so, please enter the following lines in your `.bashrc.ext`.

Using the GNU compiler:

```
if [ "$NERSC_HOST" == "cori" ]
then
  module load python/2.7-anaconda
  module swap PrgEnv-intel PrgEnv-gnu
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
  export PYTHONPATH=$HOME/picsar/python_libs:$HOME/picsar/python_bin:$PYTHONPATH
fi
```

Using the Intel compiler:

```
if [ "$NERSC_HOST" == "cori" ]
then
  module load python/2.7-anaconda
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
  export PYTHONPATH=$HOME/picsar/python_libs:$HOME/picsar/python_bin:$PYTHONPATH
fi
```

then, from your `$HOME`, type `source .bashrc`. The first line loads
`python`, `numpy`, `mpi4py`, etc. while the next lines prepare the
environment for the next installation steps.

## **2. Installing Forthon**

If you have already installed Warp, this step is already done.

Before creating the python module picsarpy for picsar, 
you must install the Forthon compiler. 
To do so: 

* Copy the last stable version of Forthon by typing:
```
git clone https://github.com/dpgrote/Forthon.git
```

* Then `cd` into the directory `Forthon` and run python `setup.py install --home=$PATH`

Here, `PATH` is where you want the `bin` and `lib` folder to be created 
and the Forthon files to be located.

NB: **On the cluster Edison at NERSC**:  
Simply type `pip install Forthon --user`

## **3. Makefile_Forthon configuration**

### **On MacOs**

Clone `picsar` where you want to install it.

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

### **On the cluster Edison at NERSC**

Clone `picsar` where you want to install it. 
If you have already installed Warp, 
we recommend to put it in the same installation directory: `$SCRATCH/warp_install`. 

Then modify the `Makefile_Forthon` to have the correct configuration.

To use the GNU compiler:
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

To use the Intel compiler
```
SRCDIR= src
BINDIR = python_bin
APPNAME=picsar
PYTHON_NAME=picsarpy
UTIL=utils
FC=Forthon
FCARGS=-v --no2underscores  --nowritemodules
FCOMP=intel
FCOMPEXEC=ftn
FARGS="-O3 -qopenmp -xAVX"
LIBDIR=
LIBS=
TESTDIR=example_scripts_python
```

### **On the cluster Cori at NERSC**

Clone `picsar` where you want to install it. 
If you have already installed Warp, 
we recommend to put it in the same installation directory: `$SCRATCH/warp_install`. 

Then modify the `Makefile_Forthon` to have the correct configuration.

To use the GNU compiler:
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

To use the Intel compiler on Haswell architecture (-xCORE-AVX2):
```
SRCDIR= src
BINDIR = python_bin
APPNAME=picsar
PYTHON_NAME=picsarpy
UTIL=utils
FC=Forthon
FCARGS=-v --no2underscores  --nowritemodules
FCOMP=intel
FCOMPEXEC=ftn
FARGS="-O3 -qopenmp -xCORE-AVX2"
LIBDIR=
LIBS=
TESTDIR=example_scripts_python
```

To use the Intel compiler on KNL architecture (-xMIC-AVX512):
```
SRCDIR= src
BINDIR = python_bin
APPNAME=picsar
PYTHON_NAME=picsarpy
UTIL=utils
FC=Forthon
FCARGS=-v --no2underscores  --nowritemodules
FCOMP=intel
FCOMPEXEC=ftn
FARGS="-O3 -qopenmp  -xMIC-AVX512"
LIBDIR=
LIBS=
TESTDIR=example_scripts_python
```

## **4. Compiling and installing**


To compile and test, invoke the rule "all": 
```
make -f Makefile_Forthon all
```
Then make sure that the folders `python_libs` and `python_bin` are in
your `$PYTHONPATH`.

On Edison and Cori, this is ensured by the setup of your `~/.bashrc.ext`.

## **5. Testing**

### **On MacOs**

Testing the code after compilation is highly recommended.

To test the compilation/execution, you can use the makefile (py.test is required):

  For all test:
  - `make -f Makefile_Forthon test`

  For each test one by one
  - Langmuir wave:           `make -f Makefile_Forthon test2`


# **Installation of the PICSAR Python module on Cori**

Cori is a system equipped of two kinds of processors with different architectures.
Cori users have access to Intel Haswell and Intel KNL partitions.
In each case, the code has to be compiled differently.
A code compiled on Haswell will work on KNL but not efficiently.
A code compiled for KNL does not work on Haswell.

## **1. Installing WARP**

In order to install WARP, you can use the WARP installation documentation.
Download WARP from the [bitbucket repository](https://bitbucket.org/berkeleylab/warp).
Instructions are located in the sources in the folder `doc`.

## **2. Installing python and packages**

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

## **3. Installing Forthon**

If you have already installed Warp, this step is already done.
If not, simply type `pip install Forthon --user`

## **4. Makefile_Forthon configuration**

### **On the cluster Cori for Intel Haswell at NERSC**

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
FARGS="-O3 -qopenmp -xCORE-AVX2 -align array64byte"
LIBDIR=
LIBS=
TESTDIR=example_scripts_python
```

### **On the cluster Cori for Intel KNL at NERSC**

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
FARGS="-O3 -fopenmp -march=knl -ffree-line-length-none -ftree-vectorize -ftree-vectorizer-verbose=0"
LIBDIR=
LIBS= -lgomp
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


## **5. Compiling and installing**

To compile and test, invoke the rule "all":
```
make -f Makefile_Forthon all
```
Then `cd` into the directory `python_module` and run `python setup.py install --user`.


## **6. Python import**

To test the import of the python module, try in a `python` shell:
```
from picsar_python import picsarpy as pxrpy
```

(Note that it must be used in parallel on NERSC machines.)

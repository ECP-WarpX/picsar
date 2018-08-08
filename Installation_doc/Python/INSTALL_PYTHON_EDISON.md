# **Installation of the PICSAR Python module on Edison**

## **1. Installing WARP**

In order to install WARP, you can use the WARP installation documentation.
Download WARP from the [bitbucket repository](https://bitbucket.org/berkeleylab/warp).
Instructions are located in the sources in the folder `doc`.

## **2. Installing python and packages**

These packages are available via the modules.
To do so, please enter the following lines in your `.bashrc.ext`.

Using the GNU compiler:

```
if [ "$NERSC_HOST" == "edison" ]
then
  module load python/2.7-anaconda
  module swap PrgEnv-intel PrgEnv-gnu
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
fi
```

Using the Intel compiler:

```
if [ "$NERSC_HOST" == "edison" ]
then
  module load python/2.7-anaconda
  export PATH=$HOME/.local/edison/2.7-anaconda/bin:$PATH
fi
```

then, from your `$HOME`, type `source .bashrc`. The first line loads
`python`, `numpy`, `mpi4py`, etc. while the next lines prepare the
environment for the next installation steps.


## **3. Installing Forthon**

If you have already installed Warp, this step is already done.
If not, simply type `pip install Forthon --user`.

## **4. Makefile_Forthon configuration**

Clone `picsar` where you want to install it.
If you have already installed Warp,
we recommend to put it in the same installation directory: `$SCRATCH/warp_install`.

After `cd` in the `picsar` directory, use the command `./configure` or 'python configure --pxr_spectral_hybrid True'  to prepare
the `Makefile`.  After this step, a file `Makefile_Forthon` must have been generated
and must look like:

```
# Source directory
SRCDIR= src
# Binary directory (.so) after compilation
BINDIR = python_module/picsar_python
# Application name
APPNAME=picsar
# Python binary name
PYTHON_NAME=picsarpy
# Path where the parser is located
UTIL=utils/forthon_parser
# We use Forthon to interface Fortran and Python
FC=Forthon
FCARGS=-v --no2underscores  --nowritemodules
# Fortran compiler
FCOMP=gfortran
FCOMPEXEC=ftn
# Fortan compilation arguments
FARGS="-O3 -fopenmp -ffree-line-length-none -ftree-vectorize -ftree-vectorizer-verbose=0"
# Library directory
LIBDIR=
# Library names
LIBS=-lgomp
# Location of the test scripts
TESTDIR=examples/example_scripts_python
```

## **5. Compiling and installing**

To compile and test, invoke the rule "all":
```
make -f Makefile_Forthon all
```
Then `cd` into the directory `python_module` and run `python setup.py install --user`.


## **6. Python import**

To test the import of the python module, try in parallel within an interactive
shell:

```
salloc -N 1 --qos premium -t 00:30:00 -L SCRATCH
```

Then:

```
srun -n 2 python-mpi -c 'from picsar_python import picsarpy as pxrpy'
```

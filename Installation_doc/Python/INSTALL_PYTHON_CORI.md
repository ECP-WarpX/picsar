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
  export PATH=$HOME/.local/cori/2.7-anaconda/bin:$PATH
fi
```

Using the Intel compiler:

```
if [ "$NERSC_HOST" == "cori" ]
then
  module load python/2.7-anaconda
  export PATH=$HOME/.local/cori/2.7-anaconda/bin:$PATH
fi
```

then, from your `$HOME`, type `source .bashrc`. The first line loads
`python`, `numpy`, `mpi4py`, etc. while the next lines prepare the
environment for the next installation steps.

## **3. Installing Forthon**

If you have already installed Warp, this step is already done.
If not, simply type `pip install Forthon --user`.

## **4. Makefile_Forthon configuration**

### **On the cluster Cori for Intel Haswell at NERSC**

Clone `picsar` where you want to install it.
If you have already installed Warp,
we recommend to put it in the same installation directory: `$SCRATCH/warp_install`.

After `cd` in the `picsar` directory, use the command `./configure` or 'python configure --pxr_spectral_hybrid True' (to use Hybrid PSATD) to prepare
the `Makefile`.
In addition, you have the possibility to use either the GNU or
the Intel compiler for an optimized compilation for either Haswell or
knl architecture. The corresponding flags are `--compiler gnu`/`--compiler intel`
and `--architecture cpu`/`--architecture knl`.

For example, if you want to use `picsar` on the knl architecture with the intel
compiler, juste type `./configure --compiler intel --architecture knl`. Make sure you have the `craype-haswell` module loaded to compile on haswell, and `craype-mic-knl` module loaded to compile on KNL.
This step should return 'Configure succeeded'. Then, a file  `Makefile_Forthon` must have been generated and must look like:

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
FCOMP=intel
FCOMPEXEC=ftn
# Fortan compilation arguments
FARGS="-O3 -qopenmp  -xMIC-AVX512"
# Library directory
LIBDIR=
# Library names
LIBS=
# Location of the test scripts
TESTDIR=examples/example_scripts_python
```

Note that by default, the flags are respectively set to `gnu` and `cpu`.

For good performances on knl architecture, we recommend to use the Intel compiler.

## **5. Compiling and installing**

To compile and test, invoke the rule "all":
```
make -f Makefile_Forthon all
```
Then `cd` into the directory `python_module` and run `python setup.py install --user` to install the PICSAR python library into your home. If you want to install it in a custom directory, use `python setup.py install --home path/to/my/directory` instead.


## **6. Python import**

To test the import of the python module, try in parallel within an interactive
shell:

`salloc -N 1 --qos interactive -t 00:30:00 -C knl -L SCRATCH` for KNL

`salloc -N 1 --qos interactive -t 00:30:00 -C haswell -L SCRATCH` for Haswell

Then:

```
srun -n 2 python-mpi -c 'from picsar_python import picsarpy as pxrpy'
```

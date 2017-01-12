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


## **3. Installing Forthon**

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

## **4. Makefile_Forthon configuration**

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
FARGS="-O3 -qopenmp -xAVX -align array64byte"
LIBDIR=
LIBS=
TESTDIR=example_scripts_python
```

## **5. Compiling and installing**

To compile and test, invoke the rule "all": 
```
make -f Makefile_Forthon all
```
Then make sure that the folders `python_libs` and `python_bin` are in
your `$PYTHONPATH`.

On Edison and Cori, this is ensured by the setup of your `~/.bashrc.ext`.


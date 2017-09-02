# **Installation of the PICSAR Python module on MIRA**

## **1. Installing WARP**

In order to install WARP, you can use the WARP installation documentation.
Download WARP from the [bitbucket repository](https://bitbucket.org/berkeleylab/warp).
Instructions are located in the sources in the folder `doc/INSTALL_MIRA/INSTALL_MIRA.md`.

## **2. Installing python and packages**

For environment configuration and Forthon installation, refer to the `doc/INSTALL_MIRA/INSTALL_MIRA.md` 
documentation file on the WARP bitbucket repo.

## **3. Installing Forthon**

For environment configuration and Forthon installation, refer to the `doc/INSTALL_MIRA/INSTALL_MIRA.md`
documentation file on the WARP bitbucket repo.

## **4. Makefile_Forthon configuration**

To use gcc compiler on powerpc arch and compile the python version of PXR use the following options in your own  `Makefile_Forthon` duplicated from `Makefile_Forthon.in`:
```
SRCDIR= src
BINDIR = python_bin
APPNAME=picsar
PYTHON_NAME=picsarpy
UTIL=utils
FC=Forthon
FCARGS=-v --no2underscores  --nowritemodules
CCOMPILER=/bgsys/drivers/V1R2M2/ppc64/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-gcc -mcmodel=medium
LDCOM=mpicc -shared
FCOMP=gfortran
FCOMPEXEC=mpif90
FARGS="-O3 -fPIC -ffree-line-length-none -ftree-vectorize -ftree-vectorizer-verbose=0 -fopenmp"
LIBDIR=
LIBS=-lgomp -lmpifort#-lgomp -lmpich-gcc -lopa-gcc -lmpl-gcc -lpami-gcc #-lmpifort -lgfortran #-lmpi -lmpi_usempif08 -lmpi_mpifh -lgomp
TESTDIR=example_scripts_python
export CC=$(CCOMPILER)
export LDSHARED=$(LDCOM)
```

## **5. Compiling and installing**

To compile and test, invoke the rule "all":
```
make -f Makefile_Forthon all
```
Then make sure that the folders `python_libs` and `python_bin` are in
your `$PYTHONPATH`.

On Edison and Cori, this is ensured by the setup of your `~/.bashrc.ext`.

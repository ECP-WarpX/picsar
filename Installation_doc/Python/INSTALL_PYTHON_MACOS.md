# **Installation of the PICSAR Python module on MacOS**

The Python version of PICSAR has been developed for the PIC code WARP.

## **1. Installing WARP**

In order to install WARP, you can use the WARP installation documentation.
Download WARP from the [bitbucket repository](https://bitbucket.org/berkeleylab/warp).
Instructions are located in the sources in the folder `doc`.

## **2. Installing python and packages for PICSAR**

If you have already installed Warp,
this step (installing python and useful package) is already done.
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

## **3. Installing Forthon**

If you have already installed Warp, this step is already done.
If not, simply type `pip install Forthon`

## **4. Makefile_Forthon configuration**

Clone `picsar` where you want to install it.

In order to automatically configure the installation, type
```
./configure
```
In order to use highly scalable hybrid PSATD with warp you need to configure python with --pxr_spectral_hybrid flag
python configure --pxr_spectral_hybrid True


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


## **5. Compiling and installing**


To compile and test, invoke the rule "all":
```
make -f Makefile_Forthon all
```
Then `cd` into the directory `python_module` and run `python setup.py install`.


## **6. Testing**

Testing the code after compilation is highly recommended.

To test the compilation/execution, you can use the makefile (needs the python package pytest):

  For all test:
  - `make -f Makefile_Forthon test`

  For each test one by one
  - Propagating beams:       `make -f Makefile_Forthon test1`
  - Langmuir wave:           `make -f Makefile_Forthon test2`

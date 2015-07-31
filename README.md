**Particle-In-Cell Scalable Application Resource (PICSAR)**
============================================================

**Overview**
------------

The PICSAR code is a "mini-app" standalone Particle-In-Cell (PIC) code that includes
the key functionalities of the WARP code main PIC loop. It is a 
compact **self-contained proxy** that adequately portrays the computational loads
and dataflow of the more complex WARP code. 

Since WARP is a very large code written in a mix of FORTRAN95, C and Python 
PICSAR will be essential for studying multi-level parallelization on the next
generation of exascale computers. 

Here are some of the specific algorithmic features of the PICSAR code :  

* The Maxwell solver uses arbitrary order finite-difference scheme (staggered/centered), 
* The particle pusher uses the Boris algorithm,
* The field gathering routine is energy conserving, 
* The current deposition and field gathering routines include high order particle shape factors.

Here are some high performance features of the PICSAR code :

* MPI parallelization for internode parallelism, 
* OpenMP parallelization for intranode parallelism,
* MPI-IO for fast parallel outputs.

**Compiling**
-------------

To build the code just type "make" in the command prompt. Make file options can be changed by editing the `Makefile`. For the gfortran compiler, simply use the flag -fopenmp to add openMP features. To set "x" OpenMP threads per MPI task, use "export OMP_NUM_THREADS=x" before starting the simulation (default will be x=1)

**Running simulations**
-----------------------

PICSAR input parameters must be provided in an input file named "input_file.pxr" in the folder where the code is ran. An example of input file is provided in this repository.

To run the executable on n MPI processes: "mpirun -np n a.out". Notice that n must be equal to nprocx*nprocy*nprocz with nprocx, nprocy, nprocz the number of processors along x,y,z directions. 

For the moment, the code outputs binary matrix files with extensions ".pxr" that can be read using python scripts. Examples of such scripts are in the folder `postproc/`
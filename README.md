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

* Particle tiling to help increase memory locality
* MPI parallelization for internode parallelism (blocking, non-blocking and Remote memory access MPI), 
* OpenMP parallelization for intranode parallelism,
* MPI-IO for fast parallel outputs.

**Compiling**
-------------

To build the code just type "make" in the command prompt. Make file options can be changed by editing the `Makefile`. For the gfortran compiler, simply use the flag -fopenmp to add openMP features. For the intel compiler, use -openmp instead. 

**Running simulations**
-----------------------

PICSAR input parameters must be provided in an input file named "input_file.pxr" in the folder where the code is ran. An example of input file is provided in this repository.

To run the executable on n MPI processes: "mpirun -np n ./picsar". Notice that if nprocx, nprocy and nprocz are provided in the input file as part of the "cpusplit" section, then n must be equal to nprocx x nprocy x nprocz with nprocx, nprocy, nprocz the number of processors along x,y,z directions. Otherwise, if nprocx, nprocy and nprocz are not defined, the code performs automatic CPU split in each direction. 

To set "x" OpenMP threads per MPI task, use "export OMP_NUM_THREADS=x" before starting the simulation (default will be x=1). OpenMP scheduling for balancing loads between tiles can be adjusted at runtime by setting the environment variable OMP_SCHEDULE to either static, guided or dynamic. To ensure that threads have enough memory space on the stack, set OMP_STACKSIZE to high enough value. In practice, export OMP_STACKSIZE=32M should be sufficient for most of test cases. 

User can specify some arguments in the command line. For the moments this feature supports only the number of tiles in each dimension and the init of particle distribution. Ex: 

mpirun -np 1 ./picsar -ntilex ntx -ntiley nty -ntilez ntz -distr 1

with ntx, nty and ntz the number of tiles in each dimension (default is one) and distr the type of particle init ("1" for init on the x-axis of the grid and "2" for Random). 

**Outputs**
-----------------------
For the moment, the code outputs binary matrix files with extensions ".pxr" that can be read using python scripts. Examples of such scripts are in the folder `postproc/`. Output frequency is controlled by setting the flag output_frequency in the output section of the input_file.pixr. Use output_frequency=-1 to disable outputs. The code places output files in a "RESULTS" directory where the code is ran. This directory has to be created before running the code in your submission script. 
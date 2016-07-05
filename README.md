**Particle-In-Cell Scalable Application Resource (PICSAR)**
============================================================

**1. Overview**
------------

The PICSAR code is a "mini-app" standalone Particle-In-Cell (PIC) code that includes
the key functionalities of the WARP code main PIC loop. It is a 
compact **self-contained proxy** that adequately portrays the computational loads
and dataflow of the more complex WARP code. 

Since WARP is a very large code written in a mix of FORTRAN95, C and Python 
PICSAR will be essential for studying multi-level parallelization on the next
generation of exascale computers. 

####A.  Here are some of the specific algorithmic features of the PICSAR code :  

* The Maxwell solver uses arbitrary order finite-difference scheme (staggered/centered), 
* The particle pusher uses the Boris algorithm,
* The field gathering routine is energy conserving, 
* The current deposition and field gathering routines include high order particle shape factors.

####B.  Here are some high performance features of the PICSAR code :

* Particle tiling to help increase memory locality
* MPI parallelization for internode parallelism (blocking, non-blocking and Remote memory access MPI), 
* OpenMP parallelization for intranode parallelism,
* MPI-IO for fast parallel outputs.
* Vectorized subroutines (field gathering, classical current deposition)

####C.  Python glue: 

* We created a Forthon parser that read Fortran source files of PICSAR and parse them to create a `picsar.v` file used by the Forthon compiler to generate a Python module for PICSAR. The Forthon parser is available in the folder `utils`. 
* Thanks to Forthon, we are able to access all high performance routines of PICSAR from python. This allows us to use PICSAR routines from WARP and vice-versa. 


**2. Compiling**
-------------

####A.  Python installation 

In order to install picsar in the form of a Python module, read detailed instructions in the file `INSTALL_PYTHON.md`

####B.  Fortran installation 

To build the code in full Fotran 90 read instructions from the file  `INSTALL_FORTRAN.md` 

**3. Running simulations**
-----------------------

####A.  Python mode

An example of python script `test.py` is provided in `example_scripts_python`. To run this script in parallel, simply type :
```
mpirun -np NMPI python test.py 
```
with NMPI the number of MPI processes. 

####B.  Fortran mode

PICSAR input parameters must be provided in an input file named `input_file.pxr` in the folder where the code is ran. An example (`test.pxr`) of input file is provided in `example_decks_fortran/`. To run the executable on n MPI processes:
```
mpirun -np n ./picsar
```
Notice that if `nprocx`, `nprocy` and `nprocz` are provided in the input file as part of the `cpusplit` section, then n must be equal to `nprocx x nprocy x nprocz` with `nprocx`, `nprocy`, `nprocz` the number of processors along x,y,z directions. Otherwise, if `nprocx`, `nprocy` and `nprocz` are not defined, the code performs automatic CPU split in each direction. User can specify some arguments in the command line. For the moments this feature supports only the number of tiles in each dimension and the init of particle distribution. Ex: `mpirun -np 1 ./picsar -ntilex ntx -ntiley nty -ntilez ntz -distr 1` with `ntx`, `nty` and `ntz` the number of tiles in each dimension (default is one) and distr the type of particle init ("1" for init on the x-axis of the grid and "2" for Random).

####C.  OpenMP

If the code (Fortran or python version) was compiled with OpenMP, you can set "x" OpenMP threads per MPI task by defining `export OMP_NUM_THREADS=x` before starting the simulation (default varies with the OS). OpenMP scheduling for balancing loads between tiles can be adjusted at runtime by setting the environment variable `OMP_SCHEDULE` to either `static`, `guided` or `dynamic`. To ensure that threads have enough memory space on the stack, set `OMP_STACKSIZE` to high enough value. In practice, `export OMP_STACKSIZE=32M` should be sufficient for most of test cases.   

**4. Outputs**
-----------------------

####A.  Field diagnostics

For the moment, the code outputs binary matrix files with extensions ".pxr" that can be read using python scripts. Examples of such scripts are in the folder `postproc/`. In the Fortran version, the output frequency is controlled by setting the flag `output_frequency` in the output section of the `input_file.pixr`. Use `output_frequency=-1` to disable outputs. The code places output files in a `RESULTS` directory where the code is ran. This directory has to be created before running the code in your submission script. 

####B.  Temporal diagnostics

Temporal diagnosctics enable to outpout the time evolution of several physical quantities such as the species kinetic energies, the field energies and the L2 norm of divergence of the field minus the charge. Temporal diagnostics can be configurated using the section `temporal`. Output format can be controled using the keyword `format`.
Temporal output files can be binary files (`format=1`) or ascii files (`format=1`).
The output frequency can be controled via the keyword `frequency`.
The different diagnostics can be activated with their flags:

* `kinE=1` for the kinetic energies
* `exE=1`, `eyE=1`, `ezE=1`: electric field energies
* `bxE=1`, `byE=1`, `bzE=1`: magnetic field energies
* `divE-rho=1`: the L2 norm of divergence of the field minus the charge

####C.  Time statistics

Time statistics refer to the computation of the simulation times spend in each significant part of the code. A final time survey is automatically provided at the end of a simulation with the stand-alone code.
The time statisctics function enables to outpout in files the time spend in the main subroutines for each iteration.
It corresponds to the section named `timestat`.

**5. Configuration of the input file**
-----------------------

####A. cpusplit section

This section enables to configure the MPI decomposition and other parameters:

* `nprocx`, `nprocy`, `nprocz`: number of processors in each direction x, y, z
* `topology`: the MPI topology, 0 corresponds to cartesian


# **Particle-In-Cell Scalable Application Resource lite (PICSARlite)**


## **1. Overview**

The Particle-In-Cell Scalable Application Resource lite (PICSARlite) is a proxy of the [PICSAR](http://picsar.net) library.

#### A.  Here are some of the specific algorithmic features of the PICSARlite algorithms:  

* The Maxwell solver use the second-order finite-difference staggered 'Yee' scheme.
* The particle pusher uses the Boris algorithm.
* The field gathering routine is 'energy conserving' (i.e. gathers the field components directly from the staggered 'Yee' mesh).
* The current deposition and field gathering routines include high-order particle shape factors (up to 3rd order).

#### B.  Here are some high performance features of the PICSARlite code :

* Particle tiling to help increase memory locality.
* OpenMP parallelization for intranode parallelism.
* MPI parallelization for internode parallelism (blocking, non-blocking and Remote memory access MPI).
* MPI-IO for fast parallel outputs.
* Minimal vectorization.

## **2. Installation**

Edit the Makefile to adjust compilation variables.
Type `make`.

## **3. Running simulations**

In serial mode, type:
```
PICSARlite/bin/picsar [input filename]
``` 

  - Note: the input filename is optional; the default is 'input_file.pxr'.

In parallel mode, type:
```
mpirun -np n PICSARlite/bin/picsar [input filename]
```

Examples of input files are provided in the `example` subdirectory.

Notice that if `nprocx`, `nprocy` and `nprocz` are provided in the input file as part of the `cpusplit` section, then n must be equal to `nprocx x nprocy x nprocz` with `nprocx`, `nprocy`, `nprocz` the number of processors along x,y,z directions. 

If `nprocx`, `nprocy` and `nprocz` are not defined, the code performs automatic CPU split in each direction. 

User can specify some arguments in the command line. For the moments this feature supports only the number of tiles in each dimension and the init of particle distribution. Ex: `mpirun -np 1 ./picsar -ntilex ntx -ntiley nty -ntilez ntz -distr 1` with `ntx`, `nty` and `ntz` the number of tiles in each dimension (default is one) and distr the type of particle init ("1" for init on the x-axis of the grid and "2" for Random).

##### - Configuration of the input file

Input files are read by PICSAR at the beginning of a run and contain all the required simulation information.

The structure of an input file is a division into sections.
Sections start by `section::name` where `name` is the section name and end with `end::name`.

These sections contain keywords and associated values to be specified.

#### B.  OpenMP

If the code was compiled with OpenMP, you can set "x" OpenMP threads per MPI task by defining `export OMP_NUM_THREADS=x` before starting the simulation (default varies with the OS). OpenMP scheduling for balancing loads between tiles can be adjusted at runtime by setting the environment variable `OMP_SCHEDULE` to either `static`, `guided` or `dynamic`. To ensure that threads have enough memory space on the stack, set `OMP_STACKSIZE` to high enough value. In practice, `export OMP_STACKSIZE=32M` should be sufficient for most of test cases.   

## **5. Outputs**


#### A.  Field diagnostics

For the moment, the code outputs binary matrix files with extensions ".pxr" that can be read using python scripts. Examples of such scripts are in the folder `postproc/`.
(Note: In order for these scripts to work, you need to add the folder `postproc`
to your `$PYTHONPATH`.)

The output frequency is controlled by setting the flag `output_frequency` in the output section of the `input_file.pixr`. Use `output_frequency=-1` to disable outputs. The code places output files in a `RESULTS` directory where the code is ran. This directory has to be created before running the code in your submission script.


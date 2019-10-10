# **Particle-In-Cell Scalable Application Resource (PICSAR)**


## **1. Overview**

The Particle-In-Cell Scalable Application Resource (PICSAR) is a high performance repository intended to help scientists porting their Particle-In-Cell (PIC) codes to the next generation of exascale computers.

PICSAR exploits the three levels of parallelism that will be required to achieve good performances on future architectures: distributed memory parallelization (internode), shared memory parallelization (intranode) and vectorization.

 

PICSAR includes:

- A high performance library of highly optimized versions of the key functionalities of the PIC loop.
- A compact "mini-app" standalone code, which serves as a self-contained proxy that adequately portrays the computational loads and dataflow of more complex PIC codes.
- A Python wrapper for using PICSAR optimized routines with Python-driven codes.

#### A.  Here are some of the specific algorithmic features of the PICSAR algorithms:  

* The Maxwell solver uses arbitrary order finite-difference scheme (staggered/centered).
* The particle pusher uses the Boris or Vay algorithm.
* The field gathering routine is energy conserving.
* The current deposition and field gathering routines include high order particle shape factors (up to 3rd order).

#### B.  Here are some high performance features of the PICSAR code :

* Vectorization.
* Particle tiling to help increase memory locality.
* OpenMP parallelization for intranode parallelism.
* MPI parallelization for internode parallelism (blocking, non-blocking and Remote memory access MPI).
* MPI-IO for fast parallel outputs.

#### C. Build Picsar as Library
* Picsar can be built as a dynamic or static library to provide tools for PIC codes.
To compile picsar as a library you need to swith $(MODE) variable in Makefile to "library", then run make with target lib:
make lib both static and dynamic picsar lib are generated in lib/ file.

#### D.  Python glue:

* We created a Forthon parser that read Fortran source files of PICSAR and parse them to create a `picsar.v` file used by the Forthon compiler to generate a Python module for PICSAR. The Forthon parser is available in the folder `utils`.
* Forthon gives access to all the high performance routines of PICSAR from python.

#### E. PICSARlite

* PICSARlite is a subset of PICSAR that can be used for testing. It is built using a Python script (`utils/generate_miniapp.py`) that accepts a number of options:
	* `--solver`: Maxwell solver method(s) to include [choices: 'all' (default), 'fdtd', 'spectral'].
	* `--pusher`: Particle pusher method(s) to include [choices: 'all' (default), 'Boris', 'Vay'].
	* `--depos`: Type(s) of charge/current deposition to include [choices: 'all' (default), 'direct','Esirkepov'].
	* `--optimization`: Flag to include the optimized versions [choices: 'on' (default), 'off'].

* For example, to create a version with the Maxwell FDTD solver, the Boris pusher, direct charge/current deposition and no optimization, type:
```
python utils/generate_miniapp.py --pusher boris --depos direct --solver fdtd --optimization off
```
* an example input script that runs with this particle configuration is given in 
	```
	picsar/PICSARlite/examples/example_decks_fortran/homogeneous_plasma_lite.pixr
	``` 

	that can be run with 

	```
	mpirun -np 8 ../../fortran_bin/picsar homogeneous_plasma_lite.pixr
	```

## **2. Installation**

Detailed installation instructions are available for various platforms for the FORTRAN and Python libraries in the `PICSAR/Installation_doc` folder.

## **3. Doxygen documentation**

PICSAR offers self-documentation of the source code via Doxygen.

To generate the code html documentation:

  - Download and install [Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html).

  - In the `PICSAR/Doxygen` folder, type:

```
doxygen Doxyfile
```

The html documentation is then accessible in `Doxygen/html/index.html`.

You can also use Doxygen's GUI frontend. Open `Doxygen/Doxyfile`, go to `Run`, click on `Run Doxygen` and finally `Show HTML output`.

For more information on Doxygen, ([see this page](https://www.stack.nl/~dimitri/doxygen/manual/doxygen_usage.html)).


## **4. Running simulations**

PICSAR can be run in two modes:

- In **pure-Fortran mode**: in this case, the code is run as a
  stand-alone application.

- In **Python mode**: in this case, PICSAR is used as a Python module. It can be used
 with existing code (e.g. Warp) to **accelerate** simulations by rerouting the calls
 to the **low-level kernels** (current deposition, field advance, particle pusher).
 More precisely, in this case, instead of calling Warp's regular kernels, the simulation
 will call PICSAR's highly-optimized kernels.

![warp_and_pxr](Doxygen/images/warp_and_picsar.png)

#### A.  Fortran mode

PICSAR input parameters must be provided in an input file named `input_file.pxr`
in the folder where the code is ran. An example (`test.pxr`) of input file is
provided in `examples/example_decks_fortran/`.
To run the executable on n MPI processes:
```
mpirun -np n ./picsar
```
Notice that if `nprocx`, `nprocy` and `nprocz` are provided in the input file as part of the `cpusplit` section, then n must be equal to `nprocx x nprocy x nprocz` with `nprocx`, `nprocy`, `nprocz` the number of processors along x,y,z directions. Otherwise, if `nprocx`, `nprocy` and `nprocz` are not defined, the code performs automatic CPU split in each direction. User can specify some arguments in the command line. For the moments this feature supports only the number of tiles in each dimension and the init of particle distribution. Ex: `mpirun -np 1 ./picsar -ntilex ntx -ntiley nty -ntilez ntz -distr 1` with `ntx`, `nty` and `ntz` the number of tiles in each dimension (default is one) and distr the type of particle init ("1" for init on the x-axis of the grid and "2" for Random).

##### - Configuration of the input file

Input files are read by PICSAR at the beginning of a run and contain all the required simulation information.
They should be in the same repository and renamed `input_file.pixr`.

Examples of input files can be found in the directory `example_decks_fortran`.

The structure of an input file is a division into sections.
Sections start by `section::name` where `name` is the section name and end with `end::name`.
Then these sections contain keywords and values to be specified according to what you want.
In order to learn how to create your own input file and what are the available sections, use the Doxygen documentation.
A page called input file configuration describes the sections and the keywords to set up a correct input file.

#### B.  Python mode

An example of python script `test.py` is provided in `examples/example_scripts_python`.
To run this script in parallel, simply type :
```
mpirun -np NMPI python test.py
```
with NMPI the number of MPI processes.

#### C.  OpenMP

If the code (Fortran or python version) was compiled with OpenMP, you can set "x" OpenMP threads per MPI task by defining `export OMP_NUM_THREADS=x` before starting the simulation (default varies with the OS). OpenMP scheduling for balancing loads between tiles can be adjusted at runtime by setting the environment variable `OMP_SCHEDULE` to either `static`, `guided` or `dynamic`. To ensure that threads have enough memory space on the stack, set `OMP_STACKSIZE` to high enough value. In practice, `export OMP_STACKSIZE=32M` should be sufficient for most of test cases.   

## **5. Outputs**


#### A.  Field diagnostics

For the moment, the code outputs binary matrix files with extensions ".pxr" that can be read using python scripts. Examples of such scripts are in the folder `postproc/`.
(Note: In order for these scripts to work, you need to add the folder `postproc`
to your `$PYTHONPATH`.)

In the Fortran version, the output frequency is controlled by setting the flag `output_frequency` in the output section of the `input_file.pixr`. Use `output_frequency=-1` to disable outputs. The code places output files in a `RESULTS` directory where the code is ran. This directory has to be created before running the code in your submission script.

#### B.  Temporal diagnostics

Temporal diagnosctics enable to outpout the time evolution of several physical quantities such as the species kinetic energies, the field energies and the L2 norm of divergence of the field minus the charge. Temporal diagnostics can be configurated using the section `temporal`. Output format can be controled using the keyword `format`.
Temporal output files can be binary files (`format=1`) or ascii files (`format=1`).
The output frequency can be controled via the keyword `frequency`.
The different diagnostics can be activated with their flags:

* `kinE=1` for the kinetic energies
* `exE=1`, `eyE=1`, `ezE=1`: electric field energies
* `bxE=1`, `byE=1`, `bzE=1`: magnetic field energies
* `divE-rho=1`: the L2 norm of divergence of the field minus the charge

#### C.  Time statistics

Time statistics refer to the computation of the simulation times spend in each significant part of the code. A final time survey is automatically provided at the end of a simulation with the stand-alone code.
The time statisctics function enables to outpout in files the time spend in the main subroutines for each iteration.
It corresponds to the section named `timestat`.

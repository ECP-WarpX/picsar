# **Installation of PICSAR in full Fortran 90 on Cori Phase II**

## Makefile configuration and compilation

### Requirement

Cori phase II is a cluster equipped of Intel Xeon Phi KNL processors.
Specific compiler flags have to be used for this architecture.

We recommend to use the Intel compiler to install PICSAR.
PICSAR also works fine with GNU.

### Compiling with Intel

To compile with Intel in production mode (optimization flags, `MODE=prod`), you can enter:
```
make SYS=cori2 MODE=prod
```

You can compile in debug mode by entering:
```
make SYS=cori2 MODE=debug
```

The makefile will use these compiler flags:
```
-O3 -xMIC-AVX512 -qopenmp -align array64byte -qopt-streaming-stores auto -qopt-report:5
```

Here, `-xMIC-AVX512` specifies that we use AVX512 vector register like on KNL.
This is important for vectorization.


### Other compilers

In this case, you have to edit the file Makefile and indicate
the following environment variables:

* FC: your MPI Fortran compiler wrapper (e.g mpif90, mpiifort, ftn etc.),

* FARGS: your compiler arguments (optimization flags etc.).
To get OpenMP version of PICSAR use the flag -fopenmp (with gfortran) and -openmp (Cray, Intel).
NB: this version of PICSAR requires at least **OpenMP 4.0**.

For instance, to use the GNU compiler:
```
make CC=cc FC=ftn FARGS='-fPIC -O3 -march=knl -ffast-math -funroll-loops'
```

The argument `-march=knl` specifies that we want to compile for KNL
architecture which is important for vectorization.

## Testing the code

To test the code on the cluster, we recommend to use an interactive session
on a single KNL node. In this example, we use a knl configured in quadrant
cache mode (default mode on Cori Phase II).

```
salloc -N 1 -t 00:30:00 -p knl -C knl,quad,cache
```

Then you will have to directly go in the test folders.
Go to `Acceptance_testing/Fortran_tests`.

Then, you can run each test one by one entering the different folders:
* `test_plasma_drift`: drifting beams are propagating in every directions
* `test_homogeneous_plasma`: an homogeneous plasma with a thermal temperature
* `test_Langmuir_wave`: The Langmuir oscillation in a hydrogen plasma

Copy the Picsar executable:

```
cp ../../../fortran_bin/picsar_cori2 .
```

Then, run each test using 4 mpi processes. Use as many OpenMP as you want.
You can have until 64 threads per MPI process. Here we can use 32 threads by doing:

```
export OMP_NUM_THREADS=32
```

Specify the OpenMP scheduler:

```
export OMP_SCHEDULE=dynamic
```

Increase the stack size to avoid unexpected segmentation fault:

```
export OMP_STACKSIZE=128M
```

Then, run the test by doing:

```
srun -n 4 -c 32 ./picsar_cori2
```

## Unit tests

PICSAR also contains unit tests to run and validate specific algorithms
and subroutines of the code.

To compile the unit tests, first clean your installation by doing:
```
make clean
```

You can compile and run each unit test one by one but here we will show you
how to compile and run all of them.

Unit tests are located in the directory `Acceptance_testing/Gcov_tests`.

To compile them with Intel on Cori phase II, enter:
```
make build_test SYS=cori2
```

To run them, go directly in the unit test folder:

```
cd Acceptance_testing/Gcov_tests
```

Use `srun` to run the test of your choice. For instance:
```
srun -n 1 -c 1 ./tile_curr_depo_3d_test
```

Some tests can be run with several OpenMP threads, others need just one MPI
rank and a single thread since it only tests vectorized subroutines:
* current_deposition_3d_test: test 3D classical current deposition subroutines
* esirkepov_2d_test: test 2D Esirkepov current deposition subroutines
* esirkepov_3d_test: test 3D Esirkepov current deposition subroutines
* field_gathering_2d_test: test 2D field gathering subroutines
* field_gathering_test: test 3D field gathering subroutines
* rho_deposition_3d_test: test 3D charge deposition subroutines
* tile_curr_depo_3d_test: test the 3D current deposition with the tiling. (you can use several OpenMP threads for this test)
* tile_field_gathering_3d_test: test the 3D field gathering with the tiling. (you can use several OpenMP threads for this test)
* tile_mpi_part_com_test.F90: test the particle exchanges between tiles and MPI domain boundaries (you can use several OpenMP threads for this test)
* tile_particle_push_3d_test: test the 3D particle pusher with the tiling (you can use several OpenMP threads for this test)
* tile_rho_depo_3d_test: test the 3D charge deposition with the tiling (you can use several OpenMP threads for this test)

For each test, a message will tell you if you pass or fail.

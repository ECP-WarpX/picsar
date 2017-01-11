# **Installation of PICSAR in full Fortran 90 on MacOs**

## Makefile configuration and compilation

### Requirement

You just need MPI and OpenMP to run PICSAR.

### Compiling with Gfortran

To compile with Gfortran in production mode (optimization flags, `MODE=prod`), just enter:
```
make
```
 at the root of the PICSAR directory. The default compiler is the GNU one (`COMP=gnu`).

You can also rapidly compile in debug mode by entering:
```
make MODE=debug
```

### Compiling with Intel

To compile with Intel in production mode (optimization flags, `MODE=prod`), you can enter:
```
make COMP=intel
```

You can compile in debug mode by entering:
```
make COMP=intel MODE=debug
```

### Other compilers

In this case, you have to edit the file Makefile and indicate 
the following environment variables:

* FC: your MPI Fortran compiler wrapper (e.g mpif90, mpiifort, ftn etc.),

* FARGS: your compiler arguments (optimization flags etc.). 
To get OpenMP version of PICSAR use the flag -fopenmp (with gfortran) and -openmp (Cray, Intel). 
NB: this version of PICSAR requires at least **OpenMP 4.0**. 

## Testing the code

You can test the code after compilation in different ways.

To test the compilation/execution, you can use the makefile (py.test is required):

  For all tests:
  > make -f Makefile test

  For each test one by one
  - Drifted plasma:     make -f Makefile test1
  - Homogeneous plasma: make -f Makefile test2
  - Langmuir wave:      make -f Makefile test3    

To test the compilation/execution, you can run test python scripts manually:

* Go the folder Acceptance_testing/Fortran_tests.

* You can perform each test one by one entering the different folders:
  - test_plasma_drift
  - test_homogeneous_plasma
  - test_Langmuir_wave
    
* You can run scripts using py.test: 
  > py.test -s --trun=1 --ttest=1
  
  --trun=0/1: this option enables/disables simulation run
  --ttest=0/1: this option enables/disables assert tests
  
* You can run scripts without py.test:
  > python <pyhthon_script> -r 1 -t 1
  
  -r 0/1: this option enables/disables simulation run
  -t 0/1: this option enables/disables simulation assert tests
  
## Unit tests

PICSAR also contains unit tests to run and validate specific algorithms 
and subroutines of the code.

To compile the unit tests, first clean your installation by doing:
```
make clean
```

You can compile and run each unit test one by one but here we will show you 
how to compile and run all of them.

Unit tests are located in the directory `Acceptance_testing/Gcov_test`.

To compile them with GNU, enter:
```
make build_test
```

To run them, enter:
```
make test_gcov
```

For each test, a message will tell you if you pass or fail.


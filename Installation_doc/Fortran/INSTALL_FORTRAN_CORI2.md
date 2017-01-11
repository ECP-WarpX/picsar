**INSTALLING PICSAR IN FULL FORTRAN 90**
========================================


**Makefile config**
-------------------------

First edit the file Makefile and indicate the following environment variables:

* FC: your MPI Fortran compiler wrapper (e.g mpif90, mpiifort, ftn etc.),

* FARGS: your compiler arguments (optimization flags etc.). To get OpenMP version of PICSAR use the flag -fopenmp (with gfortran) and -openmp (Cray, Intel). NB: this version of PICSAR requires at least **OpenMP 4.0**. 

**Compiling**
-------------------------

To compile, invoke the rule "all": 

* Make -f Makefile all

**Testing**
_________________________

Testing the code after compilation is highly recommended.

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

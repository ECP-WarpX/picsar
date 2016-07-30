FC=mpif90
FARGS= -O3 -fopenmp -g -ftree-vectorize -ftree-vectorizer-verbose=2 
#-ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -fopt-info
#FARGS=-g
SRCDIR= src
BINDIR = fortran_bin
APPNAME=picsar

$(SRCDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FARGS) -c -o $@ $<

all: clean build
test: test1 test2 test3


build:$(SRCDIR)/modules.o $(SRCDIR)/maxwell.o $(SRCDIR)/tiling.o $(SRCDIR)/sorting.o $(SRCDIR)/particles_push_2d.o $(SRCDIR)/particles_push.o $(SRCDIR)/current_deposition_2d.o $(SRCDIR)/current_deposition.o $(SRCDIR)/field_gathering_2d.o $(SRCDIR)/field_gathering.o $(SRCDIR)/mpi_derived_types.o $(SRCDIR)/boundary.o $(SRCDIR)/diags.o $(SRCDIR)/simple_io.o $(SRCDIR)/mpi_routines.o $(SRCDIR)/submain.o $(SRCDIR)/control_file.o $(SRCDIR)/main.o 
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
	
clean:
	rm -rf $(SRCDIR)/*.o *.mod $(BINDIR) RESULTS

test1:
	cd Acceptance_testing/Fortran_tests/test_plasma_drift && \
	py.test -s --ttest=1 --trun=1
test2:
	cd Acceptance_testing/Fortran_tests/test_homogeneous_plasma && \
	py.test -s --ttest=0 --trun=1	
test3:
	cd Acceptance_testing/Fortran_tests/test_Langmuir_wave && \
	py.test -s --ttest=0 --trun=1
# Previous version
#mkdir -p RESULTS
#cp $(BINDIR)/$(APPNAME) RESULTS
#cp example_decks_fortran/test.pixr RESULTS/input_file.pixr
#cd RESULTS && \
#mpirun -np 8 ./picsar -ntilex 10 -ntiley 10 -ntilez 10 -distr 1



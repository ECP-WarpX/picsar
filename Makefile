FC=mpif90
FARGS= -O3 -fopenmp
SRCDIR= src
BINDIR = fortran_bin
APPNAME=picsar


$(SRCDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FARGS) -c -o $@ $<

all: clean build test 


build:$(SRCDIR)/modules.o $(SRCDIR)/maxwell.o $(SRCDIR)/tiling.o $(SRCDIR)/particles_push.o $(SRCDIR)/current_deposition.o $(SRCDIR)/field_gathering.o $(SRCDIR)/mpi_derived_types.o $(SRCDIR)/boundary.o $(SRCDIR)/simple_io.o $(SRCDIR)/diags.o $(SRCDIR)/submain.o $(SRCDIR)/mpi_routines.o $(SRCDIR)/control_file.o  $(SRCDIR)/main.o
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
	
clean:
	rm -rf *.o *.mod $(APPNAME) *.pxr

test:
	mkdir -p RESULTS
	cp $(BINDIR)/$(APPNAME) RESULTS
	cp example_decks_fortran/test.pixr RESULTS/input_file.pixr
	cd RESULTS && \
	mpirun -np 8 ./picsar -ntilex 10 -ntiley 10 -ntilez 10 -distr 1


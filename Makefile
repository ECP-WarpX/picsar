FC=mpif90
FARGS= -O3 -fopenmp
SRCDIR= src
BINDIR = fortran_bin
APPNAME=picsar


$(SRCDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FARGS) -c -o $@ $<

all:$(SRCDIR)/modules.o $(SRCDIR)/maxwell.o $(SRCDIR)/tiling.o $(SRCDIR)/particles_push.o $(SRCDIR)/current_deposition.o $(SRCDIR)/field_gathering.o $(SRCDIR)/mpi_derived_types.o $(SRCDIR)/boundary.o $(SRCDIR)/simple_io.o $(SRCDIR)/diags.o $(SRCDIR)/submain.o $(SRCDIR)/mpi_routines.o $(SRCDIR)/control_file.o  $(SRCDIR)/main.o
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
	
clean:
	rm -rf *.o *.mod $(APPNAME) *.pxr

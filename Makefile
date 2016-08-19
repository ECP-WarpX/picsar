# ________________________________________________________________________________________
#
# PICSAR FORTRAN Makefile
#
# This makefile contains many compiling options and possibilities
#
# Mira/Vesta: FC=mpixlf90_r FARGS='-O3 -qsmp=omp -g'
# Mathieu Lobet, 2016
# ________________________________________________________________________________________

# ________________________________________________________________________________________
# Configuration (user and default)

# Compiler type
# - gnu
# - intel
# - user
COMP=gnu

# Mode
# - prod: production mode
# - debug: debug mode
# - vtune: vtune profiling
# - sde: sde profiling
# - map: Allinea Map profiling
MODE=prod

# System
# - cori2
# - cori1
# - edison
# - default
SYS=default

# Architecture
# - knl
# - ivy
# - hsw
# - host
ARCH=


# User Compiler Parameters (that can be tuned by the user)
# Fortan compiler arguments
FC=mpif90
# C compiler
CC=mpicc
# Fortran compiler arguments
FARGS= -O3 -fopenmp -JModules

# Source directory
SRCDIR=src
# Binary directory
BINDIR=fortran_bin
# Application name
APPNAME=picsar
# Module (.mod) directory
MODDIR=Modules

# ________________________________________________________________________________________
# Preparation of the flags


# GNU compiler ______________________________________________
ifeq ($(COMP),gnu)

  ifeq ($(MODE),prod)
	  FC=mpif90
	  FARGS= -O3 -fopenmp -JModules -ftree-vectorize -ftree-vectorizer-verbose=2
	  #-ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -fopt-info
	  #FARGS=-g
	else ifeq ($(MODE),debug)
	  FC=mpif90
	  FARGS= -O3 -fopenmp -g -JModules -fcheck=bound -ftree-vectorize -ftree-vectorizer-verbose=2	
	endif
	
	# ___ Architecture ________
	ifeq ($(ARCH),hsw)
    ARCH=
  endif
	
# Intel compiler  ______________________________________________
else ifeq ($(COMP),intel)

  # ___ Mode ______________
  ifeq ($(MODE),prod)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -align array64byte -qopt-streaming-stores
	else ifeq ($(MODE),debug)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1
	else ifeq ($(MODE),vtune)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1	
	else ifeq ($(MODE),sde)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1	
	else ifeq ($(MODE),map)			  
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1	
  endif
  
  # ___ Architecture ________
  ifeq ($(ARCH),host)
    LARCH= -xHOST 
  else ifeq ($(ARCH),knl)
    LARCH= -xMIC-AVX512
  else ifeq ($(ARCH),ivy)
    LARCH= -xAVX  
  else ifeq ($(ARCH),hsw)   
    LARCH= -xCORE-AVX2   
  endif
  
endif

FARGS+= $(LARCH)

# ________________________________________________________
# Not used for the moment
FSOURCE= $(wildcard $(SRCDIR)/*.F90)
FOBJS= $(FSOURCE:.F90=.o)
FDEPT= $(FSOURCE:.F90=.d)
-include $(FDEPT)
# ________________________________________________________

$(SRCDIR)/%.o $(SRCDIR)/%.mod:$(SRCDIR)/%.F90
	$(FC) -c $(FARGS) -o $@ $<

all: clean echo createdir build
test: test1 test2 test3

build:$(SRCDIR)/modules.o $(SRCDIR)/maxwell.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/sorting.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/diags.o \
	$(SRCDIR)/simple_io.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/control_file.o \
	$(SRCDIR)/main.o 
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
	
clean:
	rm -rf $(SRCDIR)/*.o *.mod $(MODDIR)/*.mod $(BINDIR) RESULTS Acceptance_testing/Gcov_tests/*.o

createdir:
	mkdir -p $(MODDIR)

echo:
	@echo ' Compiler $(COMP)'
	@echo ' Fortran wrapper $(FC)'
	@echo ' Fortran arguments $(FARGS)'

# Compiler type
# - gnu
# - intel
# - user
COMP=gnu

help:
	@echo ' ______________________________________ '
	@echo ' Makefile information'
	@echo 
	@echo ' Targets:'
	@echo ' - build'
	@echo ' - clean'
	@echo ' - buildtest'
	@echo ' - test_gcov'
	@echo	
	@echo ' COMP= Compiler type:'
	@echo ' - gnu: gnu compiler'
	@echo ' - intel: intel compiler'
	@echo ' - user: user defined makefile'	
	@echo
	@echo ' MODE= Compilation mode'
	@echo ' - prod: production compilation'
	@echo ' - debug: debug compilation'
	@echo ' - vtune: vtune analysis'	
	@echo ' ______________________________________ '
		
# ________________________________________________________________________________________
# make tests

Acceptance_testing/Gcov_tests/%.o:Acceptance_testing/Gcov_tests/%.F90
	$(FC) -c $(FARGS) -o $@ $<
	
cleantest:
	rm Acceptance_testing/Gcov_tests/*.o
	
buildtest: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering.o \
	Acceptance_testing/Gcov_tests/field_gathering_test.o \
	Acceptance_testing/Gcov_tests/current_deposition_3d_test.o \
	Acceptance_testing/Gcov_tests/esirkepov_3d_test.o \
	Acceptance_testing/Gcov_tests/esirkepov_2d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/field_gathering_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/field_gathering_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/current_deposition_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/current_deposition_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/esirkepov_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/esirkepov_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/esirkepov_2d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/esirkepov_2d_test.o
#	$(FC) -g -O0 -ftest-coverage -JModules -o Acceptance_testing/Gcov_tests/field_gathering_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/field_gathering_test.o			

test1:
	cd Acceptance_testing/Fortran_tests/test_plasma_drift && \
	py.test -s --ttest=1 --trun=1
	
test2:
	cd Acceptance_testing/Fortran_tests/test_homogeneous_plasma && \
	py.test -s --ttest=0 --trun=1	
	
test3:
	cd Acceptance_testing/Fortran_tests/test_Langmuir_wave && \
	py.test -s --ttest=0 --trun=1

test_physics:
	test1 test2 test3
	
test_gcov:
	./Acceptance_testing/Gcov_tests/field_gathering_3d_test
	./Acceptance_testing/Gcov_tests/current_deposition_3d_test
	./Acceptance_testing/Gcov_tests/esirkepov_3d_test
	./Acceptance_testing/Gcov_tests/esirkepov_2d_test

esirkepov_3d_test:	
	./Acceptance_testing/Gcov_tests/esirkepov_3d_test	
esirkepov_2d_test:	
	./Acceptance_testing/Gcov_tests/esirkepov_2d_test
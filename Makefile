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

# Compiler type (COMP)
# - gnu
# - intel
# - user
COMP=gnu

# Mode (MODE)
# - prod: production mode
# - debug: debug mode
# - vtune: vtune profiling
# - sde: sde profiling
# - map: Allinea Map profiling
MODE=prod

# System (SYS)
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

# If a system is specified
# Cori phase 1
ifeq ($(SYS),cori1)
	FC=ftn
	CC=cc
	APPNAME=picsar_cori
  ifeq ($(MODE),prod)
		COMP=none
		FARGS= -O3 -xCORE-AVX2 -qopenmp -align array64byte -qopt-streaming-stores auto
		# -qopt-report:5
		LARCH=
	else ifeq ($(MODE),debug)
		COMP=none
		FARGS= -g -O3 -xCORE-AVX2 -qopenmp -qopt-report:5 -debug inline-debug-info
		LARCH=
	else ifeq ($(MODE),vtune)
		APPNAME=picsar_cori_vtune
		COMP=none
		FARGS= -D VTUNE=1 -O3 -g -dynamic -debug inline-debug-info -qopenmp -xCORE-AVX2 -align array64byte
		CARGS= -D VTUNE=1 -O3 -g -dynamic -qopenmp -xCORE-AVX2 -I $(VTUNE_AMPLIFIER_XE_2016_DIR)/include
		LDFLAGS= $(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
		LARCH= 	
	else ifeq ($(MODE),novec)
		COMP=none
		FARGS= -g -O0 -no-simd -no-vec
		LARCH=	
	endif
# Edison
else ifeq ($(SYS),edison)
	FC=ftn
	CC=cc
	APPNAME=picsar_edison
  ifeq ($(MODE),prod)
		COMP=none
		FARGS= -O3 -xAVX -qopenmp -align array64byte -qopt-streaming-stores auto
		# -qopt-report:5
		LARCH=
	else ifeq ($(MODE),debug)
		COMP=none
		FARGS= -g -O3 -xAVX -qopenmp -qopt-report:5 -debug inline-debug-info -traceback
		LARCH=	
	else ifeq ($(MODE),novec)
		COMP=none
		FARGS= -g -O0 -no-simd -no-vec
		LARCH=	
	endif
# Cori phase 2 at NERSC
else ifeq ($(SYS),cori2)
	FC=ftn
	CC=cc
	APPNAME=picsar_cori2
  ifeq ($(MODE),prod)
		COMP=none
		FARGS= -O3 -xMIC-AVX512 -qopenmp -align array64byte -qopt-streaming-stores auto -qopt-report:5
		LARCH=
	else ifeq ($(MODE),debug)
		COMP=none
		FARGS= -g -O3 -xMIC-AVX512 -qopenmp -debug inline-debug-info -traceback -qopt-report:5
		LARCH=	
	else ifeq ($(MODE),novec)
		COMP=none
		FARGS= -g -O0 -no-simd -no-vec
		LARCH=	
	endif
# ___ Carl KNL whitebox at NERSC _____________________
else ifeq ($(SYS),carl)
	FC=mpiifort
	CC=icc
	APPNAME=picsar_carl
  ifeq ($(MODE),prod)
		COMP=none
		FARGS= -O3 -xMIC-AVX512 -qopenmp -align array64byte -qopt-streaming-stores auto -qopt-report:5
		LARCH=
	else ifeq ($(MODE),debug)
		COMP=none
		FARGS= -g -O3 -xMIC-AVX512 -qopenmp -debug inline-debug-info -traceback -qopt-report:5
		LARCH=
	else ifeq ($(MODE),vtune)
		APPNAME=picsar_carl_vtune
		COMP=none
		FARGS= -D VTUNE=1	-g -dynamic -O3 -xMIC-AVX512 -qopenmp -debug inline-debug-info -qopt-streaming-stores auto
		CARGS= -D VTUNE=1 -g -dynamic -O3 -qopenmp -xMIC-AVX512 -I $(VTUNE_AMPLIFIER_XE_2016_DIR)/include
		LDFLAGS= $(VTUNE_AMPLIFIER_XE_2016_DIR)/lib64/libittnotify.a
		LARCH= 	
	else ifeq ($(MODE),advisor)
		APPNAME=picsar_carl_advisor
		COMP=none
		FARGS= -g -O3 -xMIC-AVX512 -qopenmp -Bdynamic -debug inline-debug-info -align array64byte -qopt-streaming-stores auto
		LARCH=	
	else ifeq ($(MODE),novec)
		APPNAME=picsar_carl_novec
		COMP=none
		FARGS= -D NOVEC=0 -g -O3 -xMIC-AVX512 -qopenmp -no-simd -no-vec  -align array64byte -qopt-streaming-stores auto
		LARCH=
	else ifeq ($(MODE),nofma)		
		APPNAME=picsar_carl_nofma
		COMP=none
		FARGS=  -g -O3 -xMIC-AVX512 -qopenmp -no-fma  -align array64byte -qopt-streaming-stores auto
		LARCH=
	endif
endif

# GNU compiler ______________________________________________
ifeq ($(COMP),gnu)

  ifeq ($(MODE),prod)
	  FC=mpif90
	  FARGS= -O3 -fopenmp -JModules -ftree-vectorize 
	  #-ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -fopt-info
	  #FARGS=-g
	else ifeq ($(MODE),debug)
	  FC=mpif90
	  FARGS= -O3 -fopenmp -g -JModules -fcheck=bound -ftree-vectorize
	else ifeq ($(MODE),novec)
	  FC=mpif90
	  FARGS= -D NOVEC=0 -O3 -fopenmp -JModules
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
	  FARGS= -O3 -qopenmp -JModules -align array64byte -qopt-streaming-stores auto -qopt-report:5
	else ifeq ($(MODE),debug)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1 -qopt-report:5
	else ifeq ($(MODE),vtune)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1 -qopt-report:5
	else ifeq ($(MODE),sde)
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1 -qopt-report:5
	else ifeq ($(MODE),map)			  
	  FC=mpif90
	  FARGS= -O3 -qopenmp -JModules -check bounds -D DEBUG=1 -qopt-report:5
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
#FSOURCE= $(wildcard $(SRCDIR)/*.F90)
#FOBJS= $(FSOURCE:.F90=.o)
#FDEPT= $(FSOURCE:.F90=.d)
#-include $(FDEPT)
# ________________________________________________________


$(SRCDIR)/%.o $(SRCDIR)/%.mod $(MODDIR)/%.mod:$(SRCDIR)/%.F90
	$(FC) $(FARGS) -c -o $@ $<

$(SRCDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CARGS) -c -o $@ $<

all: echo createdir build
test: test1 test2 test3

ifeq ($(MODE),vtune)
build:$(SRCDIR)/modules.o \
	$(SRCDIR)/api_fortran_itt.o \
	$(SRCDIR)/itt_fortran.o \
	$(SRCDIR)/maxwell.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/sorting.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/charge_deposition.o \
	$(SRCDIR)/diags.o \
	$(SRCDIR)/simple_io.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/control_file.o \
	$(SRCDIR)/main.o 
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o $(LDFLAGS)
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
else ifeq ($(MODE),sde)
build:$(SRCDIR)/modules.o \
	$(SRCDIR)/api_fortran_sde.o \
	$(SRCDIR)/sde_fortran.o \
	$(SRCDIR)/maxwell.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/sorting.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/charge_deposition.o \
	$(SRCDIR)/diags.o \
	$(SRCDIR)/simple_io.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/control_file.o \
	$(SRCDIR)/main.o 
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o $(LDFLAGS)
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
else
build:$(SRCDIR)/modules.o \
	$(SRCDIR)/maxwell.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/sorting.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/charge_deposition.o \
	$(SRCDIR)/diags.o \
	$(SRCDIR)/simple_io.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/control_file.o \
	$(SRCDIR)/main.o 
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
endif
	
clean: clean_test
	rm -rf $(SRCDIR)/*.o
	rm -f *.mod
	rm -f $(BINDIR)/$(APPNAME)
	rm -rf RESULTS
	rm -r $(MODDIR)
	rm -f $(SRCDIR)/*.mod
	rm -rf *.dSYM
	rm -f Doxygen/*.tmp
	
createdir:
	mkdir -p $(MODDIR)

echo:
	@echo ' Compiler $(COMP)'
	@echo ' MPI wrapper $(FC)'
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
	@echo ' - cleantest' 
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
	@echo ' - novec: disable vectorization'	
	@echo ' - vtune: vtune analysis'	
	@echo
	@echo ' SYS= System'		
	@echo ' - edison: Edison NERSC'	
	@echo ' - carl:   NERSC KNL whitebox'
	@echo ' - cori1:  Cori phase 1 NERSC'	
	@echo ' - cori2:  Cori phase 2 NERSC'
	@echo ' ______________________________________ '
		
# ________________________________________________________________________________________
# make tests

# To be used to compile the test codes
Acceptance_testing/Gcov_tests/%.o:Acceptance_testing/Gcov_tests/%.F90
	$(FC) -c $(FARGS) -o $@ $<

# Clean files related to the tests	
clean_test:
	rm -f Acceptance_testing/Fortran_tests/*/picsar
	rm -rf Acceptance_testing/Fortran_tests/*/RESULTS	
	rm -f Acceptance_testing/Python_tests/*/*.cgm
	rm -f Acceptance_testing/Python_tests/*/*.cgmlog
	rm -f Acceptance_testing/Gcov_tests/*.o
	rm -f Acceptance_testing/Gcov_tests/*_test
	rm -rf Acceptance_testing/Gcov_tests/*.dSYM

build_tile_field_gathering_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/control_file.o \
	Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test.o	

build_field_gathering_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	Acceptance_testing/Gcov_tests/field_gathering_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/field_gathering_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/field_gathering_test.o	

build_field_gathering_2d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	Acceptance_testing/Gcov_tests/field_gathering_2d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/field_gathering_2d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/field_gathering_2d_test.o

build_current_deposition_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	Acceptance_testing/Gcov_tests/current_deposition_3d_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/current_deposition_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/current_deposition_3d_test.o

build_tile_particle_push_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/control_file.o \
	Acceptance_testing/Gcov_tests/tile_particle_push_3d_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_particle_push_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/tile_particle_push_3d_test.o

build_tile_mpi_part_com_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/particles_push_2d.o \
	$(SRCDIR)/particles_push.o \
	$(SRCDIR)/field_gathering_2d.o \
	$(SRCDIR)/field_gathering_3d_o1.o \
	$(SRCDIR)/field_gathering_3d_o2.o \
	$(SRCDIR)/field_gathering_3d_o3.o \
	$(SRCDIR)/field_gathering.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/control_file.o \
	Acceptance_testing/Gcov_tests/tile_mpi_part_com_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_mpi_part_com_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/tile_mpi_part_com_test.o

build_rho_deposition_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/charge_deposition.o \
	Acceptance_testing/Gcov_tests/rho_deposition_3d_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/rho_deposition_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/rho_deposition_3d_test.o
	
build_tile_rho_depo_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/charge_deposition.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/control_file.o \
	Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test.o

build_tile_curr_depo_3d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	$(SRCDIR)/mpi_derived_types.o \
	$(SRCDIR)/boundary.o \
	$(SRCDIR)/mpi_routines.o \
	$(SRCDIR)/control_file.o \
	Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test.o 
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test.o	

build_esirkepov_3d_test:$(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	Acceptance_testing/Gcov_tests/esirkepov_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/esirkepov_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/esirkepov_3d_test.o
	
build_esirkepov_2d_test: $(SRCDIR)/modules.o \
	$(SRCDIR)/tiling.o \
	$(SRCDIR)/current_deposition_2d.o \
	$(SRCDIR)/current_deposition.o \
	Acceptance_testing/Gcov_tests/esirkepov_2d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/esirkepov_2d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/esirkepov_2d_test.o
		
# Compilation of all the tests	
build_test: createdir \
	build_tile_field_gathering_3d_test \
	build_field_gathering_3d_test \
	build_field_gathering_2d_test \
	build_rho_deposition_3d_test \
	build_current_deposition_3d_test \
	build_tile_particle_push_3d_test \
	build_tile_rho_depo_3d_test \
	build_tile_curr_depo_3d_test \
	build_esirkepov_3d_test \
	build_esirkepov_2d_test 

#	$(FC) -g -O0 -ftest-coverage -JModules -o Acceptance_testing/Gcov_tests/field_gathering_3d_test $(SRCDIR)/*.o Acceptance_testing/Gcov_tests/field_gathering_test.o			

# __ Execute Pytest ____________________________________________________
test1:
	cd Acceptance_testing/Fortran_tests/test_plasma_drift && \
	py.test -s --ttest=1 --trun=1
	
test2:
	cd Acceptance_testing/Fortran_tests/test_homogeneous_plasma && \
	py.test -s --ttest=0 --trun=1	
	
test3:
	cd Acceptance_testing/Fortran_tests/test_Langmuir_wave && \
	py.test -s --ttest=0 --trun=1

test_pytest:
	test1 test2 test3

# __ Execute Gcov test ____________________________________________________
	
test_gcov: field_gathering_2d_test \
	field_gathering_3d_test \
	field_gathering_2d_test \
	rho_deposition_3d_test \
	tile_field_gathering_3d_test \
	tile_particle_push_3d_test \
	tile_curr_depo_3d_test \
	tile_rho_depo_3d_test \
	current_deposition_3d_test \
	esirkepov_3d_test \
	esirkepov_2d_test

current_deposition_3d_test:
	export OMP_NUM_THREADS=1
	./Acceptance_testing/Gcov_tests/current_deposition_3d_test
	
field_gathering_2d_test:
	export OMP_NUM_THREADS=1
	./Acceptance_testing/Gcov_tests/field_gathering_2d_test
	
field_gathering_3d_test:
	export OMP_NUM_THREADS=1
	mpirun -n 1 ./Acceptance_testing/Gcov_tests/field_gathering_3d_test
	
esirkepov_2d_test:	
	export OMP_NUM_THREADS=1
	./Acceptance_testing/Gcov_tests/esirkepov_2d_test
	
esirkepov_3d_test:
	export OMP_NUM_THREADS=1	
	./Acceptance_testing/Gcov_tests/esirkepov_3d_test	

rho_deposition_3d_test:
	export OMP_NUM_THREADS=1
	mpirun -n 1 ./Acceptance_testing/Gcov_tests/rho_deposition_3d_test

tile_field_gathering_3d_test:
	export OMP_NUM_THREADS=4
	mpirun -n 1 ./Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test 

tile_particle_push_3d_test:
	export OMP_NUM_THREADS=4
	mpirun -n 1 ./Acceptance_testing/Gcov_tests/tile_particle_push_3d_test
	
tile_mpi_part_com_test:
	export OMP_NUM_THREADS=2
	mpirun -n 4 ./Acceptance_testing/Gcov_tests/tile_mpi_part_com_test

tile_curr_depo_3d_test:
	export OMP_NUM_THREADS=4
	mpirun -n 1 ./Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test
	
tile_rho_depo_3d_test:
	export OMP_NUM_THREADS=4
	mpirun -n 1 ./Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test
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
# - fujitsu
COMP=gnu

# Mode (MODE)
# - prod: production mode (without FFTW)
# - prod_spectral: production mode with spectral solver and FFTW
# - debug: debug mode
# - vtune: vtune profiling
# - sde: sde profiling
# - map: Allinea Map profiling
# - library: create static and dynamic library
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
# - a64fx (only supported for fujitsu compiler)
ARCH=


# User Compiler Parameters (that can be tuned by the user)
# Fortan compiler arguments
FC=mpif90
# C compiler
CC=mpicc
# C++ compiler
CPP=mpic++

# Fortran compiler arguments
FARGS= -g -fbounds-check -O3 -fopenmp -JModules
# C compiler Arguments
CARGS= -g -O3 -fopenmp
# C++ compiler Arguments
CPPARGS= -g -O3 -fopenmp

# External libs
FFTW3_LIB=/usr/lib/x86_64-linux-gnu
FFTW3_INCLUDE=/usr/include
VTUNEDIR=/opt/intel/vtune_amplifier_xe_2017.2.0.499904

P3DFFT_INCLUDE=
P3DFFT_LIB=
IS_P3DFFT = false



# Source directory
SRCDIR=src
# Binary directory
BINDIR=fortran_bin
# Lib directory
LIBDIR=lib
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
	APPNAME=picsar_cori1
  ifeq ($(MODE),prod)
		APPNAME=picsar_cori1
		COMP=none
		FARGS= -O3 -xCORE-AVX2 -qopenmp -align array64byte -qopt-streaming-stores auto
		# -qopt-report:5
		LARCH=
	else ifeq ($(MODE),debug)
		APPNAME=picsar_cori1_debug
		COMP=none
		FARGS= -g -O3 -xCORE-AVX2 -qopenmp -qopt-report:5 -debug inline-debug-info
		LARCH=
	else ifeq ($(MODE),vtune)
		APPNAME=picsar_cori1_vtune
		COMP=none
		FARGS= -D VTUNE=1 -O3 -g -dynamic -debug inline-debug-info -qopenmp -xCORE-AVX2 -align array64byte
		CARGS= -D VTUNE=1 -O3 -g -dynamic -qopenmp -xCORE-AVX2 -I $(VTUNEDIR)/include
		LDFLAGS= $(VTUNEDIR)/lib64/libittnotify.a
		LARCH=
	else ifeq ($(MODE),sde)
		APPNAME=picsar_cori1_sde
		COMP=none
		FARGS= -D SDE=1	-g -O3 -xCORE-AVX2  -qopenmp -debug inline-debug-info -qopt-streaming-stores auto
		CARGS= -D SDE=1 -g -O3 -qopenmp -xCORE-AVX2
		LARCH=
	else ifeq ($(MODE),advisor)
		APPNAME=picsar_cori1_advisor
		COMP=none
		FARGS= -g -O3 -xCORE-AVX2 -qopenmp -dynamic -debug inline-debug-info -align array64byte
		#-qopt-streaming-stores auto
		LARCH=
	else ifeq ($(MODE),novec)
		APPNAME=picsar_cori1_novec
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
		FARGS= -O3 -xAVX -align array64byte -qopt-streaming-stores auto
		# -qopt-report:5
		LARCH=
	else ifeq ($(MODE),debug)
		APPNAME=picsar_edison_debug
		COMP=none
		FARGS= -g -O3 -xAVX -qopt-report:5 -debug inline-debug-info -traceback
		LARCH=
	else ifeq ($(MODE),sde)
		APPNAME=picsar_edison_sde
		COMP=none
		FARGS= -D SDE=1	-g -O3 -xAVX  -qopenmp -debug inline-debug-info -qopt-streaming-stores auto
		CARGS= -D SDE=1 -g -O3 -qopenmp -xAVX
		LARCH=
	else ifeq ($(MODE),novec)
		APPNAME=picsar_edison_novec
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
  	else ifeq ($(MODE),prod_spectral)
		COMP=none
		FARGS= -O3 -xMIC-AVX512 -qopenmp -align array64byte -qopt-streaming-stores auto -qopt-report:5
		LARCH=
  	else ifeq ($(MODE),debug_spectral)
		APPNAME=picsar_cori2_debug
		COMP=none
		FARGS= -g -O3 -D DEBUG=0 -xMIC-AVX512 -qopenmp -debug inline-debug-info -traceback
		LARCH=
	else ifeq ($(MODE),debug)
		APPNAME=picsar_cori2_debug
		COMP=none
		FARGS= -g -O3 -D DEBUG=0 -xMIC-AVX512 -qopenmp -debug inline-debug-info -traceback -qopt-report:5
		LARCH=
	else ifeq ($(MODE),dev)
		COMP=none
		FARGS= -O3 -D DEV=0 -xMIC-AVX512 -qopenmp -align array64byte
		# -qopt-streaming-stores auto
		LARCH=
	else ifeq ($(MODE),vtune)
		APPNAME=picsar_cori2_vtune
		COMP=none
		#FARGS= -D VTUNE=1 -O3 -g -Bdynamic -qopenmp -xMIC-AVX512 -fPIE -fPIC -align array64byte -debug inline-debug-info
		#CARGS= -D VTUNE=1 -O3 -g -Bdynamic -qopenmp -xMIC-AVX512 -fPIE -fPIC -debug inline-debug-info -I /opt/intel/vtune_amplifier_xe_2017.1.0.486011/include
		#LDFLAGS= /opt/intel/vtune_amplifier_xe_2017.1.0.486011/lib64/libittnotify.a -pie
		FARGS= -D VTUNE=1 -O3 -g -dynamic -qopenmp -xMIC-AVX512 -align array64byte -debug inline-debug-info
		CARGS= -D VTUNE=1 -O3 -g -dynamic -qopenmp -xMIC-AVX512 -debug inline-debug-info -I $(VTUNEDIR)/include
		LDFLAGS= $(VTUNEDIR)/lib64/libittnotify.a
		LARCH=
	else ifeq ($(MODE),sde)
		APPNAME=picsar_cori2_sde
		COMP=none
		FARGS= -D SDE=1 -g -O3 -xMIC-AVX512 -qopenmp -debug inline-debug-info -align array64byte
		#-qopt-streaming-stores auto
		CARGS= -D SDE=1 -g -O3 -qopenmp -xMIC-AVX512
		LARCH=
	else ifeq ($(MODE),advisor)
		APPNAME=picsar_cori2_advisor
		COMP=none
		FARGS= -g -O3 -xMIC-AVX512 -qopenmp -dynamic -debug inline-debug-info -align array64byte
		#-qopt-streaming-stores auto
		LARCH=
	else ifeq ($(MODE),novec)
		APPNAME=picsar_cori2_novec
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
		FARGS= -O3 -xMIC-AVX512 -qopenmp -align array64byte -qopt-streaming-stores auto
		LARCH=
	else ifeq ($(MODE),debug)
		APPNAME=picsar_carl_debug
		COMP=none
		FARGS= -g -O3 -D DEBUG=1 -xMIC-AVX512 -qopenmp -debug inline-debug-info -heap-arrays -fp-stack-check -traceback -qopt-report:5
		LARCH=
	else ifeq ($(MODE),vtune)
		APPNAME=picsar_carl_vtune
		COMP=none
		FARGS= -D VTUNE=1	-g -Bdynamic -O3 -xMIC-AVX512 -qopenmp -debug inline-debug-info -qopt-streaming-stores auto
		CARGS= -D VTUNE=1 -g -Bdynamic -O3 -qopenmp -xMIC-AVX512 -I $(VTUNEDIR)/include
		LDFLAGS= $(VTUNEDIR)/lib64/libittnotify.a
		LARCH=
	else ifeq ($(MODE),sde)
		APPNAME=picsar_carl_sde
		COMP=none
		FARGS= -D SDE=1	-g -O3 -xMIC-AVX512 -qopenmp -debug inline-debug-info -qopt-streaming-stores auto
		CARGS= -D SDE=1 -g -O3 -qopenmp -xMIC-AVX512
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
	  FARGS= -O3 -fopenmp -g -JModules -Wunused-variable -fcheck=bound -ftree-vectorize
	else ifeq ($(MODE),prod_spectral)
	  FC=mpif90
	  FARGS= -O3 -fopenmp -JModules -ftree-vectorize
	else ifeq ($(MODE),debug_spectral)
	  FC=mpif90
	  FARGS= -O3 -fopenmp -JModules -Wunused-variable -ftree-vectorize
	else ifeq ($(MODE),dev)
	  FC=mpif90
	  FARGS= -O3 -D DEV=1 -fopenmp -JModules -ftree-vectorize
	  #-ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -fopt-info
	  #FARGS=-g
	else ifeq ($(MODE),devdebug)
	  FC=mpif90
	  FARGS= -O3 -D DEV=1 -fopenmp -g -JModules  -Wunused-variable -fcheck=bound -ftree-vectorize
	  #-ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -fopt-info
	  #FARGS=-g
	else ifeq ($(MODE),novec)
	  FC=mpif90
	  FARGS= -D NOVEC=0 -O3 -fopenmp -JModules
        else ifeq ($(MODE),library)
          FC=mpif90
          FARGS= -O3 -fopenmp -JModules -ftree-vectorize -fPIC
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

  # Fujitsu compiler  ______________________________________________
else ifeq ($(COMP),fujitsu)

  # ___ Mode ______________
  ifeq ($(MODE),prod)
	  FC=mpifrtpx
	  FARGS= -O3 -Kfast -Kopenmp -Kopenmp_simd -Ksimd=2 -Nlst=t -Koptmsg=2
  endif

    # ___ Architecture ________
  ifeq ($(ARCH),a64fx)
    LARCH= -KA64FX -KSVE -SSL2
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

ifeq ($(MODE),$(filter $(MODE),prod_spectral debug_spectral))
	FARGS += -I$(FFTW3_INCLUDE) -D FFTW=1
	LDFLAGS += -L$(FFTW3_LIB) -lfftw3_mpi -lfftw3  -lfftw3_omp
endif
ifeq ($(IS_P3DFFT),true)
        FARGS += -I$(P3DFFT_INCLUDE)  -D P3DFFT
        LDFLAGS += $(P3DFFT_LIB)/libp3dfft.a
endif

ifeq ($(MODE),library)
        FARGS += -fPIC -I$(FFTW3_INCLUDE) -D LIBRARY=1  -D  FFTW=1
        LDFLAGS += -L$(FFTW3_LIB) -lfftw3_mpi -lfftw3  -lfftw3_omp
endif

$(SRCDIR)/%.o $(SRCDIR)/*/%.o $(SRCDIR)/*/*/%.o $(SRCDIR)/*/*/*/%.o $(SRCDIR)/%.mod $(MODDIR)/%.mod:$(SRCDIR)/%.F90
	$(FC) $(FARGS) -c -o $@ $<
$(SRCDIR)/profiling/%.o:$(SRCDIR)/profiling/%.c
	$(CC) $(CARGS) -c -o $@ $<

all: echo createdir build
test: test1 test2 test3
lib: echo createdir build_lib
build_lib:$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/GPSTD.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fastfft.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/load_balancing.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/init_external.o
	ar rcs libpxr.a $(SRCDIR)/*.o $(SRCDIR)/*/*.o  $(SRCDIR)/*/*/*.o $(SRCDIR)/*/*/*/*.o
	$(FC) $(FARGS) -shared -o libpxr.so $(SRCDIR)/*.o  $(SRCDIR)/*/*.o $(SRCDIR)/*/*/*.o  $(SRCDIR)/*/*/*/*.o
	mv libpxr.a $(LIBDIR)
	mv libpxr.so $(LIBDIR)


ifeq ($(MODE),vtune)
build:$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/profiling/api_fortran_itt.o \
	$(SRCDIR)/profiling/itt_fortran.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_circ.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_circ.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	$(SRCDIR)/main.o
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o $(SRCDIR)/*/*.o $(SRCDIR)/*/*/*.o $(SRCDIR)/*/*/*/*.o $(LDFLAGS)
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
else ifeq ($(MODE),sde)
build:$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/profiling/api_fortran_sde.o \
	$(SRCDIR)/profiling/sde_fortran.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_circ.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_circ.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	$(SRCDIR)/main.o
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o $(SRCDIR)/*/*.o $(SRCDIR)/*/*/*.o $(SRCDIR)/*/*/*/*.o $(LDFLAGS)
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
else ifeq ($(MODE),$(filter $(MODE),prod_spectral debug_spectral))
build:$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fastfft.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/GPSTD.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/kin_energy.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_circ.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_circ.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
        $(SRCDIR)/housekeeping/load_balancing.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	$(SRCDIR)/main.o
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o $(SRCDIR)/*/*.o $(SRCDIR)/*/*/*.o $(SRCDIR)/*/*/*/*.o $(LDFLAGS)
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
else
build:$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/kin_energy.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_circ.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_circ.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	$(SRCDIR)/main.o
	$(FC) $(FARGS) -o $(APPNAME) $(SRCDIR)/*.o $(SRCDIR)/*/*.o $(SRCDIR)/*/*/*.o $(SRCDIR)/*/*/*/*.o
	mkdir -p $(BINDIR)
	mv $(APPNAME) $(BINDIR)
endif

clean: clean_test
	rm -rf $(SRCDIR)/*.o
	rm -rf $(SRCDIR)/*/*.o
	rm -rf $(SRCDIR)/*/*/*.o
	rm -rf $(SRCDIR)/*/*/*/*.o
	rm -f *.mod
	rm -f $(BINDIR)/$(APPNAME)
	rm -rf RESULTS
	rm -rf $(MODDIR)
	rm -f $(SRCDIR)/*.mod
	rm -rf *.dSYM
	rm -f Doxygen/*.tmp

	rm -rf $(LIBDIR)/*

createdir:
	mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)

echo:
	@echo	''
	@echo ' MPI wrapper $(FC)'
	@echo ' Fortran arguments $(FARGS)'
	@echo	''

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
	@echo ' - build_test'
	@echo ' - clean_test'
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
	@echo ' - vtune: Vtune analysis'
	@echo ' - sde: Intel SDE analysis'
	@echo ' - advisor: Intel Advisor analysis'
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

build_tile_field_gathering_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_field_gathering_3d_test.o

build_field_gathering_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	Acceptance_testing/Gcov_tests/field_gathering_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/field_gathering_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	Acceptance_testing/Gcov_tests/field_gathering_test.o

build_field_gathering_2d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	Acceptance_testing/Gcov_tests/field_gathering_2d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/field_gathering_2d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	Acceptance_testing/Gcov_tests/field_gathering_2d_test.o

build_current_deposition_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	Acceptance_testing/Gcov_tests/current_deposition_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/current_deposition_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	Acceptance_testing/Gcov_tests/current_deposition_3d_test.o

build_tile_particle_push_3d_test: createdir \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_particle_push_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_particle_push_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_particle_push_3d_test.o

build_tile_mpi_part_com_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_mpi_part_com_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_mpi_part_com_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_mpi_part_com_test.o

build_rho_deposition_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	Acceptance_testing/Gcov_tests/rho_deposition_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/rho_deposition_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	Acceptance_testing/Gcov_tests/rho_deposition_3d_test.o

build_tile_rho_depo_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_rho_depo_3d_test.o

build_tile_curr_depo_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/tile_curr_depo_3d_test.o

build_esirkepov_3d_test:$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	Acceptance_testing/Gcov_tests/esirkepov_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/esirkepov_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	Acceptance_testing/Gcov_tests/esirkepov_3d_test.o

build_esirkepov_2d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	Acceptance_testing/Gcov_tests/esirkepov_2d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/esirkepov_2d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	Acceptance_testing/Gcov_tests/esirkepov_2d_test.o


build_maxwell_2d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fastfft.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/GPSTD.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/kin_energy.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/housekeeping/load_balancing.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/maxwell_2d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/maxwell_2d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fastfft.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/GPSTD.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/kin_energy.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/housekeeping/load_balancing.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/maxwell_2d_test.o $(LDFLAGS)

build_maxwell_3d_test: $(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fastfft.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/GPSTD.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/kin_energy.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/housekeeping/load_balancing.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/maxwell_3d_test.o
	$(FC) $(FARGS) -o Acceptance_testing/Gcov_tests/maxwell_3d_test \
	$(SRCDIR)/modules/modules.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fastfft.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/GPSTD.o \
	$(SRCDIR)/field_solvers/Maxwell/yee_solver/yee.o \
	$(SRCDIR)/field_solvers/Maxwell/karkkainen_solver/karkkainen.o \
	$(SRCDIR)/parallelization/tiling/tiling.o \
	$(SRCDIR)/housekeeping/sorting.o \
	$(SRCDIR)/particle_pushers/vay_pusher/vay_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_3d.o \
	$(SRCDIR)/particle_pushers/boris_pusher/boris_2d.o \
	$(SRCDIR)/particle_pushers/kin_energy.o \
	$(SRCDIR)/particle_pushers/laser_pusher_manager_3d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_2d.o \
	$(SRCDIR)/particle_pushers/particle_pusher_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/current_deposition_manager_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/direct/direct_current_deposition_3d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_2d.o \
	$(SRCDIR)/particle_deposition/current_deposition/esirkepov/esirkepov_3d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_2d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_2d.o \
	$(SRCDIR)/field_gathering/field_gathering_manager_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_on_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o1_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o2_3d.o \
	$(SRCDIR)/field_gathering/energy_conserving/field_gathering_o3_3d.o \
	$(SRCDIR)/parallelization/mpi/mpi_derived_types.o \
	$(SRCDIR)/housekeeping/load_balancing.o \
	$(SRCDIR)/field_solvers/Maxwell/maxwell_solver_manager.o \
	$(SRCDIR)/boundary_conditions/field_boundaries.o \
	$(SRCDIR)/boundary_conditions/particle_boundaries.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.o \
	$(SRCDIR)/field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_manager.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_2d.o \
	$(SRCDIR)/particle_deposition/charge_deposition/charge_deposition_3d.o \
	$(SRCDIR)/diags/diags.o \
	$(SRCDIR)/ios/simple_io.o \
	$(SRCDIR)/parallelization/mpi/mpi_routines.o \
	$(SRCDIR)/submain.o \
	$(SRCDIR)/initialization/control_file.o \
	Acceptance_testing/Gcov_tests/maxwell_3d_test.o $(LDFLAGS)

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
	build_esirkepov_2d_test \
	build_tile_mpi_part_com_test

build_test_spectral_3d: createdir \
	build_maxwell_3d_test
build_test_spectral_2d: createdir \
	build_maxwell_2d_test

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
	esirkepov_2d_test \
	tile_mpi_part_com_test

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

# __ Execute Maxwell solver tests ____________________________________________________

test_maxwell_solver_3d: test_plane_wave_fdtd_3d \
		test_plane_wave_psatd_3d

test_maxwell_solver_2d:	test_plane_wave_fdtd_2d \
		test_plane_wave_psatd_2d


test_plane_wave_fdtd_2d:
	cp examples/example_decks_fortran/plane_wave_test_2d.pixr \
	Acceptance_testing/Gcov_tests/input_file.pixr
	# 1 OpenMP vary number of MPIs
	export OMP_NUM_THREADS=1
	cd Acceptance_testing/Gcov_tests && \
	mpirun -np 1 ./maxwell_2d_test input_file.pixr --nsteps 161 --l_spectral .FALSE. &&\
	mpirun -np 2 ./maxwell_2d_test input_file.pixr --nsteps 161 --l_spectral .FALSE. --nprocz 2 &&\
	mpirun -np 4 ./maxwell_2d_test input_file.pixr --nsteps 161 --l_spectral .FALSE. --nprocz 2 --nprocx 2 
test_plane_wave_psatd_2d: 
	cp examples/example_decks_fortran/plane_wave_test_2d.pixr \
	Acceptance_testing/Gcov_tests/input_file.pixr
	# 1 OpenMP vary number of MPIs
	export OMP_NUM_THREADS=1
	cd Acceptance_testing/Gcov_tests && \
	mpirun -np 1 ./maxwell_2d_test input_file.pixr --nsteps 81 --l_spectral .TRUE. &&\
	mpirun -np 2 ./maxwell_2d_test input_file.pixr --nsteps 81 --l_spectral .TRUE. --nprocz 2 &&\
	mpirun -np 4 ./maxwell_2d_test input_file.pixr --nsteps 81 --l_spectral .TRUE. --nprocz 2 --nprocx 2

test_plane_wave_fdtd_3d:
	cp examples/example_decks_fortran/plane_wave_test.pixr \
	Acceptance_testing/Gcov_tests/input_file.pixr
	# 1 OpenMP vary number of MPIs
	cd Acceptance_testing/Gcov_tests && \
	export OMP_NUM_THREADS=1 && \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .FALSE. --nsteps 106 && \
	mpirun -np 4 ./maxwell_3d_test input_file.pixr --l_spectral .FALSE. --nsteps 106 && \
	mpirun -np 8 ./maxwell_3d_test input_file.pixr --l_spectral .FALSE. --nsteps 106
	# 4 OpenMP vary number of MPIs 
	cd Acceptance_testing/Gcov_tests && \
	export OMP_NUM_THREADS=4 && \
	mpirun -np 1 ./maxwell_3d_test input_file.pixr --l_spectral .FALSE. --nsteps 106 && \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .FALSE. --nsteps 106 && \
	mpirun -np 4 ./maxwell_3d_test input_file.pixr --l_spectral .FALSE. --nsteps 106
test_plane_wave_psatd_3d:
	cd Acceptance_testing/Gcov_tests && \
	export OMP_NUM_THREADS=1 && \
	mpirun -np 1 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 && \
	mpirun -np 1 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 --norderz 8 --nordery 8 --norderx 8  && \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 && \
	mpirun -np 4 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 && \
	mpirun -np 8 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61
test_plane_wave_psatd_global_3d:
	cd Acceptance_testing/Gcov_tests && \
	export OMP_NUM_THREADS=4 && \
	mpirun -np 1 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --fftw_with_mpi .TRUE. --nsteps 61  && \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 --fftw_with_mpi .TRUE. && \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 --fftw_with_mpi .TRUE. && \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 --fftw_with_mpi .TRUE. --fftw_hybrid .TRUE. --nb_group 1 \
	mpirun -np 2 ./maxwell_3d_test input_file.pixr --l_spectral .TRUE. --nsteps 61 --fftw_with_mpi .TRUE. --fftw_hybrid .TRUE. --nb_group 1 \
	--fftw_mpi_tr .TRUE. 

test_lb:
	cd Acceptance_testing/Gcov_tests && \
	export OMP_NUM_THREADS=1 &&\
	mpirun -np 8 ./maxwell_2d_test input_file.pixr --l_spectral .TRUE. --nsteps 81 --nprocz 8 --fftw_with_mpi .TRUE. --fftw_hybrid .TRUE. --nb_group 1 \
	is_lb_grp .TRUE.	

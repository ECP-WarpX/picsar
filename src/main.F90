!===============================================================================
! PICSAR 3D version 2.0 with tiling, H. VINCENTI 09/15/2015
! INCLUDES: 
! - Arbitrary order field solver (Maxwell.F90)
! - High order current deposition/field gathering routines (current_deposition.F90, field_gathering.F90)
! - MPI-parallelization (mpi_subtype_control.F90, mpi_routines.F90, boundary.F90)
! - MPI-IO outputs
! - OpenMP Parallelization (current_deposition.F90, field_gathering.F90, particle_push.F90, Maxwell.F90)
! - Tiling of particles for better memory locality
!===============================================================================

PROGRAM main
!===============================================================================
USE constants
USE fields
USE particles
USE params
USE shared_data
USE mpi_routines
USE control_file

IMPLICIT NONE
INTEGER :: i,ierror,j,l

! --- default init
  CALL default_init

! --- reads input_file
  CALL read_input_file

! --- reads from command line
  CALL read_from_cl

! --- mpi init communicator
  CALL mpi_minimal_init

! --- Check domain decomposition / Create Cartesian communicator / Allocate grid arrays
  CALL mpi_initialise

! --- allocates and inits particle distributions (on each subdomain)
  CALL initall

!----------------------------------------------
! THIS IS THE PIC ALGORITHM TIME LOOP
!----------------------------------------------
IF (rank .EQ. 0) startsim=MPI_WTIME()
CALL step(nsteps)
IF (rank .EQ. 0) endsim=MPI_WTIME()
IF (rank .EQ. 0) WRITE(0,*)  "Total runtime on ",nproc," CPUS =", endsim-startsim
CALL mpi_close

END PROGRAM main

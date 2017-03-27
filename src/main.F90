! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016,
! The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a paid-up,
! nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
! copies to the public, prepare derivative works, and perform publicly and display
! publicly, and to permit other to do so.
!
! PICSAR
!
! version 2.0
! Creation date: 09/15/2015
!
! Developers:
!  Henri Vincenti
!  Mathieu Lobet
!  Remi Lehe
!  Jean-Luc Vay
!  Guillaume Blaclard
!
! INCLUDES:
! - Arbitrary order field solver (Maxwell.F90)
! - High order current deposition/field gathering routines (current_deposition.F90, field_gathering.F90)
! - MPI-domain decomposition (mpi_subtype_control.F90, mpi_routines.F90, boundary.F90)
! - Tiling of particles for better memory locality (tiling.F90)
! - OpenMP Hybrid Parallelization (current_deposition.F90, field_gathering.F90, particle_push.F90, Maxwell.F90)
! - MPI-IO outputs

! ________________________________________________________________________________________


PROGRAM main

  USE constants
  USE fields
  USE particles
  USE params
  USE shared_data
  USE mpi_routines
  USE control_file
  USE time_stat
  USE diagnostics

! Vtune profiling
#if (defined(VTUNE) && VTUNE>0)
  USE ITT_FORTRAN
#endif
! SDE profiling
#if (defined(SDE) && SDE>0)||(defined(DFP))
  USE SDE_FORTRAN
#endif

  IMPLICIT NONE

! Intel Design Forward project
#if defined(DFP)
  CALL DFP_INIT_START
#endif

! --- default init
  CALL default_init

! --- reads input_file
  CALL read_input_file

#if (defined(VTUNE) || defined(SDE) || defined(DFP) || defined(ALLINEA))
! No command line
#else
! --- reads from command line
  CALL read_from_cl
#endif

! --- mpi init communicator
  IF (fftw_with_mpi) THEN 
    CALL mpi_minimal_init_fftw
  ELSE
    CALL mpi_minimal_init
  ENDIF 
  IF (rank .EQ. 0) THEN
    write(0,*) "_________________________________________________________________"
    write(0,*) ""
    write(0,*) " PICSAR"
    write(0,*) "_________________________________________________________________"
  ENDIF

! --- Check domain decomposition / Create Cartesian communicator / Allocate grid arrays
  CALL mpi_initialise

! --- allocates and inits particle distributions (on each subdomain)
  CALL initall

! --- Diagnostics
  CALL init_diags

! Intel Design Forward project
#if defined(DFP)
   CALL DFP_INIT_STOP
#endif

  !----------------------------------------------
  ! THIS IS THE PIC ALGORITHM TIME LOOP
  !----------------------------------------------
  IF (rank .EQ. 0) startsim=MPI_WTIME()
  CALL step(nsteps)
  IF (rank .EQ. 0) endsim=MPI_WTIME()
  IF (rank .EQ. 0) WRITE(0,*)  "Total runtime on ",nproc," CPUS =", endsim-startsim

  CALL time_statistics

  CALL mpi_close

! Intel Design Forward project
#if defined(DFP)
   CALL DFP_FINAL_STOP
#endif

END PROGRAM main

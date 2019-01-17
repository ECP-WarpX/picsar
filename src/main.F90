! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c) 2016,
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
! - High order current deposition/field gathering routines (current_deposition.F90, 
! - field_gathering.F90)
! - MPI-domain decomposition (mpi_subtype_control.F90, mpi_routines.F90, boundary.F90)
! - Tiling of particles for better memory locality (tiling.F90)
! - OpenMP Hybrid Parallelization (current_deposition.F90, field_gathering.F90, 
! - particle_push.F90, Maxwell.F90)
! - MPI-IO outputs
! ________________________________________________________________________________________


PROGRAM main
  USE PICSAR_precision
  USE constants
  USE fields
  USE particles
  USE params
  USE shared_data
  USE mpi_routines
  USE control_file
  USE time_stat
  USE diagnostics
  USE mem_status, ONLY : global_grid_mem, global_grid_tiles_mem, global_part_tiles_mem
#if defined(FFTW)
  USE mpi_fftw3
  USE fourier
  USE fastfft
  USE fftw3_fortran
#endif

! Vtune profiling
#if (defined(VTUNE) && VTUNE>0)
  USE ITT_FORTRAN
#endif
! SDE profiling
#if (defined(SDE) && SDE>0)||(defined(DFP))
  USE SDE_FORTRAN
#endif

  IMPLICIT NONE
  LOGICAL(lp) :: exist
  CHARACTER(len=250) :: str1, str2, str3
  CHARACTER(len=250) :: str4, str5, str7
  CHARACTER(len=500) :: str6
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
#if defined(FFTW)
  IF (fftw_with_mpi) THEN 
    CALL mpi_minimal_init_fftw
  ELSE
#endif
    CALL mpi_minimal_init
#if defined(FFTW)
  ENDIF 
#endif
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
  IF (rank .EQ. 0) WRITE(0,*)  "Total runtime on ",nproc," CPUS =",                   &
  endsim-startsim,"CPU AVERG TIME PER IT",(endsim-startsim)/nsteps


  ! Time statistics for the different processes of the PIC step
  CALL time_statistics

  IF (rank .EQ. 0) THEN 
	  INQUIRE(file="output_statistics.out", exist=exist)
  	IF (exist) THEN 
  		OPEN (unit=12,file="output_statistics.out", &
  		action="write",position="append", status="old")
  	ELSE
  		OPEN (unit=12,file="output_statistics.out",  &
  		  action="write",status="new")
  	ENDIF 
  	WRITE(str1,*) nx_global; WRITE(str2,*) ny_global
  	WRITE(str3,*) nz_global; WRITE(str4,*) nproc
  	! total simulation time
  	WRITE(str5,*) endsim-startsim
  	! Average time spent in different steps of the PIC loop
  	WRITE(str6,'(22(E12.5))') avetimes(1),avetimes(14),avetimes(2),avetimes(11),      &
  								avetimes(3),avetimes(4),avetimes(5),                  &
  								avetimes(6),avetimes(7),avetimes(21),                 &
  								avetimes(22),avetimes(23),avetimes(24), avetimes(25), &
  								avetimes(8),avetimes(10),                             &
  								avetimes(12), avetimes(13), avetimes(9),              &
  								avetimes(18), avetimes(19), avetimes(20)
  								
  	! Total memory used in the case (in GB)
  	WRITE(str7,'(4(E12.5))') global_grid_mem/1e9, global_grid_tiles_mem/1e9,          &
  	global_part_tiles_mem/1e9

  	! All time are put in the file on a single line
  	WRITE(12, '(512A)')  trim(adjustl(str1))//" "//trim(adjustl(str2))//" "//         &
  				  trim(adjustl(str3))//" "//trim(adjustl(str4))//" "//                &
  				  trim(adjustl(str5))//" "//trim(adjustl(str6))//                     &
  				  " "//trim(adjustl(str7))
  	CLOSE(12)
  ENDIF 

#if defined(FFTW)
  IF(l_spectral) THEN
    IF(fftw_with_mpi) THEN
      CALL DFFTW_DESTROY_PLAN(plan_r2c_mpi)
      CALL DFFTW_DESTROY_PLAN(plan_c2r_mpi)
    ELSE
      CALL fast_fftw_destroy_plan_dft(plan_r2c)
      CALL fast_fftw_destroy_plan_dft(plan_c2r)
    ENDIF
  ENDIF
#endif
  CALL mpi_close

! Intel Design Forward project
#if defined(DFP)
   CALL DFP_FINAL_STOP
#endif

END PROGRAM main

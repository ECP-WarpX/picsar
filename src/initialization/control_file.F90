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
! CONTROL_FILE.F90
!
! This file contains subroutines for reading input files and interpreting
! command line arguments.
!
! Author
! Henri Vincenti
! Mathieu Lobet
!
! date
! Creation 2015
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Module containing routines for the default initialization, for reading
!> the input file and the command line arguments.
!
!> @author
!> Henri Vincenti, !> Mathieu Lobet
!
!> @date
!> 2015-2016
!
! ________________________________________________________________________________________
MODULE control_file
  USE shared_data
  USE params
  USE fields
  USE particles
  USE params
  USE output_data
  USE time_stat
#if defined(FFTW)
  USE group_parameters
#endif
  IMPLICIT NONE
  INTEGER(idp) :: ios=0
  INTEGER(idp), PARAMETER :: fh_input = 15
  CHARACTER(LEN=string_length) :: buffer
  CHARACTER(LEN=string_length) :: section_name

  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that proceeds to default init of parameters of the PIC loop.
  !> This init is performed first in the main program.  It is called by ::main
  !>
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE default_init
    ! --- Dimension
    c_dim = 3

    ! --- Init particle tiling split
    ntilex = 1
    ntiley = 1
    ntilez = 1

    ! --- Order of Maxwell field solver (default is 2 in x, y, z)
    norderx = 2
    nordery = 2
    norderz = 2
    nx_pml = 0_idp 
    ny_pml = 0_idp 
    nz_pml = 0_idp
    absorbing_bcs = .FALSE.
    absorbing_bcs_x = .FALSE.
    absorbing_bcs_y = .FALSE.
    absorbing_bcs_z = .FALSE.
    l_nodalgrid = .FALSE.
    l_spectral = .FALSE.! (no spectral solver by default)
    g_spectral = .FALSE.! (no spectral sovler by default)
    l_staggered = .TRUE.! (staggered scheme by default )
#if defined(FFTW)
    nb_group_x = 1
    nb_group_y = 1
    nb_group_z = 1
#endif
    ! --- Order of current deposition/ field gathering
    ! (default is 1 in x, y, z)
    nox = 1
    noy = 1
    noz = 1
    nxguards=MAX(nox, 2_idp)
    nyguards=MAX(noy, 2_idp)
    nzguards=MAX(noz, 2_idp)
    nxjguards=MAX(nox, 2_idp)
    nyjguards=MAX(noy, 2_idp)
    nzjguards=MAX(noz, 2_idp)

#if defined(FFTW)
    nxg_group=max(nox, 2_idp)
    nyg_group=max(noy, 2_idp)
    nzg_group=max(noz, 2_idp)
#endif
    ! Topology
    topology = 0

    ! MPI communication
    mpicom_curr = 1

    ! Current deposition algorithm
    currdepo = 0

    ! Charge deposition algorithm
    rhodepo = 0

    ! Field gathering algorithm
    fieldgathe = 0

    ! Particle communication routine
    partcom = 0

    ! Particle pusher algorithm
    particle_pusher = 0

    ! Field gathering + part. pusher
    fg_p_pp_separated = 0

    ! Vector length current deposition
    lvec_curr_depo = 8

    ! Vector length charge deposition
    LVEC_charge_depo = 64

    ! Vector length field gathering
    LVEC_fieldgathe = 256

    ! Size of the particle mpi buffer
    mpi_buf_size = 2000

    ! Sorting activation (not activated by default)
    sorting_activated = 0
    sorting_dx = 1.
    sorting_dy = 1.
    sorting_dz = 1.
    sorting_shiftx = 0.
    sorting_shifty = 0.
    sorting_shiftz = 0.
    sorting_verbose = .TRUE.

    ! Time stats output activation
    timestat_activated = 0
    timestat_period = 0
    timestat_itstart = 0
    timestat_perit = 0
    nbuffertimestat = 1

    l_lower_order_in_v = .TRUE.

    ! --- sets coefficient multiplying Courant time step
    dtcoef = 0.7_num
    ! --- smoothing
    npass = 0
    alpha = 0.5_num
    ! --- sets max time in the simulation (in 1/w0)
    tmax = 0._num
    nsteps = 0

    !-------------------------------------------------------------------------------
    ! plasma parameters (cold plasma)
    l_particles_weight = .FALSE.! .TRUE. if particles have different weights

    ! --- quantities in plasma (or lab) frame
    !-------------------------------------------------------------------------------
    nlab  = 1.e23_num! plasma density in lab frame
    g0    = 1.0_num! initial gamma
    b0    = sqrt(1.0_num-1.0_num/g0**2)
    nc    = nlab*g0! density (in the simulation frame)
    ! plasma frequency (in the lab frame)
    wlab  = echarge*sqrt(nlab/(emass*eps0))
    ! "longitudinal" plasma frequency (in the lab frame)
    w0_l  = echarge*sqrt(nc/(g0*emass*eps0))
    ! "transverse" plasma frequency (in the lab frame)
    w0_t  = echarge*sqrt(nc/(g0**3*emass*eps0))
    w0    = w0_l
    ! --- Init number of species
    nspecies=0
    ! --- Init number of particle dumps
    npdumps = 0
    ! --- l_plasma
    l_plasma= .TRUE.
    ! --- Particle distribution
    pdistr=1
    ! Init species array
    IF (.NOT. l_species_allocated) THEN
      nspecies=0
      ALLOCATE(species_parray(1:nspecies_max))
      l_species_allocated=.TRUE.
    ENDIF

    ! Particle boundaries (0 - periodic by default)
    pbound_x_min=0
    pbound_y_min=0
    pbound_z_min=0
    pbound_x_max=0
    pbound_y_max=0
    pbound_z_max=0

    ! Temporal output
    temdiag_frequency = 0
    temdiag_format = 0

    ! SET FFTW WITH MPI FLAG
    fftw_with_mpi = .FALSE.
    fftw_hybrid = .FALSE.
    fftw_mpi_transpose = .FALSE.
    p3dfft_flag = .FALSE.
#if defined(P3DFFT) 
    p3dfft_stride = .TRUE.
#endif
  END SUBROUTINE default_init

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads command line arguments.
  !
  !> @details
  !> This routine allow to pass command line arguments to PXR
  !> and is very useful when performing parametric studies.
  !> This init is performed after every other init and thus overwrites parameters
  !> defined in default_init or read_input_file. It is called by ::main
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_from_cl
    INTEGER :: i
    DO i = 1, COMMAND_ARGUMENT_COUNT()-1, 2
      CALL GETARG(i, buffer)
      IF (INDEX(buffer, 'ntilex') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') ntilex
      ELSE IF (INDEX(buffer, 'ntiley') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') ntiley
      ELSE IF (INDEX(buffer, 'ntilez') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') ntilez
      ELSE IF (INDEX(buffer, 'distr') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') pdistr
      ELSE IF (INDEX(buffer, 'nprocx') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') nprocx
      ELSE IF (INDEX(buffer, 'nprocy') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') nprocy
      ELSE IF (INDEX(buffer, 'nprocz') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') nprocz
      ELSE IF (INDEX(buffer, 'nox') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') nox
      ELSE IF (INDEX(buffer, 'noy') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') noy
      ELSE IF (INDEX(buffer, 'noz') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') noz
      ELSE IF (INDEX(buffer, 'tmax') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) tmax
      ELSE IF (INDEX(buffer, 'nsteps') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nsteps
      ELSE IF (INDEX(buffer, 'dtcoef') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) dtcoef
      ELSE IF (INDEX(buffer, 'nx') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nx_global_grid
      ELSE IF (INDEX(buffer, 'ny') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) ny_global_grid
      ELSE IF (INDEX(buffer, 'nz') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nz_global_grid
      ELSE IF (INDEX(buffer, 'fieldgathe') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) fieldgathe
      ELSE IF (INDEX(buffer, 'currdepo') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) currdepo
      ELSE IF (INDEX(buffer, 'rhodepo') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) rhodepo
      ELSE IF (INDEX(buffer, 'partcom') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) partcom
      ELSE IF (INDEX(buffer, 'particle_pusher') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) particle_pusher
      ELSE IF (INDEX(buffer, 'sorting') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) sorting_activated
      ELSE IF (INDEX(buffer, 'lvec_curr_depo') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) lvec_curr_depo
      ELSE IF (INDEX(buffer, 'l_plasma') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) l_plasma
      ELSE IF (INDEX(buffer, 'l_spectral') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) l_spectral
      ELSE IF (INDEX(buffer, 'absorbing_bcs_x') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) absorbing_bcs_x
      ELSE IF (INDEX(buffer, 'absorbing_bcs_y') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) absorbing_bcs_y
      ELSE IF (INDEX(buffer, 'absorbing_bcs_z') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) absorbing_bcs_z
      ELSE IF (INDEX(buffer, 'nxpml') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') nx_pml
      ELSE IF (INDEX(buffer, 'nypml') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') ny_pml
      ELSE IF (INDEX(buffer, 'nzpml') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, '(i10)') nz_pml
      ELSE IF (INDEX(buffer, 'g_spectral') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) g_spectral
      ELSE IF (INDEX(buffer, 'fftw_with_mpi') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) fftw_with_mpi
      ELSE IF (INDEX(buffer, 'fftw_hybrid') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) fftw_hybrid
      ELSE IF (INDEX(buffer, 'p3dfft_flag') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) p3dfft_flag
      ELSE IF (INDEX(buffer, 'p3dfft_stride') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) p3dfft_stride
      ELSE IF (INDEX(buffer, 'fftw_mpi_tr') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) fftw_mpi_transpose
      ELSE IF (INDEX(buffer, 'lvec_charge_depo') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) lvec_charge_depo
      ELSE IF (INDEX(buffer, 'lvec_fieldgathe') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) lvec_fieldgathe
      ELSE IF (INDEX(buffer, 'mpicom_curr') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) mpicom_curr
      ELSE IF (INDEX(buffer, 'mpi_buf_size') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) mpi_buf_size
      ELSE IF (INDEX(buffer, 'nsteps') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nsteps
      ELSE IF (INDEX(buffer, 'nguardsx') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nxguards
        nxjguards=nxguards
      ELSE IF (INDEX(buffer, 'nguardsy') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nyguards
        nyjguards=nyguards
      ELSE IF (INDEX(buffer, 'nguardsz') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nzguards
        nzjguards=nzguards
#if defined(FFTW)
      ELSE IF (INDEX(buffer, 'ngguards_x') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nxg_group
      ELSE IF (INDEX(buffer, 'ngguards_y') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nyg_group
      ELSE IF (INDEX(buffer, 'ngguards_z') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nzg_group
#endif
#if defined(FFTW)
      ELSE IF (INDEX(buffer, 'nb_group_z') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nb_group_z
#if defined(P3DFFT) 
      ELSE IF (INDEX(buffer, 'nb_group_y') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) nb_group_y
#endif
#endif
      ELSE IF (INDEX(buffer, 'c_dim') .GT. 0) THEN
        CALL GETARG(i+1, buffer)
        READ(buffer, *) c_dim 
        IF(c_dim == 2) ny = 1_idp
      END IF
    END DO
 
    ! > IF absorbing_bcs_i == .TRUE. then set absorbing_bcs = .TRUE. 
    ! > to init splitted fields block matri
    IF(absorbing_bcs_x .OR. absorbing_bcs_y .OR. absorbing_bcs_z) absorbing_bcs = .TRUE.

    RETURN
  END SUBROUTINE read_from_cl

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads simulation parameters from an input file (of name input_file.pixr)
  !> This init is performed after the call to default_init and before the call to
  !> to read_from_cl. It is called by ::main
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_input_file
    INTEGER :: ix = 0
    CHARACTER(len=32) :: input_filename, input_file_default='input_file.pixr'
    
    CALL get_command_argument(1, input_filename)
    IF (LEN_TRIM(input_filename) == 0) input_filename = input_file_default
    
    ! --- OPENS INPUT FILE
    OPEN(fh_input, file=TRIM(input_filename), status='OLD', iostat=ios)
    
    IF (ios/=0) then
        write(0,*) '#############################'
        write(0,*) '# Error opening input file. # '
        write(0,*) '# Check input filename.     # '
        write(0,*) '#############################'
        call abort()
    END IF
       
    DO WHILE(ios==0)
      READ(fh_input, '(A)', iostat=ios) buffer
      ix=INDEX(buffer, 'section::')
      IF (ix .GT. 0) THEN
        section_name=buffer(ix:string_length)
        !write(0, *) TRIM(ADJUSTL(section_name))
        SELECT CASE(TRIM(ADJUSTL(section_name)))
        CASE('section::main')
          CALL read_main_section
        CASE('section::species')
          CALL read_species_section
        CASE('section::antenna')
          CALL read_antenna_section
        CASE('section::output')
          CALL read_output_section
        CASE('section::cpusplit')
          CALL read_cpusplit_section
        CASE('section::plasma')
          CALL read_plasma_section
        CASE('section::temporal')
          CALL read_temporal_output_section
        CASE('section::solver')
          CALL read_solver_section
        CASE('section::timestat')
          CALL read_timestat_section
        CASE('section::sorting')
          CALL read_sorting_section
        CASE('section::particle_dump')
          CALL read_particle_dumps_section
        END SELECT
      END IF
    END DO
    RETURN
  END SUBROUTINE read_input_file

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the cpu section in the input file
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_cpusplit_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS CPUSPLIT SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'nprocx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nprocx
      ELSE IF (INDEX(buffer, 'nprocy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nprocy
      ELSE IF (INDEX(buffer, 'nprocz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nprocz
      ELSE IF (INDEX(buffer, 'topology') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') topology
      ELSE IF (INDEX(buffer, 'mpicom_curr') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') mpicom_curr
      ELSE IF (INDEX(buffer, 'partcom') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') partcom
      ELSE IF (INDEX(buffer, 'fg_p_pp_separated') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') fg_p_pp_separated
      ELSE IF (INDEX(buffer, 'particle_pusher') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') particle_pusher
      ELSE IF (INDEX(buffer, 'lvec_curr_depo') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') lvec_curr_depo
      ELSE IF (INDEX(buffer, 'lvec_charge_depo') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') lvec_charge_depo
      ELSE IF (INDEX(buffer, 'lvec_fieldgathe') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') lvec_fieldgathe
      ELSE IF (INDEX(buffer, 'mpi_buf_size') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') mpi_buf_size
      ELSE IF (INDEX(buffer, 'c_dim') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_dim
        IF(c_dim == 2) ny = 1_idp
      ELSE IF (INDEX(buffer, 'end::cpusplit') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_cpusplit_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the plasma main properties section in the input file
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_plasma_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS CPUSPLIT SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'nlab') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) nlab
      ELSE IF (INDEX(buffer, 'gamma0') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) g0
      ELSE IF (INDEX(buffer, 'pdistr') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) pdistr
      ELSE IF (INDEX(buffer, 'particle_pusher') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') particle_pusher
      ELSE IF (INDEX(buffer, 'end::plasma') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_plasma_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the solver parameters section in the input file
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_solver_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS CPUSPLIT SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'norderx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') norderx
      ELSE IF (INDEX(buffer, 'nordery') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nordery
      ELSE IF (INDEX(buffer, 'norderz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') norderz
      ELSE IF (INDEX(buffer, 'nox') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nox
      ELSE IF (INDEX(buffer, 'nx_pml') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nx_pml
      ELSE IF (INDEX(buffer, 'ny_pml') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') ny_pml
      ELSE IF (INDEX(buffer, 'nz_pml') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nz_pml
      ELSE IF (INDEX(buffer, 'noy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') noy
      ELSE IF (INDEX(buffer, 'noz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') noz
      ELSE IF (INDEX(buffer, 'currdepo') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') currdepo
      ELSE IF (INDEX(buffer, 'fieldgathe') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') fieldgathe
      ELSE IF (INDEX(buffer, 'rhodepo') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') rhodepo
      ELSE IF (INDEX(buffer, 'partcom') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') partcom
      ELSE IF (INDEX(buffer, 'mpi_buf_size') .GT. 0) THEN
        ix = INDEX(buffer, "=")
      ELSE IF (INDEX(buffer, 'l_spectral') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) l_spectral
      ELSE IF (INDEX(buffer, 'absorbing_bcs_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) absorbing_bcs_x
      ELSE IF (INDEX(buffer, 'absorbing_bcs_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) absorbing_bcs_y
      ELSE IF (INDEX(buffer, 'absorbing_bcs_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) absorbing_bcs_z
      ELSE IF (INDEX(buffer, 'g_spectral') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) g_spectral
      ELSE IF (INDEX(buffer, 'l_staggered') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) l_staggered
      ELSE IF (INDEX(buffer, 'fftw_with_mpi') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) fftw_with_mpi
      ELSE IF (INDEX(buffer, 'fftw_mpi_tr') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) fftw_mpi_transpose
      ELSE IF (INDEX(buffer, 'fftw_hybrid') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) fftw_hybrid
      ELSE IF (INDEX(buffer, 'p3dfft_flag') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) p3dfft_flag
      ELSE IF (INDEX(buffer, 'p3dfft_stride') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) p3dfft_stride
#if defined(FFTW)
      ELSE IF (INDEX(buffer, 'nb_group_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nb_group_z
#if defined(P3DFFT)
      ELSE IF (INDEX(buffer, 'nb_group_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nb_group_y
#endif
#endif
      ELSE IF (INDEX(buffer, 'fg_p_pp_separated') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') fg_p_pp_separated
      ELSE IF (INDEX(buffer, 'end::solver') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO

    ! > IF absorbing_bcs_i == .TRUE. then set absorbing_bcs = .TRUE. 
    ! > to init splitted fields block matrixes
    IF(absorbing_bcs_x .OR. absorbing_bcs_y .OR. absorbing_bcs_z) absorbing_bcs = .TRUE.

    RETURN
  END SUBROUTINE read_solver_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the particle sorting parameters section in the input file
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_sorting_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS CPUSPLIT SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'activation') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') sorting_activated
      ELSE IF (INDEX(buffer, 'dx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_dx
      ELSE IF (INDEX(buffer, 'dy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_dy
      ELSE IF (INDEX(buffer, 'dz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_dz
      ELSE IF (INDEX(buffer, 'shiftx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_shiftx
      ELSE IF (INDEX(buffer, 'shifty') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_shifty
      ELSE IF (INDEX(buffer, 'shiftz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_shiftz
      ELSE IF (INDEX(buffer, 'verbose') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) sorting_verbose
      ELSE IF (INDEX(buffer, 'end::sorting') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_sorting_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the time statistics parameters section in the input file
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_timestat_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS CPUSPLIT SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'activation') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') timestat_activated
      ELSE IF (INDEX(buffer, 'period') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') timestat_period
      ELSE IF (INDEX(buffer, 'it_start') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') timestat_itstart
      ELSE IF (INDEX(buffer, 'per_it') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') timestat_perit
      ELSE IF (INDEX(buffer, 'buffersize') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nbuffertimestat
      ELSE IF (INDEX(buffer, 'end::timestat') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_timestat_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the general parameters section in the input file
  !> including the domain extension, the discretization, the tiling, the guard cells
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_main_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS GRID SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'c_dim') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_dim
      ELSE IF (INDEX(buffer, 'nx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nx_global_grid
        nx_global=nx_global_grid-1
      ELSE IF (INDEX(buffer, 'ny') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') ny_global_grid
        ny_global=ny_global_grid-1
      ELSE IF (INDEX(buffer, 'nz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nz_global_grid
        nz_global=nz_global_grid-1
      ELSE IF (INDEX(buffer, 'ntilex') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') ntilex
      ELSE IF (INDEX(buffer, 'ntiley') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') ntiley
      ELSE IF (INDEX(buffer, 'ntilez') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') ntilez
      ELSEIF (INDEX(buffer, 'dx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dx
      ELSE IF (INDEX(buffer, 'dy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dy
      ELSE IF (INDEX(buffer, 'dz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dz
      ELSEIF (INDEX(buffer, 'xmin') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) xmin
      ELSE IF (INDEX(buffer, 'ymin') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) ymin
      ELSE IF (INDEX(buffer, 'zmin') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) zmin
      ELSEIF (INDEX(buffer, 'xmax') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) xmax
      ELSE IF (INDEX(buffer, 'ymax') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) ymax
      ELSE IF (INDEX(buffer, 'zmax') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) zmax
      ELSE IF (INDEX(buffer, 't_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) tmax
      ELSE IF (INDEX(buffer, 'nsteps') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) nsteps
      ELSE IF (INDEX(buffer, 'dtcoef') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dtcoef
      ELSE IF (INDEX(buffer, 'nguardsx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nxguards
#if defined(FFTW)
        nxg_group=nxguards
#endif
      ELSE IF (INDEX(buffer, 'nguardsy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nyguards
#if defined(FFTW)
        nyg_group=nyguards
#endif
      ELSE IF (INDEX(buffer, 'nguardsz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nzguards
#if defined(FFTW)
        nzg_group=nzguards
#endif
      ELSE IF (INDEX(buffer, 'njguardsx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nxjguards
      ELSE IF (INDEX(buffer, 'njguardsy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nyjguards
      ELSE IF (INDEX(buffer, 'njguardsz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nzjguards
#if defined(FFTW)
      ELSE IF (INDEX(buffer, 'ngguards_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nzg_group
      ELSE IF (INDEX(buffer, 'ngguards_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nyg_group
      ELSE IF (INDEX(buffer, 'ngguards_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') nxg_group
#endif
      ELSE IF (INDEX(buffer, 'l_plasma') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) l_plasma
      ELSE IF (INDEX(buffer, 'end::main') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_main_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the species properties section in the input file
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_species_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section
    TYPE(particle_species), POINTER :: curr
    ! READS SPECIES SECTION OF INPUT FILE
    IF (.NOT. l_species_allocated) THEN
      nspecies=0
      ALLOCATE(species_parray(1:nspecies_max))
      l_species_allocated=.TRUE.
    ENDIF
    nspecies = nspecies+1
    curr => species_parray(nspecies)
    ! minimal init for species attributes
    curr%charge = -echarge
    curr%mass = emass
    curr%nppcell = 0
    curr%x_min = 0._num
    curr%x_max = 0._num
    curr%y_min = 0._num
    curr%y_max = 0._num
    curr%z_min = 0._num
    curr%z_max = 0._num
    curr%vdrift_x =0._num
    curr%vdrift_y =0._num
    curr%vdrift_z =0._num
    curr%vth_x =0._num
    curr%vth_y =0._num
    curr%vth_z =0._num
    curr%sorting_period = 0
    curr%sorting_start = 0
    curr%species_npart=0
    end_section=.FALSE.
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'name') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%name
      ELSE IF (INDEX(buffer, 'mass') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%mass
        curr%mass=curr%mass*emass
      ELSE IF (INDEX(buffer, 'charge') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%charge
        curr%charge=curr%charge*echarge
      ELSEIF (INDEX(buffer, 'nppcell') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') curr%nppcell
      ELSE IF (INDEX(buffer, 'x_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%x_min
      ELSE IF (INDEX(buffer, 'x_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%x_max
      ELSEIF (INDEX(buffer, 'y_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%y_min
      ELSE IF (INDEX(buffer, 'y_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%y_max
      ELSEIF (INDEX(buffer, 'z_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%z_min
      ELSE IF (INDEX(buffer, 'z_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%z_max
      ELSE IF (INDEX(buffer, 'vdrift_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%vdrift_x
        curr%vdrift_x=curr%vdrift_x*clight
      ELSE IF (INDEX(buffer, 'vdrift_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%vdrift_y
        curr%vdrift_y=curr%vdrift_y*clight
      ELSE IF (INDEX(buffer, 'vdrift_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%vdrift_z
        curr%vdrift_z=curr%vdrift_z*clight
      ELSE IF (INDEX(buffer, 'vth_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%vth_x
        curr%vth_x=curr%vth_x*clight
      ELSE IF (INDEX(buffer, 'vth_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%vth_y
        curr%vth_y=curr%vth_y*clight
      ELSE IF (INDEX(buffer, 'vth_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%vth_z
        curr%vth_z=curr%vth_z*clight
      ELSE IF (INDEX(buffer, 'sorting_period') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%sorting_period
      ELSE IF (INDEX(buffer, 'sorting_start') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%sorting_start
      ELSE IF (INDEX(buffer, 'end::species') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_species_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads the particle dump parameters section
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_particle_dumps_section
    INTEGER                      :: ix = 0, ispecies
    LOGICAL(lp)                       :: end_section
    TYPE(particle_dump), POINTER :: dp
    CHARACTER(LEN=string_length) :: dump_name
    ! READS SPECIES SECTION OF INPUT FILE
    IF (.NOT. l_pdumps_allocated) THEN
      npdumps=0
      ALLOCATE(particle_dumps(1:nspecies_max))
      l_pdumps_allocated=.TRUE.
    ENDIF
    npdumps=npdumps+1
    dp => particle_dumps(npdumps)
    ! minimal init for filter attributes
    dp%ispecies   = -1
    dp%dump_x_min = xmin
    dp%dump_x_max = xmax
    dp%dump_y_min = ymin
    dp%dump_y_max = ymax
    dp%dump_z_min = zmin
    dp%dump_z_max = zmax
    dp%dump_ux_min = 1e7
    dp%dump_ux_max = 1e9
    dp%dump_uy_min = 1e7
    dp%dump_uy_max = 1e9
    dp%dump_uz_min = 1e7
    dp%dump_uz_max = 1e9
    dp%diag_period = -1
    end_section=.FALSE.
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'dump_x_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_x_min
      ELSE IF (INDEX(buffer, 'dump_x_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_x_max
      ELSE IF (INDEX(buffer, 'species_name') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dump_name
        DO ispecies=1, nspecies
          IF (INDEX(dump_name, species_parray(ispecies)%name) .GT. 0) THEN
            dp%ispecies=ispecies
            EXIT
          ENDIF
        END DO
        IF (dp%ispecies .EQ. -1) THEN
          WRITE(0, *) "ERROR IN SPECIES NAME PARTICLE DUMP SECTION"

        ENDIF
      ELSE IF (INDEX(buffer, 'dump_y_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_y_min
      ELSEIF (INDEX(buffer, 'dump_y_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_y_max
      ELSE IF (INDEX(buffer, 'dump_z_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_z_min
      ELSEIF (INDEX(buffer, 'dump_z_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_z_max
      ELSE IF (INDEX(buffer, 'dump_ux_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_ux_min
      ELSE IF (INDEX(buffer, 'dump_ux_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_ux_max
      ELSE IF (INDEX(buffer, 'dump_uy_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_uy_min
      ELSEIF (INDEX(buffer, 'dump_uy_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_uy_max
      ELSE IF (INDEX(buffer, 'dump_uz_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_uz_min
      ELSEIF (INDEX(buffer, 'dump_uz_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) dp%dump_uz_max
      ELSEIF (INDEX(buffer, 'diag_period') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') dp%diag_period
      ELSE IF (INDEX(buffer, 'end::particle_dump') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_particle_dumps_section


  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that field output parameters section
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_output_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.
    ! READS GRID SECTION OF INPUT FILE
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'output_frequency') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') output_frequency
      ELSE IF (INDEX(buffer, 'output_step_min') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') output_step_min
      ELSE IF (INDEX(buffer, 'output_step_max') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') output_step_max
      ELSEIF (INDEX(buffer, 'ex') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_ex
      ELSE IF (INDEX(buffer, 'ey') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_ey
      ELSE IF (INDEX(buffer, 'ez') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_ez
      ELSEIF (INDEX(buffer, 'bx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_bx
      ELSE IF (INDEX(buffer, 'by') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_by
      ELSE IF (INDEX(buffer, 'bz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_bz
      ELSEIF (INDEX(buffer, 'jx') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_jx
      ELSE IF (INDEX(buffer, 'jy') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_jy
      ELSE IF (INDEX(buffer, 'jz') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_jz
      ELSEIF (INDEX(buffer, 'rho') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_rho
      ELSE IF (INDEX(buffer, 'dive') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_dive
      ELSE IF (INDEX(buffer, 'divj') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_divj
      ELSE IF (INDEX(buffer, 'divb') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') c_output_divb
      ELSE IF (INDEX(buffer, 'end::output') .GT. 0) THEN
        end_section =.TRUE.
      END IF
    END DO
    RETURN
  END SUBROUTINE read_output_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine that reads parameters for temporal diagnistics in the input file
  !> Temporal diagnostics are the temporal evolution of some quantities
  !> such as the particle energies and the field energies.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE read_temporal_output_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section = .FALSE.

    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      !WRITE(0, *), TRIM(ADJUSTL(buffer))
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'frequency') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_frequency
      ELSE IF (INDEX(buffer, 'format') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_format
      ELSE IF (INDEX(buffer, 'kinE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(1)
      ELSE IF (INDEX(buffer, 'exE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(2)
      ELSE IF (INDEX(buffer, 'eyE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(3)
      ELSE IF (INDEX(buffer, 'ezE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(4)
      ELSE IF (INDEX(buffer, 'bxE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(5)
      ELSE IF (INDEX(buffer, 'byE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(6)
      ELSE IF (INDEX(buffer, 'bzE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(7)
      ELSE IF (INDEX(buffer, 'divE-rho') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(8)
      ELSE IF (INDEX(buffer, 'rho') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(9)
      ELSE IF (INDEX(buffer, 'divE') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(10)
      ENDIF
    ENDDO
    RETURN
  END SUBROUTINE read_temporal_output_section

  ! ______________________________________________________________________________________
  !> @brief
  !> Initialization Laser (and antenna) section.
  !> @author
  !> Haithem Kallala
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE read_antenna_section
    INTEGER :: ix = 0
    LOGICAL(lp)  :: end_section 
    TYPE(particle_species), POINTER :: curr, curr2

    IF (.NOT. l_species_allocated) THEN
      nspecies=0
      ALLOCATE(species_parray(1:nspecies_max))
      l_species_allocated=.TRUE.
    ENDIF
    nspecies = nspecies+1
    curr => species_parray(nspecies)
    ! minimal init for species attributes
    curr%charge = 1._num
    curr%mass = emass
    curr%name='laser_antenna'
    curr%nppcell = 1
    curr%x_min = 0._num
    curr%x_max = 0._num
    curr%y_min = 0._num
    curr%y_max = 0._num
    curr%z_min = 0._num
    curr%z_max = 0._num
    curr%vdrift_x =0._num
    curr%vdrift_y =0._num
    curr%vdrift_z =0._num
    curr%vth_x =0._num
    curr%vth_y =0._num
    curr%vth_z =0._num
    curr%sorting_period = 0
    curr%sorting_start = 0
    curr%species_npart=0
    curr%ldodepos = .TRUE.
    ! --- Init default value for antenna params
    curr%is_antenna=.TRUE.
    curr%antenna_params%polangle = 0._num
    curr%antenna_params%spot_x = 0._num
    curr%antenna_params%spot_y = 0._num
    curr%antenna_params%spot_z = 0._num
    curr%antenna_params%lambda_laser = 0._num
    curr%antenna_params%laser_ctau = 0._num
    curr%antenna_params%laser_a_1 = 0._num
    curr%antenna_params%laser_a_2 = 0._num
    curr%antenna_params%laser_w0 = 0._num
    curr%antenna_params%temporal_order = 2
    curr%antenna_params%t_peak = 0._num
    curr%antenna_params%time_window = 0_idp
    curr%antenna_params%vector = 0._num 
    curr%antenna_params%polvector1 = 0._num
    curr%antenna_params%polvector2 = 0._num
     
    end_section = .FALSE.
    DO WHILE((.NOT. end_section) .AND. (ios==0))
      READ(fh_input, '(A)', iostat=ios) buffer
      IF (INDEX(buffer, '#') .GT. 0) THEN
        CYCLE
      ENDIF
      IF (INDEX(buffer, 'vector_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%vector(1)
      ELSE IF (INDEX(buffer, 'vector_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%vector(2)
      ELSE IF (INDEX(buffer, 'vector_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%vector(3)
      ELSE IF (INDEX(buffer, 'spot_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%spot_x
      ELSE IF (INDEX(buffer, 'spot_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%spot_y
      ELSE IF (INDEX(buffer, 'spot_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%spot_z
      ELSE IF (INDEX(buffer, 'lambda_laser') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%lambda_laser
      ELSE IF (INDEX(buffer, 'pvec_x') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%polvector1(1)
      ELSE IF (INDEX(buffer, 'pvec_y') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%polvector1(2)
      ELSE IF (INDEX(buffer, 'pvec_z') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%polvector1(3)
      ELSE IF (INDEX(buffer, 'laser_ctau') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%laser_ctau
      ELSE IF (INDEX(buffer, 't_peak') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%t_peak
      ELSE IF (INDEX(buffer, 'laser_a_1') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%laser_a_1
      ELSE IF (INDEX(buffer, 'laser_a_2') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%laser_a_2
      ELSE IF (INDEX(buffer, 'laser_w0') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%laser_w0
      ELSE IF (INDEX(buffer, 'polangle') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), *) curr%antenna_params%polangle
      ELSE IF (INDEX(buffer, 'temporal_order') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)') curr%antenna_params%temporal_order
      ELSE IF (INDEX(buffer, 'window') .GT. 0) THEN
        ix = INDEX(buffer, "=")
        READ(buffer(ix+1:string_length), '(i10)')curr%antenna_params%time_window
      ELSE IF (INDEX(buffer, 'end::antenna') .GT. 0) THEN
        end_section =.TRUE.
      ENDIF
    ENDDO

    ! Allocate second species for laser antenna 
    ! This species is equivalent to the first one, except it has an 
    ! oppisite charge
    nspecies = nspecies+1
    curr2 => species_parray(nspecies)
    ! minimal init for species attributes
    curr2%charge = -1._num
    curr2%mass = emass
    curr2%name='laser_antenna'
    curr2%nppcell = 1
    curr2%x_min = curr%x_min
    curr2%x_max = curr%x_max
    curr2%y_min = curr%y_min
    curr2%y_max = curr%y_max
    curr2%z_min = curr%z_min
    curr2%z_max = curr%z_max
    curr2%vdrift_x = curr%vdrift_x
    curr2%vdrift_y = curr%vdrift_y 
    curr2%vdrift_z =curr%vdrift_z
    curr2%vth_x =curr%vth_x
    curr2%vth_y =curr%vth_y
    curr2%vth_z =curr%vth_z
    curr2%sorting_period = curr%sorting_period
    curr2%sorting_start = curr%sorting_start
    curr2%species_npart=curr%species_npart
    curr2%ldodepos = .TRUE.
    ! --- Init default value for antenna params
    curr2%is_antenna=.TRUE.
    curr2%antenna_params%polangle = curr%antenna_params%polangle
    curr2%antenna_params%spot_x = curr%antenna_params%spot_x
    curr2%antenna_params%spot_y = curr%antenna_params%spot_y
    curr2%antenna_params%spot_z = curr%antenna_params%spot_z
    curr2%antenna_params%lambda_laser = curr%antenna_params%lambda_laser
    curr2%antenna_params%laser_ctau = curr%antenna_params%laser_ctau
    curr2%antenna_params%laser_a_1 = curr%antenna_params%laser_a_1
    curr2%antenna_params%laser_a_2 = curr%antenna_params%laser_a_2
    curr2%antenna_params%laser_w0 = curr%antenna_params%laser_w0
    curr2%antenna_params%temporal_order = curr%antenna_params%temporal_order
    curr2%antenna_params%t_peak = curr%antenna_params%t_peak
    curr2%antenna_params%time_window =  curr%antenna_params%time_window 
    curr2%antenna_params%vector = curr%antenna_params%vector
    curr2%antenna_params%polvector1 = curr%antenna_params%polvector1
    curr2%antenna_params%polvector2 =  curr%antenna_params%polvector2
    RETURN
  END SUBROUTINE read_antenna_section
  ! ______________________________________________________________________________________
  !> @brief
  !> Initialization of the species section and arrays.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE init_species_section

    IF (.NOT. l_species_allocated) THEN
      nspecies=0
      ALLOCATE(species_parray(1:nspecies_max))
      l_species_allocated=.TRUE.
    ENDIF
  END SUBROUTINE init_species_section

END MODULE control_file

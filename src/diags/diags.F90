! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016, ! The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! If you have questions about your rights to use or distribute this software, ! please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a paid-up, ! nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
! copies to the public, prepare derivative works, and perform publicly and display
! publicly, and to permit other to do so.
!
! DIAGS.F90
!
! Purpose:
! This file contains routines for processing the diagnostics
!
! Authors
! Henri Vincenti, ! Mathieu Lobet
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> This module contains useful diagnostics to test code correctness.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
MODULE diagnostics
  ! ________________________________________________________________________________________
  
  USE constants
  USE mpi
  IMPLICIT NONE
  
  
  CONTAINS
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Computes derived physical quantities from simulation
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE calc_diags
    ! ____________________________________________________________________________________
    USE fields
    USE field_boundary
    USE particle_boundary
    USE particles
    USE params
    USE shared_data
    USE tiling
    USE time_stat
    IMPLICIT NONE
    
    REAL(num) :: tmptime
    
#if defined(DEBUG)
    WRITE(0, *) "Calc_diags: start"
#endif
    IF (l_plasma) THEN
      ! - Computes total charge density
      CALL pxrdepose_rho_on_grid()
      
      ! - Charge boundary conditions
      CALL charge_bcs()
    ENDIF
    ! - Computes electric field divergence on grid at n+1
    IF (.False.) THEN
      
      IF (it.ge.timestat_itstart) THEN
        tmptime = MPI_WTIME()
      ENDIF
      
      IF (.not.(divE_computed))  then
        CALL calc_field_div(dive, ex, ey, ez, nx, ny, nz, nxguards, nyguards,         &
        nzguards, dx, dy, dz)
        divE_computed = .true.
      ENDIF
      
      ! Get the total number of particles
      !CALL get_tot_number_of_particles(ntot)
      
      IF (it.ge.timestat_itstart) THEN
        localtimes(9) = localtimes(9) + (MPI_WTIME() - tmptime)
      ENDIF
      
    ENDIF
    
#if defined(DEBUG)
    WRITE(0, *) "Calc_diags: stop"
#endif
    
  END SUBROUTINE calc_diags
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Computes field divergence.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE calc_field_div(divee, eex, eey, eez, nx, ny, nz, nxguard, nyguard,       &
  nzguard, dx, dy, dz)
    ! ____________________________________________________________________________________
    IMPLICIT NONE
    INTEGER(idp) ::  j, k, l
    INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                    &
    -nzguard:nz+nzguard), intent(in) :: eex, eey, eez
    REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,                    &
    -nzguard:nz+nzguard), intent(in out) :: divee
    REAL(num)    :: dx, dy, dz, invdx, invdy, invdz
    
    invdx=1.0_num/dx
    invdy=1.0_num/dy
    invdz=1.0_num/dz
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, j, k)
    !$OMP DO COLLAPSE(3)
    DO l = 0, nz
      DO k = 0, ny
        DO j = 0, nx
          divee(j, k, l) = invdx*(eex(j, k, l)-eex(j-1, k, l))+ invdy*(eey(j, k,      &
          l)-eey(j, k-1, l))+invdz*(eez(j, k, l)-eez(j, k, l-1)) 
        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
  END SUBROUTINE calc_field_div
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Initialization of the different diags.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE init_diags
    ! ____________________________________________________________________________________
    USE shared_data
    IMPLICIT NONE
    
    IF (rank.eq.0) THEN
      WRITE(0, *)
      WRITE(0, *) ' Initilization of the diags'
    ENDIF
    
    ! Creation of the result folder
#if (defined(VTUNE) || defined(SDE) || defined(DFP) || defined(ALLINEA))
#else
    IF (rank.eq.0) CALL system('mkdir RESULTS')
    IF (rank.eq.0) CALL system('rm RESULTS/*')
#endif
    ! Initialization of the temporal diags
    CALL init_temp_diags
    
    ! Init time statistics
    CALL init_time_stat_output
    
    IF (rank.eq.0) WRITE(0, *)
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Init temporal diags
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE init_temp_diags
    ! ____________________________________________________________________________________
    USE output_data
    USE particle_properties
    USE shared_data
    USE params
    
    IMPLICIT NONE
    
    INTEGER :: i
    
    IF (temdiag_frequency.gt.0) THEN
      
      i = 1
      temdiag_nb_part = 0
      temdiag_nb_field = 0
      temdiag_nb = 0
      
      ! Determine the number of diags for the particles
      ! Kinetic energy
      if (temdiag_act_list(1).gt.0) then
        temdiag_i_list(1) = i
        temdiag_nb_values(1) = nspecies
        temdiag_name_list(1) = 'kinE'
        i = i + nspecies
        temdiag_nb_part = temdiag_nb_part + nspecies
        temdiag_nb = temdiag_nb + 1
      end if
      ! Determine the number of diags for the fields
      ! Ex energy
      if (temdiag_act_list(2).gt.0) then
        temdiag_i_list(2) = i
        temdiag_name_list(2) = 'exE'
        temdiag_nb_values(2) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Ey energy
      if (temdiag_act_list(3).gt.0) then
        temdiag_i_list(3) = i
        temdiag_name_list(3) = 'eyE'
        temdiag_nb_values(3) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Ez energy
      if (temdiag_act_list(4).gt.0) then
        temdiag_i_list(4) = i
        temdiag_name_list(4) = 'ezE'
        temdiag_nb_values(4) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Bx energy
      if (temdiag_act_list(5).gt.0) then
        temdiag_i_list(5) = i
        temdiag_name_list(5) = 'bxE'
        temdiag_nb_values(5) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! By energy
      if (temdiag_act_list(6).gt.0) then
        temdiag_i_list(6) = i
        temdiag_name_list(6) = 'byE'
        temdiag_nb_values(6) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Bz energy
      if (temdiag_act_list(7).gt.0) then
        temdiag_i_list(7) = i
        temdiag_name_list(7) = 'bzE'
        temdiag_nb_values(7) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Norn divE*eps0 - rho
      if (temdiag_act_list(8).gt.0) then
        temdiag_i_list(8) = i
        temdiag_name_list(8) = 'divE-rho'
        temdiag_nb_values(8) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Norn divE*eps0 - rho
      if (temdiag_act_list(9).gt.0) then
        temdiag_i_list(9) = i
        temdiag_name_list(9) = 'rho'
        temdiag_nb_values(9) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      ! Norn divE*eps0 - rho
      if (temdiag_act_list(10).gt.0) then
        temdiag_i_list(10) = i
        temdiag_name_list(10) = 'divE'
        temdiag_nb_values(10) = 1
        temdiag_nb_field = temdiag_nb_field + 1
        temdiag_nb = temdiag_nb + 1
        i = i+1
      end if
      temdiag_totvalues = i - 1
      
      if (temdiag_nb.gt.0) then
        if (rank.eq.0) then
          write(0, '(" Initialization of the temporal diagnostics")')
        end if
      end if
      
      ! Each mpi task will write in a given file according to their rank
      IF (nproc.ge.temdiag_nb) then
        IF ((rank.ge.0).and.(rank.le.temdiag_nb)) then
          if (temdiag_act_list(rank+1).gt.0) then
            write(0, '(" Rank ", I3, ", creation of the file ", A30)') rank,          &
            "./RESULTS/"//trim(adjustl(temdiag_name_list(rank+1)))
            if (temdiag_format.eq.1) then
              open(unit=42,                                                           &
              file="./RESULTS/"//trim(adjustl(temdiag_name_list(rank+1))),            &
              ACTION='write', STATUS='new')
              write(42, *) temdiag_nb_values(rank+1), temdiag_frequency*dt
            else
              open(unit=42,                                                           &
              file="./RESULTS/"//trim(adjustl(temdiag_name_list(rank+1))),            &
              FORM="unformatted", ACCESS='stream', STATUS='new') 
              write(42) temdiag_nb_values(rank+1), temdiag_frequency*dt
            end if
            
          end if
        ENDIF
        ! If there is not enough proc, only the first one deals with it
      else
        IF (rank.eq.0) then
          do i= 1, temdiag_nb
            if (temdiag_act_list(i) .gt.0) then
              write(0, '(" Creation of the file ", A30)')                             &
              "./RESULTS/"//trim(adjustl(temdiag_name_list(i)))
              IF (temdiag_format.eq.1) THEN
                open(unit=42+i,                                                       &
                file="./RESULTS/"//trim(adjustl(temdiag_name_list(i))),               &
                ACTION='write', STATUS='new')
                write(42+i, *) temdiag_nb_values(i), temdiag_frequency*dt
              else
                open(unit=42+i,                                                       &
                file="./RESULTS/"//trim(adjustl(temdiag_name_list(i))),               &
                FORM="unformatted", ACCESS='stream', STATUS='NEW') 
                write(42+i) temdiag_nb_values(i), temdiag_frequency*dt
              end if
            end if
          end do
        end if
      end if
      
      CALL MPI_BARRIER(comm, errcode)
      
    ENDIF
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Initialize outputs of the time statistics
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE init_time_stat_output
    ! ____________________________________________________________________________________
    USE time_stat
    USE shared_data
    USE params
    IMPLICIT NONE
    
    INTEGER :: nb_timestat
    
    IF (timestat_activated.gt.0) THEN
      
      nb_timestat = 13
      
      OPEN(unit=41, file="./RESULTS/time_stat", FORM="unformatted", ACCESS='stream')
      WRITE(41) nb_timestat
      WRITE(41) timestat_period*dt
      
      IF (rank.eq.0) WRITE(0, *) ' Initialization of the time statistics output: ',   &
      nb_timestat
      
      itimestat = 1! index in the buffer
      
      ALLOCATE(buffer_timestat(nb_timestat, nbuffertimestat))! Buffer before output
      
    ELSE
      timestat_period = 0
    ENDIF
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> This subroutine determine the total number of particles in the domain
  !> from species of index is.
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @creation
  !> May 2016
  SUBROUTINE get_tot_number_of_particles_from_species(is, nptot)
    ! ____________________________________________________________________________________
    USE particle_tilemodule
    USE particle_speciesmodule
    USE tile_params
    USE particles
    USE mpi_derived_types
    USE shared_data
    USE tiling
    IMPLICIT NONE
    
    INTEGER(idp), INTENT(INOUT) :: nptot
    INTEGER(idp), INTENT(IN) :: is
    INTEGER(idp) :: nptot_loc
    INTEGER(idp) :: ix, iy, iz, ip, n
    TYPE(particle_tile), POINTER :: curr_tile
    TYPE(particle_species), POINTER :: curr
    
    nptot = 0
    
    ! Current species
    curr=>species_parray(is)
    
    ! Loop over the tiles
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) SHARED(curr,        &
    !$OMP ntilex, ntiley, ntilez) PRIVATE(ix, iy, iz, n, is, curr_tile, ip)           &
    !$OMP reduction(+:nptot_loc)   
    DO iz=1, ntilez
      DO iy=1, ntiley
        DO ix=1, ntilex
          
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          nptot_loc = nptot_loc + curr_tile%np_tile(1)
          
        End do
      End do
    End do
    !$OMP END PARALLEL DO
    
    ! All MPI reduction
    call MPI_ALLREDUCE(nptot_loc, nptot, 1_isp, mpidbl, MPI_SUM, comm, errcode)
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> This subroutine determine the total number of particles all species included
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation: May 2016
  SUBROUTINE get_tot_number_of_particles(nptot)
    ! ____________________________________________________________________________________
    USE particle_tilemodule
    USE particle_speciesmodule
    USE tile_params
    USE particles
    USE mpi_derived_types
    USE shared_data
    USE tiling
    IMPLICIT NONE
    
    INTEGER(idp), INTENT(OUT) :: nptot
    INTEGER(idp)              :: is, nptottmp
    
    nptot = 0
    
    DO is=1, nspecies
      nptottmp = 0
      CALL get_tot_number_of_particles_from_species(is, nptottmp)
      nptot = nptot + nptottmp
    ENDDO
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Determine the local kinetic energy for the species ispecies.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE get_loc_kinetic_energy(ispecies, kinetic_energy_loc)
    ! ____________________________________________________________________________________
    USE particle_tilemodule
    USE particle_speciesmodule
    USE tile_params
    USE particles
    USE mpi_derived_types
    USE shared_data
    USE tiling
    IMPLICIT NONE
    INTEGER(idp) :: ispecies
    
    REAL(num) :: kinetic_energy_loc
    INTEGER(idp) :: ix, iy, iz, np, ip, n
    TYPE(particle_tile), POINTER :: curr_tile
    TYPE(particle_species), POINTER :: curr
    INTEGER(idp), parameter :: LVEC2=16
    REAL(num)                          :: partgam
    REAL(num), dimension(:), allocatable :: gaminv
    
    kinetic_energy_loc = 0
    
    ! current species
    curr=>species_parray(ispecies)
    
    ! Loop over the tiles
    !$OMP PARALLEL DEFAULT(NONE) SHARED(curr, ntilex, ntiley, ntilez) PRIVATE(ix, iy, &
    !$OMP iz, ispecies, curr_tile, ip, np, n, gaminv, partgam)                        &
    !$OMP reduction(+:kinetic_energy_loc)    
    
    !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
    DO iz=1, ntilez
      DO iy=1, ntiley
        DO ix=1, ntilex
          
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np = curr_tile%np_tile(1)
          
          allocate(gaminv(np))
          gaminv(1:np) = curr_tile%part_gaminv(1:np)
          
          !write(0, *) " Tile total kinetic energy for tile:", ix, iy, iz
          !write(0, *) " Number of particles:", np
          !write(0, *) " partx:", curr_tile%part_x(1:10)
          !write(0, *) " gaminv:", 1./curr_tile%part_gaminv(1:10)
          
          ! Loop over the particles
          DO ip=1, np, LVEC2
#if defined __INTEL_COMPILER
            !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64, gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
            !$OMP SIMD SAFELEN(LVEC2)
#elif defined __IBMBGQ__
            !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
            !$DIR SIMD
#endif
            DO n=1, MIN(LVEC2, np-ip+1)
              
              partgam = 1./gaminv(ip+(n-1))
              
              kinetic_energy_loc = kinetic_energy_loc +                               &
              (partgam-1.)*curr_tile%pid(ip+(n-1), wpid)
              
            end do
            
          End do
          
          deallocate(gaminv)
          !write(0, *) " Total Local kinetic energy", total_kinetic_energy_loc
          
        End do
      End do
    End do
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Determine the total kinetic energy for the species ispecies.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE get_kinetic_energy(ispecies, total_kinetic_energy)
    ! ____________________________________________________________________________________
    USE particle_tilemodule
    USE particle_speciesmodule
    USE tile_params
    USE particles
    USE mpi_derived_types
    USE shared_data
    USE tiling
    IMPLICIT NONE
    INTEGER(idp) :: ispecies
    REAL(num) :: total_kinetic_energy
    
    REAL(num) :: total_kinetic_energy_loc
    INTEGER(idp) :: ix, iy, iz, np, ip, n
    TYPE(particle_tile), POINTER :: curr_tile
    TYPE(particle_species), POINTER :: curr
    INTEGER(idp), parameter :: LVEC2=16
    REAL(num) :: partgam
    
    ! current species
    curr=>species_parray(ispecies)
    
    ! Loop over the tiles
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) SHARED(curr,        &
    !$OMP ntilex, ntiley, ntilez) PRIVATE(ix, iy, iz, n, ispecies, curr_tile, ip, np, &
    !$OMP partgam) reduction(+:total_kinetic_energy_loc)   
    DO iz=1, ntilez
      DO iy=1, ntiley
        DO ix=1, ntilex
          
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          np = curr_tile%np_tile(1)
          
          !write(0, *) " Tile total kinetic energy for tile:", ix, iy, iz
          !write(0, *) " Number of particles:", np
          !write(0, *) " partx:", curr_tile%part_x(1:10)
          !write(0, *) " gaminv:", 1./curr_tile%part_gaminv(1:10)
          
          ! Loop over the particles (cache)
          DO ip=1, np, LVEC2
#if defined _OPENMP && _OPENMP>=201307
            !$OMP SIMD SAFELEN(LVEC2)
#elif defined __IBMBGQ__
            !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
            !$DIR SIMD
            !DIR ASSUME_ALIGNED partgaminv:64
#endif
            DO n=1, MIN(LVEC2, np-ip+1)
              
              partgam = 1./curr_tile%part_gaminv(ip+(n-1))
              total_kinetic_energy_loc = total_kinetic_energy_loc +                   &
              (partgam-1.)*curr_tile%pid(ip+(n-1), wpid)
              
            end do
            
          End do
          
          !write(0, *) " Total Local kinetic energy", total_kinetic_energy_loc
          
        End do
      End do
    End do
    !$OMP END PARALLEL DO
    
    ! All MPI reduction
    call MPI_ALLREDUCE(total_kinetic_energy_loc, total_kinetic_energy, 1_isp, mpidbl, &
    MPI_SUM, comm, errcode)
    
    !print*, total_kinetic_energy
    
  END SUBROUTINE get_kinetic_energy
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Determine the local field energy for the given field in 2d.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE get_loc_field_energy_2d(field, nx2, nz2, dx2, dz2, nxguard, nzguard,     &
  field_energy)
    ! ____________________________________________________________________________________
    USE constants
    IMPLICIT NONE
    
    ! __ Parameters ________________________________
    INTEGER(idp)                :: nx2, nz2
    INTEGER(idp)                :: nxguard, nzguard
    REAL(num)                   :: field_energy
    INTEGER(idp)                :: j, k, l
    REAL(num)                   :: dx2, dz2
    REAL(num), dimension(-nxguard:nx2+nxguard, 1, -nzguard:nz2+nzguard), intent(in)   &
    :: field
    field_energy = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(2) REDUCTION(+:field_energy)
    do l = 1, nz2
      do j = 1, nx2
        
        field_energy = field_energy + field(j, 1, l)*field(j, 1, l)*0.5
        
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    field_energy = field_energy*dx2*dz2
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Determine the local field energy for the given field.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE get_loc_field_energy(field, nx2, ny2, nz2, dx2, dy2, dz2, nxguard,       &
  nyguard, nzguard, field_energy)
    ! ____________________________________________________________________________________
    USE constants
    IMPLICIT NONE
    INTEGER(idp)     :: nx2, ny2, nz2
    INTEGER(idp)     :: nxguard, nyguard, nzguard
    REAL(num)                   :: field_energy
    INTEGER(idp)                :: j, k, l
    REAL(num)                   :: dx2, dy2, dz2
    REAL(num), dimension(-nxguard:nx2+nxguard, -nyguard:ny2+nyguard,                  &
    -nzguard:nz2+nzguard), intent(in) :: field
    field_energy = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(3) REDUCTION(+:field_energy)
    
    do l = 1, nz2
      do k = 1, ny2
        do j = 1, nx2
          
          field_energy = field_energy + field(j, k, l)*field(j, k, l)*0.5
          
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    field_energy = field_energy*dx2*dy2*dz2
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Determine the total field energy for the given field
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE get_field_energy_2d(field, nx2, nz2, dx2, dz2, nxguard, nzguard,         &
  field_energy)
    ! ____________________________________________________________________________________
    USE constants
    USE mpi_derived_types
    USE mpi_type_constants
    USE shared_data
    
    ! __ Parameters _____________________________________
    IMPLICIT NONE
    INTEGER(idp)                :: nx2, nz2
    INTEGER(idp)                :: nxguard, nzguard
    REAL(num), dimension(-nxguard:nx2+nxguard, 1, -nzguard:nz2+nzguard) :: field
    REAL(num)                   :: field_energy, field_energy_loc
    INTEGER(idp)                :: j, k, l
    REAL(num)                   :: dx2, dz2
    
    field_energy = 0
    field_energy_loc = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(2) REDUCTION(+:field_energy_loc)
    do l = 1, nz2
      do j = 1, nx2
        field_energy_loc = field_energy_loc + field(j, 1, l)*field(j, 1, l)*0.5
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! All MPI reduction
    call MPI_ALLREDUCE(field_energy_loc, field_energy, 1_isp, mpidbl, MPI_SUM, comm,  &
    errcode)
    
    field_energy = field_energy*dx2*dz2
    
  END SUBROUTINE
  
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Get the energy of the field component field
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 09.29.2016
  !
  !> @param[in] field array for a given field component
  !> @param[in] nx2, ny2, nz2 number of cells
  !> @param[in] dx2, dy2, dz2 space step in each direction
  !> @param[in] nxguard, nyguard, nzguard number of guard cells
  !> @param[out] field_energy energy og the corresponding field
  !
  SUBROUTINE get_field_energy(field, nx2, ny2, nz2, dx2, dy2, dz2, nxguard, nyguard,  &
  nzguard, field_energy)
    ! ____________________________________________________________________________________
    USE constants
    USE mpi_derived_types
    USE mpi_type_constants
    USE shared_data
    IMPLICIT NONE
    INTEGER(idp)                :: nx2, ny2, nz2
    INTEGER(idp)                :: nxguard, nyguard, nzguard
    REAL(num), dimension(-nxguard:nx2+nxguard, -nyguard:ny2+nyguard,                  &
    -nzguard:nz2+nzguard) :: field
    REAL(num)                   :: field_energy, field_energy_loc
    INTEGER(idp)                :: j, k, l
    REAL(num)                   :: dx2, dy2, dz2
    
    field_energy = 0
    field_energy_loc = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(3) REDUCTION(+:field_energy_loc)
    do l = 1, nz2
      do k = 1, ny2
        do j = 1, nx2
          
          field_energy_loc = field_energy_loc + field(j, k, l)*field(j, k, l)*0.5
          
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! All MPI reduction
    call MPI_ALLREDUCE(field_energy_loc, field_energy, 1_isp, mpidbl, MPI_SUM, comm,  &
    errcode)
    
    field_energy = field_energy*dx2*dy2*dz2
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Compute norm of dF/dt = divE -rho/eps0 local to the MPI domain (parallel function)
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  SUBROUTINE get_loc_norm_divErho(divee2, rho2, nx2, ny2, nz2, nxguard, nyguard,      &
  nzguard, norm)
    ! ____________________________________________________________________________________
    USE mpi_derived_types
    USE mpi_type_constants
    USE shared_data
    USE constants
    IMPLICIT NONE
    
    INTEGER(idp)                :: j, k, l
    INTEGER(idp)                :: nx2, ny2, nz2, nxguard, nyguard, nzguard
    REAL(num), dimension(-nxguard:nx2+nxguard, -nyguard:ny2+nyguard,                  &
    -nzguard:nz2+nzguard), intent(in) :: divee2, rho2
    REAL(num), intent(out)      :: norm
    
    norm = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(3) REDUCTION(+:norm)
    do l = 1, nz2
      do k = 1, ny2
        do j = 1, nx2
          
          norm = norm + (divee2(j, k, l)*eps0 - rho2(j, k, l))**2
          
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Compute the  square of the norm of array local to the MPI domain (OpenMP parallel function)
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  SUBROUTINE get_loc_norm_2(array, nx2, ny2, nz2, nxguard, nyguard, nzguard, norm)
    ! ____________________________________________________________________________________
    USE mpi_derived_types
    USE mpi_type_constants
    USE shared_data
    USE constants
    IMPLICIT NONE
    
    INTEGER(idp)                :: j, k, l
    INTEGER(idp)                :: nx2, ny2, nz2, nxguard, nyguard, nzguard
    REAL(num), dimension(-nxguard:nx2+nxguard, -nyguard:ny2+nyguard,                  &
    -nzguard:nz2+nzguard), intent(in) :: array
    REAL(num), intent(out)      :: norm
    
    norm = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(3) REDUCTION(+:norm)
    do l = 1, nz2
      do k = 1, ny2
        do j = 1, nx2
          
          norm = norm+ (array(j, k, l))**2
          
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
  END SUBROUTINE
  
  ! ____________________________________________________________________________________
  !> @brief
  !> Compute norm of dF/dt = divE -rho/eps0 (parallel function)
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  SUBROUTINE get_norm_divErho(divee2, rho2, nx2, ny2, nz2, nxguard, nyguard, nzguard, &
  norm)
    ! ____________________________________________________________________________________
    USE mpi_derived_types
    USE mpi_type_constants
    USE shared_data
    USE constants
    IMPLICIT NONE
    
    INTEGER(idp)                :: j, k, l
    INTEGER(idp)                :: nx2, ny2, nz2, nxguard, nyguard, nzguard
    REAL(num), dimension(-nxguard:nx2+nxguard, -nyguard:ny2+nyguard,                  &
    -nzguard:nz2+nzguard), intent(in) :: divee2, rho2
    REAL(num)                   :: norm_loc
    REAL(num), intent(out)      :: norm
    
    norm_loc = 0
    norm = 0
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
    !$OMP DO COLLAPSE(3) REDUCTION(+:norm_loc)
    do l = 1, nz2
      do k = 1, ny2
        do j = 1, nx2
          
          norm_loc = norm_loc + (divee2(j, k, l)*eps0 - rho2(j, k, l))**2
          
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! All MPI reduction
    call MPI_ALLREDUCE(norm_loc, norm, 1_isp, mpidbl, MPI_SUM, comm, errcode)
    
    norm = sqrt(norm)
    
  END SUBROUTINE
END MODULE diagnostics

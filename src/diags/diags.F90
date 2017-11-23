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
! ________________________________________________________________________________________
MODULE diagnostics
  USE constants
  USE mpi
  IMPLICIT NONE

  CONTAINS

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE calc_diags
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE calc_field_div(divee, eex, eey, eez, nx, ny, nz, nxguard, nyguard,       &
    nzguard, dx, dy, dz)
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


  ! ______________________________________________________________________________________
  !> @brief
  !> Computes B field divergence.
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  !
  ! ______________________________________________________________________________________
  SUBROUTINE calc_field_divB(divbb, bbx, bby, bbz, nx, ny, nz, nxguard, nyguard,&
    nzguard, dx, dy, dz)
    IMPLICIT NONE
    INTEGER(idp) ::  j, k, l
    INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,&
    -nzguard:nz+nzguard), intent(in) :: bbx, bby, bbz
    REAL(num), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,&
    -nzguard:nz+nzguard), intent(in out) :: divbb
    REAL(num)    :: dx, dy, dz, invdx, invdy, invdz

    invdx=1.0_num/dx
    invdy=1.0_num/dy
    invdz=1.0_num/dz
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, j, k)
    !$OMP DO COLLAPSE(3)
    DO l = 1, nz
      DO k = 1, ny
        DO j = 1, nx
          divbb(j, k, l) = invdx*(bbx(j+1, k, l)-bbx(j, k, l))+ invdy*(bby(j, k+1,& 
          l)-bby(j, k, l))+invdz*(bbz(j, k, l+1)-bbz(j, k, l))
        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

  END SUBROUTINE calc_field_divB



  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE init_diags
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE init_temp_diags
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
      CALL MPI_BARRIER(comm,errcode)
      IF (nproc.ge.temdiag_nb) then
        IF ((rank.ge.0).and.(rank.lt.temdiag_nb)) then
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE init_time_stat_output
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

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine determine the total number of particles in one MPI domain
  !> from species of index is.
  !
  !> @author
  !> Guillaume Blaclard
  !
  !> @creation
  !> June 2017
  ! ______________________________________________________________________________________
  SUBROUTINE get_local_number_of_particles_from_species(is, nptot_loc)
    USE particle_tilemodule
    USE particle_speciesmodule
    USE tile_params
    USE particles
    USE mpi_derived_types
    USE shared_data
    USE tiling
    IMPLICIT NONE

    INTEGER(idp), DIMENSION(1), INTENT(INOUT) :: nptot_loc
    INTEGER(idp), INTENT(IN) :: is
    INTEGER(idp) :: ix,iy,iz
    TYPE(particle_tile), POINTER :: curr_tile
    TYPE(particle_species), POINTER :: curr

    ! Current species
    curr=>species_parray(is)
    nptot_loc(1) = 0_idp
    ! Loop over the tiles
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
    !$OMP SHARED(curr, ntilex,ntiley,ntilez) &
    !$OMP PRIVATE(ix,iy,iz,is,curr_tile) &
    !$OMP reduction(+:nptot_loc)
    DO iz=1, ntilez
      DO iy=1, ntiley
        DO ix=1, ntilex

          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          nptot_loc(1) = nptot_loc(1) + curr_tile%np_tile(1)

        End do
      End do
    End do
    !$OMP END PARALLEL DO

  END SUBROUTINE get_local_number_of_particles_from_species


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine determine the total number of particles in the MPI domain
  !> from species of index is.
  !
  !> @author
  !> Mathieu Lobet
  !> Guillaume Blaclard
  !
  !> @date
  !> Creation: May 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_tot_number_of_particles_from_species(is, nptot)
    USE particle_tilemodule
    USE particle_speciesmodule
    USE tile_params
    USE particles
    USE mpi_derived_types
    USE shared_data
    USE tiling
    IMPLICIT NONE

    INTEGER(idp), DIMENSION(1), INTENT(INOUT) :: nptot
    INTEGER(idp), INTENT(IN) :: is
    INTEGER(idp), DIMENSION(1) :: nptot_loc

    nptot(1) = 0_idp
    nptot_loc(1) = 0_idp
    CALL get_local_number_of_particles_from_species(is, nptot_loc)

    ! All MPI reduction
    CALL MPI_ALLREDUCE(nptot_loc(1), nptot(1), 1_isp, MPI_INTEGER8, MPI_SUM,          &
    comm,errcode)

  END SUBROUTINE get_tot_number_of_particles_from_species

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine returns a given position, a momentum or a field of a particle in
  !> the domain from species of index ispecies.
  !> For x, y, z, ux, uy, uz, quantity must be egal respectively to 1, 2, 3, 4, 5, 6.
  !> For ex, ey, ez, bx, by, bz, quantity must be egal respectively to 7, 8, 9, 10, 11
  !> and 12.
  !
  !> @author
  !> Guillaume Blaclard
  !
  !> @creation
  !> June 2017
  ! ______________________________________________________________________________________
  SUBROUTINE getquantity(ispecies, quantity, nptot, quantityarray)
      USE particle_tilemodule
      USE particle_speciesmodule
      USE tile_params
      USE particles
      USE tiling
      IMPLICIT NONE

      INTEGER(idp), INTENT(IN) :: ispecies
      INTEGER(idp), INTENT(IN) :: nptot
      REAL(num), dimension(nptot), INTENT(OUT) :: quantityarray
      INTEGER(idp), INTENT(IN) :: quantity
      INTEGER(idp) :: ix, iy, iz, np
      INTEGER(idp) :: compt
      TYPE(particle_tile), POINTER :: curr_tile
      TYPE(particle_species), POINTER :: curr

      curr=>species_parray(ispecies)
      compt = 1

      ! Loop over the tiles
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
            np = curr_tile%np_tile(1)

            SELECT CASE (quantity)
              CASE  (1)
                quantityarray(compt:compt+np-1) = curr_tile%part_x(1:np)
              CASE  (2)
                quantityarray(compt:compt+np-1) = curr_tile%part_y(1:np)
              CASE  (3)
                quantityarray(compt:compt+np-1) = curr_tile%part_z(1:np)
              CASE  (4)
                quantityarray(compt:compt+np-1) = curr_tile%part_ux(1:np)
              CASE  (5)
                quantityarray(compt:compt+np-1) = curr_tile%part_uy(1:np)
              CASE  (6)
                quantityarray(compt:compt+np-1) = curr_tile%part_uz(1:np)
              CASE  (7)
                quantityarray(compt:compt+np-1) = curr_tile%part_ex(1:np)
              CASE  (8)
                quantityarray(compt:compt+np-1) = curr_tile%part_ey(1:np)
              CASE  (9)
                quantityarray(compt:compt+np-1) = curr_tile%part_ez(1:np)
              CASE (10)
                quantityarray(compt:compt+np-1) = curr_tile%part_bx(1:np)
              CASE (11)
                quantityarray(compt:compt+np-1) = curr_tile%part_by(1:np)
              CASE (12)
                quantityarray(compt:compt+np-1) = curr_tile%part_bz(1:np)
            END SELECT

            compt = compt + np
          END DO
        END DO
      END DO

  END SUBROUTINE getquantity

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine returns the given variable from pid arrays of particles in
  !> the domain from species of index ispecies.
  !> Then quantity_pid corresponds to something like xoldpid, wpid etc..
  !
  !> @author
  !> Guillaume Blaclard
  !
  !> @creation
  !> June 2017
  ! ______________________________________________________________________________________
  SUBROUTINE getquantity_pid(ispecies, quantitypid, nptot, quantityarray)
      USE particle_tilemodule
      USE particle_speciesmodule
      USE tile_params
      USE particles
      USE tiling
      IMPLICIT NONE

      INTEGER(idp), INTENT(IN) :: ispecies
      INTEGER(idp), INTENT(IN) :: nptot
      REAL(num), dimension(nptot), INTENT(OUT) :: quantityarray
      INTEGER(idp), INTENT(IN) :: quantitypid
      INTEGER(idp) :: ix, iy, iz, np
      INTEGER(idp) :: compt
      TYPE(particle_tile), POINTER :: curr_tile
      TYPE(particle_species), POINTER :: curr

      curr=>species_parray(ispecies)
      compt = 1

      ! Loop over the tiles
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            np = curr_tile%np_tile(1)
            quantityarray(compt:compt+np-1) = curr_tile%pid(1:np, quantitypid)
            compt = compt + np
          END DO
        END DO
      END DO

  END SUBROUTINE getquantity_pid

  ! ______________________________________________________________________________________
  !> @brief
  !> Load particle quantity and select particles which cross a moving plane and return
  !> their indexes and the number of True. Assume that the inital total number of particles
  !> is known.
  !
  !> @author
  !> Guillaume Blaclard
  !
  !> @creation
  !> June 2017
  ! ______________________________________________________________________________________
  SUBROUTINE load_test_plane(nptot, plane_position, plane_position_old,                 &
  plane_normal_vector, particle_x, particle_y, particle_z, particle_x_old,              &
  particle_y_old, particle_z_old, npnew, index_particle, particle_relative_position,    &
  particle_relative_position_old)
      IMPLICIT NONE

      INTEGER(idp), INTENT(IN) :: nptot
      REAL(num), DIMENSION(3), INTENT(IN) :: plane_position
      REAL(num), DIMENSION(3), INTENT(IN) :: plane_position_old
      REAL(num), DIMENSION(3), INTENT(IN) :: plane_normal_vector
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_x
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_y
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_z
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_x_old
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_y_old
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_z_old
      INTEGER(idp), INTENT(OUT) :: npnew
      INTEGER(idp), DIMENSION(nptot), INTENT(OUT) :: index_particle
      REAL(num), DIMENSION(nptot), INTENT(OUT) :: particle_relative_position
      REAL(num), DIMENSION(nptot), INTENT(OUT) :: particle_relative_position_old

      INTEGER(idp) :: i
      npnew = 0_idp
      DO i=1, nptot
        particle_relative_position_old(i) =                                             &
        plane_normal_vector(1)*(particle_x_old(i) - plane_position_old(1))              &
        + plane_normal_vector(2)*(particle_y_old(i) - plane_position_old(2))            &
        + plane_normal_vector(3)*(particle_z_old(i) - plane_position_old(3))

        particle_relative_position(i) =                                                 &
        plane_normal_vector(1)*(particle_x(i) - plane_position(1))                      &
        + plane_normal_vector(2)*(particle_y(i) - plane_position(2))                    &
        + plane_normal_vector(3)*(particle_z(i) - plane_position(3))

        IF ( particle_relative_position_old(i) <= 0 .AND.                               &
        particle_relative_position(i) >0) THEN
          index_particle(i) = 1
          npnew = npnew +1

        END IF
      END DO

  END SUBROUTINE load_test_plane

  ! ______________________________________________________________________________________
  !> @brief
  !> Select a particle quantity thanks to the mask index_particle.
  !
  !> @author
  !> Guillaume Blaclard
  !
  !> @creation
  !> June 2017
  ! ______________________________________________________________________________________
  SUBROUTINE select_quantity(nptot, npnew, index_particle, quantity_array_tot,          &
  quantity_array_new)
      IMPLICIT NONE

      INTEGER(idp), INTENT(IN) :: nptot
      INTEGER(idp), INTENT(IN) :: npnew
      INTEGER(idp), DIMENSION(nptot), INTENT(IN) :: index_particle
      REAL(num), DIMENSION(nptot), INTENT(IN) :: quantity_array_tot
      REAL(num), DIMENSION(npnew), INTENT(OUT) :: quantity_array_new
      INTEGER(idp) :: compt
      INTEGER(idp) :: i

      compt = 1
      DO i=1, nptot
        IF ( index_particle(i) == 1 ) THEN
          quantity_array_new(compt) = quantity_array_tot(i)
          compt = compt + 1
        END IF
      END DO
  END SUBROUTINE select_quantity

  ! ______________________________________________________________________________________
  !> @brief
  !> Return the value of captured_quantity, interpolated value between previous_quantity
  !> and current_quantity on the plane.
  !
  !> @author
  !> Guillaume Blaclard
  !
  !> @creation
  !> June 2017
  ! ______________________________________________________________________________________
  SUBROUTINE interpolate_quantity( index_particle, nptot, npnew, previous_quantity,     &
  current_quantity, particle_relative_position, particle_relative_position_old,         &
  captured_quantity)

      IMPLICIT NONE

      INTEGER(idp), INTENT(IN) :: nptot
      INTEGER(idp), INTENT(IN) :: npnew
      INTEGER(idp), DIMENSION(nptot), INTENT(IN) :: index_particle
      REAL(num), DIMENSION(nptot), INTENT(IN) :: previous_quantity
      REAL(num), DIMENSION(nptot), INTENT(IN) :: current_quantity
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_relative_position
      REAL(num), DIMENSION(nptot), INTENT(IN) :: particle_relative_position_old
      REAL(num), DIMENSION(npnew), INTENT(OUT) :: captured_quantity
      REAL(num), DIMENSION(:), ALLOCATABLE :: previous_quantity_sbs
      REAL(num), DIMENSION(:), ALLOCATABLE :: current_quantity_sbs
      REAL(num), DIMENSION(:), ALLOCATABLE :: particle_relative_position_sbs
      REAL(num), DIMENSION(:), ALLOCATABLE :: particle_relative_position_old_sbs
      REAL(num), DIMENSION(:), ALLOCATABLE :: norm_factor
      REAL(num), DIMENSION(:), ALLOCATABLE :: interp_current
      REAL(num), DIMENSION(:), ALLOCATABLE :: interp_previous
      INTEGER(idp) :: i

      IF (npnew > 0) THEN
        ALLOCATE (particle_relative_position_sbs(npnew))
        ALLOCATE (particle_relative_position_old_sbs(npnew))
        ALLOCATE (previous_quantity_sbs(npnew))
        ALLOCATE (current_quantity_sbs(npnew))
        ALLOCATE (interp_current(npnew))
        ALLOCATE (interp_previous(npnew))
        ALLOCATE (norm_factor(npnew))

        ! Take the index of particles where index_particle is true
        CALL select_quantity(nptot, npnew, index_particle, previous_quantity,           &
        previous_quantity_sbs)
        CALL select_quantity(nptot, npnew, index_particle, current_quantity,            &
        current_quantity_sbs)
        CALL select_quantity(nptot, npnew, index_particle,                              &
        particle_relative_position, particle_relative_position_sbs)
        CALL select_quantity(nptot, npnew, index_particle,                              &
        particle_relative_position_old, particle_relative_position_old_sbs)

        ! Interpolate particle quantity to the time when they cross the plane
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(npnew, norm_factor,           &
        !$OMP particle_relative_position_old_sbs, interp_current, interp_previous,      &
        !$OMP particle_relative_position_sbs, previous_quantity_sbs,                    &
        !$OMP current_quantity_sbs, captured_quantity)
        DO i=1, npnew
          norm_factor(i) = 1 / ( ABS(particle_relative_position_old_sbs(i))             &
                      + particle_relative_position_sbs(i) )
          interp_current(i) = ABS(particle_relative_position_old_sbs(i)) * norm_factor(i)
          interp_previous(i) = particle_relative_position_sbs(i) * norm_factor(i)

          captured_quantity(i) = interp_current(i) * current_quantity_sbs(i)            &
                               + interp_previous(i) * previous_quantity_sbs(i)
        END DO
        !$OMP END PARALLEL DO
        DEALLOCATE (particle_relative_position_sbs)
        DEALLOCATE (particle_relative_position_old_sbs)
        DEALLOCATE (previous_quantity_sbs)
        DEALLOCATE (current_quantity_sbs)
        DEALLOCATE (interp_current)
        DEALLOCATE (interp_previous)
        DEALLOCATE (norm_factor)
      END IF
  END SUBROUTINE interpolate_quantity


  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE get_loc_kinetic_energy(ispecies, kinetic_energy_loc)
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

            END DO

          END DO

          deallocate(gaminv)
          !write(0, *) " Total Local kinetic energy", total_kinetic_energy_loc

        End do
      End do
    End do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE get_kinetic_energy(ispecies, total_kinetic_energy)
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE get_loc_field_energy_2d(field, nx2, nz2, dx2, dz2, nxguard, nzguard,     &
    field_energy)
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE get_loc_field_energy(field, nx2, ny2, nz2, dx2, dy2, dz2, nxguard,       &
    nyguard, nzguard, field_energy)
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE get_field_energy_2d(field, nx2, nz2, dx2, dz2, nxguard, nzguard,         &
    field_energy)
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE get_field_energy(field, nx2, ny2, nz2, dx2, dy2, dz2, nxguard, nyguard,  &
    nzguard, field_energy)
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

  ! ______________________________________________________________________________________
  !> @brief
  !> Compute norm of dF/dt = divE -rho/eps0 local to the MPI domain (parallel function)
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_loc_norm_divErho(divee2, rho2, nx2, ny2, nz2, nxguard, nyguard,      &
    nzguard, norm)
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

  ! ______________________________________________________________________________________
  !> @brief
  !> Compute the  square of the norm of array local to the MPI domain
  !> (OpenMP parallel function)
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_loc_norm_2(array, nx2, ny2, nz2, nxguard, nyguard, nzguard, norm)
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

  ! ______________________________________________________________________________________
  !> @brief
  !> Compute norm of dF/dt = divE -rho/eps0 (parallel function)
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_norm_divErho(divee2, rho2, nx2, ny2, nz2, nxguard, nyguard, nzguard, &
    norm)
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

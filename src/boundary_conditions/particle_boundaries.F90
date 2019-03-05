! ______________________________________________________________________________________
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
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! PARTICLE_BOUNDARIES.F90
!
! Brief description:
! This file contains routines for the particle boundary conditions.
!
! For the moment this module handles:
! periodic external boundary conditions for particles and fields
!
! authors:
! Henri Vincenti
! Mathieu Lobet
!
!
! List of suboutines:
! - particle_bcs_tiles_and_mpi_3d
!
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Module that contains subroutines for boundary conditions on particles.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
MODULE particle_boundary

  USE shared_data
  USE fields
  USE particles
  USE tiling
  USE mpi_derived_types
  USE PICSAR_precision
  USE constants
  USE time_stat
  USE params

  IMPLICIT NONE

  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition routine on particles in 3d
  !
  !> @details
  !> This subroutine is the main one to manage the particle boundary conditions
  !> between MPI domains and between tiles.
  !> The different model of exchanges can be used via the variable partcom
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation: 2015
  !> last modified: 09/12/2016
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs
    USE mpi
#ifdef _OPENMP
    USE omp_lib
#endif
    USE picsar_precision, ONLY: num
    USE time_stat, ONLY: localtimes, timestat_itstart
    IMPLICIT NONE

    REAL(num) :: tdeb, tend
    REAL(num) :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    tdeb=MPI_WTIME()

    ! ___________________________________________
    ! Tile and MPI exchanges are done separately without OpenMP
    IF (partcom.eq.2) THEN

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_tiles: start"
#endif

      SELECT CASE (c_dim)
        ! __________________________
        ! 2D
      CASE(2)
        CALL particle_bcs_tiles_2d()
        ! __________________________
        ! 3D
      CASE DEFAULT
        CALL particle_bcs_tiles()
      END SELECT

      IF (it.ge.timestat_itstart) THEN
        localtimes(11) = localtimes(11) + (MPI_WTIME() - tmptime)
        tmptime = MPI_WTIME()
      ENDIF
      tend = MPI_WTIME()
      local_time_part=local_time_part+(tend-tdeb)

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_mpi: start"
#endif

      IF (mpicom_curr .EQ. 1) THEN
        ! Then exchange particle between MPI domains
        CALL particle_bcs_mpi_blocking()
      ELSE
        CALL particle_bcs_mpi_non_blocking()
      ENDIF

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_mpi: stop"
#endif

      IF (it.ge.timestat_itstart) THEN
        localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)
      ENDIF

      ! ___________________________________________
      ! OpenMP and MPI exchanges are done separately
    ELSE IF (partcom.eq.1) THEN

      ! First exchange particles between tiles (NO MPI at that point)
#if defined(DEBUG)
#ifdef _OPENMP
      WRITE(0, *) "particle_bcs_tiles_openmp: start"
#else
      WRITE(0, *) "particle_bcs_tiles: start"
#endif
#endif
      SELECT CASE (c_dim)
        ! __________________________
        ! 2D
      CASE(2)
#ifdef _OPENMP
        CALL particle_bcs_tiles_2d_openmp()
#else
        CALL particle_bcs_tiles_2d()
#endif
        ! __________________________
        ! 3D
      CASE DEFAULT
#ifdef _OPENMP
        CALL particle_bcs_tiles_openmp()
#else
        CALL particle_bcs_tiles()
#endif
      END SELECT
#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_tiles: stop"
#endif

      IF (it.ge.timestat_itstart) THEN
        localtimes(11) = localtimes(11) + (MPI_WTIME() - tmptime)
        tmptime = MPI_WTIME()
      ENDIF
      tend = MPI_WTIME()
      local_time_part=local_time_part+(tend-tdeb)

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_mpi: start"
#endif
      IF (mpicom_curr .EQ. 1) THEN
        ! Then exchange particle between MPI domains
        CALL particle_bcs_mpi_blocking()
      ELSE
        CALL particle_bcs_mpi_non_blocking()
      ENDIF

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_mpi: stop"
#endif

      IF (it.ge.timestat_itstart) THEN
        localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)
      ENDIF

      ! _____________________________________________
      ! Tile and MPI com in one
      ! Default
    ELSE

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_tiles_and_mpi: start"
#endif

      IF (it.ge.timestat_itstart) THEN
        tmptime = MPI_WTIME()
      ENDIF

      SELECT CASE (c_dim)
        ! __________________________
        ! 2D
      CASE(2)

#ifdef _OPENMP
        CALL particle_bcs_tiles_2d_openmp()
#else
        CALL particle_bcs_tiles_2d()
#endif

        CALL particle_bcs_mpi_non_blocking_2d()

        ! __________________________
        ! 3D
      CASE DEFAULT

        CALL particle_bcs_tiles_and_mpi_3d

      END SELECT

      IF (it.ge.timestat_itstart) THEN
        localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)
      ENDIF

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_tiles_and_mpi: stop"
#endif
    ENDIF

  END SUBROUTINE particle_bcs

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition routine on particles in 2d
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs_2d
    USE mpi
#ifdef _OPENMP
    USE omp_lib
#endif
    USE picsar_precision, ONLY: num
    USE tiling
    USE time_stat, ONLY: localtimes

    IMPLICIT NONE
    REAL(num) :: tdeb, tend
    REAL(num) :: tmptime
    tmptime = MPI_WTIME()
    tdeb=MPI_WTIME()

    ! ___________________________________________
    ! OpenMP and MPI exchanges are done separately
    IF (partcom.eq.1) THEN

      ! First exchange particles between tiles (NO MPI at that point)
#if defined(DEBUG)
#ifdef _OPENMP
      WRITE(0, *) "particle_bcs_tiles_openmp: start"
#else
      WRITE(0, *) "particle_bcs_tiles: start"
#endif
#endif

#ifdef _OPENMP
      CALL particle_bcs_tiles_2d_openmp()
#else
      CALL particle_bcs_tiles_2d()
#endif

      localtimes(11) = localtimes(11) + (MPI_WTIME() - tmptime)
      tmptime = MPI_WTIME()
      tend = MPI_WTIME()
      local_time_part=local_time_part+(tend-tdeb)

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_mpi: start"
#endif

      CALL particle_bcs_mpi_non_blocking_2d()

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_mpi: stop"
#endif

      localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)

      ! _____________________________________________
      ! Tile and MPI com in one
      ! Default
    ELSE

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_tiles_and_mpi: start"
#endif

      tmptime = MPI_WTIME()

#ifdef _OPENMP
      CALL particle_bcs_tiles_2d_openmp()
#else
      CALL particle_bcs_tiles_2d()
#endif

      CALL particle_bcs_mpi_non_blocking_2d()

      localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)

#if defined(DEBUG)
      WRITE(0, *) "particle_bcs_tiles_and_mpi: stop"
#endif
    ENDIF

  END SUBROUTINE particle_bcs_2d


  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition on tiles - 3D version Wihtout OpenMP.
  !
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs_tiles
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(:), ALLOCATABLE :: partpid

    ALLOCATE(partpid(1:npid))
    DO ispecies=1, nspecies! LOOP ON SPECIES
      curr=> species_parray(ispecies)
      ! Get first tiles dimensions (may be different from last tile)
      nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
      ny0_grid_tile = curr%array_of_tiles(1, 1, 1)%ny_grid_tile
      nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile
      DO iz=1, ntilez! LOOP ON TILES
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            nptile=curr_tile%np_tile(1)
            DO i=nptile, 1, -1! LOOP ON PARTICLES
              partx=curr_tile%part_x(i)
              party=curr_tile%part_y(i)
              partz=curr_tile%part_z(i)
              partux=curr_tile%part_ux(i)
              partuy=curr_tile%part_uy(i)
              partuz=curr_tile%part_uz(i)
              gaminv=curr_tile%part_gaminv(i)
              partpid=curr_tile%pid(i, 1:npid)

              ! Case 1: if particle did not leave tile nothing to do
              IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT.                &
              curr_tile%x_tile_max)) .AND. ((party .GE. curr_tile%y_tile_min) .AND.   &
              (party .LT. curr_tile%y_tile_max)) .AND. ((partz .GE.                   &
              curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) CYCLE

              ! Case 2: if particle left MPI domain nothing to do now
              IF ((partx .LT. x_min_local_part) .OR. (partx .GE. x_max_local_part))   &
              CYCLE
              IF ((party .LT. y_min_local_part) .OR. (party .GE. y_max_local_part))   &
              CYCLE
              IF ((partz .LT. z_min_local_part) .OR. (partz .GE. z_max_local_part))   &
              CYCLE

              ! Case 3: particles changed tile. Tranfer particle to new tile
              ! Get new indexes of particle in array of tiles
              indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),       &
              idp)+1, ntilex)
              indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy),       &
              idp)+1, ntiley)
              indz = MIN(FLOOR((partz-(z_min_local)+dz/2_num)/(nz0_grid_tile*dz),     &
              idp)+1, ntilez)
              CALL rm_particle_at_tile(curr, ix, iy, iz, i)
              CALL add_particle_at_tile(curr, indx, indy, indz, partx, party, partz,  &
              partux, partuy, partuz, gaminv, partpid)
            END DO!END LOOP ON PARTICLES
          END DO
        END DO
      END DO! END LOOP ON TILES
    END DO! END LOOP ON SPECIES
    DEALLOCATE(partpid)
  END SUBROUTINE particle_bcs_tiles

  ! ______________________________________________________________________________________
  !> Boundary condition on tiles - 3D version with OpenMP
  !> @brief
  !
  !> @details
  !> This version is efficient when the number of tiles is large
  !> compared to the number of threads.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> 2016
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs_tiles_openmp()
#ifdef _OPENMP
    USE omp_lib
#endif
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(:), ALLOCATABLE :: partpid
    INTEGER(idp) :: nthreads_tot, nthreads_loop1, nthreads_loop2

#ifdef _OPENMP
    nthreads_tot=OMP_GET_MAX_THREADS()
    !nthreads_tot=OMP_GET_NUM_THREADS()
    CALL OMP_SET_NESTED(.TRUE.)
#else
    nthreads_tot=1
#endif

    IF (nthreads_tot .GT. 1) THEN
      nthreads_loop1=MIN(nspecies, nthreads_tot)
      nthreads_loop2=MAX(1_idp, nthreads_tot/nthreads_loop1)
    ELSE
      nthreads_loop1=1
      nthreads_loop2=1
    ENDIF

    ALLOCATE(partpid(npid))
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr, ispecies, nx0_grid_tile,            &
    !$OMP ny0_grid_tile, nz0_grid_tile, ipx, ipy, ipz, partx, party, partz, partux,   &
    !$OMP partuy, partuz, gaminv, partpid, indx, indy, indz, nptile, curr_tile)       &
    !$OMP SHARED(nspecies, npid, nthreads_loop2, species_parray, ntilex, ntiley,      &
    !$OMP ntilez, x_min_local, y_min_local, z_min_local, x_min_local_part,            &
    !$OMP y_min_local_part, z_min_local_part, x_max_local_part, y_max_local_part,     &
    !$OMP z_max_local_part, dx, dy, dz) NUM_THREADS(nthreads_loop1)
    DO ispecies=1, nspecies! LOOP ON SPECIES
      curr=> species_parray(ispecies)
      ! Get first tiles dimensions (may be different from last tile)
      nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
      ny0_grid_tile = curr%array_of_tiles(1, 1, 1)%ny_grid_tile
      nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile
      DO ipz=1, 3
        DO ipy=1, 3
          DO ipx=1, 3
            !$OMP PARALLEL DO DEFAULT(NONE) SHARED(curr, npid, ntilex, ntiley,        &
            !$OMP ntilez, x_min_local, y_min_local, z_min_local, x_min_local_part,    &
            !$OMP y_min_local_part, z_min_local_part, x_max_local_part,               &
            !$OMP y_max_local_part, z_max_local_part, dx, dy, dz, nx0_grid_tile,      &
            !$OMP ny0_grid_tile, nz0_grid_tile) FIRSTPRIVATE(ipx, ipy, ipz)           &
            !$OMP PRIVATE(ix, iy, iz, i, curr_tile, nptile, partx, party, partz,      &
            !$OMP partux, partuy, partuz, gaminv, partpid, indx, indy, indz)          &
            !$OMP COLLAPSE(3) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
            DO iz=ipz, ntilez, 3! LOOP ON TILES
              DO iy=ipy, ntiley, 3
                DO ix=ipx, ntilex, 3
                  curr_tile=>curr%array_of_tiles(ix, iy, iz)
                  nptile=curr_tile%np_tile(1)
                  DO i=nptile, 1, -1! LOOP ON PARTICLES
                    partx=curr_tile%part_x(i)
                    party=curr_tile%part_y(i)
                    partz=curr_tile%part_z(i)
                    partux=curr_tile%part_ux(i)
                    partuy=curr_tile%part_uy(i)
                    partuz=curr_tile%part_uz(i)
                    gaminv=curr_tile%part_gaminv(i)
                    partpid=curr_tile%pid(i, 1:npid)

                    ! Case 1: if particle did not leave tile nothing to do
                    IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT.          &
                    curr_tile%x_tile_max)) .AND. ((party .GE. curr_tile%y_tile_min)   &
                    .AND. (party .LT. curr_tile%y_tile_max)) .AND. ((partz .GE.       &
                    curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max)))   &
                    CYCLE

                    ! Case 2: if particle left MPI domain nothing to do now
                    IF ((partx .LT. x_min_local_part) .OR. (partx .GE.                &
                    x_max_local_part)) CYCLE
                    IF ((party .LT. y_min_local_part) .OR. (party .GE.                &
                    y_max_local_part)) CYCLE
                    IF ((partz .LT. z_min_local_part) .OR. (partz .GE.                &
                    z_max_local_part)) CYCLE

                    ! Case 3: particles changed tile. Tranfer particle to new tile
                    ! Get new indexes of particle in array of tiles
                    indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx), &
                    idp)+1, ntilex)
                    indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy), &
                    idp)+1, ntiley)
                    indz =                                                            &
                    MIN(FLOOR((partz-(z_min_local)+dz/2_num)/(nz0_grid_tile*dz),      &
                    idp)+1, ntilez)
                    CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                    CALL add_particle_at_tile(curr, indx, indy, indz, partx, party,   &
                    partz, partux, partuy, partuz, gaminv, partpid)
                  END DO!END LOOP ON PARTICLES
                END DO
              END DO
            END DO! END LOOP ON TILES
            !$OMP END PARALLEL DO
          END DO
        END DO
      END DO
    END DO! END LOOP ON SPECIES
    !$OMP END PARALLEL DO
    DEALLOCATE(partpid)
  END SUBROUTINE particle_bcs_tiles_openmp

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition on tiles - 2D version without OpenMP
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> 2016
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs_tiles_2d
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indz
    INTEGER(idp) :: nptile, nx0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(:), ALLOCATABLE :: partpid

    ALLOCATE(partpid(npid))
    iy=1
    DO ispecies=1, nspecies! LOOP ON SPECIES
      curr=> species_parray(ispecies)
      ! Get first tiles dimensions (may be different from last tile)
      nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
      nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile
      DO iz=1, ntilez! LOOP ON TILES
        DO ix=1, ntilex
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          nptile=curr_tile%np_tile(1)
          DO i=nptile, 1, -1! LOOP ON PARTICLES
            partx=curr_tile%part_x(i)
            party=curr_tile%part_y(i)
            partz=curr_tile%part_z(i)
            partux=curr_tile%part_ux(i)
            partuy=curr_tile%part_uy(i)
            partuz=curr_tile%part_uz(i)
            gaminv=curr_tile%part_gaminv(i)
            partpid=curr_tile%pid(i, 1:npid)

            ! Case 1: if particle did not leave tile nothing to do
            IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT.                  &
            curr_tile%x_tile_max)) .AND. ((partz .GE. curr_tile%z_tile_min) .AND.     &
            (partz .LT. curr_tile%z_tile_max))) CYCLE

            ! Case 2: if particle left MPI domain nothing to do now
            IF ((partx .LT. x_min_local_part) .OR. (partx .GE. x_max_local_part))     &
            CYCLE
            IF ((partz .LT. z_min_local_part) .OR. (partz .GE. z_max_local_part))     &
            CYCLE

            ! Case 3: particles changed tile. Tranfer particle to new tile
            ! Get new indexes of particle in array of tiles
            indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx), idp)+1, &
            ntilex)
            indz = MIN(FLOOR((partz-(z_min_local)+dz/2_num)/(nz0_grid_tile*dz),       &
            idp)+1, ntilez)
            CALL rm_particle_at_tile(curr, ix, iy, iz, i)
            CALL add_particle_at_tile(curr, indx, iy, indz, partx, party, partz,      &
            partux, partuy, partuz, gaminv, partpid)
          END DO!END LOOP ON PARTICLES
        END DO
      END DO! END LOOP ON TILES
    END DO! END LOOP ON SPECIES
    DEALLOCATE(partpid)
  END SUBROUTINE particle_bcs_tiles_2d

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition on tiles - 2D version with OpenMP
  !
  !> @details
  !> This version is efficient when the number of tiles is large
  !> compared to the number of threads
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> 2016
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs_tiles_2d_openmp
#ifdef _OPENMP
    USE omp_lib
#endif
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz, ipx, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(:), ALLOCATABLE :: partpid
    INTEGER(idp) :: nthreads_tot, nthreads_loop1, nthreads_loop2

#ifdef _OPENMP
    nthreads_tot=OMP_GET_MAX_THREADS()
    CALL OMP_SET_NESTED(.TRUE.)
#else
    nthreads_tot=1
#endif

    IF (nthreads_tot .GT. 1) THEN
      nthreads_loop1=MIN(nspecies, nthreads_tot)
      nthreads_loop2=MAX(1_idp, nthreads_tot/nthreads_loop1)
    ELSE
      nthreads_loop1=1
      nthreads_loop2=1
    ENDIF
    ALLOCATE(partpid(npid))
    iy=1
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr, ispecies, nx0_grid_tile,            &
    !$OMP ny0_grid_tile, nz0_grid_tile, ipx, ipz, partx, party, partz, partux,        &
    !$OMP partuy, partuz, gaminv, partpid, indx, indy, indz, curr_tile, nptile)       &
    !$OMP SHARED(iy, nspecies, npid, nthreads_loop2, species_parray, ntilex, ntiley,  &
    !$OMP ntilez, x_min_local, y_min_local, z_min_local, x_min_local_part,            &
    !$OMP y_min_local_part, z_min_local_part, x_max_local_part, y_max_local_part,     &
    !$OMP z_max_local_part, dx, dy, dz) NUM_THREADS(nthreads_loop1)
    DO ispecies=1, nspecies! LOOP ON SPECIES
      curr=> species_parray(ispecies)
      ! Get first tiles dimensions (may be different from last tile)
      nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
      ny0_grid_tile = curr%array_of_tiles(1, 1, 1)%ny_grid_tile
      nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile
      DO ipz=1, 3
        DO ipx=1, 3
          !$OMP PARALLEL DO DEFAULT(NONE) SHARED(iy, curr, npid, ntilex, ntiley,      &
          !$OMP ntilez, x_min_local_part, y_min_local_part, z_min_local_part,         &
          !$OMP x_max_local_part, y_max_local_part, z_max_local_part, x_min_local,    &
          !$OMP y_min_local, z_min_local, dx, dy, dz, nx0_grid_tile, ny0_grid_tile,   &
          !$OMP nz0_grid_tile) FIRSTPRIVATE(ipx, ipz) PRIVATE(ix, iz, i, curr_tile,   &
          !$OMP nptile, partx, party, partz, partux, partuy, partuz, gaminv, partpid, &
          !$OMP indx, indy, indz) COLLAPSE(2) SCHEDULE(runtime)                       &
          !$OMP NUM_THREADS(nthreads_loop2)
          DO iz=ipz, ntilez, 3! LOOP ON TILES
            DO ix=ipx, ntilex, 3
              curr_tile=>curr%array_of_tiles(ix, iy, iz)
              nptile=curr_tile%np_tile(1)
              DO i=nptile, 1, -1! LOOP ON PARTICLES
                partx=curr_tile%part_x(i)
                !party=curr_tile%part_y(i)
                partz=curr_tile%part_z(i)
                partux=curr_tile%part_ux(i)
                partuy=curr_tile%part_uy(i)
                partuz=curr_tile%part_uz(i)
                gaminv=curr_tile%part_gaminv(i)
                partpid=curr_tile%pid(i, 1:npid)

                ! Case 1: if particle did not leave tile nothing to do
                IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT.              &
                curr_tile%x_tile_max)) .AND. ((partz .GE. curr_tile%z_tile_min) .AND. &
                (partz .LT. curr_tile%z_tile_max))) CYCLE

                ! Case 2: if particle left MPI domain nothing to do now
                IF ((partx .LT. x_min_local_part) .OR. (partx .GE. x_max_local_part)) &
                CYCLE
                IF ((partz .LT. z_min_local_part) .OR. (partz .GE. z_max_local_part)) &
                CYCLE

                ! Case 3: particles changed tile. Tranfer particle to new tile
                ! Get new indexes of particle in array of tiles
                indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),     &
                idp)+1, ntilex)
                indz = MIN(FLOOR((partz-(z_min_local)+dz/2_num)/(nz0_grid_tile*dz),   &
                idp)+1, ntilez)
                !CALL rm_particle_at_tile_2d(curr, ix, iy, iz, i)
                CALL rm_particle_at_tile_2d(curr, ix, iz, i)
                !CALL add_particle_at_tile(curr, indx, iy, indz, &
                !   partx, party, partz, partux, partuy, partuz, gaminv, partpid)
                CALL add_particle_at_tile_2d(curr, indx, indz, partx, partz, partux,  &
                partuy, partuz, gaminv, partpid)
              END DO!END LOOP ON PARTICLES
            END DO
          END DO! END LOOP ON TILES
          !$OMP END PARALLEL DO
        END DO
      END DO
    END DO! END LOOP ON SPECIES
    !$OMP END PARALLEL DO
    DEALLOCATE(partpid)
  END SUBROUTINE particle_bcs_tiles_2d_openmp

  ! Experimental subroutine
#if defined(DEV)
  ! ______________________________________________________________________________________
  ! ********** EXPERIMENTAL subroutine ****************************
  ! => should not be used
  ! This subroutine process the particle boundary conditions
  ! For this aim, the particles are reordered into buckets in the particle arrays
  ! Each bucket corresponds to the group of particle to be exchanged in a common direction
  ! Exchanges between tiles and between MPI domains is combined for better efficiency
  ! This subroutine is less efficient that particle_bcs_tiles_openmp
  ! Mathieu Lobet, 2016
  ! ______________________________________________________________________________________
  SUBROUTINE particle_bsc_openmp_reordering
USE communications, ONLY: part_com_buffer
#ifdef _OPENMP
USE omp_lib
#endif
USE picsar_precision, ONLY: idp, num

    IMPLICIT NONE

    INTEGER(idp) :: i, is, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER    :: curr_tile, curr_tile_add
    REAL(num)    :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(:), ALLOCATABLE :: partpid
    INTEGER(idp) :: test =0, nthreads_tot, nthreads_loop1, nthreads_loop2
    INTEGER(idp) :: dirx, diry, dirz
    INTEGER(idp) :: k, ib
    INTEGER(idp) :: ipmin, ipmax
    type(part_com_buffer), ALLOCATABLE, DIMENSION(:, :, :, :) :: buffer


#ifdef _OPENMP
    nthreads_tot=OMP_GET_MAX_THREADS()
    CALL OMP_SET_NESTED(.TRUE.)
#else
    nthreads_tot=1
#endif


    IF (nthreads_tot .GT. 1) THEN
      nthreads_loop1=MIN(nspecies, nthreads_tot)
      nthreads_loop2=MAX(1_idp, nthreads_tot/nthreads_loop1)
    ELSE
      nthreads_loop1=1
      nthreads_loop2=1
    ENDIF

    ALLOCATE(buffer(ntilex, ntiley, ntilez, nspecies))
    ALLOCATE(partpid(npid))

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr, is, nx0_grid_tile, ny0_grid_tile,   &
    !$OMP nz0_grid_tile, ipx, ipy, ipz, indx, indy, indz, partpid, ib, k, dirx, diry, &
    !$OMP dirz, partx, party, partz, curr_tile, nptile, partux, partuy, partuz)       &
    !$OMP SHARED(nspecies, npid, nthreads_loop2, species_parray, ntilex, ntiley,      &
    !$OMP ntilez, x_min_local, y_min_local, z_min_local, x_min_local_part,            &
    !$OMP y_min_local_part, z_min_local_part, x_max_local_part, y_max_local_part,     &
    !$OMP z_max_local_part, dx, dy, dz, buffer) NUM_THREADS(nthreads_loop1)
    DO is=1, nspecies! LOOP ON SPECIES
      curr=> species_parray(is)
      ! Get first tiles dimensions (may be different from last tile)
      nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
      ny0_grid_tile = curr%array_of_tiles(1, 1, 1)%ny_grid_tile
      nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile
      DO ipz=1, 3
        DO ipy=1, 3
          DO ipx=1, 3
            !$OMP PARALLEL DO DEFAULT(NONE) SHARED(curr, npid, ntilex, ntiley,        &
            !$OMP ntilez, x_min_local_part, y_min_local_part, z_min_local_part,       &
            !$OMP x_max_local_part, y_max_local_part, z_max_local_part, x_min_local,  &
            !$OMP y_min_local, z_min_local, dx, dy, dz, buffer, nx0_grid_tile,        &
            !$OMP ny0_grid_tile, nz0_grid_tile) FIRSTPRIVATE(ipx, ipy, ipz, is)       &
            !$OMP PRIVATE(ix, iy, iz, i, curr_tile, nptile, partx, party, partz,      &
            !$OMP partux, partuy, partuz, gaminv, partpid, indx, indy, indz, dirx,    &
            !$OMP diry, dirz, ib, k) COLLAPSE(3) SCHEDULE(runtime)                    &
            !$OMP NUM_THREADS(nthreads_loop2)
            DO iz=ipz, ntilez, 3! LOOP ON TILES
              DO iy=ipy, ntiley, 3
                DO ix=ipx, ntilex, 3
                  curr_tile=>curr%array_of_tiles(ix, iy, iz)
                  nptile=curr_tile%np_tile(1)

                  ! Temporary array
                  ALLOCATE(buffer(ix, iy, iz, is)%part_x(nptile), buffer(ix, iy, iz,  &
                  is)%part_y(nptile), buffer(ix, iy, iz, is)%part_z(nptile),          &
                  buffer(ix, iy, iz, is)%part_ux(nptile), buffer(ix, iy, iz,          &
                  is)%part_uy(nptile), buffer(ix, iy, iz, is)%part_uz(nptile),        &
                  buffer(ix, iy, iz, is)%part_gaminv(nptile), buffer(ix, iy, iz,      &
                  is)%pid(nptile, 1), buffer(ix, iy, iz, is)%boundid(nptile),         &
                  buffer(ix, iy, iz, is)%bin_npart(0:27), buffer(ix, iy, iz,          &
                  is)%bin_pos(0:27))

                  buffer(ix, iy, iz, is)%bin_npart(:) = 0

                  ! Particle reordering
                  ! Step 1 - determine the number of particles in each bin
                  DO i=1, nptile
                    partx=curr_tile%part_x(i)
                    party=curr_tile%part_y(i)
                    partz=curr_tile%part_z(i)
                    partux=curr_tile%part_ux(i)
                    partuy=curr_tile%part_uy(i)
                    partuz=curr_tile%part_uz(i)
                    gaminv=curr_tile%part_gaminv(i)
                    partpid=curr_tile%pid(i, 1:npid)


                    ! Case 1: if particle did not leave tile nothing to do
                    IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT.          &
                    curr_tile%x_tile_max)) .AND. ((party .GE. curr_tile%y_tile_min)   &
                    .AND. (party .LT. curr_tile%y_tile_max)) .AND. ((partz .GE.       &
                    curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max)))   &
                    THEN

                    buffer(ix, iy, iz, is)%bin_npart(0) = buffer(ix, iy, iz,          &
                    is)%bin_npart(0) + 1

                    buffer(ix, iy, iz, is)%boundid(i) = 0

                    CYCLE

                    ! Case 2: if particle left MPI domain nothing to do now
                  ELSE IF ((partx .LT. x_min_local_part) .OR. (partx .GE.             &
                    x_max_local_part).AND. (party .LT. y_min_local_part) .OR. (party    &
                    .GE. y_max_local_part).AND. (partz .LT. z_min_local_part) .OR.      &
                    (partz .GE. z_max_local_part)) THEN

                    buffer(ix, iy, iz, is)%bin_npart(0) = buffer(ix, iy, iz,          &
                    is)%bin_npart(0) + 1

                    buffer(ix, iy, iz, is)%boundid(i) = 0

                    CYCLE

                  ENDIF

                  indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),   &
                  idp)+1, ntilex)
                  indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy),   &
                  idp)+1, ntiley)
                  indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz),   &
                  idp)+1, ntilez)

                  ! Direction of the particle
                  dirx = indx - ix
                  diry = indy - iy
                  dirz = indz - iz

                  buffer(ix, iy, iz, is)%boundid(i) = (2+dirx) + (1+diry)*3 +         &
                  (1+dirz)*9

                  !if (buffer(ix, iy, iz, is)%boundid(i)>27) then
                  !print*, indx, indy, indz
                  !print*, dirx, diry, dirz
                  !end if

                  buffer(ix, iy, iz, is)%bin_npart(buffer(ix, iy, iz, is)%boundid(i)) &
                  = buffer(ix, iy, iz, is)%bin_npart(buffer(ix, iy, iz,               &
                  is)%boundid(i)) + 1

                ENDDO

                ! Particle reordering
                ! Step 2 - Determine bin positions
                buffer(ix, iy, iz, is)%bin_pos(0) = 1
                Do i=1, 27


                  buffer(ix, iy, iz, is)%bin_pos(i) = buffer(ix, iy, iz,              &
                  is)%bin_pos(i-1) + buffer(ix, iy, iz, is)%bin_npart(i-1)
                ENDDO

                ! Particle reordering
                ! Step 3 - reorder the particles
                DO i=1, nptile
                  ib = buffer(ix, iy, iz, is)%boundid(i)
                  k = buffer(ix, iy, iz, is)%bin_pos(ib)

                  buffer(ix, iy, iz, is)%part_x(k) = curr_tile%part_x(i)
                  buffer(ix, iy, iz, is)%part_y(k) = curr_tile%part_y(i)
                  buffer(ix, iy, iz, is)%part_z(k) = curr_tile%part_z(i)
                  buffer(ix, iy, iz, is)%part_ux(k) = curr_tile%part_ux(i)
                  buffer(ix, iy, iz, is)%part_uy(k) = curr_tile%part_uy(i)
                  buffer(ix, iy, iz, is)%part_uz(k) = curr_tile%part_uz(i)
                  buffer(ix, iy, iz, is)%part_gaminv(k) = curr_tile%part_gaminv(i)
                  buffer(ix, iy, iz, is)%pid(k, 1:npid) = curr_tile%pid(i, 1:npid)

                  buffer(ix, iy, iz, is)%bin_pos(ib) = buffer(ix, iy, iz,             &
                  is)%bin_pos(ib) + 1
                ENDDO

                buffer(ix, iy, iz, is)%bin_pos(0) = 1
                Do i=1, 27
                  buffer(ix, iy, iz, is)%bin_pos(i) = buffer(ix, iy, iz,              &
                  is)%bin_pos(i-1) + buffer(ix, iy, iz, is)%bin_npart(i-1)
                ENDDO

                ! Particle that will stay in the domain
                curr_tile%part_x(1:buffer(ix, iy, iz, is)%bin_npart(0)) = buffer(ix,  &
                iy, iz, is)%part_x(1:buffer(ix, iy, iz, is)%bin_npart(0))
                curr_tile%part_y(1:buffer(ix, iy, iz, is)%bin_npart(0)) = buffer(ix,  &
                iy, iz, is)%part_y(1:buffer(ix, iy, iz, is)%bin_npart(0))
                curr_tile%part_z(1:buffer(ix, iy, iz, is)%bin_npart(0)) = buffer(ix,  &
                iy, iz, is)%part_z(1:buffer(ix, iy, iz, is)%bin_npart(0))
                curr_tile%part_ux(1:buffer(ix, iy, iz, is)%bin_npart(0)) = buffer(ix, &
                iy, iz, is)%part_ux(1:buffer(ix, iy, iz, is)%bin_npart(0))
                curr_tile%part_uy(1:buffer(ix, iy, iz, is)%bin_npart(0)) = buffer(ix, &
                iy, iz, is)%part_uy(1:buffer(ix, iy, iz, is)%bin_npart(0))
                curr_tile%part_uz(1:buffer(ix, iy, iz, is)%bin_npart(0)) = buffer(ix, &
                iy, iz, is)%part_uz(1:buffer(ix, iy, iz, is)%bin_npart(0))
                curr_tile%part_gaminv(1:buffer(ix, iy, iz, is)%bin_npart(0)) =        &
                buffer(ix, iy, iz, is)%part_gaminv(1:buffer(ix, iy, iz,               &
                is)%bin_npart(0))
                curr_tile%pid(1:buffer(ix, iy, iz, is)%bin_npart(0), 1:npid) =        &
                buffer(ix, iy, iz, is)%pid(1:buffer(ix, iy, iz, is)%bin_npart(0),     &
                1:npid)

                curr_tile%np_tile(1) = buffer(ix, iy, iz, is)%bin_npart(0)

              END DO
            END DO
          END DO! END LOOP ON TILES
          !$OMP END PARALLEL DO
        END DO
      END DO
    END DO
  END DO! END LOOP ON SPECIES
  !$OMP END PARALLEL DO


  ! ______________________________________________________________________________________
  ! Exchange between tiles
  ! Reordered buffer are copied in parallel

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr, is, nx0_grid_tile, ny0_grid_tile,     &
  !$OMP nz0_grid_tile, ipx, ipy, ipz, ib, indx, indy, indz, nptile, ipmin, ipmax,     &
  !$OMP curr_tile) SHARED(nspecies, npid, nthreads_loop2, species_parray, ntilex,     &
  !$OMP ntiley, ntilez, x_min_local, y_min_local, z_min_local, x_max_local,           &
  !$OMP y_max_local, z_max_local, dx, dy, dz, buffer) NUM_THREADS(nthreads_loop1)
  DO is=1, nspecies! LOOP ON SPECIES
    curr=> species_parray(is)
    ! Get first tiles dimensions (may be different from last tile)
    nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
    ny0_grid_tile = curr%array_of_tiles(1, 1, 1)%ny_grid_tile
    nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile
    DO ipz=1, 3
      DO ipy=1, 3
        DO ipx=1, 3
          !$OMP PARALLEL DO DEFAULT(NONE) SHARED(curr, npid, ntilex, ntiley, ntilez,  &
          !$OMP x_min_local, y_min_local, z_min_local, x_max_local, y_max_local,      &
          !$OMP z_max_local, dx, dy, dz, buffer, nx0_grid_tile, ny0_grid_tile,        &
          !$OMP nz0_grid_tile) FIRSTPRIVATE(ipx, ipy, ipz, is) PRIVATE(ix, iy, iz, i, &
          !$OMP curr_tile, nptile, partx, party, partz, partux, partuy, partuz,       &
          !$OMP gaminv, indx, indy, indz, ipmin, ipmax, ib) COLLAPSE(3)               &
          !$OMP SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
          DO iz=ipz, ntilez, 3! LOOP ON TILES
            DO iy=ipy, ntiley, 3
              DO ix=ipx, ntilex, 3
                curr_tile=>curr%array_of_tiles(ix, iy, iz)
                nptile=curr_tile%np_tile(1)

                ! Copy of the buffers in every directions
                DO dirz = -1, 1
                  DO diry = -1, 1
                    DO dirx = -1, 1

                      indx = ix+dirx
                      indy = iy+diry
                      indz = iz+dirz

                      IF((indx>=1).and.(indx<=ntilex).and.                            &
                      (indy>=1).and.(indy<=ntiley).and. (indz>=1).and.(indz<=ntilez)) &
                      THEN

                      ib = 2+dirx + (1+diry)*3 + (1+dirz)*9

                      ipmin = buffer(ix, iy, iz, is)%bin_pos(ib)
                      ipmax = ipmin+buffer(ix, iy, iz, is)%bin_npart(ib)-1

                      CALL add_group_of_particles_at_tile(curr, indx, indy, indz,     &
                      buffer(ix, iy, iz, is)%bin_npart(ib), npid, buffer(ix, iy, iz,  &
                      is)%part_x(ipmin:ipmax), buffer(ix, iy, iz,                     &
                      is)%part_y(ipmin:ipmax), buffer(ix, iy, iz,                     &
                      is)%part_z(ipmin:ipmax), buffer(ix, iy, iz,                     &
                      is)%part_ux(ipmin:ipmax), buffer(ix, iy, iz,                    &
                      is)%part_uy(ipmin:ipmax), buffer(ix, iy, iz,                    &
                      is)%part_uz(ipmin:ipmax), buffer(ix, iy, iz,                    &
                      is)%part_gaminv(ipmin:ipmax), buffer(ix, iy, iz,                &
                      is)%pid(ipmin:ipmax, 1:npid))

                    ENDIF

                  ENDDO
                ENDDO
              ENDDO

            END DO
          END DO
        END DO! END LOOP ON TILES
        !$OMP END PARALLEL DO
      END DO
    END DO
  END DO
END DO! END LOOP ON SPECIES
!$OMP END PARALLEL DO
DEALLOCATE(partpid)
END SUBROUTINE particle_bsc_openmp_reordering
#endif

! ________________________________________________________________________________________
!> @brief
!> MPI Boundary condition routine on particles
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!
! ________________________________________________________________________________________
SUBROUTINE particle_bcs_mpi_blocking
USE buff_exchange_part, ONLY: buff_part
USE mpi
USE picsar_precision, ONLY: idp, isp, lp, num
INTEGER(isp) :: nvar! Simple implementation
INTEGER(idp) :: nold, nnew
INTEGER(isp), DIMENSION(-1:1, -1:1, -1:1) :: nptoexch
INTEGER(isp), DIMENSION(:, :, :), ALLOCATABLE :: out_dir
TYPE(buff_part), ALLOCATABLE, DIMENSION(:, :, :)    :: sendbuf
REAL(num), ALLOCATABLE, DIMENSION(:) :: recvbuf
LOGICAL(lp), ALLOCATABLE, DIMENSION(:) :: mask
INTEGER(isp) :: ibuff, nout, nbuff
INTEGER(isp) :: xbd, ybd, zbd
INTEGER(isp) :: ixp, iyp, izp
INTEGER(isp) :: nsend_buf, nrecv_buf
INTEGER(isp) :: dest, src
LOGICAL(lp)  :: out_of_bounds
INTEGER(idp) :: ispecies, i, ix, iy, iz
INTEGER(idp) :: ixtile, iytile, iztile
REAL(num) :: part_xyz, tdeb, tend
TYPE(particle_species), POINTER :: currsp
TYPE(particle_tile), POINTER :: curr

nvar=npid+7

! - LOOP ON SPECIES
DO ispecies=1, nspecies!
  tdeb=MPI_WTIME()
  ! Init send recv buffers
  currsp => species_parray(ispecies)
  nptoexch=0_isp
  nsend_buf=0_isp
  nout=0_isp
  nrecv_buf=0
  nbuff=0_isp
  ibuff=1_isp
  ! ALLOCATE and init sendbuf attributes for particle exchanges
  ALLOCATE(sendbuf(-1:1,-1:1,-1:1))
  DO iz=-1,1
    DO iy=-1,1
      DO ix=-1,1
        sendbuf(ix,iy,iz)%ibuff=1_idp
        sendbuf(ix,iy,iz)%nbuff=1_idp
        ALLOCATE(sendbuf(ix,iy,iz)%buff_arr(sendbuf(ix,iy,iz)%nbuff))
      ENDDO
    ENDDO
  ENDDO 
  ! - LOOP ON TILES 
  DO iztile=1, ntilez!LOOP ON TILES
    DO iytile=1, ntiley
      DO ixtile=1, ntilex
        curr=>currsp%array_of_tiles(ixtile, iytile, iztile)
        ! If not subdomain border, nothing to do
        IF (.NOT. curr%subdomain_bound) CYCLE
        ! Else, search for outbound particles
        ALLOCATE(mask(1:curr%np_tile(1)))
        mask=.TRUE.
        part_xyz=0.
        ! Identify outbounds particles
        DO i = 1, curr%np_tile(1)!LOOP ON PARTICLES
          xbd = 0
          ybd = 0
          zbd = 0
          out_of_bounds = .FALSE.
          part_xyz = curr%part_x(i)
          ! Particle has left this processor -x
          IF (part_xyz .LT. x_min_local_part) THEN
            xbd = -1
            IF (x_min_boundary_part) THEN
              SELECT CASE (pbound_x_min)
              CASE (1_idp)! absorbing
                mask(i)=.FALSE.
                CYCLE
              CASE (2_idp)! Reflecting
                curr%part_x(i) = part_xyz + dx
                curr%part_ux(i) = - curr%part_ux(i)
                xbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_x(i)=2._num*xmin_part-curr%part_x(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_x(i) = part_xyz + length_x_part
              END SELECT
            ENDIF
          ENDIF

          ! Particle has left this processor +x
          IF (part_xyz .GE. x_max_local_part) THEN
            xbd = 1
            IF (x_max_boundary_part) THEN
              SELECT CASE (pbound_x_max)
              CASE (1_idp)! absorbing
                mask(i)=.FALSE.
                CYCLE
              CASE (2_idp)! Reflecting
                curr%part_x(i) = part_xyz - dx
                curr%part_ux(i) = - curr%part_ux(i)
                xbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_x(i)=2._num*xmax_part-curr%part_x(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_x(i) = part_xyz - length_x_part
              END SELECT
            ENDIF
          ENDIF

          part_xyz = curr%part_y(i)
          ! Particle has left this processor -y
          IF ((part_xyz .LT. y_min_local_part) .AND. (c_dim .EQ. 3)) THEN
            ybd = -1
            IF (y_min_boundary_part) THEN
              SELECT CASE (pbound_y_min)! absorbing
              CASE (1_idp)
                mask(i)=.FALSE.
                CYCLE
              CASE (2_idp)! Reflecting
                curr%part_y(i) = part_xyz + dy
                curr%part_uy(i) = - curr%part_uy(i)
                ybd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_y(i)=2._num*ymin_part-curr%part_y(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_y(i) = part_xyz + length_y_part
              END SELECT
            ENDIF
          ENDIF

          ! Particle has left this processor
          IF ((part_xyz .GE. y_max_local_part) .AND. (c_dim .EQ. 3)) THEN
            ybd = 1
            IF (y_max_boundary_part) THEN
              SELECT CASE (pbound_y_max)
              CASE (1_idp)! absorbing
                mask(i)=.FALSE.
                CYCLE
              CASE (2_idp)! Reflecting
                curr%part_y(i) = part_xyz - dy
                curr%part_uy(i) = - curr%part_uy(i)
                ybd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_y(i)=2._num*ymax_part-curr%part_y(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_y(i) = part_xyz - length_y_part
              END SELECT
            ENDIF
          ENDIF

          part_xyz = curr%part_z(i)
          ! Particle has left this processor
          IF (part_xyz .LT. z_min_local_part) THEN
            zbd = -1
            IF (z_min_boundary_part) THEN
              SELECT CASE (pbound_z_min)
              CASE (1_idp)! absorbing
                mask(i)=.FALSE.
                CYCLE
              CASE (2_idp)! Reflecting
                curr%part_z(i) = part_xyz + dz
                curr%part_uz(i) = - curr%part_uz(i)
                zbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_z(i)=2._num*zmin_part-curr%part_z(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_z(i) = part_xyz + length_z_part
              END SELECT
            ENDIF
          ENDIF

          ! Particle has left this processor
          IF (part_xyz .GE. z_max_local_part) THEN
            zbd = 1
            ! Particle has left the system
            IF (z_max_boundary_part) THEN
              SELECT CASE (pbound_z_max)
              CASE (1_idp)! absorbing
                mask(i)=.FALSE.
                CYCLE
              CASE (2_idp)! Reflecting
                curr%part_z(i) = part_xyz - dz
                curr%part_uz(i) = - curr%part_uz(i)
                zbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_z(i)=2._num*zmax_part-curr%part_z(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_z(i) = part_xyz - length_z_part
              END SELECT
            ENDIF
          ENDIF

          IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
            ! Particle has left processor, send it to its neighbour
            mask(i)=.FALSE.
            nout=nout+1
            ibuff=sendbuf(xbd, ybd, zbd)%ibuff
            IF (ibuff+nvar .GT. sendbuf(xbd, ybd, zbd)%nbuff) THEN
                nold = sendbuf(xbd, ybd, zbd)%nbuff
                nnew = 2_idp*(sendbuf(xbd, ybd, zbd)%nbuff+nvar+1_idp) ! ARRAY LIST TYPE
				CALL resize_1D_array_real(sendbuf(xbd, ybd, zbd)%buff_arr, nold, nnew)
				sendbuf(xbd, ybd, zbd)%nbuff=nnew
            ENDIF 
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff)    = curr%part_x(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+1)  = curr%part_y(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+2)  = curr%part_z(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+3)  = curr%part_ux(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+4)  = curr%part_uy(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+5)  = curr%part_uz(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+6)  = curr%part_gaminv(i)
            sendbuf(xbd, ybd, zbd)%buff_arr(ibuff+7:ibuff+6+npid)  = curr%pid(i, 1:npid)
            nptoexch(xbd, ybd, zbd) = nptoexch(xbd, ybd, zbd)+1
            sendbuf(xbd, ybd, zbd)%ibuff  = nptoexch(xbd, ybd, zbd)*nvar+1
          ENDIF
        ENDDO!END LOOP ON PARTICLES
        ! Remove outbound particles from current tile
        CALL rm_particles_from_species_with_mask(currsp, ixtile, iytile, iztile,      &
        mask)
        DEALLOCATE(mask)
      ENDDO
    ENDDO
  ENDDO! END LOOP ON TILES
  tend=MPI_WTIME()
  local_time_part=local_time_part+(tend-tdeb)
  ! SEND/RECEIVE PARTICLES TO/FROM ADJACENT SUBDOMAINS
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1
        IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
        ixp = -ix
        iyp = -iy
        izp = -iz
        ! SEND - RECEIVE PARTICLES IN BUFFERS
        !- Get number of particles in recvbuff
        nsend_buf=nptoexch(ix, iy, iz)*nvar
        nrecv_buf=0
        dest = INT(neighbour(ix, iy, iz), isp)
        src  = INT(neighbour(ixp, iyp, izp), isp)
        CALL MPI_SENDRECV(nsend_buf, 1_isp, MPI_INTEGER, dest, tag, nrecv_buf, 1_isp, &
        MPI_INTEGER, src, tag, comm, status, errcode)
        ALLOCATE(recvbuf(1:nrecv_buf))
        CALL MPI_SENDRECV(sendbuf(ix, iy, iz)%buff_arr(1:nsend_buf), nsend_buf,       & 
        mpidbl, dest, tag, recvbuf, nrecv_buf, mpidbl, src, tag, comm, status, errcode)
        ! Add received particles to particle arrays
        DO i =1, nrecv_buf, nvar
          CALL add_particle_to_species(currsp, recvbuf(i), recvbuf(i+1),              &
          recvbuf(i+2), recvbuf(i+3), recvbuf(i+4), recvbuf(i+5), recvbuf(i+6),       &
          recvbuf(i+7:i+6+npid))
        END DO
        DEALLOCATE(recvbuf)
      ENDDO
    ENDDO
  ENDDO
  ! DEALLOCATE sendbuf and its attributes 
  DO iz=-1,1
    DO iy=-1,1
      DO ix=-1,1
        DEALLOCATE(sendbuf(ix,iy,iz)%buff_arr)
      ENDDO
    ENDDO
  ENDDO 
  DEALLOCATE(sendbuf)
END DO! End loop on species
END SUBROUTINE particle_bcs_mpi_blocking

! ________________________________________________________________________________________
!> @brief
!> MPI Boundary condition routine on particles with non blocking mpi
!> communication subroutine
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
!> Revision 11/2016
!>
!> Modifications:
!> Mathieu - 11/09/2016 - add use mpi for compatibility with my version of mpi
! ________________________________________________________________________________________
SUBROUTINE particle_bcs_mpi_non_blocking
USE buff_exchange_part, ONLY: buff_part
USE mpi
USE picsar_precision, ONLY: idp, isp, lp, num
IMPLICIT NONE
INTEGER(isp) :: nvar! Simple implementation
INTEGER(isp), DIMENSION(-1:1, -1:1, -1:1) :: nptoexch
TYPE(buff_part), ALLOCATABLE, DIMENSION(:, :, :)    :: sendbuff, recvbuff
LOGICAL(lp)  :: remove_from_sim
INTEGER(isp) :: ibuff, nbuff
INTEGER(isp) :: xbd, ybd, zbd
INTEGER(isp) :: mpitag, count
INTEGER(idp), ALLOCATABLE, DIMENSION(:, :, :, :) :: npart_recv, npart_send
INTEGER(isp) :: dest, src, ireq
INTEGER(isp), DIMENSION(:), ALLOCATABLE :: requests
INTEGER(idp) :: ispecies, i, ix, iy, iz, npcurr, ipart, nold, nnew
INTEGER(idp) :: ixtile, iytile, iztile, ispec, nmax
REAL(num) :: part_xyz
TYPE(particle_species), POINTER :: currsp
TYPE(particle_tile), POINTER :: curr

mpitag=0_isp
ALLOCATE(npart_send(1:nspecies, -1:1, -1:1, -1:1))
ALLOCATE(npart_recv(1:nspecies, -1:1, -1:1, -1:1))
ALLOCATE(requests(27*2))

npart_recv=0
npart_send=0
nvar=7+npid
! POST IRECV PARTICLES FROM ADJACENT DOMAINS
! ----- POST ISEND FOR THE NUMBER OF PARTICLES
ireq=1
DO iz = -1, 1
  DO iy = -1, 1
    DO ix = -1, 1
      IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
      count=nspecies
      dest = neighbour(ix, iy, iz)
      CALL MPI_IRECV(npart_recv(1:count, ix, iy, iz), count, MPI_INTEGER8, dest,      &
      mpitag, comm, requests(ireq), errcode)
      ireq=ireq+1
    END DO
  END DO
END DO

! ----- ALLOCATE/INIT SEND BUFFER
ALLOCATE(sendbuff(-1:1,-1:1,-1:1))
DO iz=-1,1
  DO iy=-1,1
    DO ix=-1,1
      sendbuff(ix,iy,iz)%ibuff=1_idp
      sendbuff(ix,iy,iz)%nbuff=1_idp
      ALLOCATE(sendbuff(ix,iy,iz)%buff_arr(sendbuff(ix,iy,iz)%nbuff))
    ENDDO
  ENDDO
ENDDO 

! ----- PUT PARTICLES TO BE SENT IN SENDBUFF BUFFER
nptoexch=0
DO ispecies=1, nspecies!LOOP ON SPECIES
  ! Init send recv buffers
  currsp => species_parray(ispecies)
  DO iztile=1, ntilez!LOOP ON TILES
    DO iytile=1, ntiley
      DO ixtile=1, ntilex
        curr=>currsp%array_of_tiles(ixtile, iytile, iztile)
        ! If not subdomain border, nothing to do
        IF (.NOT. curr%subdomain_bound) CYCLE
        ! Else, search for outbound particles
        part_xyz=0.
        ! Identify outbounds particles
        npcurr=curr%np_tile(1)
        DO i = npcurr, 1, -1!LOOP ON PARTICLES
          xbd = 0
          ybd = 0
          zbd = 0
          remove_from_sim = .FALSE.
          part_xyz = curr%part_x(i)
          ! Particle has left this processor
          IF (part_xyz .LT. x_min_local_part) THEN
            xbd = -1
            IF (x_min_boundary_part) THEN
              SELECT CASE (pbound_x_min)
              CASE (1_idp)! absorbing
                remove_from_sim=.TRUE.
              CASE (2_idp)! Reflecting
                curr%part_x(i) = part_xyz + dx
                curr%part_ux(i) = - curr%part_ux(i)
                xbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_x(i)=2*xmin_part-curr%part_x(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_x(i) = part_xyz + length_x_part
              END SELECT
            ENDIF
          ENDIF
          ! Particle has left this processor
          IF (part_xyz .GE. x_max_local_part) THEN
            xbd = 1
            IF (x_max_boundary_part) THEN
              SELECT CASE (pbound_x_max)
              CASE (1_idp)! absorbing
                remove_from_sim=.TRUE.
              CASE (2_idp)! Reflecting
                curr%part_x(i) = part_xyz - dx
                curr%part_ux(i) = - curr%part_ux(i)
                xbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_x(i)=2*xmax_part-curr%part_x(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_x(i) = part_xyz - length_x_part
              END SELECT
            ENDIF
          ENDIF

          part_xyz = curr%part_y(i)
          ! Particle has left this processor
          IF ((part_xyz .LT. y_min_local_part) .AND. (c_dim .EQ. 3)) THEN
            ybd = -1
            IF (y_min_boundary_part) THEN
              SELECT CASE (pbound_y_min)! absorbing
              CASE (1_idp)
                remove_from_sim=.TRUE.
              CASE (2_idp)! Reflecting
                curr%part_y(i) = part_xyz + dy
                curr%part_uy(i) = - curr%part_uy(i)
                ybd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_y(i)=2*ymin_part-curr%part_y(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_y(i) = part_xyz + length_y_part
              END SELECT
            ENDIF
          ENDIF

          ! Particle has left this processor
          IF ((part_xyz .GE. y_max_local_part) .AND. (c_dim .EQ. 3)) THEN
            ybd = 1
            IF (y_max_boundary_part) THEN
              SELECT CASE (pbound_y_max)
              CASE (1_idp)! absorbing
                remove_from_sim=.TRUE.
              CASE (2_idp)! Reflecting
                curr%part_y(i) = part_xyz - dy
                curr%part_uy(i) = - curr%part_uy(i)
                ybd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_y(i)=2*ymax_part-curr%part_y(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_y(i) = part_xyz - length_y_part
              END SELECT
            ENDIF
          ENDIF

          part_xyz = curr%part_z(i)
          ! Particle has left this processor
          IF (part_xyz .LT. z_min_local_part) THEN
            zbd = -1
            IF (z_min_boundary_part) THEN
              SELECT CASE (pbound_z_min)
              CASE (1_idp)! absorbing
                remove_from_sim=.TRUE.
              CASE (2_idp)! Reflecting
                curr%part_z(i) = part_xyz + dz
                curr%part_uz(i) = - curr%part_uz(i)
                zbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_z(i)=2*zmin_part-curr%part_z(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_z(i) = part_xyz + length_z_part
              END SELECT
            ENDIF
          ENDIF

          ! Particle has left this processor
          IF (part_xyz .GE. z_max_local_part) THEN
            zbd = 1
            ! Particle has left the system
            IF (z_max_boundary_part) THEN
              SELECT CASE (pbound_z_max)
              CASE (1_idp)! absorbing
                remove_from_sim=.TRUE.
              CASE (2_idp)! Reflecting
                curr%part_z(i) = part_xyz - dz
                curr%part_uz(i) = - curr%part_uz(i)
                zbd=0
              CASE (3_idp)! Reinjecting (not thermal for now)
                curr%part_z(i)=2*zmax_part-curr%part_z(i)
                curr%part_ux(i)=0.
                curr%part_uy(i)=0.
                curr%part_uz(i)=0.
                ! Sanity check (keep particle in same tile)
                curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
                curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
                curr%part_y(i)=MAX(curr%part_y(i), curr%y_tile_min+0.5_num*dy)
                curr%part_y(i)=MIN(curr%part_y(i), curr%y_tile_max-0.5_num*dy)
                curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
                curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
                ! Don't remove particle
                CYCLE
              CASE DEFAULT! periodic
                curr%part_z(i) = part_xyz - length_z_part
              END SELECT
            ENDIF
          ENDIF

          IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
            ! Particle has left processor, send it to its neighbour
            IF (.NOT. remove_from_sim) THEN
              !ibuff=nptoexch(xbd, ybd, zbd)*nvar+1
              ibuff=sendbuff(xbd, ybd, zbd)%ibuff
              IF (ibuff+nvar .GT. sendbuff(xbd, ybd, zbd)%nbuff) THEN
                nold = sendbuff(xbd, ybd, zbd)%nbuff
                nnew = 2_idp*(sendbuff(xbd, ybd, zbd)%nbuff+nvar+1_idp) ! ARRAY LIST TYPE
				CALL resize_1D_array_real(sendbuff(xbd, ybd, zbd)%buff_arr, nold, nnew)
				sendbuff(xbd, ybd, zbd)%nbuff=nnew
              ENDIF 
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff)    = curr%part_x(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+1)   = curr%part_y(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+2)   = curr%part_z(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+3)  = curr%part_ux(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+4)  = curr%part_uy(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+5)  = curr%part_uz(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+6)  = curr%part_gaminv(i)
              sendbuff(xbd, ybd, zbd)%buff_arr(ibuff+7:ibuff+6+npid) = curr%pid(i, 1:npid)
              npart_send(ispecies, xbd, ybd, zbd)=npart_send(ispecies, xbd, ybd,      &
              zbd)+1
              nptoexch(xbd, ybd, zbd) = nptoexch(xbd, ybd, zbd)+1
              sendbuff(xbd, ybd, zbd)%ibuff  = nptoexch(xbd, ybd, zbd)*nvar+1
              ! Remove particle of current species from current tile
            ENDIF
            CALL rm_particles_from_species(currsp, ixtile, iytile, iztile, i)
          ENDIF
        ENDDO!END LOOP ON PARTICLES
      ENDDO
    ENDDO
  ENDDO! END LOOP ON TILES

ENDDO! END LOOP ON SPECIES

! ----- POST ISEND FOR THE NUMBER OF PARTICLES
DO iz = -1, 1
  DO iy = -1, 1
    DO ix = -1, 1
      IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
      count=nspecies
      src = INT(neighbour(ix, iy, iz), isp)
      CALL MPI_ISEND(npart_send(1:count, ix, iy, iz), count, MPI_INTEGER8, src,       &
      mpitag, comm, requests(ireq), errcode)
      ireq=ireq+1
    END DO
  END DO
END DO

CALL MPI_WAITALL(ireq-1_isp, requests, MPI_STATUSES_IGNORE, errcode)
requests=0_isp
ireq=1

! ----- ALLOCATION OF RECV BUFFER 
ALLOCATE(recvbuff(-1:1,-1:1,-1:1))
DO iz = -1, 1
  DO iy = -1, 1
    DO ix = -1, 1
      recvbuff(ix,iy,iz)%ibuff=1_idp
      recvbuff(ix,iy,iz)%nbuff=SUM(npart_recv(:, ix, iy, iz))*nvar
      ALLOCATE(recvbuff(ix,iy,iz)%buff_arr(recvbuff(ix,iy,iz)%nbuff))
    END DO
  END DO
END DO

! ----- POST IRECV FOR PARTICLE DATA
DO iz = -1, 1
  DO iy = -1, 1
    DO ix = -1, 1
      count=nvar*SUM(npart_recv(:, ix, iy, iz))
      IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
      IF (count .GT. 0) THEN
        src = INT(neighbour(ix, iy, iz), isp)
        CALL MPI_IRECV(recvbuff(ix, iy, iz)%buff_arr(1:count), count,                  &
        MPI_DOUBLE_PRECISION, src, MPI_ANY_TAG, comm, requests(ireq), errcode)
        ireq=ireq+1
      ENDIF
    END DO
  END DO
END DO

! ----- POST ISEND FOR PARTICLE DATA
DO iz = -1, 1
  DO iy = -1, 1
    DO ix = -1, 1
      count=nvar*SUM(npart_send(:, ix, iy, iz))
      IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
      IF (count .GT. 0) THEN
        dest = INT(neighbour(ix, iy, iz), isp)
        CALL MPI_ISEND(sendbuff(ix, iy, iz)%buff_arr(1:count), count,                  & 
        MPI_DOUBLE_PRECISION, dest, mpitag, comm, requests(ireq), errcode)
        ireq=ireq+1
      ENDIF
    END DO
  END DO
END DO

! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA
count=ireq-1
IF (count .GT. 0_isp) THEN
  CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)
ENDIF

! ----- ADD PARTICLES FROM RECV BUFF TO SPECIES ARRAY
DO iz = -1, 1
  DO iy = -1, 1
    DO ix = -1, 1
      IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
      ispec=0
      DO ispecies=1, nspecies
        currsp=> species_parray(ispecies)
        DO ipart=1, nvar*npart_recv(ispecies, ix, iy, iz), nvar
          ibuff=ispec+ipart
          CALL add_particle_to_species(currsp,                  & 
          recvbuff(ix, iy, iz)%buff_arr(ibuff),                 &
          recvbuff(ix, iy, iz)%buff_arr(ibuff+1),               & 
          recvbuff(ix, iy, iz)%buff_arr(ibuff+2),               &
          recvbuff(ix, iy, iz)%buff_arr(ibuff+3),               & 
          recvbuff(ix, iy, iz)%buff_arr(ibuff+4),               &
          recvbuff(ix, iy, iz)%buff_arr(ibuff+5),               & 
          recvbuff(ix, iy, iz)%buff_arr(ibuff+6),               &
          recvbuff(ix, iy, iz)%buff_arr(ibuff+7:ibuff+6+npid))
        END DO
        ispec=ispec+nvar*npart_recv(ispecies, ix, iy, iz)
      END DO
    END DO
  END DO
END DO

! ----- DEALLOCATE SENDBUFF/RECVBUFF
DO iz=-1,1
  DO iy=-1,1
    DO ix=-1,1
      DEALLOCATE(sendbuff(ix,iy,iz)%buff_arr,recvbuff(ix,iy,iz)%buff_arr)
      ENDDO
  ENDDO
ENDDO 
DEALLOCATE(sendbuff,recvbuff)

DEALLOCATE(npart_send, npart_recv, requests)
END SUBROUTINE particle_bcs_mpi_non_blocking

! ________________________________________________________________________________________
!
!> @brief
!> MPI Boundary condition routine for particles in 2D x, z geometry
!
!> @author
!> Mathieu Lobet
!> Henri Vincenti
!
!> @date
!> Creation 2015
!
!> @warning
!> Need to add reflecting boundary conditions to be consistent with
!> other MPI particle exchange routines
! ________________________________________________________________________________________
SUBROUTINE particle_bcs_mpi_non_blocking_2d
USE mpi
INTEGER(isp)                  :: nvar! Simple implementation
INTEGER(isp), DIMENSION(-1:1, -1:1)       :: nptoexch
REAL(num), ALLOCATABLE, DIMENSION(:, :, :) :: sendbuff, recvbuff
INTEGER(isp)                             :: ibuff, nbuff
INTEGER(isp) :: xbd, zbd
INTEGER(isp) :: mpitag, count
INTEGER(idp), ALLOCATABLE, DIMENSION(:, :, :) :: npart_recv, npart_send
INTEGER(isp) :: dest, src, ireq
INTEGER(isp), DIMENSION(:), ALLOCATABLE :: requests
LOGICAL(lp)  :: out_of_bounds
INTEGER(idp) :: ispecies, i, ix, iz, npcurr, ipart
INTEGER(idp) :: ixtile, iztile, ispec, nmax
REAL(num) :: part_xyz
TYPE(particle_species), POINTER :: currsp
TYPE(particle_tile), POINTER :: curr

mpitag=0_isp
ALLOCATE(npart_send(1:nspecies, -1:1, -1:1))
ALLOCATE(npart_recv(1:nspecies, -1:1, -1:1))
ALLOCATE(requests(27*2))

npart_recv=0
npart_send=0
nvar=7+npid

! POST IRECV PARTICLES FROM ADJACENT DOMAINS
! ----- POST ISEND FOR THE NUMBER OF PARTICLES
ireq=1
DO iz = -1, 1
  DO ix = -1, 1
    IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
    count=nspecies
    dest = neighbour(ix, 0, iz)
    CALL MPI_IRECV(npart_recv(1:count, ix, iz), count, MPI_INTEGER8, dest,            &
    MPI_ANY_TAG, comm, requests(ireq), errcode)
    ireq=ireq+1
  END DO
END DO



! GET NUMBER OF PARTICLES AT BORDER OF CURRENT DOMAIN (INIT SEND BUFFER)
nbuff=0
DO ispecies=1, nspecies
  currsp => species_parray(ispecies)
  DO iztile=1, ntilez!LOOP ON TILES
    DO ixtile=1, ntilex
      curr=>currsp%array_of_tiles(ixtile, 1, iztile)
      IF (.NOT. curr%subdomain_bound) CYCLE
      nbuff=nbuff+ curr%np_tile(1)
    END DO
  END DO
END DO


ALLOCATE(sendbuff(1:nbuff*nvar, -1:1, -1:1))
! PUT PARTICLES TO BE SENT IN BUFFER
nptoexch=0
DO ispecies=1, nspecies!LOOP ON SPECIES
  ! Init send recv buffers
  currsp => species_parray(ispecies)
  DO iztile=1, ntilez!LOOP ON TILES
    DO ixtile=1, ntilex
      curr=>currsp%array_of_tiles(ixtile, 1, iztile)
      ! If not subdomain border, nothing to do
      IF (.NOT. curr%subdomain_bound) CYCLE
      ! Else, search for outbound particles
      part_xyz=0.
      ! Identify outbounds particles
      npcurr=curr%np_tile(1)
      DO i = npcurr, 1, -1!LOOP ON PARTICLES
        xbd = 0
        zbd = 0
        out_of_bounds = .FALSE.
        part_xyz = curr%part_x(i)
        ! Particle has left this processor
        IF (part_xyz .LT. x_min_local_part) THEN
          xbd = -1
          IF (x_min_boundary_part) THEN
            SELECT CASE (pbound_x_min)
            CASE (1_idp)! absorbing
              CALL rm_particles_from_species_2d(currsp, ixtile, iztile, i)
              CYCLE
            CASE (3_idp)! Reinjecting (not thermal for now)
              curr%part_x(i)=2.0_num*xmin_part-curr%part_x(i)
              curr%part_ux(i)=0.
              curr%part_uy(i)=0.
              curr%part_uz(i)=0.
              ! Sanity check (keep particle in same tile)
              curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
              curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
              curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
              curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
              ! Don't remove particle
              CYCLE
            CASE DEFAULT! periodic
              curr%part_x(i) = part_xyz + length_x_part
            END SELECT
          ENDIF
        ENDIF
        ! Particle has left this processor
        IF (part_xyz .GE. x_max_local_part) THEN
          xbd = 1
          IF (x_max_boundary_part) THEN
            SELECT CASE (pbound_x_max)
            CASE (1_idp)! absorbing
              CALL rm_particles_from_species_2d(currsp, ixtile, iztile, i)
              CYCLE
            CASE (3_idp)! Reinjecting (not thermal for now)
              curr%part_x(i)= 2.0_num*xmax_part-curr%part_x(i)
              curr%part_ux(i)=0.
              curr%part_uy(i)=0.
              curr%part_uz(i)=0.
              ! Sanity check (keep particle in same tile)
              curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
              curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
              curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
              curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
              ! Don't remove particle
              CYCLE
            CASE DEFAULT! periodic
              curr%part_x(i) = part_xyz - length_x_part
            END SELECT
          ENDIF
        ENDIF

        part_xyz = curr%part_z(i)
        ! Particle has left this processor
        IF (part_xyz .LT. z_min_local_part) THEN
          zbd = -1
          IF (z_min_boundary_part) THEN
            SELECT CASE (pbound_z_min)
            CASE (1_idp)! absorbing
              CALL rm_particles_from_species_2d(currsp, ixtile, iztile, i)
              CYCLE
            CASE (3_idp)! Reinjecting (not thermal for now)
              curr%part_z(i)= 2.0_num*zmin_part-curr%part_z(i)
              curr%part_ux(i)=0.
              curr%part_uy(i)=0.
              curr%part_uz(i)=0.
              ! Sanity check (keep particle in same tile)
              ! Sanity check (keep particle in same tile)
              curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
              curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
              curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
              curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
              ! Don't remove particle
              CYCLE
            CASE DEFAULT! periodic
              curr%part_z(i) = part_xyz + length_z_part
            END SELECT
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_xyz .GE. z_max_local_part) THEN
          zbd = 1
          ! Particle has left the system
          IF (z_max_boundary_part) THEN
            SELECT CASE (pbound_z_max)
            CASE (1_idp)! absorbing
              CALL rm_particles_from_species_2d(currsp, ixtile, iztile, i)
              CYCLE
            CASE (3_idp)! Reinjecting (not thermal for now)
              curr%part_z(i)= 2.0_num*zmax_part-curr%part_z(i)
              curr%part_ux(i)=0.
              curr%part_uy(i)=0.
              curr%part_uz(i)=0.
              ! Sanity check (keep particle in same tile)
              curr%part_x(i)=MAX(curr%part_x(i), curr%x_tile_min+0.5_num*dx)
              curr%part_x(i)=MIN(curr%part_x(i), curr%x_tile_max-0.5_num*dx)
              curr%part_z(i)=MAX(curr%part_z(i), curr%z_tile_min+0.5_num*dz)
              curr%part_z(i)=MIN(curr%part_z(i), curr%z_tile_max-0.5_num*dz)
              ! Don't remove particle
              CYCLE
            CASE DEFAULT! periodic
              curr%part_z(i) = part_xyz - length_z_part
            END SELECT
          ENDIF
        ENDIF

        IF (ABS(xbd) + ABS(zbd) .GT. 0) THEN
          ! Particle has left processor, send it to its neighbour
          ibuff=nptoexch(xbd, zbd)*nvar+1
          sendbuff(ibuff, xbd, zbd)    = curr%part_x(i)
          sendbuff(ibuff+1, xbd, zbd)  = curr%part_z(i)
          sendbuff(ibuff+2, xbd, zbd)  = curr%part_ux(i)
          sendbuff(ibuff+3, xbd, zbd)  = curr%part_uy(i)
          sendbuff(ibuff+4, xbd, zbd)  = curr%part_uz(i)
          sendbuff(ibuff+5, xbd, zbd)  = curr%part_gaminv(i)
          sendbuff(ibuff+6:ibuff+5+npid, xbd, zbd)  = curr%pid(i, 1:npid)
          npart_send(ispecies, xbd, zbd)=npart_send(ispecies, xbd, zbd)+1
          nptoexch(xbd, zbd) = nptoexch(xbd, zbd)+1
          ! Remove particle of current species from current tile
          CALL rm_particles_from_species_2d(currsp, ixtile, iztile, i)
        ENDIF
      ENDDO!END LOOP ON PARTICLES
    ENDDO
  ENDDO! END LOOP ON TILES
ENDDO! END LOOP ON SPECIES


! ----- POST ISEND FOR THE NUMBER OF PARTICLES
DO iz = -1, 1
  DO ix = -1, 1
    IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
    count=nspecies
    src = INT(neighbour(ix, 0, iz), isp)
    CALL MPI_ISEND(npart_send(1:count, ix, iz), count, MPI_INTEGER8, src, mpitag,     &
    comm, requests(ireq), errcode)
    ireq=ireq+1
  END DO
END DO


CALL MPI_WAITALL(ireq-1_isp, requests, MPI_STATUSES_IGNORE, errcode)
requests=0_isp
ireq=1

! ----- POST IRECV FOR PARTICLE DATA
nmax=nvar*MAXVAL(SUM(npart_recv, 1))
ALLOCATE(recvbuff(nmax, -1:1, -1:1))
DO iz = -1, 1
  DO ix = -1, 1
    count=nvar*SUM(npart_recv(:, ix, iz))
    IF (ABS(ix)  + ABS(iz) .EQ. 0) CYCLE
    IF (count .GT. 0) THEN
      src = INT(neighbour(ix, 0, iz), isp)
      CALL MPI_IRECV(recvbuff(1:count, ix, iz), count, MPI_DOUBLE_PRECISION, src,     &
      MPI_ANY_TAG, comm, requests(ireq), errcode)
      ireq=ireq+1
    ENDIF
  END DO
END DO


! ----- POST ISEND FOR PARTICLE DATA
DO iz = -1, 1
  DO ix = -1, 1
    count=nvar*SUM(npart_send(:, ix, iz))
    IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
    IF (count .GT. 0) THEN
      dest = INT(neighbour(ix, 0, iz), isp)
      CALL MPI_ISEND(sendbuff(1:count, ix, iz), count, MPI_DOUBLE_PRECISION, dest,    &
      mpitag, comm, requests(ireq), errcode)
      ireq=ireq+1
    ENDIF
  END DO
END DO


! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA
count=ireq-1
IF (count .GT. 0_isp) THEN
  CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)
ENDIF


! ----- ADD PARTICLES FROM RECV BUFF TO SPECIES ARRAY
DO iz = -1, 1
  DO ix = -1, 1
    IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
    ispec=0
    DO ispecies=1, nspecies
      currsp=> species_parray(ispecies)
      DO ipart=1, nvar*npart_recv(ispecies, ix, iz), nvar
        ibuff=ispec+ipart
        CALL add_particle_to_species_2d(currsp, recvbuff(ibuff, ix, iz),              &
        recvbuff(ibuff+1, ix, iz), recvbuff(ibuff+2, ix, iz), recvbuff(ibuff+3, ix,   &
        iz), recvbuff(ibuff+4, ix, iz), recvbuff(ibuff+5, ix, iz),                    &
        recvbuff(ibuff+6:ibuff+5+npid, ix, iz))
      END DO
      ispec=ispec+nvar*npart_recv(ispecies, ix, iz)
    END DO
  END DO
END DO


DEALLOCATE(sendbuff, recvbuff, npart_send, npart_recv, requests)
END SUBROUTINE particle_bcs_mpi_non_blocking_2d


! ________________________________________________________________________________________
!> @brief
!> This subroutine combined in a single routine the particle communications
!> between tiles and between MPI domains for 3D
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation: May 2016
!
! ________________________________________________________________________________________
SUBROUTINE particle_bcs_tiles_and_mpi_3d
USE communications, ONLY: mpi_tile_buffer
USE mpi
#ifdef _OPENMP
USE omp_lib
#endif
USE params, ONLY: mpi_buf_size, resize_factor
USE picsar_precision, ONLY: idp, isp, num
IMPLICIT NONE
INTEGER(idp)                    :: is, ix, iy, iz
INTEGER(idp)                    :: i, i2, i3
INTEGER(idp)                    :: indx, indy, indz, ipx, ipy, ipz
INTEGER(idp)                    :: xbd, ybd, zbd
INTEGER(idp)                    :: k, j, ib, ibs
INTEGER(isp)                    :: ireq
INTEGER(isp)                    :: dest, src
INTEGER(idp)                    :: nptile, nx0_grid_tile, ny0_grid_tile,              &
nz0_grid_tile
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER    :: curr_tile
REAL(num)                       :: partx, party, partz, partux, partuy, partuz,       &
gaminv
REAL(num), DIMENSION(:), ALLOCATABLE      :: partpid
INTEGER(idp)                                      :: nthreads_tot
INTEGER(idp)                                      :: nthreads_loop1, nthreads_loop2
INTEGER(idp), dimension(:, :), ALLOCATABLE         :: mpi_npart
REAL(num), dimension(:, :, :, :), ALLOCATABLE        :: bufsend
REAL(num), dimension(:, :), ALLOCATABLE            :: recvbuf
TYPE(mpi_tile_buffer), dimension(:, :, :, :), ALLOCATABLE :: tilebuf
INTEGER(isp), DIMENSION(:, :), ALLOCATABLE         :: nrecv_buf
INTEGER(isp), DIMENSION(:), ALLOCATABLE           :: reqs
INTEGER(isp)                                      :: nrecv_buf_tot, npos, typebuffer
INTEGER(idp)                                      :: new_mpi_buf_size,                &
old_mpi_buf_size
REAL(num)                                         :: nx0_grid_tile_dx
REAL(num)                                         :: ny0_grid_tile_dy,                &
nz0_grid_tile_dz
INTEGER(isp)                                      :: stats(2)
INTEGER(idp)                                      :: recvbuf_index(27)
INTEGER(idp)                                      :: lvect
REAL(num)                                         :: dxs2, dys2, dzs2

! _________________________________________________________
! Determine number of threads to be used for nested parallel region

#ifdef _OPENMP
nthreads_tot=OMP_GET_MAX_THREADS()
!nthreads_tot=OMP_GET_NUM_THREADS()
CALL OMP_SET_NESTED(.TRUE.)
#else
nthreads_tot=1
#endif

IF (nthreads_tot .GT. 1) THEN
  nthreads_loop1=MIN(nspecies, nthreads_tot)
  nthreads_loop2=MAX(1_idp, nthreads_tot/nthreads_loop1)
ELSE
  nthreads_loop1=1
  nthreads_loop2=1
ENDIF

lvect = 64

dxs2 = dx*0.5_num
dys2 = dy*0.5_num
dzs2 = dz*0.5_num

! ___________________________________________________________
! Part 1 - Determine the particle to be exchanged with other tiles or with other MPI domains
!
! In the loop on particles, when a particle has moved to another tile, this particle is
! directly transfer to this tile by memory copy
!
! When a particle has moved to another MPI domain, the particle is copied
! into a local buffer. These local buffer enables to treat exchanges between mpi domains
! and tiles in the same loop
!
! This part os the most time consuming for homogeneous plasmas

#if defined(DEBUG) && (DEBUG==3)
write(0, *) "Part 1 - Determine the particle to be exchanged with other tiles or with &
other MPI domains"
#endif

ALLOCATE(mpi_npart(27, nspecies))
ALLOCATE(tilebuf(ntilex, ntiley, ntilez, nspecies))
ALLOCATE(partpid(npid))

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr, is, ib, k, nx0_grid_tile,               &
!$OMP ny0_grid_tile, nz0_grid_tile, ipx, ipy, ipz, nx0_grid_tile_dx,                  &
!$OMP ny0_grid_tile_dy, nz0_grid_tile_dz, xbd, ybd, zbd, gaminv, partpid, indx, indy, &
!$OMP indz, partx, party, partz, curr_tile, nptile, partux, partuy, partuz, j)        &
!$OMP SHARED(npid, nspecies, nthreads_loop2, species_parray, ntilex, ntiley, ntilez,  &
!$OMP x_min_local_part, y_min_local_part, z_min_local_part, length_x_part,            &
!$OMP length_y_part, length_z_part, dxs2, dys2, dzs2, xmin_part, xmax_part,           &
!$OMP ymin_part, ymax_part, zmin_part, zmax_part, x_min_boundary_part,                &
!$OMP x_max_boundary_part, y_min_boundary_part, y_max_boundary_part,                  &
!$OMP z_min_boundary_part, z_max_boundary_part, pbound_x_min, pbound_x_max,           &
!$OMP pbound_y_min, pbound_y_max, pbound_z_min, pbound_z_max, x_min_local,            &
!$OMP y_min_local, z_min_local, x_max_local_part, y_max_local_part, z_max_local_part, &
!$OMP dx, dy, dz, mpi_npart, tilebuf, mpi_buf_size, lvect)                            &
!$OMP NUM_THREADS(nthreads_loop1)
! LOOP ON SPECIES
DO is=1, nspecies

  curr=> species_parray(is)
  ! Get first tiles dimensions (may be different from last tile)
  nx0_grid_tile = curr%array_of_tiles(1, 1, 1)%nx_grid_tile
  ny0_grid_tile = curr%array_of_tiles(1, 1, 1)%ny_grid_tile
  nz0_grid_tile = curr%array_of_tiles(1, 1, 1)%nz_grid_tile

  nx0_grid_tile_dx = 1._num/(nx0_grid_tile*dx)
  ny0_grid_tile_dy = 1._num/(ny0_grid_tile*dy)
  nz0_grid_tile_dz = 1._num/(nz0_grid_tile*dz)

  mpi_npart(:, is) = 0
  ! LOOP ON TILES
  DO ipz=1, 3
    DO ipy=1, 3
      DO ipx=1, 3


        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(npid, curr, ntilex, ntiley, ntilez,    &
        !$OMP x_min_local_part, y_min_local_part, z_min_local_part, x_max_local_part, &
        !$OMP y_max_local_part, z_max_local_part, dx, dy, dz, nx0_grid_tile,          &
        !$OMP ny0_grid_tile, nz0_grid_tile, pbound_x_min, pbound_x_max, pbound_y_min, &
        !$OMP pbound_y_max, pbound_z_min, pbound_z_max, x_min_local, y_min_local,     &
        !$OMP z_min_local, length_x_part, length_y_part, length_z_part, tilebuf,      &
        !$OMP mpi_npart, xmin_part, xmax_part, ymin_part, ymax_part, zmin_part,       &
        !$OMP zmax_part, x_min_boundary_part, x_max_boundary_part,                    &
        !$OMP y_min_boundary_part, y_max_boundary_part, z_min_boundary_part,          &
        !$OMP z_max_boundary_part, nx0_grid_tile_dx, ny0_grid_tile_dy,                &
        !$OMP nz0_grid_tile_dz, dxs2, dys2, dzs2, mpi_buf_size, lvect)                &
        !$OMP FIRSTPRIVATE(ipx, ipy, ipz, is) PRIVATE(ix, iy, iz, i, ib, k,           &
        !$OMP curr_tile, nptile, partx, party, partz, partux, partuy, partuz, gaminv, &
        !$OMP partpid, indx, indy, indz, xbd, ybd, zbd, i2, i3, new_mpi_buf_size,     &
        !$OMP old_mpi_buf_size) COLLAPSE(3) SCHEDULE(runtime)                         &
        !$OMP NUM_THREADS(nthreads_loop2)
        ! LOOP ON TILES
        DO iz=ipz, ntilez, 3
          DO iy=ipy, ntiley, 3
            DO ix=ipx, ntilex, 3
              curr_tile=>curr%array_of_tiles(ix, iy, iz)
              nptile=curr_tile%np_tile(1)

              ! Allocation of the buffer
              IF (curr_tile%subdomain_bound) THEN
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_x(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_y(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_z(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_ux(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_uy(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_uz(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%part_gaminv(mpi_buf_size, 27))
                ALLOCATE(tilebuf(ix, iy, iz, is)%pid(mpi_buf_size, npid, 27))
              ENDIF
              tilebuf(ix, iy, iz, is)%npart(1:27) = 0

#if defined(DEBUG) && (DEBUG==3)
              write(0, '("Loop on particles inside tiles: ", I2, X, I2, X, I2)')ix,   &
              iy, iz
#endif

              ! Loop on particles inside tiles
              DO i=nptile, 1, -1

                !                 DO i2=nptile, 1, -lvect
                !                   DO i3 = MIN(lvect, i2), 1, -1
                !                   i = i2 - i3 + 1

                partx=curr_tile%part_x(i)
                party=curr_tile%part_y(i)
                partz=curr_tile%part_z(i)
                partux=curr_tile%part_ux(i)
                partuy=curr_tile%part_uy(i)
                partuz=curr_tile%part_uz(i)
                gaminv=curr_tile%part_gaminv(i)
                partpid=curr_tile%pid(i, 1:npid)

                ! Case 1: if particle did not leave tile nothing to do
                IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT.              &
                curr_tile%x_tile_max)) .AND. ((party .GE. curr_tile%y_tile_min) .AND. &
                (party .LT. curr_tile%y_tile_max)) .AND. ((partz .GE.                 &
                curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) CYCLE

                ! Case 2: if particle left MPI domain
                IF (((partx .LT. x_min_local_part) .OR. (partx .GE.                   &
                x_max_local_part)) .OR. ((party .LT. y_min_local_part) .OR. (party    &
                .GE. y_max_local_part)) .OR. ((partz .LT. z_min_local_part) .OR.      &
                (partz .GE. z_max_local_part))) THEN

                ! Then we determine in which domain this particle is going
                xbd = 0
                ybd = 0
                zbd = 0

                ! Particle has left this processor -x
                IF (partx .LT. x_min_local_part) THEN
                  xbd = -1
                  IF (x_min_boundary_part) THEN
                    SELECT CASE (pbound_x_min)
                    CASE (1_idp)! absorbing
                      CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                      CYCLE! If particle has been removed CYCLE to next particle
                    CASE (3_idp)! Reinjecting (not thermal for now)
                      curr_tile%part_x(i)=2*xmin_part-curr_tile%part_x(i)
                      curr_tile%part_ux(i)=0.
                      curr_tile%part_uy(i)=0.
                      curr_tile%part_uz(i)=0.
                      ! Sanity check (keep particle in same tile)
                      curr_tile%part_x(i)=MAX(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_min+0.5_num*dx)
                      curr_tile%part_x(i)=MIN(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_max-0.5_num*dx)
                      curr_tile%part_y(i)=MAX(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_min+0.5_num*dy)
                      curr_tile%part_y(i)=MIN(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_max-0.5_num*dy)
                      curr_tile%part_z(i)=MAX(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_min+0.5_num*dz)
                      curr_tile%part_z(i)=MIN(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_max-0.5_num*dz)
                      ! Don't remove particle
                      CYCLE
                    CASE DEFAULT! periodic
                      curr_tile%part_x(i) = partx + length_x_part
                    END SELECT
                  ENDIF

                  ! Particle has left this processor +x
                ELSE IF (partx .GE. x_max_local_part) THEN
                  xbd = 1
                  IF (x_max_boundary_part) THEN
                    SELECT CASE (pbound_x_max)
                    CASE (1_idp)! absorbing
                      CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                      CYCLE! If particle has been removed CYCLE to next particle
                    CASE (3_idp)! Reinjecting (not thermal for now)
                      curr_tile%part_x(i)=2*xmax_part-curr_tile%part_x(i)
                      curr_tile%part_ux(i)=0.
                      curr_tile%part_uy(i)=0.
                      curr_tile%part_uz(i)=0.
                      ! Sanity check (keep particle in same tile)
                      curr_tile%part_x(i)=MAX(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_min+0.5_num*dx)
                      curr_tile%part_x(i)=MIN(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_max-0.5_num*dx)
                      curr_tile%part_y(i)=MAX(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_min+0.5_num*dy)
                      curr_tile%part_y(i)=MIN(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_max-0.5_num*dy)
                      curr_tile%part_z(i)=MAX(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_min+0.5_num*dz)
                      curr_tile%part_z(i)=MIN(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_max-0.5_num*dz)
                      ! Don't remove particle
                      CYCLE
                    CASE DEFAULT! periodic
                      curr_tile%part_x(i) = partx - length_x_part
                    END SELECT
                  ENDIF
                ENDIF

                ! Particle has left this processor -y
                IF ((party .LT. y_min_local_part)) THEN
                  ybd = -1
                  IF (y_min_boundary_part) THEN
                    SELECT CASE (pbound_y_min)! absorbing
                    CASE (1_idp)
                      CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                      CYCLE! If particle has been removed CYCLE to next particle
                    CASE (3_idp)! Reinjecting (not thermal for now)
                      curr_tile%part_y(i)=2*ymin_part-curr_tile%part_y(i)
                      curr_tile%part_ux(i)=0.
                      curr_tile%part_uy(i)=0.
                      curr_tile%part_uz(i)=0.
                      ! Sanity check (keep particle in same tile)
                      curr_tile%part_x(i)=MAX(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_min+0.5_num*dx)
                      curr_tile%part_x(i)=MIN(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_max-0.5_num*dx)
                      curr_tile%part_y(i)=MAX(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_min+0.5_num*dy)
                      curr_tile%part_y(i)=MIN(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_max-0.5_num*dy)
                      curr_tile%part_z(i)=MAX(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_min+0.5_num*dz)
                      curr_tile%part_z(i)=MIN(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_max-0.5_num*dz)
                      ! Don't remove particle
                      CYCLE
                    CASE DEFAULT! periodic
                      curr_tile%part_y(i) = party + length_y_part
                    END SELECT
                  ENDIF
                  ! Particle has left this processor +y
                ELSE IF ((party .GE. y_max_local_part)) THEN
                  ybd = 1
                  IF (y_max_boundary_part) THEN
                    SELECT CASE (pbound_y_max)
                    CASE (1_idp)! absorbing
                      CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                      CYCLE! If particle has been removed CYCLE to next particle
                    CASE (3_idp)! Reinjecting (not thermal for now)
                      curr_tile%part_y(i)=2*ymax_part-curr_tile%part_y(i)
                      curr_tile%part_ux(i)=0.
                      curr_tile%part_uy(i)=0.
                      curr_tile%part_uz(i)=0.
                      ! Sanity check (keep particle in same tile)
                      curr_tile%part_x(i)=MAX(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_min+0.5_num*dx)
                      curr_tile%part_x(i)=MIN(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_max-0.5_num*dx)
                      curr_tile%part_y(i)=MAX(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_min+0.5_num*dy)
                      curr_tile%part_y(i)=MIN(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_max-0.5_num*dy)
                      curr_tile%part_z(i)=MAX(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_min+0.5_num*dz)
                      curr_tile%part_z(i)=MIN(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_max-0.5_num*dz)
                      ! Don't remove particle
                      CYCLE
                    CASE DEFAULT! periodic
                      curr_tile%part_y(i) = party - length_y_part
                    END SELECT
                  ENDIF
                ENDIF

                ! Particle has left this processor -z
                IF (partz .LT. z_min_local_part) THEN
                  zbd = -1
                  IF (z_min_boundary_part) THEN
                    SELECT CASE (pbound_z_min)
                    CASE (1_idp)! absorbing
                      CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                      CYCLE! If particle has been removed CYCLE to next particle
                    CASE (3_idp)! Reinjecting (not thermal for now)
                      curr_tile%part_z(i)=2*zmin_part-curr_tile%part_z(i)
                      curr_tile%part_ux(i)=0.
                      curr_tile%part_uy(i)=0.
                      curr_tile%part_uz(i)=0.
                      ! Sanity check (keep particle in same tile)
                      curr_tile%part_x(i)=MAX(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_min+0.5_num*dx)
                      curr_tile%part_x(i)=MIN(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_max-0.5_num*dx)
                      curr_tile%part_y(i)=MAX(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_min+0.5_num*dy)
                      curr_tile%part_y(i)=MIN(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_max-0.5_num*dy)
                      curr_tile%part_z(i)=MAX(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_min+0.5_num*dz)
                      curr_tile%part_z(i)=MIN(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_max-0.5_num*dz)
                      ! Don't remove particle
                      CYCLE
                    CASE DEFAULT! periodic
                      curr_tile%part_z(i) = partz + length_z_part
                    END SELECT
                  ENDIF
                  ! Particle has left this processor +z
                ELSE IF (partz .GE. z_max_local_part) THEN
                  zbd = 1
                  ! Particle has left the system
                  IF (z_max_boundary_part) THEN
                    SELECT CASE (pbound_z_max)
                    CASE (1_idp)! absorbing
                      CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                      CYCLE! If particle has been removed CYCLE to next particle
                    CASE (3_idp)! Reinjecting (not thermal for now)
                      curr_tile%part_z(i)=2*zmax_part-curr_tile%part_z(i)
                      curr_tile%part_ux(i)=0.
                      curr_tile%part_uy(i)=0.
                      curr_tile%part_uz(i)=0.
                      ! Sanity check (keep particle in same tile)
                      curr_tile%part_x(i)=MAX(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_min+0.5_num*dx)
                      curr_tile%part_x(i)=MIN(curr_tile%part_x(i),                    &
                      curr_tile%x_tile_max-0.5_num*dx)
                      curr_tile%part_y(i)=MAX(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_min+0.5_num*dy)
                      curr_tile%part_y(i)=MIN(curr_tile%part_y(i),                    &
                      curr_tile%y_tile_max-0.5_num*dy)
                      curr_tile%part_z(i)=MAX(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_min+0.5_num*dz)
                      curr_tile%part_z(i)=MIN(curr_tile%part_z(i),                    &
                      curr_tile%z_tile_max-0.5_num*dz)
                      ! Don't remove particle
                      CYCLE
                    CASE DEFAULT! periodic
                      curr_tile%part_z(i) = partz - length_z_part
                    END SELECT
                  ENDIF
                ENDIF


                ! Particle has left processor, we put it in a local buffer
                ! Boundary index 1-27
                ib = 2+xbd + (1+ybd)*3 + (1+zbd)*9
                ! Current particle available index
                k = tilebuf(ix, iy, iz, is)%npart(ib) + 1

                ! Particle properties are placed in a buffer for each com direction
                tilebuf(ix, iy, iz, is)%part_x(k, ib) = curr_tile%part_x(i)
                tilebuf(ix, iy, iz, is)%part_y(k, ib) = curr_tile%part_y(i)
                tilebuf(ix, iy, iz, is)%part_z(k, ib) = curr_tile%part_z(i)
                tilebuf(ix, iy, iz, is)%part_ux(k, ib) = curr_tile%part_ux(i)
                tilebuf(ix, iy, iz, is)%part_uy(k, ib) = curr_tile%part_uy(i)
                tilebuf(ix, iy, iz, is)%part_uz(k, ib) = curr_tile%part_uz(i)
                tilebuf(ix, iy, iz, is)%part_gaminv(k, ib) = curr_tile%part_gaminv(i)
                tilebuf(ix, iy, iz, is)%pid(k, 1:npid, ib) = curr_tile%pid(i, 1:npid)

                IF (k.eq.SIZE(tilebuf(ix, iy, iz, is)%part_x, 1)) THEN
                  old_mpi_buf_size = SIZE(tilebuf(ix, iy, iz, is)%part_x, 1)
                  new_mpi_buf_size = SIZE(tilebuf(ix, iy, iz, is)%part_x, 1)*2
                  WRITE(0, '(" WARNING: Tile buffer array has been resized: nbpart =  &
                  ", I7, " new buffer size = ", I7)') k, new_mpi_buf_size
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_x,           &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_y,           &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_z,           &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_ux,          &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_uy,          &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_uz,          &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_2D_array_real(tilebuf(ix, iy, iz, is)%part_gaminv,      &
                  old_mpi_buf_size, new_mpi_buf_size, 27_idp, 27_idp)
                  CALL resize_3D_array_real(tilebuf(ix, iy, iz, is)%pid,              &
                  old_mpi_buf_size, new_mpi_buf_size, npid, npid, 27_idp, 27_idp)
                ENDIF

                ! Update of the number of particles
                tilebuf(ix, iy, iz, is)%npart(ib) = k

                ! The particle is deleted
                CALL rm_particle_at_tile(curr, ix, iy, iz, i)
                CYCLE

              ENDIF!end if particle left MPI domain


              ! Case 3: particles changed tile. Tranfer particle to new tile
              ! Get new indexes of particle in array of tiles
              indx = MIN(FLOOR((partx-x_min_local+dxs2)*(nx0_grid_tile_dx), idp)+1,   &
              ntilex)
              indy = MIN(FLOOR((party-y_min_local+dys2)*(ny0_grid_tile_dy), idp)+1,   &
              ntiley)
              indz = MIN(FLOOR((partz-(z_min_local)+dzs2)*(nz0_grid_tile_dz), idp)+1, &
              ntilez)
              CALL rm_particle_at_tile(curr, ix, iy, iz, i)


              CALL add_particle_at_tile(curr, indx, indy, indz, partx, party, partz,  &
              partux, partuy, partuz, gaminv, partpid)


              !ENDDO
              !End loop on particles
            END DO

            ! Reduction of the total number of particle to be
            ! exchanged in every direction for MPI
            ! Here we use a critical region since we expect the thread
            ! to not be perfectly synchronized
            ! In this part appears to have negative consequences on the parallelization,             ! this process can be performed in sequencial outside the species loop
            IF (curr_tile%subdomain_bound) THEN
              !OMP CRITICAL
              mpi_npart(1:27, is) = mpi_npart(1:27, is) + tilebuf(ix, iy, iz,         &
              is)%npart(1:27)
              !OMP END CRITICAL

            ENDIF

          END DO
        END DO
      END DO! END LOOP ON TILES
      !$OMP END PARALLEL DO

    END DO
  END DO
END DO
END DO! END LOOP ON SPECIES
!$OMP END PARALLEL DO


#if defined(DEBUG) && (DEBUG==3)
write(0, *) "Part 2 - Creation of the send buffer for the MPI communications"
#endif

! ________________________________________________________________________________________
! Part 2 - Creation of the send buffer for the MPI communications
!
! Here, a global buffer is created for the MPI communications. We copy in parallel
! the local buffers of the tiles into the global buffer using nested parallel regions, ! one for the species and the other on the communication directions
!
! There are 27 directions in 3D but only the 8 faces of the cubic domain represent
! a significant number of particles to exchange when the plasma kinetic
! distribution is uniform
!
! For drifted plasmas, this method is not the most efficient

ALLOCATE(bufsend(MAXVAL(mpi_npart(:, :)), 7+npid, 27, nspecies))

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr, is, nx0_grid_tile, ny0_grid_tile,       &
!$OMP nz0_grid_tile, ipx, ipy, ipz, ib, k, j) SHARED(nspecies, npid, nthreads_loop2,  &
!$OMP species_parray, ntilex, ntiley, ntilez, dx, dy, dz, mpi_npart, bufsend,         &
!$OMP tilebuf) NUM_THREADS(nthreads_loop1)
DO is=1, nspecies! LOOP ON SPECIES

mpi_npart(:, is) = 0

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(npid, is, tilebuf, mpi_npart, bufsend, ntilez, &
!$OMP ntiley, ntilex) PRIVATE(ix, iy, iz, k, xbd, ybd, zbd, j, ib) COLLAPSE(3)        &
!$OMP SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
DO xbd = -1, 1
  DO ybd = -1, 1
    DO zbd = -1, 1

      IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .EQ. 0) CYCLE

      ib = 2+xbd + (1+ybd)*3 + (1+zbd)*9

      DO iz=1, ntilez! LOOP ON TILES
        DO iy=1, ntiley
          DO ix=1, ntilex
            k = tilebuf(ix, iy, iz, is)%npart(ib)
            !print*, 'iz, iy, ix, ib:', iz, iy, ix, ib, 'k', k

            IF (k.eq.0) CYCLE

            j = mpi_npart(ib, is)
            bufsend(j+1:j+k, 1, ib, is) = tilebuf(ix, iy, iz, is)%part_x(1:k, ib)
            bufsend(j+1:j+k, 2, ib, is) = tilebuf(ix, iy, iz, is)%part_y(1:k, ib)
            bufsend(j+1:j+k, 3, ib, is) = tilebuf(ix, iy, iz, is)%part_z(1:k, ib)
            bufsend(j+1:j+k, 4, ib, is) = tilebuf(ix, iy, iz, is)%part_ux(1:k, ib)
            bufsend(j+1:j+k, 5, ib, is) = tilebuf(ix, iy, iz, is)%part_uy(1:k, ib)
            bufsend(j+1:j+k, 6, ib, is) = tilebuf(ix, iy, iz, is)%part_uz(1:k, ib)
            bufsend(j+1:j+k, 7, ib, is) = tilebuf(ix, iy, iz, is)%part_gaminv(1:k,    &
            ib)
            bufsend(j+1:j+k, 8:7+npid, ib, is) = tilebuf(ix, iy, iz, is)%pid(1:k,     &
            1:npid, ib)
            mpi_npart(ib, is) = j + k
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(tilebuf)

#if defined(DEBUG) && (DEBUG==3)
write(0, *) "Part 3 - MPI Communications"
#endif

! _______________________________________
! Part 3 - MPI Communications
! First, we determine the number of particles to receive

ALLOCATE(nrecv_buf(27, nspecies))

nrecv_buf(:, :)=0

! Thread version
IF (.FALSE.) THEN

ALLOCATE(reqs(2))

DO is=1, nspecies! LOOP ON SPECIES
  !curr=> species_parray(is)

  !$OMP PARALLEL DO DEFAULT(NONE) SHARED(is, mpi_npart, comm, neighbour, nrecv_buf)   &
  !$OMP FIRSTPRIVATE(tag, status, stats, reqs, MPI_STATUSES_IGNORE) PRIVATE(ix, iy,   &
  !$OMP iz, k, ipx, ipy, ipz, dest, src, ib, ibs, errcode) COLLAPSE(3)                &
  !$OMP SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)

  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1
        IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
        ! index of the communication direction
        ib = 2+ix + (1+iy)*3 + (1+iz)*9

        ! indexes of the source/destination mpi task
        ipx = -ix
        ipy = -iy
        ipz = -iz

        !ibs = 2+ipx + (1+ipy)*3 + (1+ipz)*9

        dest = INT(neighbour(ix, iy, iz), isp)
        src  = INT(neighbour(ipx, ipy, ipz), isp)

        ! Number of particle
        k = mpi_npart(ib, is)

        ! Exchange
        CALL MPI_Irecv( nrecv_buf(ib, is), 1_isp, MPI_INTEGER, src, INT(ib, isp),     &
        comm, reqs(1), errcode)

        CALL MPI_Isend(k, 1_isp, MPI_INTEGER, dest, INT(ib, isp), comm, reqs(2),      &
        errcode)

        CALL MPI_Waitall(2_isp, reqs, MPI_STATUSES_IGNORE, errcode)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ENDDO

ELSE

! _______________________________________________________
! Sequential version of the previous block

! Waitall is done after the loops
ALLOCATE(reqs(54*nspecies))
DO is=1, nspecies! LOOP ON SPECIES
  ireq = 0
  !curr=> species_parray(is)
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1
        IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
        ! index of the communication direction
        ib = 2+ix + (1+iy)*3 + (1+iz)*9

        ! indexes of the source/destination mpi task
        ipx = -ix
        ipy = -iy
        ipz = -iz

        !ibs = 2+ipx + (1+ipy)*3 + (1+ipz)*9

        dest = INT(neighbour(ix, iy, iz), isp)
        src  = INT(neighbour(ipx, ipy, ipz), isp)

        ! Number of particle
        k = mpi_npart(ib, is)

        ! Exchange
        ireq = ireq + 1
        CALL MPI_Irecv( nrecv_buf(ib, is), 1_isp, MPI_INTEGER, src, INT(ib, isp),     &
        comm, reqs(ireq), errcode)
        ireq = ireq + 1
        CALL MPI_Isend(k, 1_isp, MPI_INTEGER, dest, INT(ib, isp), comm, reqs(ireq),   &
        errcode)

      ENDDO
    ENDDO
  ENDDO
  CALL MPI_Waitall(ireq, reqs, MPI_STATUSES_IGNORE, errcode)
ENDDO

ENDIF
! ________________________________________________________________________________________
! Then we exchange the buffers
! The communication can be performed in sequential or with multiple thread


! LOOP ON SPECIES
DO is=1, nspecies

nrecv_buf_tot = SUM(nrecv_buf(:, is))
curr=> species_parray(is)
ALLOCATE(recvbuf(1:nrecv_buf_tot+1, 7+npid))

! Multithread version
IF (.FALSE.) THEN

  ! ________________________________________________
  ! Multiple thread version
  ! Determine the position of each received buffer in recvbuf
  npos=1
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1
        ib = 2+ix + (1+iy)*3 + (1+iz)*9
        recvbuf_index(ib) = npos
        npos = npos + nrecv_buf(ib, is)
      ENDDO
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(NONE) SHARED(is, mpi_npart, npid, comm, neighbour,        &
  !$OMP nrecv_buf, nrecv_buf_tot, recvbuf_index, MPI_STATUSES_IGNORE, bufsend,        &
  !$OMP mpidbl, recvbuf) FIRSTPRIVATE(tag, status, stats, reqs) PRIVATE(ix, iy, iz,   &
  !$OMP k, j, ipx, ipy, ipz, dest, src, ib, ibs, errcode, typebuffer) COLLAPSE(3)     &
  !$OMP SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1

        IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE

        ib = 2+ix + (1+iy)*3 + (1+iz)*9

        k = mpi_npart(ib, is)
        j = nrecv_buf(ib, is)

        dest = INT(neighbour(ix, iy, iz), isp)
        src  = INT(neighbour(-ix, -iy, -iz), isp)

        CALL MPI_TYPE_VECTOR(7_isp+INT(npid, isp), INT(j, isp), INT(nrecv_buf_tot+1,  &
        isp), MPI_DOUBLE_PRECISION, typebuffer, errcode)
        call MPI_TYPE_COMMIT(typebuffer, errcode)

        ! Exchange
        CALL MPI_Irecv(recvbuf(recvbuf_index(ib), 1), 1_isp, typebuffer, src, INT(ib, &
        isp), comm, reqs(1), errcode)

        CALL MPI_Isend(bufsend(1:k, 1:7+npid, ib, is), INT((7_idp+npid)*k, isp),      &
        mpidbl, dest, INT(ib, isp), comm, reqs(2), errcode)

        CALL MPI_Waitall(2_isp, reqs, MPI_STATUSES_IGNORE, errcode)

        !CALL MPI_SENDRECV(bufsend(1:k, 1:8, ib, is), 8_isp*k, mpidbl, dest, INT(ib, isp), &
        !recvbuf(recvbuf_index(ib), 1), 1_isp, typebuffer, src, INT(ib, isp), comm, status, errcode)

        call MPI_TYPE_FREE(typebuffer, errcode)

        !npos = npos + j
        !print*, ib, npos, j, iz, iy, ix

      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! ________________________________________________

ELSE

  ! Wait all is done inside the loop
  npos=1
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1

        IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE

        ib = 2+ix + (1+iy)*3 + (1+iz)*9

        k = mpi_npart(ib, is)
        j = nrecv_buf(ib, is)

        dest = INT(neighbour(ix, iy, iz), isp)
        src  = INT(neighbour(-ix, -iy, -iz), isp)

        CALL MPI_TYPE_VECTOR(7_isp+INT(npid, isp), INT(j, isp), INT(nrecv_buf_tot+1,  &
        isp), MPI_DOUBLE_PRECISION, typebuffer, errcode)
        call MPI_TYPE_COMMIT(typebuffer, errcode)
        ! Exchange
        CALL MPI_Irecv(recvbuf(npos, 1), 1_isp, typebuffer, src, INT(ib, isp), comm,  &
        reqs(1), errcode)
        CALL MPI_Isend(bufsend(1:k, 1:7+npid, ib, is), INT((7_idp+npid)*k, isp),      &
        mpidbl, dest, INT(ib, isp), comm, reqs(2), errcode)
        CALL MPI_Waitall(2_isp, reqs, MPI_STATUSES_IGNORE, errcode)

        call MPI_TYPE_FREE(typebuffer, errcode)
        npos = npos + j
        !print*, ib, npos, j, iz, iy, ix

      ENDDO
    ENDDO
  ENDDO

  ! ________________________________________________


ENDIF

#if defined(DEBUG) && (DEBUG==3)
write(0, *) " Copy of the buffers in the tile particle arrays"
#endif

! If the buffer is not empty...
IF (nrecv_buf_tot.gt.0) THEN

  ! Get first tile dimensions (may be different from last tile)
  nx0_grid_tile_dx = 1._num/(curr%array_of_tiles(1, 1, 1)%nx_grid_tile*dx)
  ny0_grid_tile_dy = 1._num/(curr%array_of_tiles(1, 1, 1)%ny_grid_tile*dy)
  nz0_grid_tile_dz = 1._num/(curr%array_of_tiles(1, 1, 1)%nz_grid_tile*dz)

  ! Parallelization on the tiles:
  ! Each tile will read the buffer and copy particles which belong
  ! to it into its own particle array

  !$OMP PARALLEL DO DEFAULT(NONE) SHARED(is, npid, ntilez, ntiley, ntilex,            &
  !$OMP nrecv_buf_tot, recvbuf, curr, dxs2, dys2, dzs2, nx0_grid_tile_dx,             &
  !$OMP ny0_grid_tile_dy, nz0_grid_tile_dz, x_min_local, y_min_local, z_min_local)    &
  !$OMP PRIVATE(ix, iy, iz, i, k, indx, indy, indz, curr_tile, nptile) COLLAPSE(3)    &
  !$OMP SCHEDULE(runtime)
  DO iz=1, ntilez! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex

        curr_tile=>curr%array_of_tiles(ix, iy, iz)
        !nptile=curr_tile%np_tile(1)

        ! If the tile is not at the boundary, no particle will be transfer
        IF (.NOT.curr_tile%subdomain_bound) CYCLE

        ! Add received particles to particle arrays
        DO i = 1, nrecv_buf_tot

          ! Get particle index in array of tile
          indx = MIN(FLOOR((recvbuf(i, 1)-x_min_local+dxs2)*(nx0_grid_tile_dx),       &
          idp)+1, ntilex)
          indy = MIN(FLOOR((recvbuf(i, 2)-y_min_local+dys2)*(ny0_grid_tile_dy),       &
          idp)+1, ntiley)
          indz = MIN(FLOOR((recvbuf(i, 3)-(z_min_local)+dzs2)*(nz0_grid_tile_dz),     &
          idp)+1, ntilez)

          ! If the particle is in the current tile
          IF ((indx.eq.ix).AND.(indy.eq.iy).AND.(indz.eq.iz)) THEN

            ! Sanity check for max number of particles in tile
            nptile = curr_tile%np_tile(1)+1! Current number of particles in the tile
            k  = curr_tile%npmax_tile! Max number of particles in the tile
            IF (nptile .GT. k) THEN
              ! Resize particle tile arrays if tile is full
              curr%are_tiles_reallocated(ix, iy, iz)=1
              CALL resize_particle_arrays(curr_tile, k, NINT(resize_factor*k+1, idp))
            ENDIF

            ! Finally, add particle to tile
            curr_tile%np_tile(1)=nptile
            curr_tile%part_x(nptile)  = recvbuf(i, 1)
            curr_tile%part_y(nptile)  = recvbuf(i, 2)
            curr_tile%part_z(nptile)  = recvbuf(i, 3)
            curr_tile%part_ux(nptile) = recvbuf(i, 4)
            curr_tile%part_uy(nptile) = recvbuf(i, 5)
            curr_tile%part_uz(nptile) = recvbuf(i, 6)
            curr_tile%part_gaminv(nptile) = recvbuf(i, 7)
            curr_tile%pid(nptile, 1:npid) = recvbuf(i, 8:7+npid)
            curr_tile%part_ex(nptile)  = 0._num
            curr_tile%part_ey(nptile)  = 0._num
            curr_tile%part_ez(nptile)  = 0._num
            curr_tile%part_bx(nptile)  = 0._num
            curr_tile%part_by(nptile)  = 0._num
            curr_tile%part_bz(nptile)  = 0._num

          ENDIF

        ENDDO

      ENDDO
    ENDDO
  ENDDO
  !!$OMP END PARALLEL DO

  ! Update total number of particle species
  curr%species_npart=curr%species_npart+nrecv_buf_tot

ENDIF

DEALLOCATE(recvbuf)

ENDDO

DEALLOCATE(mpi_npart)
DEALLOCATE(reqs)

! ________________________________________________________________________________________
! Checking
#if defined(DEBUG)
WRITE(0, *) " Checking after particle boundary conditions: start"
DO is=1, nspecies
curr=> species_parray(is)
DO iz=1, ntilez
  DO iy=1, ntiley
    DO ix=1, ntilex

      curr_tile=>curr%array_of_tiles(ix, iy, iz)
      nptile=curr_tile%np_tile(1)

      DO i=1, nptile
        partx=curr_tile%part_x(i)
        party=curr_tile%part_y(i)
        partz=curr_tile%part_z(i)
        partux=curr_tile%part_ux(i)
        partuy=curr_tile%part_uy(i)
        partuz=curr_tile%part_uz(i)
        gaminv=curr_tile%part_gaminv(i)
        partpid=curr_tile%pid(i, 1:npid)

        IF ((partx .LT. curr_tile%x_tile_min) .OR. (partx .GE. curr_tile%x_tile_max)  &
        .OR. (party .LT. curr_tile%y_tile_min) .OR. (party .GE. curr_tile%y_tile_max) &
        .OR. (partz .LT. curr_tile%z_tile_min) .OR. (partz .GE.                       &
        curr_tile%z_tile_max)) THEN


        WRITE(0, '("ERROR: particle outside the domain")')
        WRITE(0, '("Particle id:", I7, " of species ", I7)') i, is
        WRITE(0, '("In tile: ", I3, X, I3, X, I3)') ix, iy, iz
        WRITE(0, '("x:", E12.5, X, E12.5, X, E12.5)') curr_tile%x_tile_min, partx,    &
        curr_tile%x_tile_max
        WRITE(0, '("y:", E12.5, X, E12.5, X, E12.5)') curr_tile%y_tile_min, party,    &
        curr_tile%y_tile_max
        WRITE(0, '("z:", E12.5, X, E12.5, X, E12.5)') curr_tile%z_tile_min, partz,    &
        curr_tile%z_tile_max
        WRITE(0, *)

      ENDIF

    ENDDO

  ENDDO
ENDDO
ENDDO
ENDDO
WRITE(0, *) " Checking after particle boundary conditions: stop"
#endif

DEALLOCATE(partpid)
END SUBROUTINE particle_bcs_tiles_and_mpi_3d

END MODULE particle_boundary

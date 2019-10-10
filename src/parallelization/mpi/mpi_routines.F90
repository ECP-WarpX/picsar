! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! MPI_ROUTINES.F90
!
! This file contains functions for MPI.
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!
!> @brief
!> This module contains subroutines for MPI
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
MODULE mpi_routines
  USE shared_data
  USE fields
  USE mpi
  USE params
  IMPLICIT NONE
  REAL(num) :: start_time, end_time

  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine performs a minimal initialization for MPI.
  !> It calls the MPI init subroutine and determine the ranks
  !> and the number of processes.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE mpi_minimal_init()
    LOGICAL(isp) :: isinitialized
    INTEGER(isp) :: nproc_comm, rank_in_comm
    CALL MPI_INITIALIZED(isinitialized, errcode)
    IF (.NOT. isinitialized) THEN
      CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided, errcode)
    ENDIF
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
    nproc=INT(nproc_comm, idp)
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
    rank=INT(rank_in_comm, idp)
  END SUBROUTINE mpi_minimal_init


  SUBROUTINE mpi_minimal_init_fftw()
#if defined(FFTW)
    USE iso_c_binding
    USE mpi_fftw3, ONLY: fftw_mpi_init
    USE picsar_precision, ONLY: idp, isp
#endif
    LOGICAL(isp) :: isinitialized
    INTEGER(isp) :: nproc_comm, rank_in_comm
    INTEGER(isp) :: iret
    CALL MPI_INITIALIZED(isinitialized, errcode)
    IF (.NOT. isinitialized) THEN
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, errcode)
#if defined(FFTW)
      IF (provided >= MPI_THREAD_FUNNELED) THEN
        CALL DFFTW_INIT_THREADS(iret)
        fftw_threads_ok = .TRUE.
      ELSE
        fftw_threads_ok=.FALSE.
      ENDIF
      CALL FFTW_MPI_INIT()
#endif
    ENDIF
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
    nproc=INT(nproc_comm, idp)
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
    rank=INT(rank_in_comm, idp)
  END SUBROUTINE mpi_minimal_init_fftw

  ! ______________________________________________________________________________________
  !> @brief
  !> Minimal initialization when using Python.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  SUBROUTINE mpi_minimal_init_python(comm_in)
#if defined(FFTW)
    USE iso_c_binding
    USE mpi_fftw3, ONLY: fftw_mpi_init
    USE picsar_precision, ONLY: idp, isp
#endif
    LOGICAL(isp) :: isinitialized
    INTEGER(isp) :: nproc_comm, rank_in_comm
    INTEGER(idp), OPTIONAL, INTENT(IN) :: comm_in
    INTEGER(isp) :: iret

    CALL MPI_INITIALIZED(isinitialized, errcode)
    IF (.NOT. isinitialized) THEN
      CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided,  errcode)
    ENDIF

    IF (present(comm_in) .AND. comm_in .GE. 0) THEN
      CALL MPI_COMM_DUP(INT(comm_in, isp), comm, errcode)
    ELSE
      CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    ENDIF
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
    nproc=INT(nproc_comm, idp)

    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
    rank=INT(rank_in_comm, idp)
#if defined(FFTW)
    IF (.NOT. p3dfft_flag) THEN
      CALL DFFTW_INIT_THREADS(iret)
      IF( fftw_with_mpi) THEN
        CALL FFTW_MPI_INIT()
      ENDIF
    ENDIF
#endif


  END SUBROUTINE mpi_minimal_init_python


  ! ______________________________________________________________________________________
  !> @brief
  !> Get mpi topology when using warpp.
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2018
  ! ______________________________________________________________________________________
  SUBROUTINE get_neighbours_python


  
    proc_x_min = MODULO(x_coords-1,nprocx) + y_coords*nprocx +z_coords*nprocx*nprocy
    proc_x_max = MODULO(x_coords+1,nprocx) + y_coords*nprocx +z_coords*nprocx*nprocy
    proc_y_min = x_coords + MODULO(y_coords-1,nprocy)*nprocx +z_coords*nprocx*nprocy
    proc_y_max = x_coords + MODULO(y_coords+1,nprocy)*nprocx +z_coords*nprocx*nprocy
    proc_z_min = x_coords + y_coords*nprocx +MODULO(z_coords-1,nprocz)*nprocx*nprocy
    proc_z_max = x_coords + y_coords*nprocx +MODULO(z_coords+1,nprocz)*nprocx*nprocy

  END SUBROUTINE get_neighbours_python

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine creates the MPI process decomposition
  !> according to the specified topology and compite the related parameters.
  !
  !> @author
  !> Henri Vincenti
  !> Mathieu Lobet
  !
  !> @date
  !> Creation 2015
  !>
  !> Modification:
  !> Mathieu Lobet - March 2016: implementation of a random and home-made cartesian
  !>                             topology for efficiency comparison with the MPI
  !>                             Cartesian topology.
  ! ______________________________________________________________________________________
  SUBROUTINE setup_communicator
    INTEGER(isp), PARAMETER :: ndims = 3
    INTEGER(idp) :: idim
    INTEGER(isp) :: nproc_comm, dims(ndims), old_comm, ierr, neighb
    INTEGER(isp) :: proc_x_minsp, proc_x_maxsp, proc_y_minsp, proc_y_maxsp
    INTEGER(isp) :: proc_z_minsp, proc_z_maxsp
    LOGICAL(isp) :: periods(ndims), reorder, op
    INTEGER(isp) :: test_coords(ndims), rank_in_comm
    INTEGER(idp) :: ix, iy, iz
    INTEGER(idp) :: x_coords_neight, y_coords_neight, z_coords_neight
    INTEGER(idp) :: nxsplit, nysplit, nzsplit
    INTEGER(idp) :: rankyz
    INTEGER(idp) :: new_rank, old_rank
    INTEGER(idp), dimension(:, :, :), allocatable :: topo_array
    REAL(num) :: r


    nx_global=nx_global_grid-1
    ny_global=ny_global_grid-1
    nz_global=nz_global_grid-1

    !!! --- NB: CPU Split performed on number of grid points (not cells)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc_comm, ierr)
    nproc=INT(nproc_comm, idp)

    ! With fftw can only do CPU split with respect to Z direction (X in C-order)
    IF (fftw_with_mpi .AND. .NOT. fftw_hybrid) THEN
      nprocx=1
      nprocy=1
      nprocz=nproc
    ENDIF
    ! check if c_dim = 2 and fftw_mpi_transpose = false
    IF(c_dim ==2 .AND. fftw_mpi_transpose .EQV. .TRUE.) THEN
      fftw_mpi_transpose = .FALSE.
      WRITE(0, *) 'WARING fftw_mpi_transpose flag unavailable in dim 2'
      WRITE(0, *) 'fftw_mpi_transpose set to .FALSE.'
      fftw_mpi_transpose = .FALSE.
    ENDIF

    dims = (/nprocz, nprocy, nprocx/)

    ! Initial CPU split sanity check
    IF ((nprocx .EQ. 0) .AND. (nprocy .EQ. 0) .AND. (nprocz .EQ. 0)) THEN
      CALL MPI_DIMS_CREATE(nproc_comm, ndims, dims, errcode)
      nprocx = INT(dims(3), idp)
      nprocy = INT(dims(2), idp)
      nprocz = INT(dims(1), idp)
    ENDIF

  IF (nproc .NE. nprocx*nprocy*nprocz) THEN
      IF (rank .EQ. 0) THEN
        WRITE(0, *) '*** ERROR ***'
        WRITE(0, *) 'nprocx*nprocy*nprocz =/ # of MPI processes'
        WRITE(0, *) ' Check input file '
        WRITE(0, *) ' Total number of processors:', nproc
        WRITE(0, *) ' Number of processors in each direction:', nprocx, nprocy,       &
        nprocz
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF
  ENDIF

  IF(c_dim == 3) THEN
    IF (nx_global_grid .LT. nxguards .OR. ny_global_grid .LT. nyguards .OR.           &
      nz_global_grid .LT. nzguards) THEN
      IF (rank .EQ. 0) THEN
        WRITE(0, *) '*** ERROR ***'
        WRITE(0, *) 'Simulation domain is too small.'
        WRITE(0, *) 'nx_global_grid', nx_global_grid, 'nxguards', nxguards
        WRITE(0, *) 'ny_global_grid', ny_global_grid, 'nyguards', nyguards
        WRITE(0, *) 'nz_global_grid', nz_global_grid, 'nzguards', nzguards
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF
  ELSE IF(c_dim == 2) THEN
    IF (nx_global_grid .LT. nxguards .OR. nz_global_grid .LT. nzguards) THEN
      IF (rank .EQ. 0) THEN
        WRITE(0, *) '*** ERROR ***'
        WRITE(0, *) 'Simulation domain is too small.'
        WRITE(0, *) 'nx_global_grid', nx_global_grid, 'nxguards', nxguards
        WRITE(0, *) 'nz_global_grid', nz_global_grid, 'nzguards', nzguards
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF
  ENDIF
  IF (nprocx * nprocy * nprocz .GT. 0) THEN
    nxsplit = nx_global_grid / nprocx
    nysplit = ny_global_grid / nprocy
    nzsplit = nz_global_grid / nprocz
    IF(c_dim == 3 ) THEN
      IF (nxsplit .LT. nxguards .OR. nysplit .LT. nyguards .OR. nzsplit .LT. &
         nzguards)  THEN
        IF (rank .EQ. 0) THEN
          WRITE(0, *) 'WRONG CPU SPLIT nlocal<nguards 1'
          WRITE(0, *) nxsplit, nysplit, nzsplit
          WRITE(0, *) nx_global_grid, ny_global_grid, nz_global_grid
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
      ENDIF
    ELSE IF(c_dim ==2) THEN
      IF (nxsplit .LT. nxguards .OR. nzsplit .LT.  nzguards)  THEN
        IF (rank .EQ. 0) THEN
          WRITE(0, *) 'WRONG CPU SPLIT nlocal<nguards 2'
          WRITE(0, *) nxsplit, nzsplit
          WRITE(0, *) nx_global_grid, nz_global_grid
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
      ENDIF
    ENDIF
  ENDIF

periods = .FALSE.
reorder = .TRUE.

! Random topology ----------------------------------------
IF (topology == 2) THEN

  IF (rank .EQ. 0) THEN
    WRITE(0, *) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
  ENDIF

  allocate(topo_array(0:nprocx-1, 0:nprocy-1, 0:nprocz-1))

  IF (rank.EQ.0) THEN

    ! Creation of a random array for the topology
    DO iz=0, nprocz-1
      DO iy=0, nprocy-1
        DO ix=0, nprocx-1
          topo_array(ix, iy, iz) = INT(ix + (iy)*nprocx + (iz)*nprocy*nprocx, idp)
          !print*, ix + (iy)*nprocx + (iz)*nprocy*nprocx
        ENDDO
      ENDDO
    ENDDO
    DO iz =0, nprocz-1
      DO iy=0, nprocy-1
        DO ix=0, nprocx-1
          CALL RANDOM_NUMBER(r)
          x_coords = INT(r*(nprocx-1))
          ! x_coords = INT(RAN(seed)*(nprocx-1))
          CALL RANDOM_NUMBER(r)
          y_coords = INT(r*(nprocy-1))
          CALL RANDOM_NUMBER(r)
          z_coords = INT(r*(nprocz-1))
          rankyz = topo_array(ix, iy, iz)
          topo_array(ix, iy, iz) = topo_array(x_coords, y_coords, z_coords)
          topo_array(x_coords, y_coords, z_coords) = rankyz
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Sharing of the topology
  CALL MPI_BCAST(topo_array, INT(nprocx*nprocy*nprocz, isp), MPI_INTEGER8, 0_isp,     &
  comm, errcode)

  !IF (rank==3) print*, 'rank', rank, 'topo_array', topo_array

  ! Each processor determine their coordinates
  DO iz =0, nprocz-1
    DO iy=0, nprocy-1
      DO ix=0, nprocx-1

        IF (rank .EQ. topo_array(ix, iy, iz)) THEN

          x_coords = ix
          y_coords = iy
          z_coords = iz

        ENDIF

      ENDDO
    ENDDO
  ENDDO

  x_min_boundary = .FALSE.
  x_max_boundary = .FALSE.
  IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
  IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.

  y_min_boundary = .FALSE.
  y_max_boundary = .FALSE.
  IF (y_coords .EQ. 0) y_min_boundary = .TRUE.
  IF (y_coords .EQ. nprocy - 1) y_max_boundary = .TRUE.

  z_min_boundary = .FALSE.
  z_max_boundary = .FALSE.
  IF (z_coords .EQ. 0) z_min_boundary = .TRUE.
  IF (z_coords .EQ. nprocz - 1) z_max_boundary = .TRUE.

  neighbour = MPI_PROC_NULL
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1

        IF (x_coords + ix .LT. 0) THEN
          x_coords_neight = x_coords + nprocx-1
        ELSE IF (x_coords + ix .GT. nprocx-1) THEN
          x_coords_neight = x_coords - (nprocx-1)
        ELSE
          x_coords_neight = x_coords + ix
        END IF

        IF (y_coords + iy .LT. 0) THEN
          y_coords_neight = y_coords + nprocy-1
        ELSE IF (y_coords + iy .GT. nprocy-1) THEN
          y_coords_neight = y_coords - (nprocy-1)
        ELSE
          y_coords_neight = y_coords + iy
        END IF

        IF (z_coords + iz .LT. 0) THEN
          z_coords_neight = z_coords + nprocz-1
        ELSE IF (z_coords + iz .GT. nprocz-1) THEN
          z_coords_neight = z_coords - (nprocz-1)
        ELSE
          z_coords_neight = z_coords + iz
        END IF

        neighbour(ix, iy, iz) = topo_array(x_coords_neight, y_coords_neight,          &
        z_coords_neight)
      ENDDO
    ENDDO
  ENDDO

  proc_x_max = neighbour(1, 0, 0)
  proc_x_min = neighbour(-1, 0, 0)
  proc_y_max = neighbour(0, 1, 0)
  proc_y_min = neighbour(0, -1, 0)
  proc_z_max = neighbour(0, 0, 1)
  proc_z_min = neighbour(0, 0, -1)

  ! Topology - processors are ordered according to their ranks ---------------
ELSE IF (topology == 1) THEN

  IF (rank .EQ. 0) THEN
    WRITE(0, *) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
  ENDIF

  ! We first fill the x direction, then y and finally z
  x_coords = MOD(rank, nprocx)
  y_coords = MOD((rank-x_coords)/nprocx, nprocy)
  z_coords = (rank-x_coords - y_coords*nprocx)/(nprocx*nprocy)

  x_min_boundary = .FALSE.
  x_max_boundary = .FALSE.
  IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
  IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.

  y_min_boundary = .FALSE.
  y_max_boundary = .FALSE.
  IF (y_coords .EQ. 0) y_min_boundary = .TRUE.
  IF (y_coords .EQ. nprocy - 1) y_max_boundary = .TRUE.

  z_min_boundary = .FALSE.
  z_max_boundary = .FALSE.
  IF (z_coords .EQ. 0) z_min_boundary = .TRUE.
  IF (z_coords .EQ. nprocz - 1) z_max_boundary = .TRUE.

  neighbour = MPI_PROC_NULL
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1

        IF (x_coords + ix .LT. 0) THEN
          x_coords_neight = x_coords + nprocx-1
        ELSE IF (x_coords + ix .GT. nprocx-1) THEN
          x_coords_neight = x_coords - (nprocx-1)
        ELSE
          x_coords_neight = x_coords + ix
        END IF

        IF (y_coords + iy .LT. 0) THEN
          y_coords_neight = y_coords + nprocy-1
        ELSE IF (y_coords + iy .GT. nprocy-1) THEN
          y_coords_neight = y_coords - (nprocy-1)
        ELSE
          y_coords_neight = y_coords + iy
        END IF

        IF (z_coords + iz .LT. 0) THEN
          z_coords_neight = z_coords + nprocz-1
        ELSE IF (z_coords + iz .GT. nprocz-1) THEN
          z_coords_neight = z_coords - (nprocz-1)
        ELSE
          z_coords_neight = z_coords + iz
        END IF

        neighbour(ix, iy, iz) = x_coords_neight + (nprocx)*(y_coords_neight) +        &
        (nprocx*nprocy)*(z_coords_neight)

        print*, 'rank', rank, ix, iy, iz, neighbour(ix, iy, iz)

      ENDDO
    ENDDO
  ENDDO

  proc_x_max = neighbour(1, 0, 0)
  proc_x_min = neighbour(-1, 0, 0)
  proc_y_max = neighbour(0, 1, 0)
  proc_y_min = neighbour(0, -1, 0)
  proc_z_max = neighbour(0, 0, 1)
  proc_z_min = neighbour(0, 0, -1)

  ! Default topology - Cartesian topology
ELSE
  ! Set boundary to be periodic in x, y, z for particles and fields by default
  periods(c_ndims) = .TRUE.
  periods(c_ndims-1) = .TRUE.
  periods(c_ndims-2) = .TRUE.

  ! Create an array of rank corresponding to the current topology
  !ALLOCATE(rank_array(nprocs))
  old_rank = rank
  !old_rank_array(old_rank) = rank

  ! Creation of a Cartesian topology by using the original communicator
  old_comm = comm
  CALL MPI_CART_CREATE(old_comm, ndims, dims, periods, reorder, comm, errcode)
  CALL MPI_COMM_FREE(old_comm, errcode)
  CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
  CALL MPI_CART_COORDS(comm, rank_in_comm, ndims, coordinates, errcode)
  rank=INT(rank_in_comm, idp)
  ! Determine MPI ranks at the boundaries
  CALL MPI_CART_SHIFT(comm, 2_isp, 1_isp, proc_x_minsp, proc_x_maxsp, errcode)
  CALL MPI_CART_SHIFT(comm, 1_isp, 1_isp, proc_y_minsp, proc_y_maxsp, errcode)
  CALL MPI_CART_SHIFT(comm, 0_isp, 1_isp, proc_z_minsp, proc_z_maxsp, errcode)

  proc_x_min=INT(proc_x_minsp, idp)
  proc_x_max=INT(proc_x_maxsp, idp)
  proc_y_min=INT(proc_y_minsp, idp)
  proc_y_max=INT(proc_y_maxsp, idp)
  proc_z_min=INT(proc_z_minsp, idp)
  proc_z_max=INT(proc_z_maxsp, idp)

  IF (rank .EQ. 0) THEN
    WRITE(0, *) 'Cartesian topology'
    WRITE(0, *) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
    WRITE(0, *)
  ENDIF

  ! ------------------------------------------------------------
  ! Checking of the new topology
  CALL MPI_BARRIER(comm, errcode)
  new_rank = rank
  ! We first fill the x direction, then y and finally z
  x_coords = MOD(rank, nprocx)
  y_coords = MOD((rank-x_coords)/nprocx, nprocy)
  z_coords = (rank-x_coords - y_coords*nprocx)/(nprocx*nprocy)
  CALL MPI_BARRIER(comm, errcode)
  ! -------------------------------------------------------------

  x_coords = coordinates(c_ndims)
  x_min_boundary = .FALSE.
  x_max_boundary = .FALSE.
  IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
  IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.

  y_coords = coordinates(c_ndims-1)
  y_min_boundary = .FALSE.
  y_max_boundary = .FALSE.
  IF (y_coords .EQ. 0) y_min_boundary = .TRUE.
  IF (y_coords .EQ. nprocy - 1) y_max_boundary = .TRUE.

  z_coords = coordinates(c_ndims-2)
  z_min_boundary = .FALSE.
  z_max_boundary = .FALSE.
  IF (z_coords .EQ. 0) z_min_boundary = .TRUE.
  IF (z_coords .EQ. nprocz - 1) z_max_boundary = .TRUE.

  neighbour = MPI_PROC_NULL
  DO iz = -1, 1
    DO iy = -1, 1
      DO ix = -1, 1
        test_coords = coordinates
        test_coords(1) = test_coords(1)+iz
        test_coords(2) = test_coords(2)+iy
        test_coords(3) = test_coords(3)+ix
        op = .TRUE.
        DO idim = 1, ndims
          IF ((test_coords(idim) .LT. 0 .OR. test_coords(idim) .GE. dims(idim)) .AND. &
          .NOT. periods(idim)) op = .FALSE.
        ENDDO
        IF (op) THEN
          ! Then determine the rank of the neighbor processor
          CALL MPI_CART_RANK(comm, test_coords, neighb, errcode)
          neighbour(ix, iy, iz)=INT(neighb, idp)
        ENDIF
      ENDDO
    ENDDO
  ENDDO

ENDIF

IF(absorbing_bcs) THEN
  CALL get_non_periodic_mpi_bcs()
ENDIF

END SUBROUTINE setup_communicator

! ______________________________________________________________________________________
!> @brief
!> This routine subroutine sets proc_x/y/z_min/max to  MPI_PROC_NULL for mpi at the edge 
!> of the domain when using absorbing_bcs
!> @ author
!> H. Kallala
!> H. Vincenti
!> 2018
! ______________________________________________________________________________________
SUBROUTINE get_non_periodic_mpi_bcs()
  IF(absorbing_bcs_x) THEN
    IF(x_min_boundary) proc_x_min = MPI_PROC_NULL
    IF(x_max_boundary) proc_x_max = MPI_PROC_NULL
  ENDIF
  IF(absorbing_bcs_y) THEN
    IF(y_min_boundary) proc_y_min = MPI_PROC_NULL
    IF(y_max_boundary) proc_y_max = MPI_PROC_NULL
  ENDIF
  IF(absorbing_bcs_z) THEN
    IF(z_min_boundary) proc_z_min = MPI_PROC_NULL
    IF(z_max_boundary) proc_z_max = MPI_PROC_NULL
  ENDIF
END SUBROUTINE get_non_periodic_mpi_bcs

! ______________________________________________________________________________________
!> @brief
!> This routine creates mpi groups and corresponding communicators while using
!> fftw_hybrid=.TRUE. or hybrid=.TRUE.
!> @ author
!> H. Kallala
!> H. Vincenti
!> 2018
! ______________________________________________________________________________________
SUBROUTINE setup_groups
#if defined(FFTW)
USE group_parameters, ONLY: cell_x_max_g, cell_x_min_g, cell_y_max_g, cell_y_min_g,  &
  cell_z_max_g, cell_z_min_g, group_sizes, group_y_max_boundary,                     &
  group_y_min_boundary, group_z_max_boundary, group_z_min_boundary,                  &
  is_group_x_boundary_max, is_group_x_boundary_min, is_group_y_boundary_max,         &
  is_group_y_boundary_min, is_group_z_boundary_max, is_group_z_boundary_min,         &
  is_on_boundary_group_y, is_on_boundary_group_z, local_rank, local_size,            &
  mpi_comm_group_id, mpi_group_id, mpi_ordered_comm_world, mpi_root_comm,            &
  mpi_root_group, mpi_world_group, nx_group, nx_group_global, nx_group_global_array, &
  nx_group_global_grid, nx_group_grid, ny_group, ny_group_global,                    &
  ny_group_global_array, ny_group_global_grid, ny_group_grid, nz_group,              &
  nz_group_global, nz_group_global_array, nz_group_global_grid, nz_group_grid,       &
  p3d_fend, p3d_fsize, p3d_fstart, p3d_iend, p3d_isize, p3d_istart, root_rank,       &
  root_size, which_group, x_group_coords, x_max_group, x_min_group, y_group_coords,  &
  y_max_group, y_min_group, z_group_coords, z_max_group, z_min_group
USE iso_c_binding
#endif
USE mpi
#if defined(FFTW)
USE mpi_fftw3, ONLY: alloc_local, fftw_mpi_local_size_2d, fftw_mpi_local_size_3d,    &
  fftw_mpi_local_size_3d_transposed, local_nx, local_nx_tr, local_ny, local_ny_tr,   &
  local_nz, local_nz_tr, local_x0, local_x0_tr, local_y0, local_y0_tr, local_z0,     &
  local_z0_tr
#endif
USE mpi_type_constants, ONLY: mpidbl
#if defined(P3DFFT)
USE p3dfft
#endif
USE picsar_precision, ONLY: idp, isp
USE shared_data, ONLY: absorbing_bcs, c_dim, cell_x_max, cell_x_min, comm, dx, dy,   &
  dz, errcode, fftw_mpi_transpose, ix_max_r, ix_min_r, iy_max_r, iy_min_r, iz_max_r, &
  iz_min_r, nb_group, nb_group_x, nb_group_y, nb_group_z, nproc, nprocx, nprocy,     &
  nprocz, nx, nxg_group, ny, ny_global, nyg_group, nz, nz_global, nzg_group,         &
  p3dfft_flag, rank, x, x_coords, xmax, xmin, y, y_coords, y_max_boundary,           &
  y_min_boundary, ymax, ymin, z, z_coords, z_max_boundary, z_min_boundary, zmax,     &
  zmin
INTEGER(isp) :: ierr
INTEGER(idp) :: group_size
INTEGER(isp), ALLOCATABLE, DIMENSION(:)    ::  grp_comm
INTEGER(isp)                                :: roots_comm
INTEGER(idp)  :: i,j,temp
INTEGER(idp), DIMENSIOn(:), ALLOCATABLE :: all_iz_min_global, all_iz_max_global
INTEGER(idp)                            :: iz_min_global, iz_max_global,iy_min_global, &
iy_max_global
#if defined(P3DFFT)
INTEGER(isp)                            :: pdims(2)
#endif
INTEGER(idp) , ALLOCATABLE, DIMENSION(:)  :: all_iy_min_global,all_iy_max_global
INTEGER(isp) :: color
INTEGER(isp)                              :: key, key_roots,color_roots
LOGICAL(isp)                               :: is_in_place

#if defined(FFTW)
#if defined(DEBUG)
  WRITE(0, *) "setup_groups : start"
#endif

  IF(c_dim == 2) THEN
    ny_global = 1_idp
    nyguards = 0_idp
    nyjguards = 0_idp
  ENDIF

  ! - Create the base group upon which all other groups are defined
  CALL MPI_COMM_GROUP(comm, mpi_world_group, errcode)

  ! - Computes number of groups along x, y and z
  IF(.NOT. p3dfft_flag) THEN
    nb_group_y = nprocy
    nyg_group = nyguards
  ENDIF
  nb_group_x = nprocx

  ! - If 2D case set p3dfft to .FALSE. by default  (p3dfft is a 3D library)
  IF(c_dim==2 .AND. p3dfft_flag)  THEN
     IF(rank==0) THEN
       WRITE(*,*)"ERROR cant run with p3dfft in 2d case"
       WRITE(*,*)"p3dfft = .FALSE."
     ENDIF
     p3dfft_flag = .FALSE.
  ENDIF

  ! - If p3dfft is enabled, fftw_mpi_transpose is set to .FALSE.
  ! - by default (all MPI communications are handled by p3dfft in that case)
  IF(p3dfft_flag) THEN
    fftw_mpi_transpose = .FALSE.
  ENDIF

  ! - Sanity check: number of groups has to divide number of procs along each direction
  ! - If not --> abort simulation
  IF(MODULO(nprocx,nb_group_x) .NE. 0 .OR. MODULO(nprocy,nb_group_y) .NE. 0  &
  .OR. MODULO(nprocz,nb_group_z) .NE. 0 ) THEN
    IF(rank==0)   THEN
      WRITE(*,*)"ERROR please set a number of groups multiple of nproc along each direction"
    ENDIF
    CALL MPI_ABORT(comm,errcode,ierr)
  ENDIF

  ! - Computes total number of groups
  nb_group = nb_group_z*nb_group_y*nb_group_x
  ! - Computes number of ranks per group
  group_size = (nprocx*nprocy*nprocz)/nb_group
  ! - Computes number of ranks per group along x
  group_sizes(1) =nprocx/nb_group_x
  ! - Computes number of ranks per group along y
  group_sizes(2)=nprocy/nb_group_y
  ! - Computes number of ranks per group along z
  group_sizes(3) = nprocz/nb_group_z
  ! - Computes coordinates of current group in the global set of groups along z
  z_group_coords = z_coords/group_sizes(3)
  ! - Computes coordinates of current group in the global set of groups along y
  y_group_coords = y_coords/group_sizes(2)
  ! - Computes coordinates of current group in the global set of groups along x
  x_group_coords = x_coords/group_sizes(1)
  ! - Init the group id of current group based on its coordinates in the global
  ! set of  groups
  which_group =x_group_coords+y_group_coords*nb_group_x+ &
  z_group_coords*nb_group_x*nb_group_y


  ! - Array allocations
  ALLOCATE(grp_comm(nb_group))
  ALLOCATE(mpi_group_id(nb_group),mpi_comm_group_id(nb_group))


  !-- Computing the rank of the process in its MPI group
  !-- Key is initialized so that processes in the new MPI group communicators
  !-- have the same neighbors than the ones in the global communicator comm
  key = MODULO(z_coords,group_sizes(3))*group_sizes(1)*group_sizes(2) + &
        MODULO(y_coords,group_sizes(2))*group_sizes(1) + &
        MODULO(x_coords,group_sizes(1))

  ! - Init a local root rank for each group identified by color_roots
  ! - Local root ranks are used to communicate information between MPI groups
  ! - By default, the local root rank of a MPI group is
  ! - the minimum of all ranks contained in that group which is spoted by key==0 here
  color_roots=MPI_UNDEFINED
  IF(key==0) color_roots = 1_isp  ! 1 != mpi_undefined

  ! ranks of procs inside roots_comm will be in the same order as groups
  key_roots = z_group_coords*nb_group_y*nb_group_x + y_group_coords*nb_group_x &
  + x_group_coords

  roots_comm = MPI_COMM_NULL
  root_rank = MPI_PROC_NULL
  mpi_root_comm = MPI_COMM_NULL

  ! -- Create a MPI communicator associated to the local roots
  CALL MPI_COMM_SPLIT(comm,color_roots, key_roots,roots_comm,errcode)

  ! - For each process in the root_communicator get:
  ! - (i) The size of the root communicator
  ! - (ii) The rank of current process in the communicator
  ! - (iii) Init mpi_root_comm and mpi_root_group variables
  IF( roots_comm .NE. MPI_COMM_NULL) THEN
    ! - Duplicate roots_comm into mpi_root_com
    CALL MPI_COMM_DUP(roots_comm,mpi_root_comm,errcode)
    ! - Create a MPI group associated to the MPI communicatorof local root ranks
    CALL MPI_COMM_GROUP(mpi_root_comm,mpi_root_group,errcode)
    ! - Get mpi ranks in mpi_root_comm
    CALL MPI_COMM_RANK(mpi_root_comm, root_rank, errcode)
    ! - Get mpi size in mpi_root_comm
    CALL MPI_COMM_SIZE(mpi_root_comm, root_size, errcode)
    ! - Free old roots_comm and keep mpi_root_comm
    CALL MPI_COMM_FREE(roots_comm, errcode)
  ENDIF

  ! - Init group variables and communicators
  DO i=1, nb_group
    mpi_comm_group_id(i) = MPI_COMM_NULL
    mpi_group_id(i) = MPI_GROUP_NULL
    grp_comm(i)  = MPI_COMM_NULL
  ENDDO


  ! -- Create groups communicators and groups
  ! -- mpi taks in the same group have the same value of which_group which is
  ! -- unique to each group
  ! -- Mpi processes of the same group share the same value of
  ! -- colors(which_group) so mpi_comm_world will be split into nb_group
  ! -- key enables to arrange ranks in the group the same way they are in the
  ! -- ordered grid
  color = which_group
  ! -- Create communicator associated to each group
  CALL MPI_COMM_SPLIT(comm,color,key,grp_comm(which_group+1),errcode)
  ! -- Duplicate old communicators to the variables mpi_comm_group_id,
  ! -- mpi_group_id
  DO i = 1, nb_group

    IF (grp_comm(i) .NE. MPI_COMM_NULL) THEN
      CALL MPI_COMM_DUP(grp_comm(i),mpi_comm_group_id(i),errcode)
      ! -- Create a MPI group associated to the local MPI communicator
      CALL MPI_COMM_GROUP(mpi_comm_group_id(i), mpi_group_id(i), errcode)
      ! - Get mpi_comm_group_id(i) communicator size
      CALL MPI_COMM_SIZE(mpi_comm_group_id(i), local_size, errcode)
      ! - get rank of current process in mpi_comm_group_id(i)
      CALL MPI_COMM_RANK(mpi_comm_group_id(i), local_rank, errcode)
      ! - For each group, release old (& temporary) communicators
      ! - And keep only new communicators and groups
      CALL MPI_COMM_FREE(grp_comm(i), errcode)
    ENDIF
  ENDDO
  ! - Create and ordered comm world  to compute communication scheduling
  ! -- key is the topological rank of current mpi if mpi_comm_world would have
  ! -- been cartesian
  key = z_coords*nprocx*nprocy+ y_coords*nprocx + x_coords
  color_roots = MPI_UNDEFINED
  color_roots = 1_isp
  CALL MPI_COMM_SPLIT(comm,color_roots,key,mpi_ordered_comm_world,errcode)

  ! - Detect if current rank is on -z boundary of its group
  ! - Will be replaced soon by calls to MPI_CART
  IF(local_rank/group_sizes(2)==0) THEN
    is_on_boundary_group_z = .TRUE.
    group_z_min_boundary = .TRUE.
  ENDIF
  ! - Detect if current rank is on +z boundary of its group
  ! - Will be replaced soon by calls to MPI_CART
  IF(local_rank/group_sizes(2)==group_sizes(3)-1) THEN
    is_on_boundary_group_z = .TRUE.
    group_z_max_boundary = .TRUE.
  ENDIF

  ! - Detect if current rank is on +y boundary of its group
  ! - Will be replaced soon by calls to MPI_CART
  IF(MODULO(local_rank,int(group_sizes(2),isp))==0) THEN
    is_on_boundary_group_y = .TRUE.
    group_y_min_boundary = .TRUE.
  ENDIF

  ! - Detect if current rank is on -y boundary of its group
  ! - Will be replaced soon by calls to MPI_CART
  IF(MODULO(local_rank,int(group_sizes(2),isp))==group_sizes(2)-1) THEN
    is_on_boundary_group_y = .TRUE.
    group_y_max_boundary = .TRUE.
  ENDIF

  ! - Init number of cells/group along each direction (EXCLUDING guard cells)
  nx_group_global = nx ! By default there is no group along X in current version
  ny_group_global = ny_global/nb_group_y
  nz_group_global = nz_global/(nb_group_z)

  ! - Sanity check: if the number of groups does not exactly divide the number of cells
  ! - per group  along z then add one cell to the  remaining cells
  ! - of each one of the the  nz_global - nb_group_z*nz_group_global last groups
  IF(nz_global .NE. nz_group_global*(nb_group_z)) THEN
    temp = nz_global - nb_group_z*nz_group_global
    IF(nb_group_z-1-z_group_coords .LT. temp) THEN
      nz_group_global=nz_group_global+1
    ENDIF
  ENDIF
  ! - Sanity check: if the number of groups does not exactly divide the number of cells
  ! - per group  along y then add one cell to the remaining cells
  ! - of each one of the ny_global - nb_group_y*ny_group_global last groups
  IF(ny_global .NE. ny_group_global*(nb_group_y)) THEN
    temp = ny_global - nb_group_y*ny_group_global
    IF(nb_group_y-1-y_group_coords .LT. temp) THEN
      ny_group_global=ny_group_global+1
    ENDIF
  ENDIF

  ! - Array allocations and init
  ALLOCATE(nz_group_global_array(nb_group))
  ALLOCATE(ny_group_global_array(nb_group))
  ALLOCATE(nx_group_global_array(nb_group))
  y_min_group=ymin
  y_max_group=ymax
  x_min_group=xmin
  x_max_group=xmax
  z_min_group=zmin
  z_max_group=zmax

  ! - GATHER number of cells per groups along each direction from each group local root
  IF(mpi_root_comm .NE. MPI_COMM_NULL) THEN ! If == MPI_COMM_NULL current rank is not root
    CALL MPI_ALLGATHER(nz_group_global, 1_isp, MPI_INTEGER8,nz_group_global_array, &
    1_isp,MPI_INTEGER8, mpi_root_comm, errcode)
    CALL MPI_ALLGATHER(ny_group_global, 1_isp, MPI_INTEGER8,ny_group_global_array, &
    1_isp,MPI_INTEGER8, mpi_root_comm, errcode)
    CALL MPI_ALLGATHER(nx_group_global, 1_isp, MPI_INTEGER8,nx_group_global_array, &
    1_isp,MPI_INTEGER8, mpi_root_comm, errcode)

    ! - Computes z_min_group/z_max_group for the group containing current local root rank
    DO i=1, z_group_coords
      j = (i-1)*nb_group_x*nb_group_y+y_group_coords*nb_group_x+x_group_coords+1
      z_min_group = z_min_group + nz_group_global_array(j)*dz
    ENDDO
    z_max_group = z_min_group + nz_group_global*dz

    ! - Computes y_min_group/y_max_group for the group containing current local root rank
    DO i=1, y_group_coords
      j = z_group_coords*nb_group_x*nb_group_y+(i-1)*nb_group_x+x_group_coords+1
      y_min_group = y_min_group +ny_group_global_array(j)*dy
    ENDDO
    y_max_group = y_min_group + ny_group_global*dy

    ! - Computes x_min_group/x_max_group for the group containing current local root rank
    DO i=1, x_group_coords
      j = z_group_coords*nb_group_x*nb_group_y+y_group_coords*nb_group_x+i
      x_min_group = x_min_group +nx_group_global_array(j)*dx
    ENDDO
    x_max_group = x_min_group + nx_group_global*dx
  ENDIF

  ! - In each group, broadcast group dimensions and positions from its local root rank
  ! - to all other ranks
  DO i=1, nb_group
    IF(mpi_comm_group_id(i)  .NE. MPI_COMM_NULL) THEN
      CALL MPI_BCAST(z_min_group, 1_isp, mpidbl, 0_isp, mpi_comm_group_id(i),     &
      errcode)
      CALL MPI_BCAST(z_max_group, 1_isp, mpidbl, 0_isp, mpi_comm_group_id(i),     &
      errcode)
      CALL MPI_BCAST(y_min_group, 1_isp, mpidbl, 0_isp,mpi_comm_group_id(i),      &
      errcode)
      CALL MPI_BCAST(y_max_group, 1_isp, mpidbl, 0_isp,mpi_comm_group_id(i),      &
      errcode)
      CALL MPI_BCAST(x_min_group, 1_isp, mpidbl, 0_isp,mpi_comm_group_id(i),      &
      errcode)
      CALL MPI_BCAST(x_max_group, 1_isp, mpidbl, 0_isp,mpi_comm_group_id(i),      &
      errcode)
      CALL MPI_BCAST(nz_group_global_array,INT(nb_group,isp),MPI_INTEGER8,        &
      0_isp,mpi_comm_group_id(i),errcode)
      CALL MPI_BCAST(nx_group_global_array,INT(nb_group,isp),MPI_INTEGER8,       &
      0_isp,mpi_comm_group_id(i),errcode)
      CALL MPI_BCAST(ny_group_global_array,INT(nb_group,isp),MPI_INTEGER8,        &
      0_isp,mpi_comm_group_id(i),errcode)
    ENDIF
  ENDDO

  ! - For 2D case, number of guard cells per group is 0 along y
  IF(c_dim ==  2 ) THEN
    nyg_group = 0
  ENDIF

  ! - Update number of grid points per group along each direction
  nx_group_global_grid = nx_group_global+1
  ny_group_global_grid = ny_group_global+1
  nz_group_global_grid = nz_group_global+1

  ! - Get total number of cells per group
  nx_group = nx_group_global + 2*nxg_group
  ny_group = ny_group_global + 2*nyg_group
  nz_group = nz_group_global + 2*nzg_group

  ! - Get total number of grid points per group
  nx_group_grid = nx_group_global_grid + 2*nxg_group
  ny_group_grid = ny_group_global_grid + 2*nyg_group
  nz_group_grid = nz_group_global_grid + 2*nzg_group

  ! - For each group
  DO i=1, nb_group
    IF(mpi_comm_group_id(i)  .NE. MPI_COMM_NULL) THEN
      ! Case 1: FFTW-MPI is used for global FFT transpositions among groups
      IF(.NOT. p3dfft_flag) THEN
        ! - 3D case
        IF(c_dim == 3) THEN
          ! Fourier arrays have transposed dimensions
          IF(fftw_mpi_transpose) THEN
            ! - Init starting index/local size of FFT input arrays along Z
            ! - and starting index/local size of FFT output arrays along Y 
            ! - (transposed case)
            alloc_local = fftw_mpi_local_size_3d_transposed(nz_group, ny_group,         &
            nx_group/2+1, mpi_comm_group_id(i), local_nz, local_z0, local_nz_tr,        &
            local_z0_tr)
            ! - Init starting index/ local size of FFT input arrays along Y
            local_y0=0 
            local_ny=ny_group
            ! - Sanity check:  if transposition leads to local_ny > nprocy in group
            ! - in that case, local_ny can be 0 --> Abort
            ! - Same if local_nz<nprocz in group. This can lead to local_nz=0
            ! - In that case --> Abort as well
            IF(local_nz .EQ. 0_idp .OR. local_ny .EQ. 0_idp) THEN
              WRITE(0,*) 'ERROR local_ny or local_nz = 0 in rank ',rank
              CALL MPI_ABORT(comm,errcode,ierr)
            ENDIF
            ! - Init starting index and local size of output FFT arrays along Z
            local_y0_tr=0
            local_ny_tr=nz_group
          ! Fourier arrays have same dimensions than real arrays
          ELSE  
            alloc_local = FFTW_MPI_LOCAL_SIZE_3D(nz_group, ny_group, nx_group/2+1,      &
            mpi_comm_group_id(i), local_nz, local_z0)
            IF(local_nz .EQ. 0_idp ) THEN
              WRITE(0,*) 'ERROR local_nz = 0 in rank ',rank
              CALL MPI_ABORT(comm,errcode,ierr)
            ENDIF
            local_ny = ny_group 
            local_y0 =0_idp
            ! - Init starting indexes and local sizes of output FFT arrays along Y and Z
            local_ny_tr=local_ny
            local_y0_tr=local_y0
            local_nz_tr=local_nz
            local_z0_tr=local_z0
          ENDIF
          ! - Init starting index and local size of local input FFT arrays along X
          local_nx = 2_idp*(nx_group/2+1) 
          local_x0 =0_idp
          ! - Init starting index and local size of local output FFT arrays along X
          local_nx_tr=nx_group/2+1
          local_x0_tr=local_x0
        ! - 2D case 
        ELSE IF(c_dim == 2_idp) THEN
          ! - Init FFT
          alloc_local = FFTW_MPI_LOCAL_SIZE_2D(nz_group, nx_group/2+1,                  &
          MPI_COMM_GROUP_ID(i), local_nz, local_z0)

          ! - Sanity check: if local_nz<nprocz in group, this can lead to
          ! - local_nz=0 --> Abort
          IF(local_nz .EQ. 0_idp ) THEN
            WRITE(0,*) 'ERROR local_nz = 0 in rank ',rank
            CALL MPI_ABORT(comm,errcode,ierr)
          ENDIF
          local_ny = ny_group
          local_y0 = 0_idp
          local_nx = 2_idp*(nx_group/2_idp+1_idp)
          local_x0 = 0_idp
          ! - Init starting indexes/ local sizes of FFT output arrays along X,Y,Z
          ! - NB: At present, 2D geometry does not support fftw_transpose mode 
          ! - So dimensions and starting indexes of local input/output FFT arrays are 
          ! - identical 
          local_nz_tr=local_nz
          local_z0_tr=local_z0
          local_ny_tr=local_ny
          local_y0_tr=local_y0
          local_nx_tr=local_nx/2_idp
          local_x0_tr=local_x0  
        ENDIF
      ! - Case 2: P3DFFT is used for global FFT transpositions among groups
      ELSE
#if defined(P3DFFT)
       pdims(1) = INT(nprocy/nb_group_y,isp)
       pdims(2) = INT(nprocz/nb_group_z,isp)
       ! Set up P3DFFT plans and decomp
       is_in_place = .TRUE.
       CALL p3dfft_setup(pdims,INT(nx_group,isp),INT(ny_group,isp),INT(nz_group,isp),&
          mpi_comm_group_id(i),INT(nx_group,isp),INT(ny_group,isp),INT(nz_group,isp),&
          is_in_place)
       ! - Get local dimensions/starting indices of FFT arrays in real space
       CALL p3dfft_get_dims(p3d_istart,p3d_iend,p3d_isize,1_isp)
       ! - Get local dimensions/starting indices of FFT arrays in Fourier space
       CALL p3dfft_get_dims(p3d_fstart,p3d_fend,p3d_fsize,2_isp)
       ! - Init starting indexes/ local sizes of FFT input arrays along X,Y,Z
       local_nz = p3d_isize(3) ! Local size of FFT array along Z 
       local_ny = p3d_isize(2) ! Local size of FFT array along Y
       local_nx = p3d_isize(1) ! Local size of FFT array along X
       local_z0 = p3d_istart(3) - 1_idp ! Min global X-index boundary of local FFT array 
       local_y0 = p3d_istart(2) - 1_idp ! Min global Y-index boundary of local FFT array 
       local_x0 = p3d_istart(1) - 1_idp ! Min global Z-index boundary of local FFT array 
       ! - Init starting indexes/ local sizes of FFT output array along X,Y,Z
       ! - NB: At present, p3dfft does not support fftw_transpose mode 
       ! - So dimensions and starting indexes of local input/output FFT arrays are 
       ! - identical 
       local_nz_tr=p3d_fsize(3)
       local_ny_tr=p3d_fsize(2)
       local_nx_tr=p3d_fsize(1)
       local_z0_tr=p3d_fstart(3) - 1_idp
       local_y0_tr=p3d_fstart(3) - 2_idp
       local_x0_tr=p3d_fstart(3) - 3_idp
#endif
      ENDIF
    ENDIF
  ENDDO

  ! - Sanity check on the number of guard cells compared to the global size of the domain
  ! - In z and y directions
  IF(nz_global .LE. 2*nzg_group) THEN
    WRITE(*,*) '**ERROR **, nz_global too small compared to nzg_group'
    CALL MPI_ABORT(comm,errcode,ierr)
  ELSE IF (ny_global .LE. 2*nyg_group) THEN
    WRITE(*,*) '**ERROR **, ny_global too small compared to nyg_group'
    CALL MPI_ABORT(comm,errcode,ierr)
  ENDIF

  ! -- Gets min and max indices of the global FFT array without
  ! -- the group guard cells
  ! - Along Z 
  iz_min_r = 1_idp
  iz_max_r = local_nz
  IF(group_z_min_boundary) iz_min_r = nzg_group+1_idp
  IF(group_z_max_boundary) iz_max_r = local_nz-nzg_group
  ! - Along Y 
  iy_min_r = 1_idp
  iy_max_r = local_nz
  IF(group_y_min_boundary) iy_min_r = nyg_group+1_idp
  IF(group_y_max_boundary) iy_max_r = local_ny-nyg_group


  ! -- Computes global index of the lower z-boundary of the distributed FFT array
  ! -- (The one that includes group guard cells).
  ! -- 0 index corresponds to origin of the global simulation grid (WITHOUT guard cells)
  j=0
  DO i = 1 ,z_group_coords
     j = j +  nz_group_global_array(x_group_coords+                                  &
     y_group_coords*nb_group_x+(i-1_idp)*nb_group_x*nb_group_y+1_idp)
  ENDDO
  iz_min_global = local_z0-nzg_group + j !
  ! -- Computes global index of the upper z-boundary of the distributed FFT array
  ! -- (The one that includes group guard cells) 
  iz_max_global = iz_min_global + local_nz - 1_idp
  

  ! -- Array allocation for performing an MPI allgather operation
  ALLOCATE(cell_z_min_g(nprocz),cell_z_max_g(nprocz))
  ALLOCATE(all_iz_min_global(nproc),all_iz_max_global(nproc))

  ! - Gather all min indexes along z from all other procs
  CALL MPI_ALLGATHER(iz_min_global,1_isp, MPI_INTEGER8, all_iz_min_global,           &
  INT(1,isp), MPI_INTEGER8, mpi_ordered_comm_world, errcode)

  ! - Gather all max indexes along z from all other procs
  CALL MPI_ALLGATHER(iz_max_global,1_isp, MPI_INTEGER8, all_iz_max_global,           &
  INT(1,isp), MPI_INTEGER8, mpi_ordered_comm_world, errcode)

  ! -- Store min and max indices along z in 1D arrays
  DO i=1, nprocz
    cell_z_min_g(i) = all_iz_min_global(x_coords+y_coords*nprocx+(i-1)*nprocx*nprocy+1)
    cell_z_max_g(i) = all_iz_max_global(x_coords+y_coords*nprocx+(i-1)*nprocx*nprocy+1)
  ENDDO
  DEALLOCATE(all_iz_max_global,all_iz_min_global)

  ! -- Computes global index of the lower y-boundary of the distributed FFT array
  ! -- (The one that includes group guard cells). 
  ! -- 0 index corresponds to origin of the global simulation grid (WITHOUT guard cells) 
  j=0_idp
  DO i = 1_idp,y_group_coords
     j = j +ny_group_global_array(x_group_coords+(i-1_idp)*nb_group_x+               & 
     z_group_coords*nb_group_x*nb_group_y+1_idp)
  ENDDO
  iy_min_global = local_y0-nyg_group + j
  ! -- Computes global index of the upper y-boundary of the distributed FFT array
  ! -- (The one that includes group guard cells) 
  iy_max_global = iy_min_global + local_ny - 1_idp

  ALLOCATE(cell_y_min_g(nprocy),cell_y_max_g(nprocy))
  ALLOCATE(all_iy_min_global(nproc),all_iy_max_global(nproc))

  ! - Gather all min indexes along y from all other procs
  CALL MPI_ALLGATHER(iy_min_global,1_isp, MPI_INTEGER8, all_iy_min_global,&
  INT(1,isp), MPI_INTEGER8, mpi_ordered_comm_world, errcode)

  ! - Gather all max indexes along y from all other procs
  CALL MPI_ALLGATHER(iy_max_global,1_isp, MPI_INTEGER8, all_iy_max_global,&
  INT(1,isp), MPI_INTEGER8, mpi_ordered_comm_world, errcode)

  ! -- Store min and max indices along y in 1D arrays
  IF (c_dim==3) THEN
  DO i=1, nprocy
    cell_y_min_g(i) = all_iy_min_global(x_coords+(i-1)*nprocx+ & 
    z_coords*nprocx*nprocy+1_idp)
    cell_y_max_g(i) = all_iy_max_global(x_coords+(i-1)*nprocx+ &
    z_coords*nprocx*nprocy+1_idp)
  ENDDO
  ELSE
    cell_y_min_g(1)=0
    cell_y_max_g(1)=0 
  ENDIF
  
  DEALLOCATE(all_iy_max_global,all_iy_min_global)
  
  ! -- upper/lower x-boundaries of group 
  ix_min_r = 1_idp
  ! -- upper/lower x-boundaries of group
  ix_max_r = nx + 2*nxg_group

  ! -- Deallocate used arrays
  DEALLOCATE(grp_comm )
  !> Checks  if current group on domain boundary or not
  is_group_x_boundary_min = .FALSE.
  is_group_x_boundary_max = .FALSE.
  is_group_y_boundary_min = .FALSE.
  is_group_y_boundary_max = .FALSE.
  is_group_z_boundary_min = .FALSE.
  is_group_z_boundary_max = .FALSE.

  IF(x_group_coords==0_idp) is_group_x_boundary_min = .TRUE.
  IF(x_group_coords==nb_group_x-1_idp)  is_group_x_boundary_max = .TRUE.
  IF(y_group_coords==0_idp) is_group_y_boundary_min = .TRUE.
  IF(y_group_coords==nb_group_y-1_idp)  is_group_y_boundary_max = .TRUE.
  IF(z_group_coords==0_idp) is_group_z_boundary_min = .TRUE.
  IF(z_group_coords==nb_group_z-1_idp)  is_group_z_boundary_max = .TRUE.

  !> If current group is on one boundary and sorbing_bcs .TRUE. then current
  !> group contains a pml region
  IF(absorbing_bcs) THEN
    IF(is_group_x_boundary_min .OR. is_group_x_boundary_max .OR. &
       is_group_y_boundary_min .OR. is_group_y_boundary_max .OR. &
       is_group_z_boundary_min .OR. is_group_z_boundary_max) THEN
     ENDIF
  ENDIF
  ALLOCATE(cell_x_min_g(nprocx),cell_x_max_g(nprocx))
  cell_x_min_g = cell_x_min - nxguards
  cell_x_max_g = cell_x_max + nxguards
     

#if defined(DEBUG)
  WRITE(0, *) "setup_groups : end"
#endif

#endif
END SUBROUTINE setup_groups

! ______________________________________________________________________________________
!> @brief
!> This routine gathers on all ranks, the sizes of the distributed FFT array along z 
!> on each rank. This is useful to obtain the grid decomposition performed by the 
!> FFTW distributed FFT library - This is used when fftw_with_mpi=.TRUE. 
!> and fftw_hybrid=.FALSE.
!> @ author 
!> H. Kallala
!> 2018
! ______________________________________________________________________________________
SUBROUTINE adjust_grid_mpi_global
#if defined(FFTW)
USE iso_c_binding
USE mpi
USE mpi_fftw3, ONLY: local_nx, local_nx_tr, local_ny, local_ny_tr, local_nz,         &
  local_nz_tr, local_x0, local_x0_tr, local_y0, local_y0_tr, local_z0_tr
#endif
#if defined(FFTW)
USE picsar_precision, ONLY: idp, isp
USE shared_data, ONLY: cell_z_max, cell_z_min, comm, errcode, fftw_hybrid,           &
  fftw_mpi_transpose, fftw_with_mpi, ix_max_r, ix_min_r, iy_max_r, iy_min_r,         &
  iz_max_r, iz_min_r, nprocz, nx, nx_global, ny, ny_global, nz, nz_global,           &
  nz_global_grid_max, nz_global_grid_min, nz_grid, z_coords
#endif

#if defined(FFTW)
  INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: all_nz
  INTEGER(idp)  :: idim
  nz = local_nz
  nz_grid = nz + 1
  ALLOCATE(all_nz(1:nprocz))

  CALL MPI_ALLGATHER(nz, 1_isp, MPI_INTEGER8, all_nz, 1_isp, MPI_INTEGER8,    &
  comm, errcode)

  cell_z_min(1) = 0
  cell_z_max(1) = all_nz(1)-1

  DO idim=2, nprocz
    cell_z_min(idim) = cell_z_max(idim-1)+1
    cell_z_max(idim) = cell_z_min(idim) + all_nz(idim)-1
  ENDDO

  nz_global_grid_min = cell_z_min(z_coords+1)
  nz_global_grid_max = cell_z_max(z_coords+1)+1

  DEALLOCATE(all_nz)

  ix_min_r = 1
  ix_max_r = nx_global

  iy_min_r = 1
  iy_max_r = ny_global

  iz_min_r = 1
  iz_max_r = local_nz
  
  ! - Get local sizes of FFT input/output arrays along directions orthogonal to 
  ! - the direction of MPI CPU-split (FFTW-MPI only)
  IF ((.NOT. fftw_hybrid) .AND. (fftw_with_mpi)) THEN 
    IF (fftw_mpi_transpose) THEN 
      ! - Starting indexes/local sizes of FFT real arrays along Y
      local_ny=ny_global
      local_y0=0
      ! - Starting indexes/local sizes of FFT Fourier arrays along Y
      local_ny_tr=nz_global
      local_y0_tr=0
    ELSE
      ! - Starting indexes/local sizes of FFT input arrays along Y
      local_ny=ny_global
      local_y0=0
      ! - Starting indexes/local sizes of FFT output arrays along Y and Z
      local_ny_tr=local_ny
      local_y0_tr=local_y0
      local_nz_tr=local_nz
      local_z0_tr=0
    ENDIF 
    ! - Starting indexes/local sizes of FFT input arrays along X
    local_nx=(nx_global/2_idp+1)*2_idp
    local_x0=0
    ! - Starting indexes/local sizes of FFT output arrays along X
    local_nx_tr=local_nx/2_idp
    local_x0_tr=local_x0
  ENDIF 
#endif

END SUBROUTINE adjust_grid_mpi_global

! ______________________________________________________________________________________
!> @brief
!> This subroutine computes the space domain decomposition and
!> related parameters for each MPI process.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> H. Kallala
!
!> @date
!> Creation 2015
! ______________________________________________________________________________________
SUBROUTINE mpi_initialise
USE iso_c_binding
#if defined(FFTW)
USE load_balance
#endif
USE mpi
#if defined(FFTW)
USE mpi_fftw3, ONLY: alloc_local, fftw_mpi_local_size_3d,                            &
  fftw_mpi_local_size_3d_transposed, local_nz, local_nz_tr, local_z0, local_z0_tr
#endif
#if defined(FFTW)
USE picsar_precision, ONLY: idp, isp, num
USE shared_data, ONLY: absorbing_bcs, absorbing_bcs_x, absorbing_bcs_y,              &
  absorbing_bcs_z, c_dim, cell_x_max, cell_x_min, cell_y_max, cell_y_min,            &
  cell_z_max, cell_z_min, comm, dx, dy, dz, fftw_hybrid, fftw_mpi_transpose,         &
  fftw_with_mpi, length_x, length_x_part, length_y, length_y_part, length_z,         &
  length_z_part, nprocx, nprocy, nprocz, nx, nx_global, nx_global_grid_max,          &
  nx_global_grid_min, nx_grid, nxg_group, ny, ny_global, ny_global_grid_max,         &
  ny_global_grid_min, ny_grid, nyg_group, nz, nz_global, nz_global_grid_max,         &
  nz_global_grid_min, nz_grid, nzg_group, offset_grid_part_x_max,                    &
  offset_grid_part_x_min, offset_grid_part_y_max, offset_grid_part_y_min,            &
  offset_grid_part_z_max, offset_grid_part_z_min, p3dfft_flag, pbound_x_max,         &
  pbound_x_min, pbound_y_max, pbound_y_min, pbound_z_max, pbound_z_min, rank, x,     &
  x_coords, x_grid_max, x_grid_max_local, x_grid_maxs, x_grid_min, x_grid_min_local, &
  x_grid_mins, x_max_boundary, x_max_boundary_part, x_max_local, x_max_local_part,   &
  x_min_boundary, x_min_boundary_part, x_min_local, x_min_local_part, xmax,          &
  xmax_part, xmin, xmin_part, y, y_coords, y_grid_max, y_grid_max_local,             &
  y_grid_maxs, y_grid_min, y_grid_min_local, y_grid_mins, y_max_boundary,            &
  y_max_boundary_part, y_max_local, y_max_local_part, y_min_boundary,                &
  y_min_boundary_part, y_min_local, y_min_local_part, ymax, ymax_part, ymin,         &
  ymin_part, z, z_coords, z_grid_max, z_grid_max_local, z_grid_maxs, z_grid_min,     &
  z_grid_min_local, z_grid_mins, z_max_boundary, z_max_boundary_part, z_max_local,   &
  z_max_local_part, z_min_boundary, z_min_boundary_part, z_min_local,                &
  z_min_local_part, zmax, zmax_part, zmin, zmin_part
#endif
INTEGER(isp) :: idim
INTEGER(isp) :: nx0, nxp
INTEGER(isp) :: ny0, nyp
INTEGER(isp) :: nz0, nzp
#if defined(FFTW)
INTEGER(C_INTPTR_T):: kx, ly, mz
#endif

#if defined(DEBUG)
  WRITE(0, *) "mpi_initialize : start"
#endif

! Init number of guard cells of subdomains in each dimension
IF (l_smooth_compensate) THEN
  nxguards = nxguards + 1
  nyguards = nyguards + 1
  nzguards = nzguards + 1
END IF

! Init boundary conditions
IF(absorbing_bcs_x .OR. absorbing_bcs_y .OR. absorbing_bcs_z) absorbing_bcs = .TRUE.

CALL setup_communicator
ALLOCATE(x_grid_mins(1:nprocx), x_grid_maxs(1:nprocx))
ALLOCATE(y_grid_mins(1:nprocy), y_grid_maxs(1:nprocy))
ALLOCATE(z_grid_mins(1:nprocz), z_grid_maxs(1:nprocz))
ALLOCATE(cell_x_min(1:nprocx), cell_x_max(1:nprocx))
ALLOCATE(cell_y_min(1:nprocy), cell_y_max(1:nprocy))
ALLOCATE(cell_z_min(1:nprocz), cell_z_max(1:nprocz))

#if defined(FFTW)
IF(fftw_hybrid) THEN
  fftw_with_mpi = .TRUE.
  ! - Properly sets group guard cells along X 
  ! - As there is no group along x, impose nxg_group=nxguards 
  ! - (Only local FFTs along X at present)
  nxg_group=nxguards 
  ! - Properly sets group guard cells along Y
  IF (nyg_group .EQ. 0_idp) THEN
    nyg_group=nyguards ! Case when nyg_group has not been defined: use nyguards
  ENDIF
  ! - Properly sets group guard cells along Z 
  IF (nzg_group .EQ. 0_idp) THEN
    nzg_group=nzguards ! Case when nzg_group has not been defined: use nzguards
  ! - Properly sets group guard cells along X
  ! - As there is no groups along X, nxg_group and nxguards must be identical
  ENDIF
ENDIF
IF(p3dfft_flag) THEN
  fftw_with_mpi = .TRUE.
  fftw_hybrid = .TRUE.
  fftw_mpi_transpose = .FALSE.
ENDIF
! With fftw_with_mpi CPU split is performed along z only
IF (fftw_with_mpi .AND. .NOT. fftw_hybrid) THEN
  mz=INT(nz_global,C_INTPTR_T); ly=INT(ny_global,C_INTPTR_T); kx=INT(nx_global,C_INTPTR_T)
  !   get local data size and allocate (note dimension reversal)
  IF(.NOT. fftw_mpi_transpose) THEN
    alloc_local = fftw_mpi_local_size_3d(mz, ly, kx/2+1, comm, local_nz, local_z0)
  ELSE
    alloc_local = fftw_mpi_local_size_3d_transposed(mz, ly, kx/2+1, comm, local_nz,   &
    local_z0, local_nz_tr, local_z0_tr)
  ENDIF
  ! Regular CPU split
  ! If the total number of gridpoints cannot be exactly subdivided then fix
  ! The first nxp processors have nx0 cells
  ! The remaining processors have nx0+1 cells
  nx0 = nx_global / nprocx
  ny0 = ny_global / nprocy
  nz0 = nz_global /nprocz
ELSE
#endif
  ! Split is done on the total number of cells as in WARP
  ! Initial WARP split is used with each processor boundary
  ! being shared by two adjacent MPI processes
  nx0 = nx_global / nprocx
  ny0 = ny_global / nprocy
  nz0 = nz_global / nprocz
#if defined(FFTW)
ENDIF
#endif
IF (nx0 * nprocx .NE. nx_global) THEN
  nxp = (nx0 + 1) * nprocx - nx_global
ELSE
  nxp = nprocx
ENDIF

IF (ny0 * nprocy .NE. ny_global) THEN
  nyp = (ny0 + 1) * nprocy - ny_global
ELSE
  nyp = nprocy
ENDIF

IF (nz0 * nprocz .NE. nz_global) THEN
  nzp = (nz0 + 1) * nprocz - nz_global
ELSE
  nzp = nprocz
ENDIF

cell_x_min(1)=0
cell_x_max(1)=nx0-1
DO idim = 2, nxp
  cell_x_min(idim) = cell_x_max(idim-1)+1
  cell_x_max(idim) = cell_x_min(idim)+nx0-1
ENDDO
DO idim = nxp+1, nprocx
  cell_x_min(idim) = cell_x_max(idim-1)+1
  cell_x_max(idim) = cell_x_min(idim)+nx0
ENDDO

cell_y_min(1)=0
cell_y_max(1)=ny0-1
DO idim = 2, nyp
  cell_y_min(idim) = cell_y_max(idim-1)+1
  cell_y_max(idim) = cell_y_min(idim)+ny0-1
ENDDO
DO idim = nyp+1, nprocy
  cell_y_min(idim) = cell_y_max(idim-1)+1
  cell_y_max(idim) = cell_y_min(idim)+ny0
ENDDO

cell_z_min(1)=0
cell_z_max(1)=nz0-1
DO idim = 2, nzp
  cell_z_min(idim) = cell_z_max(idim-1)+1
  cell_z_max(idim) = cell_z_min(idim)+nz0-1
ENDDO
DO idim = nzp+1, nprocz
  cell_z_min(idim) = cell_z_max(idim-1)+1
  cell_z_max(idim) = cell_z_min(idim)+nz0
ENDDO

nx_global_grid_min = cell_x_min(x_coords+1)
nx_global_grid_max = cell_x_max(x_coords+1)+1

ny_global_grid_min = cell_y_min(y_coords+1)
ny_global_grid_max = cell_y_max(y_coords+1)+1

nz_global_grid_min = cell_z_min(z_coords+1)
nz_global_grid_max = cell_z_max(z_coords+1)+1


!!! --- number of gridpoints of each subdomain
nx_grid = nx_global_grid_max - nx_global_grid_min + 1
ny_grid = ny_global_grid_max - ny_global_grid_min + 1
nz_grid = nz_global_grid_max - nz_global_grid_min + 1

!!! --- number of cells of each subdomain
nx=nx_grid-1
ny=ny_grid-1
nz=nz_grid-1

!!! --- Set up global grid limits
length_x = xmax - xmin
dx = length_x / REAL(nx_global, num)
x_grid_min = xmin
x_grid_max = xmax

length_y = ymax - ymin
dy = length_y / REAL(ny_global, num)
y_grid_min = ymin
y_grid_max = ymax

length_z = zmax - zmin
dz = length_z / REAL(nz_global, num)
z_grid_min = zmin
z_grid_max = zmax

IF(c_dim == 2) THEN
  ny = 1_idp
  ny_global = 1_idp
!  cell_y_min(y_coords+1) = 1
!  cell_y_max(y_coords+1) = 1
  ny_grid = 2_idp
  nyguards = 0_idp
  nyjguards = 0_idp
ENDIF

#if defined(FFTW)
! -- Case of distributed FFT
IF(fftw_with_mpi) THEN
  ! -- Case of distributed FFT per MPI group only
  IF(fftw_hybrid) THEN
    CALL setup_groups
  ! -- Case of totally global FFT
  ELSE
    ! -- adjust nz to be equal to local_nz (since
    !!! --- it is computed differently by distributed FFT libraries)
    CALL adjust_grid_mpi_global
    IF(nz .NE. cell_z_max(z_coords+1) - cell_z_min(z_coords+1)+1) THEN
      WRITE(0, *) 'ERROR IN AJUSTING THE GRID 1'
      STOP
    ENDIF
  ENDIF
ENDIF

! -- When using distributed FFT, set
IF(fftw_hybrid) THEN
    ! Computes intersection of FFT distributed arrays and regular grid arrays
    CALL get2D_intersection_group_mpi()
    ! Load balancing is not used ,grid decomposition is that of
    ! fftw_local_size_3d
ENDIF
#endif
IF(absorbing_bcs .AND. l_spectral) THEN
  g_spectral = .TRUE. ! absorbing_bcs push only available with mult_mat_vec
                        ! routine
ENDIF
IF(absorbing_bcs .AND. .NOT. l_spectral) THEN
 IF(rank==0)  WRITE(0, *)'ERROR , pmls are not available yet with FDTD'
 STOP
ENDIF

!!! --- Set up global grid limits

CALL compute_simulation_axis()

x_min_local = x_grid_mins(x_coords+1)
x_max_local = x_grid_maxs(x_coords+1)
y_min_local = y_grid_mins(y_coords+1)
y_max_local = y_grid_maxs(y_coords+1)
z_min_local = z_grid_mins(z_coords+1)
z_max_local = z_grid_maxs(z_coords+1)

x_grid_min_local=x_min_local
y_grid_min_local=y_min_local
z_grid_min_local=z_min_local
x_grid_max_local=x_max_local
y_grid_max_local=y_max_local
z_grid_max_local=z_max_local

CALL allocate_grid_quantities()
start_time = MPI_WTIME()

! ---- set up global particle domain boundaries
! ----- Set up particle domain external boundaries flags
! ----- Set up local domain boundaries
! -- Xmin
IF  ((pbound_x_min .EQ. 3) .OR. (pbound_x_min .EQ. 1)) THEN
  xmin_part=xmin+offset_grid_part_x_min
  IF ((xmin_part .GE. x_min_local) .AND. (xmin_part .LT. x_max_local)) THEN
    x_min_boundary_part= .TRUE.
    x_min_local_part = xmin_part
  ELSE
    x_min_boundary_part= .FALSE.
    x_min_local_part = x_min_local
  ENDIF
ELSE
  xmin_part=xmin
  x_min_boundary_part=x_min_boundary
  x_min_local_part=x_min_local
ENDIF
! -- Xmax
IF  ((pbound_x_max .EQ. 3) .OR. (pbound_x_max .EQ. 1)) THEN
  xmax_part=xmax+offset_grid_part_x_max
  IF ((xmax_part .GE. x_min_local) .AND. (xmax_part .LT. x_max_local)) THEN
    x_max_boundary_part= .TRUE.
    x_max_local_part = xmax_part
  ELSE
    x_max_boundary_part= .FALSE.
    x_max_local_part = x_max_local
  ENDIF
ELSE
  xmax_part=xmax
  x_max_boundary_part=x_max_boundary
  x_max_local_part=x_max_local
ENDIF
! -- Ymin
IF  ((pbound_y_min .EQ. 3) .OR. (pbound_y_min .EQ. 1)) THEN
  ymin_part=ymin+offset_grid_part_y_min
  IF ((ymin_part .GE. y_min_local) .AND. (ymin_part .LT. y_max_local)) THEN
    y_min_boundary_part= .TRUE.
    y_min_local_part = ymin_part
  ELSE
    y_min_boundary_part= .FALSE.
    y_min_local_part = y_min_local
  ENDIF
ELSE
  ymin_part=ymin
  y_min_boundary_part=y_min_boundary
  y_min_local_part=y_min_local
ENDIF
! -- Ymax
IF  ((pbound_y_max .EQ. 3) .OR. (pbound_y_max .EQ. 1)) THEN
  ymax_part=ymax+offset_grid_part_y_max
  IF ((ymax_part .GE. y_min_local) .AND. (ymax_part .LT. y_max_local)) THEN
    y_max_boundary_part= .TRUE.
    y_max_local_part = ymax_part
  ELSE
    y_max_boundary_part= .FALSE.
    y_max_local_part = y_max_local
  ENDIF
ELSE
  ymax_part=ymax
  y_max_boundary_part=y_max_boundary
  y_max_local_part=y_max_local
ENDIF
! -- zmin
IF  ((pbound_z_min .EQ. 3) .OR. (pbound_z_min .EQ. 1)) THEN
  zmin_part=zmin+offset_grid_part_z_min
  IF ((zmin_part .GE. z_min_local) .AND. (zmin_part .LT. z_max_local)) THEN
    z_min_boundary_part= .TRUE.
    z_min_local_part = zmin_part
  ELSE
    z_min_boundary_part= .FALSE.
    z_min_local_part = z_min_local
  ENDIF
ELSE
  zmin_part=zmin
  z_min_boundary_part=z_min_boundary
  z_min_local_part=z_min_local
ENDIF
! -- zmax
IF  ((pbound_z_max .EQ. 3) .OR. (pbound_z_max .EQ. 1)) THEN
  zmax_part=zmax+offset_grid_part_z_max
  IF ((zmax_part .GE. z_min_local) .AND. (zmax_part .LT. z_max_local)) THEN
    z_max_boundary_part= .TRUE.
    z_max_local_part = zmax_part
  ELSE
    z_max_boundary_part= .FALSE.
    z_max_local_part = z_max_local
  ENDIF
ELSE
  zmax_part=zmax
  z_max_boundary_part=z_max_boundary
  z_max_local_part=z_max_local
ENDIF

! Set particle domain extent
length_x_part = xmax_part - xmin_part
length_y_part = ymax_part - ymin_part
length_z_part = zmax_part - zmin_part

#if defined(DEBUG)
  WRITE(0, *) "mpi_initialize : end"
#endif
END SUBROUTINE mpi_initialise

! ______________________________________________________________________________________
!> @brief
!> This subroutine allocates and computes the MPI process local and global
!> grid minima and maxima.
!
!> This subroutine is called in mpi_initialise().
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ______________________________________________________________________________________
SUBROUTINE compute_simulation_axis()
IMPLICIT NONE
INTEGER(idp) :: ix, iy, iz, iproc

IF (.NOT. l_axis_allocated) THEN
  ! Allocate arrays of axis
  ALLOCATE(x(-nxguards:nx+nxguards))
  ALLOCATE(y(-nyguards:ny+nyguards))
  ALLOCATE(z(-nzguards:nz+nzguards))
  ALLOCATE(x_global(-nxguards:nx_global+nxguards))
  ALLOCATE(y_global(-nyguards:ny_global+nyguards))
  ALLOCATE(z_global(-nzguards:nz_global+nzguards))
  l_axis_allocated=.TRUE.
ENDIF
!!! --- Set up global grid
DO ix = -nxguards, nx_global+nxguards
  x_global(ix) = x_grid_min + ix * dx
ENDDO
DO iy = -nyguards, ny_global+nyguards
  y_global(iy) = y_grid_min + iy * dy
ENDDO
DO iz = -nzguards, nz_global+nzguards
  z_global(iz) = z_grid_min + iz * dz
ENDDO

!!! --- Set up local grid maxima and minima
DO iproc = 1, nprocx
  x_grid_mins(iproc) = x_global(cell_x_min(iproc))
  x_grid_maxs(iproc) = x_global(cell_x_max(iproc)+1)
ENDDO
IF(c_dim == 3) THEN
  DO iproc = 1, nprocy
    y_grid_mins(iproc) = y_global(cell_y_min(iproc))
    y_grid_maxs(iproc) = y_global(cell_y_max(iproc)+1)
  ENDDO
ELSE
 y_grid_mins(1) = 0
 y_grid_maxs(1) = 0
ENDIF 
DO iproc = 1, nprocz
  z_grid_mins(iproc) = z_global(cell_z_min(iproc))
  z_grid_maxs(iproc) = z_global(cell_z_max(iproc)+1)
ENDDO

END SUBROUTINE compute_simulation_axis

! ______________________________________________________________________________________
!> @brief
!> This subroutine allocates grid quantities such as fields, currents, charge
!> and electric field divergence.
!
!> This subroutine is called at the end of mpi_initialise().
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ______________________________________________________________________________________
SUBROUTINE allocate_grid_quantities()
USE iso_c_binding
#if defined(FFTW)
USE mpi_fftw3, ONLY: alloc_local, fftw_alloc_complex, fftw_alloc_real, local_nx,     &
  local_nx_tr, local_ny, local_ny_tr, local_nz, local_nz_tr
USE picsar_precision, ONLY: idp
#endif
IMPLICIT NONE
#if defined(FFTW)
TYPE(C_PTR) :: cdata, cin
INTEGER(idp) :: imn, imx, jmn, jmx, kmn, kmx
INTEGER(idp) :: nxx, nyy, nzz
#endif
! --- Allocate regular grid quantities (in real space)
ALLOCATE(ex(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(ey(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(ez(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(bx(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(by(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(bz(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
! > When using absorbing_bcs , allocate splitted fields 
IF(absorbing_bcs) THEN
  ALLOCATE(exy(-nxguards:nx+nxguards, -nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(exz(-nxguards:nx+nxguards, -nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(eyx(-nxguards:nx+nxguards, -nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(eyz(-nxguards:nx+nxguards, -nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(ezx(-nxguards:nx+nxguards, -nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(ezy(-nxguards:nx+nxguards, -nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(bxy(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(bxz(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(byx(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(byz(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(bzx(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
  ALLOCATE(bzy(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
ENDIF
ALLOCATE(jx(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards,                     &
-nzjguards:nz+nzjguards))
ALLOCATE(jy(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards,                     &
-nzjguards:nz+nzjguards))
ALLOCATE(jz(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards,                     &
-nzjguards:nz+nzjguards))
ALLOCATE(rho(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards,                    &
-nzjguards:nz+nzjguards))
ALLOCATE(rhoold(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards,                 &
-nzjguards:nz+nzjguards))
ALLOCATE(dive(-nxguards:nx+nxguards, -nyguards:ny+nyguards,                       &
-nzguards:nz+nzguards))
ALLOCATE(divj(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(divb(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
! --- Initialize auxiliary field arrays for gather to particles
ex_p => ex
ey_p => ey
ez_p => ez
bx_p => bx
by_p => by
bz_p => bz


#if defined(FFTW)
! ---  Allocate grid quantities in Fourier space
IF (l_spectral) THEN

  ! - Case when fftw_with_mpi is .TRUE. (distributed FFT)
  IF (fftw_with_mpi) THEN
    ! - FFT arrays dimensions in Fourier space along X,Y,Z
    nkx=local_nx_tr
    nky=local_ny_tr
    nkz=local_nz_tr
    IF(c_dim == 2) THEN
      nky = 1_idp
    ENDIF
    ! - Allocate complex FFT arrays
    ! - Case when p3dfft_flag is .TRUE. (p3dfft is used for distributed FFT)
    IF(p3dfft_flag) THEN
      IF(.NOT. g_spectral) THEN
        ALLOCATE(exf(nkx,nky,nkz))
        ALLOCATE(eyf(nkx,nky,nkz))
        ALLOCATE(ezf(nkx,nky,nkz))
        ALLOCATE(bxf(nkx,nky,nkz))
        ALLOCATE(byf(nkx,nky,nkz))
        ALLOCATE(bzf(nkx,nky,nkz))
        ALLOCATE(jxf(nkx,nky,nkz))
        ALLOCATE(jyf(nkx,nky,nkz))
        ALLOCATE(jzf(nkx,nky,nkz))
        ALLOCATE(rhof(nkx,nky,nkz))
        ALLOCATE(rhooldf(nkx,nky,nkz))
      ENDIF
    ! - Case when FFTW is used for the distributed FFT
    ELSE IF(.NOT. p3dfft_flag) THEN
      IF(.NOT. g_spectral) THEN
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, exf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, eyf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, ezf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, bxf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, byf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, bzf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, jxf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, jyf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, jzf, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, rhof, [nkx, nky, nkz])
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, rhooldf, [nkx, nky, nkz])
      ENDIF
    ENDIF
    ! - Allocate real FFT arrays 
    nxx = local_nx
    nyy = local_ny
    nzz = local_nz
    IF(c_dim == 2) THEN
      nyy = 1_idp
    ENDIF
    ! - Case when p3dfft_flag is .TRUE. (p3dfft is used for distributed FFT)
    IF(.NOT. p3dfft_flag) THEN
    ! - When using absorbing_bcs, merged fields are not allocated in fourier space
    ! - neither ex_r,ey_r ... components
    ! - In this case only splitted fields are allocated  
    ! - The merge is done using local fields (ex = exy+exz )
     IF(.NOT. absorbing_bcs) THEN
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, ex_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, ey_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, ez_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, bx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, by_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, bz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, jx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, jy_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, jz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, rho_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, rhoold_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
     ELSE IF(absorbing_bcs) THEN
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, exy_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, exz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, eyx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, eyz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, ezx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, ezy_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, bxy_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, bxz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, byx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, byz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, bzx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, bzy_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, jx_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, jy_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, jz_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, rho_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
       CALL c_f_pointer(cin, rhoold_r, [nxx, nyy, nzz])
       cin = fftw_alloc_real(2 * alloc_local);
      ENDIF
    ! - Case when FFTW is used for the distributed FFT
    ELSE IF(p3dfft_flag) THEN
    ! - When using absorbing_bcs, merged fields are not allocated in fourier space
    ! - neither ex_r,ey_r ... components
    ! - In this case only splitted fields are allocated  
    ! - The merge is done using local fields (ex = exy+exz )
      IF(.NOT. absorbing_bcs) THEN
        ALLOCATE(ex_r(nxx,nyy,nzz))
        ALLOCATE(ey_r(nxx,nyy,nzz))
        ALLOCATE(ez_r(nxx,nyy,nzz))
        ALLOCATE(bx_r(nxx,nyy,nzz))
        ALLOCATE(by_r(nxx,nyy,nzz))
        ALLOCATE(bz_r(nxx,nyy,nzz))
        ALLOCATE(jx_r(nxx,nyy,nzz))
        ALLOCATE(jy_r(nxx,nyy,nzz))
        ALLOCATE(jz_r(nxx,nyy,nzz))
        ALLOCATE(rho_r(nxx,nyy,nzz))
        ALLOCATE(rhoold_r(nxx,nyy,nzz))
      ELSE IF(absorbing_bcs) THEN
        ALLOCATE(exy_r(nxx,nyy,nzz))
        ALLOCATE(exz_r(nxx,nyy,nzz))
        ALLOCATE(eyx_r(nxx,nyy,nzz))
        ALLOCATE(eyz_r(nxx,nyy,nzz))
        ALLOCATE(ezx_r(nxx,nyy,nzz))
        ALLOCATE(ezy_r(nxx,nyy,nzz))
        ALLOCATE(bxy_r(nxx,nyy,nzz))
        ALLOCATE(bxz_r(nxx,nyy,nzz))
        ALLOCATE(byx_r(nxx,nyy,nzz))
        ALLOCATE(byz_r(nxx,nyy,nzz))
        ALLOCATE(bzx_r(nxx,nyy,nzz))
        ALLOCATE(bzy_r(nxx,nyy,nzz))
        ALLOCATE(jx_r(nxx,nyy,nzz))
        ALLOCATE(jy_r(nxx,nyy,nzz))
        ALLOCATE(jz_r(nxx,nyy,nzz))
        ALLOCATE(rho_r(nxx,nyy,nzz))
        ALLOCATE(rhoold_r(nxx,nyy,nzz))
      ENDIF
    ENDIF      
  ! Case of local FFTs (purely local pseudo-spectral solver)
  ELSE IF(.NOT. fftw_with_mpi) THEN
    nkx=(2*nxguards+nx)/2+1! Real To Complex Transform
    nky=(2*nyguards+ny)
    nkz=(2*nzguards+nz)
    IF(c_dim == 2) THEN
      nky =1_idp
    ENDIF
    IF(.NOT. g_spectral) THEN
    ! - Allocate complex FFT arrays
      ALLOCATE(exf(nkx, nky, nkz))
      ALLOCATE(eyf(nkx, nky, nkz))
      ALLOCATE(ezf(nkx, nky, nkz))
      ALLOCATE(bxf(nkx, nky, nkz))
      ALLOCATE(byf(nkx, nky, nkz))
      ALLOCATE(bzf(nkx, nky, nkz))
      ALLOCATE(jxf(nkx, nky, nkz))
      ALLOCATE(jyf(nkx, nky, nkz))
      ALLOCATE(jzf(nkx, nky, nkz))
      ALLOCATE(rhof(nkx, nky, nkz))
      ALLOCATE(rhooldf(nkx, nky, nkz))
    ENDIF
    ! - Allocate real FFT arrays 
    imn=-nxguards; imx=nx+nxguards-1
    jmn=-nyguards;jmx=ny+nyguards-1
    kmn=-nzguards;kmx=nz+nzguards-1
    IF(c_dim == 2) THEN
      jmn = 0
      jmx = 0
    ENDIF
    IF (.NOT. absorbing_bcs) THEN
    ! - When using absorbing_bcs, merged fields are not allocated in fourier space
    ! - neither ex_r,ey_r ... components
    ! - In this case only splitted fields are allocated  
    ! - The merge is done using local fields (ex = exy+exz )
      ALLOCATE(ex_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(ey_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(ez_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(bx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(by_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(bz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(jx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(jy_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(jz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(rho_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(rhoold_r(imn:imx, jmn:jmx, kmn:kmx))
    ELSE IF(absorbing_bcs) THEN
      ALLOCATE(exy_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(exz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(eyx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(eyz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(ezx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(ezy_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(bxy_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(bxz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(byx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(byz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(bzx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(bzy_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(jx_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(jy_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(jz_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(rho_r(imn:imx, jmn:jmx, kmn:kmx))
      ALLOCATE(rhoold_r(imn:imx, jmn:jmx, kmn:kmx))
    ENDIF
  ENDIF
ENDIF
#endif

! --- Quantities used by the dynamic load balancer for distributed FFTs
ALLOCATE(new_cell_x_min(1:nprocx), new_cell_x_max(1:nprocx))
ALLOCATE(new_cell_y_min(1:nprocy), new_cell_y_max(1:nprocy))
ALLOCATE(new_cell_z_min(1:nprocz), new_cell_z_max(1:nprocz))

END SUBROUTINE allocate_grid_quantities

! ______________________________________________________________________________________
!> This subroutine finalizes MPI with some time information.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ______________________________________________________________________________________
SUBROUTINE mpi_close
USE mpi
INTEGER :: seconds, minutes, hours, total

IF (rank .EQ. 0) THEN
  end_time = MPI_WTIME()
  total = INT(end_time - start_time)
  seconds = MOD(total, 60)
  minutes = MOD(total / 60, 60)
  hours = total / 3600
ENDIF

CALL MPI_FINALIZE(errcode)

END SUBROUTINE mpi_close

! ______________________________________________________________________________________
!> @brief
!> Subroutine dedicated to the time report at the end of the simulation.
!> This subroutines gathers time statistics from all processes and print a summary
!> of the time spent in each important step of the main PIC loop.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
! ______________________________________________________________________________________
SUBROUTINE time_statistics
#ifdef _OPENMP
USE omp_lib
#endif
USE params, ONLY: currdepo, fg_p_pp_separated, fieldgathe, it, nsteps, rhodepo
USE picsar_precision, ONLY: idp, isp, num
USE time_stat, ONLY: avetimes, init_avetimes, init_localtimes, init_maxtimes,        &
  init_mintimes, localtimes, maxtimes, mintimes
IMPLICIT NONE
REAL(num), DIMENSION(26) :: percenttimes
INTEGER(idp)             :: nthreads_tot

! Total times
localtimes(20) = sum(localtimes(1:14)) + localtimes(24)
localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + localtimes(8) +      &
localtimes(11) + localtimes(13)
localtimes(24)  = localtimes(24)  + localtimes(26) 
localtimes(18) = localtimes(5) + localtimes(6) + localtimes(7) + localtimes(8) +      &
localtimes(24)

init_localtimes(5) = sum(init_localtimes(1:4))

! Reductions
! Maximun times
CALL MPI_REDUCE(localtimes, maxtimes, 26_isp, mpidbl, MPI_MAX, 0_isp, comm, errcode)
CALL MPI_REDUCE(init_localtimes, init_maxtimes, 5_isp, mpidbl, MPI_MAX, 0_isp, comm,  &
errcode)
! Minimum times
CALL MPI_REDUCE(localtimes, mintimes, 26_isp, mpidbl, MPI_MIN, 0_isp, comm, errcode)
CALL MPI_REDUCE(init_localtimes, init_mintimes, 5_isp, mpidbl, MPI_MIN, 0_isp, comm,  &
errcode)
! Average
CALL MPI_REDUCE(localtimes, avetimes, 26_isp, mpidbl, MPI_SUM, 0_isp, comm, errcode)
CALL MPI_REDUCE(init_localtimes, init_avetimes, 5_isp, mpidbl, MPI_SUM, 0_isp, comm,  &
errcode)
avetimes = avetimes / nproc
init_avetimes = init_avetimes / nproc

! Percentage
percenttimes = avetimes / avetimes(20) * 100

IF (rank .EQ. 0) THEN
  WRITE(0, *)                                                                         &
  '___________________________________________________________________________'
  WRITE(0, '(X, A40, X)') "Time statistics Initialization:"
  WRITE(0, *) ""
  WRITE(0, '(X, A25, X, A8, X, A8, X, A8, X, A8, X, A8)') "Step part", "min (s)",     &
  "ave (s)", "max (s)"
  WRITE(0, *)                                                                         &
  "---------------------------------------------------------------------------"
  WRITE(0, '(X, A25, 5(X, F8.2))') "Tiling and part. load:", init_mintimes(1),        &
  init_avetimes(1), init_maxtimes(1)
  WRITE(0, *) ""
  WRITE(0, *)                                                                         &
  '___________________________________________________________________________'
  WRITE(0, '(X, A40, X)') "Time statistics main loop:"
  WRITE(0, *) ""
  WRITE(0, '(X, A25, X, A8, X, A8, X, A8, X, A8, X, A8)') "Step part", "min (s)",     &
  "ave (s)", "max (s)", "per (%)", "/it (ms)"
  WRITE(0, *)                                                                         &
  "---------------------------------------------------------------------------"
  IF (fg_p_pp_separated.le.1) THEN
    WRITE(0, '(X, A25, 5(X, F8.2))') "Particle pusher + field gathering:",            &
    mintimes(1), avetimes(1), maxtimes(1), percenttimes(1), avetimes(1)/nsteps*1e3
  ELSE
    WRITE(0, '(X, A25, 5(X, F8.2))') "Field gathering:", mintimes(14), avetimes(14),  &
    maxtimes(14), percenttimes(14), avetimes(14)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Particle pusher:", mintimes(1), avetimes(1),    &
    maxtimes(1), percenttimes(1), avetimes(1)/nsteps*1e3
  ENDIF
  WRITE(0, '(X, A25, 5(X, F8.2))') "Particle MPI bound. cond.:", mintimes(2),         &
  avetimes(2), maxtimes(2), percenttimes(2), avetimes(2)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Particle OpenMP bound. cond.:", mintimes(11),     &
  avetimes(11), maxtimes(11), percenttimes(11), avetimes(11)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Current deposition:", mintimes(3), avetimes(3),   &
  maxtimes(3), percenttimes(3), avetimes(3)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Current bound. cond.:", mintimes(4), avetimes(4), &
  maxtimes(4), percenttimes(4), avetimes(4)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Push bfield fdtd:", mintimes(5), avetimes(5),     &
  maxtimes(5), percenttimes(5), avetimes(5)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "B field bound. cond.:", mintimes(6), avetimes(6), &
  maxtimes(6), percenttimes(6), avetimes(6)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Push efield fdtd:", mintimes(7), avetimes(7),     &
  maxtimes(7), percenttimes(7), avetimes(7)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Get fields psatd:", mintimes(21), avetimes(21),   &
  maxtimes(21), percenttimes(21), avetimes(21)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Execute fftw:", mintimes(22), avetimes(22),       &
  maxtimes(22), percenttimes(22), avetimes(22)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Multiply f_fields:", mintimes(23), avetimes(23),  &
  maxtimes(23), percenttimes(23), avetimes(23)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Absorbing_bcs:", mintimes(26),avetimes(26),  &
  maxtimes(26), percenttimes(23), avetimes(26)/nsteps*1e3

  WRITE(0, '(X, A25, 5(X, F8.2))') "Total solve max psatd:", mintimes(24),            &
  avetimes(24), maxtimes(24), percenttimes(24), avetimes(24)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Group mpi comm:", mintimes(25), avetimes(25),     &
  maxtimes(25), percenttimes(25), avetimes(25)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "E field bound. cond.:", mintimes(8), avetimes(8), &
  maxtimes(8), percenttimes(8), avetimes(8)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Sorting:", mintimes(10), avetimes(10),            &
  maxtimes(10), percenttimes(10), avetimes(10)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Charge deposition:", mintimes(12), avetimes(12),  &
  maxtimes(12), percenttimes(12), avetimes(12)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Charge bound. cond.:", mintimes(13),              &
  avetimes(13), maxtimes(13), percenttimes(13), avetimes(13)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Diags:", mintimes(9), avetimes(9), maxtimes(9),   &
  percenttimes(9), avetimes(9)/nsteps*1e3
  WRITE(0, *) ""
  WRITE(0, '(X, A25, 5(X, F8.2))') "Total time Maxwell solver:", mintimes(18),        &
  avetimes(18), maxtimes(18), percenttimes(18), avetimes(18)/nsteps*1e3
  WRITE(0, '(X, A25, 5(X, F8.2))') "Total time bound. cond.:", mintimes(19),          &
  avetimes(19), maxtimes(19), percenttimes(19), avetimes(19)/nsteps*1e3
  WRITE(0, '(X, A25, X, F8.2, X, F8.2, X, F8.2)') "Total time:", mintimes(20),        &
  avetimes(20), maxtimes(20)


#ifdef _OPENMP
  nthreads_tot=OMP_GET_MAX_THREADS()
  CALL OMP_SET_NESTED(.TRUE.)
#else
  nthreads_tot=1
#endif

  WRITE(0, *) ''
  IF (fg_p_pp_separated.le.1) THEN
    WRITE(0, *) 'For lib_performance python class:'
    WRITE(0, '("(nmpi=", I5, ", nomp=", I5, ", name='''', kernel=", F6.2, ",          &
    fieldgathe=", F6.2, ", part_mpi_com=", F6.2, ", part_omp_com=", F6.2, ",          &
    currdepo=", F6.2, ", currcom=", F6.2, ", maxwell=", F6.2, ", maxwellcom=", F6.2,  &
    ", sorting=", F6.2, ", rhodepo=", F6.2, ", rhocom=", F6.2, ", diags=", F6.2,      &
    ")")') nproc, nthreads_tot, avetimes(20), avetimes(1), avetimes(2), avetimes(11), &
    avetimes(3), avetimes(4), avetimes(5)+avetimes(7), avetimes(6)+avetimes(8),       &
    avetimes(10), avetimes(12), avetimes(13), avetimes(9)
  ELSE
    WRITE(0, *) 'For lib_performance python class:'
    WRITE(0, '("(nmpi=", I5, ", nomp=", I5, ", name='''', kernel=", F6.2, ",          &
    fieldgathe=", F6.2, ", partpusher=", F6.2, ", part_mpi_com=", F6.2, ",            &
    part_omp_com=", F6.2, ", currdepo=", F6.2, ", currcom=", F6.2, ", maxwell=",      &
    F6.2, ", maxwellcom=", F6.2, ", sorting=", F6.2, ", rhodepo=", F6.2, ", rhocom=", &
    F6.2, ", diags=", F6.2, ")")') nproc, nthreads_tot, avetimes(20), avetimes(14),   &
    avetimes(1), avetimes(2), avetimes(11), avetimes(3), avetimes(4),                 &
    avetimes(5)+avetimes(7), avetimes(6)+avetimes(8), avetimes(10), avetimes(12),     &
    avetimes(13), avetimes(9)
  ENDIF

ENDIF

END SUBROUTINE time_statistics

! ____________________________________________________________________________
!> @brief
!> Subroutine dedicated to the time statistics for one iteration
!
!> This subroutine prints time statistics during the simulation
!> every iterations if timestat_perit > 0.
!>
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2016
! ______________________________________________________________________________________
SUBROUTINE time_statistics_per_iteration
#ifdef _OPENMP
USE omp_lib
#endif
USE params, ONLY: fg_p_pp_separated, it, nsteps
USE picsar_precision, ONLY: isp, num
USE time_stat, ONLY: avetimes, init_avetimes, init_localtimes, init_maxtimes,        &
  init_mintimes, localtimes, maxtimes, mintimes, timestat_perit
IMPLICIT NONE

REAL(num), DIMENSION(26) :: percenttimes

! Time stats per iteration activated
IF (timestat_perit.gt.0) THEN

  ! Total times
  localtimes(20) = sum(localtimes(1:14)) + localtimes(24)
  localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + localtimes(8) +    &
  localtimes(11) + localtimes(13)
  localtimes(24) = localtimes(24) + localtimes(26) 
  localtimes(18) = localtimes(5) + localtimes(6) + localtimes(7) + localtimes(8) +    &
  localtimes(24)

  init_localtimes(5) = sum(init_localtimes(1:4))

  ! Reductions
  ! Maximun times
  CALL MPI_REDUCE(localtimes, maxtimes, 26_isp, mpidbl, MPI_MAX, 0_isp, comm,         &
  errcode)
  CALL MPI_REDUCE(init_localtimes, init_maxtimes, 5_isp, mpidbl, MPI_MAX, 0_isp,      &
  comm, errcode)
  ! Minimum times
  CALL MPI_REDUCE(localtimes, mintimes, 26_isp, mpidbl, MPI_MIN, 0_isp, comm,         &
  errcode)
  CALL MPI_REDUCE(init_localtimes, init_mintimes, 5_isp, mpidbl, MPI_MAX, 0_isp,      &
  comm, errcode)
  ! Average
  CALL MPI_REDUCE(localtimes, avetimes, 26_isp, mpidbl, MPI_SUM, 0_isp, comm,         &
  errcode)
  CALL MPI_REDUCE(init_localtimes, init_avetimes, 5_isp, mpidbl, MPI_MAX, 0_isp,      &
  comm, errcode)
  avetimes = avetimes / nproc

  ! Percentage
  percenttimes = avetimes / avetimes(20) * 100

  IF (rank .EQ. 0) THEN
    WRITE(0, *)                                                                       &
    '___________________________________________________________________________'
    WRITE(0, '(X, A25, X, A8, X, A8, X, A8, X, A8, X, A8)') "Step part", "min (s)",   &
    "ave (s)", "max (s)", "per (%)", "/it (ms)"
    WRITE(0, *)                                                                       &
    "---------------------------------------------------------------------------"
    IF (fg_p_pp_separated.le.1) THEN
      WRITE(0, '(X, A25, 5(X, F8.2))') "Particle pusher + field gathering:",          &
      mintimes(1), avetimes(1), maxtimes(1), percenttimes(1), avetimes(1)/nsteps*1e3
    ELSE
      WRITE(0, '(X, A25, 5(X, F8.2))') "Field gathering:", mintimes(14),              &
      avetimes(14), maxtimes(14), percenttimes(14), avetimes(14)/nsteps*1e3
      WRITE(0, '(X, A25, 5(X, F8.2))') "Particle pusher:", mintimes(1), avetimes(1),  &
      maxtimes(1), percenttimes(1), avetimes(1)/nsteps*1e3
    ENDIF
    WRITE(0, '(X, A25, 5(X, F8.2))') "Particle MPI bound. cond.:", mintimes(2),       &
    avetimes(2), maxtimes(2), percenttimes(2), avetimes(2)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Particle OpenMP bound. cond.:", mintimes(11),   &
    avetimes(11), maxtimes(11), percenttimes(11), avetimes(11)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Current deposition:", mintimes(3), avetimes(3), &
    maxtimes(3), percenttimes(3), avetimes(3)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Current bound. cond.:", mintimes(4),            &
    avetimes(4), maxtimes(4), percenttimes(4), avetimes(4)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Push bfield fdtd:", mintimes(5), avetimes(5),   &
    maxtimes(5), percenttimes(5), avetimes(5)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "B field bound. cond.:", mintimes(6),            &
    avetimes(6), maxtimes(6), percenttimes(6), avetimes(6)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Push efield fdtd:", mintimes(7), avetimes(7),   &
    maxtimes(7), percenttimes(7), avetimes(7)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Get fields psatd:", mintimes(21), avetimes(21), &
    maxtimes(21), percenttimes(21), avetimes(21)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Execute fftw:", mintimes(22), avetimes(22),     &
    maxtimes(22), percenttimes(22), avetimes(22)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Multiply f_fields:", mintimes(23),              &
    avetimes(23), maxtimes(23), percenttimes(23), avetimes(23)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Absorbing_bcs:", mintimes(26),&
    avetimes(26), maxtimes(26), percenttimes(26), avetimes(26)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Total solve max psatd:", mintimes(24),          &
    avetimes(24), maxtimes(24), percenttimes(24), avetimes(24)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Group mpi comm:", mintimes(25), avetimes(25),   &
    maxtimes(25), percenttimes(25), avetimes(25)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "E field bound. cond.:", mintimes(8),            &
    avetimes(8), maxtimes(8), percenttimes(8), avetimes(8)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Sorting:", mintimes(10), avetimes(10),          &
    maxtimes(10), percenttimes(10), avetimes(10)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Charge deposition:", mintimes(12),              &
    avetimes(12), maxtimes(12), percenttimes(12), avetimes(12)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Charge bound. cond.:", mintimes(13),            &
    avetimes(13), maxtimes(13), percenttimes(13), avetimes(13)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Diags:", mintimes(9), avetimes(9), maxtimes(9), &
    percenttimes(9), avetimes(9)/nsteps*1e3
    WRITE(0, *) ""
    WRITE(0, '(X, A25, 5(X, F8.2))') "Total time Maxwell solver:", mintimes(18),      &
    avetimes(18), maxtimes(18), percenttimes(18), avetimes(18)/nsteps*1e3
    WRITE(0, '(X, A25, 5(X, F8.2))') "Total time bound. cond.:", mintimes(19),        &
    avetimes(19), maxtimes(19), percenttimes(19), avetimes(19)/nsteps*1e3
    WRITE(0, '(X, A25, X, F8.2, X, F8.2, X, F8.2)') "Total time:", mintimes(20),      &
    avetimes(20), maxtimes(20)
  ENDIF
ENDIF
END SUBROUTINE time_statistics_per_iteration

! ________________________________________________________________________________________
!> @brief
!> Subroutine that gets total memory size in Bytes occupied by tile structures grid
!> arrays on local rank. These memory sizes are stored in the mem_status module 
!> variable local_grid_mem
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________
SUBROUTINE get_local_grid_mem()
  USE mem_status, ONLY: local_grid_mem
  USE mpi
  USE mpi_type_constants, ONLY: status
  USE picsar_precision, ONLY: idp, num
  USE shared_data, ONLY: absorbing_bcs, divb, dive, divj, rho, rhoold
  IMPLICIT NONE 
  INTEGER(idp) :: field_size_f, field_size_r, cc_mat_size
  local_grid_mem = 0._num
  
  ! -- Update of memory size occupied by real-space arrays 
  local_grid_mem = local_grid_mem + SIZEOF(ex)
  local_grid_mem = local_grid_mem + SIZEOF(ey)
  local_grid_mem = local_grid_mem + SIZEOF(ez)
  local_grid_mem = local_grid_mem + SIZEOF(bx)
  local_grid_mem = local_grid_mem + SIZEOF(by)
  local_grid_mem = local_grid_mem + SIZEOF(bz)
  local_grid_mem = local_grid_mem + SIZEOF(jx)
  local_grid_mem = local_grid_mem + SIZEOF(jy)
  local_grid_mem = local_grid_mem + SIZEOF(jz)
  local_grid_mem = local_grid_mem + SIZEOF(rho)
  local_grid_mem = local_grid_mem + SIZEOF(rhoold)
  local_grid_mem = local_grid_mem + SIZEOF(dive)
  local_grid_mem = local_grid_mem + SIZEOF(divj)
  local_grid_mem = local_grid_mem + SIZEOF(divb)
  IF(absorbing_bcs) THEN
    local_grid_mem = local_grid_mem + SIZEOF(sigma_x_e)*6.0_num
    local_grid_mem = local_grid_mem + SIZEOF(exy) * 12.0_num
  ENDIF

#if defined(FFTW)
  ! -- Update of memory size occupied by Fourier arrays 
  IF(g_spectral) THEN
    IF(absorbing_bcs) THEN
      field_size_f = SIZEOF(exf)*(17.0_num+12.0_num) 
    ELSE 
      field_size_f = SIZEOF(exf)*(11.0_num+6.0_num)
    ENDIF
  ELSE 
    field_size_f = SIZEOF(exf)*11.0_num
  ENDIF
  local_grid_mem = local_grid_mem + field_size_f
  IF(absorbing_bcs) THEN
    field_size_r = SIZEOF(ex_r)*17.0_num  
    ! -- Memory occupied by blocs
    cc_mat_size = SIZEOF(exf)*(34.0_num) 
  ELSE 
    field_size_r = SIZEOF(ex_r)*11.0_num 
    ! -- Memory occupued by blocs
    cc_mat_size = SIZEOF(exf)*(33.0_num)
  ENDIF
  local_grid_mem =  local_grid_mem + field_size_r + cc_mat_size
  
#endif
END SUBROUTINE get_local_grid_mem 

! ________________________________________________________________________________________
!> @brief
!> Subroutine that computes the total memory size in Bytes occupied by grid arrays 
!> on all ranks. These memory sizes are reduced on rank 0 in the mem_status module 
!> variable: global_grid_mem
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________
SUBROUTINE get_global_grid_mem()
  USE mem_status, ONLY: global_grid_mem, local_grid_mem
  USE picsar_precision, ONLY: isp
  IMPLICIT NONE 

  ! - Estimate total grid arrays memory (reduce on proc 0)
  CALL MPI_REDUCE(local_grid_mem, global_grid_mem,                                     &
  1_isp, mpidbl, MPI_SUM, 0_isp,comm,errcode)
END SUBROUTINE get_global_grid_mem 

END MODULE mpi_routines

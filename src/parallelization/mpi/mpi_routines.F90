! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software, ! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the
! U.S.g,'cooc'! Government has been granted for itself and others acting on its behalf a
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
    INTEGER(isp) :: iret
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
    USE mpi_fftw3
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
    LOGICAL(isp) :: isinitialized
    INTEGER(isp) :: nproc_comm, rank_in_comm
    INTEGER(idp), OPTIONAL, INTENT(IN) :: comm_in

    CALL MPI_INITIALIZED(isinitialized, errcode)
    IF (.NOT. isinitialized) CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided,        &
    errcode)
    IF (present(comm_in) .AND. comm_in .GE. 0) THEN
      CALL MPI_COMM_DUP(INT(comm_in, isp), comm, errcode)
    ELSE
      CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    ENDIF
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
    nproc=INT(nproc_comm, idp)

    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
    rank=INT(rank_in_comm, idp)

  END SUBROUTINE mpi_minimal_init_python

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
#if defined(FFTW) 
    USE group_parameters , ONLY : nb_group
#endif
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
          WRITE(0, *) 'WRONG CPU SPLIT nlocal<nguards'
          WRITE(0, *) nxsplit, nysplit, nzsplit
          WRITE(0, *) nx_global_grid, ny_global_grid, nz_global_grid
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
      ENDIF
    ELSE IF(c_dim ==2) THEN 
      IF (nxsplit .LT. nxguards .OR. nzsplit .LT.  nzguards)  THEN
        IF (rank .EQ. 0) THEN
          WRITE(0, *) 'WRONG CPU SPLIT nlocal<nguards'
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
END SUBROUTINE setup_communicator

!> @brief
!> This routine creates mpi groups and communicators while using
!fftw_hybrid=.TRUE.
!or hybrid=.TRUE.
!> Author : Haithem Kallala
!> 2017
SUBROUTINE setup_groups
#if defined(FFTW)
USE  group_parameters
USE mpi_fftw3
#endif
USE picsar_precision
USE shared_data
USE mpi

INTEGER(isp), PARAMETER :: ndims = 3
LOGICAL(isp) :: periods(ndims), reorder
INTEGER(isp) :: dims(ndims), ierr
INTEGER(idp) :: group_size
INTEGER(isp), ALLOCATABLE, DIMENSION(:)    :: grp_id, grp_comm, local_roots_rank
INTEGER(isp), ALLOCATABLE, DIMENSION(:, :) :: grp_ranks
INTEGER(isp)                                :: roots_grp, roots_comm
INTEGER(idp)  :: i,j,temp
INTEGER(idp), DIMENSIOn(:), ALLOCATABLE :: all_nz_group, all_iz_max_r, all_iz_min_r,  &
all_cells, all_nz_lb, all_nzp, all_iz_min_lbg, all_iz_max_lbg
INTEGER(idp)                            :: iz_min_lbg, iz_max_lbg
#if defined(FFTW)
#if defined(DEBUG)
  WRITE(0, *) "setup_groups : start"
#endif

  CALL MPI_COMM_GROUP(comm, MPI_WORLD_GROUP, errcode)
  
  nb_group = nprocx*nprocy*nb_group
  
  group_size = nproc/nb_group
  
  ALLOCATE(grp_id(nb_group), grp_comm(nb_group), local_roots_rank(nb_group))
  ALLOCATE(grp_ranks(group_size, nb_group))

  ALLOCATE(MPI_GROUP_ID(nb_group),MPI_COMM_GROUP_ID(nb_group))

  DO i=1, nb_group
    MPI_COMM_GROUP_ID(i) = MPI_COMM_NULL
    MPI_GROUP_ID(i) = MPI_GROUP_NULL
    grp_comm(i)  = MPI_COMM_NULL
    grp_id(i) = MPI_GROUP_NULL
  ENDDO

  DO j=1, nprocx*nprocy
    local_roots_rank(j) = INT(j-1, isp)
  ENDDO
  DO j=nprocx*nprocy+1, nb_group
    i = j-nprocx*nprocy
    local_roots_rank(j) = local_roots_rank(i)+group_size*nprocx*nprocy
  ENDDO
  DO j = 1, nb_group
    grp_ranks(1, j) = local_roots_rank(j)
    DO i = 2, group_size
      grp_ranks(i, j) = grp_ranks(i-1, j) + nprocx*nprocy
    ENDDO
  ENDDO
 
  DO i= 1, nb_group
    CALL MPI_GROUP_INCL(MPI_WORLD_GROUP, INT(group_size, isp), grp_ranks(:, i),       &
    grp_id(i), errcode)
    CALL MPI_COMM_CREATE_GROUP(comm, grp_id(i), 0, grp_comm(i), errcode)
  ENDDO
  
  roots_grp = MPI_GROUP_NULL
  roots_comm = MPI_COMM_NULL
  
  CALL MPI_GROUP_INCL(MPI_WORLD_GROUP, INT(nb_group, isp), local_roots_rank,          &
  roots_grp,  errcode)
  
  CALL MPI_COMM_CREATE_GROUP(comm, roots_grp, 0, roots_comm, errcode)
  MPI_ROOT_GROUP = MPI_GROUP_NULL
  MPI_ROOT_COMM  = MPI_COMM_NULL
  root_rank = MPI_PROC_NULL
  
  IF( roots_comm .NE. MPI_COMM_NULL) THEN
    CALL MPI_COMM_RANK(roots_comm, root_rank, errcode)
    CALL MPI_COMM_SIZE(roots_comm, root_size, errcode)
    MPI_ROOT_COMM = roots_comm
    MPI_ROOT_GROUP = roots_grp
  ENDIF
  
  !  ALLOCATE(MPI_GROUP_ID(INT(nb_group, isp)), MPI_COMM_GROUP_ID(INT(nb_group, isp)))
  periods = (/.FALSE., .FALSE., .FALSE./)
  dims = (/INT(group_size, isp), 1_isp, 1_isp/)
  reorder = .TRUE.
  
  DO i = 1, nb_group
    IF (grp_comm(i) .NE. MPI_COMM_NULL) THEN
     
      CALL MPI_CART_CREATE(grp_comm(i), ndims, dims, periods, reorder,                &
      MPI_COMM_GROUP_ID(i), errcode)
      CALL MPI_COMM_GROUP(MPI_COMM_GROUP_ID(i), MPI_GROUP_ID(i), errcode)
      CALL MPI_COMM_SIZE(MPI_COMM_GROUP_ID(i), local_size, errcode)
      CALL MPI_COMM_RANK(MPI_COMM_GROUP_ID(i), local_rank, errcode)
      CALL MPI_CART_COORDS(MPI_COMM_GROUP_ID(i), local_rank, ndims, group_coordinates,&
      errcode)
      
      which_group = i-1
      
      z_group_coords = which_group/nprocx/nprocy
      y_group_coords = (which_group-z_group_coords*nprocx*nprocy)/nprocx
      x_group_coords = (which_group-z_group_coords*nprocx*nprocy) -                   &
      nprocx*y_group_coords
  
    ENDIF
  ENDDO
  
  DO i = 1, nb_group
    IF(grp_comm(i) .NE. MPI_COMM_NULL) THEN
      CALL MPI_COMM_FREE(grp_comm(i), errcode)
      CALL MPI_GROUP_FREE(grp_id(i), errcode)
    ENDIF
  ENDDO
  
  is_on_boarder = .FALSE.
  group_z_min_boundary = .FALSE.
  group_z_max_boundary = .FALSE.
  
  IF(local_rank  .EQ. 0_idp) THEN
    is_on_boarder = .TRUE.
    group_z_min_boundary = .TRUE.
  ENDIF
  
  IF(local_rank .EQ. local_size-1) THEN
    is_on_boarder = .TRUE.
    group_z_max_boundary = .TRUE.
  ENDIF
  
  nx_group_global = nx
  ny_group_global = ny
  
  nz_group_global = nz_global/(nb_group/(nprocx*nprocy))
   
  IF(nz_global .NE. nz_group_global*(nb_group/(nprocx*nprocy))) THEN
    temp = INT(nz_global - nb_group*nz_group_global/nprocx/nprocy, idp)
    IF(INT(nb_group/nprocx/nprocy-1-z_group_coords, idp) .LT. temp) THEN
      nz_group_global=nz_group_global+1
    ENDIF
  ENDIF
  ! THE NEXT ERROR COULD BE REMOVED WITH ADVANCED LOAD BALANCING 
  IF(nz_group_global .LE. 2*nzg_group) THEN
    WRITE(*,*) '**ERROR **, nz_group too small compared to nzg_group'
    CALL MPI_ABORT(comm,errcode,ierr)
  ENDIF
  
  y_min_group=ymin
  y_max_group=ymax
  x_min_group=xmin
  x_max_group=xmax
  z_min_group=zmin
  z_max_group=zmax
  
  IF(MPI_ROOT_COMM .NE. MPI_COMM_NULL) THEN
  
    ALLOCATE(all_nz_group(nb_group))
    CALL MPI_ALLGATHER(nz_group_global, 1_isp, MPI_LONG_LONG_INT, all_nz_group, 1_isp,&
    MPI_LONG_LONG_INT, MPI_ROOT_COMM, errcode)
    CALL MPI_BARRIER(MPI_ROOT_COMM, errcode)
    IF(root_rank .NE. 0_isp) THEN
      DO i=1, z_group_coords
        z_min_group = z_min_group + all_nz_group(i)*dz
      ENDDO
    ENDIF
    z_max_group = z_min_group + nz_group_global*dz
    DEALLOCATE(all_nz_group)
  
  ENDIF
  
  DO i=1, nb_group
  
    IF(MPI_COMM_GROUP_ID(i)  .NE. MPI_COMM_NULL) THEN
      CALL MPI_BCAST(z_min_group, 1_isp, MPI_DOUBLE, 0_isp, MPI_COMM_GROUP_ID(i),     &
      errcode)
      CALL MPI_BCAST(z_max_group, 1_isp, MPI_DOUBLE, 0_isp, MPI_COMM_GROUP_ID(i),     &
      errcode)
    ENDIF
  
  ENDDO

  IF(c_dim ==  2 ) THEN
    nyg_group = 0 
  ENDIF 
  nx_group_global_grid = nx_group_global+1
  ny_group_global_grid = ny_group_global+1
  nz_group_global_grid = nz_group_global+1
  
  nx_group = nx_group_global + 2*nxg_group
  ny_group = ny_group_global + 2*nyg_group
  nz_group = nz_group_global + 2*nzg_group
  
  nx_group_grid = nx_group_global_grid + 2*nxg_group
  ny_group_grid = ny_group_global_grid + 2*nyg_group
  nz_group_grid = nz_group_global_grid + 2*nzg_group
   
  DO i=1, nb_group
    IF(MPI_COMM_GROUP_ID(i)  .NE. MPI_COMM_NULL) THEN
      IF(c_dim == 3) THEN
        IF(fftw_mpi_transpose) THEN
          alloc_local = fftw_mpi_local_size_3d_transposed(nz_group, ny_group,         &
          nx_group, MPI_COMM_GROUP_ID(i), local_nz, local_z0, local_ny, local_y0)
          IF(local_nz .EQ. 0_idp .OR. local_ny .EQ. 0_idp) THEN
            WRITE(0,*) 'ERROR local_ny or local_nz = 0 in rank',rank
            CALL MPI_ABORT(comm,errcode,ierr)
          ENDIF
        ELSE
          alloc_local = FFTW_MPI_LOCAL_SIZE_3D(nz_group, ny_group, nx_group,          &
          MPI_COMM_GROUP_ID(i), local_nz, local_z0)
          IF(local_nz .EQ. 0_idp ) THEN
            WRITE(0,*) 'ERROR local_nz = 0 in rank',rank
            CALL MPI_ABORT(comm,errcode,ierr)
          ENDIF
        ENDIF
      ELSE IF(c_dim == 2) THEN
        alloc_local = FFTW_MPI_LOCAL_SIZE_2D(nz_group, nx_group,                      &
        MPI_COMM_GROUP_ID(i), local_nz, local_z0)
        IF(local_nz .EQ. 0_idp ) THEN
          WRITE(0,*) 'ERROR local_nz = 0 in rank',rank
          CALL MPI_ABORT(comm,errcode,ierr)
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  
  iz_min_r = 1
  iz_max_r = local_nz
        
  IF(group_z_min_boundary) iz_min_r = nzg_group+1
  IF(group_z_max_boundary) iz_max_r = local_nz-nzg_group

  ! case where some procs solve maxwell in ghost region exclusively is only
  ! supported by load balancing
  IF(.NOT. is_lb_grp) THEN 
    IF(iz_min_r .GT. iz_max_r) THEN
      WRITE(0, *) '*** ERROR ***'
      WRITE(0, *) '*** iz_max_r < iz_min_r ***'
      CALL MPI_ABORT(comm,errcode,ierr)
    ENDIF
  ENDIF
  IF(is_lb_grp) THEN
    IF(0 .NE. MODULO(nz_global,nb_group/(nprocx*nprocy))) THEN
     WRITE(0,*) 'PROBLEM HERE, PLEASE SET nz as a multiple of nb_group'
     CALL MPI_ABORT(comm,errcode,ierr)
    ENDIF
    iz_min_lbg = local_z0-nzg_group + nz_group_global*z_group_coords
    iz_max_lbg = iz_min_lbg + local_nz - 1 
    ALLOCATE(cell_z_min_lbg(nprocz),cell_z_max_lbg(nprocz))
    ALLOCATE(all_iz_min_lbg(nproc),all_iz_max_lbg(nproc))

    CALL MPI_ALLGATHER(iz_min_lbg,1_isp, MPI_LONG_LONG_INT, all_iz_min_lbg,           &
    INT(1,isp), MPI_LONG_LONG_INT, comm, errcode)

    CALL MPI_ALLGATHER(iz_max_lbg,1_isp, MPI_LONG_LONG_INT, all_iz_max_lbg,           &
    INT(1,isp), MPI_LONG_LONG_INT, comm, errcode)

    DO i=1, nprocz
      cell_z_min_lbg(i) = all_iz_min_lbg((i-1)*nprocx*nprocy+1)
      cell_z_max_lbg(i) = all_iz_max_lbg((i-1)*nprocx*nprocy+1)
    ENDDO
    DEALLOCATE(all_iz_max_lbg,all_iz_min_lbg)
    DO i = 1, nprocz
      IF(i-1 == z_coords) THEN
        IF(cell_z_max_lbg(i) - cell_z_min_lbg(i) + 1  .NE. local_nz) THEN
          WRITE(*,*) 'ERROR in load balancing indexes'
          CALL MPI_ABORT(comm, errcode,ierr)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  !to take into account odd cases
  nz_lb = MAX(iz_max_r - iz_min_r +1,0)
  nz_grid_lb = nz_lb+1

  ALLOCATE(all_nz_lb(nproc), all_nzp(nprocz))
  
  CALL MPI_ALLGATHER(nz_lb, 1_isp, MPI_LONG_LONG_INT, all_nz_lb, 1_isp      ,         &
  MPI_LONG_LONG_INT, comm, errcode)
  
  DO i=1, nprocz
    all_nzp(i) = all_nz_lb((i-1)*nprocx*nprocy+1)
  ENDDO

  z_min_local_lb = zmin
  z_max_local_lb = zmax
  z_min_local_lb = z_min_local_lb + dz*sum(all_nzp(1:z_coords))
  z_max_local_lb = z_min_local_lb +nz_lb*dz
  cell_z_min_f(1) = 0
  cell_z_max_f(1) = all_nzp(1) - 1
  
  DO i =2, nprocz
    cell_z_min_f(i) = cell_z_max_f(i-1) + 1
    cell_z_max_f(i) = cell_z_min_f(i)-1 + all_nzp(i)
  ENDDO
  
  ix_min_r = 1
  ix_max_r = nx + 2*nxguards
  
  iy_min_r = 1
  iy_max_r = ny + 2*nyguards
  
  nz_global_grid_min_lb = cell_z_min_f(z_coords+1)
  nz_global_grid_max_lb = cell_z_max_f(z_coords+1)+1
  
  DEALLOCATE(all_nz_lb, all_nzp)
  DEALLOCATE(grp_id, grp_comm, local_roots_rank, grp_ranks)

#if defined(DEBUG)
  WRITE(0, *) "setup_groups : end"
#endif

#endif

END SUBROUTINE setup_groups

!> Computes r_local and g_local arrays (corresponding indexes between field_r(at
!mpi task rank) and field(at mpi task rank)
!> Haithem kallala 2017
SUBROUTINE compute_load_balancing_local()
#if defined(FFTW)
  USE mpi_fftw3
  USE shared_data
  USE group_parameters
  USE picsar_precision
  USE params
  USE mpi 
  INTEGER(idp)  :: index_r_first, index_r_last, index_g_first, index_g_last
  INTEGER(idp)  :: g_f,g_l,r_f,r_l,i
  INTEGER(isp)  :: ierr

  index_r_first = cell_z_min_r(z_coords+1)
  index_g_first = cell_z_min_f(z_coords+1)
  index_r_last = cell_z_max_r(z_coords+1) 
  index_g_last = cell_z_max_f(z_coords+1)
  g_f = MAX(index_g_first,index_r_first) - index_g_first
  r_f = MAX(index_g_first,index_r_first) - index_r_first
  g_l = MIN(index_g_last,index_r_last) - index_g_first
  r_l = MIN(index_g_last,index_r_last) - index_r_first
!CHECK NOT WRONG 
  IF((r_l - r_f) .NE. (g_l - g_f)) THEN
    WRITE(0,*) 'error in load balancing 1 '
    CALL MPI_ABORT(comm,errcode,ierr)
  ENDIF
!  IF((r_l - r_f .LT. 0_idp) .OR. (g_l - g_f .LT. 0_idp)) THEN
!    WRITE(0,*) 'error in load balancing 2' 
!    CALL MPI_ABORT(comm,errcode,ierr)
!  ENDIF
  size_local = MIN(r_l - r_f +1,0_idp)
  ALLOCATE(g_local(size_local),r_local(size_local))
  IF(size_local .GE. 1_idp) THEN
    g_local(1) = g_f + 1
    IF(group_z_min_boundary) g_local(1) = g_local(1) + nzg_group
    r_local(1) = r_f
    DO i = 2,size_local
      g_local(i) = g_local(i-1) +1
      r_local(i) = r_local(i-1) +1
    ENDDO
  ENDIF
#endif
END SUBROUTINE compute_load_balancing_local

!>Computes r_left and g_left (corresponding indexes between field_r(at mpi task
!>rank) and field(at mpi task rank-1)
!>Haithem Kallala 2017
SUBROUTINE compute_load_balancing_from_left()
#if defined(FFTW)
  USE mpi_fftw3
  USE shared_data
  USE group_parameters
  USE picsar_precision
  USE params
  USE mpi

  INTEGER(idp)  :: index_r_first, index_r_last, index_g_first, index_g_last
  INTEGER(idp)  :: g_f,g_l,r_f,r_l,i
  INTEGER(isp)  :: ierr
  IF(z_coords == 0) THEN
   size_left = 0_idp
    ELSE
    index_r_first = cell_z_min_r(z_coords)
    index_g_first = cell_z_min_f(z_coords+1)
    index_r_last = cell_z_max_r(z_coords)
    index_g_last = cell_z_max_f(z_coords+1)

    IF(index_r_last .LT. index_g_first)  THEN
     size_left = 0_idp
      ELSE
      g_f = MAX(index_g_first,index_r_first) - index_g_first
      r_f = MAX(index_g_first,index_r_first) - index_r_first
      g_l = MIN(index_g_last,index_r_last) - index_g_first
      r_l = MIN(index_g_last,index_r_last) - index_r_first
    !CHECK NOT WRONG 
      IF((r_l - r_f) .NE. (g_l - g_f)) THEN
        WRITE(0,*) 'error in load balancing 11 '
        CALL MPI_ABORT(comm,errcode,ierr)
      ENDIF
 !     IF((r_l - r_f .LT. 0_idp) .OR. (g_l - g_f .LT. 0_idp)) THEN
 !       WRITE(0,*) 'error in load balancing 21'
 !       CALL MPI_ABORT(comm,errcode,ierr)
 !     ENDIF
      size_left = MIN(r_l - r_f +1,0_idp)
      ALLOCATE(g_left(size_left),r_left(size_left))
      IF(size_left .GE. 1_idp) THEN
        g_left(1) = g_f + 1
        r_left(1) = r_f
        DO i = 2,size_left
          g_left(i) = g_left(i-1) +1
          r_left(i) = r_left(i-1) +1
        ENDDO
      ENDIF
    ENDIF
  ENDIF
#endif
END SUBROUTINE compute_load_balancing_from_left

!>Computes corresponding r_right and g_right indexes between field_r(at mpi task
!rank) and field(at mpi task rank+1)
SUBROUTINE compute_load_balancing_from_right()
#if defined(FFTW)
  USE mpi_fftw3
  USE shared_data
  USE group_parameters
  USE picsar_precision
  USE params
  USE mpi

  INTEGER(idp)  :: index_r_first, index_r_last, index_g_first, index_g_last
  INTEGER(idp)  :: g_f,g_l,r_f,r_l,i
  INTEGER(isp)  :: ierr

  IF(z_coords == nprocz-1) THEN
   size_right = 0_idp
  ELSE
    index_r_first = cell_z_min_r(z_coords+2)
    index_g_first = cell_z_min_f(z_coords+1)
    index_r_last = cell_z_max_r(z_coords+2)
    index_g_last = cell_z_max_f(z_coords+1)
    IF(index_g_last .LT. index_r_first) THEN
      size_right = 0_idp
    ELSE
      g_f = MAX(index_g_first,index_r_first) - index_g_first
      r_f = MAX(index_g_first,index_r_first) - index_r_first
      g_l = MIN(index_g_last,index_r_last) - index_g_first
      r_l = MIN(index_g_last,index_r_last) - index_r_first
    !CHECK NOT WRONG 
      IF((r_l - r_f) .NE. (g_l - g_f)) THEN
        WRITE(0,*) 'error in load balancing 12 '
        CALL MPI_ABORT(comm,errcode,ierr)
      ENDIF
   !   IF((r_l - r_f .LT. 0_idp) .OR. (g_l - g_f .LT. 0_idp)) THEN
   !     WRITE(0,*) 'error in load balancing 22',rank,r_l,r_f,index_r_first,index_r_last,index_g_last
   !     CALL MPI_ABORT(comm,errcode,ierr)
   !   ENDIF
      size_right = MIN(r_l - r_f +1,0_idp)
      ALLOCATE(g_right(size_right),r_right(size_right))
      IF(size_right .GE. 1_idp) THEN
        g_right(1) = g_f + 1
        IF(group_z_min_boundary) g_local(1) = g_local(1) + nzg_group
        r_right(1) = r_f
        DO i = 2,size_right
          g_right(i) = g_right(i-1) +1
          r_right(i) = r_right(i-1) +1
        ENDDO
      ENDIF
    ENDIF
  ENDIF
#endif
END SUBROUTINE compute_load_balancing_from_right
!>Send r_right and g_right to rr_left and rg_left of rank+1  respectively
SUBROUTINE Sync_exchange_load_balancing_arrays_1
#if defined(FFTW)
  USE group_parameters
  USE shared_data
  USE mpi
  INTEGER(isp)  :: requests(1)
 IF(z_coords .NE. nprocz-1) THEN
   CALL MPI_ISEND(size_right,1_isp,MPI_LONG_LONG_INT,INT(proc_z_max,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. 0) THEN
   CALL MPI_IRECV(rsize_left,1_isp,MPI_LONG_LONG_INT,INT(proc_z_min,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF

 ALLOCATE(rr_left(rsize_left))
 ALLOCATE(rg_left(rsize_left))

 IF(z_coords .NE. nprocz-1) THEN
   CALL MPI_ISEND(r_right,size_right,MPI_LONG_LONG_INT,INT(proc_z_max,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. 0) THEN
   CALL MPI_IRECV(rr_left,rsize_left,MPI_LONG_LONG_INT,INT(proc_z_min,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF  
 IF(z_coords .NE. nprocz-1) THEN
   CALL MPI_ISEND(g_right,size_right,MPI_LONG_LONG_INT,INT(proc_z_max,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. 0) THEN
   CALL MPI_IRECV(rg_left,rsize_left,MPI_LONG_LONG_INT,INT(proc_z_min,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
#endif
END SUBROUTINE Sync_exchange_load_balancing_arrays_1

!>Sends r_left and g_left to rr_right and rg_right of rank-1 respectively
SUBROUTINE Sync_exchange_load_balancing_arrays_2
#if defined(FFTW)
  USE group_parameters
  USE shared_data
  USE mpi
  INTEGER(isp)  :: requests(1)
 IF(z_coords .NE. 0) THEN
   CALL MPI_ISEND(size_left,1_isp,MPI_LONG_LONG_INT,INT(proc_z_min,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. nprocz-1) THEN
   CALL MPI_IRECV(rsize_right,1_isp,MPI_LONG_LONG_INT,INT(proc_z_max,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF

 ALLOCATE(rr_right(rsize_right))
 ALLOCATE(rg_right(rsize_right))

 IF(z_coords .NE. 0) THEN
   CALL MPI_ISEND(r_left,size_left,MPI_LONG_LONG_INT,INT(proc_z_min,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. nprocz-1) THEN
   CALL MPI_IRECV(rr_right,rsize_right,MPI_LONG_LONG_INT,INT(proc_z_max,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. 0) THEN
   CALL MPI_ISEND(g_left,size_left,MPI_LONG_LONG_INT,INT(proc_z_min,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
 IF(z_coords .NE. nprocz-1) THEN
   CALL MPI_IRECV(rg_right,rsize_right,MPI_LONG_LONG_INT,INT(proc_z_max,isp),tag,comm,requests(1),errcode)
   CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
 ENDIF
#endif
END SUBROUTINE Sync_exchange_load_balancing_arrays_2


SUBROUTINE adjust_grid_mpi_global

#if defined(FFTW)
  USE mpi_fftw3
  USE shared_data
  USE mpi
  USE picsar_precision
  INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: all_nz
  INTEGER(idp)  :: idim
  
  nz = local_nz
  nz_grid = nz + 1
  
  ALLOCATE(all_nz(1:nprocz))
  
  CALL MPI_ALLGATHER(nz, 1_isp, MPI_LONG_LONG_INT, all_nz, 1_isp, MPI_LONG_LONG_INT,    &
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
!
!> @date
!> Creation 2015
! ______________________________________________________________________________________
SUBROUTINE mpi_initialise
#if defined(FFTW)
USE mpi_fftw3
USE group_parameters
USE load_balance
#endif
INTEGER(isp) :: idim,ierr
INTEGER(isp) :: nx0, nxp
INTEGER(isp) :: ny0, nyp
INTEGER(isp) :: nz0, nzp
INTEGER(idp) :: iz
#if defined(FFTW)
INTEGER(C_INTPTR_T) :: kx, ly, mz
INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: nz_procs, all_nz
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

CALL setup_communicator

ALLOCATE(x_grid_mins(1:nprocx), x_grid_maxs(1:nprocx))
ALLOCATE(y_grid_mins(1:nprocy), y_grid_maxs(1:nprocy))
ALLOCATE(z_grid_mins(1:nprocz), z_grid_maxs(1:nprocz))
ALLOCATE(cell_x_min(1:nprocx), cell_x_max(1:nprocx))
ALLOCATE(cell_y_min(1:nprocy), cell_y_max(1:nprocy))
ALLOCATE(cell_z_min(1:nprocz), cell_z_max(1:nprocz))

#if defined(FFTW)
! With fftw_with_mpi CPU split is performed only along z
IF (fftw_with_mpi) THEN
  mz=nz_global; ly=ny_global; kx=nx_global
  !   get local data size and allocate (note dimension reversal)
  IF(.NOT. fftw_mpi_transpose) THEN
    alloc_local = fftw_mpi_local_size_3d(mz, ly, kx/2+1, comm, local_nz, local_z0)
  ELSE
    alloc_local = fftw_mpi_local_size_3d_transposed(mz, ly, kx/2+1, comm, local_nz,   &
    local_z0, local_ny, local_y0)
  ENDIF
  ALLOCATE(nz_procs(nproc))
  CALL MPI_ALLGATHER(INT(local_nz, idp), 1_isp, MPI_INTEGER8, nz_procs, INT(1, isp),  &
  MPI_INTEGER8, comm, errcode)
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

!!! --- if fftw_with_mpi = true then adjust nz to b equal to local_nz (since the
!two are computed differently

#if defined(FFTW)

IF(fftw_with_mpi ) THEN
  IF(.NOT. fftw_hybrid) THEN
    CALL adjust_grid_mpi_global
    IF(nz .NE. cell_z_max(z_coords+1) - cell_z_min(z_coords+1)+1) THEN
      WRITE(*, *), 'ERROR IN AJUSTING THE GRID 1'
      STOP
    ENDIF
  ELSE
    ALLOCATE(cell_z_min_r(1:nprocz),cell_z_max_r(1:nprocz))
    cell_z_min_r = cell_z_min
    cell_z_max_r = cell_z_max
    ALLOCATE(cell_z_min_f(1:nprocz),cell_z_max_f(1:nprocz))
    CALL setup_groups
    IF(nz_lb .NE. cell_z_max_f(z_coords+1) - cell_z_min_f(z_coords+1)+1) THEN
      WRITE(*, *), 'ERROR IN AJUSTING THE GRID 2'
      STOP
    ENDIF
  ENDIF
ENDIF
IF(fftw_hybrid) THEN 
  IF(is_lb_grp) THEN
    call get1D_intersection_group_mpi()
call mpi_barrier(comm,errcode)
    !> Computes field array indexes
   ! CALL compute_load_balancing_local
   ! CALL compute_load_balancing_from_left
   ! CALL compute_load_balancing_from_right
   ! !> Synchronizes array indexes between adjacent procs
   ! CALL Sync_exchange_load_balancing_arrays_2
   ! CALL Sync_exchange_load_balancing_arrays_1
!    IF(size_local + size_left + size_right .NE. nz_lb) THEN
!      WRITE(*,*), 'ERROR IN LOAD BALANCING'
!      CALL MPI_ABORT(comm,errcode,ierr) 
!    ENDIF
!    IF(size_local + rsize_left+ rsize_right .NE. nz) THEN
!      WRITE(*,*)  'EROOR IN LOAD BALANCING(REAL Indice)'
!      CALL MPI_ABORT(comm,errcode,ierr)
!    ENDIF
  ! IF is_lb_grp == .FALSE. THEN set grid DD params to that of
  ! fftw_mpi_local_size (Unbalanced grid for particles)
  ELSE IF(.NOT. is_lb_grp) THEN
    cell_z_min = cell_z_min_f
    cell_z_max = cell_z_max_f
print*,cell_z_min_f,cell_z_max_f
    nz = nz_lb
    nz_grid = nz_grid_lb 
    nz_global_grid_min = nz_global_grid_min_lb
    nz_global_grid_max = nz_global_grid_max_lb
    z_min_local = z_min_local_lb
    z_max_local = z_max_local_lb
  ENDIF
ENDIF
#endif
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
print*,"cocococo",size(y_global),cell_y_max
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
DO iproc = 1, nprocy
  y_grid_mins(iproc) = y_global(cell_y_min(iproc))
  y_grid_maxs(iproc) = y_global(cell_y_max(iproc)+1)
ENDDO
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
#if defined(FFTW)
USE fourier
USE mpi_fftw3
USE group_parameters
#endif
IMPLICIT NONE
#if defined(FFTW)
TYPE(C_PTR) :: cdata, cin
INTEGER(idp) :: imn, imx, jmn, jmx, kmn, kmx
#endif
! --- Allocate regular grid quantities (in real space)
ALLOCATE(ex(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(ey(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(ez(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(bx(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(by(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(bz(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
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
  IF (fftw_with_mpi) THEN
    nkx=(nx_global)/2+1! Real To Complex Transform
    nky=ny_global
    nkz=local_nz
    IF(fftw_mpi_transpose) THEN
      nkx=(nx_global)/2+1! Real To Complex Transform
      nky=nz_global
      nkz=local_ny
    ENDIF
    IF(fftw_hybrid)  THEN
      nkx = nx_group/2+1
      nky = ny_group
      nkz = local_nz
      IF(fftw_mpi_transpose) THEN
        nkx = nx_group/2+1
        nky = nz_group
        nkz = local_ny
      ENDIF
    ENDIF
    ! - Allocate complex arrays
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
    cdata = fftw_alloc_complex(alloc_local)
    ! - Allocate real arrays
    IF(fftw_mpi_transpose .AND. fftw_hybrid) THEN
      nkx = nx_group/2+1
      nky = ny_group
      nkz = local_nz
    ENDIF
    IF(fftw_mpi_transpose .AND. .NOT. fftw_hybrid) THEN
      nkx = nx_global/2+1
      nky = ny_global
      nkz = local_nz
    ENDIF
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, ex_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, ey_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, ez_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, bx_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, by_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, bz_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, jx_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, jy_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, jz_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, rho_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
    CALL c_f_pointer(cin, rhoold_r, [2*nkx, nky, nkz])
    cin = fftw_alloc_real(2 * alloc_local);
  ELSE
    nkx=(2*nxguards+nx)/2+1! Real To Complex Transform
    nky=(2*nyguards+ny)
    nkz=(2*nzguards+nz)
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
    imn=-nxguards; imx=nx+nxguards-1
    jmn=-nyguards;jmx=ny+nyguards-1
    kmn=-nzguards;kmx=nz+nzguards-1
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
  ENDIF
ENDIF
#endif
! --- Quantities used by the dynamic load balancer
ALLOCATE(new_cell_x_min(1:nprocx), new_cell_x_max(1:nprocx))
ALLOCATE(new_cell_y_min(1:nprocy), new_cell_y_max(1:nprocy))
ALLOCATE(new_cell_z_min(1:nprocz), new_cell_z_max(1:nprocz))

ALLOCATE(dive(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(divj(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
ALLOCATE(divb(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))

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
USE time_stat
USE params
#ifdef _OPENMP
USE omp_lib
#endif
IMPLICIT NONE
REAL(num), DIMENSION(25) :: percenttimes
INTEGER(idp)             :: nthreads_tot

! Total times
localtimes(20) = sum(localtimes(1:14)) + localtimes(24)
localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + localtimes(8) +      &
localtimes(11) + localtimes(13)
localtimes(18) = localtimes(5) + localtimes(6) + localtimes(7) + localtimes(8) +      &
localtimes(24)

init_localtimes(5) = sum(init_localtimes(1:4))

! Reductions
! Maximun times
CALL MPI_REDUCE(localtimes, maxtimes, 25_isp, mpidbl, MPI_MAX, 0_isp, comm, errcode)
CALL MPI_REDUCE(init_localtimes, init_maxtimes, 5_isp, mpidbl, MPI_MAX, 0_isp, comm,  &
errcode)
! Minimum times
CALL MPI_REDUCE(localtimes, mintimes, 25_isp, mpidbl, MPI_MIN, 0_isp, comm, errcode)
CALL MPI_REDUCE(init_localtimes, init_mintimes, 5_isp, mpidbl, MPI_MIN, 0_isp, comm,  &
errcode)
! Average
CALL MPI_REDUCE(localtimes, avetimes, 25_isp, mpidbl, MPI_SUM, 0_isp, comm, errcode)
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
USE time_stat
USE params
#ifdef _OPENMP
USE omp_lib
#endif
IMPLICIT NONE

REAL(num), DIMENSION(25) :: percenttimes

! Time stats per iteration activated
IF (timestat_perit.gt.0) THEN

  ! Total times
  localtimes(20) = sum(localtimes(1:14)) + localtimes(24)
  localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + localtimes(8) +    &
  localtimes(11) + localtimes(13)
  localtimes(18) = localtimes(5) + localtimes(6) + localtimes(7) + localtimes(8) +    &
  localtimes(24)

  init_localtimes(5) = sum(init_localtimes(1:4))

  ! Reductions
  ! Maximun times
  CALL MPI_REDUCE(localtimes, maxtimes, 25_isp, mpidbl, MPI_MAX, 0_isp, comm,         &
  errcode)
  CALL MPI_REDUCE(init_localtimes, init_maxtimes, 5_isp, mpidbl, MPI_MAX, 0_isp,      &
  comm, errcode)
  ! Minimum times
  CALL MPI_REDUCE(localtimes, mintimes, 25_isp, mpidbl, MPI_MIN, 0_isp, comm,         &
  errcode)
  CALL MPI_REDUCE(init_localtimes, init_mintimes, 5_isp, mpidbl, MPI_MAX, 0_isp,      &
  comm, errcode)
  ! Average
  CALL MPI_REDUCE(localtimes, avetimes, 25_isp, mpidbl, MPI_SUM, 0_isp, comm,         &
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
END MODULE mpi_routines

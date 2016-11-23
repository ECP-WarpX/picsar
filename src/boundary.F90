! ________________________________________________________________________________________
!
! BOUNDARY.F90
!
!>@author
!>Henri Vincenti,
!>Mathieu Lobet
!>
! Brief description:
!> Module containing routines for boundary conditions on fields, currents and particles
!> @brief
!>
!> For the moment this module handles:
!> periodic external boundary conditions for particles and fields
!>
! List of suboutines:
! - field_bc
! - exchange_mpi_3d_grid_array_with_guards
! - exchange_mpi_3d_grid_array_with_guards_nonblocking
! - summation_bcs
! - summation_bcs_nonblocking
! - particle_bcs_tiles_and_mpi_3d
!
! ________________________________________________________________________________________

MODULE boundary

  USE shared_data
  USE fields
  USE particles
  USE tiling
  USE mpi_derived_types
  USE constants
  USE time_stat
  USE params

  IMPLICIT NONE

  CONTAINS

  !!! --- Exchange field values at processor boundaries and apply field
  !!! --- boundary conditions
  SUBROUTINE field_bc(field, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: field
    IF (mpicom_curr.EQ.1) THEN
    	CALL exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, nx, ny, nz)
    ELSE
		CALL exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, nx, ny, nz)
    ENDIF

  END SUBROUTINE field_bc

!!! --- ROUTINE EXCHANGING GUARD REGIONS BETWEEN SUBDOMAINS (BLOCKING VERSION+DIAGONAL TRICK)
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)

    USE constants
    USE mpi
    IMPLICIT NONE

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpidbl

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    szmax = sizes(1) * sizes(2) * nzg
    sz = sizes(1) * sizes(3) * nyg
    IF (sz .GT. szmax) szmax = sz
    sz = sizes(2) * sizes(3) * nxg
    IF (sz .GT. szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

	IF (is_dtype_init(1)) THEN
    	mpi_dtypes(1) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
    	is_dtype_init(1) = .FALSE.
    ENDIF

    ! MOVE EDGES ALONG X
    CALL MPI_SENDRECV(field(1,-nyg,-nzg), 1_isp, mpi_dtypes(1), INT(proc_x_min,isp), &
        tag, temp, sz, basetype, INT(proc_x_max,isp), tag, comm, status, errcode)

    IF (proc_x_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(nx_local-nxg,-nyg,-nzg), 1_isp, mpi_dtypes(1), INT(proc_x_max,isp), &
        tag, temp, sz, basetype, INT(proc_x_min,isp), tag, comm, status, errcode)

    IF (proc_x_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

	IF (is_dtype_init(2)) THEN
    	mpi_dtypes(2) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
    	is_dtype_init(2) = .FALSE.
    ENDIF

    ! MOVE EDGES ALONG Y
    CALL MPI_SENDRECV(field(-nxg,1,-nzg), 1_isp, mpi_dtypes(2), INT(proc_y_min,isp), &
        tag, temp, sz, basetype, INT(proc_y_max,isp), tag, comm, status, errcode)

    IF (proc_y_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-nxg,ny_local-nyg,-nzg), 1_isp, mpi_dtypes(2), INT(proc_y_max,isp), &
        tag, temp, sz, basetype, INT(proc_y_min,isp), tag, comm, status, errcode)

    IF (proc_y_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF


    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg

    sz = subsizes(1) * subsizes(2) * subsizes(3)

	IF (is_dtype_init(3)) THEN
    	mpi_dtypes(3) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
    	is_dtype_init(3) = .FALSE.
    ENDIF


    ! MOVE EDGES ALONG Z
    CALL MPI_SENDRECV(field(-nxg,-nyg,1), 1_isp, mpi_dtypes(3), INT(proc_z_min,isp), &
        tag, temp, sz, basetype, INT(proc_z_max,isp), tag, comm, status, errcode)

    IF (proc_z_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-nxg,-nyg,nz_local-nzg), 1_isp, mpi_dtypes(3), INT(proc_z_max,isp), &
        tag, temp, sz, basetype, INT(proc_z_min,isp), tag, comm, status, errcode)

    IF (proc_z_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    DEALLOCATE(temp)

  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards

! ________________________________________________________________________________________
!> @brief
!> ROUTINE EXCHANGING GUARD REGIONS BETWEEN SUBDOMAINS (NON-BLOCKING VERSION+ DIAGONAL TRICK)
!>
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)
! ________________________________________________________________________________________

    USE mpi

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: basetype, sz
    INTEGER(isp):: requests(4)

    basetype = mpidbl

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    ! MOVE EDGES ALONG X
    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

	IF (is_dtype_init(4)) THEN
    	mpi_dtypes(4) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
    	is_dtype_init(4) = .FALSE.
    ENDIF

    CALL MPI_ISEND(field(1,-nyg,-nzg), 1_isp, mpi_dtypes(4), INT(proc_x_min,isp), tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(nx_local+1,-nyg,-nzg), 1_isp, mpi_dtypes(4), INT(proc_x_max,isp), tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(nx_local-nxg,-nyg,-nzg), 1_isp, mpi_dtypes(4), INT(proc_x_max,isp), tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, mpi_dtypes(4), INT(proc_x_min,isp), tag, &
        comm, requests(4), errcode)


    ! NEED TO WAIT BEFORE EXCHANGING ALONG Y (DIAGONAL TERMS)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! MOVE EDGES ALONG Y
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)

	IF (is_dtype_init(5)) THEN
    	mpi_dtypes(5) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
    	is_dtype_init(5) = .FALSE.
    ENDIF

    CALL MPI_ISEND(field(-nxg,1,-nzg), 1_isp, mpi_dtypes(5), INT(proc_y_min,isp), tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg,ny_local+1,-nzg), 1_isp, mpi_dtypes(5), INT(proc_y_max,isp), tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(-nxg,ny_local-nyg,-nzg), 1_isp, mpi_dtypes(5), INT(proc_y_max,isp), tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, mpi_dtypes(5), INT(proc_y_min,isp), tag, &
        comm, requests(4), errcode)

    ! NEED TO WAIT BEFORE EXCHANGING ALONG Z (DIAGONAL TERMS)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

   ! MOVE EDGES ALONG Z
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg

    sz = subsizes(1) * subsizes(2) * subsizes(3)

	IF (is_dtype_init(6)) THEN
    	mpi_dtypes(6) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
    	is_dtype_init(6) = .FALSE.
    ENDIF

    CALL MPI_ISEND(field(-nxg,-nyg,1), 1_isp, mpi_dtypes(6), INT(proc_z_min,isp), tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,nz_local+1), 1_isp, mpi_dtypes(6), INT(proc_z_max,isp), tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(-nxg,-nyg,nz_local-nzg), 1_isp, mpi_dtypes(6), INT(proc_z_max,isp), tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, mpi_dtypes(6), INT(proc_z_min,isp), tag, &
        comm, requests(4), errcode)

    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)


  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking


! ________________________________________________________________________________________
!> @brief
!> Routine for adding current contributions fron adjacent subdomains
  SUBROUTINE summation_bcs(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
! ________________________________________________________________________________________

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: nn, sz

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

	IF (is_dtype_init(7)) THEN
    	mpi_dtypes(7) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
    	is_dtype_init(7) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn,-nyg,-nzg), 1_isp, mpi_dtypes(7), &
        INT(neighbour( 1,0,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(-1,0,0),isp), tag, comm, status, errcode)
    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg+1,-nyg,-nzg), 1_isp, mpi_dtypes(7), &
        INT(neighbour(-1,0,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour( 1,0,0),isp), tag, comm, status, errcode)
    array(nn-nxg+1:nn,:,:) = array(nn-nxg+1:nn,:,:) + temp

    DEALLOCATE(temp)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

	IF (is_dtype_init(8)) THEN
    	mpi_dtypes(8) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
    	is_dtype_init(8) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,nn,-nzg), 1_isp, mpi_dtypes(8), &
        INT(neighbour(0, 1,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0,-1,0),isp), tag, comm, status, errcode)
    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg+1,-nzg), 1_isp, mpi_dtypes(8), &
        INT(neighbour(0,-1,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0, 1,0),isp), tag, comm, status, errcode)
    array(:,nn-nyg+1:nn,:) = array(:,nn-nyg+1:nn,:) + temp

    DEALLOCATE(temp)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

	IF (is_dtype_init(9)) THEN
    	mpi_dtypes(9) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
    	is_dtype_init(9) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,nn), 1_isp, mpi_dtypes(9) , &
        INT(neighbour(0,0, 1),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0,0,-1),isp), tag, comm, status, errcode)
    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg+1), 1_isp, mpi_dtypes(9) , &
        INT(neighbour(0,0,-1),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0,0, 1),isp), tag, comm, status, errcode)
    array(:,:,nn-nzg+1:nn) = array(:,:,nn-nzg+1:nn) + temp

    DEALLOCATE(temp)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

  END SUBROUTINE summation_bcs

!!! --- Routine for adding current contributions fron adjacent subdomains nonblocking version
  SUBROUTINE summation_bcs_nonblocking(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    USE mpi

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: nn, sz
    INTEGER(isp) :: requests(4)

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

	IF (is_dtype_init(10)) THEN
    	mpi_dtypes(10) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
    	is_dtype_init(10) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(nn,-nyg,-nzg), 1_isp, mpi_dtypes(10), INT(proc_x_max,isp), tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_x_min,isp), tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg+1,-nyg,-nzg), 1_isp, mpi_dtypes(10), INT(proc_x_min,isp), tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_x_max,isp), tag, &
    comm, requests(4), errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp1
    array(nn-nxg+1:nn,:,:) = array(nn-nxg+1:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

	IF (is_dtype_init(11)) THEN
    	mpi_dtypes(11) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
    	is_dtype_init(11) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(-nxg,nn,-nzg), 1_isp, mpi_dtypes(11), INT(proc_y_max,isp), tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_y_min,isp), tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg+1,-nzg), 1_isp, mpi_dtypes(11), INT(proc_y_min,isp), tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_y_max,isp), tag, &
    comm, requests(4), errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp1
    array(:,nn-nyg+1:nn,:) = array(:,nn-nyg+1:nn,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

	IF (is_dtype_init(12)) THEN
    	mpi_dtypes(12) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
    	is_dtype_init(12) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)),temp2(subsizes(1), subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(-nxg,-nyg,nn), 1_isp, mpi_dtypes(12), INT(proc_z_max,isp), tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_z_min,isp), tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg+1), 1_isp, mpi_dtypes(12), INT(proc_z_min,isp), tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_z_max,isp), tag, &
    comm, requests(4), errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp1
    array(:,:,nn-nzg+1:nn) = array(:,:,nn-nzg+1:nn) + temp2

    DEALLOCATE(temp1,temp2)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

  END SUBROUTINE summation_bcs_nonblocking

!!! --- Routine for adding current contributions fron adjacent subdomains persistent version
  SUBROUTINE summation_bcs_persistent_jx(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
    USE communications
    USE mpi

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarrayx, subarrayy, subarrayz, nn, sz
    INTEGER(isp) :: proc_x_min_mpisp, proc_x_max_mpisp
    INTEGER(isp) :: proc_y_min_mpisp, proc_y_max_mpisp
    INTEGER(isp) :: proc_z_min_mpisp, proc_z_max_mpisp

    proc_x_min_mpisp = proc_x_min
    proc_x_max_mpisp = proc_x_max
    proc_y_min_mpisp = proc_y_min
    proc_y_max_mpisp = proc_y_max
    proc_z_min_mpisp = proc_z_min
    proc_z_max_mpisp = proc_z_max

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    ! Initialization of the persistent communication
    IF (it .EQ. 0) THEN

      ! Init X
      subsizes(1) = nxg
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max_mpisp, tag, comm, reqperjxx(1), errcode)
      call MPI_SEND_INIT(array(-nxg+1,-nyg,-nzg), 1_isp, subarrayx, proc_x_min_mpisp, tag, comm, reqperjxx(3), errcode)

      call MPI_RECV_INIT(temp1, sz, mpidbl, proc_x_min_mpisp, tag, comm, reqperjxx(2), errcode)
      call MPI_RECV_INIT(temp2, sz, mpidbl, proc_x_max_mpisp, tag, comm, reqperjxx(4), errcode)

      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max_mpisp, tag, comm, reqperjxy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg+1,-nzg), 1_isp, subarrayy, proc_y_min_mpisp, tag, comm, reqperjxy(3), errcode)

      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max_mpisp, tag, comm, reqperjxz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg+1), 1_isp, subarrayz, proc_z_min_mpisp, tag, comm, reqperjxz(3), errcode)

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)
      CALL MPI_TYPE_FREE(subarrayz, errcode)

    ENDIF

    !! -- Summation along X- direction

    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'it',it
    !if (rank.eq.0) print*, 'summation along X'

    temp1  = 0.0_num
    temp2 = 0.0_num
    !call MPI_STARTALL(4_isp,reqperjxx,errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min_mpisp, tag, &
    comm, reqperjxx(2), errcode)
    !call MPI_START(reqperjxx(2),errcode)
    call MPI_START(reqperjxx(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max_mpisp, tag, &
    comm, reqperjxx(4), errcode)
    !call MPI_START(reqperjxx(4),errcode)
    call MPI_START(reqperjxx(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjxx, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp1
    array(nn-nxg+1:nn,:,:) = array(nn-nxg+1:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Y- direction

    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'summation along Y'

    temp1  = 0.0_num
    temp2 = 0.0_num
    call MPI_START(reqperjxy(1),errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min_mpisp, tag, &
    comm, reqperjxy(2), errcode)
    call MPI_START(reqperjxy(3),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max_mpisp, tag, &
    comm, reqperjxy(4), errcode)
    CALL MPI_WAITALL(4_isp, reqperjxy, MPI_STATUSES_IGNORE, errcode)

    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp1
    array(:,nn-nyg+1:nn,:) = array(:,nn-nyg+1:nn,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)),temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'summation along Z'

    temp1  = 0.0_num
    temp2 = 0.0_num

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min_mpisp, tag, &
    comm, reqperjxz(2), errcode)
    call MPI_START(reqperjxz(1),errcode)
    !if (rank.eq.0) print*, 'irecv 1'

    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max_mpisp, tag, &
    comm, reqperjxz(4), errcode)
    call MPI_START(reqperjxz(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjxz, MPI_STATUSES_IGNORE, errcode)
    !if (rank.eq.0) print*, 'irecv 2'

    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp1
    array(:,:,nn-nzg+1:nn) = array(:,:,nn-nzg+1:nn) + temp2
    !if (rank.eq.0) print*, 'array sum'

    DEALLOCATE(temp1,temp2)

    !if (rank.eq.0) print*, 'deallocate'

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    !if (rank.eq.0) print*, 'field_bc'

    !call MPI_REQUEST_FREE(reqperx(1),errcode)
    !call MPI_REQUEST_FREE(reqperx(3),errcode)
    !call MPI_REQUEST_FREE(reqpery(1),errcode)
    !call MPI_REQUEST_FREE(reqpery(3),errcode)
    !call MPI_REQUEST_FREE(reqperz(1),errcode)
    !call MPI_REQUEST_FREE(reqperz(3),errcode)
  END SUBROUTINE

!!! --- Routine for adding current contributions fron adjacent subdomains persistent version
  SUBROUTINE summation_bcs_persistent_jy(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
    USE communications
    USE mpi

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarrayx, subarrayy, subarrayz, nn, sz
    INTEGER(isp) :: proc_x_min_mpisp, proc_x_max_mpisp
    INTEGER(isp) :: proc_y_min_mpisp, proc_y_max_mpisp
    INTEGER(isp) :: proc_z_min_mpisp, proc_z_max_mpisp

    proc_x_min_mpisp = proc_x_min
    proc_x_max_mpisp = proc_x_max
    proc_y_min_mpisp = proc_y_min
    proc_y_max_mpisp = proc_y_max
    proc_z_min_mpisp = proc_z_min
    proc_z_max_mpisp = proc_z_max

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    ! Initialization of the persistent communication
    IF (it .EQ. 0) THEN

      ! Init X
      subsizes(1) = nxg
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max_mpisp, tag, comm, reqperjyx(1), errcode)
      call MPI_SEND_INIT(array(-nxg+1,-nyg,-nzg), 1_isp, subarrayx, proc_x_min_mpisp, tag, comm, reqperjyx(3), errcode)

      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max_mpisp, tag, comm, reqperjyy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg+1,-nzg), 1_isp, subarrayy, proc_y_min_mpisp, tag, comm, reqperjyy(3), errcode)

      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max_mpisp, tag, comm, reqperjyz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg+1), 1_isp, subarrayz, proc_z_min_mpisp, tag, comm, reqperjyz(3), errcode)

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)
      CALL MPI_TYPE_FREE(subarrayz, errcode)

    ENDIF

    !! -- Summation along X- direction

    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'it',it
    !if (rank.eq.0) print*, 'summation along X'

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min_mpisp, tag, &
    comm, reqperjyx(2), errcode)
    call MPI_START(reqperjyx(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max_mpisp, tag, &
    comm, reqperjyx(4), errcode)
    call MPI_START(reqperjyx(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjyx, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp1
    array(nn-nxg+1:nn,:,:) = array(nn-nxg+1:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Y- direction

    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'summation along Y'

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min_mpisp, tag, &
    comm, reqperjyy(2), errcode)
    CALL MPI_START(reqperjyy(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max_mpisp, tag, &
    comm, reqperjyy(4), errcode)
    call MPI_START(reqperjyy(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjyy, MPI_STATUSES_IGNORE, errcode)

    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp1
    array(:,nn-nyg+1:nn,:) = array(:,nn-nyg+1:nn,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)),temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'summation along Z'

    temp1  = 0.0_num
    temp2 = 0.0_num

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min_mpisp, tag, &
    comm, reqperjyz(2), errcode)
    call MPI_START(reqperjyz(1),errcode)
    !if (rank.eq.0) print*, 'irecv 1'

    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max_mpisp, tag, &
    comm, reqperjyz(4), errcode)
    call MPI_START(reqperjyz(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjyz, MPI_STATUSES_IGNORE, errcode)
    !if (rank.eq.0) print*, 'irecv 2'

    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp1
    array(:,:,nn-nzg+1:nn) = array(:,:,nn-nzg+1:nn) + temp2
    !if (rank.eq.0) print*, 'array sum'

    DEALLOCATE(temp1,temp2)

    !if (rank.eq.0) print*, 'deallocate'

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    !if (rank.eq.0) print*, 'field_bc'

    !call MPI_REQUEST_FREE(reqperx(1),errcode)
    !call MPI_REQUEST_FREE(reqperx(3),errcode)
    !call MPI_REQUEST_FREE(reqpery(1),errcode)
    !call MPI_REQUEST_FREE(reqpery(3),errcode)
    !call MPI_REQUEST_FREE(reqperz(1),errcode)
    !call MPI_REQUEST_FREE(reqperz(3),errcode)

  END SUBROUTINE

!!! --- Routine for adding current contributions fron adjacent subdomains persistent version
  SUBROUTINE summation_bcs_persistent_jz(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
    USE communications
    USE mpi

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarrayx, subarrayy, subarrayz, nn, sz
    INTEGER(isp) :: proc_x_min_mpisp, proc_x_max_mpisp
    INTEGER(isp) :: proc_y_min_mpisp, proc_y_max_mpisp
    INTEGER(isp) :: proc_z_min_mpisp, proc_z_max_mpisp

    proc_x_min_mpisp = proc_x_min
    proc_x_max_mpisp = proc_x_max
    proc_y_min_mpisp = proc_y_min
    proc_y_max_mpisp = proc_y_max
    proc_z_min_mpisp = proc_z_min
    proc_z_max_mpisp = proc_z_max

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    ! Initialization of the persistent communication
    IF (it .EQ. 0) THEN

      ! Init X
      subsizes(1) = nxg
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max_mpisp, tag, comm, reqperjzx(1), errcode)
      call MPI_SEND_INIT(array(-nxg+1,-nyg,-nzg), 1_isp, subarrayx, proc_x_min_mpisp, tag, comm, reqperjzx(3), errcode)

      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max_mpisp, tag, comm, reqperjzy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg+1,-nzg), 1_isp, subarrayy, proc_y_min_mpisp, tag, comm, reqperjzy(3), errcode)

      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max_mpisp, tag, comm, reqperjzz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg+1), 1_isp, subarrayz, proc_z_min_mpisp, tag, comm, reqperjzz(3), errcode)

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)
      CALL MPI_TYPE_FREE(subarrayz, errcode)

    ENDIF

    !! -- Summation along X- direction

    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'it',it
    !if (rank.eq.0) print*, 'summation along X'

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min_mpisp, tag, &
    comm, reqperjzx(2), errcode)
    call MPI_START(reqperjzx(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max_mpisp, tag, &
    comm, reqperjzx(4), errcode)
    call MPI_START(reqperjzx(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjzx, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp1
    array(nn-nxg+1:nn,:,:) = array(nn-nxg+1:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Y- direction

    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'summation along Y'

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min_mpisp, tag, &
    comm, reqperjzy(2), errcode)
    CALL MPI_START(reqperjzy(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max_mpisp, tag, &
    comm, reqperjzy(4), errcode)
    call MPI_START(reqperjzy(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjzy, MPI_STATUSES_IGNORE, errcode)

    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp1
    array(:,nn-nyg+1:nn,:) = array(:,nn-nyg+1:nn,:) + temp2

    DEALLOCATE(temp1,temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)),temp2(subsizes(1), subsizes(2), subsizes(3)))

    !if (rank.eq.0) print*, 'summation along Z'

    temp1  = 0.0_num
    temp2 = 0.0_num

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min_mpisp, tag, &
    comm, reqperjzz(2), errcode)
    call MPI_START(reqperjzz(1),errcode)
    !if (rank.eq.0) print*, 'irecv 1'

    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max_mpisp, tag, &
    comm, reqperjzz(4), errcode)
    call MPI_START(reqperjzz(3),errcode)
    CALL MPI_WAITALL(4_isp, reqperjzz, MPI_STATUSES_IGNORE, errcode)
    !if (rank.eq.0) print*, 'irecv 2'

    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp1
    array(:,:,nn-nzg+1:nn) = array(:,:,nn-nzg+1:nn) + temp2
    !if (rank.eq.0) print*, 'array sum'

    DEALLOCATE(temp1,temp2)

    !if (rank.eq.0) print*, 'deallocate'

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    !if (rank.eq.0) print*, 'field_bc'

    !call MPI_REQUEST_FREE(reqperx(1),errcode)
    !call MPI_REQUEST_FREE(reqperx(3),errcode)
    !call MPI_REQUEST_FREE(reqpery(1),errcode)
    !call MPI_REQUEST_FREE(reqpery(3),errcode)
    !call MPI_REQUEST_FREE(reqperz(1),errcode)
    !call MPI_REQUEST_FREE(reqperz(3),errcode)

  END SUBROUTINE


!!! --- Boundary condition routine for electric field
  SUBROUTINE efield_bcs
  
    REAL(num) :: tmptime
#if defined(DEBUG)
  WRITE(0,*) "efield_bcs: start"
#endif

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Electric field MPI exchange between subdomains
    CALL field_bc(ex, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(ey, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(ez, nxguards, nyguards, nzguards, nx, ny, nz)
    IF (it.ge.timestat_itstart) THEN
    localtimes(8) = localtimes(8) + (MPI_WTIME() - tmptime)
    ENDIF
#if defined(DEBUG)
  WRITE(0,*) "efield_bcs: stop"
#endif
  END SUBROUTINE efield_bcs

!!! --- Boundary condition routine for magnetic field
  SUBROUTINE bfield_bcs
    REAL(num) :: tmptime
#if defined(DEBUG)
  WRITE(0,*) "bfield_bcs: start"
#endif
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Magnetic field MPI exchange between subdomains
    CALL field_bc(bx, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(by, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(bz, nxguards, nyguards, nzguards, nx, ny, nz)
    IF (it.ge.timestat_itstart) THEN
    localtimes(6) = localtimes(6) + (MPI_WTIME() - tmptime)
    ENDIF
#if defined(DEBUG)
  WRITE(0,*) "bfield_bcs: stop"
#endif
  END SUBROUTINE bfield_bcs

!!! --- Boundary conditions routine for currents
  SUBROUTINE current_bcs
    REAL(num) :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Add current contribution from adjacent subdomains

    IF (mpicom_curr.EQ.2) THEN

      CALL summation_bcs_persistent_jx(jx, nxjguards, nyjguards, nzjguards, nx, ny, nz)
      CALL summation_bcs_persistent_jy(jy, nxjguards, nyjguards, nzjguards, nx, ny, nz)
      CALL summation_bcs_persistent_jz(jz, nxjguards, nyjguards, nzjguards, nx, ny, nz)

    ELSE IF (mpicom_curr.EQ.1) THEN

      CALL summation_bcs(jx, nxjguards, nyjguards, nzjguards, nx, ny, nz)
      CALL summation_bcs(jy, nxjguards, nyjguards, nzjguards, nx, ny, nz)
      CALL summation_bcs(jz, nxjguards, nyjguards, nzjguards, nx, ny, nz)

    ELSE

      CALL summation_bcs_nonblocking(jx, nxjguards, nyjguards, nzjguards, nx, ny, nz)
      CALL summation_bcs_nonblocking(jy, nxjguards, nyjguards, nzjguards, nx, ny, nz)
      CALL summation_bcs_nonblocking(jz, nxjguards, nyjguards, nzjguards, nx, ny, nz)

    ENDIF
    IF (it.ge.timestat_itstart) THEN
    localtimes(4) = localtimes(4) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE current_bcs

!!! --- Boundary conditions routine for charge density
SUBROUTINE charge_bcs
! Add charge contribution from adjacent subdomains
! __________________________________________________
    USE time_stat


    REAL(num) :: tmptime
    IF (it.ge.timestat_itstart) THEN
    tmptime = MPI_WTIME()
    ENDIF

    IF (mpicom_curr.EQ.1) THEN
    	CALL summation_bcs(rho, nxjguards, nyjguards, nzjguards, nx, ny, nz)
    ELSE
    	CALL summation_bcs_nonblocking(rho, nxjguards, nyjguards, nzjguards, nx, ny, nz)
    ENDIF
    IF (it.ge.timestat_itstart) THEN
    localtimes(13) = localtimes(13) + (MPI_WTIME() - tmptime)
    ENDIF

END SUBROUTINE charge_bcs


! ________________________________________________________________________________________
!> Boundary condition routine on particles in 3d
!> @brief
!
!> This subroutine is the main one to manage the particle boundary conditions
!> between MPI domains and between tiles.
!> The different model of exchanges can be used via the variable partcom
!> @details
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> last modified: 09/12/2016
  SUBROUTINE particle_bcs
! ________________________________________________________________________________________
	USE omp_lib
	USE time_stat
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
	WRITE(0,*) "particle_bcs_tiles: start"
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
	WRITE(0,*) "particle_bcs_mpi: start"
#endif

	IF (mpicom_curr .EQ. 1) THEN
		! Then exchange particle between MPI domains
		CALL particle_bcs_mpi_blocking()
	ELSE
		CALL particle_bcs_mpi_non_blocking()
	ENDIF

#if defined(DEBUG)
	WRITE(0,*) "particle_bcs_mpi: stop"
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
  WRITE(0,*) "particle_bcs_tiles_openmp: start"
#else
  WRITE(0,*) "particle_bcs_tiles: start"
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
  WRITE(0,*) "particle_bcs_tiles: stop"
#endif

    IF (it.ge.timestat_itstart) THEN
      localtimes(11) = localtimes(11) + (MPI_WTIME() - tmptime)
      tmptime = MPI_WTIME()
    ENDIF
    tend = MPI_WTIME()
    local_time_part=local_time_part+(tend-tdeb)

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_mpi: start"
#endif
    IF (mpicom_curr .EQ. 1) THEN
    	! Then exchange particle between MPI domains
    	CALL particle_bcs_mpi_blocking()
    ELSE
    	CALL particle_bcs_mpi_non_blocking()
    ENDIF

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_mpi: stop"
#endif

    IF (it.ge.timestat_itstart) THEN
      localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)
    ENDIF

    ! _____________________________________________
    ! Tile and MPI com in one
    ! Default
    ELSE

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_tiles_and_mpi: start"
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
  WRITE(0,*) "particle_bcs_tiles_and_mpi: stop"
#endif
    ENDIF

  END SUBROUTINE particle_bcs

!!! Boundary condition routine on particles
  SUBROUTINE particle_bcs_2d
  	USE omp_lib
    USE time_stat
    USE tiling

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
  WRITE(0,*) "particle_bcs_tiles_openmp: start"
#else
  WRITE(0,*) "particle_bcs_tiles: start"
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
  WRITE(0,*) "particle_bcs_mpi: start"
#endif

    	CALL particle_bcs_mpi_non_blocking_2d()

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_mpi: stop"
#endif

    localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)

    ! _____________________________________________
    ! Tile and MPI com in one
    ! Default
    ELSE

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_tiles_and_mpi: start"
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
  WRITE(0,*) "particle_bcs_tiles_and_mpi: stop"
#endif
    ENDIF

  END SUBROUTINE particle_bcs_2d

!!! Boundary condition on tiles
  SUBROUTINE particle_bcs_tiles
  
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv

    DO ispecies=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(ispecies)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO iz=1, ntilez! LOOP ON TILES
            DO iy=1, ntiley
                DO ix=1, ntilex
                    curr_tile=>curr%array_of_tiles(ix,iy,iz)
                    nptile=curr_tile%np_tile(1)
                    DO i=nptile, 1, -1! LOOP ON PARTICLES
                        partx=curr_tile%part_x(i)
                        party=curr_tile%part_y(i)
                        partz=curr_tile%part_z(i)
                        partux=curr_tile%part_ux(i)
                        partuy=curr_tile%part_uy(i)
                        partuz=curr_tile%part_uz(i)
                        gaminv=curr_tile%part_gaminv(i)
                        partw=curr_tile%pid(i,wpid)

                        ! Case 1: if particle did not leave tile nothing to do
                        IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    			   &
                        .AND. ((party .GE. curr_tile%y_tile_min) .AND. (party .LT. curr_tile%y_tile_max))  			   &
                        .AND. ((partz .GE. curr_tile%z_tile_min+zgrid) .AND. (partz .LT. curr_tile%z_tile_max+zgrid))) &
                        CYCLE

                        ! Case 2: if particle left MPI domain nothing to do now
                        IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) 			CYCLE
                        IF ((party .LT. y_min_local) .OR. (party .GE. y_max_local)) 			CYCLE
                        IF ((partz .LT. z_min_local+zgrid) .OR. (partz .GE. z_max_local+zgrid)) CYCLE

                        ! Case 3: particles changed tile. Tranfer particle to new tile
                        ! Get new indexes of particle in array of tiles
                        indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),idp)+1,ntilex)
                        indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy),idp)+1,ntiley)
                        indz = MIN(FLOOR((partz-(z_min_local+zgrid)+dz/2_num)/(nz0_grid_tile*dz),idp)+1,ntilez)
                        CALL rm_particle_at_tile(curr,ix,iy,iz,i)
                        CALL add_particle_at_tile(curr, indx,indy,indz, &
                             partx, party, partz, partux, partuy, partuz, gaminv, partw)
                    END DO !END LOOP ON PARTICLES
                END DO
            END DO
        END DO ! END LOOP ON TILES
    END DO ! END LOOP ON SPECIES
  END SUBROUTINE particle_bcs_tiles

! ________________________________________________________________________________________
!> Boundary condition on tiles - 3D version
!> @brief
!
!> This version is efficient when the number of tiles is large
!> compared to the number of threads
!
!> @author
!> Henri Vincenti
!
!> @date
!> 2016
	SUBROUTINE particle_bcs_tiles_openmp()
! ________________________________________________________________________________________

		USE omp_lib
		IMPLICIT NONE
		INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
		INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
		TYPE(particle_species), POINTER :: curr
		TYPE(particle_tile), POINTER :: curr_tile
		REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
		INTEGER(idp) :: nthreads_tot, nthreads_loop1, nthreads_loop2

#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	!nthreads_tot=OMP_GET_NUM_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif

		IF (nthreads_tot .GT. 1) THEN
			nthreads_loop1=MIN(nspecies,nthreads_tot)
			nthreads_loop2=MAX(1_idp,nthreads_tot/nthreads_loop1)
		ELSE
			nthreads_loop1=1
			nthreads_loop2=1
		ENDIF

	!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,ispecies, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz, &
	!$OMP partx,party,partz,partux,partuy,partuz,gaminv,partw,indx,indy,indz,nptile,curr_tile) &
	!$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
    !$OMP x_max_local,y_max_local,z_max_local,zgrid,dx,dy,dz) NUM_THREADS(nthreads_loop1)
    DO ispecies=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(ispecies)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO ipz=1,3
        	DO ipy=1,3
        		DO ipx=1,3
					!$OMP PARALLEL DO DEFAULT(NONE) SHARED(curr,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
					!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,zgrid)  &
					!$OMP FIRSTPRIVATE(ipx,ipy,ipz) &
					!$OMP PRIVATE(ix,iy,iz,i,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
					!$OMP indx,indy,indz) COLLAPSE(3) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
					DO iz=ipz, ntilez,3! LOOP ON TILES
						DO iy=ipy, ntiley,3
							DO ix=ipx, ntilex,3
								curr_tile=>curr%array_of_tiles(ix,iy,iz)
								nptile=curr_tile%np_tile(1)
								DO i=nptile, 1, -1! LOOP ON PARTICLES
									partx=curr_tile%part_x(i)
									party=curr_tile%part_y(i)
									partz=curr_tile%part_z(i)
									partux=curr_tile%part_ux(i)
									partuy=curr_tile%part_uy(i)
									partuz=curr_tile%part_uz(i)
									gaminv=curr_tile%part_gaminv(i)
									partw=curr_tile%pid(i,wpid)

									! Case 1: if particle did not leave tile nothing to do
									IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    &
									.AND. ((party .GE. curr_tile%y_tile_min) .AND. (party .LT. curr_tile%y_tile_max))  &
									.AND. ((partz .GE. curr_tile%z_tile_min+zgrid) .AND. (partz .LT. curr_tile%z_tile_max+zgrid))) &
									CYCLE

									! Case 2: if particle left MPI domain nothing to do now
									IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
									IF ((party .LT. y_min_local) .OR. (party .GE. y_max_local)) CYCLE
									IF ((partz .LT. z_min_local+zgrid) .OR. (partz .GE. z_max_local+zgrid)) CYCLE

									! Case 3: particles changed tile. Tranfer particle to new tile
									! Get new indexes of particle in array of tiles
									indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),idp)+1,ntilex)
									indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy),idp)+1,ntiley)
									indz = MIN(FLOOR((partz-(z_min_local+zgrid)+dz/2_num)/(nz0_grid_tile*dz),idp)+1,ntilez)
									CALL rm_particle_at_tile(curr,ix,iy,iz,i)
									CALL add_particle_at_tile(curr, indx,indy,indz, &
										 partx, party, partz, partux, partuy, partuz, gaminv, partw)
								END DO !END LOOP ON PARTICLES
							END DO
						END DO
					END DO ! END LOOP ON TILES
					!$OMP END PARALLEL DO
        		END DO
        	END DO
        END DO
    END DO ! END LOOP ON SPECIES
    !$OMP END PARALLEL DO
  END SUBROUTINE particle_bcs_tiles_openmp

! ________________________________________________________________________________________
!> @brief
!> Boundary condition on tiles
  SUBROUTINE particle_bcs_tiles_2d
! ________________________________________________________________________________________
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indz
    INTEGER(idp) :: nptile, nx0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv

    iy=1
    DO ispecies=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(ispecies)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO iz=1, ntilez! LOOP ON TILES
			DO ix=1, ntilex
				curr_tile=>curr%array_of_tiles(ix,iy,iz)
				nptile=curr_tile%np_tile(1)
				DO i=nptile, 1, -1! LOOP ON PARTICLES
					partx=curr_tile%part_x(i)
					party=curr_tile%part_y(i)
					partz=curr_tile%part_z(i)
					partux=curr_tile%part_ux(i)
					partuy=curr_tile%part_uy(i)
					partuz=curr_tile%part_uz(i)
					gaminv=curr_tile%part_gaminv(i)
					partw=curr_tile%pid(i,wpid)

					! Case 1: if particle did not leave tile nothing to do
					IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    &
					.AND. ((partz .GE. curr_tile%z_tile_min+zgrid) .AND. (partz .LT. curr_tile%z_tile_max+zgrid))) &
					CYCLE

					! Case 2: if particle left MPI domain nothing to do now
					IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
					IF ((partz .LT. z_min_local+zgrid) .OR. (partz .GE. z_max_local+zgrid)) CYCLE

					! Case 3: particles changed tile. Tranfer particle to new tile
					! Get new indexes of particle in array of tiles
					indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),idp)+1,ntilex)
					indz = MIN(FLOOR((partz-(z_min_local+zgrid)+dz/2_num)/(nz0_grid_tile*dz),idp)+1,ntilez)
					CALL rm_particle_at_tile(curr,ix,iy,iz,i)
					CALL add_particle_at_tile(curr, indx,iy,indz, &
						 partx, party, partz, partux, partuy, partuz, gaminv, partw)
				END DO !END LOOP ON PARTICLES
			END DO
        END DO ! END LOOP ON TILES
    END DO ! END LOOP ON SPECIES
  END SUBROUTINE particle_bcs_tiles_2d

! ________________________________________________________________________________________
!> Boundary condition on tiles - 2D version
!> This version is efficient when the number of tiles is large
!> compared to the number of threads
  SUBROUTINE particle_bcs_tiles_2d_openmp()
! ________________________________________________________________________________________

    USE omp_lib
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz, ipx, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: nthreads_tot, nthreads_loop1, nthreads_loop2

#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif

	IF (nthreads_tot .GT. 1) THEN
		nthreads_loop1=MIN(nspecies,nthreads_tot)
		nthreads_loop2=MAX(1_idp,nthreads_tot/nthreads_loop1)
	ELSE
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF
	iy=1
	!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,ispecies, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipz, &
	!$OMP partx,party,partz,partux,partuy,partuz,gaminv,partw,indx,indy,indz,curr_tile,nptile)    &
	!$OMP SHARED(iy,nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
    !$OMP x_max_local,y_max_local,z_max_local,zgrid,dx,dy,dz) NUM_THREADS(nthreads_loop1)
    DO ispecies=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(ispecies)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO ipz=1,3
        	DO ipx=1,3
				!$OMP PARALLEL DO DEFAULT(NONE) SHARED(iy,curr,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
				!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,zgrid)  &
				!$OMP FIRSTPRIVATE(ipx,ipz) &
				!$OMP PRIVATE(ix,iz,i,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
				!$OMP indx,indy,indz) COLLAPSE(2) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
				DO iz=ipz, ntilez,3! LOOP ON TILES
					DO ix=ipx, ntilex,3
						curr_tile=>curr%array_of_tiles(ix,iy,iz)
						nptile=curr_tile%np_tile(1)
						DO i=nptile, 1, -1! LOOP ON PARTICLES
							partx=curr_tile%part_x(i)
							!party=curr_tile%part_y(i)
							partz=curr_tile%part_z(i)
							partux=curr_tile%part_ux(i)
							partuy=curr_tile%part_uy(i)
							partuz=curr_tile%part_uz(i)
							gaminv=curr_tile%part_gaminv(i)
							partw=curr_tile%pid(i,wpid)

							! Case 1: if particle did not leave tile nothing to do
							IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    &
							.AND. ((partz .GE. curr_tile%z_tile_min+zgrid) .AND. (partz .LT. curr_tile%z_tile_max+zgrid))) &
							CYCLE

							! Case 2: if particle left MPI domain nothing to do now
							IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
							IF ((partz .LT. z_min_local+zgrid) .OR. (partz .GE. z_max_local+zgrid)) CYCLE

							! Case 3: particles changed tile. Tranfer particle to new tile
							! Get new indexes of particle in array of tiles
							indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),idp)+1,ntilex)
							indz = MIN(FLOOR((partz-(z_min_local+zgrid)+dz/2_num)/(nz0_grid_tile*dz),idp)+1,ntilez)
							!CALL rm_particle_at_tile_2d(curr,ix,iy,iz,i)
							CALL rm_particle_at_tile_2d(curr,ix,iz,i)
							!CALL add_particle_at_tile(curr, indx,iy,indz, &
							!	 partx, party, partz, partux, partuy, partuz, gaminv, partw)
							CALL add_particle_at_tile_2d(curr, indx,indz, &
								 partx, partz, partux, partuy, partuz, gaminv, partw)
						END DO !END LOOP ON PARTICLES
					END DO
				END DO ! END LOOP ON TILES
				!$OMP END PARALLEL DO
        	END DO
        END DO
    END DO ! END LOOP ON SPECIES
    !$OMP END PARALLEL DO
  END SUBROUTINE particle_bcs_tiles_2d_openmp

! Experimental subroutine
#if defined(DEV)
  ! __________________________________________________________________
  SUBROUTINE particle_bsc_openmp_reordering
  ! ********** EXPERIMENTAL subroutine ****************************
  ! => should not be used
  ! This subroutine process the particle boundary conditions
  ! For this aim, the particles are reordered into buckets in the particle arrays
  ! Each bucket corresponds to the group of particle to be exchanged in a common direction
  ! Exchanges between tiles and between MPI domains is combined for better efficiency
  ! This subroutine is less efficient that particle_bcs_tiles_openmp
  ! Mathieu Lobet, 2016
  ! __________________________________________________________________

    USE omp_lib
    USE communications
    IMPLICIT NONE

    INTEGER(idp) :: i, is, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER    :: curr_tile, curr_tile_add
    REAL(num)    :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: test =0, nthreads_tot, nthreads_loop1, nthreads_loop2
    INTEGER(idp) :: dirx,diry,dirz
    INTEGER(idp) :: k,ib
    INTEGER(idp) :: ipmin,ipmax
    type(part_com_buffer), ALLOCATABLE, DIMENSION(:,:,:,:) :: buffer


#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif


	IF (nthreads_tot .GT. 1) THEN
		nthreads_loop1=MIN(nspecies,nthreads_tot)
		nthreads_loop2=MAX(1_idp,nthreads_tot/nthreads_loop1)
	ELSE
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF

	ALLOCATE(buffer(ntilex,ntiley,ntilez,nspecies))

	  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,is, nx0_grid_tile,ny0_grid_tile,&
	  !$OMP nz0_grid_tile,ipx,ipy,ipz,indx,indy,indz,partw,ib,k,dirx,diry,dirz,&
	  !$OMP partx, party, partz, curr_tile,nptile,partux,partuy,partuz) &
	  !$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
    !$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz,buffer) NUM_THREADS(nthreads_loop1)
    DO is=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(is)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO ipz=1,3
        	DO ipy=1,3
        		DO ipx=1,3
					!$OMP PARALLEL DO DEFAULT(NONE) SHARED(curr,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
					!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, &
					!$OMP buffer, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile)  &
					!$OMP FIRSTPRIVATE(ipx,ipy,ipz,is) &
					!$OMP PRIVATE(ix,iy,iz,i,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
					!$OMP indx,indy,indz,dirx,diry,dirz,ib,k) &
					!$OMP COLLAPSE(3) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
					DO iz=ipz, ntilez,3! LOOP ON TILES
						DO iy=ipy, ntiley,3
							DO ix=ipx, ntilex,3
								curr_tile=>curr%array_of_tiles(ix,iy,iz)
								nptile=curr_tile%np_tile(1)

								! Temporary array
								ALLOCATE(buffer(ix,iy,iz,is)%part_x(nptile), &
								         buffer(ix,iy,iz,is)%part_y(nptile), &
								         buffer(ix,iy,iz,is)%part_z(nptile), &
								         buffer(ix,iy,iz,is)%part_ux(nptile), &
								         buffer(ix,iy,iz,is)%part_uy(nptile), &
								         buffer(ix,iy,iz,is)%part_uz(nptile),&
								         buffer(ix,iy,iz,is)%part_gaminv(nptile),&
								         buffer(ix,iy,iz,is)%pid(nptile,1), &
								         buffer(ix,iy,iz,is)%boundid(nptile), &
								         buffer(ix,iy,iz,is)%bin_npart(0:27), &
								         buffer(ix,iy,iz,is)%bin_pos(0:27))

								 buffer(ix,iy,iz,is)%bin_npart(:) = 0

								!buffer(ix,iy,iz,is)%part_x(1:nptile) = curr_tile%part_x(1:nptile)
								!buffer(ix,iy,iz,is)%part_y(1:nptile) = curr_tile%part_y(1:nptile)
								!buffer%part_z(1:nptile) = curr_tile%part_z(1:nptile)
								!buffer%part_ux(1:nptile) = curr_tile%part_ux(1:nptile)
								!buffer%part_uy(1:nptile) = curr_tile%part_uy(1:nptile)
								!buffer%part_uz(1:nptile) = curr_tile%part_uz(1:nptile)
								!buffer%part_gaminv(1:nptile) = curr_tile%part_gaminv(1:nptile)
								!buffer%pid(1:nptile,1) = curr_tile%pid(1:nptile,wpid)

								! Particle reordering
								! Step 1 - determine the number of particles in each bin
								DO i=1,nptile
									partx=curr_tile%part_x(i)
									party=curr_tile%part_y(i)
									partz=curr_tile%part_z(i)
									partux=curr_tile%part_ux(i)
									partuy=curr_tile%part_uy(i)
									partuz=curr_tile%part_uz(i)
									gaminv=curr_tile%part_gaminv(i)
									partw=curr_tile%pid(i,wpid)


									! Case 1: if particle did not leave tile nothing to do
									IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    &
									.AND. ((party .GE. curr_tile%y_tile_min) .AND. (party .LT. curr_tile%y_tile_max))  &
									.AND. ((partz .GE. curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) &
									THEN

										buffer(ix,iy,iz,is)%bin_npart(0) = &
								    buffer(ix,iy,iz,is)%bin_npart(0) + 1

								    buffer(ix,iy,iz,is)%boundid(i) = 0

								    CYCLE

									! Case 2: if particle left MPI domain nothing to do now
                  ELSE IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local).AND. &
									         (party .LT. y_min_local) .OR. (party .GE. y_max_local).AND. &
									         (partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) THEN

										buffer(ix,iy,iz,is)%bin_npart(0) = &
								    buffer(ix,iy,iz,is)%bin_npart(0) + 1

								    buffer(ix,iy,iz,is)%boundid(i) = 0

								    CYCLE

									ENDIF

									indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),idp)+1,ntilex)
									indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy),idp)+1,ntiley)
									indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz),idp)+1,ntilez)

                  ! Direction of the particle
									dirx = indx - ix
									diry = indy - iy
									dirz = indz - iz

								  buffer(ix,iy,iz,is)%boundid(i) = (2+dirx) + (1+diry)*3 + (1+dirz)*9

								    !if (buffer(ix,iy,iz,is)%boundid(i)>27) then
								    !print*, indx,indy,indz
								    !print*, dirx,diry,dirz
								    !end if

								  buffer(ix,iy,iz,is)%bin_npart(buffer(ix,iy,iz,is)%boundid(i)) = &
								  buffer(ix,iy,iz,is)%bin_npart(buffer(ix,iy,iz,is)%boundid(i)) + 1

								ENDDO

								! Particle reordering
								! Step 2 - Determine bin positions
								buffer(ix,iy,iz,is)%bin_pos(0) = 1
								Do i=1,27


								  buffer(ix,iy,iz,is)%bin_pos(i) = buffer(ix,iy,iz,is)%bin_pos(i-1) &
								  + buffer(ix,iy,iz,is)%bin_npart(i-1)
								ENDDO

								! Particle reordering
								! Step 3 - reorder the particles
								DO i=1,nptile
								  ib = buffer(ix,iy,iz,is)%boundid(i)
								  k = buffer(ix,iy,iz,is)%bin_pos(ib)

                  buffer(ix,iy,iz,is)%part_x(k) = curr_tile%part_x(i)
                  buffer(ix,iy,iz,is)%part_y(k) = curr_tile%part_y(i)
                  buffer(ix,iy,iz,is)%part_z(k) = curr_tile%part_z(i)
                  buffer(ix,iy,iz,is)%part_ux(k) = curr_tile%part_ux(i)
                  buffer(ix,iy,iz,is)%part_uy(k) = curr_tile%part_uy(i)
                  buffer(ix,iy,iz,is)%part_uz(k) = curr_tile%part_uz(i)
                  buffer(ix,iy,iz,is)%part_gaminv(k) = curr_tile%part_gaminv(i)
                  buffer(ix,iy,iz,is)%pid(k,1) = curr_tile%pid(i,wpid)

                  buffer(ix,iy,iz,is)%bin_pos(ib) = &
                  buffer(ix,iy,iz,is)%bin_pos(ib) + 1
								ENDDO

								buffer(ix,iy,iz,is)%bin_pos(0) = 1
								Do i=1,27
								  buffer(ix,iy,iz,is)%bin_pos(i) = buffer(ix,iy,iz,is)%bin_pos(i-1) &
								  + buffer(ix,iy,iz,is)%bin_npart(i-1)
								ENDDO

                ! Particle that will stay in the domain
                curr_tile%part_x(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_x(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%part_y(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_y(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%part_z(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_z(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%part_ux(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_ux(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%part_uy(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_uy(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%part_uz(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_uz(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%part_gaminv(1:buffer(ix,iy,iz,is)%bin_npart(0)) = &
                buffer(ix,iy,iz,is)%part_gaminv(1:buffer(ix,iy,iz,is)%bin_npart(0))
                curr_tile%pid(1:buffer(ix,iy,iz,is)%bin_npart(0),wpid) = &
                buffer(ix,iy,iz,is)%pid(1:buffer(ix,iy,iz,is)%bin_npart(0),wpid)

                curr_tile%np_tile(1) = buffer(ix,iy,iz,is)%bin_npart(0)

							END DO
						END DO
					END DO ! END LOOP ON TILES
					!$OMP END PARALLEL DO
        		END DO
        	END DO
        END DO
    END DO ! END LOOP ON SPECIES
    !$OMP END PARALLEL DO


   ! ___________________________________________________
   ! Exchange between tiles
   ! Reordered buffer are copied in parallel

	  !$OMP PARALLEL DO DEFAULT(NONE) &
	  !$OMP PRIVATE(curr,is, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,&
	  !$OMP         ipx,ipy,ipz,ib,indx,indy,indz,nptile,ipmin,ipmax,curr_tile) &
	  !$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
    !$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz,buffer) NUM_THREADS(nthreads_loop1)
    DO is=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(is)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO ipz=1,3
        	DO ipy=1,3
        		DO ipx=1,3
					!$OMP PARALLEL DO DEFAULT(NONE) SHARED(curr,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
					!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, &
					!$OMP buffer, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile)  &
					!$OMP FIRSTPRIVATE(ipx,ipy,ipz,is) &
					!$OMP PRIVATE(ix,iy,iz,i,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
					!$OMP indx,indy,indz,ipmin,ipmax,ib) COLLAPSE(3) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
					DO iz=ipz, ntilez,3! LOOP ON TILES
						DO iy=ipy, ntiley,3
							DO ix=ipx, ntilex,3
								curr_tile=>curr%array_of_tiles(ix,iy,iz)
								nptile=curr_tile%np_tile(1)

								! Copy of the buffers in every directions
                DO dirz = -1, 1
                  DO diry = -1, 1
                    DO dirx = -1, 1

                     indx = ix+dirx
                     indy = iy+diry
                     indz = iz+dirz

                     IF((indx>=1).and.(indx<=ntilex).and. &
                        (indy>=1).and.(indy<=ntiley).and. &
                        (indz>=1).and.(indz<=ntilez)) THEN

                       ib = 2+dirx + (1+diry)*3 + (1+dirz)*9

                       ipmin = buffer(ix,iy,iz,is)%bin_pos(ib)
                       ipmax = ipmin+buffer(ix,iy,iz,is)%bin_npart(ib)-1

									     CALL add_group_of_particles_at_tile(curr,indx,indy,indz, &
									     buffer(ix,iy,iz,is)%bin_npart(ib),&
										   buffer(ix,iy,iz,is)%part_x(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%part_y(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%part_z(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%part_ux(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%part_uy(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%part_uz(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%part_gaminv(ipmin:ipmax), &
										   buffer(ix,iy,iz,is)%pid(ipmin:ipmax,1))

										  ENDIF

                    ENDDO
                  ENDDO
                ENDDO

							END DO
						END DO
					END DO ! END LOOP ON TILES
					!$OMP END PARALLEL DO
        		END DO
        	END DO
        END DO
    END DO ! END LOOP ON SPECIES
    !$OMP END PARALLEL DO

  END SUBROUTINE particle_bsc_openmp_reordering
#endif

! ________________________________________________________________________________________
!> @brief
!> MPI Boundary condition routine on particles
  SUBROUTINE particle_bcs_mpi_blocking
! ________________________________________________________________________________________

    USE mpi

    INTEGER(isp), PARAMETER :: nvar=8 ! Simple implementation
    INTEGER(isp), DIMENSION(-1:1,-1:1,-1:1) :: nptoexch
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: recvbuf
    LOGICAL(lp) , ALLOCATABLE, DIMENSION(:) :: mask
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

    DO ispecies=1, nspecies !LOOP ON SPECIES
        tdeb=MPI_WTIME()
        ! Init send recv buffers
        currsp => species_parray(ispecies)
        nptoexch=0
        nsend_buf=0
        nout=0
        nrecv_buf=0
        nbuff=0
        ibuff=1
        ! GET NUMBER OF PARTICLES IN BORDER TILES
        DO iztile=1, ntilez !LOOP ON TILES
            DO iytile=1, ntiley
                DO ixtile=1, ntilex
                    curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
                    IF (.NOT. curr%subdomain_bound) CYCLE
                    	nbuff=nbuff+curr%np_tile(1)
            	END DO
            END DO
        END DO
        ALLOCATE(sendbuf(-1:1,-1:1,-1:1,1:nbuff*nvar))
        DO iztile=1, ntilez !LOOP ON TILES
            DO iytile=1, ntiley
                DO ixtile=1, ntilex
                    curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
                    ! If not subdomain border, nothing to do
                    IF (.NOT. curr%subdomain_bound) CYCLE
                    ! Else, search for outbound particles
                    ALLOCATE(mask(1:curr%np_tile(1)))
                    mask=.TRUE.
                    part_xyz=0.
                    ! Identify outbounds particles
                    DO i = 1, curr%np_tile(1) !LOOP ON PARTICLES
                        xbd = 0
                        ybd = 0
                        zbd = 0
                        out_of_bounds = .FALSE.
                        part_xyz = curr%part_x(i)
                        ! Particle has left this processor -x
                        IF (part_xyz .LT. x_min_local) THEN
                            xbd = -1
                            IF (x_min_boundary) THEN
                            	SELECT CASE (pbound_x_min)
                            	CASE (1_idp) ! absorbing
                            		mask(i)=.FALSE.
                            		CYCLE
                            	CASE (2_idp) ! Reflecting
                            	    curr%part_x(i) = part_xyz + dx
									                curr%part_ux(i) = - curr%part_ux(i)
                                  xbd=0
                            	CASE DEFAULT ! periodic
                                	curr%part_x(i) = part_xyz + length_x
                                END SELECT
                            ENDIF
                        ENDIF

                        ! Particle has left this processor +x
                        IF (part_xyz .GE. x_max_local) THEN
                            xbd = 1
                            IF (x_max_boundary) THEN
                            	SELECT CASE (pbound_x_max)
                            	CASE (1_idp) ! absorbing
                            		mask(i)=.FALSE.
                            		CYCLE
                            	CASE (2_idp) ! Reflecting
                            	    curr%part_x(i) = part_xyz - dx
									                curr%part_ux(i) = - curr%part_ux(i)
                                  xbd=0
                            	CASE DEFAULT ! periodic
                                	curr%part_x(i) = part_xyz - length_x
                                END SELECT
                            ENDIF
                        ENDIF

                        part_xyz = curr%part_y(i)
                        ! Particle has left this processor -y
                        IF ((part_xyz .LT. y_min_local) .AND. (c_dim .EQ. 3)) THEN
                            ybd = -1
                            IF (y_min_boundary) THEN
                            	SELECT CASE (pbound_y_min)! absorbing
                            	CASE (1_idp)
                            		mask(i)=.FALSE.
                            		CYCLE
                            	CASE (2_idp) ! Reflecting
                            	    curr%part_y(i) = part_xyz + dy
									                curr%part_uy(i) = - curr%part_uy(i)
                                  ybd=0
                            	CASE DEFAULT ! periodic
                               		curr%part_y(i) = part_xyz + length_y
                                END SELECT
                            ENDIF
                        ENDIF

                        ! Particle has left this processor
                        IF ((part_xyz .GE. y_max_local) .AND. (c_dim .EQ. 3)) THEN
                            ybd = 1
                            IF (y_max_boundary) THEN
                            	SELECT CASE (pbound_y_max)
                            	CASE (1_idp) ! absorbing
                            		mask(i)=.FALSE.
                            		CYCLE
                            	CASE (2_idp) ! Reflecting
                            	    curr%part_y(i) = part_xyz - dy
									                curr%part_uy(i) = - curr%part_uy(i)
									                ybd=0
                            	CASE DEFAULT ! periodic
                                	curr%part_y(i) = part_xyz - length_y
                                END SELECT
                            ENDIF
                        ENDIF

                        part_xyz = curr%part_z(i)
                        ! Particle has left this processor
                        IF (part_xyz .LT. z_min_local+zgrid) THEN
                            zbd = -1
                            IF (z_min_boundary) THEN
								              SELECT CASE (pbound_z_min)
								              CASE (1_idp) ! absorbing
									              mask(i)=.FALSE.
									              CYCLE
                              CASE (2_idp) ! Reflecting
                            	  curr%part_z(i) = part_xyz + dz
									              curr%part_uz(i) = - curr%part_uz(i)
									              zbd=0
								              CASE DEFAULT ! periodic
									              curr%part_z(i) = part_xyz + length_z
								              END SELECT
                            ENDIF
                        ENDIF

                        ! Particle has left this processor
                        IF (part_xyz .GE. z_max_local+zgrid) THEN
                            zbd = 1
                            ! Particle has left the system
                            IF (z_max_boundary) THEN
                            	SELECT CASE (pbound_z_max)
                            	CASE (1_idp) ! absorbing
                            		mask(i)=.FALSE.
                            		CYCLE
                            	CASE (2_idp) ! Reflecting
                            	    curr%part_z(i) = part_xyz - dz
									                curr%part_uz(i) = - curr%part_uz(i)
                                  zbd=0
                            	CASE DEFAULT ! periodic
                                	curr%part_z(i) = part_xyz - length_z
                                END SELECT
                            ENDIF
                        ENDIF

                        IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
                        ! Particle has left processor, send it to its neighbour
                            mask(i)=.FALSE.
                            nout=nout+1
                            ibuff=nptoexch(xbd,ybd,zbd)*nvar+1
                            sendbuf(xbd,ybd,zbd,ibuff)    = curr%part_x(i)
                            sendbuf(xbd,ybd,zbd,ibuff+1)  = curr%part_y(i)
                            sendbuf(xbd,ybd,zbd,ibuff+2)  = curr%part_z(i)
                            sendbuf(xbd,ybd,zbd,ibuff+3)  = curr%part_ux(i)
                            sendbuf(xbd,ybd,zbd,ibuff+4)  = curr%part_uy(i)
                            sendbuf(xbd,ybd,zbd,ibuff+5)  = curr%part_uz(i)
                            sendbuf(xbd,ybd,zbd,ibuff+6)  = curr%part_gaminv(i)
                            sendbuf(xbd,ybd,zbd,ibuff+7)  = curr%pid(i,wpid)
                            nptoexch(xbd,ybd,zbd) = nptoexch(xbd,ybd,zbd)+1
                        ENDIF
                    ENDDO !END LOOP ON PARTICLES
                    ! Remove outbound particles from current tile
                    CALL rm_particles_from_species_with_mask(currsp, ixtile,iytile,iztile, mask)
                    DEALLOCATE(mask)
                  ENDDO
               ENDDO
            ENDDO ! END LOOP ON TILES
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
                            nsend_buf=nptoexch(ix,iy,iz)*nvar
                            nrecv_buf=0
                            dest = INT(neighbour(ix,iy,iz),isp)
                            src  = INT(neighbour(ixp,iyp,izp),isp)
                            CALL MPI_SENDRECV(nsend_buf, 1_isp, MPI_INTEGER, dest, tag, nrecv_buf, 1_isp, &
                            MPI_INTEGER, src, tag, comm, status, errcode)
                            ALLOCATE(recvbuf(1:nrecv_buf))
                            CALL MPI_SENDRECV(sendbuf(ix,iy,iz,1:nsend_buf), nsend_buf, mpidbl, dest, tag, &
                            recvbuf, nrecv_buf, mpidbl, src, tag, comm, status, errcode)
                            ! Add received particles to particle arrays
                            DO i =1, nrecv_buf, nvar
                                CALL add_particle_to_species(currsp, recvbuf(i), recvbuf(i+1), recvbuf(i+2), &
                                recvbuf(i+3), recvbuf(i+4), recvbuf(i+5), recvbuf(i+6),recvbuf(i+7))
                            END DO
                            DEALLOCATE(recvbuf)
                    ENDDO
                ENDDO
            ENDDO
        DEALLOCATE(sendbuf)
    END DO ! End loop on species
  END SUBROUTINE particle_bcs_mpi_blocking

! ________________________________________________________________________________________
!> @brief
!> MPI Boundary condition routine on particles with non blocking mpi communication subroutine
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
  SUBROUTINE particle_bcs_mpi_non_blocking
! ________________________________________________________________________________________

		USE mpi

		IMPLICIT NONE

		INTEGER(isp), PARAMETER :: nvar=8 ! Simple implementation
		INTEGER(isp), DIMENSION(-1:1,-1:1,-1:1) :: nptoexch
		REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuff, recvbuff
		LOGICAL(lp)  :: remove_from_sim
		INTEGER(isp) :: ibuff, nbuff
		INTEGER(isp) :: xbd, ybd, zbd
		INTEGER(isp) :: mpitag, count
		INTEGER(idp), ALLOCATABLE, DIMENSION(:,:,:,:) :: npart_recv, npart_send
		INTEGER(isp) :: dest, src, ireq
		INTEGER(isp), DIMENSION(:), ALLOCATABLE :: requests
		INTEGER(idp) :: ispecies, i, ix, iy, iz, npcurr, ipart
		INTEGER(idp) :: ixtile, iytile, iztile, ispec, nmax
		REAL(num) :: part_xyz
		TYPE(particle_species), POINTER :: currsp
		TYPE(particle_tile), POINTER :: curr

		mpitag=0_isp
		ALLOCATE(npart_send(1:nspecies,-1:1,-1:1,-1:1))
		ALLOCATE(npart_recv(1:nspecies,-1:1,-1:1,-1:1))
		ALLOCATE(requests(27*2))

		npart_recv=0
		npart_send=0

	! POST IRECV PARTICLES FROM ADJACENT DOMAINS
	! ----- POST ISEND FOR THE NUMBER OF PARTICLES
	ireq=1
	DO iz = -1, 1
		DO iy = -1, 1
			DO ix = -1, 1
				IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
				count=nspecies
				dest = neighbour(ix,iy,iz)
				CALL MPI_IRECV(npart_recv(1:count,ix,iy,iz), count,  MPI_INTEGER8, dest, mpitag,    &
								comm, requests(ireq), errcode)
				ireq=ireq+1
			END DO
		END DO
	END DO


		! GET NUMBER OF PARTICLES AT BORDER OF CURRENT DOMAIN (INIT SEND BUFFER)
		nbuff=0
		DO ispecies=1, nspecies
			currsp => species_parray(ispecies)
		DO iztile=1, ntilez !LOOP ON TILES
			DO iytile=1, ntiley
				DO ixtile=1, ntilex
					curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
					IF (.NOT. curr%subdomain_bound) CYCLE
					nbuff=nbuff+ curr%np_tile(1)
				END DO
			END DO
		END DO
	END DO

	ALLOCATE(sendbuff(1:nbuff*nvar,-1:1,-1:1,-1:1))
	! PUT PARTICLES TO BE SENT IN BUFFER
	nptoexch=0
	DO ispecies=1, nspecies !LOOP ON SPECIES
			! Init send recv buffers
			currsp => species_parray(ispecies)
			DO iztile=1, ntilez !LOOP ON TILES
					DO iytile=1, ntiley
							DO ixtile=1, ntilex
									curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
									! If not subdomain border, nothing to do
									IF (.NOT. curr%subdomain_bound) CYCLE
										! Else, search for outbound particles
										part_xyz=0.
										! Identify outbounds particles
										npcurr=curr%np_tile(1)
										DO i = npcurr,1,-1 !LOOP ON PARTICLES
												xbd = 0
												ybd = 0
												zbd = 0
												remove_from_sim = .FALSE.
												part_xyz = curr%part_x(i)
												! Particle has left this processor
												IF (part_xyz .LT. x_min_local) THEN
														xbd = -1
														IF (x_min_boundary) THEN
															SELECT CASE (pbound_x_min)
															CASE (1_idp) ! absorbing
																remove_from_sim=.TRUE.
															CASE (2_idp) ! Reflecting
																	curr%part_x(i) = part_xyz + dx
																	curr%part_ux(i) = - curr%part_ux(i)
																	xbd=0
															CASE DEFAULT ! periodic
																	curr%part_x(i) = part_xyz + length_x
																END SELECT
														ENDIF
												ENDIF
												! Particle has left this processor
												IF (part_xyz .GE. x_max_local) THEN
														xbd = 1
														IF (x_max_boundary) THEN
															SELECT CASE (pbound_x_max)
															CASE (1_idp) ! absorbing
																remove_from_sim=.TRUE.
															CASE (2_idp) ! Reflecting
																	curr%part_x(i) = part_xyz - dx
																	curr%part_ux(i) = - curr%part_ux(i)
																	xbd=0
															CASE DEFAULT ! periodic
																	curr%part_x(i) = part_xyz - length_x
																END SELECT
														ENDIF
												ENDIF

												part_xyz = curr%part_y(i)
												! Particle has left this processor
												IF ((part_xyz .LT. y_min_local) .AND. (c_dim .EQ. 3)) THEN
														ybd = -1
														IF (y_min_boundary) THEN
															SELECT CASE (pbound_y_min)! absorbing
															CASE (1_idp)
																remove_from_sim=.TRUE.
															CASE (2_idp) ! Reflecting
																	curr%part_y(i) = part_xyz + dy
																	curr%part_uy(i) = - curr%part_uy(i)
																	ybd=0
															CASE DEFAULT ! periodic
																	curr%part_y(i) = part_xyz + length_y
																END SELECT
														ENDIF
												ENDIF

												! Particle has left this processor
												IF ((part_xyz .GE. y_max_local) .AND. (c_dim .EQ. 3)) THEN
														ybd = 1
														IF (y_max_boundary) THEN
															SELECT CASE (pbound_y_max)
															CASE (1_idp) ! absorbing
																remove_from_sim=.TRUE.
															CASE (2_idp) ! Reflecting
																	curr%part_y(i) = part_xyz - dy
																	curr%part_uy(i) = - curr%part_uy(i)
																	ybd=0
															CASE DEFAULT ! periodic
																	curr%part_y(i) = part_xyz - length_y
																END SELECT
														ENDIF
												ENDIF

												part_xyz = curr%part_z(i)
												! Particle has left this processor
												IF (part_xyz .LT. z_min_local+zgrid) THEN
														zbd = -1
														IF (z_min_boundary) THEN
																SELECT CASE (pbound_z_min)
																CASE (1_idp) ! absorbing
																	remove_from_sim=.TRUE.
																CASE (2_idp) ! Reflecting
																	curr%part_z(i) = part_xyz + dz
																	curr%part_uz(i) = - curr%part_uz(i)
																	zbd=0
																CASE DEFAULT ! periodic
																	curr%part_z(i) = part_xyz + length_z
																END SELECT
														ENDIF
												ENDIF

												! Particle has left this processor
												IF (part_xyz .GE. z_max_local+zgrid) THEN
														zbd = 1
														! Particle has left the system
														IF (z_max_boundary) THEN
															SELECT CASE (pbound_z_max)
															CASE (1_idp) ! absorbing
																remove_from_sim=.TRUE.
															CASE (2_idp) ! Reflecting
																	curr%part_z(i) = part_xyz - dz
																	curr%part_uz(i) = - curr%part_uz(i)
																	zbd=0
															CASE DEFAULT ! periodic
																	curr%part_z(i) = part_xyz - length_z
																END SELECT
														ENDIF
												ENDIF

												IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
												! Particle has left processor, send it to its neighbour
													IF (.NOT. remove_from_sim) THEN
																		ibuff=nptoexch(xbd,ybd,zbd)*nvar+1
																		sendbuff(ibuff,xbd,ybd,zbd)    = curr%part_x(i)
																		sendbuff(ibuff+1,xbd,ybd,zbd)  = curr%part_y(i)
																		sendbuff(ibuff+2,xbd,ybd,zbd)  = curr%part_z(i)
																		sendbuff(ibuff+3,xbd,ybd,zbd)  = curr%part_ux(i)
																		sendbuff(ibuff+4,xbd,ybd,zbd)  = curr%part_uy(i)
																		sendbuff(ibuff+5,xbd,ybd,zbd)  = curr%part_uz(i)
																		sendbuff(ibuff+6,xbd,ybd,zbd)  = curr%part_gaminv(i)
																		sendbuff(ibuff+7,xbd,ybd,zbd)  = curr%pid(i,wpid)
																		npart_send(ispecies, xbd,ybd,zbd)=npart_send(ispecies,xbd,ybd,zbd)+1
																		nptoexch(xbd,ybd,zbd) = nptoexch(xbd,ybd,zbd)+1
																		! Remove particle of current species from current tile
													ENDIF
													CALL rm_particles_from_species(currsp, ixtile, iytile, iztile, i)
												ENDIF
										ENDDO !END LOOP ON PARTICLES
									ENDDO
							 ENDDO
						ENDDO ! END LOOP ON TILES
		ENDDO ! END LOOP ON SPECIES

		! ----- POST ISEND FOR THE NUMBER OF PARTICLES
		DO iz = -1, 1
			DO iy = -1, 1
				DO ix = -1, 1
				    IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
					count=nspecies
					src = INT(neighbour(ix,iy,iz),isp)
					CALL MPI_ISEND(npart_send(1:count,ix,iy,iz), count,  MPI_INTEGER8, src, mpitag,    &
									comm, requests(ireq), errcode)
					ireq=ireq+1
				END DO
			END DO
	    END DO

		CALL MPI_WAITALL(ireq-1_isp,requests, MPI_STATUSES_IGNORE, errcode)
		requests=0_isp
		ireq=1

		! ----- POST IRECV FOR PARTICLE DATA
		nmax=nvar*MAXVAL(SUM(npart_recv,1))
		ALLOCATE(recvbuff(nmax,-1:1,-1:1,-1:1))
		DO iz = -1, 1
			DO iy = -1, 1
				DO ix = -1, 1
					count=nvar*SUM(npart_recv(:,ix,iy,iz))
				    IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
					IF (count .GT. 0) THEN
						src = INT(neighbour(ix,iy,iz),isp)
						CALL MPI_IRECV(recvbuff(1:count,ix,iy,iz),count, MPI_DOUBLE_PRECISION,src,MPI_ANY_TAG, &
										comm, requests(ireq),errcode)
						ireq=ireq+1
					ENDIF
				END DO
			END DO
	    END DO

		! ----- POST ISEND FOR PARTICLE DATA
		DO iz = -1, 1
			DO iy = -1, 1
				DO ix = -1, 1
					count=nvar*SUM(npart_send(:,ix,iy,iz))
				    IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
					IF (count .GT. 0) THEN
						dest = INT(neighbour(ix,iy,iz),isp)
						CALL MPI_ISEND(sendbuff(1:count,ix,iy,iz),count, MPI_DOUBLE_PRECISION,dest,mpitag, &
										comm, requests(ireq),errcode)
						ireq=ireq+1
					ENDIF
				END DO
			END DO
	    END DO

	    ! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA
		count=ireq-1
		IF (count .GT. 0_isp) THEN
			CALL MPI_WAITALL(count,requests, MPI_STATUSES_IGNORE, errcode)
		ENDIF

		! ----- ADD PARTICLES FROM RECV BUFF TO SPECIES ARRAY
		DO iz = -1, 1
			DO iy = -1, 1
				DO ix = -1, 1
					IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
					ispec=0
					DO ispecies=1,nspecies
						currsp=> species_parray(ispecies)
						DO ipart=1,nvar*npart_recv(ispecies,ix,iy,iz),nvar
							ibuff=ispec+ipart
							CALL add_particle_to_species(currsp, recvbuff(ibuff,ix,iy,iz), &
							recvbuff(ibuff+1,ix,iy,iz), recvbuff(ibuff+2,ix,iy,iz), 	   &
							recvbuff(ibuff+3,ix,iy,iz), recvbuff(ibuff+4,ix,iy,iz), 	   &
							recvbuff(ibuff+5,ix,iy,iz), recvbuff(ibuff+6,ix,iy,iz),        &
							recvbuff(ibuff+7,ix,iy,iz))
						END DO
						ispec=ispec+nvar*npart_recv(ispecies,ix,iy,iz)
					END DO
				END DO
			END DO
	    END DO


  DEALLOCATE(sendbuff,recvbuff,npart_send,npart_recv,requests)
  END SUBROUTINE particle_bcs_mpi_non_blocking

! ________________________________________________________________________________________
!
!> @brief
!> MPI Boundary condition routine for particles in 2D x,z geometry
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
  SUBROUTINE particle_bcs_mpi_non_blocking_2d
! ________________________________________________________________________________________

		USE mpi

    INTEGER(isp), PARAMETER                  :: nvar=8 ! Simple implementation
    INTEGER(isp), DIMENSION(-1:1,-1:1)       :: nptoexch
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: sendbuff, recvbuff
    INTEGER(isp)                             :: ibuff, nbuff
    INTEGER(isp) :: xbd, zbd
    INTEGER(isp) :: mpitag, count
    INTEGER(idp), ALLOCATABLE, DIMENSION(:,:,:) :: npart_recv, npart_send
    INTEGER(isp) :: dest, src, ireq
    INTEGER(isp), DIMENSION(:), ALLOCATABLE :: requests
    LOGICAL(lp)  :: out_of_bounds
    INTEGER(idp) :: ispecies, i, ix, iz, npcurr, ipart
    INTEGER(idp) :: ixtile, iztile, ispec, nmax
    REAL(num) :: part_xyz
    TYPE(particle_species), POINTER :: currsp
    TYPE(particle_tile), POINTER :: curr

    mpitag=0_isp
    ALLOCATE(npart_send(1:nspecies,-1:1,-1:1))
    ALLOCATE(npart_recv(1:nspecies,-1:1,-1:1))
    ALLOCATE(requests(27*2))

    npart_recv=0
    npart_send=0

	! POST IRECV PARTICLES FROM ADJACENT DOMAINS
	! ----- POST ISEND FOR THE NUMBER OF PARTICLES
	ireq=1
	DO iz = -1, 1
			DO ix = -1, 1
				IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
				count=nspecies
				dest = neighbour(ix,0,iz)
				CALL MPI_IRECV(npart_recv(1:count,ix,iz), count,  MPI_INTEGER8, dest, MPI_ANY_TAG,    &
								comm, requests(ireq), errcode)
				ireq=ireq+1
			END DO
	END DO



    ! GET NUMBER OF PARTICLES AT BORDER OF CURRENT DOMAIN (INIT SEND BUFFER)
    nbuff=0
    DO ispecies=1, nspecies
    	currsp => species_parray(ispecies)
		DO iztile=1, ntilez !LOOP ON TILES
				DO ixtile=1, ntilex
					curr=>currsp%array_of_tiles(ixtile,1,iztile)
					IF (.NOT. curr%subdomain_bound) CYCLE
					nbuff=nbuff+ curr%np_tile(1)
				END DO
		END DO
	END DO


	ALLOCATE(sendbuff(1:nbuff,-1:1,-1:1))
	! PUT PARTICLES TO BE SENT IN BUFFER
	nptoexch=0
    DO ispecies=1, nspecies !LOOP ON SPECIES
        ! Init send recv buffers
        currsp => species_parray(ispecies)
        DO iztile=1, ntilez !LOOP ON TILES
                DO ixtile=1, ntilex
                    curr=>currsp%array_of_tiles(ixtile,1,iztile)
                    ! If not subdomain border, nothing to do
                    IF (.NOT. curr%subdomain_bound) CYCLE
                    ! Else, search for outbound particles
                    part_xyz=0.
                    ! Identify outbounds particles
                	npcurr=curr%np_tile(1)
                    DO i = npcurr,1,-1 !LOOP ON PARTICLES
                        xbd = 0
                        zbd = 0
                        out_of_bounds = .FALSE.
                        part_xyz = curr%part_x(i)
                        ! Particle has left this processor
                        IF (part_xyz .LT. x_min_local) THEN
                            xbd = -1
                            IF (x_min_boundary) THEN
                            	SELECT CASE (pbound_x_min)
                            	CASE (1_idp) ! absorbing
                            		CALL rm_particles_from_species_2d(currsp, &
                                ixtile, iztile, i)
                            		CYCLE
                            	CASE DEFAULT ! periodic
                                	curr%part_x(i) = part_xyz + length_x
                                END SELECT
                            ENDIF
                        ENDIF
                        ! Particle has left this processor
                        IF (part_xyz .GE. x_max_local) THEN
                            xbd = 1
                            IF (x_max_boundary) THEN
                            	SELECT CASE (pbound_x_max)
                            	CASE (1_idp) ! absorbing
                            		CALL rm_particles_from_species_2d(currsp, &
                                ixtile, iztile, i)
                            		CYCLE
                            	CASE DEFAULT ! periodic
                                	curr%part_x(i) = part_xyz - length_x
                                END SELECT
                            ENDIF
                        ENDIF

                        part_xyz = curr%part_z(i)
                        ! Particle has left this processor
                        IF (part_xyz .LT. z_min_local) THEN
                            zbd = -1
                            IF (z_min_boundary) THEN
              								SELECT CASE (pbound_z_min)
              								CASE (1_idp) ! absorbing
              									CALL rm_particles_from_species_2d(currsp, &
                                 ixtile, iztile, i)
              									CYCLE
              								CASE DEFAULT ! periodic
              									curr%part_z(i) = part_xyz + length_z
              								END SELECT
                            ENDIF
                        ENDIF

                        ! Particle has left this processor
                        IF (part_xyz .GE. z_max_local) THEN
                            zbd = 1
                            ! Particle has left the system
                            IF (z_max_boundary) THEN
                            	SELECT CASE (pbound_z_max)
                            	CASE (1_idp) ! absorbing
                            		CALL rm_particles_from_species_2d(currsp, &
                                ixtile, iztile, i)
                            		CYCLE
                            	CASE DEFAULT ! periodic
                                	curr%part_z(i) = part_xyz - length_z
                                END SELECT
                            ENDIF
                        ENDIF

                        IF (ABS(xbd) + ABS(zbd) .GT. 0) THEN
                        ! Particle has left processor, send it to its neighbour
                            ibuff=nptoexch(xbd,zbd)*nvar+1
                            sendbuff(ibuff,xbd,zbd)    = curr%part_x(i)
                            sendbuff(ibuff+1,xbd,zbd)  = curr%part_z(i)
                            sendbuff(ibuff+2,xbd,zbd)  = curr%part_ux(i)
                            sendbuff(ibuff+3,xbd,zbd)  = curr%part_uy(i)
                            sendbuff(ibuff+4,xbd,zbd)  = curr%part_uz(i)
                            sendbuff(ibuff+5,xbd,zbd)  = curr%part_gaminv(i)
                            sendbuff(ibuff+6,xbd,zbd)  = curr%pid(i,wpid)
                        	  npart_send(ispecies, xbd,zbd)=npart_send(ispecies,xbd,zbd)+1
                            nptoexch(xbd,zbd) = nptoexch(xbd,zbd)+1
                            ! Remove particle of current species from current tile
                       	    CALL rm_particles_from_species_2d(currsp, ixtile, iztile, i)
                        ENDIF
                    ENDDO !END LOOP ON PARTICLES
                  ENDDO
            ENDDO ! END LOOP ON TILES
		ENDDO ! END LOOP ON SPECIES


		! ----- POST ISEND FOR THE NUMBER OF PARTICLES
		DO iz = -1, 1
				DO ix = -1, 1
				    IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
					count=nspecies
					src = INT(neighbour(ix,0,iz),isp)
					CALL MPI_ISEND(npart_send(1:count,ix,iz), count,  MPI_INTEGER8, src, mpitag,    &
									comm, requests(ireq), errcode)
					ireq=ireq+1
				END DO
	    END DO


		CALL MPI_WAITALL(ireq-1_isp,requests, MPI_STATUSES_IGNORE, errcode)
		requests=0_isp
		ireq=1

		! ----- POST IRECV FOR PARTICLE DATA
		nmax=nvar*MAXVAL(SUM(npart_recv,1))
		ALLOCATE(recvbuff(nmax,-1:1,-1:1))
		DO iz = -1, 1
				DO ix = -1, 1
					count=nvar*SUM(npart_recv(:,ix,iz))
				    IF (ABS(ix)  + ABS(iz) .EQ. 0) CYCLE
					IF (count .GT. 0) THEN
						src = INT(neighbour(ix,0,iz),isp)
						CALL MPI_IRECV(recvbuff(1:count,ix,iz),count, MPI_DOUBLE_PRECISION,src,MPI_ANY_TAG, &
										comm, requests(ireq),errcode)
						ireq=ireq+1
					ENDIF
				END DO
	    END DO


		! ----- POST ISEND FOR PARTICLE DATA
		DO iz = -1, 1
				DO ix = -1, 1
					count=nvar*SUM(npart_send(:,ix,iz))
				    IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
					IF (count .GT. 0) THEN
						dest = INT(neighbour(ix,0,iz),isp)
						CALL MPI_ISEND(sendbuff(1:count,ix,iz),count, MPI_DOUBLE_PRECISION,dest,mpitag, &
										comm, requests(ireq),errcode)
						ireq=ireq+1
					ENDIF
				END DO
	    END DO


	    ! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA
		count=ireq-1
		IF (count .GT. 0_isp) THEN
			CALL MPI_WAITALL(count,requests, MPI_STATUSES_IGNORE, errcode)
		ENDIF


		! ----- ADD PARTICLES FROM RECV BUFF TO SPECIES ARRAY
		DO iz = -1, 1
				DO ix = -1, 1
					IF (ABS(ix) + ABS(iz) .EQ. 0) CYCLE
					ispec=0
					DO ispecies=1,nspecies
						currsp=> species_parray(ispecies)
						DO ipart=1,nvar*npart_recv(ispecies,ix,iz),nvar
							ibuff=ispec+ipart
							CALL add_particle_to_species_2d(currsp, recvbuff(ibuff,ix,iz), &
							recvbuff(ibuff+1,ix,iz), recvbuff(ibuff+2,ix,iz), 	   &
							recvbuff(ibuff+3,ix,iz), recvbuff(ibuff+4,ix,iz), 	   &
							recvbuff(ibuff+5,ix,iz), recvbuff(ibuff+6,ix,iz))
						END DO
						ispec=ispec+nvar*npart_recv(ispecies,ix,iz)
					END DO
				END DO
	    END DO


  DEALLOCATE(sendbuff,recvbuff,npart_send,npart_recv,requests)
  END SUBROUTINE particle_bcs_mpi_non_blocking_2d


! ______________________________________________________________________________________
!> @brief
!> This subroutine combined in a single routine the particle communications between tiles
!> and between MPI domains for 3D
!
!> @author
!> Mathieu Lobet
!
!> @date
!> May 2016
!
  SUBROUTINE particle_bcs_tiles_and_mpi_3d
! ______________________________________________________________________________________
	USE omp_lib
	USE communications
	USE precomputed
	USE params
	USE mpi
	IMPLICIT NONE

	INTEGER(idp)                    :: is, ix, iy, iz
	INTEGER(idp)                    :: i,i2,i3
	INTEGER(idp)                    :: indx, indy, indz, ipx, ipy, ipz
	INTEGER(idp)                    :: xbd,ybd,zbd
	INTEGER(idp)                    :: k,j,ib,ibs
	INTEGER(isp)                    :: ireq
	INTEGER(isp)                    :: dest, src
	INTEGER(idp)                    :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
	TYPE(particle_species), POINTER :: curr
	TYPE(particle_tile), POINTER    :: curr_tile
	REAL(num)                       :: partx, party, partz, partux, partuy, partuz, partw, gaminv
	INTEGER(idp)                                      :: nthreads_tot
	INTEGER(idp)                                      :: nthreads_loop1, nthreads_loop2
	INTEGER(idp), dimension(:,:), ALLOCATABLE         :: mpi_npart
	REAL(num), dimension(:,:,:,:), ALLOCATABLE        :: bufsend
	REAL(num), dimension(:,:), ALLOCATABLE            :: recvbuf
	TYPE(mpi_tile_buffer), dimension(:,:,:,:), ALLOCATABLE :: tilebuf
	INTEGER(isp), DIMENSION(:,:), ALLOCATABLE         :: nrecv_buf
	INTEGER(isp), DIMENSION(:), ALLOCATABLE           :: reqs
	INTEGER(isp)                                      :: nrecv_buf_tot,npos, typebuffer
	INTEGER(idp)                                      :: new_mpi_buf_size,old_mpi_buf_size
	REAL(num)                                         :: nx0_grid_tile_dx
	REAL(num)                                         :: ny0_grid_tile_dy,nz0_grid_tile_dz
	INTEGER(isp)                                      :: stats(2)
	INTEGER(idp)                                      :: recvbuf_index(27)
	INTEGER(idp)                                      :: lvect


! ________________________________________
! Checking
#if defined(DEBUG)
		DO iz=1, ntilez
			DO iy=1, ntiley
				DO ix=1, ntilex
		
					curr_tile=>curr%array_of_tiles(ix,iy,iz)
					nptile=curr_tile%np_tile(1)
					
					DO i=1,nptile
						partx=curr_tile%part_x(i)
						party=curr_tile%part_y(i)
						partz=curr_tile%part_z(i)
						partux=curr_tile%part_ux(i)
						partuy=curr_tile%part_uy(i)
						partuz=curr_tile%part_uz(i)
						gaminv=curr_tile%part_gaminv(i)
						partw=curr_tile%pid(i,wpid)

						IF ((partx .LT. curr_tile%x_tile_min) .OR.   &
						    (partx .GE. curr_tile%x_tile_max) .OR.   &
						    (party .LT. curr_tile%y_tile_min) .OR.   &
						    (party .GE. curr_tile%y_tile_max) .OR.  &
						    (partz .LT. curr_tile%z_tile_min+zgrid) .OR. &
						    (partz .GE. curr_tile%z_tile_max+zgrid)) THEN
						    
						    
						    WRITE(0,'("ERROR: particle outside the domain")')
						    WRITE(0,'("x:",E12.5,X,E12.5,X,E12.5)') curr_tile%x_tile_min,partx,curr_tile%x_tile_max
						    WRITE(0,'("y:",E12.5,X,E12.5,X,E12.5)') curr_tile%y_tile_min,party,curr_tile%y_tile_max
						    WRITE(0,'("z:",E12.5,X,E12.5,X,E12.5)') curr_tile%z_tile_min,partz,curr_tile%z_tile_max
						    WRITE(0,*)
						    
						ENDIF
						
					ENDDO
		
				ENDDO
			ENDDO
		ENDDO
		WRITE(0,*) " Checking after particle boundary condition ok"
#endif

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
		nthreads_loop1=MIN(nspecies,nthreads_tot)
		nthreads_loop2=MAX(1_idp,nthreads_tot/nthreads_loop1)
	ELSE
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF

	lvect = 64

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
	write(0,*) "Part 1 - Determine the particle to be exchanged with other tiles or with other MPI domains"
#endif

	ALLOCATE(mpi_npart(27,nspecies))
	ALLOCATE(tilebuf(ntilex,ntiley,ntilez,nspecies))

! 	DO is=1, nspecies
!   	curr=> species_parray(is)
! 		DO iz=1, ntilez
! 			DO iy=1, ntiley
! 				DO ix=1, ntilex
!
! 					curr_tile=>curr%array_of_tiles(ix,iy,iz)
! 					nptile=curr_tile%np_tile(1)
!
! 					! Allocation of the buffer
! 					IF (curr_tile%subdomain_bound) THEN
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_x(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_y(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_z(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_ux(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_uy(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_uz(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%part_gaminv(mpi_buf_size,27))
! 						ALLOCATE(tilebuf(ix,iy,iz,is)%pid(mpi_buf_size,27))
! 					ENDIF
! 					tilebuf(ix,iy,iz,is)%npart(1:27) = 0
!
! 				ENDDO
! 			ENDDO
! 		ENDDO
! 	ENDDO


	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP PRIVATE(curr,is,ib,k,nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz,&
	!$OMP nx0_grid_tile_dx,ny0_grid_tile_dy,nz0_grid_tile_dz,xbd,ybd,zbd,gaminv,&
	!$OMP partw,indx,indy,indz,partx,party,partz,curr_tile,nptile,partux,partuy,partuz,&
	!$OMP j) &
	!$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
	!$OMP length_x,length_y,length_z,dxs2,dys2,dzs2, &
	!$OMP x_min_boundary,x_max_boundary,y_min_boundary,y_max_boundary,z_min_boundary,z_max_boundary,  &
	!$OMP pbound_x_min,pbound_x_max,pbound_y_min,pbound_y_max,pbound_z_min,pbound_z_max, &
	!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz,mpi_npart,tilebuf,mpi_buf_size,lvect,zgrid) &
	!$OMP NUM_THREADS(nthreads_loop1)
	! LOOP ON SPECIES
	DO is=1, nspecies

		curr=> species_parray(is)
		! Get first tiles dimensions (may be different from last tile)
		nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
		ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
		nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile

		nx0_grid_tile_dx = 1._num/(nx0_grid_tile*dx)
		ny0_grid_tile_dy = 1._num/(ny0_grid_tile*dy)
		nz0_grid_tile_dz = 1._num/(nz0_grid_tile*dz)

		mpi_npart(:,is) = 0
		! LOOP ON TILES
		DO ipz=1,3
			DO ipy=1,3
				DO ipx=1,3
				
					!$OMP PARALLEL DO DEFAULT(NONE) &
					!$OMP SHARED(curr,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, &
					!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile, &
					!$OMP pbound_x_min,pbound_x_max,pbound_y_min,pbound_y_max,pbound_z_min,pbound_z_max, &
					!$OMP length_x,length_y,length_z,tilebuf,mpi_npart, &
					!$OMP x_min_boundary,x_max_boundary,y_min_boundary,y_max_boundary,z_min_boundary,z_max_boundary, &
					!$OMP nx0_grid_tile_dx,ny0_grid_tile_dy,nz0_grid_tile_dz,dxs2,dys2,dzs2,mpi_buf_size,lvect,zgrid)  &
					!$OMP FIRSTPRIVATE(ipx,ipy,ipz,is) &
					!$OMP PRIVATE(ix,iy,iz,i,ib,k,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
					!$OMP indx,indy,indz,xbd,ybd,zbd,i2,i3,new_mpi_buf_size,old_mpi_buf_size)  &
					!$OMP COLLAPSE(3) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)
					! LOOP ON TILES
					DO iz=ipz, ntilez,3
						DO iy=ipy, ntiley,3
							DO ix=ipx, ntilex,3
								curr_tile=>curr%array_of_tiles(ix,iy,iz)
								nptile=curr_tile%np_tile(1)

								! Allocation of the buffer
								IF (curr_tile%subdomain_bound) THEN
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_x(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_y(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_z(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_ux(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_uy(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_uz(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_gaminv(mpi_buf_size,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%pid(mpi_buf_size,27))
								ENDIF
								tilebuf(ix,iy,iz,is)%npart(1:27) = 0

#if defined(DEBUG) && (DEBUG==3)
	write(0,'("Loop on particles inside tiles: ",I2,X,I2,X,I2)')ix,iy,iz
#endif

								! Loop on particles inside tiles
								DO i=nptile, 1, -1

! 								DO i2=nptile, 1, -lvect
! 									DO i3 = MIN(lvect,i2),1,-1
! 								  i = i2 - i3 + 1

									partx=curr_tile%part_x(i)
									party=curr_tile%part_y(i)
									partz=curr_tile%part_z(i)
									partux=curr_tile%part_ux(i)
									partuy=curr_tile%part_uy(i)
									partuz=curr_tile%part_uz(i)
									gaminv=curr_tile%part_gaminv(i)
									partw=curr_tile%pid(i,wpid)

									! Case 1: if particle did not leave tile nothing to do
									IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    &
									.AND. ((party .GE. curr_tile%y_tile_min) .AND. (party .LT. curr_tile%y_tile_max))  &
									.AND. ((partz .GE. curr_tile%z_tile_min+zgrid) .AND. (partz .LT. curr_tile%z_tile_max+zgrid))) &
									CYCLE

									! Case 2: if particle left MPI domain
									IF (((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) .OR. &
									 ((party .LT. y_min_local) .OR. (party .GE. y_max_local)) .OR. &
									 ((partz .LT. z_min_local+zgrid) .OR. (partz .GE. z_max_local+zgrid))) THEN

									! Then we determine in which domain this particle is going
									xbd = 0
									ybd = 0
									zbd = 0

									! Particle has left this processor -x
									IF (partx .LT. x_min_local) THEN
										xbd = -1
										IF (x_min_boundary) THEN
											SELECT CASE (pbound_x_min)
												CASE (1_idp) ! absorbing
												CALL rm_particle_at_tile(curr,ix,iy,iz,i)
												CYCLE
												CASE DEFAULT ! periodic
													curr_tile%part_x(i) = partx + length_x
												END SELECT
											ENDIF

										! Particle has left this processor +x
										ELSE IF (partx .GE. x_max_local) THEN
											xbd = 1
											IF (x_max_boundary) THEN
												SELECT CASE (pbound_x_max)
													CASE (1_idp) ! absorbing
													CALL rm_particle_at_tile(curr,ix,iy,iz,i)
														CYCLE
													 CASE DEFAULT ! periodic
													curr_tile%part_x(i) = partx - length_x
												END SELECT
											ENDIF
										ENDIF

										! Particle has left this processor -y
										IF ((party .LT. y_min_local)) THEN
											ybd = -1
											 IF (y_min_boundary) THEN
															SELECT CASE (pbound_y_min)! absorbing
															CASE (1_idp)
													CALL rm_particle_at_tile(curr,ix,iy,iz,i)
													CYCLE
												CASE DEFAULT ! periodic
													curr_tile%part_y(i) = party + length_y
												END SELECT
											ENDIF
											
										! Particle has left this processor +y
										ELSE IF ((party .GE. y_max_local)) THEN
											ybd = 1
											IF (y_max_boundary) THEN
												SELECT CASE (pbound_y_max)
												CASE (1_idp) ! absorbing
																CALL rm_particle_at_tile(curr,ix,iy,iz,i)
													CYCLE
												CASE DEFAULT ! periodic
													curr_tile%part_y(i) = party - length_y
												END SELECT
											ENDIF
										ENDIF

										! Particle has left this processor -z
										IF (partz .LT. z_min_local+zgrid) THEN
											zbd = -1
											IF (z_min_boundary) THEN
												SELECT CASE (pbound_z_min)
													CASE (1_idp) ! absorbing
														CALL rm_particle_at_tile(curr,ix,iy,iz,i)
													CYCLE
													CASE DEFAULT ! periodic
														curr_tile%part_z(i) = partz + length_z
													END SELECT
											ENDIF
										
										! Particle has left this processor +z
										ELSE IF (partz .GE. z_max_local+zgrid) THEN
											zbd = 1
											! Particle has left the system
											IF (z_max_boundary) THEN
												SELECT CASE (pbound_z_max)
												CASE (1_idp) ! absorbing
													CALL rm_particle_at_tile(curr,ix,iy,iz,i)
													CYCLE
												CASE DEFAULT ! periodic
														curr_tile%part_z(i) = partz - length_z
												END SELECT
											ENDIF
										ENDIF

                  ! Particle has left processor, we put it in a local buffer
                  IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN

                    ! Boundary index 1-27
                    ib = 2+xbd + (1+ybd)*3 + (1+zbd)*9
                    ! Current particle available index
                    k = tilebuf(ix,iy,iz,is)%npart(ib) + 1

                    ! Particle properties are placed in a buffer for each com direction
                    tilebuf(ix,iy,iz,is)%part_x(k,ib) = curr_tile%part_x(i)
                    tilebuf(ix,iy,iz,is)%part_y(k,ib) = curr_tile%part_y(i)
                    tilebuf(ix,iy,iz,is)%part_z(k,ib) = curr_tile%part_z(i)
                    tilebuf(ix,iy,iz,is)%part_ux(k,ib) = curr_tile%part_ux(i)
                    tilebuf(ix,iy,iz,is)%part_uy(k,ib) = curr_tile%part_uy(i)
                    tilebuf(ix,iy,iz,is)%part_uz(k,ib) = curr_tile%part_uz(i)
                    tilebuf(ix,iy,iz,is)%part_gaminv(k,ib) = curr_tile%part_gaminv(i)
                    tilebuf(ix,iy,iz,is)%pid(k,ib) = curr_tile%pid(i,wpid)

                    IF (k.eq.SIZE(tilebuf(ix,iy,iz,is)%part_x,1)) THEN
                      old_mpi_buf_size = SIZE(tilebuf(ix,iy,iz,is)%part_x,1)
                      new_mpi_buf_size = SIZE(tilebuf(ix,iy,iz,is)%part_x,1)*2
                      WRITE(0,'(" WARNING: Tile buffer array has been resized: nbpart = ",I7,&
" new buffer size = ",I7)') k,new_mpi_buf_size
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_x, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_y, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_z, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_ux, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_uy, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_uz, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%part_gaminv, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                      CALL resize_2D_array_real(tilebuf(ix,iy,iz,is)%pid, old_mpi_buf_size,new_mpi_buf_size,27_idp,27_idp)
                    ENDIF
                    
                    ! Update of the number of particles
                    tilebuf(ix,iy,iz,is)%npart(ib) = k

                    ! The particle is deleted
                    CALL rm_particle_at_tile(curr,ix,iy,iz,i)

									ENDIF


									  CYCLE

                  ENDIF !end if particle left MPI domain


									! Case 3: particles changed tile. Tranfer particle to new tile
									! Get new indexes of particle in array of tiles
									indx = MIN(FLOOR((partx-x_min_local+dxs2)*(nx0_grid_tile_dx),idp)+1,ntilex)
									indy = MIN(FLOOR((party-y_min_local+dys2)*(ny0_grid_tile_dy),idp)+1,ntiley)
									indz = MIN(FLOOR((partz-(z_min_local+zgrid)+dzs2)*(nz0_grid_tile_dz),idp)+1,ntilez)
									!if ((indx.le.0).or.(indy.le.0).or.(indz.le.0)) THEN
                   !print*,'xmin',x_min_local,'xmax',x_max_local,'x',partx,xbd
                   !print*,'ymin',y_min_local,'ymax',y_max_local,'y',party,ybd
                   !print*,'zmin',z_min_local,'zmax',z_max_local,'z',partz,zbd
                  !endif
									CALL rm_particle_at_tile(curr,ix,iy,iz,i)

									!IF ((indx < 1).OR.(indy < 1).OR.(indz < 1)) THEN
									!  Write(0,*) indx,indy,indz
									!  Write(0,*) partx,party,partz
									!  Write(0,*) i,nptile
									!ENDIF

									CALL add_particle_at_tile(curr, indx,indy,indz, &
										 partx, party, partz, partux, partuy, partuz, gaminv, partw)


                  !ENDDO
								!End loop on particles
								END DO

								! Reduction of the total number of particle to be
								! exchanged in every direction for MPI
								! Here we use a critical region since we expect the thread
								! to not be perfectly synchronized
								! In this part appears to have negative consequences on the parallelization,
								! this process can be performed in sequencial outside the species loop
								IF (curr_tile%subdomain_bound) THEN
								!OMP CRITICAL
                  mpi_npart(1:27,is) = mpi_npart(1:27,is) + tilebuf(ix,iy,iz,is)%npart(1:27)
								!OMP END CRITICAL

								ENDIF

							END DO
						END DO
					END DO ! END LOOP ON TILES
					!$OMP END PARALLEL DO

        		END DO
        	END DO
        END DO
    END DO ! END LOOP ON SPECIES
    !$OMP END PARALLEL DO


#if defined(DEBUG) && (DEBUG==3)
	write(0,*) "Part 2 - Creation of the send buffer for the MPI communications"
#endif

    ! ____________________________________________________________________________________
    ! Part 2 - Creation of the send buffer for the MPI communications
    !
    ! Here, a global buffer is created for the MPI communications. We copy in parallel
    ! the local buffers of the tiles into the global buffer using nested parallel regions,
    ! one for the species and the other on the communication directions
    !
    ! There are 27 directions in 3D but only the 8 faces of the cubic domain represent
    ! a significant number of particles to exchange when the plasma kinetic
    ! distribution is uniform
    !
    ! For drifted plasmas, this method is not the most efficient

    ALLOCATE(bufsend(MAXVAL(mpi_npart(:,:)),8,27,nspecies))
    !print*,'max part',MAXVAL(mpi_npart(:,:)),MINVAL(mpi_npart(:,:))
    !print*,mpi_npart(:,1)
    !print*
    !print*,mpi_npart(:,2)
    !print*


    !print*,'Creation of the send buffer'
	  !$OMP PARALLEL DO DEFAULT(NONE) &
	  !$OMP PRIVATE(curr,is, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz,&
	  !$OMP ib,k,j) &
	  !$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez, &
    !$OMP dx,dy,dz,mpi_npart,bufsend,tilebuf) &
    !$OMP NUM_THREADS(nthreads_loop1)
    DO is=1, nspecies ! LOOP ON SPECIES

     mpi_npart(:,is) = 0

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(is,tilebuf,mpi_npart,bufsend,ntilez,ntiley,ntilex) &
     !$OMP PRIVATE(ix,iy,iz,k,xbd,ybd,zbd,j,ib) &
     !$OMP COLLAPSE(3) SCHEDULE(runtime) &
     !$OMP NUM_THREADS(nthreads_loop2)
      DO xbd = -1, 1
        DO ybd = -1, 1
          DO zbd = -1, 1

            IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE

            ib = 2+xbd + (1+ybd)*3 + (1+zbd)*9

					  DO iz=1, ntilez! LOOP ON TILES
						DO iy=1, ntiley
						DO ix=1, ntilex

						  k = tilebuf(ix,iy,iz,is)%npart(ib)
					    !print*,'iz,iy,ix, ib:',iz,iy,ix,ib,'k',k

						  IF (k.eq.0) CYCLE

						  j = mpi_npart(ib,is)
						  bufsend(j+1:j+k,1,ib,is) = tilebuf(ix,iy,iz,is)%part_x(1:k,ib)
						  bufsend(j+1:j+k,2,ib,is) = tilebuf(ix,iy,iz,is)%part_y(1:k,ib)
						  bufsend(j+1:j+k,3,ib,is) = tilebuf(ix,iy,iz,is)%part_z(1:k,ib)
						  bufsend(j+1:j+k,4,ib,is) = tilebuf(ix,iy,iz,is)%part_ux(1:k,ib)
						  bufsend(j+1:j+k,5,ib,is) = tilebuf(ix,iy,iz,is)%part_uy(1:k,ib)
						  bufsend(j+1:j+k,6,ib,is) = tilebuf(ix,iy,iz,is)%part_uz(1:k,ib)
						  bufsend(j+1:j+k,7,ib,is) = tilebuf(ix,iy,iz,is)%part_gaminv(1:k,ib)
						  bufsend(j+1:j+k,8,ib,is) = tilebuf(ix,iy,iz,is)%pid(1:k,ib)
						  mpi_npart(ib,is) = j + k

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
	write(0,*) "Part 3 - MPI Communications"
#endif

    ! _______________________________________
    ! Part 3 - MPI Communications

!     Old sequential version
!     print*,'MPI Communications'
!     DO is=1, nspecies ! LOOP ON SPECIES
!       curr=> species_parray(is)
!       DO iz = -1, 1
!         DO iy = -1, 1
!           DO ix = -1, 1
!             IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
!                ib = 2+ix + (1+iy)*3 + (1+iz)*9
!                ipx = -ix
!                ipy = -iy
!                ipz = -iz
!                ! SEND - RECEIVE PARTICLES IN BUFFERS
!                !- Get number of particles in recvbuff
!                nrecv_buf=0
!                dest = INT(neighbour(ix,iy,iz),isp)
!                src  = INT(neighbour(ipx,ipy,ipz),isp)
!                k = mpi_npart(ib,is)
!                CALL MPI_SENDRECV(k, 1_isp, MPI_INTEGER, dest, tag, nrecv_buf, 1_isp, &
!                     MPI_INTEGER, src, tag, comm, status, errcode)
!                ALLOCATE(recvbuf(1:nrecv_buf,8))
!                CALL MPI_SENDRECV(bufsend(1:k,1:8,ib,is), 8*k, mpidbl, dest, tag, &
!                     recvbuf, nrecv_buf*8, mpidbl, src, tag, comm, status, errcode)
!               ! Add received particles to particle arrays
!               DO i = 1,nrecv_buf
!               CALL add_particle_to_species(curr, &
!               recvbuf(i,1),  &
!               recvbuf(i,2),  &
!               recvbuf(i,3),  &
!               recvbuf(i,4),  &
!               recvbuf(i,5),  &
!               recvbuf(i,6),  &
!               recvbuf(i,7),  &
!               recvbuf(i,8))
!               ENDDO
!
!               DEALLOCATE(recvbuf)
!           ENDDO
!         ENDDO
!       ENDDO
!     ENDDO


    ! First, we determine the number of particles to receive

    ALLOCATE(nrecv_buf(27,nspecies))

    nrecv_buf(:,:)=0

    ! Thread version
    IF (.FALSE.) THEN

    ALLOCATE(reqs(2))

    DO is=1, nspecies ! LOOP ON SPECIES
      !curr=> species_parray(is)

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(is,mpi_npart,comm,neighbour,nrecv_buf) &
     !$OMP FIRSTPRIVATE(tag,status,stats,reqs,MPI_STATUSES_IGNORE) &
     !$OMP PRIVATE(ix,iy,iz,k,ipx,ipy,ipz,dest,src,ib,ibs,errcode) &
     !$OMP COLLAPSE(3) SCHEDULE(runtime) &
     !$OMP NUM_THREADS(nthreads_loop2)

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

               dest = INT(neighbour(ix,iy,iz),isp)
               src  = INT(neighbour(ipx,ipy,ipz),isp)

               ! Number of particle
               k = mpi_npart(ib,is)

               ! Exchange
               CALL MPI_Irecv( nrecv_buf(ib,is), 1_isp, MPI_INTEGER, src, &
               INT(ib,isp), comm, reqs(1), errcode)

               CALL MPI_Isend(k, 1_isp, MPI_INTEGER, dest, INT(ib,isp), &
               comm, reqs(2), errcode)

               CALL MPI_Waitall(2_isp,reqs,MPI_STATUSES_IGNORE,errcode)

               !CALL MPI_SENDRECV(k, 1_isp, MPI_INTEGER, dest, ib, nrecv_buf(ib,is), 1_isp, &
               !     MPI_INTEGER, src, ib, comm, status, errcode)

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ENDDO

    ELSE

    ! _______________________________________________________
    ! Sequential version of the previous block
!
! Waitall is done after the loops
    ALLOCATE(reqs(54*nspecies))
    DO is=1, nspecies ! LOOP ON SPECIES
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

               dest = INT(neighbour(ix,iy,iz),isp)
               src  = INT(neighbour(ipx,ipy,ipz),isp)

               ! Number of particle
               k = mpi_npart(ib,is)

               ! Exchange
                 ireq = ireq + 1
                 CALL MPI_Irecv( nrecv_buf(ib,is), 1_isp, MPI_INTEGER, src, &
                 INT(ib,isp), comm, reqs(ireq), errcode)
                 ireq = ireq + 1
                 CALL MPI_Isend(k, 1_isp, MPI_INTEGER, dest, INT(ib,isp), &
                 comm, reqs(ireq), errcode)

               !CALL MPI_SENDRECV(k, 1_isp, MPI_INTEGER, dest, ib, nrecv_buf(ib,is), 1_isp, &
               !     MPI_INTEGER, src, ib, comm, status, errcode)
          ENDDO
        ENDDO
      ENDDO
      CALL MPI_Waitall(ireq,reqs,MPI_STATUSES_IGNORE,errcode)
    ENDDO

! First irecv and then isend
!    ALLOCATE(reqs(54*nspecies))
!    DO is=1, nspecies ! LOOP ON SPECIES
!      ireq = 0
!       DO iz = -1, 1
!         DO iy = -1, 1
!           DO ix = -1, 1
!             IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
!                ! index of the communication direction
!                ib = 2+ix + (1+iy)*3 + (1+iz)*9
!                ! Number of particle
!                k = mpi_npart(ib,is)
!                ! Exchange
!               src  = INT(neighbour(-ix,-iy,-iz),isp)
!               ireq = ireq + 1
!               CALL MPI_Irecv( nrecv_buf(ib,is), 1_isp, MPI_INTEGER, src, &
!               INT(ib,isp), comm, reqs(ireq), errcode)
!
!           ENDDO
!         ENDDO
!       ENDDO
!       DO iz = -1, 1
!         DO iy = -1, 1
!           DO ix = -1, 1
!             IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
!                ! index of the communication direction
!                ib = 2+ix + (1+iy)*3 + (1+iz)*9
!                ! Number of particle
!                k = mpi_npart(ib,is)
!                ! Exchange
!                dest = INT(neighbour(ix,iy,iz),isp)
!                ireq = ireq + 1
!                CALL MPI_Isend(k, 1_isp, MPI_INTEGER, dest, INT(ib,isp), &
!                  comm, reqs(ireq), errcode)
!           ENDDO
!         ENDDO
!       ENDDO
!      CALL MPI_Waitall(ireq,reqs,MPI_STATUSES_IGNORE,errcode)
!    ENDDO

! Wait all is done inside the loops
!     ALLOCATE(reqs(2))
!     DO is=1, nspecies ! LOOP ON SPECIES
!       !curr=> species_parray(is)
!       DO iz = -1, 1
!         DO iy = -1, 1
!           DO ix = -1, 1
!             IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
!                ! index of the communication direction
!                ib = 2+ix + (1+iy)*3 + (1+iz)*9
!
!                ! indexes of the source/destination mpi task
!                ipx = -ix
!                ipy = -iy
!                ipz = -iz
!                !ibs = 2+ipx + (1+ipy)*3 + (1+ipz)*9
!                dest = INT(neighbour(ix,iy,iz),isp)
!                src  = INT(neighbour(ipx,ipy,ipz),isp)
!
!                ! Number of particle
!                k = mpi_npart(ib,is)
!
!                ! Exchange
!                CALL MPI_Irecv( nrecv_buf(ib,is), 1_isp, MPI_INTEGER, src, &
!                INT(ib,isp), comm, reqs(1), errcode)
!                CALL MPI_Isend(k, 1_isp, MPI_INTEGER, dest, INT(ib,isp), &
!                comm, reqs(2), errcode)
!                CALL MPI_Waitall(2_isp,reqs,MPI_STATUSES_IGNORE,errcode)
!                !CALL MPI_SENDRECV(k, 1_isp, MPI_INTEGER, dest, ib, nrecv_buf(ib,is), 1_isp, &
!                !     MPI_INTEGER, src, ib, comm, status, errcode)
!
!           ENDDO
!         ENDDO
!       ENDDO
!     ENDDO


    ENDIF
    ! ______________________________________________



    ! Then we exchange the buffers
    ! The communication can be performed in sequential or with multiple thread


    ! LOOP ON SPECIES
    DO is=1, nspecies

      nrecv_buf_tot = SUM(nrecv_buf(:,is))
      curr=> species_parray(is)
      ALLOCATE(recvbuf(1:nrecv_buf_tot+1,8))

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
          npos = npos + nrecv_buf(ib,is)
          ENDDO
        ENDDO
      ENDDO

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(is,mpi_npart,comm,neighbour,nrecv_buf,nrecv_buf_tot,recvbuf_index,MPI_STATUSES_IGNORE, &
     !$OMP bufsend,mpidbl,recvbuf) &
     !$OMP FIRSTPRIVATE(tag,status,stats,reqs) &
     !$OMP PRIVATE(ix,iy,iz,k,j,ipx,ipy,ipz,dest,src,ib,ibs,errcode,typebuffer) &
     !$OMP COLLAPSE(3) SCHEDULE(runtime) &
     !$OMP NUM_THREADS(nthreads_loop2)
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1

            IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE

            ib = 2+ix + (1+iy)*3 + (1+iz)*9

            k = mpi_npart(ib,is)
            j = nrecv_buf(ib,is)

            dest = INT(neighbour(ix,iy,iz),isp)
            src  = INT(neighbour(-ix,-iy,-iz),isp)

            CALL MPI_TYPE_VECTOR(8_isp,INT(j,isp),INT(nrecv_buf_tot+1,isp),MPI_DOUBLE_PRECISION,typebuffer,errcode)
            call MPI_TYPE_COMMIT(typebuffer,errcode)

           ! Exchange
           CALL MPI_Irecv(recvbuf(recvbuf_index(ib),1),1_isp, typebuffer,src, &
           INT(ib,isp),comm,reqs(1),errcode)

           CALL MPI_Isend(bufsend(1:k,1:8,ib,is),INT(8_idp*k,isp),mpidbl,dest,INT(ib,isp), &
           comm, reqs(2), errcode)

           CALL MPI_Waitall(2_isp,reqs,MPI_STATUSES_IGNORE,errcode)

            !CALL MPI_SENDRECV(bufsend(1:k,1:8,ib,is), 8_isp*k, mpidbl, dest, INT(ib,isp), &
            !recvbuf(recvbuf_index(ib),1), 1_isp, typebuffer, src, INT(ib,isp), comm, status, errcode)

            call MPI_TYPE_FREE(typebuffer,errcode)

            !npos = npos + j
            !print*,ib,npos,j,iz,iy,ix

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      ! ________________________________________________

      ELSE

      ! ________________________________________________
      ! Sequential version

! Wait all done outside
! I do not know why this version skips particles
!       npos=1
!       ireq = 0
!       DO iz = -1, 1
!         DO iy = -1, 1
!           DO ix = -1, 1
!
!             IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
!
!             ib = 2+ix + (1+iy)*3 + (1+iz)*9
!
!             k = mpi_npart(ib,is)
!             j = nrecv_buf(ib,is)
!
!             dest = INT(neighbour(ix,iy,iz),isp)
!             src  = INT(neighbour(-ix,-iy,-iz),isp)
!
!
!             !IF (j.gt.0) THEN
!               CALL MPI_TYPE_VECTOR(8_isp,INT(j,isp),INT(nrecv_buf_tot+1,isp),mpidbl,typebuffer,errcode)
!               call MPI_TYPE_COMMIT(typebuffer,errcode)
!               ! Exchange
!               ireq = ireq + 1
!               CALL MPI_Irecv(recvbuf(npos,1),1_isp, typebuffer,src, &
!               INT(ib,isp),comm,reqs(ireq),errcode)
!             !ENDIF
!
!             !IF (k.gt.0) THEN
!               ireq = ireq + 1
!               CALL MPI_Isend(bufsend(1:k,1:8,ib,is),8_isp*k,mpidbl,dest,INT(ib,isp), &
!               comm, reqs(ireq), errcode)
!               !CALL MPI_SENDRECV(bufsend(1:k,1:8,ib,is), 8_isp*k, mpidbl, dest, tag, &
!               !recvbuf(npos,1), 1_isp, typebuffer, src, tag, comm, status, errcode)
!
!             !ENDIF
!
!             CALL MPI_TYPE_FREE(typebuffer,errcode)
!
!             npos = npos + j
!             !print*,ib,npos,j,iz,iy,ix
!
!           ENDDO
!         ENDDO
!       ENDDO
!       CALL MPI_Waitall(ireq,reqs,MPI_STATUSES_IGNORE,errcode)

! Wait all is done inside the loop
      npos=1
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1

            IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE

            ib = 2+ix + (1+iy)*3 + (1+iz)*9

            k = mpi_npart(ib,is)
            j = nrecv_buf(ib,is)

            dest = INT(neighbour(ix,iy,iz),isp)
            src  = INT(neighbour(-ix,-iy,-iz),isp)

            CALL MPI_TYPE_VECTOR(8_isp,INT(j,isp),INT(nrecv_buf_tot+1,isp),MPI_DOUBLE_PRECISION,typebuffer,errcode)
            call MPI_TYPE_COMMIT(typebuffer,errcode)
            ! Exchange
            CALL MPI_Irecv(recvbuf(npos,1),1_isp, typebuffer,src, &
            INT(ib,isp),comm,reqs(1),errcode)
            CALL MPI_Isend(bufsend(1:k,1:8,ib,is),INT(8_idp*k,isp),mpidbl,dest,INT(ib,isp), &
            comm, reqs(2), errcode)
            CALL MPI_Waitall(2_isp,reqs,MPI_STATUSES_IGNORE,errcode)

            !CALL MPI_SENDRECV(bufsend(1:k,1:8,ib,is), 8_isp*k, mpidbl, dest, tag, &
            !recvbuf(npos,1), 1_isp, typebuffer, src, tag, comm, status, errcode)

            call MPI_TYPE_FREE(typebuffer,errcode)
            npos = npos + j
            !print*,ib,npos,j,iz,iy,ix

          ENDDO
        ENDDO
      ENDDO

      ! ________________________________________________

      ENDIF

#if defined(DEBUG) && (DEBUG==3)
	write(0,*) " Copy of the buffers in the tile particle arrays"
#endif

      ! If the buffer is not empty...
      IF (nrecv_buf_tot.gt.0) THEN

      ! Get first tile dimensions (may be different from last tile)
      nx0_grid_tile_dx = 1._num/(curr%array_of_tiles(1,1,1)%nx_grid_tile*dx)
      ny0_grid_tile_dy = 1._num/(curr%array_of_tiles(1,1,1)%ny_grid_tile*dy)
      nz0_grid_tile_dz = 1._num/(curr%array_of_tiles(1,1,1)%nz_grid_tile*dz)

      ! Parallelization on the tiles:
      ! Each tile will read the buffer and copy particles which belong
      ! to it into its own particle array

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(is,ntilez,ntiley,ntilex,nrecv_buf_tot,recvbuf,curr,dxs2,dys2,dzs2, &
      !$OMP nx0_grid_tile_dx,ny0_grid_tile_dy,nz0_grid_tile_dz, &
      !$OMP x_min_local,y_min_local,z_min_local,zgrid) &
      !$OMP PRIVATE(ix,iy,iz,i,k,indx,indy,indz,curr_tile,nptile) &
      !$OMP COLLAPSE(3) SCHEDULE(runtime)
      DO iz=1, ntilez! LOOP ON TILES
			DO iy=1, ntiley
			DO ix=1, ntilex

				curr_tile=>curr%array_of_tiles(ix,iy,iz)
				!nptile=curr_tile%np_tile(1)

        ! If the tile is not at the boundary, no particle will be transfer
				IF (.NOT.curr_tile%subdomain_bound) CYCLE

        ! Add received particles to particle arrays
        DO i = 1,nrecv_buf_tot

        ! Get particle index in array of tile
		    indx = MIN(FLOOR((recvbuf(i,1)-x_min_local+dxs2)*(nx0_grid_tile_dx),idp)+1,ntilex)
		    indy = MIN(FLOOR((recvbuf(i,2)-y_min_local+dys2)*(ny0_grid_tile_dy),idp)+1,ntiley)
		    indz = MIN(FLOOR((recvbuf(i,3)-(z_min_local+zgrid)+dzs2)*(nz0_grid_tile_dz),idp)+1,ntilez)

! 				IF (((recvbuf(i,1) .LT. x_min_local) .OR. (recvbuf(i,1) .GE. x_max_local)) .OR. &
! 					 ((recvbuf(i,2) .LT. y_min_local) .OR. (recvbuf(i,2) .GE. y_max_local)) .OR. &
! 					((recvbuf(i,3) .LT. z_min_local) .OR. (recvbuf(i,3) .GE. z_max_local))) THEN
!
! 					print*,i,nrecv_buf_tot,indx,indy,indz, recvbuf(i,1),recvbuf(i,2),recvbuf(i,3),recvbuf(i,4)
!
! 				ENDIF

		    ! If the particle is in the current tile
		    IF ((indx.eq.ix).AND.(indy.eq.iy).AND.(indz.eq.iz)) THEN

          ! Sanity check for max number of particles in tile
          nptile = curr_tile%np_tile(1)+1 ! Current number of particles in the tile
          k  = curr_tile%npmax_tile  ! Max number of particles in the tile
          IF (nptile .GT. k) THEN
          ! Resize particle tile arrays if tile is full
        	curr%are_tiles_reallocated(ix,iy,iz)=1
            CALL resize_particle_arrays(curr_tile, k, NINT(resize_factor*k+1,idp))
          ENDIF

          ! Finally, add particle to tile
          curr_tile%np_tile(1)=nptile
          curr_tile%part_x(nptile)  = recvbuf(i,1)
          curr_tile%part_y(nptile)  = recvbuf(i,2)
          curr_tile%part_z(nptile)  = recvbuf(i,3)
          curr_tile%part_ux(nptile) = recvbuf(i,4)
          curr_tile%part_uy(nptile) = recvbuf(i,5)
          curr_tile%part_uz(nptile) = recvbuf(i,6)
          curr_tile%part_gaminv(nptile) = recvbuf(i,7)
          curr_tile%pid(nptile,wpid) = recvbuf(i,8)
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

	END SUBROUTINE particle_bcs_tiles_and_mpi_3d


END MODULE boundary

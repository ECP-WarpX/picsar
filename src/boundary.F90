!!!! --- MODULE COUNTAINING ROUTINE FOR BOUNDARY CONDITIONS ON FIELDS, CURRENTS AND PARTICLES
!!!! --- For the moment this module handles:
!!!! ---  periodic external boundary conditions for particles and fields
! List of suboutines:
! - field_bc
! - exchange_mpi_3d_grid_array_with_guards
! - exchange_mpi_3d_grid_array_with_guards_nonblocking
! - summation_bcs
! - summation_bcs_nonblocking
! - particle_bcs_tiles_and_mpi_3d

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

    CALL exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, nx, ny, nz)

  END SUBROUTINE field_bc

!!! --- ROUTINE EXCHANGING GUARD REGIONS BETWEEN SUBDOMAINS (BLOCKING VERSION+DIAGONAL TRICK)
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarray, basetype, sz, szmax, i, j, k, n
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

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    ! MOVE EDGES ALONG X
    CALL MPI_SENDRECV(field(0,-nyg,-nzg), 1_isp, subarray, INT(proc_x_min,isp), &
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

    CALL MPI_SENDRECV(field(nx_local+1-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_x_max,isp), &
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

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    ! MOVE EDGES ALONG Y
    CALL MPI_SENDRECV(field(-nxg,0,-nzg), 1_isp, subarray, INT(proc_y_min,isp), &
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

    CALL MPI_SENDRECV(field(-nxg,ny_local+1-nyg,-nzg), 1_isp, subarray, INT(proc_y_max,isp), &
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

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    ! MOVE EDGES ALONG Z
    CALL MPI_SENDRECV(field(-nxg,-nyg,0), 1_isp, subarray, INT(proc_z_min,isp), &
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

    CALL MPI_SENDRECV(field(-nxg,-nyg,nz_local+1-nzg), 1_isp, subarray, INT(proc_z_max,isp), &
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

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards


!!! --- ROUTINE EXCHANGING GUARD REGIONS BETWEEN SUBDOMAINS (NON-BLOCKING VERSION+ DIAGONAL TRICK)
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)
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

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_ISEND(field(0,-nyg,-nzg), 1_isp, subarray, INT(proc_x_min,isp), tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(nx_local+1,-nyg,-nzg), 1_isp, subarray, INT(proc_x_max,isp), tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(nx_local+1-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_x_max,isp), tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_x_min,isp), tag, &
        comm, requests(4), errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

    ! NEED TO WAIT BEFORE EXCHANGING ALONG Y (DIAGONAL TERMS)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! MOVE EDGES ALONG Y
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_ISEND(field(-nxg,0,-nzg), 1_isp, subarray, INT(proc_y_min,isp), tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg,ny_local+1,-nzg), 1_isp, subarray, INT(proc_y_max,isp), tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(-nxg,ny_local+1-nyg,-nzg), 1_isp, subarray, INT(proc_y_max,isp), tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_y_min,isp), tag, &
        comm, requests(4), errcode)

    ! NEED TO WAIT BEFORE EXCHANGING ALONG Z (DIAGONAL TERMS)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)


    CALL MPI_TYPE_FREE(subarray, errcode)

   ! MOVE EDGES ALONG Z
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_ISEND(field(-nxg,-nyg,0), 1_isp, subarray, INT(proc_z_min,isp), tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,nz_local+1), 1_isp, subarray, INT(proc_z_max,isp), tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(-nxg,-nyg,nz_local+1-nzg), 1_isp, subarray, INT(proc_z_max,isp), tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_z_min,isp), tag, &
        comm, requests(4), errcode)

    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)
    CALL MPI_TYPE_FREE(subarray, errcode)

  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking




!!! --- Routine for adding current contributions fron adjacent subdomains
  SUBROUTINE summation_bcs(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarray, nn, sz, i

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn+1,-nyg,-nzg), 1_isp, subarray, &
        INT(neighbour( 1,0,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(-1,0,0),isp), tag, comm, status, errcode)
    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg), 1_isp, subarray, &
        INT(neighbour(-1,0,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour( 1,0,0),isp), tag, comm, status, errcode)
    array(nn+1-nxg:nn,:,:) = array(nn+1-nxg:nn,:,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,nn+1,-nzg), 1_isp, subarray, &
        INT(neighbour(0, 1,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0,-1,0),isp), tag, comm, status, errcode)
    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg), 1_isp, subarray, &
        INT(neighbour(0,-1,0),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0, 1,0),isp), tag, comm, status, errcode)
    array(:,nn+1-nyg:nn,:) = array(:,nn+1-nyg:nn,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,nn+1), 1_isp, subarray, &
        INT(neighbour(0,0, 1),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0,0,-1),isp), tag, comm, status, errcode)
    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg), 1_isp, subarray, &
        INT(neighbour(0,0,-1),isp), tag, temp, sz, mpidbl, &
        INT(neighbour(0,0, 1),isp), tag, comm, status, errcode)
    array(:,:,nn+1-nzg:nn) = array(:,:,nn+1-nzg:nn) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

  END SUBROUTINE summation_bcs

!!! --- Routine for adding current contributions fron adjacent subdomains nonblocking version
  SUBROUTINE summation_bcs_nonblocking(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarray, nn, sz, i
    INTEGER(isp) :: requests(4)

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(nn,-nyg,-nzg), 1_isp, subarray, INT(proc_x_max,isp), tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_x_min,isp), tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_x_min,isp), tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_x_max,isp), tag, &
    comm, requests(4), errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg,:,:) = array(0:nxg,:,:) + temp1
    array(nn-nxg:nn,:,:) = array(nn-nxg:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)
    nn = ny

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(-nxg,nn,-nzg), 1_isp, subarray, INT(proc_y_max,isp), tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_y_min,isp), tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_y_min,isp), tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_y_max,isp), tag, &
    comm, requests(4), errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(:,0:nyg,:) = array(:,0:nyg,:) + temp1
    array(:,nn-nyg:nn,:) = array(:,nn-nyg:nn,:) + temp2

    DEALLOCATE(temp1,temp2)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
    nn = nz

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)),temp2(subsizes(1), subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(-nxg,-nyg,nn), 1_isp, subarray, INT(proc_z_max,isp), tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_z_min,isp), tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg), 1_isp, subarray, INT(proc_z_min,isp), tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_z_max,isp), tag, &
    comm, requests(4), errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(:,:,0:nzg) = array(:,:,0:nzg) + temp1
    array(:,:,nn-nzg:nn) = array(:,:,nn-nzg:nn) + temp2

    DEALLOCATE(temp1,temp2)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

  END SUBROUTINE summation_bcs_nonblocking

!!! --- Routine for adding current contributions fron adjacent subdomains persistent version
  SUBROUTINE summation_bcs_persistent_jx(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
    USE communications

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarrayx, subarrayy, subarrayz, nn, sz, i
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
      subsizes(1) = nxg+1
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max_mpisp, tag, comm, reqperjxx(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayx, proc_x_min_mpisp, tag, comm, reqperjxx(3), errcode)     

      call MPI_RECV_INIT(temp1, sz, mpidbl, proc_x_min_mpisp, tag, comm, reqperjxx(2), errcode)
      call MPI_RECV_INIT(temp2, sz, mpidbl, proc_x_max_mpisp, tag, comm, reqperjxx(4), errcode)  
      
      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max_mpisp, tag, comm, reqperjxy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayy, proc_y_min_mpisp, tag, comm, reqperjxy(3), errcode) 
          
      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max_mpisp, tag, comm, reqperjxz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayz, proc_z_min_mpisp, tag, comm, reqperjxz(3), errcode) 

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)  
      CALL MPI_TYPE_FREE(subarrayz, errcode)
    
    ENDIF
     
    !! -- Summation along X- direction    

    subsizes(1) = nxg+1
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

    array(0:nxg,:,:) = array(0:nxg,:,:) + temp1
    array(nn-nxg:nn,:,:) = array(nn-nxg:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)   

    !! -- Summation along Y- direction    

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
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

    array(:,0:nyg,:) = array(:,0:nyg,:) + temp1
    array(:,nn-nyg:nn,:) = array(:,nn-nyg:nn,:) + temp2

    DEALLOCATE(temp1,temp2)
    
    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
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
    
    array(:,:,0:nzg) = array(:,:,0:nzg) + temp1
    array(:,:,nn-nzg:nn) = array(:,:,nn-nzg:nn) + temp2
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

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarrayx, subarrayy, subarrayz, nn, sz, i
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
      subsizes(1) = nxg+1
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max_mpisp, tag, comm, reqperjyx(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayx, proc_x_min_mpisp, tag, comm, reqperjyx(3), errcode)     
      
      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max_mpisp, tag, comm, reqperjyy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayy, proc_y_min_mpisp, tag, comm, reqperjyy(3), errcode) 
          
      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max_mpisp, tag, comm, reqperjyz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayz, proc_z_min_mpisp, tag, comm, reqperjyz(3), errcode) 

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)  
      CALL MPI_TYPE_FREE(subarrayz, errcode)
    
    ENDIF
     
    !! -- Summation along X- direction    

    subsizes(1) = nxg+1
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

    array(0:nxg,:,:) = array(0:nxg,:,:) + temp1
    array(nn-nxg:nn,:,:) = array(nn-nxg:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)   

    !! -- Summation along Y- direction    

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
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

    array(:,0:nyg,:) = array(:,0:nyg,:) + temp1
    array(:,nn-nyg:nn,:) = array(:,nn-nyg:nn,:) + temp2

    DEALLOCATE(temp1,temp2)
    
    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
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
    
    array(:,:,0:nzg) = array(:,:,0:nzg) + temp1
    array(:,:,nn-nzg:nn) = array(:,:,nn-nzg:nn) + temp2
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

    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: subarrayx, subarrayy, subarrayz, nn, sz, i
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
      subsizes(1) = nxg+1
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max_mpisp, tag, comm, reqperjzx(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayx, proc_x_min_mpisp, tag, comm, reqperjzx(3), errcode)     
      
      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max_mpisp, tag, comm, reqperjzy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayy, proc_y_min_mpisp, tag, comm, reqperjzy(3), errcode) 
          
      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max_mpisp, tag, comm, reqperjzz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayz, proc_z_min_mpisp, tag, comm, reqperjzz(3), errcode) 

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)  
      CALL MPI_TYPE_FREE(subarrayz, errcode)
    
    ENDIF
     
    !! -- Summation along X- direction    

    subsizes(1) = nxg+1
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

    array(0:nxg,:,:) = array(0:nxg,:,:) + temp1
    array(nn-nxg:nn,:,:) = array(nn-nxg:nn,:,:) + temp2

    DEALLOCATE(temp1,temp2)   

    !! -- Summation along Y- direction    

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
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

    array(:,0:nyg,:) = array(:,0:nyg,:) + temp1
    array(:,nn-nyg:nn,:) = array(:,nn-nyg:nn,:) + temp2

    DEALLOCATE(temp1,temp2)
    
    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
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
    
    array(:,:,0:nzg) = array(:,:,0:nzg) + temp1
    array(:,:,nn-nzg:nn) = array(:,:,nn-nzg:nn) + temp2
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
    tmptime = MPI_WTIME()
    ! Electric field MPI exchange between subdomains
    CALL field_bc(ex, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(ey, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(ez, nxguards, nyguards, nzguards, nx, ny, nz)
    localtimes(8) = localtimes(8) + (MPI_WTIME() - tmptime)
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
    tmptime = MPI_WTIME()
    ! Magnetic field MPI exchange between subdomains
    CALL field_bc(bx, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(by, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(bz, nxguards, nyguards, nzguards, nx, ny, nz)
    localtimes(6) = localtimes(6) + (MPI_WTIME() - tmptime)
#if defined(DEBUG)
  WRITE(0,*) "bfield_bcs: stop"
#endif    
  END SUBROUTINE bfield_bcs

!!! --- Boundary conditions routine for currents
  SUBROUTINE current_bcs
    REAL(num) :: tmptime
    tmptime = MPI_WTIME()
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
    
    localtimes(4) = localtimes(4) + (MPI_WTIME() - tmptime)
  END SUBROUTINE current_bcs

!!! --- Boundary conditions routine for charge density
SUBROUTINE charge_bcs
! Add charge contribution from adjacent subdomains
! __________________________________________________
    USE time_stat


    REAL(num) :: tmptime
    tmptime = MPI_WTIME() 

    CALL summation_bcs_nonblocking(rho, nxjguards, nyjguards, nzjguards, nx, ny, nz)
    
    localtimes(13) = localtimes(13) + (MPI_WTIME() - tmptime) 
    
END SUBROUTINE charge_bcs


!!! Boundary condition routine on particles
  SUBROUTINE particle_bcs
  	USE omp_lib
    USE time_stat
    IMPLICIT NONE
	  REAL(num) :: tdeb, tend
    REAL(num) :: tmptime
    tmptime = MPI_WTIME()
    tdeb=MPI_WTIME()
    
    
    IF (partcom.eq.1) THEN
    
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

    localtimes(11) = localtimes(11) + (MPI_WTIME() - tmptime)
    tmptime = MPI_WTIME()
    tend = MPI_WTIME()
    local_time_part=local_time_part+(tend-tdeb)

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_mpi: start"
#endif

    ! Then exchange particle between MPI domains
    CALL particle_bcs_mpi_blocking
    
#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_mpi: stop"
#endif

    localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)
    
    ! Tile and MPI com in one
    ! Default
    ELSE

#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_tiles_and_mpi: start"
#endif
    
      tmptime = MPI_WTIME()
    
      CALL particle_bcs_tiles_and_mpi_3d
      
      localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)
    
#if defined(DEBUG)
  WRITE(0,*) "particle_bcs_tiles_and_mpi: stop"
#endif    
    ENDIF
    
  END SUBROUTINE particle_bcs

!!! Boundary condition on tiles
  SUBROUTINE particle_bcs_tiles
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile, curr_tile_add
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: test =0

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
                        IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max))    &
                        .AND. ((party .GE. curr_tile%y_tile_min) .AND. (party .LT. curr_tile%y_tile_max))  &
                        .AND. ((partz .GE. curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) &
                        CYCLE

                        ! Case 2: if particle left MPI domain nothing to do now
                        IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
                        IF ((party .LT. y_min_local) .OR. (party .GE. y_max_local)) CYCLE
                        IF ((partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) CYCLE

                        ! Case 3: particles changed tile. Tranfer particle to new tile
                        ! Get new indexes of particle in array of tiles
                        indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex)
                        indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy))+1,ntiley)
                        indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez)
                        CALL rm_particle_at_tile(curr,ix,iy,iz,i)
                        CALL add_particle_at_tile(curr, indx,indy,indz, &
                             partx, party, partz, partux, partuy, partuz, gaminv, partw)
                    END DO !END LOOP ON PARTICLES
                END DO
            END DO
        END DO ! END LOOP ON TILES
    END DO ! END LOOP ON SPECIES
  END SUBROUTINE particle_bcs_tiles

!!! Boundary condition on tiles - 3D version 
!!! This version is efficient when the number of tiles is large 
!!! compared to the number of threads 
  SUBROUTINE particle_bcs_tiles_openmp()
    USE omp_lib 
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile, curr_tile_add
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: test =0, nthreads_tot, nthreads_loop1, nthreads_loop2
	
#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif
	
	IF (nthreads_tot .GT. 1) THEN 
		nthreads_loop1=MIN(nspecies,nthreads_tot)
		nthreads_loop2=MAX(1,nthreads_tot/nthreads_loop1)
	ELSE 
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF 
	
	!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,ispecies, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz) &
	!$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, & 
    !$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz) NUM_THREADS(nthreads_loop1) 
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
					!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile)  &
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
									.AND. ((partz .GE. curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) &
									CYCLE

									! Case 2: if particle left MPI domain nothing to do now
									IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
									IF ((party .LT. y_min_local) .OR. (party .GE. y_max_local)) CYCLE
									IF ((partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) CYCLE

									! Case 3: particles changed tile. Tranfer particle to new tile
									! Get new indexes of particle in array of tiles
									indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex)
									indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy))+1,ntiley)
									indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez)
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


!!! Boundary condition on tiles
  SUBROUTINE particle_bcs_tiles_2d
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile, curr_tile_add
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: test =0

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
					.AND. ((partz .GE. curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) &
					CYCLE

					! Case 2: if particle left MPI domain nothing to do now
					IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
					IF ((partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) CYCLE

					! Case 3: particles changed tile. Tranfer particle to new tile
					! Get new indexes of particle in array of tiles
					indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex)
					indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez)
					CALL rm_particle_at_tile(curr,ix,iy,iz,i)
					CALL add_particle_at_tile(curr, indx,iy,indz, &
						 partx, party, partz, partux, partuy, partuz, gaminv, partw)
				END DO !END LOOP ON PARTICLES
			END DO
        END DO ! END LOOP ON TILES
    END DO ! END LOOP ON SPECIES
  END SUBROUTINE particle_bcs_tiles_2d
  
!!! Boundary condition on tiles - 2D version 
!!! This version is efficient when the number of tiles is large 
!!! compared to the number of threads 
  SUBROUTINE particle_bcs_tiles_2d_openmp()
    USE omp_lib 
    IMPLICIT NONE
    INTEGER(idp):: i, ispecies, ix, iy, iz, indx, indy, indz, ipx, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile, curr_tile_add
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: test =0, nthreads_tot, nthreads_loop1, nthreads_loop2
	
#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif
	
	IF (nthreads_tot .GT. 1) THEN 
		nthreads_loop1=MIN(nspecies,nthreads_tot)
		nthreads_loop2=MAX(1,nthreads_tot/nthreads_loop1)
	ELSE 
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF 
	iy=1
	!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,ispecies, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipz) &
	!$OMP SHARED(iy,nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, & 
    !$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz) NUM_THREADS(nthreads_loop1) 
    DO ispecies=1, nspecies ! LOOP ON SPECIES
        curr=> species_parray(ispecies)
        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
        DO ipz=1,3
        	DO ipx=1,3
				!$OMP PARALLEL DO DEFAULT(NONE) SHARED(iy,curr,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, & 
				!$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile)  &
				!$OMP FIRSTPRIVATE(ipx,ipz) &
				!$OMP PRIVATE(ix,iz,i,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
				!$OMP indx,indy,indz) COLLAPSE(2) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2)  
				DO iz=ipz, ntilez,3! LOOP ON TILES
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
							.AND. ((partz .GE. curr_tile%z_tile_min) .AND. (partz .LT. curr_tile%z_tile_max))) &
							CYCLE

							! Case 2: if particle left MPI domain nothing to do now
							IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) CYCLE
							IF ((partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) CYCLE

							! Case 3: particles changed tile. Tranfer particle to new tile
							! Get new indexes of particle in array of tiles
							indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex)
							indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez)
							CALL rm_particle_at_tile(curr,ix,iy,iz,i)
							CALL add_particle_at_tile(curr, indx,iy,indz, &
								 partx, party, partz, partux, partuy, partuz, gaminv, partw)
						END DO !END LOOP ON PARTICLES
					END DO
				END DO ! END LOOP ON TILES
				!$OMP END PARALLEL DO 
        	END DO 
        END DO
    END DO ! END LOOP ON SPECIES
    !$OMP END PARALLEL DO 
  END SUBROUTINE particle_bcs_tiles_2d_openmp

  ! __________________________________________________________________
  SUBROUTINE particle_bsc_openmp_reordering
  ! Experimental subroutine
  ! This subroutine process the particle boundary conditions
  ! For this aim, the particles are reordered into buckets in the particle arrays
  ! Each bucket corresponds to the group of particle to be exchanged in a common direction
  ! Exchanges between tiles and between MPI domains is combined for better efficiency
  ! This subroutine is less efficient that particle_bcs_tiles_openmp
  ! __________________________________________________________________

    USE omp_lib 
    USE communications
    IMPLICIT NONE
    
    INTEGER(idp):: i, is, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile, curr_tile_add
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
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
		nthreads_loop2=MAX(1,nthreads_tot/nthreads_loop1)
	ELSE 
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF 
	
	ALLOCATE(buffer(ntilex,ntiley,ntilez,nspecies))
	
	!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,is, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz) &
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

									indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex)
									indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy))+1,ntiley)
									indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez)

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

	  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(curr,is, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz) &
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
  
  
!!! MPI Boundary condition routine on particles
  SUBROUTINE particle_bcs_mpi_blocking
    INTEGER(isp), PARAMETER :: nvar=8 ! Simple implementation
    INTEGER(isp), DIMENSION(-1:1,-1:1,-1:1) :: nptoexch
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: temp
    LOGICAL(idp), ALLOCATABLE, DIMENSION(:) :: mask
    INTEGER(isp) :: ibuff, isend, nout, nbuff, ninit
    INTEGER(isp) :: xbd, ybd, zbd
    INTEGER(isp) :: ixp, iyp, izp
    INTEGER(isp) :: nsend_buf, nrecv_buf, npart_curr
    INTEGER(isp) :: dest, src
    LOGICAL(idp) :: out_of_bounds
    INTEGER(idp) :: ispecies, i, ip, ix, iy, iz
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
        nbuff=currsp%species_npart*nvar
        ibuff=1
        ALLOCATE(sendbuf(-1:1,-1:1,-1:1,1:nbuff))
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
                            	CASE DEFAULT ! periodic 
                                	curr%part_y(i) = part_xyz - length_y
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
                            mask(i)=.FALSE.
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
                            		mask(i)=.FALSE.
                            		CYCLE
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


  ! ______________________________________________________________________________________
  SUBROUTINE particle_bcs_tiles_and_mpi_3d
  ! This subroutine combined in a single routine the particle communications between tiles
  ! and between MPI domains for 3D
  ! 
  ! Mathieu Lobet, May 2016
  ! ______________________________________________________________________________________
    USE omp_lib 
    USE communications
    USE precomputed
    IMPLICIT NONE
    INTEGER(idp) :: i, is, ix, iy, iz, indx, indy, indz, ipx, ipy, ipz
    INTEGER(idp) :: xbd,ybd,zbd
    INTEGER(idp) :: k,j,ib,ibs
    INTEGER(isp) :: dest, src
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER    :: curr_tile, curr_tile_add
    REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
    INTEGER(idp) :: test =0, nthreads_tot, nthreads_loop1, nthreads_loop2
    INTEGER(idp), dimension(:,:), ALLOCATABLE :: mpi_npart
    REAL(num), dimension(:,:,:,:), ALLOCATABLE :: bufsend
    REAL(num), dimension(:,:), ALLOCATABLE     :: recvbuf
    TYPE(mpi_buffer), dimension(:,:,:,:), ALLOCATABLE :: tilebuf
    INTEGER(isp), DIMENSION(:,:), ALLOCATABLE :: nrecv_buf
    INTEGER(isp) :: nrecv_buf_tot,npos, typebuffer
    REAL(num) :: nx0_grid_tile_dx, ny0_grid_tile_dy, nz0_grid_tile_dz
    INTEGER(isp) :: reqs(2),stats(2)
    INTEGER(idp) :: recvbuf_index(27)
	
	  ! _________________________________________________________
	  ! Determine number of threads to be used for nested parallel region
	
#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif
	
	IF (nthreads_tot .GT. 1) THEN 
		nthreads_loop1=MIN(nspecies,nthreads_tot)
		nthreads_loop2=MAX(1,nthreads_tot/nthreads_loop1)
	ELSE 
		nthreads_loop1=1
		nthreads_loop2=1
	ENDIF 
	
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
		
	ALLOCATE(mpi_npart(27,nspecies))
	ALLOCATE(tilebuf(ntilex,ntiley,ntilez,nspecies))
	
	  !$OMP PARALLEL DO DEFAULT(NONE) &
	  !$OMP PRIVATE(curr,is,ib,k,nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz,&
	  !$OMP nx0_grid_tile_dx,ny0_grid_tile_dy,nz0_grid_tile_dz) &
	  !$OMP SHARED(nspecies,nthreads_loop2,species_parray,ntilex,ntiley,ntilez,x_min_local,y_min_local,z_min_local, & 
	  !$OMP length_x,length_y,length_z,dxs2,dys2,dzs2, &
	  !$OMP x_min_boundary,x_max_boundary,y_min_boundary,y_max_boundary,z_min_boundary,z_max_boundary,  &
	  !$OMP pbound_x_min,pbound_x_max,pbound_y_min,pbound_y_max,pbound_z_min,pbound_z_max, &
    !$OMP x_max_local,y_max_local,z_max_local,dx,dy,dz,mpi_npart,tilebuf) &
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
					!$OMP nx0_grid_tile_dx,ny0_grid_tile_dy,nz0_grid_tile_dz,dxs2,dys2,dzs2)  &
					!$OMP FIRSTPRIVATE(ipx,ipy,ipz,is) &
					!$OMP PRIVATE(ix,iy,iz,i,ib,k,curr_tile,nptile,partx,party,partz,partux,partuy,partuz,gaminv,partw, &
					!$OMP indx,indy,indz,xbd,ybd,zbd)  &
					!$OMP COLLAPSE(3) SCHEDULE(runtime) NUM_THREADS(nthreads_loop2) 
					! LOOP ON TILES
					DO iz=ipz, ntilez,3
						DO iy=ipy, ntiley,3
							DO ix=ipx, ntilex,3
								curr_tile=>curr%array_of_tiles(ix,iy,iz)
								nptile=curr_tile%np_tile(1)
								
								! Allocation of the buffer
								IF (curr_tile%subdomain_bound) THEN
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_x(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_y(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_z(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_ux(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_uy(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_uz(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%part_gaminv(2000,27))
								  ALLOCATE(tilebuf(ix,iy,iz,is)%pid(2000,27))
								ENDIF
								tilebuf(ix,iy,iz,is)%npart(1:27) = 0
								
								! LOOP ON PARTICLES
								DO i=nptile, 1, -1
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
									CYCLE

									! Case 2: if particle left MPI domain
									IF (((partx .LT. x_min_local) .OR. (partx .GE. x_max_local)) .OR. &
									 ((party .LT. y_min_local) .OR. (party .GE. y_max_local)) .OR. &
									 ((partz .LT. z_min_local) .OR. (partz .GE. z_max_local))) THEN
									
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
                    IF (partz .LT. z_min_local) THEN
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
                    ELSE IF (partz .GE. z_max_local) THEN
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
                    
                    ! Update of the number of particles
                    tilebuf(ix,iy,iz,is)%npart(ib) = k
                    
                    ! The particle is deleted
                    CALL rm_particle_at_tile(curr,ix,iy,iz,i)

									  ENDIF
									  
									  CYCLE
                    
                  ENDIF

									! Case 3: particles changed tile. Tranfer particle to new tile
									! Get new indexes of particle in array of tiles
									indx = MIN(FLOOR((partx-x_min_local+dxs2)*(nx0_grid_tile_dx))+1,ntilex)
									indy = MIN(FLOOR((party-y_min_local+dys2)*(ny0_grid_tile_dy))+1,ntiley)
									indz = MIN(FLOOR((partz-z_min_local+dzs2)*(nz0_grid_tile_dz))+1,ntilez)
									!if ((indx.eq.0).or.(indy.eq.0).or.(indz.eq.0)) THEN
                   !print*,'xmin',x_min_local,'xmax',x_max_local,'x',partx,xbd
                   !print*,'ymin',y_min_local,'ymax',y_max_local,'y',party,ybd
                   !print*,'zmin',z_min_local,'zmax',z_max_local,'z',partz,zbd
                  !endif
									CALL rm_particle_at_tile(curr,ix,iy,iz,i)
									CALL add_particle_at_tile(curr, indx,indy,indz, &
										 partx, party, partz, partux, partuy, partuz, gaminv, partw)
								END DO !END LOOP ON PARTICLES
								
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
	  !$OMP PRIVATE(curr,is, nx0_grid_tile,ny0_grid_tile,nz0_grid_tile,ipx,ipy,ipz) &
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
    ! Multiple thread version
    
    ALLOCATE(nrecv_buf(27,nspecies))
    nrecv_buf(:,:)=0
    
    IF (.FALSE.) THEN
    DO is=1, nspecies ! LOOP ON SPECIES
      !curr=> species_parray(is)
      
     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(is,mpi_npart,comm,neighbour,nrecv_buf) &
     !$OMP FIRSTPRIVATE(tag,status,stats,reqs) &
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
                        
               CALL MPI_Waitall(2_isp,reqs,stats,errcode)            
                        
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
    DO is=1, nspecies ! LOOP ON SPECIES
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
    ENDDO
    
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
              
           CALL MPI_Isend(bufsend(1:k,1:8,ib,is),8_isp*k,mpidbl,dest,INT(ib,isp), &
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
               
            CALL MPI_Isend(bufsend(1:k,1:8,ib,is),8_isp*k,mpidbl,dest,INT(ib,isp), &
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
      !$OMP x_min_local,y_min_local,z_min_local) &
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
		    indx = MIN(FLOOR((recvbuf(i,1)-x_min_local+dxs2)*(nx0_grid_tile_dx))+1,ntilex)
		    indy = MIN(FLOOR((recvbuf(i,2)-y_min_local+dys2)*(ny0_grid_tile_dy))+1,ntiley)
		    indz = MIN(FLOOR((recvbuf(i,3)-z_min_local+dzs2)*(nz0_grid_tile_dz))+1,ntilez)
		    
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
    
  END SUBROUTINE particle_bcs_tiles_and_mpi_3d



END MODULE boundary

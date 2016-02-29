!!!! --- MODULE COUNTAINING ROUTINE FOR BOUNDARY CONDITIONS ON FIELDS, CURRENTS AND PARTICLES
!!!! --- For the moment this module handles:
!!!! ---  periodic external boundary conditions for particles and fields
!!!! Suboutines:
!!!! field_bc
!!!! exchange_mpi_3d_grid_array_with_guards
!!!! exchange_mpi_3d_grid_array_with_guards_nonblocking
!!!! summation_bcs
!!!! summation_bcs_nonblocking

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
    CALL MPI_SENDRECV(field(0,-nyg,-nzg), 1_isp, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

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

    CALL MPI_SENDRECV(field(nx_local+1-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

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
    CALL MPI_SENDRECV(field(-nxg,0,-nzg), 1_isp, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

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

    CALL MPI_SENDRECV(field(-nxg,ny_local+1-nyg,-nzg), 1_isp, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

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
    CALL MPI_SENDRECV(field(-nxg,-nyg,0), 1_isp, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

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

    CALL MPI_SENDRECV(field(-nxg,-nyg,nz_local+1-nzg), 1_isp, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

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

    CALL MPI_ISEND(field(0,-nyg,-nzg), 1_isp, subarray, proc_x_min, tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(nx_local+1,-nyg,-nzg), 1_isp, subarray, proc_x_max, tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(nx_local+1-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_max, tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_min, tag, &
        comm, requests(4), errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

    ! NEED TO WAIT BEFORE EXCHANGING ALONG Y (DIAGONAL TERMS)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! MOVE EDGES ALONG Y
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_ISEND(field(-nxg,0,-nzg), 1_isp, subarray, proc_y_min, tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg,ny_local+1,-nzg), 1_isp, subarray, proc_y_max, tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(-nxg,ny_local+1-nyg,-nzg), 1_isp, subarray, proc_y_max, tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, subarray, proc_y_min, tag, &
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

    CALL MPI_ISEND(field(-nxg,-nyg,0), 1_isp, subarray, proc_z_min, tag, &
         comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,nz_local+1), 1_isp, subarray, proc_z_max, tag, &
        comm, requests(2), errcode)

    CALL MPI_ISEND(field(-nxg,-nyg,nz_local+1-nzg), 1_isp, subarray, proc_z_max, tag, &
         comm, requests(3), errcode)
    CALL MPI_IRECV(field(-nxg,-nyg,-nzg), 1_isp, subarray, proc_z_min, tag, &
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
    CALL MPI_ISEND(array(nn,-nyg,-nzg), 1_isp, subarray, proc_x_max, tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min, tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_min, tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max, tag, &
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
    CALL MPI_ISEND(array(-nxg,nn,-nzg), 1_isp, subarray, proc_y_max, tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min, tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg), 1_isp, subarray, proc_y_min, tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max, tag, &
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
    CALL MPI_ISEND(array(-nxg,-nyg,nn), 1_isp, subarray, proc_z_max, tag, &
    comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min, tag, &
    comm, requests(2), errcode)
    CALL MPI_ISEND(array(-nxg,-nyg,-nzg), 1_isp, subarray, proc_z_min, tag, &
    comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max, tag, &
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
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max, tag, comm, reqperjxx(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayx, proc_x_min, tag, comm, reqperjxx(3), errcode)     

      call MPI_RECV_INIT(temp1, sz, mpidbl, proc_x_min, tag, comm, reqperjxx(2), errcode)
      call MPI_RECV_INIT(temp2, sz, mpidbl, proc_x_max, tag, comm, reqperjxx(4), errcode)  
      
      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max, tag, comm, reqperjxy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayy, proc_y_min, tag, comm, reqperjxy(3), errcode) 
          
      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max, tag, comm, reqperjxz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayz, proc_z_min, tag, comm, reqperjxz(3), errcode) 

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
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min, tag, &
    comm, reqperjxx(2), errcode)  
    !call MPI_START(reqperjxx(2),errcode)  
    call MPI_START(reqperjxx(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max, tag, &
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
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min, tag, &
    comm, reqperjxy(2), errcode)   
    call MPI_START(reqperjxy(3),errcode) 
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max, tag, &
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

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min, tag, &
    comm, reqperjxz(2), errcode)
    call MPI_START(reqperjxz(1),errcode)
    !if (rank.eq.0) print*, 'irecv 1'
    
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max, tag, &
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
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max, tag, comm, reqperjyx(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayx, proc_x_min, tag, comm, reqperjyx(3), errcode)     
      
      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max, tag, comm, reqperjyy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayy, proc_y_min, tag, comm, reqperjyy(3), errcode) 
          
      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max, tag, comm, reqperjyz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayz, proc_z_min, tag, comm, reqperjyz(3), errcode) 

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
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min, tag, &
    comm, reqperjyx(2), errcode)    
    call MPI_START(reqperjyx(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max, tag, &
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
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min, tag, &
    comm, reqperjyy(2), errcode)   
    CALL MPI_START(reqperjyy(1),errcode) 
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max, tag, &
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

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min, tag, &
    comm, reqperjyz(2), errcode)
    call MPI_START(reqperjyz(1),errcode)
    !if (rank.eq.0) print*, 'irecv 1'
    
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max, tag, &
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
      call MPI_SEND_INIT(array(nn,-nyg,-nzg), 1_isp, subarrayx,  proc_x_max, tag, comm, reqperjzx(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayx, proc_x_min, tag, comm, reqperjzx(3), errcode)     
      
      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,nn,-nzg), 1_isp, subarrayy,  proc_y_max, tag, comm, reqperjzy(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayy, proc_y_min, tag, comm, reqperjzy(3), errcode) 
          
      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg,-nyg,nn), 1_isp, subarrayz,  proc_z_max, tag, comm, reqperjzz(1), errcode)
      call MPI_SEND_INIT(array(-nxg,-nyg,-nzg), 1_isp, subarrayz, proc_z_min, tag, comm, reqperjzz(3), errcode) 

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
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min, tag, &
    comm, reqperjzx(2), errcode)    
    call MPI_START(reqperjzx(1),errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max, tag, &
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
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min, tag, &
    comm, reqperjzy(2), errcode)   
    CALL MPI_START(reqperjzy(1),errcode) 
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max, tag, &
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

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min, tag, &
    comm, reqperjzz(2), errcode)
    call MPI_START(reqperjzz(1),errcode)
    !if (rank.eq.0) print*, 'irecv 1'
    
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max, tag, &
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
    tmptime = MPI_WTIME()
    ! Electric field MPI exchange between subdomains
    CALL field_bc(ex, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(ey, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(ez, nxguards, nyguards, nzguards, nx, ny, nz)
    localtimes(8) = localtimes(8) + (MPI_WTIME() - tmptime)
  END SUBROUTINE efield_bcs

!!! --- Boundary condition routine for magnetic field
  SUBROUTINE bfield_bcs
    REAL(num) :: tmptime
    tmptime = MPI_WTIME()
    ! Magnetic field MPI exchange between subdomains
    CALL field_bc(bx, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(by, nxguards, nyguards, nzguards, nx, ny, nz)
    CALL field_bc(bz, nxguards, nyguards, nzguards, nx, ny, nz)
    localtimes(6) = localtimes(6) + (MPI_WTIME() - tmptime)
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
    CALL summation_bcs_nonblocking(rho, nxjguards, nyjguards, nzjguards, nx, ny, nz)
END SUBROUTINE charge_bcs


!!! Boundary condition routine on particles
  SUBROUTINE particle_bcs
    USE time_stat
    IMPLICIT NONE

    REAL(num) :: tmptime
    
    tmptime = MPI_WTIME()

    ! First exchange particles between tiles (NO MPI at that point)
    CALL particle_bcs_tiles

    ! Then exchange particle between MPI domains
    CALL particle_bcs_mpi_blocking

    localtimes(2) = localtimes(2) + (MPI_WTIME() - tmptime)

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
                        IF (((partx .GE. curr_tile%x_tile_min) .AND. (partx .LT. curr_tile%x_tile_max)) &
                        .AND. ((party .GE. curr_tile%y_tile_min) .AND. (party .LT. curr_tile%y_tile_max)) &
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
                        CALL rm_particle_at_tile(curr_tile,i)
                        CALL add_particle_at_tile(curr, indx,indy,indz, &
                             partx, party, partz, partux, partuy, partuz, gaminv, partw)
                    END DO !END LOOP ON PARTICLES
                END DO
            END DO
        END DO ! END LOOP ON TILES
    END DO ! END LOOP ON SPECIES
  END SUBROUTINE particle_bcs_tiles


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
    INTEGER(isp) :: ispecies, i, ip, ix, iy, iz
    INTEGER(isp) :: ixtile, iytile, iztile
    REAL(num) :: part_xyz
    TYPE(particle_species), POINTER :: currsp
    TYPE(particle_tile), POINTER :: curr

    DO ispecies=1, nspecies !LOOP ON SPECIES
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
                    xbd = 0
                    ybd = 0
                    zbd = 0
                    part_xyz=0.
                    ! Identify outbounds particles
                    DO i = 1, curr%np_tile(1) !LOOP ON PARTICLES
                        xbd = 0
                        ybd = 0
                        zbd = 0
                        out_of_bounds = .FALSE.
                        part_xyz = curr%part_x(i)
                        ! Particle has left this processor
                        IF (part_xyz .LT. x_min_local) THEN
                            xbd = -1
                            IF (x_min_boundary) THEN
                                curr%part_x(i) = part_xyz + length_x
                            ENDIF
                        ENDIF
                        ! Particle has left this processor
                        IF (part_xyz .GE. x_max_local) THEN
                            xbd = 1
                            IF (x_max_boundary) THEN
                                curr%part_x(i) = part_xyz - length_x
                            ENDIF
                        ENDIF

                        part_xyz = curr%part_y(i)
                        ! Particle has left this processor
                        IF (part_xyz .LT. y_min_local) THEN
                            ybd = -1
                            IF (y_min_boundary) THEN
                                curr%part_y(i) = part_xyz + length_y
                            ENDIF
                        ENDIF

                        ! Particle has left this processor
                        IF (part_xyz .GE. y_max_local) THEN
                            ybd = 1
                            IF (y_max_boundary) THEN
                                curr%part_y(i) = part_xyz - length_y
                            ENDIF
                        ENDIF

                        part_xyz = curr%part_z(i)
                        ! Particle has left this processor
                        IF (part_xyz .LT. z_min_local) THEN
                            zbd = -1
                            IF (z_min_boundary) THEN
                                curr%part_z(i) = part_xyz + length_z
                            ENDIF
                        ENDIF

                        ! Particle has left this processor
                        IF (part_xyz .GE. z_max_local) THEN
                            zbd = 1
                            ! Particle has left the system
                            IF (z_max_boundary) THEN
                                curr%part_z(i) = part_xyz - length_z
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
                    CALL rm_particles_from_species(currsp, curr, mask)
                    DEALLOCATE(mask)
                  ENDDO
               ENDDO
            ENDDO ! END LOOP ON TILES
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

END MODULE boundary

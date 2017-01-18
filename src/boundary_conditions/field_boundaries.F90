! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)  
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
! FIELD_BOUNDARIES.F90
!
! Brief description:
! This file contains subroutines for boundary conditions on fields and currents.
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
! - field_bc
! - exchange_mpi_3d_grid_array_with_guards
! - exchange_mpi_3d_grid_array_with_guards_nonblocking
! - summation_bcs
! - summation_bcs_nonblocking
!
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Module that contains subroutines for boundary conditions on fields and currents.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
MODULE field_boundary
! ______________________________________________________________________________

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

  ! ____________________________________________________________________________
  !> @brief 
  !> Exchange field values at processor boundaries and apply field
  !> boundary conditions
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE field_bc(field, nxg, nyg, nzg, nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________
  
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: field
    IF (mpicom_curr.EQ.1) THEN
      CALL exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, nx, ny, nz)
    ELSE
    CALL exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, nx, ny, nz)
    ENDIF

  END SUBROUTINE field_bc

  ! ____________________________________________________________________________
  !> @brief
  !> Routine exchanging guard regions between subdomains (blocking, version+diagonal trick)
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________
  
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

  ! ____________________________________________________________________________
  !> @brief
  !> Routine exchanging guard regions between subdomains (non-blocking version 
  !> + diagonal trick)
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________

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


  ! ____________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE summation_bcs(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________

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



  ! ____________________________________________________________________________
  !> Routine for adding current contributions fron adjacent subdomains 
  ! nonblocking version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE summation_bcs_nonblocking(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________
  
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

  ! ____________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains persistent version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE summation_bcs_persistent_jx(array, nxg, nyg, nzg, &
                                         nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________
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

  ! ____________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains persistent version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE summation_bcs_persistent_jy(array, nxg, nyg, nzg, &
                                         nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________
  
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

  ! ____________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains persistent version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE summation_bcs_persistent_jz(array, nxg, nyg, nzg, &
                                         nx_local, ny_local, nz_local)
  ! ____________________________________________________________________________
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

  ! ____________________________________________________________________________
  !> @brief
  !> Boundary condition routine for electric field
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE efield_bcs
  ! ____________________________________________________________________________  
  
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

  ! ____________________________________________________________________________
  !> @brief
  !> Boundary condition routine for magnetic field
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE bfield_bcs
  ! ____________________________________________________________________________
  
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

  ! ____________________________________________________________________________
  !> @brief
  !> Boundary conditions routine for currents
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  SUBROUTINE current_bcs
  ! ____________________________________________________________________________
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

  ! ____________________________________________________________________________
  !> @brief
  !> Boundary conditions routine for charge density.
  !
  !> @details
  !> Add charge contribution from adjacent subdomains.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
SUBROUTINE charge_bcs
  ! ____________________________________________________________________________
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



  ! ______________________________________________________________________________
  !> @brief
  !> Move the global, local and tile boundaries when moving window
  !
  !> author
  !> Guillaume Blaclard
  !
  !> @date
  !> Creation: 2016
  !
  !> @param[in] dxx,dyy,dzz moving window displacement 
    SUBROUTINE pxr_move_sim_boundaries(dxx,dyy,dzz)
  ! ______________________________________________________________________________
    IMPLICIT NONE 
  
    REAL(num), INTENT(IN) :: dxx, dyy, dzz
    INTEGER(idp) :: ispecies, ix, iy, iz
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER    :: curr_tile

    ! MOVE GLOBAL BOUNDARIES 
    xmin=xmin+dxx
          xmax=xmax+dxx
          ! Update along Y
    ymin=ymin+dyy
          ymax=ymax+dyy
          ! Update along Z
    zmin=zmin+dzz
    zmax=zmax+dzz

    ! MOVE GLOBAL BOUNDARIES 
    x_grid_min=x_grid_min+dxx
          x_grid_max=x_grid_max+dxx
          ! Update along Y
    y_grid_min=y_grid_min+dyy
          y_grid_max=y_grid_max+dyy
          ! Update along Z
    z_grid_min=z_grid_min+dzz
    z_grid_max=z_grid_max+dzz

    ! MOVE LOCAL BOUNDARIES - MPI 
          ! Update along X
    x_min_local=x_min_local+dxx
          x_max_local=x_max_local+dxx
          ! Update along Y
    y_min_local=y_min_local+dyy
          y_max_local=y_max_local+dyy
          ! Update along Z
    z_min_local=z_min_local+dzz
          z_max_local=z_max_local+dzz

    ! MOVE LOCAL BOUNDARIES 
    x_grid_min_local=x_grid_min_local+dxx
          x_grid_max_local=x_grid_max_local+dxx
          ! Update along Y
    y_grid_min_local=y_grid_min_local+dyy
          y_grid_max_local=y_grid_max_local+dyy
          ! Update along Z
    z_grid_min_local=z_grid_min_local+dzz
    z_grid_max_local=z_grid_max_local+dzz

  ! MOVE TILE BOUNDARIES ALONG X,Y,Z 
    DO ispecies =1, nspecies
    curr=> species_parray(ispecies)
    DO iz=1, ntilez
           DO iy=1,ntiley
                  DO ix=1,ntilex
                          curr_tile=> curr%array_of_tiles(ix,iy,iz)
                          ! Update along X
                          curr_tile%x_grid_tile_min=curr_tile%x_grid_tile_min+dxx
                          curr_tile%x_grid_tile_max=curr_tile%x_grid_tile_max+dxx
                          ! Update along Y
                          curr_tile%y_grid_tile_min=curr_tile%y_grid_tile_min+dyy
                          curr_tile%y_grid_tile_max=curr_tile%y_grid_tile_max+dyy
                          ! Update along Z
                          curr_tile%z_grid_tile_min=curr_tile%z_grid_tile_min+dzz
                          curr_tile%z_grid_tile_max=curr_tile%z_grid_tile_max+dzz

                          ! Update along X
                          curr_tile%x_tile_min=curr_tile%x_tile_min+dxx
                          curr_tile%x_tile_max=curr_tile%x_tile_max+dxx
                          ! Update along Y
                          curr_tile%y_tile_min=curr_tile%y_tile_min+dyy
                          curr_tile%y_tile_max=curr_tile%y_tile_max+dyy
                          ! Update along Z
                          curr_tile%z_tile_min=curr_tile%z_tile_min+dzz
                          curr_tile%z_tile_max=curr_tile%z_tile_max+dzz
                    END DO 
                END DO
          END DO
      END DO

  END SUBROUTINE pxr_move_sim_boundaries


END MODULE field_boundary

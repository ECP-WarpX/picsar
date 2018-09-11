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
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Module that contains subroutines for boundary conditions on fields and currents.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
MODULE field_boundary
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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE field_bc(field, nxg, nyg, nzg, nx_local, ny_local, nz_local)
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: field
    IF (mpicom_curr.EQ.1) THEN
      CALL exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, nx, ny, nz)
    ELSE
      CALL exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg,   &
      nx, ny, nz)
    ENDIF

  END SUBROUTINE field_bc

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine exchanging guard regions between subdomains
  !> (blocking, version+diagonal trick)
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, nx_local,   &
    ny_local, nz_local)
    USE constants
    USE mpi
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpidbl

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    szmax = sizes(1) * sizes(2) * (nzg+1)
    sz = sizes(1) * sizes(3) * (nyg+1)
    IF (sz .GT. szmax) szmax = sz
    sz = sizes(2) * sizes(3) * (nxg+1)
    IF (sz .GT. szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    IF (is_dtype_init(1)) THEN
      mpi_dtypes(1) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
      is_dtype_init(1) = .FALSE.
    ENDIF

    ! MOVE EDGES ALONG X
    CALL MPI_SENDRECV(field(0, -nyg, -nzg), 1_isp, mpi_dtypes(1), INT(proc_x_min,     &
    isp), tag, temp, sz, basetype, INT(proc_x_max, isp), tag, comm, status, errcode)

    IF (proc_x_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
        DO j = -nyg, subsizes(2)-nyg-1
          DO i = nx_local, subsizes(1)+nx_local-1
            field(i, j, k) = temp(n)
            n = n + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(nx_local-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(1),          &
    INT(proc_x_max, isp), tag, temp, sz, basetype, INT(proc_x_min, isp), tag, comm,   &
    status, errcode)

    IF (proc_x_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
        DO j = -nyg, subsizes(2)-nyg-1
          DO i = -nxg, subsizes(1)-nxg-1
            field(i, j, k) = temp(n)
            n = n + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    IF (is_dtype_init(2)) THEN
      mpi_dtypes(2) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
      is_dtype_init(2) = .FALSE.
    ENDIF

    ! MOVE EDGES ALONG Y
    CALL MPI_SENDRECV(field(-nxg, 0, -nzg), 1_isp, mpi_dtypes(2), INT(proc_y_min,     &
    isp), tag, temp, sz, basetype, INT(proc_y_max, isp), tag, comm, status, errcode)

    IF (proc_y_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
        DO j = ny_local, subsizes(2)+ny_local-1
          DO i = -nxg, subsizes(1)-nxg-1
            field(i, j, k) = temp(n)
            n = n + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-nxg, ny_local-nyg, -nzg), 1_isp, mpi_dtypes(2),          &
    INT(proc_y_max, isp), tag, temp, sz, basetype, INT(proc_y_min, isp), tag, comm,   &
    status, errcode)

    IF (proc_y_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
        DO j = -nyg, subsizes(2)-nyg-1
          DO i = -nxg, subsizes(1)-nxg-1
            field(i, j, k) = temp(n)
            n = n + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF


    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    IF (is_dtype_init(3)) THEN
      mpi_dtypes(3) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
      is_dtype_init(3) = .FALSE.
    ENDIF


    ! MOVE EDGES ALONG Z
    CALL MPI_SENDRECV(field(-nxg, -nyg, 0), 1_isp, mpi_dtypes(3), INT(proc_z_min,     &
    isp), tag, temp, sz, basetype, INT(proc_z_max, isp), tag, comm, status, errcode)

    IF (proc_z_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = nz_local, subsizes(3)+nz_local-1
        DO j = -nyg, subsizes(2)-nyg-1
          DO i = -nxg, subsizes(1)-nxg-1
            field(i, j, k) = temp(n)
            n = n + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-nxg, -nyg, nz_local-nzg), 1_isp, mpi_dtypes(3),          &
    INT(proc_z_max, isp), tag, temp, sz, basetype, INT(proc_z_min, isp), tag, comm,   &
    status, errcode)

    IF (proc_z_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
        DO j = -nyg, subsizes(2)-nyg-1
          DO i = -nxg, subsizes(1)-nxg-1
            field(i, j, k) = temp(n)
            n = n + 1
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DEALLOCATE(temp)

  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, &
    nx_local, ny_local, nz_local)
    USE mpi
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: basetype, sz
    INTEGER(isp):: requests(2)

    basetype = mpidbl

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    ! MOVE EDGES ALONG X
    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    IF (is_dtype_init(4)) THEN
      mpi_dtypes(4) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
      is_dtype_init(4) = .FALSE.
    ENDIF

    ! --- +X
    CALL MPI_ISEND(field(0, -nyg, -nzg), 1_isp, mpi_dtypes(4), INT(proc_x_min, isp),  &
    tag, comm, requests(1), errcode)
    CALL MPI_IRECV(field(nx_local, -nyg, -nzg), 1_isp, mpi_dtypes(4), INT(proc_x_max, &
    isp), tag, comm, requests(2), errcode)

    ! --- Need to wait here to avoid modifying buffer at field(nx_local) and field(0)
    CALL MPI_WAITALL(2_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! --- -X
    CALL MPI_ISEND(field(nx_local-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(4),             &
    INT(proc_x_max, isp), tag, comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(4), INT(proc_x_min,     &
    isp), tag, comm, requests(2), errcode)


    ! NEED TO WAIT BEFORE EXCHANGING ALONG Y (DIAGONAL TERMS)
    CALL MPI_WAITALL(2_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! MOVE EDGES ALONG Y
    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)

    IF (is_dtype_init(5)) THEN
      mpi_dtypes(5) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
      is_dtype_init(5) = .FALSE.
    ENDIF

    ! --- +Y
    CALL MPI_ISEND(field(-nxg, 0, -nzg), 1_isp, mpi_dtypes(5), INT(proc_y_min, isp),  &
    tag, comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg, ny_local, -nzg), 1_isp, mpi_dtypes(5), INT(proc_y_max, &
    isp), tag, comm, requests(2), errcode)

    ! --- Need to wait here to avoid modifying buffer at field(nx_local) and field(0)
    CALL MPI_WAITALL(2_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! --- -Y
    CALL MPI_ISEND(field(-nxg, ny_local-nyg, -nzg), 1_isp, mpi_dtypes(5),             &
    INT(proc_y_max, isp), tag, comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(5), INT(proc_y_min,     &
    isp), tag, comm, requests(2), errcode)

    ! NEED TO WAIT BEFORE EXCHANGING ALONG Z (DIAGONAL TERMS)
    CALL MPI_WAITALL(2_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! MOVE EDGES ALONG Z
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    IF (is_dtype_init(6)) THEN
      mpi_dtypes(6) = create_3d_array_derived_type(basetype, subsizes, sizes, starts)
      is_dtype_init(6) = .FALSE.
    ENDIF

    ! --- +Z
    CALL MPI_ISEND(field(-nxg, -nyg, 0), 1_isp, mpi_dtypes(6), INT(proc_z_min, isp),  &
    tag, comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg, -nyg, nz_local), 1_isp, mpi_dtypes(6), INT(proc_z_max, &
    isp), tag, comm, requests(2), errcode)

    ! --- Need to wait here to avoid modifying buffer at field(nx_local) and field(0)
    CALL MPI_WAITALL(2_isp, requests, MPI_STATUSES_IGNORE, errcode)

    ! --- -Z
    CALL MPI_ISEND(field(-nxg, -nyg, nz_local-nzg), 1_isp, mpi_dtypes(6),             &
    INT(proc_z_max, isp), tag, comm, requests(1), errcode)
    CALL MPI_IRECV(field(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(6), INT(proc_z_min,     &
    isp), tag, comm, requests(2), errcode)

    CALL MPI_WAITALL(2_isp, requests, MPI_STATUSES_IGNORE, errcode)


  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards_nonblocking


  ! ______________________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE summation_bcs(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: array
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: nn, sz

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn =  nx_local

    IF (is_dtype_init(7)) THEN
      mpi_dtypes(7) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      is_dtype_init(7) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn, -nyg, -nzg), 1_isp, mpi_dtypes(7), INT(neighbour( 1,  &
    0, 0), isp), tag, temp, sz, mpidbl, INT(neighbour(-1, 0, 0), isp), tag, comm,     &
    status, errcode)
    array(0:nxg, :, :) = array(0:nxg, :, :) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(7),                  &
    INT(neighbour(-1, 0, 0), isp), tag, temp, sz, mpidbl, INT(neighbour( 1, 0, 0),    &
    isp), tag, comm, status, errcode)
    array(nn-nxg:nn, :, :) = array(nn-nxg:nn, :, :) + temp

    DEALLOCATE(temp)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)
    nn =  ny_local

    IF (is_dtype_init(8)) THEN
      mpi_dtypes(8) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      is_dtype_init(8) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg, nn, -nzg), 1_isp, mpi_dtypes(8), INT(neighbour(0,   &
    1, 0), isp), tag, temp, sz, mpidbl, INT(neighbour(0, -1, 0), isp), tag, comm,     &
    status, errcode)
    array(:, 0:nyg, :) = array(:, 0:nyg, :) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(8), INT(neighbour(0, &
    -1, 0), isp), tag, temp, sz, mpidbl, INT(neighbour(0, 1, 0), isp), tag, comm,     &
    status, errcode)
    array(:, nn-nyg:nn, :) = array(:, nn-nyg:nn, :) + temp

    DEALLOCATE(temp)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
    nn =  nz_local

    IF (is_dtype_init(9)) THEN
      mpi_dtypes(9) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      is_dtype_init(9) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg, -nyg, nn), 1_isp, mpi_dtypes(9), INT(neighbour(0,   &
    0, 1), isp), tag, temp, sz, mpidbl, INT(neighbour(0, 0, -1), isp), tag, comm,     &
    status, errcode)
    array(:, :, 0:nzg) = array(:, :, 0:nzg) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(9), INT(neighbour(0, &
    0, -1), isp), tag, temp, sz, mpidbl, INT(neighbour(0, 0, 1), isp), tag, comm,     &
    status, errcode)
    array(:, :, nn-nzg:nn) = array(:, :, nn-nzg:nn) + temp

    DEALLOCATE(temp)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

  END SUBROUTINE summation_bcs

  ! ______________________________________________________________________________________
  !> Routine for adding current contributions fron adjacent subdomains
  ! nonblocking version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE summation_bcs_nonblocking(array, nxg, nyg, nzg, nx_local, ny_local,      &
    nz_local)
    USE mpi
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: array
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp1, temp2
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: nn, sz
    INTEGER(isp) :: requests(4)

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx_local

    IF (is_dtype_init(10)) THEN
      mpi_dtypes(10) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      is_dtype_init(10) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(nn, -nyg, -nzg), 1_isp, mpi_dtypes(10), INT(proc_x_max,      &
    isp), tag, comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_x_min, isp), tag, comm, requests(2),   &
    errcode)
    CALL MPI_ISEND(array(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(10), INT(proc_x_min,    &
    isp), tag, comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_x_max, isp), tag, comm, requests(4),   &
    errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg, :, :) = array(0:nxg, :, :) + temp1
    array(nn-nxg:nn, :, :) = array(nn-nxg:nn, :, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)
    nn = ny_local

    IF (is_dtype_init(11)) THEN
      mpi_dtypes(11) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      is_dtype_init(11) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(-nxg, nn, -nzg), 1_isp, mpi_dtypes(11), INT(proc_y_max,      &
    isp), tag, comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_y_min, isp), tag, comm, requests(2),   &
    errcode)
    CALL MPI_ISEND(array(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(11), INT(proc_y_min,    &
    isp), tag, comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_y_max, isp), tag, comm, requests(4),   &
    errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(:, 0:nyg, :) = array(:, 0:nyg, :) + temp1
    array(:, nn-nyg:nn, :) = array(:, nn-nyg:nn, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
    nn = nz_local

    IF (is_dtype_init(12)) THEN
      mpi_dtypes(12) = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      is_dtype_init(12) = .FALSE.
    ENDIF

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_ISEND(array(-nxg, -nyg, nn), 1_isp, mpi_dtypes(12), INT(proc_z_max,      &
    isp), tag, comm, requests(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, INT(proc_z_min, isp), tag, comm, requests(2),   &
    errcode)
    CALL MPI_ISEND(array(-nxg, -nyg, -nzg), 1_isp, mpi_dtypes(12), INT(proc_z_min,    &
    isp), tag, comm, requests(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, INT(proc_z_max, isp), tag, comm, requests(4),   &
    errcode)
    CALL MPI_WAITALL(4_isp, requests, MPI_STATUSES_IGNORE, errcode)

    array(:, :, 0:nzg) = array(:, :, 0:nzg) + temp1
    array(:, :, nn-nzg:nn) = array(:, :, nn-nzg:nn) + temp2

    DEALLOCATE(temp1, temp2)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

  END SUBROUTINE summation_bcs_nonblocking

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains persistent version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE summation_bcs_persistent_jx(array, nxg, nyg, nzg, nx_local, ny_local,    &
    nz_local)
    USE communications
    USE mpi
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: array
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp1, temp2
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

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    ! Initialization of the persistent communication
    IF (it .EQ. 0) THEN

      ! Init X
      subsizes(1) = nxg+1
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx_local
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn, -nyg, -nzg), 1_isp, subarrayx, proc_x_max_mpisp,   &
      tag, comm, reqperjxx(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayx, proc_x_min_mpisp, &
      tag, comm, reqperjxx(3), errcode)

      call MPI_RECV_INIT(temp1, sz, mpidbl, proc_x_min_mpisp, tag, comm,              &
      reqperjxx(2), errcode)
      call MPI_RECV_INIT(temp2, sz, mpidbl, proc_x_max_mpisp, tag, comm,              &
      reqperjxx(4), errcode)

      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny_local
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg, nn, -nzg), 1_isp, subarrayy, proc_y_max_mpisp,   &
      tag, comm, reqperjxy(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayy, proc_y_min_mpisp, &
      tag, comm, reqperjxy(3), errcode)

      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz_local
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg, -nyg, nn), 1_isp, subarrayz, proc_z_max_mpisp,   &
      tag, comm, reqperjxz(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayz, proc_z_min_mpisp, &
      tag, comm, reqperjxz(3), errcode)

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)
      CALL MPI_TYPE_FREE(subarrayz, errcode)

    ENDIF

    !! -- Summation along X- direction

    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min_mpisp, tag, comm, reqperjxx(2),      &
    errcode)
    call MPI_START(reqperjxx(1), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max_mpisp, tag, comm, reqperjxx(4),      &
    errcode)
    call MPI_START(reqperjxx(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjxx, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg, :, :) = array(0:nxg, :, :) + temp1
    array(nn-nxg:nn, :, :) = array(nn-nxg:nn, :, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Y- direction

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)
    nn = ny_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    call MPI_START(reqperjxy(1), errcode)
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min_mpisp, tag, comm, reqperjxy(2),      &
    errcode)
    call MPI_START(reqperjxy(3), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max_mpisp, tag, comm, reqperjxy(4),      &
    errcode)
    CALL MPI_WAITALL(4_isp, reqperjxy, MPI_STATUSES_IGNORE, errcode)

    array(:, 0:nyg, :) = array(:, 0:nyg, :) + temp1
    array(:, nn-nyg:nn, :) = array(:, nn-nyg:nn, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
    nn = nz_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min_mpisp, tag, comm, reqperjxz(2),      &
    errcode)
    call MPI_START(reqperjxz(1), errcode)

    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max_mpisp, tag, comm, reqperjxz(4),      &
    errcode)
    call MPI_START(reqperjxz(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjxz, MPI_STATUSES_IGNORE, errcode)

    array(:, :, 0:nzg) = array(:, :, 0:nzg) + temp1
    array(:, :, nn-nzg:nn) = array(:, :, nn-nzg:nn) + temp2

    DEALLOCATE(temp1, temp2)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
  END SUBROUTINE

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains persistent version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE summation_bcs_persistent_jy(array, nxg, nyg, nzg, nx_local, ny_local,    &
    nz_local)
    USE communications
    USE mpi
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: array
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp1, temp2
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

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    ! Initialization of the persistent communication
    IF (it .EQ. 0) THEN

      ! Init X
      subsizes(1) = nxg+1
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx_local
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn, -nyg, -nzg), 1_isp, subarrayx, proc_x_max_mpisp,   &
      tag, comm, reqperjyx(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayx, proc_x_min_mpisp, &
      tag, comm, reqperjyx(3), errcode)

      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny_local
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg, nn, -nzg), 1_isp, subarrayy, proc_y_max_mpisp,   &
      tag, comm, reqperjyy(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayy, proc_y_min_mpisp, &
      tag, comm, reqperjyy(3), errcode)

      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz_local
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg, -nyg, nn), 1_isp, subarrayz, proc_z_max_mpisp,   &
      tag, comm, reqperjyz(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayz, proc_z_min_mpisp, &
      tag, comm, reqperjyz(3), errcode)

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)
      CALL MPI_TYPE_FREE(subarrayz, errcode)

    ENDIF

    !! -- Summation along X- direction

    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min_mpisp, tag, comm, reqperjyx(2),      &
    errcode)
    call MPI_START(reqperjyx(1), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max_mpisp, tag, comm, reqperjyx(4),      &
    errcode)
    call MPI_START(reqperjyx(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjyx, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg, :, :) = array(0:nxg, :, :) + temp1
    array(nn-nxg:nn, :, :) = array(nn-nxg:nn, :, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Y- direction

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)
    nn = ny_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min_mpisp, tag, comm, reqperjyy(2),      &
    errcode)
    CALL MPI_START(reqperjyy(1), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max_mpisp, tag, comm, reqperjyy(4),      &
    errcode)
    call MPI_START(reqperjyy(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjyy, MPI_STATUSES_IGNORE, errcode)

    array(:, 0:nyg, :) = array(:, 0:nyg, :) + temp1
    array(:, nn-nyg:nn, :) = array(:, nn-nyg:nn, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
    nn = nz_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))


    temp1  = 0.0_num
    temp2 = 0.0_num

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min_mpisp, tag, comm, reqperjyz(2),      &
    errcode)
    call MPI_START(reqperjyz(1), errcode)

    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max_mpisp, tag, comm, reqperjyz(4),      &
    errcode)
    call MPI_START(reqperjyz(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjyz, MPI_STATUSES_IGNORE, errcode)

    array(:, :, 0:nzg) = array(:, :, 0:nzg) + temp1
    array(:, :, nn-nzg:nn) = array(:, :, nn-nzg:nn) + temp2

    DEALLOCATE(temp1, temp2)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
  END SUBROUTINE

  ! ______________________________________________________________________________________
  !> @brief
  !> Routine for adding current contributions fron adjacent subdomains persistent version
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE summation_bcs_persistent_jz(array, nxg, nyg, nzg, nx_local, ny_local,    &
    nz_local)
    USE communications
    USE mpi
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg, -nyg:ny_local+nyg, -nzg:nz_local+nzg),    &
    INTENT(INOUT) :: array
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp1, temp2
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

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    ! Initialization of the persistent communication
    IF (it .EQ. 0) THEN

      ! Init X
      subsizes(1) = nxg+1
      subsizes(2) = sizes(2)
      subsizes(3) = sizes(3)
      nn = nx_local
      subarrayx = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(nn, -nyg, -nzg), 1_isp, subarrayx, proc_x_max_mpisp,   &
      tag, comm, reqperjzx(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayx, proc_x_min_mpisp, &
      tag, comm, reqperjzx(3), errcode)

      ! Init Y
      subsizes(1) = sizes(1)
      subsizes(2) = nyg+1
      subsizes(3) = sizes(3)
      nn = ny_local
      subarrayy = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg, nn, -nzg), 1_isp, subarrayy, proc_y_max_mpisp,   &
      tag, comm, reqperjzy(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayy, proc_y_min_mpisp, &
      tag, comm, reqperjzy(3), errcode)

      ! Init Z
      subsizes(1) = sizes(1)
      subsizes(2) = sizes(2)
      subsizes(3) = nzg+1
      nn = nz_local
      subarrayz = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)
      call MPI_SEND_INIT(array(-nxg, -nyg, nn), 1_isp, subarrayz, proc_z_max_mpisp,   &
      tag, comm, reqperjzz(1), errcode)
      call MPI_SEND_INIT(array(-nxg, -nyg, -nzg), 1_isp, subarrayz, proc_z_min_mpisp, &
      tag, comm, reqperjzz(3), errcode)

      CALL MPI_TYPE_FREE(subarrayx, errcode)
      CALL MPI_TYPE_FREE(subarrayy, errcode)
      CALL MPI_TYPE_FREE(subarrayz, errcode)

    ENDIF

    !! -- Summation along X- direction

    subsizes(1) = nxg+1
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_x_min_mpisp, tag, comm, reqperjzx(2),      &
    errcode)
    call MPI_START(reqperjzx(1), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_x_max_mpisp, tag, comm, reqperjzx(4),      &
    errcode)
    call MPI_START(reqperjzx(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjzx, MPI_STATUSES_IGNORE, errcode)

    array(0:nxg, :, :) = array(0:nxg, :, :) + temp1
    array(nn-nxg:nn, :, :) = array(nn-nxg:nn, :, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Y- direction

    subsizes(1) = sizes(1)
    subsizes(2) = nyg+1
    subsizes(3) = sizes(3)
    nn = ny_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num
    CALL MPI_IRECV(temp1, sz, mpidbl, proc_y_min_mpisp, tag, comm, reqperjzy(2),      &
    errcode)
    CALL MPI_START(reqperjzy(1), errcode)
    CALL MPI_IRECV(temp2, sz, mpidbl, proc_y_max_mpisp, tag, comm, reqperjzy(4),      &
    errcode)
    call MPI_START(reqperjzy(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjzy, MPI_STATUSES_IGNORE, errcode)

    array(:, 0:nyg, :) = array(:, 0:nyg, :) + temp1
    array(:, nn-nyg:nn, :) = array(:, nn-nyg:nn, :) + temp2

    DEALLOCATE(temp1, temp2)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg+1
    nn = nz_local

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1),         &
    subsizes(2), subsizes(3)))

    temp1  = 0.0_num
    temp2 = 0.0_num

    CALL MPI_IRECV(temp1, sz, mpidbl, proc_z_min_mpisp, tag, comm, reqperjzz(2),      &
    errcode)
    call MPI_START(reqperjzz(1), errcode)

    CALL MPI_IRECV(temp2, sz, mpidbl, proc_z_max_mpisp, tag, comm, reqperjzz(4),      &
    errcode)
    call MPI_START(reqperjzz(3), errcode)
    CALL MPI_WAITALL(4_isp, reqperjzz, MPI_STATUSES_IGNORE, errcode)

    array(:, :, 0:nzg) = array(:, :, 0:nzg) + temp1
    array(:, :, nn-nzg:nn) = array(:, :, nn-nzg:nn) + temp2

    DEALLOCATE(temp1, temp2)

    CALL field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
  END SUBROUTINE

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition routine for electric field
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE efield_bcs
    REAL(num) :: tmptime
#if defined(DEBUG)
    WRITE(0, *) "efield_bcs: start"
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
    WRITE(0, *) "efield_bcs: stop"
#endif
  END SUBROUTINE efield_bcs

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary condition routine for magnetic field
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE bfield_bcs

    REAL(num) :: tmptime
#if defined(DEBUG)
    WRITE(0, *) "bfield_bcs: start"
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
    WRITE(0, *) "bfield_bcs: stop"
#endif
  END SUBROUTINE bfield_bcs

  ! ______________________________________________________________________________________
  !> @brief
  !> Boundary conditions routine for currents
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE current_bcs
    REAL(num) :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Add current contribution from adjacent subdomains

    IF (mpicom_curr.EQ.2) THEN

      CALL summation_bcs_persistent_jx(jx, nxjguards, nyjguards, nzjguards, nx, ny,   &
      nz)
      CALL summation_bcs_persistent_jy(jy, nxjguards, nyjguards, nzjguards, nx, ny,   &
      nz)
      CALL summation_bcs_persistent_jz(jz, nxjguards, nyjguards, nzjguards, nx, ny,   &
      nz)

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

  ! ______________________________________________________________________________________
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
  ! ______________________________________________________________________________________
  SUBROUTINE charge_bcs
    USE time_stat

    REAL(num) :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    IF (mpicom_curr.EQ.1) THEN
      CALL summation_bcs(rho, nxjguards, nyjguards, nzjguards, nx, ny, nz)
    ELSE
      CALL summation_bcs_nonblocking(rho, nxjguards, nyjguards, nzjguards, nx, ny,    &
      nz)
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(13) = localtimes(13) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE charge_bcs

END MODULE field_boundary

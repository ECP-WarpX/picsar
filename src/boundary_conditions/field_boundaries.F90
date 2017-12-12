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
  !> <@brief
  !> Boundary condition routine for E and B between mpi_groups
  !
  !> @author
  !> Haithem Kallala
  !>date
  !> Creation 2017
  !
  !
  ! ______________________________________________________________________________________
  SUBROUTINE ebj_field_bcs_groups()
#if defined(FFTW)
    USE group_parameters
    USE mpi_fftw3
#endif
    USE shared_data
    REAL(num) :: tmptime
    INTEGER(idp)     ::  size_nx
#if defined(DEBUG)
    WRITE(0, *) "efield_bcs_group: start"
#endif
#if defined(FFTW)
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    size_nx = 2*(nx_group/2+1)
    IF(mpicom_curr .EQ. 0) THEN
      CALL field_bc_group_non_blocking(ex_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_non_blocking(ey_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_non_blocking(ez_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_non_blocking(bx_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_non_blocking(by_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_non_blocking(bz_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_non_blocking(jx_r, size_nx, ny_group, local_nz,           &
      nzg_group)
      CALL field_bc_group_non_blocking(jy_r, size_nx, ny_group, local_nz,           &
      nzg_group)
      CALL field_bc_group_non_blocking(jz_r, size_nx, ny_group, local_nz,           &
      nzg_group)
      CALL field_bc_group_non_blocking(rho_r, size_nx, ny_group, local_nz,          &
      nzg_group)
      CALL field_bc_group_non_blocking(rhoold_r, size_nx, ny_group, local_nz,       &
      nzg_group)
    ELSE
      CALL field_bc_group_blocking(ex_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(ey_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(ez_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(bx_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(by_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(bz_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(jx_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(jy_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(jz_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(rho_r, size_nx, ny_group, local_nz, nzg_group)
      CALL field_bc_group_blocking(rhoold_r, size_nx, ny_group, local_nz,           &
      nzg_group)

    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
    ENDIF
#endif
  END SUBROUTINE ebj_field_bcs_groups


!This subroutine is still not working correctly so please use mpicom_curr=0 with
!groups 
  SUBROUTINE field_bc_group_blocking(field, nxx, nyy, nzz, ngroupz)
#if defined(FFTW)
    USE group_parameters
#endif
    USE shared_data
    INTEGER(idp), INTENT(IN)  :: nxx, nyy, nzz, ngroupz
    REAL(num), INTENT(INOUT), DIMENSION(1:nxx, 1:nyy, 1:nzz)  :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: basetype
    REAL(num)   , ALLOCATABLE , DIMENSION(:) :: temp
    INTEGER(idp)                             :: sz,ix,iy,iz,n
#if defined(FFTW)
    basetype = mpidbl
    sizes(1) = nxx
    sizes(2) = nyy
    sizes(3) = nzz
    starts=1
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ngroupz

    IF (is_dtype_init(20)) THEN
      mpi_dtypes(20) = create_3d_array_derived_type(basetype, subsizes, sizes,        &
      starts)
      is_dtype_init(20) = .FALSE.
    ENDIF

    IF(XOR(group_z_min_boundary,group_z_max_boundary)) THEN
      IF(group_z_min_boundary) THEN
        CALL MPI_SEND(field(1, 1, iz_min_r), 1_isp, mpi_dtypes(20), INT(proc_z_min,   &
        isp), tag, comm, errcode)
        CALL MPI_RECV(field(1, 1, 1), 1_isp, mpi_dtypes(20), INT(proc_z_min,isp), tag,&
        comm,MPI_STATUS_IGNORE, errcode)
      ENDIF
      IF(group_z_max_boundary) THEN
        CALL MPI_RECV(field(1, 1, iz_max_r+1), 1_isp, mpi_dtypes(20), INT(proc_z_max, &
        isp), tag, comm,MPI_STATUS_IGNORE, errcode)
        CALL MPI_SEND(field(1, 1, iz_max_r-ngroupz +1), 1_isp, mpi_dtypes(20),        &
        INT(proc_z_max, isp), tag, comm, errcode)
      ENDIF
    ELSE 
      sz  = subsizes(1)*subsizes(2)*subsizes(3)
      ALLOCATE(temp(sz))
      CALL MPI_SENDRECV(field(1,1,iz_min_r), 1_isp, mpi_dtypes(20),INT(proc_z_min,    &
      isp), tag, temp, sz, basetype, INT(proc_z_max, isp), tag, comm, status,errcode)
      IF(proc_z_max .NE. MPI_PROC_NULL) THEN
        n=1
        DO iz= iz_max_r+1,nzz
          DO iy=1,nyy
            DO ix=1,nxx
               field(ix,iy,iz) = temp(n)
               n = n+1
            ENDDO 
          ENDDO
        ENDDO
      ENDIF
      CALL MPI_SENDRECV(field(1,1,iz_max_r-ngroupz+1), 1_isp,mpi_dtypes(20),          &
      INT(proc_z_max,isp), tag, temp, sz, basetype, INT(proc_z_min, isp), tag, comm,  &
      status,errcode)
      IF(proc_z_min .NE. MPI_PROC_NULL) THEN
        n=1
        DO iz= 1,ngroupz
          DO iy=1,nyy
            DO ix=1,nxx
               field(ix,iy,iz) = temp(n)
               n = n+1
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(temp)
    ENDIF
#endif
  END SUBROUTINE field_bc_group_blocking



  SUBROUTINE field_bc_group_non_blocking(field, nxx, nyy, nzz, ngroupz)
#if defined(FFTW)
    USE group_parameters
#endif
    USE shared_data
    INTEGER(idp), INTENT(IN)  :: nxx, nyy, nzz, ngroupz
    REAL(num), INTENT(INOUT), DIMENSION(1:nxx, 1:nyy, 1:nzz)  :: field
    INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER(isp) :: basetype
    INTEGER(isp):: requests_1(2), requests_2(2)

    basetype = mpidbl
    sizes(1) = nxx
    sizes(2) = nyy
    sizes(3) = nzz
    starts=1
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ngroupz
    IF (is_dtype_init(20)) THEN
      mpi_dtypes(20) = create_3d_array_derived_type(basetype, subsizes, sizes,        &
      starts)
      is_dtype_init(20) = .FALSE.
    ENDIF
#if defined(FFTW)
!CASE Where each group has more than one mpi task
    IF (XOR(group_z_min_boundary,group_z_max_boundary))  THEN
      IF(group_z_min_boundary) THEN
        CALL MPI_IRECV(field(1, 1, 1), 1_isp,mpi_dtypes(20),INT(proc_z_min,isp),      &
        tag, comm, requests_1(1), errcode)
        CALL MPI_ISEND(field(1, 1, iz_min_r), 1_isp, mpi_dtypes(20),INT(proc_z_min,   &
        isp), tag, comm, requests_1(2), errcode)
        CALL MPI_WAITALL(2_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
      ENDIF
      IF(group_z_max_boundary) THEN
        CALL MPI_IRECV(field(1, 1, iz_max_r+1), 1_isp, mpi_dtypes(20),INT(proc_z_max, &
        isp), tag, comm, requests_2(1), errcode)
        CALL MPI_ISEND(field(1, 1, iz_max_r-ngroupz +1), 1_isp, mpi_dtypes(20),       &
        INT(proc_z_max, isp), tag, comm, requests_2(2), errcode)
        CALL MPI_WAITALL(2_isp, requests_2, MPI_STATUSES_IGNORE, errcode)
      ENDIF
    ENDIF
!case where each group has one mpi task 
    IF(group_z_min_boundary .AND. group_z_max_boundary) THEN
      CALL MPI_IRECV(field(1, 1, iz_max_r+1), 1_isp,mpi_dtypes(20),INT(proc_z_max,   &
      isp), tag, comm, requests_1(2), errcode)
      CALL MPI_ISEND(field(1, 1, iz_min_r), 1_isp, mpi_dtypes(20),INT(proc_z_min,     &
      isp), tag, comm, requests_1(1), errcode)
      CALL MPI_WAITALL(2_isp, requests_1, MPI_STATUSES_IGNORE, errcode)

      CALL MPI_IRECV(field(1, 1, 1), 1_isp, mpi_dtypes(20), INT(proc_z_min,isp),      &
      tag, comm, requests_2(2), errcode)
      CALL MPI_ISEND(field(1, 1, iz_max_r-ngroupz +1), 1_isp, mpi_dtypes(20),         &
      INT(proc_z_max, isp), tag, comm, requests_2(1), errcode)
      CALL MPI_WAITALL(2_isp, requests_2, MPI_STATUSES_IGNORE, errcode)
    ENDIF
#endif
  END SUBROUTINE field_bc_group_non_blocking
  ! ______________________________________________________________________________________
  !> Routine for backward communications between groups
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017  
  !
  ! ______________________________________________________________________________________
!  SUBROUTINE group_communication_backward(coeff_norm)
!#if defined(FFTW)
!  USE group_parameters
!  USE mpi_fftw3
!#endif
!  USE shared_data
!  USE mpi  
!  USE omp_lib
!  USE time_stat
!  REAL(num)  ,INTENT(IN)  :: coeff_norm
!  INTEGER(idp) , DIMENSION(c_ndims) :: sizes, subsizes_left, subsizes_right , starts
!  INTEGER(isp) :: basetype
!  INTEGER(isp) :: requests_1(1)
!  INTEGER(idp)  :: nxx,nyy,nzz
!  REAL(num)     :: tmptime
!
!#if defined(FFTW)
!  IF (it.ge.timestat_itstart) THEN
!    tmptime = MPI_WTIME()
!  ENDIF
!  basetype = mpidbl
!  nxx = 2*(nx_group/2+1)
!  nyy = ny_group
!  nzz = local_nz
!
!  sizes(1) =nxx 
!  sizes(2) =nyy 
!  sizes(3) =nzz 
!  starts=1
!
!  subsizes_right(1) = sizes(1)
!  subsizes_right(2) = sizes(2) 
!  subsizes_right(3) = size_right
!
!  subsizes_left(1) = sizes(1) 
!  subsizes_left(2) = sizes(2)  
!  subsizes_left(3) = rsize_left
!
!  IF (is_dtype_init(40)) THEN
!    mpi_dtypes(40) = create_3d_array_derived_type(basetype, subsizes_right,sizes,  &
!    starts)
!    is_dtype_init(40) = .FALSE.
!  ENDIF
!
!  IF (is_dtype_init(41)) THEN
!    mpi_dtypes(41) = create_3d_array_derived_type(basetype,subsizes_left,sizes,  &
!    starts)
!    is_dtype_init(41) = .FALSE.
!  ENDIF
!
!  subsizes_left(3) = size_left
!  subsizes_right(3) = rsize_right
!
!  IF (is_dtype_init(42)) THEN
!
!    mpi_dtypes(42) = create_3d_array_derived_type(basetype,subsizes_left,sizes,&
!    starts)
!    is_dtype_init(42) = .FALSE.
!  ENDIF
!  IF (is_dtype_init(43)) THEN
!    mpi_dtypes(43) = create_3d_array_derived_type(basetype,subsizes_right,sizes,&
!    starts)
!    is_dtype_init(43) = .FALSE.
!  ENDIF
!
!
!  !Send ex_r to ex  of proc_z_max
!  CALL SEND_TO_RIGHT_f2r(ex_r,ex,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_RIGHT_f2r(ey_r,ey,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_RIGHT_f2r(ez_r,ez,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_RIGHT_f2r(bx_r,bx,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_RIGHT_f2r(by_r,by,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_RIGHT_f2r(bz_r,bz,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!
!  
!  !> Send ex_r to ex of proc_z_min
!  CALL SEND_TO_LEFT_f2r(ex_r,ex,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_LEFT_f2r(ey_r,ey,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_LEFT_f2r(ez_r,ez,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_LEFT_f2r(bx_r,bx,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_LEFT_f2r(by_r,by,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!  CALL SEND_TO_LEFT_f2r(bz_r,bz,nxx,nyy,nzz,coeff_norm,nx,ny,nz,nxguards,nyguards,nzguards)
!
!  IF (it.ge.timestat_itstart) THEN
!    localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
!  ENDIF
!#endif
!  END SUBROUTINE group_communication_backward
!
!  SUBROUTINE SEND_TO_RIGHT_f2r(field_in,field_out,nxx,nyy,nzz,coeff_norm,nx1,ny1,nz1,nxg,nyg,nzg)
!#if defined(FFTW)
!   USE group_parameters
!#endif
!   USE shared_data
!   USE mpi
!   REAL(num) , INTENT(IN) :: coeff_norm
!   INTEGER(idp) , INTENT(IN)  :: nxx,nyy,nzz,nxg,nyg,nzg,nx1,ny1,nz1
!   REAL(num) , INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg) ::field_out
!   REAL(num) , INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz) :: field_in
!   INTEGER(isp)   :: requests_1(1)
!   INTEGER(idp)   :: ix,iy,iz
!   REAL(num), ALLOCATABLE, DIMENSION(:,:,:) ::  temp_from_right
!
!#if defined(FFTW)
!    IF(z_coords .NE. nprocz-1 .AND. size_right .NE. 0_isp) THEN
!      CALL MPI_ISEND(field_in(1,1,g_right(1)),1_isp,mpi_dtypes(40),&
!      INT(proc_z_max,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!    ENDIF
!
!    IF(z_coords .NE. 0 .AND. rsize_left .NE. 0_idp) THEN
!      ALLOCATE(temp_from_right(nxx,nyy,rsize_left))
!
!      CALL MPI_IRECV(temp_from_right,1_isp,mpi_dtypes(41),Int(proc_z_min,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!
!      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
!      DO iy=iy_min_r, iy_max_r
!        DO ix=ix_min_r, ix_max_r
!          DO iz=1,rsize_left
!             field_out(ix-ix_min_r-nxg,iy-iy_min_r-nyg,rr_left(iz)) =&
!                 coeff_norm*temp_from_right(ix-ix_min_r+1,iy-iy_min_r+1,iz)
!          ENDDO
!        ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO
!
!      DEALLOCATE(temp_from_right)
!    ENDIF
!#endif
!  END SUBROUTINE SEND_TO_RIGHT_f2r
!
!
!  SUBROUTINE SEND_TO_LEFT_f2r(field_in,field_out,nxx,nyy,nzz,coeff_norm,nx1,ny1,nz1,nxg,nyg,nzg)
!#if defined(FFTW)
!   USE group_parameters
!#endif
!   USE shared_data
!   USE mpi
!   REAL(num) , INTENT(IN) :: coeff_norm
!   INTEGER(idp) , INTENT(IN)  :: nxx,nyy,nzz,nx1,ny1,nz1,nxg,nyg,nzg
!   REAL(num) , INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg) ::field_out
!   REAL(num) , INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz) :: field_in
!   INTEGER(isp)   :: requests_1(1)
!   INTEGER(idp)   :: ix,iy,iz
!   REAL(num), ALLOCATABLE, DIMENSION(:,:,:) ::  temp_from_left
!
!#if defined(FFTW)
!    IF(z_coords .NE. 0 .AND. size_left .NE. 0_isp) THEN
!      CALL MPI_ISEND(field_in(1,1,g_left(1)),1_isp,mpi_dtypes(42),&
!      INT(proc_z_min,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!    ENDIF
!
!    IF(z_coords .NE. nprocz-1 .AND. rsize_right .NE. 0_idp) THEN
!      ALLOCATE(temp_from_left(nxx,nyy,rsize_right))
!
!      CALL MPI_IRECV(temp_from_left,1_isp,mpi_dtypes(43),Int(proc_z_max,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!
!      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
!      DO iy=iy_min_r, iy_max_r
!        DO ix=ix_min_r, ix_max_r
!          DO iz=1,rsize_right
!             field_out(ix-ix_min_r-nxg,iy-iy_min_r-nyg,rr_right(iz)) =&
!                 coeff_norm*temp_from_left(ix-ix_min_r+1,iy-iy_min_r+1,iz)
!          ENDDO
!        ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO
!
!      DEALLOCATE(temp_from_left)
!    ENDIF
!#endif
!  END SUBROUTINE SEND_TO_LEFT_f2r
!
!
!  ! ______________________________________________________________________________________
!  !> Routine for forward comms between groups
!  !
!  !> @author
!  !> Haithem Kallala
!  !
!  !> @date
!  !> Creation 2015
!  !
!  ! ______________________________________________________________________________________
!  SUBROUTINE group_communication_forward
!#if defined(FFTW)
!    USE group_parameters
!    USE mpi_fftw3
!#endif
!    USE shared_data
!    USE mpi 
!    USE omp_lib
!    USE time_stat
!    REAL(num) , ALLOCATABLE, DIMENSION(:,:,:) :: temp_from_right
!    INTEGER(idp) , DIMENSION(c_ndims) :: sizes, subsizes_left, subsizes_right , starts
!    INTEGER(isp) :: basetype
!    INTEGER(isp) :: requests_1(1)
!    INTEGER(idp)  :: nxx,nyy,nzz
!    REAL(num)     :: tmptime
!
!#if defined(FFTW)
!    IF (it.ge.timestat_itstart) THEN
!      tmptime = MPI_WTIME()
!    ENDIF
!    basetype = mpidbl
!    sizes(1) = nx+2*nxguards+1
!    sizes(2) = ny+2*nyguards+1
!    sizes(3) = nz+2*nzguards+1
!
!    nxx = 2*(nx_group/2+1)
!    nyy = ny_group
!    nzz = local_nz
!    starts=1
!
!    subsizes_left(1) = sizes(1)
!    subsizes_left(2) = sizes(2)
!    subsizes_left(3) = size_left
!
!    subsizes_right(1) = sizes(1)
!    subsizes_right(2) = sizes(2)
!    subsizes_right(3) = rsize_right
!
!    IF (is_dtype_init(30)) THEN
!      mpi_dtypes(30) = create_3d_array_derived_type(basetype, subsizes_right, sizes,  &
!      starts)
!      is_dtype_init(30) = .FALSE.
!    ENDIF
!
!    IF (is_dtype_init(31)) THEN
!      mpi_dtypes(31) = create_3d_array_derived_type(basetype, subsizes_left, sizes,   &
!      starts)
!      is_dtype_init(31) = .FALSE.
!    ENDIF
!
!    subsizes_right(3) = size_right
!    subsizes_left(3) = rsize_left
!    IF (is_dtype_init(32)) THEN
!      mpi_dtypes(32) = create_3d_array_derived_type(basetype,subsizes_left,sizes,  &
!      starts)
!      is_dtype_init(32) = .FALSE.
!    ENDIF
!    IF (is_dtype_init(33)) THEN
!      mpi_dtypes(33) = create_3d_array_derived_type(basetype,subsizes_right,sizes,   &
!      starts)
!      is_dtype_init(33) = .FALSE.
!    ENDIF
!
!    !SEND ex to  ex_r  at(proc_z_max)
!    CALL SEND_TO_RIGHT_r2f(ex,ex_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(ey,ey_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(ez,ez_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(bx,bx_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(by,by_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(bz,bz_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(jx,jx_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(jy,jy_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(jz,jz_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(rho,rhoold_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_RIGHT_r2f(rhoold,rhoold_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
! 
!    !> Send ex to ex_r at (proc_z_min)
!    CALL SEND_TO_LEFT_r2f(ex,ex_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(ey,ey_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(ez,ez_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(bx,bx_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(by,by_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(bz,bz_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(jx,jx_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(jy,jy_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(jz,jz_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(rho,rhoold_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!    CALL SEND_TO_LEFT_r2f(rhoold,rhoold_r,nxx,nyy,nzz,nx,ny,nz,nxguards,nyguards,nzguards)
!
!    IF (it.ge.timestat_itstart) THEN
!      localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
!    ENDIF
!#endif
!  END SUBROUTINE group_communication_forward
!   
!  SUBROUTINE SEND_TO_RIGHT_r2f(field_in,field_out,nxx,nyy,nzz,nx1,ny1,nz1,nxg,nyg,nzg)
!#if defined(FFTW)
!   USE group_parameters
!#endif
!   USE shared_data
!   USE mpi
!   INTEGER(idp) , INTENT(IN)  :: nxx,nyy,nzz ,nx1,ny1,nz1,nxg,nyg,nzg
!   REAL(num) , INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg) :: field_in
!   REAL(num) , INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz) :: field_out
!   INTEGER(isp)   :: requests_1(1)
!   INTEGER(idp)   :: ix,iy,iz   
!   REAL(num), ALLOCATABLE, DIMENSION(:,:,:) ::  temp_from_right
!
!#if defined(FFTW)
!    IF(z_coords .NE. nprocz-1 .AND. rsize_right .NE. 0_isp) THEN
!      CALL MPI_ISEND(field_in(-nxg,-nyg,rr_right(1)),1_isp,mpi_dtypes(30),           &
!      INT(proc_z_max,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!    ENDIF
!
!    IF(z_coords .NE. 0 .AND. size_left .NE. 0_idp) THEN
!      ALLOCATE(temp_from_right(nx1+2*nxg+1,ny1+2*nyg+1,size_left))
!      CALL MPI_IRECV(temp_from_right(1,1,1),1_isp,mpi_dtypes(31),Int(proc_z_min,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!
!      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
!      DO iz=1,size_left
!        DO iy=iy_min_r, iy_max_r
!          DO ix=ix_min_r, ix_max_r
!           field_out(ix,iy,g_left(iz)) = temp_from_right(ix-ix_min_r+1,iy-iy_min_r+1,iz)
!          ENDDO
!        ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO
!
!      DEALLOCATE(temp_from_right)
!    ENDIF
!#endif
!  END SUBROUTINE SEND_TO_RIGHT_r2f 
!
!
!  SUBROUTINE SEND_TO_LEFT_r2f(field_in,field_out,nxx,nyy,nzz,nx1,ny1,nz1,nxg,nyg,nzg)
!#if defined(FFTW)
!   USE group_parameters
!#endif
!   USE shared_data
!   USE mpi
!   INTEGER(idp) , INTENT(IN)  :: nxx,nyy,nzz ,nx1,ny1,nz1,nxg,nyg,nzg
!   REAL(num) , INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg) :: field_in
!   REAL(num) , INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz) :: field_out
!   INTEGER(isp)   :: requests_1(1)
!   INTEGER(idp)   :: ix,iy,iz   
!   REAL(num), ALLOCATABLE, DIMENSION(:,:,:) ::  temp_from_left
!
!#if defined(FFTW)
!    IF(z_coords .NE. 0 .AND. rsize_left .NE. 0_isp) THEN
!      CALL MPI_ISEND(field_in(-nxguards,-nyguards,rr_left(1)),1_isp,mpi_dtypes(32),           &
!      INT(proc_z_min,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!    ENDIF
!
!    IF(z_coords .NE. nprocz-1 .AND. size_right .NE. 0_idp) THEN
!
!      ALLOCATE(temp_from_left(nx1+2*nxg+1,ny1+2*nyg+1,size_right))
!      CALL MPI_IRECV(temp_from_left,1_isp,mpi_dtypes(33),Int(proc_z_max,isp),tag,comm,requests_1(1),errcode)
!      CALL MPI_WAITALL(1_isp, requests_1, MPI_STATUSES_IGNORE, errcode)
!
!      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
!      DO iz=1,size_right
!        DO iy=iy_min_r, iy_max_r
!          DO ix=ix_min_r, ix_max_r
!             field_out(ix,iy,g_right(iz)) = temp_from_left(ix-ix_min_r+1,iy-iy_min_r+1,iz)
!          ENDDO
!        ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO
!
!      DEALLOCATE(temp_from_left)
!    ENDIF
!#endif
!  END SUBROUTINE SEND_TO_LEFT_r2f 
!



  SUBROUTINE generalized_comms_group_r2f()
#if defined(FFTW) 
    USE group_parameters
    USE mpi_fftw3
#endif
    USE time_stat
    USE shared_data
    USE load_balance  
    USE fields
#if defined(FFTW)
    INTEGER(idp)  :: nxx, nyy, nzz
    REAL(num)                                   :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    nxx = 2*(nx_group/2+1)
    nyy = ny_group
    nzz = local_nz
    CALL sendrecv_rf_generalized(ex,nx,nxguards,ny,nyguards,nz,nzguards,ex_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(ey,nx,nxguards,ny,nyguards,nz,nzguards,ey_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(ez,nx,nxguards,ny,nyguards,nz,nzguards,ez_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(bx,nx,nxguards,ny,nyguards,nz,nzguards,bx_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(by,nx,nxguards,ny,nyguards,nz,nzguards,by_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(bz,nx,nxguards,ny,nyguards,nz,nzguards,bz_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(jx,nx,nxguards,ny,nyguards,nz,nzguards,jx_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(jy,nx,nxguards,ny,nyguards,nz,nzguards,jy_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(jz,nx,nxguards,ny,nyguards,nz,nzguards,jz_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(rho,nx,nxguards,ny,nyguards,nz,nzguards,rho_r,nxx,nyy,nzz)
    CALL sendrecv_rf_generalized(rhoold,nx,nxguards,ny,nyguards,nz,nzguards,rhoold_r,nxx,nyy,nzz)
    IF (it.ge.timestat_itstart) THEN
      localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
    ENDIF
#endif

  END SUBROUTINE generalized_comms_group_r2f
  
  SUBROUTINE sendrecv_rf_generalized(field,nx1,nxg,ny1,nyg,nz1,nzg,field_f,nxx,nyy,nzz)

    USE load_balance
#if defined(FFTW)
    USE group_parameters
#endif
    INTEGER(idp), INTENT(IN)                    ::  nx1,nxg,ny1,nyg,nz1,nzg,nxx,nyy,nzz
    REAL(num)    ,INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg)  :: field
    REAL(num)    ,INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz)  :: field_f
    INTEGER(idp)                                ::  i , j , k , ix , iy , iz
    REAL(num) , ALLOCATABLE , DIMENSION(:,:,:)  :: temp_recv
    INTEGER(isp)                                :: requests(2)
    INTEGER(isp)                                :: rank_to_send_to, rank_to_recv_from
    INTEGER(idp)  :: n
#if defined(FFTW)
    !  i-1 = mpi test in z direction for exchanges
    ! example :
    ! i-1=0 exchange between current rank and itself
    ! i-1=1 send to right (proc_z_max)
    ! i-1=nprocz-1 send to left(proc_z_min)
    DO i=1,nprocz
        n=0
      IF(i == 1 ) CYCLE
      j = MODULO(z_coords+i-1,nprocz) +1
      ! j corresponds to the z_coords(+1) of mpi task to which the send is done 
      ! k corresponds to the z_coords(+1) of mpi task from which the recv is
      ! done
      k = MODULO(z_coords-(i-1),nprocz) +1

      rank_to_send_to = INT(array_of_ranks_to_send_to(i),isp)
      rank_to_recv_from = INT(array_of_ranks_to_recv_from(i),isp)
      IF(sizes_to_exchange_f_to_recv(k) == 0) rank_to_recv_from = MPI_PROC_NULL
      IF(sizes_to_exchange_r_to_send(j) == 0) rank_to_send_to = MPI_PROC_NULL

      ALLOCATE(temp_recv(-nxg:nxg+nx1,-nyg:nyg+ny1,1:sizes_to_exchange_f_to_recv(k)))
!      IF(sizes_to_exchange_r_to_send(j) .GT. 0) THEN
!        CALL MPI_SEND(field(1,1,r_first_cell_to_send(j)),1_isp,send_type_r(j), &
!             rank_to_send_to,tag,comm,errcode)
!                endif
!        if(sizes_to_exchange_f_to_recv(k) .GT. 0) then
!         call mpi_recv(temp_recv(1,1,1),1_isp,recv_type_f(k),rank_to_recv_from,tag,comm,status,errcode)
!        endif

      CALL MPI_SENDRECV(field(-nxg,-nyg,r_first_cell_to_send(j)),1_isp,send_type_r(j),&
      rank_to_send_to,tag,temp_recv(-nxg,-nyg,1),1_isp,recv_type_f(k),rank_to_recv_from,tag,comm,status,errcode)
    !  ALLOCATE(temp_recv(2*nxg+nx1+1,2*nyg+ny1+1,sizes_to_exchange_f_to_recv(k)))
    !  IF(send_type_r(j) .NE.  MPI_DATATYPE_NULL) THEN
    !    print*,'send      ',int(z_coords),int(rank_to_send_to),int(i),int(sizes_to_exchange_r_to_send(j))
    !    n=n+1
    !    CALL MPI_ISEND(field(1,1,r_first_cell_to_send(j)),1_isp,send_type_r(j),      &
    !    rank_to_send_to,tag,comm,requests(n),errcode)

    !  !  CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
    !  ENDIF

    !  IF(recv_type_f(k) .NE. MPI_DATATYPE_NULL) THEN
    !    print*,'recv      ',int(z_coords),int(rank_to_recv_from),int(i),int(sizes_to_exchange_f_to_recv(k))
    !    n=n+1
    !    CALL MPI_IRecv(temp_recv(1,1,1),1_isp,recv_type_f(k),                &
    !    rank_to_recv_from,tag,comm,requests(n),errcode)
    !   ! CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
    ! ENDIF
    !    IF(n .GT. 0_idp)        CALL MPI_WAITALL(INT(n,isp), requests(1), MPI_STATUSES_IGNORE, errcode)
      IF(sizes_to_exchange_f_to_recv(k) .GT. 0_idp) THEN
!        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
        DO iz=1,sizes_to_exchange_f_to_recv(k)
          DO iy=iy_min_r, iy_max_r
            DO ix=ix_min_r, ix_max_r
               field_f(ix,iy,iz-1+f_first_cell_to_recv(k)) = temp_recv(ix-ix_min_r-nxg,iy-iy_min_r-nyg,iz)
            ENDDO
          ENDDO
        ENDDO
 !     !$OMP END PARALLEL DO
      ENDIF
      DEALLOCATE(temp_recv)
      CALL MPI_BARRIER(comm,errcode)
    ENDDO
#endif
  END SUBROUTINE  sendrecv_rf_generalized

  SUBROUTINE generalized_comms_group_f2r(coeff_norm)
#if defined(FFTW) 
    USE group_parameters
    USE mpi_fftw3
#endif
    USE time_stat
    USE shared_data
    USE load_balance  
    USE fields
#if defined(FFTW)
    REAL(num)  , INTENT(IN)  :: coeff_norm
    INTEGER(idp)  :: nxx, nyy, nzz
    REAL(num)     :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    nxx = 2*(nx_group/2+1)
    nyy = ny_group
    nzz = local_nz
    CALL sendrecv_fr_generalized(coeff_norm,ex,nx,nxguards,ny,nyguards,nz,nzguards,ex_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,ey,nx,nxguards,ny,nyguards,nz,nzguards,ey_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,ez,nx,nxguards,ny,nyguards,nz,nzguards,ez_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,bx,nx,nxguards,ny,nyguards,nz,nzguards,bx_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,by,nx,nxguards,ny,nyguards,nz,nzguards,by_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,bz,nx,nxguards,ny,nyguards,nz,nzguards,bz_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,jx,nx,nxguards,ny,nyguards,nz,nzguards,jx_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,jy,nx,nxguards,ny,nyguards,nz,nzguards,jy_r,nxx,nyy,nzz)
    CALL sendrecv_fr_generalized(coeff_norm,jz,nx,nxguards,ny,nyguards,nz,nzguards,jz_r,nxx,nyy,nzz)

    IF (it.ge.timestat_itstart) THEN
      localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
    ENDIF
#endif

  END SUBROUTINE generalized_comms_group_f2r

  SUBROUTINE sendrecv_fr_generalized(coeff_norm,field,nx1,nxg,ny1,nyg,nz1,nzg,field_f,nxx,nyy,nzz)

    USE load_balance
#if defined(FFTW)
    USE group_parameters
#endif
    REAL(num)   , INTENT(IN)                      :: coeff_norm
    INTEGER(idp), INTENT(IN)                      ::  nx1,nxg,ny1,nyg,nz1,nzg,nxx,nyy,nzz
    REAL(num)   , INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg)  :: field
    REAL(num)   , INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz)  :: field_f
    INTEGER(idp)        ::  i , j , k , ix , iy , iz
    REAL(num)   , ALLOCATABLE , DIMENSION(:,:,:)  :: temp_recv
    INTEGER(isp)                                  :: requests(2)
    INTEGER(isp)                                  :: rank_to_send_to, rank_to_recv_from
#if defined(FFTW)
    !  i-1 = mpi test in z direction for exchanges
    ! example :
    ! i-1=0 exchange between current rank and itself
    ! i-1=1 send to right (proc_z_max)
    ! i-1=nprocz-1 send to left(proc_z_min)
    DO i=1,nprocz
      IF(i == 1 ) CYCLE
      j = MODULO(z_coords+i-1,nprocz) +1
      ! j corresponds to the z_coords(+1) of mpi task to which the send is done 
      ! k corresponds to the z_coords(+1) of mpi task from which the recv is
      ! done
      k = MODULO(z_coords-(i-1),nprocz) +1

      rank_to_send_to = INT(array_of_ranks_to_send_to(i),isp)
      rank_to_recv_from = INT(array_of_ranks_to_recv_from(i),isp)
      IF(sizes_to_exchange_r_to_recv(k) == 0) rank_to_recv_from = MPI_PROC_NULL
      IF(sizes_to_exchange_f_to_send(j) == 0) rank_to_send_to = MPI_PROC_NULL
      ALLOCATE(temp_recv(nxx,nyy,sizes_to_exchange_r_to_recv(k)))
      CALL MPI_SENDRECV(field_f(1,1,f_first_cell_to_send(j)),1_isp,send_type_f(j),&
      rank_to_send_to,tag,temp_recv(1,1,1),1_isp,recv_type_r(k),rank_to_recv_from,tag,comm,status,errcode)
call mpi_barrier(comm,errcode)

!      IF(send_type_f(j) .NE. MPI_DATATYPE_NULL) THEN
!        CALL MPI_ISEND(field_f(1,1,f_first_cell_to_send(j)),1_isp,send_type_f(j),      &
!        rank_to_send_to,tag,comm,requests(1),errcode)
!        
!        CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
!      ENDIF
!      IF(recv_type_r(k) .NE. MPI_DATATYPE_NULL) THEN
!        CALL MPI_IRecv(temp_recv(1,1,1),1_isp,recv_type_r(k),                &
!        rank_to_recv_from,tag,comm,requests(1),errcode)
!        CALL MPI_WAITALL(1_isp, requests, MPI_STATUSES_IGNORE, errcode)
        IF(sizes_to_exchange_r_to_recv(k) .GT. 0) THEN
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
        DO iz=1,sizes_to_exchange_r_to_recv(k)
          DO iy=iy_min_r, iy_max_r
            DO ix=ix_min_r, ix_max_r
     field(ix-ix_min_r-nxg,iy-iy_min_r-nyg,iz-1+r_first_cell_to_recv(k)) = &
                 temp_recv(ix,iy,iz)*coeff_norm  !+&
        
      !  field(ix-ix_min_r-nxg,iy-iy_min_r-nyg,iz-nzguards+r_first_cell_to_recv(k))
            ENDDO
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
      DEALLOCATE(temp_recv)
    ENDDO
#endif
  END SUBROUTINE  sendrecv_fr_generalized


  
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



  ! ______________________________________________________________________________________
  !> @brief
  !> Move the global, local and tile boundaries when moving window
  !
  !> author
  !> Guillaume Blaclard
  !
  !> @date
  !> Creation: 2016
  !
  !> @param[in] dxx, dyy, dzz moving window displacement
  ! ______________________________________________________________________________________
  SUBROUTINE pxr_move_sim_boundaries(dxx, dyy, dzz)
    IMPLICIT NONE

    REAL(num), INTENT(IN) :: dxx, dyy, dzz
    INTEGER(idp) :: ispecies, ix, iy, iz
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER    :: curr_tile

    ! MOVE GLOBAL BOUNDARIES
    xmin=xmin+dxx
    xmax=xmax+dxx
    xmin_part=xmin_part+dxx
    xmax_part=xmax_part+dxx
    ! Update along Y
    ymin=ymin+dyy
    ymax=ymax+dyy
    ymin_part=ymin_part+dyy
    ymax_part=ymax_part+dyy
    ! Update along Z
    zmin=zmin+dzz
    zmax=zmax+dzz
    zmin_part=zmin_part+dzz
    zmax_part=zmax_part+dzz

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
    x_min_local_part=x_min_local_part+dxx
    x_max_local_part=x_max_local_part+dxx
    ! Update along Y
    y_min_local=y_min_local+dyy
    y_max_local=y_max_local+dyy
    y_min_local_part=y_min_local_part+dyy
    y_max_local_part=y_max_local_part+dyy
    ! Update along Z
    z_min_local=z_min_local+dzz
    z_max_local=z_max_local+dzz
    z_min_local_part=z_min_local_part+dzz
    z_max_local_part=z_max_local_part+dzz

    ! MOVE LOCAL BOUNDARIES
    x_grid_min_local=x_grid_min_local+dxx
    x_grid_max_local=x_grid_max_local+dxx
    ! Update along Y
    y_grid_min_local=y_grid_min_local+dyy
    y_grid_max_local=y_grid_max_local+dyy
    ! Update along Z
    z_grid_min_local=z_grid_min_local+dzz
    z_grid_max_local=z_grid_max_local+dzz

    ! MOVE TILE BOUNDARIES ALONG X, Y, Z
    DO ispecies =1, nspecies
      curr=> species_parray(ispecies)
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile=> curr%array_of_tiles(ix, iy, iz)
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

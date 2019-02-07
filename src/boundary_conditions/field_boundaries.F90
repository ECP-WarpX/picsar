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
  USE PICSAR_precision
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
    USE constants, ONLY: c_ndims
    USE mpi
    USE picsar_precision, ONLY: idp, isp, num
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

    IF(c_dim == 3) THEN
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
    ENDIF  
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
    IF(c_dim == 3) THEN
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
    ENDIF
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
  !
  !> @author
  !> Haithem Kallala
  !> 
  !> @date
  !> Creation 2017
  !> 
  !> @ brief 
  !> This routine is copying values from the local  field arrays to the 
  !> FFT field arrays. FFT arrays are used to perform distributed FFTs on each MPI 
  !> groups. Local field arrays are used to perform field gathering operations during 
  !> the PIC cycle. Keeping two different grids allows for efficient load balancing 
  !> of computations in the field gathering and Maxwell solvers stages of the PIC cycle. 
  ! ______________________________________________________________________________________
  SUBROUTINE generalized_comms_group_l2g()
    USE fields, ONLY: bx, bx_r, bxy, bxy_r, bxz, bxz_r, by, by_r, byx, byx_r, byz,   &
      byz_r, bz, bz_r, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exy, exy_r, exz, exz_r, ey, &
      ey_r, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezx, ezx_r, ezy, ezy_r, jx, jx_r, jy,  &
      jy_r, jz, jz_r, nxguards, nyguards, nzguards, rho_r, rhoold_r
#if defined(FFTW) 
    USE iso_c_binding
    USE load_balance
#endif
    USE mpi
#if defined(FFTW) 
    USE mpi_fftw3, ONLY: local_nx, local_ny, local_nz
#endif
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs, nx, ny, nz, rho, rhoold
    USE time_stat, ONLY: localtimes, timestat_itstart
#if defined(FFTW)
    INTEGER(idp)  :: nxx, nyy, nzz
    REAL(num)                                   :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    nxx = local_nx 
    nyy = local_ny
    nzz = local_nz
    IF(mpicom_curr == 1) THEN
      IF(.NOT. absorbing_bcs) THEN

        !> When using period bcs for fields, standard EM equations are solved
        !> Standard EM fields are exchanged between hybrid and local fields

        CALL sendrecv_l2g_generalized(ex,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ex_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(ey,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ey_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(ez,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ez_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(bx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(by,nx,nxguards,ny,nyguards,nz,nzguards,          &
        by_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(bz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(jx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        jx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(jy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        jy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(jz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        jz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(rho,nx,nxguards,ny,nyguards,nz,nzguards,         &
        rho_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(rhoold,nx,nxguards,ny,nyguards,nz,nzguards,      &
        rhoold_r,nxx,nyy,nzz)
      ELSE IF (absorbing_bcs) THEN

        !> When using absorbing bcs for fields,splitted fields EM equations are
        !solved
        !> Splitted EM fields are exchanged between hybrid and local fields

        CALL sendrecv_l2g_generalized(exy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        exy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(eyx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        eyx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(ezx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ezx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(bxy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bxy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(byx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        byx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(bzx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bzx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(exz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        exz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(eyz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        eyz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(ezy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ezy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(bxz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bxz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(byz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        byz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(bzy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bzy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(jx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        jx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(jy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        jy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(jz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        jz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(rho,nx,nxguards,ny,nyguards,nz,nzguards,         &
        rho_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized(rhoold,nx,nxguards,ny,nyguards,nz,nzguards,      &
        rhoold_r,nxx,nyy,nzz)
      ENDIF
    ELSE 
      IF(.NOT. absorbing_bcs) THEN

        !> When using period bcs for fields, standard EM equations are solved
        !> Standard EM fields are exchanged between hybrid and local fields

        CALL sendrecv_l2g_generalized_non_blocking(ex,nx,nxguards,ny,nyguards,nz,      &
        nzguards,ex_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(ey,nx,nxguards,ny,nyguards,nz,      &
        nzguards,ey_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(ez,nx,nxguards,ny,nyguards,nz,      &
        nzguards,ez_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(bx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,bx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(by,nx,nxguards,ny,nyguards,nz,      &
        nzguards,by_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(bz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,bz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(jx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,jx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(jy,nx,nxguards,ny,nyguards,nz,      &
        nzguards,jy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(jz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,jz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(rho,nx,nxguards,ny,nyguards,nz,     &
        nzguards,rho_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(rhoold,nx,nxguards,ny,nyguards,nz,  &
        nzguards,rhoold_r,nxx,nyy,nzz)
      ELSE IF (absorbing_bcs) THEN

        !> When using absorbing bcs for fields,splitted fields EM equations are solved
        !> Splitted EM fields are exchanged between hybrid and local fields

        CALL sendrecv_l2g_generalized_non_blocking(exy,nx,nxguards,ny,nyguards,nz,      &
        nzguards,exy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(eyx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,eyx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(ezx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,ezx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(bxy,nx,nxguards,ny,nyguards,nz,      &
        nzguards,bxy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(byx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,byx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(bzx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,bzx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(exz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,exz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(eyz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,eyz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(ezy,nx,nxguards,ny,nyguards,nz,      &
        nzguards,ezy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(bxz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,bxz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(byz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,byz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(bzy,nx,nxguards,ny,nyguards,nz,      &
        nzguards,bzy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(jx,nx,nxguards,ny,nyguards,nz,      &
        nzguards,jx_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(jy,nx,nxguards,ny,nyguards,nz,      &
        nzguards,jy_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(jz,nx,nxguards,ny,nyguards,nz,      &
        nzguards,jz_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(rho,nx,nxguards,ny,nyguards,nz,     &
        nzguards,rho_r,nxx,nyy,nzz)
        CALL sendrecv_l2g_generalized_non_blocking(rhoold,nx,nxguards,ny,nyguards,nz,  &
        nzguards,rhoold_r,nxx,nyy,nzz)
      ENDIF

    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
    ENDIF
#endif
  END SUBROUTINE generalized_comms_group_l2g

  ! ______________________________________________________________________________________
  !> @brief 
  !> Routine for exchanging data from the local grid arrays to the FFT grid arrays 
  !> This routine uses blocking SENDRECV MPI exchanges
  !> @author
  !> Haithem Kallala
  !> @params[in,out] 3D REAL(num) array field_l - local field array (e.g. ex)
  !> @params[in,out] 3D REAL(num) array field_g - FFT field array (e.g. ex_r)
  !> @params[in] nx1 - number of cells along dimension X of field_l
  !> @params[in] nxg - number of guard cells along dimension X of field_l
  !> @params[in] ny1 - number of cells along dimension Y of field_l
  !> @params[in] nyg - number of guard cells along dimension Y of field_l
  !> @params[in] nz1 - number of cells along dimension Z of field_l
  !> @params[in] nzg - number of guard cells along dimension Z of field_l
  !> @params[in] nxx - size of array field_g along X 
  !> @params[in] nyy - size of array field_g along Y
  !> @params[in] nyz - size of array field_g along Y
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE sendrecv_l2g_generalized(field_l,nx1,nxg,ny1,nyg,nz1,nzg,field_g,nxx,nyy,nzz)
#if defined(FFTW)
    USE group_parameters, ONLY: array_of_ranks_to_recv_from_l2g,                     &
      array_of_ranks_to_send_to_l2g, g_first_cell_to_recv_y, g_first_cell_to_recv_z, &
      l_first_cell_to_send_y, l_first_cell_to_send_z, nb_comms_l2g, recv_type_g,     &
      send_type_l, size_exchanges_l2g_recv_y, size_exchanges_l2g_recv_z,             &
      size_exchanges_l2g_send_y, size_exchanges_l2g_send_z
    USE load_balance
#endif
    USE mpi
#if defined(FFTW)
    USE mpi_type_constants, ONLY: status
    USE picsar_precision, ONLY: idp, isp, num
#endif
    INTEGER(idp), INTENT(IN)                    ::  nx1,nxg,ny1,nyg,nz1,nzg,nxx,nyy,nzz
    REAL(num)    ,INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg)   & 
    :: field_l
    REAL(num)    ,INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz)  :: field_g
    INTEGER(idp)                                ::  ii
    INTEGER(isp)                                :: rank_to_send_to, rank_to_recv_from
#if defined(FFTW)
    DO ii=1,nb_comms_l2g
     ! ii==1 corresponds to target_rank = my_rank 
     ! and recv_rank = my_rank 
     ! this copy is done locally in fourier_psaod.F90

      IF (ii==1) CYCLE  
                         
      rank_to_send_to = INT(array_of_ranks_to_send_to_l2g(ii),isp)
      rank_to_recv_from = INT(array_of_ranks_to_recv_from_l2g(ii),isp)

      IF(size_exchanges_l2g_recv_z(ii) == 0) rank_to_recv_from = MPI_PROC_NULL
      IF(size_exchanges_l2g_send_z(ii) == 0) rank_to_send_to = MPI_PROC_NULL
      IF(size_exchanges_l2g_recv_y(ii) == 0) rank_to_recv_from =MPI_PROC_NULL
      IF(size_exchanges_l2g_send_y(ii) == 0) rank_to_send_to = MPI_PROC_NULL
      CALL MPI_SENDRECV(field_l(-nxg, l_first_cell_to_send_y(ii),                      &
      l_first_cell_to_send_z(ii)), 1_isp, send_type_l(ii),                             &
      rank_to_send_to,tag,field_g(1, g_first_cell_to_recv_y(ii),                       &
      g_first_cell_to_recv_z(ii)), 1_isp, recv_type_g(ii),                             &
      rank_to_recv_from ,tag ,comm ,status ,errcode)
    ENDDO
#endif
  END SUBROUTINE  sendrecv_l2g_generalized

  ! ______________________________________________________________________________________
  !> @brief 
  !> Routine for exchanging data from the local grid arrays to the FFT grid arrays 
  !> This routine uses ISEND and IRECV MPI exchanges (non-blocking exchanges) 
  !> @author
  !> Haithem Kallala
  !> @params[in,out] 3D REAL(num) array field_l - local field array (e.g. ex)
  !> @params[in,out] 3D REAL(num) array field_g - FFT field array (e.g. ex_r)
  !> @params[in] nx1 - number of cells along dimension X of field_l
  !> @params[in] nxg - number of guard cells along dimension X of field_l
  !> @params[in] ny1 - number of cells along dimension Y of field_l
  !> @params[in] nyg - number of guard cells along dimension Y of field_l
  !> @params[in] nz1 - number of cells along dimension Z of field_l
  !> @params[in] nzg - number of guard cells along dimension Z of field_l
  !> @params[in] nxx - size of array field_g along X 
  !> @params[in] nyy - size of array field_g along Y
  !> @params[in] nyz - size of array field_g along Y
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE sendrecv_l2g_generalized_non_blocking(field_l,nx1,nxg,ny1,nyg,nz1,nzg,    &
  field_g,nxx,nyy,nzz)
#if defined(FFTW)
USE group_parameters, ONLY: array_of_ranks_to_recv_from_l2g,                         &
  array_of_ranks_to_send_to_l2g, g_first_cell_to_recv_y, g_first_cell_to_recv_z,     &
  l_first_cell_to_send_y, l_first_cell_to_send_z, nb_comms_l2g, recv_type_g,         &
  requests_l2g, send_type_l, size_exchanges_l2g_recv_y, size_exchanges_l2g_recv_z,   &
  size_exchanges_l2g_send_y, size_exchanges_l2g_send_z
USE load_balance
#endif
USE mpi
#if defined(FFTW)
USE picsar_precision, ONLY: idp, isp, num
#endif

    INTEGER(idp), INTENT(IN)                    ::  nx1,nxg,ny1,nyg,nz1,nzg,nxx,nyy,nzz
    REAL(num)    ,INTENT(INOUT)  , DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg)   &
    :: field_l
    REAL(num)    ,INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz)  :: field_g
    INTEGER(idp)                                ::  ii
    INTEGER(isp)                                :: rank_to_send_to, rank_to_recv_from
    INTEGER(idp)                                :: n

#if defined(FFTW)
    requests_l2g = 0
    n=0
    DO ii=1,nb_comms_l2g
     ! ii==1 corresponds to target_rank = my_rank 
     ! and recv_rank = my_rank 
     ! this copy is done locally in fourier_psaod.F90
      IF(ii==1) CYCLE
      rank_to_send_to = INT(array_of_ranks_to_send_to_l2g(ii),isp)
      rank_to_recv_from = INT(array_of_ranks_to_recv_from_l2g(ii),isp)

      IF(size_exchanges_l2g_recv_z(ii) .GT. 0 .AND.                                    &
      size_exchanges_l2g_recv_y(ii) .GT. 0)   THEN                                     
        n=n+1
        CALL MPI_IRECV(field_g(1, g_first_cell_to_recv_y(ii),                          &
        g_first_cell_to_recv_z(ii)), 1_isp,recv_type_g(ii),                            &
        rank_to_recv_from,tag,comm,requests_l2g(n),errcode)
      ENDIF
      IF(size_exchanges_l2g_send_z(ii) .GT. 0                                          &
      .AND. size_exchanges_l2g_send_y(ii) .GT. 0) THEN
        n=n+1
        CALL MPI_ISEND(field_l(-nxg, l_first_cell_to_send_y(ii),                       &
         l_first_cell_to_send_z(ii)), 1_isp, send_type_l(ii)                           &
        ,rank_to_send_to,tag,comm,requests_l2g(n),errcode)
      ENDIF
    ENDDO
    CALL MPI_WAITALL(INT(n,isp),requests_l2g, MPI_STATUSES_IGNORE, errcode)
#endif
  END SUBROUTINE  sendrecv_l2g_generalized_non_blocking
  
  ! ______________________________________________________________________________________
  !
  !> @author
  !> Haithem Kallala
  !> 
  !> @date
  !> Creation 2017
  !> 
  !> @ brief 
  !> This routine is copying values from the FFT field arrays to the 
  !> local field arrays. FFT arrays are used to perform distributed FFTs on each MPI 
  !> groups. Local field arrays are used to perform field gathering operations during 
  !> the PIC cycle. Keeping two different grids allows for efficient load balancing 
  !> of computations in the field gathering and Maxwell solvers stages of the PIC cycle. 
  ! ______________________________________________________________________________________
  SUBROUTINE generalized_comms_group_g2l()
    USE fields, ONLY: bx, bx_r, bxy, bxy_r, bxz, bxz_r, by, by_r, byx, byx_r, byz,   &
      byz_r, bz, bz_r, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exy, exy_r, exz, exz_r, ey, &
      ey_r, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezx, ezx_r, ezy, ezy_r, nxguards,      &
      nyguards, nzguards
#if defined(FFTW) 
    USE iso_c_binding
    USE load_balance
#endif
    USE mpi
#if defined(FFTW) 
    USE mpi_fftw3, ONLY: local_nx, local_ny, local_nz
#endif
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs, nx, ny, nz
    USE time_stat, ONLY: localtimes, timestat_itstart
    INTEGER(idp)  :: nxx, nyy, nzz
    REAL(num)     :: tmptime
#if defined(FFTW)
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    nxx = local_nx 
    nyy = local_ny
    nzz = local_nz
    IF(mpicom_curr ==1) THEN
      IF(.NOT. absorbing_bcs) THEN

        !> When using period bcs for fields, standard EM equations are solved
        !> Standard EM fields are exchanged between hybrid and local fields

        CALL sendrecv_g2l_generalized(ex,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ex_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(ey,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ey_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(ez,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ez_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(bx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(by,nx,nxguards,ny,nyguards,nz,nzguards,          &
        by_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(bz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bz_r,nxx,nyy,nzz)
      ELSE IF(absorbing_bcs) THEN

        !> When using absorbing bcs for fields,splitted fields EM equations are
        !solved
        !> Splitted EM fields are exchanged between hybrid and local fields

        CALL sendrecv_g2l_generalized(exy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        exy_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(eyx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        eyx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(ezx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ezx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(bxy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bxy_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(byx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        byx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(bzx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bzx_r,nxx,nyy,nzz) 
        CALL sendrecv_g2l_generalized(exz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        exz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(eyz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        eyz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(ezy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ezy_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(bxz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bxz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(byz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        byz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized(bzy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bzy_r,nxx,nyy,nzz)
      ENDIF
    ELSE 
      IF (.NOT. absorbing_bcs) THEN

        !> When using period bcs for fields, standard EM equations are solved
        !> Standard EM fields are exchanged between hybrid and local fields

        CALL sendrecv_g2l_generalized_non_blocking(ex,nx,nxguards,ny,nyguards,         &
        nz,nzguards,ex_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(ey,nx,nxguards,ny,nyguards,         &
        nz,nzguards,ey_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(ez,nx,nxguards,ny,nyguards,         &
        nz,nzguards,ez_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(bx,nx,nxguards,ny,nyguards,         &
        nz,nzguards,bx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(by,nx,nxguards,ny,nyguards,         &
        nz,nzguards,by_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(bz,nx,nxguards,ny,nyguards,         &
        nz,nzguards,bz_r,nxx,nyy,nzz)
      ELSE IF(absorbing_bcs) THEN

        !> When using absorbing bcs for fields,splitted fields EM equations are
        !solved
        !> Splitted EM fields are exchanged between hybrid and local fields

        CALL sendrecv_g2l_generalized_non_blocking(exy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        exy_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(eyx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        eyx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(ezx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ezx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(bxy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bxy_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(byx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        byx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(bzx,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bzx_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(exz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        exz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(eyz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        eyz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(ezy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        ezy_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(bxz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bxz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(byz,nx,nxguards,ny,nyguards,nz,nzguards,          &
        byz_r,nxx,nyy,nzz)
        CALL sendrecv_g2l_generalized_non_blocking(bzy,nx,nxguards,ny,nyguards,nz,nzguards,          &
        bzy_r,nxx,nyy,nzz)

      ENDIF
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(25) = localtimes(25) + (MPI_WTIME() - tmptime)
    ENDIF
#endif

  END SUBROUTINE generalized_comms_group_g2l

  ! ______________________________________________________________________________________
  !> @brief 
  !> Routine for exchanging data from the FFT grid arrays to the local grid arrays 
  !> This routine uses ISEND and IRECV MPI exchanges (non-blocking exchanges) 
  !> @author
  !> Haithem Kallala
  !> @params[in,out] 3D REAL(num) array field_l - local field array (e.g. ex)
  !> @params[in,out] 3D REAL(num) array field_g - FFT field array (e.g. ex_r)
  !> @params[in] nx1 - number of cells along dimension X of field_l
  !> @params[in] nxg - number of guard cells along dimension X of field_l
  !> @params[in] ny1 - number of cells along dimension Y of field_l
  !> @params[in] nyg - number of guard cells along dimension Y of field_l
  !> @params[in] nz1 - number of cells along dimension Z of field_l
  !> @params[in] nzg - number of guard cells along dimension Z of field_l
  !> @params[in] nxx - size of array field_g along X 
  !> @params[in] nyy - size of array field_g along Y
  !> @params[in] nyz - size of array field_g along Y
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE sendrecv_g2l_generalized_non_blocking(field_l,nx1,nxg,ny1,nyg,nz1,nzg,    &
  field_g,nxx,nyy,nzz)
#if defined(FFTW)
USE group_parameters, ONLY: array_of_ranks_to_recv_from_g2l,                         &
  array_of_ranks_to_send_to_g2l, g_first_cell_to_send_y, g_first_cell_to_send_z,     &
  l_first_cell_to_recv_y, l_first_cell_to_recv_z, nb_comms_g2l, recv_type_l,         &
  requests_g2l, send_type_g, size_exchanges_g2l_recv_y, size_exchanges_g2l_recv_z,   &
  size_exchanges_g2l_send_y, size_exchanges_g2l_send_z
USE load_balance
#endif
USE mpi
#if defined(FFTW)
USE picsar_precision, ONLY: idp, isp, num
#endif

    INTEGER(idp), INTENT(IN)                      ::  nx1,nxg,ny1,nyg,nz1,nzg,nxx,nyy,nzz
    REAL(num)   , INTENT(INOUT)  ,                                                     &
    DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg)  :: field_l
    REAL(num)   , INTENT(INOUT)  ,                                                     &
    DIMENSION(nxx,nyy,nzz)  :: field_g
    INTEGER(idp)        :: ii
    INTEGER(isp)                                  :: rank_to_send_to, rank_to_recv_from
    INTEGER(idp)                                :: n
   
#if defined(FFTW)
    requests_g2l = 0 
    n = 0
    DO ii=1,nb_comms_g2l
     ! ii==1 corresponds to target_rank = my_rank 
     ! and recv_rank = my_rank 
     ! this copy is done locally in fourier_psaod.F90

      IF(ii==1) CYCLE
      rank_to_send_to = INT(array_of_ranks_to_send_to_g2l(ii),isp)
      rank_to_recv_from = INT(array_of_ranks_to_recv_from_g2l(ii),isp)
 
      IF(size_exchanges_g2l_recv_z(ii) .GT. 0 .AND.                                    &
      size_exchanges_g2l_recv_y(ii) .GT. 0) THEN                                       
        n=n+1
        CALL MPI_IRECV(field_l(-nxg, l_first_cell_to_recv_y(ii),                       &
        l_first_cell_to_recv_z(ii)), 1_isp,recv_type_l(ii)                             &
        , rank_to_recv_from,tag,comm,requests_g2l(n),errcode)
      ENDIF
      IF(size_exchanges_g2l_send_z(ii) .GT. 0 .AND.                                    &
      size_exchanges_g2l_send_y(ii) .GT. 0 ) THEN
        n=n+1
        CALL MPI_ISEND(field_g(1,g_first_cell_to_send_y(ii),                           &
        g_first_cell_to_send_z(ii)) ,1_isp,send_type_g(ii),                            &
        rank_to_send_to,tag,comm,requests_g2l(n),errcode)   
      ENDIF

    ENDDO
    CALL MPI_WAITALL(INT(n,isp),requests_g2l, MPI_STATUSES_IGNORE, errcode)
#endif
  END SUBROUTINE  sendrecv_g2l_generalized_non_blocking

  ! ______________________________________________________________________________________
  !> @brief 
  !> Routine for exchanging data from the FFT grid arrays to the local grid arrays 
  !> This routine uses blocking SENDRECV MPI exchanges 
  !> @author
  !> Haithem Kallala
  !> @params[in,out] 3D REAL(num) array field_l - local field array (e.g. ex)
  !> @params[in,out] 3D REAL(num) array field_g - FFT field array (e.g. ex_r)
  !> @params[in] nx1 - number of cells along dimension X of field_l
  !> @params[in] nxg - number of guard cells along dimension X of field_l
  !> @params[in] ny1 - number of cells along dimension Y of field_l
  !> @params[in] nyg - number of guard cells along dimension Y of field_l
  !> @params[in] nz1 - number of cells along dimension Z of field_l
  !> @params[in] nzg - number of guard cells along dimension Z of field_l
  !> @params[in] nxx - size of array field_g along X 
  !> @params[in] nyy - size of array field_g along Y
  !> @params[in] nyz - size of array field_g along Y
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE sendrecv_g2l_generalized(field_l,nx1,nxg,ny1,nyg,nz1,nzg,field_g,nxx,nyy,nzz)
#if defined(FFTW)
    USE group_parameters, ONLY: array_of_ranks_to_recv_from_g2l,                     &
      array_of_ranks_to_send_to_g2l, g_first_cell_to_send_y, g_first_cell_to_send_z, &
      l_first_cell_to_recv_y, l_first_cell_to_recv_z, nb_comms_g2l, recv_type_l,     &
      send_type_g, size_exchanges_g2l_recv_y, size_exchanges_g2l_recv_z,             &
      size_exchanges_g2l_send_y, size_exchanges_g2l_send_z
    USE load_balance
#endif
    USE mpi
#if defined(FFTW)
    USE mpi_type_constants, ONLY: status
    USE picsar_precision, ONLY: idp, isp, num
#endif
    INTEGER(idp), INTENT(IN)                      :: nx1,nxg,ny1,nyg,nz1,nzg,nxx,nyy,nzz
    REAL(num), INTENT(INOUT),                                                          &
    DIMENSION(-nxg:nx1+nxg,-nyg:ny1+nyg,-nzg:nz1+nzg)  :: field_l
    REAL(num)   , INTENT(INOUT)  , DIMENSION(nxx,nyy,nzz)  :: field_g
    INTEGER(idp)        ::ii
    INTEGER(isp)        :: rank_to_send_to,rank_to_recv_from
#if defined(FFTW)

    DO ii=1,nb_comms_g2l
     ! ii==1 corresponds to target_rank = my_rank 
     ! and recv_rank = my_rank 
     ! this copy is done locally in fourier_psaod.F90


      IF(ii==1) CYCLE
      rank_to_send_to = INT(array_of_ranks_to_send_to_g2l(ii),isp)
      rank_to_recv_from = INT(array_of_ranks_to_recv_from_g2l(ii),isp)

      IF(size_exchanges_g2l_recv_z(ii) == 0) rank_to_recv_from = MPI_PROC_NULL
      IF(size_exchanges_g2l_send_z(ii) == 0) rank_to_send_to = MPI_PROC_NULL
      IF(size_exchanges_g2l_recv_y(ii) == 0) rank_to_recv_from = MPI_PROC_NULL
      IF(size_exchanges_g2l_send_y(ii) == 0) rank_to_send_to = MPI_PROC_NULL

      CALL MPI_SENDRECV(field_g(1,g_first_cell_to_send_y(ii),                          &
      g_first_cell_to_send_z(ii)) ,1_isp,send_type_g(ii),                              &
      rank_to_send_to,tag,field_l(-nxg, l_first_cell_to_recv_y(ii),                    &
      l_first_cell_to_recv_z(ii)), 1_isp,&
      recv_type_l(ii), rank_to_recv_from, tag, comm, status, errcode)

    ENDDO
#endif
  END SUBROUTINE  sendrecv_g2l_generalized
    
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
    IF(c_dim == 3) THEN
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
    ENDIF
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
    USE communications, ONLY: reqperjxx, reqperjxy, reqperjxz
    USE constants, ONLY: c_ndims
    USE mpi
    USE picsar_precision, ONLY: idp, isp, num
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
  END SUBROUTINE summation_bcs_persistent_jx

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
    USE communications, ONLY: reqperjyx, reqperjyy, reqperjyz
    USE constants, ONLY: c_ndims
    USE mpi
    USE picsar_precision, ONLY: idp, isp, num
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
  END SUBROUTINE summation_bcs_persistent_jy

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
    USE communications, ONLY: reqperjzx, reqperjzy, reqperjzz
    USE constants, ONLY: c_ndims
    USE mpi
    USE picsar_precision, ONLY: idp, isp, num
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
  END SUBROUTINE summation_bcs_persistent_jz

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
    USE mpi
    REAL(num) :: tmptime
#if defined(DEBUG)
    WRITE(0, *) "efield_bcs: start"
#endif

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Electric field MPI exchange between subdomains
    IF(.NOT. absorbing_bcs) THEN
    
      !> When using periodic bcs, exchange standard EM fields

      CALL field_bc(ex, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(ey, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(ez, nxguards, nyguards, nzguards, nx, ny, nz)
    ELSE IF(absorbing_bcs) THEN

      !> When using absorbing bcs, exchange splitted EM fields
      CALL field_bc(exy, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(exz, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(eyx, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(eyz, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(ezx, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(ezy, nxguards, nyguards, nzguards, nx, ny, nz)
      !> When using absorbing bcs, the electric field is merged here
      !> This is done here in order not to call merge_e_fields from warp 
      CALL merge_e_fields()

    ENDIF
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
USE mpi

    REAL(num) :: tmptime
#if defined(DEBUG)
    WRITE(0, *) "bfield_bcs: start"
#endif
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Magnetic field MPI exchange between subdomains
    IF(.NOT. absorbing_bcs) THEN 

      !> When using periodic bcs, exchange standard EM fields

      CALL field_bc(bx, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(by, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(bz, nxguards, nyguards, nzguards, nx, ny, nz)
    ELSE IF(absorbing_bcs) THEN

      !> When using absorbing bcs, exchange splitted EM fields

      CALL field_bc(bxy, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(bxz, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(byx, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(byz, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(bzx, nxguards, nyguards, nzguards, nx, ny, nz)
      CALL field_bc(bzy, nxguards, nyguards, nzguards, nx, ny, nz)

      !> When using absorbing bcs, the magnetic field is merged here
      !> This is done here in order not to call merge_fields from warp 

      CALL merge_b_fields()
    ENDIF
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
    USE mpi
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
    USE mpi
    USE picsar_precision, ONLY: num
    USE time_stat, ONLY: localtimes, timestat_itstart

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

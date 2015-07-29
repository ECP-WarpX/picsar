MODULE mpi_routines
  USE shared_data
  USE fields
  USE mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(num) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init

    CALL MPI_INIT(errcode)
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    CALL MPI_COMM_SIZE(comm, nproc, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)

  END SUBROUTINE mpi_minimal_init



  SUBROUTINE setup_communicator

    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: dims(ndims), idim, old_comm, ierr
    LOGICAL :: periods(ndims), reorder, op, reset
    INTEGER :: test_coords(ndims)
    INTEGER :: ix, iy, iz
    INTEGER :: nxsplit, nysplit, nzsplit
    INTEGER :: area, minarea, nprocyz
    INTEGER :: ranges(3,1), nproc_orig, oldgroup, newgroup
    CHARACTER(LEN=11) :: str



    !!! --- NB: CPU Split performed on number of grid points (not cells)
    nx_global_grid = nx_global+1
    ny_global_grid = ny_global+1
    nz_global_grid = nz_global+1

    nproc_orig = nproc

    IF (nx_global_grid .LT. nxguards .OR. ny_global_grid .LT. nyguards .OR. nz_global_grid .LT. nzguards) THEN
      IF (rank .EQ. 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Simulation domain is too small.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

    reset = .FALSE.
    IF (MAX(nprocx,1) * MAX(nprocy,1) * MAX(nprocz,1) .GT. nproc) THEN
      reset = .TRUE.
    ELSE IF (nprocx * nprocy * nprocz .GT. 0) THEN
      ! Sanity check
      nxsplit = nx_global_grid / nprocx
      nysplit = ny_global_grid / nprocy
      nzsplit = nz_global_grid / nprocz
      IF (nxsplit .LT. nxguards .OR. nysplit .LT. nyguards .OR. nzsplit .LT. nzguards) &
          reset = .TRUE.
    ENDIF

    IF (reset) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, 'Unable to use requested processor subdivision. Using ' &
            // 'default division.'
      ENDIF
      nprocx = 0
      nprocy = 0
      nprocz = 0
    ENDIF

    IF (nprocx * nprocy * nprocz .EQ. 0) THEN
      DO WHILE (nproc .GT. 1)
        ! Find the processor split which minimizes surface area of
        ! the resulting domain

        minarea = nx_global_grid * ny_global_grid + ny_global_grid * nz_global_grid &
            + nz_global_grid * nx_global_grid

        DO ix = 1, nproc
          nprocyz = nproc / ix
          IF (ix * nprocyz .NE. nproc) CYCLE

          nxsplit = nx_global_grid / ix
          ! Actual domain must be bigger than the number of ghostcells
          IF (nxsplit .LT. nxguards) CYCLE

          DO iy = 1, nprocyz
            iz = nprocyz / iy
            IF (iy * iz .NE. nprocyz) CYCLE

            nysplit = ny_global_grid / iy
            nzsplit = nz_global_grid / iz
            ! Actual domain must be bigger than the number of ghostcells
            IF (nysplit .LT. nyguards .OR. nzsplit .LT. nzguards) CYCLE

            area = nxsplit * nysplit + nysplit * nzsplit + nzsplit * nxsplit
            IF (area .LT. minarea) THEN
              nprocx = ix
              nprocy = iy
              nprocz = iz
              minarea = area
            ENDIF
          ENDDO
        ENDDO

        IF (nprocx .GT. 0) EXIT

        ! If we get here then no suitable split could be found. Decrease the
        ! number of processors and try again.

        nproc = nproc - 1
      ENDDO
    ENDIF

    IF (nproc_orig .NE. nproc) THEN
      IF (.NOT.allow_cpu_reduce) THEN
        IF (rank .EQ. 0) THEN
          PRINT*,'*** ERROR ***'
          PRINT*,'Cannot split the domain using the requested number of CPUs.'
          PRINT*,'Try reducing the number of CPUs'
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        STOP
      ENDIF
      IF (rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Cannot split the domain using the requested number of CPUs.'
        PRINT*,'Reducing the number of CPUs'
      ENDIF
      ranges(1,1) = nproc
      ranges(2,1) = nproc_orig - 1
      ranges(3,1) = 1
      old_comm = comm
      CALL MPI_COMM_GROUP(old_comm, oldgroup, errcode)
      CALL MPI_GROUP_RANGE_EXCL(oldgroup, 1, ranges, newgroup, errcode)
      CALL MPI_COMM_CREATE(old_comm, newgroup, comm, errcode)
      IF (comm .EQ. MPI_COMM_NULL) THEN
        CALL MPI_FINALIZE(errcode)
        STOP
      ENDIF
      CALL MPI_GROUP_FREE(oldgroup, errcode)
      CALL MPI_GROUP_FREE(newgroup, errcode)
      CALL MPI_COMM_FREE(old_comm, errcode)
    ENDIF

    dims = (/nprocz, nprocy, nprocx/)
    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    periods = .FALSE.
    reorder = .TRUE.

    ! Set boundary to be periodic in x,y,z for particles and fields

    periods(c_ndims) = .TRUE.
    periods(c_ndims-1) = .TRUE.
    periods(c_ndims-2) = .TRUE.

    old_comm = comm
    CALL MPI_CART_CREATE(old_comm, ndims, dims, periods, reorder, comm, errcode)
    CALL MPI_COMM_FREE(old_comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_y_min, proc_y_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_z_min, proc_z_max, errcode)

    nprocx = dims(3)
    nprocy = dims(2)
    nprocz = dims(1)
    nprocdir = dims

    IF (rank .EQ. 0) THEN
      PRINT *, 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
    ENDIF

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
          ! MPI_CART_RANK returns an error rather than
          ! MPI_PROC_NULL if the coords are out of range.
          DO idim = 1, ndims
            IF ((test_coords(idim) .LT. 0 &
                .OR. test_coords(idim) .GE. dims(idim)) &
                .AND. .NOT. periods(idim)) op = .FALSE.
          ENDDO
          IF (op) THEN
            CALL MPI_CART_RANK(comm, test_coords, neighbour(ix,iy,iz), errcode)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: idim
    INTEGER :: nx0, nxp
    INTEGER :: ny0, nyp
    INTEGER :: nz0, nzp
    INTEGER :: iproc, ix, iy, iz

    ! Init number of guard cells of subdomains in each dimension
    nxguards = MAX(nox,norderx)+npass(1)
    nyguards = MAX(noy,nordery)+npass(2)
    nzguards = MAX(noz,norderz)+npass(3)
    IF (l_smooth_compensate) THEN
        nxguards = nxguards + 1
        nyguards = nyguards + 1
        nzguards = nzguards + 1
    END IF

    CALL setup_communicator

    ALLOCATE(npart_each_rank(nproc))
    ALLOCATE(x_grid_mins(0:nprocx-1), x_grid_maxs(0:nprocx-1))
    ALLOCATE(y_grid_mins(0:nprocy-1), y_grid_maxs(0:nprocy-1))
    ALLOCATE(z_grid_mins(0:nprocz-1), z_grid_maxs(0:nprocz-1))
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))
    ALLOCATE(cell_y_min(nprocy), cell_y_max(nprocy))
    ALLOCATE(cell_z_min(nprocz), cell_z_max(nprocz))

    nx0 = nx_global_grid / nprocx
    ny0 = ny_global_grid / nprocy
    nz0 = nz_global_grid / nprocz

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx .NE. nx_global_grid) THEN
        nxp = (nx0 + 1) * nprocx - nx_global_grid
    ELSE
        nxp = nprocx
    ENDIF

    IF (ny0 * nprocy .NE. ny_global_grid) THEN
        nyp = (ny0 + 1) * nprocy - ny_global_grid
    ELSE
        nyp = nprocy
    ENDIF

    IF (nz0 * nprocz .NE. nz_global_grid) THEN
        nzp = (nz0 + 1) * nprocz - nz_global_grid
    ELSE
        nzp = nprocz
    ENDIF

    DO idim = 1, nxp
        cell_x_min(idim) = (idim - 1) * nx0 + 1
        cell_x_max(idim) = idim * nx0
    ENDDO
    DO idim = nxp + 1, nprocx
        cell_x_min(idim) = nxp * nx0 + (idim - nxp - 1) * (nx0 + 1) + 1
        cell_x_max(idim) = nxp * nx0 + (idim - nxp) * (nx0 + 1)
    ENDDO

    DO idim = 1, nyp
        cell_y_min(idim) = (idim - 1) * ny0 + 1
        cell_y_max(idim) = idim * ny0
    ENDDO
    DO idim = nyp + 1, nprocy
        cell_y_min(idim) = nyp * ny0 + (idim - nyp - 1) * (ny0 + 1) + 1
        cell_y_max(idim) = nyp * ny0 + (idim - nyp) * (ny0 + 1)
    ENDDO

    DO idim = 1, nzp
        cell_z_min(idim) = (idim - 1) * nz0 + 1
        cell_z_max(idim) = idim * nz0
    ENDDO
    DO idim = nzp + 1, nprocz
        cell_z_min(idim) = nzp * nz0 + (idim - nzp - 1) * (nz0 + 1) + 1
        cell_z_max(idim) = nzp * nz0 + (idim - nzp) * (nz0 + 1)

    ENDDO

    nx_global_grid_min = cell_x_min(x_coords+1)
    nx_global_grid_max = cell_x_max(x_coords+1)
    n_global_grid_min(1) = nx_global_grid_min
    n_global_grid_max(1) = nx_global_grid_max

    ny_global_grid_min = cell_y_min(y_coords+1)
    ny_global_grid_max = cell_y_max(y_coords+1)
    n_global_grid_min(2) = ny_global_grid_min
    n_global_grid_max(2) = ny_global_grid_max

    nz_global_grid_min = cell_z_min(z_coords+1)
    nz_global_grid_max = cell_z_max(z_coords+1)
    n_global_grid_min(3) = nz_global_grid_min
    n_global_grid_max(3) = nz_global_grid_max

    nx_grid = nx_global_grid_max - nx_global_grid_min + 1
    ny_grid = ny_global_grid_max - ny_global_grid_min + 1
    nz_grid = nz_global_grid_max - nz_global_grid_min + 1

    !!! --- number of cells of each subdomain
    nx=nx_grid-1
    ny=ny_grid-1
    nz=nz_grid-1

    ALLOCATE(x(-nxguards:nx+nxguards), y(-nyguards:ny+nyguards), z(-nzguards:nz+nzguards))
    ALLOCATE(x_global(-nxguards:nx_global+nxguards))
    ALLOCATE(y_global(-nyguards:ny_global+nyguards))
    ALLOCATE(z_global(-nzguards:nz_global+nzguards))

    !!! -- sets xmax, ymax, zmax
    xmax = nx_global*dx
    ymax = ny_global*dy
    zmax = nz_global*dz

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
    DO iproc = 0, nprocx-1
        x_grid_mins(iproc) = x_global(cell_x_min(iproc+1)-1)
        x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocy-1
        y_grid_mins(iproc) = y_global(cell_y_min(iproc+1)-1)
        y_grid_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocz-1
        z_grid_mins(iproc) = z_global(cell_z_min(iproc+1)-1)
        z_grid_maxs(iproc) = z_global(cell_z_max(iproc+1))
    ENDDO

    x_min_local = x_grid_mins(x_coords)
    x_max_local = x_grid_maxs(x_coords)
    y_min_local = y_grid_mins(y_coords)
    y_max_local = y_grid_maxs(y_coords)
    z_min_local = z_grid_mins(z_coords)
    z_max_local = z_grid_maxs(z_coords)


    ! --- Allocate grid quantities
    ALLOCATE(ex(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(ey(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(ez(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(exsm(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(eysm(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(ezsm(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bxsm(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bysm(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bzsm(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bx(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(by(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bz(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(jx(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(jy(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(jz(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(rho(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(dive(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))

    start_time = MPI_WTIME()

  END SUBROUTINE mpi_initialise



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


END MODULE mpi_routines

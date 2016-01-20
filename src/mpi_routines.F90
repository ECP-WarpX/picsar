MODULE mpi_routines
  USE shared_data
  USE fields
  USE mpi
  IMPLICIT NONE
  !PRIVATE
  !PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(num) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init()
    LOGICAL(isp) :: isinitialized
	INTEGER(isp) :: nproc_comm, rank_in_comm
    CALL MPI_INITIALIZED(isinitialized,errcode)
    IF (.NOT. isinitialized) CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,provided,errcode)
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
	nproc=INT(nproc_comm,idp)
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
	rank=INT(rank_in_comm,idp)
  END SUBROUTINE mpi_minimal_init



  SUBROUTINE setup_communicator

    INTEGER(isp), PARAMETER :: ndims = 3
    INTEGER(idp) :: idim
	INTEGER(isp) :: nproc_comm,dims(ndims), old_comm, ierr
    LOGICAL(isp) :: periods(ndims), reorder, op, reset
    INTEGER(isp) :: test_coords(ndims), rank_in_comm
    INTEGER(idp) :: ix, iy, iz
    INTEGER(idp) :: nxsplit, nysplit, nzsplit
    INTEGER(isp) :: ranges(3,1), nproc_orig, oldgroup, newgroup
    CHARACTER(LEN=11) :: str

    nx_global=nx_global_grid-1
    ny_global=ny_global_grid-1
    nz_global=nz_global_grid-1

    !!! --- NB: CPU Split performed on number of grid points (not cells)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc_comm, ierr)
	nproc=INT(nproc_comm,idp)

    dims = (/nprocz, nprocy, nprocx/)

    ! Initial CPU split sanity check
    IF ((nprocx .EQ. 0) .AND. (nprocy .EQ. 0) .AND. (nprocz .EQ. 0)) THEN
        CALL MPI_DIMS_CREATE(nproc_comm, ndims, dims, errcode)
        nprocx = INT(dims(3),idp)
        nprocy = INT(dims(2),idp)
        nprocz = INT(dims(1),idp)
    ENDIF

    IF (nproc .NE. nprocx*nprocy*nprocz) THEN
        IF (rank .EQ. 0) THEN
            WRITE(0,*) '*** ERROR ***'
            WRITE(0,*) 'nprocx*nprocy*nprocz =/ # of MPI processes'
            WRITE(0,*) ' Check input file '
            CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
    ENDIF

    IF (nx_global_grid .LT. nxguards .OR. ny_global_grid .LT. nyguards .OR. nz_global_grid .LT. nzguards) THEN
      IF (rank .EQ. 0) THEN
        WRITE(0,*) '*** ERROR ***'
        WRITE(0,*) 'Simulation domain is too small.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

    IF (nprocx * nprocy * nprocz .GT. 0) THEN
      nxsplit = nx_global_grid / nprocx
      nysplit = ny_global_grid / nprocy
      nzsplit = nz_global_grid / nprocz
      IF (nxsplit .LT. nxguards .OR. nysplit .LT. nyguards .OR. nzsplit .LT. nzguards) THEN
        IF (rank .EQ. 0) THEN
            WRITE(0,*) 'WRONG CPU SPLIT nlocal<nguards'
            CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
      ENDIF
    ENDIF

    periods = .FALSE.
    reorder = .TRUE.

    ! Set boundary to be periodic in x,y,z for particles and fields

    periods(c_ndims) = .TRUE.
    periods(c_ndims-1) = .TRUE.
    periods(c_ndims-2) = .TRUE.

    old_comm = comm
    CALL MPI_CART_CREATE(old_comm, ndims, dims, periods, reorder, comm, errcode)
    CALL MPI_COMM_FREE(old_comm, errcode)
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
    CALL MPI_CART_COORDS(comm, rank_in_comm, ndims, coordinates, errcode)
	rank=INT(rank_in_comm,idp)
    CALL MPI_CART_SHIFT(comm, 2_isp, 1_isp, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 1_isp, 1_isp, proc_y_min, proc_y_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0_isp, 1_isp, proc_z_min, proc_z_max, errcode)

    nprocdir = dims

    IF (rank .EQ. 0) THEN
      WRITE(0,*) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
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

    INTEGER(isp) :: idim
    INTEGER(isp) :: nx0, nxp, nx_tot_grid
    INTEGER(isp) :: ny0, nyp, ny_tot_grid
    INTEGER(isp) :: nz0, nzp, nz_tot_grid
    INTEGER(isp) :: iproc, ix, iy, iz

    ! Init number of guard cells of subdomains in each dimension
	! Init in Python 
    !nxguards = MAX(nox,norderx)+npass(1)
    !nyguards = MAX(noy,nordery)+npass(2)
    !nzguards = MAX(noz,norderz)+npass(3)
    !nxjguards = MAX(nox,2)
    !nyjguards = MAX(noy,2)
    !nzjguards = MAX(noz,2)

    IF (l_smooth_compensate) THEN
        nxguards = nxguards + 1
        nyguards = nyguards + 1
        nzguards = nzguards + 1
    END IF

    CALL setup_communicator

	ALLOCATE(x_grid_mins(0:nprocx-1), x_grid_maxs(0:nprocx-1))
    ALLOCATE(y_grid_mins(0:nprocy-1), y_grid_maxs(0:nprocy-1))
    ALLOCATE(z_grid_mins(0:nprocz-1), z_grid_maxs(0:nprocz-1))
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))
    ALLOCATE(cell_y_min(nprocy), cell_y_max(nprocy))
    ALLOCATE(cell_z_min(nprocz), cell_z_max(nprocz))

	! Split is done on the total number of grid nodes
	! Initial WARP split is used with each processor boundary 
	! being shared by two adjacent MPI processes 
	nx_tot_grid = nx_global_grid+nprocx-1 ! Total number of grid points along x
	ny_tot_grid = ny_global_grid+nprocy-1 ! Total number of grid points along y
	nz_tot_grid = nz_global_grid+nprocz-1 ! Total number of grid points along z

    nx0 = nx_tot_grid / nprocx
    ny0 = ny_tot_grid / nprocy
    nz0 = nz_tot_grid / nprocz

    ! If the total number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx .NE. nx_tot_grid) THEN
        nxp = (nx0 + 1) * nprocx - nx_tot_grid
    ELSE
        nxp = nprocx
    ENDIF

    IF (ny0 * nprocy .NE. ny_tot_grid) THEN
        nyp = (ny0 + 1) * nprocy - ny_tot_grid
    ELSE
        nyp = nprocy
    ENDIF

    IF (nz0 * nprocz .NE. nz_tot_grid) THEN
        nzp = (nz0 + 1) * nprocz - nz_tot_grid
    ELSE
        nzp = nprocz
    ENDIF

	cell_x_min(1)=1
	cell_x_max(1)=1+nx0-1
    DO idim = 2, nxp
        cell_x_min(idim) = cell_x_max(idim-1)
        cell_x_max(idim) = cell_x_min(idim)+nx0-1
    ENDDO
    DO idim = nxp + 1, nprocx
        cell_x_min(idim) = cell_x_max(idim-1)
        cell_x_max(idim) = cell_x_min(idim)+nx0
    ENDDO

	cell_y_min(1)=1
	cell_y_max(1)=1+ny0-1
    DO idim = 2, nyp
        cell_y_min(idim) = cell_y_max(idim-1)
        cell_y_max(idim) = cell_y_min(idim)+ny0-1
    ENDDO
    DO idim = nyp + 1, nprocy
        cell_y_min(idim) = cell_y_max(idim-1)
        cell_y_max(idim) = cell_y_min(idim)+ny0
    ENDDO

	cell_z_min(1)=1
	cell_z_max(1)=1+nz0-1
    DO idim = 2, nzp
        cell_z_min(idim) = cell_z_max(idim-1)
        cell_z_max(idim) = cell_z_min(idim)+nz0-1
    ENDDO
    DO idim = nzp + 1, nprocz
        cell_z_min(idim) = cell_z_max(idim-1)
        cell_z_max(idim) = cell_z_min(idim)+nz0
    ENDDO


    nx_global_grid_min = cell_x_min(x_coords+1)
    nx_global_grid_max = cell_x_max(x_coords+1)


    ny_global_grid_min = cell_y_min(y_coords+1)
    ny_global_grid_max = cell_y_max(y_coords+1)

    nz_global_grid_min = cell_z_min(z_coords+1)
    nz_global_grid_max = cell_z_max(z_coords+1)


    !!! --- number of gridpoints of each subdomain
    nx_grid = nx_global_grid_max - nx_global_grid_min + 1
    ny_grid = ny_global_grid_max - ny_global_grid_min + 1
    nz_grid = nz_global_grid_max - nz_global_grid_min + 1

    !!! --- number of cells of each subdomain
    nx=nx_grid-1
    ny=ny_grid-1
    nz=nz_grid-1

	!print *,"X split: RANK, nx_global_grid, nx_grid, nx_global_grid_min, nx_global_grid_maxx",rank,nx_global_grid,nx_grid,nx_global_grid_min,nx_global_grid_max
	!print *, "Y split: RANK, ny_global_grid, ny_grid, ny_global_grid_min, ny_global_grid_maxx",rank,ny_global_grid,ny_grid,ny_global_grid_min,ny_global_grid_max
	!print *, "Z split: RANK, nz_global_grid, nz_grid, nz_global_grid_min, nz_global_grid_maxx",rank,nz_global_grid,nz_grid,nz_global_grid_min,nz_global_grid_max

	! Allocate arrays of axis
    ALLOCATE(x(-nxguards:nx+nxguards), y(-nyguards:ny+nyguards), z(-nzguards:nz+nzguards))
    ALLOCATE(x_global(-nxguards:nx_global+nxguards))
    ALLOCATE(y_global(-nyguards:ny_global+nyguards))
    ALLOCATE(z_global(-nzguards:nz_global+nzguards))

    !!! -- sets xmax, ymax, zmax
    !xmax = nx_global*dx
    !ymax = ny_global*dy
    !zmax = nz_global*dz

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
        x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1)-1)
    ENDDO
    DO iproc = 0, nprocy-1
        y_grid_mins(iproc) = y_global(cell_y_min(iproc+1)-1)
        y_grid_maxs(iproc) = y_global(cell_y_max(iproc+1)-1)
    ENDDO
    DO iproc = 0, nprocz-1
        z_grid_mins(iproc) = z_global(cell_z_min(iproc+1)-1)
        z_grid_maxs(iproc) = z_global(cell_z_max(iproc+1)-1)
    ENDDO

    x_min_local = x_grid_mins(x_coords)
    x_max_local = x_grid_maxs(x_coords)
    y_min_local = y_grid_mins(y_coords)
    y_max_local = y_grid_maxs(y_coords)
    z_min_local = z_grid_mins(z_coords)
    z_max_local = z_grid_maxs(z_coords)

    x_grid_min_local=x_min_local
    y_grid_min_local=y_min_local
    z_grid_min_local=z_min_local
    x_grid_max_local=x_max_local
    y_grid_max_local=y_max_local
    z_grid_max_local=z_max_local


	!print *, "X split: RANK, xmin, xmax, dx, x_min_local, x_max_local", rank, xmin, xmax, dx, x_min_local, x_max_local
	!print *, "Y split: RANK, ymin, ymax, dy, y_min_local, y_max_local", rank, ymin, ymax, dy, y_min_local, y_max_local
	!print *, "Z split: RANK, zmin, zmax, dz, z_min_local, z_max_local", rank, zmin, zmax, dz, z_min_local, z_max_local

    ! --- Allocate grid quantities
	PRINT *, "nxjguards,nyjguards,nzjguards", nxjguards,nyjguards,nzjguards
	PRINT *, "nxguards,nyguards,nzguards", nxguards,nyguards,nzguards
    ALLOCATE(ex(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(ey(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(ez(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bx(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(by(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(bz(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
    ALLOCATE(jx(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
    ALLOCATE(jy(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
    ALLOCATE(jz(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
    ALLOCATE(rho(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
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

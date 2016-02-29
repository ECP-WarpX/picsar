MODULE mpi_routines
  USE shared_data
  USE fields
  USE mpi
  IMPLICIT NONE
  !PRIVATE
  !PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(num) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init(comm_in)
    LOGICAL(isp) :: isinitialized
	INTEGER(isp) :: nproc_comm, rank_in_comm
	INTEGER(idp), OPTIONAL, INTENT(IN) :: comm_in 
	
    CALL MPI_INITIALIZED(isinitialized,errcode)
    IF (.NOT. isinitialized) CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,provided,errcode)
    IF (present(comm_in)) THEN 
    	CALL MPI_COMM_DUP(INT(comm_in,isp), comm, errcode)    	
    ELSE
    	CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    ENDIF
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
	nproc=INT(nproc_comm,idp)
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
	rank=INT(rank_in_comm,idp)
  END SUBROUTINE mpi_minimal_init



  SUBROUTINE setup_communicator

    INTEGER(isp), PARAMETER :: ndims = 3
    INTEGER(idp) :: idim
	INTEGER(isp) :: nproc_comm,dims(ndims), old_comm, ierr, neighb
	INTEGER(isp) :: proc_x_minsp, proc_x_maxsp,proc_y_minsp, proc_y_maxsp, & 
					proc_z_minsp, proc_z_maxsp
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
    CALL MPI_CART_SHIFT(comm, 2_isp, 1_isp, proc_x_minsp, proc_x_maxsp, errcode)
    CALL MPI_CART_SHIFT(comm, 1_isp, 1_isp, proc_y_minsp, proc_y_maxsp, errcode)
    CALL MPI_CART_SHIFT(comm, 0_isp, 1_isp, proc_z_minsp, proc_z_maxsp, errcode)
    
	proc_x_min=INT(proc_x_minsp,idp)
	proc_x_max=INT(proc_x_maxsp,idp)
	proc_y_min=INT(proc_y_minsp,idp)
	proc_y_max=INT(proc_y_maxsp,idp)
	proc_z_min=INT(proc_z_minsp,idp)
	proc_z_max=INT(proc_z_maxsp,idp)
	
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
            CALL MPI_CART_RANK(comm, test_coords, neighb, errcode)
            neighbour(ix,iy,iz)=INT(neighb,idp)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER(isp) :: idim
    INTEGER(isp) :: nx0, nxp
    INTEGER(isp) :: ny0, nyp
    INTEGER(isp) :: nz0, nzp
    INTEGER(isp) :: iproc, ix, iy, iz

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

	! Split is done on the total number of cells as in WARP 
	! Initial WARP split is used with each processor boundary 
	! being shared by two adjacent MPI processes 
    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nz0 = nz_global / nprocz

    ! If the total number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 cells
    ! The remaining processors have nx0+1 cells
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
    DO idim = nzp+1 , nprocz
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

	! Allocate arrays of axis
    ALLOCATE(x(-nxguards:nx+nxguards), y(-nyguards:ny+nyguards), z(-nzguards:nz+nzguards))
    ALLOCATE(x_global(-nxguards:nx_global+nxguards))
    ALLOCATE(y_global(-nyguards:ny_global+nyguards))
    ALLOCATE(z_global(-nzguards:nz_global+nzguards))

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

  END SUBROUTINE mpi_initialise

  SUBROUTINE allocate_grid_quantities()
      ! --- Allocate grid quantities
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
    ! --- Quantities used by the dynamic load balancer 
    ALLOCATE(new_cell_x_min(1:nprocx), new_cell_x_max(1:nprocx))
    ALLOCATE(new_cell_y_min(1:nprocy), new_cell_y_max(1:nprocy))
    ALLOCATE(new_cell_z_min(1:nprocz), new_cell_z_max(1:nprocz))
  END SUBROUTINE 
  
  
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

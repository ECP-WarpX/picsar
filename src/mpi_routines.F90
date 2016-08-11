! ______________________________________________________________________________
! 
! MPI_ROUTINES.F90
! Contains subroutines for MPI
! 

MODULE mpi_routines
  USE shared_data
  USE fields
  USE mpi
  USE params
  IMPLICIT NONE
  !PRIVATE
  !PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(num) :: start_time, end_time

CONTAINS

SUBROUTINE mpi_minimal_init()
    LOGICAL(isp) :: isinitialized
	  INTEGER(isp) :: nproc_comm, rank_in_comm
	
	  !print*,'start mpi_minimal_init()'
	  !print*,'MPI_INITIALIZED'
    CALL MPI_INITIALIZED(isinitialized,errcode)
    IF (.NOT. isinitialized) THEN
      !print*, 'MPI_INIT_THREAD'
      CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,provided,errcode)
      !CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,provided,errcode)
    ENDIF
    !print*,'MPI_COMM_DUP'
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    !print*,'MPI_COMM_SIZE'
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
    nproc=INT(nproc_comm,idp)
    !print*,'MPI_COMM_RANK'
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
	  rank=INT(rank_in_comm,idp)	  	  
	  !print*, 'end mpi_minimal_init'
  END SUBROUTINE mpi_minimal_init

  SUBROUTINE mpi_minimal_init_python(comm_in)
    LOGICAL(isp) :: isinitialized
	INTEGER(isp) :: nproc_comm, rank_in_comm
	INTEGER(idp), OPTIONAL, INTENT(IN) :: comm_in 
	
    CALL MPI_INITIALIZED(isinitialized,errcode)
    IF (.NOT. isinitialized) CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,provided,errcode)
    IF (present(comm_in) .AND. comm_in .GE. 0) THEN 
    	CALL MPI_COMM_DUP(INT(comm_in,isp), comm, errcode)    	
    ELSE
    	CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    ENDIF
    CALL MPI_COMM_SIZE(comm, nproc_comm, errcode)
    nproc=INT(nproc_comm,idp)
    
    CALL MPI_COMM_RANK(comm, rank_in_comm, errcode)
	rank=INT(rank_in_comm,idp)	  	  

  END SUBROUTINE mpi_minimal_init_python

  SUBROUTINE setup_communicator

    INTEGER(isp), PARAMETER :: ndims = 3
    INTEGER(idp) :: idim
	INTEGER(isp) :: nproc_comm,dims(ndims), old_comm, ierr, neighb
	INTEGER(isp) :: proc_x_minsp, proc_x_maxsp,proc_y_minsp, proc_y_maxsp, & 
					proc_z_minsp, proc_z_maxsp
    LOGICAL(isp) :: periods(ndims), reorder, op, reset
    INTEGER(isp) :: test_coords(ndims), rank_in_comm
    INTEGER(idp) :: ix, iy, iz
    INTEGER(idp) :: x_coords_neight, y_coords_neight, z_coords_neight
    INTEGER(idp) :: nxsplit, nysplit, nzsplit
    INTEGER(isp) :: ranges(3,1), nproc_orig, oldgroup, newgroup
    INTEGER(idp) :: rankyz
    INTEGER(idp) :: new_rank,old_rank
    INTEGER(idp), dimension(:,:,:),allocatable :: topo_array
    REAL(num) :: r
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
            WRITE(0,*) ' Total number of processors:',nproc
            WRITE(0,*) ' Number of processors in each direction:',nprocx,nprocy,nprocz
            CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
    ENDIF

    IF (nx_global_grid .LT. nxguards .OR. ny_global_grid .LT. nyguards .OR. nz_global_grid .LT. nzguards) THEN
      IF (rank .EQ. 0) THEN
        WRITE(0,*) '*** ERROR ***'
        WRITE(0,*) 'Simulation domain is too small.'
        WRITE(0,*) 'nx_global_grid',nx_global_grid,'nxguards',nxguards
        WRITE(0,*) ny_global_grid,nyguards
        WRITE(0,*) nz_global_grid,nzguards
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
            WRITE(0,*) nxsplit,nysplit,nzsplit
            WRITE(0,*) nx_global_grid,ny_global_grid,nz_global_grid
            CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
      ENDIF
    ENDIF

    periods = .FALSE.
    reorder = .TRUE.
	
    ! Random topology ----------------------------------------
    IF (topology == 2) THEN

      IF (rank .EQ. 0) THEN
        WRITE(0,*) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
      ENDIF

      allocate(topo_array(0:nprocx-1,0:nprocy-1,0:nprocz-1))
    
      IF (rank.EQ.0) THEN
      
        ! Creation of a random array for the topology
        DO iz=0,nprocz-1
          DO iy=0,nprocy-1
            DO ix=0,nprocx-1
              topo_array(ix,iy,iz) = INT(ix + (iy)*nprocx + (iz)*nprocy*nprocx,idp) 
              !print*, ix + (iy)*nprocx + (iz)*nprocy*nprocx
            ENDDO
          ENDDO
        ENDDO            
        DO iz =0,nprocz-1
          DO iy=0,nprocy-1
            DO ix=0,nprocx-1
              CALL RANDOM_NUMBER(r)
              x_coords = INT(r*(nprocx-1))
              ! x_coords = INT(RAN(seed)*(nprocx-1))
              CALL RANDOM_NUMBER(r)
              y_coords = INT(r*(nprocy-1))
              CALL RANDOM_NUMBER(r)
              z_coords = INT(r*(nprocz-1))
              rankyz = topo_array(ix,iy,iz)
              topo_array(ix,iy,iz) = topo_array(x_coords,y_coords,z_coords)
              topo_array(x_coords,y_coords,z_coords) = rankyz
            ENDDO
          ENDDO
        ENDDO
      
        !print*, 'rank',rank,'topo_array',topo_array
      
      ENDIF
      
      ! Sharing of the topology
      CALL MPI_BCAST(topo_array,INT(nprocx*nprocy*nprocz,isp),MPI_INTEGER8,0_isp,comm,errcode)
      
      !IF (rank==3) print*, 'rank',rank,'topo_array',topo_array
      
      ! Each processor determine their coordinates
      DO iz =0,nprocz-1
        DO iy=0,nprocy-1
          DO ix=0,nprocx-1      

            IF (rank .EQ. topo_array(ix,iy,iz)) THEN
            
              x_coords = ix
              y_coords = iy
              z_coords = iz
            
            ENDIF     

          ENDDO
        ENDDO
      ENDDO

      x_min_boundary = .FALSE.
      x_max_boundary = .FALSE.
      IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
      IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.
      
      y_min_boundary = .FALSE.
      y_max_boundary = .FALSE.
      IF (y_coords .EQ. 0) y_min_boundary = .TRUE.
      IF (y_coords .EQ. nprocy - 1) y_max_boundary = .TRUE.

      z_min_boundary = .FALSE.
      z_max_boundary = .FALSE.
      IF (z_coords .EQ. 0) z_min_boundary = .TRUE.
      IF (z_coords .EQ. nprocz - 1) z_max_boundary = .TRUE.   

      neighbour = MPI_PROC_NULL     
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
          
            IF (x_coords + ix .LT. 0) THEN            
              x_coords_neight = x_coords + nprocx-1           
            ELSE IF (x_coords + ix .GT. nprocx-1) THEN              
              x_coords_neight = x_coords - (nprocx-1)
            ELSE
              x_coords_neight = x_coords + ix          
            END IF 

            IF (y_coords + iy .LT. 0) THEN            
              y_coords_neight = y_coords + nprocy-1           
            ELSE IF (y_coords + iy .GT. nprocy-1) THEN              
              y_coords_neight = y_coords - (nprocy-1) 
            ELSE
              y_coords_neight = y_coords + iy             
            END IF
            
            IF (z_coords + iz .LT. 0) THEN            
              z_coords_neight = z_coords + nprocz-1           
            ELSE IF (z_coords + iz .GT. nprocz-1) THEN              
              z_coords_neight = z_coords - (nprocz-1)  
            ELSE
              z_coords_neight = z_coords + iz                              
            END IF            
                 
            neighbour(ix,iy,iz) = topo_array(x_coords_neight,y_coords_neight,z_coords_neight)
                                
           !print*,'rank',rank, ix,iy,iz, neighbour(ix,iy,iz) 
            
          ENDDO
        ENDDO
      ENDDO  

      proc_x_max = neighbour(1,0,0)            
      proc_x_min = neighbour(-1,0,0)
      proc_y_max = neighbour(0,1,0)            
      proc_y_min = neighbour(0,-1,0) 
      proc_z_max = neighbour(0,0,1)            
      proc_z_min = neighbour(0,0,-1)  
    
    ! Topology - processors are ordered according to their ranks ---------------
    ELSE IF (topology == 1) THEN

      IF (rank .EQ. 0) THEN
        WRITE(0,*) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
      ENDIF
    
      ! We first fill the x direction, then y and finally z
      x_coords = MOD(rank,nprocx)
      y_coords = MOD((rank-x_coords)/nprocx,nprocy)
      z_coords = (rank-x_coords - y_coords*nprocx)/(nprocx*nprocy)
      
      !print*,'rank:',rank,rankyz,x_coords,y_coords,z_coords
         
      x_min_boundary = .FALSE.
      x_max_boundary = .FALSE.
      IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
      IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.
      
      y_min_boundary = .FALSE.
      y_max_boundary = .FALSE.
      IF (y_coords .EQ. 0) y_min_boundary = .TRUE.
      IF (y_coords .EQ. nprocy - 1) y_max_boundary = .TRUE.

      z_min_boundary = .FALSE.
      z_max_boundary = .FALSE.
      IF (z_coords .EQ. 0) z_min_boundary = .TRUE.
      IF (z_coords .EQ. nprocz - 1) z_max_boundary = .TRUE.      

      neighbour = MPI_PROC_NULL     
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
          
            IF (x_coords + ix .LT. 0) THEN            
              x_coords_neight = x_coords + nprocx-1           
            ELSE IF (x_coords + ix .GT. nprocx-1) THEN              
              x_coords_neight = x_coords - (nprocx-1)
            ELSE
              x_coords_neight = x_coords + ix          
            END IF 

            IF (y_coords + iy .LT. 0) THEN            
              y_coords_neight = y_coords + nprocy-1           
            ELSE IF (y_coords + iy .GT. nprocy-1) THEN              
              y_coords_neight = y_coords - (nprocy-1) 
            ELSE
              y_coords_neight = y_coords + iy             
            END IF
            
            IF (z_coords + iz .LT. 0) THEN            
              z_coords_neight = z_coords + nprocz-1           
            ELSE IF (z_coords + iz .GT. nprocz-1) THEN              
              z_coords_neight = z_coords - (nprocz-1)  
            ELSE
              z_coords_neight = z_coords + iz                              
            END IF            
                 
            neighbour(ix,iy,iz) = x_coords_neight + (nprocx)*(y_coords_neight) &
                                + (nprocx*nprocy)*(z_coords_neight)
                                
           print*,'rank',rank, ix,iy,iz, neighbour(ix,iy,iz) 
            
          ENDDO
        ENDDO
      ENDDO  

      proc_x_max = neighbour(1,0,0)            
      proc_x_min = neighbour(-1,0,0)
      proc_y_max = neighbour(0,1,0)            
      proc_y_min = neighbour(0,-1,0) 
      proc_z_max = neighbour(0,0,1)            
      proc_z_min = neighbour(0,0,-1)     
          
    ! Default topology - Cartesian topology
    ELSE
   	  ! Set boundary to be periodic in x,y,z for particles and fields by default 
   	  periods(c_ndims) = .TRUE.
   	  periods(c_ndims-1) = .TRUE.
   	  periods(c_ndims-2) = .TRUE.
   	  
      ! Create an array of rank corresponding to the current topology
      !ALLOCATE(rank_array(nprocs))
      old_rank = rank
      !old_rank_array(old_rank) = rank
        
      ! Creation of a Cartesian topology by using the original communicator
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

      IF (rank .EQ. 0) THEN
        WRITE(0,*) 'Cartesian topology'
        WRITE(0,*) 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
        WRITE(0,*)
      ENDIF 

      ! ------------------------------------------------------------
      ! Checking of the new topology
      CALL MPI_BARRIER(comm,errcode)
      new_rank = rank     
      !new_rank_array(old_rank) = rank
      !IF (new_rank.NE.old_rank) WRITE(0,'(A,I5,A,I5)') 'Rank switched from ',old_rank,' to ',new_rank
      !WRITE(0,'(X,A,I5,A,I5)') 'Rank switched from ',old_rank,' to ',new_rank
      ! We first fill the x direction, then y and finally z
      x_coords = MOD(rank,nprocx)
      y_coords = MOD((rank-x_coords)/nprocx,nprocy)
      z_coords = (rank-x_coords - y_coords*nprocx)/(nprocx*nprocy)
      !WRITE(0,'(X,A30,3(X,I5))') 'Coordinates from MPI_CART_COORDS:',coordinates(1:3)
      !WRITE(0,'(X,A30,3(X,I5))') 'Theoretical coordinates:',z_coords,y_coords,x_coords
      CALL MPI_BARRIER(comm,errcode)
      ! -------------------------------------------------------------

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
              IF ((test_coords(idim) .LT. 0 &            ! If outside the grid
                  .OR. test_coords(idim) .GE. dims(idim)) & ! If outside the grid
                  .AND. .NOT. periods(idim)) op = .FALSE. ! And periodic or not in this direction
            ENDDO
            IF (op) THEN
              ! Then determine the rank of the neighbor processor
              CALL MPI_CART_RANK(comm, test_coords, neighb, errcode)
              neighbour(ix,iy,iz)=INT(neighb,idp)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
    ENDIF

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


	CALL compute_simulation_axis()

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


  SUBROUTINE compute_simulation_axis()
  	IMPLICIT NONE
  	INTEGER(idp) :: ix, iy, iz, iproc 
  	IF (.NOT. l_axis_allocated) THEN 
  		! Allocate arrays of axis
    	ALLOCATE(x(-nxguards:nx+nxguards), y(-nyguards:ny+nyguards), z(-nzguards:nz+nzguards))
    	ALLOCATE(x_global(-nxguards:nx_global+nxguards))
    	ALLOCATE(y_global(-nyguards:ny_global+nyguards))
    	ALLOCATE(z_global(-nzguards:nz_global+nzguards))
    	l_axis_allocated=.TRUE.
    ENDIF

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
    
  END SUBROUTINE


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

! ________________________________________________________________________________________
  SUBROUTINE time_statistics
!
!  Subroutine dedicated to the time report at the end of the simulation
! ________________________________________________________________________________________    
    USE time_stat
    USE params
    USE omp_lib 
    IMPLICIT NONE
        
    REAL(num), DIMENSION(20) :: mintimes, init_mintimes
    REAL(num), DIMENSION(20) :: maxtimes, init_maxtimes
    REAL(num), DIMENSION(20) :: avetimes, init_avetimes
    REAL(num), DIMENSION(20) :: percenttimes
    INTEGER(idp)             :: nthreads_tot
    
    ! Total times
    localtimes(20) = sum(localtimes(1:14))
    localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + &
                     localtimes(8) + localtimes(11) + localtimes(13)
    localtimes(18) = localtimes(5) + localtimes(6) + localtimes(7) + localtimes(8)                 
    
    init_localtimes(5) = sum(init_localtimes(1:4))
    
    ! Reductions
    ! Maximun times
    CALL MPI_REDUCE(localtimes,maxtimes,20_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)
    CALL MPI_REDUCE(init_localtimes,init_maxtimes,5_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)    
    ! Minimum times
    CALL MPI_REDUCE(localtimes,mintimes,20_isp,mpidbl,MPI_MIN,0_isp,comm,errcode)
    CALL MPI_REDUCE(init_localtimes,init_mintimes,5_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)    
    ! Average
    CALL MPI_REDUCE(localtimes,avetimes,20_isp,mpidbl,MPI_SUM,0_isp,comm,errcode)
    CALL MPI_REDUCE(init_localtimes,init_avetimes,5_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)    
    avetimes = avetimes / nproc
   
    ! Percentage 
    percenttimes = avetimes / avetimes(20) * 100
     
  IF (rank .EQ. 0) THEN
    WRITE(0,*) '___________________________________________________________________________'
    WRITE(0,'(X,A40,X)') "Time statistics Initialization:"
    WRITE(0,*) ""
    WRITE(0,'(X,A25,X,A8,X,A8,X,A8,X,A8,X,A8)') "Step part","min (s)","ave (s)","max (s)"
    WRITE(0,*) "---------------------------------------------------------------------------"
    WRITE(0,'(X,A25,5(X,F8.2))') "Tiling and part. load:", init_mintimes(1), init_avetimes(1), init_maxtimes(1)
    WRITE(0,*) ""
    WRITE(0,*) '___________________________________________________________________________'
    WRITE(0,'(X,A40,X)') "Time statistics main loop:"
    WRITE(0,*) ""
    WRITE(0,'(X,A25,X,A8,X,A8,X,A8,X,A8,X,A8)') "Step part","min (s)","ave (s)","max (s)","per (%)","/it (ms)"
    WRITE(0,*) "---------------------------------------------------------------------------"
    IF (fg_p_pp_seperated.le.1) THEN
      WRITE(0,'(X,A25,5(X,F8.2))') "Particle pusher + field gathering:", mintimes(1), avetimes(1), maxtimes(1),&
      percenttimes(1), avetimes(1)/nsteps*1e3
    ELSE
      WRITE(0,'(X,A25,5(X,F8.2))') "Field gathering:",mintimes(14), avetimes(14), maxtimes(14),&
      percenttimes(14), avetimes(14)/nsteps*1e3    
      WRITE(0,'(X,A25,5(X,F8.2))') "Particle pusher:", mintimes(1), avetimes(1), maxtimes(1),&
      percenttimes(1), avetimes(1)/nsteps*1e3
    ENDIF
    WRITE(0,'(X,A25,5(X,F8.2))') "Particle MPI bound. cond.:", mintimes(2), avetimes(2), maxtimes(2),&
    percenttimes(2), avetimes(2)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "Particle OpenMP bound. cond.:", mintimes(11), avetimes(11), maxtimes(11),&
    percenttimes(11), avetimes(11)/nsteps*1e3    
    WRITE(0,'(X,A25,5(X,F8.2))') "Current deposition:", mintimes(3), avetimes(3), maxtimes(3),&
    percenttimes(3), avetimes(3)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "Current bound. cond.:", mintimes(4), avetimes(4), maxtimes(4),&
    percenttimes(4), avetimes(4)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "Push bfield:", mintimes(5), avetimes(5), maxtimes(5),&
    percenttimes(5), avetimes(5)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "B field bound. cond.:",mintimes(6), avetimes(6), maxtimes(6),&
    percenttimes(6), avetimes(6)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "Push efield:", mintimes(7), avetimes(7), maxtimes(7),&
    percenttimes(7), avetimes(7)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "E field bound. cond.:",mintimes(8), avetimes(8), maxtimes(8),&
    percenttimes(8), avetimes(8)/nsteps*1e3
    WRITE(0,'(X,A25,5(X,F8.2))') "Sorting:",mintimes(10), avetimes(10), maxtimes(10),&
    percenttimes(10), avetimes(10)/nsteps*1e3    
    WRITE(0,'(X,A25,5(X,F8.2))') "Charge deposition:", mintimes(12), avetimes(12), maxtimes(12),&
    percenttimes(12), avetimes(12)/nsteps*1e3     
    WRITE(0,'(X,A25,5(X,F8.2))') "Charge bound. cond.:", mintimes(13), avetimes(13), maxtimes(13),&
    percenttimes(13), avetimes(13)/nsteps*1e3          
    WRITE(0,'(X,A25,5(X,F8.2))') "Diags:",mintimes(9), avetimes(9), maxtimes(9),&
    percenttimes(9), avetimes(9)/nsteps*1e3
    
    WRITE(0,*) ""
    WRITE(0,'(X,A25,5(X,F8.2))') "Total time Maxwell solver:",mintimes(18), avetimes(18), maxtimes(18), &
    percenttimes(18), avetimes(18)/nsteps*1e3    
    WRITE(0,'(X,A25,5(X,F8.2))') "Total time bound. cond.:",mintimes(19), avetimes(19), maxtimes(19), &
    percenttimes(19), avetimes(19)/nsteps*1e3
    WRITE(0,'(X,A25,X,F8.2,X,F8.2,X,F8.2)') "Total time:",mintimes(20), avetimes(20), maxtimes(20)

	
#ifdef _OPENMP
	nthreads_tot=OMP_GET_MAX_THREADS()
	CALL OMP_SET_NESTED(.TRUE.)
#else
	nthreads_tot=1
#endif

    WRITE(0,*) ''
    IF (fg_p_pp_seperated.le.1) THEN    
      WRITE(0,*) 'For lib_performance python class:' 
      WRITE(0,'("(nmpi=",I5,",nomp=",I5,",name='''',kernel=",F6.2,",fieldgave=",F6.2,",part_mpi_com=",F6.2,&
      ",part_omp_com=",F6.2,",currdepo=",F6.2,",currcom=",F6.2,",maxwell=",F6.2,&
      ",maxwellcom=",F6.2,",sorting=",F6.2,",rhodepo=",F6.2,",rhocom=",F6.2,",diags=",F6.2,")")')&
      nproc,nthreads_tot,avetimes(20), &
      avetimes(1),avetimes(2),avetimes(11),avetimes(3),avetimes(4),&
      avetimes(5)+avetimes(7),avetimes(6)+avetimes(8),avetimes(10),&
      avetimes(12),avetimes(13),avetimes(9)
    ELSE
      WRITE(0,*) 'For lib_performance python class:' 
      WRITE(0,'("(nmpi=",I5,",nomp=",I5,",name='''',kernel=",F6.2,",fieldgave=",F6.2,",partpusher=",F6.2,",part_mpi_com=",F6.2,&
      ",part_omp_com=",F6.2,",currdepo=",F6.2,",currcom=",F6.2,",maxwell=",F6.2,&
      ",maxwellcom=",F6.2,",sorting=",F6.2,",rhodepo=",F6.2,",rhocom=",F6.2,",diags=",F6.2,")")') &
      nproc,nthreads_tot,avetimes(20), &
      avetimes(14),avetimes(1),avetimes(2),avetimes(11),avetimes(3),avetimes(4),&
      avetimes(5)+avetimes(7),avetimes(6)+avetimes(8),avetimes(10),&
      avetimes(12),avetimes(13),avetimes(9)
    ENDIF

  ENDIF    
  
  END SUBROUTINE

! ________________________________________________________________________________________
  SUBROUTINE time_statistics_per_iteration
!
!  Subroutine dedicated to the time statistics for one iteration
! ________________________________________________________________________________________    

    USE time_stat
    USE params
    USE omp_lib 
    IMPLICIT NONE
        
    REAL(num), DIMENSION(20) :: mintimes, init_mintimes
    REAL(num), DIMENSION(20) :: maxtimes, init_maxtimes
    REAL(num), DIMENSION(20) :: avetimes, init_avetimes
    REAL(num), DIMENSION(20) :: percenttimes
    INTEGER(idp)             :: nthreads_tot

    ! Time stats per iteration activated
    IF (timestat_perit.gt.0) THEN

      ! Total times
      localtimes(20) = sum(localtimes(1:14))
      localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + &
                       localtimes(8) + localtimes(11) + localtimes(13)
      localtimes(18) = localtimes(5) + localtimes(6) + localtimes(7) + localtimes(8)                 
    
      init_localtimes(5) = sum(init_localtimes(1:4))
    
      ! Reductions
      ! Maximun times
      CALL MPI_REDUCE(localtimes,maxtimes,20_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)
      CALL MPI_REDUCE(init_localtimes,init_maxtimes,5_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)    
      ! Minimum times
      CALL MPI_REDUCE(localtimes,mintimes,20_isp,mpidbl,MPI_MIN,0_isp,comm,errcode)
      CALL MPI_REDUCE(init_localtimes,init_mintimes,5_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)    
      ! Average
      CALL MPI_REDUCE(localtimes,avetimes,20_isp,mpidbl,MPI_SUM,0_isp,comm,errcode)
      CALL MPI_REDUCE(init_localtimes,init_avetimes,5_isp,mpidbl,MPI_MAX,0_isp,comm,errcode)    
      avetimes = avetimes / nproc
   
      ! Percentage 
      percenttimes = avetimes / avetimes(20) * 100

      IF (rank .EQ. 0) THEN
        WRITE(0,*) '___________________________________________________________________________'
        WRITE(0,'(X,A25,X,A8,X,A8,X,A8,X,A8,X,A8)') "Step part","min (s)","ave (s)","max (s)","per (%)","/it (ms)"
        WRITE(0,*) "---------------------------------------------------------------------------"
        IF (fg_p_pp_seperated.le.1) THEN
          WRITE(0,'(X,A25,5(X,F8.2))') "Particle pusher + field gathering:", mintimes(1), avetimes(1), maxtimes(1),&
          percenttimes(1), avetimes(1)/nsteps*1e3
        ELSE
          WRITE(0,'(X,A25,5(X,F8.2))') "Field gathering:",mintimes(14), avetimes(14), maxtimes(14),&
          percenttimes(14), avetimes(14)/nsteps*1e3    
          WRITE(0,'(X,A25,5(X,F8.2))') "Particle pusher:", mintimes(1), avetimes(1), maxtimes(1),&
          percenttimes(1), avetimes(1)/nsteps*1e3
        ENDIF
        WRITE(0,'(X,A25,5(X,F8.2))') "Particle MPI bound. cond.:", mintimes(2), avetimes(2), maxtimes(2),&
        percenttimes(2), avetimes(2)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "Particle OpenMP bound. cond.:", mintimes(11), avetimes(11), maxtimes(11),&
        percenttimes(11), avetimes(11)/nsteps*1e3    
        WRITE(0,'(X,A25,5(X,F8.2))') "Current deposition:", mintimes(3), avetimes(3), maxtimes(3),&
        percenttimes(3), avetimes(3)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "Current bound. cond.:", mintimes(4), avetimes(4), maxtimes(4),&
        percenttimes(4), avetimes(4)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "Push bfield:", mintimes(5), avetimes(5), maxtimes(5),&
        percenttimes(5), avetimes(5)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "B field bound. cond.:",mintimes(6), avetimes(6), maxtimes(6),&
        percenttimes(6), avetimes(6)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "Push efield:", mintimes(7), avetimes(7), maxtimes(7),&
        percenttimes(7), avetimes(7)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "E field bound. cond.:",mintimes(8), avetimes(8), maxtimes(8),&
        percenttimes(8), avetimes(8)/nsteps*1e3
        WRITE(0,'(X,A25,5(X,F8.2))') "Sorting:",mintimes(10), avetimes(10), maxtimes(10),&
        percenttimes(10), avetimes(10)/nsteps*1e3    
        WRITE(0,'(X,A25,5(X,F8.2))') "Charge deposition:", mintimes(12), avetimes(12), maxtimes(12),&
        percenttimes(12), avetimes(12)/nsteps*1e3     
        WRITE(0,'(X,A25,5(X,F8.2))') "Charge bound. cond.:", mintimes(13), avetimes(13), maxtimes(13),&
        percenttimes(13), avetimes(13)/nsteps*1e3          
        WRITE(0,'(X,A25,5(X,F8.2))') "Diags:",mintimes(9), avetimes(9), maxtimes(9),&
        percenttimes(9), avetimes(9)/nsteps*1e3
    
        WRITE(0,*) ""
        WRITE(0,'(X,A25,5(X,F8.2))') "Total time Maxwell solver:",mintimes(18), avetimes(18), maxtimes(18), &
        percenttimes(18), avetimes(18)/nsteps*1e3    
        WRITE(0,'(X,A25,5(X,F8.2))') "Total time bound. cond.:",mintimes(19), avetimes(19), maxtimes(19), &
        percenttimes(19), avetimes(19)/nsteps*1e3
        WRITE(0,'(X,A25,X,F8.2,X,F8.2,X,F8.2)') "Total time:",mintimes(20), avetimes(20), maxtimes(20)
      ENDIF
      
    ENDIF

  END SUBROUTINE

END MODULE mpi_routines

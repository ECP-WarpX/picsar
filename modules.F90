!===============================================================================
! Contains shared data
!===============================================================================
MODULE constants
!===============================================================================
INTEGER, PARAMETER :: num = KIND(1.d0)
REAL(num), parameter :: emass   = 9.1093818800000006e-31_num,      &
                        pmass   = 1.6726231000000001e-27_num,      &
                        echarge = 1.6021764620000001e-19_num,      &
                        clight  = 299792458.0_num,                 &
                        mu0     = 1.2566370614359173e-06_num,      &
                        eps0    = 8.8541878176203892e-12_num,      &
                        pi      = 3.141592653589793_num
INTEGER, PARAMETER :: c_ndims = 3
! direction parameters
INTEGER, PARAMETER :: c_dir_x = 1
INTEGER, PARAMETER :: c_dir_y = 2
INTEGER, PARAMETER :: c_dir_z = 3
LOGICAL:: l_smooth_compensate
INTEGER, PARAMETER :: string_length = 264
! Error handling
INTEGER, PARAMETER :: c_err_bad_value = 2**4
END MODULE constants

!===============================================================================
MODULE fields
!===============================================================================
USE constants
LOGICAL:: l_lower_order_in_v, l_nodalgrid
INTEGER:: nxs=0, nys=0, nzs=0
INTEGER:: norderx, nordery, norderz
INTEGER:: nxguards,nyguards, nzguards, nox, noy, noz, npass(3)
REAL(num):: alpha(3)
REAL(num), POINTER, DIMENSION(:,:,:) :: ex,ey,ez,bx,by,bz,jx,jy,jz,xx,yy,zz, &
                                           exsm,eysm,ezsm,bxsm,bysm,bzsm
REAL(num), POINTER, DIMENSION(:) :: xcoeffs, ycoeffs, zcoeffs ! Fonberg coefficients
END MODULE fields

!===============================================================================
MODULE particles
!===============================================================================
USE constants
LOGICAL :: l_initongrid = .FALSE.
LOGICAL :: l_particles_weight = .FALSE.
LOGICAL :: l4symtry = .FALSE.
INTEGER :: nspecies
INTEGER, PARAMETER :: nspecies_max=4 ! Max number of particle species
REAL(num) :: fdxrand=0.0_num,fdzrand=0.0_num,vthx=0.0_num,vthy=0.0_num,vthz=0.0_num
LOGICAL :: l_species_allocated=.FALSE.
! # of particle tiles in each dimension
INTEGER :: ntilex, ntiley, ntilez
! Fortran object representing a particle tile
TYPE particle_tile
    LOGICAL :: l_arrays_allocated= .FALSE.
    ! Current number of particles in tile
    INTEGER :: np_tile, npmax_tile
    INTEGER :: nx_grid_tile, ny_grid_tile, nz_grid_tile
    INTEGER :: nx_cells_tile, ny_cells_tile, nz_cells_tile
    INTEGER :: nx_tile_min, nx_tile_max, ny_tile_min, ny_tile_max, &
               nz_tile_min, nz_tile_max
    ! Tile position
    REAL(num) :: x_tile_min, y_tile_min, z_tile_min
    REAL(num) :: x_tile_max, y_tile_max, z_tile_max
    ! Local grid quantities in the tile
    REAL(num), POINTER, DIMENSION(:,:,:) :: jx_tile, jy_tile, jz_tile, rho_tile
    REAL(num), POINTER, DIMENSION(:,:,:) :: ex_tile, ey_tile, ez_tile
    REAL(num), POINTER, DIMENSION(:,:,:) :: bx_tile, by_tile, bz_tile
    ! Particle arrays
    REAL(num), POINTER, DIMENSION(:) :: part_x, part_y, part_z
    REAL(num), POINTER, DIMENSION(:) :: part_ux, part_uy, part_uz
    REAL(num), POINTER, DIMENSION(:) :: part_ex, part_ey, part_ez
    REAL(num), POINTER, DIMENSION(:) :: part_bx, part_by, part_bz
    REAL(num), POINTER, DIMENSION(:) :: weight
END TYPE

! Fortran object representing a particle species
TYPE particle_species
    ! Attributes of particle species object
    CHARACTER(LEN=string_length) :: name
    REAL(num) :: charge
    REAL(num) :: mass
    REAL(num) :: x_min
    REAL(num) :: x_max
    REAL(num) :: y_min
    REAL(num) :: y_max
    REAL(num) :: z_min
    REAL(num) :: z_max
    REAL(num) :: vdrift_x
    REAL(num) :: vdrift_y
    REAL(num) :: vdrift_z
    REAL(num) :: vth_x
    REAL(num) :: vth_y
    REAL(num) :: vth_z
    INTEGER   :: species_npart
    INTEGER   :: nppspecies_max
    INTEGER   :: nppcell
    LOGICAL   :: l_arrayoftiles_allocated =.FALSE.
    ! For some stupid reason, cannot use ALLOCATABLE in derived types
    ! in Fortran 90 - Need to use POINTER instead
    TYPE(particle_tile), DIMENSION(:,:,:), POINTER :: array_of_tiles    
END TYPE
! Array of pointers to particle species objects
TYPE(particle_species), POINTER, DIMENSION(:):: species_parray
END MODULE particles

!===============================================================================
MODULE params
!===============================================================================
USE constants
INTEGER :: it,nsteps
REAL(num) :: g0,b0,dt,w0,dtcoef,tmax
REAL(num) :: theta,nlab,wlab,nc,w0_l,w0_t
LOGICAL :: l_coeffs_allocated= .FALSE., l_ck=.FALSE.
REAL(num), PARAMETER :: resize_factor=1.5_num
END MODULE params

!===============================================================================
MODULE shared_data
!===============================================================================
USE mpi
USE constants

!---------------------------------------------------------------------------
! CUSTOM DATATYPES
!---------------------------------------------------------------------------
TYPE aofreals
    REAL(num), POINTER, DIMENSION(:) :: array
END TYPE

TYPE aofint
    INTEGER, POINTER, DIMENSION(:) :: array
END TYPE
!----------------------------------------------------------------------------
! MPI data
!----------------------------------------------------------------------------
INTEGER  :: mpidbl = MPI_DOUBLE_PRECISION
INTEGER :: coordinates(c_ndims), neighbour(-1:1, -1:1, -1:1)
INTEGER :: x_coords, proc_x_min, proc_x_max
INTEGER :: y_coords, proc_y_min, proc_y_max
INTEGER :: z_coords, proc_z_min, proc_z_max
INTEGER :: errcode, comm, tag, rank
INTEGER :: nproc, nprocx, nprocy, nprocz
INTEGER :: nprocdir(c_ndims)
INTEGER :: status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank, ny_each_rank, nz_each_rank
INTEGER, ALLOCATABLE, DIMENSION(:) :: npart_each_rank
LOGICAL :: x_min_boundary, x_max_boundary
LOGICAL :: y_min_boundary, y_max_boundary
LOGICAL :: z_min_boundary, z_max_boundary
! The location of the processors
INTEGER, DIMENSION(:), ALLOCATABLE :: cell_x_min, cell_x_max
INTEGER, DIMENSION(:), ALLOCATABLE :: cell_y_min, cell_y_max
INTEGER, DIMENSION(:), ALLOCATABLE :: cell_z_min, cell_z_max
INTEGER, DIMENSION(:), ALLOCATABLE :: old_x_max, old_y_max, old_z_max
INTEGER :: nx_global_grid_min, nx_global_grid_max
INTEGER :: ny_global_grid_min, ny_global_grid_max
INTEGER :: nz_global_grid_min, nz_global_grid_max
INTEGER :: n_global_grid_min(c_ndims), n_global_grid_max(c_ndims)
! domain and loadbalancing
LOGICAL :: allow_cpu_reduce = .FALSE.
REAL(num), DIMENSION(:), ALLOCATABLE :: x_global, y_global, z_global
REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global
REAL(num), DIMENSION(:), ALLOCATABLE :: xb_offset_global
REAL(num), DIMENSION(:), ALLOCATABLE :: yb_offset_global
REAL(num), DIMENSION(:), ALLOCATABLE :: zb_offset_global
! domain limits and size
INTEGER  :: nx, ny, nz ! local number of cells
INTEGER  :: nx_grid, ny_grid, nz_grid ! local number of grid points
INTEGER  :: nx_global, ny_global, nz_global ! global number of cells
INTEGER  :: nx_global_grid, ny_global_grid, nz_global_grid ! global number of grid points
REAL(num):: dx, xmin, xmax, length_x
REAL(num):: x_min_local, x_max_local
REAL(num):: dy, ymin, ymax,length_y
REAL(num):: y_min_local, y_max_local
REAL(num):: dz, zmin, zmax,length_z
REAL(num):: z_min_local, z_max_local
! Derived types (MPI exchange)
INTEGER :: derived_type_grid
INTEGER :: derived_subarray_grid
! Axis
REAL(num), pointer, dimension(:) :: x, y, z
REAL(num), DIMENSION(:), ALLOCATABLE :: x_grid_mins, x_grid_maxs
REAL(num), DIMENSION(:), ALLOCATABLE :: y_grid_mins, y_grid_maxs
REAL(num), DIMENSION(:), ALLOCATABLE :: z_grid_mins, z_grid_maxs
REAL(num) ::  x_grid_min, x_grid_max
REAL(num) :: x_grid_min_local, x_grid_max_local
REAL(num) ::  y_grid_min, y_grid_max
REAL(num) :: y_grid_min_local, y_grid_max_local
REAL(num) :: z_grid_min, z_grid_max
REAL(num) :: z_grid_min_local, z_grid_max_local

! Total charge density
REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: rho
! Electric Field divergence
REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: dive


! Simulation time statistics
REAL(num) :: startsim =0.0_num
REAL(num) :: endsim =0.0_num

! output frequency
INTEGER :: output_frequency = -1 !(Default is no output)
INTEGER :: output_step_min = 0
INTEGER :: output_step_max = 0

! output quantity flag (Default=False)
INTEGER :: c_output_ex = 0
INTEGER :: c_output_ey = 0
INTEGER :: c_output_ez = 0
INTEGER :: c_output_bx = 0
INTEGER :: c_output_by = 0
INTEGER :: c_output_bz = 0
INTEGER :: c_output_jx = 0
INTEGER :: c_output_jy = 0
INTEGER :: c_output_jz = 0
INTEGER :: c_output_rho = 0
INTEGER :: c_output_dive = 0

! File names for output dumps
CHARACTER(LEN=string_length) :: fileex   ='ex'
CHARACTER(LEN=string_length) :: fileey   ='ey'
CHARACTER(LEN=string_length) :: fileez   ='ez'
CHARACTER(LEN=string_length) :: filebx   ='bx'
CHARACTER(LEN=string_length) :: fileby   ='by'
CHARACTER(LEN=string_length) :: filebz   ='bz'
CHARACTER(LEN=string_length) :: filejx   ='jx'
CHARACTER(LEN=string_length) :: filejy   ='jy'
CHARACTER(LEN=string_length) :: filejz   ='jz'
CHARACTER(LEN=string_length) :: filedive ='dive'
CHARACTER(LEN=string_length) :: filerho  ='rho'
END MODULE shared_data



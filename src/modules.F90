!===============================================================================
! Contains shared data
!===============================================================================
MODULE constants
!===============================================================================
INTEGER, PARAMETER :: num = 8
INTEGER, PARAMETER :: isp = 4
INTEGER, PARAMETER :: idp = 8
REAL(num), PARAMETER :: emass   = 9.10938291e-31_num,      &
                        pmass   = 1.6726231000000001e-27_num,      &
                        echarge = 1.6021764620000001e-19_num,      &
                        clight  = 2.99792458e8_num,                 &
                        mu0     = 1.2566370614359173e-06_num,      &
                        eps0    = 8.854187817620389e-12_num,      &
                        imu0    = 795774.715459,               &
                        pi      = 3.14159265358979323_num
INTEGER(isp), PARAMETER :: c_ndims = 3
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
INTEGER(idp):: nxs=0, nys=0, nzs=0
INTEGER(idp):: norderx, nordery, norderz
INTEGER(idp):: nxguards,nyguards, nzguards, nox, noy, noz, npass(3)
INTEGER(idp):: nxjguards,nyjguards, nzjguards
REAL(num):: alpha(3)
REAL(num), POINTER, DIMENSION(:,:,:) :: ex,ey,ez,bx,by,bz,jx,jy,jz
! Fonberg coefficients
REAL(num), POINTER, DIMENSION(:) :: xcoeffs, ycoeffs, zcoeffs
END MODULE fields

MODULE grid_tilemodule !#do not parse 
USE constants 
TYPE grid_tile
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile, &
 								bxtile, bytile, bztile,  jxtile, jytile, jztile,rhotile
END TYPE
TYPE(grid_tile), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: aofgrid_tiles
END MODULE grid_tilemodule

! Fortran object representing a particle tile
!===============================================================================
MODULE particle_tilemodule !#do not parse 
!===============================================================================
USE constants
TYPE particle_tile
    LOGICAL :: l_arrays_allocated= .FALSE.
    ! Current number of particles in tile
    INTEGER(idp), DIMENSION(1) :: np_tile
    INTEGER(idp) :: npmax_tile
    INTEGER(idp) :: nxg_tile, nyg_tile, nzg_tile
    INTEGER(idp) :: nx_grid_tile, ny_grid_tile, nz_grid_tile
    INTEGER(idp) :: nx_cells_tile, ny_cells_tile, nz_cells_tile
    INTEGER(idp) :: nx_tile_min, nx_tile_max, ny_tile_min, ny_tile_max, &
               nz_tile_min, nz_tile_max
    ! Tile position
    REAL(num) :: x_tile_min, y_tile_min, z_tile_min
    REAL(num) :: x_tile_max, y_tile_max, z_tile_max
    REAL(num) :: x_grid_tile_min, y_grid_tile_min, z_grid_tile_min
    REAL(num) :: x_grid_tile_max, y_grid_tile_max, z_grid_tile_max
    ! Subdomain border flags
    LOGICAL :: subdomain_bound = .FALSE.
    ! Particle arrays
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_x
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_y
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_z
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ux
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uy
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uz
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_gaminv
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ex
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ey
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ez
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_bx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_by
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_bz
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
END TYPE
END MODULE particle_tilemodule

!===============================================================================
MODULE particle_speciesmodule !#do not parse 
!===============================================================================
use particle_tilemodule
use constants
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
    INTEGER(idp)   :: species_npart
    INTEGER(idp)   :: nppspecies_max
    INTEGER(idp)   :: nppcell
    INTEGER(idp)   :: sorting_period
    LOGICAL(idp)   :: l_arrayoftiles_allocated =.FALSE.
    ! For some stupid reason, cannot use ALLOCATABLE in derived types
    ! in Fortran 90 - Need to use POINTER instead
    TYPE(particle_tile), DIMENSION(:,:,:), ALLOCATABLE :: array_of_tiles
    ! Array indicating if a tile has been reallocated
    ! Used for interfacing WARP and PXR
    INTEGER(idp), DIMENSION(:,:,:), ALLOCATABLE :: are_tiles_reallocated
END TYPE
END MODULE particle_speciesmodule

!===============================================================================
MODULE tile_params
!===============================================================================
! # of particle tiles in each dimension
INTEGER :: ntilex, ntiley, ntilez
END MODULE tile_params


!===============================================================================
MODULE particle_properties
!===============================================================================
USE constants
INTEGER(idp), PARAMETER :: npid=1
INTEGER(idp), PARAMETER :: wpid=1
LOGICAL :: l_initongrid = .FALSE.
LOGICAL :: l_particles_weight = .FALSE.
LOGICAL :: l4symtry = .FALSE.
INTEGER(idp) :: pdistr
INTEGER(idp) :: nspecies
INTEGER(idp) :: ntot ! total number of particles (all species, all subdomains -> useful for stat)
!INTEGER(idp), PARAMETER :: nspecies_max=4 ! Max number of particle species
INTEGER(idp) :: nspecies_max=4 ! Max number of particle species
REAL(num) :: fdxrand=0.0_num,fdzrand=0.0_num,vthx=0.0_num,vthy=0.0_num,vthz=0.0_num
LOGICAL :: l_species_allocated=.FALSE.
END MODULE particle_properties



!===============================================================================
MODULE particles
!===============================================================================
USE constants
USE tile_params
USE particle_tilemodule
USE particle_speciesmodule
USE particle_properties
USE grid_tilemodule

! Array of  particle species objects
TYPE(particle_species), ALLOCATABLE, TARGET, DIMENSION(:):: species_parray

END MODULE particles

!===============================================================================
MODULE params
!===============================================================================
USE constants
INTEGER(idp) :: it,nsteps
REAL(num) :: g0,b0,dt,w0,dtcoef,tmax
REAL(num) :: theta,nlab,wlab,nc,w0_l,w0_t,lambdalab
LOGICAL :: l_coeffs_allocated= .FALSE., l_ck=.FALSE.
REAL(num), PARAMETER :: resize_factor=1.5_num
INTEGER(idp) :: topology
INTEGER(idp) :: mpicom_curr
INTEGER(isp) :: seed
INTEGER(idp) :: currdepo                            ! Current deposition method
INTEGER(idp) :: fieldgave                           ! Field gathering method
END MODULE params

!===============================================================================
MODULE mpi_type_constants !#do not parse
!===============================================================================
use mpi
use constants
INTEGER(isp)  :: mpidbl = MPI_DOUBLE_PRECISION
INTEGER(isp) :: status(MPI_STATUS_SIZE)
! Derived types (MPI exchange)
INTEGER(isp) :: derived_type_grid
INTEGER(isp) :: derived_subarray_grid
END MODULE mpi_type_constants

!===============================================================================
MODULE communications
!===============================================================================
use constants
INTEGER(isp) :: reqperjxx(4),reqperjxy(4),reqperjxz(4)
INTEGER(isp) :: reqperjyx(4),reqperjyy(4),reqperjyz(4)
INTEGER(isp) :: reqperjzx(4),reqperjzy(4),reqperjzz(4)
END MODULE communications

!===============================================================================
MODULE time_stat
!===============================================================================
use constants

INTEGER(idp) :: timestat_activated
INTEGER(idp) :: timestat_period
REAL(num), dimension(20) :: localtimes

END MODULE

!===============================================================================
MODULE output_data !#do not parse
!===============================================================================
use constants

! Simulation time statistics
REAL(num) :: startsim =0.0_num
REAL(num) :: endsim =0.0_num
REAL(num) :: startit, timeit
REAL(num) :: pushtime

! output frequency
INTEGER(idp) :: output_frequency = -1 !(Default is no output)
INTEGER(idp) :: output_step_min = 0
INTEGER(idp) :: output_step_max = 0

! output quantity flag (Default=False)
INTEGER(KIND=4) :: c_output_ex = 0
INTEGER(KIND=4) :: c_output_ey = 0
INTEGER(KIND=4) :: c_output_ez = 0
INTEGER(KIND=4) :: c_output_bx = 0
INTEGER(KIND=4) :: c_output_by = 0
INTEGER(KIND=4) :: c_output_bz = 0
INTEGER(KIND=4) :: c_output_jx = 0
INTEGER(KIND=4) :: c_output_jy = 0
INTEGER(KIND=4) :: c_output_jz = 0
INTEGER(KIND=4) :: c_output_rho = 0
INTEGER(KIND=4) :: c_output_dive = 0

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

! temporal diagnostics
INTEGER(isp), dimension(10) :: temdiag_act_list                  ! Array of activation flags
CHARACTER(len=string_length), dimension(10) :: temdiag_name_list ! Filename for the different temporal diags
INTEGER(isp), dimension(10) :: temdiag_i_list                    ! Array of index to locate the value in the big array
INTEGER(isp), dimension(10) :: temdiag_nb_values                 ! Array containing the number of values in the big array

INTEGER :: temdiag_nb
INTEGER :: temdiag_nb_part
INTEGER :: temdiag_nb_field
INTEGER :: temdiag_totvalues
INTEGER :: temdiag_frequency
INTEGER :: temdiag_format
REAL(num), dimension(:),allocatable :: temdiag_array             ! Big array containing all the temporal diag at a given iteration

END MODULE output_data

!===============================================================================
MODULE timing
!===============================================================================
use constants 
REAL(num) :: dep_curr_time=0._num

END MODULE timing


!===============================================================================
MODULE shared_data
!===============================================================================
use mpi_type_constants
USE output_data
!----------------------------------------------------------------------------
! MPI subdomain data
!----------------------------------------------------------------------------
INTEGER(isp) :: errcode, provided, comm, tag
INTEGER(idp) :: rank
INTEGER(isp) :: coordinates(3) 
INTEGER (idp) :: neighbour(-1:1, -1:1, -1:1)
INTEGER(isp) :: x_coords, proc_x_min, proc_x_max
INTEGER(isp):: y_coords, proc_y_min, proc_y_max
INTEGER(isp) :: z_coords, proc_z_min, proc_z_max
INTEGER(idp) :: nproc, nprocx, nprocy, nprocz
INTEGER(isp) :: nprocdir(3)
INTEGER(idp), POINTER, DIMENSION(:) :: nx_each_rank, ny_each_rank, nz_each_rank
LOGICAL(idp) :: x_min_boundary, x_max_boundary
LOGICAL(idp) :: y_min_boundary, y_max_boundary
LOGICAL(idp) :: z_min_boundary, z_max_boundary
! The location of the processors
INTEGER(idp), DIMENSION(:), POINTER :: cell_x_min, cell_x_max
INTEGER(idp), DIMENSION(:), POINTER :: cell_y_min, cell_y_max
INTEGER(idp), DIMENSION(:), POINTER :: cell_z_min, cell_z_max
INTEGER(idp), DIMENSION(:), POINTER :: old_x_max, old_y_max, old_z_max
INTEGER(idp) :: nx_global_grid_min, nx_global_grid_max
INTEGER(idp) :: ny_global_grid_min, ny_global_grid_max
INTEGER(idp) :: nz_global_grid_min, nz_global_grid_max
! domain and loadbalancing
LOGICAL :: allow_cpu_reduce = .FALSE.
REAL(num), DIMENSION(:), POINTER :: x_global, y_global, z_global
REAL(num), DIMENSION(:), POINTER :: xb_global, yb_global, zb_global
REAL(num), DIMENSION(:), POINTER :: xb_offset_global
REAL(num), DIMENSION(:), POINTER :: yb_offset_global
REAL(num), DIMENSION(:), POINTER :: zb_offset_global
! domain limits and size
INTEGER(idp)  :: nx, ny, nz ! local number of cells
INTEGER(idp)  :: nx_grid, ny_grid, nz_grid ! local number of grid points
INTEGER(idp)  :: nx_global, ny_global, nz_global ! global number of cells
INTEGER(idp)  :: nx_global_grid, ny_global_grid, nz_global_grid ! global number of grid points
REAL(num):: dx, xmin, xmax, length_x
REAL(num):: x_min_local, x_max_local
REAL(num):: dy, ymin, ymax,length_y
REAL(num):: y_min_local, y_max_local
REAL(num):: dz, zmin, zmax,length_z
REAL(num):: z_min_local, z_max_local

! Sorting
INTEGER(idp) :: sorting_activated                   ! Activation of soting
REAL(NUM)    :: sorting_dx, sorting_dy, sorting_dz  ! Bin space steps
REAL(NUM)    :: sorting_shiftx, sorting_shifty, sorting_shiftz

! Axis
REAL(num), POINTER, DIMENSION(:) :: x, y, z
REAL(num), DIMENSION(:), POINTER :: x_grid_mins, x_grid_maxs
REAL(num), DIMENSION(:), POINTER :: y_grid_mins, y_grid_maxs
REAL(num), DIMENSION(:), POINTER :: z_grid_mins, z_grid_maxs
REAL(num) ::  x_grid_min, x_grid_max
REAL(num) :: x_grid_min_local, x_grid_max_local
REAL(num) ::  y_grid_min, y_grid_max
REAL(num) :: y_grid_min_local, y_grid_max_local
REAL(num) :: z_grid_min, z_grid_max
REAL(num) :: z_grid_min_local, z_grid_max_local

! Total charge density
REAL(num), POINTER, DIMENSION(:,:,:) :: rho
! Electric Field divergence
REAL(num), POINTER, DIMENSION(:,:,:) :: dive

END MODULE shared_data


MODULE python_pointers
USE constants
INTEGER(idp), POINTER :: partn(:)
INTEGER(idp) :: partnmax
INTEGER(idp) :: nxtg, nytg, nztg
INTEGER(idp) :: nxgt, nygt, nzgt
INTEGER(idp) :: nxct, nyct, nzct
INTEGER(idp) :: nxmin, nxmax, nymin, nymax, &
nzmin, nzmax
! Tile position
REAL(num) :: xtmin, ytmin, ztmin
REAL(num) :: xtmax, ytmax, ztmax
REAL(num) :: xgtmin, ygtmin, zgtmin
REAL(num) :: xgtmax, ygtmax, zgtmax
REAL(num), DIMENSION(:), POINTER :: partx
REAL(num), DIMENSION(:), POINTER :: party
REAL(num), DIMENSION(:), POINTER :: partz
REAL(num), DIMENSION(:), POINTER :: partux
REAL(num), DIMENSION(:), POINTER :: partuy
REAL(num), DIMENSION(:), POINTER :: partuz
REAL(num), DIMENSION(:), POINTER :: partgaminv
REAL(num), DIMENSION(:,:), POINTER :: pid
REAL(num), DIMENSION(:), POINTER :: partex
REAL(num), DIMENSION(:), POINTER :: partey
REAL(num), DIMENSION(:), POINTER :: partez
REAL(num), DIMENSION(:), POINTER :: partbx
REAL(num), DIMENSION(:), POINTER :: partby
REAL(num), DIMENSION(:), POINTER :: partbz
END MODULE python_pointers





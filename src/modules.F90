!===============================================================================
! Contains shared data
!===============================================================================
!> Module containing the Picsar constant parameters
MODULE constants
!===============================================================================
!>
!> Float precision
INTEGER, PARAMETER :: num = 8
!> Integer 4 byte precision
INTEGER, PARAMETER :: isp = 4
!> integer double precision
INTEGER, PARAMETER :: idp = 8
INTEGER, PARAMETER :: cpx = 8
!> Electron mass
REAL(num), PARAMETER :: emass   = 9.10938291e-31_num
!> Proton mass
REAL(num), PARAMETER :: pmass   = 1.6726231000000001e-27_num
!> Electron charge
REAL(num), PARAMETER :: echarge = 1.6021764620000001e-19_num
!> Speed of light in vacuum
REAL(num), PARAMETER :: clight  = 2.99792458e8_num
!> Magnetic constant
REAL(num), PARAMETER :: mu0     = 1.2566370614359173e-06_num
!> Vacuum permeability
REAL(num), PARAMETER :: eps0    = 8.854187817620389e-12_num
REAL(num), PARAMETER :: imu0    = 795774.715459
!> The famous pi value
REAL(num), PARAMETER :: pi      = 3.14159265358979323_num
!> Dimension of the cartesian topology
INTEGER(isp), PARAMETER :: c_ndims = 3
! direction parameters
!> x direction index parameter
INTEGER, PARAMETER :: c_dir_x = 1
!> y direction index parameter
INTEGER, PARAMETER :: c_dir_y = 2
!> z direction index parameter
INTEGER, PARAMETER :: c_dir_z = 3
LOGICAL:: l_smooth_compensate
!> string length parameter for some outputs
INTEGER, PARAMETER :: string_length = 264
! Error handling
INTEGER, PARAMETER :: c_err_bad_value = 2**4
!> Vector length for some vectorized loops
INTEGER(idp), PARAMETER :: LVEC = 8
END MODULE constants

!=========================================================================================
!> Module containing useful pre-computed parameters for some subroutines
MODULE precomputed
!=========================================================================================
USE constants
!> Inverse of the space discretization:
!> \f$ 1/dx \f$
REAL(num) :: dxi
!> Inverse of the space discretization:
!> \f$ 1/dy \f$
REAL(num) :: dyi
!> Inverse of the space discretization:
!> \f$ 1/dz \f$
REAL(num) :: dzi
REAL(num) :: invvol
REAL(num) :: dts2dx,dts2dy,dts2dz
REAL(num) :: dtsdx0,dtsdy0,dtsdz0
REAL(num) :: dxs2,dys2,dzs2
REAL(num) :: clightsq
END MODULE precomputed

!=========================================================================================
!> Module containing parameters and data structures for the fields
MODULE fields
!=========================================================================================
USE constants
LOGICAL:: l_lower_order_in_v, l_nodalgrid, l4symtry
INTEGER(idp):: nxs=0, nys=0, nzs=0
INTEGER(idp):: norderx, nordery, norderz
INTEGER(idp):: nxguards,nyguards, nzguards, nox, noy, noz, npass(3)
INTEGER(idp):: nxjguards,nyjguards, nzjguards    ! Guard cells for current arrays
REAL(num):: alpha(3)
REAL(num), POINTER, DIMENSION(:,:,:) :: ex,ey,ez,bx,by,bz,jx,jy,jz
! Fonberg coefficients
REAL(num), POINTER, DIMENSION(:) :: xcoeffs, ycoeffs, zcoeffs
END MODULE fields

! ________________________________________________________________________________________
!> Module containing the field tile data structure
MODULE grid_tilemodule !#do not parse
!
!
! ________________________________________________________________________________________
USE constants
TYPE grid_tile
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile, eytile, eztile, &
 								bxtile, bytile, bztile,  jxtile, jytile, jztile,rhotile
#if defined __INTEL_COMPILER
    !dir$ attributes align:64 :: extile
    !dir$ attributes align:64 :: eytile
    !dir$ attributes align:64 :: eztile
    !dir$ attributes align:64 :: bxtile
    !dir$ attributes align:64 :: bytile
    !dir$ attributes align:64 :: bztile
    !dir$ attributes align:64 :: jxtile
    !dir$ attributes align:64 :: jytile
    !dir$ attributes align:64 :: jztile
    !dir$ attributes align:64 :: rhotile
#endif
		! FastMEM attributes to manage where the data are allocated
    !DIR ATTRIBUTES FASTMEM  :: extile
    !DIR ATTRIBUTES FASTMEM  :: eytile
    !DIR ATTRIBUTES FASTMEM  :: eztile
    !DIR ATTRIBUTES FASTMEM  :: bxtile
    !DIR ATTRIBUTES FASTMEM  :: bytile
    !DIR ATTRIBUTES FASTMEM  :: bztile
    !DIR ATTRIBUTES FASTMEM  :: jxtile
    !DIR ATTRIBUTES FASTMEM  :: jytile
    !DIR ATTRIBUTES FASTMEM  :: jztile
    !DIR ATTRIBUTES FASTMEM  :: rhotile

END TYPE
TYPE(grid_tile), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: aofgrid_tiles
END MODULE grid_tilemodule

!===============================================================================
!> Module containing the Fortran object descriptor representing a particle tile
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
    !DIR ATTRIBUTES FASTMEM  :: part_x
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_y
    !DIR ATTRIBUTES FASTMEM  :: part_y
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_z
    !DIR ATTRIBUTES FASTMEM  :: part_z
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ux
    !DIR ATTRIBUTES FASTMEM  :: part_ux
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uy
    !DIR ATTRIBUTES FASTMEM  :: part_uy
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uz
    !DIR ATTRIBUTES FASTMEM  :: part_uz
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_gaminv
    !DIR ATTRIBUTES FASTMEM  :: part_gaminv
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ex
    !DIR ATTRIBUTES FASTMEM  :: part_ex
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ey
    !DIR ATTRIBUTES FASTMEM  :: part_ey
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ez
    !DIR ATTRIBUTES FASTMEM  :: part_ez
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_bx
    !DIR ATTRIBUTES FASTMEM  :: part_bx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_by
    !DIR ATTRIBUTES FASTMEM  :: part_by
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_bz
    !DIR ATTRIBUTES FASTMEM  :: part_bz
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
    !DIR ATTRIBUTES FASTMEM  :: pid
#if defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_x
    !dir$ attributes align:64 :: part_y
    !dir$ attributes align:64 :: part_z
    !dir$ attributes align:64 :: part_ux
    !dir$ attributes align:64 :: part_uy
    !dir$ attributes align:64 :: part_uz
    !dir$ attributes align:64 :: part_gaminv
    !dir$ attributes align:64 :: part_ex
    !dir$ attributes align:64 :: part_ey
    !dir$ attributes align:64 :: part_ez
    !dir$ attributes align:64 :: part_bx
    !dir$ attributes align:64 :: part_by
    !dir$ attributes align:64 :: part_bz
    !dir$ attributes align:64 :: pid
#endif
END TYPE
END MODULE particle_tilemodule

!===============================================================================
!> Module containing the Fortran object descriptor representing a particle species
MODULE particle_speciesmodule !#do not parse
!===============================================================================
use particle_tilemodule
use constants
!> Fortran object representing a particle species
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
    INTEGER(idp)   :: sorting_start     ! Sorting start iteration
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
!> Module for the tile parameters
MODULE tile_params
!===============================================================================
! # of particle tiles in each dimension
USE constants
!> Number of tile in the x direction
INTEGER(idp) :: ntilex
!> Number of tile in the y direction
INTEGER(idp) :: ntiley
!> Number of tile in the z direction
INTEGER(idp) :: ntilez
END MODULE tile_params


!===============================================================================
!> Module containing useful properties for the particles
MODULE particle_properties
!===============================================================================
	USE constants
	INTEGER(idp), PARAMETER :: npid=1
	INTEGER(idp), PARAMETER :: wpid=1
	LOGICAL :: l_initongrid = .FALSE.
	LOGICAL :: l_particles_weight = .FALSE.
  !> Particle pusher type (0: Boris, 1: Vay, Default: 0)
  INTEGER(idp) :: particle_pusher = 0
	!> Particle initial distribution
	INTEGER(idp) :: pdistr
	!> Number of species
	INTEGER(idp) :: nspecies
	!> total number of particles (all species, all subdomains -> useful for stat)
	INTEGER(idp) :: ntot
	!> Max number of particle species
	INTEGER(idp) :: nspecies_max=6
	REAL(num) :: fdxrand=0.0_num,fdzrand=0.0_num,vthx=0.0_num,vthy=0.0_num,vthz=0.0_num
	LOGICAL :: l_species_allocated=.FALSE.,  l_pdumps_allocated=.FALSE.
END MODULE particle_properties



!===============================================================================
!> Module containing the array of species
MODULE particles !#do not parse
!===============================================================================
USE constants
USE tile_params
USE particle_tilemodule
USE particle_speciesmodule
USE particle_properties
USE grid_tilemodule

	!> Array of  particle species objects
	TYPE(particle_species), ALLOCATABLE, TARGET, DIMENSION(:) :: species_parray

END MODULE particles

!=========================================================================================
!> Module containing useful configuration and simulation parameters
MODULE params
!=========================================================================================
USE constants

	!> iteration number
	INTEGER(idp)         :: it
	!> Total number of steps
	INTEGER(idp)         :: nsteps
	REAL(num)            :: g0,b0,dt,w0,dtcoef,tmax
	REAL(num)            :: theta,nlab,wlab,nc,w0_l,w0_t,lambdalab
	LOGICAL              :: l_coeffs_allocated= .FALSE., l_ck=.FALSE.
	REAL(num), PARAMETER :: resize_factor=2._num
	REAL(num), PARAMETER :: downsize_factor=0.5_num
	REAL(num), PARAMETER :: downsize_threshold=0.4_num
	INTEGER(idp) :: topology
	INTEGER(idp) :: mpicom_curr
	INTEGER(isp) :: seed
	!> Current deposition method
	INTEGER(idp) :: currdepo
	!> Charge deposition method
	INTEGER(idp) :: rhodepo
	!> Field gathering method
	INTEGER(idp) :: fieldgathe
	!> Type of comm routine to use for particles
	INTEGER(idp) :: partcom
	!> Field gathering + part. pusher seperated flag
	INTEGER(idp) :: fg_p_pp_separated
	!> Vector size for the current deposition
	INTEGER(idp) :: LVEC_curr_depo
	!> Vector size for the charge deposition
	INTEGER(idp) :: LVEC_charge_depo
	!> Vector size for the field gathering
	INTEGER(idp) :: LVEC_fieldgathe
	!> MPI buffer size
	INTEGER(isp) :: mpi_buf_size

END MODULE params

!===============================================================================
!> Module containing MPI parameters
MODULE mpi_type_constants !#do not parse
!===============================================================================
use mpi
use constants
INTEGER(isp)  :: mpidbl = MPI_DOUBLE_PRECISION
INTEGER(isp) :: status(MPI_STATUS_SIZE)
! Derived types (MPI exchange)
INTEGER(isp) :: derived_type_grid
INTEGER(isp) :: derived_subarray_grid
INTEGER(isp), DIMENSION(100) :: mpi_dtypes
LOGICAL(isp), DIMENSION(100) :: is_dtype_init = .TRUE.
END MODULE mpi_type_constants

!===============================================================================
!> Module for the communications
MODULE communications  !#do not parse
!===============================================================================
use constants
INTEGER(isp) :: reqperjxx(4),reqperjxy(4),reqperjxz(4)
INTEGER(isp) :: reqperjyx(4),reqperjyy(4),reqperjyz(4)
INTEGER(isp) :: reqperjzx(4),reqperjzy(4),reqperjzz(4)

  TYPE part_com_buffer
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_x
      !dir$ attributes align:64 :: part_x
      !DIR ATTRIBUTES FASTMEM  :: part_x
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_y
      !dir$ attributes align:64 :: part_y
      !DIR ATTRIBUTES FASTMEM  :: part_y
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_z
      !dir$ attributes align:64 :: part_z
      !DIR ATTRIBUTES FASTMEM  :: part_z
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ux
      !dir$ attributes align:64 :: part_ux
      !DIR ATTRIBUTES FASTMEM  :: part_ux
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uy
      !dir$ attributes align:64 :: part_uy
      !DIR ATTRIBUTES FASTMEM  :: part_uy
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uz
      !dir$ attributes align:64 :: part_uz
      !DIR ATTRIBUTES FASTMEM  :: part_uz
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_gaminv
      !dir$ attributes align:64 :: part_gaminv
      !DIR ATTRIBUTES FASTMEM  :: part_gaminv
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
      !dir$ attributes align:64 :: pid
      !DIR ATTRIBUTES FASTMEM :: pid
      INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: boundid
      INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: bin_npart
      INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: bin_pos
  END TYPE

  TYPE mpi_buffer
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_x
      !dir$ attributes align:64 :: part_x
      !DIR ATTRIBUTES FASTMEM  :: part_x
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_y
      !dir$ attributes align:64 :: part_y
      !DIR ATTRIBUTES FASTMEM  :: part_y
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_z
      !dir$ attributes align:64 :: part_z
      !DIR ATTRIBUTES FASTMEM  :: part_z
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_ux
      !dir$ attributes align:64 :: part_ux
      !DIR ATTRIBUTES FASTMEM  :: part_ux
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_uy
      !dir$ attributes align:64 :: part_uy
      !DIR ATTRIBUTES FASTMEM  :: part_uy
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_uz
      !dir$ attributes align:64 :: part_uz
      !DIR ATTRIBUTES FASTMEM  :: part_uz
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_gaminv
      !dir$ attributes align:64 :: part_gaminv
      !DIR ATTRIBUTES FASTMEM  :: part_gaminv
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
      !dir$ attributes align:64 :: pid
      !DIR ATTRIBUTES FASTMEM  :: pid
      INTEGER(idp), dimension(27) :: npart
  END TYPE


END MODULE communications

!===============================================================================
!> Module for the time statistics
MODULE time_stat !#do not parse
!===============================================================================
use constants

	!> Activation of the outputs
	INTEGER(idp)                           :: timestat_activated
	!> Period for the outputs
	INTEGER(idp)                           :: timestat_period
	!> First iteration for the time statistics
	INTEGER(idp)                           :: timestat_itstart
	!> ! Flag to activate the time statistics per iteration
	INTEGER(idp)                           :: timestat_perit
	!> MPI local times for the initialization
	REAL(num), dimension(5)                :: init_localtimes

	!> MPI local times for the main loop
	REAL(num), dimension(20)               :: localtimes
	!> Buffer for the output
	REAL(num), DIMENSION(:,:), POINTER     :: buffer_timestat
	INTEGER(idp)                           :: itimestat
	!> Number of entries in the buffer
	INTEGER(idp)                           :: nbuffertimestat

END MODULE

!===============================================================================
!> Module for the outputs
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
	!> Array of activation flags
	INTEGER(isp), dimension(15) :: temdiag_act_list
	!> Filename for the different temporal diags
	CHARACTER(len=string_length), dimension(10) :: temdiag_name_list
	!> Array of index to locate the value in the big array
	INTEGER(isp), dimension(15) :: temdiag_i_list
	!> Array containing the number of values in the big array
	INTEGER(isp), dimension(15) :: temdiag_nb_values

  !> Number of temoral diags
	INTEGER(idp) :: temdiag_nb
	!> Number of particle temporal diagnostics
	INTEGER(idp) :: temdiag_nb_part
	INTEGER(idp) :: temdiag_nb_field
	INTEGER(idp) :: temdiag_totvalues
	INTEGER(idp) :: temdiag_frequency
	INTEGER(idp) :: temdiag_format
	!> Big array containing all the temporal diag at a given iteration
	REAL(num), dimension(:),allocatable :: temdiag_array

	! Computation flags
	LOGICAL      :: divE_computed

	INTEGER(idp) :: npdumps
	TYPE particle_dump
		INTEGER(idp) :: ispecies
		INTEGER(idp) :: diag_period
		REAL(num) :: dump_x_min, dump_x_max
		REAL(num) :: dump_y_min, dump_y_max
		REAL(num) :: dump_z_min, dump_z_max
		REAL(num) :: dump_ux_min, dump_ux_max
		REAL(num) :: dump_uy_min, dump_uy_max
		REAL(num) :: dump_uz_min, dump_uz_max
	END TYPE particle_dump

  !> Object for the particle dumping
	TYPE(particle_dump), ALLOCATABLE, TARGET, DIMENSION(:) :: particle_dumps
END MODULE output_data

!===============================================================================
!> Module for the timing. This module should be merged with the time statistics
MODULE timing
!===============================================================================
	use constants
	REAL(num) :: dep_curr_time=0._num

END MODULE timing


!===============================================================================
!> Module for the data shared with Python.
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
	INTEGER(idp) :: x_coords, proc_x_min, proc_x_max
	INTEGER(idp):: y_coords, proc_y_min, proc_y_max
	INTEGER(idp)                        :: z_coords, proc_z_min, proc_z_max
	INTEGER(idp)                        :: nproc, nprocx, nprocy, nprocz
	INTEGER(isp)                        :: nprocdir(3)
	INTEGER(idp), POINTER, DIMENSION(:) :: nx_each_rank, ny_each_rank, nz_each_rank
	! Boundary data
	LOGICAL(idp)                        :: x_min_boundary, x_max_boundary
	LOGICAL(idp)                        :: y_min_boundary, y_max_boundary
	LOGICAL(idp)                        :: z_min_boundary, z_max_boundary
	INTEGER(idp)                        :: pbound_x_min, pbound_x_max
	INTEGER(idp)                        :: pbound_y_min, pbound_y_max
	INTEGER(idp)                        :: pbound_z_min, pbound_z_max

	! The location of the processors
	INTEGER(idp), DIMENSION(:), POINTER :: cell_x_min, cell_x_max
	INTEGER(idp), DIMENSION(:), POINTER :: cell_y_min, cell_y_max
	INTEGER(idp), DIMENSION(:), POINTER :: cell_z_min, cell_z_max
	INTEGER(idp), DIMENSION(:), POINTER :: new_cell_x_min, new_cell_x_max
	INTEGER(idp), DIMENSION(:), POINTER :: new_cell_y_min, new_cell_y_max
	INTEGER(idp), DIMENSION(:), POINTER :: new_cell_z_min, new_cell_z_max
	INTEGER(idp)                        :: nx_global_grid_min, nx_global_grid_max
	INTEGER(idp)                        :: ny_global_grid_min, ny_global_grid_max
	INTEGER(idp)                        :: nz_global_grid_min, nz_global_grid_max
	! Domain axis
	LOGICAL(idp)                        :: l_axis_allocated=.FALSE.
	REAL(num), DIMENSION(:), POINTER    :: x_global, y_global, z_global
	REAL(num), DIMENSION(:), POINTER    :: xb_global, yb_global, zb_global
	REAL(num), DIMENSION(:), POINTER    :: xb_offset_global
	REAL(num), DIMENSION(:), POINTER    :: yb_offset_global
	REAL(num), DIMENSION(:), POINTER    :: zb_offset_global
	! domain limits and size
	!> local number of cells
	INTEGER(idp)                        :: nx, ny, nz
	!> local number of grid points
	INTEGER(idp)                        :: nx_grid, ny_grid, nz_grid
	!> global number of cells
	INTEGER(idp)                        :: nx_global, ny_global, nz_global
	!> global number of grid points
	INTEGER(idp)                        :: nx_global_grid, ny_global_grid, nz_global_grid
	REAL(num)                           :: dx, xmin, xmax, length_x
	REAL(num)                           :: x_min_local, x_max_local
	REAL(num)                           :: dy, ymin, ymax,length_y
	REAL(num)                           :: y_min_local, y_max_local
	REAL(num)                           :: dz, zmin, zmax,length_z
	REAL(num)                           :: z_min_local, z_max_local

	! Sorting
	!> Activation of the sorting
	INTEGER(idp) :: sorting_activated
	!> Bin space steps
	REAL(NUM)    :: sorting_dx, sorting_dy, sorting_dz
	!> Shift of the sorting grid in respect of the origin
	REAL(NUM)    :: sorting_shiftx, sorting_shifty, sorting_shiftz
	!> verbose for the sorting (depreciated)
	LOGICAL      :: sorting_verbose

	! Axis
	!> Space dimension
	INTEGER(idp) :: c_dim = 3
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

  !> Moving window offset z-position
  REAL(num) :: zgrid =0.

	!> Total charge density
	REAL(num), POINTER, DIMENSION(:,:,:) :: rho
	!> Electric Field divergence
	REAL(num), POINTER, DIMENSION(:,:,:) :: dive

	! Values used for load balancing
	REAL(num) :: mpitime_per_it, max_time_per_it, min_time_per_it
	REAL(num) :: global_time_per_cell, global_time_per_part
	REAL(num) :: local_time_cell, local_time_part
	INTEGER(idp) :: npart_local, npart_global
END MODULE shared_data

!===============================================================================
!> Module for the Maxwell Solver coefficients
MODULE kyee_em3d
!===============================================================================
	USE constants
	REAL(num) :: alphax = 0.58333333333333337  ! 7./12.
	REAL(num) :: betaxy = 0.083333333333333329 ! 1./12.
	REAL(num) :: betaxz = 0.083333333333333329 ! 1./12.
	REAL(num) :: gammax = 0.020833333333333332 ! 1./48.
	REAL(num) :: alphay = 0.58333333333333337  ! 7./12.
	REAL(num) :: betayx = 0.083333333333333329 ! 1./12.
	REAL(num) :: betayz = 0.083333333333333329 ! 1./12.
	REAL(num) :: gammay = 0.020833333333333332 ! 1./48.
	REAL(num) :: alphaz = 0.58333333333333337  ! 7./12.
	REAL(num) :: betazx = 0.083333333333333329 ! 1./12.
	REAL(num) :: betazy = 0.083333333333333329 ! 1./12.
	REAL(num) :: gammaz = 0.020833333333333332 ! 1./48.
	REAL(num) :: deltaz = 0.000000000000000000 ! for the lehe solver
END MODULE kyee_em3d

!=========================================================================================
!> Module containing pointer to the python arrays
MODULE python_pointers
!=========================================================================================
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

  !> array for particle x position
	REAL(num), DIMENSION(:), POINTER :: partx
	!dir$ attributes align:64 :: partx
	!DIR ATTRIBUTES FASTMEM  :: partx
  !> array for particle y position
	REAL(num), DIMENSION(:), POINTER :: party
	!dir$ attributes align:64 :: party
	!DIR ATTRIBUTES FASTMEM  :: party
  !> array for particle z position
	REAL(num), DIMENSION(:), POINTER :: partz
	!dir$ attributes align:64 :: partz
	!DIR ATTRIBUTES FASTMEM  :: partz
  !> array for particle x momentum
	REAL(num), DIMENSION(:), POINTER :: partux
	!dir$ attributes align:64 :: partux
	!DIR ATTRIBUTES FASTMEM  :: partux
  !> array for particle y momentum
	REAL(num), DIMENSION(:), POINTER :: partuy
	!dir$ attributes align:64 :: partuy
	!DIR ATTRIBUTES FASTMEM  :: partuy
  !> array for particle z momentum
	REAL(num), DIMENSION(:), POINTER :: partuz
	!dir$ attributes align:64 :: partuz
	!DIR ATTRIBUTES FASTMEM  :: partuz
  !> array for the inverse of the particle gamma factor
	REAL(num), DIMENSION(:), POINTER :: partgaminv
	!dir$ attributes align:64 :: partgaminv
	!DIR ATTRIBUTES FASTMEM  :: partgaminv
	REAL(num), DIMENSION(:,:), POINTER :: pid
	!dir$ attributes align:64 :: pid
	!DIR ATTRIBUTES FASTMEM  :: pid
	REAL(num), DIMENSION(:), POINTER :: partex
	!dir$ attributes align:64 :: partex
	!DIR ATTRIBUTES FASTMEM  :: partex
	REAL(num), DIMENSION(:), POINTER :: partey
	!dir$ attributes align:64 :: partey
	!DIR ATTRIBUTES FASTMEM  :: partey
	REAL(num), DIMENSION(:), POINTER :: partez
	!dir$ attributes align:64 :: partez
	!DIR ATTRIBUTES FASTMEM  :: partez
	REAL(num), DIMENSION(:), POINTER :: partbx
	!dir$ attributes align:64 :: partbx
	!DIR ATTRIBUTES FASTMEM  :: partbx
	REAL(num), DIMENSION(:), POINTER :: partby
	!dir$ attributes align:64 :: partby
	!DIR ATTRIBUTES FASTMEM  :: partby
	REAL(num), DIMENSION(:), POINTER :: partbz
	!dir$ attributes align:64 :: partbz
	!DIR ATTRIBUTES FASTMEM  :: partbz
END MODULE python_pointers

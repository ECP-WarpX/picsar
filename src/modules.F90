! ==============================================================================
!
! MODULES.F90
!
! This file contains several modules with parameters for different purposes.
!
! ==============================================================================


!===============================================================================
! Contains shared data
!===============================================================================
!> @brief
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
!> logical precision
INTEGER, PARAMETER :: lp = 8

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
LOGICAL(lp):: l_smooth_compensate
!> string length parameter for some outputs
INTEGER, PARAMETER :: string_length = 264
!> Error handling
INTEGER, PARAMETER :: c_err_bad_value = 2**4
!> Vector length for some vectorized loops
INTEGER(idp), PARAMETER :: LVEC = 8
END MODULE constants

!=========================================================================================
!> @brief
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
!> @brief
!> Module containing parameters and data structures for the fields
MODULE fields
!=========================================================================================
  USE constants
  !> Flag: interpolation at a lower order for the field gathering
  LOGICAL(lp) :: l_lower_order_in_v
  !> Flag: use of nodal grids
  LOGICAL(lp) :: l_nodalgrid
  !> Flag: this flag needs a description, used in field gathering routines
  LOGICAL(lp) :: l4symtry
  INTEGER(idp):: nxs=0, nys=0, nzs=0
  !> order in x of the FDTD Maxwell solver
  INTEGER(idp):: norderx
  !> order in y of the FDTD Maxwell solver
  INTEGER(idp):: nordery
  !> order in z of the FDTD Maxwell solver
  INTEGER(idp):: norderz
  !> Number of guard cells in x
  INTEGER(idp):: nxguards
  !> Number of guard cells in y
  INTEGER(idp):: nyguards
  !> Number of guard cells in z
  INTEGER(idp):: nzguards
  !> interpolation order in x for the field gathering
  INTEGER(idp):: nox
  !> interpolation order in y for the field gathering
  INTEGER(idp):: noy
  !> interpolation order in z for the field gathering
  INTEGER(idp):: noz
  !> this parameter needs a description
  INTEGER(idp):: npass(3)
  !> Number of guard cells in x for the current deposition
  INTEGER(idp):: nxjguards
  !> Number of guard cells in y for the current deposition
  INTEGER(idp):: nyjguards
  !> Number of guard cells in z for the current deposition
  INTEGER(idp):: nzjguards
  !> Coefficient for the Maxwell solver, exact purpose needs a better explanation
  REAL(num):: alpha(3)
  !> MPI-domain electric field grid in x
  REAL(num), POINTER, DIMENSION(:,:,:) :: ex
  !> MPI-domain electric field grid in y
  REAL(num), POINTER, DIMENSION(:,:,:) :: ey
  !> MPI-domain electric field grid in z
  REAL(num), POINTER, DIMENSION(:,:,:) :: ez
  !> MPI-domain magnetic field grid in x
  REAL(num), POINTER, DIMENSION(:,:,:) :: bx
  !> MPI-domain magnetic field grid in y
  REAL(num), POINTER, DIMENSION(:,:,:) :: by
  !> MPI-domain magnetic field grid in z
  REAL(num), POINTER, DIMENSION(:,:,:) :: bz
  !> MPI-domain current grid in x
  REAL(num), POINTER, DIMENSION(:,:,:) :: jx
  !> MPI-domain current grid in y
  REAL(num), POINTER, DIMENSION(:,:,:) :: jy
  !> MPI-domain current grid in z
  REAL(num), POINTER, DIMENSION(:,:,:) :: jz
  !> Fonberg coefficients in x
  REAL(num), POINTER, DIMENSION(:) :: xcoeffs
  !> Fonberg coefficients in y
  REAL(num), POINTER, DIMENSION(:) :: ycoeffs
  !> Fonberg coefficients in z
  REAL(num), POINTER, DIMENSION(:) :: zcoeffs
  
END MODULE fields

! ________________________________________________________________________________________
!> @brief
!> Module containing the field tile data structure.
MODULE grid_tilemodule !#do not parse
! ________________________________________________________________________________________

  USE constants
  
  !> This object contains 3D field grids for one tile
  TYPE grid_tile
      !> Tile Electric field grid in x
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: extile
      !> Tile Electric field grid in y
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: eytile
      !> Tile Electric field grid in z
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: eztile
      !> Tile Magnetic field grid in x
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bxtile
      !> Tile Magnetic field grid in y
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bytile
      !> Tile Magnetic field grid in z
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bztile
      !> Tile Current grid in x
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jxtile
      !> Tile Current grid in y
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jytile
      !> Tile Current grid in z
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jztile
      !> Tile Charge grid
      REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: rhotile
      
! We declare arrays aligned for vectorization efficiency.
! These directives are only understood by the Intel compiler.
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
    ! FastMEM attributes to manage where the data are allocated (DDR/MCDRAM)
    ! This is not activated.
    !!DIR ATTRIBUTES FASTMEM  :: extile
    !!DIR ATTRIBUTES FASTMEM  :: eytile
    !!DIR ATTRIBUTES FASTMEM  :: eztile
    !!DIR ATTRIBUTES FASTMEM  :: bxtile
    !!DIR ATTRIBUTES FASTMEM  :: bytile
    !!DIR ATTRIBUTES FASTMEM  :: bztile
    !!DIR ATTRIBUTES FASTMEM  :: jxtile
    !!DIR ATTRIBUTES FASTMEM  :: jytile
    !!DIR ATTRIBUTES FASTMEM  :: jztile
    !!DIR ATTRIBUTES FASTMEM  :: rhotile

  END TYPE
  
  !> This array contains a list of tile field grids.
  !> This array is local to each MPI domain.
  !> Tile grids are contained in the object grid_tile (see extile,eytile...).
  TYPE(grid_tile), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: aofgrid_tiles
  
END MODULE grid_tilemodule

!===============================================================================
!> @brief
!> Module containing the Fortran object descriptor representing a particle tile.
!> Also see tiling.F90 for the definition of the tile properties.
MODULE particle_tilemodule !#do not parse
!===============================================================================
USE constants

!> Object that contains tile particle arrays and particle tile properties.
TYPE particle_tile
    !> Flag: tile arrays are allocated
    LOGICAL(lp) :: l_arrays_allocated= .FALSE.
    !> Current number of particles in tile
    INTEGER(idp), DIMENSION(1) :: np_tile
    !> Max number of particles per tile: size of the arrays
    INTEGER(idp) :: npmax_tile
    !> Number of guard cells in x
    INTEGER(idp) :: nxg_tile
    !> Number of guard cells in y
    INTEGER(idp) :: nyg_tile
    !> Number of guard cells in z
    INTEGER(idp) :: nzg_tile
    !> Number of nodes in x
    INTEGER(idp) :: nx_grid_tile
    !> Number of nodes in y
    INTEGER(idp) :: ny_grid_tile
    !> Number of nodes in z
    INTEGER(idp) :: nz_grid_tile
    !> Number of cells in x
    INTEGER(idp) :: nx_cells_tile
    !> Number of cells in y
    INTEGER(idp) :: ny_cells_tile
    !> Number of cells in z
    INTEGER(idp) :: nz_cells_tile
    !> Minimal cell index in x
    INTEGER(idp) :: nx_tile_min
    !> Maximal cell index in x
    INTEGER(idp) :: nx_tile_max
    !> Minimal cell index in y
    INTEGER(idp) :: ny_tile_min
    !> Maximal cell index in y
    INTEGER(idp) :: ny_tile_max
    !> Minimal cell index in z
    INTEGER(idp) :: nz_tile_min
    !> Maximal cell index in z
    INTEGER(idp) :: nz_tile_max
    ! Tile position
    !> Minimal tile boundary in x
    REAL(num) :: x_tile_min
    !> Minimal tile boundary in y
    REAL(num) :: y_tile_min
    !> Minimal tile boundary in z
    REAL(num) :: z_tile_min
    !> Maximal tile boundary in x
    REAL(num) :: x_tile_max
    !> Maximal tile boundary in y
    REAL(num) :: y_tile_max
    !> Maximal tile boundary in z
    REAL(num) :: z_tile_max
    
    !> Minimal grid tile boundary in x
    REAL(num) :: x_grid_tile_min
    !> Minimal grid tile boundary in y
    REAL(num) :: y_grid_tile_min
    !> Minimal grid tile boundary in z
    REAL(num) :: z_grid_tile_min
    
    !> Maximal grid tile boundary in x
    REAL(num) :: x_grid_tile_max
    !> Maximal grid tile boundary in y
    REAL(num) :: y_grid_tile_max
    !> Maximal grid tile boundary in z
    REAL(num) :: z_grid_tile_max
    !> Subdomain border flags
    LOGICAL(lp) :: subdomain_bound = .FALSE.
    ! Particle arrays
    !> Particle x position array in the tile
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_x
    !> Particle y position array in the tile
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_y
    !> Particle z position array in the tile
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_z
    !> Particle normalized momentum array in x (px/mc)
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ux
    !> Particle normalized momentum array in y (py/mc)
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uy
    !> Particle normalized momentum array in z (pz/mc)
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uz
    !> Particle Lorentz factor array
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_gaminv
    !> Particle electric field array in x
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ex
    !> Particle electric field array in y
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ey
    !> Particle electric field array in z
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ez
    !> Particle magnetic field array in x
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_bx
    !> Particle magnetic field array in y
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_by
    !> Particle magnetic field array in z
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_bz
    !> Particle weight array
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
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
    ! FastMEM attributes to manage where the data are allocated (DDR/MCDRAM)
    ! This is not activated.
    !DIR ATTRIBUTES FASTMEM  :: part_x
    !DIR ATTRIBUTES FASTMEM  :: part_y
    !DIR ATTRIBUTES FASTMEM  :: part_z
    !DIR ATTRIBUTES FASTMEM  :: part_ux
    !DIR ATTRIBUTES FASTMEM  :: part_uy
    !DIR ATTRIBUTES FASTMEM  :: part_uz
    !DIR ATTRIBUTES FASTMEM  :: part_gaminv
    !DIR ATTRIBUTES FASTMEM  :: part_ex
    !DIR ATTRIBUTES FASTMEM  :: part_ey
    !DIR ATTRIBUTES FASTMEM  :: part_ez
    !DIR ATTRIBUTES FASTMEM  :: part_bx
    !DIR ATTRIBUTES FASTMEM  :: part_by
    !DIR ATTRIBUTES FASTMEM  :: part_bz
    !DIR ATTRIBUTES FASTMEM  :: pid
END TYPE
END MODULE particle_tilemodule

!===============================================================================
!> @brief
!> Module containing the Fortran object descriptor representing a particle species
MODULE particle_speciesmodule !#do not parse
!===============================================================================
  use particle_tilemodule
  use constants
  !> Fortran object representing a particle species
  TYPE particle_species
      ! Attributes of particle species object
      !> Particle species name
      CHARACTER(LEN=string_length) :: name
      !> Particle species charge
      REAL(num) :: charge
      !> Particle species mass
      REAL(num) :: mass
      !> Particle minimal x position for initialization
      REAL(num) :: x_min
      !> Particle maximal x position for initialization
      REAL(num) :: x_max
      !> Particle minimal y position for initialization
      REAL(num) :: y_min
      !> Particle maximal y position for initialization
      REAL(num) :: y_max
      !> Particle minimal z position for initialization
      REAL(num) :: z_min
      !> Particle maximal z position for initialization
      REAL(num) :: z_max
      !> Particle drift velocity in x
      REAL(num) :: vdrift_x
      !> Particle drift velocity in y
      REAL(num) :: vdrift_y
      !> Particle drift velocity in z
      REAL(num) :: vdrift_z
      !> Particle thermic velocity in x (for distribution drawing)
      REAL(num) :: vth_x
      !> Particle thermic velocity in y (for distribution drawing)
      REAL(num) :: vth_y
      !> Particle thermic velocity in z (for distribution drawing)
      REAL(num) :: vth_z
      !> Number of particles for this species
      INTEGER(idp)   :: species_npart
      !> Maximal number of particles that arrays can contain
      INTEGER(idp)   :: nppspecies_max
      !> Number of particles per cell
      INTEGER(idp)   :: nppcell
      !> Sorting period
      INTEGER(idp)   :: sorting_period
      !> Sorting start iteration
      INTEGER(idp)   :: sorting_start 
      !> Flag indicating of the array array_of_tile has been allocated
      LOGICAL(lp)   :: l_arrayoftiles_allocated =.FALSE.
      ! For some stupid reason, cannot use ALLOCATABLE in derived types
      ! in Fortran 90 - Need to use POINTER instead
      !> List of tiles (of objects particle_tile) in the MPI domain for the 
      !> particles of this species.
      TYPE(particle_tile), DIMENSION(:,:,:), ALLOCATABLE :: array_of_tiles
      !> Array indicating if a tile has been reallocated.
      !> Used for interfacing WARP and PXR
      INTEGER(idp), DIMENSION(:,:,:), ALLOCATABLE :: are_tiles_reallocated
  END TYPE
END MODULE particle_speciesmodule

!===============================================================================
!> @brief
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
  !> Number of elements per particle in the pid particle array
  INTEGER(idp), PARAMETER :: npid=1
  !> Dimension corresponding to the weight in pid
  INTEGER(idp), PARAMETER :: wpid=1
  !> This flag seems to be unused
  LOGICAL(lp) :: l_initongrid = .FALSE.
  !> Flag to activate the use of weight for the particles
  LOGICAL(lp) :: l_particles_weight = .FALSE.
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
  !> this parameter it not used
  REAL(num) :: fdxrand=0.0_num
  !> this parameter it not used
  REAL(num) :: fdzrand=0.0_num
  REAL(num) :: vthx=0.0_num
  REAL(num) :: vthy=0.0_num
  REAL(num) :: vthz=0.0_num
  !> Flag for the allocation of the species array
  LOGICAL(lp) :: l_species_allocated=.FALSE.
  !> Flag for the allocation of the particle dump array
  LOGICAL(lp) :: l_pdumps_allocated=.FALSE.
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

  !> Array of particle species objects
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
  !> Initial frame Gamma factor (in the case of a moving frame)
  REAL(num)            :: g0
  !> Initial normalized frame velocity (sqrt(1.0_num-1.0_num/g0**2))
  REAL(num)            :: b0
  !> Time step
  REAL(num)            :: dt
  !> "longitudinal" plasma frequency in the frame of reference
  REAL(num)            :: w0
  !> Factor on the time step (> 0 and <= 1 to respect the CFL condition)
  REAL(num)            :: dtcoef
  !> Final simulation time
  REAL(num)            :: tmax
  !> The purpose of this variable us unknown and this variable is not used
  REAL(num)            :: theta
  !> Plasma density in lab frame
  REAL(num)            :: nlab
  !> Plasma frequency in the lab frame (echarge*sqrt(nlab/(emass*eps0)))
  REAL(num)            :: wlab
  !> Density in the simulation frame (nlab*g0)
  REAL(num)            :: nc
  !> "longitudinal" plasma frequency in the lab frame
  REAL(num)            :: w0_l
  !> transverse plasma frequency
  REAL(num)            :: w0_t
  !> Cold plasma wavelength (2*pi*clight/wlab)
  REAL(num)            :: lambdalab
  !> Flag True when stencil coefficients for Maxwell have been allocated
  LOGICAL(lp)          :: l_coeffs_allocated= .FALSE.
  !> This purpose of this parameter is unknown and it is not used
  LOGICAL(lp)          :: l_ck=.FALSE.
  !> This factor is used to resize particle arrays when full
  REAL(num), PARAMETER :: resize_factor=2._num
  !> This factor is used to diminish the size of particle arrays
  REAL(num), PARAMETER :: downsize_factor=0.5_num
  !> This factor is used to downsize the maximal number of particles
  REAL(num), PARAMETER :: downsize_threshold=0.4_num
  !> Type of MPI topology
  INTEGER(idp) :: topology
  !> Type if current MPI communication
  INTEGER(idp) :: mpicom_curr
  !> Seed for random drawings 
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
  INTEGER(idp) :: mpi_buf_size

END MODULE params

!===============================================================================
!> Module containing MPI parameters
MODULE mpi_type_constants !#do not parse
!===============================================================================
  use mpi
  use constants
  !> Variable with a short name that contains the double size 
  !> parameter MPI_DOUBLE_PRECISION
  INTEGER(isp)  :: mpidbl = MPI_DOUBLE_PRECISION
  !> Status parameter for MPI
  INTEGER(isp) :: status(MPI_STATUS_SIZE)
  ! Derived types (MPI exchange)
  !> Unused derived type 
  INTEGER(isp) :: derived_type_grid
  !> Unused derived type 
  INTEGER(isp) :: derived_subarray_grid
  !> Array containing a list of different derived types
  INTEGER(isp), DIMENSION(100) :: mpi_dtypes
  !> List of flags that tells if the corresponding derived type in 
  !> mpi_dtypes has been initialized
  LOGICAL(lp), DIMENSION(100) :: is_dtype_init = .TRUE.
END MODULE mpi_type_constants

!===============================================================================
!> Module for the communications
MODULE communications  !#do not parse
!===============================================================================
  use constants
  
  INTEGER(isp) :: reqperjxx(4)
  INTEGER(isp) :: reqperjxy(4)
  INTEGER(isp) :: reqperjxz(4)
  INTEGER(isp) :: reqperjyx(4)
  INTEGER(isp) :: reqperjyy(4)
  INTEGER(isp) :: reqperjyz(4)
  INTEGER(isp) :: reqperjzx(4)
  INTEGER(isp) :: reqperjzy(4)
  INTEGER(isp) :: reqperjzz(4)

  !> Structure that contains buffer arrays for particle 
  !> exchange between tiles.
  !> This structure is used by the subroutine 
  !> particle_bsc_openmp_reordering() in boundary.F90.
  TYPE part_com_buffer
      !> particle x position buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_x
      !dir$ attributes align:64 :: part_x
      !DIR ATTRIBUTES FASTMEM  :: part_x
      !> particle y position buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_y
      !dir$ attributes align:64 :: part_y
      !DIR ATTRIBUTES FASTMEM  :: part_y
      !> particle z position buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_z
      !dir$ attributes align:64 :: part_z
      !DIR ATTRIBUTES FASTMEM  :: part_z
      !> particle x momentum buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ux
      !dir$ attributes align:64 :: part_ux
      !DIR ATTRIBUTES FASTMEM  :: part_ux
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uy
      !dir$ attributes align:64 :: part_uy
      !DIR ATTRIBUTES FASTMEM  :: part_uy
      !> particle z momentum buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uz
      !dir$ attributes align:64 :: part_uz
      !DIR ATTRIBUTES FASTMEM  :: part_uz
      !> particle gamma factor buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:) :: part_gaminv
      !dir$ attributes align:64 :: part_gaminv
      !DIR ATTRIBUTES FASTMEM  :: part_gaminv
      !> particle weight buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
      !dir$ attributes align:64 :: pid
      !DIR ATTRIBUTES FASTMEM :: pid
      INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: boundid
      INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: bin_npart
      INTEGER(idp), ALLOCATABLE, DIMENSION(:) :: bin_pos
  END TYPE

  !> Structure that contains buffer arrays for particle 
  !> exchange between tiles and MPI domains.
  !> This structure is used by the subroutine 
  !> particle_bcs_tiles_and_mpi_3d() in boundary.F90.
  TYPE mpi_tile_buffer
      !> particle x position buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_x
      !dir$ attributes align:64 :: part_x
      !DIR ATTRIBUTES FASTMEM  :: part_x
      !> particle y position buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_y
      !dir$ attributes align:64 :: part_y
      !DIR ATTRIBUTES FASTMEM  :: part_y
      !> particle z position buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_z
      !dir$ attributes align:64 :: part_z
      !DIR ATTRIBUTES FASTMEM  :: part_z
      !> particle x momentum buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_ux
      !dir$ attributes align:64 :: part_ux
      !DIR ATTRIBUTES FASTMEM  :: part_ux
      !> particle y momentum buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_uy
      !dir$ attributes align:64 :: part_uy
      !DIR ATTRIBUTES FASTMEM  :: part_uy
      !> particle z momentum buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_uz
      !dir$ attributes align:64 :: part_uz
      !DIR ATTRIBUTES FASTMEM  :: part_uz
      !> particle gamma factor buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: part_gaminv
      !dir$ attributes align:64 :: part_gaminv
      !DIR ATTRIBUTES FASTMEM  :: part_gaminv
      !> particle weight buffer array
      REAL(num), ALLOCATABLE, DIMENSION(:,:) :: pid
      !dir$ attributes align:64 :: pid
      !DIR ATTRIBUTES FASTMEM  :: pid
      !> Number of particles to be exchanged in each direction
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
  !> start time
  REAL(num) :: startsim =0.0_num
  !> End time
  REAL(num) :: endsim =0.0_num
  !> Simulation time at the beginning of the iteration (used in submain.F90)
  REAL(num) :: startit
  !> Simulation time at the end of the iteration so that timeit-startit 
  !> is the iteration time (used in submain.F90)
  REAL(num) :: timeit
  !> Time spent in the pusher
  REAL(num) :: pushtime

  !> Output frequency of the field diagnostics
  INTEGER(idp) :: output_frequency = -1 !(Default is no output)
  !> First step for the field diagnostics
  INTEGER(idp) :: output_step_min = 0
  !> Last step for the field diagnostics  
  INTEGER(idp) :: output_step_max = 0

  ! output quantity flag (Default=False)
  !> Activation of the Ex electric field output
  INTEGER(KIND=4) :: c_output_ex = 0
  !> Activation of the Ey electric field output
  INTEGER(KIND=4) :: c_output_ey = 0
  !> Activation of the Ez electric field output
  INTEGER(KIND=4) :: c_output_ez = 0
  !> Activation of the Bx magnetic field output
  INTEGER(KIND=4) :: c_output_bx = 0
  !> Activation of the By magnetic field output
  INTEGER(KIND=4) :: c_output_by = 0
  !> Activation of the Bz magnetic field output
  INTEGER(KIND=4) :: c_output_bz = 0
  !> Activation of the Jx current field output
  INTEGER(KIND=4) :: c_output_jx = 0
  !> Activation of the Jy current field output
  INTEGER(KIND=4) :: c_output_jy = 0
  !> Activation of the Jz current field output
  INTEGER(KIND=4) :: c_output_jz = 0
  !> Activation of the density output
  INTEGER(KIND=4) :: c_output_rho = 0
  !> Activation of the electric field divergence output
  INTEGER(KIND=4) :: c_output_dive = 0

  ! File names for output dumps
  !> File name for the Ex electric field output
  CHARACTER(LEN=string_length) :: fileex   ='ex'
  !> File name for the Ey electric field output
  CHARACTER(LEN=string_length) :: fileey   ='ey'
  !> File name for the Ez electric field output
  CHARACTER(LEN=string_length) :: fileez   ='ez'
  !> File name for the Bx magnetic field output
  CHARACTER(LEN=string_length) :: filebx   ='bx'
  !> File name for the By magnetic field output
  CHARACTER(LEN=string_length) :: fileby   ='by'
  !> File name for the Bz magnetic field output
  CHARACTER(LEN=string_length) :: filebz   ='bz'
  !> File name for the Jx current output
  CHARACTER(LEN=string_length) :: filejx   ='jx'
  !> File name for the Jy current output
  CHARACTER(LEN=string_length) :: filejy   ='jy'
  !> File name for the Jz current output
  CHARACTER(LEN=string_length) :: filejz   ='jz'
  !> File name for the electric field divergence output
  CHARACTER(LEN=string_length) :: filedive ='dive'
  !> File name for the density output
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

  !> Number of temporal diags
  INTEGER(idp) :: temdiag_nb
  !> Number of particle temporal diagnostics
  INTEGER(idp) :: temdiag_nb_part
  !> Number of field temporal diagnostics
  INTEGER(idp) :: temdiag_nb_field
  !> Total number of temporal diagnostics
  INTEGER(idp) :: temdiag_totvalues
  !> Output frequency of the temporal diagnostics
  INTEGER(idp) :: temdiag_frequency
  !> Output format of the temporal diagnostics
  INTEGER(idp) :: temdiag_format
  !> Big array containing all the temporal diag at a given iteration
  REAL(num), dimension(:),allocatable :: temdiag_array

  ! Computation flags
  !> Flag true if the divergence of the electric field has been 
  !> calculated for the current iteration.
  LOGICAL(lp)  :: divE_computed

  ! Particle dump
  !> Flag true if the particle dumping is activated
  LOGICAL(lp) :: particle_dump_activated
  !> Number of species to dump
  INTEGER(lp) :: npdumps
  !> Structure for particle dumping
  TYPE particle_dump
    !> Corresponding species
    INTEGER(idp) :: ispecies
    !> Dumping period
    INTEGER(idp) :: diag_period
    !> Minimal particle position in x for dumping
    REAL(num) :: dump_x_min
    !> Maximal particle position in x for dumping
    REAL(num) :: dump_x_max
    REAL(num) :: dump_y_min
    REAL(num) :: dump_y_max
    REAL(num) :: dump_z_min
    REAL(num) :: dump_z_max
    REAL(num) :: dump_ux_min
    REAL(num) :: dump_ux_max
    REAL(num) :: dump_uy_min
    REAL(num) :: dump_uy_max
    REAL(num) :: dump_uz_min
    REAL(num) :: dump_uz_max
  END TYPE particle_dump

  !> Object for the particle dumping
  TYPE(particle_dump), ALLOCATABLE, TARGET, DIMENSION(:) :: particle_dumps
END MODULE output_data


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
  LOGICAL(lp)                        :: x_min_boundary, x_max_boundary
  LOGICAL(lp)                        :: y_min_boundary, y_max_boundary
  LOGICAL(lp)                        :: z_min_boundary, z_max_boundary
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
  LOGICAL(lp)                        :: l_axis_allocated=.FALSE.
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
  LOGICAL(lp)      :: sorting_verbose

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

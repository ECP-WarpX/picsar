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
! MODULES.F90
!
! This file contains several modules with parameters for different purposes.
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
! Contains shared data
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Module containing precision parameters (double, single, complex etc.)
!> for declaration of ints, floats etc.
! ________________________________________________________________________________________
MODULE PICSAR_precision
  !> Float precision
  INTEGER, PARAMETER :: num = 8
  !> Integer 4 byte precision
  INTEGER, PARAMETER :: isp = 4
  !> integer double precision
  INTEGER, PARAMETER :: idp = 8
  !> logical precision
  INTEGER, PARAMETER :: lp = 8
  !> Complex precision
  INTEGER, PARAMETER :: cpx = 8
END MODULE PICSAR_precision

! ________________________________________________________________________________________
!> @brief
!> Module containing the Picsar constant parameters
! ________________________________________________________________________________________
MODULE constants
  USE PICSAR_precision
  !> Electron mass
  REAL(num), PARAMETER :: emass   = 9.10938291e-31_num
  !> Proton mass
  REAL(num), PARAMETER :: pmass   = 1.6726231000000001e-27_num
  !> Electron charge
  REAL(num), PARAMETER :: echarge = 1.6021764620000001e-19_num
#if defined(LIBRARY)
  !> Speed of light in vacuum
  REAL(num), PARAMETER :: clight  = 1.0_num
  !> Magnetic constant
  REAL(num), PARAMETER :: mu0     = 1.0_num
  !> Vacuum permeability
  REAL(num), PARAMETER :: eps0    = 1.0_num
  REAL(num), PARAMETER :: imu0    = 1.0_num
  !> The famous pi value
#else
  !> Speed of light in vacuum
  REAL(num), PARAMETER :: clight  = 2.99792458e8_num
  !> Magnetic constant
  REAL(num), PARAMETER :: mu0     = 1.2566370614359173e-06_num
  !> Vacuum permeability
  REAL(num), PARAMETER :: eps0    = 8.854187817620389e-12_num
  REAL(num), PARAMETER :: imu0    = 795774.715459_num
#endif
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

! ________________________________________________________________________________________
!> @brief
!> Module containing useful pre-computed parameters for some subroutines
! ________________________________________________________________________________________
MODULE precomputed
  USE PICSAR_precision
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
  REAL(num) :: dts2dx, dts2dy, dts2dz
  REAL(num) :: dtsdx0, dtsdy0, dtsdz0
  REAL(num) :: dxs2, dys2, dzs2
  REAL(num) :: clightsq
END MODULE precomputed

! ________________________________________________________________________________________
!> @brief
!> Module containing parameters and data structures for the fields
! ________________________________________________________________________________________
MODULE fields
  USE PICSAR_precision
  USE constants
  !> Flag: interpolation at a lower order for the field gathering
  LOGICAL(lp) :: l_lower_order_in_v
  !> Flag: use of nodal grids
  LOGICAL(lp) :: l_nodalgrid
  !> Flag: use of PSAOTD spectral solver
  LOGICAL(lp) :: l_spectral
  !> Flag: use psatd with multiply_mat_vector_routine (not suited for prod)
  LOGICAL(lp) :: g_spectral = .FALSE.
  !> Flag: use of staggered grid
  LOGICAL(lp) :: l_staggered = .TRUE.
  !> Flag: this flag needs a description, used in field gathering routines
  LOGICAL(lp) :: l4symtry
  INTEGER(idp):: nxs=0, nys=0, nzs=0
  !> order in x of the FDTD Maxwell solver
  INTEGER(idp):: norderx
  !> order in y of the FDTD Maxwell solver
  INTEGER(idp):: nordery
  !> order in z of the FDTD Maxwell solver
  INTEGER(idp):: norderz
  !> n_pml in x direction 
  INTEGER(idp):: nx_pml
  !> n_pml in y direction 
  INTEGER(idp):: ny_pml
  !> n_pml in z direction 
  INTEGER(idp):: nz_pml
  !> shift_pml :: number of guardcells forced to 0 when using PMLS.
  !> Available option only for local solver.
  !> When using hybrid solver all guardcells are forced to 0 because the  guardcells 
  !> of the local fields (ex, ey ....) are not filled from the FFT-grid fields (ex_r,ey_r ...)
  INTEGER(idp) :: shift_x_pml ,shift_y_pml, shift_z_pml
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
  REAL(num), POINTER, DIMENSION(:, :, :) :: ex
  !> MPI-domain electric field grid in y
  REAL(num), POINTER, DIMENSION(:, :, :) :: ey
  !> MPI-domain electric field grid in z
  REAL(num), POINTER, DIMENSION(:, :, :) :: ez
  !> MPI-domain magnetic field grid in x
  REAL(num), POINTER, DIMENSION(:, :, :) :: bx
  !> MPI-domain magnetic field grid in y
  REAL(num), POINTER, DIMENSION(:, :, :) :: by
  !> MPI-domain magnetic field grid in z
  REAL(num), POINTER, DIMENSION(:, :, :) :: bz
  !> MPI-domain current grid in x
  REAL(num), POINTER, DIMENSION(:, :, :) :: jx
  !> MPI-domain current grid in y
  REAL(num), POINTER, DIMENSION(:, :, :) :: jy
  !> MPI-domain current grid in z
  REAL(num), POINTER, DIMENSION(:, :, :) :: jz
  !> MPI-domain electric field grid in x (auxiliary array for gather to particles)
  REAL(num), POINTER, DIMENSION(:, :, :) :: ex_p
  !> MPI-domain electric field grid in y (auxiliary array for gather to particles)
  REAL(num), POINTER, DIMENSION(:, :, :) :: ey_p
  !> MPI-domain electric field grid in z (auxiliary array for gather to particles)
  REAL(num), POINTER, DIMENSION(:, :, :) :: ez_p
  !> MPI-domain magnetic field grid in x (auxiliary array for gather to particles)
  REAL(num), POINTER, DIMENSION(:, :, :) :: bx_p
  !> MPI-domain magnetic field grid in y (auxiliary array for gather to particles)
  REAL(num), POINTER, DIMENSION(:, :, :) :: by_p
  !> MPI-domain magnetic field grid in z (auxiliary array for gather to particles)
  REAL(num), POINTER, DIMENSION(:, :, :) :: bz_p
  !> MPI-domain current grid in x
  !> MPI-domain electric field grid in x
  REAL(num), POINTER, DIMENSION(:, :, :) :: ex_r
  !> MPI-domain electric field grid in y
  REAL(num), POINTER, DIMENSION(:, :, :) :: ey_r
  !> MPI-domain electric field grid in z
  REAL(num), POINTER, DIMENSION(:, :, :) :: ez_r
  !> MPI-domain magnetic field grid in x
  REAL(num), POINTER, DIMENSION(:, :, :) :: bx_r
  !> MPI-domain magnetic field grid in y
  REAL(num), POINTER, DIMENSION(:, :, :) :: by_r
  !> MPI-domain magnetic field grid in z
  REAL(num), POINTER, DIMENSION(:, :, :) :: bz_r
  !> MPI-domain current grid in x
  REAL(num), POINTER, DIMENSION(:, :, :) :: jx_r
  !> MPI-domain current grid in y
  REAL(num), POINTER, DIMENSION(:, :, :) :: jy_r
  !> MPI-domain current grid in z
  REAL(num), POINTER, DIMENSION(:, :, :) :: jz_r
  !> MPI-domain current grid in z - Fourier space
  REAL(num), POINTER, DIMENSION(:, :, :) :: rho_r
  !> MPI-domain current grid in z - Fourier space
  REAL(num), POINTER, DIMENSION(:, :, :) :: rhoold_r

  !> MPI-domain splitted EM fields for PML
  REAL(num), POINTER, DIMENSIOn(:,:,:) :: exy_r, exz_r, eyx_r, eyz_r, ezx_r,&
  ezy_r, bxy_r, bxz_r, byx_r, byz_r, bzx_r, bzy_r
  REAL(num) , POINTER, DIMENSION(:,:,:) :: exy,exz,eyx,eyz,ezx,ezy, &
        bxy,bxz,byx,byz,bzx,bzy
  REAL(num) , POINTER, DIMENSION(:) :: sigma_x_e, sigma_y_e, sigma_z_e, &
        sigma_x_b, sigma_y_b, sigma_z_b
  !> MPI-domain electric field grid in x - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: exf
  !> MPI-domain electric field grid in y - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: eyf
  !> MPI-domain electric field grid in z - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: ezf
  !> MPI-domain magnetic field grid in x - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: bxf
  !> MPI-domain magnetic field grid in y - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: byf
  !> MPI-domain magnetic field grid in z - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: bzf
  !> MPI-domain current grid in x - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: jxf
  !> MPI-domain current grid in y - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: jyf
  !> MPI-domain current grid in z - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: jzf
  !> MPI-domain current grid in z - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: rhof
  !> MPI-domain current grid in z - Fourier space
  COMPLEX(cpx), POINTER, DIMENSION(:, :, :) :: rhooldf
  !> Fonberg coefficients in x
  REAL(num), POINTER, DIMENSION(:) :: xcoeffs
  !> Fonberg coefficients in y
  REAL(num), POINTER, DIMENSION(:) :: ycoeffs
  !> Fonberg coefficients in z
  REAL(num), POINTER, DIMENSION(:) :: zcoeffs

  !> Electric energy withi mpi domain
  REAL(num)                        :: electro_energy_mpi
  !> Magnetic energy withi mpi domain
  REAL(num)                        :: magnetic_energy_mpi
  !> ElectroMagnetic energy withi mpi domain
  REAL(num)                        :: electromagn_energy_mpi
  !> Total Electric energy 
  REAL(num)                        :: electro_energy_total
  !> Total Magnetic energy 
  REAL(num)                        :: magneto_energy_total
  !> Total ElectroMagnetic energy
  REAL(num)                        :: electromagn_energy_total
END MODULE fields

! ________________________________________________________________________________________
!> @brief
!> Module containing the current/charge tile data structure.
! ________________________________________________________________________________________
MODULE grid_tilemodule!#do not parse
  USE PICSAR_precision
  USE constants
  !> This object contains 3D field grids for one tile 
  TYPE grid_tile
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: arr1 ! For X current component 
                                                       ! or charge 
    !> Tile Current grid in y
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: arr2 ! For Y current component 
    !> Tile Current grid in z
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: arr3 ! For Z current component
    ! We declare arrays aligned for vectorization efficiency.
    ! These directives are only understood by the Intel compiler.
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: arr1
    !dir$ attributes align:64 :: arr2
    !dir$ attributes align:64 :: arr3
#endif
  END TYPE

  !> This array contains a list of current grid tiles
  !> This array is local to each MPI domain.
  !> This structure is used to avoid contentions during the parallel threaded 
  !> reductions of currents by threads in the large current array (used in the 
  !> the Maxwell solver). In particular, to avoid contention, copies of guard cells 
  !> from grid tiles to the large current arrays are done alternatively in X, Y, Z 
  !> directions with implicit synchronization between each direction. This imposes to 
  !> store grid tile arrays in shared memory via this structure aofgrid_tiles. 
  !> This moderately increases memory footprint but allows 
  !> for very good parallel efficiency of the current deposition. 
  !> N.B: for field gathering fields are copied from the large field arrays (used 
  !> in the Maxwell solver) to private grid tile arrays to each thread before the 
  !> field gathering step on the particles. As this step does not involve any contention, 
  !> (copy operation by threads) these private arrays can thus
  !> be allocated/de-allocated when needed and do not require to be saved in 
  !> aofgrid_tiles (therefore saving a significant amount of memory). 
  TYPE(grid_tile), ALLOCATABLE, TARGET, DIMENSION(:, :, :) :: aofgrid_tiles

END MODULE grid_tilemodule

! ________________________________________________________________________________________
!> @brief
!> Module containing the Fortran object descriptor representing a particle tile.
!> Also see tiling.F90 for the definition of the tile properties.
! ________________________________________________________________________________________
MODULE particle_tilemodule!#do not parse
  USE PICSAR_precision
  USE constants
  !> Object that contains tile particle arrays and particle tile properties.
  TYPE particle_tile
    !> Flag:  current tile kin energy
    REAL(num)   :: kin_energy_tile
    !> Flag: tile arrays are allocated
    LOGICAL(lp) :: l_arrays_allocated= .FALSE.
    !> Current number of particles in tile
    INTEGER(idp), DIMENSION(1) :: np_tile= (/0_idp/)
    !> Max number of particles per tile: size of the arrays
    INTEGER(idp) :: npmax_tile = 0_idp
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
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: pid
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
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

! ________________________________________________________________________________________
!> Module used for particle exchanges in MPI routines
! ________________________________________________________________________________________
MODULE buff_exchange_part!#do not parse
USE PICSAR_precision
USE constants
  TYPE buff_part
    INTEGER(idp) :: nbuff ! % curent size of buffer array
    INTEGER(idp) :: ibuff ! % curent position in buffer array
    REAL(num), DIMENSION(:), ALLOCATABLE  :: buff_arr
  END TYPE buff_part
END MODULE buff_exchange_part

! ________________________________________________________________________________________
!> Module defining particle_antenna type
! ________________________________________________________________________________________
MODULE antenna!#do not parse
  USE PICSAR_precision
  USE constants
  TYPE particle_antenna
    REAL(num)         ::  spot_x
    REAL(num)         ::  spot_y
    REAL(num)         ::  spot_z
    REAL(num), DIMENSION(3) ::  polvector1, polvector2, vector
    REAL(num)         :: laser_a_1!laser particle max_v_1 at focus (in clight unit)
    REAL(num)         :: laser_a_2!laser particle max_v_2 at focus (in clight unit)
    REAL(num)         :: Emax_laser_1
    REAL(num)         :: Emax_laser_2
    REAL(num)         :: Emax
    REAL(num)         :: laser_w0!laser waist at focus
    REAL(num)         :: inv_w02! 1./w0**2
    COMPLEX(cpx)      :: q_z! complex curv on the plan
    REAL(num)         :: laser_ctau! length of the pulse --->
    ! ---> (length from the peak to 1/e*pick= c*time_duration_of_the_pulse)
    REAL(num)         :: laser_tau! time duration of the pulse
    REAL(num)         :: t_peak
    REAL(num)         :: zr! rayleigh length of the laser
    REAL(num)         :: inv_zr!1/zr
    REAL(num)         :: polangle! phase shift between laser along povector2 -->
    ! --> and polvector1
    REAL(num)         :: lambda_laser
    REAL(num)         :: k0_laser
    INTEGER(idp)      :: temporal_order
    INTEGER(idp)      :: time_window! 0 for Gaussian 1 Hanning Window
  END TYPE particle_antenna
END MODULE antenna

! ________________________________________________________________________________________
!> @brief
!> Module containing the Fortran object descriptor representing a particle species
! ________________________________________________________________________________________
MODULE particle_speciesmodule!#do not parse
  USE particle_tilemodule
  USE PICSAR_precision
  USE constants
  USE antenna
  REAL(num)   :: kin_energy_mpi
  REAL(num)   :: kin_energy_total

  !> Fortran object representing a particle species
  TYPE particle_species
    ! Species kinetic energy
    REAL(num)   :: kin_energy_sp
    ! Attributes of particle species object
    !> Particle antenna flag (.FALSE. by default)
    LOGICAL(lp) :: is_antenna = .FALSE.
    !> Antenna params (init if is_antenna true)
    TYPE(particle_antenna) :: antenna_params
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
    !> Flag indicating if this particle species deposit current/charge on
    !> the grid (useful for test particles). Default is TRUE
    LOGICAL(lp)   :: ldodepos =.TRUE.
    !> Flag indicating if current particle species is freezed (no push, no field gathering)
    !> To completely stop particle routines set ldodepos to False and lfreeze to True  
    LOGICAL(lp)   :: lfreeze  =.FALSE.
    ! For some stupid reason, cannot use ALLOCATABLE in derived types
    ! in Fortran 90 - Need to use POINTER instead
    !> List of tiles (of objects particle_tile) in the MPI domain for the
    !> particles of this species.
    TYPE(particle_tile), DIMENSION(:, :, :), ALLOCATABLE :: array_of_tiles
    !> Array indicating if a tile has been reallocated.
    !> Used for interfacing WARP and PXR
    INTEGER(idp), DIMENSION(:, :, :), ALLOCATABLE :: are_tiles_reallocated
  END TYPE
END MODULE particle_speciesmodule

! ________________________________________________________________________________________
!> @brief
!> Module for the tile parameters
! ________________________________________________________________________________________
MODULE tile_params
  ! # of particle tiles in each dimension
  USE PICSAR_precision
  USE constants
  !> Number of tile in the x direction
  INTEGER(idp) :: ntilex
  !> Number of tile in the y direction
  INTEGER(idp) :: ntiley
  !> Number of tile in the z direction
  INTEGER(idp) :: ntilez
END MODULE tile_params


! ________________________________________________________________________________________
!> @brief
!> Module containing useful properties for the particles
! ________________________________________________________________________________________
MODULE particle_properties
  USE PICSAR_precision
  USE constants
  !> Number of elements per particle in the pid particle array
  !> Default is 1 i.e only particle weights are recorded
  INTEGER(idp)  :: npid=3
  !> Index in pid array corresponding to particle weights
  !> Beware: default is wpid=1. Use same in WARP when coupling WARP+PXR
  INTEGER(idp), PARAMETER :: wpid=1
  !> Index in pid array corresponding to particle ids
  INTEGER(idp)  :: ssnpid
  !> Index in pid array corresponding to old x positions of particles
  INTEGER(idp)  :: xoldpid
  !> Index in pid array corresponding to old y positions of particles
  INTEGER(idp) :: yoldpid
  !> Index in pid array corresponding to old x positions of particles
  INTEGER(idp) :: zoldpid
  !> Index in pid array corresponding to old x momentum of particles
  INTEGER(idp) :: uxoldpid
  !> Index in pid array corresponding to old y momentum of particles
  INTEGER(idp) :: uyoldpid
  !> Index in pid array corresponding to old z momentum of particles
  INTEGER(idp) :: uzoldpid 
  !> Index in pid array corresponding to old x positions of particles
  INTEGER(idp)  :: exoldpid
  !> Index in pid array corresponding to old y positions of particles
  INTEGER(idp) :: eyoldpid
  !> Index in pid array corresponding to old x positions of particles
  INTEGER(idp) :: ezoldpid
  !> Index in pid array corresponding to old x momentum of particles
  INTEGER(idp) :: bxoldpid
  !> Index in pid array corresponding to old y momentum of particles
  INTEGER(idp) :: byoldpid
  !> Index in pid array corresponding to old z momentum of particles
  INTEGER(idp) :: bzoldpid  
  !> This flag seems to be unused
  LOGICAL(lp) :: l_initongrid = .FALSE.
  !> Flag to activate the use of weight for the particles
  LOGICAL(lp) :: l_particles_weight = .FALSE.
  !> Particle pusher type (0: Boris, 1: Vay, Default: 0)
  INTEGER(idp) :: particle_pusher = 0
  !> Particle initial distribution
  INTEGER(idp) :: pdistr
  !> Number of species
  INTEGER(idp) :: nspecies = 0 
  !> total number of particles (all species, all subdomains -> useful for stat)
  INTEGER(idp) :: ntot
  !> Max number of particle species
  INTEGER(idp) :: nspecies_max = 40_idp
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
  !> Flag for the allocation of the grid tile arrays
  LOGICAL(lp) :: l_aofgrid_tiles_allocated=.FALSE.
  LOGICAL(lp) :: l_aofgrid_tiles_array_allocated=.FALSE.
  !> Flag for plasma init/push
  LOGICAL(lp) :: l_plasma = .TRUE.
END MODULE particle_properties

! ________________________________________________________________________________________
!> Module containing the array of species
! ________________________________________________________________________________________
MODULE particles!#do not parse
  USE PICSAR_precision
  USE constants
  USE tile_params
  USE particle_tilemodule
  USE particle_speciesmodule
  USE particle_properties
  USE grid_tilemodule

  !> Array of particle species objects
  TYPE(particle_species), ALLOCATABLE, TARGET, DIMENSION(:) :: species_parray

END MODULE particles

! ________________________________________________________________________________________
!> Module containing useful configuration and simulation parameters
! ________________________________________________________________________________________
MODULE params
  USE PICSAR_precision
  USE constants
  !> iteration number
  INTEGER(idp)         :: it=0_idp
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
  INTEGER(idp) :: currdepo = 0
  !> Charge deposition method
  INTEGER(idp) :: rhodepo = 0
  !> Field gathering method
  INTEGER(idp) :: fieldgathe = 0
  !> Type of comm routine to use for particles
  INTEGER(idp) :: partcom = 0
  !> Field gathering + part. pusher seperated flag
  INTEGER(idp) :: fg_p_pp_separated = 0
  !> Vector size for the current deposition
  INTEGER(idp) :: lvec_curr_depo = 16
  !> Vector size for the charge deposition
  INTEGER(idp) :: lvec_charge_depo = 16
  !> Vector size for the field gathering
  INTEGER(idp) :: lvec_fieldgathe = 16
  !> MPI buffer size
  INTEGER(idp) :: mpi_buf_size

END MODULE params

! ________________________________________________________________________________________
!> Module containing MPI parameters
! ________________________________________________________________________________________
MODULE mpi_type_constants!#do not parse
  use mpi
  USE PICSAR_precision
  USE constants
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

! ________________________________________________________________________________________
!> Module for the communications
! ________________________________________________________________________________________
MODULE communications!#do not parse
  USE PICSAR_precision
  USE constants
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
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_x
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_x
    !> particle y position buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_y
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_y
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_y
    !> particle z position buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_z
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_z
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_z
    !> particle x momentum buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_ux
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_ux
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_ux
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uy
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_uy
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_uy
    !> particle z momentum buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_uz
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_uz
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_uz
    !> particle gamma factor buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:) :: part_gaminv
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_gaminv
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_gaminv
    !> particle weight buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: pid
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: pid
#endif
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
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_x
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_x
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_x
    !> particle y position buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_y
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_y
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_y
    !> particle z position buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_z
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_z
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_z
    !> particle x momentum buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_ux
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_ux
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_ux
    !> particle y momentum buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_uy
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_uy
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_uy
    !> particle z momentum buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_uz
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_uz
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_uz
    !> particle gamma factor buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :) :: part_gaminv
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: part_gaminv
#endif
    !DIR ATTRIBUTES FASTMEM  :: part_gaminv
    !> particle weight buffer array
    REAL(num), ALLOCATABLE, DIMENSION(:, :, :) :: pid
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
    !dir$ attributes align:64 :: pid
#endif
    !DIR ATTRIBUTES FASTMEM  :: pid
    !> Number of particles to be exchanged in each direction
    INTEGER(idp), dimension(27) :: npart
  END TYPE
END MODULE communications

! ________________________________________________________________________________________
!> Module for the time statistics
! ________________________________________________________________________________________
MODULE time_stat!#do not parse
  USE PICSAR_precision
  USE constants
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
  REAL(num), dimension(26)               :: localtimes
  REAL(num), DIMENSION(26)               :: mintimes, init_mintimes
  REAL(num), DIMENSION(26)               :: maxtimes, init_maxtimes
  REAL(num), DIMENSION(26)               :: avetimes, init_avetimes
  !> Buffer for the output
  REAL(num), DIMENSION(:, :), POINTER     :: buffer_timestat
  INTEGER(idp)                           :: itimestat
  !> Number of entries in the buffer
  INTEGER(idp)                           :: nbuffertimestat
END MODULE time_stat

! ________________________________________________________________________________________
!> Module for the outputs
! ________________________________________________________________________________________
MODULE output_data!#do not parse
  USE PICSAR_precision
  USE constants

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
  INTEGER(idp) :: output_frequency = -1!(Default is no output)
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
  !> Activation of div J field divergence output
  INTEGER(KIND=4) :: c_output_divj = 0
  !> Activation of div B field divergence output
  INTEGER(KIND=4) :: c_output_divb = 0

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
  !> File name for the current field divergence output
  CHARACTER(LEN=string_length) :: filedivj ='divj'
  !> File name for the magnetic field divergence output
  CHARACTER(LEN=string_length) :: filedivb ='divb'


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
  REAL(num), dimension(:), allocatable :: temdiag_array

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
    !> Minimal particle position in y for dumping
    REAL(num) :: dump_y_min
    !> Maximal particle position in y for dumping
    REAL(num) :: dump_y_max
    !> Minimal particle position in z for dumping
    REAL(num) :: dump_z_min
    !> Maximal particle position in z for dumping
    REAL(num) :: dump_z_max
    !> Minimal particle momentum in x for dumping
    REAL(num) :: dump_ux_min
    !> Maximal particle momentum in x for dumping
    REAL(num) :: dump_ux_max
    !> Minimal particle momentum in y for dumping
    REAL(num) :: dump_uy_min
    !> Maximal particle momentum in y for dumping
    REAL(num) :: dump_uy_max
    !> Minimal particle momentum in z for dumping
    REAL(num) :: dump_uz_min
    !> Maximal particle momentum in z for dumping
    REAL(num) :: dump_uz_max
  END TYPE particle_dump

  !> Object for the particle dumping
  TYPE(particle_dump), ALLOCATABLE, TARGET, DIMENSION(:) :: particle_dumps
END MODULE output_data

!MODULE FOR GROUP params
#if defined(FFTW)
MODULE group_parameters !#do not parse
  USE mpi_type_constants
  USE picsar_precision

  !> group sizes of of all groups
  INTEGER(idp), DIMENSION(3) :: group_sizes
  !> To which group this mpi task belongs
  INTEGER(idp)    ::  which_group
  !> x y z coordinats of the group
  INTEGER(idp)    :: x_group_coords, y_group_coords, z_group_coords
  !> coordinates of each proc in its group
  INTEGER(isp)  ::  group_coordinates(3)
  !> local_size and local rank
  INTEGER(isp)    :: local_size, local_rank
  !> MPI_GROUP associated to mpi_comm_world
  INTEGER(isp)    :: mpi_world_group
  !> ARRAY of  MPI_GROUP associated to each mpi task (!= mpi_group_null or
  !mpi_comm_null if and  only if i == which group + 1
  INTEGER(isp) ,DIMENSION(:), POINTER :: mpi_group_id, mpi_comm_group_id
  !>  MPI_COMM for local roots group and MPI_GROUP for local  roots and roots
  !ranks in the mpi_root_comm
  INTEGER(isp)  :: mpi_root_comm, mpi_root_group, root_rank, root_size
  !> Ordered comm world (for computing group splitting)
  INTEGER(isp)  :: mpi_ordered_comm_world
  !> Field cell  sizes in groups without guardcells
  INTEGER(idp)  :: nx_group_global, ny_group_global, nz_group_global
  INTEGER(idp) , DIMENSION(:), POINTER ::  nz_group_global_array, ny_group_global_array, nx_group_global_array
  !> Field grid sizes in groups whithout guardcells
  INTEGER(idp)  :: nx_group_global_grid, ny_group_global_grid, nz_group_global_grid
  !> Field cell  sizes in groups with guardcells
  INTEGER(idp)  :: nx_group, ny_group, nz_group
  !> Field grid sizes in groups with guardcells
  INTEGER(idp)  :: nx_group_grid, ny_group_grid, nz_group_grid

  !> Nz grid min max group index
  INTEGER(idp)  ::   nz_grid_min_grp, nz_grid_max_grp, nz_grid_grp
  !> This flag is true if MPI task is on the edge of its group (so need
  !additional comm
  LOGICAL(lp)  ::  is_on_boundary_group_z = .FALSE.
  !> This flag is true if the MPI rank is at the inferior z group boundary
  LOGICAL(lp)  :: group_z_min_boundary = .FALSE.
  !> This flag is true if the MPI rank is at the superior z group boundary
  LOGICAL(lp)  :: group_z_max_boundary = .FALSE.
  LOGICAL(lp)  :: group_y_min_boundary = .FALSE.
  LOGICAL(lp)  :: group_y_max_boundary = .FALSE. 
  LOGICAL(lp)  :: is_on_boundary_group_y = .FALSE.
  !> minimum and maximum cell numbers in each group :
  INTEGER(idp), DIMENSION(:), POINTER  :: cell_z_min_group, cell_z_max_group
  !> physical limits of group domains
  REAL(num)                                 :: z_min_group, z_max_group
  REAL(num)                                 :: y_min_group, y_max_group
  REAL(num)                                 :: x_min_group, x_max_group

  REAL(num)                                  :: z_min_local_lb, z_max_local_lb 
  INTEGER(idp)                               :: nz_global_grid_min_lb , nz_global_grid_max_lb

  !> Cell domain for load balancing general case (taking into account
  !--guardcells)
  INTEGER(idp)  , DIMENSION(:) , POINTER :: cell_z_min_g, cell_z_max_g, &
  cell_x_min_g
  INTEGER(idp)  , DIMENSION(:) , POINTER :: cell_y_min_g, cell_y_max_g, &
  cell_x_max_g
  INTEGER(idp)  , DIMENSION(:) , POINTER :: size_exchanges_l2g_recv_z, size_exchanges_g2l_recv_z
  INTEGER(idp)  , DIMENSION(:) , POINTER :: size_exchanges_g2l_send_z,size_exchanges_l2g_send_z 
  INTEGER(idp)  , DIMENSION(:) , POINTER :: g_first_cell_to_recv_z,l_first_cell_to_recv_z
  INTEGER(idp)  , DIMENSION(:) , POINTER :: g_first_cell_to_send_z,l_first_cell_to_send_z

  INTEGER(idp)  , DIMENSION(:,:) , POINTER :: size_exchanges_l2g_recv,size_exchanges_g2l_recv
  INTEGER(idp)  , DIMENSION(:,:) , POINTER :: size_exchanges_g2l_send,size_exchanges_l2g_send


  INTEGER(idp)  , DIMENSION(:) , POINTER :: size_exchanges_l2g_recv_y,size_exchanges_g2l_recv_y
  INTEGER(idp)  , DIMENSION(:) , POINTER :: size_exchanges_g2l_send_y,size_exchanges_l2g_send_y
  INTEGER(idp)  , DIMENSION(:) , POINTER :: g_first_cell_to_recv_y,l_first_cell_to_recv_y
  INTEGER(idp)  , DIMENSION(:) , POINTER :: g_first_cell_to_send_y,l_first_cell_to_send_y


  !> TYPE IN WHICH ex_r will be recieving
  INTEGER(isp)  , DIMENSION(:) , POINTER :: recv_type_g   
  !> TYPE IN WHICH ex_r will be sending
  INTEGER(isp)  , DIMENSION(:) , POINTER :: send_type_g
  !> TYPE IN WHICH ex will be recieving
  INTEGER(isp)  , DIMENSION(:) , POINTER :: recv_type_l
  !> TYPE IN WHICH ex will be sending
  INTEGER(isp)  , DIMENSION(:) , POINTER :: send_type_l
  INTEGER(isp)  , DIMENSION(:) , POINTER :: array_of_ranks_to_send_to,array_of_ranks_to_send_to_l2g,    &
        array_of_ranks_to_send_to_g2l
  INTEGER(isp)  , DIMENSION(:) , POINTER :: array_of_ranks_to_recv_from,array_of_ranks_to_recv_from_l2g,&
  array_of_ranks_to_recv_from_g2l
  INTEGER(isp)  , DIMENSION(:) , POINTER :: requests_l2g, requests_g2l
  !> Work_array_g2l and Work_array_l2g contain are arrays of sizes nb_comms_g2l and nb_comms_g2l
  !> respectively.
  !> each cell of arrays encode the localization of array with which
  !> communication is performed 
  INTEGER(idp)  , DIMENSION(:) , POINTER :: work_array_g2l, work_array_l2g
  !> Nb_comms_g2l and nb_comms_l2g are equal to the number of send + recv calls
  !done by each mpi in mpi comms group during l->g and g->l communications respectively
  INTEGER(idp)  :: nb_comms_g2l,nb_comms_l2g
  INTEGER(isp) , DIMENSION(3) :: p3d_istart, p3d_iend , p3d_fstart,p3d_fend, p3d_fsize, p3d_isize

  !> Tells if current group is on z axis domain boundary
  LOGICAL(lp)    :: is_group_z_boundary_max, is_group_z_boundary_min
  !> Tells if current group is on y axis domain boundary
  LOGICAL(lp)    :: is_group_y_boundary_max, is_group_y_boundary_min
  !> Tells if current group is on x axis domain boundary
  LOGICAL(lp)    :: is_group_x_boundary_max, is_group_x_boundary_min

END MODULE group_parameters

#endif


! ________________________________________________________________________________________
!> Module for the data shared with Python.
! ________________________________________________________________________________________
MODULE shared_data
  use mpi_type_constants
  USE output_data
  !----------------------------------------------------------------------------
  ! MPI subdomain data
  !----------------------------------------------------------------------------
  !> FFTW distributed
  LOGICAL(lp)  :: fftw_with_mpi, fftw_mpi_transpose, fftw_threads_ok, fftw_hybrid
  !> Number of groups (this is a parameter in the input file)
  INTEGER(idp)    ::  nb_group
  !> Number of groups in each direction
  INTEGER(idp) :: nb_group_z, nb_group_y, nb_group_x
  !> Group guard cells in : (only nzg_group  and nyg_group are relevant for now)
  INTEGER(idp)  :: nzg_group, nyg_group, nxg_group
  LOGICAL(lp)   :: p3dfft_flag     = .FALSE.
  LOGICAL(lp)   :: p3dfft_stride   = .FALSE.
  LOGICAL(lp)   :: absorbing_bcs   =   .FALSE.
  LOGICAL(lp)   :: absorbing_bcs_x = .FALSE.
  LOGICAL(lp)   :: absorbing_bcs_y = .FALSE.
  LOGICAL(lp)   :: absorbing_bcs_z = .FALSE.
  LOGICAL(lp)  :: fftw_plan_measure=.TRUE.
  !> First and last indexes of real data in group (only z is relevant for now)
  INTEGER(idp)  ::   iz_min_r, iz_max_r, iy_min_r, iy_max_r, ix_min_r, ix_max_r

  !> Error code for MPI
  INTEGER(isp) :: errcode
  !> Variable used by MPI
  INTEGER(isp) :: provided
  !> Communicator used by MPI
  INTEGER(isp) :: comm
  !> Tag variable used by MPI
  INTEGER(isp) :: tag
  !> MPI process rank
  INTEGER(idp) :: rank
  !> MPI process coordinates in a Cartesian topology
  INTEGER(isp) :: coordinates(3)
  !> MI process neighbors in a Cartesian topology
  INTEGER(idp) :: neighbour(-1:1, -1:1, -1:1)
  !> MPI process x coordinate, equivalent of x_coords = coordinates(c_ndims)
  !> with a Cartesian topology
  INTEGER(idp) :: x_coords
  !> Minimal MPI process topology index in x
  INTEGER(idp) :: proc_x_min
  !> Maximal MPI process topology index in x
  INTEGER(idp) :: proc_x_max
  !> MPI process y coordinate, equivalent of x_coords = coordinates(c_ndims-1)
  !> with a Cartesian topology
  INTEGER(idp) :: y_coords
  !> Minimal MPI process topology index in y
  INTEGER(idp) :: proc_y_min
  !> Maximal MPI process topology index in y
  INTEGER(idp) :: proc_y_max
  !> MPI process y coordinate, equivalent of x_coords = coordinates(c_ndims-1)
  !> with a Cartesian topology
  INTEGER(idp) :: z_coords
  !> Minimal MPI process topology index in z
  INTEGER(idp) :: proc_z_min
  !> Maximal MPI process topology index in z
  INTEGER(idp) :: proc_z_max
  !> Number of MPI processes
  INTEGER(idp) :: nproc
  !> Number of MPI processes in x
  INTEGER(idp) :: nprocx
  !> Number of MPI processes in y
  INTEGER(idp) :: nprocy
  !> Number of MPI processes in z
  INTEGER(idp) :: nprocz
  !> This parameter is not used
  INTEGER(isp) :: nprocdir(3)
  ! Boundary data
  !> This flag is true if the MPI rank is at the inferior x grid boundary
  LOGICAL(lp)  :: x_min_boundary
  !> This flag is true if the MPI rank is at the superior x grid boundary
  LOGICAL(lp)  :: x_max_boundary
  !> This flag is true if the MPI rank is at the inferior y grid boundary
  LOGICAL(lp)  :: y_min_boundary
  !> This flag is true if the MPI rank is at the superior y grid boundary
  LOGICAL(lp)  :: y_max_boundary
  !> This flag is true if the MPI rank is at the inferior z grid boundary
  LOGICAL(lp)  :: z_min_boundary
  !> This flag is true if the MPI rank is at the superior z grid boundary
  LOGICAL(lp)  :: z_max_boundary
  !> This flag is true if the MPI rank is at the inferior x part boundary
  LOGICAL(lp)  :: x_min_boundary_part
  !> This flag is true if the MPI rank is at the superior x part boundary
  LOGICAL(lp)  :: x_max_boundary_part
  !> This flag is true if the MPI rank is at the inferior y part boundary
  LOGICAL(lp)  :: y_min_boundary_part
  !> This flag is true if the MPI rank is at the superior y part boundary
  LOGICAL(lp)  :: y_max_boundary_part
  !> This flag is true if the MPI rank is at the inferior z part boundary
  LOGICAL(lp)  :: z_min_boundary_part
  !> This flag is true if the MPI rank is at the superior z part boundary
  LOGICAL(lp)  :: z_max_boundary_part
  !> Type of particle boundary condition at the inferior x boundary
  INTEGER(idp) :: pbound_x_min
  !> Type of boundary condition at the superior x boundary
  INTEGER(idp) :: pbound_x_max
  !> Type of boundary condition at the inferior y boundary
  INTEGER(idp) :: pbound_y_min
  !> Type of boundary condition at the superior y boundary
  INTEGER(idp) :: pbound_y_max
  !> Type of boundary condition at the inferior z boundary
  INTEGER(idp) :: pbound_z_min
  !> Type of boundary condition at the superior z boundary
  INTEGER(idp) :: pbound_z_max
  ! The location of the processors
  !> Minimum cell number in x for each MPI process
  INTEGER(idp), DIMENSION(:), POINTER :: cell_x_min
  !> Maximal cell number in x for each MPI process
  INTEGER(idp), DIMENSION(:), POINTER :: cell_x_max
  !> Minimum cell number in y for each MPI process
  INTEGER(idp), DIMENSION(:), POINTER :: cell_y_min
  !> Maximal cell number in y for each MPI process
  INTEGER(idp), DIMENSION(:), POINTER :: cell_y_max
  !> Minimum cell number in z for each MPI process
  INTEGER(idp), DIMENSION(:), POINTER :: cell_z_min, cell_z_min_r, cell_z_min_f
  !> Maximal cell number in z for each MPI process
  INTEGER(idp), DIMENSION(:), POINTER :: cell_z_max, cell_z_max_r, cell_z_max_f
  INTEGER(idp), DIMENSION(:), POINTER :: cell_y_min_r,cell_y_max_r,cell_y_min_f,cell_y_max_f
  !> Used in em3dsolverPXR.py
  INTEGER(idp), DIMENSION(:), POINTER :: new_cell_x_min
  !> Used in em3dsolverPXR.py
  INTEGER(idp), DIMENSION(:), POINTER :: new_cell_x_max
  !> Used in em3dsolverPXR.py
  INTEGER(idp), DIMENSION(:), POINTER :: new_cell_y_min
  !> Used in em3dsolverPXR.py
  INTEGER(idp), DIMENSION(:), POINTER :: new_cell_y_max
  !> Used in em3dsolverPXR.py
  INTEGER(idp), DIMENSION(:), POINTER :: new_cell_z_min
  !> Used in em3dsolverPXR.py
  INTEGER(idp), DIMENSION(:), POINTER :: new_cell_z_max
  !> Minimum node number in x in the current MPI process
  INTEGER(idp)                        :: nx_global_grid_min
  !> Maximal node number in x in the current MPI process
  INTEGER(idp)                        :: nx_global_grid_max
  !> Minimum node number in y in the current MPI process
  INTEGER(idp)                        :: ny_global_grid_min
  !> Maximal node number in y in the current MPI process
  INTEGER(idp)                        :: ny_global_grid_max
  !> Minimum node number in z in the current MPI process
  INTEGER(idp)                        :: nz_global_grid_min
  !> Maximal node number in z in the current MPI process
  INTEGER(idp)                        :: nz_global_grid_max
  ! Domain axis
  !> Flag true when arrays of axis (x, y, z, x_global, y_global, z_global)
  !> are allocated (see mpi_routines.F90)
  LOGICAL(lp)                         :: l_axis_allocated=.FALSE.
  !> Global x axis cell array with guard cells (-nxguards:nx_global+nxguards)
  REAL(num), DIMENSION(:), POINTER    :: x_global
  !> Global y axis cell array with guard cells (-nyguards:ny_global+nyguards)
  REAL(num), DIMENSION(:), POINTER    :: y_global
  !> Global z axis cell array with guard cells (-nzguards:nz_global+nzguards)
  REAL(num), DIMENSION(:), POINTER    :: z_global
  ! domain limits and size
  !> local number of cells in x
  INTEGER(idp)                        :: nx
  !> local number of cells in y
  INTEGER(idp)                        :: ny
  !> local number of cells in z
  INTEGER(idp)                        :: nz
  !> local number of cells in z for  mpi groups load balancing 
  INTEGER(idp)                        :: nz_lb
  !> local number of grid points in z for  mpi groups load balancing 
  INTEGER(idp)                        :: nz_grid_lb
  !> local number of grid points in x
  INTEGER(idp)                        :: nx_grid
  !> local number of grid points in y
  INTEGER(idp)                        :: ny_grid
  !> local number of grid points in z
  INTEGER(idp)                        :: nz_grid
  !> local number of grid points in kx (Fourier Space)
  INTEGER(idp)                        :: nkx
  !> local number of grid points in ky (Fourier Space)
  INTEGER(idp)                        :: nky
  !> local number of grid points in kz (Fourier Space)"
  INTEGER(idp)                        :: nkz
  !> global number of cells in x
  INTEGER(idp)                        :: nx_global
  !> global number of cells in y
  INTEGER(idp)                        :: ny_global
  !> global number of cells in z
  INTEGER(idp)                        :: nz_global
  !> global number of grid points in x
  INTEGER(idp)                        :: nx_global_grid
  !> global number of grid points in y
  INTEGER(idp)                        :: ny_global_grid
  !> global number of grid points in z
  INTEGER(idp)                        :: nz_global_grid
  !> Space step in x
  REAL(num)                           :: dx
  !> Global grid minimal limit in x
  REAL(num)                           :: xmin
  !> Global grid maximal limit in x
  REAL(num)                           :: xmax
  !> global length in x: xmax - xmin
  REAL(num)                           :: length_x
  !> Local minimal grid limit in x
  REAL(num)                           :: x_min_local
  !> Local maximal grid limit in x
  REAL(num)                           :: x_max_local
  !> Space step in y
  REAL(num)                           :: dy
  !> Global grid minimal limit in y
  REAL(num)                           :: ymin
  !> Global grid maximal limit in y
  REAL(num)                           :: ymax
  !> global length in y: ymax - ymin
  REAL(num)                           :: length_y
  !> Local minimal grid limit in y
  REAL(num)                           :: y_min_local
  !> Local maximal grid limit in y
  REAL(num)                           :: y_max_local
  !> Space step in z
  REAL(num)                           :: dz
  !> Global grid minimal limit in z
  REAL(num)                           :: zmin
  !> Global grid maximal limit in z
  REAL(num)                           :: zmax
  !> global length in z: zmax - zmin
  REAL(num)                           :: length_z
  !> Local minimal grid limit in z
  REAL(num)                           :: z_min_local
  !> Local maximal grid limit in z
  REAL(num)                           :: z_max_local
  !> Local minimal particle domain limit in x
  REAL(num)                           :: x_min_local_part
  !> Local maximal particle domain limit in x
  REAL(num)                           :: x_max_local_part
  !> Local minimal particle domain limit in x
  REAL(num)                           :: y_min_local_part
  !> Local maximal particle domain limit in x
  REAL(num)                           :: y_max_local_part
  !> Local minimal particle domain limit in x
  REAL(num)                           :: z_min_local_part
  !> Local maximal particle domain limit in x
  REAL(num)                           :: z_max_local_part
  !> Global minimal particle domain limit in x
  REAL(num)                           :: xmin_part
  !> Global maximal particle domain limit in x
  REAL(num)                           :: xmax_part
  !> Global minimal particle domain limit in x
  REAL(num)                           :: ymin_part
  !> Global maximal particle domain limit in x
  REAL(num)                           :: ymax_part
  !> Global minimal particle domain limit in x
  REAL(num)                           :: zmin_part
  !> Global  maximal particle domain limit in x
  REAL(num)                           :: zmax_part
  !> global particle domain length in x: xmax_part - xmin_part
  REAL(num)                           :: length_x_part
  !> global particle domain length in y: ymax_part - ymin_part
  REAL(num)                           :: length_y_part
  !> global particle domain length in z: zmax_part - zmin_part
  REAL(num)                           :: length_z_part
  !> Offset between grid and particle limits (x min bound)
  !> Default is 0. NB:at present this can be used only for absorbing/reinjecting
  !> particle boundary conditions. Forced to 0 otherwise
  REAL(num)                           :: offset_grid_part_x_min =0.0_num
  !> Default is 0. NB:at present this can be used only for absorbing/reinjecting
  !> particle boundary conditions. Forced to 0 otherwise
  REAL(num)                           :: offset_grid_part_x_max =0.0_num
  !> Default is 0. NB:at present this can be used only for absorbing/reinjecting
  !> particle boundary conditions. Forced to 0 otherwise
  REAL(num)                           :: offset_grid_part_y_min =0.0_num
  !> Default is 0. NB:at present this can be used only for absorbing/reinjecting
  !> particle boundary conditions. Forced to 0 otherwise
  REAL(num)                           :: offset_grid_part_y_max =0.0_num
  !> Default is 0. NB:at present this can be used only for absorbing/reinjecting
  !> particle boundary conditions. Forced to 0 otherwise
  REAL(num)                           :: offset_grid_part_z_min =0.0_num
  !> Default is 0. NB:at present this can be used only for absorbing/reinjecting
  !> particle boundary conditions. Forced to 0 otherwise
  REAL(num)                           :: offset_grid_part_z_max =0.0_num
  ! Sorting
  !> Activation of the sorting
  INTEGER(idp) :: sorting_activated
  !> Bin space step in x for the sorting
  REAL(NUM)    :: sorting_dx
  !> Bin space step in y for the sorting
  REAL(NUM)    :: sorting_dy
  !> Bin space step in z for the sorting
  REAL(NUM)    :: sorting_dz
  !> Shift of the sorting grid in respect of the origin in x
  REAL(NUM)    :: sorting_shiftx
  !> Shift of the sorting grid in respect of the origin in x
  REAL(NUM)    :: sorting_shifty
  !> Shift of the sorting grid in respect of the origin in x
  REAL(NUM)    :: sorting_shiftz
  !> verbose for the sorting (depreciated)
  LOGICAL(lp)  :: sorting_verbose
  ! Axis
  !> Space dimension
  INTEGER(idp) :: c_dim = 3
  !> Local x axis cell array with guard cells (-nxguards:nx+nxguards)
  REAL(num), POINTER, DIMENSION(:) :: x
  !> Local y axis cell array with guard cells (-nyguards:ny+nyguards)
  REAL(num), POINTER, DIMENSION(:) :: y
  !> Local z axis cell array with guard cells (-nzguards:nz+nzguards)
  REAL(num), POINTER, DIMENSION(:) :: z
  !> Local grid minimum in x for each processor (1:nprocx)
  REAL(num), DIMENSION(:), POINTER :: x_grid_mins
  !> Local grid maximum in x for each processor (1:nprocx)
  REAL(num), DIMENSION(:), POINTER :: x_grid_maxs
  !> Local grid minimum in y for each processor (1:nprocy)
  REAL(num), DIMENSION(:), POINTER :: y_grid_mins
  !> Local grid maximum in y for each processor (1:nprocy)
  REAL(num), DIMENSION(:), POINTER :: y_grid_maxs
  !> Local grid minimum in z for each processor (1:nprocz)
  REAL(num), DIMENSION(:), POINTER :: z_grid_mins
  !> Local grid maximum in z for each processor (1:nprocz)
  REAL(num), DIMENSION(:), POINTER :: z_grid_maxs
  !> Minimal global grid limit in x
  REAL(num) :: x_grid_min
  !> Maximal global grid limit in x
  REAL(num) :: x_grid_max
  !> Minimal local grid limit in x
  REAL(num) :: x_grid_min_local
  !> Maximal local grid limit in x
  REAL(num) :: x_grid_max_local
  !> Minimal global grid limit in y
  REAL(num) :: y_grid_min
  !> Maximal global grid limit in y
  REAL(num) :: y_grid_max
  !> Minimal local grid limit in y
  REAL(num) :: y_grid_min_local
  !> Maximal local grid limit in y
  REAL(num) :: y_grid_max_local
  !> Minimal global grid limit in z
  REAL(num) :: z_grid_min
  !> Maximal global grid limit in z
  REAL(num) :: z_grid_max
  !> Minimal local grid limit in z
  REAL(num) :: z_grid_min_local
  !> Maximal local grid limit in z
  REAL(num) :: z_grid_max_local
  !> Total charge density
  REAL(num), POINTER, DIMENSION(:, :, :) :: rho, rhoold
  !> Electric Field divergence
  REAL(num), POINTER, DIMENSION(:, :, :) :: dive
  !> Current divergence
  REAL(num), POINTER, DIMENSION(:, :, :) :: divj
  !> Magnetic Field divergence
  REAL(num), POINTER, DIMENSIOn(:, :, :) :: divb

  ! Values used for load balancing
  REAL(num) :: mpitime_per_it
  REAL(num) :: max_time_per_it
  REAL(num) :: min_time_per_it
  REAL(num) :: global_time_per_cell
  REAL(num) :: global_time_per_part
  REAL(num) :: local_time_cell
  REAL(num) :: local_time_part
  INTEGER(idp) :: npart_local
  INTEGER(idp) :: npart_global
END MODULE shared_data

#if defined(FFTW)
MODULE fourier!#do not parse
  USE PICSAR_precision
  USE constants
  INTEGER(idp), DIMENSION(1) :: plan_r2c, plan_c2r
END MODULE fourier
#endif

! ________________________________________________________________________________________
!> Module for the Maxwell Solver coefficients
! ________________________________________________________________________________________
MODULE kyee_em3d
  USE PICSAR_precision
  USE constants
  !> alphax Maxwell coefficient = 7./12.
  REAL(num) :: alphax = 0.58333333333333337_num
  !> alphax Maxwell coefficient = 1./12.
  REAL(num) :: betaxy = 0.083333333333333329_num
  !> alphax Maxwell coefficient = 1./12.
  REAL(num) :: betaxz = 0.083333333333333329_num
  !> gammax Maxwell coefficient = 1./48.
  REAL(num) :: gammax = 0.020833333333333332_num
  !> alphay Maxwell coefficient = 7./12.
  REAL(num) :: alphay = 0.58333333333333337_num
  !> betayx Maxwell coefficient = 1./12.
  REAL(num) :: betayx = 0.083333333333333329_num
  !> betayz Maxwell coefficient = 1./12.
  REAL(num) :: betayz = 0.083333333333333329_num
  !> gammay Maxwell coefficient = 1./48.
  REAL(num) :: gammay = 0.020833333333333332_num
  !> alphaz Maxwell coefficient = 7./12.
  REAL(num) :: alphaz = 0.58333333333333337_num
  !> betazx Maxwell coefficient = 1./12.
  REAL(num) :: betazx = 0.083333333333333329_num
  !> betazy Maxwell coefficient = 1./12.
  REAL(num) :: betazy = 0.083333333333333329_num
  !> gammaz Maxwell coefficient = 1./48.
  REAL(num) :: gammaz = 0.020833333333333332_num
  !> Coefficient for the lehe solver
  REAL(num) :: deltaz = 0.000000000000000000_num
END MODULE kyee_em3d

! ________________________________________________________________________________________
!> Module containing pointer to the python arrays (used in em3dsolverPXR.py)
! ________________________________________________________________________________________
MODULE python_pointers
  USE PICSAR_precision
  USE constants
  !> Equivalent of pg.nps, the number of particles for each species
  INTEGER(idp), POINTER :: partn(:)
  !> Maximal number of particles, equivalent of pg.npmax in em3dsolverPXR
  INTEGER(idp) :: partnmax
  !> Number of guard cells in x, equivalent of particle_tile%nxg_tile
  INTEGER(idp) :: nxtg
  !> Number of guard cells in y, equivalent of particle_tile%nyg_tile
  INTEGER(idp) :: nytg
  !> Number of guard cells in z, equivalent of particle_tile%nzg_tile
  INTEGER(idp) :: nztg
  !> Number of nodes in x, equivalent of particle_tile%nx_grid_tile
  INTEGER(idp) :: nxgt
  !> Number of nodes in y, equivalent of particle_tile%ny_grid_tile
  INTEGER(idp) :: nygt
  !> Number of nodes in z, equivalent of particle_tile%nz_grid_tile
  INTEGER(idp) :: nzgt
  !> Number of cells in x, equivalent of particle_tile%nx_cells_tile
  INTEGER(idp) :: nxct
  !> Number of cells in y, equivalent of particle_tile%ny_cells_tile
  INTEGER(idp) :: nyct
  !> Number of cells in z, equivalent of particle_tile%nz_cells_tile
  INTEGER(idp) :: nzct
  !> Minimal cell index in x, equivalent of particle_tile%nx_tile_min
  INTEGER(idp) :: nxmin
  !> Maximal cell index in x, equivalent of particle_tile%nx_tile_max
  INTEGER(idp) :: nxmax
  !> Minimal cell index in y, equivalent of particle_tile%ny_tile_min
  INTEGER(idp) :: nymin
  !> Maximal cell index in y, equivalent of particle_tile%ny_tile_max
  INTEGER(idp) :: nymax
  !> Minimal cell index in z, equivalent of particle_tile%nz_tile_min
  INTEGER(idp) :: nzmin
  !> Maximal cell index in z, equivalent of particle_tile%nz_tile_max
  INTEGER(idp) :: nzmax
  ! Tile position
  !> Minimal tile limit in x: equivalent of particle_tile%x_tile_min
  REAL(num) :: xtmin
  !> Minimal tile limit in y: equivalent of particle_tile%y_tile_min
  REAL(num) :: ytmin
  !> Minimal tile limit in z: equivalent of particle_tile%z_tile_min
  REAL(num) :: ztmin
  !> Maximal tile limit in x: equivalent of particle_tile%x_tile_max
  REAL(num) :: xtmax
  !> Maximal tile limit in y: equivalent of particle_tile%y_tile_max
  REAL(num) :: ytmax
  !> Maximal tile limit in z: equivalent of particle_tile%z_tile_max
  REAL(num) :: ztmax
  !> Minimal grid tile boundary in x
  REAL(num) :: xgtmin
  !> Minimal grid tile boundary in y
  REAL(num) :: ygtmin
  !> Minimal grid tile boundary in z
  REAL(num) :: zgtmin
  !> Maximal grid tile boundary in x
  REAL(num) :: xgtmax
  !> Maximal grid tile boundary in y
  REAL(num) :: ygtmax
  !> Maximal grid tile boundary in z
  REAL(num) :: zgtmax
  !> array for particle x position
  REAL(num), DIMENSION(:), POINTER :: partx
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partx
#endif
  !DIR ATTRIBUTES FASTMEM  :: partx
  !> array for particle y position
  REAL(num), DIMENSION(:), POINTER :: party
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: party
#endif
  !DIR ATTRIBUTES FASTMEM  :: party
  !> array for particle z position
  REAL(num), DIMENSION(:), POINTER :: partz
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partz
#endif
  !DIR ATTRIBUTES FASTMEM  :: partz
  !> array for particle x momentum
  REAL(num), DIMENSION(:), POINTER :: partux
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partux
#endif
  !DIR ATTRIBUTES FASTMEM  :: partux
  !> array for particle y momentum
  REAL(num), DIMENSION(:), POINTER :: partuy
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partuy
#endif
  !DIR ATTRIBUTES FASTMEM  :: partuy
  !> array for particle z momentum
  REAL(num), DIMENSION(:), POINTER :: partuz
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partuz
#endif
  !DIR ATTRIBUTES FASTMEM  :: partuz
  !> array for the inverse of the particle gamma factor
  REAL(num), DIMENSION(:), POINTER :: partgaminv
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partgaminv
#endif
  !DIR ATTRIBUTES FASTMEM  :: partgaminv
  !> Array for particle weights and ids
  REAL(num), DIMENSION(:, :), POINTER :: pid
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: pid
#endif
  !DIR ATTRIBUTES FASTMEM  :: pid
  !> Particle Ex electric field
  REAL(num), DIMENSION(:), POINTER :: partex
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partex
#endif
  !DIR ATTRIBUTES FASTMEM  :: partex
  !> Particle Ey electric field
  REAL(num), DIMENSION(:), POINTER :: partey
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partey
#endif
  !DIR ATTRIBUTES FASTMEM  :: partey
  !> Particle Ez electric field
  REAL(num), DIMENSION(:), POINTER :: partez
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partez
#endif
  !DIR ATTRIBUTES FASTMEM  :: partez
  !> Particle Bx magnetic field
  REAL(num), DIMENSION(:), POINTER :: partbx
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partbx
#endif
  !DIR ATTRIBUTES FASTMEM  :: partbx
  !> Particle By magnetic field
  REAL(num), DIMENSION(:), POINTER :: partby
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partby
#endif
  !DIR ATTRIBUTES FASTMEM  :: partby
  !> Particle Bz magnetic field
  REAL(num), DIMENSION(:), POINTER :: partbz
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !dir$ attributes align:64 :: partbz
#endif
  !DIR ATTRIBUTES FASTMEM  :: partbz
END MODULE python_pointers

! ________________________________________________________________________________________
!> @brief
!> Module that stores memory sizes of Fortran arrays/derived types allocated on local rank 
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________
MODULE mem_status
  USE PICSAR_precision
  USE constants
  ! Memory size (in Bytes) occupied by grid arrays on local rank 
  REAL(num) :: local_grid_mem = 0._num
  ! Memory size (in Bytes) occupied by tiles grid arrays on local rank 
  ! (arrays used for field gathering on each tile)
  REAL(num) :: local_grid_tiles_mem = 0._num
  ! Memory size (in Bytes) occupied by particle arrays (tile structure) on local rank 
  REAL(num) :: local_part_tiles_mem = 0._num
  ! Total (all ranks) of memory size (in Bytes) occupied by grid arrays 
  REAL(num) :: global_grid_mem = 0._num
  ! Total (all ranks) of memory size (in Bytes) occupied by tile grid arrays 
  REAL(num) :: global_grid_tiles_mem = 0._num
  ! Total (all ranks) of memory size (in Bytes) occupied by particle
  ! arrays/structures on local rank 
  REAL(num) :: global_part_tiles_mem = 0._num
END MODULE mem_status

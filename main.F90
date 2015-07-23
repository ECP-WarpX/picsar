!===============================================================================
! WARPCORE 3D version 1.0, H. VINCENTI&J-L Vay 05/07/2015
! INCLUDES: 
! - Arbitrary order field solver (Maxwell.F90)
! - High order current deposition/field gathering (current_deposition.F90, field_gathering.F90)
! - MPI-parallelization (mpi_subtype_control.F90, mpi_routines.F90, boundary.F90)
! - OpenMP Paralelization (current_deposition.F90, field_gathering.F90, particle_push.F90, Maxwell.F90) 
!===============================================================================

PROGRAM main
!===============================================================================
USE constants
USE fields
USE particles
USE params
USE shared_data
USE mpi_routines

IMPLICIT NONE
INTEGER :: i,ierror,j,l

!----------------------------------------------
! INPUT PARAMETERS - EQUIVALENT TO INPUT FILE
!----------------------------------------------
! --- sets order of current deposition (between 1 and 3)
nox = 3
noy = 3
noz = 3

! --- sets order of Maxwell's solver in the three dimensions
norderx = 2
nordery = 2
norderz = 2

! --- number of particles per cell 
nppcell=1

! --- smoothing
npass = 0
alpha = 0.5_num

! --- sets coefficient multiplying Courant time step
dtcoef = 0.7_num
l_nodalgrid = .FALSE.


! --- sets max time in the simulation (in 1/w0)
tmax = 40.0_num

!-------------------------------------------------------------------------------
! plasma parameters (cold plasma)
l_particles_weight = .FALSE. ! particles have the same weight
vthx   = 0.0_num*clight      ! initial velocity spread on x(electrons only)
vthy   = 0.0_num*clight	! initial velocity spread on y (electrons only)
vthz   = 0.0_num*clight      ! initial velocity spread on z (electrons only)
theta = 0.0_num*pi           ! initial angle

! --- quantities in plasma (or lab) frame
!-------------------------------------------------------------------------------
g0    = 13.0_num          ! initial gamma
b0    = sqrt(1.0_num-1.0_num/g0**2)
nlab  = 1.e25_num            ! density in lab frame
nc    = nlab*g0          ! density (in the simulation frame)
wlab  = echarge*sqrt(nlab/(emass*eps0)) ! plasma frequency (in the lab frame)
w0_l  = echarge*sqrt(nc/(g0*emass*eps0))    ! "longitudinal" plasma frequency (in the lab frame)
w0_t  = echarge*sqrt(nc/(g0**3*emass*eps0)) ! "transverse" plasma frequency (in the lab frame)
w0    = w0_l

! --- sets # of cells
nx_global = 90
ny_global = 90
nz_global = 90

! --- sets resolution and grid size (in SI units)
dx = 6.5e-7_num
dy = 6.5e-7_num
dz = 6.5e-7_num

! --- sets simulation domain boundaries
xmin = 0.0_num
ymin = 0.0_num
zmin = 0.0_num
xmax = nx_global*dx
ymax = ny_global*dy
zmax = nz_global*dz

! --- sets domain decomposition 
nprocx = 1
nprocy = 1
nprocz = 1


!----------------------------------------------
! SET-UP COMMUNICATOR, ALLOCATE/INIT ARRAYS
!----------------------------------------------

! --- sets mpi communicator
  CALL mpi_minimal_init

! --- sets domain decompositions and allocate field arrays
  CALL mpi_initialise

! --- allocates and inits particle distributions (on each subdomain)
  CALL initall

!----------------------------------------------
! THIS IS THE PIC ALGORITHM TIME LOOP
!----------------------------------------------
IF (rank .EQ. 0) startsim=MPI_WTIME()
CALL step(nsteps)
IF (rank .EQ. 0) endsim=MPI_WTIME()
CALL mpi_close

END PROGRAM main

!===============================================================================
! Contains shared data
!===============================================================================
MODULE constants
!===============================================================================
INTEGER, PARAMETER :: num = KIND(1.d0)
REAL(num), parameter :: emass   = 9.1093818800000006e-31,      &
                        pmass   = 1.6726231000000001e-27,      &
                        echarge = 1.6021764620000001e-19,      &
                        clight  = 299792458.0,                 &
                        mu0     = 1.2566370614359173e-06,      &
                        eps0    = 8.8541878176203892e-12,      &
                        pi      = 3.141592653589793
INTEGER, PARAMETER :: c_ndims = 3
! direction parameters
INTEGER, PARAMETER :: c_dir_x = 1
INTEGER, PARAMETER :: c_dir_y = 2
INTEGER, PARAMETER :: c_dir_z = 3
LOGICAL:: l_smooth_compensate

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
INTEGER :: nparte, npartp, nppcell
REAL(num) :: fdxrand=0.0_num,fdzrand=0.0_num,vthx=0.0_num,vthy=0.0_num,vthz=0.0_num
REAL(num), POINTER, DIMENSION(:) :: xe, ye, ze, uxe, uye, uze, &
                                       exe, eye, eze, bxe, bye, bze
REAL(num), POINTER, DIMENSION(:) :: xp, yp, zp, uxp, uyp, uzp, &
                                       exp, eyp, ezp, bxp, byp, bzp
REAL(num), POINTER, DIMENSION(:) :: we, wp
END MODULE particles

!===============================================================================
MODULE params
!===============================================================================
USE constants
INTEGER :: it,nsteps
REAL(num) :: g0,b0,dt,w0,dtcoef,tmax
REAL(num) :: theta,nlab,wlab,nc,w0_l,w0_t,r,th
LOGICAL :: l_arrays_allocated= .FALSE., l_ck=.FALSE.
END MODULE params

!===============================================================================
MODULE shared_data
!===============================================================================
USE mpi
USE constants
!----------------------------------------------------------------------------
! MPI data
!----------------------------------------------------------------------------
INTEGER  :: mpireal = MPI_DOUBLE_PRECISION
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
INTEGER :: nx_global_min, nx_global_max
INTEGER :: ny_global_min, ny_global_max
INTEGER :: nz_global_min, nz_global_max
INTEGER :: n_global_min(c_ndims), n_global_max(c_ndims)
! domain and loadbalancing
LOGICAL :: allow_cpu_reduce = .FALSE.
REAL(num), DIMENSION(:), ALLOCATABLE :: x_global, y_global, z_global
REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global
REAL(num), DIMENSION(:), ALLOCATABLE :: xb_offset_global
REAL(num), DIMENSION(:), ALLOCATABLE :: yb_offset_global
REAL(num), DIMENSION(:), ALLOCATABLE :: zb_offset_global
! domain limits and size
INTEGER  :: nx, ny, nz ! local sizes
INTEGER  :: nx_global, ny_global, nz_global ! global sizes
REAL(num):: dx, xmin, xmax, length_x
REAL(num):: x_min_local, x_max_local
REAL(num):: dy, ymin, ymax,length_y
REAL(num):: y_min_local, y_max_local
REAL(num):: dz, zmin, zmax,length_z
REAL(num):: z_min_local, z_max_local
! Subtypes (field exchange)
INTEGER :: subtype_field, subtype_field_r4
INTEGER :: subarray_field, subarray_field_r4
INTEGER :: subarray_field_big, subarray_field_big_r4

! Axis
real(kind=8), pointer, dimension(:) :: x, y, z
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
REAL(num), POINTER, DIMENSION(:,:,:) :: rho

! IO/ STATISTICS
REAL(num) :: starttime=0., startsim=0., endsim=0.
REAL(num) :: runtime=0.
REAL(num), ALLOCATABLE, DIMENSION(:) :: pushb, bcs_pushb, pushe, bcs_pushe, &
                                        push_part, bcs_part, cs, bcs_cs, field_gath
INTEGER, PARAMETER :: string_length = 264
CHARACTER(LEN=string_length) :: fileex='ex'
CHARACTER(LEN=string_length) :: fileey='ey'
CHARACTER(LEN=string_length) :: fileez='ez'
CHARACTER(LEN=string_length) :: filebx='bx'
CHARACTER(LEN=string_length) :: fileby='by'
CHARACTER(LEN=string_length) :: filebz='bz'
CHARACTER(LEN=string_length) :: filejx='jx'
CHARACTER(LEN=string_length) :: filedive='dive'
CHARACTER(LEN=string_length) :: filerho='rho'
END MODULE shared_data



!------------------
! THIS IS A MODULE 
!------------------
module constants

integer, parameter :: num = 8
integer, parameter :: isp = 4
integer, parameter :: idp = 8
real(kind=8), parameter :: emass   = 9.1093818800000006e-31,pmass   = 1.6726231000000001e-27_num,echarge = 1.6021764620000001e-19_num,clight  = 299792458.0_num,mu0     = 1.2566370614359173e-06_num,eps0    = 8.8541878176203892e-12_num,pi      = 3.141592653589793_num
integer(isp), parameter :: c_ndims = 3

integer, parameter :: c_dir_x = 1
integer, parameter :: c_dir_y = 2
integer, parameter :: c_dir_z = 3
logical:: l_smooth_compensate
integer, parameter :: string_length = 264

integer, parameter :: c_err_bad_value = 2**4
end module constants
!------------------
! THIS IS A MODULE 
!------------------
module fields

use constants
logical:: l_lower_order_in_v, l_nodalgrid
integer(idp):: nxs=0, nys=0, nzs=0
integer(idp):: norderx, nordery, norderz
integer(idp):: nxguards,nyguards, nzguards, nox, noy, noz, npass(3)
integer(idp):: nxjguards,nyjguards, nzjguards
real(kind=8):: alpha(3)
real(kind=8), pointer, dimension(:,:,:) :: ex,ey,ez,bx,by,bz,jx,jy,jz

real(kind=8), pointer, dimension(:) :: xcoeffs, ycoeffs, zcoeffs
end module fields
!------------------
! THIS IS A MODULE 
!------------------
module particle_tilemodule !#do not parse
use constants
type particle_tile
logical :: l_arrays_allocated= .false.

integer(idp) :: np_tile, npmax_tile
integer(idp) :: nx_grid_tile, ny_grid_tile, nz_grid_tile
integer(idp) :: nx_cells_tile, ny_cells_tile, nz_cells_tile
integer(idp) :: nx_tile_min, nx_tile_max, ny_tile_min, ny_tile_max,nz_tile_min, nz_tile_max

real(kind=8) :: x_tile_min, y_tile_min, z_tile_min
real(kind=8) :: x_tile_max, y_tile_max, z_tile_max
real(kind=8) :: x_grid_tile_min, y_grid_tile_min, z_grid_tile_min
real(kind=8) :: x_grid_tile_max, y_grid_tile_max, z_grid_tile_max

logical :: subdomain_bound = .false.

real(kind=8), pointer, dimension(:) :: part_x
real(kind=8), pointer, dimension(:) :: part_y
real(kind=8), pointer, dimension(:) :: part_z
real(kind=8), pointer, dimension(:) :: part_ux
real(kind=8), pointer, dimension(:) :: part_uy
real(kind=8), pointer, dimension(:) :: part_uz
real(kind=8), pointer, dimension(:) :: part_ex
real(kind=8), pointer, dimension(:) :: part_ey
real(kind=8), pointer, dimension(:) :: part_ez
real(kind=8), pointer, dimension(:) :: part_bx
real(kind=8), pointer, dimension(:) :: part_by
real(kind=8), pointer, dimension(:) :: part_bz
real(kind=8), pointer, dimension(:) :: weight

real(kind=8), pointer, dimension(:,:,:) :: jxtile
real(kind=8), pointer, dimension(:,:,:) :: jytile
real(kind=8), pointer, dimension(:,:,:) :: jztile
real(kind=8), dimension(:,:,:), pointer :: rhotile
end type
end module particle_tilemodule !#do not parse
!------------------
! THIS IS A MODULE 
!------------------
module particle_speciesmodule !#do not parse #do not parse
use particle_tilemodule
use constants

type particle_species

character(len=string_length) :: name
real(kind=8) :: charge
real(kind=8) :: mass
real(kind=8) :: x_min
real(kind=8) :: x_max
real(kind=8) :: y_min
real(kind=8) :: y_max
real(kind=8) :: z_min
real(kind=8) :: z_max
real(kind=8) :: vdrift_x
real(kind=8) :: vdrift_y
real(kind=8) :: vdrift_z
real(kind=8) :: vth_x
real(kind=8) :: vth_y
real(kind=8) :: vth_z
integer(idp)   :: species_npart
integer(idp)   :: nppspecies_max
integer(idp)   :: nppcell
logical(idp)   :: l_arrayoftiles_allocated =.false.


type(particle_tile), dimension(:,:,:), pointer :: array_of_tiles
end type
end module particle_speciesmodule !#do not parse #do not parse
!------------------
! THIS IS A MODULE 
!------------------
module tile_params

integer :: ntilex, ntiley, ntilez

end module tile_params
!------------------
! THIS IS A MODULE 
!------------------
module particle_properties
use constants
integer(idp), parameter  :: nthreads_tile=1
logical :: l_initongrid = .false.
logical :: l_particles_weight = .false.
logical :: l4symtry = .false.
integer(idp) :: pdistr
integer(idp) :: nspecies
integer(idp) :: ntot
integer(idp), parameter :: nspecies_max=4
real(kind=8) :: fdxrand=0.0,fdzrand=0.0,vthx=0.0,vthy=0.0,vthz=0.0
logical :: l_species_allocated=.false.
end module particle_properties
!------------------
! THIS IS A MODULE 
!------------------
module particles #do not parse

use constants
use tile_params
use particle_tilemodule
use particle_speciesmodule
use particle_properties


type(particle_species), pointer, dimension(:):: species_parray
end module particles #do not parse
!------------------
! THIS IS A MODULE 
!------------------
module params

use constants
integer(idp) :: it,nsteps
real(kind=8) :: g0,b0,dt,w0,dtcoef,tmax
real(kind=8) :: theta,nlab,wlab,nc,w0_l,w0_t
logical :: l_coeffs_allocated= .false., l_ck=.false.
real(kind=8), parameter :: resize_factor=1.5
end module params
!------------------
! THIS IS A MODULE 
!------------------
module mpi_type_constants !#do not parse

use mpi
use constants
integer(isp)  :: mpidbl = mpi_double_precision
integer(isp) :: status(mpi_status_size)

integer(isp) :: derived_type_grid
integer(isp) :: derived_subarray_grid
end module mpi_type_constants !#do not parse
!------------------
! THIS IS A MODULE 
!------------------
module output_data !#do not parse

use constants


real(kind=8) :: startsim =0.0
real(kind=8) :: endsim =0.0
real(kind=8) :: startit, timeit
real(kind=8) :: pushtime




integer(idp) :: output_frequency = -1
integer(idp) :: output_step_min = 0
integer(idp) :: output_step_max = 0


integer(kind=4) :: c_output_ex = 0
integer(kind=4) :: c_output_ey = 0
integer(kind=4) :: c_output_ez = 0
integer(kind=4) :: c_output_bx = 0
integer(kind=4) :: c_output_by = 0
integer(kind=4) :: c_output_bz = 0
integer(kind=4) :: c_output_jx = 0
integer(kind=4) :: c_output_jy = 0
integer(kind=4) :: c_output_jz = 0
integer(kind=4) :: c_output_rho = 0
integer(kind=4) :: c_output_dive = 0


character(len=string_length) :: fileex   ='ex'
character(len=string_length) :: fileey   ='ey'
character(len=string_length) :: fileez   ='ez'
character(len=string_length) :: filebx   ='bx'
character(len=string_length) :: fileby   ='by'
character(len=string_length) :: filebz   ='bz'
character(len=string_length) :: filejx   ='jx'
character(len=string_length) :: filejy   ='jy'
character(len=string_length) :: filejz   ='jz'
character(len=string_length) :: filedive ='dive'
character(len=string_length) :: filerho  ='rho'

end module output_data !#do not parse
!------------------
! THIS IS A MODULE 
!------------------
module shared_data

use mpi_type_constants
use output_data



integer(isp) :: errcode, provided, comm, tag, rank
integer(isp) :: coordinates(3), neighbour(-1:1, -1:1, -1:1)
integer(isp) :: x_coords, proc_x_min, proc_x_max
integer(isp):: y_coords, proc_y_min, proc_y_max
integer(isp) :: z_coords, proc_z_min, proc_z_max
integer(isp) :: nproc, nprocx, nprocy, nprocz
integer(isp) :: nprocdir(3)
integer(idp), pointer, dimension(:) :: nx_each_rank, ny_each_rank, nz_each_rank
integer(idp), pointer, dimension(:) :: npart_each_rank
logical :: x_min_boundary, x_max_boundary
logical :: y_min_boundary, y_max_boundary
logical :: z_min_boundary, z_max_boundary

integer(idp), dimension(:), pointer :: cell_x_min, cell_x_max
integer(idp), dimension(:), pointer :: cell_y_min, cell_y_max
integer(idp), dimension(:), pointer :: cell_z_min, cell_z_max
integer(idp), dimension(:), pointer :: old_x_max, old_y_max, old_z_max
integer(idp) :: nx_global_grid_min, nx_global_grid_max
integer(idp) :: ny_global_grid_min, ny_global_grid_max
integer(idp) :: nz_global_grid_min, nz_global_grid_max
integer(idp) :: n_global_grid_min(3), n_global_grid_max(3)

logical :: allow_cpu_reduce = .false.
real(kind=8), dimension(:), pointer :: x_global, y_global, z_global
real(kind=8), dimension(:), pointer :: xb_global, yb_global, zb_global
real(kind=8), dimension(:), pointer :: xb_offset_global
real(kind=8), dimension(:), pointer :: yb_offset_global
real(kind=8), dimension(:), pointer :: zb_offset_global

integer(idp)  :: nx, ny, nz
integer(idp)  :: nx_grid, ny_grid, nz_grid
integer(idp)  :: nx_global, ny_global, nz_global
integer(idp)  :: nx_global_grid, ny_global_grid, nz_global_grid
real(kind=8):: dx, xmin, xmax, length_x
real(kind=8):: x_min_local, x_max_local
real(kind=8):: dy, ymin, ymax,length_y
real(kind=8):: y_min_local, y_max_local
real(kind=8):: dz, zmin, zmax,length_z
real(kind=8):: z_min_local, z_max_local


real(kind=8), pointer, dimension(:) :: x, y, z
real(kind=8), dimension(:), pointer :: x_grid_mins, x_grid_maxs
real(kind=8), dimension(:), pointer :: y_grid_mins, y_grid_maxs
real(kind=8), dimension(:), pointer :: z_grid_mins, z_grid_maxs
real(kind=8) ::  x_grid_min, x_grid_max
real(kind=8) :: x_grid_min_local, x_grid_max_local
real(kind=8) ::  y_grid_min, y_grid_max
real(kind=8) :: y_grid_min_local, y_grid_max_local
real(kind=8) :: z_grid_min, z_grid_max
real(kind=8) :: z_grid_min_local, z_grid_max_local


real(kind=8), pointer, dimension(:,:,:) :: rho

real(kind=8), pointer, dimension(:,:,:) :: dive

end module shared_data
!------------------
! THIS IS A MODULE 
!------------------
module python_pointers
use constants
integer(idp) :: partn
real(kind=8), dimension(:), pointer :: partx
real(kind=8), dimension(:), pointer :: party
real(kind=8), dimension(:), pointer :: partz
real(kind=8), dimension(:), pointer :: partux
real(kind=8), dimension(:), pointer :: partuy
real(kind=8), dimension(:), pointer :: partuz
real(kind=8), dimension(:), pointer :: partw
end module python_pointers
!------------------
! THIS IS A MODULE 
!------------------
module tiling #do not parse

use constants
use particles
use shared_data
use fields
use params
implicit none

end module tiling
!------------------
! THIS IS A MODULE 
!------------------
module mpi_derived_types !#do not parse #do not parse






use shared_data

implicit none

contains






function create_current_grid_derived_type()

integer(isp) :: create_current_grid_derived_type

create_current_grid_derived_type =create_grid_derived_type(mpidbl, nx, ny, nz, nx_global_grid_min,ny_global_grid_min, nz_global_grid_min)

end function create_current_grid_derived_type






function create_current_grid_subarray(ngx,ngy,ngz)

integer(isp) :: create_current_grid_subarray
integer(idp), intent(in) :: ngx, ngy, ngz
integer(idp) :: nxloc, nyloc, nzloc
nxloc=nx+1
nyloc=ny+1
nzloc=nz+1


create_current_grid_subarray =create_grid_subarray(mpidbl, ngx, ngy, ngz,nxloc, nyloc, nzloc)

end function create_current_grid_subarray







function create_grid_derived_type(mpitype, nx_local, ny_local, nz_local,cell_start_x_local, cell_start_y_local, cell_start_z_local)

integer(isp), intent(in) :: mpitype
integer(idp), intent(in) :: nx_local
integer(idp), intent(in) :: ny_local
integer(idp), intent(in) :: nz_local
integer(idp), intent(in) :: cell_start_x_local
integer(idp), intent(in) :: cell_start_y_local
integer(idp), intent(in) :: cell_start_z_local
integer(isp) :: create_grid_derived_type
integer(idp), dimension(c_ndims) :: n_local, n_global, start

n_local = (/nx_local+1, ny_local+1, nz_local+1/)
n_global = (/nx_global+1, ny_global+1, nz_global+1/)
start = (/cell_start_x_local, cell_start_y_local, cell_start_z_local/)

create_grid_derived_type =create_3d_array_derived_type(mpitype, n_local, n_global, start)

end function create_grid_derived_type






function create_3d_array_derived_type(mpitype, n_local, n_global, start)result(vec3d_sub)

integer(isp), intent(in) :: mpitype
integer(idp), dimension(c_ndims), intent(in) :: n_local
integer(idp), dimension(c_ndims), intent(in) :: n_global
integer(idp), dimension(c_ndims), intent(in) :: start
integer(isp), dimension(c_ndims) :: lengths, types
integer(kind=mpi_address_kind) :: disp(3), starts(3)
integer(isp) :: vec2d, vec2d_sub
integer(isp) :: vec3d, vec3d_sub, typesize

vec2d = mpi_datatype_null
call mpi_type_vector(int(n_local(2),isp), int(n_local(1),isp), int(n_global(1),isp), mpitype,vec2d, errcode)
call mpi_type_commit(vec2d, errcode)

call mpi_type_size(mpitype, typesize, errcode)
starts = start - 1
lengths = 1

disp(1) = 0
disp(2) = typesize * (starts(1) + n_global(1) * starts(2))
disp(3) = typesize * n_global(1) * n_global(2)
types(1) = mpi_lb
types(2) = vec2d
types(3) = mpi_ub

vec2d_sub = mpi_datatype_null
call mpi_type_create_struct(c_ndims, lengths, disp, types, vec2d_sub, errcode)
call mpi_type_commit(vec2d_sub, errcode)

vec3d = mpi_datatype_null
call mpi_type_contiguous(int(n_local(3),isp), vec2d_sub, vec3d, errcode)
call mpi_type_commit(vec3d, errcode)

disp(1) = 0
disp(2) = typesize * n_global(1) * n_global(2) * starts(3)
disp(3) = typesize * n_global(1) * n_global(2) * n_global(3)
types(1) = mpi_lb
types(2) = vec3d
types(3) = mpi_ub

vec3d_sub = mpi_datatype_null
call mpi_type_create_struct(c_ndims, lengths, disp, types, vec3d_sub, errcode)
call mpi_type_commit(vec3d_sub, errcode)

call mpi_type_free(vec2d, errcode)
call mpi_type_free(vec2d_sub, errcode)
call mpi_type_free(vec3d, errcode)

end function create_3d_array_derived_type



function create_grid_subarray(mpitype, ng1, ng2, ng3, n1, n2, n3)

integer(isp), intent(in) :: mpitype
integer(idp),intent(in) :: ng1,ng2,ng3, n1, n2, n3
integer(idp), dimension(3) :: n_local, ng, n_global, start
integer(isp) :: i, ndim, create_grid_subarray

n_local(1) = n1
n_local(2) = n2
n_local(3) = n3
ng(1)= ng1
ng(2)= ng2
ng(3)= ng3
ndim = 3
do i = 1, ndim
start(i) = 1 + ng(i)
n_global(i) = n_local(i) + 2 * ng(i)
enddo

create_grid_subarray =create_3d_array_derived_type(mpitype, n_local, n_global, start)

end function create_grid_subarray

end module mpi_derived_types !#do not parse
!------------------
! THIS IS A MODULE 
!------------------
module boundary #do not parse

use shared_data
use fields
use particles
use tiling
use mpi_derived_types
use constants

implicit none

end module boundary
!------------------
! THIS IS A MODULE 
!------------------
module simple_io #do not parse

use mpi_derived_types
use fields
use shared_data
implicit none

end module simple_io
!------------------
! THIS IS A MODULE 
!------------------
module diagnostics #do not parse

use constants
implicit none


end module diagnostics
!------------------
! THIS IS A MODULE 
!------------------
module mpi_routines #do not parse
use shared_data
use fields
use mpi
implicit none



real(kind=8) :: start_time, end_time

end module mpi_routines
!------------------
! THIS IS A MODULE 
!------------------
module control_file #do not parse

use shared_data
use params
use fields
use particles
implicit none

integer(idp) :: ios=0
integer(idp), parameter :: fh_input = 15
character(len=string_length) :: buffer
character(len=string_length) :: section_name

end module control_file
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_add_particle_to_species #do not parse
interface intef0205361
subroutine add_particle_to_species(currsp, partx, party, partz,partux, partuy, partuz, partw)

use tiling
implicit none
real(kind=8) :: partx, party, partz, partux, partuy, partuz, partw
type(particle_species), pointer, intent(in out) :: currsp
type(particle_tile), pointer :: curr
integer :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile, nptile
integer :: ixtile, iytile, iztile
end subroutine add_particle_to_species

end interface intef0205361
end module interf_add_particle_to_species
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_add_particle_at_tile #do not parse
interface intef0205686
subroutine add_particle_at_tile(curr, partx, party, partz,partux, partuy, partuz, partw)

use tiling
implicit none
integer :: count, nmax
real(kind=8) :: partx, party, partz, partux, partuy, partuz, partw
type(particle_tile), pointer, intent(in out) :: curr
end subroutine add_particle_at_tile

end interface intef0205686
end module interf_add_particle_at_tile
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_rm_particles_from_species #do not parse
interface intef0205867
subroutine rm_particles_from_species(currsp, curr, mask)

use tiling
type(particle_species), pointer, intent(in out) :: currsp
type(particle_tile), pointer, intent(in out) :: curr
logical, dimension (:), intent(in) :: mask
integer :: ninit, i
end subroutine rm_particles_from_species

end interface intef0205867
end module interf_rm_particles_from_species
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_rm_particle_at_tile #do not parse
interface intef0206051
subroutine rm_particle_at_tile(curr, index)

use tiling
implicit none
integer :: index
type(particle_tile), pointer, intent(in out) :: curr
end subroutine rm_particle_at_tile

end interface intef0206051
end module interf_rm_particle_at_tile
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_allocate_tile_arrays #do not parse
interface intef0206215
subroutine allocate_tile_arrays(curr_tile)

use tiling
type(particle_tile), pointer, intent(in out) :: curr_tile
integer :: nmax, nxc, nyc, nzc
end subroutine allocate_tile_arrays

end interface intef0206215
end module interf_allocate_tile_arrays
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_resize_particle_arrays #do not parse
interface intef0207542
subroutine resize_particle_arrays(curr, old_size, new_size)

use tiling
implicit none
type(particle_tile), pointer, intent(in out) :: curr
integer :: old_size, new_size
end subroutine resize_particle_arrays

end interface intef0207542
end module interf_resize_particle_arrays
!----------------------------------
!This is an interface module block 
!----------------------------------
module interf_resize_array_real #do not parse
interface intef0207716
subroutine resize_array_real(arr, old_size, new_size)

use tiling
implicit none
real(kind=8), dimension(:),pointer, intent(in out) :: arr
real(kind=8), dimension(:),pointer :: temp
integer :: old_size, new_size
end subroutine

end interface intef0207716
end module interf_resize_array_real

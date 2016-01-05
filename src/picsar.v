picsar


***** constants:
num integer /8/ 
isp integer /4/ 
idp integer /8/ 
emass real /9.1093818800000006e-31/ 
pmass real /1.6726231000000001e-27/ 
echarge real /1.6021764620000001e-19/ 
clight real /299792458.0/ 
mu0 real /1.2566370614359173e-06/ 
eps0 real /8.8541878176203892e-12/ 
pi real /3.141592653589793/ 
c_ndims integer /3/ 
c_dir_x integer /1/ 
c_dir_y integer /2/ 
c_dir_z integer /3/ 
l_smooth_compensate logical
string_length integer /264/ 
c_err_bad_value integer /2**4/ 

***** fields:
l_lower_order_in_v logical
l_nodalgrid logical
nxs integer /0/ 
nys integer /0/ 
nzs integer /0/ 
norderx integer
nordery integer
norderz integer
nxguards integer
nyguards integer
nzguards integer
nox integer
noy integer
noz integer
npass(3) integer
nxjguards integer
nyjguards integer
nzjguards integer
alpha(3) real
ex(:,:,:) _real
ey(:,:,:) _real
ez(:,:,:) _real
bx(:,:,:) _real
by(:,:,:) _real
bz(:,:,:) _real
jx(:,:,:) _real
jy(:,:,:) _real
jz(:,:,:) _real
xcoeffs(:) _real
ycoeffs(:) _real
zcoeffs(:) _real



***** tile_params:
ntilex integer
ntiley integer
ntilez integer

***** particle_properties:
nthreads_tile integer /1/ 
l_initongrid logical /.false./ 
l_particles_weight logical /.false./ 
l4symtry logical /.false./ 
pdistr integer
nspecies integer
ntot integer
nspecies_max integer /4/ 
fdxrand real /0.0/ 
fdzrand real /0.0/ 
vthx real /0.0/ 
vthy real /0.0/ 
vthz real /0.0/ 
l_species_allocated logical /.false./ 


***** params:
it integer
nsteps integer
g0 real
b0 real
dt real
w0 real
dtcoef real
tmax real
theta real
nlab real
wlab real
nc real
w0_l real
w0_t real
l_coeffs_allocated logical /.false./ 
l_ck logical /.false./ 
resize_factor real /1.5/ 



***** shared_data:
errcode integer
provided integer
comm integer
tag integer
rank integer
coordinates(3) integer
neighbour(-1:1,-1:1,-1:1) integer
x_coords integer
proc_x_min integer
proc_x_max integer
y_coords integer
proc_y_min integer
proc_y_max integer
z_coords integer
proc_z_min integer
proc_z_max integer
nproc integer
nprocx integer
nprocy integer
nprocz integer
nprocdir(3) integer
nx_each_rank(:) _integer
ny_each_rank(:) _integer
nz_each_rank(:) _integer
npart_each_rank(:) _integer
x_min_boundary logical
x_max_boundary logical
y_min_boundary logical
y_max_boundary logical
z_min_boundary logical
z_max_boundary logical
cell_x_min(:) _integer
cell_x_max(:) _integer
cell_y_min(:) _integer
cell_y_max(:) _integer
cell_z_min(:) _integer
cell_z_max(:) _integer
old_x_max(:) _integer
old_y_max(:) _integer
old_z_max(:) _integer
nx_global_grid_min integer
nx_global_grid_max integer
ny_global_grid_min integer
ny_global_grid_max integer
nz_global_grid_min integer
nz_global_grid_max integer
n_global_grid_min(3) integer
n_global_grid_max(3) integer
allow_cpu_reduce logical /.false./ 
x_global(:) _real
y_global(:) _real
z_global(:) _real
xb_global(:) _real
yb_global(:) _real
zb_global(:) _real
xb_offset_global(:) _real
yb_offset_global(:) _real
zb_offset_global(:) _real
nx integer
ny integer
nz integer
nx_grid integer
ny_grid integer
nz_grid integer
nx_global integer
ny_global integer
nz_global integer
nx_global_grid integer
ny_global_grid integer
nz_global_grid integer
dx real
xmin real
xmax real
length_x real
x_min_local real
x_max_local real
dy real
ymin real
ymax real
length_y real
y_min_local real
y_max_local real
dz real
zmin real
zmax real
length_z real
z_min_local real
z_max_local real
x(:) _real
y(:) _real
z(:) _real
x_grid_mins(:) _real
x_grid_maxs(:) _real
y_grid_mins(:) _real
y_grid_maxs(:) _real
z_grid_mins(:) _real
z_grid_maxs(:) _real
x_grid_min real
x_grid_max real
x_grid_min_local real
x_grid_max_local real
y_grid_min real
y_grid_max real
y_grid_min_local real
y_grid_max_local real
z_grid_min real
z_grid_max real
z_grid_min_local real
z_grid_max_local real
rho(:,:,:) _real
dive(:,:,:) _real

***** python_pointers:
partx(:) _real
party(:) _real
partz(:) _real
partux(:) _real
partuy(:) _real
partuz(:) _real
partw(:) _real















***** Subroutines:
push_bfield() subroutine
push_efield() subroutine
push_em3d_evec_norder(ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,by(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,mudt:real,dtsdx(norderx/2):real,dtsdy(nordery/2):real,dtsdz(norderz/2):real,nx:integer,ny:integer,nz:integer,norderx:integer,nordery:integer,norderz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nxs:integer,nys:integer,nzs:integer,l_nodalgrid:logical) subroutine
push_em3d_evec(ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,by(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,mudt:real,dtsdx:real,dtsdy:real,dtsdz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nxs:integer,nys:integer,nzs:integer,l_nodalgrid:logical) subroutine
push_em3d_bvec_norder(ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,by(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,dtsdx(norderx/2):real,dtsdy(nordery/2):real,dtsdz(norderz/2):real,nx:integer,ny:integer,nz:integer,norderx:integer,nordery:integer,norderz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nxs:integer,nys:integer,nzs:integer,l_nodalgrid:logical) subroutine
push_em3d_bvec(ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,by(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,dtsdx:real,dtsdy:real,dtsdz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nxs:integer,nys:integer,nzs:integer,l_nodalgrid:logical) subroutine
set_tile_split() subroutine
add_particle_to_species(currsp:_particle_species,partx:real,party:real,partz:real,partux:real,partuy:real,partuz:real,partw:real) subroutine
add_particle_at_tile(curr:_particle_tile,partx:real,party:real,partz:real,partux:real,partuy:real,partuz:real,partw:real) subroutine
rm_particles_from_species(currsp:_particle_species,curr:_particle_tile,mask:) subroutine
rm_particle_at_tile(curr:_particle_tile,index:integer) subroutine
allocate_tile_arrays(curr_tile:_particle_tile) subroutine
init_tile_arrays() subroutine
load_particles() subroutine
resize_particle_arrays(curr:_particle_tile,old_size:integer,new_size:integer) subroutine
resize_array_real(arr(:):_real,old_size:integer,new_size:integer) subroutine
point_to_tile(ispecies:integer,ix:integer,iy:integer,iz:integer) subroutine
set_particle_species_properties(nsp:integer,sname:string,mss:real,chrg:real,nppc:integer,xsmin:real,ysmin:real,zsmin:real,xsmax:real,ysmax:real,zsmax:real,vdxs:real,vdys:real,vdzs:real,vthxs:real,vthys:real,vthzs:real) subroutine
push_particles() subroutine
pushxyz(np:integer,xp(np):real,yp(np):real,zp(np):real,uxp(np):real,uyp(np):real,uzp(np):real,dt:real) subroutine
epush_v(np:integer,uxp(np):real,uyp(np):real,uzp(np):real,ex(np):real,ey(np):real,ez(np):real,q:real,m:real,dt:real) subroutine
bpush_v(np:integer,uxp(np):real,uyp(np):real,uzp(np):real,bx(np):real,by(np):real,bz(np):real,q:real,m:real,dt:real) subroutine
depose_currents_on_grid_jxjyjz() subroutine
depose_jxjyjz_scalar_1_1_1(jx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,np:integer,xp(np):real,yp(np):real,zp(np):real,uxp(np):real,uyp(np):real,uzp(np):real,w(np):real,q:real,xmin:real,ymin:real,zmin:real,dt:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer) subroutine
depose_jxjyjz_esirkepov_1_1_1(jx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,np:integer,xp(np):real,yp(np):real,zp(np):real,uxp(np):real,uyp(np):real,uzp(np):real,w(np):real,q:real,xmin:real,ymin:real,zmin:real,dt:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer) subroutine
depose_jxjyjz_esirkepov_n(jx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,jz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,np:integer,xp(np):real,yp(np):real,zp(np):real,uxp(np):real,uyp(np):real,uzp(np):real,w(np):real,q:real,xmin:real,ymin:real,zmin:real,dt:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nox:integer,noy:integer,noz:integer,l_particles_weight:logical,l4symtry:logical) subroutine
gete3d_energy_conserving_1_1_1(np:integer,xp(np):real,yp(np):real,zp(np):real,ex(np):real,ey(np):real,ez(np):real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,exg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,eyg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ezg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real) subroutine
getb3d_energy_conserving_1_1_1(np:integer,xp(np):real,yp(np):real,zp(np):real,bx(np):real,by(np):real,bz(np):real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,bxg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,byg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bzg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real) subroutine
gete3d_n_energy_conserving(np:integer,xp(np):real,yp(np):real,zp(np):real,ex(np):real,ey(np):real,ez(np):real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nox:integer,noy:integer,noz:integer,exg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,eyg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,ezg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,l_lower_order_in_v:logical) subroutine
getb3d_n_energy_conserving(np:integer,xp(np):real,yp(np):real,zp(np):real,bx(np):real,by(np):real,bz(np):real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nox:integer,noy:integer,noz:integer,bxg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,byg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,bzg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,l_lower_order_in_v:logical) subroutine
field_bc(field(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg):real,nxg:integer,nyg:integer,nzg:integer,nx_local:integer,ny_local:integer,nz_local:integer) subroutine
exchange_mpi_3d_grid_array_with_guards(field(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg):real,nxg:integer,nyg:integer,nzg:integer,nx_local:integer,ny_local:integer,nz_local:integer) subroutine
exchange_mpi_3d_grid_array_with_guards_nonblocking(field(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg):real,nxg:integer,nyg:integer,nzg:integer,nx_local:integer,ny_local:integer,nz_local:integer) subroutine
summation_bcs(array(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg):real,nxg:integer,nyg:integer,nzg:integer,nx_local:integer,ny_local:integer,nz_local:integer) subroutine
summation_bcs_nonblocking(array(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg):real,nxg:integer,nyg:integer,nzg:integer,nx_local:integer,ny_local:integer,nz_local:integer) subroutine
efield_bcs() subroutine
bfield_bcs() subroutine
current_bcs() subroutine
charge_bcs() subroutine
particle_bcs() subroutine
particle_bcs_tiles() subroutine
particle_bcs_mpi_blocking() subroutine
output_routines() subroutine
write_single_array_to_file(filename:string,array(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg):real,nxg:integer,nyg:integer,nzg:integer,nx_local:integer,ny_local:integer,nz_local:integer,offset:integer,err:integer) subroutine
calc_diags() subroutine
depose_rho_scalar_1_1_1(rho(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,np:integer,xp(np):real,yp(np):real,zp(np):real,w(np):real,q:real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer) subroutine
depose_rho_vechv_1_1_1(rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard)):real,np:integer,xp(np):real,yp(np):real,zp(np):real,w(np):real,q:real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer) subroutine
depose_rho_n(rho(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,np:integer,xp(np):real,yp(np):real,zp(np):real,w(np):real,q:real,xmin:real,ymin:real,zmin:real,dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,nox:integer,noy:integer,noz:integer,l_particles_weight:logical,l4symtry:logical) subroutine
calc_field_div(divee(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,eex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,eey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,eez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard):real,nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,dx:real,dy:real,dz:real) subroutine
step(nst:integer) subroutine
initall() subroutine
fd_weights(coeffs(norder/2):real,norder:integer,l_nodal:logical) subroutine
mpi_minimal_init() subroutine
setup_communicator() subroutine
mpi_initialise() subroutine
mpi_close() subroutine
default_init() subroutine
read_from_cl() subroutine
read_input_file() subroutine
read_cpusplit_section() subroutine
read_main_section() subroutine
read_species_section() subroutine
read_output_section() subroutine

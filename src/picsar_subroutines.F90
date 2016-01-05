!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_bfield
use constants
use params
use fields
use shared_data
implicit none

call push_em3d_bvec_norder(ex,ey,ez,bx,by,bz,0.5_num*dt/dx*xcoeffs,0.5_num*dt/dy*ycoeffs,0.5_num*dt/dz*zcoeffs,nx,ny,nz, norderx,nordery,norderz,nxguards,nyguards,nzguards,nxs,nys,nzs,l_nodalgrid)

end subroutine push_bfield
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_efield
use constants
use params
use fields
use shared_data
implicit none

call push_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,clight**2*mu0*dt,clight**2*dt/dx*xcoeffs,clight**2*dt/dy*ycoeffs,clight**2*dt/dz*zcoeffs,nx,ny,nz,norderx,nordery,norderz,nxguards,nyguards,nzguards,nxs,nys,nzs,l_nodalgrid)

end subroutine push_efield
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz,norderx,nordery,norderz,nxguard,nyguard,nzguard,nxs,nys,nzs,l_nodalgrid)

use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(kind=8), intent(in out), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(in), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: jx, jy, jz
real(kind=8), intent(in) :: mudt,dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
ist = 0
else
ist = 1
end if

!!$omp parallel default(shared) private(l,k,j,i)
!!$omp do collapse(3)

do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
ex(j,k,l) = ex(j,k,l) - mudt  * jx(j,k,l)
do i = 1, nordery/2
ex(j,k,l) = ex(j,k,l) + dtsdy(i) * (bz(j,k+i-ist,l)   - bz(j,k-i,l  ))
end do
do i = 1, norderz/2
ex(j,k,l) = ex(j,k,l) - dtsdz(i) * (by(j,k,l+i-ist)   - by(j,k  ,l-i))
end do
end do
end do
end do
!!$omp end do
!!$omp do collapse(3)

do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
ey(j,k,l) = ey(j,k,l) - mudt  * jy(j,k,l)
do i = 1, norderx/2
ey(j,k,l) = ey(j,k,l) - dtsdx(i) * (bz(j+i-ist,k,l)   - bz(j-i,k,l))
end do
do i = 1, norderz/2
ey(j,k,l) = ey(j,k,l) + dtsdz(i) * (bx(j,k,l+i-ist)   - bx(j,k,l-i))
end do
end do
end do
end do
!!$omp end do
!!$omp do collapse(3)

do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
ez(j,k,l) = ez(j,k,l) - mudt  * jz(j,k,l)
do i = 1, norderx/2
ez(j,k,l) = ez(j,k,l) + dtsdx(i) * (by(j+i-ist,k,l) - by(j-i,k  ,l))
end do
do i = 1, nordery/2
ez(j,k,l) = ez(j,k,l) - dtsdy(i) * (bx(j,k+i-ist,l) - bx(j  ,k-i,l))
end do
end do
end do
end do
!!$omp end do
!!$omp end parallel
return
end subroutine push_em3d_evec_norder
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_em3d_evec(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,l_nodalgrid)
use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs
real(kind=8), intent(in out), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(in), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: jx,jy,jz
real(kind=8), intent(in) :: mudt,dtsdx,dtsdy,dtsdz
integer(idp):: j,k,l
logical :: l_nodalgrid


do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
ex(j,k,l) = ex(j,k,l) + dtsdy * (bz(j,k,l)   - bz(j,k-1,l  ))- dtsdz * (by(j,k,l)   - by(j,k  ,l-1))- mudt  * jx(j,k,l)
end do
end do
end do


do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
ey(j,k,l) = ey(j,k,l) - dtsdx * (bz(j,k,l)   - bz(j-1,k,l))+ dtsdz * (bx(j,k,l)   - bx(j,k,l-1))- mudt  * jy(j,k,l)
end do
end do
end do


do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
ez(j,k,l) = ez(j,k,l) + dtsdx * (by(j,k,l) - by(j-1,k  ,l))- dtsdy * (bx(j,k,l) - bx(j  ,k-1,l))- mudt  * jz(j,k,l)
end do
end do
end do

return
end subroutine push_em3d_evec
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_em3d_bvec_norder(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,nx,ny,nz,norderx,nordery,norderz,nxguard,nyguard,nzguard,nxs,nys,nzs,l_nodalgrid)

use constants
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(kind=8), intent(in out), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(in) :: dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
integer(idp) :: i,j,k,l,ist
logical :: l_nodalgrid

if (l_nodalgrid) then
ist = 0
else
ist = 1
end if

!!$omp parallel default(shared) private(l,k,j,i)
!!$omp do collapse(3)

do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
do i = 1, nordery/2
bx(j,k,l) = bx(j,k,l) - dtsdy(i) * (ez(j,k+i,l  ) - ez(j,k-i+ist,l))
end do
do i = 1, norderz/2
bx(j,k,l) = bx(j,k,l) + dtsdz(i) * (ey(j,k,  l+i) - ey(j,k,l-i+ist))
end do
end do
end do
end do
!!$omp end do
!!$omp do collapse(3)

do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
do i = 1, norderx/2
by(j,k,l) = by(j,k,l) + dtsdx(i) * (ez(j+i,k,l  ) - ez(j-i+ist,k,l))
end do
do i = 1, norderz/2
by(j,k,l) = by(j,k,l) - dtsdz(i) * (ex(j  ,k,l+i) - ex(j,k,l-i+ist))
end do
end do
end do
end do
!!$omp end do
!!$omp do collapse(3)

do l = -nzs, nz+nzs
do k = -nys, ny+nys
do j = -nxs, nx+nxs
do i = 1, norderx/2
bz(j,k,l) = bz(j,k,l) - dtsdx(i) * (ey(j+i,k,l) - ey(j-i+ist,k,l))
end do
do i = 1, nordery/2
bz(j,k,l) = bz(j,k,l) + dtsdy(i) * (ex(j,k+i,l) - ex(j,k-i+ist,l))
end do
end do
end do
end do
!!$omp end do
!!$omp end parallel
return

end subroutine push_em3d_bvec_norder
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_em3d_bvec(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,l_nodalgrid)
use constants
integer(idp):: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs
real(kind=8), intent(in out), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(in) :: dtsdx,dtsdy,dtsdz
integer(idp) :: j,k,l
logical :: l_nodalgrid


do l = -nzs, nz+nzs-1
do k = -nys, ny+nys-1
do j = -nxs, nx+nxs
bx(j,k,l) = bx(j,k,l) - dtsdy * (ez(j,k+1,l  ) - ez(j,k,l))+ dtsdz * (ey(j,k,  l+1) - ey(j,k,l))
end do
end do
end do


do l = -nzs, nz+nzs-1
do k = -nys, ny+nys
do j = -nxs, nx+nxs-1
by(j,k,l) = by(j,k,l) + dtsdx * (ez(j+1,k,l  ) - ez(j,k,l))- dtsdz * (ex(j  ,k,l+1) - ex(j,k,l))
end do
end do
end do


do l = -nzs, nz+nzs
do k = -nys, ny+nys-1
do j = -nxs, nx+nxs-1
bz(j,k,l) = bz(j,k,l) - dtsdx * (ey(j+1,k,l) - ey(j,k,l))+ dtsdy * (ex(j,k+1,l) - ex(j,k,l))
end do
end do
end do
end subroutine push_em3d_bvec
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine set_tile_split
use tiling
implicit none
integer(kind=4) :: ix, iy, iz, ispecies
integer(kind=4) :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
integer(kind=4) :: nx0_last_tile, ny0_last_tile, nz0_last_tile
type(particle_species), pointer :: curr_sp
type(particle_tile), pointer :: curr


nx0_grid_tile = nx_grid / ntilex
ny0_grid_tile = ny_grid / ntiley
nz0_grid_tile = nz_grid / ntilez


if (nx0_grid_tile .lt. 4) then
if (rank .eq. 0) print *, "number of tiles in x to high, settting back to default value 1"
ntilex=1
end if
if (ny0_grid_tile .lt. 4) then
if (rank .eq. 0) print *, "number of tiles in y to high, setting back to default value 1"
ntiley=1
end if
if (nz0_grid_tile .lt. 4) then
if (rank .eq. 0) print *, "number of tiles in z to high, setting back to default value 1"
ntilez=1
end if


nx0_last_tile= nx0_grid_tile+(nx_grid-nx0_grid_tile*ntilex)
ny0_last_tile= ny0_grid_tile+(ny_grid-ny0_grid_tile*ntiley)
nz0_last_tile= nz0_grid_tile+(nz_grid-nz0_grid_tile*ntilez)


do ispecies =1, nspecies
curr_sp => species_parray(ispecies)
if (.not. curr_sp%l_arrayoftiles_allocated) then
allocate(curr_sp%array_of_tiles(ntilex,ntiley,ntilez))
curr_sp%l_arrayoftiles_allocated = .true.
end if

do iz=1, ntilez
do iy=1,ntiley
do ix=1,ntilex
curr=> curr_sp%array_of_tiles(ix,iy,iz)

if (ix .eq. 1) curr%subdomain_bound = .true.
if (ix .lt. ntilex) then
curr%nx_grid_tile=nx0_grid_tile
curr%nx_cells_tile=curr%nx_grid_tile-1
curr%x_tile_min=x_min_local+(ix-1)*nx0_grid_tile*dx
curr%x_grid_tile_min=curr%x_tile_min+dx/2.0_num
curr%x_tile_max=curr%x_tile_min+nx0_grid_tile*dx
curr%x_grid_tile_max=curr%x_tile_max-dx/2.0_num
curr%nx_tile_min = (ix-1)*nx0_grid_tile
curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
else
curr%subdomain_bound= .true.
curr%nx_grid_tile=nx0_last_tile
curr%nx_cells_tile=curr%nx_grid_tile-1
curr%x_tile_min=x_min_local+(ix-1)*nx0_grid_tile*dx
curr%x_grid_tile_min=curr%x_tile_min+dx/2.0_num
curr%x_tile_max=curr%x_tile_min+nx0_last_tile*dx
curr%x_grid_tile_max=curr%x_tile_max-dx/2.0_num
curr%nx_tile_min = (ix-1)*nx0_grid_tile
curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
endif

if (iy .eq. 1) curr%subdomain_bound = .true.
if (iy .lt. ntiley) then
curr%ny_grid_tile=ny0_grid_tile
curr%ny_cells_tile=curr%ny_grid_tile-1
curr%y_tile_min=y_min_local+(iy-1)*ny0_grid_tile*dy
curr%y_grid_tile_min=curr%y_tile_min+dy/2.0_num
curr%y_tile_max=curr%y_tile_min+ny0_grid_tile*dy
curr%y_grid_tile_max=curr%y_tile_max-dy/2.0_num
curr%ny_tile_min = (iy-1)*ny0_grid_tile
curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
else
curr%subdomain_bound= .true.
curr%ny_grid_tile=ny0_last_tile
curr%ny_cells_tile=curr%ny_grid_tile-1
curr%y_tile_min=y_min_local+(iy-1)*ny0_grid_tile*dy
curr%y_grid_tile_min=curr%y_tile_min+dy/2.0_num
curr%y_tile_max=curr%y_tile_min+ny0_last_tile*dy
curr%y_grid_tile_max=curr%y_tile_max-dy/2.0_num
curr%ny_tile_min = (iy-1)*ny0_grid_tile
curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
endif

if (iz .eq. 1) curr%subdomain_bound = .true.
if (iz .lt. ntilez) then
curr%nz_grid_tile=nz0_grid_tile
curr%nz_cells_tile=curr%nz_grid_tile-1
curr%z_tile_min=z_min_local+(iz-1)*nz0_grid_tile*dz
curr%z_grid_tile_min=curr%z_tile_min+dz/2.0_num
curr%z_tile_max=curr%z_tile_min+nz0_grid_tile*dz
curr%z_grid_tile_max=curr%z_tile_max-dz/2.0_num
curr%nz_tile_min = (iz-1)*nz0_grid_tile
curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
else
curr%subdomain_bound= .true.
curr%nz_grid_tile=nz0_last_tile
curr%nz_cells_tile=curr%nz_grid_tile-1
curr%z_tile_min=z_min_local+(iz-1)*nz0_grid_tile*dz
curr%z_grid_tile_min=curr%z_tile_min+dz/2.0_num
curr%z_tile_max=curr%z_tile_min+nz0_last_tile*dz
curr%z_grid_tile_max=curr%z_tile_max-dz/2.0_num
curr%nz_tile_min = (iz-1)*nz0_grid_tile
curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
endif
end do
end do
end do
end do
end subroutine
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine add_particle_to_species(currsp, partx, party, partz,partux, partuy, partuz, partw)
use interf_add_particle_at_tile
use tiling
implicit none
real(kind=8) :: partx, party, partz, partux, partuy, partuz, partw
type(particle_species), pointer, intent(in out) :: currsp
type(particle_tile), pointer :: curr
integer :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile, nptile
integer :: ixtile, iytile, iztile



nx0_grid_tile = currsp%array_of_tiles(1,1,1)%nx_grid_tile
ny0_grid_tile = currsp%array_of_tiles(1,1,1)%ny_grid_tile
nz0_grid_tile = currsp%array_of_tiles(1,1,1)%nz_grid_tile


ixtile = min(floor((partx-x_min_local)/(nx0_grid_tile*dx))+1,ntilex)
iytile = min(floor((party-y_min_local)/(ny0_grid_tile*dy))+1,ntiley)
iztile = min(floor((partz-z_min_local)/(nz0_grid_tile*dz))+1,ntilez)


curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
call add_particle_at_tile(curr, partx, party, partz,partux, partuy, partuz, partw)


currsp%species_npart=currsp%species_npart+1
end subroutine add_particle_to_species
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine add_particle_at_tile(curr, partx, party, partz,partux, partuy, partuz, partw)
use interf_resize_particle_arrays
use interf_allocate_tile_arrays
use tiling
implicit none
integer :: count, nmax
real(kind=8) :: partx, party, partz, partux, partuy, partuz, partw
type(particle_tile), pointer, intent(in out) :: curr

if (.not. curr%l_arrays_allocated) then
call allocate_tile_arrays(curr)
endif


count = curr%np_tile+1
nmax  = curr%npmax_tile
if (count .gt. nmax) then

call resize_particle_arrays(curr, nmax, nint(resize_factor*nmax+1))
endif

curr%np_tile=count
curr%part_x(count)  = partx
curr%part_y(count)  = party
curr%part_z(count)  = partz
curr%part_ux(count) = partux
curr%part_uy(count) = partuy
curr%part_uz(count) = partuz
curr%weight(count)  = partw
curr%part_ex(count)  = 0._num
curr%part_ey(count)  = 0._num
curr%part_ez(count)  = 0._num
curr%part_bx(count)  = 0._num
curr%part_by(count)  = 0._num
curr%part_bz(count)  = 0._num
end subroutine add_particle_at_tile
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine rm_particles_from_species(currsp, curr, mask)
use interf_rm_particle_at_tile
use tiling
type(particle_species), pointer, intent(in out) :: currsp
type(particle_tile), pointer, intent(in out) :: curr
logical, dimension (:), intent(in) :: mask
integer :: ninit, i
ninit= curr%np_tile
do i = ninit,1,-1
if (.not. mask(i)) then
call rm_particle_at_tile(curr,i)
currsp%species_npart=currsp%species_npart-1
endif
enddo
end subroutine rm_particles_from_species
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine rm_particle_at_tile(curr, index)
use tiling
implicit none
integer :: index
type(particle_tile), pointer, intent(in out) :: curr
if (index .eq. curr%np_tile) then


curr%np_tile=curr%np_tile-1
else


curr%part_x(index)=curr%part_x(curr%np_tile)
curr%part_y(index)=curr%part_y(curr%np_tile)
curr%part_z(index)=curr%part_z(curr%np_tile)
curr%part_ux(index)=curr%part_ux(curr%np_tile)
curr%part_uy(index)=curr%part_uy(curr%np_tile)
curr%part_uz(index)=curr%part_uz(curr%np_tile)
curr%weight(index)=curr%weight(curr%np_tile)
curr%np_tile=curr%np_tile-1
end if
end subroutine rm_particle_at_tile
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine allocate_tile_arrays(curr_tile)
use tiling
type(particle_tile), pointer, intent(in out) :: curr_tile
integer :: nmax, nxc, nyc, nzc

nmax = curr_tile%npmax_tile
allocate(curr_tile%part_x(1:nmax), curr_tile%part_y(1:nmax),curr_tile%part_z(1:nmax), curr_tile%part_ux(1:nmax),curr_tile%part_uy(1:nmax), curr_tile%part_uz(1:nmax),curr_tile%weight(1:nmax), curr_tile%part_ex(1:nmax),curr_tile%part_ey(1:nmax), curr_tile%part_ez(1:nmax),curr_tile%part_bx(1:nmax), curr_tile%part_by(1:nmax),curr_tile%part_bz(1:nmax))

nxc=curr_tile%nx_cells_tile
nyc=curr_tile%ny_cells_tile
nzc=curr_tile%nz_cells_tile
allocate(curr_tile%jxtile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),curr_tile%jytile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),curr_tile%jztile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),curr_tile%rhotile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards))
curr_tile%l_arrays_allocated = .true.

end subroutine allocate_tile_arrays
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine init_tile_arrays
use interf_allocate_tile_arrays
use tiling
implicit none
integer(kind=4) :: ispecies, ix, iy, iz
integer(kind=4) :: n1, n2, n3
type(particle_tile), pointer :: curr_tile
type(particle_species), pointer :: curr


do ispecies=1,nspecies
curr=>species_parray(ispecies)
curr%species_npart=0
do iz=1, ntilez
do iy=1, ntiley
do ix=1,ntilex
curr_tile=>curr%array_of_tiles(ix,iy,iz)

n1=curr_tile%nx_cells_tile
n2=curr_tile%ny_cells_tile
n3=curr_tile%nz_cells_tile
curr_tile%npmax_tile=n1*n2*n3*curr%nppcell
curr_tile%np_tile=0

call allocate_tile_arrays(curr_tile)
end do
end do
end do
end do





!$omp parallel do collapse(3) schedule(runtime) default(none) &
!$omp shared(species_parray, ntilex, ntiley, ntilez, nspecies) &
!$omp private(ix,iy,iz,ispecies,curr,curr_tile)
do iz=1, ntilez
do iy=1, ntiley
do ix=1, ntilex
do ispecies=1, nspecies
curr=>species_parray(ispecies)
curr_tile=>curr%array_of_tiles(ix,iy,iz)

curr_tile%part_x=0.0_num
curr_tile%part_y=0.0_num
curr_tile%part_z=0.0_num
curr_tile%part_ux=0.0_num
curr_tile%part_uy=0.0_num
curr_tile%part_uz=0.0_num
curr_tile%part_ex=0.0_num
curr_tile%part_ey=0.0_num
curr_tile%part_ez=0.0_num
curr_tile%part_bx=0.0_num
curr_tile%part_by=0.0_num
curr_tile%part_bz=0.0_num
curr_tile%weight=0.0_num
end do
end do
end do
end do
!$omp end parallel do

end subroutine init_tile_arrays
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine load_particles
use interf_add_particle_to_species
use interf_add_particle_to_species
use tiling
implicit none
type(particle_species), pointer :: curr
integer(kind=4) :: ispecies, l, k, j, ipart
integer(kind=4) :: jmin, jmax, kmin, kmax, lmin, lmax
real(kind=8) :: partx, party, partz, partux, partuy, partuz, partw
real(kind=8) :: phi, th, v
integer(kind=4) :: err, npart
real(kind=8), dimension(6) :: rng=0


if (pdistr .eq. 1) then
do ispecies=1,nspecies
curr=>species_parray(ispecies)
jmin = nint(max(curr%x_min-x_min_local,0.0_num)/dx)
jmax = nint(min(curr%x_max-x_min_local,x_max_local-x_min_local)/dx)
kmin = nint(max(curr%y_min-y_min_local,0.0_num)/dy)
kmax = nint(min(curr%y_max-y_min_local,y_max_local-y_min_local)/dy)
lmin = nint(max(curr%z_min-z_min_local,0.0_num)/dz)
lmax = nint(min(curr%z_max-z_min_local,z_max_local-z_min_local)/dz)
do l=lmin,lmax-1
do k=kmin,kmax-1
do j=jmin,jmax-1
do ipart=1,curr%nppcell

partx = x_min_local+j*dx+dx/curr%nppcell*(ipart-1)
party = y_min_local+k*dy+dy/curr%nppcell*(ipart-1)
partz = z_min_local+l*dz+dz/curr%nppcell*(ipart-1)
partw = nc*dx*dy*dz/(curr%nppcell)

call random_number(rng(1:3))
v=max(1e-10_num,rng(1))
th=2*pi*rng(2)
phi=2*pi*rng(3)
partux= curr%vdrift_x + curr%vth_x*sqrt(-2.*log(v))*cos(th)*cos(phi)
partuy= curr%vdrift_y + curr%vth_y*sqrt(-2.*log(v))*cos(th)*sin(phi)
partuz= curr%vdrift_z + curr%vth_z*sqrt(-2.*log(v))*sin(th)

call add_particle_to_species(curr, partx, party, partz,partux, partuy, partuz, partw)
end do
end do
end do
end do
end do
endif

if (pdistr .eq. 2) then
do ispecies=1,nspecies
curr=>species_parray(ispecies)
do j=0,nx-1
do k=0,ny-1
do l=0,nz-1
do ipart=1,curr%nppcell
call random_number(rng(1:6))

partx = x_min_local+min(rng(1),0.999)*(x_max_local-x_min_local)
party = y_min_local+min(rng(2),0.999)*(y_max_local-y_min_local)
partz = z_min_local+min(rng(3),0.999)*(z_max_local-z_min_local)
partw = nc*dx*dy*dz/(curr%nppcell)

v=max(1e-10_num,rng(4))
th=2*pi*rng(5)
phi=2*pi*rng(6)
partux= curr%vdrift_x + curr%vth_x*sqrt(-2.*log(v))*cos(th)*cos(phi)
partuy= curr%vdrift_y + curr%vth_y*sqrt(-2.*log(v))*cos(th)*sin(phi)
partuz= curr%vdrift_z + curr%vth_z*sqrt(-2.*log(v))*sin(th)

call add_particle_to_species(curr, partx, party, partz,partux, partuy, partuz, partw)
end do
end do
end do
end do
end do
endif


ntot=0
do ispecies=1,nspecies
curr=>species_parray(ispecies)
call mpi_allreduce(curr%species_npart,npart,1_isp, mpi_integer,mpi_sum,comm, err)
ntot=ntot+npart
if (rank .eq. 0) then
write (0,*) 'loaded npart = ', npart,' particles of species ',trim(adjustl(curr%name))
end if
end do

return
end subroutine load_particles
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine resize_particle_arrays(curr, old_size, new_size)
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use interf_resize_array_real
use tiling
implicit none

type(particle_tile), pointer, intent(in out) :: curr
integer :: old_size, new_size

curr%npmax_tile=new_size
call resize_array_real(curr%part_x, old_size, new_size)
call resize_array_real(curr%part_y, old_size, new_size)
call resize_array_real(curr%part_z, old_size, new_size)
call resize_array_real(curr%part_ux, old_size, new_size)
call resize_array_real(curr%part_uy, old_size, new_size)
call resize_array_real(curr%part_uz, old_size, new_size)
call resize_array_real(curr%weight, old_size, new_size)
call resize_array_real(curr%part_ex, old_size, new_size)
call resize_array_real(curr%part_ey, old_size, new_size)
call resize_array_real(curr%part_ez, old_size, new_size)
call resize_array_real(curr%part_bx, old_size, new_size)
call resize_array_real(curr%part_by, old_size, new_size)
call resize_array_real(curr%part_bz, old_size, new_size)
end subroutine resize_particle_arrays
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine resize_array_real(arr, old_size, new_size)
use tiling
implicit none
real(kind=8), dimension(:),pointer, intent(in out) :: arr
real(kind=8), dimension(:),pointer :: temp
integer :: old_size, new_size

allocate(temp(1:new_size))

temp(1:old_size)=arr(1:old_size)
deallocate(arr)
allocate(arr(1:new_size))
arr(1:old_size) = temp(1:old_size)
deallocate(temp)
end subroutine
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine point_to_tile(ispecies, ix, iy, iz)
use tiling
use python_pointers
implicit none
integer(idp), intent(in) :: ix,iy,iz,ispecies
type(particle_species), pointer  :: currsp
type(particle_tile), pointer ::curr_tile

currsp=> species_parray(ispecies)
curr_tile=>currsp%array_of_tiles(ix,iy,iz)
partx=>curr_tile%part_x
party=>curr_tile%part_y
partz=>curr_tile%part_z
partux=>curr_tile%part_ux
partuy=>curr_tile%part_uy
partuz=>curr_tile%part_uz
partw=>curr_tile%weight

end subroutine point_to_tile
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine set_particle_species_properties(nsp,sname,mss,chrg,nppc,xsmin,ysmin,zsmin,xsmax,ysmax,zsmax,vdxs,vdys,vdzs,vthxs,vthys,vthzs)
use tiling
implicit none
integer(idp), intent(in) :: nsp, nppc
real(kind=8), intent(in) :: mss, chrg,xsmin,ysmin,zsmin,xsmax,ysmax,zsmax,vdxs,vdys,vdzs,vthxs,vthys,vthzs
character(len=*), intent(in) :: sname
type(particle_species), pointer  :: currsp

currsp=> species_parray(nsp)
currsp%charge=chrg
currsp%mass=mss
currsp%x_min=xsmin
currsp%y_min=ysmin
currsp%z_min=zsmin
currsp%x_max=xsmax
currsp%y_max=ysmax
currsp%z_max=zsmax
currsp%vdrift_x=vdxs
currsp%vdrift_y=vdys
currsp%vdrift_z=vdzs
currsp%vth_x=vthxs
currsp%vth_y=vthys
currsp%vth_z=vthzs
currsp%nppcell=nppc


end subroutine set_particle_species_properties
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine push_particles
use particles
use constants
use fields
use params
use shared_data
use tiling
implicit none
integer(idp) :: ispecies, ix, iy, iz, count
integer(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
type(particle_species), pointer :: curr
type(particle_tile), pointer :: curr_tile
real(kind=8) :: tdeb, tend
integer(idp) :: nxc, nyc, nzc, ipmin,ipmax, np,ip
integer(idp) :: nblk=1000000

tdeb=mpi_wtime()
!$omp parallel do collapse(3) schedule(runtime) default(none) &
!$omp shared(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$omp nxjguards,nyjguards,nzjguards,ex,ey,ez,bx,by,bz,dx,dy,dz,dt,nblk) &
!$omp private(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$omp lmax,nxc,nyc,nzc, ipmin,ipmax,ip,np)
do iz=1, ntilez
do iy=1, ntiley
do ix=1, ntilex
do ispecies=1, nspecies


curr=>species_parray(ispecies)
curr_tile=>curr%array_of_tiles(ix,iy,iz)
count=curr_tile%np_tile
jmin=curr_tile%nx_tile_min-nxjguards
jmax=curr_tile%nx_tile_max+nxjguards
kmin=curr_tile%ny_tile_min-nyjguards
kmax=curr_tile%ny_tile_max+nyjguards
lmin=curr_tile%nz_tile_min-nzjguards
lmax=curr_tile%nz_tile_max+nzjguards
nxc=curr_tile%nx_cells_tile
nyc=curr_tile%ny_cells_tile
nzc=curr_tile%nz_cells_tile
curr_tile%part_ex = 0.0_num
curr_tile%part_ey = 0.0_num
curr_tile%part_ez = 0.0_num
curr_tile%part_bx=0.0_num
curr_tile%part_by=0.0_num
curr_tile%part_bz=0.0_num

do ip=1,count,nblk
np=min(count-ip+1,nblk)
ipmin=ip
ipmax=ip+np-1

call gete3d_energy_conserving_1_1_1(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),curr_tile%part_z(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax),curr_tile%part_ey(ipmin:ipmax),curr_tile%part_ez(ipmin:ipmax),curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,curr_tile%z_grid_tile_min, dx,dy,dz,curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,nzjguards,ex(jmin:jmax,kmin:kmax,lmin:lmax),ey(jmin:jmax,kmin:kmax,lmin:lmax),ez(jmin:jmax,kmin:kmax,lmin:lmax))

call getb3d_energy_conserving_1_1_1(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),curr_tile%part_z(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax),curr_tile%part_by(ipmin:ipmax),curr_tile%part_bz(ipmin:ipmax),curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,curr_tile%z_grid_tile_min, dx,dy,dz,curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,nzjguards,bx(jmin:jmax,kmin:kmax,lmin:lmax),by(jmin:jmax,kmin:kmax,lmin:lmax),bz(jmin:jmax,kmin:kmax,lmin:lmax))

call bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax),curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dt*0.5_num)

call epush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),curr_tile%part_uz(ipmin:ipmax), curr_tile%part_ex(ipmin:ipmax), curr_tile%part_ey(ipmin:ipmax),curr_tile%part_ez(ipmin:ipmax), curr%charge,curr%mass,dt)

call bpush_v(np,curr_tile%part_ux(ipmin:ipmax), curr_tile%part_uy(ipmin:ipmax),curr_tile%part_uz(ipmin:ipmax), curr_tile%part_bx(ipmin:ipmax), curr_tile%part_by(ipmin:ipmax),curr_tile%part_bz(ipmin:ipmax), curr%charge,curr%mass,dt*0.5_num)

call pushxyz(np,curr_tile%part_x(ipmin:ipmax),curr_tile%part_y(ipmin:ipmax),curr_tile%part_z(ipmin:ipmax), curr_tile%part_ux(ipmin:ipmax),curr_tile%part_uy(ipmin:ipmax),curr_tile%part_uz(ipmin:ipmax),dt)
end do
end do
end do
end do
end do
!$omp end parallel do
tend=mpi_wtime()
pushtime=pushtime+(tend-tdeb)
end subroutine push_particles
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine pushxyz(np,xp,yp,zp,uxp,uyp,uzp,dt)

use constants
use omp_lib
implicit none
integer(idp)   :: np
real(kind=8) :: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np)
real(kind=8) :: dt,gaminv,clghtisq,usq
integer(idp)  :: ip

clghtisq = 1.0_num/clight**2
!!$omp parallel do private(ip, usq, gaminv)
do ip=1,np
usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
gaminv = 1.0_num/sqrt(1.0_num + usq)
xp(ip) = xp(ip) + uxp(ip)*gaminv*dt
yp(ip) = yp(ip) + uyp(ip)*gaminv*dt
zp(ip) = zp(ip) + uzp(ip)*gaminv*dt
enddo
!!$omp end parallel do

return
end subroutine pushxyz
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine epush_v(np,uxp,uyp,uzp,ex,ey,ez,q,m,dt)


use constants
implicit none
integer(idp) :: np
real(kind=8):: uxp(np),uyp(np),uzp(np)
real(kind=8):: ex(np),ey(np),ez(np)
real(kind=8):: q,m,dt

integer(idp) :: ip
real(kind=8):: const

const = q*dt/m
!!$omp parallel do private(ip)
do ip=1,np
uxp(ip) = uxp(ip) + ex(ip)*const
uyp(ip) = uyp(ip) + ey(ip)*const
uzp(ip) = uzp(ip) + ez(ip)*const
enddo
!!$omp end parallel do

return
end subroutine epush_v
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine bpush_v(np,uxp,uyp,uzp,bx,by,bz,q,m,dt)


use constants
implicit none
integer(idp)   :: np
real(kind=8) :: uxp(np), uyp(np), uzp(np)
real(kind=8) :: bx(np), by(np), bz(np)
real(kind=8) :: q,m,dt,gaminv
integer(idp)   :: ip
real(kind=8) :: const,clghtisq,sx,sy,sz,tx,ty,tz,tsqi,uxppr,uyppr,uzppr,usq

const = q*dt*0.5_num/m
clghtisq = 1.0_num/clight**2

!!$omp parallel do private(ip, tx, ty, tz, tsqi, sx, sy, sz, uxppr, uyppr, uzppr, usq, gaminv)
do ip=1,np
usq = (uxp(ip)**2 + uyp(ip)**2+ uzp(ip)**2)*clghtisq
gaminv = 1.0_num/sqrt(1.0_num + usq)
tx = gaminv*bx(ip)*const
ty = gaminv*by(ip)*const
tz = gaminv*bz(ip)*const
tsqi = 2.0_num/(1.0_num + tx**2 + ty**2 + tz**2)
sx = tx*tsqi
sy = ty*tsqi
sz = tz*tsqi
uxppr = uxp(ip) + uyp(ip)*tz - uzp(ip)*ty
uyppr = uyp(ip) + uzp(ip)*tx - uxp(ip)*tz
uzppr = uzp(ip) + uxp(ip)*ty - uyp(ip)*tx
uxp(ip) = uxp(ip) + uyppr*sz - uzppr*sy
uyp(ip) = uyp(ip) + uzppr*sx - uxppr*sz
uzp(ip) = uzp(ip) + uxppr*sy - uyppr*sx
enddo
!!$omp end parallel do

return

end subroutine bpush_v
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_currents_on_grid_jxjyjz
use particles
use constants
use fields
use params
use shared_data
use tiling
use omp_lib
implicit none
integer(idp) :: ispecies, ix, iy, iz, count
integer(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
integer(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
type(particle_species), pointer :: curr
type(particle_tile), pointer :: curr_tile
real(kind=8) :: tdeb, tend
real(kind=8), dimension(:,:,:), pointer :: jx_tile,jy_tile,jz_tile
integer(idp) :: nxc, nyc, nzc

jx = 0.0_num
jy = 0.0_num
jz = 0.0_num

tdeb=mpi_wtime()
!$omp parallel default(none) &
!$omp shared(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguards,nyjguards,nzjguards,dx,dy,dz,dt,jx,jy,jz) &
!$omp private(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$omp lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,jx_tile,jy_tile,jz_tile,nxc,nyc,nzc)

!$omp do collapse(3) schedule(runtime)
do iz=1,ntilez
do iy=1,ntiley
do ix=1,ntilex
do ispecies=1, nspecies
curr => species_parray(ispecies)
curr_tile=>curr%array_of_tiles(ix,iy,iz)
count=curr_tile%np_tile
jmin=curr_tile%nx_tile_min
jmax=curr_tile%nx_tile_max
kmin=curr_tile%ny_tile_min
kmax=curr_tile%ny_tile_max
lmin=curr_tile%nz_tile_min
lmax=curr_tile%nz_tile_max
nxc=curr_tile%nx_cells_tile
nyc=curr_tile%ny_cells_tile
nzc=curr_tile%nz_cells_tile
curr_tile%jxtile = 0.0_num
curr_tile%jytile = 0.0_num
curr_tile%jztile = 0.0_num

call depose_jxjyjz_scalar_1_1_1(curr_tile%jxtile,curr_tile%jytile,curr_tile%jztile,count,curr_tile%part_x(1:count),curr_tile%part_y(1:count),curr_tile%part_z(1:count),curr_tile%part_ux(1:count),curr_tile%part_uy(1:count),curr_tile%part_uz(1:count),curr_tile%weight(1:count),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,curr_tile%z_grid_tile_min,dt,dx,dy,dz,curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,nzjguards)

jx(jmin:jmax,kmin:kmax,lmin:lmax) = jx(jmin:jmax,kmin:kmax,lmin:lmax) + curr_tile%jxtile(0:nxc,0:nyc,0:nzc)
jy(jmin:jmax,kmin:kmax,lmin:lmax) = jy(jmin:jmax,kmin:kmax,lmin:lmax) + curr_tile%jytile(0:nxc,0:nyc,0:nzc)
jz(jmin:jmax,kmin:kmax,lmin:lmax) = jz(jmin:jmax,kmin:kmax,lmin:lmax) + curr_tile%jztile(0:nxc,0:nyc,0:nzc)
end do
end do
end do
end do
!$omp end do

!$omp do collapse(3) schedule(runtime)
do iz=1,ntilez
do iy=1,ntiley
do ix=1,ntilex
do ispecies=1, nspecies
curr => species_parray(ispecies)
curr_tile=>curr%array_of_tiles(ix,iy,iz)
count=curr_tile%np_tile
jmin=curr_tile%nx_tile_min
jmax=curr_tile%nx_tile_max
kmin=curr_tile%ny_tile_min
kmax=curr_tile%ny_tile_max
lmin=curr_tile%nz_tile_min
lmax=curr_tile%nz_tile_max
jminc=jmin-nxjguards
jmaxc=jmax+nxjguards
kminc=kmin-nyjguards
kmaxc=kmax+nyjguards
lminc=lmin-nzjguards
lmaxc=lmax+nzjguards
nxc=curr_tile%nx_cells_tile
nyc=curr_tile%ny_cells_tile
nzc=curr_tile%nz_cells_tile



jx(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jx(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+curr_tile%jxtile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
jx(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jx(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+curr_tile%jxtile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)

jx(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jx(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+curr_tile%jxtile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
jx(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jx(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+curr_tile%jxtile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)

jx(jmin:jmax,kmin:kmax,lminc:lmin-1) = jx(jmin:jmax,kmin:kmax,lminc:lmin-1)+curr_tile%jxtile(0:nxc, 0:nyc,-nzjguards:-1)
jx(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jx(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+curr_tile%jxtile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)


jy(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jy(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+curr_tile%jytile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
jy(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jy(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+curr_tile%jytile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)

jy(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jy(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+curr_tile%jytile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
jy(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jy(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+curr_tile%jytile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)

jy(jmin:jmax,kmin:kmax,lminc:lmin-1) = jy(jmin:jmax,kmin:kmax,lminc:lmin-1)+curr_tile%jytile(0:nxc, 0:nyc,-nzjguards:-1)
jy(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jy(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+curr_tile%jytile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)


jz(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jz(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+curr_tile%jztile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
jz(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jz(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+curr_tile%jztile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)

jz(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jz(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+curr_tile%jztile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
jz(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jz(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+curr_tile%jztile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)

jz(jmin:jmax,kmin:kmax,lminc:lmin-1) = jz(jmin:jmax,kmin:kmax,lminc:lmin-1)+curr_tile%jztile(0:nxc, 0:nyc,-nzjguards:-1)
jz(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jz(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+curr_tile%jztile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)
end do
end do
end do
end do
!$omp end do
!$omp end parallel
tend=mpi_wtime()
pushtime=pushtime+(tend-tdeb)

end subroutine depose_currents_on_grid_jxjyjz
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_jxjyjz_scalar_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin,dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp, w
real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
real(kind=8) :: dxi,dyi,dzi,xint,yint,zint,oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
real(kind=8) :: x,y,z,xmid,ymid,zmid,vx,vy,vz,invvol, dts2dx, dts2dy, dts2dz
real(kind=8) :: wq, wqx, wqy, wqz, gaminv, usq, clightsq
real(kind=8), dimension(2) :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0
integer(idp) :: j,k,l,ip
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
invvol = dxi*dyi*dzi
dts2dx = 0.5_num*dt*dxi
dts2dy = 0.5_num*dt*dyi
dts2dz = 0.5_num*dt*dzi
clightsq = 1.0_num/clight**2
sx=0.0_num;sy=0.0_num;sz=0.0_num;
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;


do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi


usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clightsq
gaminv = 1.0_num/sqrt(1.0_num + usq)
vx = uxp(ip)*gaminv
vy = uyp(ip)*gaminv
vz = uzp(ip)*gaminv


wq=q*w(ip)
wqx=wq*invvol*vx
wqy=wq*invvol*vy
wqz=wq*invvol*vz


xmid=x-dts2dx*vx
ymid=y-dts2dy*vy
zmid=z-dts2dz*vz


j=floor(xmid)
k=floor(ymid)
l=floor(zmid)


xint = xmid-j
yint = ymid-k
zint = zmid-l
sx( 0) = 1.0_num-xint
sx( 1) = xint
sy( 0) = 1.0_num-yint
sy( 1) = yint
sz( 0) = 1.0_num-zint
sz( 1) = zint


xint = xmid-j-0.5_num
yint = ymid-k-0.5_num
zint = zmid-l-0.5_num
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
sz0( 0) = 1.0_num-zint
sz0( 1) = zint



jx(j,k,l)      = jx(j,k,l)+sx0(0)*sy(0)*sz(0)*wqx
jx(j+1,k,l)    = jx(j+1,k,l)+sx0(1)*sy(0)*sz(0)*wqx
jx(j,k+1,l)    = jx(j,k+1,l)+sx0(0)*sy(1)*sz(0)*wqx
jx(j+1,k+1,l)  = jx(j+1,k+1,l)+sx0(1)*sy(1)*sz(0)*wqx
jx(j,k,l+1)    = jx(j,k,l+1)+sx0(0)*sy(0)*sz(1)*wqx
jx(j+1,k,l+1)  = jx(j+1,k,l+1)+sx0(1)*sy(0)*sz(1)*wqx
jx(j,k+1,l+1)  = jx(j,k+1,l+1)+sx0(0)*sy(1)*sz(1)*wqx
jx(j+1,k+1,l+1)= jx(j+1,k+1,l+1)+sx0(1)*sy(1)*sz(1)*wqx


jy(j,k,l)      = jy(j,k,l)+sx(0)*sy0(0)*sz(0)*wqy
jy(j+1,k,l)    = jy(j+1,k,l)+sx(1)*sy0(0)*sz(0)*wqy
jy(j,k+1,l)    = jy(j,k+1,l)+sx(0)*sy0(1)*sz(0)*wqy
jy(j+1,k+1,l)  = jy(j+1,k+1,l)+sx(1)*sy0(1)*sz(0)*wqy
jy(j,k,l+1)    = jy(j,k,l+1)+sx(0)*sy0(0)*sz(1)*wqy
jy(j+1,k,l+1)  = jy(j+1,k,l+1)+sx(1)*sy0(0)*sz(1)*wqy
jy(j,k+1,l+1)  = jy(j,k+1,l+1)+sx(0)*sy0(1)*sz(1)*wqy
jy(j+1,k+1,l+1)= jy(j+1,k+1,l+1)+sx(1)*sy0(1)*sz(1)*wqy


jz(j,k,l)      = jz(j,k,l)+sx(0)*sy(0)*sz0(0)*wqz
jz(j+1,k,l)    = jz(j+1,k,l)+sx(1)*sy(0)*sz0(0)*wqz
jz(j,k+1,l)    = jz(j,k+1,l)+sx(0)*sy(1)*sz0(0)*wqz
jz(j+1,k+1,l)  = jz(j+1,k+1,l)+sx(1)*sy(1)*sz0(0)*wqz
jz(j,k,l+1)    = jz(j,k,l+1)+sx(0)*sy(0)*sz0(1)*wqz
jz(j+1,k,l+1)  = jz(j+1,k,l+1)+sx(1)*sy(0)*sz0(1)*wqz
jz(j,k+1,l+1)  = jz(j,k+1,l+1)+sx(0)*sy(1)*sz0(1)*wqz
jz(j+1,k+1,l+1)= jz(j+1,k+1,l+1)+sx(1)*sy(1)*sz0(1)*wqz
end do
return
end subroutine depose_jxjyjz_scalar_1_1_1
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_jxjyjz_esirkepov_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin,dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)

use omp_lib
use constants
implicit none
integer(idp):: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
real(kind=8), dimension(:,:,:), pointer:: jx1, jy1, jz1
real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp, w
real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
real(kind=8), dimension(:,:,:), pointer :: sdx,sdy,sdz
real(kind=8) :: clghtisq,usq,gaminv,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq,dtsdx0,dtsdy0,dtsdz0
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0
real(kind=8), dimension(:), pointer:: sx, sx0, dsx
real(kind=8), dimension(:), pointer :: sy, sy0, dsy
real(kind=8), dimension(:), pointer :: sz, sz0, dsz
integer(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc,ixmin, ixmax, iymin, iymax, izmin, izmax


dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
dtsdx0 = dt*dxi
dtsdy0 = dt*dyi
dtsdz0 = dt*dzi
invvol = 1.0_num/(dx*dy*dz)
invdtdx = 1.0_num/(dt*dy*dz)
invdtdy = 1.0_num/(dt*dx*dz)
invdtdz = 1.0_num/(dt*dx*dy)
allocate(sdx(-1:2,-1:2,-1:2),sdy(-1:2,-1:2,-1:2),sdz(-1:2,-1:2,-1:2))
allocate(sx(-1:2), sx0(-1:2), dsx(-1:2))
allocate(sy(-1:2), sy0(-1:2), dsy(-1:2))
allocate(sz(-1:2), sz0(-1:2), dsz(-1:2))
clghtisq = 1.0_num/clight**2
dtsdz0 = dt*dzi
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num
allocate(jx1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard),jy1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard),jz1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
!!$omp parallel default(none) private(ip,x,y,z,usq,vx,vy,vz,gaminv,xold,yold,zold, &
!!$omp wq,wqx,wqy,wqz,iixp0,ijxp0,ikxp0, xint,yint,zint, oxint,xintsq, oxintsq,dix,diy,diz, &
!!$omp dsx, dsy, dsz, oyint,yintsq, oyintsq, ozint,zintsq, ozintsq,ixmin, ixmax, iymin, iymax, izmin, izmax,  &
!!$omp k,j,i,kc,jc,ic, iixp, ijxp, ikxp,sx,sy,sz, sx0,sy0,sz0,sdx,sdy,sdz,jx1,jy1,jz1) &
!!$omp shared(np,xp,yp,zp,uxp,uyp,uzp,w,dxi,dyi,dzi,invdtdx,invdtdy,invdtdz,xmin,ymin,zmin,clghtisq,dtsdx0,dtsdy0,dtsdz0,q,jx,jy,jz)
jx1=0.0_num; jy1=0.0_num; jz1=0.0_num
!!$omp do
do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi

usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq
gaminv = 1.0_num/sqrt(1.0_num + usq)
vx = uxp(ip)*gaminv
vy = uyp(ip)*gaminv
vz = uzp(ip)*gaminv

xold=x-dtsdx0*vx
yold=y-dtsdy0*vy
zold=z-dtsdz0*vz

wq=q*w(ip)
wqx = wq*invdtdx
wqy = wq*invdtdy
wqz = wq*invdtdz

iixp0=floor(x)
ijxp0=floor(y)
ikxp0=floor(z)

xint=x-iixp0
yint=y-ijxp0
zint=z-ikxp0


sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
sz0( 0) = 1.0_num-zint
sz0( 1) = zint

iixp=floor(xold)
ijxp=floor(yold)
ikxp=floor(zold)

xint = xold-iixp
yint = yold-ijxp
zint = zold-ikxp

dix = iixp-iixp0
diy = ijxp-ijxp0
diz = ikxp-ikxp0

sx=0.0_num;sy=0.0_num;sz=0.0_num

sx( 0+dix) = 1.0_num-xint
sx( 1+dix) = xint
sy( 0+diy) = 1.0_num-yint
sy( 1+diy) = yint
sz( 0+diz) = 1.0_num-zint
sz( 1+diz) = zint

dsx = sx - sx0
dsy = sy - sy0
dsz = sz - sz0

ixmin = min(0,dix)
ixmax = max(0,dix)+1
iymin = min(0,diy)
iymax = max(0,diy)+1
izmin = min(0,diz)
izmax = max(0,diz)+1


do k=izmin, izmax
do j=iymin, iymax
do i=ixmin, ixmax
ic = iixp0+i
jc = ijxp0+j
kc = ikxp0+k
if(i<ixmax) then
sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) +(0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
if (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
jx1(ic,jc,kc) = jx1(ic,jc,kc) + sdx(i,j,k)
end if
if(j<iymax) then
sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) +(0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
if (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
jy1(ic,jc,kc) = jy1(ic,jc,kc) + sdy(i,j,k)
end if
if(k<izmax) then
sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) +(0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
if (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
jz1(ic,jc,kc) = jz1(ic,jc,kc) + sdz(i,j,k)
end if
end do
end do
end do
end do
!!$omp end do
!!$omp critical
jx=jx+jx1
jy=jy+jy1
jz=jz+jz1
!!$omp end critical
!!$omp end parallel
deallocate(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz,jx1,jy1,jz1)
return
end subroutine depose_jxjyjz_esirkepov_1_1_1
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_jxjyjz_esirkepov_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin,dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz,l_particles_weight,l4symtry)


use constants
implicit none
integer(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
real(kind=8), dimension(:,:,:), pointer:: jx1, jy1, jz1
real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp, w
real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
real(kind=8), dimension(:,:,:), pointer :: sdx,sdy,sdz
real(kind=8) :: clghtisq,usq,gaminv,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz,s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq,dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0
real(kind=8), dimension(:), pointer :: sx, sx0, dsx
real(kind=8), dimension(:), pointer :: sy, sy0, dsy
real(kind=8), dimension(:), pointer :: sz, sz0, dsz
integer(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc,ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx, ndtody, ndtodz,xl,xu,yl,yu,zl,zu
logical :: l_particles_weight,l4symtry


ndtodx = int(clight*dt/dx)
ndtody = int(clight*dt/dy)
ndtodz = int(clight*dt/dz)
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
dtsdx0 = dt*dxi
dtsdy0 = dt*dyi
dtsdz0 = dt*dzi
dts2dx0 = 0.5_num*dtsdx0
dts2dy0 = 0.5_num*dtsdy0
dts2dz0 = 0.5_num*dtsdz0
invvol = 1.0_num/(dx*dy*dz)
invdtdx = 1.0_num/(dt*dy*dz)
invdtdy = 1.0_num/(dt*dx*dz)
invdtdz = 1.0_num/(dt*dx*dy)

xl = -int(nox/2)-1-ndtodx
xu = int((nox+1)/2)+1+ndtodx
yl = -int(noy/2)-1-ndtody
yu = int((noy+1)/2)+1+ndtody
zl = -int(noz/2)-1-ndtodz
zu = int((noz+1)/2)+1+ndtodz
allocate(sdx(xl:xu,yl:yu,zl:zu),sdy(xl:xu,yl:yu,zl:zu),sdz(xl:xu,yl:yu,zl:zu))
allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
allocate(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))
allocate(jx1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard),jy1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard),jz1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
clghtisq = 1.0_num/clight**2
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num
jx1=0.0_num;jy1=0.0_num;jz1=0.0_num

!$omp parallel private(ip,x,y,z,usq,vx,vy,vz,gaminv,xold,yold,zold,ncells,dtsdx,dtsdy,dtsdz,dts2dx,dts2dy,dts2dz, &
!$omp icell, wq,wqx,wqy,wqz,iixp0,ijxp0,ikxp0, xint,yint,zint, oxint,xintsq, oxintsq,dix,diy,diz, &
!$omp dsx, dsy, dsz, oyint,yintsq, oyintsq, ozint,zintsq, ozintsq,ixmin, ixmax, iymin, iymax, izmin, izmax,  &
!$omp k,j,i,kc,jc,ic, iixp, ijxp, ikxp,sx,sy,sz) firstprivate(sx0,sy0,sz0,sdx,sdy,sdz,jx1,jy1,jz1)
!$omp do
do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi

usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq
gaminv = 1.0_num/sqrt(1.0_num + usq)
vx = uxp(ip)*gaminv
vy = uyp(ip)*gaminv
vz = uzp(ip)*gaminv

xold=x-dtsdx0*vx
yold=y-dtsdy0*vy
zold=z-dtsdz0*vz

if (l4symtry) then
x=abs(x)
y=abs(y)
xold=abs(xold)
yold=abs(yold)
vx = (x-xold)/dtsdx0
vy = (y-yold)/dtsdy0
end if

ncells = 1

dtsdx = dtsdx0/ncells
dtsdy = dtsdy0/ncells
dtsdz = dtsdz0/ncells
dts2dx = dts2dx0/ncells
dts2dy = dts2dy0/ncells
dts2dz = dts2dz0/ncells

x=xold
y=yold
z=zold

do icell = 1,ncells
xold = x
yold = y
zold = z
x = x+dtsdx*vx
y = y+dtsdy*vy
z = z+dtsdz*vz


if (l_particles_weight) then
wq=q*w(ip)
else
wq=q*w(1)
end if
wqx = wq*invdtdx
wqy = wq*invdtdy
wqz = wq*invdtdz



if (nox==2*(nox/2)) then
iixp0=nint(x)
else
iixp0=floor(x)
end if
if (noy==2*(noy/2)) then
ijxp0=nint(y)
else
ijxp0=floor(y)
end if
if (noz==2*(noz/2)) then
ikxp0=nint(z)
else
ikxp0=floor(z)
end if

xint=x-iixp0
yint=y-ijxp0
zint=z-ikxp0


select case(nox)
case(0)
sx0( 0) = 1.0_num
case(1)
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
case(2)
xintsq = xint*xint
sx0(-1) = 0.5_num*(0.5_num-xint)**2
sx0( 0) = 0.75_num-xintsq
sx0( 1) = 0.5_num*(0.5_num+xint)**2
case(3)
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx0(-1) = onesixth*oxintsq*oxint
sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx0( 2) = onesixth*xintsq*xint
end select

select case(noy)
case(0)
sy0( 0) = 1.0_num
case(1)
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
case(2)
yintsq = yint*yint
sy0(-1) = 0.5_num*(0.5_num-yint)**2
sy0( 0) = 0.75_num-yintsq
sy0( 1) = 0.5_num*(0.5_num+yint)**2
case(3)
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy0(-1) = onesixth*oyintsq*oyint
sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy0( 2) = onesixth*yintsq*yint
end select

select case(noz)
case(0)
sz0( 0) = 1.0_num
case(1)
sz0( 0) = 1.0_num-zint
sz0( 1) = zint
case(2)
zintsq = zint*zint
sz0(-1) = 0.5_num*(0.5_num-zint)**2
sz0( 0) = 0.75_num-zintsq
sz0( 1) = 0.5_num*(0.5_num+zint)**2
case(3)
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz0(-1) = onesixth*ozintsq*ozint
sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz0( 2) = onesixth*zintsq*zint
end select



if (nox==2*(nox/2)) then
iixp=nint(xold)
else
iixp=floor(xold)
end if
if (noy==2*(noy/2)) then
ijxp=nint(yold)
else
ijxp=floor(yold)
end if
if (noz==2*(noz/2)) then
ikxp=nint(zold)
else
ikxp=floor(zold)
end if


xint = xold-iixp
yint = yold-ijxp
zint = zold-ikxp


dix = iixp-iixp0
diy = ijxp-ijxp0
diz = ikxp-ikxp0


sx=0.0_num;sy=0.0_num;sz=0.0_num


select case(nox)
case(0)
sx( 0+dix) = 1.0_num
case(1)
sx( 0+dix) = 1.0_num-xint
sx( 1+dix) = xint
case(2)
xintsq = xint*xint
sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
sx( 0+dix) = 0.75_num-xintsq
sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
case(3)
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx(-1+dix) = onesixth*oxintsq*oxint
sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx( 2+dix) = onesixth*xintsq*xint
end select

select case(noy)
case(0)
sy( 0+diy) = 1.0_num
case(1)
sy( 0+diy) = 1.0_num-yint
sy( 1+diy) = yint
case(2)
yintsq = yint*yint
sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
sy( 0+diy) = 0.75_num-yintsq
sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
case(3)
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy(-1+diy) = onesixth*oyintsq*oyint
sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy( 2+diy) = onesixth*yintsq*yint
end select

select case(noz)
case(0)
sz( 0+diz) = 1.0_num
case(1)
sz( 0+diz) = 1.0_num-zint
sz( 1+diz) = zint
case(2)
zintsq = zint*zint
sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
sz( 0+diz) = 0.75_num-zintsq
sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
case(3)
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz(-1+diz) = onesixth*ozintsq*ozint
sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz( 2+diz) = onesixth*zintsq*zint
end select


dsx = sx - sx0
dsy = sy - sy0
dsz = sz - sz0


ixmin = min(0,dix)-int(nox/2)
ixmax = max(0,dix)+int((nox+1)/2)
iymin = min(0,diy)-int(noy/2)
iymax = max(0,diy)+int((noy+1)/2)
izmin = min(0,diz)-int(noz/2)
izmax = max(0,diz)+int((noz+1)/2)


do k=izmin, izmax
do j=iymin, iymax
do i=ixmin, ixmax
ic = iixp0+i
jc = ijxp0+j
kc = ikxp0+k
if(i<ixmax) then
sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) +(0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
if (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
jx1(ic,jc,kc) = jx1(ic,jc,kc) + sdx(i,j,k)
end if
if(j<iymax) then
sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) +(0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
if (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
jy1(ic,jc,kc) = jy1(ic,jc,kc) + sdy(i,j,k)
end if
if(k<izmax) then
sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) +(0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
if (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
jz1(ic,jc,kc) = jz1(ic,jc,kc) + sdz(i,j,k)
end if
end do
end do
end do
end do
end do
!$omp end do

!$omp critical
jx=jx+jx1
jy=jy+jy1
jz=jz+jz1
!$omp end critical
!$omp end parallel
deallocate(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz,jx1,jy1,jz1)

return
end subroutine depose_jxjyjz_esirkepov_n
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine gete3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,exg,eyg,ezg)


use omp_lib
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax,ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint,xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
real(kind=8), dimension(0:1) :: sx
real(kind=8), dimension(0:1) :: sy
real(kind=8), dimension(0:1) :: sz
real(kind=8), dimension(:), pointer :: sx0,sy0,sz0
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
allocate(sx0(0:1),sy0(0:1),sz0(0:1))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$omp parallel do private(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$omp sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi


j=floor(x)
j0=floor(x-0.5_num)
k=floor(y)
k0=floor(y-0.5_num)
l=floor(z)
l0=floor(z-0.5_num)
xint=x-j
yint=y-k
zint=z-l


sx( 0) = 1.0_num-xint
sx( 1) = xint
sy( 0) = 1.0_num-yint
sy( 1) = yint
sz( 0) = 1.0_num-zint
sz( 1) = zint
xint=x-0.5_num-j0
yint=y-0.5_num-k0
zint=z-0.5_num-l0
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
sz( 0) = 1.0_num-zint
sz( 1) = zint


ex(ip) = ex(ip) + sx(0)*sy0(0)*sz0(0)*exg(j0,k,l)
ex(ip) = ex(ip) + sx(1)*sy0(0)*sz0(0)*exg(j0+1,k,l)
ex(ip) = ex(ip) + sx(0)*sy0(1)*sz0(0)*exg(j0,k+1,l)
ex(ip) = ex(ip) + sx(1)*sy0(1)*sz0(0)*exg(j0+1,k+1,l)
ex(ip) = ex(ip) + sx(0)*sy0(0)*sz0(1)*exg(j0,k,l+1)
ex(ip) = ex(ip) + sx(1)*sy0(0)*sz0(1)*exg(j0+1,k,l+1)
ex(ip) = ex(ip) + sx(0)*sy0(1)*sz0(1)*exg(j0,k+1,l+1)
ex(ip) = ex(ip) + sx(1)*sy0(1)*sz0(1)*exg(j0+1,k+1,l+1)


ey(ip) = ey(ip) + sx0(0)*sy(0)*sz0(0)*eyg(j,k0,l)
ey(ip) = ey(ip) + sx0(1)*sy(0)*sz0(0)*eyg(j+1,k0,l)
ey(ip) = ey(ip) + sx0(0)*sy(1)*sz0(0)*eyg(j,k0+1,l)
ey(ip) = ey(ip) + sx0(1)*sy(1)*sz0(0)*eyg(j+1,k0+1,l)
ey(ip) = ey(ip) + sx0(0)*sy(0)*sz0(1)*eyg(j,k0,l+1)
ey(ip) = ey(ip) + sx0(1)*sy(0)*sz0(1)*eyg(j+1,k0,l+1)
ey(ip) = ey(ip) + sx0(0)*sy(1)*sz0(1)*eyg(j,k0+1,l+1)
ey(ip) = ey(ip) + sx0(1)*sy(1)*sz0(1)*eyg(j+1,k0+1,l+1)


ez(ip) = ez(ip) + sx0(0)*sy0(0)*sz(0)*ezg(j,k,l0)
ez(ip) = ez(ip) + sx0(1)*sy0(0)*sz(0)*ezg(j+1,k,l0)
ez(ip) = ez(ip) + sx0(0)*sy0(1)*sz(0)*ezg(j,k+1,l0)
ez(ip) = ez(ip) + sx0(1)*sy0(1)*sz(0)*ezg(j+1,k+1,l0)
ez(ip) = ez(ip) + sx0(0)*sy0(0)*sz(1)*ezg(j,k,l0+1)
ez(ip) = ez(ip) + sx0(1)*sy0(0)*sz(1)*ezg(j+1,k,l0+1)
ez(ip) = ez(ip) + sx0(0)*sy0(1)*sz(1)*ezg(j,k+1,l0+1)
ez(ip) = ez(ip) + sx0(1)*sy0(1)*sz(1)*ezg(j+1,k+1,l0+1)
end do
!!$omp end parallel do
deallocate(sx0,sz0)
return
end subroutine gete3d_energy_conserving_1_1_1
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine getb3d_energy_conserving_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,bxg,byg,bzg)


use omp_lib
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(np) :: xp,yp,zp,bx,by,bz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax,ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint,xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
real(kind=8), dimension(0:1) :: sx
real(kind=8), dimension(0:1) :: sy
real(kind=8), dimension(0:1) :: sz
real(kind=8), dimension(:), pointer :: sx0,sy0,sz0
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
allocate(sx0(0:1),sy0(0:1),sz0(0:1))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$omp parallel do private(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$omp sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi


j=floor(x)
j0=floor(x-0.5_num)
k=floor(y)
k0=floor(y-0.5_num)
l=floor(z)
l0=floor(z-0.5_num)
xint=x-j
yint=y-k
zint=z-l


sx( 0) = 1.0_num-xint
sx( 1) = xint
sy( 0) = 1.0_num-yint
sy( 1) = yint
sz( 0) = 1.0_num-zint
sz( 1) = zint
xint=x-0.5_num-j0
yint=y-0.5_num-k0
zint=z-0.5_num-l0
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
sz( 0) = 1.0_num-zint
sz( 1) = zint


bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)


by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)


bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)
end do
!!$omp end parallel do
deallocate(sx0,sz0)
return
end subroutine getb3d_energy_conserving_1_1_1
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz,exg,eyg,ezg,l_lower_order_in_v)

use omp_lib
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
logical :: l4symtry,l_lower_order_in_v
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax,ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint,xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
real(kind=8), dimension(-int(nox/2):int((nox+1)/2)) :: sx
real(kind=8), dimension(-int(noy/2):int((noy+1)/2)) :: sy
real(kind=8), dimension(-int(noz/2):int((noz+1)/2)) :: sz
real(kind=8), dimension(:), pointer :: sx0,sy0,sz0
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
iymin = -int(noy/2)
iymax =  int((noy+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1

if (l_lower_order_in_v) then
ixmin0 = -int((nox-1)/2)
ixmax0 =  int((nox)/2)
iymin0 = -int((noy-1)/2)
iymax0 =  int((noy)/2)
izmin0 = -int((noz-1)/2)
izmax0 =  int((noz)/2)
else
ixmin0 = -int((nox)/2)
ixmax0 =  int((nox+1)/2)
iymin0 = -int((noy)/2)
iymax0 =  int((noy+1)/2)
izmin0 = -int((noz)/2)
izmax0 =  int((noz+1)/2)
end if

allocate(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

signx = 1.0_num
signy = 1.0_num
!$omp parallel do private(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint, &
!$omp   sx,sy,sz,sx0,sy0,sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq,ozint,zintsq,ozintsq)
do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi

if (l_lower_order_in_v) then
if (nox==2*(nox/2)) then
j=nint(x)
j0=floor(x-0.5_num)
else
j=floor(x)
j0=floor(x)
end if
if (noy==2*(noy/2)) then
k=nint(y)
k0=floor(y-0.5_num)
else
k=floor(y)
k0=floor(y)
end if
if (noz==2*(noz/2)) then
l=nint(z)
l0=floor(z-0.5_num)
else
l=floor(z)
l0=floor(z)
end if
else
if (nox==2*(nox/2)) then
j=nint(x)
j0=floor(x)
else
j=floor(x)
j0=floor(x-0.5_num)
end if
if (noy==2*(noy/2)) then
k=nint(y)
k0=floor(y)
else
k=floor(y)
k0=floor(y-0.5_num)
end if
if (noz==2*(noz/2)) then
l=nint(z)
l0=floor(z)
else
l=floor(z)
l0=floor(z-0.5_num)
end if
end if

xint=x-j
yint=y-k
zint=z-l

if (nox==1) then
sx( 0) = 1.0_num-xint
sx( 1) = xint
elseif (nox==2) then
xintsq = xint*xint
sx(-1) = 0.5_num*(0.5_num-xint)**2
sx( 0) = 0.75_num-xintsq
sx( 1) = 0.5_num*(0.5_num+xint)**2
elseif (nox==3) then
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx(-1) = onesixth*oxintsq*oxint
sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx( 2) = onesixth*xintsq*xint
end if

if (noy==1) then
sy( 0) = 1.0_num-yint
sy( 1) = yint
elseif (noy==2) then
yintsq = yint*yint
sy(-1) = 0.5_num*(0.5_num-yint)**2
sy( 0) = 0.75_num-yintsq
sy( 1) = 0.5_num*(0.5_num+yint)**2
elseif (noy==3) then
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy(-1) = onesixth*oyintsq*oyint
sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy( 2) = onesixth*yintsq*yint
end if

if (noz==1) then
sz( 0) = 1.0_num-zint
sz( 1) = zint
elseif (noz==2) then
zintsq = zint*zint
sz(-1) = 0.5_num*(0.5_num-zint)**2
sz( 0) = 0.75_num-zintsq
sz( 1) = 0.5_num*(0.5_num+zint)**2
elseif (noz==3) then
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz(-1) = onesixth*ozintsq*ozint
sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz( 2) = onesixth*zintsq*zint
end if

xint=x-0.5_num-j0
yint=y-0.5_num-k0
zint=z-0.5_num-l0

if (l_lower_order_in_v) then

if (nox==1) then
sx0( 0) = 1.0_num
elseif (nox==2) then
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
elseif (nox==3) then
xintsq = xint*xint
sx0(-1) = 0.5_num*(0.5_num-xint)**2
sx0( 0) = 0.75_num-xintsq
sx0( 1) = 0.5_num*(0.5_num+xint)**2
end if

if (noy==1) then
sy0( 0) = 1.0_num
elseif (noy==2) then
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
elseif (noy==3) then
yintsq = yint*yint
sy0(-1) = 0.5_num*(0.5_num-yint)**2
sy0( 0) = 0.75_num-yintsq
sy0( 1) = 0.5_num*(0.5_num+yint)**2
end if

if (noz==1) then
sz0( 0) = 1.0_num
elseif (noz==2) then
sz0( 0) = 1.0_num-zint
sz0( 1) = zint
elseif (noz==3) then
zintsq = zint*zint
sz0(-1) = 0.5_num*(0.5_num-zint)**2
sz0( 0) = 0.75_num-zintsq
sz0( 1) = 0.5_num*(0.5_num+zint)**2
end if

else

if (nox==1) then
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
elseif (nox==2) then
xintsq = xint*xint
sx0(-1) = 0.5_num*(0.5_num-xint)**2
sx0( 0) = 0.75_num-xintsq
sx0( 1) = 0.5_num*(0.5_num+xint)**2
elseif (nox==3) then
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx0(-1) = onesixth*oxintsq*oxint
sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx0( 2) = onesixth*xintsq*xint
end if

if (noy==1) then
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
elseif (noy==2) then
yintsq = yint*yint
sy0(-1) = 0.5_num*(0.5_num-yint)**2
sy0( 0) = 0.75_num-yintsq
sy0( 1) = 0.5_num*(0.5_num+yint)**2
elseif (noy==3) then
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy0(-1) = onesixth*oyintsq*oyint
sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy0( 2) = onesixth*yintsq*yint
end if

if (noz==1) then
sz0( 0) = 1.0_num-zint
sz0( 1) = zint
elseif (noz==2) then
zintsq = zint*zint
sz0(-1) = 0.5_num*(0.5_num-zint)**2
sz0( 0) = 0.75_num-zintsq
sz0( 1) = 0.5_num*(0.5_num+zint)**2
elseif (noz==3) then
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz0(-1) = onesixth*ozintsq*ozint
sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz0( 2) = onesixth*zintsq*zint
end if

end if

do ll = izmin, izmax+1
do kk = iymin, iymax+1
do jj = ixmin0, ixmax0
ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)*signx
end do
end do
end do

do ll = izmin, izmax+1
do kk = iymin0, iymax0
do jj = ixmin, ixmax+1
ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)*signy
end do
end do
end do

do ll = izmin0, izmax0
do kk = iymin, iymax+1
do jj = ixmin, ixmax+1
ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj,k+kk,l0+ll)
end do
end do
end do

end do
!$omp end parallel do
deallocate(sx0,sy0,sz0)

return
end subroutine gete3d_n_energy_conserving
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz,bxg,byg,bzg,l_lower_order_in_v)

use omp_lib
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
real(kind=8), dimension(np) :: xp,yp,zp,bx,by,bz
logical :: l_lower_order_in_v
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
integer(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax,ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint,xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
real(kind=8), dimension(-int(nox/2):int((nox+1)/2)) :: sx
real(kind=8), dimension(-int(noy/2):int((noy+1)/2)) :: sy
real(kind=8), dimension(-int(noz/2):int((noz+1)/2)) :: sz
real(kind=8), dimension(:), pointer :: sx0,sy0,sz0
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
iymin = -int(noy/2)
iymax =  int((noy+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1


if (l_lower_order_in_v) then
ixmin0 = -int((nox-1)/2)
ixmax0 =  int((nox)/2)
iymin0 = -int((noy-1)/2)
iymax0 =  int((noy)/2)
izmin0 = -int((noz-1)/2)
izmax0 =  int((noz)/2)
else
ixmin0 = -int((nox)/2)
ixmax0 =  int((nox+1)/2)
iymin0 = -int((noy)/2)
iymax0 =  int((noy+1)/2)
izmin0 = -int((noz)/2)
izmax0 =  int((noz+1)/2)
end if
allocate(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

signx = 1.0_num
signy = 1.0_num

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!$omp parallel do private(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!$omp sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
do ip=1,np
x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi

if (l_lower_order_in_v) then
if (nox==2*(nox/2)) then
j=nint(x)
j0=floor(x-0.5_num)
else
j=floor(x)
j0=floor(x)
end if
if (noy==2*(noy/2)) then
k=nint(y)
k0=floor(y-0.5_num)
else
k=floor(y)
k0=floor(y)
end if
if (noz==2*(noz/2)) then
l=nint(z)
l0=floor(z-0.5_num)
else
l=floor(z)
l0=floor(z)
end if
else
if (nox==2*(nox/2)) then
j=nint(x)
j0=floor(x)
else
j=floor(x)
j0=floor(x-0.5_num)
end if
if (noy==2*(noy/2)) then
k=nint(y)
k0=floor(y)
else
k=floor(y)
k0=floor(y-0.5_num)
end if
if (noz==2*(noz/2)) then
l=nint(z)
l0=floor(z)
else
l=floor(z)
l0=floor(z-0.5_num)
end if
end if

xint=x-j
yint=y-k
zint=z-l

if (nox==1) then
sx( 0) = 1.0_num-xint
sx( 1) = xint
elseif (nox==2) then
xintsq = xint*xint
sx(-1) = 0.5_num*(0.5_num-xint)**2
sx( 0) = 0.75_num-xintsq
sx( 1) = 0.5_num*(0.5_num+xint)**2
elseif (nox==3) then
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx(-1) = onesixth*oxintsq*oxint
sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx( 2) = onesixth*xintsq*xint
end if

if (noy==1) then
sy( 0) = 1.0_num-yint
sy( 1) = yint
elseif (noy==2) then
yintsq = yint*yint
sy(-1) = 0.5_num*(0.5_num-yint)**2
sy( 0) = 0.75_num-yintsq
sy( 1) = 0.5_num*(0.5_num+yint)**2
elseif (noy==3) then
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy(-1) = onesixth*oyintsq*oyint
sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy( 2) = onesixth*yintsq*yint
end if

if (noz==1) then
sz( 0) = 1.0_num-zint
sz( 1) = zint
elseif (noz==2) then
zintsq = zint*zint
sz(-1) = 0.5_num*(0.5_num-zint)**2
sz( 0) = 0.75_num-zintsq
sz( 1) = 0.5_num*(0.5_num+zint)**2
elseif (noz==3) then
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz(-1) = onesixth*ozintsq*ozint
sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz( 2) = onesixth*zintsq*zint
end if

xint=x-0.5_num-j0
yint=y-0.5_num-k0
zint=z-0.5_num-l0

if (l_lower_order_in_v) then
if (nox==1) then
sx0( 0) = 1.0_num
elseif (nox==2) then
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
elseif (nox==3) then
xintsq = xint*xint
sx0(-1) = 0.5_num*(0.5_num-xint)**2
sx0( 0) = 0.75_num-xintsq
sx0( 1) = 0.5_num*(0.5_num+xint)**2
end if

if (noy==1) then
sy0( 0) = 1.0_num
elseif (noy==2) then
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
elseif (noy==3) then
yintsq = yint*yint
sy0(-1) = 0.5_num*(0.5_num-yint)**2
sy0( 0) = 0.75_num-yintsq
sy0( 1) = 0.5_num*(0.5_num+yint)**2
end if

if (noz==1) then
sz0( 0) = 1.0_num
elseif (noz==2) then
sz0( 0) = 1.0_num-zint
sz0( 1) = zint
elseif (noz==3) then
zintsq = zint*zint
sz0(-1) = 0.5_num*(0.5_num-zint)**2
sz0( 0) = 0.75_num-zintsq
sz0( 1) = 0.5_num*(0.5_num+zint)**2
end if
else

if (nox==1) then
sx0( 0) = 1.0_num-xint
sx0( 1) = xint
elseif (nox==2) then
xintsq = xint*xint
sx0(-1) = 0.5_num*(0.5_num-xint)**2
sx0( 0) = 0.75_num-xintsq
sx0( 1) = 0.5_num*(0.5_num+xint)**2
elseif (nox==3) then
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx0(-1) = onesixth*oxintsq*oxint
sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx0( 2) = onesixth*xintsq*xint
end if

if (noy==1) then
sy0( 0) = 1.0_num-yint
sy0( 1) = yint
elseif (noy==2) then
yintsq = yint*yint
sy0(-1) = 0.5_num*(0.5_num-yint)**2
sy0( 0) = 0.75_num-yintsq
sy0( 1) = 0.5_num*(0.5_num+yint)**2
elseif (noy==3) then
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy0(-1) = onesixth*oyintsq*oyint
sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy0( 2) = onesixth*yintsq*yint
end if

if (noz==1) then
sz0( 0) = 1.0_num-zint
sz0( 1) = zint
elseif (noz==2) then
zintsq = zint*zint
sz0(-1) = 0.5_num*(0.5_num-zint)**2
sz0( 0) = 0.75_num-zintsq
sz0( 1) = 0.5_num*(0.5_num+zint)**2
elseif (noz==3) then
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz0(-1) = onesixth*ozintsq*ozint
sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz0( 2) = onesixth*zintsq*zint
end if
end if

do ll = izmin0, izmax0
do kk = iymin0, iymax0
do jj = ixmin, ixmax+1
bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)*signx
end do
end do
end do

do ll = izmin0, izmax0
do kk = iymin, iymax+1
do jj = ixmin0, ixmax0
by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)*signy
end do
end do
end do

do ll = izmin, izmax+1
do kk = iymin0, iymax0
do jj = ixmin0, ixmax0
bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj,k0+kk,l+ll)
end do
end do
end do
end do

deallocate(sx0,sz0)

return
end subroutine getb3d_n_energy_conserving
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine field_bc(field, nxg, nyg, nzg, nx_local, ny_local, nz_local)
use boundary

integer(idp), intent(in) :: nxg, nyg, nzg
integer(idp), intent(in) :: nx_local, ny_local, nz_local
real(kind=8), dimension(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), intent(inout) :: field

call exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg, nx, ny, nz)

end subroutine field_bc
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg,nx_local, ny_local, nz_local)
use boundary

integer(idp), intent(in) :: nxg, nyg, nzg
integer(idp), intent(in) :: nx_local, ny_local, nz_local
real(kind=8), dimension(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), intent(inout) :: field
integer(idp), dimension(c_ndims) :: sizes, subsizes, starts
integer(isp) :: subarray, basetype, sz, szmax, i, j, k, n
real(kind=8), pointer :: temp(:)

basetype = mpidbl

sizes(1) = nx_local + 1 + 2 * nxg
sizes(2) = ny_local + 1 + 2 * nyg
sizes(3) = nz_local + 1 + 2 * nzg
starts = 1

szmax = sizes(1) * sizes(2) * nzg
sz = sizes(1) * sizes(3) * nyg
if (sz .gt. szmax) szmax = sz
sz = sizes(2) * sizes(3) * nxg
if (sz .gt. szmax) szmax = sz

allocate(temp(szmax))

subsizes(1) = nxg
subsizes(2) = sizes(2)
subsizes(3) = sizes(3)

sz = subsizes(1) * subsizes(2) * subsizes(3)

subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)


call mpi_sendrecv(field(0,-nyg,-nzg), 1_isp, subarray, proc_x_min,tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

if (proc_x_max .ne. mpi_proc_null) then
n = 1
do k = -nzg, subsizes(3)-nzg-1
do j = -nyg, subsizes(2)-nyg-1
do i = nx_local+1, subsizes(1)+nx_local
field(i,j,k) = temp(n)
n = n + 1
enddo
enddo
enddo
endif

call mpi_sendrecv(field(nx_local+1-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_max,tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

if (proc_x_min .ne. mpi_proc_null) then
n = 1
do k = -nzg, subsizes(3)-nzg-1
do j = -nyg, subsizes(2)-nyg-1
do i = -nxg, subsizes(1)-nxg-1
field(i,j,k) = temp(n)
n = n + 1
enddo
enddo
enddo
endif

call mpi_type_free(subarray, errcode)

subsizes(1) = sizes(1)
subsizes(2) = nyg
subsizes(3) = sizes(3)

sz = subsizes(1) * subsizes(2) * subsizes(3)

subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)


call mpi_sendrecv(field(-nxg,0,-nzg), 1_isp, subarray, proc_y_min,tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

if (proc_y_max .ne. mpi_proc_null) then
n = 1
do k = -nzg, subsizes(3)-nzg-1
do j = ny_local+1, subsizes(2)+ny_local
do i = -nxg, subsizes(1)-nxg-1
field(i,j,k) = temp(n)
n = n + 1
enddo
enddo
enddo
endif

call mpi_sendrecv(field(-nxg,ny_local+1-nyg,-nzg), 1_isp, subarray, proc_y_max,tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

if (proc_y_min .ne. mpi_proc_null) then
n = 1
do k = -nzg, subsizes(3)-nzg-1
do j = -nyg, subsizes(2)-nyg-1
do i = -nxg, subsizes(1)-nxg-1
field(i,j,k) = temp(n)
n = n + 1
enddo
enddo
enddo
endif

call mpi_type_free(subarray, errcode)

subsizes(1) = sizes(1)
subsizes(2) = sizes(2)
subsizes(3) = nzg

sz = subsizes(1) * subsizes(2) * subsizes(3)

subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)


call mpi_sendrecv(field(-nxg,-nyg,0), 1_isp, subarray, proc_z_min,tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

if (proc_z_max .ne. mpi_proc_null) then
n = 1
do k = nz_local+1, subsizes(3)+nz_local
do j = -nyg, subsizes(2)-nyg-1
do i = -nxg, subsizes(1)-nxg-1
field(i,j,k) = temp(n)
n = n + 1
enddo
enddo
enddo
endif

call mpi_sendrecv(field(-nxg,-nyg,nz_local+1-nzg), 1_isp, subarray, proc_z_max,tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

if (proc_z_min .ne. mpi_proc_null) then
n = 1
do k = -nzg, subsizes(3)-nzg-1
do j = -nyg, subsizes(2)-nyg-1
do i = -nxg, subsizes(1)-nxg-1
field(i,j,k) = temp(n)
n = n + 1
enddo
enddo
enddo
endif

call mpi_type_free(subarray, errcode)

deallocate(temp)

end subroutine exchange_mpi_3d_grid_array_with_guards
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine exchange_mpi_3d_grid_array_with_guards_nonblocking(field, nxg, nyg, nzg,nx_local, ny_local, nz_local)
use boundary

integer(idp), intent(in) :: nxg, nyg, nzg
integer(idp), intent(in) :: nx_local, ny_local, nz_local
real(kind=8), dimension(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), intent(inout) :: field
integer(idp), dimension(c_ndims) :: sizes, subsizes, starts
integer(isp) :: subarray, basetype, sz, szmax, i, j, k, n
real(kind=8), pointer :: temp(:)
integer(isp):: requests(4)

basetype = mpidbl

sizes(1) = nx_local + 1 + 2 * nxg
sizes(2) = ny_local + 1 + 2 * nyg
sizes(3) = nz_local + 1 + 2 * nzg
starts = 1


subsizes(1) = nxg
subsizes(2) = sizes(2)
subsizes(3) = sizes(3)

subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

call mpi_isend(field(0,-nyg,-nzg), 1_isp, subarray, proc_x_min, tag,comm, requests(1), errcode)
call mpi_irecv(field(nx_local+1,-nyg,-nzg), 1_isp, subarray, proc_x_max, tag,comm, requests(2), errcode)

call mpi_isend(field(nx_local+1-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_max, tag,comm, requests(3), errcode)
call mpi_irecv(field(-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_min, tag,comm, requests(4), errcode)

call mpi_type_free(subarray, errcode)


call mpi_waitall(4_isp, requests, mpi_statuses_ignore, errcode)


subsizes(1) = sizes(1)
subsizes(2) = nyg
subsizes(3) = sizes(3)

subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

call mpi_isend(field(-nxg,0,-nzg), 1_isp, subarray, proc_y_min, tag,comm, requests(1), errcode)
call mpi_irecv(field(-nxg,ny_local+1,-nzg), 1_isp, subarray, proc_y_max, tag,comm, requests(2), errcode)

call mpi_isend(field(-nxg,ny_local+1-nyg,-nzg), 1_isp, subarray, proc_y_max, tag,comm, requests(3), errcode)
call mpi_irecv(field(-nxg,-nyg,-nzg), 1_isp, subarray, proc_y_min, tag,comm, requests(4), errcode)


call mpi_waitall(4_isp, requests, mpi_statuses_ignore, errcode)


call mpi_type_free(subarray, errcode)


subsizes(1) = sizes(1)
subsizes(2) = sizes(2)
subsizes(3) = nzg

sz = subsizes(1) * subsizes(2) * subsizes(3)

subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

call mpi_isend(field(-nxg,-nyg,0), 1_isp, subarray, proc_z_min, tag,comm, requests(1), errcode)
call mpi_irecv(field(-nxg,-nyg,nz_local+1), 1_isp, subarray, proc_z_max, tag,comm, requests(2), errcode)

call mpi_isend(field(-nxg,-nyg,nz_local+1-nzg), 1_isp, subarray, proc_z_max, tag,comm, requests(3), errcode)
call mpi_irecv(field(-nxg,-nyg,-nzg), 1_isp, subarray, proc_z_min, tag,comm, requests(4), errcode)

call mpi_waitall(4_isp, requests, mpi_statuses_ignore, errcode)
call mpi_type_free(subarray, errcode)

end subroutine exchange_mpi_3d_grid_array_with_guards_nonblocking
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine summation_bcs(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
use boundary

integer(idp), intent(in) :: nxg, nyg, nzg
integer(idp), intent(in) :: nx_local, ny_local, nz_local
real(kind=8), dimension(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), intent(inout) :: array
real(kind=8), dimension(:,:,:), pointer :: temp
integer(idp), dimension(c_ndims) :: sizes, subsizes, starts
integer(isp) :: subarray, nn, sz, i

sizes(1) = nx + 1 + 2 * nxg
sizes(2) = ny + 1 + 2 * nyg
sizes(3) = nz + 1 + 2 * nzg
starts = 1


subsizes(1) = nxg
subsizes(2) = sizes(2)
subsizes(3) = sizes(3)
nn = nx

subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

sz = subsizes(1) * subsizes(2) * subsizes(3)
allocate(temp(subsizes(1), subsizes(2), subsizes(3)))

temp = 0.0_num
call mpi_sendrecv(array(nn+1,-nyg,-nzg), 1_isp, subarray,neighbour( 1,0,0), tag, temp, sz, mpidbl,neighbour(-1,0,0), tag, comm, status, errcode)
array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp

temp = 0.0_num
call mpi_sendrecv(array(-nxg,-nyg,-nzg), 1_isp, subarray,neighbour(-1,0,0), tag, temp, sz, mpidbl,neighbour( 1,0,0), tag, comm, status, errcode)
array(nn+1-nxg:nn,:,:) = array(nn+1-nxg:nn,:,:) + temp

deallocate(temp)
call mpi_type_free(subarray, errcode)


subsizes(1) = sizes(1)
subsizes(2) = nyg
subsizes(3) = sizes(3)
nn = ny

subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

sz = subsizes(1) * subsizes(2) * subsizes(3)
allocate(temp(subsizes(1), subsizes(2), subsizes(3)))

temp = 0.0_num
call mpi_sendrecv(array(-nxg,nn+1,-nzg), 1_isp, subarray,neighbour(0, 1,0), tag, temp, sz, mpidbl,neighbour(0,-1,0), tag, comm, status, errcode)
array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp

temp = 0.0_num
call mpi_sendrecv(array(-nxg,-nyg,-nzg), 1_isp, subarray,neighbour(0,-1,0), tag, temp, sz, mpidbl,neighbour(0, 1,0), tag, comm, status, errcode)
array(:,nn+1-nyg:nn,:) = array(:,nn+1-nyg:nn,:) + temp

deallocate(temp)
call mpi_type_free(subarray, errcode)


subsizes(1) = sizes(1)
subsizes(2) = sizes(2)
subsizes(3) = nzg
nn = nz

subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

sz = subsizes(1) * subsizes(2) * subsizes(3)
allocate(temp(subsizes(1), subsizes(2), subsizes(3)))

temp = 0.0_num
call mpi_sendrecv(array(-nxg,-nyg,nn+1), 1_isp, subarray,neighbour(0,0, 1), tag, temp, sz, mpidbl,neighbour(0,0,-1), tag, comm, status, errcode)
array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp

temp = 0.0_num
call mpi_sendrecv(array(-nxg,-nyg,-nzg), 1_isp, subarray,neighbour(0,0,-1), tag, temp, sz, mpidbl,neighbour(0,0, 1), tag, comm, status, errcode)
array(:,:,nn+1-nzg:nn) = array(:,:,nn+1-nzg:nn) + temp

deallocate(temp)
call mpi_type_free(subarray, errcode)

call field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

end subroutine summation_bcs
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine summation_bcs_nonblocking(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)
use boundary

integer(idp), intent(in) :: nxg, nyg, nzg
integer(idp), intent(in) :: nx_local, ny_local, nz_local
real(kind=8), dimension(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), intent(inout) :: array
real(kind=8), dimension(:,:,:), pointer :: temp1, temp2
integer(idp), dimension(c_ndims) :: sizes, subsizes, starts
integer(isp) :: subarray, nn, sz, i
integer(isp) :: requests(4)

sizes(1) = nx + 1 + 2 * nxg
sizes(2) = ny + 1 + 2 * nyg
sizes(3) = nz + 1 + 2 * nzg
starts = 1


subsizes(1) = nxg
subsizes(2) = sizes(2)
subsizes(3) = sizes(3)
nn = nx

subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

sz = subsizes(1) * subsizes(2) * subsizes(3)
allocate(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

temp1  = 0.0_num
temp2 = 0.0_num
call mpi_isend(array(nn+1,-nyg,-nzg), 1_isp, subarray, proc_x_max, tag,comm, requests(1), errcode)
call mpi_irecv(temp1, sz, mpidbl, proc_x_min, tag,comm, requests(2), errcode)
call mpi_isend(array(-nxg,-nyg,-nzg), 1_isp, subarray, proc_x_min, tag,comm, requests(3), errcode)
call mpi_irecv(temp2, sz, mpidbl, proc_x_max, tag,comm, requests(4), errcode)
call mpi_waitall(4_isp, requests, mpi_statuses_ignore, errcode)

array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp1
array(nn+1-nxg:nn,:,:) = array(nn+1-nxg:nn,:,:) + temp2

deallocate(temp1,temp2)
call mpi_type_free(subarray, errcode)


subsizes(1) = sizes(1)
subsizes(2) = nyg
subsizes(3) = sizes(3)
nn = ny

subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

sz = subsizes(1) * subsizes(2) * subsizes(3)
allocate(temp1(subsizes(1), subsizes(2), subsizes(3)), temp2(subsizes(1), subsizes(2), subsizes(3)))

temp1  = 0.0_num
temp2 = 0.0_num
call mpi_isend(array(-nxg,nn+1,-nzg), 1_isp, subarray, proc_y_max, tag,comm, requests(1), errcode)
call mpi_irecv(temp1, sz, mpidbl, proc_y_min, tag,comm, requests(2), errcode)
call mpi_isend(array(-nxg,-nyg,-nzg), 1_isp, subarray, proc_y_min, tag,comm, requests(3), errcode)
call mpi_irecv(temp2, sz, mpidbl, proc_y_max, tag,comm, requests(4), errcode)
call mpi_waitall(4_isp, requests, mpi_statuses_ignore, errcode)

array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp1
array(:,nn+1-nyg:nn,:) = array(:,nn+1-nyg:nn,:) + temp2

deallocate(temp1,temp2)
call mpi_type_free(subarray, errcode)


subsizes(1) = sizes(1)
subsizes(2) = sizes(2)
subsizes(3) = nzg
nn = nz

subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

sz = subsizes(1) * subsizes(2) * subsizes(3)
allocate(temp1(subsizes(1), subsizes(2), subsizes(3)),temp2(subsizes(1), subsizes(2), subsizes(3)))

temp1  = 0.0_num
temp2 = 0.0_num
call mpi_isend(array(-nxg,-nyg,nn+1), 1_isp, subarray, proc_z_max, tag,comm, requests(1), errcode)
call mpi_irecv(temp1, sz, mpidbl, proc_z_min, tag,comm, requests(2), errcode)
call mpi_isend(array(-nxg,-nyg,-nzg), 1_isp, subarray, proc_z_min, tag,comm, requests(3), errcode)
call mpi_irecv(temp2, sz, mpidbl, proc_z_max, tag,comm, requests(4), errcode)
call mpi_waitall(4_isp, requests, mpi_statuses_ignore, errcode)

array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp1
array(:,:,nn+1-nzg:nn) = array(:,:,nn+1-nzg:nn) + temp2

deallocate(temp1,temp2)
call mpi_type_free(subarray, errcode)

call field_bc(array, nxg, nyg, nzg, nx_local, ny_local, nz_local)

end subroutine summation_bcs_nonblocking
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine efield_bcs
use boundary

call field_bc(ex, nxguards, nyguards, nzguards, nx, ny, nz)
call field_bc(ey, nxguards, nyguards, nzguards, nx, ny, nz)
call field_bc(ez, nxguards, nyguards, nzguards, nx, ny, nz)
end subroutine efield_bcs
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine bfield_bcs
use boundary

call field_bc(bx, nxguards, nyguards, nzguards, nx, ny, nz)
call field_bc(by, nxguards, nyguards, nzguards, nx, ny, nz)
call field_bc(bz, nxguards, nyguards, nzguards, nx, ny, nz)

end subroutine bfield_bcs
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine current_bcs
use boundary

call summation_bcs_nonblocking(jx, nxjguards, nyjguards, nzjguards, nx, ny, nz)
call summation_bcs_nonblocking(jy, nxjguards, nyjguards, nzjguards, nx, ny, nz)
call summation_bcs_nonblocking(jz, nxjguards, nyjguards, nzjguards, nx, ny, nz)
end subroutine current_bcs
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine charge_bcs
use boundary

call summation_bcs_nonblocking(rho, nxjguards, nyjguards, nzjguards, nx, ny, nz)
end subroutine charge_bcs
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine particle_bcs
use boundary
implicit none


call particle_bcs_tiles


call particle_bcs_mpi_blocking

end subroutine particle_bcs
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine particle_bcs_tiles
use interf_add_particle_at_tile
use interf_rm_particle_at_tile
use boundary
implicit none
integer:: i, ispecies, ix, iy, iz, indx, indy, indz
integer(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
type(particle_species), pointer :: curr
type(particle_tile), pointer :: curr_tile, curr_tile_add
real(kind=8) :: partx, party, partz, partux, partuy, partuz, partw
integer(idp) :: test =0

do ispecies=1, nspecies
curr=> species_parray(ispecies)

nx0_grid_tile = curr%array_of_tiles(1,1,1)%nx_grid_tile
ny0_grid_tile = curr%array_of_tiles(1,1,1)%ny_grid_tile
nz0_grid_tile = curr%array_of_tiles(1,1,1)%nz_grid_tile
do iz=1, ntilez
do iy=1, ntiley
do ix=1, ntilex
curr_tile=>curr%array_of_tiles(ix,iy,iz)
nptile=curr_tile%np_tile
do i=nptile, 1, -1
partx=curr_tile%part_x(i)
party=curr_tile%part_y(i)
partz=curr_tile%part_z(i)
partux=curr_tile%part_ux(i)
partuy=curr_tile%part_uy(i)
partuz=curr_tile%part_uz(i)
partw=curr_tile%weight(i)


if (((partx .ge. curr_tile%x_tile_min) .and. (partx .lt. curr_tile%x_tile_max)).and. ((party .ge. curr_tile%y_tile_min) .and. (party .lt. curr_tile%y_tile_max)).and. ((partz .ge. curr_tile%z_tile_min) .and. (partz .lt. curr_tile%z_tile_max)))cycle


if ((partx .lt. x_min_local) .or. (partx .ge. x_max_local)) cycle
if ((party .lt. y_min_local) .or. (party .ge. y_max_local)) cycle
if ((partz .lt. z_min_local) .or. (partz .ge. z_max_local)) cycle



indx = min(floor((partx-x_min_local)/(nx0_grid_tile*dx))+1,ntilex)
indy = min(floor((party-y_min_local)/(ny0_grid_tile*dy))+1,ntiley)
indz = min(floor((partz-z_min_local)/(nz0_grid_tile*dz))+1,ntilez)
call rm_particle_at_tile(curr_tile,i)
curr_tile_add=>curr%array_of_tiles(indx,indy,indz)
call add_particle_at_tile(curr_tile_add,partx, party, partz, partux, partuy, partuz, partw)
end do
end do
end do
end do
end do
end subroutine particle_bcs_tiles
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine particle_bcs_mpi_blocking
use interf_add_particle_to_species
use interf_rm_particles_from_species
use boundary
integer(isp), parameter :: nvar=7
integer(isp), dimension(-1:1,-1:1,-1:1) :: nptoexch
real(kind=8), pointer, dimension(:,:,:,:) :: sendbuf
real(kind=8), pointer, dimension(:) :: recvbuf
real(kind=8), pointer, dimension(:) :: temp
logical, pointer, dimension(:) :: mask
integer(isp) :: ibuff, isend, nout, nbuff, ninit
integer(isp) :: xbd, ybd, zbd
integer(isp) :: ixp, iyp, izp
integer(isp) :: nsend_buf, nrecv_buf, npart_curr
integer(isp) :: dest, src
logical :: out_of_bounds
integer(isp) :: ispecies, i, ip, ix, iy, iz
integer(isp) :: ixtile, iytile, iztile
real(kind=8) :: part_xyz
type(particle_species), pointer :: currsp
type(particle_tile), pointer :: curr

do ispecies=1, nspecies

currsp => species_parray(ispecies)
nptoexch=0
nsend_buf=0
nout=0
nrecv_buf=0
nbuff=currsp%species_npart*nvar
ibuff=1
allocate(sendbuf(-1:1,-1:1,-1:1,1:nbuff))
do iztile=1, ntilez
do iytile=1, ntiley
do ixtile=1, ntilex
curr=>currsp%array_of_tiles(ixtile,iytile,iztile)

if (.not. curr%subdomain_bound) cycle

allocate(mask(1:curr%np_tile))
mask=.true.
xbd = 0
ybd = 0
zbd = 0
part_xyz=0.

do i = 1, curr%np_tile
xbd = 0
ybd = 0
zbd = 0
out_of_bounds = .false.
part_xyz = curr%part_x(i)

if (part_xyz .lt. x_min_local) then
xbd = -1
if (x_min_boundary) then
curr%part_x(i) = part_xyz + length_x
endif
endif

if (part_xyz .ge. x_max_local) then
xbd = 1
if (x_max_boundary) then
curr%part_x(i) = part_xyz - length_x
endif
endif

part_xyz = curr%part_y(i)

if (part_xyz .lt. y_min_local) then
ybd = -1
if (y_min_boundary) then
curr%part_y(i) = part_xyz + length_y
endif
endif


if (part_xyz .ge. y_max_local) then
ybd = 1
if (y_max_boundary) then
curr%part_y(i) = part_xyz - length_y
endif
endif

part_xyz = curr%part_z(i)

if (part_xyz .lt. z_min_local) then
zbd = -1
if (z_min_boundary) then
curr%part_z(i) = part_xyz + length_z
endif
endif


if (part_xyz .ge. z_max_local) then
zbd = 1

if (z_max_boundary) then
curr%part_z(i) = part_xyz - length_z
endif
endif

if (abs(xbd) + abs(ybd) + abs(zbd) .gt. 0) then

mask(i)=.false.
nout=nout+1
ibuff=nptoexch(xbd,ybd,zbd)*nvar+1
sendbuf(xbd,ybd,zbd,ibuff)    = curr%part_x(i)
sendbuf(xbd,ybd,zbd,ibuff+1)  = curr%part_y(i)
sendbuf(xbd,ybd,zbd,ibuff+2)  = curr%part_z(i)
sendbuf(xbd,ybd,zbd,ibuff+3)  = curr%part_ux(i)
sendbuf(xbd,ybd,zbd,ibuff+4)  = curr%part_uy(i)
sendbuf(xbd,ybd,zbd,ibuff+5)  = curr%part_uz(i)
sendbuf(xbd,ybd,zbd,ibuff+6)  = curr%weight(i)
nptoexch(xbd,ybd,zbd) = nptoexch(xbd,ybd,zbd)+1
endif
enddo

call rm_particles_from_species(currsp, curr, mask)
deallocate(mask)
enddo
enddo
enddo

do iz = -1, 1
do iy = -1, 1
do ix = -1, 1
if (abs(ix) + abs(iy) + abs(iz) .eq. 0) cycle
ixp = -ix
iyp = -iy
izp = -iz


nsend_buf=nptoexch(ix,iy,iz)*nvar
nrecv_buf=0
dest = neighbour(ix,iy,iz)
src  = neighbour(ixp,iyp,izp)
call mpi_sendrecv(nsend_buf, 1_isp, mpi_integer, dest, tag, nrecv_buf, 1_isp,mpi_integer, src, tag, comm, status, errcode)
allocate(recvbuf(1:nrecv_buf))
call mpi_sendrecv(sendbuf(ix,iy,iz,1:nsend_buf), nsend_buf, mpidbl, dest, tag,recvbuf, nrecv_buf, mpidbl, src, tag, comm, status, errcode)

do i =1, nrecv_buf, nvar
call add_particle_to_species(currsp, recvbuf(i), recvbuf(i+1), recvbuf(i+2),recvbuf(i+3), recvbuf(i+4), recvbuf(i+5), recvbuf(i+6))
end do
deallocate(recvbuf)
enddo
enddo
enddo
deallocate(sendbuf)
end do
end subroutine particle_bcs_mpi_blocking
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine output_routines()
use simple_io
use shared_data
use params
implicit none
character(len=string_length) :: strtemp
integer(kind=mpi_offset_kind) :: offset=0
integer(kind=4) :: err=0
if (output_frequency .lt. 1) return
if ((it .ge. output_step_min) .and. (it .le. output_step_max) .and.(mod(it-output_step_min,output_frequency) .eq. 0)) then

write(strtemp,'(i5)') it

if (c_output_ex .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(fileex))//trim(adjustl(strtemp))//'.pxr', ex, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_ey .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(fileey))//trim(adjustl(strtemp))//'.pxr', ey, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_ez .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(fileez))//trim(adjustl(strtemp))//'.pxr', ez, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_bx .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filebx))//trim(adjustl(strtemp))//'.pxr', bx, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_by .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(fileby))//trim(adjustl(strtemp))//'.pxr', by, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_bz .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filebz))//trim(adjustl(strtemp))//'.pxr', bz, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_jx .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filejx))//trim(adjustl(strtemp))//'.pxr', jx, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
endif
if (c_output_jy .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filejy))//trim(adjustl(strtemp))//'.pxr', jy, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
endif
if (c_output_jz .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filejz))//trim(adjustl(strtemp))//'.pxr', jz, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
endif
if (c_output_dive .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filedive))//trim(adjustl(strtemp))//'.pxr', dive, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
endif
if (c_output_rho .eq. 1) then

call write_single_array_to_file('./results/'//trim(adjustl(filerho))//trim(adjustl(strtemp))//'.pxr', rho, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
endif
endif
end subroutine output_routines
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine write_single_array_to_file(filename, array, nxg, nyg, nzg, nx_local, ny_local, nz_local, offset, err)
use simple_io

character(len=*), intent(in) :: filename
integer(idp), intent(in) :: nxg, nyg, nzg
integer(idp), intent(in) :: nx_local, ny_local, nz_local
real(kind=8), dimension(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), intent(inout) :: array
integer(kind=mpi_offset_kind), intent(in) :: offset
integer(isp), intent(inout) :: err
integer(isp) :: subt, suba, fh, i

call mpi_file_open(comm, trim(filename), mpi_mode_create + mpi_mode_wronly,mpi_info_null, fh, errcode)

if (errcode .ne. 0) then
if (rank .eq. 0) print *, 'file ', trim(filename), ' could not be created - check disk space'
err = ior(err, c_err_bad_value)
return
endif

subt = create_current_grid_derived_type()
suba = create_current_grid_subarray(nxguards, nyguards, nzguards)
call mpi_file_set_view(fh, offset, mpi_byte, subt, 'native',mpi_info_null, errcode)

call mpi_file_write_all(fh, array, 1_isp, suba, mpi_status_ignore, errcode)

call mpi_file_close(fh, errcode)
call mpi_type_free(subt, errcode)

end subroutine write_single_array_to_file
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine calc_diags
use diagnostics
use fields
use boundary
use particles
use params
use shared_data
use tiling
implicit none
integer(idp) :: ispecies, ix, iy, iz, count
integer(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
integer(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
type(particle_species), pointer :: curr
type(particle_tile), pointer :: curr_tile
integer(idp) :: nxc, nyc, nzc


dive=0.0_num
call calc_field_div(dive, ex, ey, ez, nx, ny, nz, nxguards, nyguards, nzguards, dx, dy, dz)


rho=0.0_num
!$omp parallel default(none) &
!$omp shared(rho,ntilex,ntiley,ntilez,nspecies,species_parray,nxjguards, &
!$omp nyjguards,nzjguards,dx,dy,dz) &
!$omp private(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin,lmax,nxc,nyc,nzc, &
!$omp jminc,jmaxc,kminc,kmaxc,lminc,lmaxc)
!$omp do collapse(3) schedule(runtime)
do iz=1, ntilez
do iy=1, ntiley
do ix=1, ntilex
do ispecies=1, nspecies
curr => species_parray(ispecies)
curr_tile=>curr%array_of_tiles(ix,iy,iz)
count= curr_tile%np_tile
jmin=curr_tile%nx_tile_min
jmax=curr_tile%nx_tile_max
kmin=curr_tile%ny_tile_min
kmax=curr_tile%ny_tile_max
lmin=curr_tile%nz_tile_min
lmax=curr_tile%nz_tile_max
nxc=curr_tile%nx_cells_tile
nyc=curr_tile%ny_cells_tile
nzc=curr_tile%nz_cells_tile
curr_tile%rhotile = 0.0_num

call depose_rho_vechv_1_1_1(curr_tile%rhotile, count,curr_tile%part_x(1:count),curr_tile%part_y(1:count),curr_tile%part_z(1:count),curr_tile%weight(1:count), curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min,dx,dy,dz,curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,nzjguards)

rho(jmin:jmax,kmin:kmax,lmin:lmax) = rho(jmin:jmax,kmin:kmax,lmin:lmax)+curr_tile%rhotile(0:nxc,0:nyc,0:nzc)
end do
end do
end do
end do
!$omp end do

!$omp do collapse(3) schedule(runtime)
do iz=1,ntilez
do iy=1,ntiley
do ix=1,ntilex
do ispecies=1, nspecies
curr => species_parray(ispecies)
curr_tile=>curr%array_of_tiles(ix,iy,iz)
count=curr_tile%np_tile
jmin=curr_tile%nx_tile_min
jmax=curr_tile%nx_tile_max
kmin=curr_tile%ny_tile_min
kmax=curr_tile%ny_tile_max
lmin=curr_tile%nz_tile_min
lmax=curr_tile%nz_tile_max
jminc=jmin-nxjguards
jmaxc=jmax+nxjguards
kminc=kmin-nyjguards
kmaxc=kmax+nyjguards
lminc=lmin-nzjguards
lmaxc=lmax+nzjguards
nxc=curr_tile%nx_cells_tile
nyc=curr_tile%ny_cells_tile
nzc=curr_tile%nz_cells_tile


rho(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = rho(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+curr_tile%rhotile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
rho(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = rho(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+curr_tile%rhotile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)

rho(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = rho(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+curr_tile%rhotile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
rho(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = rho(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+curr_tile%rhotile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)

rho(jmin:jmax,kmin:kmax,lminc:lmin-1) = rho(jmin:jmax,kmin:kmax,lminc:lmin-1)+curr_tile%rhotile(0:nxc, 0:nyc,-nzjguards:-1)
rho(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = rho(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+curr_tile%rhotile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)
end do
end do
end do
end do
!$omp end do
!$omp end parallel
call charge_bcs

end subroutine calc_diags
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_rho_scalar_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
use diagnostics
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: rho
real(kind=8) :: xp(np), yp(np), zp(np), w(np)
real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
real(kind=8) :: dxi,dyi,dzi,xint,yint,zint,oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
real(kind=8) :: x,y,z,wq,invvol
real(kind=8), dimension(2) :: sx(0:1), sy(0:1), sz(0:1)
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0
integer(idp) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
invvol = dxi*dyi*dzi
do ip=1,np

x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi

j=floor(x)
k=floor(y)
l=floor(z)

xint = x-j
yint = y-k
zint = z-l

wq=q*w(ip)*invvol

sx( 0) = 1.0_num-xint
sx( 1) = xint
sy( 0) = 1.0_num-yint
sy( 1) = yint
sz( 0) = 1.0_num-zint
sz( 1) = zint

rho(j,k,l)      = rho(j,k,l)+sx(0)*sy(0)*sz(0)*wq
rho(j+1,k,l)    = rho(j+1,k,l)+sx(1)*sy(0)*sz(0)*wq
rho(j,k+1,l)    = rho(j,k+1,l)+sx(0)*sy(1)*sz(0)*wq
rho(j+1,k+1,l)  = rho(j+1,k+1,l)+sx(1)*sy(1)*sz(0)*wq
rho(j,k,l+1)    = rho(j,k,l+1)+sx(0)*sy(0)*sz(1)*wq
rho(j+1,k,l+1)  = rho(j+1,k,l+1)+sx(1)*sy(0)*sz(1)*wq
rho(j,k+1,l+1)  = rho(j,k+1,l+1)+sx(0)*sy(1)*sz(1)*wq
rho(j+1,k+1,l+1)= rho(j+1,k+1,l+1)+sx(1)*sy(1)*sz(1)*wq
end do
return
end subroutine depose_rho_scalar_1_1_1
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_rho_vechv_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
use diagnostics
use constants
implicit none
integer(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8),intent(in out) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
real(kind=8), dimension(:,:), pointer:: rhocells
integer(idp), parameter :: lvec=8
integer(idp), dimension(lvec) :: icell
real(kind=8) :: ww
integer(idp) :: ncells
real(kind=8) :: xp(np), yp(np), zp(np), w(np)
real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
real(kind=8) :: dxi,dyi,dzi
real(kind=8) :: xint,yint,zint
real(kind=8) :: x,y,z,invvol
real(kind=8) :: sx(lvec), sy(lvec), sz(lvec), wq(lvec)
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0
integer(kind=4) :: ic,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
integer(kind=4) :: nnx, nnxy
integer(kind=4) :: moff(1:8)
real(kind=8):: mx(1:8),my(1:8),mz(1:8), sgn(1:8)


dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
invvol = dxi*dyi*dzi
ncells=(2*nxguard+nx)*(2*nyguard+ny)*(2*nzguard+nz)
allocate(rhocells(8,ncells))
rhocells=0.0_num
nnx = nx + 1 + 2*nxguard
nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
moff(1) = 0
moff(2) = 1
moff(3) = nnx
moff(4) = nnx+1
moff(5) = nnxy
moff(6) = nnxy+1
moff(7) = nnxy+nnx
moff(8) = nnxy+nnx+1
mx=(/0_num,1_num,0_num,1_num,0_num,1_num,0_num,1_num/)
my=(/0_num,0_num,1_num,1_num,0_num,0_num,1_num,1_num/)
mz=(/0_num,0_num,0_num,0_num,1_num,1_num,1_num,1_num/)
sgn=(/1_num,-1_num,-1_num,1_num,-1_num,1_num,1_num,-1_num/)


do ip=1,np,lvec




do n=1,min(lvec,np-ip+1)
nn=ip+n-1


x= (xp(nn)-xmin)*dxi
y = (yp(nn)-ymin)*dyi
z = (zp(nn)-zmin)*dzi

j=floor(x)
k=floor(y)
l=floor(z)
icell(n)=1+j+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)

sx(n) = x-j
sy(n) = y-k
sz(n) = z-l

wq(n)=q*w(nn)*invvol
end do

do n=1,min(lvec,np-ip+1)

ic=icell(n)



do nv=1,8
ww=(-mx(nv)+sx(n))*(-my(nv)+sy(n))*(-mz(nv)+sz(n))*wq(n)*sgn(nv)
rhocells(nv,ic)=rhocells(nv,ic)+ww
end do
end do
end do
do nv=1,8



do ic=1,ncells
rho(ic+moff(nv))=rho(ic+moff(nv))+rhocells(nv,ic)
end do
end do
deallocate(rhocells)
return
end subroutine depose_rho_vechv_1_1_1
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz,l_particles_weight, l4symtry)
use diagnostics
implicit none
integer(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: rho
real(kind=8) :: xp(np), yp(np), zp(np), w(np)
real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
logical :: l_particles_weight, l4symtry

real(kind=8) :: dxi,dyi,dzi,xint,yint,zint,oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
real(kind=8) :: x,y,z,wq,invvol
real(kind=8) :: sx(-int(nox/2):int((nox+1)/2)),sy(-int(noy/2):int((noy+1)/2)),sz(-int(noz/2):int((noz+1)/2))
real(kind=8), parameter :: onesixth=1.0/6.0,twothird=2.0/3.0
integer(idp) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
invvol = dxi*dyi*dzi

ixmin = -int(nox/2)
ixmax = int((nox+1)/2)
iymin = -int(noy/2)
iymax = int((noy+1)/2)
izmin = -int(noz/2)
izmax = int((noz+1)/2)

do ip=1,np


x = (xp(ip)-xmin)*dxi
y = (yp(ip)-ymin)*dyi
z = (zp(ip)-zmin)*dzi


if (l4symtry) then
x=abs(x)
y=abs(y)
end if



if (nox==2*(nox/2)) then
j=nint(x)
else
j=floor(x)
end if
if (noy==2*(noy/2)) then
k=nint(y)
else
k=floor(y)
end if
if(noz==2*(noz/2)) then
l=nint(z)
else
l=floor(z)
end if


xint = x-j
yint = y-k
zint = z-l


if (l_particles_weight) then
wq=q*w(ip)*invvol
else
wq=q*invvol*w(1)
endif


select case(nox)
case(0)
sx( 0) = 1.0_num
case(1)
sx( 0) = 1.0_num-xint
sx( 1) = xint
case(2)
xintsq = xint*xint
sx(-1) = 0.5_num*(0.5_num-xint)**2
sx( 0) = 0.75_num-xintsq
sx( 1) = 0.5_num*(0.5_num+xint)**2
case(3)
oxint = 1.0_num-xint
xintsq = xint*xint
oxintsq = oxint*oxint
sx(-1) = onesixth*oxintsq*oxint
sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
sx( 2) = onesixth*xintsq*xint
end select

select case(noy)
case(0)
sy( 0) = 1.0_num
case(1)
sy( 0) = 1.0_num-yint
sy( 1) = yint
case(2)
yintsq = yint*yint
sy(-1) = 0.5_num*(0.5_num-yint)**2
sy( 0) = 0.75_num-yintsq
sy( 1) = 0.5_num*(0.5_num+yint)**2
case(3)
oyint = 1.0_num-yint
yintsq = yint*yint
oyintsq = oyint*oyint
sy(-1) = onesixth*oyintsq*oyint
sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
sy( 2) = onesixth*yintsq*yint
end select

select case(noz)
case(0)
sz( 0) = 1.0_num
case(1)
sz( 0) = 1.0_num-zint
sz( 1) = zint
case(2)
zintsq = zint*zint
sz(-1) = 0.5_num*(0.5_num-zint)**2
sz( 0) = 0.75_num-zintsq
sz( 1) = 0.5_num*(0.5_num+zint)**2
case(3)
ozint = 1.0_num-zint
zintsq = zint*zint
ozintsq = ozint*ozint
sz(-1) = onesixth*ozintsq*ozint
sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
sz( 2) = onesixth*zintsq*zint
end select


do ll = izmin, izmax
do kk = iymin, iymax
do jj = ixmin, ixmax
rho(j+jj,k+kk,l+ll)=rho(j+jj,k+kk,l+ll)+sx(jj)*sy(kk)*sz(ll)*wq
end do
end do
end do
end do
return
end subroutine depose_rho_n
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine calc_field_div(divee, eex, eey, eez, nx, ny, nz, nxguard, nyguard, nzguard, dx, dy, dz)
use diagnostics
implicit none
integer(idp) ::  j,k,l
integer(idp) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: eex,eey,eez
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: divee
real(kind=8) :: dx, dy, dz, invdx, invdy, invdz

invdx=1.0_num/dx
invdy=1.0_num/dy
invdz=1.0_num/dz

do l = 0, nz
do k = 0, ny
do j = 0, nx
divee(j,k,l) = invdx*(eex(j,k,l)-eex(j-1,k,l))+invdy*(eey(j,k,l)-eey(j,k-1,l))+invdz*(eez(j,k,l)-eez(j,k,l-1))
end do
end do
end do

end subroutine calc_field_div
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine step(nst)

use constants
use fields
use particles
use params
use shared_data
use boundary
use omp_lib
use diagnostics
use simple_io

implicit none
integer :: nst,i


if (rank .eq. 0) then
write (0,*), "nsteps = ", nst
end if
do i=1,nst
if (rank .eq. 0) startit=mpi_wtime()
pushtime=0._num


call push_particles


call particle_bcs


call depose_currents_on_grid_jxjyjz


call current_bcs


call push_bfield


call bfield_bcs


call push_efield


call efield_bcs


call push_bfield


call bfield_bcs


call calc_diags




it = it+1
timeit=mpi_wtime()

if (rank .eq. 0) then
write(0,*) 'it = ',it,' || time = ',it*dt, " || push/part (ns)= ", pushtime*1e9_num/ntot," || tot/part (ns)= ", (timeit-startit)*1e9_num/ntot
end if
end do

end subroutine step
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine initall

use constants
use params
use fields
use particles
use shared_data
use tiling


implicit none
integer :: i,ierror,j,k,l, ispecies, ipart, count
integer :: jmin, jmax, lmin, lmax, kmin, kmax
integer :: ix, iy, iz
integer :: npartemp, ncurr
real(kind=8) :: v, th, phi
type(particle_species), pointer :: curr
type(particle_tile), pointer :: curr_tile



dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))
it = 0


if (.not. l_coeffs_allocated) then
allocate(xcoeffs(norderx/2),ycoeffs(nordery/2),zcoeffs(norderz/2))
end if


call fd_weights(xcoeffs, norderx, l_nodalgrid)
call fd_weights(ycoeffs, nordery, l_nodalgrid)
call fd_weights(zcoeffs, norderz, l_nodalgrid)




call set_tile_split


call init_tile_arrays


call load_particles




ex=0.0_num;ey=0.0_num;ez=0.0_num
bx=0.0_num;by=0.0_num;bz=0.0_num
jx=0.0_num;jy=0.0_num;jz=0.0_num


nsteps = nint(tmax/(w0_l*dt))

end subroutine initall
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine fd_weights(coeffs, norder, l_nodal)







use constants
implicit none
integer(kind=4) :: norder, n, m, mn, i, j, k
logical :: l_nodal
real(kind=8) :: z, fact, c1, c2, c3, c4, c5
real(kind=8), intent(in out), dimension(norder/2) :: coeffs
real(kind=8), pointer, dimension(:) :: x
real(kind=8), pointer, dimension(:,:) :: c

if (l_nodal) then
z=0.0_num
fact=1.0_num
else
z=0.5_num
fact=0.5_num
end if
m=1
n=norder+1

allocate(x(0:n-1))
allocate(c(0:m,0:n-1))

do i=0, n-1
x(i)=(i-n/2+1)*1.0_num
end do

c=0.0_num; c1=1.0_num; c4=x(0)-z; c(0,0)=1.0_num
do i=1, n-1
mn=min(i+1,m+1)
c2=1.0_num
c5=c4
c4=x(i)-z
do j=0, i-1
c3=x(i)-x(j)
c2=c2*c3
if (j .eq. (i-1)) then
do k=1, mn-1
c(k,i)=c1*(k*c(k-1,i-1)-c5*c(k,i-1))/c2
end do
c(0,i)=-c1*c5*c(0,i-1)/c2
do k=1, mn-1
c(k,j)=(c4*c(k,j)-k*c(k-1,j))/c3
end do
c(0,j)=c4*c(0,j)/c3
end if

end do
c1=c2
end do

do i=1, norder/2
coeffs(i)=c(m,norder/2+i-1)
end do
return
end subroutine fd_weights
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine mpi_minimal_init()
use mpi_routines
logical(isp) :: isinitialized
call mpi_initialized(isinitialized,errcode)
if (.not. isinitialized) call mpi_init_thread(mpi_thread_single,provided,errcode)
call mpi_comm_dup(mpi_comm_world, comm, errcode)
call mpi_comm_size(comm, nproc, errcode)
call mpi_comm_rank(comm, rank, errcode)
end subroutine mpi_minimal_init
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine setup_communicator
use mpi_routines

integer(isp), parameter :: ndims = 3
integer(isp) :: dims(ndims), idim, old_comm, ierr
logical(isp) :: periods(ndims), reorder, op, reset
integer(isp) :: test_coords(ndims)
integer(isp) :: ix, iy, iz
integer(idp) :: nxsplit, nysplit, nzsplit
integer(isp) :: ranges(3,1), nproc_orig, oldgroup, newgroup
character(len=11) :: str



call mpi_comm_size(mpi_comm_world, nproc, ierr)
dims = (/nprocz, nprocy, nprocx/)


if ((nprocx .eq. 0) .and. (nprocy .eq. 0) .and. (nprocz .eq. 0)) then
call mpi_dims_create(nproc, ndims, dims, errcode)
nprocx = dims(3)
nprocy = dims(2)
nprocz = dims(1)
endif

if (nproc .ne. nprocx*nprocy*nprocz) then
if (rank .eq. 0) then
write(0,*) '*** error ***'
write(0,*) 'nprocx*nprocy*nprocz =/ # of mpi processes'
write(0,*) ' check input file '
call mpi_abort(mpi_comm_world, errcode, ierr)
endif
endif

if (nx_global_grid .lt. nxguards .or. ny_global_grid .lt. nyguards .or. nz_global_grid .lt. nzguards) then
if (rank .eq. 0) then
write(0,*) '*** error ***'
write(0,*) 'simulation domain is too small.'
endif
call mpi_abort(mpi_comm_world, errcode, ierr)
endif

if (nprocx * nprocy * nprocz .gt. 0) then
nxsplit = nx_global_grid / nprocx
nysplit = ny_global_grid / nprocy
nzsplit = nz_global_grid / nprocz
if (nxsplit .lt. nxguards .or. nysplit .lt. nyguards .or. nzsplit .lt. nzguards) then
if (rank .eq. 0) then
write(0,*) 'wrong cpu split nlocal<nguards'
call mpi_abort(mpi_comm_world, errcode, ierr)
endif
endif
endif

periods = .false.
reorder = .true.



periods(c_ndims) = .true.
periods(c_ndims-1) = .true.
periods(c_ndims-2) = .true.

old_comm = comm
call mpi_cart_create(old_comm, ndims, dims, periods, reorder, comm, errcode)
call mpi_comm_free(old_comm, errcode)
call mpi_comm_rank(comm, rank, errcode)
call mpi_cart_coords(comm, rank, ndims, coordinates, errcode)
call mpi_cart_shift(comm, 2_isp, 1_isp, proc_x_min, proc_x_max, errcode)
call mpi_cart_shift(comm, 1_isp, 1_isp, proc_y_min, proc_y_max, errcode)
call mpi_cart_shift(comm, 0_isp, 1_isp, proc_z_min, proc_z_max, errcode)

nprocdir = dims

if (rank .eq. 0) then
write(0,*) 'processor subdivision is ', (/nprocx, nprocy, nprocz/)
endif

x_coords = coordinates(c_ndims)
x_min_boundary = .false.
x_max_boundary = .false.
if (x_coords .eq. 0) x_min_boundary = .true.
if (x_coords .eq. nprocx - 1) x_max_boundary = .true.

y_coords = coordinates(c_ndims-1)
y_min_boundary = .false.
y_max_boundary = .false.
if (y_coords .eq. 0) y_min_boundary = .true.
if (y_coords .eq. nprocy - 1) y_max_boundary = .true.

z_coords = coordinates(c_ndims-2)
z_min_boundary = .false.
z_max_boundary = .false.
if (z_coords .eq. 0) z_min_boundary = .true.
if (z_coords .eq. nprocz - 1) z_max_boundary = .true.

neighbour = mpi_proc_null
do iz = -1, 1
do iy = -1, 1
do ix = -1, 1
test_coords = coordinates
test_coords(1) = test_coords(1)+iz
test_coords(2) = test_coords(2)+iy
test_coords(3) = test_coords(3)+ix
op = .true.


do idim = 1, ndims
if ((test_coords(idim) .lt. 0.or. test_coords(idim) .ge. dims(idim)).and. .not. periods(idim)) op = .false.
enddo
if (op) then
call mpi_cart_rank(comm, test_coords, neighbour(ix,iy,iz), errcode)
endif
enddo
enddo
enddo

end subroutine setup_communicator
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine mpi_initialise
use mpi_routines

integer(isp) :: idim
integer(isp) :: nx0, nxp
integer(isp) :: ny0, nyp
integer(isp) :: nz0, nzp
integer(isp) :: iproc, ix, iy, iz


nxguards = max(nox,norderx)+npass(1)
nyguards = max(noy,nordery)+npass(2)
nzguards = max(noz,norderz)+npass(3)
nxjguards = max(nox,2)
nyjguards = max(noy,2)
nzjguards = max(noz,2)

if (l_smooth_compensate) then
nxguards = nxguards + 1
nyguards = nyguards + 1
nzguards = nzguards + 1
end if

call setup_communicator

allocate(npart_each_rank(nproc))
allocate(x_grid_mins(0:nprocx-1), x_grid_maxs(0:nprocx-1))
allocate(y_grid_mins(0:nprocy-1), y_grid_maxs(0:nprocy-1))
allocate(z_grid_mins(0:nprocz-1), z_grid_maxs(0:nprocz-1))
allocate(cell_x_min(nprocx), cell_x_max(nprocx))
allocate(cell_y_min(nprocy), cell_y_max(nprocy))
allocate(cell_z_min(nprocz), cell_z_max(nprocz))

nx0 = nx_global_grid / nprocx
ny0 = ny_global_grid / nprocy
nz0 = nz_global_grid / nprocz




if (nx0 * nprocx .ne. nx_global_grid) then
nxp = (nx0 + 1) * nprocx - nx_global_grid
else
nxp = nprocx
endif

if (ny0 * nprocy .ne. ny_global_grid) then
nyp = (ny0 + 1) * nprocy - ny_global_grid
else
nyp = nprocy
endif

if (nz0 * nprocz .ne. nz_global_grid) then
nzp = (nz0 + 1) * nprocz - nz_global_grid
else
nzp = nprocz
endif

do idim = 1, nxp
cell_x_min(idim) = (idim - 1) * nx0 + 1
cell_x_max(idim) = idim * nx0
enddo
do idim = nxp + 1, nprocx
cell_x_min(idim) = nxp * nx0 + (idim - nxp - 1) * (nx0 + 1) + 1
cell_x_max(idim) = nxp * nx0 + (idim - nxp) * (nx0 + 1)
enddo

do idim = 1, nyp
cell_y_min(idim) = (idim - 1) * ny0 + 1
cell_y_max(idim) = idim * ny0
enddo
do idim = nyp + 1, nprocy
cell_y_min(idim) = nyp * ny0 + (idim - nyp - 1) * (ny0 + 1) + 1
cell_y_max(idim) = nyp * ny0 + (idim - nyp) * (ny0 + 1)
enddo

do idim = 1, nzp
cell_z_min(idim) = (idim - 1) * nz0 + 1
cell_z_max(idim) = idim * nz0
enddo
do idim = nzp + 1, nprocz
cell_z_min(idim) = nzp * nz0 + (idim - nzp - 1) * (nz0 + 1) + 1
cell_z_max(idim) = nzp * nz0 + (idim - nzp) * (nz0 + 1)

enddo

nx_global_grid_min = cell_x_min(x_coords+1)
nx_global_grid_max = cell_x_max(x_coords+1)
n_global_grid_min(1) = nx_global_grid_min
n_global_grid_max(1) = nx_global_grid_max

ny_global_grid_min = cell_y_min(y_coords+1)
ny_global_grid_max = cell_y_max(y_coords+1)
n_global_grid_min(2) = ny_global_grid_min
n_global_grid_max(2) = ny_global_grid_max

nz_global_grid_min = cell_z_min(z_coords+1)
nz_global_grid_max = cell_z_max(z_coords+1)
n_global_grid_min(3) = nz_global_grid_min
n_global_grid_max(3) = nz_global_grid_max


nx_grid = nx_global_grid_max - nx_global_grid_min + 1
ny_grid = ny_global_grid_max - ny_global_grid_min + 1
nz_grid = nz_global_grid_max - nz_global_grid_min + 1


nx=nx_grid-1
ny=ny_grid-1
nz=nz_grid-1

allocate(x(-nxguards:nx+nxguards), y(-nyguards:ny+nyguards), z(-nzguards:nz+nzguards))
allocate(x_global(-nxguards:nx_global+nxguards))
allocate(y_global(-nyguards:ny_global+nyguards))
allocate(z_global(-nzguards:nz_global+nzguards))


xmax = nx_global*dx
ymax = ny_global*dy
zmax = nz_global*dz


length_x = xmax - xmin +dx
dx = length_x / real(nx_global+1, num)
x_grid_min = xmin
x_grid_max = xmax

length_y = ymax - ymin +dy
dy = length_y / real(ny_global+1, num)
y_grid_min = ymin
y_grid_max = ymax

length_z = zmax - zmin +dz
dz = length_z / real(nz_global+1, num)
z_grid_min = zmin
z_grid_max = zmax


do ix = -nxguards, nx_global+nxguards
x_global(ix) = x_grid_min + ix * dx
enddo
do iy = -nyguards, ny_global+nyguards
y_global(iy) = y_grid_min + iy * dy
enddo
do iz = -nzguards, nz_global+nzguards
z_global(iz) = z_grid_min + iz * dz
enddo


do iproc = 0, nprocx-1
x_grid_mins(iproc) = x_global(cell_x_min(iproc+1)-1)
x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1)-1)
enddo
do iproc = 0, nprocy-1
y_grid_mins(iproc) = y_global(cell_y_min(iproc+1)-1)
y_grid_maxs(iproc) = y_global(cell_y_max(iproc+1)-1)
enddo
do iproc = 0, nprocz-1
z_grid_mins(iproc) = z_global(cell_z_min(iproc+1)-1)
z_grid_maxs(iproc) = z_global(cell_z_max(iproc+1)-1)
enddo

x_min_local = x_grid_mins(x_coords)-dx/2
x_max_local = x_grid_maxs(x_coords)+dx/2
y_min_local = y_grid_mins(y_coords)-dx/2
y_max_local = y_grid_maxs(y_coords)+dx/2
z_min_local = z_grid_mins(z_coords)-dx/2
z_max_local = z_grid_maxs(z_coords)+dx/2

x_grid_min_local=x_min_local+dx/2
y_grid_min_local=y_min_local+dy/2
z_grid_min_local=z_min_local+dz/2
x_grid_max_local=x_max_local-dx/2
y_grid_max_local=y_max_local-dy/2
z_grid_max_local=z_max_local-dz/2


allocate(ex(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
allocate(ey(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
allocate(ez(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
allocate(bx(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
allocate(by(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
allocate(bz(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))
allocate(jx(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
allocate(jy(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
allocate(jz(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
allocate(rho(-nxjguards:nx+nxjguards, -nyjguards:ny+nyjguards, -nzjguards:nz+nzjguards))
allocate(dive(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards))

start_time = mpi_wtime()

end subroutine mpi_initialise
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine mpi_close
use mpi_routines

integer :: seconds, minutes, hours, total

if (rank .eq. 0) then
end_time = mpi_wtime()
total = int(end_time - start_time)
seconds = mod(total, 60)
minutes = mod(total / 60, 60)
hours = total / 3600
endif

call mpi_finalize(errcode)

end subroutine mpi_close
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine default_init
use control_file

ntilex = 1
ntiley = 1
ntilez = 1


norderx = 2
nordery = 2
norderz = 2
l_nodalgrid = .false.


nox = 1
noy = 1
noz = 1
l_lower_order_in_v = .false.


dtcoef = 0.7_num

npass = 0
alpha = 0.5_num

tmax = 40.0_num



l_particles_weight = .false.



nlab  = 1.e23_num
g0    = 130.0_num
b0    = sqrt(1.0_num-1.0_num/g0**2)
nc    = nlab*g0
wlab  = echarge*sqrt(nlab/(emass*eps0))
w0_l  = echarge*sqrt(nc/(g0*emass*eps0))
w0_t  = echarge*sqrt(nc/(g0**3*emass*eps0))
w0    = w0_l

nspecies=0


pdistr=1

if (.not. l_species_allocated) then
nspecies=0
allocate(species_parray(1:nspecies_max))
l_species_allocated=.true.
endif
end subroutine default_init
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine read_from_cl
use control_file
integer :: i, ix
do i = 1, iargc()-1,2
call getarg(i, buffer)
if (index(buffer,'ntilex') .gt. 0) then
call getarg(i+1, buffer)
read(buffer, '(i10)') ntilex
else if (index(buffer,'ntiley') .gt. 0) then
call getarg(i+1, buffer)
read(buffer, '(i10)') ntiley
else if (index(buffer,'ntilez') .gt. 0) then
call getarg(i+1, buffer)
read(buffer, '(i10)') ntilez
else if (index(buffer,'distr') .gt. 0) then
call getarg(i+1, buffer)
read(buffer, '(i10)') pdistr
end if
end do
return
end subroutine read_from_cl
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine read_input_file
use control_file
integer :: ix = 0

open(fh_input, file='input_file.pixr')
do while(ios==0)
read(fh_input, '(a)', iostat=ios) buffer
ix=index(buffer,'section::')
if (ix .gt. 0) then
section_name=buffer(ix:string_length)
select case(trim(adjustl(section_name)))
case('section::main')
call read_main_section
case('section::species')
call read_species_section
case('section::output')
call read_output_section
case('section::cpusplit')
call read_cpusplit_section
end select
end if
end do
return
end subroutine read_input_file
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine read_cpusplit_section
use control_file
integer :: ix = 0
logical :: end_section = .false.

do while((.not. end_section) .and. (ios==0))
read(fh_input, '(a)', iostat=ios) buffer
if (index(buffer,'nprocx') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') nprocx
else if (index(buffer,'nprocy') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') nprocy
else if (index(buffer,'nprocz') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') nprocz
else if (index(buffer,'end::cpusplit') .gt. 0) then
end_section =.true.
end if
end do
return
end subroutine read_cpusplit_section
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine read_main_section
use control_file
integer :: ix = 0
logical :: end_section = .false.

do while((.not. end_section) .and. (ios==0))
read(fh_input, '(a)', iostat=ios) buffer
if (index(buffer,'nx') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') nx_global_grid
nx_global=nx_global_grid-1
else if (index(buffer,'ny') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') ny_global_grid
ny_global=ny_global_grid-1
else if (index(buffer,'nz') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') nz_global_grid
nz_global=nz_global_grid-1
else if (index(buffer,'ntilex') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') ntilex
else if (index(buffer,'ntiley') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') ntiley
else if (index(buffer,'ntilez') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') ntilez
elseif (index(buffer,'dx') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) dx
else if (index(buffer,'dy') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) dy
else if (index(buffer,'dz') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) dz
elseif (index(buffer,'x_grid_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) xmin
else if (index(buffer,'y_grid_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) ymin
else if (index(buffer,'z_grid_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) zmin
else if (index(buffer,'t_max') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) tmax
else if (index(buffer,'end::main') .gt. 0) then
end_section =.true.
end if
end do
return
end subroutine read_main_section
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine read_species_section
use control_file
integer :: ix = 0
logical :: end_section
type(particle_species), pointer :: curr

if (.not. l_species_allocated) then
nspecies=0
allocate(species_parray(1:nspecies_max))
l_species_allocated=.true.
endif
nspecies = nspecies+1
curr => species_parray(nspecies)

curr%charge = -echarge
curr%mass = emass
curr%nppcell = 0
curr%x_min = 0._num
curr%x_max = 0._num
curr%y_min = 0._num
curr%y_max = 0._num
curr%z_min = 0._num
curr%z_max = 0._num
curr%vdrift_x =0._num
curr%vdrift_y =0._num
curr%vdrift_z =0._num
curr%vth_x =0._num
curr%vth_y =0._num
curr%vth_z =0._num
curr%species_npart=0
end_section=.false.
do while((.not. end_section) .and. (ios==0))
read(fh_input, '(a)', iostat=ios) buffer
if (index(buffer,'name') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%name
else if (index(buffer,'mass') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%mass
curr%mass=curr%mass*emass
else if (index(buffer,'charge') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%charge
curr%charge=curr%charge*echarge
elseif (index(buffer,'nppcell') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length),'(i10)') curr%nppcell
else if (index(buffer,'x_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%x_min
else if (index(buffer,'x_max') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%x_max
elseif (index(buffer,'y_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%y_min
else if (index(buffer,'y_max') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%y_max
elseif (index(buffer,'z_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%z_min
else if (index(buffer,'z_max') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%z_max
else if (index(buffer,'vdrift_x') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%vdrift_x
curr%vdrift_x=curr%vdrift_x*clight
else if (index(buffer,'vdrift_y') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%vdrift_y
curr%vdrift_y=curr%vdrift_y*clight
else if (index(buffer,'vdrift_z') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%vdrift_z
curr%vdrift_z=curr%vdrift_z*clight
else if (index(buffer,'vth_x') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%vth_x
curr%vth_x=curr%vth_x*clight
else if (index(buffer,'vth_y') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%vth_y
curr%vth_y=curr%vth_y*clight
else if (index(buffer,'vth_z') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), *) curr%vth_z
curr%vth_z=curr%vth_z*clight
else if (index(buffer,'end::species') .gt. 0) then
end_section =.true.
end if
end do
return
end subroutine read_species_section
!----------------------
! THIS IS A SUBROUTINE 
!----------------------
subroutine read_output_section
use control_file
integer :: ix = 0
logical :: end_section = .false.

do while((.not. end_section) .and. (ios==0))
read(fh_input, '(a)', iostat=ios) buffer
if (index(buffer,'output_frequency') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') output_frequency
else if (index(buffer,'output_step_min') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') output_step_min
else if (index(buffer,'output_step_max') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') output_step_max
elseif (index(buffer,'ex') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_ex
else if (index(buffer,'ey') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_ey
else if (index(buffer,'ez') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_ez
elseif (index(buffer,'bx') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_bx
else if (index(buffer,'by') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_by
else if (index(buffer,'bz') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_bz
elseif (index(buffer,'jx') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_jx
else if (index(buffer,'jy') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_jy
else if (index(buffer,'jz') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_jz
elseif (index(buffer,'rho') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_rho
else if (index(buffer,'dive') .gt. 0) then
ix = index(buffer, "=")
read(buffer(ix+1:string_length), '(i10)') c_output_dive
else if (index(buffer,'end::output') .gt. 0) then
end_section =.true.
end if
end do
return
end subroutine read_output_section

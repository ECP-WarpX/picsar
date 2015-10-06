!===============================================================================
! Gather electric field on particles
!===============================================================================
SUBROUTINE gather_ebfields_on_particles
USE particles
USE constants
USE fields
USE params
USE shared_data
USE tiling
IMPLICIT NONE
INTEGER :: ispecies, ix, iy, iz, count
INTEGER :: jmin, jmax, kmin, kmax, lmin, lmax
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend
REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ex_tile, ey_tile, ez_tile, bx_tile, by_tile, bz_tile
INTEGER :: nxc, nyc, nzc

tdeb=MPI_WTIME()
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray, &
!$OMP nxjguards,nyjguards,nzjguards,ex,ey,ez,bx,by,bz,dx,dy,dz) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,ex_tile,ey_tile,ez_tile,bx_tile,by_tile,bz_tile,nxc,nyc,nzc) 
DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
        DO ix=1, ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
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
                ALLOCATE(ex_tile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),  &
                         ey_tile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),  &
                         ez_tile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),  &
                         bx_tile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),  &
                         by_tile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards),  &
                         bz_tile(-nxjguards:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards))
                !!! --- Gather electric fields on particles
                curr_tile%part_ex = 0.0_num
                curr_tile%part_ey = 0.0_num
                curr_tile%part_ez = 0.0_num
                ex_tile = ex(jmin:jmax,kmin:kmax,lmin:lmax)
                ey_tile = ey(jmin:jmax,kmin:kmax,lmin:lmax)
                ez_tile = ez(jmin:jmax,kmin:kmax,lmin:lmax)
                CALL gete3d_energy_conserving_vec_1_1_1(count,curr_tile%part_x(1:count),curr_tile%part_y(1:count), &
                                      curr_tile%part_z(1:count), curr_tile%part_ex(1:count),                   &
                                      curr_tile%part_ey(1:count),curr_tile%part_ez(1:count),                   &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                     &
                                      curr_tile%z_grid_tile_min, dx,dy,dz,curr_tile%nx_cells_tile,             &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,     &
                                      nzjguards,ex_tile,ey_tile,                           &
                                      ez_tile)

                !!! --- Gather magnetic fields on particles
                curr_tile%part_bx=0.0_num
                curr_tile%part_by=0.0_num
                curr_tile%part_bz=0.0_num
                bx_tile = bx(jmin:jmax,kmin:kmax,lmin:lmax)
                by_tile = by(jmin:jmax,kmin:kmax,lmin:lmax)
                bz_tile = bz(jmin:jmax,kmin:kmax,lmin:lmax)
                CALL getb3d_energy_conserving_1_1_1(count,curr_tile%part_x(1:count),curr_tile%part_y(1:count), &
                                      curr_tile%part_z(1:count), curr_tile%part_bx(1:count),                   &
                                      curr_tile%part_by(1:count),curr_tile%part_bz(1:count),                   &
                                      curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,                     &
                                      curr_tile%z_grid_tile_min, dx,dy,dz,curr_tile%nx_cells_tile,             &
                                      curr_tile%ny_cells_tile,curr_tile%nz_cells_tile,nxjguards,nyjguards,     &
                                      nzjguards,bx_tile,by_tile,                           &
                                      bz_tile)
                DEALLOCATE(ex_tile, ey_tile, ez_tile, bx_tile,by_tile,bz_tile)
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO! END LOOP ON TILES
!$OMP END PARALLEL DO
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)

END SUBROUTINE gather_ebfields_on_particles


!=================================================================================
! Gathering of electric field from Yee grid ("energy conserving") on particles
! at order 1
SUBROUTINE gete3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg)
!=================================================================================

USE omp_lib
USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,ex,ey,ez
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1) :: sx
REAL(num), DIMENSION(0:1) :: sy
REAL(num), DIMENSION(0:1) :: sz
REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
ALLOCATE(sx0(0:1),sy0(0:1),sz0(0:1))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np
    
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=floor(x)
    j0=floor(x-0.5_num)
    k=floor(y)
    k0=floor(y-0.5_num)
    l=floor(z)
    l0=floor(z-0.5_num)
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
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
    
    ! Compute Ex on particle
    ex(ip) = ex(ip) + sx(0)*sy0(0)*sz0(0)*exg(j0,k,l)
    ex(ip) = ex(ip) + sx(1)*sy0(0)*sz0(0)*exg(j0+1,k,l)
    ex(ip) = ex(ip) + sx(0)*sy0(1)*sz0(0)*exg(j0,k+1,l)
    ex(ip) = ex(ip) + sx(1)*sy0(1)*sz0(0)*exg(j0+1,k+1,l)
    ex(ip) = ex(ip) + sx(0)*sy0(0)*sz0(1)*exg(j0,k,l+1)
    ex(ip) = ex(ip) + sx(1)*sy0(0)*sz0(1)*exg(j0+1,k,l+1)
    ex(ip) = ex(ip) + sx(0)*sy0(1)*sz0(1)*exg(j0,k+1,l+1)
    ex(ip) = ex(ip) + sx(1)*sy0(1)*sz0(1)*exg(j0+1,k+1,l+1)
    
    ! Compute Ey on particle
    ey(ip) = ey(ip) + sx0(0)*sy(0)*sz0(0)*eyg(j,k0,l)
    ey(ip) = ey(ip) + sx0(1)*sy(0)*sz0(0)*eyg(j+1,k0,l)
    ey(ip) = ey(ip) + sx0(0)*sy(1)*sz0(0)*eyg(j,k0+1,l)
    ey(ip) = ey(ip) + sx0(1)*sy(1)*sz0(0)*eyg(j+1,k0+1,l)
    ey(ip) = ey(ip) + sx0(0)*sy(0)*sz0(1)*eyg(j,k0,l+1)
    ey(ip) = ey(ip) + sx0(1)*sy(0)*sz0(1)*eyg(j+1,k0,l+1)
    ey(ip) = ey(ip) + sx0(0)*sy(1)*sz0(1)*eyg(j,k0+1,l+1)
    ey(ip) = ey(ip) + sx0(1)*sy(1)*sz0(1)*eyg(j+1,k0+1,l+1)
    
    ! Compute Ez on particle
    ez(ip) = ez(ip) + sx0(0)*sy0(0)*sz(0)*ezg(j,k,l0)
    ez(ip) = ez(ip) + sx0(1)*sy0(0)*sz(0)*ezg(j+1,k,l0)
    ez(ip) = ez(ip) + sx0(0)*sy0(1)*sz(0)*ezg(j,k+1,l0)
    ez(ip) = ez(ip) + sx0(1)*sy0(1)*sz(0)*ezg(j+1,k+1,l0)
    ez(ip) = ez(ip) + sx0(0)*sy0(0)*sz(1)*ezg(j,k,l0+1)
    ez(ip) = ez(ip) + sx0(1)*sy0(0)*sz(1)*ezg(j+1,k,l0+1)
    ez(ip) = ez(ip) + sx0(0)*sy0(1)*sz(1)*ezg(j,k+1,l0+1)
    ez(ip) = ez(ip) + sx0(1)*sy0(1)*sz(1)*ezg(j+1,k+1,l0+1)
END DO
!!$OMP END PARALLEL DO
DEALLOCATE(sx0,sz0)
RETURN
END SUBROUTINE gete3d_energy_conserving_1_1_1

!=================================================================================
! Gathering of electric field from Yee grid ("energy conserving") on particles
! at order 1
SUBROUTINE gete3d_energy_conserving_vec_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg)
!=================================================================================

USE omp_lib
USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,ex,ey,ez
REAL(num), DIMENSION(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard)) :: exg,eyg,ezg
REAL(num), DIMENSION(:,:), ALLOCATABLE :: excells, eycells, ezcells
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, wwx,wwy,wwz
INTEGER, PARAMETER :: LVEC=8
REAL(num) :: sx(LVEC), sy(LVEC), sz(LVEC), sx0(LVEC), sy0(LVEC), sz0(LVEC), wq(LVEC)
INTEGER :: ICELLSx(LVEC), ICELLSy(LVEC), ICELLSz(LVEC)
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
INTEGER :: nnx, nnxy, ic, NCELLS, icx,icy,icz, n, nv, nn
INTEGER   :: moff(1:8)
REAL(num) :: mx(1:8),my(1:8),mz(1:8), sgn(1:8)
REAL(num), DIMENSION(8,LVEC) :: extemp,eytemp,eztemp

NCELLS=(2*nxguard+nx)*(2*nyguard+ny)*(2*nzguard+nz)
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
ALLOCATE(excells(8,NCELLS),eycells(8,NCELLS),ezcells(8,NCELLS))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
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
mx=(/0,1,0,1,0,1,0,1/)
my=(/0,0,1,1,0,0,1,1/)
mz=(/0,0,0,0,1,1,1,1/)
sgn=(/1,-1,-1,1,-1,1,1,-1/)

! Reduction of ex,ey,ez into excells,eycells,ezcells
DO ic=1,NCELLS  !!! VECTOR
    !DIR$ ASSUME_ALIGNED excells:64, eycells:64, ezcells:64
    !DIR$ IVDEP
    DO nv=1,8
        excells(nv,ic)=excells(nv,ic)+exg(ic+moff(nv))
        eycells(nv,ic)=eycells(nv,ic)+eyg(ic+moff(nv))
        ezcells(nv,ic)=ezcells(nv,ic)+ezg(ic+moff(nv))
    END DO
END DO

!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np, LVEC
    !DIR$ ASSUME_ALIGNED sx:64, sy:64, sz:64
    !DIR$ ASSUME_ALIGNED sx0:64, sy0:64, sz0:64
    !DIR$ ASSUME_ALIGNED xp:64, yp:64, zp:64
    !DIR$ ASSUME_ALIGNED ICELLSx:64, ICELLSy:64, ICELLSz:64
    !DIR$ IVDEP
    DO n=1, MIN(LVEC,np-ip+1)
        nn=ip+n-1
        x = (xp(nn)-xmin)*dxi
        y = (yp(nn)-ymin)*dyi
        z = (zp(nn)-zmin)*dzi
    
        ! Compute index of particle
        j=floor(x)
        j0=floor(x-0.5_num)
        k=floor(y)
        k0=floor(y-0.5_num)
        l=floor(z)
        l0=floor(z-0.5_num)
        ICELLSx(n)=1+j0+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)
        ICELLSy(n)=1+j+nxguard+(k0+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)
        ICELLSz(n)=1+j+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l0+nzguard+1)*(ny+2*nyguard)

        ! Particle relative position to nodes
        sx(n)=x-j
        sy(n)=y-k
        sz(n)=z-l
        sx0(n)=x-0.5_num-j0
        sy0(n)=y-0.5_num-k0
        sz0(n)=z-0.5_num-l0
    END DO
    DO n=1, MIN(LVEC,np-ip+1)
        icx=ICELLSx(n)
        icy=ICELLSy(n)
        icz=ICELLSz(n)
        !DIR$ ASSUME_ALIGNED mx:64, my:64, mz:64,sgn:64
        !DIR$ ASSUME_ALIGNED excells:64, eycells:64, ezcells:64
        !DIR$ SIMD
        DO nv=1,8
            ! Compute weight for vertex nv of cells icx,icy,icz
            wwx=(-mx(nv)+sx(n))*(-my(nv)+sy0(n))* &
            (-mz(nv)+sz0(n))*sgn(nv)
            wwy=(-mx(nv)+sx0(n))*(-my(nv)+sy(n))* &
            (-mz(nv)+sz0(n))*sgn(nv)
            wwz=(-mx(nv)+sx0(n))*(-my(nv)+sy0(n))* &
            (-mz(nv)+sz(n))*sgn(nv)
            ! Compute Ex on particle
            extemp(nv,n) = extemp(nv,n) + wwx*excells(nv,icx)
    
            ! Compute Ey on particle
            eytemp(nv,n) = eytemp(nv,n) + wwy*eycells(nv,icy)
    
            ! Compute Ez on particle
            eztemp(nv,n) = eztemp(nv,n) + wwz*ezcells(nv,icz)
        END DO
    END DO 
    DO nv=1,8
        !DIR$ ASSUME_ALIGNED ex:64, ey:64, ez:64
        !DIR$ ASSUME_ALIGNED extemp:64,eytemp:64,eztemp:64
        DO n=1,MIN(LVEC,np-ip+1)
           nn=ip+n-1
           ex(nn)=ex(nn)+extemp(nv,n)
           ey(nn)=ey(nn)+eytemp(nv,n)
           ez(nn)=ez(nn)+eztemp(nv,n)
        END DO
    END DO
END DO
!!$OMP END PARALLEL DO
DEALLOCATE(excells,eycells,ezcells)
RETURN
END SUBROUTINE gete3d_energy_conserving_vec_1_1_1


!=================================================================================
! Gathering of Magnetic field from Yee grid ("energy conserving") on particles
! at order 1
SUBROUTINE getb3d_energy_conserving_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg)
!=================================================================================

USE omp_lib
USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,bx,by,bz
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
REAL(num), DIMENSION(0:1) :: sx
REAL(num), DIMENSION(0:1) :: sy
REAL(num), DIMENSION(0:1) :: sz
REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
ALLOCATE(sx0(0:1),sy0(0:1),sz0(0:1))
sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, &
!!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np
    
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=floor(x)
    j0=floor(x-0.5_num)
    k=floor(y)
    k0=floor(y-0.5_num)
    l=floor(z)
    l0=floor(z-0.5_num)
    xint=x-j
    yint=y-k
    zint=z-l
    
    ! Compute shape factors
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
    
    ! Compute Bx on particle
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(0)*bxg(j,k0,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(0)*bxg(j+1,k0,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(0)*bxg(j,k0+1,l0)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(0)*bxg(j+1,k0+1,l0)
    bx(ip) = bx(ip) + sx(0)*sy0(0)*sz0(1)*bxg(j,k0,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(0)*sz0(1)*bxg(j+1,k0,l0+1)
    bx(ip) = bx(ip) + sx(0)*sy0(1)*sz0(1)*bxg(j,k0+1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sy0(1)*sz0(1)*bxg(j+1,k0+1,l0+1)
    
    ! Compute By on particle
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(0)*byg(j0,k,l0)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(0)*byg(j0+1,k,l0)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(0)*byg(j0,k+1,l0)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(0)*byg(j0+1,k+1,l0)
    by(ip) = by(ip) + sx0(0)*sy(0)*sz0(1)*byg(j0,k,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(0)*sz0(1)*byg(j0+1,k,l0+1)
    by(ip) = by(ip) + sx0(0)*sy(1)*sz0(1)*byg(j0,k+1,l0+1)
    by(ip) = by(ip) + sx0(1)*sy(1)*sz0(1)*byg(j0+1,k+1,l0+1)
    
    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(0)*bzg(j0,k0,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(0)*bzg(j0+1,k0,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(0)*bzg(j0,k0+1,l)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(0)*bzg(j0+1,k0+1,l)
    bz(ip) = bz(ip) + sx0(0)*sy0(0)*sz(1)*bzg(j0,k0,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(0)*sz(1)*bzg(j0+1,k0,l+1)
    bz(ip) = bz(ip) + sx0(0)*sy0(1)*sz(1)*bzg(j0,k0+1,l+1)
    bz(ip) = bz(ip) + sx0(1)*sy0(1)*sz(1)*bzg(j0+1,k0+1,l+1)
END DO
!!$OMP END PARALLEL DO
DEALLOCATE(sx0,sz0)
RETURN
END SUBROUTINE getb3d_energy_conserving_1_1_1


!=================================================================================
! Gathering of electric field from Yee grid ("energy conserving") on particles
! At arbitrary order. WARNING: Highly unoptimized routine
SUBROUTINE gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,       &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,exg,eyg,ezg,l_lower_order_in_v)
!=================================================================================
USE omp_lib
USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), dimension(np) :: xp,yp,zp,ex,ey,ez
LOGICAL :: l4symtry,l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
REAL(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
REAL(num), dimension(:), allocatable :: sx0,sy0,sz0
REAL(num), parameter :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
iymin = -int(noy/2)
iymax =  int((noy+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1

IF (l_lower_order_in_v) THEN
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
ELSE
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
END IF

ALLOCATE(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

signx = 1.0_num
signy = 1.0_num
!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint, &
!$OMP   sx,sy,sz,sx0,sy0,sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq,ozint,zintsq,ozintsq)
DO ip=1,np

    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    IF (l_lower_order_in_v) THEN
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x-0.5_num)
        ELSE
            j=floor(x)
            j0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y-0.5_num)
        ELSE
            k=floor(y)
            k0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z-0.5_num)
        ELSE
            l=floor(z)
            l0=floor(z)
        END IF
    ELSE
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x)
        ELSE
            j=floor(x)
            j0=floor(x-0.5_num)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y)
        ELSE
            k=floor(y)
            k0=floor(y-0.5_num)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z)
        ELSE
            l=floor(z)
            l0=floor(z-0.5_num)
        END IF
    END IF

    xint=x-j
    yint=y-k
    zint=z-l

    IF (nox==1) THEN
        sx( 0) = 1.0_num-xint
        sx( 1) = xint
    ELSEIF (nox==2) THEN
        xintsq = xint*xint
        sx(-1) = 0.5_num*(0.5_num-xint)**2
        sx( 0) = 0.75_num-xintsq
        sx( 1) = 0.5_num*(0.5_num+xint)**2
    ELSEIF (nox==3) THEN
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx( 2) = onesixth*xintsq*xint
    END IF

    IF (noy==1) THEN
        sy( 0) = 1.0_num-yint
        sy( 1) = yint
    ELSEIF (noy==2) THEN
        yintsq = yint*yint
        sy(-1) = 0.5_num*(0.5_num-yint)**2
        sy( 0) = 0.75_num-yintsq
        sy( 1) = 0.5_num*(0.5_num+yint)**2
    ELSEIF (noy==3) THEN
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy( 2) = onesixth*yintsq*yint
    END IF

    IF (noz==1) THEN
        sz( 0) = 1.0_num-zint
        sz( 1) = zint
    ELSEIF (noz==2) THEN
        zintsq = zint*zint
        sz(-1) = 0.5_num*(0.5_num-zint)**2
        sz( 0) = 0.75_num-zintsq
        sz( 1) = 0.5_num*(0.5_num+zint)**2
    ELSEIF (noz==3) THEN
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz( 2) = onesixth*zintsq*zint
    END IF

    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0

    IF (l_lower_order_in_v) THEN

        IF (nox==1) THEN
            sx0( 0) = 1.0_num
        ELSEIF (nox==2) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==3) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num
        ELSEIF (noy==2) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==3) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        END IF

        IF (noz==1) THEN
            sz0( 0) = 1.0_num
        ELSEIF (noz==2) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==3) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        END IF

    ELSE

        IF (nox==1) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==2) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        ELSEIF (nox==3) THEN
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==2) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        ELSEIF (noy==3) THEN
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END IF

        IF (noz==1) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==2) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        ELSEIF (noz==3) THEN
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END IF

    END IF

    DO ll = izmin, izmax+1
        DO kk = iymin, iymax+1
            DO jj = ixmin0, ixmax0
                ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj,k+kk,l+ll)*signx
            END DO
        END DO
    END DO

    DO ll = izmin, izmax+1
        DO kk = iymin0, iymax0
            DO jj = ixmin, ixmax+1
                ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k0+kk,l+ll)*signy
            END DO
        END DO
    END DO

    DO ll = izmin0, izmax0
        DO kk = iymin, iymax+1
            DO jj = ixmin, ixmax+1
                ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj,k+kk,l0+ll)
            END DO
        END DO
    END DO

END DO
!$OMP END PARALLEL DO
DEALLOCATE(sx0,sy0,sz0)

RETURN
END SUBROUTINE gete3d_n_energy_conserving

!=================================================================================
! Gathering of Magnetic field from Yee grid ("energy conserving") on particles
! At arbitrary order. WARNING: Highly unoptimized routine
SUBROUTINE getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,       &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,bxg,byg,bzg,l_lower_order_in_v)
!=================================================================================
USE omp_lib
USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(np) :: xp,yp,zp,bx,by,bz
LOGICAL :: l_lower_order_in_v
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
REAL(num) :: xmin,ymin,zmin,dx,dy,dz
INTEGER :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
              ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0
REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
              xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
REAL(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
REAL(num), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
REAL(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num

dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz

ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
iymin = -int(noy/2)
iymax =  int((noy+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1


IF (l_lower_order_in_v) THEN
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    iymin0 = -int((noy-1)/2)
    iymax0 =  int((noy)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
ELSE
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    iymin0 = -int((noy)/2)
    iymax0 =  int((noy+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)
END IF
ALLOCATE(sx0(ixmin0:ixmax0),sy0(iymin0:iymax0),sz0(izmin0:izmax0))

signx = 1.0_num
signy = 1.0_num

sx=0.0_num
sy=0.0_num
sz=0.0_num
sx0=0.0_num
sy0=0.0_num
sz0=0.0_num
!$OMP PARALLEL DO PRIVATE(ip,ll,jj,kk,x,y,z,j,k,l,j0,k0,l0,xint,yint,zint,sx,sy,sz,sx0,sy0, & 
!$OMP sz0,oxint,xintsq,oxintsq,oyint,yintsq,oyintsq, ozint,zintsq,ozintsq)
DO ip=1,np
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi

    IF (l_lower_order_in_v) THEN
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x-0.5_num)
        ELSE
            j=floor(x)
            j0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y-0.5_num)
        ELSE
            k=floor(y)
            k0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z-0.5_num)
        ELSE
            l=floor(z)
            l0=floor(z)
        END IF
    ELSE
        IF (nox==2*(nox/2)) THEN
            j=nint(x)
            j0=floor(x)
        ELSE
            j=floor(x)
            j0=floor(x-0.5_num)
        END IF
        IF (noy==2*(noy/2)) THEN
            k=nint(y)
            k0=floor(y)
        ELSE
            k=floor(y)
            k0=floor(y-0.5_num)
        END IF
        IF (noz==2*(noz/2)) THEN
            l=nint(z)
            l0=floor(z)
        ELSE
            l=floor(z)
            l0=floor(z-0.5_num)
        END IF
    END IF

    xint=x-j
    yint=y-k
    zint=z-l

    IF (nox==1) THEN
        sx( 0) = 1.0_num-xint
        sx( 1) = xint
    ELSEIF (nox==2) THEN
        xintsq = xint*xint
        sx(-1) = 0.5_num*(0.5_num-xint)**2
        sx( 0) = 0.75_num-xintsq
        sx( 1) = 0.5_num*(0.5_num+xint)**2
    ELSEIF (nox==3) THEN
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx( 2) = onesixth*xintsq*xint
    END IF

    IF (noy==1) THEN
        sy( 0) = 1.0_num-yint
        sy( 1) = yint
    ELSEIF (noy==2) THEN
        yintsq = yint*yint
        sy(-1) = 0.5_num*(0.5_num-yint)**2
        sy( 0) = 0.75_num-yintsq
        sy( 1) = 0.5_num*(0.5_num+yint)**2
    ELSEIF (noy==3) THEN
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy( 2) = onesixth*yintsq*yint
    END IF

    IF (noz==1) THEN
        sz( 0) = 1.0_num-zint
        sz( 1) = zint
    ELSEIF (noz==2) THEN
        zintsq = zint*zint
        sz(-1) = 0.5_num*(0.5_num-zint)**2
        sz( 0) = 0.75_num-zintsq
        sz( 1) = 0.5_num*(0.5_num+zint)**2
    ELSEIF (noz==3) THEN
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz( 2) = onesixth*zintsq*zint
    END IF

    xint=x-0.5_num-j0
    yint=y-0.5_num-k0
    zint=z-0.5_num-l0

    IF (l_lower_order_in_v) THEN
        IF (nox==1) THEN
            sx0( 0) = 1.0_num
        ELSEIF (nox==2) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==3) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num
        ELSEIF (noy==2) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==3) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        END IF

        IF (noz==1) THEN
            sz0( 0) = 1.0_num
        ELSEIF (noz==2) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==3) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        END IF
    ELSE

        IF (nox==1) THEN
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        ELSEIF (nox==2) THEN
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        ELSEIF (nox==3) THEN
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END IF

        IF (noy==1) THEN
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        ELSEIF (noy==2) THEN
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        ELSEIF (noy==3) THEN
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END IF
 
        IF (noz==1) THEN
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        ELSEIF (noz==2) THEN
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        ELSEIF (noz==3) THEN
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END IF
    END IF

    DO ll = izmin0, izmax0
        DO kk = iymin0, iymax0
            DO jj = ixmin, ixmax+1
                bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k0+kk,l0+ll)*signx
            END DO
        END DO
    END DO

    DO ll = izmin0, izmax0
        DO kk = iymin, iymax+1
            DO jj = ixmin0, ixmax0
                by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj,k+kk,l0+ll)*signy
            END DO
        END DO
    END DO

    DO ll = izmin, izmax+1
        DO kk = iymin0, iymax0
            DO jj = ixmin0, ixmax0
                bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj,k0+kk,l+ll)
            END DO
        END DO
    END DO
END DO
!OMP END PARALLEL DO
DEALLOCATE(sx0,sz0)

RETURN
END SUBROUTINE getb3d_n_energy_conserving

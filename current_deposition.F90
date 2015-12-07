!===============================================================================
! Deposit current in each tile
!===============================================================================
SUBROUTINE depose_currents_on_grid_jxjyjz
USE particles
USE constants
USE fields
USE params
USE shared_data
USE tiling
USE omp_lib
IMPLICIT NONE
INTEGER :: ispecies, ix, iy, iz, count
INTEGER :: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
REAL(num) :: tdeb, tend
REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jx_tile,jy_tile,jz_tile
INTEGER :: nxc, nyc, nzc

jx = 0.0_num
jy = 0.0_num
jz = 0.0_num

tdeb=MPI_WTIME()
!!! ------------------------------------------------------------------------------------------------------
!!! --- Adding currents from tiles to global arrays (THIS ALGO AVOIDS REDUCTION OPERATION)
!!! ------------------------------------------------------------------------------------------------------

!! **** STEP 1:  Each tile adds current contribution to their local cells (not guardcells)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguards,nyjguards,nzjguards,dx,dy,dz,dt,jx,jy,jz,nox,noy,noz) &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,count,jmin,jmax,kmin,kmax,lmin, &
!$OMP lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,jx_tile,jy_tile,jz_tile,nxc,nyc,nzc)
!! Current deposition
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                curr_tile%jxtile = 0.0_num; curr_tile%jytile = 0.0_num
                curr_tile%jztile = 0.0_num
                ! Depose current in jtile
                CALL depose_jxjyjz_esirkepov_n(curr_tile%jxtile,curr_tile%jytile,curr_tile%jztile,count,        &
                curr_tile%part_x(1:count),curr_tile%part_y(1:count),curr_tile%part_z(1:count),                  &
                curr_tile%part_ux(1:count),curr_tile%part_uy(1:count),curr_tile%part_uz(1:count),               &
                curr_tile%weight(1:count),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,      &
                curr_tile%z_grid_tile_min,dt,dx,dy,dz,curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,          &
                curr_tile%nz_cells_tile,nxjguards,nyjguards,nzjguards,nox,noy,noz,.TRUE.,.FALSE.)
                ! Reduce jtile in j
                jx(jmin:jmax,kmin:kmax,lmin:lmax) = jx(jmin:jmax,kmin:kmax,lmin:lmax) + curr_tile%jxtile(0:nxc,0:nyc,0:nzc)
                jy(jmin:jmax,kmin:kmax,lmin:lmax) = jy(jmin:jmax,kmin:kmax,lmin:lmax) + curr_tile%jytile(0:nxc,0:nyc,0:nzc)
                jz(jmin:jmax,kmin:kmax,lmin:lmax) = jz(jmin:jmax,kmin:kmax,lmin:lmax) + curr_tile%jztile(0:nxc,0:nyc,0:nzc)
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO

!! **** STEP 2:  Each tile adds charge contribution to adjacent tiles
!! **** This step requires sync for each dimension to avoid thread contention
!----  +/- X direction
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                jminc=jmin-nxjguards; jmaxc=jmax+nxjguards
                kminc=kmin-nyjguards; kmaxc=kmax+nyjguards
                lminc=lmin-nzjguards; lmaxc=lmax+nzjguards
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! - FACES +/- X
                ! --- JX
                jx(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jx(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                curr_tile%jxtile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                jx(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jx(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                curr_tile%jxtile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                ! --- JY
                jy(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jy(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                curr_tile%jytile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                jy(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jy(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                curr_tile%jytile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                ! --- JZ
                jz(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jz(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                curr_tile%jztile(-nxjguards:-1,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                jz(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jz(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                curr_tile%jztile(nxc+1:nxc+nxjguards,-nyjguards:nyc+nyjguards,-nzjguards:nzc+nzjguards)
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!----  +/- Y direction
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                jminc=jmin-nxjguards; jmaxc=jmax+nxjguards
                kminc=kmin-nyjguards; kmaxc=kmax+nyjguards
                lminc=lmin-nzjguards; lmaxc=lmax+nzjguards
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! - FACES +/- Y
                ! --- JX
                jx(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jx(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                curr_tile%jxtile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
                jx(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jx(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                curr_tile%jxtile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                ! --- JY
                jy(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jy(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                curr_tile%jytile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
                jy(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jy(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                curr_tile%jytile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)
                ! --- JZ
                jz(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jz(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                curr_tile%jztile(0:nxc,-nyjguards:-1,-nzjguards:nzc+nzjguards)
                jz(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jz(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                curr_tile%jztile(0:nxc,nyc+1:nyc+nyjguards,-nzjguards:nzc+nzjguards)
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!----  +/- Z direction
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                jminc=jmin-nxjguards; jmaxc=jmax+nxjguards
                kminc=kmin-nyjguards; kmaxc=kmax+nyjguards
                lminc=lmin-nzjguards; lmaxc=lmax+nzjguards
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! - FACES +/- Z
                ! --- JX
                jx(jmin:jmax,kmin:kmax,lminc:lmin-1) = jx(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                curr_tile%jxtile(0:nxc, 0:nyc,-nzjguards:-1)
                jx(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jx(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                curr_tile%jxtile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)
                ! --- JY
                jy(jmin:jmax,kmin:kmax,lminc:lmin-1) = jy(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                curr_tile%jytile(0:nxc, 0:nyc,-nzjguards:-1)
                jy(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jy(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                curr_tile%jytile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)
                ! --- JZ
                jz(jmin:jmax,kmin:kmax,lminc:lmin-1) = jz(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                curr_tile%jztile(0:nxc, 0:nyc,-nzjguards:-1)
                jz(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jz(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                curr_tile%jztile(0:nxc, 0:nyc,nzc+1:nzc+nzjguards)
            END DO! END LOOP ON SPECIES
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!$OMP END PARALLEL
tend=MPI_WTIME()
pushtime=pushtime+(tend-tdeb)

END SUBROUTINE depose_currents_on_grid_jxjyjz


!!! --- Order 1 3D scalar current deposition routine (rho*v)
!!! This version does not vectorize on SIMD architectures
SUBROUTINE depose_jxjyjz_scalar_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    IMPLICIT NONE
    INTEGER :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,vx,vy,vz,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) :: wq, wqx, wqy, wqz, gaminv, usq, clightsq
    REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER :: j,k,l,ip
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

    ! LOOP ON PARTICLES
    DO ip=1,np
        ! --- computes position in  grid units at (n+1)
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        ! Computes velocity
        usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clightsq
        gaminv = 1.0_num/sqrt(1.0_num + usq)
        vx = uxp(ip)*gaminv
        vy = uyp(ip)*gaminv
        vz = uzp(ip)*gaminv

        ! --- computes particles weights
        wq=q*w(ip)
        wqx=wq*invvol*vx
        wqy=wq*invvol*vy
        wqz=wq*invvol*vz

        ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
        xmid=x-dts2dx*vx
        ymid=y-dts2dy*vy
        zmid=z-dts2dz*vz

        ! --- finds node of cell containing particles for current positions
        j=floor(xmid)
        k=floor(ymid)
        l=floor(zmid)

        ! --- computes set of coefficients for node centered quantities
        xint = xmid-j
        yint = ymid-k
        zint = zmid-l
        sx( 0) = 1.0_num-xint
        sx( 1) = xint
        sy( 0) = 1.0_num-yint
        sy( 1) = yint
        sz( 0) = 1.0_num-zint
        sz( 1) = zint

        ! --- computes set of coefficients for staggered quantities
        xint = xmid-j-0.5_num
        yint = ymid-k-0.5_num
        zint = zmid-l-0.5_num
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint

        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! - JX
        jx(j,k,l)      = jx(j,k,l)+sx0(0)*sy(0)*sz(0)*wqx
        jx(j+1,k,l)    = jx(j+1,k,l)+sx0(1)*sy(0)*sz(0)*wqx
        jx(j,k+1,l)    = jx(j,k+1,l)+sx0(0)*sy(1)*sz(0)*wqx
        jx(j+1,k+1,l)  = jx(j+1,k+1,l)+sx0(1)*sy(1)*sz(0)*wqx
        jx(j,k,l+1)    = jx(j,k,l+1)+sx0(0)*sy(0)*sz(1)*wqx
        jx(j+1,k,l+1)  = jx(j+1,k,l+1)+sx0(1)*sy(0)*sz(1)*wqx
        jx(j,k+1,l+1)  = jx(j,k+1,l+1)+sx0(0)*sy(1)*sz(1)*wqx
        jx(j+1,k+1,l+1)= jx(j+1,k+1,l+1)+sx0(1)*sy(1)*sz(1)*wqx

        ! - JY
        jy(j,k,l)      = jy(j,k,l)+sx(0)*sy0(0)*sz(0)*wqy
        jy(j+1,k,l)    = jy(j+1,k,l)+sx(1)*sy0(0)*sz(0)*wqy
        jy(j,k+1,l)    = jy(j,k+1,l)+sx(0)*sy0(1)*sz(0)*wqy
        jy(j+1,k+1,l)  = jy(j+1,k+1,l)+sx(1)*sy0(1)*sz(0)*wqy
        jy(j,k,l+1)    = jy(j,k,l+1)+sx(0)*sy0(0)*sz(1)*wqy
        jy(j+1,k,l+1)  = jy(j+1,k,l+1)+sx(1)*sy0(0)*sz(1)*wqy
        jy(j,k+1,l+1)  = jy(j,k+1,l+1)+sx(0)*sy0(1)*sz(1)*wqy
        jy(j+1,k+1,l+1)= jy(j+1,k+1,l+1)+sx(1)*sy0(1)*sz(1)*wqy

        ! - JZ
        jz(j,k,l)      = jz(j,k,l)+sx(0)*sy(0)*sz0(0)*wqz
        jz(j+1,k,l)    = jz(j+1,k,l)+sx(1)*sy(0)*sz0(0)*wqz
        jz(j,k+1,l)    = jz(j,k+1,l)+sx(0)*sy(1)*sz0(0)*wqz
        jz(j+1,k+1,l)  = jz(j+1,k+1,l)+sx(1)*sy(1)*sz0(0)*wqz
        jz(j,k,l+1)    = jz(j,k,l+1)+sx(0)*sy(0)*sz0(1)*wqz
        jz(j+1,k,l+1)  = jz(j+1,k,l+1)+sx(1)*sy(0)*sz0(1)*wqz
        jz(j,k+1,l+1)  = jz(j,k+1,l+1)+sx(0)*sy(1)*sz0(1)*wqz
        jz(j+1,k+1,l+1)= jz(j+1,k+1,l+1)+sx(1)*sy(1)*sz0(1)*wqz
    END DO
    RETURN
END SUBROUTINE depose_jxjyjz_scalar_1_1_1


!===========================================================================================
! Esirkepov current deposition algorithm at order 1 in x, y, z (nox=noy=noz=1)
SUBROUTINE depose_jxjyjz_esirkepov_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin, &
                                      dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
!===========================================================================================
USE omp_lib
USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
REAL(num), DIMENSION(:,:,:), ALLOCATABLE:: jx1, jy1, jz1
REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w
REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
REAL(num) :: clghtisq,usq,gaminv,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz, &
                                      s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,         &
                                      oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                                      dtsdx0,dtsdy0,dtsdz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx
REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
INTEGER :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                                      ixmin, ixmax, iymin, iymax, izmin, izmax

! PARAMETER INIT
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
ALLOCATE(sdx(-1:2,-1:2,-1:2),sdy(-1:2,-1:2,-1:2),sdz(-1:2,-1:2,-1:2))
ALLOCATE(sx(-1:2), sx0(-1:2), dsx(-1:2))
ALLOCATE(sy(-1:2), sy0(-1:2), dsy(-1:2))
ALLOCATE(sz(-1:2), sz0(-1:2), dsz(-1:2))
clghtisq = 1.0_num/clight**2
dtsdz0 = dt*dzi
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num
ALLOCATE(jx1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), &
         jy1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), &
         jz1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
!!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ip,x,y,z,usq,vx,vy,vz,gaminv,xold,yold,zold, &
!!$OMP wq,wqx,wqy,wqz,iixp0,ijxp0,ikxp0, xint,yint,zint, oxint,xintsq, oxintsq,dix,diy,diz, &
!!$OMP dsx, dsy, dsz, oyint,yintsq, oyintsq, ozint,zintsq, ozintsq,ixmin, ixmax, iymin, iymax, izmin, izmax,  &
!!$OMP k,j,i,kc,jc,ic, iixp, ijxp, ikxp,sx,sy,sz, sx0,sy0,sz0,sdx,sdy,sdz,jx1,jy1,jz1) &
!!$OMP SHARED(np,xp,yp,zp,uxp,uyp,uzp,w,dxi,dyi,dzi,invdtdx,invdtdy,invdtdz,xmin,ymin,zmin,clghtisq,dtsdx0,dtsdy0,dtsdz0,q,jx,jy,jz)
jx1=0.0_num; jy1=0.0_num; jz1=0.0_num
!!$OMP DO
DO ip=1,np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq
    gaminv = 1.0_num/sqrt(1.0_num + usq)
    vx = uxp(ip)*gaminv
    vy = uyp(ip)*gaminv
    vz = uzp(ip)*gaminv
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqy = wq*invdtdy
    wqz = wq*invdtdz
    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ijxp0=floor(y)
    ikxp0=floor(z)
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    yint=y-ijxp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
    sy0( 0) = 1.0_num-yint
    sy0( 1) = yint
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint
    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ijxp=floor(yold)
    ikxp=floor(zold)
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    yint = yold-ijxp
    zint = zold-ikxp
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diy = ijxp-ijxp0
    diz = ikxp-ikxp0
    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx=0.0_num;sy=0.0_num;sz=0.0_num
    ! --- computes coefficients for quantities centered between nodes
    sx( 0+dix) = 1.0_num-xint
    sx( 1+dix) = xint
    sy( 0+diy) = 1.0_num-yint
    sy( 1+diy) = yint
    sz( 0+diz) = 1.0_num-zint
    sz( 1+diz) = zint
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsy = sy - sy0
    dsz = sz - sz0
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)
    ixmax = max(0,dix)+1
    iymin = min(0,diy)
    iymax = max(0,diy)+1
    izmin = min(0,diz)
    izmax = max(0,diz)+1
    
    ! --- add current contributions
    DO k=izmin, izmax
      DO j=iymin, iymax
          DO i=ixmin, ixmax
              ic = iixp0+i
              jc = ijxp0+j
              kc = ikxp0+k
              IF(i<ixmax) THEN
                  sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
                  (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
                  IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                  jx1(ic,jc,kc) = jx1(ic,jc,kc) + sdx(i,j,k)
              END IF
              IF(j<iymax) THEN
                  sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
                  (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
                  IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                  jy1(ic,jc,kc) = jy1(ic,jc,kc) + sdy(i,j,k)
              END IF
              IF(k<izmax) THEN
                  sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
                  (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
                  IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                  jz1(ic,jc,kc) = jz1(ic,jc,kc) + sdz(i,j,k)
              END IF
          END DO
      END DO
    END DO
END DO
!!$OMP END DO
!!$OMP CRITICAL
jx=jx+jx1
jy=jy+jy1
jz=jz+jz1
!!$OMP END CRITICAL
!!$OMP END PARALLEL
DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz,jx1,jy1,jz1)
RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_1_1_1

!===========================================================================================
! ! Esirkepov current deposition algorithm for linear, quadratic or cubic splines
! WARNING: Highly unoptimized routine ---> USE INLINED ROUTINE
SUBROUTINE depose_jxjyjz_esirkepov_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin, &
dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
nox,noy,noz,l_particles_weight,l4symtry)
!===========================================================================================

USE constants
IMPLICIT NONE
INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
REAL(num), DIMENSION(:,:,:), ALLOCATABLE:: jx1, jy1, jz1
REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w
REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
REAL(num) :: clghtisq,usq,gaminv,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
REAL(num), DIMENSION(:), ALLOCATABLE :: sx, sx0, dsx
REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
INTEGER :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx, ndtody, ndtodz, &
xl,xu,yl,yu,zl,zu
LOGICAL :: l_particles_weight,l4symtry

! PARAMETER INIT
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
ALLOCATE(sdx(xl:xu,yl:yu,zl:zu),sdy(xl:xu,yl:yu,zl:zu),sdz(xl:xu,yl:yu,zl:zu))
ALLOCATE(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
ALLOCATE(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
ALLOCATE(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))
ALLOCATE(jx1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), &
         jy1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), &
         jz1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
clghtisq = 1.0_num/clight**2
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num
jx1=0.0_num;jy1=0.0_num;jz1=0.0_num

!$OMP PARALLEL PRIVATE(ip,x,y,z,usq,vx,vy,vz,gaminv,xold,yold,zold,ncells,dtsdx,dtsdy,dtsdz,dts2dx,dts2dy,dts2dz, &
!$OMP icell, wq,wqx,wqy,wqz,iixp0,ijxp0,ikxp0, xint,yint,zint, oxint,xintsq, oxintsq,dix,diy,diz, &
!$OMP dsx, dsy, dsz, oyint,yintsq, oyintsq, ozint,zintsq, ozintsq,ixmin, ixmax, iymin, iymax, izmin, izmax,  &
!$OMP k,j,i,kc,jc,ic, iixp, ijxp, ikxp,sx,sy,sz) FIRSTPRIVATE(sx0,sy0,sz0,sdx,sdy,sdz,jx1,jy1,jz1)
!$OMP DO
DO ip=1,np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq
    gaminv = 1.0_num/sqrt(1.0_num + usq)
    vx = uxp(ip)*gaminv
    vy = uyp(ip)*gaminv
    vz = uzp(ip)*gaminv
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- applies 4-fold symmetry
    IF (l4symtry) THEN
        x=abs(x)
        y=abs(y)
        xold=abs(xold)
        yold=abs(yold)
        vx = (x-xold)/dtsdx0
        vy = (y-yold)/dtsdy0
    END IF
    !computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdy = dtsdy0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dy = dts2dy0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    y=yold
    z=zold

    DO icell = 1,ncells
        xold = x
        yold = y
        zold = z
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz

        ! --- computes particles "weights"
        IF (l_particles_weight) THEN
            wq=q*w(ip)
        ELSE
            wq=q*w(1)
        END IF
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

! --- finds node of cell containing particles for current positions
! --- (different for odd/even spline orders)
        IF (nox==2*(nox/2)) THEN
            iixp0=nint(x)
        ELSE
            iixp0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            ijxp0=nint(y)
        ELSE
            ijxp0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            ikxp0=nint(z)
        ELSE
            ikxp0=floor(z)
        END IF
! --- computes distance between particle and node for current positions
        xint=x-iixp0
        yint=y-ijxp0
        zint=z-ikxp0

! --- computes coefficients for node centered quantities
        SELECT CASE(nox)
        CASE(0)
            sx0( 0) = 1.0_num
        CASE(1)
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        CASE(2)
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        CASE(3)
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END SELECT

        SELECT CASE(noy)
        CASE(0)
            sy0( 0) = 1.0_num
        CASE(1)
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        CASE(2)
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        CASE(3)
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END SELECT

        SELECT CASE(noz)
        CASE(0)
            sz0( 0) = 1.0_num
        CASE(1)
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        CASE(2)
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        CASE(3)
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END SELECT

! --- finds node of cell containing particles for old positions
! --- (different for odd/even spline orders)
        IF (nox==2*(nox/2)) THEN
            iixp=nint(xold)
        ELSE
            iixp=floor(xold)
        END IF
        IF (noy==2*(noy/2)) THEN
            ijxp=nint(yold)
        ELSE
            ijxp=floor(yold)
        END IF
        IF (noz==2*(noz/2)) THEN
            ikxp=nint(zold)
        ELSE
            ikxp=floor(zold)
        END IF

! --- computes distance between particle and node for old positions
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.0_num;sy=0.0_num;sz=0.0_num

! --- computes coefficients for quantities centered between nodes
        SELECT CASE(nox)
        CASE(0)
            sx( 0+dix) = 1.0_num
        CASE(1)
            sx( 0+dix) = 1.0_num-xint
            sx( 1+dix) = xint
        CASE(2)
            xintsq = xint*xint
            sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
            sx( 0+dix) = 0.75_num-xintsq
            sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
        CASE(3)
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx(-1+dix) = onesixth*oxintsq*oxint
            sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx( 2+dix) = onesixth*xintsq*xint
        END SELECT

        SELECT CASE(noy)
        CASE(0)
            sy( 0+diy) = 1.0_num
        CASE(1)
            sy( 0+diy) = 1.0_num-yint
            sy( 1+diy) = yint
        CASE(2)
            yintsq = yint*yint
            sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
            sy( 0+diy) = 0.75_num-yintsq
            sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
        CASE(3)
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy(-1+diy) = onesixth*oyintsq*oyint
            sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy( 2+diy) = onesixth*yintsq*yint
        END SELECT

        SELECT CASE(noz)
        CASE(0)
            sz( 0+diz) = 1.0_num
        CASE(1)
            sz( 0+diz) = 1.0_num-zint
            sz( 1+diz) = zint
        CASE(2)
            zintsq = zint*zint
            sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
            sz( 0+diz) = 0.75_num-zintsq
            sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
        CASE(3)
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz(-1+diz) = onesixth*ozintsq*ozint
            sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz( 2+diz) = onesixth*zintsq*zint
        END SELECT

! --- computes coefficients difference
        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0

! --- computes min/max positions of current contributions
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)

! --- add current contributions
        DO k=izmin, izmax
            DO j=iymin, iymax
                DO i=ixmin, ixmax
                    ic = iixp0+i
                    jc = ijxp0+j
                    kc = ikxp0+k
                    IF(i<ixmax) THEN
                        sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
                        (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
                        IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                        jx1(ic,jc,kc) = jx1(ic,jc,kc) + sdx(i,j,k)
                    END IF
                    IF(j<iymax) THEN
                        sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
                        (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
                        IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                        jy1(ic,jc,kc) = jy1(ic,jc,kc) + sdy(i,j,k)
                    END IF
                    IF(k<izmax) THEN
                        sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
                        (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
                        IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                        jz1(ic,jc,kc) = jz1(ic,jc,kc) + sdz(i,j,k)
                    END IF
                END DO
            END DO
        END DO
    END DO
END DO
!$OMP END DO

!$OMP CRITICAL
jx=jx+jx1
jy=jy+jy1
jz=jz+jz1
!$OMP END CRITICAL
!$OMP END PARALLEL
DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz,jx1,jy1,jz1)

RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_n

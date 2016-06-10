MODULE diagnostics
!!!! --- This module contains useful diagnostics to test code correctness
    USE constants
    IMPLICIT NONE


CONTAINS

    !!! --- Computes derived physical quantities from simulation
    SUBROUTINE calc_diags
        USE fields
        USE boundary
        USE particles
        USE params
        USE shared_data
        USE tiling
        USE time_stat
        IMPLICIT NONE

        REAL(num) :: tmptime
        
#if defined(DEBUG)
        WRITE(0,*) "Calc_diags: start"
#endif
        
        tmptime = MPI_WTIME()        

        ! - Computes electric field divergence on grid at n+1
        dive=0.0_num
        CALL calc_field_div(dive, ex, ey, ez, nx, ny, nz, nxguards, nyguards, nzguards, dx, dy, dz)

        ! - Computes total charge density
        CALL pxrdepose_rho_on_grid()
        
        CALL charge_bcs()

        localtimes(9) = localtimes(9) + (MPI_WTIME() - tmptime)

#if defined(DEBUG)
        WRITE(0,*) "Calc_diags: stop"
#endif

    END SUBROUTINE calc_diags

!===============================================================================
! Charge deposition algorithms
!===============================================================================
SUBROUTINE pxrdepose_rho_on_grid
USE fields
USE shared_data
USE params
IMPLICIT NONE 
INTEGER(idp) :: c_rho_old
c_rho_old = 0
rho = 0.0_num
! DEPOSIT Charge on the grid 
CALL pxrdepose_rho_on_grid_sub_openmp(rho,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	nox,noy,noz,dx,dy,dz,dt,c_rho_old)

END SUBROUTINE pxrdepose_rho_on_grid

!===============================================================================
! Deposit rho in each tile
! OpenMP version. Avoids conflict while reducing tile charge in the global 
! charge array. 
!===============================================================================
SUBROUTINE pxrdepose_rho_on_grid_sub_openmp(rhog,nxx,nyy,nzz,nxjguard,nyjguard,nzjguard, &
	noxx,noyy,nozz,dxx,dyy,dzz,dtt, c_rho_old)
USE particles
USE constants
USE tiling
USE omp_lib
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz, c_rho_old
REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
REAL(num), INTENT(IN OUT) :: rhog(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
INTEGER(idp) :: ispecies, ix, iy, iz, count
INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
TYPE(grid_tile), POINTER :: currg
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
LOGICAL(idp) :: isdeposited=.FALSE.

!$OMP PARALLEL DEFAULT(NONE)                                                              &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguard,nyjguard,              &
!$OMP nzjguard,dxx,dyy,dzz,dtt,rhog,noxx,noyy,nozz,aofgrid_tiles, c_dim, c_rho_old)       &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,currg, curr_tile,count,jmin,jmax,kmin,kmax,lmin,     &
!$OMP lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,nxc,nyc,nzc, nxjg, nyjg, nzjg, isdeposited)
!! Current deposition
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
        	curr => species_parray(1)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
            nxjg=curr_tile%nxg_tile
            nyjg=curr_tile%nyg_tile
            nzjg=curr_tile%nzg_tile
            jmin=curr_tile%nx_tile_min
        	jmax=curr_tile%nx_tile_max
            kmin=curr_tile%ny_tile_min
            kmax=curr_tile%ny_tile_max
            lmin=curr_tile%nz_tile_min
            lmax=curr_tile%nz_tile_max
            nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
            nzc=curr_tile%nz_cells_tile         
			currg=>aofgrid_tiles(ix,iy,iz)
            currg%rhotile=0._num
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
           	    curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .EQ. 0) THEN 
                	CYCLE
                ELSE 
                	isdeposited=.TRUE.
                ENDIF 
                ! Depose charge in rhotile
                SELECT CASE (c_dim) 
				CASE (2)
					SELECT CASE (c_rho_old)
					CASE(1) ! Rho at older time 
						CALL pxr_depose_rhoold_n_2dxz(currg%rhotile(:,0,:),count,       &
						curr_tile%part_x,curr_tile%part_z,     							&
						curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,     		&
						curr_tile%part_gaminv,curr_tile%pid(1,wpid),curr%charge,  		&
						curr_tile%x_grid_tile_min,  curr_tile%z_grid_tile_min,  		&
						dtt,dxx,dzz,nxc,nzc,                          					&
						nxjg,nzjg,noxx,nozz,.TRUE._idp,.FALSE._idp)
					CASE DEFAULT  ! Rho at current time 
						CALL pxr_depose_rho_n_2dxz(currg%rhotile(:,0,:),count,              &
						curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     			&
						curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,     	&
						curr_tile%z_grid_tile_min,dxx,dzz,nxc,nzc,                          &
						nxjg,nzjg,noxx,nozz,.TRUE._idp,.FALSE._idp,.FALSE._idp,0_idp)
					END SELECT 
				CASE DEFAULT 
					CALL pxr_depose_rho_n(currg%rhotile,count,                                                 &
					curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						           &
					curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,     &
					curr_tile%z_grid_tile_min,dxx,dyy,dzz,nxc,nyc,nzc,                                         &
					nxjg,nyjg,nzjg,noxx,noyy,nozz,.TRUE.,.FALSE.) 
			    END SELECT 
            END DO! END LOOP ON SPECIES
            IF (isdeposited) THEN
            	rhog(jmin:jmax,kmin:kmax,lmin:lmax)=rhog(jmin:jmax,kmin:kmax,lmin:lmax)+currg%rhotile(0:nxc,0:nyc,0:nzc)
            ENDIF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!! Adding charge from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
!+/- X
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
        	isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .GT. 0) isdeposited=.TRUE.  
            END DO
            IF (isdeposited) THEN 
            	currg=>aofgrid_tiles(ix,iy,iz)
                curr => species_parray(1)
           		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxjg=curr_tile%nxg_tile
                nyjg=curr_tile%nyg_tile
                nzjg=curr_tile%nzg_tile
                jminc=jmin-nxjg; jmaxc=jmax+nxjg
                kminc=kmin-nyjg; kmaxc=kmax+nyjg
                lminc=lmin-nzjg; lmaxc=lmax+nzjg
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! --- RHO
                ! - FACES +/- X
                rhog(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = rhog(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                currg%rhotile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                rhog(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = rhog(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                currg%rhotile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
            ENDIF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!+/- Y
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .GT. 0) isdeposited=.TRUE.  
            END DO
            IF (isdeposited) THEN 
            	currg=>aofgrid_tiles(ix,iy,iz)
                curr => species_parray(1)
           		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxjg=curr_tile%nxg_tile
                nyjg=curr_tile%nyg_tile
                nzjg=curr_tile%nzg_tile
                jminc=jmin-nxjg; jmaxc=jmax+nxjg
                kminc=kmin-nyjg; kmaxc=kmax+nyjg
                lminc=lmin-nzjg; lmaxc=lmax+nzjg
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! --- RHO
                ! - FACES +/- Y
                rhog(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = rhog(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                currg%rhotile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                rhog(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = rhog(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                currg%rhotile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
            END IF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
! +/-Z
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .GT. 0) isdeposited=.TRUE.  
            END DO
            IF (isdeposited) THEN 
            	currg=>aofgrid_tiles(ix,iy,iz)
                curr => species_parray(1)
           		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxjg=curr_tile%nxg_tile
                nyjg=curr_tile%nyg_tile
                nzjg=curr_tile%nzg_tile
                jminc=jmin-nxjg; jmaxc=jmax+nxjg
                kminc=kmin-nyjg; kmaxc=kmax+nyjg
                lminc=lmin-nzjg; lmaxc=lmax+nzjg
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! --- RHO
                ! - FACES +/- Z
                rhog(jmin:jmax,kmin:kmax,lminc:lmin-1) = rhog(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                currg%rhotile(0:nxc, 0:nyc,-nzjg:-1)
                rhog(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = rhog(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                currg%rhotile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
            END IF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE pxrdepose_rho_on_grid_sub_openmp

  !!! --- Order 1 3D scalar charge deposition routine
    !!! This version does not vectorize on SIMD architectures
    SUBROUTINE depose_rho_scalar_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), INTENT(IN OUT) :: rho
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol
        REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        DO ip=1,np
            ! --- computes current position in grid units
            x = (xp(ip)-xmin)*dxi
            y = (yp(ip)-ymin)*dyi
            z = (zp(ip)-zmin)*dzi
            ! --- finds node of cell containing particles for current positions
            j=floor(x)
            k=floor(y)
            l=floor(z)
            ! --- computes distance between particle and node for current positions
            xint = x-j
            yint = y-k
            zint = z-l
            ! --- computes particles weights
            wq=q*w(ip)*invvol
            ! --- computes coefficients for node centered quantities
            sx( 0) = 1.0_num-xint
            sx( 1) = xint
            sy( 0) = 1.0_num-yint
            sy( 1) = yint
            sz( 0) = 1.0_num-zint
            sz( 1) = zint
            ! --- add charge density contributions
            rho(j,k,l)      = rho(j,k,l)+sx(0)*sy(0)*sz(0)*wq
            rho(j+1,k,l)    = rho(j+1,k,l)+sx(1)*sy(0)*sz(0)*wq
            rho(j,k+1,l)    = rho(j,k+1,l)+sx(0)*sy(1)*sz(0)*wq
            rho(j+1,k+1,l)  = rho(j+1,k+1,l)+sx(1)*sy(1)*sz(0)*wq
            rho(j,k,l+1)    = rho(j,k,l+1)+sx(0)*sy(0)*sz(1)*wq
            rho(j+1,k,l+1)  = rho(j+1,k,l+1)+sx(1)*sy(0)*sz(1)*wq
            rho(j,k+1,l+1)  = rho(j,k+1,l+1)+sx(0)*sy(1)*sz(1)*wq
            rho(j+1,k+1,l+1)= rho(j+1,k+1,l+1)+sx(1)*sy(1)*sz(1)*wq
        END DO
        RETURN
    END SUBROUTINE depose_rho_scalar_1_1_1

!!! --- Order 1 3D vector charge deposition routine
!!! --- Computes charge density on grid vectorized with Schwarzmeier and Hewitt scheme
!!! --- This routine does vectorize on SIMD architecture but poor performances
    SUBROUTINE depose_rho_vecSH_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol, sx0, sy0, sz0, sx1, sy1, sz1
        REAL(num), ALLOCATABLE :: ww(:,:)
        INTEGER(idp), ALLOCATABLE :: ll(:,:)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: j,k,l,nn,ip,n,m,ixmin, ixmax, iymin, iymax, izmin, izmax
        INTEGER(idp) :: nblk
        INTEGER(idp) :: nnx, nnxy, ind0
        INTEGER(idp) :: moff(1:8)

        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        ! Vectorization parameters init
        nblk=4 ! multiple of vector instruction length (32 bytes on Avx)
        ALLOCATE(ww(1:8,1:nblk),ll(1:8,1:nblk))

        nnx = nx + 1 + 2*nxguard
        nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
        moff(1) = 0_idp
        moff(2) = 1_idp
        moff(3) = nnx
        moff(4) = nnx+1_idp
        moff(5) = nnxy
        moff(6) = nnxy+1_idp
        moff(7) = nnxy+nnx
        moff(8) = nnxy+nnx+1_idp

        DO ip=1,np,nblk
#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
            !DIR$ ASSUME_ALIGNED w:64,ww:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,ww,xp,yp,zp,w)
#endif 
            !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
            !DIR$ ASSUME_ALIGNED w:64,ww:64
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
            DO n=ip,MIN(ip+nblk-1,np) !!!! vector
                nn=n-ip+1
                !- Computations relative to particle n
            	! --- computes current position in grid units
            	x = (xp(n)-xmin)*dxi
            	y = (yp(n)-ymin)*dyi
            	z = (zp(n)-zmin)*dzi
            	! --- finds node of cell containing particles for current positions
            	j=floor(x)
            	k=floor(y)
            	l=floor(z)
            	! --- computes distance between particle and node for current positions
            	xint = x-j
            	yint = y-k
            	zint = z-l
            	! --- computes particles weights
            	wq=q*w(n)*invvol
                ! --- computes coefficients for node centered quantities
            	sx0 = 1.0_num-xint
            	sx1 = xint
            	sy0 = 1.0_num-yint
            	sy1 = yint
            	sz0 = 1.0_num-zint
            	sz1 = zint
            	! --- computes weight for each of the 8-vertices surrounding particle n
                ind0 = (j+nxguard+1) + (k+nyguard+1)*nnx + (l+nzguard+1)*nnxy
            	ww(1,nn) = sx0*sy0*sz0*wq
                ll(1,nn) = ind0+moff(1)
                ww(2,nn) = sx1*sy0*sz0*wq
                ll(2,nn) = ind0+moff(2)
                ww(3,nn) = sx0*sy1*sz0*wq
                ll(3,nn) = ind0+moff(3)
                ww(4,nn) = sx1*sy1*sz0*wq
                ll(4,nn) = ind0+moff(4)
                ww(5,nn) = sx0*sy0*sz1*wq
                ll(5,nn) = ind0+moff(5)
                ww(6,nn) = sx1*sy0*sz1*wq
                ll(6,nn) = ind0+moff(6)
                ww(7,nn) = sx0*sy1*sz1*wq
                ll(7,nn) = ind0+moff(7)
                ww(8,nn) = sx1*sy1*sz1*wq
                ll(8,nn) = ind0+moff(8)
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            ! --- add charge density contributions
            DO m= 1,MIN(nblk,np-ip+1)
#if defined __INTEL_COMPILER 
                !DIR$ ASSUME_ALIGNED ww:64
#elif defined __IBMBGQ__
                !IBM* ALIGN(64,ww)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO l=1,8  !!!! Vector
            		rho(ll(l,m)) = rho(ll(l,m))+ww(l,m)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	        !$OMP END SIMD
#endif
            END DO
        END DO
        DEALLOCATE(ww,ll)
        RETURN
    END SUBROUTINE depose_rho_vecSH_1_1_1

    !!! --- Order 1 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized with Nishiguchi, Orii and Yabe scheme (NOY)
    !!! --- This routine does vectorize on SIMD architecture but poor performances
    SUBROUTINE depose_rho_vecNOY_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), INTENT(IN OUT) :: rho
        REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: rho1
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol
        REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: j,k,l,vv,n,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
        INTEGER(idp), PARAMETER :: LVEC2=4
        REAL(num), DIMENSION(LVEC2,8) :: ww
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        ALLOCATE(rho1(1:LVEC2,-nxguard:nx+nxguard,-nyguard:ny+nyguard, &
                -nzguard:nz+nzguard))
        rho1=0.0_num
        DO ip=1,np, LVEC2
            !DIR$ ASSUME_ALIGNED xp:32
            !DIR$ ASSUME_ALIGNED yp:32
            !DIR$ ASSUME_ALIGNED zp:32
            !DIR$ ASSUME_ALIGNED w:32
            !DIR$ ASSUME_ALIGNED ww:32
            !DIR$ IVDEP
            DO vv=1, MIN(LVEC2,np-ip+1) !!! Vector
                n=vv+ip-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x = (xp(n)-xmin)*dxi
                y = (yp(n)-ymin)*dyi
                z = (zp(n)-zmin)*dzi
                ! --- finds node of cell containing particles for current positions
                j=floor(x)
                k=floor(y)
                l=floor(z)
                ! --- computes distance between particle and node for current positions
                xint = x-j
                yint = y-k
                zint = z-l
                ! --- computes particles weights
                wq=q*w(n)*invvol
                ! --- computes coefficients for node centered quantities
                sx( 0) = 1.0_num-xint
                sx( 1) = xint
                sy( 0) = 1.0_num-yint
                sy( 1) = yint
                sz( 0) = (1.0_num-zint)*wq
                sz( 1) = zint*wq
                ww(vv,1) = sx(0)*sy(0)*sz(0)
                ww(vv,2) = sx(1)*sy(0)*sz(0)
                ww(vv,3) = sx(0)*sy(1)*sz(0)
                ww(vv,4) = sx(1)*sy(1)*sz(0)
                ww(vv,5) = sx(0)*sy(0)*sz(1)
                ww(vv,6) = sx(1)*sy(0)*sz(1)
                ww(vv,7) = sx(0)*sy(1)*sz(1)
                ww(vv,8) = sx(1)*sy(1)*sz(1)
            END DO
!            j=1;k=1;l=1
!            !DIR$ ASSUME_ALIGNED rho1:32
!            !DIR$ ASSUME_ALIGNED ww:32
!            !DIR$ IVDEP
!            DO vv=1, MIN(LVEC2,np-ip+1) !!! Vector
!                ! --- add charge density contributions
!                rho1(vv,j,k,l)      = rho1(vv,j,k,l)+ww(vv,1)
!                rho1(vv,j+1,k,l)    = rho1(vv,j+1,k,l)+ww(vv,2)
!                rho1(vv,j,k+1,l)    = rho1(vv,j,k+1,l)+ww(vv,3)
!                rho1(vv,j+1,k+1,l)  = rho1(vv,j+1,k+1,l)+ww(vv,4)
!                rho1(vv,j,k,l+1)    = rho1(vv,j,k,l+1)+ww(vv,5)
!                rho1(vv,j+1,k,l+1)  = rho1(vv,j+1,k,l+1)+ww(vv,6)
!                rho1(vv,j,k+1,l+1)  = rho1(vv,j,k+1,l+1)+ww(vv,7)
!                rho1(vv,j+1,k+1,l+1)= rho1(vv,j+1,k+1,l+1)+ww(vv,8)
!            END DO
        END DO

        DO jj=-nxguard,nxguard+nx !!! Vector
            DO vv=1,LVEC2
                rho(jj,:,:)=rho(jj,:,:)+rho1(vv,jj,:,:)
            END DO
        END DO
        DEALLOCATE(rho1)
        RETURN
    END SUBROUTINE depose_rho_vecNOY_1_1_1
    !!! --- Order 1 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized (HV-SCHEME v2)
    !!! --- This routine does vectorize on SIMD architecture but poor performances
    SUBROUTINE depose_rho_vecHV_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num), DIMENSION(:,:), ALLOCATABLE:: rhocells
        INTEGER(idp), PARAMETER :: LVEC2=8
        INTEGER(idp), DIMENSION(LVEC2) :: ICELL
        REAL(num) :: wq
        INTEGER(idp) :: NCELLS
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi
        REAL(num) :: xint,yint,zint
        REAL(num) :: x,y,z,invvol
        REAL(num) :: sx(0:1), sy(0:1), sz(0:1)
        REAL(num) :: ww(1:LVEC2,8)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: ic,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
        INTEGER(idp) :: nnx, nnxy
        INTEGER(idp) :: moff(1:8)

        ! Init parameters
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        NCELLS=(2*nxguard+nx)*(2*nyguard+ny)*(2*nzguard+nz)
        ALLOCATE(rhocells(8,NCELLS))
        rhocells=0.0_num
        nnx = nx + 1 + 2*nxguard
        nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
        moff(1) = 0_idp
        moff(2) = 1_idp
        moff(3) = nnx
        moff(4) = nnx+1_idp
        moff(5) = nnxy
        moff(6) = nnxy+1_idp
        moff(7) = nnxy+nnx
        moff(8) = nnxy+nnx+1_idp

        ! FIRST LOOP: computes cell index of particle and their weight on vertices
        DO ip=1,np,LVEC2
            !DIR$ ASSUME_ALIGNED xp:64
            !DIR$ ASSUME_ALIGNED yp:64
            !DIR$ ASSUME_ALIGNED zp:64
            !DIR$ ASSUME_ALIGNED w:64
            !DIR$ ASSUME_ALIGNED ww:64
            !DIR$ ASSUME_ALIGNED ICELL:64
            !DIR$ IVDEP
            DO n=1,MIN(LVEC2,np-ip+1)
                nn=ip+n-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x= (xp(nn)-xmin)*dxi
                z = (zp(nn)-zmin)*dzi
                y = (yp(nn)-ymin)*dyi
                ! --- finds cell containing particles for current positions
                j=floor(x)
                k=floor(y)
                l=floor(z)
                ICELL(n)=1+j+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)
                ! --- computes distance between particle and node for current positions
                xint = x-j
                yint = y-k
                zint = z-l
                ! --- computes particles weights
                wq=q*w(nn)*invvol
                ! --- Compute weight
                sx( 0) = 1.0_num-xint
                sx( 1) = xint
                sy( 0) = 1.0_num-yint
                sy( 1) = yint
                sz( 0) = (1.0_num-zint)*wq
                sz( 1) = zint*wq
                ww(n,1) = sx(0)*sy(0)*sz(0)
                ww(n,2) = sx(1)*sy(0)*sz(0)
                ww(n,3) = sx(0)*sy(1)*sz(0)
                ww(n,4) = sx(1)*sy(1)*sz(0)
                ww(n,5) = sx(0)*sy(0)*sz(1)
                ww(n,6) = sx(1)*sy(0)*sz(1)
                ww(n,7) = sx(0)*sy(1)*sz(1)
                ww(n,8) = sx(1)*sy(1)*sz(1)
            END DO
            ! Current deposition on vertices
            DO n=1,MIN(LVEC2,np-ip+1)
                ! --- add charge density contributions to vertices of the current cell
                ic=ICELL(n)
                !DIR$ ASSUME_ALIGNED rhocells:64
                !DIR$ ASSUME_ALIGNED ww:64
                !!DIR$ NOUNROLL
                DO nv=1,8 !!! - VECTOR
                    rhocells(nv,ic)=rhocells(nv,ic)+ww(n,nv)
                END DO
            END DO
        END DO
        DO nv=1,8
            !DIR$ ASSUME_ALIGNED rhocells:64
            !DIR$ ASSUME_ALIGNED rho:64
            !DIR$ IVDEP
            DO ic=1,NCELLS  !!! VECTOR
                rho(ic+moff(nv))=rho(ic+moff(nv))+rhocells(nv,ic)
            END DO
        END DO
        DEALLOCATE(rhocells)
        RETURN
    END SUBROUTINE depose_rho_vecHV_1_1_1

    !!! --- Order 1 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized (HV-SCHEME v2)
    !!! --- This routine does vectorize on SIMD architecture but poor performances
    SUBROUTINE depose_rho_vecHVv2_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp), INTENT (IN) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num), DIMENSION(:,:), ALLOCATABLE:: rhocells
        INTEGER(idp), PARAMETER :: LVEC2=64
        INTEGER(idp), DIMENSION(LVEC2) :: ICELL
        REAL(num) :: ww
        INTEGER(idp) :: NCELLS
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi
        REAL(num) :: xint,yint,zint
        REAL(num) :: x,y,z,invvol
        REAL(num) :: sx(LVEC2), sy(LVEC2), sz(LVEC2), wq(LVEC2)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: ic,igrid,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
        INTEGER(idp) :: nnx, nnxy
        INTEGER(idp) :: moff(1:8)
        REAL(num):: mx(1:8),my(1:8),mz(1:8), sgn(1:8)
        INTEGER(idp) :: orig, jorig, korig, lorig
        INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy

        ! Init parameters
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard;
        ncx=nx+2;ncy=ny+2;ncz=nz+2
        NCELLS=ncx*ncy*ncz
        ALLOCATE(rhocells(8,NCELLS))
        rhocells=0.0_num
        nnx = ngridx
        nnxy = nnx*ngridy
        moff = (/0_idp,1_idp,nnx,nnx+1_idp,nnxy,nnxy+1_idp,nnxy+nnx,nnxy+nnx+1_idp/)
        mx=(/1_num,0_num,1_num,0_num,1_num,0_num,1_num,0_num/)
        my=(/1_num,1_num,0_num,0_num,1_num,1_num,0_num,0_num/)
        mz=(/1_num,1_num,1_num,1_num,0_num,0_num,0_num,0_num/)
        sgn=(/-1_num,1_num,1_num,-1_num,1_num,-1_num,-1_num,1_num/)
        jorig=-1; korig=-1;lorig=-1
        orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
        ngx=(ngridx-ncx)
        ngxy=(ngridx*ngridy-ncx*ncy)
        ncxy=ncx*ncy
        ! FIRST LOOP: computes cell index of particle and their weight on vertices
        DO ip=1,np,LVEC2
#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED xp:64
            !DIR$ ASSUME_ALIGNED yp:64
            !DIR$ ASSUME_ALIGNED zp:64
            !DIR$ ASSUME_ALIGNED w:64
            !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,xp,yp,zp,w,ICELL)
#endif
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif
            DO n=1,MIN(LVEC2,np-ip+1)
                nn=ip+n-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x= (xp(nn)-xmin)*dxi
                y = (yp(nn)-ymin)*dyi
                z = (zp(nn)-zmin)*dzi
                ! --- finds cell containing particles for current positions
                j=floor(x)
                k=floor(y)
                l=floor(z)
                ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
                ! --- computes distance between particle and node for current positions
                sx(n) = x-j
                sy(n) = y-k
                sz(n) = z-l
                ! --- computes particles weights
                wq(n)=q*w(nn)*invvol
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            ! Current deposition on vertices
            DO n=1,MIN(LVEC2,np-ip+1)
                ! --- add charge density contributions to vertices of the current cell
                ic=ICELL(n)
#if defined __INTEL_COMPILER 
                ic=ICELL(n)
                !DIR$ ASSUME_ALIGNED rhocells:64
                !DIR$ ASSUME_ALIGNED sx:64
                !DIR$ ASSUME_ALIGNED sy:64
                !DIR$ ASSUME_ALIGNED sz:64
#elif defined __IBMBGQ__
                !IBM* ALIGN(64,rhocells,sx,sy,sz)
#endif
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif
                DO nv=1,8 !!! - VECTOR
                    ww=(-mx(nv)+sx(n))*(-my(nv)+sy(n))* &
                        (-mz(nv)+sz(n))*wq(n)*sgn(nv)
                    rhocells(nv,ic)=rhocells(nv,ic)+ww
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            END DO
        END DO
        ! - reduction of rhocells in rho
        DO iz=1, ncz
            DO iy=1,ncy
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif
                DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                    ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                    igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                    rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1,ic)
                    rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2,ic)
                    rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3,ic)
                    rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4,ic)
                    rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5,ic)
                    rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6,ic)
                    rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7,ic)
                    rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8,ic)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	    !$OMP END SIMD
#endif
            END DO
        END DO
        DEALLOCATE(rhocells)
        RETURN
    END SUBROUTINE depose_rho_vecHVv2_1_1_1

    !!! --- Order 2 3D scalar charge deposition routine
    !!! This version does not vectorize on SIMD architectures
    SUBROUTINE depose_rho_scalar_2_2_2(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), INTENT(IN OUT) :: rho
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol,sx1,sx2,sx3,sx4,sx5,sx6,sx7,sx8,sx9
        REAL(num) :: sx(-1:1), sy(-1:1), sz(-1:1)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        DO ip=1,np
            ! --- computes current position in grid units
            x = (xp(ip)-xmin)*dxi
            y = (yp(ip)-ymin)*dyi
            z = (zp(ip)-zmin)*dzi
            ! --- finds node of cell containing particles for current positions
            j=nint(x)
            k=nint(y)
            l=nint(z)
            ! --- computes distance between particle and node for current positions
            xint = x-j
            yint = y-k
            zint = z-l
            ! --- computes particles weights
            wq=q*w(ip)*invvol
            ! --- computes coefficients for node centered quantities
            xintsq = xint*xint
            sx(-1) = 0.5_num*(0.5_num-xint)**2
            sx( 0) = 0.75_num-xintsq
            sx( 1) = 0.5_num*(0.5_num+xint)**2
            yintsq = yint*yint
            sy(-1) = 0.5_num*(0.5_num-yint)**2
            sy( 0) = 0.75_num-yintsq
            sy( 1) = 0.5_num*(0.5_num+yint)**2
            zintsq = zint*zint
            sz(-1) = 0.5_num*(0.5_num-zint)**2*wq
            sz( 0) = (0.75_num-zintsq)*wq
            sz( 1) = 0.5_num*(0.5_num+zint)**2*wq
            sx1=sx(-1)*sy(-1)
            sx2=sx(0)*sy(-1)
            sx3=sx(1)*sy(-1)
            sx4=sx(-1)*sy(0)
            sx5=sx(0)*sy(0)
            sx6=sx(1)*sy(0)
            sx7=sx(-1)*sy(1)
            sx8=sx(0)*sy(1)
            sx9=sx(1)*sy(1)
 	    ! --- add charge density contributions to the 27 
	    ! --- nearest vertices
            rho(j-1,k-1,l-1) = rho(j-1,k-1,l-1)   + sx1*sz(-1)
            rho(j,k-1,l-1)   = rho(j,k-1,l-1)     + sx2*sz(-1)
            rho(j+1,k-1,l-1) = rho(j+1,k-1,l-1)   + sx3*sz(-1)
            rho(j-1,k,l-1)   = rho(j-1,k,l-1)     + sx4*sz(-1)
            rho(j,k,l-1)     = rho(j,k,l-1)       + sx5*sz(-1)
            rho(j+1,k,l-1)   = rho(j+1,k,l-1)     + sx6*sz(-1)
            rho(j-1,k+1,l-1) = rho(j-1,k+1,l-1)   + sx7*sz(-1)
            rho(j,k+1,l-1)   = rho(j,k+1,l-1)     + sx8*sz(-1)
            rho(j+1,k+1,l-1) = rho(j+1,k+1,l-1)   + sx9*sz(-1)
            rho(j-1,k-1,l)   = rho(j-1,k-1,l)     + sx1*sz(0)
            rho(j,k-1,l)     = rho(j,k-1,l)       + sx2*sz(0)
            rho(j+1,k-1,l)   = rho(j+1,k-1,l)     + sx3*sz(0)
            rho(j-1,k,l)     = rho(j-1,k,l)       + sx4*sz(0)
            rho(j,k,l)       = rho(j,k,l)         + sx5*sz(0)
            rho(j+1,k,l)     = rho(j+1,k,l)       + sx6*sz(0)
            rho(j-1,k+1,l)   = rho(j-1,k+1,l)     + sx7*sz(0)
            rho(j,k+1,l)     = rho(j,k+1,l)       + sx8*sz(0)
            rho(j+1,k+1,l)   = rho(j+1,k+1,l)     + sx9*sz(0)
            rho(j-1,k-1,l+1) = rho(j-1,k-1,l+1)   + sx1*sz(1)
            rho(j,k-1,l+1)   = rho(j,k-1,l+1)     + sx2*sz(1)
            rho(j+1,k-1,l+1) = rho(j+1,k-1,l+1)   + sx3*sz(1)
            rho(j-1,k,l+1)   = rho(j-1,k,l+1)     + sx4*sz(1)
            rho(j,k,l+1)     = rho(j,k,l+1)       + sx5*sz(1)
            rho(j+1,k,l+1)   = rho(j+1,k,l+1)     + sx6*sz(1)
            rho(j-1,k+1,l+1) = rho(j-1,k+1,l+1)   + sx7*sz(1)
            rho(j,k+1,l+1)   = rho(j,k+1,l+1)     + sx8*sz(1)
            rho(j+1,k+1,l+1) = rho(j+1,k+1,l+1)   + sx9*sz(1)
        END DO
        RETURN
    END SUBROUTINE depose_rho_scalar_2_2_2
    !!! --- Order 1 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized (HV-SCHEME)
    !!! --- This routine does vectorize on SIMD architecture with good performances
    !!! ---  Speedup>2 on AVX 256 bits
    SUBROUTINE depose_rho_vecHVv2_2_2_2(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num), DIMENSION(:,:), ALLOCATABLE:: rhocells
        INTEGER(idp), PARAMETER :: LVEC2=64
        INTEGER(idp), DIMENSION(LVEC2) :: ICELL, IG
        REAL(num) :: ww, wwx,wwy,wwz
        INTEGER(idp) :: NCELLS
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi
        REAL(num) :: xint,yint,zint,xintsq,yintsq,zintsq
        REAL(num) :: x,y,z,invvol, wq0, wq, szy, syy0,syy1,syy2,szz0,szz1,szz2
        REAL(num) :: sx0(LVEC2), sx1(LVEC2), sx2(LVEC2)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: ic,igrid,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
        INTEGER(idp) :: nnx, nnxy, off0, ind0
        INTEGER(idp) :: moff(1:8)
        REAL(num):: ww0(1:LVEC2,1:8),www(1:LVEC2,1:8)
        INTEGER(idp) :: orig, jorig, korig, lorig
        INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy

        ! Init parameters
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        wq0=q*invvol
        ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
        ncx=nx+3;ncy=ny+3;ncz=nz+3
        NCELLS=ncx*ncy*ncz
        ALLOCATE(rhocells(8,NCELLS))
        rhocells=0.0_num
        nnx = nx + 1 + 2*nxguard
        nnxy = nnx*(ny+1+2*nyguard)
        moff = (/-nnx-nnxy,-nnxy,nnx-nnxy,-nnx,nnx,-nnx+nnxy,nnxy,nnx+nnxy/)
	    ww0=0.0_num
        jorig=-1; korig=-1;lorig=-1
        orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
        ngx=(ngridx-ncx)
        ngxy=(ngridx*ngridy-ncx*ncy)
        ncxy=ncx*ncy
        ! FIRST LOOP: computes cell index of particle and their weight on vertices
        DO ip=1,np,LVEC2
#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
            !DIR$ ASSUME_ALIGNED w:64, sx0:64,sx1:64,sx2:64
            !DIR$ ASSUME_ALIGNED ICELL:64, IG:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,xp,yp,zp,w,sx0,sx1,sx2,ICELL,IG)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
            DO n=1,MIN(LVEC2,np-ip+1)
                nn=ip+n-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x= (xp(nn)-xmin)*dxi
                y = (yp(nn)-ymin)*dyi
                z = (zp(nn)-zmin)*dzi
                ! --- finds cell containing particles for current positions
                j=nint(x)
                k=nint(y)
                l=nint(z)
                ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
                IG(n)=ICELL(n)+(k-korig)*ngx+(l-lorig)*ngxy
                ! --- computes distance between particle and node for current positions
                xint = x-j
                yint = y-k
                zint = z-l
                xintsq=xint**2
                yintsq=yint**2
                zintsq=zint**2
                ! --- computes particles weights
                wq=w(nn)*wq0
                sx0(n)=0.5_num*(0.5_num-xint)**2
                sx1(n)=(0.75_num-xintsq)
                sx2(n)=0.5_num*(0.5_num+xint)**2
                syy0=0.5_num*(0.5_num-yint)**2
                syy1=(0.75_num-yintsq)
                syy2=0.5_num*(0.5_num+yint)**2
                szz0=0.5_num*(0.5_num-zint)**2*wq
                szz1=(0.75_num-zintsq)*wq
                szz2=0.5_num*(0.5_num+zint)**2*wq
                www(n,1) = syy0*szz0
                www(n,2) = syy1*szz0
                www(n,3) = syy2*szz0
                www(n,4) = syy0*szz1
                www(n,5) = syy2*szz1
                www(n,6) = syy0*szz2
                www(n,7) = syy1*szz2
                www(n,8) = syy2*szz2
                szy=syy1*szz1 ! central point
                ww0(n,1)=szy*sx0(n)
                ww0(n,2)=szy*sx1(n)
                ww0(n,3)=szy*sx2(n)
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            ! Current deposition on vertices
            DO n=1,MIN(LVEC2,np-ip+1)
                ! --- add charge density contributions to vertices of the current cell
                !DIR$ ASSUME_ALIGNED rhocells:64
#if defined __INTEL_COMPILER 
                !DIR$ ASSUME_ALIGNED rhocells:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,rhocells)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO nv=1,8 !!! - VECTOR
                    ww=www(n,nv)
                    ! Loop on (i=-1,j,k)
                    rhocells(nv,ICELL(n)-1)=rhocells(nv,ICELL(n)-1)+ww*sx0(n)
                    ! Loop on (i=0,j,k)
                    rhocells(nv,ICELL(n))=rhocells(nv,ICELL(n))+ww*sx1(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ICELL(n)+1)=rhocells(nv,ICELL(n)+1)+ww*sx2(n)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
                DO nv=1,4
                    rho(orig+IG(n)+nv-2)=rho(orig+IG(n)+nv-2)+ww0(n,nv)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            END DO
        END DO
        ! - reduction of rhocells in rho
        DO iz=1, ncz
            DO iy=1,ncy
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif
                DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                    ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                    igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                    rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1,ic)
                    rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2,ic)
                    rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3,ic)
                    rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4,ic)
                    rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5,ic)
                    rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6,ic)
                    rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7,ic)
                    rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8,ic)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            END DO
        END DO
        DEALLOCATE(rhocells)
        RETURN
    END SUBROUTINE depose_rho_vecHVv2_2_2_2


    !!! --- Order 3 3D scalar charge deposition routine
    !!! This version does not vectorize on SIMD architectures
    SUBROUTINE depose_rho_scalar_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), INTENT(IN OUT) :: rho
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol,sx1,sx2,sx3,sx4,sx5,sx6,sx7,sx8,sx9
        REAL(num) :: sx(-1:2), sy(-1:2), sz(-1:2)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        DO ip=1,np
            ! --- computes current position in grid units
            x = (xp(ip)-xmin)*dxi
            y = (yp(ip)-ymin)*dyi
            z = (zp(ip)-zmin)*dzi
            ! --- finds node of cell containing particles for current positions
            j=floor(x)
            k=floor(y)
            l=floor(z)
            ! --- computes distance between particle and node for current positions
            xint = x-j
            yint = y-k
            zint = z-l
            ! --- computes particles weights
            wq=q*w(ip)*invvol
            ! --- computes coefficients for node centered quantities
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx(-1) = onesixth*oxintsq*oxint
            sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx( 2) = onesixth*xintsq*xint
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy(-1) = onesixth*oyintsq*oyint
            sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy( 2) = onesixth*yintsq*yint
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz(-1) = onesixth*ozintsq*ozint
            sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz( 2) = onesixth*zintsq*zint
            ! --- add charge density contributions to the 64
            ! --- nearest vertices
            ! --- add charge density contributions
            rho(j-1,k-1,l-1)=rho(j-1,k-1,l-1)+sx(-1)*sy(-1)*sz(-1)*wq
            rho(j,k-1,l-1)=rho(j,k-1,l-1)+sx(0)*sy(-1)*sz(-1)*wq
            rho(j+1,k-1,l-1)=rho(j+1,k-1,l-1)+sx(1)*sy(-1)*sz(-1)*wq
            rho(j+2,k-1,l-1)=rho(j+2,k-1,l-1)+sx(2)*sy(-1)*sz(-1)*wq
            rho(j-1,k,l-1)=rho(j-1,k,l-1)+sx(-1)*sy(0)*sz(-1)*wq
            rho(j,k,l-1)=rho(j,k,l-1)+sx(0)*sy(0)*sz(-1)*wq
            rho(j+1,k,l-1)=rho(j+1,k,l-1)+sx(1)*sy(0)*sz(-1)*wq
            rho(j+2,k,l-1)=rho(j+2,k,l-1)+sx(2)*sy(0)*sz(-1)*wq
            rho(j-1,k+1,l-1)=rho(j-1,k+1,l-1)+sx(-1)*sy(1)*sz(-1)*wq
            rho(j,k+1,l-1)=rho(j,k+1,l-1)+sx(0)*sy(1)*sz(-1)*wq
            rho(j+1,k+1,l-1)=rho(j+1,k+1,l-1)+sx(1)*sy(1)*sz(-1)*wq
            rho(j+2,k+1,l-1)=rho(j+2,k+1,l-1)+sx(2)*sy(1)*sz(-1)*wq
            rho(j-1,k+2,l-1)=rho(j-1,k+2,l-1)+sx(-1)*sy(2)*sz(-1)*wq
            rho(j,k+2,l-1)=rho(j,k+2,l-1)+sx(0)*sy(2)*sz(-1)*wq
            rho(j+1,k+2,l-1)=rho(j+1,k+2,l-1)+sx(1)*sy(2)*sz(-1)*wq
            rho(j+2,k+2,l-1)=rho(j+2,k+2,l-1)+sx(2)*sy(2)*sz(-1)*wq
            rho(j-1,k-1,l)=rho(j-1,k-1,l)+sx(-1)*sy(-1)*sz(0)*wq
            rho(j,k-1,l)=rho(j,k-1,l)+sx(0)*sy(-1)*sz(0)*wq
            rho(j+1,k-1,l)=rho(j+1,k-1,l)+sx(1)*sy(-1)*sz(0)*wq
            rho(j+2,k-1,l)=rho(j+2,k-1,l)+sx(2)*sy(-1)*sz(0)*wq
            rho(j-1,k,l)=rho(j-1,k,l)+sx(-1)*sy(0)*sz(0)*wq
            rho(j,k,l)=rho(j,k,l)+sx(0)*sy(0)*sz(0)*wq
            rho(j+1,k,l)=rho(j+1,k,l)+sx(1)*sy(0)*sz(0)*wq
            rho(j+2,k,l)=rho(j+2,k,l)+sx(2)*sy(0)*sz(0)*wq
            rho(j-1,k+1,l)=rho(j-1,k+1,l)+sx(-1)*sy(1)*sz(0)*wq
            rho(j,k+1,l)=rho(j,k+1,l)+sx(0)*sy(1)*sz(0)*wq
            rho(j+1,k+1,l)=rho(j+1,k+1,l)+sx(1)*sy(1)*sz(0)*wq
            rho(j+2,k+1,l)=rho(j+2,k+1,l)+sx(2)*sy(1)*sz(0)*wq
            rho(j-1,k+2,l)=rho(j-1,k+2,l)+sx(-1)*sy(2)*sz(0)*wq
            rho(j,k+2,l)=rho(j,k+2,l)+sx(0)*sy(2)*sz(0)*wq
            rho(j+1,k+2,l)=rho(j+1,k+2,l)+sx(1)*sy(2)*sz(0)*wq
            rho(j+2,k+2,l)=rho(j+2,k+2,l)+sx(2)*sy(2)*sz(0)*wq
            rho(j-1,k-1,l+1)=rho(j-1,k-1,l+1)+sx(-1)*sy(-1)*sz(1)*wq
            rho(j,k-1,l+1)=rho(j,k-1,l+1)+sx(0)*sy(-1)*sz(1)*wq
            rho(j+1,k-1,l+1)=rho(j+1,k-1,l+1)+sx(1)*sy(-1)*sz(1)*wq
            rho(j+2,k-1,l+1)=rho(j+2,k-1,l+1)+sx(2)*sy(-1)*sz(1)*wq
            rho(j-1,k,l+1)=rho(j-1,k,l+1)+sx(-1)*sy(0)*sz(1)*wq
            rho(j,k,l+1)=rho(j,k,l+1)+sx(0)*sy(0)*sz(1)*wq
            rho(j+1,k,l+1)=rho(j+1,k,l+1)+sx(1)*sy(0)*sz(1)*wq
            rho(j+2,k,l+1)=rho(j+2,k,l+1)+sx(2)*sy(0)*sz(1)*wq
            rho(j-1,k+1,l+1)=rho(j-1,k+1,l+1)+sx(-1)*sy(1)*sz(1)*wq
            rho(j,k+1,l+1)=rho(j,k+1,l+1)+sx(0)*sy(1)*sz(1)*wq
            rho(j+1,k+1,l+1)=rho(j+1,k+1,l+1)+sx(1)*sy(1)*sz(1)*wq
            rho(j+2,k+1,l+1)=rho(j+2,k+1,l+1)+sx(2)*sy(1)*sz(1)*wq
            rho(j-1,k+2,l+1)=rho(j-1,k+2,l+1)+sx(-1)*sy(2)*sz(1)*wq
            rho(j,k+2,l+1)=rho(j,k+2,l+1)+sx(0)*sy(2)*sz(1)*wq
            rho(j+1,k+2,l+1)=rho(j+1,k+2,l+1)+sx(1)*sy(2)*sz(1)*wq
            rho(j+2,k+2,l+1)=rho(j+2,k+2,l+1)+sx(2)*sy(2)*sz(1)*wq
            rho(j-1,k-1,l+2)=rho(j-1,k-1,l+2)+sx(-1)*sy(-1)*sz(2)*wq
            rho(j,k-1,l+2)=rho(j,k-1,l+2)+sx(0)*sy(-1)*sz(2)*wq
            rho(j+1,k-1,l+2)=rho(j+1,k-1,l+2)+sx(1)*sy(-1)*sz(2)*wq
            rho(j+2,k-1,l+2)=rho(j+2,k-1,l+2)+sx(2)*sy(-1)*sz(2)*wq
            rho(j-1,k,l+2)=rho(j-1,k,l+2)+sx(-1)*sy(0)*sz(2)*wq
            rho(j,k,l+2)=rho(j,k,l+2)+sx(0)*sy(0)*sz(2)*wq
            rho(j+1,k,l+2)=rho(j+1,k,l+2)+sx(1)*sy(0)*sz(2)*wq
            rho(j+2,k,l+2)=rho(j+2,k,l+2)+sx(2)*sy(0)*sz(2)*wq
            rho(j-1,k+1,l+2)=rho(j-1,k+1,l+2)+sx(-1)*sy(1)*sz(2)*wq
            rho(j,k+1,l+2)=rho(j,k+1,l+2)+sx(0)*sy(1)*sz(2)*wq
            rho(j+1,k+1,l+2)=rho(j+1,k+1,l+2)+sx(1)*sy(1)*sz(2)*wq
            rho(j+2,k+1,l+2)=rho(j+2,k+1,l+2)+sx(2)*sy(1)*sz(2)*wq
            rho(j-1,k+2,l+2)=rho(j-1,k+2,l+2)+sx(-1)*sy(2)*sz(2)*wq
            rho(j,k+2,l+2)=rho(j,k+2,l+2)+sx(0)*sy(2)*sz(2)*wq
            rho(j+1,k+2,l+2)=rho(j+1,k+2,l+2)+sx(1)*sy(2)*sz(2)*wq
            rho(j+2,k+2,l+2)=rho(j+2,k+2,l+2)+sx(2)*sy(2)*sz(2)*wq
        END DO
        RETURN
    END SUBROUTINE depose_rho_scalar_3_3_3

    !!! --- Order 3 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized (HV-SCHEME)
    !!! --- This routine does vectorize on SIMD architecture with good performances
    !!! ---  Speedup>2 on AVX 256 bits
    SUBROUTINE depose_rho_vecHVv2_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num), DIMENSION(:,:), ALLOCATABLE:: rhocells
        INTEGER(idp), PARAMETER :: LVEC2=8
        INTEGER(idp), DIMENSION(LVEC2) :: ICELL
        REAL(num) :: ww, wwx,wwy,wwz
        INTEGER(idp) :: NCELLS
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,invvol, wq0, wq
        REAL(num) :: sx1(LVEC2), sx2(LVEC2), sx3(LVEC2),sx4(LVEC2)
        REAL(num) :: sy(-1:2), sz(-1:2)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: ic,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
        INTEGER(idp) :: nnx, nnxy, off0, ind0
        INTEGER(idp) :: moff(1:16)
        REAL(num):: www(1:LVEC2,1:16)

        ! Init parameters
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        wq0=q*invvol
        NCELLS=(2*nxguard+nx)*(2*nyguard+ny)*(2*nzguard+nz)
        ALLOCATE(rhocells(16,NCELLS))
        rhocells=0.0_num
        nnx = nx + 1 + 2*nxguard
        nnxy = (nx+1+2*nxguard)*(ny+1+2*nyguard)
        moff = (/-1_idp-nnx-nnxy,-nnx-nnxy,1_idp-nnx-nnxy,-1_idp-nnxy,-nnxy,1_idp-nnxy,-1_idp+nnx-nnxy,nnx-nnxy, &
                -1_idp-nnx-nnxy,-nnx-nnxy,1_idp-nnx-nnxy,-1_idp-nnxy,-nnxy,1_idp-nnxy,-1_idp+nnx-nnxy,nnx-nnxy/)
        off0=1+nnx+nnxy

        ! FIRST LOOP: computes cell index of particle and their weight on vertices
        DO ip=1,np,LVEC2

#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
            !DIR$ ASSUME_ALIGNED w:64
            !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,xp,yp,zp,w,ICELL)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
            DO n=1,MIN(LVEC2,np-ip+1)
                nn=ip+n-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x= (xp(nn)-xmin)*dxi
                y = (yp(nn)-ymin)*dyi
                z = (zp(nn)-zmin)*dzi
                ! --- finds cell containing particles for current positions
                j=floor(x)
                k=floor(y)
                l=floor(z)
                ICELL(n)=1+j+nxguard+(k+nyguard+1)*(nx+2*nxguard)+(l+nzguard+1)*(ny+2*nyguard)
                wq=w(nn)*wq0
                ! --- computes distance between particle and node for current positions
                xint = x-j
                yint = y-k
                zint = z-l
                ! --- computes coefficients for node centered quantities
                oxint = 1.0_num-xint
                xintsq = xint*xint
                oxintsq = oxint*oxint
                sx1(n) = onesixth*oxintsq*oxint
                sx2(n) = twothird-xintsq*(1.0_num-xint*0.5_num)
                sx3(n) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
                sx4(n) = onesixth*xintsq*xint
                oyint = 1.0_num-yint
                yintsq = yint*yint
                oyintsq = oyint*oyint
                sy(-1) = onesixth*oyintsq*oyint
                sy( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)
                sy( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)
                sy( 2) = onesixth*yintsq*yint
                ozint = 1.0_num-zint
                zintsq = zint*zint
                ozintsq = ozint*ozint
                sz(-1) = onesixth*ozintsq*ozint*wq
                sz( 0) = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
                sz( 1) = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
                sz( 2) = onesixth*zintsq*zint*wq
                www(n,1)=sy(-1)*sz(-1)
                www(n,2)=sy(0)*sz(-1)
                www(n,3)=sy(1)*sz(-1)
                www(n,4)=sy(2)*sz(-1)
                www(n,5)=sy(-1)*sz(0)
                www(n,6)=sy(0)*sz(0)
                www(n,7)=sy(1)*sz(0)
                www(n,8)=sy(2)*sz(0)
                www(n,9)=sy(-1)*sz(1)
                www(n,10)=sy(0)*sz(1)
                www(n,11)=sy(1)*sz(1)
                www(n,12)=sy(2)*sz(1)
                www(n,13)=sy(-1)*sz(2)
                www(n,14)=sy(0)*sz(2)
                www(n,15)=sy(1)*sz(2)
                www(n,16)=sy(2)*sz(2)
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            ! Current deposition on vertices
            DO n=1,MIN(LVEC2,np-ip+1)
                ! --- add charge density contributions to vertices of the current cell
                ic=ICELL(n)
#if defined __INTEL_COMPILER 
                !DIR$ ASSUME_ALIGNED rhocells:64
#elif defined __IBMBGQ__
                !IBM* ALIGN(64,rhocells)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
                DO nv=1,16 !!! - VECTOR
                    ww=www(n,nv)
                    !ww=0.1_num
                    ! Loop on (i=-1,j,k)
                    rhocells(nv,ic-1)=rhocells(nv,ic-1)+ww*sx1(n)
                    ! Loop on (i=0,j,k)
                    rhocells(nv,ic)=rhocells(nv,ic)+ww*sx2(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic+1)=rhocells(nv,ic+1)+ww*sx3(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic+2)=rhocells(nv,ic+2)+ww*sx4(n)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            END DO
        END DO
        DO nv=1,16
            ind0=off0+moff(nv)
#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED rhocells:64
            !DIR$ ASSUME_ALIGNED rho:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,rhocells,rho)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
            DO ic=1,NCELLS  !!! VECTOR
                rho(ic+ind0)=rho(ic+ind0)+rhocells(nv,ic)
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
        END DO
        DEALLOCATE(rhocells)
        RETURN
    END SUBROUTINE depose_rho_vecHVv2_3_3_3


    !!! --- Order 3 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized (HV-SCHEME)
    !!! --- This routine does vectorize on SIMD architecture with good performances
    !!! ---  Speedup>2 on AVX 256 bits
    SUBROUTINE depose_rho_vecHVv3_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num), DIMENSION(:,:), ALLOCATABLE:: rhocells
        INTEGER(idp), PARAMETER :: LVEC2=16
        INTEGER(idp), DIMENSION(LVEC2) :: ICELL
        REAL(num) :: ww, wwx,wwy,wwz
        INTEGER(idp) :: NCELLS
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint(1:LVEC2), &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,invvol, wq0, wq
        REAL(num) :: sx1(LVEC2), sx2(LVEC2), sx3(LVEC2),sx4(LVEC2), sy1(LVEC2), sy2(LVEC2), sy3(LVEC2),sy4(LVEC2), &
                     sz1, sz2, sz3,sz4
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: ic, igrid, ic0,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
        INTEGER(idp) :: nnx, nnxy, off0, ind0
        INTEGER(idp) :: moff(1:8)
        REAL(num):: www(1:16,1:LVEC2), zdec(1:8), h1(1:8), h11(1:8), h12(1:8), sgn(1:8), szz(1:8)
        INTEGER(idp) :: orig, jorig, korig, lorig
        INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy

        ! Init parameters
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        wq0=q*invvol
        ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
        ncx=nx+5; ncy=ny+4; ncz=nz+2
        NCELLS=ncx*ncy*ncz
        ALLOCATE(rhocells(8,NCELLS))
        rhocells=0_num; www=0.0_num
        nnx = ngridx
        nnxy = ngridx*ngridy
        moff = (/-nnxy,0_idp,nnxy,2_idp*nnxy,nnx-nnxy,nnx,nnx+nnxy,nnx+2_idp*nnxy/)
        jorig=-2_idp; korig=-2_idp;lorig=-1_idp
        orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
        ngx=(ngridx-ncx)
        ngxy=(ngridx*ngridy-ncx*ncy)
        ncxy=ncx*ncy

        h1(1:4)=(/1_num,0_num,1_num,0_num/); sgn(1:4)=(/1_num,-1_num,1_num,-1_num/)
        h11(1:4)=(/0_num,1_num,1_num,0_num/); h12(1:4)=(/1_num,0_num,0_num,1_num/)

        ! FIRST LOOP: computes cell index of particle and their weight on vertices
        DO ip=1,np,LVEC2
#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
            !DIR$ ASSUME_ALIGNED w:64,sx1:64,sx2:64,sx3:64,sx4:64
            !DIR$ ASSUME_ALIGNED w:64,sy1:64,sy2:64,sy3:64,sy4:64
            !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,xp,yp,zp,w,sx1,sx2,sx3,sx4,sy1,sy2,sy3,sy4,ICELL)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
            DO n=1,MIN(LVEC2,np-ip+1)
                nn=ip+n-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x= (xp(nn)-xmin)*dxi
                y = (yp(nn)-ymin)*dyi
                z = (zp(nn)-zmin)*dzi
                ! --- finds cell containing particles for current positions
                j=floor(x)
                k=floor(y)
                l=floor(z)
                ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
                wq=w(nn)*wq0
                ! --- computes distance between particle and node for current positions
                xint = x-j
                yint= y-k
                zint(n) = z-l
                ! --- computes coefficients for node centered quantities
                oxint = 1.0_num-xint
                xintsq = xint*xint
                oxintsq = oxint*oxint
                sx1(n) = onesixth*oxintsq*oxint
                sx2(n) = twothird-xintsq*(1.0_num-xint*0.5_num)
                sx3(n) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
                sx4(n) = onesixth*xintsq*xint
                oyint = 1.0_num-yint
                yintsq = yint*yint
                oyintsq = oyint*oyint
                sy1(n) = onesixth*oyintsq*oyint*wq
                sy2(n) = (twothird-yintsq*(1.0_num-yint*0.5_num))*wq
                sy3(n) = (twothird-oyintsq*(1.0_num-oyint*0.5_num))*wq
                sy4(n) = onesixth*yintsq*yint*wq
!                ozint = 1.0_num-zint(n)
!                zintsq = zint(n)*zint(n)
!                ozintsq = ozint*ozint
!                sz1 = onesixth*ozintsq*ozint*wq
!                sz2 = (twothird-zintsq*(1.0_num-zint(n)*0.5_num))*wq
!                sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
!                sz4 = onesixth*zintsq*zint(n)*wq
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            DO n=1,MIN(LVEC2,np-ip+1)

#if defined __INTEL_COMPILER 
                !$DIR ASSUME_ALIGNED www:64, h1:64, h11:64, h12:64, zdec:64, sgn:64, szz:64
#elif defined __IBMBGQ__
                !IBM* ALIGN(64,www,h1,h11,h12,zdec,sgn,szz)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO nv=1,4 !!! Vector
                    zdec(nv)     = (h1(nv)-zint(n))*sgn(nv)
                    szz(nv)      = (twothird-zdec(nv)**2*(1.0_num-zdec(nv)*0.5_num))*h11(nv)+ &
                                    onesixth*zdec(nv)**3*h12(nv)
                    www(nv,n)    = szz(nv)*sy1(n)
                    www(nv+4,n)  = szz(nv)*sy2(n)
                    www(nv+8,n)  = szz(nv)*sy3(n)
                    www(nv+12,n) = szz(nv)*sy4(n)
                ENDDO
            END DO
            ! Current deposition on vertices
            DO n=1,MIN(LVEC2,np-ip+1)
                ! --- add charge density contributions to vertices of the current cell
                ic=ICELL(n)

#if defined __INTEL_COMPILER 
                !DIR$ ASSUME_ALIGNED rhocells:64, www:64
#elif defined __IBMBGQ__
                !IBM* ALIGN(64,rhocells,www)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO nv=1,8 !!! - VECTOR
                    ! Loop on (i=-1,j,k)
                    rhocells(nv,ic-ncx-1) = rhocells(nv,ic-ncx-1) + www(nv,n)*sx1(n)
                    ! Loop on (i=0,j,k)
                    rhocells(nv,ic-ncx)   = rhocells(nv,ic-ncx)   + www(nv,n)*sx2(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic-ncx+1) = rhocells(nv,ic-ncx+1) + www(nv,n)*sx3(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic-ncx+2) = rhocells(nv,ic-ncx+2) + www(nv,n)*sx4(n)
                    ! Loop on (i=-1,j,k)
                    rhocells(nv,ic+ncx-1) = rhocells(nv,ic+ncx-1) + www(nv+8,n)*sx1(n)
                    ! Loop on (i=0,j,k)
                    rhocells(nv,ic+ncx)   = rhocells(nv,ic+ncx)   + www(nv+8,n)*sx2(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic+ncx+1) = rhocells(nv,ic+ncx+1) + www(nv+8,n)*sx3(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic+ncx+2) = rhocells(nv,ic+ncx+2) + www(nv+8,n)*sx4(n)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	    !$OMP END SIMD
#endif
            END DO
        END DO
        ! - reduction of rhocells in rho
        DO iz=1, ncz
            DO iy=1,ncy
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                    ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                    igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                    rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1,ic)
                    rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2,ic)
                    rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3,ic)
                    rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4,ic)
                    rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5,ic)
                    rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6,ic)
                    rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7,ic)
                    rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8,ic)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	     !$OMP END SIMD
#endif
            END DO
        END DO
        DEALLOCATE(rhocells)
        RETURN
    END SUBROUTINE depose_rho_vecHVv3_3_3_3

    !!! --- Order 3 3D vector charge deposition routine
    !!! --- Computes charge density on grid vectorized (HV-SCHEME)
    !!! --- This routine does vectorize on SIMD architecture with good performances
    !!! ---  Speedup>2 on AVX 256 bits
    SUBROUTINE depose_rho_vecHVv4_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: rho(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num), DIMENSION(:,:), ALLOCATABLE:: rhocells
        INTEGER(idp), PARAMETER :: LVEC2=16
        INTEGER(idp), DIMENSION(LVEC2) :: ICELL
        REAL(num) :: ww, wwx,wwy,wwz
        INTEGER(idp) :: NCELLS
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint(1:LVEC2), &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,invvol, wq0, wq
        REAL(num) :: sx1(LVEC2), sx2(LVEC2), sx3(LVEC2),sx4(LVEC2), sy1, sy2, sy3,sy4, &
                     sz1, sz2, sz3,sz4, w1,w2
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp):: ic, igrid, ic0,j,k,l,vv,n,ip,jj,kk,ll,nv,nn
        INTEGER(idp) :: nnx, nnxy, off0, ind0
        INTEGER(idp) :: moff(1:8)
        REAL(num):: www1(LVEC2,8),www2(LVEC2,8), zdec(1:8), h1(1:8), h11(1:8), h12(1:8), sgn(1:8), szz(1:8)
        INTEGER(idp) :: orig, jorig, korig, lorig
        INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy

        ! Init parameters
        dxi = 1.0_num/dx
        dyi = 1.0_num/dy
        dzi = 1.0_num/dz
        invvol = dxi*dyi*dzi
        wq0=q*invvol
        ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
        ncx=nx+5; ncy=ny+4; ncz=nz+2
        NCELLS=ncx*ncy*ncz
        ALLOCATE(rhocells(8,NCELLS))
        rhocells=0_num
        nnx = ngridx
        nnxy = ngridx*ngridy
        moff = (/-nnxy,0_idp,nnxy,2_idp*nnxy,nnx-nnxy,nnx,nnx+nnxy,nnx+2_idp*nnxy/)
        jorig=-2_idp; korig=-2_idp;lorig=-1_idp
        orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
        ngx=(ngridx-ncx)
        ngxy=(ngridx*ngridy-ncx*ncy)
        ncxy=ncx*ncy

        ! FIRST LOOP: computes cell index of particle and their weight on vertices
        DO ip=1,np,LVEC2

#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
            !DIR$ ASSUME_ALIGNED w:64,sx1:64,sx2:64,sx3:64,sx4:64
            !DIR$ ASSUME_ALIGNED ICELL:64, www1:64, www2:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,xp,yp,zp,w,sx1,sx2,sx3,sx4,ICELL,www1,www2)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!$DIR SIMD 
#endif 
            DO n=1,MIN(LVEC2,np-ip+1)
                nn=ip+n-1
                ! Calculation relative to particle n
                ! --- computes current position in grid units
                x= (xp(nn)-xmin)*dxi
                y = (yp(nn)-ymin)*dyi
                z = (zp(nn)-zmin)*dzi
                ! --- finds cell containing particles for current positions
                j=floor(x)
                k=floor(y)
                l=floor(z)
                ICELL(n)=1+(j-jorig)+(k-korig)*(ncx)+(l-lorig)*ncxy
                wq=w(nn)*wq0
                ! --- computes distance between particle and node for current positions
                xint = x-j
                yint= y-k
                zint(n) = z-l
                ! --- computes coefficients for node centered quantities
                oxint = 1.0_num-xint
                xintsq = xint*xint
                oxintsq = oxint*oxint
                sx1(n) = onesixth*oxintsq*oxint
                sx2(n) = twothird-xintsq*(1.0_num-xint*0.5_num)
                sx3(n) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
                sx4(n) = onesixth*xintsq*xint
                oyint = 1.0_num-yint
                yintsq = yint*yint
                oyintsq = oyint*oyint
                sy1 = onesixth*oyintsq*oyint
                sy2 = (twothird-yintsq*(1.0_num-yint*0.5_num))
                sy3 = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
                sy4 = onesixth*yintsq*yint
                ozint = 1.0_num-zint(n)
                zintsq = zint(n)*zint(n)
                ozintsq = ozint*ozint
                sz1 = onesixth*ozintsq*ozint*wq
                sz2 = (twothird-zintsq*(1.0_num-zint(n)*0.5_num))*wq
                sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
                sz4 = onesixth*zintsq*zint(n)*wq
                www1(n,1)=sz1*sy1
                www1(n,2)=sz2*sy1
                www1(n,3)=sz3*sy1
                www1(n,4)=sz4*sy1
                www1(n,5)=sz1*sy2
                www1(n,6)=sz2*sy2
                www1(n,7)=sz3*sy2
                www1(n,8)=sz4*sy2
                www2(n,1)=sz1*sy3
                www2(n,2)=sz2*sy3
                www2(n,3)=sz3*sy3
                www2(n,4)=sz4*sy3
                www2(n,5)=sz1*sy4
                www2(n,6)=sz2*sy4
                www2(n,7)=sz3*sy4
                www2(n,8)=sz4*sy4
            END DO
#if defined _OPENMP && _OPENMP>=201307
       	   !$OMP END SIMD
#endif
            ! Current deposition on vertices
            DO n=1,MIN(LVEC2,np-ip+1)
                ! --- add charge density contributions to vertices of the current cell
                ic=ICELL(n)
#if defined __INTEL_COMPILER 
                !DIR$ ASSUME_ALIGNED rhocells:64, www1:64, www2:64
#elif defined __IBMBGQ__
                !IBM* ALIGN(64,rhocells,www1,www2)
#endif 
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO nv=1,8 !!! - VECTOR
                    w1=www1(n,nv)
                    ! Loop on (i=-1,j,k)
                    rhocells(nv,ic-ncx-1) = rhocells(nv,ic-ncx-1) + w1*sx1(n)
                    ! Loop on (i=0,j,k)
                    rhocells(nv,ic-ncx)   = rhocells(nv,ic-ncx)   + w1*sx2(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic-ncx+1) = rhocells(nv,ic-ncx+1) + w1*sx3(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic-ncx+2) = rhocells(nv,ic-ncx+2) + w1*sx4(n)

                    w2=www2(n,nv)
                    ! Loop on (i=-1,j,k)
                    rhocells(nv,ic+ncx-1) = rhocells(nv,ic+ncx-1) + w2*sx1(n)
                    ! Loop on (i=0,j,k)
                    rhocells(nv,ic+ncx)   = rhocells(nv,ic+ncx)   + w2*sx2(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic+ncx+1) = rhocells(nv,ic+ncx+1) + w2*sx3(n)
                    !Loop on (i=1,j,k)
                    rhocells(nv,ic+ncx+2) = rhocells(nv,ic+ncx+2) + w2*sx4(n)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	        !$OMP END SIMD
#endif
            END DO
        END DO
        ! - reduction of rhocells in rho
        DO iz=1, ncz
            DO iy=1,ncy
#if defined _OPENMP && _OPENMP>=201307
			    !$OMP SIMD 
#elif defined __IBMBGQ__
			    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			    !$DIR SIMD 
#endif 
                DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                    ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                    igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                    rho(orig+igrid+moff(1))=rho(orig+igrid+moff(1))+rhocells(1,ic)
                    rho(orig+igrid+moff(2))=rho(orig+igrid+moff(2))+rhocells(2,ic)
                    rho(orig+igrid+moff(3))=rho(orig+igrid+moff(3))+rhocells(3,ic)
                    rho(orig+igrid+moff(4))=rho(orig+igrid+moff(4))+rhocells(4,ic)
                    rho(orig+igrid+moff(5))=rho(orig+igrid+moff(5))+rhocells(5,ic)
                    rho(orig+igrid+moff(6))=rho(orig+igrid+moff(6))+rhocells(6,ic)
                    rho(orig+igrid+moff(7))=rho(orig+igrid+moff(7))+rhocells(7,ic)
                    rho(orig+igrid+moff(8))=rho(orig+igrid+moff(8))+rhocells(8,ic)
                END DO
#if defined _OPENMP && _OPENMP>=201307
       	        !$OMP END SIMD
#endif
            END DO
        END DO
        DEALLOCATE(rhocells)
        RETURN
    END SUBROUTINE depose_rho_vecHVv4_3_3_3

    !!! --- General charge deposition routine (Warning: Highly unoptimized routine)
    !!! Computes charge density on grid at arbitrary orders nox, noy and noz
    SUBROUTINE pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz, &
                        l_particles_weight, l4symtry)
        IMPLICIT NONE
        INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: rho
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        LOGICAL :: l_particles_weight, l4symtry

        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol
        REAL(num) :: sx(-int(nox/2):int((nox+1)/2)), &
                     sy(-int(noy/2):int((noy+1)/2)), &
                     sz(-int(noz/2):int((noz+1)/2))
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER(idp) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
   
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

        DO ip=1,np
        
            ! --- computes current position in grid units
            x = (xp(ip)-xmin)*dxi
            y = (yp(ip)-ymin)*dyi
            z = (zp(ip)-zmin)*dzi
        
            ! --- applies 4-fold symmetry
            IF (l4symtry) THEN
                x=abs(x)
                y=abs(y)
            END IF
      
            ! --- finds node of cell containing particles for current positions
            ! --- (different for odd/even spline orders)
            IF (nox==2*(nox/2)) THEN
                j=nint(x)
            ELSE
                j=floor(x)
            END IF
            IF (noy==2*(noy/2)) THEN
                k=nint(y)
            ELSE
                k=floor(y)
            END IF
            IF(noz==2*(noz/2)) THEN
                l=nint(z)
            ELSE
                l=floor(z)
            END IF

            ! --- computes distance between particle and node for current positions
            xint = x-j
            yint = y-k
            zint = z-l

            ! --- computes particles "weights"
            IF (l_particles_weight) THEN
                wq=q*w(ip)*invvol
            ELSE
                wq=q*invvol*w(1)
            ENDIF
      
            ! --- computes coefficients for node centered quantities
            SELECT CASE(nox)
            CASE(0)
                sx( 0) = 1.0_num
            CASE(1)
                sx( 0) = 1.0_num-xint
                sx( 1) = xint
            CASE(2)
                xintsq = xint*xint
                sx(-1) = 0.5_num*(0.5_num-xint)**2
                sx( 0) = 0.75_num-xintsq
                sx( 1) = 0.5_num*(0.5_num+xint)**2
            CASE(3)
                oxint = 1.0_num-xint
                xintsq = xint*xint
                oxintsq = oxint*oxint
                sx(-1) = onesixth*oxintsq*oxint
                sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
                sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
                sx( 2) = onesixth*xintsq*xint
            END SELECT

            SELECT CASE(noy)
            CASE(0)
                sy( 0) = 1.0_num
            CASE(1)
                sy( 0) = 1.0_num-yint
                sy( 1) = yint
            CASE(2)
                yintsq = yint*yint
                sy(-1) = 0.5_num*(0.5_num-yint)**2
                sy( 0) = 0.75_num-yintsq
                sy( 1) = 0.5_num*(0.5_num+yint)**2
            CASE(3)
                oyint = 1.0_num-yint
                yintsq = yint*yint
                oyintsq = oyint*oyint
                sy(-1) = onesixth*oyintsq*oyint
                sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
                sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
                sy( 2) = onesixth*yintsq*yint
            END SELECT

            SELECT CASE(noz)
            CASE(0)
                sz( 0) = 1.0_num
            CASE(1)
                sz( 0) = 1.0_num-zint
                sz( 1) = zint
            CASE(2)
                zintsq = zint*zint
                sz(-1) = 0.5_num*(0.5_num-zint)**2
                sz( 0) = 0.75_num-zintsq
                sz( 1) = 0.5_num*(0.5_num+zint)**2
            CASE(3)
                ozint = 1.0_num-zint
                zintsq = zint*zint
                ozintsq = ozint*ozint
                sz(-1) = onesixth*ozintsq*ozint
                sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
                sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
                sz( 2) = onesixth*zintsq*zint
            END SELECT

            ! --- add charge density contributions
            DO ll = izmin, izmax
                DO kk = iymin, iymax
                        DO jj = ixmin, ixmax
                            rho(j+jj,k+kk,l+ll)=rho(j+jj,k+kk,l+ll)+sx(jj)*sy(kk)*sz(ll)*wq
                        END DO
                END DO
            END DO
        END DO
        RETURN
    END SUBROUTINE pxr_depose_rho_n
    
    ! 2D Charge deposition at nth order 
	subroutine pxr_depose_rho_n_2dxz(rho,np,xp,yp,zp,w,q,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard,nox,noz, &
							l_particles_weight,l4symtry,l_2drz, type_rz_depose)
	   use constants 
	   implicit none
	   integer(idp) :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
	   real(num), dimension(-nxguard:nx+nxguard,0:0,-nzguard:nz+nzguard), intent(in out) :: rho
	   real(num), dimension(np) :: xp,yp,zp,w
	   real(num) :: q,dt,dx,dz,xmin,zmin
	   logical(idp) :: l_particles_weight,l4symtry,l_2drz

	   real(num) :: dxi,dzi,xint,zint, &
					   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
	   real(num) :: x,z,r,wq,invvol
	   real(num) :: sx(-int(nox/2):int((nox+1)/2)), &
					   sz(-int(noz/2):int((noz+1)/2))
	   real(num), parameter :: onesixth=1./6.,twothird=2./3.
	   integer(idp) :: j,l,ip,jj,ll,ixmin, ixmax, izmin, izmax
   
		  dxi = 1./dx
		  dzi = 1./dz
		  invvol = dxi*dzi

		  ! Davoine method : limited to order 1 in r
		  if (type_rz_depose==2) then
			 nox = 1
		  endif

		  ixmin = -int(nox/2)
		  ixmax = int((nox+1)/2)
		  izmin = -int(noz/2)
		  izmax = int((noz+1)/2)

		  do ip=1,np
		
			! --- computes current position in grid units
			if (l_2drz) then
			  r = sqrt(xp(ip)*xp(ip)+yp(ip)*yp(ip))
			  x = (r-xmin)*dxi
			  z = (zp(ip)-zmin)*dzi
			else
			  x = (xp(ip)-xmin)*dxi
			  z = (zp(ip)-zmin)*dzi
			end if
		
			! --- applies 4-fold symmetry
			if (l4symtry) then
			  x=abs(x)
			end if
		
			! --- finds node of cell containing particles for current positions 
			! --- (different for odd/even spline orders)
			if (nox==2*(nox/2)) then
			  j=nint(x)
			else
			  j=floor(x)
			end if
			if (noz==2*(noz/2)) then
			  l=nint(z)
			else
			  l=floor(z)
			end if

			! --- computes distance between particle and node for current positions
			xint = x-j
			zint = z-l

			! --- computes particles "weights"
			if (l_particles_weight) then
			  wq=q*w(ip)*invvol
			else
			  wq=q*invvol
			end if
	  
			! --- computes coefficients for node centered quantities
			if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
			   sx(0) = 1. - xint  + 1./(4*j+2)*( -xint + xint**2 )
			   sx(1) = 1. - sx(0)
			else                          ! Standard method, canonical shapes in r
			   select case(nox)
			   case(0)
				  sx( 0) = 1.
			   case(1)
				  sx( 0) = 1.-xint
				  sx( 1) = xint
			   case(2)
				  xintsq = xint*xint
				  sx(-1) = 0.5*(0.5-xint)**2
				  sx( 0) = 0.75-xintsq
				  sx( 1) = 0.5*(0.5+xint)**2
			   case(3)
				  oxint = 1.-xint
				  xintsq = xint*xint
				  oxintsq = oxint*oxint
				  sx(-1) = onesixth*oxintsq*oxint
				  sx( 0) = twothird-xintsq*(1.-xint/2)
				  sx( 1) = twothird-oxintsq*(1.-oxint/2)
				  sx( 2) = onesixth*xintsq*xint
			   end select
			endif

			select case(noz)
			 case(0)
			  sz( 0) = 1.
			 case(1)
			  sz( 0) = 1.-zint
			  sz( 1) = zint
			 case(2)
			  zintsq = zint*zint
			  sz(-1) = 0.5*(0.5-zint)**2
			  sz( 0) = 0.75-zintsq
			  sz( 1) = 0.5*(0.5+zint)**2
			 case(3)
			  ozint = 1.-zint
			  zintsq = zint*zint
			  ozintsq = ozint*ozint
			  sz(-1) = onesixth*ozintsq*ozint
			  sz( 0) = twothird-zintsq*(1.-zint/2)
			  sz( 1) = twothird-ozintsq*(1.-ozint/2)
			  sz( 2) = onesixth*zintsq*zint
			end select        

			! --- add charge density contributions
			 do ll = izmin, izmax
				do jj = ixmin, ixmax
				  rho(j+jj,0,l+ll)=rho(j+jj,0,l+ll)+sx(jj)*sz(ll)*wq
				end do
			end do

		end do

	  return
	end subroutine pxr_depose_rho_n_2dxz

	subroutine pxr_depose_rhoold_n_2dxz(rhoold,np,xp,zp,ux,uy,uz,gaminv,w,q,xmin,zmin,dt,dx,dz,nx,nz,nxguard,nzguard,nox,noz, &
							l_particles_weight,l4symtry)
	   use constants 
	   implicit none
	   integer(idp) :: np,nx,nz,nox,noz,nxguard,nzguard
	   real(num), dimension(-nxguard:nx+nxguard,0:0,-nzguard:nz+nzguard), intent(in out) :: rhoold
	   real(num), dimension(np) :: xp,zp,w,ux,uy,uz,gaminv
	   real(num) :: q,dt,dx,dz,xmin,zmin
	   logical(idp) :: l_particles_weight,l4symtry

	   real(num) :: dxi,dzi,xint,zint, &
					   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
	   real(num) :: xintold,zintold, &
					   oxintold,ozintold
	   real(num) :: x,z,xold,zold,wq,invvol,vx,vy,vz
	   real(num) :: sx(-int(nox/2):int((nox+1)/2)), &
					   sz(-int(noz/2):int((noz+1)/2))
	   real(num) :: sxold(-int(nox/2):int((nox+1)/2)), &
					   szold(-int(noz/2):int((noz+1)/2))
	   real(num), parameter :: onesixth=1./6.,twothird=2./3.
	   integer(idp) :: j,l,ip,jj,ll,jold,lold,ixmin, ixmax, izmin, izmax, istep, ndt,idt
	   real(num) :: dxp,dzp,x0,z0,x1,z1
	  
		  dxi = 1./dx
		  dzi = 1./dz
		  invvol = dxi*dzi

		  ixmin = -int(nox/2)
		  ixmax = int((nox+1)/2)
		  izmin = -int(noz/2)
		  izmax = int((noz+1)/2)
		  ndt = 1

		  do ip=1,np
	  
		   vx = ux(ip)*gaminv(ip)
		   vy = uy(ip)*gaminv(ip)
		   vz = uz(ip)*gaminv(ip)
				
		   x1 = (xp(ip)-xmin)*dxi
		   z1 = (zp(ip)-zmin)*dzi
		   x0 = x1 - vx*dt*dxi
		   z0 = z1 - vz*dt*dzi

		   dxp=(x1-x0)/ndt
		   dzp=(z1-z0)/ndt
	   
		   xold=x0
		   zold=z0

		   do idt=1,ndt
	   
			if (idt>1) then
			  xold=x
			  zold=z
			end if
			x=xold+dxp
			z=zold+dzp
		
			if (l4symtry) then
			  x=abs(x)
			  xold=abs(xold)
			end if
		
			! --- finds node of cell containing particles for current positions 
			! --- (different for odd/even spline orders)
			if (nox==2*(nox/2)) then
			  j=nint(x)
			else
			  j=floor(x)
			end if
			if (noz==2*(noz/2)) then
			  l=nint(z)
			else
			  l=floor(z)
			end if

			if (nox==2*(nox/2)) then
			  jold=nint(xold)
			else
			  jold=floor(xold)
			end if
			if (noz==2*(noz/2)) then
			  lold=nint(zold)
			else
			  lold=floor(zold)
			end if

			xint = x-j
			zint = z-l
			xintold = xold-jold
			zintold = zold-lold

			if (l_particles_weight) then
			  wq=q*w(ip)*invvol
			else
			  wq=q*w(1)*invvol
			end if
	  
			select case(nox)
			 case(0)
			  sxold( 0) = 1.
			 case(1)
			  sxold( 0) = 1.-xintold
			  sxold( 1) = xintold
			 case(2)
			  xintsq = xintold*xintold
			  sxold(-1) = 0.5*(0.5-xintold)**2
			  sxold( 0) = 0.75-xintsq
			  sxold( 1) = 0.5*(0.5+xintold)**2
			 case(3)
			  oxintold = 1.-xintold
			  xintsq = xintold*xintold
			  oxintsq = oxintold*oxintold
			  sxold(-1) = onesixth*oxintsq*oxintold
			  sxold( 0) = twothird-xintsq*(1.-xintold/2)
			  sxold( 1) = twothird-oxintsq*(1.-oxintold/2)
			  sxold( 2) = onesixth*xintsq*xintold
			end select        

			select case(noz)
			 case(0)
			  szold( 0) = 1.
			 case(1)
			  szold( 0) = 1.-zintold
			  szold( 1) = zintold
			 case(2)
			  zintsq = zintold*zintold
			  szold(-1) = 0.5*(0.5-zintold)**2
			  szold( 0) = 0.75-zintsq
			  szold( 1) = 0.5*(0.5+zintold)**2
			 case(3)
			  ozintold = 1.-zintold
			  zintsq = zintold*zintold
			  ozintsq = ozintold*ozintold
			  szold(-1) = onesixth*ozintsq*ozintold
			  szold( 0) = twothird-zintsq*(1.-zintold/2)
			  szold( 1) = twothird-ozintsq*(1.-ozintold/2)
			  szold( 2) = onesixth*zintsq*zintold
			end select 
		
			select case(nox)
			 case(0)
			  sx( 0) = 1.
			 case(1)
			  sx( 0) = 1.-xint
			  sx( 1) = xint
			 case(2)
			  xintsq = xint*xint
			  sx(-1) = 0.5*(0.5-xint)**2
			  sx( 0) = 0.75-xintsq
			  sx( 1) = 0.5*(0.5+xint)**2
			 case(3)
			  oxint = 1.-xint
			  xintsq = xint*xint
			  oxintsq = oxint*oxint
			  sx(-1) = onesixth*oxintsq*oxint
			  sx( 0) = twothird-xintsq*(1.-xint/2)
			  sx( 1) = twothird-oxintsq*(1.-oxint/2)
			  sx( 2) = onesixth*xintsq*xint
			end select        

			select case(noz)
			 case(0)
			  sz( 0) = 1.
			 case(1)
			  sz( 0) = 1.-zint
			  sz( 1) = zint
			 case(2)
			  zintsq = zint*zint
			  sz(-1) = 0.5*(0.5-zint)**2
			  sz( 0) = 0.75-zintsq
			  sz( 1) = 0.5*(0.5+zint)**2
			 case(3)
			  ozint = 1.-zint
			  zintsq = zint*zint
			  ozintsq = ozint*ozint
			  sz(-1) = onesixth*ozintsq*ozint
			  sz( 0) = twothird-zintsq*(1.-zint/2)
			  sz( 1) = twothird-ozintsq*(1.-ozint/2)
			  sz( 2) = onesixth*zintsq*zint
			end select        

			 do ll = izmin, izmax
				do jj = ixmin, ixmax

				  rhoold(jold+jj,0,lold+ll) = rhoold(jold+jj,0,lold+ll) + sxold(jj)*szold(ll)*wq

				end do
			end do
		  end do
		end do

	  return
	end subroutine pxr_depose_rhoold_n_2dxz


    !!! --- Computes field divergence
    SUBROUTINE calc_field_div(divee, eex, eey, eez, nx, ny, nz, nxguard, nyguard, nzguard, dx, dy, dz)
        IMPLICIT NONE
        INTEGER(idp) ::  j,k,l
        INTEGER(idp) :: nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: eex,eey,eez
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: divee
        REAL(num) :: dx, dy, dz, invdx, invdy, invdz

        invdx=1.0_num/dx
        invdy=1.0_num/dy
        invdz=1.0_num/dz
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,j,k)
        !$OMP DO COLLAPSE(3)
        DO l = 0, nz
            DO k = 0, ny
                DO j = 0, nx
                    divee(j,k,l) = invdx*(eex(j,k,l)-eex(j-1,k,l))+ &
                                invdy*(eey(j,k,l)-eey(j,k-1,l))+invdz*(eez(j,k,l)-eez(j,k,l-1))
                END DO
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE calc_field_div
    
    SUBROUTINE init_diags
        USE shared_data
        IMPLICIT NONE    
    
        IF (rank.eq.0) THEN 
          WRITE(0,*) 
          WRITE(0,*) ' Initilization of the diags' 
        ENDIF
    
        ! Creation of the result folder
        IF (rank.eq.0) CALL system('mkdir RESULTS')
        IF (rank.eq.0) CALL system('rm RESULTS/*')        
        ! Initialization of the temporal diags 
        CALL init_temp_diags
        
        ! Init time statistics
        CALL init_time_stat_output
        
        IF (rank.eq.0) WRITE(0,*)
        
    END SUBROUTINE
    
    !!! --- Init temporal diags
    SUBROUTINE init_temp_diags
        USE output_data
        USE particle_properties
        USE shared_data
        USE params
        
        IMPLICIT NONE
        
        INTEGER :: i
        
        IF (temdiag_frequency.gt.0) THEN
        
        i = 1
        temdiag_nb_part = 0
        temdiag_nb_field = 0
        temdiag_nb = 0
        
        ! Determine the number of diags for the particles
        ! Kinetic energy
        if (temdiag_act_list(1).gt.0) then
          temdiag_i_list(1) = i
          temdiag_nb_values(1) = nspecies
          temdiag_name_list(1) = 'kinE'
          i = i + nspecies
          temdiag_nb_part = temdiag_nb_part + nspecies
          temdiag_nb = temdiag_nb + 1
        end if
        ! Determine the number of diags for the fields
        ! Ex energy
        if (temdiag_act_list(2).gt.0) then        
          temdiag_i_list(2) = i
          temdiag_name_list(2) = 'exE'
          temdiag_nb_values(2) = 1
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1       
          i = i+1   
        end if
        ! Ey energy
        if (temdiag_act_list(3).gt.0) then        
          temdiag_i_list(3) = i
          temdiag_name_list(3) = 'eyE'
          temdiag_nb_values(3) = 1
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1    
          i = i+1      
        end if   
        ! Ez energy
        if (temdiag_act_list(4).gt.0) then        
          temdiag_i_list(4) = i
          temdiag_name_list(4) = 'ezE'   
          temdiag_nb_values(4) = 1       
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1 
          i = i+1         
        end if 
        ! Bx energy
        if (temdiag_act_list(5).gt.0) then        
          temdiag_i_list(5) = i
          temdiag_name_list(5) = 'bxE'
          temdiag_nb_values(5) = 1
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1       
          i = i+1   
        end if
        ! By energy
        if (temdiag_act_list(6).gt.0) then        
          temdiag_i_list(6) = i
          temdiag_name_list(6) = 'byE'
          temdiag_nb_values(6) = 1
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1       
          i = i+1   
        end if
        ! Bz energy
        if (temdiag_act_list(7).gt.0) then        
          temdiag_i_list(7) = i
          temdiag_name_list(7) = 'bzE'
          temdiag_nb_values(7) = 1
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1       
          i = i+1   
        end if
        ! Norn divE*eps0 - rho
        if (temdiag_act_list(8).gt.0) then        
          temdiag_i_list(8) = i
          temdiag_name_list(8) = 'divE-rho'
          temdiag_nb_values(8) = 1
          temdiag_nb_field = temdiag_nb_field + 1
          temdiag_nb = temdiag_nb + 1       
          i = i+1   
        end if        
                                
        temdiag_totvalues = i - 1
        
        if (temdiag_nb.gt.0) then
          if (rank.eq.0) then 
            write(0,'(" Initialization of the temporal diagnostics")')
          end if
        end if

        ! Each mpi task will write in a given file according to their rank
        IF (nproc.ge.temdiag_nb) then
          IF ((rank.ge.0).and.(rank.le.temdiag_nb)) then
            if (temdiag_act_list(rank+1).gt.0) then
              write(0,'(" Rank ",I3,", creation of the file ",A30)') rank,"./RESULTS/"//trim(adjustl(temdiag_name_list(rank+1)))      
              if (temdiag_format.eq.1) then
                open(unit=42,file="./RESULTS/"//trim(adjustl(temdiag_name_list(rank+1))),ACTION='write')
                write(42,*) temdiag_nb_values(rank+1), temdiag_frequency*dt
              else
                open(unit=42,file="./RESULTS/"//trim(adjustl(temdiag_name_list(rank+1))),FORM="unformatted",ACCESS='stream')
                write(42) temdiag_nb_values(rank+1), temdiag_frequency*dt
              end if
              
            end if
          ENDIF
        ! If there is not enough proc, only the first one deals with it
        else
          IF (rank.eq.0) then
            do i= 1,temdiag_nb
              if (temdiag_act_list(i) .gt.0) then
                write(0,'(" Creation of the file ",A30)') "./RESULTS/"//trim(adjustl(temdiag_name_list(i)))     
                IF (temdiag_format.eq.1) THEN 
                  open(unit=42+i,file="./RESULTS/"//trim(adjustl(temdiag_name_list(i))),ACTION='write')    
                  write(42+i,*) temdiag_nb_values(i),temdiag_frequency*dt
                else
                  open(unit=42+i,file="./RESULTS/"//trim(adjustl(temdiag_name_list(i))),FORM="unformatted",ACCESS='stream')
                  write(42+i) temdiag_nb_values(i),temdiag_frequency*dt
                end if          
              end if
            end do
          end if
        end if
            
        CALL MPI_BARRIER(comm,errcode) 
        
        ENDIF    
    END SUBROUTINE

    ! _____________________________________________________
    SUBROUTINE init_time_stat_output
    ! Initialize outputs of the time statistics
    !
    ! _____________________________________________________
      USE time_stat
      USE shared_data
      USE params
      IMPLICIT NONE 
      
      INTEGER :: nb_timestat
         
      IF (timestat_activated.gt.0) THEN
      
        nb_timestat = 11
      
        OPEN(unit=41,file="./RESULTS/time_stat",FORM="unformatted",ACCESS='stream')
        WRITE(41) nb_timestat
        WRITE(41) timestat_period*dt
        
        IF (rank.eq.0) WRITE(0,*) ' Initialization of the time statistics output: ',nb_timestat
        
        itimestat = 1 ! index in the buffer
        
        ALLOCATE(buffer_timestat(nb_timestat,nbuffertimestat)) ! Buffer before output
      
      ELSE
        timestat_period = 0       
      ENDIF
    
    END SUBROUTINE

    !!! --- Determine the local kinetic energy for the species ispecies
    SUBROUTINE get_loc_kinetic_energy(ispecies,kinetic_energy_loc)
        USE particle_tilemodule
        USE particle_speciesmodule
        USE tile_params
        USE particles
        USE mpi_derived_types
        USE shared_data
        USE tiling
        IMPLICIT NONE
        INTEGER(idp) :: ispecies 
        
        REAL(num) :: kinetic_energy_loc
        INTEGER(idp) :: ix,iy,iz,np,ip,n
        TYPE(particle_tile), POINTER :: curr_tile
        TYPE(particle_species), POINTER :: curr
        INTEGER(idp),parameter :: LVEC2=16
        REAL(num)                          :: partgam
        REAL(num),dimension(:),allocatable :: gaminv

        kinetic_energy_loc = 0
        
        ! current species    
        curr=>species_parray(ispecies)
            
        ! Loop over the tiles
        !$OMP PARALLEL &
        !$OMP DEFAULT(NONE) &
        !$OMP SHARED(curr, ntilex, ntiley, ntilez) &
        !$OMP PRIVATE(ix,iy,iz,ispecies,curr_tile,ip,np,n,gaminv,partgam) &
        !$OMP reduction(+:kinetic_energy_loc)
                        
        !$OMP DO COLLAPSE(3) SCHEDULE(runtime) 
        DO iz=1, ntilez
            DO iy=1, ntiley
                DO ix=1, ntilex
                
                    curr_tile=>curr%array_of_tiles(ix,iy,iz)
                    np = curr_tile%np_tile(1)
                    
                    allocate(gaminv(np))                  
                    gaminv(1:np) = curr_tile%part_gaminv(1:np)
                    
                    !write(0,*) " Tile total kinetic energy for tile:",ix,iy,iz
                    !write(0,*) " Number of particles:",np
                    !write(0,*) " partx:",curr_tile%part_x(1:10) 
                    !write(0,*) " gaminv:",1./curr_tile%part_gaminv(1:10)
                
                    ! Loop over the particles
                    DO ip=1,np,LVEC2
#if defined __INTEL_COMPILER 
                        !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
                        !IBM* ALIGN(64,gaminv)
#endif 
#if defined _OPENMP && _OPENMP>=201307
                        !$OMP SIMD SAFELEN(LVEC2) 
#elif defined __IBMBGQ__
			            !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			            !$DIR SIMD 
#endif 
                        DO n=1,MIN(LVEC2,np-ip+1)
            
                            partgam = 1./gaminv(ip+(n-1)) 
                            
                            kinetic_energy_loc = kinetic_energy_loc + (partgam-1.)*curr_tile%pid(ip+(n-1),wpid)
                            
                        end do                       
                    
                    End do
                       
                    deallocate(gaminv)                      
                    !write(0,*) " Total Local kinetic energy",total_kinetic_energy_loc                  
                    
                End do
            End do 
        End do
        !$OMP END DO
        !$OMP END PARALLEL
    
    end subroutine

    !!! --- Determine the total kinetic energy for the species ispecies
    SUBROUTINE get_kinetic_energy(ispecies,total_kinetic_energy)
        USE particle_tilemodule
        USE particle_speciesmodule
        USE tile_params
        USE particles
        USE mpi_derived_types
        USE shared_data
        USE tiling
        IMPLICIT NONE
        INTEGER(idp) :: ispecies 
        REAL(num) :: total_kinetic_energy
        
        REAL(num) :: total_kinetic_energy_loc 
        INTEGER(idp) :: ix,iy,iz,np,ip,n
        TYPE(particle_tile), POINTER :: curr_tile
        TYPE(particle_species), POINTER :: curr
        INTEGER(idp),parameter :: LVEC2=16
        REAL(num) :: partgam
        
        ! current species    
        curr=>species_parray(ispecies)
            
        ! Loop over the tiles
        !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
        !$OMP SHARED(curr, ntilex, ntiley, ntilez) &
        !$OMP PRIVATE(ix,iy,iz,n,ispecies,curr_tile,ip,np,partgam) &
        !$OMP reduction(+:total_kinetic_energy_loc)
        DO iz=1, ntilez
            DO iy=1, ntiley
                DO ix=1, ntilex
                
                    curr_tile=>curr%array_of_tiles(ix,iy,iz)
                    np = curr_tile%np_tile(1)
                    
                    !write(0,*) " Tile total kinetic energy for tile:",ix,iy,iz
                    !write(0,*) " Number of particles:",np
                    !write(0,*) " partx:",curr_tile%part_x(1:10) 
                    !write(0,*) " gaminv:",1./curr_tile%part_gaminv(1:10)
                
                    ! Loop over the particles (cache)
                    DO ip=1,np,LVEC2
                        !!DIR$ ASSUME_ALIGNED partgaminv:64
#if defined _OPENMP && _OPENMP>=201307
                        !$OMP SIMD SAFELEN(LVEC2) 
#elif defined __IBMBGQ__
			            !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			            !$DIR SIMD 
#endif 
                        DO n=1,MIN(LVEC2,np-ip+1)
            
                            partgam = 1./curr_tile%part_gaminv(ip+(n-1)) 
                            total_kinetic_energy_loc = total_kinetic_energy_loc + (partgam-1.)*curr_tile%pid(ip+(n-1),wpid)
                            
                        end do                       
                    
                    End do
                      
                    !write(0,*) " Total Local kinetic energy",total_kinetic_energy_loc                  
                    
                End do
            End do 
        End do
        !$OMP END PARALLEL DO
        
        ! All MPI reduction
        call MPI_ALLREDUCE(total_kinetic_energy_loc,total_kinetic_energy,1_isp,mpidbl,MPI_SUM,comm,errcode)
        
        !print*,total_kinetic_energy
    
    END SUBROUTINE get_kinetic_energy

    !!! --- Determine the local field energy for the given field
    SUBROUTINE get_loc_field_energy(field,nx2,ny2,nz2,dx2,dy2,dz2,nxguard,nyguard,nzguard,field_energy)
        USE constants
        IMPLICIT NONE
        INTEGER(idp)     :: nx2,ny2,nz2
        INTEGER(idp)     :: nxguard,nyguard,nzguard        
        REAL(num)                   :: field_energy
        INTEGER(idp)                :: j,k,l
        REAL(num)                   :: dx2,dy2,dz2
        REAL(num), dimension(-nxguard:nx2+nxguard,-nyguard:ny2+nyguard,-nzguard:nz2+nzguard), intent(in) :: field        
        field_energy = 0
 
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j)
        !$OMP DO COLLAPSE(3) REDUCTION(+:field_energy)

        do l = 1, nz2
            do k = 1, ny2
                do j = 1, nx2
                
                  field_energy = field_energy + field(j,k,l)*field(j,k,l)*0.5
                
                end do
            end do
        end do
      !$OMP END DO
      !$OMP END PARALLEL   
      
      field_energy = field_energy*dx2*dy2*dz2
    
    END SUBROUTINE
    
    !!! --- Determine the total field energy for the given field
    SUBROUTINE get_field_energy(field,nx2,ny2,nz2,dx2,dy2,dz2,nxguard,nyguard,nzguard,field_energy)
        USE constants
        USE mpi_derived_types
        USE mpi_type_constants
        USE shared_data
        IMPLICIT NONE
        INTEGER(idp)                :: nx2,ny2,nz2
        INTEGER(idp)                :: nxguard,nyguard,nzguard
        REAL(num),dimension(-nxguard:nx2+nxguard,-nyguard:ny2+nyguard,-nzguard:nz2+nzguard) :: field
        REAL(num)                   :: field_energy,field_energy_loc
        INTEGER(idp)                :: j,k,l
        REAL(num)                   :: dx2,dy2,dz2

        field_energy = 0
        field_energy_loc = 0
 
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j)
        !$OMP DO COLLAPSE(3) REDUCTION(+:field_energy_loc)
        do l = 1, nz2
            do k = 1, ny2
                do j = 1, nx2
                
                  field_energy_loc = field_energy_loc + field(j,k,l)*field(j,k,l)*0.5
                
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL   
      
        ! All MPI reduction
        call MPI_ALLREDUCE(field_energy_loc,field_energy,1_isp,mpidbl,MPI_SUM,comm,errcode)
        
        field_energy = field_energy*dx2*dy2*dz2
                   
    END SUBROUTINE

    !!! --- Compute norm of dF/dt = divE -rho/eps0 (parallel function)
    SUBROUTINE get_loc_norm_divErho(divee2,rho2,nx2, ny2, nz2, nxguard, nyguard, nzguard,norm_loc)
        USE mpi_derived_types
        USE mpi_type_constants
        USE shared_data
        USE constants
        IMPLICIT NONE  

        INTEGER(idp)                :: j,k,l
        INTEGER(idp),intent(in)     :: nx2,ny2,nz2,nxguard,nyguard,nzguard        
        REAL(num), dimension(-nxguard:nx2+nxguard,-nyguard:ny2+nyguard,-nzguard:nz2+nzguard),intent(in) :: divee2,rho2  
        REAL(num), intent(out)      :: norm_loc
        
        norm_loc = 0

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j)
        !$OMP DO COLLAPSE(3) REDUCTION(+:norm_loc)
        do l = 1, nz2
            do k = 1, ny2
                do j = 1, nx2

                    norm_loc = norm_loc + (divee2(j,k,l)*eps0 - rho2(j,k,l))**2

                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
    END SUBROUTINE


    !!! --- Compute norm of dF/dt = divE -rho/eps0 (parallel function)
    SUBROUTINE get_norm_divErho(divee2,rho2,nx2, ny2, nz2, nxguard, nyguard, nzguard,norm)
        USE mpi_derived_types
        USE mpi_type_constants
        USE shared_data
        USE constants
        IMPLICIT NONE  

        INTEGER(idp)                :: j,k,l
        INTEGER(idp)                :: nx2,ny2,nz2,nxguard,nyguard,nzguard        
        REAL(num), dimension(-nxguard:nx2+nxguard,-nyguard:ny2+nyguard,-nzguard:nz2+nzguard),intent(in) :: divee2,rho2  
        REAL(num)                   :: norm_loc
        REAL(num), intent(out)      :: norm
        
        norm_loc = 0
        norm = 0

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,k,j)
        !$OMP DO COLLAPSE(3) REDUCTION(+:norm_loc)
        do l = 1, nz2
            do k = 1, ny2
                do j = 1, nx2

                    norm_loc = norm_loc + (divee2(j,k,l)*eps0 - rho2(j,k,l))**2

                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL  

        ! All MPI reduction
        call MPI_ALLREDUCE(norm_loc,norm,1_isp,mpidbl,MPI_SUM,comm,errcode)
        
        norm = sqrt(norm)
        
    END SUBROUTINE
END MODULE diagnostics

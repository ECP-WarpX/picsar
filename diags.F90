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
        IMPLICIT NONE
        INTEGER :: ispecies, ix, iy, iz, count
        INTEGER :: jmin, jmax, kmin, kmax, lmin, lmax
        TYPE(particle_species), POINTER :: curr
        TYPE(particle_tile), POINTER :: curr_tile

        ! - Computes electric field divergence on grid at n+1
        dive=0.0_num
        CALL calc_field_div(dive, ex, ey, ez, nx, ny, nz, nxguards, nyguards, nzguards, dx, dy, dz)

        ! - Computes total charge density
        rho=0.0_num
        DO ispecies=1, nspecies! LOOP ON SPECIES
            curr => species_parray(ispecies)
            DO iz=1, ntilez ! LOOP ON TILES
                DO iy=1, ntiley
                    DO ix=1, ntilex
                        curr_tile=>curr%array_of_tiles(ix,iy,iz)
                        count= curr_tile%np_tile
                        curr_tile%rho_tile=0.0_num
                        jmin=curr_tile%nx_tile_min-nxjguards
                        jmax=curr_tile%nx_tile_max+nxjguards
                        kmin=curr_tile%ny_tile_min-nyjguards
                        kmax=curr_tile%ny_tile_max+nyjguards
                        lmin=curr_tile%nz_tile_min-nzjguards
                        lmax=curr_tile%nz_tile_max+nzjguards
                        ! Depose charge in rho_tile
                        CALL depose_rho_scalar_1_1_1(curr_tile%rho_tile, count,curr_tile%part_x(1:count), &
                             curr_tile%part_y(1:count),curr_tile%part_z(1:count),              &
                             curr_tile%weight(1:count), curr%charge,curr_tile%x_grid_tile_min, &
                             curr_tile%y_grid_tile_min, curr_tile%z_grid_tile_min,dx,dy,dz,    &
                             curr_tile%nx_cells_tile,curr_tile%ny_cells_tile,                  &
                             curr_tile%nz_cells_tile,nxjguards,nyjguards,nzjguards)
                        ! Reduce rho_tile in rho
                        rho(jmin:jmax,kmin:kmax,lmin:lmax) = rho(jmin:jmax,kmin:kmax,lmin:lmax)+ curr_tile%rho_tile
                    END DO
                END DO
            END DO !END LOOP ON TILES
        END DO !END LOOP ON SPECIES
        CALL charge_bcs

    END SUBROUTINE calc_diags

    !!! --- Order 1 3D scalar charge deposition routine
    !!! This version does not vectorize on SIMD architectures
    SUBROUTINE depose_rho_scalar_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
        USE constants
        IMPLICIT NONE
        INTEGER :: np,nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), INTENT(IN OUT) :: rho
        REAL(num) :: xp(np), yp(np), zp(np), w(np)
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
        REAL(num) :: x,y,z,wq,invvol
        REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1)
        REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
        INTEGER :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
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

    !!! --- General charge deposition routine (Warning: Highly unoptimized routine)
    !!! Computes charge density on grid at arbitrary orders nox, noy and noz
    SUBROUTINE depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz, &
                        l_particles_weight, l4symtry)
        IMPLICIT NONE
        INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
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
        INTEGER :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
   
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
    END SUBROUTINE depose_rho_n

    !!! --- Computes field divergence
    SUBROUTINE calc_field_div(divee, eex, eey, eez, nx, ny, nz, nxguard, nyguard, nzguard, dx, dy, dz)
        IMPLICIT NONE
        INTEGER ::  j,k,l
        INTEGER :: nx,ny,nz,nxguard,nyguard,nzguard
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: eex,eey,eez
        REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: divee
        REAL(num) :: dx, dy, dz, invdx, invdy, invdz

        invdx=1.0_num/dx
        invdy=1.0_num/dy
        invdz=1.0_num/dz

        DO l = 0, nz
            DO k = 0, ny
                DO j = 0, nx
                    divee(j,k,l) = invdx*(eex(j,k,l)-eex(j-1,k,l))+ &
                                invdy*(eey(j,k,l)-eey(j,k-1,l))+invdz*(eez(j,k,l)-eez(j,k,l-1))
                END DO
            END DO
        END DO

    END SUBROUTINE calc_field_div

END MODULE diagnostics
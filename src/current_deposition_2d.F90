! ________________________________________________________________________________________
! 
! CURRENT_DEPOSITION_2D.F90
!
! List of subroutines
! - pxrdepose_currents_on_grid_jxjyjz_2d
! ________________________________________________________________________________________

! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d
! 
! Main subroutine for the current deposition called in submain for the 2d
! ________________________________________________________________________________________
  USE fields
  USE shared_data
  USE params
  USE time_stat
#if defined(PROFILING) && PROFILING==2      
  USE ITT_SDE_FORTRAN                       
#endif                                   
  IMPLICIT NONE 
  
  ! __ Parameter declaration __________________________________________________
  REAL(num) :: tdeb, tend
  
  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE

    SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & 
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)!#do not parse
      USE omp_lib
      USE constants
      implicit none
  
      integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(inout) :: jx,jy,jz
      real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
      real(num)                             :: q,dt,dx,dz,xmin,zmin
      logical(idp)                          :: l_particles_weight,l4symtry,l_2drz
      real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
    END SUBROUTINE

    SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_2_2(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & 
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)!#do not parse
      USE omp_lib
      USE constants
      implicit none
  
      integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(inout) :: jx,jy,jz
      real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
      real(num)                             :: q,dt,dx,dz,xmin,zmin
      logical(idp)                          :: l_particles_weight,l4symtry,l_2drz
      real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
    END SUBROUTINE

    SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & 
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)!#do not parse
      USE omp_lib
      USE constants
      implicit none
  
      integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(inout) :: jx,jy,jz
      real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
      real(num)                             :: q,dt,dx,dz,xmin,zmin
      logical(idp)                          :: l_particles_weight,l4symtry,l_2drz
      real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
    END SUBROUTINE

  subroutine pxr_depose_jxjyjz_esirkepov2d_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                   dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                   nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose) !#do not parse
     use constants
     implicit none
     integer(idp)                           :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
     real(num), dimension(np)               :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
     real(num)                              :: q,dt,dx,dz,xmin,zmin
     logical(idp)                           :: l_particles_weight,l4symtry,l_2drz     
     real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  END SUBROUTINE
    
  END INTERFACE
  ! ___________________________________________________________________________

! For debugging    
#if defined(DEBUG)
  WRITE(0,*) "Depose_currents_on_grid: start"
#endif

  ! For time statistics
  tdeb=MPI_WTIME()

! For profiling with Vtune/SDE
#if PROFILING==2              
  CALL start_collection()     
#endif                        

  jx = 0.0_num
  jy = 0.0_num
  jz = 0.0_num

  ! __ Current deposition ________________________________________________________________

  ! _______________________________________________________
  ! Esirkepov OpenMP/tiling version non-vectorized but more optimized than the general order subroutine
  IF (currdepo.EQ.1) THEN

    ! Order 1
    IF ((nox.eq.1).AND.(noz.eq.1)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_1_1,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    ! Order 2
    ELSE IF ((nox.eq.2).AND.(noz.eq.2)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_2_2,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    ! Order 3
    ELSE IF ((nox.eq.3).AND.(noz.eq.3)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_3_3,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    ! Order n
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)
	  ENDIF

  ! _______________________________________________________
	! Default - Esirkepov parallel version with OPENMP/tiling and optimizations
  ELSE

    ! Order 1
    IF ((nox.eq.1).AND.(noz.eq.1)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_1_1,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    ! Order 2
    ELSE IF ((nox.eq.2).AND.(noz.eq.2)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_2_2,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    ! Order 3
    ELSE IF ((nox.eq.3).AND.(noz.eq.3)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_3_3,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    ! Order n
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)
	  ENDIF
	    
  ENDIF

! Stop Vtune/SDE analysis
#if PROFILING==2                     
  CALL stop_collection()            
#endif                               

  ! For time statistics
  tend = MPI_WTIME()
  localtimes(3)=localtimes(3)+(tend-tdeb)
  
! For debugging   
#if defined(DEBUG)
  WRITE(0,*) "Depose_current_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d


! ________________________________________________________________________________________
subroutine pxr_depose_jxjyjz_esirkepov2d_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
! 2D Current deposition esirkepov n order (from 0 to 3)
! This subroutine is adapted from the version of WARP
!
!
! Input parameters:
! - jx,jy,jz: current arrays
! - np: number of particles
! - xp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse of the gamma factor
! - w: weight
! - q: charge
! - xmin,zmin: minimal boundaries of the tile
! - dt, dx, dz: time and space discretization
! - nx,nz: tile grid size
! - nxguard, nzguard: guard cell numbers
! - nox, noz: shape factor order (useless here bur kept for common interface)
! - l_particles_weight: to take into account the particle weight
! - l4symtry (useless here bur kept for common interface)
! - l_2drz  (useless here bur kept for common interface)
! - type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________
   use constants
   implicit none
   integer(idp) :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
   real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
   real(num), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(num) :: q,dt,dx,dz,xmin,zmin
   logical(idp) :: l_particles_weight,l4symtry,l_2drz

   real(num) :: dxi,dzi,dtsdx,dtsdz,xint,yint,zint
   real(num),dimension(:,:), allocatable :: sdx,sdz
   real(num) :: xold,yold,zold,rold,xmid,zmid,x,y,z,r,c,s,wq,wqx,wqz, &
                   tmp,vx,vy,vz,dts2dx,dts2dz, &
                   s1x,s2x,s1z,s2z,invvol,invdtdx,invdtdz, &
                   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq, &
                   dtsdx0,dtsdz0,dts2dx0,dts2dz0
   real(num), parameter :: onesixth=1./6.
   real(num), parameter :: twothird=2./3.
   real(num), dimension(:), allocatable :: sx, sx0, dsx, sz, sz0, dsz
   integer(idp) :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc, &
                   ixmin, ixmax, izmin, izmax, icell, ncells, ndtodx, ndtodz, &
                   xl,xu,zl,zu

    ndtodx = int(clight*dt/dx)
    ndtodz = int(clight*dt/dz)
    xl = -int(nox/2)-1-ndtodx
    xu = int((nox+1)/2)+1+ndtodx
    zl = -int(noz/2)-1-ndtodz
    zu = int((noz+1)/2)+1+ndtodz
    allocate(sdx(xl:xu,zl:zu),sdz(xl:xu,zl:zu))
    allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
    allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

    sx0=0.;sz0=0.
    sdx=0.;sdz=0.
      
    ! Davoine method : limited to order 1 in r
    if (type_rz_depose==2) then
       nox = 1
    endif

      dxi = 1./dx
      dzi = 1./dz
      invvol = 1./(dx*dz)
      dtsdx0 = dt*dxi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dz0 = 0.5*dtsdz0
      invdtdx = 1./(dt*dz)
      invdtdz = 1./(dt*dx)

      do ip=1,np
      
        ! --- computes current position in grid units
        x = xp(ip)
        if (l_2drz) then
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-10) then
            c = x/r 
            s = y/r
          else
            c = 1.
            s = 0.
          end if
          x = r
        end if
        x=x*dxi
        z = zp(ip)*dzi
          
        ! --- computes velocity
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)

        ! --- computes old position in grid units
        if (l_2drz) then
          xold = xp(ip)-dt*vx
          yold = yp(ip)-dt*vy
          rold = sqrt(xold*xold+yold*yold)
          xold=rold*dxi
          vy = -vx*s+vy*c
          vx = (x-xold)/dtsdx0
        else
          xold=x-dtsdx0*vx
        end if
        zold=z-dtsdz0*vz
 
        ! --- applies 4-fold symmetry
        if (l4symtry) then
          x=abs(x)
          xold=abs(xold)
          vx = (x-xold)/dtsdx0
        end if

        ! --- sets positions relative to grid  start
        x = x-xmin*dxi
        z = z-zmin*dzi
        xold = xold-xmin*dxi
        zold = zold-zmin*dzi
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1!+max( int(abs(x-xold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        z=zold
        
        do icell = 1,ncells

        xold = x
        zold = z

        x = x+dtsdx*vx
        z = z+dtsdz*vz

        ! --- computes particles "weights"
        if (l_particles_weight) then
           wq=q*w(ip)
        else
           wq=q*w(1)
        end if
        wqx = wq*invdtdx
        wqz = wq*invdtdz

        ! --- finds node of cell containing particles for current positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
          iixp0=nint(x)
        else
          iixp0=floor(x)
        end if
        if (noz==2*(noz/2)) then
          ikxp0=nint(z)
        else
          ikxp0=floor(z)
        end if

        ! --- computes distance between particle and node for current positions
        xint=x-iixp0
        zint=z-ikxp0

        ! --- computes coefficients for node centered quantities
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sx0(0) = 1. - xint  + 1./(4*iixp0+2)*( -xint + xint**2 )
           sx0(1) = 1. - sx0(0)
        else! Standard method, canonical shapes in r 
           select case(nox)
           case(0)
              sx0( 0) = 1.
           case(1)
              sx0( 0) = 1.-xint
              sx0( 1) = xint
           case(2)
              xintsq = xint*xint
              sx0(-1) = 0.5*(0.5-xint)**2
              sx0( 0) = 0.75-xintsq
              sx0( 1) = 0.5*(0.5+xint)**2
           case(3)
              oxint = 1.-xint
              xintsq = xint*xint
              oxintsq = oxint*oxint
              sx0(-1) = onesixth*oxintsq*oxint
              sx0( 0) = twothird-xintsq*(1.-xint/2)
              sx0( 1) = twothird-oxintsq*(1.-oxint/2)
              sx0( 2) = onesixth*xintsq*xint
           end select
        endif

        select case(noz)
        case(0)
           sz0( 0) = 1.
        case(1)
           sz0( 0) = 1.-zint
           sz0( 1) = zint
        case(2)
           zintsq = zint*zint
           sz0(-1) = 0.5*(0.5-zint)**2
           sz0( 0) = 0.75-zintsq
           sz0( 1) = 0.5*(0.5+zint)**2
        case(3)
           ozint = 1.-zint
           zintsq = zint*zint
           ozintsq = ozint*ozint
           sz0(-1) = onesixth*ozintsq*ozint
           sz0( 0) = twothird-zintsq*(1.-zint/2)
           sz0( 1) = twothird-ozintsq*(1.-ozint/2)
           sz0( 2) = onesixth*zintsq*zint
        end select

        ! --- finds node of cell containing particles for old positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
           iixp=nint(xold)
        else
           iixp=floor(xold)
        end if
        if (noz==2*(noz/2)) then
           ikxp=nint(zold)
        else
           ikxp=floor(zold)
        end if

        ! --- computes distance between particle and node for old positions
        xint = xold-iixp
        zint = zold-ikxp

        ! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diz = ikxp-ikxp0

        ! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.;sz=0.

        ! --- computes coefficients for quantities centered between nodes
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sx(0) = 1. - xint  + 1./(4*iixp+2)*( -xint + xint**2 )
           sx(1) = 1. - sx(0)
        else! Standard method, canonical shapes in r 
           select case(nox)
           case(0)
              sx( 0+dix) = 1.
           case(1)
              sx( 0+dix) = 1.-xint
              sx( 1+dix) = xint
           case(2)
              xintsq = xint*xint
              sx(-1+dix) = 0.5*(0.5-xint)**2
              sx( 0+dix) = 0.75-xintsq
              sx( 1+dix) = 0.5*(0.5+xint)**2
           case(3)
              oxint = 1.-xint
              xintsq = xint*xint
              oxintsq = oxint*oxint
              sx(-1+dix) = onesixth*oxintsq*oxint
              sx( 0+dix) = twothird-xintsq*(1.-xint/2)
              sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
              sx( 2+dix) = onesixth*xintsq*xint
           end select
        endif
        
        select case(noz)
        case(0)
           sz( 0+diz) = 1.
        case(1)
           sz( 0+diz) = 1.-zint
           sz( 1+diz) = zint
        case(2)
           zintsq = zint*zint
           sz(-1+diz) = 0.5*(0.5-zint)**2
           sz( 0+diz) = 0.75-zintsq
           sz( 1+diz) = 0.5*(0.5+zint)**2
        case(3)
           ozint = 1.-zint
           zintsq = zint*zint
           ozintsq = ozint*ozint
           sz(-1+diz) = onesixth*ozintsq*ozint
           sz( 0+diz) = twothird-zintsq*(1.-zint/2)
           sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
           sz( 2+diz) = onesixth*zintsq*zint
        end select

        ! --- computes coefficients difference
        dsx = sx - sx0
        dsz = sz - sz0

        ! --- computes min/max positions of current contributions
        ixmin = min(0_idp,dix)-int(nox/2)
        ixmax = max(0_idp,dix)+int((nox+1)/2)
        izmin = min(0_idp,diz)-int(noz/2)
        izmax = max(0_idp,diz)+int((noz+1)/2)
        
        ! --- add current contributions
        ! --- NB : the current is later divided by the cylindrical cell volume in applybc_j
        do k=izmin, izmax
           do i=ixmin, ixmax
              ic = iixp0+i
              kc = ikxp0+k

              ! -- Jx
              if(i<ixmax) then
                 sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
                 if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x
                 jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
              end if
              
              ! -- Jy (2D Esirkepov scheme)
              jy(ic,kc) = jy(ic,kc) + wq*vy*invvol/ncells* &
                   ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )

              ! -- Jz
              if(k<izmax) then
                 sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov
                 if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
                 jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
              end if
           end do
        end do
     end do
    end do
    
    deallocate(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)

  return
end subroutine pxr_depose_jxjyjz_esirkepov2d_n

! ________________________________________________________________________________________
SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
!
! 2D Current deposition with the method of Esirkepov at order 1
! This function is not optimized but provides better performances than 
! using the abitrary order function
!
! Input parameters:
! - jx,jy,jz: current arrays
! - np: number of particles
! - xp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse of the gamma factor
! - w: weight
! - q: charge
! - xmin,zmin: minimal boundaries of the tile
! - dt, dx, dz: time and space discretization
! - nx,nz: tile grid size
! - nxguard, nzguard: guard cell numbers
! - nox, noz: shape factor order (useless here bur kept for common interface)
! - l_particles_weight: to take into account the particle weight
! - l4symtry (useless here bur kept for common interface)
! - l_2drz  (useless here bur kept for common interface)
! - type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none

  ! __ Parameter declaration ________________________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  logical(idp)                          :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), parameter                  :: onethird=1./3.  
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc
  integer(idp)                          :: ixmin, ixmax, izmin, izmax, icell, ndtodx, ndtodz
  integer(idp)                          :: xl,xu,zl,zu

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(-1:2,-1:2),sdz(-1:2,-1:2))
  ALLOCATE(sx(-1:2), sx0(-1:2), dsx(-1:2))
  ALLOCATE(sz(-1:2), sz0(-1:2), dsz(-1:2))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num
  
  DO ip=1,np
  
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
  
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
  
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz
  
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz
  
    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ikxp0=floor(z)
  
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0
  
    ! --- computes coefficients for node centered quantities
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
  
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint
  
    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ikxp=floor(zold)
  
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp
  
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0
  
    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
  
    ! --- computes coefficients for quantities centered between nodes
    sx( 0+dix) = 1.0_num-xint
    sx( 1+dix) = xint
  
    sz( 0+diz) = 1.0_num-zint
    sz( 1+diz) = zint
  
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0
    
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)-0
    ixmax = max(0,dix)+1
    izmin = min(0,diz)-0
    izmax = max(0,diz)+1
    
    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k
        
        ! --- Jx
        IF(i<ixmax) THEN
          sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
          if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x 
          jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
        END IF
        
        ! -- Jy (2D Esirkepov scheme)
        jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &
        ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )
        
        ! --- Jz
        IF(k<izmax) THEN
          sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov&
          if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
          jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
        END IF
      END DO
    END DO

  END DO
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)
  RETURN

END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1

! ________________________________________________________________________________________
SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_2_2(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
!
! 2D Current deposition with the method of Esirkepov at order 2
! This function is not optimized but provides better performances than 
! using the abitrary order function
!
! Input parameters:
! - jx,jy,jz: current arrays
! - np: number of particles
! - xp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse of the gamma factor
! - w: weight
! - q: charge
! - xmin,zmin: minimal boundaries of the tile
! - dt, dx, dz: time and space discretization
! - nx,nz: tile grid size
! - nxguard, nzguard: guard cell numbers
! - nox, noz: shape factor order (useless here bur kept for common interface)
! - l_particles_weight: to take into account the particle weight
! - l4symtry (useless here bur kept for common interface)
! - l_2drz  (useless here bur kept for common interface)
! - type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none
  
  ! __ Parameter declaration ________________________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  logical(idp)                          :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), parameter                  :: onethird=1./3.  
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc
  integer(idp)                          :: ixmin, ixmax, izmin, izmax, icell, ndtodx, ndtodz
  integer(idp)                          :: xl,xu,zl,zu

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(-2:2,-2:2),sdz(-2:2,-2:2))
  ALLOCATE(sx(-2:2), sx0(-2:2), dsx(-2:2))
  ALLOCATE(sz(-2:2), sz0(-2:2), dsz(-2:2))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num
  
  DO ip=1,np
  
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
  
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
  
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz
  
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz
  
    ! --- finds node of cell containing particles for current positions
    iixp0=nint(x)
    ikxp0=nint(z)
  
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0
  
    ! --- computes coefficients for node centered quantities
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2
  
    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = 0.75_num-zintsq
    sz0( 1) = 0.5_num*(0.5_num+zint)**2
  
    ! --- finds node of cell containing particles for old positions
    iixp=nint(xold)
    ikxp=nint(zold)
  
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp
  
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0
  
    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx(-2)=0.0_num
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sz(-2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
  
    ! --- computes coefficients for quantities centered between nodes
    xintsq = xint*xint
    sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
    sx( 0+dix) = 0.75_num-xintsq
    sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
  
    zintsq = zint*zint
    sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
    sz( 0+diz) = 0.75_num-zintsq
    sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
  
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0
    
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)-1
    ixmax = max(0,dix)+1
    izmin = min(0,diz)-1
    izmax = max(0,diz)+1
    
    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k
        
        ! --- Jx
        IF(i<ixmax) THEN
          sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
          if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x 
          jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
        END IF
        
        ! -- Jy (2D Esirkepov scheme)
        jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &
        ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )
        
        ! --- Jz
        IF(k<izmax) THEN
          sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov&
          if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
          jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
        END IF
      END DO
    END DO

  END DO
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)
  RETURN

End subroutine pxr_depose_jxjyjz_esirkepov2d_2_2

! ________________________________________________________________________________________
subroutine pxr_depose_jxjyjz_esirkepov2d_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
!
! 2D Current deposition with the method of Esirkepov at order 3
! This function is not optimized but provides better performances than 
! using the abitrary order function
!
! Input parameters:
! - jx,jy,jz: current arrays
! - np: number of particles
! - xp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse of the gamma factor
! - w: weight
! - q: charge
! - xmin,zmin: minimal boundaries of the tile
! - dt, dx, dz: time and space discretization
! - nx,nz: tile grid size
! - nxguard, nzguard: guard cell numbers
! - nox, noz: shape factor order (useless here bur kept for common interface)
! - l_particles_weight: to take into account the particle weight
! - l4symtry (useless here bur kept for common interface)
! - l_2drz  (useless here bur kept for common interface)
! - type_rz_depose (useless here bur kept for common interface)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none
  
  ! __ Parameter declaration _______________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  logical(idp)                          :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc
  integer(idp)                          :: ixmin, ixmax, izmin, izmax, icell, ndtodx, ndtodz
  integer(idp)                          :: xl,xu,zl,zu

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(-2:3,-2:3),sdz(-2:3,-2:3))
  ALLOCATE(sx(-2:3), sx0(-2:3), dsx(-2:3))
  ALLOCATE(sz(-2:3), sz0(-2:3), dsz(-2:3))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num

  DO ip=1,np
  
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
  
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
  
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz
  
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz  

    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ikxp0=floor(z)
  
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx0(-1) = onesixth*oxintsq*oxint
    sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx0( 2) = onesixth*xintsq*xint
  
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz0(-1) = onesixth*ozintsq*ozint
    sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz0( 2) = onesixth*zintsq*zint

    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ikxp=floor(zold)

    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp

    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0

    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx(-2)=0.0_num
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sx(3)=0.0_num
    sz(-2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
    sz(3)=0.0_num

    ! --- computes coefficients for quantities centered between nodes
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1+dix) = onesixth*oxintsq*oxint
    sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx( 2+dix) = onesixth*xintsq*xint
  
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1+diz) = onesixth*ozintsq*ozint
    sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz( 2+diz) = onesixth*zintsq*zint

    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0
    
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)-1
    ixmax = max(0,dix)+2
    izmin = min(0,diz)-1
    izmax = max(0,diz)+2

    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k
        
        ! --- Jx
        IF(i<ixmax) THEN
          sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
          if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x 
          jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
        END IF
        
        ! -- Jy (2D Esirkepov scheme)
        jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &
        ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )
        
        ! --- Jz
        IF(k<izmax) THEN
          sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov&
          if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
          jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
        END IF

! __ DEbug _______________________________
!  print*,'sum',sum(jx),sum(jy),sum(jz)
!  print*,'j',jy(ic,kc)
!  print*,'sdx',sdx(i,k),sdz(i,k)
!  print*,'wq',wqx,wqz,wq
!  print*,'s',sz0(k),sx0(i)
!  print*,'dsx',dsx(i),dsz(k)
!  read*        
        
      END DO
    END DO

  END DO
  
  
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)

End subroutine pxr_depose_jxjyjz_esirkepov2d_3_3
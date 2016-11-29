!_________________________________________________________________________________________
!
! FIELD_GATHERING_2D.F90
!
! List of subroutines:
! - geteb2dxz_energy_conserving
!
! - pxr_gete2dxz_n_energy_conserving
! - pxr_getb2dxz_n_energy_conserving
! - pxr_gete2dxz_energy_conserving_3_3
! - pxr_getb2dxz_energy_conserving_3_3
!_________________________________________________________________________________________


!_________________________________________________________________________________________
!> General subroutines for the 2D cartesian field gathering
!> @brief
!
!> @author
!> Mathieu Lobet
!
!> @date
!> 2016
SUBROUTINE geteb2dxz_energy_conserving(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,exg,eyg,ezg,bxg,byg,bzg,l4symtry,l_lower_order_in_v)
!_________________________________________________________________________________________

  USE constants
  USE params
  implicit none
  
  integer(idp)                  :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
  logical(idp), intent(in)      :: l4symtry,l_lower_order_in_v
  real(num), dimension(np)      :: xp,yp,zp,ex,ey,ez,bx,by,bz
  real(num), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg    
  real(num), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
  real(num)                     :: xmin,ymin,zmin,dx,dy,dz
              
  ! ______________________________________________                                  
  ! Arbitrary order, non-optimized subroutines
  IF (fieldgathe.eq.1) THEN
  

    !!! --- Gather electric field on particles
    CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)
  
  ! ________________________________________
  ! Optimized subroutines, default  
  ELSE

  
    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
  
    !!! --- Gather electric field on particles
    CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)
                                      
    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

    !!! --- Gather electric field on particles
    CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_3_3(np,xp,zp,ex,ey,ez,xmin,zmin,   &
                                        dx,dz,nx,nz,nxguard,nzguard, &
                                        exg,eyg,ezg,l_lower_order_in_v)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_3_3(np,xp,zp,bx,by,bz,xmin,zmin,   &
                                        dx,dz,nx,nz,nxguard,nzguard, &
                                        bxg,byg,bzg,l_lower_order_in_v) 
    
    ! Arbitrary order             
    ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
      !!! --- Gather magnetic fields on particles
     CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                 dx,dz,nx,nz,nxguard,nzguard, &
                                 nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)		  
    
    ENDIF				  	    
  ENDIF
END SUBROUTINE geteb2dxz_energy_conserving

!_________________________________________________________________________________________
!> 2D electric field gathering routine 
subroutine pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard, &
                                       nox,noz,exg,eyg,ezg,l4symtry,l_2drz,l_lower_order_in_v)
!_________________________________________________________________________________________                                       
      use constants
      implicit none
      
      integer(idp)             :: np,nx,nz,nox,noz,nxguard,nzguard
      real(num), dimension(np) :: xp,yp,zp,ex,ey,ez
      logical(idp)             :: l4symtry,l_2drz,l_lower_order_in_v
      real(num), dimension(-nxguard:nx+nxguard,1,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(num)                :: xmin,zmin,dx,dz,costheta,sintheta
      integer(idp)             :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                      ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0, l0
      real(num) :: dxi, dzi, x, y, z, r, xint, zint, &
                      xintsq,oxint,zintsq,ozint,oxintsq,ozintsq,signx
      real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(num), dimension(:), allocatable :: sx0,sz0
      real(num), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox+1)/2)-1
      izmin = -int(noz/2)
      izmax =  int((noz+1)/2)-1

      if (l_lower_order_in_v) then
        ixmin0 = -int((nox-1)/2)
        ixmax0 =  int((nox)/2)
        izmin0 = -int((noz-1)/2)
        izmax0 =  int((noz)/2)
      else
        ixmin0 = -int((nox)/2)
        ixmax0 =  int((nox+1)/2)
        izmin0 = -int((noz)/2)
        izmax0 =  int((noz+1)/2)
      end if
      allocate(sx0(ixmin0:ixmax0),sz0(izmin0:izmax0))

      signx = 1.

      do ip=1,np

        if (l_2drz) then
          x = xp(ip)
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-20) then
            costheta=x/r
            sintheta=y/r
          else  
            costheta=1.
            sintheta=0.
          end if
          x = (r-xmin)*dxi
        else
          x = (xp(ip)-xmin)*dxi
        end if

        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
        end if
        
        if (l_lower_order_in_v) then
          if (nox==2*(nox/2)) then
            j=nint(x)
            j0=floor(x-0.5)
          else
            j=floor(x)
            j0=floor(x)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z-0.5)
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
            j0=floor(x-0.5)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z)
          else
            l=floor(z)
            l0=floor(z-0.5)
          end if
        end if
        
        xint=x-j
        zint=z-l

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

        xint=x-0.5-j0
        zint=z-0.5-l0

        if (l_lower_order_in_v) then
        
         if (nox==1) then
          sx0( 0) = 1.
         elseif (nox==2) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         end if

         if (noz==1) then
          sz0( 0) = 1.
         elseif (noz==2) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         end if

        else

         if (nox==1) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==2) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
         end if

         if (noz==1) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==2) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
         end if

        end if

        if (l_2drz) then
       
!          write(0,*) 'field gathering needs to be done for fstype=4 in EM-RZ'
!          stop
          do ll = izmin, izmax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sz(ll)*sx0(jj)*(exg(j0+jj,1,l+ll)*costheta-eyg(j0+jj,1,l+ll)*sintheta)
              ey(ip) = ey(ip) + sz(ll)*sx0(jj)*(exg(j0+jj,1,l+ll)*sintheta+eyg(j0+jj,1,l+ll)*costheta)
            end do
          end do

        else

          do ll = izmin, izmax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j0+jj,1,l+ll)*signx
            end do
          end do

          do ll = izmin, izmax+1
            do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj,1,l+ll)
            end do
          end do

        end if

          do ll = izmin0, izmax0
            do jj = ixmin, ixmax+1
              ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj,1,l0+ll)
            end do
          end do
                     
     end do
     deallocate(sx0,sz0)
     
   return
 end subroutine pxr_gete2dxz_n_energy_conserving


! ________________________________________________________________________________________
!> 2D magnetic field gathering routine 
subroutine pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard, &
                                       nox,noz,bxg,byg,bzg,l4symtry,l_2drz,l_lower_order_in_v)
! ________________________________________________________________________________________
      use constants
      implicit none
      
      integer(idp) :: np,nx,nz,nox,noz,nxguard,nzguard
      real(num), dimension(np) :: xp,yp,zp,bx,by,bz
      logical(idp) :: l4symtry,l_2drz,l_lower_order_in_v
      real(num), dimension(-nxguard:nx+nxguard,1,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(num) :: xmin,zmin,dx,dz
      integer(idp) :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                      ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0, l0
      real(num) :: dxi, dzi, x, y, z, xint, zint, &
                      xintsq,oxint,zintsq,ozint,oxintsq,ozintsq,signx, &
                      r, costheta, sintheta
      real(num), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(num), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(num), dimension(:), allocatable :: sx0,sz0
      real(num), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox+1)/2)-1
      izmin = -int(noz/2)
      izmax =  int((noz+1)/2)-1

      if (l_lower_order_in_v) then
        ixmin0 = -int((nox-1)/2)
        ixmax0 =  int((nox)/2)
        izmin0 = -int((noz-1)/2)
        izmax0 =  int((noz)/2)
      else
        ixmin0 = -int((nox)/2)
        ixmax0 =  int((nox+1)/2)
        izmin0 = -int((noz)/2)
        izmax0 =  int((noz+1)/2)
      end if
      allocate(sx0(ixmin0:ixmax0),sz0(izmin0:izmax0))

      signx = 1.

      sx=0
      sz=0.
      sx0=0.
      sz0=0.

      do ip=1,np

        if (l_2drz) then
          x = xp(ip)
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-20) then
            costheta=x/r
            sintheta=y/r
          else  
            costheta=1.
            sintheta=0.
          end if
          x = (r-xmin)*dxi
        else
          x = (xp(ip)-xmin)*dxi
        end if

        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
        end if

        if (l_lower_order_in_v) then
          if (nox==2*(nox/2)) then
            j=nint(x)
            j0=floor(x-0.5)
          else
            j=floor(x)
            j0=floor(x)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z-0.5)
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
            j0=floor(x-0.5)
          end if
          if (noz==2*(noz/2)) then
            l=nint(z)
            l0=floor(z)
          else
            l=floor(z)
            l0=floor(z-0.5)
          end if
        end if

        xint=x-j
        zint=z-l

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

        xint=x-0.5-j0
        zint=z-0.5-l0

        if (l_lower_order_in_v) then
        
         if (nox==1) then
          sx0( 0) = 1.
         elseif (nox==2) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         end if

         if (noz==1) then
          sz0( 0) = 1.
         elseif (noz==2) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         end if

        else

         if (nox==1) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         elseif (nox==2) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
         end if

         if (noz==1) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         elseif (noz==2) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
         end if

        end if

        if (l_2drz) then

          do ll = izmin0, izmax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sz0(ll)*(bxg(j+jj,1,l0+ll)*costheta-byg(j+jj,1,l0+ll)*sintheta)
              by(ip) = by(ip) + sx(jj)*sz0(ll)*(bxg(j+jj,1,l0+ll)*sintheta+byg(j+jj,1,l0+ll)*costheta)
            end do
          end do

        else

          do ll = izmin0, izmax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj,1,l0+ll)*signx
            end do
          end do

          do ll = izmin0, izmax0
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj,1,l0+ll)
            end do
          end do

        end if

        do ll = izmin, izmax+1
            do jj = ixmin0, ixmax0
              bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j0+jj,1,l+ll)
            end do
        end do
                 
     end do
     deallocate(sx0,sz0)

   return
 end subroutine pxr_getb2dxz_n_energy_conserving



!_________________________________________________________________________________________
!> Field gathering cartesian in 2D for the electric field
!> @brief
!
!> This function is vectorized
!> @details
subroutine pxr_gete2dxz_energy_conserving_3_3(np,xp,zp,ex,ey,ez,xmin,zmin,dx,dz,nx,nz, &
													nxguard,nzguard,exg,eyg,ezg,l_lower_order_in_v)
!                                       
!
! Inputs:
! - np: particle number
! - xp, zp: particle position arrays
! - ex,ey,ez: electric field particle arrays
! - xmin, zmin: tile boundaries
! - dx, dz: space steps
! - nx,nz: space discretization
! - nxguard, nzguard: number of guard cells
! - exg, eyg,ezg: field arrays
! - l_lower_order_in_v:
!
! output:
! - ex,ey,ez: electric field particle arrays
!_________________________________________________________________________________________
	use constants
	implicit none
	
	
	integer(idp)                  :: np,nx,nz,nox,noz,nxguard,nzguard
	real(num), dimension(np)      :: xp,zp,ex,ey,ez
	logical(idp)                  :: l_lower_order_in_v
	real(num), dimension(-nxguard:nx+nxguard,1,-nzguard:nz+nzguard) :: exg,eyg,ezg
	real(num)                     :: xmin,zmin,dx,dz
	integer(isp)                  :: ip, j, l, jj, ll, j0, l0
	real(num)                     :: dxi, dzi, x, z, r, xint, zint, &
																		xintsq,oxint,zintsq,ozint,oxintsq,ozintsq,signx
	real(num), DIMENSION(-1:2)    :: sx,sx0
	real(num), DIMENSION(-1:2)    :: sz,sz0
	real(num), parameter          :: onesixth=1./6.,twothird=2./3.

	dxi = 1./dx
	dzi = 1./dz

	! ___________________________   
	IF (l_lower_order_in_v) THEN
      
#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
		DO ip=1,np
	
			x = (xp(ip)-xmin)*dxi
			z = (zp(ip)-zmin)*dzi

			! Compute index of particle
			j=floor(x)
			j0=floor(x)

			l=floor(z)
			l0=floor(z)

			xint=x-j
			zint=z-l

			! Compute shape factors
			oxint = 1.0_num-xint
			xintsq = xint*xint
			oxintsq = oxint*oxint
			sx(-1) = onesixth*oxintsq*oxint
			sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
			sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
			sx( 2) = onesixth*xintsq*xint

			ozint = 1.0_num-zint
			zintsq = zint*zint
			ozintsq = ozint*ozint
			sz(-1) = onesixth*ozintsq*ozint
			sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
			sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
			sz( 2) = onesixth*zintsq*zint

			xint=x-0.5_num-j0
			zint=z-0.5_num-l0

			xintsq = xint*xint
			sx0(-1) = 0.5_num*(0.5_num-xint)**2
			sx0( 0) = 0.75_num-xintsq
			sx0( 1) = 0.5_num*(0.5_num+xint)**2

			zintsq = zint*zint
			sz0(-1) = 0.5_num*(0.5_num-zint)**2
			sz0( 0) = 0.75_num-zintsq
			sz0( 1) = 0.5_num*(0.5_num+zint)**2

			! Compute Ex on particle
			ex(ip) = ex(ip) + sx0(-1)*sz(-1)*exg(j0-1,1,l-1)
			ex(ip) = ex(ip) + sx0(0)*sz(-1)*exg(j0,1,l-1)
			ex(ip) = ex(ip) + sx0(1)*sz(-1)*exg(j0+1,1,l-1)
			ex(ip) = ex(ip) + sx0(-1)*sz(0)*exg(j0-1,1,l)
			ex(ip) = ex(ip) + sx0(0)*sz(0)*exg(j0,1,l)
			ex(ip) = ex(ip) + sx0(1)*sz(0)*exg(j0+1,1,l)
			ex(ip) = ex(ip) + sx0(-1)*sz(1)*exg(j0-1,1,l+1)
			ex(ip) = ex(ip) + sx0(0)*sz(1)*exg(j0,1,l+1)
			ex(ip) = ex(ip) + sx0(1)*sz(1)*exg(j0+1,1,l+1)
			ex(ip) = ex(ip) + sx0(-1)*sz(2)*exg(j0-1,1,l+2)
			ex(ip) = ex(ip) + sx0(0)*sz(2)*exg(j0,1,l+2)
			ex(ip) = ex(ip) + sx0(1)*sz(2)*exg(j0+1,1,l+2)

			! Compute Ey on particle
			ey(ip) = ey(ip) + sx(-1)*sz(-1)*eyg(j-1,1,l-1)
			ey(ip) = ey(ip) + sx(0)*sz(-1)*eyg(j,1,l-1)
			ey(ip) = ey(ip) + sx(1)*sz(-1)*eyg(j+1,1,l-1)
			ey(ip) = ey(ip) + sx(2)*sz(-1)*eyg(j+2,1,l-1)
			ey(ip) = ey(ip) + sx(-1)*sz(0)*eyg(j-1,1,l)
			ey(ip) = ey(ip) + sx(0)*sz(0)*eyg(j,1,l)
			ey(ip) = ey(ip) + sx(1)*sz(0)*eyg(j+1,1,l)
			ey(ip) = ey(ip) + sx(2)*sz(0)*eyg(j+2,1,l)
			ey(ip) = ey(ip) + sx(-1)*sz(1)*eyg(j-1,1,l+1)
			ey(ip) = ey(ip) + sx(0)*sz(1)*eyg(j,1,l+1)
			ey(ip) = ey(ip) + sx(1)*sz(1)*eyg(j+1,1,l+1)
			ey(ip) = ey(ip) + sx(2)*sz(1)*eyg(j+2,1,l+1)
			ey(ip) = ey(ip) + sx(-1)*sz(2)*eyg(j-1,1,l+2)
			ey(ip) = ey(ip) + sx(0)*sz(2)*eyg(j,1,l+2)
			ey(ip) = ey(ip) + sx(1)*sz(2)*eyg(j+1,1,l+2)
			ey(ip) = ey(ip) + sx(2)*sz(2)*eyg(j+2,1,l+2)

			! Compute Ez on particle
			ez(ip) = ez(ip) + sx(-1)*sz0(-1)*ezg(j-1,1,l0-1)
			ez(ip) = ez(ip) + sx(0)*sz0(-1)*ezg(j,1,l0-1)
			ez(ip) = ez(ip) + sx(1)*sz0(-1)*ezg(j+1,1,l0-1)
			ez(ip) = ez(ip) + sx(2)*sz0(-1)*ezg(j+2,1,l0-1)
			ez(ip) = ez(ip) + sx(-1)*sz0(0)*ezg(j-1,1,l0)
			ez(ip) = ez(ip) + sx(0)*sz0(0)*ezg(j,1,l0)
			ez(ip) = ez(ip) + sx(1)*sz0(0)*ezg(j+1,1,l0)
			ez(ip) = ez(ip) + sx(2)*sz0(0)*ezg(j+2,1,l0)
			ez(ip) = ez(ip) + sx(-1)*sz0(1)*ezg(j-1,1,l0+1)
			ez(ip) = ez(ip) + sx(0)*sz0(1)*ezg(j,1,l0+1)
			ez(ip) = ez(ip) + sx(1)*sz0(1)*ezg(j+1,1,l0+1)
			ez(ip) = ez(ip) + sx(2)*sz0(1)*ezg(j+2,1,l0+1)
		
		END DO
#if defined _OPENMP && _OPENMP>=201307
			!$OMP END SIMD 
#endif        
        
	! ___________________________         
	! l_lower_order_in_v is false      
	ELSE

#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
		DO ip=1,np

			x = (xp(ip)-xmin)*dxi
			z = (zp(ip)-zmin)*dzi

			! Compute index of particle
			j=floor(x)
			j0=floor(x-0.5_num)
			l=floor(z)
			l0=floor(z-0.5_num)
		
			xint=x-j
			zint=z-l

			oxint = 1.0_num-xint
			xintsq = xint*xint
			oxintsq = oxint*oxint
			sx(-1) = onesixth*oxintsq*oxint
			sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
			sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
			sx( 2) = onesixth*xintsq*xint

			ozint = 1.0_num-zint
			zintsq = zint*zint
			ozintsq = ozint*ozint
			sz(-1) = onesixth*ozintsq*ozint
			sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
			sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
			sz( 2) = onesixth*zintsq*zint
		
			xint=x-0.5_num-j0
			zint=z-0.5_num-l0

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

			! Compute Ex on particle
			ex(ip) = ex(ip) + sx0(-1)*sz(-1)*exg(j0-1,1,l-1)
			ex(ip) = ex(ip) + sx0(0)*sz(-1)*exg(j0,1,l-1)
			ex(ip) = ex(ip) + sx0(1)*sz(-1)*exg(j0+1,1,l-1)
			ex(ip) = ex(ip) + sx0(2)*sz(-1)*exg(j0+2,1,l-1)
			ex(ip) = ex(ip) + sx0(-1)*sz(0)*exg(j0-1,1,l)
			ex(ip) = ex(ip) + sx0(0)*sz(0)*exg(j0,1,l)
			ex(ip) = ex(ip) + sx0(1)*sz(0)*exg(j0+1,1,l)
			ex(ip) = ex(ip) + sx0(2)*sz(0)*exg(j0+2,1,l)
			ex(ip) = ex(ip) + sx0(-1)*sz(1)*exg(j0-1,1,l+1)
			ex(ip) = ex(ip) + sx0(0)*sz(1)*exg(j0,1,l+1)
			ex(ip) = ex(ip) + sx0(1)*sz(1)*exg(j0+1,1,l+1)
			ex(ip) = ex(ip) + sx0(2)*sz(1)*exg(j0+2,1,l+1)
			ex(ip) = ex(ip) + sx0(-1)*sz(2)*exg(j0-1,1,l+2)
			ex(ip) = ex(ip) + sx0(0)*sz(2)*exg(j0,1,l+2)
			ex(ip) = ex(ip) + sx0(1)*sz(2)*exg(j0+1,1,l+2)
			ex(ip) = ex(ip) + sx0(2)*sz(2)*exg(j0+2,1,l+2)

			! Compute Ey on particle
			ey(ip) = ey(ip) + sx(-1)*sz(-1)*eyg(j-1,1,l-1)
			ey(ip) = ey(ip) + sx(0)*sz(-1)*eyg(j,1,l-1)
			ey(ip) = ey(ip) + sx(1)*sz(-1)*eyg(j+1,1,l-1)
			ey(ip) = ey(ip) + sx(2)*sz(-1)*eyg(j+2,1,l-1)
			ey(ip) = ey(ip) + sx(-1)*sz(0)*eyg(j-1,1,l)
			ey(ip) = ey(ip) + sx(0)*sz(0)*eyg(j,1,l)
			ey(ip) = ey(ip) + sx(1)*sz(0)*eyg(j+1,1,l)
			ey(ip) = ey(ip) + sx(2)*sz(0)*eyg(j+2,1,l)
			ey(ip) = ey(ip) + sx(-1)*sz(1)*eyg(j-1,1,l+1)
			ey(ip) = ey(ip) + sx(0)*sz(1)*eyg(j,1,l+1)
			ey(ip) = ey(ip) + sx(1)*sz(1)*eyg(j+1,1,l+1)
			ey(ip) = ey(ip) + sx(2)*sz(1)*eyg(j+2,1,l+1)
			ey(ip) = ey(ip) + sx(-1)*sz(2)*eyg(j-1,1,l+2)
			ey(ip) = ey(ip) + sx(0)*sz(2)*eyg(j,1,l+2)
			ey(ip) = ey(ip) + sx(1)*sz(2)*eyg(j+1,1,l+2)
			ey(ip) = ey(ip) + sx(2)*sz(2)*eyg(j+2,1,l+2)

			! Compute Ez on particle
			ez(ip) = ez(ip) + sx(-1)*sz0(-1)*ezg(j-1,1,l0-1)
			ez(ip) = ez(ip) + sx(0)*sz0(-1)*ezg(j,1,l0-1)
			ez(ip) = ez(ip) + sx(1)*sz0(-1)*ezg(j+1,1,l0-1)
			ez(ip) = ez(ip) + sx(2)*sz0(-1)*ezg(j+2,1,l0-1)
			ez(ip) = ez(ip) + sx(-1)*sz0(0)*ezg(j-1,1,l0)
			ez(ip) = ez(ip) + sx(0)*sz0(0)*ezg(j,1,l0)
			ez(ip) = ez(ip) + sx(1)*sz0(0)*ezg(j+1,1,l0)
			ez(ip) = ez(ip) + sx(2)*sz0(0)*ezg(j+2,1,l0)
			ez(ip) = ez(ip) + sx(-1)*sz0(1)*ezg(j-1,1,l0+1)
			ez(ip) = ez(ip) + sx(0)*sz0(1)*ezg(j,1,l0+1)
			ez(ip) = ez(ip) + sx(1)*sz0(1)*ezg(j+1,1,l0+1)
			ez(ip) = ez(ip) + sx(2)*sz0(1)*ezg(j+2,1,l0+1)
			ez(ip) = ez(ip) + sx(-1)*sz0(2)*ezg(j-1,1,l0+2)
			ez(ip) = ez(ip) + sx(0)*sz0(2)*ezg(j,1,l0+2)
			ez(ip) = ez(ip) + sx(1)*sz0(2)*ezg(j+1,1,l0+2)
			ez(ip) = ez(ip) + sx(2)*sz0(2)*ezg(j+2,1,l0+2)

		end do
#if defined _OPENMP && _OPENMP>=201307
			!$OMP END SIMD 
#endif
	ENDIF
	return
 end subroutine pxr_gete2dxz_energy_conserving_3_3


!_________________________________________________________________________________________
subroutine pxr_getb2dxz_energy_conserving_3_3(np,xp,zp,bx,by,bz,xmin,zmin,dx,dz,nx,nz,&
                    nxguard,nzguard,bxg,byg,bzg,l_lower_order_in_v)
!                                       
! Field gathering cartesian in 2D for the magnetic field
! This function is vectorized                                    
!_________________________________________________________________________________________                                       
      use constants
      implicit none

      ! __ Parameter declaration ___________________________________________
      integer(idp)                       :: np,nx,nz,nox,noz,nxguard,nzguard
      real(num), dimension(np)           :: xp,zp,bx,by,bz
      logical(idp)                       :: l_lower_order_in_v
      real(num), dimension(-nxguard:nx+nxguard,1,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(num)                          :: xmin,zmin,dx,dz
      integer(idp)                       :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                                            ixmin0, ixmax0, izmin0, izmax0, jj, ll, j0, l0
      real(num)                          :: dxi, dzi, x, y, z, xint, zint, &
                                            xintsq,oxint,zintsq,ozint,oxintsq,ozintsq
      real(num), DIMENSION(-1:2)         :: sx, sx0
      real(num), DIMENSION(-1:2)         :: sz, sz0
      real(num), parameter               :: onesixth=1./6.,twothird=2./3.

      ! ___________________________     
      ! Compute parameters

      dxi = 1./dx
      dzi = 1./dz

      sx=0
      sz=0.
      sx0=0.
      sz0=0.

      ! ___________________________   
      IF (l_lower_order_in_v) THEN

#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
      DO ip=1,np
    
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
    
    ! Compute index of particle
    j=floor(x)
    j0=floor(x)
    
    l=floor(z)
    l0=floor(z)
    
    xint=x-j
    zint=z-l
    
    ! Compute shape factors
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1) = onesixth*oxintsq*oxint
    sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
    sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
    sx( 2) = onesixth*xintsq*xint
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1) = onesixth*ozintsq*ozint
    sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
    sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
    sz( 2) = onesixth*zintsq*zint
    
    xint=x-0.5_num-j0
    zint=z-0.5_num-l0
    
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2

    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = 0.75_num-zintsq
    sz0( 1) = 0.5_num*(0.5_num+zint)**2
    
    ! Compute Bx on particle
    bx(ip) = bx(ip) + sx(-1)*sz0(-1)*bxg(j-1,1,l0-1)
    bx(ip) = bx(ip) + sx(0)*sz0(-1)*bxg(j,1,l0-1)
    bx(ip) = bx(ip) + sx(1)*sz0(-1)*bxg(j+1,1,l0-1)
    bx(ip) = bx(ip) + sx(2)*sz0(-1)*bxg(j+2,1,l0-1)
    bx(ip) = bx(ip) + sx(-1)*sz0(0)*bxg(j-1,1,l0)
    bx(ip) = bx(ip) + sx(0)*sz0(0)*bxg(j,1,l0)
    bx(ip) = bx(ip) + sx(1)*sz0(0)*bxg(j+1,1,l0)
    bx(ip) = bx(ip) + sx(2)*sz0(0)*bxg(j+2,1,l0)
    bx(ip) = bx(ip) + sx(-1)*sz0(1)*bxg(j-1,1,l0+1)
    bx(ip) = bx(ip) + sx(0)*sz0(1)*bxg(j,1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sz0(1)*bxg(j+1,1,l0+1)
    bx(ip) = bx(ip) + sx(2)*sz0(1)*bxg(j+2,1,l0+1)
    
    ! Compute By on particle
    by(ip) = by(ip) + sx0(-1)*sz0(-1)*byg(j0-1,1,l0-1)
    by(ip) = by(ip) + sx0(0)*sz0(-1)*byg(j0,1,l0-1)
    by(ip) = by(ip) + sx0(1)*sz0(-1)*byg(j0+1,1,l0-1)
    by(ip) = by(ip) + sx0(-1)*sz0(0)*byg(j0-1,1,l0)
    by(ip) = by(ip) + sx0(0)*sz0(0)*byg(j0,1,l0)
    by(ip) = by(ip) + sx0(1)*sz0(0)*byg(j0+1,1,l0)
    by(ip) = by(ip) + sx0(-1)*sz0(1)*byg(j0-1,1,l0+1)
    by(ip) = by(ip) + sx0(0)*sz0(1)*byg(j0,1,l0+1)
    by(ip) = by(ip) + sx0(1)*sz0(1)*byg(j0+1,1,l0+1)
    
    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(-1)*sz(-1)*bzg(j0-1,1,l-1)
    bz(ip) = bz(ip) + sx0(0)*sz(-1)*bzg(j0,1,l-1)
    bz(ip) = bz(ip) + sx0(1)*sz(-1)*bzg(j0+1,1,l-1)
    bz(ip) = bz(ip) + sx0(-1)*sz(0)*bzg(j0-1,1,l)
    bz(ip) = bz(ip) + sx0(0)*sz(0)*bzg(j0,1,l)
    bz(ip) = bz(ip) + sx0(1)*sz(0)*bzg(j0+1,1,l)
    bz(ip) = bz(ip) + sx0(-1)*sz(1)*bzg(j0-1,1,l+1)
    bz(ip) = bz(ip) + sx0(0)*sz(1)*bzg(j0,1,l+1)
    bz(ip) = bz(ip) + sx0(1)*sz(1)*bzg(j0+1,1,l+1)
    bz(ip) = bz(ip) + sx0(-1)*sz(2)*bzg(j0-1,1,l+2)
    bz(ip) = bz(ip) + sx0(0)*sz(2)*bzg(j0,1,l+2)
    bz(ip) = bz(ip) + sx0(1)*sz(2)*bzg(j0+1,1,l+2)      
      
      end do      
      
      ! ___________________________         
      ! l_lower_order_in_v is false
      ELSE


#if defined _OPENMP && _OPENMP>=201307
			!$OMP SIMD 
#elif defined __IBMBGQ__
			!IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
			!DIR$ SIMD 
#endif
#if defined __INTEL_COMPILER 
!DIR$ IVDEP
!DIR$ DISTRIBUTE POINT
#endif
      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        z = (zp(ip)-zmin)*dzi
    
        ! Compute index of particle
        j=floor(x)
        j0=floor(x-0.5_num)
        
        l=floor(z)
        l0=floor(z-0.5_num)
    
        xint=x-j
        zint=z-l

        ! Compute shape factors
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
        sx( 2) = onesixth*xintsq*xint
    
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
        sz( 2) = onesixth*zintsq*zint

        xint=x-0.5-j0
        zint=z-0.5-l0

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

    ! Compute Bx on particle
    bx(ip) = bx(ip) + sx(-1)*sz0(-1)*bxg(j-1,1,l0-1)
    bx(ip) = bx(ip) + sx(0)*sz0(-1)*bxg(j,1,l0-1)
    bx(ip) = bx(ip) + sx(1)*sz0(-1)*bxg(j+1,1,l0-1)
    bx(ip) = bx(ip) + sx(2)*sz0(-1)*bxg(j+2,1,l0-1)
    bx(ip) = bx(ip) + sx(-1)*sz0(0)*bxg(j-1,1,l0)
    bx(ip) = bx(ip) + sx(0)*sz0(0)*bxg(j,1,l0)
    bx(ip) = bx(ip) + sx(1)*sz0(0)*bxg(j+1,1,l0)
    bx(ip) = bx(ip) + sx(2)*sz0(0)*bxg(j+2,1,l0)
    bx(ip) = bx(ip) + sx(-1)*sz0(1)*bxg(j-1,1,l0+1)
    bx(ip) = bx(ip) + sx(0)*sz0(1)*bxg(j,1,l0+1)
    bx(ip) = bx(ip) + sx(1)*sz0(1)*bxg(j+1,1,l0+1)
    bx(ip) = bx(ip) + sx(2)*sz0(1)*bxg(j+2,1,l0+1)
    bx(ip) = bx(ip) + sx(-1)*sz0(2)*bxg(j-1,1,l0+2)
    bx(ip) = bx(ip) + sx(0)*sz0(2)*bxg(j,1,l0+2)
    bx(ip) = bx(ip) + sx(1)*sz0(2)*bxg(j+1,1,l0+2)
    bx(ip) = bx(ip) + sx(2)*sz0(2)*bxg(j+2,1,l0+2)

    ! Compute By on particle
    by(ip) = by(ip) + sx0(-1)*sz0(-1)*byg(j0-1,1,l0-1)
    by(ip) = by(ip) + sx0(0)*sz0(-1)*byg(j0,1,l0-1)
    by(ip) = by(ip) + sx0(1)*sz0(-1)*byg(j0+1,1,l0-1)
    by(ip) = by(ip) + sx0(2)*sz0(-1)*byg(j0+2,1,l0-1)
    by(ip) = by(ip) + sx0(-1)*sz0(0)*byg(j0-1,1,l0)
    by(ip) = by(ip) + sx0(0)*sz0(0)*byg(j0,1,l0)
    by(ip) = by(ip) + sx0(1)*sz0(0)*byg(j0+1,1,l0)
    by(ip) = by(ip) + sx0(2)*sz0(0)*byg(j0+2,1,l0)
    by(ip) = by(ip) + sx0(-1)*sz0(1)*byg(j0-1,1,l0+1)
    by(ip) = by(ip) + sx0(0)*sz0(1)*byg(j0,1,l0+1)
    by(ip) = by(ip) + sx0(1)*sz0(1)*byg(j0+1,1,l0+1)
    by(ip) = by(ip) + sx0(2)*sz0(1)*byg(j0+2,1,l0+1)
    by(ip) = by(ip) + sx0(-1)*sz0(2)*byg(j0-1,1,l0+2)
    by(ip) = by(ip) + sx0(0)*sz0(2)*byg(j0,1,l0+2)
    by(ip) = by(ip) + sx0(1)*sz0(2)*byg(j0+1,1,l0+2)
    by(ip) = by(ip) + sx0(2)*sz0(2)*byg(j0+2,1,l0+2)

    ! Compute Bz on particle
    bz(ip) = bz(ip) + sx0(-1)*sz(-1)*bzg(j0-1,1,l-1)
    bz(ip) = bz(ip) + sx0(0)*sz(-1)*bzg(j0,1,l-1)
    bz(ip) = bz(ip) + sx0(1)*sz(-1)*bzg(j0+1,1,l-1)
    bz(ip) = bz(ip) + sx0(2)*sz(-1)*bzg(j0+2,1,l-1)
    bz(ip) = bz(ip) + sx0(-1)*sz(0)*bzg(j0-1,1,l)
    bz(ip) = bz(ip) + sx0(0)*sz(0)*bzg(j0,1,l)
    bz(ip) = bz(ip) + sx0(1)*sz(0)*bzg(j0+1,1,l)
    bz(ip) = bz(ip) + sx0(2)*sz(0)*bzg(j0+2,1,l)
    bz(ip) = bz(ip) + sx0(-1)*sz(1)*bzg(j0-1,1,l+1)
    bz(ip) = bz(ip) + sx0(0)*sz(1)*bzg(j0,1,l+1)
    bz(ip) = bz(ip) + sx0(1)*sz(1)*bzg(j0+1,1,l+1)
    bz(ip) = bz(ip) + sx0(2)*sz(1)*bzg(j0+2,1,l+1)
    bz(ip) = bz(ip) + sx0(-1)*sz(2)*bzg(j0-1,1,l+2)
    bz(ip) = bz(ip) + sx0(0)*sz(2)*bzg(j0,1,l+2)
    bz(ip) = bz(ip) + sx0(1)*sz(2)*bzg(j0+1,1,l+2)
    bz(ip) = bz(ip) + sx0(2)*sz(2)*bzg(j0+2,1,l+2)
                 
     end do
#if defined _OPENMP && _OPENMP>=201307
			!$OMP END SIMD 
#endif

   ENDIF
   return
 end subroutine pxr_getb2dxz_energy_conserving_3_3
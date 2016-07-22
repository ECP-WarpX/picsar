#! /usr/bin/python
from numpy import *

# Parameters
nox=3 # order of gathering
noz=3
l4symtry=False
l_particles_weight=True
final_loop_lin=2
filename="depose_jxjyjz_esirkepov_"+str(nox)+"_"+str(noz)+".F90"
subroutine_deposej="depose_jxjyjz_esirkepov2d_lin_"+str(nox)+"_"+str(noz)
indent_cont_1="                                      "
indent_cont_2="              "
indent_1="  "
indent_2=indent_1+indent_1
indent_3=indent_2+indent_1
indent_4=indent_3+indent_1
indent_5=indent_4+indent_1
indent_6=indent_5+indent_1


# Routine Variables explicitly written
ndtodx = 0
ndtodz = 0
xl = -int(nox/2)-1-ndtodx
xu = int((nox+1)/2)+1+ndtodx
zl = -int(noz/2)-1-ndtodz
zu = int((noz+1)/2)+1+ndtodz
ixmin = -1-int(nox/2)
ixmax = 1+int((nox+1)/2)
izmin = -1-int(noz/2)
izmax = 1+int((noz+1)/2)

# Function to handle "+ int" when int=0 and int =/0
def plusi(i):
    if (i == 0):
        return ""
    else:
        if (i>0):
            return "+"+str(abs(i))
        else:
            return "-"+str(abs(i))

# Check
print 'ixmin',ixmin,'ixmax',ixmax
print 'izmin',izmin,'izmax',izmax

# Open file
fh=open(filename,"w")



# ------- Current deposition

fh.write("SUBROUTINE "+subroutine_deposej+"(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &\n");
fh.write(indent_cont_1+"dt,dx,dz,nx,nz,nxguard,nzguard, &\n");
fh.write(indent_cont_1+"nox,noz,l_particles_weight)\n");

fh.write(indent_1+"USE omp_lib\n")
fh.write(indent_1+"USE constants\n")
fh.write(indent_1+"implicit none\n")
fh.write(indent_1+"\n")
fh.write(indent_1+"! Parameter declaration\n");
fh.write(indent_1+"integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard\n")
fh.write(indent_1+"real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz\n")
fh.write(indent_1+"real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w\n")
fh.write(indent_1+"real(num)                             :: q,dt,dx,dz,xmin,zmin\n")
fh.write(indent_1+"logical(idp)                          :: l_particles_weight\n")

fh.write(indent_1+"real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint\n")
fh.write(indent_1+"real(num),dimension(:,:), allocatable :: sdx,sdz\n")
fh.write(indent_1+"real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz\n")
fh.write(indent_1+"real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz\n")
fh.write(indent_1+"real(num)                             :: invvol,invdtdx,invdtdz\n")
fh.write(indent_1+"real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq\n")
fh.write(indent_1+"real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0\n")
fh.write(indent_1+"real(num), parameter                  :: onesixth=1./6.,twothird=2./3.\n")
fh.write(indent_1+"real(num), parameter                  :: onethird=1./3.\n")
fh.write(indent_1+"real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz\n")
fh.write(indent_1+"integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc\n")
fh.write(indent_1+"integer(idp)                          :: ixmin, ixmax, izmin, izmax, icell, ndtodx, ndtodz\n")
fh.write(indent_1+"integer(idp)                          :: xl,xu,zl,zu\n")

fh.write(indent_1+"\n");
fh.write(indent_1+"! Parameter initialization\n");
fh.write(indent_1+"dxi = 1.0_num/dx\n");
fh.write(indent_1+"dzi = 1.0_num/dz\n");
fh.write(indent_1+"dtsdx0 = dt*dxi\n");
fh.write(indent_1+"dtsdz0 = dt*dzi\n");
fh.write(indent_1+"invvol = 1.0_num/(dx*dz)\n");
fh.write(indent_1+"invdtdx = 1.0_num/(dt*dz)\n");
fh.write(indent_1+"invdtdz = 1.0_num/(dt*dx)\n");
fh.write(indent_1+"dtsdz0 = dt*dzi\n");

fh.write(indent_1+"allocate(sdx(%d:%d,%d:%d),sdz(%d:%d,%d:%d))\n"%(xl,xu,zl,zu,xl,xu,zl,zu))
fh.write(indent_1+"ALLOCATE(sx("+str(xl)+":"+str(xu)+"), sx0("+str(xl)+":"+str(xu)+"), dsx("+str(xl)+":"+str(xu)+"))\n");
fh.write(indent_1+"ALLOCATE(sz("+str(zl)+":"+str(zu)+"), sz0("+str(zl)+":"+str(zu)+"), dsz("+str(zl)+":"+str(zu)+"))\n");

fh.write(indent_1+"sx0=0.0_num;sz0=0.0_num\n");
fh.write(indent_1+"sdx=0.0_num;sdz=0.0_num\n");

fh.write(indent_1+"\n");
fh.write(indent_1+"DO ip=1,np\n");

fh.write(indent_1+"\n");
fh.write(indent_2+"! --- computes current position in grid units\n");
fh.write(indent_2+"x = (xp(ip)-xmin)*dxi\n");
fh.write(indent_2+"z = (zp(ip)-zmin)*dzi\n");

fh.write(indent_1+"\n");
fh.write(indent_2+"! --- computes velocity\n");
fh.write(indent_2+"vx = uxp(ip)*gaminv(ip)\n");
fh.write(indent_2+"vy = uyp(ip)*gaminv(ip)\n");
fh.write(indent_2+"vz = uzp(ip)*gaminv(ip)\n");

fh.write(indent_1+"\n");
fh.write(indent_2+"! --- computes old position in grid units\n");
fh.write(indent_2+"xold=x-dtsdx0*vx\n");
fh.write(indent_2+"zold=z-dtsdz0*vz\n");
if(l4symtry):
    fh.write(indent_2+"! --- applies 4-fold symmetry\n");
    fh.write(indent_2+"x=abs(x)\n");
    fh.write(indent_2+"xold=abs(xold)\n");
    fh.write(indent_2+"vx = (x-xold)/dtsdx0\n");

fh.write(indent_1+"\n");    
fh.write(indent_2+"! --- computes particles weights\n");
if (l_particles_weight):
    fh.write(indent_2+"wq=q*w(ip)\n");
else:
    fh.write(indent_2+"wq=q*w(1)\n");
fh.write(indent_2+"wqx = wq*invdtdx\n");
fh.write(indent_2+"wqz = wq*invdtdz\n");

fh.write(indent_1+"\n");
fh.write(indent_2+"! --- finds node of cell containing particles for current positions\n");
if(nox==2*(nox/2)):
    fh.write(indent_2+"iixp0=nint(x)\n");
else:
    fh.write(indent_2+"iixp0=floor(x)\n");
if(noz==2*(noz/2)):
    fh.write(indent_2+"ikxp0=nint(z)\n");
else:
    fh.write(indent_2+"ikxp0=floor(z)\n");

fh.write(indent_1+"\n");    
fh.write(indent_2+"! --- computes distance between particle and node for current positions\n");
fh.write(indent_2+"xint=x-iixp0\n");
fh.write(indent_2+"zint=z-ikxp0\n")

fh.write(indent_1+"\n");   
fh.write(indent_2+"! --- computes coefficients for node centered quantities\n")
if(nox==0):
    fh.write(indent_2+"sx0( 0) = 1.0_num\n");
if(nox==1):
    fh.write(indent_2+"sx0( 0) = 1.0_num-xint\n");
    fh.write(indent_2+"sx0( 1) = xint\n");
if(nox==2):
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"sx0(-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_2+"sx0( 0) = 0.75_num-xintsq\n");
    fh.write(indent_2+"sx0( 1) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_2+"oxint = 1.0_num-xint\n");
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"oxintsq = oxint*oxint\n");
    fh.write(indent_2+"sx0(-1) = onesixth*oxintsq*oxint\n");
    fh.write(indent_2+"sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
    fh.write(indent_2+"sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
    fh.write(indent_2+"sx0( 2) = onesixth*xintsq*xint\n");

fh.write(indent_1+"\n");   
if(noz==0):
    fh.write(indent_2+"sz0( 0) = 1.0_num\n");
if(noz==1):
    fh.write(indent_2+"sz0( 0) = 1.0_num-zint\n");
    fh.write(indent_2+"sz0( 1) = zint\n");
if(noz==2):
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"sz0(-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_2+"sz0( 0) = 0.75_num-zintsq\n");
    fh.write(indent_2+"sz0( 1) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_2+"ozint = 1.0_num-zint\n");
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"ozintsq = ozint*ozint\n");
    fh.write(indent_2+"sz0(-1) = onesixth*ozintsq*ozint\n");
    fh.write(indent_2+"sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
    fh.write(indent_2+"sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
    fh.write(indent_2+"sz0( 2) = onesixth*zintsq*zint\n");

fh.write(indent_1+"\n");   
fh.write(indent_2+"! --- finds node of cell containing particles for old positions\n");
if(nox==2*(nox/2)):
    fh.write(indent_2+"iixp=nint(xold)\n");
else:
    fh.write(indent_2+"iixp=floor(xold)\n");
if(noz==2*(noz/2)):
    fh.write(indent_2+"ikxp=nint(zold)\n");
else:
    fh.write(indent_2+"ikxp=floor(zold)\n");

fh.write(indent_1+"\n");   
fh.write(indent_2+"! --- computes distance between particle and node for old positions\n");
fh.write(indent_2+"xint = xold-iixp\n");
fh.write(indent_2+"zint = zold-ikxp\n");

fh.write(indent_1+"\n");   
fh.write(indent_2+"! --- computes node separation between old and current positions\n");
fh.write(indent_2+"dix = iixp-iixp0\n");
fh.write(indent_2+"diz = ikxp-ikxp0\n");

fh.write(indent_1+"\n");   
fh.write(indent_2+"! --- zero out coefficients (needed because of different dix and diz for each particle)\n");
for i in range(xl,xu+1):
  fh.write(indent_2+"sx(%d)=0.0_num\n"%(i));
for i in range(zl,zu+1):  
  fh.write(indent_2+"sz(%d)=0.0_num\n"%(i));

fh.write(indent_1+"\n");       
fh.write(indent_2+"! --- computes coefficients for quantities centered between nodes\n");
if(nox==0):
    fh.write(indent_2+"sx( 0+dix) = 1.0_num\n");
if(nox==1):
    fh.write(indent_2+"sx( 0+dix) = 1.0_num-xint\n");
    fh.write(indent_2+"sx( 1+dix) = xint\n");
if(nox==2):
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"sx(-1+dix) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_2+"sx( 0+dix) = 0.75_num-xintsq\n");
    fh.write(indent_2+"sx( 1+dix) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_2+"oxint = 1.0_num-xint\n");
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"oxintsq = oxint*oxint\n");
    fh.write(indent_2+"sx(-1+dix) = onesixth*oxintsq*oxint\n");
    fh.write(indent_2+"sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
    fh.write(indent_2+"sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
    fh.write(indent_2+"sx( 2+dix) = onesixth*xintsq*xint\n");

fh.write(indent_1+"\n");   
if(noz==0):
    fh.write(indent_2+"sz( 0+diz) = 1.0_num\n");
if(noz==1):
    fh.write(indent_2+"sz( 0+diz) = 1.0_num-zint\n");
    fh.write(indent_2+"sz( 1+diz) = zint\n");
if(noz==2):
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"sz(-1+diz) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_2+"sz( 0+diz) = 0.75_num-zintsq\n");
    fh.write(indent_2+"sz( 1+diz) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_2+"ozint = 1.0_num-zint\n");
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"ozintsq = ozint*ozint\n");
    fh.write(indent_2+"sz(-1+diz) = onesixth*ozintsq*ozint\n");
    fh.write(indent_2+"sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
    fh.write(indent_2+"sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
    fh.write(indent_2+"sz( 2+diz) = onesixth*zintsq*zint\n");

fh.write(indent_1+"\n");   
fh.write(indent_2+"! --- computes coefficients difference\n");
for i in range(xl,xu+1):
  fh.write(indent_2+"dsx(%d) = sx(%d) - sx0(%d)\n"%(i,i,i));
  fh.write(indent_2+"dsz(%d) = sz(%d) - sz0(%d)\n"%(i,i,i));
fh.write(indent_2+"\n")

# The final loop ic completely linearized
# This is for the vectorized version
if final_loop_lin==1:
  fh.write(indent_2+"! --- add current contributions\n")

  # For x
  fh.write(indent_1+"\n");  
  l = 0
  for k in range(izmin,izmax+1):  
    for i in range(ixmin,8+ixmin-1+1):
      l += 1
      if(i<ixmax):
        #print i,k,l      
        fh.write(indent_2+"sdx(n,%d)  = wqx*dsx(%d)*( sz0(%d) + 0.5*dsz(%d) )\n"%(l,i,k,k));
        if (i>ixmin):
          fh.write(indent_2+"sdx(n,%d) = sdx(n,%d)+sdx(n,%d)\n"%(l,l,l-1));       
  # For y
  fh.write(indent_1+"\n");    
  l = 0
  for k in range(izmin,izmax+1):  
    for i in range(ixmin,8+ixmin-1+1):
       l += 1   
       if (i<=ixmax):
         fh.write(indent_2+"sdy(n,%d) = wq*vy*invvol* &\n"%(l)) 
         fh.write(indent_2+"( (sz0(%d)+0.5*dsz(%d))*sx0(%d) + (0.5*sz0(%d)+onethird*dsz(%d))*dsx(%d))\n"%(k,k,i,k,k,i))   
     
  # For z
  fh.write(indent_1+"\n");  
  l = 0
  for k in range(izmin,izmax+1):  
    for i in range(ixmin,8+ixmin-1+1):
      l+=1
      if((k<izmax)and(i<=ixmax)):
        fh.write(indent_2+"sdz(n,%d)  = wqz*dsz(%d)*(sx0(%d)+0.5*dsx(%d))\n"%(l,k,i,i)); 
        if(k>izmin):
          fh.write(indent_2+"sdz(n,%d) = sdz(n,%d)+sdz(n,%d)\n"%(l,l,l-1));             
                            
# The final loop is not linearized
elif (final_loop_lin==0):

  fh.write(indent_2+"! --- computes min/max positions of current contributions\n")
  fh.write(indent_2+"ixmin = min(0,dix)-"+str(int(nox/2))+"\n")
  fh.write(indent_2+"ixmax = max(0,dix)+"+str(int((nox+1)/2))+"\n")
  fh.write(indent_2+"izmin = min(0,diz)-"+str(int(noz/2))+"\n")
  fh.write(indent_2+"izmax = max(0,diz)+"+str(int((noz+1)/2))+"\n")
  
  fh.write(indent_2+"\n")
  fh.write(indent_2+"! --- add current contributions\n")
  fh.write(indent_2+"DO k=izmin, izmax\n")
  fh.write(indent_3+"DO i=ixmin, ixmax\n")
  fh.write(indent_4+"ic = iixp0+i\n")
  fh.write(indent_4+"kc = ikxp0+k\n")

  fh.write(indent_4+"\n")
  fh.write(indent_4+"! --- Jx\n")  
  fh.write(indent_4+"IF(i<ixmax) THEN\n")
  fh.write(indent_5+"sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov\n")
  fh.write(indent_5+"if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x \n")
  fh.write(indent_5+"jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current\n")
  fh.write(indent_4+"END IF\n")

  fh.write(indent_4+"\n")
  fh.write(indent_4+"! -- Jy (2D Esirkepov scheme)\n")  
  fh.write(indent_4+"jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &\n") 
  fh.write(indent_4+"( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )\n")                

  fh.write(indent_4+"\n") 
  fh.write(indent_4+"! --- Jz\n") 
  fh.write(indent_4+"IF(k<izmax) THEN\n")
  fh.write(indent_5+"sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov\n")
  fh.write(indent_5+"if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z\n")
  fh.write(indent_5+"jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current\n")
  fh.write(indent_4+"END IF\n")
 
  fh.write(indent_3+"END DO\n")
  fh.write(indent_2+"END DO\n")

# The final loop is not linearized
elif (final_loop_lin==2):
    
  # For x
  fh.write(indent_1+"\n");    
  l = 0
  for k in range(izmin,izmax+1):  
    for i in range(ixmin,ixmax+1):
      if (i<ixmax):
        fh.write(indent_4+"sdx(n,%d,%d)  = wqx*dsx(%d)*( sz0(%d) + 0.5*dsz(%d) )  \n"%(i,k,i,k,k))
        if (i>ixmin):
          fh.write(indent_4+"sdx(n,%d,%d)=sdx(n,%d,%d)+sdx(n,%d-1,%d)    \n"%(i,k,i,k,i,k))  
              
  # For y
  fh.write(indent_1+"\n");    
  l = 0
  for k in range(izmin,izmax+1):  
    for i in range(ixmin,ixmax+1):
        fh.write(indent_4+"sdy(n,%d,%d) = wq*vy*invvol* &\n"%(i,k)) 
        fh.write(indent_4+"( (sz0(%d)+0.5*dsz(%d))*sx0(%d) + (0.5*sz0(%d)+onethird*dsz(%d))*dsx(%d) )\n"%(k,k,i,k,k,i))     
    
  # For z
  fh.write(indent_1+"\n");    
  l = 0
  for k in range(izmin,izmax+1):  
    for i in range(ixmin,ixmax+1):
      if (k<izmax):
        fh.write(indent_5+"sdz(n,%d,%d)=wqz*dsz(%d)*(sx0(%d)+0.5*dsx(%d))    \n"%(i,k,k,i,i))
        if (k>izmin):
          fh.write(indent_5+"sdz(n,%d,%d)=sdz(n,%d,%d)+sdz(n,%d,%d-1)    \n"%(i,k,i,k,i,k))
                                                                
# End do particles  
fh.write("\n")
fh.write(indent_1+"END DO\n")


fh.write(indent_1+"DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)\n")

fh.write(indent_1+"RETURN\n");
fh.write("END SUBROUTINE "+subroutine_deposej+"\n");

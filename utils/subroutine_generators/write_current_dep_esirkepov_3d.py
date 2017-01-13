#! /usr/bin/python
from numpy import *

# Parameters
nox=1 # order of gathering
noy=1
noz=1
l4symtry=False
l_particles_weight=True
final_loop_lin=0
filename="depose_jxjyjz_esirkepov_"+str(nox)+"_"+str(noy)+"_"+str(noz)+".F90"
subroutine_deposej="depose_jxjyjz_esirkepov_lin_"+str(nox)+"_"+str(noy)+"_"+str(noz)
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
ndtody = 0
ndtodz = 0
xl = -int(nox/2)-1-ndtodx
xu = int((nox+1)/2)+1+ndtodx
yl = -int(noy/2)-1-ndtody
yu = int((noy+1)/2)+1+ndtody
zl = -int(noz/2)-1-ndtodz
zu = int((noz+1)/2)+1+ndtodz
ixmin = -1-int(nox/2)
ixmax = 1+int((nox+1)/2)
iymin = -1-int(noy/2)
iymax = 1+int((noy+1)/2)
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

# Open file
fh=open(filename,"w")



# ------- Current deposition

fh.write("SUBROUTINE "+subroutine_deposej+"(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &\n");
fh.write(indent_cont_1+"dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &\n");
fh.write(indent_cont_1+"nox,noy,noz,l_particles_weight,l4symtry)\n");
fh.write(indent_1+"USE omp_lib\n");
fh.write(indent_1+"USE constants\n");
fh.write(indent_1+"IMPLICIT NONE\n");
fh.write(indent_1+"INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard\n");
fh.write(indent_1+"REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz\n");
#fh.write(indent_1+"REAL(num), DIMENSION(:,:,:), ALLOCATABLE:: jx1, jy1, jz1\n");
fh.write(indent_1+"REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv\n");
fh.write(indent_1+"REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin\n");
fh.write(indent_1+"REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint\n");
fh.write(indent_1+"REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz\n");
fh.write(indent_1+"REAL(num) :: clghtisq,usq,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz, &\n");
fh.write(indent_cont_1+"s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,         &\n");
fh.write(indent_cont_1+"oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &\n");
fh.write(indent_cont_1+"dtsdx0,dtsdy0,dtsdz0\n");
fh.write(indent_1+"REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz\n");
fh.write(indent_1+"INTEGER(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &\n");
fh.write(indent_cont_1+"ixmin, ixmax, iymin, iymax, izmin, izmax\n");
fh.write(indent_1+"LOGICAL(idp) :: l_particles_weight,l4symtry\n");
fh.write("\n");

fh.write("! PARAMETER INIT\n");
fh.write("dxi = 1.0_num/dx\n");
fh.write("dyi = 1.0_num/dy\n");
fh.write("dzi = 1.0_num/dz\n");
fh.write("dtsdx0 = dt*dxi\n");
fh.write("dtsdy0 = dt*dyi\n");
fh.write("dtsdz0 = dt*dzi\n");
fh.write("invvol = 1.0_num/(dx*dy*dz)\n");
fh.write("invdtdx = 1.0_num/(dt*dy*dz)\n");
fh.write("invdtdy = 1.0_num/(dt*dx*dz)\n");
fh.write("invdtdz = 1.0_num/(dt*dx*dy)\n");
fh.write("ALLOCATE(sdx("+str(xl)+":"+str(xu)+","+str(yl)+":"+str(yu)+","+str(zl)+":"+str(zu)+"),sdy("\
        +str(xl)+":"+str(xu)+","+str(yl)+":"+str(yu)+","+str(zl)+":"+str(zu)+"),sdz("+str(xl)+":"+str(xu)+\
        ","+str(yl)+":"+str(yu)+","+str(zl)+":"+str(zu)+"))\n");
fh.write("ALLOCATE(sx("+str(xl)+":"+str(xu)+"), sx0("+str(xl)+":"+str(xu)+"), dsx("+str(xl)+":"+str(xu)+"))\n");
fh.write("ALLOCATE(sy("+str(yl)+":"+str(yu)+"), sy0("+str(yl)+":"+str(yu)+"), dsy("+str(yl)+":"+str(yu)+"))\n");
fh.write("ALLOCATE(sz("+str(zl)+":"+str(zu)+"), sz0("+str(zl)+":"+str(zu)+"), dsz("+str(zl)+":"+str(zu)+"))\n");
#fh.write("ALLOCATE(jx1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), &\n");
#fh.write("         jy1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), &\n");
#fh.write("         jz1(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))\n");
fh.write("clghtisq = 1.0_num/clight**2\n");
fh.write("sx0=0.0_num;sy0=0.0_num;sz0=0.0_num\n");
fh.write("sdx=0.0_num;sdy=0.0_num;sdz=0.0_num\n");
#fh.write("jx1=0.0_num;jy1=0.0_num;jz1=0.0_num\n");
fh.write("dtsdz0 = dt*dzi\n");
#fh.write("!!$OMP PARALLEL PRIVATE(ip,x,y,z,usq,vx,vy,vz,gaminv,xold,yold,zold, &\n");
#fh.write("!!$OMP wq,wqx,wqy,wqz,iixp0,ijxp0,ikxp0, xint,yint,zint, oxint,xintsq, oxintsq,dix,diy,diz, &\n");
#fh.write("!!$OMP dsx, dsy, dsz, oyint,yintsq, oyintsq, ozint,zintsq, ozintsq,ixmin, ixmax, iymin, iymax, izmin, izmax,  &\n");
#fh.write("!!$OMP k,j,i,kc,jc,ic, iixp, ijxp, ikxp,sx,sy,sz) FIRSTPRIVATE(sx0,sy0,sz0,sdx,sdy,sdz,jx1,jy1,jz1)&\n");
#fh.write("!!$OMP DO\n");
fh.write("DO ip=1,np\n");
fh.write(indent_1+"! --- computes current position in grid units\n");
fh.write(indent_1+"x = (xp(ip)-xmin)*dxi\n");
fh.write(indent_1+"y = (yp(ip)-ymin)*dyi\n");
fh.write(indent_1+"z = (zp(ip)-zmin)*dzi\n");
fh.write(indent_1+"! --- computes velocity\n");
#fh.write(indent_1+"usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq\n");
#fh.write(indent_1+"gaminv = 1.0_num/sqrt(1.0_num + usq)\n");
fh.write(indent_1+"vx = uxp(ip)*gaminv(ip)\n");
fh.write(indent_1+"vy = uyp(ip)*gaminv(ip)\n");
fh.write(indent_1+"vz = uzp(ip)*gaminv(ip)\n");
fh.write(indent_1+"! --- computes old position in grid units\n");
fh.write(indent_1+"xold=x-dtsdx0*vx\n");
fh.write(indent_1+"yold=y-dtsdy0*vy\n");
fh.write(indent_1+"zold=z-dtsdz0*vz\n");
if(l4symtry):
    fh.write(indent_1+"! --- applies 4-fold symmetry\n");
    fh.write(indent_1+"x=abs(x)\n");
    fh.write(indent_1+"y=abs(y)\n");
    fh.write(indent_1+"xold=abs(xold)\n");
    fh.write(indent_1+"yold=abs(yold)\n");
    fh.write(indent_1+"vx = (x-xold)/dtsdx0\n");
    fh.write(indent_1+"vy = (y-yold)/dtsdy0\n");
fh.write(indent_1+"! --- computes particles weights\n");
if (l_particles_weight):
    fh.write(indent_1+"wq=q*w(ip)\n");
else:
    fh.write(indent_1+"wq=q*w(1)\n");
fh.write(indent_1+"wqx = wq*invdtdx\n");
fh.write(indent_1+"wqy = wq*invdtdy\n");
fh.write(indent_1+"wqz = wq*invdtdz\n");

fh.write(indent_1+"! --- finds node of cell containing particles for current positions\n");
if(nox==2*(nox/2)):
    fh.write(indent_1+"iixp0=nint(x)\n");
else:
    fh.write(indent_1+"iixp0=floor(x)\n");
if(noy==2*(noy/2)):
    fh.write(indent_1+"ijxp0=nint(y)\n");
else:
    fh.write(indent_1+"ijxp0=floor(y)\n");
if(noz==2*(noz/2)):
    fh.write(indent_1+"ikxp0=nint(z)\n");
else:
    fh.write(indent_1+"ikxp0=floor(z)\n");
    
fh.write(indent_1+"! --- computes distance between particle and node for current positions\n");
fh.write(indent_1+"xint=x-iixp0\n");
fh.write(indent_1+"yint=y-ijxp0\n");
fh.write(indent_1+"zint=z-ikxp0\n")

fh.write(indent_1+"! --- computes coefficients for node centered quantities\n")
if(nox==0):
    fh.write(indent_1+"sx0( 0) = 1.0_num\n");
if(nox==1):
    fh.write(indent_1+"sx0( 0) = 1.0_num-xint\n");
    fh.write(indent_1+"sx0( 1) = xint\n");
if(nox==2):
    fh.write(indent_1+"xintsq = xint*xint\n");
    fh.write(indent_1+"sx0(-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_1+"sx0( 0) = 0.75_num-xintsq\n");
    fh.write(indent_1+"sx0( 1) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_1+"oxint = 1.0_num-xint\n");
    fh.write(indent_1+"xintsq = xint*xint\n");
    fh.write(indent_1+"oxintsq = oxint*oxint\n");
    fh.write(indent_1+"sx0(-1) = onesixth*oxintsq*oxint\n");
    fh.write(indent_1+"sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)\n");
    fh.write(indent_1+"sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)\n");
    fh.write(indent_1+"sx0( 2) = onesixth*xintsq*xint\n");
    
if(noy==0):
    fh.write(indent_1+"sy0( 0) = 1.0_num\n");
if (noy==1):
    fh.write(indent_1+"sy0( 0) = 1.0_num-yint\n");
    fh.write(indent_1+"sy0( 1) = yint\n");
if(noy==2):
    fh.write(indent_1+"yintsq = yint*yint\n");
    fh.write(indent_1+"sy0(-1) = 0.5_num*(0.5_num-yint)**2\n");
    fh.write(indent_1+"sy0( 0) = 0.75_num-yintsq\n");
    fh.write(indent_1+"sy0( 1) = 0.5_num*(0.5_num+yint)**2\n");
if(noy==3):
    fh.write(indent_1+"oyint = 1.0_num-yint\n");
    fh.write(indent_1+"yintsq = yint*yint\n");
    fh.write(indent_1+"oyintsq = oyint*oyint\n");
    fh.write(indent_1+"sy0(-1) = onesixth*oyintsq*oyint\n");
    fh.write(indent_1+"sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)\n");
    fh.write(indent_1+"sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)\n");
    fh.write(indent_1+"sy0( 2) = onesixth*yintsq*yint\n");

if(noz==0):
    fh.write(indent_1+"sz0( 0) = 1.0_num\n");
if(noz==1):
    fh.write(indent_1+"sz0( 0) = 1.0_num-zint\n");
    fh.write(indent_1+"sz0( 1) = zint\n");
if(noz==2):
    fh.write(indent_1+"zintsq = zint*zint\n");
    fh.write(indent_1+"sz0(-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_1+"sz0( 0) = 0.75_num-zintsq\n");
    fh.write(indent_1+"sz0( 1) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_1+"ozint = 1.0_num-zint\n");
    fh.write(indent_1+"zintsq = zint*zint\n");
    fh.write(indent_1+"ozintsq = ozint*ozint\n");
    fh.write(indent_1+"sz0(-1) = onesixth*ozintsq*ozint\n");
    fh.write(indent_1+"sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)\n");
    fh.write(indent_1+"sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)\n");
    fh.write(indent_1+"sz0( 2) = onesixth*zintsq*zint\n");

fh.write(indent_1+"! --- finds node of cell containing particles for old positions\n");
if(nox==2*(nox/2)):
    fh.write(indent_1+"iixp=nint(xold)\n");
else:
    fh.write(indent_1+"iixp=floor(xold)\n");
if(noy==2*(noy/2)):
    fh.write(indent_1+"ijxp=nint(yold)\n");
else:
    fh.write(indent_1+"ijxp=floor(yold)\n");
if(noz==2*(noz/2)):
    fh.write(indent_1+"ikxp=nint(zold)\n");
else:
    fh.write(indent_1+"ikxp=floor(zold)\n");

fh.write(indent_1+"! --- computes distance between particle and node for old positions\n");
fh.write(indent_1+"xint = xold-iixp\n");
fh.write(indent_1+"yint = yold-ijxp\n");
fh.write(indent_1+"zint = zold-ikxp\n");

fh.write(indent_1+"! --- computes node separation between old and current positions\n");
fh.write(indent_1+"dix = iixp-iixp0\n");
fh.write(indent_1+"diy = ijxp-ijxp0\n");
fh.write(indent_1+"diz = ikxp-ikxp0\n");

fh.write(indent_1+"! --- zero out coefficients (needed because of different dix and diz for each particle)\n");
fh.write(indent_1+"sx=0.0_num;sy=0.0_num;sz=0.0_num\n");

fh.write(indent_1+"! --- computes coefficients for quantities centered between nodes\n");

if(nox==0):
    fh.write(indent_1+"sx( 0+dix) = 1.0_num\n");
if(nox==1):
    fh.write(indent_1+"sx( 0+dix) = 1.0_num-xint\n");
    fh.write(indent_1+"sx( 1+dix) = xint\n");
if(nox==2):
    fh.write(indent_1+"xintsq = xint*xint\n");
    fh.write(indent_1+"sx(-1+dix) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_1+"sx( 0+dix) = 0.75_num-xintsq\n");
    fh.write(indent_1+"sx( 1+dix) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_1+"oxint = 1.0_num-xint\n");
    fh.write(indent_1+"xintsq = xint*xint\n");
    fh.write(indent_1+"oxintsq = oxint*oxint\n");
    fh.write(indent_1+"sx(-1+dix) = onesixth*oxintsq*oxint\n");
    fh.write(indent_1+"sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)\n");
    fh.write(indent_1+"sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)\n");
    fh.write(indent_1+"sx( 2+dix) = onesixth*xintsq*xint\n");
    
if(noy==0):
    fh.write(indent_1+"sy( 0+diy) = 1.0_num\n");
if (noy==1):
    fh.write(indent_1+"sy( 0+diy) = 1.0_num-yint\n");
    fh.write(indent_1+"sy( 1+diy) = yint\n");
if(noy==2):
    fh.write(indent_1+"yintsq = yint*yint\n");
    fh.write(indent_1+"sy(-1+diy) = 0.5_num*(0.5_num-yint)**2\n");
    fh.write(indent_1+"sy( 0+diy) = 0.75_num-yintsq\n");
    fh.write(indent_1+"sy( 1+diy) = 0.5_num*(0.5_num+yint)**2\n");
if(noy==3):
    fh.write(indent_1+"oyint = 1.0_num-yint\n");
    fh.write(indent_1+"yintsq = yint*yint\n");
    fh.write(indent_1+"oyintsq = oyint*oyint\n");
    fh.write(indent_1+"sy(-1+diy) = onesixth*oyintsq*oyint\n");
    fh.write(indent_1+"sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)\n");
    fh.write(indent_1+"sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)\n");
    fh.write(indent_1+"sy( 2+diy) = onesixth*yintsq*yint\n");

if(noz==0):
    fh.write(indent_1+"sz( 0+diz) = 1.0_num\n");
if(noz==1):
    fh.write(indent_1+"sz( 0+diz) = 1.0_num-zint\n");
    fh.write(indent_1+"sz( 1+diz) = zint\n");
if(noz==2):
    fh.write(indent_1+"zintsq = zint*zint\n");
    fh.write(indent_1+"sz(-1+diz) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_1+"sz( 0+diz) = 0.75_num-zintsq\n");
    fh.write(indent_1+"sz( 1+diz) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_1+"ozint = 1.0_num-zint\n");
    fh.write(indent_1+"zintsq = zint*zint\n");
    fh.write(indent_1+"ozintsq = ozint*ozint\n");
    fh.write(indent_1+"sz(-1+diz) = onesixth*ozintsq*ozint\n");
    fh.write(indent_1+"sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)\n");
    fh.write(indent_1+"sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)\n");
    fh.write(indent_1+"sz( 2+diz) = onesixth*zintsq*zint\n");

fh.write(indent_1+"! --- computes coefficients difference\n");
fh.write(indent_1+"dsx = sx - sx0\n");
fh.write(indent_1+"dsy = sy - sy0\n");
fh.write(indent_1+"dsz = sz - sz0\n");
fh.write(indent_1+"\n")
# The final loop ic completely linearized
if final_loop_lin==0:
  fh.write(indent_1+"! --- add current contributions\n")
#   for k in range(izmin, izmax+1):
#        for j in range(iymin, iymax+1):
#            for i  in range(ixmin, ixmax+1):
#                 ic="iixp0"+plusi(i)
#                 jc="ijxp0"+plusi(j)
#                 kc="ikxp0"+plusi(k)
#                 if(i<ixmax):
#                     fh.write(indent_1+"sdx("+str(i)+","+str(j)+","+str(k)+")  = wqx*dsx("+str(i)+")*((sy0("+str(j)+")+0.5_num*dsy("+str(j)+"))*sz0("+str(k)+") + &\n");
#                     fh.write(indent_1+"(0.5_num*sy0("+str(j)+")+1.0_num/3.0_num*dsy("+str(j)+"))*dsz("+str(k)+"))\n");
#                     if (i>ixmin):
#                         fh.write(indent_1+"sdx("+str(i)+","+str(j)+","+str(k)+")=sdx("+str(i)+","+str(j)+","+str(k)+")+sdx("+str(i-1)+","+str(j)+","+str(k)+")\n");
#                     fh.write(indent_1+"jx("+ic+","+jc+","+kc+") = jx("+ic+","+jc+","+kc+") + sdx("+str(i)+","+str(j)+","+str(k)+")\n");
#                 if(j<iymax):
#                     fh.write(indent_1+"sdy("+str(i)+","+str(j)+","+str(k)+")  = wqy*dsy("+str(j)+")*((sz0("+str(k)+")+0.5_num*dsz("+str(k)+"))*sx0("+str(i)+") + &\n");
#                     fh.write(indent_1+"(0.5_num*sz0("+str(k)+")+1.0_num/3.0_num*dsz("+str(k)+"))*dsx("+str(i)+"))\n");
#                     if(j>iymin):
#                         fh.write(indent_1+"sdy("+str(i)+","+str(j)+","+str(k)+")=sdy("+str(i)+","+str(j)+","+str(k)+")+sdy("+str(i)+","+str(j-1)+","+str(k)+")\n");
#                     fh.write(indent_1+"jy("+ic+","+jc+","+kc+") = jy("+ic+","+jc+","+kc+") + sdy("+str(i)+","+str(j)+","+str(k)+")\n");
#                 if(k<izmax):
#                     fh.write(indent_1+"sdz("+str(i)+","+str(j)+","+str(k)+")  = wqz*dsz("+str(k)+")*((sx0("+str(i)+")+0.5_num*dsx("+str(i)+"))*sy0("+str(j)+") + &\n");
#                     fh.write(indent_1+"(0.5_num*sx0("+str(i)+")+1.0_num/3.0_num*dsx("+str(i)+"))*dsy("+str(j)+"))\n");
#                     if(k>izmin):
#                         fh.write(indent_1+"sdz("+str(i)+","+str(j)+","+str(k)+")=sdz("+str(i)+","+str(j)+","+str(k)+")+sdz("+str(i)+","+str(j)+","+str(k-1)+")\n");
#                     fh.write(indent_1+"jz("+ic+","+jc+","+kc+") = jz("+ic+","+jc+","+kc+") + sdz("+str(i)+","+str(j)+","+str(k)+")\n");

  fh.write(indent_1+"\n")
  for k in range(izmin, izmax+1):
       for j in range(iymin, iymax+1):
           for i  in range(ixmin, ixmax+1):
                if(i<ixmax):
                    fh.write(indent_1+"sdx("+str(i)+","+str(j)+","+str(k)+")  = wqx*dsx("+str(i)+")*((sy0("+str(j)+\
                    ")+0.5_num*dsy("+str(j)+"))*sz0("+str(k)+") + &\n");
                    fh.write(indent_1+"(0.5_num*sy0("+str(j)+")+1.0_num/3.0_num*dsy("+str(j)+"))*dsz("+str(k)+"))\n");
                    if (i>ixmin):
                        fh.write(indent_1+"sdx("+str(i)+","+str(j)+","+str(k)+")=sdx("+str(i)+","+str(j)+\
                        ","+str(k)+")+sdx("+str(i-1)+","+str(j)+","+str(k)+")\n");

  fh.write(indent_1+"\n")                   
  for k in range(izmin, izmax+1):
       for j in range(iymin, iymax+1):
           for i  in range(ixmin, ixmax+1):                
                if(j<iymax):
                    fh.write(indent_1+"sdy("+str(i)+","+str(j)+","+str(k)+")  = wqy*dsy("+str(j)+")*((sz0("+str(k)+")+0.5_num*dsz("+str(k)+"))*sx0("+str(i)+") + &\n");
                    fh.write(indent_1+"(0.5_num*sz0("+str(k)+")+1.0_num/3.0_num*dsz("+str(k)+"))*dsx("+str(i)+"))\n");
                    if(j>iymin):
                        fh.write(indent_1+"sdy("+str(i)+","+str(j)+","+str(k)+")=sdy("+str(i)+","+\
                        str(j)+","+str(k)+")+sdy("+str(i)+","+str(j-1)+","+str(k)+")\n");

  fh.write(indent_1+"\n")
  for k in range(izmin, izmax+1):
       for j in range(iymin, iymax+1):
           for i  in range(ixmin, ixmax+1):                  
                if(k<izmax):
                    fh.write(indent_1+"sdz("+str(i)+","+str(j)+","+str(k)+")  = wqz*dsz("+str(k)+")*((sx0("+str(i)+")+0.5_num*dsx("+str(i)+"))*sy0("+str(j)+") + &\n");
                    fh.write(indent_1+"(0.5_num*sx0("+str(i)+")+1.0_num/3.0_num*dsx("+str(i)+"))*dsy("+str(j)+"))\n");
                    if(k>izmin):
                        fh.write(indent_1+"sdz("+str(i)+","+str(j)+","+str(k)+")=sdz("+str(i)+","\
                        +str(j)+","+str(k)+")+sdz("+str(i)+","+str(j)+","+str(k-1)+")\n");                        

  fh.write(indent_1+"\n")                                            
  for k in range(izmin, izmax+1):
       for j in range(iymin, iymax+1):
           for i  in range(ixmin, ixmax+1):
                ic="iixp0"+plusi(i)
                jc="ijxp0"+plusi(j)
                kc="ikxp0"+plusi(k)
                if(i<ixmax):                    
                    fh.write(indent_1+"jx("+ic+","+jc+","+kc+") = jx("+ic+","+jc+","+kc+") + sdx("+str(i)+\
                    ","+str(j)+","+str(k)+")\n");

  fh.write(indent_1+"\n")  
  for k in range(izmin, izmax+1):
       for j in range(iymin, iymax+1):
           for i  in range(ixmin, ixmax+1):
                ic="iixp0"+plusi(i)
                jc="ijxp0"+plusi(j)
                kc="ikxp0"+plusi(k)                    
                if(j<iymax):
                    fh.write(indent_1+"jy("+ic+","+jc+","+kc+") = jy("+ic+","+jc+","+kc+") + sdy("+str(i)+\
                    ","+str(j)+","+str(k)+")\n");

  fh.write(indent_1+"\n")
  for k in range(izmin, izmax+1):
       for j in range(iymin, iymax+1):
           for i  in range(ixmin, ixmax+1):
                ic="iixp0"+plusi(i)
                jc="ijxp0"+plusi(j)
                kc="ikxp0"+plusi(k)                    
                if(k<izmax):
                    fh.write(indent_1+"jz("+ic+","+jc+","+kc+") = jz("+ic+","+jc+","+kc+") + sdz("+str(i)+","+str(j)+","+str(k)+")\n");


# The final loop is not linearized
elif (final_loop_lin==1):

  fh.write(indent_1+"! --- computes min/max positions of current contributions\n")
  fh.write(indent_1+"ixmin = min(0,dix)-"+str(int(nox/2))+"\n")
  fh.write(indent_1+"ixmax = max(0,dix)+"+str(int((nox+1)/2))+"\n")
  fh.write(indent_1+"iymin = min(0,diy)-"+str(int(noy/2))+"\n")
  fh.write(indent_1+"iymax = max(0,diy)+"+str(int((noy+1)/2))+"\n")
  fh.write(indent_1+"izmin = min(0,diz)-"+str(int(noz/2))+"\n")
  fh.write(indent_1+"izmax = max(0,diz)+"+str(int((noz+1)/2))+"\n")

  fh.write(indent_1+"! --- add current contributions\n")
  fh.write(indent_1+"DO k=izmin, izmax\n")
  fh.write(indent_2+"DO j=iymin, iymax\n")
  fh.write(indent_3+"DO i=ixmin, ixmax\n")
  fh.write(indent_4+"ic = iixp0+i\n")
  fh.write(indent_4+"jc = ijxp0+j\n")
  fh.write(indent_4+"kc = ikxp0+k\n")
  
  fh.write(indent_4+"IF(i<ixmax) THEN\n")
  fh.write(indent_5+"sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &\n")
  fh.write(indent_5+"(0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k)) \n")
  fh.write(indent_5+"IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k) \n")
  fh.write(indent_5+"jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)\n")
  fh.write(indent_4+"END IF\n")
  
  fh.write(indent_4+"IF(j<iymax) THEN\n")
  fh.write(indent_5+"sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &\n")
  fh.write(indent_5+"(0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))\n")
  fh.write(indent_5+"IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)\n")
  fh.write(indent_5+"jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)\n")
  fh.write(indent_4+"END IF\n")
 
  fh.write(indent_4+"IF(k<izmax) THEN\n")
  fh.write(indent_5+"sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &\n")
  fh.write(indent_5+"(0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))\n")
  fh.write(indent_5+"IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)\n")
  fh.write(indent_5+"jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)\n")
  fh.write(indent_4+"END IF\n")
 
  fh.write(indent_3+"END DO\n")
  fh.write(indent_2+"END DO\n")
  fh.write(indent_1+"END DO\n")

elif (final_loop_lin==2):

  fh.write(indent_1+"! --- computes min/max positions of current contributions\n")
  fh.write(indent_1+"ixmin = min(0,dix)-"+str(int(nox/2))+"\n")
  fh.write(indent_1+"ixmax = max(0,dix)+"+str(int((nox+1)/2))+"\n")
  fh.write(indent_1+"iymin = min(0,diy)-"+str(int(noy/2))+"\n")
  fh.write(indent_1+"iymax = max(0,diy)+"+str(int((noy+1)/2))+"\n")
  fh.write(indent_1+"izmin = min(0,diz)-"+str(int(noz/2))+"\n")
  fh.write(indent_1+"izmax = max(0,diz)+"+str(int((noz+1)/2))+"\n")

  fh.write(indent_1+"! --- add current contributions\n")
  fh.write(indent_1+"DO k=izmin, izmax\n")
  fh.write(indent_1+"!$OMP SIMD aligned(sy0,dyz:64)\n")
  fh.write(indent_2+"DO j=iymin, iymax\n")
  fh.write(indent_3+"jc = ijxp0+j\n")
  fh.write(indent_3+"kc = ikxp0+k\n")

  for i  in range(ixmin, ixmax+1):
    ic="iixp0"+plusi(i)
    if(i<ixmax):
      fh.write(indent_3+"sdx("+str(i)+",j,k)  = wqx*dsx("+str(i)+")*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &\n");
      fh.write(indent_3+"(0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))\n");
      if (i>ixmin):
        fh.write(indent_3+"sdx("+str(i)+",j,k)=sdx("+str(i)+",j,k)+sdx("+str(i-1)+",j,k)\n");
      fh.write(indent_3+"jx("+ic+",jc,kc) = jx("+ic+",jc,kc) + sdx("+str(i)+",j,k)\n");  
                  
  fh.write(indent_2+"END DO\n")
  fh.write(indent_1+"END DO\n")

  fh.write(indent_1+"DO k=izmin, izmax\n")
  fh.write(indent_1+"!$OMP SIMD aligned(sx0,dsx:64)\n")
  fh.write(indent_2+"DO i=ixmin, ixmax\n")
  fh.write(indent_3+"ic = iixp0+i\n")
  fh.write(indent_3+"kc = ikxp0+k\n")

  for j in range(iymin, iymax+1):
    jc="ijxp0"+plusi(j)
    if(j<iymax):
      fh.write(indent_3+"sdy(i,"+str(j)+",k)  = wqy*dsy("+str(j)+")*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &\n");
      fh.write(indent_3+"(0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))\n");
      if(j>iymin):
        fh.write(indent_3+"sdy(i,"+str(j)+",k)=sdy(i,"+\
        str(j)+",k)+sdy(i,"+str(j-1)+",k)\n");
      fh.write(indent_3+"jy(ic,"+jc+",kc) = jy(ic,"+jc+",kc) + sdy(i,"+str(j)+",k)\n");

  fh.write(indent_2+"END DO\n")
  fh.write(indent_1+"END DO\n")

  fh.write(indent_2+"DO j=iymin, iymax\n")
  fh.write(indent_1+"!$OMP SIMD aligned(sx0,dsx:64)\n")
  fh.write(indent_3+"DO i=ixmin, ixmax\n")
  fh.write(indent_3+"ic = iixp0+i\n")
  fh.write(indent_3+"jc = ijxp0+j\n")
 
  for k in range(izmin, izmax+1):
    kc="ikxp0"+plusi(k) 
    if(k<izmax):
      fh.write(indent_3+"sdz(i,j,"+str(k)+")  = wqz*dsz("+str(k)+")*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &\n");
      fh.write(indent_3+"(0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))\n");
      if(k>izmin):
        fh.write(indent_3+"sdz(i,j,"+str(k)+")=sdz(i,"\
        +str(j)+","+str(k)+")+sdz(i,j,"+str(k-1)+")\n");   
      fh.write(indent_3+"jz(ic,jc,"+kc+") = jz(ic,jc,"+kc+") + sdz(i,j,"+str(k)+")\n");
                     
  fh.write(indent_2+"END DO\n")
  fh.write(indent_1+"END DO\n")

                                                            
# End do particles  
fh.write("\n")
fh.write("END DO\n")
#fh.write("!!$OMP END DO\n")

#fh.write("!$OMP CRITICAL\n")
#fh.write("jx=jx+jx1\n")
#fh.write("jy=jy+jy1\n")
#fh.write("jz=jz+jz1\n")
#fh.write("!$OMP END CRITICAL\n")
#fh.write("!!$OMP END PARALLEL\n")
fh.write("DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)\n")

fh.write("RETURN\n");
fh.write("END SUBROUTINE "+subroutine_deposej+"\n");
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
subroutine_deposej="depose_jxjyjz_esirkepov_gpu_"+str(nox)+"_"+str(noy)+"_"+str(noz)
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
fh.write("SUBROUTINE "+subroutine_deposej+"(jx, jx_nguard, jx_nvalid, jy, jy_nguard, &\n");
fh.write(indent_cont_1+"jy_nvalid, jz, jz_nguard, jz_nvalid, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &\n");
fh.write(indent_cont_1+"xmin, ymin, zmin, dt, dx, dy, dz)\n");
fh.write(indent_1+"USE constants\n");
fh.write(indent_1+"USE constants, ONLY: clight\n");
fh.write(indent_1+"USE picsar_precision, ONLY: idp, num\n");
fh.write(indent_1+"IMPLICIT NONE\n");
fh.write(indent_1+"INTEGER(idp) :: np\n");
fh.write(indent_1+"INTEGER(idp), intent(in) :: jx_nguard(3), jx_nvalid(3), jy_nguard(3), jy_nvalid(3), &\n");
fh.write(indent_1+"jz_nguard(3), jz_nvalid(3)\n");
fh.write(indent_1+"REAL(num), intent(IN OUT):: jx(-jx_nguard(1):jx_nvalid(1)+jx_nguard(1)-1,&\n");
fh.write(indent_1+"                            -jx_nguard(2):jx_nvalid(2)+jx_nguard(2)-1,&\n");
fh.write(indent_1+"                            -jx_nguard(3):jx_nvalid(3)+jx_nguard(3)-1 )\n");
fh.write(indent_1+"REAL(num), intent(IN OUT):: jy(-jy_nguard(1):jy_nvalid(1)+jy_nguard(1)-1,&\n");
fh.write(indent_1+"                            -jy_nguard(2):jy_nvalid(2)+jy_nguard(2)-1,&\n");
fh.write(indent_1+"                            -jy_nguard(3):jy_nvalid(3)+jy_nguard(3)-1 )\n");
fh.write(indent_1+"REAL(num), intent(IN OUT):: jz(-jz_nguard(1):jz_nvalid(1)+jz_nguard(1)-1,&\n");
fh.write(indent_1+"                            -jz_nguard(2):jz_nvalid(2)+jz_nguard(2)-1,&\n");
fh.write(indent_1+"                            -jz_nguard(3):jz_nvalid(3)+jz_nguard(3)-1 )\n");
fh.write(indent_1+"REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv\n");
fh.write(indent_1+"REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin\n");
fh.write(indent_1+"REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint\n");
fh.write(indent_1+"REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz\n");
fh.write(indent_1+"REAL(num) :: clghtisq,usq,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz, &\n");
fh.write(indent_1+"             s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,&\n");
fh.write(indent_1+"             oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &\n");
fh.write(indent_1+"             dtsdx0,dtsdy0,dtsdz0\n");
fh.write(indent_1+"REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num\n");
fh.write(indent_1+"REAL(num), PARAMETER :: onethird=1.0_num/3.0_num\n");
fh.write(indent_1+"REAL(num), PARAMETER :: twothird=2.0_num/3.0_num\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz\n");
fh.write(indent_1+"INTEGER(idp) :: iixp0,ijxp0,ikxp0\n");
fh.write(indent_1+"INTEGER(idp) :: iixp,ijxp,ikxp,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &\n");
fh.write(indent_1+"                ixmin, ixmax, iymin, iymax, izmin, izmax\n");
fh.write(indent_1+"INTEGER(isp) :: n,nn\n");
fh.write(indent_1+"LOGICAL(idp) :: l_particles_weight,l4symtry\n");
fh.write("\n");

fh.write(indent_1+"! PARAMETER INIT\n");
fh.write(indent_1+"dxi = 1.0_num/dx\n");
fh.write(indent_1+"dyi = 1.0_num/dy\n");
fh.write(indent_1+"dzi = 1.0_num/dz\n");
fh.write(indent_1+"dtsdx0 = dt*dxi\n");
fh.write(indent_1+"dtsdy0 = dt*dyi\n");
fh.write(indent_1+"dtsdz0 = dt*dzi\n");
fh.write(indent_1+"invvol = 1.0_num/(dx*dy*dz)\n");
fh.write(indent_1+"invdtdx = 1.0_num/(dt*dy*dz)\n");
fh.write(indent_1+"invdtdy = 1.0_num/(dt*dx*dz)\n");
fh.write(indent_1+"invdtdz = 1.0_num/(dt*dx*dy)\n");
fh.write(indent_1+"ALLOCATE(sdx("+str(xl)+":"+str(xu)+","+str(yl)+":"+str(yu)+","+str(zl)+":"+str(zu)+"),sdy("\
        +str(xl)+":"+str(xu)+","+str(yl)+":"+str(yu)+","+str(zl)+":"+str(zu)+"),sdz("+str(xl)+":"+str(xu)+\
        ","+str(yl)+":"+str(yu)+","+str(zl)+":"+str(zu)+"))\n");
fh.write(indent_1+"ALLOCATE(sx("+str(xl)+":"+str(xu)+"), sx0("+str(xl)+":"+str(xu)+"), dsx("+str(xl)+":"+str(xu)+"))\n");
fh.write(indent_1+"ALLOCATE(sy("+str(yl)+":"+str(yu)+"), sy0("+str(yl)+":"+str(yu)+"), dsy("+str(yl)+":"+str(yu)+"))\n");
fh.write(indent_1+"ALLOCATE(sz("+str(zl)+":"+str(zu)+"), sz0("+str(zl)+":"+str(zu)+"), dsz("+str(zl)+":"+str(zu)+"))\n");
fh.write(indent_1+"clghtisq = 1.0_num/clight**2\n");
fh.write(indent_1+"sx0=0.0_num;sy0=0.0_num;sz0=0.0_num\n");
fh.write(indent_1+"sdx=0.0_num;sdy=0.0_num;sdz=0.0_num\n");
fh.write(indent_1+"dtsdz0 = dt*dzi\n");
fh.write(indent_1+"!$acc parallel deviceptr(jx, jy, jz, xp, yp, zp, uxp, uyp, uzp, w, gaminv)\n");
fh.write(indent_1+"!$acc loop gang vector private(sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1))\n");
fh.write(indent_1+"DO ip=1,np\n");
fh.write(indent_2+"! --- computes current position in grid units\n");
fh.write(indent_2+"x = (xp(ip)-xmin)*dxi\n");
fh.write(indent_2+"y = (yp(ip)-ymin)*dyi\n");
fh.write(indent_2+"z = (zp(ip)-zmin)*dzi\n");
fh.write(indent_2+"! --- computes velocity\n");
#fh.write(indent_1+"gaminv = 1.0_num/sqrt(1.0_num + usq)\n");
fh.write(indent_2+"vx = uxp(ip)*gaminv(ip)\n");
fh.write(indent_2+"vy = uyp(ip)*gaminv(ip)\n");
fh.write(indent_2+"vz = uzp(ip)*gaminv(ip)\n");
fh.write(indent_2+"! --- computes old position in grid units\n");
fh.write(indent_2+"xold=x-dtsdx0*vx\n");
fh.write(indent_2+"yold=y-dtsdy0*vy\n");
fh.write(indent_2+"zold=z-dtsdz0*vz\n");
if(l4symtry):
    fh.write(indent_2+"! --- applies 4-fold symmetry\n");
    fh.write(indent_2+"x=abs(x)\n");
    fh.write(indent_2+"y=abs(y)\n");
    fh.write(indent_2+"xold=abs(xold)\n");
    fh.write(indent_2+"yold=abs(yold)\n");
    fh.write(indent_2+"vx = (x-xold)/dtsdx0\n");
    fh.write(indent_2+"vy = (y-yold)/dtsdy0\n");
fh.write(indent_2+"! --- computes particles weights\n");
if (l_particles_weight):
    fh.write(indent_2+"wq=q*w(ip)\n");
else:
    fh.write(indent_2+"wq=q*w(1)\n");
fh.write(indent_2+"wqx = wq*invdtdx\n");
fh.write(indent_2+"wqy = wq*invdtdy\n");
fh.write(indent_2+"wqz = wq*invdtdz\n");

fh.write(indent_2+"! --- finds node of cell containing particles for current positions\n");
if(nox==2*(nox/2)):
    fh.write(indent_2+"iixp0=nint(x)\n");
else:
    fh.write(indent_2+"iixp0=floor(x)\n");
if(noy==2*(noy/2)):
    fh.write(indent_2+"ijxp0=nint(y)\n");
else:
    fh.write(indent_2+"ijxp0=floor(y)\n");
if(noz==2*(noz/2)):
    fh.write(indent_2+"ikxp0=nint(z)\n");
else:
    fh.write(indent_2+"ikxp0=floor(z)\n");
    
fh.write(indent_2+"! --- computes distance between particle and node for current positions\n");
fh.write(indent_2+"xint=x-iixp0\n");
fh.write(indent_2+"yint=y-ijxp0\n");
fh.write(indent_2+"zint=z-ikxp0\n")

fh.write(indent_2+"! --- computes coefficients for node centered quantities\n")
fh.write(indent_2+"sx0=0.0_num;sy0=0.0_num;sz0=0.0_num\n");
fh.write(indent_2+"sdx=0.0_num;sdy=0.0_num;sdz=0.0_num\n");
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
    fh.write(indent_2+"sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)\n");
    fh.write(indent_2+"sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)\n");
    fh.write(indent_2+"sx0( 2) = onesixth*xintsq*xint\n");
    
if(noy==0):
    fh.write(indent_2+"sy0( 0) = 1.0_num\n");
if (noy==1):
    fh.write(indent_2+"sy0( 0) = 1.0_num-yint\n");
    fh.write(indent_2+"sy0( 1) = yint\n");
if(noy==2):
    fh.write(indent_2+"yintsq = yint*yint\n");
    fh.write(indent_2+"sy0(-1) = 0.5_num*(0.5_num-yint)**2\n");
    fh.write(indent_2+"sy0( 0) = 0.75_num-yintsq\n");
    fh.write(indent_2+"sy0( 1) = 0.5_num*(0.5_num+yint)**2\n");
if(noy==3):
    fh.write(indent_2+"oyint = 1.0_num-yint\n");
    fh.write(indent_2+"yintsq = yint*yint\n");
    fh.write(indent_2+"oyintsq = oyint*oyint\n");
    fh.write(indent_2+"sy0(-1) = onesixth*oyintsq*oyint\n");
    fh.write(indent_2+"sy0( 0) = twothird-yintsq*(1.0_num-yint*0.5_num)\n");
    fh.write(indent_2+"sy0( 1) = twothird-oyintsq*(1.0_num-oyint*0.5_num)\n");
    fh.write(indent_2+"sy0( 2) = onesixth*yintsq*yint\n");

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
    fh.write(indent_2+"sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)\n");
    fh.write(indent_2+"sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)\n");
    fh.write(indent_2+"sz0( 2) = onesixth*zintsq*zint\n");

fh.write(indent_2+"! --- finds node of cell containing particles for old positions\n");
if(nox==2*(nox/2)):
    fh.write(indent_2+"iixp=nint(xold)\n");
else:
    fh.write(indent_2+"iixp=floor(xold)\n");
if(noy==2*(noy/2)):
    fh.write(indent_2+"ijxp=nint(yold)\n");
else:
    fh.write(indent_2+"ijxp=floor(yold)\n");
if(noz==2*(noz/2)):
    fh.write(indent_2+"ikxp=nint(zold)\n");
else:
    fh.write(indent_2+"ikxp=floor(zold)\n");

fh.write(indent_2+"! --- computes distance between particle and node for old positions\n");
fh.write(indent_2+"xint = xold-iixp\n");
fh.write(indent_2+"yint = yold-ijxp\n");
fh.write(indent_2+"zint = zold-ikxp\n");

fh.write(indent_2+"! --- computes node separation between old and current positions\n");
fh.write(indent_2+"dix = iixp-iixp0\n");
fh.write(indent_2+"diy = ijxp-ijxp0\n");
fh.write(indent_2+"diz = ikxp-ikxp0\n");

fh.write(indent_2+"! --- zero out coefficients (needed because of different dix and diz for each particle)\n");
if(noz==1):
  for i in range(-1,3):
    fh.write(indent_2+"sx(%d)=0.0_num;sy(%d)=0.0_num;sz(%d)=0.0_num\n"%(i,i,i));
if(noz==2):
  for i in range(-2,3):
    fh.write(indent_2+"sx(%d)=0.0_num;sy(%d)=0.0_num;sz(%d)=0.0_num\n"%(i,i,i));
if(noz==3):
  for i in range(-2,3):
    fh.write(indent_2+"sx(%d)=0.0_num;sy(%d)=0.0_num;sz(%d)=0.0_num\n"%(i,i,i));
      
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
    fh.write(indent_2+"sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)\n");
    fh.write(indent_2+"sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)\n");
    fh.write(indent_2+"sx( 2+dix) = onesixth*xintsq*xint\n");
    
if(noy==0):
    fh.write(indent_2+"sy( 0+diy) = 1.0_num\n");
if (noy==1):
    fh.write(indent_2+"sy( 0+diy) = 1.0_num-yint\n");
    fh.write(indent_2+"sy( 1+diy) = yint\n");
if(noy==2):
    fh.write(indent_2+"yintsq = yint*yint\n");
    fh.write(indent_2+"sy(-1+diy) = 0.5_num*(0.5_num-yint)**2\n");
    fh.write(indent_2+"sy( 0+diy) = 0.75_num-yintsq\n");
    fh.write(indent_2+"sy( 1+diy) = 0.5_num*(0.5_num+yint)**2\n");
if(noy==3):
    fh.write(indent_2+"oyint = 1.0_num-yint\n");
    fh.write(indent_2+"yintsq = yint*yint\n");
    fh.write(indent_2+"oyintsq = oyint*oyint\n");
    fh.write(indent_2+"sy(-1+diy) = onesixth*oyintsq*oyint\n");
    fh.write(indent_2+"sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)\n");
    fh.write(indent_2+"sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)\n");
    fh.write(indent_2+"sy( 2+diy) = onesixth*yintsq*yint\n");

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
    fh.write(indent_2+"sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)\n");
    fh.write(indent_2+"sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)\n");
    fh.write(indent_2+"sz( 2+diz) = onesixth*zintsq*zint\n");

fh.write(indent_2+"! --- computes coefficients difference\n");
fh.write(indent_2+"dsx = sx - sx0\n");
fh.write(indent_2+"dsy = sy - sy0\n");
fh.write(indent_2+"dsz = sz - sz0\n");
fh.write(indent_2+"\n")

# End of loop with parameters
#fh.write("!!$OMP END DO\n")

# The final loop ic completely linearized
if final_loop_lin==0:
  fh.write(indent_2+"! --- add current contributions\n")
  for k in range(izmin, izmax+1):
        for j in range(iymin, iymax+1):
            for i  in range(ixmin, ixmax+1):
                 ic="iixp0"+plusi(i)
                 jc="ijxp0"+plusi(j)
                 kc="ikxp0"+plusi(k)
                 if(i<ixmax):
                     fh.write(indent_2+"sdx("+str(i)+","+str(j)+","+str(k)+")  = wqx*dsx("+str(i)+")*((sy0("+str(j)+")+0.5_num*dsy("+str(j)+"))*sz0("+str(k)+") + &\n");
                     fh.write(indent_2+"(0.5_num*sy0("+str(j)+")+onethird*dsy("+str(j)+"))*dsz("+str(k)+"))\n");
                     if (i>ixmin):
                         fh.write(indent_2+"sdx("+str(i)+","+str(j)+","+str(k)+")=sdx("+str(i)+","+str(j)+","+str(k)+")+sdx("+str(i-1)+","+str(j)+","+str(k)+")\n");
                     fh.write(indent_2+"!$acc atomic update\n");
                     fh.write(indent_2+"jx("+ic+","+jc+","+kc+") = jx("+ic+","+jc+","+kc+") + sdx("+str(i)+","+str(j)+","+str(k)+")\n");
                 if(j<iymax):
                     fh.write(indent_2+"sdy("+str(i)+","+str(j)+","+str(k)+")  = wqy*dsy("+str(j)+")*((sz0("+str(k)+")+0.5_num*dsz("+str(k)+"))*sx0("+str(i)+") + &\n");
                     fh.write(indent_2+"(0.5_num*sz0("+str(k)+")+onethird*dsz("+str(k)+"))*dsx("+str(i)+"))\n");
                     if(j>iymin):
                         fh.write(indent_2+"sdy("+str(i)+","+str(j)+","+str(k)+")=sdy("+str(i)+","+str(j)+","+str(k)+")+sdy("+str(i)+","+str(j-1)+","+str(k)+")\n");
                     fh.write(indent_2+"!$acc atomic update\n");
                     fh.write(indent_2+"jy("+ic+","+jc+","+kc+") = jy("+ic+","+jc+","+kc+") + sdy("+str(i)+","+str(j)+","+str(k)+")\n");
                 if(k<izmax):
                     fh.write(indent_2+"sdz("+str(i)+","+str(j)+","+str(k)+")  = wqz*dsz("+str(k)+")*((sx0("+str(i)+")+0.5_num*dsx("+str(i)+"))*sy0("+str(j)+") + &\n");
                     fh.write(indent_2+"(0.5_num*sx0("+str(i)+")+onethird*dsx("+str(i)+"))*dsy("+str(j)+"))\n");
                     if(k>izmin):
                         fh.write(indent_2+"sdz("+str(i)+","+str(j)+","+str(k)+")=sdz("+str(i)+","+str(j)+","+str(k)+")+sdz("+str(i)+","+str(j)+","+str(k-1)+")\n");
                     fh.write(indent_2+"!$acc atomic update\n");
                     fh.write(indent_2+"jz("+ic+","+jc+","+kc+") = jz("+ic+","+jc+","+kc+") + sdz("+str(i)+","+str(j)+","+str(k)+")\n");

  # Original order
  if False:
    l=0
    nshift=1
    fh.write(indent_1+"\n")
    fh.write(indent_1+"!Original order x\n")
    fh.write(indent_1+"\n")
    for k in range(izmin, izmax+1):
        for j in range(iymin, iymax+1):
            for i  in range(ixmin, ixmax+1):
               if(i<ixmax):
                 l += 1
                 fh.write(indent_1+"sdx(n,"+str(l)+")  = wqx*dsx("+str(i)+")*((sy0("+str(j)+\
                      ")+0.5_num*dsy("+str(j)+"))*sz0("+str(k)+") + &\n");
                 fh.write(indent_1+"(0.5_num*sy0("+str(j)+")+onethird*dsy("+str(j)+"))*dsz("+str(k)+"))\n");
                 if (i>ixmin):
                   fh.write(indent_1+"sdx(n,"+str(l)+")=sdx(n,"+str(l)+")+sdx(n,"+str(l-1)+")\n");

    l = 0
    nshift=(iymax+1-iymin)
    fh.write(indent_1+"!Original order y\n")     
    fh.write(indent_1+"\n")                   
    for k in range(izmin, izmax+1):
      for j in range(iymin, iymax+1):
          for i  in range(ixmin, ixmax+1):                
               if(j<iymax):
                   l += 1
                   fh.write(indent_1+"sdy(n,"+str(l)+")  = wqy*dsy("+str(j)+")*((sz0("+str(k)+")+0.5_num*dsz("+str(k)+"))*sx0("+str(i)+") + &\n");
                   fh.write(indent_1+"(0.5_num*sz0("+str(k)+")+onethird*dsz("+str(k)+"))*dsx("+str(i)+"))\n");
                   if(j>iymin):
                       fh.write(indent_1+"sdy(n,"+str(l)+")=sdy(n,"+str(l)+")+sdy(n,"+str(l-nshift)+")\n");

    l=0
    nshift = (izmax+1-izmin)*(iymax+1-iymin)
    fh.write(indent_1+"!Original order z\n")  
    fh.write(indent_1+"\n")
    for k in range(izmin, izmax+1):
          for j in range(iymin, iymax+1):
              for i  in range(ixmin, ixmax+1):                  
                   if(k<izmax):
                       l += 1
                       fh.write(indent_1+"sdz(n,"+str(l)+")  = wqz*dsz("+str(k)+")*((sx0("+str(i)+")+0.5_num*dsx("+str(i)+"))*sy0("+str(j)+") + &\n");
                       fh.write(indent_1+"(0.5_num*sx0("+str(i)+")+onethird*dsx("+str(i)+"))*dsy("+str(j)+"))\n");
                       if(k>izmin):
                           fh.write(indent_1+"sdz(n,"+str(l)+")=sdz(n,"+str(l)+")+sdz(n,"+str(l-nshift)+")\n");  

# Henri's method developped for Esirkepov: reordering of the weights
  # For the x direction
  if False:
    if(noz==1):
      l=0
      nyz=(izmax+1-izmin)*(iymax+1-iymin)
      fh.write(indent_2+"\n")
      for i  in range(ixmin, ixmax+1):
        if(i<ixmax):
          for k in range(izmin, izmax+1):
             for j in range(iymin, iymax+1):
                  l += 1
                  fh.write(indent_2+"sdx(n,"+str(l)+")  = wqx*dsx("+str(i)+")*((sy0("+str(j)+\
                        ")+0.5_num*dsy("+str(j)+"))*sz0("+str(k)+") + &\n");
                  fh.write(indent_2+"(0.5_num*sy0("+str(j)+")+onethird*dsy("+str(j)+"))*dsz("+str(k)+"))\n");
                  if (i>ixmin):
                    fh.write(indent_2+"sdx(n,"+str(l)+")=sdx(n,"+str(l)+")+sdx(n,"+str(l-nyz)+")\n");
    if(noz==2):
      l=0
      nyz=(izmax+1-izmin)*(iymax+1-iymin)  
      fh.write(indent_2+"\n")    
      for k in range(izmin, izmax+1):
          for j in range(iymin, iymax+2):     # +1 en y
              for i  in range(ixmin, ixmax+1):
                 if(i<ixmax):
                   l += 1
                   if j < iymax+1:
                     fh.write(indent_2+"sdx(n,"+str(l)+")  = wqx*dsx("+str(i)+")*((sy0("+str(j)+\
                          ")+0.5_num*dsy("+str(j)+"))*sz0("+str(k)+") + &\n");
                     fh.write(indent_2+"(0.5_num*sy0("+str(j)+")+onethird*dsy("+str(j)+"))*dsz("+str(k)+"))\n");
                     if (i>ixmin):
                       fh.write(indent_2+"sdx(n,"+str(l)+")=sdx(n,"+str(l)+")+sdx(n,"+str(l-1)+")\n");
                   else:
                      fh.write(indent_2+"sdx(n,"+str(l)+")  = 0.\n");                     

    # For the y direction
    if(noz==1):
      l=0
      nxz=(izmax+1-izmin)*4
      print ixmin,ixmax
      print iymin,iymax
      print izmin,izmax
      print nxz  
      fh.write(indent_2+"\n")        
      for j in range(iymin, iymax+1):    
        if(j<iymax):         
          for k in range(izmin, izmax+1):
               for i  in range(ixmin, ixmax+1):  
                  l += 1              
                  fh.write(indent_2+"sdy(n,"+str(l)+")  = wqy*dsy("+str(j)+")*((sz0("+str(k)+")+0.5_num*dsz("+str(k)+"))*sx0("+str(i)+") + &\n");
                  fh.write(indent_2+"(0.5_num*sz0("+str(k)+")+onethird*dsz("+str(k)+"))*dsx("+str(i)+"))\n");
                  if(j>iymin):
                     fh.write(indent_2+"sdy(n,"+str(l)+")=sdy(n,"+str(l)+")+sdy(n,"+str(l-nxz)+")\n");
    if(noz==2):
      l=0
      nshift=(ixmax+2-ixmin)
      print ixmin,ixmax
      print iymin,iymax
      print izmin,izmax
      print nshift
      fh.write(indent_2+"\n")
      for k in range(izmin, izmax+1): 
         for i  in range(ixmin, ixmax+2): # +1 en x            
            for j in range(iymin, iymax+1): 
               if(j<iymax): 
                  l += 1
                  if (i < ixmax+1):             
                    fh.write(indent_2+"sdy(n,"+str(l)+")  = wqy*dsy("+str(j)+")*((sz0("+str(k)+")+0.5_num*dsz("+str(k)+"))*sx0("+str(i)+") + &\n");
                    fh.write(indent_2+"(0.5_num*sz0("+str(k)+")+onethird*dsz("+str(k)+"))*dsx("+str(i)+"))\n");
                    if(j>iymin):
                      fh.write(indent_2+"sdy(n,"+str(l)+")=sdy(n,"+str(l)+")+sdy(n,"+str(l-1)+")\n");      
                  else:
                    fh.write(indent_2+"sdy(n,"+str(l)+")  = 0.\n");  

    # For the z direction                                      
    if(noz==1):
      l=0
      nxy=(ixmax+1-ixmin)*(iymax+1-iymin)
      print ixmin,ixmax
      print iymin,iymax
      print izmin,izmax
      print nxy
      fh.write(indent_1+"\n")
      for k in range(izmin, izmax+1):
        if(k<izmax):
           for j in range(iymin, iymax+1):
               for i  in range(ixmin, ixmax+1):                  
                  l += 1
                  fh.write(indent_1+"sdz(n,"+str(l)+")  = wqz*dsz("+str(k)+")*((sx0("+str(i)+")+0.5_num*dsx("+str(i)+"))*sy0("+str(j)+") + &\n");
                  fh.write(indent_1+"(0.5_num*sx0("+str(i)+")+onethird*dsx("+str(i)+"))*dsy("+str(j)+"))\n");
                  if(k>izmin):
                    fh.write(indent_1+"sdz(n,"+str(l)+")=sdz(n,"+str(l)+")+sdz(n,"+str(l-nxy)+")\n"); 
    if(noz==2):
      l=0
      nshift=(ixmax+2-ixmin)
      print ixmin,ixmax
      print iymin,iymax
      print izmin,izmax
      print nshift
      fh.write(indent_2+"\n")
      for j in range(iymin, iymax+1):       
         for i  in range(ixmin, ixmax+2):  # +1 en x
              for k in range(izmin, izmax+1): 
                 if(k<izmax):                
                   l += 1
                   if (i < ixmax+1):
                     fh.write(indent_2+"sdz(n,"+str(l)+")  = wqz*dsz("+str(k)+")*((sx0("+str(i)+")+0.5_num*dsx("+str(i)+"))*sy0("+str(j)+") + &\n");
                     fh.write(indent_2+"(0.5_num*sx0("+str(i)+")+onethird*dsx("+str(i)+"))*dsy("+str(j)+"))\n");
                     if(k>izmin):
                       fh.write(indent_2+"sdz(n,"+str(l)+")=sdz(n,"+str(l)+")+sdz(n,"+str(l-1)+")\n");     
                   else:
                     fh.write(indent_2+"sdz(n,"+str(l)+")  = 0.\n"); 
                      
  # Linearization of the current
  if False:
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

fh.write(indent_1+"END DO\n")
fh.write(indent_1+"!$acc end loop\n")
fh.write(indent_1+"!$acc end parallel\n")

#fh.write("!$OMP CRITICAL\n")
#fh.write("jx=jx+jx1\n")
#fh.write("jy=jy+jy1\n")
#fh.write("jz=jz+jz1\n")
#fh.write("!$OMP END CRITICAL\n")
#fh.write("!!$OMP END PARALLEL\n")
fh.write("DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)\n")

fh.write("RETURN\n");
fh.write("END SUBROUTINE "+subroutine_deposej+"\n");

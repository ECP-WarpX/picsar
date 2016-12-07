# ______________________________________________________________________________
#
# WRITE_GATHERING_ROUTINES.PY
#
# This python script enables to inline the original field gathering 
# subroutine in 3D
#
#


#! /usr/bin/python
from numpy import *

# Parameters
nox=1 # order of gathering
noy=1
noz=1

l_lower_order_in_v=True
  
filename="gathering_routines_2d_"+str(nox)+"_"+str(noz)+".F90"
subroutine_b_field="getb2dxz_energy_conserving_"+str(nox)+"_"+str(noz)
subroutine_e_field="gete2dxz_energy_conserving_"+str(nox)+"_"+str(noz)
indent_cont_1="                                      "
indent_cont_2="              "
indent_1="  "
indent_2=indent_1+indent_1
indent_3=indent_2+indent_1
indent_4=indent_3+indent_1


# Routine Variables explicitly written
ixmin = -int(nox/2)
ixmax =  int((nox+1)/2)-1
izmin = -int(noz/2)
izmax =  int((noz+1)/2)-1
if (l_lower_order_in_v):
    ixmin0 = -int((nox-1)/2)
    ixmax0 =  int((nox)/2)
    izmin0 = -int((noz-1)/2)
    izmax0 =  int((noz)/2)
else:
    ixmin0 = -int((nox)/2)
    ixmax0 =  int((nox+1)/2)
    izmin0 = -int((noz)/2)
    izmax0 =  int((noz+1)/2)


# Function to handle "+ int" when int=0 and int =/0
def plusi(i):
    if (i == 0):
        return ""
    elif (i <0):
        return str(i)    
    else:
        return "+"+str(i)

# Open file
fh=open(filename,"w")



# ------- BFIELD
# Write order n scalar gathering routine for BFIELD
# ------ Write electric field gathering/ energy conserving scheme  with OpenMP parallelization
fh.write("SUBROUTINE "+subroutine_b_field+"(np,xp,yp,zp,bx,by,bz,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard,&\n");
fh.write(indent_cont_1+"bxg,byg,bzg,l4symtry,l_2drz,l_lower_order_in_v) &\n");
fh.write(indent_1+"USE omp_lib\n");
fh.write(indent_1+"USE constants\n");
fh.write(indent_1+"IMPLICIT NONE\n");
fh.write(indent_1+"integer(idp) :: np,nx,nz,nox,noz,nxguard,nzguard\n");
fh.write(indent_1+"REAL(num), DIMENSION(np) :: xp,zp,bx,by,bz\n");
fh.write(indent_1+"REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg\n");
fh.write(indent_1+"REAL(num) :: xmin,ymin,zmin,dx,dy,dz\n");
fh.write(indent_1+"INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &\n");
fh.write(indent_cont_2+"ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0\n");
fh.write(indent_1+"REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &\n");
fh.write(indent_cont_2+"xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq\n");
fh.write(indent_1+"REAL(num), DIMENSION("+str(-int(nox/2))+":"+str(int((nox+1)/2))+") :: sx\n");
fh.write(indent_1+"REAL(num), DIMENSION("+str(-int(noy/2))+":"+str(int((noy+1)/2))+") :: sy\n");
fh.write(indent_1+"REAL(num), DIMENSION("+str(-int(noz/2))+":"+str(int((noz+1)/2))+") :: sz\n");
fh.write(indent_1+"REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0\n");
fh.write(indent_1+"REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num\n");
fh.write("\n");
fh.write(indent_1+"dxi = 1.0_num/dx\n");
fh.write(indent_1+"dzi = 1.0_num/dz\n");
fh.write(indent_1+"ALLOCATE(sx0("+str(ixmin0)+":"+str(ixmax0)+"),sz0("+str(izmin0)+":"+str(izmax0)+"))\n");
fh.write(indent_1+"sx=0.0_num\n");
fh.write(indent_1+"sz=0.0_num\n");
fh.write(indent_1+"sx0=0.0_num\n");
fh.write(indent_1+"sz0=0.0_num\n");

# _______________________________
# lower_order_in_v true

fh.write(indent_1+"IF (l_lower_order_in_v) THEN\n")
fh.write(indent_2+"DO ip=1,np,lvect\n")
fh.write(indent_3+"DO n=1,MIN(lvect,np-ip+1)\n")
fh.write(indent_3+"\n")

fh.write(indent_4+"nn=ip+n-1\n")

fh.write(indent_4+"x = (xp(nn)-xmin)*dxi\n")
fh.write(indent_4+"z = (zp(nn)-zmin)*dzi\n")
fh.write(indent_4+"\n")
fh.write(indent_4+"! Compute index of particle\n")
if(nox==2*(nox/2)):
      fh.write(indent_4+"j=nint(x)\n");
      fh.write(indent_4+"j0=floor(x-0.5_num)\n");
else:
      fh.write(indent_4+"j=floor(x)\n");
      fh.write(indent_4+"j0=floor(x)\n");
fh.write(indent_2+"\n")          
if(noz==2*(noz/2)):
      fh.write(indent_4+"l=nint(z)\n");
      fh.write(indent_4+"l0=floor(z-0.5_num)\n");
else:
      fh.write(indent_4+"l=floor(z)\n");
      fh.write(indent_4+"l0=floor(z)\n");

fh.write(indent_2+"\n")        
fh.write(indent_2+"xint=x-j\n");
fh.write(indent_2+"zint=z-l\n");
fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute shape factors\n")

if(nox==1):
    fh.write(indent_2+"sx(n, 0) = 1.0_num-xint\n");
    fh.write(indent_2+"sx(n, 1) = xint\n");
if(nox==2):
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"sx(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_2+"sx(n, 0) = 0.75_num-xintsq\n");
    fh.write(indent_2+"sx(n, 1) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_2+"oxint = 1.0_num-xint\n");
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"oxintsq = oxint*oxint\n");
    fh.write(indent_2+"sx(n,-1) = onesixth*oxintsq*oxint\n");
    fh.write(indent_2+"sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
    fh.write(indent_2+"sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
    fh.write(indent_2+"sx(n, 2) = onesixth*xintsq*xint\n");

if(noz==1):
    fh.write(indent_2+"sz(n, 0) = 1.0_num-zint\n");
    fh.write(indent_2+"sz(n, 1) = zint\n");
if(noz==2):
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"sz(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_2+"sz(n, 0) = 0.75_num-zintsq\n");
    fh.write(indent_2+"sz(n, 1) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_2+"ozint = 1.0_num-zint\n");
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"ozintsq = ozint*ozint\n");
    fh.write(indent_2+"sz(n,-1) = onesixth*ozintsq*ozint\n");
    fh.write(indent_2+"sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
    fh.write(indent_2+"sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
    fh.write(indent_2+"sz(n, 2) = onesixth*zintsq*zint\n");

fh.write(indent_2+"\n");
fh.write(indent_2+"xint=x-0.5_num-j0\n");
fh.write(indent_2+"zint=z-0.5_num-l0\n");

fh.write(indent_2+"\n");
if(nox==1):
    fh.write(indent_2+"sx0(n, 0) = 1.0_num\n");
if (nox==2):
    fh.write(indent_2+"sx0(n, 0) = 1.0_num-xint\n");
    fh.write(indent_2+"sx0(n, 1) = xint\n");
if (nox==3):
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"sx0(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_2+"sx0(n, 0) = 0.75_num-xintsq\n");
    fh.write(indent_2+"sx0(n, 1) = 0.5_num*(0.5_num+xint)**2\n");

if(noz==1):
    fh.write(indent_2+"sz0(n, 0) = 1.0_num\n");
if (noz==2):
    fh.write(indent_2+"sz0(n, 0) = 1.0_num-zint\n");
    fh.write(indent_2+"sz0(n, 1) = zint\n");
if (noz==3):
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"sz0(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_2+"sz0(n, 0) = 0.75_num-zintsq\n");
    fh.write(indent_2+"sz0(n, 1) = 0.5_num*(0.5_num+zint)**2\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute Bx on particle\n")
for ll in range(izmin0, izmax0+1):
        for jj in range(ixmin, ixmax+2):
            fh.write(indent_2+"bx(nn) = bx(nn) + sx(n,"+str(jj)+")*sz0(n,"+str(ll)+")*bxg(j"+plusi(jj)+",1,l0"+plusi(ll)+")\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute By on particle\n")
for ll in range(izmin0, izmax0+1):
        for jj in range(ixmin0, ixmax0+1):
            fh.write(indent_2+"by(nn) = by(nn) + sx0(n,"+str(jj)+")*sz0(n,"+str(ll)+")*byg(j0"+plusi(jj)+",1,l0"+plusi(ll)+")\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute Bz on particle\n")
for ll in range(izmin, izmax+2):
        for jj in range(ixmin0, ixmax0+1):
            fh.write(indent_2+"bz(nn) = bz(nn) + sx0(n,"+str(jj)+")*sz(n,"+str(ll)+")*bzg(j0"+plusi(jj)+",1,l"+plusi(ll)+")\n");

fh.write(indent_2+"END DO\n");
fh.write(indent_1+"END DO\n");
        
# _______________________________
# lower_order_in_v false
        
fh.write(indent_1+"ELSE\n");

fh.write(indent_2+"DO ip=1,np,lvect\n")
fh.write(indent_3+"DO n=1,MIN(lvect,np-ip+1)\n")
fh.write(indent_3+"\n")

fh.write(indent_4+"nn=ip+n-1\n")
fh.write(indent_4+"\n")

fh.write(indent_4+"x = (xp(nn)-xmin)*dxi\n");
fh.write(indent_4+"z = (zp(nn)-zmin)*dzi\n");
fh.write(indent_4+"\n")
fh.write(indent_4+"! Compute index of particle\n")
if(nox==2*(nox/2)):
      fh.write(indent_4+"j=nint(x)\n");
      fh.write(indent_4+"j0=floor(x)\n");
else:
      fh.write(indent_4+"j=floor(x)\n");
      fh.write(indent_4+"j0=floor(x-0.5_num)\n");
fh.write(indent_2+"\n")         
if(noz==2*(noz/2)):
      fh.write(indent_4+"l=nint(z)\n");
      fh.write(indent_4+"l0=floor(z)\n");
else:
      fh.write(indent_4+"l=floor(z)\n");
      fh.write(indent_4+"l0=floor(z-0.5_num)\n");

fh.write(indent_4+"\n")        
fh.write(indent_4+"xint=x-j\n");
fh.write(indent_4+"zint=z-l\n");
fh.write(indent_4+"\n")
fh.write(indent_4+"! Compute shape factors\n")

if(nox==1):
    fh.write(indent_4+"sx(n, 0) = 1.0_num-xint\n");
    fh.write(indent_4+"sx(n, 1) = xint\n");
if(nox==2):
    fh.write(indent_4+"xintsq = xint*xint\n");
    fh.write(indent_4+"sx(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_4+"sx(n, 0) = 0.75_num-xintsq\n");
    fh.write(indent_4+"sx(n, 1) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_4+"oxint = 1.0_num-xint\n");
    fh.write(indent_4+"xintsq = xint*xint\n");
    fh.write(indent_4+"oxintsq = oxint*oxint\n");
    fh.write(indent_4+"sx(n,-1) = onesixth*oxintsq*oxint\n");
    fh.write(indent_4+"sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
    fh.write(indent_4+"sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
    fh.write(indent_4+"sx(n, 2) = onesixth*xintsq*xint\n");

if(noz==1):
    fh.write(indent_4+"sz(n, 0) = 1.0_num-zint\n");
    fh.write(indent_4+"sz(n, 1) = zint\n");
if(noz==2):
    fh.write(indent_4+"zintsq = zint*zint\n");
    fh.write(indent_4+"sz(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_4+"sz(n, 0) = 0.75_num-zintsq\n");
    fh.write(indent_4+"sz(n, 1) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_4+"ozint = 1.0_num-zint\n");
    fh.write(indent_4+"zintsq = zint*zint\n");
    fh.write(indent_4+"ozintsq = ozint*ozint\n");
    fh.write(indent_4+"sz(n,-1) = onesixth*ozintsq*ozint\n");
    fh.write(indent_4+"sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
    fh.write(indent_4+"sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
    fh.write(indent_4+"sz(n, 2) = onesixth*zintsq*zint\n");

fh.write(indent_4+"\n");
fh.write(indent_4+"xint=x-0.5_num-j0\n");
fh.write(indent_4+"zint=z-0.5_num-l0\n");

fh.write(indent_4+"\n");
if(nox==1):
    fh.write(indent_4+"sx0( 0) = 1.0_num-xint\n");
    fh.write(indent_4+"sx0( 1) = xint\n");
if(nox==2):
    fh.write(indent_4+"xintsq = xint*xint\n");
    fh.write(indent_4+"sx0(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_4+"sx0(n, 0) = 0.75_num-xintsq\n");
    fh.write(indent_4+"sx0(n, 1) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_4+"oxint = 1.0_num-xint\n");
    fh.write(indent_4+"xintsq = xint*xint\n");
    fh.write(indent_4+"oxintsq = oxint*oxint\n");
    fh.write(indent_4+"sx0(n,-1) = onesixth*oxintsq*oxint\n");
    fh.write(indent_4+"sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
    fh.write(indent_4+"sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
    fh.write(indent_4+"sx0(n, 2) = onesixth*xintsq*xint\n");

fh.write(indent_4+"\n");
if(noz==1):
    fh.write(indent_4+"sz0(n, 0) = 1.0_num-zint\n");
    fh.write(indent_4+"sz0(n, 1) = zint\n");
if(noz==2):
    fh.write(indent_4+"zintsq = zint*zint\n");
    fh.write(indent_4+"sz0(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_4+"sz0(n, 0) = 0.75_num-zintsq\n");
    fh.write(indent_4+"sz0(n, 1) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_4+"ozint = 1.0_num-zint\n");
    fh.write(indent_4+"zintsq = zint*zint\n");
    fh.write(indent_4+"ozintsq = ozint*ozint\n");
    fh.write(indent_4+"sz0(n,-1) = onesixth*ozintsq*ozint\n");
    fh.write(indent_4+"sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
    fh.write(indent_4+"sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
    fh.write(indent_4+"sz0(n, 2) = onesixth*zintsq*zint\n");

fh.write(indent_4+"\n")
fh.write(indent_4+"! Compute Bx on particle\n")
for ll in range(izmin0, izmax0+1):
        for jj in range(ixmin, ixmax+2):
            fh.write(indent_4+"bx(nn) = bx(nn) + sx(n,"+str(jj)+")*sz0(n,"+str(ll)+")*bxg(j"+plusi(jj)+",1,l0"+plusi(ll)+")\n");

fh.write(indent_4+"\n")
fh.write(indent_4+"! Compute By on particle\n")
for ll in range(izmin0, izmax0+1):
        for jj in range(ixmin0, ixmax0+1):
            fh.write(indent_4+"by(nn) = by(nn) + sx0(n,"+str(jj)+")*sz0(n,"+str(ll)+")*byg(j0"+plusi(jj)+",1,l0"+plusi(ll)+")\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute Bz on particle\n")
for ll in range(izmin, izmax+2):
        for jj in range(ixmin0, ixmax0+1):
            fh.write(indent_2+"bz(nn) = bz(nn) + sx0(n,"+str(jj)+")*sz(n,"+str(ll)+")*bzg(j0"+plusi(jj)+",1,l"+plusi(ll)+")\n");

fh.write(indent_2+"END DO\n");
fh.write(indent_1+"END DO\n");

fh.write(indent_1+"END IF\n");

fh.write(indent_1+"RETURN\n");
fh.write("END SUBROUTINE "+subroutine_b_field+"\n");


fh.write(indent_1+"\n")
fh.write(indent_1+"\n")
fh.write(indent_1+"\n")
fh.write(indent_1+"\n")
fh.write(indent_1+"\n")
fh.write(indent_1+"\n")


# ------- EFIELD
# Write order n scalar gathering routine for BFIELD
# ------ Write electric field gathering/ energy conserving scheme  with OpenMP parallelization
fh.write("SUBROUTINE "+subroutine_e_field+"(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,       &\n");
fh.write(indent_cont_1+"dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &\n");
fh.write(indent_cont_1+"exg,eyg,ezg)\n");
fh.write("USE omp_lib\n");
fh.write("USE constants\n");
fh.write("IMPLICIT NONE\n");
fh.write("INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard\n");
fh.write("REAL(num), DIMENSION(np) :: xp,yp,zp,ex,ey,ez\n");
fh.write("REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg\n");
fh.write("REAL(num) :: xmin,ymin,zmin,dx,dy,dz\n");
fh.write("INTEGER(idp) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &\n");
fh.write(indent_cont_2+"ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll, j0, k0, l0\n");
fh.write("REAL(num) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &\n");
fh.write(indent_cont_2+"xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq\n");
fh.write("REAL(num), DIMENSION("+str(-int(nox/2))+":"+str(int((nox+1)/2))+") :: sx\n");
fh.write("REAL(num), DIMENSION("+str(-int(noy/2))+":"+str(int((noy+1)/2))+") :: sy\n");
fh.write("REAL(num), DIMENSION("+str(-int(noz/2))+":"+str(int((noz+1)/2))+") :: sz\n");
fh.write("REAL(num), DIMENSION(:), ALLOCATABLE :: sx0,sy0,sz0\n");
fh.write("REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num\n");
fh.write("\n");
fh.write("dxi = 1.0_num/dx\n");
fh.write("dyi = 1.0_num/dy\n");
fh.write("dzi = 1.0_num/dz\n");
fh.write("ALLOCATE(sx0("+str(ixmin0)+":"+str(ixmax0)+"),sz0("+str(izmin0)+":"+str(izmax0)+"))\n");
fh.write("sx=0.0_num\n");
fh.write("sz=0.0_num\n");
fh.write("sx0=0.0_num\n");
fh.write("sz0=0.0_num\n");

fh.write(indent_1+"IF (l_lower_order_in_v) THEN\n")
fh.write(indent_2+"DO ip=1,np,lvect\n")
fh.write(indent_3+"DO n=1,MIN(lvect,np-ip+1)\n")
fh.write(indent_3+"\n")

fh.write(indent_4+"nn=ip+n-1\n")

fh.write(indent_2+"\n")
fh.write(indent_2+"x = (xp(nn)-xmin)*dxi\n");
fh.write(indent_2+"z = (zp(nn)-zmin)*dzi\n");
fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute index of particle\n")
if (l_lower_order_in_v):
    if(nox==2*(nox/2)):
        fh.write(indent_2+"j=nint(x)\n");
        fh.write(indent_2+"j0=floor(x-0.5_num)\n");
    else:
        fh.write(indent_2+"j=floor(x)\n");
        fh.write(indent_2+"j0=floor(x)\n");
    if(noz==2*(noz/2)):
        fh.write(indent_2+"l=nint(z)\n");
        fh.write(indent_2+"l0=floor(z-0.5_num)\n");
    else:
        fh.write(indent_2+"l=floor(z)\n");
        fh.write(indent_2+"l0=floor(z)\n");
else:
    if(nox==2*(nox/2)):
        fh.write(indent_2+"j=nint(x)\n");
        fh.write(indent_2+"j0=floor(x)\n");
    else:
        fh.write(indent_2+"j=floor(x)\n");
        fh.write(indent_2+"j0=floor(x-0.5_num)\n");
    if(noz==2*(noz/2)):
        fh.write(indent_2+"l=nint(z)\n");
        fh.write(indent_2+"l0=floor(z)\n");
    else:
        fh.write(indent_2+"l=floor(z)\n");
        fh.write(indent_2+"l0=floor(z-0.5_num)\n");
fh.write(indent_2+"xint=x-j\n");
fh.write(indent_2+"zint=z-l\n");
fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute shape factors\n")
if(nox==1):
    fh.write(indent_2+"sx(n, 0) = 1.0_num-xint\n");
    fh.write(indent_2+"sx(n, 1) = xint\n");
if(nox==2):
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"sx(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
    fh.write(indent_2+"sx(n, 0) = 0.75_num-xintsq\n");
    fh.write(indent_2+"sx(n, 1) = 0.5_num*(0.5_num+xint)**2\n");
if (nox==3):
    fh.write(indent_2+"oxint = 1.0_num-xint\n");
    fh.write(indent_2+"xintsq = xint*xint\n");
    fh.write(indent_2+"oxintsq = oxint*oxint\n");
    fh.write(indent_2+"sx(n,-1) = onesixth*oxintsq*oxint\n");
    fh.write(indent_2+"sx(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
    fh.write(indent_2+"sx(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
    fh.write(indent_2+"sx(n, 2) = onesixth*xintsq*xint\n");

if(noz==1):
    fh.write(indent_2+"sz(n, 0) = 1.0_num-zint\n");
    fh.write(indent_2+"sz(n, 1) = zint\n");
if(noz==2):
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"sz(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
    fh.write(indent_2+"sz(n, 0) = 0.75_num-zintsq\n");
    fh.write(indent_2+"sz(n, 1) = 0.5_num*(0.5_num+zint)**2\n");
if(noz==3):
    fh.write(indent_2+"ozint = 1.0_num-zint\n");
    fh.write(indent_2+"zintsq = zint*zint\n");
    fh.write(indent_2+"ozintsq = ozint*ozint\n");
    fh.write(indent_2+"sz(n,-1) = onesixth*ozintsq*ozint\n");
    fh.write(indent_2+"sz(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
    fh.write(indent_2+"sz(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
    fh.write(indent_2+"sz(n, 2) = onesixth*zintsq*zint\n");
    
fh.write(indent_2+"\n");
fh.write(indent_2+"xint=x-0.5_num-j0\n");
fh.write(indent_2+"zint=z-0.5_num-l0\n");

if(l_lower_order_in_v):
    fh.write(indent_2+"\n"); 
    if(nox==1):
        fh.write(indent_2+"sx0(n, 0) = 1.0_num\n");
    if (nox==2):
        fh.write(indent_2+"sx0(n, 0) = 1.0_num-xint\n");
        fh.write(indent_2+"sx0(n, 1) = xint\n");
    if (nox==3):
        fh.write(indent_2+"xintsq = xint*xint\n");
        fh.write(indent_2+"sx0(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
        fh.write(indent_2+"sx0(n, 0) = 0.75_num-xintsq\n");
        fh.write(indent_2+"sx0(n, 1) = 0.5_num*(0.5_num+xint)**2\n");


    fh.write(indent_2+"\n");     
    if(noz==1):
        fh.write(indent_2+"sz0(n, 0) = 1.0_num\n");
    if (noz==2):
        fh.write(indent_2+"sz0(n, 0) = 1.0_num-zint\n");
        fh.write(indent_2+"sz0(n, 1) = zint\n");
    if (noz==3):
        fh.write(indent_2+"zintsq = zint*zint\n");
        fh.write(indent_2+"sz0(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
        fh.write(indent_2+"sz0(n, 0) = 0.75_num-zintsq\n");
        fh.write(indent_2+"sz0(n, 1) = 0.5_num*(0.5_num+zint)**2\n");   
        
else:
    fh.write(indent_2+"\n"); 
    if(nox==1):
        fh.write(indent_2+"sx0(n, 0) = 1.0_num-xint\n");
        fh.write(indent_2+"sx0(n, 1) = xint\n");
    if(nox==2):
        fh.write(indent_2+"xintsq = xint*xint\n");
        fh.write(indent_2+"sx0(n,-1) = 0.5_num*(0.5_num-xint)**2\n");
        fh.write(indent_2+"sx0(n, 0) = 0.75_num-xintsq\n");
        fh.write(indent_2+"sx0(n, 1) = 0.5_num*(0.5_num+xint)**2\n");
    if (nox==3):
        fh.write(indent_2+"oxint = 1.0_num-xint\n");
        fh.write(indent_2+"xintsq = xint*xint\n");
        fh.write(indent_2+"oxintsq = oxint*oxint\n");
        fh.write(indent_2+"sx0(n,-1) = onesixth*oxintsq*oxint\n");
        fh.write(indent_2+"sx0(n, 0) = twothird-xintsq*(1.0_num-xint*0.5_num)\n");
        fh.write(indent_2+"sx0(n, 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)\n");
        fh.write(indent_2+"sx0(n, 2) = onesixth*xintsq*xint\n");
        
    fh.write(indent_2+"\n");    
    if(noz==1):
        fh.write(indent_2+"sz0(n, 0) = 1.0_num-zint\n");
        fh.write(indent_2+"sz0(n, 1) = zint\n");
    if(noz==2):
        fh.write(indent_2+"zintsq = zint*zint\n");
        fh.write(indent_2+"sz0(n,-1) = 0.5_num*(0.5_num-zint)**2\n");
        fh.write(indent_2+"sz0(n, 0) = 0.75_num-zintsq\n");
        fh.write(indent_2+"sz0(n, 1) = 0.5_num*(0.5_num+zint)**2\n");
    if(noz==3):
        fh.write(indent_2+"ozint = 1.0_num-zint\n");
        fh.write(indent_2+"zintsq = zint*zint\n");
        fh.write(indent_2+"ozintsq = ozint*ozint\n");
        fh.write(indent_2+"sz0(n,-1) = onesixth*ozintsq*ozint\n");
        fh.write(indent_2+"sz0(n, 0) = twothird-zintsq*(1.0_num-zint*0.5_num)\n");
        fh.write(indent_2+"sz0(n, 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)\n");
        fh.write(indent_2+"sz0(n, 2) = onesixth*zintsq*zint\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute Ex on particle\n")
for ll in range(izmin, izmax+2):
        for jj in range(ixmin0, ixmax0+1):
            fh.write(indent_2+"ex(nn) = ex(nn) + sx0(n,"+str(jj)+")*sz(n,"+str(ll)+")*exg(j0"+plusi(jj)+",1,l"+plusi(ll)+")\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute Ey on particle\n")
for ll in range(izmin, izmax+2):
        for jj in range(ixmin, ixmax+2):
            fh.write(indent_2+"ey(nn) = ey(nn) + sx(n,"+str(jj)+")*sz(n,"+str(ll)+")*eyg(j"+plusi(jj)+",1,l"+plusi(ll)+")\n");

fh.write(indent_2+"\n")
fh.write(indent_2+"! Compute Ez on particle\n")
for ll in range(izmin0, izmax0+1):
        for jj in range(ixmin, ixmax+2):
            fh.write(indent_2+"ez(nn) = ez(nn) + sx(n,"+str(jj)+")*sz0(n,"+str(ll)+")*ezg(j"+plusi(jj)+",1,l0"+plusi(ll)+")\n");          

fh.write("END DO\n");
fh.write("END DO\n");
fh.write("DEALLOCATE(sx0,sz0)\n");
fh.write("RETURN\n");
fh.write("END SUBROUTINE "+subroutine_e_field+"\n");

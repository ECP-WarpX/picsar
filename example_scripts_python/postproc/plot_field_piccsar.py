#! /usr/bin/python
from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import glob 
import os
from PyLoadArrayPicssar import LoadBinNumPyArray3D


# Initialisation (chargement champ)
workdir='/Users/henrivincenti/warpcore3D_git_version/RESULTS';
os.chdir(workdir);
filey='ey0.pxr';
Ey=LoadBinNumPyArray3D(filey,91,91,91);

Eyplaneyz=abs(Ey[45,:,:]);
Eyplanexz=abs(Ey[:,45,:]);
Eyplanexy=abs(Ey[:,:,45]);

extent = [0., 91., 0., 91.]
fig = plt.figure(num=1,figsize=(3,3),dpi=80,facecolor='w',edgecolor='w')
plt.subplot(221)
imgplot=plt.imshow(Eyplaneyz,extent=extent,origin='lower',interpolation='nearest')
plt.gca().set_aspect('auto', adjustable='box')
plt.gca().set_title("y-z plane")
#imgplot.set_clim(-4,12)

plt.subplot(222)
imgplot=plt.imshow(Eyplanexz,extent=extent,origin='lower',interpolation='nearest')
plt.gca().set_aspect('auto', adjustable='box')
plt.gca().set_title("x-z plane")

plt.subplot(223)
imgplot=plt.imshow(Eyplanexy,extent=extent,origin='lower',interpolation='nearest')
plt.gca().set_aspect('auto', adjustable='box')
plt.gca().set_title("x-y plane")

Eyline1=Ey[45,45,:];
plt.subplot(224)
plt.plot(Eyline1)
#plt.ylim([0., 1000000000000])
plt.show()
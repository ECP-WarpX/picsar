#! /usr/bin/python
from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import glob 
import os
from PyLoadArrayPicssar import LoadBinNumPyArray3D
from numpy import linalg as LA

# Paramaters -constants
eps0=8.85e-12

# Initialisation (chargement champ)
it=1
workdir='/Users/henrivincenti/warpcore3D_git_version/RESULTS/';
#workdir='/Users/henrivincenti/PICSAR_JLMOD/RESULTS/';
os.chdir(workdir);
filedive=workdir+'dive'+str(it)+'.pxr';
filerho=workdir+'rho'+str(it)+'.pxr';

dive=LoadBinNumPyArray3D(filedive,91,91,91);
rho=LoadBinNumPyArray3D(filerho,91,91,91);

rhoc1=rho[:,45,:];
divec1=dive[:,45,:];
extent = [0., 91., 0., 91.]
fig = plt.figure(num=1,figsize=(3,3),dpi=80,facecolor='w',edgecolor='w')
plt.subplot(121)
imgplot=plt.imshow(rhoc1,extent=extent,origin='lower',interpolation='nearest')
plt.gca().set_aspect('auto', adjustable='box')
plt.colorbar()

epsilon=rhoc1-divec1*eps0
plt.subplot(122)
imgplot=plt.imshow(epsilon,extent=extent,origin='lower',interpolation='nearest')
plt.gca().set_aspect('auto', adjustable='box')
plt.colorbar()
#imgplot.set_clim(-19,-14)
#plt.ylim([0., 1000000000000])



print("Differences norme L2 ||rho-dive|| iteration it = " + str(it))
print(LA.norm((dive*eps0-rho)))


print("Total charge " + str(it))
print(np.sum(rho))
print("Total dive" + str(it))
print(np.sum(dive*eps0))
plt.show()
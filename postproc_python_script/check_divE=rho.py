"""
 _______________________________________________________________________________

 *** Copyright Notice ***

 "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)  
 2016, The Regents of the University of California, through Lawrence Berkeley 
 National Laboratory (subject to receipt of any required approvals from the 
 U.S. Dept. of Energy). All rights reserved.

 If you have questions about your rights to use or distribute this software, 
 please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

 NOTICE.
 This Software was developed under funding from the U.S. Department of Energy 
 and the U.S. Government consequently retains certain rights. As such, the U.S. 
 Government has been granted for itself and others acting on its behalf a  
 paid-up, nonexclusive, irrevocable, worldwide license in the Software to 
 reproduce, distribute copies to the public, prepare derivative works, and 
 perform publicly and display publicly, and to permit other to do so. 


 Script to check divE=rho using PICSAR field outputs (divE and rho).
 
 Developers:
 Henri Vincenti
 
 Date:
 Creation 2015

 _______________________________________________________________________________

"""

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
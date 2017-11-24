from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import sys,getopt
import sys
#from PyLoadArrayPicsar import LoadBinNumPyArray3D

from Field import *

# ______________________________________________________________________________
# RCparams

mpl.rcParams['font.size'] = 14
mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.subplot.wspace'] = 0.3
mpl.rcParams['figure.subplot.hspace'] = 0.3
# ______________________________________________________________________________
# Main


 # ________________________________________________________________________________
 # Reading of the output arguments

xslice = None
yslice = None
zslice = None

#  print 'Input file is ', file
file = sys.argv[1]
 # ________________________________________________________________________________
  # Creation of an object field from the file
field=Field(file)
plt.figure()
plt.grid(True)
nx=field.nx
nz=field.nz
xx=np.linspace(field.xmin,field.xmax,nx)

zz=np.linspace(field.zmin,field.zmax,nz)
X,Z=np.meshgrid(xx,zz)
ff=field.f[:,0,:]
plt.pcolormesh(X,Z,ff)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('z')
plt.xlim(field.xmin,field.xmax)
plt.ylim(field.zmin,field.zmax)
plt.title(file)
plt.ion()
plt.show()
	 


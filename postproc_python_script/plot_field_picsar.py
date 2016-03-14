#! /usr/bin/python
from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import glob 
import os
import sys,getopt
from PyLoadArrayPicssar import LoadBinNumPyArray3D


# Initialisation (chargement champ)
#workdir='/Users/henrivincenti/warpcore3D_git_version/RESULTS';
#os.chdir(workdir);
#filey='ey0.pxr';

# Reading of the output arguments

def main(argv):
  try:
    opts, args = getopt.getopt(argv,"hf:n:",["file="])
  except getopt.GetoptError:
    print 'Specify file name'
    sys.exit(2)
    
  filey = argv[0]
    
  for opt, arg in opts:
    if opt == '-h':
      print 'Help page'
      print '-f file name'
      print '-n dimension nx=ny=nz=n of the case'
      sys.exit()
    elif opt in ("-f", "--file"):
      filey = arg
    elif opt in "-n":
      print arg
      nx = int(arg)
      ny = int(arg)
      nz = int(arg)
            
  print 'Input file is ', filey

  Ey=LoadBinNumPyArray3D(filey,nx,ny,nz);

  x = linspace(0,100.,nx)
  y = linspace(0,100.,ny)
  z = linspace(0,100.,nz)

  Eyplaneyz=abs(Ey[nx/2,:,:]);
  Eyplanexz=abs(Ey[:,ny/2,:]);
  Eyplanexy=abs(Ey[:,:,nz/2]);

  extent = [0., 100., 0., 100.]
  fig = plt.figure(num=1,figsize=(10,8),dpi=80,facecolor='w',edgecolor='w')
  plt.subplot(221)
  
  im = plt.pcolormesh(y,z,Eyplaneyz)
  plt.title("y-z plane")
  plt.xlabel('y')
  plt.ylabel('z')
  cb = plt.colorbar()
  
  #imgplot=plt.imshow(Eyplaneyz,extent=extent,origin='lower',interpolation='nearest')
  #plt.gca().set_aspect('auto', adjustable='box')
  #plt.gca().set_title("y-z plane")
  #imgplot.set_clim(-4,12)

  plt.subplot(222)
  im = plt.pcolormesh(x,z,Eyplanexz)
  plt.title("x-z plane")
  cb = plt.colorbar()
  
  #imgplot=plt.imshow(Eyplanexz,extent=extent,origin='lower',interpolation='nearest')
  #plt.gca().set_aspect('auto', adjustable='box')
  #plt.gca().set_title("x-z plane")

  plt.subplot(223)
  im = plt.pcolormesh(x,y,Eyplanexz)
  plt.title("x-y plane")
  cb = plt.colorbar()
  
  #imgplot=plt.imshow(Eyplanexy,extent=extent,origin='lower',interpolation='nearest')
  #plt.gca().set_aspect('auto', adjustable='box')
  #plt.gca().set_title("x-y plane")

  Eyline1=Ey[nx/2,ny/2,:];
  plt.subplot(224)
  plt.plot(Eyline1)
  #plt.ylim([0., 1000000000000])
  plt.xlabel("z")
  plt.show()
  
if __name__ == "__main__":
   main(sys.argv[1:])  

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
 
 
 Field.py
 
 This file contains the class Field to read and plot Picsar field diagnostics.
 
 Developers:
 Mathieu Lobet
 Remi Lehe
 
 Date:
 Creation 2016
 _______________________________________________________________________________
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import struct
import matplotlib.gridspec as gridspec

# ______________________________________________________________________________

class Field:
  """
  Class Field
  
  This class enable to manage field outputs from the mini App
  
  """
  def __init__(self,filename):  

    with open(filename,'rb') as file:
      fileContent = file.read()
  
    l = 0
    self.xmin = struct.unpack("d", fileContent[l:l+8])[0]; l+=8
    self.xmax = struct.unpack("d", fileContent[l:l+8])[0]; l+=8  
    self.nx = struct.unpack("i", fileContent[l:l+4])[0]; l+=4

    self.ymin = struct.unpack("d", fileContent[l:l+8])[0]; l+=8
    self.ymax = struct.unpack("d", fileContent[l:l+8])[0]; l+=8  
    self.ny = struct.unpack("i", fileContent[l:l+4])[0]; l+=4
    
    self.zmin = struct.unpack("d", fileContent[l:l+8])[0]; l+=8
    self.zmax = struct.unpack("d", fileContent[l:l+8])[0]; l+=8  
    self.nz = struct.unpack("i", fileContent[l:l+4])[0]; l+=4  
    
    print ('Openning of ',filename)
    print ('xmin: %f, xmax: %f, nx: %d'%(self.xmin,self.xmax,self.nx))
    print ('ymin: %f, ymax: %f, ny: %d'%(self.ymin,self.ymax,self.ny))
    print ('zmin: %f, zmax: %f, nz: %d'%(self.zmin,self.zmax,self.nz))

    self.x = np.linspace(self.xmin,self.xmax,self.nx, endpoint=False)
    self.y = np.linspace(self.ymin,self.ymax,self.ny, endpoint=False)
    self.z = np.linspace(self.zmin,self.zmax,self.nz, endpoint=False)
    
    self.dx = self.x[1] - self.x[0]
    if(self.ny > 1):
	self.dy = self.y[1] - self.y[0]
    else:
	#2d case
	self.dy = 1
    self.dz = self.z[1] - self.z[0]
            
    self.f = np.zeros([self.nz,self.ny,self.nx])  
    
    ncells = self.nx*self.ny*self.nz
    print ('Total number of cells:',ncells)
    print ('Size of file:',len(fileContent))
    
    for iz in range(self.nz):
      for iy in range(self.ny):      
        for ix in range(self.nx):      
          self.f[iz,iy,ix] = struct.unpack("d", fileContent[l:l+8])[0]; l+=8

  def plot_3d(self,ax,ixslices=[],iyslices=[],izslices=[]):
    
    levels = np.linspace(-1, 1, 40)
    
    if len(ixslices)>0:
      for i in range(len(ixslices)):
        ixs = ixslices[i]
        X, Y = np.meshgrid(self.y, self.z)
        ax.contourf(X,Y,self.x[ixs] + self.f[:,:,ixs], xdir='x')

    if len(iyslices)>0:
      for i in range(len(iyslices)):
        iys = iyslices[i]
        X, Y = np.meshgrid(self.x, self.z)
        ax.contourf(X,Y,self.y[iys] + self.f[:,iys,:], ydir='y')

    if len(izslices)>0:
      for i in range(len(izslices)):
        izs = izslices[i]
        X, Y = np.meshgrid(self.x, self.y)
        ax.contourf(X,Y,self.z[izs] + self.f[:,izs,:], zdir='z')

    ax.set_xlim3d(self.xmin, self.xmax)
    ax.set_ylim3d(self.ymin, self.ymax)
    ax.set_zlim3d(self.zmin, self.zmax)
    
  def plot(self,fig,ax,slice='x',slice_index=0,slice_pos=None):

    if slice_pos==None:
      index = slice_index
    else:
      if slice=='x':
        index = int((slice_pos-self.xmin)/self.dx)
      if slice=='y':
        index = int((slice_pos-self.ymin)/self.dy)
      if slice=='z':
        index = int((slice_pos-self.zmin)/self.dz)
    
      print ('index:',index)
    
    if slice=='x':

      im = ax.pcolormesh(self.y,self.z,self.f[:,:,index])
      ax.set_xlabel('y')
      ax.set_ylabel('z')   
      ax.set_title('slice at x=%f'%(self.x[index]))   
      ax.set_xlim([self.ymin,self.ymax])
      ax.set_ylim([self.zmin,self.zmax])
      cb = plt.colorbar(im,ax=ax)

    if slice=='y':

      im = ax.pcolormesh(self.x,self.z,self.f[:,index,:])
      cb = plt.colorbar(im,ax=ax)
      ax.set_xlabel('x')
      ax.set_ylabel('z')  
      ax.set_xlim([self.xmin,self.xmax])
      ax.set_ylim([self.zmin,self.zmax])       
      ax.set_title('slice at y=%f'%(self.y[index]))        
      
    if slice=='z':

      im = ax.pcolormesh(self.x,self.y,self.f[index,:,:])
      cb = plt.colorbar(im,ax=ax) 
      ax.set_xlabel('x')
      ax.set_ylabel('y')   
      ax.set_xlim([self.xmin,self.xmax])
      ax.set_ylim([self.ymin,self.ymax])      
      ax.set_title('slice at z=%f'%(self.z[index])) 

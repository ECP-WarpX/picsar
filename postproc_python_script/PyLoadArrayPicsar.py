#! /usr/bin/python
from numpy import *
import numpy as np
import struct

def LoadBinNumPyArray3D(filename,dimx,dimy,dimz):
    """
    Read Picsar output arrays
    """
    fd=open(filename,'rb');
    dim=np.zeros(3,dtype=np.int32)
    dim[0]=dimx;
    dim[1]=dimy;
    dim[2]=dimz;
    Mat3D = np.transpose(np.reshape(np.fromfile(fd,dtype=np.float64,count=np.prod(dim)),dim));
    return Mat3D


def read_picsar_temporal_diags(filename):
  """
  Read Picsar temporal files
  """
  with open(filename,'rb') as file:
    fileContent = file.read()
  
  l = 0
  nc = struct.unpack("i", fileContent[l:l+4])[0]; l+=4
  dt = struct.unpack("d", fileContent[l:l+8])[0]; l+=8
  nl = (len(fileContent)-12)/nc/8
  
  print
  print ' Reading',filename
  print ' Total number of bytes',len(fileContent)
  print ' Number of columns:',nc
  print ' Number of iterations:',nl
  print ' Time step:',dt,'s'
  
  if nc > 1:
    array = zeros([nc,nl])
    t = linspace(0,nl*dt,nl)
    
    for i in range(nl):
      for j in range(nc):
        array[j,i] = struct.unpack("d", fileContent[l:l+8])[0]; l+=8
  else:
    array = zeros(nl)
    t = linspace(0,nl*dt,nl)
    
    for i in range(nl):
        array[i] = struct.unpack("d", fileContent[l:l+8])[0]; l+=8   
        
  return t,array   



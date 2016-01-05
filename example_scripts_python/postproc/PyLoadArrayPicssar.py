#! /usr/bin/python
from numpy import *
import numpy as np

def LoadBinNumPyArray3D(filename,dimx,dimy,dimz):
    fd=open(filename,'rb');
    dim=np.zeros(3,dtype=np.int32)
    dim[0]=dimx;
    dim[1]=dimy;
    dim[2]=dimz;
    Mat3D = np.transpose(np.reshape(np.fromfile(fd,dtype=np.float64,count=np.prod(dim)),dim));
    return Mat3D




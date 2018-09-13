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


 plot_field_picsar.py

 Script that can be called to plot field outputs.
 
 This script uses command arguments:
 -h help information
 -f file path
 -n dimension of the case
 --xslice= to make a slice at the position x=
 --yslice= to make a slice at the position y=
 --zslice= to make a slice at the position z=
  
 Enter: `python plot_field_picsar.py -h`
 to have the help information.
 
 Developers:
 Mathieu Lobet
 
 Date:
 Creation 2016

 _______________________________________________________________________________

"""

#! /usr/bin/python
from numpy import *
import matplotlib.pyplot as plt
import numpy as np
import glob 
import os
import sys,getopt
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

def main(argv):

 # ________________________________________________________________________________
 # Reading of the output arguments
  try:
    opts, args = getopt.getopt(argv,"hf:n:",["file=","xslice=","yslice=","zslice="])
  except getopt.GetoptError:
    print 'Specify file name'
    sys.exit(2)
    
  filey = argv[0]
  xslice = None
  yslice = None
  zslice = None
    
  for opt, arg in opts:
    if opt == '-h':
      print 'Help page'
      print '-f file name'
      print '-n dimension nx=ny=nz=n of the case'
      sys.exit()
    elif opt in ("-f", "--file"):
      file = arg
    elif opt in ("--xslice"): 
      xslice = float(arg)     
    elif opt in ("--yslice"): 
      yslice = float(arg) 
    elif opt in ("--zslice"): 
      zslice = float(arg)                      
  print 'Input file is ', file

  if xslice == None:
    xslice = 0
  if yslice == None:
    yslice = 0    
  if zslice == None:
    zslice = 0
 # ________________________________________________________________________________
  # Creation of an object field from the file
  field=Field(file)

  fig = plt.figure(figsize=(13,9))
  gs = gridspec.GridSpec(2, 2)
  ax0 = plt.subplot(gs[0, 0])
  ax1 = plt.subplot(gs[1, 0])
  ax2 = plt.subplot(gs[0, 1])

  field.plot(fig,ax0,slice='x',slice_pos=xslice)
  field.plot(fig,ax1,slice='y',slice_pos=yslice)
  field.plot(fig,ax2,slice='z',slice_pos=zslice)

  fig.tight_layout()
  
  plt.show()
  
if __name__ == "__main__":
   main(sys.argv[1:])  

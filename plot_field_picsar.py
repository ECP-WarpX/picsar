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
#  print 'Input file is ', file
  file = "RESULTS/"+filey+".pxr"
  if xslice == None:
    xslice = 0e-6
  if yslice == None:
    yslice = 0
  if zslice == None:
    zslice = 0e-7
 # ________________________________________________________________________________
  # Creation of an object field from the file
  field=Field(file)

  fig = plt.figure(figsize=(13,9))
  gs = gridspec.GridSpec(2, 2)
  ax0 = plt.subplot(gs[0, 0])
  ax1 = plt.subplot(gs[1, 0])
  ax2 = plt.subplot(gs[0, 1])
  ax3 = plt.subplot(gs[1,1])
  field.plot(fig,ax0,slice='x',slice_pos=xslice)
  field.plot(fig,ax1,slice='y',slice_pos=yslice)
  field.plot(fig,ax2,slice='z',slice_pos=zslice)
  field.plot1d(fig,ax3,'x','y',yslice,zslice)
  fig.tight_layout()
  plt.ion()
  plt.show()
	 
if __name__ == "__main__":
   main(sys.argv[1:])


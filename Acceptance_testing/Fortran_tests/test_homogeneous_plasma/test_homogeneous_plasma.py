#! /usr/bin/python
# ____________________________________________________________________
#
# Test Langmuir wave
# Fortran version
# ____________________________________________________________________


from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import glob 
import os
import sys,getopt
import struct
from PyLoadArrayPicsar import *
from numpy import linalg as LA
from subprocess import call

def test_homogeneous_plasma(trun,ttest,tpath):
  """
  
  """

  print ' _________________________________________'
  print ''
  print ' Test homogeneous plasma'  
  print ' _________________________________________'

  trun=int(trun)
  ttest = int(ttest)

  print
  print ' Running simulation:',trun
  print ' Using assert:',ttest
  print ' Run in path:',tpath
  print

  # ____________________________________________________________________
  # RCparams
  mpl.rcParams['font.size'] = 15
  mpl.rcParams['legend.fontsize'] = 12
  mpl.rcParams['figure.facecolor'] = 'white'
  
  # ____________________________________________________________________
  # Computation

  current_path = os.getcwd()
  file_path = ''
  for arg in sys.argv:
    if 'test_langmuir_wave.py' in arg:
      file_path = os.path.dirname(os.path.realpath(arg))
  
  print ' Current directory:',current_path
  print ' Script directory:',file_path
  print 

  if (trun):
    
    # Change working directory for the simulation one 
    if len(tpath)>0:
      os.chdir(tpath)
    
    input_file = '../../../example_decks_fortran/homogeneous_plasma.pixr'
    picsar_exe = '../../../fortran_bin/picsar'
    
    # Get the input file from the example
    call(["cp", "-p",input_file,'input_file.pixr']) 
    call(["cp", "-p",picsar_exe,'.']) 
    
    # Run picsar
    #omp_num_threads = 2
    #call(["export", "OMP_NUM_THREADS=2"]) 
    os.putenv('OMP_NUM_THREADS','2')
    call(["rm","RESULTS/*"])
    call(["mpirun","-n","4","./picsar"])
    #call(["sh","launcher"])
  
  # Parameters
  echarge = 1.60217662E-19
  emass = 9.10938356E-31
  eps0 = 8.85418782E-12
  n = 1e25
  v = 0.1
  gamma = 1./sqrt(1.-v**2)
  wplasma = echarge*sqrt(n/(emass*eps0))
  Tplasma = 2.*pi*sqrt(gamma)/wplasma
  print ' Gamma:',gamma 
  print ' Plasma frequency:',wplasma
  print ' Plasma period:',Tplasma,'s'
  
  # Opening of the temporal files
  print
  print ' _________________________________________'
  print ' Checking energy balance'
  # Kinetic energy
  t,kinE = read_picsar_temporal_diags('RESULTS/kinE')

  # Opening of the ezfield
  t,ezE = read_picsar_temporal_diags('RESULTS/ezE')

  # Opening of the eyfield
  t,eyE = read_picsar_temporal_diags('RESULTS/eyE')

  # Opening of the exfield
  t,exE = read_picsar_temporal_diags('RESULTS/exE')

  # Opening of the exfield
  t,bzE = read_picsar_temporal_diags('RESULTS/bzE')

  # Opening of the exfield
  t,byE = read_picsar_temporal_diags('RESULTS/byE')

  # Opening of the exfield
  t,bxE = read_picsar_temporal_diags('RESULTS/bxE')
  
  total_energy = sum(kinE,axis=0) + ezE + exE + eyE + bzE + bxE + byE
      
  min_totalE = min(total_energy)
  max_totalE = max(total_energy)
  diffrel = (max_totalE - min_totalE)/max_totalE
  print ' Relative error on the total energy:',diffrel
  if ttest: assert diffrel < 1e-2
      
  # Plotting      
  fig = plt.figure(figsize=(12,8))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :])
  
  ax.plot(t,kinE[0,:],label='Electron kinetic energy')
  ax.plot(t,kinE[1,:],label='Proton kinetic energy')
  ax.plot(t,ezE,label='Ez energy')
  ax.plot(t,eyE,label='Ey energy')  
  ax.plot(t,exE,label='Ex energy')
  ax.plot(t,bzE,label='Bz energy')
  ax.plot(t,byE,label='By energy',color='brown')  
  ax.plot(t,bxE,label='Bx energy',color='orange')  
  ax.plot(t,total_energy,label='tot',color='k',ls=':')
  ax.plot([Tplasma,Tplasma],[0.,max(total_energy)],ls='--',color='k')
  
  ax.set_xlabel('t (s)')
  ax.set_ylabel('Energy (J)')
  
  plt.annotate('', xy=(0, kinE[0,0]*0.5), xycoords='data',xytext=(Tplasma, kinE[0,0]*0.5), textcoords='data',arrowprops={'arrowstyle': '<->'})
  
  ax.legend(loc='upper center',ncol=3,borderaxespad=-3)
  
  # Test divergence
  print ' ________________________________________ '
  print ' Check DivE = rho/eps0'
  
  # Opening of divE-rho
  t,diverho = read_picsar_temporal_diags('RESULTS/divE-rho')

  fig1 = plt.figure(figsize=(12,8))
  gs1 = gridspec.GridSpec(2, 2)
  ax1 = plt.subplot(gs1[:, :])
  
  ax1.plot(t,diverho,label=r'$\nabla E \times \varepsilon_0 - \rho$',lw=2) 
  ax1.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
  
  ax1.set_xlabel('t (s)')
  
  print ' _________________________________ '
  print ' Check DivE = rho/eps0'
  for it in range(0,50,10):
    dive=LoadBinNumPyArray3D('RESULTS/dive' + str(it) + '.pxr',100,100,100);
    rho=LoadBinNumPyArray3D('RESULTS/rho'+ str(it) + '.pxr',100,100,100);  
    norm = LA.norm((dive*eps0-rho)) 
    print
    print(" Differences norme L2 ||rho-divE|| iteration it = " + str(it))
    print "",LA.norm((dive*eps0-rho))
    print " Total charge "
    print "",np.sum(rho)
    print " Total divergence at "
    print "",np.sum(dive*eps0)
    if ttest: assert (norm < 1E-5),"L2 norm||DivE - rho/eps0|| too high"

  if ttest: assert (max(diverho) < 1E-5),"L2 norm||DivE - rho/eps0|| too high"

  # ____________________________________________________
  # Advice
  
  print
  print ' _______________________________________'
  print ' Advice for users:' 
  print ' - Check that the energy is constant with time'
  print ' - Check that divE = rho/eps0 for each tests'
  print
  
  plt.show()

if __name__ == "__main__":

  argv = sys.argv[1:]
  run = 1
  test = 1
  path = ''

  try:
    opts, args = getopt.getopt(argv,"hr:t:p:",["test=","run=",'path='])
  except getopt.GetoptError:
    help()
    sys.exit(2)
    
  for opt, arg in opts:
    if opt == '-h':
      help()
      sys.exit()
    elif opt in ("-t", "--test"):
      test = int(arg)
    elif opt in ("-r", "--run"):
      run = int(arg)
    elif opt in ("-p", "--path"):
      path = arg
  test_homogeneous_plasma(run,test,path)  



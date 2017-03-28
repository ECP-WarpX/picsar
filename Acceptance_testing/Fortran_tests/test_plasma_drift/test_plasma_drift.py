#! /usr/bin/python

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

 test_plasma_drift.py

 Script to launch the plasma drift test using Fortran
 - analysis of the energy
 - analysis of the beam trajectories
 - analysis of divE = rho/eps0, total divergence and total charge
 
 Developer:
 Mathieu Lobet
 
 Date
 Creation 2016
 _______________________________________________________________________________
"""


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
from Field import *

# _______________________________________________________________________
# Functions

def test_plasma_drift(tpath,trun,ttest,tshow):
  """
  Main function to manage the simulation and the analysis.
  """

  print (' _________________________________________')
  print ('')
  print (' Test plasma drift')
  print (' _________________________________________')

  trun=int(trun)
  ttest = int(ttest)
  tshow= int(tshow)

  print
  print (' Running simulation:',trun)
  print (' Using assert:',ttest)
  print (' Run in path:',tpath)
  print (' Show the results in figures:',tshow)
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
  
  print(' Current directory:',current_path)
  print(' Script directory:',file_path)
  print
  
  if (trun):
    
    # Change working directory for the simulation one 
    if len(tpath)>0:
      os.chdir(path)
    
    input_file = '../../../examples/example_decks_fortran/drifted_plasma.pixr'
    picsar_exe = '../../../fortran_bin/picsar'
    
    # Get the input file from the example
    call(["cp", "-p",input_file,'input_file.pixr']) 
    call(["cp", "-p",picsar_exe,'.']) 
    
    # Run picsar
    #omp_num_threads = 2
    #call(["export", "OMP_NUM_THREADS=2"]) 
    os.putenv('OMP_NUM_THREADS','2')
    call(["rm","RESULTS/*"])
    call(["mkdir","-p","RESULTS"])
    call(["mpirun","-n","4","./picsar"])
    #call(["sh","launcher"])
    # Run picsar
    
    #call(["sh","launcher"])
  
  # ____________________________________________________________________
  # Analysis
  
  # Parameters
  echarge = 1.60217662E-19
  emass = 9.10938356E-31
  eps0 = 8.85418782E-12
  n = 1e25
  v = 0.1
  gamma = 1./sqrt(1.-v**2)
  wplasma = echarge*sqrt(n/(emass*eps0))
  Tplasma = 2.*pi*sqrt(gamma)/wplasma
  
  print(' Gamma:',gamma)
  print(' Plasma frequency:',wplasma)
  print(' Plasma period:',Tplasma,'s')
  
  # Test energy
  # Opening of the temporal files
  print
  print (' _________________________________________')
  print (' Checking energy balance')
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
  print
  print(' Relative error on the total energy:',diffrel)
  if ttest: assert diffrel < 1e-2
      
  # Plotting 
  if tshow:     
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2)
    ax = plt.subplot(gs[:, :])
  
    for sp in range(len(kinE)):
      ax.plot(t,kinE[sp,:],label='Species %d'%sp)
    
    ax.plot(t,ezE,label='Ez energy')
    ax.plot(t,eyE,label='Ey energy')  
    ax.plot(t,exE,label='Ex energy')
    ax.plot(t,bzE,label='Bz energy')
    ax.plot(t,byE,label='By energy',color='brown')  
    ax.plot(t,bxE,label='Bx energy',color='orange')  
    ax.plot(t,total_energy,label='tot',color='k',ls=':')
    ax.plot([Tplasma,Tplasma],[0.,max(total_energy)],ls='--',color='k')
  
    ax.set_ylim([0,max(total_energy)*1.1])
  
    ax.set_xlabel('t (s)')
    ax.set_ylabel('Energy (J)')
  
    plt.annotate('', xy=(0, kinE[0,0]*0.5), xycoords='data',xytext=(Tplasma, kinE[0,0]*0.5), textcoords='data',arrowprops={'arrowstyle': '<->'})
  
    ax.legend(loc='upper center',ncol=5,borderaxespad=-3)
 
  # Test divergence
  print (' ________________________________________ ')
  print (' Check DivE = rho/eps0')
  
  # Opening of divE-rho
  t,diverho = read_picsar_temporal_diags('RESULTS/divE-rho')
  t,rho = read_picsar_temporal_diags('RESULTS/rho')
  t,dive= read_picsar_temporal_diags('RESULTS/divE')
  
  print()
  print ('max(||diverho||):',max(diverho))
  print ('max(||rho||(t)):',max(rho))
  print ('min(||rho||(t)):',min(rho)) 
  print ('max(||divE||(t)):',max(dive)*eps0)
  print ('min(||divE||(t)):',min(dive)*eps0)
  print()
  
  if tshow:
    fig1 = plt.figure(figsize=(12,8))
    gs1 = gridspec.GridSpec(7, 4)
    ax1 = plt.subplot(gs1[0:3, :])
    ax2 = plt.subplot(gs1[4:7, :])
  
    ax1.plot(t,diverho,label=r'$||\nabla E \times \varepsilon_0 - \rho||$',lw=2) 
    ax1.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
    ax1.set_xlabel('t (s)')

    ax2.plot(t,dive*eps0,label=r'$|| \nabla E \times \varepsilon_0 ||$',lw=2) 
    ax2.plot(t,rho,label=r'$|| \rho ||$',color='r',lw=2,ls='--')
    ax2.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
    ax2.set_xlabel('t (s)')

  #if ttest: assert (max(diverho) < 1E-3),"L2 norm||DivE - rho/eps0|| too high"

  # Analyse of the files
  if 1: # Temporarily removed due to MPI-IO issues (plateform dependent)
      for it in range(0,21,10):
        dive=Field('RESULTS/dive' + str(it) + '.pxr')
        rho=Field('RESULTS/rho'+ str(it) + '.pxr')  
        F = ((dive.f*eps0-rho.f)) 
        min_F = amin(abs(F))
        max_F = amax(abs(F))
        ave_F = average(abs(F))
        
        print()
        print(" Iteration it = " + str(it))
        print(" Total charge:",sqrt(sum(rho.f**2)))
        print(" Total divergence:",sqrt(sum(dive.f**2))*eps0)
        print(" min(divE*eps0-rho)):",min_F)
        print(" max(divE*eps0-rho)):",max_F)
        print(" ave(divE*eps0-rho)):",ave_F)
        
        #if ttest: assert (max_F < 1E-3),"L2 norm||DivE*eps0 - rho|| too high"

  # ____________________________________________________
  # Advice
  
  print
  print (' _______________________________________')
  print (' Advice for users:' )
  print (' - Check that the energy keeps constant with time')
  print (' - Check that divE = rho/eps0 for each tests')
  print
  
  if tshow: plt.show()

if __name__ == "__main__":

  argv = sys.argv[1:]
  run = 1
  test = 1
  show = 1
  path = ''

  try:
    opts, args = getopt.getopt(argv,"hr:t:p:w:",["test=","run=",'path=','show='])
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
    elif opt in ("-w", "--show"):
      show = int(arg)
  test_langmuir_wave(path,run,test,show)  


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

 test_langmuir_wave.py

 Script to launch the Langmuir wave test using Fortran
 
 Developer:
 Mathieu Lobet
 
 Date
 Creation 2016
 _______________________________________________________________________________
"""

# _______________________________________________________________________
# Functions

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
import pytest
from Field import *

def help():
    print (' Help:')
    print (' -t --test: if 1 use assert')
    print (' -r --run: if 1 run the code before analysis')
    print (' -w --show: if 1 show the results with Matplotlib')
    print (' -f --folder: path to the test folder')


def test_langmuir_wave(tpath,trun,ttest,tshow):
  """
  Function to launch and analyse the Langmuir Wave test case
  """
  
  print
  print (' _________________________________________')
  print ('')
  print (' Test Langmuir Wave')
  print (' _________________________________________')

  trun=int(trun)
  ttest = int(ttest)
  tshow=int(tshow)

  print
  print (' Running simulation:',trun)
  print (' Using assert:',ttest)
  print (' Run in path:',tpath)
  print (' Show results with Matplotlib:',tshow)
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
  
  print (' Current directory:',current_path)
  print (' Script directory:',file_path)
  print

  if (trun):
    
    # Change working directory for the simulation one 
    if len(tpath)>0:
      os.chdir(path)
    
    input_file = '../../../examples/example_decks_fortran/langmuir_wave.pixr'
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

  # ____________________________________________________________________
  # Analysis
  
  # Parameters
  echarge = 1.60217662E-19
  emass = 9.10938356E-31
  eps0 = 8.85418782E-12
  n = 1.1e25
  v = 0.1
  gamma = 1./sqrt(1.-v**2)
  wplasma = echarge*sqrt(n/(emass*eps0))
  Tplasma = 2.*pi*sqrt(gamma)/wplasma
  print (' Gamma:',gamma)
  print (' Plasma frequency:',wplasma)
  print (' Plasma period:',Tplasma,'s')
  
  # Opening of the temporal files
  print
  print (' _________________________________________')
  print (' Checking energy balance')
  # Kinetic energy
  t,kinE = read_picsar_temporal_diags('RESULTS/kinE')
  kinE = kinE[:,1:]

  # Opening of the ezfield
  t,ezE = read_picsar_temporal_diags('RESULTS/ezE')
  ezE = (ezE[1:] + ezE[0:-1])*0.5

  # Opening of the eyfield
  t,eyE = read_picsar_temporal_diags('RESULTS/eyE')
  eyE = (eyE[1:] + eyE[0:-1])*0.5

  # Opening of the exfield
  t,exE = read_picsar_temporal_diags('RESULTS/exE')
  exE = (exE[1:] + exE[0:-1])*0.5
  
  # Opening of the exfield
  t,bzE = read_picsar_temporal_diags('RESULTS/bzE')
  bzE = (bzE[1:] + bzE[0:-1])*0.5

  # Opening of the exfield
  t,byE = read_picsar_temporal_diags('RESULTS/byE')
  byE = (byE[1:] + byE[0:-1])*0.5
  
  # Opening of the exfield
  t,bxE = read_picsar_temporal_diags('RESULTS/bxE')
  bxE = (bxE[1:] + bxE[0:-1])*0.5
  
  total_energy = sum(kinE,axis=0) + ezE + exE + eyE + bzE + bxE + byE

  min_totalE = min(total_energy)
  max_totalE = max(total_energy)
  diffrel = (max_totalE - min_totalE)/max_totalE
  print (' Relative error on the total energy:',diffrel)
  if ttest: assert diffrel < 1e-2
      
  # Plotting 
  if tshow:
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2)
    ax = plt.subplot(gs[:, :])
  
    ax.plot(t[1:],kinE[0,:],label='Electron kinetic energy')
    ax.plot(t[1:],kinE[1,:],label='Proton kinetic energy')
    ax.plot(t[1:],ezE,label='Ez energy')
    ax.plot(t[1:],eyE,label='Ey energy')  
    ax.plot(t[1:],exE,label='Ex energy')
    ax.plot(t[1:],bzE,label='Bz energy')
    ax.plot(t[1:],byE,label='By energy')  
    ax.plot(t[1:],bxE,label='Bx energy',color='orange')  
    ax.plot(t[1:],total_energy,label='tot',color='k',ls=':')
    ax.plot([Tplasma,Tplasma],[0.,max(total_energy)],ls='--',color='k')

    ax.set_xlabel('t (s)')
    ax.set_ylabel('Energy (J)')
  
    ax.set_ylim([0,max(total_energy)*1.2])

    plt.annotate('', xy=(0, kinE[0,0]*0.5), xycoords='data',xytext=(Tplasma, kinE[0,0]*0.5), textcoords='data',arrowprops={'arrowstyle': '<->'})
  
    plt.text(Tplasma*0.5, 0.3*max(total_energy),'Plasma period = %g s'%Tplasma, horizontalalignment='center',verticalalignment='center')
  
  # Theoretical oscilations
  ekinth = zeros(len(kinE[0,:]))
  A = 0.5 * max(kinE[0,:])
  ekinth = A + A*cos(2*pi*t[1:]/(Tplasma*0.5))  

  if tshow:  
    ax.plot(t[1:],ekinth,ls='--',label='Th. Elec. kin. energy')
    ax.legend(loc='upper center',ncol=4,borderaxespad=-2)
  

  # Test oscillations
  diffth = abs(kinE[0,:] - ekinth)
  print
  print (' Maximum difference between theory and simulation:',max(diffth))
  if ttest: assert (max(diffth) < 1E-1*max(kinE[0,:]))," The difference between simulation and theory is too high (> %g)"%(1E-1*max(kinE[0,:]))
  
  # Test divergence
  print (' ________________________________________ ')
  print (' Check DivE = rho/eps0')
  
  # Opening of divE-rho
  t,diverho = read_picsar_temporal_diags('RESULTS/divE-rho')
  t,rho = read_picsar_temporal_diags('RESULTS/rho')
  t,dive= read_picsar_temporal_diags('RESULTS/divE')
  
  print()
  print ('max(||diverho||(t)):',max(diverho))
  print ('max(||diverho||/||rho||):',max(diverho/rho))
  print ('max(||rho||(t)):',max(rho))
  print ('min(||rho||(t)):',min(rho)) 
  print ('max(||divE||(t)):',max(dive)*eps0)
  print ('min(||divE||(t)):',min(dive)*eps0) 
  print()

  #if ttest: assert (max(diverho/rho) < 1E-3),"L2 norm ||DivE*eps0 - rho||/||rho|| too high  (> 1E-3)"

  if tshow:
    fig1 = plt.figure(figsize=(12,8))
    gs1 = gridspec.GridSpec(7, 4)
    ax1 = plt.subplot(gs1[0:3, :])
    ax2 = plt.subplot(gs1[4:7, :])
  
    ax1.plot(t,diverho/rho,label=r'$||\nabla E \times \varepsilon_0 - \rho||$',lw=2) 
    ax1.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
    ax1.set_xlabel('t (s)')
    ax1.set_yscale('log')

    ax2.plot(t,dive*eps0,label=r'$|| \nabla E \times \varepsilon_0 ||$',lw=2) 
    ax2.plot(t,rho,label=r'$|| \rho ||$',color='r',lw=2,ls='--')
    ax2.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
    ax2.set_xlabel('t (s)')
    ax2.set_yscale('log')
  
  # Analyse of the files
  if 1: # Temporarily removed due to MPI-IO issues (plateform dependent)
      for it in range(0,30,10):
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
        
        #if ttest: assert (max_F < 1E-3),"L2 norm ||DivE*eps0 - rho|| too high  (> 1E-3)"
        

  # ____________________________________________________
  # Advice
  
  print ()
  print (' _______________________________________')
  print (' Advice for users:' )
  print (' - Check that the energy is constant with time')
  print (' - Check that divE = rho/eps0 for each tests')
  print (' - Check the energy oscillating behavior')
  print ()
  
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


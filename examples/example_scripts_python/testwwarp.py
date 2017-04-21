from mpi4py import MPI
import sys
import os
currentdir=os.getcwd()
sys.path.append(currentdir+'/python_bin/')
print(currentdir+'/python_bin/')
from picsar_python import picsarpy as pxrpy
import numpy as np
import warp as wp

pxr = pxrpy.picsar

#### Input parameters (replace input file)
pxr.default_init()
##Simulation box size and duration
pxr.nx_global_grid=10
pxr.ny_global_grid=10
pxr.nz_global_grid=10
pxr.xmin=0.0
pxr.ymin=0.0
pxr.zmin=0.0
pxr.tmax=0.5

## Resolution
pxr.dx=6.5e-7
pxr.dy=6.5e-7
pxr.dz=6.5e-7

## Maxwell solver
pxr.norderx=2
pxr.nordery=2
pxr.norderz=2


## Current deposition order
pxr.nox=1
pxr.noy=1
pxr.noz=1

## Tiling parameters
pxr.ntilex=5
pxr.ntiley=5
pxr.ntilez=5

## Particles species
pxr.nspecies=2 # Set number of particle species
pxr.pdistr=2 # Random distribution for particles
# Set particle species 1 properties
pxr.set_particle_species_properties(1,"electron",1.,-1.,10,0.,0.,0.,1.3,1.3,1.3,0.,0.,0.1,0.,0.,0.)
# Set particle species 2 properties
pxr.set_particle_species_properties(2,"proton",1836.,-1.,10,0.,0.,0.,1.3,1.3,1.3,0.,0.,0.,0.,0.,0.)



#### MPI INIT
pxr.mpi_minimal_init()
pxr.mpi_initialise()
mpirank= pxr.rank


#### Init particle distributions  for each species (with tiles)
pxr.initall() #Fortran routine
#initall_py() # Python init (To be written)


#### PIC LOOP with intrinsic step function
ntsteps=10
#start=MPI.Wtime()
#pxr.step(ntsteps)
#endt=MPI.Wtime()
#if (mpirank==0):
#    print("Total simulation time with step in Fortran 90 (s) ="+str(endt-start))

def steppy(nt):
    for i in range(0,nt):
        # Push particles
        pxr.push_particles()
        pxr.particle_bcs()
        
        # Deposit currents on the grid
        pxr.depose_currents_on_grid_jxjyjz()
        pxr.current_bcs()

        # Push EM fields
        pxr.push_bfield()
        pxr.bfield_bcs()
        
        pxr.push_efield()
        pxr.efield_bcs()
        pxr.push_bfield()
        pxr.bfield_bcs()

        # Compute diags
        pxr.calc_diags()

        # print diags
        if (mpirank==0):
            print("Iteration number "+str(i)+"|| time (s) ="+str(pxr.dt*i))


def pntfieldswarp():
    wp.w3d.nx = pxr.nx_global
    wp.w3d.ny = pxr.ny_global
    wp.w3d.nz = pxr.nz_global
    wp.w3d.xmmin = pxr.xmin
    wp.w3d.xmmax = pxr.xmax
    wp.w3d.ymmin = pxr.ymin
    wp.w3d.ymmax = pxr.ymax
    wp.w3d.zmmin = pxr.zmin
    wp.w3d.zmmax = pxr.zmax
    
    wp.em = wp.EM3D(nxguard=pxr.nxguards,
                    nyguard=pxr.nyguards,
                    nzguard=pxr.nzguards)

    wp.em.finalize()
    del wp.em.fields.Ex; wp.em.fields.Ex = pxr.ex
    del wp.em.fields.Ey; wp.em.fields.Ey = pxr.ey
    del wp.em.fields.Ez; wp.em.fields.Ez = pxr.ez
    del wp.em.fields.Bx; wp.em.fields.Ex = pxr.bx
    del wp.em.fields.By; wp.em.fields.Ey = pxr.by
    del wp.em.fields.Bz; wp.em.fields.Ez = pxr.bz
    
pntfieldswarp()

#### PIC LOOP written in python


#start=MPI.Wtime()
steppy(ntsteps)
#endt=MPI.Wtime()
#if (mpirank==0):
#    print("Total simulation time with step in python (s) ="+str(endt-start))

wp.em.pfezg(ncolor=21)

#wp.winon()


## MPICLOSE

#pxr.mpi_close()


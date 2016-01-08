from warp import *
from mpi4py import MPI
import sys
import os
currentdir=os.getcwd()
sys.path.append(currentdir+'/python_bin/')
print(currentdir+'/python_bin/')
import picsarpy as pxr
import numpy as np


#### Input parameters (replace input file)
pxr.default_init()
##Simulation box size and duration
pxr.picsar.nx_global_grid=100
pxr.picsar.ny_global_grid=100
pxr.picsar.nz_global_grid=100
pxr.picsar.xmin=0.0
pxr.picsar.ymin=0.0
pxr.picsar.zmin=0.0
pxr.picsar.tmax=0.5

## Resolution
pxr.picsar.dx=6.5e-7
pxr.picsar.dy=6.5e-7
pxr.picsar.dz=6.5e-7

## Maxwell solver
pxr.picsar.norderx=2
pxr.picsar.nordery=2
pxr.picsar.norderz=2


## Current deposition order
pxr.picsar.nox=1
pxr.picsar.noy=1
pxr.picsar.noz=1

## Tiling parameters
pxr.picsar.ntilex=5
pxr.picsar.ntiley=5
pxr.picsar.ntilez=5

#### Init particle distribution

## Properties of particle species
name=["electron","proton"]
charge=[-1.,1.]
mass=[1.,1836.]
nppcell=[10,10]
xmin=[0.,0.];ymin=[0.,0.];zmin=[0.,0.]
xmax=[1.3,1.3];ymax=[1.3,1.3];zmax=[1.3,1.3]
vdriftx=[0.,0.]; vdrifty=[0.,0.]; vdriftz=[0.,0.]
vthx=[0.,0.]; vthy=[0.,0.]; vthz=[0.,0.]

## Set number of particle species
pxr.picsar.nspecies=len(name)

## Set particle species  properties in PICSAR
for i in range(0,pxr.picsar.nspecies):
    pxr.set_particle_species_properties(i+1,name[i],mass[i],charge[i],nppcell[i],xmin[i],ymin[i],zmin[i],xmax[i],ymax[i],zmax[i],vdriftx[i],vdrifty[i],vdriftz[i],vthx[i],vthy[i],vthz[i])

#### MPI INIT
pxr.mpi_minimal_init()
pxr.mpi_initialise()
mpirank= pxr.picsar.rank


#### Routine that Init particle distributions  for each species (with tiles)
def density_profile_global(x,y,z):
    # Use defined plasma profile goes here
    # e.g: plasma disk of radius Rd
    r=np.sqrt(x**2+y**2+z**2)
    Rd=20*pxr.picsar.dx
    return np.where(r>Rd,0.,1.)

#### Routine that load particles in PICSAR with the density profile
#### specified in routine density_profile_global
def py_load_particles():
        for ispecies in range(0,pxr.picsar.nspecies):
                x,y,z=np.mgrid[0:pxr.picsar.nx,0:pxr.picsar.ny,0:pxr.picsar.nz]
                x=x.flatten()
                y=y.flatten()
                z=z.flatten()
                x=x+pxr.picsar.x_grid_min_local+(x-0.5)*pxr.picsar.dx;
                y=y+pxr.picsar.y_grid_min_local+(y-0.5)*pxr.picsar.dy;
                z=z+pxr.picsar.z_grid_min_local+(z-0.5)*pxr.picsar.dz;
                dens=density_profile_global(x,y,z)
                x=np.extract(dens>0.,x)
                y=np.extract(dens>0.,y)
                z=np.extract(dens>0.,z)
                dens=np.extract(dens>0.,dens)
                partw = dens*pxr.picsar.nc*pxr.picsar.dx*pxr.picsar.dy*pxr.picsar.dz/(nppcell[ispecies])
                for i in range(0,nppcell[ispecies]):
                    pxr.py_add_particles_to_species(ispecies+1,len(x),x+pxr.picsar.dx*i/nppcell[ispecies],\
                    y+pxr.picsar.dy*i/nppcell[ispecies], z+pxr.picsar.dz*i/nppcell[ispecies], 0., 0., 0., partw)


def initallpy():
    # Set time step
    pxr.picsar.dtcoef=0.7
    pxr.picsar.dt = pxr.picsar.dtcoef/(pxr.picsar.clight* \
    np.sqrt(1.0/pxr.picsar.dx**2+1.0/pxr.picsar.dy**2+1.0/pxr.picsar.dz**2))
    # Set tile split
    pxr.set_tile_split()
    # Allocate and init array of tiles
    pxr.init_tile_arrays()
    # Load particles
    py_load_particles()
    # Init stencil coefficients
    pxr.init_stencil_coefficients()

#### PIC LOOP with intrinsic step function
ntsteps=10
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
            print("Iteration number "+str(i)+"|| time (s) ="+str(pxr.picsar.dt*i))


#### INIT+ PIC LOOP written in python
initallpy()
start=MPI.Wtime()
steppy(ntsteps)
endt=MPI.Wtime()
if (mpirank==0):
    print("Total simulation time with step in python (s) ="+str(endt-start))


## MPICLOSE

pxr.mpi_close()


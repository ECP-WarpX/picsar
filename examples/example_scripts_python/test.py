from mpi4py import MPI
import sys
import os
from picsar_python import picsarpy as pxrpy
import numpy as np
home=os.getenv('HOME')
#currentdir=os.getcwd()
#sys.path.append(currentdir+'/python_bin/')
#print(currentdir+'/python_bin/')

pxr=pxrpy.picsar

#### Input parameters (replace input file)
pxr.default_init()

pxr.topology = 0

pxr.nprocx = 1
pxr.nprocy = 1
pxr.nprocz = 1

##Simulation box size and duration
pxr.nx_global_grid=100
pxr.ny_global_grid=100
pxr.nz_global_grid=100
pxr.xmin=-5e-6
pxr.xmax=5e-6
pxr.ymin=-5e-6
pxr.ymax=5e-6
pxr.zmin=-5e-6
pxr.zmax=5e-6
pxr.dtcoef=0.7

## Maxwell solver
pxr.norderx=2
pxr.nordery=2
pxr.norderz=2

## Current deposition order
pxr.nox=1
pxr.noy=1
pxr.noz=1

## Number of guard cells
pxr.nxguards=2
pxr.nyguards=2
pxr.nzguards=2
pxr.nxjguards=2
pxr.nyjguards=2
pxr.nzjguards=2

## Tiling parameters
pxr.ntilex=5
pxr.ntiley=5
pxr.ntilez=5

## Sorting
pxr.sorting_activated=1
pxr.sorting_dx = 1.
pxr.sorting_dy = 1.
pxr.sorting_dz = 1.
pxr.sorting_shiftx = 0.
pxr.sorting_shifty = 0.
pxr.sorting_shiftz = 0.
        
#### Init particle distribution

## Properties of particle species
name=["electron","proton"]
charge=[-1.,1.]
mass=[1.,1836.]
nppcell=[5,5]
xmin=[-1.3,-1.3];ymin=[-1.3,-1.3];zmin=[-1.3,-1.3]
xmax=[1.3,1.3];ymax=[1.3,1.3];zmax=[1.3,1.3]
vdriftx=[0.,0.]; vdrifty=[0.,0.]; vdriftz=[0.,0.]
vthx=[0.,0.]; vthy=[0.,0.]; vthz=[0.,0.]
sorting_period = [5,5]
sorting_start = [0,0]

## Set number of particle species
pxr.nspecies=len(name)

## Set particle species  properties in PICSAR
print " pxr.set_particle_species_properties"
for i in range(0,pxr.nspecies):
    pxr.set_particle_species_properties(i+1,name[i],mass[i],charge[i],nppcell[i],\
    xmin[i],ymin[i],zmin[i],xmax[i],ymax[i],zmax[i],vdriftx[i],vdrifty[i],vdriftz[i],\
    vthx[i],vthy[i],vthz[i],sorting_period[i],sorting_start[i])

#### MPI INIT
print " MPI init"
pxr.mpi_minimal_init(-1)
pxr.mpi_initialise()
mpirank= pxr.rank


#### Routine that Init particle distributions  for each species (with tiles)
def density_profile_global(x,y,z):
    # Use defined plasma profile goes here
    # e.g: plasma disk of radius Rd
    r=np.sqrt(x**2+y**2+z**2)
    Rd=(pxr.xmax-pxr.xmin)
    return np.where(r>Rd,0.,1.)

#### Routine that load particles in PICSAR with the density profile
#### specified in routine density_profile_global
def py_load_particles():
        for ispecies in range(0,pxr.nspecies):
                x,y,z=np.mgrid[0:pxr.nx,0:pxr.ny,0:pxr.nz]
                x=x.flatten()
                y=y.flatten()
                z=z.flatten()
                x=pxr.x_grid_min_local+x*pxr.dx
                y=pxr.y_grid_min_local+y*pxr.dy
                z=pxr.z_grid_min_local+z*pxr.dz
                dens=density_profile_global(x,y,z)
                x=np.extract(dens>0.,x)
                y=np.extract(dens>0.,y)
                z=np.extract(dens>0.,z)
                dens=np.extract(dens>0.,dens)
                partw = dens*pxr.nc*pxr.dx*pxr.dy*pxr.dz/(nppcell[ispecies])
                for i in range(0,nppcell[ispecies]):
                    pxr.py_add_particles_to_species(ispecies+1,len(x),x+pxr.dx*i/nppcell[ispecies],\
                    y+pxr.dy*i/nppcell[ispecies], z+pxr.dz*i/nppcell[ispecies], 
                    np.zeros(np.size(x)), np.zeros(np.size(x)), np.zeros(np.size(x)),np.ones(np.size(x)), partw)

def initallpy():
    pxr.dt = pxr.dtcoef/(pxr.clight* \
    np.sqrt(1.0/pxr.dx**2+1.0/pxr.dy**2+1.0/pxr.dz**2))

    pxr.sorting_dx *= pxr.dx
    pxr.sorting_dy *= pxr.dy
    pxr.sorting_dz *= pxr.dz

    # Set tile split
    print " pxr.set_tile_split()"
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
        
        print(" Starting iteration "+str(i))
    
        # Push particles
        pxr.push_particles()
        pxr.particle_bcs()
        
        # Sorting
        pxr.particle_sorting_sub()
        
        # Deposit currents on the grid
        pxr.pxrdepose_currents_on_grid_jxjyjz()
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


#### INIT+ PIC LOOP written in python
print " initallpy()"
initallpy()
start=MPI.Wtime()
steppy(ntsteps)
endt=MPI.Wtime()
if (mpirank==0):
    print("Total simulation time with step in python (s) ="+str(endt-start))


## MPICLOSE

pxr.mpi_close()


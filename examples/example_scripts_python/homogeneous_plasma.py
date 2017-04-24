from warp import *
from warp.field_solvers.em3dsolverPXR import *
import os
from warp.data_dumping.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from mpi4py import MPI
EnableAll()
home=os.getenv('HOME')

l_pxr=1

# --- flags turning off unnecessary diagnostics (ignore for now)
top.ifzmmnt = 0
top.itmomnts = 0
top.itplps = 0
top.itplfreq = 0
top.zzmomnts = 0
top.zzplps = 0
top.zzplfreq = 0
top.nhist = top.nt
top.iflabwn = 0
w3d.lrhodia3d = false
w3d.lgetese3d = false
w3d.lgtlchg3d = false
#top.nxproc = 2
#top.nyproc = 2
#top.nzproc = 2

# ----------
# Parameters
# ----------

dfact = 1
dxfact = 16
dtfact = 8*4
N_step = 20000/dtfact

#Two-layer foil:
#Carbon layer
carbon_layer_start     = 2
carbon_layer_width     = 10
carbon_layer_thickness = 0.075
carbon_layer_e_density = 400.
nppcell_carbon         = 25
#Hydrogen layer
hydrogen_layer_width     = 6
hydrogen_layer_thickness = 0.05
hydrogen_layer_e_density = 21.
nppcell_hydrogen         = 16

#Laser at the left border:
a0             = 10.
laser_duration = 10
laser_width    = 4/2 #(=2w_0)
#S-pol
#Transverse profile - Gaussian
#Longitudinal profile Super-Gaussian 12
#focused at x=6 from the left border

#Mesh: 
dt=0.005
dx=dy=dz=0.01

# --- scaling
carbon_layer_e_density/=dfact
hydrogen_layer_e_density/=dfact
dx*=dxfact
dy=dz=dx
dt*=dtfact

#-------------------------------------------------------------------------------
# main parameters
#-------------------------------------------------------------------------------
dim = "3d"                 # 3D calculation
#dim = "2d"                 # 2D calculation 
#dim = "1d"                 # 1D calculation 
dpi=100                     # graphics resolution
l_test             = 0      # Will open output window on screen 
                            # and stop before entering main loop.
l_gist             = 0      # Turns gist plotting on/off
l_restart          = false  # To restart simulation from an old run (works?)
restart_dump       = ""     # dump file to restart from (works?)
l_moving_window    = 0      # on/off (Galilean) moving window
l_plasma           = 1      # on/off plasma
l_usesavedist      = 0      # if on, uses dump of beam particles distribution
l_smooth           = 1      # on/off smoothing of current density
l_laser            = 0      # on/off laser
l_pdump            = 0      # on/off regular dump of beam data
stencil            = 0      # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F 
                            # use 0 or 1; 2 does not verify Gauss Law
if dim=="1d":stencil=0
dtcoef             = dt/dx  # coefficient to multiply default time step that is set at the EM solver CFL
top.depos_order    = 1      # particles deposition order (1=linear, 2=quadratic, 3=cubic)
top.efetch         = 4      # field gather type (1=from nodes "momentum conserving"; 4=from Yee mesh "energy conserving")

top.runid          = "homogeneous_plasma"                         # run name
top.pline1         = "basic lpa"                         # comment line on plots
top.runmaker       = "J.-L. Vay, M. Lobet"                        # run makers
top.lrelativ       = true                                # on/off relativity (for particles push)
top.pgroup.lebcancel_pusher=0                         # flag for particle pusher (0=Boris pusher; 1=Vay PoP 08 pusher)
#top.ibpush=2
l_verbose          = 0                                   # verbosity level (0=off; 1=on)

#-------------------------------------------------------------------------------
# diagnostics parameters + a few other settings
#-------------------------------------------------------------------------------
live_plot_freq     = 10   # frequency (in time steps) of live plots (off is l_test is off)

fielddiag_period   = 500/dtfact
partdiag_period    = 500/dtfact

#-------------------------------------------------------------------------------
# laser parameters
#-------------------------------------------------------------------------------
# --- in lab frame
lambda_laser       = 0.8e-6                      # wavelength 
laser_width       *= lambda_laser
laser_waist        = laser_width/2  
laser_radius       = laser_waist/2.354820        # FWHM -> radius
laser_duration     = laser_duration*lambda_laser/clight # laser duration 
laser_polangle     = 0#pi/2                       # polarization (0=aligned with x; pi/2=aligned with y)
k0                 = 2.*pi/lambda_laser
w0                 = k0*clight
ZR                 = 0.5*k0*(laser_waist**2)   # Rayleigh length
Eamp               = a0*w0*emass*clight/echarge
Bamp               = Eamp/clight 
if l_laser==0:Eamp=Bamp=0.

#-------------------------------------------------------------------------------
# plasma layers
#-------------------------------------------------------------------------------
dfact             = 1.                                  # coefficient factor for plasma density (for scaled simulations)
densc             = emass*eps0*w0**2/echarge**2         # critical density

#-------------------------------------------------------------------------------
# Carbon plasma 
#-------------------------------------------------------------------------------
dens0_C           = dfact*carbon_layer_e_density*densc    # plasma density
wp_C              = sqrt(dens0_C*echarge**2/(eps0*emass)) # plasma frequency
kp_C              = wp_C/clight                           # plasma wavenumber
lambda_plasma_C   = 2.*pi/kp_C                            # plasma wavelength

#-------------------------------------------------------------------------------
# Carbon plasma 
#-------------------------------------------------------------------------------
dens0_H           = dfact*hydrogen_layer_e_density*densc  # plasma density
wp_H              = sqrt(dens0_H*echarge**2/(eps0*emass)) # plasma frequency
kp_H              = wp_H/clight                            # plasma wavenumber
lambda_plasma_H   = 2.*pi/kp_H                            # plasma wavelength

#-------------------------------------------------------------------------------
# print some plasma parameters to the screen
#-------------------------------------------------------------------------------
print("the laser spot size is: ")
print laser_waist
print("the Rayleigh length is: ")
print ZR
print("the laser wavelength is: ")
print lambda_laser
print("the Carbon plasma wavelength is: ")
print lambda_plasma_C
print("the Hydrogen plasma wavelength is: ")
print lambda_plasma_H

#-------------------------------------------------------------------------------
# number of plasma macro-particles/cell
#-------------------------------------------------------------------------------
nppcellx_C = 1#5
nppcelly_C = 1#5
nppcellz_C = 1#5

nppcellx_H = 1#4
nppcelly_H = 1#4
nppcellz_H = 1#4

if dim=="2d":
  nppcelly_C = nppcelly_H = 1
if dim=="1d":
  nppcellx_C = nppcellx_H = 1
  nppcelly_C = nppcelly_H = 1

#-------------------------------------------------------------------------------
# grid dimensions, nb cells and BC
#-------------------------------------------------------------------------------
w3d.zmmax = 5*lambda_laser
w3d.zmmin = -5.*lambda_laser
w3d.xmmin = -5*lambda_laser
w3d.xmmax = 5*lambda_laser
w3d.ymmin = -5*lambda_laser
w3d.ymmax = 5*lambda_laser

w3d.nx = nint((w3d.xmmax-w3d.xmmin)/(dx*lambda_laser))
w3d.ny = nint((w3d.ymmax-w3d.ymmin)/(dy*lambda_laser))
w3d.nz = nint((w3d.zmmax-w3d.zmmin)/(dz*lambda_laser))

if dim in ["1d"]:
    w3d.nx = 2
    w3d.xmmin = -float(w3d.nx)/2
    w3d.xmmax = float(w3d.nx)/2
if dim in ["1d","2d"]:
    w3d.ny = 2
    w3d.ymmin = -float(w3d.ny)/2
    w3d.ymmax = float(w3d.ny)/2

w3d.dx = (w3d.xmmax-w3d.xmmin)/w3d.nx
w3d.dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
w3d.dz = (w3d.zmmax-w3d.zmmin)/w3d.nz

# --- sets field boundary conditions
# --- longitudinal
w3d.bound0  = w3d.boundnz = openbc
# --- transverse
w3d.boundxy = periodic

# --- sets particles boundary conditions
# --- longitudinal
top.pbound0  = absorb
top.pboundnz = absorb
# --- transverse
top.pboundxy = periodic

#-------------------------------------------------------------------------------
# set graphics
#-------------------------------------------------------------------------------
if l_gist:
 if l_test:
  winon(0,dpi=dpi)
 else:
  setup()
else:
 setup()
 
#-------------------------------------------------------------------------------
# set particles weights
#-------------------------------------------------------------------------------
weight_C   = dens0_C*w3d.dx*w3d.dy*w3d.dz/(nppcellx_C*nppcelly_C*nppcellz_C) 
weight_H   = dens0_H*w3d.dx*w3d.dy*w3d.dz/(nppcellx_H*nppcelly_H*nppcellz_H) 
top.wpid = nextpid() # Activate variable weights in the method addpart

# --- create plasma species
elec_C = Species(type=Electron,weight=weight_C,name='elec_C')
elec_H = Species(type=Electron,weight=weight_H,name='elec_H')
ions_C = Species(type=Carbon,weight=weight_C/6,charge_state=6.,name='ions_C')
ions_H = Species(type=Proton,weight=weight_H,name='ions_H')

top.depos_order[...] = top.depos_order[0,0] # sets deposition order of all species = those of species 0
top.efetch[...] = top.efetch[0] # same for field gathering
if dim in ["1d","2d"]:
  top.depos_order[1,:]=1
if dim=="1d":
  top.depos_order[0,:]=1

#-------------------------------------------------------------------------------
# set smoothing of current density
#-------------------------------------------------------------------------------
if l_smooth:
  # --- 1 time nilinear (0.25,0.5,0.25) + 1 time relocalization (-1, 3/2,-1.)
  npass_smooth = [[ 1 , 1 ],[ 0 , 0 ],[ 1 , 1 ]]
  alpha_smooth = [[ 0.5, 3./2],[ 0.5, 3.],[0.5, 3./2]]
  stride_smooth = [[ 1 , 1 ],[ 1 , 1 ],[ 1 , 1 ]]
  if dim=='1d':
    for i in range(len(npass_smooth[0])):
      npass_smooth[0][i]=0
  if dim in ['1d','2d']:
    for i in range(len(npass_smooth[0])):
      npass_smooth[1][i]=0
else:
  npass_smooth = [[ 0 ],[ 0 ],[ 0 ]]
  alpha_smooth = [[ 1.],[ 1.],[ 1.]]
  stride_smooth = [[ 1 ],[ 1 ],[ 1 ]]

#-------------------------------------------------------------------------------
# initializes WARP
#-------------------------------------------------------------------------------
top.fstype = -1 # sets field solver to None (desactivates electrostatic solver)
package('w3d');
generate()
#-------------------------------------------------------------------------------
# set a few shortcuts
#-------------------------------------------------------------------------------
pg = top.pgroup

if l_plasma:

    zmin = w3d.zmmin
    zmax = w3d.zmmax
    xmin = w3d.xmmin
    xmax = w3d.xmmax
    if dim=='3d':
        ymin= w3d.ymmin
        ymax= w3d.ymmax
        np = w3d.nx*w3d.ny*nint((zmax-zmin)/w3d.dz)*nppcellx_H*nppcelly_H*nppcellz_H
    else:
        ymin=ymax=0.
        np = w3d.nx*nint((zmax-zmin)/w3d.dz)*nppcellx_H*nppcelly_H*nppcellz_H
    
    elec_H.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,
                       vthx=0.2*clight,vthy=0.,vthz=0.,
                       spacing='uniform')

    ions_H.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,
                       vthx=0.,vthy=0.,vthz=0.,
                       spacing='uniform')


laser_total_duration=1.25*laser_duration
#-------------------------------------------------------------------------------
# set laser pulse shape
#-------------------------------------------------------------------------------
def laser_amplitude(time):
 global laser_total_duration,Eamp
 #fx=Exp[-((x-t)*2/xblen)]^12
 #xblen=10 - duration of the pulse
 return Eamp*exp(-(2*(time-0.5*laser_total_duration)/laser_duration)**12)

def laser_profile(x,y,z):
  global laser_waist, ZR
  #fy=Exp[-(y*2/yblen)^2]
  #yblen=4
  r2 = x**2 + y**2
  W0 = laser_waist
  Wz = W0*sqrt(1.+(z/ZR)**2) # Beam radius varies along propagation direction
  return (W0/Wz)*exp(-r2/Wz**2)

#-------------------------------------------------------------------------------
# set laser amplitude by combining the pulse shape, laser profile, and laser phase
#-------------------------------------------------------------------------------

if l_laser:
  laser_source_z=10*w3d.dz
  
  #-------------------------------------------------------------------------------
  def laser_func(x,y,t):
    global laser_amplitude,laser_phase,laser_profile,k0,w0,ZR
    z0 = (carbon_layer_start*lambda_laser-laser_source_z)
    Rz = z0*(1.+(ZR/z0)**2)
    E = laser_amplitude(t)*laser_profile(x,y,z0)
    r = sqrt(x*x+y*y)
    angle = w0*t+k0*z0-arctan(z0/ZR)+k0*r**2/(2.*Rz)   
    #return [E*cos(angle),E*sin(angle)]   # for circularly polarized laser
    return [0,E*sin(angle)]               # for linearly polarized laser
    #  return [E*sin(angle),0]  
  
  
else:
  laser_func=laser_source_z=Eamp=laser_polangle=None
             # for linearly polarized laser

#-------------------------------------------------------------------------------
# initializes main field solver block
#-------------------------------------------------------------------------------
if l_pxr:
    ntilex = 1#max(1,w3d.nx/30)
    ntiley = 1#max(1,w3d.ny/30)
    ntilez = 1#max(1,w3d.nz/30)
#    pg.sw=0.
    em = EM3DPXR(laser_func=laser_func,
                 laser_source_z=laser_source_z,
                 laser_polangle=laser_polangle,
                 laser_emax=Eamp,
                 stencil=stencil,
                 npass_smooth=npass_smooth,
                 alpha_smooth=alpha_smooth,
                 stride_smooth=stride_smooth,
                 l_2dxz=dim=="2d",
                 l_1dz=dim=="1d",
                 dtcoef=dtcoef,
                 l_getrho=1,
                 spectral=0,
                 current_cor=0,
                 listofallspecies=listofallspecies,
                 ntilex=ntilex,
                 ntiley=ntiley,
                 ntilez=ntilez,
                 # Guard cells
                 #nxguard=6,
                 #nyguard=6,
                 #nzguard=6,
                 currdepo=currdepo,   
                 mpicom_curr=mpicom_curr,
                 fieldgathe=fieldgathe,
                 sorting=sort,
                 partcom=partcom,
                 l_verbose=l_verbose)
    step = em.step
else:
    em = EM3D(   laser_func=laser_func,
                 laser_source_z=laser_source_z,
                 laser_polangle=laser_polangle,
                 laser_emax=Eamp,
                 stencil=stencil,
                 npass_smooth=npass_smooth,
                 alpha_smooth=alpha_smooth,
                 stride_smooth=stride_smooth,
                 l_2dxz=dim=="2d",
                 l_1dz=dim=="1d",
                 dtcoef=dtcoef,
                 l_getrho=1,
                 l_verbose=l_verbose)

#-------------------------------------------------------------------------------
# restarts from dump file
#-------------------------------------------------------------------------------
if l_restart:
  restore(dump_file)

#-------------------------------------------------------------------------------
# register solver
#-------------------------------------------------------------------------------
print 'register solver'
registersolver(em)
em.finalize()
loadrho()

print 'done'

# _______________________________________________
# Install After Step:
  
def liveplots():
  print "liveplots"

installafterstep(liveplots)


# Load additional OpenPMD diagnostic
diag_f = FieldDiagnostic( period=fielddiag_period, top=top, w3d=w3d, em=em,
                          comm_world=comm_world, lparallel_output=lparallel )
diag_elec_H = ParticleDiagnostic( period=partdiag_period, top=top, w3d=w3d,
            species = {"elec_H" : elec_H},
            comm_world=comm_world, lparallel_output=lparallel )
diag_ions_H = ParticleDiagnostic( period=partdiag_period, top=top, w3d=w3d,
            species = {"ions_H" : ions_H},
            comm_world=comm_world, lparallel_output=lparallel )
#installafterstep( diag_f.write )
#installafterstep( diag_elec_C.write )
#installafterstep( diag_elec_H.write )
#installafterstep( diag_ions_C.write )
#installafterstep( diag_ions_H.write )

# Time 
tottime = AppendableArray()
def accuttime():
  global tottime
  tottime.append(time.clock())
  if me==0 and top.it%200==0:
    f=PW.PW('tottime.pdb')
    f.time=tottime[:]
    f.close()

installafterstep(accuttime)

if me==0:
    f=PW.PW('sim_params.pck')
    f.weight_H=weight_H
    f.dt=top.dt
    f.dens0_H=dens0_H
    f.close()

print '\nInitialization complete\n'

# if this is a test, then stop, else execute main loop
if l_test:
  print '<<< To execute n steps, type "step(n)" at the prompt >>>'
#  raise('')
else:
  step(N_step,1,1)
  


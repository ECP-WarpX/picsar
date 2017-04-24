from warp import *
from warp.field_solvers.em3dsolverPXR import *
import os
from warp.data_dumping.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from mpi4py import MPI
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
top.lcomm_cartesian=1
w3d.lrhodia3d = false
w3d.lgetese3d = false
w3d.lgtlchg3d = false

# Flags turning off auto decomp and let user specify its decomp
top.lautodecomp = 1 # Particles
top.lfsautodecomp = 1 # fields 

# Flags turning on/off load balancing
load_balance=1
dlb_freq=11
dlb_threshold=5 # dynamic load balancing threshold in % 
dlb_at_init=0 # Do a load balancing of the simulation at init 

# ----------
# Parameters
# ----------

dfact = 1
dxfact = 8
dtfact = 8
N_step = 20000/dtfact

#Two-layer foil:
#Carbon layer
carbon_layer_start     = 2
carbon_layer_width     = 6
carbon_layer_thickness = 0.075
carbon_layer_e_density = 4.
nppcell_carbon         = 2500
#Hydrogen layer
hydrogen_layer_width     = 6
hydrogen_layer_thickness = 0.05
hydrogen_layer_e_density = 2.
nppcell_hydrogen         = 1600

#Laser at the left border:
a0             = 100.
laser_duration = 10
laser_width    = 4/2 #(=2w_0)
#S-pol
#Transverse profile - Gaussian
#Longitudinal profile Super-Gaussian 12
#focused at x=6 from the left border

#Mesh: 
dt=0.0015
dx=dy=0.0035
dz=0.002

# --- scaling
carbon_layer_e_density/=dfact
hydrogen_layer_e_density/=dfact
dx*=dxfact
dz=dy=dx
dt*=dtfact

#-------------------------------------------------------------------------------
# main parameters
#-------------------------------------------------------------------------------
dim = "3d"                 # 3D calculation
#dim = "2d"                 # 2D calculation 
#dim = "1d"                 # 1D calculation 
dpi=100                     # graphics resolution
l_test             = 1     # Will open output window on screen
                            # and stop before entering main loop.
l_gist             = 1      # Turns gist plotting on/off
l_restart          = false  # To restart simulation from an old run (works?)
restart_dump       = ""     # dump file to restart from (works?)
l_moving_window    = 1      # on/off (Galilean) moving window
l_plasma           = 1    # on/off plasma
l_usesavedist      = 0      # if on, uses dump of beam particles distribution
l_smooth           = 1      # on/off smoothing of current density
l_laser            = 1      # on/off laser
l_pdump            = 0      # on/off regular dump of beam data
stencil            = 0      # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F 
                            # use 0 or 1; 2 does not verify Gauss Law
if dim=="1d":stencil=0
dtcoef             = dt/dx  # coefficient to multiply default time step that is set at the EM solver CFL
top.depos_order    = 1      # particles deposition order (1=linear, 2=quadratic, 3=cubic)
top.efetch         = 4      # field gather type (1=from nodes "momentum conserving"; 4=from Yee mesh "energy conserving")

top.runid          = "ion acceleration"                         # run name
top.pline1         = "basic lpa"                         # comment line on plots
top.runmaker       = "H. Vincenti,"                        # run makers
top.lrelativ       = true                                # on/off relativity (for particles push)
top.pgroup.lebcancel_pusher=0                         # flag for particle pusher (0=Boris pusher; 1=Vay PoP 08 pusher)
#top.ibpush=2
l_verbose          = 0                                   # verbosity level (0=off; 1=on)

#-------------------------------------------------------------------------------
# diagnostics parameters + a few other settings
#-------------------------------------------------------------------------------
live_plot_freq     = 1000  # frequency (in time steps) of live plots (off is l_test is off)

fielddiag_period   = 150#200000/dtfact
partdiag_period    = 150#200000/dtfact
l_parallelo=False

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
nppcellx_C = 2#5
nppcelly_C = 2#5
nppcellz_C = 2#5

nppcellx_H = 2#4
nppcelly_H = 2#4
nppcellz_H = 2#4

if dim=="2d":
  nppcelly_C = nppcelly_H = 1
if dim=="1d":
  nppcellx_C = nppcellx_H = 1
  nppcelly_C = nppcelly_H = 1

#-------------------------------------------------------------------------------
# grid dimensions, nb cells and BC
#-------------------------------------------------------------------------------
w3d.zmmax = carbon_layer_start*2.*lambda_laser
w3d.zmmin = 0.
w3d.xmmin = -0.5*carbon_layer_width*lambda_laser
w3d.xmmax = -w3d.xmmin
w3d.ymmin = -0.5*carbon_layer_width*lambda_laser
w3d.ymmax = -w3d.ymmin

w3d.nx = nint((w3d.xmmax-w3d.xmmin)/(dx*lambda_laser))
w3d.ny = nint((w3d.ymmax-w3d.ymmin)/(dy*lambda_laser))
w3d.nz = nint((w3d.zmmax-w3d.zmmin)/(dz*lambda_laser))

if dim in ["1d"]:
    w3d.nx = 2
    w3d.xmmin = -float(w3d.nx)/2.
    w3d.xmmax = float(w3d.nx)/2.
if dim in ["1d","2d"]:
    w3d.ny = 2
    w3d.ymmin = -float(w3d.ny)/2.
    w3d.ymmax = float(w3d.ny)/2.

w3d.dx = (w3d.xmmax-w3d.xmmin)/w3d.nx
w3d.dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
w3d.dz = (w3d.zmmax-w3d.zmmin)/w3d.nz

# --- sets field boundary conditions
# --- longitudinal
w3d.bound0  = w3d.boundnz = periodic
# --- transverse
w3d.boundxy = periodic

# --- sets particles boundary conditions
# --- longitudinal
top.pbound0  = periodic
top.pboundnz = periodic
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
weight_H   = dens0_H*w3d.dx*w3d.dy*w3d.dz/(nppcellx_C*nppcelly_C*nppcellz_C)
top.wpid = nextpid() # Activate variable weights in the method addpart

# --- create plasma species
elec_C = Species(type=Electron,weight=weight_C,name='elec_C')
elec_H = Species(type=Electron,weight=weight_H,name='elec_H')
ions_C = Species(type=Carbon,weight=weight_C/6.,charge_state=6.,name='ions_C')
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
# User defined decomp goes here 
#top.nxprocs=2
#top.nyprocs=2
#top.nzprocs=2
#top.fsdecomp.nxprocs=top.nxprocs
#top.fsdecomp.nyprocs=top.nyprocs
#top.fsdecomp.nzprocs=top.nzprocs
#top.ppdecomp.nxprocs=top.nxprocs
#top.ppdecomp.nyprocs=top.nyprocs
#top.ppdecomp.nzprocs=top.nzprocs
#top.fsdecomp.nx=[20,87]
#top.fsdecomp.nx=[107]
#top.fsdecomp.ny=[107]
#top.fsdecomp.ny=[20,87]
#top.fsdecomp.nz=[10,61]
#top.fsdecomp.nz=[71]
#top.ppdecomp.nx=[20,87]
#top.ppdecomp.nx=[107]
#top.ppdecomp.ny=[20,87]
#top.ppdecomp.ny=[107]
#top.ppdecomp.nz=[10,61]
#top.ppdecomp.nz=[71]
#top.userdecompx = top.fsdecomp.nx
#top.userdecompy = top.fsdecomp.ny
#top.userdecompz = top.fsdecomp.nz
    
top.fstype = -1 # sets field solver to None (desactivates electrostatic solver)
package('w3d'); generate()

#-------------------------------------------------------------------------------
# set a few shortcuts
#-------------------------------------------------------------------------------
pg = top.pgroup

if l_plasma:

    zmin = carbon_layer_start*lambda_laser
    zmax = zmin+carbon_layer_thickness*lambda_laser
    xmin = w3d.xmmin
    xmax = w3d.xmmax
    if dim=='3d':
        ymin=xmin
        ymax=xmax
        np = w3d.nx*w3d.ny*nint((zmax-zmin)/w3d.dz)*nppcellx_C*nppcelly_C*nppcellz_C
    else:
        ymin=ymax=0.
        np = w3d.nx*nint((zmax-zmin)/w3d.dz)*nppcellx_C*nppcelly_C*nppcellz_C
    
    elec_C.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,
                       vthx=0.,vthy=0.,vthz=0.,
                       spacing='uniform')

    ions_C.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,
                       vthx=0.,vthy=0.,vthz=0.,
                       spacing='uniform')
 
    zmin=zmax+0.
    zmax+=hydrogen_layer_thickness*lambda_laser                    
    if dim=='3d':
        ymin=xmin
        ymax=xmax
        np = w3d.nx*w3d.ny*nint((zmax-zmin)/w3d.dz)*nppcellx_H*nppcelly_H*nppcellz_H
    else:
        ymin=ymax=0.
        np = w3d.nx*nint((zmax-zmin)/w3d.dz)*nppcellx_H*nppcelly_H*nppcellz_H

    elec_H.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,
                       vthx=0.,vthy=0.,vthz=0.,
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
else:
  laser_func=laser_source_z=Eamp=None

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
#  return [E*sin(angle),0]               # for linearly polarized laser

#-------------------------------------------------------------------------------
# initializes main field solver block
#-------------------------------------------------------------------------------
if l_pxr:
    ntilex = max(1,w3d.nxlocal/10)
    ntiley = max(1,w3d.nylocal/10)
    ntilez = max(1,w3d.nzlocal/10)
#    pg.sw=0.
    em = EM3DPXR(       laser_func=laser_func,
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
                 l_getrho=0,
                 spectral=0,
                 current_cor=0,
                 listofallspecies=listofallspecies,
                 ntilex=ntilex,
                 ntiley=ntiley,
                 ntilez=ntilez,
                 l_verbose=l_verbose,
                 dload_balancing=load_balance, 
                 dlb_freq=dlb_freq, 
                 dlb_threshold=dlb_threshold, 
                 dlb_at_init=dlb_at_init)
    step = em.step
else:
    em = EM3D(       laser_func=laser_func,
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
                 l_getrho=0,
                 l_verbose=l_verbose)

print(em.nxlocal)
print(em.nylocal)
print(em.nzlocal)

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
  
def liveplots():
  if top.it%live_plot_freq==0:
    fma(0);


    if dim=="1d":
      em.pfex(view=3,titles=0,gridscale=1.e-12,direction=1)
      ptitles('Ex [V/m]','z [um]')
      density=elec_C.get_density()+elec_H.get_density()
      plsys(4)
      pla(density)
      ptitles('n','z [um]','X [um]')
      pzxez(view=5,msize=2,titles=0,gridscale=1.e-9)
      ptitles('Ez [V/m]','z [um]')
      refresh()
    else:    
      em.pfey(view=3,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-12,l_transpose=1,direction=1)
      ptitles('Ey [TV/m]','z [um]','X [um]')
      density=elec_C.get_density()+elec_H.get_density()
      #density=0
      if dim=='3d':density=density[:,w3d.ny/2,:]
      ppg(transpose(density),view=4,titles=0, \
          xmin=w3d.zmmin+top.zgrid,xmax=w3d.zmmax+top.zgrid, \
          ymin=w3d.xmmin,ymax=w3d.xmmax,
          xscale=1e6,yscale=1.e6)#,gridscale=1./dens0)

      ptitles('n','z [um]','X [um]')
      em.pfez(view=5,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-9,l_transpose=1,direction=1)
      ptitles('Ez [GV/m]','z [um]','X [um]')
      #ions_C.ppzx(view=6)
      #elec_C.ppzx(color=red,view=6)
      #ions_H.ppzx(color=blue,view=6)
      #elec_H.ppzx(color=cyan,view=6)

installafterstep(liveplots)


# Load additional OpenPMD diagnostic
diag_f = FieldDiagnostic( period=fielddiag_period, top=top, w3d=w3d, em=em,
                          comm_world=comm_world, lparallel_output=l_parallelo )
diag_elec_C = ParticleDiagnostic( period=partdiag_period, top=top, w3d=w3d,
            species = {"elec_C" : elec_C},
            comm_world=comm_world, lparallel_output=l_parallelo )
diag_elec_H = ParticleDiagnostic( period=partdiag_period, top=top, w3d=w3d,
            species = {"elec_H" : elec_H},
            comm_world=comm_world, lparallel_output=l_parallelo )
diag_ions_C = ParticleDiagnostic( period=partdiag_period, top=top, w3d=w3d,
            species = {"ions_C" : ions_C},
            comm_world=comm_world, lparallel_output=l_parallelo )
diag_ions_H = ParticleDiagnostic( period=partdiag_period, top=top, w3d=w3d,
            species = {"ions_H" : ions_H},
            comm_world=comm_world, lparallel_output=l_parallelo )

installafterstep( diag_f.write )
installafterstep( diag_elec_C.write )
installafterstep( diag_elec_H.write )
installafterstep( diag_ions_C.write )
installafterstep( diag_ions_H.write )

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
    f.weight_C=weight_C
    f.weight_H=weight_H
    f.dt=top.dt
    f.dens0_C=dens0_C
    f.dens0_H=dens0_H
    f.close()

print '\nInitialization complete\n'

# if this is a test, then stop, else execute main loop
if l_test:
  print '<<< To execute n steps, type "step(n)" at the prompt >>>'
  tdeb=MPI.Wtime()
  em.step(50,1,1)
  tend=MPI.Wtime()
  print("Final runtime (s): "+str(tend-tdeb))
#  raise('')
else:
  em.step(1000,1,1)
  
#pxr.point_to_tile(1,1,1,1)

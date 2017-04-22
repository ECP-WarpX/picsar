"""
Launches laser at z=-2.5 microns.
"""
from mpi4py import MPI
from warp import *
from warp.field_solvers.em3dsolverPXR import *
EnableAll()

l_2d    = 0

l_test  = 1 # --- open window on screen if true, save on disk in cgm file otherwise
ncells = 16*2 # --- nb cells in x
nzfact = 10 # --- multiplication factor for nb cells in z

# --- grid dimensions and nb cells
w3d.xmmin = -5.e-6     # --- min x of simulation box
w3d.xmmax = 5.e-6      # --- max x of simulation box
w3d.zmmin = -5.e-6     # --- min z of simulation box
w3d.zmmax = 5.e-6      # --- max z of simulation box
w3d.nx = ncells        # --- nb cells in x
w3d.nz = ncells*nzfact # --- nb cells in z
if l_2d:
  # --- sets y min/max and ny so that dy=1 in 2D
  w3d.ymmin = -1.
  w3d.ymmax = 1.
  w3d.ny = 2
else:
  # --- sets y min/max and ny to x values in 3D
  w3d.ymmin = w3d.xmmin
  w3d.ymmax = w3d.xmmax
  w3d.ny = ncells
w3d.dx = dx = (w3d.xmmax-w3d.xmmin)/w3d.nx
w3d.dy = dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
w3d.dz = dz = (w3d.zmmax-w3d.zmmin)/w3d.nz

# --- field boundary conditions
bounds = zeros(6)
bounds[:] = openbc
bounds[:] = periodic
  
# --- laser shape
lambda_laser       = 1.e-6
k0                 = 2.*pi/lambda_laser
w0                 = k0*clight
a0                 = 3.                         # normalized potential vector (amplitude)
laser_spotsize     = 0.2*(w3d.xmmax-w3d.xmmin)  # spot size
laser_risetime     = 25./w0
laser_polangle     = pi/2      # polarization (0=aligned with x; pi/2=aligned with y)
W0                 = laser_spotsize/2 # Radis of beam at waist
zr                 = (pi*W0**2)/lambda_laser # Rayleigh length  
E0 = a0*w0*emass*clight/echarge

# --- open graphics window (on screen if l_test=1, on disk otherwise)
if l_test:
  window(0,dpi=100)
else:
  setup()
  
palette('rainbow.gp')  
palette('bluewhitered.gp')

ions = Species(type=Proton) # --- unused but necessary to declare at least one particle species
top.fstype=-1               # --- turns off electrostatic solver
generate()                  # initializes internal arrays

#-------------------------------------------------------------------------------
# --- set laser pulse shape
#-------------------------------------------------------------------------------
# --- set position of laser injection
Zc = w3d.zmmin+0.25*(w3d.zmmax-w3d.zmmin)
# --- set position of laser focus (relative to Zc)
Zf = 1.5*zr

def laser_amplitude(time):
 return E0*exp(-0.5*((time-4*laser_risetime)/laser_risetime)**2)

def laser_amplitude(time):
 return E0 # assume no risetime for this example


def laser_profile(x,y,z):
  global laser_waist
  r2 = x**2 + y**2
  Wz = W0*sqrt(1.+(z/zr)**2) # Beam radius varies along propagation direction
  return (W0/Wz)*exp(-r2/Wz**2)
  
#-------------------------------------------------------------------------------
# set laser amplitude by combining the pulse shape, laser profile, and laser phase
#-------------------------------------------------------------------------------
def laser_func(x,y,t):
  global laser_amplitude,laser_phase,laser_profile,k0,w0
  z0 = Zf
  Rz = z0*(1.+(zr/z0)**2)
  E = laser_amplitude(t)*laser_profile(x,y,z0)
  r = sqrt(x*x+y*y)
  angle = w0*t+k0*z0-arctan(z0/zr)+k0*r**2/(2.*Rz)   
  return [E*cos(angle),E*sin(angle)]   # for circularly polarized laser
#  return [0,E*sin(angle)]               # for linearly polarized laser
#  return [E*sin(angle),0]               # for linearly polarized laser

#-------------------------------------------------------------------------------
# --- initialization of electromagnetic solver
#-------------------------------------------------------------------------------
if 0:
    em = EM3D(bounds=bounds,
          l_2dxz=l_2d,
          laser_func=laser_func,
          laser_emax=E0,
          dtcoef=0.995,
          stencil=1,
          l_pushf=0,
          l_setcowancoefs=1,
          laser_source_z=Zc)
else:
    em = EM3DPXR(spectral=0,
          bounds=bounds,
          l_2dxz=l_2d,
          laser_func=laser_func,
          laser_emax=E0,
          dtcoef=0.995,
          stencil=0,
          l_pushf=0,
          current_cor=0,
          l_setcowancoefs=0,
          laser_source_z=Zc)
  
registersolver(em)
em.finalize()

#-------------------------------------------------------------------------------
# --- definition and installation of plotting routine
#-------------------------------------------------------------------------------
def mkplots():
    if top.it%10==0: # --- plots every 10 time steps
      window(0);
      # --- 2D plot of Ey at top of frame
      fma();em.pfey(direction=1,l_transpose=1,view=9,cmin=-E0,cmax=E0);
      pldj([Zc],[w3d.xmmin],[Zc],[w3d.xmmax],color=green)
      pldj([Zc+Zf],[w3d.xmmin],[Zc+Zf],[w3d.xmmax],color=magenta)
      # --- get ey and bx (on processor 0)
      if l_2d:
        ey=em.gatherey()
        bx=em.gatherbx()
      else:
        ey=em.gatherey(direction=1,slice=w3d.ny/2)
        bx=em.gatherbx(direction=1,slice=w3d.ny/2)
      if me==0: # --- only on processor 0
        # --- get ey and bx on axis
        ey=ey[w3d.nx/2,:]
        bx=bx[w3d.nx/2,:]
        # --- scale by E0
        ey/=E0;bx/=E0
        # --- compute z locations
        z=w3d.zmmin+arange(w3d.nz+1)*w3d.dz
        z*=1.e6 # --- converts to microns
        # --- plot ey,bx and ey-bx*c on axis
#        plsys(10) # selects lower plot
#        pla(ey,z,width=4)     
#        pla(bx*clight,z,color=blue,width=4)     
#        pla(0.5*(ey-bx*clight),z,color=red)     
        # --- sets plot limits
#        limits(w3d.zmmin*1.e6,w3d.zmmax*1.e6,-1.5,1.5)
        # --- add title
#        ptitles('','z (microns)','E/E0')
        # --- refresh window
        refresh()
#installafterstep(mkplots)
 
# inject laser for some time
t0=time.clock()
start=MPI.Wtime()

em.step(200*ncells/32)  
#step(200*ncells/32)  

endt=MPI.Wtime()
print 'time = ',time.clock()-t0,endt-start
# now plot electric field from analytical formula under paraxial approximation 
# from http://www.rp-photonics.com/gaussian_beams.html

#winon(1)
#palette('bluewhitered.gp')
if l_2d:
    x0,z0 = getmesh2d(w3d.xmmin,w3d.dx,w3d.nx-1,
                    w3d.zmmin,w3d.dz,w3d.nz-1)
    z0-=Zc+Zf
    Rz = z0*(1.+(zr/z0)**2)
    E = E0*laser_profile(x0,x0*0.,z0)*cos(k0*z0-arctan(z0/zr)+k0*x0**2/(2.*Rz))
    ppg(transpose(E),view=10,cmin=-E0,cmax=E0,xmin=w3d.zmmin,xmax=w3d.zmmax,ymin=w3d.xmmin,ymax=w3d.xmmax)
    pldj([Zc],[w3d.xmmin],[Zc],[w3d.xmmax],color=green)
    pldj([Zc+Zf],[w3d.xmmin],[Zc+Zf],[w3d.xmmax],color=magenta)
    ptitles('E_y ^^analytic','z','x')
    refresh()

else:
    em.pfey()

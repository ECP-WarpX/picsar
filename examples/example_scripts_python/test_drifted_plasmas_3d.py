# ______________________________________________________________________________
#
# Execution test: Drifted plasmas
# We advice to not delete or modify this script, else make a copy
#
# In this script, electron-positron beams (drifted plasmas) are sent in each 
# direction of the domain: x,y,z.
#
# This script enables to test that neutral plasmas propagate correctly 
# and that particle communications work properly
#
# ______________________________________________________________________________

from warp import *
from warp.field_solvers.em3dsolverPXR import *
from warp.init_tools import *
from warp.data_dumping.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from mpi4py import MPI
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy import linalg as LA

def test_drifted_plasmas():
  """
  main function for the case test_drifted_plasmas
  
  """

  # Picsar flag: 0 warp routines, 1 picsar routines
  l_pxr=1
  l_temporal_diags=0
  l_pytest=0
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

  #-------------------------------------------------------------------------------
  # Parameters
  #-------------------------------------------------------------------------------

  # Number of iterations
  Nsteps = 100

  #Mesh: normalized at the plasma frequency
  dx=0.04
  dy=0.04
  dz=0.04
  dt=0.9*1./sqrt(1./dx**2+1./dy**2+1./dz**2)

  # Plasma properties
  plasma_elec_density = 0.01
  nppcell         = 5


  #Laser at the left border:
  a0             = 1.
  laser_duration = 10
  laser_width    = 4/2 #(=2w_0)
  #S-pol
  #Transverse profile - Gaussian
  #Longitudinal profile Super-Gaussian 12
  #focused at x=6 from the left border

  # --- scaling

  #-------------------------------------------------------------------------------
  # main parameters
  #-------------------------------------------------------------------------------
  dim = "3d"                 # 3D calculation
  #dim = "2d"                 # 2D calculation 
  #dim = "1d"                 # 1D calculation 
  dpi=100                     # graphics resolution
  l_test             = 1      # Will open output window on screen
                              # and stop before entering main loop.
  l_gist             = 0      # Turns gist plotting on/off
  l_matplotlib       = 1
  l_restart          = false  # To restart simulation from an old run (works?)
  restart_dump       = ""     # dump file to restart from (works?)
  l_moving_window    = 1      # on/off (Galilean) moving window
  l_plasma           = 1      # on/off plasma
  l_usesavedist      = 0      # if on, uses dump of beam particles distribution
  l_smooth           = 0      # on/off smoothing of current density
  l_laser            = 0      # on/off laser
  l_pdump            = 0      # on/off regular dump of beam data
  stencil            = 0      # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F 
                              # use 0 or 1; 2 does not verify Gauss Law
  if dim=="1d":stencil=0
  dtcoef             = 1.  # coefficient to multiply default time step that is set at the EM solver CFL
  top.depos_order    = 3      # particles deposition order (1=linear, 2=quadratic, 3=cubic)
  top.efetch         = 4      # field gather type (1=from nodes "momentum conserving"; 4=from Yee mesh "energy conserving")

  top.runid          = "Test_drifted_plasma"                         # run name
  top.pline1         = "Test"                         # comment line on plots
  top.runmaker       = "J.L. Vay, M. Lobet"                        # run makers
  top.lrelativ       = true                                # on/off relativity (for particles push)
  top.pgroup.lebcancel_pusher=0                         # flag for particle pusher (0=Boris pusher; 1=Vay PoP 08 pusher)
  #top.ibpush=2
  l_verbose          = 0                                   # verbosity level (0=off; 1=on)

  #-------------------------------------------------------------------------------
  # diagnostics parameters + a few other settings
  #-------------------------------------------------------------------------------
  live_plot_freq     = 1  # frequency (in time steps) of live plots (off is l_test is off)

  fielddiag_period   = 500
  partdiag_period    = 500

  #-------------------------------------------------------------------------------
  # laser parameters
  #-------------------------------------------------------------------------------
  # --- in lab frame
  lambda_laser       = 1.e-6                      # wavelength 
  laser_width       *= lambda_laser
  laser_waist        = laser_width/2  
  laser_radius       = laser_waist/2.354820        # FWHM -> radius
  laser_duration     = laser_duration*lambda_laser/clight # laser duration 
  laser_polangle     = 0#pi/2                       # polarization (0=aligned with x; pi/2=aligned with y)
  k0                 = 2.*pi/lambda_laser
  w0                 = k0*clight
  ZR                 = 0.5*k0*(laser_waist**2)            # Rayleigh length
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
  dens0             = plasma_elec_density*densc             # plasma density
  wp                = sqrt(dens0*echarge**2/(eps0*emass)) # plasma frequency
  kp                = wp/clight                           # plasma wavenumber
  lambda_plasma     = 2.*pi/kp                            # plasma wavelength
  Tplasma           = 2.*pi/wp                            # plasma period

  #-------------------------------------------------------------------------------
  # print some plasma parameters to the screen
  #-------------------------------------------------------------------------------
  print " dx: ",dx,"dy: ",dy,"dz: ",dz,"dt: ",dt
  print " Plasma elec. density: ",plasma_elec_density,'nc',dens0 
  print " Plasma wavelength: ",lambda_plasma
  print " Plasma frequency: ",wp
  print ' Plasma period:',Tplasma

  #-------------------------------------------------------------------------------
  # number of plasma macro-particles/cell
  #-------------------------------------------------------------------------------
  nppcellx = 1#5
  nppcelly = 1#5
  nppcellz = 1#5

  #-------------------------------------------------------------------------------
  # Algorithm choices
  # See Doxygen doc for more information
  #-------------------------------------------------------------------------------
  # Optional: current deposition algorithm, 
  currdepo=0
  # Optional: mpi com for the current desposition
  mpicom_curr=0
  # Field gathering method
  fieldgathe=0
  # Type of particle communication
  partcom =0
  # Field gathering and particle pusher together
  fg_p_pp_separated=0

  #-------------------------------------------------------------------------------
  # grid dimensions, nb cells and BC
  #-------------------------------------------------------------------------------
  w3d.zmmin = -2.*lambda_plasma
  w3d.zmmax = -w3d.zmmin
  w3d.xmmin = -2.*lambda_plasma
  w3d.xmmax = -w3d.xmmin
  w3d.ymmin = -2.*lambda_plasma
  w3d.ymmax = -w3d.ymmin

  w3d.nx = nint((w3d.xmmax-w3d.xmmin)/(dx*lambda_plasma))
  w3d.ny = nint((w3d.ymmax-w3d.ymmin)/(dy*lambda_plasma))
  w3d.nz = nint((w3d.zmmax-w3d.zmmin)/(dz*lambda_plasma))

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
  weight   = dens0 *w3d.dx*w3d.dy*w3d.dz/(nppcellx*nppcelly*nppcellz)
  top.wpid = nextpid() # Activate variable weights in the method addpart

  # --- create plasma species
  electron0 = Species(type=Electron,weight=weight,name='electron0')
  positron0 = Species(type=Positron,weight=weight,name='positron0')
  electron1 = Species(type=Electron,weight=weight,name='electron1')
  positron1 = Species(type=Positron,weight=weight,name='positron1')
  electron2 = Species(type=Electron,weight=weight,name='electron2')
  positron2 = Species(type=Positron,weight=weight,name='positron2')

  # --- Init the sorting
  sort = Sorting(periods=[10,10,10,10],starts=[0,0,0,0],activated=1,dx=0.5,dy=0.5,dz=0.5,xshift=-0.5,yshift=-0.5,zshift=-0.5)

  top.depos_order[...] = top.depos_order[0,0] # sets deposition order of all species = those of species 0
  top.efetch[...] = top.efetch[0] # same for field gathering

  #-------------------------------------------------------------------------------
  # set smoothing of current density
  #-------------------------------------------------------------------------------
  if l_smooth:
    # --- 1 time nilinear (0.25,0.5,0.25) + 1 time relocalization (-1, 3/2,-1.)
    npass_smooth = [[ 1 , 1 ],[ 0 , 0 ],[ 1 , 1 ]]
    alpha_smooth = [[ 0.5, 3./2],[ 0.5, 3.],[0.5, 3./2]]
    stride_smooth = [[ 1 , 1 ],[ 1 , 1 ],[ 1 , 1 ]]
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

      spatial_extent = 0.05
      speed = 0.9*clight

      zmin = w3d.zmmin*spatial_extent
      zmax = w3d.zmmax*spatial_extent
      xmin = w3d.xmmin*spatial_extent
      xmax = w3d.xmmax*spatial_extent
      ymin = w3d.ymmin*spatial_extent
      ymax = w3d.ymmax*spatial_extent
      
      np = w3d.nx*w3d.ny*nint((zmax-zmin)/w3d.dz)*nppcellx*nppcelly*nppcellz
      
      electron0.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,vzmean=speed,vthx=0.,vthy=0.,vthz=0.,spacing='uniform')
      positron0.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,vzmean=speed,vthx=0.,vthy=0.,vthz=0.,spacing='uniform')

      electron1.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,vymean=speed,vthx=0.,vthy=0.,vthz=0.,spacing='uniform')
      positron1.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,vymean=speed,vthx=0.,vthy=0.,vthz=0.,spacing='uniform')

      electron2.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,vxmean=speed,vthx=0.,vthy=0.,vthz=0.,spacing='uniform')
      positron2.add_uniform_box(np,xmin,xmax,ymin,ymax,zmin,zmax,vxmean=speed,vthx=0.,vthy=0.,vthz=0.,spacing='uniform')

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
  # initializes main field solver block
  #-------------------------------------------------------------------------------
  if l_pxr:
      ntilex = max(1,w3d.nxlocal/10)
      ntiley = max(1,w3d.nylocal/10)
      ntilez = max(1,w3d.nzlocal/10)
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
                   fg_p_pp_separated=fg_p_pp_separated,
                   l_verbose=l_verbose,
		   l_debug = 0)
      step = em.step
  else:
      print ' em=EM3D'
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
                   l_pushf=1,
                   l_verbose=l_verbose)

  #-------------------------------------------------------------------------------
  # restarts from dump file
  #-------------------------------------------------------------------------------
  if l_restart:
    restore(dump_file)

  #-------------------------------------------------------------------------------
  # register solver
  #-------------------------------------------------------------------------------
  print ' register solver'
  registersolver(em)
  em.finalize()
  loadrho()

  print 'done'

  # ______________________________________________________________________________
  # Diagnostics
  # ______________________________________________________________________________

  t = 0
  imu0 = 12.566370614E-7
  t_array = []
  e0kinE_array = []
  p0kinE_array = []
  e1kinE_array = []
  p1kinE_array = []
  e2kinE_array = []
  p2kinE_array = []
  ezE_array = []
  eyE_array = []
  exE_array = []
  bzE_array = []
  byE_array = []
  bxE_array = []
  div_array=[]
  etot_array = []
    
  # ______________________________________________________________________________
  # User actions at each iteration
  def liveplots():
    if top.it%live_plot_freq==0:
      fma(0); 
      em.pfez(view=3,titles=0,xscale=1./lambda_plasma,yscale=1./lambda_plasma,gridscale=1.e-12,l_transpose=1,direction=1)
      ptitles('Ez [TV/m]','Z [lambda]','X [lambda]')
      
      density=electron0.get_density()
      density=density[:,w3d.ny/2,:]
      ppg(transpose(density),view=4,titles=0,xmin=w3d.zmmin+top.zgrid,xmax=w3d.zmmax+top.zgrid,ymin=w3d.xmmin,ymax=w3d.xmmax,xscale=1e6,yscale=1.e6)#,gridscale=1./dens0)

      density=positron0.get_density()
      density=density[:,w3d.ny/2,:]
      ppg(transpose(density),view=5,titles=0,xmin=w3d.zmmin+top.zgrid,xmax=w3d.zmmax+top.zgrid,ymin=w3d.xmmin,ymax=w3d.xmmax,xscale=1e6,yscale=1.e6)#,gridscale=1./dens0)
   
   
    if (l_temporal_diags):
     
      # With Picsar
      if l_pxr:
        # Kinetic energy  
        # Using python
        if False:
          e0kinE = 0
	  e1kinE = 0
          e2kinE = 0
          p0kinE = 0
          p1kinE = 0
          p2kinE = 0
          for i in range(em.ntilez):
            for j in range(em.ntiley):
              for k in range(em.ntilex):
                  e0kinE += sum((1./electron0.pgroups[i][j][k].gaminv[0:electron0.pgroups[i][j][k].nps]-1.)*electron0.pgroups[i][j][k].pid[0:electron0.pgroups[i][j][k].nps,0])*9.10938215E-31*299792458.**2
                  e1kinE += sum((1./electron1.pgroups[i][j][k].gaminv[0:electron1.pgroups[i][j][k].nps]-1.)*electron1.pgroups[i][j][k].pid[0:electron1.pgroups[i][j][k].nps,0])*9.10938215E-31*299792458.**2
                  e2kinE += sum((1./electron2.pgroups[i][j][k].gaminv[0:electron2.pgroups[i][j][k].nps]-1.)*electron2.pgroups[i][j][k].pid[0:electron2.pgroups[i][j][k].nps,0])*9.10938215E-31*299792458.**2
                  p0kinE += sum((1./positron0.pgroups[i][j][k].gaminv[0:positron0.pgroups[i][j][k].nps]-1.)*positron0.pgroups[i][j][k].pid[0:positron0.pgroups[i][j][k].nps,0])*9.10938215E-31*299792458.**2
                  p1kinE += sum((1./positron1.pgroups[i][j][k].gaminv[0:positron1.pgroups[i][j][k].nps]-1.)*positron1.pgroups[i][j][k].pid[0:positron1.pgroups[i][j][k].nps,0])*9.10938215E-31*299792458.**2
                  p2kinE += sum((1./positron2.pgroups[i][j][k].gaminv[0:positron2.pgroups[i][j][k].nps]-1.)*positron2.pgroups[i][j][k].pid[0:positron2.pgroups[i][j][k].nps,0])*9.10938215E-31*299792458.**2
        else:
          # Using the fortran subroutine
          e0kinE = em.get_kinetic_energy(1)*9.10938215E-31*299792458.**2 
          e1kinE = em.get_kinetic_energy(3)*9.10938215E-31*299792458.**2 
          e2kinE = em.get_kinetic_energy(5)*9.10938215E-31*299792458.**2 
          p0kinE = em.get_kinetic_energy(2)*9.10938215E-31*299792458.**2 
          p1kinE = em.get_kinetic_energy(4)*9.10938215E-31*299792458.**2
          p2kinE = em.get_kinetic_energy(6)*9.10938215E-31*299792458.**2

        print ' Electron kinetic energy (J)',e0kinE,e1kinE,e2kinE
        print ' Positron kinetic energy (J)',p0kinE,p1kinE,p2kinE    
        
        # Electric energy
        # Using python
        if False:  
          ezE = 0.5*sum(em.fields.Ez[0:em.nxlocal,0:em.nylocal,0:em.nzlocal]**2)*em.dx*em.dy*em.dz*8.85418782E-12
          exE = 0.5*sum(em.fields.Ex[0:em.nxlocal,0:em.nylocal,0:em.nzlocal]**2)*em.dx*em.dy*em.dz*8.85418782E-12
          eyE = 0.5*sum(em.fields.Ey[0:em.nxlocal,0:em.nylocal,0:em.nzlocal]**2)*em.dx*em.dy*em.dz*8.85418782E-12
          bzE = 0.5*sum(em.fields.Bz[0:em.nxlocal,0:em.nylocal,0:em.nzlocal]**2)*em.dx*em.dy*em.dz*imu0
          bxE = 0.5*sum(em.fields.Bx[0:em.nxlocal,0:em.nylocal,0:em.nzlocal]**2)*em.dx*em.dy*em.dz*imu0
          byE = 0.5*sum(em.fields.By[0:em.nxlocal,0:em.nylocal,0:em.nzlocal]**2)*em.dx*em.dy*em.dz*imu0
          div = 0
          print ' Electric energy (J):',ezE,eyE,exE
          print ' Magnetic energy (J):',bzE,byE,bxE
        # using picsar fortran subroutine
        else:
          exE = em.get_field_energy('ex')*8.85418782E-12 
          eyE = em.get_field_energy('ey')*8.85418782E-12 
          ezE = em.get_field_energy('ez')*8.85418782E-12 
          bxE = em.get_field_energy('bx')*8.85418782E-12 
          byE = em.get_field_energy('by')*8.85418782E-12 
          bzE = em.get_field_energy('bz')*8.85418782E-12   
          div = em.get_normL2_divEeps0_rho()  
          print ' Electric energy (J):',ezE,eyE,exE
          print ' Magnetic energy (J):',bzE,byE,bxE
          print ' NormL2 of DivE*eps0 - rho:',div
        
        etot = ezE + exE + eyE 
        etot += e0kinE + e1kinE + e2kinE
        etot += p0kinE + p1kinE + p2kinE 
        etot += bzE + bxE + byE
        print ' Total energy (J)',etot
        
        # Put in arrays
        e0kinE_array.append(e0kinE)
        p0kinE_array.append(p0kinE)
        e1kinE_array.append(e1kinE)
        p1kinE_array.append(p1kinE)
        e2kinE_array.append(e2kinE)
        p2kinE_array.append(p2kinE)
        ezE_array.append(ezE)
        eyE_array.append(eyE)
        exE_array.append(exE)
        bzE_array.append(bzE)
        byE_array.append(byE)
        bxE_array.append(bxE)  
        div_array.append(div)
        etot_array.append(etot)
        t_array.append(top.time)

      #With warp      
      else:
      
        # Get divergence
        divE = em.getdive()
        rho = em.getrho()
        F = em.getf()
        divarray = divE*eps0 - rho
        div = LA.norm((divarray))
        maxdiv = abs(divarray).max()
        mindiv = abs(divarray).min()
          
        print ' NormL2 of DivE*eps0 - rho:',div,LA.norm(F),'Max',maxdiv,'Min',mindiv
      
        #electron0.getn()      # selectron macro-particle number
        #x = electron0.getx()      # selectron macro-particle x-coordinates 
        #y = electron0.gety()      # selectron macro-particle y-coordinates      
        #z = electron0.getz()      # selectron macro-particle x-coordinates
        #pid = electron0.getpid()    
        
        div_array.append(div)
        t_array.append(top.time)
      
  installafterstep(liveplots)

  # _____________________________________________________________________________
  #
  # Outputs
  # _____________________________________________________________________________

  # Setup the diagnostics
  diag_period = 10

  remove_existing_directory( ['diags'] )

  diag1 = FieldDiagnostic( period=diag_period, top=top, w3d=w3d, em=em,\
    comm_world=comm_world , write_dir='diags')

  installafterstep( diag1.write )

  diag2 = ParticleDiagnostic( period=diag_period, top=top, w3d=w3d, \
    species={ species.name : species for species in listofallspecies }, \
    comm_world=comm_world , write_dir='diags')

  installafterstep( diag2.write )

  tottime = AppendableArray()
  def accuttime():
    global tottime
    tottime.append(time.clock())
    if me==0 and top.it%200==0:
      f=PW.PW('tottime.pdb')
      f.time=tottime[:]
      f.close()

  #installafterstep(accuttime)


  print '\nInitialization complete\n'

  # ______________________________________________________________________________
  # running
  em.step(Nsteps,1,1)
    
  # ______________________________________________________________________________
  #
  # Checkings and analysis
  # ______________________________________________________________________________

  if me==0:
 
    if (l_temporal_diags):
   
      # ______________________________________________________________________________
      # RCparams
      if l_matplotlib:
        mpl.rcParams['font.size'] = 15
        mpl.rcParams['legend.fontsize'] = 12
        mpl.rcParams['figure.facecolor'] = 'white'

      t_array=array(t_array)

      print ' Plasma period:',Tplasma

      # Plotting    
      if l_matplotlib:  
        fig = plt.figure(figsize=(16,8))
        gs = gridspec.GridSpec(2, 2)
        ax = plt.subplot(gs[:, :])

        ax.plot(t_array,e0kinE_array,label='Elec. kin. energy')
        ax.plot(t_array,p0kinE_array,label='Posi. kin. energy')
        ax.plot(t_array,e1kinE_array,label='Elec. kin. energy')
        ax.plot(t_array,p1kinE_array,label='Posi. kin. energy')
        ax.plot(t_array,e2kinE_array,label='Elec. kin. energy')
        ax.plot(t_array,p2kinE_array,label='Posi. kin. energy')
        ax.plot(t_array,ezE_array,label='Ez energy')
        ax.plot(t_array,eyE_array,label='Ey energy')
        ax.plot(t_array,exE_array,label='Ex energy')
        ax.plot(t_array,bzE_array,label='Bz energy')
        ax.plot(t_array,byE_array,label='By energy')
        ax.plot(t_array,bxE_array,label='Bx energy')
        ax.plot(t_array,etot_array,label='Total energy',color='k',ls=':')

        ax.set_ylim([0,max(etot_array)*1.1])

        ax.set_xlabel('t (s)')
        ax.set_ylabel('Energy (J)')

        ax.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
      
      # _____________________________________________________________________
      # DivE*eps0 - rho
      if l_matplotlib:
        fig1 = plt.figure(figsize=(12,8))
        gs1 = gridspec.GridSpec(2, 2)
        ax1 = plt.subplot(gs[:, :])    

        ax1.plot(t_array,div_array,label=r'$\nabla E \times \varepsilon_0 - \rho$',lw=2) 
        ax1.legend(loc='upper center',ncol=4,borderaxespad=-2,fontsize=20)
      
        ax1.set_xlabel('t (s)')
      
      if l_pytest: assert (max(div_array) < 1E-5),"L2 norm||DivE - rho/eps0|| too high"
    
    # _____________________________________________________________________
    # Advice
    
    print
    print ' _______________________________________'
    print ' Advice for users:' 
    print ' - Check that the energy is constant with time'
    print ' - Check that divE = rho/eps0 for each tests'
    print ' - Check the energy oscillating behavior'
    print    
    
    if l_matplotlib: plt.show()
  
if __name__ == "__main__":

	# Launch the test
  test_drifted_plasmas() 

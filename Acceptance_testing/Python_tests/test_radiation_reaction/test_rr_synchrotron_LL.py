# --- Input script to test for RR losses 
# --- in a case of 'synchrotron radiation'

# --- Import needed python modules (WARP, PXR, etc.) 
from warp import *
from warp.field_solvers.em3dsolverPXR import *
import os
from warp.data_dumping.openpmd_diag import FieldDiagnostic, ParticleDiagnostic, \
    ParticleAccumulator, ProbeParticleDiagnostic
from mpi4py import MPI

def test_rr_synchrotron_LL():
	# ----------------------------------------------------------------------------------------
	# --- CODE PARAMETERS ---- (USER: PLEASE IGNORE THIS PART) 
	# ----------------------------------------------------------------------------------------
	# --- flags turning off unnecessary diagnostics (USER: please ignore)
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
	top.lcomm_cartesian=1
	# Flags turning off auto decomp and let user specify its decomp
	top.lautodecomp = 1 # Particles
	top.lfsautodecomp = 1 # fields
	# Flags turning on/off load balancing
	load_balance=0
	dlb_freq=100
	dlb_threshold=10 # dynamic load balancing threshold in %
	dlb_at_init=1 # Do a load balancing of the simulation at init

	# ----------------------------------------------------------------------------------------
	# --- USER PARAMETERS ------- 
	# ----------------------------------------------------------------------------------------

	# Turn PXR booster (ON/OFF) -- ON by default 
	l_pxr=1

	# Case dimension 
	dim = "2d"  # ("1d", "2d" or "3d" geometries)

	# Laser wavelength (use for normalization of input params)
	lambda_laser       = 1e-6     # laser wavelength (in meters)
	b0 				   = 100.      #  Normalized amplitude of the constant magnetic field

								  # (a0>1: relativistic motion)

	# Maxwell solver parameters 
	spectral = 1 # (1 - PSATD pseudo-spectral solver, 0 - FDTD solver)
	stencil = 1  # (if spectral is 0 : 0 - FDTD-Yee, 1 - FDTD KK)

	# Particle pusher type 
	particle_pusher = 4  # (0 : Boris, 1: Vay, 2: Boris with RR (S09), 3: Boris with RR (B08), 4: Boris with RR (LL),)

	# Scale factor (spatial, density) to scale down simulation 
	dfact = 1. # Scale factor in density
	dxfact = 1. # Scale factor in space
	dtfact = 1. # Scale factor in time 

	# Mesh resolution 
	coeff_yee=0.995 # Coefficients to apply to the ratio c delta t/ delta x
	dx=dy=dz=0.01   # Grid resolutions in laser wavelength
	dt=dx
	epsi=1e-3  #relative error tolerance

	# Axis origin and simulation box size (in laser wavelength)
	x_min=-0.5
	y_min=-0.5
	z_min=-0.5
	L_x=2.5
	L_y=2.5
	L_z=2.5

	# Field and particle boundary conditions 
	field_bc_xy = periodic
	field_bc_z  = periodic
	part_bc_xy  = periodic
	part_bc_z   = periodic

	# Simulation duration (in lase periods)
	t_max=16

	# Number of particles per cell (here e- only)
	nppcellx_e = 1 # Along X
	nppcelly_e = 1 # Along Y 
	nppcellz_e = 1 # Along Z

	# Initial conditions for electron 
	px0 = 100.
	py0 = 0.
	pz0 = 0.
	x0  = 0. 
	y0  = 0. 
	z0  = 0. 

	# Current deposition and field gathering 
	current_cor = 0 # Flag for current correction (spectral mode)
	depos_order = 3 # particles deposition order (Particle shape factor)
					# (1=linear, 2=quadratic, 3=cubic)

	# -- Control flags 
	l_test             = 0      # Do code init without starting PIC cycle 
	l_plasma           = 1      # Turns on/off plasma initialization 
	l_smooth           = 1      # on/off smoothing of current density
	l_laser            = 1      # turns on/off laser
	l_gist             = 0      # Turns gist plotting on/off
	l_restart          = false  # To restart simulation from an old run (works?)
	l_moving_window    = 1      # on/off (Galilean) moving window
	l_usesavedist      = 0      # if on, uses dump of beam particles distribution
	l_pdump            = 0      # on/off regular dump of beam data
	l_verbose          = 0      # verbosity level (0=off; 1=on)
	restart_dump       = ""     # dump file to restart from (works?)
	dpi=100                     # graphics resolution
	# Flags controlling I/Os parameters
	live_plot_freq     = 1000000000 # frequency (in time steps) of live plots 
									# (off is l_test is off)
	fielddiag_period      = 10.      # In laser period
	partdiag_period       = 0.0001
	lparallelo = False
	lpar_metadata= False  
	onefile_per_flush = True 


	# ----------------------------------------------------------------------------------------
	# --- CODE INIT (USER: PLEASE IGNORE THIS PART) ------- 
	# ----------------------------------------------------------------------------------------
	# Init Maxwell solver parameters 
	if spectral:
		stencil=1
		norderx=100
		nordery=100
		norderz=100
		ntsub=inf
		current_cor=0
		nxguard=8
		nyguard=8
		nzguard=8
	else:
		norderx=2
		nordery=2
		norderz=2
		ntsub=1
		current_cor=0
		nxguard=3
		nyguard=3
		nzguard=3
	
	# --- scaling of numerical parameters
	dx*=dxfact
	dy=dz=dx
	dt*=dtfact
	if dim=="1d":stencil=0
	dtcoef             = dt/dx  # coefficient to multiply default time step 
								# that is set at the EM solver CFL
	if spectral == 0 and stencil ==0 : dtcoef*=dtcoef*coeff_yee
	top.depos_order    = depos_order  
	top.efetch         = 4      # field gather type 
								#(1=from nodes "momentum conserving"; 
								# 4=from Yee mesh "energy conserving")
	top.runid          = "Test of classical radiation reaction "   # run name
	top.pline1         = "RR effect "                              # comment line on plots
	top.runmaker       = "H. Vincenti"                             # run makers
	top.lrelativ       = true                                      # on/off relativity 
																   #(for particles push)
	top.pgroup.lebcancel_pusher=particle_pusher                    

	# --- Laser parameters 
	las_period         = lambda_laser/clight
	k0                 = 2.*pi/lambda_laser
	w0                 = k0*clight

	# --- Plasma parameters 
	densc             = emass*eps0*w0**2/echarge**2         # critical density

	if dim=="2d":
	  nppcelly_e= nppcelly_e = 1
	if dim=="1d":
	  nppcellx_e = nppcellx_e = 1
	  nppcelly_e = nppcelly_e = 1

	# --- Grid mesh 
	w3d.xmmin = x_min*lambda_laser
	w3d.xmmax = w3d.xmmin+L_x*lambda_laser
	w3d.ymmin = y_min*lambda_laser
	w3d.ymmax = w3d.ymmin+L_y*lambda_laser
	w3d.zmmin = z_min*lambda_laser
	w3d.zmmax = w3d.zmmin+L_z*lambda_laser

	w3d.dx=dx*lambda_laser
	w3d.dy=dy*lambda_laser
	w3d.dz=dz*lambda_laser

	# --- Adjustements to get a power of 2 cells in each direction 
	# --- (only is spectral flag activated)
	if (spectral): 
		w3d.nx = nint((w3d.xmmax-w3d.xmmin)/w3d.dx)
		if (w3d.nx/64. is not nint(w3d.nx/64)):
			w3d.nx=64*nint(w3d.nx/64)
		w3d.nz = nint((w3d.zmmax-w3d.zmmin)/w3d.dz)
		if (w3d.nz/64. is not nint(w3d.nz/64)): 
			w3d.nz=64*nint(w3d.nz/64)
		w3d.ny = nint((w3d.ymmax-w3d.ymmin)/w3d.dy)
		if (w3d.ny/64. is not nint(w3d.ny/64)):
			w3d.ny=64*max(1,nint(w3d.ny/64))
		w3d.xmmax = w3d.xmmin+w3d.nx*w3d.dx
		w3d.ymmax = w3d.ymmin+w3d.ny*w3d.dy
		w3d.zmmax = w3d.zmmin+w3d.nz*w3d.dz
	if dim in ["1d"]:
		w3d.nx = 2
		w3d.xmmin = -float(w3d.nx)/2.
		w3d.xmmax = float(w3d.nx)/2.
		w3d.dx=1.
		nxguard=0
	if dim in ["1d","2d"]:
		w3d.ny = 2
		w3d.ymmin = -float(w3d.ny)/2.
		w3d.ymmax = float(w3d.ny)/2.
		w3d.dy=1.
		nyguard=0

	print 'w3d.dx, w3d.dy, w3d.dz =', w3d.dx, w3d.dy, w3d.dz
	print 'w3d.nx, w3d.ny, w3d.nz =', w3d.nx, w3d.ny, w3d.nz
	print 'w3d.xmmax, w3d.ymmax, w3d.zmmax =', w3d.xmmax, w3d.ymmax, w3d.zmmax

	# --- Sets field boundary conditions
	# longitudinal
	w3d.bound0  = w3d.boundnz = field_bc_z
	# transverse
	w3d.boundxy = field_bc_xy

	# --- sets particles boundary conditions
	# longitudinal
	top.pbound0  = part_bc_z
	top.pboundnz = part_bc_z
	# transverse
	top.pboundxy = part_bc_xy

	top.wpid = top.nextpid() # Activate variable weights in the method addpart
        if(l_pxr):
          top.nextpid()
          top.nextpid()
	# --- create plasma species
	electron = Species(type=Electron,weight=1.,name='electron')

	# --- sets deposition order of all species = those of species 0
	top.depos_order[...] = top.depos_order[0,0] 

	# --- same for field gathering
	top.efetch[...] = top.efetch[0] 
	if dim in ["1d","2d"]:
	  top.depos_order[1,:]=1
	if dim=="1d":
	  top.depos_order[0,:]=1

	# --- Sets smoothing of current densities
	if l_smooth:
	  # - 1 time bilinear (0.25,0.5,0.25) + 1 time relocalization (-1, 3/2,-1.)
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

	# --- Initializes WARP 
	top.fstype = -1 # sets field solver to None (desactivates electrostatic solver)
	package('w3d');generate()
	pg = top.pgroup


	# --- Sets offset between particle and grid (for PMLs)
	offset_x_part_grid=[0.,0.]
	offset_y_part_grid=[0.,0.]
	offset_z_part_grid=[0.,0.]

	if (top.pboundxy==absorb):
		offset_x_part_grid=[max(1,int(0.2*lambda_laser/w3d.dx))*w3d.dx,
		-max(1,int(0.2*lambda_laser/w3d.dx))*w3d.dx]
		if dim=="3d":
		  offset_y_part_grid=[max(1,int(0.2*lambda_laser/w3d.dy))*w3d.dy,
		  -max(1,int(0.2*lambda_laser/w3d.dy))*w3d.dy]
	if (top.pboundnz==absorb):
		offset_z_part_grid[1]=-max(1,int(0.2*lambda_laser/w3d.dz))*w3d.dz
	if (top.pbound0==absorb):
		offset_z_part_grid[0]=max(1,int(0.2*lambda_laser/w3d.dz))*w3d.dz

	# --- Initializes electron position 
	if l_plasma:
		# Computes electron gamma and velocities
		gam0=np.sqrt(1.+px0**2+py0**2+pz0**2)
		vx0=(px0/gam0)*clight
		vy0=(py0/gam0)*clight
		vz0=(pz0/gam0)*clight
		# Normalizes position 
		x0=x0*lambda_laser
		y0=y0*lambda_laser
		z0=z0*lambda_laser
		# Add 1 electron to "electron" species object 
		electron.addpart(x=x0,y=y0,z=z0,       \
						 vx=vx0,vy=vy0,vz=vz0, \
						 w=1.)

	# --- Initializes electromagnetic solver object em 
	if l_pxr:
		if (dim=='3d'):
			ntilex = max(1,w3d.nxlocal/8)
			ntiley = max(1,w3d.nylocal/8)
			ntilez = max(1,w3d.nzlocal/8)
		else: # dim is 2D
			ntilex = max(1,w3d.nxlocal/30)
			ntiley = 1
			ntilez = max(1,w3d.nzlocal/30)
		em = EM3DPXR(laser_func=None,
					 stencil=stencil,
					 npass_smooth=npass_smooth,
					 alpha_smooth=alpha_smooth,
					 stride_smooth=stride_smooth,
					 l_2dxz=dim=="2d",
					 l_1dz=dim=="1d",
					 dtcoef=dtcoef,
					 l_getrho=0,
					 spectral=spectral,
					 norderx=norderx,
					 nordery=nordery,
					 norderz=norderz,
					 nxguard=nxguard,
					 nyguard=nyguard,
					 nzguard=nzguard,
					 ntsub=ntsub,
					 current_cor=current_cor,
					 listofallspecies=listofallspecies,
					 ntilex=ntilex,
					 ntiley=ntiley,
					 ntilez=ntilez,
					 offset_x_part_grid=offset_x_part_grid,
					 offset_y_part_grid=offset_y_part_grid,
					 offset_z_part_grid=offset_z_part_grid,
					 #l_reinject=[1.,1.,1.,1.,1.,1.],
					 l_verbose=l_verbose,
					 dload_balancing=load_balance,
					 dlb_freq=dlb_freq,
					 dlb_threshold=dlb_threshold,
					 dlb_at_init=dlb_at_init,
					 l_setcowancoefs=True,
					 mpi_buf_size=80,
					 partcom=1)
		step = em.step
	else:
		em = EM3DFFT(laser_func=None,
					 stencil=stencil,
					 npass_smooth=npass_smooth,
					 alpha_smooth=alpha_smooth,
					 stride_smooth=stride_smooth,
					 l_2dxz=dim=="2d",
					 l_1dz=dim=="1d",
					 dtcoef=dtcoef,
					 l_getrho=0,
					 spectral=spectral,
					 norderx=norderx,
					 nordery=nordery,
					 norderz=norderz,
					 nyguard=nyguard, 
					 nzguard=nzguard, 
					 ntsub=ntsub,
					 current_cor=current_cor,
					 l_verbose=l_verbose,
					 l_setcowancoefs=True)


	# --- Restarts from dump file 
	if l_restart:
	  restore(dump_file)

	# --- Register em solver 
	print 'register solver'
	registersolver(em)
	em.finalize()
	loadrho()
	print 'done'

	# Initializes constant magnetic field 
	B_max=b0*emass*w0/echarge
	em.fields.By[...] = B_max       

	# --- Normalize output freq quantities  
	fielddiag_period       = max(int(fielddiag_period*las_period /top.dt),1)
	partdiag_period        = max(int(partdiag_period*las_period /top.dt),1)

	# --- Compute number of time steps 
	N_step = int(t_max*las_period/top.dt)

	# Comparison with analytical formula 
	# Calculation of important quantities and comparison with formula
	def calc_data():
		# Constants
		c=299792458 
		me=9.11e-31
		mi=1836.*me
		eps0= 8.854187817620e-12
		mu0=1.25663706e-6
		e=1.6e-19;

		las_lambda=1e-6
		las_time=las_lambda/c
		las_omega=2*np.pi/las_time
		cnob=1/(me*las_omega/e)                   # normalisation B --> unites a0
		cnop=1/(me*c)
	
		B=b0/cnob

		K=mu0*e**4*B**2/(6*np.pi*me**3*c)

		# box dimensions in laser wavelength	
		N=np.int(np.ceil(t_max/dz)) # number of iterations

		# Electron x,u init
		t=np.linspace(0,t_max*las_time,N)
		gam_a=np.zeros(N)
		ux0=px0*c # in si units
		uz0=pz0*c # in si units
		gam_a[0]=np.sqrt(1+(ux0/c)**2+(uz0/c)**2)

		A=(gam_a[0]+1)/(gam_a[0]-1)

		gam_a=(A+np.exp(-K*t*2))/(A-np.exp(-K*t*2))
		return(gam_a)
	gam_a=calc_data()  
	
	failure=False
	def check_model(): 
		#global gam_a, epsi, failure
		c=299792458 
		t=top.it
		gam=np.sqrt(1+(electron.getux()**2+electron.getuy()**2+electron.getuz()**2)/c**2)
		#print "time %s, gam_a, %s,  gam, %s."%(t,gam_a[t-1],gam)
		if np.abs((gam_a[t-1]-gam)/gam) > epsi : #relative error
			failure=True
		

	# Comparison with analytical formula    
	installafterstep(check_model)   


	print '\nInitialization complete\n'

	# if this is a test, then stop, else execute main loop
	if l_test:
	  if(me==0):
		print '<<< To execute n steps, type "step(n)" at the prompt >>>'
	else:
	  if (me==0):
		  print("Total number of steps: ",N_step)
		  tdeb=MPI.Wtime()
	  em.step(N_step)
	  tend=MPI.Wtime()
	  if(me==0):
		  print("Final runtime (s): "+str(tend-tdeb))

  	# Assert failure 
  	assert failure == False


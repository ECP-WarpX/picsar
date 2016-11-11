"""Class for 2D & 3D FFT-based electromagnetic solver """
from warp.field_solvers.em3dsolverFFT import *
from warp.particles.species import *

try:
    #import warp.field_solvers.GPSTD as gpstd
	import GPSTDPXR as gpstd
except:
	#import GPSTDPXR as gpstd
    import warp.field_solvers.GPSTD as gpstd
try:
    import picsarpy as pxrpy
    pxr = pxrpy.picsar
    l_pxr=True
except:
    l_pxr=False
try:
    from mpi4py import MPI
except:
    print 'Error cannot import mpi4py'
try:
    import os as os
except:
    print 'Error cannot import os'
try:
    # Try to import fortran wrapper of FFTW
    # import pyfftw
    # fft = pyfftw.interfaces.numpy_fft
    import fastfftforpy as fftpy
    import fastfftpy as fstpy
    fst=fstpy.fastfft
    fft=fftpy
    l_fftw=True
except:
    fft = np.fft
    l_fftw=False



class EM3DPXR(EM3DFFT):

    __em3dpxrinputs__ = []
    __flaginputs__ = {'ntilex':1,
                      'ntiley':1,
                      'ntilez':1,
                      'listofallspecies':[],
                      'dload_balancing':0,
                      'dlb_freq':1,
                      'dlb_threshold':20,
                      'dlb_at_init':1,
                      'it_dlb_init':11,
                      'l_output_grid':0,
                      'l_output_freq':1,
                      'rhodepo':0,      # Charge deposition method
                      'currdepo':0,     # Current deposition method
                      'mpicom_curr':1,  # Com type Current deposition
                      'fieldgathe':0,   # Field gathering method
                      'partcom':0,      # Particle communication
                      'fg_p_pp_separated':0,
                      'lvec_curr_depo':8,
                      'lvec_charge_depo':64,
                      'lvec_fieldgathe':512,
                      'sorting':None,
                      'l_debug':0
                      }

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass


        self.processdefaultsfromdict(EM3DPXR.__flaginputs__,kw)

        if (self.l_debug):
          print("Call __init__")
          print(' Debug prints activated')

        EM3DFFT.__init__(self,kwdict=kw)

        self.l_pxr = l_pxr
        self.l_fftw = l_fftw

        # If sorting undefined
        if self.sorting==None:
          self.sorting = Sorting([],[],activated=0,dx=1.,dy=1.,dz=1.,xshift=-0.5,yshift=-0.5,zshift=-0.5)

        if (self.l_debug): print("End __init__")

    def finalize(self,lforce=False):
        if self.finalized and not lforce: return
        if self.l_pxr:
        	EM3D.finalize(self)
        	self.allocatefieldarraysFFT()
        	self.allocatefieldarraysPXR()
        else:
          EM3DFFT.finalize(self)

    def convertindtoproc(self,ix,iy,iz,nx,ny,nz):
      ixt = ix
      iyt = iy
      izt = iz

      if (ixt < 0   ): ixt = nx - 1
      if (ixt > nx-1): ixt = 0
      if (iyt < 0   ): iyt = ny - 1
      if (iyt > ny-1): iyt = 0
      if (izt < 0   ): izt = nz - 1
      if (izt > nz-1): izt = 0

      convertindextoproc = ixt + iyt*nx + izt*nx*ny

      return convertindextoproc


    def allocatefieldarraysPXR(self):

        if (self.l_debug): print("allocatefieldarraysPXR")

        # Set up case dimensionality
        if (self.l_2dxz):
          pxr.c_dim=2
        else:
          pxr.c_dim=3
          
        # Set up PXR MPI Data
        if (self.l_debug): print(" Setup PXR MPI Data")
        pxr.nprocx=top.fsdecomp.nxprocs
        pxr.nprocy=top.fsdecomp.nyprocs
        pxr.nprocz=top.fsdecomp.nzprocs
        ixcpu=top.fsdecomp.ixproc
        iycpu=top.fsdecomp.iyproc
        izcpu=top.fsdecomp.izproc
        for iz in range(-1,2):
            for iy in range(-1,2):
                for ix in range(-1,2):
                    indtoproc=self.convertindtoproc(ixcpu+ix,iycpu+iy,izcpu+iz,pxr.nprocx,pxr.nprocy,pxr.nprocz)
                    pxr.neighbour[ix+1,iy+1,iz+1]=indtoproc

        if (ixcpu==0):
            pxr.x_min_boundary=1
        if (ixcpu==pxr.nprocx-1):
            pxr.x_max_boundary=1
        if (iycpu==0):
            pxr.y_min_boundary=1
        if (iycpu==pxr.nprocy-1):
            pxr.y_max_boundary=1
        if (izcpu==0):
            pxr.z_min_boundary=1
        if (izcpu==pxr.nprocz-1):
            pxr.z_max_boundary=1

        pxr.x_coords=ixcpu
        pxr.y_coords=iycpu
        pxr.z_coords=izcpu

        # MPI boundaries index in global array
        if (self.l_debug): print(" MPI boundaries index in global array")
        pxr.cell_x_min=top.fsdecomp.ix
        pxr.cell_y_min=top.fsdecomp.iy
        pxr.cell_z_min=top.fsdecomp.iz
        pxr.cell_x_max=pxr.cell_x_min+(top.fsdecomp.nx-1)
        pxr.cell_y_max=pxr.cell_y_min+(top.fsdecomp.ny-1)
        pxr.cell_z_max=pxr.cell_z_min+(top.fsdecomp.nz-1)
        pxr.x_grid_mins=top.fsdecomp.xmin
        pxr.x_grid_maxs=top.fsdecomp.xmax
        pxr.y_grid_mins=top.fsdecomp.ymin
        pxr.y_grid_maxs=top.fsdecomp.ymax
        pxr.z_grid_mins=top.fsdecomp.zmin
        pxr.z_grid_maxs=top.fsdecomp.zmax

        # Particle boundaries for PXR
        if (self.l_debug): print(" Setup particle boundaries for PXR")
        if (top.pbound0 == absorb):
        	pxr.pbound_z_min=1
        elif(top.pbound0 == reflect):
        	pxr.pbound_z_min=2
        else: # Default is periodic
        	pxr.pbound_z_min=0

        if (top.pboundnz == absorb):
            pxr.pbound_z_max=1
        elif(top.pboundnz == reflect):
            pxr.pbound_z_max=2
        else: # Default is periodic
            pxr.pbound_z_max=0

        if (top.pboundxy == absorb):
            pxr.pbound_x_min=1
            pxr.pbound_x_max=1
            pxr.pbound_y_min=1
            pxr.pbound_y_max=1
        elif(top.pboundxy == reflect):
            pxr.pbound_x_min=2
            pxr.pbound_x_max=2
            pxr.pbound_y_min=2
            pxr.pbound_y_max=2
        else: # Default is periodic
            pxr.pbound_x_min=0
            pxr.pbound_x_max=0
            pxr.pbound_y_min=0
            pxr.pbound_y_max=0

        # --- number of grid cells
        if (self.l_debug): (" Setup number of grid cells for PXR")
        pxr.nx_global = w3d.nx
        pxr.ny_global = w3d.ny
        pxr.nz_global = w3d.nz
        pxr.nx_global_grid = pxr.nx_global+1
        pxr.ny_global_grid = pxr.ny_global+1
        pxr.nz_global_grid = pxr.nz_global+1

        pxr.nx = self.nxlocal
        pxr.ny = self.nylocal
        pxr.nz = self.nzlocal
        pxr.nx_grid=pxr.nx+1
        pxr.ny_grid=pxr.ny+1
        pxr.nz_grid=pxr.nz+1


        # --- number of guard cells
        if (self.l_debug): print(" Setup number of guard cells for PXR")
        pxr.nxguards = self.nxguard
        pxr.nyguards = self.nyguard
        pxr.nzguards = self.nzguard
        pxr.nxjguards = self.nxguard
        pxr.nyjguards = self.nyguard
        pxr.nzjguards = self.nzguard

        # --- grid dimensions
        pxr.xmin = w3d.xmmin
        pxr.ymin = w3d.ymmin
        pxr.zmin = w3d.zmmin
        pxr.xmax = w3d.xmmax
        pxr.ymax = w3d.ymmax
        pxr.zmax = w3d.zmmax
        pxr.x_grid_min=pxr.xmin
        pxr.x_grid_max=pxr.xmax
        pxr.y_grid_min=pxr.ymin
        pxr.y_grid_max=pxr.ymax
        pxr.z_grid_min=pxr.zmin
        pxr.z_grid_max=pxr.zmax

        pxr.x_min_local = self.fields.xmin
        pxr.x_max_local = self.fields.xmax
        pxr.y_min_local = self.fields.ymin
        pxr.y_max_local = self.fields.ymax
        pxr.z_min_local = self.fields.zmin
        pxr.z_max_local = self.fields.zmax
        pxr.x_grid_min_local=pxr.x_min_local
        pxr.x_grid_max_local=pxr.x_max_local
        pxr.y_grid_min_local=pxr.y_min_local
        pxr.y_grid_max_local=pxr.y_max_local
        pxr.z_grid_min_local=pxr.z_min_local
        pxr.z_grid_max_local=pxr.z_max_local
        pxr.zgrid=top.zgrid

        pxr.length_x = pxr.xmax-pxr.xmin
        pxr.length_y = pxr.ymax-pxr.ymin
        pxr.length_z = pxr.zmax-pxr.zmin

        # INIT MPI_DATA FOR PICSAR
        # Init communicator variable in picsar
        if (self.l_debug): print(" Init communicator variable in PXR")
        pxr.mpi_minimal_init_python(top.fsdecomp.mpi_comm)

        # allocate grid quantities
        if (self.l_debug): print(" Allocate grid quantities in PXR")
        pxr.allocate_grid_quantities()
        pxr.compute_simulation_axis()

        # set time step
        pxr.dt = top.dt

        # --- Resolution
        if (self.l_debug): print(" Setup resolution and related variables in PXR")
        pxr.dx = self.dx
        pxr.dy = self.dy
        pxr.dz = self.dz
        pxr.dxi = 1./self.dx
        pxr.dyi = 1./self.dy
        pxr.dzi = 1./self.dz
        pxr.invvol = pxr.dxi*pxr.dyi*pxr.dzi
        pxr.dts2dx = 0.5*pxr.dt*pxr.dxi
        pxr.dts2dy = 0.5*pxr.dt*pxr.dyi
        pxr.dts2dz = 0.5*pxr.dt*pxr.dzi
        pxr.clightsq = 1.0/pxr.clight**2

        # --- Maxwell solver
        pxr.norderx = self.norderx
        pxr.nordery = self.nordery
        pxr.norderz = self.norderz

        pxr.xcoeffs = self.fields.xcoefs
        pxr.ycoeffs = self.fields.ycoefs
        pxr.zcoeffs = self.fields.zcoefs

        # Set coefficient for Maxwell solver
        if (self.l_debug): print(" Set coefficient for Maxwell solver")
        pxr.alphax = em3d.alphax
        pxr.alphay = em3d.alphay
        pxr.alphaz = em3d.alphaz
        pxr.betaxy = em3d.betaxy
        pxr.betayx = em3d.betayx
        pxr.betaxz = em3d.betaxz
        pxr.betazx = em3d.betazx
        pxr.betayz = em3d.betayz
        pxr.betazy = em3d.betazy
        pxr.gammax = em3d.gammax
        pxr.gammay = em3d.gammay
        pxr.gammaz = em3d.gammaz
        pxr.deltaz = em3d.deltaz

        pxr.ex = self.fields.Ex
        pxr.ey = self.fields.Ey
        pxr.ez = self.fields.Ez
        pxr.bx = self.fields.Bx
        pxr.by = self.fields.By
        pxr.bz = self.fields.Bz
        pxr.jx = self.fields.Jx
        pxr.jy = self.fields.Jy
        pxr.jz = self.fields.Jz

        pxr.l_nodalgrid = self.l_nodalgrid

        pxr.nxs = 0
        pxr.nys = 0
        pxr.nzs = 0

        # Current deposition
        pxr.nox = top.depos_order[0][0]
        pxr.noy = top.depos_order[1][0]
        pxr.noz = top.depos_order[2][0]

        if (self.l_debug): print(" Set up algorithms in PXR")

        # Charge deposition algorithm
        pxr.rhodepo=self.rhodepo
        # Current deposition algorithm
        pxr.currdepo=self.currdepo
        # Tye of MPI communication for the current
        pxr.mpicom_curr=self.mpicom_curr
        # Field gathering method
        pxr.fieldgathe=self.fieldgathe
        # Particle communication
        pxr.partcom=self.partcom
        # Field gathering and PArticle pusher separated
        pxr.fg_p_pp_separated=self.fg_p_pp_separated
        # Particle pusher type
        pxr.particle_pusher = top.pgroup.lebcancel_pusher
        # lvec size for the current deposition
        pxr.lvec_curr_depo = self.lvec_curr_depo
        # lvec size for the charge deposition
        pxr.lvec_charge_depo = self.lvec_charge_depo
        # lvec size for the field gathering
        if (self.lvec_fieldgathe==0):
          if ((pxr.nox==3)and(pxr.noy==3)and(pxr.noz==3)):
            pxr.lvec_fieldgathe = 64
          else:
            pxr.lvec_fieldgathe = 512
        else:
          pxr.lvec_fieldgathe = self.lvec_fieldgathe


        #Type of field gathering
        pxr.l4symtry=w3d.l4symtry
        pxr.l_lower_order_in_v = self.l_lower_order_in_v

        # --- Tiling parameters
        pxr.ntilex = self.ntilex
        pxr.ntiley = self.ntiley
        pxr.ntilez = self.ntilez

        # --- Sorting parameters
        if (self.l_debug): print(" Setup sorting parameters in PXR")
        pxr.sorting_activated = self.sorting.activated
        pxr.sorting_dx = self.sorting.dx*pxr.dx
        pxr.sorting_dy = self.sorting.dy*pxr.dy
        pxr.sorting_dz = self.sorting.dz*pxr.dz
        pxr.sorting_shiftx = self.sorting.xshift*pxr.dx
        pxr.sorting_shifty = self.sorting.yshift*pxr.dy
        pxr.sorting_shiftz = self.sorting.zshift*pxr.dz
        pxr.sorting_verbose = self.sorting.verbose

        # --- time statistics
        self.time_stat_loc_array = zeros([20])

        # --- species section
        pxr.nspecies_max=top.pgroup.ns

        # --- allocates array of species
        if (self.l_debug): print(" Allocates array of species")
        pxr.init_species_section()

        for i,s in enumerate(self.listofallspecies):
            # Check for sorting
            if (i >= len(self.sorting.periods)):
              self.sorting.periods.append(0)
              self.sorting.starts.append(0)
            # initialize species in pxr
            pxr.set_particle_species_properties(i+1,s.name,s.mass,s.charge,0, \
                                                0.,0.,0.,0.,0.,0., \
                                                0.,0.,0.,0.,0.,0., \
                                                self.sorting.periods[i],self.sorting.starts[i])
            pxr.nspecies+=1

        pxr.set_tile_split()
        pxr.init_tile_arrays()

        for i,s in enumerate(self.listofallspecies):
            pxr.py_add_particles_to_species(i+1, s.nps,
                                            s.getx(bcast=0,gather=0),
                                            s.gety(bcast=0,gather=0),
                                            s.getz(bcast=0,gather=0),
                                            s.getux(bcast=0,gather=0),
                                            s.getuy(bcast=0,gather=0),
                                            s.getuz(bcast=0,gather=0),
                                            s.getgaminv(bcast=0,gather=0),
                                            s.getweights(bcast=0,gather=0)  )

        top.pgroup.npmax=0
        top.pgroup.ns=1
        top.pgroup.nps=0
        top.pgroup.gchange()

        # --- mirror PXR tile structure in Warp with list of pgroups
        if (self.l_debug): print(" Mirror PXR tile structure in Warp with list of pgroups")
        for i,s in enumerate(self.listofallspecies):
            s.pgroups = []
            s.jslist = [0]
            s.sw=1.
            for iz in range(1,self.ntilez+1):
                xygroup=[]
                for iy in range(1,self.ntiley+1):
                    xgroup=[]
                    for ix in range(1,self.ntilex+1):
                        pg = ParticleGroup()
                        xgroup.append(pg)
                        pxr.point_to_tile(i+1, ix, iy, iz)
                        pg.npmax = 0
                        pxr.partnmax
                        pg.ns=1
                        pg.npid=top.npid
                        pg.gchange()
                        pg.sq = s.charge
                        pg.sm = s.mass
                        pg.sw = s.sw
                        pg.npmax = pxr.partnmax
                        pg.nps = pxr.partn
                        pg.ins[0] = 1
                        pg.sid[0]=0
                        pg.xp = pxr.partx
                        pg.yp = pxr.party
                        pg.zp = pxr.partz
                        pg.uxp = pxr.partux
                        pg.uyp = pxr.partuy
                        pg.uzp = pxr.partuz
                        pg.pid = fzeros([pg.npmax,top.npid])
                        pg.pid = pxr.pid
                        pg.gaminv = pxr.partgaminv
                        pg.ex = pxr.partex
                        pg.ey = pxr.partey
                        pg.ez = pxr.partez
                        pg.bx = pxr.partbx
                        pg.by = pxr.partby
                        pg.bz = pxr.partbz
                        pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                    xygroup.append(xgroup)
                s.pgroups.append(xygroup)
            pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))
#        for i,s in enumerate(self.listofallspecies):
#            def ppzx(self,**kw):
#                for pg in self.pgroups:
#                   self._callppfunc(ppzx,pgroup=pg,**kw)

        if (self.l_debug): print("End allocatefieldarraysPXR")

#            s.ppzx = ppzx
    def aliasparticlearrays(self):
        # --- Detect if tile arrays have been reallocated in PXR
        # --- and make proper aliasing in WARP

        isrealloc=zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8'))
        for i,s in enumerate(self.listofallspecies):
            pxr.get_are_tiles_reallocated(i+1, self.ntilex, self.ntiley, self.ntilez,isrealloc)
            ix,iy,iz=where(isrealloc==1)
            for il in range(0,len(ix)):
                pg = s.pgroups[iz[il]][iy[il]][ix[il]]
                pxr.point_to_tile(i+1, ix[il]+1, iy[il]+1, iz[il]+1)
                pg.npmax = 0
                pxr.partnmax
                pg.ns=1
                pg.npid=top.npid
                pg.gchange()
                pg.sq = s.charge
                pg.sm = s.mass
                pg.sw = s.sw
                pg.npmax = pxr.partnmax
                pg.nps = pxr.partn
                pg.ins[0] = 1
                pg.sid[0]=0
                pg.xp = pxr.partx
                pg.yp = pxr.party
                pg.zp = pxr.partz
                pg.uxp = pxr.partux
                pg.uyp = pxr.partuy
                pg.uzp = pxr.partuz
                pg.pid = fzeros([pg.npmax,top.npid])
                pg.pid = pxr.pid
                pg.gaminv = pxr.partgaminv
                pg.ex = pxr.partex
                pg.ey = pxr.partey
                pg.ez = pxr.partez
                pg.bx = pxr.partbx
                pg.by = pxr.partby
                pg.bz = pxr.partbz
            pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))

    def print_nptiles(self,ispecies):
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    print ix,iy,iz,pxr.partn[0], pxr.partnmax
    def print_nptiles_sp0(self):
        s=self.listofallspecies[0]
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(1, ix, iy, iz)
                    print ix,iy,iz,pxr.partn, pxr.partnmax
                    print ix,iy,iz,s.pgroups[iz-1][iy-1][ix-1].nps, s.pgroups[iz-1][iy-1][ix-1].npmax
    def ppzx_ptiles(self,ispecies,ppg,colors=['black','blue','red','green'],msize=2):
        ncolor = len(colors)
        ic=0
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    ppg(pxr.partx[:pxr.partn[0]],pxr.partz[:pxr.partn[0]],color=colors[ic%ncolor],msize=msize)
                    ic+=1

    def ppzx_ptiles_v2(self,ispecies,ppg,**kw):
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    ppg(pxr.partx[:pxr.partn[0]],pxr.partz[:pxr.partn[0]],kwdict=kw)

    def push_e(self,dir=1.):
        """
        Electric field Maxwell solver
        """
        
        tdeb=MPI.Wtime()
        
        dt = dir*top.dt/self.ntsub
        if self.novercycle==1:
            if dir>0.:
                doit=True
            else:
                doit=False
        else:
            if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
                doit=True
            else:
                doit=False
        if doit:
            if self.laser_mode==1:
                self.add_laser(self.fields)
                if dir<0.:
                    self.fields.Ex_inz*=-1.
                    self.fields.Ey_inz*=-1.
            if self.l_verbose:print 'push_e',self,dt,top.it,self.icycle

            if self.l_pxr:
                f=self.fields
                l_pushe=False
                tdebcell=MPI.Wtime()
                if self.l_2dxz:
                    if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
                        pxr.pxrpush_em2d_evec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                                  f.Jx,f.Jy,f.Jz,
                                                  clight**2*mu0*dt,
                                                  clight**2*dt/f.dx*f.xcoefs[0],
                                                  clight**2*dt/f.dy*f.ycoefs[0],
                                                  clight**2*dt/f.dz*f.zcoefs[0],
                                                  f.nx,f.ny,f.nz,
                                                  f.nxguard,f.nyguard,f.nzguard,
                                                  0,0,0,f.l_nodalgrid)

                    else:
                        pxr.pxrpush_em2d_evec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                                  f.Jx,f.Jy,f.Jz,
                                                  clight**2*mu0*dt,
                                                  clight**2*dt/f.dx*f.xcoefs,
                                                  clight**2*dt/f.dy*f.ycoefs,
                                                  clight**2*dt/f.dz*f.zcoefs,
                                                  f.nx,f.ny,f.nz,
                                                  f.norderx,f.nordery,f.norderz,
                                                  f.nxguard,f.nyguard,f.nzguard,
                                                  0,0,0,f.l_nodalgrid)
                else:
                    if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
                        pxr.pxrpush_em3d_evec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              f.Jx,f.Jy,f.Jz,
                                              clight**2*mu0*dt,
                                              clight**2*dt/f.dx*f.xcoefs[0],
                                              clight**2*dt/f.dy*f.ycoefs[0],
                                              clight**2*dt/f.dz*f.zcoefs[0],
                                              f.nx,f.ny,f.nz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
                    else:
                        pxr.pxrpush_em3d_evec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              f.Jx,f.Jy,f.Jz,
                                              clight**2*mu0*dt,
                                              clight**2*dt/f.dx*f.xcoefs,
                                              clight**2*dt/f.dy*f.ycoefs,
                                              clight**2*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
                tendcell=MPI.Wtime()
                pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
            else:
                l_pushe=True
            push_em3d_eef(self.block,dt,0,self.l_pushf,self.l_pushpot,l_pushe)
            if self.laser_mode==1 and dir<0.:
                self.fields.Ex_inz*=-1.
                self.fields.Ey_inz*=-1.
        if self.refinement is not None:
            self.__class__.__bases__[1].push_e(self.field_coarse,dir)
            
        tend=MPI.Wtime()
        self.time_stat_loc_array[7] += (tend-tdeb)

    def push_b_part_1(self,dir=1.):
      """
      Magnetic field solver
      """
    
      tdeb=MPI.Wtime()
      
      dt = dir*top.dt/self.ntsub
      if self.novercycle==1:
        if dir>0.:
          doit=True
        else:
          doit=False
      else:
        if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
          doit=True
        else:
          doit=False
      if doit:
        if self.l_verbose:print 'push_b part 1',self,dt,top.it,self.icycle,dir
        if self.l_pxr:
          tdebcell=MPI.Wtime()
          f=self.fields
          l_pushb=False
          if self.l_2dxz:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                pxr.pxrpush_em2d_bvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx*f.xcoefs[0],
                            0.5*dt/f.dy*f.ycoefs[0],
                            0.5*dt/f.dz*f.zcoefs[0],
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            0,0,0,f.l_nodalgrid)
              elif (f.stencil==1): # Karkainnen solver
                pxr.pxr_push_em3d_kyeebvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx,
                            0.5*dt/f.dy,
                            0.5*dt/f.dz,
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            f.l_2dxz)
            else: #nth order solver >  2
              pxr.pxrpush_em2d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          else:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                pxr.pxrpush_em3d_bvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx*f.xcoefs[0],
                            0.5*dt/f.dy*f.ycoefs[0],
                            0.5*dt/f.dz*f.zcoefs[0],
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            0,0,0,f.l_nodalgrid)
              elif (f.stencil==1): # Karkainnen solver
                pxr.pxr_push_em3d_kyeebvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx,
                            0.5*dt/f.dy,
                            0.5*dt/f.dz,
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            f.l_2dxz)
            else: #nth order solver >  2
              pxr.pxrpush_em3d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
        else:
          l_pushb=True
        push_em3d_bf(self.block,dt,1,self.l_pushf,self.l_pushpot,l_pushb)
      if self.refinement is not None:
        self.__class__.__bases__[1].push_b_part_1(self.field_coarse,dir)

      tend=MPI.Wtime()
      self.time_stat_loc_array[5] += (tend-tdeb)

    def push_b_part_2(self):
      """
      Magnetic field solver
      """

      tdeb=MPI.Wtime()
    
      if top.efetch[0] != 4 and (self.refinement is None):self.node2yee3d()
      dt = top.dt/self.ntsub
      if self.ntsub<1.:
        self.novercycle = nint(1./self.ntsub)
        self.icycle = (top.it-1)%self.novercycle
      else:
        self.novercycle = 1
        self.icycle = 0
      if self.icycle==0:
        if self.l_verbose:print 'push_b part 2',self,dt,top.it,self.icycle
        if self.l_pxr:
          f=self.fields
          l_pushb=False
          tdebcell=MPI.Wtime()
          if self.l_2dxz:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                pxr.pxrpush_em2d_bvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx*f.xcoefs[0],
                            0.5*dt/f.dy*f.ycoefs[0],
                            0.5*dt/f.dz*f.zcoefs[0],
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            0,0,0,f.l_nodalgrid)
              elif (f.stencil==1): # Karkainnen solver
                pxr.pxr_push_em3d_kyeebvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx,
                            0.5*dt/f.dy,
                            0.5*dt/f.dz,
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            f.l_2dxz)
            else: #nth order solver >  2
              pxr.pxrpush_em2d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          else:
            if (f.norderx==2) & (f.nordery==2) & (f.norderz==2):
              if (f.stencil==0): # Yee solver
                pxr.pxrpush_em3d_bvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx*f.xcoefs[0],
                            0.5*dt/f.dy*f.ycoefs[0],
                            0.5*dt/f.dz*f.zcoefs[0],
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            0,0,0,f.l_nodalgrid)
              elif (f.stencil==1): # Karkainnen solver
                pxr.pxr_push_em3d_kyeebvec(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                            0.5*dt/f.dx,
                            0.5*dt/f.dy,
                            0.5*dt/f.dz,
                            f.nx,f.ny,f.nz,
                            f.nxguard,f.nyguard,f.nzguard,
                            f.l_2dxz)
            else:  #nth order solver >  2
              pxr.pxrpush_em3d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                          0.5*dt/f.dx*f.xcoefs,
                          0.5*dt/f.dy*f.ycoefs,
                          0.5*dt/f.dz*f.zcoefs,
                          f.nx,f.ny,f.nz,
                          f.norderx,f.nordery,f.norderz,
                          f.nxguard,f.nyguard,f.nzguard,
                          0,0,0,f.l_nodalgrid)
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
        else:
          l_pushb=True
        push_em3d_bf(self.block,dt,2,self.l_pushf,self.l_pushpot,l_pushb)
      if self.refinement is not None:
        self.__class__.__bases__[1].push_b_part_2(self.field_coarse)

      tend=MPI.Wtime()
      self.time_stat_loc_array[5] += (tend-tdeb)

    def push_spectral_psaotd(self):
        """
        PSAOTD Maxwell solver
        """
        if self.l_pxr:
          tdebcell=MPI.Wtime()

        if top.it%100==0:print 'push PSAOTD',top.it
        if top.efetch[0] != 4 and (self.refinement is None) and not self.l_nodalgrid:self.node2yee3d()

        if self.ntsub==inf:
          self.GPSTDMaxwell.fields['rhoold']=self.fields.Rhoold
          self.fields.Rho=self.fields.Rhoarray[...,0]
          self.GPSTDMaxwell.fields['rhonew']=self.fields.Rho
        else:
          if self.l_pushf:
        #                self.fields.Rho=self.fields.Rhoarray[...,0]
            self.GPSTDMaxwell.fields['rhoold']=self.fields.Rhoold.copy()
            self.GPSTDMaxwell.fields['rhonew']=self.fields.Rho.copy()
            self.GPSTDMaxwell.fields['drho']=self.fields.Rho-self.fields.Rhoold

        self.GPSTDMaxwell.fields['jx']=self.fields.Jx
        self.GPSTDMaxwell.fields['jy']=self.fields.Jy
        self.GPSTDMaxwell.fields['jz']=self.fields.Jz

        self.GPSTDMaxwell.push_fields()

        b=self.block

        # --- sides
        if b.xlbnd==openbc:self.xlPML.push()
        if b.xrbnd==openbc:self.xrPML.push()
        if b.ylbnd==openbc:self.ylPML.push()
        if b.yrbnd==openbc:self.yrPML.push()
        if b.zlbnd==openbc:self.zlPML.push()
        if b.zrbnd==openbc:self.zrPML.push()

        # --- edges
        if(b.xlbnd==openbc and b.ylbnd==openbc):self.xlylPML.push()
        if(b.xrbnd==openbc and b.ylbnd==openbc):self.xrylPML.push()
        if(b.xlbnd==openbc and b.yrbnd==openbc):self.xlyrPML.push()
        if(b.xrbnd==openbc and b.yrbnd==openbc):self.xryrPML.push()
        if(b.xlbnd==openbc and b.zlbnd==openbc):self.xlzlPML.push()
        if(b.xrbnd==openbc and b.zlbnd==openbc):self.xrzlPML.push()
        if(b.xlbnd==openbc and b.zrbnd==openbc):self.xlzrPML.push()
        if(b.xrbnd==openbc and b.zrbnd==openbc):self.xrzrPML.push()
        if(b.ylbnd==openbc and b.zlbnd==openbc):self.ylzlPML.push()
        if(b.yrbnd==openbc and b.zlbnd==openbc):self.yrzlPML.push()
        if(b.ylbnd==openbc and b.zrbnd==openbc):self.ylzrPML.push()
        if(b.yrbnd==openbc and b.zrbnd==openbc):self.yrzrPML.push()

        # --- corners
        if(b.xlbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc):self.xlylzlPML.push()
        if(b.xrbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc):self.xrylzlPML.push()
        if(b.xlbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc):self.xlyrzlPML.push()
        if(b.xrbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc):self.xryrzlPML.push()
        if(b.xlbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc):self.xlylzrPML.push()
        if(b.xrbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc):self.xrylzrPML.push()
        if(b.xlbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc):self.xlyrzrPML.push()
        if(b.xrbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc):self.xryrzrPML.push()

        #    if em.pml_method==2:
        #      self.fields.spectral=0
        #      scale_em3d_bnd_fields(self.block,top.dt,self.l_pushf)
        #      self.fields.spectral=1

        if self.boris_cor:
          self.boris_correction()
        if self.l_pxr:
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
          self.time_stat_loc_array[7] += (tendcell-tdebcell)

    def current_cor_spectral(self):
        """
        Current spectral correction
        """
    
        if self.l_pxr:
            tdebcell=MPI.Wtime()
            
        j=1j      # imaginary number
        emK = self.FSpace
        em = self
        f = self.fields
        ixl,ixu,iyl,iyu,izl,izu = emK.get_ius()

        fields_shape = [ixu-ixl,iyu-iyl,izu-izl]

        if emK.planj_rfftn is None:
          emK.planj_rfftn= emK.create_plan_rfftn(np.asarray(fields_shape))
        if emK.planj_irfftn is None:
          emK.planj_irfftn= emK.create_plan_irfftn(np.asarray(fields_shape))

        self.wrap_periodic_BC([f.Rho,f.Rhoold_local,f.Jx,f.Jy,f.Jz])

        if emK.nx>1:JxF = emK.rfftn(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)
        if emK.ny>1:JyF = emK.rfftn(squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)
        if emK.nz>1:JzF = emK.rfftn(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]),plan=emK.planj_rfftn)

        em.dRhoodtF = emK.rfftn(squeeze((f.Rho-f.Rhoold_local)[ixl:ixu,iyl:iyu,izl:izu]/top.dt),plan=emK.planj_rfftn)

        # --- get longitudinal J
        divJ = 0.
        if emK.nx>1:divJ += emK.kxmn*JxF
        if emK.ny>1:divJ += emK.kymn*JyF
        if emK.nz>1:divJ += emK.kzmn*JzF

        if emK.nx>1:
          Jxl = emK.kxpn*divJ
        if emK.ny>1:
          Jyl = emK.kypn*divJ
        if emK.nz>1:
          Jzl = emK.kzpn*divJ

        # --- get transverse J
        if emK.nx>1:
          Jxt = JxF-Jxl
        if emK.ny>1:
          Jyt = JyF-Jyl
        if emK.nz>1:
          Jzt = JzF-Jzl

        if emK.nx>1:
          Jxl = j*em.dRhoodtF*emK.kxpn/emK.kmag
        if emK.ny>1:
          Jyl = j*em.dRhoodtF*emK.kypn/emK.kmag
        if emK.nz>1:
          Jzl = j*em.dRhoodtF*emK.kzpn/emK.kmag

        if emK.nx>1:
          JxF = Jxt+Jxl
        if emK.ny>1:
          JyF = Jyt+Jyl
        if emK.nz>1:
          JzF = Jzt+Jzl


        if emK.nx>1:
          Jx = emK.irfftn(JxF, np.asarray(np.shape(squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn, field_out=squeeze(f.Jx[ixl:ixu,iyl:iyu,izl:izu]))
          Jx.resize(fields_shape)
          f.Jx[ixl:ixu,iyl:iyu,izl:izu] = Jx.real
        if emK.ny>1:
          Jy = emK.irfftn(JyF, np.asarray(np.shape(squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn, field_out=squeeze(f.Jy[ixl:ixu,iyl:iyu,izl:izu]))
          Jy.resize(fields_shape)
          f.Jy[ixl:ixu,iyl:iyu,izl:izu] = Jy.real
        if emK.nz>1:
          Jz = emK.irfftn(JzF, np.asarray(np.shape(squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]))), plan=emK.planj_irfftn, field_out=squeeze(f.Jz[ixl:ixu,iyl:iyu,izl:izu]))
          Jz.resize(fields_shape)
          f.Jz[ixl:ixu,iyl:iyu,izl:izu] = Jz.real

        # Time statistics
        if self.l_pxr:
          tendcell=MPI.Wtime()
          pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
          self.time_stat_loc_array[16] += (tendcell-tdebcell)
          

    def exchange_e(self,dir=1.):
        """
        Electric field boundary conditions
        """
        
        t0 = MPI.Wtime()
        
        if self.novercycle==1:
            if dir>0.:
                doit=True
            else:
                doit=False
        else:
            if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
                doit=True
            else:
                doit=False
        if doit:
            if 0:#self.l_pxr:
                print 'exchange e pxr'
                pxr.efield_bcs()
            else:
                em3d_exchange_e(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_e(self.field_coarse)
            
        t1 = MPI.Wtime()
        self.time_stat_loc_array[8] += (t1-t0)

    def exchange_b(self,dir=1.):
        """
        Magnetic field boundary conditions
        """
        
        t0 = MPI.Wtime()
        
        if self.novercycle==1:
            if dir>0.:
                doit=True
            else:
                doit=False
        else:
            if self.icycle==0 or (self.icycle==self.novercycle-1 and dir>0.):
                doit=True
            else:
                doit=False
        if doit:
            if self.l_verbose:print 'exchange_b',self,top.it,self.icycle
            if 0:#self.l_pxr:
                print 'exchange b pxr'
                pxr.bfield_bcs()
            else:
                em3d_exchange_b(self.block)
        if self.refinement is not None:
            self.__class__.__bases__[1].exchange_b(self.field_coarse,dir)

        t1 = MPI.Wtime()
        self.time_stat_loc_array[6] += (t1-t0)


    def step(self,n=1,freq_print=10,lallspecl=0):
      """
      This function performs a range of Particle-In-Cell iterations
      
      Inputs:
      - n: number of iterations
      - freq_print: print frequency
      """

      if (self.l_debug): print("Call step")

      stdout_stat=10
      tdeb=MPI.Wtime()
      for i in range(n):
          if(me==0):
              if top.it%freq_print==0:print 'it = %g time = %g'%(top.it,top.time)
          if lallspecl:
              l_first=l_last=1
          else:
              if i==0:
                  l_first=1
              else:
                  l_first=0
              if i==n-1:
                  l_last=1
              else:
                  l_last=0
                  
          self.onestep(l_first,l_last)
          
          if(l_pxr & (top.it%stdout_stat==0) & (pxr.rank==0)):
              tend=MPI.Wtime()
              mpi_time_per_stat=(tend-tdeb)
              tdeb=MPI.Wtime()
              print("time/stdout_stat (s)",mpi_time_per_stat)

      if (self.l_debug): print("End step")



    def onestep(self,l_first,l_last):
        """
        Perform a single particle-in-cell step
        """

        if (self.l_debug): print("Call onestep")

        # --- Iteration number
        pxr.it = top.it

        # --- call beforestep functions
        if (self.l_debug): print("Call beforestep functions")
        callbeforestepfuncs.callfuncsinlist()
        top.zgrid+=top.vbeamfrm*top.dt
        top.zbeam=top.zgrid
        # --- gather fields from grid to particles
        if (self.l_debug): print("Call Field gathering and particle push")
#        w3d.pgroupfsapi = top.pgroup
#        for js in range(top.pgroup.ns):
#          self.fetcheb(js)
        if l_pxr:
            tdebpart=0.
            tendpart=0.
            tdebfield=0.
            tendfield=0.
            tdebcell=0.
            tendcell=0.
            tdeb=MPI.Wtime()
            pxr.local_time_part=0.
            pxr.local_time_cell=0.
        # --- push
        if l_first:
            if l_pxr:
                # Particle pusher
                tdebpart=MPI.Wtime()
                pxr.pxrpush_particles_part2()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                self.time_stat_loc_array[0] += (tendpart-tdebpart)

                # Particle boundary consitions
                #pxr.particle_bcs_2d()
                tdebpart=MPI.Wtime()
                pxr.particle_bcs()
                tendpart=MPI.Wtime()
                self.time_stat_loc_array[1] += (tendpart-tdebpart)
                
                #for i,s in enumerate(self.listofallspecies):
                #    for pg in s.flatten(s.pgroups):
                #        particleboundaries3d(pg,-1,False)
                #pxr.particle_bcs_tiles()
                self.aliasparticlearrays()
            else:
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.push_velocity_second_half(0,pg)
                        self.push_positions(0,pg)
                        particleboundaries3d(pg,-1,False)
        else:
            if l_pxr:
                # Particle pusher
                tdebpart=MPI.Wtime()
                #pxr.push_particles()
                pxr.field_gathering_plus_particle_pusher()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                self.time_stat_loc_array[0] += (tendpart-tdebpart)

                # Particle boundary conditions
                tdebpart=MPI.Wtime()
                pxr.particle_bcs()
                tendpart=MPI.Wtime()
                self.time_stat_loc_array[1] += (tendpart-tdebpart)
                
                
                #for i,s in enumerate(self.listofallspecies):
                #    for pg in s.flatten(s.pgroups):
                #        particleboundaries3d(pg,-1,False)
                #pxr.particle_bcs_tiles()
                self.aliasparticlearrays()

            else:
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        tendpart=MPI.Wtime()
                        self.push_velocity_full(0,pg)
                        self.push_positions(0,pg)
                        tendpart=MPI.Wtime()
                        self.time_stat_loc_array[0] += (tendpart-tdebpart)

                        # Particle boundary conditions
                        tdebpart=MPI.Wtime()
                        particleboundaries3d(pg,-1,False)
                        tendpart=MPI.Wtime()
                        self.time_stat_loc_array[1] += (tendpart-tdebpart)
                        
        # --- Particle sorting
        if (self.l_debug): print("Call Particle Sorting")
        if l_pxr:
          if ((self.sorting.activated)and(top.it>=0)):
            tdebpart=MPI.Wtime()
            pxr.particle_sorting_sub()
            tendpart=MPI.Wtime()
            self.time_stat_loc_array[10] += (tendpart-tdebpart)

        # --- call beforeloadrho functions
        if (self.l_debug): print("Call beforeloadrho functions")
        beforeloadrho.callfuncsinlist()
        pgroups = []
        for i,s in enumerate(self.listofallspecies):
            pgroups+=s.flatten(s.pgroups)
        self.pgroups = pgroups
#        self.loadsource(pgroups=pgroups)
        #tdebpart=MPI.Wtime()

        # Call user-defined injection routines
        if (self.l_debug): print("Call user-defined injection routines")
        userinjection.callfuncsinlist()
        if (self.l_pxr):
            pxr.zgrid=self.zgrid
        self.loadrho(pgroups=pgroups)
        self.loadj(pgroups=pgroups)
        # Moving window

        #tendpart=MPI.Wtime()
        #pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
#        self.solve2ndhalf()

        #tdebcell=MPI.Wtime()
        # --- dosolve
        # Current deposition + Maxwell

        self.dosolve()

        #tendcell=MPI.Wtime()
        #pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)

        if l_pxr:
            if l_last:
                tdebpart=MPI.Wtime()
                pxr.pxrpush_particles_part1()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                self.time_stat_loc_array[0] += (tendpart-tdebpart)
        else:
            t0=MPI.Wtime()
            for i,s in enumerate(self.listofallspecies):
                for pg in s.flatten(s.pgroups):
                    w3d.pgroupfsapi = pg
                    self.fetcheb(0,pg)
                    if l_last:
                        self.push_velocity_first_half(0,pg)
            t1 = MPI.Wtime()
            self.time_stat_loc_array[0] += (t1-t0)

        # --- update time, time counter
        top.time+=top.dt
        if top.it%top.nhist==0:
#           zmmnt()
           minidiag(top.it,top.time,top.lspecial)
        top.it+=1

        # Load balance every dlb_freq time step
        if (self.l_debug): print("Call Load balance")
        if (l_pxr & (self.dload_balancing & (top.it%self.dlb_freq==0))):
            pxr.mpitime_per_it=pxr.local_time_part+pxr.local_time_cell
            pxr.get_max_time_per_it()
            pxr.get_min_time_per_it()
            ## --- Compute time per part and per cell
            pxr.compute_time_per_part()
            pxr.compute_time_per_cell()
            imbalance=(pxr.max_time_per_it-pxr.min_time_per_it)/pxr.min_time_per_it*100.
            if (imbalance>self.dlb_threshold):
            	if (self.l_2dxz):
                	self.load_balance_2d(str(imbalance)+"%")
                else:
                    self.load_balance_3d(str(imbalance)+"%")
        # Try to Load balance at init
        if ((top.it==self.it_dlb_init) & self.dlb_at_init & self.dload_balancing):
          pxr.mpitime_per_it=pxr.local_time_part+pxr.local_time_cell
          pxr.get_max_time_per_it()
          pxr.get_min_time_per_it()
          ## --- Compute time per part and per cell
          pxr.compute_time_per_part()
          pxr.compute_time_per_cell()
          if (self.l_2dxz):
              self.load_balance_2d('Init')
          else:
            self.load_balance_3d('Init')

        # PXr custom outputs mpi-io
        if (self.l_debug): print("Call PXR custom outputs mpi-io")
        if(l_pxr & self.l_output_grid & (top.it % self.l_output_freq ==0)):
          self.output_pxr(top.it)

        # --- call afterstep functions
        callafterstepfuncs.callfuncsinlist()

    def load_balance_3d(self,imbalance):
        """
        Load balance between MPI domains in 3D
        """
        if (l_pxr):
        
            tdeb = MPIWTIME()
        
            ## --- Compute time per part and per cell
            pxr.compute_time_per_part()
            pxr.compute_time_per_cell()

            ## --- Compute new split along each dimension
            pxr.compute_new_split(pxr.global_time_per_part,pxr.global_time_per_cell,pxr.nx_global,pxr.ny_global,pxr.nz_global,
                              pxr.new_cell_x_min,pxr.new_cell_x_max,pxr.new_cell_y_min,pxr.new_cell_y_max,
                              pxr.new_cell_z_min,pxr.new_cell_z_max,pxr.nprocx,pxr.nprocy,pxr.nprocz)
            isnewsplit=sum(pxr.cell_x_min-pxr.new_cell_x_min)+sum(pxr.cell_x_max-pxr.new_cell_x_max)+ \
                	   sum(pxr.cell_y_min-pxr.new_cell_y_min)+sum(pxr.cell_y_max-pxr.new_cell_y_max)+ \
                	   sum(pxr.cell_z_min-pxr.new_cell_z_min)+sum(pxr.cell_z_max-pxr.new_cell_z_max)
            if (isnewsplit==0):
            	if(pxr.rank==0):
                	print("Optimal load balancing already achieved by current implementation")
            else:
            	if(pxr.rank==0):
                	print("trying to load balance the simulation, imbalance=", imbalance)
                ## --- Compute limits for all procs
                ix1old=np.zeros(pxr.nproc,dtype="i8"); ix2old=np.zeros(pxr.nproc,dtype="i8")
                iy1old=np.zeros(pxr.nproc,dtype="i8"); iy2old=np.zeros(pxr.nproc,dtype="i8")
                iz1old=np.zeros(pxr.nproc,dtype="i8"); iz2old=np.zeros(pxr.nproc,dtype="i8")
                ix1new=np.zeros(pxr.nproc,dtype="i8"); ix2new=np.zeros(pxr.nproc,dtype="i8")
                iy1new=np.zeros(pxr.nproc,dtype="i8"); iy2new=np.zeros(pxr.nproc,dtype="i8")
                iz1new=np.zeros(pxr.nproc,dtype="i8"); iz2new=np.zeros(pxr.nproc,dtype="i8")

                pxr.get_1darray_proclimits(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,
                                        pxr.cell_x_min,pxr.cell_y_min,pxr.cell_z_min,
                                        pxr.cell_x_max,pxr.cell_y_max,pxr.cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,
                                        top.lcomm_cartesian)
                pxr.get_1darray_proclimits(ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,
                                        pxr.new_cell_x_min,pxr.new_cell_y_min,pxr.new_cell_z_min,
                                        pxr.new_cell_x_max,pxr.new_cell_y_max,pxr.new_cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,top.lcomm_cartesian)
                ## --- Compute new sizes for grid arrays
                nx_new=pxr.new_cell_x_max[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+1
                ny_new=pxr.new_cell_y_max[pxr.y_coords]-pxr.new_cell_y_min[pxr.y_coords]+1
                nz_new=pxr.new_cell_z_max[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+1

                ## --- Remap field arrays
                # -- Ex
                ex_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(ex_new,nx_new,ny_new,nz_new,
                                        	pxr.ex,pxr.nx,pxr.ny,pxr.nz,
                                        	pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                        	ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                        	ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.ex=ex_new
                # -- Ey
                ey_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(ey_new,nx_new,ny_new,nz_new,
                                        	pxr.ey,pxr.nx,pxr.ny,pxr.nz,
                                        	pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                        	ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                        	ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.ey=ey_new
                # -- Ez
                ez_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(ez_new,nx_new,ny_new,nz_new,
                                        	pxr.ez,pxr.nx,pxr.ny,pxr.nz,
                                        	pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                        	ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                        	ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.ez=ez_new
                # -- Bx
                bx_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(bx_new,nx_new,ny_new,nz_new,
                                        	pxr.bx,pxr.nx,pxr.ny,pxr.nz,
                                        	pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                        	ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                        	ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.bx=bx_new
                # -- By
                by_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(by_new,nx_new,ny_new,nz_new,
                                        	pxr.by,pxr.nx,pxr.ny,pxr.nz,
                                        	pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                        	ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                        	ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.by=by_new
                # -- Bz
                bz_new=zeros((nx_new+2*pxr.nxguards+1,ny_new+2*pxr.nyguards+1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_3d_field_component(bz_new,nx_new,ny_new,nz_new,
                                        	pxr.bz,pxr.nx,pxr.ny,pxr.nz,
                                        	pxr.nxguards,pxr.nyguards,pxr.nzguards,
                                        	ix1old, ix2old, iy1old, iy2old, iz1old, iz2old,
                                        	ix1new, ix2new, iy1new, iy2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.bz=bz_new
                ## -- Reallocate current arrays
                # Currents are recomputed each iteration so no need to exchange them
                jx_new=zeros((nx_new+2*pxr.nxjguards+1,ny_new+2*pxr.nyjguards+1,nz_new+2*pxr.nzjguards+1),order='F')
                jy_new=zeros((nx_new+2*pxr.nxjguards+1,ny_new+2*pxr.nyjguards+1,nz_new+2*pxr.nzjguards+1),order='F')
                jz_new=zeros((nx_new+2*pxr.nxjguards+1,ny_new+2*pxr.nyjguards+1,nz_new+2*pxr.nzjguards+1),order='F')
                pxr.jx=jx_new
                pxr.jy=jy_new
                pxr.jz=jz_new

                # Update pxr new array dimensions
                pxr.nx=nx_new
                pxr.ny=ny_new
                pxr.nz=nz_new
                pxr.nx_grid=pxr.nx+1
                pxr.ny_grid=pxr.ny+1
                pxr.nz_grid=pxr.nz+1

                # Test if domain has been resized - used for particle remaping
                isnewdom=pxr.cell_x_min[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+pxr.cell_x_max[pxr.x_coords]-pxr.new_cell_x_max[pxr.x_coords]+ \
                pxr.cell_y_min[pxr.y_coords]-pxr.new_cell_y_min[pxr.y_coords]+pxr.cell_y_max[pxr.y_coords]-pxr.new_cell_y_max[pxr.y_coords]+ \
                pxr.cell_z_min[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+pxr.cell_z_max[pxr.z_coords]-pxr.new_cell_z_max[pxr.z_coords]

                # Update new subdomain index arrays
                pxr.cell_x_min=pxr.new_cell_x_min
                pxr.cell_x_max=pxr.new_cell_x_max
                pxr.cell_y_min=pxr.new_cell_y_min
                pxr.cell_y_max=pxr.new_cell_y_max
                pxr.cell_z_min=pxr.new_cell_z_min
                pxr.cell_z_max=pxr.new_cell_z_max
                pxr.nx_global_grid_min = pxr.cell_x_min[pxr.x_coords]
                pxr.nx_global_grid_max = pxr.cell_x_max[pxr.x_coords]+1
                pxr.ny_global_grid_min = pxr.cell_y_min[pxr.y_coords]
                pxr.ny_global_grid_max = pxr.cell_y_max[pxr.y_coords]+1
                pxr.nz_global_grid_min = pxr.cell_z_min[pxr.z_coords]
                pxr.nz_global_grid_max = pxr.cell_z_max[pxr.z_coords]+1


                # Update global simulation axis
                pxr.compute_simulation_axis()

                # Set new min and max for local domain
                pxr.x_min_local = pxr.x_grid_mins[pxr.x_coords]
                pxr.x_max_local = pxr.x_grid_maxs[pxr.x_coords]
                pxr.y_min_local = pxr.y_grid_mins[pxr.y_coords]
                pxr.y_max_local = pxr.y_grid_maxs[pxr.y_coords]
                pxr.z_min_local = pxr.z_grid_mins[pxr.z_coords]
                pxr.z_max_local = pxr.z_grid_maxs[pxr.z_coords]
                pxr.x_grid_min_local=pxr.x_min_local
                pxr.x_grid_max_local=pxr.x_max_local
                pxr.y_grid_min_local=pxr.y_min_local
                pxr.y_grid_max_local=pxr.y_max_local
                pxr.z_grid_min_local=pxr.z_min_local
                pxr.z_grid_max_local=pxr.z_max_local

                ##--- Alias WARP grid arrays on pxr new arrays
                self.nxlocal=pxr.nx
                self.nylocal=pxr.ny
                self.nzlocal=pxr.nz
                self.ymminlocal = pxr.y_min_local
                self.zmminlocal = pxr.z_min_local
                self.fields.xmin = pxr.x_min_local
                self.fields.xmax = pxr.x_max_local
                self.fields.ymin = pxr.y_min_local
                self.fields.ymax = pxr.y_max_local
                self.fields.zmin = pxr.z_min_local
                self.fields.zmax = pxr.z_max_local
                
                # Udpate domain decomposition in WARP
                top.fsdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.fsdecomp.ny=pxr.cell_y_max-pxr.cell_y_min+1
                top.fsdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.fsdecomp.ix=pxr.cell_x_min
                top.fsdecomp.iy=pxr.cell_y_min
                top.fsdecomp.iz=pxr.cell_z_min
                top.fsdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.fsdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.fsdecomp.ymin=pxr.cell_y_min*pxr.dy
                top.fsdecomp.ymax=(pxr.cell_y_max+1)*pxr.dy
                top.fsdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.fsdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz
                top.ppdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.ppdecomp.ny=pxr.cell_y_max-pxr.cell_y_min+1
                top.ppdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.ppdecomp.ix=pxr.cell_x_min
                top.ppdecomp.iy=pxr.cell_y_min
                top.ppdecomp.iz=pxr.cell_z_min
                top.ppdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.ppdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.ppdecomp.ymin=pxr.cell_y_min*pxr.dy
                top.ppdecomp.ymax=(pxr.cell_y_max+1)*pxr.dy
                top.ppdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.ppdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz

                # Reallocate warp arrays
                self.allocatefieldarrays()
                # Alias newly allocated arrays on WARP structure
                self.fields.Ex=pxr.ex
                self.fields.Ey=pxr.ey
                self.fields.Ez=pxr.ez
                self.fields.Bx=pxr.bx
                self.fields.By=pxr.by
                self.fields.Bz=pxr.bz
                self.fields.Exp=pxr.ex
                self.fields.Eyp=pxr.ey
                self.fields.Ezp=pxr.ez
                self.fields.Bxp=pxr.bx
                self.fields.Byp=pxr.by
                self.fields.Bzp=pxr.bz
                self.fields.Jx=pxr.jx
                self.fields.Jy=pxr.jy
                self.fields.Jz=pxr.jz


                em3d_exchange_e(self.block)
                em3d_exchange_b(self.block)

                # If domain has been resized, do a new tile split and exchange particles
                if 1:#((isnewdom != 0)):
                   # Now exchanging particles
                    pxr.create_new_tile_split()

                    self.ntilex = pxr.ntilex
                    self.ntiley = pxr.ntiley
                    self.ntilez = pxr.ntilez

                  # Alias PXR tiles to WARP pgroups
                    for i,s in enumerate(self.listofallspecies):
                      s.pgroups = []
                      s.jslist = [0]
                      for iz in range(1,self.ntilez+1):
                        xygroup=[]
                        for iy in range(1,self.ntiley+1):
                          xgroup=[]
                          for ix in range(1,self.ntilex+1):
                            pg = ParticleGroup()
                            xgroup.append(pg)
                            pxr.point_to_tile(i+1, ix, iy, iz)
                            pg.npmax = 0
                            pxr.partnmax
                            pg.ns=1
                            pg.npid=top.npid
                            pg.gchange()
                            pg.sq = s.charge
                            pg.sm = s.mass
                            pg.sw = s.sw
                            pg.npmax = pxr.partnmax
                            pg.nps = pxr.partn
                            pg.ins[0] = 1
                            pg.sid[0]=0
                            pg.xp = pxr.partx
                            pg.yp = pxr.party
                            pg.zp = pxr.partz
                            pg.uxp = pxr.partux
                            pg.uyp = pxr.partuy
                            pg.uzp = pxr.partuz
                            pg.pid = fzeros([pg.npmax,top.npid])
                            pg.pid = pxr.pid
                            pg.gaminv = pxr.partgaminv
                            pg.ex = pxr.partex
                            pg.ey = pxr.partey
                            pg.ez = pxr.partez
                            pg.bx = pxr.partbx
                            pg.by = pxr.partby
                            pg.bz = pxr.partbz
                            pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                          xygroup.append(xgroup)
                        s.pgroups.append(xygroup)
                      pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))
#                pxr.particle_bcs_mpi_blocking()
                pxr.remap_particles(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,
                            ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,
                            pxr.cell_x_min,pxr.cell_x_max,pxr.cell_y_min,pxr.cell_y_max,
                            pxr.cell_z_min,pxr.cell_z_max,
                            pxr.rank, pxr.nproc, pxr.nprocx, pxr.nprocy,pxr.nprocz,top.lcomm_cartesian)

            # Time statistics
            tend=MPI.Wtime()
            self.time_stat_loc_array[15] += (tend-tdeb)
                
    def load_balance_2d(self,imbalance):
        if (l_pxr):
            ## --- Compute time per part and per cell
            pxr.compute_time_per_part()
            pxr.compute_time_per_cell()

            ## --- Compute new split along each dimension
            pxr.compute_new_split_2d(pxr.global_time_per_part,pxr.global_time_per_cell,pxr.nx_global,pxr.nz_global,
                              pxr.new_cell_x_min,pxr.new_cell_x_max,
                              pxr.new_cell_z_min,pxr.new_cell_z_max,pxr.nprocx,pxr.nprocz)
            isnewsplit=sum(pxr.cell_x_min-pxr.new_cell_x_min)+sum(pxr.cell_x_max-pxr.new_cell_x_max)+ \
                	   sum(pxr.cell_z_min-pxr.new_cell_z_min)+sum(pxr.cell_z_max-pxr.new_cell_z_max)
            if (isnewsplit==0):
                if(pxr.rank==0):
                  print("Optimal load balancing already achieved by current implementation")
            else:
                if(pxr.rank==0):
                  print("trying to load balance the simulation, imbalance=", imbalance)

                ## --- Compute limits for all procs
                ix1old=np.zeros(pxr.nproc,dtype="i8"); ix2old=np.zeros(pxr.nproc,dtype="i8")
                iy1old=np.zeros(pxr.nproc,dtype="i8"); iy2old=np.zeros(pxr.nproc,dtype="i8")
                iz1old=np.zeros(pxr.nproc,dtype="i8"); iz2old=np.zeros(pxr.nproc,dtype="i8")
                ix1new=np.zeros(pxr.nproc,dtype="i8"); ix2new=np.zeros(pxr.nproc,dtype="i8")
                iy1new=np.zeros(pxr.nproc,dtype="i8"); iy2new=np.zeros(pxr.nproc,dtype="i8")
                iz1new=np.zeros(pxr.nproc,dtype="i8"); iz2new=np.zeros(pxr.nproc,dtype="i8")

                pxr.get_1darray_proclimits(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,
                                        pxr.cell_x_min,pxr.cell_y_min,pxr.cell_z_min,
                                        pxr.cell_x_max,pxr.cell_y_max,pxr.cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,
                                        top.lcomm_cartesian)
                pxr.get_1darray_proclimits(ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,
                                        pxr.new_cell_x_min,pxr.new_cell_y_min,pxr.new_cell_z_min,
                                        pxr.new_cell_x_max,pxr.new_cell_y_max,pxr.new_cell_z_max,
                                        pxr.nprocx, pxr.nprocy, pxr.nprocz, pxr.nproc,top.lcomm_cartesian)
                ## --- Compute new sizes for grid arrays
                nx_new=pxr.new_cell_x_max[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+1
                nz_new=pxr.new_cell_z_max[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+1

                ## --- Remap field arrays
                # -- Ex
                ex_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(ex_new,nx_new,nz_new,
                                        	pxr.ex,pxr.nx,pxr.nz,
                                        	pxr.nxguards,pxr.nzguards,
                                        	ix1old, ix2old, iz1old, iz2old,
                                        	ix1new, ix2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.ex=ex_new
                # -- Ey
                ey_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(ey_new,nx_new,nz_new,
                                        	pxr.ey,pxr.nx,pxr.nz,
                                        	pxr.nxguards,pxr.nzguards,
                                        	ix1old, ix2old, iz1old, iz2old,
                                        	ix1new, ix2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.ey=ey_new
                # -- Ez
                ez_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(ez_new,nx_new,nz_new,
                                        	pxr.ez,pxr.nx,pxr.nz,
                                        	pxr.nxguards,pxr.nzguards,
                                        	ix1old, ix2old, iz1old, iz2old,
                                        	ix1new, ix2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.ez=ez_new
                # -- Bx
                bx_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(bx_new,nx_new,nz_new,
                                        	pxr.bx,pxr.nx,pxr.nz,
                                        	pxr.nxguards,pxr.nzguards,
                                        	ix1old, ix2old, iz1old, iz2old,
                                        	ix1new, ix2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.bx=bx_new
                # -- By
                by_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(by_new,nx_new,nz_new,
                                        	pxr.by,pxr.nx,pxr.nz,
                                        	pxr.nxguards,pxr.nzguards,
                                        	ix1old, ix2old, iz1old, iz2old,
                                        	ix1new, ix2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.by=by_new
                # -- Bz
                bz_new=zeros((nx_new+2*pxr.nxguards+1,1,nz_new+2*pxr.nzguards+1),order='F')
                pxr.mpi_remap_2d_field_component(bz_new,nx_new,nz_new,
                                        	pxr.bz,pxr.nx,pxr.nz,
                                        	pxr.nxguards,pxr.nzguards,
                                        	ix1old, ix2old, iz1old, iz2old,
                                        	ix1new, ix2new, iz1new, iz2new,
                                        	pxr.rank, pxr.nproc)
                pxr.bz=bz_new
                ## -- Reallocate current arrays
                # Currents are recomputed each iteration so no need to exchange them
                jx_new=zeros((nx_new+2*pxr.nxjguards+1,1,nz_new+2*pxr.nzjguards+1),order='F')
                jy_new=zeros((nx_new+2*pxr.nxjguards+1,1,nz_new+2*pxr.nzjguards+1),order='F')
                jz_new=zeros((nx_new+2*pxr.nxjguards+1,1,nz_new+2*pxr.nzjguards+1),order='F')
                pxr.jx=jx_new
                pxr.jy=jy_new
                pxr.jz=jz_new

                # Update pxr new array dimensions
                pxr.nx=nx_new
                pxr.nz=nz_new
                pxr.nx_grid=pxr.nx+1
                pxr.nz_grid=pxr.nz+1

                # Test if domain has been resized - used for particle remaping
                isnewdom=pxr.cell_x_min[pxr.x_coords]-pxr.new_cell_x_min[pxr.x_coords]+pxr.cell_x_max[pxr.x_coords]-pxr.new_cell_x_max[pxr.x_coords]+ \
                pxr.cell_z_min[pxr.z_coords]-pxr.new_cell_z_min[pxr.z_coords]+pxr.cell_z_max[pxr.z_coords]-pxr.new_cell_z_max[pxr.z_coords]

                # Update new subdomain index arrays
                pxr.cell_x_min=pxr.new_cell_x_min
                pxr.cell_x_max=pxr.new_cell_x_max
                pxr.cell_z_min=pxr.new_cell_z_min
                pxr.cell_z_max=pxr.new_cell_z_max
                pxr.nx_global_grid_min = pxr.cell_x_min[pxr.x_coords]
                pxr.nx_global_grid_max = pxr.cell_x_max[pxr.x_coords]+1
                pxr.nz_global_grid_min = pxr.cell_z_min[pxr.z_coords]
                pxr.nz_global_grid_max = pxr.cell_z_max[pxr.z_coords]+1


                # Update global simulation axis
                pxr.compute_simulation_axis()

                # Set new min and max for local domain
                pxr.x_min_local = pxr.x_grid_mins[pxr.x_coords]
                pxr.x_max_local = pxr.x_grid_maxs[pxr.x_coords]
                pxr.z_min_local = pxr.z_grid_mins[pxr.z_coords]
                pxr.z_max_local = pxr.z_grid_maxs[pxr.z_coords]
                pxr.x_grid_min_local=pxr.x_min_local
                pxr.x_grid_max_local=pxr.x_max_local
                pxr.z_grid_min_local=pxr.z_min_local
                pxr.z_grid_max_local=pxr.z_max_local

                ##--- Alias WARP grid arrays on pxr new arrays
                self.nxlocal=pxr.nx
                self.nzlocal=pxr.nz
                self.xmminlocal = pxr.x_min_local
                self.zmminlocal = pxr.z_min_local
                self.fields.xmin = pxr.x_min_local
                self.fields.xmax = pxr.x_max_local
                self.fields.zmin = pxr.z_min_local
                self.fields.zmax = pxr.z_max_local
                
                # Udpate domain decomposition in WARP
                top.fsdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.fsdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.fsdecomp.ix=pxr.cell_x_min
                top.fsdecomp.iz=pxr.cell_z_min
                top.fsdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.fsdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.fsdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.fsdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz
                top.ppdecomp.nx=pxr.cell_x_max-pxr.cell_x_min+1
                top.ppdecomp.nz=pxr.cell_z_max-pxr.cell_z_min+1
                top.ppdecomp.ix=pxr.cell_x_min
                top.ppdecomp.iz=pxr.cell_z_min
                top.ppdecomp.xmin=pxr.cell_x_min*pxr.dx
                top.ppdecomp.xmax=(pxr.cell_x_max+1)*pxr.dx
                top.ppdecomp.zmin=pxr.cell_z_min*pxr.dz
                top.ppdecomp.zmax=(pxr.cell_z_max+1)*pxr.dz

                # Reallocate warp arrays
                self.allocatefieldarrays()
                if (self.spectral==1):
                  self.allocatefieldarraysFFT()
                # Alias newly allocated arrays on WARP structure
                self.fields.Ex=pxr.ex
                self.fields.Ey=pxr.ey
                self.fields.Ez=pxr.ez
                self.fields.Bx=pxr.bx
                self.fields.By=pxr.by
                self.fields.Bz=pxr.bz
                self.fields.Exp=pxr.ex
                self.fields.Eyp=pxr.ey
                self.fields.Ezp=pxr.ez
                self.fields.Bxp=pxr.bx
                self.fields.Byp=pxr.by
                self.fields.Bzp=pxr.bz
                self.fields.Jx=pxr.jx
                self.fields.Jy=pxr.jy
                self.fields.Jz=pxr.jz


                em3d_exchange_e(self.block)
                em3d_exchange_b(self.block)
                #If domain has been resized, do a new tile split and exchange particles
                if 1:#((isnewdom != 0)):
                # Now exchanging particles
                    pxr.create_new_tile_split()
                    pxr.remap_particles_2d(ix1old,ix2old,iz1old,iz2old,
                                            ix1new,ix2new,iz1new,iz2new,
                                            pxr.cell_x_min,pxr.cell_x_max,
                                            pxr.cell_z_min,pxr.cell_z_max,
                                            pxr.rank, pxr.nproc, pxr.nprocx,pxr.nprocz,top.lcomm_cartesian)
                    self.ntilex = pxr.ntilex
                    self.ntilez = pxr.ntilez

					# Alias PXR tiles to WARP pgroups
                    for i,s in enumerate(self.listofallspecies):
                        s.pgroups = []
                        s.jslist = [0]
                        for iz in range(1,self.ntilez+1):
                          xygroup=[]
                          for iy in range(1,self.ntiley+1):
                            xgroup=[]
                            for ix in range(1,self.ntilex+1):
                              pg = ParticleGroup()
                              xgroup.append(pg)
                              pxr.point_to_tile(i+1, ix, iy, iz)
                              pg.npmax = 0
                              pxr.partnmax
                              pg.ns=1
                              pg.npid=top.npid
                              pg.gchange()
                              pg.sq = s.charge
                              pg.sm = s.mass
                              pg.sw = s.sw
                              pg.npmax = pxr.partnmax
                              pg.nps = pxr.partn
                              pg.ins[0] = 1
                              pg.sid[0]=0
                              pg.xp = pxr.partx
                              pg.yp = pxr.party
                              pg.zp = pxr.partz
                              pg.uxp = pxr.partux
                              pg.uyp = pxr.partuy
                              pg.uzp = pxr.partuz
                              pg.pid = fzeros([pg.npmax,top.npid])
                              pg.pid = pxr.pid
                              pg.gaminv = pxr.partgaminv
                              pg.ex = pxr.partex
                              pg.ey = pxr.partey
                              pg.ez = pxr.partez
                              pg.bx = pxr.partbx
                              pg.by = pxr.partby
                              pg.bz = pxr.partbz
                              pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                            xygroup.append(xgroup)
                          s.pgroups.append(xygroup)
                        pxr.set_are_tiles_reallocated(i+1, self.ntilex,self.ntiley,self.ntilez,zeros((self.ntilex,self.ntiley,self.ntilez),dtype=dtype('i8')))


    def output_pxr(self,iter):
      pxr.py_mpi_output_grid_quantity('ez',pxr.ez,pxr.nx,pxr.ny,pxr.nz,pxr.nxguards,pxr.nyguards,pxr.nzguards,iter)


    def fetcheb(self,js,pg=None):
        if self.l_verbose:print me,'enter fetcheb'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        w3d.ipminfsapi=pg.ins[js]
        w3d.npfsapi=pg.nps[js]
        pg.ex[il:iu]=0.
        pg.ey[il:iu]=0.
        pg.ez[il:iu]=0.
        pg.bx[il:iu]=0.
        pg.by[il:iu]=0.
        pg.bz[il:iu]=0.
        self.fetche()
        self.fetchb()

    def push_velocity_full(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_velocity_full'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        if pg.lebcancel_pusher:
          ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                            pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                            pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                            pg.sq[js],pg.sm[js],top.dt,0)
        else:
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)
          # --- push velocity from magnetic field
          bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                      pg.sq[js],pg.sm[js],top.dt, top.ibpush)
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)

        if self.l_verbose:print me,'exit push_ions_velocity_first_half'

    def push_velocity_first_half(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_velocity_first_half'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        if pg.lebcancel_pusher:
          ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                            pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                            pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                            pg.sq[js],pg.sm[js],top.dt,1)
        else:
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)
          # --- push velocity from magnetic field
          bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                      pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)

        if self.l_verbose:print me,'exit push_ions_velocity_first_half'

    def push_velocity_second_half(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_velocity_second_half'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        if pg.lebcancel_pusher:
          ebcancelpush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                            pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                            pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                            pg.sq[js],pg.sm[js],top.dt,2)
        else:
          # --- push velocity from magnetic field
          bpush3d (np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                      pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                      pg.sq[js],pg.sm[js],0.5*top.dt, top.ibpush)
          # --- push velocity from electric field (half step)
          epush3d(np,pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                     pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                     pg.sq[js],pg.sm[js],0.5*top.dt)
        # --- update gamma
        self.set_gamma(js,pg)

        if self.l_verbose:print me,'exit push_ions_velocity_second_half'

    def set_gamma(self,js,pg=None):
        if self.l_verbose:print me,'enter set_gamma'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        # --- update gamma
        gammaadv(np,pg.gaminv[il:iu],pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                 top.gamadv,top.lrelativ)

        if self.l_verbose:print me,'exit push_ions_velocity_second_half'

    def push_positions(self,js,pg=None):
        if self.l_verbose:print me,'enter push_ions_positions'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        xpush3d(np,pg.xp[il:iu],pg.yp[il:iu],pg.zp[il:iu],
                       pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],
                       pg.gaminv[il:iu],top.dt)

        if self.l_verbose:print me,'exit push_ions_positions'

    def loadsource(self,lzero=None,lfinalize_rho=None,pgroups=None,**kw):
        '''
        Current and charge deposition, uses particles from top directly.
        
        Inputs:
           - lzero
           - lfinalize_rho
           - pgroups
        '''
        
        # --- Note that the grid location is advanced even if no field solve
        # --- is being done.
        self.advancezgrid()
        # --- If ldosolve is false, then skip the gather of rho, unless
        # --- lzero is also false, in which case the solver is assumed to
        # --- be gathering the source (for example during an EGUN iteration).
        if not self.ldosolve and lzero: return
        if lzero is None: lzero = w3d.lzerorhofsapi
        if lfinalize_rho is None: lfinalize_rho = w3d.lfinalizerhofsapi

        self.setparticledomains()
        self.allocatedataarrays()
        if lzero: self.zerosourcep()

        if pgroups is None: pgroups = [top.pgroup]

        if  l_pxr:
            # --- PICSAR current deposition
            # --- js = 0
             f=self.fields

             for pgroup in pgroups:
                if w3d.js1fsapi >= 0: js1 = w3d.js1fsapi
                else:                 js1 = 0
                if w3d.js2fsapi >= 0: js2 = w3d.js2fsapi+1
                else:                 js2 = pgroup.ns

                jslist = kw.get('jslist',None)
                if jslist is None: jslist = range(js1,js2)

                for js in jslist:
                    n = pgroup.nps[js]
                    if n == 0: continue
                    if pgroup.ldts[js]:
                        indts = top.ndtstorho[pgroup.ndts[js]-1]
                        iselfb = pgroup.iselfb[js]
                        self.setsourcepforparticles(0,indts,iselfb)

                        if self.debug:
                            i1 = pgroup.ins[js]-1
                            i2 = pgroup.ins[js]+pgroup.nps[js]-1
                            if self.nxlocal > 0:
                                x = pgroup.xp[i1:i2]
                                if self.l4symtry: x = abs(x)
                                if self.solvergeom == w3d.RZgeom:
                                    y = pgroup.yp[i1:i2]
                                    x = sqrt(x**2 + y**2)
                                assert x.min() >= self.xmminp,\
                                       "Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,x.min())
                                assert x.max() < self.xmmaxp,\
                                       "Particles in species %d have x above the grid when depositing the source, max x = %e"%(js,x.max())
                            if self.nylocal > 0:
                                y = pgroup.yp[i1:i2]
                                if self.l4symtry or self.l2symtry: y = abs(y)
                                assert y.min() >= self.ymminp,\
                                       "Particles in species %d have y below the grid when depositing the source, min y = %e"%(js,y.min())
                                assert y.max() < self.ymmaxp,\
                                       "Particles in species %d have y above the grid when depositing the source, max y = %e"%(js,y.max())
                            if self.nzlocal > 0:
                                z = pgroup.zp[i1:i2]
                                assert z.min() >= self.zmminp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z below the grid when depositing the source, min z = %e"%(js,z.min())
                                assert z.max() < self.zmmaxp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z above the grid when depositing the source, max z = %e"%(js,z.max())
                                       
             # ___________________________________
             # Depose currents in PXR

             t0 = MPI.Wtime()

             pxr.jx = self.fields.Jx
             pxr.jy = self.fields.Jy
             pxr.jz = self.fields.Jz

             if pxr.c_dim == 2:

               pxr.pxrdepose_currents_on_grid_jxjyjz_2d()

               #pxr.pxrdepose_currents_on_grid_jxjyjz_sub_openmp(f.Jx,f.Jy,f.Jz,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,
               #pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt)

             elif pxr.c_dim ==3:

               pxr.pxrdepose_currents_on_grid_jxjyjz()

             # Time statistics
             t1 = MPI.Wtime()
             self.time_stat_loc_array[2] += (t1-t0)

             # ___________________________________
              # Depose charge density in PXR if required
              
             if self.l_getrho : # Depose Rho in PXR
             
               t0 = MPI.Wtime()
             
               if pxr.c_dim == 2:

                 pxr.pxrdepose_rho_on_grid_sub_openmp_2d(f.Rho,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt,0)

               elif pxr.c_dim ==3:

                 pxr.rho = self.fields.Rho
                 pxr.pxrdepose_rho_on_grid()

               # Time statistics
               t1 = MPI.Wtime()
               self.time_stat_loc_array[12] += (t1-t0)

             #pxr.pxrdepose_rho_on_grid_sub_openmp_3d(f.Rho,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt,0)
             if self.current_cor: # Depose Rhoold_local in PXR
                 t0 = MPI.Wtime()
                 pxr.pxrdepose_rho_on_grid_sub_openmp_3d(f.Rhoold_local,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt,1)
                 t1 = MPI.Wtime()
                 self.time_stat_loc_array[12] += (t1-t0)
               
        else:

            for pgroup in pgroups:

                if w3d.js1fsapi >= 0: js1 = w3d.js1fsapi
                else:                 js1 = 0
                if w3d.js2fsapi >= 0: js2 = w3d.js2fsapi+1
                else:                 js2 = pgroup.ns

                jslist = kw.get('jslist',None)
                if jslist is None: jslist = range(js1,js2)

                for js in jslist:
                    n = pgroup.nps[js]
                    if n == 0: continue
                    if pgroup.ldts[js]:
                        indts = top.ndtstorho[pgroup.ndts[js]-1]
                        iselfb = pgroup.iselfb[js]
                        self.setsourcepforparticles(0,indts,iselfb)

                        if self.debug:
                            i1 = pgroup.ins[js]-1
                            i2 = pgroup.ins[js]+pgroup.nps[js]-1
                            if self.nxlocal > 0:
                                x = pgroup.xp[i1:i2]
                                if self.l4symtry: x = abs(x)
                                if self.solvergeom == w3d.RZgeom:
                                    y = pgroup.yp[i1:i2]
                                    x = sqrt(x**2 + y**2)
                                assert x.min() >= self.xmminp,\
                                       "Particles in species %d have x below the grid when depositing the source, min x = %e"%(js,x.min())
                                assert x.max() < self.xmmaxp,\
                                       "Particles in species %d have x above the grid when depositing the source, max x = %e"%(js,x.max())
                            if self.nylocal > 0:
                                y = pgroup.yp[i1:i2]
                                if self.l4symtry or self.l2symtry: y = abs(y)
                                assert y.min() >= self.ymminp,\
                                       "Particles in species %d have y below the grid when depositing the source, min y = %e"%(js,y.min())
                                assert y.max() < self.ymmaxp,\
                                       "Particles in species %d have y above the grid when depositing the source, max y = %e"%(js,y.max())
                            if self.nzlocal > 0:
                                z = pgroup.zp[i1:i2]
                                assert z.min() >= self.zmminp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z below the grid when depositing the source, min z = %e"%(js,z.min())
                                assert z.max() < self.zmmaxp+self.getzgridndts()[indts],\
                                       "Particles in species %d have z above the grid when depositing the source, max z = %e"%(js,z.max())

                        self.setsourcep(js,pgroup,self.getzgridndts()[indts])

        # --- Only finalize the source if lzero is true, which means the this
        # --- call to loadsource should be a complete operation.
        self.sourcepfinalized = False
        if lzero and lfinalize_rho: self.finalizesourcep()


    def apply_bndconditions(self,js,pg=None):
        if self.l_verbose:print me,'enter apply_ions_bndconditions'
        # --- apply boundary conditions
        if pg is None:
            pg = top.pgroup
        if pg.nps[js]==0:return
        self.apply_bnd_conditions(js,pg)
        if self.l_verbose:print me,'exit apply_ions_bndconditions'

    def apply_bnd_conditions(self,js,pg=None):
        if self.l_verbose:print me,'enter apply_bnd_conditions'
        if pg is None:
            pg = top.pgroup
        if pg.nps[js]==0:return
        il = pg.ins[js]-1
        iu = il+pg.nps[js]
        #stckxy3d(pg.nps[js],pg.xp[il:iu],w3d.xmmax,w3d.xmmin,w3d.dx,
        #              pg.yp[il:iu],w3d.ymmax,w3d.ymmin,w3d.dy,
        #              pg.zp[il:iu],w3d.zmminlocal,w3d.dz,
        #              pg.uxp[il:iu],pg.uyp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
        #              top.zgrid,top.zbeam,w3d.l2symtry,w3d.l4symtry,top.pboundxy,true)
        stckxy3d(pg,js,top.zbeam,true)
        partbndwithdata(pg.nps[js],pg.xp[il:iu],pg.uxp[il:iu],pg.gaminv[il:iu],
                        w3d.xmmaxlocal,w3d.xmminlocal,w3d.dx,0.,
                        top.pboundxy,top.pboundxy)
        partbndwithdata(pg.nps[js],pg.yp[il:iu],pg.uyp[il:iu],pg.gaminv[il:iu],
                        w3d.ymmaxlocal,w3d.ymminlocal,w3d.dy,0.,
                        top.pboundxy,top.pboundxy)
        if js==0 or js==w3d.nzp-1:
          if js==0:top.pboundnz=-1
          if js==w3d.nzp-1:top.pbound0=-1
          partbndwithdata(pg.nps[js],pg.zp[il:iu],pg.uzp[il:iu],pg.gaminv[il:iu],
                          w3d.zmmaxlocal,w3d.zmminlocal,w3d.dz,top.zgrid,
                          top.pbound0,top.pboundnz)
          if js==0:top.pboundnz=0
          if js==w3d.nzp-1:top.pbound0=0
        if self.scraper is not None:self.scraper.scrape(js)
        processlostpart(pg,js+1,top.clearlostpart,top.time+top.dt*pg.ndts[js],top.zbeam)
        if self.l_verbose:print me,'enter apply_bnd_conditions'

    def get_kinetic_energy(self,sp,**kw):
        """
        Get the total kinetic energy of the species sp using PICSAR fortran subroutines

        input:
        - sp: species number
        """
        total_kinetic_energy = zeros(1)
        if self.l_verbose:print me,'compute kinetic energy on species',sp
        pxr.get_kinetic_energy(sp,total_kinetic_energy)
        #print total_kinetic_energy,sp
        return total_kinetic_energy[0]

    def get_field_energy(self,field,**kw):
        """
        Get the total field energy for the given component.
        The field energy is calculated in parallel with a picsar fortran subroutine.

        input:
        - field: field component
        """
        field_energy = zeros(1)

        if pxr.c_dim==2:

          if field=='ex':
            pxr.get_field_energy_2d(self.fields.Ex,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='ey':
            pxr.get_field_energy_2d(self.fields.Ey,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='ez':
            pxr.get_field_energy_2d(self.fields.Ez,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='bx':
            pxr.get_field_energy_2d(self.fields.Bx,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='by':
            pxr.get_field_energy_2d(self.fields.By,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          elif field=='bz':
            pxr.get_field_energy_2d(self.fields.Bz,pxr.nx,pxr.nz,pxr.dx,pxr.dz,pxr.nxguards,pxr.nzguards,field_energy)
          return field_energy[0]

        else:

          if field=='ex':
            pxr.get_field_energy(self.fields.Ex,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='ey':
            pxr.get_field_energy(self.fields.Ey,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='ez':
            pxr.get_field_energy(self.fields.Ez,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='bx':
            pxr.get_field_energy(self.fields.Bx,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='by':
            pxr.get_field_energy(self.fields.By,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          elif field=='bz':
            pxr.get_field_energy(self.fields.Bz,pxr.nx,pxr.ny,pxr.nz,pxr.dx,pxr.dy,pxr.dz,pxr.nxguards,pxr.nyguards,pxr.nzguards,field_energy)
          return field_energy[0]

    def get_normL2_divEeps0_rho(self):
        """
        Compute the L2 norm of divE*eps0 - rho
        Computation of rho has to be activated

        """
        div = zeros(1)

        pxr.calc_field_div(pxr.dive,pxr.ex, pxr.ey, pxr.ez, pxr.nx,pxr.ny,pxr.nz,pxr.nxguards,pxr.nyguards,pxr.nzguards,pxr.dx,pxr.dy,pxr.dz)

        pxr.get_norm_diverho(pxr.dive,pxr.rho,pxr.nx,pxr.ny,pxr.nz,pxr.nxguards,pxr.nyguards,pxr.nzguards,div)

        return div[0]

    def display_picsar_time_statistics(self):
        """
        Display the Picsar time statistics
        """
        pxr.time_statistics()

    def display_time_statistics(self):
        """
        Display the time statistics
        """
        self.time_stat_ave_array = zeros([20])
        self.time_stat_min_array = zeros([20])
        self.time_stat_max_array = zeros([20])
        nproc = pxr.nprocx*pxr.nprocy*pxr.nprocz

        MPI.COMM_WORLD.Reduce([self.time_stat_loc_array,MPI.DOUBLE], [self.time_stat_ave_array,MPI.DOUBLE], op=MPI.SUM, root=0)
        MPI.COMM_WORLD.Reduce([self.time_stat_loc_array,MPI.DOUBLE], [self.time_stat_min_array,MPI.DOUBLE], op=MPI.MIN, root=0)
        MPI.COMM_WORLD.Reduce([self.time_stat_loc_array,MPI.DOUBLE], [self.time_stat_max_array,MPI.DOUBLE], op=MPI.MAX, root=0)

        self.time_stat_ave_array[:] /= nproc

        if me==0:

          print ' _____________________________________________________________'
          print
          print '  Time statisctics'
          print ' _____________________________________________________________'

          print ' Parts                              {:^8} {:^8} {:^8}'.format('min', 'ave', 'max')
          print ' -------------------------------------------------------------'
          print ' Particle pusher + field gathering: {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[0],self.time_stat_ave_array[0],self.time_stat_max_array[0])
          print ' Particle boundary conditions:      {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[1],self.time_stat_ave_array[1],self.time_stat_max_array[1])
          print ' Current deposition:                {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[2],self.time_stat_ave_array[2],self.time_stat_max_array[2])
          print ' Current bound. cond.:              {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[3],self.time_stat_ave_array[3],self.time_stat_max_array[3])
          print ' Magnetic field solver:             {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[5],self.time_stat_ave_array[5],self.time_stat_max_array[5])
          print ' Magnetic field bound. cond.:       {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[6],self.time_stat_ave_array[6],self.time_stat_max_array[6])
          print ' Electric field solver:             {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[7],self.time_stat_ave_array[7],self.time_stat_max_array[7])
          print ' Electric field bound. cond.:       {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[8],self.time_stat_ave_array[8],self.time_stat_max_array[8])
          print ' Particle sorting:                  {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[10],self.time_stat_ave_array[10],self.time_stat_max_array[10])
          print ' Charge deposition:                 {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[12],self.time_stat_ave_array[12],self.time_stat_max_array[12])
          print ' Charge bound. cond.:               {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[13],self.time_stat_ave_array[13],self.time_stat_max_array[13])
          print ' Load balancing:                    {:8.3f} {:8.3f} {:8.3f}'.format(self.time_stat_min_array[15],self.time_stat_ave_array[15],self.time_stat_max_array[15])
          print


    def allocatefieldarraysFFT(self):
        def fc(x,norder):
            fact1 = 1
            fact2 = 1
            result = 0
            for i in range(abs(norder)/2):
              fact1 *= max(i,1)
              fact2 *= max(2*i,1)*max(2*i-1,1)
              result += x**(2*i+1)*fact2/float(2**(2*i)*fact1**2*(2*i+1))
            return result


        f=self.fields
        b=self.block
        s=self
        f.spectral = (self.spectral > 0)
        bc_periodic = [self.bounds[0]==periodic,
                       self.bounds[2]==periodic,
                       self.bounds[4]==periodic]
        if self.current_cor:
            f.nxdrho = f.nx
            f.nydrho = f.ny
            f.nzdrho = f.nz
            f.nxdrhoguard = f.nxguard
            f.nydrhoguard = f.nyguard
            f.nzdrhoguard = f.nzguard
            f.gchange()

        if self.spectral:

            kwGPSTD = {'l_staggered':s.l_spectral_staggered,\
                     'l_staggered_a_la_brendan':s.l_staggered_a_la_brendan, \
                     'spectral':s.spectral,\
                     'norderx':s.norderx,\
                     'nordery':s.nordery,\
                     'norderz':s.norderz,\
                     'nxguard':s.nxguard,\
                     'nyguard':s.nyguard,\
                     'nzguard':s.nzguard,\
                     'dt':top.dt,\
                     'dx':w3d.dx,\
                     'dy':w3d.dy,\
                     'dz':w3d.dz,\
                     'ntsub':s.ntsub,\
                     'l_pushf':s.l_pushf,\
                     'l_pushg':s.l_pushg,\
                     'l_getrho':s.l_getrho,\
                     'clight':clight}

            if s.ntsub is np.inf:
                if not self.l_getrho:
                    self.l_getrho = True
                    f.nxr = f.nx
                    f.nyr = f.ny
                    f.nzr = f.nz
                    f.gchange()

                self.GPSTDMaxwell = gpstd.PSATD_Maxwell(yf=self.fields,
                                                  eps0=eps0,
                                                  bc_periodic=bc_periodic,
                                                  **kwGPSTD)
            else:
                if self.l_pushf and not self.l_getrho:
                    self.l_getrho = True
                    f.nxr = f.nx
                    f.nyr = f.ny
                    f.nzr = f.nz
                    f.gchange()
                self.GPSTDMaxwell = gpstd.GPSTD_Maxwell(yf=self.fields,
                                                  eps0=eps0,
                                                  bc_periodic=bc_periodic,
                                                  **kwGPSTD)

            self.FSpace = self.GPSTDMaxwell
        else:
            kwFS = {'l_staggered':s.l_spectral_staggered,\
                     'l_staggered_a_la_brendan':s.l_staggered_a_la_brendan, \
                     'spectral':s.spectral,\
                     'norderx':s.norderx,\
                     'nordery':s.nordery,\
                     'norderz':s.norderz,\
                     'nxguard':s.nxguard,\
                     'nyguard':s.nyguard,\
                     'nzguard':s.nzguard,\
                     'dt':top.dt,\
                     'dx':w3d.dx,\
                     'dy':w3d.dy,\
                     'nx':max([1,self.fields.nx]),\
                     'ny':max([1,self.fields.ny]),\
                     'nz':max([1,self.fields.nz]),\
                     'dz':w3d.dz}
            self.FSpace = Fourier_Space(bc_periodic=bc_periodic,**kwFS)

        # --- computes Brendan's Jz,Jx multipliers
        if self.Jmult and self.GPSTDMaxwell.nz>1:
                k = self.GPSTDMaxwell.k
                if self.GPSTDMaxwell.nx>1:kxvzdto2 = 0.5*self.GPSTDMaxwell.kx*clight*top.dt
                if self.GPSTDMaxwell.ny>1:kyvzdto2 = 0.5*self.GPSTDMaxwell.ky*clight*top.dt
                kzvzdto2 = 0.5*self.GPSTDMaxwell.kz*clight*top.dt
                sinkzvzdto2 = sin(kzvzdto2)
                coskzvzdto2 = cos(kzvzdto2)
                kdto2 = 0.5*k*clight*top.dt
                sinkdto2 = sin(kdto2)
                coskdto2 = cos(kdto2)
                numer = clight*top.dt*k*self.kz*(self.sinkzvzdto2**2-self.sinkdto2**2)
                denom = 2*sinkdto2*sinkzvzdto2 \
                      * (self.GPSTDMaxwell.kz*sinkzvzdto2*coskdto2-k*coskzvzdto2*sinkdto2)
                denomno0 = where(denom==0.,0.0001,self.denom)

                raise Exception("What is the 3-D version of Brendan's correction?")

                ktest=where((pi/2-kxvzdto2**2/(2*pi))>0,(pi/2-kxvzdto2**2/(2*pi)),0)

                Jmultiplier = where(abs(self.kzvzdto2)<ktest,numer/denomno0,0)

                self.Jmultiplier[0,:]=self.Jmultiplier[1,:]
                self.Jmultiplier[:,0]=self.Jmultiplier[:,1]

        # --- set Ex,By multipliers (ebcor=0,1,2)
        if self.l_correct_num_Cherenkov and self.spectral:
              emK = self.FSpace
#              k = emK.k
              k = sqrt(emK.kx_unmod*emK.kx_unmod+emK.ky_unmod*emK.ky_unmod+emK.kz_unmod*emK.kz_unmod)
              if top.boost_gamma==1.:
                  raise Exception('Error: l_correct_num_Cherenkov=True with top.boost_gamma=1.')

              b0 = sqrt(1.-1./top.boost_gamma**2)
              self.b0=b0
              self.ebcor = 2

              if 0:

              # --- old coefs
                  # --- set Ex,By multipliers (ebcor=0,1,2)
                  if self.ebcor==2:
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      self.Exmultiplier = self.kzvzdto2*self.coskzvzdto2/self.sinkzvzdto2

                      self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      self.sinkdto2 = sin(self.kdto2)
                      self.coskdto2 = cos(self.kdto2)
                      self.Bymultiplier = self.kdto2*self.coskdto2/self.sinkdto2

                  if self.ebcor==1:
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      self.sinkdto2 = sin(self.kdto2)
                      self.coskdto2 = cos(self.kdto2)
                      self.Exmultiplier = self.kdto2*self.sinkdto2**2*self.sinkzvzdto2*self.coskzvzdto2/ \
                        (self.kzvzdto2*(self.kdto2*self.sinkdto2**2+ \
                        (self.sinkdto2*self.coskdto2-self.kdto2)*self.sinkzvzdto2**2))

              else:
              # --- new cooefs
                  if self.ebcor==2:
                      # --- set Ex multiplier
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      self.Exmultiplier = self.kzvzdto2*self.coskzvzdto2/self.sinkzvzdto2
                      # --- set By multiplier
                      if self.norderx is inf:
                          self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      else:
                          self.kdto2 = sqrt((fc(sin(emK.kx_unmod*0.5*self.dx),self.norderx)/(0.5*self.dx))**2+ \
                              (fc(sin(emK.kz_unmod*0.5*self.dz),self.norderz)/(0.5*self.dz))**2)
                          self.kdto2 = where(self.kdto2==0,0.0001,0.5*self.kdto2*clight*top.dt)
                      if 0:#self.solver==PSATD:
                          self.Bymultiplier = self.kdto2/tan(self.kdto2)
                      else:
                          self.thetadto2=self.ntsub*arcsin(self.kdto2/self.ntsub)
                          self.Bymultiplier = self.kdto2/(tan(self.thetadto2)*cos(self.thetadto2/self.ntsub))

                  if self.ebcor==1:
                      self.kzvzdto2 = where(emK.kz_unmod==0,0.0001,0.5*emK.kz_unmod*b0*clight*top.dt)
                      self.sinkzvzdto2 = sin(self.kzvzdto2)
                      self.coskzvzdto2 = cos(self.kzvzdto2)
                      if self.norderx is None:
                          self.kdto2 = where(k==0,0.0001,0.5*k*clight*top.dt)
                      else:
                          self.kdto2 = sqrt((fc(sin(emK.kx_unmod*0.5*self.dx),self.norderx)/(0.5*self.dx))**2+ \
                              (fc(sin(emK.kz_unmod*0.5*self.dz),self.norderz)/(0.5*self.dz))**2)
                          self.kdto2 = where(self.kdto2==0,0.0001,0.5*self.kdto2*clight*top.dt)
                          self.kzvzdto2 = fc(sin(emK.kz_unmod*0.5*self.dz),self.norderz)/(0.5*self.dz)
                          self.kzvzdto2 = where(self.kzvzdto2==0,0.0001,0.5*self.kzvzdto2*b0*clight*top.dt)
                      if 0:#:self.solver==PSATD:
                          self.sinkdto2 = sin(self.kdto2)
                          self.coskdto2 = cos(self.kdto2)
                          self.Exmultiplier = self.kdto2*self.sinkdto2**2*self.sinkzvzdto2*self.coskzvzdto2/ \
                           (self.kzvzdto2*(self.kdto2*self.sinkdto2**2+ \
                           (self.sinkdto2*self.coskdto2-self.kdto2)*self.sinkzvzdto2**2))
                      else:
                          self.thetadto2=self.ntsub*arcsin(self.kdto2/self.ntsub)
                          self.Exmultiplier = self.ntsub*self.sinkzvzdto2*self.coskzvzdto2*sin(self.thetadto2)**2/ \
                           (self.kzvzdto2*(self.ntsub*sin(self.thetadto2)**2-self.sinkzvzdto2**2* \
                           (self.ntsub-sin(2*self.thetadto2)/sin(2*self.thetadto2/self.ntsub))))


        if 0:#self.spectral:
                  emK = self.FSpace
                  b0 = sqrt(1.-1./top.boost_gamma**2)
                  self.cut = 0.6
                  k = sqrt(emK.kx_unmod*emK.kx_unmod+emK.kz_unmod*emK.kz_unmod)
                  self.k_source_filter = where(k*self.dz/pi>self.cut*min(1.,self.dz/(b0*clight*top.dt)),0.,1.)
                  if self.l_getrho:emK.add_Sfilter('rho',self.k_source_filter)
                  emK.add_Sfilter('jx',self.k_source_filter)
                  emK.add_Sfilter('jy',self.k_source_filter)
                  emK.add_Sfilter('jz',self.k_source_filter)

        if self.spectral:
            kwPML = kwGPSTD

            if s.ntsub==inf:
                GPSTD_PML = gpstd.PSATD_Maxwell_PML
            else:
                GPSTD_PML = gpstd.GPSTD_Maxwell_PML

            # --- sides
            if b.xlbnd==openbc: s.xlPML = GPSTD_PML(syf=b.sidexl.syf,**kwPML)
            if b.xrbnd==openbc: s.xrPML = GPSTD_PML(syf=b.sidexr.syf,**kwPML)
            if b.ylbnd==openbc: s.ylPML = GPSTD_PML(syf=b.sideyl.syf,**kwPML)
            if b.yrbnd==openbc: s.yrPML = GPSTD_PML(syf=b.sideyr.syf,**kwPML)
            if b.zlbnd==openbc: s.zlPML = GPSTD_PML(syf=b.sidezl.syf,**kwPML)
            if b.zrbnd==openbc: s.zrPML = GPSTD_PML(syf=b.sidezr.syf,**kwPML)

            # --- edges
            if(b.xlbnd==openbc and b.ylbnd==openbc): s.xlylPML = GPSTD_PML(syf=b.edgexlyl.syf,**kwPML)
            if(b.xrbnd==openbc and b.ylbnd==openbc): s.xrylPML = GPSTD_PML(syf=b.edgexryl.syf,**kwPML)
            if(b.xlbnd==openbc and b.yrbnd==openbc): s.xlyrPML = GPSTD_PML(syf=b.edgexlyr.syf,**kwPML)
            if(b.xrbnd==openbc and b.yrbnd==openbc): s.xryrPML = GPSTD_PML(syf=b.edgexryr.syf,**kwPML)
            if(b.xlbnd==openbc and b.zlbnd==openbc): s.xlzlPML = GPSTD_PML(syf=b.edgexlzl.syf,**kwPML)
            if(b.xrbnd==openbc and b.zlbnd==openbc): s.xrzlPML = GPSTD_PML(syf=b.edgexrzl.syf,**kwPML)
            if(b.xlbnd==openbc and b.zrbnd==openbc): s.xlzrPML = GPSTD_PML(syf=b.edgexlzr.syf,**kwPML)
            if(b.xrbnd==openbc and b.zrbnd==openbc): s.xrzrPML = GPSTD_PML(syf=b.edgexrzr.syf,**kwPML)
            if(b.ylbnd==openbc and b.zlbnd==openbc): s.ylzlPML = GPSTD_PML(syf=b.edgeylzl.syf,**kwPML)
            if(b.yrbnd==openbc and b.zlbnd==openbc): s.yrzlPML = GPSTD_PML(syf=b.edgeyrzl.syf,**kwPML)
            if(b.ylbnd==openbc and b.zrbnd==openbc): s.ylzrPML = GPSTD_PML(syf=b.edgeylzr.syf,**kwPML)
            if(b.yrbnd==openbc and b.zrbnd==openbc): s.yrzrPML = GPSTD_PML(syf=b.edgeyrzr.syf,**kwPML)

            # --- corners
            if(b.xlbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc): s.xlylzlPML = GPSTD_PML(syf=b.cornerxlylzl.syf,**kwPML)
            if(b.xrbnd==openbc and b.ylbnd==openbc and b.zlbnd==openbc): s.xrylzlPML = GPSTD_PML(syf=b.cornerxrylzl.syf,**kwPML)
            if(b.xlbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc): s.xlyrzlPML = GPSTD_PML(syf=b.cornerxlyrzl.syf,**kwPML)
            if(b.xrbnd==openbc and b.yrbnd==openbc and b.zlbnd==openbc): s.xryrzlPML = GPSTD_PML(syf=b.cornerxryrzl.syf,**kwPML)
            if(b.xlbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc): s.xlylzrPML = GPSTD_PML(syf=b.cornerxlylzr.syf,**kwPML)
            if(b.xrbnd==openbc and b.ylbnd==openbc and b.zrbnd==openbc): s.xrylzrPML = GPSTD_PML(syf=b.cornerxrylzr.syf,**kwPML)
            if(b.xlbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc): s.xlyrzrPML = GPSTD_PML(syf=b.cornerxlyrzr.syf,**kwPML)
            if(b.xrbnd==openbc and b.yrbnd==openbc and b.zrbnd==openbc): s.xryrzrPML = GPSTD_PML(syf=b.cornerxryrzr.syf,**kwPML)


class Sorting:
  """
    Class Sorting

    Used to setup the sorting with picsars

    - activated: >0 sorting is activated
    - periods: list containing the sorting periods for each species
    - starts: first iteration before the start of the sorting
    - dx, dy, dz: the bin size normalized to the cell size. For instance, a dx of 1 corresponds to the cell dx.
    - xshift,yshift,zshift: shift of the sorting grid. The shift is normalized to dx,dy,dz. For instance a shift of 1 corresponds of 1 space step.

  """
  def __init__(self,periods,starts,activated=1,dx=1.,dy=1.,dz=1.,xshift=0.,yshift=0,zshift=0,verbose=False):
    self.activated = activated
    self.periods = periods
    self.starts = starts
    self.dx = dx
    self.dy = dy
    self.dz = dz
    self.xshift = xshift
    self.yshift = yshift
    self.zshift = zshift
    self.verbose = verbose
"""Class for 2D & 3D FFT-based electromagnetic solver """
from warp.field_solvers.em3dsolverFFT import *
try:
    import picsarpy as pxrpy
    pxr = pxrpy.picsar
    l_pxr=True
    print 'PICSAR package found and loaded.'
except:
    l_pxr=False
    print 'PICSAR package not found.'
try: 
    from mpi4py import MPI 
except: 
    print 'Error cannot import mpi4py'  
    
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
                      'l_output_freq':1
                      }

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(EM3DPXR.__flaginputs__,kw)
        EM3DFFT.__init__(self,kwdict=kw)
        
        self.l_pxr = l_pxr

    def finalize(self,lforce=False):
        if self.finalized and not lforce: return
        EM3DFFT.finalize(self)
        if l_pxr:self.allocatefieldarraysPXR()
        loadrho()
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
    
    	# Set up case dimensionality 
    	if (self.l_2dxz): 
    		pxr.c_dim=2
    	else: 
    		pxr.c_dim=3
        # Set up PXR MPI Data
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
        if (top.pbound0 == absorb): 
        	pxr.pbound_z_min=1
        else: # Default is periodic 
        	pxr.pbound_z_min=0
        	
        if (top.pboundnz == absorb): 
        	pxr.pbound_z_max=1
        else: # Default is periodic 
        	pxr.pbound_z_max=0
        	
        if (top.pboundxy == absorb): 
        	pxr.pbound_x_min=1
        	pxr.pbound_x_max=1
        	pxr.pbound_y_min=1
        	pxr.pbound_y_max=1
        else: # Default is periodic 
        	pxr.pbound_x_min=0
        	pxr.pbound_x_max=0
        	pxr.pbound_y_min=0
        	pxr.pbound_y_max=0
        
        # --- number of grid cells
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
        
        pxr.length_x=pxr.xmax-pxr.xmin
        pxr.length_y=pxr.ymax-pxr.ymin
        pxr.length_z=pxr.zmax-pxr.zmin
        
		# INIT MPI_DATA FOR PICSAR
        # Init communicator variable in picsar 
        pxr.mpi_minimal_init(top.fsdecomp.mpi_comm)
        
        # allocate grid quantities 
        pxr.allocate_grid_quantities()
        pxr.compute_simulation_axis()
        
        # set time step 
        pxr.dt = top.dt

        # --- Resolution
        pxr.dx = self.dx
        pxr.dy = self.dy
        pxr.dz = self.dz

        # --- Maxwell solver
        pxr.norderx = self.norderx
        pxr.nordery = self.nordery
        pxr.norderz = self.norderz
        
        pxr.xcoeffs = self.fields.xcoefs
        pxr.ycoeffs = self.fields.ycoefs
        pxr.zcoeffs = self.fields.zcoefs
        
        # Set coefficient for Maxwell solver
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

        # --- Tiling parameters
        pxr.ntilex = self.ntilex
        pxr.ntiley = self.ntiley
        pxr.ntilez = self.ntilez

        # --- species section
        pxr.nspecies_max=top.pgroup.ns
        
            # --- allocates array of species
        pxr.init_species_section() 
        
        for i,s in enumerate(self.listofallspecies):
            pxr.set_particle_species_properties(i+1,s.name,s.mass,s.charge,0, \
                                                0.,0.,0.,0.,0.,0., \
                                                0.,0.,0.,0.,0.,0.)
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

    def push_b_part_1(self,dir=1.):
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
												  0.5*dt/f.dx*f.xcoefs[0],
												  0.5*dt/f.dy*f.ycoefs[0],
												  0.5*dt/f.dz*f.zcoefs[0],
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
												  0.5*dt/f.dx*f.xcoefs[0],
												  0.5*dt/f.dy*f.ycoefs[0],
												  0.5*dt/f.dz*f.zcoefs[0],
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

    def push_b_part_2(self):
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
												  0.5*dt/f.dx*f.xcoefs[0],
												  0.5*dt/f.dy*f.ycoefs[0],
												  0.5*dt/f.dz*f.zcoefs[0],
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
												  0.5*dt/f.dx*f.xcoefs[0],
												  0.5*dt/f.dy*f.ycoefs[0],
												  0.5*dt/f.dz*f.zcoefs[0],
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


    def exchange_e(self,dir=1.):
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

    def exchange_b(self,dir=1.):
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
         

    def step(self,n=1,freq_print=10,lallspecl=0):
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
            	       


    def onestep(self,l_first,l_last):

        # --- call beforestep functions
        callbeforestepfuncs.callfuncsinlist()
    
        top.zgrid+=top.vbeamfrm*top.dt
        top.zbeam=top.zgrid
        # --- gather fields from grid to particles
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
                tdebpart=MPI.Wtime()
                pxr.pxrpush_particles_part2()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                pxr.particle_bcs()
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
                tdebpart=MPI.Wtime()
                pxr.push_particles()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
                pxr.particle_bcs()
                #for i,s in enumerate(self.listofallspecies):
                #    for pg in s.flatten(s.pgroups):
                #        particleboundaries3d(pg,-1,False)
                #pxr.particle_bcs_tiles()
                self.aliasparticlearrays()

            else:
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.push_velocity_full(0,pg)
                        self.push_positions(0,pg)
                        particleboundaries3d(pg,-1,False)

        # --- call beforeloadrho functions
        beforeloadrho.callfuncsinlist()
        pgroups = []
        for i,s in enumerate(self.listofallspecies):
            pgroups+=s.flatten(s.pgroups)
        self.pgroups = pgroups
#        self.loadsource(pgroups=pgroups)
        #tdebpart=MPI.Wtime()
        self.loadrho(pgroups=pgroups)
        self.loadj(pgroups=pgroups)
        #tendpart=MPI.Wtime()
        #pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
#        self.solve2ndhalf()
        #tdebcell=MPI.Wtime()
        self.dosolve()
        #tendcell=MPI.Wtime()
        #pxr.local_time_cell=pxr.local_time_cell+(tendcell-tdebcell)
    
        if l_pxr:
            if l_last:
                tdebpart=MPI.Wtime()
                pxr.pxrpush_particles_part1()
                tendpart=MPI.Wtime()
                pxr.local_time_part=pxr.local_time_part+(tendpart-tdebpart)
        else:
            for i,s in enumerate(self.listofallspecies):
                for pg in s.flatten(s.pgroups):
                    w3d.pgroupfsapi = pg
                    self.fetcheb(0,pg)
                    if l_last:
                        self.push_velocity_first_half(0,pg)

        # --- update time, time counter
        top.time+=top.dt
        if top.it%top.nhist==0:
#           zmmnt()
           minidiag(top.it,top.time,top.lspecial)
        top.it+=1
        
        # Load balance every dlb_freq time step
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
        if(l_pxr & self.l_output_grid & (top.it % self.l_output_freq ==0)): 
        	self.output_pxr(top.it)

        # --- call afterstep functions
        callafterstepfuncs.callfuncsinlist()

    def load_balance_3d(self,imbalance): 
        if (l_pxr): 
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
             	self.xmminlocal = pxr.x_min_local
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
        '''Charge deposition, uses particles from top directly
          - jslist: option list of species to load'''
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
             tdeb = MPI.Wtime()
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
             pxr.pxrdepose_currents_on_grid_jxjyjz_sub_openmp(f.Jx,f.Jy,f.Jz,pxr.nx,pxr.ny,pxr.nz,pxr.nxjguards,
             pxr.nyjguards,pxr.nzjguards,pxr.nox,pxr.noy,pxr.noz,pxr.dx,pxr.dy,pxr.dz,pxr.dt)  
             tend=MPI.Wtime()
             pxr.local_time_part+=(tend-tdeb)
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
            
                        

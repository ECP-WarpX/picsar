"""Class for 2D & 3D FFT-based electromagnetic solver """
from em3dsolverFFT import *
try:
    import picsarpy as pxrpy
    pxr = pxrpy.picsar
    l_pxr=True
    print 'PICSAR package found and loaded.'
except:
    l_pxr=False
    print 'PICSAR package not found.'

class EM3DPXR(EM3DFFT):

    __em3dpxrinputs__ = []
    __flaginputs__ = {'ntilex':1,
                      'ntiley':1,
                      'ntilez':1,
                      'listofallspecies':[],
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

    def allocatefieldarraysPXR(self):
        # --- number of grid cells
        pxr.nx_global = w3d.nx
        pxr.ny_global = w3d.ny
        pxr.nz_global = w3d.nz
        pxr.nx_global_grid = pxr.nx_global+1
        pxr.ny_global_grid = pxr.ny_global+1
        pxr.nz_global_grid = pxr.nz_global+1
        
        pxr.nx = self.nx
        pxr.ny = self.ny
        pxr.nz = self.nz

        pxr.nprocx=0
        pxr.nprocy=0
        pxr.nprocz=0
		
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
        
        pxr.x_min_local = self.fields.xmin
        pxr.x_max_local = self.fields.xmax
        pxr.y_min_local = self.fields.ymin
        pxr.y_max_local = self.fields.ymax
        pxr.z_min_local = self.fields.zmin
        pxr.z_max_local = self.fields.zmax
		
		
        pxr.mpi_minimal_init()
        pxr.mpi_initialise()
#        print 'nprocx',pxr.nprocx
        

        
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
        
        pxr.ex = self.fields.Ex
        pxr.ey = self.fields.Ey
        pxr.ez = self.fields.Ez
        pxr.bx = self.fields.Bx
        pxr.by = self.fields.By
        pxr.bz = self.fields.Bz
        
        pxr.l_nodalgrid = self.l_nodalgrid
        
        pxr.nxs = 0
        pxr.nys = 0
        pxr.nzs = 0


		# Current deposition
        pxr.nox=1
        pxr.noy=1
        pxr.noz=1

        # --- Tiling parameters
        pxr.ntilex = self.ntilex
        pxr.ntiley = self.ntiley
        pxr.ntilez = self.ntilez

        # --- species section
        pxr.nspecies_max=top.pgroup.ns
        
            # --- allocates array of species
        pxr.init_species_section() 
        
        for i,s in enumerate(self.listofallspecies):
            pxr.set_particle_species_properties(i+1,s.name,s.mass,s.charge,1, \
                                                0.,0.,0.,0.,0.,0., \
                                                0.,0.,0.,0.,0.,0.)
            pxr.nspecies+=1

        pxr.set_tile_split()
        pxr.init_tile_arrays()
        
        for i,s in enumerate(self.listofallspecies):
            print 'add particles %g to pxr'%i
            pxr.py_add_particles_to_species(i+1, s.nps, 
                                            s.getx(bcast=0,gather=0), 
                                            s.gety(bcast=0,gather=0), 
                                            s.getz(bcast=0,gather=0), 
                                            s.getux(bcast=0,gather=0), 
                                            s.getuy(bcast=0,gather=0), 
                                            s.getuz(bcast=0,gather=0), 
                                            s.getweights(bcast=0,gather=0)  )      

        top.pgroup.npmax=0
        top.pgroup.ns=1
        top.pgroup.nps=0
        top.pgroup.gchange()
        
        # --- mirror PXR tile structure in Warp with list of pgroups
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
                        pg.sw = 1.
                        pg.npmax = pxr.partnmax
                        pg.nps = pxr.partn
                        pg.ins[0] = 1
                        pg.xp = pxr.partx
                        pg.yp = pxr.party
                        pg.zp = pxr.partz
                        pg.uxp = pxr.partux
                        pg.uyp = pxr.partuy
                        pg.uzp = pxr.partuz
                        pg.pid = fzeros([pg.npmax,top.npid])
                        pg.pid[:,top.wpid-1] = pxr.partw
                        pg.gaminv = zeros(pg.npmax)
                        self.set_gamma(js=0,pg=pg)
                        pg.ex = pxr.partex
                        pg.ey = pxr.partey
                        pg.ez = pxr.partez
                        pg.bx = pxr.partbx
                        pg.by = pxr.partby
                        pg.bz = pxr.partbz
                        pg.lebcancel_pusher=top.pgroup.lebcancel_pusher
                    xygroup.append(xgroup)
                s.pgroups.append(xygroup)

#        for i,s in enumerate(self.listofallspecies):
#            def ppzx(self,**kw):
#                for pg in self.pgroups:
#                   self._callppfunc(ppzx,pgroup=pg,**kw)
                    
#            s.ppzx = ppzx
            
    def print_nptiles(self,ispecies):
        for iz in range(1,self.ntilez+1):
            for iy in range(1,self.ntiley+1):
                for ix in range(1,self.ntilex+1):
                    pxr.point_to_tile(ispecies, ix, iy, iz)
                    print ix,iy,iz,pxr.partn[0]
    
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
                if self.l_2dxz:
                    pxr.pxrpush_em2d_evec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              f.J[...,0],f.J[...,1],f.J[...,2],
                                              clight**2*mu0*dt,        
                                              clight**2*dt/f.dx*f.xcoefs,
                                              clight**2*dt/f.dy*f.ycoefs,
                                              clight**2*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
                else:
                    pxr.pxrpush_em3d_evec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              f.J[...,0],f.J[...,1],f.J[...,2],
                                              clight**2*mu0*dt,        
                                              clight**2*dt/f.dx*f.xcoefs,
                                              clight**2*dt/f.dy*f.ycoefs,
                                              clight**2*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
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
                f=self.fields
                l_pushb=False
                if self.l_2dxz:
                    pxr.pxrpush_em2d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              0.5*dt/f.dx*f.xcoefs,
                                              0.5*dt/f.dy*f.ycoefs,
                                              0.5*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
                else:
                    pxr.pxrpush_em3d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              0.5*dt/f.dx*f.xcoefs,
                                              0.5*dt/f.dy*f.ycoefs,
                                              0.5*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
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
                if self.l_2dxz:
                    pxr.pxrpush_em2d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              0.5*dt/f.dx*f.xcoefs,
                                              0.5*dt/f.dy*f.ycoefs,
                                              0.5*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
                else:
                    pxr.pxrpush_em3d_bvec_norder(f.Ex,f.Ey,f.Ez,f.Bx,f.By,f.Bz,
                                              0.5*dt/f.dx*f.xcoefs,
                                              0.5*dt/f.dy*f.ycoefs,
                                              0.5*dt/f.dz*f.zcoefs,
                                              f.nx,f.ny,f.nz,
                                              f.norderx,f.nordery,f.norderz,
                                              f.nxguard,f.nyguard,f.nzguard,
                                              0,0,0,f.l_nodalgrid)
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
        for i in range(n):
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
       
    def onestep(self,l_first,l_last):
        # --- call beforestep functions
        callbeforestepfuncs.callfuncsinlist()
    
        top.zgrid+=top.vbeamfrm*top.dt
        top.zbeam=top.zgrid

        # --- gather fields from grid to particles
#        w3d.pgroupfsapi = top.pgroup
#        for js in range(top.pgroup.ns):
#          self.fetcheb(js)

        # --- push 
        if l_first:
            if l_pxr:
                pxr.pxrpush_particles_part2()
                #pxr.particle_bcs()
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.set_gamma(0,pg)
                        particleboundaries3d(pg,-1,False)
            else:
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.push_velocity_second_half(0,pg)
                        self.push_positions(0,pg)
                        particleboundaries3d(pg,-1,False)
        else:        
            if l_pxr:
                pxr.push_particles()
                #pxr.particle_bcs()
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.set_gamma(0,pg)
                        particleboundaries3d(pg,-1,False)
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
        self.loadrho(pgroups=pgroups)
        self.loadj(pgroups=pgroups)
                      
#        self.solve2ndhalf()
        self.dosolve()
    
        if l_pxr:
            if l_last:
                print top.it,'part 1'
                pxr.pxrpush_particles_part1()
                for i,s in enumerate(self.listofallspecies):
                    for pg in s.flatten(s.pgroups):
                        self.set_gamma(0,pg)
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

        # --- call afterstep functions
        callafterstepfuncs.callfuncsinlist()

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

"""Class for 2D & 3D FFT-based electromagnetic solver """
from warp_init_tools.plasma_initialization import * 
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
import numpy as np 

class PlasmaInjectorPXR(PlasmaInjector):
    def continuous_injection(self):
        """
        Routine which is called by warp at each timestep
        """
        # Move the injection position with the moving window
        self.z_inject += self.top.vbeamfrm * self.top.dt
        # Move the position of the end of the plasma by its mean velocity
        self.z_end_plasma += self.v_plasma * self.top.dt
        # Move the position of the limit beyond which no plasma is injected
        # (It moves along with the plasma because the user gives p_zmax at t=0)
        self.p_zmax += self.v_plasma * self.top.dt

        # Add slices filled with plasma
        while (self.z_end_plasma < self.z_inject) and \
            (self.z_end_plasma < self.p_zmax) :

            # Add one slice
            self.load_plasmaPXR( self.z_end_plasma + self.w3d.dz,
                              self.z_end_plasma,
                              self.z_end_plasma + self.w3d.dz )

            # One slice has been added ; increment the position of plasma end
            self.z_end_plasma += self.w3d.dz

    def load_plasmaPXR( self, z_end_plasma, zmin, zmax ):
        """
        Load plasma between zmin and zmax in PXR.
        
        The positions of the particles along the z axis are of the form
        z = z_end_plasma - i*dz - 0.5*dz
        and satisfy zmin <= z < zmax
        where i is an integer, and dz is the spacing
        between particles

        Parameters
        ----------
        z_end_plasma : float
           Position of the global end of the plasma
        
        zmin, zmax : floats (meters)
           Positions between which the plasma is to be loaded, in the
           local domain
        """
        # Get 1d array of evenly-spaced positions for the particles along z
        dz = self.w3d.dz / self.p_nz
        # Get the min and max indices i for z = z_end_plasma - i*dz - 0.5*dz
        i_min = int( (z_end_plasma-zmin)/dz - 0.5 ) 
        i_max = int( (z_end_plasma-zmax)/dz + 0.5 )
        i_max = max( i_max, 0 )
        z_reg = z_end_plasma - dz*( np.arange( i_max, i_min+1 ) + 0.5 )

        # Get the corresponding particle positions at injection
        if self.dim == "3d":
            zp, xp, yp = np.meshgrid( z_reg, self.x_reg, self.y_reg )
            z0 = zp.flatten()
            x0 = xp.flatten()
            y0 = yp.flatten()
            r0 = np.sqrt( x0**2 + y0**2 )
        elif self.dim == "circ":
            # (copy=True is important here, since it allows to
            # change the theta angles individually)
            zp, rp, thetap = np.meshgrid( z_reg, self.x_reg,
                                          self.theta_reg, copy=True )
            # Prevent the particles from being aligned along any direction
            theta = unalign_angles( thetap )
            r0 = rp.flatten()
            z0 = zp.flatten()
            x0 = r0 * np.cos( theta )
            y0 = r0 * np.sin( theta )
        elif self.dim == "2d":
            zp, xp = np.meshgrid( z_reg, self.x_reg )
            z0 = zp.flatten()
            x0 = xp.flatten()
            y0 = np.zeros_like(x0)
            r0 = abs(x0)
        elif self.dim == "1d":
            z0 = z_reg
            x0 = y0 = np.zeros_like(z0)
            r0 = abs(x0)

        # Modulate the weights according to the density
        # (Take into account the motion of the plasma to retrieve the
        # the position where this slice of plasma was at t=0, this
        # is because the dens_func is given at t=0)
        if self.dens_func is not None:
            w = self.dens_func( x0, y0, z0 - self.v_plasma*self.top.time )
        else:
            w = np.ones_like( z0 - self.v_plasma*self.top.time )
        # In circ, the particles have larger weight at higher radius
        if self.dim == "circ":
            w = w * r0 / self.w3d.xmmax

        # Add random momenta for the particles
        ux = self.ux_m
        uy = self.uy_m
        uz = self.uz_m
        gamma_inv = self.gamma_inv_m
        n_part = len(z0)
        if self.ux_th != 0:
            ux += self.ux_th * np.random.normal( size=n_part )
        if self.uy_th != 0:
            uy += self.uy_th * np.random.normal( size=n_part )
        if self.uz_th != 0:
            uz += self.uz_th * np.random.normal( size=n_part )
        if (self.ux_th !=0) or (self.uy_th !=0) or (self.uy_th !=0):
            gamma_inv = 1./np.sqrt( 1 + ux**2 + uy**2 + uz**2 )

        # Load the electrons and ions (on top of each other)
        if self.elec is not None:
            # Use the random momenta

            # Filter particles outside the box 
            cond=(x0>=pxr.x_min_local) & (x0<pxr.x_max_local) &\
                 (y0>=pxr.y_min_local) & (y0<pxr.y_max_local) &\
                 (z0>=pxr.z_min_local+pxr.zgrid) & (z0<pxr.z_max_local+pxr.zgrid) 
            x0=x0[cond]
            y0=y0[cond]
            z0=z0[cond]
            nps0=np.size(x0)
            #print("size",nps0, np.amax(np.abs(uz)),np.amax(np.abs(uy)),np.amax(np.abs(ux)), \
            #np.amax(z0),np.amax(y0),np.amax(x0),np.amin(z0),np.amin(y0),np.amin(x0))
            pxr.py_add_particles_to_species(1, nps0, 
                                            x0, 
                                            y0, 
                                            z0, 
                                            ux*np.ones(nps0), 
                                            uy*np.ones(nps0), 
                                            uz*np.ones(nps0),
                                            gamma_inv*np.ones(nps0), 
                                            w)     
        if self.ions is not None:
            # For each element, only add particles to the lowest charge state
            for element in self.ions.keys():
                # TO BE DONE SEARCH FOR INDEX SPECIES IN LISTOFALLSPECIES
                # Use only the mean momenta
                lowest_state_species = self.ions[ element ][0]
                lowest_state_species.addpart( x=x0, y=y0, z=z0,
                    vx=c*self.ux_m*self.gamma_inv_m,
                    vy=c*self.uy_m*self.gamma_inv_m,
                    vz=c*self.uz_m*self.gamma_inv_m,
                    gi=self.gamma_inv_m, w=w,
                    lallindomain=False )
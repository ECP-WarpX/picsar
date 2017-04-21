"""
 _______________________________________________________________________________

 *** Copyright Notice ***

 "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
 2016, The Regents of the University of California, through Lawrence Berkeley
 National Laboratory (subject to receipt of any required approvals from the
 U.S. Dept. of Energy). All rights reserved.

 If you have questions about your rights to use or distribute this software,
 please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

 NOTICE.
 This Software was developed under funding from the U.S. Department of Energy
 and the U.S. Government consequently retains certain rights. As such, the U.S.
 Government has been granted for itself and others acting on its behalf a
 paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 reproduce, distribute copies to the public, prepare derivative works, and
 perform publicly and display publicly, and to permit other to do so.


 Class for 2D & 3D FFT-based electromagnetic solver

 Developers:
 Henri Vincenti

 Date:
 Creation 2016

 _______________________________________________________________________________

"""
from warp.field_solvers.GPSTD import *
try:
    from picsar_python import picsarpy as pxrpy
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

class GPSTDPXR(GPSTD):

    def create_fortran_matrix_blocks(self):
        mymat = self.mymat
        nrow=len(mymat)
        # Allocate new block matrix in Fortran
        # And corresponding vector blocks
        pxr.allocate_new_matrix_vector(nrow)
        # Get Matrix index
        self.matrix_index = pxr.nmatrixes
        # Alias Fortran blocks to mymat elements
        for i in range(1,nrow+1):
            for j in range(1,nrow+1):
                mymat[i-1][j-1] = np.asarray(mymat[i-1][j-1])
                if (mymat[i-1][j-1].ndim == 2):
                    ki = self.fields_name[i-1]
                    n1=self.fields[ki].shape[0]
                    n2=self.fields[ki].shape[1]
                    n3=self.fields[ki].shape[2]
                    n1mymat=mymat[i-1][j-1].shape[0]
                    n2mymat=mymat[i-1][j-1].shape[1]
                    if (n1==1):
                        mymat[i-1][j-1]=np.reshape(mymat[i-1][j-1],(1,n1mymat,n2mymat))
                    elif(n2==1):
                        mymat[i-1][j-1]=np.reshape(mymat[i-1][j-1],(n1mymat,1,n2mymat))
                    elif(n3==1):
                        mymat[i-1][j-1]=np.reshape(mymat[i-1][j-1],(n1mymat,n2mymat,1))
                else:
                    if (np.size(mymat[i-1][j-1])==1):
                        mymat[i-1][j-1]=np.reshape(mymat[i-1][j-1],(1,1,1))
                if mymat[i-1][j-1].dtype is not np.dtype('complex128'):
                    mymat[i-1][j-1]=mymat[i-1][j-1].astype(np.complex128)
                mymat[i-1][j-1]=np.asfortranarray(mymat[i-1][j-1])
                n1 = mymat[i-1][j-1].shape[0]
                n2 = mymat[i-1][j-1].shape[1]
                n3 = mymat[i-1][j-1].shape[2]
                pxr.point_to_matrix_block_p2f(self.mymat[i-1][j-1],n1,n2,n3,i,j,self.matrix_index)

    def get_Ffields(self):
        ixl,ixu,iyl,iyu,izl,izu = self.get_ius()
        if self.Ffields=={}:
            self.fields_shape = [ixu-ixl,iyu-iyl,izu-izl]
            for k in self.fields.keys():
                self.plan_rfftn[k] = self.create_plan_rfftn(np.asarray(self.fields_shape))
                self.Ffields[k]    =self.rfftn(self.fields[k][ixl:ixu,iyl:iyu,izl:izu],plan=self.plan_rfftn[k])
        else:
            for k in self.fields.keys():
                self.Ffields[k]=self.rfftn(self.fields[k][ixl:ixu,iyl:iyu,izl:izu],field_out=self.Ffields[k],plan=self.plan_rfftn[k])

    def get_fields(self):
        ixl,ixu,iyl,iyu,izl,izu = self.get_ius()
        if (self.plan_irfftn=={}):
            for k in self.fields.keys():
                    self.plan_irfftn[k] = self.create_plan_irfftn(np.asarray(self.fields_shape))
        for k in self.fields.keys():
            if not self.LSource[k]:
                shapek = np.asarray(np.shape(self.fields[k][ixl:ixu,iyl:iyu,izl:izu]))
                self.fields[k][ixl:ixu,iyl:iyu,izl:izu] = self.irfftn(self.Ffields[k], shapek, field_out=self.fields[k][ixl:ixu,iyl:iyu,izl:izu], plan=self.plan_irfftn[k])

    def push_fields(self):

        # --- Fourier transforming fields
        self.get_Ffields()
        # --- filter sources before push
        for k in self.Sfilters.keys():
           self.Ffields[k]*=self.Sfilters[k]
        mymat = self.mymat
        n = len(mymat)
        # --- set dictionary of field values before time step
        oldfields = {}
        for k in self.Ffields.keys():
            oldfields[k] = self.Ffields[k].copy(order='F')
        # --- set dictionary of field flags for update
        updated_fields = {}
        for k in self.Ffields.keys():
            updated_fields[k] = False
        # --- Alias block vectors in Fortran
        for i in range(1,n+1):
            ki = self.fields_name[i-1]
            n1r=self.fields[ki].shape[0]
            n2r=self.fields[ki].shape[1]
            n3r=self.fields[ki].shape[2]
            nfs = self.Ffields[ki].ndim
            if (nfs < 3):
                n1=self.Ffields[ki].shape[0]
                n2=self.Ffields[ki].shape[1]
                if(n1r==1):
                    self.Ffields[ki]=np.reshape(self.Ffields[ki],(1,n1,n2))
                    oldfields[ki]=np.reshape(oldfields[ki],(1,n1,n2))
                elif(n2r==1):
                    self.Ffields[ki]=np.reshape(self.Ffields[ki],(n1,1,n2))
                    oldfields[ki]=np.reshape(oldfields[ki],(n1,1,n2))
                elif(n3r==1):
                    self.Ffields[ki]=np.reshape(self.Ffields[ki],(n1,n2,1))
                    oldfields[ki]=np.reshape(oldfields[ki],(n1,n2,1))
            nn1=self.Ffields[ki].shape[0]
            nn2=self.Ffields[ki].shape[1]
            nn3=self.Ffields[ki].shape[2]
            pxr.point_to_vector_block_p2f(self.Ffields[ki],nn1,nn2,nn3,i, \
                                        self.matrix_index,False,self.LSource[ki])
            pxr.point_to_vector_block_p2f(oldfields[ki],nn1,nn2,nn3,i, \
                                        self.matrix_index,True,self.LSource[ki])
        # --- fields update in FORTRAN
        pxr.multiply_mat_vector(self.matrix_index)

        # Deleting old copies of the fields
        del oldfields

        # --- filter fields after push
        for k in self.Ffilters.keys():
           self.Ffields[k]*=self.Ffilters[k]
        # Fourier transforming back fields
        self.get_fields()
        # --- set periodic BC
        if self.bc_periodic[0]:
            ngx = self.nxguard
        else:
            ngx = 0
        if self.bc_periodic[1]:
            ngy = self.nyguard
        else:
            ngy = 0
        if self.bc_periodic[2]:
            ngz = self.nzguard
        else:
            ngz = 0

        if self.bc_periodic[0]:
            for k in self.fields.keys():
                if updated_fields[k]:
                    f = self.fields[k]
                    f[-ngx-1:,...]=f[ngx:2*ngx+1,...]
                    f[:ngx,...]=f[-2*ngx-1:-ngx-1,...]
        if self.bc_periodic[1]:
            for k in self.fields.keys():
                if updated_fields[k]:
                    f = self.fields[k]
                    f[:,-ngy-1:,:]=f[:,ngy:2*ngy+1,:]
                    f[:,:ngy,:]=f[:,-2*ngy-1:-ngy-1,:]
        if self.bc_periodic[2]:
            for k in self.fields.keys():
                if updated_fields[k]:
                    f = self.fields[k]
                    f[...,-ngz-1:]=f[...,ngz:2*ngz+1]
                    f[...,:ngz]=f[...,-2*ngz-1:-ngz-1]

        del updated_fields

class PSATD_Maxwell_PML(GPSTDPXR):

    __flaginputs__ = {'syf':None,'l_pushf':False,'l_pushg':False,'clight':299792458.0}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(GPSTD_Maxwell_PML.__flaginputs__,kw)
        syf=self.syf
        nx = np.max([1,syf.nx])
        ny = np.max([1,syf.ny])
        nz = np.max([1,syf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        dt=self.dt
        cdt=dt*self.clight
        self.wdt = self.k*cdt
        self.coswdt=np.cos(self.wdt)
        self.sinwdt=np.sin(self.wdt)

        j = 1j

        if self.l_pushf:
            self.add_fields({"exx":syf.exx, \
                             "exy":syf.exy, \
                             "exz":syf.exz, \
                             "eyx":syf.eyx, \
                             "eyy":syf.eyy, \
                             "eyz":syf.eyz, \
                             "ezx":syf.ezx, \
                             "ezy":syf.ezy, \
                             "ezz":syf.ezz})
        else:
            self.add_fields({"exy":syf.exy, \
                             "exz":syf.exz, \
                             "eyx":syf.eyx, \
                             "eyz":syf.eyz, \
                             "ezx":syf.ezx, \
                             "ezy":syf.ezy})
        if self.l_pushg:
            self.add_fields({"bxx":syf.bxx, \
                             "bxy":syf.bxy, \
                             "bxz":syf.bxz, \
                             "byx":syf.byx, \
                             "byy":syf.byy, \
                             "byz":syf.byz, \
                             "bzx":syf.bzx, \
                             "bzy":syf.bzy, \
                             "bzz":syf.bzz})
        else:
            self.add_fields({"bxy":syf.bxy, \
                             "bxz":syf.bxz, \
                             "byx":syf.byx, \
                             "byz":syf.byz, \
                             "bzx":syf.bzx, \
                             "bzy":syf.bzy})
        if self.l_pushf:
            self.add_fields({"fx":syf.fx, \
                             "fy":syf.fy, \
                             "fz":syf.fz})

        if self.l_pushg:
            self.add_fields({"gx":syf.gx, \
                             "gy":syf.gy, \
                             "gz":syf.gz})

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt
        cdt=dt*self.clight
        C=self.coswdt
        S=self.sinwdt

        if self.nx>1:
            axm = j*S*self.kxmn
            axp = j*S*self.kxpn
        else:
            axm = axp = 0.

        if self.ny>1:
            aym = j*S*self.kymn
            ayp = j*S*self.kypn
        else:
            aym = ayp = 0.

        if self.nz>1:
            azm = j*S*self.kzmn
            azp = j*S*self.kzpn
        else:
            azm = azp = 0.

        self.mymat = self.getmaxwellmat_pml(C,S,axp,ayp,azp,axm,aym,azm)
        self.create_fortran_matrix_blocks()

    def getmaxwellmat_pml(self,C,S,axp,ayp,azp,axm,aym,azm):
        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            # --- bx
            if self.l_pushg:mymat.add_op('bxx',{'bxx':C,'gx':axm,'gy':axm,'gz':axm})
            mymat.add_op('bxy',{'bxy':C,'ezx':-ayp,'ezy':-ayp,'ezz':-ayp})
            mymat.add_op('bxz',{'bxz':C,'eyx': azp,'eyy': azp,'eyz': azp})
            # --- by
            mymat.add_op('byx',{'byx':C,'ezx': axp,'ezy': axp,'ezz': axp})
            if self.l_pushg:mymat.add_op('byy',{'byy':C,'gx':aym,'gy':aym,'gz':aym})
            mymat.add_op('byz',{'byz':C,'exx':-azp,'exy':-azp,'exz':-azp})
            # --- bz
            mymat.add_op('bzx',{'bzx':C,'eyx':-axp,'eyy':-axp,'eyz':-axp})
            mymat.add_op('bzy',{'bzy':C,'exx': ayp,'exy': ayp,'exz': ayp})
            if self.l_pushg:mymat.add_op('bzz',{'bzz':C,'gx':azm,'gy':azm,'gz':azm})
        else:
            # --- bx
            if self.l_pushg:mymat.add_op('bxx',{'bxx':C,'gx':axm,'gy':axm,'gz':axm})
            mymat.add_op('bxy',{'bxy':C,'ezx':-ayp,'ezy':-ayp})
            mymat.add_op('bxz',{'bxz':C,'eyx': azp,'eyz': azp})
            # --- by
            mymat.add_op('byx',{'byx':C,'ezx': axp,'ezy': axp})
            if self.l_pushg:mymat.add_op('byy',{'byy':C,'gx':aym,'gy':aym,'gz':aym})
            mymat.add_op('byz',{'byz':C,'exy':-azp,'exz':-azp})
            # --- bz
            mymat.add_op('bzx',{'bzx':C,'eyx':-axp,'eyz':-axp})
            mymat.add_op('bzy',{'bzy':C,'exy': ayp,'exz': ayp})
            if self.l_pushg:mymat.add_op('bzz',{'bzz':C,'gx':azm,'gy':azm,'gz':azm})

        if self.l_pushg:
            # --- ex
            if self.l_pushf:mymat.add_op('exx',{'exx':C,'fx':axp,'fy':axp,'fz':axp})
            mymat.add_op('exy',{'exy':C,'bzx': aym,'bzy': aym,'bzz': aym})
            mymat.add_op('exz',{'exz':C,'byx':-azm,'byy':-azm,'byz':-azm})
            # --- ey
            mymat.add_op('eyx',{'eyx':C,'bzx':-axm,'bzy':-axm,'bzz':-axm})
            if self.l_pushf:mymat.add_op('eyy',{'eyy':C,'fx':ayp,'fy':ayp,'fz':ayp})
            mymat.add_op('eyz',{'eyz':C,'bxx': azm,'bxy': azm,'bxz': azm})
            # --- ez
            mymat.add_op('ezx',{'ezx':C,'byx': axm,'byy': axm,'byz': axm})
            mymat.add_op('ezy',{'ezy':C,'bxx':-aym,'bxy':-aym,'bxz':-aym})
            if self.l_pushf:mymat.add_op('ezz',{'ezz':C,'fx':azp,'fy':azp,'fz':azp})
        else:
            # --- ex
            if self.l_pushf:mymat.add_op('exx',{'exx':C,'fx':axp,'fy':axp,'fz':axp})
            mymat.add_op('exy',{'exy':C,'bzx': aym,'bzy': aym})
            mymat.add_op('exz',{'exz':C,'byx':-azm,'byz':-azm})
            # --- ey
            mymat.add_op('eyx',{'eyx':C,'bzx':-axm,'bzy':-axm})
            if self.l_pushf:mymat.add_op('eyy',{'eyy':C,'fx':ayp,'fy':ayp,'fz':ayp})
            mymat.add_op('eyz',{'eyz':C,'bxy': azm,'bxz': azm})
            # --- ez
            mymat.add_op('ezx',{'ezx':C,'byx': axm,'byz': axm})
            mymat.add_op('ezy',{'ezy':C,'bxy':-aym,'bxz':-aym})
            if self.l_pushf:mymat.add_op('ezz',{'ezz':C,'fx':azp,'fy':azp,'fz':azp})

        if self.l_pushf:
            mymat.add_op('fx',{'fx':C,'exx':axm,'exy':axm,'exz':axm})
            mymat.add_op('fy',{'fy':C,'eyx':aym,'eyy':aym,'eyz':aym})
            mymat.add_op('fz',{'fz':C,'ezx':azm,'ezy':azm,'ezz':azm})

        if self.l_pushg:
            mymat.add_op('gx',{'gx':C,'bxx':axp,'bxy':axp,'bxz':axp})
            mymat.add_op('gy',{'gy':C,'byx':ayp,'byy':ayp,'byz':ayp})
            mymat.add_op('gz',{'gz':C,'bzx':azp,'bzy':azp,'bzz':azp})

        return mymat.mat

    def push(self):
        syf = self.syf

        self.push_fields()

        for f in self.fields.values():
            if self.nx>1:
                f[:self.nxguard,...]=0.
                f[-self.nxguard/2:,...]=0.
            if self.ny>1:
                f[:,:self.nyguard,:]=0.
                f[:,-self.nyguard/2:,:]=0.
            if self.nz>1:
                f[...,:self.nzguard/2]=0.
                f[...,-self.nzguard/2:]=0.

#      scale_em3d_split_fields(syf,top.dt,self.l_pushf)

        return

class GPSTD_Maxwell(GPSTDPXR):

    __flaginputs__ = {'yf':None,
                      'l_pushf':False,
                      'l_pushg':False,
                      'clight':299792458.0,
                      'eps0':8.854187817620389e-12,
                       'V_galilean':np.array([0.,0.,0.]),
                      'V_pseudogalilean':np.array([0.,0.,0.])}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(GPSTD_Maxwell.__flaginputs__,kw)
        yf=self.yf
        nx = np.max([1,yf.nx])
        ny = np.max([1,yf.ny])
        nz = np.max([1,yf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        j = 1j

        self.add_fields({"bx":yf.Bx,"by":yf.By,"bz":yf.Bz, \
                         "ex":yf.Ex,"ey":yf.Ey,"ez":yf.Ez})
        self.add_fields({"rho":yf.Rho},True)

        if self.l_pushf:
            self.add_fields({"f":yf.F})
            self.add_fields({"rhoold":yf.Rhoold},True)
            self.add_fields({"rhonew":yf.Rho},True)
            self.add_fields({"drho":yf.Rho},True)
        if self.l_pushg:
            self.add_fields({"g":yf.G})
#        if self.l_pushf or self.l_getrho:
#            self.add_fields({"rho":yf.Rho},True)
        self.add_fields({"rho":yf.Rho},True)
        self.add_fields({"drho":yf.Rhoold},True)
        self.add_fields({"jx":yf.Jx,"jy":yf.Jy,"jz":yf.Jz},True)

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt/self.ntsub
        cdt=dt*self.clight

        if self.nx>1:
            axm = j*dt*self.clight*self.kxm
            axp = j*dt*self.clight*self.kxp
        else:
            axm = axp = 0.

        if self.ny>1:
            aym = j*dt*self.clight*self.kym
            ayp = j*dt*self.clight*self.kyp
        else:
            aym = ayp = 0.

        if self.nz>1:
            azm = j*dt*self.clight*self.kzm
            azp = j*dt*self.clight*self.kzp
        else:
            azm = azp = 0.

        if self.nx>1:
            axp0 = 0.5/self.ntsub
            axm0 = 0.65/self.ntsub
        else:
            axm0 = axp0 = 0.

        if self.ny>1:
            ayp0 = 0.55/self.ntsub
            aym0 = 0.45/self.ntsub
        else:
            aym0 = ayp0 = 0.

        if self.nz>1:
            azp0 = 0.35/self.ntsub
            azm0 = 0.25/self.ntsub
        else:
            azm0 = azp0 = 0.

        self.mymatref = self.getmaxwellmat(axp0,ayp0,azp0,axm0,aym0,azm0, \
                            0.1/self.ntsub,0.11/self.ntsub,m0,m1, \
                            0.5*self.dx,0.5*self.dy,0.5*self.dz,l_matref=1)

        matcompress = getmatcompress(self.mymatref)

        self.mymatref = exp_by_squaring_matrixlist(self.mymatref, self.ntsub, matcompress=matcompress)
        if np.all(self.V_galilean==0.):
            self.mymat = self.getmaxwellmat(axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,\
                         self.kx_unmod,self.ky_unmod,self.kz_unmod,l_matref=0,matcompress=matcompress)
        else:
            self.mymat = self.getmaxwellmat_galilean(self.kxpn,self.kypn,self.kzpn,\
                        self.kxmn,self.kymn,self.kzmn,dt,cdt,self.V_galilean)

        self.create_fortran_matrix_blocks()

    def getmaxwellmat(self,axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,
                      kx_unmod,ky_unmod,kz_unmod,l_matref=0,
                      matcompress=None):
        c=self.clight

        if self.l_pushf:
            matpushrho = GPSTD_Matrix(self.fields)
#            matpushrho.add_op('rho',{'rho':1.,'jx':-axm/c,'jy':-aym/c,'jz':-azm/c})
            matpushrho.add_op('rho',{'rho':1.,'drho':1./self.ntsub})

        matpushb = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            matpushb.add_op('bx',{'bx':1.,'ey': azp/(2*c),'ez':-ayp/(2*c),'g':axm/(2*c)})
            matpushb.add_op('by',{'by':1.,'ex':-azp/(2*c),'ez': axp/(2*c),'g':aym/(2*c)})
            matpushb.add_op('bz',{'bz':1.,'ex': ayp/(2*c),'ey':-axp/(2*c),'g':azm/(2*c)})
        else:
            matpushb.add_op('bx',{'bx':1.,'ey': azp/(2*c),'ez':-ayp/(2*c)})
            matpushb.add_op('by',{'by':1.,'ex':-azp/(2*c),'ez': axp/(2*c)})
            matpushb.add_op('bz',{'bz':1.,'ex': ayp/(2*c),'ey':-axp/(2*c)})

        matpushe = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            matpushe.add_op('ex',{'ex':1.,'by':-azm*c,'bz': aym*c,'f':axp,'jx':-dt/self.eps0})
            matpushe.add_op('ey',{'ey':1.,'bx': azm*c,'bz':-axm*c,'f':ayp,'jy':-dt/self.eps0})
            matpushe.add_op('ez',{'ez':1.,'bx':-aym*c,'by': axm*c,'f':azp,'jz':-dt/self.eps0})
        else:
            matpushe.add_op('ex',{'ex':1.,'by':-azm*c,'bz': aym*c,'jx':-dt/self.eps0})
            matpushe.add_op('ey',{'ey':1.,'bx': azm*c,'bz':-axm*c,'jy':-dt/self.eps0})
            matpushe.add_op('ez',{'ez':1.,'bx':-aym*c,'by': axm*c,'jz':-dt/self.eps0})

        if self.l_pushf:
            matpushf = GPSTD_Matrix(self.fields)
            matpushf.add_op('f',{'f':1.,'ex':axm/2,'ey':aym/2,'ez':azm/2,'rho':-0.5*cdt/self.eps0})

        if self.l_pushg:
            matpushg = GPSTD_Matrix(self.fields)
            matpushg.add_op('g',{'g':1.,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        if self.l_pushf:
            mymat_init = GPSTD_Matrix(self.fields)
            mymat_init.add_op('rho',{'rho':0.,'rhoold':1.})
            mymat_init.add_op('drho',{'drho':0.,'rhonew':1.,'rhoold':-1.})

        if self.l_pushf:
            mymat = multmat(matpushf.mat,matpushb.mat,matcompress=matcompress)
            if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushe.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushrho.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushf.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)
        else:
            mymat = multmat(matpushb.mat,matpushe.mat,matcompress=matcompress)
            if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)

        self.mymat = exp_by_squaring_matrixlist(mymat, self.ntsub, matcompress=matcompress)
        if self.l_pushf:
            self.mymat = multmat(mymat_init.mat,self.mymat)

        if l_matref:
            self.matpushb=matpushb
            self.matpushe=matpushe
            if self.l_pushf:self.matpushf=matpushf
            if self.l_pushg:self.matpushg=matpushg
            if self.l_pushf:self.matpushrho=matpushrho

        return self.mymat

    def getmaxwellmat_galilean(self,axp,ayp,azp,axm,aym,azm,dt,cdt,m0,m1,
                      kx_unmod,ky_unmod,kz_unmod,l_matref=0,
                      matcompress=None,V_galilean=[0.,0.,0.]):

        # --- equivalent to J constant in lab frame grid when 1
        # --- equivalent to J constant in moving frame grid when 0
        l_matpushj = 0

        j=1j
        c=self.clight
        V0  = np.sqrt(V_galilean[0]*V_galilean[0]+V_galilean[1]*V_galilean[1]+V_galilean[2]*V_galilean[2])
        w=self.k*c
        kV=self.kx_unmod*V_galilean[0]+self.ky_unmod*V_galilean[1]+self.kz_unmod*V_galilean[2]
        T=np.exp(j*kV*dt)
        T2=np.exp(j*kV*dt/2)

        Theta=np.exp(j*kV*self.dt)
        coef = -j*self.divsetorig(j*kV,1.-Theta,-1./self.dt)

        is_singular=(w!=0) & (kV==0)
        coef[is_singular] = j/self.dt
        self.cr=coef.copy()
        self.cro=coef.copy()*Theta
        coef /= self.kmag
        self.JxCorRhomult = coef*self.kxpn
        self.JyCorRhomult = coef*self.kypn
        self.JzCorRhomult = coef*self.kzpn

        self.JxCorRhooldmult = coef*self.kxpn*Theta
        self.JyCorRhooldmult = coef*self.kypn*Theta
        self.JzCorRhooldmult = coef*self.kzpn*Theta

        axp = axp
        ayp = ayp
        azp = azp
        axm = axm
        aym = aym
        azm = azm

        if self.l_pushf:
            matpushrho = GPSTD_Matrix(self.fields)
    #            matpushrho.add_op('rho',{'rho':T,'jx':-axm/c,'jy':-aym/c,'jz':-azm/c})
    #            matpushrho.add_op('rho',{'rho':T,'drho':1./self.ntsub})
            alpha = 0.5*j*kV*dt
            matpushrho.add_op('rho',{'rho':(1.+alpha)/(1.-alpha),'drho':1./(self.ntsub*(1.-alpha))})
            matpushdrho = GPSTD_Matrix(self.fields)
            matpushdrho.add_op('drho',{'drho':T})

        matpushb = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            matpushb.add_op('bx',{'bx':T2,'ey': azp*T2/(2*c),'ez':-ayp*T2/(2*c),'g':axm*T2/(2*c)})
            matpushb.add_op('by',{'by':T2,'ex':-azp*T2/(2*c),'ez': axp*T2/(2*c),'g':aym*T2/(2*c)})
            matpushb.add_op('bz',{'bz':T2,'ex': ayp*T2/(2*c),'ey':-axp*T2/(2*c),'g':azm*T2/(2*c)})
        else:
            matpushb.add_op('bx',{'bx':T2,'ey': azp*T2/(2*c),'ez':-ayp*T2/(2*c)})
            matpushb.add_op('by',{'by':T2,'ex':-azp*T2/(2*c),'ez': axp*T2/(2*c)})
            matpushb.add_op('bz',{'bz':T2,'ex': ayp*T2/(2*c),'ey':-axp*T2/(2*c)})

        matpushe = GPSTD_Matrix(self.fields)
        if self.l_pushf:
            matpushe.add_op('ex',{'ex':T,'by':-azm*T2*c,'bz': aym*T2*c,'f':axp*T2,'jx':-dt*T/self.eps0})
            matpushe.add_op('ey',{'ey':T,'bx': azm*T2*c,'bz':-axm*T2*c,'f':ayp*T2,'jy':-dt*T/self.eps0})
            matpushe.add_op('ez',{'ez':T,'bx':-aym*T2*c,'by': axm*T2*c,'f':azp*T2,'jz':-dt*T/self.eps0})
        else:
            matpushe.add_op('ex',{'ex':T,'by':-azm*T2*c,'bz': aym*T2*c,'jx':-dt*T/self.eps0})
            matpushe.add_op('ey',{'ey':T,'bx': azm*T2*c,'bz':-axm*T2*c,'jy':-dt*T/self.eps0})
            matpushe.add_op('ez',{'ez':T,'bx':-aym*T2*c,'by': axm*T2*c,'jz':-dt*T/self.eps0})

        if self.l_pushf:
            matpushf = GPSTD_Matrix(self.fields)
            matpushf.add_op('f',{'f':T2,'ex':axm*T2/2,'ey':aym*T2/2,'ez':azm*T2/2,'rho':-0.5*cdt/self.eps0})

        if self.l_pushg:
            matpushg = GPSTD_Matrix(self.fields)
            matpushg.add_op('g',{'g':T,'bx':axp*T2*c,'by':ayp*T2*c,'bz':azp*T2*c})

        mymat_init = GPSTD_Matrix(self.fields)
        if l_matpushj:
            shift_init = np.exp((-0.5+0*0.5/self.ntsub)*j*kV*self.dt)
        else:
            shift_init = 1.
        mymat_init.add_op('jx',{'jx':shift_init})
        mymat_init.add_op('jy',{'jy':shift_init})
        mymat_init.add_op('jz',{'jz':shift_init})
        if self.l_pushf:
            mymat_init.add_op('rho',{'rhoold':1.})
            if l_matpushj:
                mymat_init.add_op('drho',{'rhonew':np.exp(-j*kV*self.dt),'rhoold':-1.})
            else:
                mymat_init.add_op('drho',{'rhonew':np.exp(-j*kV*self.dt),'rhoold':-1.})

        matpushj = GPSTD_Matrix(self.fields)
        matpushj.add_op('jx',{'jx':T})
        matpushj.add_op('jy',{'jy':T})
        matpushj.add_op('jz',{'jz':T})

        if self.l_pushf:
            mymat = multmat(matpushf.mat,matpushb.mat,matcompress=matcompress)
            mymat = multmat(mymat,matpushe.mat,matcompress=matcompress)
        else:
            mymat = multmat(matpushb.mat,matpushe.mat,matcompress=matcompress)
        if self.l_pushg:mymat = multmat(mymat,matpushg.mat,matcompress=matcompress)
        if self.l_pushf:mymat = multmat(mymat,matpushdrho.mat,matcompress=matcompress)
        if self.l_pushf:mymat = multmat(mymat,matpushrho.mat,matcompress=matcompress)
        if self.l_pushf:mymat = multmat(mymat,matpushf.mat,matcompress=matcompress)
        mymat = multmat(mymat,matpushb.mat,matcompress=matcompress)

        if l_matpushj:mymat = multmat(mymat,matpushj.mat,matcompress=matcompress)

        mymat = exp_by_squaring_matrixlist(mymat, self.ntsub, matcompress=matcompress)

        self.mymat = multmat(mymat_init.mat,mymat)

        if l_matref:
            self.matpushb=matpushb
            self.matpushe=matpushe
            if self.l_pushf:self.matpushf=matpushf
            if self.l_pushg:self.matpushg=matpushg
            if self.l_pushf:self.matpushrho=matpushrho

        return self.mymat

    def push(self):

        self.push_fields()

        return

class PSATD_Maxwell(GPSTDPXR):

    __flaginputs__ = {'yf':None,
                      'l_pushf':False,
                      'l_pushg':False,
                      'clight':299792458.0,
                      'eps0':8.854187817620389e-12,
                      'V_galilean':np.array([0.,0.,0.]),
                      'V_pseudogalilean':np.array([0.,0.,0.])}

    def __init__(self,**kw):
        try:
            kw['kwdict'].update(kw)
            kw = kw['kwdict']
            del kw['kwdict']
        except KeyError:
            pass

        self.processdefaultsfromdict(PSATD_Maxwell.__flaginputs__,kw)

        yf=self.yf

        nx = np.max([1,yf.nx])
        ny = np.max([1,yf.ny])
        nz = np.max([1,yf.nz])
        kw['nx']=nx
        kw['ny']=ny
        kw['nz']=nz

        GPSTD.__init__(self,kwdict=kw)

        j = 1j

        self.add_fields({"bx":yf.Bx,"by":yf.By,"bz":yf.Bz, \
                         "ex":yf.Ex,"ey":yf.Ey,"ez":yf.Ez})
        if self.l_pushf:
            self.add_fields({"f":yf.F})
        if self.l_pushg:
            self.add_fields({"g":yf.G})
        self.add_fields({"rhoold":yf.Rhoold},True)
        self.add_fields({"rhonew":yf.Rho},True)
        self.add_fields({"jx":yf.Jx,"jy":yf.Jy,"jz":yf.Jz},True)

        self.get_Ffields()

        m0 = 0.
        m1 = 1.
        dt=self.dt
        cdt=dt*self.clight
        self.wdt = self.k*cdt
        self.coswdt=np.cos(self.wdt)
        self.sinwdt=np.sin(self.wdt)
        C=self.coswdt
        S=self.sinwdt

        if np.all(self.V_galilean==0.) and np.all(self.V_pseudogalilean==0.):
            self.mymat = self.getmaxwellmat(self.kxpn,self.kypn,self.kzpn,\
                     self.kxmn,self.kymn,self.kzmn,dt,cdt)
        else:
            if np.any(self.V_galilean<>0.):
                self.mymat = self.getmaxwellmat_galilean(self.kxpn,self.kypn, \
                            self.kzpn, self.kxmn,self.kymn,self.kzmn,dt,cdt, \
                            self.V_galilean)

            if np.any(self.V_pseudogalilean<>0.):
                self.mymat = self.getmaxwellmat_pseudogalilean(self.kxpn, \
                             self.kypn, self.kzpn, self.kxmn,self.kymn, \
                             self.kzmn,dt,cdt,self.V_pseudogalilean)

        self.create_fortran_matrix_blocks()

    def getmaxwellmat(self,kxpn,kypn,kzpn,kxmn,kymn,kzmn,dt,cdt):

        j = 1j
        c=self.clight
        C=self.coswdt
        S=self.sinwdt

        Soverk = self.divsetorig(S,self.kmag,self.dt*self.clight)
        Jmult = 1./(self.kmag*self.clight*self.eps0)

        EJmult = -self.divsetorig(S,self.kmag*self.clight*self.eps0,self.dt/self.eps0)

        ERhomult = j*(-EJmult/dt-1./self.eps0)/self.kmag
        ERhooldmult = j*(C/self.eps0+EJmult/dt) /self.kmag

        BJmult = j*(C-1.)*Jmult/self.clight

        FJmult = j*(C-1.)*Jmult
        FRhomult = self.divsetorig(C-1.,dt*self.kmag**2*self.clight*self.eps0,-0.5*self.dt*self.clight/self.eps0)

        if self.nx>1:
            axm = j*S*self.kxmn
            axp = j*S*self.kxpn
            kxpn = self.kxpn
            kxmn = self.kxmn
        else:
            axm = axp = 0.
            bxm = bxp = 0.
            kxpn = kxmn = 0.

        if self.ny>1:
            aym = j*S*self.kymn
            ayp = j*S*self.kypn
            kypn = self.kypn
            kymn = self.kymn
        else:
            aym = ayp = 0.
            bym = byp = 0.
            kypn = kymn = 0.

        if self.nz>1:
            azm = j*S*self.kzmn
            azp = j*S*self.kzpn
            kzpn = self.kzpn
            kzmn = self.kzmn
        else:
            azm = azp = 0.
            bzm = bzp = 0.
            kzpn = kzmn = 0.

        self.BJmult = BJmult
        self.EJmult = EJmult
        self.ERhomult = ERhomult
        self.ERhooldmult = ERhooldmult
        self.Jmult = Jmult
        self.Soverk = Soverk

        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'g':axm/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'g':aym/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'g':azm/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
        else:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})

        if self.l_pushf:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'f':axp,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'f':ayp,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'f':azp,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})
        else:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})

        if self.l_pushf:
            mymat.add_op('f',{'f':C,'ex':axm,'ey':aym,'ez':azm, \
                                    'jx': kxmn*FJmult,'jy': kymn*FJmult,'jz': kzmn*FJmult, \
                                    'rhonew':FRhomult,\
                                    'rhoold':-FRhomult - Soverk/self.eps0})

        if self.l_pushg:
            mymat.add_op('g',{'g':C,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        return mymat.mat
    def getmaxwellmat_galilean(self,kxpn,kypn,kzpn,kxmn,kymn,kzmn,dt,cdt,V_galilean=np.array([0.,0.,0.])):

        j = 1j
        V0 = np.linalg.norm(V_galilean)
        c=self.clight
        C=self.coswdt
        S=self.sinwdt
        kV=self.kx_unmod*V_galilean[0]+self.ky_unmod*V_galilean[1]+self.kz_unmod*V_galilean[2]
        Theta=T=np.exp(j*kV*dt)
        CT = C*Theta
        ST = S*Theta
        w = self.k*c
        kVow = self.divsetorig(kV,self.kmag*c,V0/c)
        So1mT = self.divsetorig(S,1.-T,j*c/V0)
        onemCo1mT = self.divsetorig(1.-C,1.-T,0.)

        denom = (w*w-kV*kV)
        self.denom=denom

        X1 = self.divsetorig(1.-CT+j*kVow*ST, denom, dt**2*0.5)
        X2 = self.divsetorig(1.+j*kVow*T*So1mT+kVow**2*T*onemCo1mT, denom, dt**2/6)
        X3 = T*self.divsetorig(C+j*kVow*T*So1mT+kVow**2*onemCo1mT, denom, -dt**2/3)

        #Apply a special limit when kV=0 but w!=0
        is_singular=(w!=0) & (kV==0)
        X2[is_singular]=(1.-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2
        X3[is_singular] = T[is_singular]*(C[is_singular]-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2

        Soverk = self.divsetorig(S,self.kmag,self.dt*self.clight)
        Jmult = 1./(self.kmag*self.clight*self.eps0)

        EJmult = -self.divsetorig(ST,self.kmag*self.clight*self.eps0,dt/self.eps0)+j*X1*kV/self.eps0

        ERhomult = -j*c**2*X2*self.k/self.eps0
        ERhooldmult = j*c**2*X3*self.k/self.eps0

        BJmult = -self.k*j*X1/self.eps0

        FJmult = j*(C-1.)*Jmult
        FRhomult = (C-1.)/(dt*self.kmag**2*self.clight*self.eps0)

        coef = -j*self.divsetorig(j*kV,1.-T,-1./dt)
        coef[is_singular]=j/dt
        self.CDcoef = coef.copy()*np.exp(0.5*j*kV*dt)

        self.cr=coef.copy()
        self.cro=coef.copy()*T
        coef /= self.kmag
        self.JxCorRhomult = coef*self.kxpn
        self.JyCorRhomult = coef*self.kypn
        self.JzCorRhomult = coef*self.kzpn

        self.JxCorRhooldmult = coef*self.kxpn*T
        self.JyCorRhooldmult = coef*self.kypn*T
        self.JzCorRhooldmult = coef*self.kzpn*T

        if len(self.dims)==1:
            FRhomult[0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==2:
            FRhomult[0,0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==3:
            FRhomult[0,0,0] = -0.5*self.dt*self.clight/self.eps0

        if self.nx>1:
            axm = j*ST*self.kxmn
            axp = j*ST*self.kxpn
            kxpn = self.kxpn
            kxmn = self.kxmn
        else:
            axm = axp = 0.
            bxm = bxp = 0.
            kxpn = kxmn = 0.

        if self.ny>1:
            aym = j*ST*self.kymn
            ayp = j*ST*self.kypn
            kypn = self.kypn
            kymn = self.kymn
        else:
            aym = ayp = 0.
            bym = byp = 0.
            kypn = kymn = 0.

        if self.nz>1:
            azm = j*ST*self.kzmn
            azp = j*ST*self.kzpn
            kzpn = self.kzpn
            kzmn = self.kzmn
        else:
            azm = azp = 0.
            bzm = bzp = 0.
            kzpn = kzmn = 0.

        self.BJmult = BJmult
        self.EJmult = EJmult
        self.ERhomult = ERhomult
        self.ERhooldmult = ERhooldmult
        self.Jmult = Jmult
        self.Soverk = Soverk
        self.X1=X1
        self.X2=X2
        self.X3=X3
        self.kVow = kVow
        self.So1mT = So1mT
        self.onemCo1mT = onemCo1mT
        self.kV=kV
        self.T=T
        self.CT = CT
        self.ST = ST

        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            mymat.add_op('bx',{'bx':CT,'ey': azp/c,'ez':-ayp/c,'g':axm/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':CT,'ex':-azp/c,'ez': axp/c,'g':aym/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':CT,'ex': ayp/c,'ey':-axp/c,'g':azm/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
        else:
            mymat.add_op('bx',{'bx':CT,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':CT,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':CT,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})

        if self.l_pushf:
            mymat.add_op('ex',{'ex':CT,'by':-azm*c,'bz': aym*c,'jx':EJmult,'f':axp,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':CT,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'f':ayp,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':CT,'bx':-aym*c,'by': axm*c,'jz':EJmult,'f':azp,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})
        else:
            mymat.add_op('ex',{'ex':CT,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':CT,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':CT,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})

        if self.l_pushf:
            print 'l_pushf not yet implemented in PSATD Galilean'
            raise
            mymat.add_op('f',{'f':CT,'ex':axm,'ey':aym,'ez':azm, \
                                    'jx': kxmn*FJmult,'jy': kymn*FJmult,'jz': kzmn*FJmult, \
                                    'rhonew':FRhomult,\
                                    'rhoold':-FRhomult - Soverk/self.eps0})

        if self.l_pushg:
            print 'l_pushg not yet implemented in PSATD Galilean'
            raise
            mymat.add_op('g',{'g':CT,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        return mymat.mat

    def getmaxwellmat_pseudogalilean(self,kxpn,kypn,kzpn,kxmn,kymn,kzmn,dt,cdt,V_galilean=np.array([0.,0.,0.])):

        j = 1j
        V0 = np.sqrt(V_galilean[0]*V_galilean[0]+V_galilean[1]*V_galilean[1]+V_galilean[2]*V_galilean[2])
        c=self.clight
        C=self.coswdt
        S=self.sinwdt
        kV=self.kx_unmod*V_galilean[0]+self.ky_unmod*V_galilean[1]+self.kz_unmod*V_galilean[2]
        Theta=T=np.exp(j*kV*dt)
        invT = np.exp(-j*kV*dt)
        CT = C*Theta
        ST = S*Theta
        w = self.k*c
        kVow = self.divsetorig(kV,self.kmag*c,V0/c)
        So1mT = self.divsetorig(S,1.-T,j*c/V0)
        onemCo1mT = self.divsetorig(1.-C,1.-T,0.)

        denom = (w*w-kV*kV)
        self.denom=denom

        X1 = self.divsetorig(1.-CT+j*kVow*ST, denom, dt**2*0.5)
        X2 = self.divsetorig(1.+j*kVow*T*So1mT+kVow**2*T*onemCo1mT, denom, dt**2/6)
        X3 = self.divsetorig(C+j*kVow*T*So1mT+kVow**2*onemCo1mT, denom, -dt**2/3)

        #Apply a special limit when kV=0 but w!=0
        is_singular=(w!=0.) & (kV==0.)
        X2[is_singular]=(1.-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2
        X3[is_singular] = (C[is_singular]-S[is_singular]/(w[is_singular]*dt))/w[is_singular]**2

        X1 *= np.exp(-0.5*j*kV*dt)

        Soverk = self.divsetorig(S,self.kmag,self.dt*self.clight)
        Jmult = 1./(self.kmag*self.clight*self.eps0)

        EJmult = -np.exp(-0.5*j*kV*dt)*self.divsetorig(ST,self.kmag*self.clight*self.eps0,dt/self.eps0)+j*X1*kV/self.eps0

        ERhomult = -j*c**2*X2*self.k/self.eps0
        ERhooldmult = j*c**2*X3*self.k/self.eps0

        BJmult = -self.k*j*X1/self.eps0

        FJmult = j*(C-1.)*Jmult
        FRhomult = (C-1.)/(dt*self.kmag**2*self.clight*self.eps0)

        coef = self.divsetorig(j*kV*dt,T-1.,1.)
        coef[is_singular]=1.

        self.CDcoef = coef.copy()*np.exp(0.5*j*kV*dt)
        coef *= j/(dt*self.kmag)

        self.JxCorRhomult = coef*self.kxpn*np.exp(0.5*j*kV*dt)
        self.JyCorRhomult = coef*self.kypn*np.exp(0.5*j*kV*dt)
        self.JzCorRhomult = coef*self.kzpn*np.exp(0.5*j*kV*dt)

        self.JxCorRhooldmult = coef*self.kxpn*np.exp(0.5*j*kV*dt)
        self.JyCorRhooldmult = coef*self.kypn*np.exp(0.5*j*kV*dt)
        self.JzCorRhooldmult = coef*self.kzpn*np.exp(0.5*j*kV*dt)

        if len(self.dims)==1:
            FRhomult[0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==2:
            FRhomult[0,0] = -0.5*self.dt*self.clight/self.eps0
        if len(self.dims)==3:
            FRhomult[0,0,0] = -0.5*self.dt*self.clight/self.eps0

        if self.nx>1:
            axm = j*S*self.kxmn
            axp = j*S*self.kxpn
            kxpn = self.kxpn
            kxmn = self.kxmn
        else:
            axm = axp = 0.
            bxm = bxp = 0.
            kxpn = kxmn = 0.

        if self.ny>1:
            aym = j*S*self.kymn
            ayp = j*S*self.kypn
            kypn = self.kypn
            kymn = self.kymn
        else:
            aym = ayp = 0.
            bym = byp = 0.
            kypn = kymn = 0.

        if self.nz>1:
            azm = j*S*self.kzmn
            azp = j*S*self.kzpn
            kzpn = self.kzpn
            kzmn = self.kzmn
        else:
            azm = azp = 0.
            bzm = bzp = 0.
            kzpn = kzmn = 0.

        self.BJmult = BJmult
        self.EJmult = EJmult
        self.ERhomult = ERhomult
        self.ERhooldmult = ERhooldmult
        self.Jmult = Jmult
        self.Soverk = Soverk
        self.X1=X1
        self.X2=X2
        self.X3=X3
        self.kVow = kVow
        self.So1mT = So1mT
        self.onemCo1mT = onemCo1mT
        self.kV=kV
        self.T=T
        self.CT = CT
        self.ST = ST

        mymat = GPSTD_Matrix(self.fields)
        if self.l_pushg:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'g':axm/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'g':aym/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'g':azm/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
        else:
            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})

        if self.l_pushf:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'f':axp,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'f':ayp,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'f':azp,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})
        else:
            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})

        if self.l_pushf:
            print 'l_pushf not yet implemented in PSATD PseudoGalilean'
            raise
            mymat.add_op('f',{'f':C,'ex':axm,'ey':aym,'ez':azm, \
                                    'jx': kxmn*FJmult,'jy': kymn*FJmult,'jz': kzmn*FJmult, \
                                    'rhonew':FRhomult,\
                                    'rhoold':-FRhomult - Soverk/self.eps0})

        if self.l_pushg:
            mymat.add_op('g',{'g':C,'bx':axp*c,'by':ayp*c,'bz':azp*c})

        return mymat.mat

    def push(self):

        self.push_fields()

        return

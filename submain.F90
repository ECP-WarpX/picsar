!===============================================================================
! TIME LOOP OF THE PIC ALGORITHM
SUBROUTINE step(nst)
!===============================================================================
USE constants
USE fields
USE particles
USE params
USE shared_data
USE boundary
USE omp_lib
USE diagnostics
USE simple_io

IMPLICIT NONE
INTEGER :: nst,i, k
INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0
INTEGER :: err
CHARACTER(LEN=string_length) :: strtemp
INTEGER :: ixsource, iysource, izsource
REAL(num) :: xsource, ysource, zsource
REAL(num), DIMENSION(-nxguards:nx+nxguards, -nyguards:ny+nyguards, -nzguards:nz+nzguards) :: dive

DO i=1,nst

!!!! --- Input external field in the simulation box
!xsource=xmin+length_x/2.0
!ysource=ymin+length_y/2.0
!zsource=zmin+length_z/2.0
!IF (((xsource .GT. x_min_local) .AND.(xsource .LT. x_max_local)) .AND. &
!((ysource .GT. y_min_local) .AND.(ysource .LT. y_max_local)) .AND. &
!((zsource .GT. z_min_local) .AND.(zsource .LT. z_max_local)) .AND. (it .EQ. 2)) THEN
!ixsource =NINT((xsource-x_min_local)/(x_max_local-x_min_local)*nx)
!iysource =NINT((ysource-y_min_local)/(y_max_local-y_min_local)*ny)
!izsource =NINT((zsource-z_min_local)/(z_max_local-z_min_local)*nz)
!ey(ixsource,iysource,izsource) = 0.0_num
!!CALL smooth3d_121(ey,nx,ny,nz,npass,alpha)
!ENDIF

!!! --- Push velocity with B half step
CALL bpush_v(npartp,uxp(1:npartp),uyp(1:npartp),uzp(1:npartp),        &
bxp(1:npartp),byp(1:npartp),bzp(1:npartp),echarge,pmass,dt*0.5_num)
CALL bpush_v(nparte,uxe(1:nparte),uye(1:nparte),uze(1:nparte),        &
bxe(1:nparte),bye(1:nparte),bze(1:nparte),-echarge,emass,dt*0.5_num)

!!! --- Push velocity with E half step
CALL epush_v(npartp,uxp(1:npartp),uyp(1:npartp),uzp(1:npartp),          &
exp(1:npartp),eyp(1:npartp),ezp(1:npartp),echarge,pmass,dt*0.5_num)
CALL epush_v(nparte,uxe(1:nparte),uye(1:nparte),uze(1:nparte),          &
exe(1:nparte),eye(1:nparte),eze(1:nparte),-echarge,emass,dt*0.5_num)


!!! --- Push B field half a time step
CALL push_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                                                  &
                           0.5_num*dt/dx*xcoeffs,0.5_num*dt/dy*ycoeffs,0.5_num*dt/dz*zcoeffs,  &
                           nx,ny,nz, norderx,nordery,norderz,                                  &
                           nxguards,nyguards,nzguards,nxs,nys,nzs,                             &
                           l_nodalgrid)

!!! --- Boundary conditions for B
CALL bfield_bcs

!!!! --- push electrons positions a time step
CALL pushxyz(nparte,xe(1:nparte),ye(1:nparte),ze(1:nparte),        &
uxe(1:nparte),uye(1:nparte),uze(1:nparte),dt)
!!! --- push protons positions a time step
CALL pushxyz(npartp,xp(1:npartp),yp(1:npartp),zp(1:npartp),        &
uxp(1:npartp),uyp(1:npartp),uzp(1:npartp),dt)


!!! --- Apply BC on particles
CALL particle_bcs_e
CALL particle_bcs_p

jx = 0.0_num
jy = 0.0_num
jz = 0.0_num


!!! --- Deposit current for electrons
CALL depose_jxjyjz_esirkepov_n(jx,jy,jz,nparte,xe(1:nparte),ye(1:nparte),ze(1:nparte),                &
                               uxe(1:nparte),uye(1:nparte),uze(1:nparte),                             &
                               we(1:nparte),-echarge,x_min_local,y_min_local,z_min_local,dt,dx,dy,dz, &
                               nx,ny,nz,nxguards,nyguards,nzguards,                                   &
                               nox,noy,noz,l_particles_weight,l4symtry)
!!! --- Deposit current for protons
CALL depose_jxjyjz_esirkepov_n(jx,jy,jz,npartp,xp(1:npartp),yp(1:npartp),zp(1:npartp),                &
                               uxp(1:npartp),uyp(1:npartp),uzp(1:npartp),                             &
                               wp(1:npartp),echarge,x_min_local,y_min_local,z_min_local,dt,dx,dy,dz,  &
                               nx,ny,nz,nxguards,nyguards,nzguards,                                   &
                               nox,noy,noz,l_particles_weight,l4symtry)


!!! --- Boundary conditions for currents
CALL current_bcs


!!! --- push E field  a full time step
CALL push_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,clight**2*mu0*dt,    &
                          clight**2*dt/dx*xcoeffs,clight**2*dt/dy*ycoeffs, &
                          clight**2*dt/dz*zcoeffs,nx,ny,nz,                &
                          norderx,nordery,norderz,                         &
                          nxguards,nyguards,nzguards,nxs,nys,nzs,          &
                          l_nodalgrid)

!!! --- Boundary conditions for E
CALL efield_bcs

! - Compute electric field divergence on grid at n+1
dive=0.0_num
CALL calc_field_div(dive, ex, ey, ez, nx, ny, nz, nxguards, nyguards, nzguards, dx, dy, dz)


!!! --- push B field half a time step
CALL push_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                                                  &
                           0.5_num*dt/dx*xcoeffs,0.5_num*dt/dy*ycoeffs,0.5_num*dt/dz*zcoeffs,  &
                           nx,ny,nz, norderx,nordery,norderz,                                  &
                           nxguards,nyguards,nzguards,nxs,nys,nzs,                             &
                           l_nodalgrid)

!!! --- Boundary conditions for B
CALL bfield_bcs

exsm = ex
eysm = ey
ezsm = ez
bxsm = bx
bysm = by
bzsm = bz

!!! --- Gather electric fields on protons
exp=0.0_num;eyp=0.0_num;ezp=0.0_num
CALL gete3d_n_energy_conserving(npartp,xp(1:npartp),yp(1:npartp),zp(1:npartp),             &
                                      exp(1:npartp),eyp(1:npartp),ezp(1:npartp),           &
                                      x_min_local,y_min_local,z_min_local,                 &
                                      dx,dy,dz,nx,ny,nz,nxguards,nyguards,nzguards,        &
                                      nox,noy,noz,exsm,eysm,ezsm,l_lower_order_in_v)
!!! --- Gather electric fields on electrons
exe=0.0_num;eye=0.0_num;eze=0.0_num
CALL gete3d_n_energy_conserving(nparte,xe(1:nparte),ye(1:nparte),ze(1:nparte),             &
                                      exe(1:nparte),eye(1:nparte),eze(1:nparte),           &
                                      x_min_local,y_min_local,z_min_local,                 &
                                      dx,dy,dz,nx,ny,nz,nxguards,nyguards,nzguards,        &
                                      nox,noy,noz,exsm,eysm,ezsm,l_lower_order_in_v)
!!! --- Gather magnetic fields on protons
bxp=0.0_num;byp=0.0_num;bzp=0.0_num
CALL getb3d_n_energy_conserving(npartp,xp(1:npartp),yp(1:npartp),zp(1:npartp),             &
                                      bxp(1:npartp),byp(1:npartp),bzp(1:npartp),           &
                                      x_min_local,y_min_local,z_min_local,                 &
                                      dx,dy,dz,nx,ny,nz,nxguards,nyguards,nzguards,        &
                                      nox,noy,noz,bxsm,bysm,bzsm,l_lower_order_in_v)
!!! --- Gather magnetic fields on electrons
bxe=0.0_num;bye=0.0_num;bze=0.0_num
CALL getb3d_n_energy_conserving(nparte,xe(1:nparte),ye(1:nparte),ze(1:nparte),             &
                                      bxe(1:nparte),bye(1:nparte),bze(1:nparte),           &
                                      x_min_local,y_min_local,z_min_local,                 &
                                      dx,dy,dz,nx,ny,nz,nxguards,nyguards,nzguards,        &
                                      nox,noy,noz,bxsm,bysm,bzsm,l_lower_order_in_v)


!!! --- Push velocity with E half step
CALL epush_v(npartp,uxp(1:npartp),uyp(1:npartp),uzp(1:npartp),        &
             exp(1:npartp),eyp(1:npartp),ezp(1:npartp),echarge,pmass,dt*0.5_num)
CALL epush_v(nparte,uxe(1:nparte),uye(1:nparte),uze(1:nparte),        &
            exe(1:nparte),eye(1:nparte),eze(1:nparte),-echarge,emass,dt*0.5_num)

!!! --- Push velocity with B full step
CALL bpush_v(npartp,uxp(1:npartp),uyp(1:npartp),uzp(1:npartp),        &
            bxp(1:npartp),byp(1:npartp),bzp(1:npartp),echarge,pmass,dt*0.5_num)
CALL bpush_v(nparte,uxe(1:nparte),uye(1:nparte),uze(1:nparte),        &
            bxe(1:nparte),bye(1:nparte),bze(1:nparte),-echarge,emass,dt*0.5_num)


!!! --- Compute charge density on grid at n+1
!! - Compute proton density on grid
rho=0.0_num
CALL depose_rho_n(rho, npartp,xp(1:npartp),yp(1:npartp),zp(1:npartp),wp(1:npartp),  &
echarge,x_min_local,y_min_local,z_min_local,dx,dy,dz,nx,ny,nz,nxguards,    &
nyguards,nzguards,nox,noy,noz,l_particles_weight,l4symtry)

!! - Compute electron density on grid
CALL depose_rho_n(rho, nparte,xe(1:nparte),ye(1:nparte),ze(1:nparte),we(1:nparte),  &
-echarge,x_min_local,y_min_local,z_min_local,dx,dy,dz,nx,ny,nz,nxguards,   &
nyguards,nzguards,nox,noy,noz,l_particles_weight,l4symtry)

CALL charge_bcs

!!! --- Write output to disk
WRITE(strtemp,'(I5)') it
!! -- Write grid quantities
! - Write current density ex
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileex))// &
TRIM(ADJUSTL(strtemp))//'.pxr', ex, offset, err)
! - Write current density jx
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileey))// &
TRIM(ADJUSTL(strtemp))//'.pxr', ey, offset, err)
! - Write current density jx
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileez))// &
TRIM(ADJUSTL(strtemp))//'.pxr', ez, offset, err)
! - Write current density jx
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filebx))// &
TRIM(ADJUSTL(strtemp))//'.pxr', bx, offset, err)
! - Write current density jx
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileby))// &
TRIM(ADJUSTL(strtemp))//'.pxr', by, offset, err)
! - Write current density jx
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filebz))// &
TRIM(ADJUSTL(strtemp))//'.pxr', bz, offset, err)
! - Write current density jx
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filejx))// &
TRIM(ADJUSTL(strtemp))//'.pxr', jx, offset, err)
! - Write electric field divergence div E
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filedive))// &
TRIM(ADJUSTL(strtemp))//'.pxr', dive, offset, err)
! - Write total charge density rho
CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filerho))// &
TRIM(ADJUSTL(strtemp))//'.pxr', rho, offset, err)

it = it+1

IF (rank .EQ. 0) THEN
	write(0,*) 'it = ',it,'time = ',it*dt
END IF

END DO

END SUBROUTINE step

!===============================================================================
! ALLOCATE/INIT PLASMA AND FIELD ARRAYS
SUBROUTINE initall
!===============================================================================
USE constants
USE params
USE fields
USE particles
USE shared_data
!use IFPORT ! uncomment if using the intel compiler (for rand)
IMPLICIT NONE
INTEGER:: i,ierror,j,k,l
INTEGER :: npartemp
INTEGER :: ixsource, iysource, izsource
REAL(num) :: xsource, ysource, zsource
!real(8) :: rand

!!! --- Set time step
dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))

!!! --- Set number of particles
!nparte = nx*ny*nz*nppcell
!npartp = nx*ny*nz*nppcell
nparte=0
npartp=0
IF ((y_max_local .GT. ymax/2.0_num) .AND. (y_min_local .LT. ymax/2.0_num)) THEN
    nparte=1
    npartp=1
    npartemp=2*nppcell*nx*ny*nz ! large fixed table size to avoid frequent alloc-de-alloc to increase size at runtime
ELSE
    nparte=0
    npartp=0
    npartemp=2*nppcell*nx*ny*nz
END IF

!!! --- Allocate particle arrays
IF (.NOT. l_arrays_allocated) THEN
  ! --- allocate field arrays -> See MPI_ROUTINES

  ! --- allocate electron parameters
  ALLOCATE(xe(npartemp),ye(npartemp),ze(npartemp),uxe(npartemp),uye(npartemp),uze(npartemp),we(npartemp), &
           exe(npartemp),eye(npartemp),eze(npartemp),bxe(npartemp),bye(npartemp),bze(npartemp))
  ! --- allocate proton parameters
  ALLOCATE(xp(npartemp),yp(npartemp),zp(npartemp),uxp(npartemp),uyp(npartemp),uzp(npartemp),wp(npartemp), &
           exp(npartemp),eyp(npartemp),ezp(npartemp),bxp(npartemp),byp(npartemp),bzp(npartemp))
  l_arrays_allocated=.true.
  ! --- allocate coefficient arrays
  ALLOCATE(xcoeffs(norderx/2),ycoeffs(nordery/2),zcoeffs(norderz/2))
END IF

!!! --- Initialize particle and field arrays
ex=0.0_num;ey=0.0_num;ez=0.0_num
bx=0.0_num;by=0.0_num;bz=0.0_num
exsm=0.0_num;eysm=0.0_num;ezsm=0.0_num
bxsm=0.0_num;bysm=0.0_num;bzsm=0.0_num
jx=0.0_num;jy=0.0_num;jz=0.0_num
exe=0.0_num;eye=0.0_num;eze=0.0_num
bxe=0.0_num;bye=0.0_num;bze=0.0_num
uxe=0.0_num;uye=0.0_num;uze=0.0_num
exp=0.0_num;eyp=0.0_num;ezp=0.0_num
bxp=0.0_num;byp=0.0_num;bzp=0.0_num
uxp=0.0_num;uyp=0.0_num;uzp=0.0_num
we = nc*dx*dy*dz/(nppcell)
wp = we
it = 0

!!! --- Initialize stencil coefficients array for Maxwell field solver
CALL FD_weights(xcoeffs, norderx, l_nodalgrid)
CALL FD_weights(ycoeffs, nordery, l_nodalgrid)
CALL FD_weights(zcoeffs, norderz, l_nodalgrid)

!!! --- Initialize particle positions  in each subdomain
! Species 1 : Protons
DO i = 1, npartp
    xp(i) = x_min_local+(x_max_local-x_min_local)/2.0_num
    yp(i) = y_min_local+(y_max_local-y_min_local)/2.0_num
    zp(i) = z_min_local+(z_max_local-z_min_local)/2.0_num
END DO

IF ((y_max_local .GT. ymax/2.0_num) .AND. (y_min_local .LT. ymax/2.0_num)) THEN
    PRINT *, "xp", xp(1)/dx
    PRINT *, "yp", yp(1)/dy
    PRINT *, "zp", zp(1)/dz
END IF
PRINT *, "rank", rank, "x_min_local", x_min_local, "y_min_local", y_min_local, "z_min_local", z_min_local
PRINT *, "rank", rank, "x_max_local", x_max_local, "y_max_local", y_max_local, "z_max_local", z_max_local

! Species 2 : Electrons
xe = xp
ye = yp
ze = zp

!!! --- Initialize particles velocities (cold plasma here -simple case)
! Species 1 : Protons
uxp=0.0_num*clight
uyp=0.8_num*clight
uzp=0.0_num*clight
! Species 2 : Electrons
uxe=0.0_num*clight
uye=-0.8_num*clight
uze=0.0_num*clight

!!! --- Initialize external field at t=0
xsource=xmin+length_x/2.0_num
ysource=ymin+length_y/2.0_num
zsource=zmin+length_z/2.0_num

IF (((xsource .GE. x_min_local) .AND.(xsource .LE. x_max_local)) .AND. &
    ((ysource .GT. y_min_local) .AND.(ysource .LE. y_max_local)) .AND. &
    ((zsource .GE. z_min_local) .AND.(zsource .LE. z_max_local))) THEN
    ixsource =NINT((xsource-x_min_local)/(x_max_local-x_min_local)*nx)
    iysource =NINT((ysource-y_min_local)/(y_max_local-y_min_local)*ny)
    izsource =NINT((zsource-z_min_local)/(z_max_local-z_min_local)*nz)
    ey(ixsource,iysource,izsource) = 0.0_num
ENDIF

!!! --- set number of time steps
nsteps = nint(tmax/(w0_l*dt))
IF (rank .EQ. 0) THEN
    WRITE (0,*), "nsteps = ", nsteps
    WRITE (0,*), "Loaded npart= ", nparte*nprocx*nprocy*nprocz, " of species electron"
    WRITE (0,*), "Loaded npart= ", npartp*nprocx*nprocy*nprocz, " of species  proton "
END IF

!!! --- Allocate stats tables
ALLOCATE(pushb(1:nsteps))
ALLOCATE(bcs_pushb(1:nsteps))
ALLOCATE(pushe(1:nsteps))
ALLOCATE(bcs_pushe(1:nsteps))
ALLOCATE(push_part(1:nsteps))
ALLOCATE(bcs_part(1:nsteps))
ALLOCATE(cs(1:nsteps))
ALLOCATE(bcs_cs(1:nsteps))
ALLOCATE(field_gath(1:nsteps))
pushb=0.0_num
bcs_pushb=0.0_num
pushe=0.0_num
bcs_pushe=0.0_num
push_part=0.0_num
bcs_part=0.0_num
cs=0.0_num
bcs_cs=0.0_num
END SUBROUTINE initall

!===============================================================================
! COMPUTE STENCIL COEFFICIENTS FOR MAXWELL FIELD SOLVER
!===============================================================================
SUBROUTINE FD_weights(coeffs, norder, l_nodal)
!adapted from Matlab code from Fornberg (1998)
!Calculates FD weights. The parameters are:
!z   location where approximations are to be accurate,
!n   number of grid points,
!m   highest derivative that we want to find weights for
!c   array size m+1,length(x) containing (as output) in
!successive rows the weights for derivatives 0,1,...,m.
USE constants
IMPLICIT NONE
INTEGER :: norder, n, m, mn, i, j, k
LOGICAL :: l_nodal
REAL(num) :: z, fact, c1, c2, c3, c4, c5
REAL(num), INTENT(IN OUT), DIMENSION(norder/2) :: coeffs
REAL(num), ALLOCATABLE, DIMENSION(:) :: x
REAL(num), ALLOCATABLE, DIMENSION(:,:) :: c

IF (l_nodal) THEN
    z=0.0_num
    fact=1.0_num
ELSE
    z=0.5_num
    fact=0.5_num
END IF
m=1
n=norder+1

ALLOCATE(x(0:n-1))
ALLOCATE(c(0:m,0:n-1))

DO i=0, n-1
    x(i)=(i-n/2+1)*1.0_num
END DO

c=0.0_num; c1=1.0_num; c4=x(0)-z; c(0,0)=1.0_num
DO i=1, n-1
    mn=min(i+1,m+1)
    c2=1.0_num
    c5=c4
    c4=x(i)-z
    DO j=0, i-1
        c3=x(i)-x(j)
        c2=c2*c3
        IF (j .EQ. (i-1)) THEN
            DO k=1, mn-1
                c(k,i)=c1*(k*c(k-1,i-1)-c5*c(k,i-1))/c2
            END DO
            c(0,i)=-c1*c5*c(0,i-1)/c2
            DO k=1, mn-1
                c(k,j)=(c4*c(k,j)-k*c(k-1,j))/c3
            END DO
            c(0,j)=c4*c(0,j)/c3
        END IF

    END DO
    c1=c2
END DO

DO i=1, norder/2
    coeffs(i)=c(m,norder/2+i-1)
END DO
RETURN
END SUBROUTINE FD_weights

!===============================================================================
!     --- performs a smoothing sequence of n passes of a three points stencil
!     --- with coefficients [(1-alpha)/2, alpha, (1-alpha)/2] along x, y and/or z.
SUBROUTINE smooth3d_121(q,nx,ny,nz,npass,alpha)
!===============================================================================
USE constants
IMPLICIT NONE
INTEGER :: nx,ny,nz,i,j,k,l,npass(3)

REAL(num), DIMENSION(1:nx,1:ny,1:nz) :: q
REAL(num) :: temp0, temp1, alpha(3), a, b

!!! ---  x smoothing
IF (nx>0) THEN
    a = alpha(1)
    b = (1.0_num-a)/2.0_num
    DO i=1,npass(1)
        DO  l=1,nz
            DO  k=1,ny
                temp0 = q(1,k,l)
                j=1; q(j,k,l) = q(j,k,l)*(1.-b)+b*q(j+1,k,l)
                DO  j=1,nx-1
                    temp1 = q(j,k,l)
                    q(j,k,l) = a*q(j,k,l)+b*(temp0+q(j+1,k,l))
                    temp0 = temp1
                END DO
                j=nx;q(j,k,l) = q(j,k,l)*(1.0_num-b)+b*temp0
            END DO
        END DO
    END DO
END IF

!!! ---   y smoothing
IF (ny>0) THEN
    a = alpha(2)
    b = (1.0_num-a)/2.0_num
        DO i=1,npass(2)
            DO  l=1,nz
                DO  j=1,nx
                    temp0 = q(j,1,l)
                    k=1; q(j,k,l) = q(j,k,l)*(1.0_num-b)+b*q(j,k+1,l)
                    DO  k=1,ny-1
                        temp1 = q(j,k,l)
                        q(j,k,l) = a*q(j,k,l)+b*(temp0+q(j,k+1,l))
                        temp0 = temp1
                    END DO
                    k=ny;q(j,k,l) = q(j,k,l)*(1.0_num-b)+b*temp0
                END DO
            END DO
        END DO
END IF

!!! ---   z smoothing
IF (nz>0) THEN
    a = alpha(3)
    b = (1.0_num-a)/2.0_num
        DO i=1,npass(3)
            DO  k=1,ny
                DO  j=1,nx
                    temp0 = q(j,k,1)
                    l=1; q(j,k,l) = q(j,k,l)*(10_num-b)+b*q(j,k,l+1)
                    DO  l=1,nz-1
                        temp1 = q(j,k,l)
                        q(j,k,l) = a*q(j,k,l)+b*(temp0+q(j,k,l+1))
                        temp0 = temp1
                    END DO
                    l=nz;q(j,k,l) = q(j,k,l)*(1.0_num-b)+b*temp0
                END DO
            END DO
        END DO
END IF

RETURN
END SUBROUTINE smooth3d_121

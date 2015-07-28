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
INTEGER :: ixsource, iysource, izsource
REAL(num) :: xsource, ysource, zsource

!!! --- This is the PIC LOOP
DO i=1,nst

!!! --- External field
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

!!! --- Advance velocity half a time step
CALL push_particles_v

!!! --- Push B field half a time step
CALL push_bfield

!!! --- Boundary conditions for B
CALL bfield_bcs

!!! --- Push particles a full time step
CALL push_particles_xyz

!!! --- Apply BC on particles
CALL particle_bcs

!!! --- Deposit current of particle species on the grid
CALL depose_currents_on_grid_jxjyjz

!!! --- Boundary conditions for currents
CALL current_bcs

!!! --- push E field  a full time step
CALL push_efield

!!! --- Boundary conditions for E
CALL efield_bcs

!!! --- push B field half a time step
CALL push_bfield

!!! --- Boundary conditions for B
CALL bfield_bcs

!!! --- Gather electromagnetic fields from the grid to particle species
CALL gather_ebfields_on_particles

!!! --- Advance velocity half a time step
CALL push_particles_v

!!! --- derived quantitites
CALL calc_diags

!!! --- output 
CALL output_routines

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
INTEGER:: i,ierror,j,k,l, ispecies
INTEGER :: npartemp
INTEGER :: ixsource, iysource, izsource
REAL(num) :: xsource, ysource, zsource
!real(8) :: rand

!!! --- Set time step
dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))
it = 0

!!! --- Allocate arrays
IF (.NOT. l_arrays_allocated) THEN
    nppspecies_max=2*nppcell*nx*ny*nz
    !! -- Allocate species arrays
ALLOCATE(species_parray(1:nspecies_max))
    DO ispecies=1,nspecies
    ALLOCATE(species_parray(ispecies)%part_x(1:nppspecies_max),    &
             species_parray(ispecies)%part_y(1:nppspecies_max),    &
             species_parray(ispecies)%part_z(1:nppspecies_max),    &
             species_parray(ispecies)%part_ux(1:nppspecies_max),   &
             species_parray(ispecies)%part_uy(1:nppspecies_max),   &
             species_parray(ispecies)%part_uz(1:nppspecies_max),   &
             species_parray(ispecies)%part_ex(1:nppspecies_max),   &
             species_parray(ispecies)%part_ey(1:nppspecies_max),   &
             species_parray(ispecies)%part_ez(1:nppspecies_max),   &
             species_parray(ispecies)%part_bx(1:nppspecies_max),   &
             species_parray(ispecies)%part_by(1:nppspecies_max),   &
             species_parray(ispecies)%part_bz(1:nppspecies_max),   &
             species_parray(ispecies)%weight(1:nppspecies_max))
    END DO
    l_arrays_allocated=.TRUE.
    ! --- allocate Maxwell solver coefficient arrays
    ALLOCATE(xcoeffs(norderx/2),ycoeffs(nordery/2),zcoeffs(norderz/2))
END IF
!!! --- Initialize particle and field arrays
! - Init grid arrays
ex=0.0_num;ey=0.0_num;ez=0.0_num
bx=0.0_num;by=0.0_num;bz=0.0_num
exsm=0.0_num;eysm=0.0_num;ezsm=0.0_num
bxsm=0.0_num;bysm=0.0_num;bzsm=0.0_num
jx=0.0_num;jy=0.0_num;jz=0.0_num

! - Init particle species attributes/arrays (2 species here)
! - Species 1: electron
species_parray(1)%name='electron'
species_parray(1)%mass=emass
species_parray(1)%charge=-echarge
species_parray(1)%species_npart=0
species_parray(1)%part_x=0.0_num
species_parray(1)%part_y=0.0_num
species_parray(1)%part_z=0.0_num
species_parray(1)%part_ux=0.0_num
species_parray(1)%part_uy=0.0_num
species_parray(1)%part_uz=0.0_num
species_parray(1)%part_ex=0.0_num
species_parray(1)%part_ey=0.0_num
species_parray(1)%part_ez=0.0_num
species_parray(1)%part_bx=0.0_num
species_parray(1)%part_by=0.0_num
species_parray(1)%part_bz=0.0_num
species_parray(1)%weight=nc*dx*dy*dz/(nppcell)

! - Species 2: proton
species_parray(2)%name='proton'
species_parray(2)%mass=pmass
species_parray(2)%charge=echarge
species_parray(2)%species_npart=0
species_parray(2)%part_x=0.0_num
species_parray(2)%part_y=0.0_num
species_parray(2)%part_z=0.0_num
species_parray(2)%part_ux=0.0_num
species_parray(2)%part_uy=0.0_num
species_parray(2)%part_uz=0.0_num
species_parray(2)%part_ex=0.0_num
species_parray(2)%part_ey=0.0_num
species_parray(2)%part_ez=0.0_num
species_parray(2)%part_bx=0.0_num
species_parray(2)%part_by=0.0_num
species_parray(2)%part_bz=0.0_num
species_parray(2)%weight=nc*dx*dy*dz/(nppcell)

!!! --- Initialize stencil coefficients array for Maxwell field solver
CALL FD_weights(xcoeffs, norderx, l_nodalgrid)
CALL FD_weights(ycoeffs, nordery, l_nodalgrid)
CALL FD_weights(zcoeffs, norderz, l_nodalgrid)

!!! --- Sets-up particle space/ velocity distribution
IF ((y_max_local .GT. ymax/2.0_num) .AND. (y_min_local .LT. ymax/2.0_num)) THEN
    ! - 1 electron
    species_parray(1)%species_npart=1
    species_parray(1)%part_x(1)=x_min_local+(x_max_local-x_min_local)/2.0_num
    species_parray(1)%part_y(1)=y_min_local+(y_max_local-y_min_local)/2.0_num
    species_parray(1)%part_z(1)=z_min_local+(z_max_local-z_min_local)/2.0_num
    species_parray(1)%part_uy(1)=-0.8_num*clight
    ! - 1 proton
    species_parray(2)%species_npart=1
    species_parray(2)%part_x(1)=x_min_local+(x_max_local-x_min_local)/2.0_num
    species_parray(2)%part_y(1)=y_min_local+(y_max_local-y_min_local)/2.0_num
    species_parray(2)%part_z(1)=z_min_local+(z_max_local-z_min_local)/2.0_num
    species_parray(2)%part_uy(1)=0.8_num*clight
END IF


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

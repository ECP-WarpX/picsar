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
INTEGER :: i,ierror,j,k,l, ispecies, ipart
INTEGER :: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER :: npartemp
INTEGER :: ixsource, iysource, izsource
REAL(num) :: xsource, ysource, zsource
TYPE(particle_species), POINTER :: curr
!real(8) :: rand

!!! --- Set time step
dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))
it = 0

!!! --- Allocate species arrays
IF (.NOT. l_arrays_allocated) THEN
    DO ispecies=1,nspecies
    curr=>species_parray(ispecies)
    curr%nppspecies_max=2*(curr%nppcell)*nx*ny*nz
    ALLOCATE(curr%part_x(1:curr%nppspecies_max),    &
             curr%part_y(1:curr%nppspecies_max),    &
             curr%part_z(1:curr%nppspecies_max),    &
             curr%part_ux(1:curr%nppspecies_max),   &
             curr%part_uy(1:curr%nppspecies_max),   &
             curr%part_uz(1:curr%nppspecies_max),   &
             curr%part_ex(1:curr%nppspecies_max),   &
             curr%part_ey(1:curr%nppspecies_max),   &
             curr%part_ez(1:curr%nppspecies_max),   &
             curr%part_bx(1:curr%nppspecies_max),   &
             curr%part_by(1:curr%nppspecies_max),   &
             curr%part_bz(1:curr%nppspecies_max),   &
             curr%weight(1:curr%nppspecies_max))
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

! - Init particle species arrays (2 species here)
DO ispecies=1,nspecies
    species_parray(ispecies)%species_npart=0
    species_parray(ispecies)%part_x=0.0_num
    species_parray(ispecies)%part_y=0.0_num
    species_parray(ispecies)%part_z=0.0_num
    species_parray(ispecies)%part_ux=0.0_num
    species_parray(ispecies)%part_uy=0.0_num
    species_parray(ispecies)%part_uz=0.0_num
    species_parray(ispecies)%part_ex=0.0_num
    species_parray(ispecies)%part_ey=0.0_num
    species_parray(ispecies)%part_ez=0.0_num
    species_parray(ispecies)%part_bx=0.0_num
    species_parray(ispecies)%part_by=0.0_num
    species_parray(ispecies)%part_bz=0.0_num
    species_parray(ispecies)%weight=nc*dx*dy*dz/(curr%nppcell) !uniform density for the moment
END DO

!!! --- Sets-up particle space distribution
DO ispecies=1,nspecies
    curr=>species_parray(ispecies)
    IF ((((curr%x_min .GE. x_min_local) .AND. (curr%x_min .LT. x_max_local)) &
        .OR. ((curr%x_max .GE. x_min_local) .AND. (curr%x_max .LT. x_max_local))) .AND. &
        (((curr%x_min .GE. x_min_local) .AND. (curr%x_min .LT. x_max_local)) &
        .OR. ((curr%x_max .GE. x_min_local) .AND. (curr%x_max .LT. x_max_local))) .AND. &
        (((curr%x_min .GE. x_min_local) .AND. (curr%x_min .LT. x_max_local)) &
        .OR. ((curr%x_max .GE. x_min_local) .AND. (curr%x_max .LT. x_max_local))))  THEN
        jmin = NINT(MAX(curr%x_min-x_min_local,0.0_num)/dx)
        jmax = NINT(MIN(curr%x_max-x_min_local,x_max_local-x_min_local)/dx)
        kmin = NINT(MAX(curr%y_min-y_min_local,0.0_num)/dy)
        kmax = NINT(MIN(curr%y_max-y_min_local,y_max_local-y_min_local)/dy)
        lmin = NINT(MAX(curr%z_min-z_min_local,0.0_num)/dz)
        lmax = NINT(MIN(curr%z_max-z_min_local,z_max_local-z_min_local)/dz)
        DO l=lmin+1,lmax-1
            DO k=kmin+1,kmax-1
                DO j=jmin+1,jmax-1
                    DO ipart=1,curr%nppcell
                        curr%species_npart=curr%species_npart+1
                        curr%part_x(curr%species_npart) = x_min_local+(j-1)*dx+dx/curr%nppcell*ipart
                        curr%part_y(curr%species_npart) = y_min_local+(k-1)*dy+dy/curr%nppcell*ipart
                        curr%part_z(curr%species_npart) = z_min_local+(l-1)*dz+dz/curr%nppcell*ipart
                        curr%part_ux(curr%species_npart)= curr%vdrift_x
                        curr%part_uy(curr%species_npart)= curr%vdrift_y
                        curr%part_uz(curr%species_npart)= curr%vdrift_z
                    END DO
                END DO
            END DO
        END DO
    END IF
END DO
!!! --- Initialize stencil coefficients array for Maxwell field solver
CALL FD_weights(xcoeffs, norderx, l_nodalgrid)
CALL FD_weights(ycoeffs, nordery, l_nodalgrid)
CALL FD_weights(zcoeffs, norderz, l_nodalgrid)


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

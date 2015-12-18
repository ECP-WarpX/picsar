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
INTEGER :: nst,i

!!! --- This is the main PIC LOOP
IF (rank .EQ. 0) THEN
WRITE (0,*), "nsteps = ", nst
END IF
DO i=1,nst
    IF (rank .EQ. 0) startit=MPI_WTIME()
    pushtime=0._num
!
!    !!! --- Field gather & particle push
     CALL push_particles
!
!    !!! --- Apply BC on particles
    CALL particle_bcs
!
!    !!! --- Deposit current of particle species on the grid
    CALL depose_currents_on_grid_jxjyjz
!
!    !!! --- Boundary conditions for currents
    CALL current_bcs
!
!    !!! --- Push B field half a time step
    CALL push_bfield
!
!    !!! --- Boundary conditions for B
    CALL bfield_bcs
!
!    !!! --- Push E field  a full time step
    CALL push_efield
!
!    !!! --- Boundary conditions for E
    CALL efield_bcs
!
!    !!! --- push B field half a time step
    CALL push_bfield
!
!    !!! --- Boundary conditions for B
    CALL bfield_bcs
!
!    !!! --- Computes derived quantities
    CALL calc_diags
!
!    !!! --- Output simulation results
!    CALL output_routines

    it = it+1
    timeit=MPI_WTIME()

   IF (rank .EQ. 0) THEN
        WRITE(0,*) 'it = ',it,' || time = ',it*dt, " || push/part (ns)= ", pushtime*1e9_num/ntot, &
        " || tot/part (ns)= ", (timeit-startit)*1e9_num/ntot
    END IF
END DO

END SUBROUTINE step

!===============================================================================
! INIT PLASMA AND FIELD ARRAYS AT it=0
SUBROUTINE initall
!===============================================================================
USE constants
USE params
USE fields
USE particles
USE shared_data
USE tiling

!use IFPORT ! uncomment if using the intel compiler (for rand)
IMPLICIT NONE
INTEGER :: i,ierror,j,k,l, ispecies, ipart, count
INTEGER :: jmin, jmax, lmin, lmax, kmin, kmax
INTEGER :: ix, iy, iz
INTEGER :: npartemp, ncurr
REAL(num) :: v, th, phi
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
!real(8) :: rand

!!! --- Set time step/ it
dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))
it = 0

!!! --- Allocate coefficient arrays for Maxwell solver
IF (.NOT. l_coeffs_allocated) THEN
    ALLOCATE(xcoeffs(norderx/2),ycoeffs(nordery/2),zcoeffs(norderz/2))
END IF

!!! --- Initialize stencil coefficients array for Maxwell field solver
CALL FD_weights(xcoeffs, norderx, l_nodalgrid)
CALL FD_weights(ycoeffs, nordery, l_nodalgrid)
CALL FD_weights(zcoeffs, norderz, l_nodalgrid)


! ------ INIT PARTICLE DISTRIBUTIONS
!!! --- Set tile split for particles
CALL set_tile_split

! - Allocate particle arrays for each tile of each species
CALL init_tile_arrays

! - Load particle distribution on each tile
CALL load_particles

! ----- INIT FIELD ARRAYS
!!! --- Initialize field/currents arrays
! - Init grid arrays
ex=0.0_num;ey=0.0_num;ez=0.0_num
bx=0.0_num;by=0.0_num;bz=0.0_num
jx=0.0_num;jy=0.0_num;jz=0.0_num

!!! --- set number of time steps
nsteps = nint(tmax/(w0_l*dt))

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
INTEGER(KIND=4) :: norder, n, m, mn, i, j, k
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
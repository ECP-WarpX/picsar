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
USE sorting
USE mpi_routines
#if (defined(VTUNE) && VTUNE>0)
USE ITT_FORTRAN                     
#endif                                  
#if (defined(SDE) && SDE>0)||(defined(DFP))
USE SDE_FORTRAN                     
#endif          

IMPLICIT NONE
INTEGER(idp) :: nst,i

!!! --- This is the main PIC LOOP
IF (rank .EQ. 0) THEN
WRITE (0,*) "nsteps = ", nst
END IF

!!! --- Start Vtune/SDE analysis
#if VTUNE==1                      
CALL start_vtune_collection()                
#endif
#if SDE==1                      
CALL start_sde_collection()                
#endif
#if ALLINEA==1      
CALL ALLINEA_START_SAMPLING
#endif                         
! Intel Design Forward project
#if defined(DFP)
 CALL DFP_MAIN_START()
#endif

! ________________________________________________________________________________________
! 
! Main loop
! ________________________________________________________________________________________

! ___________________________________________
! Loop in 3D
IF (c_dim.eq.3) THEN

  DO i=1,nst
      IF (rank .EQ. 0) startit=MPI_WTIME()
    
      !!! --- Init iteration variables
      pushtime=0._num
      divE_computed = .False.
    
      !!! --- Field gather & particle push
      !IF (rank .EQ. 0) PRINT *, "#1"
      CALL field_gathering_plus_particle_pusher
      !IF (rank .EQ. 0) PRINT *, "#2"
      !!! --- Apply BC on particles
      CALL particle_bcs
      !IF (rank .EQ. 0) PRINT *, "#3"
      !!! --- Particle Sorting
      !write(0,*),'Sorting'
      CALL pxr_particle_sorting
      !IF (rank .EQ. 0) PRINT *, "#4"
      !!! --- Deposit current of particle species on the grid
      !write(0,*),'Depose currents'
      CALL pxrdepose_currents_on_grid_jxjyjz  
      !IF (rank .EQ. 0) PRINT *, "#5"
      !!! --- Boundary conditions for currents
      !write(0,*),'Current_bcs'
      CALL current_bcs
      !IF (rank .EQ. 0) PRINT *, "#6"
      !!! --- Push B field half a time step
      !write(0,*),'push_bfield'
      CALL push_bfield
      !IF (rank .EQ. 0) PRINT *, "#7"
      !!! --- Boundary conditions for B
      CALL bfield_bcs
      !IF (rank .EQ. 0) PRINT *, "#8"
      !!! --- Push E field  a full time step
      CALL push_efield
      !IF (rank .EQ. 0) PRINT *, "#9"
      !!! --- Boundary conditions for E
      CALL efield_bcs
      !IF (rank .EQ. 0) PRINT *, "#10"
      !!! --- push B field half a time step
      CALL push_bfield
      !IF (rank .EQ. 0) PRINT *, "#11"
      !!! --- Boundary conditions for B
      CALL bfield_bcs
      !IF (rank .EQ. 0) PRINT *, "#12"
      !!! --- Computes derived quantities
      CALL calc_diags
      !IF (rank .EQ. 0) PRINT *, "#13"
      !!! --- Output simulation results
      CALL output_routines
      !IF (rank .EQ. 0) PRINT *, "#14"
    
      it = it+1
      timeit=MPI_WTIME()

      CALL time_statistics_per_iteration

     IF (rank .EQ. 0)  THEN
          WRITE(0,*) 'it = ',it,' || time = ',it*dt, " || push/part (ns)= ", pushtime*1e9_num/ntot, &
          " || tot/part (ns)= ", (timeit-startit)*1e9_num/ntot
      END IF
  END DO

! ___________________________________________
! Loop in 2D
ELSE IF (c_dim.eq.2) THEN

  DO i=1,nst
      IF (rank .EQ. 0) startit=MPI_WTIME()

      !!! --- Init iteration variables
      pushtime=0._num
      divE_computed = .False.

      !!! --- Field gather & particle push
      CALL field_gathering_plus_particle_pusher

      !!! --- Apply BC on particles
      CALL particle_bcs_2d

      !!! --- Deposit current of particle species on the grid
      CALL pxrdepose_currents_on_grid_jxjyjz_2d

      !!! --- Boundary conditions for currents
      CALL current_bcs

      !!! --- Push B field half a time step
      CALL push_bfield_2d

      it = it+1
      timeit=MPI_WTIME()

     IF (rank .EQ. 0)  THEN
          WRITE(0,*) 'it = ',it,' || time = ',it*dt, " || push/part (ns)= ", pushtime*1e9_num/ntot, &
          " || tot/part (ns)= ", (timeit-startit)*1e9_num/ntot
      END IF
  END DO
  
ENDIF

!!! --- Stop Vtune analysis
#if VTUNE==1            
CALL stop_vtune_collection()      
#endif
#if SDE==1            
CALL stop_sde_collection()      
#endif 
#if ALLINEA==1      
CALL ALLINEA_STOP_SAMPLING
#endif   
! Intel Design Forward project
#if defined(DFP)
 CALL DFP_MAIN_STOP
#endif

! Intel Design Forward project
#if defined(DFP)
 CALL DFP_FINAL_START
#endif

!!! --- Output time statistics
CALL final_output_time_statistics

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
	USE time_stat
	USE precomputed

	!use IFPORT ! uncomment if using the intel compiler (for rand)
	IMPLICIT NONE
	INTEGER(idp)                    :: ispecies,i
	REAL(num)                       :: tdeb
	TYPE(particle_species), POINTER :: curr
	TYPE(particle_dump), POINTER    :: dp 

	! Time statistics
	init_localtimes(:) = 0
	localtimes(:)=0

	! Dimension parameter check
	IF (c_dim.eq.2) THEN
		dy = 1.
		ntiley = 1
		ymin = 0
		ymax = 0
		ny = 1
	ENDIF

	! Few calculations and updates      
	nc    = nlab*g0          ! density (in the simulation frame)
	wlab  = echarge*sqrt(nlab/(emass*eps0)) ! plasma frequency (in the lab frame)
	lambdalab = 2*pi*clight/wlab
	w0_l  = echarge*sqrt(nc/(g0*emass*eps0))    ! "longitudinal" plasma frequency (in the lab frame)
	w0_t  = echarge*sqrt(nc/(g0**3*emass*eps0)) ! "transverse" plasma frequency (in the lab frame)
	w0    = w0_l 

	!!! --- Set time step/ it
	IF (c_dim.eq.3) THEN
		dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))
	ELSE IF (c_dim.eq.2) THEN
		dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dz**2))
	ENDIF
	it = 0

	!!! --- set number of time steps or total time
	if (tmax.eq.0) then
		tmax = nsteps*w0_l*dt
	else
		nsteps = nint(tmax/(w0_l*dt))
	endif

	!!! --- Sorting

	sorting_dx = sorting_dx*dx
	sorting_dy = sorting_dy*dy
	sorting_dz = sorting_dz*dz

	sorting_shiftx = sorting_shiftx*dx
	sorting_shifty = sorting_shifty*dy
	sorting_shiftz = sorting_shiftz*dz

	!!! --- Precomputed parameters
	dxi = 1.0_num/dx
	dyi = 1.0_num/dy
	dzi = 1.0_num/dz
	invvol = dxi*dyi*dzi
	dts2dx = 0.5_num*dt*dxi
	dts2dy = 0.5_num*dt*dyi
	dts2dz = 0.5_num*dt*dzi
	dtsdx0 = dt*dxi
	dtsdy0 = dt*dyi
	dtsdz0 = dt*dzi
	clightsq = 1.0_num/clight**2
	dxs2 = dx*0.5_num
	dys2 = dy*0.5_num
	dzs2 = dz*0.5_num

	!- Init stencil coefficients
	CALL init_stencil_coefficients()

	! Summary
	IF (rank .EQ. 0) THEN
		write(0,*) ''
		write(0,*) 'SIMULATION PARAMETERS:'
		write(0,*) 'Dimension:',c_dim  
		write(0,*) 'dx, dy, dz:',dx,dy,dz
		write(0,*) 'dt:',dt,'s',dt*1e15,'fs'
		write(0,*) 'Total time:',tmax,'plasma periods:',tmax/w0_l,'s'
		write(0,*) 'Number of steps:',nsteps
		write(0,*) 'Tiles:',ntilex,ntiley,ntilez
		write(0,*) 'MPI com current:',mpicom_curr
		write(0,*) 'Current deposition method:',currdepo
		write(0,*) 'Charge deposition algo:',rhodepo
		write(0,*) 'Field gathering method:',fieldgathe
		write(0,*) 'Field gathering plus particle pusher seperated:',fg_p_pp_separated
		write(0,*) 'Current/field gathering order:',nox,noy,noz
		write(0,'(" Particle communication: (partcom=",I1,")")') partcom
		IF (particle_pusher.eq.1) THEN
			write(0,'(" Pusher: Jean-Luc Vay algorithm, particle_pusher=",I1)') particle_pusher
		ELSE
			write(0,'(" Pusher: Boris algorithm (particle_pusher=",I1,")")') particle_pusher
		ENDIF
		write(0,*) 'Maxwell derivative coeff:',xcoeffs
		write(0,*) 'MPI buffer size:',mpi_buf_size 
		WRITE(0,*) ''
		WRITE(0,*) 'Vector length current deposition',lvec_curr_depo
		WRITE(0,*) 'Vector length charge deposition',lvec_charge_depo
		WRITE(0,*) 'Vector length field gathering',lvec_fieldgathe  
		write(0,*) ''
		write(0,*) 'PLASMA PROPERTIES:'
		write(0,*) 'Distribution:',pdistr
		write(0,*) 'Density in the lab frame:',nlab,'m^-3'
		write(0,*) 'Density in the simulation frame:',nc,'m^-3'
		write(0,*) 'Cold plasma frequency in the lab frame:',wlab,'s^-1'
		write(0,*) 'cold plasma wavelength:',lambdalab,'m',lambdalab*1e6,'um'
		write(0,*) ''
	
		write(0,'(" MPI domain decomposition")')  
		write(0,*) 'Topology:',topology
		write(0,'(" Local number of cells:",I5,X,I5,X,I5)') nx,ny,nz
		write(0,'(" Local number of grid point:",I5,X,I5,X,I5)') nx_grid,ny_grid,nz_grid  
		write(0,'(" Guard cells:",I5,X,I5,X,I5)') nxguards,nyguards,nzguards
		write(0,*) ''  
		
		! Sorting
		IF (sorting_activated.gt.0) THEN 
			write(0,*) 'Particle sorting activated'
			write(0,*) 'dx:',sorting_dx
			write(0,*) 'dy:',sorting_dy
			write(0,*) 'dz:',sorting_dz 
			write(0,*) 'shiftx:',sorting_shiftx
			write(0,*) 'shifty:',sorting_shifty
			write(0,*) 'shiftz:',sorting_shiftz    
			write(0,*) ''    
		ELSE
			write(0,*) 'Particle sorting non-activated'
			write(0,*) ''    
		ENDIF
	
		! Species properties
		write(0,*)  'Number of species:',nspecies
		DO ispecies=1,nspecies
			curr => species_parray(ispecies)
			write(0,*) trim(adjustl(curr%name))
			write(0,*) 'Drift velocity:',curr%vdrift_x,curr%vdrift_y,curr%vdrift_z
			write(0,*) 'Sorting period:',curr%sorting_period 
			write(0,*) 'Sorting start:',curr%sorting_start     
			write(0,*) ''
		end do
	
		! Diags
		IF (timestat_activated.gt.0) THEN  
			write(0,*) 'Output of time statistics activated'
			write(0,*) 'Computation of the time statistics starts at',timestat_itstart
			write(0,*) 'Buffer size:',nbuffertimestat
		ELSE
			write(0,*) 'Output of time statistics non-activated'
			write(0,'(X,"Computation of the time statistics starts at iteration:",I5)') timestat_itstart    
		ENDIF
		write(0,*) 
		
		! Particle Dump
		IF (npdumps.gt.0) THEN
			DO i = 1, npdumps
				dp => particle_dumps(i)
				WRITE(0,'(" Dump number: ",I2)') i
				WRITE(0,'(" species name: ",A10)') species_parray(dp%ispecies)%name
				WRITE(0,*) 
			ENDDO
		ELSE
			WRITE(0,'(" No particle dump (",I2,")")') npdumps
			WRITE(0,*) 
		ENDIF
	
	end if

	! ------ INIT PARTICLE DISTRIBUTIONS

	tdeb=MPI_WTIME()

	!!! --- Set tile split for particles
	CALL set_tile_split

	IF (rank .EQ. 0) write(0,*), "Set tile split: done"

	! - Allocate particle arrays for each tile of each species
	CALL init_tile_arrays

	IF (rank .EQ. 0) write(0,*), "Initialization of the tile arrays: done"

	! - Load particle distribution on each tile
	CALL load_particles

	IF (rank .EQ. 0) write(0,*), "Creation of the particles: done"

	init_localtimes(1) = MPI_WTIME() - tdeb

	! - Estimate tile size 
	CALL estimate_memory_consumption

! ----- INIT FIELD ARRAYS
!!! --- Initialize field/currents arrays
! - Init grid arrays
ex=0.0_num;ey=0.0_num;ez=0.0_num
bx=0.0_num;by=0.0_num;bz=0.0_num
jx=0.0_num;jy=0.0_num;jz=0.0_num

END SUBROUTINE initall

!===============================================================================
! INIT STENCIL COEFFICIENTS
!===============================================================================

SUBROUTINE init_stencil_coefficients()
USE constants
USE params
USE fields
USE particles
USE shared_data
USE tiling
IMPLICIT NONE

!!! --- Allocate coefficient arrays for Maxwell solver
IF (.NOT. l_coeffs_allocated) THEN
    ALLOCATE(xcoeffs(norderx/2),ycoeffs(nordery/2),zcoeffs(norderz/2))
    l_coeffs_allocated=.TRUE.
END IF

!!! --- Initialize stencil coefficients array for Maxwell field solver
CALL FD_weights(xcoeffs, norderx, l_nodalgrid)
CALL FD_weights(ycoeffs, nordery, l_nodalgrid)
CALL FD_weights(zcoeffs, norderz, l_nodalgrid)

END SUBROUTINE init_stencil_coefficients


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
INTEGER(idp) :: norder, n, m, mn, i, j, k
LOGICAL      :: l_nodal
REAL(num)    :: z, fact, c1, c2, c3, c4, c5
REAL(num), INTENT(IN OUT), DIMENSION(norder/2) :: coeffs
REAL(num), ALLOCATABLE, DIMENSION(:)           :: x
REAL(num), ALLOCATABLE, DIMENSION(:,:)         :: c

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

! ____________________________
! For debugging


!!! --- Test current
SUBROUTINE current_debug
  USE fields
  USE shared_data
  IMPLICIT NONE
  
  INTEGER :: i
  
  !jx(1:nx,1:ny,1:nz) = 1.
  !jx(1,1:ny,1:nz) = 0.5
  !jx(nx,ny,nz) = 0.5
  i = 0
  jy(i:nx,i:ny,i:nz) = 1.

  jy(i,i:ny,i:nz) = 0.5
  jy(nx,i:ny,i:nz) = 0.5
  jy(i:nx,i,i:nz) = 0.5
  jy(i:nx,ny,i:nz) = 0.5
  jy(i:nx,i:ny,i) = 0.5
  jy(i:nx,i:ny,nz) = 0.5

  jy(i:nx,i,i) = 0.25
  jy(i,i:ny,i) = 0.25
  jy(i,i,i:nz) = 0.25
  jy(i:nx,ny,nz) = 0.25
  jy(nx,i:ny,nz) = 0.25
  jy(nx,ny,i:nz) = 0.25

  jy(i,i,i) = 0.125
  jy(nx,i,i) = 0.125
  jy(i,ny,i) = 0.125
  jy(i,i,nz) = 0.125
  jy(nx,ny,i) = 0.125
  jy(nx,i,nz) = 0.125
  jy(i,ny,nz) = 0.125
  jy(nx,ny,nz) = 0.125
  !jz(1:nx,1:ny,1:nz) = 1.
  !jz(1,1,1) = 0.5
  !jz(nx,ny,nz) = 0.5
  !!! --- End debug
END SUBROUTINE


! ________________________________

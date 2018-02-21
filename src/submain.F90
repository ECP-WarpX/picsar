! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software, ! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! SUBMAIN.F90
!
! This file is for the main PIC loop without initialization.
!
! Developers:
! Henri Vincenti
! Mathieu Lobet
!
! Date:
! Creation 2015
!
! Modifications:
! Mathieu Lobet - 2016 - Summary of the main parameters printed
!                        at the beginning of the simulation.
! Mathieu Lobet - 2016 - Creation of a partial 2D loop (does not include all 2d steps)
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Subroutine that performs the time loop of the PIC algorithm.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
!
!> @param[in] nst number of time steps
!
! ________________________________________________________________________________________
SUBROUTINE step(nst)

  USE constants
  USE fields
  USE particles
  USE params
  USE shared_data
  USE field_boundary
  USE particle_boundary
  USE diagnostics
  USE simple_io
  USE sorting
  USE mpi_routines
#if defined(FFTW)
  USE gpstd_solver
  USE mpi_fftw3
  USE group_parameters
#endif
#if (defined(VTUNE) && VTUNE>0)
  USE ITT_FORTRAN
#endif
#if (defined(SDE) && SDE>0)||(defined(DFP))
  USE SDE_FORTRAN
#endif

  IMPLICIT NONE
  INTEGER(idp) :: nst, i

  !!! --- This is the main PIC LOOP
  IF (rank .EQ. 0) THEN
    WRITE (0, *) "nsteps = ", nst
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

  ! ______________________________________________________________________________________
  !
  ! Main loop
  ! ______________________________________________________________________________________

  ! ___________________________________________
  ! Loop in 3D
  IF (c_dim.eq.3) THEN
    rhoold=0.0_num
    rho = 0.0_num
    DO i=1, nst
      IF (rank .EQ. 0) startit=MPI_WTIME()

      !!! --- Init iteration variables
      pushtime=0._num
      divE_computed = .False.
      IF (l_plasma) THEN
        !!! --- Field gather & particle push
        !IF (rank .EQ. 0) PRINT *, "#1"
        CALL field_gathering_plus_particle_pusher
        !IF (rank .EQ. 0) PRINT *, "#2"
        !!! --- Push virtual laser particles
        CALL push_laser_particles
        !!! --- Apply BC on particles
        CALL particle_bcs
        !IF (rank .EQ. 0) PRINT *, "#3"
#if defined(FFTW)
        IF (l_spectral) THEN
          CALL  copy_field(rhoold, nx+2*nxguards+1, ny+2*nyguards+1,      &
                nz+2*nzguards+1, rho, nx+2*nxguards+1, ny+2*nyguards+1,   &
                nz+2*nzguards+1)
          CALL pxrdepose_rho_on_grid
          CALL charge_bcs
        ENDIF
#endif
        !!! --- Particle Sorting
        !write(0, *), 'Sorting'
        CALL pxr_particle_sorting
        !IF (rank .EQ. 0) PRINT *, "#4"
        !!! --- Deposit current of particle species on the grid
        !write(0, *), 'Depose currents'
        CALL pxrdepose_currents_on_grid_jxjyjz
        !IF (rank .EQ. 0) PRINT *, "#5"
        !!! --- Boundary conditions for currents
        !write(0, *), 'Current_bcs'
        CALL current_bcs
      ENDIF
#if defined(FFTW)
      IF (l_spectral) THEN

        !!! --- FFTW FORWARD - FIELD PUSH - FFTW BACKWARD
        CALL push_psatd_ebfield_3d
        !IF (rank .EQ. 0) PRINT *, "#0"
        !!! --- Boundary conditions for E AND B
        CALL efield_bcs
        CALL bfield_bcs
        IF(absorbing_bcs) THEN
          CALL field_damping_bcs()
          CALL merge_fields()
        ENDIF
      ELSE
#endif
        !IF (rank .EQ. 0) PRINT *, "#6"
        !!! --- Push B field half a time step
        !write(0, *), 'push_bfield'
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
#if defined(FFTW)
      ENDIF
#endif
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
        WRITE(0, *) 'it = ', it, ' || time = ', it*dt, " || push/part (ns)= ",        &
        pushtime*1e9_num/ntot, " || tot/part (ns)= ", (timeit-startit)*1e9_num/ntot
      END IF
    END DO

    ! ___________________________________________
    ! Loop in 2D
  ELSE IF (c_dim.eq.2) THEN

    DO i=1, nst
      IF (rank .EQ. 0) startit=MPI_WTIME()

      !!! --- Init iteration variables
      pushtime=0._num
      divE_computed = .False.

      !!! --- Field gather & particle push
      CALL field_gathering_plus_particle_pusher
      call push_laser_particles()

      !!! --- Apply BC on particles
      CALL particle_bcs_2d

      !!! --- Deposit current of particle species on the grid
      CALL pxr_particle_sorting

      CALL pxrdepose_currents_on_grid_jxjyjz_2d

      !!! --- Boundary conditions for currents
      CALL current_bcs
#if defined(FFTW)
        IF (l_spectral) THEN
          CALL  copy_field(rhoold, nx+2*nxguards+1, ny+2*nyguards+1,      &
                nz+2*nzguards+1, rho, nx+2*nxguards+1, ny+2*nyguards+1,   &
                nz+2*nzguards+1)
          CALL pxrdepose_rho_on_grid
          CALL charge_bcs
        ENDIF
#endif

#if defined(FFTW)
      IF (l_spectral) THEN
        !!! --- FFTW FORWARD - FIELD PUSH - FFTW BACKWARD
        CALL push_psatd_ebfield_2d
        !IF (rank .EQ. 0) PRINT *, "#0"
        !!! --- Boundary conditions for E AND B
        CALL efield_bcs
        CALL bfield_bcs
        IF (absorbing_bcs) THEN
          CALL field_damping_bcs()
          CALL merge_fields()
        ENDIF
      ELSE
#endif
        !IF (rank .EQ. 0) PRINT *, "#6"
        !!! --- Push B field half a time step
        !write(0, *), 'push_bfield'
        CALL push_bfield_2d
        !IF (rank .EQ. 0) PRINT *, "#7"
        !!! --- Boundary conditions for B
        CALL bfield_bcs
        !IF (rank .EQ. 0) PRINT *, "#8"
        !!! --- Push E field  a full time step
        CALL push_efield_2d
        !IF (rank .EQ. 0) PRINT *, "#9"
        !!! --- Boundary conditions for E
        CALL efield_bcs
        !IF (rank .EQ. 0) PRINT *, "#10"
        !!! --- push B field half a time step
        CALL push_bfield_2d
        !IF (rank .EQ. 0) PRINT *, "#11"
        !!! --- Boundary conditions for B
        CALL bfield_bcs

#if defined(FFTW)
      ENDIF
#endif
      !IF (rank .EQ. 0) PRINT *, "#12"
      !!! --- Computes derived quantities
      CALL calc_diags
      !IF (rank .EQ. 0) PRINT *, "#13"
      !!! --- Output simulation results
      CALL output_routines
      !IF (rank .EQ. 0) PRINT *, "#14"
      it = it +1
      timeit=MPI_WTIME()

      CALL time_statistics_per_iteration

      IF (rank .EQ. 0)  THEN
        WRITE(0, *) 'it = ', it, ' || time = ', it*dt, " || push/part (ns)= ",        &
        pushtime*1e9_num/ntot, " || tot/part (ns)= ", (timeit-startit)*1e9_num/ntot
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
! ________________________________________________________________________________________
!> @brief
!> init pml arrays 
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2018

SUBROUTINE init_pml_arrays
  USE field_boundary
  USE fields
  USE shared_data
  USE constants  
  USE params, ONLY : dt 
  USE mpi
  USE mpi_derived_types

  LOGICAL(lp)  :: is_intersection_x, is_intersection_y, is_intersection_z
  INTEGER(idp) :: ix,iy,iz,pow
  REAL(num)    :: coeff,b_offset, e_offset
  INTEGER(idp) :: type_id  
  REAL(num)    , ALLOCATABLE, DIMENSION(:) :: temp
  coeff = 1.80_num
  b_offset = 0.5_num
  e_offset = 0._num
  pow = 4_idp

  !> Inits pml arrays of the same size as ex fields in the daming direction!
  !> sigmas are 1d arrray to economize memory
  
  !> Pml arrays are initializd to 0.0_num because exp(0.0_num) = 1.0!
  !> So in case of no pml region inside the mpi domain, pml acts as vaccum

  !> first each proc allocates sigma as an array of size n_global + 2*ng +1
  !> then each proc will compute sigma in the whole domain and will then extract
  !> the relevent part to its domain decomposition
  !> Note that a pml region can overlap many procs even with local psatd !
  !> In this case (when pml overlaps many procs) field guardcells can act as a
  !> pml region too (imagine nx_pml = 15 and  nx_local = 10, then 5 cells of
  !> proc number 2 will be pml regions. Say that we have 3 nxguards, then the 3
  !> guardcells and 2 physical cells are pmls

  ALLOCATE(sigma_x_e(-nxguards:nx_global+nxguards));
  sigma_x_e = 0.0_num
  ALLOCATE(sigma_x_b(-nxguards:nx_global+nxguards));
  sigma_x_b = 0.0_num
  ALLOCATE(sigma_y_e(-nyguards:ny_global+nyguards));
  sigma_y_e = 0.0_num
  ALLOCATE(sigma_y_b(-nyguards:ny_global+nyguards));
  sigma_y_b = 0.0_num
  ALLOCATE(sigma_z_e(-nzguards:nz_global+nzguards));
  sigma_z_e = 0.0_num
  ALLOCATE(sigma_z_b(-nzguards:nz_global+nzguards));
  sigma_z_b = 0.0_num

  !> Inits sigma_x_e and sigma_x_b in the lower bound of the domain along x
  !> axis
  !> first, each proc will compute sigma in the whole domain
  !> Need more straightforward way to do this
  DO ix = 0,nx_pml-1
    sigma_x_e(ix) =  coeff*clight/dx*(nx_pml-ix-e_offset)**pow
    sigma_x_b(ix) =  coeff*clight/dx*(nx_pml-ix-b_offset)**pow
  ENDDO

  ! > Inits sigma_x_e and sigma_x_b in the upper bound of the domain along x
  !> axis
  !> first, each proc will compute sigma in the whole domain
  !> Need more straightforward way to do this
  DO ix = nx_global-nx_pml, nx_global-1
    sigma_x_e(ix) = coeff*clight/dx *(ix-(nx_global-nx_pml-1)+e_offset)**pow
    sigma_x_b(ix) = coeff*clight/dx *(ix-(nx_global-nx_pml-1)+b_offset)**pow
  ENDDO

  !> Each proc extracts the relevent part of sigma 
  ALLOCATE(temp(-nxguards:nx_global+nxguards)) 
  temp = sigma_x_e
  DEALLOCATE(sigma_x_e); ALLOCATE(sigma_x_e(-nxguards:nx+nxguards))
  sigma_x_e = temp(cell_x_min(x_coords+1)-nxguards:cell_x_max(x_coords+1)+1+nxguards)
 
  temp = sigma_x_b
  DEALLOCATE(sigma_x_b); ALLOCATE(sigma_x_b(-nxguards:nx+nxguards))
  sigma_x_b = temp(cell_x_min(x_coords+1)-nxguards:cell_x_max(x_coords+1)+1+nxguards)
  DEALLOCATE(temp)
 
  ! > Inits sigma_y_e and sigma_y_b in the lower bound of the domain along y
  !> axis
  !> first, each proc will compute sigma in the whole domain
  !> Need more straightforward way to do this
  DO iy = 0 , ny_pml-1
    sigma_y_e(iy) =  coeff*clight/dy*(ny_pml-iy-e_offset)**pow
    sigma_y_b(iy) =  coeff*clight/dy*(ny_pml-iy-b_offset)**pow
  ENDDO

  ! > Inits sigma_y_e and sigma_y_b in the upper bound of the domain along y
  !> axis
  !> first, each proc will compute sigma in the whole domain
  !> Need more straightforward way to do this
  DO iy = ny_global-ny_pml,ny_global-1
    sigma_y_e(iy) = coeff*clight/dy*(iy-(ny_global-ny_pml-1)+e_offset)**pow
    sigma_y_b(iy) = coeff*clight/dy*(iy-(ny_global-ny_pml-1)+b_offset)**pow
  ENDDO
  
  !> Each proc extracts the relevent part of sigma 
  ALLOCATE(temp(-nyguards:ny_global+nyguards)) 
  temp = sigma_y_e
  DEALLOCATE(sigma_y_e); ALLOCATE(sigma_y_e(-nyguards:ny+nyguards))
  sigma_y_e = temp(cell_y_min(y_coords+1)-nyguards:cell_y_max(y_coords+1)+1+nyguards)
  temp = sigma_y_b
  DEALLOCATE(sigma_y_b); ALLOCATE(sigma_y_b(-nyguards:ny+nyguards))
  sigma_y_b = temp(cell_y_min(y_coords+1)-nyguards:cell_y_max(y_coords+1)+1+nyguards)
  DEALLOCATE(temp)

  ! > Inits sigma_z_e and sigma_z_b in the lower bound of the domain along z
  !> axis
  !> first, each proc will compute sigma in the whole domain
  !> Need more straightforward way to do this
  DO iz =0 , nz_pml-1
    sigma_z_e(iz) =  coeff*clight/dz*(nz_pml-iz-e_offset)**pow
    sigma_z_b(iz) =  coeff*clight/dz*(nz_pml-iz-b_offset)**pow
  ENDDO

  ! > Inits sigma_z_e and sigma_z_b in the upper bound of the domain along z
  !> axis
  !> first, each proc will compute sigma in the whole domain
  !> Need more straightforward way to do this
  DO iz =  nz_global-nz_pml,nz_global-1 
     sigma_z_e(iz) = coeff*clight/dz*(iz-(nz_global-nz_pml-1)+e_offset)**pow
     sigma_z_b(iz) = coeff*clight/dz*(iz-(nz_global-nz_pml-1)+b_offset)**pow
  ENDDO

  !> Each proc extracts the relevent part of sigma 
  ALLOCATE(temp(-nzguards:nz_global+nzguards))
  temp = sigma_z_e
  DEALLOCATE(sigma_z_e); ALLOCATE(sigma_z_e(-nzguards:nz+nzguards))
  sigma_z_e = temp(cell_z_min(z_coords+1)-nzguards:cell_z_max(z_coords+1)+1+nzguards)

  temp = sigma_z_b 
  DEALLOCATE(sigma_z_b); ALLOCATE(sigma_z_b(-nzguards:nz+nzguards))
  sigma_z_b = temp(cell_z_min(z_coords+1)-nzguards:cell_z_max(z_coords+1)+1+nzguards)
  DEALLOCATE(temp)
 
  !> sigma=exp(-sigma*dt) 
  !> Uses an exact formulation to damp fields :
  !> dE/dt = -sigma * E => E(n)=exp(-sigma*dt)*E(n-1)
  !> Note that fdtd pml solving requires field time centering
  sigma_x_e = EXP(-sigma_x_e*dt)
  sigma_y_e = EXP(-sigma_y_e*dt)
  sigma_z_e = EXP(-sigma_z_e*dt)
  sigma_x_b = EXP(-sigma_x_b*dt)
  sigma_y_b = EXP(-sigma_y_b*dt)
  sigma_z_b = EXP(-sigma_z_b*dt)


END SUBROUTINE init_pml_arrays

! ________________________________________________________________________________________
!> @brief
!> Initialize the plasma and field arrays at it=0.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
SUBROUTINE initall
  ! ________________________________________________________________________________________
  USE constants
  USE params
  USE fields
  USE particles
  USE shared_data
  USE tiling
  USE time_stat
#if defined(FFTW)
  USE group_parameters
  USE fourier_psaotd
  USE gpstd_solver
#endif
  USE precomputed

  !use IFPORT ! uncomment if using the intel compiler (for rand)
  IMPLICIT NONE
  INTEGER(idp)                    :: ispecies, i
  REAL(num)                       :: tdeb
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_dump), POINTER    :: dp

  ! Time statistics
  init_localtimes(:) = 0
  localtimes(:)=0

  ! Dimension parameter check
  IF (c_dim.eq.2) THEN
    dy = HUGE(1.0_num)
    ntiley = 1_idp
    ymin = 0_idp
    ymax = 0_idp
    ny = 1_idp
    nyguards = 0_idp
    ymin = -HUGE(1.0_num) 
    ymax = HUGE(1.0_num) 
    y_min_local = -HUGE(1.0_num)
    y_max_local = HUGE(1.0_num)
    cell_y_min = 0_idp
    cell_y_max = 0_idp
  ENDIF

  ! Few calculations and updates
  nc    = nlab*g0! density (in the simulation frame)
  wlab  = echarge*sqrt(nlab/(emass*eps0))! plasma frequency (in the lab frame)
  lambdalab = 2*pi*clight/wlab
  w0_l  = echarge*sqrt(nc/(g0*emass*eps0))! "longitudinal" plasma frequency (in the lab frame)
  w0_t  = echarge*sqrt(nc/(g0**3*emass*eps0))! "transverse" plasma frequency (in the lab frame)
  w0    = w0_l

  !!! --- Set time step/ it
  IF (c_dim.eq.3) THEN
    IF (l_spectral) THEN
      dt=MIN(dx, dy, dz)/clight
    ELSE
      dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))
    ENDIF
  ELSE IF (c_dim.eq.2) THEN
    IF (l_spectral) THEN 
      dt=MIN(dx,dz)/clight
    ELSE
      dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dz**2))
    ENDIF
  ENDIF
  it = 0

  !!! --- set number of time steps or total time
  if (nsteps .eq. 0) then
    nsteps = nint(tmax/(w0_l*dt))
  else
    tmax = nsteps*w0_l*dt
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

  ! - Init stencil coefficients for FDTD

  IF(.NOT. l_spectral) THEN
    CALL init_stencil_coefficients()
  ENDIF
  ! Summary
  IF (rank .EQ. 0) THEN
    write(0, *) ''
    write(0, *) 'SIMULATION PARAMETERS:'
    write(0, *) 'Dimension:', c_dim
    write(0, *) 'dx, dy, dz:', dx, dy, dz
    write(0, *) 'dt:', dt, 's', dt*1e15, 'fs'
    write(0, '(" Coefficient on dt determined via the CFL (dtcoef): ", F12.5)')       &
    dtcoef
    write(0, *) 'Total time:', tmax, 'plasma periods:', tmax/w0_l, 's'
    write(0, *) 'Number of steps:', nsteps
    write(0, *) 'Tiles:', ntilex, ntiley, ntilez
    write(0, *) 'MPI com current:', mpicom_curr
    write(0, *) 'Current deposition method:', currdepo
    write(0, *) 'Charge deposition algo:', rhodepo
    write(0, *) 'Field gathering method:', fieldgathe
    write(0, *) 'Field gathering plus particle pusher seperated:', fg_p_pp_separated
    write(0, *) 'Current/field gathering order:', nox, noy, noz
    write(0, '(" Particle communication: (partcom=", I1, ")")') partcom
    IF (particle_pusher.eq.1) THEN
      write(0, '(" Pusher: Jean-Luc Vay algorithm, particle_pusher=", I1)')           &
      particle_pusher
    ELSE
      write(0, '(" Pusher: Boris algorithm (particle_pusher=", I1, ")")')             &
      particle_pusher
    ENDIF
    IF(.NOT. l_spectral) THEN
      write(0, *) 'Maxwell derivative coeff:', xcoeffs
    ENDIF
    write(0, *) 'MPI buffer size:', mpi_buf_size
    WRITE(0, *) ''
    WRITE(0, *) 'Vector length current deposition', lvec_curr_depo
    WRITE(0, *) 'Vector length charge deposition', lvec_charge_depo
    WRITE(0, *) 'Vector length field gathering', lvec_fieldgathe
    write(0, *) ''
    write(0, *) 'PLASMA PROPERTIES:'
    write(0, *) 'Distribution:', pdistr
    write(0, *) 'Density in the lab frame:', nlab, 'm^-3'
    write(0, *) 'Density in the simulation frame:', nc, 'm^-3'
    write(0, *) 'Cold plasma frequency in the lab frame:', wlab, 's^-1'
    write(0, *) 'cold plasma wavelength:', lambdalab, 'm', lambdalab*1e6, 'um'
    write(0, *) ''

    write(0, '(" MPI domain decomposition")')
    write(0, *) 'Topology:', topology
    write(0, '(" Local number of cells:", I5, X, I5, X, I5)') nx, ny, nz
    write(0, '(" Local number of grid point:", I5, X, I5, X, I5)') nx_grid, ny_grid,  &
    nz_grid
    write(0, '(" Guard cells:", I5, X, I5, X, I5)') nxguards, nyguards, nzguards
    write(0, *) ''
    IF(absorbing_bcs) THEN
      write(0, '(" Absorbing field bcs with n_pml =:", I5, X, I5, X, I5)') nx_pml, ny_pml,nz_pml
    ELSE 
      write(0, '(" Periodic field bcs")')
    ENDIF
    write(0, '(" FFTW - parameters ")')
#if defined(FFTW)
    IF (l_spectral)    write(0, '(" PSATD Maxwell Solver")')
    IF (fftw_with_mpi) write(0, '(" FFTW distributed version - MPI ")')
    IF (fftw_hybrid)   write(0, '(" FFTW distributed version")')
    IF (p3dfft_flag)   write(0, '(" USING PDFFT")')
    IF(p3dfft_stride .AND. p3dfft_flag) write(0, '(" USING STRIDED  PDFFT")')
    IF(p3dfft_stride .EQV. .FALSE. .AND. p3dfft_flag) write(0, '(" USING UNSTRIDED PDFFT")')
    IF(fftw_hybrid) write(0, '(" nb_groups :", I5, X, I5, X, I5)') nb_group_x,nb_group_y,nb_group_z
    IF (fftw_threads_ok) write(0, '(" FFTW MPI - Threaded support enabled ")')
    IF (fftw_mpi_transpose) write(0, '(" FFTW MPI Transpose plans enabled ")')
#endif
    ! Sorting
    IF (sorting_activated.gt.0) THEN
      write(0, *) 'Particle sorting activated'
      write(0, *) 'dx:', sorting_dx
      write(0, *) 'dy:', sorting_dy
      write(0, *) 'dz:', sorting_dz
      write(0, *) 'shiftx:', sorting_shiftx
      write(0, *) 'shifty:', sorting_shifty
      write(0, *) 'shiftz:', sorting_shiftz
      write(0, *) ''
    ELSE
      write(0, *) 'Particle sorting non-activated'
      write(0, *) ''
    ENDIF

    ! Species properties
    write(0, *)  'Number of species:', nspecies
    DO ispecies=1, nspecies
      curr => species_parray(ispecies)
      write(0, *) trim(adjustl(curr%name))
      write(0, *) 'Charge:', curr%charge
      write(0, *) 'Drift velocity:', curr%vdrift_x, curr%vdrift_y, curr%vdrift_z
      write(0, *) 'Sorting period:', curr%sorting_period
      write(0, *) 'Sorting start:', curr%sorting_start
      write(0, *) ''
    end do

    ! Diags
    IF (timestat_activated.gt.0) THEN
      write(0, *) 'Output of time statistics activated'
      write(0, *) 'Computation of the time statistics starts at', timestat_itstart
      write(0, *) 'Buffer size:', nbuffertimestat
    ELSE
      write(0, *) 'Output of time statistics non-activated'
      write(0, '(X, "Computation of the time statistics starts at iteration:", I5)')  &
      timestat_itstart
    ENDIF
    write(0, *)

    ! Particle Dump
    IF (npdumps.gt.0) THEN
      DO i = 1, npdumps
        dp => particle_dumps(i)
        WRITE(0, '(" Dump number: ", I2)') i
        WRITE(0, '(" species name: ", A10)') species_parray(dp%ispecies)%name
        WRITE(0, *)
      ENDDO
    ELSE
      WRITE(0, '(" No particle dump (", I2, ")")') npdumps
      WRITE(0, *)
    ENDIF

  end if

  ! ------ INIT PARTICLE DISTRIBUTIONS

  tdeb=MPI_WTIME()

  !!! --- Set tile split for particles
  CALL set_tile_split

  IF (rank .EQ. 0) write(0, *) "Set tile split: done"

  ! - Allocate particle arrays for each tile of each species
  CALL init_tile_arrays

  IF (rank .EQ. 0) write(0, *) "Initialization of the tile arrays: done"

  ! - Load particle distribution on each tile
  CALL load_particles

  ! - Load laser antenna particles
  CALL load_laser

  IF (rank .EQ. 0) write(0, *) "Creation of the particles: done"

  init_localtimes(1) = MPI_WTIME() - tdeb

  ! - If absorbing bcs then init pml arrays
  IF(absorbing_bcs) THEN
    CALL init_pml_arrays
  ENDIF

#if defined(FFTW)
  ! -Init Fourier
  IF (l_spectral) THEN
    CALL init_plans_blocks
  ENDIF
#endif
  ! - Estimate tile size
  CALL estimate_total_memory_consumption

  ! ----- INIT FIELD ARRAYS
  !!! --- Initialize field/currents arrays
  ! - Init grid arrays
  ex=0.0_num;ey=0.0_num;ez=0.0_num
  bx=0.0_num;by=0.0_num;bz=0.0_num
  jx=0.0_num;jy=0.0_num;jz=0.0_num
  IF(absorbing_bcs) THEN
    exy = 0.5_num * ex
    exz = 0.5_num * ex
    eyx = 0.5_num * ey
    eyz = 0.5_num * ey
    ezx = 0.5_num * ez
    ezy = 0.5_num * ez
    bxy = 0.5_num * bx
    bxz = 0.5_num * bx
    byx = 0.5_num * by
    byz = 0.5_num * by
    bzx = 0.5_num * bz
    bzy = 0.5_num * bz
  ENDIF
END SUBROUTINE initall

! ________________________________________________________________________________________
!> @brief
!> Subroutine that computes total memory consumption (particle/grid tile structures and 
!> regular grid arrays used in the Maxwell solver)
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________
SUBROUTINE estimate_total_memory_consumption
  USE tiling
  USE mpi_routines
  USE mem_status
  IMPLICIT NONE 
  REAL(num) :: total_memory=0._num, avg_per_mpi=0._num
  ! - Get local/global memory occupied by tile arrays (grid and particles)
  CALL get_local_tile_mem()
  CALL get_global_tile_mem()

  ! - Get local/global memory occupied by grid arrays (Maxwell solver)
  CALL get_local_grid_mem()
  CALL get_global_grid_mem()

  ! - Output results on standard output 
  IF (rank .EQ. 0) THEN 
    total_memory=global_grid_mem+global_grid_tiles_mem+global_part_tiles_mem
    avg_per_mpi=total_memory/nproc
    WRITE(0, *) 'Total memory (GB) for grid arrays ', global_grid_mem/1e9
    WRITE(0, *) 'Total memory (GB) for grid tile arrays ', global_grid_tiles_mem/1e9
    WRITE(0, *) 'Total memory (GB) for particle tile arrays ', global_part_tiles_mem/1e9
    WRITE(0, *) 'Total memory (GB)', total_memory/1e9
    WRITE(0, *) 'Avg memory (GB) /MPI process ', avg_per_mpi/1e9
  ENDIF 
END SUBROUTINE estimate_total_memory_consumption



! ________________________________________________________________________________________
!> @brief
!> Initialize stencil coefficients.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
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
    ALLOCATE(xcoeffs(norderx/2), ycoeffs(nordery/2), zcoeffs(norderz/2))
    l_coeffs_allocated=.TRUE.
  END IF

  !!! --- Initialize stencil coefficients array for Maxwell field solver
  CALL FD_weights(xcoeffs, norderx, l_nodalgrid)
  CALL FD_weights(ycoeffs, nordery, l_nodalgrid)
  CALL FD_weights(zcoeffs, norderz, l_nodalgrid)

END SUBROUTINE init_stencil_coefficients


! ________________________________________________________________________________________
!> @brief
!> Compute stencil coefficients for Maxwell field solver
!
!> Adapted from Matlab code from Fornberg (1998)
!> Calculates FD weights. The parameters are:
!> @param[in] z location where approximations are to be accurate.
!> n   number of grid points, !> m   highest derivative that we want to find weights for
!> c   array size m+1, length(x) containing (as output) in
!> successive rows the weights for derivatives 0, 1, ..., m.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
SUBROUTINE FD_weights(coeffs, norder, l_nodal)
  ! ________________________________________________________________________________________

  USE constants
  IMPLICIT NONE
  INTEGER(idp) :: norder, n, m, mn, i, j, k
  LOGICAL(lp)  :: l_nodal
  REAL(num)    :: z, fact, c1, c2, c3, c4, c5
  REAL(num), INTENT(IN OUT), DIMENSION(norder/2) :: coeffs
  REAL(num), ALLOCATABLE, DIMENSION(:)           :: x
  REAL(num), ALLOCATABLE, DIMENSION(:, :)         :: c

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
  ALLOCATE(c(0:m, 0:n-1))

  DO i=0, n-1
    x(i)=(i-n/2+1)*1.0_num
  END DO

  c=0.0_num; c1=1.0_num; c4=x(0)-z; c(0, 0)=1.0_num
  DO i=1, n-1
    mn=min(i+1, m+1)
    c2=1.0_num
    c5=c4
    c4=x(i)-z
    DO j=0, i-1
      c3=x(i)-x(j)
      c2=c2*c3
      IF (j .EQ. (i-1)) THEN
        DO k=1, mn-1
          c(k, i)=c1*(k*c(k-1, i-1)-c5*c(k, i-1))/c2
        END DO
        c(0, i)=-c1*c5*c(0, i-1)/c2
        DO k=1, mn-1
          c(k, j)=(c4*c(k, j)-k*c(k-1, j))/c3
        END DO
        c(0, j)=c4*c(0, j)/c3
      END IF

    END DO
    c1=c2
  END DO

  DO i=1, norder/2
    coeffs(i)=c(m, norder/2+i-1)
  END DO
  RETURN
END SUBROUTINE FD_weights

! ______________________________________
! For debugging


! ________________________________________________________________________________________
!> @brief
!> Subroutine to test current
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
SUBROUTINE current_debug
  ! ________________________________________________________________________________________
  USE fields
  USE shared_data
  IMPLICIT NONE

  INTEGER :: i

  !jx(1:nx, 1:ny, 1:nz) = 1.
  !jx(1, 1:ny, 1:nz) = 0.5
  !jx(nx, ny, nz) = 0.5
  i = 0
  jy(i:nx, i:ny, i:nz) = 1.

  jy(i, i:ny, i:nz) = 0.5
  jy(nx, i:ny, i:nz) = 0.5
  jy(i:nx, i, i:nz) = 0.5
  jy(i:nx, ny, i:nz) = 0.5
  jy(i:nx, i:ny, i) = 0.5
  jy(i:nx, i:ny, nz) = 0.5

  jy(i:nx, i, i) = 0.25
  jy(i, i:ny, i) = 0.25
  jy(i, i, i:nz) = 0.25
  jy(i:nx, ny, nz) = 0.25
  jy(nx, i:ny, nz) = 0.25
  jy(nx, ny, i:nz) = 0.25

  jy(i, i, i) = 0.125
  jy(nx, i, i) = 0.125
  jy(i, ny, i) = 0.125
  jy(i, i, nz) = 0.125
  jy(nx, ny, i) = 0.125
  jy(nx, i, nz) = 0.125
  jy(i, ny, nz) = 0.125
  jy(nx, ny, nz) = 0.125
  !jz(1:nx, 1:ny, 1:nz) = 1.
  !jz(1, 1, 1) = 0.5
  !jz(nx, ny, nz) = 0.5
  !!! --- End debug
END SUBROUTINE current_debug

! ______________________________________

! ______________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! TILE_FIELD_GATHERING_3D_TEST.F90
!
! Test code for the field gathering with tiling in 3D
!
! Mathieu Lobet, 2016.08
! ______________________________________________________________________________

PROGRAM tile_field_gathering_3d_test
  USE constants
  USE fields
  USE particles
  USE params
  USE shared_data
  USE mpi_routines
  USE control_file
  USE tiling

  ! ______________________________________________________________________________________
  ! Parameters

  TYPE(particle_species), POINTER :: curr
  INTEGER(idp)                    :: jmin, jmax, kmin, kmax, lmin, lmax
  REAL(num)                       :: partx, party, partz, partux, partuy, partuz, gaminv
  REAL(num), DIMENSION(:), ALLOCATABLE :: partpid
  REAL(num)                       :: phi, th, up, clightsq
  REAL(num)                                :: Ef, Bf
  REAL(num), DIMENSION(6)                  :: rng=0_num
  REAL(num)                                :: epsilon
  REAL(num)                                :: t0
  LOGICAL(lp)                              :: passed
  REAL(num), dimension(10)                 :: t
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumex,tilesumey,tilesumez
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumbx,tilesumby,tilesumbz
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  CHARACTER(len=64)                        :: title

  ! ______________________________________________________________________________________
  ! Initialization
  ! --- default init
  CALL default_init

  ! --- Dimension
  c_dim = 3

  ! --- Number of processors
  nprocx=1
  nprocy=1
  nprocz=1

  ! --- Domain size
  nx_global_grid=50
  ny_global_grid=50
  nz_global_grid=50

  ! --- Domain extension
  xmin=0
  ymin=0
  zmin=0
  xmax=1e-6
  ymax=1e-6
  zmax=1e-6

  ! --- Init particle tiling split
  ntilex = 6
  ntiley = 6
  ntilez = 6

  ! --- Vector length field gathering
  LVEC_fieldgathe = 256

  ! --- Interpolation
  l_lower_order_in_v = .TRUE.

  ! --- Field properties
  Ef = 1E12 !V/m
  Bf = 1E4  !T

  ! --- Relative error to pass the test
  epsilon = 1e-6_num
  passed=.TRUE.

  ! --- Init number of species
  nspecies=2

  ! --- Species creation
  IF (.NOT. l_species_allocated) THEN
      ALLOCATE(species_parray(1:nspecies_max))
      l_species_allocated=.TRUE.
  ENDIF

  ! Electrons
  curr => species_parray(1)
  curr%charge = -echarge
  curr%mass = emass
  curr%nppcell = 40
  curr%x_min = xmin
  curr%x_max = xmax
  curr%y_min = ymin
  curr%y_max = ymax
  curr%z_min = zmin
  curr%z_max = zmax
  curr%vdrift_x =0._num
  curr%vdrift_y =0._num
  curr%vdrift_z =0._num
  curr%vth_x =0._num
  curr%vth_y =0._num
  curr%vth_z =0._num
  curr%sorting_period = 0
  curr%sorting_start = 0
  curr%species_npart=0

  ! Positrons
  curr => species_parray(2)
  curr%charge = echarge
  curr%mass = emass
  curr%nppcell = 40
  curr%x_min = xmin
  curr%x_max = xmax
  curr%y_min = ymin
  curr%y_max = ymax
  curr%z_min = zmin
  curr%z_max = zmax
  curr%vdrift_x =0._num
  curr%vdrift_y =0._num
  curr%vdrift_z =0._num
  curr%vth_x =0._num
  curr%vth_y =0._num
  curr%vth_z =0._num
  curr%sorting_period = 0
  curr%sorting_start = 0
  curr%species_npart=0

  ! --- mpi init communicator
  CALL mpi_minimal_init

  ! --- Check domain decomposition / Create Cartesian communicator / Allocate grid arrays
  CALL mpi_initialise

  ! --- Set tile split for particles
  CALL set_tile_split

  ! --- Allocate particle arrays for each tile of each species
  CALL init_tile_arrays

  ! --- Init "by hand" of the particle properties
  ALLOCATE(partpid(npid))
  clightsq=1/clight**2
  DO ispecies=1,nspecies
      curr=>species_parray(ispecies)
      jmin = NINT(MAX(curr%x_min-x_min_local,0.0_num)/dx)
      jmax = NINT(MIN(curr%x_max-x_min_local,x_max_local-x_min_local)/dx)
      kmin = NINT(MAX(curr%y_min-y_min_local,0.0_num)/dy)
      kmax = NINT(MIN(curr%y_max-y_min_local,y_max_local-y_min_local)/dy)
      lmin = NINT(MAX(curr%z_min-z_min_local,0.0_num)/dz)
      lmax = NINT(MIN(curr%z_max-z_min_local,z_max_local-z_min_local)/dz)
      DO l=lmin,lmax-1
          DO k=kmin,kmax-1
              DO j=jmin,jmax-1
                  DO ipart=1,curr%nppcell
                      CALL RANDOM_NUMBER(rng(1:6))
                      ! Sets positions and weight
                      partx = x_min_local+MIN(rng(1),0.999_num)*(x_max_local-x_min_local)
                      party = y_min_local+MIN(rng(2),0.999_num)*(y_max_local-y_min_local)
                      partz = z_min_local+MIN(rng(3),0.999_num)*(z_max_local-z_min_local)
                      partpid(wpid) = nc*dx*dy*dz/(curr%nppcell)
                      ! Sets velocity
                      up=rng(4)*200

                      th = -0.5*pi + rng(5)*pi
                      phi = rng(6)*2*pi

                      partux = up*cos(th)*cos(phi)
                      partuy = up*cos(th)*sin(phi)
                      partuz = up*sin(th)

                      gaminv = 1./sqrt(1.0_num + (partux**2 + partuy**2 + partuz**2)*clightsq)
                      ! Adds particle to array of tiles of current species
                      CALL add_particle_to_species(curr, partx, party, partz, &
                      partux, partuy, partuz, gaminv, partpid)
                  END DO
              END DO
          END DO
      END DO
  END DO ! END LOOP ON SPECIES

  DEALLOCATE(partpid)
  ! Allocate array to check the results
  ALLOCATE(tilesumex(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumey(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumez(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumbx(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumby(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumbz(ntilez,ntiley,ntilex))

  ! --- Field init
  CALL RANDOM_NUMBER(ex)
  CALL RANDOM_NUMBER(ey)
  CALL RANDOM_NUMBER(ez)
  CALL RANDOM_NUMBER(bx)
  CALL RANDOM_NUMBER(by)
  CALL RANDOM_NUMBER(bz)

  ex = ex*Ef
  ey = ey*Ef
  ez = ez*Ef

  bx = bx*Bf
  by = by*Bf
  bz = bz*Bf

  dtcoef = 0.5
  dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))

  write(0,*) 'dx',dx,'dy',dy,'dz',dz,'dt',dt
  write(0,*) 'sum(ex)',sum(ex),'sum(ey)',sum(ey),'sum(ez)',sum(ez)
  ! ______________________________________________________________________________________
  ! Test of the different subroutines with tiling

  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  sumex = 0

  i = 1
  name(i) = 'pxr_gete3d_n_energy_conserving'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 2 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'gete3d_energy_conserving_scalar_1_1_1'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 1 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)


  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV4_1_1_1'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 0 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  ! _________________________________________
  ! Test of extra developer's functions
#if defined(DEV)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV1_1_1_1'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 6 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV3_1_1_1'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 7 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

#endif
  ! End test of extra developer's functions
  ! _________________________________________

  ! Computation of the relative error
  CALL compute_err(i,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,epsilon,passed)

  title = 'Results order 1'
  CALL display_statistics(title,i,name,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,t)

  ! __ Order 2 __________________
  write(0,*) ''

  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  sumex = 0

  i = 1
  name(i) = 'pxr_gete3d_n_energy_conserving'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 2 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'gete3d_energy_conserving_scalar_2_2_2'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 1 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV4_2_2_2'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 0 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  ! _________________________________________
  ! Test of extra developer's functions
#if defined(DEV)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV1_2_2_2'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 6 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV3_2_2_2'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 7 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

#endif
  ! End test of extra developer's functions
  ! _________________________________________

  ! Computation of the relative error
  CALL compute_err(i,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,epsilon,passed)

  title = 'Results order 2'
  CALL display_statistics(title,i,name,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,t)

  ! __ Order 3 __________________
  write(0,*) ''

  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  sumex = 0

  i = 1
  name(i) = 'pxr_gete3d_n_energy_conserving'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 2 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'gete3d_energy_conserving_scalar_3_3_3'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 1 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV3_3_3_3'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 0 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  ! _________________________________________
  ! Test of extra developer's functions
#if defined(DEV)

  i = i+1
  name(i) = 'gete3d_energy_conserving_vec2_3_3_3'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 5 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_blockvec2_3_3_3'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 4 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vec_3_3_3'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 6 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

  i = i+1
  name(i) = 'geteb3d_energy_conserving_vecV2_3_3_3'
  write(0,*) 'Computation of ',name(i)
  fieldgathe = 7 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  t(i) = MPI_WTIME() - t0
  CALL check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)

#endif
  ! End test of extra developer's functions
  ! _________________________________________

  ! Computation of the relative error
  CALL compute_err(i,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,epsilon,passed)

  title = 'Results order 3'
  CALL display_statistics(title,i,name,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,t)

  ! ___ Final exam ____________________________________________
  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST TILING FIELD GATHERING 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST TILING FIELD GATHERING 3D FAILED **********  \e[0m \n"')
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')

  ! ______________________________________________________________________________________

  CALL mpi_close

END PROGRAM

! ______________________________________________________________________________________
! External subroutines


SUBROUTINE check_field_gathering(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz)
  USE particles
  USE constants
  USE tiling
  IMPLICIT NONE

  ! ___ Parameter declaration ________________________________________
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumex,tilesumey,tilesumez
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumbx,tilesumby,tilesumbz
  INTEGER(idp)                    :: ispecies, ix, iy, iz, count
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile

  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr=>species_parray(1)
        curr_tile=>curr%array_of_tiles(ix,iy,iz)

        DO ispecies=1, nspecies ! LOOP ON SPECIES
            ! - Get current tile properties
            ! - Init current tile variables
          curr=>species_parray(ispecies)
          curr_tile=>curr%array_of_tiles(ix,iy,iz)
          count=curr_tile%np_tile(1)

          tilesumex(iz,iy,ix)=SUM(curr_tile%part_ex(1:count))
          tilesumey(iz,iy,ix)=SUM(curr_tile%part_ey(1:count))
          tilesumez(iz,iy,ix)=SUM(curr_tile%part_ez(1:count))
          tilesumbx(iz,iy,ix)=SUM(curr_tile%part_bx(1:count))
          tilesumby(iz,iy,ix)=SUM(curr_tile%part_by(1:count))
          tilesumbz(iz,iy,ix)=SUM(curr_tile%part_bz(1:count))

        END DO! END LOOP ON SPECIES
      END DO
    END DO
  END DO! END LOOP ON TILES
END SUBROUTINE

SUBROUTINE compute_err(n,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,epsilon,passed)
  USE constants
  IMPLICIT NONE

  INTEGER(isp)                             :: n
  INTEGER(isp)                             :: i
  REAL(num)                                :: epsilon
  LOGICAL(lp), INTENT(INOUT)                   :: passed
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10), INTENT(INOUT)  :: errex,errey,errez
  REAL(num), dimension(10), INTENT(INOUT)  :: errbx,errby,errbz

  IF (n.gt.1) THEN
    DO i = 2,n
      errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
      errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
      errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
      errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
      errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
      errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)

      IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))

    ENDDO
  ENDIF

END SUBROUTINE


SUBROUTINE display_statistics(title,n,name,sumex,sumey,sumez,sumbx,sumby,sumbz, &
           errex,errey,errez,errbx,errby,errbz,t)
  USE constants
  IMPLICIT NONE

  CHARACTER(len=64)                        :: title
  INTEGER(isp)                             :: n
  REAL(num), dimension(10)                 :: t
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  INTEGER(isp)                             :: i

  write(0,*)
  write(0,'(A40)') title
  write(0,'(A40)') 'Electric field'
  write(0,'(A40, 7(A13))') "Subroutines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,7(X,E12.5))') name(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i), t(i)
  ENDDO
  write(0,'(A40)') 'Magnetic field'
  write(0,'(A40, 7(A13))') "Subroutines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,7(X,E12.5))') name(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i), t(i)
  ENDDO
END SUBROUTINE

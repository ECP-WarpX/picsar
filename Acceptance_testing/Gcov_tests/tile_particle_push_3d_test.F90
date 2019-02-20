! ________________________________________________________________________________________
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
! TILE_PARTICLE_PUSH_3D_TEST.F90
!
! Test code for the field gathering + particle pusher with tiling in 3D
!
! Mathieu Lobet, 2016.08
! ________________________________________________________________________________________

PROGRAM tile_particle_push_3d_test
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

  TYPE(particle_species), POINTER          :: curr
  TYPE(particle_species)                   :: curr0
  INTEGER(idp)                             :: jmin, jmax, kmin, kmax, lmin, lmax
  REAL(num)                                :: partx, party, partz
  REAL(num)                                :: partux, partuy, partuz, gaminv
  REAL(num), DIMENSION(:), ALLOCATABLE     :: partpid
  REAL(num)                                :: phi, th, up, clightsq
  REAL(num)                                :: Ef, Bf
  REAL(num), DIMENSION(6)                  :: rng=0_num
  REAL(num)                                :: epsilon
  REAL(num)                                :: t0
  LOGICAL(lp)                              :: passed
  REAL(num), dimension(10)                 :: tfg,tpp
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumex,tilesumey,tilesumez
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumbx,tilesumby,tilesumbz
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumx,tilesumy,tilesumz
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumpx,tilesumpy,tilesumpz
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: sumx,sumy,sumz
  REAL(num), dimension(10)                 :: sumpx,sumpy,sumpz
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  REAL(num), dimension(10)                 :: errx,erry,errz
  REAL(num), dimension(10)                 :: errpx,errpy,errpz
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
  nx_global_grid=25
  ny_global_grid=25
  nz_global_grid=25

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

  ! --- Order (these values will be changed later in the code, 
  ! this is the maximum in thi test)
  nox=3
  noy=3
  noz=3

  ! --- Init guard cells
  ! (Since we will use the same arrays for all orders, 
  ! we put the maximum of guard cells)
  nxguards=3
  nyguards=3
  nzguards=3
  nxjguards=3
  nyjguards=3
  nzjguards=3

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
  ALLOCATE(tilesumx(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumy(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumz(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumpx(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumpy(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumpz(ntilez,ntiley,ntilex))
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

  curr0 = curr
  ! ______________________________________________________________________________________
  ! Test of the different subroutines with tiling

  ! ______________
  ! Order 1

  nox=1 ; noy=1 ; noz=1

  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  errx = 0
  erry = 0
  errz = 0
  errpx = 0
  errpy = 0
  errpz = 0
  tfg = 0
  tpp = 0

  i = 1
  name(i) = 'field_gathering_sub + particle_pusher_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 2
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tfg(i) = MPI_WTIME() - t0
  ! Particle pusher then
  t0 = MPI_WTIME()
  CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_sub + particle_pusher_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0
  nox=1 ; noy=1 ; noz=1
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tfg(i) = MPI_WTIME() - t0
  ! Particle pusher then
  t0 = MPI_WTIME()
  CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 0
  nox=1 ; noy=1 ; noz=1
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0
  nox=1 ; noy=1 ; noz=1
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_cacheblock_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 0
  nox=1 ; noy=1 ; noz=1
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_cacheblock_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0
  nox=1 ; noy=1 ; noz=1
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  ! Computation of the relative error
  CALL compute_err(i,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  epsilon,passed)

  title = 'Results order 1'
  CALL display_statistics(title,i,name,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  tfg,tpp)

  ! __ Order 2 __________________
  write(0,*) ''

  nox=2 ; noy=2 ; noz=2

  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  errx = 0
  erry = 0
  errz = 0
  errpx = 0
  errpy = 0
  errpz = 0
  tfg = 0
  tpp = 0
  i = 0
  curr = curr0

  i = 1
  name(i) = 'field_gathering_sub + particle_pusher_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 2
  nox=2 ; noy=2 ; noz=2
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tfg(i) = MPI_WTIME() - t0
  ! Particle pusher then
  t0 = MPI_WTIME()
  CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_sub + particle_pusher_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0
  nox=2 ; noy=2 ; noz=2
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tfg(i) = MPI_WTIME() - t0
  ! Particle pusher then
  t0 = MPI_WTIME()
  CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 0
  nox=2 ; noy=2 ; noz=2
  ! field gathering and particle pusher
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0
  nox=2 ; noy=2 ; noz=2
  ! field gathering and particle pusher
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 0
  nox=2 ; noy=2 ; noz=2
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0
  nox=2 ; noy=2 ; noz=2
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  ! Computation of the relative error
  CALL compute_err(i,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  epsilon,passed)

  title = 'Results order 2'
  CALL display_statistics(title,i,name,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  tfg,tpp)

  ! __ Order 3 __________________
  write(0,*) ''

  nox=3 ; noy=3 ; noz=3

  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  errx = 0
  erry = 0
  errz = 0
  errpx = 0
  errpy = 0
  errpz = 0
  i = 0
  curr = curr0

  i = 1
  name(i) = 'field_gathering_sub + particle_pusher_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 2  
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tfg(i) = MPI_WTIME() - t0
  ! Particle pusher then
  t0 = MPI_WTIME()
  CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_sub + particle_pusher_sub(JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 1
  nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tfg(i) = MPI_WTIME() - t0
  ! Particle pusher then
  t0 = MPI_WTIME()
  CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub, scalar (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 1 ; nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub, scalar (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 1 ; nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub, vectorized (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 0 ; nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_sub, vectorized (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0 ; nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_cacheblock_sub (Boris)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 0
  fieldgathe = 0 ; nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  curr = curr0

  i = i+1
  name(i) = 'field_gathering_plus_particle_pusher_cacheblock_sub (JLVay)'
  write(0,*) 'Computation of ',name(i)
  particle_pusher = 1
  fieldgathe = 0 ; nox=3 ; noy=3 ; noz=3
  ! field gathering first
  t0 = MPI_WTIME()
  CALL field_gathering_plus_particle_pusher_cacheblock_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
  nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)
  tpp(i) = MPI_WTIME() - t0
  ! Checking
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz,&
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)

  ! Computation of the relative error
  CALL compute_err(i,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  epsilon,passed)

  title = 'Results order 3'
  CALL display_statistics(title,i,name,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  tfg,tpp)

  ! ___ Final exam ____________________________________________
  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST TILING FIELD GATHERING + PARTICLE PUSHER 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST TILING FIELD GATHERING + PARTICLE PUSHER 3D FAILED **********  \e[0m \n"')
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')

  ! ______________________________________________________________________________________

  CALL mpi_close

END PROGRAM

! ______________________________________________________________________________________
! External subroutines


SUBROUTINE check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz)
  USE particles
  USE constants
  USE tiling
  IMPLICIT NONE

  ! ___ Parameter declaration ________________________________________
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumex,tilesumey,tilesumez
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumbx,tilesumby,tilesumbz
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumx,tilesumy,tilesumz
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumpx,tilesumpy,tilesumpz
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
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

          tilesumx(iz,iy,ix)=SUM(curr_tile%part_x(1:count))
          tilesumy(iz,iy,ix)=SUM(curr_tile%part_y(1:count))
          tilesumz(iz,iy,ix)=SUM(curr_tile%part_z(1:count))
          tilesumpx(iz,iy,ix)=SUM(curr_tile%part_ux(1:count))
          tilesumpy(iz,iy,ix)=SUM(curr_tile%part_uy(1:count))
          tilesumpz(iz,iy,ix)=SUM(curr_tile%part_uz(1:count))

        END DO! END LOOP ON SPECIES
      END DO
    END DO
  END DO! END LOOP ON TILES
END SUBROUTINE

SUBROUTINE compute_err(n,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  epsilon,passed)
  USE constants
  IMPLICIT NONE

  INTEGER(isp)                             :: n
  INTEGER(isp)                             :: i
  REAL(num)                                :: epsilon
  LOGICAL(lp), INTENT(INOUT)                   :: passed
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: sumx,sumy,sumz
  REAL(num), dimension(10)                 :: sumpx,sumpy,sumpz
  REAL(num), dimension(10),INTENT(INOUT)   :: errex,errey,errez
  REAL(num), dimension(10),INTENT(INOUT)   :: errbx,errby,errbz
  REAL(num), dimension(10),INTENT(INOUT)   :: errx,erry,errz
  REAL(num), dimension(10),INTENT(INOUT)   :: errpx,errpy,errpz

  IF (n.gt.1) THEN
    DO i = 2,n
      errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
      errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
      errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
      errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
      errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
      errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)

      errx(i) = abs((sumx(i) - sumx(1)))/sumx(1)
      erry(i) = abs((sumy(i) - sumy(1)))/sumy(1)
      errz(i) = abs((sumz(i) - sumz(1)))/sumz(1)
      errpx(i) = abs((sumpx(i) - sumpx(1)))/sumpx(1)
      errpy(i) = abs((sumpy(i) - sumpy(1)))/sumpy(1)
      errpz(i) = abs((sumpz(i) - sumpz(1)))/sumpz(1)

      IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (erry(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errz(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errpx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errpy(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errpz(i) .gt. epsilon) passed = (passed.and.(.false.))

    ENDDO
  ENDIF

END SUBROUTINE


SUBROUTINE display_statistics(title,n,name,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz, &
  tfg,tpp)
  USE constants
  IMPLICIT NONE

  CHARACTER(len=64)                        :: title
  INTEGER(isp)                             :: n
  REAL(num), dimension(10)                 :: tfg,tpp
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: sumx,sumy,sumz
  REAL(num), dimension(10)                 :: sumpx,sumpy,sumpz
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  REAL(num), dimension(10)                 :: errx,erry,errz
  REAL(num), dimension(10)                 :: errpx,errpy,errpz
  INTEGER(isp)                             :: i

  write(0,*)
  write(0,'(A40)') title
  write(0,'(A60)') 'Field gathering, electric field'
  write(0,'(A60, 7(A13))') "Subroutines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A60,7(X,E12.5))') name(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i), tfg(i)
  ENDDO

  write(0,'(A60)') 'Field gathering, electric field, magnetic field'
  write(0,'(A60, 7(A13))') "Subroutines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
  write(0,'(A60,7(X,E12.5))') name(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i), tfg(i)
  ENDDO

  write(0,'(A40)') 'Particle pusher'
  write(0,'(A60, 7(A13))') "Subroutines", "sum(x)", "sum(y)", "sum(z)", "err x", "err y", "err z", "time (s)"
  write(0,'(A60, 7(A13))') "Subroutines", "sum(px)", "sum(py)", "sum(pz)", "err px", "err py", "err pz"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A60,7(X,E12.5))') name(i), sumx(i), sumy(i), sumz(i), errx(i), erry(i), errz(i), tpp(i)
    write(0,'(A60,7(X,E12.5))') '', sumpx(i), sumpy(i), sumpz(i), errpx(i), errpy(i), errpz(i), tpp(i)
  ENDDO

END SUBROUTINE

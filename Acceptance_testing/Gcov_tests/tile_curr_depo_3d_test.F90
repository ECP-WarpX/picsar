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
! TILE_CURR_DEPO_3D_TEST.F90
!
! Test code for the current deposition with tiling in 3D
!
! Mathieu Lobet, 2016.08
! ______________________________________________________________________________

PROGRAM tile_curr_depo_3d_test
  USE constants
  USE fields
  USE particles
  USE params
  USE shared_data
  USE mpi_routines
  USE control_file
  USE tiling

  ! ____________________________________________________________________________
  ! Parameters

  TYPE(particle_species), POINTER :: curr
  INTEGER(idp)                    :: jmin, jmax, kmin, kmax, lmin, lmax
  REAL(num)                                :: partx, party, partz
  REAL(num)                                :: partux, partuy, partuz, gaminv
  REAL(num), DIMENSION(:), ALLOCATABLE     :: partpid
  REAL(num)                                :: phi, th, up, clightsq
  REAL(num)                                :: Ef, Bf
  REAL(num), DIMENSION(6)                  :: rng=0_num
  REAL(num)                                :: epsilon
  REAL(num)                                :: t0
  LOGICAL(lp)                              :: passed
  REAL(num), dimension(10)                 :: t
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumjx,tilesumjy,tilesumjz
  REAL(num), dimension(10)                 :: sumjx,sumjy,sumjz
  REAL(num), dimension(10)                 :: errjx,errjy,errjz
  CHARACTER(len=64)                        :: title
  ! ____________________________________________________________________________
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
  lvec_curr_depo = 8

  ! --- Interpolation
  l_lower_order_in_v = .TRUE.

  ! --- Field properties
  Ef = 1E12 !V/m
  Bf = 1E4  !T

  ! --- Relative error to pass the test
  epsilon = 1e-5_num
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
  ALLOCATE(tilesumjx(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumjy(ntilez,ntiley,ntilex))
  ALLOCATE(tilesumjz(ntilez,ntiley,ntilex))

  dtcoef = 0.5
  dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))

  write(0,*) 'dx',dx,'dy',dy,'dz',dz,'dt',dt
  ! ____________________________________________________________________________
  ! Test of the different subroutines with tiling

  ! __ Order 1 __________________
  write(0,*) ''

  errjx = 0
  errjy = 0
  errjz = 0

  i = 1
  name(i) = 'Esirkepov order n sequential version'
  currdepo = 2 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical scalar order 1 sequential version'
  currdepo = 5 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical scalar order 1 openmp version'
  currdepo = 4 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical vectorized order 1 openmp version'
  currdepo = 3 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Esirkepov order 1 openmp version'
  currdepo = 1 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

!   i = i+1
!   name(i) = 'Esirkepov vectorized order 1 openmp version'
!   currdepo = 0 ; nox=1 ; noy=1 ; noz=1
!   t0 = MPI_WTIME()
!   CALL pxrdepose_currents_on_grid_jxjyjz
!   t(i) = MPI_WTIME() - t0
!   CALL check(tilesumjx,tilesumjy,tilesumjz)
!   sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  ! Computation of the relative error
  CALL compute_err(i,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,epsilon,passed)

  title = 'Results order 1'
  CALL display_statistics(title,i,name,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,t)

  ! __ Order 2 __________________
  write(0,*) ''

  errjx = 0
  errjy = 0
  errjz = 0
  tilesumjx = 0
  tilesumjy = 0
  tilesumjz = 0

  i = 1
  name(i) = 'Esirkepov order n sequential version'
  currdepo = 2 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical scalar order 2 sequential version'
  currdepo = 5 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical scalar order 2 openmp version'
  currdepo = 4 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical vectorized order 2 openmp version'
  currdepo = 3 ;  nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Esirkepov order 2 openmp version'
  currdepo = 1 ;  nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

!   i = i+1
!   name(i) = 'Esirkepov vectorized order 2 openmp version'
!   currdepo = 0 ;  nox=2 ; noy=2 ; noz=2
!   t0 = MPI_WTIME()
!   CALL pxrdepose_currents_on_grid_jxjyjz
!   t(i) = MPI_WTIME() - t0
!   CALL check(tilesumjx,tilesumjy,tilesumjz)
!   sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  ! Computation of the relative error
  CALL compute_err(i,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,epsilon,passed)

  title = 'Results order 2'
  CALL display_statistics(title,i,name,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,t)

  ! __ Order 3 __________________
  write(0,*) ''

  errjx = 0
  errjy = 0
  errjz = 0
  tilesumjx = 0
  tilesumjy = 0
  tilesumjz = 0

  i = 1
  name(i) = 'Esirkepov order n sequential version'
  write(0,*) 'Computation of ',name(i)
  currdepo = 2 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical scalar order 3 sequential version'
  write(0,*) 'Computation of ',name(i)
  currdepo = 5 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical scalar order 3 openmp version'
  write(0,*) 'Computation of ',name(i)
  currdepo = 4 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Classical vectorized order 3 openmp version'
  write(0,*) 'Computation of ',name(i)
  currdepo = 3 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  i = i+1
  name(i) = 'Esirkepov order 3 openmp version'
  write(0,*) 'Computation of ',name(i)
  currdepo = 1 ;  nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_currents_on_grid_jxjyjz
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumjx,tilesumjy,tilesumjz)
  sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

!   i = i+1
!   name(i) = 'Esirkepov vectorized order 2 openmp version'
!   currdepo = 0 ;  nox=2 ; noy=2 ; noz=2
!   t0 = MPI_WTIME()
!   CALL pxrdepose_currents_on_grid_jxjyjz
!   t(i) = MPI_WTIME() - t0
!   CALL check(tilesumjx,tilesumjy,tilesumjz)
!   sumjx(i) = SUM(tilesumjx) ; sumjy(i) = SUM(tilesumjy) ; sumjz(i) = SUM(tilesumjz)

  ! Computation of the relative error
  CALL compute_err(i,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,epsilon,passed)

  title = 'Results order 3'
  CALL display_statistics(title,i,name,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,t)

  ! ___ Final exam ____________________________________________
  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST TILING CURRENT DEPOSITION 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST TILING CURRENT DEPOSITION 3D FAILED **********  \e[0m \n"')
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')

  CALL MPI_FINALIZE(errcode)

  ! ____________________________________________________________________________

END PROGRAM

! ______________________________________________________________________________
! External subroutines

SUBROUTINE depose_currents_on_grid_jxjyjz
! ______________________________________________________________________________
  USE fields
  USE shared_data
  USE params
  USE time_stat

  IMPLICIT NONE
  REAL(num) :: tdeb


  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE

      ! ____________________________________________________________________________________
      ! Generic current deposition routine
       SUBROUTINE depose_jxjyjz(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
             dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
             nox,noy,noz,current_depo_algo) !#do not parse
         USE constants
         IMPLICIT NONE
         INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard,current_depo_algo
         REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
         REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
         REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
       END SUBROUTINE

       ! Interface for subroutine with no reduction - classical deposition order 1
        SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
             uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
             dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz,lvect) !#do not parse
          USE constants
          IMPLICIT NONE
          INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
          INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
          REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
          REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
          REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
          INTEGER(idp) :: ncx, ncy, ncz
          INTEGER(idp) :: lvect
        END SUBROUTINE

       ! Interface for subroutine with no reduction - classical deposition order 2
        SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
             uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
             dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz,lvect) !#do not parse
          USE constants
          IMPLICIT NONE
          INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
          INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
          REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
          REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
          REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
          INTEGER(idp) :: ncx, ncy, ncz
          INTEGER(idp) :: lvect
        END SUBROUTINE

       ! Interface for subroutine with no reduction - classical deposition order 3
        SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
             uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
             dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz,lvect) !#do not parse
          USE constants
          IMPLICIT NONE
          INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
          INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
          REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
          REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
          REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
          INTEGER(idp) :: ncx, ncy, ncz
          INTEGER(idp) :: lvect
        END SUBROUTINE

      SUBROUTINE current_reduction_1_1_1(jx,jy,jz,jxcells,jycells,jzcells,ncells, &
                 nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
        USE constants
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                 :: nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
        INTEGER(idp), INTENT(IN)                 :: nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN), DIMENSION(8,ncells):: jxcells,jycells,jzcells
      END SUBROUTINE

      SUBROUTINE current_reduction_2_2_2(jx,jy,jz,jxcells,jycells,jzcells,ncells,&
                 nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
        USE constants
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                 :: nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
        INTEGER(idp), INTENT(IN)                 :: nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN), DIMENSION(8,ncells):: jxcells,jycells,jzcells
      END SUBROUTINE

      SUBROUTINE current_reduction_3_3_3(jx,jy,jz,jxcells,jycells,jzcells,ncells,&
                 nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
        USE constants
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                 :: nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                 :: ncx, ncy, ncz
        INTEGER(idp), INTENT(IN)                 :: nxguard,nyguard,nzguard
        REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
        REAL(num),INTENT(IN), DIMENSION(8,ncells):: jxcells,jycells,jzcells
      END SUBROUTINE

    END INTERFACE
  ! ___________________________________________________________________________

#if defined(DEBUG)
  WRITE(0,*) "Depose_currents_on_grid: start"
#endif

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==2
  CALL start_vtune_collection()
#endif
#if SDE==2
  CALL start_vtune_collection()
#endif

  jx = 0.0_num
  jy = 0.0_num
  jz = 0.0_num

  ! Current deposition branches

  ! _______________________________________________________
  ! Classical current deposition, non-optimized/no tiling


    ! _______________________________________________________
    ! Classical current deposition, non-optimized/no tiling
    IF (currdepo.EQ.5) THEN

       IF ((nox.eq.noy).AND.(noy.eq.noz)) THEN
          CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(depose_jxjyjz, &
               jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt, 3_idp)
          ! The last argument is 3: this means the scalar routines will be used inside `depose_jxjyjz`
       ELSE
          CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, &
               jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
               nox,noy,noz,dx,dy,dz,dt,1_idp)
          ! The last argument is 1: this means the generic esirkepov routine will be used inside `depose_jxjyjz`
       ENDIF


    ! _______________________________________________________
    ! Classical current deposition, non-optimized/tiling
    ELSE IF (currdepo.EQ.4) THEN

       IF ((nox.eq.noy).AND.(noy.eq.noz)) THEN
          CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz, &
               jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt, 3_idp)
          ! The last argument is 3: this means the scalar routines will be used inside `depose_jxjyjz`
       ELSE
          CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, &
               jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
               nox,noy,noz,dx,dy,dz,dt,1_idp)
          ! The last argument is 1: this means the generic esirkepov routine will be used inside `depose_jxjyjz`
       ENDIF

    ! _______________________________________________________
    ! Classical current deposition, parallel, vectorized
    ELSE IF (currdepo.EQ.3) THEN

      IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
        ! Old version with reduction for each species
        !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_vecHVv3_3_3_3, &
        !jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
        CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( &
             depose_jxjyjz_vecHV_vnr_3_3_3, current_reduction_3_3_3,&
             jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
      ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN
        ! Old version with reduction for each species
        !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_vecHVv2_2_2_2, &
        !jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
        CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( &
             depose_jxjyjz_vecHV_vnr_2_2_2, current_reduction_2_2_2,&
             jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
      ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN
        ! Old version with reduction for each species
        !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_vecHVv2_1_1_1, &
        !jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
        CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( &
             depose_jxjyjz_vecHV_vnr_1_1_1, current_reduction_1_1_1,&
             jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
        !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v3( &
        !     depose_jxjyjz_vecHV_vnr_1_1_1, current_reduction_1_1_1,&
        !     jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
      ELSE
        CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, &
             jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
             nox,noy,noz,dx,dy,dz,dt,1_idp)
      ENDIF
    ! _______________________________________________________
    ! Esirkepov sequential version
    ELSE IF (currdepo.EQ.2) THEN

      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
           nox,noy,noz,dx,dy,dz,dt)

    ! _______________________________________________________
    ! Esirkepov tiling version
    ELSE IF (currdepo.EQ.1) THEN

       CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, &
            jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
            nox,noy,noz,dx,dy,dz,dt, 0_idp)
          ! The last argument is 0: this means the optimized esirkepov routine will be used inside `depose_jxjyjz`

    ! _______________________________________________________
    ! Default - Esirkepov parallel version with OPENMP/tiling and optimizations
    ELSE IF (currdepo .EQ. 0) THEN

       CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, &
            jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
            nox,noy,noz,dx,dy,dz,dt, 0_idp)
       ! The last argument is 1: this means the optimized esirkepov routine will be used inside `depose_jxjyjz`

      ! Arbitrary order
    ELSE

       CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz, &
            jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
            nox,noy,noz,dx,dy,dz,dt, 1_idp)
       ! The last argument is 1: this means the generic esirkepov routine will be used inside `depose_jxjyjz`
    ENDIF

END SUBROUTINE depose_currents_on_grid_jxjyjz

SUBROUTINE check(tilesumjx,tilesumjy,tilesumjz)
  USE particles
  USE constants
  USE tiling
  IMPLICIT NONE

  ! ___ Parameter declaration ________________________________________
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumjx,tilesumjy,tilesumjz
  INTEGER(idp)                    :: ix, iy, iz
  TYPE(grid_tile), POINTER        :: currg

  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex

        currg=>aofgrid_tiles(ix,iy,iz)

        tilesumjx(iz,iy,ix)=SUM(currg%arr1)
        tilesumjy(iz,iy,ix)=SUM(currg%arr2)
        tilesumjz(iz,iy,ix)=SUM(currg%arr3)

      END DO
    END DO
  END DO! END LOOP ON TILES
END SUBROUTINE

SUBROUTINE compute_err(n,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,epsilon,passed)
  USE constants
  IMPLICIT NONE

  INTEGER(isp)                             :: n
  INTEGER(isp)                             :: i
  REAL(num)                                :: epsilon
  LOGICAL(lp), INTENT(INOUT)                   :: passed
  REAL(num), dimension(10)                 :: sumjx,sumjy,sumjz
  REAL(num), dimension(10), INTENT(INOUT)  :: errjx,errjy,errjz

  IF (n.gt.1) THEN
    DO i = 2,n
      errjx(i) = abs((sumjx(i) - sumjx(1)))/sumjx(1)
      errjy(i) = abs((sumjy(i) - sumjy(1)))/sumjy(1)
      errjz(i) = abs((sumjz(i) - sumjz(1)))/sumjz(1)

      IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))

    ENDDO
  ENDIF

END SUBROUTINE


SUBROUTINE display_statistics(title,n,name,sumjx,sumjy,sumjz, &
           errjx,errjy,errjz,t)
  USE constants
  IMPLICIT NONE

  CHARACTER(len=64)                        :: title
  INTEGER(isp)                             :: n
  REAL(num), dimension(10)                 :: t
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), dimension(10)                 :: sumjx,sumjy,sumjz
  REAL(num), dimension(10)                 :: errjx,errjy,errjz
  INTEGER(isp)                             :: i

  write(0,*)
  write(0,'(A40)') title
  write(0,'(A40)') 'Current field'
  write(0,'(A50, 7(A13))') "Subroutines", "sum(jx)", "sum(jy)", "sum(jz)", "err jx", "err jy", "err jz", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A50,7(X,E12.5))') name(i), sumjx(i), sumjy(i), sumjz(i), errjx(i), errjy(i), errjz(i), t(i)
  ENDDO
END SUBROUTINE

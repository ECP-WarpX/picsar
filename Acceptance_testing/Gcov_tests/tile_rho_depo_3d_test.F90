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
! TILE_RHO_DEPO_3D_TEST.F90
!
! Test code for the charge deposition with tiling in 3D
!
! Mathieu Lobet, 2016.08
! ______________________________________________________________________________

PROGRAM tile_rho_depo_3d_test
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

  TYPE(particle_species), POINTER          :: curr
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
  REAL(num), dimension(10)                 :: t
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumrho
  REAL(num), dimension(10)                 :: sumrho
  REAL(num), dimension(10)                 :: errrho
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

  ! --- Domain extension
  xmin=0
  ymin=0
  zmin=0
  xmax=1e-6
  ymax=1e-6
  zmax=1e-6

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

  ! --- Init particle tiling split
  ntilex = 6
  ntiley = 6
  ntilez = 6

  ! --- Vector length field gathering
  lvec_rho_depo = 64

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
  curr%x_max = xmax*0.5
  curr%y_min = ymin
  curr%y_max = ymax*0.5
  curr%z_min = zmin
  curr%z_max = zmax*0.5
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
  END DO ! end loop on speciess

  DEALLOCATE(partpid)
  ! Allocate array to check the results
  ALLOCATE(tilesumrho(ntilez,ntiley,ntilex))

  ! Time step using the CFL
  dtcoef = 0.5
  dt = dtcoef/(clight*sqrt(1.0_num/dx**2+1.0_num/dy**2+1.0_num/dz**2))

  write(0,*) 'dx',dx,'dy',dy,'dz',dz,'dt',dt
  ! ____________________________________________________________________________
  ! Test of the different subroutines with tiling

  ! __ Order 1 __________________
  write(0,*) ''

  errrho = 0

  ! Reference
  i = 1
  name(i) = 'Charge deposition order n'
  rhodepo = 2 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  i = i+1
  name(i) = 'Charge deposition scalar'
  rhodepo = 1 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  i = i+1
  name(i) = 'Charge deposition vectorized'
  rhodepo = 0 ; nox=1 ; noy=1 ; noz=1
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  ! Computation of the relative error
  CALL compute_err(i,sumrho, &
           errrho,epsilon,passed)

  title = 'Results order 1'
  CALL display_statistics(title,i,name,sumrho, &
           errrho,t)

  ! __ Order 2 __________________
  write(0,*) ''

  errrho = 0
  tilesumrho = 0

  i = 1
  name(i) = 'Charge deposition order n'
  rhodepo = 2 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  i = i+1
  name(i) = 'Charge deposition scalar'
  rhodepo = 1 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  i = i+1
  name(i) = 'Charge deposition vectorized'
  rhodepo = 0 ; nox=2 ; noy=2 ; noz=2
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)


  ! Computation of the relative error
  CALL compute_err(i,sumrho, &
           errrho,epsilon,passed)

  title = 'Results order 2'
  CALL display_statistics(title,i,name,sumrho, &
           errrho,t)

  ! __ Order 3 __________________
  write(0,*) ''

  errrho = 0
  tilesumrho = 0


  i = 1
  name(i) = 'Charge deposition order n'
  write(0,*) 'Computation of ',name(i)
  rhodepo = 2 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  i = i+1
  name(i) = 'Charge deposition scalar'
  write(0,*) 'Computation of ',name(i)
  rhodepo = 1 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)

  i = i+1
  name(i) = 'Charge deposition vectorized'
  write(0,*) 'Computation of ',name(i)
  rhodepo = 0 ; nox=3 ; noy=3 ; noz=3
  t0 = MPI_WTIME()
  CALL pxrdepose_rho_on_grid
  t(i) = MPI_WTIME() - t0
  CALL check(tilesumrho)
  sumrho(i) = SUM(tilesumrho)


  ! Computation of the relative error
  CALL compute_err(i,sumrho, &
           errrho,epsilon,passed)

  title = 'Results order 3'
  CALL display_statistics(title,i,name,sumrho, &
           errrho,t)

  ! ___ Final exam ____________________________________________
  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST TILING CHARGE DEPOSITION 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST TILING CHARGE DEPOSITION 3D FAILED **********  \e[0m \n"')
    ! Failure exit
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')

  CALL MPI_FINALIZE(errcode)
  ! ____________________________________________________________________________

END PROGRAM

! ______________________________________________________________________________
! External subroutines



SUBROUTINE check(tilesumrho)
  USE particles
  USE constants
  USE tiling
  IMPLICIT NONE

  ! ___ Parameter declaration ________________________________________
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumrho
  INTEGER(idp)             :: ix, iy, iz
  TYPE(grid_tile), POINTER        :: currg

  DO iz=1, ntilez ! LOOP ON TILES
    DO iy=1, ntiley
      DO ix=1, ntilex

        currg=>aofgrid_tiles(ix,iy,iz)

        tilesumrho(iz,iy,ix)=SUM(currg%arr1)

      END DO
    END DO
  END DO! END LOOP ON TILES
END SUBROUTINE

SUBROUTINE compute_err(n,sumrho, &
           errrho,epsilon,passed)
  USE constants
  IMPLICIT NONE

  INTEGER(isp)                             :: n
  INTEGER(isp)                             :: i
  REAL(num)                                :: epsilon
  LOGICAL(lp), INTENT(INOUT)                   :: passed
  REAL(num), dimension(10)                 :: sumrho
  REAL(num), dimension(10), INTENT(INOUT)  :: errrho

  IF (n.gt.1) THEN
    DO i = 2,n
      errrho(i) = abs((sumrho(i) - sumrho(1)))/abs(sumrho(1))

      IF (errrho(i) .gt. epsilon) passed = (passed.and.(.false.))

    ENDDO
  ENDIF

END SUBROUTINE


SUBROUTINE display_statistics(title,n,name,sumrho, &
           errrho,t)
  USE constants
  IMPLICIT NONE

  CHARACTER(len=64)                        :: title
  INTEGER(isp)                             :: n
  REAL(num), dimension(10)                 :: t
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), dimension(10)                 :: sumrho
  REAL(num), dimension(10)                 :: errrho
  INTEGER(isp)                             :: i

  write(0,*)
  write(0,'(A40)') title
  write(0,'(A40)') 'Rho field'
  write(0,'(A50, 7(A13))') "Subroutines", "sum(rho)", "err rho", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A50,7(X,E12.5))') name(i), sumrho(i), errrho(i), t(i)
  ENDDO
END SUBROUTINE

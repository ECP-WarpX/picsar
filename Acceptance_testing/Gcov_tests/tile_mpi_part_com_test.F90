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
! TILE_MPI_PART_COM_TEST.F90
!
! Test code for the particle communications
!
! Mathieu Lobet, 2016.08
! ______________________________________________________________________________

PROGRAM tile_mpi_part_com_test
	USE constants
	USE fields
	USE particles
	USE params
	USE shared_data
	USE mpi_routines
	USE control_file
	USE tiling
	USE particle_boundary

	! ____________________________________________________________________________
	! Parameters

	TYPE(particle_species), POINTER          :: curr
	TYPE(particle_species)                   :: curr0
	REAL(num)                                :: partx, party, partz
	REAL(num)                                :: partux, partuy, partuz, gaminv
	REAL(num), DIMENSION(:), ALLOCATABLE		 :: partpid
	REAL(num)                                :: phi, th, up, clightsq
	REAL(num)                                :: Ef, Bf
	REAL(num), DIMENSION(6)                  :: rng=0_num
  REAL(num)                                :: epsilon
  REAL(num)                                :: t0
  LOGICAL(lp)                              :: passed
  REAL(num), dimension(10)                 :: ttile,tmpi
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumex,tilesumey,tilesumez
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumbx,tilesumby,tilesumbz
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumx,tilesumy,tilesumz
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumpx,tilesumpy,tilesumpz
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: tilesumga
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: sumx,sumy,sumz
  REAL(num), dimension(10)                 :: sumpx,sumpy,sumpz
  REAL(num), dimension(10)                 :: sumga
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  REAL(num), dimension(10)                 :: errx,erry,errz
  REAL(num), dimension(10)                 :: errpx,errpy,errpz
  REAL(num), dimension(10)                 :: errga

	! ____________________________________________________________________________
	! Initialization
	! --- default init
	CALL default_init

	! --- Dimension
	c_dim = 3

	! --- Number of processors
	nprocx=1
	nprocy=2
	nprocz=2

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
	ntilex = 10
	ntiley = 5
	ntilez = 5

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

	! --- Number of iterations
	nsteps = 30

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
	curr%nppcell = 20
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
	curr%nppcell = 20
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
											up=rng(4)*500

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
  ALLOCATE(tilesumga(ntilez,ntiley,ntilex))

  errx=0
  erry=0
  errz=0
  errpx=0
  errpy=0
  errpz=0
  errga=0

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

  ! Time set to 0
  tmpi = 0
  ttile = 0

  ! Check initial
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz,tilesumga)
  sumex(1) = SUM(tilesumex) ; sumey(1) = SUM(tilesumey) ; sumez(1) = SUM(tilesumez)
  sumbx(1) = SUM(tilesumbx) ; sumby(1) = SUM(tilesumby) ; sumbz(1) = SUM(tilesumbz)
  sumx(1) = SUM(tilesumx) ; sumy(1) = SUM(tilesumy) ; sumz(1) = SUM(tilesumz)
  sumpx(1) = SUM(tilesumpx) ; sumpy(1) = SUM(tilesumpy) ; sumpz(1) = SUM(tilesumpz)
  sumga(1) = SUM(tilesumga)

  IF (rank.eq.0) THEN
    write(0,*) 'dx',dx,'dy',dy,'dz',dz,'dt',dt
    write(0,*)
    write(0,*) 'Initial State'
    write(0,'(X,"sum(x**2): ",E12.5)') sumx(1)
    write(0,'(X,"sum(y**2): ",E12.5)') sumy(1)
    write(0,'(X,"sum(z**2): ",E12.5)') sumz(1)
    write(0,'(X,"sum(px**2): ",E12.5)') sumpx(1)
    write(0,'(X,"sum(py**2): ",E12.5)') sumpy(1)
    write(0,'(X,"sum(pz**2): ",E12.5)') sumpz(1)
    write(0,'(X,"sum(gamma**2): ",E12.5)') sumga(1)
    write(0,*)
  ENDIF
  curr0 = curr
  ! ____________________________________________________________________________
  ! Test of the different subroutines with tiling

  i = 1
  name(i) = 'particle_bcs_tiles + particle_bcs_mpi_blocking'
  IF (rank.eq.0) write(0,*) 'Computation of ',name(i)

  nox=1 ; noy=1 ; noz=1


  DO it=1,nsteps

    ! Particle pusher
    CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
    nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

    ! Tile communications
    t0 = MPI_WTIME()
    CALL particle_bcs_tiles()
    ttile(i) = ttile(i) + MPI_WTIME() - t0
    ! MPI communications
    t0 = MPI_WTIME()
    CALL particle_bcs_mpi_blocking()
    tmpi(i) = tmpi(i) + MPI_WTIME() - t0

  ENDDO
  ! Check
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz,tilesumga)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)
  sumga(i) = SUM(tilesumga)



  curr = curr0
  i = i+1
  name(i) = 'particle_bcs_tiles + particle_bcs_mpi_nonblocking'
  IF (rank.eq.0) write(0,*) 'Computation of ',name(i)
  nox=1 ; noy=1 ; noz=1


  DO it=1,nsteps

    ! Particle pusher
    CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
    nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

    ! Tile communications
    t0 = MPI_WTIME()
    CALL particle_bcs_tiles()
    ttile(i) = ttile(i) + MPI_WTIME() - t0
    ! MPI communications
    t0 = MPI_WTIME()
    CALL particle_bcs_mpi_non_blocking()
    tmpi(i) = tmpi(i) + MPI_WTIME() - t0

  ENDDO
  ! Check
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz,tilesumga)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)
  sumga(i) = SUM(tilesumga)


  curr = curr0
  i = i+1
  name(i) = 'particle_bcs_tiles_openmp + particle_bcs_mpi_nonblocking'
  IF (rank.eq.0) write(0,*) 'Computation of ',name(i)
  nox=1 ; noy=1 ; noz=1


  DO it=1,nsteps

    ! Particle pusher
    CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
    nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

    ! Tile communications
    t0 = MPI_WTIME()
    CALL particle_bcs_tiles_openmp()
    ttile(i) = ttile(i) + MPI_WTIME() - t0
    ! MPI communications
    t0 = MPI_WTIME()
    CALL particle_bcs_mpi_non_blocking()
    tmpi(i) = tmpi(i) + MPI_WTIME() - t0

  ENDDO
  ! Check
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz,tilesumga)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)
  sumga(i) = SUM(tilesumga)

  curr = curr0
  i = i+1
  name(i) = 'particle_bcs_tiles_and_mpi_3d'
  IF (rank.eq.0) write(0,*) 'Computation of ',name(i)
  nox=1 ; noy=1 ; noz=1


  DO it=1,nsteps

    ! Particle pusher
    CALL particle_pusher_sub(ex,ey,ez,bx,by,bz,nx,ny,nz,nxguards,nyguards, &
    nzguards,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,l_lower_order_in_v)

    ! Tile and MPI communications
    t0 = MPI_WTIME()
    CALL particle_bcs_tiles_and_mpi_3d()
    ttile(i) = ttile(i) + MPI_WTIME() - t0

  ENDDO
  ! Check
  CALL check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz,tilesumga)
  sumex(i) = SUM(tilesumex) ; sumey(i) = SUM(tilesumey) ; sumez(i) = SUM(tilesumez)
  sumbx(i) = SUM(tilesumbx) ; sumby(i) = SUM(tilesumby) ; sumbz(i) = SUM(tilesumbz)
  sumx(i) = SUM(tilesumx) ; sumy(i) = SUM(tilesumy) ; sumz(i) = SUM(tilesumz)
  sumpx(i) = SUM(tilesumpx) ; sumpy(i) = SUM(tilesumpy) ; sumpz(i) = SUM(tilesumpz)
  sumga(i) = SUM(tilesumga)

  ! ________________________________________________________
  ! Computation of the relative error
  CALL compute_err(i,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz,sumga, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz,errga, &
  epsilon,passed)

  ! Display the results
  CALL display_statistics(i,name,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz,sumga, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz,errga, &
  tmpi,ttile)



  ! ___ Final exam ____________________________________________
  IF (rank.eq.0) THEN
    write(0,*)
    IF (passed) THEN
      !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
      !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
      CALL system('printf "\e[32m ********** TEST TILE AND MPI PARTICLE COMMUNICATIONS PASSED **********  \e[0m \n"')
    ELSE
      !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
      !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
      CALL system('printf "\e[31m ********** TEST TILE AND MPI PARTICLE COMMUNICATIONS FAILED **********  \e[0m \n"')
      CALL EXIT(9)
    ENDIF
  ENDIF

  IF (rank.eq.0) write(0,'(" ____________________________________________________________________________")')
  ! ____________________________________________________________________________

  CALL mpi_close

END PROGRAM

! ______________________________________________________________________________
! External subroutines


SUBROUTINE check(tilesumex,tilesumey,tilesumez,tilesumbx,tilesumby,tilesumbz, &
  tilesumx,tilesumy,tilesumz,tilesumpx,tilesumpy,tilesumpz,tilesumga)
  USE particles
  USE constants
  USE tiling
  IMPLICIT NONE

  ! ___ Parameter declaration ________________________________________
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumex,tilesumey,tilesumez
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumbx,tilesumby,tilesumbz
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumx,tilesumy,tilesumz
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumpx,tilesumpy,tilesumpz
  REAL(num), DIMENSION(ntilez,ntiley,ntilex) :: tilesumga
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

          tilesumx(iz,iy,ix)=SUM(curr_tile%part_x(1:count)**2)
          tilesumy(iz,iy,ix)=SUM(curr_tile%part_y(1:count)**2)
          tilesumz(iz,iy,ix)=SUM(curr_tile%part_z(1:count)**2)
          tilesumpx(iz,iy,ix)=SUM(curr_tile%part_ux(1:count)**2)
          tilesumpy(iz,iy,ix)=SUM(curr_tile%part_uy(1:count)**2)
          tilesumpz(iz,iy,ix)=SUM(curr_tile%part_uz(1:count)**2)

          tilesumga(iz,iy,ix)=SUM((1./curr_tile%part_gaminv(1:count))**2)

        END DO! END LOOP ON SPECIES
      END DO
    END DO
  END DO! END LOOP ON TILES
END SUBROUTINE

SUBROUTINE compute_err(n,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz,sumga, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz,errga, &
  epsilon,passed)
  USE constants
  USE mpi_routines
  USE shared_data
  IMPLICIT NONE

  INTEGER(isp)                             :: n
  INTEGER(isp)                             :: i
  REAL(num)                                :: epsilon
  LOGICAL(lp), INTENT(INOUT)                   :: passed
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: sumx,sumy,sumz
  REAL(num), dimension(10)                 :: sumpx,sumpy,sumpz,sumga
  REAL(num), dimension(10),INTENT(INOUT)   :: errex,errey,errez
  REAL(num), dimension(10),INTENT(INOUT)   :: errbx,errby,errbz
  REAL(num), dimension(10),INTENT(INOUT)   :: errx,erry,errz
  REAL(num), dimension(10),INTENT(INOUT)   :: errpx,errpy,errpz,errga

  REAL(num), dimension(10)                 :: sumxtot,sumytot,sumztot
  REAL(num), dimension(10)                 :: sumpxtot,sumpytot,sumpztot,sumgatot

  ! MPI REDUCTION
  CALL MPI_REDUCE(sumx(1),sumxtot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(sumy(1),sumytot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(sumz(1),sumztot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(sumpx(1),sumpxtot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(sumpy(1),sumpytot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(sumpz(1),sumpztot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(sumga(1),sumgatot(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)

  IF (n.gt.1) THEN
    DO i = 2,n
      !errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
      !errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
      !errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
      !errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
      !errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
      !errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)

      errx(i) = abs((sumxtot(i) - sumxtot(1)))/sumxtot(1)
      erry(i) = abs((sumytot(i) - sumytot(1)))/sumytot(1)
      errz(i) = abs((sumztot(i) - sumztot(1)))/sumztot(1)
      errpx(i) = abs((sumpxtot(i) - sumpxtot(1)))/sumpxtot(1)
      errpy(i) = abs((sumpytot(i) - sumpytot(1)))/sumpytot(1)
      errpz(i) = abs((sumpztot(i) - sumpztot(1)))/sumpztot(1)

      errga(i) = abs((sumgatot(i) - sumgatot(1)))/sumgatot(1)


      !IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
      !IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
      !IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))

      !IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
      !IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
      !IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (erry(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errz(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errpx(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errpy(i) .gt. epsilon) passed = (passed.and.(.false.))
      IF (errpz(i) .gt. epsilon) passed = (passed.and.(.false.))

      IF (errga(i) .gt. epsilon) passed = (passed.and.(.false.))

    ENDDO
  ENDIF

END SUBROUTINE


SUBROUTINE display_statistics(n,name,&
  sumex,sumey,sumez,sumbx,sumby,sumbz, &
  sumx,sumy,sumz,sumpx,sumpy,sumpz,sumga, &
  errex,errey,errez,errbx,errby,errbz, &
  errx,erry,errz,errpx,errpy,errpz,errga, &
  tmpi,ttile)
  USE constants
  USE shared_data
  IMPLICIT NONE

  INTEGER(isp)                             :: n
  REAL(num), dimension(10)                 :: ttile,tmpi
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: sumx,sumy,sumz
  REAL(num), dimension(10)                 :: sumpx,sumpy,sumpz,sumga
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  REAL(num), dimension(10)                 :: errx,erry,errz
  REAL(num), dimension(10)                 :: errpx,errpy,errpz,errga
  REAL(num), dimension(10)                 :: ttileave,tmpiave
  REAL(num), dimension(10)                 :: ttilemax,tmpimax
  REAL(num), dimension(10)                 :: ttilemin,tmpimin
  INTEGER(isp)                             :: i

  CALL MPI_REDUCE(ttile(1),ttileave(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(tmpi(1),tmpiave(1),n,mpidbl,MPI_SUM,0_isp,comm,errcode)
  CALL MPI_REDUCE(ttile(1),ttilemax(1),n,mpidbl,MPI_MAX,0_isp,comm,errcode)
  CALL MPI_REDUCE(tmpi(1),tmpimax(1),n,mpidbl,MPI_MAX,0_isp,comm,errcode)
  CALL MPI_REDUCE(ttile(1),ttilemin(1),n,mpidbl,MPI_MIN,0_isp,comm,errcode)
  CALL MPI_REDUCE(tmpi(1),tmpimin(1),n,mpidbl,MPI_MIN,0_isp,comm,errcode)

  IF (rank.eq.0) THEN

    ttileave(:) = ttileave(:) / (nproc)
    tmpiave(:) = tmpiave(:) / (nproc)

    write(0,*)
    write(0,'(A40)') 'tile and mpi particle communications'
    write(0,'(" _________________________________________________________________________")')
    DO i = 1,n
      write(0,*)
      write(0,'(X,A60)') name(i)
      write(0,'(X,7(A13))')   "min time tile (s)", " ave time tile (s)", " max time tile (s)"
      write(0,'(6(X,E12.5))') ttilemin(i),ttileave(i),ttilemax(i)
      write(0,'(X,7(A18))')   "min time mpi (s)", " ave time mpi (s)", " max time mpi (s)"
      write(0,'(6(X,E12.5))') tmpimin(i),tmpiave(i),tmpimax(i)
      write(0,'(7(A13))')  "sum(x)", "sum(y)", "sum(z)", "err x", "err y", "err z"
      write(0,'(7(X,E12.5))') sumx(i), sumy(i), sumz(i), errx(i), erry(i), errz(i)
      write(0,'(7(A13))') "sum(px)", "sum(py)", "sum(pz)", "err px", "err py", "err pz"
      write(0,'(7(X,E12.5))') sumpx(i), sumpy(i), sumpz(i), errpx(i), errpy(i), errpz(i)
      write(0,'(2(A13))') "sum(ga)", "err gamma"
      write(0,'(2(X,E12.5))') sumga(i), errga(i)
    ENDDO

  ENDIF

END SUBROUTINE

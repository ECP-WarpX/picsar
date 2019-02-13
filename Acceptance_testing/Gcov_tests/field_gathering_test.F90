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
! FIELD_GATHERING_3D_TEST.F90
!
! Test code for the field gathering in 3D
!
! Mathieu Lobet, 2016.08
! ________________________________________________________________________________________

PROGRAM field_gathering_3d_test

  USE mpi
  USE constants
  IMPLICIT NONE

  ! __________________________________________________________________
  ! Parameter declaration

  INTEGER(idp)                             :: i,n
  INTEGER(idp)                             :: np
  INTEGER(idp)                             :: lvect
  INTEGER(idp)                             :: nx,ny,nz
  INTEGER(idp)                             :: nxguard,nyguard,nzguard
  INTEGER                                  :: ierr,provided
  LOGICAL(lp)                              :: l_lower_order_in_v, l_nodal
  LOGICAL(lp)                              :: passed
  REAL(num)                                :: xmin,ymin,zmin
  REAL(num)                                :: xmax,ymax,zmax
  REAL(num)                                :: q,m
  REAL(num)                                :: Ef,Bf
  REAL(num)                                :: epsilon
  REAL(num)                                :: t0
  REAL(num)                                :: dx,dy,dz,dt
  REAL(num), dimension(:), allocatable     :: xp,yp,zp
  REAL(num), dimension(:), allocatable     :: uxp,uyp,uzp
  REAL(num), dimension(:), allocatable     :: gaminv,w
  REAL(num), dimension(:), allocatable     :: ex,ey,ez
  REAL(num), dimension(:), allocatable     :: bx,by,bz
  REAL(num), dimension(:,:,:), allocatable :: exg,eyg,ezg
  REAL(num), dimension(:,:,:), allocatable :: bxg,byg,bzg
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz
  REAL(num), dimension(10)                 :: te,tb
  CHARACTER(len=64), dimension(10)         :: namee,nameb

  INTEGER(idp)                             :: nguard(3), nvalid(3)

  write(0,'(" ____________________________________________________________________________")')
  write(0,*) 'TEST: field gathering 3D'

  ! __ MPI init for the time _____________________________________________
  CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,provided,ierr)

  ! _________________________________________________________________
  ! Parameter initialization

  !$OMP PARALLEL &
  !$OMP DEFAULT(private)

  np = 100

  nx = 100
  ny = 100
  nz = 100

  lvect = 512

  nxguard = 3
  nyguard = 3
  nzguard = 3

  xmin = 0.
  ymin = 0.
  zmin = 0.

  Ef = 1E12 !V/m
  Bf = 1E4  !T

  q = -1._num
  m = 1._num

  dx = 1.E-6
  dy = 1.E-6
  dz = 1.E-6

  dt = 0.5_num * 1._num / sqrt(1._num / dx**2 + 1._num / dy**2 + 1._num / dz**2 )

  epsilon = 1E-6

  l_lower_order_in_v = .TRUE.
  l_nodal = .FALSE.
  passed = .TRUE.

  xmax = xmin + (nx-1)*dx
  ymax = ymin + (ny-1)*dy
  zmax = zmin + (nz-1)*dz

  write(0,*) 'l_lower_order_in_v:',l_lower_order_in_v
  write(0,'(" xmax:",F12.5," ymax:",F12.5," zmax",F12.5 )') xmax,ymax,zmax

  ! __ array allocations ________________________________________________

  ALLOCATE(exg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(eyg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(ezg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  ALLOCATE(bxg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(byg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(bzg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  ALLOCATE(xp(np),yp(np),zp(np))

  ALLOCATE(uxp(np),uyp(np),uzp(np))
  ALLOCATE(w(np),gaminv(np))

  ALLOCATE(ex(np),ey(np),ez(np))
  ALLOCATE(bx(np),by(np),bz(np))

  ! __ Initialization of the arrays _________________________________

  call random_seed()
  !CALL init_random_seed()


  CALL RANDOM_NUMBER(xp(1:np))
  CALL RANDOM_NUMBER(yp(1:np))
  CALL RANDOM_NUMBER(zp(1:np))

  CALL RANDOM_NUMBER(exg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  CALL RANDOM_NUMBER(eyg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  CALL RANDOM_NUMBER(ezg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  CALL RANDOM_NUMBER(bxg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  CALL RANDOM_NUMBER(byg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  CALL RANDOM_NUMBER(bzg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  CALL RANDOM_NUMBER(uxp(1:np))
  CALL RANDOM_NUMBER(uyp(1:np))
  CALL RANDOM_NUMBER(uzp(1:np))

  CALL RANDOM_NUMBER(w(1:np))

  exg = exg*Ef
  eyg = eyg*Ef
  ezg = ezg*Ef

  bxg = bxg*Bf
  byg = byg*Bf
  bzg = bzg*Bf

  xp = xmin + xp*(xmax - xmin)
  yp = ymin + yp*(ymax - ymin)
  zp = zmin + zp*(zmax - zmin)

  gaminv = 1._num / sqrt(1 - (uxp**2 + uyp**2 + uxp**2))

  uxp = uxp * gaminv
  uyp = uyp * gaminv
  uzp = uzp * gaminv

  sumex = 0
  sumey = 0
  sumez = 0
  sumbx = 0
  sumby = 0
  sumbz = 0
  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0

  ! ___________________________________________________
  ! Test of the functions

  ! __ Electric and Magnetic field Order 1 ________________________________________

  i = 1
  !write(0,*) 'test reference: pxr_gete3d_n_energy_conserving'
  namee(i) = 'pxr_gete3d_n_energy_conserving'
  nameb(i) = 'pxr_getb3d_n_energy_conserving'
  t0 = MPI_WTIME()
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL pxr_gete3d_n_energy_conserving(np,xp,yp,zp, &
      ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,1_idp,1_idp,1_idp, &
      exg,nguard,nvalid, &
      eyg,nguard,nvalid, &
      ezg,nguard,nvalid, &
      .FALSE._idp,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() - t0

  t0 = MPI_WTIME()
  CALL pxr_getb3d_n_energy_conserving(np,xp,yp,zp, &
      bx,by,bz,xmin,ymin,zmin,dx,dy,dz,1_idp,1_idp,1_idp, &
      bxg,nguard,nvalid, &
      byg,nguard,nvalid, &
      bzg,nguard,nvalid, &
      .FALSE._idp,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME()-t0

  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) ;
  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) ;
  !write(0,*) sum(ex),sum(ey),sum(ez)

  !write(0,*) 'test 2: gete3d_energy_conserving_scalar_1_1_1'
  i = i + 1
  namee(i) = 'gete3d_energy_conserving_scalar_1_1_1'
  nameb(i) = 'getb3d_energy_conserving_scalar_1_1_1'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  t0 = MPI_WTIME()
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL gete3d_energy_conserving_scalar_1_1_1(np,xp,yp,zp, &
    ex,ey,ez,xmin,ymin,zmin,dx,dy,dz, &
    exg,nguard,nvalid, &
    eyg,nguard,nvalid, &
    ezg,nguard,nvalid, &
    l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() - t0
  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_scalar_1_1_1(np,xp,yp,zp, &
    bx,by,bz,xmin,ymin,zmin,dx,dy,dz, &
    bxg,nguard,nvalid, &
    byg,nguard,nvalid, &
    bzg,nguard,nvalid, &
    l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV2_1_1_1'
  nameb(i) = 'geteb3d_energy_conserving_vecV2_1_1_1'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV2_1_1_1(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV3_1_1_1'
  nameb(i) = 'geteb3d_energy_conserving_vecV3_1_1_1'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV3_1_1_1(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0
  tb(i) = 0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV4_1_1_1'
  nameb(i) = 'geteb3d_energy_conserving_vecV4_1_1_1'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV4_1_1_1(np,xp,yp,zp, &
         ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz, &
         exg,nguard,nvalid, &
	 eyg,nguard,nvalid, &
	 ezg,nguard,nvalid, &
	 bxg,nguard,nvalid, &
	 byg,nguard,nvalid, &
	 bzg,nguard,nvalid, &
	 lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0
  tb(i) = 0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  ! _________________________________________
  ! Test of extra developer's functions
#if defined(DEV)

  i = i + 1
  namee(i) = 'gete3d_energy_conserving_vec_1_1_1'
  nameb(i) = 'getb3d_energy_conserving_vec_1_1_1'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  t0 = MPI_WTIME()
  CALL gete3d_energy_conserving_vec_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_vec_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV1_1_1_1'
  nameb(i) = 'geteb3d_energy_conserving_vecV1_1_1_1'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV1_1_1_1(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

#endif

  n = i

  write(0,*)
  write(0,'(" Results Electric and Magnetic field gathering order 1")')
  write(0,'(A40, 7(A13))') "Subroutines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez", "time (s)"
  write(0,'(A40, 7(A13))') "Subroutines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,7(X,E12.5))') namee(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i), te(i)
    write(0,'(A40,7(X,E12.5))') nameb(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i), tb(i)
    write(0,'(X,"Total time: ",(X,E12.5))') tb(i) + te(i)
    write(0,*)
  ENDDO

  ! ______________________________________________________________________________________
  ! Electric and Magnetic field Order 2
  sumex = 0 ; sumey = 0 ; sumez = 0
  sumbx = 0 ; sumby = 0 ; sumbz = 0
  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  bx = 0
  by = 0
  bz = 0
  ex = 0
  ey = 0
  ez = 0
  te = 0
  tb = 0
  write(0,*)


  i = 1
  !write(0,*) 'test reference: pxr_gete3d_n_energy_conserving'
  namee(i) = 'pxr_gete3d_n_energy_conserving'
  nameb(i) = 'pxr_getb3d_n_energy_conserving'
  t0 = MPI_WTIME()
  CALL pxr_gete3d_n_energy_conserving(np,xp,yp,zp, &
      ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,2_idp,2_idp,2_idp, &
      exg,nguard,nvalid, &
      eyg,nguard,nvalid, &
      ezg,nguard,nvalid, &
      .FALSE._idp,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL pxr_getb3d_n_energy_conserving(np,xp,yp,zp, &
      bx,by,bz,xmin,ymin,zmin,dx,dy,dz,2_idp,2_idp,2_idp, &
      bxg,nguard,nvalid, &
      byg,nguard,nvalid, &
      bzg,nguard,nvalid, &
      .FALSE._idp,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) ;
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'gete3d_energy_conserving_scalar_2_2_2'
  nameb(i) = 'getb3d_energy_conserving_scalar_2_2_2'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL gete3d_energy_conserving_scalar_2_2_2(np,xp,yp,zp, &
    ex,ey,ez,xmin,ymin,zmin,dx,dy,dz, &
    exg,nguard,nvalid, &
    eyg,nguard,nvalid, &
    ezg,nguard,nvalid, &
    l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0
  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_scalar_2_2_2(np,xp,yp,zp, &
    bx,by,bz,xmin,ymin,zmin,dx,dy,dz, &
    bxg,nguard,nvalid, &
    byg,nguard,nvalid, &
    bzg,nguard,nvalid, &
    l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV3_2_2_2'
  nameb(i) = 'geteb3d_energy_conserving_vecV3_2_2_2'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV3_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV4_2_2_2'
  nameb(i) = 'geteb3d_energy_conserving_vecV4_2_2_2'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL geteb3d_energy_conserving_vecV4_2_2_2(np,xp,yp,zp, &
     ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz, &
     exg,nguard,nvalid, &
	 eyg,nguard,nvalid, &
	 ezg,nguard,nvalid, &
	 bxg,nguard,nvalid, &
	 byg,nguard,nvalid, &
	 bzg,nguard,nvalid, &
	 lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  ! _________________________________________
  ! Test of extra developer's functions
#if defined(DEV)

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV1_2_2_2'
  nameb(i) = 'geteb3d_energy_conserving_vecV1_2_2_2'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV1_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV2_2_2_2'
  nameb(i) = 'geteb3d_energy_conserving_vecV2_2_2_2'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV2_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'gete3d_energy_conserving_vec_2_2_2'
  nameb(i) = 'getb3d_energy_conserving_vec_2_2_2'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  t0 = MPI_WTIME()
  CALL gete3d_energy_conserving_vec_2_2_2(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0
  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_vec_2_2_2(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

#endif

  n=i

  write(0,*)
  write(0,'(" Results Electric and Magnetic field gathering order 1")')
  write(0,'(A40, 7(A13))') "Subroutines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez", "time (s)"
  write(0,'(A40, 7(A13))') "Subroutines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,7(X,E12.5))') namee(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i), te(i)
    write(0,'(A40,7(X,E12.5))') nameb(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i), tb(i)
    write(0,'(X,"Total time: ",(X,E12.5))') tb(i) + te(i)
    write(0,*)
  ENDDO


  ! __ Electric and Magnetic field Order 3 ______________________________________

  sumex = 0
  sumey = 0
  sumez = 0
  sumbx = 0
  sumby = 0
  sumbz = 0
  errex = 0
  errey = 0
  errez = 0
  errbx = 0
  errby = 0
  errbz = 0
  bx = 0
  by = 0
  bz = 0
  ex = 0
  ey = 0
  ez = 0
  lvect = 64
  write(0,*)
  i = 1

  !write(0,*) 'test reference: pxr_gete3d_n_energy_conserving'
  namee(i) = 'pxr_gete3d_n_energy_conserving'
  nameb(i) = 'pxr_getb3d_n_energy_conserving'
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  t0 = MPI_WTIME()
  CALL pxr_gete3d_n_energy_conserving(np,xp,yp,zp, &
      ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,3_idp,3_idp,3_idp, &
      exg,nguard,nvalid, &
      eyg,nguard,nvalid, &
      ezg,nguard,nvalid, &
      .FALSE._idp,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  CALL pxr_getb3d_n_energy_conserving(np,xp,yp,zp, &
      bx,by,bz,xmin,ymin,zmin,dx,dy,dz,3_idp,3_idp,3_idp, &
      bxg,nguard,nvalid, &
      byg,nguard,nvalid, &
      bzg,nguard,nvalid, &
      .FALSE._idp,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)

  !write(0,*) 'test 2: gete3d_energy_conserving_scalar_1_1_1'
  i = i + 1
  namee(i) = 'gete3d_energy_conserving_scalar_3_3_3'
  nameb(i) = 'getb3d_energy_conserving_scalar_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  te(i) = MPI_WTIME() -t0
  CALL gete3d_energy_conserving_scalar_3_3_3(np,xp,yp,zp, &
    ex,ey,ez,xmin,ymin,zmin,dx,dy,dz, &
    exg,nguard,nvalid, &
    eyg,nguard,nvalid, &
    ezg,nguard,nvalid, &
    l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_scalar_3_3_3(np,xp,yp,zp, &
    bx,by,bz,xmin,ymin,zmin,dx,dy,dz, &
    bxg,nguard,nvalid, &
    byg,nguard,nvalid, &
    bzg,nguard,nvalid, &
    l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vecV2_3_3_3'
  nameb(i) = 'geteb3d_energy_conserving_vecV2_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vecV2_3_3_3(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  ! _________________________________________
  ! Test of extra developer's functions
#if defined(DEV)

  !write(0,*) 'test 2: gete3d_energy_conserving_scalar_1_1_1'
  i = i + 1
  namee(i) = 'gete3d_energy_conserving_linear_3_3_3'
  nameb(i) = 'getb3d_energy_conserving_linear_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0
  te(i) = MPI_WTIME() -t0
  CALL gete3d_energy_conserving_linear_3_3_3(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_linear_3_3_3(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'gete3d_energy_conserving_vec_3_3_3'
  nameb(i) = 'getb3d_energy_conserving_vec_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  CALL gete3d_energy_conserving_vec_3_3_3(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_vec_3_3_3(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_vec_3_3_3'
  nameb(i) = 'geteb3d_energy_conserving_vec_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_vec_3_3_3(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'gete3d_energy_conserving_vec2_3_3_3'
  nameb(i) = 'getb3d_energy_conserving_vec2_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  CALL gete3d_energy_conserving_vec2_3_3_3(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0

  t0 = MPI_WTIME()
  CALL getb3d_energy_conserving_vec2_3_3_3(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  tb(i) = MPI_WTIME() -t0

  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

  i = i + 1
  namee(i) = 'geteb3d_energy_conserving_blockvec2_3_3_3'
  nameb(i) = 'getwb3d_energy_conserving_blockvec2_3_3_3'
  ex = 0 ; ey = 0 ; ez = 0
  bx = 0 ; by = 0 ; bz = 0

  t0 = MPI_WTIME()
  CALL geteb3d_energy_conserving_blockvec2_3_3_3(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v, l_nodal)
  te(i) = MPI_WTIME() -t0


  sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
  sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz)
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

#endif

  ! End test of extra developer's functions
  ! _________________________________________

  n = i

  write(0,*)
  write(0,'(" Results Electric and Magnetic field gathering order 1")')
  write(0,'(A40, 7(A13))') "Subroutines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez", "time (s)"
  write(0,'(A40, 7(A13))') "Subroutines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz", "time (s)"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,7(X,E12.5))') namee(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i), te(i)
    write(0,'(A40,7(X,E12.5))') nameb(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i), tb(i)
    write(0,'(X,"Total time: ",(X,E12.5))') tb(i) + te(i)
    write(0,*)
  ENDDO

  ! __ Check results ___________

  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST FIELD GATHERING 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST FIELD GATHERING 3D FAILED **********  \e[0m \n"')
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')

  !$OMP END PARALLEL

  ! finalyze MPI
  call MPI_FINALIZE(ierr)

END PROGRAM

! ________________________________________________________________________________________
! External functions

! ________________________________________________________________________________________
! subroutine init_random_seed()
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
! ________________________________________________________________________________________
!   use iso_fortran_env, only: int64
!   implicit none
!   integer, allocatable :: seed(:)
!   integer :: i, n, un, istat, dt(8), pid
!   integer(int64) :: t
!
!   call random_seed(size = n)
!   allocate(seed(n))
  ! First try if the OS provides a random number generator
!   open(newunit=un, file="/dev/urandom", access="stream", &
!        form="unformatted", action="read", status="old", iostat=istat)
!   if (istat == 0) then
!      read(un) seed
!      close(un)
!   else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
!      call system_clock(t)
!      if (t == 0) then
!         call date_and_time(values=dt)
!         t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
!              + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
!              + dt(3) * 24_int64 * 60 * 60 * 1000 &
!              + dt(5) * 60 * 60 * 1000 &
!              + dt(6) * 60 * 1000 + dt(7) * 1000 &
!              + dt(8)
!      end if
!      pid = getpid()
!      t = ieor(t, int(pid, kind(t)))
!      do i = 1, n
!         seed(i) = lcg(t)
!      end do
!   end if
!   call random_seed(put=seed)
! contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
!   function lcg(s)
!     integer :: lcg
!     integer(int64) :: s
!     if (s == 0) then
!        s = 104729
!     else
!        s = mod(s, 4294967296_int64)
!     end if
!     s = mod(s * 279470273_int64, 4294967291_int64)
!     lcg = int(mod(s, int(huge(0), int64)), kind(0))
!   end function lcg
! end subroutine init_random_seed

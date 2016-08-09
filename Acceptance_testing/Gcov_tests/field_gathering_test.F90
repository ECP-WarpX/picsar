! ________________________________________________________________________________________
!
! FIELD_GATHERING_3D_TEST.F90
! Test code for the field gathering in 3D
!
! This program does not need MPI
! ________________________________________________________________________________________

PROGRAM field_gathering_3d_test

  USE constants
  IMPLICIT NONE

  INTEGER(idp)                             :: i,n
  INTEGER(idp)                             :: np
  INTEGER(idp)                             :: lvect  
  INTEGER(idp)                             :: nx,ny,nz
  INTEGER(idp)                             :: nxguard,nyguard,nzguard
  LOGICAL(idp)                             :: l_lower_order_in_v
  LOGICAL                                  :: passed  
  REAL(num)                                :: xmin,ymin,zmin
  REAL(num)                                :: xmax,ymax,zmax
  REAL(num)                                :: Ef,Bf   
  REAL(num)                                :: epsilon   
  REAL(num)                                :: dx,dy,dz  
  REAL(num), dimension(:), allocatable     :: xp,yp,zp 
  REAL(num), dimension(:), allocatable     :: ex,ey,ez  
  REAL(num), dimension(:), allocatable     :: bx,by,bz   
  REAL(num), dimension(:,:,:), allocatable :: exg,eyg,ezg 
  REAL(num), dimension(:,:,:), allocatable :: bxg,byg,bzg   
  REAL(num), dimension(10)                 :: sumex,sumey,sumez
  REAL(num), dimension(10)                 :: sumbx,sumby,sumbz  
  REAL(num), dimension(10)                 :: errex,errey,errez
  REAL(num), dimension(10)                 :: errbx,errby,errbz  
  CHARACTER(len=64), dimension(10)         :: name
  CHARACTER(len=512)                       :: line  
    
  write(0,'(" ____________________________________________________________________________")')
  write(0,*) 'TEST: field gathering 3D'

  ! _____________________________________________
  ! Parameter initialization
  
  np = 50000
  
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
  
  dx = 1.E-6
  dy = 1.E-6
  dz = 1.E-6
  
  epsilon = 1E-6
  
  l_lower_order_in_v = .FALSE.
  
  passed = .TRUE.
  
  xmax = xmin + (nx-1)*dx
  ymax = ymin + (ny-1)*dy  
  zmax = zmin + (nz-1)*dz  

  write(0,*) 'l_lower_order_in_v:',l_lower_order_in_v
  write(0,'(" xmax:",F12.5," ymax:",F12.5," zmax",F12.5 )') xmax,ymax,zmax
    
  ALLOCATE(exg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(eyg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))  
  ALLOCATE(ezg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  ALLOCATE(bxg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(byg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))  
  ALLOCATE(bzg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  
  ALLOCATE(xp(np),yp(np),zp(np)) 
    
  ALLOCATE(ex(np),ey(np),ez(np))    
  ALLOCATE(bx(np),by(np),bz(np))   
  ! Initialization of the arrays
  CALL init_random_seed()
  
  CALL RANDOM_NUMBER(xp(1:np))
  CALL RANDOM_NUMBER(yp(1:np))  
  CALL RANDOM_NUMBER(zp(1:np)) 
  
  CALL RANDOM_NUMBER(exg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))   
  CALL RANDOM_NUMBER(eyg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))  
  CALL RANDOM_NUMBER(ezg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))  

  CALL RANDOM_NUMBER(bxg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))   
  CALL RANDOM_NUMBER(byg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))  
  CALL RANDOM_NUMBER(bzg(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))  

  exg = exg*Ef
  eyg = eyg*Ef
  ezg = ezg*Ef

  bxg = bxg*Bf
  byg = byg*Bf
  bzg = bzg*Bf

  xp = xmin + xp*(xmax - xmin)
  yp = ymin + yp*(ymax - ymin)
  zp = zmin + zp*(zmax - zmin)
  
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
  
  ! __ Electric field Order 1 ________________________________________
  
  i = 1
  !write(0,*) 'test reference: pxr_gete3d_n_energy_conserving'
  name(i) = 'pxr_gete3d_n_energy_conserving'
  CALL pxr_gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,&
                                 dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                 1_idp,1_idp,1_idp,exg,eyg,ezg,.FALSE.,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) ; i = i + 1
  !write(0,*) sum(ex),sum(ey),sum(ez)  

  !write(0,*) 'test 1: gete3d_energy_conserving_1_1_1'
  name(i) = 'gete3d_energy_conserving_1_1_1'
  ex = 0
  ey = 0
  ez = 0
  CALL gete3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) 
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1
  !write(0,*) sum(ex),sum(ey),sum(ez)

  !write(0,*) 'test 2: gete3d_energy_conserving_scalar_1_1_1'
  name(i) = 'gete3d_energy_conserving_scalar_1_1_1'
  ex = 0
  ey = 0
  ez = 0
	CALL gete3d_energy_conserving_scalar_1_1_1(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1	

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  name(i) = 'geteb3d_energy_conserving_1_1_1'
  ex = 0
  ey = 0
  ez = 0
	CALL geteb3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1	

	n = i-1
	
	write(0,*)
  write(0,'(" Results Electric field order 1")')	
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
		write(0,'(A40,6(X,E12.5))') name(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i)
	ENDDO

  ! __ Magnetic field Order 1 ________________________________________
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
	write(0,*)
	
  i = 1
  !write(0,*) 'test reference: pxr_getb3d_n_energy_conserving'
  name(i) = 'pxr_getb3d_n_energy_conserving'
  CALL pxr_getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,&
                                 dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                 1_idp,1_idp,1_idp,bxg,byg,bzg,.FALSE.,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) ; i = i + 1

  !write(0,*) 'test 1: getb3d_energy_conserving_1_1_1'
  name(i) = 'getb3d_energy_conserving_1_1_1'
  bx = 0
  by = 0
  bz = 0
  CALL getb3d_energy_conserving_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)
	IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1

  !write(0,*) 'test 2: getb3d_energy_conserving_scalar_1_1_1'
  name(i) = 'getb3d_energy_conserving_scalar_1_1_1'
  bx = 0
  by = 0
  bz = 0
  CALL getb3d_energy_conserving_scalar_1_1_1(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)		
	i = i + 1

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  name(i) = 'geteb3d_energy_conserving_1_1_1'
  bx = 0
  by = 0
  bz = 0
	CALL geteb3d_energy_conserving_1_1_1(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)		
	i = i + 1

	n = i-1
	write(0,*)
  write(0,'(" Results Magnetic field order 1")')	
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
		write(0,'(A40,6(X,E12.5))') name(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i)
	ENDDO
  ! __ Electric field Order 2 ______________________________________
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
	write(0,*)


  i = 1
  !write(0,*) 'test reference: pxr_gete3d_n_energy_conserving'
  name(i) = 'pxr_gete3d_n_energy_conserving'
  CALL pxr_gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,&
                                 dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                 2_idp,2_idp,2_idp,exg,eyg,ezg,.FALSE.,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) ; i = i + 1

  !write(0,*) 'test 1: gete3d_energy_conserving_2_2_2'
  name(i) = 'gete3d_energy_conserving_2_2_2'
  ex = 0
  ey = 0
  ez = 0
  CALL gete3d_energy_conserving_2_2_2(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) 
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1

  !write(0,*) 'test 2: geteb3d_energy_conserving_2_2_2'
  name(i) = 'geteb3d_energy_conserving_2_2_2'
  ex = 0
  ey = 0
  ez = 0
	CALL geteb3d_energy_conserving_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1	

	n = i-1
	write(0,*)
  write(0,'(" Results Electric field order 2")')	
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
		write(0,'(A40,6(X,E12.5))') name(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i)
	ENDDO

  ! __ Magnetic field Order 2 ______________________________________
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
	write(0,*)
	
  i = 1
  !write(0,*) 'test reference: pxr_getb3d_n_energy_conserving'
  name(i) = 'pxr_getb3d_n_energy_conserving'
  CALL pxr_getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,&
                                 dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                 2_idp,2_idp,2_idp,bxg,byg,bzg,.FALSE.,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) ; i = i + 1

  !write(0,*) 'test 1: getb3d_energy_conserving_2_2_2'
  name(i) = 'getb3d_energy_conserving_2_2_2'
  bx = 0
  by = 0
  bz = 0
  CALL getb3d_energy_conserving_2_2_2(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)
	IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1

  !write(0,*) 'test 2: geteb3d_energy_conserving_2_2_2'
  name(i) = 'geteb3d_energy_conserving_2_2_2'
  bx = 0
  by = 0
  bz = 0
  ex = 0
  ey = 0
  ez = 0
	CALL geteb3d_energy_conserving_2_2_2(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)		
	i = i + 1

	n = i-1
	write(0,*)
  write(0,'(" Results Magnetic field order 2")')	
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
		write(0,'(A40,6(X,E12.5))') name(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i)
	ENDDO

  ! __ Electric field Order 3 ______________________________________
  
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
	write(0,*)
  i = 1

  !write(0,*) 'test reference: pxr_gete3d_n_energy_conserving'
  name(i) = 'pxr_gete3d_n_energy_conserving'
  CALL pxr_gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,&
                                 dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                 3_idp,3_idp,3_idp,exg,eyg,ezg,.FALSE.,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) ; i = i + 1

  !write(0,*) 'test 1: gete3d_energy_conserving_2_2_2'
  name(i) = 'gete3d_energy_conserving_3_3_3'
  ex = 0
  ey = 0
  ez = 0
  CALL gete3d_energy_conserving_3_3_3(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez) 
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1

  !write(0,*) 'test 2: gete3d_energy_conserving_scalar_1_1_1'
  name(i) = 'gete3d_energy_conserving_scalar_3_3_3'
  ex = 0
  ey = 0
  ez = 0
	CALL gete3d_energy_conserving_scalar_3_3_3(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1	

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  name(i) = 'geteb3d_energy_conserving_3_3_3'
  ex = 0
  ey = 0
  ez = 0
	CALL geteb3d_energy_conserving_3_3_3(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
	sumex(i)=sum(ex) ; sumey(i) = sum(ey) ; sumez(i) = sum(ez)
	errex(i) = abs((sumex(i) - sumex(1)))/sumex(1)
	errey(i) = abs((sumey(i) - sumey(1)))/sumey(1)
	errez(i) = abs((sumez(i) - sumez(1)))/sumez(1)
	IF (errex(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errey(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errez(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1	

	n = i-1
	write(0,*)
  write(0,'(" Results Electric field order 3")')	
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(ex)", "sum(ey)", "sum(ez)", "err ex", "err ey", "err ez"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
		write(0,'(A40,6(X,E12.5))') name(i), sumex(i), sumey(i), sumez(i), errex(i), errey(i), errez(i)
	ENDDO

  ! __ Magnetic field Order 3 ______________________________________
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
	write(0,*)

  i = 1
  !write(0,*) 'test reference: pxr_getb3d_n_energy_conserving'
  name(i) = 'pxr_getb3d_n_energy_conserving'
  CALL pxr_getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,&
                                 dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                 3_idp,3_idp,3_idp,bxg,byg,bzg,.FALSE.,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) ; i = i + 1

  name(i) = 'getb3d_energy_conserving_3_3_3'
  bx = 0
  by = 0
  bz = 0
  CALL getb3d_energy_conserving_3_3_3(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)
	IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))
	i = i + 1

  !write(0,*) 'test 2: getb3d_energy_conserving_scalar_1_1_1'
  name(i) = 'getb3d_energy_conserving_scalar_3_3_3'
  bx = 0
  by = 0
  bz = 0
  CALL getb3d_energy_conserving_scalar_3_3_3(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      bxg,byg,bzg,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)
	IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))			
	i = i + 1

  !write(0,*) 'test 3: geteb3d_energy_conserving_1_1_1'
  name(i) = 'geteb3d_energy_conserving_3_3_3'
  bx = 0
  by = 0
  bz = 0
	CALL geteb3d_energy_conserving_3_3_3(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,   &
                                      dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      exg,eyg,ezg,bxg,byg,bzg,lvect,l_lower_order_in_v)
	sumbx(i)=sum(bx) ; sumby(i) = sum(by) ; sumbz(i) = sum(bz) 
	errbx(i) = abs((sumbx(i) - sumbx(1)))/sumbx(1)
	errby(i) = abs((sumby(i) - sumby(1)))/sumby(1)
	errbz(i) = abs((sumbz(i) - sumbz(1)))/sumbz(1)
	IF (errbx(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errby(i) .gt. epsilon) passed = (passed.and.(.false.))
	IF (errbz(i) .gt. epsilon) passed = (passed.and.(.false.))	
	i = i + 1

	n = i-1
	write(0,*)
  write(0,'(" Results Magnetic field order 3")')	
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(bx)", "sum(by)", "sum(bz)", "err bx", "err by", "err bz"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    !write(line,'("printf ""\e[32m",A40,6(X,E12.5),"\e[0m \n """)') name(i), sumbx(i), sumby(i), &
    !sumbz(i), errbx(i), errby(i), errbz(i)
    !CALL system(trim(adjustl(line)))
		write(0,'(A40,6(X,E12.5))') name(i), sumbx(i), sumby(i), sumbz(i), errbx(i), errby(i), errbz(i)
	ENDDO

  ! __ Check results ___________
  
	write(0,*)
  IF (passed) THEN
		!write(0,'("\033[32m **** TEST PASSED **** \033[0m")')	
		!CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')  
		CALL system('printf "\e[32m ********** TEST FIELD GATHERING PASSED **********  \e[0m \n"')
  ELSE
		!write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
		!CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'") 		
		CALL system('printf "\e[31m ********** TEST FIELD GATHERING FAILED **********  \e[0m \n"')
  ENDIF
  
  write(0,'(" ____________________________________________________________________________")')

END PROGRAM

! ________________________________________________________________________________________
! External functions

! ________________________________________________________________________________________
subroutine init_random_seed()
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
! ________________________________________________________________________________________
	use iso_fortran_env, only: int64
	implicit none
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid
	integer(int64) :: t

	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(newunit=un, file="/dev/urandom", access="stream", &
			 form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
		 read(un) seed
		 close(un)
	else
		 ! Fallback to XOR:ing the current time and pid. The PID is
		 ! useful in case one launches multiple instances of the same
		 ! program in parallel.
		 call system_clock(t)
		 if (t == 0) then
				call date_and_time(values=dt)
				t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
						 + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
						 + dt(3) * 24_int64 * 60 * 60 * 1000 &
						 + dt(5) * 60 * 60 * 1000 &
						 + dt(6) * 60 * 1000 + dt(7) * 1000 &
						 + dt(8)
		 end if
		 pid = getpid()
		 t = ieor(t, int(pid, kind(t)))
		 do i = 1, n
				seed(i) = lcg(t)
		 end do
	end if
	call random_seed(put=seed)
contains
	! This simple PRNG might not be good enough for real work, but is
	! sufficient for seeding a better PRNG.
	function lcg(s)
		integer :: lcg
		integer(int64) :: s
		if (s == 0) then
			 s = 104729
		else
			 s = mod(s, 4294967296_int64)
		end if
		s = mod(s * 279470273_int64, 4294967291_int64)
		lcg = int(mod(s, int(huge(0), int64)), kind(0))
	end function lcg
end subroutine init_random_seed

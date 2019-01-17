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
! CURRENT_DEPOSITION_3D_TEST.F90
!
! Test code for the classical current deposition in 3D
!
! This program does not need MPI or OpenMP but some subroutines are vectorized
! Mathieu Lobet, 2016.08
! ________________________________________________________________________________________

PROGRAM current_deposition_3d_test

  USE constants

  IMPLICIT NONE

  ! __________________________________________________________________
  ! Parameter declaration

  INTEGER(idp)                             :: i,n
  INTEGER(idp)                             :: np
  INTEGER(idp)                             :: lvect
  INTEGER(idp)                             :: nx,ny,nz
  INTEGER(idp)                             :: ncx,ncy,ncz
  INTEGER(idp)                             :: nxguard,nyguard,nzguard
  INTEGER(idp)                             :: ncells
  LOGICAL(lp)                              :: passed
  REAL(num)                                :: xmin,ymin,zmin
  REAL(num)                                :: xmax,ymax,zmax
  REAL(num)                                :: q,m
  REAL(num)                                :: epsilon
  REAL(num)                                :: dx,dy,dz,dt
  REAL(num), dimension(:), allocatable     :: xp,yp,zp
  REAL(num), dimension(:), allocatable     :: uxp,uyp,uzp
  REAL(num), dimension(:), allocatable     :: up,theta,phi
  REAL(num), dimension(:), allocatable     :: gaminv,w
  REAL(num), dimension(:,:,:), allocatable :: jx,jy,jz
  REAL(num), dimension(:,:), allocatable :: jxcells,jycells,jzcells
  REAL(num), dimension(10)                 :: sumjx,sumjy,sumjz
  REAL(num), dimension(10)                 :: errjx,errjy,errjz
  CHARACTER(len=64), dimension(10)         :: name

  INTEGER(idp)                             :: nguard(3), nvalid(3)

  write(0,'(" ____________________________________________________________________________")')
  write(0,*) 'TEST: current deposition 3D'

  ! _________________________________________________________________
  ! Parameter initialization

  np = 50000

  nx = 100
  ny = 100
  nz = 100

  q = -1.
  m = 1.

  nxguard = 3
  nyguard = 3
  nzguard = 3

  xmin = 0.
  ymin = 0.
  zmin = 0.

  dx = 1.E-6
  dy = 1.E-6
  dz = 1.E-6

  dt = 0.5_num * 1._num / ((clight *sqrt(1._num / dx**2 + 1._num / dy**2 + 1._num / dz**2 )))

  epsilon = 1E-3

  passed = .TRUE.

  xmax = xmin + (nx-1)*dx
  ymax = ymin + (ny-1)*dy
  zmax = zmin + (nz-1)*dz

  lvect = 8

  write(0,'(" xmax:",F12.5," ymax:",F12.5," zmax",F12.5 )') xmax,ymax,zmax
  write(0,'(" dt:",E12.5)') dt


  ! __ array allocations ________________________________________________

  ALLOCATE(jx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(jy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(jz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  ALLOCATE(xp(np),yp(np),zp(np))

  ALLOCATE(up(np),theta(np),phi(np))

  ALLOCATE(uxp(np),uyp(np),uzp(np))
  ALLOCATE(w(np),gaminv(np))

  ! __ Initialization of the arrays _________________________________
  CALL random_seed()
  !CALL init_random_seed()

  CALL RANDOM_NUMBER(xp(1:np))
  CALL RANDOM_NUMBER(yp(1:np))
  CALL RANDOM_NUMBER(zp(1:np))

  CALL RANDOM_NUMBER(up(1:np))
  CALL RANDOM_NUMBER(theta(1:np))
  CALL RANDOM_NUMBER(phi(1:np))

  CALL RANDOM_NUMBER(w(1:np))

  xp = xmin + xp*(xmax - xmin)
  yp = ymin + yp*(ymax - ymin)
  zp = zmin + zp*(zmax - zmin)

  theta = -0.5*pi + theta*pi
  phi = phi*2*pi

  up = up * 200.

  uxp = up*cos(theta)*cos(phi)
  uyp = up*cos(theta)*sin(phi)
  uzp = up*sin(theta)

  gaminv = 1./sqrt(1._num + (uxp**2 + uyp**2 + uzp**2))

  errjx = 0
  errjy = 0
  errjz = 0

  DEALLOCATE(up,theta,phi)

  write(0,'(" Total particle energy:",E12.5)') sum(w/gaminv)
  write(0,'(" Min/Max energy:",2(X,E12.5))') minval(1/gaminv),maxval(1/gaminv)
  write(0,'(" Min velocity:",3(X,E12.5))') minval(uxp*gaminv),minval(uyp*gaminv),minval(uzp*gaminv)
  write(0,'(" Min momentum:",3(X,E12.5))') minval(uxp),minval(uyp),minval(uzp)
  write(0,'(" Max momentum:",3(X,E12.5))') maxval(uxp),maxval(uyp),maxval(uzp)
  ! _____________________________________________________________
  ! Test of the functions

  ! __ Order 1 __________________________________________________

  i = 1
  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'pxr_depose_jxjyjz_esirkepov_n'
  !print*,trim(adjustl(name(i)))
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL pxr_depose_jxjyjz_esirkepov_n( &
          jx,nguard,nvalid,     &
          jy,nguard,nvalid,     &
          jz,nguard,nvalid,     &
          np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
          dt,dx,dy,dz,1_idp,1_idp,1_idp,.TRUE._idp,.FALSE._idp)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.

  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  name(i) = 'depose_jxjyjz_scalar_1_1_1'
  CALL depose_jxjyjz_scalar_1_1_1( &
       jx,nguard,nvalid,     &
       jy,nguard,nvalid,     &
       jz,nguard,nvalid,     &
       np,xp,yp,zp,uxp,uyp,uzp,    &
       gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz)

  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  name(i) = 'depose_jxjyjz_vecHVv2_1_1_1'
  CALL depose_jxjyjz_vecHVv2_1_1_1( &
       jx,nguard,nvalid,     &
       jy,nguard,nvalid,     &
       jz,nguard,nvalid,     &
       np,xp,yp,zp,uxp,uyp,uzp,     &
       gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_vecHV_vnr_1_1_1'
  ! Determine ncells
  ncx=nx+3; ncy=ny+3; ncz=nz+3
  ncells = ncx*ncy*ncz
  ALLOCATE(jxcells(8,ncells),jycells(8,ncells),jzcells(8,ncells))
  jxcells = 0 ; jycells = 0; jzcells = 0
  CALL depose_jxjyjz_vecHV_vnr_1_1_1(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz,lvect)
  CALL current_reduction_1_1_1(jx,jy,jz,jxcells,jycells,jzcells,ncells,nx,ny,nz,&
                               nxguard,nyguard,nzguard,ncx,ncy,ncz)
  WRITE(0,*)'Sum:',sum(jxcells),sum(jycells),sum(jzcells)
  DEALLOCATE(jxcells,jycells,jzcells)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  n = i-1
  write(0,*)
  write(0,'(" Current deposition order 1")')
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(jx)", "sum(jy)", "sum(jz)", "err jx", "err jy", "err jz"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,6(X,E12.5))') name(i), sumjx(i), sumjy(i), sumjz(i), errjx(i), errjy(i), errjz(i)
  ENDDO

  ! __ Order 2 __________________________________________________

  i = 1
  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'pxr_depose_jxjyjz_esirkepov_n'
  !print*,trim(adjustl(name(i)))
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL pxr_depose_jxjyjz_esirkepov_n( &
        jx,nguard,nvalid,     &
        jy,nguard,nvalid,     &
        jz,nguard,nvalid,     &
        np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
        dt,dx,dy,dz,2_idp,2_idp,2_idp,.TRUE._idp,.FALSE._idp)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  i = i + 1


  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_scalar_2_2_2'
  CALL depose_jxjyjz_scalar_2_2_2( &
       jx,nguard,nvalid,     &
       jy,nguard,nvalid,     &
       jz,nguard,nvalid,     &
       np,xp,yp,zp,uxp,uyp,uzp,    &
       gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz) ;
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_vecHVv2_2_2_2'
  CALL depose_jxjyjz_vecHVv2_2_2_2( &
       jx,nguard,nvalid,     &
       jy,nguard,nvalid,     &
       jz,nguard,nvalid,     &
       np,xp,yp,zp,uxp,uyp,uzp,     &
       gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_vecHV_vnr_2_2_2'
  ! Determine ncells
  ! Originally was ncx=nxc+4; ncy=nyc+4; ncz=nzc+4
  ! But we need one more for the algorithm at order 2
  ncx=nx+5; ncy=ny+5; ncz=nz+5
  ncells = ncx*ncy*ncz
  ALLOCATE(jxcells(8,ncells),jycells(8,ncells),jzcells(8,ncells))
  jxcells = 0 ; jycells = 0; jzcells = 0
  CALL depose_jxjyjz_vecHV_vnr_2_2_2(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz,lvect)
  CALL current_reduction_2_2_2(jx,jy,jz,jxcells,jycells,jzcells,ncells,nx,ny,nz,&
                               nxguard,nyguard,nzguard,ncx,ncy,ncz)
  WRITE(0,*)'Sum:',sum(jxcells),sum(jycells),sum(jzcells)
  DEALLOCATE(jxcells,jycells,jzcells)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  n = i-1
  write(0,*)
  write(0,'(" Current deposition order 1")')
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(jx)", "sum(jy)", "sum(jz)", "err jx", "err jy", "err jz"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,6(X,E12.5))') name(i), sumjx(i), sumjy(i), sumjz(i), errjx(i), errjy(i), errjz(i)
  ENDDO

  ! __ Order 3 __________________________________________________

  i = 1
  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'pxr_depose_jxjyjz_esirkepov_n'
  !print*,trim(adjustl(name(i)))
  nguard = (/ nxguard, nyguard, nzguard /)
  nvalid = (/ nx+1, ny+1, nz+1 /)
  CALL pxr_depose_jxjyjz_esirkepov_n( &
          jx,nguard,nvalid,     &
          jy,nguard,nvalid,     &
          jz,nguard,nvalid,     &
          np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
          dt,dx,dy,dz,3_idp,3_idp,3_idp,.TRUE._idp,.FALSE._idp)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_scalar_3_3_3'
  !print*,trim(adjustl(name(i)))
  CALL depose_jxjyjz_scalar_3_3_3( &
       jx,nguard,nvalid,     &
       jy,nguard,nvalid,     &
       jz,nguard,nvalid,     &
       np,xp,yp,zp,uxp,uyp,uzp,    &
       gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_vecHVv3_3_3_3'
  !print*,trim(adjustl(name(i)))
  CALL depose_jxjyjz_vecHVv3_3_3_3( &
       jx,nguard,nvalid,     &
       jy,nguard,nvalid,     &
       jz,nguard,nvalid,     &
       np,xp,yp,zp,uxp,uyp,uzp,     &
       gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  jx(:,:,:) = 0.
  jy(:,:,:) = 0.
  jz(:,:,:) = 0.
  name(i) = 'depose_jxjyjz_vecHV_vnr_3_3_3'
  !print*,trim(adjustl(name(i)))
  ! Determine ncells
  ! Originally was ncx=nxc+4; ncy=nyc+4; ncz=nzc+4
  ! But we need one more for the algorithm at order 2
  ncx=nx+5; ncy=ny+4; ncz=nz+3
  ncells = ncx*ncy*ncz
  ALLOCATE(jxcells(8,ncells),jycells(8,ncells),jzcells(8,ncells))
  jxcells = 0 ; jycells = 0; jzcells = 0
  CALL depose_jxjyjz_vecHV_vnr_3_3_3(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz,lvect)
  CALL current_reduction_3_3_3(jx,jy,jz,jxcells,jycells,jzcells,ncells,nx,ny,nz,&
                               nxguard,nyguard,nzguard,ncx,ncy,ncz)
  WRITE(0,*)'Sum jxcells:',sum(jxcells),'Sum jycells:',sum(jycells),'Sum jzcells',sum(jzcells)
  DEALLOCATE(jxcells,jycells,jzcells)
  sumjx(i)=sum(jx) ; sumjy(i) = sum(jy) ; sumjz(i) = sum(jz)
  errjx(i) = abs((sumjx(i) - sumjx(1)))/abs(sumjx(1))
  errjy(i) = abs((sumjy(i) - sumjy(1)))/abs(sumjy(1))
  errjz(i) = abs((sumjz(i) - sumjz(1)))/abs(sumjz(1))
  IF (errjx(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjy(i) .gt. epsilon) passed = (passed.and.(.false.))
  IF (errjz(i) .gt. epsilon) passed = (passed.and.(.false.))
  i = i + 1

  n = i-1
  write(0,*)
  write(0,'(" Current deposition order 3")')
  write(0,'(A40, 6(A10))') "Subrtouines", "sum(jx)", "sum(jy)", "sum(jz)", "err jx", "err jy", "err jz"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,6(X,E12.5))') name(i), sumjx(i), sumjy(i), sumjz(i), errjx(i), errjy(i), errjz(i)
  ENDDO

  ! __ Check results ___________

  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST CURRENT DEPOSITION 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST CURRENT DEPOSITION 3D FAILED **********  \e[0m \n"')
  ENDIF

  write(0,'(" ____________________________________________________________________________")')


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

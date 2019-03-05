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
! RHO_DEPOSITION_3D_TEST.F90
!
! Test code for the charge deposition in 3D
!
! This program is not parallel but some subroutines are vectorized
! Mathieu Lobet, 2016.08
! ______________________________________________________________________________

PROGRAM rho_deposition_3d_test

  USE constants
  IMPLICIT NONE

  ! __________________________________________________________________
  ! Parameter declaration

  INTEGER(idp)                             :: i,n
  INTEGER(idp)                             :: np,nc
  INTEGER(idp)                             :: lvect
  INTEGER(idp)                             :: nx,ny,nz
  INTEGER(idp)                             :: nxguard,nyguard,nzguard
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
  REAL(num), dimension(:,:,:), allocatable :: rho, rho_ref
  REAL(num), dimension(10)                 :: sum_rho
  REAL(num), dimension(10)                 :: ave_err_rho
  REAL(num), dimension(10)                 :: min_err_rho
  REAL(num), dimension(10)                 :: max_err_rho
  REAL(num), dimension(10)                 :: err_sum_rho
  CHARACTER(len=64), dimension(10)         :: name
  CHARACTER(len=64)                        :: title

  write(0,'(" ____________________________________________________________________________")')
  write(0,*) 'TEST: current deposition 3D'

  ! _________________________________________________________________
  ! Parameter initialization

  ! Number of particles
  np = 50000

  ! Number of cells
  nx = 100
  ny = 100
  nz = 100

  ! Charge
  q = -1.

  ! Mass
  m = 1.

  ! Guard cell number
  nxguard = 3
  nyguard = 3
  nzguard = 3

  ! Origin of the domain
  xmin = 0.
  ymin = 0.
  zmin = 0.

  ! Space steps
  dx = 1.E-6
  dy = 1.E-6
  dz = 1.E-6

  ! Time step computed with the CFL
  dt = 0.5_num * 1._num / ((clight *sqrt(1._num / dx**2 + 1._num / dy**2 + 1._num / dz**2 )))

  ! Max relative error required to pass the tests
  epsilon = 1E-3

  ! Variable that becomes FALSE in case of error
  passed = .TRUE.

  ! Domain max
  xmax = xmin + (nx)*dx
  ymax = ymin + (ny)*dy
  zmax = zmin + (nz)*dz

  ! Vector length for vectorized subroutines
  lvect = 8

  write(0,'(" xmax:",F12.5," ymax:",F12.5," zmax",F12.5 )') xmax,ymax,zmax
  write(0,'(" dt:",E12.5)') dt


  ! __ array allocations ________________________________________________

  nc = (nx+1+2*nxguard)*(ny+1+2*nyguard)*(nz+1+2*nzguard)

  ALLOCATE(rho(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))
  ALLOCATE(rho_ref(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard))

  ALLOCATE(xp(np),yp(np),zp(np))

  ALLOCATE(up(np),theta(np),phi(np))

  ALLOCATE(uxp(np),uyp(np),uzp(np))
  ALLOCATE(w(np),gaminv(np))

  ! __ Initialization of the arrays _________________________________
  call random_seed
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

  ave_err_rho = 0
  err_sum_rho = 0

  DEALLOCATE(up,theta,phi)

  write(0,'(" Total particle energy:",E12.5)') sum(w/gaminv)
  write(0,'(" Min/Max energy:",2(X,E12.5))') minval(1/gaminv),maxval(1/gaminv)
  write(0,'(" Min velocity:",3(X,E12.5))') minval(uxp*gaminv),minval(uyp*gaminv),minval(uzp*gaminv)
  write(0,'(" Min momentum:",3(X,E12.5))') minval(uxp),minval(uyp),minval(uzp)
  write(0,'(" Max momentum:",3(X,E12.5))') maxval(uxp),maxval(uyp),maxval(uzp)
  ! _____________________________________________________________
  ! Test of the functions

  ! __ Order 1 __________________________________________________

  ! Reference
  i = 1
  rho(:,:,:) = 0.
  name(i) = 'pxr_depose_rho_n'
  !print*,trim(adjustl(name(i)))
  CALL pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                    nxguard,nyguard,nzguard,1_idp,1_idp,1_idp, &
                    .TRUE.,.FALSE.)
  sum_rho(i)=sum(rho)
  rho_ref = rho


  i = i + 1
  rho(:,:,:) = 0.
  name(i) = 'depose_rho_scalar_1_1_1'
  !print*,trim(adjustl(name(i)))
  CALL depose_rho_scalar_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
          nxguard,nyguard,nzguard,lvect)
  CALL check(rho_ref,rho,nc,sum_rho(1),sum_rho(i),min_err_rho(i),ave_err_rho(i),&
       max_err_rho(i),err_sum_rho(i),epsilon, passed)

  i = i + 1
  rho(:,:,:) = 0.
  name(i) = 'depose_rho_vecHVv2_1_1_1'
  !print*,trim(adjustl(name(i)))
  CALL depose_rho_vecHVv2_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
          nxguard,nyguard,nzguard,lvect)
  CALL check(rho_ref,rho,nc,sum_rho(1),sum_rho(i),min_err_rho(i),ave_err_rho(i),&
       max_err_rho(i),err_sum_rho(i),epsilon, passed)

  n = i
  title = "Rho deposition order 1"
  CALL output_check(title,INT(n,isp),name,sum_rho,ave_err_rho,err_sum_rho)

  ! __ Order 2 __________________________________________________

  sum_rho = 0
  ave_err_rho = 0
  err_sum_rho = 0

  i = 1
  rho(:,:,:) = 0.
  name(i) = 'pxr_depose_rho_n'
  !print*,trim(adjustl(name(i)))
  CALL pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                    nxguard,nyguard,nzguard,2_idp,2_idp,2_idp, &
                    .TRUE.,.FALSE.)
  sum_rho(i)=sum(rho)
  rho_ref = rho

  i = i + 1
  rho(:,:,:) = 0.
  name(i) = 'depose_rho_scalar_2_2_2'
  !print*,trim(adjustl(name(i)))
  CALL depose_rho_scalar_2_2_2(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
          nxguard,nyguard,nzguard,lvect)
  CALL check(rho_ref,rho,nc,sum_rho(1),sum_rho(i),min_err_rho(i),ave_err_rho(i),&
       max_err_rho(i),err_sum_rho(i),epsilon, passed)

  i = i + 1
  rho(:,:,:) = 0.
  name(i) = 'depose_rho_vecHVv4_2_2_2'
  !print*,trim(adjustl(name(i)))
  CALL depose_rho_vecHVv2_2_2_2(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
          nxguard,nyguard,nzguard,lvect)
  CALL check(rho_ref,rho,nc,sum_rho(1),sum_rho(i),min_err_rho(i),ave_err_rho(i),&
       max_err_rho(i),err_sum_rho(i),epsilon, passed)


  n = i
  title = "Rho deposition order 2"
  CALL output_check(title,INT(n,isp),name,sum_rho,ave_err_rho,err_sum_rho)

  ! __ Order 3 __________________________________________________

  i = 1
  rho(:,:,:) = 0.
  name(i) = 'pxr_depose_rho_n'
  !print*,trim(adjustl(name(i)))
  CALL pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                    nxguard,nyguard,nzguard,3_idp,3_idp,3_idp, &
                    .TRUE.,.FALSE.)
  sum_rho(i)=sum(rho)
  rho_ref = rho

  i = i + 1
  rho(:,:,:) = 0.
  name(i) = 'depose_rho_scalar_3_3_3'
  !print*,trim(adjustl(name(i)))
  CALL depose_rho_scalar_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
          nxguard,nyguard,nzguard,lvect)
  CALL check(rho_ref,rho,nc,sum_rho(1),sum_rho(i),min_err_rho(i),ave_err_rho(i),&
       max_err_rho(i),err_sum_rho(i),epsilon, passed)

  i = i + 1
  rho(:,:,:) = 0.
  name(i) = 'depose_rho_vecHVv4_3_3_3'
  !print*,trim(adjustl(name(i)))
  CALL depose_rho_vecHVv4_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
          nxguard,nyguard,nzguard,lvect)
  CALL check(rho_ref,rho,nc,sum_rho(1),sum_rho(i),min_err_rho(i),ave_err_rho(i),&
       max_err_rho(i),err_sum_rho(i),epsilon, passed)

  title = "Rho deposition order 3"
  CALL output_check(title,INT(i,isp),name,sum_rho,ave_err_rho,err_sum_rho)


  ! __ Check results ___________

  write(0,*)
  IF (passed) THEN
    !write(0,'("\033[32m **** TEST PASSED **** \033[0m")')
    !CALL system('echo -e "\e[32m **** TEST PASSED **** \e[0m"')
    CALL system('printf "\e[32m ********** TEST CHARGE DEPOSITION 3D PASSED **********  \e[0m \n"')
  ELSE
    !write(0,'("\033[31m **** TEST FAILED **** \033[0m")')
    !CALL system("echo -e '\e[31m **********  TEST FAILED ********** \e[0m'")
    CALL system('printf "\e[31m ********** TEST CHARGE DEPOSITION 3D FAILED **********  \e[0m \n"')
    CALL EXIT(9)
  ENDIF

  write(0,'(" ____________________________________________________________________________")')

END PROGRAM

! ________________________________________________________________________________________
! External functions

! ________________________________________________________________
!
! Check the results by comparison with the original cases
!
SUBROUTINE check(rho_ref,rho,nc,sum_rho_ref,sum_rho,min_err_rho,ave_err_rho,max_err_rho,err_sum_rho,epsilon,passed)
! ________________________________________________________________

  USE constants
  IMPLICIT NONE

  INTEGER(idp)                   :: nc
  REAL(num)                      :: ave_err_rho
  REAL(num)                      :: err_sum_rho
  REAL(num)                      :: sum_rho, sum_rho_ref
  REAL(num)                      :: epsilon
  REAL(num), DIMENSION(nc)       :: rho_ref, rho
  LOGICAL(lp), INTENT(INOUT)         :: passed

  REAL(num), DIMENSION(nc)       :: diff_rho
  REAL(num)                      :: max_err_rho,min_err_rho


  diff_rho = rho - rho_ref

  ! Sum
  sum_rho = sum(rho)

  ! Error point by point
  ave_err_rho = abs(sum(diff_rho/rho_ref,MASK=((rho_ref.ne.0).and.(rho.ne.0))))/nc !/abs(sum_rho_ref)
  max_err_rho = maxval(abs(diff_rho/rho_ref),MASK=((rho_ref.ne.0).and.(rho.ne.0)))
  min_err_rho = minval(abs(diff_rho/rho_ref),MASK=((rho_ref.ne.0).and.(rho.ne.0)))

  ! Error sum
  err_sum_rho = abs((sum_rho - sum_rho))/abs(sum_rho)

  !
  IF (err_sum_rho .gt. epsilon) passed = (passed.and.(.false.))
  IF (ave_err_rho .gt. epsilon) passed = (passed.and.(.false.))
  IF (min_err_rho .gt. epsilon) passed = (passed.and.(.false.))
  IF (max_err_rho .gt. epsilon) passed = (passed.and.(.false.))

END SUBROUTINE
! ________________________________________________________________


! ________________________________________________________________
!
! This subroutine output the results of the checks
!
SUBROUTINE output_check(title,n,name,sum_rho,err_rho,err_sum_rho)
! ________________________________________________________________

  USE constants
  IMPLICIT NONE

  CHARACTER(len=64), INTENT(IN)            :: title
  INTEGER(isp), INTENT(IN)                 :: n
  CHARACTER(len=64), dimension(10)         :: name
  REAL(num), dimension(10)                 :: sum_rho
  REAL(num), dimension(10)                 :: err_rho, err_sum_rho
  INTEGER(isp)                             :: i

  write(0,*)
  write(0,'(A40)') title
  write(0,'(A40, 6(X,A13))') "Subroutines", "sum(rho)", "err sum ", "ave err point"
  write(0,'(" _____________________________________________________")')
  DO i = 1,n
    write(0,'(A40,6(X,E13.5))') name(i), sum_rho(i), err_sum_rho(i), err_rho(i)
  ENDDO

END SUBROUTINE


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

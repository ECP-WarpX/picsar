! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016,
! The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a paid-up,
! nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
! copies to the public, prepare derivative works, and perform publicly and display
! publicly, and to permit other to do so.
!
!
! init_Fourier.F90
!
! Purpose:
! This file contains subroutines for Fourier init (allocation and init of Fourier k-vectors, 
! creation of planes for FFTW etc.)
!
! Author:
! Henri Vincenti
!
! Date:
! Creation March, 2017
!
! ________________________________________________________________________________________

! ---  

MODULE fourier_psaotd
IMPLICIT NONE 
CONTAINS 

SUBROUTINE rfftfreq(nxx,kxx,dxx)
USE constants 
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nxx
REAL(num), INTENT(IN) :: dxx
REAL(num), DIMENSION(nxx/2+1), INTENT(OUT) :: kxx
INTEGER(idp) :: i, l,m,n
REAL(num) :: fe

fe=1_num/dxx

! Even case
IF (MOD(nxx,2) .EQ. 0) THEN 
	n=nxx/2
! Odd case 
ELSE
	n=(nxx-1)/2
ENDIF 

kxx(1)=0
DO i=1,n
	kxx(i+1)=kxx(i)+1
END DO

kxx=kxx/(dxx*nxx)
END SUBROUTINE rfftfreq


! ---  fftfreq as in numpy 
SUBROUTINE fftfreq(nxx,kxx, dxx)
USE constants 
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nxx
REAL(num), INTENT(IN) :: dxx
REAL(num), DIMENSION(nxx), INTENT(OUT) :: kxx
INTEGER(idp) :: i,n
REAL(num) :: fe

fe=1_num/dxx

n=nxx

kxx(1)=0
IF (MOD(n,2) .EQ. 0)THEN
  ! First part of k [0,...,n/2-1]
	DO i=1,n/2_idp-1_idp
		kxx(i+1)=kxx(i)+1_idp 
	END DO 
  ! Second part of k [-n/2,-1]
	kxx(n/2_idp+1)=-n/2_idp
	DO i=n/2_idp+1,n-1
		kxx(i+1)=kxx(i)+1_idp
	END DO 
ELSE
  ! First part of k [0,...,(n-1)/2]
	DO i=1,(n-1_idp)/2_idp
		kxx(i+1)=kxx(i)+1_idp 
	END DO 
  ! Second part of k [-(n-1)/2,-1]
	kxx((n-1_idp)/2_idp+2_idp)=-(n-1_idp)/2_idp
	DO i=(n-1_idp)/2_idp+2_idp,n-1
		kxx(i+1)=kxx(i)+1_idp
	END DO
ENDIF  
kxx=kxx/(dxx*nxx)
END SUBROUTINE fftfreq

! ---  fftfreq as in numpy 
SUBROUTINE fftfreq2(nxx,kxx, dxx)
USE constants 
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nxx
REAL(num), INTENT(IN) :: dxx
REAL(num), DIMENSION(nxx), INTENT(OUT) :: kxx
INTEGER(idp) :: i,n
REAL(num) :: fe

fe=1_num/dxx

n=nxx

kxx(1)=0
DO i=1,n-1
	kxx(i+1)=kxx(i)+1
END DO

kxx=kxx/(dxx*n)
END SUBROUTINE fftfreq2


SUBROUTINE init_fourier
USE shared_data 
USE fastfft 
Use fourier 
USE fftw3_fortran
USE fields
USE omp_lib  
IMPLICIT NONE 
REAL(num) :: xi 
INTEGER(idp) :: i, nopenmp, l, m, n
REAL(num), DIMENSION(:), ALLOCATABLE :: wx,wy,wz
REAL(num), DIMENSION(:), POINTER :: temp  
REAL(num), DIMENSION(:), ALLOCATABLE :: xcoefs,ycoefs,zcoefs
REAL(num), DIMENSION(:), ALLOCATABLE :: kxtemp, kytemp, kztemp
COMPLEX(cpx) :: ii 
INTEGER(idp) :: nfftx,nffty,nfftz 
ALLOCATE(kxtemp(nkx),kytemp(nky),kztemp(nkz))
ALLOCATE(xcoefs(nkx),ycoefs(nky),zcoefs(nkz))

IF (fftw_with_mpi) THEN 
	nfftx=nx_global
	nffty=ny_global
	nfftz=nz_global
ELSE
	nfftx=nx+2*nxguards+1
	nffty=ny+2*nyguards+1
	nfftz=nz+2*nzguards+1
ENDIF 

! Init k-vectors 
CALL rfftfreq(nfftx,kxunit,1.0_num)!2._num/dx*pi*(rfftfreq(nx,1.0_num))
CALL fftfreq(nffty,kyunit, 1.0_num)!2._num/dx*pi*(rfftfreq(ny,1.0_num))
IF (fftw_with_mpi) THEN 
	CALL initkvectors_mpi(nfftx,nffty,nfftz)
ELSE
	CALL initkvectors(nfftx,nffty,nfftz)
ENDIF 

kxunit=2._num/dx*pi*kxunit
kyunit=2._num/dy*pi*kyunit
kzunit=2._num/dz*pi*kzunit

! 
IF (l_staggered) THEN
	xi = 2.0_num
ELSE
	xi = 1.0_num 
ENDIF 
ALLOCATE(wx(norderx/2),wy(nordery/2),wz(norderz/2))
CALL FD_weights_hvincenti(norderx,wx,l_staggered)
CALL FD_weights_hvincenti(nordery,wy,l_staggered)
CALL FD_weights_hvincenti(norderz,wz,l_staggered)
IF (rank .EQ. 0) THEN 
	PRINT *, 'COEFFICIENTS HVINCENTI', wx(:)
ENDIF 

! - Init xcoefs 
xcoefs=0._num
DO i=1, norderx/2
	xcoefs=xcoefs+wx(i)*2.0_num*SIN(kxunit*(xi*(i-1)+1_idp)*dx/xi)
ENDDO 
! - Init ycoefs
ycoefs=0._num 
DO i=1, nordery/2
	ycoefs=ycoefs+wy(i)*2.0_num*SIN(kyunit*(xi*(i-1)+1_idp)*dy/xi)
ENDDO 
! - Init zcoefs 
zcoefs=0._num
DO i=1, norderz/2
	zcoefs=zcoefs+wz(i)*2.0_num*SIN(kzunit*(xi*(i-1)+1_idp)*dz/xi)
ENDDO 

! Init kxunit_mod 
kxtemp=kxunit*dx
WHERE(kxtemp==0.) kxtemp=1.0_num!kxtemp(1)=1. ! Put 1 where kx==0
kxunit_mod=kxunit*xcoefs/(kxtemp)
! Init kyunit_mod 
kytemp=kyunit*dy
WHERE(kytemp==0.) kytemp=1.0_num ! Put 1 where ky==0
kyunit_mod=kyunit*ycoefs/(kytemp)
! Init kzunit_mod 
kztemp=kzunit*dz
WHERE(kztemp==0.) kztemp=1.0_num! kztemp(1)=1. ! Put 1 where kz==0
kzunit_mod=kzunit*zcoefs/(kztemp)
! Init kxn, kx_unmod 
DO n=1,nkz
	DO m=1,nky
		kxn(:,m,n) = kxunit_mod
		kx_unmod(:,m,n) = kxunit
	END DO 
END DO 
! Init kyn, ky_unmod 
DO n=1,nkz
	DO l=1,nkx
        kyn(l,:,n) = kyunit_mod
        ky_unmod(l,:,n) = kyunit
	END DO 
END DO 
! Init kzn, kz_unmod 
DO m=1,nky
	DO l=1,nkx
        kzn(l,m,:) = kzunit_mod
        kz_unmod(l,m,:) = kzunit
	END DO 
END DO 
! - Init kx, ky, kz, k , kmag 
kx=kxn
ky=kyn
kz=kzn
k = ABS(SQRT(kxn**2+kyn**2+kzn**2))
kmag = k 
WHERE(kmag==0.) kmag=1._num
!kmag(1,1,1)=1._num ! Remove k=0 
kxn = kxn/kmag
kyn = kyn/kmag
kzn = kzn/kmag

ii=(0,1.)
IF (l_staggered) THEN 
	kxmn = kxn*EXP(-ii*kx_unmod*dx/2.0_num)
	kxpn = kxn*EXP( ii*kx_unmod*dx/2.0_num)
	kymn = kyn*EXP(-ii*ky_unmod*dy/2.0_num)
	kypn = kyn*EXP( ii*ky_unmod*dy/2.0_num)
	kzmn = kzn*EXP(-ii*kz_unmod*dz/2.0_num)
	kzpn = kzn*EXP( ii*kz_unmod*dz/2.0_num)
	kxm = kx*EXP(-ii*kx_unmod*dx/2.0_num)
	kxp = kx*EXP( ii*kx_unmod*dx/2.0_num)
	kym = ky*EXP(-ii*ky_unmod*dy/2.0_num)
	kyp = ky*EXP( ii*ky_unmod*dy/2.0_num)
	kzm = kz*EXP(-ii*kz_unmod*dz/2.0_num)
	kzp = kz*EXP( ii*kz_unmod*dz/2.0_num)
ELSE 
	kxmn = kxn
	kxpn = kxn
	kymn = kyn
	kypn = kyn
	kzmn = kzn
	kzpn = kzn
	kxm = kx
	kxp = kx
	kym = ky
	kyp = ky
	kzm = kz
	kzp = kz
END IF 


#ifdef _OPENMP
      nopenmp=OMP_GET_MAX_THREADS()
#else
      nopenmp=1
#endif

! - Init PSAOTD coefficients 
CALL init_psaotd()
! Compute plans for forward and backward 
! Plan is computed once with one field component and reused with other 
! components 
IF (fftw_with_mpi) THEN 
	CALL init_plans_fourier_mpi(nopenmp) 
ELSE 
	CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp,nfftx,nffty,nfftz,ex, &
		exf,plan_r2c,INT(FFTW_MEASURE,idp),INT(FFTW_FORWARD,idp))
	CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp,nfftx,nffty,nfftz,exf, &
		ex,plan_c2r,INT(FFTW_MEASURE,idp),INT(FFTW_BACKWARD,idp))
ENDIF 

END SUBROUTINE init_fourier 

SUBROUTINE initkvectors_mpi(nfftx,nffty,nfftz)
USE mpi_fftw3
USE PICSAR_precision 
USE fourier 
USE shared_data
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nfftx, nffty, nfftz
REAL(num), DIMENSION(:), ALLOCATABLE :: kxtemp, kytemp, kztemp, kzfftfreq_temp 

CALL rfftfreq(nfftx,kxunit,1.0_num)!2._num/dx*pi*(rfftfreq(nx,1.0_num))
CALL fftfreq(nffty,kyunit, 1.0_num)!2._num/dx*pi*(rfftfreq(ny,1.0_num))
ALLOCATE(kzfftfreq_temp(nz_global+1))
CALL fftfreq(nfftz,kzfftfreq_temp, 1.0_num)!2._num/dx*pi*(rfftfreq(nz,1.0_num))
PRINT *, "rank",  local_z0, local_nz, nz_global+1
kzunit(1:nkz)=kzfftfreq_temp(local_z0+1:local_z0+1+local_nz-1)
DEALLOCATE(kzfftfreq_temp)

END SUBROUTINE initkvectors_mpi 


SUBROUTINE initkvectors(nfftx, nffty, nfftz)
USE PICSAR_precision 
USE fourier 
USE shared_data
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nfftx, nffty, nfftz
CALL rfftfreq(nfftx,kxunit,1.0_num)!2._num/dx*pi*(rfftfreq(nx,1.0_num))
CALL fftfreq(nffty,kyunit, 1.0_num)!2._num/dx*pi*(rfftfreq(ny,1.0_num))
CALL fftfreq(nfftz,kzunit, 1.0_num)!2._num/dx*pi*(rfftfreq(nz,1.0_num))
END SUBROUTINE initkvectors

SUBROUTINE init_plans_fourier_mpi(nopenmp)
USE PICSAR_precision
USE shared_data
USE fields 
USE mpi_fftw3
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nopenmp 
INTEGER(C_INT) :: nopenmp_cint
INTEGER(C_INTPTR_T) :: nx_cint, ny_cint, nz_cint
nopenmp_cint=nopenmp 

IF  (fftw_threads_ok) THEN 
	CALL  DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
ENDIF 
nz_cint=nz_global+1
ny_cint=ny_global+1
nx_cint=nx_global+1

plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint,ny_cint,nx_cint, &
				ex_r, exf, comm, FFTW_MEASURE);
plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(nz_cint,ny_cint,nx_cint, &
				exf, ex_r, comm, FFTW_MEASURE);

END SUBROUTINE init_plans_fourier_mpi 


! - Computes Fourier domain coefficients for an order-p stencil 
SUBROUTINE FD_weights_hvincenti(p,w, is_staggered)
USE picsar_PRECISION  
IMPLICIT NONE 
LOGICAL(idp), INTENT(IN) :: is_staggered 
INTEGER(idp), INTENT(IN) :: p 
REAL(num), DIMENSION(p/2), INTENT(OUT) :: w
INTEGER(idp) :: i, l 
REAL(num) :: lognumer, logdenom 

DO i=1,p/2
	l=i
	IF (is_staggered) THEN 
		lognumer = LOG(16.0_num)*(1.0_num-p/2.0_num)+logfactorial(p-1_idp)*2.0_num
		logdenom = LOG(2.0_num*l-1.0_num)*2.0_num+ &
		logfactorial(p/2_idp+l-1_idp)+logfactorial(p/2_idp-l)+ &
		2.0_num*logfactorial(p/2_idp-1_idp)
	ELSE
		lognumer = logfactorial(p/2_idp)*2.0_num
		logdenom = logfactorial(p/2_idp+l)+ &
		logfactorial(p/2_idp-l)+LOG(1.0_num*l)
	ENDIF 
	w(i) = (-1.0_num)**(l+1)*EXP(lognumer-logdenom)
END DO 
END SUBROUTINE FD_weights_hvincenti

! - Computes factorial of n 
FUNCTION factorial(n)
USE constants 
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: n 
INTEGER(idp) :: factorial 
INTEGER(idp) :: i, Ans 
Ans = 1 
DO i = 1, n
	Ans = Ans * i 
END DO 
factorial = Ans
END FUNCTION factorial

FUNCTION logfactorial(n) ! returns log(n!)
use PICSAR_PRECISION
INTEGER(idp), INTENT(IN)  :: n
REAL(num)                 :: logfactorial,x
INTEGER(idp)              :: k 
IF(n.eq.0) THEN 
    logfactorial=0.
ELSE 
     x=log(1.*n)
     logfactorial=x
      DO k=2,n-1
          x=log(1.*k)
          logfactorial=logfactorial+x
      ENDDO
ENDIF
RETURN 
END FUNCTION logfactorial

SUBROUTINE get_Ffields
USE shared_data
USE fields 
USE fourier 
USE fastfft 
IMPLICIT NONE
INTEGER(idp) :: nfftx,nffty,nfftz 
nfftx=nx+2*nxguards+1
nffty=ny+2*nyguards+1
nfftz=nz+2*nzguards+1

! Get Fourier transform of all fields components and currents 
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ex, exf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ey, eyf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ez, ezf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bx, bxf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,by, byf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bz, bzf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jx, jxf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jy, jyf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jz, jzf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rhoold, rhooldf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rho, rhof, plan_r2c)

END SUBROUTINE get_Ffields


SUBROUTINE get_Ffields_mpi
USE shared_data
USE fields 
USE mpi_fftw3
IMPLICIT NONE
INTEGER(idp) :: ix,iy,iz


! Copy array values before FFT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=1,nz
	DO iy=1,ny_global
		DO ix=1,nx_global
			ex_r(ix,iy,iz)=ex(ix-1,iy-1,iz-1)
			ey_r(ix,iy,iz)=ey(ix-1,iy-1,iz-1)
			ez_r(ix,iy,iz)=ez(ix-1,iy-1,iz-1)
			bx_r(ix,iy,iz)=bx(ix-1,iy-1,iz-1)
			by_r(ix,iy,iz)=by(ix-1,iy-1,iz-1)
			bz_r(ix,iy,iz)=bz(ix-1,iy-1,iz-1)
			jx_r(ix,iy,iz)=jx(ix-1,iy-1,iz-1)
			jy_r(ix,iy,iz)=jy(ix-1,iy-1,iz-1)
			jz_r(ix,iy,iz)=jz(ix-1,iy-1,iz-1)
			rho_r(ix,iy,iz)=rho(ix-1,iy-1,iz-1)
			rhoold_r(ix,iy,iz)=rhoold(ix-1,iy-1,iz-1)
		END DO
	END DO
END DO 
!$OMP END PARALLEL DO 


! Get global Fourier transform of all fields components and currents 

CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ex_r, exf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ey_r, eyf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ez_r, ezf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bx_r, exf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, by_r, eyf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bz_r, ezf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r, jxf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r, jyf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r, jzf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r, rhof)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r, rhooldf)

END SUBROUTINE get_Ffields_mpi


SUBROUTINE get_fields
USE params
USE shared_data
USE fields 
USE fourier 
USE fastfft 
IMPLICIT NONE
REAL(num) :: coeff_norm 
INTEGER(idp) :: ix,iy,iz
! Get Inverse Fourier transform of all fields components and currents 
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,exf, ex, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,eyf, ey, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,ezf, ez, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,bxf, bx, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,byf, by, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,bzf, bz, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,jxf, jx, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,jyf, jy, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,jzf, jz, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,rhof, rho, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nx,ny,nz,rhooldf, rhoold, plan_c2r)
coeff_norm= SIZE(ex)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=-nzguards,nz+nzguards 
	DO iy=-nyguards,ny+nyguards 
		DO ix=-nxguards,nx+nxguards
				ex(ix,iy,iz)=ex(ix,iy,iz)/coeff_norm
				ey(ix,iy,iz)=ey(ix,iy,iz)/coeff_norm 
				ez(ix,iy,iz)=ez(ix,iy,iz)/coeff_norm 
				bx(ix,iy,iz)=bx(ix,iy,iz)/coeff_norm
				by(ix,iy,iz)=by(ix,iy,iz)/coeff_norm
				bz(ix,iy,iz)=bz(ix,iy,iz)/coeff_norm 
				jx(ix,iy,iz)=jx(ix,iy,iz)/coeff_norm 
				jy(ix,iy,iz)=jy(ix,iy,iz)/coeff_norm 
				jz(ix,iy,iz)=jz(ix,iy,iz)/coeff_norm 
				rho(ix,iy,iz)=rho(ix,iy,iz)/coeff_norm 
				rhoold(ix,iy,iz)=rhoold(ix,iy,iz)/coeff_norm
		END DO 
	END DO 
END DO 
!$OMP END PARALLEL DO 
!ex=0.
!ex(1:nkx,:,:)=k*clight*dt
END SUBROUTINE get_fields

SUBROUTINE get_fields_mpi
USE shared_data
USE fields 
USE mpi_fftw3
IMPLICIT NONE
REAL(num) :: coeff_norm 
INTEGER(idp) :: ix,iy,iz

! Get global Fourier transform of all fields components and currents 

CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, exf, ex_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, eyf, ey_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, ezf, ez_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, bxf, ex_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, byf, ey_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, bzf, ez_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jxf, jx_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jyf, jy_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jzf, jz_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, rhof, rho_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, rhooldf, rhoold_r)

coeff_norm=1.!(nx_global)*(ny_global)*(nz)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=1,nz
	DO iy=1,ny_global
		DO ix=1,nx_global
			ex(ix-1,iy-1,iz-1)=ex_r(ix,iy,iz)/coeff_norm
			ey(ix-1,iy-1,iz-1)=ey_r(ix,iy,iz)/coeff_norm
			ez(ix-1,iy-1,iz-1)=ez_r(ix,iy,iz)/coeff_norm
			bx(ix-1,iy-1,iz-1)=bx_r(ix,iy,iz)/coeff_norm
			by(ix-1,iy-1,iz-1)=by_r(ix,iy,iz)/coeff_norm
			bz(ix-1,iy-1,iz-1)=bz_r(ix,iy,iz)/coeff_norm
			jx(ix-1,iy-1,iz-1)=jx_r(ix,iy,iz)/coeff_norm
			jy(ix-1,iy-1,iz-1)=jy_r(ix,iy,iz)/coeff_norm
			jz(ix-1,iy-1,iz-1)=jz_r(ix,iy,iz)/coeff_norm
			rho(ix-1,iy-1,iz-1)=rho_r(ix,iy,iz)/coeff_norm
			rhoold(ix-1,iy-1,iz-1)=rhoold_r(ix,iy,iz)/coeff_norm
		END DO
	END DO
END DO 
!$OMP END PARALLEL DO 

END SUBROUTINE get_fields_mpi




SUBROUTINE init_psaotd 
USE shared_data
USE fourier 
USE params 
IMPLICIT NONE 
COMPLEX(cpx) :: jj
REAL(num) ::  cdt
COMPLEX(cpx), DIMENSION(:,:,:), ALLOCATABLE :: temp 

ALLOCATE(temp(nkx,nky,nkz))
ALLOCATE(EJmult(nkx,nky,nkz))
ALLOCATE(BJmult(nkx,nky,nkz))
ALLOCATE(ERhomult(nkx,nky,nkz))
ALLOCATE(ERhooldmult(nkx,nky,nkz))
ALLOCATE(coswdt(nkx,nky,nkz))
ALLOCATE(sinwdt(nkx,nky,nkz))
ALLOCATE(axm(nkx,nky,nkz),axp(nkx,nky,nkz))
ALLOCATE(aym(nkx,nky,nkz),ayp(nkx,nky,nkz))
ALLOCATE(azm(nkx,nky,nkz),azp(nkx,nky,nkz))


jj=(0.,1.)
cdt=clight*dt 
temp= k*cdt 
coswdt=COS(ABS(temp))
sinwdt=SIN(ABS(temp))

EJmult=-sinwdt/(kmag*clight*eps0)
EJmult(1,1,1)=dt/eps0
ERhomult=jj*(-EJmult/dt-1.0_num/eps0)/kmag 
ERhooldmult = jj*(coswdt/eps0+EJmult/dt)/kmag 
temp=1.0_num/(kmag*clight*eps0) ! Jmult 
BJmult=jj*(coswdt-1.0_num)*temp/clight 

axm = jj*sinwdt*kxmn
axp = jj*sinwdt*kxpn
aym = jj*sinwdt*kymn
ayp = jj*sinwdt*kypn
azm = jj*sinwdt*kzmn
azp = jj*sinwdt*kzpn

END SUBROUTINE init_psaotd 


SUBROUTINE push_psaotd_ebfielfs
USE shared_data
USE fields 
USE fourier 
IMPLICIT NONE 
INTEGER(idp) ::  ix, iy, iz 
REAL(num) :: invc
invc=1.0_num/clight 
! - Push B a full time step 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=1,nkz 
	DO iy=1,nky 
		DO ix=1,nkx
			! - Bx 
			bxf(ix,iy,iz) = coswdt(ix,iy,iz)*bxf(ix,iy,iz) + &
			invc*azp(ix,iy,iz)*eyf(ix,iy,iz) -             &
			invc*ayp(ix,iy,iz)*ezf(ix,iy,iz) +             &
		  kzpn(ix,iy,iz)*BJmult(ix,iy,iz)*jyf(ix,iy,iz) -  &
			kypn(ix,iy,iz)*BJmult(ix,iy,iz)*jzf(ix,iy,iz)

			! - By 
			byf(ix,iy,iz) = coswdt(ix,iy,iz)*byf(ix,iy,iz) - & 
			invc*azp(ix,iy,iz)*exf(ix,iy,iz) 			         + & 
			invc*axp(ix,iy,iz)*ezf(ix,iy,iz) 			         - &
			kzpn(ix,iy,iz)*BJmult(ix,iy,iz)*jxf(ix,iy,iz)  +  &
			kxpn(ix,iy,iz)*BJmult(ix,iy,iz)*jzf(ix,iy,iz) 

			! - Bz 
			bzf(ix,iy,iz) = coswdt(ix,iy,iz)*bzf(ix,iy,iz) + &
			invc*ayp(ix,iy,iz)*exf(ix,iy,iz)   -             &
			invc*axp(ix,iy,iz)*eyf(ix,iy,iz)   +             &
			kypn(ix,iy,iz)*BJmult(ix,iy,iz)*jxf(ix,iy,iz) -  &
			kxpn(ix,iy,iz)*BJmult(ix,iy,iz)*jyf(ix,iy,iz)

			! Push E a full time step 
			! - Ex 
			exf(ix,iy,iz) = coswdt(ix,iy,iz)*exf(ix,iy,iz) -   &
			azm(ix,iy,iz)*clight*byf(ix,iy,iz) +               &
		  aym(ix,iy,iz)*clight*bzf(ix,iy,iz) + 			         &
		  EJmult(ix,iy,iz)*jxf(ix,iy,iz)     +               &
		  kxpn(ix,iy,iz)*ERhomult(ix,iy,iz)*rhof(ix,iy,iz)   &
			+ kxpn(ix,iy,iz)*ERhooldmult(ix,iy,iz)*rhooldf(ix,iy,iz) 

			! - Ey 
			eyf(ix,iy,iz) = coswdt(ix,iy,iz)*eyf(ix,iy,iz) +    &
			azm(ix,iy,iz)*clight*bxf(ix,iy,iz) -                &
			axm(ix,iy,iz)*clight*bzf(ix,iy,iz) +                &
			EJmult(ix,iy,iz)*jyf(ix,iy,iz) +                    &
			kypn(ix,iy,iz)*ERhomult(ix,iy,iz)*rhof(ix,iy,iz)    &
			+ kypn(ix,iy,iz)*ERhooldmult(ix,iy,iz)*rhooldf(ix,iy,iz)

			! - Ez 
			ezf(ix,iy,iz) = coswdt(ix,iy,iz)*ezf(ix,iy,iz) -     &
			aym(ix,iy,iz)*clight*bxf(ix,iy,iz) +				 &
		  axm(ix,iy,iz)*clight*byf(ix,iy,iz) +                 &
			EJmult(ix,iy,iz)*jzf(ix,iy,iz) +                     &
		  kzpn(ix,iy,iz)*ERhomult(ix,iy,iz)*rhof(ix,iy,iz) + &
			kzpn(ix,iy,iz)*ERhooldmult(ix,iy,iz)*rhooldf(ix,iy,iz)
		END DO 
	END DO 
END DO 
!$OMP END PARALLEL DO 

END SUBROUTINE push_psaotd_ebfielfs




! - Push b 
!            mymat.add_op('bx',{'bx':C,'ey': azp/c,'ez':-ayp/c,'jy': kzpn*BJmult,'jz':-kypn*BJmult})
!            mymat.add_op('by',{'by':C,'ex':-azp/c,'ez': axp/c,'jx':-kzpn*BJmult,'jz': kxpn*BJmult})
!            mymat.add_op('bz',{'bz':C,'ex': ayp/c,'ey':-axp/c,'jx': kypn*BJmult,'jy':-kxpn*BJmult})
! - Push e 
!            mymat.add_op('ex',{'ex':C,'by':-azm*c,'bz': aym*c,'jx':EJmult,'rhonew':kxpn*ERhomult,'rhoold':kxpn*ERhooldmult})
!            mymat.add_op('ey',{'ey':C,'bx': azm*c,'bz':-axm*c,'jy':EJmult,'rhonew':kypn*ERhomult,'rhoold':kypn*ERhooldmult})
!            mymat.add_op('ez',{'ez':C,'bx':-aym*c,'by': axm*c,'jz':EJmult,'rhonew':kzpn*ERhomult,'rhoold':kzpn*ERhooldmult})


END MODULE fourier_psaotd
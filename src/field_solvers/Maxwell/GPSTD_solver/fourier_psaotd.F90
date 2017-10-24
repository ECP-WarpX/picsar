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
! fourier_psaotd.F90
!
! Purpose:
! This file contains subroutines for:
! (i)  Fourier init (allocation and init of Fourier k-vectors,
! creation of planes for FFTW etc.),
! (ii) Forward/Backward Fourier transform of EM fields using ,
! FFTW (Distributed/Shared version),
! (iii) Maxwell push in the spectral space using the Pseudo-Spectral Arbitrary Order
! Analytical Time Domain (PSAOTD) Maxwell solver.
!
! Authors:
! Henri Vincenti
! Jean-Luc Vay
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
REAL(num) :: df

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

df=1_num/(dxx*nxx)

kxx=kxx*df

END SUBROUTINE rfftfreq


! ---  fftfreq as in numpy
SUBROUTINE fftfreq(nxx,kxx, dxx)
USE constants
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nxx
REAL(num), INTENT(IN) :: dxx
REAL(num), DIMENSION(nxx), INTENT(OUT) :: kxx
INTEGER(idp) :: i,n
REAL(num) :: df

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

df=1_num/(dxx*nxx)

kxx=kxx*df
END SUBROUTINE fftfreq

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
COMPLEX(cpx) :: ii
INTEGER(idp) :: nfftx,nffty,nfftz, imn,imx,jmn,jmx,kmn,kmx
ALLOCATE(xcoefs(nkx),ycoefs(nky),zcoefs(nkz))

IF (fftw_with_mpi) THEN
	nfftx=nx_global
	nffty=ny_global
	nfftz=nz_global
ELSE
	nfftx=nx+2*nxguards
	nffty=ny+2*nyguards
	nfftz=nz+2*nzguards
ENDIF

! Init k-vectors
IF (fftw_with_mpi) THEN
	CALL initkvectors_mpi(nfftx,nffty,nfftz)
ELSE
	CALL initkvectors(nfftx,nffty,nfftz)
ENDIF

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

! - Init xcoefs
xcoefs=0._num
DO i=1, norderx/2
	xcoefs=xcoefs+wx(i)*2.0_num*SIN(kxunit*dx*((i-1)+1_idp/xi))
ENDDO
! - Init ycoefs
ycoefs=0._num
DO i=1, nordery/2
	ycoefs=ycoefs+wy(i)*2.0_num*SIN(kyunit*dy*((i-1)+1_idp/xi))
ENDDO
! - Init zcoefs
zcoefs=0._num
DO i=1, norderz/2
	zcoefs=zcoefs+wz(i)*2.0_num*SIN(kzunit*dz*((i-1)+1_idp/xi))
ENDDO

! Init kxunit_mod
kxunit_mod=xcoefs/dx
! Init kyunit_mod
kyunit_mod=ycoefs/dy
! Init kzunit_mod
kzunit_mod=zcoefs/dz

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
k =SQRT(kxn**2+kyn**2+kzn**2)
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
	CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp,nfftx,nffty,nfftz,ex_r, &
		exf,plan_r2c,INT(FFTW_MEASURE,idp),INT(FFTW_FORWARD,idp))
	CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp,nfftx,nffty,nfftz,exf, &
		ex_r,plan_c2r,INT(FFTW_MEASURE,idp),INT(FFTW_BACKWARD,idp))
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

CALL rfftfreq(nfftx,kxunit,dx/(2_num*pi))
CALL fftfreq(nffty,kyunit, dy/(2_num*pi))
ALLOCATE(kzfftfreq_temp(nz_global+1))
CALL fftfreq(nfftz,kzfftfreq_temp, dz/(2_num*pi))
kzunit(1:nkz)=kzfftfreq_temp(local_z0+1:local_z0+1+local_nz-1)
DEALLOCATE(kzfftfreq_temp)

END SUBROUTINE initkvectors_mpi


SUBROUTINE initkvectors(nfftx, nffty, nfftz)
USE PICSAR_precision
USE fourier
USE shared_data
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nfftx, nffty, nfftz
CALL rfftfreq(nfftx,kxunit,dx/(2_num*pi))
CALL fftfreq(nffty,kyunit,dy/(2_num*pi))
CALL fftfreq(nfftz,kzunit,dz/(2_num*pi))
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
nz_cint=nz_global
ny_cint=ny_global
nx_cint=nx_global

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
USE time_stat
USE params

IMPLICIT NONE
INTEGER(idp) :: nfftx,nffty,nfftz, nxx,nyy,nzz
REAL(num)    :: tmptime
nfftx=nx+2*nxguards
nffty=ny+2*nyguards
nfftz=nz+2*nzguards
nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;

IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF

! Init fourier fields fields
call normalize_Fourier(ex_r,nfftx,nffty,nfftz,ex,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(ey_r,nfftx,nffty,nfftz,ey,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(ez_r,nfftx,nffty,nfftz,ez,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(bx_r,nfftx,nffty,nfftz,bx,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(by_r,nfftx,nffty,nfftz,by,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(bz_r,nfftx,nffty,nfftz,bz,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(jx_r,nfftx,nffty,nfftz,jx,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(jy_r,nfftx,nffty,nfftz,jy,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(jz_r,nfftx,nffty,nfftz,jz,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(rho_r,nfftx,nffty,nfftz,rho,nxx,nyy,nzz,1.0_num)
call normalize_Fourier(rhoold_r,nfftx,nffty,nfftz,rhoold,nxx,nyy,nzz,1.0_num)
IF (it.ge.timestat_itstart) THEN
  localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
ENDIF

! Do Fourier Transform
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ex_r, exf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ey_r, eyf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ez_r, ezf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bx_r, bxf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,by_r, byf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bz_r, bzf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jx_r, jxf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jy_r, jyf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jz_r, jzf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rhoold_r, rhooldf, plan_r2c)
CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rho_r, rhof, plan_r2c)
IF (it.ge.timestat_itstart) THEN
  localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
ENDIF

END SUBROUTINE get_Ffields


SUBROUTINE get_Ffields_mpi
USE shared_data
USE fields
USE mpi_fftw3
USE time_stat
USE params
IMPLICIT NONE
INTEGER(idp) :: ix,iy,iz
REAL(num)    :: tmptime
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
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
IF (it.ge.timestat_itstart) THEN
  localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
ENDIF


! Get global Fourier transform of all fields components and currents
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF

CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ex_r, exf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ey_r, eyf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ez_r, ezf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bx_r, bxf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, by_r, byf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bz_r, bzf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r, jxf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r, jyf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r, jzf)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r, rhof)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r, rhooldf)
IF (it.ge.timestat_itstart) THEN
  localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
ENDIF

END SUBROUTINE get_Ffields_mpi


SUBROUTINE get_fields
USE params
USE shared_data
USE fields
USE fourier
USE fastfft
USE time_stat
IMPLICIT NONE
REAL(num) :: coeff_norm,tmptime
INTEGER(idp) :: ix,iy,iz,nxx,nyy,nzz,nfftx,nffty,nfftz
nfftx=nx+2*nxguards
nffty=ny+2*nyguards
nfftz=nz+2*nzguards
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF

! Get Inverse Fourier transform of all fields components and currents
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,exf, ex_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,eyf, ey_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,ezf, ez_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,bxf, bx_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,byf, by_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,bzf, bz_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,jxf, jx_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,jyf, jy_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,jzf, jz_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,rhof, rho_r, plan_c2r)
CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,rhooldf, rhoold_r, plan_c2r)

IF (it.ge.timestat_itstart) THEN
  localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
ENDIF

coeff_norm= 1._num/SIZE(ex_r)
nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;

IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
call normalize_Fourier(ex,nxx,nyy,nzz,ex_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(ey,nxx,nyy,nzz,ey_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(ez,nxx,nyy,nzz,ez_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(bx,nxx,nyy,nzz,bx_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(by,nxx,nyy,nzz,by_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(bz,nxx,nyy,nzz,bz_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(jx,nxx,nyy,nzz,jx_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(jy,nxx,nyy,nzz,jy_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(jz,nxx,nyy,nzz,jz_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(rho,nxx,nyy,nzz,rho_r,nfftx,nffty,nfftz,coeff_norm)
call normalize_Fourier(rhoold,nxx,nyy,nzz,rhoold_r,nfftx,nffty,nfftz,coeff_norm)

IF (it.ge.timestat_itstart) THEN
  localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
ENDIF

END SUBROUTINE get_fields

SUBROUTINE normalize_Fourier(ex_out,n1,n2,n3,ex_in,nxx,nyy,nzz,coeff_norm)
USE PICSAR_precision
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz,n1,n2,n3
REAL(num), INTENT(IN) :: coeff_norm
REAL(num), DIMENSION(nxx,nyy,nzz), INTENT(IN OUT) :: ex_in
REAL(num), DIMENSION(n1,n2,n3), INTENT(IN OUT) :: ex_out
INTEGER(idp) :: ix,iy,iz
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=1,MIN(nzz,n3)
	DO iy=1,MIN(nyy,n2)
		DO ix=1,MIN(nxx,n1)
				ex_out(ix,iy,iz)=ex_in(ix,iy,iz)*coeff_norm
		END DO
	END DO
END DO
!$OMP END PARALLEL DO
!ex=0.
END SUBROUTINE normalize_Fourier

SUBROUTINE get_fields_mpi
USE shared_data
USE fields
USE mpi_fftw3
USE time_stat
USE params
IMPLICIT NONE
REAL(num) :: coeff_norm,tmptime
INTEGER(idp) :: ix,iy,iz

! Get global Fourier transform of all fields components and currents
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF

CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, exf, ex_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, eyf, ey_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, ezf, ez_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, exf, ex_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, eyf, ey_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, ezf, ez_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jxf, jx_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jyf, jy_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, jzf, jz_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, rhof, rho_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, rhooldf, rhoold_r)
IF (it.ge.timestat_itstart) THEN
  localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
ENDIF

coeff_norm=(nx_global)*(ny_global)*(nz_global)

IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
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
IF (it.ge.timestat_itstart) THEN
  localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
ENDIF

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
EJmult(1,1,1)=-dt/eps0
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
USE time_stat
USE params
IMPLICIT NONE
INTEGER(idp) ::  ix, iy, iz
REAL(num) :: invc,tmptime
invc=1.0_num/clight
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF

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
IF (it.ge.timestat_itstart) THEN
  localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime)
ENDIF

END SUBROUTINE push_psaotd_ebfielfs

END MODULE fourier_psaotd

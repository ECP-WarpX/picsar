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
INTEGER(idp) :: nfftx,nffty,nfftz, imn,imx,jmn,jmx,kmn,kmx
ALLOCATE(kxtemp(nkx),kytemp(nky),kztemp(nkz))
ALLOCATE(xcoefs(nkx),ycoefs(nky),zcoefs(nkz))

! Init k-vectors 
CALL rfftfreq(nx+2*nxguards,kxunit,1.0_num)!2._num/dx*pi*(rfftfreq(nx,1.0_num))
CALL fftfreq(ny+2*nyguards,kyunit, 1.0_num)!2._num/dx*pi*(rfftfreq(ny,1.0_num))
CALL fftfreq(nz+2*nzguards,kzunit, 1.0_num)!2._num/dx*pi*(rfftfreq(nz,1.0_num))
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
nfftx=nx+2*nxguards
nffty=ny+2*nyguards
nfftz=nz+2*nzguards
CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp,nfftx,nffty,nfftz,ex_r, &
    exf,plan_r2c,INT(FFTW_MEASURE,idp),INT(FFTW_FORWARD,idp))
CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp,nfftx,nffty,nfftz,exf, &
    ex_r,plan_c2r,INT(FFTW_MEASURE,idp),INT(FFTW_BACKWARD,idp))

END SUBROUTINE init_fourier 

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
INTEGER(idp) :: nfftx,nffty,nfftz, nxx,nyy,nzz
nfftx=nx+2*nxguards
nffty=ny+2*nyguards
nfftz=nz+2*nzguards
nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1; 
! Get Fourier transform of all fields components and currents 
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ex(imn:imx,jmn:jmx,kmn:kmx), exf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ey(imn:imx,jmn:jmx,kmn:kmx), eyf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ez(imn:imx,jmn:jmx,kmn:kmx), ezf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bx(imn:imx,jmn:jmx,kmn:kmx), bxf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,by(imn:imx,jmn:jmx,kmn:kmx), byf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bz(imn:imx,jmn:jmx,kmn:kmx), bzf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jx(imn:imx,jmn:jmx,kmn:kmx), jxf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jy(imn:imx,jmn:jmx,kmn:kmx), jyf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jz(imn:imx,jmn:jmx,kmn:kmx), jzf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rhoold(imn:imx,jmn:jmx,kmn:kmx), rhooldf, plan_r2c)
!CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rho(imn:imx,jmn:jmx,kmn:kmx), rhof, plan_r2c)


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

! Do Fourier Transform 
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

END SUBROUTINE get_Ffields

SUBROUTINE get_fields
USE params
USE shared_data
USE fields 
USE fourier 
USE fastfft 
IMPLICIT NONE
REAL(num) :: coeff_norm 
INTEGER(idp) :: ix,iy,iz,nxx,nyy,nzz,nfftx,nffty,nfftz
nfftx=nx+2*nxguards
nffty=ny+2*nyguards
nfftz=nz+2*nzguards
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

!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,exf, ex(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,eyf, ey(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,ezf, ez(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,bxf, bx(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,byf, by(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,bzf, bz(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,jxf, jx(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,jyf, jy(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,jzf, jz(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,rhof, rho(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,rhooldf, rhoold(imn:imx,jmn:jmx,kmn:kmx), plan_c2r)
!coeff_norm= SIZE(ex(imn:imx,jmn:jmx,kmn:kmx))
coeff_norm= 1._num/SIZE(ex_r)
nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1; 
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
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
!DO iz=-nzguards,nz+nzguards-1 
!	DO iy=-nyguards,ny+nyguards-1 
!		DO ix=-nxguards,nx+nxguards-1
!				ex(ix,iy,iz)=ex_r(ix,iy,iz)/coeff_norm
!				ey(ix,iy,iz)=ey_r(ix,iy,iz)/coeff_norm 
!				ez(ix,iy,iz)=ez_r(ix,iy,iz)/coeff_norm 
!				bx(ix,iy,iz)=bx_r(ix,iy,iz)/coeff_norm
!				by(ix,iy,iz)=by_r(ix,iy,iz)/coeff_norm
!				bz(ix,iy,iz)=bz_r(ix,iy,iz)/coeff_norm 
!				jx(ix,iy,iz)=jx_r(ix,iy,iz)/coeff_norm 
!				jy(ix,iy,iz)=jy_r(ix,iy,iz)/coeff_norm 
!				jz(ix,iy,iz)=jz_r(ix,iy,iz)/coeff_norm 
!				rho(ix,iy,iz)=rho_r(ix,iy,iz)/coeff_norm 
!				rhoold(ix,iy,iz)=rhoold_r(ix,iy,iz)/coeff_norm
!		END DO 
!	END DO 
!END DO 
!!$OMP END PARALLEL DO 
!ex=0.
!ex(1:nkx,:,:)=k*clight*dt
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
		DO ix=1,MIN(nyy,n1)
				ex_out(ix,iy,iz)=ex_in(ix,iy,iz)*coeff_norm
		END DO 
	END DO 
END DO 
!$OMP END PARALLEL DO 
!ex=0.
END SUBROUTINE normalize_Fourier

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

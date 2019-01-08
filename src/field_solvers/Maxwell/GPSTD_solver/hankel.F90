
MODULE Matrix_H
  USE PICSAR_precision
  REAL(num), POINTER, DIMENSION(:, :, :) :: Ma_tot
  REAL(num), POINTER, DIMENSION(:, :, :) :: invM_tot


END Matrix_H


MODULE Hankel !# do not parse

IMPLICIT NONE

CONTAINS
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
FUNCTION inv(A) RESULT(Ainv)
  real(num), dimension(:,:), intent(in) :: A
  real(num), dimension(size(A,1),size(A,2)) :: Ainv

  real(num), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer  :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  write(*,*) "DGETRF"
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  write(*,*) "DGETRI"
  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
END FUNCTION inv


! Matrix multiplication for RZ

SUBROUTINE dgemm_example( a,b,c, nfftr)
!!!!!!!!!!!!!!!!
!Brock Palen
!brockp@umich.edu
!
!Example Callding DGEMM from Fortran
!See:
!
!http://www.netlib.org/blas/dgemm.f
!!!!!!!!!!!!!!!!
USE PICSAR_precision
#ifdef _OPENMP
  use omp_lib
#endif
implicit none
EXTERNAL DGEMM
INTEGER (idp), INTENT(in) :: nfftr
INTEGER  (idp):: I, J
REAL(num) ,  INTENT(in) :: a(:,:)
REAL(num) , INTENT(in) :: b(:,:)
REAL(num) , INTENT(out) :: c(:,:)
INTEGER (idp) :: M    ! Number of rows of matrix op(a)
INTEGER (idp) :: N    ! Number of columns of matrix op(b)
INTEGER (idp) :: K    ! Number of columns of matrix op(a) and rows op(b)
INTEGER (idp) :: lda  ! n entry, LDA specifies the first dimension of A as declared
             ! in the calling (sub) program. When  TRANSA = 'N' or 'n' then
             ! LDA must be at least  max( 1, m ), otherwise  LDA must be at
             ! least  max( 1, k ).  
INTEGER (idp) :: ldb  ! On entry, LDB specifies the first dimension of B as declared
             ! in the calling (sub) program. When  TRANSB = 'N' or 'n' then
             ! LDB must be at least  max( 1, k ), otherwise  LDB must be at
             ! least  max( 1, n ). 
INTEGER (idp) :: ldc  ! On entry, LDC specifies the first dimension of C as declared
             ! in  the  calling  (sub)  program.   LDC  must  be  at  least
             ! max( 1, m ).
REAL (num) ::  alpha, beta

#ifdef _OPENMP
 write(*,*)"Will use: ", omp_get_num_procs(), " of ", &
            omp_get_max_threads(), " available"
#endif
M= nfftr
N=M
K=M
lda=M
ldb=M
ldc=M
alpha=1.0
beta=0.0


CALL DGEMM('N','N',M,N,K,alpha,a,lda,b,ldb,beta,c,ldc)

END SUBROUTINE dgemm_example

SUBROUTINE hankel_matrix_init(nfftr,imode,p,rmax,invM)

  USE PICSAR_precision
  !> the matrix initialisation depend on two important variables 
  !> the considered mode imode and the order of hankel transform p
  !> Very important!! : p should be in {m-1, m, m+1}
  INTEGER (idp), intent(in) :: nfftr, imode, p 
  REAL (num), intent(in) :: rmax
  REAL (num), dimension(:,:),allocatable, intent(out):: invM
  REAL (num) :: PI  = 4 * atan (1.0_8)
  !REAL (num) ::  rmax,dr
  !USE shared_data, ONLY: dx, nx
  INTEGER (idp) :: k, j,i, m,p_denom
  REAL (num) , dimension (:), allocatable :: x,r,nu,alphas,bjn,djn,fjn,byn,dyn,fyn,denom
  REAL (num) , dimension (:,:), allocatable :: denom1
  !REAL (num) , dimension (:,:) , allocatable:: alphas1
  REAL (num) , dimension (:,:), allocatable :: y,numer, jn_d, jn_f,yn,yn_d,yn_f
  !REAL (num) , dimension (:,:), allocatable :: invM, Ma
  !REAL (kind =4),dimension (:,:), allocatable :: Ma
  !dr=1.
  !p= 0
  !> For writing simplicity we use m instead of imode 
  m= imode
  !nfftr=10
  !rmax=nx*dx
  ALLOCATE (r(nfftr))
  ALLOCATE (x(nfftr))
  ALLOCATE (nu(nfftr))
  ALLOCATE (alphas(nfftr))
  ALLOCATE (bjn(nfftr))
  ALLOCATE (djn(nfftr))
  ALLOCATE (fjn(nfftr))
  ALLOCATE (byn(nfftr))
  ALLOCATE (dyn(nfftr))
  ALLOCATE (fyn(nfftr))
  ALLOCATE (denom(nfftr))
  !ALLOCATE (r1(nfftr,1))
  ALLOCATE (denom1(nfftr,1))
  !ALLOCATE (alphas1(1,nfftr))
  ALLOCATE (y(nfftr,nfftr))
  ALLOCATE (numer(nfftr,nfftr))
  ALLOCATE (invM(nfftr,nfftr))
  !ALLOCATE (Ma(nfftr,nfftr))
  ALLOCATE (jn_d(nfftr,nfftr))
  ALLOCATE (jn_f(nfftr,nfftr))
  ALLOCATE (yn(nfftr,nfftr))
  ALLOCATE (yn_d(nfftr,nfftr))
  ALLOCATE (yn_f(nfftr,nfftr))
!  ALLOCATE (b(nfftr,nfftr))
!  ALLOCATE (c(nfftr,nfftr))
!  ALLOCATE (d(nfftr,nfftr))
  Do i=1, nfftr
    Do j=1, nfftr
      b(i,j)=1.0
    end do
  end do


  IF (p .eq.m ) THEN
    p_denom= p+1
  ELSE
    p_denom =p
  ENDIF
  CALL  jyzo  (m,nfftr,nu)
  if (m .ne.0) then
    alphas(1)= 0._num
    alphas(2:nfftr)=nu(1:nfftr-1)
  else
    alphas=nu
  endif
 ! Do i =1,nfftr
 !   write(*,*),"alphas",  alphas(i)
 ! end do
  Do i =1, nfftr
    if (alphas(i) .eq. 0._num) then
      bjn(i)=0.
    else
      Call jyndd ( p_denom, alphas(i), bjn(i),djn(i),fjn(i),byn(i),dyn(i),fyn(i) )
    endif
    write (*,*), "alpha",alphas(i), "bjn", bjn(i)
  end do
  Do i=1, nfftr
    denom(i)= PI*rmax**2*bjn(i)**2
    !write (*,*), "denom", denom(i)
  end do
  ! work has stopped here
  Do k=1, nfftr
    r(k)= (k-1)*dr+0.5_num
    !write (*,*), "r", r(k)
  End do
  Do i =1,nfftr
    Do j=1, nfftr
      y(i,j)= alphas(i)*r(j)/rmax
      !write(*,*), "alphas", alphas(i), "r", r(j), "yij", y(i,j)
    end do
  end do
  Do i=1, nfftr
    Do j=1, nfftr
      if (y(i,j) .eq. 0.) then
        numer(i,j)=0.
      else
        Call jyndd ( p, y(i,j), numer(i,j), jn_d(i,j), jn_f(i,j),yn(i,j),yn_d(i,j),yn_f(i,j) )
      end if
      !write (*,*), "y",y(i,j), "numer", numer(i,j)
    end do
  end do
  denom1(:,1)=denom
  if (m .ne. 0) then
    Do i=1, nfftr
      invM(2:,i)= numer(2:,i)/ denom1(2:,1)
    end do
    if (p == m-1) then
      invM(1,:)= r**(m-1)/(PI *rmax**(m+1))
    else
      invM(1,:)=0_num
    end if
  else
    Do i =1, nfftr
      invM(1:,i)= numer(1:,i)/denom1(1:,1)
    end do

  end if

!  Do i= 1, nfftr
!    Do j=1,nfftr
!      WRITE(*,*) invM(i,j)
!    end do
!  end do
  DEALLOCATE (x,r, nu, alphas, bjn, djn, fjn, byn, dyn, fyn, denom, denom1,)
  DEALLOCATE (y,numer,invM,jn_d,jn_f,yn,yn_d,yn_f)
END SUBROUTINE hankel_matrix_init                          
! end of subroutine init matrix
! this should be a new subroutine

 
SUBROUTINE Hankel_Mtot_and_invMtot(invM_tot, Ma_tot)
  USE PICSAR_precision
  USE shared_data, ONLY: dx, nx, nxguards, nmodes
  !REAL(num), dimension(:,:,:), allocatable, intent (inout) :: invM_tot, Ma_tot
  REAL(num), dimension(:,:), allocatable :: invM , Ma
  INTEGER (idp) :: imode 
#if defined(LIBRARY)
   nfftr = nx+2*nxguards+1
#else
   !> When using picsar
   nfftr = nx+2*nxguards
#endif
  
  ALLOCATE (invM_tot(nfftr,nfftr,nmodes))
  ALLOCATE (Ma_tot(nfftr,nfftr,nmodes))
  ALLOCATE (Ma(nfftr,nfftr))
  ALLOCATE (invM(nfftr,nfftr))
  DO imode=1, nmodes 
    CALL hankel_matrix_init(nfftr,imode,p,rmax,invM)
    Ma=inv1(invM)
    invM_tot(:,:,imode)=invM
    Ma_tot(:,:,imode)= Ma    
    invM(:,:)=0._num
    Ma(:,:)=0._num
  END DO

END SUBROUTINE Hankel_Mtot_and_invMtot

  !Ma=inv1(invM)
  

  !write(*,*) "Ma"
  !Do i= 1, nfftr
  !  Do j=1,nfftr
  !    WRITE(*,*) Ma(i,j)
  !  end do
  !end do
  !write (*,*) "code stopped here"

SUBROUTINE get_Hfields(Ma_tot,nfftr)
  USE PICSAR_precision
  USE shared_data, ONLY: nmodes
  USE fields, ONLY:  el_h, ep_h, em_h, bl_h, bp_h, bm_h, jl_h, jp_h, jm_h, rhoold_h, rho_h
  USE fields, ONLY : el_f, er_f, et_f, bl_f, br_f, bt_f, jl_f, jr_f, jt_f, rho_f, rhoold_h
  USE Matrix_H, ONLY:  Ma_tot
  INTEGER (idp), intent(in) :: nfftr
  REAL(num), dimension(:,:), allocatable :: invM , Ma
  INTEGER (idp) :: imode 
   
  ALLOCATE (Ma_tot(nfftr,nfftr,nmodes))
  DO imode=1, nmodes
    Ma = Ma_tot(:,:,imode)
    call cpu_time ( t1 )
    Call  dgemm_example(el_f(:,:imode), Ma, el_h(:,:,imode), nfftr)
    Call  dgemm_example(ep_f(:,:imode), Ma, em_h(:,:,imode), nfftr)
    Call  dgemm_example(em_f(:,:imode), Ma, ep_h(:,:,imode), nfftr)
    Call  dgemm_example(bl_f(:,:imode), Ma, bl_h(:,:,imode), nfftr)
    Call  dgemm_example(bp_f(:,:imode), Ma, bm_h(:,:,imode), nfftr)
    Call  dgemm_example(bm_f(:,:imode), Ma, bp_h(:,:,imode), nfftr)
    Call  dgemm_example(jl_f(:,:imode), Ma, jl_h(:,:,imode), nfftr)
    Call  dgemm_example(jm_f(:,:imode), Ma, jp_h(:,:,imode), nfftr)
    Call  dgemm_example(jp_f(:,:imode), Ma, jm_h(:,:,imode), nfftr)
    Call  dgemm_example(rho_f(:,:imode), Ma, rho_h(:,:,imode), nfftr)
    Call  dgemm_example(rhoold_f(:,:imode), Ma, rhoold_h(:,:,imode), nfftr)
    call cpu_time ( t2 )
  write ( *, * ) 'Elapsed CPU time = ', t2 - t1

  !write (*,*) "dgemm successful"

  !DO i=1,nfftr
  !  Do j=1, nfftr
  !    write(*,*) c(i,j)
  !  end do
  !end do

  !Call  dgemm_example(c, invM ,d, nfftr)
  !write (*,*) "dgemm inv successful"

  !DO i=1,nfftr
  !  Do j=1, nfftr
  !    write(*,*) d(i,j)
  !  end do
  !end do

END SUBROUTINE get_Hfields


END MODULE Hankel

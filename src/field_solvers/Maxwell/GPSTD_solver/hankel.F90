
!MODULE Matrix_H
!  USE PICSAR_precision
!  REAL(num), POINTER, DIMENSION(:, :) :: Ma_1
!  REAL(num), POINTER, DIMENSION(:, :) :: invM_1
!  REAL(num), POINTER, DIMENSION(:, :) :: Ma
!  REAL(num), POINTER, DIMENSION(:, :) :: invM
!  REAL(num), POINTER, DIMENSION(:, :) :: Ma1
!  REAL(num), POINTER, DIMENSION(:, :) :: invM1
!END Matrix_H


MODULE Hankel !# do not parse

IMPLICIT NONE
CONTAINS

!SUBROUTINE allocate_matrix
!  USE fields, ONLY :: nx, nxguards
!  
!  nfftr= nx+2*nxguards
!  ALLOCATE (Ma_1(nfftr,nfftr))
!  ALLOCATE (invM_1(nfftr,nfftr))
!  ALLOCATE (Ma(nfftr,nfftr))
!  ALLOCATE (invM(nfftr,nfftr))
!  ALLOCATE (Ma1(nfftr,nfftr))
!  ALLOCATE (invM1(nfftr,nfftr))

!END SUBROUTINE allocate_matrix
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
FUNCTION inv(A) RESULT(Ainv)
  USE PICSAR_precision
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
COMPLEX(cpx) ,  INTENT(in) :: a(:,:)
REAL(num) , INTENT(in) :: b(:,:)
COMPLEX(cpx) , INTENT(out) :: c(:,:)
REAL (num) , dimension(size(a,1),size(a,2)) :: a_r, a_im , c_r, c_im
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
COMPLEX(cpx) :: ii

ii=DCMPLX(0._num,1_num)
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

a_r = DREAL(a)
a_im = DIMAG (a)
CALL DGEMM('N','N',M,N,K,alpha,a_r,lda,b,ldb,beta,c_r,ldc)
CALL DGEMM('N','N',M,N,K,alpha,a_im,lda,b,ldb,beta,c_im,ldc)

c(:,:)=c_r(:,:)+ ii* c_im(:,:)
END SUBROUTINE dgemm_example

SUBROUTINE hankel_matrix_init(nfftr,imode,p,rmax,invM)
  USE math_tools
  USE PICSAR_precision
  !> the matrix initialisation depend on two important variables 
  !> the considered mode imode and the order of hankel transform p
  !> Very important!! : p should be in {m-1, m, m+1}
  USE shared_data, ONLY : dx
  INTEGER (idp), intent(in) :: nfftr, imode, p 
  REAL (num), intent(in) :: rmax
  REAL (num), dimension(:,:),intent(inout):: invM
  REAL (num) :: PI  = 4 * atan (1.0_8)
  REAL (num) ::  dr
  !USE shared_data, ONLY: dx, nx
  INTEGER (idp) :: k, j,i, m,p_denom
  REAL (num) , dimension (:), allocatable :: x,r,nu,alphas,bjn,djn,fjn,byn,dyn,fyn,denom, rj1,ry0,ry1
  REAL (num) , dimension (:,:), allocatable :: denom1
  !REAL (num) , dimension (:,:) , allocatable:: alphas1
  REAL (num) , dimension (:,:), allocatable :: y,numer, jn_d, jn_f,yn,yn_d,yn_f
  !REAL (num) , dimension (:,:), allocatable :: invM, Ma
  !REAL (kind =4),dimension (:,:), allocatable :: Ma
  !dr=1.
  !p= 0
  !> For writing simplicity we use m instead of imode 
  m= imode
  dr= dx
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
  ALLOCATE (rj1(nfftr))
  ALLOCATE (ry0(nfftr))
  ALLOCATE (ry1(nfftr))
  !ALLOCATE (r1(nfftr,1))
  ALLOCATE (denom1(nfftr,1))
  !ALLOCATE (alphas1(1,nfftr))
  ALLOCATE (y(nfftr,nfftr))
  ALLOCATE (numer(nfftr,nfftr))
  !ALLOCATE (invM(nfftr,nfftr))
  !ALLOCATE (Ma(nfftr,nfftr))
  ALLOCATE (jn_d(nfftr,nfftr))
  ALLOCATE (jn_f(nfftr,nfftr))
  ALLOCATE (yn(nfftr,nfftr))
  ALLOCATE (yn_d(nfftr,nfftr))
  ALLOCATE (yn_f(nfftr,nfftr))
!  ALLOCATE (b(nfftr,nfftr))
!  ALLOCATE (c(nfftr,nfftr))
!  ALLOCATE (d(nfftr,nfftr))
!  Do i=1, nfftr
!    Do j=1, nfftr
!      b(i,j)=1.0
!    end do
!  end do


  IF (p .eq.m ) THEN
    p_denom= p+1
  ELSE
    p_denom =p
  ENDIF
  CALL  jyzo ( m, nfftr, nu, rj1, ry0, ry1 )
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
  DEALLOCATE (y,numer,jn_d,jn_f,yn,yn_d,yn_f)
  DEALLOCATE (rj1,ry0, ry1)
END SUBROUTINE hankel_matrix_init                          
! end of subroutine init matrix
! this should be a new subroutine

 
SUBROUTINE Hankel_M_and_invM()
  USE PICSAR_precision
  USE fields , ONLY : Ma, Ma_1, Ma1, invM, invM_1, invM1
  USE shared_data, ONLY: dx, nx, nmodes, xmax
  USE fields, ONLY : nxguards
  !REAL(num), dimension(:,:,:), allocatable, intent (inout) :: invM_tot, Ma_tot
  !REAL(num), dimension(:,:), allocatable :: invM , Ma
  INTEGER (idp) :: imode,p,nfftr  
#if defined(LIBRARY)
   nfftr = nx+2*nxguards+1
#else
   !> When using picsar
   nfftr = nx+2*nxguards
#endif
  
  !ALLOCATE (invM_tot(nfftr,nfftr,nmodes))
  !ALLOCATE (Ma_tot(nfftr,nfftr,nmodes))
  !ALLOCATE (Ma(nfftr,nfftr))
  !ALLOCATE (invM(nfftr,nfftr))
  DO imode=1, nmodes 
    DO p=imode-1, imode+1
      CALL hankel_matrix_init(nfftr,imode,p,xmax,invM_1)
      CALL hankel_matrix_init(nfftr,imode,p,xmax,invM)
      CALL hankel_matrix_init(nfftr,imode,p,xmax,invM1)
    END DO
    Ma=inv(invM)
    Ma_1 =inv(invM_1)
    Ma1 = inv(invM1)
    !invM_tot(:,:,imode)=invM
    !Ma_tot(:,:,imode)= Ma    
    !invM(:,:)=0._num
    !Ma(:,:)=0._num
  END DO

END SUBROUTINE Hankel_M_and_invM

  !Ma=inv1(invM)
  

  !write(*,*) "Ma"
  !Do i= 1, nfftr
  !  Do j=1,nfftr
  !    WRITE(*,*) Ma(i,j)
  !  end do
  !end do
  !write (*,*) "code stopped here"

SUBROUTINE get_Hfields()
  USE PICSAR_precision
  USE fourier_psaotd
  USE shared_data, ONLY: nmodes
  USE fields, ONLY:  el_h, ep_h, em_h, bl_h, bp_h, bm_h, jl_h, jp_h, jm_h, rhoold_h, rho_h
  !USE fields, ONLY : el_f, er_f, et_f, bl_f, br_f, bt_f, jl_f, jr_f, jt_f, rho_f, rhoold_h
  USE fields, ONLY : el_f, ep_f, em_f, bl_f, bp_f, bm_f, jl_f, jp_f, jm_f, rho_f, rhoold_f
  USE fields, ONLY : Ma, Ma1, Ma_1
  USE shared_data, ONLY:  nx
  USE fields, ONLY:  nxguards
  !INTEGER (idp), intent(in) :: nfftr
  !REAL(num), dimension(:,:), allocatable :: Ma
  INTEGER(idp) :: nfftr
  INTEGER (idp) :: imode
  REAL(num) :: t1, t2
  COMPLEX(cpx),POINTER, DIMENSION(:, :) :: el_h_ptr, ep_h_ptr, em_h_ptr, bl_h_ptr, bp_h_ptr, &
                            bm_h_ptr, jl_h_ptr, jp_h_ptr, jm_h_ptr, rhoold_h_ptr, rho_h_ptr, &
                                 el_f_ptr, ep_f_ptr, em_f_ptr, bl_f_ptr, bp_f_ptr, bm_f_ptr, &
                                    jl_f_ptr, jp_f_ptr, jm_f_ptr, rho_f_ptr, rhoold_f_ptr 
  !REAL (num), dimension (:,:,:), allocatable :: el_h_r, ep_h_r, em_h_r, bl_h_r, bp_h_r, bm_h_r, jl_h_r, &
  !                            jp_h_r, jm_h_r, rhoold_h_r, rho_h_r, el_h_im, ep_h_im, em_h_im, bl_h_im,  &
  !                            bp_h_im, bm_h_im, jl_h_im, jp_h_im, jm_h_im, rhoold_h_im, rho_h_im
 
  !REAL (num), dimension (:,:), allocatable :: f_h_r, f_f_r, f_h_im, f_f_im
  !ALLOCATE (Ma(nfftr,nfftr))
  !ALLOCATE (f_h_r(size(el_h,1), size(el_h,2)))
  !ALLOCATE (f_f_r(size(el_f,1), size(el_f,2)))
  !ALLOCATE (f_h_im(size(el_h,1), size(el_h,2)))
  !ALLOCATE (f_f_im(size(el_f,1), size(el_f,2)))
#if defined(LIBRARY)
    nfftr=nx+2*nxguards+1
#else
    nfftr=nx+2*nxguards
#endif
  Call get_Ffields_AM_rz()
  
  DO imode=1, nmodes
    !Ma = Ma_tot(:,:,imode)
    Call Hankel_M_and_invM()
    call cpu_time ( t1 )
    el_h_ptr=>el_h(:,:,imode)
    ep_h_ptr=>ep_h(:,:,imode)
    em_h_ptr=>em_h(:,:,imode)
    bl_h_ptr=>bl_h(:,:,imode)
    bp_h_ptr=>bp_h(:,:,imode)
    bm_h_ptr=>bm_h(:,:,imode)
    jl_h_ptr=>bm_h(:,:,imode)
    jp_h_ptr=>jp_h(:,:,imode)
    jm_h_ptr=>jm_h(:,:,imode)
    rhoold_h_ptr=>rhoold_h(:,:,imode)
    rho_h_ptr=>rho_h(:,:,imode)
    el_f_ptr=>el_f(:,:,imode)
    ep_f_ptr=>ep_f(:,:,imode)
    em_f_ptr=>em_f(:,:,imode)
    bl_f_ptr=>bl_f(:,:,imode)
    bp_f_ptr=>bp_f(:,:,imode)
    bm_f_ptr=>bm_f(:,:,imode)
    jl_f_ptr=>jl_f(:,:,imode)
    jp_f_ptr=>jp_f(:,:,imode)
    jm_f_ptr=>jm_f(:,:,imode)
    rho_f_ptr=>rho_f(:,:,imode)
    rhoold_f_ptr=>rhoold_f(:,:,imode)
    Call  dgemm_example(el_f_ptr, Ma, el_h_ptr, nfftr)
    Call  dgemm_example(ep_f_ptr, Ma_1, em_h_ptr, nfftr)
    Call  dgemm_example(em_f_ptr, Ma1, ep_h_ptr, nfftr)
    Call  dgemm_example(bl_f_ptr, Ma, bl_h_ptr, nfftr)
    Call  dgemm_example(bp_f_ptr, Ma_1, bm_h_ptr, nfftr)
    Call  dgemm_example(bm_f_ptr, Ma1, bp_h_ptr, nfftr)
    Call  dgemm_example(jl_f_ptr, Ma, jl_h_ptr, nfftr)
    Call  dgemm_example(jm_f_ptr, Ma1, jp_h_ptr, nfftr)
    Call  dgemm_example(jp_f_ptr, Ma_1, jm_h_ptr, nfftr)
    Call  dgemm_example(rho_f_ptr, Ma, rho_h_ptr, nfftr)
    Call  dgemm_example(rhoold_f_ptr, Ma, rhoold_h_ptr, nfftr)
    call cpu_time ( t2 )
  write ( *, * ) 'Elapsed CPU time = ', t2 - t1
  
  END DO
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


SUBROUTINE get_Hfields_inv()
  USE PICSAR_precision
  USE shared_data, ONLY: nmodes
  USE fields, ONLY:  el_h, ep_h, em_h, bl_h, bp_h, bm_h, jl_h, jp_h, jm_h, rhoold_h, rho_h
  !USE fields, ONLY : el_f, er_f, et_f, bl_f, br_f, bt_f, jl_f, jr_f, jt_f, rho_f, rhoold_h
  USE fields, ONLY:  el_h_inv, ep_h_inv, em_h_inv, bl_h_inv, bp_h_inv, bm_h_inv 
  !                   jl_h_inv, jp_h_inv, jm_h_inv, rhoold_h_inv, rho_h_inv
  USE fields, ONLY : invM, invM1, invM_1
  USE shared_data, ONLY:  nx
  USE fields, ONLY:  nxguards
  COMPLEX(cpx),POINTER, DIMENSION(:, :) :: el_h_ptr, ep_h_ptr, em_h_ptr, bl_h_ptr, bp_h_ptr, bm_h_ptr,&
                                  el_h_inv_ptr, ep_h_inv_ptr, em_h_inv_ptr, bl_h_inv_ptr, bp_h_inv_ptr, bm_h_inv_ptr
  !INTEGER(idp), INTENT(IN):: nfftr
  !REAL(num), dimension(:,:), allocatable :: invM 
  INTEGER (idp) :: imode, nfftr
  REAL (num) :: t1, t2
  !ALLOCATE (invM(nfftr,nfftr))
#if defined(LIBRARY)
    nfftr=nx+2*nxguards+1
#else
    nfftr=nx+2*nxguards
#endif


  DO imode=1, nmodes
   ! Ma = Ma_tot(:,:,imode)
    Call Hankel_M_and_invM()
    call cpu_time ( t1 )
    el_h_ptr=>el_h(:,:,imode)
    ep_h_ptr=>ep_h(:,:,imode)
    em_h_ptr=>em_h(:,:,imode)
    bl_h_ptr=>bl_h(:,:,imode)
    bp_h_ptr=>bp_h(:,:,imode)
    bm_h_ptr=>bm_h(:,:,imode)
    el_h_inv_ptr=>el_h_inv(:,:,imode)
    ep_h_inv_ptr=>ep_h_inv(:,:,imode)
    em_h_inv_ptr=>em_h_inv(:,:,imode)
    bl_h_inv_ptr=>bl_h_inv(:,:,imode)
    bp_h_inv_ptr=>bp_h_inv(:,:,imode)
    bm_h_inv_ptr=>bm_h_inv(:,:,imode)
    Call  dgemm_example(el_h_ptr, invM, el_h_inv_ptr, nfftr)
    Call  dgemm_example(em_h_ptr, invM_1, em_h_inv_ptr, nfftr)
    Call  dgemm_example(ep_h_ptr, invM1, ep_h_inv_ptr, nfftr)
    Call  dgemm_example(bl_h_ptr, invM, bl_h_inv_ptr, nfftr)
    Call  dgemm_example(bm_h_ptr, invM_1, bm_h_inv_ptr, nfftr)
    Call  dgemm_example(bp_h_ptr, invM1, bp_h_inv_ptr, nfftr)
!    Call  dgemm_example(jl_h(:,:,imode), invM, jl_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(jp_h(:,:,imode), invM1, jp_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(jm_h(:,:,imode), invM_1, jm_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(rho_h(:,:,imode), invM, rho_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(rhoold_h(:,:,imode), invM, rho_hold_inv(:,:,imode), nfftr)
    call cpu_time ( t2 )
  write ( *, * ) 'Elapsed CPU time = ', t2 - t1
  END DO

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

END SUBROUTINE get_Hfields_inv

END MODULE Hankel


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
USE math_tools
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
  IMPLICIT NONE
  real(num), dimension(:,:), intent(in) :: A
  real(num), dimension(size(A,1),size(A,2)) :: Ainv

  real(num), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer(isp), dimension(size(A,1)) :: ipiv   ! pivot indices
  integer(isp)  :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)
  !write (*,*) "n =", n
  !write(*,*) "DGETRF"
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  !write(*,*) "DGETRI"
  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
END FUNCTION inv


! Matrix multiplication for RZ

SUBROUTINE dgemm_example( a,b,c, nfftr, nffty)
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
INTEGER (idp), INTENT(in) :: nfftr,nffty
INTEGER  (idp):: I, J
COMPLEX(cpx) , dimension(:,:), INTENT(in) :: a
REAL(num) ,  dimension(:,:), INTENT(in) :: b
COMPLEX(cpx) , dimension(:,:), INTENT(out) :: c
COMPLEX (cpx), DIMENSION(size(c,2), size(c,1)) :: d 
REAL (num) , dimension(size(a,1),size(a,2)) :: a_r, a_im 
REAL (num) , dimension(size(c,2),size(c,1)) :: d_r, d_im
INTEGER (isp) :: M    ! Number of rows of matrix op(a)
INTEGER  (isp):: N    ! Number of columns of matrix op(b)
INTEGER  (isp):: K    ! Number of columns of matrix op(a) and rows op(b)
INTEGER (isp):: lda  ! n entry, LDA specifies the first dimension of A as declared
             ! in the calling (sub) program. When  TRANSA = 'N' or 'n' then
             ! LDA must be at least  max( 1, m ), otherwise  LDA must be at
             ! least  max( 1, k ).  
INTEGER (isp) :: ldb  ! On entry, LDB specifies the first dimension of B as declared
             ! in the calling (sub) program. When  TRANSB = 'N' or 'n' then
             ! LDB must be at least  max( 1, k ), otherwise  LDB must be at
             ! least  max( 1, n ). 
INTEGER (isp) :: ldc  ! On entry, LDC specifies the first dimension of C as declared
             ! in  the  calling  (sub)  program.   LDC  must  be  at  least
             ! max( 1, m ).
REAL (num) ::  alpha, beta
COMPLEX(cpx) :: ii
ii=DCMPLX(0._num,1_num)

a_r= 0.0_num
a_im =0.0_num
d=0.0_num
d_im=0.0_num
d_r =0.0_num
#ifdef _OPENMP
 !write(*,*)"Will use: ", omp_get_num_procs(), " of ", &
 !           omp_get_max_threads(), " available"
#endif
M= nffty
N= nfftr
K=nfftr
lda=K
ldb=K
ldc=M
alpha=1.0
beta=0.0

a_r = REAL(a)
a_im = IMAG (a)
CALL DGEMM('T','N',M,N,K,alpha,a_r,lda,b,ldb,beta,d_r,ldc)
CALL DGEMM('T','N',M,N,K,alpha,a_im,lda,b,ldb,beta,d_im,ldc)
d=d_r+ ii* d_im
c= transpose (d)
!c(10,10)= 100.
!write (0,*) "c" , c(10,10)
END SUBROUTINE dgemm_example
SUBROUTINE hankel_matrix_init(nfftr,imode,p,rmax,invM)
  USE math_tools
  USE PICSAR_precision
  !> the matrix initialisation depend on two important variables 
  !> the considered mode imode and the order of hankel transform p
  !> Very important!! : p should be in {m-1, m, m+1}
  USE shared_data, ONLY : dx
  IMPLICIT NONE
  INTEGER (idp), intent(in) :: imode, p 
  INTEGER (idp), intent(in) :: nfftr
  REAL (num), intent(in) :: rmax
  REAL (num), dimension(:,:),intent(inout):: invM
  REAL (num) :: PI  = 4 * atan (1.0_8)
  REAL (num) ::  dr
  !USE shared_data, ONLY: dx, nx
  INTEGER (idp) :: k, j,i, m,p_denom
  REAL (num) , dimension (:), allocatable :: x,r,nu, nu1,alphas,bjn,djn,fjn,byn,dyn,fyn,denom, rj1,ry0,ry1
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
  !rmax=nfftr*dx
  ALLOCATE (r(nfftr))
  ALLOCATE (x(nfftr))
  ALLOCATE (nu(nfftr))
  ALLOCATE (nu1(nfftr-1))
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

  invM=0.
  IF (p .eq.m ) THEN
    p_denom= p+1
  ELSE
    p_denom =p
  ENDIF
  if (m .ne.0) then
    alphas(1)= 0._num
    CALL  jyzo ( m, nfftr-1, nu1, rj1, ry0, ry1 )
    alphas(2:nfftr)=nu1(1:nfftr-1)
  else
    CALL  jyzo ( m, nfftr, nu, rj1, ry0, ry1 )
    alphas=nu
  endif
 ! Do i =1,nfftr
 !   write(*,*),"alphas",  alphas(i)
 ! end do
  Do i =1, nfftr
    !if (alphas(i) .eq. 0._num) then
    !  bjn(i)=0.
    !else
    if (p_denom .ge. 0) then
      Call JYNDD ( p_denom, alphas(i), bjn(i), djn(i),fjn(i),byn(i),dyn(i),fyn(i) )
    else
      Call JYNDD ( -p_denom, alphas(i), bjn(i), djn(i),fjn(i),byn(i),dyn(i),fyn(i) )
      bjn(i)=(-1.)**(-(p_denom))*bjn(i)
    end if
    !Call jyndd ( p_denom, alphas(i), bjn(i),djn(i),fjn(i),byn(i),dyn(i),fyn(i) )
    !endif
    !write (*,*), "alpha",alphas(i), "bjn", bjn(i)
  end do
  Do i=1, nfftr
    denom(i)= PI*rmax**2*bjn(i)**2
    !write (*,*), "denom", denom(i)
  end do
  ! work has stopped here
  Do k=1, nfftr
    !r(k)= (k-1)*dr+0.5_num
    r(k)= dr*((k-1)+0.5_num)
    !write (*,*), "r", r(k)
  End do
  Do j =1,nfftr
    Do i=1, nfftr
      y(i,j)= alphas(i)*r(j)/rmax
      !write(*,*), "alphas", alphas(i), "r", r(j), "yij", y(i,j)
    end do
  end do
  Do j=1, nfftr
    Do i=1, nfftr
      !if (y(i,j) .eq. 0.) then
      !  numer(i,j)=0.
      !else
      if (p .ge. 0) then
        Call jyndd ( p, y(i,j), numer(i,j), jn_d(i,j), jn_f(i,j),yn(i,j),yn_d(i,j),yn_f(i,j) )
      else
        Call jyndd ( -p, y(i,j), numer(i,j), jn_d(i,j), jn_f(i,j),yn(i,j),yn_d(i,j),yn_f(i,j) )
        !write (0,*) "jn_d" , jn_d(i,j) , (-1.)**(abs(p))
        numer(i,j)=(-1.)**(-(p))*numer(i,j)
        !write (0,*) "jn_d1" , jn_d(i,j)
      endif
      !Call jyndd ( p, y(i,j), numer(i,j), jn_d(i,j), jn_f(i,j),yn(i,j),yn_d(i,j),yn_f(i,j) )
      !end if
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
      invM(1,:)=0.0_num
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
  DEALLOCATE (x,r, nu,nu1, alphas, bjn, djn, fjn, byn, dyn, fyn, denom, denom1)
  DEALLOCATE (y,numer,jn_d,jn_f,yn,yn_d,yn_f)
  DEALLOCATE (rj1,ry0, ry1)
END SUBROUTINE hankel_matrix_init                          
! end of subroutine init matrix
! this should be a new subroutine

subroutine rmat_svd_lapack ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! R8MAT_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    Moreover, the nonzero entries of S are positive, and appear
!    in order, from largest magnitude to smallest.
!
!    This routine calls the LAPACK routine DGESVD to compute the
!    factorization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
  use picsar_precision
  implicit none

  integer ( isp ), intent(in):: m
  integer ( isp ), intent (in) :: n

  real ( num ), dimension (:,:), intent(in) ::  a(m,n)
  real ( num ), dimension (:,:):: a_copy(m,n)
  integer ( isp ) ::i
  integer ( isp ) ::info
  integer ( isp ) ::lda
  integer ( isp ) ::ldu
  integer ( isp ) ::ldv
  character ::jobu
  character ::jobv
  integer ( isp ) ::lwork
  real ( num ), dimension (:) :: sdiag(min(m,n))
  real ( num ), dimension (:,:), intent(out) ::  s(m,n)
  real ( num ), dimension (:,:), intent(out) ::  u(m,m)
  real ( num ), dimension (:,:), intent(out) ::  v(n,n)
  real ( num ), allocatable, dimension ( : ) :: work

  lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

  allocate ( work(1:lwork) )
!
!  Compute the eigenvalues and eigenvectors.
!
  jobu = 'A'
  jobv = 'A'
  lda = m
  ldu = m
  ldv = n
!
!  The input matrix is destroyed by the routine.  Since we need to keep
!  it around, we only pass a copy to the routine.
!
  a_copy(1:m,1:n) = a(1:m,1:n)

  call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_SVD_LAPACK - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    return
  end if
!
!  Make the MxN matrix S from the diagonal values in SDIAG.
!
  s(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    s(i,i) = sdiag(i)
  end do
!
!  Transpose V.
!
  v = transpose ( v )

  deallocate ( work )

  return
end subroutine rmat_svd_lapack

 
subroutine pseudo_inverse ( m, n, u, s, v, a_pseudo )

!*****************************************************************************80
!
!! PSEUDO_INVERSE computes the pseudoinverse.
!
!  Discussion:
!
!    Given the singular value decomposition of a real MxN matrix A:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    the pseudo inverse is the NxM matrix A+ with the form
!
!      A+ = V * S+ * U'

!    where
!
!      S+ is the NxM matrix whose nonzero diagonal elements are
!      the inverses of the corresponding diagonal elements of S.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
!    Output, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
!

  use picsar_precision
  implicit none

  integer ( isp ), intent(in):: m
  integer ( isp ) , intent (in) ::n

  real ( num ), dimension (:,:), intent (out) :: a_pseudo(n,m)
  integer ( idp ) :: i
  real ( num ), dimension (:,:), intent (in):: s(m,n)
  real ( num ) , dimension (:,:)::sp(n,m)
  real ( num ), dimension (:,:), intent (in):: u(m,m)
  real ( num ), dimension (:,:), intent (in):: v(n,n)

  sp(1:n,1:m) = 0.0D+00
  do i = 1, min ( m, n )
    if ( s(i,i) /= 0.0D+00 ) then
      sp(i,i) = 1.0D+00 / s(i,i)
    end if
  end do
  a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )

  return
end subroutine pseudo_inverse



SUBROUTINE Hankel_M_and_invM(imode)
  USE PICSAR_precision
  USE fields , ONLY : Ma, Ma_1, Ma1, invM, invM_1, invM1
  USE shared_data, ONLY: dx, nx, nmodes, xmax
  USE fields, ONLY : nxguards
  IMPLICIT NONE
  !REAL(num), dimension(:,:,:), allocatable, intent (inout) :: invM_tot, Ma_tot
  !REAL(num), dimension(:,:), allocatable :: invM , Ma
  INTEGER (idp), INTENT (IN) :: imode
  INTEGER (idp) :: p, k,m
  INTEGER (idp)  :: nfftr
  REAL (NUM) :: rmax 
  real ( num ), allocatable, dimension ( :, : ) :: a
  real ( num ), allocatable, dimension ( :, : ) :: a_pseudo
  real ( num ), allocatable, dimension ( :, : ) :: s
  real ( num ), allocatable, dimension ( :, : ) :: u
  real ( num ), allocatable, dimension ( :, : ) :: v
  real ( num ), allocatable, dimension ( :, : ) :: a1
  real ( num ), allocatable, dimension ( :, : ) :: a_pseudo1
  real ( num ), allocatable, dimension ( :, : ) :: s1
  real ( num ), allocatable, dimension ( :, : ) :: u1
  real ( num ), allocatable, dimension ( :, : ) :: v1
#if defined(LIBRARY)
   nfftr = nx+2*nxguards+1
#else
   !> When using picsar
   nfftr = nx+2*nxguards
   write (*,*) "using picsar"
#endif
  

  allocate ( a(1:nfftr-1,1:nfftr) )
  allocate ( a_pseudo(1:nfftr,1:nfftr-1) )
  allocate ( u(1:nfftr-1,1:nfftr) )
 ! allocate ( u2(1:m,1:m) )
  allocate ( s(1:nfftr-1,1:nfftr) )
 ! allocate ( s2(1:m,1:n) )
  allocate ( v(1:nfftr,1:nfftr) )
  allocate ( a1(1:nfftr-1,1:nfftr) )
  allocate ( a_pseudo1(1:nfftr,1:nfftr-1) )
  allocate ( u1(1:nfftr-1,1:nfftr) )
 ! allocate ( u2(1:m,1:m) )
  allocate ( s1(1:nfftr-1,1:nfftr) )
 ! allocate ( s2(1:m,1:n) )
  allocate ( v1(1:nfftr,1:nfftr) )
  !ALLOCATE (invM_tot(nfftr,nfftr,nmodes))
  !ALLOCATE (Ma_tot(nfftr,nfftr,nmodes))
  !ALLOCATE (Ma(nfftr,nfftr))
  !ALLOCATE (invM(nfftr,nfftr))
  !DO imode=0, nmodes-1 
  !a=0.0_num
  !a_pseudo=0.0_num
  !u= 0.0_num
  !s=0.0_num
  !v=0.0_num
  !a1=0.0_num
  !a_pseudo1=0.0_num
  !u1= 0.0_num
  !s1=0.0_num
  !v1= 0.0_num
  rmax= dx *nfftr
  !rmax=xmax
  DO p=imode-1, imode+1
    SELECT CASE(p-imode)
      CASE(-1)
        CALL hankel_matrix_init(nfftr,imode,p,rmax,invM_1)
        Ma_1 =inv(invM_1)
      CASE(0)
        CALL hankel_matrix_init(nfftr,imode,p,rmax,invM)
        IF (imode .NE. 0) THEN
          a(1:nfftr-1,1:nfftr)= invM(2:nfftr,1:nfftr)
          call rmat_svd_lapack ( INT((nfftr-1),isp), INT(nfftr,isp), a, u, s, v )
          call pseudo_inverse ( INT((nfftr-1),isp), INT(nfftr,isp), u, s, v, a_pseudo )
          Ma(1:nfftr,2:nfftr)= a_pseudo(1:nfftr,1:nfftr-1)
          Ma(:,1)=0.0_num
          !write (*,*) "pseudo"
        ELSE
          Ma=inv(invM)
        END IF 
      CASE(1)
        CALL hankel_matrix_init(nfftr,imode,p,rmax,invM1)
        IF (imode .NE. 0) THEN
          !write (0,*), "invM1=======",invM1 
          a1(1:nfftr-1,1:nfftr)= invM1(2:nfftr,1:nfftr)
          call rmat_svd_lapack ( INT((nfftr-1),isp), INT(nfftr,isp), a1, u1, s1, v1 )
          call pseudo_inverse ( INT((nfftr-1),isp), INT(nfftr,isp), u1, s1, v1, a_pseudo1 )
          Ma1(1:nfftr,2:nfftr)= a_pseudo1(1:nfftr,1:nfftr-1)
          Ma1(:,1)=0.0_num
          !write (*,*) "pseudo"
        ELSE
          Ma1=inv(invM1)
        END IF
    END SELECT
  END DO

  DEALLOCATE (a, a_pseudo, u, s, v, a1, a_pseudo1, u1, s1, v1)

!    DO p=imode-1, imode+1
!      CALL hankel_matrix_init(nfftr,imode,p,xmax,invM_1)
!      CALL hankel_matrix_init(nfftr,imode,p,xmax,invM)
!      CALL hankel_matrix_init(nfftr,imode,p,xmax,invM1)
!    END DO
!    write(*,*) "doing inv"

!  if ((imode .ne. 0) .and. (p .ne. (imode-1))) then
!    a(1:nfftr-1,1:nfftr)= invM(2:nfftr,1:nfftr)
!    call rmat_svd_lapack ( (nfftr-1), nfftr, a, u, s, v )
!    call pseudo_inverse ( (nfftr-1), nfftr, u, s, v, a_pseudo )
!    Ma(1:nfftr,2:nfftr)= a_pseudo(1:nfftr,1:nfftr-1)
!    Ma(:,1)=0.
!    write (*,*) "pseudo"
!  else
!    Ma=inv1(invM)
!  end if
!  write(*,*) "Ma"
!  Do i= 1, nfftr
!    Do j=1,nfftr
!      WRITE(*,*) "Ma",  Ma(i,j) , "i=", i, "j=", j
!    end do
!  end do
!  write (*,*) "code stopped here"


!    Ma=inv(invM)
!    write (*,*) "inv done"
!    Ma_1 =inv(invM_1)
!    Ma1 = inv(invM1)
    !invM_tot(:,:,imode)=invM
    !Ma_tot(:,:,imode)= Ma    
    !invM(:,:)=0._num
    !Ma(:,:)=0._num
 ! END DO


       ! write (0,*), "Ma =================="
       ! DO k=1,nx
       !   DO m=1, nx
       !     write (0,*) , "k= ", k, "m= ", m, "Ma", Ma(k,m)
       !   END DO
       ! END DO
       ! write (0,*), "Ma_1 =================="
       ! DO k=1,nx
       !   DO m=1,nx
       !     write (0,*) , "k= ", k, "m= ", m, "Ma_1", Ma_1(k,m)
       !   END DO
       ! END DO
       ! write (0,*), "Ma1 =================="
       ! DO k=1,nx
       !   DO m=1,nx
       !     write (0,*) , "k= ", k, "m= ", m, "Ma1", Ma1(k,m)
       !   END DO
       ! END DO
       ! write (0,*), "invM =================="
       ! DO k=1,nx
       !   DO m=1,nx
       !     write (0,*) , "k= ", k, "m= ", m, "invM", invM(k,m)
       !   END DO
       ! END DO
       ! write (0,*), "invM_1 =================="
       ! DO k=1,nx
       !   DO m=1, nx
       !     write (0,*) , "k= ", k, "m= ", m, "invM_1", invM_1(k,m)
       !   END DO
       ! END DO
       ! write (0,*), "invM1 =================="
       ! DO k=1,nx
       !   DO m=1,nx
       !     write (0,*) , "k= ", k, "m= ", m, "invM1", invM1(k,m)
       !   END DO
       ! END DO

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
  USE fields !, ONLY:  el_h, ep_h, em_h, bl_h, bp_h, bm_h, jl_h, jp_h, jm_h, rhoold_h, rho_h
  !USE fields, ONLY : el_f, er_f, et_f, bl_f, br_f, bt_f, jl_f, jr_f, jt_f, rho_f, rhoold_h
  !USE fields, ONLY : el_f, ep_f, em_f, bl_f, bp_f, bm_f, jl_f, jp_f, jm_f, rho_f, rhoold_f
  !USE fields, ONLY : Ma, Ma1, Ma_1
  USE shared_data, ONLY:  nx, ny, nmodes 
  USE fields, ONLY:  nxguards, nyguards
  !INTEGER (idp), intent(in) :: nfftr
  !REAL(num), dimension(:,:), allocatable :: Ma
  IMPLICIT NONE
  INTEGER(idp) :: nfftr, nffty
  INTEGER (idp) :: imode, i,j 
  REAL(num) :: t1, t2
  complex(cpx)::ii
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
    nffty=ny+2*nyguards+1
#else
    nfftr=nx+2*nxguards
    nffty=ny+2*nyguards
#endif
  ii= DCMPLX(0.0_NUM,1.0_NUM)
  !el_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !em_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !ep_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !bl_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !bm_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !bp_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !jl_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !jm_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !jp_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !rho_h = DCMPLX(0.0_NUM, 0.0_NUM)
  !rhoold_h = DCMPLX(0.0_NUM, 0.0_NUM)
  Call get_Ffields_AM_rz()
  !el_f=el_c
  !em_f= (er_c+ii*et_c)/2.0_num
  !ep_f=(er_c-ii*et_c)/2.0_num
  !bm_f= (br_c+ii*bt_c)/2.0_num
  !bp_f=(br_c-ii*bt_c)/2.0_num  
  !write (*,*) "==========START GET H FIELD HERE AFTER GET F FIELD =============================="
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*)"er_c", er_c(i,j,:)
  ! end do
  !end do
    
  !write (*,*) "get_Ffields_AM_rz  em_f" !, MAXVAL(abs(em_f))
  !DO i=1, nfftr
  ! DO j=1, ny
     !write (0,*), em_f(i,j,:)
  ! end do
  !end do

  !write (*,*) "get_Ffields_AM_rz  ep_f", MAXVAL(abs(ep_f))
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_f(i,j,:)
  ! end do
  !end do
  DO imode=1, nmodes
    !Ma = Ma_tot(:,:,imode)
    Call Hankel_M_and_invM(imode-1)
  !  !write (0,*), "Ma ======" , Ma
  !  !write (0,*), "Ma_1 ======" , Ma_1
  !  !write (0,*), "Ma1 ======" , Ma1
  !  !call cpu_time ( t1 )
  !  el_h_ptr=>el_h(:,:,imode)
  !  ep_h_ptr=>ep_h(:,:,imode)
  !  em_h_ptr=>em_h(:,:,imode)
  !  bl_h_ptr=>bl_h(:,:,imode)
  !  bp_h_ptr=>bp_h(:,:,imode)
  !  bm_h_ptr=>bm_h(:,:,imode)
  !  jl_h_ptr=>jl_h(:,:,imode)
  !  jp_h_ptr=>jp_h(:,:,imode)
  !  jm_h_ptr=>jm_h(:,:,imode)
  !  rhoold_h_ptr=>rhoold_h(:,:,imode)
  !  rho_h_ptr=>rho_h(:,:,imode)
  !  el_f_ptr=>el_f(:,:,imode)
  !  ep_f_ptr=>ep_f(:,:,imode)
  !  em_f_ptr=>em_f(:,:,imode)
  !  bl_f_ptr=>bl_f(:,:,imode)
  !  bp_f_ptr=>bp_f(:,:,imode)
  !  bm_f_ptr=>bm_f(:,:,imode)
  !  jl_f_ptr=>jl_f(:,:,imode)
  !  jp_f_ptr=>jp_f(:,:,imode)
  !  jm_f_ptr=>jm_f(:,:,imode)
  !  rho_f_ptr=>rho_f(:,:,imode)
  !  rhoold_f_ptr=>rhoold_f(:,:,imode)
  !  Call  dgemm_example(el_f_ptr, Ma, el_h_ptr, nfftr, nffty)
  !  Call  dgemm_example(ep_f_ptr, Ma_1, em_h_ptr, nfftr, nffty)
  !  Call  dgemm_example(em_f_ptr, Ma1, ep_h_ptr, nfftr, nffty)
  !  Call  dgemm_example(bl_f_ptr, Ma, bl_h_ptr, nfftr, nffty)
  !  Call  dgemm_example(bp_f_ptr, Ma_1, bm_h_ptr, nfftr, nffty)
  !  Call  dgemm_example(bm_f_ptr, Ma1, bp_h_ptr, nfftr, nffty)
    !Call  dgemm_example(jl_f_ptr, Ma, jl_h_ptr, nfftr, nffty)
    !Call  dgemm_example(jm_f_ptr, Ma1, jp_h_ptr, nfftr, nffty)
    !Call  dgemm_example(jp_f_ptr, Ma_1, jm_h_ptr, nfftr, nffty)
    !Call  dgemm_example(rho_f_ptr, Ma, rho_h_ptr, nfftr, nffty)
    !Call  dgemm_example(rhoold_f_ptr, Ma, rhoold_h_ptr, nfftr, nffty)
    !call cpu_time ( t2 )

    Call  dgemm_example(el_f(:,:,imode), Ma, el_h(:,:,imode), nfftr, nffty )
    Call  dgemm_example(em_f(:,:,imode), Ma_1, em_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(ep_f(:,:,imode), Ma1, ep_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bl_f(:,:,imode), Ma, bl_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bm_f(:,:,imode), Ma_1, bm_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bp_f(:,:,imode), Ma1, bp_h(:,:,imode), nfftr, nffty)
    !Call  dgemm_example(jl_f(:,:,imode), Ma, jl_h(:,:,imode), nfftr, nffty)
    !Call  dgemm_example(jm_f(:,:,imode), Ma1, jp_h(:,:,imode), nfftr, nffty)
    !Call  dgemm_example(jp_f(:,:,imode), Ma_1, jm_h(:,:,imode), nfftr, nffty)
    !Call  dgemm_example(rho_f(:,:,imode), Ma, rho_h(:,:,imode), nfftr, nffty)
    !Call  dgemm_example(rhoold_f(:,:,imode), Ma, rhoold_h(:,:,imode), nfftr, nffty)
  !write ( *, * ) 'Elapsed CPU time = ', t2 - t1
  !write (0,*) "el_h(10,10)", el_h(10,10,:),el_h_ptr(10,10)
  !write (0,*) "em_f"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_f(i,j,imode)
  ! end do
  !end do

  !write (0,*) "ep_f"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_f(i,j,imode)
  ! end do
  !end do
 
  !write (0,*) "ep_h", MAXVAL(abs(ep_h))

  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_h(i,j,imode)
  ! end do
  !end do

  !write (0,*) "em_h", MAXVAL(abs(em_h))
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_h(i,j,imode)
  ! end do
  !end do 
  END DO
  !el_h= el_f
  !em_h= ep_f
  !ep_h= em_f
  !bl_h= bl_f
  !bm_h= bp_f
  !bp_h= bm_f

  !write (0,*) "em_f"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_f(i,j,:)
  ! end do
  !end do

  !write (0,*) "ep_f"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_f(i,j,:)
  ! end do
  !end do

  !write (0,*) "ep_h", MAXVAL(abs(ep_h))
  !write (0,*) "el_h", MAXVAL(abs(el_h))
  !write (0,*) "bp_h", MAXVAL(abs(bp_h))
  !write (0,*) "bm_h", MAXVAL(abs(bm_h))
  !write (0,*) "bl_h", MAXVAL(abs(bl_h))
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_h(i,j,:)
  ! end do
  !end do

  !write (0,*) "em_h", MAXVAL(abs(em_h))
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_h(i,j,:)
  ! end do
  !end do
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
  !write (*,*) "====================== END OF  GET H FIELD  =============================="
END SUBROUTINE get_Hfields


SUBROUTINE get_Hfields_inv()
  USE PICSAR_precision
  USE shared_data, ONLY: nmodes
  USE fields !, ONLY:  el_h, ep_h, em_h, bl_h, bp_h, bm_h, jl_h, jp_h, jm_h, rhoold_h, rho_h
  !USE fields, ONLY : el_f, er_f, et_f, bl_f, br_f, bt_f, jl_f, jr_f, jt_f, rho_f, rhoold_h
  !USE fields, ONLY:  el_h_inv, ep_h_inv, em_h_inv, bl_h_inv, bp_h_inv, bm_h_inv 
  !                   jl_h_inv, jp_h_inv, jm_h_inv, rhoold_h_inv, rho_h_inv
  !USE fields, ONLY : invM, invM1, invM_1
  USE shared_data, ONLY:  nx,ny
  USE fields, ONLY:  nxguards, nyguards
  IMPLICIT NONE
  COMPLEX(cpx),POINTER, DIMENSION(:, :) :: el_h_ptr, ep_h_ptr, em_h_ptr, bl_h_ptr, bp_h_ptr, bm_h_ptr,&
                                  el_h_inv_ptr, ep_h_inv_ptr, em_h_inv_ptr, bl_h_inv_ptr, bp_h_inv_ptr, bm_h_inv_ptr
  !INTEGER(idp), INTENT(IN):: nfftr
  !REAL(num), dimension(:,:), allocatable :: invM 
  INTEGER (idp) :: imode, nfftr, i,j, nffty
  REAL (num) :: t1, t2
  !ALLOCATE (invM(nfftr,nfftr))
#if defined(LIBRARY)
    nfftr=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
#else
    nfftr=nx+2*nxguards
    nffty=ny+2*nyguards
#endif
  !ii = DCMPLX(0.0_NUM,1.0_NUM)
  !el_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
  !ep_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
  !em_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
  !bl_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
  !bm_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
  !bp_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
  !write (*,*) "====================== START OF  GET H FIELD INV   =============================="
  DO imode=1, nmodes
  ! ! Ma = Ma_tot(:,:,imode)
    Call Hankel_M_and_invM(imode-1)
  !  !call cpu_time ( t1 )
  !  el_h_ptr=>el_h(:,:,imode)
  !  ep_h_ptr=>ep_h(:,:,imode)
  !  em_h_ptr=>em_h(:,:,imode)
  !  bl_h_ptr=>bl_h(:,:,imode)
  !  bp_h_ptr=>bp_h(:,:,imode)
  !  bm_h_ptr=>bm_h(:,:,imode)
  !  el_h_inv_ptr=>el_h_inv(:,:,imode)
  !  ep_h_inv_ptr=>ep_h_inv(:,:,imode)
  !  em_h_inv_ptr=>em_h_inv(:,:,imode)
  !  bl_h_inv_ptr=>bl_h_inv(:,:,imode)
  !  bp_h_inv_ptr=>bp_h_inv(:,:,imode)
  !  bm_h_inv_ptr=>bm_h_inv(:,:,imode)
  !  !write (0,*) "ep_h", MAXVAL(ABS(ep_h))
  !  Call  dgemm_example(el_h_ptr, invM, el_h_inv_ptr, nfftr, nffty)
  !  Call  dgemm_example(em_h_ptr, invM_1, em_h_inv_ptr, nfftr, nffty)
  !  Call  dgemm_example(ep_h_ptr, invM1, ep_h_inv_ptr, nfftr, nffty)
  !  Call  dgemm_example(bl_h_ptr, invM, bl_h_inv_ptr, nfftr, nffty)
  !  Call  dgemm_example(bm_h_ptr, invM_1, bm_h_inv_ptr, nfftr, nffty)
  !  Call  dgemm_example(bp_h_ptr, invM1, bp_h_inv_ptr, nfftr, nffty)
 !   Call  dgemm_example(jl_h(:,:,imode), invM, jl_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(jp_h(:,:,imode), invM1, jp_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(jm_h(:,:,imode), invM_1, jm_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(rho_h(:,:,imode), invM, rho_h_inv(:,:,imode), nfftr)
 !   Call  dgemm_example(rhoold_h(:,:,imode), invM, rho_hold_inv(:,:,imode), nfftr)
    !call cpu_time ( t2 )

    Call  dgemm_example(el_h(:,:,imode), invM, el_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(em_h(:,:,imode), invM_1, em_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(ep_h(:,:,imode), invM1, ep_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bl_h(:,:,imode), invM, bl_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bm_h(:,:,imode), invM_1, bm_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bp_h(:,:,imode), invM1, bp_h_inv(:,:,imode), nfftr, nffty)
 ! write ( *, * ) 'Elapsed CPU time = ', t2 - t1

  !write (0,*) "el_h(10,10)", el_h_inv(10,10,:), "-------", el_h_inv_ptr(10,10)
  !write (0,*) "ep_h"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_h(i,j,imode)
  ! end do
  !end do
  !write (0,*) "ep_h_inv"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_h_inv(i,j,imode)
  ! end do
  !end do
  !write (0,*) "em_h"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_h(i,j,imode)
  ! end do
  !end do
  !write (0,*) "em_h_inv"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_h_inv(i,j,imode)
  ! end do
  !end do

  END DO
  !el_h_inv=el_h
  !em_h_inv=em_h
  !ep_h_inv= ep_h
  !bl_h_inv=bl_h
  !bm_h_inv=bm_h
  !bp_h_inv=bp_h


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

  !write (0,*) "ep_h"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_h(i,j,:)
  ! end do
  !end do
  !write (0,*) "ep_h_inv"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), ep_h_inv(i,j,:)
  ! end do
  !end do
  !write (0,*) "em_h"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_h(i,j,:)
  ! end do
  !end do
  !write (0,*) "em_h_inv"
  !DO i=1, nfftr
  ! DO j=1, ny
  !   write (0,*), em_h_inv(i,j,:)
  ! end do
  !end do


write (*,*) "====================== END OF  GET H FIELD INV =============================="
END SUBROUTINE get_Hfields_inv


  SUBROUTINE fft_backward_c2c_local_AM_rz(nfftx,nffty,nfftz)
    USE PICSAR_precision
    USE fastfft
    USE fields !, ONLY : el_c, er_c, et_c, bl_c, br_c, bt_c, jl_c, jr_c, jt_c,rho_c,   &
                !       rhoold_c
    !USE fields, ONLY : el_h_inv, em_h_inv, ep_h_inv, bl_h_inv, bp_h_inv, bm_h_inv
    ! jl_h_inv, jp_h_inv, jm_h_inv, rho_h_inv,   &
    !                   rhoold_h_inv

    USE fourier, ONLY: plan_rz_f_inv
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE time_stat, ONLY: timestat_itstart, localtimes
    IMPLICIT NONE
    REAL(num)   :: tmptime
    INTEGER(idp), INTENT(IN)     :: nfftx,nffty,nfftz
    COMPLEX(cpx), dimension(:,:,:), allocatable :: er_h_inv, et_h_inv, br_h_inv,bt_h_inv
    COMPLEX (cpx) :: ii
    INTEGER :: i,j
    ii = DCMPLX(0._num,1._num)

    ALLOCATE (er_h_inv(nfftx,nffty,nfftz))
    ALLOCATE (et_h_inv(nfftx,nffty,nfftz))
    ALLOCATE (br_h_inv(nfftx,nffty,nfftz))
    ALLOCATE (bt_h_inv(nfftx,nffty,nfftz))

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    
    !write (*,*) "BEFORE get_Hfields_inv  em_h_inv", MAXVAL(abs(em_h_inv)) , "ep_h_inv", MAXVAL(abs(ep_h_inv)) 
    CALL get_Hfields_inv()
    !write (*,*) "AFTER get_Hfields_inv em_h_inv",  MAXVAL(abs(em_h_inv)) , "ep_h_inv", MAXVAL(abs(ep_h_inv))
    er_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
    et_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
    br_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)
    bt_h_inv = DCMPLX(0.0_NUM, 0.0_NUM)

    er_h_inv = em_h_inv+ep_h_inv
    et_h_inv = ep_h_inv-em_h_inv
    br_h_inv = bm_h_inv+bp_h_inv
    bt_h_inv = bp_h_inv-bm_h_inv
    
    !write (*,*) "AFTER get_Hfields_inv et_h_inv"
  !DO i=1, nfftx
  ! DO j=1, nffty
  !   write (0,*), et_h_inv(i,j,:)
  ! end do
  !end do    
  !  write (*,*) "AFTER get_Hfields_inv er_h_inv"

  !DO i=1, nfftx
  ! DO j=1, nffty
  !   write (0,*), er_h_inv(i,j,:)
  ! end do
  !end do
  !el_c = DCMPLX(0.0_NUM, 0.0_NUM)
  !er_c = DCMPLX(0.0_NUM, 0.0_NUM)
  !et_c = DCMPLX(0.0_NUM, 0.0_NUM)
  !bl_c = DCMPLX(0.0_NUM, 0.0_NUM)
  !br_c = DCMPLX(0.0_NUM, 0.0_NUM)
  !bt_c = DCMPLX(0.0_NUM, 0.0_NUM)
  !el_c =el_h_inv
  !er_c=er_h_inv
  !et_c=et_h_inv
  !bl_c=bl_h_inv
  !br_c=br_h_inv
  !bt_c=-bt_h_inv
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, el_h_inv, el_c,plan_rz_f_inv)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, er_h_inv, er_c,plan_rz_f_inv)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, et_h_inv, et_c,plan_rz_f_inv)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bl_h_inv, bl_c,plan_rz_f_inv)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, br_h_inv, br_c,plan_rz_f_inv)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bt_h_inv, bt_c,plan_rz_f_inv)
   
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, el_f, el_c,plan_rz_f_inv)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, ep_f, er_c,plan_rz_f_inv)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, em_f, et_c,plan_rz_f_inv)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bl_f, bl_c,plan_rz_f_inv)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bp_f, br_c,plan_rz_f_inv)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bm_f, bt_c,plan_rz_f_inv)

    !write (*,*) "AFTER get_fields_inv  er_c", MAXVAL(ABS(er_c))
    !write (*,*) "AFTER get_fields_inv  el_c", MAXVAL(ABS(el_c))
    !write (*,*) "AFTER get_fields_inv  et_c", MAXVAL(ABS(et_c))
    !write (*,*) "AFTER get_fields_inv  br_c", MAXVAL(ABS(br_c))
    !write (*,*) "AFTER get_fields_inv  bl_c", MAXVAL(ABS(bl_c))
    !write (*,*) "AFTER get_fields_inv  bt_c", MAXVAL(ABS(bt_c))
    !DO i=1,nfftx
    !  DO j=1,nffty
    !    write (0,*) "i=", i, "j=", j, "er_c", er_c(i,j,:)
    !    write (0,*) "i=", i, "j=", j, "el_c", el_c(i,j,:)
    !    write (0,*) "i=", i, "j=", j, "et_c", et_c(i,j,:)
    !    write (0,*) "i=", i, "j=", j, "br_c", br_c(i,j,:)
    !    write (0,*) "i=", i, "j=", j, "bl_c", bl_c(i,j,:)
    !    write (0,*) "i=", i, "j=", j, "bt_c", bt_c(i,j,:)
    !  END DO
    !END DO 

    et_c = ii*et_c
    bt_c= ii*bt_c
    !write (*,*) "AFTER get_fields_inv et_c after ii ",  MAXVAL(abs(et_c)) , "er_c", MAXVAL(abs(er_c))
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

    DEALLOCATE (er_h_inv)
    DEALLOCATE (et_h_inv)
    DEALLOCATE (br_h_inv)
    DEALLOCATE (bt_h_inv)

  END SUBROUTINE fft_backward_c2c_local_AM_rz
  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes backward C2R local FFTs only in one direction the longitudinal one
  !>  - concerns only the local 
  !> pseudo-spectral solver
  !
  !> @author
  !> Imen Zemzemi
  !
  !> @date
  !> Creation 2018
  ! ______________________________________________________________________________________

  SUBROUTINE get_fields_AM_rz()
    USE fastfft
    USE fields, ONLY:  nxguards, nyguards
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY:  ny, nx, nmodes
    USE time_stat, ONLY: timestat_itstart, localtimes
    USE gpstd_solver
    IMPLICIT NONE
    REAL(num) :: tmptime
    INTEGER(idp) :: nfftx, nffty, nfftz
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
#endif
    nfftz=nmodes
    ! Get Inverse Fourier transform C2C of all fields components 
    CALL fft_backward_c2c_local_AM_rz(nfftx,nffty,nfftz)

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
#if !defined(LIBRARY)
     !CALL copy_field_backward_AM_rz
#endif
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE get_fields_AM_rz
END MODULE Hankel

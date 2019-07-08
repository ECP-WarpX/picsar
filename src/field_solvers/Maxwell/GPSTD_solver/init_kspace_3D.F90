! ________________________________________________________________________________________
!> @brief
!> This module contains subroutines that init Fourier domain arrays and variables
!> required in the spectral Maxwell solver step of the PIC cycle
!> This includes init of block matrixes for the GPSTD at t=0 as well as k-vectors
!> (local, global, or semi-global/hybrid depending on the FFT algorithm)
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________

MODULE math_tools !#do not parse 

  USE PICSAR_PRECISION
  IMPLICIT NONE
  CONTAINS 
  ! ______________________________________________________________________________________
  !> @brief
  !> This function computes SINC value of an array of real
  !
  !> @author
  !> H. Kallala
  !
  !> @params[in] block - array of REAL(num)
  !> @params[in] n1 - INTEGER(idp) - size of array block along dimension 1
  !> @params[in] n2 - INTEGER(idp) - size of array block along dimension 2
  !> @params[in] n3 - INTEGER(idp) - size of array block along dimension 3
  !> @params[out] sinc_block - array of REAL(num) - SINC of input array
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  FUNCTION sinc_block(n1, n2, n3, block)
    INTEGER(idp), INTENT(IN)                     :: n1, n2, n3
    REAL(num), DIMENSION(:, :, :), INTENT(in)  :: block
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: sinc_block
    INTEGER(idp)       :: i, j, k

    ALLOCATE(sinc_block(n1, n2, n3))
    DO k=1, n3
      DO j = 1, n2
        DO i = 1, n1
          sinc_block(i, j, k)=sinc(block(i, j, k))
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END FUNCTION sinc_block
  ! - Computes factorial of n
  FUNCTION factorial(n)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: n
    INTEGER(idp) :: factorial
    INTEGER(idp) :: i, ans

    ans = 1
    DO i = 1, n
      ans = ans * i
    END DO
    factorial = ans
  END FUNCTION factorial

  FUNCTION logfactorial(n)! returns log(n!)
    INTEGER(idp), INTENT(IN)  :: n
    REAL(num)                 :: logfactorial, x
    INTEGER(idp)              :: k

    IF(n.EQ.0_idp) THEN
      logfactorial=0.0_num
    ELSE
      x=log(1.0_num*n)
      logfactorial=x
      DO k=2, n-1
        x=log(1.0_num*k)
        logfactorial=logfactorial+x
      ENDDO
    ENDIF
    RETURN
  END FUNCTION logfactorial
  ! ______________________________________________________________________________________
  !> @brief
  !> This function computes SINC value of a REAL(num)
  !  sinc(x) = sin(x)/x if x != 0 else sinc(x) = 1.0
  !> @author
  !> H. Kallala
  !
  !> @params[in] x -  REAL(num)
  !> @params[out] sinc - REAL(num) - returns SINC of input variable
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________

  FUNCTION sinc (x)
    USE picsar_precision
    IMPLICIT NONE
    REAL(num) :: sinc
    REAL(num), INTENT(IN) ::x

    IF (x .ne. 0.0_num) THEN
      sinc=sin(x)/x
    ELSE
      sinc=1.0_num
    ENDIF
    RETURN
  END FUNCTION sinc


!  SUBROUTINE jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
!  USE picsar_precision
!  !*****************************************************************************80
!  !
!  !! JYNDD: Bessel functions Jn(x) and Yn(x), first and second
!  !derivatives.
!  !
!  !  Licensing:
!  !
!  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.
!  !    However, 
!  !    they give permission to incorporate this routine into a user
!  !    program 
!  !    provided that the copyright is acknowledged.
!  !
!  !  Modified:
!  !
!  !    02 August 2012
!  !
!  !  Author:
!  !
!  !    Shanjie Zhang, Jianming Jin
!  !
!  !  Reference:
!  !
!  !    Shanjie Zhang, Jianming Jin,
!  !    Computation of Special Functions,
!  !    Wiley, 1996,
!  !    ISBN: 0-471-11963-6,
!  !    LC: QA351.C45.
!  !
!  !  Parameters:
!  !
!  !    Input, integer ( kind = 4 ) N, the order.
!  !
!  !    Input, real ( kind = 8 ) X, the argument.
!  !
!  !    Output, real ( kind = 8 ) BJN, DJN, FJN, BYN, DYN, FYN, the values
!  !    of
!  !    Jn(x), Jn'(x), Jn"(x), Yn(x), Yn'(x), Yn"(x).
!  !
!    IMPLICIT NONE
!  
!    REAL ( num ) bj(102)
!    REAL ( num ) bjn
!    REAL ( num ) byn
!    REAL ( num ) bs
!    REAL ( num ) by(102)
!    REAL ( num ) djn
!    REAL ( num ) dyn
!    REAL ( num ) e0
!    REAL ( num ) ec
!    REAL ( num ) f
!    REAL ( num ) f0
!    REAL ( num ) f1
!    REAL ( num ) fjn
!    REAL ( num ) fyn
!    INTEGER ( idp ) k
!    INTEGER ( idp ) m
!    INTEGER ( idp ) mt
!    INTEGER ( idp ) n
!    INTEGER ( idp ) nt
!    REAL ( idp ) s1
!    REAL ( idp ) su
!    REAL ( idp ) x
!  
!    DO nt = 1, 900
!      mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt ) &
!        - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
!      IF ( 20 < mt ) THEN
!        EXIT
!      END IF 
!    END DO
!  
!    m = nt
!    bs = 0.0D+00
!    f0 = 0.0D+00
!    f1 = 1.0D-35
!    su = 0.0D+00
!    DO k = m, 0, -1
!      f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
!      IF ( k <= n + 1 ) THEN
!        bj(k+1) = f
!      END IF
!      IF ( k == 2 * int ( k / 2 ) ) THEN
!        bs = bs + 2.0D+00 * f
!        IF ( k /= 0 ) THEN
!          su = su + ( -1.0D+00 ) ** ( k / 2 ) * f / k
!        END IF
!      END IF
!      f0 = f1
!      f1 = f
!    END DO
!  
!    DO k = 0, n + 1
!      bj(k+1) = bj(k+1) / ( bs - f )
!    END DO
!  
!    bjn = bj(n+1)
!    ec = 0.5772156649015329D+00
!    e0 = 0.3183098861837907D+00
!    s1 = 2.0D+00 * e0 * ( log ( x / 2.0D+00 ) + ec ) * bj(1)
!    f0 = s1 - 8.0D+00 * e0 * su / ( bs - f )
!    f1 = ( bj(2) * f0 - 2.0D+00 * e0 / x ) / bj(1)
!  
!    by(1) = f0
!    by(2) = f1
!    DO k = 2, n + 1 
!      f = 2.0D+00 * ( k - 1.0D+00 ) * f1 / x - f0
!      by(k+1) = f
!      f0 = f1
!      f1 = f
!    END DO
!  
!    byn = by(n+1)
!    djn = - bj(n+2) + n * bj(n+1) / x
!    dyn = - by(n+2) + n * by(n+1) / x
!    fjn = ( n * n / ( x * x ) - 1.0D+00 ) * bjn - djn / x
!    fyn = ( n * n / ( x * x ) - 1.0D+00 ) * byn - dyn / x
!  
!    RETURN
!  END SUBROUTINE jyndd
!  
!  
!  
!  SUBROUTINE jyzo ( n, nt, rj0 )
!  USE picsar_precision 
!  !*****************************************************************************80
!  !
!  !! JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
!  !
!  !  Licensing:
!  !
!  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!  !    they give permission to incorporate this routine into a user program 
!  !    provided that the copyright is acknowledged.
!  !
!  !  Modified:
!  !
!  !    28 July 2012
!  !
!  !  Author:
!  !
!  !    Shanjie Zhang, Jianming Jin
!  !
!  !  Reference:
!  !
!  !    Shanjie Zhang, Jianming Jin,
!  !    Computation of Special Functions,
!  !    Wiley, 1996,
!  !    ISBN: 0-471-11963-6,
!  !    LC: QA351.C45.
!  !
!  !  Parameters:
!  !
!  !    Input, integer ( kind = 4 ) N, the order of the Bessel functions.
!  !
!  !    Input, integer ( kind = 4 ) NT, the number of zeros.
!  !
!  !    Output, real ( kind = 8 ) RJ0(NT), RJ1(NT), RY0(NT), RY1(NT), the zeros 
!  !    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!  !
!    IMPLICIT NONE
!  
!    INTEGER ( idp ) nt
!  
!    REAL ( num ) bjn
!    REAL ( num ) byn
!    REAL ( num ) djn
!    REAL ( num ) dyn
!    REAL ( num ) fjn
!    REAL ( num ) fyn
!    INTEGER ( idp ) l
!    INTEGER ( idp ) n
!    REAL ( num ) n_r8
!    REAL ( num ) rj0(nt)
!    REAL ( num ) x
!    REAL ( num ) x0
!  
!    n_r8 = REAL ( n, num )
!  
!    IF ( n <= 20 ) THEN
!      x = 2.82141D+00 + 1.15859D+00 * n_r8 
!    ELSE
!      x = n + 1.85576D+00 * n_r8 ** 0.33333D+00 &
!        + 1.03315D+00 / n_r8 ** 0.33333D+00
!    END IF
!  
!    l = 0
!  
!    DO
!  
!      x0 = x
!      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
!      x = x - bjn / djn
!  
!      IF ( 1.0D-09 < abs ( x - x0 ) ) THEN
!        CYCLE
!      END IF
!  
!      l = l + 1
!      rj0(l) = x
!      x = x + 3.1416D+00 + ( 0.0972D+00 + 0.0679D+00 * n_r8 &
!        - 0.000354D+00 * n_r8 ** 2 ) / l
!  
!      IF ( nt <= l ) THEN
!        EXIT
!      END IF
!  
!    END DO 
!  
!    RETURN
!  END SUBROUTINE jyzo 

       SUBROUTINE JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
!
!       ===========================================================
!       Purpose: Compute Bessel functions Jn(x) and Yn(x), and
!                their first and second derivatives
!       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
!                n   ---  Order of Jn(x) and Yn(x)
!       Output:  BJN ---  Jn(x)
!                DJN ---  Jn'(x)
!                FJN ---  Jn"(x)
!                BYN ---  Yn(x)
!                DYN ---  Yn'(x)
!                FYN ---  Yn"(x)
!       Routines called:
!                JYNBH to compute Jn and Yn
!       ===========================================================
!       
        USE PICSAR_precision
        IMPLICIT NONE
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER (IDP), INTENT(IN) :: N
        REAL (NUM) , INTENT(IN) :: X
        REAL (NUM) , INTENT (OUT) :: BJN,DJN,FJN,BYN,DYN,FYN
        REAL(NUM), DIMENSION(:) :: BJ(2),BY(2)
        integer (idp):: NM
        CALL JYNBH(N+1,N,X,NM,BJ,BY)
!       Compute derivatives by differentiation formulas
        BJN=BJ(1)
        BYN=BY(1)
        DJN=-BJ(2)+N*BJ(1)/X
        DYN=-BY(2)+N*BY(1)/X
        FJN=(N*N/(X*X)-1.0D0)*BJN-DJN/X
        FYN=(N*N/(X*X)-1.0D0)*BYN-DYN/X
        RETURN
        END

        SUBROUTINE JYNBH(N,NMIN,X,NM,BJ,BY)
!
!       =====================================================
!       Purpose: Compute Bessel functions Jn(x), Yn(x)
!       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
!                n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
!                nmin -- Lowest order computed  ( nmin ≥ 0 )
!       Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
!                BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =====================================================
!
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        USE PICSAR_precision
        IMPLICIT NONE
        INTEGER (IDP), INTENT(IN) :: N, NMIN
        INTEGER (IDP), INTENT(OUT) :: NM
        REAL (NUM), INTENT(IN) :: X
        REAL (NUM), DIMENSION(:), intent(out) :: BJ(0:N-NMIN),BY(0:N-NMIN)
        REAL (NUM), DIMENSION(:) ::A(4),B(4),A1(4),B1(4)
        REAL (NUM):: PI=3.141592653589793D0
        REAL (NUM):: R2P=.63661977236758D0
        REAL (NUM):: BS,SU,SV,F2,F1,F, S0, BJ0, BJ1, BJK, BY0, BY1, BYK, EC, T1, T2, P0,Q0,P1, Q1, CU
        INTEGER (IDP):: K, M, KY
        NM=N
        IF (X.LT.1.0D-100) THEN
           DO 10 K=NMIN,N
              BJ(K-NMIN)=0.0D0
10            BY(K-NMIN)=-1.0D+300
           IF (NMIN.EQ.0) BJ(0)=1.0D0
           RETURN
        ENDIF
        IF (X.LE.300.0.OR.N.GT.INT(0.9*X, IDP)) THEN
!          Backward recurrence for Jn
           IF (N.EQ.0) NM=1
           M=MSTA1(X,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(X,NM,15)
           ENDIF
           BS=0.0D0
           SU=0.0D0
           SV=0.0D0
           F2=0.0D0
           F1=1.0D-100
           F=0.0D0
           DO 15 K=M,0,-1
              F=2.0D0*(K+1.0D0)/X*F1-F2
              IF (K.LE.NM .AND. K.GE.NMIN) THEN
                BJ(K-NMIN)=F
              END IF
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 BS=BS+2.0D0*F
                 SU=SU+(-1)**(K/2)*F/K
              ELSE IF (K.GT.1) THEN
                 SV=SV+(-1)**(K/2)*K/(K*K-1.0D0)*F
              ENDIF
              F2=F1
15            F1=F
           S0=BS+F
           DO 20 K=NMIN,NM
20            BJ(K-NMIN)=BJ(K-NMIN)/S0
!          Estimates for Yn at start of recurrence
           BJ0 = F1 / S0
           BJ1 = F2 / S0
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           BY0=R2P*(EC*BJ0-4.0D0*SU/S0)
           BY1=R2P*((EC-1.0D0)*BJ1-BJ0/X-4.0D0*SV/S0)
           IF (0.GE.NMIN) BY(0-NMIN)=BY0
           IF (1.GE.NMIN) BY(1-NMIN)=BY1
           KY=2
        ELSE
!          Hankel expansion
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,  &
                  -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00, &
                  .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00, &
                  .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00, &
                  -.1993531733751297D+01,.2724882731126854D+02/
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO 25 K=1,4
              P0=P0+A(K)*X**(-2*K)
25            Q0=Q0+B(K)*X**(-2*K-1)
           CU=DSQRT(R2P/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           IF (0.GE.NMIN) BJ(0-NMIN)=BJ0
           IF (0.GE.NMIN) BY(0-NMIN)=BY0
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO 30 K=1,4
              P1=P1+A1(K)*X**(-2*K)
30            Q1=Q1+B1(K)*X**(-2*K-1)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
           IF (1.GE.NMIN) BJ(1-NMIN)=BJ1
           IF (1.GE.NMIN) BY(1-NMIN)=BY1
           DO 35 K=2,NM
              BJK=2.0D0*(K-1.0D0)/X*BJ1-BJ0
              IF (K.GE.NMIN) BJ(K-NMIN)=BJK
              BJ0=BJ1
35            BJ1=BJK
           KY=2
        ENDIF
!       Forward recurrence for Yn
        DO 45 K=KY,NM
           BYK=2.0D0*(K-1.0D0)*BY1/X-BY0
           IF (K.GE.NMIN) BY(K-NMIN)=BYK
           BY0=BY1
45         BY1=BYK
        RETURN
        END
        INTEGER FUNCTION MSTA1(X,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================
!
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        USE PICSAR_PRECISION
        IMPLICIT NONE
        REAL (NUM), INTENT (IN) :: X
        INTEGER (ISP), INTENT (IN) :: MP
        INTEGER (IDP) :: N0, IT, NN, N1
        REAL (NUM) :: A0,F0,F, F1
        A0=DABS(X)
        N0=INT(1.1D0*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END
       INTEGER FUNCTION MSTA2(X,N,MP)
!
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        USE PICSAR_PRECISION
        IMPLICIT NONE
        REAL (NUM), INTENT (IN) :: X
        INTEGER (ISP), INTENT (IN):: MP
        INTEGER (IDP), INTENT (IN) :: N
        REAL (NUM) :: A0, HMP, EJN, OBJ,F0, F1, F
        INTEGER (IDP) ::  N0, IT, NN, N1
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)+1
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL FUNCTION ENVJ(N,X)
        USE PICSAR_PRECISION
        REAL (NUM), INTENT (IN):: X
        INTEGER (IDP) , INTENT (IN) :: N
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
        SUBROUTINE JYZO(N,NT,RJ0,RJ1,RY0,RY1)
!
!       ======================================================
!       Purpose: Compute the zeros of Bessel functions Jn(x),
!                Yn(x), and their derivatives
!       Input :  n  --- Order of Bessel functions  (n >= 0)
!                NT --- Number of zeros (roots)
!       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
!                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
!                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
!                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
!       Routine called: JYNDD for computing Jn(x), Yn(x), and
!                       their first and second derivatives
!       ======================================================
!
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        USE PICSAR_PRECISION
        IMPLICIT NONE
        INTEGER (IDP), INTENT (IN):: N, NT
        REAL (NUM), DIMENSION(:), INTENT(OUT):: RJ0(NT),RJ1(NT),RY0(NT),RY1(NT)
        INTEGER (IDP):: L
        REAL (NUM):: X, XGUESS, X0, PI, BJN, BYN, DJN, DYN,FJN, FYN
        PI=3.141592653589793D0
!       -- Newton method for j_{N,L}
!       1) initial guess for j_{N,1}
        IF (N.LE.20) THEN
           X=2.82141+1.15859*N
        ELSE
!          Abr & Stg (9.5.14)
           X=N+1.85576*N**0.33333+1.03315/N**0.33333
        ENDIF
        L=0
!       2) iterate
        XGUESS=X
10      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-BJN/DJN
        IF (X-X0.LT.-1) X=X0-1
        IF (X-X0.GT.1) X=X0+1
        IF (DABS(X-X0).GT.1.0D-11) GO TO 10
!       3) initial guess for j_{N,L+1}
        IF (L.GE.1)THEN
           IF (X.LE.RJ0(L)+0.5) THEN
              X=XGUESS+PI
              XGUESS=X
              GO TO 10
           ENDIF
        END IF
        L=L+1
        RJ0(L)=X
!       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.0972d0+0.0679*N-0.000354*N**2)/L, 0d0)
        IF (L.LT.NT) GO TO 10
!       -- Newton method for j_{N,L}'
        IF (N.LE.20) THEN
           X=0.961587+1.07703*N
        ELSE
           X=N+0.80861*N**0.33333+0.07249/N**0.33333
        ENDIF
        IF (N.EQ.0) X=3.8317
        L=0
        XGUESS=X
15      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-DJN/FJN
        IF (X-X0.LT.-1) X=X0-1
        IF (X-X0.GT.1) X=X0+1
        IF (DABS(X-X0).GT.1.0D-11) GO TO 15
        IF (L.GE.1)THEN
           IF (X.LE.RJ1(L)+0.5) THEN
              X=XGUESS+PI
              XGUESS=X
              GO TO 15
           ENDIF
        END IF
        L=L+1
        RJ1(L)=X
!       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.4955d0+0.0915*N-0.000435*N**2)/L, 0d0)
        IF (L.LT.NT) GO TO 15
!       -- Newton method for y_{N,L}
        IF (N.LE.20) THEN
           X=1.19477+1.08933*N
        ELSE
           X=N+0.93158*N**0.33333+0.26035/N**0.33333
        ENDIF
        L=0
        XGUESS=X
20      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-BYN/DYN
        IF (X-X0.LT.-1) X=X0-1
        IF (X-X0.GT.1) X=X0+1
        IF (DABS(X-X0).GT.1.0D-11) GO TO 20
        IF (L.GE.1)THEN
           IF (X.LE.RY0(L)+0.5) THEN
              X=XGUESS+PI
              XGUESS=X
              GO TO 20
           END IF
        END IF
        L=L+1
        RY0(L)=X
!       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.312d0+0.0852*N-0.000403*N**2)/L,0d0)
        IF (L.LT.NT) GO TO 20
!       -- Newton method for y_{N,L}'
        IF (N.LE.20) THEN
           X=2.67257+1.16099*N
        ELSE
           X=N+1.8211*N**0.33333+0.94001/N**0.33333
        ENDIF
        L=0
        XGUESS=X
25      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-DYN/FYN
        IF (DABS(X-X0).GT.1.0D-11) GO TO 25
        IF (L.GE.1) THEN
           IF (X.LE.RY1(L)+0.5) THEN
              X=XGUESS+PI
              XGUESS=X
              GO TO 25
           END IF
        END IF
        L=L+1
        RY1(L)=X
!       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.197d0+0.0643*N-0.000286*N**2)/L,0d0)
        IF (L.LT.NT) GO TO 25
        RETURN
        END




END MODULE math_tools


MODULE gpstd_solver
  USE math_tools
  USE PICSAR_PRECISION
  IMPLICIT NONE
  COMPLEX(cpx), DIMENSION(:), ALLOCATABLE :: kxc, kxb, kxf, kyc, kyb, kyf, kzc, kzb,  &
  kzf
  COMPLEX(cpx), DIMENSION(:,:), ALLOCATABLE :: krc

  ! Flattened array of matrix_blocks indices that are usefull for the PSATD push 
  ! in the case of periodic boundary conditions
  ! if is_usefull_per(i) == 1 then the cc_mat(1)%matrix_blocks(k,j) is allocated 
  ! and corresponds to a required block for the PSATD computation
  ! else if is_usefull_per(i) == 0 then cc_mat(1)%matrix_blocks(k,j)  is 1x1x1 null block
  ! with i = (k-1) * 11 + (j-1) + 1
 ! INTEGER(isp)  :: is_usefull_per(121) = &
 !    (/1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1 ,&
 !      0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0 ,&
 !      1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0 ,&
 !      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,&
 !      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,&
 !      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  INTEGER(idp) :: is_usefull_per(33) = &
     (/ 0,  4,  5,  6,  9, 10, 12, 14, 16, 18, 20, 21, 24,               & 
       25, 26, 30, 31, 32, 34, 35, 36, 40, 41, 44, 46, 48,               & 
       50, 52, 55, 56, 60, 61, 62/)

  ! Flattened array of matrix_blocks indices that are usefull for the PSATD push 
  ! in the case of absorbing boundary conditions
  ! if is_usefull_abs(i) == 1 then the cc_mat(1)%matrix_blocks(k,j) is allocated 
  ! and corresponds to a required block for the PSATD computation
  ! else if is_usefull_abs(i) == 0 then cc_mat(1)%matrix_blocks(k,j)  is 1x1x1 null block
  ! with i = (k-1) * 17 + (j-1) + 1
  !INTEGER(isp) :: is_usefull_abs(289) = &
  !   (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, &
  !      0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, &
  !      1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, &
  !      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, &
  !      0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, &
  !      0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, &
  !      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, &
  !      1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, &
  !      0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  !      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  !      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  !      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  !      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
  !      0, 0, 0/)
  INTEGER(idp) :: is_usefull_abs(51) = &
        (/0,  10,  11,  12,  15,  16,  18,  25,  26,  36,  44,            &
        45,  47,  49,  50,  54,  57,  58,  72,  76,  77,  82,             &
        83,  84,  90,  91,  92, 106, 107, 108, 115, 116, 121,             &
        122, 126, 140, 141, 144, 148, 150, 153, 154, 162, 172,            &
        173, 180, 182, 183, 187, 188, 198/)

  CONTAINS
#if defined(FFTW)
  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes dimensions on which the FFT is performed depending on
  !> the type of FFT performed (local FFT or distributed FFT but without MPI groups).
  !> It includes 2D and 3D dimensions.
  !> N.B: this routine is deprecated and will be removed in upcoming versions.
  !
  !> IMPORTANT : this routine has been adjusted to the spectral azimuthal RZ
  !> without major changes to avoid complexity and duplication of routines
  !> Here: nx, nxguards are respectively nr and nrguards 
  !>       ny, nyguards are respectively nl and nlguards
  !>       nfftz is nmodes as actually fields are 2D but we chose a 3D structure
  !>       for better performance  
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  !
  !> @params[in,out] nfftx INTEGER(idp) - Number of points in spectral space along X
  !> @params[in,out] nffty INTEGER(idp) - Number of points in spectral space along Y
  !> @params[in,out] nfftz INTEGER(idp) - Number of points in spectral space along Z
  ! ______________________________________________________________________________________
  SUBROUTINE select_case_dims_local(nfftx, nffty, nfftz)
    USE fields, ONLY: nzguards, nxguards, nyguards, l_AM_rz
    USE group_parameters, ONLY: nx_group, ny_group, nz_group, p3d_fsize
    USE iso_c_binding
    USE mpi_fftw3, ONLY: local_nz, local_nz_tr
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: c_dim, fftw_hybrid, fftw_mpi_transpose, fftw_with_mpi,    &
      nx, nx_global, ny, ny_global, nz, nz_global, p3dfft_flag, nmodes
    INTEGER(idp), INTENT(INOUT) :: nfftx, nffty, nfftz

    IF(fftw_with_mpi) THEN
      !> Global pseudo spectral solver, Only periodic bcs, no guard cells
      IF(.NOT. fftw_hybrid) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx=nx_global
          nffty=ny_global
          nfftz=local_nz
        ELSE IF(fftw_mpi_transpose) THEN
          nfftx=nx_global
          nffty=nz_global
          nfftz=local_nz_tr
        ENDIF
      !> Hybrid pseudo spectral solver
      ELSE IF(fftw_hybrid) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx = nx_group
          nffty = ny_group
          nfftz = local_nz
        ELSE IF(fftw_mpi_transpose) THEN
          nfftx = nx_group
          nffty = nz_group
          nfftz = local_nz_tr
        ENDIF
      ENDIF
#if defined(P3DFFT)
      !> If p3d_fft to perform ffts and create groups in y and z directions
      IF(p3dfft_flag) THEN
        nfftx = p3d_fsize(1)
        nffty = p3d_fsize(2)
        nfftz = p3d_fsize(3)
      ENDIF
#endif
    !> Local pseudo spectral solver
    ELSE IF(.NOT. fftw_with_mpi) THEN
      !> When using picsar with smilei
#if defined(LIBRARY)
      nfftx = nx+2*nxguards+1
      nffty = ny+2*nyguards+1
      IF (l_AM_rz) THEN
        nfftz = nmodes
      ELSE
        nfftz = nz+2*nzguards+1
      END IF
      !> if l_AM_rz is not precised than it's the cartesian geometry 
      !> that is true by default otherwise nfftz is the number of modes 
#else
      !> When using picsar
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      IF (l_AM_rz) THEN
        nfftz = nmodes
      ELSE
        nfftz = nz+2*nzguards
      END IF 
      !> if l_AM_rz is not precised than it's the cartesian geometry 
      !> that is true by default otherwise nfftz is the number of modes
#endif
      IF(c_dim ==2) THEN
        nffty = 1
      ENDIF
    ENDIF
  END SUBROUTINE select_case_dims_local

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes dimensions on which the FFT is performed depending on
  !> the type of FFT performed (local FFT, distributed FFT with and without MPI groups).
  !> N.B: this subroutine will fully replace deprecated select_case_dims_local subroutine
  !> in future release of PICSAR.
  !
  !> IMPORTANT : this routine has been adjusted to the spectral azimuthal RZ
  !> without major changes to avoid complexity and duplication of routines
  !> Here: nx, nxguards are respectively nr and nrguards 
  !>       ny, nyguards are respectively nl and nlguards
  !>       nfftz is nmodes as actually fields are 2D but we chose a 3D structure
  !>       for better performance  
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  !
  !> @params[in,out] nfftx INTEGER(idp) - Number of points in spectral space along X
  !> @params[in,out] nffty INTEGER(idp) - Number of points in spectral space along Y
  !> @params[in,out] nfftz INTEGER(idp) - Number of points in spectral space along Z
  ! ______________________________________________________________________________________
  SUBROUTINE select_case_dims_global(nfftx, nffty, nfftz)
    USE fields, ONLY: nzguards, nxguards, nyguards, l_AM_rz
    USE group_parameters, ONLY: nx_group, ny_group, nz_group
    USE iso_c_binding
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: c_dim, fftw_hybrid, fftw_mpi_transpose, fftw_with_mpi,    &
      nx, nx_global, ny, ny_global, nz, nz_global, p3dfft_flag, p3dfft_stride, nmodes
    INTEGER(idp), INTENT(INOUT) :: nfftx, nffty, nfftz
    !> When using global or hybrid pseudo spectral solver
    IF( fftw_with_mpi) THEN
      !> When using global pseudo spectral solver with periodic bcs and no guard
      !> cells
      IF(.NOT. fftw_hybrid) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx=nx_global
          nffty=ny_global
          nfftz=nz_global
        ELSE IF(fftw_mpi_transpose) THEN
          nfftx=nx_global
          nffty=nz_global
          nfftz=ny_global
        ENDIF
      ELSE IF(fftw_hybrid) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx = nx_group
          nffty = ny_group
          nfftz = nz_group
        ELSE IF(fftw_mpi_transpose) THEN
          nfftx = nx_group
          nffty = nz_group
          nfftz = ny_group
        ENDIF
      ENDIF
    !> local pseudo spectral solver
    ELSE IF(.NOT. fftw_with_mpi) THEN
#if defined(LIBRARY)
      !> When using Smilei with picsar
      nfftx = nx+2*nxguards+1
      nffty = ny+2*nyguards+1
      IF (l_AM_rz) THEN
        nfftz = nmodes
      ELSE
        nfftz = nz+2*nzguards+1
      END IF 
      !> if l_AM_rz is not precised than it's the cartesian geometry 
      !> that is true by default otherwise nfftz is the number of modes

#else
      !> When using only picsar
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      IF (l_AM_rz) THEN
        nfftz = nmodes
      ELSE
        nfftz = nz+2*nzguards
      !> if l_AM_rz is not precised than it's the cartesian geometry 
      !> that is true by default otherwise nfftz is the number of modes
      ENDIF
#endif
    ENDIF
    IF(c_dim ==2) THEN
      nffty=1
    ENDIF

#if defined(P3DFFT)
    !> When using P3DFFT library to perform ffts
    IF(p3dfft_flag) THEN
      IF(p3dfft_stride) THEN
       nfftx =nz_group
       nffty =ny_group
       nfftz =nx_group
      !> When P3D is compiled with stride flag on
      ELSE IF( .NOT. p3dfft_stride) THEN
       nfftx = nx_group
       nffty = ny_group
       nfftz = nz_group
     ENDIF
    ENDIF
#endif
  END SUBROUTINE select_case_dims_global


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
  USE PICSAR_precision
  USE math_tools
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
  !#if defined (DEBUG=1)
  !Do i= 1, nfftr
  !  Do j=1,nfftr
  !    if  (invM(i,j) .eq. 0.) then 
  !    write (0,*) "i,j=", i, j
  !    endif
  !  end do
  !end do
  !#endif
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
        !write(0,*), "nfftr=", nfftr, "imode=", imode, "p=", p, "rmax=", rmax
        CALL hankel_matrix_init(nfftr,imode,p,rmax,invM_1)
        Ma_1 =inv(invM_1)
        !write(*,*) "case -1"
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
        !write (*,*) "case 0"
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
        !write (*,*) "case 1"
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



END SUBROUTINE Hankel_M_and_invM

  !Ma=inv1(invM)
  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes block matrixes of fourier space wavelength vector
  !> As well as different other blocks usefull to compute block matrixes for psatd
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE init_kspace
    USE constants, ONLY: clight
    USE fields, ONLY: l_staggered, norderx, nordery, norderz, nxguards, nyguards,    &
      nzguards
    USE iso_c_binding
    USE matrix_coefficients, ONLY: at_op, kspace
    USE matrix_data, ONLY: nmatrixes2, ns_max
    USE omp_lib
    USE fields, ONLY : nxguards, nyguards, nzguards, l_staggered, l_AM_rz
    USE fields, ONLY : norderx, nordery, norderz
    USE params, ONLY: dt
    USE picsar_precision, ONLY: idp, num, lp, cpx
    USE shared_data, ONLY: c_dim, fftw_mpi_transpose, p3dfft_flag, p3dfft_stride

    REAL(num), ALLOCATABLE, DIMENSION(:, :, :)    :: temp, temp2
    INTEGER(idp)                                  :: i, j, k
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: nfftx, nffty, nfftz, nfftxr
    LOGICAL(lp)                                   :: switch
    LOGICAL (lp), DIMENSION(:), ALLOCATABLE :: switch_rz
    !> kspace is the (inf)finite order wave vector block matrix in 
    !> different directions
    !> It has 10 blocks: 
    !> i=3,6,9 => kspace(nmatrixes2)%block_vector(i) = 
    !> centered  derivative  operator
    !> along x,y,z directions respectively.
    !> i=1,4,7 => kspace(nmatrixes2)%block_vector(i) =
    !> derivative operator from dual to primal meshgrid (involved in Maxwell
    !> Ampere equation) along x,y,z respectively
    !> i=2,5,8 => kspace(nmatrixes2)%block_vector(i) =
    !> derivative operator from primal to dual meshgrid (involved in Maxwell
    !> Faraday equation) along x,y,z respectively
    !> kspace(nmatrixes2)%block_vector(10) = Absolute value of wave vector


    !> PS: If using fftw_mpi_transpose or strided p3dfft then
    !> y and z axis (or x and z for strided p3dfft) are transposed inside kspace
    !> blocks. 
    !> But the convention above remains  identical
  
    !> at_op is the block matrix for different operators involved
    !> in psatd block matrixes computations.
    !> It has 4 blocks
  
    !> The first component (for null frequency) of each block
    !> except block 2 is computed using Taylor expansion 

    !> at_op(nmatrixes2)%block_vector(1) = sin(|K|*c*dt)/(|K|)
    !> at_op(nmatrixes2)%block_vector(2) = cos(|K|*c*dt)
    !> at_op(nmatrixes2)%block_vector(3) = (1-cos(|K|*c*dt))/|K|**2
    !> at_op(nmatrixes2)%block_vector(4) = (sin(|K|*c*dt)/|K|-c*dt)/|K|**2

    nmatrixes2=nmatrixes2+1

    IF(.NOT. ASSOCIATED(kspace)) THEN
      ALLOCATE(kspace(ns_max))
    ENDIF
    ALLOCATE(kspace(nmatrixes2)%block_vector(10_idp))

    IF(.NOT. ASSOCIATED(at_op)) THEN
      ALLOCATE(at_op(ns_max))
    ENDIF
    ALLOCATE(at_op(nmatrixes2)%block_vector(4_idp))

    CALL select_case_dims_local(nfftx, nffty, nfftz)
    !> if l_AM_rz is true then nfftxr=nfftx because it's complex to complex
    !> transform in this case 
    IF (l_AM_rz) THEN 
      nfftxr = nfftx
    ELSE 
      nfftxr = nfftx/2+1
    ENDIF 
    IF(p3dfft_flag) nfftxr = nfftx
    DO i = 1_idp, 10_idp
      ALLOCATE(kspace(nmatrixes2)%block_vector(i)%block3dc(nfftxr, nffty, nfftz))
    ENDDO
    DO i = 1_idp, 4_idp
      ALLOCATE(at_op(nmatrixes2)%block_vector(i)%block3dc(nfftxr, nffty, nfftz))
    ENDDO
    !construct kspace
    ii=DCMPLX(0.0_num, 1.0_num)
    !> computes wave vector for a staggered or an unstaggered grid 
    !> takes into account norderx, nordery, norderz
    !> if norder == 0 then compute wave vector for an infinite order stencil
    CALL compute_k_vec(l_staggered)
    !> the loop for spectral RZ is interpreted as follows:
    !> nfftz is the number of azimuthal modes in theta == nmodes
    !> nffty is the number of point in the spectral grid along l direction
    !> nfftxr is the number of points in the spectral grid along r direction
    !!!!!!! Those values are computed in the select_case_dims_local routine
    DO k = 1, nfftz
      DO j = 1, nffty
        DO i = 1, nfftxr
          !> Here we separate the case of spectral RZ from the cartesian one
          !> we only have two k one in the longitudinal component
          !> the other one is in the radial direction noted in the article k_|_ 
          !> k_|_  depends on the mode which is passed by nfftz here
          IF (l_AM_rz) THEN
            kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kyf(j)
            kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kyb(j)
            kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kyc(j)
            kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = krc(i,k)
            kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = krc(i,k)
            kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = krc(i,k)
            kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = krc(i,k)
            kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = krc(i,k)
            kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = krc(i,k)
            !write (0,*) "i=", i, "j=", j, "k=", k,  "krc=", krc(i,k), "kyc=", kyc(j)
            
          ELSE  
            IF(.NOT. p3dfft_flag) THEN
              IF(.NOT. fftw_mpi_transpose) THEN
                kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
                kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
                kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
                IF(c_dim == 3) THEN
                  !IF (l_AM_rz) THEN
                  !  kspace(nmatrixes2)%block_vector(4)%block3dc(i,j, k) = krc(j,k)
                  !  kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = krc(j,k)
                  !  kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = krc(j,k)
                  !ELSE
                    kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                    kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                    kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
                  !ENDIF 
                ELSE IF(c_dim == 2) THEN
                  !> If c_dim == 2 Then y derivative is null
                  !> c_dim = 2 cannot be used with p3dfft or fftw_mpi_transpose
                  !> flag
                  kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k)= (0.0_num,0.0_num) 
                  kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k)= (0.0_num,0.0_num)
                  kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k)= (0.0_num,0.0_num) 
                ENDIF
                !> If we do spectral with azimutal modes in cylindrical then
                !> ky=kz=kr 
                !IF (l_AM_rz) THEN
                !  kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = krc(j,k)
                !  kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = krc(j,k)
                !  kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = krc(j,k)
                !ENDIF         
                kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kzf(k)
                kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kzb(k)
                kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kzc(k)
              ELSE IF(fftw_mpi_transpose) THEN
                !> If fftw_mpi_transpose kyc is the derivative operator along z and 
                !> kzc is the derivative  operator along y
                kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
                kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
                kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
                kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kzf(k)
                kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kzb(k)
                kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kzc(k)
                kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kyf(j)
                kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kyb(j)
                kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kyc(j)
              ENDIF
            ELSE IF(p3dfft_flag) THEN
              IF(p3dfft_stride) THEN
                  !> If p3dfft_stride x and z axis are transposed: 
                  !> kzc is the derivative along x and kxc is the derivative
                  !> along z
                  kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kzf(k)
                  kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kzb(k)
                  kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kzc(k)
                  kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                  kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                  kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
                  kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kxf(i)
                  kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kxb(i)
                  kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kxc(i)
              ELSE IF(.NOT. p3dfft_stride) THEN
                  kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
                  kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
                  kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
                  kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                  kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                  kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
                  kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kzf(k)
                  kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kzb(k)
                  kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kzc(k)
              ENDIF
            ENDIF
          ENDIF  
        ENDDO
      ENDDO
    ENDDO

    !> Computes the norm of wave vector in fourier space 

    IF (l_AM_rz) THEN 
      kspace(nmatrixes2)%block_vector(10)%block3dc=                                  &
      SQRT(ABS(kspace(nmatrixes2)%block_vector(1)%block3dc)**2 +                     &
        ABS(kspace(nmatrixes2)%block_vector(4)%block3dc)**2) 
    ELSE 
      kspace(nmatrixes2)%block_vector(10)%block3dc=                                    &
      SQRT(ABS(kspace(nmatrixes2)%block_vector(9)%block3dc)**2 +                       &
          ABS(kspace(nmatrixes2)%block_vector(6)%block3dc)**2 +                        &
          ABS(kspace(nmatrixes2)%block_vector(3)%block3dc)**2)

    END IF
    switch = .FALSE.
    !write (0,*), "max of absolute value of k", MAXVAL(ABS(kspace(nmatrixes2)%block_vector(10)%block3dc))
    ALLOCATE(temp(nfftxr, nffty, nfftz))
    ALLOCATE(temp2(nfftxr, nffty, nfftz))

    temp=dt*clight*REAL(kspace(nmatrixes2)%block_vector(10)%block3dc, num)
    temp2=sinc_block(nfftxr, nffty, nfftz, temp)

    at_op(nmatrixes2)%block_vector(1)%block3dc = DCMPLX(temp2, 0._num)

    !write (0,*) "block_vector(1) ", MAXVAL(ABS(at_op(nmatrixes2)%block_vector(1)%block3dc))
    at_op(nmatrixes2)%block_vector(1)%block3dc =                                      &
    clight*dt*at_op(nmatrixes2)%block_vector(1)%block3dc
    temp2=COS(temp)

    at_op(nmatrixes2)%block_vector(2)%block3dc = DCMPLX(temp2, 0._num)
    temp=0.5_num*temp

    !write (0,*) "block_vector(2) ", MAXVAL(ABS(at_op(nmatrixes2)%block_vector(2)%block3dc))
    !at_op(nmatrixes2)%block_vector(3)%block3dc = 2._num*(clight*dt/2.0_num)**2        &
    !*sinc_block(nfftxr, nffty, nfftz, temp)*sinc_block(nfftxr, nffty, nfftz,    &
    !temp)
    ALLOCATE (switch_rz(nfftz))
    switch_rz= .FALSE.	
    !> if current mpi task contains the null frequency then this processor it
    !> tagged by switch = .TRUE. in order perform Taylor expansion
    !> for at_op...(i)(1,1,1) only in this mpi task
    IF (l_AM_rz) THEN
      DO i=1, nfftz
        IF(ABS(kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, i)) .EQ. 0.0_num) THEN
          kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, i) = DCMPLX(1.0_num,         &
          0.0_num)
          switch_rz(i) = .TRUE.
        ENDIF
      END DO
    ELSE 
      IF(ABS(kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)) .EQ. 0.0_num) THEN
        kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = DCMPLX(1.0_num,         &
        0.0_num)
        switch = .TRUE.
      ENDIF
     END IF 


    at_op(nmatrixes2)%block_vector(3)%block3dc = (DCMPLX(1.0_num, 0.0_num) -          &
    at_op(nmatrixes2)%block_vector(2)%block3dc)                                       &
    /kspace(nmatrixes2)%block_vector(10)%block3dc**2
    !(1-C)/k^2
    at_op(nmatrixes2)%block_vector(4)%block3dc =                                      &
    (at_op(nmatrixes2)%block_vector(1)%block3dc-clight*dt) /                          &
    kspace(nmatrixes2)%block_vector(10)%block3dc/                                     &
    kspace(nmatrixes2)%block_vector(10)%block3dc

    !> Performs Taylor expansion for
    !> at_op(nmatrixes2)%block_vector(3-4)%block3dc(1, 1, 1)
    !> Taylor expansion for
    !> at_op(nmatrixes2)%block_vector(1)%block3dc(1, 1, 1) is performed inside
    !> sinc function

    IF (l_AM_rz) THEN
      DO i=1, nfftz
        IF (switch_rz(i)) THEN
          at_op(nmatrixes2)%block_vector(3)%block3dc(1, 1, i) = (clight*dt)**2/2.0_num
          at_op(nmatrixes2)%block_vector(4)%block3dc(1, 1,                                &
          i)=DCMPLX(-(clight*dt)**3/6.0_num, 0.0_num)
          kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, i)=DCMPLX(0._num, 0._num)
        END IF
      END DO
    ELSE
      IF (switch) THEN
        at_op(nmatrixes2)%block_vector(3)%block3dc(1, 1, 1) = (clight*dt)**2/2.0_num
        at_op(nmatrixes2)%block_vector(4)%block3dc(1, 1,                                &
        1)=DCMPLX(-(clight*dt)**3/6.0_num, 0.0_num)
        kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)=DCMPLX(0._num, 0._num)
      END IF 
    END IF 

    !write (0,*) "block_vector(3) ", MAXVAL(ABS(at_op(nmatrixes2)%block_vector(3)%block3dc))
    !write (0,*) "block_vector(4) ", MAXVAL(abs(at_op(nmatrixes2)%block_vector(4)%block3dc))
    DEALLOCATE(temp, temp2, switch_rz)
  END SUBROUTINE init_kspace

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine deallocates all block matrixes already initialized
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE delete_k_space
    USE matrix_coefficients, ONLY: at_op, kspace
    USE matrix_data, ONLY: nmatrixes2
    USE picsar_precision, ONLY: idp
    INTEGER(idp)  :: i
    DO i = 1,10
       DEALLOCATE(kspace(nmatrixes2)%block_vector(i)%block3dc)
    ENDDO
    DO i=1,4
       DEALLOCATE(at_op(nmatrixes2)%block_vector(i)%block3dc)
    ENDDO
    DEALLOCATE(kxc,kxb,kxf,kyc,kyb,kyf,kzc,kzb,kzf,krc)

  END SUBROUTINE delete_k_space

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes K-vectors along X,Y,Z
  !
  !> @author
  !> Haithem Kallala
  !
  !> @params[in] l_stg LOGICAL(lp) - Assumes staggered grid for l_stg==.TRUE.
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE compute_k_vec(l_stg)
    USE fields, ONLY: nordery, norderz, norderx, l_AM_rz
    USE group_parameters, ONLY: p3d_fend, p3d_fsize, p3d_fstart
    USE iso_c_binding
    USE mpi_fftw3, ONLY: local_nz, local_nz_tr, local_z0, local_z0_tr
    USE picsar_precision, ONLY: idp, num, lp, cpx
    USE shared_data, ONLY: dx, dy, dz, fftw_mpi_transpose, fftw_with_mpi, nz,        &
      p3dfft_flag, p3dfft_stride, nmodes
    IMPLICIT NONE
    LOGICAL(lp), INTENT(IN)                     :: l_stg
    COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     ::                                     &
     kxct,kxbt,kxft,kyct,kybt,kyft,kzct,kzbt,kzft, k_temp
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: i, nfftx, nffty, nfftz
    REAL(num)                                     :: sd
    INTEGER(idp)                                  :: temp_order

#if defined(LIBRARY)
    !Due to different staggering in PICSAR and SMILEI
    ii = DCMPLX(0.0_num, -1.0_num)
#else
    ii = DCMPLX(0.0_num, 1.0_num)
#endif
    !> If fftw_mpi_transpose then use FFTW_MPI_TRANSPOSED_OUT/IN plans
    !> fftw_mpi_transpose avoids spurious mpi_alltoall call for each
    !> fftw_mpi_exec call. (initially fftw_mpi_exec call mpi_alltoall two
    !> times to perform global data transposition along y and z axis)
    !> Hence fftw_mpi_exec is faster when using transposed plans
    !> But the user should keep in mind that fourier fields are then transposed
    !> in memory and data splitting along mpi procs is done along y axis instead
    !of z 
    !> axis  with regular fftw_mpi.
    !> block matrixes are also transposed conveniently  during init_gpstd when
    !> using transposed plans
    !> A similar optimization is possible when using p3dfft (p3dfft_stride =
    !.TRUE.)  but z and x axis are then transposed 


    !> If fftw_mpi_transpose, y and z are transposed so nordery, norderz, and  dy 
    !> dz are switched
    IF(fftw_mpi_transpose) THEN
      sd=dz
      dz=dy
      dy=sd
      temp_order = norderz
      norderz=nordery
      nordery=temp_order
    ENDIF
    !> if p3dfft_stride x and z are transposed, so norderz, norderx and dz, dy
    !> are switched
    IF(p3dfft_flag) THEN
      IF(p3dfft_stride) THEN
        sd = dz
        dz = dx
        dx = sd
        temp_order = norderz
        norderz = norderx
        norderx = temp_order
      ENDIF
    ENDIF

    !> computes fourier space size for all the group (if hybrid)  
    !> or only locally (if local psatd) or for the whole domain(if gloal psatd)
    CALL select_case_dims_global(nfftx, nffty, nfftz)

    !> computes wave vector components in each direction
    !> for the case of spectral AM RZ kr is always centered 
    !> nfftx == nfftr and dx= =dr , kr depend on the mode 
    !> nffty == nfftl this one can be staggered ky == kl and dy==dl
    IF (.NOT. l_AM_rz) THEN
      CALL compute_k_1d( nfftx,kxc,kxf,kxb,norderx,dx,l_stg)
      CALL compute_k_1d( nffty,kyc,kyf,kyb,nordery,dy,l_stg)
      CALL compute_k_1d( nfftz,kzc,kzf,kzb,norderz,dz,l_stg)
    ELSE IF (l_AM_rz) THEN
      CALL  compute_kr_1d(nfftx,krc,dx,nmodes) 
      CALL  compute_k_1d( nffty,kyc,kyf,kyb,nordery,dy,l_stg)
      !write (0,*) "dy = ", dy, "nordery= ", nordery
    END IF
    
    !write (0,*), "===========================kr=============================="
     
    ! Selects only haf of  kx because r2c and c2r ffts
    ! the case not l_AM_rz is added to this one because we perform a complex to
    ! complex transform 
    IF((.NOT. p3dfft_flag) .AND. (.NOT. l_AM_rz) ) THEN
      ALLOCATE(k_temp(nfftx));
      k_temp = kxc;
      DEALLOCATE(kxc); ALLOCATE(kxc(nfftx/2+1)) ; kxc = k_temp(1:nfftx/2+1)
      k_temp = kxb;
      DEALLOCATE(kxb); ALLOCATE(kxb(nfftx/2+1)) ; kxb = k_temp(1:nfftx/2+1)
      k_temp = kxf;
      DEALLOCATE(kxf); ALLOCATE(kxf(nfftx/2+1)) ; kxf = k_temp(1:nfftx/2+1)
      DEALLOCATE(k_temp)
    ENDIF

    !> Selects only relevent wave vector components for each processor    
    IF(fftw_with_mpi) THEN
      IF( .NOT. p3dfft_flag) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          ALLOCATE(k_temp(nfftz))
          k_temp = kzc
          DEALLOCATE(kzc);ALLOCATE(kzc(local_nz))
          kzc = k_temp(local_z0+1:local_z0+local_nz)
          k_temp = kzf
          DEALLOCATE(kzf);ALLOCATE(kzf(local_nz))
          kzf = k_temp(local_z0+1:local_z0+local_nz)
          k_temp = kzb
          DEALLOCATE(kzb);ALLOCATE(kzb(local_nz))
          kzb = k_temp(local_z0+1:local_z0+local_nz)
        ELSE IF(fftw_mpi_transpose) THEN
          ALLOCATE(k_temp(nfftz))
          k_temp = kzc
          DEALLOCATE(kzc);ALLOCATE(kzc(local_nz_tr))
          kzc = k_temp(local_z0_tr+1:local_z0_tr+local_nz_tr)
          k_temp = kzf
          DEALLOCATE(kzf);ALLOCATE(kzf(local_nz_tr))
          kzf = k_temp(local_z0_tr+1:local_z0_tr+local_nz_tr)
          k_temp = kzb
          DEALLOCATE(kzb);ALLOCATE(kzb(local_nz_tr))
          kzb = k_temp(local_z0_tr+1:local_z0_tr+local_nz_tr)
        ENDIF
        DEALLOCATE(k_temp)
      ELSE IF(p3dfft_flag) THEN
          ALLOCATE(kxct(nfftx),kxbt(nfftx),kxft(nfftx),kyct(nffty),kybt(nffty),        &
          kyft(nffty),kzct(nfftz),kzbt(nfftz),kzft(nfftz))
          kxct = kxc; kxbt = kxb ; kxft = kxf ;
          kyct = kyc; kybt = kyb ; kyft = kyf ;
          kzct = kzc; kzbt = kzb ; kzft = kzf ;
          DEALLOCATE(kxc,kxf,kxb,kyc,kyf,kyb,kzc,kzf,kzb)

          ALLOCATE(kxc(p3d_fsize(1)),kxf(p3d_fsize(1)),kxb(p3d_fsize(1)))
          ALLOCATE(kyc(p3d_fsize(2)),kyf(p3d_fsize(2)),kyb(p3d_fsize(2)))
          ALLOCATE(kzc(p3d_fsize(3)),kzf(p3d_fsize(3)),kzb(p3d_fsize(3)))
          kxc = kxct(p3d_fstart(1):p3d_fend(1))
          kxb = kxbt(p3d_fstart(1):p3d_fend(1))
          kxf = kxft(p3d_fstart(1):p3d_fend(1))
          kyc = kyct(p3d_fstart(2):p3d_fend(2))
          kyb = kybt(p3d_fstart(2):p3d_fend(2))
          kyf = kyft(p3d_fstart(2):p3d_fend(2))
          kzc = kzct(p3d_fstart(3):p3d_fend(3))
          kzb = kzbt(p3d_fstart(3):p3d_fend(3))
          kzf = kzft(p3d_fstart(3):p3d_fend(3))
          DEALLOCATE(kxct,kxbt,kxft,kyct,kybt,kyft,kzct,kzbt,kzft)
       ENDIF
    ENDIF

    !> If fftw_mpi_transpose , p3dfft_stride reswitch parameters
    IF(fftw_mpi_transpose) THEN
      sd=dz
      dz=dy
      dy=sd
      temp_order = norderz
      norderz=nordery
      nordery=temp_order
    ENDIF
    IF(p3dfft_flag) THEN
      IF(p3dfft_stride) THEN
        sd = dz
        dz = dx
        dx = sd
        temp_order = norderz
        norderz = norderx
        norderx = temp_order
      ENDIF
    ENDIF
  END SUBROUTINE compute_k_vec

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes a 1D k-vector along a given direction
  !
  !> @author
  !> Haithem Kallala
  !
  !> @params[in] nfft INTEGER(idp) - number of points on which the FFT is performed on
  !> current axis
  !> @params[in] norder - INTEGER(idp) - stencil spatial order
  !> @params[in] d  - REAL(num) - sampling period in real space
  !> @params[in] l_stg - LOGICAL(lp) - Assumes staggered grid for l_stg==.TRUE.
  !> @params[in,out] kvec - array of REAL(num) - kvector
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE compute_k_1d(nfft,kvec,kvecf,kvecb,norder,d,l_stg)
     USE constants, ONLY: pi
     USE picsar_precision, ONLY: cpx, idp, lp, num
     REAL(num) , INTENT(IN)  :: d
     INTEGER(idp) , INTENT(IN) :: norder,nfft
     COMPLEX(cpx) , DIMENSION(:) , ALLOCATABLE , INTENT(INOUT) :: kvec,kvecf,kvecb
     COMPLEX(cpx) , DIMENSION(:) , ALLOCATABLE :: k_true
     LOGICAL(lp)  , INTENT(IN)             :: l_stg
     COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     ::  ones, onesp
     REAL(num), ALLOCATABLE, DIMENSION(:)        :: FD
     INTEGER(idp)                                ::j,i
     COMPLEX(cpx)                                ::  ii

     ii = (0.0_num,1.0_num)

     ALLOCATE(ones(nfft), onesp(nfft))
     ALLOCATE(kvec(nfft),kvecf(nfft),kvecb(nfft), k_true(nfft))
     kvec=(0._num, 0._num)
     kvecb=(0._num, 0._num)
     kvecf=(0._num, 0._num)
     k_true=(0._num, 0._num)
     DO j=1_idp, nfft
       ones(j)  = DCMPLX(j-1.0_num, 0.0_num)
       onesp(j) = DCMPLX(j-1.0_num, 0.0_num)
       IF(j .GT. nfft/2_idp +1) THEN
         ones(j)  =DCMPLX(-ones(j))
         onesp(j) =DCMPLX( nfft + onesp(j))
       ENDIF
     ENDDO
     ! > By convention, if norder == 0 then computes derivative operatior with
     ! > infinite order stencil

     IF (norder .ne. 0_idp) THEN
       ALLOCATE(FD(norder/2))
       !> Computes finite difference coefficients for a staggered or an
       !> unstaggered grid
       CALL FD_weights_hvincenti(norder, FD, l_stg)
       DO i=1_idp, norder/2
         IF(l_stg) THEN
            kvec=kvec+2.0_num/d*FD(i)*SIN((i*2.0_num-1.0_num)*PI*ones/nfft)
         ELSE
            CALL fftfreq(nfft, k_true,  d)
            kvec=kvec+2.0_num/d*FD(i)*SIN(i*k_true*d)
            !kvec=kvec+2.0_num/d*FD(i)*SIN(i*2.0_num*PI*onesp/nfft)
         ENDIF
       ENDDO

       DEALLOCATE(FD)
     ELSE
       !> If norder == 0 then computes the exact wave vector
       !> with an infinite stencil
       CALL fftfreq(nfft, kvec,  d)
     ENDIF
     !> If staggered grid then computes staggered derivative operator  
     !> (dual to primal and primal to dual)  

     IF(l_stg) THEN
       kvecf=kvec*EXP(-ii*PI*onesp/nfft)
       kvecb=kvec*EXP(ii*PI*onesp/nfft)
     ELSE
     !> If unstaggered grid 
       kvecb=kvec
       kvecf=kvec
     ENDIF
     !write (0,*) "nfft", nfft, "size", SIZE(kvec), "delta = ", d,  "k_order=", kvec
     DEALLOCATE(onesp,ones)
     DEALLOCATE(k_true)
  END SUBROUTINE compute_k_1d

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes a 1D k-vector along r  direction using the zeros
  !> of Bessel  functions
  !> PS : kr is always centred for the moment and we calculate the exact value !  
  !> for the moment l_stg is always false and norder should be only 0 
  !> @author
  !> Imen Zemzemi
  !
  !> @params[in] nfft INTEGER(idp) - number of points on which the FFT is performed on
  !> current axis
  !> @params[in] norder - INTEGER(idp) - stencil spatial order
  !> @params[in] d  - REAL(num) - sampling period in real space
  !> @params[in] l_stg - LOGICAL(lp) - Assumes staggered grid for l_stg==.TRUE.
  !> @params[in,out] kvec - array of REAL(num) - kvector
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE compute_kr_1d(nfft,kvec,d,nmodes)
     USE shared_data, ONLY: nx
     USE picsar_precision, ONLY: idp, cpx
     REAL(num) , INTENT(IN)  :: d
     INTEGER(idp) , INTENT(IN) :: nfft,nmodes
     COMPLEX(cpx) , DIMENSION(:,:) , ALLOCATABLE , INTENT(INOUT) :: kvec
     COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     ::  ones
     REAL ( idp ), ALLOCATABLE, DIMENSION(:) :: nu, RJ1,RY0,RY1, nu1
     INTEGER(idp) ::  i,k
     ALLOCATE (nu(nfft))
     ALLOCATE (nu1(nfft-1))
     ALLOCATE (RJ1(nfft))
     ALLOCATE (RY0(nfft))
     ALLOCATE (RY1(nfft))
     ALLOCATE(ones(nfft))
     ALLOCATE(kvec(nfft,nmodes))
     kvec=(0._num, 0._num)
     !kr = 2*np.pi * self.trans[m].dht0.get_nu()
     DO k=0,nmodes-1
       IF (k .eq. 0) THEN
          CALL JYZO(k,nfft,nu,RJ1,RY0,RY1)
         DO i=1,nfft
           kvec(i,k+1) = nu(i)/(nfft*d) 
         END DO
       ELSE
         CALL JYZO(k,nfft-1,nu1,RJ1,RY0,RY1)
         kvec(1,k+1)=0._num 
         DO i=2,nfft
           kvec(i,k+1)=nu1(i-1)/(nfft*d)
         END DO
       END IF
     END DO

  DEALLOCATE (nu, nu1, RJ1, RY0, RY1, ones)

  END SUBROUTINE compute_kr_1d

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes a 1D k-vector along a given direction
  !
  !> @author
  !> H. Vincenti
  !
  !> @params[in] nxx - INTEGER(idp) - number of points on which the FFT is performed on
  !> current axis
  !> @params[in] dxx - REAL(num) - sampling period along current axis
  !> @params[in,out] kxx - array of REAL(num) - kvector along current durection
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________

  SUBROUTINE fftfreq(nxx, kxx, dxx)
    USE constants, ONLY: pi
    USE picsar_precision, ONLY: cpx, idp, num
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)                    :: nxx
    REAL(num), INTENT(IN)                    :: dxx
    COMPLEX(cpx), INTENT(OUT), DIMENSION(nxx)  :: kxx
    INTEGER(idp) :: i, n
    REAL(num) :: fe

    fe=1_num/dxx
    n=nxx
    kxx(1)=DCMPLX(0.0_num, 0.0_num)
    IF (MOD(n, 2) .EQ. 0)THEN
      ! First part of k [0, ..., n/2-1]
      DO i=1, n/2_idp-1_idp
        kxx(i+1)=kxx(i)+(1., 0.)
      END DO
      ! Second part of k [-n/2, -1]
      kxx(n/2_idp+1)=-n/2_idp
      DO i=n/2_idp+1, n-1
        kxx(i+1)=kxx(i)+DCMPLX(1.0_num, 0.0_num)
      END DO
    ELSE
      ! First part of k [0, ..., (n-1)/2]
      DO i=1, (n-1_idp)/2_idp
        kxx(i+1)=kxx(i)+(1., 0.)
      END DO
      ! Second part of k [-(n-1)/2, -1]
      kxx((n-1_idp)/2_idp+2_idp)=-DCMPLX((n-1_idp)/2_idp, 0.0_num)
      DO i=(n-1_idp)/2_idp+2_idp, n-1
        kxx(i+1)=kxx(i)+(1.0_num, 0.0_num)
      END DO
    ENDIF
    kxx=kxx/(dxx*nxx)*2.0_num*PI
  END SUBROUTINE fftfreq

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine allocated block matrixes  
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE init_gpstd()
    USE constants, ONLY: clight, eps0, mu0
    USE fields, ONLY: bxf, byf, bzf, exf, eyf, ezf, g_spectral, jxf, jyf, jzf, rhof, &
      rhooldf, l_AM_rz,  invM, invM_1, invM1, Ma, Ma_1, Ma1
    USE iso_c_binding
    USE matrix_coefficients, ONLY: at_op, cc_mat, kspace, vnew, vold, hankel_mat
    USE matrix_data, ONLY: nmatrixes, nmatrixes2, nmatrixes_h
    USE mpi_fftw3, ONLY: alloc_local, fftw_alloc_complex
    USE omp_lib
    USE params, ONLY: dt
    USE picsar_precision, ONLY: idp, num, lp, cpx
    USE shared_data, ONLY: absorbing_bcs, c_dim, fftw_with_mpi, nkx, nky, nkz, nx,   &
      ny, nz, p3dfft_flag, nmodes

    INTEGER(idp)           :: i, j, imode
    COMPLEX(cpx)           :: ii
    INTEGER(idp)           :: nfftx, nffty, nfftz,nfftxr, nbloc_ccmat, nbloc_vnew, nblock_hankel
    LOGICAL(lp)            :: switch
    REAL(num)              :: coeff_norm
    TYPE(C_PTR)            :: cdata
 
    IF(absorbing_bcs) THEN
      !> When using pmls, cc_mat is a 12x17 matrix
      nbloc_ccmat = 17_idp
      nbloc_vnew = 12_idp
    ELSE IF(.NOT. absorbing_bcs) THEN
      !> When using peridic bcs, cc_mat is a 6x11 matrix
      nbloc_ccmat = 11_idp
      nbloc_vnew = 6_idp
    ENDIF
  
    CALL select_case_dims_local(nfftx, nffty, nfftz)
    ii=DCMPLX(0.0_num, 1.0_num)
    CALL allocate_new_matrix_vector(nbloc_ccmat)
    nblock_hankel  = 6_idp
    CALL allocate_new_hankel_matrix_vector(nfftz, nblock_hankel)
    nfftxr = nfftx/2+1
    IF(p3dfft_flag) nfftxr = nfftx
    IF (l_AM_rz) THEN
      nfftxr = nfftx
    ENDIF  
    CALL init_kspace
    !> Here again in the call of init_space we specify if we are using the spectral cylindrical 
    !> and nfftxr should be nfftr, nffty is nfftl and nfftz is nmodes
    nkx = nfftxr
    nky = nffty
    IF (l_AM_rz) THEN
      nkz=nmodes
    ELSE
      nkz = nfftz
    END IF
    !> Allocates cc_mat  block matrix
    !> cc_mat blocks are initally as an nbloc_ccmat x nbloc_ccmat block matrix 
    !> At the end of the routine, useless blcoks are deleted  
    DO i=1_idp, nbloc_ccmat
      DO j=1_idp, nbloc_ccmat
        ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i, j)%block3dc(nkx, nky,    &
        nkz))
        cc_mat(nmatrixes)%block_matrix2d(i, j)%nx = nkx
        cc_mat(nmatrixes)%block_matrix2d(i, j)%ny = nky
        cc_mat(nmatrixes)%block_matrix2d(i, j)%nz = nkz
      ENDDO
    ENDDO

    DO i=1_idp, nfftz
      DO j=1_idp, nblock_hankel
        ALLOCATE(hankel_mat(nmatrixes_h)%mode_block_matrix2d(i, j)%block2dc(nkx, nkx))    
        hankel_mat(nmatrixes_h)%mode_block_matrix2d(i, j)%nx=nkx
      ENDDO
    ENDDO

    

    !> If g_spectral then psatd uses multiply_mat_vec routine in GPSTD.F90 
    !> to perform the maxwell push in Fourier space
    !> So we need to allocate vold/vnew vector blocks
    !> else if g_spectral == false these arrays are not allocated, and
    !> push_psaotd_ebfields_3d/2d is used to perform the maxwell push in Fourier space

    !> When using absorbing_bcs, g_spectral = .TRUE. is needed
    IF(g_spectral) THEN
      IF(p3dfft_flag) THEN  ! hybrid with p3dfft
        DO i = 1,nbloc_ccmat
          ALLOCATE(vold(nmatrixes)%block_vector(i)%block3dc(nkx,nky,nkz))
        ENDDO
        DO i = 1,nbloc_vnew
          ALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc(nkx,nky,nkz))
        ENDDO
      ELSE IF(fftw_with_mpi) THEN ! hybrid or global with fftw
        DO i =1,nbloc_ccmat
          cdata = fftw_alloc_complex(alloc_local)
          CALL c_f_pointer(cdata, vold(nmatrixes)%block_vector(i)%block3dc, [nkx, nky, nkz])
        ENDDO
        DO i=1,nbloc_vnew
          cdata = fftw_alloc_complex(alloc_local)
          CALL c_f_pointer(cdata, vnew(nmatrixes)%block_vector(i)%block3dc,[nkx, nky, nkz])
        ENDDO
      ELSE IF(.NOT. fftw_with_mpi) THEN ! local psatd
        DO i = 1,nbloc_ccmat
          ALLOCATE(vold(nmatrixes)%block_vector(i)%block3dc(nkx,nky,nkz))
        ENDDO
        DO i = 1,nbloc_vnew
          ALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc(nkx,nky,nkz))
        ENDDO
      ENDIF
      DO i = 1,nbloc_ccmat
        vold(nmatrixes)%block_vector(i)%nx = nfftxr
        vold(nmatrixes)%block_vector(i)%ny = nffty
        vold(nmatrixes)%block_vector(i)%nz = nfftz
      ENDDO
      DO i=1,nbloc_vnew
        vnew(nmatrixes)%block_vector(i)%nx = nfftxr
        vnew(nmatrixes)%block_vector(i)%ny = nffty
        vnew(nmatrixes)%block_vector(i)%nz = nfftz
      ENDDO
      DO i=nbloc_vnew + 1_idp,nbloc_ccmat
        ALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc(1,1,1))
        vnew(nmatrixes)%block_vector(i)%nx = 1
        vnew(nmatrixes)%block_vector(i)%ny = 1
        vnew(nmatrixes)%block_vector(i)%nz = 1
      ENDDO
    ENDIF

    !> Init all blocks to 0.0
    DO i = 1,nbloc_ccmat
      DO j=1,nbloc_ccmat
        cc_mat(nmatrixes)%block_matrix2d(i, j)%block3dc = CMPLX(0.0_num,0.0_num)
      ENDDO
    ENDDO

    DO i = 1,nfftz
      DO j=1,nblock_hankel
        hankel_mat(nmatrixes_h)%mode_block_matrix2d(i, j)%block2dc = CMPLX(0.0_num,0.0_num)
      ENDDO
    ENDDO
 

    IF (absorbing_bcs) THEN
      !> When using pmls, splitted fields EM equations are solved
      !> The following routine solved these pushes fourier splitted fields in
      !> fourier space
      !> ex = exy + exz
      !> ey = eyx + eyz  ...

      !> Current contribution to EM acts only on the first half component of the
      !> field component
      !> Thus current and density contributions only acts influences
      !> exy , eyx , ezx, bxy, byx, bzx
      !> Hence, contribution of J to exz is null
      CALL compute_cc_mat_splitted_fields()
    ELSE IF(.NOT. absorbing_bcs) THEN  
      !> When not using absorbing bcs, standard EM equations are solved in
      ! fourier space
      IF (l_AM_rz) THEN
        DO imode =1,nfftz
          Call Hankel_M_and_invM(imode-1)  
          hankel_mat(nmatrixes_h)%mode_block_matrix2d(imode, 1)%block2dc = Ma
          hankel_mat(nmatrixes_h)%mode_block_matrix2d(imode, 2)%block2dc = Ma_1
          hankel_mat(nmatrixes_h)%mode_block_matrix2d(imode, 3)%block2dc = Ma1
          hankel_mat(nmatrixes_h)%mode_block_matrix2d(imode, 4)%block2dc = invM
          hankel_mat(nmatrixes_h)%mode_block_matrix2d(imode, 5)%block2dc = invM_1
          hankel_mat(nmatrixes_h)%mode_block_matrix2d(imode, 6)%block2dc = invM1

        ENDDO
        DEALLOCATE(Ma,Ma_1,Ma1,invM,invM_1,invM1)
        !write (0,*) " computing merged fields in RZ"
        CALL compute_cc_mat_merged_fields_AM_rz()
      ELSE
        CALL compute_cc_mat_merged_fields()
      ENDIF
    ENDIF
  
 
    !> Renormalize cc_mat blocks
    !> Because fftw_r2c and followed by fftw_c2r multiplies fields by 
    !> nfftx*nffty*nfftz 
    !> This way, no need to normalize fields in a separate step


    ! Introduce fft normalisation factor in mat bloc mult
    CALL select_case_dims_global(nfftx,nffty,nfftz)
    !> Here i'm not sure if it shoud be divided by nfftx* nffty or just nffty ....
    IF (l_AM_rz) THEN
      coeff_norm = 1.0_num
    ELSE 
      coeff_norm = 1.0_num/(nfftx*nffty*nfftz)
    ENDIF
    DO i=1,nbloc_ccmat
      DO j=1,nbloc_ccmat
          cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc =                            &
          coeff_norm*cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc
      ENDDO
    ENDDO

    !> Delete uninitialized blocks
    DO i=1,nbloc_ccmat
      DO j=1,nbloc_ccmat
        IF(sum(abs(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc))                   &
        /size(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc)  == 0.0_num) THEN
          DEALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc)
          ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc(1,1,1))
          cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc(1,1,1) = (0._num,0._num)
          cc_mat(nmatrixes)%block_matrix2d(i,j)%nx = 1
          cc_mat(nmatrixes)%block_matrix2d(i,j)%ny = 1
          cc_mat(nmatrixes)%block_matrix2d(i,j)%nz = 1
        ENDIF
      ENDDO
    ENDDO

    !> Delete kspace and at_op blocks
    !> Might not delete these blocks if current filtering or field correction is
    !> needed in Fourier space
    !CALL delete_k_space
  END SUBROUTINE init_gpstd

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine inits block matrixes with splitted fields EM equations when
  !> using PMLs
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________


  SUBROUTINE compute_cc_mat_splitted_fields()
    USE constants, ONLY: clight, eps0, mu0
    USE matrix_coefficients, ONLY: at_op, cc_mat, kspace
    USE matrix_data, ONLY: nmatrixes, nmatrixes2
    USE params, ONLY: dt
    USE picsar_precision, ONLY: cpx, idp, lp, num
    INTEGER(idp) :: i,j,k
    COMPLEX(cpx) ::  ii
    LOGICAL(lp)  :: switch

    !> cc_mat_(nmatrixes)block_matrix2d(i,j) components are sorted using the
    !>following nomenclature 
    !> In this case cc_mat is a 12x17 block matrix
    !> 1-> exyf; 2->exzf; 3->eyxf; 4->eyzf; 5->ezxf; 6->ezyf
    !> 7-> bxyf; 8->bxzf; 9->byxf; 10->byzf; 11->bzxf; 12->bzyf
    !> 13-> jxf; 14->jyf; 15->jzf; 16->rhooldf; 17->rhof
    !> cc_mat_(nmatrixes)block_matrix2d(i,j) is the contribution of the j-th
    !> scalar field to the i-th scalar field 

    ii=DCMPLX(0.0_num, 1.0_num)

    !> Contribution of E field to E field and B field to B field
    DO i=1, 12
     cc_mat(nmatrixes)%block_matrix2d(i, i)%block3dc =&
     AT_OP(nmatrixes2)%block_vector(2)%block3dc
    ENDDO
     ! contribution of current to elec field 
     ! by convention, j will only contribute to exy, eyx, ezx 
    DO i = 1, 3
      j = 2*(i-1) + 1
      k = i+12
      cc_mat(nmatrixes)%block_matrix2d(j,k)%block3dc = &
      (-1._num)*clight*mu0*AT_OP(nmatrixes2)%block_vector(1)%block3dc
    ENDDO

    !contribution rho old by convention only contributes to exy ,eyx ezx
    switch = .FALSE.

    !> Spots mpis that contain null frequency to perform Taylor expansion later
    IF(ABS(Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)) .EQ. 0.0_num)    THEN
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = (1.0_num, 0.0_num)
      switch = .TRUE.
    ENDIF
    ! Contribution of rhooldf to E
    !> rhoold only contributes to bxy, byx, bzx

    DO i = 1, 3
      j = 2*(i-1)+1
      cc_mat(nmatrixes)%block_matrix2d(j, 16_idp)%block3dc = DCMPLX(0.,&
      1.)*(AT_OP(nmatrixes2)%block_vector(2)%block3dc&
      -1./(clight*dt)*AT_OP(nmatrixes2)%block_vector(1)%block3dc)&
      /Kspace(nmatrixes2)%block_vector(10)%block3dc**2
      cc_mat(nmatrixes)%block_matrix2d(j, 16_idp)%block3dc =&
      cc_mat(nmatrixes)%block_matrix2d(j, 16_idp)%block3dc&
      *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc  
     IF(switch) THEN
        cc_mat(nmatrixes)%block_matrix2d(j, 16_idp)%block3dc(1, 1, 1) =&
        -1.0_num/3.0_num*(0.0_num, 1.0_num)*(clight*dt)**2    
     ENDIF

      !> If current mpi task contains null frequency then performs Taylor
      !expansion for cc_mat(nmatrixes)%block_matrix2d(i, 16_idp)%block3dc(1, 1,
      !1)
      cc_mat(nmatrixes)%block_matrix2d(j, 16_idp)%block3dc = 1.0_num/eps0&
      *cc_mat(nmatrixes)%block_matrix2d(j, 16_idp)%block3dc
    ENDDO
    !> End contribution rhooldf to E
    
    !> Begin contribution rhof to E
    !> rho only contributes to bxy, byx, bzx

    DO i = 1, 3
      j = 2*(i-1)+1
      cc_mat(nmatrixes)%block_matrix2d(j, 17_idp)%block3dc = DCMPLX(0.,&
      1.)*(1./(clight*dt)* AT_OP(nmatrixes2)%block_vector(1)%block3dc-DCMPLX(1.,     &
      0.))/Kspace(nmatrixes2)%block_vector(10)%block3dc**2
      cc_mat(nmatrixes)%block_matrix2d(j, 17_idp)%block3dc =&
      cc_mat(nmatrixes)%block_matrix2d(j, 17_idp)%block3dc&
      *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
      IF(switch) THEN
        cc_mat(nmatrixes)%block_matrix2d(j, 17_idp)%block3dc(1, 1, 1) =&
        -1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
      ENDIF

      !> If current mpi task contains null frequency then performs Taylor
      !expansion for cc_mat(nmatrixes)%block_matrix2d(i, 17_idp)%block3dc(1, 1,
      !1)
      cc_mat(nmatrixes)%block_matrix2d(j, 17_idp)%block3dc = 1.0_num/eps0 *&
      cc_mat(nmatrixes)%block_matrix2d(j, 17_idp)%block3dc
    ENDDO
    !> END contribution rhof to E    
    IF(switch) THEN
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)   = DCMPLX(0., 0.)
    ENDIF

    !>  Begin Contribution of j to B
    !> j only contributes to bxy, byx, bzx
    cc_mat(nmatrixes)%block_matrix2d(7, 14)%block3dc = - mu0*&
    ii*Kspace(nmatrixes2)%block_vector(8)%block3dc*&
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(7, 15)%block3dc = -&
    mu0*(-ii)*Kspace(nmatrixes2)%block_vector(5)%block3dc*&
    AT_OP(nmatrixes2)%block_vector(3)%block3dc



    cc_mat(nmatrixes)%block_matrix2d(9, 13)%block3dc = -&
    mu0*(-ii)*Kspace(nmatrixes2)%block_vector(8)%block3dc*&
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(9, 15)%block3dc = - mu0*&
    ii*Kspace(nmatrixes2)%block_vector(2)%block3dc*&
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(11, 13)%block3dc = - mu0*&
    ii*Kspace(nmatrixes2)%block_vector(5)%block3dc*&
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(11, 14)%block3dc =&
    -mu0*(-ii)*Kspace(nmatrixes2)%block_vector(2)%block3dc*&
    AT_OP(nmatrixes2)%block_vector(3)%block3dc
    !> End contribution J to B

    !> Begin contribution of B to E

    cc_mat(nmatrixes)%block_matrix2d(2,9)%block3dc = -&
    ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(2,10)%block3dc = -&
    ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(1, 11)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(1, 12)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(4, 7)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(3, 11)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 12)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc



    cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc = -&
    ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc = -&
    ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 10)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    !> End contribution of B to E
 
    !> Begin contribution E to B
 
    cc_mat(nmatrixes)%block_matrix2d(8, 4)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(8, 3)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(7, 5)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(7, 6)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(10, 1)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
   
    cc_mat(nmatrixes)%block_matrix2d(10, 2)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
    cc_mat(nmatrixes)%block_matrix2d(9, 5)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
    cc_mat(nmatrixes)%block_matrix2d(9, 6)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
   
    cc_mat(nmatrixes)%block_matrix2d(12, 1)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
    cc_mat(nmatrixes)%block_matrix2d(12, 2)%block3dc =&
    ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
    cc_mat(nmatrixes)%block_matrix2d(11, 3)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
   
    cc_mat(nmatrixes)%block_matrix2d(11, 4)%block3dc =&
    -ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight&
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
    !> End contribution E to B
    END SUBROUTINE compute_cc_mat_splitted_fields

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine inits block matrixes with standard fields EM equations when
  !> NOT using PMLs
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________


  SUBROUTINE compute_cc_mat_merged_fields()
    USE constants, ONLY: clight, eps0, mu0
    USE matrix_coefficients, ONLY: at_op, cc_mat, kspace
    USE matrix_data, ONLY: nmatrixes, nmatrixes2
    USE params, ONLY: dt
    USE picsar_precision, ONLY: cpx, idp, lp, num
    INTEGER(idp)  :: i,j   
    COMPLEX(cpx) ::  ii
    LOGICAL(lp)            :: switch
    ii=DCMPLX(0.0_num, 1.0_num)

    !> cc_mat_(nmatrixes)block_matrix2d(i,j) components are sorted using the
    !>following nomenclature 
    !> In this case cc_mat is a 6x11 block matrix
    !> 1-> exf; 2->eyf; 3->ezf; 4->bxf; 5->byf; 6->bzf
    !> 7->jxf; 8->jyf; 9->jzf; 10-> rhooldf; 11->rhof
    !> cc_mat_(nmatrixes)block_matrix2d(i,j) is the contribution of the j-th
    !> scalar field to the i-th scalar field 
    

  
    !> Contribution of B field to E field update
    cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc = -                               &
    ii*kspace(nmatrixes2)%block_vector(7)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc
    cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(7)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(2, 6)%block3dc =                                 &
    -ii*kspace(nmatrixes2)%block_vector(1)%block3dc*clight                            &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc = -                               &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 5)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(1)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc
    
    !> End contribution B field to E field
    
    !> Contribution of E field to B field
    cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(8)%block3dc/clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc =                                 &
    -ii*kspace(nmatrixes2)%block_vector(5)%block3dc/clight                            &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc =                                 &
    -ii*kspace(nmatrixes2)%block_vector(8)%block3dc/clight                            &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 3)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(2)%block3dc/clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(5)%block3dc/clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6, 2)%block3dc =                                 &
    -ii*kspace(nmatrixes2)%block_vector(2)%block3dc/clight                            &
    *at_op(nmatrixes2)%block_vector(1)%block3dc
   
    !> End contribiton E field to B field
 
    !> Contribution of E field to E field and B field to B field
    DO i=1, 6
      cc_mat(nmatrixes)%block_matrix2d(i, i)%block3dc =                               &
      at_op(nmatrixes2)%block_vector(2)%block3dc
    ENDDO
    !> End contribution of E field To E field and B field to B field    

    !> Contribution of J field to E field
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, i+6)%block3dc =                             &
      (-1._num)*clight*mu0*at_op(nmatrixes2)%block_vector(1)%block3dc
    ENDDO
    ! End contribution of J field to E field

    !> Contribution of J field to B field
    cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc = - mu0*                          &
    ii*kspace(nmatrixes2)%block_vector(8)%block3dc*                                   &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc = -                               &
    mu0*(-ii)*kspace(nmatrixes2)%block_vector(5)%block3dc*                            &
    at_op(nmatrixes2)%block_vector(3)%block3dc



    cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc = -                               &
    mu0*(-ii)*kspace(nmatrixes2)%block_vector(8)%block3dc*                            &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc = - mu0*                          &
    ii*kspace(nmatrixes2)%block_vector(2)%block3dc*                                   &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc = - mu0*                          &
    ii*kspace(nmatrixes2)%block_vector(5)%block3dc*                                   &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc =                                 &
    -mu0*(-ii)*kspace(nmatrixes2)%block_vector(2)%block3dc*                           &
    at_op(nmatrixes2)%block_vector(3)%block3dc
    
    !> End contribution of J field to B field

    !> Contribution of rhoold field to E field

    !> if current mpi task contains the null frequency then this processor it
    !> tagged by switch = .TRUE. in order perform Taylor expansion
    !> for certain blocks

    switch = .FALSE.
    !> Spots mpis that contain null frequency to perform Taylor expansion later
    IF(ABS(kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)) .EQ. 0.0_num) THEN
      kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = (1.0_num, 0.0_num)
      switch = .TRUE.
    ENDIF
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc = DCMPLX(0.,               &
      1.)*(at_op(nmatrixes2)%block_vector(2)%block3dc                                 &
      -1./(clight*dt)*at_op(nmatrixes2)%block_vector(1)%block3dc)                     &
      /kspace(nmatrixes2)%block_vector(10)%block3dc**2

      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc =                          &
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc                            &
      *kspace(nmatrixes2)%block_vector(3*i-1)%block3dc

      !> If current mpi task contains null frequency then performs Taylor
      !> expansion for cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc(1, 1,
      !1)
      IF(switch) THEN
        cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc(1, 1, 1) =               &
        -1.0_num/3.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
      ENDIF
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc = 1.0_num/eps0             &
      *cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc
    ENDDO
    !> End contribution of rhoold field to E field
  
    !> Contribution of rho field to E field
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc = DCMPLX(0.,               &
      1.)*(1./(clight*dt)* at_op(nmatrixes2)%block_vector(1)%block3dc -DCMPLX(1.,     &
      0.))/kspace(nmatrixes2)%block_vector(10)%block3dc**2

      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc =                          &
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc                            &
      *kspace(nmatrixes2)%block_vector(3*i-1)%block3dc

      !> If current mpi task contains null frequency then performs Taylor
      !expansion for cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc(1, 1,
      !1)
      IF(switch) THEN
        cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc(1, 1, 1) =               &
        -1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
      ENDIF
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc = 1.0_num/eps0 *           &
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc
    ENDDO
    IF(switch) THEN
      kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)   = DCMPLX(0., 0.)
    ENDIF
    !> End contribution of rho field to E field   
  END SUBROUTINE compute_cc_mat_merged_fields
  
  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine inits block matrixes with fields EM equations in cylindrical azimuthal modes when
  !> NOT using PMLs
  !
  !> @author
  !> Imene Zemzemi
  !
  !> @date
  !> Creation September 2018
  ! ______________________________________________________________________________________


  SUBROUTINE compute_cc_mat_merged_fields_AM_rz()
    USE shared_data
    USE matrix_coefficients
    USE constants
    USE params, ONLY : dt
    IMPLICIT NONE
    INTEGER(idp)  :: i,j, imode   
    COMPLEX(cpx) ::  ii
    LOGICAL(lp), DIMENSION(:), ALLOCATABLE      :: switch
    ii=DCMPLX(0.0_num, 1.0_num)

    !> cc_mat_(nmatrixes)block_matrix2d(i,j) components are sorted using the
    !>following nomenclature 
    !> In this case cc_mat is a 6x11 block matrix
    !> 1-> elf; 2->e+f; 3->e-f; 4->blf; 5->b+f; 6->b-f
    !> 7->jlf; 8->j+f; 9->j-f; 10-> rhooldf; 11->rhof
    !> cc_mat_(nmatrixes)block_matrix2d(i,j) is the contribution of the j-th
    !> scalar field to the i-th scalar field 
    

  
    !> Contribution of B field to E field update
    cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc = -                               &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *0.5*at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(2, 5)%block3dc =                                 &
    kspace(nmatrixes2)%block_vector(1)%block3dc*clight                                &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc = -                               &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *0.5*at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 6)%block3dc = -                               &
    kspace(nmatrixes2)%block_vector(1)%block3dc*clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc
    
    !> End contribution B field to E field
    !> Contribution of E field to B field
    cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc = -                               &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc/clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc = -                               &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc/clight                             &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc/(2.*clight)                        &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 2)%block3dc = -                               &
    kspace(nmatrixes2)%block_vector(1)%block3dc/clight                                &
    *at_op(nmatrixes2)%block_vector(1)%block3dc               

    cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc =                                 &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc/(2.*clight)                        &
    *at_op(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6, 3)%block3dc =                                 &
    kspace(nmatrixes2)%block_vector(1)%block3dc/clight                                &
    *at_op(nmatrixes2)%block_vector(1)%block3dc
   
    !> End contribiton E field to B field
 
    !> Contribution of E field to E field and B field to B field
    DO i=1, 6
      cc_mat(nmatrixes)%block_matrix2d(i, i)%block3dc =                               &
      at_op(nmatrixes2)%block_vector(2)%block3dc
    ENDDO
    !> End contribution of E field To E field and B field to B field    

    !> Contribution of J field to E field
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, i+6)%block3dc =                              &
      (-1._num)*clight*mu0*at_op(nmatrixes2)%block_vector(1)%block3dc
    ENDDO
    ! End contribution of J field to E field

    !> Contribution of J field to B field
    cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc =  mu0*                             &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*                                    &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc =  mu0*                             &
    kspace(nmatrixes2)%block_vector(4)%block3dc*                                       &
    at_op(nmatrixes2)%block_vector(3)%block3dc



    cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc = - mu0*                           &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*0.5*                                 &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(5, 8)%block3dc = mu0*                            &
    kspace(nmatrixes2)%block_vector(1)%block3dc*                                       &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc = - mu0*                           &
    ii*kspace(nmatrixes2)%block_vector(4)%block3dc*0.5*                                 &
    at_op(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(6, 9)%block3dc =  - mu0*                            &
    kspace(nmatrixes2)%block_vector(1)%block3dc*                                        &
    at_op(nmatrixes2)%block_vector(3)%block3dc
    
    !> End contribution of J field to B field

    !> Contribution of rhoold field to E field

    !> if current mpi task contains the null frequency then this processor it
    !> tagged by switch = .TRUE. in order perform Taylor expansion
    !> for certain blocks
    ALLOCATE(switch(nmodes))
    switch = .FALSE.
    DO imode=1, nmodes
      !> Spots mpis that contain null frequency to perform Taylor expansion later
      IF(ABS(kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, imode)) .EQ. 0.0_num) THEN
        kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, imode) = (1.0_num, 0.0_num)
        switch(imode) = .TRUE.
      ENDIF
    END DO

    
    cc_mat(nmatrixes)%block_matrix2d(1, 10_idp)%block3dc = DCMPLX(0.,1.)            &
    *kspace(nmatrixes2)%block_vector(1)%block3dc                                    &
    *(at_op(nmatrixes2)%block_vector(2)%block3dc                                    &
    -1./(clight*dt)*at_op(nmatrixes2)%block_vector(1)%block3dc)                     &
    /kspace(nmatrixes2)%block_vector(10)%block3dc**2


    cc_mat(nmatrixes)%block_matrix2d(2, 10_idp)%block3dc = -                        &
    kspace(nmatrixes2)%block_vector(4)%block3dc                                     &
    *(at_op(nmatrixes2)%block_vector(2)%block3dc                                    &
    -1./(clight*dt)*at_op(nmatrixes2)%block_vector(1)%block3dc)                     &
    /(2.*kspace(nmatrixes2)%block_vector(10)%block3dc**2)

    cc_mat(nmatrixes)%block_matrix2d(3, 10_idp)%block3dc =                          &
    kspace(nmatrixes2)%block_vector(4)%block3dc                                     &
    *(at_op(nmatrixes2)%block_vector(2)%block3dc                                    &
    -1./(clight*dt)*at_op(nmatrixes2)%block_vector(1)%block3dc)                     &
    /(2.*kspace(nmatrixes2)%block_vector(10)%block3dc**2)

     !> IMPORTANT TO REMEMBER TAYLOR EXPANSION SHOULD BE VERIFIED FOR AM_rz
    !> expansion for cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc(1, 1,
    !1)

    DO imode=1, nmodes
      IF (switch(imode)) THEN
        cc_mat(nmatrixes)%block_matrix2d(1, 10_idp)%block3dc(1, 1, imode) =             &
           -1.0_num/3.0_num*DCMPLX(0.0_num, 1.0_num)*(clight*dt)**2

       cc_mat(nmatrixes)%block_matrix2d(2, 10_idp)%block3dc(1, 1, imode) =             &
           1.0_num/3.0_num*(clight*dt)**2

        cc_mat(nmatrixes)%block_matrix2d(3, 10_idp)%block3dc(1, 1, imode) =             &
           -1.0_num/3.0_num*(clight*dt)**2  
      END IF
    END DO
    
    DO i = 1, 3
      DO imode=1,nmodes
        IF(switch(imode)) THEN
          cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc(1, 1, 1) =             &
           -1.0_num/3.0_num*DCMPLX(0.0_num, 1.0_num)*(clight*dt)**2
        ENDIF
      END DO
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc = 1.0_num/eps0           &
      *cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc
    ENDDO
    !> End contribution of rhoold field to E field
 
 
    !> Contribution of rho field to E field
    
    cc_mat(nmatrixes)%block_matrix2d(1, 11_idp)%block3dc = - DCMPLX(0.,1.)          &
    *kspace(nmatrixes2)%block_vector(1)%block3dc                                    & 
    *(DCMPLX(1.,0.)- 1./(clight*dt)* at_op(nmatrixes2)%block_vector(1)%block3dc)    &
    /kspace(nmatrixes2)%block_vector(10)%block3dc**2

    cc_mat(nmatrixes)%block_matrix2d(2, 11_idp)%block3dc =                           &
    kspace(nmatrixes2)%block_vector(4)%block3dc                                     & 
    *(DCMPLX(1.,0.)- 1./(clight*dt)* at_op(nmatrixes2)%block_vector(1)%block3dc)    &
    /(2.*kspace(nmatrixes2)%block_vector(10)%block3dc**2)
    
    cc_mat(nmatrixes)%block_matrix2d(3, 11_idp)%block3dc =   -                        &
    kspace(nmatrixes2)%block_vector(4)%block3dc                                     & 
    *(DCMPLX(1.,0.)- 1./(clight*dt)* at_op(nmatrixes2)%block_vector(1)%block3dc)    &
    /(2.*kspace(nmatrixes2)%block_vector(10)%block3dc**2)


    !> If current mpi task contains null frequency then performs Taylor
    !expansion for cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc(1, 1,
    !1)
    
    DO imode=1, nmodes
      IF (switch(imode)) THEN
        cc_mat(nmatrixes)%block_matrix2d(1, 11_idp)%block3dc(1, 1, imode) =               &
           1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2

        cc_mat(nmatrixes)%block_matrix2d(2, 11_idp)%block3dc(1, 1, imode) =               &
          -DCMPLX(0.,1.)* 1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2

        cc_mat(nmatrixes)%block_matrix2d(3, 11_idp)%block3dc(1, 1, imode) =               &
           DCMPLX(0.,1.)*1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
      END IF
    END DO
  
    !> IMPORTANT TO REMEMBER TAYLOR EXPANSION SHOULD BE VERIFIED FOR AM_rz
    Do i = 1, 3
      DO imode=1,nmodes
        IF(switch(imode)) THEN
          cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc(1, 1, 1) =               &
          -1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
        ENDIF
      END DO
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc = 1.0_num/eps0 *           &
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc
    ENDDO

    DO imode=1, nmodes
      IF(switch(imode)) THEN
        kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, imode)   = DCMPLX(0., 0.)
      ENDIF
    END DO
    !> End contribution of rho field to E field   
  END SUBROUTINE compute_cc_mat_merged_fields_AM_rz

  ! ______________________________________________________________________________________
  !> @brief
  !> This function computes coefficients of order p stencil for centered/staggered
  !> scheme - Taken from H. Vincenti and J-L Vay, CPC, 200, 147 (2016).
  !
  !> @author
  !> H. Vincenti
  !> H. Kallala
  !
  !> @params[in] is_staggered - LOGICAL(lp) - assumes staggered grid if
  !> is_staggered==.TRUE.
  !> @params[in] p - INTEGER(idp) - spatial order p of the stencil
  !> @params[out] w - array of REAL(num) of size p/2 - array containing
  !> stencil coefficients
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE FD_weights_hvincenti(p, w, is_staggered)
    USE picsar_precision, ONLY: idp, lp, num
    IMPLICIT NONE
    LOGICAL(lp) , INTENT(IN) :: is_staggered
    INTEGER(idp), INTENT(IN) :: p
    REAL(num), DIMENSION(p/2), INTENT(OUT) :: w
    INTEGER(idp) :: i, l
    REAL(num) :: lognumer, logdenom

    DO i=1, p/2
      l=i
      IF (is_staggered) THEN
        lognumer =LOG(16.0_num)*(1.0_num-p/2.0_num)+logfactorial(p-1_idp)*2.0_num
        logdenom = LOG(2.0_num*l-1.0_num)*2.0_num+                                    &
        logfactorial(p/2_idp+l-1_idp)+logfactorial(p/2_idp-l)+                        &
        2.0_num*logfactorial(p/2_idp-1_idp)
      ELSE
        lognumer = logfactorial(p/2_idp)*2.0_num
        logdenom = logfactorial(p/2_idp+l)+ logfactorial(p/2_idp-l)+LOG(1.0_num*l)
      ENDIF
      w(i) = (-1.0_num)**(l+1)*EXP(lognumer-logdenom)
    END DO
  END SUBROUTINE FD_weights_hvincenti

  SUBROUTINE copy_field(ex_out, n1, n2, n3, ex_in, nxx, nyy, nzz)
    USE omp_lib
    USE picsar_precision, ONLY: idp, num
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, n1, n2, n3
    REAL(num), DIMENSION(nxx, nyy, nzz), INTENT(IN OUT) :: ex_in
    REAL(num), DIMENSION(n1, n2, n3), INTENT(IN OUT) :: ex_out
    INTEGER(idp) :: ix, iy, iz

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=1, MIN(nzz, n3)
      DO iy=1, MIN(nyy, n2)
        DO ix=1, MIN(nxx, n1)
          ex_out(ix, iy, iz)=ex_in(ix, iy, iz)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE copy_field


! TO finish 
  SUBROUTINE copy_field_forward_AM_rz()
    USE fields, ONLY : el, er, et, bl, br, bt, jl, jr, jt
    USE fields, ONLY : el_c, er_c, et_c, bl_c, br_c, bt_c, jl_c, jr_c, jt_c,rho_c,   &
                       rhoold_c
    USE omp_lib
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: rho, rhoold
    IMPLICIT NONE
    INTEGER(idp) :: il, ir, imode, ill, irr, irrr, illl
    INTEGER(idp) , dimension(3) :: lbound_r, ubound_r, lbound_p,ubound_p,lbound_s, ubound_s

    lbound_r = LBOUND(el_c)
    lbound_p = LBOUND(el)
    ubound_r = UBOUND(el_c)
    ubound_p = UBOUND(el)
    lbound_s = LBOUND(jl)
    ubound_s = UBOUND(jl)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ir, il, imode, irr, ill,irrr,illl) COLLAPSE(3)
    DO imode=lbound_r(3),ubound_r(3)
      DO il=lbound_r(2),ubound_r(2)
        DO ir=lbound_r(1),ubound_r(1)
          irr = ir - lbound_r(1) +lbound_p(1)
          ill = il - lbound_r(2) +lbound_p(2)
          irrr = ir - lbound_r(1) +lbound_s(1)
          illl = il - lbound_r(2) +lbound_s(2)

          el_c(ir, il, imode)=el(irr, ill, imode)
          er_c(ir, il, imode)=er(irr, ill, imode)
          et_c(ir, il, imode)=et(irr, ill, imode)
          bl_c(ir, il, imode)=bl(irr, ill, imode)
          br_c(ir, il, imode)=br(irr, ill, imode)
          bt_c(ir, il, imode)=bt(irr, ill, imode)
          jl_c(ir, il, imode)=jl(irrr, illl, imode)
          jr_c(ir, il, imode)=jr(irrr, illl, imode)
          jt_c(ir, il, imode)=jt(irrr, illl, imode)
          rho_c(ir, il, imode)=rho(irrr, illl, imode)
          rhoold_c(ir, il, imode)=rhoold(irrr, illl, imode)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE copy_field_forward_AM_rz

  SUBROUTINE copy_field_forward()
    USE fields, ONLY: bx, bx_r, bxy, bxy_r, bxz, bxz_r, by, by_r, byx, byx_r, byz,   &
      byz_r, bz, bz_r, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exy, exy_r, exz, exz_r, ey, &
      ey_r, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezx, ezx_r, ezy, ezy_r, jx, jx_r, jy,  &
      jy_r, jz, jz_r, rho_r, rhoold_r
    USE omp_lib
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: absorbing_bcs, nx, ny, nz, rho, rhoold
    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz, ixx, iyy, izz, ixxx, iyyy, izzz
    INTEGER(idp) , dimension(3) :: lbound_r, ubound_r, lbound_p ,ubound_p,lbound_s, ubound_s

    IF(absorbing_bcs) THEN 
       lbound_r = LBOUND(exy_r)
       lbound_p = LBOUND(exy)
       ubound_r = UBOUND(exy_r)
       ubound_p = UBOUND(exy)
       lbound_s = LBOUND(jx)
       ubound_s = UBOUND(jx)
    ELSE
       lbound_r = LBOUND(ex_r)
       lbound_p = LBOUND(ex)
       ubound_r = UBOUND(ex_r)
       ubound_p = UBOUND(ex)
       lbound_s = LBOUND(jx)
       ubound_s = UBOUND(jx)
    ENDIF
    ! When using periodic bcs, standard EM fields are communicated 
    ! Else, when using absorbing bcs, splitted EM fields are communicated 
    IF(.NOT. absorbing_bcs) THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, ixx, iyy ,izz,ixxx,iyyy,izzz) COLLAPSE(3)
      DO iz=lbound_r(3),ubound_r(3)
        DO iy=lbound_r(2),ubound_r(2)
          DO ix=lbound_r(1),ubound_r(1)
            ixx = ix - lbound_r(1) +lbound_p(1)
            iyy = iy - lbound_r(2) +lbound_p(2)
            izz = iz - lbound_r(3) +lbound_p(3)
            ixxx = ix - lbound_r(1) +lbound_s(1)
            iyyy = iy - lbound_r(2) +lbound_s(2)
            izzz = iz - lbound_r(3) +lbound_s(3)

            ex_r(ix, iy, iz)=ex(ixx, iyy, izz)
            ey_r(ix, iy, iz)=ey(ixx, iyy, izz)
            ez_r(ix, iy, iz)=ez(ixx, iyy, izz)
            bx_r(ix, iy, iz)=bx(ixx, iyy, izz)
            by_r(ix, iy, iz)=by(ixx, iyy, izz)
            bz_r(ix, iy, iz)=bz(ixx, iyy, izz)
            jx_r(ix, iy, iz)=jx(ixxx, iyyy, izzz)
            jy_r(ix, iy, iz)=jy(ixxx, iyyy, izzz)
            jz_r(ix, iy, iz)=jz(ixxx, iyyy, izzz)
            rho_r(ix, iy, iz)=rho(ixxx, iyyy, izzz)
            rhoold_r(ix, iy, iz)=rhoold(ixxx, iyyy, izzz)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ELSE IF(absorbing_bcs) THEN

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, ixx, iyy , izz, ixxx, iyyy , izzz) COLLAPSE(3)
      DO iz=lbound_r(3),ubound_r(3)
        DO iy=lbound_r(2),ubound_r(2)
          DO ix=lbound_r(1),ubound_r(1)
            ixx = ix - lbound_r(1) + lbound_p(1)
            iyy = iy - lbound_r(2) + lbound_p(2)
            izz = iz - lbound_r(3) + lbound_p(3)
            ixxx = ix - lbound_r(1) + lbound_s(1)
            iyyy = iy - lbound_r(2) + lbound_s(2)
            izzz = iz - lbound_r(3) + lbound_s(3)
            exy_r(ix, iy, iz)=exy(ixx, iyy, izz)
            eyx_r(ix, iy, iz)=eyx(ixx, iyy, izz)
            ezx_r(ix, iy, iz)=ezx(ixx, iyy, izz)
            bxy_r(ix, iy, iz)=bxy(ixx, iyy, izz)
            byx_r(ix, iy, iz)=byx(ixx, iyy, izz)
            bzx_r(ix, iy, iz)=bzx(ixx, iyy, izz)
            exz_r(ix, iy, iz)=exz(ixx, iyy, izz)
            eyz_r(ix, iy, iz)=eyz(ixx, iyy, izz)
            ezy_r(ix, iy, iz)=ezy(ixx, iyy, izz)
            bxz_r(ix, iy, iz)=bxz(ixx, iyy, izz)
            byz_r(ix, iy, iz)=byz(ixx, iyy, izz)
            bzy_r(ix, iy, iz)=bzy(ixx, iyy, izz)
            jx_r(ix, iy, iz)=jx(ixxx, iyyy, izzz)
            jy_r(ix, iy, iz)=jy(ixxx, iyyy, izzz)
            jz_r(ix, iy, iz)=jz(ixxx, iyyy, izzz)
            rho_r(ix, iy, iz)=rho(ixxx, iyyy, izzz)
            rhoold_r(ix, iy, iz)=rhoold(ixxx, iyyy, izzz)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO  
    ENDIF
  END SUBROUTINE copy_field_forward

!  SUBROUTINE copy_field_backward_AM_rz()
!    USE fields, ONLY : el, er, et, bl, br, bt, jl, jr, jt
!    USE fields, ONLY : el_c, er_c, et_c, bl_c, br_c, bt_c, jl_c, jr_c, jt_c
!    USE omp_lib
!    USE picsar_precision, ONLY: idp
!
!    IMPLICIT NONE
!    INTEGER(idp) :: ir, il, imode, irr ,ill 
!    INTEGER(idp) , dimension(3) :: lbound_r, ubound_r, lbound_p ,ubound_p
!
!    lbound_r = LBOUND(el_c)
!    lbound_p = LBOUND(el)
!    ubound_r = UBOUND(el_c)
!    ubound_p = UBOUND(el)
!    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ir, il, imode, irr , ill) COLLAPSE(3)
!    DO imode=lbound_r(3),ubound_r(3)
!      DO il=lbound_r(2),ubound_r(2)
!        DO ir=lbound_r(1),ubound_r(1)
!          irr = ir - lbound_r(1) +lbound_p(1)
!          ill = il - lbound_r(2) +lbound_p(2)
!          el(irr, ill, imode)=el_c(ir, il, imode)
!          er(irr, ill, imode)=er_c(ir, il, imode)
!          et(irr, ill, imode)=et_c(ir, il, imode)
!          bl(irr, ill, imode)=bl_c(ir, il, imode)
!          br(irr, ill, imode)=br_c(ir, il, imode)
!          bt(irr, ill, imode)=bt_c(ir, il, imode)
!        END DO
!      END DO
!    END DO
!    !$OMP END PARALLEL DO
!  END SUBROUTINE copy_field_backward_AM_rz


  SUBROUTINE copy_field_backward()
    USE fields, ONLY: bx, bx_r, bxy, bxy_r, bxz, bxz_r, by, by_r, byx, byx_r, byz,   &
      byz_r, bz, bz_r, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exy, exy_r, exz, exz_r, ey, &
      ey_r, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezx, ezx_r, ezy, ezy_r, jx, jx_r, jy,  &
      jy_r, jz, jz_r
    USE omp_lib
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: absorbing_bcs, nx, ny, nz

    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz, ixx ,iyy , izz
    INTEGER(idp) , dimension(3) :: lbound_r, ubound_r, lbound_p ,ubound_p

    IF(absorbing_bcs) THEN
       lbound_r = LBOUND(exy_r)
       lbound_p = LBOUND(exy)
       ubound_r = UBOUND(exy_r)
       ubound_p = UBOUND(exy)
    ELSE
       lbound_r = LBOUND(ex_r)
       lbound_p = LBOUND(ex)
       ubound_r = UBOUND(ex_r)
       ubound_p = UBOUND(ex)
    ENDIF

    ! When using periodic bcs, standard EM fields are communicated 
    ! Else, when using absorbing bcs, splitted EM fields are communicated 
    IF(.NOT. absorbing_bcs) THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, ixx , iyy ,izz) COLLAPSE(3)
      DO iz=lbound_r(3),ubound_r(3)
        DO iy=lbound_r(2),ubound_r(2)
          DO ix=lbound_r(1),ubound_r(1)
            ixx = ix - lbound_r(1) +lbound_p(1)
            iyy = iy - lbound_r(2) +lbound_p(2)
            izz = iz - lbound_r(3) +lbound_p(3)
            ex(ixx, iyy, izz)=ex_r(ix, iy, iz)
            ey(ixx, iyy, izz)=ey_r(ix, iy, iz)
            ez(ixx, iyy, izz)=ez_r(ix, iy, iz)
            bx(ixx, iyy, izz)=bx_r(ix, iy, iz)
            by(ixx, iyy, izz)=by_r(ix, iy, iz)
            bz(ixx, iyy, izz)=bz_r(ix, iy, iz)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ELSE IF(absorbing_bcs) THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, ixx, iyy , izz) COLLAPSE(3)
      DO iz=lbound_r(3),ubound_r(3)
        DO iy=lbound_r(2),ubound_r(2)
          DO ix=lbound_r(1),ubound_r(1)
            ixx = ix - lbound_r(1) +lbound_p(1)
            iyy = iy - lbound_r(2) +lbound_p(2)
            izz = iz - lbound_r(3) +lbound_p(3)
            exy(ixx, iyy, izz)=exy_r(ix, iy, iz)
            eyx(ixx, iyy, izz)=eyx_r(ix, iy, iz)
            ezx(ixx, iyy, izz)=ezx_r(ix, iy, iz)
            bxy(ixx, iyy, izz)=bxy_r(ix, iy, iz)
            byx(ixx, iyy, izz)=byx_r(ix, iy, iz)
            bzx(ixx, iyy, izz)=bzx_r(ix, iy, iz)
            exz(ixx, iyy, izz)=exz_r(ix, iy, iz)
            eyz(ixx, iyy, izz)=eyz_r(ix, iy, iz)
            ezy(ixx, iyy, izz)=ezy_r(ix, iy, iz)
            bxz(ixx, iyy, izz)=bxz_r(ix, iy, iz)
            byz(ixx, iyy, izz)=byz_r(ix, iy, iz)
            bzy(ixx, iyy, izz)=bzy_r(ix, iy, iz)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ENDIF
  END SUBROUTINE copy_field_backward

#endif
END MODULE gpstd_solver

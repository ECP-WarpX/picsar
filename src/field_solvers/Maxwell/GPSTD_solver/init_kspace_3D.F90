!> @brief
!> This module contains subroutines for psatd_with_blocks
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________

MODULE gpstd_solver
  IMPLICIT NONE 
  CONTAINS 
!> @brief
!> This subroutine calculates spectral space and other usefull blocks
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________

  SUBROUTINE init_kspace(nx,ny,nz,dx,dy,dz,dt,porderx,pordery,porderz)
    USE matrix_coefficients
    USE CONSTANTS
    USE mpi_fftw3
    USE omp_lib
    USE fourier_psaotd , ONLY : FD_weights_hvincenti
    INTEGER(idp) , INTENT(IN)                     :: nx,ny,nz,porderx,pordery,porderz
    REAL(num)    , INTENT(IN)                     :: dx,dy,dz,dt
    LOGICAL(lp)                                   :: l_stg
    REAL(num)    , ALLOCATABLE , DIMENSION(:)     :: FD_x,FD_y,FD_z,onesx,onesy,onesz
    REAL(num)    , ALLOCATABLE , DIMENSION(:)     :: onesxp,onesyp,oneszp
    COMPLEX(cpx) , ALLOCATABLE , DIMENSION(:)     :: kxf,kyf,kzf,kxb,kyb,kzb,kxc,kyc,kzc,kxff,kxbb,kxcc
    REAL(num)    , ALLOCATABLE , DIMENSION(:,:,:) :: temp,temp2
    INTEGER(idp)                                  :: i,j,k
    COMPLEX(cpx)                                  :: ii
   
    IF(.NOT. ASSOCIATED(Kspace)) THEN
      ALLOCATE(KSPACE(ns_max))
    ENDIF
    nmatrixes2=nmatrixes2+1
    ALLOCATE(Kspace(nmatrixes2)%block_vector(10_idp))   ! 3 forward 3 backward
    IF(.NOT. ASSOCIATED(AT_OP)) THEN
      ALLOCATE(AT_OP(ns_max))
    ENDIF
    ALLOCATE(AT_OP(nmatrixes2)%block_vector(4_idp))  !S/k,C,(1-C)/k^2
    DO i = 1_idp , 10_idp
      ALLOCATE(Kspace(nmatrixes2)%block_vector(i)%block3dc(nx/2+1,ny,nz))
    ENDDO
    DO i = 1_idp , 4_idp
      ALLOCATE(AT_OP(nmatrixes2)%block_vector(i)%block3dc(nx/2+1,ny,nz))
    ENDDO
    !construct kspace
    l_stg = .TRUE.
    ii=(0.,1.)
    ALLOCATE(onesx(nx/2+1),onesxp(nx/2+1))
    ALLOCATE(onesy(ny),onesyp(ny))
    ALLOCATE(onesz(nz),oneszp(nz))
    DO i=1_idp,nx/2+1
      onesx(i)  = i-1_idp
      onesxp(i) = i-1_idp
    ENDDO
    DO j=1_idp,ny
      onesy(j)  = j-1_idp
      onesyp(j) = j-1_idp
      IF(j .GT. ny/2_idp +1) THEN
        onesy(j) = - onesy(j)
        onesyp(j) =  ny + onesyp(j)
      ENDIF
    ENDDO
    DO k=1_idp,nz
      onesz(k)  = k-1_idp
      oneszp(k) = k-1_idp
      IF(k .GT. nz/2_idp +1) THEN 
        onesz(k) = - onesz(k)
        oneszp(k) = ny + oneszp(k)
      ENDIF
    ENDDO
    IF (porderx .ne. 0_idp ) THEN  ! if 0 then infinite order
      ALLOCATE(FD_x(porderx/2))
      CALL FD_weights_hvincenti(porderx,FD_x,l_stg)
      ALLOCATE(kxf(nx/2+1),kxb(nx/2+1),kxc(nx/2+1))
      kxf=(0._num,0._num)*kxf
      kxb=(0._num,0._num)*kxb
      kxc=(0._num,0._num)*kxc
      DO i=1_idp,porderx/2
        kxc=kxc+2.0_num/dx*FD_x(i)*SIN((i*2.0_num-1.0_num)*PI*onesx*1.0_num/nx)
      ENDDO
    ELSE 
      ALLOCATE(kxff(nx),kxbb(nx),kxcc(nx))
      CALL fftfreq(nx,kxff,dx)
      CALL fftfreq(nx,kxbb,dx)
      CALL fftfreq(nx,kxcc,dx)
      kxf=kxff(1:nx/2+1)
      kxb=kxbb(1:nx/2+1) 
      kxc=kxcc(1:nx/2+1)
      DEALLOCATE(kxff,kxbb,kxcc)
    ENDIF
      kxf=kxc*EXP(-ii*PI*onesxp/nx)
      kxb=kxc*EXP(ii*PI*onesxp/nx)
    IF (pordery .ne. 0_idp) THEN 
      ALLOCATE(FD_y(pordery/2))
      CALL FD_weights_hvincenti(pordery,FD_y,l_stg)
      ALLOCATE(kyf(ny),kyb(ny),kyc(ny))
      kyf=(0._num,0._num)*kyf
      kyb=(0._num,0._num)*kyb
      kyc=(0._num,0._num)*kyc
      DO i=1_idp,pordery/2
        kyc=kyc+2.0_num/dy*FD_y(i)*SIN((i*2.0_num-1.0_num)*PI*onesy/ny)
      ENDDO
    ELSE
      CALL fftfreq(ny,kyf,dy)
      CALL fftfreq(ny,kyb,dy)
      CALL fftfreq(ny,kyc,dy)
    ENDIF
      kyf=kyc*EXP(-ii*PI*onesyp/ny)
      kyb=kyc*EXP(ii*PI*onesyp/ny)
    IF (porderz .ne. 0_idp) THEN
      ALLOCATE(FD_z(porderz/2))
      CALL FD_weights_hvincenti(porderz,FD_z,l_stg)
      ALLOCATE(kzf(nz),kzb(nz),kzc(nz))
      kzf=(0._num,0._num)*kzf
      kzb=(0._num,0._num)*kzb
      kzc=(0._num,0._num)*kzc
      DO i=1_idp,porderz/2
            kzc=kzc+2.0_num/dz*FD_z(i)*SIN((i*2.0_num-1.0_num)*PI*onesz/nz)
      ENDDO 
    ELSE
      CALL fftfreq(nz,kzf,dz)
      CALL fftfreq(nz,kzb,dz)
      CALL fftfreq(nz,kzc,dz)
    ENDIF
      kzf=kzc*EXP(-ii*PI*oneszp/nz)
      kzb=kzc*EXP(ii*PI*oneszp/nz)
    DO i = 1,nx/2+1
      DO j = 1,ny
        DO k = 1,nz
          Kspace(nmatrixes2)%block_vector(1)%block3dc(i,j,k) = kxf(i)
          Kspace(nmatrixes2)%block_vector(2)%block3dc(i,j,k) = kxb(i)
          Kspace(nmatrixes2)%block_vector(3)%block3dc(i,j,k) = kxc(i)
          Kspace(nmatrixes2)%block_vector(4)%block3dc(i,j,k) = kyf(j)
          Kspace(nmatrixes2)%block_vector(5)%block3dc(i,j,k) = kyb(j)
          Kspace(nmatrixes2)%block_vector(6)%block3dc(i,j,k) = kyc(j)
          Kspace(nmatrixes2)%block_vector(7)%block3dc(i,j,k) = kzf(k)
          Kspace(nmatrixes2)%block_vector(8)%block3dc(i,j,k) = kzb(k)
          Kspace(nmatrixes2)%block_vector(9)%block3dc(i,j,k) = kzc(k)
          Kspace(nmatrixes2)%block_vector(10)%block3dc(i,j,k) = SQRT(kxc(i)**2+kyc(j)**2+kzc(k)**2) 
        ENDDO
      ENDDO
    ENDDO
    
    ALLOCATE(temp(nx/2+1,ny,nz))
    ALLOCATE(temp2(nx/2+1,ny,nz))
    temp=dt*clight*REAL(Kspace(nmatrixes2)%block_vector(10)%block3dc,num)
    temp2=sinc_block(nx/2+1_idp,ny,nz,temp) 
    AT_OP(nmatrixes2)%block_vector(1)%block3dc = CMPLX(temp2,0._num)
    AT_OP(nmatrixes2)%block_vector(1)%block3dc = clight*dt*AT_OP(nmatrixes2)%block_vector(1)%block3dc
    temp2=COS(temp)
    AT_OP(nmatrixes2)%block_vector(2)%block3dc = CMPLX(temp2,0._num)
    temp=0.5_num*temp
    AT_OP(nmatrixes2)%block_vector(3)%block3dc = 2._num*(clight*dt/2.)**2&
    *sinc_block(nx/2+1,ny,nz,temp)*sinc_block(nx/2+1,ny,nz,temp)
   
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1) = (1.,0.)
  
    AT_OP(nmatrixes2)%block_vector(3)%block3dc = ((1.,0.) -  AT_OP(nmatrixes2)%block_vector(2)%block3dc) &
    /Kspace(nmatrixes2)%block_vector(10)%block3dc**2
    AT_OP(nmatrixes2)%block_vector(3)%block3dc(1,1,1) = (clight*dt)**2/2.0_num 
      
  
  
             !(1-C)/k^2
    temp=2._num*temp
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)=(1.0_num,0.0_num)
    AT_OP(nmatrixes2)%block_vector(4)%block3dc = (AT_OP(nmatrixes2)%block_vector(1)%block3dc-clight*dt)&
    / Kspace(nmatrixes2)%block_vector(10)%block3dc/ Kspace(nmatrixes2)%block_vector(10)%block3dc
    AT_OP(nmatrixes2)%block_vector(4)%block3dc(1,1,1)=CMPLX(-(clight*dt)**3/6.0_num,0.0_num)  
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)=(0._num,0._num)
    DEALLOCATE(kxf,kxb,kxc,kyf,kyb,kyc,kzf,kzb,kzc,temp,temp2,onesx,onesy,onesz) 
  END SUBROUTINE init_kspace
  
  SUBROUTINE fftfreq(nxx,kxx, dxx)
  USE constants
    IMPLICIT NONE
    INTEGER(idp)  , INTENT(IN)                    :: nxx
    REAL(num)     , INTENT(IN)                    :: dxx
    COMPLEX(cpx)  , INTENT(OUT) , DIMENSION(nxx)  :: kxx
    INTEGER(idp) :: i,n
    REAL(num) :: fe
  
    fe=1_num/dxx
    n=nxx
    kxx(1)=(0.,0.)
    IF (MOD(n,2) .EQ. 0)THEN
    ! First part of k [0,...,n/2-1]
      DO i=1,n/2_idp-1_idp
         kxx(i+1)=kxx(i)+(1.,0.)
      END DO
    ! Second part of k [-n/2,-1]
      kxx(n/2_idp+1)=-n/2_idp
      DO i=n/2_idp+1,n-1
        kxx(i+1)=kxx(i)+(1.,0.)
      END DO
    ELSE
    ! First part of k [0,...,(n-1)/2]
      DO i=1,(n-1_idp)/2_idp
        kxx(i+1)=kxx(i)+(1.,0.)
      END DO
    ! Second part of k [-(n-1)/2,-1]
      kxx((n-1_idp)/2_idp+2_idp)=-cmplx((n-1_idp)/2_idp,0.)
      DO i=(n-1_idp)/2_idp+2_idp,n-1
        kxx(i+1)=kxx(i)+(1.0_num,0.0_num)
      END DO
    ENDIF
    kxx=kxx/(dxx*nxx)/2.0_num*PI
  END SUBROUTINE fftfreq
  
  FUNCTION sinc_block(n1,n2,n3,block)
    USE  constants 
    INTEGER(idp), INTENT(IN)                     :: n1,n2,n3
    REAL(num)   , DIMENSION(:,:,:) , INTENT(in)  :: block
    REAL(num)   , DIMENSION(:,:,:) , ALLOCATABLE :: sinc_block
    INTEGER(idp)       :: i,j,k
    ALLOCATE(sinc_block(n1,n2,n3))
    DO k=1,n3
      DO j = 1,n2
        DO i = 1 , n1
          sinc_block(i,j,k)=sinc(block(i,j,k))
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END FUNCTION sinc_block
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
!> @brief
!> This subroutine inits fftw plans for gpstd
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________

 
  SUBROUTINE init_plans_gpstd(n1,n2,n3)
  USE fields
  USE matrix_coefficients
  USE shared_data
  USE fastfft
  Use fourier
  USE fftw3_fortran
  
#ifdef _OPENMP
   USE omp_lib
#endif
  INTEGER(idp), INTENT(in):: n1,n2,n3
  INTEGER(idp)           ::  nopenmp
#ifdef _OPENMP
    nopenmp=OMP_GET_MAX_THREADS()
    CALL OMP_SET_NESTED(.TRUE.)
#else
    nopenmp=1_idp
#endif
  CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp,n1,n2,n3,ex_r,vold(1)%block_vector(1)%block3dc &
       ,plan_r2c,INT(FFTW_MEASURE,idp),INT(FFTW_FORWARD,idp))
  CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp,n1,n2,n3,vold(1)%block_vector(1)%block3dc,&
      ex_r,plan_c2r,INT(FFTW_MEASURE,idp),INT(FFTW_BACKWARD,idp))
  END SUBROUTINE

!> @brief
!> This subroutine destroys fftw plans at the end of the program
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________

  SUBROUTINE destroy_plans_gpstd()
    USE  fastfft
    USE fourier
    CALL fast_fftw_destroy_plan_dft(plan_r2c)
    CALL fast_fftw_destroy_plan_dft(plan_c2r)
  END SUBROUTINE

!> @brief
!> This subroutine frees memory from useless blocks
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________


  SUBROUTINE delete_arrays
    USE matrix_coefficients
    IMPLICIT NONE
    INTEGER(idp)     ::  i,j
    LOGICAL(lp)      :: needed
    DO i = 1_idp,10_idp
       DEALLOCATE(Kspace(nmatrixes2)%block_vector(i)%block3dc)
    ENDDO
    DO i = 1_idp,4_idp
       DEALLOCATE(AT_OP(nmatrixes2)%block_vector(i)%block3dc)
    ENDDO
    DEALLOCATE(Kspace(nmatrixes2)%block_vector)
    DEALLOCATE(AT_OP(nmatrixes2)%block_vector)
    DO i = 1,11
      DO j=1,11
        CALL is_calculation_needed(i,j,needed)
        IF(needed .EQV. .FALSE.) THEN 
          IF(associated(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc)) THEN
            DEALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE
!> @brief
!> This subroutine constructs gpstd_blocks and puts them into cc_mat operator
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________


  SUBROUTINE init_gpstd(nx,ny,nz,dx,dy,dz,dt,porderx,pordery,porderz)
    USE matrix_coefficients
    USE PICSAR_PRECISION
    USE CONSTANTS
    USE mpi_fftw3
    USE omp_lib
    INTEGER(idp),INTENT(IN):: nx,ny,nz,porderx,pordery,porderz
    REAL(num),INTENT(IN)   :: dx,dy,dz,dt
    INTEGER(idp)           :: i,j,k,p
    COMPLEX(cpx)           :: ii
    LOGICAL(lp)            :: needed
    ii=(0.,1.)
    CALL allocate_new_matrix_vector(11_idp)
    CALL init_kspace(nx,ny,nz,dx,dy,dz,dt,porderx,pordery,porderz)
    DO i=1_idp,6_idp
      DO j=1_idp,11_idp
        CALL is_calculation_needed(i,j,needed) 
        IF(.NOT. needed) CYCLE
        ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc(nx/2+1,ny,nz))
        cc_mat(nmatrixes)%block_matrix2d(i,j)%nx = nx/2+1
        cc_mat(nmatrixes)%block_matrix2d(i,j)%ny = ny
        cc_mat(nmatrixes)%block_matrix2d(i,j)%nz = nz
      ENDDO
    ENDDO
    DO i=1_idp,11_idp
      ALLOCATE(vold(nmatrixes)%block_vector(i)%block3dc(nx/2+1,ny,nz))
      !vold(nmatrixes)%block_vector(i)%block3dc = CMPLX(0.,0.)
      vold(nmatrixes)%block_vector(i)%nx = nx/2+1
      vold(nmatrixes)%block_vector(i)%ny = ny 
      vold(nmatrixes)%block_vector(i)%nz = nz
    ENDDO
    DO i=1_idp,6_idp
      ALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc(nx/2+1,ny,nz))
      vnew(nmatrixes)%block_vector(i)%nx = nx/2+1
      vnew(nmatrixes)%block_vector(i)%ny = ny
      vnew(nmatrixes)%block_vector(i)%nz = nz
    ENDDO
    cc_mat(nmatrixes)%block_matrix2d(1,5)%block3dc =         &
    - ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight  &
      *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(1,6)%block3dc =         &
      ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight  &
      *AT_OP(nmatrixes2)%block_vector(1)%block3dc
    cc_mat(nmatrixes)%block_matrix2d(2,4)%block3dc =         &
      ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight  &
      *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(2,6)%block3dc =         &
     -ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight  &
      *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3,4)%block3dc =         &
    - ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight  &
      *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3,5)%block3dc =         &
      ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight  &
      *AT_OP(nmatrixes2)%block_vector(1)%block3dc
 
    cc_mat(nmatrixes)%block_matrix2d(4,2)%block3dc =         &
      ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight &
       *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(4,3)%block3dc =         &
      -ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight &
       *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5,1)%block3dc =         &      
      -ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight &
       *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5,3)%block3dc =         &
      ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight &
       *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6,1)%block3dc =         & 
      ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight &
       *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6,2)%block3dc =         &
      -ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight &
       *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    DO i=1,6
       cc_mat(nmatrixes)%block_matrix2d(i,i)%block3dc =        &
          AT_OP(nmatrixes2)%block_vector(2)%block3dc
    ENDDO
    DO i = 1 , 3 
       cc_mat(nmatrixes)%block_matrix2d(i,i+6)%block3dc =       &
         (-1._num)*clight*mu0*AT_OP(nmatrixes2)%block_vector(1)%block3dc
    ENDDO
    cc_mat(nmatrixes)%block_matrix2d(4,8)%block3dc = &
    - mu0* ii*Kspace(nmatrixes2)%block_vector(8)%block3dc*AT_OP(nmatrixes2)%block_vector(3)%block3dc
    cc_mat(nmatrixes)%block_matrix2d(4,9)%block3dc = &
    - mu0*(-ii)*Kspace(nmatrixes2)%block_vector(5)%block3dc*AT_OP(nmatrixes2)%block_vector(3)%block3dc
    cc_mat(nmatrixes)%block_matrix2d(5,7)%block3dc = &
    - mu0*(-ii)*Kspace(nmatrixes2)%block_vector(8)%block3dc*AT_OP(nmatrixes2)%block_vector(3)%block3dc
    cc_mat(nmatrixes)%block_matrix2d(5,9)%block3dc = &
    - mu0* ii*Kspace(nmatrixes2)%block_vector(2)%block3dc*AT_OP(nmatrixes2)%block_vector(3)%block3dc
     cc_mat(nmatrixes)%block_matrix2d(6,7)%block3dc = &
        - mu0* ii*Kspace(nmatrixes2)%block_vector(5)%block3dc*AT_OP(nmatrixes2)%block_vector(3)%block3dc
     cc_mat(nmatrixes)%block_matrix2d(6,8)%block3dc = &
       -mu0*(-ii)*Kspace(nmatrixes2)%block_vector(2)%block3dc*AT_OP(nmatrixes2)%block_vector(3)%block3dc
    !contribution rho old ( vold(10))
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1) = (1.0_num,0.0_num)
    DO i = 1,3
      cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc =        &
        CMPLX(0.,1.)*(AT_OP(nmatrixes2)%block_vector(2)%block3dc   &
        -1./(clight*dt)*AT_OP(nmatrixes2)%block_vector(1)%block3dc)&
        /Kspace(nmatrixes2)%block_vector(10)%block3dc**2
  
      cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc =        &
         cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc &
         *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
  
      cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc(1,1,1) = &
         -1.0_num/3.0_num*(0.0_num,1.0_num)*(clight*dt)**2

      cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc = 1./eps0 &
        *cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc
    ENDDO
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1) = (0.,0.)
    !contribution rho new (vold(11)) 
    DO i = 1 , 3
      cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc =  &
      CMPLX(0.,1.)*(1./(clight*dt)*                          &
      AT_OP(nmatrixes2)%block_vector(1)%block3dc             &
      -CMPLX(1.,0.))/Kspace(nmatrixes2)%block_vector(10)%block3dc**2
  
      cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc = &
         cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc &
         *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
      cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc(1,1,1) = &
       -1.0_num/6.0_num*(0.0_num,1.0_num)*(clight*dt)**2  

       cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc = 1./eps0 &
       *  cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc     
    ENDDO
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)   = (0.,0.)
    !CALL delete_arrays 
  END SUBROUTINE init_gpstd
!> @brief
!> This subroutine executes forward fftw on relevent fields
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________



  SUBROUTINE execute_fftw_gpstd_r2c
    USE fields
    USE matrix_coefficients
    USE shared_data
    USE time_stat
    USE params      
    USE fourier_psaotd
    USE fourier
    USE fastfft  
    INTEGER(idp) :: nfftx,nffty,nfftz, nxx,nyy,nzz
    REAL(num)    :: tmptime
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nz+2*nzguards
    nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;
    
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
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
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ex_r,vold(nmatrixes)%block_vector(1)%block3dc, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ey_r,vold(nmatrixes)%block_vector(2)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,ez_r,vold(nmatrixes)%block_vector(3)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bx_r,vold(nmatrixes)%block_vector(4)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,by_r,vold(nmatrixes)%block_vector(5)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,bz_r,vold(nmatrixes)%block_vector(6)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jx_r,vold(nmatrixes)%block_vector(7)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jy_r,vold(nmatrixes)%block_vector(8)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,jz_r,vold(nmatrixes)%block_vector(9)%block3dc , plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rhoold_r,vold(nmatrixes)%block_vector(10)%block3dc, plan_r2c)
    CALL fast_fftw3d_r2c_with_plan(nfftx,nffty,nfftz,rho_r,vold(nmatrixes)%block_vector(11)%block3dc , plan_r2c)
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE


!> @brief
!> This subroutine executes backward fftw on relevent fields
!
!> @author
!> Haithem Kallala
!
!> @date
!> Creation 2017
! ________________________________________________________________________________________

  SUBROUTINE execute_fftw_gpstd_c2r
    USE fields
    USE matrix_coefficients
    USE shared_data
    USE time_stat
    USE params
    USE fourier_psaotd
    Use fourier
    USE fastfft
    INTEGER(idp) :: nfftx,nffty,nfftz, nxx,nyy,nzz
    REAL(num)    :: tmptime,coeff_norm
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nz+2*nzguards
    nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,vnew(nmatrixes)%block_vector(1)%block3dc,ex_r,plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,vnew(nmatrixes)%block_vector(2)%block3dc,ey_r,plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,vnew(nmatrixes)%block_vector(3)%block3dc,ez_r,plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,vnew(nmatrixes)%block_vector(4)%block3dc,bx_r,plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,vnew(nmatrixes)%block_vector(5)%block3dc,by_r,plan_c2r)
    CALL fast_fftw3d_c2r_with_plan(nfftx,nffty,nfftz,vnew(nmatrixes)%block_vector(6)%block3dc,bz_r,plan_c2r)
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF
    coeff_norm = 1._num/(nfftx*nffty*nfftz)
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
 
    CALL normalize_Fourier(ex,nxx,nyy,nzz,ex_r,nfftx,nffty,nfftz,coeff_norm)
    CALL normalize_Fourier(ey,nxx,nyy,nzz,ey_r,nfftx,nffty,nfftz,coeff_norm)
    CALL normalize_Fourier(ez,nxx,nyy,nzz,ez_r,nfftx,nffty,nfftz,coeff_norm)
    CALL normalize_Fourier(bx,nxx,nyy,nzz,bx_r,nfftx,nffty,nfftz,coeff_norm)
    CALL normalize_Fourier(by,nxx,nyy,nzz,by_r,nfftx,nffty,nfftz,coeff_norm)
    CALL normalize_Fourier(bz,nxx,nyy,nzz,bz_r,nfftx,nffty,nfftz,coeff_norm)
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
 
  
  END SUBROUTINE 
END MODULE














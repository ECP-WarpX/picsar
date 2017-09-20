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
  USE PICSAR_PRECISION
  IMPLICIT NONE 
  COMPLEX(cpx), DIMENSION(:)      , ALLOCATABLE :: kxc,kxb,kxf,kyc,kyb,kyf,kzc,kzb,kzf
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
#if defined(FFTW)
  SUBROUTINE select_case_dims_local(nfftx,nffty,nfftz)
    USE shared_data
    USE mpi_fftw3
    USE group_parameters
    USE params
    USE fields      , ONLY : nxguards,nyguards,nzguards

    INTEGER(idp)  , INTENT(INOUT) :: nfftx,nffty,nfftz
    IF( fftw_with_mpi) THEN
      IF(.NOT. fftw_hybrid) THEN
        nfftx=nx_global
        nffty=ny_global
        nfftz=local_nz
      ELSE
        nfftx = nx_group
        nffty = ny_group
        nfftz = local_nz
      ENDIF
    ELSE
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      nfftz = nz+2*nzguards
    ENDIF
  END SUBROUTINE
  SUBROUTINE select_case_dims_global(nfftx,nffty,nfftz)
    USE shared_data
    USE mpi_fftw3
    USE group_parameters
    USE params
    USE fields      , ONLY : nxguards,nyguards,nzguards

    INTEGER(idp)  , INTENT(INOUT) :: nfftx,nffty,nfftz
    IF( fftw_with_mpi) THEN
      IF(.NOT. fftw_hybrid) THEN
        nfftx=nx_global
        nffty=ny_global
        nfftz=nz_global
      ELSE
        nfftx = nx_group
        nffty = ny_group
        nfftz = nz_group
      ENDIF
    ELSE
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      nfftz = nz+2*nzguards
    ENDIF
  END SUBROUTINE



  SUBROUTINE init_kspace
    USE matrix_coefficients
    USE CONSTANTS
    USE mpi_fftw3
    USE omp_lib
    USE shared_data !:, ONLY : dx,dy,dz,nx,ny,nz
    USE fields      , ONLY : nxguards,nyguards,nzguards
    USE fields , ONLY : norderx,nordery,norderz
    USE params , ONLY : dt
    USE mpi_fftw3
    LOGICAL(lp)                                   :: l_stg
    REAL(num)    , ALLOCATABLE , DIMENSION(:,:,:) :: temp,temp2
    INTEGER(idp)                                  :: i,j,k
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: nfftx,nffty,nfftz
    LOGICAL(lp)                                   :: switch
    IF(.NOT. ASSOCIATED(Kspace)) THEN
      ALLOCATE(KSPACE(ns_max))
    ENDIF
    IF(.NOT. ASSOCIATED(AT_OP)) THEN
      ALLOCATE(AT_OP(ns_max))
    ENDIF
    nmatrixes2=nmatrixes2+1
    ALLOCATE(Kspace(nmatrixes2)%block_vector(10_idp))   ! 3 forward 3 backward
    IF(.NOT. ASSOCIATED(AT_OP)) THEN
      ALLOCATE(AT_OP(ns_max))
    ENDIF
    ALLOCATE(AT_OP(nmatrixes2)%block_vector(4_idp))  !S/k,C,(1-C)/k^2
  !  IF( fftw_with_mpi) THEN
  !    IF(.NOT. fftw_mpi_transpose) THEN
  !      nfftx=nx_global
  !      nffty=ny_global
  !      nfftz=local_nz
  !    ELSE
  !      nfftx = nx_global
  !      nffty = local_ny
  !      nfftz = nz_global
  !    ENDIF
  !  ENDIF
  !  IF(.NOT. fftw_with_mpi) THEN
  !    nfftx = nx+2*nxguards
  !    nffty = ny+2*nyguards
  !    nfftz = nz+2*nzguards
  !  ENDIF
  !  IF(fftw_hybrid) THEN 
  !     nfftx = nx_group_grid
  !     nffty = ny_group_grid 
  !     nfftz = local_nz
  !  ENDIF
    CALL select_case_dims_local(nfftx,nffty,nfftz)
    DO i = 1_idp , 10_idp
      ALLOCATE(Kspace(nmatrixes2)%block_vector(i)%block3dc(nfftx/2+1,nffty,nfftz))
    ENDDO
    DO i = 1_idp , 4_idp
      ALLOCATE(AT_OP(nmatrixes2)%block_vector(i)%block3dc(nfftx/2+1,nffty,nfftz))
    ENDDO
    !construct kspace
    l_stg = .TRUE.
    ii=DCMPLX(0.0_num,1.0_num)
    CALL compute_k_vec_nompi(l_stg)
    DO i = 1,nfftx/2+1
      DO j = 1,nffty
        DO k = 1,nfftz
          Kspace(nmatrixes2)%block_vector(1)%block3dc(i,j,k) = kxf(i)
          Kspace(nmatrixes2)%block_vector(2)%block3dc(i,j,k) = kxb(i)
          Kspace(nmatrixes2)%block_vector(3)%block3dc(i,j,k) = kxc(i)
          Kspace(nmatrixes2)%block_vector(4)%block3dc(i,j,k) = kyf(j)
          Kspace(nmatrixes2)%block_vector(5)%block3dc(i,j,k) = kyb(j)
          Kspace(nmatrixes2)%block_vector(6)%block3dc(i,j,k) = kyc(j)
          Kspace(nmatrixes2)%block_vector(7)%block3dc(i,j,k) = kzf(k)
          Kspace(nmatrixes2)%block_vector(8)%block3dc(i,j,k) = kzb(k)
          Kspace(nmatrixes2)%block_vector(9)%block3dc(i,j,k) = kzc(k)
          Kspace(nmatrixes2)%block_vector(10)%block3dc(i,j,k) = SQRT((kxc(i)**2+kyc(j)**2+kzc(k)**2))
        ENDDO
      ENDDO
    ENDDO
    switch = .FALSE.
    ALLOCATE(temp(nfftx/2+1,nffty,nfftz))
    ALLOCATE(temp2(nfftx/2+1,nffty,nfftz))
    temp=dt*clight*REAL(Kspace(nmatrixes2)%block_vector(10)%block3dc,num)
    temp2=sinc_block(nfftx/2+1_idp,nffty,nfftz,temp) 
    AT_OP(nmatrixes2)%block_vector(1)%block3dc = DCMPLX(temp2,0._num)
    AT_OP(nmatrixes2)%block_vector(1)%block3dc = clight*dt*AT_OP(nmatrixes2)%block_vector(1)%block3dc
    temp2=COS(temp)
    AT_OP(nmatrixes2)%block_vector(2)%block3dc = DCMPLX(temp2,0._num)
    temp=0.5_num*temp
    AT_OP(nmatrixes2)%block_vector(3)%block3dc = 2._num*(clight*dt/2.0_num)**2&
    *sinc_block(nfftx/2+1,nffty,nfftz,temp)*sinc_block(nfftx/2+1,nffty,nfftz,temp)

    IF(ABS(Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)) .EQ. 0.0_num) THEN 
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1) = DCMPLX(1.0_num,0.0_num)
      switch = .TRUE.
    ENDIF
      AT_OP(nmatrixes2)%block_vector(3)%block3dc = (DCMPLX(1.0_num,0.0_num) -  AT_OP(nmatrixes2)%block_vector(2)%block3dc) &
      /Kspace(nmatrixes2)%block_vector(10)%block3dc**2
               !(1-C)/k^2
      AT_OP(nmatrixes2)%block_vector(4)%block3dc = (AT_OP(nmatrixes2)%block_vector(1)%block3dc-clight*dt)&
      / Kspace(nmatrixes2)%block_vector(10)%block3dc/ Kspace(nmatrixes2)%block_vector(10)%block3dc
    IF(switch) THEN
      AT_OP(nmatrixes2)%block_vector(3)%block3dc(1,1,1) = (clight*dt)**2/2.0_num
      AT_OP(nmatrixes2)%block_vector(4)%block3dc(1,1,1)=DCMPLX(-(clight*dt)**3/6.0_num,0.0_num)  
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)=DCMPLX(0._num,0._num)
    ENDIF
    DEALLOCATE(temp,temp2) 
  END SUBROUTINE init_kspace

  SUBROUTINE compute_k_vec_nompi(l_stg)
    USE constants
    USE shared_data !, ONLY : dx,dy,dz,nx,ny,nz
    USE fields 
    USE mpi_fftw3
    IMPLICIT NONE 
    LOGICAL(lp)  , INTENT(IN)  :: l_stg
    COMPLEX(cpx) , ALLOCATABLE , DIMENSION(:)     :: kxff,kxbb,kxcc,k_temp
    REAL(num)    , ALLOCATABLE , DIMENSION(:)     :: FD_x,FD_y,FD_z
    COMPLEX(cpx) , ALLOCATABLE , DIMENSION(:)     :: onesxp,onesyp,oneszp,onesx,onesy,onesz
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: i,j,k
    INTEGER(idp)                                  :: nfftx,nffty,nfftz
    ii = DCMPLX(0.0_num,1.0_num)
    !IF(.NOT. fftw_with_mpi) THEN
    !  nfftx=nx+2*nxguards
    !  nffty=ny+2*nyguards
    !  nfftz=nz+2*nzguards
    !ENDIF
    !IF(fftw_with_mpi) THEN
    !  IF(.NOT. fftw_mpi_transpose) THEN
    !    nfftx = nx_global 
    !    nffty = ny_global
    !    nfftz = nz_global 
    !  ELSE
    !    nfftx = nx_global
    !    nffty = ny_global
    !    nfftz = nz_global
    !  ENDIF
    !ENDIF
    CALL select_case_dims_global(nfftx,nffty,nfftz)
    ALLOCATE(onesx(nfftx/2+1),onesxp(nfftx/2+1))
    ALLOCATE(onesy(nffty),onesyp(nffty))
    ALLOCATE(onesz(nfftz),oneszp(nfftz))
    DO i=1_idp,nfftx/2+1
      onesx(i)  = DCMPLX(i-1.0_num,0.0_num)
      onesxp(i) = DCMPLX(i-1.0_num,0.0_num)
    ENDDO
    DO j=1_idp,nffty
      onesy(j)  = DCMPLX(j-1.0_num,0.0_num)
      onesyp(j) = DCMPLX(j-1.0_num,0.0_num)
      IF(j .GT. nffty/2_idp +1) THEN
        onesy(j)  =DCMPLX(-onesy(j))
        onesyp(j) =DCMPLX( nffty + onesyp(j))
      ENDIF
    ENDDO
    DO k=1_idp,nfftz
      onesz(k)  = DCMPLX(k-1.0_num)
      oneszp(k) = DCMPLX(k-1.0_num)
      IF(k .GT. nfftz/2_idp +1) THEN
        onesz(k) = DCMPLX(- onesz(k))
        oneszp(k) = DCMPLX(nfftz + oneszp(k))
      ENDIF
    ENDDO
    IF(.NOT. ALLOCATED(kxf)) THEN
      ALLOCATE(kxf(nfftx/2+1),kxb(nfftx/2+1),kxc(nfftx/2+1))
      ALLOCATE(kyf(nffty),kyb(nffty),kyc(nffty))
      ALLOCATE(kzf(nfftz),kzb(nfftz),kzc(nfftz))
    ELSE
      DEALLOCATE(kxf,kyf,kzf,kxb,kyb,kzb,kxc,kyc,kzc)
      ALLOCATE(kxf(nfftx/2+1),kxb(nfftx/2+1),kxc(nfftx/2+1))
      ALLOCATE(kyf(nffty),kyb(nffty),kyc(nffty))
      ALLOCATE(kzf(nfftz),kzb(nfftz),kzc(nfftz))
    ENDIF
    IF (norderx .ne. 0_idp ) THEN  ! if 0 then infinite order
      ALLOCATE(FD_x(norderx/2))
      CALL FD_weights_hvincenti(norderx,FD_x,l_stg)
      kxf=(0._num,0._num)*kxf
      kxb=(0._num,0._num)*kxb
      kxc=(0._num,0._num)*kxc
      DO i=1_idp,norderx/2
        kxc=kxc+2.0_num/dx*FD_x(i)*SIN((i*2.0_num-1.0_num)*PI*onesx/nfftx)
      ENDDO
      DEALLOCATE(FD_x)
    ELSE
      ALLOCATE(kxff(nfftx),kxbb(nfftx),kxcc(nfftx))
      CALL fftfreq(nfftx,kxff,dx)
      CALL fftfreq(nfftx,kxbb,dx)
      CALL fftfreq(nfftx,kxcc,dx)
      kxf=kxff(1:nfftx/2+1)
      kxb=kxbb(1:nfftx/2+1)
      kxc=kxcc(1:nfftx/2+1)
      DEALLOCATE(kxff,kxbb,kxcc)
    ENDIF
    IF (nordery .ne. 0_idp) THEN
      ALLOCATE(FD_y(nordery/2))
      CALL FD_weights_hvincenti(nordery,FD_y,l_stg)
      kyf=(0._num,0._num)*kyf
      kyb=(0._num,0._num)*kyb
      kyc=(0._num,0._num)*kyc
      DO i=1_idp,nordery/2
        kyc=kyc+2.0_num/dy*FD_y(i)*SIN((i*2.0_num-1.0_num)*PI*onesy/nffty)
      ENDDO
      DEALLOCATE(FD_y)
    ELSE
      CALL fftfreq(nffty,kyf,dy)
      CALL fftfreq(nffty,kyb,dy)
      CALL fftfreq(nffty,kyc,dy)
    ENDIF
    IF (norderz .ne. 0_idp) THEN
      ALLOCATE(FD_z(norderz/2))
      CALL FD_weights_hvincenti(norderz,FD_z,l_stg)
      kzf=(0._num,0._num)*kzf
      kzb=(0._num,0._num)*kzb
      kzc=(0._num,0._num)*kzc
      DO i=1_idp,norderz/2
        kzc=kzc+2.0_num/dz*FD_z(i)*SIN((i*2.0_num-1.0_num)*PI*onesz/nfftz)
      ENDDO
      DEALLOCATE(FD_z)
    ELSE
      CALL fftfreq(nfftz,kzf,dz)
      CALL fftfreq(nfftz,kzb,dz)
      CALL fftfreq(nfftz,kzc,dz)
    ENDIF
    IF(l_stg) THEN
    kxf=kxc*EXP(-ii*PI*onesxp/nfftx)
    kxb=kxc*EXP(ii*PI*onesxp/nfftx)

    kyf=kyc*EXP(-ii*PI*onesyp/nffty)
    kyb=kyc*EXP(ii*PI*onesyp/nffty)

    kzf=kzc*EXP(-ii*PI*oneszp/nfftz)
    kzb=kzc*EXP(ii*PI*oneszp/nfftz)
    ELSE
      kxf=kxc
      kxb=kxc
      kyf=kyc
      kzb=kzc
    ENDIF
    IF(fftw_with_mpi) THEN
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
      ELSE 
        ALLOCATE(k_temp(nffty))
        k_temp = kyc
        DEALLOCATE(kyc);ALLOCATE(kyc(local_ny))
        kyc = k_temp(local_y0+1:local_y0+local_ny)
        k_temp = kyf
        DEALLOCATE(kyf);ALLOCATE(kyf(local_ny))
        kyf = k_temp(local_y0+1:local_y0+local_ny)
        k_temp = kyb
        DEALLOCATE(kyb);ALLOCATE(kyb(local_ny))
        kyb = k_temp(local_y0+1:local_y0+local_ny)
      ENDIF
      DEALLOCATE(k_temp)
    ENDIF
    DEALLOCATE(onesx,onesy,onesz,onesxp,onesyp,oneszp)
  END SUBROUTINE
  SUBROUTINE fftfreq(nxx,kxx, dxx)
  USE constants
    IMPLICIT NONE
    INTEGER(idp)  , INTENT(IN)                    :: nxx
    REAL(num)     , INTENT(IN)                    :: dxx
    COMPLEX(cpx)  , INTENT(OUT) , DIMENSION(:)  :: kxx
    INTEGER(idp) :: i,n
    REAL(num) :: fe
  
    fe=1_num/dxx
    n=nxx
    kxx(1)=dcmplx(0.,0.)
    IF (MOD(n,2) .EQ. 0)THEN
    ! First part of k [0,...,n/2-1]
      DO i=1,n/2_idp-1_idp
         kxx(i+1)=kxx(i)+(1.,0.)
      END DO
    ! Second part of k [-n/2,-1]
      kxx(n/2_idp+1)=-n/2_idp
      DO i=n/2_idp+1,n-1
        kxx(i+1)=kxx(i)+dcmplx(1.,0.)
      END DO
    ELSE
    ! First part of k [0,...,(n-1)/2]
      DO i=1,(n-1_idp)/2_idp
        kxx(i+1)=kxx(i)+(1.,0.)
      END DO
    ! Second part of k [-(n-1)/2,-1]
      kxx((n-1_idp)/2_idp+2_idp)=-dcmplx((n-1_idp)/2_idp,0.)
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

 
  SUBROUTINE init_plans_gpstd
  USE fields
  USE matrix_coefficients
  USE shared_data
  USE fastfft
  Use fourier, ONLY: plan_r2c,plan_c2r
!  USE fftw3_fortran
  USE mpi_fftw3 
#ifdef _OPENMP
   USE omp_lib
#endif
  INTEGER(idp)           :: n1,n2,n3
  INTEGER(idp)           ::  nopenmp
  INTEGER(C_INT) :: nopenmp_cint
  INTEGER(C_INTPTR_T) :: nx_cint, ny_cint, nz_cint

#ifdef _OPENMP
    nopenmp=OMP_GET_MAX_THREADS()
#else
    nopenmp=1_idp
#endif
  nopenmp_cint=nopenmp
  IF(fftw_threads_ok) THEN
    CALL  DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
  ENDIF
  IF(.NOT. fftw_with_mpi) THEN 
    n1=nx+2*nxguards
    n2=ny+2*nyguards
    n3=nz+2*nzguards
    CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp,n1,n2,n3,ex_r,vold(1)%block_vector(1)%block3dc &
       ,plan_r2c,INT(FFTW_MEASURE,idp),INT(FFTW_FORWARD,idp))
    CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp,n1,n2,n3,vnew(1)%block_vector(1)%block3dc,&
      ex_r,plan_c2r,INT(FFTW_MEASURE,idp),INT(FFTW_BACKWARD,idp))
  ENDIF
  IF(fftw_with_mpi) THEN
     nx_cint = nx_global
     ny_cint = ny_global
     nz_cint = nz_global 
     IF(.NOT. fftw_mpi_transpose) THEN
       plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint,ny_cint,nx_cint, &
                      ex_r, vold(nmatrixes)%block_vector(1)%block3dc, comm, FFTW_MEASURE)
       plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(nz_cint,ny_cint,nx_cint, &
                      vnew(nmatrixes)%block_vector(1)%block3dc , ex_r, comm, FFTW_MEASURE)
     ELSE
       plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint,ny_cint,nx_cint, &
                      ex_r, vold(nmatrixes)%block_vector(1)%block3dc,comm,FFTW_MPI_TRANSPOSED_OUT)
       plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(ny_cint,nz_cint,nx_cint, &
                      vnew(nmatrixes)%block_vector(1)%block3dc , ex_r,comm,FFTW_MPI_TRANSPOSED_IN)
     ENDIF
  ENDIF
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
    USE shared_data
    USE fields , ONLY : l_spectral
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
    DO i = 1,6
      DO j=1,11
        CALL is_calculation_needed(i,j,needed)
        IF(needed .EQV. .FALSE.) THEN 
          IF(ASSOCIATED(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc)) THEN
            DEALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    IF(l_spectral) THEN
      DO i=1,11
         DEALLOCATE(vold(nmatrixes)%block_vector(i)%block3dc)
         IF(i .LE.  6) DEALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc)
      ENDDO
    ENDIF
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


SUBROUTINE init_gpstd()
  USE matrix_coefficients
  USE PICSAR_PRECISION
  USE CONSTANTS
  USE mpi_fftw3
  USE omp_lib
  USE shared_data !, ONLY : dx,dy,dz,nx,ny,nz
  USE fields , ONLY : norderx,nordery,norderz,nxguards,nyguards,nzguards
  USE params , ONLY : dt

  INTEGER(idp)           :: i,j,k,p
  COMPLEX(cpx)           :: ii
  LOGICAL(lp)            :: needed
  INTEGER(idp)           :: nfftx,nffty,nfftz
  LOGICAL(lp)            :: switch
!  IF( fftw_with_mpi) THEN
!    IF(.NOT. fftw_mpi_transpose) THEN
!      nfftx=nx_global
!      nffty=ny_global
!      nfftz=local_nz
!    ELSE
!      nfftx = nx_global
!      nffty = local_ny 
!      nfftz = nz_global
!    ENDIF
!  ENDIF
!  IF(.NOT. fftw_with_mpi) THEN
!    nfftx = nx+2*nxguards
!    nffty = ny+2*nyguards
!    nfftz = nz+2*nzguards
!  ENDIF
  CALL select_case_dims_local(nfftx,nffty,nfftz)
        !print*,int(nfftx,isp),int(nffty,isp),int(nfftz,isp),"nfft"
  ii=DCMPLX(0.,1.)
  CALL allocate_new_matrix_vector(11_idp)
  CALL init_kspace
  DO i=1_idp,6_idp
    DO j=1_idp,11_idp
      CALL is_calculation_needed(i,j,needed) 
      IF(.NOT. needed) CYCLE
      ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc(nfftx/2+1,nffty,nfftz))
      cc_mat(nmatrixes)%block_matrix2d(i,j)%nx = nfftx/2+1
      cc_mat(nmatrixes)%block_matrix2d(i,j)%ny = nffty
      cc_mat(nmatrixes)%block_matrix2d(i,j)%nz = nfftz
    ENDDO
  ENDDO
  DO i=1_idp,11_idp
    ALLOCATE(vold(nmatrixes)%block_vector(i)%block3dc(nfftx/2+1,nffty,nfftz))
    !vold(nmatrixes)%block_vector(i)%block3dc = DCMPLX(0.,0.)
    vold(nmatrixes)%block_vector(i)%nx = nfftx/2+1
    vold(nmatrixes)%block_vector(i)%ny = nffty 
    vold(nmatrixes)%block_vector(i)%nz = nfftz
  ENDDO
  DO i=1_idp,6_idp
    ALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc(nfftx/2+1,nffty,nfftz))
    vnew(nmatrixes)%block_vector(i)%nx = nfftx/2+1
    vnew(nmatrixes)%block_vector(i)%ny = nffty
    vnew(nmatrixes)%block_vector(i)%nz = nfftz
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
  switch = .FALSE.
  IF(ABS(Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)) .EQ. 0.0) THEN
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1) = (1.0_num,0.0_num)
    switch = .TRUE.
  ENDIF
  DO i = 1,3
    cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc =        &
      DCMPLX(0.,1.)*(AT_OP(nmatrixes2)%block_vector(2)%block3dc   &
      -1./(clight*dt)*AT_OP(nmatrixes2)%block_vector(1)%block3dc)&
      /Kspace(nmatrixes2)%block_vector(10)%block3dc**2

    cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc =        &
       cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc &
       *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
    IF(switch) THEN
      cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc(1,1,1) = &
       -1.0_num/3.0_num*(0.0_num,1.0_num)*(clight*dt)**2
    ENDIF
    cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc = 1.0_num/eps0 &
      *cc_mat(nmatrixes)%block_matrix2d(i,10_idp)%block3dc
  ENDDO
  IF(switch) THEN
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1) = DCMPLX(0.,0.)
  ENDIF
  !contribution rho new (vold(11)) 
  DO i = 1 , 3
    cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc =  &
    DCMPLX(0.,1.)*(1./(clight*dt)*                          &
    AT_OP(nmatrixes2)%block_vector(1)%block3dc             &
    -DCMPLX(1.,0.))/Kspace(nmatrixes2)%block_vector(10)%block3dc**2

    cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc = &
       cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc &
       *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
    IF(switch) THEN
      cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc(1,1,1) = &
     -1.0_num/6.0_num*(0.0_num,1.0_num)*(clight*dt)**2  
    ENDIF
    cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc = 1.0_num/eps0 &
     *  cc_mat(nmatrixes)%block_matrix2d(i,11_idp)%block3dc     
  ENDDO
  IF(switch) THEN
    Kspace(nmatrixes2)%block_vector(10)%block3dc(1,1,1)   = DCMPLX(0.,0.)
  ENDIF
!if(rank .EQ. 0) THEN
!do i=1,local_nz
!print*,kzc(i),rank,int(local_nz,isp)
!enddo
!ENDIF
!call mpi_barrier(comm,errcode)
!if(rank .EQ. 1) THEN
!do i=1,local_nz
!print*,kzc(i),rank,int(local_nz,isp)
!enddo
!ENDIF


  CALL delete_arrays 
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
!  USE fourier_psaotd
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

SUBROUTINE execute_fftw_r2c_mpi
USE shared_data
USE fields
USE mpi_fftw3
USE time_stat
USE params
USE matrix_coefficients
USE fourier
USE fastfft

INTEGER(idp)  :: nfftx,nffty,nfftz
REAL(num)     :: tmptime
INTEGER(idp)  :: ix,iy,iz
nfftx = nx_global
nffty = ny_global
nfftz = local_nz


IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
! Copy array values before FFT
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=1,local_nz
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
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,ex_r,vold(nmatrixes)%block_vector(1)%block3dc)
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,ey_r,vold(nmatrixes)%block_vector(2)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,ez_r,vold(nmatrixes)%block_vector(3)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,bx_r,vold(nmatrixes)%block_vector(4)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,by_r,vold(nmatrixes)%block_vector(5)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,bz_r,vold(nmatrixes)%block_vector(6)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,jx_r,vold(nmatrixes)%block_vector(7)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,jy_r,vold(nmatrixes)%block_vector(8)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,jz_r,vold(nmatrixes)%block_vector(9)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,rhoold_r,vold(nmatrixes)%block_vector(10)%block3dc )
CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi,rho_r,vold(nmatrixes)%block_vector(11)%block3dc )
IF (it.ge.timestat_itstart) THEN
  localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
ENDIF
END SUBROUTINE


SUBROUTINE execute_fftw_mpi_c2r
USE params
USE shared_data
USE fields
USE fourier
USE fastfft
USE time_stat
USE mpi_fftw3
USE matrix_coefficients
IMPLICIT NONE
REAL(num) :: coeff_norm,tmptime
INTEGER(idp) :: ix,iy,iz,nfftx,nffty,nfftz

IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,vnew(nmatrixes)%block_vector(1)%block3dc,ex_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,vnew(nmatrixes)%block_vector(2)%block3dc,ey_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,vnew(nmatrixes)%block_vector(3)%block3dc,ez_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,vnew(nmatrixes)%block_vector(4)%block3dc,bx_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,vnew(nmatrixes)%block_vector(5)%block3dc,by_r)
CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,vnew(nmatrixes)%block_vector(6)%block3dc,bz_r)
IF (it.ge.timestat_itstart) THEN
  localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
ENDIF
coeff_norm=1.0_num/((nx_global)*(ny_global)*(nz_global))
IF (it.ge.timestat_itstart) THEN
  tmptime = MPI_WTIME()
ENDIF
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz) COLLAPSE(3)
DO iz=1,local_nz
  DO iy=1,ny_global
    DO ix=1,nx_global
      ex(ix-1,iy-1,iz-1)=ex_r(ix,iy,iz)*coeff_norm
      ey(ix-1,iy-1,iz-1)=ey_r(ix,iy,iz)*coeff_norm
      ez(ix-1,iy-1,iz-1)=ez_r(ix,iy,iz)*coeff_norm
      bx(ix-1,iy-1,iz-1)=bx_r(ix,iy,iz)*coeff_norm
      by(ix-1,iy-1,iz-1)=by_r(ix,iy,iz)*coeff_norm
      bz(ix-1,iy-1,iz-1)=bz_r(ix,iy,iz)*coeff_norm
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
IF (it.ge.timestat_itstart) THEN
  localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
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
  !  USE fourier_psaotd
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
                lognumer =LOG(16.0_num)*(1.0_num-p/2.0_num)+logfactorial(p-1_idp)*2.0_num
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
IF(n.EQ.0_idp) THEN
    logfactorial=0.
ELSE
     x=log(1.0_num*n)
     logfactorial=x
      DO k=2,n-1
          x=log(1.0_num*k)
          logfactorial=logfactorial+x
      ENDDO
ENDIF
RETURN
END FUNCTION logfactorial

SUBROUTINE normalize_Fourier(ex_out,n1,n2,n3,ex_in,nxx,nyy,nzz,coeff_norm)
USE PICSAR_precision
USE omp_lib
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
END SUBROUTINE normalize_Fourier
#endif
END MODULE














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
  COMPLEX(cpx), DIMENSION(:), ALLOCATABLE :: kxc, kxb, kxf, kyc, kyb, kyf, kzc, kzb,  &
  kzf
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
  SUBROUTINE select_case_dims_local(nfftx, nffty, nfftz)
    USE shared_data
    USE mpi_fftw3
    USE group_parameters
    USE params
    USE fields, ONLY : nxguards, nyguards, nzguards

    INTEGER(idp), INTENT(INOUT) :: nfftx, nffty, nfftz
    IF( fftw_with_mpi) THEN
      IF(.NOT. fftw_hybrid) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx=nx_global
          nffty=ny_global
          nfftz=local_nz
        ELSE
          nfftx=nx_global
          nffty=nz_global
          nfftz=local_ny
        ENDIF
      ELSE
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx = nx_group
          nffty = ny_group
          nfftz = local_nz
        ELSE
          nfftx = nx_group
          nffty = nz_group
          nfftz = local_ny
        ENDIF
      ENDIF
#if defined(P3DFFT)
      IF(p3dfft) THEN
        nfftx = p3d_fsize(1)
        nffty = p3d_fsize(2)
        nfftz = p3d_fsize(3)
      ENDIF
#endif

    ELSE
#if defined(LIBRARY)
      nfftx = nx+2*nxguards+1
      nffty = ny+2*nyguards+1
      nfftz = nz+2*nzguards+1
#else
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      nfftz = nz+2*nzguards
#endif
      IF(c_dim ==2) THEN
        nffty = 1
      ENDIF
    ENDIF
  END SUBROUTINE select_case_dims_local

  SUBROUTINE select_case_dims_global(nfftx, nffty, nfftz)
    USE shared_data
    USE mpi_fftw3
    USE group_parameters
    USE params
    USE fields, ONLY : nxguards, nyguards, nzguards

    INTEGER(idp), INTENT(INOUT) :: nfftx, nffty, nfftz
    IF( fftw_with_mpi) THEN
      IF(.NOT. fftw_hybrid) THEN
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx=nx_global
          nffty=ny_global
          nfftz=nz_global
        ELSE
          nfftx=nx_global
          nffty=nz_global
          nfftz=ny_global
        ENDIF
      ELSE
        IF(.NOT. fftw_mpi_transpose) THEN
          nfftx = nx_group
          nffty = ny_group
          nfftz = nz_group
        ELSE
          nfftx = nx_group
          nffty = nz_group
          nfftz = ny_group
        ENDIF
      ENDIF
    ELSE
#if defined(LIBRARY)
      nfftx = nx+2*nxguards+1
      nffty = ny+2*nyguards+1
      nfftz = nz+2*nzguards+1
#else
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      nfftz = nz+2*nzguards
#endif
    ENDIF
    IF(c_dim ==2) THEN
      nffty=1
    ENDIF

#if defined(P3DFFT) 
    IF(p3dfft) THEN
      IF(p3dfft_stride) THEN
       nfftx =nz_group
       nffty =ny_group
       nfftz =nx_group
      ELSE
       nfftx = nx_group 
       nffty = ny_group
       nfftz = nz_group 
     ENDIF
    ENDIF
#endif
  END SUBROUTINE select_case_dims_global


  SUBROUTINE init_kspace
    USE matrix_coefficients
    USE CONSTANTS
    USE mpi_fftw3
    USE omp_lib
    USE shared_data
    USE fields, ONLY : nxguards, nyguards, nzguards
    USE fields, ONLY : norderx, nordery, norderz, l_staggered
    USE params, ONLY : dt
    REAL(num), ALLOCATABLE, DIMENSION(:, :, :)    :: temp, temp2
    INTEGER(idp)                                  :: i, j, k
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: nfftx, nffty, nfftz, nfftxr
    LOGICAL(lp)                                   :: switch

    IF(.NOT. ASSOCIATED(Kspace)) THEN
      ALLOCATE(KSPACE(ns_max))
    ENDIF
    IF(.NOT. ASSOCIATED(AT_OP)) THEN
      ALLOCATE(AT_OP(ns_max))
    ENDIF
    nmatrixes2=nmatrixes2+1
    ALLOCATE(Kspace(nmatrixes2)%block_vector(10_idp))! 3 forward 3 backward
    IF(.NOT. ASSOCIATED(AT_OP)) THEN
      ALLOCATE(AT_OP(ns_max))
    ENDIF
    ALLOCATE(AT_OP(nmatrixes2)%block_vector(4_idp))!S/k, C, (1-C)/k^2
    CALL select_case_dims_local(nfftx, nffty, nfftz)
    nfftxr = nfftx/2+1
    IF(p3dfft) nfftxr = nfftx
    DO i = 1_idp, 10_idp
      ALLOCATE(Kspace(nmatrixes2)%block_vector(i)%block3dc(nfftxr, nffty, nfftz))
    ENDDO
    DO i = 1_idp, 4_idp
      ALLOCATE(AT_OP(nmatrixes2)%block_vector(i)%block3dc(nfftxr, nffty, nfftz))
    ENDDO
    !construct kspace
    ii=DCMPLX(0.0_num, 1.0_num)
    CALL compute_k_vec(l_staggered)
    DO i = 1, nfftxr
      DO j = 1, nffty
        DO k = 1, nfftz
          IF(.NOT. p3dfft) THEN
            IF(.NOT. fftw_mpi_transpose) THEN
              Kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
              Kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
              Kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
              IF(c_dim == 3) THEN
                Kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                Kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                Kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
              ELSE IF(c_dim == 2) THEN
                Kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k)= (0.0_num,0.0_num) 
                Kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k)= (0.0_num,0.0_num)
                Kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k)= (0.0_num,0.0_num) 
              ENDIF
              Kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kzf(k)
              Kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kzb(k)
              Kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kzc(k)
            ELSE
              Kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
              Kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
              Kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
              Kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kzf(k)
              Kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kzb(k)
              Kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kzc(k)
              Kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kyf(j)
              Kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kyb(j)
              Kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kyc(j)
            ENDIF
          ELSE
            IF(p3dfft_stride) THEN

                Kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kzf(k)
                Kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kzb(k)
                Kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kzc(k)
                Kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                Kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                Kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
                Kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kxf(i)
                Kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kxb(i)
                Kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kxc(i)

            ELSE 
                Kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
                Kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
                Kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
                Kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                Kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                Kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
                Kspace(nmatrixes2)%block_vector(7)%block3dc(i, j, k) = kzf(k)
                Kspace(nmatrixes2)%block_vector(8)%block3dc(i, j, k) = kzb(k)
                Kspace(nmatrixes2)%block_vector(9)%block3dc(i, j, k) = kzc(k)

            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    Kspace(nmatrixes2)%block_vector(10)%block3dc= SQRT(ABS(Kspace(nmatrixes2)%block_vector(9)%block3dc)**2 + &
        ABS(Kspace(nmatrixes2)%block_vector(6)%block3dc)**2 + &
        ABS(Kspace(nmatrixes2)%block_vector(3)%block3dc)**2)
    switch = .FALSE.
    ALLOCATE(temp(nfftxr, nffty, nfftz))
    ALLOCATE(temp2(nfftxr, nffty, nfftz))
    temp=dt*clight*REAL(Kspace(nmatrixes2)%block_vector(10)%block3dc, num)
    temp2=sinc_block(nfftxr, nffty, nfftz, temp)
    AT_OP(nmatrixes2)%block_vector(1)%block3dc = DCMPLX(temp2, 0._num)
    AT_OP(nmatrixes2)%block_vector(1)%block3dc =                                      &
    clight*dt*AT_OP(nmatrixes2)%block_vector(1)%block3dc
    temp2=COS(temp)
    AT_OP(nmatrixes2)%block_vector(2)%block3dc = DCMPLX(temp2, 0._num)
    temp=0.5_num*temp
    AT_OP(nmatrixes2)%block_vector(3)%block3dc = 2._num*(clight*dt/2.0_num)**2        &
    *sinc_block(nfftxr, nffty, nfftz, temp)*sinc_block(nfftxr, nffty, nfftz,    &
    temp)

    IF(ABS(Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)) .EQ. 0.0_num) THEN
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = DCMPLX(1.0_num,         &
      0.0_num)
      switch = .TRUE.

    ENDIF
    AT_OP(nmatrixes2)%block_vector(3)%block3dc = (DCMPLX(1.0_num, 0.0_num) -          &
    AT_OP(nmatrixes2)%block_vector(2)%block3dc)                                       &
    /Kspace(nmatrixes2)%block_vector(10)%block3dc**2
    !(1-C)/k^2
    AT_OP(nmatrixes2)%block_vector(4)%block3dc =                                      &
    (AT_OP(nmatrixes2)%block_vector(1)%block3dc-clight*dt) /                          &
    Kspace(nmatrixes2)%block_vector(10)%block3dc/                                     &
    Kspace(nmatrixes2)%block_vector(10)%block3dc
    IF(switch) THEN
      AT_OP(nmatrixes2)%block_vector(3)%block3dc(1, 1, 1) = (clight*dt)**2/2.0_num
      AT_OP(nmatrixes2)%block_vector(4)%block3dc(1, 1,                                &
      1)=DCMPLX(-(clight*dt)**3/6.0_num, 0.0_num)
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)=DCMPLX(0._num, 0._num)
    ENDIF
    DEALLOCATE(temp, temp2)
  END SUBROUTINE init_kspace

  SUBROUTINE compute_k_vec(l_stg)
    USE constants
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE group_parameters
    IMPLICIT NONE
    LOGICAL(lp), INTENT(IN)                     :: l_stg
    COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     :: kxct,kxbt,kxft,kyct,kybt,kyft,kzct,kzbt,kzft, k_temp
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: i, j, k
    INTEGER(idp)                                  :: nfftx, nffty, nfftz
    REAL(num)                                     :: sd
    INTEGER(idp)                                  :: temp_order

#if defined(LIBRARY)
    !Due to different staggering in PICSAR and SMILEI
    ii = DCMPLX(0.0_num, -1.0_num)
#else
    ii = DCMPLX(0.0_num, 1.0_num)
#endif
    IF(fftw_mpi_transpose) THEN
     ! switch  z and y
      sd=dz
      dz=dy
      dy=sd
      temp_order = norderz
      norderz=nordery
      nordery=temp_order
    ENDIF
    IF(p3dfft) THEN
      ! switch z and x
      IF(p3dfft_stride) THEN
        sd = dz
        dz = dx
        dx = sd
        temp_order = norderz
        norderz = norderx
        norderx = temp_order
      ENDIF
    ENDIF       
    CALL select_case_dims_global(nfftx, nffty, nfftz)
    CALL compute_k_1d( nfftx,kxc,kxf,kxb,norderx,dx,l_stg)
    CALL compute_k_1d( nffty,kyc,kyf,kyb,nordery,dy,l_stg)
    CALL compute_k_1d( nfftz,kzc,kzf,kzb,norderz,dz,l_stg)
    ! delete second part of kx because real
    IF(.NOT. p3dfft) THEN
      ALLOCATE(k_temp(nfftx));
      k_temp = kxc;
      DEALLOCATE(kxc); ALLOCATE(kxc(nfftx/2+1)) ; kxc = k_temp(1:nfftx/2+1)
      k_temp = kxb;
      DEALLOCATE(kxb); ALLOCATE(kxb(nfftx/2+1)) ; kxb = k_temp(1:nfftx/2+1)
      k_temp = kxf;
      DEALLOCATE(kxf); ALLOCATE(kxf(nfftx/2+1)) ; kxf = k_temp(1:nfftx/2+1)
      DEALLOCATE(k_temp)
    ENDIF

    IF(fftw_with_mpi) THEN
      IF( .NOT. p3dfft) THEN
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
          ALLOCATE(k_temp(nfftz))
          k_temp = kzc
          DEALLOCATE(kzc);ALLOCATE(kzc(local_ny))
          kzc = k_temp(local_y0+1:local_y0+local_ny)
          k_temp = kzf
          DEALLOCATE(kzf);ALLOCATE(kzf(local_ny))
          kzf = k_temp(local_y0+1:local_y0+local_ny)
          k_temp = kzb
          DEALLOCATE(kzb);ALLOCATE(kzb(local_ny))
          kzb = k_temp(local_y0+1:local_y0+local_ny)
        ENDIF
        DEALLOCATE(k_temp)
      ELSE IF(p3dfft) THEN
        ALLOCATE(kxct(nfftx),kxbt(nfftx),kxft(nfftx),kyct(nffty),kybt(nffty),kyft(nffty),&
          kzct(nfftz),kzbt(nfftz),kzft(nfftz))
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
    IF(fftw_mpi_transpose) THEN
      sd=dz
      dz=dy
      dy=sd
      temp_order = norderz
      norderz=nordery
      nordery=temp_order
    ENDIF
    IF(p3dfft) THEN
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
  SUBROUTINE compute_k_1d(nfft,kvec,kvecf,kvecb,norder,d,l_stg)
     USE picsar_precision
     USE constants
     REAL(num) , INTENT(IN)  :: d
     INTEGER(idp) , INTENT(IN) :: norder,nfft
     COMPLEX(cpx) , DIMENSION(:) , ALLOCATABLE , INTENT(INOUT) :: kvec,kvecf,kvecb
     LOGICAL(lp)  , INTENT(IN)             :: l_stg
     COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     ::  ones, onesp
     REAL(num), ALLOCATABLE, DIMENSION(:)        :: FD
     INTEGER(idp)                                ::j,i
     COMPLEX(cpx)                                ::  ii
    
     ii = (0.0_num,1.0_num)

     ALLOCATE(ones(nfft), onesp(nfft))
     ALLOCATE(kvec(nfft),kvecf(nfft),kvecb(nfft))
     DO j=1_idp, nfft
       ones(j)  = DCMPLX(j-1.0_num, 0.0_num)
       onesp(j) = DCMPLX(j-1.0_num, 0.0_num)
       IF(j .GT. nfft/2_idp +1) THEN
         ones(j)  =DCMPLX(-ones(j))
         onesp(j) =DCMPLX( nfft + onesp(j))
       ENDIF
     ENDDO
     IF (norder .ne. 0_idp) THEN
       ALLOCATE(FD(norder/2))
       CALL FD_weights_hvincenti(norder, FD, l_stg)
       kvec=(0._num, 0._num)*kvec
       kvecb=(0._num, 0._num)*kvecb
       kvecf=(0._num, 0._num)*kvecf
       DO i=1_idp, norder/2
         kvec=kvec+2.0_num/d*FD(i)*SIN((i*2.0_num-1.0_num)*PI*ones/nfft)
       ENDDO
     ELSE
       CALL fftfreq(nfft, kvec,  d)
     ENDIF
     IF(l_stg) THEN
       kvecf=kvec*EXP(-ii*PI*onesp/nfft)
       kvecb=kvec*EXP(ii*PI*onesp/nfft)
     ELSE
       kvecb=kvec
       kvecf=kvec
     ENDIF
     DEALLOCATE(onesp,ones,FD)
  END SUBROUTINE compute_k_1d
  SUBROUTINE fftfreq(nxx, kxx, dxx)
    USE constants
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)                    :: nxx
    REAL(num), INTENT(IN)                    :: dxx
    COMPLEX(cpx), INTENT(OUT), DIMENSION(:)  :: kxx
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
    kxx=kxx/(dxx*nxx)/2.0_num*PI
  END SUBROUTINE fftfreq

  FUNCTION sinc_block(n1, n2, n3, block)
    USE  constants
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
  !> This subroutine constructs gpstd_blocks and puts them into cc_mat operator
  !
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ________________________________________________________________________________________


  SUBROUTINE init_gpstd() bind(C, name='init_gpstd_pxr')
    USE matrix_coefficients
    USE PICSAR_PRECISION
    USE CONSTANTS
    USE mpi_fftw3
    USE omp_lib
    USE shared_data!, ONLY : dx, dy, dz, nx, ny, nz
    USE fields, ONLY : g_spectral, norderx, nordery, norderz, nxguards, nyguards,     &
         nzguards, exf, eyf, ezf, bxf, byf, bzf, jxf, jyf, jzf, rhooldf, rhof
    USE params, ONLY : dt

    INTEGER(idp)           :: i, j, k, p
    COMPLEX(cpx)           :: ii
    LOGICAL(lp)            :: needed
    INTEGER(idp)           :: nfftx, nffty, nfftz,nfftxr
    LOGICAL(lp)            :: switch
    REAL(num)              :: coeff_norm

    CALL select_case_dims_local(nfftx, nffty, nfftz)
    ii=DCMPLX(0.0_num, 1.0_num)
    CALL allocate_new_matrix_vector(11_idp)
    nfftxr = nfftx/2+1
    IF(p3dfft) nfftxr = nfftx 
    CALL init_kspace
    nkx = nfftxr 
    nky = nffty
    nkz = nfftz
    DO i=1_idp, 11_idp
      DO j=1_idp, 11_idp
          CALL is_calculation_needed(i, j, needed)
          IF(g_spectral .OR. needed) THEN
            ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i, j)%block3dc(nfftxr, nffty,    &
            nfftz))
            cc_mat(nmatrixes)%block_matrix2d(i, j)%nx = nfftxr
            cc_mat(nmatrixes)%block_matrix2d(i, j)%ny = nffty
            cc_mat(nmatrixes)%block_matrix2d(i, j)%nz = nfftz
          ENDIF
      ENDDO
    ENDDO
    IF(g_spectral) THEN
        vold(nmatrixes)%block_vector(1)%block3dc => exf
        vold(nmatrixes)%block_vector(2)%block3dc => eyf
        vold(nmatrixes)%block_vector(3)%block3dc => ezf       
        vold(nmatrixes)%block_vector(4)%block3dc => bxf
        vold(nmatrixes)%block_vector(5)%block3dc => byf
        vold(nmatrixes)%block_vector(6)%block3dc => bzf
        vold(nmatrixes)%block_vector(7)%block3dc => jxf
        vold(nmatrixes)%block_vector(8)%block3dc => jyf
        vold(nmatrixes)%block_vector(9)%block3dc => jzf
        vold(nmatrixes)%block_vector(10)%block3dc => rhooldf
        vold(nmatrixes)%block_vector(11)%block3dc => rhof
      DO i=1_idp,11_idp
        ALLOCATE(vnew(nmatrixes)%block_vector(i)%block3dc(nfftxr, nffty,&
        nfftz))
        vnew(nmatrixes)%block_vector(i)%nx = nfftxr
        vnew(nmatrixes)%block_vector(i)%ny = nffty
        vnew(nmatrixes)%block_vector(i)%nz = nfftz
        vold(nmatrixes)%block_vector(i)%nx = nfftxr
        vold(nmatrixes)%block_vector(i)%ny = nffty
        vold(nmatrixes)%block_vector(i)%nz = nfftz
      ENDDO
      DO i = 1,11
        DO j=1,11 
          cc_mat(nmatrixes)%block_matrix2d(i, j)%block3dc = CMPLX(0.0_num,0.0_num)
        ENDDO
          cc_mat(nmatrixes)%block_matrix2d(i, i)%block3dc = CMPLX(1.0_num,0.0_num) 
      ENDDO
    ENDIF


    cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc = -                               &
    ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc =                                 &
    ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc
    cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc =                                 &
    ii*Kspace(nmatrixes2)%block_vector(7)%block3dc*clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(2, 6)%block3dc =                                 &
    -ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight                            &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc = -                               &
    ii*Kspace(nmatrixes2)%block_vector(4)%block3dc*clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(3, 5)%block3dc =                                 &
    ii*Kspace(nmatrixes2)%block_vector(1)%block3dc*clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc =                                 &
    ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc =                                 &
    -ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight                            &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc =                                 &
    -ii*Kspace(nmatrixes2)%block_vector(8)%block3dc/clight                            &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(5, 3)%block3dc =                                 &
    ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc =                                 &
    ii*Kspace(nmatrixes2)%block_vector(5)%block3dc/clight                             &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    cc_mat(nmatrixes)%block_matrix2d(6, 2)%block3dc =                                 &
    -ii*Kspace(nmatrixes2)%block_vector(2)%block3dc/clight                            &
    *AT_OP(nmatrixes2)%block_vector(1)%block3dc

    DO i=1, 6
      cc_mat(nmatrixes)%block_matrix2d(i, i)%block3dc =                               &
      AT_OP(nmatrixes2)%block_vector(2)%block3dc
    ENDDO
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, i+6)%block3dc =                             &
      (-1._num)*clight*mu0*AT_OP(nmatrixes2)%block_vector(1)%block3dc
    ENDDO
    cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc = - mu0*                          &
    ii*Kspace(nmatrixes2)%block_vector(8)%block3dc*                                   &
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc = -                               &
    mu0*(-ii)*Kspace(nmatrixes2)%block_vector(5)%block3dc*                            &
    AT_OP(nmatrixes2)%block_vector(3)%block3dc



    cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc = -                               &
    mu0*(-ii)*Kspace(nmatrixes2)%block_vector(8)%block3dc*                            &
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc = - mu0*                          &
    ii*Kspace(nmatrixes2)%block_vector(2)%block3dc*                                   &
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc = - mu0*                          &
    ii*Kspace(nmatrixes2)%block_vector(5)%block3dc*                                   &
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc =                                 &
    -mu0*(-ii)*Kspace(nmatrixes2)%block_vector(2)%block3dc*                           &
    AT_OP(nmatrixes2)%block_vector(3)%block3dc


    !contribution rho old
    switch = .FALSE.
    IF(ABS(Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)) .EQ. 0.0_num) THEN
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = (1.0_num, 0.0_num)
      switch = .TRUE.
    ENDIF
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc = DCMPLX(0.,               &
      1.)*(AT_OP(nmatrixes2)%block_vector(2)%block3dc                                 &
      -1./(clight*dt)*AT_OP(nmatrixes2)%block_vector(1)%block3dc)                     &
      /Kspace(nmatrixes2)%block_vector(10)%block3dc**2

      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc =                          &
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc                            &
      *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
      IF(switch) THEN
        cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc(1, 1, 1) =               &
        -1.0_num/3.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
      ENDIF
      cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc = 1.0_num/eps0             &
      *cc_mat(nmatrixes)%block_matrix2d(i, 10_idp)%block3dc
    ENDDO
    IF(switch) THEN
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = DCMPLX(0.0_num,         &
      0.0_num)
    ENDIF
    !contribution rho new 
    DO i = 1, 3
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc = DCMPLX(0.,               &
      1.)*(1./(clight*dt)* AT_OP(nmatrixes2)%block_vector(1)%block3dc -DCMPLX(1.,     &
      0.))/Kspace(nmatrixes2)%block_vector(10)%block3dc**2

      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc =                          &
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc                            &
      *Kspace(nmatrixes2)%block_vector(3*i-1)%block3dc
      IF(switch) THEN
        cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc(1, 1, 1) =               &
        -1.0_num/6.0_num*(0.0_num, 1.0_num)*(clight*dt)**2
      ENDIF
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc = 1.0_num/eps0 *           &
      cc_mat(nmatrixes)%block_matrix2d(i, 11_idp)%block3dc
    ENDDO
    IF(switch) THEN
      Kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)   = DCMPLX(0., 0.)
    ENDIF
    CALL select_case_dims_global(nfftx,nffty,nfftz)
    coeff_norm = 1.0_num/(nfftx*nffty*nfftz)  
    DO i=1,11
      DO j=1,11
        CALL is_calculation_needed(i, j, needed)
        IF(needed .OR. g_spectral) THEN
          cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc = coeff_norm*cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE init_gpstd

  SUBROUTINE FD_weights_hvincenti(p, w, is_staggered)
    USE picsar_PRECISION
    IMPLICIT NONE
    LOGICAL(idp), INTENT(IN) :: is_staggered
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

  ! - Computes factorial of n
  FUNCTION factorial(n)
    USE constants
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
    use PICSAR_PRECISION
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


  SUBROUTINE copy_field(ex_out, n1, n2, n3, ex_in, nxx, nyy, nzz)
    USE PICSAR_precision
    USE omp_lib
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

  SUBROUTINE copy_field_forward()
    USE PICSAR_precision
    USE omp_lib
    USE fields
    USE shared_data
    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=-nzguards, nz+nzguards-1
      DO iy=-nyguards,ny+nyguards-1
        DO ix=-nxguards,nx+nxguards-1
          ex_r(ix, iy, iz)=ex(ix, iy, iz)
          ey_r(ix, iy, iz)=ey(ix, iy, iz)
          ez_r(ix, iy, iz)=ez(ix, iy, iz)
          bx_r(ix, iy, iz)=bx(ix, iy, iz)
          by_r(ix, iy, iz)=by(ix, iy, iz)
          bz_r(ix, iy, iz)=bz(ix, iy, iz)
          jx_r(ix, iy, iz)=jx(ix, iy, iz)
          jy_r(ix, iy, iz)=jy(ix, iy, iz)
          jz_r(ix, iy, iz)=jz(ix, iy, iz)
          rho_r(ix, iy, iz)=rho(ix, iy, iz)
          rhoold_r(ix, iy, iz)=rhoold(ix, iy, iz)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE copy_field_forward

  SUBROUTINE copy_field_backward()
    USE PICSAR_precision
    USE omp_lib
    USE fields
    USE shared_data
    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=-nzguards, nz+nzguards-1
      DO iy=-nyguards,ny+nyguards-1
        DO ix=-nxguards,nx+nxguards-1
          ex(ix, iy, iz)=ex_r(ix, iy, iz)
          ey(ix, iy, iz)=ey_r(ix, iy, iz)
          ez(ix, iy, iz)=ez_r(ix, iy, iz)
          bx(ix, iy, iz)=bx_r(ix, iy, iz)
          by(ix, iy, iz)=by_r(ix, iy, iz)
          bz(ix, iy, iz)=bz_r(ix, iy, iz)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE copy_field_backward


  SUBROUTINE is_calculation_needed(irow, icol, needed)
    USE picsar_precision
    INTEGER(idp), INTENT(IN)    :: irow, icol
    LOGICAL(lp), INTENT(INOUT) :: needed

    needed = .TRUE.
    IF(irow .LE. 3_idp) THEN
      IF((icol .LE. 3_idp) .AND. (icol .NE. irow)) THEN
        needed = .FALSE.
        RETURN
      ENDIF
      IF(icol  .EQ. irow+3_idp) THEN
        needed = .FALSE.
        RETURN
      ENDIF
      IF((icol .GE. 7_idp) .AND. (icol .LE. 9_idp) .AND. (icol .NE. irow + 6_idp))    &
      THEN
      needed = .FALSE.
      RETURN
    ENDIF
  ENDIF
  IF ((irow .LE. 6_idp) .AND. (irow .GE. 4_idp)) THEN
    IF ((icol .LE. 3_idp) .AND. icol .EQ. irow - 3_idp) THEN
      needed = .FALSE.
      RETURN
    ENDIF
    IF((icol .LE. 6_idp) .AND. (icol .GE. 4_idp) .AND. icol .NE. irow) THEN
      needed = .FALSE.
      RETURN
    ENDIF
    IF((icol .GE. 7_idp) .AND. (icol .EQ. irow + 3_idp)) THEN
      needed = .FALSE.
      RETURN
    ENDIF
    IF((icol .EQ. 10_idp) .OR. (icol .EQ. 11_idp)) THEN
      needed = .FALSE.
      RETURN
    ENDIF
  ENDIF
  IF (irow .GE. 7_idp) THEN
    needed = .FALSE.
    RETURN
  ENDIF
END SUBROUTINE is_calculation_needed
#endif
END MODULE














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
END MODULE math_tools


MODULE gpstd_solver
  USE math_tools
  USE PICSAR_PRECISION
  IMPLICIT NONE
  COMPLEX(cpx), DIMENSION(:), ALLOCATABLE :: kxc, kxb, kxf, kyc, kyb, kyf, kzc, kzb,  &
  kzf

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
    USE fields, ONLY: nxguards, nyguards, nzguards
    USE group_parameters, ONLY: nx_group, ny_group, nz_group, p3d_fsize
    USE iso_c_binding
    USE mpi_fftw3, ONLY: local_nz, local_nz_tr
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: c_dim, fftw_hybrid, fftw_mpi_transpose, fftw_with_mpi,    &
      nx, nx_global, ny, ny_global, nz, nz_global, p3dfft_flag
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
      nfftz = nz+2*nzguards+1
#else
      !> When using picsar
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      nfftz = nz+2*nzguards
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
    USE fields, ONLY: nxguards, nyguards, nzguards
    USE group_parameters, ONLY: nx_group, ny_group, nz_group
    USE iso_c_binding
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: c_dim, fftw_hybrid, fftw_mpi_transpose, fftw_with_mpi,    &
      nx, nx_global, ny, ny_global, nz, nz_global, p3dfft_flag, p3dfft_stride
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
      nfftz = nz+2*nzguards+1
#else
      !> When using only picsar
      nfftx = nx+2*nxguards
      nffty = ny+2*nyguards
      nfftz = nz+2*nzguards
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
    USE params, ONLY: dt
    USE picsar_precision, ONLY: cpx, idp, lp, num
    USE shared_data, ONLY: c_dim, fftw_mpi_transpose, p3dfft_flag, p3dfft_stride

    REAL(num), ALLOCATABLE, DIMENSION(:, :, :)    :: temp, temp2
    INTEGER(idp)                                  :: i, j, k
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: nfftx, nffty, nfftz, nfftxr
    LOGICAL(lp)                                   :: switch

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
    nfftxr = nfftx/2+1
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
    DO k = 1, nfftz
      DO j = 1, nffty
        DO i = 1, nfftxr
          IF(.NOT. p3dfft_flag) THEN
            IF(.NOT. fftw_mpi_transpose) THEN
              kspace(nmatrixes2)%block_vector(1)%block3dc(i, j, k) = kxf(i)
              kspace(nmatrixes2)%block_vector(2)%block3dc(i, j, k) = kxb(i)
              kspace(nmatrixes2)%block_vector(3)%block3dc(i, j, k) = kxc(i)
              IF(c_dim == 3) THEN
                kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k) = kyf(j)
                kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k) = kyb(j)
                kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k) = kyc(j)
              ELSE IF(c_dim == 2) THEN
                !> If c_dim == 2 Then y derivative is null
                !> c_dim = 2 cannot be used with p3dfft or fftw_mpi_transpose
                !> flag
                kspace(nmatrixes2)%block_vector(4)%block3dc(i, j, k)= (0.0_num,0.0_num) 
                kspace(nmatrixes2)%block_vector(5)%block3dc(i, j, k)= (0.0_num,0.0_num)
                kspace(nmatrixes2)%block_vector(6)%block3dc(i, j, k)= (0.0_num,0.0_num) 
              ENDIF
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
        ENDDO
      ENDDO
    ENDDO

    !> Computes the norm of wave vector in fourier space 
    kspace(nmatrixes2)%block_vector(10)%block3dc=                                    &
    SQRT(ABS(kspace(nmatrixes2)%block_vector(9)%block3dc)**2 +                       &
        ABS(kspace(nmatrixes2)%block_vector(6)%block3dc)**2 +                        &
        ABS(kspace(nmatrixes2)%block_vector(3)%block3dc)**2)
    switch = .FALSE.

    ALLOCATE(temp(nfftxr, nffty, nfftz))
    ALLOCATE(temp2(nfftxr, nffty, nfftz))

    temp=dt*clight*REAL(kspace(nmatrixes2)%block_vector(10)%block3dc, num)
    temp2=sinc_block(nfftxr, nffty, nfftz, temp)

    at_op(nmatrixes2)%block_vector(1)%block3dc = DCMPLX(temp2, 0._num)

    at_op(nmatrixes2)%block_vector(1)%block3dc =                                      &
    clight*dt*at_op(nmatrixes2)%block_vector(1)%block3dc
    temp2=COS(temp)

    at_op(nmatrixes2)%block_vector(2)%block3dc = DCMPLX(temp2, 0._num)
    temp=0.5_num*temp

    at_op(nmatrixes2)%block_vector(3)%block3dc = 2._num*(clight*dt/2.0_num)**2        &
    *sinc_block(nfftxr, nffty, nfftz, temp)*sinc_block(nfftxr, nffty, nfftz,    &
    temp)

    !> if current mpi task contains the null frequency then this processor it
    !> tagged by switch = .TRUE. in order perform Taylor expansion
    !> for at_op...(i)(1,1,1) only in this mpi task
    IF(ABS(kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)) .EQ. 0.0_num) THEN
      kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1) = DCMPLX(1.0_num,         &
      0.0_num)
      switch = .TRUE.
    ENDIF
    
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

    IF(switch) THEN
      at_op(nmatrixes2)%block_vector(3)%block3dc(1, 1, 1) = (clight*dt)**2/2.0_num
      at_op(nmatrixes2)%block_vector(4)%block3dc(1, 1,                                &
      1)=DCMPLX(-(clight*dt)**3/6.0_num, 0.0_num)
      kspace(nmatrixes2)%block_vector(10)%block3dc(1, 1, 1)=DCMPLX(0._num, 0._num)
    ENDIF
    DEALLOCATE(temp, temp2)
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
    !DEALLOCATE(kxc,kxb,kxf,kyc,kyb,kyf,kzc,kzb,kzf)

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
    USE fields, ONLY: norderx, nordery, norderz
    USE group_parameters, ONLY: p3d_fend, p3d_fsize, p3d_fstart
    USE iso_c_binding
    USE mpi_fftw3, ONLY: local_nz, local_nz_tr, local_z0, local_z0_tr
    USE picsar_precision, ONLY: cpx, idp, lp, num
    USE shared_data, ONLY: dx, dy, dz, fftw_mpi_transpose, fftw_with_mpi, nz,        &
      p3dfft_flag, p3dfft_stride
    IMPLICIT NONE
    LOGICAL(lp), INTENT(IN)                     :: l_stg
    COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     ::                                     &
     kxct,kxbt,kxft,kyct,kybt,kyft,kzct,kzbt,kzft, k_temp
    COMPLEX(cpx)                                  :: ii
    INTEGER(idp)                                  :: nfftx, nffty, nfftz
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
    CALL compute_k_1d( nfftx,kxc,kxf,kxb,norderx,dx,l_stg)
    CALL compute_k_1d( nffty,kyc,kyf,kyb,nordery,dy,l_stg)
    CALL compute_k_1d( nfftz,kzc,kzf,kzb,norderz,dz,l_stg)

    ! Selects only haf of  kx because r2c and c2r ffts
    IF(.NOT. p3dfft_flag) THEN
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
     LOGICAL(lp)  , INTENT(IN)             :: l_stg
     COMPLEX(cpx), ALLOCATABLE, DIMENSION(:)     ::  ones, onesp
     REAL(num), ALLOCATABLE, DIMENSION(:)        :: FD
     INTEGER(idp)                                ::j,i
     COMPLEX(cpx)                                ::  ii

     ii = (0.0_num,1.0_num)

     ALLOCATE(ones(nfft), onesp(nfft))
     ALLOCATE(kvec(nfft),kvecf(nfft),kvecb(nfft))
     kvec=(0._num, 0._num)
     kvecb=(0._num, 0._num)
     kvecf=(0._num, 0._num)
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
              kvec=kvec+2.0_num/d*FD(i)*SIN(i*2.0_num*PI*onesp/nfft)
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

     DEALLOCATE(onesp,ones)
  END SUBROUTINE compute_k_1d

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
    kxx=kxx/(dxx*nxx)/2.0_num*PI
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
      rhooldf
    USE iso_c_binding
    USE matrix_coefficients, ONLY: at_op, cc_mat, kspace, vnew, vold
    USE matrix_data, ONLY: nmatrixes, nmatrixes2
    USE mpi_fftw3, ONLY: alloc_local, fftw_alloc_complex
    USE omp_lib
    USE params, ONLY: dt
    USE picsar_precision, ONLY: cpx, idp, lp, num
    USE shared_data, ONLY: absorbing_bcs, c_dim, fftw_with_mpi, nkx, nky, nkz, nx,   &
      ny, nz, p3dfft_flag

    INTEGER(idp)           :: i, j,incr, lin_ind
    COMPLEX(cpx)           :: ii
    INTEGER(idp)           :: nfftx, nffty, nfftz,nfftxr, nbloc_ccmat, nbloc_vnew
    LOGICAL(lp)            :: switch
    REAL(num)              :: coeff_norm
    TYPE(C_PTR)            :: cdata
    INTEGER(idp) , ALLOCATABLE, DIMENSION(:) :: is_usefull
 
    IF(absorbing_bcs) THEN
      !> When using pmls, cc_mat is a 12x17 matrix
      nbloc_ccmat = 17_idp
      nbloc_vnew = 12_idp
      ALLOCATE(is_usefull(51))
      is_usefull=is_usefull_abs
    ELSE IF(.NOT. absorbing_bcs) THEN
      !> When using peridic bcs, cc_mat is a 6x11 matrix
      nbloc_ccmat = 11_idp
      nbloc_vnew = 6_idp
      ALLOCATE(is_usefull(33))
      is_usefull=is_usefull_per
    ENDIF
  
    CALL select_case_dims_local(nfftx, nffty, nfftz)
    ii=DCMPLX(0.0_num, 1.0_num)
    CALL allocate_new_matrix_vector(nbloc_ccmat)
    nfftxr = nfftx/2+1
    IF(p3dfft_flag) nfftxr = nfftx
    CALL init_kspace
    nkx = nfftxr
    nky = nffty
    nkz = nfftz
    !> Allocates cc_mat  block matrix
    !> cc_mat blocks are initally as an nbloc_ccmat x nbloc_ccmat block matrix 
    !> At the end of the routine, useless blcoks are deleted  
    incr = 1
    DO i=1_idp, nbloc_ccmat
      DO j=1_idp, nbloc_ccmat
        lin_ind = (i-1)*nbloc_ccmat + (j-1) 
        IF (is_usefull(incr) == lin_ind) THEN
          ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i, j)%block3dc(nfftxr, nffty,    &
          nfftz))
          cc_mat(nmatrixes)%block_matrix2d(i, j)%nx = nfftxr
          cc_mat(nmatrixes)%block_matrix2d(i, j)%ny = nffty
          cc_mat(nmatrixes)%block_matrix2d(i, j)%nz = nfftz
          IF (incr<SIZE(is_usefull)) THEN
             incr = incr + 1
          ENDIF
        ELSE  
          ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(i, j)%block3dc(1,1,1))
          cc_mat(nmatrixes)%block_matrix2d(i, j)%nx = 1_idp
          cc_mat(nmatrixes)%block_matrix2d(i, j)%ny = 1_idp
          cc_mat(nmatrixes)%block_matrix2d(i, j)%nz = 1_idp
          
        ENDIF
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
      CALL compute_cc_mat_merged_fields()
    ENDIF
  
 
    !> Renormalize cc_mat blocks
    !> Because fftw_r2c and followed by fftw_c2r multiplies fields by 
    !> nfftx*nffty*nfftz 
    !> This way, no need to normalize fields in a separate step


    ! Introduce fft normalisation factor in mat bloc mult
    CALL select_case_dims_global(nfftx,nffty,nfftz)
    coeff_norm = 1.0_num/(nfftx*nffty*nfftz)

    DO i=1,nbloc_ccmat
      DO j=1,nbloc_ccmat
          cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc =                            &
          coeff_norm*cc_mat(nmatrixes)%block_matrix2d(i,j)%block3dc
      ENDDO
    ENDDO



    !> Delete kspace and at_op blocks
    !> Might not delete these blocks if current filtering or field correction is
    !> needed in Fourier space
    CALL delete_k_space
    DEALLOCATE(is_usefull)
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

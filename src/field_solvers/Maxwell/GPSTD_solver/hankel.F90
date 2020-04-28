
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


SUBROUTINE get_Hfields()
  USE PICSAR_precision
  USE gpstd_solver
  USE fourier_psaotd
  USE shared_data, ONLY: nmodes
  USE matrix_coefficients, ONLY : hankel_mat
  USE matrix_data, ONLY:  nmatrixes_h
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
  !COMPLEX(cpx),POINTER, DIMENSION(:, :) :: el_h_ptr, ep_h_ptr, em_h_ptr, bl_h_ptr, bp_h_ptr, &
  !                          bm_h_ptr, jl_h_ptr, jp_h_ptr, jm_h_ptr, rhoold_h_ptr, rho_h_ptr, &
  !                               el_f_ptr, ep_f_ptr, em_f_ptr, bl_f_ptr, bp_f_ptr, bm_f_ptr, &
  !                                  jl_f_ptr, jp_f_ptr, jm_f_ptr, rho_f_ptr, rhoold_f_ptr 
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
  !#if defined(DEBUG) 
  !WRITE(0, *) "before fourier transform"
  !#endif 
  Call get_Ffields_AM_rz()
  !#if defined(DEBUG) 
  !WRITE(0, *) "after fourier transform"
  !#endif 
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
      !#if defined(DEBUG)
      !WRITE(0, *) "Hankel matrix not through yet"
      !#endif 
    !Call Hankel_M_and_invM(imode-1)
      !#if defined(DEBUG)
      !WRITE(0, *) "Hankel matrix through"
      ! #endif 

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

    Call  dgemm_example(el_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 1)%block2dc, el_h(:,:,imode), nfftr, nffty )
    Call  dgemm_example(em_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 2)%block2dc, em_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(ep_f(:,:,imode), hankel_mat(nmatrixes_h)% & 
      mode_block_matrix2d(imode, 3)%block2dc, ep_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bl_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 1)%block2dc, bl_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bm_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 2)%block2dc, bm_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bp_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 3)%block2dc, bp_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(jl_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 1)%block2dc, jl_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(jm_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 2)%block2dc, jm_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(jp_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 3)%block2dc, jp_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(rho_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 1)%block2dc, rho_h(:,:,imode), nfftr, nffty)
    Call  dgemm_example(rhoold_f(:,:,imode), hankel_mat(nmatrixes_h)% &
      mode_block_matrix2d(imode, 1)%block2dc, rhoold_h(:,:,imode), nfftr, nffty)
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
  USE gpstd_solver
  USE shared_data, ONLY: nmodes
  USE matrix_coefficients, ONLY : hankel_mat
  USE matrix_data, ONLY:  nmatrixes_h
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
  !  Call Hankel_M_and_invM(imode-1)
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

    Call  dgemm_example(el_h(:,:,imode), hankel_mat(nmatrixes_h)% & 
       mode_block_matrix2d(imode, 4)%block2dc, el_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(em_h(:,:,imode), hankel_mat(nmatrixes_h)% &
       mode_block_matrix2d(imode, 5)%block2dc, em_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(ep_h(:,:,imode), hankel_mat(nmatrixes_h)% &
       mode_block_matrix2d(imode, 6)%block2dc, ep_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bl_h(:,:,imode), hankel_mat(nmatrixes_h)% &
       mode_block_matrix2d(imode, 4)%block2dc, bl_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bm_h(:,:,imode), hankel_mat(nmatrixes_h)% &
       mode_block_matrix2d(imode, 5)%block2dc, bm_h_inv(:,:,imode), nfftr, nffty)
    Call  dgemm_example(bp_h(:,:,imode), hankel_mat(nmatrixes_h)% &
       mode_block_matrix2d(imode, 6)%block2dc, bp_h_inv(:,:,imode), nfftr, nffty)
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


!write (*,*) "====================== END OF  GET H FIELD INV =============================="
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


  SUBROUTINE copy_field_backward_AM_rz()
    USE fields, ONLY : el, er, et, bl, br, bt, jl, jr, jt
    USE fields, ONLY : el_c, er_c, et_c, bl_c, br_c, bt_c, jl_c, jr_c, jt_c
    USE omp_lib
    USE picsar_precision, ONLY: idp

    IMPLICIT NONE
    INTEGER(idp) :: ir, il, imode, irr ,ill
    INTEGER(idp) , dimension(3) :: lbound_r, ubound_r, lbound_p ,ubound_p

    lbound_r = LBOUND(el_c)
    lbound_p = LBOUND(el)
    ubound_r = UBOUND(el_c)
    ubound_p = UBOUND(el)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ir, il, imode, irr , ill) COLLAPSE(3)
    DO imode=lbound_r(3),ubound_r(3)
      DO il=lbound_r(2),ubound_r(2)
        DO ir=lbound_r(1),ubound_r(1)
          irr = ir - lbound_r(1) +lbound_p(1)
          ill = il - lbound_r(2) +lbound_p(2)
          el(irr, ill, imode)=el_c(ir, il, imode)
          er(irr, ill, imode)=er_c(ir, il, imode)
          et(irr, ill, imode)=et_c(ir, il, imode)
          bl(irr, ill, imode)=bl_c(ir, il, imode)
          br(irr, ill, imode)=br_c(ir, il, imode)
          bt(irr, ill, imode)=bt_c(ir, il, imode)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE copy_field_backward_AM_rz






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
   ! USE gpstd_solver
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
     CALL copy_field_backward_AM_rz
#endif
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE get_fields_AM_rz
END MODULE Hankel

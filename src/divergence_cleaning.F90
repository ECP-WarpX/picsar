
SUBROUTINE laser_field_correction
  USE PICSAR_precision
  USE shared_data , ONLY : nx, ny, dy, nmodes
  USE constants, ONLY: clight
  USE fields !, ONLY: l_staggered, em_h, ep_h, el_h, bm_h, bp_h, bl_h, nxguards, nyguards
  USE gpstd_solver
  USE hankel
  IMPLICIT NONE
  COMPLEX (cpx) :: ii
  INTEGER (idp) :: i,j,k, nfftx, nffty
  COMPLEX(cpx) , dimension (:) , allocatable :: filter_array
  COMPLEX(cpx), dimension (:) , allocatable :: kyc_inv, ky_true
  COMPLEX (cpx), dimension(:,:,:), allocatable :: w, inv_w
  REAL (num) :: prop_dir

#if defined(LIBRARY)
   nfftx = nx+2*nxguards+1
   nffty = ny+2*nyguards+1
#else
   !> When using picsar
   nfftx = nx+2*nxguards
   nffty = ny+2* nyguards
#endif
  ALLOCATE (ky_true(nffty))
  ALLOCATE (kyc_inv(nffty))
  ALLOCATE (w (nfftx, nffty, nmodes))
  ALLOCATE (inv_w (nfftx, nffty, nmodes))
  ALLOCATE (filter_array (nffty))
  ii= DCMPLX(0.0_num, 1.0_num)
  prop_dir= 1._num

  CALL get_Hfields
  !CALL compute_k_vec(l_staggered)
  CALL fftfreq(nffty, ky_true, dy)
  write(0,*) "max ", maxval(abs(el_h)) , maxval(abs(em_h)), maxval(abs(ep_h))  
  filter_array= (1._num- sin(0.5_num*ky_true*dy**2))*(1._num+sin(0.5_num*ky_true*dy**2))
  
  write(0,*) "kyc" , maxval(abs(kyc))
  write(0,*) "kyc" , maxval(abs(ky_true))
  !DO j= 1, nffty
  !  DO i= 1, nfftx
  !    em_h(i,j,:)= em_h(i,j,:)*filter_array(j)
  !    ep_h(i,j,:)= ep_h(i,j,:)* filter_array(j)
  !   END DO
  !END DO
  
  DO i=1, nffty
    IF (ky_true(i) .eq. 0.) THEN
     kyc_inv=0._num
    ELSE
     kyc_inv(i)=1./ky_true(i)
    END IF
  END DO
  
  
  !DO j=1, nffty
  !  DO i=1, nfftx
  !    el_h(i,j,:)=ii* krc(i,:)*(ep_h(i,j,:)-em_h(i,j,:))*kyc_inv(i)
  !  END DO
  !END DO  
  !w = clight* sqrt(kyc**2+krc**2)
  !w *= sign(kyc)* prop_dir
  
  DO k=1 , nmodes
    DO j=1, nffty
      DO i= 1, nfftx
        el_h(i,j,k)=ii* krc(i,k)*(ep_h(i,j,k)-em_h(i,j,k))*kyc_inv(j)
        w(i,j,k) =  prop_dir *clight* sqrt(ky_true(j)**2+krc(i,k)**2)
        w(i,j,k)= dsign (REAL(w(i,j,k), KIND=8), REAL(ky_true(j), KIND=8))
        IF (w(i,j,k) .eq. 0.) THEN     
          inv_w(i,j,k)= 0._num
        ELSE
          inv_w(i,j,k)= 1./w(i,j,k)
        END IF
        bp_h(i,j,k) = -ii*inv_w(i,j,k)*(ky_true(j)*ep_h(i,j,k)-0.5_num*ii*krc(i,k)*el_h(i,j,k))
        bm_h(i,j,k) = -ii*inv_w(i,j,k)*(-ky_true(j)*em_h(i,j,k) -0.5_num*ii*krc(i,k)*el_h(i,j,k))
        bl_h(i,j,k) = inv_w(i,j,k)* krc(i,k) *(ep_h(i,j,k)+em_h(i,j,k))
      END DO
    END DO 
  END DO
  
  write(0, *) "max b ", maxval(abs(bl_h)), maxval(abs(bm_h)), maxval(abs(bp_h))
  !bp_h  = -ii*inv_w*(kyc*ep_h-0.5_num*krc*el_h)
  !bm_h = -ii*inv_w*(-kyc*em_h -0.5_num*krc*el_h)
  !bl_h = inv_w* krc *(ep_h+em_h)
  
  
  CALL get_fields_AM_rz()
  
#if defined(LIBRARY)
  write(0, *) "max b ", maxval(abs(el_c)), maxval(abs(er_c)), maxval(abs(et_c)), maxval(abs(bl_c)), maxval(abs(br_c)), maxval(abs(bt_c))
#else
  write(0, *) "max b ", maxval(abs(el)), maxval(abs(er)), maxval(abs(et)), maxval(abs(bl)), maxval(abs(br)), maxval(abs(bt))
#endif
END SUBROUTINE laser_field_correction


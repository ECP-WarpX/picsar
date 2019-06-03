! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software, ! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
MODULE link_external_tools
  USE iso_c_binding
  CONTAINS 
  SUBROUTINE init_params_external(n1,n2,n3,d1,d2,d3,dtt,ng1,ng2,ng3,nor1,nor2,nor3,is_spec,&
      field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11, cdim) &
      BIND(C,name='init_params_picsar') 
    USE fastfft
    USE fields, ONLY: bx, bx_r, bxf, by, by_r, byf, bz, bz_r, bzf, ex, ex_r, exf,    &
      ey, ey_r, eyf, ez, ez_r, ezf, jx, jx_r, jxf, jy, jy_r, jyf, jz, jz_r, jzf,     &
      l_spectral, l_staggered, norderx, nordery, norderz, nxguards, nyguards,        &
      nzguards, rho_r, rhof, rhoold_r, rhooldf, xcoeffs, ycoeffs, zcoeffs
#if defined(FFTW)
    USE fourier_psaotd
#endif
    USE iso_c_binding
    USE params, ONLY: dt
    USE picsar_precision, ONLY: idp, isp, lp
    USE shared_data, ONLY: c_dim, dx, dy, dz, fftw_hybrid, fftw_mpi_transpose,       &
      fftw_threads_ok, fftw_with_mpi, nkx, nky, nkz, nx, ny, nz, p3dfft_flag,        &
      p3dfft_stride, rank
    IMPLICIT NONE 
    INTEGER(C_INT) , INTENT(IN) :: n1,n2,n3,ng1,ng2,ng3,nor1,nor2,nor3,cdim
    REAL(C_DOUBLE) , INTENT(INOUT), TARGET , DIMENSION(-ng3:n3+ng3,-ng2:n2+ng2,-ng1:n1+ng1) :: &
        field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11
    REAL(C_DOUBLE) , INTENT(IN) ::d1,d2,d3,dtt
    INTEGER(idp) :: imn, imx, jmn, jmx, kmn, kmx
    LOGICAL(C_BOOL)   , INTENT(IN)   :: is_spec
    LOGICAL(lp)                      :: l_stg
    INTEGER(isp)                     :: iret

    IF(rank==0) PRINT*, 'BEGIN INIT EXTERNAL'
    l_spectral  = LOGICAL(is_spec,lp) 
    fftw_with_mpi = .FALSE. 
    fftw_hybrid = .FALSE.
    fftw_mpi_transpose = .FALSE.
    l_staggered = .TRUE.
    fftw_threads_ok = .FALSE.
    CALL DFFTW_INIT_THREADS(iret)
    fftw_threads_ok = .TRUE.
    p3dfft_flag = .FALSE.
    p3dfft_stride = .FALSE.
    c_dim = INT(cdim,idp)
    nx = INT(n3,idp)
    ny = INT(n2,idp)
    nz = INT(n1,idp)
    nxguards = INT(ng3,idp) 
    nyguards = INT(ng2,idp)
    nzguards = INT(ng1,idp)
    dx = d3
    dy = d2
    dz = d1
    dt = dtt

    norderx = INT(nor3,idp)
    nordery = INT(nor2,idp)
    norderz = INT(nor1,idp) 

    IF(.NOT. l_spectral) THEN
      ex => field3
      ey => field2
      ez => field1
      bx => field6
      by => field5
      bz => field4

      jx => field9 
      jy => field8
      jz => field7
    ELSE
      ex_r => field3
      ey_r => field2
      ez_r => field1
      bx_r => field6
      by_r => field5
      bz_r => field4

      jx_r => field9
      jy_r => field8
      jz_r => field7
      rho_r =>field10
      rhoold_r =>field11

      nkx=(2*nxguards+nx+1)/2+1! Real To Complex Transform
      nky=(2*nyguards+ny+1)
      nkz=(2*nzguards+nz+1)

      IF(.NOT. ASSOCIATED(exf)) ALLOCATE(exf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(eyf)) ALLOCATE(eyf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(ezf)) ALLOCATE(ezf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(bxf)) ALLOCATE(bxf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(byf)) ALLOCATE(byf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(bzf)) ALLOCATE(bzf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(jxf)) ALLOCATE(jxf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(jyf)) ALLOCATE(jyf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(jzf)) ALLOCATE(jzf(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(rhof)) ALLOCATE(rhof(nkx, nky, nkz))
      IF(.NOT. ASSOCIATED(rhooldf)) ALLOCATE(rhooldf(nkx, nky, nkz))

    ENDIF
    IF(l_spectral) CALL init_plans_blocks
    IF(.NOT. l_spectral) THEN 
      ALLOCATE(xcoeffs(norderx/2))
      ALLOCATE(ycoeffs(nordery/2))
      ALLOCATE(zcoeffs(norderz/2))
      l_stg = .TRUE.
      CALL FD_weights_hvincenti(norderx,xcoeffs,l_stg)
      CALL FD_weights_hvincenti(nordery,ycoeffs,l_stg)
      CALL FD_weights_hvincenti(norderz,zcoeffs,l_stg)
      xcoeffs = dt/dx*xcoeffs
      ycoeffs = dt/dy*ycoeffs
      zcoeffs = dt/dz*zcoeffs
    ENDIF 
    IF(rank==0) PRINT*, 'END INIT EXTERNAL'
  END SUBROUTINE init_params_external


  SUBROUTINE init_params_external_RZ(nr, nl,nmodes_in, imode,dr,dl,dtt,ngr,ngl,norr,norl,&
      field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11) &
      BIND(C,name='init_params_picsar_AM')
    USE fastfft
    USE fields, ONLY: bl_c, br_c,bt_c, bl_f, bp_f, bm_f, bl_h, bm_h, bp_h, bl_h_inv, bm_h_inv, bp_h_inv,  &
      el_c, er_c,et_c, el_f, ep_f, em_f, el_h, em_h, ep_h, el_h_inv, em_h_inv, ep_h_inv,  &
      jl_c, jr_c,jt_c, jl_f, jp_f, jm_f, jl_h, jm_h, jp_h,  &
      rho_f, rhoold_f, rho_h, rhoold_h, invM, invM1, invM_1, Ma, Ma1, Ma_1,&
      l_spectral, l_staggered, l_AM_RZ, norderx, nordery, nxguards, nyguards,        &
      rhof, rho_c,rhoold_c, rhooldf, xcoeffs, ycoeffs, zcoeffs
#if defined(FFTW)
    USE fourier_psaotd
    USE hankel 
#endif
    USE iso_c_binding
    USE params, ONLY: dt
    USE picsar_precision, ONLY: idp, isp, lp
    USE shared_data, ONLY: c_dim, dx, dy, fftw_hybrid, fftw_mpi_transpose,       &
      fftw_threads_ok, fftw_with_mpi, nkx, nky, nx, ny, nmodes, p3dfft_flag,        &
      p3dfft_stride, rank
    IMPLICIT NONE
    INTEGER(C_INT) , INTENT(IN) :: nr,nl,nmodes_in,imode,ngr,ngl,norr,norl
    COMPLEX(C_DOUBLE) , INTENT(INOUT), TARGET , DIMENSION(-ngr:nr+ngr,-ngl:nl+ngl,0:nmodes_in-1) :: &
         field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11
    REAL(C_DOUBLE) , INTENT(IN) ::dr,dl,dtt
    INTEGER(idp) :: imn, imx, jmn, jmx, kmn, kmx
    LOGICAL(lp)                      :: l_stg
    INTEGER(isp)                     :: iret

    IF(rank==0) PRINT*, 'BEGIN INIT EXTERNAL'
    l_spectral  = .TRUE.
    fftw_with_mpi = .FALSE.
    fftw_hybrid = .FALSE.
    fftw_mpi_transpose = .FALSE.
    l_staggered = .FALSE.
    fftw_threads_ok = .FALSE.
    l_AM_RZ= .TRUE.
    l_spectral= .TRUE.
    C_DIM= 3
    CALL DFFTW_INIT_THREADS(iret)
    fftw_threads_ok = .TRUE.
    p3dfft_flag = .FALSE.
    p3dfft_stride = .FALSE.
    c_dim = INT(3,idp)
    nx = INT(nr,idp)
    ny = INT(nl,idp)
    nmodes = INT(nmodes,idp)
    nxguards = INT(ngr,idp)
    nyguards = INT(ngl,idp)
    dx = dr
    dy = dl
    dt = dtt
    norderx = INT(norr,idp)
    nordery = INT(norl,idp)

      
   IF ((l_spectral) .AND. (l_AM_RZ)) THEN 
      el_c => field3
      er_c => field2
      et_c => field1
      bl_c=> field6
      br_c => field5
      bt_c => field4

      jl_c => field9
      jr_c=> field8
      jt_c => field7
      rho_c =>field10
      rhoold_c =>field11

      nkx=(2*nxguards+nx+1)! Real To Complex Transform
      nky=(2*nyguards+ny+1)
      nmodes= nmodes_in

      IF(.NOT. ASSOCIATED(el_f)) THEN
         write(0,*) "before allocation"
         ALLOCATE(el_f(nkx, nky, nmodes))
         write(0,*) "After allocation"
      ENDIF
      IF(.NOT. ASSOCIATED(em_f)) ALLOCATE(em_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(ep_f)) ALLOCATE(ep_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bl_f)) ALLOCATE(bl_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bm_f)) ALLOCATE(bm_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bp_f)) ALLOCATE(bp_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(jl_f)) ALLOCATE(jl_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(jm_f)) ALLOCATE(jm_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(jp_f)) ALLOCATE(jp_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(rho_f)) ALLOCATE(rho_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(rhoold_f)) ALLOCATE(rhoold_f(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(el_h)) ALLOCATE(el_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(em_h)) ALLOCATE(em_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(ep_h)) ALLOCATE(ep_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bl_h)) ALLOCATE(bl_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bm_h)) ALLOCATE(bm_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bp_h)) ALLOCATE(bp_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(jl_h)) ALLOCATE(jl_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(jm_h)) ALLOCATE(jm_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(jp_h)) ALLOCATE(jp_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(rho_h)) ALLOCATE(rho_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(rhoold_h)) ALLOCATE(rhoold_h(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(el_h_inv)) ALLOCATE(el_h_inv(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(em_h_inv)) ALLOCATE(em_h_inv(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(ep_h_inv)) ALLOCATE(ep_h_inv(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bl_h_inv)) ALLOCATE(bl_h_inv(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bm_h_inv)) ALLOCATE(bm_h_inv(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(bp_h_inv)) ALLOCATE(bp_h_inv(nkx, nky, nmodes))
      IF(.NOT. ASSOCIATED(invM)) ALLOCATE(invM(nkx, nkx))
      IF(.NOT. ASSOCIATED(invM_1)) ALLOCATE(invM_1(nkx, nkx))
      IF(.NOT. ASSOCIATED(invM1)) ALLOCATE(invM1(nkx, nkx))
      IF(.NOT. ASSOCIATED(Ma)) ALLOCATE(Ma(nkx, nkx))
      IF(.NOT. ASSOCIATED(Ma1)) ALLOCATE(Ma1(nkx, nkx))
      IF(.NOT. ASSOCIATED(Ma_1)) THEN
        ALLOCATE(Ma_1(nkx, nkx))
        CALL init_plans_blocks_rz    
        CALL init_rz_fields_coupling   
      ENDIF
    ENDIF

    !CALL init_plans_blocks_RZ
    IF(rank==0) PRINT*, 'END INIT EXTERNAL'
  END SUBROUTINE init_params_external_RZ

SUBROUTINE init_rz_fields_coupling
  USE PICSAR_precision
  USE fields

  el_h= DCMPLX(0.0_num, 0.0_num)
  em_h= DCMPLX(0.0_num, 0.0_num)
  ep_h= DCMPLX(0.0_num, 0.0_num)
  bl_h= DCMPLX(0.0_num, 0.0_num)
  bm_h= DCMPLX(0.0_num, 0.0_num)
  bp_h= DCMPLX(0.0_num, 0.0_num)
  jl_h= DCMPLX(0.0_num, 0.0_num)
  jm_h= DCMPLX(0.0_num, 0.0_num)
  jp_h= DCMPLX(0.0_num, 0.0_num)
  rho_h= DCMPLX(0.0_num, 0.0_num)
  rhoold_h= DCMPLX(0.0_num, 0.0_num)
  el_h_inv= DCMPLX(0.0_num, 0.0_num)
  em_h_inv= DCMPLX(0.0_num, 0.0_num)
  ep_h_inv= DCMPLX(0.0_num, 0.0_num)
  bl_h_inv= DCMPLX(0.0_num, 0.0_num)
  bm_h_inv= DCMPLX(0.0_num, 0.0_num)
  bp_h_inv= DCMPLX(0.0_num, 0.0_num)
  el_f= DCMPLX(0.0_num, 0.0_num)
  em_f= DCMPLX(0.0_num, 0.0_num)
  ep_f=  DCMPLX(0.0_num, 0.0_num)
  bl_f= DCMPLX(0.0_num, 0.0_num)
  bm_f= DCMPLX(0.0_num, 0.0_num)
  bp_f= DCMPLX(0.0_num, 0.0_num)

END SUBROUTINE init_rz_fields_coupling
  SUBROUTINE evec3d_push_norder(ex, ey, ez, bx, by, bz, jx, jy, jz, dt, dtsdx,  &
  dtsdy, dtsdz, nx, ny, nz, norderx, nordery, norderz, nxguard, nyguard,nzguard)
  USE omp_lib
  USE picsar_precision, ONLY: idp, num
  INTEGER(idp), INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard
  INTEGER(idp), INTENT(IN) :: norderx, nordery, norderz
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: Jx, Jy, Jz
  REAL(num), INTENT(IN) :: dt, dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  INTEGER(idp) :: i, j, k, l, ist,nxs,nys,nzs
  ist = 1
  nxs =nxguard
  nys =nyguard
  nzs =nzguard

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(3)
  ! advance Ex
  DO l = -nzs, nz+nzs-1
    DO k = -nys, ny+nys-1
      DO j = -nxs, nx+nxs
        Ex(j, k, l) = Ex(j, k, l) - dt  * Jx(j, k, l)
        DO i = 1,nordery/2 ! MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          IF((k+i .GT. ny+nyguard) .OR.(k-i+ist .LT. -nyguard)) CYCLE
          Ex(j, k, l) = Ex(j, k, l) + dtsdy(i) * (Bz(j, k+i, l)   - Bz(j, k-i+ist, l  &
          ))
        ENDDO
        DO i = 1,norderz/2 ! MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          IF((l+i .GT. nz+nzguard) .OR.(l-i+ist .LT. -nzguard)) CYCLE
          Ex(j, k, l) = Ex(j, k, l) - dtsdz(i) * (By(j, k, l+i)   - By(j, k,      &
          l-i+ist))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  !$OMP DO COLLAPSE(3)
  ! advance Ey
  DO l = -nzs, nz+nzs-1
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs-1
        Ey(j, k, l) = Ey(j, k, l) - dt  * Jy(j, k, l)
        DO i = 1,MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          IF((j+i .GT. nx+nxguard) .OR.(j-i+ist .LT. -nxguard)) CYCLE
          Ey(j, k, l) = Ey(j, k, l) - dtsdx(i) * (Bz(j+i, k, l)   - Bz(j-i+ist, k,    &
          l))
        ENDDO
        DO i = 1,MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          IF((l+i .GT. nz+nzguard) .OR.(l-i+ist .LT. -nzguard)) CYCLE
          Ey(j, k, l) = Ey(j, k, l) + dtsdz(i) * (Bx(j, k, l+i)   - Bx(j, k,      &
          l-i+ist))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
  ! advance Ez
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys-1
      DO j = -nxs, nx+nxs-1
        Ez(j, k, l) = Ez(j, k, l) - dt  * Jz(j, k, l)
        DO i = 1,MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          IF((j+i .GT. nx+nxguard) .OR.(j-i+ist .LT. -nxguard)) CYCLE
          Ez(j, k, l) = Ez(j, k, l) + dtsdx(i) * (By(j+i, k, l) - By(j-i+ist, k, l))
        ENDDO
        DO i = 1,MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          IF((k+i .GT. ny+nyguard) .OR.(k-i+ist .LT. -nyguard)) CYCLE
          Ez(j, k, l) = Ez(j, k, l) - dtsdy(i) * (Bx(j, k+i, l) - Bx(j, k-i+ist, l))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN
END SUBROUTINE evec3d_push_norder

SUBROUTINE bvec3d_push_norder(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx,  &
  ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard)
  USE picsar_precision, ONLY: idp, num
  INTEGER(idp)          :: nx, ny, nz, nxguard, nyguard, nzguard,           &
  norderx, nordery, norderz
  REAL(num), INTENT(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN) :: dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  INTEGER(idp)          :: i, j, k, l,nxs,nys,nzs, ist

  ist = 1
  nxs =nxguard
  nys =nyguard
  nzs =nzguard
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(3)
  DO l = -nzs+1, nz+nzs-1
    DO k = -nys+1, ny+nys-1
      DO j = -nxs, nx+nxs-1
        DO i = 1,MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          IF((k+i-ist .GT. ny+nyguard) .OR.(k-i .LT. -nyguard)) CYCLE
          Bx(j, k, l) = Bx(j, k, l) - dtsdy(i) * (Ez(j, k+i-ist, l  ) - Ez(j, k-i,    &
          l))
        ENDDO
        DO i = 1,MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          IF((l+i-ist .GT. nz+nzguard) .OR.(l-i .LT. -nzguard)) CYCLE
          Bx(j, k, l) = Bx(j, k, l) + dtsdz(i) * (Ey(j, k, l+i-ist) - Ey(j, k, l-i))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! advance By
  !$OMP DO COLLAPSE(3)
  DO l = -nzs+1, nz+nzs-1
    DO k = -nys, ny+nys-1
      DO j = -nxs+1, nx+nxs-1
        DO i = 1,MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          IF((j+i-ist .GT. nx+nxguard) .OR.(j-i .LT. -nxguard)) CYCLE
          By(j, k, l) = By(j, k, l) + dtsdx(i) * (Ez(j+i-ist, k, l  ) - Ez(j-i, k,    &
          l))
        ENDDO
        DO i = 1,MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          IF((l+i-ist .GT. nz+nzguard) .OR.(l-i .LT. -nzguard)) CYCLE
          By(j, k, l) = By(j, k, l) - dtsdz(i) * (Ex(j, k, l+i-ist) - Ex(j, k, l-i))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! advance Bz
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs-1
    DO k = -nys+1, ny+nys-1
      DO j = -nxs+1, nx+nxs-1
        DO i = 1,MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          IF((j+i-ist .GT. nx+nxguard) .OR.(j-i .LT. -nxguard)) CYCLE
          Bz(j, k, l) = Bz(j, k, l) - dtsdx(i) * (Ey(j+i-ist, k, l) - Ey(j-i, k, l))
        ENDDO
        DO i = 1,MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          IF((k+i-ist .GT. ny+nyguard) .OR.(k-i .LT. -nyguard)) CYCLE
          Bz(j, k, l) = Bz(j, k, l) + dtsdy(i) * (Ex(j, k+i-ist, l) - Ex(j, k-i, l))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN

END SUBROUTINE bvec3d_push_norder

SUBROUTINE solve_maxwell_fdtd_pxr() bind(C,name='solve_maxwell_fdtd_pxr')
  USE fields, ONLY: bx, by, bz, ex, ey, ez, jx, jy, jz, norderx, nordery, norderz,   &
    nxguards, nyguards, nzguards, xcoeffs, ycoeffs, zcoeffs
  USE params, ONLY: dt
  USE shared_data, ONLY: nx, ny, nz

  CALL evec3d_push_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,dt,xcoeffs,ycoeffs,zcoeffs,&
        nx, ny, nz, norderx, nordery, norderz, nxguards,nyguards,nzguards)
  CALL bvec3d_push_norder(ex,ey,ez,bx,by,bz,xcoeffs,ycoeffs,zcoeffs,&
        nx, ny, nz, norderx, nordery, norderz, nxguards,nyguards,nzguards)

END SUBROUTINE solve_maxwell_fdtd_pxr
END MODULE

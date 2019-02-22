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
      field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11, cdim, &
      is_pml_x, is_pml_y, is_pml_z)
      BIND(C,name='init_params_picsar') 
    USE fastfft
    USE fields, ONLY: ezf, ez_r, ez, jx_r, jxf, nordery, ey_r, rhooldf, l_staggered, &
      ex_r, rhof, bx_r, jz, by_r, bxf, l_spectral, zcoeffs, rho_r, xcoeffs, bz,      &
      nzguards, nxguards, norderz, bz_r, nyguards, jzf, jy, jx, ex, bx, jz_r, jy_r,  &
      by, ycoeffs, eyf, rhoold_r, jyf, norderx, byf, bzf, exf, ey
#if defined(FFTW)
    USE fourier_psaotd
#endif
    USE iso_c_binding
    USE params, ONLY: dt
    USE picsar_precision, ONLY: idp, isp, lp
    USE shared_data, ONLY: nz, ny, fftw_with_mpi, nx, nkx, p3dfft_stride, nky,       &
      p3dfft_flag, fftw_threads_ok, dx, c_dim, nkz, fftw_mpi_transpose, dy, rank,    &
      fftw_hybrid, dz, absorbing_bcs_x, absorbing_bcs_y, absorbing_bcs_z,            &
      absorbing_bcs
    IMPLICIT NONE 
    INTEGER(C_INT) , INTENT(IN) :: n1,n2,n3,ng1,ng2,ng3,nor1,nor2,nor3,cdim
    REAL(C_DOUBLE) , INTENT(INOUT), TARGET , DIMENSION(-ng3:n3+ng3,-ng2:n2+ng2,-ng1:n1+ng1) :: &
        field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11
    REAL(C_DOUBLE) , INTENT(IN) ::d1,d2,d3,dtt
    INTEGER(idp) :: imn, imx, jmn, jmx, kmn, kmx
    LOGICAL(C_BOOL)   , INTENT(IN)   :: is_spec
    LOGICAL(lp)                      :: l_stg
    INTEGER(isp)                     :: iret
    LOGICAL(C_BOOL),    INTENT(IN)   :: is_pml_x, is_pml_y, is_pml_z     
    INTEGER(isp) , DIMENSIOn(1:2*cdim) :: mpi_neighbors 
    

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

    absorbing_bcs_x = LOGICAL(is_pml_z,lp)
    absorbing_bcs_y = LOGICAL(is_pml_y,lp)
    absorbing_bcs_z = LOGICAL(is_pml_x,lp)
 
    IF(.NOT. l_spectral) THEN
      absorbing_bcs_x = .FALSE._lp
      absorbing_bcs_y = .FALSE._lp
      absorbing_bcs_z = .FALSE._lp
    ENDIF
    IF(absorbing_bcs_x .OR. absorbing_bcs_y .OR. absorbing_bcs_y) THEN
      absorbing_bcs = .TRUE._lp
    ELSE
      absorbing_bcs = .FALSE._lp
    ENDIF
    g_spectral = absorbing_bcs

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
      ex => field3
      ey => field2
      ez => field1
      bx => field6
      by => field5
      bz => field4
      jx => field9
      jy => field8
      jz => field7
      rho =>field10
      rhoold =>field11


      nkx=(2*nxguards+nx+1)/2+1! Real To Complex Transform
      nky=(2*nyguards+ny+1)
      nkz=(2*nzguards+nz+1)

      IF(.NOT. g_spectral) THEN
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
    ENDIF
    dx = -dx; dy = -dy; dz = -dz
    IF(l_spectral) THEN
      CALL init_plans_blocks()
      IF(absoring_bcs) THEN
        ALLOCATE(exy(nx,ny,z))
        ALLOCATE(exz(nx,ny,z))
        ALLOCATE(eyx(nx,ny,z))
        ALLOCATE(eyz(nx,ny,z))
        ALLOCATE(ezx(nx,ny,z))
        ALLOCATE(ezy(nx,ny,z))
        ALLOCATE(bxy(nx,ny,z))
        ALLOCATE(bxz(nx,ny,z))
        ALLOCATE(byx(nx,ny,z))
        ALLOCATE(byz(nx,ny,z))
        ALLOCATE(bzx(nx,ny,z))
        ALLOCATE(bzy(nx,ny,z))

        exy_r => exy
        exz_r => exz
        eyx_r => eyx
        eyz_r => eyz
        ezx_r => ezx
        ezy_r => ezy
        bxy_r => bxy
        bxz_r => bxz
        byx_r => byx
        byz_r => byz
        bzx_r => bzx
        bzy_r => bzy

        dx = -dx; dy = -dy; dz = -dz
        CALL init_pml_arrays()
        dx = -dx; dy = -dy; dz = -dz

      ENDIF
    ELSE IF (.NOT. l_spectral) THEN 
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


  SUBROUTINE push_psatd_ebfield_picsar() BIND(C,name="push_psatd_ebfield_picsar")
    USE field_boundary
    USE shared_data , ONLY : absorbing_bcs

  
    CALL push_psatd_ebfield()
    CALL efield_bcs
    CALL bfield_bcs
    IF(absorbing_bcs) THEN
      CALL field_damping_bcs()
      CALL merge_fields()
    ENDIF



  END SUBROUTINE push_psatd_ebfield_picsar

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

SUBROUTINE solve_maxwell_fdtd_picsar() bind(C,name='solve_maxwell_fdtd_picsar')
  USE fields, ONLY: ez, nordery, jz, zcoeffs, xcoeffs, bz, nzguards, nxguards,       &
    norderz, nyguards, jy, jx, ex, bx, by, ycoeffs, norderx, ey
  USE params, ONLY: dt
  USE shared_data, ONLY: nz, ny, nx

  CALL evec3d_push_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,dt,xcoeffs,ycoeffs,zcoeffs,&
        nx, ny, nz, norderx, nordery, norderz, nxguards,nyguards,nzguards)
  CALL bvec3d_push_norder(ex,ey,ez,bx,by,bz,xcoeffs,ycoeffs,zcoeffs,&
        nx, ny, nz, norderx, nordery, norderz, nxguards,nyguards,nzguards)

END SUBROUTINE solve_maxwell_fdtd_picsar
END MODULE

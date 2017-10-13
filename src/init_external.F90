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
!
SUBROUTINE init_params_external(n1,n2,n3,d1,d2,d3,dtt,ng1,ng2,ng3,nor1,nor2,nor3,is_spec,&
        field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11) &
        bind(C,name='init_params_picsar') 
    USE iso_c_binding 
    USE params
    USE shared_data
    USE constants
    USE picsar_precision
    USE fields 
#if defined(FFTW)
    USE fourier_psaotd
    USE fourier
#endif 
    IMPLICIT NONE 
    INTEGER(C_INT) , INTENT(IN) :: n1,n2,n3,ng1,ng2,ng3,nor1,nor2,nor3
    REAL(C_DOUBLE) , INTENT(INOUT), TARGET , DIMENSION(1:2*ng3+n3+1,1:2*ng2+n2,1:2*ng1+n1) :: &
        field1,field2,field3,field4,field5,field6,field7,field8,field9,field10,field11
    REAL(C_DOUBLE) , INTENT(IN) ::d1,d2,d3,dtt
    INTEGER(idp) :: imn, imx, jmn, jmx, kmn, kmx
    LOGICAL(C_BOOL)   , INTENT(IN)   :: is_spec

    l_spectral  = LOGICAL(is_spec,lp) 
    g_spectral  = .FALSE.
    fftw_with_mpi = .FALSE. 
    fftw_hybrid = .FALSE.
    hybrid_2 = .FALSE.
    fftw_mpi_transpose = .FALSE.

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
    norderx = REAL(nor3,idp)
    nordery = REAL(nor2,idp)
    norderz = REAL(nor1,idp) 

    ex => field1
    ey => field2
    ez => field3
    bx => field4
    by => field5
    bz => field6
    jx => field7 
    jy => field8
    jz => field9
    rho => field10
    rhoold =>field11
    
    nkx=(2*nxguards+nx)/2+1! Real To Complex Transform
    nky=(2*nyguards+ny)
    nkz=(2*nzguards+nz)
    IF (l_spectral) THEN
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
      imn=-nxguards;imx=nx+nxguards-1
      jmn=-nyguards;jmx=ny+nyguards-1
      kmn=-nzguards;kmx=nz+nzguards-1
      IF(.NOT. ASSOCIATED(ex_r)) ALLOCATE(ex_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(ey_r)) ALLOCATE(ey_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(ez_r)) ALLOCATE(ez_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(bx_r)) ALLOCATE(bx_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(by_r)) ALLOCATE(by_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(bz_r)) ALLOCATE(bz_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(jx_r)) ALLOCATE(jx_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(jy_r)) ALLOCATE(jy_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(jz_r)) ALLOCATE(jz_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(rho_r)) ALLOCATE(rho_r(imn:imx, jmn:jmx, kmn:kmx))
      IF(.NOT. ASSOCIATED(rhoold_r)) ALLOCATE(rhoold_r(imn:imx, jmn:jmx, kmn:kmx))
    ENDIF
END SUBROUTINE

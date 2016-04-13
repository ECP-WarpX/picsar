! ______________________________________________________________________________
! 
! CURRENT_DEPOSITION.F90
!
! List of subroutines:
! - pxrdepose_currents_on_grid_jxjyjz
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2
! - pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq
! - pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp
! - pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq
! 
! Classical non-optimized:
! - depose_jxjyjz_scalar_1_1_1
! - depose_jxjyjz_scalar_2_2_2
! - depose_jxjyjz_scalar_3_3_3
!
! Classical parallel/vectorized
! - depose_jxjyjz_vecHVv2_1_1_1
! - depose_jxjyjz_vecHVv2_2_2_2
! - depose_jxjyjz_vecHVv2_3_3_3
! - depose_jxjyjz_vecHVv3_3_3_3
!
! Subroutines with no reduction inside
! - depose_jxjyjz_vecHV_vnr_1_1_1
! - depose_jxjyjz_vecHV_vnr_2_2_2
! - depose_jxjyjz_vecHV_vnr_3_3_3
! - depose_jxjyjz_vecHV_vnr_3_3_3
!
! Esirkepov:
! - depose_jxjyjz_esirkepov_1_1_1
! - depose_jxjyjz_esirkepov_2_2_2
! - depose_jxjyjz_esirkepov_3_3_3
! - warp_depose_jxjyjz_esirkepov_n
!
! - depose_jxjyjz_tile_vecHV_1_1_1
! - depose_jxjyjz_tile_vecHV_2_2_2
! - depose_jxjyjz_tile_vecHV_3_3_3
!
! - current_reduction
!
! ______________________________________________________________________________


! _____________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz
  ! 
  ! Main subroutine called in submain
  ! ___________________________________________________________________
  USE fields
  USE shared_data
  USE params
  USE time_stat
#if defined(PROFILING) && PROFILING==2      
  USE ITT_SDE_FORTRAN                       
#endif                                   
  
  
  IMPLICIT NONE 
  REAL(num) :: tdeb

  ! Interfaces for func_order
  INTERFACE
  
    ! ____________________________________________________________________________________
    ! Classical current deposition - non optimized - order 1
    SUBROUTINE depose_jxjyjz_scalar_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin           
    END SUBROUTINE
    
    ! Classical current deposition - non optimized - order 2
    SUBROUTINE depose_jxjyjz_scalar_2_2_2(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin           
    END SUBROUTINE

    ! Classical current deposition - non optimized - order 3
    SUBROUTINE depose_jxjyjz_scalar_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin           
    END SUBROUTINE
  
    SUBROUTINE depose_jxjyjz_vecHVv2_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    END SUBROUTINE depose_jxjyjz_vecHVv2_1_1_1    
    SUBROUTINE depose_jxjyjz_vecHVv2_2_2_2(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    END SUBROUTINE depose_jxjyjz_vecHVv2_2_2_2  
    SUBROUTINE depose_jxjyjz_vecHVv2_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    END SUBROUTINE depose_jxjyjz_vecHVv2_3_3_3
    SUBROUTINE depose_jxjyjz_vecHVv3_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    END SUBROUTINE depose_jxjyjz_vecHVv3_3_3_3 

     ! Interface for subroutine with no reduction - classical deposition order 1
      SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
      	USE constants
      	IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
        REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
        REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
        REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        INTEGER(idp) :: ncx, ncy, ncz
      END SUBROUTINE

     ! Interface for subroutine with no reduction - classical deposition order 2
      SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
      	USE constants
      	IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
        REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
        REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
        REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        INTEGER(idp) :: ncx, ncy, ncz
      END SUBROUTINE

     ! Interface for subroutine with no reduction - classical deposition order 3    
      SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
      	USE constants
      	IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
        REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
        REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
        REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        INTEGER(idp) :: ncx, ncy, ncz
      END SUBROUTINE

    
    ! Esirkepov at any order   
    SUBROUTINE pxr_depose_jxjyjz_esirkepov_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
           nox,noy,noz,l_particles_weight,l4symtry) !#do not parse

      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
      LOGICAL(idp) :: l_particles_weight,l4symtry
    END SUBROUTINE

    ! Esirkepov order 1    
    SUBROUTINE depose_jxjyjz_esirkepov_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
           nox,noy,noz,l_particles_weight,l4symtry) !#do not parse

      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
      LOGICAL(idp) :: l_particles_weight,l4symtry
    END SUBROUTINE

    ! Esirkepov order 2
    SUBROUTINE depose_jxjyjz_esirkepov_2_2_2(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
           nox,noy,noz,l_particles_weight,l4symtry) !#do not parse

      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
      LOGICAL(idp) :: l_particles_weight,l4symtry
    END SUBROUTINE

    ! Esirkepov order 3
    SUBROUTINE depose_jxjyjz_esirkepov_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
           nox,noy,noz,l_particles_weight,l4symtry) !#do not parse

      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
      LOGICAL(idp) :: l_particles_weight,l4symtry
    END SUBROUTINE
    
    ! Esirkepov optimized order 1
    SUBROUTINE depose_jxjyjz_esirkepov_vecHV_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
           nox,noy,noz,l_particles_weight,l4symtry)  !#do not parse

      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
      LOGICAL(idp) :: l_particles_weight,l4symtry
    END SUBROUTINE
    
  END INTERFACE
  
#if defined(DEBUG)
  WRITE(0,*) "Depose_currents_on_grid: start"
#endif

  tdeb=MPI_WTIME()

#if PROFILING==2              
  CALL start_collection()     
#endif                        

  jx = 0.0_num
  jy = 0.0_num
  jz = 0.0_num

  ! Current deposition branches

  ! _______________________________________________________
  ! Classical current deposition, non-optimized/no tiling
  
  IF (currdepo.EQ.5) THEN

 
    IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
  
    CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(depose_jxjyjz_scalar_3_3_3, &
    jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(depose_jxjyjz_scalar_2_2_2, &
    jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)

    ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(depose_jxjyjz_scalar_1_1_1, &
    jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ENDIF
  ! _______________________________________________________
  ! Classical current deposition, non-optimized/tiling
    
  ELSE IF (currdepo.EQ.4) THEN
 
    IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
  
    CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_scalar_3_3_3, &
    jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_scalar_2_2_2, &
    jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)

    ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_scalar_1_1_1, &
    jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ENDIF
 
  ! _______________________________________________________
  ! Classical current deposition, parallel, vectorized
  ELSE IF (currdepo.EQ.3) THEN
  
    IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN
  
      ! Old version with reduction for each species
      !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_vecHVv3_3_3_3, &
      !jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)

      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( &
           depose_jxjyjz_vecHV_vnr_3_3_3, &
           jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN

      ! Old version with reduction for each species
      !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_vecHVv2_2_2_2, &
      !jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)

      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( &
           depose_jxjyjz_vecHV_vnr_2_2_2, &
           jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ELSE IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN

      ! Old version with reduction for each species
      !CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(depose_jxjyjz_vecHVv2_1_1_1, &
      !jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
           
      CALL pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2( &
           depose_jxjyjz_vecHV_vnr_1_1_1, &
           jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt)
    
    ELSE
    
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(pxr_depose_jxjyjz_esirkepov_n, &
           jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	         nox,noy,noz,dx,dy,dz,dt)    
    
    ENDIF
  ! _______________________________________________________    
  ! Esirkepov sequential version   
  ELSE IF (currdepo.EQ.2) THEN
  
    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	nox,noy,noz,dx,dy,dz,dt)
  ! _______________________________________________________
  ! Esirkepov tiling version 
  ELSE IF (currdepo.EQ.1) THEN

    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(pxr_depose_jxjyjz_esirkepov_n, &
         jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)
  ! _______________________________________________________	
	! Default - Esirkepov parallel version with OPENMP/tiling and optimizations
  ELSE
    
    ! Order 1
    IF ((nox.eq.1).AND.(noy.eq.1).AND.(noz.eq.1)) THEN

    !CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz_esirkepov_lin_1_1_1, &
    !     jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	  !     nox,noy,noz,dx,dy,dz,dt)

    CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz_esirkepov_1_1_1, &
         jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)

    ! Order 2
    ELSE IF ((nox.eq.2).AND.(noy.eq.2).AND.(noz.eq.2)) THEN

      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz_esirkepov_2_2_2, &
         jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)
	       
    ! Order 3
    ELSE IF ((nox.eq.3).AND.(noy.eq.3).AND.(noz.eq.3)) THEN    

      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(depose_jxjyjz_esirkepov_3_3_3, &
         jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)    
    
    ! Arbitrary order
    ELSE
  
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(pxr_depose_jxjyjz_esirkepov_n, &
         jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
	       nox,noy,noz,dx,dy,dz,dt)
	  END IF
	   
  ENDIF

!!! --- Stop Vtune analysis
#if PROFILING==2                     
  CALL stop_collection()            
#endif                               

  localtimes(3)=localtimes(3)+(MPI_WTIME()-tdeb)

#if defined(DEBUG)
  WRITE(0,*) "Depose_current_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz

!===============================================================================
! Deposit current in each tile with the classical method
! OpenMP version. Avoids conflict while reducing tile currents in the global 
! current array. 
!===============================================================================
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp(func_order,jxg,jyg,jzg,nxx,nyy,nzz,nxjguard,nyjguard,nzjguard, &
	noxx,noyy,nozz,dxx,dyy,dzz,dtt)
  USE particles
  USE constants
  USE tiling
  USE omp_lib
  USE timing
  USE time_stat
  IMPLICIT NONE

    ! Interfaces for func_order
    INTERFACE
      SUBROUTINE func_order(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, & !#do not parse
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      	USE constants
      	IMPLICIT NONE
      	INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      	REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      	REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      	REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      	REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
        REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
      END SUBROUTINE func_order
    END INTERFACE
  
  ! Parameters
  INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
  REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  INTEGER(idp) :: ispecies, ix, iy, iz, count
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER :: curr_tile
  TYPE(grid_tile), POINTER :: currg
  REAL(num) :: tdeb, tend
  INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(idp) :: isdeposited=.FALSE.

  !$OMP PARALLEL DEFAULT(NONE)                                                              &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguard,nyjguard,              &
  !$OMP nzjguard,dxx,dyy,dzz,dtt,jxg,jyg,jzg,noxx,noyy,nozz,aofgrid_tiles)                  &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,currg, curr_tile,count,jmin,jmax,kmin,kmax,lmin,     &
  !$OMP lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,nxc,nyc,nzc, nxjg, nyjg, nzjg, isdeposited)
  !! Current deposition
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
          	  curr => species_parray(1)
              curr_tile=>curr%array_of_tiles(ix,iy,iz)
              nxjg=curr_tile%nxg_tile
              nyjg=curr_tile%nyg_tile
              nzjg=curr_tile%nzg_tile
              jmin=curr_tile%nx_tile_min
          	  jmax=curr_tile%nx_tile_max
              kmin=curr_tile%ny_tile_min
              kmax=curr_tile%ny_tile_max
              lmin=curr_tile%nz_tile_min
              lmax=curr_tile%nz_tile_max
              nxc=curr_tile%nx_cells_tile
              nyc=curr_tile%ny_cells_tile
              nzc=curr_tile%nz_cells_tile         
			        currg=>aofgrid_tiles(ix,iy,iz)
              currg%jxtile=0.
              currg%jytile=0.
              currg%jztile=0.!jzg(jmin:jmax,kmin:kmax,lmin:lmax)
              isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
             	    curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  count=curr_tile%np_tile(1)
                  IF (count .EQ. 0) THEN 
                  	CYCLE
                  ELSE 
                  	isdeposited=.TRUE.
                  ENDIF 
                  ! Depose current in jtile
                  CALL func_order(currg%jxtile,currg%jytile,currg%jztile,count,&
                  curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						           &
                  curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,  			   &
                  curr_tile%part_gaminv,  &
                  curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,     &
                  curr_tile%z_grid_tile_min,dtt,dxx,dyy,dzz,nxc,nyc,nzc,                                     &
                  nxjg,nyjg,nzjg) 
              END DO! END LOOP ON SPECIES
              IF (isdeposited) THEN
              	jxg(jmin:jmax,kmin:kmax,lmin:lmax)=jxg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jxtile(0:nxc,0:nyc,0:nzc)
              	jyg(jmin:jmax,kmin:kmax,lmin:lmax)=jyg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jytile(0:nxc,0:nyc,0:nzc)
              	jzg(jmin:jmax,kmin:kmax,lmin:lmax)=jzg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jztile(0:nxc,0:nyc,0:nzc)
              ENDIF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
          	isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
                  curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  count=curr_tile%np_tile(1)
                  IF (count .GT. 0) isdeposited=.TRUE.  
              END DO
              IF (isdeposited) THEN 
              	currg=>aofgrid_tiles(ix,iy,iz)
                  curr => species_parray(1)
             		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                  kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                  lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                  nxjg=curr_tile%nxg_tile
                  nyjg=curr_tile%nyg_tile
                  nzjg=curr_tile%nzg_tile
                  jminc=jmin-nxjg; jmaxc=jmax+nxjg
                  kminc=kmin-nyjg; kmaxc=kmax+nyjg
                  lminc=lmin-nzjg; lmaxc=lmax+nzjg
                  nxc=curr_tile%nx_cells_tile
                  nyc=curr_tile%ny_cells_tile
                  nzc=curr_tile%nz_cells_tile
                  ! ----- Add guardcells in adjacent tiles
                  ! --- JX
                  ! - FACES +/- X
                  jxg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jxg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jxtile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  jxg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jxg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jxtile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JY
                  ! - FACES +/- X
                  jyg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jyg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jytile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  jyg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jyg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jytile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JZ
                  ! - FACES +/- X
                  jzg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jzg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jztile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  jzg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jzg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jztile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
              ENDIF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
              isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
                  curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  count=curr_tile%np_tile(1)
                  IF (count .GT. 0) isdeposited=.TRUE.  
              END DO
              IF (isdeposited) THEN 
              	currg=>aofgrid_tiles(ix,iy,iz)
                  curr => species_parray(1)
             		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                  kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                  lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                  nxjg=curr_tile%nxg_tile
                  nyjg=curr_tile%nyg_tile
                  nzjg=curr_tile%nzg_tile
                  jminc=jmin-nxjg; jmaxc=jmax+nxjg
                  kminc=kmin-nyjg; kmaxc=kmax+nyjg
                  lminc=lmin-nzjg; lmaxc=lmax+nzjg
                  nxc=curr_tile%nx_cells_tile
                  nyc=curr_tile%ny_cells_tile
                  nzc=curr_tile%nz_cells_tile
                  ! ----- Add guardcells in adjacent tiles
                  ! --- JX
                  ! - FACES +/- Y
                  jxg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jxg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                  currg%jxtile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                  jxg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jxg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                  currg%jxtile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JY
                  ! - FACES +/- Y
                  jyg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jyg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                  currg%jytile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                  jyg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jyg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                  currg%jytile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JZ
                  ! - FACES +/- Y
                  jzg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jzg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                  currg%jztile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                  jzg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jzg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                  currg%jztile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
              END IF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
              isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
                  curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  count=curr_tile%np_tile(1)
                  IF (count .GT. 0) isdeposited=.TRUE.  
              END DO
              IF (isdeposited) THEN 
              	currg=>aofgrid_tiles(ix,iy,iz)
                  curr => species_parray(1)
             		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                  kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                  lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                  nxjg=curr_tile%nxg_tile
                  nyjg=curr_tile%nyg_tile
                  nzjg=curr_tile%nzg_tile
                  jminc=jmin-nxjg; jmaxc=jmax+nxjg
                  kminc=kmin-nyjg; kmaxc=kmax+nyjg
                  lminc=lmin-nzjg; lmaxc=lmax+nzjg
                  nxc=curr_tile%nx_cells_tile
                  nyc=curr_tile%ny_cells_tile
                  nzc=curr_tile%nz_cells_tile
                  ! ----- Add guardcells in adjacent tiles
                  ! --- JX
                  ! - FACES +/- Z
                  jxg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jxg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                  currg%jxtile(0:nxc, 0:nyc,-nzjg:-1)
                  jxg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jxg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                  currg%jxtile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
                  ! --- JY
                  ! - FACES +/- Z
                  jyg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jyg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                  currg%jytile(0:nxc, 0:nyc,-nzjg:-1)
                  jyg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jyg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                  currg%jytile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
                  ! --- JZ
                  ! - FACES +/- Z
                  jzg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jzg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                  currg%jztile(0:nxc, 0:nyc,-nzjg:-1)
                  jzg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jzg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                  currg%jztile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
              END IF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
   
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp


!===============================================================================
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2(&
func_order,jxg,jyg,jzg,nxx,nyy,nzz,nxjguard,nyjguard,nzjguard, &
	noxx,noyy,nozz,dxx,dyy,dzz,dtt)
! Deposit current in each tile with the classical method
! OpenMP version. Avoids conflict while reducing tile currents in the global 
! current array. 
! In this second version, the transient current arrays are reduced 
! after the current deposition for all species and not for each species.
!
! Inputs: 
! - func_order: represent the subroutine to be used for current deposition depending on the selected order
! - jxg,jyg,jzg: current arrays
! - nxx,nyy,nzz: cell number in each direction
! - nxjguard,nyjguard,nzjguard: guard cells
! - noxx,noyy,nozz: orders for current deposition
! - dxx,dyy,dzz,dtt: space and time steps
!===============================================================================	
  USE particles
  USE constants
  USE tiling
  USE omp_lib
  USE timing
  USE time_stat
  IMPLICIT NONE

  ! _______________________________________________________________________
  ! Interfaces for func_order
  INTERFACE
      SUBROUTINE func_order(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz) !#do not parse
      	USE constants
      	IMPLICIT NONE
        INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz,ncells
        INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
        REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
        REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
        REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
        INTEGER(idp) :: ncx, ncy, ncz
      END SUBROUTINE func_order
  END INTERFACE

  ! _______________________________________________________________________  
  ! Parameters
  INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
  REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: jxcells,jycells,jzcells
  INTEGER(idp) :: ispecies, ix, iy, iz, np, ncells
  INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  INTEGER(idp) :: ncz,ncy,ncx
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER :: curr_tile
  TYPE(grid_tile), POINTER :: currg
  INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(idp) :: isdeposited=.FALSE.

  ! _______________________________________________________________________
  !$OMP PARALLEL DEFAULT(NONE)                                                              &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguard,nyjguard,              &
  !$OMP nzjguard,dxx,dyy,dzz,dtt,jxg,jyg,jzg,noxx,noyy,nozz,aofgrid_tiles)                  &
  !$OMP PRIVATE(ix,iy,iz,ispecies,ncells,curr,currg, curr_tile,np,jmin,jmax,kmin,kmax,lmin,  &
  !$OMP lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,nxc,nyc,nzc, nxjg, nyjg, nzjg, isdeposited,  &
  !$OMP jxcells,jycells,jzcells,ncx,ncy,ncz)
  
  !! Current deposition
  ! Loop on the tiles
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
          	  curr => species_parray(1)
              curr_tile=>curr%array_of_tiles(ix,iy,iz)
              nxjg=curr_tile%nxg_tile
              nyjg=curr_tile%nyg_tile
              nzjg=curr_tile%nzg_tile
              jmin=curr_tile%nx_tile_min
          	  jmax=curr_tile%nx_tile_max
              kmin=curr_tile%ny_tile_min
              kmax=curr_tile%ny_tile_max
              lmin=curr_tile%nz_tile_min
              lmax=curr_tile%nz_tile_max
              nxc=curr_tile%nx_cells_tile
              nyc=curr_tile%ny_cells_tile
              nzc=curr_tile%nz_cells_tile         
			        currg=>aofgrid_tiles(ix,iy,iz)
              currg%jxtile=0.
              currg%jytile=0.
              currg%jztile=0.!jzg(jmin:jmax,kmin:kmax,lmin:lmax)
              isdeposited=.FALSE.
              
              ncx = nxc+noxx+2
              ncy = nyc+noyy+2
              ncz = nzc+nozz+2
              
              IF ((noxx.eq.3).and.(noyy.eq.3).and.(nozz.eq.3))  THEN
                ncx=nxc+5; ncy=nyc+4; ncz=nzc+3
              ENDIF
              
              ncells = ncx*ncy*ncz
              ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
              jxcells=0.0_num
              jycells=0.0_num
              jzcells=0.0_num
              
              DO ispecies=1, nspecies ! LOOP ON SPECIES
             	    curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  np=curr_tile%np_tile(1)
                  IF (np .EQ. 0) THEN 
                  	CYCLE
                  ELSE 
                  	isdeposited=.TRUE.
                  ENDIF 
                  
                  ! Depose current in jtile
                  CALL func_order(jxcells,jycells,jzcells,np,ncells, &
                  curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						           &
                  curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,  			   &
                  curr_tile%part_gaminv,  &
                  curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,     &
                  curr_tile%z_grid_tile_min,dtt,dxx,dyy,dzz,nxc,nyc,nzc,                                     &
                  nxjg,nyjg,nzjg,ncx,ncy,ncz) 
                                                      
                  !print*,sum(jxcells), sum(jycells), sum(jzcells)                             
                                           
              END DO! END LOOP ON SPECIES
              
              CALL current_reduction(currg%jxtile,currg%jytile,currg%jztile, &
              jxcells,jycells,jzcells,ncells,nxc,nyc,nzc,nxjg,nyjg,nzjg,ncx,ncy,ncz)
              
              DEALLOCATE(jxcells,jzcells,jycells)
              
              IF (isdeposited) THEN
              	jxg(jmin:jmax,kmin:kmax,lmin:lmax)=jxg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jxtile(0:nxc,0:nyc,0:nzc)
              	jyg(jmin:jmax,kmin:kmax,lmin:lmax)=jyg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jytile(0:nxc,0:nyc,0:nzc)
              	jzg(jmin:jmax,kmin:kmax,lmin:lmax)=jzg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jztile(0:nxc,0:nyc,0:nzc)
              ENDIF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
  !+/- X
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
          	isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
                  curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  np=curr_tile%np_tile(1)
                  IF (np .GT. 0) isdeposited=.TRUE.  
              END DO
              IF (isdeposited) THEN 
              	currg=>aofgrid_tiles(ix,iy,iz)
                  curr => species_parray(1)
             		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                  kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                  lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                  nxjg=curr_tile%nxg_tile
                  nyjg=curr_tile%nyg_tile
                  nzjg=curr_tile%nzg_tile
                  jminc=jmin-nxjg; jmaxc=jmax+nxjg
                  kminc=kmin-nyjg; kmaxc=kmax+nyjg
                  lminc=lmin-nzjg; lmaxc=lmax+nzjg
                  nxc=curr_tile%nx_cells_tile
                  nyc=curr_tile%ny_cells_tile
                  nzc=curr_tile%nz_cells_tile
                  ! ----- Add guardcells in adjacent tiles
                  ! --- JX
                  ! - FACES +/- X
                  jxg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jxg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jxtile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  jxg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jxg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jxtile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JY
                  ! - FACES +/- X
                  jyg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jyg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jytile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  jyg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jyg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jytile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JZ
                  ! - FACES +/- X
                  jzg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jzg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jztile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                  jzg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jzg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                  currg%jztile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
              ENDIF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !+/- Y
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
              isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
                  curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  np=curr_tile%np_tile(1)
                  IF (np .GT. 0) isdeposited=.TRUE.  
              END DO
              IF (isdeposited) THEN 
              	currg=>aofgrid_tiles(ix,iy,iz)
                  curr => species_parray(1)
             		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                  kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                  lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                  nxjg=curr_tile%nxg_tile
                  nyjg=curr_tile%nyg_tile
                  nzjg=curr_tile%nzg_tile
                  jminc=jmin-nxjg; jmaxc=jmax+nxjg
                  kminc=kmin-nyjg; kmaxc=kmax+nyjg
                  lminc=lmin-nzjg; lmaxc=lmax+nzjg
                  nxc=curr_tile%nx_cells_tile
                  nyc=curr_tile%ny_cells_tile
                  nzc=curr_tile%nz_cells_tile
                  ! ----- Add guardcells in adjacent tiles
                  ! --- JX
                  ! - FACES +/- Y
                  jxg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jxg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                  currg%jxtile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                  jxg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jxg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                  currg%jxtile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JY
                  ! - FACES +/- Y
                  jyg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jyg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                  currg%jytile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                  jyg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jyg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                  currg%jytile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
                  ! --- JZ
                  ! - FACES +/- Y
                  jzg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jzg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                  currg%jztile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                  jzg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jzg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                  currg%jztile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
              END IF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  ! +/-Z
  !$OMP DO COLLAPSE(3) SCHEDULE(runtime)
  DO iz=1,ntilez
      DO iy=1,ntiley
          DO ix=1,ntilex
              isdeposited=.FALSE.
              DO ispecies=1, nspecies ! LOOP ON SPECIES
                  curr => species_parray(ispecies)
                  curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  np=curr_tile%np_tile(1)
                  IF (np .GT. 0) isdeposited=.TRUE.  
              END DO
              IF (isdeposited) THEN 
              	currg=>aofgrid_tiles(ix,iy,iz)
                  curr => species_parray(1)
             		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                  jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                  kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                  lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                  nxjg=curr_tile%nxg_tile
                  nyjg=curr_tile%nyg_tile
                  nzjg=curr_tile%nzg_tile
                  jminc=jmin-nxjg; jmaxc=jmax+nxjg
                  kminc=kmin-nyjg; kmaxc=kmax+nyjg
                  lminc=lmin-nzjg; lmaxc=lmax+nzjg
                  nxc=curr_tile%nx_cells_tile
                  nyc=curr_tile%ny_cells_tile
                  nzc=curr_tile%nz_cells_tile
                  ! ----- Add guardcells in adjacent tiles
                  ! --- JX
                  ! - FACES +/- Z
                  jxg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jxg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                  currg%jxtile(0:nxc, 0:nyc,-nzjg:-1)
                  jxg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jxg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                  currg%jxtile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
                  ! --- JY
                  ! - FACES +/- Z
                  jyg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jyg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                  currg%jytile(0:nxc, 0:nyc,-nzjg:-1)
                  jyg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jyg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                  currg%jytile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
                  ! --- JZ
                  ! - FACES +/- Z
                  jzg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jzg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                  currg%jztile(0:nxc, 0:nyc,-nzjg:-1)
                  jzg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jzg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                  currg%jztile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
              END IF
          END DO
      END DO
  END DO!END LOOP ON TILES
  !$OMP END DO
  !$OMP END PARALLEL
   
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2


!===============================================================================
! Deposit current in each tile with Esirkepov method
! OpenMP version. Avoids conflict while reducing tile currents in the global 
! current array. 
!===============================================================================
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp(func_order,jxg,jyg,jzg,nxx,nyy,nzz,nxjguard,nyjguard,nzjguard, &
	noxx,noyy,nozz,dxx,dyy,dzz,dtt)
USE particles
USE constants
USE tiling
USE omp_lib
USE timing
USE time_stat
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
INTEGER(idp) :: ispecies, ix, iy, iz, count
INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
TYPE(grid_tile), POINTER :: currg
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
LOGICAL(idp) :: isdeposited=.FALSE.

! For the func_order input function
INTERFACE
  SUBROUTINE func_order(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, & !#do not parse
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, & !#do not parse
           nox,noy,noz,l_particles_weight,l4symtry) !#do not parse

    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    LOGICAL(idp) :: l_particles_weight,l4symtry
  END SUBROUTINE

END INTERFACE

tdeb=MPI_WTIME()
!$OMP PARALLEL DEFAULT(NONE)                                                              &
!$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguard,nyjguard,              &
!$OMP nzjguard,dxx,dyy,dzz,dtt,jxg,jyg,jzg,noxx,noyy,nozz,aofgrid_tiles)                  &
!$OMP PRIVATE(ix,iy,iz,ispecies,curr,currg, curr_tile,count,jmin,jmax,kmin,kmax,lmin,     &
!$OMP lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,nxc,nyc,nzc, nxjg, nyjg, nzjg, isdeposited)
!! Current deposition
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
        	curr => species_parray(1)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
            nxjg=curr_tile%nxg_tile
            nyjg=curr_tile%nyg_tile
            nzjg=curr_tile%nzg_tile
            jmin=curr_tile%nx_tile_min
        	jmax=curr_tile%nx_tile_max
            kmin=curr_tile%ny_tile_min
            kmax=curr_tile%ny_tile_max
            lmin=curr_tile%nz_tile_min
            lmax=curr_tile%nz_tile_max
            nxc=curr_tile%nx_cells_tile;
            nyc=curr_tile%ny_cells_tile
            nzc=curr_tile%nz_cells_tile         
			      currg=>aofgrid_tiles(ix,iy,iz)
            currg%jxtile=0.
            currg%jytile=0.
            currg%jztile=0.!jzg(jmin:jmax,kmin:kmax,lmin:lmax)
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
           	    curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .EQ. 0) THEN 
                	CYCLE
                ELSE 
                	isdeposited=.TRUE.
                ENDIF 
                
                ! Depose current in jtile                
                CALL func_order(currg%jxtile,currg%jytile,                              &
                currg%jztile,count,                              									       &
                curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						           &
                curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,curr_tile%part_gaminv,  			   &
                curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,     &
                curr_tile%z_grid_tile_min,dtt,dxx,dyy,dzz,nxc,nyc,nzc,                                     &
                nxjg,nyjg,nzjg,noxx,noyy,nozz,.TRUE._idp,.FALSE._idp) 
                
                !CALL pxr_depose_jxjyjz_esirkepov_n(currg%jxtile,currg%jytile,                              &
                !currg%jztile,count,                              									       &
                !curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						           &
                !curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,curr_tile%part_gaminv,  			   &
                !curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,     &
                !curr_tile%z_grid_tile_min,dtt,dxx,dyy,dzz,nxc,nyc,nzc,                                     &
                !nxjg,nyjg,nzjg,noxx,noyy,nozz,.TRUE.,.FALSE.) 
                
            END DO! END LOOP ON SPECIES
            IF (isdeposited) THEN
            	jxg(jmin:jmax,kmin:kmax,lmin:lmax)=jxg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jxtile(0:nxc,0:nyc,0:nzc)
            	jyg(jmin:jmax,kmin:kmax,lmin:lmax)=jyg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jytile(0:nxc,0:nyc,0:nzc)
            	jzg(jmin:jmax,kmin:kmax,lmin:lmax)=jzg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jztile(0:nxc,0:nyc,0:nzc)
            ENDIF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!! Adding currents from guard cells of adjacent subdomains (AVOIDS REDUCTION OPERATION)
!+/- X
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
        	isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .GT. 0) isdeposited=.TRUE.  
            END DO
            IF (isdeposited) THEN 
            	currg=>aofgrid_tiles(ix,iy,iz)
                curr => species_parray(1)
           		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxjg=curr_tile%nxg_tile
                nyjg=curr_tile%nyg_tile
                nzjg=curr_tile%nzg_tile
                jminc=jmin-nxjg; jmaxc=jmax+nxjg
                kminc=kmin-nyjg; kmaxc=kmax+nyjg
                lminc=lmin-nzjg; lmaxc=lmax+nzjg
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! --- JX
                ! - FACES +/- X
                jxg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jxg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                currg%jxtile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                jxg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jxg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                currg%jxtile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                ! --- JY
                ! - FACES +/- X
                jyg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jyg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                currg%jytile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                jyg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jyg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                currg%jytile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                ! --- JZ
                ! - FACES +/- X
                jzg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc) = jzg(jminc:jmin-1,kminc:kmaxc,lminc:lmaxc)+  &
                currg%jztile(-nxjg:-1,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
                jzg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc) = jzg(jmax+1:jmaxc,kminc:kmaxc,lminc:lmaxc)+  &
                currg%jztile(nxc+1:nxc+nxjg,-nyjg:nyc+nyjg,-nzjg:nzc+nzjg)
            ENDIF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!+/- Y
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .GT. 0) isdeposited=.TRUE.  
            END DO
            IF (isdeposited) THEN 
            	currg=>aofgrid_tiles(ix,iy,iz)
                curr => species_parray(1)
           		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxjg=curr_tile%nxg_tile
                nyjg=curr_tile%nyg_tile
                nzjg=curr_tile%nzg_tile
                jminc=jmin-nxjg; jmaxc=jmax+nxjg
                kminc=kmin-nyjg; kmaxc=kmax+nyjg
                lminc=lmin-nzjg; lmaxc=lmax+nzjg
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! --- JX
                ! - FACES +/- Y
                jxg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jxg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                currg%jxtile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                jxg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jxg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                currg%jxtile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
                ! --- JY
                ! - FACES +/- Y
                jyg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jyg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                currg%jytile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                jyg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jyg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                currg%jytile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
                ! --- JZ
                ! - FACES +/- Y
                jzg(jmin:jmax,kminc:kmin-1,lminc:lmaxc) = jzg(jmin:jmax,kminc:kmin-1,lminc:lmaxc)+  &
                currg%jztile(0:nxc,-nyjg:-1,-nzjg:nzc+nzjg)
                jzg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc) = jzg(jmin:jmax,kmax+1:kmaxc,lminc:lmaxc)+  &
                currg%jztile(0:nxc,nyc+1:nyc+nyjg,-nzjg:nzc+nzjg)
            END IF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
! +/-Z
!$OMP DO COLLAPSE(3) SCHEDULE(runtime)
DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .GT. 0) isdeposited=.TRUE.  
            END DO
            IF (isdeposited) THEN 
            	currg=>aofgrid_tiles(ix,iy,iz)
                curr => species_parray(1)
           		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                jmin=curr_tile%nx_tile_min; jmax=curr_tile%nx_tile_max
                kmin=curr_tile%ny_tile_min; kmax=curr_tile%ny_tile_max
                lmin=curr_tile%nz_tile_min; lmax=curr_tile%nz_tile_max
                nxjg=curr_tile%nxg_tile
                nyjg=curr_tile%nyg_tile
                nzjg=curr_tile%nzg_tile
                jminc=jmin-nxjg; jmaxc=jmax+nxjg
                kminc=kmin-nyjg; kmaxc=kmax+nyjg
                lminc=lmin-nzjg; lmaxc=lmax+nzjg
                nxc=curr_tile%nx_cells_tile
                nyc=curr_tile%ny_cells_tile
                nzc=curr_tile%nz_cells_tile
                ! ----- Add guardcells in adjacent tiles
                ! --- JX
                ! - FACES +/- Z
                jxg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jxg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                currg%jxtile(0:nxc, 0:nyc,-nzjg:-1)
                jxg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jxg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                currg%jxtile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
                ! --- JY
                ! - FACES +/- Z
                jyg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jyg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                currg%jytile(0:nxc, 0:nyc,-nzjg:-1)
                jyg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jyg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                currg%jytile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
                ! --- JZ
                ! - FACES +/- Z
                jzg(jmin:jmax,kmin:kmax,lminc:lmin-1) = jzg(jmin:jmax,kmin:kmax,lminc:lmin-1)+  &
                currg%jztile(0:nxc, 0:nyc,-nzjg:-1)
                jzg(jmin:jmax,kmin:kmax,lmax+1:lmaxc) = jzg(jmin:jmax,kmin:kmax,lmax+1:lmaxc)+  &
                currg%jztile(0:nxc, 0:nyc,nzc+1:nzc+nzjg)
            END IF
        END DO
    END DO
END DO!END LOOP ON TILES
!$OMP END DO
!$OMP END PARALLEL
tend=MPI_WTIME()
dep_curr_time=dep_curr_time+(tend-tdeb)  
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp

SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq(func_order,jxg,jyg,jzg,nxx,nyy,nzz,nxjguard,nyjguard,nzjguard, &
	noxx,noyy,nozz,dxx,dyy,dzz,dtt)
USE particles
USE constants
USE tiling
USE timing
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), POINTER, DIMENSION(:,:,:) :: jxp, jyp, jzp
INTEGER(idp) :: ispecies, ix, iy, iz, count
INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
TYPE(grid_tile), POINTER :: currg
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
LOGICAL(idp) :: isdeposited=.FALSE.

INTERFACE
  SUBROUTINE func_order(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, & !#do not parse
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard) !#do not parse
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
      REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp,w,gaminv
      REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
  END SUBROUTINE
END INTERFACE

DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
        	curr => species_parray(1)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
            nxjg=curr_tile%nxg_tile
            nyjg=curr_tile%nyg_tile
            nzjg=curr_tile%nzg_tile
            jmin=curr_tile%nx_tile_min-nxjg
        	jmax=curr_tile%nx_tile_max+nxjg
            kmin=curr_tile%ny_tile_min-nyjg
            kmax=curr_tile%ny_tile_max+nyjg
            lmin=curr_tile%nz_tile_min-nzjg
            lmax=curr_tile%nz_tile_max+nzjg
            nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
            nzc=curr_tile%nz_cells_tile         
			currg=>aofgrid_tiles(ix,iy,iz)
            currg%jxtile=0.
            currg%jytile=0.
            currg%jztile=0.!jzg(jmin:jmax,kmin:kmax,lmin:lmax)
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
           	    curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .EQ. 0) THEN 
                	CYCLE
                ELSE 
                	isdeposited=.TRUE.
                ENDIF 
                ! Depose current in jtile
                CALL func_order(currg%jxtile,currg%jytile,currg%jztile,count,   &
                curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						                      &
                curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,curr_tile%part_gaminv,  			      &
                curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,        &
                curr_tile%z_grid_tile_min,dtt,dxx,dyy,dzz,nxc,nyc,nzc,                                        &
                nxjg,nyjg,nzjg) 
            END DO! END LOOP ON SPECIES
            IF (isdeposited) THEN
            	jxg(jmin:jmax,kmin:kmax,lmin:lmax)=jxg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jxtile
            	jyg(jmin:jmax,kmin:kmax,lmin:lmax)=jyg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jytile
            	jzg(jmin:jmax,kmin:kmax,lmin:lmax)=jzg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jztile
            ENDIF
        END DO
    END DO
END DO!END LOOP ON TILES


END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq


!===============================================================================
! Deposit current in each tile
! Sequential version 
!===============================================================================
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq(jxg,jyg,jzg,nxx,nyy,nzz,nxjguard,nyjguard,nzjguard, &
	noxx,noyy,nozz,dxx,dyy,dzz,dtt)
USE particles
USE constants
USE tiling
USE timing
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
INTEGER(idp), INTENT(IN) :: noxx,noyy,nozz
REAL(num), INTENT(IN) :: dxx,dyy,dzz, dtt
REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
REAL(num), POINTER, DIMENSION(:,:,:) :: jxp, jyp, jzp
INTEGER(idp) :: ispecies, ix, iy, iz, count
INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER(idp) :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
TYPE(particle_species), POINTER :: curr
TYPE(particle_tile), POINTER :: curr_tile
TYPE(grid_tile), POINTER :: currg
REAL(num) :: tdeb, tend
INTEGER(idp) :: nxc, nyc, nzc, nxjg, nyjg, nzjg
LOGICAL(idp) :: isdeposited=.FALSE.

DO iz=1,ntilez
    DO iy=1,ntiley
        DO ix=1,ntilex
        	curr => species_parray(1)
            curr_tile=>curr%array_of_tiles(ix,iy,iz)
            nxjg=curr_tile%nxg_tile
            nyjg=curr_tile%nyg_tile
            nzjg=curr_tile%nzg_tile
            jmin=curr_tile%nx_tile_min-nxjg
        	jmax=curr_tile%nx_tile_max+nxjg
            kmin=curr_tile%ny_tile_min-nyjg
            kmax=curr_tile%ny_tile_max+nyjg
            lmin=curr_tile%nz_tile_min-nzjg
            lmax=curr_tile%nz_tile_max+nzjg
            nxc=curr_tile%nx_cells_tile; nyc=curr_tile%ny_cells_tile
            nzc=curr_tile%nz_cells_tile         
			currg=>aofgrid_tiles(ix,iy,iz)
            currg%jxtile=0.
            currg%jytile=0.
            currg%jztile=0.!jzg(jmin:jmax,kmin:kmax,lmin:lmax)
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
           	    curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                IF (count .EQ. 0) THEN 
                	CYCLE
                ELSE 
                	isdeposited=.TRUE.
                ENDIF 
                ! Depose current in jtile
                CALL pxr_depose_jxjyjz_esirkepov_n(currg%jxtile,currg%jytile,                                 &
                currg%jztile,count,                              									                        &
                curr_tile%part_x,curr_tile%part_y,curr_tile%part_z,     						                      &
                curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,curr_tile%part_gaminv,  			      &
                curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,curr_tile%y_grid_tile_min,        &
                curr_tile%z_grid_tile_min,dtt,dxx,dyy,dzz,nxc,nyc,nzc,                                        &
                nxjg,nyjg,nzjg,noxx,noyy,nozz,.TRUE._idp,.FALSE._idp) 
            END DO! END LOOP ON SPECIES
            IF (isdeposited) THEN
            	jxg(jmin:jmax,kmin:kmax,lmin:lmax)=jxg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jxtile
            	jyg(jmin:jmax,kmin:kmax,lmin:lmax)=jyg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jytile
            	jzg(jmin:jmax,kmin:kmax,lmin:lmax)=jzg(jmin:jmax,kmin:kmax,lmin:lmax)+currg%jztile
            ENDIF
        END DO
    END DO
END DO!END LOOP ON TILES


END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq



!!! ---------  CURRENT DEPOSITION ROUTINES 

!!! --- Order 1 3D scalar current deposition routine (rho*v)
!!! This version does not vectorize on SIMD architectures
SUBROUTINE depose_jxjyjz_scalar_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,vx,vy,vz,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) :: wq, wqx, wqy, wqz, clightsq
    REAL(num), DIMENSION(2) :: sx(0:1), sy(0:1), sz(0:1), sx0(0:1), sy0(0:1), sz0(0:1)
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip
    dxi = 1.0_num/dx
    dyi = 1.0_num/dy
    dzi = 1.0_num/dz
    invvol = dxi*dyi*dzi
    dts2dx = 0.5_num*dt*dxi
    dts2dy = 0.5_num*dt*dyi
    dts2dz = 0.5_num*dt*dzi
    clightsq = 1.0_num/clight**2
    sx=0.0_num;sy=0.0_num;sz=0.0_num;
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;

    ! LOOP ON PARTICLES
    DO ip=1,np
        ! --- computes position in  grid units at (n+1)
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        ! Computes velocity
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)

        ! --- computes particles weights
        wq=q*w(ip)
        wqx=wq*invvol*vx
        wqy=wq*invvol*vy
        wqz=wq*invvol*vz

        ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
        xmid=x-dts2dx*vx
        ymid=y-dts2dy*vy
        zmid=z-dts2dz*vz

        ! --- finds node of cell containing particles for current positions
        j=floor(xmid)
        k=floor(ymid)
        l=floor(zmid)
        j0=floor(xmid-0.5_num)
        k0=floor(ymid-0.5_num)
        l0=floor(zmid-0.5_num)
        ! --- computes set of coefficients for node centered quantities
        xint = xmid-j
        yint = ymid-k
        zint = zmid-l
        sx( 0) = 1.0_num-xint
        sx( 1) = xint
        sy( 0) = 1.0_num-yint
        sy( 1) = yint
        sz( 0) = 1.0_num-zint
        sz( 1) = zint

        ! --- computes set of coefficients for staggered quantities
        xint = xmid-j0-0.5_num
        yint = ymid-k0-0.5_num
        zint = zmid-l0-0.5_num
        sx0( 0) = 1.0_num-xint
        sx0( 1) = xint
        sy0( 0) = 1.0_num-yint
        sy0( 1) = yint
        sz0( 0) = 1.0_num-zint
        sz0( 1) = zint

        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! - JX
        jx(j0  ,k  ,l  )    = jx(j0  ,k  ,l  )  +   sx0(0)*sy(0)*sz(0)*wqx
        jx(j0+1,k  ,l  )    = jx(j0+1,k  ,l  )  +   sx0(1)*sy(0)*sz(0)*wqx
        jx(j0  ,k+1,l  )    = jx(j0  ,k+1,l  )  +   sx0(0)*sy(1)*sz(0)*wqx
        jx(j0+1,k+1,l  )    = jx(j0+1,k+1,l  )  +   sx0(1)*sy(1)*sz(0)*wqx
        jx(j0  ,k  ,l+1)    = jx(j0  ,k  ,l+1)  +   sx0(0)*sy(0)*sz(1)*wqx
        jx(j0+1,k  ,l+1)    = jx(j0+1,k  ,l+1)  +   sx0(1)*sy(0)*sz(1)*wqx
        jx(j0  ,k+1,l+1)    = jx(j0  ,k+1,l+1)  +   sx0(0)*sy(1)*sz(1)*wqx
        jx(j0+1,k+1,l+1)    = jx(j0+1,k+1,l+1)  +   sx0(1)*sy(1)*sz(1)*wqx

        ! - JY
        jy(j  ,k0  ,l  )    = jy(j  ,k0  ,l  )  +   sx(0)*sy0(0)*sz(0)*wqy
        jy(j+1,k0  ,l  )    = jy(j+1,k0  ,l  )  +   sx(1)*sy0(0)*sz(0)*wqy
        jy(j  ,k0+1,l  )    = jy(j  ,k0+1,l  )  +   sx(0)*sy0(1)*sz(0)*wqy
        jy(j+1,k0+1,l  )    = jy(j+1,k0+1,l  )  +   sx(1)*sy0(1)*sz(0)*wqy
        jy(j  ,k0  ,l+1)    = jy(j  ,k0  ,l+1)  +   sx(0)*sy0(0)*sz(1)*wqy
        jy(j+1,k0  ,l+1)    = jy(j+1,k0  ,l+1)  +   sx(1)*sy0(0)*sz(1)*wqy
        jy(j  ,k0+1,l+1)    = jy(j  ,k0+1,l+1)  +   sx(0)*sy0(1)*sz(1)*wqy
        jy(j+1,k0+1,l+1)    = jy(j+1,k0+1,l+1)  +   sx(1)*sy0(1)*sz(1)*wqy

        ! - JZ
        jz(j  ,k  ,l0  )    = jz(j  ,k  ,l0  )  +   sx(0)*sy(0)*sz0(0)*wqz
        jz(j+1,k  ,l0  )    = jz(j+1,k  ,l0  )  +   sx(1)*sy(0)*sz0(0)*wqz
        jz(j  ,k+1,l0  )    = jz(j  ,k+1,l0  )  +   sx(0)*sy(1)*sz0(0)*wqz
        jz(j+1,k+1,l0  )    = jz(j+1,k+1,l0  )  +   sx(1)*sy(1)*sz0(0)*wqz
        jz(j  ,k  ,l0+1)    = jz(j  ,k  ,l0+1)  +   sx(0)*sy(0)*sz0(1)*wqz
        jz(j+1,k  ,l0+1)    = jz(j+1,k  ,l0+1)  +   sx(1)*sy(0)*sz0(1)*wqz
        jz(j  ,k+1,l0+1)    = jz(j  ,k+1,l0+1)  +   sx(0)*sy(1)*sz0(1)*wqz
        jz(j+1,k+1,l0+1)    = jz(j+1,k+1,l0+1)  +   sx(1)*sy(1)*sz0(1)*wqz
    END DO
    RETURN
END SUBROUTINE depose_jxjyjz_scalar_1_1_1



! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHVv2_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
!!! --- Order 1 3D vector current deposition routine (rho*v)
!!! This versions have good performances on SIMD architectures
!!! Providing that OpenMP 4.0 is available (Directive SIMD)
! ________________________________________________________________________________________
    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) :: clightsq

    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic
    INTEGER(idp) :: nnx, nnxy, n,nn,nv
    INTEGER(idp) :: moff(1:8)
    REAL(num):: mx(1:8),my(1:8),mz(1:8), sgn(1:8)

    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL
    REAL(num), DIMENSION(LVEC) :: sx, sy, sz, sx0, sy0, sz0,wqx,wqy,wqz
    REAL(num) :: wwx,wwy,wwz, wq,vx,vy,vz, wx,wx0, wy,wy0, wz,wz0
    INTEGER(idp) :: orig, jorig, korig, lorig, igrid
    INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy

    dxi = 1.0_num/dx
    dyi = 1.0_num/dy
    dzi = 1.0_num/dz
    invvol = dxi*dyi*dzi
    dts2dx = 0.5_num*dt*dxi
    dts2dy = 0.5_num*dt*dyi
    dts2dz = 0.5_num*dt*dzi
    sx=0.0_num;sy=0.0_num;sz=0.0_num
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
    ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard;
    ncx=nx+3;ncy=ny+3;ncz=nz+3
    NCELLS=ncx*ncy*ncz
    ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
    jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num;
    nnx = ngridx
    nnxy = nnx*ngridy
    
    ! Relative position of the nodes in respect to the particle computed node
    moff = (/0_idp,1_idp,nnx,nnx+1_idp,nnxy,nnxy+1_idp,nnxy+nnx,nnxy+nnx+1_idp/)
    
    mx=(/1_num,0_num,1_num,0_num,1_num,0_num,1_num,0_num/)
    my=(/1_num,1_num,0_num,0_num,1_num,1_num,0_num,0_num/)
    mz=(/1_num,1_num,1_num,1_num,0_num,0_num,0_num,0_num/)
    sgn=(/-1_num,1_num,1_num,-1_num,1_num,-1_num,-1_num,1_num/)
    jorig=-2; korig=-2;lorig=-2
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy
    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
        !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64
        !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
        !DIR$ ASSUME_ALIGNED w:64, gaminv:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx = uxp(nn)*gaminv(nn)
            vy = uyp(nn)*gaminv(nn)
            vz = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol
            wqx(n)=wq*vx
            wqy(n)=wq*vy
            wqz(n)=wq*vz

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx
            ymid=y-dts2dy*vy
            zmid=z-dts2dz*vz

            ! --- finds node of cell containing particles for current positions
            j=floor(xmid)
            k=floor(ymid)
            l=floor(zmid)
            j0=floor(xmid-0.5_num)
            k0=floor(ymid-0.5_num)
            l0=floor(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

            ! --- computes set of coefficients for node centered quantities
            sx(n) = xmid-j
            sy(n) = ymid-k
            sz(n) = zmid-l

            ! --- computes set of coefficients for staggered quantities
            sx0(n) = xmid-j0-0.5_num
            sy0(n) = ymid-k0-0.5_num
            sz0(n) = zmid-l0-0.5_num
        END DO
        !$OMP END SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !DIR$ ASSUME_ALIGNED mx:64, my:64, mz:64,sgn:64
            !$OMP SIMD
            DO nv=1,8
                wx=-mx(nv)+sx(n)
                wx0=-mx(nv)+sx0(n)
                wy=-my(nv)+sy(n)
                wy0=-my(nv)+sy0(n)
                wz=-mz(nv)+sz(n)
                wz0=-mz(nv)+sz0(n)
                wwx=wx0*wy*wz*wqx(n)*sgn(nv)
                wwy=wx*wy0*wz*wqy(n)*sgn(nv)
                wwz=wx*wy*wz0*wqz(n)*sgn(nv)
                ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
                ! - JX
                jxcells(nv,ICELL(n,1))=jxcells(nv,ICELL(n,1))+wwx
                ! - JY
                jycells(nv,ICELL(n,2))=jycells(nv,ICELL(n,2))+wwy
                ! - JZ
                jzcells(nv,ICELL(n,3))=jzcells(nv,ICELL(n,3))+wwz
            END DO
            !$OMP END SIMD
        END DO
    END DO
    ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
    DO iz=1, ncz
        DO iy=1,ncy
            !$OMP SIMD
            DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                ! jx
                jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1,ic)
                jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2,ic)
                jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3,ic)
                jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4,ic)
                jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5,ic)
                jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6,ic)
                jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7,ic)
                jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8,ic)
                ! jy
                jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1,ic)
                jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2,ic)
                jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3,ic)
                jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4,ic)
                jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5,ic)
                jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6,ic)
                jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7,ic)
                jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8,ic)
                ! jz
                jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1,ic)
                jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2,ic)
                jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3,ic)
                jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4,ic)
                jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5,ic)
                jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6,ic)
                jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7,ic)
                jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8,ic)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    DEALLOCATE(jxcells,jycells,jzcells)
    RETURN
END SUBROUTINE depose_jxjyjz_vecHVv2_1_1_1

! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz)
! Order 1 3D vector current deposition routine (rho*v)
! This versions have good performances on SIMD architectures
! Providing that OpenMP 4.0 is available (Directive SIMD)
! This subroutine is similar to depose_jxjyjz_vecHVv2_1_1_1 
! without the reduction process at the end
!
! Inputs:
! - jxcells, jycells, jzcells: transient current arrays
! - np: particle number
! - ncells: number of cells in the tile
! - xp,yp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse Lorentz factor arrays
! - w: particle wight arrays
! - q: charge
! - xmin,ymin,zmin: tile minimum positions
! - dt,dx,dy,dz: time and space steps
! - nx,ny,nz: tile cell numbers in each direction
! - nxguard,nyguard,nzguard: guard cells
! - ncx,ncy,ncz: tile cell extended number (depends on the order)
!
! Ouputs:
! - jxcells, jycells, jzcells updated
! ________________________________________________________________________________________
    USE constants
    USE precomputed
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz
    INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
    REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
    REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: xint,yint,zint
    REAL(num) :: x,y,z,xmid,ymid,zmid

    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic
    INTEGER(idp) :: n,nn,nv
    REAL(num):: mx(1:8),my(1:8),mz(1:8), sgn(1:8)
    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL
    REAL(num), DIMENSION(LVEC) :: sx, sy, sz, sx0, sy0, sz0,wqx,wqy,wqz
    REAL(num) :: wwx,wwy,wwz, wq,vx,vy,vz, wx,wx0, wy,wy0, wz,wz0
    INTEGER(idp) :: jorig, korig, lorig
    INTEGER(idp) :: ncx, ncy, ncxy, ncz
    INTEGER(idp) :: ix,iy,iz
 
    ! Parameters
    ncxy=ncx*ncy
    
    sx=0.0_num;sy=0.0_num;sz=0.0_num
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num

    mx=(/1_num,0_num,1_num,0_num,1_num,0_num,1_num,0_num/)
    my=(/1_num,1_num,0_num,0_num,1_num,1_num,0_num,0_num/)
    mz=(/1_num,1_num,1_num,1_num,0_num,0_num,0_num,0_num/)
    sgn=(/-1_num,1_num,1_num,-1_num,1_num,-1_num,-1_num,1_num/)
    
    jorig=-2
    korig=-2
    lorig=-2

    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
        !DIR$ ASSUME_ALIGNED sx:64,sy:64,sz:64
        !DIR$ ASSUME_ALIGNED sx0:64,sy0:64,sz0:64
        !DIR$ ASSUME_ALIGNED w:64, gaminv:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx = uxp(nn)*gaminv(nn)
            vy = uyp(nn)*gaminv(nn)
            vz = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol
            wqx(n)=wq*vx
            wqy(n)=wq*vy
            wqz(n)=wq*vz

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx
            ymid=y-dts2dy*vy
            zmid=z-dts2dz*vz

            ! --- finds node of cell containing particles for current positions
            j=floor(xmid)
            k=floor(ymid)
            l=floor(zmid)
            j0=floor(xmid-0.5_num)
            k0=floor(ymid-0.5_num)
            l0=floor(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

            ! --- computes set of coefficients for node centered quantities
            sx(n) = xmid-j
            sy(n) = ymid-k
            sz(n) = zmid-l

            ! --- computes set of coefficients for staggered quantities
            sx0(n) = xmid-j0-0.5_num
            sy0(n) = ymid-k0-0.5_num
            sz0(n) = zmid-l0-0.5_num
        END DO
        !$OMP END SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !DIR$ ASSUME_ALIGNED mx:64, my:64, mz:64,sgn:64
            !$OMP SIMD
            DO nv=1,8
                wx=-mx(nv)+sx(n)
                wx0=-mx(nv)+sx0(n)
                wy=-my(nv)+sy(n)
                wy0=-my(nv)+sy0(n)
                wz=-mz(nv)+sz(n)
                wz0=-mz(nv)+sz0(n)
                wwx=wx0*wy*wz*wqx(n)*sgn(nv)
                wwy=wx*wy0*wz*wqy(n)*sgn(nv)
                wwz=wx*wy*wz0*wqz(n)*sgn(nv)
                ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
                ! - JX
                jxcells(nv,ICELL(n,1))=jxcells(nv,ICELL(n,1))+wwx
                ! - JY
                jycells(nv,ICELL(n,2))=jycells(nv,ICELL(n,2))+wwy
                ! - JZ
                jzcells(nv,ICELL(n,3))=jzcells(nv,ICELL(n,3))+wwz
            END DO
            !$OMP END SIMD
        END DO
    END DO

    RETURN
END SUBROUTINE depose_jxjyjz_vecHV_vnr_1_1_1


!!! --- Order 2 3D scalar current deposition routine (jx*v)
!!! This version does not vectorize on SIMD architectures
SUBROUTINE depose_jxjyjz_scalar_2_2_2(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,vx,vy,vz,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) :: wq, wqx, wqy, wqz, clightsq
    REAL(num), DIMENSION(3) :: sx(-1:1), sy(-1:1), sz(-1:1), sx0(-1:1), sy0(-1:1), sz0(-1:1)
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip
    dxi = 1.0_num/dx
    dyi = 1.0_num/dy
    dzi = 1.0_num/dz
    invvol = dxi*dyi*dzi
    dts2dx = 0.5_num*dt*dxi
    dts2dy = 0.5_num*dt*dyi
    dts2dz = 0.5_num*dt*dzi
    clightsq = 1.0_num/clight**2
    sx=0.0_num;sy=0.0_num;sz=0.0_num;
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;

    ! LOOP ON PARTICLES
    DO ip=1,np
        ! --- computes position in  grid units at (n+1)
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        ! Computes velocity
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)

        ! --- computes particles weights
        wq=q*w(ip)*invvol
        wqx=wq*vx
        wqy=wq*vy
        wqz=wq*vz

        ! Gets position in grid units at (n+1/2) for computing jx(n+1/2)
        xmid=x-dts2dx*vx
        ymid=y-dts2dy*vy
        zmid=z-dts2dz*vz

        ! --- finds node of cell containing particles for current positions
        j=nint(xmid)
        k=nint(ymid)
        l=nint(zmid)
        j0=nint(xmid-0.5_num)
        k0=nint(ymid-0.5_num)
        l0=nint(zmid-0.5_num)
        ! --- computes set of coefficients for node centered quantities
        xint = xmid-j
        yint = ymid-k
        zint = zmid-l
        xintsq = xint*xint
        sx(-1) = 0.5_num*(0.5_num-xint)**2
        sx( 0) = 0.75_num-xintsq
        sx( 1) = 0.5_num*(0.5_num+xint)**2
        yintsq = yint*yint
        sy(-1) = 0.5_num*(0.5_num-yint)**2
        sy( 0) = 0.75_num-yintsq
        sy( 1) = 0.5_num*(0.5_num+yint)**2
        zintsq = zint*zint
        sz(-1) = 0.5_num*(0.5_num-zint)**2
        sz( 0) = (0.75_num-zintsq)
        sz( 1) = 0.5_num*(0.5_num+zint)**2

        ! --- computes set of coefficients for staggered quantities
        xint = xmid-j0-0.5_num
        yint = ymid-k0-0.5_num
        zint = zmid-l0-0.5_num
        xintsq = xint*xint
        sx0(-1) = 0.5_num*(0.5_num-xint)**2
        sx0( 0) = 0.75_num-xintsq
        sx0( 1) = 0.5_num*(0.5_num+xint)**2
        yintsq = yint*yint
        sy0(-1) = 0.5_num*(0.5_num-yint)**2
        sy0( 0) = 0.75_num-yintsq
        sy0( 1) = 0.5_num*(0.5_num+yint)**2
        zintsq = zint*zint
        sz0(-1) = 0.5_num*(0.5_num-zint)**2
        sz0( 0) = (0.75_num-zintsq)
        sz0( 1) = 0.5_num*(0.5_num+zint)**2

        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! --- to the 27 nearest vertices
        ! - JX
        jx(j0-1,k-1,l-1)  = jx(j0-1,k-1,l-1)  +   sx0(-1)*sy(-1)*sz(-1)*wqx
        jx(j0  ,k-1,l-1)  = jx(j0  ,k-1,l-1)  +   sx0(0 )*sy(-1)*sz(-1)*wqx
        jx(j0+1,k-1,l-1)  = jx(j0+1,k-1,l-1)  +   sx0(1 )*sy(-1)*sz(-1)*wqx
        jx(j0-1,k  ,l-1)  = jx(j0-1,k  ,l-1)  +   sx0(-1)*sy(0 )*sz(-1)*wqx
        jx(j0  ,k  ,l-1)  = jx(j0  ,k  ,l-1)  +   sx0(0 )*sy(0 )*sz(-1)*wqx
        jx(j0+1,k  ,l-1)  = jx(j0+1,k  ,l-1)  +   sx0(1 )*sy(0 )*sz(-1)*wqx
        jx(j0-1,k+1,l-1)  = jx(j0-1,k+1,l-1)  +   sx0(-1)*sy(1 )*sz(-1)*wqx
        jx(j0  ,k+1,l-1)  = jx(j0  ,k+1,l-1)  +   sx0(0 )*sy(1 )*sz(-1)*wqx
        jx(j0+1,k+1,l-1)  = jx(j0+1,k+1,l-1)  +   sx0(1 )*sy(1 )*sz(-1)*wqx
        jx(j0-1,k-1,l  )  = jx(j0-1,k-1,l  )  +   sx0(-1)*sy(-1)*sz(0 )*wqx
        jx(j0  ,k-1,l  )  = jx(j0  ,k-1,l  )  +   sx0(0 )*sy(-1)*sz(0 )*wqx
        jx(j0+1,k-1,l  )  = jx(j0+1,k-1,l  )  +   sx0(1 )*sy(-1)*sz(0 )*wqx
        jx(j0-1,k  ,l  )  = jx(j0-1,k  ,l  )  +   sx0(-1)*sy(0 )*sz(0 )*wqx
        jx(j0  ,k  ,l  )  = jx(j0  ,k  ,l  )  +   sx0(0 )*sy(0 )*sz(0 )*wqx
        jx(j0+1,k  ,l  )  = jx(j0+1,k  ,l  )  +   sx0(1 )*sy(0 )*sz(0 )*wqx
        jx(j0-1,k+1,l  )  = jx(j0-1,k+1,l  )  +   sx0(-1)*sy(1 )*sz(0 )*wqx
        jx(j0  ,k+1,l  )  = jx(j0  ,k+1,l  )  +   sx0(0 )*sy(1 )*sz(0 )*wqx
        jx(j0+1,k+1,l  )  = jx(j0+1,k+1,l  )  +   sx0(1 )*sy(1 )*sz(0 )*wqx
        jx(j0-1,k-1,l+1)  = jx(j0-1,k-1,l+1)  +   sx0(-1)*sy(-1)*sz(1 )*wqx
        jx(j0  ,k-1,l+1)  = jx(j0  ,k-1,l+1)  +   sx0(0 )*sy(-1)*sz(1 )*wqx
        jx(j0+1,k-1,l+1)  = jx(j0+1,k-1,l+1)  +   sx0(1 )*sy(-1)*sz(1 )*wqx
        jx(j0-1,k  ,l+1)  = jx(j0-1,k  ,l+1)  +   sx0(-1)*sy(0 )*sz(1 )*wqx
        jx(j0  ,k  ,l+1)  = jx(j0  ,k  ,l+1)  +   sx0(0 )*sy(0 )*sz(1 )*wqx
        jx(j0+1,k  ,l+1)  = jx(j0+1,k  ,l+1)  +   sx0(1 )*sy(0 )*sz(1 )*wqx
        jx(j0-1,k+1,l+1)  = jx(j0-1,k+1,l+1)  +   sx0(-1)*sy(1 )*sz(1 )*wqx
        jx(j0  ,k+1,l+1)  = jx(j0  ,k+1,l+1)  +   sx0(0 )*sy(1 )*sz(1 )*wqx
        jx(j0+1,k+1,l+1)  = jx(j0+1,k+1,l+1)  +   sx0(1 )*sy(1 )*sz(1 )*wqx

!        ! - JY
        jy(j-1,k0-1,l-1)  = jy(j-1,k0-1,l-1)  +   sx(-1)*sy0(-1)*sz(-1)*wqy
        jy(j  ,k0-1,l-1)  = jy(j  ,k0-1,l-1)  +   sx(0 )*sy0(-1)*sz(-1)*wqy
        jy(j+1,k0-1,l-1)  = jy(j+1,k0-1,l-1)  +   sx(1 )*sy0(-1)*sz(-1)*wqy
        jy(j-1,k0  ,l-1)  = jy(j-1,k0  ,l-1)  +   sx(-1)*sy0(0 )*sz(-1)*wqy
        jy(j  ,k0  ,l-1)  = jy(j  ,k0  ,l-1)  +   sx(0 )*sy0(0 )*sz(-1)*wqy
        jy(j+1,k0  ,l-1)  = jy(j+1,k0  ,l-1)  +   sx(1 )*sy0(0 )*sz(-1)*wqy
        jy(j-1,k0+1,l-1)  = jy(j-1,k0+1,l-1)  +   sx(-1)*sy0(1 )*sz(-1)*wqy
        jy(j  ,k0+1,l-1)  = jy(j  ,k0+1,l-1)  +   sx(0 )*sy0(1 )*sz(-1)*wqy
        jy(j+1,k0+1,l-1)  = jy(j+1,k0+1,l-1)  +   sx(1 )*sy0(1 )*sz(-1)*wqy
        jy(j-1,k0-1,l  )  = jy(j-1,k0-1,l  )  +   sx(-1)*sy0(-1)*sz(0 )*wqy
        jy(j  ,k0-1,l  )  = jy(j  ,k0-1,l  )  +   sx(0 )*sy0(-1)*sz(0 )*wqy
        jy(j+1,k0-1,l  )  = jy(j+1,k0-1,l  )  +   sx(1 )*sy0(-1)*sz(0 )*wqy
        jy(j-1,k0  ,l  )  = jy(j-1,k0  ,l  )  +   sx(-1)*sy0(0 )*sz(0 )*wqy
        jy(j  ,k0  ,l  )  = jy(j  ,k0  ,l  )  +   sx(0 )*sy0(0 )*sz(0 )*wqy
        jy(j+1,k0  ,l  )  = jy(j+1,k0  ,l  )  +   sx(1 )*sy0(0 )*sz(0 )*wqy
        jy(j-1,k0+1,l  )  = jy(j-1,k0+1,l  )  +   sx(-1)*sy0(1 )*sz(0 )*wqy
        jy(j  ,k0+1,l  )  = jy(j  ,k0+1,l  )  +   sx(0 )*sy0(1 )*sz(0 )*wqy
        jy(j+1,k0+1,l  )  = jy(j+1,k0+1,l  )  +   sx(1 )*sy0(1 )*sz(0 )*wqy
        jy(j-1,k0-1,l+1)  = jy(j-1,k0-1,l+1)  +   sx(-1)*sy0(-1)*sz(1 )*wqy
        jy(j  ,k0-1,l+1)  = jy(j  ,k0-1,l+1)  +   sx(0 )*sy0(-1)*sz(1 )*wqy
        jy(j+1,k0-1,l+1)  = jy(j+1,k0-1,l+1)  +   sx(1 )*sy0(-1)*sz(1 )*wqy
        jy(j-1,k0  ,l+1)  = jy(j-1,k0  ,l+1)  +   sx(-1)*sy0(0 )*sz(1 )*wqy
        jy(j  ,k0  ,l+1)  = jy(j  ,k0  ,l+1)  +   sx(0 )*sy0(0 )*sz(1 )*wqy
        jy(j+1,k0  ,l+1)  = jy(j+1,k0  ,l+1)  +   sx(1 )*sy0(0 )*sz(1 )*wqy
        jy(j-1,k0+1,l+1)  = jy(j-1,k0+1,l+1)  +   sx(-1)*sy0(1 )*sz(1 )*wqy
        jy(j  ,k0+1,l+1)  = jy(j  ,k0+1,l+1)  +   sx(0 )*sy0(1 )*sz(1 )*wqy
        jy(j+1,k0+1,l+1)  = jy(j+1,k0+1,l+1)  +   sx(1 )*sy0(1 )*sz(1 )*wqy

        ! - JZ
        jz(j-1,k-1,l0-1)  = jz(j-1,k-1,l0-1)  +   sx(-1)*sy(-1)*sz0(-1)*wqz
        jz(j  ,k-1,l0-1)  = jz(j  ,k-1,l0-1)  +   sx(0 )*sy(-1)*sz0(-1)*wqz
        jz(j+1,k-1,l0-1)  = jz(j+1,k-1,l0-1)  +   sx(1 )*sy(-1)*sz0(-1)*wqz
        jz(j-1,k  ,l0-1)  = jz(j-1,k  ,l0-1)  +   sx(-1)*sy(0 )*sz0(-1)*wqz
        jz(j  ,k  ,l0-1)  = jz(j  ,k  ,l0-1)  +   sx(0 )*sy(0 )*sz0(-1)*wqz
        jz(j+1,k  ,l0-1)  = jz(j+1,k  ,l0-1)  +   sx(1 )*sy(0 )*sz0(-1)*wqz
        jz(j-1,k+1,l0-1)  = jz(j-1,k+1,l0-1)  +   sx(-1)*sy(1 )*sz0(-1)*wqz
        jz(j  ,k+1,l0-1)  = jz(j  ,k+1,l0-1)  +   sx(0 )*sy(1 )*sz0(-1)*wqz
        jz(j+1,k+1,l0-1)  = jz(j+1,k+1,l0-1)  +   sx(1 )*sy(1 )*sz0(-1)*wqz
        jz(j-1,k-1,l0  )  = jz(j-1,k-1,l0  )  +   sx(-1)*sy(-1)*sz0(0 )*wqz
        jz(j  ,k-1,l0  )  = jz(j  ,k-1,l0  )  +   sx(0 )*sy(-1)*sz0(0 )*wqz
        jz(j+1,k-1,l0  )  = jz(j+1,k-1,l0  )  +   sx(1 )*sy(-1)*sz0(0 )*wqz
        jz(j-1,k  ,l0  )  = jz(j-1,k  ,l0  )  +   sx(-1)*sy(0 )*sz0(0 )*wqz
        jz(j  ,k  ,l0  )  = jz(j  ,k  ,l0  )  +   sx(0 )*sy(0 )*sz0(0 )*wqz
        jz(j+1,k  ,l0  )  = jz(j+1,k  ,l0  )  +   sx(1 )*sy(0 )*sz0(0 )*wqz
        jz(j-1,k+1,l0  )  = jz(j-1,k+1,l0  )  +   sx(-1)*sy(1 )*sz0(0 )*wqz
        jz(j  ,k+1,l0  )  = jz(j  ,k+1,l0  )  +   sx(0 )*sy(1 )*sz0(0 )*wqz
        jz(j+1,k+1,l0  )  = jz(j+1,k+1,l0  )  +   sx(1 )*sy(1 )*sz0(0 )*wqz
        jz(j-1,k-1,l0+1)  = jz(j-1,k-1,l0+1)  +   sx(-1)*sy(-1)*sz0(1 )*wqz
        jz(j  ,k-1,l0+1)  = jz(j  ,k-1,l0+1)  +   sx(0 )*sy(-1)*sz0(1 )*wqz
        jz(j+1,k-1,l0+1)  = jz(j+1,k-1,l0+1)  +   sx(1 )*sy(-1)*sz0(1 )*wqz
        jz(j-1,k  ,l0+1)  = jz(j-1,k  ,l0+1)  +   sx(-1)*sy(0 )*sz0(1 )*wqz
        jz(j  ,k  ,l0+1)  = jz(j  ,k  ,l0+1)  +   sx(0 )*sy(0 )*sz0(1 )*wqz
        jz(j+1,k  ,l0+1)  = jz(j+1,k  ,l0+1)  +   sx(1 )*sy(0 )*sz0(1 )*wqz
        jz(j-1,k+1,l0+1)  = jz(j-1,k+1,l0+1)  +   sx(-1)*sy(1 )*sz0(1 )*wqz
        jz(j  ,k+1,l0+1)  = jz(j  ,k+1,l0+1)  +   sx(0 )*sy(1 )*sz0(1 )*wqz
        jz(j+1,k+1,l0+1)  = jz(j+1,k+1,l0+1)  +   sx(1 )*sy(1 )*sz0(1 )*wqz
    END DO
    RETURN
END SUBROUTINE depose_jxjyjz_scalar_2_2_2

!!! --- Order 2 3D vector current deposition routine (rho*v)
!!! This versions have good performances on SIMD architectures
!!! Providing that OpenMP 4.0 is available (Directive SIMD)
!!! Use with nox=3
SUBROUTINE depose_jxjyjz_vecHVv2_2_2_2(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    USE precomputed
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid
    REAL(num) ::   wqx,wqy,wqz,ww, wwx, wwy, wwz
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic
    INTEGER(idp) :: nnx, nnxy, n,nn,nv
    INTEGER(idp) :: moff(1:8)
    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL, IG
    REAL(num) :: vx,vy,vz
    REAL(num) :: ww0x(LVEC,4),ww0y(LVEC,4),ww0z(LVEC,4), wwwx(LVEC,8), &
    wwwy(LVEC,8),wwwz(LVEC,8), wq
    REAL(num) :: sx0(LVEC),sx1(LVEC),sx2(LVEC)
    REAL(num) :: sx00(LVEC),sx01(LVEC),sx02(LVEC)
    REAL(num) :: sy0,sy1,sy2,sy00,sy01,sy02
    REAL(num) :: sz0,sz1,sz2,sz00,sz01,sz02, syz
    INTEGER(idp) :: igrid,orig, jorig, korig, lorig
    INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy

    ww0x=0._num; ww0y=0._num; ww0z=0._num
    ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
    ncx=nx+4;ncy=ny+4;ncz=nz+4
    NCELLS=ncx*ncy*ncz
    ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
    jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num
    nnx = nx + 1 + 2*nxguard
    nnxy = nnx*(ny+1+2*nyguard)
    moff = (/-nnx-nnxy,-nnxy,nnx-nnxy,-nnx,nnx,-nnx+nnxy,nnxy,nnx+nnxy/)
    jorig=-2; korig=-2;lorig=-2
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy
    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED sx0:64,sx1:64,sx2:64
        !DIR$ ASSUME_ALIGNED sx00:64,sx01:64,sx02:64
        !DIR$ ASSUME_ALIGNED w:64, wwwx:64,wwwy:64,wwwz:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx = uxp(nn)*gaminv(nn)
            vy = uyp(nn)*gaminv(nn)
            vz = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol
            wqx=wq*vx
            wqy=wq*vy
            wqz=wq*vz

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx
            ymid=y-dts2dy*vy
            zmid=z-dts2dz*vz

            ! --- finds node of cell containing particles for current positions
            j=nint(xmid)
            k=nint(ymid)
            l=nint(zmid)
            j0=nint(xmid-0.5_num)
            k0=nint(ymid-0.5_num)
            l0=nint(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy
            
            IG(n,1) = ICELL(n,1) + ncx + ncxy
            IG(n,2) = ICELL(n,2) + ncx + ncxy
            IG(n,3) = ICELL(n,3) + ncx + ncxy
            
            ! Old IG when centered nodes are directly put in jx,jy,jz
            !IG(n,1)=ICELL(n,1)+(k-korig)*ngx+(l-lorig)*ngxy
            !IG(n,2)=ICELL(n,2)+(k0-korig)*ngx+(l-lorig)*ngxy
            !IG(n,3)=ICELL(n,3)+(k-korig)*ngx+(l0-lorig)*ngxy

            ! --- computes set of coefficients for node centered quantities
            xint = xmid-j
            yint = ymid-k
            zint = zmid-l
            xintsq= xint**2
            yintsq= yint**2
            zintsq= zint**2
            sx0(n)=0.5_num*(0.5_num-xint)**2
            sx1(n)=(0.75_num-xintsq)
            sx2(n)=0.5_num*(0.5_num+xint)**2
            sy0=0.5_num*(0.5_num-yint)**2
            sy1=(0.75_num-yintsq)
            sy2=0.5_num*(0.5_num+yint)**2
            sz0=0.5_num*(0.5_num-zint)**2
            sz1=(0.75_num-zintsq)
            sz2=0.5_num*(0.5_num+zint)**2

            ! --- computes set of coefficients for staggered quantities
            xint = xmid-j0-0.5_num
            yint = ymid-k0-0.5_num
            zint = zmid-l0-0.5_num
            xintsq= xint**2
            yintsq= yint**2
            zintsq= zint**2
            sx00(n)=0.5_num*(0.5_num-xint)**2
            sx01(n)=(0.75_num-xintsq)
            sx02(n)=0.5_num*(0.5_num+xint)**2
            sy00=0.5_num*(0.5_num-yint)**2
            sy01=(0.75_num-yintsq)
            sy02=0.5_num*(0.5_num+yint)**2
            sz00=0.5_num*(0.5_num-zint)**2
            sz01=(0.75_num-zintsq)
            sz02=0.5_num*(0.5_num+zint)**2

            ! -- Weights for planes of 8  vertices
            ! Weights - X
            wwwx(n,1) = sy0*sz0*wqx
            wwwx(n,2) = sy1*sz0*wqx
            wwwx(n,3) = sy2*sz0*wqx
            wwwx(n,4) = sy0*sz1*wqx
            wwwx(n,5) = sy2*sz1*wqx
            wwwx(n,6) = sy0*sz2*wqx
            wwwx(n,7) = sy1*sz2*wqx
            wwwx(n,8) = sy2*sz2*wqx

            ! Weights - Y
            wwwy(n,1) = sy00*sz0*wqy
            wwwy(n,2) = sy01*sz0*wqy
            wwwy(n,3) = sy02*sz0*wqy
            wwwy(n,4) = sy00*sz1*wqy
            wwwy(n,5) = sy02*sz1*wqy
            wwwy(n,6) = sy00*sz2*wqy
            wwwy(n,7) = sy01*sz2*wqy
            wwwy(n,8) = sy02*sz2*wqy

            ! Weights - Z
            wwwz(n,1) = sy0*sz00*wqz
            wwwz(n,2) = sy1*sz00*wqz
            wwwz(n,3) = sy2*sz00*wqz
            wwwz(n,4) = sy0*sz01*wqz
            wwwz(n,5) = sy2*sz01*wqz
            wwwz(n,6) = sy0*sz02*wqz
            wwwz(n,7) = sy1*sz02*wqz
            wwwz(n,8) = sy2*sz02*wqz

            ! -- 3 remaining central points
            syz=sz1*sy1*wqx
            ww0x(n,1)=syz*sx00(n)
            ww0x(n,2)=syz*sx01(n)
            ww0x(n,3)=syz*sx02(n)
            syz=sz1*sy01*wqy
            ww0y(n,1)=syz*sx0(n)
            ww0y(n,2)=syz*sx1(n)
            ww0y(n,3)=syz*sx2(n)
            syz=sz01*sy1*wqz
            ww0z(n,1)=syz*sx0(n)
            ww0z(n,2)=syz*sx1(n)
            ww0z(n,3)=syz*sx2(n)
        END DO
        !$OMP END SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !$OMP SIMD
            DO nv=1,8
                ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
                ! - JX
                wwx=wwwx(n,nv)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)-1) = jxcells(nv,ICELL(n,1)-1) +wwx*sx00(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1))   = jxcells(nv,ICELL(n,1))   +wwx*sx01(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+1) = jxcells(nv,ICELL(n,1)+1) +wwx*sx02(n)
                ! - JY
                wwy=wwwy(n,nv)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)-1) = jycells(nv,ICELL(n,2)-1) +wwy*sx0(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2))   = jycells(nv,ICELL(n,2))   +wwy*sx1(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+1) = jycells(nv,ICELL(n,2)+1) +wwy*sx2(n)
                ! - JZ
                wwz=wwwz(n,nv)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)-1) = jzcells(nv,ICELL(n,3)-1) +wwz*sx0(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3))   = jzcells(nv,ICELL(n,3))   +wwz*sx1(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+1) = jzcells(nv,ICELL(n,3)+1) +wwz*sx2(n)
            END DO
            !$OMP END SIMD
            !$OMP SIMD
            DO nv=1,4
                jxcells(1,IG(n,1)+nv-2) = jxcells(1,IG(n,1)+nv-2) + ww0x(n,nv)
                jycells(1,IG(n,2)+nv-2) = jycells(1,IG(n,2)+nv-2) + ww0y(n,nv)
                jzcells(1,IG(n,3)+nv-2) = jzcells(1,IG(n,3)+nv-2) + ww0z(n,nv)
                
                ! Old but correct and efficient
                !jx(orig+IG(n,1)+nv-2)=jx(orig+IG(n,1)+nv-2)+ww0x(n,nv)
                !jy(orig+IG(n,2)+nv-2)=jy(orig+IG(n,2)+nv-2)+ww0y(n,nv)
                !jz(orig+IG(n,3)+nv-2)=jz(orig+IG(n,3)+nv-2)+ww0z(n,nv)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
    DO iz=1, ncz
        DO iy=1,ncy
            !$OMP SIMD
            DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                ! jx
                jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1,ic)
                jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2,ic)
                jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3,ic)
                jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4,ic)
                jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5,ic)
                jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6,ic)
                jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7,ic)
                jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8,ic)
                ! jy
                jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1,ic)
                jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2,ic)
                jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3,ic)
                jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4,ic)
                jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5,ic)
                jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6,ic)
                jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7,ic)
                jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8,ic)
                ! jz
                jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1,ic)
                jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2,ic)
                jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3,ic)
                jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4,ic)
                jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5,ic)
                jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6,ic)
                jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7,ic)
                jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8,ic)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    DEALLOCATE(jxcells,jycells,jzcells)
    RETURN
END SUBROUTINE depose_jxjyjz_vecHVv2_2_2_2


! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz)
! Order 2 3D vector current deposition routine (rho*v)
! This versions have good performances on SIMD architectures
! Providing that OpenMP 4.0 is available (Directive SIMD)
! This subroutine is similar to depose_jxjyjz_vecHVv2_1_1_1 
! without the reduction process at the end
!
! Inputs:
! - jxcells, jycells, jzcells: transient current arrays
! - np: particle number
! - ncells: number of cells in the tile
! - xp,yp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse Lorentz factor arrays
! - w: particle wight arrays
! - q: charge
! - xmin,ymin,zmin: tile minimum positions
! - dt,dx,dy,dz: time and space steps
! - nx,ny,nz: tile cell numbers in each direction
! - nxguard,nyguard,nzguard: guard cells
! - ncx,ncy,ncz: tile cell extended number (depends on the order)
!
! Ouputs:
! - jxcells, jycells, jzcells updated
! ________________________________________________________________________________________
    USE constants
    USE precomputed
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: np,nx,ny,nz,ncells
    INTEGER(idp), INTENT(IN) :: ncx, ncy, ncz
    INTEGER(idp)             :: nxguard,nyguard,nzguard
    REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    
    REAL(num) :: xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid
    REAL(num) ::   wqx,wqy,wqz,ww, wwx, wwy, wwz
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, ic
    INTEGER(idp) :: nnx, nnxy, n,nn,nv
    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL, IG
    REAL(num) :: vx,vy,vz
    REAL(num) :: ww0x(LVEC,4),ww0y(LVEC,4),ww0z(LVEC,4), wwwx(LVEC,8), &
    wwwy(LVEC,8),wwwz(LVEC,8), wq
    REAL(num) :: sx0(LVEC),sx1(LVEC),sx2(LVEC)
    REAL(num) :: sx00(LVEC),sx01(LVEC),sx02(LVEC)
    REAL(num) :: sy0,sy1,sy2,sy00,sy01,sy02
    REAL(num) :: sz0,sz1,sz2,sz00,sz01,sz02, syz
    INTEGER(idp) :: igrid,orig, jorig, korig, lorig
    INTEGER(idp) ::  ncxy, ix,iy,iz, ngridx, ngridy, ngx, ngxy

    ! __________________________________________________________
    ! Parameters
    
    ww0x=0._num
    ww0y=0._num
    ww0z=0._num
    
    ngridx=nx+1+2*nxguard
    ngridy=ny+1+2*nyguard

    nnx = nx + 1 + 2*nxguard
    nnxy = nnx*(ny+1+2*nyguard)
    
    jorig=-2
    korig=-2
    lorig=-2
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy

    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED sx0:64,sx1:64,sx2:64
        !DIR$ ASSUME_ALIGNED sx00:64,sx01:64,sx02:64
        !DIR$ ASSUME_ALIGNED w:64, wwwx:64,wwwy:64,wwwz:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx = uxp(nn)*gaminv(nn)
            vy = uyp(nn)*gaminv(nn)
            vz = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol
            wqx=wq*vx
            wqy=wq*vy
            wqz=wq*vz

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx
            ymid=y-dts2dy*vy
            zmid=z-dts2dz*vz

            ! --- finds node of cell containing particles for current positions
            j=nint(xmid)
            k=nint(ymid)
            l=nint(zmid)
            j0=nint(xmid-0.5_num)
            k0=nint(ymid-0.5_num)
            l0=nint(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy
            IG(n,1) = ICELL(n,1) + ncx + ncxy
            IG(n,2) = ICELL(n,2) + ncx + ncxy
            IG(n,3) = ICELL(n,3) + ncx + ncxy

            ! --- computes set of coefficients for node centered quantities
            xint = xmid-j
            yint = ymid-k
            zint = zmid-l
            xintsq= xint**2
            yintsq= yint**2
            zintsq= zint**2
            sx0(n)=0.5_num*(0.5_num-xint)**2
            sx1(n)=(0.75_num-xintsq)
            sx2(n)=0.5_num*(0.5_num+xint)**2
            sy0=0.5_num*(0.5_num-yint)**2
            sy1=(0.75_num-yintsq)
            sy2=0.5_num*(0.5_num+yint)**2
            sz0=0.5_num*(0.5_num-zint)**2
            sz1=(0.75_num-zintsq)
            sz2=0.5_num*(0.5_num+zint)**2

            ! --- computes set of coefficients for staggered quantities
            xint = xmid-j0-0.5_num
            yint = ymid-k0-0.5_num
            zint = zmid-l0-0.5_num
            xintsq= xint**2
            yintsq= yint**2
            zintsq= zint**2
            sx00(n)=0.5_num*(0.5_num-xint)**2
            sx01(n)=(0.75_num-xintsq)
            sx02(n)=0.5_num*(0.5_num+xint)**2
            sy00=0.5_num*(0.5_num-yint)**2
            sy01=(0.75_num-yintsq)
            sy02=0.5_num*(0.5_num+yint)**2
            sz00=0.5_num*(0.5_num-zint)**2
            sz01=(0.75_num-zintsq)
            sz02=0.5_num*(0.5_num+zint)**2

            ! -- Weights for planes of 8  vertices
            ! Weights - X
            wwwx(n,1) = sy0*sz0*wqx
            wwwx(n,2) = sy1*sz0*wqx
            wwwx(n,3) = sy2*sz0*wqx
            wwwx(n,4) = sy0*sz1*wqx
            wwwx(n,5) = sy2*sz1*wqx
            wwwx(n,6) = sy0*sz2*wqx
            wwwx(n,7) = sy1*sz2*wqx
            wwwx(n,8) = sy2*sz2*wqx

            ! Weights - Y
            wwwy(n,1) = sy00*sz0*wqy
            wwwy(n,2) = sy01*sz0*wqy
            wwwy(n,3) = sy02*sz0*wqy
            wwwy(n,4) = sy00*sz1*wqy
            wwwy(n,5) = sy02*sz1*wqy
            wwwy(n,6) = sy00*sz2*wqy
            wwwy(n,7) = sy01*sz2*wqy
            wwwy(n,8) = sy02*sz2*wqy

            ! Weights - Z
            wwwz(n,1) = sy0*sz00*wqz
            wwwz(n,2) = sy1*sz00*wqz
            wwwz(n,3) = sy2*sz00*wqz
            wwwz(n,4) = sy0*sz01*wqz
            wwwz(n,5) = sy2*sz01*wqz
            wwwz(n,6) = sy0*sz02*wqz
            wwwz(n,7) = sy1*sz02*wqz
            wwwz(n,8) = sy2*sz02*wqz

            ! -- 3 remaining central points
            syz=sz1*sy1*wqx
            ww0x(n,1)=syz*sx00(n)
            ww0x(n,2)=syz*sx01(n)
            ww0x(n,3)=syz*sx02(n)
            syz=sz1*sy01*wqy
            ww0y(n,1)=syz*sx0(n)
            ww0y(n,2)=syz*sx1(n)
            ww0y(n,3)=syz*sx2(n)
            syz=sz01*sy1*wqz
            ww0z(n,1)=syz*sx0(n)
            ww0z(n,2)=syz*sx1(n)
            ww0z(n,3)=syz*sx2(n)
        END DO
        !$OMP END SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !$OMP SIMD
            DO nv=1,8
                ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
                ! - JX
                wwx=wwwx(n,nv)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)-1) = jxcells(nv,ICELL(n,1)-1) +wwx*sx00(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1))   = jxcells(nv,ICELL(n,1))   +wwx*sx01(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+1) = jxcells(nv,ICELL(n,1)+1) +wwx*sx02(n)
                ! - JY
                wwy=wwwy(n,nv)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)-1) = jycells(nv,ICELL(n,2)-1) +wwy*sx0(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2))   = jycells(nv,ICELL(n,2))   +wwy*sx1(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+1) = jycells(nv,ICELL(n,2)+1) +wwy*sx2(n)
                ! - JZ
                wwz=wwwz(n,nv)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)-1) = jzcells(nv,ICELL(n,3)-1) +wwz*sx0(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3))   = jzcells(nv,ICELL(n,3))   +wwz*sx1(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+1) = jzcells(nv,ICELL(n,3)+1) +wwz*sx2(n)
            END DO
            !$OMP END SIMD
            !$OMP SIMD
            DO nv=1,4
                jxcells(1,IG(n,1)+nv-2) = jxcells(1,IG(n,1)+nv-2) + ww0x(n,nv)
                jycells(1,IG(n,2)+nv-2) = jycells(1,IG(n,2)+nv-2) + ww0y(n,nv)
                jzcells(1,IG(n,3)+nv-2) = jzcells(1,IG(n,3)+nv-2) + ww0z(n,nv)
            END DO
            !$OMP END SIMD
        END DO
    END DO

    RETURN
END SUBROUTINE depose_jxjyjz_vecHV_vnr_2_2_2


!!! --- Order 3 3D scalar current deposition routine (rho*v)
!!! This version does not vectorize on SIMD architectures
SUBROUTINE depose_jxjyjz_scalar_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,vx,vy,vz,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) :: wq, wqx, wqy, wqz, clightsq
    REAL(num), DIMENSION(4) :: sx(-1:2), sy(-1:2), sz(-1:2), sx0(-1:2), sy0(-1:2), sz0(-1:2)
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip
    dxi = 1.0_num/dx
    dyi = 1.0_num/dy
    dzi = 1.0_num/dz
    invvol = dxi*dyi*dzi
    dts2dx = 0.5_num*dt*dxi
    dts2dy = 0.5_num*dt*dyi
    dts2dz = 0.5_num*dt*dzi
    clightsq = 1.0_num/clight**2
    sx=0.0_num;sy=0.0_num;sz=0.0_num;
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num;

    ! LOOP ON PARTICLES
    DO ip=1,np
        ! --- computes position in  grid units at (n+1)
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        ! Computes velocity
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)

        ! --- computes particles weights
        wq=q*w(ip)*invvol
        wqx=wq*vx
        wqy=wq*vy
        wqz=wq*vz

        ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
        xmid=x-dts2dx*vx
        ymid=y-dts2dy*vy
        zmid=z-dts2dz*vz

        ! --- finds node of cell containing particles for current positions
        j=floor(xmid)
        k=floor(ymid)
        l=floor(zmid)
        j0=floor(xmid-0.5_num)
        k0=floor(ymid-0.5_num)
        l0=floor(zmid-0.5_num)
        ! --- computes set of coefficients for node centered quantities
        xint = xmid-j
        yint = ymid-k
        zint = zmid-l
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx(-1) = onesixth*oxintsq*oxint
        sx( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx( 2) = onesixth*xintsq*xint
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy(-1) = onesixth*oyintsq*oyint
        sy( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy( 2) = onesixth*yintsq*yint
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz(-1) = onesixth*ozintsq*ozint
        sz( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz( 2) = onesixth*zintsq*zint

        ! --- computes set of coefficients for staggered quantities
        xint = xmid-j0-0.5_num
        yint = ymid-k0-0.5_num
        zint = zmid-l0-0.5_num
        oxint = 1.0_num-xint
        xintsq = xint*xint
        oxintsq = oxint*oxint
        sx0(-1) = onesixth*oxintsq*oxint
        sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
        sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
        sx0( 2) = onesixth*xintsq*xint
        oyint = 1.0_num-yint
        yintsq = yint*yint
        oyintsq = oyint*oyint
        sy0(-1) = onesixth*oyintsq*oyint
        sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
        sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
        sy0( 2) = onesixth*yintsq*yint
        ozint = 1.0_num-zint
        zintsq = zint*zint
        ozintsq = ozint*ozint
        sz0(-1) = onesixth*ozintsq*ozint
        sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
        sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
        sz0( 2) = onesixth*zintsq*zint

        ! --- add current contributions in the form rho(n+1/2)v(n+1/2)
        ! --- to the 64 nearest vertices
        ! - JX
        jx(j0-1,k-1,l-1)  = jx(j0-1,k-1,l-1)  +   sx0(-1)*sy0(-1)*sz(-1)*wqx
        jx(j0  ,k-1,l-1)  = jx(j0  ,k-1,l-1)  +   sx0(0 )*sy(-1)*sz(-1)*wqx
        jx(j0+1,k-1,l-1)  = jx(j0+1,k-1,l-1)  +   sx0(1 )*sy(-1)*sz(-1)*wqx
        jx(j0+2,k-1,l-1)  = jx(j0+2,k-1,l-1)  +   sx0(2 )*sy(-1)*sz(-1)*wqx
        jx(j0-1,k  ,l-1)  = jx(j0-1,k  ,l-1)  +   sx0(-1)*sy(0 )*sz(-1)*wqx
        jx(j0  ,k  ,l-1)  = jx(j0  ,k  ,l-1)  +   sx0(0 )*sy(0 )*sz(-1)*wqx
        jx(j0+1,k  ,l-1)  = jx(j0+1,k  ,l-1)  +   sx0(1 )*sy(0 )*sz(-1)*wqx
        jx(j0+2,k  ,l-1)  = jx(j0+2,k  ,l-1)  +   sx0(2 )*sy(0 )*sz(-1)*wqx
        jx(j0-1,k+1,l-1)  = jx(j0-1,k+1,l-1)  +   sx0(-1)*sy(1 )*sz(-1)*wqx
        jx(j0  ,k+1,l-1)  = jx(j0  ,k+1,l-1)  +   sx0(0 )*sy(1 )*sz(-1)*wqx
        jx(j0+1,k+1,l-1)  = jx(j0+1,k+1,l-1)  +   sx0(1 )*sy(1 )*sz(-1)*wqx
        jx(j0+2,k+1,l-1)  = jx(j0+2,k+1,l-1)  +   sx0(2 )*sy(1 )*sz(-1)*wqx
        jx(j0-1,k+2,l-1)  = jx(j0-1,k+2,l-1)  +   sx0(-1)*sy(2 )*sz(-1)*wqx
        jx(j0  ,k+2,l-1)  = jx(j0  ,k+2,l-1)  +   sx0(0 )*sy(2 )*sz(-1)*wqx
        jx(j0+1,k+2,l-1)  = jx(j0+1,k+2,l-1)  +   sx0(1 )*sy(2 )*sz(-1)*wqx
        jx(j0+2,k+2,l-1)  = jx(j0+2,k+2,l-1)  +   sx0(2 )*sy(2 )*sz(-1)*wqx
        jx(j0-1,k-1,l  )  = jx(j0-1,k-1,l  )  +   sx0(-1)*sy(-1)*sz(0 )*wqx
        jx(j0  ,k-1,l  )  = jx(j0  ,k-1,l  )  +   sx0(0 )*sy(-1)*sz(0 )*wqx
        jx(j0+1,k-1,l  )  = jx(j0+1,k-1,l  )  +   sx0(1 )*sy(-1)*sz(0 )*wqx
        jx(j0+2,k-1,l  )  = jx(j0+2,k-1,l  )  +   sx0(2 )*sy(-1)*sz(0 )*wqx
        jx(j0-1,k  ,l  )  = jx(j0-1,k  ,l  )  +   sx0(-1)*sy(0 )*sz(0 )*wqx
        jx(j0  ,k  ,l  )  = jx(j0  ,k  ,l  )  +   sx0(0 )*sy(0 )*sz(0 )*wqx
        jx(j0+1,k  ,l  )  = jx(j0+1,k  ,l  )  +   sx0(1 )*sy(0 )*sz(0 )*wqx
        jx(j0+2,k  ,l  )  = jx(j0+2,k  ,l  )  +   sx0(2 )*sy(0 )*sz(0 )*wqx
        jx(j0-1,k+1,l  )  = jx(j0-1,k+1,l  )  +   sx0(-1)*sy(1 )*sz(0 )*wqx
        jx(j0  ,k+1,l  )  = jx(j0  ,k+1,l  )  +   sx0(0 )*sy(1 )*sz(0 )*wqx
        jx(j0+1,k+1,l  )  = jx(j0+1,k+1,l  )  +   sx0(1 )*sy(1 )*sz(0 )*wqx
        jx(j0+2,k+1,l  )  = jx(j0+2,k+1,l  )  +   sx0(2 )*sy(1 )*sz(0 )*wqx
        jx(j0-1,k+2,l  )  = jx(j0-1,k+2,l  )  +   sx0(-1)*sy(2 )*sz(0 )*wqx
        jx(j0  ,k+2,l  )  = jx(j0  ,k+2,l  )  +   sx0(0 )*sy(2 )*sz(0 )*wqx
        jx(j0+1,k+2,l  )  = jx(j0+1,k+2,l  )  +   sx0(1 )*sy(2 )*sz(0 )*wqx
        jx(j0+2,k+2,l  )  = jx(j0+2,k+2,l  )  +   sx0(2 )*sy(2 )*sz(0 )*wqx
        jx(j0-1,k-1,l+1)  = jx(j0-1,k-1,l+1)  +   sx0(-1)*sy(-1)*sz(1 )*wqx
        jx(j0  ,k-1,l+1)  = jx(j0  ,k-1,l+1)  +   sx0(0 )*sy(-1)*sz(1 )*wqx
        jx(j0+1,k-1,l+1)  = jx(j0+1,k-1,l+1)  +   sx0(1 )*sy(-1)*sz(1 )*wqx
        jx(j0+2,k-1,l+1)  = jx(j0+2,k-1,l+1)  +   sx0(2 )*sy(-1)*sz(1 )*wqx
        jx(j0-1,k  ,l+1)  = jx(j0-1,k  ,l+1)  +   sx0(-1)*sy(0 )*sz(1 )*wqx
        jx(j0  ,k,  l+1)  = jx(j0  ,k  ,l+1)  +   sx0(0 )*sy(0 )*sz(1 )*wqx
        jx(j0+1,k  ,l+1)  = jx(j0+1,k  ,l+1)  +   sx0(1 )*sy(0 )*sz(1 )*wqx
        jx(j0+2,k  ,l+1)  = jx(j0+2,k  ,l+1)  +   sx0(2 )*sy(0 )*sz(1 )*wqx
        jx(j0-1,k+1,l+1)  = jx(j0-1,k+1,l+1)  +   sx0(-1)*sy(1 )*sz(1 )*wqx
        jx(j0  ,k+1,l+1)  = jx(j0  ,k+1,l+1)  +   sx0(0 )*sy(1 )*sz(1 )*wqx
        jx(j0+1,k+1,l+1)  = jx(j0+1,k+1,l+1)  +   sx0(1 )*sy(1 )*sz(1 )*wqx
        jx(j0+2,k+1,l+1)  = jx(j0+2,k+1,l+1)  +   sx0(2 )*sy(1 )*sz(1 )*wqx
        jx(j0-1,k+2,l+1)  = jx(j0-1,k+2,l+1)  +   sx0(-1)*sy(2 )*sz(1 )*wqx
        jx(j0  ,k+2,l+1)  = jx(j0  ,k+2,l+1)  +   sx0(0 )*sy(2 )*sz(1 )*wqx
        jx(j0+1,k+2,l+1)  = jx(j0+1,k+2,l+1)  +   sx0(1 )*sy(2 )*sz(1 )*wqx
        jx(j0+2,k+2,l+1)  = jx(j0+2,k+2,l+1)  +   sx0(2 )*sy(2 )*sz(1 )*wqx
        jx(j0-1,k-1,l+2)  = jx(j0-1,k-1,l+2)  +   sx0(-1)*sy(-1)*sz(2 )*wqx
        jx(j0  ,k-1,l+2)  = jx(j0  ,k-1,l+2)  +   sx0(0 )*sy(-1)*sz(2 )*wqx
        jx(j0+1,k-1,l+2)  = jx(j0+1,k-1,l+2)  +   sx0(1 )*sy(-1)*sz(2 )*wqx
        jx(j0+2,k-1,l+2)  = jx(j0+2,k-1,l+2)  +   sx0(2 )*sy(-1)*sz(2 )*wqx
        jx(j0-1,k  ,l+2)  = jx(j0-1,k  ,l+2)  +   sx0(-1)*sy(0 )*sz(2 )*wqx
        jx(j0  ,k  ,l+2)  = jx(j0  ,k  ,l+2)  +   sx0(0 )*sy(0 )*sz(2 )*wqx
        jx(j0+1,k  ,l+2)  = jx(j0+1,k  ,l+2)  +   sx0(1 )*sy(0 )*sz(2 )*wqx
        jx(j0+2,k  ,l+2)  = jx(j0+2,k  ,l+2)  +   sx0(2 )*sy(0 )*sz(2 )*wqx
        jx(j0-1,k+1,l+2)  = jx(j0-1,k+1,l+2)  +   sx0(-1)*sy(1 )*sz(2 )*wqx
        jx(j0  ,k+1,l+2)  = jx(j0  ,k+1,l+2)  +   sx0(0 )*sy(1 )*sz(2 )*wqx
        jx(j0+1,k+1,l+2)  = jx(j0+1,k+1,l+2)  +   sx0(1 )*sy(1 )*sz(2 )*wqx
        jx(j0+2,k+1,l+2)  = jx(j0+2,k+1,l+2)  +   sx0(2 )*sy(1 )*sz(2 )*wqx
        jx(j0-1,k+2,l+2)  = jx(j0-1,k+2,l+2)  +   sx0(-1)*sy(2 )*sz(2 )*wqx
        jx(j0  ,k+2,l+2)  = jx(j0  ,k+2,l+2)  +   sx0(0 )*sy(2 )*sz(2 )*wqx
        jx(j0+1,k+2,l+2)  = jx(j0+1,k+2,l+2)  +   sx0(1 )*sy(2 )*sz(2 )*wqx
        jx(j0+2,k+2,l+2)  = jx(j0+2,k+2,l+2)  +   sx0(2 )*sy(2 )*sz(2 )*wqx

        ! - JY
        jy(j-1,k0-1,l-1)  = jy(j-1,k0-1,l-1)  +   sx(-1)*sy0(-1)*sz(-1)*wqy
        jy(j  ,k0-1,l-1)  = jy(j  ,k0-1,l-1)  +   sx(0 )*sy0(-1)*sz(-1)*wqy
        jy(j+1,k0-1,l-1)  = jy(j+1,k0-1,l-1)  +   sx(1 )*sy0(-1)*sz(-1)*wqy
        jy(j+2,k0-1,l-1)  = jy(j+2,k0-1,l-1)  +   sx(2 )*sy0(-1)*sz(-1)*wqy
        jy(j-1,k0  ,l-1)  = jy(j-1,k0  ,l-1)  +   sx(-1)*sy0(0 )*sz(-1)*wqy
        jy(j  ,k0  ,l-1)  = jy(j  ,k0  ,l-1)  +   sx(0 )*sy0(0 )*sz(-1)*wqy
        jy(j+1,k0  ,l-1)  = jy(j+1,k0  ,l-1)  +   sx(1 )*sy0(0 )*sz(-1)*wqy
        jy(j+2,k0  ,l-1)  = jy(j+2,k0  ,l-1)  +   sx(2 )*sy0(0 )*sz(-1)*wqy
        jy(j-1,k0+1,l-1)  = jy(j-1,k0+1,l-1)  +   sx(-1)*sy0(1 )*sz(-1)*wqy
        jy(j  ,k0+1,l-1)  = jy(j  ,k0+1,l-1)  +   sx(0 )*sy0(1 )*sz(-1)*wqy
        jy(j+1,k0+1,l-1)  = jy(j+1,k0+1,l-1)  +   sx(1 )*sy0(1 )*sz(-1)*wqy
        jy(j+2,k0+1,l-1)  = jy(j+2,k0+1,l-1)  +   sx(2 )*sy0(1 )*sz(-1)*wqy
        jy(j-1,k0+2,l-1)  = jy(j-1,k0+2,l-1)  +   sx(-1)*sy0(2 )*sz(-1)*wqy
        jy(j  ,k0+2,l-1)  = jy(j  ,k0+2,l-1)  +   sx(0 )*sy0(2 )*sz(-1)*wqy
        jy(j+1,k0+2,l-1)  = jy(j+1,k0+2,l-1)  +   sx(1 )*sy0(2 )*sz(-1)*wqy
        jy(j+2,k0+2,l-1)  = jy(j+2,k0+2,l-1)  +   sx(2 )*sy0(2 )*sz(-1)*wqy
        jy(j-1,k0-1,l  )  = jy(j-1,k0-1,l  )  +   sx(-1)*sy0(-1)*sz(0 )*wqy
        jy(j  ,k0-1,l  )  = jy(j  ,k0-1,l  )  +   sx(0 )*sy0(-1)*sz(0 )*wqy
        jy(j+1,k0-1,l  )  = jy(j+1,k0-1,l  )  +   sx(1 )*sy0(-1)*sz(0 )*wqy
        jy(j+2,k0-1,l  )  = jy(j+2,k0-1,l  )  +   sx(2 )*sy0(-1)*sz(0 )*wqy
        jy(j-1,k0  ,l  )  = jy(j-1,k0  ,l  )  +   sx(-1)*sy0(0 )*sz(0 )*wqy
        jy(j  ,k0  ,l  )  = jy(j  ,k0  ,l  )  +   sx(0 )*sy0(0 )*sz(0 )*wqy
        jy(j+1,k0  ,l  )  = jy(j+1,k0  ,l  )  +   sx(1 )*sy0(0 )*sz(0 )*wqy
        jy(j+2,k0  ,l  )  = jy(j+2,k0  ,l  )  +   sx(2 )*sy0(0 )*sz(0 )*wqy
        jy(j-1,k0+1,l  )  = jy(j-1,k0+1,l  )  +   sx(-1)*sy0(1 )*sz(0 )*wqy
        jy(j  ,k0+1,l  )  = jy(j  ,k0+1,l  )  +   sx(0 )*sy0(1 )*sz(0 )*wqy
        jy(j+1,k0+1,l  )  = jy(j+1,k0+1,l  )  +   sx(1 )*sy0(1 )*sz(0 )*wqy
        jy(j+2,k0+1,l  )  = jy(j+2,k0+1,l  )  +   sx(2 )*sy0(1 )*sz(0 )*wqy
        jy(j-1,k0+2,l  )  = jy(j-1,k0+2,l  )  +   sx(-1)*sy0(2 )*sz(0 )*wqy
        jy(j  ,k0+2,l  )  = jy(j  ,k0+2,l  )  +   sx(0 )*sy0(2 )*sz(0 )*wqy
        jy(j+1,k0+2,l  )  = jy(j+1,k0+2,l  )  +   sx(1 )*sy0(2 )*sz(0 )*wqy
        jy(j+2,k0+2,l  )  = jy(j+2,k0+2,l  )  +   sx(2 )*sy0(2 )*sz(0 )*wqy
        jy(j-1,k0-1,l+1)  = jy(j-1,k0-1,l+1)  +   sx(-1)*sy0(-1)*sz(1 )*wqy
        jy(j  ,k0-1,l+1)  = jy(j  ,k0-1,l+1)  +   sx(0 )*sy0(-1)*sz(1 )*wqy
        jy(j+1,k0-1,l+1)  = jy(j+1,k0-1,l+1)  +   sx(1 )*sy0(-1)*sz(1 )*wqy
        jy(j+2,k0-1,l+1)  = jy(j+2,k0-1,l+1)  +   sx(2 )*sy0(-1)*sz(1 )*wqy
        jy(j-1,k0  ,l+1)  = jy(j-1,k0  ,l+1)  +   sx(-1)*sy0(0 )*sz(1 )*wqy
        jy(j  ,k0,  l+1)  = jy(j  ,k0  ,l+1)  +   sx(0 )*sy0(0 )*sz(1 )*wqy
        jy(j+1,k0  ,l+1)  = jy(j+1,k0  ,l+1)  +   sx(1 )*sy0(0 )*sz(1 )*wqy
        jy(j+2,k0  ,l+1)  = jy(j+2,k0  ,l+1)  +   sx(2 )*sy0(0 )*sz(1 )*wqy
        jy(j-1,k0+1,l+1)  = jy(j-1,k0+1,l+1)  +   sx(-1)*sy0(1 )*sz(1 )*wqy
        jy(j  ,k0+1,l+1)  = jy(j  ,k0+1,l+1)  +   sx(0 )*sy0(1 )*sz(1 )*wqy
        jy(j+1,k0+1,l+1)  = jy(j+1,k0+1,l+1)  +   sx(1 )*sy0(1 )*sz(1 )*wqy
        jy(j+2,k0+1,l+1)  = jy(j+2,k0+1,l+1)  +   sx(2 )*sy0(1 )*sz(1 )*wqy
        jy(j-1,k0+2,l+1)  = jy(j-1,k0+2,l+1)  +   sx(-1)*sy0(2 )*sz(1 )*wqy
        jy(j  ,k0+2,l+1)  = jy(j  ,k0+2,l+1)  +   sx(0 )*sy0(2 )*sz(1 )*wqy
        jy(j+1,k0+2,l+1)  = jy(j+1,k0+2,l+1)  +   sx(1 )*sy0(2 )*sz(1 )*wqy
        jy(j+2,k0+2,l+1)  = jy(j+2,k0+2,l+1)  +   sx(2 )*sy0(2 )*sz(1 )*wqy
        jy(j-1,k0-1,l+2)  = jy(j-1,k0-1,l+2)  +   sx(-1)*sy0(-1)*sz(2 )*wqy
        jy(j  ,k0-1,l+2)  = jy(j  ,k0-1,l+2)  +   sx(0 )*sy0(-1)*sz(2 )*wqy
        jy(j+1,k0-1,l+2)  = jy(j+1,k0-1,l+2)  +   sx(1 )*sy0(-1)*sz(2 )*wqy
        jy(j+2,k0-1,l+2)  = jy(j+2,k0-1,l+2)  +   sx(2 )*sy0(-1)*sz(2 )*wqy
        jy(j-1,k0  ,l+2)  = jy(j-1,k0  ,l+2)  +   sx(-1)*sy0(0 )*sz(2 )*wqy
        jy(j  ,k0  ,l+2)  = jy(j  ,k0  ,l+2)  +   sx(0 )*sy0(0 )*sz(2 )*wqy
        jy(j+1,k0  ,l+2)  = jy(j+1,k0  ,l+2)  +   sx(1 )*sy0(0 )*sz(2 )*wqy
        jy(j+2,k0  ,l+2)  = jy(j+2,k0  ,l+2)  +   sx(2 )*sy0(0 )*sz(2 )*wqy
        jy(j-1,k0+1,l+2)  = jy(j-1,k0+1,l+2)  +   sx(-1)*sy0(1 )*sz(2 )*wqy
        jy(j  ,k0+1,l+2)  = jy(j  ,k0+1,l+2)  +   sx(0 )*sy0(1 )*sz(2 )*wqy
        jy(j+1,k0+1,l+2)  = jy(j+1,k0+1,l+2)  +   sx(1 )*sy0(1 )*sz(2 )*wqy
        jy(j+2,k0+1,l+2)  = jy(j+2,k0+1,l+2)  +   sx(2 )*sy0(1 )*sz(2 )*wqy
        jy(j-1,k0+2,l+2)  = jy(j-1,k0+2,l+2)  +   sx(-1)*sy0(2 )*sz(2 )*wqy
        jy(j  ,k0+2,l+2)  = jy(j  ,k0+2,l+2)  +   sx(0 )*sy0(2 )*sz(2 )*wqy
        jy(j+1,k0+2,l+2)  = jy(j+1,k0+2,l+2)  +   sx(1 )*sy0(2 )*sz(2 )*wqy
        jy(j+2,k0+2,l+2)  = jy(j+2,k0+2,l+2)  +   sx(2 )*sy0(2 )*sz(2 )*wqy

        ! - JZ
        jz(j-1,k-1,l0-1)  = jz(j-1,k-1,l0-1)  +   sx(-1)*sy(-1)*sz0(-1)*wqz
        jz(j  ,k-1,l0-1)  = jz(j  ,k-1,l0-1)  +   sx(0 )*sy(-1)*sz0(-1)*wqz
        jz(j+1,k-1,l0-1)  = jz(j+1,k-1,l0-1)  +   sx(1 )*sy(-1)*sz0(-1)*wqz
        jz(j+2,k-1,l0-1)  = jz(j+2,k-1,l0-1)  +   sx(2 )*sy(-1)*sz0(-1)*wqz
        jz(j-1,k  ,l0-1)  = jz(j-1,k  ,l0-1)  +   sx(-1)*sy(0 )*sz0(-1)*wqz
        jz(j  ,k  ,l0-1)  = jz(j  ,k  ,l0-1)  +   sx(0 )*sy(0 )*sz0(-1)*wqz
        jz(j+1,k  ,l0-1)  = jz(j+1,k  ,l0-1)  +   sx(1 )*sy(0 )*sz0(-1)*wqz
        jz(j+2,k  ,l0-1)  = jz(j+2,k  ,l0-1)  +   sx(2 )*sy(0 )*sz0(-1)*wqz
        jz(j-1,k+1,l0-1)  = jz(j-1,k+1,l0-1)  +   sx(-1)*sy(1 )*sz0(-1)*wqz
        jz(j  ,k+1,l0-1)  = jz(j  ,k+1,l0-1)  +   sx(0 )*sy(1 )*sz0(-1)*wqz
        jz(j+1,k+1,l0-1)  = jz(j+1,k+1,l0-1)  +   sx(1 )*sy(1 )*sz0(-1)*wqz
        jz(j+2,k+1,l0-1)  = jz(j+2,k+1,l0-1)  +   sx(2 )*sy(1 )*sz0(-1)*wqz
        jz(j-1,k+2,l0-1)  = jz(j-1,k+2,l0-1)  +   sx(-1)*sy(2 )*sz0(-1)*wqz
        jz(j  ,k+2,l0-1)  = jz(j  ,k+2,l0-1)  +   sx(0 )*sy(2 )*sz0(-1)*wqz
        jz(j+1,k+2,l0-1)  = jz(j+1,k+2,l0-1)  +   sx(1 )*sy(2 )*sz0(-1)*wqz
        jz(j+2,k+2,l0-1)  = jz(j+2,k+2,l0-1)  +   sx(2 )*sy(2 )*sz0(-1)*wqz
        jz(j-1,k-1,l0  )  = jz(j-1,k-1,l0  )  +   sx(-1)*sy(-1)*sz0(0 )*wqz
        jz(j  ,k-1,l0  )  = jz(j  ,k-1,l0  )  +   sx(0 )*sy(-1)*sz0(0 )*wqz
        jz(j+1,k-1,l0  )  = jz(j+1,k-1,l0  )  +   sx(1 )*sy(-1)*sz0(0 )*wqz
        jz(j+2,k-1,l0  )  = jz(j+2,k-1,l0  )  +   sx(2 )*sy(-1)*sz0(0 )*wqz
        jz(j-1,k  ,l0  )  = jz(j-1,k  ,l0  )  +   sx(-1)*sy(0 )*sz0(0 )*wqz
        jz(j  ,k  ,l0  )  = jz(j  ,k  ,l0  )  +   sx(0 )*sy(0 )*sz0(0 )*wqz
        jz(j+1,k  ,l0  )  = jz(j+1,k  ,l0  )  +   sx(1 )*sy(0 )*sz0(0 )*wqz
        jz(j+2,k  ,l0  )  = jz(j+2,k  ,l0  )  +   sx(2 )*sy(0 )*sz0(0 )*wqz
        jz(j-1,k+1,l0  )  = jz(j-1,k+1,l0  )  +   sx(-1)*sy(1 )*sz0(0 )*wqz
        jz(j  ,k+1,l0  )  = jz(j  ,k+1,l0  )  +   sx(0 )*sy(1 )*sz0(0 )*wqz
        jz(j+1,k+1,l0  )  = jz(j+1,k+1,l0  )  +   sx(1 )*sy(1 )*sz0(0 )*wqz
        jz(j+2,k+1,l0  )  = jz(j+2,k+1,l0  )  +   sx(2 )*sy(1 )*sz0(0 )*wqz
        jz(j-1,k+2,l0  )  = jz(j-1,k+2,l0  )  +   sx(-1)*sy(2 )*sz0(0 )*wqz
        jz(j  ,k+2,l0  )  = jz(j  ,k+2,l0  )  +   sx(0 )*sy(2 )*sz0(0 )*wqz
        jz(j+1,k+2,l0  )  = jz(j+1,k+2,l0  )  +   sx(1 )*sy(2 )*sz0(0 )*wqz
        jz(j+2,k+2,l0  )  = jz(j+2,k+2,l0  )  +   sx(2 )*sy(2 )*sz0(0 )*wqz
        jz(j-1,k-1,l0+1)  = jz(j-1,k-1,l0+1)  +   sx(-1)*sy(-1)*sz0(1 )*wqz
        jz(j  ,k-1,l0+1)  = jz(j  ,k-1,l0+1)  +   sx(0 )*sy(-1)*sz0(1 )*wqz
        jz(j+1,k-1,l0+1)  = jz(j+1,k-1,l0+1)  +   sx(1 )*sy(-1)*sz0(1 )*wqz
        jz(j+2,k-1,l0+1)  = jz(j+2,k-1,l0+1)  +   sx(2 )*sy(-1)*sz0(1 )*wqz
        jz(j-1,k  ,l0+1)  = jz(j-1,k  ,l0+1)  +   sx(-1)*sy(0 )*sz0(1 )*wqz
        jz(j  ,k,  l0+1)  = jz(j  ,k  ,l0+1)  +   sx(0 )*sy(0 )*sz0(1 )*wqz
        jz(j+1,k  ,l0+1)  = jz(j+1,k  ,l0+1)  +   sx(1 )*sy(0 )*sz0(1 )*wqz
        jz(j+2,k  ,l0+1)  = jz(j+2,k  ,l0+1)  +   sx(2 )*sy(0 )*sz0(1 )*wqz
        jz(j-1,k+1,l0+1)  = jz(j-1,k+1,l0+1)  +   sx(-1)*sy(1 )*sz0(1 )*wqz
        jz(j  ,k+1,l0+1)  = jz(j  ,k+1,l0+1)  +   sx(0 )*sy(1 )*sz0(1 )*wqz
        jz(j+1,k+1,l0+1)  = jz(j+1,k+1,l0+1)  +   sx(1 )*sy(1 )*sz0(1 )*wqz
        jz(j+2,k+1,l0+1)  = jz(j+2,k+1,l0+1)  +   sx(2 )*sy(1 )*sz0(1 )*wqz
        jz(j-1,k+2,l0+1)  = jz(j-1,k+2,l0+1)  +   sx(-1)*sy(2 )*sz0(1 )*wqz
        jz(j  ,k+2,l0+1)  = jz(j  ,k+2,l0+1)  +   sx(0 )*sy(2 )*sz0(1 )*wqz
        jz(j+1,k+2,l0+1)  = jz(j+1,k+2,l0+1)  +   sx(1 )*sy(2 )*sz0(1 )*wqz
        jz(j+2,k+2,l0+1)  = jz(j+2,k+2,l0+1)  +   sx(2 )*sy(2 )*sz0(1 )*wqz
        jz(j-1,k-1,l0+2)  = jz(j-1,k-1,l0+2)  +   sx(-1)*sy(-1)*sz0(2 )*wqz
        jz(j  ,k-1,l0+2)  = jz(j  ,k-1,l0+2)  +   sx(0 )*sy(-1)*sz0(2 )*wqz
        jz(j+1,k-1,l0+2)  = jz(j+1,k-1,l0+2)  +   sx(1 )*sy(-1)*sz0(2 )*wqz
        jz(j+2,k-1,l0+2)  = jz(j+2,k-1,l0+2)  +   sx(2 )*sy(-1)*sz0(2 )*wqz
        jz(j-1,k  ,l0+2)  = jz(j-1,k  ,l0+2)  +   sx(-1)*sy(0 )*sz0(2 )*wqz
        jz(j  ,k  ,l0+2)  = jz(j  ,k  ,l0+2)  +   sx(0 )*sy(0 )*sz0(2 )*wqz
        jz(j+1,k  ,l0+2)  = jz(j+1,k  ,l0+2)  +   sx(1 )*sy(0 )*sz0(2 )*wqz
        jz(j+2,k  ,l0+2)  = jz(j+2,k  ,l0+2)  +   sx(2 )*sy(0 )*sz0(2 )*wqz
        jz(j-1,k+1,l0+2)  = jz(j-1,k+1,l0+2)  +   sx(-1)*sy(1 )*sz0(2 )*wqz
        jz(j  ,k+1,l0+2)  = jz(j  ,k+1,l0+2)  +   sx(0 )*sy(1 )*sz0(2 )*wqz
        jz(j+1,k+1,l0+2)  = jz(j+1,k+1,l0+2)  +   sx(1 )*sy(1 )*sz0(2 )*wqz
        jz(j+2,k+1,l0+2)  = jz(j+2,k+1,l0+2)  +   sx(2 )*sy(1 )*sz0(2 )*wqz
        jz(j-1,k+2,l0+2)  = jz(j-1,k+2,l0+2)  +   sx(-1)*sy(2 )*sz0(2 )*wqz
        jz(j  ,k+2,l0+2)  = jz(j  ,k+2,l0+2)  +   sx(0 )*sy(2 )*sz0(2 )*wqz
        jz(j+1,k+2,l0+2)  = jz(j+1,k+2,l0+2)  +   sx(1 )*sy(2 )*sz0(2 )*wqz
        jz(j+2,k+2,l0+2)  = jz(j+2,k+2,l0+2)  +   sx(2 )*sy(2 )*sz0(2 )*wqz
    END DO
    RETURN
END SUBROUTINE depose_jxjyjz_scalar_3_3_3

!!! --- Order 3 3D vector current deposition routine (rho*v)
!!! This versions have good performances on SIMD architectures
!!! Providing that OpenMP 4.0 is available (Directive SIMD)
!!! Use with nox=4
SUBROUTINE depose_jxjyjz_vecHVv2_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w,gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint, &
                   oxint,oyint,xintsq,yintsq,oxintsq,oyintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) ::   ww, wwx, wwy, wwz, usq, clightsq
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic, ix, iy, iz
    INTEGER(idp) :: nnx, nnxy,ngridx, ngridy, n,nn,nv
    INTEGER(idp) :: moff(1:8)
    INTEGER(idp), PARAMETER :: LVEC2=32
    REAL(num) :: zint(LVEC),zint0(LVEC2)
    INTEGER(idp), DIMENSION(LVEC2,3) :: ICELL
    REAL(num), DIMENSION(LVEC2) :: vx,vy,vz
    REAL(num) ::  wwwx(LVEC2,16), wwwy(LVEC2,16),wwwz(LVEC2,16), wq
    REAL(num) :: sx1(LVEC2),sx2(LVEC2),sx3(LVEC2),sx4(LVEC2)
    REAL(num) :: sx01(LVEC2),sx02(LVEC2),sx03(LVEC2),sx04(LVEC2)
    REAL(num) :: sy1(LVEC2),sy2(LVEC2),sy3(LVEC2),sy4(LVEC2)
    REAL(num) :: sy01(LVEC2),sy02(LVEC2),sy03(LVEC2),sy04(LVEC2)
    REAL(num), DIMENSION(4) :: szz, zdec, h1, h11, h12, sgn
    INTEGER(idp) :: orig, ncxy, ncx, ncy, ncz, ngx, ngxy, igrid, jorig, korig, lorig

    dxi = 1.0_num/dx
    dyi = 1.0_num/dy
    dzi = 1.0_num/dz
    invvol = dxi*dyi*dzi
    dts2dx = 0.5_num*dt*dxi
    dts2dy = 0.5_num*dt*dyi
    dts2dz = 0.5_num*dt*dzi
    clightsq = 1.0_num/clight**2
    ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
    ncx=nx+5; ncy=ny+5; ncz=nz+5
    NCELLS=ncx*ncy*ncz
    ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
    jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num;
    nnx = ngridx
    nnxy = ngridx*ngridy
    moff = (/-nnxy,0_idp,nnxy,2_idp*nnxy,nnx-nnxy,nnx,nnx+nnxy,nnx+2_idp*nnxy/)
    jorig=-3_idp; korig=-3_idp;lorig=-3_idp
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy

    h1=(/1_num,0_num,1_num,0_num/); sgn=(/1_num,-1_num,1_num,-1_num/)
    h11=(/0_num,1_num,1_num,0_num/); h12=(/1_num,0_num,0_num,1_num/)
    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC2
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED vx:64,vy:64,vz:64
        !DIR$ ASSUME_ALIGNED sx1:64,sx2:64,sx3:64,sx4:64
        !DIR$ ASSUME_ALIGNED sy1:64,sy2:64,sy3:64,sy4:64
        !DIR$ ASSUME_ALIGNED sx01:64,sx02:64,sx03:64,sx04:64
        !DIR$ ASSUME_ALIGNED sy01:64,sy02:64,sy03:64,sy04:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC2,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx(n) = uxp(nn)*gaminv(nn)
            vy(n) = uyp(nn)*gaminv(nn)
            vz(n) = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx(n)
            ymid=y-dts2dy*vy(n)
            zmid=z-dts2dz*vz(n)

            ! --- finds node of cell containing particles for current positions
            j=floor(xmid)
            k=floor(ymid)
            l=floor(zmid)
            j0=floor(xmid-0.5_num)
            k0=floor(ymid-0.5_num)
            l0=floor(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

            ! --- computes set of coefficients for node centered quantities
            xint    = xmid-j
            yint    = ymid-k
            zint(n) = zmid-l
            oxint   = 1.0_num-xint
            xintsq  = xint*xint
            oxintsq = oxint*oxint
            sx1(n)  = onesixth*oxintsq*oxint
            sx2(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
            sx3(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
            sx4(n)  = onesixth*xintsq*xint
            oyint   = 1.0_num-yint
            yintsq  = yint*yint
            oyintsq = oyint*oyint
            sy1(n)  = onesixth*oyintsq*oyint*wq
            sy2(n)  = (twothird-yintsq*(1.0_num-yint*0.5_num))*wq
            sy3(n)  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))*wq
            sy4(n)  = onesixth*yintsq*yint*wq

            ! --- computes set of coefficients for staggered quantities
            xint     = xmid-j0-0.5_num
            yint     = ymid-k0-0.5_num
            zint0(n) = zmid-l0-0.5_num
            oxint    = 1.0_num-xint
            xintsq   = xint*xint
            oxintsq  = oxint*oxint
            sx01(n)  = onesixth*oxintsq*oxint
            sx02(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
            sx03(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
            sx04(n)  = onesixth*xintsq*xint
            oyint    = 1.0_num-yint
            yintsq   = yint*yint
            oyintsq  = oyint*oyint
            sy01(n)  = onesixth*oyintsq*oyint*wq
            sy02(n)  = (twothird-yintsq*(1.0_num-yint*0.5_num))*wq
            sy03(n)  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))*wq
            sy04(n)  = onesixth*yintsq*yint*wq
        END DO
        !$OMP END SIMD

        ! Compute weights
        DO n=1,MIN(LVEC2,np-ip+1)
            !DIR$ ASSUME_ALIGNED w:64, wwwx:64,wwwy:64,wwwz:64
            !$OMP SIMD
            DO nv=1,4 !!! Vector
                ! - Weiths for jx
                zdec(nv)      = (h1(nv)-zint(n))*sgn(nv)
                szz(nv)       = (twothird-zdec(nv)**2*(1.0_num-zdec(nv)*0.5_num))*h11(nv) &
                              +onesixth*zdec(nv)**3*h12(nv)
                wwwx(nv,n)    = szz(nv)*sy1(n)*vx(n)
                wwwx(nv+4,n)  = szz(nv)*sy2(n)*vx(n)
                wwwx(nv+8,n)  = szz(nv)*sy3(n)*vx(n)
                wwwx(nv+12,n) = szz(nv)*sy4(n)*vx(n)
                ! - Weiths for jy
                wwwy(nv,n)    = szz(nv)*sy01(n)*vy(n)
                wwwy(nv+4,n)  = szz(nv)*sy02(n)*vy(n)
                wwwy(nv+8,n)  = szz(nv)*sy03(n)*vy(n)
                wwwy(nv+12,n) = szz(nv)*sy04(n)*vy(n)
                ! - Weiths for jz
                zdec(nv)      = (h1(nv)-zint0(n))*sgn(nv)
                szz(nv)       = (twothird-zdec(nv)**2*(1.0_num-zdec(nv)*0.5_num))*h11(nv) &
                              +onesixth*zdec(nv)**3*h12(nv)
                wwwz(nv,n)    = szz(nv)*sy1(n)*vz(n)
                wwwz(nv+4,n)  = szz(nv)*sy2(n)*vz(n)
                wwwz(nv+8,n)  = szz(nv)*sy3(n)*vz(n)
                wwwz(nv+12,n) = szz(nv)*sy4(n)*vz(n)

            ENDDO
            !$OMP END SIMD
        END DO

        ! Add weights to nearest vertices
        DO n=1,MIN(LVEC2,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !$OMP SIMD
            DO nv=1,8
                ! --- JX
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)-ncx-1) = jxcells(nv,ICELL(n,1)-ncx-1) + wwwx(nv,n)*sx01(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)-ncx)   = jxcells(nv,ICELL(n,1)-ncx)   + wwwx(nv,n)*sx02(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+1) = jxcells(nv,ICELL(n,1)-ncx+1) + wwwx(nv,n)*sx03(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+2) = jxcells(nv,ICELL(n,1)-ncx+2) + wwwx(nv,n)*sx04(n)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)+ncx-1) = jxcells(nv,ICELL(n,1)+ncx-1) + wwwx(nv+8,n)*sx01(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)+ncx)   = jxcells(nv,ICELL(n,1)+ncx)   + wwwx(nv+8,n)*sx02(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+1) = jxcells(nv,ICELL(n,1)+ncx+1) + wwwx(nv+8,n)*sx03(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+2) = jxcells(nv,ICELL(n,1)+ncx+2) + wwwx(nv+8,n)*sx04(n)

                ! --- JY
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)-ncx-1) = jycells(nv,ICELL(n,2)-ncx-1) + wwwy(nv,n)*sx1(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)-ncx)   = jycells(nv,ICELL(n,2)-ncx)   + wwwy(nv,n)*sx2(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+1) = jycells(nv,ICELL(n,2)-ncx+1) + wwwy(nv,n)*sx3(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+2) = jycells(nv,ICELL(n,2)-ncx+2) + wwwy(nv,n)*sx4(n)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)+ncx-1) = jycells(nv,ICELL(n,2)+ncx-1) + wwwy(nv+8,n)*sx1(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)+ncx)   = jycells(nv,ICELL(n,2)+ncx)   + wwwy(nv+8,n)*sx2(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+1) = jycells(nv,ICELL(n,2)+ncx+1) + wwwy(nv+8,n)*sx3(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+2) = jycells(nv,ICELL(n,2)+ncx+2) + wwwy(nv+8,n)*sx4(n)

                ! --- JZ
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)-ncx-1) = jzcells(nv,ICELL(n,3)-ncx-1) + wwwz(nv,n)*sx1(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)-ncx)   = jzcells(nv,ICELL(n,3)-ncx)   + wwwz(nv,n)*sx2(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+1) = jzcells(nv,ICELL(n,3)-ncx+1) + wwwz(nv,n)*sx3(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+2) = jzcells(nv,ICELL(n,3)-ncx+2) + wwwz(nv,n)*sx4(n)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)+ncx-1) = jzcells(nv,ICELL(n,3)+ncx-1) + wwwz(nv+8,n)*sx1(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)+ncx)   = jzcells(nv,ICELL(n,3)+ncx)   + wwwz(nv+8,n)*sx2(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+1) = jzcells(nv,ICELL(n,3)+ncx+1) + wwwz(nv+8,n)*sx3(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+2) = jzcells(nv,ICELL(n,3)+ncx+2) + wwwz(nv+8,n)*sx4(n)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
    DO iz=1, ncz
        DO iy=1,ncy
            !$OMP SIMD
            DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                ! jx
                jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1,ic)
                jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2,ic)
                jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3,ic)
                jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4,ic)
                jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5,ic)
                jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6,ic)
                jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7,ic)
                jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8,ic)
                ! jy
                jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1,ic)
                jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2,ic)
                jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3,ic)
                jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4,ic)
                jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5,ic)
                jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6,ic)
                jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7,ic)
                jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8,ic)
                ! jz
                jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1,ic)
                jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2,ic)
                jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3,ic)
                jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4,ic)
                jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5,ic)
                jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6,ic)
                jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7,ic)
                jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8,ic)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    DEALLOCATE(jxcells,jycells,jzcells)
    RETURN
END SUBROUTINE depose_jxjyjz_vecHVv2_3_3_3


!!! --- Order 3 3D vector current deposition routine (rho*v)
!!! This versions have good performances on SIMD architectures
!!! Providing that OpenMP 4.0 is available (Directive SIMD)
!!! Use with nox=3
SUBROUTINE depose_jxjyjz_vecHVv3_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard)
    USE constants
    IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard
    REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: dxi,dyi,dzi,xint,yint,zint, &
                 oxint,oyint,ozint,xintsq,yintsq,zintsq, oxintsq,oyintsq, ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid,invvol, dts2dx, dts2dy, dts2dz
    REAL(num) ::   ww, wwx, wwy, wwz, clightsq
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic, ix, iy, iz
    INTEGER(idp) :: nnx, nnxy,ngridx, ngridy, n,nn,nv
    INTEGER(idp) :: moff(1:8)
    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL
    REAL(num), DIMENSION(LVEC) :: vx,vy,vz
    REAL(num) ::  wwwx(LVEC,16), wwwy(LVEC,16),wwwz(LVEC,16), wq
    REAL(num) :: sx1(LVEC),sx2(LVEC),sx3(LVEC),sx4(LVEC)
    REAL(num) :: sx01(LVEC),sx02(LVEC),sx03(LVEC),sx04(LVEC)
    REAL(num) :: sy1,sy2,sy3,sy4,sz1,sz2,sz3,sz4
    REAL(num) :: sy01,sy02,sy03,sy04,sz01,sz02,sz03,sz04
    REAL(num), DIMENSION(4) :: szz, zdec, h1, h11, h12, sgn
    REAL(num):: wwwx1(LVEC,8),wwwx2(LVEC,8),wwwy1(LVEC,8),wwwy2(LVEC,8),wwwz1(LVEC,8),wwwz2(LVEC,8)
    REAL(num):: wx1,wx2,wy1,wy2,wz1,wz2
    INTEGER(idp) :: orig, ncxy, ncx, ncy, ncz, ngx, ngxy, igrid, jorig, korig, lorig

    dxi = 1.0_num/dx
    dyi = 1.0_num/dy
    dzi = 1.0_num/dz
    invvol = dxi*dyi*dzi
    dts2dx = 0.5_num*dt*dxi
    dts2dy = 0.5_num*dt*dyi
    dts2dz = 0.5_num*dt*dzi
    clightsq = 1.0_num/clight**2
    ngridx=nx+1+2*nxguard;ngridy=ny+1+2*nyguard
    ncx=nx+5; ncy=ny+4; ncz=nz+3
    NCELLS=ncx*ncy*ncz
    ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
    jxcells=0.0_num; jycells=0.0_num; jzcells=0.0_num;
    nnx = ngridx
    nnxy = ngridx*ngridy
    moff = (/-nnxy,0_idp,nnxy,2_idp*nnxy,nnx-nnxy,nnx,nnx+nnxy,nnx+2_idp*nnxy/)
    jorig=-2_idp; korig=-2_idp;lorig=-2_idp
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy

    h1=(/1_num,0_num,1_num,0_num/); sgn=(/1_num,-1_num,1_num,-1_num/)
    h11=(/0_num,1_num,1_num,0_num/); h12=(/1_num,0_num,0_num,1_num/)
    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED vx:64,vy:64,vz:64
        !DIR$ ASSUME_ALIGNED sx1:64,sx2:64,sx3:64,sx4:64
        !DIR$ ASSUME_ALIGNED sx01:64,sx02:64,sx03:64,sx04:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx(n) = uxp(nn)*gaminv(nn)
            vy(n) = uyp(nn)*gaminv(nn)
            vz(n) = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx(n)
            ymid=y-dts2dy*vy(n)
            zmid=z-dts2dz*vz(n)

            ! --- finds node of cell containing particles for current positions
            j=floor(xmid)
            k=floor(ymid)
            l=floor(zmid)
            j0=floor(xmid-0.5_num)
            k0=floor(ymid-0.5_num)
            l0=floor(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

            ! --- computes set of coefficients for node centered quantities
            xint    = xmid-j
            yint    = ymid-k
            zint    = zmid-l
            oxint   = 1.0_num-xint
            xintsq  = xint*xint
            oxintsq = oxint*oxint
            sx1(n)  = onesixth*oxintsq*oxint
            sx2(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
            sx3(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
            sx4(n)  = onesixth*xintsq*xint
            oyint   = 1.0_num-yint
            yintsq  = yint*yint
            oyintsq = oyint*oyint
            sy1  = onesixth*oyintsq*oyint
            sy2  = (twothird-yintsq*(1.0_num-yint*0.5_num))
            sy3  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
            sy4  = onesixth*yintsq*yint
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz1 = onesixth*ozintsq*ozint*wq
            sz2 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
            sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
            sz4 = onesixth*zintsq*zint*wq

            ! --- computes set of coefficients for staggered quantities
            xint     = xmid-j0-0.5_num
            yint     = ymid-k0-0.5_num
            zint     = zmid-l0-0.5_num
            oxint    = 1.0_num-xint
            xintsq   = xint*xint
            oxintsq  = oxint*oxint
            sx01(n)  = onesixth*oxintsq*oxint
            sx02(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
            sx03(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
            sx04(n)  = onesixth*xintsq*xint
            oyint    = 1.0_num-yint
            yintsq   = yint*yint
            oyintsq  = oyint*oyint
            sy01  = onesixth*oyintsq*oyint
            sy02  = (twothird-yintsq*(1.0_num-yint*0.5_num))
            sy03  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
            sy04  = onesixth*yintsq*yint
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz01 = onesixth*ozintsq*ozint*wq
            sz02 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
            sz03 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
            sz04 = onesixth*zintsq*zint*wq
            ! --- computes weights
            ! - X
            wwwx1(n,1)=sz1*sy1
            wwwx1(n,2)=sz2*sy1
            wwwx1(n,3)=sz3*sy1
            wwwx1(n,4)=sz4*sy1
            wwwx1(n,5)=sz1*sy2
            wwwx1(n,6)=sz2*sy2
            wwwx1(n,7)=sz3*sy2
            wwwx1(n,8)=sz4*sy2
            wwwx2(n,1)=sz1*sy3
            wwwx2(n,2)=sz2*sy3
            wwwx2(n,3)=sz3*sy3
            wwwx2(n,4)=sz4*sy3
            wwwx2(n,5)=sz1*sy4
            wwwx2(n,6)=sz2*sy4
            wwwx2(n,7)=sz3*sy4
            wwwx2(n,8)=sz4*sy4
            ! - Y
            wwwy1(n,1)=sz1*sy01
            wwwy1(n,2)=sz2*sy01
            wwwy1(n,3)=sz3*sy01
            wwwy1(n,4)=sz4*sy01
            wwwy1(n,5)=sz1*sy02
            wwwy1(n,6)=sz2*sy02
            wwwy1(n,7)=sz3*sy02
            wwwy1(n,8)=sz4*sy02
            wwwy2(n,1)=sz1*sy03
            wwwy2(n,2)=sz2*sy03
            wwwy2(n,3)=sz3*sy03
            wwwy2(n,4)=sz4*sy03
            wwwy2(n,5)=sz1*sy04
            wwwy2(n,6)=sz2*sy04
            wwwy2(n,7)=sz3*sy04
            wwwy2(n,8)=sz4*sy04
            ! - Z
            wwwz1(n,1)=sz01*sy1
            wwwz1(n,2)=sz02*sy1
            wwwz1(n,3)=sz03*sy1
            wwwz1(n,4)=sz04*sy1
            wwwz1(n,5)=sz01*sy2
            wwwz1(n,6)=sz02*sy2
            wwwz1(n,7)=sz03*sy2
            wwwz1(n,8)=sz04*sy2
            wwwz2(n,1)=sz01*sy3
            wwwz2(n,2)=sz02*sy3
            wwwz2(n,3)=sz03*sy3
            wwwz2(n,4)=sz04*sy3
            wwwz2(n,5)=sz01*sy4
            wwwz2(n,6)=sz02*sy4
            wwwz2(n,7)=sz03*sy4
            wwwz2(n,8)=sz04*sy4
        END DO
        !$OMP END SIMD

        ! Add weights to nearest vertices
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !$OMP SIMD
            DO nv=1,8
                ! --- JX
                wx1=wwwx1(n,nv); wx2=wwwx2(n,nv)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)-ncx-1) = jxcells(nv,ICELL(n,1)-ncx-1) + wx1*sx01(n)*vx(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)-ncx)   = jxcells(nv,ICELL(n,1)-ncx)   + wx1*sx02(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+1) = jxcells(nv,ICELL(n,1)-ncx+1) + wx1*sx03(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+2) = jxcells(nv,ICELL(n,1)-ncx+2) + wx1*sx04(n)*vx(n)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)+ncx-1) = jxcells(nv,ICELL(n,1)+ncx-1) + wx2*sx01(n)*vx(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)+ncx)   = jxcells(nv,ICELL(n,1)+ncx)   + wx2*sx02(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+1) = jxcells(nv,ICELL(n,1)+ncx+1) + wx2*sx03(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+2) = jxcells(nv,ICELL(n,1)+ncx+2) + wx2*sx04(n)*vx(n)

                ! --- JY
                wy1=wwwy1(n,nv); wy2=wwwy2(n,nv)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)-ncx-1) = jycells(nv,ICELL(n,2)-ncx-1) + wy1*sx1(n)*vy(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)-ncx)   = jycells(nv,ICELL(n,2)-ncx)   + wy1*sx2(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+1) = jycells(nv,ICELL(n,2)-ncx+1) + wy1*sx3(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+2) = jycells(nv,ICELL(n,2)-ncx+2) + wy1*sx4(n)*vy(n)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)+ncx-1) = jycells(nv,ICELL(n,2)+ncx-1) + wy2*sx1(n)*vy(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)+ncx)   = jycells(nv,ICELL(n,2)+ncx)   + wy2*sx2(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+1) = jycells(nv,ICELL(n,2)+ncx+1) + wy2*sx3(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+2) = jycells(nv,ICELL(n,2)+ncx+2) + wy2*sx4(n)*vy(n)

                ! --- JZ
                wz1=wwwz1(n,nv); wz2=wwwz2(n,nv)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)-ncx-1) = jzcells(nv,ICELL(n,3)-ncx-1) + wz1*sx1(n)*vz(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)-ncx)   = jzcells(nv,ICELL(n,3)-ncx)   + wz1*sx2(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+1) = jzcells(nv,ICELL(n,3)-ncx+1) + wz1*sx3(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+2) = jzcells(nv,ICELL(n,3)-ncx+2) + wz1*sx4(n)*vz(n)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)+ncx-1) = jzcells(nv,ICELL(n,3)+ncx-1) + wz2*sx1(n)*vz(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)+ncx)   = jzcells(nv,ICELL(n,3)+ncx)   + wz2*sx2(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+1) = jzcells(nv,ICELL(n,3)+ncx+1) + wz2*sx3(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+2) = jzcells(nv,ICELL(n,3)+ncx+2) + wz2*sx4(n)*vz(n)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
    DO iz=1, ncz
        DO iy=1,ncy
            !$OMP SIMD
            DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                ! jx
                jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1,ic)
                jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2,ic)
                jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3,ic)
                jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4,ic)
                jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5,ic)
                jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6,ic)
                jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7,ic)
                jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8,ic)
                ! jy
                jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1,ic)
                jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2,ic)
                jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3,ic)
                jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4,ic)
                jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5,ic)
                jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6,ic)
                jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7,ic)
                jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8,ic)
                ! jz
                jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1,ic)
                jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2,ic)
                jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3,ic)
                jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4,ic)
                jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5,ic)
                jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6,ic)
                jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7,ic)
                jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8,ic)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    DEALLOCATE(jxcells,jycells,jzcells)
    RETURN
END SUBROUTINE depose_jxjyjz_vecHVv3_3_3_3

! ________________________________________________________________________________________
SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3(jxcells,jycells,jzcells,np,ncells,xp,yp,zp,&
           uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
           dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz)
! Order 3 3D vector current deposition routine (rho*v)
! This versions have good performances on SIMD architectures
! Providing that OpenMP 4.0 is available (Directive SIMD)
! This subroutine is similar to depose_jxjyjz_vecHVv2_1_1_1 
! without the reduction process at the end
!
! Inputs:
! - jxcells, jycells, jzcells: transient current arrays
! - np: particle number
! - ncells: number of cells in the tile
! - xp,yp,zp: particle position arrays
! - uxp,uyp,uzp: particle momentum arrays
! - gaminv: inverse Lorentz factor arrays
! - w: particle wight arrays
! - q: charge
! - xmin,ymin,zmin: tile minimum positions
! - dt,dx,dy,dz: time and space steps
! - nx,ny,nz: tile cell numbers in each direction
! - nxguard,nyguard,nzguard: guard cells
! - ncx,ncy,ncz: tile cell extended number (depends on the order)
!
! Ouputs:
! - jxcells, jycells, jzcells updated
! ________________________________________________________________________________________
    USE constants
    USE precomputed
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)                      :: np,nx,ny,nz
    INTEGER(idp), INTENT(IN)                      :: nxguard,nyguard,nzguard
    REAL(num), DIMENSION(8,ncells), INTENT(INOUT) :: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np), INTENT(IN) :: xp,yp,zp,uxp,uyp,uzp, gaminv, w
    REAL(num), INTENT(IN) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    INTEGER(idp), INTENT(IN)                      :: ncx, ncy, ncz

    REAL(num) :: xint,yint,zint, &
                 oxint,oyint,ozint,xintsq,yintsq,zintsq, oxintsq,oyintsq, ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid
    REAL(num) ::   ww, wwx, wwy, wwz
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic, ix, iy, iz
    INTEGER(idp) :: nnx, nnxy,ngridx, ngridy, n,nn,nv
    INTEGER(idp) :: moff(1:8)
    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL
    REAL(num), DIMENSION(LVEC) :: vx,vy,vz
    REAL(num) ::  wwwx(LVEC,16), wwwy(LVEC,16),wwwz(LVEC,16), wq
    REAL(num) :: sx1(LVEC),sx2(LVEC),sx3(LVEC),sx4(LVEC)
    REAL(num) :: sx01(LVEC),sx02(LVEC),sx03(LVEC),sx04(LVEC)
    REAL(num) :: sy1,sy2,sy3,sy4,sz1,sz2,sz3,sz4
    REAL(num) :: sy01,sy02,sy03,sy04,sz01,sz02,sz03,sz04
    REAL(num), DIMENSION(4) :: szz, zdec, h1, h11, h12, sgn
    REAL(num):: wwwx1(LVEC,8),wwwx2(LVEC,8),wwwy1(LVEC,8),wwwy2(LVEC,8),wwwz1(LVEC,8),wwwz2(LVEC,8)
    REAL(num):: wx1,wx2,wy1,wy2,wz1,wz2
    INTEGER(idp) :: orig, ncxy, ngx, ngxy, igrid, jorig, korig, lorig

    ! ___________________________________________________________
    ! Parameters
    ngridx=nx+1+2*nxguard
    ngridy=ny+1+2*nyguard
    
    nnx = ngridx
    nnxy = ngridx*ngridy
    moff = (/-nnxy,0_idp,nnxy,2_idp*nnxy,nnx-nnxy,nnx,nnx+nnxy,nnx+2_idp*nnxy/)   
    jorig=-2_idp
    korig=-2_idp
    lorig=-2_idp
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy

    h1=(/1_num,0_num,1_num,0_num/)
    sgn=(/1_num,-1_num,1_num,-1_num/)
    h11=(/0_num,1_num,1_num,0_num/)
    h12=(/1_num,0_num,0_num,1_num/)

    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED vx:64,vy:64,vz:64
        !DIR$ ASSUME_ALIGNED sx1:64,sx2:64,sx3:64,sx4:64
        !DIR$ ASSUME_ALIGNED sx01:64,sx02:64,sx03:64,sx04:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! Computes velocity
            vx(n) = uxp(nn)*gaminv(nn)
            vy(n) = uyp(nn)*gaminv(nn)
            vz(n) = uzp(nn)*gaminv(nn)

            ! --- computes particles weights
            wq=q*w(nn)*invvol

            ! Gets position in grid units at (n+1/2) for computing rho(n+1/2)
            xmid=x-dts2dx*vx(n)
            ymid=y-dts2dy*vy(n)
            zmid=z-dts2dz*vz(n)

            ! --- finds node of cell containing particles for current positions
            j=floor(xmid)
            k=floor(ymid)
            l=floor(zmid)
            j0=floor(xmid-0.5_num)
            k0=floor(ymid-0.5_num)
            l0=floor(zmid-0.5_num)
            ICELL(n,1)=1+(j0-jorig)+(k-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

            ! --- computes set of coefficients for node centered quantities
            xint    = xmid-j
            yint    = ymid-k
            zint    = zmid-l
            oxint   = 1.0_num-xint
            xintsq  = xint*xint
            oxintsq = oxint*oxint
            sx1(n)  = onesixth*oxintsq*oxint
            sx2(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
            sx3(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
            sx4(n)  = onesixth*xintsq*xint
            oyint   = 1.0_num-yint
            yintsq  = yint*yint
            oyintsq = oyint*oyint
            sy1  = onesixth*oyintsq*oyint
            sy2  = (twothird-yintsq*(1.0_num-yint*0.5_num))
            sy3  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
            sy4  = onesixth*yintsq*yint
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz1 = onesixth*ozintsq*ozint*wq
            sz2 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
            sz3 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
            sz4 = onesixth*zintsq*zint*wq

            ! --- computes set of coefficients for staggered quantities
            xint     = xmid-j0-0.5_num
            yint     = ymid-k0-0.5_num
            zint     = zmid-l0-0.5_num
            oxint    = 1.0_num-xint
            xintsq   = xint*xint
            oxintsq  = oxint*oxint
            sx01(n)  = onesixth*oxintsq*oxint
            sx02(n)  = twothird-xintsq*(1.0_num-xint*0.5_num)
            sx03(n)  = twothird-oxintsq*(1.0_num-oxint*0.5_num)
            sx04(n)  = onesixth*xintsq*xint
            oyint    = 1.0_num-yint
            yintsq   = yint*yint
            oyintsq  = oyint*oyint
            sy01  = onesixth*oyintsq*oyint
            sy02  = (twothird-yintsq*(1.0_num-yint*0.5_num))
            sy03  = (twothird-oyintsq*(1.0_num-oyint*0.5_num))
            sy04  = onesixth*yintsq*yint
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz01 = onesixth*ozintsq*ozint*wq
            sz02 = (twothird-zintsq*(1.0_num-zint*0.5_num))*wq
            sz03 = (twothird-ozintsq*(1.0_num-ozint*0.5_num))*wq
            sz04 = onesixth*zintsq*zint*wq
            ! --- computes weights
            ! - X
            wwwx1(n,1)=sz1*sy1
            wwwx1(n,2)=sz2*sy1
            wwwx1(n,3)=sz3*sy1
            wwwx1(n,4)=sz4*sy1
            wwwx1(n,5)=sz1*sy2
            wwwx1(n,6)=sz2*sy2
            wwwx1(n,7)=sz3*sy2
            wwwx1(n,8)=sz4*sy2
            wwwx2(n,1)=sz1*sy3
            wwwx2(n,2)=sz2*sy3
            wwwx2(n,3)=sz3*sy3
            wwwx2(n,4)=sz4*sy3
            wwwx2(n,5)=sz1*sy4
            wwwx2(n,6)=sz2*sy4
            wwwx2(n,7)=sz3*sy4
            wwwx2(n,8)=sz4*sy4
            ! - Y
            wwwy1(n,1)=sz1*sy01
            wwwy1(n,2)=sz2*sy01
            wwwy1(n,3)=sz3*sy01
            wwwy1(n,4)=sz4*sy01
            wwwy1(n,5)=sz1*sy02
            wwwy1(n,6)=sz2*sy02
            wwwy1(n,7)=sz3*sy02
            wwwy1(n,8)=sz4*sy02
            wwwy2(n,1)=sz1*sy03
            wwwy2(n,2)=sz2*sy03
            wwwy2(n,3)=sz3*sy03
            wwwy2(n,4)=sz4*sy03
            wwwy2(n,5)=sz1*sy04
            wwwy2(n,6)=sz2*sy04
            wwwy2(n,7)=sz3*sy04
            wwwy2(n,8)=sz4*sy04
            ! - Z
            wwwz1(n,1)=sz01*sy1
            wwwz1(n,2)=sz02*sy1
            wwwz1(n,3)=sz03*sy1
            wwwz1(n,4)=sz04*sy1
            wwwz1(n,5)=sz01*sy2
            wwwz1(n,6)=sz02*sy2
            wwwz1(n,7)=sz03*sy2
            wwwz1(n,8)=sz04*sy2
            wwwz2(n,1)=sz01*sy3
            wwwz2(n,2)=sz02*sy3
            wwwz2(n,3)=sz03*sy3
            wwwz2(n,4)=sz04*sy3
            wwwz2(n,5)=sz01*sy4
            wwwz2(n,6)=sz02*sy4
            wwwz2(n,7)=sz03*sy4
            wwwz2(n,8)=sz04*sy4
        END DO
        !$OMP END SIMD

        ! Add weights to nearest vertices
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !$OMP SIMD
            DO nv=1,8
                ! --- JX
                wx1=wwwx1(n,nv); wx2=wwwx2(n,nv)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)-ncx-1) = jxcells(nv,ICELL(n,1)-ncx-1) + wx1*sx01(n)*vx(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)-ncx)   = jxcells(nv,ICELL(n,1)-ncx)   + wx1*sx02(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+1) = jxcells(nv,ICELL(n,1)-ncx+1) + wx1*sx03(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+2) = jxcells(nv,ICELL(n,1)-ncx+2) + wx1*sx04(n)*vx(n)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)+ncx-1) = jxcells(nv,ICELL(n,1)+ncx-1) + wx2*sx01(n)*vx(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)+ncx)   = jxcells(nv,ICELL(n,1)+ncx)   + wx2*sx02(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+1) = jxcells(nv,ICELL(n,1)+ncx+1) + wx2*sx03(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+2) = jxcells(nv,ICELL(n,1)+ncx+2) + wx2*sx04(n)*vx(n)

                ! --- JY
                wy1=wwwy1(n,nv); wy2=wwwy2(n,nv)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)-ncx-1) = jycells(nv,ICELL(n,2)-ncx-1) + wy1*sx1(n)*vy(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)-ncx)   = jycells(nv,ICELL(n,2)-ncx)   + wy1*sx2(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+1) = jycells(nv,ICELL(n,2)-ncx+1) + wy1*sx3(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+2) = jycells(nv,ICELL(n,2)-ncx+2) + wy1*sx4(n)*vy(n)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)+ncx-1) = jycells(nv,ICELL(n,2)+ncx-1) + wy2*sx1(n)*vy(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)+ncx)   = jycells(nv,ICELL(n,2)+ncx)   + wy2*sx2(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+1) = jycells(nv,ICELL(n,2)+ncx+1) + wy2*sx3(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+2) = jycells(nv,ICELL(n,2)+ncx+2) + wy2*sx4(n)*vy(n)

                ! --- JZ
                wz1=wwwz1(n,nv); wz2=wwwz2(n,nv)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)-ncx-1) = jzcells(nv,ICELL(n,3)-ncx-1) + wz1*sx1(n)*vz(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)-ncx)   = jzcells(nv,ICELL(n,3)-ncx)   + wz1*sx2(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+1) = jzcells(nv,ICELL(n,3)-ncx+1) + wz1*sx3(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+2) = jzcells(nv,ICELL(n,3)-ncx+2) + wz1*sx4(n)*vz(n)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)+ncx-1) = jzcells(nv,ICELL(n,3)+ncx-1) + wz2*sx1(n)*vz(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)+ncx)   = jzcells(nv,ICELL(n,3)+ncx)   + wz2*sx2(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+1) = jzcells(nv,ICELL(n,3)+ncx+1) + wz2*sx3(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+2) = jzcells(nv,ICELL(n,3)+ncx+2) + wz2*sx4(n)*vz(n)
            END DO
            !$OMP END SIMD
        END DO
    END DO

    RETURN
END SUBROUTINE depose_jxjyjz_vecHV_vnr_3_3_3



!===========================================================================================
! Esirkepov current deposition algorithm at order 1 in x, y, z (nox=noy=noz=1)
SUBROUTINE depose_jxjyjz_esirkepov_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                      dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,l_particles_weight,l4symtry)
!===========================================================================================
USE omp_lib
USE constants
IMPLICIT NONE
INTEGER(idp):: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
REAL(num) :: clghtisq,usq,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz, &
                                      s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,         &
                                      oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                                      dtsdx0,dtsdy0,dtsdz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx
REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
INTEGER(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                                      ixmin, ixmax, iymin, iymax, izmin, izmax
LOGICAL(idp) :: l_particles_weight,l4symtry

! PARAMETER INIT
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
dtsdx0 = dt*dxi
dtsdy0 = dt*dyi
dtsdz0 = dt*dzi
invvol = 1.0_num/(dx*dy*dz)
invdtdx = 1.0_num/(dt*dy*dz)
invdtdy = 1.0_num/(dt*dx*dz)
invdtdz = 1.0_num/(dt*dx*dy)
ALLOCATE(sdx(-1:2,-1:2,-1:2),sdy(-1:2,-1:2,-1:2),sdz(-1:2,-1:2,-1:2))
ALLOCATE(sx(-1:2), sx0(-1:2), dsx(-1:2))
ALLOCATE(sy(-1:2), sy0(-1:2), dsy(-1:2))
ALLOCATE(sz(-1:2), sz0(-1:2), dsz(-1:2))
clghtisq = 1.0_num/clight**2
dtsdz0 = dt*dzi
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num
!!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ip,x,y,z,usq,vx,vy,vz,gaminv,xold,yold,zold, &
!!$OMP wq,wqx,wqy,wqz,iixp0,ijxp0,ikxp0, xint,yint,zint, oxint,xintsq, oxintsq,dix,diy,diz, &
!!$OMP dsx, dsy, dsz, oyint,yintsq, oyintsq, ozint,zintsq, ozintsq,ixmin, ixmax, iymin, iymax, izmin, izmax,  &
!!$OMP k,j,i,kc,jc,ic, iixp, ijxp, ikxp,sx,sy,sz, sx0,sy0,sz0,sdx,sdy,sdz,jx1,jy1,jz1) &
!!$OMP SHARED(np,xp,yp,zp,uxp,uyp,uzp,w,dxi,dyi,dzi,invdtdx,invdtdy,invdtdz,xmin,ymin,zmin,clghtisq,dtsdx0,dtsdy0,dtsdz0,q,jx,jy,jz)
!!$OMP DO
DO ip=1,np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqy = wq*invdtdy
    wqz = wq*invdtdz
    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ijxp0=floor(y)
    ikxp0=floor(z)
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    yint=y-ijxp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
    sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
    sy0( 0) = 1.0_num-yint
    sy0( 1) = yint
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint
    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ijxp=floor(yold)
    ikxp=floor(zold)
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    yint = yold-ijxp
    zint = zold-ikxp
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diy = ijxp-ijxp0
    diz = ikxp-ikxp0
    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx=0.0_num;sy=0.0_num;sz=0.0_num
    ! --- computes coefficients for quantities centered between nodes
    sx( 0+dix) = 1.0_num-xint
    sx( 1+dix) = xint
    sy( 0+diy) = 1.0_num-yint
    sy( 1+diy) = yint
    sz( 0+diz) = 1.0_num-zint
    sz( 1+diz) = zint
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsy = sy - sy0
    dsz = sz - sz0
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)
    ixmax = max(0,dix)+1
    iymin = min(0,diy)
    iymax = max(0,diy)+1
    izmin = min(0,diz)
    izmax = max(0,diz)+1
    
    ! --- add current contributions
    DO k=izmin, izmax
      DO j=iymin, iymax
          DO i=ixmin, ixmax
              ic = iixp0+i
              jc = ijxp0+j
              kc = ikxp0+k
              IF(i<ixmax) THEN
                  sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
                  (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
                  IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                  jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)
              END IF
              IF(j<iymax) THEN
                  sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
                  (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
                  IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                  jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)
              END IF
              IF(k<izmax) THEN
                  sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
                  (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
                  IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                  jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)
              END IF
          END DO
      END DO
    END DO

ENDDO


!!$OMP END DO
!!$OMP CRITICAL
!jx=jx+jx1
!jy=jy+jy1
!jz=jz+jz1
!!$OMP END CRITICAL
!!$OMP END PARALLEL
DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)
RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_1_1_1



!===========================================================================================
! Esirkepov current deposition algorithm at order 2 in x, y, z (nox=noy=noz=2)
SUBROUTINE depose_jxjyjz_esirkepov_2_2_2(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                      dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,l_particles_weight,l4symtry)
!===========================================================================================
  USE omp_lib
  USE constants
  IMPLICIT NONE
  INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
  REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
  REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
  REAL(num) :: clghtisq,usq,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz, &
                                      s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,         &
                                      oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                                      dtsdx0,dtsdy0,dtsdz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx
REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
INTEGER :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                                      ixmin, ixmax, iymin, iymax, izmin, izmax
LOGICAL(idp) :: l_particles_weight,l4symtry

! PARAMETER INIT
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
dtsdx0 = dt*dxi
dtsdy0 = dt*dyi
dtsdz0 = dt*dzi
invvol = 1.0_num/(dx*dy*dz)
invdtdx = 1.0_num/(dt*dy*dz)
invdtdy = 1.0_num/(dt*dx*dz)
invdtdz = 1.0_num/(dt*dx*dy)
ALLOCATE(sdx(-2:2,-2:2,-2:2),sdy(-2:2,-2:2,-2:2),sdz(-2:2,-2:2,-2:2))
ALLOCATE(sx(-2:2), sx0(-2:2), dsx(-2:2))
ALLOCATE(sy(-2:2), sy0(-2:2), dsy(-2:2))
ALLOCATE(sz(-2:2), sz0(-2:2), dsz(-2:2))
clghtisq = 1.0_num/clight**2
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num

dtsdz0 = dt*dzi
DO ip=1,np
  ! --- computes current position in grid units
  x = (xp(ip)-xmin)*dxi
  y = (yp(ip)-ymin)*dyi
  z = (zp(ip)-zmin)*dzi
  ! --- computes velocity
  vx = uxp(ip)*gaminv(ip)
  vy = uyp(ip)*gaminv(ip)
  vz = uzp(ip)*gaminv(ip)
  ! --- computes old position in grid units
  xold=x-dtsdx0*vx
  yold=y-dtsdy0*vy
  zold=z-dtsdz0*vz
  ! --- computes particles weights
  wq=q*w(ip)
  wqx = wq*invdtdx
  wqy = wq*invdtdy
  wqz = wq*invdtdz
  ! --- finds node of cell containing particles for current positions
  iixp0=nint(x)
  ijxp0=nint(y)
  ikxp0=nint(z)
  ! --- computes distance between particle and node for current positions
  xint=x-iixp0
  yint=y-ijxp0
  zint=z-ikxp0
  ! --- computes coefficients for node centered quantities
  xintsq = xint*xint
  sx0(-1) = 0.5_num*(0.5_num-xint)**2
  sx0( 0) = 0.75_num-xintsq
  sx0( 1) = 0.5_num*(0.5_num+xint)**2
  yintsq = yint*yint
  sy0(-1) = 0.5_num*(0.5_num-yint)**2
  sy0( 0) = 0.75_num-yintsq
  sy0( 1) = 0.5_num*(0.5_num+yint)**2
  zintsq = zint*zint
  sz0(-1) = 0.5_num*(0.5_num-zint)**2
  sz0( 0) = 0.75_num-zintsq
  sz0( 1) = 0.5_num*(0.5_num+zint)**2
  ! --- finds node of cell containing particles for old positions
  iixp=nint(xold)
  ijxp=nint(yold)
  ikxp=nint(zold)
  ! --- computes distance between particle and node for old positions
  xint = xold-iixp
  yint = yold-ijxp
  zint = zold-ikxp
  ! --- computes node separation between old and current positions
  dix = iixp-iixp0
  diy = ijxp-ijxp0
  diz = ikxp-ikxp0
  ! --- zero out coefficients (needed because of different dix and diz for each particle)
  sx=0.0_num;sy=0.0_num;sz=0.0_num
  ! --- computes coefficients for quantities centered between nodes
  xintsq = xint*xint
  sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
  sx( 0+dix) = 0.75_num-xintsq
  sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
  yintsq = yint*yint
  sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
  sy( 0+diy) = 0.75_num-yintsq
  sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
  zintsq = zint*zint
  sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
  sz( 0+diz) = 0.75_num-zintsq
  sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
  ! --- computes coefficients difference
  dsx = sx - sx0
  dsy = sy - sy0
  dsz = sz - sz0

! --- computes min/max positions of current contributions
  ixmin = min(0,dix)-1
  ixmax = max(0,dix)+1
  iymin = min(0,diy)-1
  iymax = max(0,diy)+1
  izmin = min(0,diz)-1
  izmax = max(0,diz)+1
  
  ! --- add current contributions
  DO k=izmin, izmax
    DO j=iymin, iymax
      DO i=ixmin, ixmax
        ic = iixp0+i
        jc = ijxp0+j
        kc = ikxp0+k
        IF(i<ixmax) THEN
          sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
          (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k)) 
          IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k) 
          jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)
        END IF
        IF(j<iymax) THEN
          sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
          (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
          IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
          jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)
        END IF
        IF(k<izmax) THEN
          sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
          (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
          IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
          jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)
        END IF
      END DO
    END DO
  END DO
END DO
DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)
RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_2_2_2

! _________________________________________________________________
SUBROUTINE depose_jxjyjz_esirkepov_3_3_3(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                      dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &                                      
                                      nox,noy,noz,l_particles_weight,l4symtry)
!
!
! Esirkepov current deposition at order 3 in x, y, z (nox=noy=noz=3)
! _________________________________________________________________  
  USE omp_lib
  USE constants
  IMPLICIT NONE
  INTEGER :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
  REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
  REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
  REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
  REAL(num) :: clghtisq,usq,xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz, &
                                      s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,         &
                                      oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                                      dtsdx0,dtsdy0,dtsdz0
  REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
  REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx
  REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
  REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
  INTEGER :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                                      ixmin, ixmax, iymin, iymax, izmin, izmax
  LOGICAL(idp) :: l_particles_weight,l4symtry

! PARAMETER INIT
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
dtsdx0 = dt*dxi
dtsdy0 = dt*dyi
dtsdz0 = dt*dzi
invvol = 1.0_num/(dx*dy*dz)
invdtdx = 1.0_num/(dt*dy*dz)
invdtdy = 1.0_num/(dt*dx*dz)
invdtdz = 1.0_num/(dt*dx*dy)
ALLOCATE(sdx(-2:3,-2:3,-2:3),sdy(-2:3,-2:3,-2:3),sdz(-2:3,-2:3,-2:3))
ALLOCATE(sx(-2:3), sx0(-2:3), dsx(-2:3))
ALLOCATE(sy(-2:3), sy0(-2:3), dsy(-2:3))
ALLOCATE(sz(-2:3), sz0(-2:3), dsz(-2:3))
clghtisq = 1.0_num/clight**2
sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num
dtsdz0 = dt*dzi
DO ip=1,np
  ! --- computes current position in grid units
  x = (xp(ip)-xmin)*dxi
  y = (yp(ip)-ymin)*dyi
  z = (zp(ip)-zmin)*dzi
  ! --- computes velocity
  vx = uxp(ip)*gaminv(ip)
  vy = uyp(ip)*gaminv(ip)
  vz = uzp(ip)*gaminv(ip)
  ! --- computes old position in grid units
  xold=x-dtsdx0*vx
  yold=y-dtsdy0*vy
  zold=z-dtsdz0*vz
  ! --- computes particles weights
  wq=q*w(ip)
  wqx = wq*invdtdx
  wqy = wq*invdtdy
  wqz = wq*invdtdz
  ! --- finds node of cell containing particles for current positions
  iixp0=floor(x)
  ijxp0=floor(y)
  ikxp0=floor(z)
  ! --- computes distance between particle and node for current positions
  xint=x-iixp0
  yint=y-ijxp0
  zint=z-ikxp0
  ! --- computes coefficients for node centered quantities
  oxint = 1.0_num-xint
  xintsq = xint*xint
  oxintsq = oxint*oxint
  sx0(-1) = onesixth*oxintsq*oxint
  sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
  sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
  sx0( 2) = onesixth*xintsq*xint
  oyint = 1.0_num-yint
  yintsq = yint*yint
  oyintsq = oyint*oyint
  sy0(-1) = onesixth*oyintsq*oyint
  sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
  sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
  sy0( 2) = onesixth*yintsq*yint
  ozint = 1.0_num-zint
  zintsq = zint*zint
  ozintsq = ozint*ozint
  sz0(-1) = onesixth*ozintsq*ozint
  sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
  sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
  sz0( 2) = onesixth*zintsq*zint
  ! --- finds node of cell containing particles for old positions
  iixp=floor(xold)
  ijxp=floor(yold)
  ikxp=floor(zold)
  ! --- computes distance between particle and node for old positions
  xint = xold-iixp
  yint = yold-ijxp
  zint = zold-ikxp
  ! --- computes node separation between old and current positions
  dix = iixp-iixp0
  diy = ijxp-ijxp0
  diz = ikxp-ikxp0
  ! --- zero out coefficients (needed because of different dix and diz for each particle)
  sx=0.0_num;sy=0.0_num;sz=0.0_num
  ! --- computes coefficients for quantities centered between nodes
  oxint = 1.0_num-xint
  xintsq = xint*xint
  oxintsq = oxint*oxint
  sx(-1+dix) = onesixth*oxintsq*oxint
  sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
  sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
  sx( 2+dix) = onesixth*xintsq*xint
  oyint = 1.0_num-yint
  yintsq = yint*yint
  oyintsq = oyint*oyint
  sy(-1+diy) = onesixth*oyintsq*oyint
  sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
  sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
  sy( 2+diy) = onesixth*yintsq*yint
  ozint = 1.0_num-zint
  zintsq = zint*zint
  ozintsq = ozint*ozint
  sz(-1+diz) = onesixth*ozintsq*ozint
  sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
  sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
  sz( 2+diz) = onesixth*zintsq*zint
  ! --- computes coefficients difference
  dsx = sx - sx0
  dsy = sy - sy0
  dsz = sz - sz0

  ! --- computes min/max positions of current contributions
  ixmin = min(0,dix)-1
  ixmax = max(0,dix)+2
  iymin = min(0,diy)-1
  iymax = max(0,diy)+2
  izmin = min(0,diz)-1
  izmax = max(0,diz)+2
  
  ! --- add current contributions
  DO k=izmin, izmax
    DO j=iymin, iymax
      DO i=ixmin, ixmax
        ic = iixp0+i
        jc = ijxp0+j
        kc = ikxp0+k
        IF(i<ixmax) THEN
          sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
          (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k)) 
          IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k) 
          jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)
        END IF
        IF(j<iymax) THEN
          sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
          (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
          IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
          jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)
        END IF
        IF(k<izmax) THEN
          sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
          (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
          IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
          jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)
        END IF
      END DO
    END DO
  END DO
  
END DO
DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)
RETURN
END SUBROUTINE depose_jxjyjz_esirkepov_3_3_3


! _________________________________________________________________ 
SUBROUTINE depose_jxjyjz_esirkepov_vecHV_1_1_1(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                      dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                      nox,noy,noz,l_particles_weight,l4symtry)
!
!
! Esirkepov current deposition optimized at order 1
! _________________________________________________________________  
USE constants
USE precomputed
IMPLICIT NONE
    INTEGER(idp) :: np,nx,ny,nz,nxguard,nyguard,nzguard, nox, noy,noz
    REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: jxcells,jycells,jzcells
    REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
    REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
    REAL(num) :: xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq, oxintsq,oyintsq, ozintsq
    REAL(num) :: x,y,z,xmid,ymid,zmid
    REAL(num) ::   ww, wwx, wwy, wwz
    REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
    INTEGER(idp) :: j,k,l,j0,k0,l0,ip, NCELLS, ic, ix, iy, iz
    INTEGER(idp) :: nnx, nnxy,ngridx, ngridy, n,nn,nv
    INTEGER(idp) :: moff(1:8)

    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL
    REAL(num), DIMENSION(LVEC) :: vx,vy,vz
    REAL(num) ::  wwwx(LVEC,16), wwwy(LVEC,16),wwwz(LVEC,16), wq
    REAL(num) :: sx1(LVEC),sx2(LVEC),sx3(LVEC),sx4(LVEC)
    REAL(num) :: sx01(LVEC),sx02(LVEC),sx03(LVEC),sx04(LVEC)
    REAL(num) :: sy1,sy2,sy3,sy4,sz1,sz2,sz3,sz4
    REAL(num) :: sy01,sy02,sy03,sy04,sz01,sz02,sz03,sz04
    REAL(num), DIMENSION(4) :: szz, zdec, h1, h11, h12, sgn
    REAL(num):: wwwx1(LVEC,8),wwwx2(LVEC,8),wwwy1(LVEC,8),wwwy2(LVEC,8),wwwz1(LVEC,8),wwwz2(LVEC,8)
    REAL(num):: wx1,wx2,wy1,wy2,wz1,wz2
    INTEGER(idp) :: orig, ncxy, ncx, ncy, ncz, ngx, ngxy, igrid, jorig, korig, lorig
    LOGICAL(idp) :: l_particles_weight,l4symtry ! Useless here but need to be passed in argument to match func_order arguments
    REAL(num), DIMENSION(:), ALLOCATABLE:: sx, sx0, dsx
    REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
    REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
    REAL(num) ::  dtsdx,dtsdy,dtsdz
    INTEGER(idp) :: dix,diy,diz
    REAL(num) ::dtsdx0,dtsdy0,dtsdz0
    INTEGER(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp
    INTEGER(idp) :: iixporig,ijxporig,ikxporig
    REAL(num) :: invdtdx,invdtdy,invdtdz
    REAL(num) :: wqx,wqy,wqz
    REAL(num) :: xold,yold,zold
  

    ngridx=nx+1+2*nxguard
    ngridy=ny+1+2*nyguard
    ncx=nx+4
    ncy=ny+4
    ncz=nz+4
    NCELLS=ncx*ncy*ncz
    ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
    jxcells=0.0_num
    jycells=0.0_num
    jzcells=0.0_num
    nnx = ngridx
    nnxy = ngridx*ngridy
    moff = (/-nnxy,0_idp,nnxy,2_idp*nnxy,nnx-nnxy,nnx,nnx+nnxy,nnx+2_idp*nnxy/)
    iixporig=-2_idp
    ijxporig=-2_idp
    ikxporig=-2_idp
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ncxy=ncx*ncy
    ngxy=(ngridx*ngridy-ncxy)

    h1=(/1_num,0_num,1_num,0_num/); sgn=(/1_num,-1_num,1_num,-1_num/)
    h11=(/0_num,1_num,1_num,0_num/); h12=(/1_num,0_num,0_num,1_num/)
    ! LOOP ON PARTICLES
    DO ip=1,np, LVEC
        !DIR$ ASSUME_ALIGNED xp:64,yp:64,zp:64
        !DIR$ ASSUME_ALIGNED vx:64,vy:64,vz:64
        !DIR$ ASSUME_ALIGNED sx1:64,sx2:64,sx3:64,sx4:64
        !DIR$ ASSUME_ALIGNED sx01:64,sx02:64,sx03:64,sx04:64
        !DIR$ ASSUME_ALIGNED ICELL:64
        !$OMP SIMD
        DO n=1,MIN(LVEC,np-ip+1)
            nn=ip+n-1
            ! --- computes position in  grid units at (n+1)
            x = (xp(nn)-xmin)*dxi
            y = (yp(nn)-ymin)*dyi
            z = (zp(nn)-zmin)*dzi

            ! --- Computes velocity
            vx(n) = uxp(nn)*gaminv(nn)
            vy(n) = uyp(nn)*gaminv(nn)
            vz(n) = uzp(nn)*gaminv(nn)

            ! --- computes old position in grid units
						xold=x-dtsdx0*vx(n)
						yold=y-dtsdy0*vy(n)
						zold=z-dtsdz0*vz(n)
						
						! --- computes particles weights
						wq=q*w(ip)
						wqx = wq*invdtdx
						wqy = wq*invdtdy
						wqz = wq*invdtdz

						! --- finds node of cell containing particles for current positions
						iixp0=floor(x)
						ijxp0=floor(y)
						ikxp0=floor(z)
						! --- computes distance between particle and node for current positions
						xint=x-iixp0
						yint=y-ijxp0
						zint=z-ikxp0

						! --- computes coefficients for node centered quantities
						sx0(0) = 1.0_num-xint
						sx0(1) = xint
						sy0(0) = 1.0_num-yint
						sy0(1) = yint
						sz0(0) = 1.0_num-zint
						sz0(1) = zint

						! --- finds node of cell containing particles for old positions
						iixp=floor(xold)
						ijxp=floor(yold)
						ikxp=floor(zold)
  
						! --- computes distance between particle and node for old positions
						xint = xold-iixp
						yint = yold-ijxp
						zint = zold-ikxp

            ! --- computes node separation between old and current positions
            dix = iixp-iixp0
            diy = ijxp-ijxp0
            diz = ikxp-ikxp0
  
            ! --- zero out coefficients (needed because of different dix and diz for each particle)
            sx=0.0_num;sy=0.0_num;sz=0.0_num
  
            ! --- computes coefficients for quantities centered between nodes
            sx( 0+dix) = 1.0_num-xint
            sx( 1+dix) = xint
            sy( 0+diy) = 1.0_num-yint
            sy( 1+diy) = yint
            sz( 0+diz) = 1.0_num-zint
            sz( 1+diz) = zint

            ! --- computes coefficients difference
            dsx = sx - sx0
            dsy = sy - sy0
            dsz = sz - sz0  
  
            ICELL(n,1)=1+(iixp0-iixporig)+(ijxp0-ijxporig)*ncx+(ikxp0-ikxporig)*ncxy
            ICELL(n,2)=1+(j-jorig)+(k0-korig)*ncx+(l-lorig)*ncxy
            ICELL(n,3)=1+(j-jorig)+(k-korig)*ncx+(l0-lorig)*ncxy

        END DO
        !$OMP END SIMD

        ! Add weights to nearest vertices
        DO n=1,MIN(LVEC,np-ip+1)
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
            !$OMP SIMD
            DO nv=1,8
                ! --- JX
                wx1=wwwx1(n,nv); wx2=wwwx2(n,nv)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)-ncx-1) = jxcells(nv,ICELL(n,1)-ncx-1) + wx1*sx01(n)*vx(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)-ncx)   = jxcells(nv,ICELL(n,1)-ncx)   + wx1*sx02(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+1) = jxcells(nv,ICELL(n,1)-ncx+1) + wx1*sx03(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)-ncx+2) = jxcells(nv,ICELL(n,1)-ncx+2) + wx1*sx04(n)*vx(n)
                ! Loop on (i=-1,j,k)
                jxcells(nv,ICELL(n,1)+ncx-1) = jxcells(nv,ICELL(n,1)+ncx-1) + wx2*sx01(n)*vx(n)
                ! Loop on (i=0,j,k)
                jxcells(nv,ICELL(n,1)+ncx)   = jxcells(nv,ICELL(n,1)+ncx)   + wx2*sx02(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+1) = jxcells(nv,ICELL(n,1)+ncx+1) + wx2*sx03(n)*vx(n)
                !Loop on (i=1,j,k)
                jxcells(nv,ICELL(n,1)+ncx+2) = jxcells(nv,ICELL(n,1)+ncx+2) + wx2*sx04(n)*vx(n)

                ! --- JY
                wy1=wwwy1(n,nv); wy2=wwwy2(n,nv)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)-ncx-1) = jycells(nv,ICELL(n,2)-ncx-1) + wy1*sx1(n)*vy(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)-ncx)   = jycells(nv,ICELL(n,2)-ncx)   + wy1*sx2(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+1) = jycells(nv,ICELL(n,2)-ncx+1) + wy1*sx3(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)-ncx+2) = jycells(nv,ICELL(n,2)-ncx+2) + wy1*sx4(n)*vy(n)
                ! Loop on (i=-1,j,k)
                jycells(nv,ICELL(n,2)+ncx-1) = jycells(nv,ICELL(n,2)+ncx-1) + wy2*sx1(n)*vy(n)
                ! Loop on (i=0,j,k)
                jycells(nv,ICELL(n,2)+ncx)   = jycells(nv,ICELL(n,2)+ncx)   + wy2*sx2(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+1) = jycells(nv,ICELL(n,2)+ncx+1) + wy2*sx3(n)*vy(n)
                !Loop on (i=1,j,k)
                jycells(nv,ICELL(n,2)+ncx+2) = jycells(nv,ICELL(n,2)+ncx+2) + wy2*sx4(n)*vy(n)

                ! --- JZ
                wz1=wwwz1(n,nv); wz2=wwwz2(n,nv)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)-ncx-1) = jzcells(nv,ICELL(n,3)-ncx-1) + wz1*sx1(n)*vz(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)-ncx)   = jzcells(nv,ICELL(n,3)-ncx)   + wz1*sx2(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+1) = jzcells(nv,ICELL(n,3)-ncx+1) + wz1*sx3(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)-ncx+2) = jzcells(nv,ICELL(n,3)-ncx+2) + wz1*sx4(n)*vz(n)
                ! Loop on (i=-1,j,k)
                jzcells(nv,ICELL(n,3)+ncx-1) = jzcells(nv,ICELL(n,3)+ncx-1) + wz2*sx1(n)*vz(n)
                ! Loop on (i=0,j,k)
                jzcells(nv,ICELL(n,3)+ncx)   = jzcells(nv,ICELL(n,3)+ncx)   + wz2*sx2(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+1) = jzcells(nv,ICELL(n,3)+ncx+1) + wz2*sx3(n)*vz(n)
                !Loop on (i=1,j,k)
                jzcells(nv,ICELL(n,3)+ncx+2) = jzcells(nv,ICELL(n,3)+ncx+2) + wz2*sx4(n)*vz(n)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
    DO iz=1, ncz
        DO iy=1,ncy
            !$OMP SIMD
            DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                ! jx
                jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1,ic)
                jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2,ic)
                jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3,ic)
                jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4,ic)
                jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5,ic)
                jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6,ic)
                jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7,ic)
                jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8,ic)
                ! jy
                jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1,ic)
                jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2,ic)
                jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3,ic)
                jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4,ic)
                jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5,ic)
                jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6,ic)
                jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7,ic)
                jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8,ic)
                ! jz
                jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1,ic)
                jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2,ic)
                jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3,ic)
                jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4,ic)
                jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5,ic)
                jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6,ic)
                jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7,ic)
                jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8,ic)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    DEALLOCATE(jxcells,jycells,jzcells)
    RETURN


END SUBROUTINE

!===========================================================================================
! ! Esirkepov current deposition algorithm for linear, quadratic or cubic splines
! WARNING: Highly unoptimized routine ---> USE INLINED ROUTINE
SUBROUTINE pxr_depose_jxjyjz_esirkepov_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
nox,noy,noz,l_particles_weight,l4symtry)
!===========================================================================================

USE constants
IMPLICIT NONE
INTEGER(idp) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
REAL(num), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
REAL(num), DIMENSION(np) :: xp,yp,zp,uxp,uyp,uzp, w, gaminv
REAL(num) :: q,dt,dx,dy,dz,xmin,ymin,zmin
REAL(num) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: sdx,sdy,sdz
REAL(num) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
REAL(num), PARAMETER :: onesixth=1.0_num/6.0_num,twothird=2.0_num/3.0_num
REAL(num), DIMENSION(:), ALLOCATABLE :: sx, sx0, dsx
REAL(num), DIMENSION(:), ALLOCATABLE :: sy, sy0, dsy
REAL(num), DIMENSION(:), ALLOCATABLE :: sz, sz0, dsz
INTEGER(idp) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx, ndtody, ndtodz, &
xl,xu,yl,yu,zl,zu
LOGICAL(idp) :: l_particles_weight,l4symtry

! PARAMETER INIT
ndtodx = int(clight*dt/dx)
ndtody = int(clight*dt/dy)
ndtodz = int(clight*dt/dz)
dxi = 1.0_num/dx
dyi = 1.0_num/dy
dzi = 1.0_num/dz
dtsdx0 = dt*dxi
dtsdy0 = dt*dyi
dtsdz0 = dt*dzi
dts2dx0 = 0.5_num*dtsdx0
dts2dy0 = 0.5_num*dtsdy0
dts2dz0 = 0.5_num*dtsdz0
invvol = 1.0_num/(dx*dy*dz)
invdtdx = 1.0_num/(dt*dy*dz)
invdtdy = 1.0_num/(dt*dx*dz)
invdtdz = 1.0_num/(dt*dx*dy)

xl = -int(nox/2)-1-ndtodx
xu = int((nox+1)/2)+1+ndtodx
yl = -int(noy/2)-1-ndtody
yu = int((noy+1)/2)+1+ndtody
zl = -int(noz/2)-1-ndtodz
zu = int((noz+1)/2)+1+ndtodz


ALLOCATE(sdx(xl:xu,yl:yu,zl:zu),sdy(xl:xu,yl:yu,zl:zu),sdz(xl:xu,yl:yu,zl:zu))
ALLOCATE(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
ALLOCATE(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
ALLOCATE(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

sx0=0.0_num;sy0=0.0_num;sz0=0.0_num
sdx=0.0_num;sdy=0.0_num;sdz=0.0_num


DO ip=1,np
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    y = (yp(ip)-ymin)*dyi
    z = (zp(ip)-zmin)*dzi
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    yold=y-dtsdy0*vy
    zold=z-dtsdz0*vz
    ! --- applies 4-fold symmetry
    IF (l4symtry) THEN
        x=abs(x)
        y=abs(y)
        xold=abs(xold)
        yold=abs(yold)
        vx = (x-xold)/dtsdx0
        vy = (y-yold)/dtsdy0
    END IF
    !computes maximum number of cells traversed by particle in a given dimension
    ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))

    dtsdx = dtsdx0/ncells
    dtsdy = dtsdy0/ncells
    dtsdz = dtsdz0/ncells
    dts2dx = dts2dx0/ncells
    dts2dy = dts2dy0/ncells
    dts2dz = dts2dz0/ncells

    x=xold
    y=yold
    z=zold
    DO icell = 1,ncells
        xold = x
        yold = y
        zold = z
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz

        ! --- computes particles "weights"
        IF (l_particles_weight) THEN
            wq=q*w(ip)
        ELSE
            wq=q*w(1)
        END IF
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

! --- finds node of cell containing particles for current positions
! --- (different for odd/even spline orders)
        IF (nox==2*(nox/2)) THEN
            iixp0=nint(x)
        ELSE
            iixp0=floor(x)
        END IF
        IF (noy==2*(noy/2)) THEN
            ijxp0=nint(y)
        ELSE
            ijxp0=floor(y)
        END IF
        IF (noz==2*(noz/2)) THEN
            ikxp0=nint(z)
        ELSE
            ikxp0=floor(z)
        END IF
! --- computes distance between particle and node for current positions
        xint=x-iixp0
        yint=y-ijxp0
        zint=z-ikxp0

! --- computes coefficients for node centered quantities
        SELECT CASE(nox)
        CASE(0)
            sx0( 0) = 1.0_num
        CASE(1)
            sx0( 0) = 1.0_num-xint
            sx0( 1) = xint
        CASE(2)
            xintsq = xint*xint
            sx0(-1) = 0.5_num*(0.5_num-xint)**2
            sx0( 0) = 0.75_num-xintsq
            sx0( 1) = 0.5_num*(0.5_num+xint)**2
        CASE(3)
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx0(-1) = onesixth*oxintsq*oxint
            sx0( 0) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx0( 1) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx0( 2) = onesixth*xintsq*xint
        END SELECT

        SELECT CASE(noy)
        CASE(0)
            sy0( 0) = 1.0_num
        CASE(1)
            sy0( 0) = 1.0_num-yint
            sy0( 1) = yint
        CASE(2)
            yintsq = yint*yint
            sy0(-1) = 0.5_num*(0.5_num-yint)**2
            sy0( 0) = 0.75_num-yintsq
            sy0( 1) = 0.5_num*(0.5_num+yint)**2
        CASE(3)
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy0(-1) = onesixth*oyintsq*oyint
            sy0( 0) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy0( 1) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy0( 2) = onesixth*yintsq*yint
        END SELECT

        SELECT CASE(noz)
        CASE(0)
            sz0( 0) = 1.0_num
        CASE(1)
            sz0( 0) = 1.0_num-zint
            sz0( 1) = zint
        CASE(2)
            zintsq = zint*zint
            sz0(-1) = 0.5_num*(0.5_num-zint)**2
            sz0( 0) = 0.75_num-zintsq
            sz0( 1) = 0.5_num*(0.5_num+zint)**2
        CASE(3)
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz0(-1) = onesixth*ozintsq*ozint
            sz0( 0) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz0( 1) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz0( 2) = onesixth*zintsq*zint
        END SELECT

! --- finds node of cell containing particles for old positions
! --- (different for odd/even spline orders)
        IF (nox==2*(nox/2)) THEN
            iixp=nint(xold)
        ELSE
            iixp=floor(xold)
        END IF
        IF (noy==2*(noy/2)) THEN
            ijxp=nint(yold)
        ELSE
            ijxp=floor(yold)
        END IF
        IF (noz==2*(noz/2)) THEN
            ikxp=nint(zold)
        ELSE
            ikxp=floor(zold)
        END IF

! --- computes distance between particle and node for old positions
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.0_num;sy=0.0_num;sz=0.0_num

! --- computes coefficients for quantities centered between nodes
        SELECT CASE(nox)
        CASE(0)
            sx( 0+dix) = 1.0_num
        CASE(1)
            sx( 0+dix) = 1.0_num-xint
            sx( 1+dix) = xint
        CASE(2)
            xintsq = xint*xint
            sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
            sx( 0+dix) = 0.75_num-xintsq
            sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
        CASE(3)
            oxint = 1.0_num-xint
            xintsq = xint*xint
            oxintsq = oxint*oxint
            sx(-1+dix) = onesixth*oxintsq*oxint
            sx( 0+dix) = twothird-xintsq*(1.0_num-xint/2.0_num)
            sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint/2.0_num)
            sx( 2+dix) = onesixth*xintsq*xint
        END SELECT

        SELECT CASE(noy)
        CASE(0)
            sy( 0+diy) = 1.0_num
        CASE(1)
            sy( 0+diy) = 1.0_num-yint
            sy( 1+diy) = yint
        CASE(2)
            yintsq = yint*yint
            sy(-1+diy) = 0.5_num*(0.5_num-yint)**2
            sy( 0+diy) = 0.75_num-yintsq
            sy( 1+diy) = 0.5_num*(0.5_num+yint)**2
        CASE(3)
            oyint = 1.0_num-yint
            yintsq = yint*yint
            oyintsq = oyint*oyint
            sy(-1+diy) = onesixth*oyintsq*oyint
            sy( 0+diy) = twothird-yintsq*(1.0_num-yint/2.0_num)
            sy( 1+diy) = twothird-oyintsq*(1.0_num-oyint/2.0_num)
            sy( 2+diy) = onesixth*yintsq*yint
        END SELECT

        SELECT CASE(noz)
        CASE(0)
            sz( 0+diz) = 1.0_num
        CASE(1)
            sz( 0+diz) = 1.0_num-zint
            sz( 1+diz) = zint
        CASE(2)
            zintsq = zint*zint
            sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
            sz( 0+diz) = 0.75_num-zintsq
            sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
        CASE(3)
            ozint = 1.0_num-zint
            zintsq = zint*zint
            ozintsq = ozint*ozint
            sz(-1+diz) = onesixth*ozintsq*ozint
            sz( 0+diz) = twothird-zintsq*(1.0_num-zint/2.0_num)
            sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint/2.0_num)
            sz( 2+diz) = onesixth*zintsq*zint
        END SELECT

! --- computes coefficients difference
        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0

! --- computes min/max positions of current contributions
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)
! --- add current contributions
        DO k=izmin, izmax
            DO j=iymin, iymax
                DO i=ixmin, ixmax
                    ic = iixp0+i
                    jc = ijxp0+j
                    kc = ikxp0+k
                    IF(i<ixmax) THEN
                        sdx(i,j,k)  = wqx*dsx(i)*((sy0(j)+0.5_num*dsy(j))*sz0(k) + &
                        (0.5_num*sy0(j)+1.0_num/3.0_num*dsy(j))*dsz(k))
                        IF (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                        jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)
                    END IF
                    IF(j<iymax) THEN
                        sdy(i,j,k)  = wqy*dsy(j)*((sz0(k)+0.5_num*dsz(k))*sx0(i) + &
                        (0.5_num*sz0(k)+1.0_num/3.0_num*dsz(k))*dsx(i))
                        IF (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                        jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)
                    END IF
                    IF(k<izmax) THEN
                        sdz(i,j,k)  = wqz*dsz(k)*((sx0(i)+0.5_num*dsx(i))*sy0(j) + &
                        (0.5_num*sx0(i)+1.0_num/3.0_num*dsx(i))*dsy(j))
                        IF (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                        jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)
                    END IF
                END DO
            END DO
        END DO
    END DO
END DO

DEALLOCATE(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)

RETURN
END SUBROUTINE pxr_depose_jxjyjz_esirkepov_n


subroutine warp_depose_jxjyjz_esirkepov_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 nox,noy,noz,l_particles_weight,l4symtry)
! ===========================================
! warp_depose_jxjyjz_esirkepov_n
! 
! Warp fonction for esirkepov
! =========================================== 
   use constants
   implicit none
   integer(8) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jx
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jy
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: jz
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,v,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin, clghtisq,usq,gaminv
   logical(8) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(:,:,:),allocatable :: sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                   dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(:),allocatable :: sx, sx0, dsx
   real(kind=8), DIMENSION(:),allocatable :: sy, sy0, dsy
   real(kind=8), DIMENSION(:),allocatable :: sz, sz0, dsz
   integer(8) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx, ndtody, ndtodz, &
                   xl,xu,yl,yu,zl,zu
    PRINT *, "warp depose_jxjyjz_esirkepov_n in PXR: np,l_particles_wight,l4symtry", np,l_particles_weight,l4symtry
    ndtodx = int(clight*dt/dx)
    ndtody = int(clight*dt/dy)
    ndtodz = int(clight*dt/dz)
    xl = -int(nox/2)-1-ndtodx
    xu = int((nox+1)/2)+1+ndtodx
    yl = -int(noy/2)-1-ndtody
    yu = int((noy+1)/2)+1+ndtody
    zl = -int(noz/2)-1-ndtodz
    zu = int((noz+1)/2)+1+ndtodz
    allocate(sdx(xl:xu,yl:yu,zl:zu),sdy(xl:xu,yl:yu,zl:zu),sdz(xl:xu,yl:yu,zl:zu))
    allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
    allocate(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
    allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))
	clghtisq = 1.0_num/clight**2
    sx0=0.;sy0=0.;sz0=0.
    sdx=0.;sdy=0.;sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx0 = dt*dxi
      dtsdy0 = dt*dyi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dy0 = 0.5*dtsdy0
      dts2dz0 = 0.5*dtsdz0
      invvol = 1./(dx*dy*dz)
      invdtdx = 1./(dt*dy*dz)
      invdtdy = 1./(dt*dx*dz)
      invdtdz = 1./(dt*dx*dy)

      do ip=1,np
      
        ! --- computes current position in grid units
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi
        usq = (uxp(ip)**2 + uyp(ip)**2+uzp(ip)**2)*clghtisq
        gaminv = 1.0_num/sqrt(1.0_num + usq)
        ! --- computes velocity
        vx = uxp(ip)*gaminv
        vy = uyp(ip)*gaminv
        vz = uzp(ip)*gaminv
        
        ! --- computes old position in grid units
        xold=x-dtsdx0*vx
        yold=y-dtsdy0*vy
        zold=z-dtsdz0*vz
 
        ! --- applies 4-fold symmetry
        if (l4symtry) then
          x=abs(x)
          y=abs(y)
          xold=abs(xold)
          yold=abs(yold)
          vx = (x-xold)/dtsdx0
          vy = (y-yold)/dtsdy0
        end if
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdy = dtsdy0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dy = dts2dy0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        y=yold
        z=zold
        
        do icell = 1,ncells

        xold = x
        yold = y
        zold = z
        
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz
        
        ! --- computes particles "weights"
        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q*w(1)
        end if
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

        ! --- finds node of cell containing particles for current positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
          iixp0=nint(x)
        else
          iixp0=floor(x)
        end if
        if (noy==2*(noy/2)) then
          ijxp0=nint(y)
        else
          ijxp0=floor(y)
        end if
        if (noz==2*(noz/2)) then
          ikxp0=nint(z)
        else
          ikxp0=floor(z)
        end if

        ! --- computes distance between particle and node for current positions
        xint=x-iixp0
        yint=y-ijxp0
        zint=z-ikxp0

        ! --- computes coefficients for node centered quantities
        select case(nox)
         case(0)
          sx0( 0) = 1.
         case(1)
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         case(2)
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy0( 0) = 1.
         case(1)
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         case(2)
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz0( 0) = 1.
         case(1)
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         case(2)
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
        end select        

        ! --- finds node of cell containing particles for old positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
          iixp=nint(xold)
        else
          iixp=floor(xold)
        end if
        if (noy==2*(noy/2)) then
          ijxp=nint(yold)
        else
          ijxp=floor(yold)
        end if
        if (noz==2*(noz/2)) then
          ikxp=nint(zold)
        else
          ikxp=floor(zold)
        end if

        ! --- computes distance between particle and node for old positions
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

        ! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

        ! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.;sy=0.;sz=0.

        ! --- computes coefficients for quantities centered between nodes
        select case(nox)
         case(0)
          sx( 0+dix) = 1.
         case(1)
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
         case(2)
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy( 0+diy) = 1.
         case(1)
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
         case(2)
          yintsq = yint*yint
          sy(-1+diy) = 0.5*(0.5-yint)**2
          sy( 0+diy) = 0.75-yintsq
          sy( 1+diy) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1+diy) = onesixth*oyintsq*oyint
          sy( 0+diy) = twothird-yintsq*(1.-yint/2)
          sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
          sy( 2+diy) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz( 0+diz) = 1.
         case(1)
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint
         case(2)
          zintsq = zint*zint
          sz(-1+diz) = 0.5*(0.5-zint)**2
          sz( 0+diz) = 0.75-zintsq
          sz( 1+diz) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1+diz) = onesixth*ozintsq*ozint
          sz( 0+diz) = twothird-zintsq*(1.-zint/2)
          sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
          sz( 2+diz) = onesixth*zintsq*zint
        end select        

        ! --- computes coefficients difference
        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0
        
        ! --- computes min/max positions of current contributions
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)

        ! --- add current contributions
        do k=izmin, izmax
          do j=iymin, iymax
            do i=ixmin, ixmax
              ic = iixp0+i
              jc = ijxp0+j
              kc = ikxp0+k
              if(i<ixmax) then
                sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
                if (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                jx(ic,jc,kc) = jx(ic,jc,kc) + sdx(i,j,k)
              end if
              if(j<iymax) then
                sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
                if (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                jy(ic,jc,kc) = jy(ic,jc,kc) + sdy(i,j,k)
              end if
              if(k<izmax) then
                sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
                if (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                jz(ic,jc,kc) = jz(ic,jc,kc) + sdz(i,j,k)
              end if
            end do        
          end do        
        end do        

      end do
 
    end do

    deallocate(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)

  return
end subroutine warp_depose_jxjyjz_esirkepov_n

! ======================================================
subroutine picsar_depose_jxjyjz_esirkepov_n(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 nox,noy,noz,l_particles_weight,l4symtry)
   use constants
   implicit none
   integer(8) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(8) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(:,:,:),allocatable :: sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                   dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(:),allocatable :: sx, sx0, dsx
   real(kind=8), DIMENSION(:),allocatable :: sy, sy0, dsy
   real(kind=8), DIMENSION(:),allocatable :: sz, sz0, dsz
   integer(8) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ndtodx, ndtody, ndtodz, &
                   xl,xu,yl,yu,zl,zu
    PRINT *, "pxr depose_jxjyjz_esirkepov_n: np,l_particles_wight,l4symtry", np,l_particles_weight,l4symtry
    ndtodx = int(clight*dt/dx)
    ndtody = int(clight*dt/dy)
    ndtodz = int(clight*dt/dz)
    xl = -int(nox/2)-1-ndtodx
    xu = int((nox+1)/2)+1+ndtodx
    yl = -int(noy/2)-1-ndtody
    yu = int((noy+1)/2)+1+ndtody
    zl = -int(noz/2)-1-ndtodz
    zu = int((noz+1)/2)+1+ndtodz
    allocate(sdx(xl:xu,yl:yu,zl:zu),sdy(xl:xu,yl:yu,zl:zu),sdz(xl:xu,yl:yu,zl:zu))
    allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
    allocate(sy(yl:yu), sy0(yl:yu), dsy(yl:yu))
    allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

    sx0=0.;sy0=0.;sz0=0.
    sdx=0.;sdy=0.;sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx0 = dt*dxi
      dtsdy0 = dt*dyi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dy0 = 0.5*dtsdy0
      dts2dz0 = 0.5*dtsdz0
      invvol = 1./(dx*dy*dz)
      invdtdx = 1./(dt*dy*dz)
      invdtdy = 1./(dt*dx*dz)
      invdtdz = 1./(dt*dx*dy)

      do ip=1,np
      
        ! --- computes current position in grid units
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi
        
        ! --- computes velocity
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)
        
        ! --- computes old position in grid units
        xold=x-dtsdx0*vx
        yold=y-dtsdy0*vy
        zold=z-dtsdz0*vz
 
        ! --- applies 4-fold symmetry
        if (l4symtry) then
          x=abs(x)
          y=abs(y)
          xold=abs(xold)
          yold=abs(yold)
          vx = (x-xold)/dtsdx0
          vy = (y-yold)/dtsdy0
        end if
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1!+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdy = dtsdy0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dy = dts2dy0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        y=yold
        z=zold
        
        do icell = 1,ncells

        xold = x
        yold = y
        zold = z
        
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz
        
        ! --- computes particles "weights"
        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q*w(1)
        end if
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

        ! --- finds node of cell containing particles for current positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
          iixp0=nint(x)
        else
          iixp0=floor(x)
        end if
        if (noy==2*(noy/2)) then
          ijxp0=nint(y)
        else
          ijxp0=floor(y)
        end if
        if (noz==2*(noz/2)) then
          ikxp0=nint(z)
        else
          ikxp0=floor(z)
        end if

        ! --- computes distance between particle and node for current positions
        xint=x-iixp0
        yint=y-ijxp0
        zint=z-ikxp0

        ! --- computes coefficients for node centered quantities
        select case(nox)
         case(0)
          sx0( 0) = 1.
         case(1)
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         case(2)
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy0( 0) = 1.
         case(1)
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         case(2)
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz0( 0) = 1.
         case(1)
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         case(2)
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
        end select        

        ! --- finds node of cell containing particles for old positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
          iixp=nint(xold)
        else
          iixp=floor(xold)
        end if
        if (noy==2*(noy/2)) then
          ijxp=nint(yold)
        else
          ijxp=floor(yold)
        end if
        if (noz==2*(noz/2)) then
          ikxp=nint(zold)
        else
          ikxp=floor(zold)
        end if

        ! --- computes distance between particle and node for old positions
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

        ! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

        ! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.;sy=0.;sz=0.

        ! --- computes coefficients for quantities centered between nodes
        select case(nox)
         case(0)
          sx( 0+dix) = 1.
         case(1)
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
         case(2)
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy( 0+diy) = 1.
         case(1)
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
         case(2)
          yintsq = yint*yint
          sy(-1+diy) = 0.5*(0.5-yint)**2
          sy( 0+diy) = 0.75-yintsq
          sy( 1+diy) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1+diy) = onesixth*oyintsq*oyint
          sy( 0+diy) = twothird-yintsq*(1.-yint/2)
          sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
          sy( 2+diy) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz( 0+diz) = 1.
         case(1)
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint
         case(2)
          zintsq = zint*zint
          sz(-1+diz) = 0.5*(0.5-zint)**2
          sz( 0+diz) = 0.75-zintsq
          sz( 1+diz) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1+diz) = onesixth*ozintsq*ozint
          sz( 0+diz) = twothird-zintsq*(1.-zint/2)
          sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
          sz( 2+diz) = onesixth*zintsq*zint
        end select        

        ! --- computes coefficients difference
        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0
        
        ! --- computes min/max positions of current contributions
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)

        ! --- add current contributions
        do k=izmin, izmax
          do j=iymin, iymax
            do i=ixmin, ixmax
              ic = iixp0+i
              jc = ijxp0+j
              kc = ikxp0+k
              if(i<ixmax) then
                sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
                if (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                cj(ic,jc,kc,1) = cj(ic,jc,kc,1) + sdx(i,j,k)
              end if
              if(j<iymax) then
                sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
                if (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                cj(ic,jc,kc,2) = cj(ic,jc,kc,2) + sdy(i,j,k)
              end if
              if(k<izmax) then
                sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
                if (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                cj(ic,jc,kc,3) = cj(ic,jc,kc,3) + sdz(i,j,k)
              end if
            end do        
          end do        
        end do        

      end do
 
    end do

    deallocate(sdx,sdy,sdz,sx,sx0,dsx,sy,sy0,dsy,sz,sz0,dsz)

  return
end subroutine picsar_depose_jxjyjz_esirkepov_n


! ________________________________________________________________________________________
SUBROUTINE current_reduction(jx,jy,jz,jxcells,jycells,jzcells,ncells,nx,ny,nz,nxguard,nyguard,nzguard,ncx,ncy,ncz)
! This subroutine performs the reduction of jxcellx, jycells and jzcells into jx,jy and jz.
! This subroutine is called after the loop on particles where 
! depose_jxjyjz_vecHV_vnr_1_1_1 isperformed for each species
!
! Inputs:
! - jx,jy,jz
! - jxcells,jycells,jzcells: transient current arrays
! - ncells: tile cell numbers
! - nx,ny,nz: tile cell numbers in each direction
! - nxguard,nyguard,nzguard
! - ncx,ncy,ncz
!
! Outputs:
! - jx,jy,jz updated
! ________________________________________________________________________________________
    USE constants
    USE precomputed
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nx,ny,nz,ncells
    INTEGER(idp), INTENT(IN) :: nxguard,nyguard,nzguard
    REAL(num),INTENT(IN OUT) :: jx(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jy(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN OUT) :: jz(1:(1+nx+2*nxguard)*(1+ny+2*nyguard)*(1+nz+2*nzguard))
    REAL(num),INTENT(IN), DIMENSION(8,ncells):: jxcells,jycells,jzcells


    INTEGER(idp) :: nnx, nnxy, nn,nv
    INTEGER(idp) :: moff(1:8)

    INTEGER(idp), DIMENSION(LVEC,3) :: ICELL
    INTEGER(idp) :: orig, jorig, korig, lorig, igrid,ic
    INTEGER(idp) :: ncx, ncy, ncxy, ncz,ix,iy,iz, ngridx, ngridy, ngx, ngxy


    ! _____________________________________________________________
    ! Parameters
    
    ngridx=nx+1+2*nxguard
    ngridy=ny+1+2*nyguard

    nnx = ngridx
    nnxy = nnx*ngridy
    moff = (/0_idp,1_idp,nnx,nnx+1_idp,nnxy,nnxy+1_idp,nnxy+nnx,nnxy+nnx+1_idp/)

    jorig=-2
    korig=-2
    lorig=-2
    orig=jorig+nxguard+nnx*(korig+nyguard)+(lorig+nzguard)*nnxy
    ngx=(ngridx-ncx)
    ngxy=(ngridx*ngridy-ncx*ncy)
    ncxy=ncx*ncy
    
    ! ____________________________________________________________
    ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
    DO iz=1, ncz
        DO iy=1,ncy
            !$OMP SIMD
            DO ix=1,ncx !! VECTOR (take ncx multiple of vector length)
                ic=ix+(iy-1)*ncx+(iz-1)*ncxy
                igrid=ic+(iy-1)*ngx+(iz-1)*ngxy
                ! jx
                jx(orig+igrid+moff(1))=jx(orig+igrid+moff(1))+jxcells(1,ic)
                jx(orig+igrid+moff(2))=jx(orig+igrid+moff(2))+jxcells(2,ic)
                jx(orig+igrid+moff(3))=jx(orig+igrid+moff(3))+jxcells(3,ic)
                jx(orig+igrid+moff(4))=jx(orig+igrid+moff(4))+jxcells(4,ic)
                jx(orig+igrid+moff(5))=jx(orig+igrid+moff(5))+jxcells(5,ic)
                jx(orig+igrid+moff(6))=jx(orig+igrid+moff(6))+jxcells(6,ic)
                jx(orig+igrid+moff(7))=jx(orig+igrid+moff(7))+jxcells(7,ic)
                jx(orig+igrid+moff(8))=jx(orig+igrid+moff(8))+jxcells(8,ic)
                ! jy
                jy(orig+igrid+moff(1))=jy(orig+igrid+moff(1))+jycells(1,ic)
                jy(orig+igrid+moff(2))=jy(orig+igrid+moff(2))+jycells(2,ic)
                jy(orig+igrid+moff(3))=jy(orig+igrid+moff(3))+jycells(3,ic)
                jy(orig+igrid+moff(4))=jy(orig+igrid+moff(4))+jycells(4,ic)
                jy(orig+igrid+moff(5))=jy(orig+igrid+moff(5))+jycells(5,ic)
                jy(orig+igrid+moff(6))=jy(orig+igrid+moff(6))+jycells(6,ic)
                jy(orig+igrid+moff(7))=jy(orig+igrid+moff(7))+jycells(7,ic)
                jy(orig+igrid+moff(8))=jy(orig+igrid+moff(8))+jycells(8,ic)
                ! jz
                jz(orig+igrid+moff(1))=jz(orig+igrid+moff(1))+jzcells(1,ic)
                jz(orig+igrid+moff(2))=jz(orig+igrid+moff(2))+jzcells(2,ic)
                jz(orig+igrid+moff(3))=jz(orig+igrid+moff(3))+jzcells(3,ic)
                jz(orig+igrid+moff(4))=jz(orig+igrid+moff(4))+jzcells(4,ic)
                jz(orig+igrid+moff(5))=jz(orig+igrid+moff(5))+jzcells(5,ic)
                jz(orig+igrid+moff(6))=jz(orig+igrid+moff(6))+jzcells(6,ic)
                jz(orig+igrid+moff(7))=jz(orig+igrid+moff(7))+jzcells(7,ic)
                jz(orig+igrid+moff(8))=jz(orig+igrid+moff(8))+jzcells(8,ic)
            END DO
            !$OMP END SIMD
        END DO
    END DO
    RETURN
END SUBROUTINE current_reduction

SUBROUTINE add_pxrjxjyjz_towarp_j(jwarp,cdims,nxx,nyy,nzz,nxg,nyg,nzg)
USE fields
USE shared_data
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nxx, nyy, nzz, nxg, nyg, nzg, cdims
REAL(num), INTENT (IN OUT), DIMENSION(-nxg:nxx+nxg,-nyg:nyy+nyg,-nzg:nzz+nzg,cdims) :: jwarp
INTEGER(idp) :: idim

jwarp(:,:,:,1)=jwarp(:,:,:,1)+jx
jwarp(:,:,:,2)=jwarp(:,:,:,2)+jy
jwarp(:,:,:,3)=jwarp(:,:,:,3)+jz

END SUBROUTINE add_pxrjxjyjz_towarp_j



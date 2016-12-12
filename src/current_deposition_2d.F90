! ________________________________________________________________________________________
! 
! CURRENT_DEPOSITION_2D.F90
!
! List of subroutines
! - pxrdepose_currents_on_grid_jxjyjz_2d
!
! - pxr_depose_jxjyjz_esirkepov2d_n
! - pxr_depose_jxjyjz_esirkepov2d_1_1
! - pxr_depose_jxjyjz_esirkepov2d_2_2
! - pxr_depose_jxjyjz_esirkepov2d_3_3
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Main subroutine for the current deposition called in submain for the 2d
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d
! ________________________________________________________________________________________
  USE fields
  USE shared_data
  USE params
  USE time_stat
#if defined(PROFILING) && PROFILING==2      
  USE ITT_SDE_FORTRAN                       
#endif                                   
  IMPLICIT NONE 
  
  ! __ Parameter declaration __________________________________________________
  REAL(num) :: tdeb, tend
  
  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE

    SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & 
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)!#do not parse
      USE omp_lib
      USE constants
      implicit none
      integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      integer(idp)                          :: lvect      
      real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(inout) :: jx,jy,jz
      real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
      real(num)                             :: q,dt,dx,dz,xmin,zmin
      LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
      real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
    END SUBROUTINE

    SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_2_2(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & 
                                            dt,dx,dz,nx,nz,nxguard,nzguard, &
                                            nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)!#do not parse
      USE omp_lib
      USE constants
      implicit none
      integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      integer(idp)                          :: lvect
      real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(inout) :: jx,jy,jz
      real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
      real(num)                             :: q,dt,dx,dz,xmin,zmin
      LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
      real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
    END SUBROUTINE

    SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & 
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                    nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)!#do not parse
      USE omp_lib
      USE constants
      implicit none
      integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      integer(idp)                          :: lvect
      real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(inout) :: jx,jy,jz
      real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
      real(num)                             :: q,dt,dx,dz,xmin,zmin
      LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
      real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
    END SUBROUTINE

  subroutine pxr_depose_jxjyjz_esirkepov2d_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                   dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                   nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose) !#do not parse
     use constants
     implicit none
     integer(idp)                           :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
     real(num), dimension(np)               :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
     real(num)                              :: q,dt,dx,dz,xmin,zmin
     LOGICAL(lp)                            :: l_particles_weight,l4symtry,l_2drz     
     real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  END SUBROUTINE

#if defined(DEV)
  subroutine pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,&
                                                   gaminv,w,q,xmin,zmin, &
                                                   dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                   nox,noz,lvect,l_particles_weight,&
                                                   l4symtry,l_2drz,type_rz_depose) !#do not parse
    USE omp_lib
    USE constants
    implicit none
    integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
    integer(idp)                          :: lvect
    real(num), dimension((1+nx+2*nxguard)*(1+nz+2*nzguard)), intent(in out) :: jx,jy,jz
    real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
    real(num)                             :: q,dt,dx,dz,xmin,zmin
    LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
  END SUBROUTINE
#endif
  
  END INTERFACE
  ! ___________________________________________________________________________

! For debugging    
#if defined(DEBUG)
  WRITE(0,*) "Depose_currents_on_grid: start"
#endif

  ! For time statistics
  tdeb=MPI_WTIME()

! For profiling with Vtune/SDE
#if PROFILING==2              
  CALL start_collection()     
#endif                        

  jx = 0.0_num
  jy = 0.0_num
  jz = 0.0_num

  ! __ Current deposition ________________________________________________________________

  ! _______________________________________________________
  ! Esirkepov general order subroutine
  IF (currdepo.EQ.2) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
         nox,noy,noz,dx,dy,dz,dt)
  ! _______________________________________________________
  ! Esirkepov OpenMP/tiling version non-vectorized but more optimized than the general order subroutine
  ELSE IF (currdepo.EQ.1) THEN

    ! Order 1
    IF ((nox.eq.1).AND.(noz.eq.1)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_1_1,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
    ! Order 2
    ELSE IF ((nox.eq.2).AND.(noz.eq.2)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_2_2,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
    ! Order 3
    ELSE IF ((nox.eq.3).AND.(noz.eq.3)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_3_3,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
    ! Order n
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
         nox,noy,noz,dx,dy,dz,dt)
    ENDIF

  ! _______________________________________________________
  ! Default - Esirkepov parallel version with OPENMP/tiling and optimizations
  ELSE

    ! Order 1
    IF ((nox.eq.1).AND.(noz.eq.1)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_1_1,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
    ! Order 2
    ELSE IF ((nox.eq.2).AND.(noz.eq.2)) THEN
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_2_2,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
    ! Order 3
    ELSE IF ((nox.eq.3).AND.(noz.eq.3)) THEN 
      CALL pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(pxr_depose_jxjyjz_esirkepov2d_3_3,jx,jy,jz,&
                     nx,ny,nz,nxjguards,nyjguards,nzjguards,nox,noy,noz,dx,dy,dz,dt,lvec_curr_depo)
    ! Order n
    ELSE
      CALL pxrdepose_currents_on_grid_jxjyjz_sub_openmp(jx,jy,jz,nx,ny,nz,nxjguards,nyjguards,nzjguards, &
         nox,noy,noz,dx,dy,dz,dt)
    ENDIF
      
  ENDIF

! Stop Vtune/SDE analysis
#if PROFILING==2
  CALL stop_collection() 
#endif

  ! For time statistics
  tend = MPI_WTIME()
  localtimes(3)=localtimes(3)+(tend-tdeb)
  
! For debugging   
#if defined(DEBUG)
  WRITE(0,*) "Depose_current_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d



! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile with Esirkepov method in 2D
!
!> @details
!> This subroutine is called from Fortran main program and contains an interface argument 
!> OpenMP version. Avoids conflict while reducing tile currents in the global 
!> current array.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(curr_depo_sub,jxg,jyg,jzg,&
                nxx,nyy,nzz,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt,lvect)
! ________________________________________________________________________________________
  USE particles
  USE constants
  USE tiling
  USE omp_lib
  USE timing
  USE time_stat
  IMPLICIT NONE

  ! __ Parameter declaration _______________________________________________________________
  INTEGER(idp), INTENT(IN)  :: nxx,nyy,nzz,nxjguard,nyjguard,nzjguard
  INTEGER(idp), INTENT(IN)  :: noxx,noyy,nozz
  INTEGER(idp), INTENT(IN)  :: lvect
  REAL(num), INTENT(IN)     :: dxx,dyy,dzz, dtt
  REAL(num), INTENT(IN OUT) :: jxg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jyg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  REAL(num), INTENT(IN OUT) :: jzg(-nxjguard:nxx+nxjguard,-nyjguard:nyy+nyjguard,-nzjguard:nzz+nzjguard)
  INTEGER(idp)              :: ispecies, ix, iy, iz, count
  INTEGER(idp)              :: jmin, jmax, kmin, kmax, lmin, lmax
  INTEGER(idp)              :: jminc, jmaxc, kminc, kmaxc, lminc, lmaxc
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER  :: currg
  REAL(num)                 :: tdeb, tend
  INTEGER(idp)              :: nxc, nyc, nzc, nxjg, nyjg, nzjg
  LOGICAL(lp)               :: isdeposited=.FALSE.

  ! ___ Interface _________________________________________________
  ! For the func_order input function
  INTERFACE
    SUBROUTINE curr_depo_sub(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, & !#do not parse
             dt,dx,dz,nx,nz,nxguard,nzguard, & !#do not parse
             nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose) !#do not parse

      USE constants
      IMPLICIT NONE
      INTEGER(idp)             :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
      INTEGER(idp)             :: lvect
      REAL(num), DIMENSION(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
      REAL(num), DIMENSION(np) :: xp,zp,uxp,uyp,uzp, w, gaminv
      REAL(num)                :: q,dt,dx,dz,xmin,zmin
      LOGICAL(lp)              :: l_particles_weight,l4symtry,l_2drz
    END SUBROUTINE

  END INTERFACE
  
  
  tdeb=MPI_WTIME()
  !$OMP PARALLEL DEFAULT(NONE)                                                              &
  !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,nxjguard,nyjguard,              &
  !$OMP nzjguard,dxx,dyy,dzz,dtt,jxg,jyg,jzg,noxx,noyy,nozz,aofgrid_tiles,c_dim,lvect)      &
  !$OMP PRIVATE(ix,iy,iz,ispecies,curr,currg, curr_tile,count,jmin,jmax,kmin,kmax,lmin,     &
  !$OMP lmax,jminc,jmaxc,kminc,kmaxc,lminc,lmaxc,nxc,nyc,nzc, nxjg, nyjg, nzjg, isdeposited)
  !! Current deposition
  !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
  DO iz=1,ntilez
        DO ix=1,ntilex
          curr => species_parray(1)
            curr_tile=>curr%array_of_tiles(ix,1,iz)
            nxjg=curr_tile%nxg_tile
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
            currg=>aofgrid_tiles(ix,1,iz)
            currg%jxtile=0.
            currg%jytile=0.
            currg%jztile=0.!jzg(jmin:jmax,kmin:kmax,lmin:lmax)
            isdeposited=.FALSE.
            DO ispecies=1, nspecies ! LOOP ON SPECIES
                 curr => species_parray(ispecies)
                curr_tile=>curr%array_of_tiles(ix,1,iz)
                count=curr_tile%np_tile(1)
                IF (count .EQ. 0) THEN 
                  CYCLE
                ELSE 
                  isdeposited=.TRUE.
                ENDIF 
 
                ! Depose current in jtile                
                CALL curr_depo_sub(currg%jxtile,currg%jytile,                              &
                currg%jztile,count,                                                       &
                curr_tile%part_x,curr_tile%part_z,                            &
                curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,curr_tile%part_gaminv,           &
                curr_tile%pid(1,wpid),curr%charge,curr_tile%x_grid_tile_min,     &
                curr_tile%z_grid_tile_min,dtt,dxx,dzz,nxc,nzc,                                     &
                nxjg,nzjg,noxx,nozz,lvect,.TRUE._idp,.FALSE._idp,.FALSE._idp,0_idp) 
                
                
            END DO! END LOOP ON SPECIES
            IF (isdeposited) THEN
              jxg(jmin:jmax,0,lmin:lmax)=jxg(jmin:jmax,0,lmin:lmax)+currg%jxtile(0:nxc,0,0:nzc)
              jyg(jmin:jmax,0,lmin:lmax)=jyg(jmin:jmax,0,lmin:lmax)+currg%jytile(0:nxc,0,0:nzc)
              jzg(jmin:jmax,0,lmin:lmax)=jzg(jmin:jmax,0,lmin:lmax)+currg%jztile(0:nxc,0,0:nzc)
            ENDIF
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
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp


! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition esirkepov n order (from 0 to 3)
!
!> @details
!> This subroutine is adapted from the version of WARP.
!> This subroutine is called in pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp().
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx,jy,jz current arrays
!> @param[in] np number of particles
!> @param[in] xp,zp particle position arrays
!> @param[in] uxp,uyp,uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin,zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx,nz tile grid size
!> @param[in] nxguard,nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
!
subroutine pxr_depose_jxjyjz_esirkepov2d_n(jx,jy,jz,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
! ________________________________________________________________________________________
   use constants
   implicit none
   integer(idp) :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
   real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
   real(num), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(num) :: q,dt,dx,dz,xmin,zmin
   LOGICAL(lp)  :: l_particles_weight,l4symtry,l_2drz

   real(num) :: dxi,dzi,dtsdx,dtsdz,xint,yint,zint
   real(num),dimension(:,:), allocatable :: sdx,sdz
   real(num) :: xold,yold,zold,rold,xmid,zmid,x,y,z,r,c,s,wq,wqx,wqz, &
                   tmp,vx,vy,vz,dts2dx,dts2dz, &
                   s1x,s2x,s1z,s2z,invvol,invdtdx,invdtdz, &
                   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq, &
                   dtsdx0,dtsdz0,dts2dx0,dts2dz0
   real(num), parameter :: onesixth=1./6.
   real(num), parameter :: twothird=2./3.
   real(num), dimension(:), allocatable :: sx, sx0, dsx, sz, sz0, dsz
   integer(idp) :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc, &
                   ixmin, ixmax, izmin, izmax, icell, ncells, ndtodx, ndtodz, &
                   xl,xu,zl,zu

    ndtodx = int(clight*dt/dx)
    ndtodz = int(clight*dt/dz)
    xl = -int(nox/2)-1-ndtodx
    xu = int((nox+1)/2)+1+ndtodx
    zl = -int(noz/2)-1-ndtodz
    zu = int((noz+1)/2)+1+ndtodz
    allocate(sdx(xl:xu,zl:zu),sdz(xl:xu,zl:zu))
    allocate(sx(xl:xu), sx0(xl:xu), dsx(xl:xu))
    allocate(sz(zl:zu), sz0(zl:zu), dsz(zl:zu))

    sx0=0.;sz0=0.
    sdx=0.;sdz=0.
      
    ! Davoine method : limited to order 1 in r
    if (type_rz_depose==2) then
       nox = 1
    endif

      dxi = 1./dx
      dzi = 1./dz
      invvol = 1./(dx*dz)
      dtsdx0 = dt*dxi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dz0 = 0.5*dtsdz0
      invdtdx = 1./(dt*dz)
      invdtdz = 1./(dt*dx)

      do ip=1,np
      
        ! --- computes current position in grid units
        x = xp(ip)
        if (l_2drz) then
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-10) then
            c = x/r 
            s = y/r
          else
            c = 1.
            s = 0.
          end if
          x = r
        end if
        x=x*dxi
        z = zp(ip)*dzi
          
        ! --- computes velocity
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)

        ! --- computes old position in grid units
        if (l_2drz) then
          xold = xp(ip)-dt*vx
          yold = yp(ip)-dt*vy
          rold = sqrt(xold*xold+yold*yold)
          xold=rold*dxi
          vy = -vx*s+vy*c
          vx = (x-xold)/dtsdx0
        else
          xold=x-dtsdx0*vx
        end if
        zold=z-dtsdz0*vz
 
        ! --- applies 4-fold symmetry
        if (l4symtry) then
          x=abs(x)
          xold=abs(xold)
          vx = (x-xold)/dtsdx0
        end if

        ! --- sets positions relative to grid  start
        x = x-xmin*dxi
        z = z-zmin*dzi
        xold = xold-xmin*dxi
        zold = zold-zmin*dzi
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1!+max( int(abs(x-xold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        z=zold
        
        do icell = 1,ncells

        xold = x
        zold = z

        x = x+dtsdx*vx
        z = z+dtsdz*vz

        ! --- computes particles "weights"
        if (l_particles_weight) then
           wq=q*w(ip)
        else
           wq=q*w(1)
        end if
        wqx = wq*invdtdx
        wqz = wq*invdtdz

        ! --- finds node of cell containing particles for current positions 
        ! --- (different for odd/even spline orders)
        if (nox==2*(nox/2)) then
          iixp0=nint(x)
        else
          iixp0=floor(x)
        end if
        if (noz==2*(noz/2)) then
          ikxp0=nint(z)
        else
          ikxp0=floor(z)
        end if

        ! --- computes distance between particle and node for current positions
        xint=x-iixp0
        zint=z-ikxp0

        ! --- computes coefficients for node centered quantities
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sx0(0) = 1. - xint  + 1./(4*iixp0+2)*( -xint + xint**2 )
           sx0(1) = 1. - sx0(0)
        else! Standard method, canonical shapes in r 
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
        endif

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
        if (noz==2*(noz/2)) then
           ikxp=nint(zold)
        else
           ikxp=floor(zold)
        end if

        ! --- computes distance between particle and node for old positions
        xint = xold-iixp
        zint = zold-ikxp

        ! --- computes node separation between old and current positions
        dix = iixp-iixp0
        diz = ikxp-ikxp0

        ! --- zero out coefficients (needed because of different dix and diz for each particle)
        sx=0.;sz=0.

        ! --- computes coefficients for quantities centered between nodes
        if (type_rz_depose == 2) then ! Davoine method, modified particle shapes in r
           sx(0) = 1. - xint  + 1./(4*iixp+2)*( -xint + xint**2 )
           sx(1) = 1. - sx(0)
        else! Standard method, canonical shapes in r 
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
        endif
        
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
        dsz = sz - sz0

        ! --- computes min/max positions of current contributions
        ixmin = min(0_idp,dix)-int(nox/2)
        ixmax = max(0_idp,dix)+int((nox+1)/2)
        izmin = min(0_idp,diz)-int(noz/2)
        izmax = max(0_idp,diz)+int((noz+1)/2)
        
        ! --- add current contributions
        ! --- NB : the current is later divided by the cylindrical cell volume in applybc_j
        do k=izmin, izmax
           do i=ixmin, ixmax
              ic = iixp0+i
              kc = ikxp0+k

              ! -- Jx
              if(i<ixmax) then
                 sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
                 if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x
                 jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
              end if
              
              ! -- Jy (2D Esirkepov scheme)
              jy(ic,kc) = jy(ic,kc) + wq*vy*invvol/ncells* &
                   ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )

              ! -- Jz
              if(k<izmax) then
                 sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov
                 if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
                 jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
              end if
           end do
        end do
     end do
    end do
    
    deallocate(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)

  return
end subroutine pxr_depose_jxjyjz_esirkepov2d_n

! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 1
!> This function is not optimized but provides better performances than 
!> using the abitrary order function
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx,jy,jz current arrays
!> @param[in] np number of particles
!> @param[in] xp,zp particle position arrays
!> @param[in] uxp,uyp,uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin,zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx,nz tile grid size
!> @param[in] nxguard,nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
!
SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none

  ! __ Parameter declaration ________________________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  integer(idp)                          :: lvect
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), parameter                  :: onethird=1./3.  
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,idx,idz,i,k,ic,kc
  integer(isp)                          :: dix,diz
  integer(isp)                          :: ixmin, ixmax, izmin, izmax
  integer(idp)                          :: icell, ndtodx, ndtodz
  integer(idp)                          :: xl,xu,zl,zu

  ! __ Parameter initialization ______________________________
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(-1:2,-1:2),sdz(-1:2,-1:2))
  ALLOCATE(sx(-1:2), sx0(-1:2), dsx(-1:2))
  ALLOCATE(sz(-1:2), sz0(-1:2), dsz(-1:2))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num
  
  DO ip=1,np
  
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
  
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
  
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz
  
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz
  
    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ikxp0=floor(z)
  
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0
  
    ! --- computes coefficients for node centered quantities
    sx0( 0) = 1.0_num-xint
    sx0( 1) = xint
  
    sz0( 0) = 1.0_num-zint
    sz0( 1) = zint
  
    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ikxp=floor(zold)
  
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp
  
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0
  
    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
  
    ! --- computes coefficients for quantities centered between nodes
    sx( 0+dix) = 1.0_num-xint
    sx( 1+dix) = xint
  
    sz( 0+diz) = 1.0_num-zint
    sz( 1+diz) = zint
  
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0
    
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)-0
    ixmax = max(0,dix)+1
    izmin = min(0,diz)-0
    izmax = max(0,diz)+1
    
    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k
        
        ! --- Jx
        IF(i<ixmax) THEN
          sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
          if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x 
          jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
        END IF
        
        ! -- Jy (2D Esirkepov scheme)
        jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &
        ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )
        
        ! --- Jz
        IF(k<izmax) THEN
          sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov&
          if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
          jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
        END IF
      END DO
    END DO

  END DO
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)
  RETURN

END SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_1_1

! ________________________________________________________________________________________
! @brief 
!> 2D Current deposition with the method of Esirkepov at order 2
!> This function is not optimized but provides better performances than 
!> using the abitrary order function
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx,jy,jz current arrays
!> @param[in] np number of particles
!> @param[in] xp,zp particle position arrays
!> @param[in] uxp,uyp,uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin,zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx,nz tile grid size
!> @param[in] nxguard,nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
!
SUBROUTINE pxr_depose_jxjyjz_esirkepov2d_2_2(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none
  
  ! __ Parameter declaration ________________________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  integer(idp)                          :: lvect  
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), parameter                  :: onethird=1./3.  
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,idx,idz,i,k,ic,kc
  integer(isp)                          :: dix,diz
  integer(isp)                          :: ixmin, ixmax, izmin, izmax
  integer(idp)                          :: icell, ndtodx, ndtodz
  integer(idp)                          :: xl,xu,zl,zu

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(-2:2,-2:2),sdz(-2:2,-2:2))
  ALLOCATE(sx(-2:2), sx0(-2:2), dsx(-2:2))
  ALLOCATE(sz(-2:2), sz0(-2:2), dsz(-2:2))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num
  
  DO ip=1,np
  
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
  
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
  
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz
  
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz
  
    ! --- finds node of cell containing particles for current positions
    iixp0=nint(x)
    ikxp0=nint(z)
  
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0
  
    ! --- computes coefficients for node centered quantities
    xintsq = xint*xint
    sx0(-1) = 0.5_num*(0.5_num-xint)**2
    sx0( 0) = 0.75_num-xintsq
    sx0( 1) = 0.5_num*(0.5_num+xint)**2
  
    zintsq = zint*zint
    sz0(-1) = 0.5_num*(0.5_num-zint)**2
    sz0( 0) = 0.75_num-zintsq
    sz0( 1) = 0.5_num*(0.5_num+zint)**2
  
    ! --- finds node of cell containing particles for old positions
    iixp=nint(xold)
    ikxp=nint(zold)
  
    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp
  
    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0
  
    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx(-2)=0.0_num
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sz(-2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
  
    ! --- computes coefficients for quantities centered between nodes
    xintsq = xint*xint
    sx(-1+dix) = 0.5_num*(0.5_num-xint)**2
    sx( 0+dix) = 0.75_num-xintsq
    sx( 1+dix) = 0.5_num*(0.5_num+xint)**2
  
    zintsq = zint*zint
    sz(-1+diz) = 0.5_num*(0.5_num-zint)**2
    sz( 0+diz) = 0.75_num-zintsq
    sz( 1+diz) = 0.5_num*(0.5_num+zint)**2
  
    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0
    
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)-1
    ixmax = max(0,dix)+1
    izmin = min(0,diz)-1
    izmax = max(0,diz)+1
    
    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k
        
        ! --- Jx
        IF(i<ixmax) THEN
          sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
          if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x 
          jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
        END IF
        
        ! -- Jy (2D Esirkepov scheme)
        jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &
        ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+onethird*dsz(k))*dsx(i) )
        
        ! --- Jz
        IF(k<izmax) THEN
          sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov&
          if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
          jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
        END IF
      END DO
    END DO

  END DO
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)
  RETURN

End subroutine pxr_depose_jxjyjz_esirkepov2d_2_2

! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 3
!> This function is not optimized but provides better performances than 
!> using the abitrary order function
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx,jy,jz current arrays
!> @param[in] np number of particles
!> @param[in] xp,zp particle position arrays
!> @param[in] uxp,uyp,uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin,zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx,nz tile grid size
!> @param[in] nxguard,nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
!
subroutine pxr_depose_jxjyjz_esirkepov2d_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none
  
  ! __ Parameter declaration _______________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  integer(idp)                          :: lvect
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:), allocatable :: sdx,sdz
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,idx,idz,i,k,ic,kc
  integer(idp)                          :: icell, ndtodx, ndtodz
  integer(isp)                          :: dix,diz
  integer(isp)                          :: ixmin, ixmax, izmin, izmax  
  integer(idp)                          :: xl,xu,zl,zu

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(-2:3,-2:3),sdz(-2:3,-2:3))
  ALLOCATE(sx(-2:3), sx0(-2:3), dsx(-2:3))
  ALLOCATE(sz(-2:3), sz0(-2:3), dsz(-2:3))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num

  DO ip=1,np
  
    ! --- computes current position in grid units
    x = (xp(ip)-xmin)*dxi
    z = (zp(ip)-zmin)*dzi
  
    ! --- computes velocity
    vx = uxp(ip)*gaminv(ip)
    vy = uyp(ip)*gaminv(ip)
    vz = uzp(ip)*gaminv(ip)
  
    ! --- computes old position in grid units
    xold=x-dtsdx0*vx
    zold=z-dtsdz0*vz
  
    ! --- computes particles weights
    wq=q*w(ip)
    wqx = wq*invdtdx
    wqz = wq*invdtdz  

    ! --- finds node of cell containing particles for current positions
    iixp0=floor(x)
    ikxp0=floor(z)
  
    ! --- computes distance between particle and node for current positions
    xint=x-iixp0
    zint=z-ikxp0

    ! --- computes coefficients for node centered quantities
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx0(-1) = onesixth*oxintsq*oxint
    sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx0( 2) = onesixth*xintsq*xint
  
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz0(-1) = onesixth*ozintsq*ozint
    sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz0( 2) = onesixth*zintsq*zint

    ! --- finds node of cell containing particles for old positions
    iixp=floor(xold)
    ikxp=floor(zold)

    ! --- computes distance between particle and node for old positions
    xint = xold-iixp
    zint = zold-ikxp

    ! --- computes node separation between old and current positions
    dix = iixp-iixp0
    diz = ikxp-ikxp0

    ! --- zero out coefficients (needed because of different dix and diz for each particle)
    sx(-2)=0.0_num
    sx(-1)=0.0_num
    sx(0)=0.0_num
    sx(1)=0.0_num
    sx(2)=0.0_num
    sx(3)=0.0_num
    sz(-2)=0.0_num
    sz(-1)=0.0_num
    sz(0)=0.0_num
    sz(1)=0.0_num
    sz(2)=0.0_num
    sz(3)=0.0_num

    ! --- computes coefficients for quantities centered between nodes
    oxint = 1.0_num-xint
    xintsq = xint*xint
    oxintsq = oxint*oxint
    sx(-1+dix) = onesixth*oxintsq*oxint
    sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
    sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
    sx( 2+dix) = onesixth*xintsq*xint
  
    ozint = 1.0_num-zint
    zintsq = zint*zint
    ozintsq = ozint*ozint
    sz(-1+diz) = onesixth*ozintsq*ozint
    sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
    sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
    sz( 2+diz) = onesixth*zintsq*zint

    ! --- computes coefficients difference
    dsx = sx - sx0
    dsz = sz - sz0
    
    ! --- computes min/max positions of current contributions
    ixmin = min(0,dix)-1
    ixmax = max(0,dix)+2
    izmin = min(0,diz)-1
    izmax = max(0,diz)+2

    ! --- add current contributions
    DO k=izmin, izmax
      DO i=ixmin, ixmax
        ic = iixp0+i
        kc = ikxp0+k
        
        ! --- Jx
        IF(i<ixmax) THEN
          sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k) )    ! Wx coefficient from esirkepov
          if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)         ! Integration of Wx along x 
          jx(ic,kc) = jx(ic,kc) + sdx(i,k)              ! Deposition on the current
        END IF
        
        ! -- Jy (2D Esirkepov scheme)
        jy(ic,kc) = jy(ic,kc) + wq*vy*invvol* &
        ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i) )
        
        ! --- Jz
        IF(k<izmax) THEN
          sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i))        ! Wz coefficient from esirkepov&
          if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)         ! Integration of Wz along z
          jz(ic,kc) = jz(ic,kc) + sdz(i,k)                  ! Deposition on the current
        END IF

! __ DEbug _______________________________
!  print*,'sum',sum(jx),sum(jy),sum(jz)
!  print*,'j',jy(ic,kc)
!  print*,'sdx',sdx(i,k),sdz(i,k)
!  print*,'wq',wqx,wqz,wq
!  print*,'s',sz0(k),sx0(i)
!  print*,'dsx',dsx(i),dsz(k)
!  read*        
        
      END DO
    END DO

  END DO
  
  
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)

End subroutine pxr_depose_jxjyjz_esirkepov2d_3_3


#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> 2D Current deposition with the method of Esirkepov at order 3
!> This function is semi-vectorized: only the first part
!> with the computation of indexes is vectorized
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx,jy,jz current arrays
!> @param[in] np number of particles
!> @param[in] xp,zp particle position arrays
!> @param[in] uxp,uyp,uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin,zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx,nz tile grid size
!> @param[in] nxguard,nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] lvect vector length for vectorization
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
!
subroutine pxr_depose_jxjyjz_esirkepov2d_svec_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose)
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none
  
  ! __ Parameter declaration _______________________________________________
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard, type_rz_depose
  integer(idp)                          :: lvect
  real(num), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num),dimension(:,:,:), allocatable :: sdx,sdz,sdy
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), parameter                  :: onethird=1./3.
  real(num), dimension(:), allocatable  :: sx, sx0, dsx, sz, sz0, dsz
  integer(isp)                          :: iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc
  integer(isp)                          :: icell, ndtodx, ndtodz
  integer(isp), dimension(lvect)        :: ixmin, ixmax, izmin, izmax, iixp0, ikxp0
  integer(isp)                          :: xl,xu,zl,zu,n,nn

  ! Parameter initialization
  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi
  allocate(sdx(lvect,-2:3,-2:3),sdz(lvect,-2:3,-2:3))
  ALLOCATE(sx(-2:3), sx0(-2:3), dsx(-2:3))
  ALLOCATE(sz(-2:3), sz0(-2:3), dsz(-2:3))
  sx0=0.0_num;sz0=0.0_num
  sdx=0.0_num;sdz=0.0_num

  ! Outer loop on particles with period LVEC
  DO ip=1,np, LVEC

#if defined __INTEL_COMPILER 
        !DIR$ ASSUME_ALIGNED xp:64,zp:64,gaminv:64
        !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
#elif defined __IBMBGQ__
        !IBM* ALIGN(64,xp,zp)
        !IBM* ALIGN(64,uxp,uyp,uzp)
        !IBM* ALIGN(64,ICELL)
#endif
#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD 
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
    !$DIR SIMD 
#endif

    ! Inner loop on particle
    DO n=1,MIN(LVECT,np-ip+1)
      nn=ip+n-1  
  
      ! --- computes current position in grid units
      x = (xp(nn)-xmin)*dxi
      z = (zp(nn)-zmin)*dzi
  
      ! --- computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)
  
      ! --- computes old position in grid units
      xold=x-dtsdx0*vx
      zold=z-dtsdz0*vz
  
      ! --- computes particles weights
      wq=q*w(nn)
      wqx = wq*invdtdx
      wqz = wq*invdtdz  

      ! --- finds node of cell containing particles for current positions
      iixp0(n)=floor(x)
      ikxp0(n)=floor(z)
  
      ! --- computes distance between particle and node for current positions
      xint=x-iixp0(n)
      zint=z-ikxp0(n)

      ! --- computes coefficients for node centered quantities
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx0( 2) = onesixth*xintsq*xint
  
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz0( 2) = onesixth*zintsq*zint

      ! --- finds node of cell containing particles for old positions
      iixp=floor(xold)
      ikxp=floor(zold)

      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      zint = zold-ikxp

      ! --- computes node separation between old and current positions
      dix = iixp-iixp0(n)
      diz = ikxp-ikxp0(n)

      ! --- zero out coefficients (needed because of different dix and diz for each particle)
      sx(-2)=0.0_num
      sx(-1)=0.0_num
      sx(0)=0.0_num
      sx(1)=0.0_num
      sx(2)=0.0_num
      sx(3)=0.0_num
      sz(-2)=0.0_num
      sz(-1)=0.0_num
      sz(0)=0.0_num
      sz(1)=0.0_num
      sz(2)=0.0_num
      sz(3)=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1+dix) = onesixth*oxintsq*oxint
      sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2+dix) = onesixth*xintsq*xint
  
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1+diz) = onesixth*ozintsq*ozint
      sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2+diz) = onesixth*zintsq*zint

      ! --- computes coefficients difference
      dsx(-2) = sx(-2) - sx0(-2)
      dsz(-2) = sz(-2) - sz0(-2)
      dsx(-1) = sx(-1) - sx0(-1)
      dsz(-1) = sz(-1) - sz0(-1)
      dsx(0) = sx(0) - sx0(0)
      dsz(0) = sz(0) - sz0(0)
      dsx(1) = sx(1) - sx0(1)
      dsz(1) = sz(1) - sz0(1)
      dsx(2) = sx(2) - sx0(2)
      dsz(2) = sz(2) - sz0(2)
      dsx(3) = sx(3) - sx0(3)
      dsz(3) = sz(3) - sz0(3)
    
      ! --- computes min/max positions of current contributions
      ixmin(n) = min(0,dix)-1
      ixmax(n) = max(0,dix)+2
      izmin(n) = min(0,diz)-1
      izmax(n) = max(0,diz)+2

      sdx(n,-2,-2)  = wqx*dsx(-2)*( sz0(-2) + 0.5*dsz(-2) )  
      sdx(n,-1,-2)  = wqx*dsx(-1)*( sz0(-2) + 0.5*dsz(-2) )  
      sdx(n,-1,-2)=sdx(n,-1,-2)+sdx(n,-1-1,-2)    
      sdx(n,0,-2)  = wqx*dsx(0)*( sz0(-2) + 0.5*dsz(-2) )  
      sdx(n,0,-2)=sdx(n,0,-2)+sdx(n,0-1,-2)    
      sdx(n,1,-2)  = wqx*dsx(1)*( sz0(-2) + 0.5*dsz(-2) )  
      sdx(n,1,-2)=sdx(n,1,-2)+sdx(n,1-1,-2)    
      sdx(n,2,-2)  = wqx*dsx(2)*( sz0(-2) + 0.5*dsz(-2) )  
      sdx(n,2,-2)=sdx(n,2,-2)+sdx(n,2-1,-2)    
      sdx(n,-2,-1)  = wqx*dsx(-2)*( sz0(-1) + 0.5*dsz(-1) )  
      sdx(n,-1,-1)  = wqx*dsx(-1)*( sz0(-1) + 0.5*dsz(-1) )  
      sdx(n,-1,-1)=sdx(n,-1,-1)+sdx(n,-1-1,-1)    
      sdx(n,0,-1)  = wqx*dsx(0)*( sz0(-1) + 0.5*dsz(-1) )  
      sdx(n,0,-1)=sdx(n,0,-1)+sdx(n,0-1,-1)    
      sdx(n,1,-1)  = wqx*dsx(1)*( sz0(-1) + 0.5*dsz(-1) )  
      sdx(n,1,-1)=sdx(n,1,-1)+sdx(n,1-1,-1)    
      sdx(n,2,-1)  = wqx*dsx(2)*( sz0(-1) + 0.5*dsz(-1) )  
      sdx(n,2,-1)=sdx(n,2,-1)+sdx(n,2-1,-1)    
      sdx(n,-2,0)  = wqx*dsx(-2)*( sz0(0) + 0.5*dsz(0) )  
      sdx(n,-1,0)  = wqx*dsx(-1)*( sz0(0) + 0.5*dsz(0) )  
      sdx(n,-1,0)=sdx(n,-1,0)+sdx(n,-1-1,0)    
      sdx(n,0,0)  = wqx*dsx(0)*( sz0(0) + 0.5*dsz(0) )  
      sdx(n,0,0)=sdx(n,0,0)+sdx(n,0-1,0)    
      sdx(n,1,0)  = wqx*dsx(1)*( sz0(0) + 0.5*dsz(0) )  
      sdx(n,1,0)=sdx(n,1,0)+sdx(n,1-1,0)    
      sdx(n,2,0)  = wqx*dsx(2)*( sz0(0) + 0.5*dsz(0) )  
      sdx(n,2,0)=sdx(n,2,0)+sdx(n,2-1,0)    
      sdx(n,-2,1)  = wqx*dsx(-2)*( sz0(1) + 0.5*dsz(1) )  
      sdx(n,-1,1)  = wqx*dsx(-1)*( sz0(1) + 0.5*dsz(1) )  
      sdx(n,-1,1)=sdx(n,-1,1)+sdx(n,-1-1,1)    
      sdx(n,0,1)  = wqx*dsx(0)*( sz0(1) + 0.5*dsz(1) )  
      sdx(n,0,1)=sdx(n,0,1)+sdx(n,0-1,1)    
      sdx(n,1,1)  = wqx*dsx(1)*( sz0(1) + 0.5*dsz(1) )  
      sdx(n,1,1)=sdx(n,1,1)+sdx(n,1-1,1)    
      sdx(n,2,1)  = wqx*dsx(2)*( sz0(1) + 0.5*dsz(1) )  
      sdx(n,2,1)=sdx(n,2,1)+sdx(n,2-1,1)    
      sdx(n,-2,2)  = wqx*dsx(-2)*( sz0(2) + 0.5*dsz(2) )  
      sdx(n,-1,2)  = wqx*dsx(-1)*( sz0(2) + 0.5*dsz(2) )  
      sdx(n,-1,2)=sdx(n,-1,2)+sdx(n,-1-1,2)    
      sdx(n,0,2)  = wqx*dsx(0)*( sz0(2) + 0.5*dsz(2) )  
      sdx(n,0,2)=sdx(n,0,2)+sdx(n,0-1,2)    
      sdx(n,1,2)  = wqx*dsx(1)*( sz0(2) + 0.5*dsz(2) )  
      sdx(n,1,2)=sdx(n,1,2)+sdx(n,1-1,2)    
      sdx(n,2,2)  = wqx*dsx(2)*( sz0(2) + 0.5*dsz(2) )  
      sdx(n,2,2)=sdx(n,2,2)+sdx(n,2-1,2)    
      sdx(n,-2,3)  = wqx*dsx(-2)*( sz0(3) + 0.5*dsz(3) )  
      sdx(n,-1,3)  = wqx*dsx(-1)*( sz0(3) + 0.5*dsz(3) )  
      sdx(n,-1,3)=sdx(n,-1,3)+sdx(n,-1-1,3)    
      sdx(n,0,3)  = wqx*dsx(0)*( sz0(3) + 0.5*dsz(3) )  
      sdx(n,0,3)=sdx(n,0,3)+sdx(n,0-1,3)    
      sdx(n,1,3)  = wqx*dsx(1)*( sz0(3) + 0.5*dsz(3) )  
      sdx(n,1,3)=sdx(n,1,3)+sdx(n,1-1,3)    
      sdx(n,2,3)  = wqx*dsx(2)*( sz0(3) + 0.5*dsz(3) )  
      sdx(n,2,3)=sdx(n,2,3)+sdx(n,2-1,3)   

      sdy(n,-2,-2) = wq*vy*invvol* &
      ( (sz0(-2)+0.5*dsz(-2))*sx0(-2) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-2) )
      sdy(n,-1,-2) = wq*vy*invvol* &
      ( (sz0(-2)+0.5*dsz(-2))*sx0(-1) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-1) )
      sdy(n,0,-2) = wq*vy*invvol* &
      ( (sz0(-2)+0.5*dsz(-2))*sx0(0) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(0) )
      sdy(n,1,-2) = wq*vy*invvol* &
      ( (sz0(-2)+0.5*dsz(-2))*sx0(1) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(1) )
      sdy(n,2,-2) = wq*vy*invvol* &
      ( (sz0(-2)+0.5*dsz(-2))*sx0(2) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(2) )
      sdy(n,3,-2) = wq*vy*invvol* &
      ( (sz0(-2)+0.5*dsz(-2))*sx0(3) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(3) )
      sdy(n,-2,-1) = wq*vy*invvol* &
      ( (sz0(-1)+0.5*dsz(-1))*sx0(-2) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-2) )
      sdy(n,-1,-1) = wq*vy*invvol* &
      ( (sz0(-1)+0.5*dsz(-1))*sx0(-1) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-1) )
      sdy(n,0,-1) = wq*vy*invvol* &
      ( (sz0(-1)+0.5*dsz(-1))*sx0(0) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(0) )
      sdy(n,1,-1) = wq*vy*invvol* &
      ( (sz0(-1)+0.5*dsz(-1))*sx0(1) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(1) )
      sdy(n,2,-1) = wq*vy*invvol* &
      ( (sz0(-1)+0.5*dsz(-1))*sx0(2) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(2) )
      sdy(n,3,-1) = wq*vy*invvol* &
      ( (sz0(-1)+0.5*dsz(-1))*sx0(3) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(3) )
      sdy(n,-2,0) = wq*vy*invvol* &
      ( (sz0(0)+0.5*dsz(0))*sx0(-2) + (0.5*sz0(0)+onethird*dsz(0))*dsx(-2) )
      sdy(n,-1,0) = wq*vy*invvol* &
      ( (sz0(0)+0.5*dsz(0))*sx0(-1) + (0.5*sz0(0)+onethird*dsz(0))*dsx(-1) )
      sdy(n,0,0) = wq*vy*invvol* &
      ( (sz0(0)+0.5*dsz(0))*sx0(0) + (0.5*sz0(0)+onethird*dsz(0))*dsx(0) )
      sdy(n,1,0) = wq*vy*invvol* &
      ( (sz0(0)+0.5*dsz(0))*sx0(1) + (0.5*sz0(0)+onethird*dsz(0))*dsx(1) )
      sdy(n,2,0) = wq*vy*invvol* &
      ( (sz0(0)+0.5*dsz(0))*sx0(2) + (0.5*sz0(0)+onethird*dsz(0))*dsx(2) )
      sdy(n,3,0) = wq*vy*invvol* &
      ( (sz0(0)+0.5*dsz(0))*sx0(3) + (0.5*sz0(0)+onethird*dsz(0))*dsx(3) )
      sdy(n,-2,1) = wq*vy*invvol* &
      ( (sz0(1)+0.5*dsz(1))*sx0(-2) + (0.5*sz0(1)+onethird*dsz(1))*dsx(-2) )
      sdy(n,-1,1) = wq*vy*invvol* &
      ( (sz0(1)+0.5*dsz(1))*sx0(-1) + (0.5*sz0(1)+onethird*dsz(1))*dsx(-1) )
      sdy(n,0,1) = wq*vy*invvol* &
      ( (sz0(1)+0.5*dsz(1))*sx0(0) + (0.5*sz0(1)+onethird*dsz(1))*dsx(0) )
      sdy(n,1,1) = wq*vy*invvol* &
      ( (sz0(1)+0.5*dsz(1))*sx0(1) + (0.5*sz0(1)+onethird*dsz(1))*dsx(1) )
      sdy(n,2,1) = wq*vy*invvol* &
      ( (sz0(1)+0.5*dsz(1))*sx0(2) + (0.5*sz0(1)+onethird*dsz(1))*dsx(2) )
      sdy(n,3,1) = wq*vy*invvol* &
      ( (sz0(1)+0.5*dsz(1))*sx0(3) + (0.5*sz0(1)+onethird*dsz(1))*dsx(3) )
      sdy(n,-2,2) = wq*vy*invvol* &
      ( (sz0(2)+0.5*dsz(2))*sx0(-2) + (0.5*sz0(2)+onethird*dsz(2))*dsx(-2) )
      sdy(n,-1,2) = wq*vy*invvol* &
      ( (sz0(2)+0.5*dsz(2))*sx0(-1) + (0.5*sz0(2)+onethird*dsz(2))*dsx(-1) )
      sdy(n,0,2) = wq*vy*invvol* &
      ( (sz0(2)+0.5*dsz(2))*sx0(0) + (0.5*sz0(2)+onethird*dsz(2))*dsx(0) )
      sdy(n,1,2) = wq*vy*invvol* &
      ( (sz0(2)+0.5*dsz(2))*sx0(1) + (0.5*sz0(2)+onethird*dsz(2))*dsx(1) )
      sdy(n,2,2) = wq*vy*invvol* &
      ( (sz0(2)+0.5*dsz(2))*sx0(2) + (0.5*sz0(2)+onethird*dsz(2))*dsx(2) )
      sdy(n,3,2) = wq*vy*invvol* &
      ( (sz0(2)+0.5*dsz(2))*sx0(3) + (0.5*sz0(2)+onethird*dsz(2))*dsx(3) )
      sdy(n,-2,3) = wq*vy*invvol* &
      ( (sz0(3)+0.5*dsz(3))*sx0(-2) + (0.5*sz0(3)+onethird*dsz(3))*dsx(-2) )
      sdy(n,-1,3) = wq*vy*invvol* &
      ( (sz0(3)+0.5*dsz(3))*sx0(-1) + (0.5*sz0(3)+onethird*dsz(3))*dsx(-1) )
      sdy(n,0,3) = wq*vy*invvol* &
      ( (sz0(3)+0.5*dsz(3))*sx0(0) + (0.5*sz0(3)+onethird*dsz(3))*dsx(0) )
      sdy(n,1,3) = wq*vy*invvol* &
      ( (sz0(3)+0.5*dsz(3))*sx0(1) + (0.5*sz0(3)+onethird*dsz(3))*dsx(1) )
      sdy(n,2,3) = wq*vy*invvol* &
      ( (sz0(3)+0.5*dsz(3))*sx0(2) + (0.5*sz0(3)+onethird*dsz(3))*dsx(2) )
      sdy(n,3,3) = wq*vy*invvol* &
      ( (sz0(3)+0.5*dsz(3))*sx0(3) + (0.5*sz0(3)+onethird*dsz(3))*dsx(3) )

      sdz(n,-2,-2)=wqz*dsz(-2)*(sx0(-2)+0.5*dsx(-2))    
      sdz(n,-1,-2)=wqz*dsz(-2)*(sx0(-1)+0.5*dsx(-1))    
      sdz(n,0,-2)=wqz*dsz(-2)*(sx0(0)+0.5*dsx(0))    
      sdz(n,1,-2)=wqz*dsz(-2)*(sx0(1)+0.5*dsx(1))    
      sdz(n,2,-2)=wqz*dsz(-2)*(sx0(2)+0.5*dsx(2))    
      sdz(n,3,-2)=wqz*dsz(-2)*(sx0(3)+0.5*dsx(3))    
      sdz(n,-2,-1)=wqz*dsz(-1)*(sx0(-2)+0.5*dsx(-2))    
      sdz(n,-2,-1)=sdz(n,-2,-1)+sdz(n,-2,-1-1)    
      sdz(n,-1,-1)=wqz*dsz(-1)*(sx0(-1)+0.5*dsx(-1))    
      sdz(n,-1,-1)=sdz(n,-1,-1)+sdz(n,-1,-1-1)    
      sdz(n,0,-1)=wqz*dsz(-1)*(sx0(0)+0.5*dsx(0))    
      sdz(n,0,-1)=sdz(n,0,-1)+sdz(n,0,-1-1)    
      sdz(n,1,-1)=wqz*dsz(-1)*(sx0(1)+0.5*dsx(1))    
      sdz(n,1,-1)=sdz(n,1,-1)+sdz(n,1,-1-1)    
      sdz(n,2,-1)=wqz*dsz(-1)*(sx0(2)+0.5*dsx(2))    
      sdz(n,2,-1)=sdz(n,2,-1)+sdz(n,2,-1-1)    
      sdz(n,3,-1)=wqz*dsz(-1)*(sx0(3)+0.5*dsx(3))    
      sdz(n,3,-1)=sdz(n,3,-1)+sdz(n,3,-1-1)    
      sdz(n,-2,0)=wqz*dsz(0)*(sx0(-2)+0.5*dsx(-2))    
      sdz(n,-2,0)=sdz(n,-2,0)+sdz(n,-2,0-1)    
      sdz(n,-1,0)=wqz*dsz(0)*(sx0(-1)+0.5*dsx(-1))    
      sdz(n,-1,0)=sdz(n,-1,0)+sdz(n,-1,0-1)    
      sdz(n,0,0)=wqz*dsz(0)*(sx0(0)+0.5*dsx(0))    
      sdz(n,0,0)=sdz(n,0,0)+sdz(n,0,0-1)    
      sdz(n,1,0)=wqz*dsz(0)*(sx0(1)+0.5*dsx(1))    
      sdz(n,1,0)=sdz(n,1,0)+sdz(n,1,0-1)    
      sdz(n,2,0)=wqz*dsz(0)*(sx0(2)+0.5*dsx(2))    
      sdz(n,2,0)=sdz(n,2,0)+sdz(n,2,0-1)    
      sdz(n,3,0)=wqz*dsz(0)*(sx0(3)+0.5*dsx(3))    
      sdz(n,3,0)=sdz(n,3,0)+sdz(n,3,0-1)    
      sdz(n,-2,1)=wqz*dsz(1)*(sx0(-2)+0.5*dsx(-2))    
      sdz(n,-2,1)=sdz(n,-2,1)+sdz(n,-2,1-1)    
      sdz(n,-1,1)=wqz*dsz(1)*(sx0(-1)+0.5*dsx(-1))    
      sdz(n,-1,1)=sdz(n,-1,1)+sdz(n,-1,1-1)    
      sdz(n,0,1)=wqz*dsz(1)*(sx0(0)+0.5*dsx(0))    
      sdz(n,0,1)=sdz(n,0,1)+sdz(n,0,1-1)    
      sdz(n,1,1)=wqz*dsz(1)*(sx0(1)+0.5*dsx(1))    
      sdz(n,1,1)=sdz(n,1,1)+sdz(n,1,1-1)    
      sdz(n,2,1)=wqz*dsz(1)*(sx0(2)+0.5*dsx(2))    
      sdz(n,2,1)=sdz(n,2,1)+sdz(n,2,1-1)    
      sdz(n,3,1)=wqz*dsz(1)*(sx0(3)+0.5*dsx(3))    
      sdz(n,3,1)=sdz(n,3,1)+sdz(n,3,1-1)    
      sdz(n,-2,2)=wqz*dsz(2)*(sx0(-2)+0.5*dsx(-2))    
      sdz(n,-2,2)=sdz(n,-2,2)+sdz(n,-2,2-1)    
      sdz(n,-1,2)=wqz*dsz(2)*(sx0(-1)+0.5*dsx(-1))    
      sdz(n,-1,2)=sdz(n,-1,2)+sdz(n,-1,2-1)    
      sdz(n,0,2)=wqz*dsz(2)*(sx0(0)+0.5*dsx(0))    
      sdz(n,0,2)=sdz(n,0,2)+sdz(n,0,2-1)    
      sdz(n,1,2)=wqz*dsz(2)*(sx0(1)+0.5*dsx(1))    
      sdz(n,1,2)=sdz(n,1,2)+sdz(n,1,2-1)    
      sdz(n,2,2)=wqz*dsz(2)*(sx0(2)+0.5*dsx(2))    
      sdz(n,2,2)=sdz(n,2,2)+sdz(n,2,2-1)    
      sdz(n,3,2)=wqz*dsz(2)*(sx0(3)+0.5*dsx(3))    
      sdz(n,3,2)=sdz(n,3,2)+sdz(n,3,2-1) 

    ENDDO

    ! Inner loop on particle
    DO n=1,MIN(LVECT,np-ip+1)

      ! --- add current contributions
      DO k=izmin(n), izmax(n)
        DO i=ixmin(n), ixmax(n)
          ic = iixp0(n)+i
          kc = ikxp0(n)+k
        
          ! --- Jx
          jx(ic,kc) = jx(ic,kc) + sdx(n,i,k)  ! Deposition on the current
        
          ! -- Jy (2D Esirkepov scheme)
          jy(ic,kc) = jy(ic,kc) + sdy(n,i,k)
        
          ! --- Jz
          jz(ic,kc) = jz(ic,kc) + sdz(n,i,k)  ! Deposition on the current
        
        END DO
      END DO
      
    ENDDO
  END DO
  DEALLOCATE(sdx,sdz,sx,sx0,dsx,sz,sz0,dsz)

End subroutine pxr_depose_jxjyjz_esirkepov2d_svec_3_3
#endif

#if defined(DEV)
! ________________________________________________________________________________________
!> @brief
!> Vectorized 2D Current deposition with the method of Esirkepov at order 3.
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
! Input parameters:
!> @param[inout] jx,jy,jz current arrays
!> @param[in] np number of particles
!> @param[in] xp,zp particle position arrays
!> @param[in] uxp,uyp,uzp particle momentum arrays
!> @param[in] gaminv inverse of the gamma factor
!> @param[in] w particle weight
!> @param[in] q particle charge
!> @param[in] xmin,zmin minimal boundaries of the tile
!> @param[in] dt, dx, dz time and space discretization
!> @param[in] nx,nz tile grid size
!> @param[in] nxguard,nzguard guard cell numbers
!> @param[in] nox, noz shape factor order (useless here but kept for common interface)
!> @param[in] lvect vector length for vectorization
!> @param[in] l_particles_weight to take into account the particle weight
!> @param[in] l4symtry (useless here bur kept for common interface)
!> @param[in] l_2drz  (useless here bur kept for common interface)
!> @param[in] type_rz_depose (useless here bur kept for common interface)
subroutine pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3(jx,jy,jz,np,xp,zp,uxp,uyp,uzp,&
          gaminv,w,q,xmin,zmin, &
          dt,dx,dz,nx,nz,nxguard,nzguard,&
          nox,noz,lvect,l_particles_weight,l4symtry,l_2drz,type_rz_depose) !#do not parse
! ________________________________________________________________________________________

  USE omp_lib
  USE constants
  implicit none
  
  ! __ Parameter declaration ____________________________________________________
  ! In/out parameters
  integer(idp)                          :: np,nx,nz,nox,noz,nxguard,nzguard,type_rz_depose
  integer(idp)                          :: lvect
  real(num), dimension((1+nx+2*nxguard)*(1+nz+2*nzguard)), intent(in out) :: jx,jy,jz
  real(num), dimension(np)              :: xp,zp,uxp,uyp,uzp,gaminv,w
  real(num)                             :: q,dt,dx,dz,xmin,zmin
  LOGICAL(lp)                           :: l_particles_weight,l4symtry,l_2drz
  
  ! Local parameters
  real(num)                             :: dxi,dzi,dtsdx,dtsdz,xint,zint
  real(num)                             :: xold,zold,rold,xmid,zmid,x,z,c,s,wq,wqx,wqz
  real(num)                             :: tmp,vx,vy,vz,dts2dx,dts2dz
  real(num)                             :: invvol,invdtdx,invdtdz
  real(num)                             :: oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
  real(num)                             :: dtsdx0,dtsdz0,dts2dx0,dts2dz0
  real(num), parameter                  :: onesixth=1./6.,twothird=2./3.
  real(num), parameter                  :: onethird=1./3.  
  real(num), dimension(-2:3)            :: sx, sx0, dsx, sz, sz0, dsz
  integer(idp)                          :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc
  integer(idp)                          :: ixmin, ixmax, izmin, izmax, ndtodx, ndtodz
  integer(idp)                          :: xl,xu,zl,zu
  integer(isp)                          :: ngridx, ncx,ncz,ncells
  integer(isp)                          :: n,nn,nv  
  integer(isp)                          :: iixporig, ikxporig, igrid, ix, iz, orig, nnx
  INTEGER(isp), DIMENSION(LVEC,3)       :: ICELL  
  REAL(num),DIMENSION(:,:), ALLOCATABLE :: jxcells,jycells,jzcells 
  REAL(num),DIMENSION(:,:), ALLOCATABLE :: sdx,sdy,sdz   

  ! __ Parameter initialization __________________________________________________

  dxi = 1.0_num/dx
  dzi = 1.0_num/dz
  dtsdx0 = dt*dxi
  dtsdz0 = dt*dzi
  invvol = 1.0_num/(dx*dz)
  invdtdx = 1.0_num/(dt*dz)
  invdtdz = 1.0_num/(dt*dx)
  dtsdz0 = dt*dzi

  sx0(-2:3)=0._num
  sz0(-2:3)=0._num

  ngridx=nx+1+2*nxguard
  ncx=nx+1+2*nxguard
  ncz=nz+1+2*nzguard
  NCELLS=ncx*ncz
  iixporig=-nxguard
  ikxporig=-nzguard
  ALLOCATE(jxcells(8,NCELLS),jycells(8,NCELLS),jzcells(8,NCELLS))
  
  
  jxcells = 0._num
  jycells = 0._num
  jzcells = 0._num
  
  ALLOCATE(sdx(lvect,48),sdy(lvect,48),sdz(lvect,48))
  
  sdx = 0._num
  sdz = 0._num

  nnx = ngridx
  orig=(nxguard+iixporig) + (nzguard+ikxporig)*nnx
  
  ! Outer loop on particles with period LVEC
  DO ip=1,np, LVEC

#if defined __INTEL_COMPILER 
        !DIR$ ASSUME_ALIGNED xp:64,zp:64,gaminv:64
        !DIR$ ASSUME_ALIGNED uxp:64,uyp:64,uzp:64
        !DIR$ ASSUME_ALIGNED ICELL:64
#elif defined __IBMBGQ__
        !IBM* ALIGN(64,xp,yp,zp)
        !IBM* ALIGN(64,uxp,uyp,uzp)
        !IBM* ALIGN(64,ICELL)
#endif
#if defined _OPENMP && _OPENMP>=201307
    !$OMP SIMD 
#elif defined __IBMBGQ__
    !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
    !$DIR SIMD 
#endif

    ! Inner loop on particle
    DO n=1,MIN(LVECT,np-ip+1)
      nn=ip+n-1  

      ! --- computes current position in grid units
      x = (xp(nn)-xmin)*dxi
      z = (zp(nn)-zmin)*dzi
  
      ! --- computes velocity
      vx = uxp(nn)*gaminv(nn)
      vy = uyp(nn)*gaminv(nn)
      vz = uzp(nn)*gaminv(nn)
  
      ! --- computes old position in grid units
      xold=x-dtsdx0*vx
      zold=z-dtsdz0*vz
  
      ! --- computes particles weights
      wq=q*w(ip)
      wqx = wq*invdtdx
      wqz = wq*invdtdz
  
      ! --- finds node of cell containing particles for current positions
      iixp0=floor(x)
      ikxp0=floor(z)
  
      ! --- computes distance between particle and node for current positions
      xint=x-iixp0
      zint=z-ikxp0

      ! --- computes coefficients for node centered quantities
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx0(-1) = onesixth*oxintsq*oxint
      sx0( 0) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx0( 1) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx0( 2) = onesixth*xintsq*xint
  
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz0(-1) = onesixth*ozintsq*ozint
      sz0( 0) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz0( 1) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz0( 2) = onesixth*zintsq*zint
  
      ! --- finds node of cell containing particles for old positions
      iixp=floor(xold)
      ikxp=floor(zold)
  
      ! --- computes distance between particle and node for old positions
      xint = xold-iixp
      zint = zold-ikxp
  
      ! --- computes node separation between old and current positions
      dix = iixp-iixp0
      diz = ikxp-ikxp0

      ! --- zero out coefficients (needed because of different dix and diz for each particle)
      sx(-2)=0.0_num
      sx(-1)=0.0_num
      sx(0)=0.0_num
      sx(1)=0.0_num
      sx(2)=0.0_num
      sx(3)=0.0_num
      sz(-2)=0.0_num
      sz(-1)=0.0_num
      sz(0)=0.0_num
      sz(1)=0.0_num
      sz(2)=0.0_num
      sz(3)=0.0_num

      ! --- computes coefficients for quantities centered between nodes
      oxint = 1.0_num-xint
      xintsq = xint*xint
      oxintsq = oxint*oxint
      sx(-1+dix) = onesixth*oxintsq*oxint
      sx( 0+dix) = twothird-xintsq*(1.0_num-xint*0.5_num)
      sx( 1+dix) = twothird-oxintsq*(1.0_num-oxint*0.5_num)
      sx( 2+dix) = onesixth*xintsq*xint
  
      ozint = 1.0_num-zint
      zintsq = zint*zint
      ozintsq = ozint*ozint
      sz(-1+diz) = onesixth*ozintsq*ozint
      sz( 0+diz) = twothird-zintsq*(1.0_num-zint*0.5_num)
      sz( 1+diz) = twothird-ozintsq*(1.0_num-ozint*0.5_num)
      sz( 2+diz) = onesixth*zintsq*zint

      ! --- computes coefficients difference
      dsx(-2) = sx(-2) - sx0(-2)
      dsz(-2) = sz(-2) - sz0(-2)
      dsx(-1) = sx(-1) - sx0(-1)
      dsz(-1) = sz(-1) - sz0(-1)
      dsx(0) = sx(0) - sx0(0)
      dsz(0) = sz(0) - sz0(0)
      dsx(1) = sx(1) - sx0(1)
      dsz(1) = sz(1) - sz0(1)
      dsx(2) = sx(2) - sx0(2)
      dsz(2) = sz(2) - sz0(2)
      dsx(3) = sx(3) - sx0(3)
      dsz(3) = sz(3) - sz0(3)
      
      ! --- Position of the first cell (-2,-2,-2)
      icell(n,1) = (iixp0-iixporig-1)+(ikxp0-ikxporig-2)*ncx
      
      ! --- Weight
      sdx(n,1)  = wqx*dsx(-2)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n,2)  = wqx*dsx(-1)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n,2) = sdx(n,2)+sdx(n,1)
      sdx(n,3)  = wqx*dsx(0)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n,3) = sdx(n,3)+sdx(n,2)
      sdx(n,4)  = wqx*dsx(1)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n,4) = sdx(n,4)+sdx(n,3)
      sdx(n,5)  = wqx*dsx(2)*( sz0(-2) + 0.5*dsz(-2) )
      sdx(n,5) = sdx(n,5)+sdx(n,4)
      sdx(n,9)  = wqx*dsx(-2)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n,10)  = wqx*dsx(-1)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n,10) = sdx(n,10)+sdx(n,9)
      sdx(n,11)  = wqx*dsx(0)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n,11) = sdx(n,11)+sdx(n,10)
      sdx(n,12)  = wqx*dsx(1)*( sz0(-1) + 0.5*dsz(-1) )
      sdx(n,12) = sdx(n,12)+sdx(n,11)
      sdx(n,13)  = wqx*dsx(2)*( sz0(-1) + 0.5*dsz(-1) )
    sdx(n,13) = sdx(n,13)+sdx(n,12)
    sdx(n,17)  = wqx*dsx(-2)*( sz0(0) + 0.5*dsz(0) )
    sdx(n,18)  = wqx*dsx(-1)*( sz0(0) + 0.5*dsz(0) )
    sdx(n,18) = sdx(n,18)+sdx(n,17)
    sdx(n,19)  = wqx*dsx(0)*( sz0(0) + 0.5*dsz(0) )
    sdx(n,19) = sdx(n,19)+sdx(n,18)
    sdx(n,20)  = wqx*dsx(1)*( sz0(0) + 0.5*dsz(0) )
    sdx(n,20) = sdx(n,20)+sdx(n,19)
    sdx(n,21)  = wqx*dsx(2)*( sz0(0) + 0.5*dsz(0) )
    sdx(n,21) = sdx(n,21)+sdx(n,20)
    sdx(n,25)  = wqx*dsx(-2)*( sz0(1) + 0.5*dsz(1) )
    sdx(n,26)  = wqx*dsx(-1)*( sz0(1) + 0.5*dsz(1) )
    sdx(n,26) = sdx(n,26)+sdx(n,25)
    sdx(n,27)  = wqx*dsx(0)*( sz0(1) + 0.5*dsz(1) )
    sdx(n,27) = sdx(n,27)+sdx(n,26)
    sdx(n,28)  = wqx*dsx(1)*( sz0(1) + 0.5*dsz(1) )
    sdx(n,28) = sdx(n,28)+sdx(n,27)
    sdx(n,29)  = wqx*dsx(2)*( sz0(1) + 0.5*dsz(1) )
    sdx(n,29) = sdx(n,29)+sdx(n,28)
    sdx(n,33)  = wqx*dsx(-2)*( sz0(2) + 0.5*dsz(2) )
    sdx(n,34)  = wqx*dsx(-1)*( sz0(2) + 0.5*dsz(2) )
    sdx(n,34) = sdx(n,34)+sdx(n,33)
    sdx(n,35)  = wqx*dsx(0)*( sz0(2) + 0.5*dsz(2) )
    sdx(n,35) = sdx(n,35)+sdx(n,34)
    sdx(n,36)  = wqx*dsx(1)*( sz0(2) + 0.5*dsz(2) )
    sdx(n,36) = sdx(n,36)+sdx(n,35)
    sdx(n,37)  = wqx*dsx(2)*( sz0(2) + 0.5*dsz(2) )
    sdx(n,37) = sdx(n,37)+sdx(n,36)
    sdx(n,41)  = wqx*dsx(-2)*( sz0(3) + 0.5*dsz(3) )
    sdx(n,42)  = wqx*dsx(-1)*( sz0(3) + 0.5*dsz(3) )
    sdx(n,42) = sdx(n,42)+sdx(n,41)
    sdx(n,43)  = wqx*dsx(0)*( sz0(3) + 0.5*dsz(3) )
    sdx(n,43) = sdx(n,43)+sdx(n,42)
    sdx(n,44)  = wqx*dsx(1)*( sz0(3) + 0.5*dsz(3) )
    sdx(n,44) = sdx(n,44)+sdx(n,43)
    sdx(n,45)  = wqx*dsx(2)*( sz0(3) + 0.5*dsz(3) )
    sdx(n,45) = sdx(n,45)+sdx(n,44)
  
    sdy(n,1) = wq*vy*invvol* &
    ( (sz0(-2)+0.5*dsz(-2))*sx0(-2) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-2))
    sdy(n,2) = wq*vy*invvol* &
    ( (sz0(-2)+0.5*dsz(-2))*sx0(-1) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(-1))
    sdy(n,3) = wq*vy*invvol* &
    ( (sz0(-2)+0.5*dsz(-2))*sx0(0) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(0))
    sdy(n,4) = wq*vy*invvol* &
    ( (sz0(-2)+0.5*dsz(-2))*sx0(1) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(1))
    sdy(n,5) = wq*vy*invvol* &
    ( (sz0(-2)+0.5*dsz(-2))*sx0(2) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(2))
    sdy(n,6) = wq*vy*invvol* &
    ( (sz0(-2)+0.5*dsz(-2))*sx0(3) + (0.5*sz0(-2)+onethird*dsz(-2))*dsx(3))
    sdy(n,9) = wq*vy*invvol* &
    ( (sz0(-1)+0.5*dsz(-1))*sx0(-2) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-2))
    sdy(n,10) = wq*vy*invvol* &
    ( (sz0(-1)+0.5*dsz(-1))*sx0(-1) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(-1))
    sdy(n,11) = wq*vy*invvol* &
    ( (sz0(-1)+0.5*dsz(-1))*sx0(0) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(0))
    sdy(n,12) = wq*vy*invvol* &
    ( (sz0(-1)+0.5*dsz(-1))*sx0(1) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(1))
    sdy(n,13) = wq*vy*invvol* &
    ( (sz0(-1)+0.5*dsz(-1))*sx0(2) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(2))
    sdy(n,14) = wq*vy*invvol* &
    ( (sz0(-1)+0.5*dsz(-1))*sx0(3) + (0.5*sz0(-1)+onethird*dsz(-1))*dsx(3))
    sdy(n,17) = wq*vy*invvol* &
    ( (sz0(0)+0.5*dsz(0))*sx0(-2) + (0.5*sz0(0)+onethird*dsz(0))*dsx(-2))
    sdy(n,18) = wq*vy*invvol* &
    ( (sz0(0)+0.5*dsz(0))*sx0(-1) + (0.5*sz0(0)+onethird*dsz(0))*dsx(-1))
    sdy(n,19) = wq*vy*invvol* &
    ( (sz0(0)+0.5*dsz(0))*sx0(0) + (0.5*sz0(0)+onethird*dsz(0))*dsx(0))
    sdy(n,20) = wq*vy*invvol* &
    ( (sz0(0)+0.5*dsz(0))*sx0(1) + (0.5*sz0(0)+onethird*dsz(0))*dsx(1))
    sdy(n,21) = wq*vy*invvol* &
    ( (sz0(0)+0.5*dsz(0))*sx0(2) + (0.5*sz0(0)+onethird*dsz(0))*dsx(2))
    sdy(n,22) = wq*vy*invvol* &
    ( (sz0(0)+0.5*dsz(0))*sx0(3) + (0.5*sz0(0)+onethird*dsz(0))*dsx(3))
    sdy(n,25) = wq*vy*invvol* &
    ( (sz0(1)+0.5*dsz(1))*sx0(-2) + (0.5*sz0(1)+onethird*dsz(1))*dsx(-2))
    sdy(n,26) = wq*vy*invvol* &
    ( (sz0(1)+0.5*dsz(1))*sx0(-1) + (0.5*sz0(1)+onethird*dsz(1))*dsx(-1))
    sdy(n,27) = wq*vy*invvol* &
    ( (sz0(1)+0.5*dsz(1))*sx0(0) + (0.5*sz0(1)+onethird*dsz(1))*dsx(0))
    sdy(n,28) = wq*vy*invvol* &
    ( (sz0(1)+0.5*dsz(1))*sx0(1) + (0.5*sz0(1)+onethird*dsz(1))*dsx(1))
    sdy(n,29) = wq*vy*invvol* &
    ( (sz0(1)+0.5*dsz(1))*sx0(2) + (0.5*sz0(1)+onethird*dsz(1))*dsx(2))
    sdy(n,30) = wq*vy*invvol* &
    ( (sz0(1)+0.5*dsz(1))*sx0(3) + (0.5*sz0(1)+onethird*dsz(1))*dsx(3))
    sdy(n,33) = wq*vy*invvol* &
    ( (sz0(2)+0.5*dsz(2))*sx0(-2) + (0.5*sz0(2)+onethird*dsz(2))*dsx(-2))
    sdy(n,34) = wq*vy*invvol* &
    ( (sz0(2)+0.5*dsz(2))*sx0(-1) + (0.5*sz0(2)+onethird*dsz(2))*dsx(-1))
    sdy(n,35) = wq*vy*invvol* &
    ( (sz0(2)+0.5*dsz(2))*sx0(0) + (0.5*sz0(2)+onethird*dsz(2))*dsx(0))
    sdy(n,36) = wq*vy*invvol* &
    ( (sz0(2)+0.5*dsz(2))*sx0(1) + (0.5*sz0(2)+onethird*dsz(2))*dsx(1))
    sdy(n,37) = wq*vy*invvol* &
    ( (sz0(2)+0.5*dsz(2))*sx0(2) + (0.5*sz0(2)+onethird*dsz(2))*dsx(2))
    sdy(n,38) = wq*vy*invvol* &
    ( (sz0(2)+0.5*dsz(2))*sx0(3) + (0.5*sz0(2)+onethird*dsz(2))*dsx(3))
    sdy(n,41) = wq*vy*invvol* &
    ( (sz0(3)+0.5*dsz(3))*sx0(-2) + (0.5*sz0(3)+onethird*dsz(3))*dsx(-2))
    sdy(n,42) = wq*vy*invvol* &
    ( (sz0(3)+0.5*dsz(3))*sx0(-1) + (0.5*sz0(3)+onethird*dsz(3))*dsx(-1))
    sdy(n,43) = wq*vy*invvol* &
    ( (sz0(3)+0.5*dsz(3))*sx0(0) + (0.5*sz0(3)+onethird*dsz(3))*dsx(0))
    sdy(n,44) = wq*vy*invvol* &
    ( (sz0(3)+0.5*dsz(3))*sx0(1) + (0.5*sz0(3)+onethird*dsz(3))*dsx(1))
    sdy(n,45) = wq*vy*invvol* &
    ( (sz0(3)+0.5*dsz(3))*sx0(2) + (0.5*sz0(3)+onethird*dsz(3))*dsx(2))
    sdy(n,46) = wq*vy*invvol* &
    ( (sz0(3)+0.5*dsz(3))*sx0(3) + (0.5*sz0(3)+onethird*dsz(3))*dsx(3))
  
    sdz(n,1)  = wqz*dsz(-2)*(sx0(-2)+0.5*dsx(-2))
    sdz(n,2)  = wqz*dsz(-2)*(sx0(-1)+0.5*dsx(-1))
    sdz(n,3)  = wqz*dsz(-2)*(sx0(0)+0.5*dsx(0))
    sdz(n,4)  = wqz*dsz(-2)*(sx0(1)+0.5*dsx(1))
    sdz(n,5)  = wqz*dsz(-2)*(sx0(2)+0.5*dsx(2))
    sdz(n,6)  = wqz*dsz(-2)*(sx0(3)+0.5*dsx(3))
    sdz(n,9)  = wqz*dsz(-1)*(sx0(-2)+0.5*dsx(-2))
    sdz(n,9) = sdz(n,9)+sdz(n,8)
    sdz(n,10)  = wqz*dsz(-1)*(sx0(-1)+0.5*dsx(-1))
    sdz(n,10) = sdz(n,10)+sdz(n,9)
    sdz(n,11)  = wqz*dsz(-1)*(sx0(0)+0.5*dsx(0))
    sdz(n,11) = sdz(n,11)+sdz(n,10)
    sdz(n,12)  = wqz*dsz(-1)*(sx0(1)+0.5*dsx(1))
    sdz(n,12) = sdz(n,12)+sdz(n,11)
    sdz(n,13)  = wqz*dsz(-1)*(sx0(2)+0.5*dsx(2))
    sdz(n,13) = sdz(n,13)+sdz(n,12)
    sdz(n,14)  = wqz*dsz(-1)*(sx0(3)+0.5*dsx(3))
    sdz(n,14) = sdz(n,14)+sdz(n,13)
    sdz(n,17)  = wqz*dsz(0)*(sx0(-2)+0.5*dsx(-2))
    sdz(n,17) = sdz(n,17)+sdz(n,16)
    sdz(n,18)  = wqz*dsz(0)*(sx0(-1)+0.5*dsx(-1))
    sdz(n,18) = sdz(n,18)+sdz(n,17)
    sdz(n,19)  = wqz*dsz(0)*(sx0(0)+0.5*dsx(0))
    sdz(n,19) = sdz(n,19)+sdz(n,18)
    sdz(n,20)  = wqz*dsz(0)*(sx0(1)+0.5*dsx(1))
    sdz(n,20) = sdz(n,20)+sdz(n,19)
    sdz(n,21)  = wqz*dsz(0)*(sx0(2)+0.5*dsx(2))
    sdz(n,21) = sdz(n,21)+sdz(n,20)
    sdz(n,22)  = wqz*dsz(0)*(sx0(3)+0.5*dsx(3))
    sdz(n,22) = sdz(n,22)+sdz(n,21)
    sdz(n,25)  = wqz*dsz(1)*(sx0(-2)+0.5*dsx(-2))
    sdz(n,25) = sdz(n,25)+sdz(n,24)
    sdz(n,26)  = wqz*dsz(1)*(sx0(-1)+0.5*dsx(-1))
    sdz(n,26) = sdz(n,26)+sdz(n,25)
    sdz(n,27)  = wqz*dsz(1)*(sx0(0)+0.5*dsx(0))
    sdz(n,27) = sdz(n,27)+sdz(n,26)
    sdz(n,28)  = wqz*dsz(1)*(sx0(1)+0.5*dsx(1))
    sdz(n,28) = sdz(n,28)+sdz(n,27)
    sdz(n,29)  = wqz*dsz(1)*(sx0(2)+0.5*dsx(2))
    sdz(n,29) = sdz(n,29)+sdz(n,28)
    sdz(n,30)  = wqz*dsz(1)*(sx0(3)+0.5*dsx(3))
    sdz(n,30) = sdz(n,30)+sdz(n,29)
    sdz(n,33)  = wqz*dsz(2)*(sx0(-2)+0.5*dsx(-2))
    sdz(n,33) = sdz(n,33)+sdz(n,32)
    sdz(n,34)  = wqz*dsz(2)*(sx0(-1)+0.5*dsx(-1))
    sdz(n,34) = sdz(n,34)+sdz(n,33)
    sdz(n,35)  = wqz*dsz(2)*(sx0(0)+0.5*dsx(0))
    sdz(n,35) = sdz(n,35)+sdz(n,34)
    sdz(n,36)  = wqz*dsz(2)*(sx0(1)+0.5*dsx(1))
    sdz(n,36) = sdz(n,36)+sdz(n,35)
    sdz(n,37)  = wqz*dsz(2)*(sx0(2)+0.5*dsx(2))
    sdz(n,37) = sdz(n,37)+sdz(n,36)
    sdz(n,38)  = wqz*dsz(2)*(sx0(3)+0.5*dsx(3))
    sdz(n,38) = sdz(n,38)+sdz(n,37)

    !print*,'Particle:',nn,n,ip

    ENDDO
#if defined _OPENMP && _OPENMP>=201307
    !$OMP END SIMD 
#endif

    ! Add weights to nearest vertices
    DO n=1,MIN(LVECT,np-ip+1)
#if defined __INTEL_COMPILER 
            !DIR$ ASSUME_ALIGNED jxcells:64, jycells:64, jzcells:64
#elif defined __IBMBGQ__
            !IBM* ALIGN(64,jxcells, jycells, jzcells)
#endif 
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !$DIR SIMD 
#endif
      DO nv=1,8

        !print*,nv,ICELL(n,1)
      
        jxcells(nv,ICELL(n,1)) = jxcells(nv,ICELL(n,1)) + sdx(n,nv)
        jxcells(nv,ICELL(n,1)+ncx) = jxcells(nv,ICELL(n,1)+ncx) + sdx(n,nv+8)
        jxcells(nv,ICELL(n,1)+2*ncx) = jxcells(nv,ICELL(n,1)+2*ncx) + sdx(n,nv+16)
        jxcells(nv,ICELL(n,1)+3*ncx) = jxcells(nv,ICELL(n,1)+3*ncx) + sdx(n,nv+24)
        jxcells(nv,ICELL(n,1)+4*ncx) = jxcells(nv,ICELL(n,1)+4*ncx) + sdx(n,nv+32)
        jxcells(nv,ICELL(n,1)+5*ncx) = jxcells(nv,ICELL(n,1)+5*ncx) + sdx(n,nv+40)

        jycells(nv,ICELL(n,1)) = jycells(nv,ICELL(n,1)) + sdy(n,nv)
        jycells(nv,ICELL(n,1)+ncx) = jycells(nv,ICELL(n,1)+ncx) + sdy(n,nv+8)
        jycells(nv,ICELL(n,1)+2*ncx) = jycells(nv,ICELL(n,1)+2*ncx) + sdy(n,nv+16)
        jycells(nv,ICELL(n,1)+3*ncx) = jycells(nv,ICELL(n,1)+3*ncx) + sdy(n,nv+24)
        jycells(nv,ICELL(n,1)+4*ncx) = jycells(nv,ICELL(n,1)+4*ncx) + sdy(n,nv+32)
        jycells(nv,ICELL(n,1)+5*ncx) = jycells(nv,ICELL(n,1)+5*ncx) + sdy(n,nv+40)

        jzcells(nv,ICELL(n,1)) = jzcells(nv,ICELL(n,1)) + sdz(n,nv)
        jzcells(nv,ICELL(n,1)+ncx) = jzcells(nv,ICELL(n,1)+ncx) + sdz(n,nv+8)
        jzcells(nv,ICELL(n,1)+2*ncx) = jzcells(nv,ICELL(n,1)+2*ncx) + sdz(n,nv+16)
        jzcells(nv,ICELL(n,1)+3*ncx) = jzcells(nv,ICELL(n,1)+3*ncx) + sdz(n,nv+24)
        jzcells(nv,ICELL(n,1)+4*ncx) = jzcells(nv,ICELL(n,1)+4*ncx) + sdz(n,nv+32)
        
                                      
      ENDDO
#if defined _OPENMP && _OPENMP>=201307
         !$OMP END SIMD
#endif
    ENDDO
  ENDDO

  !print*,'Reduction of jxcells,jycells,jzcells in jx,jy,jz'

  ! Reduction of jxcells,jycells,jzcells in jx,jy,jz
  DO iz=1, ncz-2
#if defined _OPENMP && _OPENMP>=201307
      !$OMP SIMD 
#elif defined __IBMBGQ__
      !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER 
      !$DIR SIMD 
#endif
      DO ix=1,ncx-8 !! VECTOR (take ncx multiple of vector length)
      
        ic=ix+(iz-1)*ncx
        igrid =orig+ix+(iz-1)*nnx

        !print*
        !print*,'ic',ic
        !print*,'ix/ncx',ix,ncx,'iz/ncz',iz,ncz
        !print*,'igrid',igrid
        
        ! jx
        jx(igrid)=jx(igrid)+jxcells(1,ic)
        jx(igrid+1)=jx(igrid+1)+jxcells(2,ic)
        jx(igrid+2)=jx(igrid+2)+jxcells(3,ic)
        jx(igrid+3)=jx(igrid+3)+jxcells(4,ic)
        jx(igrid+4)=jx(igrid+4)+jxcells(5,ic)
        jx(igrid+5)=jx(igrid+5)+jxcells(6,ic)
        jx(igrid+6)=jx(igrid+6)+jxcells(7,ic)
        jx(igrid+7)=jx(igrid+7)+jxcells(8,ic)
        ! jy
        jy(igrid)  =jy(igrid)  +jycells(1,ic)
        jy(igrid+1)=jy(igrid+1)+jycells(2,ic)
        jy(igrid+2)=jy(igrid+2)+jycells(3,ic)
        jy(igrid+3)=jy(igrid+3)+jycells(4,ic)
        jy(igrid+4)=jy(igrid+4)+jycells(5,ic)
        jy(igrid+5)=jy(igrid+5)+jycells(6,ic)
        jy(igrid+6)=jy(igrid+6)+jycells(7,ic)
        jy(igrid+7)=jy(igrid+7)+jycells(8,ic)
        ! jz
        jz(igrid)  =jz(igrid)  +jzcells(1,ic)
        jz(igrid+1)=jz(igrid+1)+jzcells(2,ic)
        jz(igrid+2)=jz(igrid+2)+jzcells(3,ic)
        jz(igrid+3)=jz(igrid+3)+jzcells(4,ic)
        jz(igrid+4)=jz(igrid+4)+jzcells(5,ic)
        jz(igrid+5)=jz(igrid+5)+jzcells(6,ic)
        jz(igrid+6)=jz(igrid+6)+jzcells(7,ic)
        jz(igrid+7)=jz(igrid+7)+jzcells(8,ic)
                
      END DO
#if defined _OPENMP && _OPENMP>=201307
         !$OMP END SIMD
#endif
  ENDDO  
  
  DEALLOCATE(jxcells,jycells,jzcells)
  DEALLOCATE(sdx,sdy,sdz)
  
  !print*,'finished'
  
End subroutine pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3
#endif

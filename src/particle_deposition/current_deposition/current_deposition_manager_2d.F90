! ______________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)  
! 2016, The Regents of the University of California, through Lawrence Berkeley 
! National Laboratory (subject to receipt of any required approvals from the 
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software, 
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy 
! and the U.S. Government consequently retains certain rights. As such, the U.S. 
! Government has been granted for itself and others acting on its behalf a  
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to 
! reproduce, distribute copies to the public, prepare derivative works, and 
! perform publicly and display publicly, and to permit other to do so. 
! 
! CURRENT_DEPOSITION_MANAGER_2D.F90
!
! Developers
! Henri Vincenti,
! Mathieu Lobet
!
! Description:
! This file contains subroutines for managing the current deposition in 2D.
!
! List of subroutines
! - pxrdepose_currents_on_grid_jxjyjz_2d
!
! ______________________________________________________________________________

! ______________________________________________________________________________
!> @brief
!> Main subroutine for the current deposition called in submain in 2d.

!> @details
!> This subroutine determines which model and algorithm to use depending on the parameter
!> currdepo and the interpolation order.

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2016


SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_2d
! ______________________________________________________________________________
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



! ______________________________________________________________________________
!> @brief
!> Deposit current in each tile with Esirkepov method in 2D
!
!> @details
!> This subroutine is called from Fortran main program and contains an interface argument 
!> OpenMP version. Avoids conflict while reducing tile currents in the global 
!> current array. 

!> @author
!> Henri Vincenti
!> Mathieu Lobet

!> @date
!> Creation 2016
SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp(curr_depo_sub,jxg,jyg,jzg,&
                nxx,nyy,nzz,nxjguard,nyjguard,nzjguard,noxx,noyy,nozz,dxx,dyy,dzz,dtt,lvect)
! ______________________________________________________________________________
  USE particles
  USE constants
  USE tiling
  USE omp_lib
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
 
END SUBROUTINE pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp


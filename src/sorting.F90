! ______________________________________________________________________________
! 
! SORTING.F90
! Contains sorting algorithms for the particles
! 
! List of subroutines:
! - pxr_particle_sorting
! - particle_sorting_sub

MODULE sorting

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE pxr_particle_sorting
    ! __________________________________________________________________________________
    ! 
    ! Main subroutine called to sort the particles in the Fortran PIC loop
    ! __________________________________________________________________________________
    USE tiling
    USE shared_data
    USE constants
    USE particles
    USE time_stat
    USE params
    IMPLICIT NONE

    !LOGICAL(idp) :: verbose=.TRUE.
    REAL(num) :: tdeb, tend
    
    IF ((sorting_activated.gt.0)) THEN
    
      tdeb=MPI_WTIME()
    
      CALL particle_sorting_sub
        
      tend=MPI_WTIME()
      localtimes(10) = localtimes(10) + (tend-tdeb)    
  
    ENDIF
  
  END SUBROUTINE

  SUBROUTINE particle_sorting_sub
    ! __________________________________________________________________________________
    ! 
    ! General subroutine for the particle sorting, used in Python PIC loop
    ! __________________________________________________________________________________
    USE tiling
    USE shared_data
    USE constants
    USE particles
    USE time_stat
    USE params
    IMPLICIT NONE

    INTEGER(idp) :: ispecies, ix, iy, iz, count

    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile
    TYPE(grid_tile), POINTER :: currg
    INTEGER(idp) :: nxc, nyc, nzc, np,ip
    INTEGER(idp) :: nxjg,nyjg,nzjg
    LOGICAL(idp) :: isgathered=.FALSE.
    LOGICAL(idp) :: verbose=.TRUE.
    REAL(num) :: sxmin,symin,szmin
    REAL(num) :: tdeb, tend
    
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
    !$OMP SHARED(ntilex,ntiley,ntilez,nspecies,species_parray,aofgrid_tiles,dx,dy,dz,it,rank, &
    !$OMP sorting_shiftx,sorting_shifty,sorting_shiftz,sorting_dx,sorting_dy,sorting_dz,sorting_verbose) &
    !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile,currg,count, &
    !$OMP nxc,nyc,nzc,nxjg,nyjg,nzjg,isgathered,sxmin,symin,szmin)
    DO iz=1, ntilez ! LOOP ON TILES
      DO iy=1, ntiley
        DO ix=1, ntilex
					curr=>species_parray(1)
			    curr_tile=>curr%array_of_tiles(ix,iy,iz)
			    nxjg=curr_tile%nxg_tile
			    nyjg=curr_tile%nyg_tile
			    nzjg=curr_tile%nzg_tile
			    nxc=curr_tile%nx_cells_tile
			    nyc=curr_tile%ny_cells_tile
			    nzc=curr_tile%nz_cells_tile
			    isgathered=.FALSE.
			    
			    
			    ! Loop over the species
			    DO ispecies=1, nspecies
			        curr=>species_parray(ispecies)
              curr_tile=>curr%array_of_tiles(ix,iy,iz)
              count=curr_tile%np_tile(1)
              IF (count .GT. 0) isgathered=.TRUE.                    	
          END DO
          
          !if (rank.eq.0) print*,ix,iy,iz
          
          IF (isgathered) THEN
            currg=>aofgrid_tiles(ix,iy,iz)

			      ! Loop over the species            
				    DO ispecies=1, nspecies

					    curr=>species_parray(ispecies)
				      
				      ! If the sorting period > 0 and the current iteration corresponds to a multiple of the period
				      IF ((it.ge.curr%sorting_start).AND.(curr%sorting_period.gt.0).AND.(MOD(it,curr%sorting_period).eq.0)) THEN

              
                IF ((sorting_verbose).and.(rank.eq.0).and. &
                    (iz.eq.1).and.(iy.eq.1).and.(ix.eq.1)) WRITE(0,*) 'Particle sorting, species',ispecies
 
					      ! - Get current tile properties
					      ! - Init current tile variables				    
					      curr_tile=>curr%array_of_tiles(ix,iy,iz)
					      count=curr_tile%np_tile(1)
					      
					      
					      IF (count .EQ. 0) CYCLE
       
                ! Sorting algorithm inside the tiles
                
                !if (rank.eq.0) print*, curr_tile%x_tile_min, curr_tile%y_tile_min, curr_tile%z_tile_min
                !if (rank.eq.0) print*, curr_tile%x_tile_max, curr_tile%y_tile_max, curr_tile%z_tile_max 
                
                sxmin = curr_tile%x_tile_min + sorting_shiftx
                symin = curr_tile%y_tile_min + sorting_shifty
                szmin = curr_tile%z_tile_min + sorting_shiftz
                
                CALL pxr_particle_bin_sorting(count,curr_tile%part_x,curr_tile%part_y,curr_tile%part_z, &
                curr_tile%part_ux,curr_tile%part_uy,curr_tile%part_uz,curr_tile%part_gaminv,&
                curr_tile%pid, wpid, sxmin,symin,szmin, &
                curr_tile%x_tile_max,curr_tile%y_tile_max,curr_tile%z_tile_max, &
                sorting_dx, sorting_dy, sorting_dz)
     
              ENDIF
     
				    END DO! END LOOP ON SPECIES
			    ENDIF
        END DO
      END DO
    END DO! END LOOP ON TILES
    !$OMP END PARALLEL DO
  
  END SUBROUTINE particle_sorting_sub


  SUBROUTINE pxr_particle_bin_sorting(np2,xp,yp,zp,ux,uy,uz,gam,pid,wpid,xmin2,ymin2,zmin2,xmax2,ymax2,zmax2,dxf,dyf,dzf)
    ! __________________________________________________________________________________
    ! 
    ! Particle bin sorting algorithm
    ! This subroutine uses a bin sorting algorithm to sort particles (including their property arrays)
    ! 
    ! np: number of particles
    ! xp,yp,zp: particle positions
    ! xp,yp,zp: particle momenta
    ! xmin,ymin,zmin: minimum point position on the local grid 
    ! xmax,ymax,zmax: maximum point position on the local grid           
    ! dxf,dyf,dzf: bin space steps
    ! ___________________________________________________________________________________
    USE constants
    implicit none

    integer(idp) :: ip,np2
    integer(idp) :: j,k,ic,nbhc
    integer(idp) :: ix,iy,iz
    integer(idp) :: nx2,ny2,nz2
    integer(idp) :: nx3,ny3,nz3
    
    integer(idp) :: wpid
    
    real(num) :: dxi,dyi,dzi
    real(num) :: dxf,dyf,dzf           
    real(num) :: x2,y2,z2
    real(num) :: xmin2,ymin2,zmin2
    real(num) :: xmax2,ymax2,zmax2
            
    real(num), dimension(np2), intent(inout)      :: xp,yp,zp
    real(num), dimension(np2), intent(inout)      :: ux,uy,uz
    real(num), dimension(np2), intent(inout)      :: gam
        
    REAL(num), DIMENSION(np2,1), intent(inout)      :: pid
        
    real(num), dimension(np2)      :: xps,yps,zps     
    real(num), dimension(np2)      :: uxs,uys,uzs 
    real(num), dimension(np2)      :: gams 
    real(num), dimension(np2,wpid)      :: pids
            
    integer(idp), dimension(np2)               :: hcnb        ! Cell number
    integer(idp), dimension(:),allocatable    :: piihc       ! Particle indexes in the grid
    integer(idp), dimension(:),allocatable    :: nbppc       ! Number of particles per cells

    ! Bin sizes
    dxi = 1./dxf
    dyi = 1./dyf
    dzi = 1./dzf
    
    nx3 = ceiling((xmax2-xmin2)*dxi)
    ny3 = ceiling((ymax2-ymin2)*dyi)
    nz3 = ceiling((zmax2-zmin2)*dzi)
    
    ! Number of bins
    nbhc = nx3*ny3*nz3
  
    allocate(piihc(nbhc))
    allocate(nbppc(nbhc))
      
    hcnb = 0
    piihc = 0
    nbppc = 0    
  
    ! Counting sort
    ! The criteria is the position in term of bin position
    !if (rank.eq.0) print*, 'Counting sort, phase 1'
    DO ip=1,np2
    
      !if (rank.eq.0) print*, 'Ip',ip,np2
      !if (rank.eq.0) print*, 'nx',nx3,ny3,nz3
      
      x2 = (xp(ip)-xmin2)*dxi
      y2 = (yp(ip)-ymin2)*dyi
      z2 = (zp(ip)-zmin2)*dzi
    
      ix = floor(x2)
      iy = floor(y2)
      iz = floor(z2)
      
      ! Bin id
      hcnb(ip) = iz*nx3*ny3 + iy*nx3 + ix+1
      !IF (hcnb(ip) > nbhc) THEN
      !  if (rank.eq.0) print*, 'Bin id',hcnb(ip),nbhc
      !  if (rank.eq.0) print*, 'Particle ix,iy,iz',ix,iy,iz
      !  if (rank.eq.0) print*, 'Particle x,y,z',xp(ip),yp(ip),zp(ip)
      !  if (rank.eq.0) print*, 'Particle x2,y2,z2',x2,y2,z2
      !  if (rank.eq.0) print*, 'Particle dx,dy,dz',dx,dy,dz  
      !  if (rank.eq.0) print*, 'Particle nx,ny,nz',nx3,ny3,nz3          
      !ENDIF

      ! We count the number of particles in each bin
      nbppc(hcnb(ip)) = nbppc(hcnb(ip))+1  
      !if (rank.eq.0) print*, 'nbppc(hcnb(ip)) = nbppc(hcnb(ip))+1',nbppc(hcnb(ip))
      
      !write(0,'( " Particle",X,I4)') ip 
    ENDDO
    
    ! Determine particle indexes in the bin grid
    !if (rank.eq.0) print*, 'Counting sort, phase 2'    
    k=0
    !!!$OMP SIMD reduction(+:k)
    DO ic = 1,nbhc
       piihc(ic) = k
       k = k+nbppc(ic)
       !write(0,'( " Cell",X,I4)') ic
    END DO
    
    ! Sorting of the particles including their properties
    !if (rank.eq.0) print*, 'Counting sort, phase 3'
    DO ip=1,np2

      !write(0,'( " Particle",X,I4)') ip   
      k = hcnb(ip)      
      piihc(k) = piihc(k) + 1
      
      !IF (k > nbhc) THEN
      !  print*,'k>nbhc'
      !ENDIF
      !IF (piihc(k)>np2) THEN
      !  print*,'piihc(k)>np2'
      !ENDIF      
      
      !write(0,*) k,np,piihc(k),nbhc

      pids(piihc(k),:) = pid(ip,:)
      
      xps(piihc(k)) = xp(ip)
      yps(piihc(k)) = yp(ip)
      zps(piihc(k)) = zp(ip)

      uxs(piihc(k)) = ux(ip)
      uys(piihc(k)) = uy(ip)
      uzs(piihc(k)) = uz(ip)
      
      gams(piihc(k)) = gam(ip)
    
    END DO
  
    ! Copy back to the original arrays

    pid = pids
    
    xp = xps
    yp = yps
    zp = zps

    ux = uxs
    uy = uys
    uz = uzs
    
    gam = gams
  
    deallocate(piihc,nbppc)
  
  END SUBROUTINE
  
END MODULE sorting

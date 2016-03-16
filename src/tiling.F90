MODULE tiling
!!!! --- This module contains useful diagnostics to test code correctness
    USE constants
    USE particles
    USE shared_data
    USE fields
    USE params
    IMPLICIT NONE

CONTAINS

	SUBROUTINE set_tile_split()
	IMPLICIT NONE 
		! Set tile split for species arrays
		CALL set_tile_split_for_species(species_parray,nspecies,ntilex,ntiley,ntilez,nx_grid,ny_grid,nz_grid, &
    									  x_min_local,y_min_local,z_min_local,x_max_local,y_max_local,z_max_local)
    								
		! ALLOCATE grid tile arrays
		ALLOCATE(aofgrid_tiles(ntilex,ntiley,ntilez))
	END SUBROUTINE set_tile_split


    !!! --- Set particle tile split in space
    SUBROUTINE set_tile_split_for_species(species_array,nspec,ntx,nty,ntz,nxgrid,nygrid,nzgrid, &
    									  xminlocal,yminlocal,zminlocal,xmaxlocal,ymaxlocal,zmaxlocal) !#do not parse 
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN) :: nspec, nxgrid, nygrid, nzgrid
        INTEGER(idp), INTENT(IN OUT) ::  ntx, nty, ntz
        REAL(num), INTENT(IN) :: xminlocal,yminlocal,zminlocal,xmaxlocal,ymaxlocal,zmaxlocal
        TYPE(particle_species), INTENT(IN OUT), TARGET, DIMENSION(nspec) :: species_array
        INTEGER(idp) :: ix, iy, iz, ispecies
        INTEGER(idp) :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
        INTEGER(idp) :: nx0_last_tile, ny0_last_tile, nz0_last_tile
        TYPE(particle_species), POINTER :: curr_sp
        TYPE(particle_tile), POINTER :: curr

        ! Tile-split
        nx0_grid_tile = nxgrid / ntx
        ny0_grid_tile = nygrid / nty
        nz0_grid_tile = nzgrid / ntz

        ! Some sanity check
        IF (nx0_grid_tile .LT. 4) THEN
            IF (rank .EQ. 0) PRINT *, "number of tiles in X is to high, settting back to default value 1"
            ntx=1
        END IF
        IF (ny0_grid_tile .LT. 4) THEN
            IF (rank .EQ. 0) PRINT *, "number of tiles in Y is to high, setting back to default value 1"
            nty=1
        END IF
        IF (nz0_grid_tile .LT. 4) THEN
            IF (rank .EQ. 0) PRINT *, "number of tiles in Z is to high, setting back to default value 1"
            ntz=1
        END IF
        !-- N.B: If the number of grid points cannot be equally divided between
        !-- tiles then give remaining points to last tile in each dimension
        nx0_last_tile= nx0_grid_tile+(nxgrid-nx0_grid_tile*ntx)
        ny0_last_tile= ny0_grid_tile+(nygrid-ny0_grid_tile*nty)
        nz0_last_tile= nz0_grid_tile+(nzgrid-nz0_grid_tile*ntz)


        !- Allocate object array of tiles for particles
        DO ispecies =1, nspecies
            curr_sp => species_array(ispecies)
            IF (.NOT. curr_sp%l_arrayoftiles_allocated) THEN
                ALLOCATE(curr_sp%array_of_tiles(ntx,nty,ntz))
                ALLOCATE(curr_sp%are_tiles_reallocated(ntx,nty,ntz))
                curr_sp%are_tiles_reallocated= 0
                curr_sp%l_arrayoftiles_allocated = .TRUE.
            END IF
            ! Sets tile spatial extents for current species
            DO iz=1, ntz
                DO iy=1,nty
                    DO ix=1,ntx
                        curr=> curr_sp%array_of_tiles(ix,iy,iz)
                        !------------- X- DIRECTION
						! FIRST TILE in X DIRECTION
                        IF (ix .EQ. 1) THEN
							curr%subdomain_bound = .TRUE.
							curr%nx_grid_tile=nx0_grid_tile
							curr%nx_cells_tile=curr%nx_grid_tile-1
							curr%x_grid_tile_min=xminlocal
							curr%x_grid_tile_max=curr%x_grid_tile_min+(nx0_grid_tile-1)*dx
							curr%x_tile_min=curr%x_grid_tile_min
							curr%x_tile_max=curr%x_grid_tile_max+dx/2.0_num
							curr%nx_tile_min = (ix-1)*nx0_grid_tile
							curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
						ENDIF
                        IF ((ix .LT. ntx) .AND. (ix .GT. 1)) THEN
                            curr%nx_grid_tile=nx0_grid_tile
                            curr%nx_cells_tile=curr%nx_grid_tile-1
                            curr%x_grid_tile_min=xminlocal+(ix-1)*nx0_grid_tile*dx
                            curr%x_grid_tile_max=curr%x_grid_tile_min+curr%nx_cells_tile*dx
                            curr%nx_tile_min = (ix-1)*nx0_grid_tile
                            curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
							curr%x_tile_min= curr%x_grid_tile_min-dx/2.0_num
							curr%x_tile_max= curr%x_grid_tile_max+dx/2.0_num
                        END IF
						! LAST TILE in X DIRECTION
						IF (ix .EQ. ntx) THEN
                            curr%subdomain_bound= .TRUE.
                            curr%nx_grid_tile=nx0_last_tile
                            curr%nx_cells_tile=curr%nx_grid_tile-1
                            curr%x_grid_tile_min=xminlocal+(ix-1)*nx0_grid_tile*dx
                            curr%x_grid_tile_max=curr%x_grid_tile_min+curr%nx_cells_tile*dx
                            curr%nx_tile_min = (ix-1)*nx0_grid_tile
                            curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
							curr%x_tile_min= curr%x_grid_tile_min-dx/2.0_num
							curr%x_tile_max= xmaxlocal
                        ENDIF
                        !------------- Y- DIRECTION
						! FIRST TILE in Y DIRECTION
                        IF (iy .EQ. 1) THEN
							curr%subdomain_bound = .TRUE.
							curr%ny_grid_tile=ny0_grid_tile
							curr%ny_cells_tile=curr%ny_grid_tile-1
							curr%y_grid_tile_min=yminlocal
							curr%y_grid_tile_max=curr%y_grid_tile_min+(ny0_grid_tile-1)*dy
							curr%y_tile_min=curr%y_grid_tile_min
							curr%y_tile_max=curr%y_grid_tile_max+dy/2.0_num
							curr%ny_tile_min = (iy-1)*ny0_grid_tile
							curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
						ENDIF
                        IF ((iy .LT. nty) .AND. (iy .GT. 1)) THEN
                            curr%ny_grid_tile=ny0_grid_tile
                            curr%ny_cells_tile=curr%ny_grid_tile-1
                            curr%y_grid_tile_min=yminlocal+(iy-1)*ny0_grid_tile*dy
                            curr%y_grid_tile_max=curr%y_grid_tile_min+curr%ny_cells_tile*dy
                            curr%ny_tile_min = (iy-1)*ny0_grid_tile
                            curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
							curr%y_tile_min= curr%y_grid_tile_min-dy/2.0_num
							curr%y_tile_max= curr%y_grid_tile_max+dy/2.0_num
                        END IF
						! LAST TILE in Y DIRECTION
						IF (iy .EQ. nty) THEN
                            curr%subdomain_bound= .TRUE.
                            curr%ny_grid_tile=ny0_last_tile
                            curr%ny_cells_tile=curr%ny_grid_tile-1
                            curr%y_grid_tile_min=yminlocal+(iy-1)*ny0_grid_tile*dy
                            curr%y_grid_tile_max=curr%y_grid_tile_min+curr%ny_cells_tile*dy
                            curr%ny_tile_min = (iy-1)*ny0_grid_tile
                            curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
							curr%y_tile_min= curr%y_grid_tile_min-dy/2.0_num
							curr%y_tile_max= ymaxlocal
                        ENDIF
                        !------------- Z- DIRECTION
						! FIRST TILE in Z DIRECTION
                        IF (iz .EQ. 1) THEN
							curr%subdomain_bound = .TRUE.
							curr%nz_grid_tile=nz0_grid_tile
							curr%nz_cells_tile=curr%nz_grid_tile-1
							curr%z_grid_tile_min=zminlocal
							curr%z_grid_tile_max=curr%z_grid_tile_min+(nz0_grid_tile-1)*dz
							curr%z_tile_min=curr%z_grid_tile_min
							curr%z_tile_max=curr%z_grid_tile_max+dz/2.0_num
							curr%nz_tile_min = (iz-1)*nz0_grid_tile
							curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
						ENDIF
                        IF ((iz .LT. ntz) .AND. (iz .GT. 1)) THEN
                            curr%nz_grid_tile=nz0_grid_tile
                            curr%nz_cells_tile=curr%nz_grid_tile-1
                            curr%z_grid_tile_min=zminlocal+(iz-1)*nz0_grid_tile*dz
                            curr%z_grid_tile_max=curr%z_grid_tile_min+curr%nz_cells_tile*dz
                            curr%nz_tile_min = (iz-1)*nz0_grid_tile
                            curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
							curr%z_tile_min= curr%z_grid_tile_min-dz/2.0_num
							curr%z_tile_max= curr%z_grid_tile_max+dz/2.0_num
                        END IF
						! LAST TILE in Z DIRECTION
						IF (iz .EQ. ntz) THEN
                            curr%subdomain_bound= .TRUE.
                            curr%nz_grid_tile=nz0_last_tile
                            curr%nz_cells_tile=curr%nz_grid_tile-1
                            curr%z_grid_tile_min=zminlocal+(iz-1)*nz0_grid_tile*dz
                            curr%z_grid_tile_max=curr%z_grid_tile_min+curr%nz_cells_tile*dz
                            curr%nz_tile_min = (iz-1)*nz0_grid_tile
                            curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
							curr%z_tile_min= curr%z_grid_tile_min-dz/2.0_num
							curr%z_tile_max= zmaxlocal
                        ENDIF
                    END DO
                END DO
            END DO
        END DO ! END DO SPECIES
    END SUBROUTINE set_tile_split_for_species

    !!! --- Add particle to array of tiles
    SUBROUTINE add_particle_to_species(currsp, partx, party, partz, &
               partux, partuy, partuz, gaminv, partw)
        IMPLICIT NONE
        REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
        TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
        TYPE(particle_tile), POINTER :: curr
        INTEGER(idp) :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile, nptile
        INTEGER(idp) :: ixtile, iytile, iztile


        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = currsp%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = currsp%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = currsp%array_of_tiles(1,1,1)%nz_grid_tile

        ! Get particle index in array of tile
		ixtile = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex)
		iytile = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy))+1,ntiley)
		iztile = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez)

        ! Point to current tile arr_of_tiles(ixtile,iytile,iztile)
        !curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
        
        CALL add_particle_at_tile(currsp,ixtile,iytile,iztile, partx, party, partz, &
                partux, partuy, partuz, gaminv, partw)

        ! Update total number of particle species
        currsp%species_npart=currsp%species_npart+1
    END SUBROUTINE add_particle_to_species

    SUBROUTINE add_particle_at_tile(currsp, ixt, iyt, izt, partx, party, partz, &
                partux, partuy, partuz, gaminv, partw)
        IMPLICIT NONE
        INTEGER(idp) :: count, nmax, ixt, iyt, izt
        REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv, partw
        TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
        TYPE(particle_tile), POINTER :: curr
        
        curr=>currsp%array_of_tiles(ixt,iyt,izt)
        ! If no particles in tile, allocate particle arrays
        IF (.NOT. curr%l_arrays_allocated) THEN
            CALL allocate_tile_arrays(curr)
        ENDIF

        ! Sanity check for max number of particles in tile
        count = curr%np_tile(1)+1
        nmax  = curr%npmax_tile
        IF (count .GT. nmax) THEN
        ! Resize particle tile arrays if tile is full
        	currsp%are_tiles_reallocated(ixt,iyt,izt)=1
            CALL resize_particle_arrays(curr, nmax, NINT(resize_factor*nmax+1,idp))
        ENDIF
        ! Finally, add particle to tile
        curr%np_tile(1)=count
        curr%part_x(count)  = partx
        curr%part_y(count)  = party
        curr%part_z(count)  = partz
        curr%part_ux(count) = partux
        curr%part_uy(count) = partuy
        curr%part_uz(count) = partuz
        curr%part_gaminv(count) = gaminv
        curr%pid(count,wpid) = partw
        curr%part_ex(count)  = 0._num
        curr%part_ey(count)  = 0._num
        curr%part_ez(count)  = 0._num
        curr%part_bx(count)  = 0._num
        curr%part_by(count)  = 0._num
        curr%part_bz(count)  = 0._num
    END SUBROUTINE add_particle_at_tile

    !!! --- Remove particles from tile using a mask variable
    !!! --- This technique avoids packing or reallocating arrays
    SUBROUTINE rm_particles_from_species(currsp, curr, mask)
        TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        LOGICAL(idp), DIMENSION (:), INTENT(IN) :: mask
        INTEGER(idp) :: ninit, i
        ninit= curr%np_tile(1)
        DO i = ninit,1,-1
            IF (.NOT. mask(i)) THEN
                CALL rm_particle_at_tile(curr,i)
                currsp%species_npart=currsp%species_npart-1
            ENDIF
        ENDDO
    END SUBROUTINE rm_particles_from_species

    SUBROUTINE rm_particle_at_tile(curr, index)
        IMPLICIT NONE
        INTEGER(idp) :: index
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        IF (index .EQ. curr%np_tile(1)) THEN
            ! If particle i is last element
            ! Simply decreases particle number
            curr%np_tile(1)=curr%np_tile(1)-1
        ELSE
            ! Particle i replaced by last element
            ! Particle number is decreased
            curr%part_x(index)=curr%part_x(curr%np_tile(1))
            curr%part_y(index)=curr%part_y(curr%np_tile(1))
            curr%part_z(index)=curr%part_z(curr%np_tile(1))
            curr%part_ux(index)=curr%part_ux(curr%np_tile(1))
            curr%part_uy(index)=curr%part_uy(curr%np_tile(1))
            curr%part_uz(index)=curr%part_uz(curr%np_tile(1))
            curr%part_gaminv(index)=curr%part_gaminv(curr%np_tile(1))
            curr%pid(index,wpid)=curr%pid(curr%np_tile(1),wpid)
            curr%np_tile=curr%np_tile(1)-1
        END IF
    END SUBROUTINE rm_particle_at_tile

    SUBROUTINE allocate_tile_arrays(curr_tile)
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr_tile
        INTEGER(idp) :: nmax, nxc, nyc, nzc
        INTEGER(idp) :: nxjg,nyjg,nzjg
        ! ALLOCATE PARTICLE ARRAYS
        nmax = curr_tile%npmax_tile
        ALLOCATE(curr_tile%part_x(1:nmax), curr_tile%part_y(1:nmax),    	   &
                 curr_tile%part_z(1:nmax), curr_tile%part_ux(1:nmax),          &
                 curr_tile%part_uy(1:nmax), curr_tile%part_uz(1:nmax),         &
                 curr_tile%pid(1:nmax,1:npid), curr_tile%part_ex(1:nmax),      &
                 curr_tile%part_ey(1:nmax), curr_tile%part_ez(1:nmax),         &
                 curr_tile%part_bx(1:nmax), curr_tile%part_by(1:nmax),         &
                 curr_tile%part_bz(1:nmax),curr_tile%part_gaminv(1:nmax))
        curr_tile%l_arrays_allocated = .TRUE.

    END SUBROUTINE allocate_tile_arrays

	SUBROUTINE init_tile_arrays()
	IMPLICIT NONE 
	
		CALL init_tile_arrays_for_species(nspecies, species_parray, aofgrid_tiles, ntilex, ntiley, ntilez)
	
	END SUBROUTINE init_tile_arrays




    SUBROUTINE init_tile_arrays_for_species(nspec, species_array, aofgtiles, ntx, nty, ntz)
        IMPLICIT NONE
        TYPE(grid_tile), DIMENSION(ntx,nty,ntz), INTENT(IN OUT) :: aofgtiles
        TYPE(particle_species), DIMENSION(nspec), TARGET, INTENT(IN OUT) :: species_array
        INTEGER(idp), INTENT(IN) :: nspec, ntx, nty, ntz
        INTEGER(idp) :: ispecies, ix, iy, iz
        INTEGER(idp) :: n1, n2, n3, ng1, ng2, ng3
        TYPE(particle_tile), POINTER :: curr_tile
        TYPE(particle_species), POINTER :: curr
        
        ! Allocate particle tile arrays 
        DO ispecies=1,nspec ! LOOP ON SPECIES
            curr=>species_array(ispecies)
            curr%species_npart=0
            DO iz=1, ntz! LOOP ON TILES
                DO iy=1, nty
                    DO ix=1,ntx
                        curr_tile=>curr%array_of_tiles(ix,iy,iz)
                        ! - Max size of particle arrays of current ile
                        n1=curr_tile%nx_cells_tile
                        n2=curr_tile%ny_cells_tile
                        n3=curr_tile%nz_cells_tile
                        curr_tile%npmax_tile=n1*n2*n3*curr%nppcell
                        curr_tile%np_tile(1)=0
                        IF ((ix .GT. 1) .AND. (ix .LT. ntx)) THEN
                        	curr_tile%nxg_tile=MAX(nox,2)
                        ELSE
                        	curr_tile%nxg_tile=nxjguards
                        END IF
                        IF ((iy .GT. 1) .AND. (iy .LT. nty)) THEN
                        	curr_tile%nyg_tile=MAX(noy,2)
                        ELSE
                        	curr_tile%nyg_tile=nyjguards
                        END IF
                        IF ((iz .GT. 1) .AND. (iz .LT. ntz)) THEN
                        	curr_tile%nzg_tile=MAX(noz,2)
                        ELSE
                        	curr_tile%nzg_tile=nzjguards
                        END IF                        
                        ! - Allocate arrays of current tile
                       CALL allocate_tile_arrays(curr_tile)
                    END DO
                END DO
            END DO
        END DO
        

        ! Init partile tile arrays in parallel - first touch policy
        ! - Init array of current tile
        ! - For some reason, don't set all values to zero?????
        ! - Have to set it manually for each element through
        ! - a DO loop see add_particle_at_tile
        !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) &
        !$OMP SHARED(species_array, ntx, nty, ntz, nspec) &
        !$OMP PRIVATE(ix,iy,iz,ispecies,curr,curr_tile)
        DO iz=1, ntz ! LOOP ON TILES
            DO iy=1, nty
                DO ix=1, ntx
                    DO ispecies=1, nspec ! LOOP ON SPECIES
                        curr=>species_array(ispecies)
                        curr_tile=>curr%array_of_tiles(ix,iy,iz)
                        !!! --- Init tile arrays
                        curr_tile%part_x=0.0_num
                        curr_tile%part_y=0.0_num
                        curr_tile%part_z=0.0_num
                        curr_tile%part_ux=0.0_num
                        curr_tile%part_uy=0.0_num
                        curr_tile%part_uz=0.0_num
                        curr_tile%part_gaminv=0.0_num
                        curr_tile%part_ex=0.0_num
                        curr_tile%part_ey=0.0_num
                        curr_tile%part_ez=0.0_num
                        curr_tile%part_bx=0.0_num
                        curr_tile%part_by=0.0_num
                        curr_tile%part_bz=0.0_num
                        curr_tile%pid=0.0_num
                    END DO! END LOOP ON SPECIES
                END DO
            END DO
        END DO! END LOOP ON TILES
        !$OMP END PARALLEL DO
        ! Allocate grid tile arrays
        curr=>species_array(1)
		DO iz=1, ntz ! LOOP ON TILES
            DO iy=1, nty
                DO ix=1, ntx
                   		curr_tile=>curr%array_of_tiles(ix,iy,iz)
                        ! - Max size of particle arrays of current ile
                        n1=curr_tile%nx_cells_tile
                        n2=curr_tile%ny_cells_tile
                        n3=curr_tile%nz_cells_tile
                        ng1=curr_tile%nxg_tile
                        ng2=curr_tile%nyg_tile
                        ng3=curr_tile%nzg_tile
                        ALLOCATE(aofgtiles(ix,iy,iz)%extile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%eytile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%eztile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%bxtile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%bytile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%bztile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%jxtile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%jytile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%jztile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                        ALLOCATE(aofgtiles(ix,iy,iz)%rhotile(-ng1:n1+ng1,-ng2:n2+ng2,-ng3:n3+ng3))
                END DO
            END DO
        END DO! END LOOP ON TILES
    END SUBROUTINE init_tile_arrays_for_species

    SUBROUTINE load_particles
        IMPLICIT NONE
        TYPE(particle_species), POINTER :: curr
        INTEGER(idp) :: ispecies, l, k, j, ipart
        INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
        REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
        REAL(num) :: phi, th, v, usq, clightsq
        INTEGER(KIND=4) :: err, npart
        REAL(num), DIMENSION(6) :: rng=0_num
		clightsq=1/clight**2
        !!! --- Sets-up particle space distribution (homogeneous case - default)
        IF (pdistr .EQ. 1) THEN
            DO ispecies=1,nspecies
                curr=>species_parray(ispecies)
                jmin = NINT(MAX(curr%x_min-x_min_local,0.0_num)/dx)
                jmax = NINT(MIN(curr%x_max-x_min_local,x_max_local-x_min_local)/dx)
                kmin = NINT(MAX(curr%y_min-y_min_local,0.0_num)/dy)
                kmax = NINT(MIN(curr%y_max-y_min_local,y_max_local-y_min_local)/dy)
                lmin = NINT(MAX(curr%z_min-z_min_local,0.0_num)/dz)
                lmax = NINT(MIN(curr%z_max-z_min_local,z_max_local-z_min_local)/dz)
                DO l=lmin,lmax-1
                    DO k=kmin,kmax-1
                        DO j=jmin,jmax-1
                            DO ipart=1,curr%nppcell
                                ! Sets positions and weight
                                partx = x_min_local+j*dx+dx/curr%nppcell*(ipart-1)
                                party = y_min_local+k*dy+dy/curr%nppcell*(ipart-1)
                                partz = z_min_local+l*dz+dz/curr%nppcell*(ipart-1)
                                partw = nc*dx*dy*dz/(curr%nppcell)
                                ! Sets velocity
                                CALL RANDOM_NUMBER(rng(1:3))
                                v=MAX(1e-10_num,rng(1))
                                th=2*pi*rng(2)
                                phi=2*pi*rng(3)
                                partux= curr%vdrift_x + curr%vth_x*sqrt(-2.*LOG(v))*COS(th)*COS(phi)
                                partuy= curr%vdrift_y + curr%vth_y*sqrt(-2.*LOG(v))*COS(th)*SIN(phi)
                                partuz= curr%vdrift_z + curr%vth_z*sqrt(-2.*LOG(v))*SIN(th)
                                usq = (partux**2 + partuy**2+partuz**2)*clightsq
       							gaminv = 1.0_num/sqrt(1.0_num + usq)
                                ! Adds particle to array of tiles of current species
                                CALL add_particle_to_species(curr, partx, party, partz, &
                                partux, partuy, partuz, gaminv, partw)
                            END DO
                        END DO
                    END DO
                END DO
            END DO ! END LOOP ON SPECIES
        ENDIF

        IF (pdistr .EQ. 2) THEN
            DO ispecies=1,nspecies
                curr=>species_parray(ispecies)
                DO j=0,nx-1
                    DO k=0,ny-1
                        DO l=0,nz-1
                            DO ipart=1,curr%nppcell
                                CALL RANDOM_NUMBER(rng(1:6))
                                ! Sets positions and weight
                                partx = x_min_local+MIN(rng(1),0.999)*(x_max_local-x_min_local)
                                party = y_min_local+MIN(rng(2),0.999)*(y_max_local-y_min_local)
                                partz = z_min_local+MIN(rng(3),0.999)*(z_max_local-z_min_local)
                                partw = nc*dx*dy*dz/(curr%nppcell)
                                ! Sets velocity
                                v=MAX(1e-10_num,rng(4))
                                th=2*pi*rng(5)
                                phi=2*pi*rng(6)
                                partux= curr%vdrift_x + curr%vth_x*sqrt(-2.*LOG(v))*COS(th)*COS(phi)
                                partuy= curr%vdrift_y + curr%vth_y*sqrt(-2.*LOG(v))*COS(th)*SIN(phi)
                                partuz= curr%vdrift_z + curr%vth_z*sqrt(-2.*LOG(v))*SIN(th)
                                usq = (partux**2 + partuy**2+partuz**2)*clightsq
       							gaminv = 1.0_num/sqrt(1.0_num + usq)
                                ! Adds particle to array of tiles of current species
                                CALL add_particle_to_species(curr, partx, party, partz, &
                                partux, partuy, partuz, gaminv, partw)
                            END DO
                        END DO
                    END DO
                END DO
            END DO ! END LOOP ON SPECIES
        ENDIF

        ! Collects total number of particles from other subdomains (useful for statistics)
        ntot=0
        DO ispecies=1,nspecies
            curr=>species_parray(ispecies)
            CALL MPI_ALLREDUCE(curr%species_npart,npart,1_isp, MPI_INTEGER,MPI_SUM,comm, err)
            ntot=ntot+npart
            IF (rank .EQ. 0) THEN
                WRITE (0,*) 'Loaded npart = ', npart,' particles of species ', &
                TRIM(ADJUSTL(curr%name))
            END IF
        END DO

        RETURN
    END SUBROUTINE load_particles

    SUBROUTINE resize_particle_arrays(curr, old_size, new_size)
        IMPLICIT NONE

        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        INTEGER(idp) :: old_size, new_size

        curr%npmax_tile=new_size
        CALL resize_1D_array_real(curr%part_x, old_size, new_size)
        CALL resize_1D_array_real(curr%part_y, old_size, new_size)
        CALL resize_1D_array_real(curr%part_z, old_size, new_size)
        CALL resize_1D_array_real(curr%part_ux, old_size, new_size)
        CALL resize_1D_array_real(curr%part_uy, old_size, new_size)
        CALL resize_1D_array_real(curr%part_uz, old_size, new_size)
        CALL resize_1D_array_real(curr%part_gaminv, old_size, new_size)
        CALL resize_2D_array_real(curr%pid, old_size, new_size,npid,npid)
        CALL resize_1D_array_real(curr%part_ex, old_size, new_size)
        CALL resize_1D_array_real(curr%part_ey, old_size, new_size)
        CALL resize_1D_array_real(curr%part_ez, old_size, new_size)
        CALL resize_1D_array_real(curr%part_bx, old_size, new_size)
        CALL resize_1D_array_real(curr%part_by, old_size, new_size)
        CALL resize_1D_array_real(curr%part_bz, old_size, new_size)
    END SUBROUTINE resize_particle_arrays

    SUBROUTINE resize_1D_array_real(arr, old_size, new_size)
        IMPLICIT NONE
        REAL(num), DIMENSION(:),ALLOCATABLE, INTENT(IN OUT) :: arr
        REAL(num), DIMENSION(:),ALLOCATABLE :: temp
        INTEGER(idp) :: old_size, new_size

        ALLOCATE(temp(1:new_size))
        ! reshape array
        temp(1:old_size)=arr(1:old_size)
        DEALLOCATE(arr)
        ALLOCATE(arr(1:new_size))
        arr(1:old_size) = temp(1:old_size)
        DEALLOCATE(temp)
    END SUBROUTINE resize_1D_array_real
    
    SUBROUTINE resize_2D_array_real(arr, nx_old,nx_new,ny_old,ny_new)
        IMPLICIT NONE
        REAL(num), DIMENSION(:,:),ALLOCATABLE, INTENT(IN OUT) :: arr
        INTEGER(idp),INTENT(IN) :: nx_old,ny_old,nx_new,ny_new
        REAL(num), DIMENSION(:,:),ALLOCATABLE :: temp

        ALLOCATE(temp(1:nx_new,1:ny_new))
        ! reshape array
        temp(1:nx_old,1:ny_old)=arr(1:nx_old,1:ny_old)
        DEALLOCATE(arr)
        ALLOCATE(arr(1:nx_new,1:ny_new))
        arr(1:nx_old,1:ny_old) = temp(1:nx_old,1:ny_old)
        DEALLOCATE(temp)
    END SUBROUTINE resize_2D_array_real

    SUBROUTINE get_local_number_of_part(npart)
    	INTEGER(idp), INTENT(IN OUT) :: npart
    	INTEGER(idp) :: ispecies
    	TYPE(particle_species), POINTER :: curr
    	
    	npart=0
		DO ispecies=1, nspecies ! LOOP ON SPECIES
			curr=> species_parray(ispecies)
			npart=npart+curr%species_npart
		END DO ! END LOOP ON SPECIES

    END SUBROUTINE get_local_number_of_part
    ! ----- SUBROUTINES DEDICATED FOR PYTHON INTERFACE
    
    !This subroutine returns pointer arrays on a given tile 
    ! of a given species (USED mainly by python interface)
    SUBROUTINE point_to_tile(ispecies, ix, iy, iz)
        USE python_pointers
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN) :: ix,iy,iz,ispecies
        TYPE(particle_species), POINTER  :: currsp
        TYPE(particle_tile), POINTER ::curr_tile

        currsp=> species_parray(ispecies)
        curr_tile=>currsp%array_of_tiles(ix,iy,iz)
		! Tile extent and dimension
		nxtg=curr_tile%nxg_tile
		nytg=curr_tile%nyg_tile
		nztg=curr_tile%nzg_tile
		nxgt=curr_tile%nx_grid_tile
		nygt=curr_tile%ny_grid_tile
		nzgt=curr_tile%nz_grid_tile
		nxct=curr_tile%nx_cells_tile
		nyct=curr_tile%ny_cells_tile
		nzct=curr_tile%nz_cells_tile
		nxmin=curr_tile%nx_tile_min
		nxmax=curr_tile%nx_tile_max
		nymin=curr_tile%ny_tile_min
		nymax=curr_tile%ny_tile_max
		nzmin=curr_tile%nz_tile_min
		nzmax=curr_tile%nz_tile_max
		xtmin=curr_tile%x_tile_min
		xtmax=curr_tile%x_tile_max
		ytmin=curr_tile%y_tile_min
		ytmax=curr_tile%y_tile_max
		ztmin=curr_tile%z_tile_min
		ztmax=curr_tile%z_tile_max
		xgtmin=curr_tile%x_grid_tile_min
		xgtmax=curr_tile%x_grid_tile_max
		ygtmin=curr_tile%y_grid_tile_min
		ygtmax=curr_tile%y_grid_tile_max
		zgtmin=curr_tile%z_grid_tile_min
		zgtmax=curr_tile%z_grid_tile_max
		! Number of particles in the tile
		partn => curr_tile%np_tile
		partnmax = curr_tile%npmax_tile
		! Particle arrays in the tile
        partx=>curr_tile%part_x
        party=>curr_tile%part_y
        partz=>curr_tile%part_z
        partux=>curr_tile%part_ux
        partuy=>curr_tile%part_uy
        partuz=>curr_tile%part_uz
        partgaminv=>curr_tile%part_gaminv
        pid=>curr_tile%pid
        partex=>curr_tile%part_ex
        partey=>curr_tile%part_ey
        partez=>curr_tile%part_ez
        partbx=>curr_tile%part_bx
        partby=>curr_tile%part_by
        partbz=>curr_tile%part_bz


    END SUBROUTINE point_to_tile

    !This subroutine returns pointer arrays on a given tile
    ! of a given species (USED mainly by python interface)
    SUBROUTINE set_particle_species_properties(nsp,sname,mss,chrg,nppc,xsmin,ysmin,zsmin,xsmax,ysmax,zsmax, &
		vdxs,vdys,vdzs,vthxs,vthys,vthzs)
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN) :: nsp, nppc
		REAL(num), INTENT(IN) :: mss, chrg,xsmin,ysmin,zsmin,xsmax,ysmax,zsmax,vdxs,vdys,vdzs,vthxs,vthys,vthzs
		CHARACTER(LEN=*), INTENT(IN) :: sname
        TYPE(particle_species), POINTER  :: currsp

        currsp=> species_parray(nsp)
		currsp%charge=chrg
		currsp%mass=mss
		currsp%x_min=xsmin
		currsp%y_min=ysmin
		currsp%z_min=zsmin
		currsp%x_max=xsmax
		currsp%y_max=ysmax
		currsp%z_max=zsmax
		currsp%vdrift_x=vdxs
		currsp%vdrift_y=vdys
		currsp%vdrift_z=vdzs
		currsp%vth_x=vthxs
		currsp%vth_y=vthys
		currsp%vth_z=vthzs
		currsp%nppcell=nppc
		currsp%name=sname
		!PRINT *, "species name", sname
		!PRINT *, "species mass", mss
		!PRINT *, "species charge", chrg

    END SUBROUTINE set_particle_species_properties

    !!! --- Add particle to array of tiles
    SUBROUTINE py_add_particles_to_species(nsp, npart, partx, party, partz, &
               partux, partuy, partuz, gaminv, partw)
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN) :: nsp, npart
        REAL(num), DIMENSION(npart), INTENT(IN) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
        TYPE(particle_species), POINTER :: currsp
        INTEGER(idp) :: i
        currsp=>species_parray(nsp)

        DO i=1,npart
            CALL add_particle_to_species(currsp, partx(i), party(i), partz(i), &
                partux(i), partuy(i), partuz(i), gaminv(i), partw(i))
        END DO
    END SUBROUTINE py_add_particles_to_species
        
    !!! --- Get logical array are_tiles_reallocated for a given species
    SUBROUTINE get_are_tiles_reallocated(nsp, ntx, nty, ntz, atrealloc)
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN) :: nsp, ntx,nty,ntz
        INTEGER(idp), DIMENSION(ntx,nty,ntz), INTENT(IN OUT) :: atrealloc
        TYPE(particle_species), POINTER :: currsp
        currsp=>species_parray(nsp)

        atrealloc=currsp%are_tiles_reallocated
        
    END SUBROUTINE get_are_tiles_reallocated
    
     !!! --- Set logical array are_tiles_reallocated for a given species
    SUBROUTINE set_are_tiles_reallocated(nsp, ntx, nty, ntz, atrealloc)
        IMPLICIT NONE
        INTEGER(idp), INTENT(IN) :: nsp, ntx,nty,ntz
        INTEGER(idp), DIMENSION(ntx,nty,ntz), INTENT(IN) :: atrealloc
        TYPE(particle_species), POINTER :: currsp
        currsp=>species_parray(nsp)

        currsp%are_tiles_reallocated=atrealloc
        
    END SUBROUTINE set_are_tiles_reallocated


END MODULE tiling
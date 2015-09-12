MODULE tiling
!!!! --- This module contains useful diagnostics to test code correctness
    USE constants
    USE particles
    USE shared_data
    USE fields
    USE params
    IMPLICIT NONE

CONTAINS

    !!! --- Set particle tile split in space
    SUBROUTINE set_tile_split
        IMPLICIT NONE
        INTEGER :: ix, iy, iz, ispecies
        INTEGER :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
        INTEGER :: nx0_last_tile, ny0_last_tile, nz0_last_tile
        TYPE(particle_species), POINTER :: curr_sp
        TYPE(particle_tile), POINTER :: curr
        ! Tile-split
        nx0_grid_tile = nx_grid / ntilex
        ny0_grid_tile = ny_grid / ntiley
        nz0_grid_tile = nz_grid / ntilez
        PRINT *, "nx_grid", nx_grid, "ntilex", ntilex,"nx_grid_tile", nx0_grid_tile
        PRINT *, "ny_grid", ny_grid, "ntiley", ntiley,"ny_grid_tile", ny0_grid_tile
        PRINT *, "nz_grid", nz_grid, "ntilez", ntilez,"nz_grid_tile", nz0_grid_tile
        !-- N.B: If the number of grid points cannot be equally divided between
        !-- tiles then give remaining points to last tile in each dimension
        nx0_last_tile= nx0_grid_tile+(nx_grid-nx0_grid_tile*ntilex)
        ny0_last_tile= ny0_grid_tile+(ny_grid-ny0_grid_tile*ntiley)
        nz0_last_tile= nz0_grid_tile+(nz_grid-nz0_grid_tile*ntilez)

        !- Allocate object array of tiles
        DO ispecies =1, nspecies
            curr_sp => species_parray(ispecies)
            IF (.NOT. curr_sp%l_arrayoftiles_allocated) THEN
                ALLOCATE(curr_sp%array_of_tiles(ntilex,ntiley,ntilez))
                curr_sp%l_arrayoftiles_allocated = .TRUE.
            END IF
            ! Sets tile spatial extents for current species
            DO iz=1, ntilez
                DO iy=1,ntiley
                    DO ix=1,ntilex
                        curr=> curr_sp%array_of_tiles(ix,iy,iz)
                        ! X- Direction
                        IF (ix .LT. ntilex) THEN
                            curr%nx_grid_tile=nx0_grid_tile
                            curr%nx_cells_tile=curr%nx_grid_tile-1
                            curr%x_tile_min=x_min_local+(ix-1)*nx0_grid_tile*dx
                            curr%x_tile_max=curr%x_tile_min+nx0_grid_tile*dx
                            curr%nx_tile_min = (ix-1)*nx0_grid_tile
                            curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
                        ELSE ! LAST TILE in X DIRECTION
                            curr%nx_grid_tile=nx0_last_tile
                            curr%nx_cells_tile=curr%nx_grid_tile-1
                            curr%x_tile_min=x_min_local+(ix-1)*nx0_grid_tile*dx
                            curr%x_tile_max=curr%x_tile_min+nx0_last_tile*dx
                            curr%nx_tile_min = (ix-1)*nx0_grid_tile
                            curr%nx_tile_max = curr%nx_tile_min+curr%nx_cells_tile
                        ENDIF
                        ! Y -DIRECTION
                        IF (iy .LT. ntiley) THEN
                            curr%ny_grid_tile=ny0_grid_tile
                            curr%ny_cells_tile=curr%ny_grid_tile-1
                            curr%y_tile_min=y_min_local+(iy-1)*ny0_grid_tile*dy
                            curr%y_tile_max=curr%y_tile_min+ny0_grid_tile*dy
                            curr%ny_tile_min = (iy-1)*ny0_grid_tile
                            curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
                        ELSE ! LAST TILE in Y DIRECTION
                            curr%ny_grid_tile=ny0_last_tile
                            curr%ny_cells_tile=curr%ny_grid_tile-1
                            curr%y_tile_min=y_min_local+(iy-1)*ny0_grid_tile*dy
                            curr%y_tile_max=curr%y_tile_min+ny0_last_tile*dy
                            curr%ny_tile_min = (iy-1)*ny0_grid_tile
                            curr%ny_tile_max = curr%ny_tile_min+curr%ny_cells_tile
                        ENDIF
                        ! Z- DIRECTION
                        IF (iz .LT. ntilez) THEN
                            curr%nz_grid_tile=nz0_grid_tile
                            curr%nz_cells_tile=curr%nz_grid_tile-1
                            curr%z_tile_min=z_min_local+(iz-1)*nz0_grid_tile*dz
                            curr%z_tile_max=curr%z_tile_min+nz0_grid_tile*dz
                            curr%nz_tile_min = (iz-1)*nz0_grid_tile
                            curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
                        ELSE ! LAST TILE in Z DIRECTION
                            curr%nz_grid_tile=nz0_last_tile
                            curr%nz_cells_tile=curr%nz_grid_tile-1
                            curr%z_tile_min=z_min_local+(iz-1)*nz0_grid_tile*dz
                            curr%z_tile_max=curr%z_tile_min+nz0_last_tile*dz
                            curr%nz_tile_min = (iz-1)*nz0_grid_tile
                            curr%nz_tile_max = curr%nz_tile_min+curr%nz_cells_tile
                        ENDIF
                    END DO
                END DO
            END DO
        END DO ! END DO SPECIES
    END SUBROUTINE

    !!! --- Add particle to array of tiles
    SUBROUTINE add_particle_to_species(currsp, partx, party, partz, &
               partux, partuy, partuz, partw)
        IMPLICIT NONE
        REAL(num) :: partx, party, partz, partux, partuy, partuz, partw
        TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
        TYPE(particle_tile), POINTER :: curr
        INTEGER :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile, nptile
        INTEGER :: ixtile, iytile, iztile
        REAL(num) :: resize_factor = 1.5_num


        ! Get first tiles dimensions (may be different from last tile)
        nx0_grid_tile = currsp%array_of_tiles(1,1,1)%nx_grid_tile
        ny0_grid_tile = currsp%array_of_tiles(1,1,1)%ny_grid_tile
        nz0_grid_tile = currsp%array_of_tiles(1,1,1)%nz_grid_tile

        ! Get particle index in array of tile
        ixtile = FLOOR((partx-x_min_local)/(nx0_grid_tile*dx))+1
        iytile = FLOOR((party-y_min_local)/(ny0_grid_tile*dy))+1
        iztile = FLOOR((partz-z_min_local)/(nz0_grid_tile*dz))+1

        ! Point to current tile arr_of_tiles(ixtile,iytile,iztile)
        curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
        CALL add_particle_at_tile(curr, partx, party, partz, &
                partux, partuy, partuz, partw)

        ! Update total number of particle species
        currsp%species_npart=currsp%species_npart+1
    END SUBROUTINE add_particle_to_species

    SUBROUTINE add_particle_at_tile(curr, partx, party, partz, &
                partux, partuy, partuz, partw)
        IMPLICIT NONE
        INTEGER :: count, nmax
        REAL(num) :: partx, party, partz, partux, partuy, partuz, partw
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        ! If no particles in tile, allocate particle arrays
        IF (.NOT. curr%l_arrays_allocated) THEN
            CALL allocate_tile_arrays(curr)
        ENDIF

        ! Sanity check for max number of particles in tile
        count = curr%np_tile+1
        nmax  = curr%npmax_tile
        IF (count .GT. nmax) THEN
        ! Resize particle tile arrays if tile is full
            CALL resize_particle_arrays(curr, nmax, NINT(resize_factor*nmax+1))
        ENDIF
        ! Finally, add particle to tile
        curr%np_tile=count
        curr%part_x(count)  = partx
        curr%part_y(count)  = party
        curr%part_z(count)  = partz
        curr%part_ux(count) = partux
        curr%part_uy(count) = partuy
        curr%part_uz(count) = partuz
        curr%weight(count)  = partw
    END SUBROUTINE add_particle_at_tile

    !!! --- Remove particles from tile using a mask variable
    !!! --- This technique avoids packing or reallocating arrays
    SUBROUTINE rm_particles_from_species(currsp, curr, mask)
        TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        LOGICAL, DIMENSION (:), INTENT(IN) :: mask
        INTEGER :: ninit, i
        ninit= curr%np_tile
        DO i = ninit,1,-1
            IF (.NOT. mask(i)) THEN
                CALL rm_particle_at_tile(curr,i)
                currsp%species_npart=currsp%species_npart-1
            ENDIF
        ENDDO
    END SUBROUTINE rm_particles_from_species

    SUBROUTINE rm_particle_at_tile(curr, index)
        IMPLICIT NONE
        INTEGER :: index
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        IF (index .EQ. curr%np_tile) THEN
            ! If particle i is last element
            ! Simply decreases particle number
            curr%np_tile=curr%np_tile-1
        ELSE
            ! Particle i replaced by last element
            ! Particle number is decreased
            curr%part_x(index)=curr%part_x(curr%np_tile)
            curr%part_y(index)=curr%part_y(curr%np_tile)
            curr%part_z(index)=curr%part_z(curr%np_tile)
            curr%part_ux(index)=curr%part_ux(curr%np_tile)
            curr%part_uy(index)=curr%part_uy(curr%np_tile)
            curr%part_uz(index)=curr%part_uz(curr%np_tile)
            curr%weight(index)=curr%weight(curr%np_tile)
            curr%np_tile=curr%np_tile-1
        END IF
    END SUBROUTINE rm_particle_at_tile

    SUBROUTINE allocate_tile_arrays(curr_tile)
        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr_tile
        INTEGER :: nmax, nxc, nyc, nzc
        ! ALLOCATE PARTICLE ARRAYS
        nmax = curr_tile%npmax_tile
        ALLOCATE(curr_tile%part_x(1:nmax))
        ALLOCATE(curr_tile%part_y(1:nmax))
        ALLOCATE(curr_tile%part_z(1:nmax))
        ALLOCATE(curr_tile%part_ux(1:nmax))
        ALLOCATE(curr_tile%part_uy(1:nmax))
        ALLOCATE(curr_tile%part_uz(1:nmax))
        ALLOCATE(curr_tile%weight(1:nmax))
        ALLOCATE(curr_tile%part_ex(1:nmax))
        ALLOCATE(curr_tile%part_ey(1:nmax))
        ALLOCATE(curr_tile%part_ez(1:nmax))
        ALLOCATE(curr_tile%part_bx(1:nmax))
        ALLOCATE(curr_tile%part_by(1:nmax))
        ALLOCATE(curr_tile%part_bz(1:nmax))

        ! ALLOCATE GRID ARRAYS
        nxc=curr_tile%nx_cells_tile
        nyc=curr_tile%ny_cells_tile
        nzc=curr_tile%nz_cells_tile
        ALLOCATE(curr_tile%rho_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz), &
                 curr_tile%jx_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%jy_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%jz_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%ex_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%ey_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%ez_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%bx_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%by_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz),  &
                 curr_tile%bz_tile(-nox:nxc+nox,-noy:nyc+noy,-noz:nzc+noz))

        curr_tile%l_arrays_allocated = .TRUE.

    END SUBROUTINE allocate_tile_arrays

    SUBROUTINE init_tile_arrays
        IMPLICIT NONE
        INTEGER :: ispecies, ix, iy, iz
        INTEGER :: n1, n2, n3
        TYPE(particle_tile), POINTER :: curr_tile
        TYPE(particle_species), POINTER :: curr

        DO ispecies=1,nspecies ! LOOP ON SPECIES
            curr=>species_parray(ispecies)
            curr%species_npart=0
            DO iz=1, ntilez! LOOP ON TILES
                DO iy=1, ntiley
                    DO ix=1,ntilex
                        curr_tile=>curr%array_of_tiles(ix,iy,iz)
                        ! - Max size of particle arrays of current ile
                        n1=curr_tile%nx_cells_tile
                        n2=curr_tile%ny_cells_tile
                        n3=curr_tile%nz_cells_tile
                        curr_tile%npmax_tile=n1*n2*n3*curr%nppcell
                        curr_tile%np_tile=0
                        ! - Allocate arrays of current tile
                        CALL allocate_tile_arrays(curr_tile)
                        ! - Init array of current tile
                        curr_tile%part_x=0.0_num
                        curr_tile%part_y=0.0_num
                        curr_tile%part_z=0.0_num
                        curr_tile%part_ux=0.0_num
                        curr_tile%part_uy=0.0_num
                        curr_tile%part_uz=0.0_num
                        curr_tile%part_ex=0.0_num
                        curr_tile%part_ey=0.0_num
                        curr_tile%part_ez=0.0_num
                        curr_tile%part_bx=0.0_num
                        curr_tile%part_by=0.0_num
                        curr_tile%part_bz=0.0_num
                        curr_tile%weight=0.0_num
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE init_tile_arrays

    SUBROUTINE load_particles
        IMPLICIT NONE
        TYPE(particle_species), POINTER :: curr
        INTEGER :: ispecies, l, k, j, ipart
        INTEGER :: jmin, jmax, kmin, kmax, lmin, lmax
        REAL(num) :: partx, party, partz, partux, partuy, partuz, partw
        REAL(num) :: phi, th, v
        !!! --- Sets-up particle space distribution (default is homogeneous)
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
                            v=MAX(1e-10_num,RAND())
                            th=2*pi*RAND()
                            phi=2*pi*RAND()
                            partux= curr%vdrift_x !+ curr%vth_x*sqrt(-2.*LOG(v))*COS(th)*COS(phi)
                            partuy= curr%vdrift_y !+ curr%vth_y*sqrt(-2.*LOG(v))*COS(th)*SIN(phi)
                            partuz= curr%vdrift_z !+ curr%vth_z*sqrt(-2.*LOG(v))*SIN(th)
                            ! Adds particle to array of tiles of current species
                            CALL add_particle_to_species(curr, partx, party, partz, &
                            partux, partuy, partuz, partw)
                        END DO
                    END DO
                END DO
            END DO
        END DO ! END LOOP ON SPECIES
    END SUBROUTINE load_particles

    SUBROUTINE resize_particle_arrays(curr, old_size, new_size)
        IMPLICIT NONE

        TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr
        INTEGER old_size, new_size

        curr%npmax_tile=new_size
        CALL resize_array_real(curr%part_x, old_size, new_size)
        CALL resize_array_real(curr%part_y, old_size, new_size)
        CALL resize_array_real(curr%part_z, old_size, new_size)
        CALL resize_array_real(curr%part_ux, old_size, new_size)
        CALL resize_array_real(curr%part_uy, old_size, new_size)
        CALL resize_array_real(curr%part_uz, old_size, new_size)
        CALL resize_array_real(curr%weight, old_size, new_size)
        CALL resize_array_real(curr%part_ex, old_size, new_size)
        CALL resize_array_real(curr%part_ey, old_size, new_size)
        CALL resize_array_real(curr%part_ez, old_size, new_size)
        CALL resize_array_real(curr%part_bx, old_size, new_size)
        CALL resize_array_real(curr%part_by, old_size, new_size)
        CALL resize_array_real(curr%part_bz, old_size, new_size)
    END SUBROUTINE resize_particle_arrays

    SUBROUTINE resize_array_real(arr, old_size, new_size)
        IMPLICIT NONE
        REAL(num), DIMENSION(:), POINTER, INTENT(IN OUT) :: arr
        REAL(num), DIMENSION(:), POINTER :: temp
        INTEGER old_size, new_size

        ALLOCATE(temp(1:new_size))
        ! reshape array
        temp(1:old_size)=arr(1:old_size)
        DEALLOCATE(arr)
        ALLOCATE(arr(1:new_size))
        arr(1:old_size) = temp(1:old_size)
        DEALLOCATE(temp)
    END SUBROUTINE

END MODULE tiling
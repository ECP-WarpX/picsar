! ________________________________________________________________________________________
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
! TILING.F90
!
! Purpose:
! This file contains useful subroutines for the tiling.
!
! Authors:
! Henri Vincenti
!
! Date
! Creation 2015
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Module that contains subroutines for the particle/grid tiling.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation: 2015
! ________________________________________________________________________________________
MODULE tiling
  USE constants
  USE particles
  USE shared_data
  USE fields
  USE params
  IMPLICIT NONE

  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> Main subroutine that sets the tile split for all species.
  !> This subroutine is called by ::initall
  !>
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE set_tile_split()
    IMPLICIT NONE

    ! Set tile split for species arrays
    CALL set_tile_split_for_species(species_parray, nspecies, ntilex, ntiley, ntilez, &
    nx_grid, ny_grid, nz_grid, x_min_local, y_min_local, z_min_local, x_max_local,    &
    y_max_local, z_max_local)

    ! ALLOCATE grid tile arrays
    ALLOCATE(aofgrid_tiles(ntilex, ntiley, ntilez))

  END SUBROUTINE set_tile_split


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine sets particle tile split in space for a given species
  !> This subroutine is called by ::set_tile_split
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE set_tile_split_for_species(species_array, nspec, ntx, nty, ntz, nxgrid,  &
    nygrid, nzgrid, xminlocal, yminlocal, zminlocal, xmaxlocal, ymaxlocal, zmaxlocal)

    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)        :: nspec, nxgrid, nygrid, nzgrid
    INTEGER(idp), INTENT(IN OUT)    ::  ntx, nty, ntz
    REAL(num), INTENT(IN)           :: xminlocal, yminlocal, zminlocal, xmaxlocal,    &
    ymaxlocal, zmaxlocal
    TYPE(particle_species), INTENT(IN OUT), TARGET, DIMENSION(nspec) :: species_array
    INTEGER(idp)                    :: ix, iy, iz, ispecies
    INTEGER(idp)                    :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    INTEGER(idp)                    :: nx0_last_tile, ny0_last_tile, nz0_last_tile
    TYPE(particle_species), POINTER :: curr_sp
    TYPE(particle_tile), POINTER    :: curr

    ! Tile-split
    nx0_grid_tile = nxgrid / ntx
    ny0_grid_tile = nygrid / nty
    nz0_grid_tile = nzgrid / ntz

    ! Some sanity check
    IF (nx0_grid_tile .LT. 4) THEN
      IF (rank .EQ. 0) PRINT *, "number of tiles in X is too high, settting back to   &
      default value 1"
      ntx=1
    END IF
    IF (ny0_grid_tile .LT. 4) THEN
      IF(c_dim .EQ. 3) THEN
        IF (rank .EQ. 0) PRINT *, "number of tiles in Y is too high, setting back to  &
        default value 1"
      ENDIF
      nty=1
    END IF
    IF (nz0_grid_tile .LT. 4) THEN
      IF (rank .EQ. 0) PRINT *, "number of tiles in Z is too high, setting back to    &
      default value 1"
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
        ALLOCATE(curr_sp%array_of_tiles(ntx, nty, ntz))
        ALLOCATE(curr_sp%are_tiles_reallocated(ntx, nty, ntz))
        curr_sp%are_tiles_reallocated= 0
        curr_sp%l_arrayoftiles_allocated = .TRUE.
      END IF
      ! Sets tile spatial extents for current species
      DO iz=1, ntz
        DO iy=1, nty
          DO ix=1, ntx
            curr=> curr_sp%array_of_tiles(ix, iy, iz)
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
              IF (c_dim .EQ. 3) curr%subdomain_bound = .TRUE.
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
              IF (c_dim .EQ. 3) curr%subdomain_bound = .TRUE.
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
    END DO! END DO SPECIES
  END SUBROUTINE set_tile_split_for_species

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine adds particle of given species to the corresponding
  !> tile particle array depending on the particle position in 3D.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE add_particle_to_species(currsp, partx, party, partz, partux, partuy,     &
    partuz, gaminv, partpid)
    IMPLICIT NONE
    REAL(num), INTENT(IN) :: partx, party, partz, partux, partuy, partuz
    REAL(num), INTENT(IN) :: gaminv
    REAL(num), DIMENSION(:), INTENT(IN) :: partpid
    TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
    INTEGER(idp) :: nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    INTEGER(idp) :: ixtile, iytile, iztile

    ! Get first tiles dimensions (may be different from last tile)
    nx0_grid_tile = currsp%array_of_tiles(1, 1, 1)%nx_grid_tile
    ny0_grid_tile = currsp%array_of_tiles(1, 1, 1)%ny_grid_tile
    nz0_grid_tile = currsp%array_of_tiles(1, 1, 1)%nz_grid_tile

    ! Get particle index in array of tile
    ixtile = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx), idp)+1,       &
    ntilex)
    iytile = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy), idp)+1,       &
    ntiley)
    iztile = MIN(FLOOR((partz-(z_min_local)+dz/2_num)/(nz0_grid_tile*dz), idp)+1,     &
    ntilez)

    CALL add_particle_at_tile(currsp, ixtile, iytile, iztile, partx, party, partz,    &
    partux, partuy, partuz, gaminv, partpid)

    ! Update total number of particle species
    currsp%species_npart=currsp%species_npart+1
  END SUBROUTINE add_particle_to_species


  ! ______________________________________________________________________________________
  !> @brief
  !> In 3D, add a particle with its properties to the list of particles
  !> inside the tile. This subroutine is called in add_particle_to_species().
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE add_particle_at_tile(currsp, ixt, iyt, izt, partx, party, partz, partux, &
    partuy, partuz, gaminv, partpid)
    IMPLICIT NONE
    INTEGER(idp) :: count, nmax, ixt, iyt, izt
    REAL(num), INTENT(IN) :: partx, party, partz, partux, partuy, partuz
    REAL(num), INTENT(IN) :: gaminv
    REAL(num), DIMENSION(:) :: partpid
    TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
    TYPE(particle_tile), POINTER :: curr

    curr=>currsp%array_of_tiles(ixt, iyt, izt)
    ! If no particles in tile, allocate particle arrays
    IF (.NOT. curr%l_arrays_allocated) THEN
      CALL allocate_tile_arrays(curr)
    ENDIF

    ! Sanity check for max number of particles in tile
    count = curr%np_tile(1)+1
    nmax  = curr%npmax_tile
    IF (count .GT. nmax) THEN
      ! Resize particle tile arrays if tile is full
      currsp%are_tiles_reallocated(ixt, iyt, izt)=1
      CALL resize_particle_arrays(curr, nmax, NINT(resize_factor*nmax+1, idp))
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
    curr%pid(count, 1:npid) = partpid
    curr%part_ex(count)  = 0._num
    curr%part_ey(count)  = 0._num
    curr%part_ez(count)  = 0._num
    curr%part_bx(count)  = 0._num
    curr%part_by(count)  = 0._num
    curr%part_bz(count)  = 0._num
  END SUBROUTINE add_particle_at_tile

  ! ______________________________________________________________________________________
  !> @brief
  !> In 3D, add a group of particles with their properties to the list of particles
  !> inside the tile.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE add_group_of_particles_at_tile(currsp, ixt, iyt, izt, np, npiid, partx,  &
    party, partz, partux, partuy, partuz, gaminv, partpid)
    IMPLICIT NONE
    INTEGER(idp) :: count, nmax, ixt, iyt, izt, np, npnew, npiid
    REAL(num), DIMENSION(np) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(np, npiid) :: partpid
    TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
    TYPE(particle_tile), POINTER :: curr

    curr=>currsp%array_of_tiles(ixt, iyt, izt)
    ! If no particles in tile, allocate particle arrays
    IF (.NOT. curr%l_arrays_allocated) THEN
      CALL allocate_tile_arrays(curr)
    ENDIF

    ! Sanity check for max number of particles in tile
    count = curr%np_tile(1)
    npnew = count + np
    nmax  = curr%npmax_tile
    IF (npnew .GT. nmax) THEN
      ! Resize particle tile arrays if tile is full
      currsp%are_tiles_reallocated(ixt, iyt, izt)=1
      CALL resize_particle_arrays(curr, nmax, NINT(resize_factor*nmax+1, idp))
    ENDIF
    ! Finally, add particle to tile
    curr%np_tile(1)=npnew
    curr%part_x(count+1:npnew)  = partx
    curr%part_y(count+1:npnew)  = party
    curr%part_z(count+1:npnew)  = partz
    curr%part_ux(count+1:npnew) = partux
    curr%part_uy(count+1:npnew) = partuy
    curr%part_uz(count+1:npnew) = partuz
    curr%part_gaminv(count+1:npnew) = gaminv
    curr%pid(count+1:npnew, 1:npiid) = partpid
    curr%part_ex(count+1:npnew)  = 0._num
    curr%part_ey(count+1:npnew)  = 0._num
    curr%part_ez(count+1:npnew)  = 0._num
    curr%part_bx(count+1:npnew)  = 0._num
    curr%part_by(count+1:npnew)  = 0._num
    curr%part_bz(count+1:npnew)  = 0._num
  END SUBROUTINE add_group_of_particles_at_tile


  ! ______________________________________________________________________________________
  !> @brief
  !> Remove particles from tile using a mask variable.
  !> This technique avoids packing or reallocating arrays.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE rm_particles_from_species_with_mask(currsp, ixt, iyt, izt, mask)
    TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
    LOGICAL(lp), DIMENSION (:), INTENT(IN) :: mask
    INTEGER(idp), INTENT(IN) :: ixt, iyt, izt
    INTEGER(idp) :: ninit, i
    TYPE(particle_tile), POINTER :: curr

    curr=>currsp%array_of_tiles(ixt, iyt, izt)
    ninit= curr%np_tile(1)
    DO i = ninit, 1, -1
      IF (.NOT. mask(i)) THEN
        CALL rm_particle_at_tile(currsp, ixt, iyt, izt, i)
        currsp%species_npart=currsp%species_npart-1
      ENDIF
    ENDDO
  END SUBROUTINE rm_particles_from_species_with_mask

  ! ______________________________________________________________________________________
  !> @brief
  !> Remove a particle in a given tile from species currsp.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE rm_particles_from_species(currsp, ixt, iyt, izt, ipart)
    TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
    INTEGER(idp), INTENT(IN) :: ipart, ixt, iyt, izt

    CALL rm_particle_at_tile(currsp, ixt, iyt, izt, ipart)
    currsp%species_npart=currsp%species_npart-1
  END SUBROUTINE rm_particles_from_species

  ! ______________________________________________________________________________________
  !> @brief
  !> Remove a particle in a given tile from species currsp.
  !> This subroutine is called in rm_particles_from_species.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE rm_particle_at_tile(currsp, ixt, iyt, izt, index)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)                        :: index, ixt, iyt, izt
    TYPE(particle_species), POINTER, INTENT(IN OUT) :: currsp
    TYPE(particle_tile), POINTER                    :: curr
    curr=>currsp%array_of_tiles(ixt, iyt, izt)

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
      curr%pid(index, 1:npid)=curr%pid(curr%np_tile(1), 1:npid)
      curr%np_tile=curr%np_tile(1)-1
    END IF

    ! Avoid memory leaks
    ! Reduce array size if # of particles in array lower than
    ! 30% of array size
    IF(curr%np_tile(1) .LT. FLOOR(downsize_threshold*curr%npmax_tile)) THEN
      IF (FLOOR(downsize_factor*curr%npmax_tile) .GT. 0) THEN
        CALL resize_particle_arrays(curr, curr%npmax_tile,                            &
        FLOOR(downsize_factor*curr%npmax_tile, idp))
        currsp%are_tiles_reallocated(ixt, iyt, izt)=1
      ENDIF
    ENDIF
  END SUBROUTINE rm_particle_at_tile


  ! ______________________________________________________________________________________
  !> @brief
  !> Allocate arrays in curr_tile for particles.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE allocate_tile_arrays(curr_tile)
    TYPE(particle_tile), POINTER, INTENT(IN OUT) :: curr_tile
    INTEGER(idp) :: nmax
    ! ALLOCATE PARTICLE ARRAYS
    nmax = curr_tile%npmax_tile
    ALLOCATE(curr_tile%part_x(1:nmax), curr_tile%part_y(1:nmax),                      &
    curr_tile%part_z(1:nmax), curr_tile%part_ux(1:nmax), curr_tile%part_uy(1:nmax),   &
    curr_tile%part_uz(1:nmax), curr_tile%pid(1:nmax, 1:npid),                         &
    curr_tile%part_ex(1:nmax), curr_tile%part_ey(1:nmax), curr_tile%part_ez(1:nmax),  &
    curr_tile%part_bx(1:nmax), curr_tile%part_by(1:nmax), curr_tile%part_bz(1:nmax),  &
    curr_tile%part_gaminv(1:nmax))
    curr_tile%l_arrays_allocated = .TRUE.

  END SUBROUTINE allocate_tile_arrays

  ! ______________________________________________________________________________________
  !> @brief
  !> Main subroutine to init arrays of tiles and species for tiling.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE init_tile_arrays()
    IMPLICIT NONE

    CALL init_tile_arrays_for_species(nspecies, species_parray, aofgrid_tiles,        &
    ntilex, ntiley, ntilez)

  END SUBROUTINE init_tile_arrays


  ! ______________________________________________________________________________________
  !
  !> This subroutine allocates and initialize (first touch) the particle property arrays
  !> for each tile.
  !> @brief
  !
  !> @details
  !> This subroutine also allocates the local field arrays contained
  !> in aofgrid_tiles for each tile.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  ! ______________________________________________________________________________________
  SUBROUTINE init_tile_arrays_for_species(nspec2, species_array, aofgtiles, ntx2,     &
    nty2, ntz2)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN)        :: nspec2, ntx2, nty2, ntz2
    TYPE(grid_tile), DIMENSION(ntx2, nty2, ntz2), INTENT(IN OUT)        :: aofgtiles
    TYPE(particle_species), DIMENSION(nspec2), TARGET, INTENT(IN OUT) ::              &
    species_array

    INTEGER(idp)                    :: ispecies, ix, iy, iz
    INTEGER(idp)                    :: n1, n2, n3, ng1, ng2, ng3
    TYPE(particle_tile), POINTER    :: curr_tile
    TYPE(particle_species), POINTER :: curr

    IF (nspec2 .EQ. 0) RETURN

    ! Allocate particle tile arrays
    DO ispecies=1, nspec2! LOOP ON SPECIES
      curr=>species_array(ispecies)
      curr%species_npart=0
      DO iz=1, ntz2! LOOP ON TILES
        DO iy=1, nty2
          DO ix=1, ntx2
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            ! - Max size of particle arrays of current ile
            n1=curr_tile%nx_cells_tile
            n2=curr_tile%ny_cells_tile
            n3=curr_tile%nz_cells_tile
            curr_tile%npmax_tile=MAX(n1*n2*n3*curr%nppcell, 1_idp)
            curr_tile%np_tile(1)=0
            ! Set number of guard cells for each tile
            IF ((ix .GT. 1) .AND. (ix .LT. ntx2)) THEN
              curr_tile%nxg_tile=MAX(nox+1, 2_idp)
            ELSE
              curr_tile%nxg_tile=nxjguards
            END IF
            IF ((iy .GT. 1) .AND. (iy .LT. nty2) .AND. (c_dim .EQ. 3)) THEN
              curr_tile%nyg_tile=MAX(noy+1, 2_idp)
            ELSE
              curr_tile%nyg_tile=nyjguards
            END IF
            IF ((iz .GT. 1) .AND. (iz .LT. ntz2)) THEN
              curr_tile%nzg_tile=MAX(noz+1, 2_idp)
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
    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE)                     &
    !$OMP SHARED(species_array, ntx2, nty2, ntz2, nspec2) PRIVATE(ix, iy, iz,         &
    !$OMP ispecies, curr, curr_tile)
    DO iz=1, ntz2! LOOP ON TILES
      DO iy=1, nty2
        DO ix=1, ntx2
          DO ispecies=1, nspec2! LOOP ON SPECIES
            curr=>species_array(ispecies)
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
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
    DO iz=1, ntz2! LOOP ON TILES
      DO iy=1, nty2
        DO ix=1, ntx2
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          ! - Max size of particle arrays of current ile
          n1=curr_tile%nx_cells_tile
          n2=curr_tile%ny_cells_tile
          n3=curr_tile%nz_cells_tile
          ng1=curr_tile%nxg_tile
          ng2=curr_tile%nyg_tile
          ng3=curr_tile%nzg_tile
          ALLOCATE(aofgtiles(ix, iy, iz)%arr1(-ng1:n1+ng1, -ng2:n2+ng2,               &
          -ng3:n3+ng3))
          ALLOCATE(aofgtiles(ix, iy, iz)%arr2(-ng1:n1+ng1, -ng2:n2+ng2,               &
          -ng3:n3+ng3))
          ALLOCATE(aofgtiles(ix, iy, iz)%arr3(-ng1:n1+ng1, -ng2:n2+ng2,               &
          -ng3:n3+ng3))
        END DO
      END DO
    END DO! END LOOP ON TILES
  END SUBROUTINE init_tile_arrays_for_species

  ! ______________________________________________________________________________________
  !> @brief
  !> Initialize the particle properties (positions and velocities) according to
  !> the specified distribution.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE load_particles
    IMPLICIT NONE
    TYPE(particle_species), POINTER :: curr
    INTEGER(idp) :: ispecies, l, k, j, ipart
    INTEGER(idp) :: jmin, jmax, kmin, kmax, lmin, lmax
    REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num) :: phi, th, v, clightsq, partvx, partvy, partvz
    REAL(num), DIMENSION(:), ALLOCATABLE :: partpid
    INTEGER(idp) :: npart
    INTEGER(isp) :: err
    REAL(num), DIMENSION(6) :: rng=0_num
    clightsq=1/clight**2
    ALLOCATE(partpid(npid))
    !!! --- Sets-up particle space distribution
    !!! --- (homogeneous case, uniform space distribution - default)
    IF (pdistr .EQ. 1) THEN
      DO ispecies=1, nspecies
        curr=>species_parray(ispecies)
        IF (curr%is_antenna) CYCLE
        IF (c_dim.eq.3) THEN
          jmin = NINT(MAX(curr%x_min-x_min_local, 0.0_num)/dx)
          jmax = NINT(MIN(curr%x_max-x_min_local, x_max_local-x_min_local)/dx)
          kmin = NINT(MAX(curr%y_min-y_min_local, 0.0_num)/dy)
          kmax = NINT(MIN(curr%y_max-y_min_local, y_max_local-y_min_local)/dy)
          lmin = NINT(MAX(curr%z_min-z_min_local, 0.0_num)/dz)
          lmax = NINT(MIN(curr%z_max-z_min_local, z_max_local-z_min_local)/dz)

          DO l=lmin, lmax-1
            DO k=kmin, kmax-1
              DO j=jmin, jmax-1
                DO ipart=1, curr%nppcell
                  ! Sets positions and weight
                  partx = x_min_local+j*dx+dx/curr%nppcell*(ipart-0.5_num)
                  party = y_min_local+k*dy+dy/curr%nppcell*(ipart-0.5_num)
                  partz = z_min_local+l*dz+dz/curr%nppcell*(ipart-0.5_num)
                  partpid(wpid) = nc*dx*dy*dz/(curr%nppcell)
                  ! Sets velocity
                  CALL RANDOM_NUMBER(rng(1:3))
                  v=MAX(1e-10_num, rng(1))
                  th=2*pi*rng(2)
                  phi=2*pi*rng(3)

                  partvx= curr%vdrift_x +                                             &
                  curr%vth_x*sqrt(-2.*LOG(v))*COS(th)*COS(phi)
                  partvy= curr%vdrift_y +                                             &
                  curr%vth_y*sqrt(-2.*LOG(v))*COS(th)*SIN(phi)
                  partvz= curr%vdrift_z + curr%vth_z*sqrt(-2.*LOG(v))*SIN(th)

                  gaminv = sqrt(1.0_num - (partvx**2 + partvy**2 +                    &
                  partvz**2)*clightsq)
                  partux = partvx /gaminv
                  partuy = partvy /gaminv
                  partuz = partvz /gaminv

                  ! Adds particle to array of tiles of current species
                  CALL add_particle_to_species(curr, partx, party, partz, partux,     &
                  partuy, partuz, gaminv, partpid)
                END DO
              END DO
            END DO
          END DO

        ENDIF
      END DO! END LOOP ON SPECIES

      !!! --- Sets-up particle space distribution (random space distribution)
    ELSE IF (pdistr .EQ. 2) THEN
      DO ispecies=1, nspecies
        curr=>species_parray(ispecies)
        IF (curr%is_antenna) CYCLE
        jmin = NINT(MAX(curr%x_min-x_min_local, 0.0_num)/dx)
        jmax = NINT(MIN(curr%x_max-x_min_local, x_max_local-x_min_local)/dx)
        kmin = NINT(MAX(curr%y_min-y_min_local, 0.0_num)/dy)
        kmax = NINT(MIN(curr%y_max-y_min_local, y_max_local-y_min_local)/dy)
        lmin = NINT(MAX(curr%z_min-z_min_local, 0.0_num)/dz)
        lmax = NINT(MIN(curr%z_max-z_min_local, z_max_local-z_min_local)/dz)
        DO l=lmin, lmax-1
          DO k=kmin, kmax-1
            DO j=jmin, jmax-1
              DO ipart=1, curr%nppcell
                CALL RANDOM_NUMBER(rng(1:6))
                ! Sets positions and weight
                partx = x_min_local+MIN(rng(1), 0.999_num)*(x_max_local-x_min_local)
                party = y_min_local+MIN(rng(2), 0.999_num)*(y_max_local-y_min_local)
                partz = z_min_local+MIN(rng(3), 0.999_num)*(z_max_local-z_min_local)
                partpid(wpid) = nc*dx*dy*dz/(curr%nppcell)
                ! Sets velocity
                v=MAX(1e-10_num, rng(4))
                th=2*pi*rng(5)
                phi=2*pi*rng(6)

                partvx= curr%vdrift_x + curr%vth_x*sqrt(-2.*LOG(v))*COS(th)*COS(phi)
                partvy= curr%vdrift_y + curr%vth_y*sqrt(-2.*LOG(v))*COS(th)*SIN(phi)
                partvz= curr%vdrift_z + curr%vth_z*sqrt(-2.*LOG(v))*SIN(th)

                gaminv = sqrt(1.0_num - (partvx**2 + partvy**2 + partvz**2)*clightsq)
                partux = partvx /gaminv
                partuy = partvy /gaminv
                partuz = partvz /gaminv

                ! Adds particle to array of tiles of current species
                CALL add_particle_to_species(curr, partx, party, partz, partux,       &
                partuy, partuz, gaminv, partpid)
              END DO
            END DO
          END DO
        END DO
      END DO! END LOOP ON SPECIES


      !!! --- Sets-up particle space distribution (random space with a given velocity)
    ELSE IF (pdistr .EQ. 3) THEN
      DO ispecies=1, nspecies
        IF (curr%is_antenna) CYCLE
        curr=>species_parray(ispecies)
        jmin = NINT(MAX(curr%x_min-x_min_local, 0.0_num)/dx)
        jmax = NINT(MIN(curr%x_max-x_min_local, x_max_local-x_min_local)/dx)
        kmin = NINT(MAX(curr%y_min-y_min_local, 0.0_num)/dy)
        kmax = NINT(MIN(curr%y_max-y_min_local, y_max_local-y_min_local)/dy)
        lmin = NINT(MAX(curr%z_min-z_min_local, 0.0_num)/dz)
        lmax = NINT(MIN(curr%z_max-z_min_local, z_max_local-z_min_local)/dz)
        DO l=lmin, lmax-1
          DO k=kmin, kmax-1
            DO j=jmin, jmax-1
              DO ipart=1, curr%nppcell
                CALL RANDOM_NUMBER(rng(1:6))
                ! Sets positions and weight
                partx = x_min_local+MIN(rng(1), 0.999_num)*(x_max_local-x_min_local)
                party = y_min_local+MIN(rng(2), 0.999_num)*(y_max_local-y_min_local)
                partz = z_min_local+MIN(rng(3), 0.999_num)*(z_max_local-z_min_local)
                partpid(wpid) = nc*dx*dy*dz/(curr%nppcell)

                ! Sets velocity
                v=MAX(1e-10_num, rng(4))
                th=2*pi*rng(5)
                phi=2*pi*rng(6)

                partvx= curr%vdrift_x + curr%vth_x*COS(th)*COS(phi)
                partvy= curr%vdrift_y + curr%vth_x*COS(th)*SIN(phi)
                partvz= curr%vdrift_z + curr%vth_x*SIN(th)

                gaminv = sqrt(1.0_num - (partvx**2 + partvy**2 + partvz**2)*clightsq)
                partux = partvx /gaminv
                partuy = partvy /gaminv
                partuz = partvz /gaminv

                ! Adds particle to array of tiles of current species
                CALL add_particle_to_species(curr, partx, party, partz, partux,       &
                partuy, partuz, gaminv, partpid)
              END DO
            END DO
          END DO
        END DO
      END DO! END LOOP ON SPECIES

    ENDIF

    ! Collects total number of particles from other subdomains (useful for statistics)
    ntot=0
    DO ispecies=1, nspecies
      curr=>species_parray(ispecies)
      IF (curr%is_antenna) CYCLE
      CALL MPI_ALLREDUCE(curr%species_npart, npart, 1_isp, MPI_INTEGER8, MPI_SUM,     &
      comm, err)
      ntot=ntot+npart
      IF (rank .EQ. 0) THEN
        WRITE (0, *) 'Loaded npart = ', npart, ' particles of species ',              &
        TRIM(ADJUSTL(curr%name))
      END IF
    END DO

    DEALLOCATE(partpid)
    RETURN
  END SUBROUTINE load_particles

  ! ______________________________________________________________________________________
  !> @brief
  !> Resize particle arrays when they reach a threshold.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  !
  ! ______________________________________________________________________________________
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
    CALL resize_2D_array_real(curr%pid, old_size, new_size, npid, npid)
    CALL resize_1D_array_real(curr%part_ex, old_size, new_size)
    CALL resize_1D_array_real(curr%part_ey, old_size, new_size)
    CALL resize_1D_array_real(curr%part_ez, old_size, new_size)
    CALL resize_1D_array_real(curr%part_bx, old_size, new_size)
    CALL resize_1D_array_real(curr%part_by, old_size, new_size)
    CALL resize_1D_array_real(curr%part_bz, old_size, new_size)
  END SUBROUTINE resize_particle_arrays

  ! ______________________________________________________________________________________
  !> @brief
  !> Resize a 1D array of reals.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE resize_1D_array_real(arr, old_size, new_size)
    IMPLICIT NONE
    REAL(num), DIMENSION(:), ALLOCATABLE, INTENT(IN OUT) :: arr
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp
    INTEGER(idp) :: old_size, new_size

    ! - Allocate temporary array for copying arr before its de-allocation/re-allocation
    ALLOCATE(temp(1:new_size))
    ! reshape array
    IF (new_size .GT. old_size) THEN
	  temp(1:old_size)=arr(1:old_size)
    ELSE
	  temp(1:new_size)=arr(1:new_size)
    ENDIF
    DEALLOCATE(arr)
    ALLOCATE(arr(1:new_size))
    arr=temp
    DEALLOCATE(temp)
  END SUBROUTINE resize_1D_array_real

  ! ______________________________________________________________________________________
  !> @brief
  !> Resize a 2D array of reals.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2015
  !
  ! ______________________________________________________________________________________
  SUBROUTINE resize_2D_array_real(arr, nx_old, nx_new, ny_old, ny_new)
    IMPLICIT NONE
    REAL(num), DIMENSION(:, :), ALLOCATABLE, INTENT(IN OUT) :: arr
    INTEGER(idp), INTENT(IN) :: nx_old, ny_old, nx_new, ny_new
    INTEGER(idp)            :: nx_temp, ny_temp
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: temp

    ! - Allocate temporary array for copying arr before its de-allocation/re-allocation
    ALLOCATE(temp(1:nx_new, 1:ny_new))
    ! reshape array
    IF (nx_new .GT. nx_old) THEN
	  nx_temp=nx_old
    ELSE
	  nx_temp=nx_new
    ENDIF
    IF (ny_new .GT. ny_old) THEN
	  ny_temp=ny_old
    ELSE
	  ny_temp=ny_new
    ENDIF
    temp(1:nx_temp, 1:ny_temp)= arr(1:nx_temp, 1:ny_temp)
    DEALLOCATE(arr)
    ALLOCATE(arr(1:nx_new, 1:ny_new))
    arr=temp
    DEALLOCATE(temp)
  END SUBROUTINE resize_2D_array_real

  ! ______________________________________________________________________________________
  !> @brief
  !> Resize a 3D array of reals.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation: 2017
  !
  ! ______________________________________________________________________________________
  SUBROUTINE resize_3D_array_real(arr, nx_old, nx_new, ny_old, ny_new, nz_old,        &
    nz_new)
    IMPLICIT NONE
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE, INTENT(IN OUT) :: arr
    INTEGER(idp), INTENT(IN) :: nx_old, ny_old, nz_old, nx_new, ny_new, nz_new
    INTEGER(idp)            :: nx_temp, ny_temp, nz_temp
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp
 
    ! - Allocate temporary array for copying arr before its de-allocation/re-allocation
    ALLOCATE(temp(1:nx_new, 1:ny_new, 1:nz_new))
    ! reshape array
    IF (nx_new .GT. nx_old) THEN
      nx_temp=nx_old
    ELSE
      nx_temp=nx_new
    ENDIF
    IF (ny_new .GT. ny_old) THEN
      ny_temp=ny_old
    ELSE
      ny_temp=ny_new
    ENDIF
    IF (nz_new .GT. nz_old) THEN
      nz_temp=nz_old
    ELSE
      nz_temp=nz_new
    ENDIF
    temp(1:nx_temp, 1:ny_temp, 1:nz_temp)= arr(1:nx_temp, 1:ny_temp, 1:nz_temp)
    DEALLOCATE(arr)
    ALLOCATE(arr(1:nx_new, 1:ny_new, 1:nz_new))
    arr=temp
    DEALLOCATE(temp)
  END SUBROUTINE resize_3D_array_real

 ! _______________________________________________________________________________________
 !> @brief
 !> Subroutine that gets total memory size in Bytes occupied by tile structures (grid/
 !> particles) on local rank. These memory sizes are stored in the mem_status module 
 !> variables: local_grid_tiles_mem (for grid) and local_part_tiles_mem (for particles)
 !
 !> @author
 !> Henri Vincenti
 !
 !> @date
 !> Creation 2018
 ! _______________________________________________________________________________________
  SUBROUTINE get_local_tile_mem()
    USE constants, ONLY: num
    USE grid_tilemodule, ONLY: aofgrid_tiles
    USE particles, ONLY: species_parray
    USE particle_properties, ONLY : nspecies
    USE tile_params
    USE mem_status, ONLY: local_grid_tiles_mem, local_part_tiles_mem
    IMPLICIT NONE 
    INTEGER(idp) :: ispecies, ix,iy,iz
    TYPE(particle_species), POINTER :: curr_sp
    TYPE(particle_tile), POINTER :: curr

    ! Get local size of species_parray structure 
    local_part_tiles_mem=0._num
    DO ispecies =1, nspecies
      curr_sp => species_parray(ispecies)
      ! Sets tile spatial extents for current species
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr=>curr_sp%array_of_tiles(ix, iy, iz)
            ! If no particles in tile, allocate particle arrays
            IF (.NOT. curr%l_arrays_allocated) THEN
              CYCLE
            ELSE 
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_x)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_y)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_z)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_ux)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_uy)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_uz)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_ex)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_ey)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_ez)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_bx)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_by)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_bz)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%pid)
              local_part_tiles_mem=local_part_tiles_mem+SIZEOF(curr%part_gaminv)
            ENDIF
          END DO
        END DO 
      END DO
    END DO

    ! Get local size of aofgrid_tiles structure 
    local_grid_tiles_mem=0._num
    DO iz=1, ntilez! LOOP ON TILES
      DO iy=1, ntiley
        DO ix=1, ntilex
          local_grid_tiles_mem=local_grid_tiles_mem+                                   &
          SIZEOF(aofgrid_tiles(ix, iy, iz)%arr1)
          local_grid_tiles_mem=local_grid_tiles_mem+                                   &
          SIZEOF(aofgrid_tiles(ix, iy, iz)%arr2)
          local_grid_tiles_mem=local_grid_tiles_mem+                                   &
          SIZEOF(aofgrid_tiles(ix, iy, iz)%arr3)
        END DO
      END DO
    END DO! END LOOP ON TILES
  END SUBROUTINE get_local_tile_mem 

 ! _______________________________________________________________________________________
 !> @brief
 !> Subroutine that computes the total memory size in Bytes occupied by tile structures 
 !> (grid/particles) on all ranks. These memory sizes are reduced on rank 0
 !>  in the mem_status module variables: global_grid_tiles_mem (for grid) and 
 !> global_part_tiles_mem (for particles)
 !
 !> @author
 !> Henri Vincenti
 !
 !> @date
 !> Creation 2018
 ! _______________________________________________________________________________________
  SUBROUTINE get_global_tile_mem()
    USE constants, ONLY: isp
    USE shared_data, ONLY: errcode, comm
    USE mpi_type_constants, ONLY: mpidbl
    USE mem_status, ONLY: local_grid_tiles_mem, local_part_tiles_mem,                  &
    global_grid_tiles_mem, global_part_tiles_mem 
    IMPLICIT NONE 

    ! - Estimate total particle arrays memory (reduce on proc 0)
    CALL MPI_REDUCE(local_part_tiles_mem, global_part_tiles_mem,                       &
    1_isp, mpidbl, MPI_SUM, 0_isp,comm,errcode)
    ! - Estimate total grid tile arrays memory (reduce on proc 0)
    CALL MPI_REDUCE(local_grid_tiles_mem, global_grid_tiles_mem,                       &
    1_isp, mpidbl, MPI_SUM, 0_isp,comm,errcode)
  END SUBROUTINE get_global_tile_mem 

END MODULE tiling

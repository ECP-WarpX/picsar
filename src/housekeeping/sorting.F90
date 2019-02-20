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
! SORTING.F90
!
! Developer:
! Mathieu Lobet
!
! This file contains subroutines for the particle sorting.
!
! List of subroutines:
! - pxr_particle_sorting
! - particle_sorting_sub
! - pxr_particle_bin_sorting
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Module for the particle sorting
!>
!> @details
!> This module contains sorting algorithms for the particles
!>
!> @author
!> Mathieu Lobet
!>
!> @date
!> Creation 04/2016
!
! ________________________________________________________________________________________
MODULE sorting
  IMPLICIT NONE

  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> Main subroutine called to sort the particles in the Fortran PIC loop
  !
  !> @author
  !> Mathieu Lobet
  !
  !> @details
  !> This subroutine is called in the main loop of the Fortran kernel.
  !> This subroutine calls particle_sorting_sub() and times it.
  !
  !> @date
  !> 2016
  ! ______________________________________________________________________________________
  SUBROUTINE pxr_particle_sorting
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: num
    USE shared_data, ONLY: sorting_activated
    USE tiling
    USE time_stat, ONLY: localtimes, timestat_itstart
    IMPLICIT NONE

    !LOGICAL(lp)  :: verbose=.TRUE.
    REAL(num) :: tdeb, tend

    IF ((sorting_activated.gt.0)) THEN
      IF (it.ge.timestat_itstart) THEN
        tdeb=MPI_WTIME()
      ENDIF
      CALL particle_sorting_sub
      IF (it.ge.timestat_itstart) THEN
        tend=MPI_WTIME()
        localtimes(10) = localtimes(10) + (tend-tdeb)
      ENDIF
    ENDIF

  END SUBROUTINE

  ! ______________________________________________________________________________________
  ! particle_sorting_sub
  !
  !> @brief
  !> General subroutine for the particle sorting, used in Python PIC loop
  !>
  !> @details
  !> This subroutine is called in pxr_particle_sorting() used in the main loop
  !>
  !> @author
  !> Mathieu Lobet
  !
  !> @date 2016
  ! ______________________________________________________________________________________
  SUBROUTINE particle_sorting_sub
    USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
    USE mpi
    USE params, ONLY: it
    USE particle_properties, ONLY: nspecies, wpid
    USE particle_speciesmodule, ONLY: particle_species
    USE particle_tilemodule, ONLY: particle_tile
    USE particles, ONLY: species_parray
    USE picsar_precision, ONLY: idp, lp, num
    USE shared_data, ONLY: c_dim, dx, dy, dz, rank, sorting_dx, sorting_dy,          &
      sorting_dz, sorting_shiftx, sorting_shifty, sorting_shiftz, sorting_verbose,   &
      x, y, z
    USE tile_params, ONLY: ntilex, ntiley, ntilez
    USE tiling
    IMPLICIT NONE

    INTEGER(idp)                    :: ispecies, ix, iy, iz, count

    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER    :: curr_tile
    TYPE(grid_tile), POINTER        :: currg
    INTEGER(idp)                    :: nxc, nyc, nzc
    INTEGER(idp)                    :: nxjg, nyjg, nzjg
    LOGICAL(lp)                     :: isgathered=.FALSE.
    REAL(num)                       :: sxmin, symin, szmin

    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) SHARED(ntilex,      &
    !$OMP ntiley, ntilez, nspecies, species_parray, aofgrid_tiles, dx, dy, dz, it,    &
    !$OMP rank, sorting_shiftx, sorting_shifty, sorting_shiftz, sorting_dx,c_dim,     &
    !$OMP sorting_dy, sorting_dz, sorting_verbose) PRIVATE(ix, iy, iz, ispecies,      &
    !$OMP curr, curr_tile, currg, count, nxc, nyc, nzc, nxjg, nyjg, nzjg, isgathered, &
    !$OMP sxmin, symin, szmin)
    DO iz=1, ntilez! LOOP ON TILES
      DO iy=1, ntiley
        DO ix=1, ntilex
          curr=>species_parray(1)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
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
            curr_tile=>curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            IF (count .GT. 0) isgathered=.TRUE.
          END DO

          !if (rank.eq.0) print*, ix, iy, iz

          IF (isgathered) THEN
            currg=>aofgrid_tiles(ix, iy, iz)

            ! Loop over the species
            DO ispecies=1, nspecies

              curr=>species_parray(ispecies)

              ! If the sorting period > 0 and the current iteration corresponds
              ! to a multiple of the period
              IF                                                                      &
              ((it.ge.curr%sorting_start).AND.(curr%sorting_period.gt.0).AND.(MOD(it, &
              curr%sorting_period).eq.0)) THEN

              ! - Get current tile properties
              ! - Init current tile variables
              curr_tile=>curr%array_of_tiles(ix, iy, iz)
              count=curr_tile%np_tile(1)

              IF (count .EQ. 0) CYCLE

              ! Sorting algorithm inside the tiles
              IF(c_dim ==3) THEN
                sxmin = curr_tile%x_tile_min + sorting_shiftx
                symin = curr_tile%y_tile_min + sorting_shifty
                szmin = curr_tile%z_tile_min + sorting_shiftz

                CALL pxr_particle_bin_sorting(count, curr_tile%part_x,                &
                curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,                &
                curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv,          &
                curr_tile%pid, wpid, sxmin, symin, szmin, curr_tile%x_tile_max +      &
                sorting_dx, curr_tile%y_tile_max + sorting_dy, curr_tile%z_tile_max + &
                sorting_dz, sorting_dx, sorting_dy, sorting_dz)
              ELSE
                sxmin = curr_tile%x_tile_min + sorting_shiftx
                szmin = curr_tile%z_tile_min + sorting_shiftz

                CALL pxr_particle_bin_sorting_2d(count, curr_tile%part_x,              &
                curr_tile%part_z, curr_tile%part_ux,                                   &
                curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv,           &
                curr_tile%pid, wpid, sxmin,  szmin, curr_tile%x_tile_max +             &
                sorting_dx                         ,curr_tile%z_tile_max +             &
                sorting_dz, sorting_dx, sorting_dz)

              ENDIF
            ENDIF
          END DO! END LOOP ON SPECIES
        ENDIF
      END DO
    END DO
  END DO! END LOOP ON TILES
  !$OMP END PARALLEL DO

END SUBROUTINE particle_sorting_sub


! ______________________________________________________________________________________
! pxr_particle_bin_sorting
!
!>
!> @brief
!> Particle cell sorting subroutine using the bin sorting algorithm
!>
!> @details
!> This subroutine uses a bin sorting algorithm to sort particles
!> (including their property arrays).
!> Here, the bins corresponds to the cells of a cartesian array.
!> The cell size is specified by the user.
!>
!> @author
!> Mathieu Lobet
!
!> @date 2016
!
!> @param[in] np2 number of particles
!> @param[inout] xp, yp, zp particle positions
!> @param[inout] ux, uy, uz particle momenta
!> @param[inout] gam particle gamma factor
!> @param[inout] pid particle id
!> @param[inout] wpid particle weight
!> @param[in] xmin2, ymin2, zmin2 minimum point position on the local grid
!> @param[in] xmax2, ymax2, zmax2 maximum point position on the local grid
!> @param[in] dxf, dyf, dzf bin space steps
!
! ______________________________________________________________________________________
SUBROUTINE pxr_particle_bin_sorting(np2, xp, yp, zp, ux, uy, uz, gam, pid, wpid,    &
  xmin2, ymin2, zmin2, xmax2, ymax2, zmax2, dxf, dyf, dzf)
  USE picsar_precision, ONLY: idp, num
  implicit none
  integer(idp) :: ip, np2
  integer(idp) :: k, ic, nbhc
  integer(idp) :: ix, iy, iz
  integer(idp) :: nx3, ny3, nz3
  integer(idp) :: wpid
  real(num)    :: dxi, dyi, dzi
  real(num)    :: dxf, dyf, dzf
  real(num)    :: x2, y2, z2
  real(num)    :: xmin2, ymin2, zmin2
  real(num)    :: xmax2, ymax2, zmax2
  real(num), dimension(np2), intent(inout)      :: xp, yp, zp
  real(num), dimension(np2), intent(inout)      :: ux, uy, uz
  real(num), dimension(np2), intent(inout)      :: gam
  REAL(num), DIMENSION(np2, 1), intent(inout)    :: pid
  real(num), dimension(np2)                     :: xps, yps, zps
  real(num), dimension(np2)                     :: uxs, uys, uzs
  real(num), dimension(np2)                     :: gams
  real(num), dimension(np2, wpid)                :: pids
  integer(idp), dimension(np2)                  :: hcnb! Cell number
  integer(idp), dimension(:), allocatable        :: piihc! Particle indexes in the grid
  integer(idp), dimension(:), allocatable        :: nbppc! Number of particles per cells

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
  DO ip=1, np2
    x2 = (xp(ip)-xmin2)*dxi
    y2 = (yp(ip)-ymin2)*dyi
    z2 = (zp(ip)-zmin2)*dzi

    ix = floor(x2)
    iy = floor(y2)
    iz = floor(z2)

    ! Bin id
    hcnb(ip) = iz*nx3*ny3 + iy*nx3 + ix+1
#if defined(DEBUG)
    IF ((hcnb(ip) > nbhc).OR.(hcnb(ip)<1)) THEN
      print*, 'Bin id', ip, hcnb(ip), nbhc
      print*, 'Particle ix, iy, iz', ix, iy, iz
      print*, 'Particle x, y, z', xp(ip), yp(ip), zp(ip)
      print*, 'Particle x2, y2, z2', x2, y2, z2
      print*, 'xmin, ymin, zmin', xmin2, ymin2, zmin2
      print*, 'xmax, ymax, zmax', xmax2, ymax2, zmax2
      print*, 'Particle dx, dy, dz', dxi, dyi, dzi
      print*, 'Particle nx, ny, nz', nx3, ny3, nz3
      stop
    ENDIF
#endif
    ! We count the number of particles in each bin
    nbppc(hcnb(ip)) = nbppc(hcnb(ip))+1

  ENDDO

  ! Determine particle indexes in the bin grid
  k=0
  DO ic = 1, nbhc
    piihc(ic) = k
    k = k+nbppc(ic)
  END DO

  ! Sorting of the particles including their properties
  !if (rank.eq.0) print*, 'Counting sort, phase 3'
  DO ip=1, np2

    k = hcnb(ip)
    piihc(k) = piihc(k) + 1
    pids(piihc(k), :) = pid(ip, :)

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

  deallocate(piihc, nbppc)

END SUBROUTINE


SUBROUTINE pxr_particle_bin_sorting_2d(np2, xp, zp, ux, uy, uz, gam, pid, wpid,    &
  xmin2,  zmin2, xmax2,  zmax2, dxf,  dzf)
  USE picsar_precision, ONLY: idp, num
  implicit none
  integer(idp) :: ip, np2
  integer(idp) :: k, ic, nbhc
  integer(idp) :: ix,  iz
  integer(idp) :: nx3,  nz3
  integer(idp) :: wpid
  real(num)    :: dxi, dzi
  real(num)    :: dxf, dzf
  real(num)    :: x2, z2
  real(num)    :: xmin2,  zmin2
  real(num)    :: xmax2,  zmax2
  real(num), dimension(np2), intent(inout)      :: xp,  zp
  real(num), dimension(np2), intent(inout)      :: ux, uy, uz
  real(num), dimension(np2), intent(inout)      :: gam
  REAL(num), DIMENSION(np2, 1), intent(inout)    :: pid
  real(num), dimension(np2)                     :: xps, zps
  real(num), dimension(np2)                     :: uxs, uys, uzs
  real(num), dimension(np2)                     :: gams
  real(num), dimension(np2, wpid)                :: pids
  integer(idp), dimension(np2)                  :: hcnb! Cell number
  integer(idp), dimension(:), allocatable        :: piihc! Particle indexes in the grid
  integer(idp), dimension(:), allocatable        :: nbppc! Number of particles per cells

  ! Bin sizes
  dxi = 1./dxf
  dzi = 1./dzf

  nx3 = ceiling((xmax2-xmin2)*dxi)
  nz3 = ceiling((zmax2-zmin2)*dzi)

  ! Number of bins
  nbhc = nx3*nz3

  allocate(piihc(nbhc))
  allocate(nbppc(nbhc))

  hcnb = 0
  piihc = 0
  nbppc = 0

  ! Counting sort
  ! The criteria is the position in term of bin position
  DO ip=1, np2
    x2 = (xp(ip)-xmin2)*dxi
    z2 = (zp(ip)-zmin2)*dzi

    ix = floor(x2)
    iz = floor(z2)

    ! Bin id
    hcnb(ip) = iz*nx3 +  ix+1
#if defined(DEBUG)
    IF ((hcnb(ip) > nbhc).OR.(hcnb(ip)<1)) THEN
      print*, 'Bin id', ip, hcnb(ip), nbhc
      print*, 'Particle ix, iz', ix,  iz
      print*, 'Particle x,  z', xp(ip), zp(ip)
      print*, 'Particle x2, z2', x2,  z2
      print*, 'xmin,  zmin', xmin2,  zmin2
      print*, 'xmax,  zmax', xmax2,  zmax2
      print*, 'Particle dx, dz', dxi, dzi
      print*, 'Particle nx, nz', nx3, nz3
      stop
    ENDIF
#endif
    ! We count the number of particles in each bin
    nbppc(hcnb(ip)) = nbppc(hcnb(ip))+1

  ENDDO

  ! Determine particle indexes in the bin grid
  k=0
  DO ic = 1, nbhc
    piihc(ic) = k
    k = k+nbppc(ic)
  END DO

  ! Sorting of the particles including their properties
  !if (rank.eq.0) print*, 'Counting sort, phase 3'
  DO ip=1, np2

    k = hcnb(ip)
    piihc(k) = piihc(k) + 1
    pids(piihc(k), :) = pid(ip, :)

    xps(piihc(k)) = xp(ip)
    zps(piihc(k)) = zp(ip)

    uxs(piihc(k)) = ux(ip)
    uys(piihc(k)) = uy(ip)
    uzs(piihc(k)) = uz(ip)

    gams(piihc(k)) = gam(ip)
  END DO

  ! Copy back to the original arrays
  pid = pids
  xp = xps
  zp = zps
  ux = uxs
  uy = uys
  uz = uzs
  gam = gams

  deallocate(piihc, nbppc)

END SUBROUTINE pxr_particle_bin_sorting_2d


END MODULE sorting

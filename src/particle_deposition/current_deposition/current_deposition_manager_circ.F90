! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
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
! CURRENT_DEPOSITION_MANAGER_circ.F90
!
! Developers
! Jean-Luc Vay
!
! Description:
! This file contains subroutines for managing the current deposition in 3D.
!
! List of subroutines:
! - pxrdepose_currents_on_grid_jrjtjl
!
! ________________________________________________________________________________________

SUBROUTINE depose_jrjtjl(jr, jt, jl, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
  xmin, zmin, dt, dx, dz, nx, nz, nmodes, nxguard, nzguard, nox, noz, current_depo_algo)
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nx, nz, nmodes, nox, noz, nxguard, nzguard, current_depo_algo
  REAL(num), DIMENSION(-nxguard:nx+nxguard, -nzguard:nz+nzguard, 0:nmodes-1), intent(in out) :: jr, jt, jl
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dx, dz, xmin, zmin
  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(3), nvalid(3)
  nguard = (/ nxguard, nzguard, 0_idp /)
  nvalid = (/ nx+1, nz+1, nmodes /)

  call depose_jrjtjl_generic( jr, nguard, nvalid, &
                              jt, nguard, nvalid, &
                              jl, nguard, nvalid, nmodes, &
                              np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
                              xmin, zmin, dt, dx, dz, nox, noz, current_depo_algo)
END SUBROUTINE


! ________________________________________________________________________________________
!> @brief
!> Generic subroutines for current deposition, adapted for field
!> arrays having different sizes depending on their nodal/cell-centered nature
!>
!> @details
!> This routine calls the relevant current deposition routine depending
!> on the order of the particle shape and the selected algorithm.
!>
! ________________________________________________________________________________________
SUBROUTINE depose_jrjtjl_generic( jr, jr_nguard, jr_nvalid, &
                                  jt, jt_nguard, jt_nvalid, &
                                  jl, jl_nguard, jl_nvalid, nmodes, &
                                  np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q, &
                                  xmin, zmin, dt, dr, dz, &
                                  nox, noz, current_depo_algo)     !#do not wrap
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp) :: np, nox, noz, current_depo_algo
  INTEGER(idp), intent(in) :: jr_nguard(2), jr_nvalid(2), &
                              jt_nguard(2), jt_nvalid(2), &
                              jl_nguard(2), jl_nvalid(2), nmodes
  REAL(num), intent(IN OUT):: jr(-jr_nguard(1):jr_nvalid(1)+jr_nguard(1)-1, &
                                 -jr_nguard(2):jr_nvalid(2)+jr_nguard(2)-1, &
                                 0:nmodes-1 )
  REAL(num), intent(IN OUT):: jt(-jt_nguard(1):jt_nvalid(1)+jt_nguard(1)-1, &
                                 -jt_nguard(2):jt_nvalid(2)+jt_nguard(2)-1, &
                                 0:nmodes-1 )
  REAL(num), intent(IN OUT):: jl(-jl_nguard(1):jl_nvalid(1)+jl_nguard(1)-1, &
                                 -jl_nguard(2):jl_nvalid(2)+jl_nguard(2)-1, &
                                 0:nmodes-1 )
  REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
  REAL(num) :: q, dt, dr, dz, xmin, zmin


  Select CASE(current_depo_algo)
    ! Scalar classical current deposition subroutines
 
   CASE DEFAULT

    IF ((nox.eq.1).and.(noz.eq.1)) THEN
      CALL depose_jxjyjz_esirkepov_1_1_1( jr, jr_nguard, jr_nvalid, &
                                          jt, jt_nguard, jt_nvalid, &
                                          jl, jl_nguard, jl_nvalid, nmodes, &
                                          np, xp, yp, zp, uxp, uyp, uzp, gaminv, w,  &
                                          q, xmin, zmin, dt, dr, dz)
 
    ENDIF

  END SELECT

END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> Main subroutine for managing the current deposition across tiles
!
!> @details
!> This subroutine is called in submain.F90 in step().
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> 2015-2016
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jrjtjl
  USE fields, ONLY: nox, noy, jr, jt, jl, nxjguards, nyjguards
  USE mpi
  USE params, ONLY: currdepo, lvec_curr_depo, dt, it
  USE particle_properties, ONLY: nspecies
  USE picsar_precision, ONLY: idp, num, cpx
  USE shared_data, ONLY: nx, ny, nmodes, xmin, ymin, dx, dy
  USE time_stat, ONLY: timestat_itstart, localtimes
  IMPLICIT NONE
  REAL(num) :: tdeb, tend

  ! ___________________________________________________________________________
  ! Interfaces for func_order
  INTERFACE
    ! ____________________________________________________________________________________
    ! Generic current deposition routine
    SUBROUTINE depose_jrjtjl(jr, jt, jl, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,     &
      xmin, zmin, dt, dx, dz, nx, nz, nmodes, nxguard, nzguard, nox, noz, current_depo_algo)
      USE picsar_precision, ONLY: idp, num, cpx
      IMPLICIT NONE
      INTEGER(idp) :: np, nx, nz, nmodes, nox, noz, nxguard, nzguard, current_depo_algo
      COMPLEX(cpx), DIMENSION(-nxguard:nx+nxguard, -nzguard:nz+nzguard, 0:nmodes-1), intent(in out) :: jr, jt, jl
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dx, dz, xmin, zmin
    END SUBROUTINE


  END INTERFACE
  ! ___________________________________________________________________________

#if defined(DEBUG)
  WRITE(0, *) "Depose_currents_on_grid: start"
#endif

  IF (it.ge.timestat_itstart) THEN
    tdeb=MPI_WTIME()
  ENDIF

#if VTUNE==2
  CALL start_vtune_collection()
#endif
#if SDE==2
  CALL start_vtune_collection()
#endif

  IF (nspecies .EQ. 0_idp) RETURN

  ! RZ - CIRC
  jr = 0.0_cpx
  jt = 0.0_cpx
  jl = 0.0_cpx

  CALL pxrdepose_currents_on_grid_jrjtjl_classical_sub_seq(depose_jrjtjl, jr, jt, &
  jl, nx, ny, nxjguards, nyjguards, nox, noy, dx, dy, dt, 3_idp)

  !!! --- Stop Vtune analysis
#if VTUNE==2
  CALL stop_vtune_collection()
#endif
#if SDE==2
  CALL stop_sde_collection()
#endif

  IF (it.ge.timestat_itstart) THEN
    tend = MPI_WTIME()
    localtimes(3)=localtimes(3)+(tend-tdeb)
  ENDIF

#if defined(DEBUG)
  WRITE(0, *) "Depose_current_on_grid: stop"
#endif

END SUBROUTINE pxrdepose_currents_on_grid_jrjtjl

! ________________________________________________________________________________________
!> @brief
!> Deposit current in each tile sequentially
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrdepose_currents_on_grid_jrjtjl_classical_sub_seq(func_order, jrg, jtg,  &
  jlg, nr, nl, nrjguard, nljguard, nor, nol, dr, dl,    &
  dtt, currrent_depo_algo )
  USE grid_tilemodule, ONLY: grid_tile, aofgrid_tiles
  USE particle_properties, ONLY: nspecies, wpid
  USE particle_speciesmodule, ONLY: particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, num, lp
  USE tile_params, ONLY: ntilez, ntilex, ntiley
  USE tiling

  IMPLICIT NONE
  INTEGER(idp), INTENT(IN) :: nr, nl, nrjguard, nljguard
  INTEGER(idp), INTENT(IN) :: nor, nol
  REAL(num), INTENT(IN) :: dr, dl, dtt
  COMPLEX(cpx), INTENT(IN OUT) :: jrg(-nrjguard:nr+nrjguard, &
                                   -nljguard:nl+nljguard, 0:nmodes-1)
  COMPLEX(cpx), INTENT(IN OUT) :: jtg(-nrjguard:nr+nrjguard, &
                                   -nljguard:nl+nljguard, 0:nmodes-1)
  COMPLEX(cpx), INTENT(IN OUT) :: jlg(-nrjguard:nr+nrjguard, &
                                   -nljguard:nl+nljguard, 0:nmodes-1)
  INTEGER(idp)                    :: ispecies, ix, iy, count
  INTEGER(idp)                    :: currrent_depo_algo
  INTEGER(idp)                    :: jmin, jmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(particle_tile), POINTER    :: curr_tile
  TYPE(grid_tile), POINTER        :: currg
  INTEGER(idp)                    :: nrc, nlc, nrjg, nljg
  LOGICAL(lp)                     :: isdeposited=.FALSE.

  ! Interfaces for func_order
  INTERFACE
    SUBROUTINE func_order(jr, jt, jl, np, xp, yp, zp, uxp, uyp, uzp, gaminv, w, q,    &
      rmin, zmin, dt, dr, dl, nr, nl, nmodes, nrguard, nlguard, nor, nol, current_depo_algo)  !#do not parse
      USE PICSAR_precision
      USE constants
      IMPLICIT NONE
      INTEGER(idp) :: np, nr, nl, nmodes, nor, nol, nrguard, nlguard,       &
      current_depo_algo
      COMPLEX(cpx), DIMENSION(-nrguard:nr+nrguard, &
                              -nlguard:nl+nlguard, &
                              0:nmodes-1), intent(in out) :: jr, jt, jl
      REAL(num), DIMENSION(np) :: xp, yp, zp, uxp, uyp, uzp, w, gaminv
      REAL(num) :: q, dt, dr, dl, rmin, zmin
    END SUBROUTINE
  END INTERFACE


  IF (nspecies .EQ. 0_idp) RETURN
  DO iy=1, ntiley
      DO ix=1, ntilex
        curr => species_parray(1)
        curr_tile=>curr%array_of_tiles(ix, iy, 1)
        nrjg=curr_tile%nxg_tile
        nljg=curr_tile%nyg_tile
        jmin=curr_tile%nx_tile_min-nrjg
        jmax=curr_tile%nx_tile_max+nrjg
        lmin=curr_tile%ny_tile_min-nljg
        lmax=curr_tile%ny_tile_max+nljg
        nrc=curr_tile%nx_cells_tile; 
        nlc=curr_tile%ny_cells_tile
        currg=>aofgrid_tiles(ix, iy, 1)
        currg%carr1=0.
        currg%carr2=0.
        currg%carr3=0.
        isdeposited=.FALSE.
        DO ispecies=1, nspecies! LOOP ON SPECIES
          curr => species_parray(ispecies)
          IF (.NOT. curr%ldodepos) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, 1)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) THEN
            CYCLE
          ELSE
            isdeposited=.TRUE.
          ENDIF
          ! Depose current in jtile
          CALL func_order(currg%carr1, currg%carr2, currg%carr3, count,            &
          curr_tile%part_x, curr_tile%part_y, curr_tile%part_z, curr_tile%part_ux,    &
          curr_tile%part_uy, curr_tile%part_uz, curr_tile%part_gaminv,                &
          curr_tile%pid(1, wpid), curr%charge, curr_tile%x_grid_tile_min,             &
          curr_tile%y_grid_tile_min, dtt, dr, dl,   &
          nrc, nlc, nmodes, nrjg, nljg, nor, nol, currrent_depo_algo)
        END DO! END LOOP ON SPECIES
        IF (isdeposited) THEN
          jrg(jmin:jmax, lmin:lmax, 1:nmodes) = & 
          jrg(jmin:jmax, lmin:lmax, 1:nmodes) + currg%arr1
          
          jtg(jmin:jmax, lmin:lmax, 1:nmodes) = &
          jtg(jmin:jmax, lmin:lmax, 1:nmodes) + currg%arr2
          
          jlg(jmin:jmax, lmin:lmax, 1:nmodes) = &
          jlg(jmin:jmax, lmin:lmax, 1:nmodes) + currg%arr3
        ENDIF
      END DO
  END DO!END LOOP ON TILES


END SUBROUTINE pxrdepose_currents_on_grid_jrjtjl_classical_sub_seq



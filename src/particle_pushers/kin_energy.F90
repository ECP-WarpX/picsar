! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
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
! PARTICLE_PUSHER_MANAGER.F90
!
! Subroutines for computing particle kin_energy.
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @author
!> Haithem Kallala

!> @date
!> Creation 2017
! ________________________________________________________________________________________
SUBROUTINE compute_kin_energy
  USE grid_tilemodule, ONLY: aofgrid_tiles, grid_tile
  USE mpi
  USE params, ONLY: fieldgathe, lvec_fieldgathe, nc
  USE particle_properties, ONLY: nspecies
  USE particle_speciesmodule, ONLY: kin_energy_mpi, kin_energy_total,                &
    particle_species
  USE particle_tilemodule, ONLY: particle_tile
  USE particles, ONLY: species_parray
  USE picsar_precision, ONLY: idp, isp, num
  USE shared_data, ONLY: comm, dx, dy, dz, errcode
  USE tile_params, ONLY: ntilex, ntiley, ntilez
  IMPLICIT NONE
  INTEGER(idp)             :: ispecies, ix, iy, iz, count
  INTEGER(idp)             :: jmin, jmax, kmin, kmax, lmin, lmax
  TYPE(particle_species), POINTER :: curr
  TYPE(grid_tile), POINTER        :: currg
  TYPE(particle_tile), POINTER    :: curr_tile
  REAL(num)                :: tdeb, tend, ppid
  INTEGER(idp)             :: nxc, nyc, nzc, ipmin, ipmax, ip
#if defined(DEBUG)
  WRITE(0, *) "compute_kin_energy : start"
#endif

  IF (nspecies .EQ. 0_idp) RETURN
  DO ispecies=1, nspecies! LOOP ON SPECIES
    curr=>species_parray(ispecies)
    curr%kin_energy_sp=0.0_num
  ENDDO

  !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime) DEFAULT(NONE) SHARED(dx, dy, dz,    &
  !$OMP nc, ntilex, ntiley, ntilez, nspecies, species_parray, aofgrid_tiles,          &
  !$OMP LVEC_fieldgathe) PRIVATE(ix, iy, iz, ispecies, curr, curr_tile, currg, count, &
  !$OMP ppid)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        DO ispecies=1, nspecies
          curr=>species_parray(ispecies)
          IF (curr%is_antenna) CYCLE
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) CYCLE
          ppid = nc*dx*dy*dz/(curr%nppcell)
          CALL compute_kin_energy_vector(count, curr_tile%part_gaminv, curr%mass,     &
          curr_tile%kin_energy_tile, ppid)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        DO ispecies=1, nspecies
          curr=>species_parray(ispecies)
          curr_tile=>curr%array_of_tiles(ix, iy, iz)
          count=curr_tile%np_tile(1)
          IF (count .EQ. 0) CYCLE
          curr%kin_energy_sp = curr%kin_energy_sp + curr_tile%kin_energy_tile
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  kin_energy_mpi = 0.0_num
  DO ispecies=1, nspecies! LOOP ON SPECIES
    curr=>species_parray(ispecies)
    kin_energy_mpi=kin_energy_mpi+curr%kin_energy_sp
    curr%kin_energy_sp = 0.0_num
  ENDDO
  kin_energy_total = 0.0_num
  CALL MPI_ALLREDUCE(kin_energy_mpi, kin_energy_total, 1_isp, MPI_DOUBLE, MPI_SUM,    &
  comm, errcode)

END SUBROUTINE compute_kin_energy

SUBROUTINE compute_kin_energy_vector(np, gaminv, mass, kin_e, ppid)
  USE constants, ONLY: clight
  USE picsar_precision, ONLY: idp, num
  IMPLICIT NONE
  INTEGER(idp), INTENT(IN)   :: np
  REAL(num), INTENT(IN)   :: gaminv(np), mass, ppid
  REAL(num), INTENT(INOUT)   :: kin_e
  REAL(num)                ::  mclightsq
  INTEGER(idp)             :: ip
  kin_e = 0.0_num
  mclightsq=mass*clight**2*ppid
#if !defined PICSAR_NO_ASSUMED_ALIGNMENT && defined __INTEL_COMPILER
  !DIR$ ASSUME_ALIGNED gaminv:64
#elif defined __IBMBGQ__
  !IBM* ALIGN(64, gaminv)
#endif
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP SIMD
#endif
#elif defined __IBMBGQ__
  !IBM* SIMD_LEVEL
#elif defined __INTEL_COMPILER
  !DIR$ SIMD
#endif
  DO ip=1, np
    kin_e = kin_e + (1.0_num/gaminv(ip)-1.0_num)*mclightsq
  ENDDO
#if defined _OPENMP && _OPENMP>=201307
#ifndef NOVEC
  !$OMP END SIMD
#endif
#endif

END SUBROUTINE compute_kin_energy_vector

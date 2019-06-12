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
! FIELD_GATHERING_MANAGER_CIRC.F90
!
! This file comtains subroutines to manage the field gathering in 2D.
!
! Developers:
! - Jean-Luc Vay
! - David Grote
!
! List of subroutines:
! - geteb2dcirc_energy_conserving
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> General subroutines for the RZ cylindrical field gathering
!
!> @author
!> Jean-Luc Vay
!> David Grote
!
!> @date
!> Creation 2016
!
!> @param[in] np Number of particles
!> @param[in] xp, yp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nx, nz space discretization
!> @param[in] nmodes number of cylindrical modes (including mode 0)
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] erg, etg, elg field arrays
!> @param[in] brg, btg, blg field arrays
!
! ________________________________________________________________________________________
SUBROUTINE geteb2dcirc_energy_conserving(np, xp, yp, zp, ex, ey, ez, bx, by, bz, &
                                         xmin, zmin, dx, dz, nx, nz, nmodes, nxguard, nzguard, nox, noz, &
                                         l_lower_order_in_v, l_nodal, &
                                         erg, etg, elg, brg, btg, blg)
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp)                  :: np, nx, nz, nmodes, nox, noz, nxguard, nzguard
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  logical(lp)                   :: l_lower_order_in_v, l_nodal
  complex(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard, 0:nmodes-1) :: erg, etg, elg
  complex(num), dimension(-nxguard:nx+nxguard, -nzguard:nz+nzguard, 0:nmodes-1) :: brg, btg, blg
  real(num)                     :: xmin, zmin, dx, dz

  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(2), nvalid(2)
  nguard = (/ nxguard, nzguard /)
  nvalid = (/ nx+1, nz+1 /)

  call geteb2dcirc_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by, bz,    &
                     xmin, zmin, dx, dz, nmodes, nox, noz, l_lower_order_in_v, l_nodal, &
                     erg, nguard, nvalid, &
                     etg, nguard, nvalid, &
                     elg, nguard, nvalid, &
                     brg, nguard, nvalid, &
                     btg, nguard, nvalid, &
                     blg, nguard, nvalid)
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> General subroutines for the RZ field gathering, adapted for field
!> arrays having different sizes depending on their nodal/cell-centered nature
!>
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] erg_nguard, etg_nguard, elg_nguard number of guard cells of the
!> erg, etg, elg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] erg_nvalid, etg_nvalid, elg_nvalid number of valid gridpoints
!>  (i.e. not guard cells) of the erg, etg, elg arrays (1d arrays containing 2 integers)
!> @param[in] brg, btg, blg magnetic field grids
!> @param[in] brg_nguard, btg_nguard, blg_nguard number of guard cells of the
!> brg, btg, blg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] brg_nvalid, btg_nvalid, blg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the brg, btg, blg arrays (1d arrays containing 2 integers)
!> @param[in] erg, etg, elg electric field grid
!>
! ________________________________________________________________________________________
SUBROUTINE geteb2dcirc_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by, bz, &
                     xmin, zmin, dx, dz, nmodes, nox, noz, l_lower_order_in_v, l_nodal, &
                     erg, erg_nguard, erg_nvalid, etg, etg_nguard, etg_nvalid, elg, elg_nguard, elg_nvalid, &
                     brg, brg_nguard, brg_nvalid, btg, btg_nguard, btg_nvalid, blg, blg_nguard, blg_nvalid) !#do not wrap
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp), intent(IN)      :: np, nmodes, nox, noz
  logical(lp), intent(IN)       :: l_lower_order_in_v, l_nodal
  integer(idp), intent(IN)      :: erg_nguard(2), erg_nvalid(2), etg_nguard(2)
  integer(idp), intent(IN)      :: etg_nvalid(2), elg_nguard(2), elg_nvalid(2)
  integer(idp), intent(IN)      :: brg_nguard(2), brg_nvalid(2), btg_nguard(2)
  integer(idp), intent(IN)      :: btg_nvalid(2), blg_nguard(2), blg_nvalid(2)
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  real(num), intent(IN)         :: xmin, zmin, dx, dz
  complex(num), intent(IN):: erg(-erg_nguard(1):erg_nvalid(1)+erg_nguard(1)-1,           &
                                 -erg_nguard(2):erg_nvalid(2)+erg_nguard(2)-1, 0:nmodes-1)
  complex(num), intent(IN):: etg(-etg_nguard(1):etg_nvalid(1)+etg_nguard(1)-1,           &
                                 -etg_nguard(2):etg_nvalid(2)+etg_nguard(2)-1, 0:nmodes-1)
  complex(num), intent(IN):: elg(-elg_nguard(1):elg_nvalid(1)+elg_nguard(1)-1,           &
                                 -elg_nguard(2):elg_nvalid(2)+elg_nguard(2)-1, 0:nmodes-1)
  complex(num), intent(IN):: brg(-brg_nguard(1):brg_nvalid(1)+brg_nguard(1)-1,           &
                                 -brg_nguard(2):brg_nvalid(2)+brg_nguard(2)-1, 0:nmodes-1)
  complex(num), intent(IN):: btg(-btg_nguard(1):btg_nvalid(1)+btg_nguard(1)-1,           &
                                 -btg_nguard(2):btg_nvalid(2)+btg_nguard(2)-1, 0:nmodes-1)
  complex(num), intent(IN):: blg(-blg_nguard(1):blg_nvalid(1)+blg_nguard(1)-1,           &
                                 -blg_nguard(2):blg_nvalid(2)+blg_nguard(2)-1, 0:nmodes-1)

  ! ______________________________________________
  ! Arbitrary order, non-optimized subroutines

  !!! --- Gather electric field on particles
  CALL pxr_gete2drz_n_energy_conserving(np, xp, yp, zp, ex, ey, ez, &
                                        xmin, zmin, dx, dz, nmodes, nox, noz, &
                                        erg, erg_nguard, erg_nvalid, &
                                        etg, etg_nguard, etg_nvalid, &
                                        elg, elg_nguard, elg_nvalid, &
                                        l_lower_order_in_v, l_nodal)

  !!! --- Gather magnetic fields on particles
  CALL pxr_getb2drz_n_energy_conserving(np, xp, yp, zp, bx, by, bz, &
                                        xmin, zmin, dx, dz, nmodes, nox, noz, &
                                        brg, brg_nguard, brg_nvalid, &
                                        btg, btg_nguard, btg_nvalid, &
                                        blg, blg_nguard, blg_nvalid, &
                                        l_lower_order_in_v, l_nodal)

END SUBROUTINE

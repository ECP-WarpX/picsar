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
! FIELD_GATHERING_MANAGER_2D.F90
!
! This file comtains subroutines to manage the field gathering in 2D.
!
! Developers:
! - Henri vincenti
! - Mathieu Lobet
!
! List of subroutines:
! - geteb2dxz_energy_conserving
!
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> General subroutines for the 2D cartesian field gathering
!
!> @author
!> Mathieu Lobet
!
!> @date
!> Creation 2016
!
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] nx, nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] exg, eyg, ezg field arrays
!> @param[in] bxg, byg, bzg field arrays
!> @param[in] l4symetry
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!> @param[in] field_gathe_algo Gathering algorithm
!> @param[in] lvect vector length
!
! ________________________________________________________________________________________
SUBROUTINE geteb2dxz_energy_conserving(np, xp, yp, zp, ex, ey, ez, bx, by, bz, xmin,  &
  ymin, zmin, dx, dy, dz, nx, ny, nz, nxguard, nyguard, nzguard, nox, noy, noz, exg,    &
  eyg, ezg, bxg, byg, bzg, l4symtry, l_lower_order_in_v, lvect, field_gathe_algo)
  USE picsar_precision, ONLY: idp, lp, num
  USE fields, ONLY: l_nodalgrid
  implicit none

  integer(idp)                  :: np, nx, ny, nz, nox, noy, noz, nxguard, nyguard,   &
  nzguard
  integer(idp)                  :: field_gathe_algo
  integer(idp)                  :: lvect
  logical(lp) , intent(in)      :: l4symtry, l_lower_order_in_v
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  real(num), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: exg, eyg, ezg
  real(num), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard, -nzguard:nz+nzguard) &
  :: bxg, byg, bzg
  real(num)                     :: xmin, ymin, zmin, dx, dy, dz

  ! Build array of guard cells and valid cells, to pass them to the generic routine
  integer(idp)                       :: nguard(2), nvalid(2)
  nguard = (/ nxguard, nzguard /)
  nvalid = (/ nx+1, nz+1 /)

  call geteb2dxz_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by, bz,    &
  xmin, ymin, zmin, dx, dy, dz, nox, noy, noz, exg, nguard, nvalid, eyg, nguard,      &
  nvalid, ezg, nguard, nvalid, bxg, nguard, nvalid, byg, nguard, nvalid, bzg, nguard, &
  nvalid, l4symtry, l_lower_order_in_v, l_nodalgrid, lvect, field_gathe_algo)
END SUBROUTINE

! ________________________________________________________________________________________
!> @brief
!> General subroutines for the 2D field gathering, adapted for field
!> arrays having different sizes depending on their nodal/cell-centered nature
!>
!> @param[in] np Number of particles
!> @param[in] xp, zp particle position arrays
!> @param[inout] ex, ey, ez electric field particle arrays
!> @param[inout] bx, by, bz magnetic field particle arrays
!> @param[in] xmin, zmin tile boundaries
!> @param[in] dx, dz space steps
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!> exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!>  (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the
!> bxg, byg, bzg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 2 integers)
!> @param[in] l4symetry
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!> @param[in] field_gathe_algo Gathering algorithm
!> @param[in] lvect vector length
!> @param[in] exg, eyg, ezg electric field grid
!>
! ________________________________________________________________________________________
SUBROUTINE geteb2dxz_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by,    &
  bz, xmin, ymin, zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg,    &
  eyg_nguard, eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid,     &
  byg, byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l4symtry,                   &
  l_lower_order_in_v, l_nodal, lvect, field_gathe_algo)            !#do not wrap
  USE params, ONLY: fieldgathe, lvec_fieldgathe
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp)                  :: field_gathe_algo
  integer(idp)                  :: np, nox, noy, noz
  integer(idp), intent(IN)      :: exg_nguard(2), exg_nvalid(2), eyg_nguard(2),       &
  eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2), bxg_nguard(2), bxg_nvalid(2),          &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  LOGICAL(lp), intent(in)      :: l4symtry, l_lower_order_in_v, l_nodal
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  real(num)                     :: xmin, ymin, zmin, dx, dy, dz
  integer(idp)                  :: lvect
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)

  IF (field_gathe_algo.lt.0) return

  ! ______________________________________________
  ! Arbitrary order, non-optimized subroutines
  IF (field_gathe_algo.eq.2) THEN


    !!! --- Gather electric field on particles
    CALL pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,    &
    dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,  &
    ezg_nguard, ezg_nvalid, l4symtry, .FALSE._idp, l_lower_order_in_v, l_nodal)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,    &
    dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,  &
    bzg_nguard, bzg_nvalid, l4symtry, .FALSE._idp, l_lower_order_in_v, l_nodal)

    ! ______________________________________________
    ! Arbitrary order, scalar subroutines
  ELSE IF (field_gathe_algo.eq.1) THEN

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_scalar_1_1( np, xp, zp, ex, ey, ez, xmin,   &
      zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
      ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_scalar_1_1( np, xp, zp, bx, by, bz, xmin,   &
      zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
      bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_scalar_2_2( np, xp, zp, ex, ey, ez, xmin,   &
      zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
      ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_scalar_2_2( np, xp, zp, bx, by, bz, xmin,   &
      zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
      bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_scalar_3_3( np, xp, zp, ex, ey, ez, xmin,   &
      zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
      ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_scalar_3_3( np, xp, zp, bx, by, bz, xmin,   &
      zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
      bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

      ! Arbitrary order
    ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,  &
      dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
      ezg, ezg_nguard, ezg_nvalid, l4symtry, .FALSE._idp, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,  &
      dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
      bzg, bzg_nguard, bzg_nvalid, l4symtry, .FALSE._idp, l_lower_order_in_v, l_nodal)
    ENDIF

    ! ________________________________________
    ! Optimized subroutines, default
  ELSE


    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_vect_1_1( np, xp, zp, ex, ey, ez, xmin,     &
      zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
      ezg_nguard, ezg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_vect_1_1( np, xp, zp, bx, by, bz, xmin,     &
      zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
      bzg_nguard, bzg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_vect_2_2( np, xp, zp, ex, ey, ez, xmin,     &
      zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
      ezg_nguard, ezg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_vect_2_2( np, xp, zp, bx, by, bz, xmin,     &
      zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
      bzg_nguard, bzg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      !!! --- Gather electric and magnetic field on particles
      CALL pxr_geteb2dxz_energy_conserving_vect_3_3( np, xp, zp, ex, ey, ez, bx, by,  &
      bz, xmin, zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,           &
      eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,      &
      byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, LVEC_fieldgathe,           &
      l_lower_order_in_v, l_nodal)

      ! Arbitrary order
    ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,  &
      dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
      ezg, ezg_nguard, ezg_nvalid, l4symtry, .FALSE._idp, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,  &
      dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
      bzg, bzg_nguard, bzg_nvalid, l4symtry, .FALSE._idp, l_lower_order_in_v, l_nodal)
    ENDIF
  ENDIF
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
!> @param[in] exg_nguard, eyg_nguard, ezg_nguard number of guard cells of the
!> exg, eyg, ezg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] exg_nvalid, eyg_nvalid, ezg_nvalid number of valid gridpoints
!>  (i.e. not guard cells) of the exg, eyg, ezg arrays (1d arrays containing 2 integers)
!> @param[in] bxg, byg, bzg magnetic field grids
!> @param[in] bxg_nguard, byg_nguard, bzg_nguard number of guard cells of the
!> bxg, byg, bzg arrays in each direction (1d arrays containing 2 integers)
!> @param[in] bxg_nvalid, byg_nvalid, bzg_nvalid number of valid gridpoints
!> (i.e. not guard cells) of the bxg, byg, bzg arrays (1d arrays containing 2 integers)
!> @param[in] l4symetry
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!> @param[in] field_gathe_algo Gathering algorithm
!> @param[in] lvect vector length
!> @param[in] exg, eyg, ezg electric field grid
!>
! ________________________________________________________________________________________
SUBROUTINE geteb2drz_energy_conserving_generic(np, xp, yp, zp, ex, ey, ez, bx, by,    &
  bz, xmin, ymin, zmin, dx, dy, dz, nox, noy, noz, exg, exg_nguard, exg_nvalid, eyg,    &
  eyg_nguard, eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid,     &
  byg, byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, l4symtry,                   &
  l_lower_order_in_v, l_nodal, lvect, field_gathe_algo)            !#do not wrap
  USE params, ONLY: fieldgathe, lvec_fieldgathe
  USE picsar_precision, ONLY: idp, lp, num
  implicit none

  integer(idp)                  :: field_gathe_algo
  integer(idp)                  :: np, nox, noy, noz
  integer(idp), intent(IN)      :: exg_nguard(2), exg_nvalid(2), eyg_nguard(2),       &
  eyg_nvalid(2), ezg_nguard(2), ezg_nvalid(2), bxg_nguard(2), bxg_nvalid(2),          &
  byg_nguard(2), byg_nvalid(2), bzg_nguard(2), bzg_nvalid(2)
  LOGICAL(lp), intent(in)      :: l4symtry, l_lower_order_in_v, l_nodal
  real(num), dimension(np)      :: xp, yp, zp, ex, ey, ez, bx, by, bz
  real(num)                     :: xmin, ymin, zmin, dx, dy, dz
  integer(idp)                  :: lvect
  REAL(num), intent(IN):: exg(-exg_nguard(1):exg_nvalid(1)+exg_nguard(1)-1,           &
  -exg_nguard(2):exg_nvalid(2)+exg_nguard(2)-1)
  REAL(num), intent(IN):: eyg(-eyg_nguard(1):eyg_nvalid(1)+eyg_nguard(1)-1,           &
  -eyg_nguard(2):eyg_nvalid(2)+eyg_nguard(2)-1)
  REAL(num), intent(IN):: ezg(-ezg_nguard(1):ezg_nvalid(1)+ezg_nguard(1)-1,           &
  -ezg_nguard(2):ezg_nvalid(2)+ezg_nguard(2)-1)
  REAL(num), intent(IN):: bxg(-bxg_nguard(1):bxg_nvalid(1)+bxg_nguard(1)-1,           &
  -bxg_nguard(2):bxg_nvalid(2)+bxg_nguard(2)-1)
  REAL(num), intent(IN):: byg(-byg_nguard(1):byg_nvalid(1)+byg_nguard(1)-1,           &
  -byg_nguard(2):byg_nvalid(2)+byg_nguard(2)-1)
  REAL(num), intent(IN):: bzg(-bzg_nguard(1):bzg_nvalid(1)+bzg_nguard(1)-1,           &
  -bzg_nguard(2):bzg_nvalid(2)+bzg_nguard(2)-1)

  IF (field_gathe_algo.lt.0) return

  ! ______________________________________________
  ! Arbitrary order, non-optimized subroutines
  IF (field_gathe_algo.eq.2) THEN


    !!! --- Gather electric field on particles
    CALL pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,    &
    dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,  &
    ezg_nguard, ezg_nvalid, l4symtry, .TRUE._idp, l_lower_order_in_v, l_nodal)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,    &
    dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,  &
    bzg_nguard, bzg_nvalid, l4symtry, .TRUE._idp, l_lower_order_in_v, l_nodal)

    ! ______________________________________________
    ! Arbitrary order, scalar subroutines
  ELSE IF (field_gathe_algo.eq.1) THEN

!   IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

!     !!! --- Gather electric field on particles
!     CALL pxr_gete2dxz_energy_conserving_scalar_1_1( np, xp, zp, ex, ey, ez, xmin,   &
!     zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
!     ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
!     !!! --- Gather magnetic fields on particles
!     CALL pxr_getb2dxz_energy_conserving_scalar_1_1( np, xp, zp, bx, by, bz, xmin,   &
!     zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
!     bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

!   ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

!     !!! --- Gather electric field on particles
!     CALL pxr_gete2dxz_energy_conserving_scalar_2_2( np, xp, zp, ex, ey, ez, xmin,   &
!     zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
!     ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
!     !!! --- Gather magnetic fields on particles
!     CALL pxr_getb2dxz_energy_conserving_scalar_2_2( np, xp, zp, bx, by, bz, xmin,   &
!     zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
!     bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

!   ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

!     !!! --- Gather electric field on particles
!     CALL pxr_gete2dxz_energy_conserving_scalar_3_3( np, xp, zp, ex, ey, ez, xmin,   &
!     zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
!     ezg_nguard, ezg_nvalid, l_lower_order_in_v, l_nodal)
!     !!! --- Gather magnetic fields on particles
!     CALL pxr_getb2dxz_energy_conserving_scalar_3_3( np, xp, zp, bx, by, bz, xmin,   &
!     zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
!     bzg_nguard, bzg_nvalid, l_lower_order_in_v, l_nodal)

!     ! Arbitrary order
!   ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,  &
      dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
      ezg, ezg_nguard, ezg_nvalid, l4symtry, .TRUE._idp, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,  &
      dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
      bzg, bzg_nguard, bzg_nvalid, l4symtry, .TRUE._idp, l_lower_order_in_v, l_nodal)
!   ENDIF

    ! ________________________________________
    ! Optimized subroutines, default
  ELSE


!   IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

!     !!! --- Gather electric field on particles
!     CALL pxr_gete2dxz_energy_conserving_vect_1_1( np, xp, zp, ex, ey, ez, xmin,     &
!     zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
!     ezg_nguard, ezg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
!     !!! --- Gather magnetic fields on particles
!     CALL pxr_getb2dxz_energy_conserving_vect_1_1( np, xp, zp, bx, by, bz, xmin,     &
!     zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
!     bzg_nguard, bzg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)

!   ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

!     !!! --- Gather electric field on particles
!     CALL pxr_gete2dxz_energy_conserving_vect_2_2( np, xp, zp, ex, ey, ez, xmin,     &
!     zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid, ezg,    &
!     ezg_nguard, ezg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)
!     !!! --- Gather magnetic fields on particles
!     CALL pxr_getb2dxz_energy_conserving_vect_2_2( np, xp, zp, bx, by, bz, xmin,     &
!     zmin, dx, dz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid, bzg,    &
!     bzg_nguard, bzg_nvalid, LVEC_fieldgathe, l_lower_order_in_v, l_nodal)

!   ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

!     !!! --- Gather electric and magnetic field on particles
!     CALL pxr_geteb2dxz_energy_conserving_vect_3_3( np, xp, zp, ex, ey, ez, bx, by,  &
!     bz, xmin, zmin, dx, dz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard,           &
!     eyg_nvalid, ezg, ezg_nguard, ezg_nvalid, bxg, bxg_nguard, bxg_nvalid, byg,      &
!     byg_nguard, byg_nvalid, bzg, bzg_nguard, bzg_nvalid, LVEC_fieldgathe,           &
!     l_lower_order_in_v, l_nodal)

!     ! Arbitrary order
!   ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving( np, xp, yp, zp, ex, ey, ez, xmin, zmin,  &
      dx, dz, nox, noz, exg, exg_nguard, exg_nvalid, eyg, eyg_nguard, eyg_nvalid,     &
      ezg, ezg_nguard, ezg_nvalid, l4symtry, .TRUE._idp, l_lower_order_in_v, l_nodal)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_n_energy_conserving( np, xp, yp, zp, bx, by, bz, xmin, zmin,  &
      dx, dz, nox, noz, bxg, bxg_nguard, bxg_nvalid, byg, byg_nguard, byg_nvalid,     &
      bzg, bzg_nguard, bzg_nvalid, l4symtry, .TRUE._idp, l_lower_order_in_v, l_nodal)
!   ENDIF
  ENDIF
END SUBROUTINE

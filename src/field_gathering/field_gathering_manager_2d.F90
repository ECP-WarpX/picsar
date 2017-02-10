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
! ______________________________________________________________________________


! ______________________________________________________________________________
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
!> @param[in] xp,zp particle position arrays
!> @param[inout] ex,ey,ez electric field particle arrays
!> @param[inout] bx,by,bz magnetic field particle arrays
!> @param[in] xmin,zmin tile boundaries
!> @param[in] dx,dz space steps
!> @param[in] nx,nz space discretization
!> @param[in] nxguard, nzguard number of guard cells
!> @param[in] exg, eyg,ezg field arrays
!> @param[in] bxg, byg,bzg field arrays
!> @param[in] l4symetry
!> @param[in] l_lower_order_in_v flag to determine if we interpolate at a lower order
!> @param[in] field_gathe_algo Gathering algorithm
!> @param[in] lvect vector length
!
SUBROUTINE geteb2dxz_energy_conserving(np,xp,yp,zp,ex,ey,ez,bx,by,bz,&
                                       xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                                       nxguard,nyguard,nzguard, &
                                       nox,noy,noz,exg,eyg,ezg,bxg,byg,bzg,&
                                       l4symtry,l_lower_order_in_v,&
                                       lvect, &
                                       field_gathe_algo)
! ______________________________________________________________________________

  USE constants
  USE params
  implicit none

  integer(idp)                  :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
  integer(idp)                  :: field_gathe_algo
  integer(idp)                  :: lvect
  logical(idp), intent(in)      :: l4symtry,l_lower_order_in_v
  real(num), dimension(np)      :: xp,yp,zp,ex,ey,ez,bx,by,bz
  real(num), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
  real(num), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
  real(num)                     :: xmin,ymin,zmin,dx,dy,dz

  IF (field_gathe_algo.lt.0) return

  ! ______________________________________________
  ! Arbitrary order, non-optimized subroutines
  IF (field_gathe_algo.eq.2) THEN


    !!! --- Gather electric field on particles
    CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                          dx,dz,nx,nz,nxguard,nzguard, &
                                          nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
    !!! --- Gather magnetic fields on particles
    CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                          dx,dz,nx,nz,nxguard,nzguard, &
                                          nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)

  ! ______________________________________________
  ! Arbitrary order, scalar subroutines
  ELSE IF (field_gathe_algo.eq.1) THEN

    ! Order 3
    IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_scalar_3_3(np,xp,zp,ex,ey,ez,xmin,zmin,dx,dz,nx,nz, &
                                                     nxguard,nzguard,exg,eyg,ezg,l_lower_order_in_v)

      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_scalar_3_3(np,xp,zp,ex,ey,ez,xmin,zmin,dx,dz,nx,nz, &
                                                     nxguard,nzguard,exg,eyg,ezg,l_lower_order_in_v)

    ! Arbitrary order
    ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                             dx,dz,nx,nz,nxguard,nzguard, &
                                             nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                            dx,dz,nx,nz,nxguard,nzguard, &
                                            nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)

    ENDIF

  ! ________________________________________
  ! Optimized subroutines, default
  ELSE


    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_vect_1_1(np,xp,zp,ex,ey,ez,xmin,zmin,   &
                                            dx,dz,nx,nz,nxguard,nzguard, &
                                            exg,eyg,ezg,LVEC_fieldgathe,l_lower_order_in_v)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_vect_1_1(np,xp,zp,bx,by,bz,xmin,zmin,   &
                                            dx,dz,nx,nz,nxguard,nzguard, &
                                            bxg,byg,bzg,LVEC_fieldgathe,l_lower_order_in_v)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_energy_conserving_vect_2_2(np,xp,zp,ex,ey,ez,xmin,zmin,   &
                                            dx,dz,nx,nz,nxguard,nzguard, &
                                            exg,eyg,ezg,lvect,l_lower_order_in_v)
      !!! --- Gather magnetic fields on particles
      CALL pxr_getb2dxz_energy_conserving_vect_2_2(np,xp,zp,bx,by,bz,xmin,zmin,   &
                                            dx,dz,nx,nz,nxguard,nzguard, &
                                            bxg,byg,bzg,lvect,l_lower_order_in_v)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      !!! --- Gather electric and magnetic field on particles
      CALL pxr_geteb2dxz_energy_conserving_vect_3_3(np,xp,zp,ex,ey,ez,bx,by,bz,xmin,zmin,   &
                                            dx,dz,nx,nz,nxguard,nzguard, &
                                            exg,eyg,ezg,bxg,byg,bzg,lvect, &
                                            l_lower_order_in_v)

    ! Arbitrary order
    ELSE

      !!! --- Gather electric field on particles
      CALL pxr_gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,&
                                             dx,dz,nx,nz,nxguard,nzguard, &
                                             nox,noz,exg,eyg,ezg,l4symtry,.FALSE._idp,l_lower_order_in_v)
      !!! --- Gather magnetic fields on particles
     CALL pxr_getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,&
                                             dx,dz,nx,nz,nxguard,nzguard, &
                                            nox,noz,bxg,byg,bzg,l4symtry,.FALSE._idp,l_lower_order_in_v)

    ENDIF
  ENDIF
END SUBROUTINE geteb2dxz_energy_conserving

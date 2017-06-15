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
! KARKAINNEN.F90
!
! Purpose:
! This file contains subroutines for the Maxwell solver of Karkainnen.
!
! Authors:
! Henri Vincenti
!
! Date:
! Creation 2015
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Karkainnen 2D/3D
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxr_push_em3d_kyeebvec(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx,    &
ny, nz, nxguard, nyguard, nzguard, l_2dxz) 
  USE kyee_em3d
  IMPLICIT NONE
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN) :: dtsdx, dtsdy, dtsdz
  INTEGER(idp) :: j, k, l
  LOGICAL(lp)  :: l_2dxz
  
  IF (.NOT.l_2dxz) THEN
    
    ! advance Bx
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(j, k, l) SHARED(Ex, Ez, Ey, Bx, By, Bz,      &
    !$OMP alphax, alphay, alphaz, betaxy, betaxz, betayx, betayz, betazx, betazy,     &
    !$OMP gammax, gammay, gammaz, dtsdx, dtsdy, dtsdz, nx, ny, nz)  
    !$OMP DO COLLAPSE(3)
    DO l = 0, nz-1
      DO k = 0, ny-1
        DO j = 0, nx
          Bx(j, k, l) = Bx(j, k, l) - alphay*dtsdy * (Ez(j, k+1, l  ) - Ez(j, k, l    &
          )) - betayx*dtsdy * (Ez(j+1, k+1, l  ) - Ez(j+1, k, l  ) +  Ez(j-1, k+1, l  &
          ) - Ez(j-1, k, l  )) - betayz*dtsdy * (Ez(j, k+1, l+1) - Ez(j, k, l+1) +    &
          Ez(j, k+1, l-1) - Ez(j, k, l-1)) - gammay*dtsdy * (Ez(j+1, k+1, l+1) -      &
          Ez(j+1, k, l+1) +  Ez(j-1, k+1, l+1) - Ez(j-1, k, l+1) +  Ez(j+1, k+1, l-1) &
          - Ez(j+1, k, l-1) +  Ez(j-1, k+1, l-1) - Ez(j-1, k, l-1)) + alphaz*dtsdz *  &
          (Ey(j, k, l+1) - Ey(j, k, l  )) + betazx*dtsdz * (Ey(j+1, k, l+1) - Ey(j+1, &
          k, l  ) +  Ey(j-1, k, l+1) - Ey(j-1, k, l  )) + betazy*dtsdz * (Ey(j, k+1,  &
          l+1) - Ey(j, k+1, l  ) +  Ey(j, k-1, l+1) - Ey(j, k-1, l  )) + gammaz*dtsdz &
          * (Ey(j+1, k+1, l+1) - Ey(j+1, k+1, l  ) +  Ey(j-1, k+1, l+1) - Ey(j-1,     &
          k+1, l  ) +  Ey(j+1, k-1, l+1) - Ey(j+1, k-1, l  ) +  Ey(j-1, k-1, l+1) -   &
          Ey(j-1, k-1, l  ))                 
        END DO
      END DO
    END DO
    !$OMP END DO
    
    ! advance By
    !$OMP DO COLLAPSE(3)
    DO l = 0, nz-1
      DO k = 0, ny
        DO j = 0, nx-1
          By(j, k, l) = By(j, k, l) + alphax*dtsdx * (Ez(j+1, k, l  ) - Ez(j, k, l    &
          )) + betaxy*dtsdx * (Ez(j+1, k+1, l  ) - Ez(j, k+1, l  ) +  Ez(j+1, k-1, l  &
          ) - Ez(j, k-1, l  )) + betaxz*dtsdx * (Ez(j+1, k, l+1) - Ez(j, k, l+1) +    &
          Ez(j+1, k, l-1) - Ez(j, k, l-1)) + gammax*dtsdx * (Ez(j+1, k+1, l+1) -      &
          Ez(j, k+1, l+1) +  Ez(j+1, k-1, l+1) - Ez(j, k-1, l+1) +  Ez(j+1, k+1, l-1) &
          - Ez(j, k+1, l-1) +  Ez(j+1, k-1, l-1) - Ez(j, k-1, l-1)) - alphaz*dtsdz *  &
          (Ex(j, k, l+1) - Ex(j, k, l  )) - betazx*dtsdz * (Ex(j+1, k, l+1) - Ex(j+1, &
          k, l  ) +  Ex(j-1, k, l+1) - Ex(j-1, k, l  )) - betazy*dtsdz * (Ex(j, k+1,  &
          l+1) - Ex(j, k+1, l  ) +  Ex(j, k-1, l+1) - Ex(j, k-1, l  )) - gammaz*dtsdz &
          * (Ex(j+1, k+1, l+1) - Ex(j+1, k+1, l  ) +  Ex(j-1, k+1, l+1) - Ex(j-1,     &
          k+1, l  ) +  Ex(j+1, k-1, l+1) - Ex(j+1, k-1, l  ) +  Ex(j-1, k-1, l+1) -   &
          Ex(j-1, k-1, l  ))                 
        END DO
      END DO
    END DO
    !$OMP END DO
    ! advance Bz
    !$OMP DO COLLAPSE(3)
    DO l = 0, nz
      DO k = 0, ny-1
        DO j = 0, nx-1
          Bz(j, k, l) = Bz(j, k, l) - alphax*dtsdx * (Ey(j+1, k, l  ) - Ey(j, k, l    &
          )) - betaxy*dtsdx * (Ey(j+1, k+1, l  ) - Ey(j, k+1, l  ) +  Ey(j+1, k-1, l  &
          ) - Ey(j, k-1, l  )) - betaxz*dtsdx * (Ey(j+1, k, l+1) - Ey(j, k, l+1) +    &
          Ey(j+1, k, l-1) - Ey(j, k, l-1)) - gammax*dtsdx * (Ey(j+1, k+1, l+1) -      &
          Ey(j, k+1, l+1) +  Ey(j+1, k-1, l+1) - Ey(j, k-1, l+1) +  Ey(j+1, k+1, l-1) &
          - Ey(j, k+1, l-1) +  Ey(j+1, k-1, l-1) - Ey(j, k-1, l-1)) + alphay*dtsdy *  &
          (Ex(j, k+1, l  ) - Ex(j, k, l  )) + betayx*dtsdy * (Ex(j+1, k+1, l  ) -     &
          Ex(j+1, k, l  ) +  Ex(j-1, k+1, l  ) - Ex(j-1, k, l  )) + betayz*dtsdy *    &
          (Ex(j, k+1, l+1) - Ex(j, k, l+1) +  Ex(j, k+1, l-1) - Ex(j, k, l-1)) +      &
          gammay*dtsdy * (Ex(j+1, k+1, l+1) - Ex(j+1, k, l+1) +  Ex(j-1, k+1, l+1) -  &
          Ex(j-1, k, l+1) +  Ex(j+1, k+1, l-1) - Ex(j+1, k, l-1) +  Ex(j-1, k+1, l-1) &
          - Ex(j-1, k, l-1))                 
        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  ELSE
    
    k=0
    ! advance Bx
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(j, l) SHARED(k, Ex, Ey, Ez, Bx, By, Bz,      &
    !$OMP alphax, alphaz, betaxz, betazx, dtsdx, dtsdz, nx, nz)  
    !$OMP DO COLLAPSE(2)
    DO l = 0, nz-1
      DO j = 0, nx
        Bx(j, k, l) = Bx(j, k, l) +    alphaz*dtsdz * (Ey(j, k, l+1) - Ey(j, k, l  )) &
        +    betazx*dtsdz * (Ey(j+1, k, l+1) - Ey(j+1, k, l  ) +  Ey(j-1, k, l+1) -   &
        Ey(j-1, k, l  ))  
      END DO
    END DO
    !$OMP END DO
    ! advance By
    !$OMP DO COLLAPSE(2)
    DO l = 0, nz-1
      DO j = 0, nx-1
        By(j, k, l) = By(j, k, l) +    alphax*dtsdx * (Ez(j+1, k, l  ) - Ez(j, k, l   &
        )) +    betaxz*dtsdx * (Ez(j+1, k, l+1) - Ez(j, k, l+1) +  Ez(j+1, k, l-1) -  &
        Ez(j, k, l-1)) -    alphaz*dtsdz * (Ex(j, k, l+1) - Ex(j, k, l  )) -          &
        betazx*dtsdz * (Ex(j+1, k, l+1) - Ex(j+1, k, l  ) +  Ex(j-1, k, l+1) -        &
        Ex(j-1, k, l  ))     
      END DO
    END DO
    !$OMP END DO
    ! advance Bz
    !$OMP DO COLLAPSE(2)
    DO l = 0, nz
      DO j = 0, nx-1
        Bz(j, k, l) = Bz(j, k, l) -    alphax*dtsdx * (Ey(j+1, k, l  ) - Ey(j, k, l   &
        )) -    betaxz*dtsdx * (Ey(j+1, k, l+1) - Ey(j, k, l+1) +  Ey(j+1, k, l-1) -  &
        Ey(j, k, l-1))  
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  END IF
  
  RETURN
END SUBROUTINE pxr_push_em3d_kyeebvec

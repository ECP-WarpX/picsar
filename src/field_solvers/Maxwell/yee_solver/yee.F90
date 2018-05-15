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
! YEE.F90
!
! Purpose:
! This file contains SUBROUTINEs for the FDTD solver of K, Yee
! and the generalization at order n.
!
! Authors:
! Henri Vincenti
! Mathieu Lobet
!
! Date:
! Creation 2015
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee 3D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em3d_evec_norder(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx,  &
  dtsdy, dtsdz, nx, ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs,  &
  nys, nzs, l_nodalgrid)
  USE constants

  INTEGER(idp), INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  INTEGER(idp), INTENT(IN) :: norderx, nordery, norderz
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: Jx, Jy, Jz
  REAL(num), INTENT(IN) :: mudt, dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  INTEGER(idp) :: i, j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid

  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(3)
  ! advance Ex
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
        DO i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Ex(j, k, l) = Ex(j, k, l) + dtsdy(i) * (Bz(j, k+i-ist, l)   - Bz(j, k-i, l  &
          ))
        END DO
        DO i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          Ex(j, k, l) = Ex(j, k, l) - dtsdz(i) * (By(j, k, l+i-ist)   - By(j, k,      &
          l-i))
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
  ! advance Ey
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        Ey(j, k, l) = Ey(j, k, l) - mudt  * Jy(j, k, l)
        DO i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          Ey(j, k, l) = Ey(j, k, l) - dtsdx(i) * (Bz(j+i-ist, k, l)   - Bz(j-i, k,    &
          l))
        END DO
        DO i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          Ey(j, k, l) = Ey(j, k, l) + dtsdz(i) * (Bx(j, k, l+i-ist)   - Bx(j, k,      &
          l-i))
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
  ! advance Ez
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
        DO i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          Ez(j, k, l) = Ez(j, k, l) + dtsdx(i) * (By(j+i-ist, k, l) - By(j-i, k, l))
        END DO
        DO i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Ez(j, k, l) = Ez(j, k, l) - dtsdy(i) * (Bx(j, k+i-ist, l) - Bx(j, k-i, l))
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  RETURN
END SUBROUTINE pxrpush_em3d_evec_norder

! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee 2D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em2d_evec_norder(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx,  &
  dtsdy, dtsdz, nx, ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs,  &
  nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, norderx,      &
  nordery, norderz
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: Jx, Jy, Jz
  REAL(num), INTENT(IN) :: mudt, dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  INTEGER(idp) :: i, j, k, l, ist
  LOGICAL(lp), INTENT(IN)  :: l_nodalgrid

  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF

  k = 0

  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
      DO i = 1, norderz/2
        Ex(j, k, l) = Ex(j, k, l) - dtsdz(i) * (By(j, k, l+i-ist)   - By(j, k, l-i))
      END DO
    END DO
  END DO
  !$OMP END DO

  ! advance Ey
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ey(j, k, l) = Ey(j, k, l) - mudt  * Jy(j, k, l)
      DO i = 1, norderx/2
        Ey(j, k, l) = Ey(j, k, l) - dtsdx(i) * (Bz(j+i-ist, k, l)   - Bz(j-i, k, l))
      END DO
      DO i = 1, norderz/2
        Ey(j, k, l) = Ey(j, k, l) + dtsdz(i) * (Bx(j, k, l+i-ist)   - Bx(j, k, l-i))
      END DO
    END DO
  END DO
  !$OMP END DO

  ! advance Ez
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
      DO i = 1, norderx/2
        Ez(j, k, l) = Ez(j, k, l) + dtsdx(i) * (By(j+i-ist, k, l) - By(j-i, k, l))
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN
END SUBROUTINE pxrpush_em2d_evec_norder


! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em2d_evec(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx, dtsdy,  &
  dtsdz, nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard)    :: Jx, Jy, Jz
  REAL(num) , INTENT(IN)  :: mudt, dtsdx, dtsdy, dtsdz
  INTEGER(idp) :: j, k, l, ist
  LOGICAL(lp) ,INTENT(IN) :: l_nodalgrid
  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF

  k = 0_idp
  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
      Ex(j, k, l) = Ex(j, k, l) - dtsdz * (By(j, k, l+1-ist)   - By(j, k, l-1))
    END DO
  END DO
  !$OMP END DO

  ! advance Ey
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ey(j, k, l) = Ey(j, k, l) - mudt  * Jy(j, k, l)
      Ey(j, k, l) = Ey(j, k, l) - dtsdx * (Bz(j+1-ist, k, l)   - Bz(j-1, k, l))
      Ey(j, k, l) = Ey(j, k, l) + dtsdz * (Bx(j, k, l+1-ist)   - Bx(j, k, l-1))
    END DO
  END DO
  !$OMP END DO

  ! advance Ez
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
      Ez(j, k, l) = Ez(j, k, l) + dtsdx * (By(j+1-ist, k, l) - By(j-1, k, l))
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN
END SUBROUTINE pxrpush_em2d_evec


! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em_pml_2d_evec(ex,ey,ez,bx,by,bz,exy,exz,eyx, eyz,ezx, ezy, bxy, bxz,byx, byz, bzx, bzy,  &
  jx, jy, jz, mudt, dtsdx, dtsdy,  &
  dtsdz, nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz,exy,exz,eyx, eyz,ezx, ezy, bxy, bxz,byx, byz, bzx, bzy
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard)    :: Jx, Jy, Jz
  REAL(num) , INTENT(IN)  :: mudt, dtsdx, dtsdy, dtsdz
  INTEGER(idp) :: j, k, l, ist
  LOGICAL(lp) ,INTENT(IN) :: l_nodalgrid
  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF

  k = 0_idp
  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
      Ex(j, k, l) = Ex(j, k, l) - dtsdz * (By(j, k, l+1-ist)   - By(j, k, l-1))
    END DO
  END DO
  !$OMP END DO

  ! advance Ey
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Eyx(j, k, l) = Eyx(j, k, l) - mudt  * Jy(j, k, l)
      Eyx(j, k, l) = Eyx(j, k, l) - dtsdx * (Bz(j+1-ist, k, l)   - Bz(j-1, k, l))
      Eyz(j, k, l) = Eyz(j, k, l) + dtsdz * (Bx(j, k, l+1-ist)   - Bx(j, k, l-1))
    END DO
  END DO
  !$OMP END DO

  ! advance Ez
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
      Ez(j, k, l) = Ez(j, k, l) + dtsdx * (By(j+1-ist, k, l) - By(j-1, k, l))
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN
END SUBROUTINE pxrpush_em_pml_2d_evec





!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em_upml_2d_evec(dex,dey,dez,hx,hy,hz,ex,ey,ez,bx,by,bz,  &
  jx, jy, jz,dtt, mudt, dtsdx, dtsdy,  &
  dtsdz, nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  USE fields , ONLY : sigma_x_e,sigma_y_e,sigma_z_e
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: dex,dey,dez,ex,ey,ez,bx,by,bz,hx,hy,hz
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard)    :: Jx, Jy, Jz
  REAL(num) , INTENT(IN)  :: mudt, dtsdx, dtsdy, dtsdz,dtt
  INTEGER(idp) :: j, k, l, ist
  LOGICAL(lp) ,INTENT(IN) :: l_nodalgrid
  REAL(num) :: dexold, deyold ,dezold
  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF
  k = 0_idp
  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j,dexold,deyold,dezold)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      dexold = dex(j,k,l) 
      dex(j,k,l) = dex(j,k,l) - dtsdz * (hy(j, k, l+1-ist)   - hy(j, k, l-1)) 



      Ex(j, k, l) = (1.0_num-sigma_z_e(l)*dtt)/(1.0_num+sigma_z_e(l)*dtt)  &
            *Ex(j, k, l) + &
             (1.0_num/(1.0_num+sigma_z_e(l)*dtt)/eps0)* &
             ((1.0_num + sigma_x_e(j)*dtt)*dex(j,k,l) - (1.0_num - sigma_x_e(j)*dtt) * dexold) - mudt*jx(j,k,l) 

    END DO
  END DO
  !$OMP END DO

  ! advance Ey
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      deyold = dey(j,k,l) 

      dey(j,k,l) = (1.0_num/dtt-sigma_z_e(l))/(1.0_num/dtt+sigma_z_e(l)) * dey(j,k,l) + &  
      1.0_num/(1.0_num/dtt+sigma_z_e(l))/dtt * & 
      (dtsdz * (hx(j, k, l+1-ist)   - hx(j, k, l-1))   - dtsdx * (hz(j+1-ist, k, l)   - hz(j-1, k, l))) 

      
      Ey(j,k,l) = (1.0_num - sigma_x_e(j)*dtt)/(1.0_num + sigma_x_e(j)*dtt) * ey(j,k,l)  + &
      (1.0_num/(1.0_num+sigma_x_e(j)*dtt)/eps0) * (dey(j,k,l) -deyold) - mudt  * Jy(j, k, l)
    END DO
  END DO
  !$OMP END DO

  ! advance Ez
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
   
      dezold = dez(j,k,l) 
      dez(j,k,l) = ( 1.0_num/dtt-sigma_x_e(j))/(1.0_num/dtt+sigma_x_e(j))  *  dez(j,k,l) &
      + 1.0_num/(1.0_num/dtt+sigma_x_e(j)) *   dtsdx/dtt * (hy(j+1-ist, k, l) - hy(j-1, k, l))  
       
   
       Ez(j,k,l ) = ez(j,k,l) +1.0_num/eps0*( (1.0_num + sigma_z_e(l)*dtt) *dez(j,k,l)  - (1.0_num - sigma_z_e(l)*dtt) * dezold) -&
       mudt  * Jz(j, k, l)
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  RETURN
END SUBROUTINE pxrpush_em_upml_2d_evec





! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee 3D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em3d_evec(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx, dtsdy,  &
  dtsdz, nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: jx, jy, jz
  REAL(num), INTENT(IN) :: mudt, dtsdx, dtsdy, dtsdz
  INTEGER(idp):: j, k, l
  LOGICAL(lp)  :: l_nodalgrid

  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        Ex(j, k, l) = Ex(j, k, l) + dtsdy * (Bz(j, k, l)   - Bz(j, k-1, l  )) - dtsdz &
        * (By(j, k, l)   - By(j, k, l-1)) - mudt  * jx(j, k, l)
      END DO
    END DO
  END DO
  !$OMP END DO
  ! advance Ey
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        Ey(j, k, l) = Ey(j, k, l) - dtsdx * (Bz(j, k, l)   - Bz(j-1, k, l)) + dtsdz * &
        (Bx(j, k, l)   - Bx(j, k, l-1)) - mudt  * jy(j, k, l)
      END DO
    END DO
  END DO
  !$OMP END DO
  ! advance Ez
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        Ez(j, k, l) = Ez(j, k, l) + dtsdx * (By(j, k, l) - By(j-1, k, l)) - dtsdy *   &
        (Bx(j, k, l) - Bx(j, k-1, l)) - mudt  * jz(j, k, l)
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN
END SUBROUTINE pxrpush_em3d_evec


! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 3D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em3d_bvec_norder(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx,  &
  ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs, nys, nzs,          &
  l_nodalgrid)
  USE constants
  INTEGER(idp), INTENT(IN)    :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs,      &
  norderx, nordery, norderz
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN) :: dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  INTEGER(idp)          :: i, j, k, l, ist
  LOGICAL(lp), INTENT(IN) :: l_nodalgrid

  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF

  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        DO i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Bx(j, k, l) = Bx(j, k, l) - dtsdy(i) * (Ez(j, k+i, l  ) - Ez(j, k-i+ist,    &
          l))
        END DO
        DO i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          Bx(j, k, l) = Bx(j, k, l) + dtsdz(i) * (Ey(j, k, l+i) - Ey(j, k, l-i+ist))
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        DO i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          By(j, k, l) = By(j, k, l) + dtsdx(i) * (Ez(j+i, k, l  ) - Ez(j-i+ist, k,    &
          l))
        END DO
        DO i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          By(j, k, l) = By(j, k, l) - dtsdz(i) * (Ex(j, k, l+i) - Ex(j, k, l-i+ist))
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs
        DO i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          Bz(j, k, l) = Bz(j, k, l) - dtsdx(i) * (Ey(j+i, k, l) - Ey(j-i+ist, k, l))
        END DO
        DO i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Bz(j, k, l) = Bz(j, k, l) + dtsdy(i) * (Ex(j, k+i, l) - Ex(j, k-i+ist, l))
        END DO
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN

END SUBROUTINE pxrpush_em3d_bvec_norder

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 2D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em2d_bvec_norder(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx,  &
  ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs, nys, nzs,          &
  l_nodalgrid)
  USE constants
  INTEGER(idp) , INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, norderx,      &
  nordery, norderz
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN) :: dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  INTEGER(idp) :: i, j, k, l, ist
  LOGICAL(lp) ,INTENT(IN) :: l_nodalgrid

  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF

  k = 0

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  ! advance Bx
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      DO i = 1, norderz/2
        Bx(j, k, l) = Bx(j, k, l) + dtsdz(i) * (Ey(j, k, l+i) - Ey(j, k, l-i+ist))
      END DO
    END DO
  END DO
  !$OMP END DO

  ! advance By
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      DO i = 1, norderx/2
        By(j, k, l) = By(j, k, l) + dtsdx(i) * (Ez(j+i, k, l  ) - Ez(j-i+ist, k, l))
      END DO
      DO i = 1, norderz/2
        By(j, k, l) = By(j, k, l) - dtsdz(i) * (Ex(j, k, l+i) - Ex(j, k, l-i+ist))
      END DO
    END DO
  END DO
  !$OMP END DO

  ! advance Bz
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      DO i = 1, norderx/2
        Bz(j, k, l) = Bz(j, k, l) - dtsdx(i) * (Ey(j+i, k, l) - Ey(j-i+ist, k, l))
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN

END SUBROUTINE pxrpush_em2d_bvec_norder

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em2d_bvec(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx, ny, nz, &
  nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN) :: dtsdx, dtsdy, dtsdz
  INTEGER(idp) :: j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid

  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF
  k = 0_idp
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Bx(j, 0, l) = Bx(j, 0, l) + dtsdz * (Ey(j, 0, l+1) - Ey(j, 0, l-1+ist))
    END DO
  END DO
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      By(j, 0, l) = By(j, 0, l) + dtsdx * (Ez(j+1, 0, l  ) - Ez(j-1+ist, 0, l))
      By(j, 0, l) = By(j, 0, l) - dtsdz * (Ex(j, 0, l+1) - Ex(j, 0, l-1+ist))
    END DO
  END DO
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Bz(j, 0, l) = Bz(j, 0, l) - dtsdx * (Ey(j+1, 0, l) - Ey(j-1+ist, 0, l))
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  RETURN

END SUBROUTINE pxrpush_em2d_bvec

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em_pml_2d_bvec(ex,ey,ez,bx,by,bz,exy,exz,eyx,eyz,ezx,ezy,bxy,bxz,byx,byz,bzx,bzy,    &
  dtsdx, dtsdy, dtsdz, nx, ny, nz, &
  nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz,exy,exz,eyx,eyz,ezx,ezy,bxy,bxz,byx,byz,bzx,bzy
  REAL(num), INTENT(IN) :: dtsdx, dtsdy, dtsdz
  INTEGER(idp) :: j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF
  k = 0_idp
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Bx(j, 0, l) = Bx(j, 0, l) + dtsdz * (Ey(j, 0, l+1) - Ey(j, 0, l-1+ist))
    END DO
  END DO
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Byx(j, 0, l) = Byx(j, 0, l) + dtsdx * (Ez(j+1, 0, l  ) - Ez(j-1+ist, 0, l))
      Byz(j, 0, l) = Byz(j, 0, l) - dtsdz * (Ex(j, 0, l+1) - Ex(j, 0, l-1+ist))
    END DO
  END DO
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      Bz(j, 0, l) = Bz(j, 0, l) - dtsdx * (Ey(j+1, 0, l) - Ey(j-1+ist, 0, l))
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  

  
  RETURN

END SUBROUTINE pxrpush_em_pml_2d_bvec




! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em_upml_2d_bvec(hx,hy,hz,ex,ey,ez,bx,by,bz,dtt,    &
  dtsdx, dtsdy, dtsdz, nx, ny, nz, &
  nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  USE fields , ONLY : sigma_x_b,sigma_y_b,sigma_z_b
  INTEGER(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz,hx,hy,hz
  REAL(num), INTENT(IN) :: dtsdx, dtsdy, dtsdz,dtt
  INTEGER(idp) :: j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  REAL(num) :: bxold,byold, bzold
  IF (l_nodalgrid) THEN
    ist = 0
  ELSE
    ist = 1
  END IF
  k = 0_idp
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j,bxold,byold,bzold)
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      bxold = bx(j,k,l)
      Bx(j, 0, l) = Bx(j, 0, l) + dtsdz * (Ey(j, 0, l+1) - Ey(j, 0, l-1+ist))
      hx(j, 0, l) = (1.0_num -sigma_z_b(l)*dtt)/(1.0_num + sigma_z_b(l)*dtt) * hx(j,0,l) + &
         (1.0_num / ((1.0_num+sigma_z_b(l)*dtt)*mu0))* &
         ((1.0_num + sigma_x_b(j)*dtt) * bx(j,k,l) -( 1.0_num - sigma_x_b(j)*dtt)*bxold)
   ! hx(j,k,l ) =exp(-sigma_z_b(l)*dtt) * hx(j,0,l) + &
   ! imu0*((1.0_num + sigma_x_b(j)*dtt) * bx(j,k,l) -( 1.0_num - sigma_x_b(j)*dtt)*bxold)
    END DO
  END DO
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      byold = by(j,0,l) 
      
     by(j,0,l) = (1.0_num/dtt - sigma_z_b(l))/(1.0_num/dtt+sigma_z_b(l))*by(j,0,l) - 1.0_num/(1.0_num/dtt + sigma_z_b(l))/dtt * &
     (-dtsdx * (Ez(j+1, 0, l  ) - Ez(j-1+ist, 0, l)) +dtsdz * (Ex(j, 0, l+1) - Ex(j, 0, l-1+ist))) 
 !        by(j,k,l ) = exp(-sigma_z_b(l) *dtt) * by(j,k,l) +dtsdx * (Ez(j+1, 0, l  ) &
 !        - Ez(j-1+ist, 0, l)) - dtsdz * (Ex(j, 0, l+1) - Ex(j, 0, l-1+ist))

 !      hy(j,0,l ) = exp(-sigma_x_b(j) *dtt) *hy(j,k,l) +(1.0_num / mu0)*(by(j,0,l) - byold)
      hy(j,0,l) = (1.0_num - sigma_x_b(j)*dtt)/(1.0_num+sigma_x_b(j)*dtt) * hy(j,0,l) + &
     (1.0_num / mu0)*(by(j,0,l) - byold)

    END DO
  END DO
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(2)
  DO l = -nzs, nz+nzs
    DO j = -nxs, nx+nxs
      bzold = bz(j,0,l)
      bz(j,0,l) = (1.0_num/dtt - sigma_x_b(j))/(1.0_num/dtt+sigma_x_b(j))*bz(j,0,l) -  &
      (1.0_num/(1.0_num/dtt+sigma_x_b(j)))/dtt *dtsdx * (Ey(j+1, 0, l) - Ey(j-1+ist, 0, l)) 
! bz(j,k,l )= exp(-sigma_x_b(j) *dtt) * bz(j,k,l) - dtsdx * (Ey(j+1, 0, l) - Ey(j-1+ist, 0, l))

      hz(j,0,l) = hz(j,0,l) + 1.0_num/mu0*((1.0_num + sigma_z_b(l)* dtt )*bz(j,0,l) - (1.0_num - sigma_z_b(l)* dtt)*bzold)

    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  
  RETURN

END SUBROUTINE pxrpush_em_upml_2d_bvec








! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 3D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em3d_bvec(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx, ny, nz, &
  nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)
  USE constants
  INTEGER(idp):: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  REAL(num), INTENT(IN OUT), DIMENSION(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  REAL(num), INTENT(IN) :: dtsdx, dtsdy, dtsdz
  INTEGER(idp) :: j, k, l
  LOGICAL(lp)  :: l_nodalgrid

  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs-1
    DO k = -nys, ny+nys-1
      DO j = -nxs, nx+nxs
        Bx(j, k, l) = Bx(j, k, l) - dtsdy * (Ez(j, k+1, l  ) - Ez(j, k, l)) + dtsdz * &
        (Ey(j, k, l+1) - Ey(j, k, l))
      END DO
    END DO
  END DO
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs-1
    DO k = -nys, ny+nys
      DO j = -nxs, nx+nxs-1
        By(j, k, l) = By(j, k, l) + dtsdx * (Ez(j+1, k, l  ) - Ez(j, k, l)) - dtsdz * &
        (Ex(j, k, l+1) - Ex(j, k, l))
      END DO
    END DO
  END DO
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(3)
  DO l = -nzs, nz+nzs
    DO k = -nys, ny+nys-1
      DO j = -nxs, nx+nxs-1
        Bz(j, k, l) = Bz(j, k, l) - dtsdx * (Ey(j+1, k, l) - Ey(j, k, l)) + dtsdy *   &
        (Ex(j, k+1, l) - Ex(j, k, l))
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
END SUBROUTINE pxrpush_em3d_bvec

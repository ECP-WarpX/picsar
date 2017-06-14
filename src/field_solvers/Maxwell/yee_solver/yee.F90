! ______________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
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
! YEE.F90
!
! Purpose:
! This file contains subroutines for the FDTD solver of K, Yee
! and the generalization at order n.
!
! Authors:
! Henri Vincenti
! Mathieu Lobet
!
! Date:
! Creation 2015
! ______________________________________________________________________________


! ______________________________________________________________________________
!> @brief
!> Push electric field Yee 3D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em3d_evec_norder(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx,  &
dtsdy, dtsdz, nx, ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs,  &
nys, nzs, l_nodalgrid)    
  ! ______________________________________________________________________________
  
  use constants
  
  integer(idp), intent(IN) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  integer(idp), intent(IN) :: norderx, nordery, norderz
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: Jx, Jy, Jz
  real(num), intent(IN) :: mudt, dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  integer(idp) :: i, j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  
  if (l_nodalgrid) then
    ist = 0
  else
    ist = 1
  end if
  
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(3)
  ! advance Ex
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
        do i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Ex(j, k, l) = Ex(j, k, l) + dtsdy(i) * (Bz(j, k+i-ist, l)   - Bz(j, k-i, l  &
          ))
        end do
        do i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          Ex(j, k, l) = Ex(j, k, l) - dtsdz(i) * (By(j, k, l+i-ist)   - By(j, k,      &
          l-i))
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
  ! advance Ey
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        Ey(j, k, l) = Ey(j, k, l) - mudt  * Jy(j, k, l)
        do i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          Ey(j, k, l) = Ey(j, k, l) - dtsdx(i) * (Bz(j+i-ist, k, l)   - Bz(j-i, k,    &
          l))
        end do
        do i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          Ey(j, k, l) = Ey(j, k, l) + dtsdz(i) * (Bx(j, k, l+i-ist)   - Bx(j, k,      &
          l-i))
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
  ! advance Ez
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
        do i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          Ez(j, k, l) = Ez(j, k, l) + dtsdx(i) * (By(j+i-ist, k, l) - By(j-i, k, l))
        end do
        do i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Ez(j, k, l) = Ez(j, k, l) - dtsdy(i) * (Bx(j, k+i-ist, l) - Bx(j, k-i, l))
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  return
end subroutine pxrpush_em3d_evec_norder

! ______________________________________________________________________________
!> @brief
!> Push electric field Yee 2D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em2d_evec_norder(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx,  &
dtsdy, dtsdz, nx, ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs,  &
nys, nzs, l_nodalgrid)    
  ! ______________________________________________________________________________
  
  use constants
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, norderx,      &
  nordery, norderz
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: Jx, Jy, Jz
  real(num), intent(IN) :: mudt, dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  integer(idp) :: i, j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  
  if (l_nodalgrid) then
    ist = 0
  else
    ist = 1
  end if
  
  k = 0
  
  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
      do i = 1, norderz/2
        Ex(j, k, l) = Ex(j, k, l) - dtsdz(i) * (By(j, k, l+i-ist)   - By(j, k, l-i))
      end do
    end do
  end do
  !$OMP END DO
  
  ! advance Ey
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ey(j, k, l) = Ey(j, k, l) - mudt  * Jy(j, k, l)
      do i = 1, norderx/2
        Ey(j, k, l) = Ey(j, k, l) - dtsdx(i) * (Bz(j+i-ist, k, l)   - Bz(j-i, k, l))
      end do
      do i = 1, norderz/2
        Ey(j, k, l) = Ey(j, k, l) + dtsdz(i) * (Bx(j, k, l+i-ist)   - Bx(j, k, l-i))
      end do
    end do
  end do
  !$OMP END DO
  ! advance Ez
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
      do i = 1, norderx/2
        Ez(j, k, l) = Ez(j, k, l) + dtsdx(i) * (By(j+i-ist, k, l) - By(j-i, k, l))
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine pxrpush_em2d_evec_norder


! ______________________________________________________________________________
!> @brief
!> Push electric field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em2d_evec(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx, dtsdy,  &
dtsdz, nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)   
  ! ______________________________________________________________________________
  
  use constants
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: Jx, Jy, Jz
  real(num), intent(IN) :: mudt, dtsdx, dtsdy, dtsdz
  integer(idp) :: j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  
  if (l_nodalgrid) then
    ist = 0
  else
    ist = 1
  end if
  
  k = 0
  
  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ex(j, k, l) = Ex(j, k, l) - mudt  * Jx(j, k, l)
      Ex(j, k, l) = Ex(j, k, l) - dtsdz * (By(j, k, l+1-ist)   - By(j, k, l-1))
    end do
  end do
  !$OMP END DO
  
  ! advance Ey
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ey(j, k, l) = Ey(j, k, l) - mudt  * Jy(j, k, l)
      Ey(j, k, l) = Ey(j, k, l) - dtsdx * (Bz(j+1-ist, k, l)   - Bz(j-1, k, l))
      Ey(j, k, l) = Ey(j, k, l) + dtsdz * (Bx(j, k, l+1-ist)   - Bx(j, k, l-1))
    end do
  end do
  !$OMP END DO
  
  ! advance Ez
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ez(j, k, l) = Ez(j, k, l) - mudt  * Jz(j, k, l)
      Ez(j, k, l) = Ez(j, k, l) + dtsdx * (By(j+1-ist, k, l) - By(j-1, k, l))
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine pxrpush_em2d_evec


! ______________________________________________________________________________
!> @brief
!> Push electric field Yee 3D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em3d_evec(ex, ey, ez, bx, by, bz, jx, jy, jz, mudt, dtsdx, dtsdy,  &
dtsdz, nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)   
  ! ______________________________________________________________________________
  
  use constants
  
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,          &
  -nzguard:nz+nzguard) :: jx, jy, jz
  real(num), intent(IN) :: mudt, dtsdx, dtsdy, dtsdz
  integer(idp):: j, k, l
  LOGICAL(lp)  :: l_nodalgrid
  
  ! advance Ex
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        Ex(j, k, l) = Ex(j, k, l) + dtsdy * (Bz(j, k, l)   - Bz(j, k-1, l  )) - dtsdz &
        * (By(j, k, l)   - By(j, k, l-1)) - mudt  * jx(j, k, l)  
      end do
    end do
  end do
  !$OMP END DO
  ! advance Ey
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        Ey(j, k, l) = Ey(j, k, l) - dtsdx * (Bz(j, k, l)   - Bz(j-1, k, l)) + dtsdz * &
        (Bx(j, k, l)   - Bx(j, k, l-1)) - mudt  * jy(j, k, l)  
      end do
    end do
  end do
  !$OMP END DO
  ! advance Ez
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        Ez(j, k, l) = Ez(j, k, l) + dtsdx * (By(j, k, l) - By(j-1, k, l)) - dtsdy *   &
        (Bx(j, k, l) - Bx(j, k-1, l)) - mudt  * jz(j, k, l)  
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
end subroutine pxrpush_em3d_evec


! ______________________________________________________________________________
!> @brief
!> Push magnetic field Yee 3D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em3d_bvec_norder(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx,  &
ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs, nys, nzs,          &
l_nodalgrid)    
  ! ______________________________________________________________________________
  
  use constants
  
  integer(idp)          :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs,      &
  norderx, nordery, norderz
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN) :: dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  integer(idp)          :: i, j, k, l, ist
  LOGICAL(lp)                :: l_nodalgrid
  
  if (l_nodalgrid) then
    ist = 0
  else
    ist = 1
  end if
  
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        do i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Bx(j, k, l) = Bx(j, k, l) - dtsdy(i) * (Ez(j, k+i, l  ) - Ez(j, k-i+ist,    &
          l))
        end do
        do i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          Bx(j, k, l) = Bx(j, k, l) + dtsdz(i) * (Ey(j, k, l+i) - Ey(j, k, l-i+ist))
        end do
      end do
    end do
  end do
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        do i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          By(j, k, l) = By(j, k, l) + dtsdx(i) * (Ez(j+i, k, l  ) - Ez(j-i+ist, k,    &
          l))
        end do
        do i = 1, MIN(MIN(norderz/2, (nz-l)+nzguard), l+nzguard)
          By(j, k, l) = By(j, k, l) - dtsdz(i) * (Ex(j, k, l+i) - Ex(j, k, l-i+ist))
        end do
      end do
    end do
  end do
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs
        do i = 1, MIN(MIN(norderx/2, (nx-j)+nxguard), j+nxguard)
          Bz(j, k, l) = Bz(j, k, l) - dtsdx(i) * (Ey(j+i, k, l) - Ey(j-i+ist, k, l))
        end do
        do i = 1, MIN(MIN(nordery/2, (ny-k)+nyguard), k+nyguard)
          Bz(j, k, l) = Bz(j, k, l) + dtsdy(i) * (Ex(j, k+i, l) - Ex(j, k-i+ist, l))
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
  
end subroutine pxrpush_em3d_bvec_norder

! ______________________________________________________________________________
!> @brief
!> Push magnetic field Yee 2D arbitrary order
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em2d_bvec_norder(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx,  &
ny, nz, norderx, nordery, norderz, nxguard, nyguard, nzguard, nxs, nys, nzs,          &
l_nodalgrid)    
  ! ______________________________________________________________________________
  
  use constants
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs, norderx,      &
  nordery, norderz
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN) :: dtsdx(norderx/2), dtsdy(nordery/2), dtsdz(norderz/2)
  integer(idp) :: i, j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  
  if (l_nodalgrid) then
    ist = 0
  else
    ist = 1
  end if
  
  k = 0
  
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j, i)
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      do i = 1, norderz/2
        Bx(j, k, l) = Bx(j, k, l) + dtsdz(i) * (Ey(j, k, l+i) - Ey(j, k, l-i+ist))
      end do
    end do
  end do
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      do i = 1, norderx/2
        By(j, k, l) = By(j, k, l) + dtsdx(i) * (Ez(j+i, k, l  ) - Ez(j-i+ist, k, l))
      end do
      do i = 1, norderz/2
        By(j, k, l) = By(j, k, l) - dtsdz(i) * (Ex(j, k, l+i) - Ex(j, k, l-i+ist))
      end do
    end do
  end do
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      do i = 1, norderx/2
        Bz(j, k, l) = Bz(j, k, l) - dtsdx(i) * (Ey(j+i, k, l) - Ey(j-i+ist, k, l))
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
  
end subroutine pxrpush_em2d_bvec_norder

! ______________________________________________________________________________
!> @brief
!> Push magnetic field Yee 2D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em2d_bvec(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx, ny, nz, &
nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)   
  ! ______________________________________________________________________________
  use constants
  integer(idp) :: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz
  real(num), intent(IN) :: dtsdx, dtsdy, dtsdz
  integer(idp) :: j, k, l, ist
  LOGICAL(lp)  :: l_nodalgrid
  
  if (l_nodalgrid) then
    ist = 0
  else
    ist = 1
  end if
  
  k = 0_idp
  
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Bx(j, 0, l) = Bx(j, 0, l) + dtsdz * (Ey(j, 0, l+1) - Ey(j, 0, l-1+ist))
    end do
  end do
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      By(j, 0, l) = By(j, 0, l) + dtsdx * (Ez(j+1, 0, l  ) - Ez(j-1+ist, 0, l))
      By(j, 0, l) = By(j, 0, l) - dtsdz * (Ex(j, 0, l+1) - Ex(j, 0, l-1+ist))
    end do
  end do
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(2)
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Bz(j, 0, l) = Bz(j, 0, l) - dtsdx * (Ey(j+1, 0, l) - Ey(j-1+ist, 0, l))
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  return
  
end subroutine pxrpush_em2d_bvec

! ______________________________________________________________________________
!> @brief
!> Push magnetic field Yee 3D order 2
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
subroutine pxrpush_em3d_bvec(ex, ey, ez, bx, by, bz, dtsdx, dtsdy, dtsdz, nx, ny, nz, &
nxguard, nyguard, nzguard, nxs, nys, nzs, l_nodalgrid)   
  ! ______________________________________________________________________________
  use constants
  integer(idp):: nx, ny, nz, nxguard, nyguard, nzguard, nxs, nys, nzs
  real(num), intent(IN OUT), dimension(-nxguard:nx+nxguard, -nyguard:ny+nyguard,      &
  -nzguard:nz+nzguard) :: ex, ey, ez, bx, by, bz 
  real(num), intent(IN) :: dtsdx, dtsdy, dtsdz
  integer(idp) :: j, k, l
  LOGICAL(lp)  :: l_nodalgrid
  
  ! advance Bx
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l, k, j)
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs-1
    do k = -nys, ny+nys-1
      do j = -nxs, nx+nxs
        Bx(j, k, l) = Bx(j, k, l) - dtsdy * (Ez(j, k+1, l  ) - Ez(j, k, l)) + dtsdz * &
        (Ey(j, k, l+1) - Ey(j, k, l)) 
      end do
    end do
  end do
  !$OMP END DO
  ! advance By
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs-1
    do k = -nys, ny+nys
      do j = -nxs, nx+nxs-1
        By(j, k, l) = By(j, k, l) + dtsdx * (Ez(j+1, k, l  ) - Ez(j, k, l)) - dtsdz * &
        (Ex(j, k, l+1) - Ex(j, k, l)) 
      end do
    end do
  end do
  !$OMP END DO
  ! advance Bz
  !$OMP DO COLLAPSE(3)
  do l = -nzs, nz+nzs
    do k = -nys, ny+nys-1
      do j = -nxs, nx+nxs-1
        Bz(j, k, l) = Bz(j, k, l) - dtsdx * (Ey(j+1, k, l) - Ey(j, k, l)) + dtsdy *   &
        (Ex(j, k+1, l) - Ex(j, k, l)) 
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine pxrpush_em3d_bvec

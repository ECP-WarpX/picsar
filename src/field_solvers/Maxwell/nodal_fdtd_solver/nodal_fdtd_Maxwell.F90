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
! nodal_fdtd_Maxwell.F90
!
! Purpose:
! This file contains SUBROUTINEs for the FDTD solver of K, Yee
! and the generalization at order n.
!
! Authors:
! Jean-Luc Vay
!
! Date:
! Creation 2018
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Push electric field nodal 3D order n
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Jean-Luc Vay
!> Weiqun Zhang
!> Maxence Thevenet
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________
subroutine pxrpush_em3d_enodal( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex, exlo, exhi,&
     ey, eylo, eyhi, &
     ez, ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     jx, jxlo, jxhi, &
     jy, jylo, jyhi, &
     jz, jzlo, jzhi, &
     norderx, nordery, norderz, &
     mudt, dtsdx, dtsdy, dtsdz)
USE picsar_precision, ONLY: idp, num, isp
! ______________________________________________________________________________

#ifdef WARPX
  integer(isp), intent(in) :: &
       xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),&
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3),&
       jxlo(3),jxhi(3),jylo(3),jyhi(3),jzlo(3),jzhi(3)
#else
  integer(idp), intent(in) :: &
       xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),&
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3),&
       jxlo(3),jxhi(3),jylo(3),jyhi(3),jzlo(3),jzhi(3)
#endif
  real(num), intent(IN OUT):: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(num), intent(IN OUT):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(num), intent(IN OUT):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))

  real(num), intent(IN):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(num), intent(IN):: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(num), intent(IN):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))

  real(num), intent(IN):: jx(jxlo(1):jxhi(1),jxlo(2):jxhi(2),jxlo(3):jxhi(3))
  real(num), intent(IN):: jy(jylo(1):jyhi(1),jylo(2):jyhi(2),jylo(3):jyhi(3))
  real(num), intent(IN):: jz(jzlo(1):jzhi(1),jzlo(2):jzhi(2),jzlo(3):jzhi(3))

  real(num), intent(IN) :: mudt,dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
  
  integer, intent(IN) :: norderx, nordery, norderz

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j, i), &
  !$OMP SHARED(norderx, nordery, norderz), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, mudt, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz, jx, jy, jz)
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ex,Bz,By,jx)
!$acc loop gang vector collapse(3)
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
        Ex(j,k,l) = Ex(j,k,l) - mudt * Jx(j,k,l)
        do i = 1, MIN(MIN(nordery/2, exhi(2)-k), k+xlo(2)-exlo(2))
          Ex(j,k,l) = Ex(j,k,l) + dtsdy(i) * (Bz(j,k+i,l) - Bz(j,k-i,l))
        end do
        do i = 1, MIN(MIN(norderz/2, exhi(3)-l), l+xlo(3)-exlo(3))
          Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i) - By(j,k,l-i))
        end do
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ey,Bz,Bx,jy)
!$acc loop gang vector collapse(3)
  do l     = ylo(3), yhi(3)
    do k   = ylo(2), yhi(2)
      do j = ylo(1), yhi(1)
        Ey(j,k,l) = Ey(j,k,l) - mudt * Jy(j,k,l)
        do i = 1, MIN(MIN(norderx/2, eyhi(1)-j), j+ylo(1)-eylo(1))
          Ey(j,k,l) = Ey(j,k,l) - dtsdx(i) * (Bz(j+i,k,l) - Bz(j-i,k,l))
        end do
        do i = 1, MIN(MIN(norderz/2, eyhi(3)-l), l+ylo(3)-eylo(3))
          Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i) - Bx(j,k,l-i))
        end do
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ez,By,Bx,jz)
!$acc loop gang vector collapse(3)
  do l     = zlo(3), zhi(3)
    do k   = zlo(2), zhi(2)
      do j = zlo(1), zhi(1)
        Ez(j,k,l) = Ez(j,k,l) - mudt * Jz(j,k,l)
        do i = 1, MIN(MIN(norderx/2, exhi(1)-j), j+xlo(1)-exlo(1))
          Ez(j,k,l) = Ez(j,k,l) + dtsdx(i) * (By(j+i,k,l) - By(j-i,k,l))
        end do
        do i = 1, MIN(MIN(nordery/2, exhi(2)-k), k+xlo(2)-exlo(2))
          Ez(j,k,l) = Ez(j,k,l) - dtsdy(i) * (Bx(j,k+i,l) - Bx(j,k-i,l))
        end do
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

end subroutine pxrpush_em3d_enodal

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field nodal 3D order n
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Jean-Luc Vay
!> Weiqun Zhang
!> Maxence Thevenet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine pxrpush_em3d_bnodal( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     norderx, nordery, norderz, &
     dtsdx, dtsdy, dtsdz)
USE picsar_precision, ONLY: idp, num, isp
! ______________________________________________________________________________


#ifdef WARPX
  integer(isp) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),&
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3)
#else
  integer(idp) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),&
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3)
#endif
  real(num), intent(IN):: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(num), intent(IN):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(num), intent(IN):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))

  real(num), intent(INOUT):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(num), intent(INOUT):: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(num), intent(INOUT):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))

  integer, intent(IN) :: norderx, nordery, norderz

  real(num), intent(IN) :: dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j, i), &
  !$OMP SHARED(norderx, nordery, norderz), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz)
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Bx,Ez,Ey)
!$acc loop gang vector collapse(3)
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
        do i = 1, MIN(MIN(nordery/2, exhi(2)-k), k+xlo(2)-exlo(2))
          Bx(j,k,l) = Bx(j,k,l) - dtsdy(i) * (Ez(j,k+i,l) - Ez(j,k-i,l))
        end do
        do i = 1, MIN(MIN(norderz/2, exhi(3)-l), l+xlo(3)-exlo(3))
          Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,l+i) - Ey(j,k,l-i))
        end do
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(By,Ez,Ex)
!$acc loop gang vector collapse(3)
  do l     = ylo(3), yhi(3)
    do k   = ylo(2), yhi(2)
      do j = ylo(1), yhi(1)
        do i = 1, MIN(MIN(norderx/2, exhi(1)-j), j+xlo(1)-exlo(1))
          By(j,k,l) = By(j,k,l) + dtsdx(i) * (Ez(j+i,k,l) - Ez(j-i,k,l))
        end do
        do i = 1, MIN(MIN(norderz/2, exhi(3)-l), l+xlo(3)-exlo(3))
          By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j,k,l+i) - Ex(j,k,l-i))
        end do
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Bz,Ey,Ex)
!$acc loop gang vector collapse(3)
  do l     = zlo(3), zhi(3)
    do k   = zlo(2), zhi(2)
      do j = zlo(1), zhi(1)
        do i = 1, MIN(MIN(norderx/2, exhi(1)-j), j+xlo(1)-exlo(1))
          Bz(j,k,l) = Bz(j,k,l) - dtsdx(i) * (Ey(j+i,k,l) - Ey(j-i,k,l))
        end do
        do i = 1, MIN(MIN(nordery/2, exhi(2)-k), k+xlo(2)-exlo(2))
          Bz(j,k,l) = Bz(j,k,l) + dtsdy(i) * (Ex(j,k+i,l) - Ex(j,k-i,l))
        end do
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

end subroutine pxrpush_em3d_bnodal


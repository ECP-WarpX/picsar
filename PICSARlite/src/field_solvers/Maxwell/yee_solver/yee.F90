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
!> Push electric field Yee 3D order 2
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine pxrpush_em3d_evec( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     jx, jxlo, jxhi, &
     jy, jylo, jyhi, &
     jz, jzlo, jzhi, &
     mudt, dtsdx,dtsdy,dtsdz)
! ______________________________________________________________________________

  use constants

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

  real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, mudt, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz, jx, jy, jz)
  !$OMP DO COLLAPSE(3)
#endif
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
        Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l) - Bz(j,k-1,l  )) &
                              - dtsdz * (By(j,k,l) - By(j,k  ,l-1)) &
                              - mudt  * jx(j,k,l)
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
  do l     = ylo(3), yhi(3)
    do k   = ylo(2), yhi(2)
      do j = ylo(1), yhi(1)
        Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                              + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                              - mudt  * jy(j,k,l)
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
  do l     = zlo(3), zhi(3)
    do k   = zlo(2), zhi(2)
      do j = zlo(1), zhi(1)
        Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                              - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                              - mudt  * jz(j,k,l)
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

end subroutine pxrpush_em3d_evec

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 3D order 2
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine pxrpush_em3d_bvec( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     dtsdx,dtsdy,dtsdz)
! ______________________________________________________________________________

  USE constants

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

  real(num), intent(IN) :: dtsdx,dtsdy,dtsdz

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz)
  !$OMP DO COLLAPSE(3)
#endif
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
           Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j  ,k+1,l  ) - Ez(j,k,l)) &
                                 + dtsdz * (Ey(j  ,k  ,l+1) - Ey(j,k,l))
       end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
  do l     = ylo(3), yhi(3)
    do k   = ylo(2), yhi(2)
      do j = ylo(1), yhi(1)
        By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k  ,l  ) - Ez(j,k,l)) &
                              - dtsdz * (Ex(j  ,k  ,l+1) - Ex(j,k,l))
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
  do l     = zlo(3), zhi(3)
    do k   = zlo(2), zhi(2)
      do j = zlo(1), zhi(1)
           Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k  ,l  ) - Ey(j,k,l)) &
                                 + dtsdy * (Ex(j  ,k+1,l  ) - Ex(j,k,l))
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

end subroutine pxrpush_em3d_bvec

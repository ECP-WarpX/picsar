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
! karkkainen.F90
!
! Purpose:
! This file contains SUBROUTINEs for the FDTD solver of Karkkainen.
!
! Authors:
! Henri Vincenti
! Mathieu Lobet
! Weiqun Zhang
! Jean-Luc Vay
! Maxence Thevenet
!
! Date:
! Creation 2015
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the magnetic field with the 3D Yee FDTD
!> scheme (order 2).
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!>
!> Also, this routine allows to specify the indices between which
!> the fields should be updated. These indices are not necessarily the same
!> as the full bounds of the field arrays. For instance, if this routine is
!> called by different OpenMP threads, each thread will handle a different part
!> of the array, and thus the corresponding indices that are passed as argument
!> to this routine will be different.
!>
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!>
!> @param[in] xlo the lowest indices (in 3D) at which to update the Bx field
!> @param[in] xhi the highest indices (in 3D) at which to update the Bx field
!> @param[in] ylo the lowest indices (in 3D) at which to update the By field
!> @param[in] yhi the highest indices (in 3D) at which to update the By field
!> @param[in] zlo the lowest indices (in 3D) at which to update the Bz field
!> @param[in] zhi the highest indices (in 3D) at which to update the Bz field
!> @param[inout] bx the array of values of the Bx field
!> @param[in] bxlo the lowest bound (in 3D) of the array `bx` 
!> (`bxlo` is always lower than `xlo`)
!> @param[in] bxhi the highest bound (in 3D) of the array `bx` 
!> (`bxlo` is always higher than `xlo`)
!> @param[inout] by the array of values of the by field
!> @param[in] bylo the lowest bound (in 3D) of the array `by` 
!> (`bylo` is always lower than `ylo`)
!> @param[in] byhi the highest bound (in 3D) of the array `by`
!>  (`bylo` is always higher than `ylo`)
!> @param[inout] bz the array of values of the bz field
!> @param[in] bzlo the lowest bound (in 3D) of the array `bz` 
!> (`bzlo` is always lower than `zlo`)
!> @param[in] bzhi the highest bound (in 3D) of the array `bz` 
!> (`bzlo` is always higher than `zlo`)
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine pxrpush_em3d_bvec_ckc( &
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

  integer :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3), &
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3)

  real(num), intent(IN):: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(num), intent(IN):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(num), intent(IN):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))

  real(num), intent(INOUT):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(num), intent(INOUT):: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(num), intent(INOUT):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))

  real(num), intent(IN) :: dtsdx,dtsdy,dtsdz

  integer :: j,k,l

  real(num) :: delta, rx, ry, rz, betaxz, betaxy, betayx, betayz, betazx, betazy 
  real(num) :: beta, alphax, alphay, alphaz, gammax, gammay, gammaz
  
  ! CKC push
  ! computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
  delta = max(dtsdx,dtsdy,dtsdz)
  rx = (dtsdx/delta)**2
  ry = (dtsdy/delta)**2
  rz = (dtsdz/delta)**2
  beta = 0.125*(1.-rx*ry*rz/(ry*rz+rz*rx+rx*ry))
  betaxy = ry*beta
  betaxz = rz*beta
  betayx = rx*beta
  betayz = rz*beta
  betazx = rx*beta
  betazy = ry*beta
  gammax = ry*rz*(1./16.-0.125*ry*rz/(ry*rz+rz*rx+rx*ry))
  gammay = rx*rz*(1./16.-0.125*rx*rz/(ry*rz+rz*rx+rx*ry))
  gammaz = rx*ry*(1./16.-0.125*rx*ry/(ry*rz+rz*rx+rx*ry))
  alphax = 1. - 2.*betaxy - 2.* betaxz - 4.*gammax
  alphay = 1. - 2.*betayx - 2.* betayz - 4.*gammay
  alphaz = 1. - 2.*betazx - 2.* betazy - 4.*gammaz

  betaxy = dtsdx*betaxy
  betaxz = dtsdx*betaxz
  betayx = dtsdy*betayx
  betayz = dtsdy*betayz
  betazx = dtsdz*betazx
  betazy = dtsdz*betazy
  alphax = dtsdx*alphax
  alphay = dtsdy*alphay
  alphaz = dtsdz*alphaz
  gammax = dtsdx*gammax
  gammay = dtsdy*gammay
  gammaz = dtsdz*gammaz

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx), &
  !$OMP SHARED(alphax, alphay, alphaz, gammax, gammay, gammaz), &
  !$OMP SHARED(betaxy, betaxz, betayx, betayz, betazx, betazy), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz)
  !$OMP DO COLLAPSE(3)
#endif
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
        Bx(j,k,l) = Bx(j,k,l) - alphay * (Ez(j  ,k+1,l  ) - Ez(j,  k  ,l  )) &
                              - betayx * (Ez(j+1,k+1,l  ) - Ez(j+1,k  ,l  )  &
                                       +  Ez(j-1,k+1,l  ) - Ez(j-1,k  ,l  )) &
                              - betayz * (Ez(j  ,k+1,l+1) - Ez(j  ,k  ,l+1)  &
                                       +  Ez(j  ,k+1,l-1) - Ez(j  ,k  ,l-1)) &
                              - gammay * (Ez(j+1,k+1,l+1) - Ez(j+1,k  ,l+1)  &
                                       +  Ez(j-1,k+1,l+1) - Ez(j-1,k  ,l+1)  &
                                       +  Ez(j+1,k+1,l-1) - Ez(j+1,k  ,l-1)  &
                                       +  Ez(j-1,k+1,l-1) - Ez(j-1,k  ,l-1)) &
                              + alphaz * (Ey(j  ,k  ,l+1) - Ey(j,  k,  l  )) &
                              + betazx * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  )  &
                                       +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  )) &
                              + betazy * (Ey(j  ,k+1,l+1) - Ey(j  ,k+1,l  )  &
                                       +  Ey(j  ,k-1,l+1) - Ey(j  ,k-1,l  )) &
                              + gammaz * (Ey(j+1,k+1,l+1) - Ey(j+1,k+1,l  )  &
                                       +  Ey(j-1,k+1,l+1) - Ey(j-1,k+1,l  )  &
                                       +  Ey(j+1,k-1,l+1) - Ey(j+1,k-1,l  )  &
                                       +  Ey(j-1,k-1,l+1) - Ey(j-1,k-1,l  ))
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
        By(j,k,l) = By(j,k,l) + alphax * (Ez(j+1,k  ,l  ) - Ez(j,  k,  l  )) &
                              + betaxy * (Ez(j+1,k+1,l  ) - Ez(j  ,k+1,l  )  &
                                       +  Ez(j+1,k-1,l  ) - Ez(j  ,k-1,l  )) &
                              + betaxz * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1)  &
                                       +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                              + gammax * (Ez(j+1,k+1,l+1) - Ez(j  ,k+1,l+1)  &
                                       +  Ez(j+1,k-1,l+1) - Ez(j  ,k-1,l+1)  &
                                       +  Ez(j+1,k+1,l-1) - Ez(j  ,k+1,l-1)  &
                                       +  Ez(j+1,k-1,l-1) - Ez(j  ,k-1,l-1)) &
                              - alphaz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                              - betazx * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  )  &
                                       +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  )) &
                              - betazy * (Ex(j  ,k+1,l+1) - Ex(j  ,k+1,l  )  &
                                       +  Ex(j  ,k-1,l+1) - Ex(j  ,k-1,l  )) &
                              - gammaz * (Ex(j+1,k+1,l+1) - Ex(j+1,k+1,l  )  &
                                       +  Ex(j-1,k+1,l+1) - Ex(j-1,k+1,l  )  &
                                       +  Ex(j+1,k-1,l+1) - Ex(j+1,k-1,l  )  &
                                       +  Ex(j-1,k-1,l+1) - Ex(j-1,k-1,l  ))
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
  do l       = zlo(3), zhi(3)
    do k    = zlo(2), zhi(2)
      do j = zlo(1), zhi(1)
        Bz(j,k,l) = Bz(j,k,l) - alphax * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  )) &
                              - betaxy * (Ey(j+1,k+1,l  ) - Ey(j  ,k+1,l  )  &
                                       +  Ey(j+1,k-1,l  ) - Ey(j  ,k-1,l  )) &
                              - betaxz * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1)  &
                                       +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) &
                              - gammax * (Ey(j+1,k+1,l+1) - Ey(j  ,k+1,l+1)  &
                                       +  Ey(j+1,k-1,l+1) - Ey(j  ,k-1,l+1)  &
                                       +  Ey(j+1,k+1,l-1) - Ey(j  ,k+1,l-1)  &
                                       +  Ey(j+1,k-1,l-1) - Ey(j  ,k-1,l-1)) &
                              + alphay * (Ex(j  ,k+1,l  ) - Ex(j  ,k  ,l  )) &
                              + betayx * (Ex(j+1,k+1,l  ) - Ex(j+1,k  ,l  )  &
                                       +  Ex(j-1,k+1,l  ) - Ex(j-1,k  ,l  )) &
                              + betayz * (Ex(j  ,k+1,l+1) - Ex(j  ,k  ,l+1)  &
                                       +  Ex(j  ,k+1,l-1) - Ex(j  ,k  ,l-1)) &
                              + gammay * (Ex(j+1,k+1,l+1) - Ex(j+1,k  ,l+1)  &
                                       +  Ex(j-1,k+1,l+1) - Ex(j-1,k  ,l+1)  &
                                       +  Ex(j+1,k+1,l-1) - Ex(j+1,k  ,l-1)  &
                                       +  Ex(j-1,k+1,l-1) - Ex(j-1,k  ,l-1))
      end do
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

end subroutine pxrpush_em3d_bvec_ckc


! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the magnetic field with the 2D CKC FDTD
!> scheme (order 2).
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Weiqun Zhang
!> Mathieu Lobet
!> Jean-Luc Vay
!> Maxence Thevenet
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine pxrpush_em2d_bvec_ckc( &
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
  integer(isp) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2),&
       bxlo(2),bxhi(2),bylo(2),byhi(2),bzlo(2),bzhi(2)
#else
  integer(idp) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2),&
       bxlo(2),bxhi(2),bylo(2),byhi(2),bzlo(2),bzhi(2)
#endif
  real(num), intent(IN):: ex(exlo(1):exhi(1),exlo(2):exhi(2))
  real(num), intent(IN):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2))
  real(num), intent(IN):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2))

  real(num), intent(INOUT):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
  real(num), intent(INOUT):: by(bylo(1):byhi(1),bylo(2):byhi(2))
  real(num), intent(INOUT):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2))

  real(num), intent(IN) :: dtsdx,dtsdy,dtsdz

  integer :: j,k
  
  real(num) :: delta, rx, rz, betaxz, betazx, alphax, alphaz

  ! dtsdy should not be used.
  ! Cole-Karkkainen-Cowan push
  ! computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
  delta = max(dtsdx,dtsdz)
  rx = (dtsdx/delta)**2
  rz = (dtsdz/delta)**2
  betaxz = 0.125*rz
  betazx = 0.125*rx
  alphax = 1. - 2.*betaxz
  alphaz = 1. - 2.*betazx

  betaxz = dtsdx*betaxz
  betazx = dtsdz*betazx
  alphax = dtsdx*alphax
  alphaz = dtsdz*alphaz

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, alphax, alphaz), &
  !$OMP SHARED(betaxz, betazx, ex, ey, ez, bx, by, bz)
  !$OMP DO COLLAPSE(2)
#endif
  do k    = xlo(2), xhi(2)
    do j  = xlo(1), xhi(1)
      Bx(j,k) = Bx(j,k) + alphaz * (Ey(j  ,k+1) - Ey(j,  k)) &
                        + betazx * (Ey(j+1,k+1) - Ey(j+1,k)  &
                                 +  Ey(j-1,K+1) - Ey(j-1,k))
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
  do k    = ylo(2), yhi(2)
    do j = ylo(1), yhi(1)
      By(j,k) = By(j,k) + alphax * (Ez(j+1,k  ) - Ez(j,k))   &
                        + betaxz * (Ez(j+1,k+1) - Ez(j,k+1)  &
                                 +  Ez(j+1,k-1) - Ez(j,k-1)) &
                        - alphaz * (Ex(j  ,k+1) - Ex(j,k  )) &
                        - betazx * (Ex(j+1,k+1) - Ex(j+1,k)  &
                                 +  Ex(j-1,k+1) - Ex(j-1,k))
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
  do k    = zlo(2), zhi(2)
    do j = zlo(1), zhi(1)
      Bz(j,k) = Bz(j,k) - alphax * (Ey(j+1,k  ) - Ey(j,k))  &
                        - betaxz * (Ey(j+1,k+1) - Ey(j,k+1) &
                                 +  Ey(j+1,k-1) - Ey(j,k-1))
    end do
  end do
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif
end subroutine pxrpush_em2d_bvec_ckc

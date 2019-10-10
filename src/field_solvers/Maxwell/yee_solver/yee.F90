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
  USE picsar_precision, ONLY: idp, lp, num

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
  USE picsar_precision, ONLY: idp, lp, num
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
!> This SUBROUTINE is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet

!> @date
!> Creation 2015
! ________________________________________________________________________________________
SUBROUTINE pxrpush_em2d_evec( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex, exlo, exhi, &
     ey, eylo, eyhi, &
     ez, ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     jx, jxlo, jxhi, &
     jy, jylo, jyhi, &
     jz, jzlo, jzhi, &
     mudt, dtsdx, dtsdy, dtsdz)
     USE picsar_precision, ONLY: idp, isp, num


#ifdef WARPX
  integer(isp), intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2),&
       bxlo(2),bxhi(2),bylo(2),byhi(2),bzlo(2),bzhi(2),&
       jxlo(2),jxhi(2),jylo(2),jyhi(2),jzlo(2),jzhi(2)
#else
  integer(idp), intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2),&
       bxlo(2),bxhi(2),bylo(2),byhi(2),bzlo(2),bzhi(2),&
       jxlo(2),jxhi(2),jylo(2),jyhi(2),jzlo(2),jzhi(2)
#endif
  real(num), intent(IN OUT):: ex(exlo(1):exhi(1),exlo(2):exhi(2))
  real(num), intent(IN OUT):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2))
  real(num), intent(IN OUT):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2))

  real(num), intent(IN):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
  real(num), intent(IN):: by(bylo(1):byhi(1),bylo(2):byhi(2))
  real(num), intent(IN):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2))

  real(num), intent(IN):: jx(jxlo(1):jxhi(1),jxlo(2):jxhi(2))
  real(num), intent(IN):: jy(jylo(1):jyhi(1),jylo(2):jyhi(2))
  real(num), intent(IN):: jz(jzlo(1):jzhi(1),jzlo(2):jzhi(2))

  real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz

  integer :: j,k

  ! dtsdy should not be used.

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, mudt, dtsdx, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz, jx, jy, jz)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(ex,By,jx)
!$acc loop gang vector collapse(2)
  do k   = xlo(2), xhi(2)
    do j = xlo(1), xhi(1)
      ex(j,k) = ex(j,k) - dtsdz * (By(j,k) - By(j,k-1)) &
                        - mudt  * jx(j,k)
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Ey,Bz,Bx,jy)
!$acc loop gang vector collapse(2)
  do k   = ylo(2), yhi(2)
    do j = ylo(1), yhi(1)
      Ey(j,k) = Ey(j,k) - dtsdx * (Bz(j,k) - Bz(j-1,k)) &
                        + dtsdz * (Bx(j,k) - Bx(j,k-1)) &
                        - mudt  * jy(j,k)
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Ez,By,jz)
!$acc loop gang vector collapse(2)
  do k   = zlo(2), zhi(2)
    do j = zlo(1), zhi(1)
      Ez(j,k) = Ez(j,k) + dtsdx * (By(j,k) - By(j-1,k  )) &
                        - mudt  * jz(j,k)
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

END SUBROUTINE pxrpush_em2d_evec

! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee RZ order 2
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!> David Grote

!> @date
!> Creation 2019
! ________________________________________________________________________________________
subroutine pxrpush_emrz_evec( &
     rlo, rhi, tlo, thi, zlo, zhi, &
     Er, erlo, erhi, &
     Et, etlo, ethi, &
     Ez, ezlo, ezhi, &
     Br, brlo, brhi, &
     Bt, btlo, bthi, &
     Bz, bzlo, bzhi, &
     Jr, jrlo, jrhi, &
     Jt, jtlo, jthi, &
     Jz, jzlo, jzhi, &
     mudt, dtsdx, dtsdy, dtsdz, rmin, dr)
     USE picsar_precision, ONLY: idp, isp, num


#ifdef WARPX
  integer(isp), intent(in) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2),&
       jrlo(2),jrhi(2),jtlo(2),jthi(2),jzlo(2),jzhi(2)
#else
  integer(idp), intent(in) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2),&
       jrlo(2),jrhi(2),jtlo(2),jthi(2),jzlo(2),jzhi(2)
#endif
  real(num), intent(IN OUT):: Er(erlo(1):erhi(1),erlo(2):erhi(2))
  real(num), intent(IN OUT):: Et(etlo(1):ethi(1),etlo(2):ethi(2))
  real(num), intent(IN OUT):: Ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2))

  real(num), intent(IN):: Br(brlo(1):brhi(1),brlo(2):brhi(2))
  real(num), intent(IN):: Bt(btlo(1):bthi(1),btlo(2):bthi(2))
  real(num), intent(IN):: Bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2))

  real(num), intent(IN):: Jr(jrlo(1):jrhi(1),jrlo(2):jrhi(2))
  real(num), intent(IN):: Jt(jtlo(1):jthi(1),jtlo(2):jthi(2))
  real(num), intent(IN):: Jz(jzlo(1):jzhi(1),jzlo(2):jzhi(2))

  real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz,rmin,dr

  integer :: j,k

  real(num) :: ru, rd

  ! dtsdy should not be used.

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j, ru, rd), &
  !$OMP SHARED(rlo, rhi, tlo, thi, zlo, zhi, mudt, dtsdx, dtsdz, rmin, dr), &
  !$OMP SHARED(Er, Et, Ez, Br, Bt, Bz, Jr, Jt, Jz)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Er,Bt,Jr)
!$acc loop gang vector collapse(2)
  do k   = rlo(2), rhi(2)
    do j = rlo(1), rhi(1)
      Er(j,k) = Er(j,k) - dtsdz * (Bt(j,k) - Bt(j,k-1)) &
                        - mudt  * Jr(j,k)
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Et,Bz,Br,Jt)
!$acc loop gang vector collapse(2)
  do k   = tlo(2), thi(2)
    do j = tlo(1), thi(1)
      if (j /= 0 .or. rmin /= 0.) then
        Et(j,k) = Et(j,k) - dtsdx * (Bz(j,k) - Bz(j-1,k)) &
                          + dtsdz * (Br(j,k) - Br(j,k-1)) &
                          - mudt  * Jt(j,k)
      endif
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Ez,Bt,Jz)
!$acc loop gang vector collapse(2)
  do k   = zlo(2), zhi(2)
    do j = zlo(1), zhi(1)
      if (j /= 0 .or. rmin /= 0.) then
        ru = 1. + 0.5/(rmin/dr + j)
        rd = 1. - 0.5/(rmin/dr + j)
        Ez(j,k) = Ez(j,k) + dtsdx * (ru*Bt(j,k) - rd*Bt(j-1,k)) &
                          - mudt  * Jz(j,k)
      else
        Ez(j,k) = Ez(j,k) + 4.*dtsdx * Bt(j,k) &
                          - mudt  * Jz(j,k)
      end if
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif
end subroutine pxrpush_emrz_evec

! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee RZ multimode order 2
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Jean-Luc Vay
!> Henri Vincenti
!> Remi Lehe
!> David Grote

!> @date
!> Creation 2019
! ________________________________________________________________________________________
subroutine pxrpush_emrz_evec_multimode( &
     rlo, rhi, tlo, thi, zlo, zhi, &
     nmodes, &
     Er, erlo, erhi, &
     Et, etlo, ethi, &
     Ez, ezlo, ezhi, &
     Br, brlo, brhi, &
     Bt, btlo, bthi, &
     Bz, bzlo, bzhi, &
     Jr, jrlo, jrhi, &
     Jt, jtlo, jthi, &
     Jz, jzlo, jzhi, &
     mudt, dtsdx, dtsdy, dtsdz, rmin, dr)
     USE picsar_precision, ONLY: idp, isp, num

  integer(idp) :: nmodes
#ifdef WARPX
  integer(isp), intent(in) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2),&
       jrlo(2),jrhi(2),jtlo(2),jthi(2),jzlo(2),jzhi(2)
#else
  integer(idp), intent(in) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2),&
       jrlo(2),jrhi(2),jtlo(2),jthi(2),jzlo(2),jzhi(2)
#endif
  complex(num), intent(IN OUT):: Er(erlo(1):erhi(1),erlo(2):erhi(2),0:nmodes-1)
  complex(num), intent(IN OUT):: Et(etlo(1):ethi(1),etlo(2):ethi(2),0:nmodes-1)
  complex(num), intent(IN OUT):: Ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),0:nmodes-1)

  complex(num), intent(IN):: Br(brlo(1):brhi(1),brlo(2):brhi(2),0:nmodes-1)
  complex(num), intent(IN):: Bt(btlo(1):bthi(1),btlo(2):bthi(2),0:nmodes-1)
  complex(num), intent(IN):: Bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),0:nmodes-1)

  complex(num), intent(IN):: Jr(jrlo(1):jrhi(1),jrlo(2):jrhi(2),0:nmodes-1)
  complex(num), intent(IN):: Jt(jtlo(1):jthi(1),jtlo(2):jthi(2),0:nmodes-1)
  complex(num), intent(IN):: Jz(jzlo(1):jzhi(1),jzlo(2):jzhi(2),0:nmodes-1)

  real(num), intent(IN) :: mudt, dtsdx, dtsdy, dtsdz, rmin, dr

  integer(idp) :: j, k, m
  real(kind=8) :: w, r, rd, ru, dt
  complex(kind=8) :: i=(0., 1.)

  ! ===============================
  !             2-D RZ multipole
  ! ===============================

  dt = dtsdx*dr

  do m = 0, nmodes-1

     ! advance Er
#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j, ru, rd, r), &
  !$OMP SHARED(rlo, rhi, tlo, thi, zlo, zhi, mudt, dtsdx, dtsdz, rmin, dr), &
  !$OMP SHARED(Er, Et, Ez, Br, Bt, Bz, Jr, Jt, Jz)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Er,Bt,Bz,Jr)
!$acc loop gang vector collapse(2)
     do k   = rlo(2), rhi(2)
       do j = rlo(1), rhi(1)
           r = rmin+j*dr+0.5*dr
           Er(j,k,m) = Er(j,k,m) - i*m*dt*Bz(j,k,m)/r &
                - dtsdz * (Bt(j,k,m)   - Bt(j  ,k-1,m)) &
                - mudt  * Jr(j,k,m)
        end do
     end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
#endif

     ! advance Etheta
#ifndef WARPX
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Et,Bz,Br,Jt)
!$acc loop gang vector collapse(2)
     do k   = tlo(2), thi(2)
       do j = tlo(1), thi(1)
         if (j /= 0) then
           ! Equation used in the bulk of the grid
           Et(j,k,m) = Et(j,k,m) - dtsdx * (Bz(j,k,m) - Bz(j-1,k,m)) &
                + dtsdz * (Br(j,k,m) - Br(j,k-1,m)) &
                - mudt * Jt(j,k,m)
         endif
       end do
       if (tlo(1) <= 0 .and. 0 <= thi(1)) then
          j = 0
          if (rmin /= 0.) then
             Et(j,k,m) = Et(j,k,m) - dtsdx * (Bz(j,k,m) - Bz(j-1,k,m)) &
                               + dtsdz * (Br(j,k,m) - Br(j,k-1,m)) &
                               - mudt  * Jt(j,k,m)
          else
             if (m == 1) then
                ! The bulk equation could in principle be used here since it does not diverge
                ! on axis. However, it typically gives poore results e.g. for the propagation
                ! of a laser pulse (The field is spuriously reduced on axis.) For this reason
                ! a modified on-axis condition is used here : we use the fact that
                ! Etheta(r=0,m=1) should equal -iEr(r=0,m=1), for the fields Er and Et to be
                ! independent of theta at r=0. Now with linear interpolation :
                ! Er(r=0,m=1) = 0.5*[Er(r=dr/2,m=1)+Er(r=-dr/2,m=1)]
                ! And using the rule applying for the guards cells (see em3d_applybc_e)
                ! Er(r=-dr/2,m=1) = Er(r=dr/2,m=1). Thus :
                Et(j,k,m) = -i*Er(j,k,m)
             else
                ! Etheta should remain 0 on axis, for modes different than m=1
                Et(j,k,m) = 0
             endif
          endif
       end if
     end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
#endif

     ! advance Ez
#ifndef WARPX
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Ez,Br,Bt,Jz)
!$acc loop gang vector collapse(2)
     do k   = zlo(2), zhi(2)
       do j = zlo(1), zhi(1)
         if (j /= 0) then
           ru = 1. + 0.5/(rmin/dr + j)
           rd = 1. - 0.5/(rmin/dr + j)
           r = rmin + j*dr
           Ez(j,k,m) = Ez(j,k,m) + dtsdx * (ru*Bt(j,k,m) - rd*Bt(j-1  ,k,m)) &
                + i*m*dt*Br(j,k,m)/r &
                - mudt  * Jz(j,k,m)
         end if
       end do
       if (zlo(1) <= 0 .and. 0 <= zhi(1)) then
         j = 0
         if (rmin == 0.) then
           if (m == 0) then
              Ez(j,k,m) = Ez(j,k,m) + 4.*dtsdx * Bt(j,k,m) &
                             - mudt  * Jz(j,k,m)
           else
              ! Ez should remain 0 on axis, for modes with m>0,
              ! but the bulk equation does not necessarily ensure this.
              Ez(j,k,m) = 0.
           endif
         else
           ru = 1. + 0.5/(rmin/dr)
           rd = 1. - 0.5/(rmin/dr)
           r = rmin + j*dr
           Ez(j,k,m) = Ez(j,k,m) + dtsdx * (ru*Bt(j,k,m) - rd*Bt(j-1  ,k,m)) &
                + i*m*dt*Br(j,k,m)/r &
                - mudt  * Jz(j,k,m)
         end if
       end if
     end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

  end do ! nmodes
return
end subroutine pxrpush_emrz_evec_multimode

! ________________________________________________________________________________________
!> @brief
!> Push electric field Yee 3D order 2
!> This SUBROUTINE is general enough to be called by AMReX.
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
SUBROUTINE pxrpush_em3d_evec( &
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
USE picsar_precision, ONLY: idp, isp, num
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

  real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, mudt, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz, jx, jy, jz)
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ex,Bz,By,jx)
!$acc loop gang vector collapse(3)
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
        Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l) - Bz(j,k-1,l  )) &
                              - dtsdz * (By(j,k,l) - By(j,k  ,l-1)) &
                              - mudt  * jx(j,k,l)
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
        Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                              + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                              - mudt  * jy(j,k,l)
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
        Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                              - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                              - mudt  * jz(j,k,l)
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

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
  USE picsar_precision, ONLY: idp, lp, num
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
  USE picsar_precision, ONLY: idp, lp, num
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
!> This SUBROUTINE pushes the magnetic field with the 2D Yee FDTD
!> scheme (order 2).
!> This SUBROUTINE is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!> regions.
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
SUBROUTINE pxrpush_em2d_bvec( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     dtsdx,dtsdy,dtsdz)
USE picsar_precision, ONLY: idp, isp, num
! ______________________________________________________________________________


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

  ! dtsdy should not be used.

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Bx,Ey)
!$acc loop gang vector collapse(2)
  do k   = xlo(2), xhi(2)
    do j = xlo(1), xhi(1)
        Bx(j,k) = Bx(j,k) + dtsdz * (Ey(j  ,k+1) - Ey(j,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(By,Ez,Ex)
!$acc loop gang vector collapse(2)
  do k   = ylo(2), yhi(2)
    do j = ylo(1), yhi(1)
        By(j,k) = By(j,k) + dtsdx * (Ez(j+1,k  ) - Ez(j,k)) &
                          - dtsdz * (Ex(j  ,k+1) - Ex(j,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Bz,Ey)
!$acc loop gang vector collapse(2)
  do k   = zlo(2), zhi(2)
    do j = zlo(1), zhi(1)
      Bz(j,k) = Bz(j,k) - dtsdx * (Ey(j+1,k  ) - Ey(j,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif
END SUBROUTINE pxrpush_em2d_bvec

! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the magnetic field with the RZ Yee FDTD
!> scheme (order 2).
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!> regions.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!> David Grote
!
!> @date
!> Creation 2019
! ________________________________________________________________________________________
subroutine pxrpush_emrz_bvec( &
     rlo, rhi, tlo, thi, zlo, zhi, &
     Er, erlo, erhi,&
     Et, etlo, ethi, &
     Ez, ezlo, ezhi, &
     Br, brlo, brhi, &
     Bt, btlo, bthi, &
     Bz, bzlo, bzhi, &
     dtsdr,dtsdt,dtsdz,rmin,dr)
USE picsar_precision, ONLY: idp, isp, num
! ______________________________________________________________________________


#ifdef WARPX
  integer(isp) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2)
#else
  integer(idp) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2)
#endif
  real(num), intent(IN):: Er(erlo(1):erhi(1),erlo(2):erhi(2))
  real(num), intent(IN):: Et(etlo(1):ethi(1),etlo(2):ethi(2))
  real(num), intent(IN):: Ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2))

  real(num), intent(INOUT):: Br(brlo(1):brhi(1),brlo(2):brhi(2))
  real(num), intent(INOUT):: Bt(btlo(1):bthi(1),btlo(2):bthi(2))
  real(num), intent(INOUT):: Bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2))

  real(num), intent(IN) :: dtsdr,dtsdt,dtsdz,rmin,dr

  integer :: j,k

  real(num) :: ru, rd

  ! dtsdt should not be used.

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j, ru, rd), &
  !$OMP SHARED(rlo, rhi, tlo, thi, zlo, zhi, dtsdr, dtsdz, rmin, dr), &
  !$OMP SHARED(Er, Et, Ez, Br, Bt, Bz)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Br,Et)
!$acc loop gang vector collapse(2)
  do k   = rlo(2), rhi(2)
    do j = rlo(1), rhi(1)
        Br(j,k) = Br(j,k) + dtsdz * (Et(j  ,k+1) - Et(j,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Bt,Ez,Er)
!$acc loop gang vector collapse(2)
  do k   = tlo(2), thi(2)
    do j = tlo(1), thi(1)
        Bt(j,k) = Bt(j,k) + dtsdr * (Ez(j+1,k  ) - Ez(j,k)) &
                          - dtsdz * (Er(j  ,k+1) - Er(j,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Bz,Et)
!$acc loop gang vector collapse(2)
  do k   = zlo(2), zhi(2)
    do j = zlo(1), zhi(1)
      ru = 1. + 0.5/(rmin/dr + j + 0.5)
      rd = 1. - 0.5/(rmin/dr + j + 0.5)
      Bz(j,k) = Bz(j,k) - dtsdr * (ru*Et(j+1,k  ) - rd*Et(j,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif
end subroutine pxrpush_emrz_bvec

! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the magnetic field with the RZ Yee FDTD multimode
!> scheme (order 2).
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!> regions.
!
!> @author
!> Jean-Luc Vay
!> Henri Vincenti
!> Remi Lehe
!> David Grote
!
!> @date
!> Creation 2019
! ________________________________________________________________________________________
subroutine pxrpush_emrz_bvec_multimode( &
     rlo, rhi, tlo, thi, zlo, zhi, &
     nmodes, &
     Er, erlo, erhi, &
     Et, etlo, ethi, &
     Ez, ezlo, ezhi, &
     Br, brlo, brhi, &
     Bt, btlo, bthi, &
     Bz, bzlo, bzhi, &
     dtsdx, dtsdy, dtsdz, rmin, dr)
     USE picsar_precision, ONLY: idp, isp, num

  integer(idp) :: nmodes
#ifdef WARPX
  integer(isp), intent(in) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2)
#else
  integer(idp), intent(in) :: rlo(2), rhi(2), tlo(2), thi(2), zlo(2), zhi(2), &
       erlo(2),erhi(2),etlo(2),ethi(2),ezlo(2),ezhi(2),&
       brlo(2),brhi(2),btlo(2),bthi(2),bzlo(2),bzhi(2)
#endif
  complex(num), intent(IN):: Er(erlo(1):erhi(1),erlo(2):erhi(2),0:nmodes-1)
  complex(num), intent(IN):: Et(etlo(1):ethi(1),etlo(2):ethi(2),0:nmodes-1)
  complex(num), intent(IN):: Ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),0:nmodes-1)

  complex(num), intent(IN OUT):: Br(brlo(1):brhi(1),brlo(2):brhi(2),0:nmodes-1)
  complex(num), intent(IN OUT):: Bt(btlo(1):bthi(1),btlo(2):bthi(2),0:nmodes-1)
  complex(num), intent(IN OUT):: Bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),0:nmodes-1)

  real(num), intent(IN) :: dtsdx, dtsdy, dtsdz, rmin, dr

  integer(idp) :: j, k, m
  real(kind=8) :: w, r, rd, ru, dt
  complex(kind=8) :: i=(0., 1.)

  dt = dtsdx*dr
  do m = 0, nmodes-1

     ! advance Br
#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j, ru, rd), &
  !$OMP SHARED(rlo, rhi, tlo, thi, zlo, zhi, dtsdr, dtsdz, rmin, dr), &
  !$OMP SHARED(Er, Et, Ez, Br, Bt, Bz)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Br,Et,Ez)
!$acc loop gang vector collapse(2)
     do k   = rlo(2), rhi(2)
       do j = rlo(1), rhi(1)
           if (j==0 .and. rmin==0) then
              ! On axis
              if (.not. m == 1) then
                 ! Br should remain 0 on axis, for modes different than m=1,
                 ! but the bulk equation does not necessarily ensure this.
                 Br(j,k,m) = 0.
              else
                 ! For the mode m = 1, the bulk equation diverges on axis
                 ! (due to the 1/r terms). The following expressions regularize
                 ! these divergences by assuming, on axis :
                 ! Ez/r = 0/r + dEz/dr
                 Br(j,k,m) = Br(j,k,m) + i*m*dt*Ez(j+1,k,m)/dr &
                      + dtsdz * (Et(j,  k+1,m) - Et(j,k,m))
              endif
           else
              ! Equations in the bulk of the grid
              r = rmin + j*dr
              Br(j,k,m) = Br(j,k,m) + i*m*dt*Ez(j,k,m)/r &
                   + dtsdz * (Et(j,  k+1,m) - Et(j,k,m))
           endif
        end do
     end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
#endif

     ! advance Btheta
#ifndef WARPX
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Bt,Ez,Er)
!$acc loop gang vector collapse(2)
     do k   = tlo(2), thi(2)
       do j = tlo(1), thi(1)
           Bt(j,k,m) = Bt(j,k,m) + dtsdx * (Ez(j+1,k  ,m) - Ez(j,k,m)) &
                - dtsdz * (Er(j,k+1,m) - Er(j,k,m))
        end do
     end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
#endif

     ! advance Bz
#ifndef WARPX
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Bz,Er,Et)
!$acc loop gang vector collapse(2)
     do k   = zlo(2), zhi(2)
       do j = zlo(1), zhi(1)
           r  = rmin + j*dr + 0.5*dr
           ru = 1. + 0.5/(rmin/dr + j + 0.5)
           rd = 1. - 0.5/(rmin/dr + j + 0.5)
           Bz(j,k,m) = Bz(j,k,m) - dtsdx * (ru*Et(j+1,k,m) - rd*Et(j,k,m)) &
                - i*m*dt*Er(j,k,m)/r
        end do
     end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

  end do ! nmodes
return
end subroutine pxrpush_emrz_bvec_multimode

! ________________________________________________________________________________________
!> @brief
!> Push magnetic field Yee 3D order 2
!> This SUBROUTINE is general enough to be called by AMReX.
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
SUBROUTINE pxrpush_em3d_bvec( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     dtsdx,dtsdy,dtsdz)
USE picsar_precision, ONLY: idp, isp, num
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

  real(num), intent(IN) :: dtsdx,dtsdy,dtsdz

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, bx, by, bz)
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Bx,Ez,Ey)
!$acc loop gang vector collapse(3)
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
           Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j  ,k+1,l  ) - Ez(j,k,l)) &
                                 + dtsdz * (Ey(j  ,k  ,l+1) - Ey(j,k,l))
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
        By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k  ,l  ) - Ez(j,k,l)) &
                              - dtsdz * (Ex(j  ,k  ,l+1) - Ex(j,k,l))
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
           Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k  ,l  ) - Ey(j,k,l)) &
                                 + dtsdy * (Ex(j  ,k+1,l  ) - Ex(j,k,l))
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

END SUBROUTINE pxrpush_em3d_bvec



! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the electric field with charge-conserving term,
!> using the 2D Yee FDTD scheme (order 2).
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!> regions.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!> Remi Lehe
!
!> @date
!> Creation 2018
! ________________________________________________________________________________________
subroutine pxrpush_em2d_evec_f( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     f, flo, fhi, &
     dtsdx, dtsdy, dtsdz)
USE picsar_precision, ONLY: idp, isp, num
! ______________________________________________________________________________


#ifdef WARPX
  integer(isp) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2), flo(2), fhi(2)
#else
  integer(idp) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2), flo(2), fhi(2)
#endif
  real(num), intent(INOUT):: ex(exlo(1):exhi(1),exlo(2):exhi(2))
  real(num), intent(INOUT):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2))
  real(num), intent(INOUT):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2))

  real(num), intent(IN):: f(flo(1):fhi(1),flo(2):fhi(2))

  real(num), intent(IN) :: dtsdx, dtsdy, dtsdz

  integer :: j,k

  ! dtsdy should not be used.

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx, dtsdz), &
  !$OMP SHARED(ex, ey, ez, f)
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Ex,F)
!$acc loop gang vector collapse(2)
  do k   = xlo(2), xhi(2)
    do j = xlo(1), xhi(1)
        Ex(j,k) = Ex(j,k) + dtsdx * (F(j+1,k) - F(j  ,k))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(2)
#endif
!$acc parallel deviceptr(Ez,F)
!$acc loop gang vector collapse(2)
  do k   = zlo(2), zhi(2)
    do j = zlo(1), zhi(1)
      Ez(j,k) = Ez(j,k) + dtsdz * (F(j,k+1) - F(j,k  ))
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif
end subroutine pxrpush_em2d_evec_f

! ________________________________________________________________________________________
!> @brief
!> This subroutine pushes the electric field with charge-conserving term,
!> using the 3D Yee FDTD scheme (order 2).
!> This subroutine is general enough to be called by AMReX.
!> OMP pragmas are ignored when compiled for WarpX.
!
!> @author
!> Henri Vincenti
!> Mathieu Lobet
!> Weiqun Zhang
!> Jean-Luc Vay
!> Maxence Thevenet
!> Remi Lehe
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
subroutine pxrpush_em3d_evec_f( &
     xlo, xhi, ylo, yhi, zlo, zhi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     f, flo, fhi, &
     dtsdx,dtsdy,dtsdz)
USE picsar_precision, ONLY: idp, isp, num
! ______________________________________________________________________________


#ifdef WARPX
  integer(isp) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),flo(3),fhi(3)
#else
  integer(idp) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),flo(3),fhi(3)
#endif
  real(num), intent(INOUT):: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(num), intent(INOUT):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(num), intent(INOUT):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))

  real(num), intent(IN):: f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

  real(num), intent(IN) :: dtsdx,dtsdy,dtsdz

  integer :: j,k,l

#ifndef WARPX
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l, k, j), &
  !$OMP SHARED(xlo, xhi, ylo, yhi, zlo, zhi, dtsdx, dtsdy, dtsdz), &
  !$OMP SHARED(ex, ey, ez, f)
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ex,F)
!$acc loop gang vector collapse(3)
  do l     = xlo(3), xhi(3)
    do k   = xlo(2), xhi(2)
      do j = xlo(1), xhi(1)
         Ex(j,k,l) = Ex(j,k,l) + dtsdx * (F(j+1,k  ,l  ) - F(j,k,l))
       end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ey,F)
!$acc loop gang vector collapse(3)
  do l     = ylo(3), yhi(3)
    do k   = ylo(2), yhi(2)
      do j = ylo(1), yhi(1)
        Ey(j,k,l) = Ey(j,k,l) + dtsdy * (F(j  ,k+1,l  ) - F(j,k,l))
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP DO COLLAPSE(3)
#endif
!$acc parallel deviceptr(Ez,F)
!$acc loop gang vector collapse(3)
  do l     = zlo(3), zhi(3)
    do k   = zlo(2), zhi(2)
      do j = zlo(1), zhi(1)
        Ez(j,k,l) = Ez(j,k,l) + dtsdz * (F(j  ,k  ,l+1) - F(j,k,l))
      end do
    end do
  end do
!$acc end loop
!$acc end parallel
#ifndef WARPX
  !$OMP END DO
  !$OMP END PARALLEL
#endif

end subroutine pxrpush_em3d_evec_f

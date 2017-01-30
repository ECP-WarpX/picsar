
! ______________________________________________________________________________
!> @brief
!> This subroutine pushes the electric field with the 3D Yee FDTD 
!> scheme (order 2).
!> This subroutine is adapted for Boxlib and do not contain $OMP parallel 
!> regions.
!
!> @author
!> Weiqun Zhang
!> Mathieu Lobet
!
!> @date
!> Creation 2015
subroutine warpx_pxr_push_em3d_evec( &
     lo, hi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     jx, jxlo, jxhi, &
     jy, jylo, jyhi, &
     jz, jzlo, jzhi, &
     mudt,    &
     dtsdx,dtsdy,dtsdz,&
     norder) bind(c) !#do not parse
! ______________________________________________________________________________

  use constants

  integer, intent(in) :: lo(3), hi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),&
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3),&
       jxlo(3),jxhi(3),jylo(3),jyhi(3),jzlo(3),jzhi(3),&
       norder
       
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

  do l       = lo(3), hi(3)
     do k    = lo(2), hi(2)
        do j = lo(1), hi(1)           
           Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l) - Bz(j,k-1,l  )) &
                                 - dtsdz * (By(j,k,l) - By(j,k  ,l-1)) &
                                 - mudt  * jx(j,k,l)
              
           Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                                 + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                                 - mudt  * jy(j,k,l)

           Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                                 - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                                 - mudt  * jz(j,k,l)
        end do
     end do
  end do

end subroutine warpx_pxr_push_em3d_evec



! ______________________________________________________________________________
!> @brief
!> This subroutine pushes the magnetic field with the 3D Yee FDTD 
!> scheme (order 2).
!> This subroutine is adapted for Boxlib and do not contain $OMP parallel 
!> regions.
!
!> @author
!> Weiqun Zhang
!> Mathieu Lobet
!
!> @date
!> Creation 2015
subroutine warpx_pxr_push_em3d_bvec( &
     lo, hi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     mudt,    &
     dtsdx,dtsdy,dtsdz,&
     norder) bind(c) !#do not parse
! ______________________________________________________________________________

  USE constants

  integer :: lo(3), hi(3), &
       exlo(3),exhi(3),eylo(3),eyhi(3),ezlo(3),ezhi(3),&
       bxlo(3),bxhi(3),bylo(3),byhi(3),bzlo(3),bzhi(3),&
       norder
       
  real(num), intent(IN):: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
  real(num), intent(IN):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
  real(num), intent(IN):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))
  
  real(num), intent(INOUT):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(num), intent(INOUT):: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(num), intent(INOUT):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
  
  real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz

  integer :: j,k,l

  do l       = lo(3), hi(3)
     do k    = lo(2), hi(2)
        do j = lo(1), hi(1)           
              Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j  ,k+1,l  ) - Ez(j,k,l)) &
                                    + dtsdz * (Ey(j  ,k  ,l+1) - Ey(j,k,l))
              By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k  ,l  ) - Ez(j,k,l)) &
                                    - dtsdz * (Ex(j  ,k  ,l+1) - Ex(j,k,l))
              Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k  ,l  ) - Ey(j,k,l)) &
                                    + dtsdy * (Ex(j  ,k+1,l  ) - Ex(j,k,l))
          end do
      end do
  end do

end subroutine warpx_pxr_push_em3d_bvec


! ______________________________________________________________________________
!> @brief
!> This subroutine pushes the electric field with the 2D Yee FDTD 
!> scheme (order 2).
!> This subroutine is adapted for Boxlib and do not contain $OMP parallel 
!> regions.
!
!> @author
!> Weiqun Zhang
!> Mathieu Lobet
!
!> @date
!> Creation 2015
subroutine warpx_pxr_push_em2d_evec( &
     lo, hi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     jx, jxlo, jxhi, &
     jy, jylo, jyhi, &
     jz, jzlo, jzhi, &
     mudt,    &
     dtsdx,dtsdy,dtsdz,&
     norder) bind(c) !#do not parse
! ______________________________________________________________________________

  use constants

  integer, intent(in) :: lo(2), hi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2),&
       bxlo(2),bxhi(2),bylo(2),byhi(2),bzlo(2),bzhi(2),&
       jxlo(2),jxhi(2),jylo(2),jyhi(2),jzlo(2),jzhi(2),&
       norder
       
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

  ! dtsdy should be be used.  It is set to nan by WarpX.

  do k    = lo(2), hi(2)
     do j = lo(1), hi(1)           
        Ex(j,k) = Ex(j,k) - dtsdz * (By(j,k) - By(j,k-1)) &
                          - mudt  * jx(j,k)
              
        Ey(j,k) = Ey(j,k) - dtsdx * (Bz(j,k) - Bz(j-1,k)) &
                          + dtsdz * (Bx(j,k) - Bx(j,k-1)) &
                          - mudt  * jy(j,k)

        Ez(j,k) = Ez(j,k) + dtsdx * (By(j,k) - By(j-1,k  )) &
                          - dtsdy * (Bx(j,k) - Bx(j  ,k-1)) &
                          - mudt  * jz(j,k)
     end do
  end do

end subroutine warpx_pxr_push_em2d_evec



! ______________________________________________________________________________
!> @brief
!> This subroutine pushes the magnetic field with the 2D Yee FDTD 
!> scheme (order 2).
!> This subroutine is adapted for Boxlib and do not contain $OMP parallel 
!> regions.
!
!> @author
!> Weiqun Zhang
!> Mathieu Lobet
!
!> @date
!> Creation 2015
subroutine warpx_pxr_push_em2d_bvec( &
     lo, hi, &
     ex,exlo,exhi,&
     ey,eylo, eyhi, &
     ez,ezlo, ezhi, &
     bx, bxlo, bxhi, &
     by, bylo, byhi, &
     bz, bzlo, bzhi, &
     mudt,    &
     dtsdx,dtsdy,dtsdz,&
     norder) bind(c) !#do not parse
! ______________________________________________________________________________

  USE constants

  integer :: lo(2), hi(2), &
       exlo(2),exhi(2),eylo(2),eyhi(2),ezlo(2),ezhi(2),&
       bxlo(2),bxhi(2),bylo(2),byhi(2),bzlo(2),bzhi(2),&
       norder
       
  real(num), intent(IN):: ex(exlo(1):exhi(1),exlo(2):exhi(2))
  real(num), intent(IN):: ey(eylo(1):eyhi(1),eylo(2):eyhi(2))
  real(num), intent(IN):: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2))
  
  real(num), intent(INOUT):: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
  real(num), intent(INOUT):: by(bylo(1):byhi(1),bylo(2):byhi(2))
  real(num), intent(INOUT):: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2))
  
  real(num), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz

  integer :: j,k

  ! dtsdy should be be used.  It is set to nan by WarpX.

  do k    = lo(2), hi(2)
     do j = lo(1), hi(1)           
        Bx(j,k) = Bx(j,k) + dtsdz * (Ey(j  ,k+1) - Ey(j,k))
        By(j,k) = By(j,k) + dtsdx * (Ez(j+1,k  ) - Ez(j,k)) &
                          - dtsdz * (Ex(j  ,k+1) - Ex(j,k))
        Bz(j,k) = Bz(j,k) - dtsdx * (Ey(j+1,k  ) - Ey(j,k))
     end do
  end do

end subroutine warpx_pxr_push_em2d_bvec

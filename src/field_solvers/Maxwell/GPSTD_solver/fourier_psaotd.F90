! ________________________________________________________________________________________

! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c) 2016,
! The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a paid-up,
! nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute
! copies to the public, prepare derivative works, and perform publicly and display
! publicly, and to permit other to do so.
!
!
! fourier_psaotd.F90
!
! Purpose:
! This file contains subroutines for:
! (i)  fftw init plans   
! (ii) Forward/Backward Fourier transform of EM fields using,
! FFTW (Distributed/Shared version) or P3DFFT
! (iii) Maxwell push in the spectral space using the Pseudo-Spectral Arbitrary Order
! Analytical Time Domain (PSAOTD) Maxwell solver.
!
! Authors:
! Henri Vincenti
! Jean-Luc Vay
!
! Date:
! Creation March, 2017
!
! ________________________________________________________________________________________
MODULE fourier_psaotd
  USE gpstd_solver
  USE matrix_coefficients

  IMPLICIT NONE
  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes FFTW plans 
  !
  !> @author
  !> Haithem Kallala
  !> H. Vincenti 
  !> @date
  !> Creation 2017
  !
  !> @params[in] nopenmp - INTEGER(idp) - number of OpenMP threads/MPI processes  
  ! ______________________________________________________________________________________
  SUBROUTINE init_plans_fourier_mpi(nopenmp)
    USE PICSAR_precision
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE group_parameters
    USE matrix_coefficients
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nopenmp
    INTEGER(C_INT) :: nopenmp_cint, iret
    INTEGER(C_INTPTR_T) :: nx_cint, ny_cint, nz_cint
    INTEGER(idp)        :: i
    INTEGER(isp)        :: planner_flag_1, planner_flag_2
    nopenmp_cint=nopenmp
    IF(.NOT. p3dfft_flag) THEN
      IF  (fftw_threads_ok) THEN
        CALL  DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
      ENDIF
    ENDIF
    !> If fftw_mpi_transpose then use FFTW_MPI_TRANSPOSED_OUT/IN plans
    !> fftw_mpi_transpose avoids spurious mpi_alltoall call for each
    !> fftw_mpi_exec call. (initially fftw_mpi_exec call mpi_alltoall two
    !> times to perform global data transposition along y and z axis back and
    !> forth)
    !> The second mpi_alltoall is ommited when using fftw_mpi_transpose.
    !> Hence fftw_mpi_exec is faster when using transposed plans
    !> But the user should keep in mind that fourier fields are not transposed
    !> back to have the same sizes as real fields since z and y axis are now
    !> switched.

    !> Block matrixes are also transposed conveniently  during init_gpstd when
    !> using transposed plans
    !> A similar optimization is possible when using p3dfft (p3dfft_stride =
    !.TRUE.) but z and x axis are then transposed 

    IF(fftw_mpi_transpose) THEN
      planner_flag_1 = IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_OUT)
      planner_flag_2 = IOR(FFTW_MEASURE,FFTW_MPI_TRANSPOSED_IN)
    ELSE
      planner_flag_1 = FFTW_MEASURE
      planner_flag_2 = FFTW_MEASURE
    ENDIF
    IF(.NOT. fftw_hybrid) THEN
      nz_cint=nz_global
      ny_cint=ny_global
      nx_cint=nx_global
      IF(c_dim == 3) THEN 
        plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint, ny_cint, nx_cint, ex_r, exf,   &
        comm, planner_flag_1)
        plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(nz_cint, ny_cint, nx_cint, exf, ex_r,   &
        comm, planner_flag_2)
      ELSE IF(c_dim == 2) THEN
        plan_r2c_mpi = fftw_mpi_plan_dft_r2c_2d(nz_cint, nx_cint,ex_r, exf,   &
        comm, planner_flag_1)
        plan_c2r_mpi = fftw_mpi_plan_dft_c2r_2d(nz_cint,nx_cint, exf,ex_r,   &
        comm, planner_flag_2)
      ENDIF
    ELSE 
      nz_cint = nz_group
      ny_cint = ny_group
      nx_cint = nx_group
      DO i=1, nb_group
        IF(MPI_COMM_GROUP_ID(i) .NE. MPI_COMM_NULL) THEN
          IF(c_dim == 3) THEN
            plan_r2c_mpi = fftw_mpi_plan_dft_r2c_3d(nz_cint, ny_cint, nx_cint, ex_r,    &
            exf, MPI_COMM_GROUP_ID(i), planner_flag_1)
            plan_c2r_mpi = fftw_mpi_plan_dft_c2r_3d(nz_cint, ny_cint, nx_cint, exf,     &
            ex_r, MPI_COMM_GROUP_ID(i), planner_flag_2)
          ELSE IF(c_dim == 2) THEN
            IF(absorbing_bcs .EQV. .FALSE.) THEN
              plan_r2c_mpi = fftw_mpi_plan_dft_r2c_2d(nz_cint, nx_cint,ex_r,    &
              exf, MPI_COMM_GROUP_ID(i), planner_flag_1)
              plan_c2r_mpi = fftw_mpi_plan_dft_c2r_2d(nz_cint, nx_cint,exf,     &
              ex_r, MPI_COMM_GROUP_ID(i), planner_flag_2)
            ELSE IF(absorbing_bcs) THEN
              plan_r2c_mpi = fftw_mpi_plan_dft_r2c_2d(nz_cint, nx_cint,exy_r,    &
              exf, MPI_COMM_GROUP_ID(i), planner_flag_1)
              plan_c2r_mpi = fftw_mpi_plan_dft_c2r_2d(nz_cint, nx_cint, exf,     &
              exy_r, MPI_COMM_GROUP_ID(i), planner_flag_2)

            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE init_plans_fourier_mpi

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes forward  R2Clocal FFTs - concerns only the local 
  !> pseudo-spectral solver
  !
  !> @author
  !> H. Vincenti 
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE get_Ffields()
    USE shared_data
    USE fields
    USE fourier
    USE fastfft
    USE time_stat
    USE params

    IMPLICIT NONE
    INTEGER(idp) :: nfftx, nffty, nfftz, nxx, nyy, nzz
    INTEGER(idp) :: ix,iy,iz
    REAL(num)    :: tmptime
    nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
    nfftz=nz+2*nzguards+1
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nz+2*nzguards
#endif
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Init fourier fields fields
#if !defined(LIBRARY)
     CALL copy_field_forward()
#endif
   ! reflective bcs after pml
   IF(absorbing_bcs) THEN 
      CALL set_fields_null_guards_cells()
   ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    
    ! Perform local forward FFTs R2C of all grid arrays 
    CALL fft_forward_r2c_local(nfftx,nffty,nfftz)

  END SUBROUTINE get_Ffields

  SUBROUTINE set_fields_null_guards_cells()
    USE fields
    USE shared_data
     
    IMPLICIT NONE
    INTEGER(idp)  :: ix, iy , iz

     IF(.NOT. u_pml) THEN
       IF(absorbing_bcs_x) THEN
         IF(x_min_boundary) THEN
           DO ix = -nxguards,-1
             exy_r(ix,:,:) = 0.0_num
             exz_r(ix,:,:) = 0.0_num
             eyx_r(ix,:,:) = 0.0_num
             eyz_r(ix,:,:) = 0.0_num
             ezx_r(ix,:,:) = 0.0_num  
             ezy_r(ix,:,:) = 0.0_num
             bxy_r(ix,:,:) = 0.0_num
             bxz_r(ix,:,:) = 0.0_num
             byx_r(ix,:,:) = 0.0_num
             byz_r(ix,:,:) = 0.0_num
             bzx_r(ix,:,:) = 0.0_num
             bzy_r(ix,:,:) = 0.0_num
           ENDDO
         ENDIF
         IF(x_max_boundary) THEN
           DO ix=nx,nx+nxguards-1 
             exy_r(ix,:,:) = 0.0_num
             exz_r(ix,:,:) = 0.0_num
             eyx_r(ix,:,:) = 0.0_num
             eyz_r(ix,:,:) = 0.0_num
             ezx_r(ix,:,:) = 0.0_num
             ezy_r(ix,:,:) = 0.0_num
             bxy_r(ix,:,:) = 0.0_num
             bxz_r(ix,:,:) = 0.0_num
             byx_r(ix-1,:,:) = 0.0_num
             byz_r(ix-1,:,:) = 0.0_num
             bzx_r(ix-1,:,:) = 0.0_num
             bzy_r(ix-1,:,:) = 0.0_num
           ENDDO
           byx_r(nx+nxguards-1,:,:) = 0.0_num
           byz_r(nx+nxguards-1,:,:) = 0.0_num
           bzx_r(nx+nxguards-1,:,:) = 0.0_num
           bzy_r(nx+nxguards-1,:,:) = 0.0_num
         ENDIF
       ENDIF
       IF(c_dim == 3) THEN 
         IF(absorbing_bcs_y) THEN
           IF(y_min_boundary) THEN
             DO iy = -nyguards,-1
               exy_r(:,iy,:) = 0.0_num
               exz_r(:,iy,:) = 0.0_num
               eyx_r(:,iy,:) = 0.0_num
               eyz_r(:,iy,:) = 0.0_num
               ezx_r(:,iy,:) = 0.0_num
               ezy_r(:,iy,:) = 0.0_num
               bxy_r(:,iy,:) = 0.0_num
               bxz_r(:,iy,:) = 0.0_num
               byx_r(:,iy,:) = 0.0_num
               byz_r(:,iy,:) = 0.0_num
               bzx_r(:,iy,:) = 0.0_num
               bzy_r(:,iy,:) = 0.0_num
             ENDDO
           ENDIF
           IF(y_max_boundary) THEN
             DO iy=ny,ny+nyguards-1  
               exy_r(:,iy,:) = 0.0_num
               exz_r(:,iy,:) = 0.0_num
               eyx_r(:,iy,:) = 0.0_num
               eyz_r(:,iy,:) = 0.0_num
               ezx_r(:,iy,:) = 0.0_num
               ezy_r(:,iy,:) = 0.0_num
               bxy_r(:,iy-1,:) = 0.0_num
               bxz_r(:,iy-1,:) = 0.0_num
               byx_r(:,iy,:) = 0.0_num
               byz_r(:,iy,:) = 0.0_num
               bzx_r(:,iy-1,:) = 0.0_num
               bzy_r(:,iy-1,:) = 0.0_num
             ENDDO
             bxy_r(:,ny+nyguards-1,:) = 0.0_num
             bxz_r(:,ny+nyguards-1,:) = 0.0_num
             bzx_r(:,ny+nyguards-1,:) = 0.0_num
             bzy_r(:,ny+nyguards-1,:) = 0.0_num
           ENDIF
         ENDIF
       ENDIF
       IF(absorbing_bcs_z) THEN
         IF(z_min_boundary) THEN
           DO iz = -nzguards,-1
             exy_r(:,:,iz) = 0.0_num
             exz_r(:,:,iz) = 0.0_num
             eyx_r(:,:,iz) = 0.0_num
             eyz_r(:,:,iz) = 0.0_num
             ezx_r(:,:,iz) = 0.0_num
             ezy_r(:,:,iz) = 0.0_num
             bxy_r(:,:,iz) = 0.0_num
             bxz_r(:,:,iz) = 0.0_num
             byx_r(:,:,iz) = 0.0_num
             byz_r(:,:,iz) = 0.0_num
             bzx_r(:,:,iz) = 0.0_num
             bzy_r(:,:,iz) = 0.0_num
           ENDDO
         ENDIF
         IF(z_max_boundary) THEN
           DO iz=nz ,nz+nzguards-1  
             exy_r(:,:,iz) = 0.0_num
             exz_r(:,:,iz) = 0.0_num
             eyx_r(:,:,iz) = 0.0_num
             eyz_r(:,:,iz) = 0.0_num
             ezx_r(:,:,iz) = 0.0_num
             ezy_r(:,:,iz) = 0.0_num
             bxy_r(:,:,iz-1) = 0.0_num
             bxz_r(:,:,iz-1) = 0.0_num
             byx_r(:,:,iz-1) = 0.0_num
             byz_r(:,:,iz-1) = 0.0_num
             bzx_r(:,:,iz) = 0.0_num
             bzy_r(:,:,iz) = 0.0_num
           ENDDO
           bxy_r(:,:,nz+nzguards-1) = 0.0_num
           bxz_r(:,:,nz+nzguards-1) = 0.0_num
           byx_r(:,:,nz+nzguards-1) = 0.0_num
           byz_r(:,:,nz+nzguards-1) = 0.0_num
         ENDIF
       ENDIF
     ELSE IF(u_pml) THEN
       IF(absorbing_bcs_x) THEN
         IF(x_min_boundary) THEN
           DO ix = -nxguards,-1
             dex_r(ix,:,:) = 0.0_num
             dey_r(ix,:,:) = 0.0_num
             dez_r(ix,:,:) = 0.0_num  
             bx_r(ix,:,:) = 0.0_num
             by_r(ix,:,:) = 0.0_num
             bz_r(ix,:,:) = 0.0_num
             ex_r(ix,:,:) = 0.0_num
             ey_r(ix,:,:) = 0.0_num
             ez_r(ix,:,:) = 0.0_num
             hx_r(ix,:,:) = 0.0_num
             hy_r(ix,:,:) = 0.0_num
             hz_r(ix,:,:) = 0.0_num

           ENDDO
         ENDIF
         IF(x_max_boundary) THEN
           DO ix=nx,nx+nxguards-1 
             dex_r(ix,:,:) = 0.0_num
             dey_r(ix,:,:) = 0.0_num
             dez_r(ix,:,:) = 0.0_num
             bx_r(ix,:,:) = 0.0_num
             by_r(ix-1,:,:) = 0.0_num
             bz_r(ix-1,:,:) = 0.0_num
             dex_r(ix,:,:) = 0.0_num
             dey_r(ix,:,:) = 0.0_num
             dez_r(ix,:,:) = 0.0_num
             hx_r(ix,:,:) = 0.0_num
             hy_r(ix-1,:,:) = 0.0_num
             hz_r(ix-1,:,:) = 0.0_num

           ENDDO
           by_r(nx+nxguards-1,:,:) = 0.0_num
           bz_r(nx+nxguards-1,:,:) = 0.0_num
           hy_r(nx+nxguards-1,:,:) = 0.0_num
           hz_r(nx+nxguards-1,:,:) = 0.0_num

         ENDIF
       ENDIF
       IF(c_dim == 3) THEN 
         IF(absorbing_bcs_y) THEN
           IF(y_min_boundary) THEN
             DO iy = -nyguards,-1
               dex_r(:,iy,:) = 0.0_num
               dey_r(:,iy,:) = 0.0_num
               dez_r(:,iy,:) = 0.0_num
               bx_r(:,iy,:) = 0.0_num
               by_r(:,iy,:) = 0.0_num
               bz_r(:,iy,:) = 0.0_num
               ex_r(:,iy,:) = 0.0_num
               ey_r(:,iy,:) = 0.0_num
               ez_r(:,iy,:) = 0.0_num
               hx_r(:,iy,:) = 0.0_num
               hy_r(:,iy,:) = 0.0_num
               hz_r(:,iy,:) = 0.0_num

             ENDDO
           ENDIF
           IF(y_max_boundary) THEN
             DO iy=ny,ny+nyguards-1  
                dex_r(:,iy,:) = 0.0_num
                dey_r(:,iy,:) = 0.0_num
                dez_r(:,iy,:) = 0.0_num
                bx_r(:,iy-1,:) = 0.0_num
                by_r(:,iy,:) = 0.0_num
                bz_r(:,iy-1,:) = 0.0_num
                ex_r(:,iy,:) = 0.0_num
                ey_r(:,iy,:) = 0.0_num
                ez_r(:,iy,:) = 0.0_num
                hx_r(:,iy-1,:) = 0.0_num
                hy_r(:,iy,:) = 0.0_num
                hz_r(:,iy-1,:) = 0.0_num

             ENDDO
             bx_r(:,ny+nyguards-1,:) = 0.0_num
             bz_r(:,ny+nyguards-1,:) = 0.0_num
             hx_r(:,ny+nyguards-1,:) = 0.0_num
             hz_r(:,ny+nyguards-1,:) = 0.0_num
           ENDIF
         ENDIF
       ENDIF
       IF(absorbing_bcs_z) THEN
         IF(z_min_boundary) THEN
           DO iz = -nzguards,-1
             dex_r(:,:,iz) = 0.0_num
             dey_r(:,:,iz) = 0.0_num
             dez_r(:,:,iz) = 0.0_num
             bx_r(:,:,iz) = 0.0_num
             by_r(:,:,iz) = 0.0_num
             bz_r(:,:,iz) = 0.0_num
             ex_r(:,:,iz) = 0.0_num
             ey_r(:,:,iz) = 0.0_num
             ez_r(:,:,iz) = 0.0_num
             hx_r(:,:,iz) = 0.0_num
             hy_r(:,:,iz) = 0.0_num
             hz_r(:,:,iz) = 0.0_num

           ENDDO
         ENDIF
         IF(z_max_boundary) THEN
           DO iz=nz ,nz+nzguards-1  
             dex_r(:,:,iz) = 0.0_num
             dey_r(:,:,iz) = 0.0_num
             dez_r(:,:,iz) = 0.0_num
             bx_r(:,:,iz-1) = 0.0_num
             by_r(:,:,iz-1) = 0.0_num
             bz_r(:,:,iz) = 0.0_num
             ex_r(:,:,iz) = 0.0_num
             ey_r(:,:,iz) = 0.0_num
             ez_r(:,:,iz) = 0.0_num
             hx_r(:,:,iz-1) = 0.0_num
             hy_r(:,:,iz-1) = 0.0_num
             hz_r(:,:,iz) = 0.0_num

           ENDDO
           bx_r(:,:,nz+nzguards-1) = 0.0_num
           by_r(:,:,nz+nzguards-1) = 0.0_num
           hx_r(:,:,nz+nzguards-1) = 0.0_num
           hy_r(:,:,nz+nzguards-1) = 0.0_num

         ENDIF
       ENDIF
     ENDIF
  END SUBROUTINE set_fields_null_guards_cells
 
  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes forward R2C distributed FFTs (with MPI groups)
  !
  !> @author
  !> Haithem Kallala
  !> H. Vincenti 
  !> 
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE get_Ffields_mpi_lb()
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    USE field_boundary

    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz,ix1,iy1,iz1, ix2,iy2,iz2, ix3,iy3,iz3
    INTEGER(idp) , DIMENSION(3) :: ubound_w , lbound_w, ubound_e,lbound_e
    REAL(num)    :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    ubound_w = UBOUND(jx)
    lbound_w = LBOUND(jx)
    IF(absorbing_bcs) THEN
      IF(.NOT.  u_pml) THEN
        ubound_e = UBOUND(exy)
        lbound_e = LBOUND(exy)
      ELSE IF(u_pml) THEN
        ubound_e = UBOUND(dex)
        lbound_e = LBOUND(dex)
      ENDIF
    ELSE
      ubound_e = UBOUND(ex)
      lbound_e = LBOUND(ex)
    ENDIF
    ! Performs copies of overlapping portions of local arrays (ex,ey,ez, etc.) 
    ! and FFT arrays (ex_r,ey_r,ez_r etc.) - non overlapping portions requires 
    ! MPI exchanges that are performed further below in this subroutine 

    ! When using periodic bcs, standard EM fields are communicated between local
    ! and hybrid grids
    ! Else, when using absorbing bcs, splitted EM fields are communicated 
    IF(.NOT. absorbing_bcs) THEN 
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, ix1, ix2, ix3, iy1, iy2,  &
      !$OMP  iz1,iz2) COLLAPSE(3)
      DO iz=1,size_exchanges_l2g_recv_z(1)
        DO iy =1,size_exchanges_l2g_recv_y(1)
          DO ix =ix_min_r,ix_max_r
             ix1 = ix
             iy1 = iy-1+g_first_cell_to_recv_y(1) 
             iz1 = iz-1+g_first_cell_to_recv_z(1)
             ix2  =ix-ix_min_r+lbound_e(1)
             iy2 = iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards
             iz2 = iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards
             ix3 = ix-ix_min_r+lbound_w(1)
             iy3 = iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards
             iz3 = iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards
             ex_r(ix1,iy1,iz1) = ex(ix2,iy2,iz2)
             ey_r(ix1,iy1,iz1) = ey(ix2,iy2,iz2)
             ez_r(ix1,iy1,iz1) = ez(ix2,iy2,iz2)
             bx_r(ix1,iy1,iz1) = bx(ix2,iy2,iz2)
             by_r(ix1,iy1,iz1) = by(ix2,iy2,iz2)
             bz_r(ix1,iy1,iz1) = bz(ix2,iy2,iz2)
             jx_r(ix1,iy1,iz1) = jx(ix3,iy3,iz3)
             jy_r(ix1,iy1,iz1) = jy(ix3,iy3,iz3)
             jz_r(ix1,iy1,iz1) = jz(ix3,iy3,iz3)
             rhoold_r(ix1,iy1,iz1) = rhoold(ix3,iy3,iz3)
             rho_r(ix1,iy1,iz1) = rho(ix3,iy3,iz3)
          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ELSE IF(absorbing_bcs) THEN
      IF(.NOT. u_pml) THEN
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
       DO iz=1,size_exchanges_l2g_recv_z(1)
          DO iy =1,size_exchanges_l2g_recv_y(1)
            DO ix =ix_min_r,ix_max_r
               ix1 = ix
               iy1 = iy-1+g_first_cell_to_recv_y(1)
               iz1 = iz-1+g_first_cell_to_recv_z(1)
               ix2 = ix-ix_min_r-nxguards
               iy2 = iy-1+l_first_cell_to_send_y(1) + lbound_e(2)+nyguards
               iz2 = iz-1+l_first_cell_to_send_z(1) + lbound_e(3)+nyguards
               ix3 = ix-ix_min_r+lbound_w(1)
               iy3 = iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards
               iz3 = iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards
               exy_r(ix1,iy1,iz1) = exy(ix2,iy2,iz2)
               exz_r(ix1,iy1,iz1) = exz(ix2,iy2,iz2)
               eyx_r(ix1,iy1,iz1) = eyx(ix2,iy2,iz2)
               eyz_r(ix1,iy1,iz1) = eyz(ix2,iy2,iz2)
               ezx_r(ix1,iy1,iz1) = ezx(ix2,iy2,iz2)
               ezy_r(ix1,iy1,iz1) = ezy(ix2,iy2,iz2)
               bxy_r(ix1,iy1,iz1) = bxy(ix2,iy2,iz2)
               bxz_r(ix1,iy1,iz1) = bxz(ix2,iy2,iz2)
               byx_r(ix1,iy1,iz1) = byx(ix2,iy2,iz2)
               byz_r(ix1,iy1,iz1) = byz(ix2,iy2,iz2)
               bzx_r(ix1,iy1,iz1) = bzx(ix2,iy2,iz2)
               bzy_r(ix1,iy1,iz1) = bzy(ix2,iy2,iz2)
               jx_r(ix1,iy1,iz1) = jx(ix3,iy3,iz3)
               jy_r(ix1,iy1,iz1) = jy(ix3,iy3,iz3)
               jz_r(ix1,iy1,iz1) = jz(ix3,iy3,iz3)
               rhoold_r(ix1,iy1,iz1) = rhoold(ix3,iy3,iz3)
               rho_r(ix1,iy1,iz1) = rho(ix3,iy3,iz3)
            ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
       ELSE IF(u_pml) THEN
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)

        DO iz=1,size_exchanges_l2g_recv_z(1)
          DO iy =1,size_exchanges_l2g_recv_y(1)
            DO ix =ix_min_r,ix_max_r
              ix1 = ix
              iy1 = iy-1+g_first_cell_to_recv_y(1)
              iz1 = iz-1+g_first_cell_to_recv_z(1)
              ix2  =ix-ix_min_r+lbound_e(1)
              iy2 = iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards
              iz2 = iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards
              ix3 = ix-ix_min_r+lbound_w(1)
              iy3 = iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards
              iz3 = iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards
              dex_r(ix1,iy1,iz1) = dex(ix2,iy2,iz2)
              dey_r(ix1,iy1,iz1) = dey(ix2,iy2,iz2)
              dez_r(ix1,iy1,iz1) = dez(ix2,iy2,iz2)
              hx_r(ix1,iy1,iz1) = hx(ix2,iy2,iz2)
              hy_r(ix1,iy1,iz1) = hy(ix2,iy2,iz2)
              hz_r(ix1,iy1,iz1) = hz(ix2,iy2,iz2)
              ex_r(ix1,iy1,iz1) = ex(ix2,iy2,iz2)
              ey_r(ix1,iy1,iz1) = ey(ix2,iy2,iz2)
              ez_r(ix1,iy1,iz1) = ez(ix2,iy2,iz2)
              bx_r(ix1,iy1,iz1) = bx(ix2,iy2,iz2)
              by_r(ix1,iy1,iz1) = by(ix2,iy2,iz2)
              bz_r(ix1,iy1,iz1) = bz(ix2,iy2,iz2)

              jx_r(ix1,iy1,iz1) = jx(ix3,iy3,iz3)
              jy_r(ix1,iy1,iz1) = jy(ix3,iy3,iz3)
              jz_r(ix1,iy1,iz1) = jz(ix3,iy3,iz3)
              rhoold_r(ix1,iy1,iz1) = rhoold(ix3,iy3,iz3)
              rho_r(ix1,iy1,iz1) = rho(ix3,iy3,iz3)
           ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
       ENDIF
     ENDIF
    ! Timers 
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF

    ! Performs MPI exchanges for non-overlapping portions of local and FFT ARRAYS 
    CALL generalized_comms_group_l2g()
    !> Set splitted fields to 0 in the guardcells next to pml region to act as a
    !> reflective mirrorr
    IF(absorbing_bcs) THEN 
      CALL set_fields_null_guardcells_groups()
    ENDIF
   
    CALL fft_forward_r2c_hybrid() 

  END SUBROUTINE get_Ffields_mpi_lb 
  

  SUBROUTINE set_fields_null_guardcells_groups()
    USE fields
    USE shared_data
    USE group_parameters
    INTEGER(idp) :: ix,iy,iz,ixx,iyy,izz
    IF(.NOT. u_pml) THEN
       ! reflective bcs after pml
       IF(absorbing_bcs_x) THEN
         IF(is_group_x_boundary_min) THEN
           DO ix = cell_x_min_g(x_coords+1),-1
             ixx = ix-cell_x_min_g(x_coords+1)+1_idp
             exy_r(ixx,:,:) = 0.0_num
             exz_r(ixx,:,:) = 0.0_num
             eyx_r(ixx,:,:) = 0.0_num
             eyz_r(ixx,:,:) = 0.0_num
             ezx_r(ixx,:,:) = 0.0_num  
             ezy_r(ixx,:,:) = 0.0_num
             bxy_r(ixx,:,:) = 0.0_num
             bxz_r(ixx,:,:) = 0.0_num
             byx_r(ixx,:,:) = 0.0_num
             byz_r(ixx,:,:) = 0.0_num
             bzx_r(ixx,:,:) = 0.0_num
             bzy_r(ixx,:,:) = 0.0_num
           ENDDO
         ENDIF
         IF(is_group_x_boundary_max) THEN
           DO ix=nx_global ,cell_x_max_g(x_coords+1) 
             ixx = ix - cell_x_min_g(x_coords+1)+1_idp 
             exy_r(ixx,:,:) = 0.0_num
             exz_r(ixx,:,:) = 0.0_num
             eyx_r(ixx,:,:) = 0.0_num
             eyz_r(ixx,:,:) = 0.0_num
             ezx_r(ixx,:,:) = 0.0_num
             ezy_r(ixx,:,:) = 0.0_num
             bxy_r(ixx,:,:) = 0.0_num
             bxz_r(ixx,:,:) = 0.0_num
             byx_r(ixx-1,:,:) = 0.0_num
             byz_r(ixx-1,:,:) = 0.0_num
             bzx_r(ixx-1,:,:) = 0.0_num
             bzy_r(ixx-1,:,:) = 0.0_num
           ENDDO
           IF(x_max_boundary) THEN
             byx_r(ixx,:,:) = 0.0_num
             byz_r(ixx,:,:) = 0.0_num
             bzx_r(ixx,:,:) = 0.0_num
             bzy_r(ixx,:,:) = 0.0_num
           ENDIF        
         ENDIF
       ENDIF
       IF(c_dim == 3_idp) THEN
         IF(absorbing_bcs_y) THEN
           IF(is_group_y_boundary_min) THEN
             DO iy = cell_y_min_g(y_coords+1),-1
               iyy = iy-cell_y_min_g(y_coords+1)+1_idp
               exy_r(:,iyy,:) = 0.0_num
               exz_r(:,iyy,:) = 0.0_num
               eyx_r(:,iyy,:) = 0.0_num
               eyz_r(:,iyy,:) = 0.0_num
               ezx_r(:,iyy,:) = 0.0_num
               ezy_r(:,iyy,:) = 0.0_num
               bxy_r(:,iyy,:) = 0.0_num
               bxz_r(:,iyy,:) = 0.0_num
               byx_r(:,iyy,:) = 0.0_num
               byz_r(:,iyy,:) = 0.0_num
               bzx_r(:,iyy,:) = 0.0_num
               bzy_r(:,iyy,:) = 0.0_num
             ENDDO
           ENDIF
           IF(is_group_y_boundary_max) THEN
             DO iy=ny_global ,cell_y_max_g(y_coords+1)  
               iyy = iy - cell_y_min_g(y_coords+1)+1_idp 
               exy_r(:,iyy,:) = 0.0_num
               exz_r(:,iyy,:) = 0.0_num
               eyx_r(:,iyy,:) = 0.0_num
               eyz_r(:,iyy,:) = 0.0_num
               ezx_r(:,iyy,:) = 0.0_num
               ezy_r(:,iyy,:) = 0.0_num
               bxy_r(:,iyy-1,:) = 0.0_num
               bxz_r(:,iyy-1,:) = 0.0_num
               byx_r(:,iyy,:) = 0.0_num
               byz_r(:,iyy,:) = 0.0_num
               bzx_r(:,iyy-1,:) = 0.0_num
               bzy_r(:,iyy-1,:) = 0.0_num
             ENDDO
             IF(y_max_boundary) THEN
               bxy_r(:,iyy,:) = 0.0_num
               bxz_r(:,iyy,:) = 0.0_num
               bzx_r(:,iyy,:) = 0.0_num
               bzy_r(:,iyy,:) = 0.0_num
             ENDIF
           ENDIF
         ENDIF
       ENDIF
       IF(absorbing_bcs_z) THEN
         IF(is_group_z_boundary_min) THEN
           DO iz = cell_z_min_g(z_coords+1),-1
             izz = iz-cell_z_min_g(z_coords+1)+1_idp
             exy_r(:,:,izz) = 0.0_num
             exz_r(:,:,izz) = 0.0_num
             eyx_r(:,:,izz) = 0.0_num
             eyz_r(:,:,izz) = 0.0_num
             ezx_r(:,:,izz) = 0.0_num
             ezy_r(:,:,izz) = 0.0_num
             bxy_r(:,:,izz) = 0.0_num
             bxz_r(:,:,izz) = 0.0_num
             byx_r(:,:,izz) = 0.0_num
             byz_r(:,:,izz) = 0.0_num
             bzx_r(:,:,izz) = 0.0_num
             bzy_r(:,:,izz) = 0.0_num
           ENDDO
         ENDIF
         IF(is_group_z_boundary_max) THEN
           DO iz=nz_global ,cell_z_max_g(z_coords+1)  
             izz = iz - cell_z_min_g(z_coords+1)+1_idp 
             exy_r(:,:,izz) = 0.0_num
             exz_r(:,:,izz) = 0.0_num
             eyx_r(:,:,izz) = 0.0_num
             eyz_r(:,:,izz) = 0.0_num
             ezx_r(:,:,izz) = 0.0_num
             ezy_r(:,:,izz) = 0.0_num
             bxy_r(:,:,izz-1) = 0.0_num
             bxz_r(:,:,izz-1) = 0.0_num
             byx_r(:,:,izz-1) = 0.0_num
             byz_r(:,:,izz-1) = 0.0_num
             bzx_r(:,:,izz) = 0.0_num
             bzy_r(:,:,izz) = 0.0_num
           ENDDO
           IF(z_max_boundary) THEN
             bxy_r(:,:,izz) = 0.0_num
             bxz_r(:,:,izz) = 0.0_num
             byx_r(:,:,izz) = 0.0_num
             byz_r(:,:,izz) = 0.0_num
           ENDIF 
         ENDIF
       ENDIF
       ELSE IF(u_pml) THEN
         ! reflective bcs after pml
         IF(absorbing_bcs_x) THEN
           IF(is_group_x_boundary_min) THEN
             DO ix = cell_x_min_g(x_coords+1),-1
               ixx = ix-cell_x_min_g(x_coords+1)+1_idp
               dex_r(ixx,:,:) = 0.0_num
               dey_r(ixx,:,:) = 0.0_num
               dez_r(ixx,:,:) = 0.0_num  
               hx_r(ixx,:,:) = 0.0_num
               hy_r(ixx,:,:) = 0.0_num
               hz_r(ixx,:,:) = 0.0_num
             ENDDO
           ENDIF
           IF(is_group_x_boundary_max) THEN
             DO ix=nx_global ,cell_x_max_g(x_coords+1) 
               ixx = ix - cell_x_min_g(x_coords+1)+1_idp 
               dex_r(ixx,:,:) = 0.0_num
               dey_r(ixx,:,:) = 0.0_num
               dez_r(ixx,:,:) = 0.0_num
               hx_r(ixx,:,:) = 0.0_num
               hy_r(ixx-1,:,:) = 0.0_num
               hz_r(ixx-1,:,:) = 0.0_num
             ENDDO
             IF(x_max_boundary) THEN
               hy_r(ixx,:,:) = 0.0_num
               hz_r(ixx,:,:) = 0.0_num
             ENDIF        
           ENDIF
         ENDIF
         IF(c_dim == 3_idp) THEN
           IF(absorbing_bcs_y) THEN
             IF(is_group_y_boundary_min) THEN
               DO iy = cell_y_min_g(y_coords+1),-1
                 iyy = iy-cell_y_min_g(y_coords+1)+1_idp
                 dex_r(:,iyy,:) = 0.0_num
                 dey_r(:,iyy,:) = 0.0_num
                 dez_r(:,iyy,:) = 0.0_num
                 hx_r(:,iyy,:) = 0.0_num
                 hy_r(:,iyy,:) = 0.0_num
                 hz_r(:,iyy,:) = 0.0_num
               ENDDO
             ENDIF
             IF(is_group_y_boundary_max) THEN
               DO iy=ny_global ,cell_y_max_g(y_coords+1)  
                 iyy = iy - cell_y_min_g(y_coords+1)+1_idp 
                 dex_r(:,iyy,:) = 0.0_num
                 dey_r(:,iyy,:) = 0.0_num
                 dez_r(:,iyy,:) = 0.0_num
                 hx_r(:,iyy-1,:) = 0.0_num
                 hy_r(:,iyy,:) = 0.0_num
                 hz_r(:,iyy-1,:) = 0.0_num
               ENDDO
               IF(y_max_boundary) THEN
                 hx_r(:,iyy,:) = 0.0_num
                 hz_r(:,iyy,:) = 0.0_num
               ENDIF
             ENDIF
           ENDIF
         ENDIF
         IF(absorbing_bcs_z) THEN
           IF(is_group_z_boundary_min) THEN
             DO iz = cell_z_min_g(z_coords+1),-1
               izz = iz-cell_z_min_g(z_coords+1)+1_idp
               dex_r(:,:,izz) = 0.0_num
               dey_r(:,:,izz) = 0.0_num
               dez_r(:,:,izz) = 0.0_num
               hx_r(:,:,izz) = 0.0_num
               hy_r(:,:,izz) = 0.0_num
               hz_r(:,:,izz) = 0.0_num
             ENDDO
           ENDIF
           IF(is_group_z_boundary_max) THEN
             DO iz=nz_global ,cell_z_max_g(z_coords+1)  
               izz = iz - cell_z_min_g(z_coords+1)+1_idp 
               dex_r(:,:,izz) = 0.0_num
               dey_r(:,:,izz) = 0.0_num
               dez_r(:,:,izz) = 0.0_num
               hx_r(:,:,izz-1) = 0.0_num
               hy_r(:,:,izz-1) = 0.0_num
               hz_r(:,:,izz) = 0.0_num
             ENDDO
             IF(z_max_boundary) THEN
               hx_r(:,:,izz) = 0.0_num
               hy_r(:,:,izz) = 0.0_num
             ENDIF 
           ENDIF
         ENDIF
       ENDIF
  END SUBROUTINE set_fields_null_guardcells_groups

 ! _______________________________________________________________________________________
 !> @brief
 !> This subroutine is used perform forward R2C distributed FFTs with fftw_with_mpi=true
 !> and without MPI groups 
 !> N.B: this routine is deprecated and will be integrally replaced by get_Ffields_mpi_lb 
 !> in the future 
 !>  
 !> @author
 !> Henri Vincenti
 !
 !> @date
 !> Creation 2017
 ! _______________________________________________________________________________________
  SUBROUTINE get_Ffields_mpi
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    USE field_boundary
    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz
    REAL(num)    :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=iz_min_r, iz_max_r
      DO iy=iy_min_r, iy_max_r
        DO ix=ix_min_r, ix_max_r
          ex_r(ix, iy, iz)=ex(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          ey_r(ix, iy, iz)=ey(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          ez_r(ix, iy, iz)=ez(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          bx_r(ix, iy, iz)=bx(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          by_r(ix, iy, iz)=by(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          bz_r(ix, iy, iz)=bz(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          jx_r(ix, iy, iz)=jx(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          jy_r(ix, iy, iz)=jy(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          jz_r(ix, iy, iz)=jz(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          rho_r(ix, iy, iz)=rho(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
          rhoold_r(ix, iy, iz)=rhoold(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    ! Get global Fourier transform of all fields components and currents
    CALL fft_forward_r2c_hybrid

  END SUBROUTINE get_Ffields_mpi

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes backward C2R local FFTs - concerns only the local 
  !> pseudo-spectral solver
  !
  !> @author
  !> H. Vincenti 
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE get_fields()
    USE params
    USE shared_data
    USE fields
    USE fourier
    USE fastfft
    USE time_stat
    IMPLICIT NONE
    REAL(num) :: tmptime
    INTEGER(idp) :: ix, iy, iz, nxx, nyy, nzz, nfftx, nffty, nfftz
    nxx=nx+2*nxguards+1; nyy=ny+2*nyguards+1; nzz=nz+2*nzguards+1;
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
    nfftz=nz+2*nzguards+1
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nz+2*nzguards
#endif


    ! Get Inverse Fourier transform of all fields components 
    CALL fft_backward_c2r_local(nfftx,nffty,nfftz)

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

#if !defined (LIBRARY)
     CALL copy_field_backward
#endif
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE get_fields

 ! _______________________________________________________________________________________
 !> @brief
 !> This subroutine is used perform backward C2R distributed FFTs with fftw_with_mpi=true
 !> and without MPI groups 
 !> N.B: this routine is deprecated and will be integrally replaced by get_fields_mpi_lb 
 !> in the future 
 !> 
 !> @author
 !> Henri Vincenti
 !
 !> @date
 !> Creation 2017
 ! _______________________________________________________________________________________
  SUBROUTINE get_fields_mpi
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    IMPLICIT NONE
    REAL(num) :: tmptime
    INTEGER(idp) :: ix, iy, iz

    ! Get global Fourier transform of all fields components 
    CALL fft_backward_c2r_hybrid
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
    DO iz=iz_min_r, iz_max_r
      DO iy=iy_min_r, iy_max_r
        DO ix=ix_min_r, ix_max_r
          ex(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=ex_r(ix, iy, iz)
          ey(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=ey_r(ix, iy, iz)
          ez(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=ez_r(ix, iy, iz)
          bx(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=bx_r(ix, iy, iz)
          by(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=by_r(ix, iy, iz)
          bz(ix-ix_min_r, iy-iy_min_r, iz-iz_min_r)=bz_r(ix, iy, iz)
        END DO
      END DO
    END DO
      !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE get_fields_mpi


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes backward C2R distributed FFTs (with MPI groups)
  !
  !> @author
  !> Haithem Kallala
  !> H. Vincenti 
  !> 
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
   SUBROUTINE get_fields_mpi_lb
    USE shared_data
    USE fields
    USE mpi_fftw3
    USE time_stat
    USE params
    USE group_parameters
    USE field_boundary
    IMPLICIT NONE
    REAL(num) ::  tmptime
    INTEGER(idp) :: ix, iy, iz,ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3
    INTEGER(idp) , DIMENSION(3) :: lbound_w , ubound_w
    
    ! Perform distributed C2R FFTs of all grid arrays (including fields, currents, charge)
    CALL fft_backward_c2r_hybrid

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    ! Performs copies of overlapping portions of FFT arrays (ex_r,ey_r,ez_r, etc.) 
    ! and local arrays (ex,ey,ez etc.) - non overlapping portions requires 
    ! MPI exchanges that are performed further below in this subroutine

    ! When using periodic bcs, standard EM fields are communicated between local
    ! and hybrid grids
    ! Else, when using absorbing bcs, splitted EM fields are communicated 
    IF(.NOT. absorbing_bcs) THEN
      ubound_w = UBOUND(ex)
      lbound_w = LBOUND(ex)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, ix1, ix2, ix3, iy1, iy2,  &
      !$OMP iy3, iz1, iz2, iz3) COLLAPSE(3)
      DO iz=1,size_exchanges_g2l_recv_z(1)
        DO iy=1,size_exchanges_g2l_recv_y(1)
          DO ix=ix_min_r, ix_max_r
            ix1 = ix-ix_min_r+lbound_w(1)
            iy1 = iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards
            iz1 =  iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards
            ix2 = ix
            iy2 = iy-1+g_first_cell_to_send_y(1)
            iz2 = iz-1+g_first_cell_to_send_z(1)
            ex(ix1,iy1,iz1) = ex_r(ix2,iy2,iz2) 
            ey(ix1,iy1,iz1) = ey_r(ix2,iy2,iz2)
            ez(ix1,iy1,iz1) = ez_r(ix2,iy2,iz2)
            bx(ix1,iy1,iz1) = bx_r(ix2,iy2,iz2)
            by(ix1,iy1,iz1) = by_r(ix2,iy2,iz2)
            bz(ix1,iy1,iz1) = bz_r(ix2,iy2,iz2)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ELSE IF(absorbing_bcs) THEN
      IF(.NOT. u_pml) THEN
        ubound_w = UBOUND(exy)
        lbound_w = LBOUND(exy)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
        DO iz=1,size_exchanges_g2l_recv_z(1)
          DO iy=1,size_exchanges_g2l_recv_y(1)
            DO ix=ix_min_r, ix_max_r
              ix1 = ix-ix_min_r+lbound_w(1)
              iy1 = iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards
              iz1 = iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards
              ix2 = ix
              iy2 = iy-1+g_first_cell_to_send_y(1)
              iz2 = iz-1+g_first_cell_to_send_z(1)
              exy(ix1,iy1,iz1) = exy_r(ix2,iy2,iz2)
              exz(ix1,iy1,iz1) = exz_r(ix2,iy2,iz2)
              eyx(ix1,iy1,iz1) = eyx_r(ix2,iy2,iz2)
              eyz(ix1,iy1,iz1) = eyz_r(ix2,iy2,iz2)
              ezx(ix1,iy1,iz1) = ezx_r(ix2,iy2,iz2)
              ezy(ix1,iy1,iz1) = ezy_r(ix2,iy2,iz2)
              bxy(ix1,iy1,iz1) = bxy_r(ix2,iy2,iz2)
              bxz(ix1,iy1,iz1) = bxz_r(ix2,iy2,iz2)
              byx(ix1,iy1,iz1) = byx_r(ix2,iy2,iz2)
              byz(ix1,iy1,iz1) = byz_r(ix2,iy2,iz2)
              bzx(ix1,iy1,iz1) = bzx_r(ix2,iy2,iz2)
              bzy(ix1,iy1,iz1) = bzy_r(ix2,iy2,iz2)
            END DO
          END DO
        END DO
        !$OMP END PARALLEL DO
      ELSE IF(u_pml) THEN

        ubound_w = UBOUND(dex)
        lbound_w = LBOUND(dex)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
        DO iz=1,size_exchanges_g2l_recv_z(1)
          DO iy=1,size_exchanges_g2l_recv_y(1)
            DO ix=ix_min_r, ix_max_r
              ix1 = ix-ix_min_r+lbound_w(1)
              iy1 = iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards
              iz1 =  iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards
              ix2 = ix
              iy2 = iy-1+g_first_cell_to_send_y(1)
              iz2 = iz-1+g_first_cell_to_send_z(1)
              dex(ix1,iy1,iz1) = dex_r(ix2,iy2,iz2)
              dey(ix1,iy1,iz1) = dey_r(ix2,iy2,iz2)
              dez(ix1,iy1,iz1) = dez_r(ix2,iy2,iz2)
              bx(ix1,iy1,iz1) = bx_r(ix2,iy2,iz2)
              by(ix1,iy1,iz1) = by_r(ix2,iy2,iz2)
              bz(ix1,iy1,iz1) = bz_r(ix2,iy2,iz2)
            END DO
          END DO
        END DO
        !$OMP END PARALLEL DO
      ENDIF
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    ! Performs MPI exchanges for non-overlapping portions of local and FFT ARRAYS 
    CALL generalized_comms_group_g2l()
  END SUBROUTINE get_fields_mpi_lb
  
  SUBROUTINE fft_forward_r2c_local(nfftx,nffty,nfftz) 
    USE fields
    USE fastfft
    USE fourier
    USE time_stat
    USE shared_data
    USE params
    REAL(num)   :: tmptime
    INTEGER(idp), INTENT(IN)    ::   nfftx,nffty,nfftz
    

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    IF(g_spectral) THEN
      IF(.NOT. absorbing_bcs) THEN
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ex_r,                      &
             vold(nmatrixes)%block_vector(1)%block3dc, plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ey_r,                      &
            vold(nmatrixes)%block_vector(2)%block3dc, plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ez_r,                      &
            vold(nmatrixes)%block_vector(3)%block3dc , plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bx_r,                      &
            vold(nmatrixes)%block_vector(4)%block3dc , plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, by_r,                      &
            vold(nmatrixes)%block_vector(5)%block3dc , plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bz_r,                      &  
            vold(nmatrixes)%block_vector(6)%block3dc , plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jx_r,                      &
            vold(nmatrixes)%block_vector(7)%block3dc , plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jy_r,                      &
            vold(nmatrixes)%block_vector(8)%block3dc, plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jz_r,                      &
            vold(nmatrixes)%block_vector(9)%block3dc, plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rhoold_r,                  &
            vold(nmatrixes)%block_vector(10)%block3dc,plan_r2c)
        CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rho_r,                     &
            vold(nmatrixes)%block_vector(11)%block3dc, plan_r2c)
      ELSE IF(absorbing_bcs) THEN
        IF(.NOT. u_pml) THEN
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, exy_r,                      & 
               vold(nmatrixes)%block_vector(1)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, exz_r,                      & 
              vold(nmatrixes)%block_vector(2)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, eyx_r,                      & 
              vold(nmatrixes)%block_vector(3)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, eyz_r,                      & 
              vold(nmatrixes)%block_vector(4)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ezx_r,                      & 
              vold(nmatrixes)%block_vector(5)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ezy_r,                      &  
              vold(nmatrixes)%block_vector(6)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bxy_r,                      & 
               vold(nmatrixes)%block_vector(7)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bxz_r,                      & 
              vold(nmatrixes)%block_vector(8)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, byx_r,                      & 
              vold(nmatrixes)%block_vector(9)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, byz_r,                      & 
              vold(nmatrixes)%block_vector(10)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bzx_r,                      & 
              vold(nmatrixes)%block_vector(11)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bzy_r,                      &  
              vold(nmatrixes)%block_vector(12)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jx_r,                      & 
              vold(nmatrixes)%block_vector(13)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jy_r,                      & 
              vold(nmatrixes)%block_vector(14)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jz_r,                      & 
              vold(nmatrixes)%block_vector(15)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rhoold_r,                  & 
              vold(nmatrixes)%block_vector(16)%block3dc,plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rho_r,                     & 
              vold(nmatrixes)%block_vector(17)%block3dc, plan_r2c) 
        ELSE IF(u_pml) THEN
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, dex_r,                      &
               vold(nmatrixes)%block_vector(1)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, dey_r,                      &
              vold(nmatrixes)%block_vector(2)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, dez_r,                      &
              vold(nmatrixes)%block_vector(3)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bx_r,                      &
              vold(nmatrixes)%block_vector(4)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, by_r,                      &
              vold(nmatrixes)%block_vector(5)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bz_r,                      &
              vold(nmatrixes)%block_vector(6)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ex_r,                      &
               vold(nmatrixes)%block_vector(7)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ey_r,                      &
              vold(nmatrixes)%block_vector(8)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ez_r,                      &
              vold(nmatrixes)%block_vector(9)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, hx_r,                      &
              vold(nmatrixes)%block_vector(10)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, hy_r,                      &
              vold(nmatrixes)%block_vector(11)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, hz_r,                      &
              vold(nmatrixes)%block_vector(12)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jx_r,                      &
              vold(nmatrixes)%block_vector(13)%block3dc , plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jy_r,                      &
              vold(nmatrixes)%block_vector(14)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jz_r,                      &
              vold(nmatrixes)%block_vector(15)%block3dc, plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rhoold_r,                  &
              vold(nmatrixes)%block_vector(16)%block3dc,plan_r2c)
          CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rho_r,                     &
              vold(nmatrixes)%block_vector(17)%block3dc, plan_r2c)
        ENDIF
      ENDIF
    ELSE IF (.NOT. g_spectral) THEN
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ex_r, exf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ey_r, eyf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, ez_r, ezf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bx_r, bxf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, by_r, byf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, bz_r, bzf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jx_r, jxf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jy_r, jyf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, jz_r, jzf, plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rhoold_r, rhooldf,plan_r2c)
      CALL fast_fftw3d_r2c_with_plan(nfftx, nffty, nfftz, rho_r, rhof, plan_r2c) 
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_forward_r2c_local


  SUBROUTINE fft_forward_r2c_hybrid()
    USE fields
    USE fastfft
    USE mpi_fftw3
    USE time_stat
    USE shared_data
    USE params
    USE group_parameters
#if defined(P3DFFT) 
    USE p3dfft
#endif
    REAL(num)   :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
!call test_fftw_mpi_from_python
    IF(g_spectral) THEN
#if defined(P3DFFT)
      IF(.NOT. p3dfft_flag) THEN
#endif
       IF(absorbing_bcs) THEN
         IF(.NOT. u_pml) THEN
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, exy_r,&
           vold(nmatrixes)%block_vector(1)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, exz_r,&
              vold(nmatrixes)%block_vector(2)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, eyx_r,&
              vold(nmatrixes)%block_vector(3)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, eyz_r,&
              vold(nmatrixes)%block_vector(4)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ezx_r,&
              vold(nmatrixes)%block_vector(5)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ezy_r,&
              vold(nmatrixes)%block_vector(6)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bxy_r,&
              vold(nmatrixes)%block_vector(7)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bxz_r,&
              vold(nmatrixes)%block_vector(8)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, byx_r,&
              vold(nmatrixes)%block_vector(9)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, byz_r,&
              vold(nmatrixes)%block_vector(10)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bzx_r,&
              vold(nmatrixes)%block_vector(11)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bzy_r,&
              vold(nmatrixes)%block_vector(12)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r,&
              vold(nmatrixes)%block_vector(13)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r,&
              vold(nmatrixes)%block_vector(14)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r,&
              vold(nmatrixes)%block_vector(15)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r,&
              vold(nmatrixes)%block_vector(16)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r,&
              vold(nmatrixes)%block_vector(17)%block3dc)
         ELSE IF(u_pml) THEN
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, dex_r,&
           vold(nmatrixes)%block_vector(1)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, dey_r,&
              vold(nmatrixes)%block_vector(2)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, dez_r,&
              vold(nmatrixes)%block_vector(3)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, hx_r,&
              vold(nmatrixes)%block_vector(4)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, hy_r,&
              vold(nmatrixes)%block_vector(5)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, hz_r,&
              vold(nmatrixes)%block_vector(6)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r,&
              vold(nmatrixes)%block_vector(7)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r,&
              vold(nmatrixes)%block_vector(8)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r,&
              vold(nmatrixes)%block_vector(9)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r,&
              vold(nmatrixes)%block_vector(10)%block3dc)
           CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r,&
              vold(nmatrixes)%block_vector(11)%block3dc)
         ENDIF
       ELSE IF(.NOT. absorbing_bcs) THEN
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ex_r,                           &
            vold(nmatrixes)%block_vector(1)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ey_r,                           &
            vold(nmatrixes)%block_vector(2)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ez_r,                           &
            vold(nmatrixes)%block_vector(3)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bx_r,                           &
            vold(nmatrixes)%block_vector(4)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, by_r,                           &
            vold(nmatrixes)%block_vector(5)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bz_r,                           &
            vold(nmatrixes)%block_vector(6)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r,                           &
            vold(nmatrixes)%block_vector(7)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r,                           &
            vold(nmatrixes)%block_vector(8)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r,                           &
            vold(nmatrixes)%block_vector(9)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r,                       &
            vold(nmatrixes)%block_vector(10)%block3dc)
         CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r,                          &
            vold(nmatrixes)%block_vector(11)%block3dc)
       ENDIF
#if defined(P3DFFT)
      ELSE IF(p3dfft_flag) THEN
        IF(.NOT. absorbing_bcs) THEN
          CALL p3dfft_ftran_r2c (ex_r,vold(nmatrixes)%block_vector(1)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (ey_r,vold(nmatrixes)%block_vector(2)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (ez_r,vold(nmatrixes)%block_vector(3)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (bx_r,vold(nmatrixes)%block_vector(4)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (by_r,vold(nmatrixes)%block_vector(5)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (bz_r,vold(nmatrixes)%block_vector(6)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (jx_r,vold(nmatrixes)%block_vector(7)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (jy_r,vold(nmatrixes)%block_vector(8)%block3dc,'fft')
          CALL p3dfft_ftran_r2c (jz_r,vold(nmatrixes)%block_vector(9)%block3dc,'fft')
          CALL p3dfft_ftran_r2c(rhoold_r,vold(nmatrixes)%block_vector(10)%block3dc,'fft')
          CALL p3dfft_ftran_r2c(rho_r,vold(nmatrixes)%block_vector(11)%block3dc,'fft')
        ELSE IF (absorbing_bcs) THEN
          IF(.NOT. u_pml) THEN
            CALL p3dfft_ftran_r2c (exy_r,vold(nmatrixes)%block_vector(1)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (exz_r,vold(nmatrixes)%block_vector(2)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (eyx_r,vold(nmatrixes)%block_vector(3)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (eyz_r,vold(nmatrixes)%block_vector(4)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (ezx_r,vold(nmatrixes)%block_vector(5)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (ezy_r,vold(nmatrixes)%block_vector(6)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (bxy_r,vold(nmatrixes)%block_vector(7)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (bxz_r,vold(nmatrixes)%block_vector(8)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (byx_r,vold(nmatrixes)%block_vector(9)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (byz_r,vold(nmatrixes)%block_vector(10)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (bzx_r,vold(nmatrixes)%block_vector(11)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (bzy_r,vold(nmatrixes)%block_vector(12)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (jx_r,vold(nmatrixes)%block_vector(13)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (jy_r,vold(nmatrixes)%block_vector(14)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (jz_r,vold(nmatrixes)%block_vector(15)%block3dc,'fft')
            CALL p3dfft_ftran_r2c(rhoold_r,vold(nmatrixes)%block_vector(16)%block3dc,'fft')
            CALL p3dfft_ftran_r2c(rho_r,vold(nmatrixes)%block_vector(17)%block3dc,'fft')
          ELSE IF(u_pml) THEN
            CALL p3dfft_ftran_r2c (dex_r,vold(nmatrixes)%block_vector(1)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (dey_r,vold(nmatrixes)%block_vector(2)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (dez_r,vold(nmatrixes)%block_vector(3)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (hx_r,vold(nmatrixes)%block_vector(4)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (hy_r,vold(nmatrixes)%block_vector(5)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (hz_r,vold(nmatrixes)%block_vector(6)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (jx_r,vold(nmatrixes)%block_vector(7)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (jy_r,vold(nmatrixes)%block_vector(8)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (jz_r,vold(nmatrixes)%block_vector(9)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (rhoold_r,vold(nmatrixes)%block_vector(10)%block3dc,'fft')
            CALL p3dfft_ftran_r2c (rho_r,vold(nmatrixes)%block_vector(11)%block3dc,'fft')
          ENDIF
        ENDIF
      ENDIF
#endif
    ELSE IF(.NOT. g_spectral) THEN
#if defined(P3DFFT)
      IF(.NOT. p3dfft_flag) THEN
#endif
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ex_r, exf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ey_r, eyf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, ez_r, ezf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bx_r, bxf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, by_r, byf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, bz_r, bzf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jx_r, jxf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jy_r, jyf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, jz_r, jzf)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rho_r, rhof)
        CALL fftw_mpi_execute_dft_r2c(plan_r2c_mpi, rhoold_r, rhooldf)
#if defined(P3DFFT)
      ELSE IF(p3dfft_flag) THEN
        CALL p3dfft_ftran_r2c (ex_r,exf,'fft')
        CALL p3dfft_ftran_r2c (ey_r,eyf,'fft')
        CALL p3dfft_ftran_r2c (ez_r,ezf,'fft')
        CALL p3dfft_ftran_r2c (bx_r,bxf,'fft')
        CALL p3dfft_ftran_r2c (by_r,byf,'fft')
        CALL p3dfft_ftran_r2c (bz_r,bzf,'fft')
        CALL p3dfft_ftran_r2c (jx_r,jxf,'fft')
        CALL p3dfft_ftran_r2c (jy_r,jyf,'fft')
        CALL p3dfft_ftran_r2c (jz_r,jzf,'fft')
        CALL p3dfft_ftran_r2c (rho_r,rhof,'fft')
        CALL p3dfft_ftran_r2c (rhoold_r,rhooldf,'fft')
      ENDIF
#endif
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_forward_r2c_hybrid

  SUBROUTINE fft_backward_c2r_local(nfftx,nffty,nfftz)
    USE fields
    USE fourier
    USE fastfft
    USE time_stat
    USE shared_data
    USE params
    REAL(num)   :: tmptime
    INTEGER(idp), INTENT(IN)     :: nfftx,nffty,nfftz



    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    IF(.NOT. g_spectral) THEN
      CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, exf, ex_r, plan_c2r)
      CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, eyf, ey_r, plan_c2r)
      CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, ezf, ez_r, plan_c2r)
      CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, bxf, bx_r, plan_c2r)
      CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, byf, by_r, plan_c2r)
      CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz, bzf, bz_r, plan_c2r)
    ELSE IF(g_spectral) THEN
      IF(absorbing_bcs) THEN
        IF(.NOT. u_pml) THEN
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           & 
          vnew(nmatrixes)%block_vector(1)%block3dc, exy_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(2)%block3dc, exz_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(3)%block3dc, eyx_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(4)%block3dc, eyz_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(5)%block3dc, ezx_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(6)%block3dc, ezy_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(7)%block3dc, bxy_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(8)%block3dc, bxz_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(9)%block3dc, byx_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(10)%block3dc, byz_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(11)%block3dc, bzx_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(12)%block3dc, bzy_r, plan_c2r)
        ELSE IF(u_pml) THEN
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(1)%block3dc, dex_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(2)%block3dc, dey_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(3)%block3dc, dez_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(4)%block3dc, bx_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(5)%block3dc, by_r, plan_c2r)
          CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
          vnew(nmatrixes)%block_vector(6)%block3dc, bz_r, plan_c2r)
        ENDIF
      ELSE IF(.NOT. absorbing_bcs) THEN
        CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
        vnew(nmatrixes)%block_vector(1)%block3dc, ex_r, plan_c2r)
        CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
        vnew(nmatrixes)%block_vector(2)%block3dc, ey_r, plan_c2r)
        CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           & 
        vnew(nmatrixes)%block_vector(3)%block3dc, ez_r, plan_c2r)
        CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
        vnew(nmatrixes)%block_vector(4)%block3dc, bx_r, plan_c2r)
        CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
        vnew(nmatrixes)%block_vector(5)%block3dc, by_r, plan_c2r)
        CALL fast_fftw3d_c2r_with_plan(nfftx, nffty, nfftz,                           &
        vnew(nmatrixes)%block_vector(6)%block3dc, bz_r, plan_c2r)
      ENDIF
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF



  END SUBROUTINE fft_backward_c2r_local

  SUBROUTINE fft_backward_c2r_hybrid()
    USE fields
    USE fastfft
    USE mpi_fftw3
    USE time_stat
    USE shared_data
    USE params
#if defined(P3DFFT)
    USE p3dfft
#endif
    REAL(num)   :: tmptime

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    IF(.NOT. g_spectral) THEN
#if defined(P3DFFT)
    IF(.NOT. p3dfft_flag) THEN
#endif

      CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, exf, ex_r)
      CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, eyf, ey_r)
      CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, ezf, ez_r)
      CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, bxf, bx_r)
      CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, byf, by_r)
      CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi, bzf, bz_r)

#if defined(P3DFFT)
    ELSE IF(p3dfft_flag) THEN
      CALL p3dfft_btran_c2r (exf,ex_r,'tff')
      CALL p3dfft_btran_c2r (eyf,ey_r,'tff')
      CALL p3dfft_btran_c2r (ezf,ez_r,'tff')
      CALL p3dfft_btran_c2r (bxf,bx_r,'tff')
      CALL p3dfft_btran_c2r (byf,by_r,'tff')
      CALL p3dfft_btran_c2r (bzf,bz_r,'tff')
    ENDIF
#endif

    ELSE IF(g_spectral) THEN
#if defined(P3DFFT)
      IF(.NOT. p3dfft_flag) THEN
#endif
        IF(.NOT. absorbing_bcs) THEN
          CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
          vnew(nmatrixes)%block_vector(1)%block3dc,ex_r)
          CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
          vnew(nmatrixes)%block_vector(2)%block3dc,ey_r)
          CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
          vnew(nmatrixes)%block_vector(3)%block3dc,ez_r)
          CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
          vnew(nmatrixes)%block_vector(4)%block3dc,bx_r)
          CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
          vnew(nmatrixes)%block_vector(5)%block3dc,by_r)
          CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
          vnew(nmatrixes)%block_vector(6)%block3dc,bz_r)
        ELSE IF (absorbing_bcs) THEN
          IF(.NOT. u_pml) THEN
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(1)%block3dc,exy_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(2)%block3dc,exz_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(3)%block3dc,eyx_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(4)%block3dc,eyz_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(5)%block3dc,ezx_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(6)%block3dc,ezy_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(7)%block3dc,bxy_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(8)%block3dc,bxz_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(9)%block3dc,byx_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(10)%block3dc,byz_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(11)%block3dc,bzx_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(12)%block3dc,bzy_r)
          ELSE IF(u_pml) THEN
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(1)%block3dc,dex_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(2)%block3dc,dey_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(3)%block3dc,dez_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(4)%block3dc,bx_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(5)%block3dc,by_r)
            CALL fftw_mpi_execute_dft_c2r(plan_c2r_mpi,                                    &
            vnew(nmatrixes)%block_vector(6)%block3dc,bz_r)
          ENDIF
        ENDIF
#if defined(P3DFFT)
      ELSE IF(p3dfft_flag) THEN
        IF (.NOT. absorbing_bcs) THEN
          CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(1)%block3dc,ex_r,'tff')
          CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(2)%block3dc,ey_r,'tff')
          CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(3)%block3dc,ez_r,'tff')
          CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(4)%block3dc,bx_r,'tff')
          CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(5)%block3dc,by_r,'tff')
          CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(6)%block3dc,bz_r,'tff')
        ELSE IF(absorbing_bcs) THEN
          IF(.NOT. u_pml) THEN
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(1)%block3dc,exy_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(2)%block3dc,exz_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(3)%block3dc,eyx_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(4)%block3dc,eyz_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(5)%block3dc,ezx_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(6)%block3dc,ezy_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(7)%block3dc,bxy_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(8)%block3dc,bxz_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(9)%block3dc,byx_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(10)%block3dc,byz_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(11)%block3dc,bzx_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(12)%block3dc,bzy_r,'tff')
          ELSE IF(u_pml) THEN
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(1)%block3dc,dex_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(2)%block3dc,dey_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(3)%block3dc,dez_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(4)%block3dc,hx_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(5)%block3dc,hy_r,'tff')
            CALL p3dfft_btran_c2r (vnew(nmatrixes)%block_vector(6)%block3dc,hz_r,'tff')
          ENDIF
        ENDIF

      ENDIF
#endif
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_backward_c2r_hybrid

  SUBROUTINE push_psaotd_ebfielfs_2d() 
    USE shared_data
    USE fields
    USE fourier
    USE time_stat
    USE params
    USE mpi_fftw3
    IMPLICIT NONE
    INTEGER(idp) ::  ix, iy, iz, nxx, nzz
    REAL(num) :: tmptime
    COMPLEX(cpx) :: bxfold, byfold, bzfold, exfold, eyfold, ezfold

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    nxx=nkx
    nzz=nkz
    iy=1_idp

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iz, exfold, eyfold, ezfold,     &
    !$OMP bxfold, byfold, bzfold) COLLAPSE(2)
    DO iz=1, nzz
        DO ix=1, nxx
          ! - Bx
          exfold=exf(ix, iy, iz)
          eyfold=eyf(ix, iy, iz)
          ezfold=ezf(ix, iy, iz)
          bxfold=bxf(ix, iy, iz)
          byfold=byf(ix, iy, iz)
          bzfold=bzf(ix, iy, iz)

          bxf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(4, 4)%block3dc(ix, iy,   &
          iz)*bxfold + cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc(ix, iy,        &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc(ix, iy,        &
          iz)*jyf(ix, iy, iz)  
          ! - By
          byf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc(ix, iy,   &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc(ix, iy,        &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(5, 3)%block3dc(ix, iy,        &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc(ix, iy,        &
          iz)*jxf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc(ix,   &
          iy, iz)*jzf(ix, iy, iz)


          ! - Bz
          bzf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc(ix, iy,   &
          iz)*bzfold + cc_mat(nmatrixes)%block_matrix2d(6, 2)%block3dc(ix, iy,        &
          iz)*eyfold +         cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc(ix,    &
          iy, iz)*jyf(ix, iy, iz)

          ! Push E a full time step
          ! - Ex
          exf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(1, 1)%block3dc(ix, iy,   &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(1, 7)%block3dc(ix, iy,        &
          iz)*jxf(ix, iy, iz)     + cc_mat(nmatrixes)%block_matrix2d(1,               &
          11)%block3dc(ix, iy, iz)*rhof(ix, iy, iz) +                                 &
          cc_mat(nmatrixes)%block_matrix2d(1, 10)%block3dc(ix, iy, iz)*rhooldf(ix,    &
          iy, iz)

          ! - Ey
          eyf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(2, 2)%block3dc(ix, iy,   &
          iz)*eyfold  + cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc(ix, iy,        &
          iz)*bxfold  + cc_mat(nmatrixes)%block_matrix2d(2, 6)%block3dc(ix, iy,       &
          iz)*bzfold  + cc_mat(nmatrixes)%block_matrix2d(2, 8)%block3dc(ix, iy,       &
          iz)*jyf(ix, iy, iz) 

          ! - Ez
          ezf(ix, iy, iz) = cc_mat(nmatrixes)%block_matrix2d(3, 3)%block3dc(ix, iy,   &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(3, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(3, 9)%block3dc(ix, iy,        &
          iz)*jzf(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(3, 11)%block3dc(ix,  &
          iy, iz)*rhof(ix, iy, iz) + cc_mat(nmatrixes)%block_matrix2d(3,              &
          10)%block3dc(ix, iy, iz)*rhooldf(ix, iy, iz)


        END DO
    END DO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE push_psaotd_ebfielfs_2d


  SUBROUTINE push_psaotd_ebfielfs_3d()
    USE shared_data
    USE fields
    USE fourier
    USE time_stat
    USE params
    USE mpi_fftw3
    USE group_parameters
    IMPLICIT NONE
    INTEGER(idp) ::  ix, iy, iz, nxx, nyy, nzz
    REAL(num) :: tmptime
    COMPLEX(cpx) :: bxfold, byfold, bzfold, exfold, eyfold, ezfold,&
        jxfold,jyfold,jzfold, rhofold,rhooldfold

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    nxx=nkx
    nyy=nky
    nzz=nkz
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, exfold, eyfold, ezfold,     &
    !$OMP bxfold, byfold, bzfold,jxfold,jyfold,jzfold,rhofold,rhooldfold) COLLAPSE(3)
    DO iz=1, nzz
      DO iy=1, nyy
        DO ix=1, nxx
          ! - Bx
          exfold=exf(ix, iy, iz)
          eyfold=eyf(ix, iy, iz)
          ezfold=ezf(ix, iy, iz)
          bxfold=bxf(ix, iy, iz)
          byfold=byf(ix, iy, iz)
          bzfold=bzf(ix, iy, iz)
          jxfold=jxf(ix, iy, iz)
          jyfold=jyf(ix, iy, iz)
          jzfold=jzf(ix, iy, iz)
          rhofold=rhof(ix, iy, iz)
          rhooldfold=rhooldf(ix, iy, iz)


          bxf(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(4, 4)%block3dc(ix, iy,                     &
          iz)*bxfold + cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc(ix, iy,        &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc(ix, iy,        &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc(ix, iy,        &
          iz)*jyfold + cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc(ix,            &
          iy, iz)*jzfold

          ! - By
          byf(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc(ix, iy,                     &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc(ix, iy,        &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(5, 3)%block3dc(ix, iy,        &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc(ix, iy,        &
          iz)*jxfold + cc_mat(nmatrixes)%block_matrix2d(5, 9)%block3dc(ix,            &
          iy, iz)*jzfold


          ! - Bz
          bzf(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc(ix, iy,                     &
          iz)*bzfold + cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc(ix, iy,        &
          iz)*exfold+ cc_mat(nmatrixes)%block_matrix2d(6, 2)%block3dc(ix, iy,         &
          iz)*eyfold+ cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc(ix, iy,         &
          iz)*jxfold+ cc_mat(nmatrixes)%block_matrix2d(6, 8)%block3dc(ix,             &
          iy, iz)*jyfold

          ! Push E a full time step
          ! - Ex
          exf(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(1, 1)%block3dc(ix, iy,                     &
          iz)*exfold + cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc(ix, iy,        &
          iz)*bzfold + cc_mat(nmatrixes)%block_matrix2d(1, 7)%block3dc(ix, iy,        &
          iz)*jxfold     + cc_mat(nmatrixes)%block_matrix2d(1,                        &
          11)%block3dc(ix, iy, iz)*rhofold       +                                    &
          cc_mat(nmatrixes)%block_matrix2d(1, 10)%block3dc(ix, iy, iz)*rhooldfold

          ! - Ey
          eyf(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(2, 2)%block3dc(ix, iy,                     &
          iz)*eyfold + cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc(ix, iy,        &
          iz)*bxfold  + cc_mat(nmatrixes)%block_matrix2d(2, 6)%block3dc(ix, iy,       &
          iz)*bzfold  + cc_mat(nmatrixes)%block_matrix2d(2, 8)%block3dc(ix, iy,       &
          iz)*jyfold + cc_mat(nmatrixes)%block_matrix2d(2, 11)%block3dc(ix,           &
          iy, iz)*rhofold + cc_mat(nmatrixes)%block_matrix2d(2,                       &
          10)%block3dc(ix, iy, iz)*rhooldfold


          ! - Ez
          ezf(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(3, 3)%block3dc(ix, iy,                     &
          iz)*ezfold + cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc(ix, iy,        &
          iz)*bxfold + cc_mat(nmatrixes)%block_matrix2d(3, 5)%block3dc(ix, iy,        &
          iz)*byfold + cc_mat(nmatrixes)%block_matrix2d(3, 9)%block3dc(ix, iy,        &
          iz)*jzfold + cc_mat(nmatrixes)%block_matrix2d(3, 11)%block3dc(ix,           &
          iy, iz)*rhofold + cc_mat(nmatrixes)%block_matrix2d(3,                       &
          10)%block3dc(ix, iy, iz)*rhooldfold
        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO
    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE push_psaotd_ebfielfs_3d

  SUBROUTINE init_plans_blocks() 
    USE shared_data
    USE fastfft
    Use fourier
    USE fftw3_fortran
    USE fields
    USE omp_lib
    USE params

    INTEGER(idp) :: nfftx, nffty, nfftz, nopenmp
#ifdef _OPENMP
    nopenmp=OMP_GET_MAX_THREADS()
#else
    nopenmp=1
#endif

    IF (fftw_with_mpi) THEN
      nfftx=nx_global
      nffty=ny_global
      nfftz=nz_global
    ELSE
#if defined(LIBRARY)
      nfftx=nx+2*nxguards+1
      nffty=ny+2*nyguards+1
      nfftz=nz+2*nzguards+1
#else
      nfftx=nx+2*nxguards
      nffty=ny+2*nyguards
      nfftz=nz+2*nzguards
#endif
    ENDIF
    !> Init matrix blocks for psatd
    CALL init_gpstd()
    IF(rank==0) WRITE(0, *) 'INIT GPSTD MATRIX DONE'


    !> If g_spectral == .TRUE. then exf is not initialized in mpi_routine.F90
    !> Instead, vector blocks structures are used to store fourier fields
    !> and multiply_mat_vector(GPSTD.F90) is used to push fields in Fourier
    !> space 
    !> exf only points to vold(nmatrixes)%block_vector(1)%block3dc to initialize
    !> fftw plans
    !> If g_spectral == .FALSE. then exf is already allocated 

    IF(g_spectral) THEN
      exf => vold(nmatrixes)%block_vector(1)%block3dc
      IF(absorbing_bcs ) THEN
        IF(.NOT. u_pml)   ex_r => exy_r
      ENDIF
    ENDIF      
    !> Init fftw plans if used
    !> NB: if p3dfft is used, then the init is performed in mpi_routines during
    !> p3dfft_setup
    IF(.NOT. p3dfft_flag ) THEN
      !> if fftw_with_mpi perform fftw_init_plans in the following routine
      IF (fftw_with_mpi) THEN
        CALL init_plans_fourier_mpi(nopenmp)
      !> If local psatd, plans are initialized here
      ELSE IF(.NOT. fftw_with_mpi) THEN
        IF(c_dim ==3) THEN
          CALL fast_fftw_create_plan_r2c_3d_dft(nopenmp, nfftx, nffty, nfftz,ex_r, exf,  &
          plan_r2c, INT(FFTW_MEASURE, idp), INT(FFTW_FORWARD, idp))
          CALL fast_fftw_create_plan_c2r_3d_dft(nopenmp, nfftx, nffty, nfftz, exf,ex_r,  &
          plan_c2r, INT(FFTW_MEASURE, idp), INT(FFTW_BACKWARD, idp))
        ELSE IF(c_dim == 2) THEN
          CALL fast_fftw_create_plan_r2c_2d_dft(nopenmp, nfftx, nfftz, ex_r, exf,        &
          plan_r2c, INT(FFTW_MEASURE, idp), INT(FFTW_FORWARD, idp))
          CALL fast_fftw_create_plan_c2r_2d_dft(nopenmp, nfftx, nfftz, exf, ex_r,        &
          plan_c2r, INT(FFTW_MEASURE, idp), INT(FFTW_BACKWARD, idp))
        ENDIF
      ENDIF
      IF(rank==0) WRITE(0, *) 'INIT GPSTD PLANS DONE'
    ENDIF
  END SUBROUTINE init_plans_blocks
END MODULE fourier_psaotd

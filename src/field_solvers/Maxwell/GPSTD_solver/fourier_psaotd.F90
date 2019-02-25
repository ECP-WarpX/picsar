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
    USE fields, ONLY: ex_r, exf, exy_r
    USE group_parameters, ONLY: mpi_comm_group_id, nx_group, ny_group, nz_group
    USE iso_c_binding
    USE mpi
    USE mpi_fftw3, ONLY: fftw_estimate, fftw_measure, fftw_mpi_plan_dft_c2r_2d,      &
      fftw_mpi_plan_dft_c2r_3d, fftw_mpi_plan_dft_r2c_2d, fftw_mpi_plan_dft_r2c_3d,  &
      fftw_mpi_transposed_in, fftw_mpi_transposed_out, plan_c2r_mpi, plan_r2c_mpi
    USE picsar_precision, ONLY: idp, isp
    USE shared_data, ONLY: absorbing_bcs, c_dim, comm, fftw_hybrid,                  &
      fftw_mpi_transpose, fftw_plan_measure, fftw_threads_ok, nb_group, nx_global,   &
      ny_global, nz_global, p3dfft_flag
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
    IF(fftw_plan_measure) THEN
       planner_flag_1 = FFTW_MEASURE
       planner_flag_2 = FFTW_MEASURE
    ELSE
       planner_flag_1 = FFTW_ESTIMATE
       planner_flag_2 = FFTW_ESTIMATE
    ENDIF
    
    IF(fftw_mpi_transpose) THEN
      planner_flag_1 = IOR(planner_flag_1,FFTW_MPI_TRANSPOSED_OUT)
      planner_flag_2 = IOR(planner_flag_2,FFTW_MPI_TRANSPOSED_IN)
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
    USE fastfft
    USE fields, ONLY: bxy_r, bxz_r, byx_r, byz_r, bzx_r, bzy_r, exy_r, exz_r, eyx_r, &
      eyz_r, ezx_r, ezy_r, nxguards, nyguards, nzguards, shift_x_pml, shift_y_pml,   &
      shift_z_pml
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs, absorbing_bcs_x, absorbing_bcs_y,          &
      absorbing_bcs_z, c_dim, nx, ny, nz, x, x_max_boundary, x_min_boundary, y,      &
      y_max_boundary, y_min_boundary, z, z_max_boundary, z_min_boundary
    USE time_stat, ONLY: localtimes, timestat_itstart

     

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
    IF(absorbing_bcs_x) THEN
      IF(x_min_boundary) THEN
        DO ix = -nxguards,-nxguards+shift_x_pml-1
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
        DO ix=nx+nxguards-shift_x_pml,nx+nxguards-1 
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
          DO iy = -nyguards,-nyguards+shift_y_pml-1
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
          DO iy=ny+nyguards-shift_y_pml,ny+nyguards-1  
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
        DO iz = -nzguards,-nzguards+shift_z_pml-1
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
        DO iz=nz + nzguards-shift_z_pml,nz+nzguards-1  
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
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    
    ! Perform local forward FFTs R2C of all grid arrays 
    CALL fft_forward_r2c_local(nfftx,nffty,nfftz)

  END SUBROUTINE get_Ffields

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
  SUBROUTINE get_Ffields_AM_rz()
    USE fastfft
    USE fields !, ONLY: nxguards, nyguards
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: nmodes, ny, nx
    USE time_stat, ONLY: timestat_itstart, localtimes

    IMPLICIT NONE
    INTEGER(idp) :: nfftx, nffty, nfftz, nxx, nyy, nzz
    REAL(num)    :: tmptime
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
#endif
    nfftz=nmodes
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ! Init fourier fields fields
    write (*,*) "BEFORE COPYING" , MAXVAL(abs(er_c))
#if !defined(LIBRARY)
     !CALL copy_field_forward_AM_rz()
#endif
    write (*,*) "AFTER COPYING" , MAXVAL(abs(er_c))
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF

    ! Perform local forward FFTs C2C of all grid arrays 
    CALL fft_forward_c2c_local_AM_rz(nfftx,nffty,nfftz)

  END SUBROUTINE get_Ffields_AM_rz
 
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
    USE field_boundary
    USE fields, ONLY: bx, bx_r, bxf, bxy, bxy_r, bxz, bxz_r, by, by_r, byf, byx,     &
      byx_r, byz, byz_r, bz, bz_r, bzf, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exf, exy,  &
      exy_r, exz, exz_r, ey, ey_r, eyf, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezf, ezx,  &
      ezx_r, ezy, ezy_r, g_spectral, jx, jx_r, jxf, jy, jy_r, jyf, jz, jz_r, jzf,    &
      nxguards, nyguards, nzguards, rho_r, rhof, rhoold_r, rhooldf
    USE group_parameters, ONLY: cell_x_max_g, cell_x_min_g, cell_y_max_g,            &
      cell_y_min_g, cell_z_max_g, cell_z_min_g, g_first_cell_to_recv_y,              &
      g_first_cell_to_recv_z, is_group_x_boundary_max, is_group_x_boundary_min,      &
      is_group_y_boundary_max, is_group_y_boundary_min, is_group_z_boundary_max,     &
      is_group_z_boundary_min, l_first_cell_to_send_y, l_first_cell_to_send_z,       &
      size_exchanges_l2g_recv_y, size_exchanges_l2g_recv_z
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs, absorbing_bcs_x, absorbing_bcs_y,          &
      absorbing_bcs_z, c_dim, ix_max_r, ix_min_r, nx_global, ny_global, nz_global,   &
      rho, rhoold, x, x_coords, x_max_boundary, y, y_coords, y_max_boundary, z,      &
      z_coords, z_max_boundary
    USE time_stat, ONLY: localtimes, timestat_itstart
    


    IMPLICIT NONE
    INTEGER(idp) :: ix, iy, iz,ixx,iyy,izz
    INTEGER(idp) , DIMENSION(3) :: ubound_w , lbound_w, ubound_e,lbound_e
    REAL(num)    :: tmptime
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF

    ubound_w = UBOUND(jx)
    lbound_w = LBOUND(jx)
    IF(absorbing_bcs) THEN
      ubound_e = UBOUND(exy)
      lbound_e = LBOUND(exy)
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
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=1,size_exchanges_l2g_recv_z(1)
        DO iy =1,size_exchanges_l2g_recv_y(1)
          DO ix =ix_min_r,ix_max_r
             ex_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   ex(ix-ix_min_r+lbound_e(1),iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards)
             ey_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   ey(ix-ix_min_r+lbound_e(1),iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards)
             ez_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   ez(ix-ix_min_r+lbound_e(1),iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards)

             bx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   bx(ix-ix_min_r+lbound_e(1),iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards)
             by_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   by(ix-ix_min_r+lbound_e(1),iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards)
             bz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   bz(ix-ix_min_r+lbound_e(1),iy-1+l_first_cell_to_send_y(1)+lbound_e(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_e(3)+nzguards)
             jx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   jx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             jy_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   jy(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             jz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   jz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             rho_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   rho(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             rhoold_r(ix,iy-1+g_first_cell_to_recv_y(1)                                &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   rhoold(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)

          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
    ELSE IF(absorbing_bcs) THEN

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
     DO iz=1,size_exchanges_l2g_recv_z(1)
       DO iy =1,size_exchanges_l2g_recv_y(1)
         DO ix =ix_min_r,ix_max_r
            exy_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  exy(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            eyx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  eyx(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            ezx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  ezx(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            exz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  exz(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            eyz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  eyz(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            ezy_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  ezy(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            bxy_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  bxy(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            byx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  byx(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            bzx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  bzx(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            bxz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  bxz(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            byz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  byz(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
            bzy_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
            ,iz-1+g_first_cell_to_recv_z(1))=                                        &
                  bzy(ix-ix_min_r-nxguards,iy-1+l_first_cell_to_send_y(1),            &
                  iz-1+l_first_cell_to_send_z(1))
             jx_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   jx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             jy_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   jy(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             jz_r(ix,iy-1+g_first_cell_to_recv_y(1)                                    &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   jz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             rho_r(ix,iy-1+g_first_cell_to_recv_y(1)                                   &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   rho(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)
             rhoold_r(ix,iy-1+g_first_cell_to_recv_y(1)                                &
             ,iz-1+g_first_cell_to_recv_z(1))=                                         &
                   rhoold(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_send_y(1)+lbound_w(2)+nyguards,             &
                   iz-1+l_first_cell_to_send_z(1)+lbound_w(3)+nzguards)

         ENDDO
       ENDDO
     ENDDO
     !$OMP END PARALLEL DO
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
    ENDIF
    CALL fft_forward_r2c_hybrid() 

  END SUBROUTINE get_Ffields_mpi_lb 

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
    USE field_boundary
    USE fields, ONLY: bx, bx_r, by, by_r, bz, bz_r, ex, ex_r, ey, ey_r, ez, ez_r,    &
      jx, jx_r, jy, jy_r, jz, jz_r, rho_r, rhoold_r
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: ix_max_r, ix_min_r, iy_max_r, iy_min_r, iz_max_r,         &
      iz_min_r, rho, rhoold
    USE time_stat, ONLY: localtimes, timestat_itstart
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
    USE fastfft
    USE fields, ONLY: nxguards, nyguards, nzguards
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: nx, ny, nz
    USE time_stat, ONLY: localtimes, timestat_itstart
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
    USE fields, ONLY: bx, bx_r, by, by_r, bz, bz_r, ex, ex_r, ey, ey_r, ez, ez_r
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: ix_max_r, ix_min_r, iy_max_r, iy_min_r, iz_max_r,         &
      iz_min_r
    USE time_stat, ONLY: localtimes, timestat_itstart
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
    USE field_boundary
    USE fields, ONLY: bx, bx_r, bxy, bxy_r, bxz, bxz_r, by, by_r, byx, byx_r, byz,   &
      byz_r, bz, bz_r, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exy, exy_r, exz, exz_r, ey, &
      ey_r, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezx, ezx_r, ezy, ezy_r, nxguards,      &
      nyguards, nzguards
    USE group_parameters, ONLY: g_first_cell_to_send_y, g_first_cell_to_send_z,      &
      l_first_cell_to_recv_y, l_first_cell_to_recv_z, size_exchanges_g2l_recv_y,     &
      size_exchanges_g2l_recv_z
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs, ix_max_r, ix_min_r, y, z
    USE time_stat, ONLY: localtimes, timestat_itstart
    IMPLICIT NONE
    REAL(num) ::  tmptime
    INTEGER(idp) :: ix, iy, iz
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
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=1,size_exchanges_g2l_recv_z(1)
        DO iy=1,size_exchanges_g2l_recv_y(1)
          DO ix=ix_min_r, ix_max_r
            ex(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards) =                                            &
                   ex_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            ey(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   ey_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            ez(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   ez_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            bx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   bx_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            by(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   by_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            bz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   bz_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                  ,iz-1+g_first_cell_to_send_z(1))
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ELSE IF(absorbing_bcs) THEN
      ubound_w = UBOUND(exy)
      lbound_w = LBOUND(exy)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz) COLLAPSE(3)
      DO iz=1,size_exchanges_g2l_recv_z(1)
        DO iy=1,size_exchanges_g2l_recv_y(1)
          DO ix=ix_min_r, ix_max_r

            exy(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards) =                                            &
                   exy_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            eyx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   eyx_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            ezx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   ezx_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            bxy(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   bxy_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            byx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   byx_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            bzx(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   bzx_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                  ,iz-1+g_first_cell_to_send_z(1))
            exz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards) =                                            &
                   exz_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            eyz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   eyz_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            ezy(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   ezy_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            bxz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   bxz_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            byz(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   byz_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                   ,iz-1+g_first_cell_to_send_z(1))

            bzy(ix-ix_min_r+lbound_w(1),iy-1+l_first_cell_to_recv_y(1) + lbound_w(2)+nyguards,                      &
            iz-1+l_first_cell_to_recv_z(1) + lbound_w(3)+nzguards)=                                             &
                   bzy_r(ix,iy-1+g_first_cell_to_send_y(1)                                &
                  ,iz-1+g_first_cell_to_send_z(1))
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(21) = localtimes(21) + (MPI_WTIME() - tmptime)
    ENDIF
    ! Performs MPI exchanges for non-overlapping portions of local and FFT ARRAYS 
    CALL generalized_comms_group_g2l()
  END SUBROUTINE get_fields_mpi_lb
 
  SUBROUTINE fft_forward_c2c_local_AM_rz(nfftx,nffty,nfftz)
    USE fastfft
    USE fields, ONLY: g_spectral
    !USE fields, ONLY : ex, ey, ez, bx, by, bz, jx, jy, jz
    USE fields, ONLY : el_c, er_c, et_c, bl_c, br_c, bt_c, jl_c, jr_c, jt_c,rho_c,   &
                       rhoold_c
    USE fields, ONLY : el_f, ep_f, em_f, bl_f, bp_f, bm_f, jl_f, jp_f, jm_f, rho_f,   &
                       rhoold_f
    USE shared_data, ONLY : absorbing_bcs

    USE fourier, ONLY: plan_rz_f
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE time_stat, ONLY: timestat_itstart, localtimes
    REAL(num)   :: tmptime
    INTEGER(idp), INTENT(IN)    ::   nfftx,nffty,nfftz
    COMPLEX(cpx), dimension(:,:,:), allocatable :: ert_p, ert_m, brt_p,brt_m,jrt_p, jrt_m
    COMPLEX (cpx) :: ii  
    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    ALLOCATE (ert_p(nfftx,nffty,nfftz))
    ALLOCATE (ert_m(nfftx,nffty,nfftz))
    ALLOCATE (brt_p(nfftx,nffty,nfftz))
    ALLOCATE (brt_m(nfftx,nffty,nfftz))
    ALLOCATE (jrt_p(nfftx,nffty,nfftz))
    ALLOCATE (jrt_m(nfftx,nffty,nfftz))
    

    ert_p=DCMPLX(0.0_num,0.0_num)
    ert_m=DCMPLX(0.0_num,0.0_num)
    brt_p=DCMPLX(0.0_num,0.0_num)
    brt_m=DCMPLX(0.0_num,0.0_num)
    jrt_p=DCMPLX(0.0_num,0.0_num)
    jrt_m=DCMPLX(0.0_num,0.0_num)
    ii=DCMPLX(0.0_num, 1.0_num)
    ert_p = (er_c+ii*et_c)/2.0_num
    ert_m = (er_c-ii*et_c)/2.0_num
    brt_p = (br_c+ii*bt_c)/2.0_num
    brt_m = (br_c-ii*bt_c)/2.0_num   
    jrt_p = (jr_c+ii*jt_c)/2.0_num
    jrt_m = (jr_c-ii*jt_c)/2.0_num
    write (*,*) "fft_forward_c2c_local_AM_rz ert_p" , MAXVAL(abs(ert_p))
    write (*,*) "fft_forward_c2c_local_AM_rz ert_m" , MAXVAL(abs(ert_m))
    write (*,*) "fft_forward_c2c_local_AM_rz  brt_p" , MAXVAL(abs(brt_p))
    write (*,*) "fft_forward_c2c_local_AM_rz brt_m" , MAXVAL(abs(brt_m))
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, el_c, el_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, ert_m, ep_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, ert_p, em_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bl_c, bl_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, brt_m, bp_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, brt_p, bm_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, jl_c, jl_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, jrt_m, jp_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, jrt_p, jm_f, plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, rhoold_c,rhoold_f,plan_rz_f)
    CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, rho_c, rho_f, plan_rz_f)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, el_c, el_f, plan_rz_f)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, er_c, ep_f, plan_rz_f)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, et_c, em_f, plan_rz_f)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bl_c, bl_f, plan_rz_f)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, br_c, bp_f, plan_rz_f)
    !CALL fast_fftw1d_3d_array_with_plan(nfftx, nffty, nfftz, bt_c, bm_f, plan_rz_f)
    !write (*,*) "et_f =" , MAXVAL(abs(em_f))
    !write (*,*) "er_f =" , abs(ep_f)
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF
    write (0,*) "nffty ==========", REAL(nffty,num)
    el_f=el_f/REAL(nffty,num)
    ep_f= ep_f/REAL(nffty,num)
    em_f=em_f /REAL(nffty,num)
    bl_f= bl_f/REAL(nffty,num)
    bp_f= bp_f/REAL(nffty,num)
    bm_f=bm_f/REAL(nffty,num)
    DEALLOCATE (ert_p)
    DEALLOCATE (ert_m)
    DEALLOCATE (brt_p)
    DEALLOCATE (brt_m)
    DEALLOCATE (jrt_p)
    DEALLOCATE (jrt_m)
    write (*,*) "em_f =" , MAXVAL(abs(em_f))
    write (*,*) "ep_f =" , MAXVAL(abs(ep_f))
    !write (*,*) "fft_forward_c2c_local_AM_rz er_c" , MAXVAL(abs(er_c))
    !write (*,*) "fft_forward_c2c_local_AM_rz et_c" , MAXVAL(abs(ii*et_c))
  END SUBROUTINE fft_forward_c2c_local_AM_rz

 
  SUBROUTINE fft_forward_r2c_local(nfftx,nffty,nfftz) 
    USE fastfft
    USE fields, ONLY: bx, bx_r, bxf, bxy, bxy_r, bxz, bxz_r, by, by_r, byf, byx,     &
      byx_r, byz, byz_r, bz, bz_r, bzf, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exf, exy,  &
      exy_r, exz, exz_r, ey, ey_r, eyf, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezf, ezx,  &
      ezx_r, ezy, ezy_r, g_spectral, jx, jx_r, jxf, jy, jy_r, jyf, jz, jz_r, jzf,    &
      rho_r, rhof, rhoold_r, rhooldf
    USE fourier, ONLY: plan_r2c
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs
    USE time_stat, ONLY: localtimes, timestat_itstart

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
    USE fastfft
    USE fields, ONLY: bx, bx_r, bxf, bxy, bxy_r, bxz, bxz_r, by, by_r, byf, byx,     &
      byx_r, byz, byz_r, bz, bz_r, bzf, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exf, exy,  &
      exy_r, exz, exz_r, ey, ey_r, eyf, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezf, ezx,  &
      ezx_r, ezy, ezy_r, g_spectral, jx, jx_r, jxf, jy, jy_r, jyf, jz, jz_r, jzf,    &
      rho_r, rhof, rhoold_r, rhooldf
    USE iso_c_binding
    USE mpi
    USE mpi_fftw3, ONLY: fftw_mpi_execute_dft_r2c, plan_r2c_mpi
#if defined(P3DFFT) 
    USE p3dfft
#endif
    USE params, ONLY: it
    USE picsar_precision, ONLY: num
    USE shared_data, ONLY: absorbing_bcs, p3dfft_flag
    USE time_stat, ONLY: localtimes, timestat_itstart

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
    USE fastfft
    USE fields, ONLY: bx, bx_r, bxf, bxy, bxy_r, bxz, bxz_r, by, by_r, byf, byx,     &
      byx_r, byz, byz_r, bz, bz_r, bzf, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exf, exy,  &
      exy_r, exz, exz_r, ey, ey_r, eyf, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezf, ezx,  &
      ezx_r, ezy, ezy_r, g_spectral, jx, jx_r, jxf, jy, jy_r, jyf, jz, jz_r, jzf,    &
      rho_r, rhof, rhoold_r, rhooldf
    USE fourier, ONLY: plan_c2r
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num
    USE shared_data, ONLY: absorbing_bcs
    USE time_stat, ONLY: localtimes, timestat_itstart
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
    USE fastfft
    USE fields, ONLY: bx, bx_r, bxf, bxy, bxy_r, bxz, bxz_r, by, by_r, byf, byx,     &
      byx_r, byz, byz_r, bz, bz_r, bzf, bzx, bzx_r, bzy, bzy_r, ex, ex_r, exf, exy,  &
      exy_r, exz, exz_r, ey, ey_r, eyf, eyx, eyx_r, eyz, eyz_r, ez, ez_r, ezf, ezx,  &
      ezx_r, ezy, ezy_r, g_spectral, jx, jx_r, jxf, jy, jy_r, jyf, jz, jz_r, jzf,    &
      rho_r, rhof, rhoold_r, rhooldf
    USE iso_c_binding
    USE mpi
    USE mpi_fftw3, ONLY: fftw_mpi_execute_dft_c2r, plan_c2r_mpi
#if defined(P3DFFT)
    USE p3dfft
#endif
    USE params, ONLY: it
    USE picsar_precision, ONLY: num
    USE shared_data, ONLY: absorbing_bcs, p3dfft_flag
    USE time_stat, ONLY: localtimes, timestat_itstart


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
        ENDIF

      ENDIF
#endif
    ENDIF
    IF (it.ge.timestat_itstart) THEN
      localtimes(22) = localtimes(22) + (MPI_WTIME() - tmptime)
    ENDIF

  END SUBROUTINE fft_backward_c2r_hybrid

  SUBROUTINE push_psaotd_ebfielfs_2d()
    USE fields, ONLY: bxf, byf, bzf, exf, eyf, ezf, jxf, jyf, jzf, rhof, rhooldf
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: cpx, idp, num
    USE shared_data, ONLY: nkx, nkz
    USE time_stat, ONLY: localtimes, timestat_itstart

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
    USE fields, ONLY: bxf, byf, bzf, exf, eyf, ezf, jxf, jyf, jzf, rhof, rhooldf
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: cpx, idp, num
    USE shared_data, ONLY: nkx, nky, nkz
    USE time_stat, ONLY: localtimes, timestat_itstart

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

  SUBROUTINE push_psaotd_ebfielfs_AM_rz()
    USE fields, ONLY: em_h, jl_h, rhoold_h, rho_h, bl_h, jm_h, ep_h, jp_h, bp_h, bm_h, el_h
    USE iso_c_binding
    USE mpi
    USE params, ONLY: it
    USE picsar_precision, ONLY: idp, num, cpx
    USE shared_data, ONLY: nkx, nky, nmodes
    USE time_stat, ONLY: timestat_itstart, localtimes

    IMPLICIT NONE
    INTEGER(idp) ::  ix, iy, iz, nxx, nyy, nzz
    REAL(num) :: tmptime
    COMPLEX(cpx) :: bl_h_old, bp_h_old, bm_h_old, el_h_old, ep_h_old, em_h_old,&
        jl_h_old,jp_h_old,jm_h_old, rho_h_old,rhoold_h_old

    IF (it.ge.timestat_itstart) THEN
      tmptime = MPI_WTIME()
    ENDIF
    nxx=nkx
    nyy=nky
    nzz=nmodes
    write (0,*) "nxx= ", nxx, "nyy = ", nyy, "nzz =", nzz
    write (0,*), "max of ep_h_old", MAXVAL(ABS (ep_h))
    write (0,*), "max of em_h_old", MAXVAL(ABS (em_h))
    write (0,*), "max of el_h_old", MAXVAL(ABS (el_h))
    write (0,*), "max of bl_h_old ", MAXVAL(ABS (bl_h))
    write (0,*), "max of bm_h_old ", MAXVAL(ABS (bm_h ))
    write (0,*), "max of bp_h_old ", MAXVAL(ABS (bp_h))
    write (0,*), "max of jp_h_old ", MAXVAL(ABS (jp_h ))
    write (0,*), "max of jl_h_old ", MAXVAL(ABS (jl_h ))
    write (0,*), "max of jm_h_old ", MAXVAL(ABS (jm_h ))
    write (0,*), "max of rho_h_old", MAXVAL(ABS (rho_h))
    write (0,*), "max of rhoold_h_old", MAXVAL(ABS (rhoold_h))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix, iy, iz, el_h_old, ep_h_old, em_h_old,     &
    !$OMP bl_h_old, bp_h_old, bm_h_old,jl_h_old,jp_h_old,jm_h_old,rho_h_old,rhoold_h_old) COLLAPSE(3)
    DO iz=1, nzz
      DO iy=1, nyy
        DO ix=1, nxx
          ! - Bx
          el_h_old=el_h(ix, iy, iz)
          ep_h_old=ep_h(ix, iy, iz)
          em_h_old=em_h(ix, iy, iz)
          bl_h_old=bl_h(ix, iy, iz)
          bp_h_old=bp_h(ix, iy, iz)
          bm_h_old=bm_h(ix, iy, iz)
          !jl_h_old=jl_h(ix, iy, iz)
          !jp_h_old=jp_h(ix, iy, iz)
          !jm_h_old=jm_h(ix, iy, iz)
          !rho_h_old=rho_h(ix, iy, iz)
          !rhoold_h_old=rhoold_h(ix, iy, iz)
          jl_h_old=0.0
          jp_h_old=0.0
          jm_h_old=0.0
          rho_h_old=0.0
          rhoold_h_old=0.0

          bl_h(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(4, 4)%block3dc(ix, iy,                     &
          iz)*bl_h_old + cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc(ix, iy,        &
          iz)*ep_h_old + cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc(ix, iy,        &
          iz)*em_h_old + cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc(ix, iy,        &
          iz)*jp_h_old + cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc(ix,            &
          iy, iz)*jm_h_old

          ! - By
          bp_h(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc(ix, iy,                     &
          iz)*bp_h_old + cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc(ix, iy,        &
          iz)*el_h_old + cc_mat(nmatrixes)%block_matrix2d(5, 2)%block3dc(ix, iy,        &
          iz)*ep_h_old + cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc(ix, iy,        &
          iz)*jl_h_old + cc_mat(nmatrixes)%block_matrix2d(5, 8)%block3dc(ix,            &
          iy, iz)*jp_h_old


          ! - Bz
          bm_h(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc(ix, iy,                     &
          iz)*bm_h_old + cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc(ix, iy,        &
          iz)*el_h_old+ cc_mat(nmatrixes)%block_matrix2d(6, 3)%block3dc(ix, iy,         &
          iz)*em_h_old+ cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc(ix, iy,         &
          iz)*jl_h_old+ cc_mat(nmatrixes)%block_matrix2d(6, 9)%block3dc(ix,             &
          iy, iz)*jm_h_old

          ! Push E a full time step
          ! - Ex
          el_h(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(1, 1)%block3dc(ix, iy,                     &
          iz)*el_h_old + cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc(ix, iy,        &
          iz)*bp_h_old + cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc(ix, iy,        &
          iz)*bm_h_old + cc_mat(nmatrixes)%block_matrix2d(1, 7)%block3dc(ix, iy,        &
          iz)*jl_h_old     + cc_mat(nmatrixes)%block_matrix2d(1,                        &
          11)%block3dc(ix, iy, iz)*rho_h_old       +                                    &
          cc_mat(nmatrixes)%block_matrix2d(1, 10)%block3dc(ix, iy, iz)*rhoold_h_old

          ! - Ey
          ep_h(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(2, 2)%block3dc(ix, iy,                     &
          iz)*ep_h_old + cc_mat(nmatrixes)%block_matrix2d(2, 4)%block3dc(ix, iy,        &
          iz)*bl_h_old  + cc_mat(nmatrixes)%block_matrix2d(2, 5)%block3dc(ix, iy,       &
          iz)*bm_h_old  + cc_mat(nmatrixes)%block_matrix2d(2, 8)%block3dc(ix, iy,       &
          iz)*jp_h_old + cc_mat(nmatrixes)%block_matrix2d(2, 11)%block3dc(ix,           &
          iy, iz)*rho_h_old + cc_mat(nmatrixes)%block_matrix2d(2,                       &
          10)%block3dc(ix, iy, iz)*rhoold_h_old


          ! - Ez
          em_h(ix, iy, iz) =                                                           &
          cc_mat(nmatrixes)%block_matrix2d(3, 3)%block3dc(ix, iy,                     &
          iz)*em_h_old + cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc(ix, iy,        &
          iz)*bl_h_old + cc_mat(nmatrixes)%block_matrix2d(3, 6)%block3dc(ix, iy,        &
          iz)*bm_h_old + cc_mat(nmatrixes)%block_matrix2d(3, 9)%block3dc(ix, iy,        &
          iz)*jm_h_old + cc_mat(nmatrixes)%block_matrix2d(3, 11)%block3dc(ix,           &
          iy, iz)*rho_h_old + cc_mat(nmatrixes)%block_matrix2d(3,                       &
          10)%block3dc(ix, iy, iz)*rhoold_h_old

          IF (cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc(ix,iy,iz)/= cc_mat(nmatrixes) &
             %block_matrix2d(5, 5)%block3dc(ix,iy,iz))THEN
            write (0,*) "(5,5)   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF 
          IF (cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc(ix,iy,iz)/= cc_mat(nmatrixes)&
              %block_matrix2d(5, 1)%block3dc(ix,iy,iz))THEN
            write (0,*) "(5,1)   ","ix=", ix, "iy=", iy, "iz=", iz
          END IF 
          IF (cc_mat(nmatrixes)%block_matrix2d(5, 2)%block3dc(ix,iy,iz)/= cc_mat(nmatrixes)&
              %block_matrix2d(5, 2)%block3dc(ix,iy,iz))THEN
            write (0,*) "(5,2)   ","ix=", ix, "iy=", iy, "iz=", iz
          END IF 
          IF (cc_mat(nmatrixes)%block_matrix2d(5, 7)%block3dc(ix,iy,iz)/= cc_mat(nmatrixes)&
              %block_matrix2d(5, 7)%block3dc(ix,iy,iz))THEN
            write (0,*) "(5,7)   ","ix=", ix, "iy=", iy, "iz=", iz
          END IF 
          IF (cc_mat(nmatrixes)%block_matrix2d(5, 8)%block3dc(ix,iy,iz)/= cc_mat(nmatrixes)&
             %block_matrix2d(5, 8)%block3dc(ix,iy,iz))THEN
            write (0,*) "(5,8)   ","ix=", ix, "iy=", iy, "iz=", iz
          END IF
          IF (cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc(ix,iy,iz)/= cc_mat(nmatrixes)&
            %block_matrix2d(6, 6)%block3dc(ix,iy,iz))THEN
            write (0,*) "(6,6)   ",  "ix=", ix, "iy=", iy, "iz=", iz
          END IF
          IF (ep_h_old /= ep_h_old) THEN 
            write (0,*) "ep_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz 
          END IF 
          IF (el_h_old /= el_h_old) THEN 
            write (0,*) "el_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
          IF (bp_h_old /= bp_h_old) THEN 
            write (0,*) "bp_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
          IF (jp_h_old /= jp_h_old) THEN 
            write (0,*) "jp_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
          IF (jl_h_old /= jl_h_old) THEN 
            write (0,*) "jl_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
           IF (bl_h_old /= bl_h_old) THEN
            write (0,*) "bl_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
           IF (em_h_old /= em_h_old) THEN
            write (0,*) "em_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
           IF (bm_h_old /= bm_h_old) THEN
            write (0,*) "bm_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
           IF (jm_h_old /= jm_h_old) THEN
            write (0,*) "jm_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
           IF (rho_h_old /= rho_h_old) THEN
            write (0,*) "bl_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
           IF (rhoold_h_old /= rhoold_h_old) THEN
            write (0,*) "rhoold_h_old   ", "ix=", ix, "iy=", iy, "iz=", iz
          END IF
        END DO
      END DO
    END DO
    
   !$OMP END PARALLEL DO 
    write (0,*), "(3,3)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(3, 3)%block3dc))
    write (0,*), "(3,4)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(3, 4)%block3dc))
    write (0,*), "(3,6)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(3, 6)%block3dc))
    write (0,*), "(3,9)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(3, 9)%block3dc))
    write (0,*), "(2,10)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(3, 10)%block3dc))
    write (0,*), "(2,11)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(2, 11)%block3dc))
    write (0,*), "(2,2)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(2, 2)%block3dc))
    write (0,*), "(2,4)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(2,4)%block3dc))
    write (0,*), "(2,5)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(2, 5)%block3dc))
    write (0,*), "(1,1)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(1, 1)%block3dc))
    write (0,*), "(1,5)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(1, 5)%block3dc))
    write (0,*), "(1,6)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(1, 6)%block3dc))
    write (0,*), "(1,7)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(1, 7)%block3dc))
    write (0,*), "(1,10)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(1, 10)%block3dc))
    write (0,*), "(1,11)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(1, 11)%block3dc))
    write (0,*), "(6,9)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(6, 9)%block3dc))
    write (0,*), "(6,7)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(6, 7)%block3dc))
    write (0,*), "(6,3)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(6, 3)%block3dc))
    write (0,*), "(6,1)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(6, 1)%block3dc))
    write (0,*), "(6,6)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(6, 6)%block3dc))
    write (0,*), "(5,5)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(5, 5)%block3dc))
    write (0,*), "(5,8)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(5, 8)%block3dc))
    write (0,*), "(5,7)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(5,7)%block3dc))
    write (0,*), "(5,2)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(5, 2)%block3dc))
    write (0,*), "(5,1)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(5, 1)%block3dc))
    write (0,*), "(4,4)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(4, 4)%block3dc))
    write (0,*), "(4,9)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(4, 9)%block3dc))
    write (0,*), "(4,8)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(4, 8)%block3dc))
    write (0,*), "(4,3)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(4, 3)%block3dc))  
    write (0,*), "(4,2)", MAXVAL(ABS (cc_mat(nmatrixes)%block_matrix2d(4, 2)%block3dc))

    !write (0,*), "ccmat2 ", MAXVAL(ABS (bl_h_old ))
    !write (0,*), "ccmat3 ", MAXVAL(ABS (bm_h_old ))
    !write (0,*), "max of jp_h_old ", MAXVAL(ABS (jp_h_old ))
    !write (0,*), "max of rho_h_old", MAXVAL(ABS (rho_h_old))
    !write (0,*), "max of rhoold_h_old", MAXVAL(ABS (rhoold_h_old))

    IF (it.ge.timestat_itstart) THEN
      localtimes(23) = localtimes(23) + (MPI_WTIME() - tmptime)
    ENDIF
  END SUBROUTINE push_psaotd_ebfielfs_AM_rz


  SUBROUTINE init_plans_blocks()
    USE fastfft
    USE fftw3_fortran, ONLY: fftw_backward, fftw_forward, fftw_measure
    USE fields, ONLY: ex_r, exf, exy_r, g_spectral, nxguards, nyguards, nzguards
    USE fourier, ONLY: plan_c2r, plan_r2c
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: absorbing_bcs, c_dim, fftw_with_mpi, nx, nx_global, ny,   &
      ny_global, nz, nz_global, p3dfft_flag, rank

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
      IF(absorbing_bcs) THEN
        ex_r => exy_r
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


  SUBROUTINE init_plans_blocks_rz()

!    USE fastfft
!    USE fftw3_fortran, ONLY: fftw_measure, fftw_backward, fftw_forward
!    USE fields, ONLY: ex_r, nzguards, nxguards, g_spectral, nyguards, exf, exy_r
!    USE fourier, ONLY: plan_c2r, plan_r2c
!    USE iso_c_binding
!    USE omp_lib
!    USE picsar_precision, ONLY: idp
!    USE shared_data, ONLY: nz, ny, fftw_with_mpi, nx, nx_global, p3dfft_flag,
!&
!     ny_global, c_dim, nz_global, rank, absorbing_bcs


    USE fastfft
    USE fftw3_fortran, ONLY: fftw_measure, fftw_backward, fftw_forward
    USE fields, ONLY: el_c, nxguards, nyguards, el_f
    USE fourier, ONLY: plan_rz_f, plan_rz_f_inv
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: idp
    USE shared_data, ONLY: nmodes, ny, nx, &
      c_dim, rank

    INTEGER(idp) :: nfftx, nffty,nfftz, nopenmp
#ifdef _OPENMP
    nopenmp=OMP_GET_MAX_THREADS()
#else
    nopenmp=1
#endif
#if defined(LIBRARY)
    nfftx=nx+2*nxguards+1
    nffty=ny+2*nyguards+1
    nfftz=nmodes
#else
    nfftx=nx+2*nxguards
    nffty=ny+2*nyguards
    nfftz=nmodes

#endif

!> Init matrix blocks for psatd
    CALL init_gpstd()
    IF(rank==0) WRITE(0, *) 'INIT GPSTD MATRIX DONE RZ'

    CALL fast_fftw_create_plan_1d_3darray_dft(nopenmp, nfftx,nffty,nfftz,el_c,el_f,  &
          plan_rz_f, INT(FFTW_MEASURE, idp), INT(FFTW_FORWARD, idp))

    CALL fast_fftw_create_plan_1d_3darray_dft(nopenmp, nfftx, nffty, nfftz,el_f,el_c,  &
          plan_rz_f_inv, INT(FFTW_MEASURE,idp), INT(FFTW_BACKWARD,idp))

    IF(rank==0) WRITE(0, *) 'INIT GPSTD PLANS DONE RZ'

  END SUBROUTINE init_plans_blocks_rz

END MODULE fourier_psaotd

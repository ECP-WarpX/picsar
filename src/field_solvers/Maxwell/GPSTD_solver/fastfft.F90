!________________________________________________________________________________________
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
! fastfft.F90 
!
! This file contains Fortran subroutines and modules to interface PICSAR 
! With the FFTW3 library and perform 3D Complex/Real to Complex/Real 1D, 2D and 3D
! FFTs/IFFTs using the shared memory version of FFTW (no MPI). It is mainly made of  
! Fortran routines written for the HPC ython module fftw_for_py written by H. Vincenti. 
!
! Developers:
! Henri Vincenti
!
! Date:
! Creation April 2017 
!
! ________________________________________________________________________________________

!**********************************************
!* SECTION 1: Module definitions
!**********************************************

! Module that defines parameters and function interfaces
! For FFTW Fortran wrapping
MODULE fftw3_fortran !#do not parse
  use, intrinsic :: iso_c_binding
  use PICSAR_precision 
  include 'fftw3.f03'
  integer(idp), parameter :: nmaxplan=100000
  ! > Fortran Integer Array where C integer pointers to plans are stored
  integer(idp), DIMENSION(nmaxplan) :: plans_cint
  integer(C_INT) :: nplan=0
END MODULE fftw3_fortran

MODULE mpi_fftw3 !#do not parse
  use, intrinsic :: iso_c_binding
  use PICSAR_precision 
  include 'fftw3-mpi.f03'
  integer(idp), parameter :: nmaxplan_mpi=100000
  ! > Fortran Integer Array where C integer pointers to plans are stored
  integer(idp), DIMENSION(nmaxplan_mpi) :: plans_cint_mpi
  integer(C_INT) :: nplan_mpi=0
  ! - Starting indexes and dimensions along X,Y,Z 
  ! - of the input FFT arrays 
  INTEGER(C_INTPTR_T) :: alloc_local ! used for memory alloc in FFTW
  INTEGER(C_INTPTR_T) :: local_x0, local_nx
  INTEGER(C_INTPTR_T) :: local_y0, local_ny
  INTEGER(C_INTPTR_T) :: local_z0, local_nz
  ! - Starting indexes and dimensions along X,Y,Z 
  ! - of the  output FFT array (tranpsosed or not)
  ! (these parameters are different from the ones 
  ! - of real arrays above, when fftw_mpi_transpose=.TRUE.)
  INTEGER(C_INTPTR_T) :: local_x0_tr, local_nx_tr
  INTEGER(C_INTPTR_T) :: local_y0_tr, local_ny_tr
  INTEGER(C_INTPTR_T) :: local_z0_tr, local_nz_tr
  TYPE(C_PTR) :: plan_r2c_mpi, plan_c2r_mpi 
END MODULE mpi_fftw3

!**********************************************
!* SECTION 2: PLAN CREATION ROUTINES 1D,2D,3D, ND
!* C2C, R2C and C2R
!**********************************************

MODULE fastfft !# do not parse

IMPLICIT NONE 

CONTAINS 

! Subroutine that creates a 3D complex to cpomplex plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_3d_dft(nopenmp,nx,ny,nz,array_in,array_out, &
    plan,plan_type,dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp

    IMPLICIT NONE

    INTEGER(idp), INTENT(IN) ::  nopenmp,nx,ny,nz
    COMPLEX(cpx), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_in, array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret, plan_type_cint, dir_cint, &
		   nx_cint, ny_cint, nz_cint, nopenmp_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    ny_cint=ny
    nz_cint=nz
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_3D(plans_cint(nplan), nx_cint,ny_cint,nz_cint, &
                            array_in,array_out,  &
                            dir_cint,plan_type_cint)
    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_3d_dft


! Subroutine that creates a 3D real to complex plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_r2c_3d_dft(nopenmp,nx,ny,nz,array_in, &
    array_out,plan,plan_type,dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx,ny,nz
    REAL(num), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx/2+1,ny,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret, nopenmp_cint, nx_cint, ny_cint, nz_cint, &
		      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    ny_cint=ny
    nz_cint=nz
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_R2C_3D(plans_cint(nplan), nx_cint,ny_cint,nz_cint, &
                            array_in,array_out,  &
                            plan_type_cint,dir_cint)
    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_r2c_3d_dft

! Subroutine that creates a 3D complex to real plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_c2r_3d_dft(nopenmp,nx,ny,nz,array_in, &
    array_out,plan,plan_type,dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx,ny,nz
    REAL(num), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_out
    COMPLEX(cpx), DIMENSION(nx/2+1,ny,nz), INTENT(IN OUT)  :: array_in
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint, ny_cint, nz_cint, &
                      plan_type_cint, dir_cint
    ! Conversion integer idp to C_INT
    nx_cint=nx
    ny_cint=ny
    nz_cint=nz
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_C2R_3D(plans_cint(nplan), nx_cint,ny_cint,nz_cint, &
                            array_in,array_out,  &
                            plan_type_cint,dir_cint)
    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_c2r_3d_dft


! Subroutine that creates a 2D plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_2d_dft(nopenmp,nx,nz,array_in, &
    array_out,plan,plan_type,dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx,nz
    COMPLEX(cpx), DIMENSION(nx,nz), INTENT(IN OUT)  :: array_in, array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint, nz_cint, &
                      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    nz_cint=nz
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_2D(plans_cint(nplan), nx_cint,nz_cint, array_in, &
                            array_out,                                    &
                            dir_cint,plan_type_cint)

    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_2d_dft

! Subroutine that creates a 2D real to complex plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_r2c_2d_dft(nopenmp,nx,nz,array_in,&
    array_out,plan,plan_type,dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx,nz
    REAL(num), DIMENSION(nx,nz), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx/2+1,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint, nz_cint, &
                      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    nz_cint=nz
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_R2C_2D(plans_cint(nplan), nx_cint,nz_cint, array_in, &
                            array_out,                                        &
                            plan_type_cint,dir_cint)
    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_r2c_2d_dft

! Subroutine that creates a 2D complex to real plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_c2r_2d_dft(nopenmp,nx,nz,array_in, &
    array_out,plan,plan_type, dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx,nz
    REAL(num), DIMENSION(nx,nz), INTENT(IN OUT)  :: array_out
    COMPLEX(cpx), DIMENSION(nx/2+1,nz), INTENT(IN OUT)  :: array_in
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint, nz_cint, &
                      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    nz_cint=nz
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_C2R_2D(plans_cint(nplan),nx_cint,nz_cint, array_in, &
                            array_out,                                       &
                            plan_type_cint,dir_cint)

    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_c2r_2d_dft

! Subroutine that creates a 1D plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_1d_dft(nopenmp,nx,array_in,array_out, &
    plan,plan_type, dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx
    COMPLEX(cpx), DIMENSION(nx), INTENT(IN OUT)  :: array_in, array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint,        &
                      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_1D(plans_cint(nplan), nx_cint,array_in,array_out,  &
                            dir_cint,plan_type_cint)

    ! return index of plan
    plan(1)=nplan
END SUBROUTINE fast_fftw_create_plan_1d_dft

! Subroutine that creates a 1D real to complex plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_r2c_1d_dft(nopenmp,nx,array_in,array_out, &
    plan,plan_type, dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx
    REAL(num), DIMENSION(nx), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx/2+1), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint,        &
                      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_R2C_1D(plans_cint(nplan), nx_cint, array_in,array_out,&
                            plan_type_cint,dir_cint)

    ! return index of plan
    plan(1)=nplan

END SUBROUTINE fast_fftw_create_plan_r2c_1d_dft

! Subroutine that creates a 1D complex to real plan using plan_type optimization
! plan_type can either be FFTW_ESTIMATE (low overhead, low optimization)
! FFTW_MEASURE (moderate to high overhead, high optimization)
! FFTW_EXHAUSTIVE (very high overhead, brute force optimization)
SUBROUTINE fast_fftw_create_plan_c2r_1d_dft(nopenmp,nx,array_in,array_out, &
    plan,plan_type,dir)
    USE fftw3_fortran, ONLY: nplan, plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nopenmp, nx
    REAL(num), DIMENSION(nx), INTENT(IN OUT)  :: array_out
    COMPLEX(cpx), DIMENSION(nx/2+1), INTENT(IN OUT)  :: array_in
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp), INTENT(IN) :: plan_type, dir
    INTEGER(C_INT) :: iret,nopenmp_cint, nx_cint,        &
                      plan_type_cint, dir_cint

    ! Conversion integer idp to C_INT
    nx_cint=nx
    nopenmp_cint=nopenmp
    dir_cint=dir
    plan_type_cint=plan_type

    ! Plan creation
    nplan=nplan+1
    CALL DFFTW_INIT_THREADS(iret)
    CALL DFFTW_PLAN_WITH_NTHREADS(nopenmp_cint)
    CALL DFFTW_PLAN_DFT_C2R_1D(plans_cint(nplan), nx_cint,   &
                            array_in, array_out,              &
                            plan_type_cint,dir_cint)
    ! return index of plan
    plan(1)=nplan

END SUBROUTINE fast_fftw_create_plan_c2r_1d_dft

!**********************************************
!* SECTION 3: FFTW EXECUTION 1D/2D/3D
!* C2C, C2R, R2C
!**********************************************

! --------- 3D ROUTINES
! Subroutine that perform 3D Complex to Complex
! DFT along previously defined strategy "plan"
SUBROUTINE fast_fftw3d_with_plan(nx,ny,nz,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp
    INTEGER(idp), INTENT(IN) ::  nx, ny,nz
    COMPLEX(cpx), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw3d_with_plan

! Subroutine that perform 3D Real 
! to Complex DFT along previously defined strategy "plan"
SUBROUTINE fast_fftw3d_r2c_with_plan(nx,ny,nz,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nx, ny,nz
    REAL(num), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx/2+1,ny,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT_R2C(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw3d_r2c_with_plan

! Subroutine that perform 3D Real to Complex DFT along
! previously defined strategy "plan"
SUBROUTINE fast_fftw3d_c2r_with_plan(nx,ny,nz,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nx, ny,nz
    COMPLEX(cpx), DIMENSION(nx/2+1,ny,nz), INTENT(IN OUT)  :: array_in
    REAL(num), DIMENSION(nx,ny,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT_C2R(plans_cint(plan), array_in, array_out)
END SUBROUTINE fast_fftw3d_c2r_with_plan



! --------- 2D ROUTINES
! Subroutine that perform 2D Complex to Complex DFT
! along previously defined strategy "plan"
SUBROUTINE fast_fftw2d_with_plan(nx,nz,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp
    INTEGER(idp), INTENT(IN) ::  nx,nz
    COMPLEX(cpx), DIMENSION(nx,nz), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw2d_with_plan


! Subroutine that perform 2D Real to Complex 
! DFT along previously defined strategy "plan"
SUBROUTINE fast_fftw2d_r2c_with_plan(nx,nz,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nx,nz
    REAL(num), DIMENSION(nx,nz), INTENT(IN)  :: array_in
    COMPLEX(cpx), DIMENSION(nx/2+1,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT_R2C(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw2d_r2c_with_plan

! Subroutine that perform 2D Real to Complex 
! DFT along previously defined strategy "plan"
SUBROUTINE fast_fftw2d_c2r_with_plan(nx,nz,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) :: nx,nz
    COMPLEX(cpx), DIMENSION(nx/2+1,nz), INTENT(IN OUT)  :: array_in
    REAL(num), DIMENSION(nx,nz), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT_C2R(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw2d_c2r_with_plan


! --------- 1D ROUTINES
! Subroutine that perform 1D Complex to Complex 
! DFT along previously defined strategy "plan"
SUBROUTINE fast_fftw1d_with_plan(nx,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp
    INTEGER(idp), INTENT(IN) ::  nx
    COMPLEX(cpx), DIMENSION(nx), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw1d_with_plan

! Subroutine that perform 1D Real to Complex DFT
! along previously defined strategy "plan"
SUBROUTINE fast_fftw1d_r2c_with_plan(nx,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) :: nx
    REAL(num), DIMENSION(nx), INTENT(IN OUT)  :: array_in
    COMPLEX(cpx), DIMENSION(nx/2+1), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT_R2C(plans_cint(iplan), array_in, array_out)

END SUBROUTINE fast_fftw1d_r2c_with_plan

! Subroutine that perform 1D Real to Complex DFT
! along previously defined strategy "plan"
SUBROUTINE fast_fftw1d_c2r_with_plan(nx,array_in, array_out, plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE omp_lib
    USE picsar_precision, ONLY: cpx, idp, num
    INTEGER(idp), INTENT(IN) ::  nx
    COMPLEX(cpx), DIMENSION(nx/2+1), INTENT(IN OUT)  :: array_in
    REAL(num), DIMENSION(nx), INTENT(IN OUT)  :: array_out
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_EXECUTE_DFT_C2R(plans_cint(iplan), array_in, array_out)
END SUBROUTINE fast_fftw1d_c2r_with_plan

!**********************************************
!* SECTION 4: plan destruction (1D,2D,3D)
!**********************************************

! Subroutine that destroys 
!previously build plan (1D/2D/3D)
SUBROUTINE fast_fftw_destroy_plan_dft(plan)
    USE fftw3_fortran, ONLY: plans_cint
    USE iso_c_binding
    USE picsar_precision, ONLY: idp
    INTEGER(idp), DIMENSION(1), INTENT(IN OUT) :: plan
    INTEGER(idp) :: iplan

    iplan=plan(1)
    CALL DFFTW_DESTROY_PLAN(plans_cint(iplan))

END SUBROUTINE fast_fftw_destroy_plan_dft

END MODULE fastfft

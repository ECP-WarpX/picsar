! ________________________________________________________________________________________
!
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
! LOAD_BALANCING.F90
!
! MODULE FOR LOAD BALANCING 3D EM PIC SIMULATIONS
!
! Author:
! Henri Vincenti
!
! Date
! v 1.0 March 16 2016
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> This module is dedicated to the load balancing. It contains subroutines to determine
!> the load imbalance and improve the load between the MPI domains.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation March 2016
! ________________________________________________________________________________________
MODULE load_balance
  USE fields
  USE shared_data
  USE tiling
  IMPLICIT NONE

  CONTAINS

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine needs a description.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE compute_currproc_array_dimensions(nnew, ncmin, ncmax, np, mpi_rank)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN OUT) :: nnew
    INTEGER(idp), INTENT(IN) :: np, mpi_rank
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) :: ncmin, ncmax

    nnew=ncmax(mpi_rank)-ncmin(mpi_rank)+1

  END SUBROUTINE compute_currproc_array_dimensions

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine needs a description.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_1Darray_proclimits(ix1, ix2, iy1, iy2, iz1, iz2, cxmin, cymin,       &
    czmin, cxmax, cymax, czmax, npx, npy, npz, np, l_cart_comm)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: npx, npy, npz, np
    LOGICAL(lp)  :: l_cart_comm
    INTEGER(idp), INTENT(IN OUT), DIMENSION(0:np-1) :: ix1, ix2, iy1, iy2, iz1, iz2
    INTEGER(idp), INTENT(IN), DIMENSION(0:npx-1) :: cxmin, cxmax
    INTEGER(idp), INTENT(IN), DIMENSION(0:npy-1) :: cymin, cymax
    INTEGER(idp), INTENT(IN), DIMENSION(0:npz-1) :: czmin, czmax
    INTEGER(idp) :: ix, iy, iz
    INTEGER(idp) :: curr_rank

    DO iz=0, npz-1
      DO iy=0, npy-1
        DO ix=0, npx-1
          CALL pxr_convertindtoproc(comm, ix, iy, iz, npx, npy, npz, curr_rank,       &
          l_cart_comm)
          ix1(curr_rank) = cxmin(ix)
          ix2(curr_rank) = cxmax(ix)+1
          iy1(curr_rank) = cymin(iy)
          iy2(curr_rank) = cymax(iy)+1
          iz1(curr_rank) = czmin(iz)
          iz2(curr_rank) = czmax(iz)+1
        END DO
      END DO
    END DO

  END SUBROUTINE get_1Darray_proclimits

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine converts indexes of a rank in a rank number.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE pxr_convertindtoproc(mpi_comm_in, ix, iy, iz, npx, npy, npz, curr_rank,  &
    l_cart_comm)
    IMPLICIT NONE
    INTEGER(isp), INTENT(IN) :: mpi_comm_in
    INTEGER(idp), INTENT(IN) :: npx, npy, npz, ix, iy, iz
    LOGICAL(lp), INTENT(IN) :: l_cart_comm
    INTEGER(idp), INTENT(IN OUT) :: curr_rank
    INTEGER(isp) :: mpi_rank, mpi_comm_isp
    INTEGER(idp) :: ixt, iyt, izt
    mpi_comm_isp=INT(mpi_comm_in, isp)

    IF (l_cart_comm) THEN
      CALL MPI_CART_RANK(mpi_comm_isp, (/INT(iz, isp), INT(iy, isp), INT(ix, isp)/),  &
      mpi_rank, errcode)
      curr_rank=INT(mpi_rank, idp)
    ELSE
      ixt=ix
      iyt=iy
      izt=iz
      IF (ixt .LT. 0)     ixt=npx-1
      IF (ixt .GT. npx-1) ixt=0
      IF (iyt .LT. 0)     iyt=npy-1
      IF (iyt .GT. npx-1) iyt=0
      IF (izt .LT. 0)     izt=npz-1
      IF (izt .GT. npx-1) izt=0
      curr_rank= ixt +iyt*npx+izt*npx*npy
    ENDIF

  END SUBROUTINE pxr_convertindtoproc

  ! ______________________________________________________________________________________
  !> @brief
  !> Remap fields based on new split
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE mpi_remap_2D_field_component(field_new, nx_new, nz_new, field_old,       &
    nx_old, nz_old, nxg, nzg, ix1old, ix2old, iz1old, iz2old, ix1new, ix2new, iz1new,   &
    iz2new, iproc, np)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: iproc, nx_new, nz_new, nx_old, nz_old, np, nxg,nzg
    REAL(num), INTENT(IN OUT), DIMENSION(-nxg:nx_new+nxg, 1, -nzg:nz_new+nzg) ::      &
    field_new
    REAL(num), INTENT(IN), DIMENSION(-nxg:nx_old+nxg, 1, -nzg:nz_old+nzg) ::          &
    field_old
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) ::  ix1old, ix2old, iz1old, iz2old
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) ::  ix1new, ix2new, iz1new, iz2new

    INTEGER(isp) :: curr_rank, ix, iy, iz

    ! ---- MAP field__new
    CALL remap_em_2Dfields(field_old, nx_old, nz_old, ix1old, ix2old, iz1old, iz2old, &
    field_new, nx_new, nz_new, nxg, nzg, ix1new, ix2new, iz1new, iz2new, iproc, np,   &
    comm, errcode)
  END SUBROUTINE mpi_remap_2D_field_component

  ! ______________________________________________________________________________________
  !> @brief
  !> Remap fields based on new split
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE mpi_remap_3D_field_component(field_new, nx_new, ny_new, nz_new,          &
    field_old, nx_old, ny_old, nz_old, nxg, nyg, nzg, ix1old, ix2old, iy1old, iy2old,   &
    iz1old, iz2old, ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, iproc, np)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: iproc, nx_new, ny_new, nz_new, nx_old, ny_old,        &
    nz_old, np, nxg, nyg, nzg
    REAL(num), INTENT(IN OUT), DIMENSION(-nxg:nx_new+nxg, -nyg:ny_new+nyg,            &
    -nzg:nz_new+nzg) :: field_new
    REAL(num), INTENT(IN), DIMENSION(-nxg:nx_old+nxg, -nyg:ny_old+nyg,                &
    -nzg:nz_old+nzg) :: field_old
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) ::  ix1old, ix2old, iy1old, iy2old,   &
    iz1old, iz2old
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) ::  ix1new, ix2new, iy1new, iy2new,   &
    iz1new, iz2new

    INTEGER(isp) :: curr_rank, ix, iy, iz

    ! ---- MAP field__new
    CALL remap_em_3Dfields(field_old, nx_old, ny_old, nz_old, ix1old, ix2old, iy1old, &
    iy2old, iz1old, iz2old, field_new, nx_new, ny_new, nz_new, nxg, nyg, nzg, ix1new, &
    ix2new, iy1new, iy2new, iz1new, iz2new, iproc, np, comm, errcode)
  END SUBROUTINE mpi_remap_3D_field_component


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine remaps emfield_old in emfield_new and
  !> takes care of all MPI exchanges between different MPI_PROCESSES
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE remap_em_3Dfields(emfield_old, nxold, nyold, nzold, ix1old, ix2old,      &
    iy1old, iy2old, iz1old, iz2old, emfield_new, nxnew, nynew, nznew, nxg, nyg, nzg,    &
    ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, iproc, nprocs, communicator,        &
    ierrcode)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxold, nyold, nzold, nxnew, nynew, nznew, iproc,      &
    nprocs
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(isp), INTENT(IN) :: communicator
    INTEGER(isp), INTENT(IN OUT) ::  ierrcode
    REAL(num), INTENT(IN), DIMENSION(-nxg:nxold+nxg, -nyg:nyold+nyg, -nzg:nzold+nzg)  &
    :: emfield_old
    REAL(num), INTENT(IN  OUT), DIMENSION(-nxg:nxnew+nxg, -nyg:nynew+nyg,             &
    -nzg:nznew+nzg) :: emfield_new
    INTEGER(idp), INTENT(IN), DIMENSION(0:nprocs-1) :: ix1old, ix2old, iy1old,        &
    iy2old, iz1old, iz2old
    INTEGER(idp), INTENT(IN), DIMENSION(0:nprocs-1) :: ix1new, ix2new, iy1new,        &
    iy2new, iz1new, iz2new
    INTEGER(idp) :: ix, iy, iz, nsubx, nsuby, nsubz
    INTEGER(idp) :: ix3min, ix3max, iy3min, iy3max, iz3min, iz3max
    INTEGER(idp) :: ix1newip, ix2newip, iy1newip, iy2newip, iz1newip, iz2newip
    INTEGER(idp) :: ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip, iz2oldip
    INTEGER(idp) :: ixmin_old, ixmax_old, iymin_old, iymax_old, izmin_old, izmax_old
    INTEGER(idp) :: ixmin_new, ixmax_new, iymin_new, iymax_new, izmin_new, izmax_new
    LOGICAL(lp)  :: l_is_intersection
    INTEGER(isp), DIMENSION(0:nprocs-1) :: sendtype, recvtype
    INTEGER(isp), DIMENSION(0:2_isp*nprocs-1) :: requests
    INTEGER(isp) :: mpitag, proc_rank, i, nsreq, nrreq, error, count
    INTEGER(isp), PARAMETER :: nd=3
    INTEGER(isp), DIMENSION(nd) :: nsub, nglob, nglob_old, start

    mpitag=0_isp
    sendtype=0_isp
    recvtype=0_isp
    requests=0_isp

    nglob(1) = nxnew+2*nxg+1
    nglob(2) = nynew+2*nyg+1
    nglob(3) = nznew+2*nzg+1
    nglob_old(1) = nxold+2*nxg+1
    nglob_old(2) = nyold+2*nyg+1
    nglob_old(3) = nzold+2*nzg+1

    ! ------ DATA TO BE RECEIVED BY OTHER PROCS
    ! Computes intersection between new proc limit and old adjacent procs limits
    ! Recvtypes are processed in this part
    ix1newip = ix1new(iproc)
    ix2newip = ix2new(iproc)
    iy1newip = iy1new(iproc)
    iy2newip = iy2new(iproc)
    iz1newip = iz1new(iproc)
    iz2newip = iz2new(iproc)
    DO i=0, nprocs-1
      CALL get_3Dintersection(ix1newip, ix2newip, iy1newip, iy2newip, iz1newip,       &
      iz2newip, ix1old(i), ix2old(i), iy1old(i), iy2old(i), iz1old(i), iz2old(i),     &
      ix3min, ix3max, iy3min, iy3max, iz3min, iz3max, l_is_intersection)
      ! If i == iproc just do a copy of emfield_old in emfield_new
      IF ((i .EQ. iproc) .AND. l_is_intersection) THEN
        ixmin_old = ix3min - ix1old(i)  ; ixmax_old = ix3max - ix1old(i)
        iymin_old = iy3min - iy1old(i)  ; iymax_old = iy3max - iy1old(i)
        izmin_old = iz3min - iz1old(i)  ; izmax_old = iz3max - iz1old(i)
        ixmin_new = ix3min - ix1newip  ; ixmax_new = ix3max - ix1newip
        iymin_new = iy3min - iy1newip  ; iymax_new = iy3max - iy1newip
        izmin_new = iz3min - iz1newip ; izmax_new = iz3max - iz1newip
        emfield_new(ixmin_new:ixmax_new, iymin_new:iymax_new, izmin_new:izmax_new) =  &
        emfield_old(ixmin_old:ixmax_old, iymin_old:iymax_old, izmin_old:izmax_old)
        CYCLE
      END IF

      ! Found intersection area between new proc and old adjacent proc
      ! Creates RECV TYPE FOR THIS Volume
      IF (l_is_intersection .AND. (i .NE. iproc)) THEN
        !--- Create recv type
        nsub(1)  = ix3max-ix3min+1
        nsub(2)  = iy3max-iy3min+1
        nsub(3)  = iz3max-iz3min+1
        ! Arrays assumed to start at index 0 in MPI_TYPE_CREATE
        start(1) = ix3min-ix1new(iproc)+1-1+(nxg)
        start(2) = iy3min-iy1new(iproc)+1-1+(nyg)
        start(3) = iz3min-iz1new(iproc)+1-1+(nzg)
        CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob, nsub, start, MPI_ORDER_FORTRAN,      &
        MPI_DOUBLE_PRECISION, recvtype(i), ierrcode)
        ! COMMIT DATA TYPE (Really important otherwise -> MPI_ERR_TYPE or Wrong results)
        CALL MPI_TYPE_COMMIT(recvtype(i), ierrcode)
      ENDIF
    END DO

    ! ------ DATA TO BE SENT TO OTHER PROCS
    ! Computes intersection between old proc limit and new adjacent procs limits
    ! Sendtypes are processed in this part
    ix1oldip = ix1old(iproc)
    ix2oldip = ix2old(iproc)
    iy1oldip = iy1old(iproc)
    iy2oldip = iy2old(iproc)
    iz1oldip = iz1old(iproc)
    iz2oldip = iz2old(iproc)
    DO i=0, nprocs-1
      CALL get_3Dintersection(ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip,       &
      iz2oldip, ix1new(i), ix2new(i), iy1new(i), iy2new(i), iz1new(i), iz2new(i),     &
      ix3min, ix3max, iy3min, iy3max, iz3min, iz3max, l_is_intersection)
      ! Case i == iproc already treated in first loop of this subroutine
      IF (i .EQ. iproc) CYCLE

      ! Found intersection area between old proc and new adjacent procs
      ! Create sendtype for this volume
      IF (l_is_intersection .AND. (i .NE. iproc)) THEN
        !--- Create send type
        nsub(1)  = ix3max-ix3min+1
        nsub(2)  = iy3max-iy3min+1
        nsub(3)  = iz3max-iz3min+1
        ! Arrays assumed to start at index 0 in MPI_TYPE_CREATE
        start(1) = ix3min-ix1oldip+1-1+(nxg)
        start(2) = iy3min-iy1oldip+1-1+(nyg)
        start(3) = iz3min-iz1oldip+1-1+(nzg)
        CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob_old, nsub, start, MPI_ORDER_FORTRAN,  &
        MPI_DOUBLE_PRECISION, sendtype(i), ierrcode)
        ! COMMIT DATA TYPE (Really important otherwise -> MPI_ERR_TYPE or Wrong results)
        CALL MPI_TYPE_COMMIT(sendtype(i), ierrcode)
      ENDIF
    END DO


    ! POST THE IRECVs if any
    nrreq=0;
    DO i=0, nprocs-1
      IF (recvtype(i) .NE. 0) THEN
        !--- Post IRECV for this area
        CALL MPI_IRECV(emfield_new(-nxg, -nyg, -nzg), 1_isp, recvtype(i), i,          &
        MPI_ANY_TAG, communicator, requests(nrreq), ierrcode)
        nrreq=nrreq+1
      ENDIF
    END DO

    !POST THE ISENDs if any
    nsreq=0;
    DO i=0, nprocs-1
      IF (sendtype(i) .NE. 0) THEN
        !--- Post ISEND for this area
        CALL MPI_ISEND(emfield_old(-nxg, -nyg, -nzg), 1_isp, sendtype(i), i, mpitag,  &
        communicator, requests(nrreq+nsreq), ierrcode)
        nsreq=nsreq+1
      ENDIF
    END DO

    ! DO SOME SYNC BEFORE GOING ON
    count=nsreq+nrreq
    CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)

    ! FREE ALL DATATYPES
    DO i=0, nprocs-1
      IF (sendtype(i) .NE. 0) CALL MPI_TYPE_FREE(sendtype(i), ierrcode)
      IF (recvtype(i) .NE. 0) CALL MPI_TYPE_FREE(recvtype(i), ierrcode)
    END DO

  END SUBROUTINE remap_em_3Dfields


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine remaps emfield_old in emfield_new and
  !> takes care of all MPI exchanges between different MPI_PROCESSES
  !> 2D CASE X-Z
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE remap_em_2Dfields(emfield_old, nxold, nzold, ix1old, ix2old, iz1old,     &
    iz2old, emfield_new, nxnew, nznew, nxg, nzg, ix1new, ix2new, iz1new, iz2new, iproc, &
    nprocs, communicator, ierrcode)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxold, nzold, nxnew, nznew, iproc, nprocs
    INTEGER(idp), INTENT(IN) :: nxg, nzg
    INTEGER(isp), INTENT(IN) :: communicator
    INTEGER(isp), INTENT(IN OUT) ::  ierrcode
    REAL(num), INTENT(IN), DIMENSION(-nxg:nxold+nxg, 1, -nzg:nzold+nzg) ::            &
    emfield_old
    REAL(num), INTENT(IN  OUT), DIMENSION(-nxg:nxnew+nxg, 1, -nzg:nznew+nzg) ::       &
    emfield_new
    INTEGER(idp), INTENT(IN), DIMENSION(0:nprocs-1) :: ix1old, ix2old, iz1old, iz2old
    INTEGER(idp), INTENT(IN), DIMENSION(0:nprocs-1) :: ix1new, ix2new, iz1new, iz2new
    INTEGER(idp) :: ix, iy, iz, nsubx, nsubz
    INTEGER(idp) :: ix3min, ix3max, iz3min, iz3max
    INTEGER(idp) :: ix1newip, ix2newip, iz1newip, iz2newip
    INTEGER(idp) :: ix1oldip, ix2oldip, iz1oldip, iz2oldip
    INTEGER(idp) :: ixmin_old, ixmax_old, izmin_old, izmax_old
    INTEGER(idp) :: ixmin_new, ixmax_new, izmin_new, izmax_new
    LOGICAL(lp)  :: l_is_intersection
    INTEGER(isp), DIMENSION(0:nprocs-1) :: sendtype, recvtype
    INTEGER(isp), DIMENSION(0:2_isp*nprocs-1) :: requests
    INTEGER(isp) :: mpitag, proc_rank, i, nsreq, nrreq, error, count
    INTEGER(isp), PARAMETER :: nd=3
    INTEGER(isp), DIMENSION(nd) :: nsub, nglob, nglob_old, start

    mpitag=0_isp

    sendtype=0_isp
    recvtype=0_isp
    requests=0_isp

    nglob(1) = nxnew+2*nxg+1
    nglob(2) = 1
    nglob(3) = nznew+2*nzg+1
    nglob_old(1) = nxold+2*nxg+1
    nglob_old(2) = 1
    nglob_old(3) = nzold+2*nzg+1

    ! ------ DATA TO BE RECEIVED BY OTHER PROCS
    ! Computes intersection between new proc limit and old adjacent procs limits
    ! Recvtypes are processed in this part
    ix1newip = ix1new(iproc)
    ix2newip = ix2new(iproc)
    iz1newip = iz1new(iproc)
    iz2newip = iz2new(iproc)
    DO i=0, nprocs-1
      CALL get_2Dintersection(ix1newip, ix2newip, iz1newip, iz2newip, ix1old(i),      &
      ix2old(i), iz1old(i), iz2old(i), ix3min, ix3max, iz3min, iz3max,                &
      l_is_intersection)
      ! If i == iproc just do a copy of emfield_old in emfield_new
      IF ((i .EQ. iproc) .AND. l_is_intersection) THEN
        ixmin_old = ix3min - ix1old(i)  ; ixmax_old = ix3max - ix1old(i)
        izmin_old = iz3min - iz1old(i)  ; izmax_old = iz3max - iz1old(i)
        ixmin_new = ix3min - ix1newip  ; ixmax_new = ix3max - ix1newip
        izmin_new = iz3min - iz1newip ; izmax_new = iz3max - iz1newip
        emfield_new(ixmin_new:ixmax_new, 1, izmin_new:izmax_new) =                    &
        emfield_old(ixmin_old:ixmax_old, 1, izmin_old:izmax_old)
        CYCLE
      END IF

      ! Found intersection area between new proc and old adjacent proc
      ! Creates RECV TYPE FOR THIS Volume
      IF (l_is_intersection .AND. (i .NE. iproc)) THEN
        !--- Create recv type
        nsub(1)  = ix3max-ix3min+1
        nsub(2)  = 1
        nsub(3)  = iz3max-iz3min+1
        ! Arrays assumed to start at index 0 in MPI_TYPE_CREATE
        start(1) = ix3min-ix1new(iproc)+1-1+(nxg)
        start(2) = 0
        start(3) = iz3min-iz1new(iproc)+1-1+(nzg)
        CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob, nsub, start, MPI_ORDER_FORTRAN,      &
        MPI_DOUBLE_PRECISION, recvtype(i), ierrcode)
        ! COMMIT DATA TYPE (Really important otherwise -> MPI_ERR_TYPE or Wrong results)
        CALL MPI_TYPE_COMMIT(recvtype(i), ierrcode)
      ENDIF
    END DO

    ! ------ DATA TO BE SENT TO OTHER PROCS
    ! Computes intersection between old proc limit and new adjacent procs limits
    ! Sendtypes are processed in this part
    ix1oldip = ix1old(iproc)
    ix2oldip = ix2old(iproc)
    iz1oldip = iz1old(iproc)
    iz2oldip = iz2old(iproc)
    DO i=0, nprocs-1
      CALL get_2Dintersection(ix1oldip, ix2oldip, iz1oldip, iz2oldip, ix1new(i),      &
      ix2new(i), iz1new(i), iz2new(i), ix3min, ix3max, iz3min, iz3max,                &
      l_is_intersection)
      ! Case i == iproc already treated in first loop of this subroutine
      IF (i .EQ. iproc) CYCLE

      ! Found intersection area between old proc and new adjacent procs
      ! Create sendtype for this volume
      IF (l_is_intersection .AND. (i .NE. iproc)) THEN
        !--- Create send type
        nsub(1)  = ix3max-ix3min+1
        nsub(2)  = 1
        nsub(3)  = iz3max-iz3min+1
        ! Arrays assumed to start at index 0 in MPI_TYPE_CREATE
        start(1) = ix3min-ix1oldip+1-1+(nxg)
        start(2) = 0
        start(3) = iz3min-iz1oldip+1-1+(nzg)
        CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob_old, nsub, start, MPI_ORDER_FORTRAN,  &
        MPI_DOUBLE_PRECISION, sendtype(i), ierrcode)
        ! COMMIT DATA TYPE (Really important otherwise -> MPI_ERR_TYPE or Wrong results)
        CALL MPI_TYPE_COMMIT(sendtype(i), ierrcode)
      ENDIF
    END DO


    ! POST THE IRECVs if any
    nrreq=0;
    DO i=0, nprocs-1
      IF (recvtype(i) .NE. 0) THEN
        !--- Post IRECV for this area
        CALL MPI_IRECV(emfield_new(-nxg, 1, -nzg), 1_isp, recvtype(i), i,             &
        MPI_ANY_TAG, communicator, requests(nrreq), ierrcode)
        nrreq=nrreq+1
      ENDIF
    END DO

    !POST THE ISENDs if any
    nsreq=0;
    DO i=0, nprocs-1
      IF (sendtype(i) .NE. 0) THEN
        !--- Post ISEND for this area
        CALL MPI_ISEND(emfield_old(-nxg, 1, -nzg), 1_isp, sendtype(i), i, mpitag,     &
        communicator, requests(nrreq+nsreq), ierrcode)
        nsreq=nsreq+1
      ENDIF
    END DO

    ! DO SOME SYNC BEFORE GOING ON
    count=nsreq+nrreq
    CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)

    ! FREE ALL DATATYPES
    DO i=0, nprocs-1
      IF (sendtype(i) .NE. 0) CALL MPI_TYPE_FREE(sendtype(i), ierrcode)
      IF (recvtype(i) .NE. 0) CALL MPI_TYPE_FREE(recvtype(i), ierrcode)
    END DO

  END SUBROUTINE remap_em_2Dfields


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine get intersection area between two 3D domains
  !> If no intersection l_is_intersection is .FALSE.
  !> Useful to determine wether to send/recv datas bases on new CPU split
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_3Dintersection(ix1min, ix1max, iy1min, iy1max, iz1min, iz1max,       &
    ix2min, ix2max, iy2min, iy2max, iz2min, iz2max, ix3min, ix3max, iy3min, iy3max,     &
    iz3min, iz3max, l_is_intersection)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) ::  ix1min, ix1max, iy1min, iy1max, iz1min, iz1max
    INTEGER(idp), INTENT(IN) ::  ix2min, ix2max, iy2min, iy2max, iz2min, iz2max
    INTEGER(idp), INTENT(IN OUT) ::  ix3min, ix3max, iy3min, iy3max, iz3min, iz3max
    LOGICAL(lp), INTENT(IN OUT) :: l_is_intersection
    LOGICAL(lp)  :: l_is_intersectionx, l_is_intersectiony, l_is_intersectionz
    ix3min=0; iy3min=0; iz3min=0
    ix3max=0; iy3max=0; iz3max=0
    l_is_intersectionx=.FALSE.
    l_is_intersectiony=.FALSE.
    l_is_intersectionz=.FALSE.
    l_is_intersection=.FALSE.

    ! - X DIRECTION
    IF (ix2min .GE. ix1min) THEN
      IF(ix2min .LE. ix1max) THEN
        l_is_intersectionx=.TRUE.
        ix3min=ix2min
        ix3max=MIN(ix2max, ix1max)
      ENDIF
    ELSE
      IF(ix1min .LE. ix2max) THEN
        l_is_intersectionx=.TRUE.
        ix3min=ix1min
        ix3max=MIN(ix1max, ix2max)
      ENDIF
    ENDIF

    ! - Y DIRECTION
    IF (iy2min .GE. iy1min) THEN
IF(iy2min .LE. iy1max) THEN
l_is_intersectiony=.TRUE.
iy3min=iy2min
iy3max=MIN(iy2max, iy1max)
ENDIF
ELSE
IF(iy1min .LE. iy2max) THEN
l_is_intersectiony=.TRUE.
iy3min=iy1min
iy3max=MIN(iy1max, iy2max)
ENDIF
ENDIF

! - Z DIRECTION
IF (iz2min .GE. iz1min) THEN
IF(iz2min .LE. iz1max) THEN
l_is_intersectionz=.TRUE.
iz3min=iz2min
iz3max=MIN(iz2max, iz1max)
ENDIF
ELSE
IF(iz1min .LE. iz2max) THEN
l_is_intersectionz=.TRUE.
iz3min=iz1min
iz3max=MIN(iz1max, iz2max)
ENDIF
ENDIF

IF (l_is_intersectionx .AND. l_is_intersectiony .AND. l_is_intersectionz)         &
l_is_intersection=.TRUE.

END SUBROUTINE get_3Dintersection


! ______________________________________________________________________________________
!> @brief
!> This subroutine get intersection area between two 2D domains
!> If no intersection l_is_intersection is .FALSE.
!> Useful to determine wether to send/recv datas bases on new CPU split
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ______________________________________________________________________________________
SUBROUTINE get_2Dintersection(ix1min, ix1max, iz1min, iz1max, ix2min, ix2max,       &
iz2min, iz2max, ix3min, ix3max, iz3min, iz3max, l_is_intersection)
IMPLICIT NONE
INTEGER(idp), INTENT(IN) ::  ix1min, ix1max, iz1min, iz1max
INTEGER(idp), INTENT(IN) ::  ix2min, ix2max, iz2min, iz2max
INTEGER(idp), INTENT(IN OUT) ::  ix3min, ix3max, iz3min, iz3max
LOGICAL(lp), INTENT(IN OUT) :: l_is_intersection
LOGICAL(lp)  :: l_is_intersectionx, l_is_intersectionz
ix3min=0; iz3min=0
ix3max=0; iz3max=0
l_is_intersectionx=.FALSE.
l_is_intersectionz=.FALSE.
l_is_intersection=.FALSE.

! - X DIRECTION
IF (ix2min .GE. ix1min) THEN
IF(ix2min .LE. ix1max) THEN
l_is_intersectionx=.TRUE.
ix3min=ix2min
ix3max=MIN(ix2max, ix1max)
ENDIF
ELSE
IF(ix1min .LE. ix2max) THEN
l_is_intersectionx=.TRUE.
ix3min=ix1min
ix3max=MIN(ix1max, ix2max)
ENDIF
ENDIF
! - Z DIRECTION
IF (iz2min .GE. iz1min) THEN
IF(iz2min .LE. iz1max) THEN
l_is_intersectionz=.TRUE.
iz3min=iz2min
iz3max=MIN(iz2max, iz1max)
ENDIF
ELSE
IF(iz1min .LE. iz2max) THEN
l_is_intersectionz=.TRUE.
iz3min=iz1min
iz3max=MIN(iz1max, iz2max)
ENDIF
ENDIF

IF (l_is_intersectionx .AND. l_is_intersectionz) l_is_intersection=.TRUE.

END SUBROUTINE get_2Dintersection

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine get intersection area between local arrays (used in the PIC loop, 
  !> index l in the subroutine)
  !> and FFT distributed arrays (used in distributed FFTs, index g) for each MPI group. 
  !> This routine is useful to determine wether to send/recv data 
  !> between local and distributed arrays (due to their different CPU split)
  !> @author
  !> Haithem Kallala
  !> @date
  !> Creation 2018
  ! ______________________________________________________________________________________

  SUBROUTINE get2D_intersection_group_mpi
    USE fields, ONLY: nxguards, nyguards, nzguards
#if defined(FFTW)
    USE group_parameters, ONLY: cell_y_max_g, cell_y_min_g, cell_z_max_g,            &
      cell_z_min_g, g_first_cell_to_recv_y, g_first_cell_to_recv_z,                  &
      g_first_cell_to_send_y, g_first_cell_to_send_z, group_y_max_boundary,          &
      group_y_min_boundary, group_z_max_boundary, group_z_min_boundary,              &
      l_first_cell_to_recv_y, l_first_cell_to_recv_z, l_first_cell_to_send_y,        &
      l_first_cell_to_send_z, ny_group, ny_group_global_array,                       &
      nz_group_global_array, size_exchanges_g2l_recv_y, size_exchanges_g2l_recv_z,   &
      size_exchanges_g2l_send_y, size_exchanges_g2l_send_z,                          &
      size_exchanges_l2g_recv_y, size_exchanges_l2g_recv_z,                          &
      size_exchanges_l2g_send_y, size_exchanges_l2g_send_z, y_group_coords,          &
      z_group_coords
    USE iso_c_binding
#endif
    USE mpi
    USE mpi_derived_types
    USE params, ONLY: mpicom_curr
    USE picsar_precision, ONLY: idp, isp, lp
    USE shared_data, ONLY: absorbing_bcs_y, absorbing_bcs_z, c_dim, cell_y_max,      &
      cell_y_min, cell_z_max, cell_z_min, nb_group_x, nb_group_y, nb_group_z,        &
      nprocy, nprocz, ny, ny_global, nyg_group, nz, nz_global, nzg_group,            &
      p3dfft_flag, rank, x, y, y_coords, y_max_boundary, y_min_boundary, z,          &
      z_coords, z_max_boundary, z_min_boundary
    IMPLICIT NONE
    INTEGER(idp)                   :: i,j
    INTEGER(idp)                   :: iz1min,iz1max,iz2min,iz2max, iy1min,              &
    iy1max,iy2min,iy2max  
    INTEGER(isp)                   :: ierr
    INTEGER(idp)                     :: nb_proc_per_group_z,nb_proc_per_group_y
    INTEGER(idp)                     :: cell_min_group_z , cell_min_group_y
    INTEGER(idp)                     :: cell_max_group_z , cell_max_group_y
    INTEGER(idp)                     :: group_rank
    INTEGER(idp)    , ALLOCATABLE, DIMENSION(:) :: nbz_group_temp, nby_group_temp
    INTEGER(idp)       :: myself_exchange_l2g_send_z , myself_exchange_l2g_recv_z ,      &
    myself_exchange_g2l_send_z , myself_exchange_g2l_recv_z     ,                        &
    myself_exchange_l2g_send_y , myself_exchange_l2g_recv_y     ,                        &
    myself_exchange_g2l_send_y , myself_exchange_g2l_recv_y

#if defined(FFTW)

    ! -- Array allocation 
    ALLOCATE(size_exchanges_l2g_recv_z(nprocz));size_exchanges_l2g_recv_z = 0_idp
    ALLOCATE(g_first_cell_to_recv_z(nprocz));g_first_cell_to_recv_z = 0_idp
    ALLOCATE(size_exchanges_g2l_send_z(nprocz));size_exchanges_g2l_send_z = 0_idp
    ALLOCATE(g_first_cell_to_send_z(nprocz));g_first_cell_to_send_z = 0_idp

    ! -- Computes number of procs per group along y and z 
    nb_proc_per_group_z = nprocz/nb_group_z
    nb_proc_per_group_y = nprocy/nb_group_y
    ! -- Get min and max cell indices of current rank 
    iz1min = cell_z_min_g(z_coords+1) 
    iz1max = cell_z_max_g(z_coords+1)
    ALLOCATE(nbz_group_temp(nb_group_z))
    DO i=1,nb_group_z
      nbz_group_temp(i) = nz_group_global_array((i-1)*nb_group_x*nb_group_y+1)
    ENDDO
    ALLOCATE(nby_group_temp(nb_group_y))
    DO i=1,nb_group_y
      nby_group_temp(i) = ny_group_global_array((i-1)*nb_group_x+1)
    ENDDO

    ! -- Boolean variables to check if current rank is at group boundaries along Z
    cell_min_group_z = 0
    cell_min_group_z = SUM(nbz_group_temp(1:z_group_coords))
    cell_max_group_z = cell_min_group_z + nbz_group_temp(1+z_group_coords) -1
    ! -- Determine Z intersection of current rank with all other ranks along Z direction 
    DO i = 1,nprocz
      iz2min = cell_z_min(i)
      iz2max = cell_z_max(i)  
      ! Computes domain intersection along z axis to compute g_indexes
      CALL compute_send_recv_sizes_and_index_g2l_copies(iz1min, iz1max, iz2min,     &
      iz2max, size_exchanges_l2g_recv_z(i), g_first_cell_to_recv_z(i),              &
      size_exchanges_g2l_send_z(i), g_first_cell_to_send_z(i),cell_min_group_z,     &
      cell_max_group_z, nz_global,nzg_group,absorbing_bcs_z)
    ENDDO

      ALLOCATE(size_exchanges_l2g_recv_y(nprocy)); size_exchanges_l2g_recv_y = 0_idp
      ALLOCATE(g_first_cell_to_recv_y(nprocy)); g_first_cell_to_recv_y = 1_idp
      ALLOCATE(size_exchanges_g2l_send_y(nprocy)); size_exchanges_g2l_send_y = 0_idp
      ALLOCATE(g_first_cell_to_send_y(nprocy)); g_first_cell_to_send_y = 1_idp
    ! Compute f index through y axis
    IF(p3dfft_flag) THEN
      ! when using p3dfft also computes domain intersection along y axis  to compute
      ! f_indexes
      iy1min = cell_y_min_g(y_coords+1) 
      iy1max = cell_y_max_g(y_coords+1)  
      cell_min_group_y = 0
      cell_min_group_y = SUM(nby_group_temp(1:y_group_coords))
      cell_max_group_y = cell_min_group_y + nby_group_temp(1+y_group_coords) -1
      DO i = 1 , nprocy
         iy2min = cell_y_min(i)
         iy2max = cell_y_max(i) 
         CALL compute_send_recv_sizes_and_index_g2l_copies(iy1min, iy1max,         & 
         iy2min,iy2max,size_exchanges_l2g_recv_y(i),                                   &
         g_first_cell_to_recv_y(i),                                                    &
         size_exchanges_g2l_send_y(i), g_first_cell_to_send_y(i),cell_min_group_y,&
         cell_max_group_y,ny_global,nyg_group,absorbing_bcs_y)
      ENDDO
    ELSE 
       !if not using p3dfft then f_indexes are set manually to contain guardcells (in
       !this case no exchanges are done in y_direction since guardcells are copyed when
       !doing communications in z direction
       size_exchanges_l2g_recv_y(y_coords+1) = MIN(2*nyguards + ny ,ny_group)
       size_exchanges_g2l_send_y(y_coords+1) = MIN(2*nyguards + ny ,ny_group)
       g_first_cell_to_send_y(y_coords+1) = 1
       g_first_cell_to_recv_y(y_coords+1) = 1
    ENDIF

    !END OF Field_f perspective, begin field perspective

    !begin field perspective by computing indexes OF ex to exchange with ex_r
    ALLOCATE(size_exchanges_g2l_recv_z(nprocz)); size_exchanges_g2l_recv_z= 0_idp
    ALLOCATE(l_first_cell_to_recv_z(nprocz)); l_first_cell_to_recv_z = 0_idp
    ALLOCATE(size_exchanges_l2g_send_z(nprocz)); size_exchanges_l2g_send_z = 0_idp
    ALLOCATE(l_first_cell_to_send_z(nprocz));l_first_cell_to_send_z = 0_idp

    !compute r index trhoug z axis
    iz1min = cell_z_min(z_coords+1)
    iz1max = cell_z_max(z_coords+1)
    DO i=1,nprocz
      iz2min = cell_z_min_g(i)
      iz2max = cell_z_max_g(i) 
      group_rank = (i-1)/nb_proc_per_group_z
      cell_min_group_z = 0
      cell_min_group_z = SUM(nbz_group_temp(1:group_rank))
      cell_max_group_z =cell_min_group_z + nbz_group_temp(1+group_rank)-1

      !-- computes domain intersection in z direction to compute r_indexes
      CALL compute_send_recv_sizes_and_index_l2g_copies(iz1min, iz1max, iz2min,     &
      iz2max,size_exchanges_g2l_recv_z(i), l_first_cell_to_recv_z(i),                  &
      size_exchanges_l2g_send_z(i), l_first_cell_to_send_z(i),                         &
      cell_min_group_z,cell_max_group_z,nz_global,nzg_group,absorbing_bcs_z)
    ENDDO
      ALLOCATE(size_exchanges_g2l_recv_y(nprocy));size_exchanges_g2l_recv_y = 0_idp
      ALLOCATE(l_first_cell_to_recv_y(nprocy));l_first_cell_to_recv_y = 0_idp
      ALLOCATE(size_exchanges_l2g_send_y(nprocy));size_exchanges_l2g_send_y = 0_idp 
      ALLOCATE(l_first_cell_to_send_y(nprocy));l_first_cell_to_send_y = 0_idp
    IF(p3dfft_flag) THEN
      iy1min = cell_y_min(y_coords+1)
      iy1max = cell_y_max(y_coords+1)
      !i-- compute r index trhou y axis
      DO i=1,nprocy
        iy2min = cell_y_min_g(i)
        iy2max = cell_y_max_g(i) 
        group_rank = (i-1)/nb_proc_per_group_y
        cell_min_group_y = 0
        cell_min_group_y = SUM(nby_group_temp(1:group_rank))
        cell_max_group_y = cell_min_group_y + nby_group_temp(1+group_rank)-1
       !-- when using p3dfft computes domain intersection in y direction to compute
       !-- r_indexes
        CALL compute_send_recv_sizes_and_index_l2g_copies(iy1min, iy1max, iy2min,   &
        iy2max,size_exchanges_g2l_recv_y(i), l_first_cell_to_recv_y(i),                &
        size_exchanges_l2g_send_y(i), l_first_cell_to_send_y(i),                       &
        cell_min_group_y,cell_max_group_y,ny_global,nyg_group,absorbing_bcs_y)
      ENDDO
    ELSE 
      !-- when not using p3dfft r_indexes in y direction are set manually to perform
      !-- guardcells communications in y directly
      size_exchanges_g2l_recv_y(y_coords+1) =MIN(2*nyguards + ny ,ny_group)
      size_exchanges_l2g_send_y(y_coords+1) = MIN(2*nyguards + ny ,ny_group)
      l_first_cell_to_recv_y(y_coords+1) = -nyguards
      l_first_cell_to_send_y(y_coords+1) = -nyguards 
    ENDIF
    IF(nprocz==nb_group_z) THEN
    !-- corrects a bug that occurs when a group contains 1 mpi in z 
    !-- direction (since each processor intersects itself twice in this case)
     size_exchanges_l2g_send_z =0 
     size_exchanges_g2l_recv_z =0 
     size_exchanges_g2l_send_z =0 
     size_exchanges_l2g_recv_z = 0
     l_first_cell_to_recv_z =1 
     g_first_cell_to_recv_z =1 
     g_first_cell_to_send_z =1 
     l_first_cell_to_send_z =1 

     size_exchanges_l2g_send_z(z_coords+1) = nz + 2*nzguards
     size_exchanges_l2g_recv_z(z_coords+1) = nz + 2*nzguards
     size_exchanges_g2l_recv_z(z_coords+1) = nz + 2*nzguards
     size_exchanges_g2l_send_z(z_coords+1) = nz + 2*nzguards
     l_first_cell_to_recv_z(z_coords+1) = -nzguards
     g_first_cell_to_recv_z(z_coords+1) = 1
     g_first_cell_to_send_z(z_coords+1) = 1
     l_first_cell_to_send_z(z_coords+1) = -nzguards
    ENDIF
    DEALLOCATE(nbz_group_temp,nby_group_temp)
    IF(nprocy==nb_group_y) THEN
    !-- coorects a bug that occurs when a group contains 1 mpi in 
    !-- y direction (since  each processor intersects itself twice
     size_exchanges_l2g_send_y =0
     size_exchanges_g2l_recv_y =0
     size_exchanges_g2l_send_y =0
     size_exchanges_l2g_recv_y = 0
     l_first_cell_to_recv_y =1
     g_first_cell_to_recv_y =1
     g_first_cell_to_send_y =1
     l_first_cell_to_send_y =1

     size_exchanges_l2g_send_y(y_coords+1) = ny + 2*nyguards
     size_exchanges_l2g_recv_y(y_coords+1) = ny + 2*nyguards
     size_exchanges_g2l_recv_y(y_coords+1) = ny + 2*nyguards
     size_exchanges_g2l_send_y(y_coords+1) = ny + 2*nyguards
     l_first_cell_to_recv_y(y_coords+1) = -nyguards
     g_first_cell_to_recv_y(y_coords+1) = 1
     g_first_cell_to_send_y(y_coords+1) = 1
     l_first_cell_to_send_y(y_coords+1) = -nyguards
    ENDIF
    IF(c_dim == 2) THEN
      size_exchanges_l2g_send_y(y_coords+1) = 1
      size_exchanges_l2g_recv_y(y_coords+1) = 1
      size_exchanges_g2l_recv_y(y_coords+1) = 1
      size_exchanges_g2l_send_y(y_coords+1) = 1
      l_first_cell_to_recv_y(y_coords+1) = 0
      g_first_cell_to_recv_y(y_coords+1) = 1 
      g_first_cell_to_send_y(y_coords+1) = 1
      l_first_cell_to_send_y(y_coords+1) = 0
    ENDIF
    !-- Creates derived types for communications
    CALL create_derived_types_groups()
    ! -- Saves value of size_exchanges_g2l_send_z(z_coords+1) and sets it to 1 to avoid a bug in case this is 0
    myself_exchange_g2l_send_z = size_exchanges_g2l_send_z(z_coords+1)
    myself_exchange_l2g_recv_z =  size_exchanges_l2g_recv_z(z_coords+1)
    size_exchanges_g2l_send_z(z_coords+1) = 1_idp
    size_exchanges_l2g_recv_z(z_coords+1) = 1_idp
    myself_exchange_g2l_recv_z = size_exchanges_g2l_recv_z(z_coords+1)
    myself_exchange_l2g_send_z =  size_exchanges_l2g_send_z(z_coords+1)
    size_exchanges_g2l_recv_z(z_coords+1) = 1_idp
    size_exchanges_l2g_send_z(z_coords+1) = 1_idp

    myself_exchange_g2l_send_y = size_exchanges_g2l_send_y(y_coords+1)
    myself_exchange_l2g_recv_y =  size_exchanges_l2g_recv_y(y_coords+1)
    size_exchanges_g2l_send_y(y_coords+1) = 1_idp
    size_exchanges_l2g_recv_y(y_coords+1) = 1_idp
    myself_exchange_g2l_recv_y = size_exchanges_g2l_recv_y(y_coords+1)
    myself_exchange_l2g_send_y =  size_exchanges_l2g_send_y(y_coords+1)
    size_exchanges_g2l_recv_y(y_coords+1) = 1_idp
    size_exchanges_l2g_send_y(y_coords+1) = 1_idp

    !-- Compresses arrays for communications (deletes useless send recv targets and
    !-- types etc ....)
    CALL compute_effective_communication_setup()
    size_exchanges_g2l_recv_y(1) = myself_exchange_g2l_recv_y
    size_exchanges_g2l_send_y(1) = myself_exchange_g2l_send_y
    size_exchanges_l2g_recv_y(1) = myself_exchange_l2g_recv_y
    size_exchanges_l2g_send_y(1) = myself_exchange_l2g_send_y
    size_exchanges_g2l_recv_z(1) = myself_exchange_g2l_recv_z
    size_exchanges_g2l_send_z(1) = myself_exchange_g2l_send_z
    size_exchanges_l2g_recv_z(1) = myself_exchange_l2g_recv_z
    size_exchanges_l2g_send_z(1) = myself_exchange_l2g_send_z

#endif
  END SUBROUTINE get2D_intersection_group_mpi

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine cleans arrays for mpi group communications (deletes useless
  !> cells of ranks, sizes, first_indexes, types ...)  
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2018
  ! ______________________________________________________________________________________
  SUBROUTINE compute_effective_communication_setup() 
#if defined(FFTW) 
   USE group_parameters, ONLY: array_of_ranks_to_recv_from,                          &
     array_of_ranks_to_recv_from_g2l, array_of_ranks_to_recv_from_l2g,               &
     array_of_ranks_to_send_to, array_of_ranks_to_send_to_g2l,                       &
     array_of_ranks_to_send_to_l2g, g_first_cell_to_recv_y, g_first_cell_to_recv_z,  &
     g_first_cell_to_send_y, g_first_cell_to_send_z, l_first_cell_to_recv_y,         &
     l_first_cell_to_recv_z, l_first_cell_to_send_y, l_first_cell_to_send_z,         &
     mpi_ordered_comm_world, nb_comms_g2l, nb_comms_l2g, recv_type_g, recv_type_l,   &
     requests_g2l, requests_l2g, send_type_g, send_type_l, size_exchanges_g2l_recv,  &
     size_exchanges_g2l_recv_y, size_exchanges_g2l_recv_z, size_exchanges_g2l_send,  &
     size_exchanges_g2l_send_y, size_exchanges_g2l_send_z, size_exchanges_l2g_recv,  &
     size_exchanges_l2g_recv_y, size_exchanges_l2g_recv_z, size_exchanges_l2g_send,  &
     size_exchanges_l2g_send_y, size_exchanges_l2g_send_z, work_array_g2l,           &
     work_array_l2g
#endif
#if defined(FFTW) 
   USE params, ONLY: mpicom_curr
   USE picsar_precision, ONLY: idp, isp
#endif
#if defined(FFTW) 
   INTEGER(idp)     :: ll,ii,i,j,k,n,jj,jjj,kk,kkk,shift_y,shift_z
   INTEGER(idp)  , ALLOCATABLE, DIMENSION(:) :: temp1,temp2,temp3,temp4,&
        temp5,temp6,temp7,temp8
   INTEGER(isp)  , ALLOCATABLE, DIMENSION(:) :: temp_rs, temp_sr
   INTEGER(isp) ,  ALLOCATABLE, DIMENSION(:,:,:):: topo_array
   INTEGER(isp)                                 :: ordered_rank
  !-- begin compute array of target and recv procs   
  !-- array_of_ranks_to_send_to and array_of_ranks_to_recv_from respectively
  !-- these arrays are 1d-nprocy*nprocz sizes 
  !-- But you can think of these arrays as 2d arrays (nprocy,nprocz)
  
  !-- array_of_ranks_to_send_to(i,j) = 
  !-- rank_of_mpi_with: y_coords_target =  y_coords + (j-1) , 
  !-- z_coords_target = z_coords + (i - 1) 
  !-- + periodic bcs on procs
  
  !-- array_of_ranks_to_recv_from(i,j) = 
  !-- rank_of_mpi_with: y_coords_target =  y_coords - (j-1) , 
  !-- z_coords_target =  z_coords - (i - 1) 
  !-- +periodic bcs on procs

  ! -- Using mpi_ordered_comm_world : computes arrays to send to and arrays of
  ! -- ranks to recieve from depending on current mpi z_coords and y_coords   
  ! -- mpi_ordered_comm_world is a cartesian grid like communicator with
  ! -- grid coordinates (x_coords,y_coords,z_coords) 
  ! -- Works for both cartesian and random topologies
  
  !-- How it works: 
  !-- The current rank communicates along y axis (with a fixed z coord
  !-- target/source) 
  !-- and after all communications are performed along y axis only increments z
  !-- coords target/source
  !-- loop over nprocz, loop over nprocy (first loop over all y and then
  !--  increment z)
  !-- Sends are posted to the i-th  right neighbor  and at the same time 
  !-- Recieves are posted to i-th left neighbor
  ordered_rank = INT(rank,isp)
  ALLOCATE(topo_array(0:nprocx-1,0:nprocy-1,0:nprocz-1))
  CALL MPI_ALLGATHER(ordered_rank, 1_isp, MPI_INTEGER,topo_array,1_isp,MPI_INTEGER,  &
       mpi_ordered_comm_world,  errcode)

  ALLOCATE(array_of_ranks_to_send_to(nprocy*nprocz)) 
  DO i=0,nprocy-1
    DO j=0,nprocz-1
      shift_y = MODULO(y_coords + i,nprocy) 
      shift_z = MODULO(z_coords + j,nprocz)
      array_of_ranks_to_send_to(i+j*nprocy+1) = topo_array(x_coords,shift_y,shift_z)
    ENDDO
  ENDDO

  ALLOCATE(array_of_ranks_to_recv_from(nprocy*nprocz))
  DO i=0,nprocy-1
    DO j=0,nprocz-1
      shift_y = MODULO(y_coords - i,nprocy) 
      shift_z = MODULO(z_coords - j,nprocz)
      array_of_ranks_to_recv_from(i+j*nprocy+1) = topo_array(x_coords,shift_y,shift_z)
    ENDDO
  ENDDO


  DEALLOCATE(topo_array)

  !-- end compute array of target and recv procs 

  !-- begin extend local and global indexes and sizes
  !-- To get an extended (nprocy*nprocz) array for sizes and indexes
  !-- To be compressed later with the same protocole
  ALLOCATE(size_exchanges_g2l_send(nprocy,nprocz), &
      size_exchanges_g2l_recv(nprocy,nprocz)    ) 
  ALLOCATE(size_exchanges_l2g_send(nprocy,nprocz), &
      size_exchanges_l2g_recv(nprocy,nprocz)    )
  !-- size_exchanges_g2l_send / size_exchanges_g2l_recv  /size_exchanges_l2g_send
  !-- size_exchanges_l2g_recv 2d arrays  of (nprocy,nprocz) sizes  
  !-- size_exchanges_g2l_send(i,j) =  
  !-- size_exchanges_g2l_send_y(i)*size_exchanges_g2l_send_z(j) 
  !-- so basically the number of double precision scalars to exchange between
  !-- current proc and proc with y_coor=i-1,z_coord=j-1 (and same x_coord)

  DO i=1,nprocy
    DO j=1,nprocz
      size_exchanges_g2l_send(i,j) = size_exchanges_g2l_send_y(i)* &
      size_exchanges_g2l_send_z(j)
      size_exchanges_l2g_send(i,j) = size_exchanges_l2g_send_y(i)* &
      size_exchanges_l2g_send_z(j)
      size_exchanges_l2g_recv(i,j) = size_exchanges_l2g_recv_y(i)* &
      size_exchanges_l2g_recv_z(j)
      size_exchanges_g2l_recv(i,j) = size_exchanges_g2l_recv_y(i)* &
      size_exchanges_g2l_recv_z(j)
    ENDDO
  ENDDO
  ALLOCATE(temp1(nprocz),temp2(nprocz),temp3(nprocz),temp4(nprocz))
  ALLOCATE(temp5(nprocz),temp6(nprocz),temp7(nprocz),temp8(nprocz))
  !-- transfrom size_exchanges_l2g/52l_recv/send_y/z from arrays of size
  !-- nprocy or nprocz to arrays of size nprocy*nprocz
  !-- you can still think of these new arrays as 2d arrays of size (nprocy,nprocz)
  !-- same thing is done for r/s_first_cell_to_r/s_y/z
  temp1=g_first_cell_to_send_z
  temp2=l_first_cell_to_send_z
  temp3=g_first_cell_to_recv_z
  temp4=l_first_cell_to_recv_z
  temp5=size_exchanges_g2l_send_z
  temp6=size_exchanges_l2g_send_z
  temp7=size_exchanges_l2g_recv_z 
  temp8=size_exchanges_g2l_recv_z

  DEALLOCATE(g_first_cell_to_send_z,l_first_cell_to_send_z,&
  g_first_cell_to_recv_z,l_first_cell_to_recv_z)

  DEALLOCATE(size_exchanges_g2l_send_z,size_exchanges_l2g_send_z,&
  size_exchanges_l2g_recv_z,size_exchanges_g2l_recv_z)

  ALLOCATE(g_first_cell_to_send_z(nprocy*nprocz));   g_first_cell_to_send_z=0
  ALLOCATE(g_first_cell_to_recv_z(nprocy*nprocz));   g_first_cell_to_recv_z=0
  ALLOCATE(l_first_cell_to_send_z(nprocy*nprocz));   l_first_cell_to_send_z=0
  ALLOCATE(l_first_cell_to_recv_z(nprocy*nprocz));   l_first_cell_to_recv_z=0
  ALLOCATE(size_exchanges_g2l_send_z(nprocy*nprocz));   size_exchanges_g2l_send_z =0
  ALLOCATE(size_exchanges_l2g_send_z(nprocy*nprocz));   size_exchanges_l2g_send_z=0
  ALLOCATE(size_exchanges_l2g_recv_z(nprocy*nprocz));   size_exchanges_l2g_recv_z=0
  ALLOCATE(size_exchanges_g2l_recv_z(nprocy*nprocz));   size_exchanges_g2l_recv_z=0
  

  DO i=1,nprocy
    DO j=1,nprocz
      g_first_cell_to_send_z(i+(j-1)*nprocy) =        temp1(j)
      l_first_cell_to_send_z(i+(j-1)*nprocy) =        temp2(j)
      g_first_cell_to_recv_z(i+(j-1)*nprocy) =        temp3(j)
      l_first_cell_to_recv_z(i+(j-1)*nprocy) =        temp4(j)
      size_exchanges_g2l_send_z(i+(j-1)*nprocy) = temp5(j)
      size_exchanges_l2g_send_z(i+(j-1)*nprocy) = temp6(j)
      size_exchanges_l2g_recv_z(i+(j-1)*nprocy) = temp7(j)
      size_exchanges_g2l_recv_z(i+(j-1)*nprocy) = temp8(j)
    ENDDO
  ENDDO
  DEALLOCATE(temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8)

  ALLOCATE(temp1(nprocy),temp2(nprocy),temp3(nprocy),temp4(nprocy))
  ALLOCATE(temp5(nprocy),temp6(nprocy),temp7(nprocy),temp8(nprocy))

  temp1=g_first_cell_to_send_y
  temp2=l_first_cell_to_send_y
  temp3=g_first_cell_to_recv_y
  temp4=l_first_cell_to_recv_y
  temp5=size_exchanges_g2l_send_y
  temp6=size_exchanges_l2g_send_y
  temp7=size_exchanges_l2g_recv_y
  temp8=size_exchanges_g2l_recv_y

  DEALLOCATE(g_first_cell_to_send_y,l_first_cell_to_send_y,      &
  g_first_cell_to_recv_y,l_first_cell_to_recv_y)
  DEALLOCATE(size_exchanges_g2l_send_y,size_exchanges_l2g_send_y,&
  size_exchanges_l2g_recv_y,size_exchanges_g2l_recv_y)


  ALLOCATE(g_first_cell_to_send_y(nprocy*nprocz)); g_first_cell_to_send_y=0
  ALLOCATE(g_first_cell_to_recv_y(nprocy*nprocz)); g_first_cell_to_recv_y=0
  ALLOCATE(l_first_cell_to_send_y(nprocy*nprocz)); l_first_cell_to_send_y=0
  ALLOCATE(l_first_cell_to_recv_y(nprocy*nprocz)); l_first_cell_to_recv_y=0
  ALLOCATE(size_exchanges_g2l_send_y(nprocy*nprocz)); size_exchanges_g2l_send_y=0
  ALLOCATE(size_exchanges_l2g_send_y(nprocy*nprocz)); size_exchanges_l2g_send_y=0
  ALLOCATE(size_exchanges_l2g_recv_y(nprocy*nprocz)); size_exchanges_l2g_recv_y=0
  ALLOCATE(size_exchanges_g2l_recv_y(nprocy*nprocz)); size_exchanges_g2l_recv_y=0

  DO i=1,nprocy
   DO j=1,nprocz
     g_first_cell_to_send_y(i+(j-1)*nprocy) =        temp1(i)
     l_first_cell_to_send_y(i+(j-1)*nprocy) =        temp2(i)
     g_first_cell_to_recv_y(i+(j-1)*nprocy) =        temp3(i)
     l_first_cell_to_recv_y(i+(j-1)*nprocy) =        temp4(i)
     size_exchanges_g2l_send_y(i+(j-1)*nprocy) = temp5(i)
     size_exchanges_l2g_send_y(i+(j-1)*nprocy) = temp6(i)
     size_exchanges_l2g_recv_y(i+(j-1)*nprocy) = temp7(i)
     size_exchanges_g2l_recv_y(i+(j-1)*nprocy) = temp8(i)
    ENDDO
  ENDDO
  DEALLOCATE(temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8)

  !-- END of extending r_arraysy/z and f_arrays_y/z
  
  
  !-- if non blocking comms are used then computes then number of called mpi_isend
  !-- mpi_irecv for each field component exchanged 
  !-- this is important to allocate a suitable mpi_requests array 
  IF(mpicom_curr .EQ. 0) THEN
    n=0
    DO i=1,nprocy
      DO j=1,nprocz
        IF(i ==y_coords+1 .AND. j==z_coords+1) CYCLE
        IF(size_exchanges_l2g_send(i,j) .GT. 0) n = n + 1
        IF(size_exchanges_l2g_recv(i,j) .GT. 0) n = n +1
      ENDDO
    ENDDO
    ! mpi request array for l->g communication
    ALLOCATE(requests_l2g(n))
    n=0 
    DO i=1,nprocy
      DO j=1,nprocz
        IF(i==y_coords+1 .AND. j==z_coords+1) CYCLE
        IF(size_exchanges_g2l_send(i,j) .GT. 0) n = n+1
        IF(size_exchanges_g2l_recv(i,j) .GT. 0) n = n+1
      ENDDO
    ENDDO
    ! mpi request array for g->l
    ALLOCATE(requests_g2l(n))
  ENDIF
  
  
  !-- end non blocking comms conditioning
  
  !-- Begin computing work_array_g2l and work_array_ arrays
  !-- first we browse all size_exchanges_g2l_send and size_exchanges_g2l_recv and 
  !-- its corresponding
  !-- then we determine if there is any communication done between current rank and
  !-- array_of_rank_to_s/r(i,j) 
  !-- if  yes then n+=1 ;  work_array_(n) = (i-1)+(j-1)*nprocy (sort of compressing
  !-- indexes of array_of_ranks) 
  !-- work_array(n) is a sort of a relative compressed 2-d distance for current proc
  !-- regarding the procs with which comms are done (wether send or recieve)
  !-- During the same iteration
  !-- Knowing the relative distance we can then determine back 
  !-- the ranks of procs to send to and the ranks of procs to recieve from
  !-- eventually (Assuming the comms protocole stated before)
  n=0
  DO j=1,nprocz
    DO i=1,nprocy

      ii =MODULO(y_coords+i-1,nprocy)+1; kk=MODULO(y_coords-(i-1),nprocy) +1
      jj = MODULO(z_coords+j-1,nprocz) +1; ll =MODULO(z_coords-(j-1),nprocz)+1
      IF(size_exchanges_l2g_send(ii,jj) .GT. 0     &
        .OR. size_exchanges_l2g_recv(kk,ll) .GT. 0)  n=n+1
    ENDDO
  ENDDO
  ALLOCATE(work_array_l2g(n))
  work_array_l2g=0
  n=0
  DO j=1,nprocz
    DO i=1,nprocy
      ii =MODULO(y_coords+i-1,nprocy)+1; kk= MODULO(y_coords-(i-1),nprocy) +1
      jj = MODULO(z_coords+j-1,nprocz) +1; ll = MODULO(z_coords-(j-1),nprocz)+1
      IF(size_exchanges_l2g_send(ii,jj) .GT. 0     &
        .OR. size_exchanges_l2g_recv(kk,ll) .GT. 0) THEN
        n=n+1
        work_array_l2g(n)=(i-1)+(j-1)*nprocy
      ENDIF
    ENDDO
  ENDDO

  n=0

  DO j=1,nprocz
    DO i=1,nprocy
      ii =MODULO(y_coords+i-1,nprocy)+1; kk= MODULO(y_coords-(i-1),nprocy) +1
      jj = MODULO(z_coords+j-1,nprocz) +1; ll = MODULO(z_coords-(j-1),nprocz) +1
      IF(size_exchanges_g2l_send(ii,jj) .GT. 0     &
        .OR. size_exchanges_g2l_recv(kk,ll) .GT. 0)  n=n+1
    ENDDO
  ENDDO
  ALLOCATE(work_array_g2l(n))
  work_array_g2l=0
  n=0

  DO j=1,nprocz
    DO i=1,nprocy
      ii =MODULO(y_coords+i-1,nprocy)+1; kk= MODULO(y_coords-(i-1),nprocy) +1
      jj =MODULO(z_coords+j-1,nprocz) +1; ll = MODULO(z_coords-(j-1),nprocz) +1
      IF(size_exchanges_g2l_send(ii,jj) .GT. 0     &
        .OR. size_exchanges_g2l_recv(kk,ll) .GT. 0) THEN  
        n=n+1
        work_array_g2l(n)=(i-1)+(j-1)*nprocy
      ENDIF
  ENDDO
  ENDDO

  nb_comms_l2g = SIZE(work_array_l2g)
  nb_comms_g2l = SIZE(work_array_g2l)

  !--end work_array computing section

  ALLOCATE(temp1(nprocy*nprocz)); temp1 = size_exchanges_g2l_send_z
  ALLOCATE(temp2(nprocy*nprocz)); temp2 = size_exchanges_g2l_recv_z
  ALLOCATE(temp3(nprocy*nprocz)); temp3 = size_exchanges_g2l_send_y
  ALLOCATE(temp4(nprocy*nprocz)); temp4 = size_exchanges_g2l_recv_y
  ALLOCATE(temp5(nprocy*nprocz)); temp5 = g_first_cell_to_send_z
  ALLOCATE(temp6(nprocy*nprocz)); temp6 = l_first_cell_to_recv_z
  ALLOCATE(temp7(nprocy*nprocz)); temp7 = g_first_cell_to_send_y
  ALLOCATE(temp8(nprocy*nprocz)); temp8 = l_first_cell_to_recv_y

  DEALLOCATE(size_exchanges_g2l_send_z,size_exchanges_g2l_recv_z,&
  size_exchanges_g2l_send_y ,size_exchanges_g2l_recv_y, g_first_cell_to_send_z,&
  l_first_cell_to_recv_z, g_first_cell_to_send_y,l_first_cell_to_recv_y) 

  ALLOCATE(size_exchanges_g2l_send_z(nb_comms_g2l))
  ALLOCATE(size_exchanges_g2l_recv_z(nb_comms_g2l))
  ALLOCATE(size_exchanges_g2l_send_y(nb_comms_g2l))
  ALLOCATE(size_exchanges_g2l_recv_y(nb_comms_g2l))
  ALLOCATE(g_first_cell_to_send_z(nb_comms_g2l))
  ALLOCATE(l_first_cell_to_recv_z(nb_comms_g2l))
  ALLOCATE(g_first_cell_to_send_y(nb_comms_g2l))
  ALLOCATE(l_first_cell_to_recv_y(nb_comms_g2l))
  ALLOCATE(array_of_ranks_to_send_to_g2l(nb_comms_g2l))
  ALLOCATE(array_of_ranks_to_recv_from_g2l(nb_comms_g2l))
  ALLOCATE(temp_rs(nprocy*nprocz));ALLOCATE(temp_sr(nprocy*nprocz))
  temp_rs = send_type_g;   temp_sr = recv_type_l
  DEALLOCATE(send_type_g); ALLOCATE(send_type_g(nb_comms_g2l))
  DEALLOCATE(recv_type_l); ALLOCATE(recv_type_l(nb_comms_g2l))


  !-- compressing arrays global to local
  !-- saving only  relevent indexes and sizes ranks and types
  DO ii=1,nb_comms_g2l
     i = work_array_g2l(ii) 
     k =   MODULO(i,nprocy)  !Y
     j =   i/nprocy          !Z
     jj  =  MODULO(z_coords + j,nprocz) 
     jjj = MODULO(z_coords -j,nprocz) 
     kk  =  MODULO(y_coords +  k,nprocy) 
     kkk =  MODULO(y_coords-k,nprocy)  
     size_exchanges_g2l_send_z(ii) = temp1(jj*nprocy+kk+1)
     size_exchanges_g2l_recv_z(ii) = temp2(jjj*nprocy+kkk+1)
     size_exchanges_g2l_send_y(ii) = temp3(jj*nprocy+kk+1)
     size_exchanges_g2l_recv_y(ii) = temp4(jjj*nprocy+kkk+1)
     g_first_cell_to_send_z(ii) = temp5(jj*nprocy+kk+1)
     l_first_cell_to_recv_z(ii) =temp6(jjj*nprocy+kkk+1)
     g_first_cell_to_send_y(ii) = temp7(jj*nprocy+kk+1)
     l_first_cell_to_recv_y(ii) = temp8(jjj*nprocy+kkk+1)
     array_of_ranks_to_send_to_g2l(ii) = array_of_ranks_to_send_to(k+j*nprocy+1)
     array_of_ranks_to_recv_from_g2l(ii) = array_of_ranks_to_recv_from(k+j*nprocy+1)
     send_type_g(ii) = temp_rs(jj*nprocy+kk+1)
     recv_type_l(ii) = temp_sr(jjj*nprocy+kkk+1)
  ENDDO

  !-- end
  
  !-- compress arrays real to fourier

  temp1 = size_exchanges_l2g_send_z
  temp2 = size_exchanges_l2g_recv_z
  temp3 = size_exchanges_l2g_send_y
  temp4 = size_exchanges_l2g_recv_y
  temp5 = l_first_cell_to_send_z
  temp6 = g_first_cell_to_recv_z
  temp7 = l_first_cell_to_send_y
  temp8 = g_first_cell_to_recv_y
  
  DEALLOCATE(size_exchanges_l2g_send_z,size_exchanges_l2g_recv_z,&
  size_exchanges_l2g_send_y ,size_exchanges_l2g_recv_y,l_first_cell_to_send_z,&
  g_first_cell_to_recv_z, l_first_cell_to_send_y,g_first_cell_to_recv_y)

  ALLOCATE(size_exchanges_l2g_send_z(nb_comms_l2g))
  ALLOCATE(size_exchanges_l2g_recv_z(nb_comms_l2g))
  ALLOCATE(size_exchanges_l2g_send_y(nb_comms_l2g))
  ALLOCATE(size_exchanges_l2g_recv_y(nb_comms_l2g))
  ALLOCATE(l_first_cell_to_send_z(nb_comms_l2g))
  ALLOCATE(g_first_cell_to_recv_z(nb_comms_l2g))
  ALLOCATE(l_first_cell_to_send_y(nb_comms_l2g))
  ALLOCATE(g_first_cell_to_recv_y(nb_comms_l2g))
  ALLOCATE(array_of_ranks_to_send_to_l2g(nb_comms_l2g))
  ALLOCATE(array_of_ranks_to_recv_from_l2g(nb_comms_l2g))
  temp_rs = send_type_l;   temp_sr = recv_type_g
  DEALLOCATE(send_type_l); ALLOCATE(send_type_l(nb_comms_l2g))
  DEALLOCATE(recv_type_g); ALLOCATE(recv_type_g(nb_comms_l2g))



  !-- compressing arrays local to global
  !-- saving onlt  relevent indexes and sizes ranks and types

  DO ii=1,nb_comms_l2g
     i = work_array_l2g(ii)
     k =   MODULO(i,nprocy)  !Y
     j =   i/nprocy          !Z
     jj  =  MODULO(z_coords + j,nprocz)
     jjj =  MODULO(z_coords -j,nprocz)
     kk  =  MODULO(y_coords +  k,nprocy)
     kkk =  MODULO(y_coords-k,nprocy) 
     size_exchanges_l2g_send_z(ii) = temp1(jj*nprocy + kk+1)
     size_exchanges_l2g_recv_z(ii) = temp2(jjj*nprocy +kkk+1)
     size_exchanges_l2g_send_y(ii) = temp3(jj*nprocy+  kk +1 )
     size_exchanges_l2g_recv_y(ii) = temp4( jjj*nprocy +kkk+1)
     l_first_cell_to_send_z(ii) = temp5(jj*nprocy +kk+1)
     g_first_cell_to_recv_z(ii) =temp6(jjj*nprocy +kkk+1)
     l_first_cell_to_send_y(ii) = temp7(jj*nprocy+kk+1)
     g_first_cell_to_recv_y(ii) = temp8(jjj*nprocy+kkk+1)
     array_of_ranks_to_send_to_l2g(ii) = array_of_ranks_to_send_to(k+j*nprocy+1)
     array_of_ranks_to_recv_from_l2g(ii) = array_of_ranks_to_recv_from(k+j*nprocy+1)
    send_type_l(ii) = temp_rs(jj*nprocy+kk+1)
    recv_type_g(ii) = temp_sr(jjj*nprocy+kkk+1)

  ENDDO
  DEALLOCATE(temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8)
  !--end

  DEALLOCATE(temp_rs,temp_sr)
  !-- end

#endif
  END SUBROUTINE compute_effective_communication_setup


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine creates mpi derived types for group comms when using pdfft
  !> @author
  !> Haithem Kallala
  !
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________

  SUBROUTINE create_derived_types_groups
  USE constants, ONLY: c_ndims
  USE fields, ONLY: nxguards, nyguards, nzguards
#if defined(FFTW)
  USE group_parameters, ONLY: nx_group, recv_type_g, recv_type_l, send_type_g,       &
    send_type_l, size_exchanges_g2l_recv_y, size_exchanges_g2l_recv_z,               &
    size_exchanges_g2l_send_y, size_exchanges_g2l_send_z, size_exchanges_l2g_recv_y, &
    size_exchanges_l2g_recv_z, size_exchanges_l2g_send_y, size_exchanges_l2g_send_z
  USE iso_c_binding
#endif
  USE mpi
  USE mpi_derived_types
#if defined(FFTW)
  USE mpi_fftw3, ONLY: local_nx, local_ny, local_nz
#endif
  USE mpi_type_constants, ONLY: mpidbl
  USE picsar_precision, ONLY: idp, isp
  USE shared_data, ONLY: nprocy, nprocz, nx, ny, nz, y, z
  INTEGER(idp)         ::  i,j
  INTEGER(idp), DIMENSION(c_ndims) :: sizes, subsizes, starts
  INTEGER(isp)                     :: basetype

#if defined(FFTW)
   basetype = mpidbl

   ALLOCATE(send_type_g(nprocz*nprocy),recv_type_g(nprocz*nprocy))

   DO i = 1,nprocz
     DO j=1,nprocy
       sizes(1) = local_nx
       sizes(2) = local_ny
       sizes(3) = local_nz

!creates mpi_type_to_recv_in_f_field from current mpi to the mpi with rank
! is determined by y_coord_target = j-1, z_coords_target = i-1 , x_coords_target
! = x_coords_curent_rank 

       subsizes(1) = MIN(2*nxguards + nx ,local_nx)
       subsizes(2) = size_exchanges_l2g_recv_y(j)
       subsizes(3) = size_exchanges_l2g_recv_z(i)
       starts = 1
       recv_type_g((i-1)*nprocy+j) = create_3d_array_derived_type(basetype, subsizes,sizes,starts)

!creates mpi_type_to_send_from_f_field from current mpi to the mpi with rank
! is determined by y_coord_target = j-1, z_coords_target = i-1 , x_coords_target
! = x_coords_curent_rank 

       subsizes(1) = MIN(2*nxguards + nx ,nx_group)
       subsizes(2) = size_exchanges_g2l_send_y(j)
       subsizes(3) = size_exchanges_g2l_send_z(i)
       send_type_g((i-1)*nprocy+j) = create_3d_array_derived_type(basetype,subsizes,sizes,starts)
     ENDDO
   ENDDO
   ALLOCATE(send_type_l(nprocz*nprocy),recv_type_l(nprocz*nprocy))
   DO i = 1,nprocz
     Do j = 1,nprocy
       sizes(1) = 2*nxguards + nx + 1
       sizes(2) = 2*nyguards + ny + 1
       sizes(3) = 2*nzguards + nz + 1

!creates mpi_type_to_send_from_r_field from current mpi to the mpi with rank
! is determined by y_coord_target = j-1, z_coords_target = i-1 , x_coords_target
! = x_coords_curent_rank 

       subsizes(1) = MIN(2*nxguards + nx ,local_nx)
       subsizes(2) = size_exchanges_g2l_recv_y(j)
       subsizes(3) = size_exchanges_g2l_recv_z(i)
       starts = 1
       recv_type_l((i-1)*nprocy+j) = create_3d_array_derived_type(basetype,subsizes,sizes,starts)

!creates mpi_type_to_send_from_r_field from current mpi to the mpi with rank
! is determined by y_coord_target = j-1, z_coords_target = i-1 , x_coords_target
! = x_coords_curent_rank 

       subsizes(1) = MIN(2*nxguards + nx ,nx_group)
       subsizes(2) = size_exchanges_l2g_send_y(j)
       subsizes(3) = size_exchanges_l2g_send_z(i)
       send_type_l((i-1)*nprocy+j) =create_3d_array_derived_type(basetype,subsizes,sizes,starts)
     ENDDO
   ENDDO
#endif
  END SUBROUTINE create_derived_types_groups

  ! ______________________________________________________________________________________
  !> @brief
  !> This routine computes for current rank, the starting indices and sizes (in the local 
  !> field arrays) of data to be sent to / received from FFT arrays of another rank. 
  !> It thus basically treats send/rev exchanges required for copies from/to the local
  !> field arrays. This is useful for the copy operation from the local PIC grid 
  !> to the FFT arrays before the global FFT operation amongst a given MPI group. 
  !> This routine computes along a given axis (X,Y or Z): 
  !> (i) intersection of the local field array (without guard cells) of current rank
  !> with the FFT field array (without group guard cells) from another rank. 
  !> Returns size and first cell of data to be received from the FFT array
  !> to the local arrays. 
  !> (ii) intersection of the local field array  (without guard cells) of current rank
  !> with the FFT field array (with group guard cells) from another rank.  
  !> Returns size and first cell of data to be sent to FFT array from the local array
  !> @param[in] iz1min INTEGER(idp) - min index of current rank along axis (X,Y or Z)
  !> @param[in] iz1max INTEGER(idp) - max index of current rank along axis (X,Y or Z)
  !> @param[in] iz2min INTEGER(idp) - min index of other rank along axis (X,Y or Z)
  !> @param[in] iz2max INTEGER8(idp) - max index of other rank along axis (X,Y or Z)
  !> @param[in] cell_g_min INTEGER(idp) - index of first cell of group (without
  !> guards cells)
  !> @param[in] cell_g_max INTEGER(idp) - index of last cell of group (without
  !> guardcells)
  !> @param[in,out] size_to_exchange_recv INTEGER(idp) - size of data to be received 
  !> from the FFT array of other rank 
  !> @param[in,out] first_cell_recv INTEGER(idp) - position in local array where data 
  !> is to be received (from the FFT array of other rank)
  !> @param[in,out] size_to_exchange_send INTEGER(idp) - size of data to be sent to 
  !> other rank 
  !> @param[in,out] first_cell_send INTEGER(idp) - starting position in local array of 
  !> data to be sent to other rank 
  !> @param[in] n_global INTEGER(idp) - number of cells of the global grid along X,Y or Z
  !> @param[in] n_guards INTEGER(idp) - number of group guard cells along X, Y or Z 
  !> @author
  !> Haithem Kallala 
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
  SUBROUTINE compute_send_recv_sizes_and_index_l2g_copies                     &
        (iz1min ,iz1max ,iz2min ,iz2max,                                      &
        size_to_exchange_recv, first_cell_recv,                               &
        size_to_exchange_send,first_cell_send,cell_g_min,cell_g_max,n_global, &
        n_guards,is_open_bc)
    USE picsar_precision, ONLY: idp, lp
    INTEGER(idp) , INTENT(IN)  :: iz1min, iz1max, iz2min, iz2max
    INTEGER(idp) , INTENT(INOUT)  ::   size_to_exchange_recv, first_cell_recv
    INTEGER(idp) , INTENT(INOUT)  ::   size_to_exchange_send, first_cell_send
    INTEGER(idp)                  :: index_rf, index_ff, index_rl, index_fl
    INTEGER(idp)                  :: select_case
    INTEGER(idp) , INTENT(IN)     :: cell_g_max,cell_g_min
    INTEGER(idp) , INTENT(IN)     :: n_global, n_guards
    LOGIcAL(lp)  , INTENT(IN)     :: is_open_bc

#if defined(FFTW)    
    index_ff = MAX(iz2min, cell_g_min)
    index_fl = MIN(iz2max, cell_g_max)
    index_rf = iz1min
    index_rl = iz1max

    ! -- Computes sizes to be recieved by local fields during global to local
    ! -- process
    IF(index_fl .GE. index_ff) THEN
      size_to_exchange_recv = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + & 
      1_idp, 0_idp)
      first_cell_recv = MAX(index_rf,index_ff) - index_rf 
    ELSE
      size_to_exchange_recv = 0_idp
    ENDIF
    ! -- first cell recv need to be set to a value contained in lower_b and
    ! -- upper_b of local field even if no exchanges are performed
    IF(size_to_exchange_recv==0_idp) first_cell_recv=0_idp

    ! -- Computes sizes to be sent by local field during local to global process
    index_ff = iz2min
    index_fl = iz2max
    size_to_exchange_send = 0_idp 
    !-- set a flag for different possible configurations for global field limits
    !-- If global field laps a physical domain boundary, max cell index and min
    !-- cell index need to be recalculated accordingly
    IF(iz2min .GE. 0_idp .AND. iz2max .GE. 0_idp .AND. iz2min .LT. n_global & 
    .AND. iz2max .LT. n_global)  THEN
       select_case = 0_idp  ! trivial case : global field is in real domain
    ELSE IF((iz2min .LT. 0_idp .AND. iz2max .LT. 0_idp) .OR.                &
    (iz2min .GE. n_global .AND. iz2max .GE. n_global)) THEN
       select_case = 1      ! global field is contained in a physical domain
                            ! boundary

    ELSE IF (iz2min .LT. 0_idp .AND. iz2max .GE. 0_idp) THEN
       select_case = 2      ! global field begins in physical domain boundary
                            ! and ends in real domain

    ELSE IF (iz2min .LT. n_global .AND. iz2max .GE. n_global) THEN
       select_case = 3      ! global field begins real domain and ends 
                            ! in physical domain boundary region
   ENDIF

   IF (select_case == 0) THEN
      ! -- Trivial case : global grid field is inside the domain 
      ! -- intersection is computed normally
      ! -- No boundary guard cells are contained in global field 
      index_ff = iz2min
      index_fl = iz2max
      size_to_exchange_send = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + &
      1_idp, 0_idp)
      first_cell_send = MAX(index_rf,index_ff) - index_rf
   ELSE IF (select_case ==  1) THEN
      ! -- The full global field is contained in a boundary ghost region 
      ! -- the limits of global fields are flipped to the corresponding 
      ! -- area near the opposite edge of the domain
      ! - In that case, we have to take into account simulation domain boundaries
      ! in  MPI exchanges (periodic, reflective, PML etc.) 
      IF(is_open_bc) THEN
        size_to_exchange_send = 0_idp ! if absorbing_bcs then don't check
                                      ! intersection outside the simulation
                                      ! domain
      ELSE IF(.NOT. is_open_bc) THEN
        index_ff = MODULO(iz2min,n_global)
        index_fl = MODULO(iz2max,n_global)
        size_to_exchange_send = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + &
        1_idp, 0_idp)
        first_cell_send = MAX(index_rf,index_ff) - index_rf
      ENDIF
   !  -- The next two cases correspond to a global field partially lapping
   !  -- domain boundary ghost region
   !  -- In this case only a portion of global field indexes are flipped to the
   !  -- corresponding area near the opposite edge of the domain
   ! - In that case, we have to take into account simulation domain boundaries
   ! in  MPI exchanges (periodic, reflective, PML etc.) 

   ELSE IF (select_case ==2) THEN
     ! -- First half intersection  
      IF(is_open_bc) THEN
        size_to_exchange_send = 0_idp ! if absorbing bcs then don't check
                                      ! intersection outside the simulation domain 
      ELSE IF(.NOT. is_open_bc) THEN
        index_ff = MODULO(iz2min,n_global)
        index_fl = n_global - 1_idp
        size_to_exchange_send = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + &
        1_idp, 0_idp)
        first_cell_send = MAX(index_rf,index_ff) - index_rf
      ENDIF
      ! -- If first half intersection returns 0 check if the rest of the global
      ! -- field intersects 
  
      ! -- Second half intersection
      IF (size_to_exchange_send .EQ. 0_idp) THEN
        index_ff = 0_idp
        index_fl = iz2max
        size_to_exchange_send = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + &
        1_idp, 0_idp)
        first_cell_send = MAX(index_rf,index_ff) - index_rf
      ENDIF

    ELSE IF (select_case == 3) THEN
      IF(is_open_bc) THEN
        size_to_exchange_send = 0_idp 
      ELSE IF(.NOT. is_open_bc) THEN  
        index_ff = 0_idp
        index_fl = MODULO(iz2max,n_global)
        size_to_exchange_send = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + &
        1_idp, 0_idp)
        first_cell_send = MAX(index_rf,index_ff) - index_rf
      ENDIF
      ! -- If first half intersection returns 0 check if the rest of the global
      ! -- field intersects 
      ! -- Second half intersection
      IF (size_to_exchange_send .EQ. 0_idp) THEN
        index_ff = iz2min
        index_fl = n_global - 1_idp
        size_to_exchange_send = MAX(MIN(index_rl,index_fl) -MAX(index_rf,index_ff) + &
        1_idp, 0_idp)
        first_cell_send = MAX(index_rf,index_ff) - index_rf
      ENDIF
    ENDIF
    ! -- first cell send need to be set to a value contained in lower_b and
    ! -- upper_b of local field even if no exchanges are performed
    IF (size_to_exchange_send==0_idp) first_cell_send=0_idp  
#endif
  END SUBROUTINE compute_send_recv_sizes_and_index_l2g_copies

  ! ______________________________________________________________________________________
  !> @brief
  !> This routine computes for current rank, the starting indices and sizes (in the FFT 
  !> field arrays) of data to be sent to / received from local arrays of another rank. 
  !> It thus basically treats send/rev exchanges required for copies from/to the FFT
  !> arrays. This is useful for the copy operation from the local PIC 
  !> grif to the global FFT grid before the global FFT operation amongst a given MPI group. 
  !> This routine computes along a given axis (X,Y or Z): 
  !> (i) intersection of the FFT field array (with group guard cells) of current rank
  !> with the local field array (without guard cells) from another rank. 
  !> Returns size and first cell of data to be received from the local array
  !> to the FFT  arrays
  !> (ii) intersection of the FFT field array  (without group guard cells) of current rank
  !> with the local field array (without guard cells) from another rank.  
  !> Returns size and first cell of data to be sent to local array from the FFT array
  !> @param[in] iz1min INTEGER(idp) - min index of current rank along axis (X,Y or Z)
  !> @param[in] iz1max INTEGER(idp) - max index of current rank along axis (X,Y or Z)
  !> @param[in] iz2min INTEGER(idp) - min index of other rank along axis (X,Y or Z)
  !> @param[in] iz2max INTEGER8(idp) - max index of other rank along axis (X,Y or Z)
  !> @param[in] cell_g_min INTEGER(idp) - index of first cell of group (without
  !> guards cells)
  !> @param[in] cell_g_max INTEGER(idp) - index of last cell of group (without
  !> guardcells)
  !> @param[in,out] size_to_exchange_recv INTEGER(idp) - size of data to be received 
  !> from the FFT array of other rank 
  !> @param[in,out] first_cell_recv INTEGER(idp) - position in local array where data 
  !> is to be received (from the FFT array of other rank)
  !> @param[in,out] size_to_exchange_send INTEGER(idp) - size of data to be sent to 
  !> other rank 
  !> @param[in,out] first_cell_send INTEGER(idp) - starting position in local array of 
  !> data to be sent to other rank 
  !> @param[in] n_global
  !> @param[in] n_guards
  !> @author
  !> Haithem Kallala 
  !> @date
  !> Creation 2017
  ! ______________________________________________________________________________________
SUBROUTINE compute_send_recv_sizes_and_index_g2l_copies                           &
        (iz1min,iz1max,iz2min,iz2max,                                             &
        size_to_exchange_recv,first_cell_recv,                                    &
        size_to_exchange_send,first_cell_send,cell_g_min,cell_g_max,n_global, &
        n_guards,is_open_bc)
    USE picsar_precision, ONLY: idp, lp
    INTEGER(idp) , INTENT(IN)  :: iz1min, iz1max, iz2min, iz2max
    INTEGER(idp)               :: iz1min_temp
    INTEGER(idp) , INTENT(INOUT)  ::   size_to_exchange_recv,first_cell_recv,     &
    size_to_exchange_send,first_cell_send
    INTEGER(idp)                  :: index_rf, index_ff, index_rl, index_fl
    INTEGER(idp)                  :: select_case
    INTEGER(idp) ,INTENT(IN)      :: cell_g_min, cell_g_max
    INTEGER(idp), INTENT(IN)      :: n_global,n_guards
    LOGICAL(lp), INTENT(IN)       :: is_open_bc   

#if defined(FFTW)
    !-- set a flag for different possible configurations for global field limits
    !-- If global field laps a physical domain boundary, max cell index and min
    !-- cell index need to be recalculated accordingly.

    ! - Most trivial case - current proc needs copies from inner simulation domain 
    IF(iz1min .GE. 0_idp .AND. iz1max .GE. 0_idp .AND. iz1min .LT. n_global      &
    .AND. iz1max .LT. n_global)  THEN
       select_case = 0_idp  
    ! - Current proc needs copies from outer guard regions of the simulation domain only 
    ! - In that case, we have to take into account simulation domain boundaries in 
    ! - MPI exchanges (periodic, reflective, PML etc.) 
    ELSE IF((iz1min .LT. 0_idp .AND. iz1max .LT. 0_idp) .OR.                     & 
    (iz1min .GE. n_global .AND. iz1max .GE. n_global)) THEN
       select_case = 1     
    ! - Current proc needs copies from outer guard regions of the simulation domain 
    ! - and from the inner simulation domain grid. 
    ! - In that case, we have to take into account simulation domain boundaries in 
    ! - MPI exchanges (periodic, reflective, PML etc.) 
    ELSE IF (iz1min .LT. 0_idp .AND. iz1max .GE. 0_idp) THEN 
       select_case = 2      
    ! - Current proc needs copies from outer guard regions of the simulation domain 
    ! - and from the inner simulation domain grid.  
    ! - In that case, we have to take into account simulation domain boundaries in 
    ! - MPI exchanges (periodic, reflective, PML etc.) 
    ELSE IF (iz1min .LT. n_global .AND. iz1max .GE. n_global) THEN
       select_case = 3 
   ENDIF
    index_rf = iz2min ! low index of local grid array (excluding its own guard cells) 
    index_rl = iz2max ! upper index of local grid array (excluding its own guard cells) 
    IF(select_case == 0_idp) THEN   ! most trivial case 
      index_ff = iz1min ! low index of FFT array (including its own guard cells)
      index_fl = iz1max ! upper index of FFT array (including its own guard cells)
      size_to_exchange_recv = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) +   &
       1_idp,0_idp)
      first_cell_recv = MAX(index_ff,index_rf) - index_ff + 1
    ELSE IF(select_case == 1_idp) THEN
      IF(is_open_bc) THEN
        size_to_exchange_recv = 0_idp ! if absorbing bcs then don't check
                                      ! intersection outside the simulation
                                      ! domain
      ELSE IF(.NOT. is_open_bc) THEN
        index_ff = MODULO(iz1min,n_global) ! This treats the case of periodic bc 
        index_fl = MODULO(iz1max,n_global) ! This treats the case of periodic bc
        size_to_exchange_recv = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) +   &
        1_idp,0_idp) ! - Yields 0 if no intersection 
        first_cell_recv = MAX(index_ff,index_rf) - index_ff + 1 
      ENDIF
    ELSE IF(select_case == 2) THEN
      IF(is_open_bc) THEN
        size_to_exchange_recv = 0_idp  ! if absorbing bcs then don't check
                                       ! intersection outside the simulation
                                       ! domain
      ELSE IF(.NOT. is_open_bc) THEN 
      ! - First compute intersection of part outside the simulation domain 
        index_ff = MODULO(iz1min,n_global) ! - This assume periodic bc 
        index_fl = n_global - 1_idp        ! - This assume periodic bc 
        size_to_exchange_recv = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) +   &
        1_idp,0_idp) !- Yields 0 if no intersection 
        first_cell_recv = MAX(index_ff,index_rf) - index_ff + 1 
      ENDIF
      ! - If intersection of this part is none --> Check intersection with the part 
      ! - located exclusively inside the simulation domain 
      IF(size_to_exchange_recv .EQ. 0_idp) THEN
        index_ff = 0_idp 
        index_fl = iz1max
        size_to_exchange_recv = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + &
         1_idp,0_idp) ! - Yields 0 if no intersection 
        first_cell_recv = MAX(index_ff,index_rf) - index_ff + 1 - iz1min 
      ENDIF
    ELSE IF(select_case == 3) THEN
      IF(is_open_bc) THEN
        size_to_exchange_recv = 0_idp  ! if absorbing bcs then don't check
                                       ! intersection outside the simulation domain 
      ELSE IF(.NOT. is_open_bc) THEN
      ! - First compute intersection outside the simulation domain 
        index_ff = 0_idp
        index_fl = MODULO(iz1max,n_global)
        size_to_exchange_recv = MAX(MIN(index_rl,index_fl) -MAX(index_rf,index_ff) +&
         1_idp,0_idp) ! - Yields 0 if no intersection 
        first_cell_recv = MAX(index_ff,index_rf) - index_ff + 1 + (n_global -iz1min)
      ENDIF
      ! - If intersection of this part is none --> Check intersection with the part 
      ! - located exclusively inside the simulation domain 
      IF(size_to_exchange_recv .EQ. 0_idp) THEN
      index_ff = iz1min
      index_fl = n_global - 1_idp
      size_to_exchange_recv = MAX(MIN(index_rl,index_fl) -MAX(index_rf,index_ff) +&
       1_idp,0_idp) ! - Yields 0 if no intersection 
      first_cell_recv = MAX(index_ff,index_rf) - index_ff + 1

      ENDIF
    ENDIF 
    ! -- first cell recv need to be set to a value contained in lower_b and
    ! -- upper_b of hybrid field even if no exchanges are performed
    IF(size_to_exchange_recv==0_idp) first_cell_recv = 1_idp
 
    size_to_exchange_send = 0_idp
    first_cell_send = 0_idp 
    index_rf = iz2min
    index_rl = iz2max 

    index_ff = MAX(iz1min,cell_g_min)
    index_fl = MIN(iz1max,cell_g_max)    
    size_to_exchange_send = MAX(MIN(index_rl,index_fl) - MAX(index_rf,index_ff) + 1_idp,0) 

    first_cell_send = MAX(index_ff,index_rf) - iz1min + 1
    IF(select_case == 1_idp) size_to_exchange_send = 0_idp
    IF(index_ff .GT. index_fl) size_to_exchange_send = 0_idp

     ! -- first cell recv need to be set to a value contained in lower_b and
     ! -- upper_b of hybrid field even if no exchanges are performed
    IF(size_to_exchange_send==0_idp) first_cell_send = 1_idp
#endif
  END SUBROUTINE compute_send_recv_sizes_and_index_g2l_copies
 
 ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes the total time per part for particle subroutines
  !> (i.e  mainly particle push, field gathering and current deposition)
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE compute_time_per_part()
    IMPLICIT NONE
    REAL(num) :: global_time_part
    CALL get_local_number_of_part(npart_local)
    global_time_part=0.
    ! Get max time per it
    IF (npart_local .EQ. 0) THEN
      local_time_part=0
    ENDIF
    CALL MPI_ALLREDUCE(local_time_part, global_time_part, 1_isp, MPI_REAL8, MPI_SUM,  &
    comm, errcode)
    CALL MPI_ALLREDUCE(npart_local, npart_global, 1_isp, MPI_REAL8, MPI_SUM, comm,    &
    errcode)
    IF (npart_global .EQ. 0_idp) THEN
      global_time_per_part=0_num
    ELSE
      global_time_per_part=global_time_part/npart_global
    ENDIF
  END SUBROUTINE compute_time_per_part


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes the total time per cell for em field subroutines
  !> (i.e field pusher)
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE compute_time_per_cell()
    IMPLICIT NONE
    REAL(num) :: global_time_cell
    global_time_cell=0.
    ! Get max time per it
    CALL MPI_ALLREDUCE(local_time_cell, global_time_cell, 1_isp, MPI_REAL8, MPI_SUM,  &
    comm, errcode)
    SELECT CASE(c_dim)
    CASE(2)
      global_time_per_cell=global_time_cell/(nx_global*nz_global)
    CASE DEFAULT! #3D Case
      global_time_per_cell=global_time_cell/(nx_global*ny_global*nz_global)
    END SELECT
  END SUBROUTINE compute_time_per_cell


  ! ______________________________________________________________________________________
  !> @brief
  !> Get the maximum time of all process per iteration.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_max_time_per_it()
    IMPLICIT NONE
    ! Get max time per it
    CALL MPI_ALLREDUCE(mpitime_per_it, max_time_per_it, 1_isp, MPI_REAL8, MPI_MAX,    &
    comm, errcode)

  END SUBROUTINE get_max_time_per_it

  ! ______________________________________________________________________________________
  !> @brief
  !> Get the minimum time of all process per iteration.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_min_time_per_it()
    IMPLICIT NONE
    ! Get max time per it
    CALL MPI_ALLREDUCE(mpitime_per_it, min_time_per_it, 1_isp, MPI_REAL8, MPI_MIN,    &
    comm, errcode)
  END SUBROUTINE get_min_time_per_it

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine needs a better description.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE compute_new_split_2D(tppart, tpcell, nx_glob, nz_glob, ncxmin, ncxmax,   &
    nczmin, nczmax, npx, npz)
    IMPLICIT NONE
    REAL(num), INTENT(IN) :: tppart, tpcell
    INTEGER(idp), INTENT(IN) :: nx_glob, nz_glob, npx, npz
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npx) :: ncxmin, ncxmax
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npz) :: nczmin, nczmax
    REAL(num), DIMENSION(:), ALLOCATABLE :: load_on_x, load_on_z
    ALLOCATE(load_on_x(0:nx_glob-1), load_on_z(0:nz_glob-1))
    load_on_x=0.
    load_on_z=0.
    ! Compute load in X and compute new split in X
    CALL get_projected_load_on_x(nx_glob, load_on_x, tppart, tpcell)
    CALL balance_in_dir(load_on_x, nx_glob, npx, ncxmin, ncxmax)
    ! Compute load in X and compute new split in Z
    CALL get_projected_load_on_z(nz_glob, load_on_z, tppart, tpcell)
    CALL balance_in_dir(load_on_z, nz_glob, npz, nczmin, nczmax)

    DEALLOCATE(load_on_x, load_on_z)

  END SUBROUTINE compute_new_split_2D

  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine needs a better description.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE compute_new_split(tppart, tpcell, nx_glob, ny_glob, nz_glob, ncxmin,     &
    ncxmax, ncymin, ncymax, nczmin, nczmax, npx, npy, npz)
    IMPLICIT NONE
    REAL(num), INTENT(IN) :: tppart, tpcell
    INTEGER(idp), INTENT(IN) :: nx_glob, ny_glob, nz_glob, npx, npy, npz
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npx) :: ncxmin, ncxmax
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npy) :: ncymin, ncymax
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npz) :: nczmin, nczmax
    REAL(num), DIMENSION(:), ALLOCATABLE :: load_on_x, load_on_y, load_on_z
    ALLOCATE(load_on_x(0:nx_glob-1), load_on_y(0:ny_glob-1), load_on_z(0:nz_glob-1))
    load_on_x=0.
    load_on_y=0.
    load_on_z=0.

    ! Compute load in X and compute new split in X
    CALL get_projected_load_on_x(nx_glob, load_on_x, tppart, tpcell)
    CALL balance_in_dir(load_on_x, nx_glob, npx, ncxmin, ncxmax)
    ! Compute load in X and compute new split in Y
    CALL get_projected_load_on_y(ny_glob, load_on_y, tppart, tpcell)
    CALL balance_in_dir(load_on_y, ny_glob, npy, ncymin, ncymax)
    ! Compute load in X and compute new split in Z
    CALL get_projected_load_on_z(nz_glob, load_on_z, tppart, tpcell)
    CALL balance_in_dir(load_on_z, nz_glob, npz, nczmin, nczmax)

    DEALLOCATE(load_on_x, load_on_y, load_on_z)

  END SUBROUTINE compute_new_split



  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes new load in direction dir based on the projected
  !> load_in_dir
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE balance_in_dir(load_in_dir, ncellmaxdir, nproc_in_dir, idirmin, idirmax)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nproc_in_dir, ncellmaxdir
    REAL(num), DIMENSION(0:ncellmaxdir-1), INTENT(IN) :: load_in_dir
    INTEGER(idp), DIMENSION(0:nproc_in_dir-1), INTENT(IN OUT) :: idirmin, idirmax
    INTEGER(idp) :: iproc, icell
    REAL(num) :: balanced_load=0_num, curr_balanced_load=0_num, curr_proc_load=0_num
    LOGICAL(lp)  :: not_balanced
    REAL(num) :: delta
    REAL(num), DIMENSION(0:nproc_in_dir-1) :: load_per_proc

    load_per_proc=0_num
    balanced_load=SUM(load_in_dir)/(nproc_in_dir)
    icell=0
    curr_balanced_load=balanced_load
    DO iproc=0, nproc_in_dir-1
      not_balanced=.TRUE.
      curr_proc_load=0_num
      idirmin(iproc)=icell
      delta=0_num
      DO WHILE((not_balanced) .AND. (icell .LT. ncellmaxdir-1))
        curr_proc_load=curr_proc_load+load_in_dir(icell)
        IF(curr_proc_load .GE. curr_balanced_load)  THEN
          not_balanced=.FALSE.
        ELSE
          icell=icell+1
        ENDIF
      END DO

      idirmax(iproc)=icell
      load_per_proc(iproc)=curr_proc_load
      delta=(iproc+1_num)*balanced_load-SUM(load_per_proc)
      curr_balanced_load=balanced_load+delta
      icell=icell+1
      ! Sanity check
      IF ((iproc .EQ. nproc_in_dir-1) .AND. (icell .LT. ncellmaxdir-1)) THEN
        idirmax(iproc)=ncellmaxdir-1
      ENDIF
    END DO


  END SUBROUTINE balance_in_dir


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes the computational load projected on X-Axis
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_projected_load_on_x(nxg, load_on_x, time_per_part, time_per_cell)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxg
    REAL(num), INTENT(IN OUT), DIMENSION(0:nxg-1) :: load_on_x
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icellx
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile

    ALLOCATE(load_part(0:nxg-1))
    ALLOCATE(load_part_sum(0:nxg-1))
    load_part=0
    load_part_sum=0

    ! Get local distribution of particles along X-axis
    DO ispecies=1, nspecies
      curr => species_parray(ispecies)
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile => curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            DO ip=1, count
              icellx=FLOOR((curr_tile%part_x(ip)-x_grid_min)/dx)
              load_part(icellx)=load_part(icellx)+1
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Get contributions on X-axis from other MPI domains
    CALL MPI_ALLREDUCE(load_part, load_part_sum, INT(nxg, isp), MPI_INTEGER8,         &
    MPI_SUM, comm, errcode)

    ! Computes load_in_x
    SELECT CASE (c_dim)
    CASE (2)! 2D CASE
      load_on_x = load_part_sum*time_per_part + nz_global*time_per_cell
    CASE DEFAULT! 3D CASE
      load_on_x = load_part_sum*time_per_part + ny_global * nz_global*time_per_cell
    END SELECT

    DEALLOCATE(load_part, load_part_sum)

  END SUBROUTINE get_projected_load_on_x


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes the computational load projected on Y-Axis
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_projected_load_on_y(nyg, load_on_y, time_per_part, time_per_cell)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nyg
    REAL(num), INTENT(IN OUT), DIMENSION(0:nyg-1) :: load_on_y
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icelly
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile

    ALLOCATE(load_part(0:nyg-1))
    ALLOCATE(load_part_sum(0:nyg-1))
    load_part=0
    load_part_sum=0
    ! Get local distribution of particles along Y-axis
    DO ispecies=1, nspecies
      curr => species_parray(ispecies)
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile => curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            DO ip=1, count
              icelly=FLOOR((curr_tile%part_y(ip)-y_grid_min)/dy)
              load_part(icelly)=load_part(icelly)+1
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Get contributions on Y-axis from other MPI domains
    CALL MPI_ALLREDUCE(load_part, load_part_sum, INT(nyg, isp), MPI_INTEGER8,         &
    MPI_SUM, comm, errcode)

    ! Computes load_in_y
    load_on_y = load_part_sum*time_per_part + nx_global * nz_global*time_per_cell

    DEALLOCATE(load_part, load_part_sum)

  END SUBROUTINE get_projected_load_on_y


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine computes the computational load projected on Z-Axis
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE get_projected_load_on_z(nzg, load_on_z, time_per_part, time_per_cell)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nzg
    REAL(num), INTENT(IN OUT), DIMENSION(0:nzg-1) :: load_on_z
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icellz
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile), POINTER :: curr_tile

    ALLOCATE(load_part(0:nzg-1))
    ALLOCATE(load_part_sum(0:nzg-1))
    load_part=0
    load_part_sum=0
    ! Get local distribution of particles along Z-axis
    DO ispecies=1, nspecies
      curr => species_parray(ispecies)
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile => curr%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            DO ip=1, count
              icellz=FLOOR((curr_tile%part_z(ip)-z_grid_min)/dz)
              load_part(icellz)=load_part(icellz)+1
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Get contributions on Z-axis from other MPI domains
    CALL MPI_ALLREDUCE(load_part, load_part_sum, INT(nzg, isp), MPI_INTEGER8,         &
    MPI_SUM, comm, errcode)

    ! Computes load_in_z
    SELECT CASE (c_dim)
    CASE (2)! 2D CASE
      load_on_z = load_part_sum*time_per_part + nx_global *time_per_cell
    CASE DEFAULT! 3D CASE
      load_on_z = load_part_sum*time_per_part + nx_global * ny_global*time_per_cell
    END SELECT
    DEALLOCATE(load_part, load_part_sum)

  END SUBROUTINE get_projected_load_on_z


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine create new array_of_tiles for each species
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2016
  ! ______________________________________________________________________________________
  SUBROUTINE create_new_tile_split()
#ifdef _OPENMP
    USE omp_lib
#endif
    IMPLICIT NONE
    TYPE(particle_species), DIMENSION(:), ALLOCATABLE, TARGET :: new_species_parray
    TYPE(particle_species), POINTER :: currsp, currsp_new
    TYPE(particle_tile), POINTER :: curr_tile
    INTEGER(idp)  :: nxmin, nxmax, nymin, nymax, nzmin, nzmax, oxmin, oxmax, oymin,   &
    oymax, ozmin, ozmax
    INTEGER(idp) :: ix, iy, iz, ip, indx, indy, indz, ispecies, count
    INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
    REAL(num) :: partx, party, partz, partux, partuy, partuz, gaminv
    REAL(num), DIMENSION(npid) :: partpid
    INTEGER(idp) :: ntilex_new, ntiley_new, ntilez_new, nthreads_tot

#ifdef _OPENMP
    nthreads_tot=OMP_GET_MAX_THREADS()
#else
    nthreads_tot=1
#endif

    IF (nthreads_tot .GT. 1) THEN
      ! Udpate optimal number of tiles
      SELECT CASE(c_dim)
      CASE (2)
        ntilex_new = MAX(1_idp, nx/35)
        ntiley_new = 1
        ntilez_new = MAX(1_idp, nz/35)
      CASE DEFAULT
        ntilex_new = MAX(1_idp, nx/10)
        ntiley_new = MAX(1_idp, ny/10)
        ntilez_new = MAX(1_idp, nz/10)
      END SELECT
    ELSE
      ntilex_new = 1
      ntiley_new = 1
      ntilez_new = 1
    ENDIF

    ! Allocate new species array
    ALLOCATE(new_species_parray(nspecies))
    ! Copy properties of each species from old array to new array
    DO ispecies=1, nspecies
      currsp=>species_parray(ispecies)
      currsp_new=>new_species_parray(ispecies)
      currsp_new%charge=currsp%charge
      currsp_new%mass=currsp%mass
      currsp_new%x_min=currsp%x_min
      currsp_new%x_max=currsp%x_max
      currsp_new%y_min=currsp%y_min
      currsp_new%y_max=currsp%y_max
      currsp_new%z_min=currsp%z_min
      currsp_new%z_max=currsp%z_max
      currsp_new%vth_x=currsp%vth_x
      currsp_new%vth_y=currsp%vth_y
      currsp_new%vth_z=currsp%vth_z
      currsp_new%vdrift_x=currsp%vdrift_x
      currsp_new%vdrift_y=currsp%vdrift_y
      currsp_new%vdrift_z=currsp%vdrift_z
      currsp_new%nppcell=currsp%nppcell
    END DO

    CALL set_tile_split_for_species(new_species_parray, nspecies, ntilex_new,         &
    ntiley_new, ntilez_new, nx_grid, ny_grid, nz_grid, x_min_local, y_min_local,      &
    z_min_local, x_max_local, y_max_local, z_max_local)


    ! Deallocate former grid tile array/ ALLOCATE new one
    DEALLOCATE(aofgrid_tiles)
    ALLOCATE(aofgrid_tiles(ntilex_new, ntiley_new, ntilez_new))

    ! Init new tile arrays of new species array
    CALL init_tile_arrays_for_species(nspecies, new_species_parray, aofgrid_tiles,    &
    ntilex_new, ntiley_new, ntilez_new)

    ! Copy particles from former species array to new species array
    DO ispecies=1, nspecies
      currsp=>species_parray(ispecies)
      currsp_new=>new_species_parray(ispecies)
      ! Get first tiles dimensions (may be different from last tile)
      nx0_grid_tile = currsp_new%array_of_tiles(1, 1, 1)%nx_grid_tile
      ny0_grid_tile = currsp_new%array_of_tiles(1, 1, 1)%ny_grid_tile
      nz0_grid_tile = currsp_new%array_of_tiles(1, 1, 1)%nz_grid_tile
      DO iz=1, ntilez
        DO iy=1, ntiley
          DO ix=1, ntilex
            curr_tile=>currsp%array_of_tiles(ix, iy, iz)
            count=curr_tile%np_tile(1)
            SELECT CASE (c_dim)
            CASE (2)
              DO ip=count, 1, -1
                partx=curr_tile%part_x(ip)
                party=curr_tile%part_y(ip)
                partz=curr_tile%part_z(ip)
                partux=curr_tile%part_ux(ip)
                partuy=curr_tile%part_uy(ip)
                partuz=curr_tile%part_uz(ip)
                gaminv=curr_tile%part_gaminv(ip)
                partpid=curr_tile%pid(ip, 1:npid)
                ! CASE 1: particle outside MPI domain temporarily put it
                ! in the first tile of new species array
                IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local) .OR.       &
                (partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) THEN
                CALL add_particle_at_tile(currsp_new, 1_idp, 1_idp, 1_idp, partx,     &
                party, partz, partux, partuy, partuz, gaminv, partpid)
                ! CASE 2: particle is in the new domain just add
                ! it to proper tile of new species array
              ELSE
                indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),     &
                idp)+1, ntilex_new)
                indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz),     &
                idp)+1, ntilez_new)
                CALL add_particle_at_tile(currsp_new, indx, 1_idp, indz, partx,       &
                party, partz, partux, partuy, partuz, gaminv, partpid)

              ENDIF
              currsp_new%species_npart=currsp_new%species_npart+1
              CALL rm_particles_from_species(currsp, ix, iy, iz, ip)
            END DO
          CASE DEFAULT
            DO ip=count, 1, -1
              partx=curr_tile%part_x(ip)
              party=curr_tile%part_y(ip)
              partz=curr_tile%part_z(ip)
              partux=curr_tile%part_ux(ip)
              partuy=curr_tile%part_uy(ip)
              partuz=curr_tile%part_uz(ip)
              gaminv=curr_tile%part_gaminv(ip)
              partpid=curr_tile%pid(ip, 1:npid)
              ! CASE 1: particle outside MPI domain temporarily put it
              ! in the first tile of new species array
              IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local) .OR. (party  &
              .LT. y_min_local) .OR. (party .GE. y_max_local) .OR. (partz .LT.        &
              z_min_local) .OR. (partz .GE. z_max_local)) THEN
              CALL add_particle_at_tile(currsp_new, 1_idp, 1_idp, 1_idp, partx,       &
              party, partz, partux, partuy, partuz, gaminv, partpid)
              ! CASE 2: particle is in the new domain just add it
              !  to proper tile of new species array
            ELSE
              indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx),       &
              idp)+1, ntilex_new)
              indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy),       &
              idp)+1, ntiley_new)
              indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz),       &
              idp)+1, ntilez_new)
              CALL add_particle_at_tile(currsp_new, indx, indy, indz, partx, party,   &
              partz, partux, partuy, partuz, gaminv, partpid)
            ENDIF
            currsp_new%species_npart=currsp_new%species_npart+1
            CALL rm_particles_from_species(currsp, ix, iy, iz, ip)
          END DO
        END SELECT
      END DO
    END DO
  END DO
END DO


! Update tile sizes
ntilex=ntilex_new
ntiley=ntiley_new
ntilez=ntilez_new


! Deallocate old array
DEALLOCATE(species_parray)

! Reallocate with new dimensions
ALLOCATE(species_parray(nspecies))

! Copy species properties
DO ispecies=1, nspecies
  currsp=>species_parray(ispecies)
  currsp_new=>new_species_parray(ispecies)
  currsp%charge=currsp_new%charge
  currsp%mass=currsp_new%mass
  currsp%x_min=currsp_new%x_min
  currsp%x_max=currsp_new%x_max
  currsp%y_min=currsp_new%y_min
  currsp%y_max=currsp_new%y_max
  currsp%z_min=currsp_new%z_min
  currsp%z_max=currsp_new%z_max
  currsp%vth_x=currsp_new%vth_x
  currsp%vth_y=currsp_new%vth_y
  currsp%vth_z=currsp_new%vth_z
  currsp%vdrift_x=currsp_new%vdrift_x
  currsp%vdrift_y=currsp_new%vdrift_y
  currsp%vdrift_z=currsp_new%vdrift_z
  currsp%nppcell=currsp_new%nppcell
END DO


! Set tile split for species_parray
CALL set_tile_split_for_species(species_parray, nspecies, ntilex, ntiley, ntilez,     &
nx_grid, ny_grid, nz_grid, x_min_local, y_min_local, z_min_local, x_max_local,        &
y_max_local, z_max_local)

DEALLOCATE(aofgrid_tiles)
ALLOCATE(aofgrid_tiles(ntilex, ntiley, ntilez))

CALL init_tile_arrays_for_species(nspecies, species_parray, aofgrid_tiles, ntilex,    &
ntiley, ntilez)

! Copy particles
! Copy particles from former tiles in first tile of new species array
DO ispecies=1, nspecies
  currsp=>species_parray(ispecies)
  currsp_new=>new_species_parray(ispecies)
  DO iz=1, ntilez
    DO iy=1, ntiley
      DO ix=1, ntilex
        curr_tile=>currsp_new%array_of_tiles(ix, iy, iz)
        count=curr_tile%np_tile(1)
        DO ip=count, 1, -1
          partx=curr_tile%part_x(ip)
          party=curr_tile%part_y(ip)
          partz=curr_tile%part_z(ip)
          partux=curr_tile%part_ux(ip)
          partuy=curr_tile%part_uy(ip)
          partuz=curr_tile%part_uz(ip)
          gaminv=curr_tile%part_gaminv(ip)
          partpid=curr_tile%pid(ip, 1:npid)
          CALL add_particle_at_tile(currsp, ix, iy, iz, partx, party, partz, partux,  &
          partuy, partuz, gaminv, partpid)
          CALL rm_particles_from_species(currsp_new, ix, iy, iz, ip)
          currsp%species_npart=currsp%species_npart+1
        END DO
      END DO
    END DO
  END DO
END DO

!Deallocate new_species_parray
DEALLOCATE(new_species_parray)

END SUBROUTINE create_new_tile_split

! ________________________________________________________________________________________
!> @brief
!> This subroutines needs a description.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ______________________________________________________________________________________
SUBROUTINE remap_particles(ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, ix1new,    &
ix2new, iy1new, iy2new, iz1new, iz2new, ncxmin, ncxmax, ncymin, ncymax, nczmin,       &
nczmax, iproc, ncpus, npx, npy, npz, l_cart_comm)
USE iso_c_binding

IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: iproc, ncpus, npx, npy, npz
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ix1old, ix2old, iy1old, iy2old,     &
iz1old, iz2old
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ix1new, ix2new, iy1new, iy2new,     &
iz1new, iz2new
INTEGER(idp), INTENT(IN), DIMENSION(0:npx-1) :: ncxmin, ncxmax
INTEGER(idp), INTENT(IN), DIMENSION(0:npy-1) :: ncymin, ncymax
INTEGER(idp), INTENT(IN), DIMENSION(0:npz-1) :: nczmin, nczmax
LOGICAL(lp), INTENT(IN) :: l_cart_comm
INTEGER(idp) :: i, ipart, ixtile, iytile, iztile, nmax, ispec, ispecies, npcurr
INTEGER(idp) :: ix3min, ix3max, iy3min, iy3max, iz3min, iz3max
INTEGER(idp) :: ix1newip, ix2newip, iy1newip, iy2newip, iz1newip, iz2newip
INTEGER(idp) :: ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip, iz2oldip
INTEGER(isp) :: mpitag, count
INTEGER(isp), ALLOCATABLE, DIMENSION(:) :: recv_rank, send_rank, requests
INTEGER(idp), ALLOCATABLE, DIMENSION(:, :) :: npart_recv, npart_send
INTEGER(idp), PARAMETER :: nmax_neighbours=10**3
INTEGER(idp) :: nvar
INTEGER(idp) :: nsend, nrecv, nsdat, nrdat, ibuff, curr_rank, iprocx, iprocy, iprocz, &
icx, icy, icz, isend
LOGICAL(lp)  :: l_is_intersection
INTEGER(idp), ALLOCATABLE, DIMENSION (:) :: nptoexch
REAL(num), ALLOCATABLE, DIMENSION (:, :) :: sendbuff
REAL(num), ALLOCATABLE, DIMENSION (:, :) :: recvbuff
REAL(num) :: part_xyz
TYPE(particle_species), POINTER :: currsp
TYPE(particle_tile), POINTER :: curr

mpitag=0_isp
nvar=7+npid

ALLOCATE(recv_rank(nmax_neighbours), send_rank(nmax_neighbours))
recv_rank=-1
send_rank=-1
nsend=0
nrecv=0

! ---- INTERSECTION OF NEW DOMAIN PROC WITH OLD DOMAINS =/ PROC
ix1newip = ix1new(iproc)
ix2newip = ix2new(iproc)
iy1newip = iy1new(iproc)
iy2newip = iy2new(iproc)
iz1newip = iz1new(iproc)
iz2newip = iz2new(iproc)
DO i=0, ncpus-1
  CALL get_3Dintersection(ix1newip, ix2newip, iy1newip, iy2newip, iz1newip, iz2newip, &
  ix1old(i), ix2old(i), iy1old(i), iy2old(i), iz1old(i), iz2old(i), ix3min, ix3max,   &
  iy3min, iy3max, iz3min, iz3max, l_is_intersection)
  ! Case i == iproc cycle
  IF (i .EQ. iproc) CYCLE
  ! Found intersection area between new proc and old adjacent proc
  ! Creates RECV TYPE FOR THIS Volume
  IF (l_is_intersection .AND. (i .NE. iproc)) THEN
    !--- Put rank in "receive from" list
    nrecv=nrecv+1
    recv_rank(nrecv)=i
  ENDIF
END DO

! ---- INTERSECTION OF OLD DOMAIN PROC WITH NEW DOMAINS =/ PROC
ix1oldip = ix1old(iproc)
ix2oldip = ix2old(iproc)
iy1oldip = iy1old(iproc)
iy2oldip = iy2old(iproc)
iz1oldip = iz1old(iproc)
iz2oldip = iz2old(iproc)
DO i=0, ncpus-1
  CALL get_3Dintersection(ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip, iz2oldip, &
  ix1new(i), ix2new(i), iy1new(i), iy2new(i), iz1new(i), iz2new(i), ix3min, ix3max,   &
  iy3min, iy3max, iz3min, iz3max, l_is_intersection)
  ! Case i == iproc cycle
  IF (i .EQ. iproc) CYCLE

  ! Found intersection area between old proc and new adjacent procs
  ! Create sendtype for this volume
  IF (l_is_intersection .AND. (i .NE. iproc)) THEN
    !--- Put rank in "receive from" list
    nsend=nsend+1
    send_rank(nsend)=i
  ENDIF
END DO


ALLOCATE(npart_send(nspecies, nsend), npart_recv(nspecies, nrecv),                    &
requests(nsend+nrecv))
npart_send=0_idp
npart_recv=0_idp
requests=0_isp

! ----- POST IRECV TO GET NUMBER OF PARTICLES
DO i=1, nrecv
  count=nspecies
  CALL MPI_IRECV(npart_recv(1:count, i), count, MPI_INTEGER8, recv_rank(i),           &
  MPI_ANY_TAG, comm, requests(i), errcode)
END DO

! ----- IDENTIFY PARTICLES TO BE SENT/ PLACE PARTICLES IN SEND BUFFER
ALLOCATE(sendbuff(nvar*npart_local, nsend), nptoexch(nsend))
nptoexch=0
DO ispecies=1, nspecies!LOOP ON SPECIES
  ! Init send recv buffers
  currsp => species_parray(ispecies)
  DO iztile=1, ntilez!LOOP ON TILES
    DO iytile=1, ntiley
      DO ixtile=1, ntilex
        curr=>currsp%array_of_tiles(ixtile, iytile, iztile)
        ! If not subdomain border, nothing to do
        IF (.NOT. curr%subdomain_bound) CYCLE
        ! search for outbound particles
        part_xyz=0.
        ! Identify outbounds particles and compute destination
        npcurr=curr%np_tile(1)
        DO i = npcurr, 1, -1!LOOP ON PARTICLES
          iprocx = x_coords
          iprocy = y_coords
          iprocz = z_coords
          part_xyz = curr%part_x(i)
          icx=(part_xyz-xmin)/dx
          ! Particle has left this processor
          IF ((part_xyz .LT. x_min_local) .OR. (part_xyz .GE. x_max_local)) THEN
            CALL get_proc_interval(iprocx, icx, ncxmin, ncxmax, npx)
          ENDIF
          part_xyz = curr%part_y(i)
          icy=(part_xyz-ymin)/dy
          ! Particle has left this processor
          IF (((part_xyz .LT. y_min_local) .OR. (part_xyz .GE. y_max_local))) THEN
            CALL get_proc_interval(iprocy, icy, ncymin, ncymax, npy)
          ENDIF
          part_xyz = curr%part_z(i)
          icz=(part_xyz-zmin)/dz
          ! Particle has left this processor
          IF ((part_xyz .LT. z_min_local) .OR. (part_xyz .GE. z_max_local)) THEN
            ! Find new proc using bissection algorithm (log(nproc))
            CALL get_proc_interval(iprocz, icz, nczmin, nczmax, npz)
          ENDIF
          ! Particles has to be sent to another proc
          IF((iprocx .NE. x_coords) .OR. (iprocy .NE. y_coords) .OR. (iprocz .NE.     &
          z_coords))  THEN
          ! Finds indices in buffer for curr_rank using a binary search algorithm
          CALL pxr_convertindtoproc(comm, iprocx, iprocy, iprocz, npx, npy, npz,    &
          curr_rank, l_cart_comm)
          CALL binary_search(isend, INT(curr_rank, isp), send_rank(1:nsend), nsend)
          ibuff=nvar*nptoexch(isend)+1
          sendbuff(ibuff, isend)    = curr%part_x(i)
          sendbuff(ibuff+1, isend)  = curr%part_y(i)
          sendbuff(ibuff+2, isend)  = curr%part_z(i)
          sendbuff(ibuff+3, isend)  = curr%part_ux(i)
          sendbuff(ibuff+4, isend)  = curr%part_uy(i)
          sendbuff(ibuff+5, isend)  = curr%part_uz(i)
          sendbuff(ibuff+6, isend)  = curr%part_gaminv(i)
          sendbuff(ibuff+7:ibuff+7+(npid-1), isend)  = curr%pid(i, 1:npid)
          npart_send(ispecies, isend)=npart_send(ispecies, isend)+1
          nptoexch(isend)=nptoexch(isend)+1
          ! Remove particle of current species from current tile
          CALL rm_particles_from_species(currsp, ixtile, iytile, iztile, i)
        ENDIF
      ENDDO!END LOOP ON PARTICLES
    ENDDO
  ENDDO
ENDDO! END LOOP ON TILES
END DO! End loop on species

! ----- POST ISEND FOR THE NUMBER OF PARTICLES
DO i=1, nsend
count=nspecies
CALL MPI_ISEND(npart_send(1:count, i), count, MPI_INTEGER8, send_rank(i), mpitag,   &
comm, requests(nrecv+i), errcode)
END DO


! ----- SYNC THE NUMBER OF PARTICLES BEFORE RECEIVING DATA
count=nsend+nrecv
IF (count .GT. 0) THEN
CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)
ENDIF
requests=0_isp
nsdat=0
nrdat=0

! ----- POST IRECV FOR PARTICLE DATA
nmax=nvar*MAXVAL(SUM(npart_recv, 1))
ALLOCATE(recvbuff(nmax, nrecv))
DO i=1, nrecv
count=nvar*SUM(npart_recv(:, i))
IF (count .GT. 0) THEN
  nrdat=nrdat+1
  CALL MPI_IRECV(recvbuff(1:count, i), count, MPI_DOUBLE_PRECISION, recv_rank(i),   &
  MPI_ANY_TAG, comm, requests(nrdat), errcode)
ENDIF
END DO

! ----- POST ISEND FOR PARTICLES DATA
DO i=1, nsend
count=nvar*SUM(npart_send(:, i))
IF (count .GT. 0) THEN
  nsdat=nsdat+1
  CALL MPI_ISEND(sendbuff(1:count, i), count, MPI_DOUBLE_PRECISION, send_rank(i),   &
  mpitag, comm, requests(nrdat+nsdat), errcode)
ENDIF
END DO

! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA
count=nrdat+nsdat
IF (count .GT. 0_isp) THEN
CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)
ENDIF

! ----- ADD PARTICLES FROM RECV BUFF TO SPECIES ARRAY
DO i =1, nrecv
ispec=0
DO ispecies=1, nspecies
  currsp=> species_parray(ispecies)
  DO ipart=1, nvar*npart_recv(ispecies, i), nvar
    ibuff=ispec+ipart
    CALL add_particle_to_species(currsp, recvbuff(ibuff, i), recvbuff(ibuff+1, i),  &
    recvbuff(ibuff+2, i), recvbuff(ibuff+3, i), recvbuff(ibuff+4, i),               &
    recvbuff(ibuff+5, i), recvbuff(ibuff+6, i), recvbuff(ibuff+7:ibuff+6+npid, i))
  END DO
  ispec=ispec+nvar*npart_recv(ispecies, i)
END DO
END DO


DEALLOCATE(sendbuff, recvbuff, nptoexch, npart_send, npart_recv, requests)

END SUBROUTINE remap_particles

! ________________________________________________________________________________________
!> @brief
!> This subroutines needs a description.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE remap_particles_2D(ix1old, ix2old, iz1old, iz2old, ix1new, ix2new, iz1new, &
iz2new, ncxmin, ncxmax, nczmin, nczmax, iproc, ncpus, npx, npz, l_cart_comm)
USE iso_c_binding


IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: iproc, ncpus, npx, npz
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ix1old, ix2old, iz1old, iz2old
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ix1new, ix2new, iz1new, iz2new
INTEGER(idp), INTENT(IN), DIMENSION(0:npx-1) :: ncxmin, ncxmax
INTEGER(idp), INTENT(IN), DIMENSION(0:npz-1) :: nczmin, nczmax
LOGICAL(lp), INTENT(IN) :: l_cart_comm
INTEGER(idp) :: i, ipart, ixtile, iytile, iztile, nmax, ispec, ispecies, npcurr, npy
INTEGER(idp) :: ix3min, ix3max, iz3min, iz3max
INTEGER(idp) :: ix1newip, ix2newip, iz1newip, iz2newip
INTEGER(idp) :: ix1oldip, ix2oldip, iz1oldip, iz2oldip
INTEGER(isp) :: mpitag, count
INTEGER(isp), ALLOCATABLE, DIMENSION(:) :: recv_rank, send_rank, requests
INTEGER(idp), ALLOCATABLE, DIMENSION(:, :) :: npart_recv, npart_send
INTEGER(idp), PARAMETER :: nmax_neighbours=10**3
INTEGER(idp) :: nvar
INTEGER(idp) :: nsend, nrecv, nsdat, nrdat, ibuff, curr_rank, iprocx, iprocy, iprocz, &
icx, icy, icz, isend
LOGICAL(lp)  :: l_is_intersection
INTEGER(idp), ALLOCATABLE, DIMENSION (:) :: nptoexch
REAL(num), ALLOCATABLE, DIMENSION (:, :) :: sendbuff
REAL(num), ALLOCATABLE, DIMENSION (:, :) :: recvbuff
REAL(num) :: part_xyz
TYPE(particle_species), POINTER :: currsp
TYPE(particle_tile), POINTER :: curr

npy=1! 2D CASE
mpitag=0_isp
nvar=npid+7

ALLOCATE(recv_rank(nmax_neighbours), send_rank(nmax_neighbours))
recv_rank=-1
send_rank=-1
nsend=0
nrecv=0

! ---- INTERSECTION OF NEW DOMAIN PROC WITH OLD DOMAINS =/ PROC
ix1newip = ix1new(iproc)
ix2newip = ix2new(iproc)
iz1newip = iz1new(iproc)
iz2newip = iz2new(iproc)
DO i=0, ncpus-1
CALL get_2Dintersection(ix1newip, ix2newip, iz1newip, iz2newip, ix1old(i),          &
ix2old(i), iz1old(i), iz2old(i), ix3min, ix3max, iz3min, iz3max, l_is_intersection)
! Case i == iproc cycle
IF (i .EQ. iproc) CYCLE
! Found intersection area between new proc and old adjacent proc
! Creates RECV TYPE FOR THIS Volume
IF (l_is_intersection .AND. (i .NE. iproc)) THEN
  !--- Put rank in "receive from" list
  nrecv=nrecv+1
  recv_rank(nrecv)=i
ENDIF
END DO

! ---- INTERSECTION OF OLD DOMAIN PROC WITH NEW DOMAINS =/ PROC
ix1oldip = ix1old(iproc)
ix2oldip = ix2old(iproc)
iz1oldip = iz1old(iproc)
iz2oldip = iz2old(iproc)
DO i=0, ncpus-1
CALL get_2Dintersection(ix1oldip, ix2oldip, iz1oldip, iz2oldip, ix1new(i),          &
ix2new(i), iz1new(i), iz2new(i), ix3min, ix3max, iz3min, iz3max, l_is_intersection)
! Case i == iproc cycle
IF (i .EQ. iproc) CYCLE

! Found intersection area between old proc and new adjacent procs
! Create sendtype for this volume
IF (l_is_intersection .AND. (i .NE. iproc)) THEN
  !--- Put rank in "receive from" list
  nsend=nsend+1
  send_rank(nsend)=i
ENDIF
END DO

ALLOCATE(npart_send(nspecies, nsend), npart_recv(nspecies, nrecv),                    &
requests(nsend+nrecv))
npart_send=0_idp
npart_recv=0_idp
requests=0_isp

! ----- POST IRECV TO GET NUMBER OF PARTICLES
DO i=1, nrecv
count=nspecies
CALL MPI_IRECV(npart_recv(1:count, i), count, MPI_INTEGER8, recv_rank(i),           &
MPI_ANY_TAG, comm, requests(i), errcode)
END DO

! ----- IDENTIFY PARTICLES TO BE SENT/ PLACE PARTICLES IN SEND BUFFER
ALLOCATE(sendbuff(nvar*npart_local, nsend), nptoexch(nsend))
nptoexch=0
DO ispecies=1, nspecies!LOOP ON SPECIES
! Init send recv buffers
currsp => species_parray(ispecies)
DO iztile=1, ntilez!LOOP ON TILES
  DO iytile=1, ntiley
    DO ixtile=1, ntilex
      curr=>currsp%array_of_tiles(ixtile, iytile, iztile)
      ! If not subdomain border, nothing to do
      IF (.NOT. curr%subdomain_bound) CYCLE
      ! search for outbound particles
      part_xyz=0.
      ! Identify outbounds particles and compute destination
      npcurr=curr%np_tile(1)
      DO i = npcurr, 1, -1!LOOP ON PARTICLES
        iprocx = x_coords
        iprocy = y_coords
        iprocz = z_coords
        part_xyz = curr%part_x(i)
        icx=(part_xyz-xmin)/dx
        ! Particle has left this processor
        IF ((part_xyz .LT. x_min_local) .OR. (part_xyz .GE. x_max_local)) THEN
          CALL get_proc_interval(iprocx, icx, ncxmin, ncxmax, npx)
        ENDIF
        part_xyz = curr%part_z(i)
        icz=(part_xyz-zmin)/dz
        ! Particle has left this processor
        IF ((part_xyz .LT. z_min_local) .OR. (part_xyz .GE. z_max_local)) THEN
          ! Find new proc using bissection algorithm (log(nproc))
          CALL get_proc_interval(iprocz, icz, nczmin, nczmax, npz)
        ENDIF
        ! Particles has to be sent to another proc
        IF((iprocx .NE. x_coords) .OR. (iprocy .NE. y_coords) .OR. (iprocz .NE.     &
        z_coords))  THEN
        ! Finds indices in buffer for curr_rank using a binary search algorithm
        CALL pxr_convertindtoproc(comm, iprocx, iprocy, iprocz, npx, npy, npz,    &
        curr_rank, l_cart_comm)
        CALL binary_search(isend, INT(curr_rank, isp), send_rank(1:nsend), nsend)
        ibuff=nvar*nptoexch(isend)+1
        sendbuff(ibuff, isend)    = curr%part_x(i)
        sendbuff(ibuff+1, isend)  = curr%part_y(i)
        sendbuff(ibuff+2, isend)  = curr%part_z(i)
        sendbuff(ibuff+3, isend)  = curr%part_ux(i)
        sendbuff(ibuff+4, isend)  = curr%part_uy(i)
        sendbuff(ibuff+5, isend)  = curr%part_uz(i)
        sendbuff(ibuff+6, isend)  = curr%part_gaminv(i)
        sendbuff(ibuff+7:ibuff+6+npid, isend)  = curr%pid(i, 1:npid)
        npart_send(ispecies, isend)=npart_send(ispecies, isend)+1
        nptoexch(isend)=nptoexch(isend)+1
        ! Remove particle of current species from current tile
        CALL rm_particles_from_species(currsp, ixtile, iytile, iztile, i)
      ENDIF
    ENDDO!END LOOP ON PARTICLES
  ENDDO
ENDDO
ENDDO! END LOOP ON TILES
END DO! End loop on species

! ----- POST ISEND FOR THE NUMBER OF PARTICLES
DO i=1, nsend
count=nspecies
CALL MPI_ISEND(npart_send(1:count, i), count, MPI_INTEGER8, send_rank(i), mpitag,   &
comm, requests(nrecv+i), errcode)
END DO


! ----- SYNC THE NUMBER OF PARTICLES BEFORE RECEIVING DATA
count=nsend+nrecv
IF (count .GT. 0) THEN
CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)
ENDIF
requests=0_isp
nsdat=0
nrdat=0

! ----- POST IRECV FOR PARTICLE DATA
nmax=nvar*MAXVAL(SUM(npart_recv, 1))
ALLOCATE(recvbuff(nmax, nrecv))
DO i=1, nrecv
count=nvar*SUM(npart_recv(:, i))
IF (count .GT. 0) THEN
nrdat=nrdat+1
CALL MPI_IRECV(recvbuff(1:count, i), count, MPI_DOUBLE_PRECISION, recv_rank(i),   &
MPI_ANY_TAG, comm, requests(nrdat), errcode)
ENDIF
END DO

! ----- POST ISEND FOR PARTICLES DATA
DO i=1, nsend
count=nvar*SUM(npart_send(:, i))
IF (count .GT. 0) THEN
nsdat=nsdat+1
CALL MPI_ISEND(sendbuff(1:count, i), count, MPI_DOUBLE_PRECISION, send_rank(i),   &
mpitag, comm, requests(nrdat+nsdat), errcode)
ENDIF
END DO

! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA
count=nrdat+nsdat
IF (count .GT. 0_isp) THEN
CALL MPI_WAITALL(count, requests, MPI_STATUSES_IGNORE, errcode)
ENDIF

! ----- ADD PARTICLES FROM RECV BUFF TO SPECIES ARRAY
DO i =1, nrecv
ispec=0
DO ispecies=1, nspecies
currsp=> species_parray(ispecies)
DO ipart=1, nvar*npart_recv(ispecies, i), nvar
  ibuff=ispec+ipart
  CALL add_particle_to_species(currsp, recvbuff(ibuff, i), recvbuff(ibuff+1, i),  &
  recvbuff(ibuff+2, i), recvbuff(ibuff+3, i), recvbuff(ibuff+4, i),               &
  recvbuff(ibuff+5, i), recvbuff(ibuff+6, i), recvbuff(ibuff+7:ibuff+6+npid, i))
END DO
ispec=ispec+nvar*npart_recv(ispecies, i)
END DO
END DO


DEALLOCATE(sendbuff, recvbuff, nptoexch, npart_send, npart_recv, requests)

END SUBROUTINE remap_particles_2D


! ________________________________________________________________________________________
!> @brief
!> Finds proc index using binary search algorithm
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE get_proc_interval(iproc, ic, ncmin, ncmax, ncpus)
USE iso_c_binding

IMPLICIT NONE
INTEGER(idp), INTENT(IN OUT) :: iproc
INTEGER(idp), INTENT(IN) :: ic, ncpus
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ncmin, ncmax
INTEGER(idp) :: imin, imax, imid
LOGICAL(lp)  :: is_not_found

imin=0
imax=ncpus-1
iproc=-1
is_not_found=.TRUE.
DO WHILE((is_not_found) .AND. (imax .GE. imin) )
imid=(imax+imin)/2
IF((ic .GE. ncmin(imid)) .AND. (ic .LE. ncmax(imid))) THEN
is_not_found=.FALSE.
iproc=imid
ELSE
IF(ic .LT. ncmin(imid)) THEN
  imax=imid-1
ELSE
  imin=imid+1
ENDIF
ENDIF
END DO

END SUBROUTINE get_proc_interval

! ________________________________________________________________________________________
!> @brief
!> Binary search. This subroutine needs a better description.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2016
! ________________________________________________________________________________________
SUBROUTINE binary_search(isend, crank, arr, narr)

IMPLICIT NONE
INTEGER(idp), INTENT(IN OUT) :: isend
INTEGER(isp), INTENT(IN) :: crank
INTEGER(idp), INTENT(IN) :: narr
INTEGER(isp), INTENT(IN), DIMENSION(narr) :: arr
INTEGER(idp) :: imin, imax, imid
LOGICAL(lp)  :: is_not_found

imin=1
imax=narr
isend=-1
is_not_found=.TRUE.
DO WHILE((is_not_found) .AND. (imax .GE. imin) )
imid=(imax+imin)/2
IF(arr(imid) .EQ. crank) THEN
is_not_found=.FALSE.
isend=imid
ELSE
IF(arr(imid) .GE. crank) THEN
  imax=imid-1
ELSE
  imin=imid+1
ENDIF
ENDIF
END DO

END SUBROUTINE binary_search

END MODULE load_balance

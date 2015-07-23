!!!! --- MODULE COUNTAINING ROUTINE FOR BOUNDARY CONDITIONS ON FIELDS, CURRENTS AND PARTICLES
!!!! --- For the moment this module handles:
!!!! ---  periodic external boundary conditions for particles and fields
MODULE boundary

  USE shared_data
  USE fields
  USE particles
  USE mpi_subtype_control
  USE constants

  IMPLICIT NONE

CONTAINS

!!! --- Exchange field values at processor boundaries and apply field
!!! --- boundary conditions
  SUBROUTINE field_bc(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(-ng:,-ng:,-ng:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, ng, nx, ny, nz)

  END SUBROUTINE field_bc

!!! --- Routine for exchanging guard regions between subdomains
  SUBROUTINE do_field_mpi_with_lengths(field, ng, nx_local, ny_local, &
      nz_local)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(-ng:,-ng:,-ng:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    sizes(1) = nx_local + 1 + 2 * ng
    sizes(2) = ny_local + 1 + 2 * ng
    sizes(3) = nz_local + 1 + 2 * ng
    starts = 1

    szmax = sizes(1) * sizes(2) * ng
    sz = sizes(1) * sizes(3) * ng
    IF (sz .GT. szmax) szmax = sz
    sz = sizes(2) * sizes(3) * ng
    IF (sz .GT. szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(0,-ng,-ng), 1, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (proc_x_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -ng, subsizes(3)-ng-1
      DO j = -ng, subsizes(2)-ng-1
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(nx_local+1-ng,-ng,-ng), 1, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

    IF (proc_x_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -ng, subsizes(3)-ng-1
      DO j = -ng, subsizes(2)-ng-1
      DO i = -ng, subsizes(1)-ng-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(-ng,0,-ng), 1, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

    IF (proc_y_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -ng, subsizes(3)-ng-1
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = -ng, subsizes(1)-ng-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-ng,ny_local+1-ng,-ng), 1, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

    IF (proc_y_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -ng, subsizes(3)-ng-1
      DO j = -ng, subsizes(2)-ng-1
      DO i = -ng, subsizes(1)-ng-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_subtype(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(-ng,-ng,0), 1, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

    IF (proc_z_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = -ng, subsizes(2)-ng-1
      DO i = -ng, subsizes(1)-ng-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-ng,-ng,nz_local+1-ng), 1, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

    IF (proc_z_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -ng, subsizes(3)-ng-1
      DO j = -ng, subsizes(2)-ng-1
      DO i = -ng, subsizes(1)-ng-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE do_field_mpi_with_lengths

!!! --- Routine for adding current contributions fron adjacent subdomains
  SUBROUTINE processor_summation_bcs(array, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(-ng:,-ng:,-ng:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, nn, sz, i

    sizes(1) = nx + 1 + 2 * ng
    sizes(2) = ny + 1 + 2 * ng
    sizes(3) = nz + 1 + 2 * ng
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = ng
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn+1,-ng,-ng), 1, subarray, &
        neighbour( 1,0,0), tag, temp, sz, mpireal, &
        neighbour(-1,0,0), tag, comm, status, errcode)
    array(0:ng-1,:,:) = array(0:ng-1,:,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-ng,-ng,-ng), 1, subarray, &
        neighbour(-1,0,0), tag, temp, sz, mpireal, &
        neighbour( 1,0,0), tag, comm, status, errcode)
    array(nn+1-ng:nn,:,:) = array(nn+1-ng:nn,:,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = ng
    subsizes(3) = sizes(3)
    nn = ny

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-ng,nn+1,-ng), 1, subarray, &
        neighbour(0, 1,0), tag, temp, sz, mpireal, &
        neighbour(0,-1,0), tag, comm, status, errcode)
    array(:,0:ng-1,:) = array(:,0:ng-1,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-ng,-ng,-ng), 1, subarray, &
        neighbour(0,-1,0), tag, temp, sz, mpireal, &
        neighbour(0, 1,0), tag, comm, status, errcode)
    array(:,nn+1-ng:nn,:) = array(:,nn+1-ng:nn,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = ng
    nn = nz

    subarray = create_3d_array_subtype(mpireal, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-ng,-ng,nn+1), 1, subarray, &
        neighbour(0,0, 1), tag, temp, sz, mpireal, &
        neighbour(0,0,-1), tag, comm, status, errcode)
    array(:,:,0:ng-1) = array(:,:,0:ng-1) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-ng,-ng,-ng), 1, subarray, &
        neighbour(0,0,-1), tag, temp, sz, mpireal, &
        neighbour(0,0, 1), tag, comm, status, errcode)
    array(:,:,nn+1-ng:nn) = array(:,:,nn+1-ng:nn) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array, ng)

  END SUBROUTINE processor_summation_bcs

!!! --- Boundary condition routine for electric field
  SUBROUTINE efield_bcs

    INTEGER :: i, ng

    ng=nxguards !we assume nxguards=nyguards=nzguards
    ! These are the MPI boundaries
    CALL field_bc(ex, ng)
    CALL field_bc(ey, ng)
    CALL field_bc(ez, ng)

  END SUBROUTINE efield_bcs

!!! --- Boundary condition routine for magnetic field
  SUBROUTINE bfield_bcs

    INTEGER :: i, ng

    ng=nxguards !we assume nxguards=nyguards=nzguards
    ! These are the MPI boundaries
    CALL field_bc(bx, ng)
    CALL field_bc(by, ng)
    CALL field_bc(bz, ng)

  END SUBROUTINE bfield_bcs

!!! --- Boundary conditions routine for currents
  SUBROUTINE current_bcs

    ! domain is decomposed. We just add currents at subdomains borders
    ! Here we also assume nxguards=nyguards=nzguards

    CALL processor_summation_bcs(jx, nxguards)
    CALL processor_summation_bcs(jy, nxguards)
    CALL processor_summation_bcs(jz, nxguards)

  END SUBROUTINE current_bcs


!!! Boundary condition routine on electrons
  SUBROUTINE particle_bcs_e
    INTEGER, PARAMETER :: nvar=7 ! Simple implementation
    INTEGER, DIMENSION(-1:1,-1:1,-1:1) :: nptoexch
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: temp
    LOGICAL, DIMENSION(nparte) :: mask
    INTEGER :: ibuff, isend, nout, nbuff
    INTEGER :: xbd, ybd, zbd
    INTEGER :: ixp, iyp, izp
    INTEGER :: nsend_buf, nrecv_buf
    INTEGER :: dest, src
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, i, ip, ix, iy, iz
    REAL(num) :: part_pos
    part_pos=0.
    nptoexch=0
    nsend_buf=0
    nout=0
    nrecv_buf=0
    mask=.TRUE.
    xbd = 0
    ybd = 0
    zbd = 0
    nbuff=nparte*nvar
    ibuff=1
    ALLOCATE(sendbuf(-1:1,-1:1,-1:1,1:nbuff))
    sendbuf=0.
    ! Identify destination of particles to send and pack them
    DO i = 1, nparte
        xbd = 0
        ybd = 0
        zbd = 0
        out_of_bounds = .FALSE.

         part_pos = xe(i)
          ! Particle has left this processor
          IF (part_pos .LT. x_min_local) THEN
            xbd = -1
            IF (x_min_boundary) THEN
                xbd = -1
                xe(i) = part_pos + length_x
            ENDIF
          ENDIF

          ! Particle has left this processor
          IF (part_pos .GE. x_max_local) THEN
            xbd = 1
            IF (x_max_boundary) THEN
                xbd = 1
                xe(i) = part_pos - length_x
            ENDIF
          ENDIF

         part_pos = ye(i)
         ! Particle has left this processor
          IF (part_pos .LT. y_min_local) THEN
            ybd = -1
            IF (y_min_boundary) THEN
                ye(i) = part_pos + length_y
            ENDIF
          ENDIF


          ! Particle has left this processor
          IF (part_pos .GE. y_max_local) THEN
                ybd = 1
            IF (y_max_boundary) THEN
                ye(i) = part_pos - length_y
            ENDIF
          ENDIF

          part_pos = ze(i)
          ! Particle has left this processor
          IF (part_pos .LT. z_min_local) THEN
            zbd = -1
            IF (z_min_boundary) THEN
                ze(i) = part_pos + length_z
            ENDIF
          ENDIF

                 ! Particle has left this processor
          IF (part_pos .GE. z_max_local) THEN
            zbd = 1
            ! Particle has left the system
            IF (z_max_boundary) THEN
                zbd = 1
                ze(i) = part_pos - length_z
            ENDIF
          ENDIF

        IF ((ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) .AND. (nproc .GT. 1)) THEN
          ! Particle has left processor, send it to its neighbour
          mask(i)=.FALSE.
          nout=nout+1
          nptoexch(xbd,ybd,zbd) = nptoexch(xbd,ybd,zbd)+1
          sendbuf(xbd,ybd,zbd,ibuff)    = xe(i)
          sendbuf(xbd,ybd,zbd,ibuff+1)  = ye(i)
          sendbuf(xbd,ybd,zbd,ibuff+2)  = ze(i)
          sendbuf(xbd,ybd,zbd,ibuff+3)  = uxe(i)
          sendbuf(xbd,ybd,zbd,ibuff+4)  = uye(i)
          sendbuf(xbd,ybd,zbd,ibuff+5)  = uze(i)
          sendbuf(xbd,ybd,zbd,ibuff+6)  = we(i)
          ibuff=ibuff+nvar
        ENDIF
      ENDDO
      ! REMOVE OUTBOUND PARTICLES FROM ARRAYS
      ! update positions and velocity arrays (fields are re-calculated)
      IF (nproc .GT. 1) THEN
        IF (nout .GT. 0) THEN
            DO i = 1, nparte
                IF (.NOT. mask(i)) THEN
                    xe(i)=xe(nparte)
                    ye(i)=ye(nparte)
                    ze(i)=ze(nparte)
                    uxe(i)=uxe(nparte)
                    uye(i)=uye(nparte)
                    uze(i)=uze(nparte)
                    we(i) =we(nparte)
                    nparte=nparte-1
                ENDIF
            ENDDO
        ENDIF
      	! swap Particles
      	DO iz = -1, 1
        	DO iy = -1, 1
          		DO ix = -1, 1
            			IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
            			ixp = -ix
            			iyp = -iy
            			izp = -iz

            			! SEND - RECEIVE PARTICLES IN BUFFERS
            			!- Get number of particles in recvbuff
            			nsend_buf=nptoexch(ix,iy,iz)*nvar
            			nrecv_buf=0
            			dest = neighbour(ix,iy,iz)
            			src  = neighbour(ixp,iyp,izp)
            			CALL MPI_SENDRECV(nsend_buf, 1, MPI_INTEGER, dest, tag, nrecv_buf, 1, &
            			MPI_INTEGER, src, tag, comm, status, errcode)
            			ALLOCATE(recvbuf(nrecv_buf))

            			!- Send/receive particles to/from neighbour
            			CALL MPI_SENDRECV(sendbuf(ix,iy,iz,1:nsend_buf), nsend_buf, mpireal, dest, tag, &
            			recvbuf, nrecv_buf, mpireal, src, tag, comm, status, errcode)

            			! Add received particles to particle arrays
            			ip=0
            			DO i =1, nrecv_buf, nvar
            				ip=ip+1
            				xe(nparte+ip)=recvbuf(i)
            				ye(nparte+ip)=recvbuf(i+1)
            				ze(nparte+ip)=recvbuf(i+2)
            				uxe(nparte+ip)=recvbuf(i+3)
            				uye(nparte+ip)=recvbuf(i+4)
            				uze(nparte+ip)=recvbuf(i+5)
                            we(nparte+ip)=recvbuf(i+6)
            			END DO
            			DEALLOCATE(recvbuf)
           			   ! Increase final number of particles
            			nparte=nparte+ip
          		ENDDO
        	ENDDO
      	ENDDO
      ENDIF
      DEALLOCATE(sendbuf)
    END SUBROUTINE particle_bcs_e

!!! Boundary condition routine on protons
  SUBROUTINE particle_bcs_p
    INTEGER, PARAMETER :: nvar=7 ! Simple implementation
    INTEGER, DIMENSION(-1:1,-1:1,-1:1) :: nptoexch
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: temp
    LOGICAL, DIMENSION(npartp) :: mask
    INTEGER :: ibuff, isend, nout
    INTEGER :: xbd, ybd, zbd
    INTEGER :: ixp, iyp, izp
    INTEGER :: nsend_buf, nrecv_buf
    INTEGER :: dest, src
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, i, ip, ix, iy, iz
    REAL(num) :: part_pos

   nptoexch=0
   nsend_buf=0
   nout=0
   nrecv_buf=0
   mask=.TRUE.
   part_pos=0.
   ibuff=1
   ALLOCATE(sendbuf(-1:1,-1:1,-1:1,1:npartp*nvar))
    ! Identify destination of particles to send and pack them
    DO i = 1, npartp
        xbd = 0
        ybd = 0
        zbd = 0
        out_of_bounds = .FALSE.

        part_pos = xp(i)
          ! Particle has left this processor
          IF (part_pos .LT. x_min_local) THEN
            xbd = -1
            IF (x_min_boundary) THEN
                xbd = -1
                xp(i) = part_pos + length_x
            ENDIF
          ENDIF

          ! Particle has left this processor
          IF (part_pos .GT. x_max_local) THEN
            xbd = 1
            IF (x_max_boundary) THEN
                xbd = 1
                xp(i) = part_pos - length_x
            ENDIF
          ENDIF

        part_pos = yp(i)
         ! Particle has left this processor
          IF (part_pos .LT. y_min_local) THEN
            ybd = -1
            IF (y_min_boundary) THEN
                yp(i) = part_pos + length_y
            ENDIF
          ENDIF


          ! Particle has left this processor
          IF (part_pos .GT. y_max_local) THEN
                ybd = 1
            IF (y_max_boundary) THEN
                yp(i) = part_pos - length_y
            ENDIF
          ENDIF

        part_pos = zp(i)
          ! Particle has left this processor
          IF (part_pos .LT. z_min_local) THEN
            zbd = -1
            IF (z_min_boundary) THEN
                zp(i) = part_pos + length_z
            ENDIF
          ENDIF

                 ! Particle has left this processor
          IF (part_pos .GT. z_max_local) THEN
            zbd = 1
            ! Particle has left the system
            IF (z_max_boundary) THEN
                zbd = 1
                zp(i) = part_pos - length_z
            ENDIF
          ENDIF

        IF ((ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) .AND. (nproc .GT. 1)) THEN
          ! Particle has left processor, send it to its neighbour
          mask(i)=.FALSE.
          nout=nout+1
          nptoexch(xbd,ybd,zbd) = nptoexch(xbd,ybd,zbd)+1
          sendbuf(xbd,ybd,zbd,ibuff)    = xp(i)
          sendbuf(xbd,ybd,zbd,ibuff+1)  = yp(i)
          sendbuf(xbd,ybd,zbd,ibuff+2)  = zp(i)
          sendbuf(xbd,ybd,zbd,ibuff+3)  = uxp(i)
          sendbuf(xbd,ybd,zbd,ibuff+4)  = uyp(i)
          sendbuf(xbd,ybd,zbd,ibuff+5)  = uzp(i)
          sendbuf(xbd,ybd,zbd,ibuff+6)  = wp(i)
          ibuff=ibuff+nvar
        ENDIF

      ENDDO
      ! REMOVE OUTBOUND PARTICLES FROM ARRAYS
      ! update positions and velocity arrays (fields are re-calculated)
      IF (nproc .GT. 1) THEN
        IF (nout .GT. 0) THEN
            DO i = 1, npartp
                IF (.NOT. mask(i)) THEN
                    xp(i)=xp(npartp)
                    yp(i)=yp(npartp)
                    zp(i)=zp(npartp)
                    uxp(i)=uxp(npartp)
                    uyp(i)=uyp(npartp)
                    uzp(i)=uzp(npartp)
                    wp(i)=wp(npartp)
                    npartp=npartp-1;
                ENDIF
            ENDDO
        ENDIF
      	! swap Particles
     	 DO iz = -1, 1
        	DO iy = -1, 1
          		DO ix = -1, 1
            			IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
            			ixp = -ix
            			iyp = -iy
            			izp = -iz

            			! SEND- RECEIVE PARTICLES IN BUFFERS
            			!- Get number of particles in recvbuff from neighbour
            			nsend_buf=nptoexch(ix,iy,iz)*nvar
           			    nrecv_buf=0
            			dest = neighbour(ix,iy,iz)
            			src  = neighbour(ixp,iyp,izp)
            			CALL MPI_SENDRECV(nsend_buf, 1, MPI_INTEGER, dest, tag, nrecv_buf, 1, &
            			MPI_INTEGER, src, tag, comm, status, errcode)
            			ALLOCATE(recvbuf(nrecv_buf))

            			!- Send/receive particles to/from neighbour

            			CALL MPI_SENDRECV(sendbuf(ix,iy,iz,1:nsend_buf), nsend_buf, mpireal, dest, tag, &
            			recvbuf, nrecv_buf, mpireal, src, tag, comm, status, errcode)

            			! Add received particles to particle arrays
           		        ip=0
           		        DO i =1, nrecv_buf, nvar
           			        ip=ip+1
                            xp(npartp+ip)=recvbuf(i)
                            yp(npartp+ip)=recvbuf(i+1)
                            zp(npartp+ip)=recvbuf(i+2)
            				uxp(npartp+ip)=recvbuf(i+3)
            				uyp(npartp+ip)=recvbuf(i+4)
            				uzp(npartp+ip)=recvbuf(i+5)
                            wp(npartp+ip)=recvbuf(i+6)
            			END DO
            			DEALLOCATE(recvbuf)
            			! Increase final number of particles
                        npartp=npartp+ip
          		ENDDO
        	ENDDO
      	ENDDO
      ENDIF
      DEALLOCATE(sendbuf)
  END SUBROUTINE particle_bcs_p


END MODULE boundary

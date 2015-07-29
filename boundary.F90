!!!! --- MODULE COUNTAINING ROUTINE FOR BOUNDARY CONDITIONS ON FIELDS, CURRENTS AND PARTICLES
!!!! --- For the moment this module handles:
!!!! ---  periodic external boundary conditions for particles and fields
MODULE boundary

  USE shared_data
  USE fields
  USE particles
  USE mpi_derived_types
  USE constants

  IMPLICIT NONE

CONTAINS

!!! --- Exchange field values at processor boundaries and apply field
!!! --- boundary conditions
  SUBROUTINE field_bc(field, nxg, nyg, nzg)

    INTEGER, INTENT(IN) :: nxg, nyg, nzg
    REAL(num), DIMENSION(-nxg:,-nyg:,-nzg:), INTENT(INOUT) :: field

    CALL exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, nx, ny, nz)

  END SUBROUTINE field_bc

!!! --- Routine for exchanging guard regions between subdomains
  SUBROUTINE exchange_mpi_3d_grid_array_with_guards(field, nxg, nyg, nzg, &
             nx_local, ny_local, nz_local)

    INTEGER, INTENT(IN) :: nxg, nyg, nzg
    REAL(num), DIMENSION(-nxg:,-nyg:,-nzg:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpidbl

    sizes(1) = nx_local + 1 + 2 * nxg
    sizes(2) = ny_local + 1 + 2 * nyg
    sizes(3) = nz_local + 1 + 2 * nzg
    starts = 1

    szmax = sizes(1) * sizes(2) * nzg
    sz = sizes(1) * sizes(3) * nyg
    IF (sz .GT. szmax) szmax = sz
    sz = sizes(2) * sizes(3) * nxg
    IF (sz .GT. szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(0,-nyg,-nzg), 1, subarray, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (proc_x_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = nx_local+1, subsizes(1)+nx_local
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(nx_local+1-nxg,-nyg,-nzg), 1, subarray, proc_x_max, &
        tag, temp, sz, basetype, proc_x_min, tag, comm, status, errcode)

    IF (proc_x_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(-nxg,0,-nzg), 1, subarray, proc_y_min, &
        tag, temp, sz, basetype, proc_y_max, tag, comm, status, errcode)

    IF (proc_y_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = ny_local+1, subsizes(2)+ny_local
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-nxg,ny_local+1-nyg,-nzg), 1, subarray, proc_y_max, &
        tag, temp, sz, basetype, proc_y_min, tag, comm, status, errcode)

    IF (proc_y_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg

    sz = subsizes(1) * subsizes(2) * subsizes(3)

    subarray = create_3d_array_derived_type(basetype, subsizes, sizes, starts)

    CALL MPI_SENDRECV(field(-nxg,-nyg,0), 1, subarray, proc_z_min, &
        tag, temp, sz, basetype, proc_z_max, tag, comm, status, errcode)

    IF (proc_z_max .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = nz_local+1, subsizes(3)+nz_local
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_SENDRECV(field(-nxg,-nyg,nz_local+1-nzg), 1, subarray, proc_z_max, &
        tag, temp, sz, basetype, proc_z_min, tag, comm, status, errcode)

    IF (proc_z_min .NE. MPI_PROC_NULL) THEN
      n = 1
      DO k = -nzg, subsizes(3)-nzg-1
      DO j = -nyg, subsizes(2)-nyg-1
      DO i = -nxg, subsizes(1)-nxg-1
        field(i,j,k) = temp(n)
        n = n + 1
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    CALL MPI_TYPE_FREE(subarray, errcode)

    DEALLOCATE(temp)

  END SUBROUTINE exchange_mpi_3d_grid_array_with_guards

!!! --- Routine for adding current contributions fron adjacent subdomains
  SUBROUTINE summation_bcs(array, nxg, nyg, nzg)

    INTEGER, INTENT(IN) :: nxg, nyg, nzg
    REAL(num), DIMENSION(-nxg:,-nyg:,-nzg:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, nn, sz, i

    sizes(1) = nx + 1 + 2 * nxg
    sizes(2) = ny + 1 + 2 * nyg
    sizes(3) = nz + 1 + 2 * nzg
    starts = 1

    !! -- Summation along X- direction
    subsizes(1) = nxg
    subsizes(2) = sizes(2)
    subsizes(3) = sizes(3)
    nn = nx

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nn+1,-nyg,-nzg), 1, subarray, &
        neighbour( 1,0,0), tag, temp, sz, mpidbl, &
        neighbour(-1,0,0), tag, comm, status, errcode)
    array(0:nxg-1,:,:) = array(0:nxg-1,:,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg), 1, subarray, &
        neighbour(-1,0,0), tag, temp, sz, mpidbl, &
        neighbour( 1,0,0), tag, comm, status, errcode)
    array(nn+1-nxg:nn,:,:) = array(nn+1-nxg:nn,:,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Y- direction
    subsizes(1) = sizes(1)
    subsizes(2) = nyg
    subsizes(3) = sizes(3)
    nn = ny

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,nn+1,-nzg), 1, subarray, &
        neighbour(0, 1,0), tag, temp, sz, mpidbl, &
        neighbour(0,-1,0), tag, comm, status, errcode)
    array(:,0:nyg-1,:) = array(:,0:nyg-1,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg), 1, subarray, &
        neighbour(0,-1,0), tag, temp, sz, mpidbl, &
        neighbour(0, 1,0), tag, comm, status, errcode)
    array(:,nn+1-nyg:nn,:) = array(:,nn+1-nyg:nn,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    !! -- Summation along Z- direction
    subsizes(1) = sizes(1)
    subsizes(2) = sizes(2)
    subsizes(3) = nzg
    nn = nz

    subarray = create_3d_array_derived_type(mpidbl, subsizes, sizes, starts)

    sz = subsizes(1) * subsizes(2) * subsizes(3)
    ALLOCATE(temp(subsizes(1), subsizes(2), subsizes(3)))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,nn+1), 1, subarray, &
        neighbour(0,0, 1), tag, temp, sz, mpidbl, &
        neighbour(0,0,-1), tag, comm, status, errcode)
    array(:,:,0:nzg-1) = array(:,:,0:nzg-1) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-nxg,-nyg,-nzg), 1, subarray, &
        neighbour(0,0,-1), tag, temp, sz, mpidbl, &
        neighbour(0,0, 1), tag, comm, status, errcode)
    array(:,:,nn+1-nzg:nn) = array(:,:,nn+1-nzg:nn) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array, nxg, nyg, nzg)

  END SUBROUTINE summation_bcs

!!! --- Boundary condition routine for electric field
  SUBROUTINE efield_bcs
    ! Electric field MPI exchange between subdomains
    CALL field_bc(ex, nxguards, nyguards, nzguards)
    CALL field_bc(ey, nxguards, nyguards, nzguards)
    CALL field_bc(ez, nxguards, nyguards, nzguards)
  END SUBROUTINE efield_bcs

!!! --- Boundary condition routine for magnetic field
  SUBROUTINE bfield_bcs
    ! Magnetic field MPI exchange between subdomains
    CALL field_bc(bx, nxguards, nyguards, nzguards)
    CALL field_bc(by, nxguards, nyguards, nzguards)
    CALL field_bc(bz, nxguards, nyguards, nzguards)

  END SUBROUTINE bfield_bcs

!!! --- Boundary conditions routine for currents
  SUBROUTINE current_bcs
    ! Add current contribution from adjacent subdomains
    CALL summation_bcs(jx, nxguards, nyguards, nzguards)
    CALL summation_bcs(jy, nxguards, nyguards, nzguards)
    CALL summation_bcs(jz, nxguards, nyguards, nzguards)
  END SUBROUTINE current_bcs

!!! --- Boundary conditions routine for charge density
SUBROUTINE charge_bcs
! Add charge contribution from adjacent subdomains
    CALL summation_bcs(rho, nxguards, nyguards, nzguards)
END SUBROUTINE charge_bcs


!!! Boundary condition routine on particles
  SUBROUTINE particle_bcs
    INTEGER, PARAMETER :: nvar=7 ! Simple implementation
    INTEGER, DIMENSION(-1:1,-1:1,-1:1) :: nptoexch
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:) :: sendbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: recvbuf
    REAL(num), ALLOCATABLE, DIMENSION(:) :: temp
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask
    INTEGER :: ibuff, isend, nout, nbuff
    INTEGER :: xbd, ybd, zbd
    INTEGER :: ixp, iyp, izp
    INTEGER :: nsend_buf, nrecv_buf, npart_curr
    INTEGER :: dest, src
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, i, ip, ix, iy, iz
    REAL(num) :: part_xyz
    TYPE(particle_species), POINTER :: curr
    ! Identify destination of particles to send and pack them
    DO ispecies=1, nspecies
        ! Init send recv buffers
        curr => species_parray(ispecies)
        part_xyz=0.
        nptoexch=0
        nsend_buf=0
        nout=0
        nrecv_buf=0
        xbd = 0
        ybd = 0
        zbd = 0
        nbuff=curr%species_npart*nvar
        ibuff=1
        ALLOCATE(sendbuf(-1:1,-1:1,-1:1,1:nbuff))
        sendbuf=0.
        ALLOCATE(mask(1:curr%species_npart))
        mask=.TRUE.
        ! Identify outbounds particles
        DO i = 1, curr%species_npart
            xbd = 0
            ybd = 0
            zbd = 0
            out_of_bounds = .FALSE.
            part_xyz = curr%part_x(i)
            ! Particle has left this processor
            IF (part_xyz .LT. x_min_local) THEN
                xbd = -1
                IF (x_min_boundary) THEN
                    xbd = -1
                    curr%part_x(i) = part_xyz + length_x
                ENDIF
            ENDIF

            ! Particle has left this processor
            IF (part_xyz .GE. x_max_local) THEN
                xbd = 1
                IF (x_max_boundary) THEN
                    xbd = 1
                    curr%part_x(i) = part_xyz - length_x
                ENDIF
            ENDIF

            part_xyz = curr%part_y(i)
            ! Particle has left this processor
            IF (part_xyz .LT. y_min_local) THEN
                ybd = -1
                IF (y_min_boundary) THEN
                    curr%part_y(i) = part_xyz + length_y
                ENDIF
            ENDIF


            ! Particle has left this processor
            IF (part_xyz .GE. y_max_local) THEN
                ybd = 1
                IF (y_max_boundary) THEN
                    curr%part_y(i) = part_xyz - length_y
                ENDIF
            ENDIF

            part_xyz = curr%part_z(i)
            ! Particle has left this processor
            IF (part_xyz .LT. z_min_local) THEN
                zbd = -1
                IF (z_min_boundary) THEN
                    curr%part_z(i) = part_xyz + length_z
                ENDIF
            ENDIF

            ! Particle has left this processor
            IF (part_xyz .GE. z_max_local) THEN
                zbd = 1
                ! Particle has left the system
                IF (z_max_boundary) THEN
                    zbd = 1
                    curr%part_z(i) = part_xyz - length_z
                ENDIF
            ENDIF

            IF ((ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) .AND. (nproc .GT. 1)) THEN
            ! Particle has left processor, send it to its neighbour
                mask(i)=.FALSE.
                nout=nout+1
                nptoexch(xbd,ybd,zbd) = nptoexch(xbd,ybd,zbd)+1
                sendbuf(xbd,ybd,zbd,ibuff)    = curr%part_x(i)
                sendbuf(xbd,ybd,zbd,ibuff+1)  = curr%part_y(i)
                sendbuf(xbd,ybd,zbd,ibuff+2)  = curr%part_z(i)
                sendbuf(xbd,ybd,zbd,ibuff+3)  = curr%part_ux(i)
                sendbuf(xbd,ybd,zbd,ibuff+4)  = curr%part_uy(i)
                sendbuf(xbd,ybd,zbd,ibuff+5)  = curr%part_uz(i)
                sendbuf(xbd,ybd,zbd,ibuff+6)  = curr%weight(i)
                ibuff=ibuff+nvar
            ENDIF
        ENDDO
        ! REMOVE OUTBOUND PARTICLES FROM ARRAYS
        ! update positions and velocity arrays (fields are re-calculated)
        IF (nproc .GT. 1) THEN
            IF (nout .GT. 0) THEN
                DO i = 1, curr%species_npart
                    IF (.NOT. mask(i)) THEN
                        curr%part_x(i)=curr%part_x(curr%species_npart)
                        curr%part_y(i)=curr%part_y(curr%species_npart)
                        curr%part_z(i)=curr%part_z(curr%species_npart)
                        curr%part_ux(i)=curr%part_ux(curr%species_npart)
                        curr%part_uy(i)=curr%part_uy(curr%species_npart)
                        curr%part_uz(i)=curr%part_uz(curr%species_npart)
                        curr%weight(i)=curr%weight(curr%species_npart)
                        curr%species_npart=curr%species_npart-1
                    ENDIF
                ENDDO
            ENDIF
            ! SEND/RECEIVE PARTICLES TO/FROM ADJACENT SUBDOMAINS
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
                            CALL MPI_SENDRECV(sendbuf(ix,iy,iz,1:nsend_buf), nsend_buf, mpidbl, dest, tag, &
                            recvbuf, nrecv_buf, mpidbl, src, tag, comm, status, errcode)

                            ! Add received particles to particle arrays
                            ip=0
                            DO i =1, nrecv_buf, nvar
                                ip=ip+1
                                curr%part_x(curr%species_npart+ip)=recvbuf(i)
                                curr%part_y(curr%species_npart+ip)=recvbuf(i+1)
                                curr%part_z(curr%species_npart+ip)=recvbuf(i+2)
                                curr%part_ux(curr%species_npart+ip)=recvbuf(i+3)
                                curr%part_uy(curr%species_npart+ip)=recvbuf(i+4)
                                curr%part_uz(curr%species_npart+ip)=recvbuf(i+5)
                                curr%weight(curr%species_npart+ip)=recvbuf(i+6)
                            END DO
                            DEALLOCATE(recvbuf)
                            ! Increase final number of particles
                            curr%species_npart=curr%species_npart+ip
                    ENDDO
                ENDDO
            ENDDO
          ENDIF
        DEALLOCATE(sendbuf)
        DEALLOCATE(mask)
    END DO ! End loop on species
  END SUBROUTINE particle_bcs

END MODULE boundary

MODULE load_balance 
USE fields 
USE shared_data
USE tiling 
IMPLICIT NONE 

CONTAINS 

!SUBROUTINE mpi_load_balance_fields 






!END SUBROUTINE mpi_load_balance_fields 




! This subroutine remaps emfield_old in emfield_new and 
! takes care of all MPI exchanges between different MPI_PROCESSES
SUBROUTINE remap_em_3Dfields(emfield_old,nxold,nyold,nzold,                 &
                           ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,     &
                           emfield_new,nxnew,nynew,nznew,                 &
                           ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,     &
                           iproc, nprocs, communicator, ierrcode)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxold, nyold, nzold, nxnew, nynew,nznew, iproc, nprocs
    INTEGER(isp), INTENT(IN) :: communicator
    INTEGER(isp), INTENT(IN OUT) ::  ierrcode
    REAL(num), INTENT(IN), DIMENSION(nxold,nyold,nzold) :: emfield_old
    REAL(num), INTENT(IN  OUT), DIMENSION(nxnew,nynew,nznew) :: emfield_new
    INTEGER, INTENT(IN), DIMENSION(nprocs) :: ix1old, ix2old, iy1old, iy2old, iz1old, iz2old
    INTEGER, INTENT(IN), DIMENSION(nprocs) :: ix1new, ix2new, iy1new, iy2new, iz1new, iz2new
    INTEGER(idp) :: ix3min,ix3max,iy3min,iy3max,iz3min,iz3max
    INTEGER(idp) :: ix1newip, ix2newip, iy1newip, iy2newip, iz1newip, iz2newip
    INTEGER(idp) :: ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip, iz2oldip
    INTEGER(idp) :: ixmin_old, ixmax_old, iymin_old, iymax_old, izmin_old, izmax_old
    INTEGER(idp) :: ixmin_new, ixmax_new, iymin_new, iymax_new, izmin_new, izmax_new
    LOGICAL(idp) :: l_is_intersection
    INTEGER(isp), DIMENSION(nprocs) :: sendtype, recvtype 
    INTEGER(isp), DIMENSION(2_isp*nprocs) :: requests 
    INTEGER(isp) :: mpitag, proc_rank,  i, nsreq, nrreq, error
    INTEGER(isp), PARAMETER :: nd=3
    INTEGER(isp), DIMENSION(nd) :: nsub, nglob, start
    
    nglob(1) = nxnew
    nglob(2) = nynew
    nglob(3) = nznew
    
    ! ------ DATA TO BE RECEIVED BY OTHER PROCS 
    ! Computes intersection between new proc limit and old adjacent procs limits 
    ! Recvtypes are processed in this part 
    ix1newip = ix1new(iproc)
    ix2newip = ix2new(iproc)
    iy1newip = iy1new(iproc)
    iy2newip = iy2new(iproc)
    iz1newip = iz1new(iproc)
    iz2newip = iz2new(iproc)
    nrreq=0
    DO i=1, nprocs 
        CALL get_3Dintersection(ix1newip, ix2newip, iy1newip, iy2newip,         &
                                iz1newip, iz2newip,                             & 
                                ix1old(i), ix2old(i), iy1old(i), iy2old(i),     &
                                iz1old(i), iz2old(i),                           & 
                                ix3min,ix3max,iy3min,iy3max,iz3min,iz3max,      &          
                                l_is_intersection)
                                
        ! If i == iproc just do a copy of emfield_old in emfield_new
        IF ((i .EQ. iproc) .AND. l_is_intersection) THEN  
            ixmin_old = ix3min - ix1old(i) + 1 ; ixmax_old = ix3max - ix1old(i) +1
            iymin_old = iy3min - iy1old(i) + 1 ; iymax_old = iy3max - iy1old(i) + 1 
            izmin_old = iz3min - iz1old(i) + 1 ; izmax_old = iz3max - iz1old(i) + 1 
            ixmin_new = ix3min - ix1newip + 1 ; ixmax_new = ix3max - ix1newip + 1
            iymin_new = iy3min - ix1newip + 1 ; iymax_new = iy3max - iy1newip + 1
            izmin_new = iz3min - iz1newip + 1; izmax_new = iz3max - iz1newip +1
            emfield_new(ixmin_new:ixmax_new,iymin_new:iymax_new,izmin_new:izmax_new) = &
            emfield_old(ixmin_old:ixmax_old,iymin_old:iymax_old,izmin_old:izmax_old)
            CYCLE
        END IF
                                
        ! Found intersection area between new proc and old adjacent proc 
        ! Post IRECV  
        IF (l_is_intersection .AND. (i .NE. iproc)) THEN 
            nrreq=nrreq+1
            !--- Create recv type 
            nsub(1)  = ix3max-ix3min+1
            nsub(2)  = iy3max-iy3min+1
            nsub(3)  = iz3max-iz3min+1
            ! Arrays assumed to start at index 0 in MPI_TYPE_CREATE
            start(1) = ix3min-ix1new(iproc)+1-1 
            start(2) = iy3min-iy1new(iproc)+1-1
            start(3) = iz3min-iz1new(iproc)+1-1
            CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob, nsub, start, MPI_ORDER_FORTRAN, &
                                         MPI_REAL8, recvtype(nrreq), ierrcode)
            !--- Post IRECV for this area 
            CALL MPI_IRECV(emfield_new, 1_isp,  recvtype(nrreq), i, mpitag,    &
                            communicator, requests(nrreq), ierrcode)
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
    nsreq=0
    DO i=1, nprocs 
        CALL get_3Dintersection(ix1oldip, ix2oldip, iy1oldip, iy2oldip,         &
                                iz1oldip, iz2oldip,                             & 
                                ix1new(i), ix2new(i), iy1new(i), iy2new(i),     &
                                iz1new(i), iz2new(i),                           & 
                                ix3min,ix3max,iy3min,iy3max,iz3min,iz3max,      &          
                                l_is_intersection)
                                
        ! Case i == iproc already treated in first loop of this subroutine  
        IF (i .EQ. iproc) CYCLE 
                                
        ! Found intersection area between old proc and new adjacent procs
        ! Post ISEND  
        IF (l_is_intersection .AND. (i .NE. iproc)) THEN 
            nsreq=nsreq+1
            !--- Create send type 
            nsub(1)  = ix3max-ix3min+1
            nsub(2)  = iy3max-iy3min+1
            nsub(3)  = iz3max-iz3min+1
            ! Arrays assumed to start at index 0 in MPI_TYPE_CREATE
            start(1) = ix3min-ix1oldip+1-1 
            start(2) = iy3min-iy1oldip+1-1
            start(3) = iz3min-iz1oldip+1-1
            CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob, nsub, start, MPI_ORDER_FORTRAN, &
                                         MPI_REAL8, sendtype(nsreq), ierrcode)
            !--- Post IRECV for this area 
            CALL MPI_IRECV(emfield_old, 1_isp,  sendtype(nsreq), i, mpitag,    &
                            communicator, requests(nrreq+nsreq), ierrcode)
        ENDIF 
    END DO   
    
    IF (nsreq .NE. nrreq) PRINT *, "HUM, SOMETHING WENT WRONG IN 3D MPI REMAPING OF FIELDS"    
    
    ! DO SOME SYNC BEFORE GOING ON  
    CALL MPI_WAITALL(nsreq+nrreq, requests, MPI_STATUSES_IGNORE, errcode)
    
    ! FREE ALL DATATYPES 
    DO i=1,nsreq
        CALL MPI_TYPE_FREE(sendtype(i), ierrcode)
        CALL MPI_TYPE_FREE(recvtype(i), ierrcode)
    END DO 
END SUBROUTINE remap_em_3Dfields


! This subroutine get intersection area between two 3D domains  
! If no intersection l_is_intersection is .FALSE.
! Useful to determine wether to send/recv datas bases on new CPU split 
SUBROUTINE get_3Dintersection(ix1min,ix1max,iy1min,iy1max,iz1min,iz1max, & 
                              ix2min,ix2max,iy2min,iy2max,iz2min,iz2max, & 
                              ix3min,ix3max,iy3min,iy3max,iz3min,iz3max, &                                 
                              l_is_intersection)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) ::  ix1min, ix1max, iy1min, iy1max, iz1min, iz1max
    INTEGER(idp), INTENT(IN) ::  ix2min, ix2max, iy2min, iy2max, iz2min, iz2max
    INTEGER(idp), INTENT(IN OUT) ::  ix3min, ix3max, iy3min, iy3max, iz3min, iz3max
    LOGICAL(idp), INTENT(IN OUT) :: l_is_intersection 
    LOGICAL(idp) :: l_is_intersectionx, l_is_intersectiony, l_is_intersectionz
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
            ix3max=MIN(ix2max,ix1max)
        ENDIF 
    ELSE 
        IF(ix1min .LE. ix2max) THEN 
            l_is_intersectionx=.TRUE.
            ix3min=ix1min
            ix3max=MIN(ix1max,ix2max)
        ENDIF
    ENDIF

    ! - Y DIRECTION
    IF (iy2min .GE. iy1min) THEN 
        IF(iy2min .LE. iy1max) THEN 
            l_is_intersectiony=.TRUE.
            iy3min=iy2min
            iy3max=MIN(iy2max,iy1max)
        ENDIF 
    ELSE 
        IF(iy1min .LE. iy2max) THEN 
            l_is_intersectiony=.TRUE.
            iy3min=iy1min
            iy3max=MIN(iy1max,iy2max)
        ENDIF
    ENDIF
    
    ! - Z DIRECTION
    IF (iz2min .GE. iz1min) THEN 
        IF(iz2min .LE. iz1max) THEN 
            l_is_intersectionz=.TRUE.
            iz3min=iz2min
            iz3max=MIN(iz2max,iz1max)
        ENDIF 
    ELSE 
        IF(iz1min .LE. iz2max) THEN 
            l_is_intersectionz=.TRUE.
            iz3min=iz1min
            iz3max=MIN(iz1max,iz2max)
        ENDIF
    ENDIF
    
    IF (l_is_intersectionx .AND. l_is_intersectiony .AND. l_is_intersectionz) l_is_intersection=.TRUE. 
        
END SUBROUTINE get_3Dintersection

SUBROUTINE remap_particles()
    IMPLICIT NONE 
    ! Remap particles between procs 
END SUBROUTINE remap_particles 


! This subroutine computes the total time per part for particle subroutines 
! (i.e  mainly particle push, field gathering and current deposition)
SUBROUTINE compute_time_per_part()
    IMPLICIT NONE 
    REAL(num) :: global_time_part
    CALL get_local_number_of_part(npart_local)
    global_time_part=0.
    ! Get max time per it 
    CALL MPI_ALLREDUCE(local_time_part, global_time_part, 1_isp, MPI_REAL8, MPI_SUM, comm, errcode)
    CALL MPI_ALLREDUCE(npart_local, npart_global, 1_isp, MPI_REAL8, MPI_SUM, comm, errcode)
    global_time_per_part=global_time_part/npart_global 

END SUBROUTINE compute_time_per_part

! This subroutine computes the total time per cell for em field subroutines 
! (i.e field pusher) 
SUBROUTINE compute_time_per_cell()
    IMPLICIT NONE 
    REAL(num) :: global_time_cell
    global_time_cell=0.
    ! Get max time per it 
    CALL MPI_ALLREDUCE(local_time_cell, global_time_cell, 1_isp, MPI_REAL8, MPI_SUM, comm, errcode)
    global_time_per_cell=global_time_cell/(nx_global*ny_global*nz_global)

END SUBROUTINE compute_time_per_cell



SUBROUTINE get_max_time_per_it()
    IMPLICIT NONE 
    ! Get max time per it 
    CALL MPI_ALLREDUCE(mpitime_per_it, max_time_per_it, 1_isp, MPI_REAL8, MPI_MAX, comm, errcode)

END SUBROUTINE get_max_time_per_it 

SUBROUTINE get_min_time_per_it()
    IMPLICIT NONE 
    ! Get max time per it 
    CALL MPI_ALLREDUCE(mpitime_per_it, min_time_per_it, 1_isp, MPI_REAL8, MPI_MIN, comm, errcode)
END SUBROUTINE get_min_time_per_it 

SUBROUTINE compute_new_split()
    IMPLICIT NONE
    REAL(num), DIMENSION(:), ALLOCATABLE :: load_on_x, load_on_y, load_on_z
    ALLOCATE(load_on_x(0:nx_global-1),load_on_y(0:ny_global-1),load_on_z(0:nz_global-1))
    load_on_x=0.
    load_on_y=0.
    load_on_z=0.
    
    ! Compute load in X and compute new split in X 
    CALL get_projected_load_on_x(nx_global,load_on_x,global_time_per_part,global_time_per_cell)
    CALL balance_in_dir(load_on_x,nx_global,nprocx,new_cell_x_min,new_cell_x_max)
    ! Compute load in Y and compute new split in Y 
    CALL get_projected_load_on_y(ny_global,load_on_y,global_time_per_part,global_time_per_cell)
    CALL balance_in_dir(load_on_y,ny_global,nprocy,new_cell_y_min,new_cell_y_max)
    ! Compute load in Z and compute new split in Z 
    CALL get_projected_load_on_z(nz_global,load_on_z,global_time_per_part,global_time_per_cell)
    CALL balance_in_dir(load_on_z,nz_global,nprocz,new_cell_z_min,new_cell_z_max)
    
    DEALLOCATE(load_on_x,load_on_y,load_on_z)
    
    END SUBROUTINE compute_new_split

! This subroutine computes new load in direction dir based on the projected 
! load_in_dir 
SUBROUTINE balance_in_dir(load_in_dir, ncellmaxdir, nproc_in_dir, idirmin, idirmax)
    IMPLICIT NONE 
    REAL(num), DIMENSION(0:ncellmaxdir-1), INTENT(IN) :: load_in_dir
    INTEGER(idp), DIMENSION(nproc_in_dir), INTENT(IN OUT) :: idirmin, idirmax
    INTEGER(idp), INTENT(IN) :: nproc_in_dir, ncellmaxdir 
    INTEGER(idp) :: iproc, icell 
    REAL(num) :: balanced_load=0_num, curr_proc_load=0_num  
    LOGICAL(idp) :: not_balanced 

    balanced_load=SUM(load_in_dir)/nproc_in_dir
    icell=0 
    DO iproc=1, nproc_in_dir
        not_balanced=.TRUE. 
        curr_proc_load=0_num
        idirmin(iproc)=icell
        DO WHILE((not_balanced) .AND. (icell .LT. ncellmaxdir))  
            curr_proc_load=curr_proc_load+load_in_dir(icell)
            icell=icell+1
            IF(curr_proc_load .GE. balanced_load) not_balanced=.FALSE. 
        END DO 
        idirmax(iproc)=icell-1
    END DO 
END SUBROUTINE balance_in_dir


! This subroutine computes the computational load projected on X-Axis 
SUBROUTINE get_projected_load_on_x(nxg,load_on_x,time_per_part,time_per_cell)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: nxg
    REAL(num), INTENT(IN OUT), DIMENSION(0:nxg-1) :: load_on_x 
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell 
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icellx
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile),  POINTER :: curr_tile

    ALLOCATE(load_part(0:nxg-1))
    ALLOCATE(load_part_sum(0:nxg-1))
    load_part=0
    load_part_sum=0

    ! Get local distribution of particles along X-axis 
    DO ispecies=1, nspecies  
        curr => species_parray(ispecies)
        DO iz=1,ntilez
            DO iy=1,ntiley
                DO ix=1,ntilex
                    curr_tile => curr%array_of_tiles(ix,iy,iz)
                    count=curr_tile%np_tile(1)
                    DO ip=1,count
                        icellx=FLOOR((curr_tile%part_x(ip)-x_grid_min)/dx)
                        load_part(icellx)=load_part(icellx)+1                
                    END DO 
                END DO
            END DO
        END DO
    END DO
    
    ! Get contributions on X-axis from other MPI domains 
    CALL MPI_ALLREDUCE(load_part, load_part_sum, INT(nxg,isp), MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Computes load_in_x
    load_on_x = load_part_sum*time_per_part + ny_global * nz_global*time_per_cell
    
    DEALLOCATE(load_part,load_part_sum)
    
END SUBROUTINE get_projected_load_on_x

! This subroutine computes the computational load projected on Y-Axis 
SUBROUTINE get_projected_load_on_y(nyg,load_on_y,time_per_part,time_per_cell)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: nyg
    REAL(num), INTENT(IN OUT), DIMENSION(0:nyg-1) :: load_on_y 
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell 
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part 
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icelly
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile),  POINTER :: curr_tile

    ALLOCATE(load_part(0:nyg-1))
    ALLOCATE(load_part_sum(0:nyg-1))
    load_part=0
    load_part_sum=0
    ! Get local distribution of particles along Y-axis 
    DO ispecies=1, nspecies  
        curr => species_parray(ispecies)
        DO iz=1,ntilez
            DO iy=1,ntiley
                DO ix=1,ntilex
                    curr_tile => curr%array_of_tiles(ix,iy,iz)
                    count=curr_tile%np_tile(1)
                    DO ip=1,count
                        icelly=FLOOR((curr_tile%part_y(ip)-y_grid_min)/dy)
                        load_part(icelly)=load_part(icelly)+1                
                    END DO 
                END DO
            END DO
        END DO
    END DO
    
    ! Get contributions on Y-axis from other MPI domains 
    CALL MPI_ALLREDUCE(load_part, load_part_sum, INT(nyg,isp), MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Computes load_in_y
    load_on_y = load_part_sum*time_per_part + nx_global * nz_global*time_per_cell
    
    
    DEALLOCATE(load_part,load_part_sum)
    
END SUBROUTINE get_projected_load_on_y

! THh subroutine computes the computational load projected on Z-Axis 
SUBROUTINE get_projected_load_on_z(nzg,load_on_z,time_per_part,time_per_cell)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: nzg
    REAL(num), INTENT(IN OUT), DIMENSION(0:nzg-1) :: load_on_z 
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell 
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part 
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icellz
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile),  POINTER :: curr_tile

    ALLOCATE(load_part(0:nzg-1))
    ALLOCATE(load_part_sum(0:nzg-1))
    load_part=0
    load_part_sum=0
    ! Get local distribution of particles along Z-axis 
    DO ispecies=1, nspecies  
        curr => species_parray(ispecies)
        DO iz=1,ntilez
            DO iy=1,ntiley
                DO ix=1,ntilex
                    curr_tile => curr%array_of_tiles(ix,iy,iz)
                    count=curr_tile%np_tile(1)
                    DO ip=1,count
                        icellz=FLOOR((curr_tile%part_z(ip)-z_grid_min)/dz)
                        load_part(icellz)=load_part(icellz)+1                
                    END DO 
                END DO
            END DO
        END DO
    END DO
    
    ! Get contributions on Z-axis from other MPI domains 
    CALL MPI_ALLREDUCE(load_part, load_part_sum, INT(nzg,isp), MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Computes load_in_z
    load_on_z = load_part_sum*time_per_part + nx_global * ny_global*time_per_cell
    
    DEALLOCATE(load_part,load_part_sum)
    
END SUBROUTINE get_projected_load_on_z


END MODULE load_balance 
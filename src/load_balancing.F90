MODULE load_balance 
USE fields 
USE shared_data
USE tiling 
IMPLICIT NONE 

CONTAINS 



! Main subroutine for load balancing 
SUBROUTINE mpi_load_balance()
    IMPLICIT NONE 
    REAL(num), PARAMETER :: imbalance_threshold = 0.15 ! 15% imbalance threshold  
    REAL(num) :: imbalance
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: ix1old, ix2old, iy1old, iy2old, iz1old, iz2old
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: ix1new, ix2new, iy1new, iy2new, iz1new, iz2new
    INTEGER(idp) :: nx_new, ny_new, nz_new
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: ex_new, ey_new, ez_new, bx_new, by_new, bz_new
    
    ALLOCATE(ix1old(0:nproc-1),ix2old(0:nproc-1),iy1old(0:nproc-1),iy2old(0:nproc-1), &
            iz1old(0:nproc-1),iz2old(0:nproc-1))
    ALLOCATE(ix1new(0:nproc-1),ix2new(0:nproc-1),iy1new(0:nproc-1),iy2new(0:nproc-1), &
            iz1new(0:nproc-1),iz2new(0:nproc-1))
    
    ! Get max time per it of all processors
    CALL get_max_time_per_it()
    ! Get min time per it of all processors 
    CALL get_min_time_per_it()
    ! Compute current imbalance 
    imbalance = (max_time_per_it-min_time_per_it)/min_time_per_it
    ! Imbalance starts to be too high --> Try load balance the simulation 
    IF (imbalance .GT. 0.15_num)  THEN 
        ! Compute time per part for particle routines 
        CALL compute_time_per_part()
        ! Compute time per cell for maxwell routines 
        CALL compute_time_per_cell() 
        ! Compute new split based on projected loads on each axis
        CALL compute_new_split(global_time_per_part,global_time_per_cell,nx_global,ny_global,nz_global, &
           new_cell_x_min,new_cell_x_max,new_cell_y_min,new_cell_y_max,new_cell_z_min,new_cell_z_max,&
           nprocx,nprocy,nprocz)
        ! Get 1D array for proc limits 
        CALL get_1Darray_proclimits(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,    &
                                    cell_x_min,cell_y_min,cell_z_min,             & 
                                    cell_x_max,cell_y_max,cell_z_max,             &
                                    nprocx, nprocy, nprocz, nproc,.TRUE.)
        CALL get_1Darray_proclimits(ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,     &
                                    new_cell_x_min,new_cell_y_min,new_cell_z_min,  & 
                                    new_cell_x_max,new_cell_y_max,new_cell_z_max,  &
                                    nprocx, nprocy, nprocz, nproc,.TRUE.)
        ! Compute new dimensions 
        CALL compute_currproc_array_dimensions(nx_new,new_cell_x_min,new_cell_x_max,nproc,x_coords)
        CALL compute_currproc_array_dimensions(ny_new,new_cell_y_min,new_cell_y_max,nproc,y_coords)
        CALL compute_currproc_array_dimensions(nz_new,new_cell_z_min,new_cell_z_max,nproc,z_coords)
        
        ! Remap field Ex
        ALLOCATE(ex_new(-nxguards:nx_new+nxguards,-nyguards:ny_new+nyguards,-nzguards:nz_new+nzguards))
        CALL mpi_remap_3D_field_component(ex_new,nx_new,ny_new,nz_new,                  &
                                        ex,nx,ny,nz,                                    &
                                        nxguards,nyguards,nzguards,                     &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        rank, nproc)
        DEALLOCATE(ex)
        ex=>ex_new
        ! Remap field Ey
        ALLOCATE(ey_new(-nxguards:nx_new+nxguards,-nyguards:ny_new+nyguards,-nzguards:nz_new+nzguards))
        CALL mpi_remap_3D_field_component(ey_new,nx_new,ny_new,nz_new,                  &
                                        ey,nx,ny,nz,                                    &
                                        nxguards,nyguards,nzguards,                     &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        rank, nproc)
        DEALLOCATE(ey)
        ey=>ey_new
        ! Remap field Ez
        ALLOCATE(ez_new(-nxguards:nx_new+nxguards,-nyguards:ny_new+nyguards,-nzguards:nz_new+nzguards))
        CALL mpi_remap_3D_field_component(ez_new,nx_new,ny_new,nz_new,                  &
                                        ez,nx,ny,nz,                                    &
                                        nxguards,nyguards,nzguards,                     &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        rank, nproc)
        DEALLOCATE(ez)
        ez=>ez_new
        ! Remap field Bx
        ALLOCATE(bx_new(-nxguards:nx_new+nxguards,-nyguards:ny_new+nyguards,-nzguards:nz_new+nzguards))
        CALL mpi_remap_3D_field_component(bx_new,nx_new,ny_new,nz_new,                  &
                                        bx,nx,ny,nz,                                    &
                                        nxguards,nyguards,nzguards,                     &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        rank, nproc)
        DEALLOCATE(bx)
        bx=>bx_new
        ! Remap field By
        ALLOCATE(by_new(-nxguards:nx_new+nxguards,-nyguards:ny_new+nyguards,-nzguards:nz_new+nzguards))
        CALL mpi_remap_3D_field_component(by_new,nx_new,ny_new,nz_new,                  &
                                        by,nx,ny,nz,                                    &
                                        nxguards,nyguards,nzguards,                     &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        rank, nproc)
        DEALLOCATE(by)
        by=>by_new
        ! Remap field Bz
        ALLOCATE(bz_new(-nxguards:nx_new+nxguards,-nyguards:ny_new+nyguards,-nzguards:nz_new+nzguards))
        CALL mpi_remap_3D_field_component(bz_new,nx_new,ny_new,nz_new,                  &
                                        bz,nx,ny,nz,                                    &
                                        nxguards,nyguards,nzguards,                     &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        rank, nproc)
        DEALLOCATE(ez)
        bz=>bz_new

    ENDIF 
    
END SUBROUTINE mpi_load_balance 


SUBROUTINE compute_currproc_array_dimensions(nnew,ncmin,ncmax,np,mpi_rank)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN OUT) :: nnew
    INTEGER(idp), INTENT(IN) :: np,mpi_rank
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) :: ncmin,ncmax

    nnew=ncmax(mpi_rank)-ncmin(mpi_rank)+1

END SUBROUTINE compute_currproc_array_dimensions

SUBROUTINE get_1Darray_proclimits(ix1,ix2,iy1,iy2,iz1,iz2,cxmin,cymin,czmin, & 
                                    cxmax,cymax,czmax,npx,npy,npz,np,l_cart_comm)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: npx, npy, npz, np
    LOGICAL(idp) :: l_cart_comm
    INTEGER(idp), INTENT(IN OUT), DIMENSION(0:np-1) :: ix1, ix2, iy1, iy2, iz1, iz2
    INTEGER(idp), INTENT(IN),  DIMENSION(0:npx-1) :: cxmin, cxmax
    INTEGER(idp), INTENT(IN), DIMENSION(0:npy-1) :: cymin, cymax
    INTEGER(idp), INTENT(IN), DIMENSION(0:npz-1) :: czmin, czmax
    INTEGER(idp) :: ix, iy, iz
    INTEGER(idp) :: curr_rank
    
    
    DO iz=0,npz-1
        DO iy=0,npy-1
            DO ix=0,npx-1
                CALL pxr_convertindtoproc(comm,ix,iy,iz,npx,npy,npz,curr_rank,l_cart_comm)
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



SUBROUTINE pxr_convertindtoproc(mpi_comm_in,ix,iy,iz,npx,npy,npz,curr_rank,l_cart_comm)
    IMPLICIT NONE 
    INTEGER(isp), INTENT(IN) :: mpi_comm_in
    INTEGER(idp), INTENT(IN) :: npx, npy, npz, ix,iy,iz
    LOGICAL(idp), INTENT(IN) :: l_cart_comm
    INTEGER(idp), INTENT(IN OUT) :: curr_rank
    INTEGER(isp) :: mpi_rank, mpi_comm_isp
    INTEGER(idp) :: ixt, iyt, izt 
    mpi_comm_isp=INT(mpi_comm_in,isp)

    IF (l_cart_comm) THEN
         CALL MPI_CART_RANK(mpi_comm_isp, (/INT(iz,isp), INT(iy,isp), INT(ix,isp)/), mpi_rank, errcode) 
         curr_rank=INT(mpi_rank,idp)
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


! Remap fields based on new split 
SUBROUTINE mpi_remap_3D_field_component(field_new,nx_new,ny_new,nz_new,                 &
                                        field_old,nx_old,ny_old,nz_old,                 &
                                        nxg,nyg,nzg,                                    &
                                        ix1old, ix2old, iy1old, iy2old, iz1old, iz2old, &
                                        ix1new, ix2new, iy1new, iy2new, iz1new, iz2new, & 
                                        iproc, np)
    IMPLICIT NONE 
    REAL(num), INTENT(IN OUT), DIMENSION(-nxg:nx_new+nxg,-nyg:ny_new+nyg,-nzg:nz_new+nzg) :: field_new
    REAL(num), INTENT(IN), DIMENSION(-nxg:nx_old+nxg,-nyg:ny_old+nyg,-nzg:nz_old+nzg) :: field_old
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) ::  ix1old, ix2old, iy1old, iy2old, iz1old, iz2old
    INTEGER(idp), DIMENSION(0:np-1), INTENT(IN) ::  ix1new, ix2new, iy1new, iy2new, iz1new, iz2new
    INTEGER(idp), INTENT(IN) :: iproc, nx_new, ny_new, nz_new, nx_old, ny_old, nz_old, np, nxg, nyg, nzg
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: ex_new, ey_new, ez_new, bx_new, by_new, bz_new
 
    INTEGER(isp) :: curr_rank, ix, iy, iz    

    ! ---- MAP field__new
    CALL remap_em_3Dfields(field_old,nx_old,ny_old,nz_old,ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,  &
                           field_new,nx_new,ny_new,nz_new,nxg,nyg,nzg,                                &
                           ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,                                 &
                           iproc,np,comm,errcode)
END SUBROUTINE mpi_remap_3D_field_component


SUBROUTINE test_dealloc_from_fortran()
IMPLICIT NONE 
REAL(num), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: ex_new

ALLOCATE(ex_new(-nxguards:nx+nxguards,-nyguards:ny+nyguards,-nzguards:nz+nzguards))
DEALLOCATE(ex)
ex=>ex_new

END SUBROUTINE test_dealloc_from_fortran


! This subroutine remaps emfield_old in emfield_new and 
! takes care of all MPI exchanges between different MPI_PROCESSES
SUBROUTINE remap_em_3Dfields(emfield_old,nxold,nyold,nzold,               &
                           ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,     &
                           emfield_new,nxnew,nynew,nznew,nxg,nyg,nzg,     &
                           ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,     &
                           iproc, nprocs, communicator, ierrcode)
    IMPLICIT NONE
    INTEGER(idp), INTENT(IN) :: nxold, nyold, nzold, nxnew, nynew,nznew, iproc, nprocs
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(isp), INTENT(IN) :: communicator
    INTEGER(isp), INTENT(IN OUT) ::  ierrcode
    REAL(num), INTENT(IN), DIMENSION(-nxg:nxold+nxg,-nyg:nyold+nyg,-nzg:nzold+nzg) :: emfield_old
    REAL(num), INTENT(IN  OUT), DIMENSION(-nxg:nxnew+nxg,-nyg:nynew+nyg,-nzg:nznew+nzg) :: emfield_new
    INTEGER(idp), INTENT(IN), DIMENSION(0:nprocs-1) :: ix1old, ix2old, iy1old, iy2old, iz1old, iz2old
    INTEGER(idp), INTENT(IN), DIMENSION(0:nprocs-1) :: ix1new, ix2new, iy1new, iy2new, iz1new, iz2new
    INTEGER(idp) :: ix, iy, iz, nsubx, nsuby, nsubz
    INTEGER(idp) :: ix3min,ix3max,iy3min,iy3max,iz3min,iz3max
    INTEGER(idp) :: ix1newip, ix2newip, iy1newip, iy2newip, iz1newip, iz2newip
    INTEGER(idp) :: ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip, iz2oldip
    INTEGER(idp) :: ixmin_old, ixmax_old, iymin_old, iymax_old, izmin_old, izmax_old
    INTEGER(idp) :: ixmin_new, ixmax_new, iymin_new, iymax_new, izmin_new, izmax_new
    LOGICAL(idp) :: l_is_intersection
    INTEGER(isp), DIMENSION(0:nprocs-1) :: sendtype, recvtype 
    INTEGER(isp), DIMENSION(0:2_isp*nprocs-1) :: requests 
    INTEGER(isp) :: mpitag, proc_rank,  i, nsreq, nrreq, error, count
    INTEGER(isp), PARAMETER :: nd=3
    INTEGER(isp), DIMENSION(nd) :: nsub, nglob, nglob_old, start


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
        CALL get_3Dintersection(ix1newip, ix2newip, iy1newip, iy2newip,         &
                                iz1newip, iz2newip,                             & 
                                ix1old(i), ix2old(i), iy1old(i), iy2old(i),     &
                                iz1old(i), iz2old(i),                           & 
                                ix3min,ix3max,iy3min,iy3max,iz3min,iz3max,      &          
                                l_is_intersection)
        !PRINT *, "rank,i",rank,i, "l_is_intersection,newvsold",l_is_intersection,iz3min,iz3max
        ! If i == iproc just do a copy of emfield_old in emfield_new
        IF ((i .EQ. iproc) .AND. l_is_intersection) THEN  
            ixmin_old = ix3min - ix1old(i)  ; ixmax_old = ix3max - ix1old(i) 
            iymin_old = iy3min - iy1old(i)  ; iymax_old = iy3max - iy1old(i)  
            izmin_old = iz3min - iz1old(i)  ; izmax_old = iz3max - iz1old(i)  
            ixmin_new = ix3min - ix1newip  ; ixmax_new = ix3max - ix1newip 
            iymin_new = iy3min - iy1newip  ; iymax_new = iy3max - iy1newip 
            izmin_new = iz3min - iz1newip ; izmax_new = iz3max - iz1newip 
            emfield_new(ixmin_new:ixmax_new,iymin_new:iymax_new,izmin_new:izmax_new) = &
            emfield_old(ixmin_old:ixmax_old,iymin_old:iymax_old,izmin_old:izmax_old)
            !PRINT *, rank, i,"Simple copy", "iz3min, iz3max, iz1newip,iz2newip,izmin_new,izmax_new",iz3min, iz3max, iz1newip,iz2newip,izmin_new,izmax_new
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
            CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob, nsub, start, MPI_ORDER_FORTRAN, &
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
        CALL get_3Dintersection(ix1oldip, ix2oldip, iy1oldip, iy2oldip,         &
                                iz1oldip, iz2oldip,                             & 
                                ix1new(i), ix2new(i), iy1new(i), iy2new(i),     &
                                iz1new(i), iz2new(i),                           & 
                                ix3min,ix3max,iy3min,iy3max,iz3min,iz3max,      &          
                                l_is_intersection)
       ! PRINT *, "rank,i",rank,i, "l_is_intersection, oldvsnew",l_is_intersection,iz3min,iz3max
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
            CALL MPI_TYPE_CREATE_SUBARRAY(nd, nglob_old, nsub, start, MPI_ORDER_FORTRAN, &
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
            CALL MPI_IRECV(emfield_new(-nxg,-nyg,-nzg), 1_isp,  recvtype(i), i, mpitag,    &
                            communicator, requests(nrreq), ierrcode)
            nrreq=nrreq+1
        ENDIF
    END DO 
    
    !POST THE ISENDs if any
    nsreq=0;
    DO i=0, nprocs-1
        IF (sendtype(i) .NE. 0) THEN 
            !--- Post ISEND for this area 
            CALL MPI_ISEND(emfield_old(-nxg,-nyg,-nzg), 1_isp,  sendtype(i), i, mpitag,    &
                           communicator, requests(nrreq+nsreq), ierrcode)
            nsreq=nsreq+1
        ENDIF
    END DO 
    
   ! DO SOME SYNC BEFORE GOING ON  
    count=nsreq+nrreq
    CALL MPI_WAITALL(count,requests, MPI_STATUSES_IGNORE, errcode)
       
    ! FREE ALL DATATYPES 
    DO i=0,nprocs-1
        IF (sendtype(i) .NE. 0) CALL MPI_TYPE_FREE(sendtype(i), ierrcode)
        IF (recvtype(i) .NE. 0) CALL MPI_TYPE_FREE(recvtype(i), ierrcode)
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


! This subroutine computes the total time per part for particle subroutines 
! (i.e  mainly particle push, field gathering and current deposition)
SUBROUTINE compute_time_per_part()
    IMPLICIT NONE 
    REAL(num) :: global_time_part
    CALL get_local_number_of_part(npart_local)
    global_time_part=0.
    ! Get max time per it
    IF (npart_local .EQ. 0) THEN 
        local_time_part=0
    ENDIF  
    CALL MPI_ALLREDUCE(local_time_part, global_time_part, 1_isp, MPI_REAL8, MPI_SUM, comm, errcode)
    CALL MPI_ALLREDUCE(npart_local, npart_global, 1_isp, MPI_REAL8, MPI_SUM, comm, errcode)
    IF (npart_global .EQ. 0_idp) THEN 
        global_time_per_part=0_num
    ELSE 
        global_time_per_part=global_time_part/npart_global 
    ENDIF 

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


SUBROUTINE compute_new_split(tppart,tpcell,nx_glob,ny_glob,nz_glob, &
           ncxmin,ncxmax,ncymin,ncymax,nczmin,nczmax, npx,npy,npz)
    IMPLICIT NONE
    REAL(num), INTENT(IN) :: tppart, tpcell
    INTEGER(idp), INTENT(IN) :: nx_glob, ny_glob, nz_glob, npx, npy, npz
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npx) :: ncxmin, ncxmax
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npy) :: ncymin, ncymax
    INTEGER(idp), INTENT(IN OUT), DIMENSION(npz) :: nczmin, nczmax
    REAL(num), DIMENSION(:), ALLOCATABLE :: load_on_x,load_on_y, load_on_z
    ALLOCATE(load_on_x(0:nx_glob-1),load_on_y(0:ny_glob-1),load_on_z(0:nz_glob-1))
    load_on_x=0.
    load_on_y=0.
    load_on_z=0.
    
    ! Compute load in X and compute new split in X 
    CALL get_projected_load_on_x(nx_glob,load_on_x,tppart,tpcell)
    CALL balance_in_dir(load_on_x,nx_glob,npx,ncxmin,ncxmax)
    ! Compute load in X and compute new split in Y 
    CALL get_projected_load_on_y(ny_glob,load_on_y,tppart,tpcell)
    CALL balance_in_dir(load_on_y,ny_glob,npy,ncymin,ncymax)
    ! Compute load in X and compute new split in Z
    CALL get_projected_load_on_z(nz_glob,load_on_z,tppart,tpcell)
    CALL balance_in_dir(load_on_z,nz_glob,npz,nczmin,nczmax)
    
    DEALLOCATE(load_on_x,load_on_y,load_on_z)
    
END SUBROUTINE compute_new_split 
 
 


! This subroutine computes new load in direction dir based on the projected 
! load_in_dir 
SUBROUTINE balance_in_dir(load_in_dir, ncellmaxdir, nproc_in_dir, idirmin, idirmax)
    IMPLICIT NONE 
    REAL(num), DIMENSION(0:ncellmaxdir-1), INTENT(IN) :: load_in_dir
    INTEGER(idp), DIMENSION(0:nproc_in_dir-1), INTENT(IN OUT) :: idirmin, idirmax
    INTEGER(idp), INTENT(IN) :: nproc_in_dir, ncellmaxdir 
    INTEGER(idp) :: iproc, icell 
    REAL(num) :: balanced_load=0_num, curr_balanced_load=0_num, curr_proc_load=0_num  
    LOGICAL(idp) :: not_balanced 
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

! This subroutine computes the computational load projected on Z-Axis 
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

!This subroutine create new array_of_tiles for each species 
SUBROUTINE create_new_tile_split()
IMPLICIT NONE 
TYPE(particle_species), DIMENSION(:), ALLOCATABLE, TARGET :: new_species_parray  
TYPE(particle_species), POINTER :: currsp, currsp_new
TYPE(particle_tile), POINTER :: curr_tile
INTEGER(idp)  :: nxmin,nxmax,nymin,nymax,nzmin,nzmax,&
                            oxmin,oxmax,oymin,oymax,ozmin,ozmax
INTEGER(idp) :: ix, iy, iz, ip, indx, indy, indz, ispecies, count
INTEGER(idp) :: nptile, nx0_grid_tile, ny0_grid_tile, nz0_grid_tile
REAL(num) :: partx, party, partz, partux, partuy, partuz, partw, gaminv
INTEGER(idp) :: ntilex_new, ntiley_new, ntilez_new

 
! Udpate optimal number of tiles 
ntilex_new = MAX(1,nx/10)
ntiley_new = MAX(1,ny/10)
ntilez_new = MAX(1,nz/10)
PRINT *, "rank, ntilexn, ntileyn, ntilezn", rank, ntilex_new , ntiley_new, ntilez_new 

! Allocate new species array 
ALLOCATE(new_species_parray(nspecies))

! Copy properties of each species from old array to new array 
DO ispecies=1,nspecies 
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


CALL set_tile_split_for_species(new_species_parray,nspecies,ntilex_new,ntiley_new,ntilez_new,nx_grid,ny_grid,nz_grid, &
                                x_min_local,y_min_local,z_min_local,x_max_local,y_max_local,z_max_local)


                              
                                
! Deallocate former grid tile array/ ALLOCATE new one 
DEALLOCATE(aofgrid_tiles)
ALLOCATE(aofgrid_tiles(ntilex_new,ntiley_new,ntilez_new))

! Init new tile arrays of new species array 
CALL init_tile_arrays_for_species(nspecies, new_species_parray, aofgrid_tiles, ntilex_new, ntiley_new, ntilez_new)


! Copy particles from former tiles in first tile of new species array
DO ispecies=1,nspecies
    currsp=>species_parray(ispecies)
    currsp_new=>new_species_parray(ispecies)
    ! Get first tiles dimensions (may be different from last tile)
    nx0_grid_tile = currsp_new%array_of_tiles(1,1,1)%nx_grid_tile
    ny0_grid_tile = currsp_new%array_of_tiles(1,1,1)%ny_grid_tile
    nz0_grid_tile = currsp_new%array_of_tiles(1,1,1)%nz_grid_tile
    DO iz=1,ntilez
        DO iy=1,ntiley
            DO ix=1,ntilex 
                curr_tile=>currsp%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                DO ip=1,count
                    partx=curr_tile%part_x(ip)
                    party=curr_tile%part_y(ip)
                    partz=curr_tile%part_z(ip)
                    partux=curr_tile%part_ux(ip)
                    partuy=curr_tile%part_uy(ip)
                    partuz=curr_tile%part_uz(ip)
                    gaminv=curr_tile%part_gaminv(ip)
                    partw=curr_tile%pid(ip,wpid)
                    
                    ! CASE 1: particle outside MPI domain temporarily put it 
                    ! in the first tile of new species array
                    IF ((partx .LT. x_min_local) .OR. (partx .GE. x_max_local) .OR. &
                        (party .LT. y_min_local) .OR. (party .GE. y_max_local) .OR. &
                        (partz .LT. z_min_local) .OR. (partz .GE. z_max_local)) THEN 
                        CALL add_particle_at_tile(currsp_new, 1_idp,1_idp,1_idp, &
                             partx, party, partz, partux, partuy, partuz, gaminv, partw)
                    ! CASE 2: particle is in the new domain just add it to proper tile of new species array
                    ELSE 
                        indx = MIN(FLOOR((partx-x_min_local+dx/2_num)/(nx0_grid_tile*dx))+1,ntilex_new)
                        indy = MIN(FLOOR((party-y_min_local+dy/2_num)/(ny0_grid_tile*dy))+1,ntiley_new)
                        indz = MIN(FLOOR((partz-z_min_local+dz/2_num)/(nz0_grid_tile*dz))+1,ntilez_new)
                        CALL add_particle_at_tile(currsp_new, indx,indy,indz, &
                             partx, party, partz, partux, partuy, partuz, gaminv, partw)
                    
                    ENDIF 
                    currsp_new%species_npart=currsp_new%species_npart+1
                END DO
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
DO ispecies=1,nspecies 
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
CALL set_tile_split_for_species(species_parray,nspecies,ntilex,ntiley,ntilez,nx_grid,ny_grid,nz_grid, &
                                x_min_local,y_min_local,z_min_local,x_max_local,y_max_local,z_max_local)
                                
DEALLOCATE(aofgrid_tiles)
ALLOCATE(aofgrid_tiles(ntilex,ntiley,ntilez))

CALL init_tile_arrays_for_species(nspecies, species_parray, aofgrid_tiles, ntilex, ntiley, ntilez)


! Copy particles 
! Copy particles from former tiles in first tile of new species array
DO ispecies=1,nspecies
    currsp=>species_parray(ispecies)
    currsp_new=>new_species_parray(ispecies)
    DO iz=1,ntilez
        DO iy=1,ntiley
            DO ix=1,ntilex 
                curr_tile=>currsp_new%array_of_tiles(ix,iy,iz)
                count=curr_tile%np_tile(1)
                DO ip=1,count
                    partx=curr_tile%part_x(ip)
                    party=curr_tile%part_y(ip)
                    partz=curr_tile%part_z(ip)
                    partux=curr_tile%part_ux(ip)
                    partuy=curr_tile%part_uy(ip)
                    partuz=curr_tile%part_uz(ip)
                    gaminv=curr_tile%part_gaminv(ip)
                    partw=curr_tile%pid(ip,wpid)
                    CALL add_particle_at_tile(currsp, ix,iy,iz, &
                         partx, party, partz, partux, partuy, partuz, gaminv, partw)
                    
                    currsp%species_npart=currsp%species_npart+1
                END DO
            END DO
        END DO
    END DO 
END DO



!Deallocate new_species_parray 
DEALLOCATE(new_species_parray)

END SUBROUTINE 

SUBROUTINE remap_particles(ix1old,ix2old,iy1old,iy2old,iz1old,iz2old,     &
                           ix1new,ix2new,iy1new,iy2new,iz1new,iz2new,     &
                           ncxmin,ncxmax,ncymin,ncymax,nczmin,nczmax,     &
                           iproc, ncpus,npx, npy,npz, l_cart_comm)
                           
IMPLICIT NONE
INTEGER(idp), INTENT(IN) :: iproc, ncpus, npx, npy, npz
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ix1old, ix2old, iy1old, iy2old, iz1old, iz2old
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ix1new, ix2new, iy1new, iy2new, iz1new, iz2new
INTEGER(idp), INTENT(IN), DIMENSION(0:npx-1) :: ncxmin, ncxmax
INTEGER(idp), INTENT(IN), DIMENSION(0:npy-1) :: ncymin, ncymax
INTEGER(idp), INTENT(IN), DIMENSION(0:npz-1) :: nczmin, nczmax
LOGICAL(idp), INTENT(IN) :: l_cart_comm
INTEGER(idp) :: i, ipart, ixtile, iytile, iztile, nmax, ispec, ispecies 
INTEGER(idp) :: ix3min,ix3max,iy3min,iy3max,iz3min,iz3max
INTEGER(idp) :: ix1newip, ix2newip, iy1newip, iy2newip, iz1newip, iz2newip
INTEGER(idp) :: ix1oldip, ix2oldip, iy1oldip, iy2oldip, iz1oldip, iz2oldip
INTEGER(isp) :: mpitag, count
INTEGER(isp), ALLOCATABLE, DIMENSION(:) :: recv_rank, send_rank, requests
INTEGER(idp), ALLOCATABLE, DIMENSION(:,:) :: npart_recv, npart_send
INTEGER(idp), PARAMETER :: nmax_neighbours=10**3, nvar=8
INTEGER(idp) :: nsend, nrecv, ibuff, curr_rank,iprocx,iprocy,iprocz,icx,icy,icz,isend
LOGICAL(idp) :: l_is_intersection
INTEGER(idp), ALLOCATABLE, DIMENSION (:) :: nptoexch
REAL(num), ALLOCATABLE, DIMENSION (:,:) :: sendbuff
REAL(num), ALLOCATABLE, DIMENSION (:,:) :: recvbuff
REAL(num) :: part_xyz
TYPE(particle_species), POINTER :: currsp
TYPE(particle_tile), POINTER :: curr

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
    CALL get_3Dintersection(ix1newip, ix2newip, iy1newip, iy2newip,         &
                            iz1newip, iz2newip,                             & 
                            ix1old(i), ix2old(i), iy1old(i), iy2old(i),     &
                            iz1old(i), iz2old(i),                           & 
                            ix3min,ix3max,iy3min,iy3max,iz3min,iz3max,      &          
                            l_is_intersection)
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
    CALL get_3Dintersection(ix1oldip, ix2oldip, iy1oldip, iy2oldip,         &
                            iz1oldip, iz2oldip,                             & 
                            ix1new(i), ix2new(i), iy1new(i), iy2new(i),     &
                            iz1new(i), iz2new(i),                           & 
                            ix3min,ix3max,iy3min,iy3max,iz3min,iz3max,      &          
                            l_is_intersection)
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


ALLOCATE(npart_send(nspecies,nsend), npart_recv(nspecies,nrecv), requests(nsend+nrecv))
npart_send=0_idp
npart_recv=0_idp
requests=0_isp

! ----- POST IRECV TO GET NUMBER OF PARTICLES 
DO i=1, nrecv
    count=nspecies
    CALL MPI_IRECV(npart_recv(1,i), count,  MPI_INTEGER8, recv_rank(i), mpitag,    &
                            comm, requests(i), errcode)    
END DO 

! ----- IDENTIFY PARTICLES TO BE SENT/ PLACE PARTICLES IN SEND BUFFER 
ALLOCATE(sendbuff(nvar*npart_local,nsend), nptoexch(nsend))
sendbuff=0_num
nptoexch=0
DO ispecies=1, nspecies !LOOP ON SPECIES
    ! Init send recv buffers
    currsp => species_parray(ispecies)
    DO iztile=1, ntilez !LOOP ON TILES
        DO iytile=1, ntiley
            DO ixtile=1, ntilex
                curr=>currsp%array_of_tiles(ixtile,iytile,iztile)
                ! If not subdomain border, nothing to do
                IF (.NOT. curr%subdomain_bound) CYCLE
                ! search for outbound particles
                part_xyz=0.
                ! Identify outbounds particles and compute destination 
                DO i = 1, curr%np_tile(1) !LOOP ON PARTICLES
                    iprocx = x_coords
                    iprocy = y_coords
                    iprocz = z_coords
                    part_xyz = curr%part_x(i)
                    ! Particle has left this processor
                    IF ((part_xyz .LT. x_min_local) .OR. (part_xyz .GE. x_max_local)) THEN
                        icx=part_xyz/dx
                        CALL get_proc_interval(iprocx,icx,ncxmin,ncxmax,ncpus)
                    ENDIF
                    part_xyz = curr%part_y(i)
                    ! Particle has left this processor
                    IF ((part_xyz .LT. y_min_local) .OR. (part_xyz .GE. y_max_local)) THEN
                        icy=part_xyz/dy
                        CALL get_proc_interval(iprocy,icy,ncymin,ncymax,ncpus)
                    ENDIF       
                    part_xyz = curr%part_z(i)
                    ! Particle has left this processor
                    IF ((part_xyz .LT. z_min_local) .OR. (part_xyz .GE. z_max_local)) THEN 
                        icz=part_xyz/dz
                        ! Find new proc using bissection algorithm (log(nproc))
                        CALL get_proc_interval(iprocz,icz,nczmin,nczmax,ncpus)
                    ENDIF 
                    ! Particles has to be sent to another proc 
                    IF((iprocx .NE. x_coords) .OR. (iprocy .NE. y_coords) .OR. (iprocz .NE. z_coords))  THEN 
                        ! Finds indices in buffer for curr_rank using a binary search algorithm
                        CALL pxr_convertindtoproc(comm,iprocx,iprocy,iprocz,npx,npy,npz,curr_rank,l_cart_comm)
                        CALL binary_search(isend,curr_rank,send_rank(1:nsend),nsend)
                        PRINT *, "curr_rank", rank, curr_rank, isend, nsend, send_rank(1:nsend)
                        nptoexch(isend)=nptoexch(isend)+1
                        ibuff=nvar*nptoexch(isend)+1
                        sendbuff(ibuff,   isend)    = curr%part_x(i)
                        sendbuff(ibuff+1, isend)  = curr%part_y(i)
                        sendbuff(ibuff+2, isend)  = curr%part_z(i)
                        sendbuff(ibuff+3, isend)  = curr%part_ux(i)
                        sendbuff(ibuff+4, isend)  = curr%part_uy(i)
                        sendbuff(ibuff+5, isend)  = curr%part_uz(i)
                        sendbuff(ibuff+6, isend)  = curr%part_gaminv(i)
                        sendbuff(ibuff+7, isend)  = curr%pid(i,wpid)
                        npart_send(i, isend)=npart_send(isend,i)+1
                        ! Remove particle from current tile 
                        CALL rm_particle_at_tile(curr,i)
                    ENDIF           
                ENDDO !END LOOP ON PARTICLES
              ENDDO
           ENDDO
        ENDDO ! END LOOP ON TILES
END DO ! End loop on species

! ----- POST ISEND FOR THE NUMBER OF PARTICLES 
DO i=1, nsend
    count=nspecies
    CALL MPI_ISEND(npart_send(1,i), count,  MPI_INTEGER8, send_rank(i), mpitag,    &
                            comm, requests(nrecv+i), errcode)    
END DO 


! ----- SYNC THE NUMBER OF PARTICLES BEFORE RECEIVING DATA
count=nsend+nrecv
CALL MPI_WAITALL(count,requests, MPI_STATUSES_IGNORE, errcode)
requests=0_isp

! ----- POST IRECV FOR PARTICLE DATA 
nmax=MAXVAL(SUM(npart_recv,1))
ALLOCATE(recvbuff(nsend,nmax))
DO i=1, nrecv
    count=SUM(npart_recv(:,i))
    IF (count .GT. 0) THEN 
        CALL MPI_IRECV(recvbuff(1,i),count, MPI_REAL8,recv_rank(i),mpitag, &
                                comm, requests(i),errcode)
    ENDIF
END DO

! ----- POST ISEND FOR PARTICLES DATA 
DO i=1, nsend
    count=SUM(npart_send(:,i))
    IF (count .GT. 0) THEN 
        CALL MPI_ISEND(sendbuff(1,i), count,  MPI_REAL8, send_rank(i), mpitag,    &
                            comm, requests(nrecv+i), errcode)    
    ENDIF
END DO 

! ----- SYNC MPI EXCHANGES FOR PARTICLE DATA 
count=nsend+nrecv
CALL MPI_WAITALL(count,requests, MPI_STATUSES_IGNORE, errcode)

! ----- ADD PARTICLES IN RECV BUFF IN SPECIES ARRAY 
DO i =1, nrecv
    ispec=1
    DO ispecies=1,nspecies
        currsp=> species_parray(ispecies) 
        DO ipart=1,npart_recv(ispecies,i)
            ibuff=ispec+ipart
            CALL add_particle_to_species(currsp, recvbuff(ibuff,i), recvbuff(ibuff+1,i), recvbuff(ibuff+2,i), &
            recvbuff(ibuff+3,i), recvbuff(ibuff+4,i), recvbuff(ibuff+5,i), recvbuff(ibuff+6,i),recvbuff(ibuff+7,i))
        END DO 
        ispec=ispec+npart_recv(ispecies,i)
    END DO 
END DO

DEALLOCATE(sendbuff,recvbuff,nptoexch,npart_send,npart_recv,requests)

END SUBROUTINE remap_particles

! Finds proc index using binary search algorithm 
SUBROUTINE get_proc_interval(iproc,ic,ncmin,ncmax,ncpus)
IMPLICIT NONE 
INTEGER(idp), INTENT(IN OUT) :: iproc 
INTEGER(idp), INTENT(IN) :: ic, ncpus
INTEGER(idp), INTENT(IN), DIMENSION(0:ncpus-1) :: ncmin, ncmax 
INTEGER(idp) :: imin, imax, imid
LOGICAL(idp) :: is_not_found=.TRUE. 
imin=0
imax=ncpus-1
iproc=-1

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

END SUBROUTINE 

SUBROUTINE binary_search(isend,crank,arr,narr)
IMPLICIT NONE 
INTEGER(idp), INTENT(IN OUT) :: isend
INTEGER(idp), INTENT(IN) :: crank 
INTEGER(idp), INTENT(IN) :: narr
INTEGER(isp), INTENT(IN), DIMENSION(narr) :: arr
INTEGER(idp) :: imin, imax, imid
LOGICAL(idp) :: is_not_found=.TRUE. 

imin=1
imax=narr
isend=-1
PRINT *, "isend, imin, imax, arr, narr",isend, imin, imax, arr, narr
DO WHILE((is_not_found) .AND. (imax .GE. imin) ) 
    imid=(imax+imin)/2
    IF(arr(imid) .EQ. crank) THEN 
        is_not_found=.FALSE. 
        isend=imid
    ELSE
        IF(arr(imid) .LT. crank) THEN 
            imax=imid-1
        ELSE 
            imin=imid+1
        ENDIF
    ENDIF
END DO 

END SUBROUTINE 


END MODULE load_balance 
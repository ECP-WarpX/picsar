MODULE load_balance 
USE fields 
USE particles 
USE params
USE shared_data
IMPLICIT NONE 

CONTAINS 

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


SUBROUTINE balance_in_dir(load_in_dir, ncellmaxdir, nproc_in_dir, idirmin, idirmax)
    IMPLICIT NONE 
    REAL(num), DIMENSION(ncellmaxdir), INTENT(IN) :: load_in_dir
    REAL(num), DIMENSION(nproc_in_dir), INTENT(IN OUT) :: idirmin, idirmax
    INTEGER(idp), INTENT(IN) :: nproc_in_dir, ncellmaxdir 
    INTEGER(idp) :: iproc, icell 
    REAL(num) :: balanced_load=0_num, curr_proc_load=0_num  
    LOGICAL(idp) :: not_balanced 

    balanced_load=SUM(load_in_dir)/nproc_in_dir
    icell=1 

    DO iproc=1, nproc_in_dir
        not_balanced=.TRUE. 
        curr_proc_load=0_num
        idirmin(iproc)=icell
        DO WHILE(not_balanced)  
            curr_proc_load=curr_proc_load+load_in_dir(icell)
            icell=icell+1
            IF(curr_proc_load .GE. balanced_load) not_balanced=.FALSE. 
        END DO 
        idirmax(iproc)=icell-1
    END DO 
END SUBROUTINE balance_in_dir



SUBROUTINE get_projected_load_on_x(nxg,load_on_x,time_per_part,time_per_cell)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: nxg
    REAL(num), INTENT(IN OUT), DIMENSION(nxg) :: load_on_x 
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell 
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icellx
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile),  POINTER :: curr_tile

    ALLOCATE(load_part(nxg))
    ALLOCATE(load_part_sum(nxg))
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
    
END SUBROUTINE get_projected_load_on_x

SUBROUTINE get_projected_load_on_y(nyg,load_on_y,time_per_part,time_per_cell)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: nyg
    REAL(num), INTENT(IN OUT), DIMENSION(nyg) :: load_on_y 
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell 
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part 
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icelly
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile),  POINTER :: curr_tile

    ALLOCATE(load_part(nyg))
    ALLOCATE(load_part_sum(nyg))
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
    
END SUBROUTINE get_projected_load_on_y

SUBROUTINE get_projected_load_on_z(nzg,load_on_z,time_per_part,time_per_cell)
    IMPLICIT NONE 
    INTEGER(idp), INTENT(IN) :: nzg
    REAL(num), INTENT(IN OUT), DIMENSION(nzg) :: load_on_z 
    REAL(num), INTENT(IN) :: time_per_part, time_per_cell 
    INTEGER(idp), DIMENSION(:), ALLOCATABLE :: load_part_sum, load_part 
    INTEGER(idp) :: ispecies, ix, iy, iz, ip, count, icellz
    TYPE(particle_species), POINTER :: curr
    TYPE(particle_tile),  POINTER :: curr_tile

    ALLOCATE(load_part(nzg))
    ALLOCATE(load_part_sum(nzg))
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
    
END SUBROUTINE get_projected_load_on_z


END MODULE load_balance 
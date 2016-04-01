MODULE simple_io

  USE mpi_derived_types
  USE fields
  USE shared_data
  IMPLICIT NONE

CONTAINS

    SUBROUTINE output_routines()
        USE shared_data
        USE params
        USE time_stat
        IMPLICIT NONE
        CHARACTER(LEN=string_length) :: strtemp
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0
        INTEGER(KIND=4) :: err=0
        REAL(num) :: tmptime
        
        IF (output_frequency .LT. 1) RETURN
        tmptime = MPI_WTIME()
        IF ((it .GE. output_step_min) .AND. (it .LE. output_step_max) .AND. &
            (MOD(it-output_step_min,output_frequency) .EQ. 0)) THEN
            !!! --- Write output to disk
            WRITE(strtemp,'(I5)') it
            !! -- Write grid quantities
            IF (c_output_ex .EQ. 1) THEN
                ! - Write current density ex
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileex))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', ex, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_ey .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileey))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', ey, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_ez .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileez))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', ez, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_bx .EQ. 1) THEN
                ! - Write magnetic field bx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filebx))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', bx, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_by .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileby))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', by, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_bz .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filebz))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', bz, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_jx .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filejx))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', jx, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_jy .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filejy))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', jy, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_jz .EQ. 1) THEN
                ! - Write current density jx
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filejz))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', jz, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_dive .EQ. 1) THEN
                ! - Write electric field divergence div E
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filedive))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', dive, nxguards, nyguards, nzguards, nx,ny,nz, offset, err)
            ENDIF
            IF (c_output_rho .EQ. 1) THEN
                ! - Write total charge density rho
                CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filerho))// &
                TRIM(ADJUSTL(strtemp))//'.pxr', rho, nxjguards, nyjguards, nzjguards, nx,ny,nz, offset, err)
            ENDIF
        ENDIF
        
        localtimes(9) = localtimes(9) + (MPI_WTIME() - tmptime)
        
    END SUBROUTINE output_routines


    ! --------------------------------------------------------------------------
    ! Output temporal diagnostics
    ! --------------------------------------------------------------------------
    SUBROUTINE output_temporal_diagnostics
        USE shared_data
        USE params
        USE time_stat
        USE diagnostics
        USE output_data
        USE particle_properties
        USE constants
        USE fields
        IMPLICIT NONE
        
        REAL(num), dimension(:), allocatable :: local_values,global_values
        INTEGER(idp) :: ispecies
        INTEGER(isp) :: i
        REAL(num) :: tmptime
        
        IF ((temdiag_frequency.gt.0).and.(MOD(it,temdiag_frequency).EQ. 0)) THEN
        
          tmptime = MPI_WTIME()
        
          Allocate(local_values(temdiag_totvalues),global_values(temdiag_totvalues))
          
          ! Kinetic energy
          if (temdiag_act_list(1).gt.0) then
            DO ispecies=1,nspecies
              CALL get_loc_kinetic_energy(ispecies,local_values(temdiag_i_list(1)+ispecies-1))
              local_values(temdiag_i_list(1)+ispecies-1) = local_values(temdiag_i_list(1)+ispecies-1)*emass*clight**2
              !if (local_values(temdiag_i_list(1)+ispecies-1) .lt. 1E-15) local_values(temdiag_i_list(1)+ispecies-1)=1e-15
            ENDDO
          end if

          ! Ex energy
          if (temdiag_act_list(2).gt.0) then        
            CALL get_loc_field_energy(ex,nx,ny,nz,dx,dy,dz,nxguards,nyguards,nzguards,local_values(temdiag_i_list(2)))
            local_values(temdiag_i_list(2)) = local_values(temdiag_i_list(2))*eps0
          end if
          
          ! Ey energy
          if (temdiag_act_list(3).gt.0) then        
            CALL get_loc_field_energy(ey,nx,ny,nz,dx,dy,dz,nxguards,nyguards,nzguards,local_values(temdiag_i_list(3)))
            local_values(temdiag_i_list(3)) = local_values(temdiag_i_list(3))*eps0
          end if

          ! Ez energy
          if (temdiag_act_list(4).gt.0) then        
            CALL get_loc_field_energy(ez,nx,ny,nz,dx,dy,dz,nxguards,nyguards,nzguards,local_values(temdiag_i_list(4)))
            local_values(temdiag_i_list(4)) = local_values(temdiag_i_list(4))*eps0
          end if

          ! Bx energy
          if (temdiag_act_list(5).gt.0) then        
            CALL get_loc_field_energy(bx,nx,ny,nz,dx,dy,dz,nxguards,nyguards,nzguards,local_values(temdiag_i_list(5)))
            local_values(temdiag_i_list(5)) = local_values(temdiag_i_list(5))*imu0
          end if
          
          ! By energy
          if (temdiag_act_list(6).gt.0) then        
            CALL get_loc_field_energy(by,nx,ny,nz,dx,dy,dz,nxguards,nyguards,nzguards,local_values(temdiag_i_list(6)))
            local_values(temdiag_i_list(6)) = local_values(temdiag_i_list(6))*imu0
          end if

          ! Bz energy
          if (temdiag_act_list(7).gt.0) then        
            CALL get_loc_field_energy(bz,nx,ny,nz,dx,dy,dz,nxguards,nyguards,nzguards,local_values(temdiag_i_list(7)))
            local_values(temdiag_i_list(7)) = local_values(temdiag_i_list(7))*imu0
          end if

          ! DivE*eps0 - rho
          if (temdiag_act_list(8).gt.0) then  
            CALL get_loc_norm_divErho(dive,rho,nx, ny, nz, nxguards, nyguards, nzguards,local_values(temdiag_i_list(8)))  
            local_values(temdiag_i_list(8)) = local_values(temdiag_i_list(8))
          end if
                    
          ! MPI all reduction
          call MPI_ALLREDUCE(local_values(1),global_values(1),temdiag_totvalues,mpidbl,MPI_SUM,comm,errcode)
          
          ! sqrt for DivE*eps0 - rho
          if (temdiag_act_list(8).gt.0) then 
            global_values(temdiag_i_list(8)) = sqrt(global_values(temdiag_i_list(8)))
          end if
          
          ! _____________
          ! Debug
          !IF (rank==0) THEN
          !  write(0,*) temdiag_totvalues
          !  write(0,*) local_values(temdiag_i_list(1):temdiag_i_list(1)+temdiag_nb_values(1)-1)
          !  write(0,*) global_values(temdiag_i_list(1):temdiag_i_list(1)+temdiag_nb_values(1)-1)
          !  write(0,*) 'exE',global_values(temdiag_i_list(2):temdiag_i_list(2)+temdiag_nb_values(2)-1)
          !  write(0,*) 'eyE',global_values(temdiag_i_list(3):temdiag_i_list(3)+temdiag_nb_values(3)-1)
          !  write(0,*) 'ezE',global_values(temdiag_i_list(4):temdiag_i_list(4)+temdiag_nb_values(4)-1)
          !ENDIF
          
          ! Output
          ! Each mpi task will write in a given file according to their rank
          IF (nproc.ge.temdiag_nb) then
            IF ((rank.ge.0).and.(rank.le.temdiag_nb)) then
              
              ! Ascii format
              IF (temdiag_format.eq.1) then
                !write(42,'(X,3(E14.10E3,X))') local_values(temdiag_i_list(rank+1):temdiag_i_list(rank+1)+temdiag_nb_values(rank+1)-1)
                write(42,*) global_values(temdiag_i_list(rank+1):temdiag_i_list(rank+1)+temdiag_nb_values(rank+1)-1)
                
              ! Binary format
              ELSE
                write(42) global_values(temdiag_i_list(rank+1):temdiag_i_list(rank+1)+temdiag_nb_values(rank+1)-1)
              ENDIF
               
            end if
          else
            if (rank.eq.0) then
              DO i=1,temdiag_nb
                IF (temdiag_format.eq.1) then
                  !write(42,'(X,3(E14.10,X))') local_values(temdiag_i_list(i):temdiag_i_list(i)+temdiag_nb_values(i)-1)  
                  write(42+i,*) global_values(temdiag_i_list(i):temdiag_i_list(i)+temdiag_nb_values(i)-1)  
                ! Binary format    
                else
                  write(42+i) global_values(temdiag_i_list(i):temdiag_i_list(i)+temdiag_nb_values(i)-1)
                end if                          
              ENDDO
            ENDIF       
          ENDIF
          
          !CALL get_kinetic_energy(ispecies,temp_diag(1))
          !if (rank==0) write(0,*) "kinetic energy",temp_diag(1)
          
          localtimes(9) = localtimes(9) + (MPI_WTIME() - tmptime)
          
        ENDIF
        
    END SUBROUTINE

  !----------------------------------------------------------------------------
  ! This subroutine writes a grid quantity (e.g EM fields, Currents)
  ! to disk using MPI-IO (H. VINCENTI)
  !----------------------------------------------------------------------------
  SUBROUTINE write_single_array_to_file(filename, array, nxg, nyg, nzg, nx_local, ny_local, nz_local, offset, err)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(idp), INTENT(IN) :: nxg, nyg, nzg
    INTEGER(idp), INTENT(IN) :: nx_local, ny_local, nz_local
    REAL(num), DIMENSION(-nxg:nx_local+nxg,-nyg:ny_local+nyg,-nzg:nz_local+nzg), INTENT(INOUT) :: array
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER(isp), INTENT(INOUT) :: err
    INTEGER(isp) :: subt, suba, fh, i

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode .NE. 0) THEN
      IF (rank .EQ. 0) PRINT *, 'file ', TRIM(filename), ' could not be created - Check disk space'
      err = IOR(err, c_err_bad_value)
      RETURN
    ENDIF

    subt = create_current_grid_derived_type()
    suba = create_current_grid_subarray(nxguards, nyguards, nzguards)
    CALL MPI_FILE_SET_VIEW(fh, offset, MPI_BYTE, subt, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_WRITE_ALL(fh, array, 1_isp, suba, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_TYPE_FREE(subt, errcode)

  END SUBROUTINE write_single_array_to_file

  SUBROUTINE output_time_statistics
  ! ______________________________________________________
  ! Output of the time statistics
  ! 
  ! ______________________________________________________
    USE time_stat
    USE params
    USE shared_data
    IMPLICIT NONE
    
    REAL(num), DIMENSION(20) :: avetimes
    
    IF ((timestat_period.gt.0).and.(MOD(it,timestat_period).eq.0)) then
    
    
      localtimes(20) = sum(localtimes(1:11))
      localtimes(19) = localtimes(2) + localtimes(4) + localtimes(6) + &
                       localtimes(8) + localtimes(11)
    
      ! Average
      CALL MPI_REDUCE(localtimes,avetimes,20,mpidbl,MPI_SUM,0,comm,errcode)
      avetimes = avetimes / nproc
    
      IF (rank.eq.0) THEN
        write(41) avetimes(1:11)
      end if
      
      
    endif
  
  END SUBROUTINE

END MODULE simple_io

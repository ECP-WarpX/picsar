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

END MODULE simple_io

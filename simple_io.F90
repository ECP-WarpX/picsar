MODULE simple_io

  USE mpi_derived_types
  USE fields

  IMPLICIT NONE

CONTAINS

    SUBROUTINE output_routines()
        USE shared_data
        USE params
        IMPLICIT NONE
        CHARACTER(LEN=string_length) :: strtemp
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0
        INTEGER :: err=0
        !!! --- Write output to disk
        WRITE(strtemp,'(I5)') it
        !! -- Write grid quantities
        ! - Write current density ex
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileex))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', ex, offset, err)
        ! - Write current density jx
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileey))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', ey, offset, err)
        ! - Write current density jx
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileez))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', ez, offset, err)
        ! - Write current density jx
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filebx))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', bx, offset, err)
        ! - Write current density jx
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(fileby))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', by, offset, err)
        ! - Write current density jx
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filebz))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', bz, offset, err)
        ! - Write current density jx
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filejx))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', jx, offset, err)
        ! - Write electric field divergence div E
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filedive))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', dive, offset, err)
        ! - Write total charge density rho
        CALL write_single_array_to_file('./RESULTS/'//TRIM(ADJUSTL(filerho))// &
        TRIM(ADJUSTL(strtemp))//'.pxr', rho, offset, err)
    END SUBROUTINE output_routines






  !----------------------------------------------------------------------------
  ! This subroutine writes a grid quantity (e.g EM fields, Currents)
  ! to disk using MPI-IO (H. VINCENTI)
  !----------------------------------------------------------------------------

  SUBROUTINE write_single_array_to_file(filename, array, offset, err)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(num), DIMENSION(:,:,:), INTENT(INOUT) :: array
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: subt, suba, fh, i

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

    CALL MPI_FILE_WRITE_ALL(fh, array, 1, suba, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_TYPE_FREE(subt, errcode)

  END SUBROUTINE write_single_array_to_file

END MODULE simple_io

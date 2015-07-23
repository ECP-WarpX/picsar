MODULE simple_io

  USE mpi_subtype_control
  USE fields

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! This subroutine writes a grid quantity (e.g EM fields, Currents)
  ! to disk using MPI-IO (H. VINCENTI)
  !----------------------------------------------------------------------------

  SUBROUTINE write_single_array_to_file(filename, array, offset, err)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(num), DIMENSION(:,:,:), INTENT(INOUT) :: array
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: ng
    INTEGER :: subtype, subarray, fh, i

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode .NE. 0) THEN
      IF (rank .EQ. 0) PRINT *, 'file ', TRIM(filename), ' could not be created - Check disk space'
      err = IOR(err, c_err_bad_value)
      RETURN
    ENDIF

    ng=nxguards ! we assume nxguards=nyguards=nzguards
    subtype = create_current_field_subtype()
    subarray = create_current_field_subarray(ng)
    CALL MPI_FILE_SET_VIEW(fh, offset, MPI_BYTE, subtype, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_WRITE_ALL(fh, array, 1, subarray, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_TYPE_FREE(subtype, errcode)

  END SUBROUTINE write_single_array_to_file

END MODULE simple_io

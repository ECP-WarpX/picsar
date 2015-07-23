MODULE mpi_subtype_control

  !----------------------------------------------------------------------------
  ! This module contains the subroutines which create the subarray used
  ! in boundary exchange for the fields
  !----------------------------------------------------------------------------

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! Frees the subtypes created by create_subtypes
  !----------------------------------------------------------------------------

  SUBROUTINE free_subtypes

    CALL MPI_TYPE_FREE(subtype_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big, errcode)

    CALL MPI_TYPE_FREE(subtype_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big_r4, errcode)

  END SUBROUTINE free_subtypes



  !----------------------------------------------------------------------------
  ! create_current_field_subtype - Creates the subtype corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subtype(basetype_in)

    INTEGER :: create_current_field_subtype
    INTEGER, OPTIONAL, INTENT(IN) :: basetype_in
    INTEGER :: basetype

    IF (PRESENT(basetype_in)) THEN
      basetype = basetype_in
    ELSE
      basetype = mpireal
    ENDIF

    create_current_field_subtype = &
        create_field_subtype(basetype, nx, ny, nz, nx_global_min, &
            ny_global_min, nz_global_min)

  END FUNCTION create_current_field_subtype



  !----------------------------------------------------------------------------
  ! create_current_field_subarray - Creates the subarray corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subarray(ng, basetype_in)

    INTEGER :: create_current_field_subarray
    INTEGER, INTENT(IN) :: ng
    INTEGER, OPTIONAL, INTENT(IN) :: basetype_in
    INTEGER :: basetype

    IF (PRESENT(basetype_in)) THEN
      basetype = basetype_in
    ELSE
      basetype = mpireal
    ENDIF

    create_current_field_subarray = &
        create_field_subarray(basetype, ng, nx+1, ny+1, nz+1)

  END FUNCTION create_current_field_subarray


  !----------------------------------------------------------------------------
  ! create_field_subtype - Creates a subtype representing the local processor
  ! for any arbitrary arrangement of an array covering the entire spatial
  ! domain. Only used directly during load balancing
  !----------------------------------------------------------------------------

  FUNCTION create_field_subtype(basetype, nx_local, ny_local, nz_local, &
      cell_start_x_local, cell_start_y_local, cell_start_z_local)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, INTENT(IN) :: nx_local
    INTEGER, INTENT(IN) :: ny_local
    INTEGER, INTENT(IN) :: nz_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER, INTENT(IN) :: cell_start_y_local
    INTEGER, INTENT(IN) :: cell_start_z_local
    INTEGER :: create_field_subtype
    INTEGER, DIMENSION(c_ndims) :: n_local, n_global, start

    n_local = (/nx_local+1, ny_local+1, nz_local+1/)
    n_global = (/nx_global+1, ny_global+1, nz_global+1/)
    start = (/cell_start_x_local, cell_start_y_local, cell_start_z_local/)

    create_field_subtype = &
        create_3d_array_subtype(basetype, n_local, n_global, start)

  END FUNCTION create_field_subtype



  !----------------------------------------------------------------------------
  ! create_1d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 1D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_1d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec1d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(1), INTENT(IN) :: n_local
    INTEGER, DIMENSION(1), INTENT(IN) :: n_global
    INTEGER, DIMENSION(1), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(1)
    INTEGER :: vec1d, vec1d_sub, typesize

    vec1d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CONTIGUOUS(n_local(1), basetype, vec1d, errcode)
    CALL MPI_TYPE_COMMIT(vec1d, errcode)

    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    starts = start - 1
    lengths = 1

    disp(1) = 0
    disp(2) = typesize * starts(1)
    disp(3) = typesize * n_global(1)
    types(1) = MPI_LB
    types(2) = vec1d
    types(3) = MPI_UB

    vec1d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec1d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec1d_sub, errcode)

    CALL MPI_TYPE_FREE(vec1d, errcode)

  END FUNCTION create_1d_array_subtype



  !----------------------------------------------------------------------------
  ! create_2d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 2D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_2d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec2d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(2), INTENT(IN) :: n_local
    INTEGER, DIMENSION(2), INTENT(IN) :: n_global
    INTEGER, DIMENSION(2), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(2)
    INTEGER :: vec2d, vec2d_sub, typesize

    vec2d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_VECTOR(n_local(2), n_local(1), n_global(1), basetype, &
        vec2d, errcode)
    CALL MPI_TYPE_COMMIT(vec2d, errcode)

    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    starts = start - 1
    lengths = 1

    disp(1) = 0
    disp(2) = typesize * (starts(1) + n_global(1) * starts(2))
    disp(3) = typesize * n_global(1) * n_global(2)
    types(1) = MPI_LB
    types(2) = vec2d
    types(3) = MPI_UB

    vec2d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec2d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec2d_sub, errcode)

    CALL MPI_TYPE_FREE(vec2d, errcode)

  END FUNCTION create_2d_array_subtype



  !----------------------------------------------------------------------------
  ! create_3d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 3D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_3d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec3d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(3), INTENT(IN) :: n_local
    INTEGER, DIMENSION(3), INTENT(IN) :: n_global
    INTEGER, DIMENSION(3), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(3)
    INTEGER :: vec2d, vec2d_sub
    INTEGER :: vec3d, vec3d_sub, typesize

    vec2d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_VECTOR(n_local(2), n_local(1), n_global(1), basetype, &
        vec2d, errcode)
    CALL MPI_TYPE_COMMIT(vec2d, errcode)

    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    starts = start - 1
    lengths = 1

    disp(1) = 0
    disp(2) = typesize * (starts(1) + n_global(1) * starts(2))
    disp(3) = typesize * n_global(1) * n_global(2)
    types(1) = MPI_LB
    types(2) = vec2d
    types(3) = MPI_UB

    vec2d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec2d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec2d_sub, errcode)

    vec3d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CONTIGUOUS(n_local(3), vec2d_sub, vec3d, errcode)
    CALL MPI_TYPE_COMMIT(vec3d, errcode)

    disp(1) = 0
    disp(2) = typesize * n_global(1) * n_global(2) * starts(3)
    disp(3) = typesize * n_global(1) * n_global(2) * n_global(3)
    types(1) = MPI_LB
    types(2) = vec3d
    types(3) = MPI_UB

    vec3d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec3d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec3d_sub, errcode)

    CALL MPI_TYPE_FREE(vec2d, errcode)
    CALL MPI_TYPE_FREE(vec2d_sub, errcode)
    CALL MPI_TYPE_FREE(vec3d, errcode)

  END FUNCTION create_3d_array_subtype



  FUNCTION create_field_subarray(basetype, ng, n1, n2, n3)

    INTEGER, INTENT(IN) :: basetype, ng, n1
    INTEGER, INTENT(IN), OPTIONAL :: n2, n3
    INTEGER, DIMENSION(3) :: n_local, n_global, start
    INTEGER :: i, ndim, create_field_subarray

    n_local(1) = n1
    ndim = 1
    IF (PRESENT(n2)) THEN
      n_local(2) = n2
      ndim = 2
    ENDIF
    IF (PRESENT(n3)) THEN
      n_local(3) = n3
      ndim = 3
    ENDIF

    DO i = 1, ndim
      start(i) = 1 + ng
      n_global(i) = n_local(i) + 2 * ng
    ENDDO

    IF (PRESENT(n3)) THEN
      create_field_subarray = &
          create_3d_array_subtype(basetype, n_local, n_global, start)
    ELSE IF (PRESENT(n2)) THEN
      create_field_subarray = &
          create_2d_array_subtype(basetype, n_local, n_global, start)
    ELSE
      create_field_subarray = &
          create_1d_array_subtype(basetype, n_local, n_global, start)
    ENDIF

  END FUNCTION create_field_subarray

END MODULE mpi_subtype_control

MODULE mpi_derived_types

  !----------------------------------------------------------------------------
  ! This module contains the subroutines which create the subarray used
  ! in boundary exchange for the fields
  !----------------------------------------------------------------------------

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! create_current_field_derived_type - Creates the derived type 
  ! corresponding to the current CPU split
  !----------------------------------------------------------------------------

  FUNCTION create_current_grid_derived_type()

    INTEGER :: create_current_grid_derived_type

    create_current_grid_derived_type = &
        create_grid_derived_type(mpidbl, nx, ny, nz, nx_global_grid_min, &
            ny_global_grid_min, nz_global_grid_min)

  END FUNCTION create_current_grid_derived_type

!----------------------------------------------------------------------------
! create_current_field_derived_type - Creates the derived type
! corresponding to the current CPU split
!----------------------------------------------------------------------------

  FUNCTION create_current_grid_subarray(ngx,ngy,ngz)

    INTEGER :: create_current_grid_subarray
    INTEGER, INTENT(IN) :: ngx, ngy, ngz


    create_current_grid_subarray = &
    create_grid_subarray(mpidbl, ngx, ngy, ngz, &
                              nx+1, ny+1, nz+1)

  END FUNCTION create_current_grid_subarray


  !----------------------------------------------------------------------------
  ! create_grid_derived_type - Creates a derived type representing the layout 
  ! of local CPU among the global simulation domain
  !----------------------------------------------------------------------------

  FUNCTION create_grid_derived_type(mpitype, nx_local, ny_local, nz_local, &
      cell_start_x_local, cell_start_y_local, cell_start_z_local)

    INTEGER, INTENT(IN) :: mpitype
    INTEGER, INTENT(IN) :: nx_local
    INTEGER, INTENT(IN) :: ny_local
    INTEGER, INTENT(IN) :: nz_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER, INTENT(IN) :: cell_start_y_local
    INTEGER, INTENT(IN) :: cell_start_z_local
    INTEGER :: create_grid_derived_type
    INTEGER, DIMENSION(c_ndims) :: n_local, n_global, start

    n_local = (/nx_local+1, ny_local+1, nz_local+1/)
    n_global = (/nx_global+1, ny_global+1, nz_global+1/)
    start = (/cell_start_x_local, cell_start_y_local, cell_start_z_local/)

    create_grid_derived_type = &
        create_3d_array_derived_type(mpitype, n_local, n_global, start)

  END FUNCTION create_grid_derived_type

  !----------------------------------------------------------------------------
  ! create_3d_array_derived_type - Creates a derived type representing the 
  ! localization of current CPU among simulation domain
  !----------------------------------------------------------------------------

  FUNCTION create_3d_array_derived_type(mpitype, n_local, n_global, start) &
      RESULT(vec3d_sub)

    INTEGER, INTENT(IN) :: mpitype
    INTEGER, DIMENSION(3), INTENT(IN) :: n_local
    INTEGER, DIMENSION(3), INTENT(IN) :: n_global
    INTEGER, DIMENSION(3), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(3)
    INTEGER :: vec2d, vec2d_sub
    INTEGER :: vec3d, vec3d_sub, typesize

    vec2d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_VECTOR(n_local(2), n_local(1), n_global(1), mpitype, &
        vec2d, errcode)
    CALL MPI_TYPE_COMMIT(vec2d, errcode)

    CALL MPI_TYPE_SIZE(mpitype, typesize, errcode)
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

  END FUNCTION create_3d_array_derived_type



  FUNCTION create_grid_subarray(mpitype, ng1, ng2, ng3, n1, n2, n3)

    INTEGER, INTENT(IN) :: mpitype, ng1,ng2,ng3, n1, n2, n3
    INTEGER, DIMENSION(3) :: n_local, ng, n_global, start
    INTEGER :: i, ndim, create_grid_subarray

    n_local(1) = n1
    n_local(2) = n2
    n_local(3) = n3
    ng(1)= ng1
    ng(2)= ng2
    ng(3)= ng3
    ndim = 3
    DO i = 1, ndim
      start(i) = 1 + ng(i)
      n_global(i) = n_local(i) + 2 * ng(i)
    ENDDO

    create_grid_subarray = &
          create_3d_array_derived_type(mpitype, n_local, n_global, start)

  END FUNCTION create_grid_subarray

END MODULE mpi_derived_types

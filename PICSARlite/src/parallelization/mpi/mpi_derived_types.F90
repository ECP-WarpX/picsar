! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! “Particle In Cell Scalable Application Resource (PICSAR) v2”, Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! MPI_DERIVED_TYPES.F90
!
! Purpose:
! This file contains subroutines to create diverse MPI derived types.
!
! Authors:
! Henri Vincenti
!
! Date:
! Creation 2015
! ________________________________________________________________________________________

! ________________________________________________________________________________________
!> @brief
!> This module contains the subroutines which create the subarray used
!> in boundary exchange for the fields.
!
!> @author
!> Henri Vincenti
!
!> @date
!> Creation 2015
! ________________________________________________________________________________________
MODULE mpi_derived_types!#do not parse
  USE shared_data
  IMPLICIT NONE

  CONTAINS
  ! ______________________________________________________________________________________
  !> @brief
  !> create_current_field_derived_type - Creates the derived type
  !> corresponding to the current CPU split.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  FUNCTION create_current_grid_derived_type()
    INTEGER(isp) :: create_current_grid_derived_type

    create_current_grid_derived_type = create_grid_derived_type(mpidbl, nx, ny, nz,   &
    nx_global_grid_min, ny_global_grid_min, nz_global_grid_min)

  END FUNCTION create_current_grid_derived_type

  ! ______________________________________________________________________________________
  !> @brief
  !> Create a subarray from the current grid according to the given parameters.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  FUNCTION create_current_grid_subarray(ngx, ngy, ngz)

    INTEGER(isp) :: create_current_grid_subarray
    INTEGER(idp), INTENT(IN) :: ngx, ngy, ngz
    INTEGER(idp) :: nxloc, nyloc, nzloc
    nxloc=nx+1
    nyloc=ny+1
    nzloc=nz+1

    create_current_grid_subarray = create_grid_subarray(mpidbl, ngx, ngy, ngz, nxloc, &
    nyloc, nzloc)
  END FUNCTION create_current_grid_subarray

  ! ______________________________________________________________________________________
  !> @brief
  !> create_grid_derived_type - Creates a derived type representing the layout
  !> of local CPU among the global simulation domain
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  FUNCTION create_grid_derived_type(mpitype, nx_local, ny_local, nz_local,            &
    cell_start_x_local, cell_start_y_local, cell_start_z_local)
    INTEGER(isp), INTENT(IN) :: mpitype
    INTEGER(idp), INTENT(IN) :: nx_local
    INTEGER(idp), INTENT(IN) :: ny_local
    INTEGER(idp), INTENT(IN) :: nz_local
    INTEGER(idp), INTENT(IN) :: cell_start_x_local
    INTEGER(idp), INTENT(IN) :: cell_start_y_local
    INTEGER(idp), INTENT(IN) :: cell_start_z_local
    INTEGER(isp) :: create_grid_derived_type, ndims
    INTEGER(isp), DIMENSION(c_ndims) :: n_local, n_global, start

    n_local = (/nx_local, ny_local, nz_local/)
    n_global = (/nx_global, ny_global, nz_global/)
    start = (/cell_start_x_local, cell_start_y_local, cell_start_z_local/)
    ! New version with MPI_TYPE_CREATE_SUBARRAY
    ndims=c_ndims
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, n_global, n_local, start, MPI_ORDER_FORTRAN, &
    mpitype, create_grid_derived_type, errcode)
    CALL MPI_TYPE_COMMIT(create_grid_derived_type, errcode)

  END FUNCTION create_grid_derived_type

  ! ______________________________________________________________________________________
  !> @brief
  !> create_3d_array_derived_type OLD VERSION - USE MPI_TYPE_CREATE_SUBARRAY
  !> instead
  !> - Creates a derived type representing the
  !> localization of current CPU among simulation domain
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  ! ______________________________________________________________________________________
  FUNCTION create_3d_array_derived_type(mpitype, n_local, n_global, start)            &
    RESULT(vec3d_sub)
    INTEGER(isp), INTENT(IN) :: mpitype
    INTEGER(idp), DIMENSION(c_ndims), INTENT(IN) :: n_local
    INTEGER(idp), DIMENSION(c_ndims), INTENT(IN) :: n_global
    INTEGER(idp), DIMENSION(c_ndims), INTENT(IN) :: start
    INTEGER(isp), DIMENSION(c_ndims) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(3)
    INTEGER(isp) :: vec2d, vec2d_sub
    INTEGER(isp) :: vec3d, vec3d_sub, typesize

    vec2d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_VECTOR(INT(n_local(2), isp), INT(n_local(1), isp), INT(n_global(1), &
    isp), mpitype, vec2d, errcode)
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
    CALL MPI_TYPE_CREATE_STRUCT(c_ndims, lengths, disp, types, vec2d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec2d_sub, errcode)

    vec3d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CONTIGUOUS(INT(n_local(3), isp), vec2d_sub, vec3d, errcode)
    CALL MPI_TYPE_COMMIT(vec3d, errcode)

    disp(1) = 0
    disp(2) = typesize * n_global(1) * n_global(2) * starts(3)
    disp(3) = typesize * n_global(1) * n_global(2) * n_global(3)
    types(1) = MPI_LB
    types(2) = vec3d
    types(3) = MPI_UB

    vec3d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(c_ndims, lengths, disp, types, vec3d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec3d_sub, errcode)

    CALL MPI_TYPE_FREE(vec2d, errcode)
    CALL MPI_TYPE_FREE(vec2d_sub, errcode)
    CALL MPI_TYPE_FREE(vec3d, errcode)

  END FUNCTION create_3d_array_derived_type


  ! ______________________________________________________________________________________
  !> @brief
  !> This subroutine creates a subarray from a given grid according
  !> to the given parameters.
  !
  !> @author
  !> Henri Vincenti
  !
  !> @date
  !> Creation 2015
  !
  !> @param[in] mpitype MPI type
  !> @param[in] ng1 guard cells in x
  !> @param[in] ng2 guard cells in y
  !> @param[in] ng3 guard cells in z
  !> @param[in] n1 number of cells in x
  !> @param[in] n2 number of cells in y
  !> @param[in] n3 number of cells in z
  !
  ! ______________________________________________________________________________________
  FUNCTION create_grid_subarray(mpitype, ng1, ng2, ng3, n1, n2, n3)
    INTEGER(isp), INTENT(IN) :: mpitype
    INTEGER(idp), INTENT(IN) :: ng1, ng2, ng3, n1, n2, n3
    INTEGER(idp), DIMENSION(3) :: n_local, ng, n_global, start
    INTEGER(isp) :: i, ndim, create_grid_subarray

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
    n_local=n_local-1! remove last point

    ! New version with  MPI_TYPE_CREATE_SUBARRAY
    CALL MPI_TYPE_CREATE_SUBARRAY(ndim, INT(n_global, isp), INT(n_local, isp),        &
    INT(start-1, isp), MPI_ORDER_FORTRAN, mpitype, create_grid_subarray, errcode)
    CALL MPI_TYPE_COMMIT(create_grid_subarray, errcode)

  END FUNCTION create_grid_subarray

END MODULE mpi_derived_types

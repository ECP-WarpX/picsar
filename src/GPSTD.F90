MODULE matrix_data
USE constants
! Maximum number of instances (matrix_blocks and vector_blocks)
INTEGER(idp), PARAMETER ::  ns_max=40
INTEGER(idp) :: nmatrixes=0
END MODULE matrix_data


! MODULE DEFINING BLOCKS AND MATRIX BLOCKS
! This can be used for interfacing GPSTD with PYTHON 
! Or directly implementing GPSTD in FORTRAN 90  
MODULE matrix_coefficients !#do not parse 
USE matrix_data
IMPLICIT NONE 
TYPE block3d
    COMPLEX(cpx), POINTER, DIMENSION(:,:,:) :: block3dc
    INTEGER(idp) :: nx, ny, nz
    LOGICAL(idp) :: is_source_variable = .FALSE. 
END TYPE block3d
TYPE matrix_blocks
    TYPE(block3d), POINTER, DIMENSION(:,:) :: block_matrix2d
    INTEGER(idp) :: nblocks
END TYPE matrix_blocks
TYPE vector_blocks
    TYPE(block3d), POINTER, DIMENSION(:) :: block_vector
    INTEGER(idp) :: nblocks
END TYPE vector_blocks

! Array of 2D block matrixes (contaning 3d blocks coefficients for GPSTD_Maxwell, GPSTD_Maxwell_PML etc.)
TYPE(matrix_blocks), POINTER, DIMENSION(:) :: cc_mat
! Arrays of 1D block vectors  (containing 3d blocks Ex, Ey, Ez etc.)
TYPE(vector_blocks), POINTER, DIMENSION(:) :: vnew, vold
END MODULE matrix_coefficients 


! SUBROUTINE for associating fortran block pointers to python numpy arrays
! (e.g mymat[][])
SUBROUTINE point_to_matrix_block_p2f(ain,n1,n2,n3,bid1,bid2,mat_index)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(8), INTENT(IN) :: n1,n2,n3, mat_index, bid1, bid2
COMPLEX(cpx), INTENT(IN), TARGET, DIMENSION(n1,n2,n3) :: ain

cc_mat(mat_index)%block_matrix2d(bid1,bid2)%block3dc=>ain
cc_mat(mat_index)%block_matrix2d(bid1,bid2)%nx = n1
cc_mat(mat_index)%block_matrix2d(bid1,bid2)%ny = n2 
cc_mat(mat_index)%block_matrix2d(bid1,bid2)%nz = n3

END SUBROUTINE point_to_matrix_block_p2f

SUBROUTINE point_to_matrix_block(ain,n1,n2,n3,mat_index, bid1, bid2)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(8), INTENT(IN) :: n1,n2,n3, mat_index, bid1, bid2
COMPLEX(cpx), INTENT(IN OUT), DIMENSION(n1,n2,n3) :: ain

ain = cc_mat(mat_index)%block_matrix2d(bid1,bid2)%block3dc(:,:,:)
END SUBROUTINE point_to_matrix_block

SUBROUTINE point_to_vec_block(ain,n1,n2,n3,mat_index, bid1, old)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(8), INTENT(IN) :: n1,n2,n3, mat_index, bid1
LOGICAL(8), INTENT(IN) :: old
COMPLEX(cpx), INTENT(IN OUT), DIMENSION(n1,n2,n3) :: ain
IF (old) THEN 
    ain = vold(mat_index)%block_vector(bid1)%block3dc(:,:,:)
ELSE
    ain = vnew(mat_index)%block_vector(bid1)%block3dc(:,:,:)
ENDIF
END SUBROUTINE point_to_vec_block

SUBROUTINE modify_vec_block(value,mat_index, bid1, old)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(8), INTENT(IN) :: mat_index, bid1
LOGICAL(8), INTENT(IN) :: old
COMPLEX(cpx), INTENT(IN) :: value
IF (old) THEN 
    vold(mat_index)%block_vector(bid1)%block3dc(:,:,:) = value
ELSE
    vnew(mat_index)%block_vector(bid1)%block3dc(:,:,:) = value
ENDIF
END SUBROUTINE modify_vec_block


! SUBROUTINE for associating fortran block pointers to python numpy arrays 
! (Ex, Ey etc.)
SUBROUTINE point_to_vector_block_p2f(ain,n1,n2,n3,iv,mat_index,old,is_source)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: n1,n2,n3, mat_index, iv
LOGICAL(idp), INTENT(IN) :: old, is_source
COMPLEX(cpx), INTENT(IN), TARGET, DIMENSION(n1,n2,n3) :: ain

IF (old) THEN 
    vold(mat_index)%block_vector(iv)%block3dc=>ain
    vold(mat_index)%block_vector(iv)%nx = n1
    vold(mat_index)%block_vector(iv)%ny = n2
    vold(mat_index)%block_vector(iv)%nz = n3
    vold(mat_index)%block_vector(iv)%is_source_variable = is_source
ELSE
    vnew(mat_index)%block_vector(iv)%block3dc=>ain
    vnew(mat_index)%block_vector(iv)%nx = n1
    vnew(mat_index)%block_vector(iv)%ny = n2
    vnew(mat_index)%block_vector(iv)%nz = n3
    vnew(mat_index)%block_vector(iv)%is_source_variable = is_source
ENDIF

END SUBROUTINE point_to_vector_block_p2f


SUBROUTINE nullify_vector_block(bid,mat_index)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(8), INTENT(IN) :: mat_index, bid

NULLIFY(vold(mat_index)%block_vector(bid)%block3dc)
NULLIFY(vnew(mat_index)%block_vector(bid)%block3dc)

END SUBROUTINE nullify_vector_block

! SUBROUTINE THAT ALLOCATE NEW BLOCK MATRIXES/VECTORS
SUBROUTINE allocate_new_matrix_vector(nvar)
USE matrix_coefficients 
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: nvar

IF (.NOT. associated(cc_mat)) THEN
    ALLOCATE(cc_mat(ns_max))
ENDIF

IF (.NOT. associated(vold)) THEN
    ALLOCATE(vold(ns_max))
ENDIF

IF (.NOT. associated(vnew)) THEN
    ALLOCATE(vnew(ns_max))
ENDIF

nmatrixes=nmatrixes+1
ALLOCATE(cc_mat(nmatrixes)%block_matrix2d(nvar,nvar), &
         vold(nmatrixes)%block_vector(nvar), &
         vnew(nmatrixes)%block_vector(nvar))
cc_mat(nmatrixes)%nblocks=nvar 
vold(nmatrixes)%nblocks=nvar 
vnew(nmatrixes)%nblocks=nvar 

END SUBROUTINE allocate_new_matrix_vector


SUBROUTINE multiply_mat_vector(matrix_index)
USE matrix_coefficients 
USE omp_lib
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: matrix_index
INTEGER(idp) :: irow, icol, nrow, ncol, nthreads_tot, nthreads_loop1, nthreads_loop2
INTEGER(idp) :: n1vec, n2vec, n3vec, n1mat, n2mat, n3mat
TYPE(block3d), POINTER :: pvec_new, pvec_old, p_mat
nrow=cc_mat(matrix_index)%nblocks
ncol=nrow
 
#ifdef _OPENMP
    nthreads_tot=OMP_GET_MAX_THREADS()
    CALL OMP_SET_NESTED(.TRUE.)
#else
    nthreads_tot=1
#endif

IF (nthreads_tot .GT. 1) THEN 
    nthreads_loop2 = 2
    nthreads_loop1=nthreads_tot/nthreads_loop2
    IF (nthreads_loop1*nthreads_loop2 .NE. nthreads_tot) THEN 
        nthreads_loop1=nthreads_tot
        nthreads_loop2=1
    ENDIF 
ELSE 
    nthreads_loop1=1
    nthreads_loop2=1
ENDIF 

!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(NONE) & 
!$OMP SHARED(nrow,ncol,nthreads_loop2,vnew,vold,cc_mat,matrix_index) &
!$OMP PRIVATE(irow,icol,pvec_new,pvec_old,p_mat,n1vec,n2vec,n3vec,n1mat,n2mat,n3mat) &
!$OMP NUM_THREADS(nthreads_loop1) 
DO irow=1,nrow
    pvec_new=>vnew(matrix_index)%block_vector(irow)
    IF (pvec_new%is_source_variable) THEN 
        CYCLE
    ENDIF
    pvec_new%block3dc=0.
    DO icol=1,ncol
        pvec_old=>vold(matrix_index)%block_vector(icol)
        p_mat=>cc_mat(matrix_index)%block_matrix2d(irow,icol)
        n1vec=pvec_old%nx
        n2vec=pvec_old%ny
        n3vec=pvec_old%nz
        n1mat=p_mat%nx
        n2mat=p_mat%ny
        n3mat=p_mat%nz
        CALL multiply_unit_blocks(pvec_new%block3dc,pvec_old%block3dc, & 
                n1vec,n2vec,n3vec, p_mat%block3dc,n1mat, n2mat, n3mat, nthreads_loop2)
    END DO 
END DO 
!$OMP END PARALLEL DO 

END SUBROUTINE 

SUBROUTINE multiply_unit_blocks(anew,block1,n1,n2,n3,coeff1,nc1,nc2,nc3,nthreads)
USE constants 
USE omp_lib
IMPLICIT NONE 
INTEGER(idp), INTENT(IN) :: n1,n2,n3, nc1,nc2,nc3, nthreads
COMPLEX(cpx), INTENT(IN OUT), DIMENSION(n1,n2,n3) :: anew, block1
COMPLEX(cpx), INTENT(IN OUT), DIMENSION(nc1,nc2,nc3) :: coeff1
INTEGER(idp) :: i,j,k 

IF (nc1*nc2*nc3 .EQ. 1) THEN 
    IF ((REALPART(coeff1(1,1,1)) .EQ. 0.) .AND. (IMAGPART(coeff1(1,1,1)) .EQ. 0.)) THEN 
       RETURN 
    ELSE IF ((REALPART(coeff1(1,1,1)) .EQ. 1.) .AND. (IMAGPART(coeff1(1,1,1)) .EQ. 0.)) THEN 
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i) & 
        !$OMP SHARED(anew,block1, n1,n2,n3) NUM_THREADS(nthreads) 
        DO k=1,n3
            DO j=1,n2
                DO i=1,n1
                    anew(i,j,k)=anew(i,j,k)+block1(i,j,k)
                END DO 
            END DO 
        END DO 
        !$OMP END PARALLEL DO
    ELSE 
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i) & 
        !$OMP SHARED(anew,block1,coeff1, n1,n2,n3) NUM_THREADS(nthreads) 
        DO k=1,n3
            DO j=1,n2
                DO i=1,n1
                    anew(i,j,k)=anew(i,j,k)+block1(i,j,k)*coeff1(1,1,1)
                END DO 
            END DO 
        END DO 
        !$OMP END PARALLEL DO 
    END IF
ELSE 
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i) & 
    !$OMP SHARED(anew,block1,coeff1, n1,n2,n3) NUM_THREADS(nthreads) 
    DO k=1,n3
        DO j=1,n2
            DO i=1,n1
                anew(i,j,k)=anew(i,j,k)+block1(i,j,k)*coeff1(i,j,k)
            END DO 
        END DO 
    END DO 
    !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE multiply_unit_blocks


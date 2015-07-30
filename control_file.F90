MODULE control_file

  USE shared_data
  USE params
  USE fields
  USE particles
  IMPLICIT NONE

  INTEGER :: ios=0
  INTEGER, PARAMETER :: fh_input = 15
  CHARACTER(LEN=string_length) :: buffer
  CHARACTER(LEN=string_length) :: section_name

CONTAINS

    SUBROUTINE read_input_file
        INTEGER :: ix = 0
        ! --- OPENS INPUT FILE
        OPEN(fh_input, file='input_file.pixr')
        DO WHILE(ios==0)
            READ(fh_input, '(A)', iostat=ios) buffer
            ix=INDEX(buffer,'section::')
            IF (ix .GT. 0) THEN
                section_name=buffer(ix:string_length)
                SELECT CASE(TRIM(ADJUSTL(section_name)))
                CASE('section::grid')
                    CALL read_grid_section
                CASE('section::species')
                    CALL read_species_section
                CASE('section::output')
                    CALL read_output_section
                CASE('section::cpusplit')
                    CALL read_cpusplit_section
                END SELECT
            END IF
        END DO
        RETURN
    END SUBROUTINE read_input_file


    SUBROUTINE read_cpusplit_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS CPUSPLIT SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            IF (INDEX(buffer,'nprocx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nprocx
            ELSE IF (INDEX(buffer,'nprocy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nprocy
            ELSE IF (INDEX(buffer,'nprocz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nprocz
            ELSE IF (INDEX(buffer,'end::cpusplit') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_cpusplit_section

    SUBROUTINE read_grid_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS GRID SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            IF (INDEX(buffer,'nx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nx_global
            ELSE IF (INDEX(buffer,'ny') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') ny_global
            ELSE IF (INDEX(buffer,'nz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nz_global
            ELSEIF (INDEX(buffer,'dx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) dx
            ELSE IF (INDEX(buffer,'dy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) dy
            ELSE IF (INDEX(buffer,'dz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) dz
            ELSEIF (INDEX(buffer,'x_grid_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) xmin
            ELSE IF (INDEX(buffer,'y_grid_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) ymin
            ELSE IF (INDEX(buffer,'z_grid_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) zmin
            ELSE IF (INDEX(buffer,'end::grid') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_grid_section

    SUBROUTINE read_species_section
        INTEGER :: ix = 0
        LOGICAL :: end_section
        TYPE(particle_species), POINTER :: curr
        ! READS SPECIES SECTION OF INPUT FILE
        IF (.NOT. l_species_allocated) THEN
            nspecies=0
            ALLOCATE(species_parray(1:nspecies_max))
            l_species_allocated=.TRUE.
        ENDIF
        nspecies = nspecies+1
        curr => species_parray(nspecies)
        end_section=.FALSE.
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            IF (INDEX(buffer,'name') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%name
            ELSE IF (INDEX(buffer,'mass') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%mass
                curr%mass=curr%mass*emass
            ELSE IF (INDEX(buffer,'charge') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%charge
                curr%charge=curr%charge*echarge
            ELSEIF (INDEX(buffer,'nppcell') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),'(i10)') curr%nppcell
            ELSE IF (INDEX(buffer,'x_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%x_min
            ELSE IF (INDEX(buffer,'x_max') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%x_max
            ELSEIF (INDEX(buffer,'y_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%y_min
            ELSE IF (INDEX(buffer,'y_max') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%y_max
            ELSEIF (INDEX(buffer,'z_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%z_min
            ELSE IF (INDEX(buffer,'z_max') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%z_max
            ELSE IF (INDEX(buffer,'vdrift_x') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vdrift_x
            ELSE IF (INDEX(buffer,'vdrift_y') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vdrift_y
            ELSE IF (INDEX(buffer,'vdrift_z') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vdrift_z
            ELSE IF (INDEX(buffer,'end::species') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_species_section

    SUBROUTINE read_output_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS GRID SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            IF (INDEX(buffer,'output_frequency') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') output_frequency
            ELSE IF (INDEX(buffer,'output_step_min') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') output_step_min
            ELSE IF (INDEX(buffer,'output_step_max') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') output_step_max
            ELSEIF (INDEX(buffer,'ex') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_ex
            ELSE IF (INDEX(buffer,'ey') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_ey
            ELSE IF (INDEX(buffer,'ez') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_ez
            ELSEIF (INDEX(buffer,'bx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_bx
            ELSE IF (INDEX(buffer,'by') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_by
            ELSE IF (INDEX(buffer,'bz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_bz
            ELSEIF (INDEX(buffer,'jx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_jx
            ELSE IF (INDEX(buffer,'jy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_jy
            ELSE IF (INDEX(buffer,'jz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_jz
            ELSEIF (INDEX(buffer,'rho') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_rho
            ELSE IF (INDEX(buffer,'dive') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') c_output_dive
            ELSE IF (INDEX(buffer,'end::output') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_output_section

END MODULE control_file

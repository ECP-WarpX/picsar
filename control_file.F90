MODULE control_file

  USE shared_data
  USE params
  USE fields
  USE particles
  IMPLICIT NONE

  INTEGER(idp) :: ios=0
  INTEGER(idp), PARAMETER :: fh_input = 15
  CHARACTER(LEN=string_length) :: buffer
  CHARACTER(LEN=string_length) :: section_name

CONTAINS
    ! Routine that proceeds to default init
    SUBROUTINE default_init
        ! --- Init particle tiling split
        ntilex = 1
        ntiley = 1
        ntilez = 1

        ! --- Order of Maxwell field solver (default is 2 in x,y,z)
        norderx = 2
        nordery = 2
        norderz = 2
        l_nodalgrid = .FALSE.
        ! --- Order of current deposition/ field gathering 
        ! (default is 2 in x,y,z)
        nox = 1
        noy = 1
        noz = 1
        l_lower_order_in_v = .FALSE.

        ! --- sets coefficient multiplying Courant time step
        dtcoef = 0.7_num
        ! --- smoothing
        npass = 0
        alpha = 0.5_num
        ! --- sets max time in the simulation (in 1/w0)
        tmax = 40.0_num

        !-------------------------------------------------------------------------------
        ! plasma parameters (cold plasma)
        l_particles_weight = .FALSE. ! .TRUE. if particles have different weights

        ! --- quantities in plasma (or lab) frame
        !-------------------------------------------------------------------------------
        nlab  = 1.e23_num            ! plasma density in lab frame
        g0    = 130.0_num          ! initial gamma
        b0    = sqrt(1.0_num-1.0_num/g0**2)
        nc    = nlab*g0          ! density (in the simulation frame)
        wlab  = echarge*sqrt(nlab/(emass*eps0)) ! plasma frequency (in the lab frame)
        w0_l  = echarge*sqrt(nc/(g0*emass*eps0))    ! "longitudinal" plasma frequency (in the lab frame)
        w0_t  = echarge*sqrt(nc/(g0**3*emass*eps0)) ! "transverse" plasma frequency (in the lab frame)
        w0    = w0_l
        ! --- Init number of species
        nspecies=0

        ! --- Particle distribution
        pdistr=1

        RETURN
    END SUBROUTINE default_init

    ! Routine that reads command line arguments
    ! Useful for parametric studies
    SUBROUTINE read_from_cl
        INTEGER :: i, ix
        DO i = 1, IARGC()-1,2
            CALL GETARG(i, buffer)
            IF (INDEX(buffer,'ntilex') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') ntilex
            ELSE IF (INDEX(buffer,'ntiley') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') ntiley
            ELSE IF (INDEX(buffer,'ntilez') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') ntilez
            ELSE IF (INDEX(buffer,'distr') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') pdistr
            END IF
        END DO
        RETURN
    END SUBROUTINE read_from_cl

    ! Routine that reads simulation parameters from input file
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
                CASE('section::main')
                    CALL read_main_section
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

    SUBROUTINE read_main_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS GRID SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            IF (INDEX(buffer,'nx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nx_global_grid
                nx_global=nx_global_grid-1
            ELSE IF (INDEX(buffer,'ny') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') ny_global_grid
                ny_global=ny_global_grid-1
            ELSE IF (INDEX(buffer,'nz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nz_global_grid
                nz_global=nz_global_grid-1
            ELSE IF (INDEX(buffer,'ntilex') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') ntilex
            ELSE IF (INDEX(buffer,'ntiley') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') ntiley
            ELSE IF (INDEX(buffer,'ntilez') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') ntilez
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
            ELSE IF (INDEX(buffer,'t_max') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) tmax
            ELSE IF (INDEX(buffer,'end::main') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_main_section

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
        ! minimal init for species attributes
        curr%charge = -echarge
        curr%mass = emass
        curr%nppcell = 0
        curr%x_min = 0._num
        curr%x_max = 0._num
        curr%y_min = 0._num
        curr%y_max = 0._num
        curr%z_min = 0._num
        curr%z_max = 0._num
        curr%vdrift_x =0._num
        curr%vdrift_y =0._num
        curr%vdrift_z =0._num
        curr%vth_x =0._num
        curr%vth_y =0._num
        curr%vth_z =0._num
        curr%species_npart=0
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
                curr%vdrift_x=curr%vdrift_x*clight
            ELSE IF (INDEX(buffer,'vdrift_y') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vdrift_y
                curr%vdrift_y=curr%vdrift_y*clight
            ELSE IF (INDEX(buffer,'vdrift_z') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vdrift_z
                curr%vdrift_z=curr%vdrift_z*clight
            ELSE IF (INDEX(buffer,'vth_x') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vth_x
                curr%vth_x=curr%vth_x*clight
            ELSE IF (INDEX(buffer,'vth_y') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vth_y
                curr%vth_y=curr%vth_y*clight
            ELSE IF (INDEX(buffer,'vth_z') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%vth_z
                curr%vth_z=curr%vth_z*clight
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

MODULE control_file

  USE shared_data
  USE params
  USE fields
  USE particles
  USE params
  USE output_data
  USE time_stat
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
        nxguards=MAX(nox,2)
        nyguards=MAX(noy,2)
        nzguards=MAX(noz,2)
        nxjguards=MAX(nox,2)
        nyjguards=MAX(noy,2)
        nzjguards=MAX(noz,2)
        
        ! Topology
        topology = 0
        
        ! MPI communication
        mpicom_curr = 0
        
        ! Current deposition algorithm
        currdepo = 0
        
        ! Field gathering algorithm 
        fieldgave = 0
        
        ! Sorting activation (not activated by default)       
        sorting_activated = 0
        sorting_dx = 1.
        sorting_dy = 1.
        sorting_dz = 1.
        sorting_shiftx = 0.
        sorting_shifty = 0.
        sorting_shiftz = 0.
         
        ! Time stats output activation 
        timestat_activated = 0
        timestat_period = 0
        
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
        g0    = 1.0_num          ! initial gamma
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
		! Init species array
		IF (.NOT. l_species_allocated) THEN
			nspecies=0
			ALLOCATE(species_parray(1:nspecies_max))
			l_species_allocated=.TRUE.
		ENDIF
		
		
		  ! Temporal output
		  temdiag_frequency = 0
		  temdiag_format = 0
		
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
            ELSE IF (INDEX(buffer,'nprocx') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') nprocx
            ELSE IF (INDEX(buffer,'nprocy') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') nprocy
            ELSE IF (INDEX(buffer,'nprocz') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') nprocz  
            ELSE IF (INDEX(buffer,'nox') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') nox
            ELSE IF (INDEX(buffer,'noy') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') noy
            ELSE IF (INDEX(buffer,'noz') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, '(i10)') noz                  
            ELSE IF (INDEX(buffer,'tmax') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, *) tmax   
            ELSE IF (INDEX(buffer,'nx') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, *) nx_global_grid   
            ELSE IF (INDEX(buffer,'ny') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, *) ny_global_grid   
            ELSE IF (INDEX(buffer,'nz') .GT. 0) THEN
                CALL GETARG(i+1, buffer)
                READ(buffer, *) nz_global_grid                                                                                   
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
                !write(0,*) TRIM(ADJUSTL(section_name))
                SELECT CASE(TRIM(ADJUSTL(section_name)))
                CASE('section::main')
                    CALL read_main_section
                CASE('section::species')
                    CALL read_species_section
                CASE('section::output')
                    CALL read_output_section
                CASE('section::cpusplit')
                    CALL read_cpusplit_section
                CASE('section::plasma')
                    CALL read_plasma_section
                CASE('section::temporal')
                    CALL read_temporal_output_section
                CASE('section::solver')                    
                    CALL read_solver_section
                CASE('section::timestat')                     
                    CALL read_timestat_section
                CASE('section::sorting')      
                    CALL read_sorting_section
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
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF             
            IF (INDEX(buffer,'nprocx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nprocx
            ELSE IF (INDEX(buffer,'nprocy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nprocy
            ELSE IF (INDEX(buffer,'nprocz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nprocz
            ELSE IF (INDEX(buffer,'topology') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') topology
            ELSE IF (INDEX(buffer,'mpicom_curr') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') mpicom_curr
            ELSE IF (INDEX(buffer,'end::cpusplit') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_cpusplit_section

    SUBROUTINE read_plasma_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS CPUSPLIT SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF             
            IF (INDEX(buffer,'nlab') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) nlab
            ELSE IF (INDEX(buffer,'gamma0') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) g0
            ELSE IF (INDEX(buffer,'pdistr') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) pdistr                
            ELSE IF (INDEX(buffer,'end::plasma') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_plasma_section

    SUBROUTINE read_solver_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS CPUSPLIT SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF             
            IF (INDEX(buffer,'norderx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') norderx
            ELSE IF (INDEX(buffer,'nordery') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nordery
            ELSE IF (INDEX(buffer,'norderz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') norderz 
            ELSE IF (INDEX(buffer,'nox') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nox  
            ELSE IF (INDEX(buffer,'noy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') noy
            ELSE IF (INDEX(buffer,'noz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') noz                                                             
             ELSE IF (INDEX(buffer,'currdepo') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') currdepo        
             ELSE IF (INDEX(buffer,'fieldgave') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') fieldgave                                                 
            ELSE IF (INDEX(buffer,'end::solver') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_solver_section

    SUBROUTINE read_sorting_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS CPUSPLIT SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF                 
             IF (INDEX(buffer,'activation') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') sorting_activated 
             ELSE IF (INDEX(buffer,'dx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),*) sorting_dx 
             ELSE IF (INDEX(buffer,'dy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),*) sorting_dy  
             ELSE IF (INDEX(buffer,'dz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),*) sorting_dz 
             ELSE IF (INDEX(buffer,'shiftx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),*) sorting_shiftx                                          
             ELSE IF (INDEX(buffer,'shifty') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),*) sorting_shifty 
             ELSE IF (INDEX(buffer,'shiftz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length),*) sorting_shiftz                                                              
            ELSE IF (INDEX(buffer,'end::sorting') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_sorting_section

    SUBROUTINE read_timestat_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS CPUSPLIT SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF             
            IF (INDEX(buffer,'activation') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') timestat_activated
            ELSE IF (INDEX(buffer,'period') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') timestat_period                                        
            ELSE IF (INDEX(buffer,'end::timestat') .GT. 0) THEN
                end_section =.TRUE.
            END IF
        END DO
        RETURN
    END SUBROUTINE read_timestat_section 

    SUBROUTINE read_main_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.
        ! READS GRID SECTION OF INPUT FILE
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF 
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
            ELSEIF (INDEX(buffer,'xmin') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) xmin
            ELSE IF (INDEX(buffer,'ymin') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) ymin
            ELSE IF (INDEX(buffer,'zmin') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) zmin
            ELSEIF (INDEX(buffer,'xmax') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) xmax
            ELSE IF (INDEX(buffer,'ymax') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) ymax
            ELSE IF (INDEX(buffer,'zmax') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) zmax
            ELSE IF (INDEX(buffer,'t_max') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) tmax
            ELSE IF (INDEX(buffer,'nguardsx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nxguards
             ELSE IF (INDEX(buffer,'nguardsy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nyguards   
            ELSE IF (INDEX(buffer,'nguardsz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nzguards
            ELSE IF (INDEX(buffer,'njguardsx') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nxjguards
             ELSE IF (INDEX(buffer,'njguardsy') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nyjguards   
            ELSE IF (INDEX(buffer,'njguardsz') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') nzjguards                                                                         
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
        curr%sorting_period = 0
        curr%species_npart=0
        end_section=.FALSE.
        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF
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
            ELSE IF (INDEX(buffer,'sorting_period') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), *) curr%sorting_period             
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
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF
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

    SUBROUTINE read_temporal_output_section
        INTEGER :: ix = 0
        LOGICAL :: end_section = .FALSE.    

        DO WHILE((.NOT. end_section) .AND. (ios==0))
            READ(fh_input, '(A)', iostat=ios) buffer
            !WRITE(0,*),TRIM(ADJUSTL(buffer))
            IF (INDEX(buffer,'#') .GT. 0) THEN
               CYCLE
            ENDIF
            IF (INDEX(buffer,'frequency') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_frequency
            ELSE IF (INDEX(buffer,'format') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_format     
            ELSE IF (INDEX(buffer,'kinE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(1) 
            ELSE IF (INDEX(buffer,'exE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(2)  
            ELSE IF (INDEX(buffer,'eyE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(3)  
            ELSE IF (INDEX(buffer,'ezE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(4)
            ELSE IF (INDEX(buffer,'bxE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(5)  
            ELSE IF (INDEX(buffer,'byE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(6)  
            ELSE IF (INDEX(buffer,'bzE') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(7)     
            ELSE IF (INDEX(buffer,'divE-rho') .GT. 0) THEN
                ix = INDEX(buffer, "=")
                READ(buffer(ix+1:string_length), '(i10)') temdiag_act_list(8)                            
            ENDIF            
        ENDDO    
        RETURN
    END SUBROUTINE

    SUBROUTINE init_species_section
        ! INIT SPECIES SECTION 
        IF (.NOT. l_species_allocated) THEN
            nspecies=0
            ALLOCATE(species_parray(1:nspecies_max))
            l_species_allocated=.TRUE.
        ENDIF
    END SUBROUTINE init_species_section
    
    
    
END MODULE control_file

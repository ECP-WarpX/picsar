import os, sys
import shutil

# The two following classes are inspired by the parser in
# picsar/utils/forthon_parser.

class Subroutines( object ):
    def __init__( self, file ):

        f = open('./src/%s'%file,"r")
        listlines=f.readlines()

        names=[]
        istart=[]
        # Init iend with -1 to avoid code duplication, this element is removed
        # afterwards
        iend=[-1]
        proceduremod=[]
        Nlines=len(listlines)
        proc_mod=True
        modulename=""
        for i in range(0,Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")
            # We are in a module
            if (("module" in curr_word_list) & (curr_line.find("use")==-1) \
                & (curr_line.find("!")==-1)    \
                & (curr_line.find("end module")==-1)):
                modulename=curr_word_list[curr_word_list.index("module")+1]
            # We are in a procedure module
            if(curr_line.find("contains")>=0 ):
                proc_mod=True
            # We are out of a module
            if ((curr_line.find("end module")>=0)):
                proc_mod=False
                modulename=""
            # This is a subroutine block
            if ((curr_line.find("subroutine")>=0) \
                & (curr_line.find("!")==-1)    \
                & (curr_line.find("end subroutine")==-1)):
                ip=curr_line.find("(")

                # Do not compt the subroutines which are in interface
                if i > iend[-1]:
                    # This is a subroutine in a procedure module
                    if (proc_mod):
                        proceduremod.append(modulename)
                    else:
                        proceduremod.append("")
                    if (ip == -1): # Case no parenthesis
                        ip=len(curr_line)
                    subname=curr_line[curr_line.find("subroutine")+\
                                                         len("subroutine"):ip]
                    names.append(subname.strip())
                    istart.append(i)
                    search_end=True
                    while True: # get end of block
                        i=i+1

                        if (i>=Nlines):
                            sys.exit("ERROR: missing end subroutine block")
                        curr_line=listlines[i].lower()

                        # If an interface is defined in the subroutine,
                        # we pass the interface block without looking for
                        # "end subroutine".
                        if (curr_line.find("interface")>=0):
                          search_end = False
                        if (curr_line.find("end interface")>=0):
                          search_end = True
                        if ((curr_line.find("end subroutine")>=0)and(search_end)):
                            break
                    iend.append(i)

        # Remove the -1 value of iend
        iend.remove(-1)

        # Copy also the docstring, we start from each istart and we go back
        # until we meet the string '!_______'
        for (i, compt) in zip(istart, range(len(istart))):
            if '! ____________________________' in listlines[i-1]:
                for j in range(i-2):
                    j = i-2-j
                    curr_line=listlines[j]
                    if '! ____________________________' in curr_line:
                        istart[compt] = j
                        break
            elif '! ____________________________' in listlines[i-2]:
                for j in range(i-3):
                    j = i-3-j
                    curr_line=listlines[j]
                    if '! ____________________________' in curr_line:
                        istart[compt] = j
                        break


        self.file         = file
        self.names        = names
        self.istart       = istart
        self.iend         = iend
        self.proceduremod = proceduremod
        self.isroutine    = True

class Modules( object ):
    def __init__( self, file ):

        f = open('./src/%s'%file,"r")
        listlines=f.readlines()

        names=[]
        istart=[]
        iend=[]
        proceduremod=[]
        Nlines=len(listlines)
        for i in range(0,Nlines):
            curr_line=listlines[i].lower()
            curr_line=self.rm_comments(curr_line)
            curr_word_list=curr_line.split(" ")
            # This is a module block
            proc_mod=False
            if (("module" in curr_word_list) & (curr_line.find("use ")==-1)   \
               & (curr_line.find("!")==-1)    \
               & (curr_line.find("end module")==-1) ):
                subname=curr_line[                                            \
                         curr_line.find("module")+len("module"):len(curr_line)]
                istart.append(i)
                while True: # get end of block
                    i=i+1
                    if (i>=Nlines):
                        sys.exit("ERROR: missing end module block statement")
                    curr_line=listlines[i].lower()
                    curr_word_list=curr_line.split(" ")
                    if (("type" in curr_word_list) &                          \
                       (not "end" in curr_word_list)):
                        typename=curr_word_list[curr_word_list.index("type")+1]
                        # This is a module type
                        if (subname.find(typename) >=0):
                            subname=subname+"moduletype"
                    if ((curr_line.find("end module")>=0)):
                        break
                    if ((curr_line.find("contains")>=0)):
                        proc_mod=True

                iend.append(i)
                names.append(subname.strip())
                proceduremod.append(proc_mod)

        # Copy also the docstring, we start from each istart and we go back
        # until we meet the string '!_______'
        for (i, compt) in zip(istart, range(len(istart))):
            if '! ____________________________' in listlines[i-1]:
                for j in range(i-2):
                    j = i-2-j
                    curr_line=listlines[j]
                    if '! ____________________________' in curr_line:
                        istart[compt] = j
                        break
            elif '! ____________________________' in listlines[i-2]:
                for j in range(i-3):
                    j = i-3-j
                    curr_line=listlines[j]
                    if '! ____________________________' in curr_line:
                        istart[compt] = j
                        break

        self.file         = file
        self.names        = names
        self.istart       = istart
        self.iend         = iend
        self.proceduremod = proceduremod
        self.isroutine    = False

    ### Remove comments from current line
    def rm_comments(self, line):
        icomm=line.find("!")
        icomm_defined=line.find("defined")
        icomm_omp=line.find("!$omp")
        #Problen in the parser with intel directive
        #icomm_intel=line.find("!DIR$")
        icomm_intel=-1
        icomm_ibm=line.find("!ibm*")
        iwrapinstr=line.find("!#do not wrap")
        if ((icomm >=0) & \
            (iwrapinstr==-1) & \
            (icomm_defined==-1) & \
            (icomm_omp==-1) & \
            (icomm_intel==-1) &\
             (icomm_ibm==-1)):
            linecopy=line[0:icomm]
        else:
            linecopy=line
        return linecopy

class MiniAppParser( object ):

    def __init__( self, type_solver, type_pusher, type_depos ):

        self.clean_folder()

        # Create the folder for the mini app
        os.makedirs('./PICSARlite')

        # To be completed
        # Solver
	generic_modules = [
                           "diagnostics",\
                           "control_file",\
                           "simple_io",\
                           "PICSAR_precision",\
                           "constants",\
                           "precomputed",\
                           "params",\
                           "mpi_type_constants",\
                           "communications",\
                           "time_stat",\
                           "output_data",\
                           "shared_data",\
                           "mem_status",\
                           "mpi_derived_types",\
                           "mpi_routines" ]
        generic_routines = [
		           "calc_diags",\
		           "calc_field_div",\
		           "calc_field_divB",\
		           "init_diags",\
		           "init_temp_diags",\
		           "init_time_stat_output",\
		           "get_loc_kinetic_energy",\
		           "get_kinetic_energy",\
		           "get_loc_field_energy_2d",\
		           "get_loc_field_energy",\
		           "get_field_energy_2d",\
		           "get_field_energy",\
		           "get_field_energy",\
		           "get_loc_norm_2",\
		           "get_norm_divErho",\
			   "output_routines",\
			   "output_temporal_diagnostics",\
			   "write_3d_field_array_to_file",\
			   "write_single_array_to_file",\
			   "mpi_minimal_init",\
			   "setup_communicator",\
			   "mpi_initialise",\
			   "compute_simulation_axis",\
			   "allocate_grid_quantities",\
			   "mpi_close",\
			   "time_statistics",\
			   "time_statistics_per_iteration",\
			   "get_local_grid_mem ",\
			   "get_global_grid_mem",\
  			   "set_tile_split"
  			   "set_tile_split_for_species"
  			   "add_particle_to_species_2d"
  			   "add_particle_to_species"
  			   "add_particle_at_tile_2d"
  			   "add_particle_at_tile"
  			   "add_group_of_particles_at_tile"
  			   "rm_particles_from_species_with_mask"
  			   "rm_particles_from_species_2d"
  			   "rm_particles_from_species"
  			   "rm_particle_at_tile_2d"
  			   "rm_particle_at_tile"
  			   "allocate_tile_arrays"
  			   "init_tile_arrays"
  			   "init_tile_arrays_for_species"
  			   "load_particles"
  			   "resize_particle_arrays"
  			   "resize_1D_array_real"
  			   "resize_2D_array_real"
  			   "resize_3D_array_real"
  			   "get_local_number_of_part"
  			   "point_to_tile"
  			   "set_particle_species_properties"
  			   "get_are_tiles_reallocated"
  			   "set_are_tiles_reallocated"
  			   "estimate_tiles_memory_consumption"
  			   "load_laser_species"
  			   "load_laser"
  			   "product_matrix_2c2"
  			   "get_local_tile_mem"
  			   "get_global_tile_mem"
			   "default_init"
			   "read_from_cl"
			   "read_input_file"
			   "read_cpusplit_section"
			   "read_plasma_section"
			   "read_solver_section"
			   "read_sorting_section"
			   "read_timestat_section"
			   "read_main_section"
			   "read_species_section"
			   "read_particle_dumps_section"
			   "read_output_section"
			   "read_temporal_output_section"
			   "read_antenna_section"
			   "init_species_section"
			   "initall"




]





        generic_modules_solver = ["fields","field_boundary"]
        generic_routines_solver = [
  "field_bc"
  "exchange_mpi_3d_grid_array_with_guards"
  "exchange_mpi_3d_grid_array_with_guards_nonblocking"
  "summation_bcs"
  "summation_bcs_nonblocking"
  "summation_bcs_persistent_jx"
  "summation_bcs_persistent_jy"
  "summation_bcs_persistent_jz"
  "efield_bcs"
  "bfield_bcs"
  "current_bcs"
  "charge_bcs"
]
        generic_routines_pusher = []
        generic_modules_pusher = [
			 "grid_tilemodule",\
			 "buff_exchange_part",\
			 "antenna",\
			 "particle_speciesmodule",\
			 "particle_properties",\
			 "particles",\
			 "buff_exchange_part",\
			 "sorting",\
			 "particle_boundary" ]

        generic_modules_depos = [
			"grid_tilemodule",\
			"tile_params",\
			"tiling",\
			"particle_tilemodule"]
        generic_routines_depos = []
	modules_spectral = [
			"math_tools",\
			"gpstd_solver",\
			"fourier_psaotd",\
			"fftw3_fortran",\
			"mpi_fftw3",\
			"fastfft",\
			"matrix_data",\
			"matrix_coefficients",\
			"fourier",\
			"group_parameters",\
			"load_balance"]
	routine_spectral  = [
    mpi_routines.F90:END SUBROUTINE get_non_periodic_mpi_bcs
    mpi_routines.F90:END SUBROUTINE setup_groups
      mpi_routines.F90:END SUBROUTINE adjust_grid_mpi_global
    mpi_routines.F90:  END SUBROUTINE mpi_minimal_init_fftw
          END SUBROUTINE generalized_comms_group_l2g
          END SUBROUTINE  sendrecv_l2g_generalized
          END SUBROUTINE  sendrecv_l2g_generalized_non_blocking
          END SUBROUTINE generalized_comms_group_g2l
          END SUBROUTINE  sendrecv_g2l_generalized_non_blocking
          END SUBROUTINE  sendrecv_g2l_generalized

	END SUBROUTINE field_damping_bcs
	END SUBROUTINE merge_fields
	END SUBROUTINE merge_e_fields
	END SUBROUTINE merge_b_fields
  END SUBROUTINE push_psatd_ebfield
	init_pml_arrays
  END SUBROUTINE select_case_dims_local
  END SUBROUTINE select_case_dims_global
  END SUBROUTINE init_kspace
  END SUBROUTINE delete_k_space
  END SUBROUTINE compute_k_vec
  END SUBROUTINE compute_k_1d
  END SUBROUTINE fftfreq
  END SUBROUTINE init_gpstd
  END SUBROUTINE compute_cc_mat_splitted_fields
  END SUBROUTINE compute_cc_mat_merged_fields
  END SUBROUTINE FD_weights_hvincenti
  END SUBROUTINE copy_field
  END SUBROUTINE copy_field_forward
  END SUBROUTINE copy_field_backward
  END SUBROUTINE init_plans_fourier_mpi
  END SUBROUTINE get_Ffields
  END SUBROUTINE get_Ffields_mpi_lb 
  END SUBROUTINE get_Ffields_mpi
  END SUBROUTINE get_fields
  END SUBROUTINE get_fields_mpi
  END SUBROUTINE get_fields_mpi_lb
  END SUBROUTINE fft_forward_r2c_local
  END SUBROUTINE fft_forward_r2c_hybrid
  END SUBROUTINE fft_backward_c2r_local
  END SUBROUTINE fft_backward_c2r_hybrid
  END SUBROUTINE push_psaotd_ebfielfs_2d
  END SUBROUTINE push_psaotd_ebfielfs_3d
  END SUBROUTINE init_plans_blocks

fast_fftw_create_plan_r2c_3d_dft
 fast_fftw_create_plan_c2r_3d_dft
 fast_fftw_create_plan_r2c_2d_dft
 fast_fftw_create_plan_c2r_2d_dft
 fast_fftw3d_c2r_with_plan
 fast_fftw3d_r2c_with_plan
allocate_new_matrix_vector
multiply_mat_vector
multiply_unit_blocks




]
	fdtd_routines= [
END SUBROUTINE push_bfield
        END SUBROUTINE push_efield
END SUBROUTINE push_bfield_2d
END SUBROUTINE push_efield_2d
END SUBROUTINE pxrpush_em3d_evec_norder
END SUBROUTINE pxrpush_em2d_evec_norder
END SUBROUTINE pxrpush_em2d_evec
END SUBROUTINE pxrpush_em3d_evec
END SUBROUTINE pxrpush_em3d_bvec_norder
END SUBROUTINE pxrpush_em2d_bvec_norder
END SUBROUTINE pxrpush_em2d_bvec
END SUBROUTINE pxrpush_em3d_bvec


]

        if type_solver == 'all':
            self.list_available_modules = []
            self.list_available_routines = []

        elif type_solver == 'fdtd':
            self.list_available_modules = []
            self.list_available_routines = []

        elif type_solver == 'spectral':
            self.list_available_modules = []
            self.list_available_routines = []

        else:
            print('#########################################################' \
                  '#########')
            print('Wrong solver argument there. The list of available '       \
                  'solver is:')
            print('- all: every routines and modules from picsar')
            print('- fdtd: every routines and modules related to the Finate ' \
                  'Domain')
            print('        Time Domain Maxwell solver in 3D.')
            print('- spectral: every routines and modules related to the '    \
                  'Pseudo-')
            print('        Spectral Analytical Time Domain Maxwell solver in '\
                  '3D.')
            print('#########################################################' \
                              '#########')

        # Pusher
        if type_pusher == 'all':
            self.list_available_modules.append([])
            self.list_available_routines.append([])

        elif type_pusher == 'boris':
            self.list_available_modules.append([])
            self.list_available_routines.append([])

        elif type_pusher == 'vay':
            self.list_available_modules.append([])
            self.list_available_routines.append([])

        else:
            print('#########################################################' \
                  '#########')
            print('Wrong pusher argument there. The list of available '       \
                  'pusher is:')
            print('- all')
            print('- boris')
            print('- vay')
            print('#########################################################' \
                              '#########')

        # Deposition
        if type_depos == 'all':
            self.list_available_modules.append([])
            self.list_available_routines.append([])

        elif type_depos == 'direct':
            self.list_available_modules.append([])
            self.list_available_routines.append([])

        elif type_depos == 'esirkepov':
            self.list_available_modules.append([])
            self.list_available_routines.append([])

        else:
            print('#########################################################' \
                  '#########')
            print('Wrong particle deposition argument there. The list of '    \
                  'available pusher is:')
            print('- all')
            print('- direct')
            print('- esirkepov')
            print('#########################################################' \
                              '#########')

        self.list_available_routines = None
        self.list_available_modules = None

        #LIST ALL .F90 or .F files in current directory
        self.listfiles = self.create_listfiles('./src')

        # Reconstruct PICSARlite
        for file in self.listfiles:
            self.write_available_routines(file)

        # Remove unavailable routines
        self.availablelistfiles = self.create_listfiles('./PICSARlite/src')
        for file in self.availablelistfiles:
            self.comment_unavailable_routine(file)

        # Copy some extra needed folders
        self.copy_utils_example()

    def clean_folder(self):

        # Delete all files
        if os.path.exists('./PICSARlite'):
            shutil.rmtree('./PICSARlite')

    def create_listfiles(self, folder):
        listfiles = []
        for root, subFolder, files in os.walk(folder):
            for file in files:
                if len(root) > len(folder):
                     listfiles.append('%s/%s'%(root[len(folder)+1:], file))
        return listfiles

    def write_available_routines(self, file):

        r = Subroutines(file)
        m = Modules(file)

        # Check if subroutines are inside modules and get the number of
        # different modules there are
        nb_modules = list(set(r.proceduremod))

        try: nb_modules.remove('')
        except (ValueError):  pass

        # Remove the unnecessary strin '\n' at the end of each modules
        for imodule in range(len(nb_modules)):
            nb_modules[imodule] = nb_modules[imodule][:-1]

        for module in m.names:
            # Loop on each module occurence, for each print the modules then
            # subroutines inside and finish by "end module"
            if module in nb_modules:
                if (self.list_available_modules is None):
                    iname = m.names.index(module)
                    self.copy_files_from_picsar(m, iname)

                    # Loop on the available routines
                    for iname in range(len(r.names)):
                        if module == r.proceduremod[iname][:-1]:
                            if self.list_available_routines is None:
                                self.copy_files_from_picsar(r, iname)

                            elif r.names[iname] in self.list_available_routines:
                                self.copy_files_from_picsar(r, iname)

                    # End the module
                    self.end_module(r, module)

                elif (module in self.list_available_modules):
                    iname = m.names.index(module)
                    self.copy_files_from_picsar(m, iname)

                    # Loop on the available routines
                    for iname in range(len(r.names)):
                        if module == r.proceduremod[iname][:-1]:
                            if self.list_available_routines is None:
                                self.copy_files_from_picsar(r, iname)

                            elif r.names[iname] in self.list_available_routines:
                                self.copy_files_from_picsar(r, iname)

                    # End the module
                    self.end_module(r, module)
            else:
                    iname = m.names.index(module)
                    self.copy_files_from_picsar(m, iname)

                    if self.list_available_modules is None:
                        self.copy_files_from_picsar(m, iname)

                    elif m.names[iname] in self.list_available_modules:
                        self.copy_files_from_picsar(m, iname)

        # Add the subroutines out of modules
        for iname in range(len(r.names)):
            if r.proceduremod[iname] == '':
                if self.list_available_routines is None:
                    self.copy_files_from_picsar(r, iname)

                elif r.names[iname] in self.list_available_routines:
                    self.copy_files_from_picsar(r, iname)

    def copy_files_from_picsar( self, routine, index ):
        # routine is an instance of the class Routines or Modules

        file = routine.file

        if not os.path.exists('./PICSARlite/src/%s'%file):

            divided_path = file.split('/')
            path_folder = './PICSARlite/src/'

            for folder in divided_path[:-1]:
                path_folder += '%s/'%folder
            if not os.path.exists(path_folder):
                os.makedirs(path_folder)

            # Write the description of the file
            f = open('./src/%s'%file, 'r')
            lines = f.readlines()

            # Copy the file header
            text_header = ''
            for line in lines:
                text_header += line

            split_str = '! ________________________________________________'\
                        '______________________________'

            # Cut the beginning if needed
            text_header_temp = text_header.split(split_str)[1].split('!')
            text_header = ''
            for line in text_header_temp[1:]:
                text_header += '!' + line

            # Remove the list of routines
            text_header = text_header.split('! List of subroutines:')[0]

            fnew = open('./PICSARlite/src/%s'%file, 'w')
            fnew.writelines(split_str)
            fnew.writelines(text_header)
            fnew.writelines(split_str+'\n\n')
            fnew.close()

        # Write the subroutine inside the file
        f = open('./src/%s'%file, 'r')
        lines = f.readlines()
        routine_txt = ''
        for line in lines[routine.istart[index]:routine.iend[index]+1]:
            routine_txt += line

        fnew = open('PICSARlite/src/%s'%file, 'aw')
        fnew.writelines(routine_txt +'\n')
        fnew.close()

    def end_module( self, routine, module ):
        file = routine.file

        f = open('./PICSARlite/src/%s'%file, 'aw')
        f.writelines('END MODULE %s'%module +'\n')
        f.close()

    def comment_unavailable_routine(self, file):
        """
            When all the files are written, check inside the files if some
            unavailable routines are unused but called. It may introduce
            a crash at compilation time.

            To avoid it, this function comments these routines and inserts
            stops.
        """


        def formatting_line(nb_blanks, str, add_ampersand=True):
            "Add the & symbol at the 87th place."
            L = len(str)+nb_blanks
            space = ' '
            strnew = ''
            # Add nb_blanks at the beginning
            for n in range(nb_blanks):
                strnew += space

            strnew += str

            for n in range(86-L):
                strnew += space
            if add_ampersand:
                strnew += '&'
            return strnew

        # Find the unwanted call
        f = open('./PICSARlite/src/%s'%file,"r")
        listlines=f.readlines()
        Nlines = len(listlines)
        istart = [-1]
        iend   = [-1]
        list_routine_name = []
        nb_blanks = []

        for i in range(0, Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")
            # We find a CALL
            if (("call" in curr_word_list) & (curr_line.find("!")==-1 )):
                # Get the name of the following routine
                indexcall = curr_word_list.index("call")
                routine_name = curr_word_list[indexcall+1].split('(')[0]

                # If the routine is in the list, everything fine, if not the
                # line and the block should be commented.
                if not routine_name in self.list_available_routines:
                    list_routine_name.append(routine_name)
                    istart.append(i)
                    curr_line_old = curr_line

                    # Store the number of blanks before the call
                    nb_blanks.append(indexcall)

                    while True:
                        if (i>=Nlines):
                            sys.exit("ERROR: missing end call block")
                        if not "&" in curr_line_old:
                            break
                        else:
                            curr_line_old = curr_line
                            i += 1
                            curr_line=listlines[i].lower()
                    iend.append(i)

        # Comment and stop the code
        fnew = open('./PICSARlite/src/%s'%file,"w")
        listlines_new = []
        compt = 0
        for i in range(0, Nlines):
            # If no more wrong subroutines, finsh the file
            if (compt == len(iend)-1) :
                listlines_new.append(listlines[i])

            elif i >= iend[compt]:
                # It is regular lines
                if i != istart[compt+1]:
                    listlines_new.append(listlines[i])

                # It is a call block
                else:
                    # Comment the lines
                    for iblock in range(iend[compt+1]-istart[compt+1]):
                        listlines_new.append('!'+listlines[i+iblock][1:])

                    # Print error message
                    error_message = \
                formatting_line( nb_blanks[compt], \
                                       "WRITE(0,*) 'ERROR:',") + '\n' \
              + formatting_line( nb_blanks[compt]+11, \
                    "'The subroutine %s', "%(list_routine_name[compt])) + '\n'\
              + formatting_line( nb_blanks[compt]+11, \
                                "'cannot be used in this configuration.'",  \
                                add_ampersand=False)

                    listlines_new.append("\n" + error_message + "\n")
                    listlines_new.append(formatting_line( nb_blanks[compt], \
                                        "STOP", add_ampersand=False)+ "\n \n")
                    if compt == len(iend)-2:
                        compt = 0
                    else:
                        compt += 1

        fnew.writelines(listlines_new)
        fnew.close()

    def copy_utils_example(self):
        os.system('cp -r ./utils ./PICSARlite/utils')
        os.system('cp -r ./examples ./PICSARlite/examples')

###############################################################################
# Main

arglist=sys.argv
#try:
type_parser = "all" #str(arglist[1])
type_pusher = "all" #str(arglist[2])
type_depos  = "all" #str(arglist[3])
miniapp = MiniAppParser(type_parser, type_pusher, type_depos)

# except(IndexError):
#     print('##################################################################')
#     print('Some arguments are needed there.')
#     print('The list of available solvers is:')
#     print('- all: every routines and modules from picsar')
#     print('- fdtd: every routines and modules related to the Finate Domain')
#     print('        Time Domain Maxwell solver in 3D.')
#     print('- spectral: every routines and modules related to the Pseudo-')
#     print('        Spectral Analytical Time Domain Maxwell solver in 3D.')
#     print('##################################################################')

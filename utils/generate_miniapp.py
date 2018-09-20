import os, sys
import shutil
from datetime import datetime

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
                            sys.exit("ERROR: missing end subroutine block in "+file)
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
        self.type_solver = type_solver
        self.type_pusher = type_pusher
        self.type_depos  = type_depos

        # Create the folder for the mini app
        os.makedirs('./PICSARlite')

        # To be completed
        # Solver
	generic_modules = ["diagnostics",\
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

        generic_routines=["calc_diags",\
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
  			    "set_tile_split",\
  			    "set_tile_split_for_species",\
  			    "add_particle_to_species_2d",\
  			    "add_particle_to_species",\
  			    "add_particle_at_tile_2d",\
  			    "add_particle_at_tile",\
  			    "add_group_of_particles_at_tile",\
  			    "rm_particles_from_species_with_mask",\
  			    "rm_particles_from_species_2d",\
  			    "rm_particles_from_species",\
  			    "rm_particle_at_tile_2d",\
  			    "rm_particle_at_tile",\
  			    "allocate_tile_arrays",\
  			    "init_tile_arrays",\
  			    "init_tile_arrays_for_species",\
  			    "load_particles",\
  			    "resize_particle_arrays",\
  			    "resize_1D_array_real",\
  			    "resize_2D_array_real",\
  			    "resize_3D_array_real",\
  			    "get_local_number_of_part",\
  			    "point_to_tile",\
  			    "set_particle_species_properties",\
  			    "get_are_tiles_reallocated",\
  			    "set_are_tiles_reallocated",\
  			    "estimate_tiles_memory_consumption",\
  			    "load_laser_species",\
  			    "load_laser",\
  			    "product_matrix_2c2",\
  			    "get_local_tile_mem",\
  			    "get_global_tile_mem",\
			    "default_init",\
			    "read_from_cl",\
			    "read_input_file",\
			    "read_cpusplit_section",\
			    "read_plasma_section",\
			    "read_solver_section",\
			    "read_sorting_section",\
			    "read_timestat_section",\
			    "read_main_section",\
			    "read_species_section",\
			    "read_particle_dumps_section",\
			    "read_output_section",\
			    "read_temporal_output_section",\
			    "read_antenna_section",\
			    "init_species_section",\
			    "initall",\
			    "pxr_gete2dxz_energy_conserving_vect_1_1",\
			    "pxr_getb2dxz_energy_conserving_vect_1_1",\
			    "pxr_gete2dxz_energy_conserving_vect_2_2",\
			    "pxr_getb2dxz_energy_conserving_vect_2_2",\
			    "pxr_gete2dxz_energy_conserving_vect_2_2",\
			    "pxr_getb2dxz_energy_conserving_vect_2_2",\
			    "pxr_gete2dxz_energy_conserving_scalar_3_3",\
			    "pxr_gete2dxz_energy_conserving_vect_3_3",\
			    "pxr_getb2dxz_energy_conserving_scalar_3_3",\
			    "pxr_getb2dxz_energy_conserving_vect_3_3",\
			    "pxr_geteb2dxz_energy_conserving_vect_3_3",\
			    "pxr_gete2dxz_n_energy_conserving",\
			    "pxr_getb2dxz_n_energy_conserving",\
			    "gete3d_energy_conserving_scalar_1_1_1",\
			    "getb3d_energy_conserving_scalar_1_1_1",\
			    "gete3d_energy_conserving_vec_1_1_1",\
			    "getb3d_energy_conserving_vec_1_1_1",\
			    "geteb3d_energy_conserving_vecV1_1_1_1",\
			    "geteb3d_energy_conserving_vecV2_1_1_1",\
			    "geteb3d_energy_conserving_vecV3_1_1_1",\
			    "geteb3d_energy_conserving_vecV4_1_1_1",\
			    "geteb3d_energy_conserving_vec_1_1_1_v2",\
			    "geteb3d_energy_conserving_vec_1_1_1_sub",\
			    "gete3d_energy_conserving_scalar_2_2_2",\
			    "getb3d_energy_conserving_scalar_2_2_2",\
			    "gete3d_energy_conserving_vec_2_2_2",\
			    "getb3d_energy_conserving_vec_2_2_2",\
			    "geteb3d_energy_conserving_vecV1_2_2_2",\
			    "geteb3d_energy_conserving_vecV2_2_2_2", \
			    "geteb3d_energy_conserving_vecV3_2_2_2", \
			    "geteb3d_energy_conserving_vecV4_2_2_2",\
			    "gete3d_energy_conserving_scalar_3_3_3",\
			    "getb3d_energy_conserving_scalar_3_3_3",\
			    "gete3d_energy_conserving_linear_3_3_3", \
			    "getb3d_energy_conserving_linear_3_3_3",\
			    "gete3d_energy_conserving_vec_3_3_3", \
			    "getb3d_energy_conserving_vec_3_3_3",\
			    "gete3d_energy_conserving_vec2_3_3_3",\
			    "getb3d_energy_conserving_vec2_3_3_3",\
			    "geteb3d_energy_conserving_vec_3_3_3",\
			    "geteb3d_energy_conserving_vecV2_3_3_3",\
			    "geteb3d_energy_conserving_vecV3_3_3_3",\
			    "geteb3d_energy_conserving_blockvec_3_3_3",\
			    "geteb3d_energy_conserving_blockvec2_3_3_3",\
			    "pxrgete3d_n_energy_conserving",\
			    "pxrgetb3d_n_energy_conserving",\
			    "pxr_getb3d_n_energy_conserving",\
			    "pxr_gete3d_n_energy_conserving",\
			    "geteb2dxz_energy_conserving",\
			    "geteb2dxz_energy_conserving_generic",\
			    "field_gathering",\
			    "field_gathering_sub",\
			    "geteb3d_energy_conserving",\
			    "geteb3d_energy_conserving_generic"
			    "mpi_send",
			    "mpi_recv",
			    "mpi_isend",
			    "mpi_irecv",
			    "mpi_sendrecv", 
			    "mpi_wait",
			    "mpi_waitall", 
			    "mpi_type_free",
			    "mpi_start mpi_type_free",
			    "mpi_send_init",
			    "mpi_type_commit",
			    "mpi_type_vector",
			    "mpi_reduce",
			    "mpi_allgather",
			    "mpi_abort",
			    "mpi_bcast mpi_barrier",
			    "mpi_comm_free",
			    "mpi_comm_rank",
			    "mpi_comm_size",
			    "mpi_comm_group", 
			    "mpi_comm_dup",
			    "mpi_comm_split", 
			    "mpi_comm_free",
			    "mpi_cart_rank",
			    "mpi_comm_group",
			    "mpi_cart_shift",
			    "mpi_cart_create",
			    "mpi_dims_create",
			    "mpi_init_thread",
			    "mpi_initialized",
			    "mpi_type_create_subarray",
			    "mpi_type_commit",
			    "mpi_type_vector",
			    "mpi_type_commit",
			    "mpi_type_create_struct",
			    "mpi_type_size"]

        generic_modules_solver = ["fields","field_boundary"]

        generic_routines_solver = [
				   "field_bc",\
 				   "exchange_mpi_3d_grid_array_with_guards",\
 				   "exchange_mpi_3d_grid_array_with_guards_nonblocking",\
 				   "summation_bcs",\
 				   "summation_bcs_nonblocking",\
 				   "summation_bcs_persistent_jx",\
 				   "summation_bcs_persistent_jy",\
 				   "summation_bcs_persistent_jz",\
 				   "efield_bcs",\
 				   "bfield_bcs",\
 				   "current_bcs",\
 				   "charge_bcs"]
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

        generic_routines_pusher = [
				   "set_tile_split",\
				   "particle_bcs",\
				   "particle_bcs_2d",\
				   "particle_bcs_tiles",\
				   "particle_bcs_tiles_openmp",\
				   "particle_bcs_tiles_2d",\
				   "particle_bcs_tiles_2d_openmp",\
				   "particle_bsc_openmp_reordering",\
				   "particle_bcs_mpi_blocking",\
				   "particle_bcs_mpi_non_blocking",\
				   "particle_bcs_mpi_non_blocking_2d",\
				   "particle_bcs_tiles_and_mpi_3d",\
				   "push_laser_particles",\
				   "laserp_pusher_gaussian",\
				   "laserp_pusher_hanning",\
				   "gaussian_profile",\
				   "hanning_profile",\
				   "field_gathering_plus_particle_pusher_sub_2d",\
				   "field_gathering_plus_particle_pusher",\
				   "field_gathering_plus_particle_pusher_sub",\
				   "field_gathering_plus_particle_pusher_cacheblock_sub",\
				   "particle_pusher_sub",\
				   "pxrpush_particles_part1",\
				   "pxrpush_particles_part1_sub",\
				   "pxrpush_particles_part2",\
				   "field_gathering_plus_particle_pusher_1_1_1",\
				   "field_gathering_plus_particle_pusher_2_2_2",\
				   "field_gathering_plus_particle_pusher_3_3_3",\
				   "pxr_pushxyz",\
				   "pxr_push2dxz",\
				   "pxr_pushxz"]


	boris_modules = []
	boris_routines = [
			  "pxr_pushxz",\
			  "pxr_push2dxz",\
			  "pxr_boris_push_u_3d",\
			  "pxr_boris_push_rr_S09_u_3d",\
			  "pxr_boris_push_rr_B08_u_3d",\
			  "pxr_boris_push_rr_LL_u_3d",\
			  "pxr_boris_push_u_3d_block",\
			  "pxr_pushxyz",\
			  "pxr_epush_v",\
			  "pxr_bpush_v",\
			  "pxr_set_gamma"]

	vay_pusher_modules = []

	vay_pusher_routines= ["pxr_ebcancelpush3d" ]

        generic_modules_depos = [
			"grid_tilemodule",\
			"tile_params",\
			"tiling",\
			"particle_tilemodule"]
        generic_routines_depos = [
			"depose_rho_scalar_1_1_1",\
			"depose_rho_scalar_2_2_2",\
			"depose_rho_scalar_3_3_3",\
			"depose_rho_vecSH_1_1_1",\
			"depose_rho_vecNOY_1_1_1",\
			"depose_rho_vecHV_1_1_1",\
			"depose_rho_vecHVv2_1_1_1",\
			"depose_rho_vecHVv2_2_2_2",\
			"depose_rho_vecHVv2_3_3_3",\
			"depose_rho_vecHVv3_3_3_3",\
			"depose_rho_vecHVv4_3_3_3",\
			"pxr_depose_rho_n",\
			"pxr_depose_rho_n_2dxz",\
			"pxr_depose_rhoold_n_2dxz",\
			"pxrdepose_rho_on_grid",\
			"pxrdepose_rho_on_grid_sub_openmp_3d_n",\
			"pxrdepose_rho_on_grid_sub_openmp_3d",\
			"pxrdepose_rho_on_grid_sub_openmp_2d",\
			"pxrdepose_rho_on_grid_sub_openmp_3d_scalar",\
			"pxrdepose_rho_on_grid_sub_openmp_3d_vecto",\
			"depose_jxjyjz",\
			"depose_jxjyjz_generic",\
			"pxrdepose_currents_on_grid_jxjyjz",\
			"depose_jxjyjz_2d depose_jxjyjz_generic_2d",\
	                "pxrdepose_currents_on_grid_jxjyjz_2d",\
			"pxrdepose_currents_on_grid_jxjyjz_2d",\
		        "curr_depo_sub",\
			"func_order",\
			]

	esirkepov_modules=[]

	esirkepov_routines=["pxrdepose_currents_on_grid_jxjyjz_sub_openmp",\
			"pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp",\
			"depose_jxjyjz_esirkepov",\
			"depose_jxjyjz_esirkepov_1_1_1",\
			"depose_jxjyjz_esirkepov_2_2_2",\
			"depose_jxjyjz_esirkepov_3_3_3",\
			"pxr_depose_jxjyjz_esirkepov_n",\
			"pxr_depose_jxjyjz_esirkepov2d_n",\
			"pxr_depose_jxjyjz_esirkepov2d_1_1",\
			"pxr_depose_jxjyjz_esirkepov2d_2_2",\
			"pxr_depose_jxjyjz_esirkepov2d_3_3",\
                        "pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3",\
			"  depose_jxjyjz_esirkepov_2d"]

	direct_modules = []

	direct_routines = ["depose_jxjyjz_scalar_1_1_1",\
			"depose_jxjyjz_vecHVv2_1_1_1",\
			"depose_jxjyjz_vecHV_vnr_1_1_1",\
			"depose_jxjyjz_scalar_2_2_2",\
			"depose_jxjyjz_vecHVv2_2_2_2",\
			"depose_jxjyjz_vecHV_vnr_2_2_2",\
			"depose_jxjyjz_scalar_3_3_3",\
			"depose_jxjyjz_vecHVv3_3_3_3",\
			"depose_jxjyjz_vecHV_vnr_3_3_3",\
			"current_reduction_1_1_1",\
			"current_reduction_2_2_2",\
			"current_reduction_3_3_3"]

	spectral_modules= [
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

	spectral_routines  = [
			"get_non_periodic_mpi_bcs",
			"setup_groups",
			"adjust_grid_mpi_global",
			"mpi_minimal_init_fftw",
			"sendrecv_l2g_generalized",
			"sendrecv_l2g_generalized_non_blocking",
                        "generalized_comms_group_l2g",
			"generalized_comms_group_g2l",
			"sendrecv_g2l_generalized_non_blocking",
			"sendrecv_g2l_generalized",
			"field_damping_bcs",
			"merge_fields",
			"merge_e_fields",
			"merge_b_fields",
			"push_psatd_ebfield",
			"init_pml_arrays",
			"select_case_dims_local",
			"select_case_dims_global",
			"init_kspace",
			"delete_k_space",
			"compute_k_vec",
			"compute_k_1d",
			"fftfreq",
			"init_gpstd",
			"compute_cc_mat_splitted_fields",
			"compute_cc_mat_merged_fields",
			"FD_weights_hvincenti",
			"copy_field",
			"copy_field_forward",
			"copy_field_backward",
			"init_plans_fourier_mpi",
			"get_Ffields",
			"get_Ffields_mpi_lb",
			"get_Ffields_mpi",
			"get_fields",
			"get_fields_mpi",
			"get_fields_mpi_lb",
			"fft_forward_r2c_local",
			"fft_forward_r2c_hybrid",
			"fft_backward_c2r_local",
			"fft_backward_c2r_hybrid",
			"push_psaotd_ebfielfs_2d",
			"push_psaotd_ebfielfs_3d",
			"init_plans_blocks",
			"fast_fftw_create_plan_r2c_3d_dft",
			"fast_fftw_create_plan_c2r_3d_dft",
			"fast_fftw_create_plan_r2c_2d_dft",
			"fast_fftw_create_plan_c2r_2d_dft",
			"fast_fftw3d_c2r_with_plan",
			"fast_fftw3d_r2c_with_plan",
			"allocate_new_matrix_vector",
			"multiply_mat_vector",
			"multiply_unit_blocks"
			"fftw_mpi_local_size_3d" ,
			"fftw_mpi_local_size_3d_transposed" ,
			"fftw_mpi_local_size_2d p3dfft_setup" ,
			"p3dfft_get_dims" ,
			"dfftw_init_threads",
			"fftw_mpi_init fftw_mpi_plan_dft_r2c_3d",
			"fftw_mpi_plan_dft_c2r_3d",
			"fftw_mpi_plan_dft_r2c_2d" ,
			"fftw_mpi_plan_dft_c2r_2d" ,
			"p3dfft_ftran_r2c" ,
			"fftw_mpi_execute_dft_r2c",
			"p3dfft_ftran_c2r" ,
			"fftw_mpi_execute_dft_c2r"]
	fdtd_modules = []

	fdtd_routines = [
		        "push_bfield",
		        "push_efield",
		        "init_stencil_coefficients",
			"FD_weights",
		        "push_bfield_2d",
		        "push_efield_2d",
		        "pxrpush_em3d_evec_norder",
		        "pxrpush_em2d_evec_norder",
		        "pxrpush_em2d_evec",
		        "pxrpush_em3d_evec",
		        "pxrpush_em3d_bvec_norder",
		        "pxrpush_em2d_bvec_norder",
		        "pxrpush_em2d_bvec",
		        "pxrpush_em3d_bvec"]

	self.list_available_modules = generic_modules + generic_modules_solver + generic_modules_pusher + generic_modules_depos 
	self.list_available_routines =  generic_routines + generic_routines_solver + generic_routines_pusher + generic_routines_depos
        if type_solver == 'all':
            self.list_available_modules+= spectral_modules + fdtd_modules
            self.list_available_routines += spectral_routines + fdtd_routines

        elif type_solver == 'fdtd':
            self.list_available_modules+= fdtd_modules
            self.list_available_routines+= fdtd_routines

        elif type_solver == 'spectral':
            self.list_available_modules += spectral_modules
            self.list_available_routines += spectral_routines 

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
            self.list_available_modules +=  vay_pusher_modules + boris_modules
            self.list_available_routines +=  vay_pusher_routines + boris_routines

        elif type_pusher == 'boris':
            self.list_available_modules += boris_modules
            self.list_available_routines += boris_routines

        elif type_pusher == 'vay':
            self.list_available_modules += vay_pusher_modules
            self.list_available_routines +=vay_pusher_routines

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
            self.list_available_modules += direct_modules + esirkepov_modules
            self.list_available_routines += direct_routines + esirkepov_routines

        elif type_depos == 'direct':
            self.list_available_modules += direct_modules
            self.list_available_routines += direct_routines

        elif type_depos == 'esirkepov':
            self.list_available_modules += esirkepov_modules
            self.list_available_routines +=esirkepov_routines

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
	print "list of available routines  ",self.list_available_routines
	print "list of available  modules  ",self.list_available_modules

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
        self.generate_main()

    def clean_folder(self):

        # Delete all files
        if os.path.exists('./PICSARlite'):
            shutil.rmtree('./PICSARlite')

    def create_listfiles(self, folder):
        listfiles = []
        for root, subFolder, files in os.walk(folder):
            for file in files:
                if len(root) > len(folder)-5:
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
                #if not routine_name in self.list_available_routines:
                if 1:
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
                                       "WRITE(0,*) 'ABORT.',") + '\n' \
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


    def generate_main( self ):
        file = 'main.F90'
        fnew = open('./PICSARlite/src/%s'%file,"w")
        listlines = []

        # Get the current date
        now = datetime.now()

        # Print the header
        header = '! ________________________________________________________' \
           + '________________________________ \n! \n'                        \
           + '! *** Copyright Notice *** \n! \n'                              \
           + '! "Particle In Cell Scalable Application Resource (PICSAR) v2",'\
           + ' Copyright (c) 2016, \n! The Regents of the University of '     \
           + 'California, through Lawrence Berkeley National \n! Laboratory ' \
           + '(subject to receipt of any required approvals from the U.S. '   \
           + 'Dept. of Energy). \n! All rights reserved. \n! \n'              \
           + '! If you have questions about your rights to use or distribute '\
           + "this software, \n! please contact Berkeley Lab's Innovation & " \
           + 'Partnerships Office at  IPO@lbl.gov. \n! \n! NOTICE. \n! '      \
           + 'This Software was developed under funding from the U.S.'        \
           + ' Department of Energy \n! and the U.S. Government consequently' \
           + 'retains certain rights. As such, the U.S. \n! Government has '  \
           + 'been granted for itself and others acting on its behalf a '     \
           + 'paid-up, \n! nonexclusive, irrevocable, worldwide license in '  \
           + 'the Software to reproduce, distribute \n! copies to the '       \
           + 'public, prepare derivative works, and perform publicly and '    \
           + 'display \n! publicly, and to permit other to do so. \n! \n'     \
           + '! PICSAR MINIAPP \n! \n'                                        \
           + '! Creation date: %s/%s/%s \n! \n'%(now.day, now.month, now.year)\
           + '! DEVELOPERS: \n! - Henri Vincenti \n! - Mathieu Lobet \n! - '  \
           + 'Remi Lehe \n! - Jean-Luc Vay \n! - Guillaume Blaclard \n! - '   \
           + 'Haithem Kallala \n! \n'                                         \
           + '! INCLUDES \n! - Maxwell solver: %s \n'%(self.type_solver)      \
           + '! - Particle pusher: %s \n'%(self.type_pusher)                  \
           + '! - Particle deposition: %s \n'%(self.type_depos)  + '! \n'     \
           + '! ____________________________________________________________' \
           + '____________________________ \n\n\n'
        listlines.append(header)


        # Print the main file
        listlines.append('PROGRAM main\n')
        listmodulesmain = ["PICSAR_precision", "constants", "fields",         \
                        "particles", "params", "shared_data", "mpi_routines", \
                        "control_file", "time_stat", "diagnostics"]

        for module in listmodulesmain:
            #if module in self.list_available_modules:
            if 1:
                listlines.append('  USE %s \n'%module)

        #if "mem_status" in self.list_available_modules:
        if 1:
            listlines.append('  USE mem_status, ONLY : global_grid_mem, '     \
                            +'global_grid_tiles_mem, global_part_tiles_mem\n\n')

        if (self.type_solver == 'all') | (self.type_solver == 'spectral'):
            listlines.append('#if defined(FFTW)\n  USE mpi_fftw3 \n  '       \
             +'USE fourier \n  USE fastfft \n  USE fftw3_fortran \n'        \
             +'#endif\n\n')

        listlines.append('  IMPLICIT NONE\n')
        listlines.append('  LOGICAL :: exist\n')
        listlines.append('  CHARACTER(len=250) :: str1, str2, str3\n')
        listlines.append('  CHARACTER(len=250) :: str4, str5, str7\n')
        listlines.append('  CHARACTER(len=500) :: str6\n\n')
        listlines.append('! --- default init\n')
        listlines.append('  CALL default_init\n\n')
        listlines.append('! --- reads input_file\n')
        listlines.append('  CALL read_input_file\n\n')
        listlines.append('! --- reads from command line\n')
        listlines.append('  CALL read_from_cl\n\n')
        listlines.append('! --- mpi init communicator\n')

        if (self.type_solver == 'all') | (self.type_solver == 'spectral'):
            listlines.append('#if defined(FFTW)\n')
            listlines.append('  IF (fftw_with_mpi) THEN \n')
            listlines.append('    CALL mpi_minimal_init_fftw\n')
            listlines.append('  ELSE\n')
            listlines.append('#endif\n')
            listlines.append('    CALL mpi_minimal_init\n')
            listlines.append('#if defined(FFTW)\n')
            listlines.append('  ENDIF \n')
            listlines.append('#endif\n\n')
        else:
            listlines.append('  CALL mpi_minimal_init\n\n')

        listlines.append('  IF (rank .EQ. 0) THEN\n')
        listlines.append('    write(0,*) "___________________________________'\
                       + '______________________________"\n')
        listlines.append('    write(0,*) ""\n')
        listlines.append('    write(0,*) " PICSAR"\n')
        listlines.append('    write(0,*) "___________________________________'\
                       + '______________________________"\n')
        listlines.append('  ENDIF\n\n')

        listlines.append('! --- Check domain decomposition / Create Cartesian'\
                       + ' communicator / Allocate grid arrays\n')
        listlines.append('  CALL mpi_initialise\n\n')
        listlines.append('! --- allocates and inits particle distributions ' \
                       + '(on each subdomain)\n')
        listlines.append('  CALL initall\n\n')
        listlines.append('! --- Diagnostics\n')
        listlines.append('  CALL init_diags\n\n')

        listlines.append('  !----------------------------------------------\n')
        listlines.append('  ! THIS IS THE PIC ALGORITHM TIME LOOP\n')
        listlines.append('  !----------------------------------------------\n')
        listlines.append('  IF (rank .EQ. 0) startsim=MPI_WTIME()\n')
        listlines.append('  CALL step(nsteps)\n\n')

        listlines.append('  IF (rank .EQ. 0) endsim=MPI_WTIME()\n')
        listlines.append('  IF (rank .EQ. 0) WRITE(0,*)  "Total runtime on ",'\
                       + 'nproc," CPUS =",                   &\n')
        listlines.append('  endsim-startsim,"CPU AVERG TIME PER IT", '\
                       + '(endsim-startsim)/nsteps\n\n')

        listlines.append('  ! Time statistics for the different processes of '\
                       + 'the PIC step\n')
        listlines.append('  CALL time_statistics\n\n')

        listlines.append('  IF (rank .EQ. 0) THEN \n')
        listlines.append('	  INQUIRE(file="output_statistics.out", '\
                       + 'exist=exist)\n')
        listlines.append('  	IF (exist) THEN \n')
        listlines.append('  		OPEN (unit=12,file="output_statistics'\
                       + '.out", &\n')
        listlines.append('  		action="write",position="append", status='\
                       + '"old")\n')
        listlines.append('  	ELSE\n')
        listlines.append('  		OPEN (unit=12,file="output_statistics.out"'\
                       + ',  &\n')
        listlines.append('  		  action="write",status="new")\n')
        listlines.append('  	ENDIF \n')
        listlines.append('  	WRITE(str1,*) nx_global; WRITE(str2,*) '\
                       + 'ny_global\n')
        listlines.append('  	WRITE(str3,*) nz_global; WRITE(str4,*) nproc\n')
        listlines.append('  	! total simulation time\n')
        listlines.append('  	WRITE(str5,*) endsim-startsim\n')
        listlines.append('  	! Average time spent in different steps of '\
                       + 'the PIC loop\n')
        listlines.append("  	WRITE(str6,'(22(E12.5))') avetimes(1),"\
                       + 'avetimes(14),avetimes(2),avetimes(11),      &\n')
        listlines.append('  								avetimes(3),'\
                       + 'avetimes(4),avetimes(5),                  &\n')
        listlines.append('  								avetimes(6),'\
                       + 'avetimes(7),avetimes(21),                 &\n')
        listlines.append('  								avetimes(22),'\
                       + 'avetimes(23),avetimes(24), avetimes(25), &\n')
        listlines.append('  								avetimes(8),'\
                       + 'avetimes(10),                             &\n')
        listlines.append('  								avetimes(12),'\
                       + ' avetimes(13), avetimes(9),              &\n')
        listlines.append('  								avetimes(18), '\
                       + 'avetimes(19), avetimes(20)\n\n')

        listlines.append('  	! Total memory used in the case (in GB)\n')
        listlines.append("  	WRITE(str7,'(4(E12.5))') global_grid_mem/1e9,"\
                       + ' global_grid_tiles_mem/1e9,          &\n')
        listlines.append('  	global_part_tiles_mem/1e9\n\n')

        listlines.append('  	! All time are put in the file on a single'\
                       + ' line\n')
        listlines.append("  	WRITE(12, '(512A)')"+'  trim(adjustl(str1))//"'\
                       + ' "//trim(adjustl(str2))//" "//         &\n')
        listlines.append('  				  trim(adjustl(str3))//" '\
                       + '"//trim(adjustl(str4))//" "//                &\n')
        listlines.append(' 		    	      trim(adjustl(str5))//" '\
                       + '"//trim(adjustl(str6))//                 &\n')
        listlines.append('  				  " "//trim(adjustl(str7))\n')
        listlines.append('  	CLOSE(12)\n')
        listlines.append('  ENDIF \n\n')

        if (self.type_solver == 'all') | (self.type_solver == 'spectral'):
            listlines.append('#if defined(FFTW)\n')
            listlines.append('  IF(l_spectral) THEN\n')
            listlines.append('    IF(fftw_with_mpi) THEN\n')
            listlines.append('      CALL DFFTW_DESTROY_PLAN(plan_r2c_mpi)\n')
            listlines.append('      CALL DFFTW_DESTROY_PLAN(plan_c2r_mpi)\n')
            listlines.append('    ELSE\n')
            listlines.append('      CALL fast_fftw_destroy_plan_dft,'\
                           + '(plan_r2c)\n')
            listlines.append('      CALL fast_fftw_destroy_plan_dft,'\
                           + '(plan_c2r)\n')
            listlines.append('    ENDIF\n')
            listlines.append('  ENDIF\n')
            listlines.append('#endif\n')
            listlines.append('  CALL mpi_close\n\n')
            listlines.append('! Intel Design Forward project\n')
            listlines.append('#if defined(DFP)\n')
            listlines.append('   CALL DFP_FINAL_STOP\n')
            listlines.append('#endif\n\n')
            listlines.append('END PROGRAM main\n')

        else:
            listlines.append('END PROGRAM main\n')

        fnew.writelines(listlines)
        fnew.close()


###############################################################################
# Main

arglist=sys.argv
#try:
type_parser = str(arglist[1])
type_pusher = str(arglist[2])
type_depos  = str(arglist[3])
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

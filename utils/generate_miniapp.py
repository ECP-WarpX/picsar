import os, sys
import shutil
from datetime import datetime

# From the picsar root directory:
# python utils/generate_miniapp.py \
# --pusher boris --depos direct --solver fdtd --optimization off \
# --charge off --laser off --geom 3d --order 1 --diags off --errchk off
# Then:
# pushd PICSARlite
# Go into src/submain.F90 and remove "USE diagnostics" and "USE simple_io".
# make
# pushd examples/example_decks_fortran
# mpirun -np 8 ../../fortran_bin/picsar homogeneous_plasma_lite.pixr

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

class Interface( object ):
    def __init__( self, file ):

        f = open('./src/%s'%file,"r")
        listlines=f.readlines()

        istart=[]
        # Init iend with -1 to avoid code duplication, this element is removed
        # afterwards
        iend=[-1]
        Nlines=len(listlines)
        proc_mod=True
        for i in range(0,Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")

            # This is a subroutine block
            if ((curr_line.find("interface")>=0) \
                & (curr_line.find("!")==-1)):
                if i > iend[-1]:
                    istart.append(i)

                    while True: # get end of block
                        i=i+1

                        if (i>=Nlines):
                            sys.exit("ERROR: missing end interface block in "+file)
                        curr_line=listlines[i].lower()

                        if (curr_line.find("end interface")>=0):
                          break
                    iend.append(i)

        # Remove the -1 value of iend
        iend.remove(-1)

        # Create a file with only the interface in it and apply the subroutine
        # class to it.
        for j in range(len(istart)):
            ftemp = open('./src/%s_temp_%s'%(file, j), "w")
            listlines_temp = []
            for i in range(istart[j]+1, iend[j], 1):
                listlines_temp.append(listlines[i])

            ftemp.writelines(listlines_temp)
            ftemp.close()

        self.file         = file
        self.istart       = istart
        self.iend         = iend

        self.routine_istart = []
        self.routine_iend   = []
        self.routine_names  = []

        for j in range(len(istart)):
            r = Subroutines('%s_temp_%s'%(file, j))
            self.routine_istart.append(r.istart)
            self.routine_iend.append(r.iend)
            self.routine_names.append(r.names)
            os.system('rm ./src/%s_temp_%s'%(file, j))


class MiniAppParser( object ):

    def __init__( self, type_solver, type_pusher, type_depos,
                  flag_optimization, flag_charge, flag_laser, flag_geom,
                  flag_order, flag_diags, flag_errchk ):

        self.clean_folder()
        self.type_solver = type_solver
        self.type_pusher = type_pusher
        self.type_depos  = type_depos
        self.flag_optimization = flag_optimization
        self.flag_charge = flag_charge
        self.flag_laser = flag_laser
        self.flag_geom = flag_geom
        self.flag_order = flag_order
        self.flag_diags = flag_diags
        self.flag_errchk = flag_errchk

        self.include_geom_2d = False
        self.include_geom_3d = False
        if flag_geom == 'all':
            self.include_geom_2d = True
            self.include_geom_3d = True
        elif flag_geom == '2d':
            self.include_geom_2d = True
        elif flag_geom == '3d':
            self.include_geom_3d = True
        else:
            print('#########################################################' \
                  '#########')
            print('Wrong geometry argument there. The list of available '     \
                  'geometry is:')
            print('- all: both 2-dimensional and 3-dimensional')
            print('- 2d: 2-dimensional')
            print('- 3d: 3-dimensional')
            print('#########################################################' \
                              '#########')

        self.include_order_1 = False
        self.include_order_2 = False
        self.include_order_3 = False
        self.include_order_n = False
        if flag_order == 'all':
            self.include_order_1 = True
            self.include_order_2 = True
            self.include_order_3 = True
            self.include_order_n = True
        elif flag_order == '1':
            self.include_order_1 = True
        elif flag_order == '2':
            self.include_order_2 = True
        elif flag_order == '3':
            self.include_order_3 = True
        elif flag_order == 'n':
            self.include_order_n = True
        else:
            print('#########################################################' \
                  '#########')
            print('Wrong order there. The list of available order is:')
            print('- all: orders 1, 2, 3, n')
            print('- 1: order 1')
            print('- 2: order 2')
            print('- 3: order 3')
            print('- n: order n')
            print('#########################################################' \
                              '#########')

        self.include_solver_fdtd = False
        self.include_solver_spectral = False
        if type_solver == 'all':
            self.include_solver_fdtd = True
            self.include_solver_spectral = True
        elif type_solver == 'fdtd':
            self.include_solver_fdtd = True
        elif type_solver == 'spectral':
            self.include_solver_spectral = True
        else:
            print('#########################################################' \
                  '#########')
            print('Wrong solver argument there. The list of available '       \
                  'solver is:')
            print('- all: every routines and modules from picsar')
            print('- fdtd: every routines and modules related to the Finite ' \
                  'Domain')
            print('        Time Domain Maxwell solver in 3D.')
            print('- spectral: every routines and modules related to the '    \
                  'Pseudo-')
            print('        Spectral Analytical Time Domain Maxwell solver in '\
                  '3D.')
            print('#########################################################' \
                              '#########')

        self.include_pusher_boris = False
        self.include_pusher_vay = False
        if type_pusher == 'all':
            self.include_pusher_boris = True
            self.include_pusher_vay = True
        elif type_pusher == 'boris':
            self.include_pusher_boris = True
        elif type_pusher == 'vay':
            self.include_pusher_vay = True
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

        self.include_depos_direct = False
        self.include_depos_esirkepov = False
        if type_depos == 'all':
            self.include_depos_direct = True
            self.include_depos_esirkepov = True
        elif type_depos == 'direct':
            self.include_depos_direct = True
        elif type_depos == 'esirkepov':
            self.include_depos_esirkepov = True
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

        # Create the folder for the mini app
        os.makedirs('./PICSARlite')

        # To be completed
        # Solver
        generic_modules = [
                        "control_file",\
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
                        "mpi_routines", \
                        "python_pointers" ]

        external_modules = [
                         "p3dfft",\
                         "mpi",\
                         "omp_lib",\
                         "cufft",\
                         "iso_c_binding" ]


        diag_modules=["diagnostics", \
                      "simple_io" ]

        diag_routines=[
                        "calc_diags",\
                        "calc_field_div",\
                        "calc_field_divB",\
                        "init_diags",\
                        "init_temp_diags",\
                        "init_time_stat_output",\
                        "get_loc_kinetic_energy",\
                        "get_kinetic_energy",\
                        "get_loc_norm_2",\
                        "system",\
                        "get_norm_divErho",\
                        "output_routines",\
                        "output_temporal_diagnostics",\
                        "write_3d_field_array_to_file",\
                        "write_single_array_to_file",\
                        "write_particles_to_file",\
                        "get_particles_to_dump",\
                        "concatenate_particle_variable",\
                        "write_particle_variable",\
                        "output_time_statistics",\
                        "final_output_time_statistics",\
                        "pxrdepose_rho_on_grid",\
                        ]

        diag_routines_3D=[
                        "get_loc_field_energy",\
                        "get_field_energy",\
                        ]

        diag_routines_2D=[
                        "get_loc_field_energy_2d",\
                        "get_field_energy_2d",\
                        ]

        generic_routines=[
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
                        "add_group_of_particles_at_tile",\
                        "rm_particles_from_species_with_mask",\
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
                        "mpi_send",
                        "mpi_recv",
                        "mpi_isend",
                        "mpi_irecv",
                        "mpi_sendrecv",
                        "mpi_wait",
                        "mpi_waitall",
                        "mpi_type_free",
                        "mpi_send_init",
                        "mpi_start mpi_type_free",
                        "mpi_type_commit",
                        "mpi_type_vector",
                        "mpi_reduce",
                        "mpi_allgather",
                        "mpi_abort",
                        "mpi_bcast mpi_barrier",
                        "mpi_bcast",
                        "mpi_comm_free",
                        "mpi_comm_rank",
                        "mpi_comm_size",
                        "mpi_comm_group",
                        "mpi_comm_split",
                        "mpi_comm_free",
                        "mpi_comm_dup",
                        "mpi_cart_rank",
                        "mpi_comm_group",
                        "mpi_cart_shift",
                        "mpi_cart_create",
                        "mpi_dims_create",
                        "mpi_init_thread",
                        "mpi_initialized",
                        "mpi_type_create_subarray",
                        "mpi_type_commit",
                        "mpi_type_contiguous",
                        "mpi_type_vector",
                        "mpi_type_commit",
                        "mpi_type_create_struct",
                        "mpi_type_size",
                        "mpi_recv_init",
                        "mpi_start",
                        "mpi_barrier",
                        "mpi_allreduce",
                        "mpi_file_open",
                        "mpi_file_set_view",
                        "mpi_file_write_all",
                        "mpi_file_close",
                        "mpi_finalize",
                        "omp_set_nested",
                        "start_vtune_collection",
                        "start_sde_collection",
                        "stop_vtune_collection",
                        "stop_sde_collection",
                        "compute_send_recv_sizes_and_index_g2l_copies",
                        "compute_send_recv_sizes_and_index_l2g_copies",
                        "create_derived_types_groups",
                        "get_local_grid_mem",
                        "pxr_particle_sorting",
                        "particle_sorting_sub",
                        "get_local_number_of_particles_from_species",
                        "start_collection",
                        "stop_collection",
                        "random_number",
                        "get_loc_norm_diverho",
                        "getarg",
                        "get_command_argument",
                        "abort",
                        "mpi_cart_coords",
                        "allinea_start_sampling",
                        "allinea_stop_sampling",
                        "dfp_main_start",
                        "dfp_main_stop",
                        "dfp_final_start",
                        "estimate_total_memory_consumption",
                        "init_splitted_fields_random",
                        "step",
#                         "select_quantity",
#                         "pxr_convertindtoproc",
#                         "remap_em_2dfields",
#                         "remap_em_3dfields",
#                         "get_2dintersection",
#                         "get_3dintersection",
#                         "get_projected_load_on_x",
#                         "get_projected_load_on_y",
#                         "get_projected_load_on_z",
#                         "balance_in_dir",
#                         "get_proc_interval",
#                         "binary_search",
#                         "compute_effective_communication_setup",
                        ]

        generic_routines_3d = [
                        "add_particle_to_species",\
                        "add_particle_at_tile",\
                        "rm_particles_from_species",\
                        "rm_particle_at_tile",\
                        "pxr_particle_bin_sorting",
                            ]

        generic_routines_2d = [
                        "add_particle_to_species_2d",\
                        "add_particle_at_tile_2d",\
                        "rm_particles_from_species_2d",\
                        "rm_particle_at_tile_2d",\
                        "pxr_particle_bin_sorting_2d",\
                            ]

        generic_modules_solver = ["fields","field_boundary"]

        solver_routines = [
                        "field_bc",\
                        "summation_bcs",\
                        "summation_bcs_nonblocking",\
                        "summation_bcs_persistent_jx",\
                        "summation_bcs_persistent_jy",\
                        "summation_bcs_persistent_jz",\
                        "efield_bcs",\
                        "bfield_bcs",\
                        "current_bcs",\
                        "charge_bcs"]

        solver_routines_3d = [
                        "exchange_mpi_3d_grid_array_with_guards",\
                        "exchange_mpi_3d_grid_array_with_guards_nonblocking",\
                                     ]

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

        pusher_routines = [
                        "set_tile_split",\
                            ]

        pusher_routines_3d = [
                        "particle_bcs",\
                        "particle_bcs_tiles",\
                        "particle_bcs_tiles_openmp",\
#                        "particle_bsc_openmp_reordering",\ not used
                        "particle_bcs_mpi_blocking",\
                        "particle_bcs_mpi_non_blocking",\
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
                        "particle_bcs_tiles_and_mpi_3d",\
                            ]

        pusher_routines_2d = [
                        "particle_bcs_2d",\
                        "particle_bcs_tiles_2d",\
                        "particle_bcs_tiles_2d_openmp",\
                        "particle_bcs_mpi_non_blocking_2d",\
                        "field_gathering_plus_particle_pusher_sub_2d",\
                        "pxr_push2dxz",\
                                     ]

        boris_modules = []

        boris_routines_scalar_2d = [
                        "pxr_pushxz",\
                        "pxr_push2dxz"]

        boris_routines_scalar_3d = [
                        "pxr_boris_push_u_3d",\
                        "pxr_pushxyz",\
                        "pxr_epush_v",\
                        "pxr_bpush_v",\
                        "pxr_set_gamma"]

        boris_routines_vector_3d = [
                        "pxr_boris_push_rr_S09_u_3d",\
                        "pxr_boris_push_rr_B08_u_3d",\
                        "pxr_boris_push_rr_LL_u_3d",\
                        "pxr_boris_push_u_3d_block"]

        vay_pusher_modules = []

        vay_pusher_routines= ["pxr_ebcancelpush3d" ]

        laser_routines = [
                        "load_laser_species",\
                        "load_laser",\
                        "push_laser_particles",\
                        "laserp_pusher_gaussian",\
                        "laserp_pusher_hanning",\
                        "gaussian_profile",\
                        "hanning_profile",\
                         ]

        generic_modules_depos = [
                        "grid_tilemodule",\
                        "tile_params",\
                        "tiling",\
                        "particle_tilemodule"]

        depos_routines_charge_3d_openmp = [
                        "pxrdepose_rho_on_grid_sub_openmp_3d",\
                        "pxrdepose_rho_on_grid_sub_openmp_3d_n",\
                        ]

        depos_routines_charge_3d_openmp_scalar = [
                        "pxrdepose_rho_on_grid_sub_openmp_3d_scalar",\
                        ]

        depos_routines_charge_3d_openmp_vector = [
                        "pxrdepose_rho_on_grid_sub_openmp_3d_vecto",\
                        ]

        depos_routines_charge_2d_openmp = [
                        "pxrdepose_rho_on_grid_sub_openmp_2d",\
                        ]

        depos_scalar_routines_charge_2d = [
                        "pxr_depose_rho_n_2dxz",\
                        "pxr_depose_rhoold_n_2dxz",\
                                                  ]

        depos_scalar_routines_charge_3d_o1 = [
                        "depose_rho_scalar_1_1_1",\
                        ]

        depos_scalar_routines_charge_3d_o2 = [
                        "depose_rho_scalar_2_2_2",\
                        ]

        depos_scalar_routines_charge_3d_o3 = [
                        "depose_rho_scalar_3_3_3",\
                        ]

        depos_scalar_routines_charge_3d_on = [
                        "pxr_depose_rho_n",\
                        ]

        depos_generic_routines_current = [
                        "func_order",\
                        "curr_depo_sub",\
                        ]

        depos_generic_routines_current_3d = [
                        "depose_jxjyjz",\
                        "depose_jxjyjz_generic",\
                        "pxrdepose_currents_on_grid_jxjyjz",\
                        "pxrdepose_currents_on_grid_jxjyjz_classical_sub_seq",
                        ]

        depos_generic_routines_current_2d = [
                        "depose_jxjyjz_2d",\
                        "depose_jxjyjz_generic_2d",\
                        "depose_jxjyjz_2d depose_jxjyjz_generic_2d",\
                        "pxrdepose_currents_on_grid_jxjyjz_2d",\
                                                   ]
        depos_vector_routines_charge_3d_o1 = [
                        "depose_rho_vecSH_1_1_1",\
                        "depose_rho_vecNOY_1_1_1",\
                        "depose_rho_vecHV_1_1_1",\
                        "depose_rho_vecHVv2_1_1_1",\
                        ]

        depos_vector_routines_charge_3d_o2 = [
                        "depose_rho_vecHVv2_2_2_2",\
                        ]

        depos_vector_routines_charge_3d_o3 = [
                        "depose_rho_vecHVv2_3_3_3",\
                        "depose_rho_vecHVv3_3_3_3",\
                        "depose_rho_vecHVv4_3_3_3",\
                        ]

        esirkepov_routines_generic_openmp = [
                        "pxrdepose_currents_on_grid_jxjyjz_esirkepov2d_sub_openmp",\
                        "pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_openmp",\
                        ]

        direct_routines_generic_openmp = [
                        "pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp",\
                        ]

        direct_routines_generic_openmp_vect = [
                        "pxrdepose_currents_on_grid_jxjyjz_classical_sub_openmp_v2",\
                        ]

        depos_vector_routines_current = [
                        "pxrdepose_currents_on_grid_jxjyjz_sub_openmp",\
                        "curr_reduc_sub",\
                        ]

        esirkepov_routines_scalar_3d=[
                        "pxrdepose_currents_on_grid_jxjyjz_esirkepov_sub_seq",\
                        "depose_jxjyjz_esirkepov"]

        esirkepov_routines_scalar_2d=[
                        "depose_jxjyjz_esirkepov_2d"]

        esirkepov_routines_scalar_o1_3d=[
                        "depose_jxjyjz_esirkepov_1_1_1",\
                                     ]
                        
        esirkepov_routines_scalar_o1_2d=[
                        "pxr_depose_jxjyjz_esirkepov2d_1_1",\
                                     ]

        esirkepov_routines_scalar_o2_3d=[
                        "depose_jxjyjz_esirkepov_2_2_2",\
                                     ]

        esirkepov_routines_scalar_o2_2d=[
                        "pxr_depose_jxjyjz_esirkepov2d_2_2",\
                                     ]

        esirkepov_routines_scalar_o3_3d=[
                        "depose_jxjyjz_esirkepov_3_3_3",\
                                     ]

        esirkepov_routines_scalar_o3_2d=[
                        "pxr_depose_jxjyjz_esirkepov2d_3_3",\
                                     ]

        esirkepov_routines_scalar_on_3d=[
                        "pxr_depose_jxjyjz_esirkepov_n",\
                                     ]

        esirkepov_routines_scalar_on_2d=[
                        "pxr_depose_jxjyjz_esirkepov2d_n",\
                                     ]

        esirkepov_routines_vector_o3_2d=["pxr_depose_jxjyjz_esirkepov2d_vecHV_3_3"]

        direct_routines_scalar_o1 = ["depose_jxjyjz_scalar_1_1_1"]
        direct_routines_scalar_o2 = ["depose_jxjyjz_scalar_2_2_2"]
        direct_routines_scalar_o3 = ["depose_jxjyjz_scalar_3_3_3"]

        direct_routines_vector_o1 = ["depose_jxjyjz_vecHVv2_1_1_1",\
                        "depose_jxjyjz_vecHV_vnr_1_1_1",\
                        "current_reduction_1_1_1",\
                        ]

        direct_routines_vector_o2 = [
                        "depose_jxjyjz_vecHVv2_2_2_2",\
                        "depose_jxjyjz_vecHV_vnr_2_2_2",\
                        "current_reduction_2_2_2",\
                        ]

        direct_routines_vector_o3 = [
                        "depose_jxjyjz_vecHVv3_3_3_3",\
                        "depose_jxjyjz_vecHV_vnr_3_3_3",\
                        "current_reduction_3_3_3"
                         ]
                       
        gather_routines_scalar_2d = [
                        "geteb2dxz_energy_conserving",\
                        "geteb2dxz_energy_conserving_generic",\
                                            ]
        gather_routines_scalar_2d_o3 = [
                        "pxr_gete2dxz_energy_conserving_scalar_3_3",\
                        "pxr_getb2dxz_energy_conserving_scalar_3_3",\
                                            ]

        gather_routines_scalar_2d_on = [
                        "pxr_gete2dxz_n_energy_conserving",\
                        "pxr_getb2dxz_n_energy_conserving",\
                                            ]

        gather_routines_scalar_3d = [
                        "field_gathering",\
                        "field_gathering_sub",\
                        "geteb3d_energy_conserving",\
                        "geteb3d_energy_conserving_generic"
                         ]

        gather_routines_scalar_3d_o1 = [
                        "gete3d_energy_conserving_scalar_1_1_1",\
                        "getb3d_energy_conserving_scalar_1_1_1",\
                                            ]

        gather_routines_scalar_3d_o2 = [
                        "gete3d_energy_conserving_scalar_2_2_2",\
                        "getb3d_energy_conserving_scalar_2_2_2",\
                                            ]

        gather_routines_scalar_3d_o3 = [
                        "gete3d_energy_conserving_scalar_3_3_3",\
                        "getb3d_energy_conserving_scalar_3_3_3",\
                        "gete3d_energy_conserving_linear_3_3_3", \
                        "getb3d_energy_conserving_linear_3_3_3",\
                                            ]

        gather_routines_scalar_3d_on = [
                        "pxrgete3d_n_energy_conserving",\
                        "pxrgetb3d_n_energy_conserving",\
                        "pxr_getb3d_n_energy_conserving",\
                        "pxr_gete3d_n_energy_conserving",\
                                            ]

        gather_routines_vector_o1 = [
                        "pxr_gete2dxz_energy_conserving_vect_1_1",\
                        "pxr_getb2dxz_energy_conserving_vect_1_1",\
                        "gete3d_energy_conserving_vec_1_1_1",\
                        "getb3d_energy_conserving_vec_1_1_1",\
                        "geteb3d_energy_conserving_vecV1_1_1_1",\
                        "geteb3d_energy_conserving_vecV2_1_1_1",\
                        "geteb3d_energy_conserving_vecV3_1_1_1",\
                        "geteb3d_energy_conserving_vecV4_1_1_1",\
                        "geteb3d_energy_conserving_vec_1_1_1_v2",\
                        "geteb3d_energy_conserving_vec_1_1_1_sub",\
                        ]

        gather_routines_vector_o2 = [
                        "pxr_gete2dxz_energy_conserving_vect_2_2",\
                        "pxr_getb2dxz_energy_conserving_vect_2_2",\
                        "pxr_gete2dxz_energy_conserving_vect_2_2",\
                        "pxr_getb2dxz_energy_conserving_vect_2_2",\
                        "gete3d_energy_conserving_vec_2_2_2",\
                        "getb3d_energy_conserving_vec_2_2_2",\
                        "geteb3d_energy_conserving_vecV1_2_2_2",\
                        "geteb3d_energy_conserving_vecV2_2_2_2", \
                        "geteb3d_energy_conserving_vecV3_2_2_2", \
                        "geteb3d_energy_conserving_vecV4_2_2_2",\
                        ]

        gather_routines_vector_o3 = [
                        "pxr_gete2dxz_energy_conserving_vect_3_3",\
                        "pxr_getb2dxz_energy_conserving_vect_3_3",\
                        "pxr_geteb2dxz_energy_conserving_vect_3_3",\
                        "gete3d_energy_conserving_vec_3_3_3", \
                        "getb3d_energy_conserving_vec_3_3_3",\
                        "gete3d_energy_conserving_vec2_3_3_3",\
                        "getb3d_energy_conserving_vec2_3_3_3",\
                        "geteb3d_energy_conserving_vec_3_3_3",\
                        "geteb3d_energy_conserving_vecV2_3_3_3",\
                        "geteb3d_energy_conserving_vecV3_3_3_3",\
                        "geteb3d_energy_conserving_blockvec_3_3_3",\
                        "geteb3d_energy_conserving_blockvec2_3_3_3",\
                        ]

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
                        "init_plans_blocks",
                        "allocate_new_matrix_vector",
                        "multiply_mat_vector",
                        "multiply_unit_blocks",
                        "get2d_intersection_group_mpi",
                        "c_f_pointer",
                        "p3dfft_get_dims" ,
                        "p3dfft_setup",
                        "dfftw_init_threads",
                        "dfftw_plan_with_nthreads",
#                         "dfftw_plan_dft_1d",
#                         "dfftw_plan_dft_2d",
#                         "dfftw_plan_dft_3d",
#                         "dfftw_plan_dft_r2c_1d",
#                         "dfftw_plan_dft_c2r_1d",
                        "dfftw_execute_dft",
                        "dfftw_execute_dft_r2c",
                        "dfftw_execute_dft_c2r",
                        "dfftw_destroy_plan",
                        "fftw_mpi_init",
                        "p3dfft_ftran_r2c" ,
                        "fftw_mpi_execute_dft_r2c",
                        "p3dfft_ftran_c2r" ,
                        "fftw_mpi_execute_dft_c2r",
                        "p3dfft_btran_c2r",
                        ]

        spectral_routines_2d = [
                        "push_psaotd_ebfielfs_2d",
                        "fast_fftw_create_plan_r2c_2d_dft",
                        "fast_fftw_create_plan_c2r_2d_dft",
                        "fftw_mpi_local_size_2d p3dfft_setup" ,
                        "dfftw_plan_dft_r2c_2d",
                        "dfftw_plan_dft_c2r_2d",
                        "fftw_mpi_plan_dft_r2c_2d" ,
                        "fftw_mpi_plan_dft_c2r_2d" ,
                               ]

        spectral_routines_3d = [
                        "push_psaotd_ebfielfs_3d",
                        "fast_fftw_create_plan_r2c_3d_dft",
                        "fast_fftw_create_plan_c2r_3d_dft",
                        "fast_fftw3d_c2r_with_plan",
                        "fast_fftw3d_r2c_with_plan",
                        "fftw_mpi_local_size_3d" ,
                        "fftw_mpi_local_size_3d_transposed" ,
                        "dfftw_plan_dft_r2c_3d",
                        "dfftw_plan_dft_c2r_3d",
                        "fftw_mpi_init fftw_mpi_plan_dft_r2c_3d",
                        "fftw_mpi_plan_dft_c2r_3d",
                               ]

        fdtd_modules = []

        fdtd_routines = [
                        "push_bfield",
                        "push_efield",
                        "init_stencil_coefficients",
                        "FD_weights"
                        ]

        fdtd_routines_2d = [
                        "push_bfield_2d",
                        "push_efield_2d",
                        "pxrpush_em2d_evec_norder",
                        "pxrpush_em2d_evec",
                        "pxrpush_em2d_bvec_norder",
                        "pxrpush_em2d_bvec"
                           ]

        fdtd_routines_3d = [
                        "pxrpush_em3d_evec_norder",
                        "pxrpush_em3d_evec",
                        "pxrpush_em3d_bvec_norder",
                        "pxrpush_em3d_bvec"
                           ]

        self.list_available_modules = generic_modules                        \
                                    + generic_modules_solver                 \
                                    + generic_modules_pusher                 \
                                    + generic_modules_depos                  \
                                    + external_modules

        self.list_available_routines = generic_routines                      \
                                     + solver_routines                       \
                                     + pusher_routines                       \
                                     + depos_generic_routines_current

        if self.include_geom_2d:
            self.list_available_routines += generic_routines_2d            \
                                          + pusher_routines_2d             \
                                          + depos_generic_routines_current_2d

        if self.include_geom_3d:
            self.list_available_routines += generic_routines_3d     \
                                          + pusher_routines_3d      \
                                          + solver_routines_3d      \
                                          + depos_generic_routines_current_3d

        if self.flag_diags == 'on':
            self.list_available_modules += diag_modules
            self.list_available_routines += diag_routines
            if self.include_geom_2d:
                self.list_available_routines += diag_routines_2D
            if self.include_geom_3d:
                self.list_available_routines += diag_routines_3D

        if self.flag_laser == 'on':
            self.list_available_routines += laser_routines


        # add deposition routines

        if self.flag_charge == 'on':
            if self.include_geom_2d:
                self.list_available_routines += depos_routines_charge_2d_openmp
            if self.include_geom_3d:
                self.list_available_routines += depos_routines_charge_3d_openmp

        if self.include_depos_esirkepov:
            self.list_available_routines += esirkepov_routines_generic_openmp
        elif self.include_depos_direct:
            self.list_available_routines += direct_routines_generic_openmp

        if self.flag_optimization == 'on':
            if self.flag_charge == 'on':
                if self.include_geom_3d:
                    self.list_available_routines += depos_routines_charge_3d_openmp_vector
                    if self.include_order_1:
                        self.list_available_routines += depos_vector_routines_charge_3d_o1
                    if self.include_order_2:
                        self.list_available_routines += depos_vector_routines_charge_3d_o2
                    if self.include_order_3:
                        self.list_available_routines += depos_vector_routines_charge_3d_o3
                    if self.include_order_n:
                        self.list_available_routines += depos_vector_routines_charge_3d_on

            self.list_available_routines += depos_vector_routines_current
            if self.include_geom_2d:
                if self.include_depos_esirkepov:
                    if self.include_order_3:
                        self.list_available_routines += esirkepov_routines_vector_2d_o3
            elif self.include_geom_3d:
                if self.include_depos_direct:
                    if self.include_order_1:
                        self.list_available_routines += direct_routines_vector_o1
                    if self.include_order_2:
                        self.list_available_routines += direct_routines_vector_o2
                    if self.include_order_3:
                        self.list_available_routines += direct_routines_vector_o3
                    if self.include_order_n:
                        self.list_available_routines += direct_routines_generic_openmp_vect

        elif self.flag_optimization == 'off':
            if self.flag_charge == 'on':
                if self.include_geom_2d:
                    self.list_available_routines += depos_scalar_routines_charge_2d
                elif self.include_geom_3d:
                    self.list_available_routines += depos_routines_charge_3d_openmp_vector
                    if self.include_order_1:
                        self.list_available_routines += depos_scalar_routines_charge_3d_o1
                    if self.include_order_2:
                        self.list_available_routines += depos_scalar_routines_charge_3d_o2
                    if self.include_order_3:
                        self.list_available_routines += depos_scalar_routines_charge_3d_o3
                    if self.include_order_n:
                        self.list_available_routines += depos_scalar_routines_charge_3d_on

            if self.include_geom_2d:
                if self.include_depos_esirkepov:
                    self.list_available_routines += esirkepov_routines_scalar_2d
                    if self.include_order_1:
                        self.list_available_routines += esirkepov_routines_scalar_o1_2d
                    if self.include_order_2:
                        self.list_available_routines += esirkepov_routines_scalar_o2_2d
                    if self.include_order_3:
                        self.list_available_routines += esirkepov_routines_scalar_o3_2d
                    if self.include_order_n:
                        self.list_available_routines += esirkepov_routines_scalar_on_2d

            elif self.include_geom_3d:
                if self.include_depos_esirkepov:
                    self.list_available_routines += esirkepov_routines_scalar_3d
                    if self.include_order_1:
                        self.list_available_routines += esirkepov_routines_scalar_o1_3d
                    if self.include_order_2:
                        self.list_available_routines += esirkepov_routines_scalar_o2_3d
                    if self.include_order_3:
                        self.list_available_routines += esirkepov_routines_scalar_o3_3d
                    if self.include_order_n:
                        self.list_available_routines += esirkepov_routines_scalar_on_3d
                elif self.include_depos_direct:
                    if self.include_order_1:
                        self.list_available_routines += direct_routines_scalar_o1
                    if self.include_order_2:
                        self.list_available_routines += direct_routines_scalar_o2
                    if self.include_order_3:
                        self.list_available_routines += direct_routines_scalar_o3

        # Add gather subroutines
        if self.include_geom_2d:
            self.list_available_routines += gather_routines_scalar_2d
        if self.include_geom_3d:
            self.list_available_routines += gather_routines_scalar_3d
            if self.include_order_1:
                self.list_available_routines += gather_routines_scalar_3d_o1
            if self.include_order_2:
                self.list_available_routines += gather_routines_scalar_3d_o2
            if self.include_order_3:
                self.list_available_routines += gather_routines_scalar_3d_o3
            if self.include_order_n:
                self.list_available_routines += gather_routines_scalar_3d_on

        if self.flag_optimization == 'on':
            if self.include_order_1:
                self.list_available_routines += gather_routines_vector_o1
            if self.include_order_2:
                self.list_available_routines += gather_routines_vector_o2
            if self.include_order_3:
                self.list_available_routines += gather_routines_vector_o3

        if self.include_solver_spectral:
            self.list_available_modules  += spectral_modules
            self.list_available_routines += spectral_routines
            if self.include_geom_2d:
                self.list_available_routines += spectral_routines_2d
            if self.include_geom_3d:
                self.list_available_routines += spectral_routines_3d

        if self.include_solver_fdtd:
            self.list_available_modules  += fdtd_modules
            self.list_available_routines += fdtd_routines
            if self.include_geom_2d:
                self.list_available_routines += fdtd_routines_2d
            if self.include_geom_3d:
                self.list_available_routines += fdtd_routines_3d

        # Pusher
        if self.include_pusher_boris:
            self.list_available_modules  += boris_modules
            if self.include_geom_2d:
                self.list_available_routines += boris_routines_scalar_2d
            if self.include_geom_3d:
                self.list_available_routines += boris_routines_scalar_3d
                if self.flag_optimization == 'on':
                    self.list_available_routines += boris_routines_vector_3d

        if self.include_pusher_vay:
            self.list_available_modules  += vay_pusher_modules
            self.list_available_routines += vay_pusher_routines

        #LIST ALL .F90 or .F files in current directory
        self.listfiles = self.create_listfiles('./src')

        # Reconstruct PICSARlite
        print "Write routines in PICSARlite"
        for file in self.listfiles:
            self.write_available_routines(file)

        # Remove unavailable routines
        self.availablelistfiles = self.create_listfiles('./PICSARlite/src')
        #self.availablelistfiles = ['submain.F90']

        for file in self.availablelistfiles:
                self.comment_unavailable_routine(file)
                self.comment_unavailable_use(file)

        # Copy some extra needed folders
        print "Add extra files"
        self.copy_extra_files()
        print "Generate main"
        self.generate_main()
        print "Generate Makefile"
        self.generate_makefile()

    def clean_folder(self):

        # Delete all files
        if os.path.exists('./PICSARlite'):
            shutil.rmtree('./PICSARlite')

    def create_listfiles(self, folder):
        listfiles = []
        for root, subFolder, files in os.walk(folder):
            for file in files:
                if len(root[len(folder)+1:]) > 0:
                    listfiles.append('%s/%s'%(root[len(folder)+1:], file))
                else:
                    listfiles.append(file)
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

        # Lower the case of the modules
        lower_list_available_modules = []
        for module in self.list_available_modules:
            lower_list_available_modules.append(module.lower())
        # Lower the case of the routines
        lower_list_available_routines = []
        for routine in self.list_available_routines:
            lower_list_available_routines.append(routine.lower())

        for module in m.names:

            iname = m.names.index(module)

            # Loop on each module occurence, for each print the modules then
            # subroutines inside and finish by "end module"
            if (module in nb_modules) &                                  \
               (module in lower_list_available_modules):

                self.copy_files_from_picsar(m, iname)

            # Modules without subroutines
            else:
                if m.names[iname] in lower_list_available_modules:
                    self.copy_files_from_picsar(m, iname)

        # Add the subroutines out of modules
        for iname in range(len(r.names)):
            if r.proceduremod[iname] == '':
                if r.names[iname] in lower_list_available_routines:
                    self.copy_files_from_picsar(r, iname)

        if os.path.exists('./PICSARlite/src/%s'%file):
            self.manage_routine_interface(file)

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

    def manage_routine_interface(self, file):

        # Reopen the file to track the interface
        interface = Interface('../PICSARlite/src/%s'%file)
        if len(interface.istart) < 1:
            return

        lower_list_available_routines = []
        for routine in self.list_available_routines:
            lower_list_available_routines.append(routine.lower())

        f = open('PICSARlite/src/%s'%file, 'r')
        listlines = f.readlines()
        Nlines = len(listlines)
        compt = 0
        istart = [-1]; istart += interface.istart
        iend   = [-1]; iend += interface.iend

        fnew =  open('./PICSARlite/src/%s'%file, 'w')
        listlines_new = []
        for i in range(Nlines):
            if (compt == len(iend)-1) & (i >= iend[compt]):
                listlines_new.append(listlines[i])

            elif i >= iend[compt]:
                if i != istart[compt+1]:
                    listlines_new.append(listlines[i])

                # Begin the interface
                else:
                    listlines_new.append(listlines[i])
                    for index in range(len(interface.routine_istart[compt])):
                        if interface.routine_names[compt][index] in \
                                                lower_list_available_routines:
                            istart_routine = \
                interface.routine_istart[compt][index] + istart[compt+1]+1
                            iend_routine = \
                interface.routine_iend[compt][index] + istart[compt+1] + 1
                            for line in listlines[istart_routine:\
                                                  iend_routine+1]:
                                listlines_new.append(line)
                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

        fnew.writelines(listlines_new)
        fnew.close()

        # Comment calls if needed
        self.clean_call_interface(file, interface )

    def clean_call_interface(self, file, interface):

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

        # Lower the case of the routines
        lower_list_available_routines = []
        for routine in self.list_available_routines:
            lower_list_available_routines.append(routine.lower())

        interface_routine_names = []
        for i in range(len(interface.routine_names)):
            interface_routine_names += interface.routine_names[i]

        for i in range(0, Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")
            # We find a CALL
            if (("call" in curr_word_list) & (curr_line.find("!")==-1 )):
                # Get the name of the following routine
                indexcall = curr_word_list.index("call")
                routine_name = curr_word_list[indexcall+1].split('(')
                routine_name = routine_name[1].split(',')[0]
                routine_name = routine_name.split('\n')[0]

                # If the routine is in the list, everything fine, if not the
                # line and the block should be commented.
                if ((not routine_name in lower_list_available_routines) \
                 & (routine_name in interface_routine_names)):
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
            # If no more wrong subroutines, finish the file
            if (compt == len(iend)-1) & (i >= iend[compt]):
                listlines_new.append(listlines[i])

            elif i >= iend[compt]:
                # It is regular lines
                if i != istart[compt+1]:
                    listlines_new.append(listlines[i])

                # It is a call block
                elif self.flag_errchk == 'off':
                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

                else:
                    # Comment the lines
                    if iend[compt+1] == istart[compt+1]:
                        listlines_new.append('!'+listlines[i][1:])
                    else:
                        for iblock in range(iend[compt+1]-istart[compt+1]):
                            line = listlines[i+iblock][1:].split('\n')[0]
                            listlines_new.append('!'+line+'\n')

                    # Print error message
                    error_message = \
                formatting_line( nb_blanks[compt], \
                                       "WRITE(0,*) 'ABORT.',") + '\n' \
              + formatting_line( nb_blanks[compt]+11, \
                    "'The routine %s', "%(list_routine_name[compt])) + '\n'\
              + formatting_line( nb_blanks[compt]+11, \
                                "' cannot be used in this configuration.'",  \
                                add_ampersand=False)

                    listlines_new.append("\n" + error_message + "\n")
                    listlines_new.append(formatting_line( nb_blanks[compt], \
                                        "STOP", add_ampersand=False)+ "\n \n")
                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

        fnew.writelines(listlines_new)
        fnew.close()

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

        # Lower the case of the routines
        lower_list_available_routines = []
        for routine in self.list_available_routines:
            lower_list_available_routines.append(routine.lower())

        for i in range(0, Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")
            # We find a CALL
            if (("call" in curr_word_list) & (curr_line.find("!")==-1 )):
                # Get the name of the following routine
                indexcall = curr_word_list.index("call")
                routine_name = curr_word_list[indexcall+1].split('(')
                routine_name = routine_name[0].split('\n')[0]

                # If the routine is in the list, everything fine, if not the
                # line and the block should be commented.
                if not routine_name in lower_list_available_routines:
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
            # If no more wrong subroutines, finish the file
            if (compt == len(iend)-1) & (i >= iend[compt]):
                listlines_new.append(listlines[i])

            elif i >= iend[compt]:
                # It is regular lines
                if i != istart[compt+1]:
                    listlines_new.append(listlines[i])

                # It is a call block
                elif self.flag_errchk == 'off':
                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

                else:
                    # Comment the lines
                    if iend[compt+1] == istart[compt+1]:
                        listlines_new.append('!'+listlines[i][1:])
                    else:
                    # elif self.flag_errchk == 'on':
                        for iblock in range(iend[compt+1]-istart[compt+1]):
                            line = listlines[i+iblock][1:].split('\n')[0]
                            listlines_new.append('!'+line+'\n')

                    # Print error message
                    error_message = \
                formatting_line( nb_blanks[compt], \
                                       "WRITE(0,*) 'ABORT.',") + '\n' \
              + formatting_line( nb_blanks[compt]+11, \
                    "'The routine %s', "%(list_routine_name[compt])) + '\n'\
              + formatting_line( nb_blanks[compt]+11, \
                                "' cannot be used in this configuration.'",  \
                                add_ampersand=False)

                    listlines_new.append("\n" + error_message + "\n")
                    listlines_new.append(formatting_line( nb_blanks[compt], \
                                        "STOP", add_ampersand=False)+ "\n \n")
                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

        fnew.writelines(listlines_new)
        fnew.close()


    def comment_unavailable_use(self, file):
        """
            When all the files are written, check inside the files if some
            unavailable module are unvailable but used. It may introduce
            a crash at compilation time.

            To avoid it, this function removes the "USE" for these modules.
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
        list_module_name = []
        nb_blanks = []

        # Lower the case of the routines
        lower_list_available_modules = []
        for module in self.list_available_modules:
            lower_list_available_modules.append(module.lower())

        for i in range(0, Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")
            # We find a CALL
            if (("use" in curr_word_list) & (curr_line.find("!")==-1 )):
                # Get the name of the following routine
                indexcall = curr_word_list.index("use")
                module_name = curr_word_list[indexcall+1].split(' ')
                module_name = module_name[0].split(',')
                module_name = module_name[0].split('\n')[0]

                # If the module is in the list, everything fine, if not the
                # line should be commented.
                if not module_name in lower_list_available_modules:
                    list_module_name.append(module_name)
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
            # If no more wrong subroutines, finish the file
            if (compt == len(iend)-1) & (i >= iend[compt]):
                listlines_new.append(listlines[i])

            elif i >= iend[compt]:
                # It is regular lines
                if i != istart[compt+1]:
                    listlines_new.append(listlines[i])

                # It is a call block
                elif self.flag_errchk == 'off':
                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

                else:
                    # elif self.flag_errchk == 'on':
                    # Comment the lines
                    if iend[compt+1] == istart[compt+1]:
                        listlines_new.append('!'+listlines[i])

                    else:
                        for iblock in range(iend[compt+1]-istart[compt+1]):
                            line = listlines[i+iblock].split('\n')[0]
                            listlines_new.append('!'+line+'\n')

                    if compt == len(iend)-1:
                        compt = 0
                    else:
                        compt += 1

        fnew.writelines(listlines_new)
        fnew.close()


    def copy_extra_files(self):
        os.system('cp -r ./utils ./PICSARlite/utils')
        os.system('cp -r ./examples ./PICSARlite/examples')
        os.system('cp ./configure ./PICSARlite/configure')
        os.system('cp ./Makefile_Forthon.in ./PICSARlite/Makefile_Forthon.in')

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

        # Low the case of the routines
        lower_list_available_modules = []
        for module in self.list_available_modules:
            lower_list_available_modules.append(module.lower())

        for module in listmodulesmain:
            if module in lower_list_available_modules:
                listlines.append('  USE %s \n'%module)

        if "mem_status" in lower_list_available_modules:
            listlines.append('  USE mem_status, ONLY : global_grid_mem, '     \
                            +'global_grid_tiles_mem, global_part_tiles_mem\n\n')

        if self.include_solver_spectral:
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

        if self.include_solver_spectral:
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
        if self.flag_diags == 'on':
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

        if self.include_solver_spectral:
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

        fnew.writelines(listlines)
        fnew.close()

    def generate_makefile( self ):
        file = 'Makefile'
        fold = open('./%s'%file, "r")
        listlines=fold.readlines()
        Nlines = len(listlines)
        listlinesnew = []

        # Low the case of the available files and remove the .F90 if needed
        lower_availablelistfiles = []
        for availablefile in self.availablelistfiles:
            if availablefile[-4:] == '.F90':
                availablefile = availablefile[:-4]
            lower_availablelistfiles.append(availablefile.lower())

        for i in range(0, Nlines):
            curr_line=listlines[i].lower()
            curr_word_list=curr_line.split(" ")

            if  "mode=" in curr_word_list:
                if self.type_solver == 'spectral':
                    line = listlines[i][:5] + "prod_spectral"
                    listlinesnew.append(line)
                else:
                    line = listlines[i][:5] + "prod"
                    listlinesnew.append(line)

            # Copy the first part of the old Makefile
            if curr_line.find("lib: echo createdir build_lib")==-1:
                listlinesnew.append(listlines[i])

            #Break before the build_lib
            else:
                listlinesnew.append(listlines[i])
                break



        # Append the existing files to the build
        flag_end = False
        while True:
            i += 1
            curr_line=listlines[i].lower()
            curr_word_list= curr_line.split(" ")
            # If clean is in the current line, it means that the build is done
            if i >= Nlines:
                break
                #sys.exit("ERROR: missing end clean")
            else:
                if ("clean" in curr_line):
                    break

            while True:
                if flag_end:
                    break

                if i >= Nlines:
                    sys.exit("ERROR: missing end build block ")

                else:
                    if "endif" in curr_line:
                        listlinesnew.append(listlines[i])
                        flag_end = True
                        break
                    elif ("ifeq" in curr_line) | ("else" in curr_line):
                        listlinesnew.append(listlines[i])
                        break
                    else:
                        curr_line=listlines[i].lower()
                        curr_word_list= curr_line.split(" ")
                        if (len(curr_word_list) ==2) & \
                            (curr_line.find('$(srcdir)') >0):
                            if "build" in curr_line:
                                if "build_lib" in curr_line:
                                    curr_word_list= \
                                            curr_line.split("build_lib:")[1]
                                else:
                                    curr_word_list= \
                                                curr_line.split("build:")
                                    curr_word_list =curr_word_list[1]
                                curr_word_list=curr_word_list.split("\t")
                                routine = curr_word_list[0].split(' ')[0]
                                routine = routine[10:-2]


                            else:
                                curr_word_list=curr_line.split("\t")
                                routine = curr_word_list[1].split(' ')[0]
                                routine = routine[10:-2]

                            # Parse the ith line of build_lib, check if the file exists
                            # If not skip the line
                            if routine in lower_availablelistfiles:
                                listlinesnew.append(listlines[i])

                        else:
                            listlinesnew.append(listlines[i])
                        i += 1

        # End the file until the acceptance tests
        while True:

            if i >= Nlines:
                sys.exit("ERROR: missing acceptance tests ")
            else:
                curr_line=listlines[i].lower()

                # --- removes 'cleantest'
                if ("clean:" in curr_line) and ("clean_test" in curr_line):
                    listlines[i] = listlines[i].replace("clean_test","")

                curr_word_list= curr_line.split(" ")
                if "# __________________________________________" in curr_line:
                    break
                else:
                    listlinesnew.append(listlines[i])
                    i += 1

        fnew = open('./PICSARlite/%s'%file, "w")
        fnew.writelines(listlinesnew)
        fnew.close()


###############################################################################
import argparse

parser = argparse.ArgumentParser(description='Parse inputs to generate_miniapp.',
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--solver', dest='type_solver', default='fdtd',
                    choices = ['all','fdtd','spectral'],
                    help = 'Maxwell solvers among: \n'\
                    +'- fdtd:     3D-Finite Domain Time Domain Maxwell '\
                    + 'solver\n'\
                    + '- spectral: 3D-Pseudo-Spectral Analytical Time Domain' \
                    + ' Maxwell solver. \n'       \
                    + '- all:      both solvers \n'    \
                    + 'Default all.\n\n')

parser.add_argument('--pusher', dest='type_pusher', default='boris',
                    choices = ['all','boris','vay'],
                    help='Particle pushers among: \n'\
                    + '- boris: Boris pusher \n- vay:   Vay pusher\n' \
                    + '- all:   both particle pushers \n'    \
                    + 'Default all.\n\n')

parser.add_argument('--depos',  dest='type_depos',  default='direct',
                    choices = ['all','direct','esirkepov'],
                    help='Particle depositions among: \n'\
                    + '- direct:    Direct deposition \n'\
                    + '- esirkepov: Esirkepov deposition\n' \
                    + '- all:       both particle depositions\n'   \
                    + 'Default all.\n\n')

parser.add_argument('--optimization',  dest='flag_optimization', default='off',
                    choices = ['on','off'],
                    help='flag optimization \n' \
                    +'Default on.\n\n')

parser.add_argument('--charge',  dest='flag_charge', default='off',
                    choices = ['on','off'],
                    help='flag charge \n' \
                    +'Default on.\n\n')

parser.add_argument('--laser',  dest='flag_laser', default='off',
                    choices = ['on','off'],
                    help='flag include laser pusher\n' \
                    +'Default on.\n\n')

parser.add_argument('--geom',  dest='flag_geom', default='3d',
                    choices = ['2d','3d','all'],
                    help='Geometry among: \n'\
                    + '- 2d:    2-dimensional\n'\
                    + '- 3d:    3-dimensional\n' \
                    + '- all:   both 2-dimensional and 3-dimensional\n'   \
                    + 'Default all.\n\n')

parser.add_argument('--order',  dest='flag_order', default='1',
                    choices = ['1','2','3','n','all'],
                    help='Field gathering order among: \n'\
                    + '- 1:    order 1 \n'\
                    + '- 2:    order 2 \n'\
                    + '- 3:    order 3 \n'\
                    + '- n:    order n \n'\
                    + '- all:  all orders 1, 2, 3, n\n'   \
                    + 'Default all.\n\n')

parser.add_argument('--diags',  dest='flag_diags', default='off',
                    choices = ['on','off'],
                    help='flag diagnostics \n' \
                    +'Default on.\n\n')

parser.add_argument('--errchk',  dest='flag_errchk', default='off',
                    choices = ['on','off'],
                    help='flag error, print warning messages at runtime when' \
                    +'the code is trying to use unavailable routines. \n' \
                    +'Note that this flag must be set to off when the flag' \
                    +'--diags is off also. \n'
                    +'Default off.\n')

args=parser.parse_args()

type_solver = args.type_solver
type_pusher = args.type_pusher
type_depos  = args.type_depos
flag_charge = args.flag_charge
flag_laser  = args.flag_laser
flag_geom   = args.flag_geom
flag_order  = args.flag_order
flag_diags  = args.flag_diags
flag_errchk = args.flag_errchk
flag_optimization = args.flag_optimization

miniapp = MiniAppParser(type_solver, type_pusher, type_depos, flag_optimization,
       flag_charge, flag_laser, flag_geom, flag_order, flag_diags, flag_errchk)

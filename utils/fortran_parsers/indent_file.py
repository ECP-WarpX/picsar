"""
 _______________________________________________________________________________

  *** Copyright Notice ***

 "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
 2016, The Regents of the University of California, through Lawrence Berkeley
 National Laboratory (subject to receipt of any required approvals from the
 U.S. Dept. of Energy). All rights reserved.

 If you have questions about your rights to use or distribute this software,
 please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

 NOTICE.
 This Software was developed under funding from the U.S. Department of Energy
 and the U.S. Government consequently retains certain rights. As such, the U.S.
 Government has been granted for itself and others acting on its behalf a
 paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 reproduce, distribute copies to the public, prepare derivative works, and
 perform publicly and display publicly, and to permit other to do so.

 PARSER FOR CORRECTING INDENTATION IN FORTRAN FILES 
 H. VINCENTI - June, 2017


 The function correct_indentation below reads Fortran 90 files and correct its indentation 
 according to guidelines detailed in CONTRIBUTING.md 

 Input arguments: filename

 Outputs:
 1. New fortran file guidelines compliant

 Developers:
 Henri Vincenti

 Date:
 Creation: June 1 2017
 _______________________________________________________________________________
"""

import re


def correct_indentation(input_file,output_file,indent_block):
	# ----- Parse input file 
	# -- Open input file 
	f_input=open(input_file,"r")
	# -- Get input file lines as a list of strings
	listlines_input=f_input.readlines()
	f_input.close()
	nlines=len(listlines_input)
	# -- Init output list of lines 
	listlines_output=[]
	# -- Block pattern 
	pattern_block_start="^\s*(DO|SUBROUTINE|MODULE|PROGRAM|IF.*THEN|FUNCTION|SELECT|INTERFACE|TYPE ).*"
	pattern_block_end="^\s*(END |ENDDO|ENDIF|END$|END\!)"
	pattern_block_intermediate="^\s*(CASE|ELSE)"
	pattern_no_trailing_spaces="^\s*(.*)(?<=\S)\s*$"
	pattern_comp_directive="^\s*#"
	# -- Sequentially go through input list of lines to identify blocks to indent 
	# - Some init first 
	curr_indent=0
	for iline,line in enumerate(listlines_input): 
		# - Process current line without spaces
		match_line_no_space=re.search(pattern_no_trailing_spaces,line)
		# - find if current line matches beginning of a new Fortran block to indent
		match_block_start=re.search(pattern_block_start,line,re.IGNORECASE)
		match_block_intermediate  =re.search(pattern_block_intermediate,line,re.IGNORECASE)
		match_block_end  =re.search(pattern_block_end,line,re.IGNORECASE)
		match_block_comp_dir=re.search(pattern_comp_directive,line,re.IGNORECASE)
		if match_line_no_space is None: # Line contains only spaces; replace with empty line
			listlines_output.append("\n")
		elif match_block_start is not None: # Current line starts a Fortran block		
			nwspaces=' ' * curr_indent
			listlines_output.append(nwspaces+match_line_no_space.group(1)+"\n")
			curr_indent+=indent_block
		elif match_block_comp_dir is not None: # Current line is compiler dir no indent
			listlines_output.append(match_line_no_space.group(1)+"\n")	
		elif match_block_end is not None: # Current line ends a Fortran block 
			curr_indent-=indent_block
			nwspaces=' ' * curr_indent
			listlines_output.append(nwspaces+match_line_no_space.group(1)+"\n")
		elif match_block_intermediate is not None:  # Current line is intermediate block (else
													# case, etc. 
			curr_indent-=indent_block
			nwspaces=' ' * curr_indent
			listlines_output.append(nwspaces+match_line_no_space.group(1)+"\n")
			curr_indent+=indent_block
		else: # Current line is regular line (not a Fortran block) 
			nwspaces=' ' * curr_indent
			listlines_output.append(nwspaces+match_line_no_space.group(1)+"\n")
		#print iline, curr_indent

	# -- Write indented output file 
	f_output=open(output_file,"w")
	f_output.writelines(listlines_output)
	f_output.close()


# Size of indent (in number of spaces) for Fortran blocks
# A Fortran block starts either with a SUBROUTINE,MODULE,PROGRAM, FUNCTION etc. key word
indent_block=2
#LIST ALL .F90 or .F files in current directory
listfiles=["modules/modules.F90", \
           "housekeeping/sorting.F90", \
           "field_solvers/Maxwell/maxwell_solver_manager.F90", \
           "field_solvers/Maxwell/yee_solver/yee.F90", \
           "field_solvers/Maxwell/karkkainen_solver/karkkainen.F90", \
           "field_solvers/Maxwell/GPSTD_solver/GPSTD.F90", \
           "field_solvers/Maxwell/GPSTD_solver/fourier_psaotd.F90", \
           "field_solvers/Maxwell/GPSTD_solver/fastfft.F90", \
           "field_solvers/Maxwell/GPSTD_solver/init_kspace_3D.F90", \
           "particle_pushers/kin_energy.F90", \
           "parallelization/tiling/tiling.F90", \
           "particle_pushers/boris_pusher/boris_2d.F90", \
           "particle_pushers/boris_pusher/boris_3d.F90", \
           "particle_pushers/vay_pusher/vay_3d.F90", \
           "particle_pushers/particle_pusher_manager_2d.F90", \
           "particle_pushers/particle_pusher_manager_3d.F90", \
           "particle_pushers/laser_pusher_manager_3d.F90", \
           "particle_deposition/current_deposition/direct/direct_current_deposition_3d.F90", \
           "particle_deposition/current_deposition/esirkepov/esirkepov_2d.F90", \
           "particle_deposition/current_deposition/esirkepov/esirkepov_3d.F90", \
           "particle_deposition/current_deposition/current_deposition_manager_2d.F90", \
           "particle_deposition/current_deposition/current_deposition_manager_3d.F90", \
           "field_gathering/field_gathering_manager_2d.F90", \
           "field_gathering/field_gathering_manager_3d.F90", \
           "field_gathering/energy_conserving/field_gathering_on_3d.F90",\
           "field_gathering/energy_conserving/field_gathering_o1_3d.F90",\
           "field_gathering/energy_conserving/field_gathering_o2_3d.F90",\
           "field_gathering/energy_conserving/field_gathering_o3_3d.F90",\
           "field_gathering/energy_conserving/field_gathering_on_2d.F90",\
           "field_gathering/energy_conserving/field_gathering_o1_2d.F90",\
           "field_gathering/energy_conserving/field_gathering_o2_2d.F90",\
           "field_gathering/energy_conserving/field_gathering_o3_2d.F90",\
           "parallelization/mpi/mpi_derived_types.F90",\
           "boundary_conditions/field_boundaries.F90", \
           "boundary_conditions/particle_boundaries.F90", \
           "ios/simple_io.F90", \
           "particle_deposition/charge_deposition/charge_deposition_2d.F90", \
           "particle_deposition/charge_deposition/charge_deposition_3d.F90", \
           "particle_deposition/charge_deposition/charge_deposition_manager.F90", \
           "diags/diags.F90", \
           "submain.F90", \
           "parallelization/mpi/mpi_routines.F90",\
           "initialization/control_file.F90", \
           "housekeeping/load_balancing.F90"]

for file in listfiles: 
	correct_indentation(file,file,indent_block)

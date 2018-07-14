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


 The function justify_file below reads Fortran 90 files and correct their justification 
 according to guidelines detailed in CONTRIBUTING.md 

 Input arguments: filename

 Outputs:
 1. New fortran file guidelines compliant

 Developers:
 Henri Vincenti

 Date:
 Creation: June 14 2017
 _______________________________________________________________________________
"""

import re

#Maximum columns of characters in a code line 
max_chars=85
# Ampersand column location 
ampersand_col=max_chars+1

"""
 --- This function pre-processes all lines of the file 
 --- to correctly insert spaces between commas. 
 --- The guideline imposes one space only after a comma and no space before the comma. 
 --- Input parameters: 
     - string input_filename containing the name of file to be pre-processed
 --- Output parameters: 
     - No outputs 
"""
def parse_commas(input_filename):
    fin=open(input_filename,'r')
    list_of_lines=fin.readlines()
    fin.close()
    fout=open(input_filename,'w')
    for line in list_of_lines: 
       output_line=re.sub("\s*,\s*",", ",line) 
       fout.write(output_line)
    fout.close()

"""
 --- This function truncates recusrively a line that exceeds the maximum line length 
 --- max_chars. 
 --- Input parameters: 
     - input_line: initial line to be truncated 
     - max_chars : maximum line length (in number of characters)
     - baseindent: base indent of the current input_line in the main program.
       This baseindent is useful to add the required number of leading spaces 
       when truncating the line into continuation line. 
     - directive : indicate if the current_line is prefixed by a directive in the main 
       code. If yes, a directive key-word has to be inserted after each continuation line.  
     - ampersand_col: position of the ampersand char & (in number of characters) 
       for line continuation 
     - split_char: string containing splitting character used to truncate source lines 
 --- Output parameters: 
     - python string containing a new line properly truncated
"""
def operate_truncation(input_line,max_chars,baseindent,directive,ampersand_col,split_char): 
   current_line=input_line
   lbindent=len(baseindent)
   l0=len(input_line)+lbindent
   # --- Line is too long - split it recursively 
   # Beware!!!!: A line has to be split at special characters and not in the middle of 
   # variables names. String names can be split. 
   if (l0>max_chars): 
      # Split current line into chunks separated by commas
      list_chunks=current_line.split(split_char)
      line_length=lbindent
      line_split_1=''
      isplit=0
      # Check for last chunk which is before ampersand comment 
      for i in range(len(list_chunks)): 
         line_length+=len(list_chunks[i])+len(split_char)
         if (line_length>ampersand_col): 
            break
         else: 
            line_split_1=line_split_1+list_chunks[i]+split_char
         isplit+=1
      if isplit==0: 
         print("Error, cannot split current line :",current_line)
         return baseindent+current_line+"\n"
      else:
         line_split_2= split_char.join(list_chunks[isplit:])
         wampersand=' ' * (ampersand_col-len(line_split_1)-lbindent)
         line_split_1=line_split_1+wampersand+"&"+"\n"
         return baseindent+line_split_1+ \
         operate_truncation(directive+line_split_2,max_chars,baseindent,directive, \
         ampersand_col, split_char)
   else: 
      return baseindent+current_line

"""
 --- This function pre-processes a line 
 --- by splitting it in code part and comment part 
 --- If current line is a continuation line, this function is called 
 --- recursively until all continuation lines are concatenated 
 --- This will then allow to simply truncate the line properly according to 
 --- our pre-defined guidelines. 
 --- Input paramaters: 
     - fh: file handler 
     - curr_line: string containing current line to be pre-processed 
 --- Output paramaters: 
     - python list of three strings ([code part,comment part,directive part] of 
       the current line) 
"""
def pre_process_line(fh,curr_line):
   # -- Test if current line is a pre-processor directive 
   # -- nothing to do in this case 
   pattern_preproc="^#"
   if re.search(pattern_preproc,curr_line): 
      return [curr_line[:-1],'','']
   # -- Test if current line is a compiler directive
   pattern_dir="^\s*(!\$[A-Z,a-z]* )"
   match_dir=re.search(pattern_dir,curr_line)
   if (match_dir is None): # This is a regular code line 
      directive=''
      split_line_omp=curr_line
   else: # This is a compiler directive
      directive=match_dir.group(1)
      split_line_omp=curr_line.split(directive)[1]
   # -- remove comments from the code part of the line
   # - Split current line (without new line char) 
   # - to separate code part from comments 
   split_line_com = split_line_omp[:-1].split('!',1)
   code_part    = split_line_com[0]
   if len(split_line_com)>1: 
      comment_part = '!'+split_line_com[1]
   else: 
      comment_part = ''
   # -- Test if code part has an ampersand continuation line
   # -- If it does, call recursively pre_process_line
   split_code_amper=code_part.split('&')
   code_part=re.sub("\s+$","",split_code_amper[0])
   code_part=re.sub("^\s+","",code_part)
   if (len(split_code_amper)>1): 
      curr_line = fh.readline()
      [code_part_sub,comment_part_sub,directive]=pre_process_line(fh,curr_line)
      return [code_part+" "+code_part_sub, comment_part+" "+comment_part_sub,directive]
   else: 
      return [code_part,comment_part,directive]

"""
 --- This function gets the current indent (leading spaces) of a line 
      of code 
 --- Input parameters: 
     - string curr_line containing the input line of code
 --- Output parameters: 
     - Leading spaces of curr_line 
"""
def get_current_indent(curr_line): 
   pattern_indent='^( *)'
   match_pat=re.search(pattern_indent,curr_line)
   return match_pat.group(0) 


def justify_file(input_file,output_file): 
	# Current code is for testing only 
	parse_commas(input_file)
	fin=open(input_file,"r")
	list_of_output_lines=[]
	curr_line=fin.readline()
	split_char=" "
	while (curr_line != ""):
	   # Get indentation lever of current line 
	   curr_indent=get_current_indent(curr_line)
	   # Pre-process current line 
	   [codepart,commentpart,directive]=pre_process_line(fin,curr_line)
	   #print(directive,codepart,commentpart)
	   # Assemble pre-processed line 
	   preproc_line= directive+codepart#+commentpart+"\n"
	   # Check pre-processed line length and truncate recursively if necessary
	   output_line = operate_truncation(preproc_line,max_chars,curr_indent,directive, \
       ampersand_col, split_char)
	   # Add potential existing comments at EOL
	   output_line = output_line+commentpart+"\n"
	   #fout.write(output_line)
	   list_of_output_lines.append(output_line)
	   curr_line=fin.readline()
	fin.close()
	fout=open(output_file,"w")
	fout.write("".join(list_of_output_lines))
	fout.close()

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
    justify_file(file,file)


import re

# Size of indent (in number of spaces) for Fortran blocks
# A Fortran block starts either with a SUBROUTINE,MODULE,PROGRAM, FUNCTION etc. key word
indent_block=2
# Input file name to parse 
input_file="/Users/henrivincenti/picsar/src/modules/modules.F90"
#Output file name 
output_file="/Users/henrivincenti/picsar/src/modules/test_modules.F90"

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
pattern_block_end="^\s*(END)"
pattern_block_intermediate="^\s*(CASE|ELSE)"
pattern_no_trailing_spaces="^\s*(.*)"
pattern_comp_directive="^#"
# -- Sequentially go through input list of lines to identify blocks to indent 
# - Some init first 
curr_indent=0
for iline,line in enumerate(listlines_input): 
    # - Process current line without spaces
	match_line_no_space=re.search(pattern_no_trailing_spaces,line)
	# - find if current line matches beginning of a new Fortran block to indent
	match_block_start=re.search(pattern_block_start,line)
	match_block_intermediate  =re.search(pattern_block_intermediate,line)
	match_block_end  =re.search(pattern_block_end,line)
	match_block_comp_dir=re.search(pattern_comp_directive,line)
	if match_block_start is not None: # Current line starts a Fortran block		
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

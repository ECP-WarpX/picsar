"""
This script replaces the syntax

USE some_module

with the syntax

USE some_module, ONLY: some_variable1, some_variable2

The variables are automatically detected.

Usage
-----
python use_module_only.py
"""
import glob
import re, sys, os

# Modules that are not defined within picsar
known_external_modules = ['mpi', 'omp_lib', 'p3dfft', 'itt_sde_fortran']


def reconstruct_lines(lines):
    """Reconstruction full lines from Fortran broken lines using &"""
    new_lines = []
    current_line = ''
    for line in lines:
        new_line = line.rstrip(' ')
        if new_line.endswith('&'):
            new_line = new_line.rstrip('& ')
            current_line += new_line
        else:
            current_line += new_line
            new_lines.append(current_line)
            current_line = ''
    return new_lines


def get_module_variables(listlines, dict_modules):
    """Modifies `dict_modules` to with keys the module name
    and value a list of existing module variables"""
    current_key = None
    inside_type = False
    for line in listlines:
        # Detect beginning of type
        if re.match('\s*type', line, re.IGNORECASE):
            if re.search('::', line): #Single-line type declaration
                continue
            else:
                inside_type = True
        # Detect module beginning
        m = re.match("^\s*module\s*(\w+)", line, flags=re.IGNORECASE)
        if m:
            if current_key is not None:
                raise ValueError('Still parsing %s' %line )
            current_key = m.group(1).lower()
            dict_modules[current_key] = []
        # Detect module end
        m = re.match("^\s*end module", line, flags=re.IGNORECASE)
        if m:
            current_key = None
        # Detect variable names
        if (not inside_type) and (current_key is not None):
            m = re.match('.*::(.*)', line)
            if m:
                for variable in re.findall('(\w+)', m.group(1)):
                    if not re.match('[\d_]', variable): # Exclude pure numbers
                        dict_modules[current_key].append( variable.lower() )
        # Detect end of type
        if re.match('\s*end type', line, re.IGNORECASE):
            inside_type = False

def get_subroutines(listlines):
    """Return dictionary with text of the subroutines"""
    dict_subroutines = {}
    current_key = None
    for line in listlines:
        # Detect beginning of subroutine
        m = re.match('\s*subroutine (\w+)', line, re.IGNORECASE)
        if m:
            current_key = m.group(1).lower()
            dict_subroutines[current_key] = ''  # Initialize empty string
        # Add line to text
        if current_key is not None:
            dict_subroutines[current_key] += line
        # Detect end of subroutine
        m = re.match('\s*end subroutine (\w+)', line, re.IGNORECASE)
        if m:
            current_key = None
    return( dict_subroutines )

def get_sub_module( dict_subs, dict_modules ):
    dict_subs_modules = {}

    # Loop over the subroutines
    for name, text in dict_subs.items():
        module_list = re.findall( '\n\s*use (\w+)', text, re.IGNORECASE )
        dict_subs_modules[name] = {}
        # Loop over the modules used in this subroutine
        for module in module_list:
            module = module.lower()
            if module not in dict_modules.keys():
                if module not in known_external_modules:
                    print('WARNING: Module %s seems to be external.' %module)
                continue
            dict_subs_modules[name][module] = []
            # Loop over the variables defined by this module
            for variable in dict_modules[module]:
                # Check whether this variable is being used
                if re.search('[\W_]%s[\W\n]'%variable, text, re.IGNORECASE):
                    dict_subs_modules[name][module].append(variable)

    return dict_subs_modules

def rewrite_subroutines( lines, dict_subs_modules ):
    """Modifies the list lines in-place by replacing 'use *' with 'use * only'"""
    current_subroutine = None
    i = 0
    N_lines = len(lines)
    while i < N_lines:
        line = lines[i]
        # Detect beginning of subroutine
        m = re.match('\s*subroutine (\w+)', line, re.IGNORECASE)
        if m:
            current_subroutine = m.group(1).lower()
        # Detect end of subroutine
        if re.match('\s*end subroutine', line, re.IGNORECASE):
            current_subroutine = None
        # Detect end of module
        if re.match('\s*end module', line, re.IGNORECASE):
            inside_module = False
        # Replace use
        if current_subroutine is not None:
            m = re.match('(\s*)use (\w+)', line, re.IGNORECASE)
            if m:
                # Find module name
                module = m.group(2).lower()
                # Collect complete the line, including line breaks
                complete_line = line
                while '&' in line:
                    lines[i] = '' # Erase line in final text
                    i += 1
                    line = lines[i]
                    complete_line = complete_line.rstrip('& \n') + ' ' + line.lstrip(' ')
                # Rewrite line
                if module in dict_subs_modules[ current_subroutine ]:
                    variable_list = dict_subs_modules[ current_subroutine ][module]
                    if variable_list == []:
                        # No used variables: the module does not need to be used
                        lines[i] = ''
                    else:
                        variables = ', '.join(variable_list)
                        indent = m.group(1)
                        new_line = 'USE %s, ONLY: %s\n' %(m.group(2), variables)
                        final_line = format_less_than_75_characters( new_line, indent )
                        lines[i] = final_line
                elif module not in known_external_modules:
                    raise RuntimeError('   Skipping %s in %s' %(module, current_subroutine))
        # Go the next line
        i += 1


def format_less_than_75_characters( line, indent ):
    words = line.split(' ')
    total_line = ''
    new_line = indent
    for word in words:
        new_line += word + ' '
        if len(new_line) > 75:
            total_line += new_line + ' &\n'
            new_line = indent
    total_line += new_line
    return total_line.rstrip(' ')


if __name__ == '__main__':


    # Go through all files and find modules and the variables defined
    dict_modules = {}
    for filename in glob.iglob('../../src/**/*.F90', recursive=True):
        print(filename)
        with open(filename) as f:
            lines = f.readlines()
        lines = reconstruct_lines(lines)
        get_module_variables(lines, dict_modules)
    print('')
    for key in dict_modules:
        print(key)

    # Go through all files and replace USE module syntax
    for filename in glob.iglob('../../src/**/*.F90', recursive=True):
        print(filename)
        with open(filename) as f:
            lines = f.readlines()
        # Extract the text of each subroutines
        dict_subs = get_subroutines(reconstruct_lines(lines))
        # Extract dictionaries with keys (subroutine, module) and
        # values list of variables used in the subroutine
        dict_subs_modules = get_sub_module( dict_subs, dict_modules )
        # Rewrite routine by using the "only" syntax ; this is done by modifying the list `lines`
        rewrite_subroutines( lines, dict_subs_modules )
        with open(filename, 'w') as f:
            f.write(''.join(lines) )
        print('')

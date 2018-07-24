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
known_external_modules = ['mpi', 'omp_lib', 'p3dfft', 'itt_sde_fortran', 'iso_c_binding']

def reconstruct_lines(lines):
    """Reconstruction full lines from Fortran broken lines using &; and without comments"""
    new_lines = []
    current_line = ''
    for line in lines:
        # Remove comments
        new_line = re.sub('!.*', '', line)
        new_line = new_line.rstrip(' ')
        if new_line.endswith('&'):
            new_line = new_line.rstrip('& ')
            current_line += new_line
        else:
            current_line += new_line
            new_lines.append(current_line)
            current_line = ''
    return new_lines


def get_module_variables(listlines, dict_modules, dict_used_modules):
    """Modifies `dict_modules` by adding entries, with keys the module name
    and value a list of existing module variables"""
    current_module = None
    inside_type = False
    inside_subroutine = False
    for line in listlines:
        # Detect beginning of type
        m = re.match('\s*type (\w+)', line, re.IGNORECASE)
        if m:
            dict_modules[current_module].append(m.group(1).lower())
            if re.search('::', line): #Single-line type declaration
                continue
            else:
                inside_type = True
        # Detect beginning of subroutine
        m = re.match('\s*subroutine (\w+)', line, re.IGNORECASE)
        if m and current_module is not None:
            dict_modules[current_module].append(m.group(1).lower())
            inside_subroutine = True
        # Detect module beginning
        m = re.match("^\s*module\s*(\w+)", line, flags=re.IGNORECASE)
        if m:
            if current_module is not None:
                raise ValueError('Still parsing %s' %line )
            current_module = m.group(1).lower()
            dict_modules[current_module] = []
            dict_used_modules[current_module] = []
        # Detect module end
        m = re.match("^\s*end module", line, flags=re.IGNORECASE)
        if m:
            current_module = None
        # Detect variable names and used modules
        if (not inside_type) and (not inside_subroutine) and (current_module is not None):
            # Detect variable name
            m = re.match('.*::(.*)', line)
            if m:
                for previous_char, variable in re.findall('(\W)(\w+)', m.group(1)):
                    # Exclude default values for variables, such .TRUE. and 'ex'
                    if previous_char not in ['=', '.', "'"]:
                        # Exclude pure numbers
                        if not re.match('[\d]', variable):
                            dict_modules[current_module].append( variable.lower() )
            # Detect used module
            m = re.match('\s*use (\w+)', line, re.IGNORECASE)
            if m:
                dict_used_modules[current_module].append(m.group(1).lower())
        # Detect end of type
        if re.match('\s*end type', line, re.IGNORECASE):
            inside_type = False
        # Detect end of subroutine
        if re.match('\s*end subroutine', line, re.IGNORECASE):
            inside_subroutine = False

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

def get_sub_module( dict_subs, dict_modules, dict_used_modules ):
    dict_subs_modules = {}
    dict_ifdef_modules = {}

    # Loop over the subroutines
    for name, text in dict_subs.items():
        # Initialize dictionaries
        dict_ifdef_modules[name] = {}
        dict_subs_modules[name] = {}
        # Loop over lines and find all modules
        lines = text.split('\n')
        ifdef = None
        module_list = []
        for line in lines:
            # Find beginning of ifdef
            m = re.match( '\s*#if defined\((\w+)\)', line, re.IGNORECASE )
            if m:
                ifdef = m.group(1)
            # Find modules
            m = re.match( '\s*use (\w+)', line, re.IGNORECASE )
            if m:
                module = m.group(1).lower()
                module_list.append( module )
                dict_ifdef_modules[name][module] = ifdef
            # Find end of ifdef
            if re.match('\s*#endif', line, re.IGNORECASE):
                ifdef = None
        # Also include all used modules
        module_set = set()
        for module in module_list:
            module_set.add(module)
            if module in known_external_modules:
                continue
            for used_module in dict_used_modules[module]:
                module_set.add( used_module )
                # Determine whether this used module should be imported outside of an ifdef
                if (dict_ifdef_modules[name][module] is None) or \
                    ((used_module in dict_ifdef_modules[name]) and (dict_ifdef_modules[name][used_module] is None)):
                    dict_ifdef_modules[name][used_module] = None
                else:
                    # Otherwise, check that it is using the same ifdef
                    if used_module in dict_ifdef_modules[name]:
                        assert dict_ifdef_modules[name][used_module] == dict_ifdef_modules[name][module]
                    else:
                        # Not yet imported: import it within the same ifdef
                        dict_ifdef_modules[name][used_module] = dict_ifdef_modules[name][module]

        # Loop over the modules used in this subroutine
        dict_subs_modules[name] = {}
        for module in module_set:
            dict_subs_modules[name][module] = []
            if module in known_external_modules:
                continue
            # Loop over the variables defined by this module
            for variable in dict_modules[module]:
                # Check whether this variable is being used
                if re.search('[\W_]%s[\W\n]'%variable, text, re.IGNORECASE):
                    dict_subs_modules[name][module].append(variable)

    return dict_subs_modules, dict_ifdef_modules

def rewrite_subroutines( lines, dict_subs_modules, dict_ifdef_modules ):
    """Modifies the list lines in-place by replacing 'use *' with 'use * only'"""
    current_subroutine = None
    i = 0
    N_lines = len(lines)
    replaced_modules = False
    while i < N_lines:
        line = lines[i]
        # Detect beginning of subroutine
        m = re.match('\s*subroutine (\w+)', line, re.IGNORECASE)
        if m:
            current_subroutine = m.group(1).lower()
            replaced_modules = False  # Did not yet replace modules
        # Detect end of subroutine
        if re.match('\s*end subroutine', line, re.IGNORECASE):
            current_subroutine = None
        # Detect end of module
        if re.match('\s*end module', line, re.IGNORECASE):
            inside_module = False
        if current_subroutine is not None:
            # Find first module
            m = re.match('(\s*)use (\w+)', line, re.IGNORECASE)
            if m:
                # Remove all modules lines
                while re.match('(\s*)use (\w+)', line, re.IGNORECASE):
                    # Remove ifdef before and after
                    if re.match('\s*#if defined', lines[i-1], re.IGNORECASE):
                        lines[i-1] = ''
                    if re.match('\s*#endif', lines[i+1], re.IGNORECASE):
                        lines[i+1] = ''
                    # Collect complete the line, including line breaks
                    while '&' in line:
                        lines[i] = '' # Erase line in final text
                        i += 1
                        line = lines[i]
                    lines[i] = ''
                    i += 1
                    line = lines[i]
                i -= 1
                # Add the used modules (if it was not yet done)
                indent = m.group(1)
                if not replaced_modules:
                    for module in dict_subs_modules[current_subroutine].keys():
                        # Add ifdef if needed
                        if dict_ifdef_modules[current_subroutine][module] is not None:
                            lines[i] += '#if defined(%s)\n' %dict_ifdef_modules[current_subroutine][module]
                        # Write external modules
                        if module in known_external_modules:
                            new_line = '%sUSE %s\n' %(indent, module)
                            lines[i] += new_line
                        # Otherwise add variables to the current line
                        else:
                            variable_list = dict_subs_modules[ current_subroutine ][module]
                            if variable_list != []:
                                variables = ', '.join(variable_list)
                                new_line = 'USE %s, ONLY: %s\n' %(module, variables)
                                final_line = format_less_than_75_characters( new_line, indent )
                                lines[i] += final_line
                        # Add ifdef if needed
                        if dict_ifdef_modules[current_subroutine][module] is not None:
                            lines[i] += '#endif %s\n' %dict_ifdef_modules[current_subroutine][module]
                    replaced_modules = True
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
            new_line = indent*2
    total_line += new_line
    return total_line.rstrip(' ')

def remove_empty_endif( lines ):
    N_lines = len(lines)
    # Erase empty if / endif
    for i in range(N_lines):
        # Check the if defined
        m = re.match('\s*#if defined\((\w+)\)', lines[i])
        if m:
            # Check if the next line is already the next endif
            if re.match('\s*#endif %s' %m.group(1), lines[i+1]):
                # In this case erase both
                lines[i] = ''
                lines[i+1] = ''
    # Erase unneeded match endif/if
    for i in range(N_lines):
        # Check the if defined
        m = re.match('\s*#if defined\((\w+)\)', lines[i])
        if m:
            # Check if the previous line is a corresponding endif
            if re.match('\s*#endif %s' %m.group(1), lines[i-1]):
                # In this case erase both
                lines[i] = ''
                lines[i-1] = ''
    # Erase names after endif
    for i in range(N_lines):
        m = re.match('(\s*#endif)', lines[i])
        if m:
            lines[i] = m.group(1) + '\n'

if __name__ == '__main__':


    # Go through all files and find modules and the variables defined
    dict_modules = {}
    dict_used_modules = {}
    print('Scanning modules...')
    for filename in glob.iglob('../../src/**/*.F90', recursive=True):
        with open(filename) as f:
            lines = f.readlines()
        lines = reconstruct_lines(lines)
        get_module_variables(lines, dict_modules, dict_used_modules)
    print('')

    # Go through all files and replace USE module syntax
    for filename in glob.iglob('../../src/**/*.F90', recursive=True):
        print('Rewriting %s' %filename)
        with open(filename) as f:
            lines = f.readlines()
        # Extract the text of each subroutines
        dict_subs = get_subroutines(reconstruct_lines(lines))
        # Extract dictionaries with keys (subroutine, module) and
        # values list of variables used in the subroutine
        dict_subs_modules, dict_ifdef_modules = \
            get_sub_module( dict_subs, dict_modules, dict_used_modules )
        # Rewrite routine by using the "only" syntax ; this is done by modifying the list `lines`
        rewrite_subroutines( lines, dict_subs_modules, dict_ifdef_modules )
        with open(filename, 'w') as f:
            f.write(''.join(lines) )
        print('')

    # Clean empty endif
    for filename in glob.iglob('../../src/**/*.F90', recursive=True):
        print('Cleaning %s' %filename)
        with open(filename) as f:
            lines = f.readlines()
        remove_empty_endif(lines)
        with open(filename, 'w') as f:
            f.write(''.join(lines) )

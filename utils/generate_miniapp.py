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
                & (curr_line.find("end subroutine")==-1) \
                & (curr_line.find("#do not parse")==-1) ):
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
                        if (subname.find("#do not parse")==-1):
                            break
                        else:
                            pass

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
        iparseinstr=line.find("!#do not parse")
        iwrapinstr=line.find("!#do not wrap")
        if ((icomm >=0) & \
            (iparseinstr==-1) & \
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

        #LIST ALL .F90 or .F files in current directory
        self.listfiles = self.create_listfiles()

        for file in self.listfiles:
            self.write_available_routines(file)

    def clean_folder(self):

        # Delete all files
        if os.path.exists('./PICSARlite'):
            shutil.rmtree('./PICSARlite')

    def create_listfiles(self):
        listfiles = []
        for root, subFolder, files in os.walk('./src'):
            for file in files:
                if len(root) > 5:
                     listfiles.append('%s/%s'%(root[6:], file))
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

        # Just indenpendant modules or subroutines
        if len(nb_modules) == 0:
            # Loop on the available modules
            for iname in range(len(m.names)):
                if self.list_available_modules is None:
                    self.copy_files_from_picsar(m, iname)

                elif m.names[iname] in self.list_available_modules:
                    self.copy_files_from_picsar(m, iname)

            for iname in range(len(r.names)):
                if self.list_available_routines is None:
                    self.copy_files_from_picsar(r, iname)

                elif r.names[iname] in self.list_available_routines:
                    self.copy_files_from_picsar(r, iname)

        # Loop on each module occurence, for each print the modules then
        # subroutines inside and finish by "end module"
        else:
            for module in nb_modules:
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

    def copy_files_from_picsar( self, routine, index ):
        # routine is an instance of the class Routines or Modules

        file = routine.file

        if not os.path.exists('PICSARlite/src/%s'%file):

            divided_path = file.split('/')
            path_folder = 'PICSARlite/src/'

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

            fnew = open('PICSARlite/src/%s'%file, 'w')
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

        f = open('PICSARlite/src/%s'%file, 'aw')
        f.writelines('END MODULE %s'%module +'\n')
        f.close()

###############################################################################
# Main

arglist=sys.argv
try:
    type_parser = str(arglist[1])
    type_pusher = str(arglist[2])
    type_depos  = str(arglist[3])
    miniapp = MiniAppParser(type_parser, type_pusher, type_depos)

except(IndexError):
    print('##################################################################')
    print('Some arguments are needed there.')
    print('The list of available solvers is:')
    print('- all: every routines and modules from picsar')
    print('- fdtd: every routines and modules related to the Finate Domain')
    print('        Time Domain Maxwell solver in 3D.')
    print('- spectral: every routines and modules related to the Pseudo-')
    print('        Spectral Analytical Time Domain Maxwell solver in 3D.')
    print('##################################################################')

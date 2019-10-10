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

 PARSER FOR FORTHON
 H. VINCENTI - DEC 7, 2015


 This function reads a group of Fortran 90 files and produces a FORTHON
 interface file .v used by the forthon compiler to produce a python
 module file .so

 Input arguments: filename

 Outputs:
 1. New fortran file filename_forthon.F90 forthon compliant
 2. interface file .v with subroutines/modules declaration
 NB :
 - User has to define variables  in a Fortran MODULE
 - Routines in procedure modules are not treated differently than routines outside modules
 - Variables in a Module are grouped inside a common "Module" of variables in the .v file
 - Derived types have to be declared in modules of same name (see FORTHON doc)

 Developers:
 Henri Vincenti
 Mathieu Lobet

 Date:
 Creation: Dec 7 2015

 REVISION:
 - Mathieu Lobet - 08.18.2016 - better management of the directives

 _______________________________________________________________________________
"""

import numpy as np
import sys
import os
import time


###### MAIN FUNCTION FOR PRE-PARSING SEVERAL .F90 FILES INTO A SINGLE FILE

def fortran_preparser(filelist,filewrite):

    ##### Open new file to write modules
    # Subroutine file
    fws=open(filewrite+"_subroutines.F90","w")
    # Module file
    fwm=open(filewrite+"_modules.F90","w")

    # Parsing modules and subroutines for all files
    listlines_modules=[]
    listlines_subroutines=[]
    for ifile in range(0,len(filelist)):
        ##### Open Fortran source file
        fr=open(filelist[ifile],"r")
        ##### Cleaning text from comments, continuations &
        listlines=fr.readlines()
        Nlines=len(listlines)

        print("Pre-formatting "+str(Nlines)+" lines of file "+filelist[ifile])
        listlines_1=(preprocess_file(listlines))
        print("Sanity check of file "+filelist[ifile])
        #listlines_2=sanity_check_1(listlines_1)
        #listlines=sanity_check_3(listlines_1)
        listlines=listlines_1
        ##### Parse modules
        print("Now parsing modules of file "+filelist[ifile])
        listlines_modules=listlines_modules+preparse_modules(listlines)

         ##### Parse subroutines
        print("Now parsing subroutines of file "+filelist[ifile])
        listlines_subroutines=listlines_subroutines+preparse_subroutines(listlines)

    # For each subroutine/function, check if an interface A is needed
    # If needed,  place interface in a module intmod
    # All modules are placed at the end of module files
    # In subroutines wher A is needed append "use intmod"
    [listlines_interfaces, listlines_subroutines]=preparse_subroutine_interfaces(listlines_subroutines)
    listlines_modules=listlines_modules+listlines_interfaces
    ###### Write all modules  to file
    fwm.writelines(listlines_modules)

    ###### Write all subroutines to file
    fws.writelines(listlines_subroutines)


def preparse_modules(listlines):
    listlines_new=[]
    # Identify module blocks
    [names,istart,iend,procmod] = get_module_blocks(listlines)
    for i in range(0,len(names)):
        listlines_new.append("!------------------\n")
        listlines_new.append("! THIS IS A MODULE \n")
        listlines_new.append("!------------------\n")
        if (procmod[i]):
            listlines_new.append("module "+names[i].replace("moduletype","")+" #do not parse"+"\n")
        else:
            listlines_new.append("module "+names[i].replace("moduletype","")+"\n")
        for iline in range(istart[i]+1,iend[i]):
            listlines_new.append(listlines[iline]+"\n")
        listlines_new.append("end module "+names[i].replace("moduletype","")+"\n")
    return listlines_new


def preparse_subroutines(listlines):
    listlines_new=[]
    # Identify subroutine blocks
    [names,istart,iend,proceduremod] = get_subroutine_blocks(listlines)
    for i in range(0,len(names)):
        listlines_new.append("!----------------------\n")
        listlines_new.append("! THIS IS A SUBROUTINE \n")
        listlines_new.append("!----------------------\n")
        listlines_new.append(listlines[istart[i]]+"\n")
        if (proceduremod[i]!=""):
                listlines_new.append("use "+proceduremod[i]+"\n")
        for iline in range(istart[i]+1,iend[i]+1):
            listlines_new.append(listlines[iline]+"\n")

    # Identify function blocks
#    [names,istart,iend,proceduremod] = get_function_blocks(listlines)
#    for i in range(0,len(names)):
#        listlines_new.append("!----------------------\n")
#        listlines_new.append("! THIS IS A FUNCTION   \n")
#        listlines_new.append("!----------------------\n")
#        listlines_new.append(listlines[istart[i]]+"\n")
#        if (proceduremod[i]!=""):
#            listlines_new.append("use "+proceduremod[i]+"\n")
#        for iline in range(istart[i]+1,iend[i]+1):
#            listlines_new.append(listlines[iline]+"\n")


    return listlines_new


def preparse_subroutine_interfaces(listlines_sub):

    # Identify subroutines that need interface and create corresponding interface block
    [names,istart,iend,proceduremod] = get_subroutine_blocks(listlines_sub)
    listlines_interfaces=[]
    isininterface=[False for _ in range(len(names))]
    for i in range(0,len(names)):
        [argsname, argstype, wrap] = parse_subroutine_args(preprocess_file(listlines_sub[istart[i]:iend[i]+1]))
        needsinterface=False
        for iargs in range(0, len(argstype)):
            currtype=argstype[iargs]
            #print(currtype.strip())
            ipoint=currtype.strip().find("_")
            # This subroutine needs an interface block
            if (ipoint==0):
                needsinterface=True
        if (needsinterface):
            isininterface[i]=True
            # Get subroutine tag
            subtag=[]
            for itag in range(istart[i]+1,iend[i]):
                if (listlines_sub[itag].find("use ")>=0):
                    subtag.append(listlines_sub[itag])
                if (listlines_sub[itag].find("::")>=0):
                    subtag.append(listlines_sub[itag])
                if (listlines_sub[itag].find("implicit none")>=0):
                    subtag.append(listlines_sub[itag])
            # Creates interface block
            interf_name=str(time.clock()).replace(".","")
            listlines_interfaces.append("!----------------------------------\n")
            listlines_interfaces.append("!This is an interface module block \n")
            listlines_interfaces.append("!----------------------------------\n")
            listlines_interfaces.append("module interf_"+names[i].strip()+" #do not parse\n")
            listlines_interfaces.append("interface intef"+interf_name+"\n")
            listlines_interfaces.append(listlines_sub[istart[i]]+"\n")
            listlines_interfaces=listlines_interfaces+subtag
            listlines_interfaces.append(listlines_sub[iend[i]]+"\n")
            listlines_interfaces.append("end interface intef"+interf_name+"\n")
            listlines_interfaces.append("end module interf_"+names[i].strip()+"\n")

    # Add "use interfacemodulename" to subroutines calling an interfaced subroutine
    listlines_sub_new=listlines_sub[:]
    for i in range(0,len(names)):
        for iline in range(istart[i],iend[i]+1):
            currline=listlines_sub[iline]
            # This is a call to a subroutine
            if ((currline.find("call "))>=0):
                calltoint=False
                interface_name=""
                for iname in range(0,len(names)):
                    # Needs to add "use interfacemodulename' in this subroutine
                    if ((currline.find(names[iname])>=0) & (isininterface[iname])):
                        calltoint=True
                        interface_name=names[iname]
                        break
                if (calltoint):
                    indx=listlines_sub_new.index(listlines_sub[istart[i]])
                    listlines_sub_new.insert(indx+1,"use interf_"+interface_name.strip()+"\n")
    return [listlines_interfaces, listlines_sub_new]


# This functions performs sanity check for Forthon compliance
# ** In particular:
# - Check 1: Fortran derived types have to be declared in separate modules with name Typename+"Module"
def sanity_check_1(listline):
    lenlist=len(listline)
    ## Identify module blocks and positions in file
    [names,istart,iend,procmod] = get_module_blocks(listline)
    ## Get type definition blocks
    [names_t,istart_t,iend_t] = get_DerivedType_def_blocks(listline)

    ## Sanity check 1: Fortran derived types definition blocks have to be placed in separate modules with names
    # Find type declaration that does not comply with Forthon rules
    compliant =[False for _ in range(len(names_t))]
    for i in range(0,len(names_t)):
        for imod in range(0,len(names)):
            # Type def block is in module imod
            if ((istart_t[i]>=istart[imod]) & (iend_t[i]<=iend[imod])):
                # Check if module name is typename+"module"
                if (names[imod].find(names_t[i].strip()+"module") >=0):
                    compliant[i]=True

    # Re-write new listline without non compliant Type declaration
    listline_new=[]
    for i in range(0,len(listline)):
        iscompliant=True
        for idecl in range(0,len(names_t)):
            if ((i>=istart_t[idecl]) & (i<=iend_t[idecl]) & (not compliant[idecl])):
                iscompliant=False
        if (iscompliant):
            listline_new.append(listline[i])
    # Append non compliant type declaration in modules at beginning of file
    new_type_mod=[]
    for i in range(0,len(names_t)):
        if (not compliant[i]):
            new_type_mod.append("module "+names_t[i]+"module")
            for iblock in range(istart_t[i],iend_t[i]+1):
                new_type_mod.append(listline[iblock])
            new_type_mod.append("end module ")
    listline_new=new_type_mod+listline_new

    return listline_new

# - Check 2: Arrays of derived types are not supported.
#    + If present in a module xxx the module won't be parsed is not parsed
def sanity_check_2(listline):
    lenlist=len(listline)
    ## Identify module blocks and positions in file
    [names,istart,iend,procmod] = get_module_blocks(listline)

    ## Get type definition blocks
    [names_t,istart_t,iend_t] = get_DerivedType_def_blocks(listline)

    ## Get array of derived type declaration lines
    ideclines=[]
    for i in range(0,len(listline)):
        curr_line=listline[i]
        idecl=curr_line.find("::")
        if ((idecl >= 0) & (curr_line[0:idecl].find("type(")>=0)):
            [argsname, argstype, initvalue]=parse_decl_line(curr_line)
            # This is an array of type declaration line
            for id in range(0,len(argsname)):
                # This is an array
                if (argsname[id].find("(")>=0):
                    ideclines.append(i)
                    break
    # Detect arrays of derived types in type declaration
    isintype=[False for _ in range(len(ideclines))]
    isinmodule=[False for _ in range(len(ideclines))]
    indtype=[-1 for _ in range(len(ideclines))]
    indmod=[-1 for _ in range(len(ideclines))]
    for i in range(0,len(ideclines)):
        for idecl in range(0,len(names_t)):
            # Array of derived types in a type def block
            if ((ideclines[i]>=istart_t[idecl]) & (ideclines[i]<=iend_t[idecl])):
                isintype[i]=True
                indtype[i]=istart_t[idecl]
        for ibl in range(0,len(names)):
            # Array of derived types outside a type def block but in a module (can be in a subroutine)
            if ((ideclines[i]>=istart[ibl]) & (ideclines[i]<=iend[ibl])):
                isinmodule[i]=True
                indmod[i]=istart[ibl]

    # Re-write new listline without non compliant Type declaration
    listline_new=[]
    for i in range(0,len(listline)):
            # Array of derived
            if ((not i in ideclines)):
                listline_new.append(listline[i])
            else:
                indx=ideclines.index(i)
                if (isintype[indx]):
                    indmodtype=listline_new.index(listline[indmod[indx]])
                    listline_new.append(listline[i])
                    listline_new[indmodtype]=listline_new[indmodtype]+" #do not parse"
                else:
                    if (not isinmodule[indx]):
                        listline_new.append(listline[i])

    # Add array of derived types defined in modules inside their own "unused" modules
    for i in range(0,len(ideclines)):
        if (isinmodule[i] and (not isintype[i])):
            # Find module in which
            moduleline=listline[indmod[i]]
            modulename="unused"+str(time.clock()).replace(".","")
            # Insert use module modulename at the beginning of the containing module
            listline_new.insert(listline_new.index(moduleline)+1,"use "+modulename)
            # Create a new module containing array of derived type
            listline_new.insert(listline_new.index(moduleline),"module "+modulename+ " #do not parse")
            listline_new.insert(listline_new.index(moduleline),listline[ideclines[i]])
            listline_new.insert(listline_new.index(moduleline),"end module "+modulename)


    return listline_new

# - Check 3: Arrays of derived types are not supported.
#    + If present in a module xxx the module won't be parsed is not parsed
def sanity_check_3(listline):
    lenlist=len(listline)
    ## Identify module blocks and positions in file
    [names,istart,iend,procmod] = get_module_blocks(listline)

    ## Get type definition blocks
    [names_t,istart_t,iend_t] = get_DerivedType_def_blocks(listline)

    ## Get array of derived type declaration lines
    ideclines=[]
    for i in range(0,len(listline)):
        curr_line=listline[i]
        idecl=curr_line.find("::")
        if ((idecl >= 0) & (curr_line[0:idecl].find("type(")>=0)):
            [argsname, argstype, initvalue]=parse_decl_line(curr_line)
            # This is an array of type declaration line
            for id in range(0,len(argsname)):
                # This is an array
                if (argsname[id].find("(")>=0):
                    ideclines.append(i)
                    break
    # Detect arrays of derived types in modules
    isinmodule=[False for _ in range(len(ideclines))]
    indmod=[-1 for _ in range(len(ideclines))]
    for i in range(0,len(ideclines)):
        for ibl in range(0,len(names)):
            # Array of derived types outside a type def block but in a module (can be in a subroutine)
            if ((ideclines[i]>=istart[ibl]) & (ideclines[i]<=iend[ibl])):
                isinmodule[i]=True
                indmod[i]=istart[ibl]

    # Re-write new listline with flag do not parse
    listline_new=[]
    for i in range(0,len(listline)):
            # Array of derived
            if ((not i in ideclines)):
                listline_new.append(listline[i])
            else:
                indx=ideclines.index(i)
                if (isinmodule[indx]):
                    indmodtype=listline_new.index(listline[indmod[indx]])
                    listline_new.append(listline[i])
                    listline_new[indmodtype]=listline_new[indmodtype]+" #do not parse"
                else:
                    listline_new.append(listline[i])

    return listline_new







def get_DerivedType_def_blocks(listlines):
    names=[]
    istart=[]
    iend=[]
    Nlines=len(listlines)
    for i in range(0,Nlines):
        curr_line=listlines[i]
        curr_line=rm_comments(curr_line)
        curr_word_list=curr_line.split(" ")
        # This is a type declaration block
        proc_mod=False
        if (("type" in curr_word_list) & (curr_line.find("(")==-1) & (curr_line.find("end type")==-1)):
            subname=curr_line[curr_line.find("type")+len("type"):len(curr_line)]
            istart.append(i)
            while True: # get end of block
                i=i+1
                if (i>=Nlines):
                    sys.exit("ERROR: missing end module block statement")
                curr_line=listlines[i]
                if ((curr_line.find("end type")>=0)):
                    break

            iend.append(i)
            names.append(subname.strip())
    return[names,istart,iend]




###### MAIN FUNCTION FOR PARSING an application FILE INTO .v file
def parse_forthon_app(appname):

    # Open Subroutines file
    fr1=open(appname+"_subroutines.F90","r")
    listlines_subroutines=fr1.readlines()

    # Open Modules file
    fr2=open(appname+"_modules.F90","r")
    listlines_modules=fr2.readlines()

    ##### Open appname.v file
    # Get file root name
    filewrite=appname+".v"
    # Open .v file
    fw=open(filewrite,"w")
    # Write down appname as python module name in .v file
    fw.write(appname+"\n")
    fw.write("\n")

    ##### Parse the Fortran source file in .v file
    ## Number of lines of the source file
    listlines_subroutines=(preprocess_file(listlines_subroutines))
    listlines_modules=(preprocess_file(listlines_modules))
    ## Identify and sparse blocks in order
    # - 1 Parse module blocks
    [namesmod,istartmod,iendmod,procmod] = get_module_blocks(listlines_modules)
    listlines_modules=parse_module_blocks(fw,listlines_modules, namesmod, istartmod, iendmod)

    # - 2 Parse subroutine blocks
    [names,istart,iend,proceduremod] = get_subroutine_blocks(listlines_subroutines)
    parse_subroutine_blocks(fw,listlines_subroutines, names, namesmod, istart, iend)

    # - 3 Parse final .F90 file used by Forthon
    # Open appname.F90 file
    fr3=open(appname+".F90","w")
    # Only write modules that were not parsed in .v
#    for i in range(0,len(namesmod)):
#        if (namesmod[i].find("#do not parse")>=0):
#            fr3.writelines(reformat_file(listlines_modules[istartmod[i]:iendmod[i]+1]))
    fr3.writelines(reformat_file(listlines_modules))
    # Write all subroutines
    fr3.writelines(reformat_file(listlines_subroutines))


# Function for parsing a file
def parse_file(fileread):
    ##### Open Fortran source file
    fr=open(fileread,"r")

    ##### Open corresponding .v file
    # Get file root name
    beg=fileread.find(".")
    end=len(fileread)
    fileroot=fileread.replace(fileread[beg:end],"")
    filewrite=fileroot+".v"
    # Open .v file
    fw=open(filewrite,"w")
    # Write down module name in .v file
    fw.write(fileroot+"\n")
    fw.write("\n")


    ##### Parse the Fortran source file in .v file
    ## Number of lines of the source file
    listlines=fr.readlines()
    Nlines=len(listlines)
    print("Pre-processing "+str(Nlines)+" lines of file "+fileread)
    listlines=(preprocess_file(listlines))

    ## Identify and parse blocks in order
    # - 1 Parse module blocks
    [namesmod,istart,iend,procmod] = get_module_blocks(listlines)
    parse_module_blocks(fw,listlines, namesmod, istart, iend)
    # - 2 Parse subroutine blocks
    [names,istart,iend,proceduremod] = get_subroutine_blocks(listlines)
    parse_subroutine_blocks(fw,listlines, names, namesmod, istart, iend)

###### ----------------------------------------------------
###### FUNCTIONS FOR MODULE BLOCK PARSING
#### - INDENTIFY MODULE BLOCKS THAT NEED TO BE PARSED
def get_module_blocks(listlines):
    names=[]
    istart=[]
    iend=[]
    proceduremod=[]
    Nlines=len(listlines)
    for i in range(0,Nlines):
        curr_line=listlines[i]
        curr_line=rm_comments(curr_line)
        curr_word_list=curr_line.split(" ")
        # This is a module block
        proc_mod=False
        if (("module" in curr_word_list) & (curr_line.find("use ")==-1) & (curr_line.find("end module")==-1)):
            subname=curr_line[curr_line.find("module")+len("module"):len(curr_line)]
            istart.append(i)
            while True: # get end of block
                i=i+1
                if (i>=Nlines):
                    sys.exit("ERROR: missing end module block statement")
                curr_line=listlines[i]
                curr_word_list=curr_line.split(" ")
                if (("type" in curr_word_list) & (not "end" in curr_word_list)):
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
    return[names,istart,iend,proceduremod]

#### - PARSE EACH SUBROUTINE BLOCK
def parse_module_blocks(fw,listlines,names,istart,iend):
    nblocks=len(names)
    listlines_new=[]
    # Get each module block and parse it
    print("Parsed "+str(nblocks)+" module block(s) in file .v file\n")
    for iblock in range(0,nblocks):
        fw.write("\n")
        # Modules with "#do not parse" are not parsed in .v file
        if (names[iblock].find("#do not parse")==-1):
            # This is a Type module block
            if(names[iblock].find("moduletype")>=0):
                modinv=[]
                # Print module in .v with % instead of #
                names[iblock]=names[iblock].replace("moduletype","")
                names[iblock]=names[iblock].replace("module","")
                fw.write("%%%%% "+names[iblock].strip()+":\n")
                modinv.append("%%%%% "+names[iblock].strip()+":\n")
                print("derived type "+names[iblock])
                # Parse module variable names and types
                [listvars,listtypes,listinitv] = parse_module_vars(listlines[istart[iblock]:iend[iblock]+1])
                # Write parsed module in .v file
                for i in range(0,len(listvars)):
                    if (listvars[i].find("cobj__")==-1):
                        if(listinitv[i]==""):
                            fw.write(listvars[i]+" "+listtypes[i]+"\n")
                            modinv.append(listvars[i]+" "+listtypes[i]+"\n")
                        else:
                            fw.write(listvars[i]+" "+listtypes[i]+" /"+listinitv[i]+"/ \n")
                            modinv.append(listvars[i]+" "+listtypes[i]+"\n")
                # Generate Fortran module procedure for derived type
                # using Forthon
                newmod=gen_forthon_dtypef90(names[iblock],modinv)
                listlines_new=listlines_new+newmod
            # This is a standard module block
            else:
                fw.write("***** "+names[iblock].strip()+":\n")
                print("module "+names[iblock])
                # Parse module variable names and types
                [listvars,listtypes,listinitv] = parse_module_vars(listlines[istart[iblock]:iend[iblock]+1])
                # Write parsed module in .v file
                for i in range(0,len(listvars)):
                    if (listvars[i].find("cobj__")==-1):
                        if(listinitv[i]==""):
                            fw.write(listvars[i]+" "+listtypes[i]+"\n")
                        else:
                            fw.write(listvars[i]+" "+listtypes[i]+" /"+listinitv[i]+"/ \n")
                listlines_new=listlines_new+listlines[istart[iblock]:iend[iblock]+1]
        else:
            listlines_new=listlines_new+listlines[istart[iblock]:iend[iblock]+1]
    return listlines_new

def gen_forthon_dtypef90(modname,modinv):
    fmod=open("tempapp"+".F90","w")
    fmod.writelines("! FORTRAN FILE")
    fmod.close()
    fv=open("tempapp"+".v","w")
    fv.write("tempapp\n")
    fv.write("")
    fv.writelines(modinv)
    fv.close()
    os.system("Forthon -v --no2underscores -g tempapp")
    dirs=[x[0] for x in os.walk("build/")]
    print(dirs)
    filedir=dirs[3]
    fnewmod=open(filedir+"/tempapp_p.F90","r")
    newmod_pre= fnewmod.readlines()
    newmod_post=[]
    # Only get procedure module (not additional subroutines)
    isinmod=True
    i=0
    while(isinmod):
      newmod_post.append(newmod_pre[i])
      i=i+1
      if (newmod_pre[i].find("END MODULE")>=0):
        newmod_post.append(newmod_pre[i])
        isinmod=False

    fnewmod.close()
    os.system("rm -rf build/")
    os.system("rm -f tempapp.v")
    os.system("rm -f tempapp.F90")
    return newmod_post


### - PARSE SUBROUTINE ARGS AND TYPES
def parse_module_vars(listlines):
    len_block=len(listlines)
    listvars=[]
    listtypes=[]
    listinitv=[]
    # Get module vars, types and init values
    for i in range(0,len_block):
        curr_line=listlines[i]
        idecl=curr_line.find("::")
        if (idecl >= 0):
            [argsname, argstype, initvalue]=parse_decl_line(curr_line)
            for i in range(0,len(argsname)):
                listvars.append(argsname[i])
                listtypes.append(argstype[i])
                listinitv.append(initvalue[i])

    return [listvars,listtypes,listinitv]
###### ----------------------------------------------------



###### ----------------------------------------------------
###### FUNCTIONS FOR SUBROUTINE BLOCK PARSING
#### - INDENTIFY SUBROUTINE BLOCKS THAT NEED TO BE PARSED
def get_subroutine_blocks(listlines):
    names=[]
    istart=[]
    iend=[]
    proceduremod=[]
    Nlines=len(listlines)
    proc_mod=True
    modulename=""
    for i in range(0,Nlines):
        curr_line=listlines[i]
        curr_word_list=curr_line.split(" ")
        # We are in a module
        if (("module" in curr_word_list) & (curr_line.find("use")==-1) & (curr_line.find("end module")==-1)):
            modulename=curr_word_list[curr_word_list.index("module")+1]
        # We are in a procedure module
        if(curr_line.find("contains")>=0 ):
            proc_mod=True
        # We are out of a module
        if ((curr_line.find("end module")>=0)):
            proc_mod=False
            modulename=""
        # This is a subroutine block
        if ((curr_line.find("subroutine")>=0) & (curr_line.find("end subroutine")==-1) \
        & (curr_line.find("#do not parse")==-1) ):
            ip=curr_line.find("(")
            # This is a subroutine in a procedure module
            if (proc_mod):
                proceduremod.append(modulename)
            else:
                proceduremod.append("")
            if (ip == -1): # Case no parenthesis
                ip=len(curr_line)
            subname=curr_line[curr_line.find("subroutine")+len("subroutine"):ip]
            names.append(subname.strip())
            istart.append(i)
            search_end=True
            while True: # get end of block
                i=i+1
                if (i>=Nlines):
                    sys.exit("ERROR: missing end subroutine block")
                curr_line=listlines[i]

                # If an interface is defined in the subroutine
                # We pass the interface block without looking for "end subroutine"
                if (curr_line.find("interface")>=0):
                  search_end = False
                if (curr_line.find("end interface")>=0):
                  search_end = True
                if ((curr_line.find("end subroutine")>=0)and(search_end)):
                    break
            iend.append(i)
    return[names,istart,iend,proceduremod]


def get_function_blocks(listlines):
    names=[]
    istart=[]
    iend=[]
    proceduremod=[]
    Nlines=len(listlines)
    proc_mod=True
    modulename=""
    for i in range(0,Nlines):
        curr_line=listlines[i]
        curr_word_list=curr_line.split(" ")
        # We are in a module
        if (("module" in curr_word_list) & (curr_line.find("use")==-1) & (curr_line.find("end module")==-1)):
            modulename=curr_word_list[curr_word_list.index("module")+1]
        # We are in a procedure module
        if(curr_line.find("contains")>=0 ):
            proc_mod=True
        # We are out of a module
        if ((curr_line.find("end module")>=0)):
            proc_mod=False
            modulename=""
        # This is a subroutine block
        if ((curr_line.find("function")>=0) & (curr_line.find("end function")==-1)):
            ip=curr_line.find("(")
            # This is a function in a procedure module
            if (proc_mod):
                proceduremod.append(modulename)
            else:
                proceduremod.append("")
            if (ip == -1): # Case no parenthesis
                ip=len(curr_line)
            subname=curr_line[curr_line.find("function")+len("function"):ip]
            names.append(subname.strip())
            istart.append(i)
            while True: # get end of block
                i=i+1
                if (i>=Nlines):
                    sys.exit("ERROR: missing end function block")
                curr_line=listlines[i]
                if (curr_line.find("end function")>=0):
                    break
            iend.append(i)

    return[names,istart,iend,proceduremod]








#### - PARSE EACH SUBROUTINE BLOCK
def parse_subroutine_blocks(fw,listlines,names,namesmod,istart,iend):
    nblocks=len(names)
    fw.write("\n")
    fw.write("***** Subroutines:\n")
    print("Parsed "+str(nblocks)+" subroutine block(s) in .v file \n")

    # Get non used derived type
    nonusedtype=[]
    for i in range(0,len(namesmod)):
        if ((namesmod[i].find("module")>=0) & (namesmod[i].find("#do not parse")>=0)):
            nonusedtype.append(namesmod[i].replace("moduletype","").replace("#do not parse","").replace("module","").replace(" ",""))

    # Get each subroutine block and parse it
    for iblock in range(0,nblocks):
        print("subroutine "+names[iblock])
        # Parse subroutine args names and types
        [argsname, argstype, wrap] = parse_subroutine_args(listlines[istart[iblock]:iend[iblock]+1])
        parse=True
        for itype in range(0,len(argsname)):
            for inon in range(0,len(nonusedtype)):
                if (argstype[itype].find(nonusedtype[inon])>=0):
                    parse=False
        if (parse and wrap):
            # Write parsed subroutine in .v file
            semicol=""
            if (argsname!=[""]):
                semicol=":"
            fw.write(names[iblock].strip()+"("+argsname[0].strip()+semicol+argstype[0].strip().replace("_istarget_",""))
            for i in range(1,len(argsname)):
                fw.write(","+argsname[i]+semicol+argstype[i].replace("_istarget_",""))
            fw.write(") subroutine\n")

### - PARSE SUBROUTINE ARGS AND TYPES
def parse_subroutine_args(listlines):
    len_block=len(listlines)

    # Get subroutine full arg list (without &, spaces and parenthesis)
    curr_line=listlines[0].replace(" ","")
    istart=curr_line.find("(")
    iend=curr_line.find(")")
    if (istart>=0):
        curr_line=curr_line[istart+1:iend]
    else: # subroutine has no arguments in brackets ()
        curr_line=""
    # Parse full arg chain
    argslist=curr_line.split(",")

    # Get args types
    argstypes=["" for _ in range(len(argslist))]
    wrap = True
    for i in range(0,len_block):
        curr_line=listlines[i]
        if curr_line.find("#do not wrap") >= 0:
            wrap = False
        idecl=curr_line.find("::")
        if (idecl >= 0):
            [argsname, argstype, initvalue]=parse_decl_line(curr_line)
            for ivar in range(0,len(argsname)):
                varname=argsname[ivar].strip()
                ideb=varname.find("(")
                if (ideb>=0): # var is an array
                    varname=varname[0:ideb] # remove array def from varname
                    ist= curr_line.find(varname+"(")
                # add this to type
                if varname in argslist:
                    idx=argslist.index(varname)
                    argslist[idx]=argsname[ivar] # in Case of array
                    argstypes[idx]=argstype[ivar].replace(" ","")

    return [argslist,argstypes,wrap]

###### ----------------------------------------------------

###### ----------------------------------------------------
######## GENERAL PURPOSE FUNCTIONS
### preprocess the file
### Remove comments and continuing files
def preprocess_file(listline):
    lenlist=len(listline)
    # First removes comments and new line
    # Make everything lower case (Fortran is case insensitive)
    for i in range(0,lenlist):
        curr_line=listline[i]
        # Line is not considered
        if (curr_line.find("!#do not consider")>=0):
          curr_line = ""
        else:
          # Compiler directive are case sensitive do not lower
          if not ((curr_line.find("#ifdef")>=0) or (curr_line.find("#if")>=0) or (curr_line.find("#else")>=0) or \
          (curr_line.find("#endif")>=0)):
              curr_line=curr_line.lower()
          curr_line=rm_comments(curr_line)
          curr_line=rm_newline(curr_line)
          curr_line=rm_spechar(curr_line)
        listline[i]=curr_line

    # Concatenate continuing lines in single lines
    # and create new list
    listline_new=[]
    nlines=0
    i=0
    while(i<lenlist):
        curr_line=listline[i]
        [nlines,concline]=concatenate(listline,i)
        i=i+nlines+1
        listline_new.append(concline)
    return listline_new

### Reformat lines before writing them to files
def reformat_file(listline):
    lenlist=len(listline)
    # First removes comments and new line
    # Make everything lower case (Fortran is case insensitive)
    for i in range(0,lenlist):
        curr_line=listline[i]
        if curr_line.find("\n"):
        	curr_line=curr_line+"\n"
        curr_line=curr_line.replace("#do not parse","")
        listline[i]=curr_line

    return listline



### Remove comments from current line
def rm_comments(line):
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

### Remove new line symbols
def rm_newline(line):
    linecopy=line.replace("\n","")
    return linecopy

### Remove special characters
def rm_spechar(line):
    idecl=line.find("::")
    if (idecl>=0):
        linecopy=line.replace("_num","")
        linecopy=linecopy.replace("(num)","(kind=8)")
        #linecopy=linecopy.replace("allocatable","pointer")
    else:
        linecopy=line
    return linecopy


### Concatenate continuing line and return
### (i) number of concatenated lines
### (ii) the concatenated line string
def concatenate(listline,istart):
    concline=""
    nlines=0
    i=istart
    curr_line=listline[i].strip()
    curr_line=rm_comments(curr_line)
    endofcont=False
    # If not a precompiler directve, since they can contain &&
    if (not(curr_line.find("#if")>=0))and(not(curr_line.find("#elif")>=0))and(not(curr_line.find("#else")>=0)):
      while (not endofcont):
          iesper=curr_line.find("&")
          iomp=curr_line.find("!$omp")
          if ((iesper>= 0) & (iomp==-1)):
            concline=concline+curr_line[0:iesper].strip()
            i=i+1
            nlines=nlines+1
            curr_line=listline[i]
            curr_line=rm_comments(curr_line)
          else:
            concline=concline+curr_line.strip()
            endofcont=True
    else:
      nlines=0
      concline = curr_line
    return [nlines, concline]

# Get type of variables
def get_type(line):
    typeargs=line.replace(" ","").split(",")
    typechain=""
    addname=""
    # First write dimensions if array
    for i in range(0,len(typeargs)):
        curr_arg=typeargs[i]
        if (curr_arg.find("dimension")>=0):
            # In this case adds dimensions to all varnames
            ideb=line.find("dimension(")
            addname=addname+line[ideb+len("dimension(")-1:len(line)].split(")")[0]+")"
        if (curr_arg.find("type(")>=0):
            ideb=curr_arg.find("(")
            ifin=curr_arg.find(")")
            typechain= typechain+curr_arg[ideb+1:ifin]
            #break
        if (curr_arg.find("integer")>=0):
            typechain=typechain+"integer"
        if (curr_arg.find("real")>=0):
            typechain=typechain+"real"
        if (curr_arg.find("complex")>=0):
            typechain=typechain+"complex"
        if (curr_arg.find("pointer")>=0):
            typechain= "_"+typechain
        if (curr_arg.find("target")>=0):
            typechain= "_istarget_"+typechain
        if (curr_arg.find("allocatable")>=0):
            typechain= "_"+typechain
        if (curr_arg.find("logical")>=0):
            typechain= typechain+"logical"
        if (curr_arg.find("character")>=0):
            typechain= typechain+"string"
    return [addname,typechain]

# Parses Fortran 90 declaration line
def parse_decl_line(line):
    idecl=line.find("::")
    argsname=[]
    argstype=[]
    initvalue=[]
    # Parse types
    [addname,currtype]=get_type(line[0:idecl])
    # Parse var names and init values
    varchain=line[idecl+2:len(line)]
    icomma=find_all(varchain,",")
    ileftp=find_all(varchain,"(")
    irightp=find_all(varchain,")")
    # Finf indexes of variables without comma in the parenthesis (i.e arrays)
    icomnew=icomma[:]
    for icom in range(0,len(icomma)):
        for ip in range(0,len(ileftp)):
            if ((icomma[icom]<=irightp[ip]) & (icomma[icom]>=ileftp[ip])):
                icomnew.pop(icomnew.index(icomma[icom]))

    # Get vars and init values
    ist=0
    for i in range(0,len(icomnew)+1):
        # Still in var list
        if (i<len(icomnew)):
            currval=""
            currvar=varchain[ist:icomnew[i]]
            ieq=currvar.find("=")
            if (ieq>=0):
                currval=currvar[ieq+1:len(currvar)]
                currvar=currvar[0:ieq]
            ipar=currvar.find("(")
            # Conflict array def- Keep variable side def
            if ((ipar>=0) & (addname != "")):
                addnamec=""
            else:
                addnamec=addname
            argsname.append((currvar+addnamec).replace(" ",""))
            argstype.append(currtype.replace(" ",""))
            initvalue.append(currval.replace(" ",""))
            ist=icomnew[i]+1
        # end of var list
        else:
            currval=""
            currvar=varchain[ist:len(varchain)]
            ieq=currvar.find("=")
            if (ieq>=0):
                currval=currvar[ieq+1:len(currvar)]
                currvar=currvar[0:ieq]
            ipar=currvar.find("(")
            # Conflict array def- Keep variable side def
            if ((ipar>=0) & (addname != "")):
                addnamec=""
            else:
                addnamec=addname
            argsname.append((currvar+addnamec).replace(" ",""))
            argstype.append(currtype.replace(" ",""))
            initvalue.append(currval.replace(" ",""))

    return [argsname, argstype, initvalue]

def find_all(line,substr):
    iocc=[]
    istart=0
    while True:
        ieq=line.find(substr)
        if (ieq==-1):
            break
        else:
            iocc.append(istart+ieq)
            line=line[ieq+1:len(line)]
            istart=istart+ieq+1
    return iocc


def fuse_subroutine_module_files(appname):

    ##### Open new file to write modules
    # Subroutine file
    fws=open(appname+"_subroutines.F90","r")
    listline_s=fws.readlines()
    # Module file
    fwm=open(appname+"_modules.F90","r")
    listline_m=fwm.readlines()

    # remove # parsing flags before compilation
    for i in range(0,len(listline_s)):
        listline_s[i]=listline_s[i].replace("#do not parse","")
    for i in range(0,len(listline_m)):
        listline_m[i]=listline_m[i].replace("#do not parse","")

    # Write all lines in a single file
    fwrite=open(appname+".F90","w")
    fwrite.writelines(listline_m)
    fwrite.writelines(listline_s)

    # Delete old files
    fws.close()
    fwm.close()
    #os.remove(appname+"_subroutines.F90")
    #os.remove(appname+"_modules.F90")
    fwrite.close()

def remove_file(filename):
    try:
            os.remove(filename)
    except OSError:
            pass
#### TO DO LIST
#- If subroutine arguments is an array of types

#### ADD in doc
# After subroutine alwas use end subroutine. This is also tru for all other declaratio tokens (type, module etc.)
# Derived types have to be placed in their own modules etc.


###### ----------------------------------------------------
###### SMALL TEST
###### ----------------------------------------------------

appname="picsar"
workdir="src/"
os.chdir(workdir)
remove_file(appname+".F90")
remove_file(appname+".v")

#LIST ALL .F90 or .F files in current directory
listfiles=["modules/modules.F90", \
           "housekeeping/sorting.F90", \
           "field_solvers/Maxwell/maxwell_solver_manager.F90", \
           "field_solvers/Maxwell/yee_solver/yee.F90", \
           "field_solvers/Maxwell/karkkainen_solver/karkkainen.F90", \
           "field_solvers/Maxwell/GPSTD_solver/GPSTD.F90", \
           "parallelization/tiling/tiling.F90", \
           "particle_pushers/boris_pusher/boris_2d.F90", \
           "particle_pushers/boris_pusher/boris_3d.F90", \
           "particle_pushers/vay_pusher/vay_3d.F90", \
           "particle_pushers/laser_pusher_manager_3d.F90", \
           "particle_pushers/particle_pusher_manager_2d.F90", \
           "particle_pushers/particle_pusher_manager_3d.F90", \
           "particle_deposition/current_deposition/direct/direct_current_deposition_2d.F90", \
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

def create_listfiles(folder):
    listfiles = []
    for root, subFolder, files in os.walk(folder):
        for file in files:
           formatted_file = '%s/%s'%(root, file)
           listfiles.append(formatted_file[2:])
    return listfiles

listfiles_miniapp = create_listfiles('.')
print listfiles_miniapp

new_listfile = []
for files in listfiles:
    if files in listfiles_miniapp:
       new_listfile.append(files)

listfiles = new_listfile
print listfiles

# Pre-parse all application files in two .F90 files
# appname_subroutines.F90 and appnam_modules.F90
fortran_preparser(listfiles,appname)

# Parse appname_subroutines.F90 and modules.F90
# Create appname.v file
parse_forthon_app(appname)

# Fuse appname_subroutines.F90 and modules.F90
# in appname.F90
#fuse_subroutine_module_files(appname)

# Compile app with forthon compiler

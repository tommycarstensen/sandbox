#!/bin/env python
#!/usr/bin/python2
# Script for finding complexes where the ligand is close to crystal contacts
#

import sys
sys.path.append('/local/chresten/pKa/trunk/pKarun')

import WHAT_IF
import os
import sys
import math
from read_write_files import *


def generate_crystal_contacts_shell_old(id):

    command='getmol '+id+'\nY\n\n\n %SOUSHL\n %makmol\n\n cryst_'+id+'.pdb\n\n tot 0\n\n'
    
    x=WHAT_IF.Whatif()
        
    x.RedirectOutput('cryst_%s.log'%id)
    x.Run(command,'log_%s.log'%id)

    return

def generate_crystal_contacts_shell(id):

    command='setwif 1052 0\n\ngetmol '+id+'\nY\n\n\n %delwat\n %SOUSHL\n \n %makmol\n\n cryst_'+id+'.pdb\n\n tot 0\n\n'
    
    x=WHAT_IF.Whatif()


    x.RedirectOutput('cryst_%s.log'%id)
    x.Run(command,'log_%s.log'%id)

    return                                      # the above two commands will be used to obtain...
                                                # the orig. pdb file + crystal contacts complex

def generate_crystal_contacts_shell_from_file(filename):

    command='getmol '+filename+'\n\n\n %delwat\n %dellig\n %SOUSHL\n \n %makmol\n\n cryst_'+filename[-8:-4]+'.pdb\n\n tot 0\n\n'
    print command
    x=WHAT_IF.Whatif()
    
    
    x.RedirectOutput('cryst_%s.log'%filename[-8:-4])
    x.Run(command,'log_%s.log'%filename[-8:-4])

    return                                      # the above two commands will be used to obtain...
                                                # the orig. pdb file + crystal contacts complex


	



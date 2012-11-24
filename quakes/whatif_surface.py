#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2007

##import sys
##sys.path.append('/software/Python-2.5/lib/python2.5/site-packages/')
##import MMTK
##from MMTK.MolecularSurface import surfaceAndVolume
##from MMTK.MolecularSurface import surfaceAtoms
##from MMTK.PDB import PDBConfiguration
##
##
##configuration = PDBConfiguration('2lzm.pdb')
##molecules = configuration.createAll()
##surf = surfaceAtoms(molecules)
##print surf
####sv = surfaceAndVolume(

import os

fd = open('script.txt','w')
fd.writelines([
    '/software/whatif_debugged/DO_WHATIF.COM <<EOF\n',
    'GETMOL 1atp.pdb\n',
    'setname\n',
    '%DELWAT\n',
    ##'WATRAD=1.4\n',
    ##'ACCTYP=1\n',
    '%SETACC\n',
    'NOWAT 0\n',
    'NOWAT 0\n',
    '%LISTA\n',
    'TOT 0\n',
    'STOP\n',
    'Y\n',
    ])
fd.close()

os.system('source script.txt > tmp.txt')

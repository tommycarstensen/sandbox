#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2009

## script for generating crystal contacts of biological units

import os, numpy, math
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import biounit,parse_pdb
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import combinatorics

path_pdb = '/data/pdb-v3.2'

d_radii_vdw = {
    'H':1.20,
    'C':1.70,
    'N':1.55,
    'O':1.52,
    'S':1.80,
    }

l_translations = combinatorics.permutation_w_rep([-1,0,1,],3)
l_translations.remove([0,0,0,])

s_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' # 0123456789

## Give the cutoff for "NEAR" contacts (Angstrom)
if '-sympar' in sys.argv:
    sympar = float(sys.argv[sys.argv.index('-sympar')+1])
else:
    sympar = 5.

def main(pdb):

    if not os.path.isfile('%s.pdb' %(pdb)):
        os.system('cp %s/%s/pdb%s.ent %s.pdb' %(path_pdb,pdb[1:3],pdb,pdb,))

    ##
    ## whatif crystal contacts (for comparison to Tommy method)
    ##
    source = 'whatif.src'
    fd = open(source,'w')
    fd.writelines([
        '/local/software/whatif/DO_WHATIF.COM <<EOF\n',
        'GETMOL %s.pdb\n' %(pdb,),
        '%s\n' %(pdb,),
        '%DELWAT\n',
        '%DELLIG\n',
        '%SYMPAR\n',
        '%f\n' %(sympar),
        '%SOUSHL\n', ## crystal contacts
        '%MAKMOL\n',
        '\n', ## The file header will be copied from a PDB file. Hit return for the default header that has no information in it.
        '%s_soushl.pdb\n' %(pdb),
        'TOT 0\n',
        '\n', ## REMARKS
        'STOP\n',
        'Y\n',
        ])
    fd.close()
    if os.path.isfile('%s_soushl.pdb' %(pdb)):
        os.remove('%s_soushl.pdb' %(pdb))
##    os.system('source %s > whatif_surface/%s.out' %(source, pdb))
    os.system('source %s/%s > %s.out' %(os.getcwd(), source, pdb,))

    os.system('rm DRG* DAVADRUG.PDB ALTERR.LOG PDBFILE FOR*.DAT WHATIF.FIG')
    if os.path.isfile('%s_soushl.pdb' %(pdb)):
        os.system('rm whatif.src')
    os.system('rm %s.out' %(pdb))
##    stop_dont_proceed_with_tommys_method
##    os.system('rm %s.pdb' %(pdb))

    return


if __name__ == '__main__':

    pdb = sys.argv[sys.argv.index('-pdb')+1]
    main(pdb)

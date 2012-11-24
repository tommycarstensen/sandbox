import os
import sys

path_pdb = '/data/pdb-v3.2'

pdb = sys.argv[sys.argv.index('-pdb')+1]

def main(pdb,):

    if not os.path.isfile('%s.pdb' %(pdb)):
        os.system('cp %s/%s/pdb%s.ent %s.pdb' %(path_pdb,pdb[1:3],pdb,pdb,))

    ##
    ## calculate surface
    ## 
    lines = [
##        '/software/whatif/DO_WHATIF.COM <<EOF\n',
        '/local/software/whatif/DO_WHATIF.COM <<EOF\n',
        'GETMOL %s.pdb\n' %(pdb),
        '%s\n' %(pdb),
        '%DELWAT\n',
        '%DELLIG\n',
        ##'WATRAD=1.4\n',
        ##'ACCTYP=1\n',
        '%SETACC\n',
        'NOWAT 0\n',
        'NOWAT 0\n',
        '%LISTA\n',
        'TOT 0\n',
        'STOP\n',
        'Y\n',
        ]

    source = 'WHATIF_SETACC.src'
        
    fd = open(source,'w')
    fd.writelines(lines)
    fd.close()

    print 'source %s/%s > %s_SETACC.out' %(os.getcwd(),source, pdb,)
    os.system('source %s/%s > %s_SETACC.out' %(os.getcwd(),source, pdb,))

    ##
    ## clean up the mess that whatif left
    ##
    l_fn = ['ALTERR.LOG','PDBFILE','DAVADRUG.PDB','WHATIF.FIG',]
    for fn in l_fn:
        if os.path.isfile(fn):
            os.remove(fn)
    l = os.listdir(os.getcwd())
    for s in l:
        if s[:3] == 'DRG' and s[6] == '.':
            os.remove(s)
    os.remove(source)

    return

if __name__ == '__main__':
    main(pdb)

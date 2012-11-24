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

def main(pdb, fn_in, fn_out,):

    run_whatif(pdb, fn_in, fn_out,)

    d_acc = parse_output(fn_out,)

    return d_acc


def run_whatif(pdb,fn_in,fn_out,bool_recalc=True,):

    if os.path.isfile(fn_out,):
        if bool_recalc == False:
            return
        elif bool_recalc == True:
            os.remove(fn_out,)
        

    ##
    ## calculate surface
    ## 
    lines = [
##        '/software/whatif/DO_WHATIF.COM <<EOF\n',
        '/local/software/whatif/DO_WHATIF.COM <<EOF\n',
        'GETMOL %s\n' %(fn_in),
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

    print 'source %s/%s > %s' %(os.getcwd(), source, fn_out,)
    os.system('source /local/tc/quakes/%s > %s' %(source, fn_out,))

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

    return


def parse_output(fn_out,):

    ##
    ## parse output
    ##
    d_acc = {}

    fd = open(fn_out,'r')
    lines = fd.readlines()
    fd.close()

    for i in range(len(lines)):
        line = lines[i]
        if line[:37] == 'New accessibilities calculated ... : ':
            k = i
            ## loop until residue
            for j in range(i+2,len(lines)):
                line = lines[j]
                if 'Residue:' in line:
                    index = line.index('Residue:')
                    chain = lines[j][index+29]
                    if lines[j][index+21:index+25].strip() == 'OXT':
                        res_no = 'OXT'
                    else:
                        res_no = int(lines[j][index+21:index+25])
                    iCode = lines[j][index+26]
                    if not chain in d_acc.keys():
                        d_acc[chain] = {}
                    if not res_no in d_acc[chain].keys():
                        d_acc[chain][res_no] = {}
                    d_acc[chain][res_no][iCode] = {}
                    res_name = lines[j][15:18]
                    if res_name == 'HOH':
                        continue
##                    if lines[j][index+21:index+25].strip() == 'OXT':
##                        continue
                    ## loop until column headers of residue
                    for k in range(j+1,len(lines)):
                        line = lines[k]
                        if 'Atom     X     Y     Z   Acc   B   WT   VdW  Colr   AtOK  Val' in line:
                            ## loop over data lines
                            for l in range(k+1,len(lines)):
                                line = lines[l]
                                if line == ' \n':
                                    continue
                                if 'Residue:' in line:
                                    break
                                if 'Option not found' in line:
                                    break
                                if 'Which in turn is attached to' in line:
                                    break
                                atom_name = line[:4].strip()
                                if atom_name == "O'":
                                    atom_name = 'O'
                                if line[24:28].strip() == '---':
                                    continue
##                                print pdb
##                                print 'i', lines[i]
##                                print 'j', lines[j]
##                                print 'k', lines[k]
##                                print 'l', lines[l]
                                acc = float(line[24:28])
                                if atom_name in d_acc[chain][res_no][iCode].keys():
                                    print chain,res_no,iCode,atom_name
                                    stop
                                d_acc[chain][res_no][iCode][atom_name] = acc
                            break ## Atom in line
                    continue ## Residue: in line
            break ## new accesibilities calculated

    return d_acc

if __name__ == '__main__':
    main()

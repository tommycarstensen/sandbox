## built-in
import urllib2
## add-on
import numpy
## my own sandbox...
import parse_mmCIF

path_pdb = 'pathtoyourpdbs'
path_pdb = 'tmp'
l_pdbs = ['listofyourpdbs']
import os
l_fns = os.listdir(path_pdb)
l_pdbs = [fn[:4] for fn in l_fns]


def main():

    d_atoms = identify_CH_bonds()

    loop_pdbs(l_pdbs,d_atoms,)

    return


def loop_pdbs(l_pdbs,d_atoms,):

    for pdb in l_pdbs:

        fp = '%s.pdb' %(os.path.join(path_pdb,pdb,))
        fd = open(fp,'r')
        lines = fd.readlines()
        fd.close()

        for line in lines:
            record = line[:6].strip()
            ## also skip MODRES...
            if record != 'ATOM':
                continue
            res_name = line[17:20].strip()
##            ## skip if not residue with pi system
##            if res_name not in [
##                ## Pi system
##                'TYR','PHE','HIS','TRP',
##                ]:
##                continue
            atom_name = line[12:16].strip()

            ## skip if not C-H bond
            bool_CH = False
            if atom_name in d_atoms[res_name]:
                bool_CH = True
            if (
                (res_name == 'TYR' and atom_name in ['CD1','CD2','CE1','CE2',])
                or
                (res_name == 'PHE' and atom_name in ['CD1','CD2','CE1','CE2',])
                or
                (res_name == 'HIS' and atom_name in ['CG','CD2','ND1','CE1','NE2',])
                or
                (res_name == 'TRP' and atom_name in ['CD2','CE2','CE3','CZ2','CZ3','CH2',])
                ):
                res_name not in ['TYR','PHE','HIS','TRP',]
                ):
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            coord = numpy.array([x,y,z,])
            print coord
            stop

    return


def identify_CH_bonds():
    
    ##
    ## identify all C-H single bonds in the standard residues
    ##
    d_atoms = {}
    for residue in [
        'ALA',
##        'ALA','CYS','ASP','GLU','PHE',
##        'GLY','HIS','ILE','LYS','LEU',
##        'MET','ASN','PRO','GLN','ARG',
##        'SER','THR','VAL','TRP','TYR',
        ]:
        lines = urllib2.urlopen('http://www.pdb.org/pdb/files/ligand/%s.cif' %(residue)).readlines()
        d = parse_mmCIF.main(residue,lines)
        d_atoms[residue] = []
        for i in range(len(d['_chem_comp_bond.comp_id'])):
            if d['_chem_comp_bond.value_order'][i] != 'SING':
                continue
            atom1 = d['_chem_comp_bond.atom_id_1'][i]
            atom2 = d['_chem_comp_bond.atom_id_2'][i]
            ## heavy element is always listed before hydrogen
            if atom1[0] != 'C' or atom2[0] != 'H':
                continue
            print residue, d['_chem_comp_bond.atom_id_1'][i], d['_chem_comp_bond.atom_id_2'][i]
            d_atoms[residue] += [atom1]

    return d_atoms


if __name__ == '__main__':
    main()

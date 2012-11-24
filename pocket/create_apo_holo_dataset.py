## 1) parse proteins with motions from molmovdb
## 2) identify those proteins that exist in an apo and holo form (few?) by het records
## 3a) calc motion between pdbs
## 3b) identify binding site and see if do better when motion described by nma mode 7

## OR

## identify binding site residues and not active site residues
## require that PDBs are from the same space group so observed motion (structural difference) is not due to space group differences

## do succes rate of bind site id for brackets of apo-holo-overlaps (0.8-0.9 etc.)

## i want a big dataset, but i also have to run it through the other programs...

## exclude 95% identity to reduce size of dataset but only compare 100% identity

import urllib2, os, numpy, math, copy
import sys
sys.path.append('C:/Users/Tommy Carstensen/Documents/svn/tc_sandbox/pdb')
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('C:/Users/Tommy Carstensen/Documents/svn/tc_sandbox/quakes')
sys.path.append('/home/tc/svn/tc_sandbox/quakes')
import smallmolecules
sys.path.append('/home/tc/svn/GoodVibes')
import NMA
sys.path.append('/home/tc/svn/Protool/')
import geometry

path_mmCIF = '/media/Elements/mmCIF'
path_mmCIF = '/home/tc/Downloads/mmCIF'
path_mmCIF = '/media/Tommy/mmCIF'
path_mmCIF = '/media/WDMyBook1TB/2TB/mmCIF'

def main():

##    d_apo2holo = {}
##    set_pdbs_apo = set()
##    for cutoff in [
####        6,8,
##        10,
##        ]:
##        fd = open('max_cutoff%i.txt' %(cutoff), 'r')
##        lines = fd.readlines()
##        fd.close()
##        l_overlaps = []
##        for line in lines:
##            index = line.index(' ')
##            cluster = line[:index]
##            l = eval(line[index+1:-1])
##            if len(l) == 6:
##                [overlap,overlap_mode7,pdb_holo,pdb_apo,mode,rmsd,] = l
##            elif len(l) == 5:
##                [overlap,pdb_holo,pdb_apo,mode,rmsd,] = l
##            elif len(l) == 7:
##                [overlap,overlap_mode7,pdb_holo,pdb_apo,mode,rmsd,ligand,] = l
##            set_pdbs_apo |= set([pdb_apo])
##            d_apo2holo[pdb_apo] = {
##                'holo':pdb_holo,'ligand':ligand,
##                'overlap':overlap_mode7,'rmsd':rmsd,
##                }
##    print set_pdbs_apo
##    print d_apo2holo
##    print len(set_pdbs_apo)
##    print len(d_apo2holo.keys())
##    stop

    print 'identify structures with large ligands'
    l_pdbs_ligands_large, d_pdbs_ligands_large, l_pdbs_ligands = find_holo()

######    print '1j8a' in l_pdbs_ligands
######
######    d = parse_mmCIF.main('1j8a',)
######    l_nonpolymers_large, l_saccharides, l_nonpolymers_all = identify_ligands(
######        '1j8a',d,
######        )
######    print l_nonpolymers_all
######
######    l_include_additional = ['_entity.type__non-polymer','_one_polysaccharide',]
######    set_pdbs = combine_sets(l_include_additional=l_include_additional)
######    l_pdbs = list(set_pdbs)
######    l_pdbs.sort()
######    print '1j8a' in l_pdbs
######
######    stop

##    ## print statistics
##    print 'pdbs with large ligands', len(l_pdbs_ligands_large)
##    stop

    print 'find apo structures from clusters containing holo structures'
    d_clusters = find_apo(l_pdbs_ligands_large,l_pdbs_ligands,)
##    ## print some statistics...
##    print 'clusters', len(d_clusters.keys())
##    apo = 0
##    holo = 0
##    for cluster in d_clusters.keys():
##        apo += len(d_clusters[cluster]['apo'])
##        holo += len(d_clusters[cluster]['holo'])
##    print 'apo and holo in clusters', apo, holo
##    stop_tmp_statistics

##    find_different_ligand_binding_sites(d_clusters,d_pdbs_ligands_large,)

##    d_clusters = find_overlap_with_other_DB(d_clusters)
####    ## print some statistics...
####    print 'clusters', len(d_clusters.keys())
####    apo = 0
####    holo = 0
####    for cluster in d_clusters.keys():
####        apo += len(d_clusters[cluster]['apo'])
####        holo += len(d_clusters[cluster]['holo'])
####    print 'apo and holo in clusters', apo, holo
####    stop_tmp_statistics

    ## calculations carried out on apo structures
    ## therefore find apo structure which best describes apo-to-holo motion by mode 7
    l_pdbs_apo_best = identify_best_apo_structure(d_clusters, d_pdbs_ligands_large)
    print l_pdbs_apo_best
    stop

    identify_sequence_similar_apo_pdbs(l_pdbs_apo_best)

    return


def parse_csa():

    ## read CSA database
    fd = open('/home/tc/UCD/GV_ligand_binding_site_identification/datasets/CSA_2_2_12.dat','r')
    lines = fd.readlines()
    fd.close()

    d_catres = {}
    for i_line in range(1,len(lines)):
        if i_line % 10000 == 0:
            print 'csa', i_line, len(lines)
        line = lines[i_line]
        l = line.strip().split(',')

        pdbID = l[0]
        site_number = int(l[1])
        res_name = l[2] ## residue type
        chainID = l[3]
##        res_no = int(l[4])
        res_no = l[4] ## use string for set comparisons with mmCIF
##        chemical_function = l[5]
##        evidence_type = l[6]
##        refID = l[7] ## literature entry

        ## ion
        if chainID == '':
            continue
        if res_name == '':
            continue
        if site_number > 0:
            continue

        if not pdbID in d_catres.keys():
            d_catres[pdbID] = {
                'chains':[],'res_names':[],'res_nos':[],
##                'refs':[],
                }
        d_catres[pdbID]['chains'] += [chainID]
        d_catres[pdbID]['res_names'] += [res_name]
        d_catres[pdbID]['res_nos'] += [res_no]
##        d_catres[pdbID]['refs'] += [refID]

    return d_catres


def find_different_ligand_binding_sites(d_clusters,d_pdbs_ligands_large,):

    print 'what does this function do?'
    ## it was probably meant to check if the ligand binding site is identical across holo structures...

    for cluster in d_clusters.keys():
        set_pdbs_holo = d_clusters[cluster]['holo']
        l_pdbs_holo = list(set_pdbs_holo)
        for i1 in range(len(l_pdbs_holo)-1):
            pdb_holo1 = l_pdbs_holo[i1]
            ligand1 = d_pdbs_ligands_large[pdb_holo1]
            coord1 = find_ligand_center(pdb_holo1,ligand1,)
            for i2 in range(i1+1,len(l_pdbs_holo)):
                pdb_holo2 = l_pdbs_holo[i2]
                ligand2 = d_pdbs_ligands_large[pdb_holo2]
                coord2 = find_ligand_center(pdb_holo2,ligand2,)
##                print math.sqrt(sum((coord2-coord1)**2))
                print cluster, pdb_holo1, pdb_holo2, ligand1, ligand2
                if [ligand1,ligand2].count('saccharide') == 1:
                    print ligand1,ligand2
                    print pdb_holo1, pdb_holo2
                    stop
##                if ligand1 != ligand2:
##                    print ligand1, ligand2
##                    print pdb_holo1, pdb_holo2
##                    stop

    stop_end           
    return


def find_ligand_center(pdb_holo,ligand,):

    d_mmCIF_holo = parse_mmCIF.main(pdb_holo,)

    l_coords_ligand = []
    for i in range(len(d_mmCIF_holo['_atom_site.id'])):
        if not d_mmCIF_holo['_atom_site.group_PDB'][i] == 'HETATM':
            continue
        if d_mmCIF_holo['_atom_site.label_comp_id'][i] == 'HOH':
            continue
        if d_mmCIF_holo['_atom_site.label_comp_id'][i] not in ligand:
            continue
        x = float(d_mmCIF_holo['_atom_site.Cartn_x'][i])
        y = float(d_mmCIF_holo['_atom_site.Cartn_y'][i])
        z = float(d_mmCIF_holo['_atom_site.Cartn_z'][i])
        coord = numpy.array([x,y,z,])
        l_coords_ligand += [coord]

    position_ligand = sum(l_coords_ligand)/len(l_coords_ligand)

    return position_ligand


def find_overlap_with_other_DB(d_clusters):

    set_pdbs = set()
    for cluster in d_clusters.keys():
        set_pdbs |= d_clusters[cluster]['holo']
        set_pdbs |= d_clusters[cluster]['apo']

    d_sets = {}
    l_fns = [
        'ligasite.csv', ## redundant, *apo* or holo, organized
##        'CSA_2_2_12.dat', ## redundancy, *apo* and holo, mixed
##        'BindingMOAD_data_for_PDB.csv', ## *holo* only, Ki values
        ]

    for fn in l_fns:

        fd = open('datasets/%s' %(fn),'r')
        lines = fd.readlines()
        fd.close()

        l_pdbs = []
        for line in lines[1:]:
            l_pdbs += [line[:4].lower()]

        print fn, len( set(l_pdbs) & set_pdbs ), len(set(l_pdbs)), len(set_pdbs)
        set_overlap = set(l_pdbs) & set_pdbs
        d_sets[fn] = set_overlap
        count = 0
        for cluster in d_clusters.keys():
            if (
                len(set_overlap & d_clusters[cluster]['apo']) > 0
                or
                len(set_overlap & d_clusters[cluster]['holo']) > 0
                ):
                count += 1
            ## not present and only one db compared (then manipulate d_clusters)
            elif len(l_fns) == 1:
                del d_clusters[cluster]
        print 'present in', count, 'clusters'

    return d_clusters


def count_heavy_atoms_in_formula(formula):

    l_atoms = formula.split()
    count_atoms_heavy = 0
    for atoms in l_atoms:

        ## charge (not atom)
        if l_atoms.index(atoms) == len(l_atoms)-1 and atoms[0] in ['-','1','2','3','4','5','6',]:
            continue

        if l_atoms.index(atoms) == len(l_atoms)-1 and atoms[0] in ['+','7','8','9','0',]:
            print formula
            print atoms
            stop

        ## parse element
        element = atoms
        for s in '0123456789':
            element = element.replace(s,'')

        ## only count heavy atoms
        if element == 'H':
            continue

        for s in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            atoms = atoms.replace(s,'')

        ## only one atom of that element type
        if len(atoms) == 0:
            count_atoms_heavy += 1
        ## multiple atoms of that element type
        else:
            count_atoms_heavy += int(atoms)

    return count_atoms_heavy


def identify_ligands(pdb,d):

    l_nonpolymers_large = []
    l_saccharides = []
    l_nonpolymers_all = []

    d_smallmolecules = smallmolecules.main()

    for i in range(len(d['_chem_comp.type'])):

        if d['_chem_comp.id'][i] in ['HOH','UNL',]:
            continue

        if d['_chem_comp.type'][i].lower() == 'non-polymer':

            hetID = d['_chem_comp.id'][i]
            bool_small = False
            for k in ['ions','metals','clusters','solutes','functional',]:
                if hetID in d_smallmolecules[k]:
                    bool_small = True
                    break
            if bool_small == True:
                continue

            l_nonpolymers_all += [hetID]

            formula = d['_chem_comp.formula'][i]
            count_atoms_heavy = count_heavy_atoms_in_formula(formula)
            if count_atoms_heavy < 10:
                continue

            l_nonpolymers_large += [hetID]

        elif 'saccharide' in d['_chem_comp.type'][i].lower():
            hetID = d['_chem_comp.id'][i]
            l_saccharides += [hetID]

    return l_nonpolymers_large, l_saccharides, l_nonpolymers_all


def find_holo():

    print 'identify structures with large ligands'

##    if os.path.isfile('pdbIDs_large_ligands.out'):
    if os.path.isfile('pdbIDs_ligands.out'):

        fd = open('pdbIDs_large_ligands.out')
        s = fd.read()
        fd.close()
        l_pdbs_ligands_large = s.split(',')

        fd = open('pdbIDs_large_ligands.dict')
        s = fd.read()
        fd.close()
        d_pdbs_ligands_large = eval(s)

        fd = open('pdbIDs_ligands.out')
        s = fd.read()
        fd.close()
        l_pdbs_ligands = s.split(',')

        return l_pdbs_ligands_large, d_pdbs_ligands_large, l_pdbs_ligands

    d_csa = parse_csa()

    ## get list of small molecules to exclude
    d_smallmolecules = smallmolecules.main()

    ## include if ligand (also mono- and poly-saccharides)
    l_include_additional = ['_entity.type__non-polymer','_one_polysaccharide',]
    set_pdbs = combine_sets(l_include_additional=l_include_additional)
    l_pdbs = list(set_pdbs)
    l_pdbs.sort()

    fd = open('%s/list_one_polysaccharide.txt' %(path_mmCIF,),'r')
    s = fd.read() ## tmp!!!
    fd.close()
    set_one_polysaccharide = set(s.lower().split('\n'))

    ##
    ## find pdbs with large ligands
    ##
    l_pdbs_ligands_large = []
    l_pdbs_ligands = []
    d_pdbs_ligands_large = {}
    for i_pdb in range(len(l_pdbs)):

        if i_pdb % 100 == 0:
            print i_pdb, len(l_pdbs)

        pdb = l_pdbs[i_pdb]

##        if not pdb in set_polysaccharide:
##            continue

        d = parse_mmCIF.main(pdb,)
        l_nonpolymers_large, l_saccharides, l_nonpolymers_all = identify_ligands(
            pdb,d,
            )

        ## pdb with one or more ligands?
        if len(l_nonpolymers_all) > 0:
            l_pdbs_ligands += [pdb]

        ## more than one ligand?
        if (
            len(l_nonpolymers_large) > 1
            or
            (
                len(l_nonpolymers_large) > 0
                and
                len(l_saccharides) > 0
                )
            ):
##            if pdb in set_polysaccharide-set(['2zhn','2z8l','2wm0','2w68',]):
##                fd = open('remediation_nonpolymer_to_saccharide.txt','a')
##                fd.write('%s %s %s\n' %(pdb,l_nonpolymers,l_saccharides,))
##                fd.close()
            continue
        if len(
            ## nonpolymers and saccharides
            set( set(l_nonpolymers_large) | set(l_saccharides) )
            -
            ## not coenzymes or prostethic groups
            set( set(d_smallmolecules['coenzymes']) | set(d_smallmolecules['prosthetic groups']) )
            ) == 0:
            continue

        if len(l_nonpolymers_large) == 1 and len(l_saccharides) == 0:
            hetID = ''.join(l_nonpolymers_large)
            ## no more than 1 molecule! otherwise floating around everywhere! e.g. 2pc2
            ## if hetID not part of larger entity (e.g. 1a85)
##            if hetID in d['_pdbx_entity_nonpoly.comp_id']:
##            _pdbx_nonpoly_scheme.entity_id
##            _pdbx_nonpoly_scheme.mon_id
##            entity_ID = d['_pdbx_entity_nonpoly.entity_id'][d['_pdbx_entity_nonpoly.comp_id'].index(hetID)]
##            ## if more than 1 molecule
##            if d['_entity.pdbx_number_of_molecules'][d['_entity.id'].index(entity_ID)] > 1:
##                continue
            ## more than 1 molecule?
            if d['_pdbx_nonpoly_scheme.mon_id'].count(hetID) > 1:
                continue

        ## only accept one mono- or poly-saccharide
        if len(l_saccharides) > 0 and len(l_nonpolymers_large) == 0:
            if pdb not in set_one_polysaccharide:
                continue

        ## only accept pdb if in catalytic site atlas
        if not pdb in d_csa.keys():
            continue

        ## only accept ligand if vicinal to catalytic site
        ## check that at least one binding residue and one catalytic residue is the same
        d_mmCIF = parse_mmCIF.main(
            pdb,
            l_data_categories = ['_struct_site_gen'],
            l_data_categories_break = ['_struct_site_gen'],
            )

        ## skip if error (_struct_site_gen missing; e.g. 547 in 2of2)
        if not '_struct_site_gen.auth_seq_id' in d_mmCIF.keys():
            fd = open('remediation_struct_site_missing.txt','a')
            fd.write('%s %s\n' %(pdb,str(l_nonpolymers_large),))
            fd.close()
            continue

        if len(
            set(d_csa[pdb]['res_nos'])
            &
            set(d_mmCIF['_struct_site_gen.auth_seq_id']) ## *not* label_seq_id cf. 1a2b
            ) == 0:
            print
            print 'x', d_csa[pdb]['res_nos']
            print 'x', d_mmCIF['_struct_site_gen.auth_seq_id']
            print pdb
            print l_nonpolymers_large
            continue

        l_pdbs_ligands_large += [pdb]
        ## append nonpolymers
        if len(l_nonpolymers_large) == 1 and len(l_saccharides) == 0:
            d_pdbs_ligands_large[pdb] = l_nonpolymers_large
            print pdb, l_pdbs.index(pdb), 'of', len(l_pdbs), 'holo', ''.join(l_nonpolymers_large)
        ## append saccharide
        elif len(l_nonpolymers_large) == 0 and len(l_saccharides) > 0:
            d_pdbs_ligands_large[pdb] = l_saccharides
            print pdb, l_pdbs.index(pdb), 'of', len(l_pdbs), 'holo', 'saccharide'
        ## not expected
        else:
            print pdb
            print l_nonpolymers_large
            print l_saccharides
            stop

    s = ','.join(l_pdbs_ligands_large)
    fd = open('pdbIDs_large_ligands.out','w')
    fd.write(s)
    fd.close()
    fd = open('pdbIDs_large_ligands.dict','w')
    fd.write(str(d_pdbs_ligands_large))
    fd.close()
    fd = open('pdbIDs_ligands.out','w')
    fd.write(str(l_pdbs_ligands))
    fd.close()

    return l_pdbs_ligands_large, d_pdbs_ligands_large, l_pdbs_ligands


def combine_sets(l_include_additional=[],l_exclude_additional=[]):

    l_include_fns = [
        ##
        '_exptl.method__X-RAY_DIFFRACTION',
        ## include if protein
        '_entity_poly.type__polypeptide(L)',
        ## include if monomer
        '_pdbx_struct_assembly.oligomeric_details__monomeric',
        ## include if monomer (ASU)
        '_entity.pdbx_number_of_molecules__1',
        ## and monomer (biounit) ... otherwise which is A and which is B? diff cryst contacts...
        '_one_polypeptide',
        ]

    l_exclude_fns = [
        ## exclude if zero occupancy residues other than at the terminals
        '_pdbx_unobs_residues__NONTERMINAL',
        ## exclude if unobs alpha carbon atoms
        '_pdbx_unobs_atoms__CA',
        ##
##        '_entity_poly.type__polysaccharide(D)',
        '_entity_poly.type__polyribonucleotide','_entity_poly.type__polydeoxyribonucleotide',
        ## risk bound ligands if included
        '_pdbx_struct_mod_residue__notMSE',
        ]

    set_include_pdbs = set()
    for fn in l_include_fns:
        fd = open('%s/list%s.txt' %(path_mmCIF,fn,),'r')
        s = fd.read()
        fd.close()
        if set_include_pdbs == set():
            set_include_pdbs = set(s.lower().split('\n'))
        else:
            set_include_pdbs &= set(s.lower().split('\n'))
    if len(l_include_additional) > 0:
        set_include_additional = set()
        for fn in l_include_additional:
            fd = open('%s/list%s.txt' %(path_mmCIF,fn,),'r')
            s = fd.read()
            fd.close()
            ## either set1 or set2
            set_include_additional |= set(s.lower().split('\n'))
        ## both set1 and set2
        set_include_pdbs &= set_include_additional

    set_exclude_pdbs = set()
    for fn in l_exclude_fns+l_exclude_additional:
        fd = open('%s/list%s.txt' %(path_mmCIF,fn,),'r')
        s = fd.read()
        fd.close()
        set_exclude_pdbs |= set(s.lower().split('\n'))

    set_pdbs = set_include_pdbs-set_exclude_pdbs

    return set_pdbs


def find_apo(l_pdbs_large_ligands,l_pdbs_ligands,):

    print 'find apo'

    if os.path.isfile('find_apo.out'):
        fd = open('find_apo.out')
        s = fd.read()
        fd.close()
        d_clusters = eval(s)
        return d_clusters

    d_clusters = {}

    ## exclude if ligand
    l_exclude_additional = [
##        '_entity.type__non-polymer', ## don't exclude small molecules
##        '_one_polysaccharide',
        ]
    l_pdbs_apo = list(combine_sets(l_exclude_additional=l_exclude_additional))

    fd = open('%s/bc-100.out' %(path_mmCIF),'r')
    lines = fd.readlines()
    fd.close()

    ## loop over clusters
    for i in range(len(lines)):
        if i % 10000 == 0:
            print 'cluster', i, 'of', len(lines)
        line = lines[i]
        l_pdbs = line.lower().split()
        for i_pdb in range(len(l_pdbs)):
            l_pdbs[i_pdb] = l_pdbs[i_pdb][:4]

        ## 1) at least one of the  pdbs with large ligands is present in the cluster
        set_pdbs_holo = set(l_pdbs) & set(l_pdbs_large_ligands)
        if len( set_pdbs_holo ) == 0:
            continue

        set_pdbs_apo = set(l_pdbs) & set(l_pdbs_apo)
        set_pdbs_apo -= set(l_pdbs_large_ligands)
##        if len(set_pdbs_apo) != len(set_pdbs_apo-set(l_pdbs_ligands)):
##            print set_pdbs_apo
##            print set_pdbs_apo&set(l_pdbs_ligands)
##            stop
        set_pdbs_apo -= set(l_pdbs_ligands)
        ## no protein without ligand...
        if len(set_pdbs_apo) == 0:
            continue

##        for pdb_apo in list(set_pdbs_apo):
##            d = parse_mmCIF.main(pdb_apo,)
##            l_nonpolymers_large, l_saccharides, l_nonpolymers_all = identify_ligands(
##                pdb_apo,d,
##                )
##            if len(l_nonpolymers_large) > 0 or len(l_saccharides) > 0:# or len(l_nonpolymers_all) > 0:
##                set_pdbs_apo.remove(pdb_apo)

        d_clusters[i] = {
            'holo':set_pdbs_holo,
            'apo':set_pdbs_apo,
            }

    fd = open('find_apo.out','w')
    fd.write(str(d_clusters))
    fd.close()

    return d_clusters


def identify_best_apo_structure(d_clusters, d_pdbs_ligands_large,):

    print 'identify best apo structure'

##    l_pdbs_apo_best = ['1f10', '1tgn', '1eie', '2cbe', '3k0n', '2hnp', '1uor', '1z1i', '2ppn', '1xqz', '4ape', '1pdb', '1jam', '2wbz', '1bbc', '1gc7', '1a3h', '1rtc', '3lsd', '1wvw', '2e3m', '2jic', '1arl', '3pte', '1yes', '1akz', '3npo', '1alb', '2exo', '1af9', '2zco', '3ly3', '1zg4', '1hfd', '1hka', '1mzl', '1erk', '1ak1', '1thu', '1p38', '1iad', '3a0x', '3ars', '2paw', '1ri5', '1cua', '2qev', '1lmn', '1kx8', '1pnz', '1boi', '2gdn', '3g4p', '1gy0', '1ey0', '1ifc', '1sll', '2gg4', '2sil', '1jcf', '1ks9', '1mla', '3blm', '1znw', '2sga', '1qmt', '2zj8', '2hbj', '1wos', '3mft', '1xix', '1qtr', '2vfb', '1tje', '1bk7', '1i04', '3hny', '1eyh', '1lz4', '1fus', '1rd6', '1cz1', '2o0k', '3aap', '1ptd', '1gcu', '1mtz', '2bv9', '1gqz', '1mri', '1zsa', '1lp8', '3pkv', '1ahc', '3dt8', '1xza', '2vfy', '1lr9', '1t7n', '1nm8', '1gnd', '1sqg', '3g6l', '2ggo', '1eur', '1h09', '2xly', '1tm2', '2ac4', '1yhv', '2plc', '1gnv', '1bue', '1xqo', '3c0e', '1ojq', '1gbs', '3ejf', '3ewq', '1syc', '1sye', '1wka', '1kqx']
##    set_pdbs_apo = set(['1ak1', '1lp8', '2exo', '1syc', '1cua', '1znw', '1gbs', '1w8v', '3a0x', '1sqg', '3npo', '2e3m', '3k0n', '3aap', '1bk7', '1tgn', '1ahc', '3mft', '1erk', '1sll', '1rtc', '1gy0', '2hbj', '1rd6', '1mzl', '1hka', '2ppn', '2vfy', '1yhv', '1eyd', '1eur', '1bbc', '1mri', '2zj8', '2ggo', '2vfb', '1f10', '3c0e', '1tje', '1wvw', '1wos', '1qtr', '1pdb', '2qev', '3g6l', '1arl', '1akz', '2paw', '1kqx', '3pte', '1iad', '2d59', '1vds', '1lmn', '1kf5', '1jcf', '2zco', '2ac4', '1ojq', '1ey0', '1yes', '1sye', '1tqo', '2gg4', '1xqz', '1ri5', '1p38', '1z1i', '4ape', '1ifb', '1jam', '1mtz', '2sil', '1xqo', '3ewq', '2sga', '3blm'])
    l_pdbs_apo_best = []
    l_clusters = d_clusters.keys()
    l_clusters.sort()

    if os.path.isfile('rmsd_v_overlap/all.txt'):
        os.remove('rmsd_v_overlap/all.txt')

    for cluster in l_clusters:

##        if cluster != 2:
##            continue ## tmp!!!

##        if cluster == 2:
##            pdb_apo = 'xxxx'
##            pdb_holo = 'xxxx'
##            l_pdbs_apo_best += [[pdb_apo,pdb_holo,]]
##            continue
        if cluster == 5:
            pdb_apo = '1j8a'
            pdb_holo = '1gi4'
            l_pdbs_apo_best += [[pdb_apo,pdb_holo,]]
            continue
        if cluster == 11:
            pdb_apo = '1rnm'
            pdb_holo = '1rno'
            l_pdbs_apo_best += [[pdb_apo,pdb_holo,]]
            continue
        if cluster == 12:
            pdb_apo = '3hku'
            pdb_holo = '3mzc'
            l_pdbs_apo_best += [[pdb_apo,pdb_holo,]]
            continue
        if cluster == 94:
            pdb_apo = '8est'
            pdb_holo = '1esa'
            l_pdbs_apo_best += [[pdb_apo,pdb_holo,]]
            continue

        if os.path.isfile('rmsd_v_overlap/cluster%i.txt' %(cluster)):
            continue
##            os.remove('rmsd_v_overlap/cluster%i.txt' %(cluster))

        l_pdbs_holo = list(d_clusters[cluster]['holo'])
        l_pdbs_apo = list(d_clusters[cluster]['apo'])
        l_pdbs_holo.sort()
        l_pdbs_apo.sort()

        print 'cluster', cluster, len(l_clusters), 'holo', len(l_pdbs_holo), 'apo', len(l_pdbs_apo)

        max_overlap = [0,'pdb_holo','pdb_apo','mode',]
        for pdb_holo in l_pdbs_holo:

            if pdb_holo in ['1aq7']: ## tmp!!! peptide ligand
                stop_tmp
                continue

            if pdb_holo in [
                '1eou', ## pdbx_poly_seq_scheme.auth_mon_id error
                '1rct','1axb','1tem','1h2f', ## _struct_ref_seq_dif missing
                '4a3h','2wf5', ## incorrect _struct_ref_seq.pdbx_db_accession
                '1n8q','1no3',  ## _struct_ref.pdbx_seq_one_letter_code wrong
                '1b0j', ## _struct_ref_seq_dif wrong
                ]: 
                continue

            if len(l_pdbs_holo) > 5 and len(l_pdbs_apo) > 25:
                print cluster, pdb_holo, l_pdbs_holo.index(pdb_holo), len(l_pdbs_holo)

            d_holo = parse_mmCIF.main(pdb_holo,)
            SG_holo = ''.join(d_holo['_symmetry.space_group_name_H-M'])

            if (
                '(' in ''.join(d_holo['_entity_poly.pdbx_seq_one_letter_code'])
                and
                '_pdbx_struct_mod_residue.id' not in d_holo.keys()
                ):
                seq = ''.join(d_holo['_entity_poly.pdbx_seq_one_letter_code'])
                modres = seq[seq.index('(')+1:seq.index(')')]
                chem_comp_type = d_holo['_chem_comp.type'][d_holo['_chem_comp.id'].index(modres)]
                fd = open('remediation_pdbx_struct_mod_residue.txt','a')
                fd.write('%s %s %s\n' %(pdb_holo,modres,chem_comp_type,))
                fd.close()
                continue

            d_coords, l_coords_alpha_holo = mmCIF2coords.main(pdb_holo, d_holo)

            l_coords_alpha_holo_original = list(l_coords_alpha_holo)

            for pdb_apo in l_pdbs_apo:

                if pdb_holo == '1o3h' and pdb_apo == '1yp9': ## some problem with matrix diagonalization...?
                    continue

                if pdb_apo in [
                    '1ps3', ## seq_dif mmCIF data category missing (LYS/GLU 907 conflict)
                    '1ee3','1ell','1elm','1m4l','1o03','1bza', ## _struct_ref_seq.pdbx_db_accession wrong
                    '2bdx', ## somehow this dimer got included...
                    ]:
                    continue

##                if (
##                    cluster == 5 and (
##                        pdb_apo != '1tgn'
####                        or
####                        pdb_holo != '1c1p'
##                        )
##                    ):
##                    continue ## tmp!!!

                d_apo = parse_mmCIF.main(pdb_apo,)
                SG_apo = ''.join(d_apo['_symmetry.space_group_name_H-M'])

                ## require space groups to be identical so structural differences are due to ligands and not crystal contacts
                if SG_holo != SG_apo:
                    continue

                if (
                    '(' in ''.join(d_apo['_entity_poly.pdbx_seq_one_letter_code'])
                    and
                    '_pdbx_struct_mod_residue.id' not in d_apo.keys()
                    ):
                    seq = ''.join(d_apo['_entity_poly.pdbx_seq_one_letter_code'])
                    modres = seq[seq.index('(')+1:seq.index(')')]
                    chem_comp_type = d_apo['_chem_comp.type'][d_apo['_chem_comp.id'].index(modres)]
                    fd = open('remediation_pdbx_struct_mod_residue.txt','a')
                    fd.write('%s %s %s\n' %(pdb_apo,modres,chem_comp_type,))
                    fd.close()
                    continue

                d_coords, l_coords_alpha_apo = mmCIF2coords.main(pdb_apo, d_apo)

                l_coords_alpha_holo = list(l_coords_alpha_holo_original)

                l_coords_alpha_apo, l_coords_alpha_holo = sequential_alignment_of_coordinates(
                    l_coords_alpha_apo,l_coords_alpha_holo,
                    d_apo,d_holo,
                    pdb_apo, pdb_holo,
                    )

                ##
                ## align apo and holo structure
                ##
                instance_geometry = geometry.geometry()
                rmsd = instance_geometry.superpose(l_coords_alpha_holo,l_coords_alpha_apo)
                tv1 = instance_geometry.fitcenter
                rm = instance_geometry.rotation
                tv2 = instance_geometry.refcenter
##                if rmsd > 2.5:
##                    print rmsd, pdb_apo, pdb_holo
##                    print
##                    print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_holo:index2_seq_holo]
##                    print 'original'
##                    print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'][:10]
##                    print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][:10]
##                    print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'][-10:]
##                    print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][-10:]
##                    print
##                    print 'seq excl question marks'
##                    print index1_seq_apo,index2_seq_apo
##                    print index1_seq_holo,index2_seq_holo
##                    print
##                    print 'excluding question marks'
##                    print l_wt_seq_apo[:10]
##                    print l_wt_seq_holo[:10]
##                    print l_wt_seq_apo[-10:]
##                    print l_wt_seq_holo[-10:]
##                    print
##                    print 'coord overlap'
##                    print index1_coord_apo, index1_coord_holo
##                    print index2_coord_apo, index2_coord_holo
##                    print
##                    print 'aligned'
##                    print l_wt_seq_apo[index1_coord_apo:index1_coord_apo+10]
##                    print l_wt_seq_apo[index1_coord_apo-10:index2_coord_apo]
##                    print l_wt_seq_holo[index1_coord_holo:index1_coord_holo+10]
##                    print l_wt_seq_holo[index2_coord_holo-10:index2_coord_holo]
##                    print
##                    print l_coords_alpha_apo[0]
##                    print l_coords_alpha_holo[0]
##                    print l_coords_alpha_apo[-1]
##                    print l_coords_alpha_holo[-1]
##                    print
##                    print rmsd
##                    print pdb_apo, pdb_holo
##                    stop

                ## structural alignment
                for i_coord in range(len(l_coords_alpha_apo)):
                    l_coords_alpha_apo[i_coord] = numpy.dot(l_coords_alpha_apo[i_coord]-tv1,rm)+tv2

                ##
                ## vector from apo to holo
                ##
                vector_apo2holo = []
                for i in range(len(l_coords_alpha_holo)):
                    vector_apo2holo += [
                        l_coords_alpha_holo[i][0]-l_coords_alpha_apo[i][0],
                        l_coords_alpha_holo[i][1]-l_coords_alpha_apo[i][1],
                        l_coords_alpha_holo[i][2]-l_coords_alpha_apo[i][2],
                        ]
                vector_apo2holo = numpy.array(vector_apo2holo)

                ##
                ## calculate normal modes of apo structure
                ##
                cutoff = 10
                try:
                    matrix_hessian = NMA.hessian_calculation(l_coords_alpha_apo, cutoff, verbose = False)
                    eigenvectors, eigenvalues = NMA.diagonalize_hessian(matrix_hessian, verbose = False)
                except:
                    print pdb_apo, pdb_holo, len(l_coords_alpha_apo)
                    matrix_hessian = NMA.hessian_calculation(l_coords_alpha_apo, cutoff, verbose = False)
                    eigenvectors, eigenvalues = NMA.diagonalize_hessian(matrix_hessian, verbose = False)

                ##
                ## calculate overlap between normal modes and difference vector
                ##
##                for mode in range(6,12+1):
                for mode in range(6,6+1):
                    try:
                        eigenvector = eigenvectors[mode]
                    except:
                        print mode, len(eigenvectors)
                        print pdb_apo, pdb_holo
                        print len(l_coords_alpha_apo), len(l_coords_alpha_holo)
                        d_coords, l_coords_alpha_apo = mmCIF2coords.main(pdb_apo, d_apo)
                        d_coords, l_coords_alpha_holo = mmCIF2coords.main(pdb_holo, d_holo)
                        print len(l_coords_alpha_apo), len(l_coords_alpha_holo)

                        print pdb_apo, pdb_holo
                        print
                        print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_holo:index2_seq_holo]
                        print 'original'
                        print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'][:10]
                        print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][:10]
                        print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'][-10:]
                        print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][-10:]
                        print
                        print 'seq excl question marks'
                        print index1_seq_apo,index2_seq_apo
                        print index1_seq_holo,index2_seq_holo
                        print
                        print 'excluding question marks'
                        print l_wt_seq_apo[:10]
                        print l_wt_seq_holo[:10]
                        print l_wt_seq_apo[-10:]
                        print l_wt_seq_holo[-10:]
                        print
                        print 'coord overlap'
                        print index1_coord_apo, index1_coord_holo
                        print index2_coord_apo, index2_coord_holo
                        print
                        print 'aligned'
                        print l_wt_seq_apo[index1_coord_apo:index1_coord_apo+10]
                        print l_wt_seq_apo[index1_coord_apo-10:index2_coord_apo]
                        print l_wt_seq_holo[index1_coord_holo:index1_coord_holo+10]
                        print l_wt_seq_holo[index2_coord_holo-10:index2_coord_holo]
                        print
                        print l_coords_alpha_apo[0]
                        print l_coords_alpha_holo[0]
                        print l_coords_alpha_apo[-1]
                        print l_coords_alpha_holo[-1]
                        print
                        print pdb_apo, pdb_holo
                        stop
                        stop

                    overlap = abs(
                        numpy.dot(eigenvector,vector_apo2holo)
                        /
                        math.sqrt(
                            numpy.dot(eigenvector,eigenvector)
                            *
                            numpy.dot(vector_apo2holo,vector_apo2holo)
                            )
                        )

                    if mode == 6:
                        overlap_mode7 = overlap

                    if overlap > max_overlap[0]:
                        max_overlap = [
                            overlap,overlap_mode7,pdb_holo,pdb_apo,
                            mode,rmsd,d_pdbs_ligands_large[pdb_holo],
                            ]
                        print overlap, overlap_mode7, pdb_apo, pdb_holo

                if len(l_pdbs_holo) > 5 and len(l_pdbs_apo) > 25:
                    print 'cluster', cluster,
                    print 'holo', pdb_holo, l_pdbs_holo.index(pdb_holo), len(l_pdbs_holo),
                    print 'apo', pdb_apo,  l_pdbs_apo.index(pdb_apo), len(l_pdbs_apo),
                    print 'rmsd', round(rmsd,2), 'overlap', round(overlap_mode7,2)
                fd = open('rmsd_v_overlap/cluster%i.txt' %(cluster),'a')
                fd.write('%s %s\n' %(rmsd,overlap_mode7,))
                fd.close()
                fd = open('rmsd_v_overlap/all.txt','a')
                fd.write('%s %s\n' %(rmsd,overlap_mode7,))
                fd.close()

        if max_overlap[0] == 0:
            continue
            
        print cluster, max_overlap
        fd = open('max_cutoff%s.txt' %(cutoff),'a')
        fd.write('%s %s %s %s\n' %(cluster, str(max_overlap), pdb_apo, pdb_holo,) )
        fd.close()
        pdb_apo = max_overlap[3]
        pdb_holo = max_overlap[2]
##        d_clusters[cluster] = pdb_apo
        l_pdbs_apo_best += [[pdb_apo,pdb_holo,]]

    print l_pdbs_apo_best
    stop_end
    return l_pdbs_apo_best


def sequential_alignment_of_coordinates(
    l_coords_alpha_apo,l_coords_alpha_holo,
    d_apo,d_holo,
    pdb_apo, pdb_holo,
    ):

    ## identical number of coordinates and identical sequences
    if (
        len(l_coords_alpha_holo) == len(l_coords_alpha_apo)
##                    ):
        and
        d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'] == d_holo['_pdbx_poly_seq_scheme.pdb_mon_id']
        ):

        l_coords_alpha_apo = l_coords_alpha_apo
        l_coords_alpha_holo = l_coords_alpha_holo

    ## different number of coordinates or different sequences
    else:

        ##
        ## to be able to compare sequences replace with db_mon_id (assume protein from same species)
        ## it's ok to make permanent change to mmCIF
        ##
        for d in [d_apo, d_holo]:
            if '_struct_ref_seq_dif.seq_num' in d.keys():
                for i_struct_ref_seq_dif in range(len(d['_struct_ref_seq_dif.seq_num'])):

                    ## continue if expression tag / cloning artifact...
                    if d['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif] == '?':
                        continue

                    try:
                        seq_num = int(d['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif])
                    except:
                        print pdb_apo, pdb_holo
                        seq_num = int(d['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif])
                    print '_struct_ref_seq_dif.seq_num', seq_num
                    print '_struct_ref_seq_dif.mon_id', d['_struct_ref_seq_dif.mon_id'][i_struct_ref_seq_dif]
                    print '_struct_ref_seq_dif.db_mon_id', d['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif]

                    print '_pdbx_poly_seq_scheme.pdb_mon_id', d['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1]
                    print '_pdbx_poly_seq_scheme.seq_id', d['_pdbx_poly_seq_scheme.seq_id'][seq_num-1]

                    ## permanent change to mmCIF accepted
                    if d['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1] != '?':
                        d['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1] = d['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif]

        ##
        ## make list *after* replacement
        ##
        l_pdbx_poly_seq_scheme__pdb_mon_id__apo = list(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])
        l_pdbx_poly_seq_scheme__pdb_mon_id__holo = list(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])

        ##
        ## replace N-terminal cloning artifacts with question marks
        ## question marks are removed from sequence in next step...
        ##
        d = {
            'apo':{
                'd_coords':d_apo,'l_coords':l_coords_alpha_apo,
                '_pdbx_poly_seq_scheme.pdb_mon_id':l_pdbx_poly_seq_scheme__pdb_mon_id__apo,
                },
            'holo':{
                'd_coords':d_holo,'l_coords':l_coords_alpha_holo,
                '_pdbx_poly_seq_scheme.pdb_mon_id':l_pdbx_poly_seq_scheme__pdb_mon_id__holo,
                },
            }
        for k in [
            'apo','holo',
            ]:

            d_mmCIF = d[k]['d_coords']

            if not '_struct_ref_seq_dif.seq_num' in d_mmCIF.keys():
                break
            
            for i_struct_ref_seq_dif in range(len(d_mmCIF['_struct_ref_seq_dif.seq_num'])):

                if (
                    ## cloning artifact / expression tag
                    d_mmCIF['_struct_ref_seq_dif.details'][i_struct_ref_seq_dif] in ['CLONING ARTIFACT','EXPRESSION TAG',]
                    and
                    ## not in uniprot
                    d_mmCIF['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif] not in ['?','.',]
                    ):
                    print d_mmCIF['_struct_ref_seq_dif.details']
                    print d_mmCIF['_struct_ref_seq_dif.db_mon_id']
                    print k, pdb_apo, pdb_holo
                    stop
                if not (
                    ## cloning artifact / expression tag
                    d_mmCIF['_struct_ref_seq_dif.details'][i_struct_ref_seq_dif] in ['CLONING ARTIFACT','EXPRESSION TAG',]
                    and
                    ## not in uniprot
                    d_mmCIF['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif] in ['?','.',]
                    ):
                    continue

                seq_num = int(d_mmCIF['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif])

                ## remove coordinate *if* present for given residue
                if d_mmCIF['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1] != '?':
                    ## remove N-terminal coordinate
                    if d_mmCIF['_struct_ref_seq_dif.seq_num'][:seq_num] == [str(i) for i in range(1,seq_num+1)]:
                        d[k]['l_coords'] = d[k]['l_coords'][1:]
                    ## remove C-terminal coordinate
                    elif d_mmCIF['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif:] == [str(i) for i in range(seq_num,int(d_mmCIF['_struct_ref_seq_dif.seq_num'][-1])+1)]:
                        d[k]['l_coords'] = d[k]['l_coords'][:-1]
                    else:
                        print d_mmCIF['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif:]
                        print [str(i) for i in range(seq_num,int(d_mmCIF['_struct_ref_seq_dif.seq_num'][-1])+1)]
                        print
                        print seq_num, d_mmCIF['_struct_ref_seq_dif.seq_num']
                        print d_mmCIF['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif:]
                        print
                        print d_mmCIF['_struct_ref_seq_dif.seq_num'][:seq_num]
                        print '!='
                        print [str(i) for i in range(1,seq_num+1)]
                        print pdb_apo, pdb_holo, k
                        stop_not_Nterminal
                ## remove from sequence in next step...
                d[k]['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1] = '?'

        l_coords_alpha_apo = d['apo']['l_coords']
        l_coords_alpha_holo = d['holo']['l_coords']

        ## solution that works in all cases (of missing terminal residues, not terminal sequence differences)
        ## also for 2d59 and 2d5a, which have residues missing at the Nterm and Cterm, respectively
        ## and 3koi, 1eou
##        index1_seq_apo = next((i for i,v in enumerate(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
##        index2_seq_apo = len(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))
##        index1_seq_holo = next((i for i,v in enumerate(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
##        index2_seq_holo = len(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))
        index1_seq_apo = next((i for i,v in enumerate(d['apo']['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
        index2_seq_apo = len(d['apo']['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d['apo']['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))
        index1_seq_holo = next((i for i,v in enumerate(d['holo']['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
        index2_seq_holo = len(d['holo']['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d['holo']['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))

        ## cap apo and holo sequence to achieve same length and same sequence
##        l_wt_seq_apo = list(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_apo:index2_seq_apo])
##        l_wt_seq_holo = list(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_holo:index2_seq_holo])
        l_wt_seq_apo = list(d['apo']['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_apo:index2_seq_apo])
        l_wt_seq_holo = list(d['holo']['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_holo:index2_seq_holo])
        if len(l_wt_seq_apo) != len(l_coords_alpha_apo) or len(l_wt_seq_holo) != len(l_coords_alpha_holo):
            print len(l_wt_seq_apo)
            print len(l_coords_alpha_apo)
            print len(l_wt_seq_holo)
            print len(l_coords_alpha_holo)
            print pdb_apo, pdb_holo

            d_coords, l_coords_alpha_holo = mmCIF2coords.main(pdb_holo, d_holo)
            print len(l_coords_alpha_holo)

            stop

        if '-'.join(l_wt_seq_apo) in '-'.join(l_wt_seq_holo):
            index = '-'.join(l_wt_seq_holo).index('-'.join(l_wt_seq_apo))
            index /= 4
            index1_coord_apo = 0
            index2_coord_apo = len(l_wt_seq_apo)
            index1_coord_holo = index
            index2_coord_holo = index+len(l_wt_seq_apo)
        elif '-'.join(l_wt_seq_holo) in '-'.join(l_wt_seq_apo):
            index = '-'.join(l_wt_seq_apo).index('-'.join(l_wt_seq_holo))
            index /= 4
            index1_coord_apo = index
            index2_coord_apo = index+len(l_wt_seq_holo)
            index1_coord_holo = 0
            index2_coord_holo = len(l_wt_seq_holo)
        else:
            x1 = commonOverlapNaive('-'.join(l_wt_seq_apo),'-'.join(l_wt_seq_holo),)
            x2 = commonOverlapNaive('-'.join(l_wt_seq_holo),'-'.join(l_wt_seq_apo),)
            if x1 > x2:
                index1_coord_apo = -((x1-3)/4+1)
                index2_coord_apo = len(l_wt_seq_apo)
                index1_coord_holo = 0
                index2_coord_holo = ((x1-3)/4+1)
            elif x2 > x1:
                index1_coord_holo = -((x2-3)/4+1) ## correct!!!
                index2_coord_holo = len(l_wt_seq_holo) ## correct!!!
                index1_coord_apo = 0 ## correct!!!
                index2_coord_apo = ((x2-3)/4+1) ## correct!!!
            else:
                print
                print 'x1', x1, 'x2', x2
                print
                print l_wt_seq_apo
                print l_wt_seq_holo
                print pdb_apo
                print pdb_holo
                for y in range(len(l_wt_seq_apo)):
                    if l_wt_seq_apo[y] != l_wt_seq_holo[y]:
                        print y
                        print l_wt_seq_apo[:y+1]
                        print l_wt_seq_holo[:y+1]
                        break
##                print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id']

                print '>%s _struct_ref.pdbx_seq_one_letter_code' %(pdb_apo)
                db_seq_apo = ''.join(d_apo['_struct_ref.pdbx_seq_one_letter_code']).upper().replace(' ','')
                print db_seq_apo
                print '>%s _struct_ref.pdbx_seq_one_letter_code' %(pdb_holo)
                db_seq_holo = ''.join(d_holo['_struct_ref.pdbx_seq_one_letter_code']).upper().replace(' ','')
                print db_seq_holo
                
                print '>%s _entity_poly.pdbx_seq_one_letter_code_can' %(pdb_apo)
                print ''.join(d_apo['_entity_poly.pdbx_seq_one_letter_code_can'])
                print '>%s _entity_poly.pdbx_seq_one_letter_code_can' %(pdb_holo)
                print ''.join(d_holo['_entity_poly.pdbx_seq_one_letter_code_can'])

                print
                print d_apo['_struct_ref_seq.pdbx_db_accession']
                print d_holo['_struct_ref_seq.pdbx_db_accession']

                import urllib2
                UNP = ''.join(d_apo['_struct_ref_seq.pdbx_db_accession'])
                if UNP[0] != 'P':
                    UNP = ''.join(d_holo['_struct_ref_seq.pdbx_db_accession'])
                url = urllib2.urlopen('http://www.uniprot.org/uniprot/%s.fasta' %(UNP))
                lines = url.readlines()
                
                seq = ''.join(lines[1:]).replace('\n','').strip()
                print
                print seq
                print seq == db_seq_apo
                print seq == db_seq_holo
                print len(seq), len(db_seq_apo), len(db_seq_holo)
                for y in range(len(seq)):
                    if seq[y] != db_seq_apo[y]:
                        print 'apo', pdb_apo
                        print seq[:y+1]
                        print db_seq_apo[:y+1]
                        break
                for y in range(len(seq)):
                    if seq[y] != db_seq_holo[y]:
                        print 'holo', pdb_holo
                        print seq[:y+1]
                        print db_seq_holo[:y+1]
                        break

                print 'x1', x1, 'x2', x2
                print pdb_apo, pdb_holo
                print db_seq_apo == db_seq_holo

                stop_no_overlap

            if l_wt_seq_apo[index1_coord_apo:index2_coord_apo] != l_wt_seq_holo[index1_coord_holo:index2_coord_holo]:
                print l_wt_seq_apo
                print l_wt_seq_holo
                print l_wt_seq_apo[index1_coord_apo:index2_coord_apo]
                print l_wt_seq_holo[index1_coord_holo:index2_coord_holo]
                print index1_coord_apo, index2_coord_apo
                print index1_coord_holo, index2_coord_holo
                print pdb_apo, pdb_holo

                print 'x', x, (x-3)/4
                print 'index_apo', index1_coord_apo, index2_coord_apo
                print 'index_holo', index1_coord_holo, index2_coord_holo
                print 'raw', l_wt_seq_apo
                print 'raw', l_wt_seq_holo
##                        print '***', l_wt_seq_apo[index1_coord_apo:index2_coord_apo], index1_coord_apo, index2_coord_apo
##                        print '***', l_wt_seq_holo[index1_coord_holo:index2_coord_holo], index1_coord_holo, index2_coord_holo
                print pdb_apo, pdb_holo
                print 'aln', l_wt_seq_holo[-(x-3)/4:len(l_wt_seq_holo)]
                print 'aln', l_wt_seq_apo[0:(x-3)/4]
                print 'text1[-x:] == text2[:x]'
                print 'x', x, (x-3)/4, len(l_wt_seq_apo), len(l_wt_seq_holo)

                for y in range(len(l_wt_seq_apo)):
                    if l_wt_seq_apo[y] != l_wt_seq_holo[y]:
                        print l_wt_seq_apo[:y+1]
                        print l_wt_seq_holo[:y+1]
                        break

                stop1

        l_coords_alpha_apo = l_coords_alpha_apo[index1_coord_apo:index2_coord_apo]
        l_coords_alpha_holo = l_coords_alpha_holo[index1_coord_holo:index2_coord_holo]
        if len(l_coords_alpha_apo) != len(l_coords_alpha_holo) or len(l_coords_alpha_apo) < 10:
            print
            print 'poly_seq'
            print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id']
            print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id']
            print
            print index1_coord_apo, index2_coord_apo
            print index1_coord_holo, index2_coord_holo

            for y in range(min(len(l_wt_seq_apo[index1_coord_apo:index2_coord_apo]),len(l_wt_seq_holo))):
                if l_wt_seq_apo[y] != l_wt_seq_holo[y]:
                    print l_wt_seq_apo[:y+1]
                    print l_wt_seq_holo[:y+1]
                    break
            print 'coords_apo', len(l_coords_alpha_apo)
            print 'coords_holo', len(l_coords_alpha_holo)

            d_coords, l_coords_alpha_apo = mmCIF2coords.main(pdb_apo, d_apo)
            d_coords, l_coords_alpha_holo = mmCIF2coords.main(pdb_holo, d_holo)
            print 'coords_apo', len(l_coords_alpha_apo)
            print 'coords_holo', len(l_coords_alpha_holo)
            l_coords_alpha_apo = l_coords_alpha_apo[index1_coord_apo:index2_coord_apo]
            l_coords_alpha_holo = l_coords_alpha_holo[index1_coord_holo:index2_coord_holo]
            print 'coords_apo', len(l_coords_alpha_apo)
            print 'coords_holo', len(l_coords_alpha_holo)

            print 'seq'
            print len(l_wt_seq_apo)
            print len(l_wt_seq_holo)
            print len(l_wt_seq_apo[index1_coord_apo:index2_coord_apo])
            print len(l_wt_seq_holo[index1_coord_holo:index2_coord_holo])

            print 'poly_seq_scheme', len(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])
            print 'poly_seq_scheme', len(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])
            print 'wt', l_wt_seq_apo[index1_coord_apo:index2_coord_apo]
            print 'wt', l_wt_seq_holo[index1_coord_holo:index2_coord_holo]
            print 'x1', x1, 'x2', x2
            print pdb_apo, pdb_holo
            stop_diff_len

    return l_coords_alpha_apo, l_coords_alpha_holo


def commonOverlapNaive(text1,text2):

##    x = commonOverlapNaive('Fire at Will','William Riker is number one',)
##    print x
    ## problem if seq1 = xxxACDExxx and seq2 = xxxACDExxx

    x = min(len(text1), len(text2))
    while x > 0:
        if text1[-x:] == text2[:x]:
            break
        x -= 1

    return x


def identify_sequence_similar_apo_pdbs(self,l_pdbs,):

    stop_already_used_clusters

    l_pdbs_sequence_unique

    if not os.path.isfile('bc-95.out'):
        urllines = urllib2.urlopen('ftp://resources.rcsb.org/sequence/clusters/bc-95.out')
        lines = urllines.readlines()
        fd = open('bc-95.out','w')
        fd.writelines(lines)
        fd.close()
    else:
        fd = open('bc-95.out','r')
        lines = fd.readlines()
        fd.close()

    for line in lines:
        l_pdbs_line = line
        print l_pdbs
        stop
        if len( set(l_pdbs) & set(l_pdbs_line) ) > 1:
            print l_pdbs_line
            print set(l_pdbs) & set(l_pdbs_line)
            stop

    return l_pdbs_sequence_unique


if __name__ == '__main__':
    main()

##                        if l_wt_seq_apo[:len(l_wt_seq_holo)] == l_wt_seq_holo:
##                            print 'a'
##                            index1_coord_apo = 0
##                            index2_coord_apo = len(l_coords_alpha_apo)-(len(l_wt_seq_apo)-len(l_wt_seq_holo))
##                            index1_coord_holo = 0
##                            index2_coord_holo = len(l_coords_alpha_holo)
##                            print index1_coord_apo, index2_coord_apo
##                            print index1_coord_holo, index2_coord_holo
##                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][0]
##                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][0]
##                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][-1]
##                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][-1]
##                            if pdb_apo == '3koi':
##                                stopa
####                            stopa
##                        elif l_wt_seq_apo[-len(l_wt_seq_holo):] == l_wt_seq_holo:
##                            print 'b'
##                            index1_coord_apo = (len(l_wt_seq_apo)-len(l_wt_seq_holo))
##                            index2_coord_apo = len(l_coords_alpha_apo)
##                            index1_coord_holo = 0
##                            index2_coord_holo = len(l_coords_alpha_holo)
####                            print index1_coord_apo, index2_coord_apo
####                            print index1_coord_holo, index2_coord_holo
####                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][0]
####                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][0]
####                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][-1]
####                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][-1]
####                            stopb
##                        elif l_wt_seq_holo[:len(l_wt_seq_apo)] == l_wt_seq_apo:
##                            print 'c'
##                            index1_coord_apo = 0
##                            index2_coord_apo = len(l_coords_alpha_apo)
##                            index1_coord_holo = 0
##                            index2_coord_holo = len(l_coords_alpha_holo)-(len(l_wt_seq_holo)-len(l_wt_seq_apo))
####                            print index1_coord_apo, index2_coord_apo
####                            print index1_coord_holo, index2_coord_holo
####                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][0]
####                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][0]
####                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][-1]
####                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][-1]
####                            stopc
##                        elif l_wt_seq_holo[-len(l_wt_seq_apo):] == l_wt_seq_apo:
##                            print 'd'
##                            index1_coord_apo = 0
##                            index2_coord_apo = len(l_coords_alpha_apo)
##                            index1_coord_holo = len(l_wt_seq_holo)-len(l_wt_seq_apo)
##                            index2_coord_holo = len(l_coords_alpha_holo)
##                            print index1_coord_apo, index2_coord_apo
##                            print index1_coord_holo, index2_coord_holo
##                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][0]
##                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][0]
##                            print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][-1]
##                            print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][-1]
####                            stopd
##                        else:
##                            
##                            print
##                            print d_apo['_pdbx_poly_seq_scheme.pdb_mon_id']
##                            print
##                            print d_holo['_pdbx_poly_seq_scheme.pdb_mon_id']
##                            print
##                            print l_wt_seq_apo
##                            print
##                            print l_wt_seq_holo
##                            if ''.join(l_wt_seq_apo) in ''.join(l_wt_seq_holo):
##                                index = ''.join(l_wt_seq_holo).index(''.join(l_wt_seq_apo))
##                                index1_coord_apo = 0
##                                index2_coord_apo = len(l_coords_alpha_apo)
##                                index1_coord_holo = index/3
##                                index2_coord_holo = index/3+len(l_coords_alpha_apo)
####                                print index1_coord_apo, index2_coord_apo
####                                print index1_coord_holo, index2_coord_holo
####                                print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][0]
####                                print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][0]
####                                print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][-1]
####                                print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][-1]
####                                print pdb_apo, pdb_holo
####                                stope
##                            elif ''.join(l_wt_seq_holo) in ''.join(l_wt_seq_apo):
##                                index = ''.join(l_wt_seq_apo).index(''.join(l_wt_seq_holo))
##                                index1_coord_apo = index/3
##                                index2_coord_apo = index/3+len(l_coords_alpha_holo)
##                                index1_coord_holo = 0
##                                index2_coord_holo = len(l_coords_alpha_holo)
####                                print index1_coord_apo, index2_coord_apo
####                                print index1_coord_holo, index2_coord_holo
####                                print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][0]
####                                print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][0]
####                                print l_coords_alpha_apo[index1_coord_apo:index2_coord_apo][-1]
####                                print l_coords_alpha_holo[index1_coord_holo:index2_coord_holo][-1]
####                                print pdb_apo, pdb_holo
####                                stopf
##                            else:
####                            print d_apo['_struct_ref_seq_dif.seq_num']
####                            print d_holo['_struct_ref_seq_dif.seq_num']
##                                print l_wt_seq_apo
##                                print l_wt_seq_holo
##                                stop
##                    else:
##                        index1_coord_apo = max(0,index1_seq_holo-index1_seq_apo)
##                        index2_coord_apo = len(l_coords_alpha_apo)+min(0,index2_seq_holo-index2_seq_apo)
##                        index1_coord_holo = max(0,index1_seq_apo-index1_seq_holo)
##                        index2_coord_holo = len(l_coords_alpha_holo)+min(0,index2_seq_apo-index2_seq_holo)


##                    s_seqres_holo = ''.join(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id']).replace('?','')
##                    s_seqres_apo = ''.join(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id']).replace('?','')
##                    ## check
##                    if len(s_seqres_holo) % 3 != 0 or len(s_seqres_apo) % 3 != 0:
##                        stopstop
##                    if s_seqres_holo in s_seqres_apo:
##                        index1 = s_seqres_apo.index(s_seqres_holo)
##                        index2 = index1+len(s_seqres_holo)
##                        index1 /= 3
##                        index2 /= 3
##                        l_coords_alpha_apo = l_coords_alpha_apo[index1:index2]
##                    elif s_seqres_apo in s_seqres_holo:
##                        index1 = s_seqres_holo.index(s_seqres_apo)
##                        index2 = index1+len(s_seqres_apo)
##                        index1 /= 3
##                        index2 /= 3
##                        l_coords_alpha_holo = l_coords_alpha_holo[index1:index2]
##                    else:
##                        stop

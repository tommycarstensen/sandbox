## built ins
import math
## add ons
import numpy
## local
import create_apo_holo_dataset
import sys
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('/home/tc/svn/Protool/')
import geometry
sys.path.append('/home/tc/svn/GoodVibes')
import NMA

path_mmCIF = '/media/WDMyBook1TB/2TB/mmCIF'

def main():

    set_pdbs = exclude_include()
    l_pdbs_remove = [
        '4a3h','2wf5','1arl','1ee3', ## incorrect _struct_ref_seq.pdbx_db_accession
        '1uyd','1uye','1uyf','2byh','2byi', ## remediation _struct_ref_seq_dif
        '2xdu','3dn8','3dna','1ps3','1ouf','1l35','2eun','1rtc','1zon', ## _struct_ref_seq_dif missing
        '1pwl','1pwm','2fz8','2fz9', ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
        ]
    set_pdbs.remove('1f92') ## remediation _struct_ref_seq_dif incorrect residue number
    set_pdbs.remove('2f6f') ## remediation _pdbx_poly_seq_scheme.auth_mon_id wrong
    set_pdbs.remove('3a5j') ## remediation _struct_ref_seq_dif.db_mon_id is ? but should be MET
    set_pdbs.remove('2rhx') ## remediation _struct_ref_seq_dif.db_mon_id is ? but should be SER
    set_pdbs.remove('2fzb') ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
    set_pdbs.remove('2fzd') ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
    set_pdbs.remove('3dn5') ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
    set_pdbs.remove('1x96') ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
    set_pdbs.remove('1x97') ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
    set_pdbs.remove('1x98') ## remediation incorrect _struct_ref.pdbx_seq_one_letter_code
    set_pdbs.remove('1z3n') ## GenBank DBref - not an error...
    set_pdbs.remove('1z8a') ## GenBank DBref - not an error...
    set_pdbs.remove('1z89') ## GenBank DBref - not an error...
    set_pdbs.remove('2pf8') ## stupid use of alt_ids (C for highest occupancy and only altloc)
    set_pdbs.remove('2pyr') ## stupid use of alt_ids (G and R)
    set_pdbs.remove('3pdn') ## stupid use of alt_ids (B and C)
    set_pdbs.remove('2v4c') ## alt_id B used for 100% occupancy atoms
    set_pdbs.remove('1jxt') ## weird alt_id microheterogeneity...
    set_pdbs.remove('1jxu') ## weird alt_id microheterogeneity...
    set_pdbs.remove('1jxw') ## weird alt_id microheterogeneity...
    set_pdbs.remove('1jxx') ## weird alt_id microheterogeneity...
    set_pdbs.remove('1jxy') ## weird alt_id microheterogeneity...
##    set_pdbs.remove('1ac4') ## multiple strains and taxonomy ids but all same organism (S. cerevisiae)...
##    set_pdbs.remove('1ac8') ## multiple strains and taxonomy ids but all same organism (S. cerevisiae)...
##    set_pdbs.remove('1aeb') ## multiple strains and taxonomy ids but all same organism (S. cerevisiae)...
##    set_pdbs.remove('2rbt') ## multiple strains and taxonomy ids but all same organism (S. cerevisiae)... UNP A7A026, TAX 307796, STRAIN YJM789
##    set_pdbs.remove('2rbu') ## multiple strains and taxonomy ids but all same organism (S. cerevisiae)... UNP A7A026, TAX 307796, STRAIN YJM789
##    set_pdbs.remove('2rbv') ## multiple strains and taxonomy ids but all same organism (S. cerevisiae)... UNP A7A026, TAX 307796, STRAIN YJM789
    for pdb in l_pdbs_remove:
        set_pdbs.remove(pdb)

    fd = open('%s/bc-100.out' %(path_mmCIF),'r')
    lines = fd.readlines()
    fd.close()

    for i_line in range(len(lines)):
        cluster = i_line
        if cluster < 4816:
            continue
##        if cluster not in [5,]:
##            continue
        line = lines[i_line]
        l_pdbs = line.lower().split()
        l_pdbs.sort()
        for i_pdb in range(len(l_pdbs)):
            l_pdbs[i_pdb] = l_pdbs[i_pdb][:4]

        for i_pdb1 in range(0,len(l_pdbs)-1):

            pdb1 = l_pdbs[i_pdb1]

##            if pdb1 != '1t49': ## tmp!!!
##                continue

            if not pdb1 in set_pdbs:
                continue

            print pdb1
            stop

            d_mmCIF1 = parse_mmCIF.main(pdb1,)

            bool_monomeric = check_monomeric(d_mmCIF1)
            if bool_monomeric == False:
                if i_pdb1 == 0:
                    break
                else:
                    continue

            bool_remediation_modres = check_modres(d_mmCIF1,pdb1,)
            if bool_remediation_modres == True:
                continue

            if '_struct_ref_seq_dif.details' in d_mmCIF1.keys():
                if 'DELETION' in d_mmCIF1['_struct_ref_seq_dif.details']:
                    continue

            for i_entity in range(len(d_mmCIF1['_entity.id'])):
                if d_mmCIF1['_entity.type'][i_entity] == 'polymer':
                    if int(d_mmCIF1['_entity.pdbx_number_of_molecules'][i_entity]) != 1:
                        print d_mmCIF1['_entity.pdbx_number_of_molecules']
                        print pdb1, cluster
                        stop

            SG1 = d_mmCIF1['_symmetry.space_group_name_H-M']

            for i_pdb2 in range(i_pdb1+1,len(l_pdbs)):

                pdb2 = l_pdbs[i_pdb2]

##                if pdb2 != '2pf8': ## tmp!!!
##                    continue

##                if pdb1 != '3fui' or pdb2 != '3fuj':
##                    continue

                if not pdb2 in set_pdbs:
                    continue

                d_mmCIF2 = parse_mmCIF.main(pdb2,)

                bool_monomeric = check_monomeric(d_mmCIF2)
                if bool_monomeric == False:
                    continue

                bool_remediation_modres = check_modres(d_mmCIF2,pdb2,)
                if bool_remediation_modres == True:
                    continue

                if '_struct_ref_seq_dif.seq_num' in d_mmCIF2.keys():
                    if 'DELETION' in d_mmCIF2['_struct_ref_seq_dif.details']:
                        continue

                ## biounit monomeric?
                for i_entity in range(len(d_mmCIF2['_entity.id'])):
                    if d_mmCIF2['_entity.type'][i_entity] == 'polymer':
                        if int(d_mmCIF2['_entity.pdbx_number_of_molecules'][i_entity]) != 1:
                            continue

                SG2 = d_mmCIF2['_symmetry.space_group_name_H-M']

                if SG1 != SG2:
                    continue

                ## parse coordinates again after being shortened in previous loop
                try:
                    d_coords1, l_coords_alpha1 = mmCIF2coords.main(pdb1, d_mmCIF1)
                except:
                    fd = open('remediation_atom_site.label_alt_id.txt','a')
                    fd.write('%s\n' %(pdb1,))
                    fd.close()
                try:
                    d_coords2, l_coords_alpha2 = mmCIF2coords.main(pdb2, d_mmCIF2)
                except:
                    fd = open('remediation_atom_site.label_alt_id.txt','a')
                    fd.write('%s\n' %(pdb2,))
                    fd.close()

                ## align sequences/coordinates
                try:
                    l_coords_alpha1, l_coords_alpha2 = create_apo_holo_dataset.sequential_alignment_of_coordinates(
                        l_coords_alpha1, l_coords_alpha2,
                        d_mmCIF1, d_mmCIF2,
                        pdb1, pdb2,
                        )
                except:
                    fd = open('remediation_struct_ref_seq_dif.txt','a')
                    fd.write(
                        '%s %s %s %s\n' %(
                            pdb1,pdb2,
                            d_mmCIF1['_struct_ref_seq.pdbx_db_accession'],
                            d_mmCIF2['_struct_ref_seq.pdbx_db_accession'],
                            )
                        )
                    fd.close()
                    continue
                if len(l_coords_alpha1) != len(l_coords_alpha2):
                    print d_mmCIF1['_pdbx_poly_seq_scheme.pdb_mon_id']
                    print d_mmCIF2['_pdbx_poly_seq_scheme.pdb_mon_id']
                    print 'coords', len(l_coords_alpha1), len(l_coords_alpha2)
                    print 'seq', len(d_mmCIF1['_pdbx_poly_seq_scheme.pdb_mon_id'])
                    print 'seq', len(d_mmCIF2['_pdbx_poly_seq_scheme.pdb_mon_id'])
                    print pdb1, pdb2
                    d_coords1, l_coords_alpha1 = mmCIF2coords.main(pdb1, d_mmCIF1)
                    d_coords1, l_coords_alpha2 = mmCIF2coords.main(pdb1, d_mmCIF2)
                    print len(l_coords_alpha1), len(l_coords_alpha2)
                    stop
                    continue

                ##
                ## align structure 1 and 2
                ##
                instance_geometry = geometry.geometry()
                rmsd = instance_geometry.superpose(l_coords_alpha1,l_coords_alpha2)
                tv1 = instance_geometry.fitcenter
                rm = instance_geometry.rotation
                tv2 = instance_geometry.refcenter

                ## structural alignment
                for i_coord in range(len(l_coords_alpha2)):
                    l_coords_alpha2[i_coord] = numpy.dot(l_coords_alpha2[i_coord]-tv1,rm)+tv2

                ##
                ## vector from structure 1 to 2
                ##
                vector = []
                for i in range(len(l_coords_alpha1)):
                    vector += [
                        l_coords_alpha1[i][0]-l_coords_alpha2[i][0],
                        l_coords_alpha1[i][1]-l_coords_alpha2[i][1],
                        l_coords_alpha1[i][2]-l_coords_alpha2[i][2],
                        ]
                vector = numpy.array(vector)

                ##
                ## calculate normal modes of structure 1
                ##
                cutoff = 10
                try:
                    matrix_hessian1 = NMA.hessian_calculation(l_coords_alpha1, cutoff, verbose = False)
                    eigenvectors1, eigenvalues1 = NMA.diagonalize_hessian(matrix_hessian1, verbose = False)
                    matrix_hessian2 = NMA.hessian_calculation(l_coords_alpha2, cutoff, verbose = False)
                    eigenvectors2, eigenvalues2 = NMA.diagonalize_hessian(matrix_hessian2, verbose = False)
                except:
                    continue

                ##
                ## calculate overlap between normal modes and difference vector
                ##
                eigenvector1 = eigenvectors1[6]
                eigenvector2 = eigenvectors2[6]

                overlap1 = calc_overlap(eigenvector1,vector)
                overlap2 = calc_overlap(eigenvector2,vector)
                overlap3a = calc_overlap(eigenvector1,eigenvector2)
                overlap3b = calc_overlap(eigenvectors1[6],eigenvectors2[7])
                overlap3c = calc_overlap(eigenvectors1[7],eigenvectors2[6])
                overlap3 = max(overlap3a,overlap3b,overlap3c)

                fd = open('rmsd_v_overlap2/cluster%i.txt' %(i_line),'a')
                fd.write('%s %s\n' %(rmsd,overlap1))
                fd.close()
                fd = open('rmsd_v_overlap2/cluster%i.txt' %(i_line),'a')
                fd.write('%s %s\n' %(rmsd,overlap2))
                fd.close()
                fd = open('rmsd_v_overlap2/cluster%i_ev_v_ev.txt' %(i_line),'a')
                fd.write('%s %s\n' %(rmsd,overlap3a))
                fd.close()
                fd = open('rmsd_v_overlap2/cluster%i_ev_v_ev_max.txt' %(i_line),'a')
                fd.write('%s %s\n' %(rmsd,overlap3))
                fd.close()
                print pdb1, pdb2, 'cluster', i_line, 'size', len(l_pdbs),
                print 'overlap', '%4.2f' %(round(overlap1,2)), '%4.2f' %(round(overlap2,2)), '%4.2f' %(round(overlap3,2)), 'rmsd', '%4.2f' %(round(rmsd,2))

    return


def calc_overlap(eigenvector,vector):

    overlap = abs(
        numpy.dot(eigenvector,vector)
        /
        math.sqrt(
            numpy.dot(eigenvector,eigenvector)
            *
            numpy.dot(vector,vector)
            )
        )

    return overlap


def exclude_include():

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

    set_exclude_pdbs = set()
    for fn in l_exclude_fns:
        fd = open('%s/list%s.txt' %(path_mmCIF,fn,),'r')
        s = fd.read()
        fd.close()
        set_exclude_pdbs |= set(s.lower().split('\n'))

    set_pdbs = set_include_pdbs-set_exclude_pdbs

    return set_pdbs


def check_modres(d_mmCIF,pdb,):

    bool_remediation_modres = False

    if (
        '(' in ''.join(d_mmCIF['_entity_poly.pdbx_seq_one_letter_code'])
        and
        '_pdbx_struct_mod_residue.id' not in d_mmCIF.keys()
        ):
        seq = ''.join(d_mmCIF['_entity_poly.pdbx_seq_one_letter_code'])
        modres = seq[seq.index('(')+1:seq.index(')')]
        chem_comp_type = d_mmCIF['_chem_comp.type'][d_mmCIF['_chem_comp.id'].index(modres)]
        fd = open('remediation_pdbx_struct_mod_residue.txt','a')
        fd.write('%s %s %s\n' %(pdb,modres,chem_comp_type,))
        fd.close()
        bool_remediation_modres = True

    return bool_remediation_modres


def check_monomeric(d_mmCIF):

    bool_monomeric = True

    if '_pdbx_struct_assembly.oligomeric_details' in d_mmCIF.keys():
        if (
            d_mmCIF['_pdbx_struct_assembly.oligomeric_details']
            !=
            len(d_mmCIF['_pdbx_struct_assembly.oligomeric_details'])*['monomeric']
            ):
            bool_monomeric = False

    return bool_monomeric


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

        ## solution that works in all cases (of missing terminal residues, not terminal sequence differences)
        ## also for 2d59 and 2d5a, which have residues missing at the Nterm and Cterm, respectively
        ## and 3koi, 1eou
        index1_seq_apo = next((i for i,v in enumerate(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
        index2_seq_apo = len(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))
        index1_seq_holo = next((i for i,v in enumerate(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
        index2_seq_holo = len(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))

        ## to be able to compare sequences replace with db_mon_id (assume protein from same species)
        for d in [d_apo, d_holo]:
            if '_struct_ref_seq_dif.seq_num' in d.keys():
                for i_struct_ref_seq_dif in range(len(d['_struct_ref_seq_dif.seq_num'])):
                    ## continue if expression tag...
                    if d['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif] == '?':
                        continue

                    seq_num = int(d['_struct_ref_seq_dif.seq_num'][i_struct_ref_seq_dif])
                    print seq_num
                    print d['_struct_ref_seq_dif.mon_id'][i_struct_ref_seq_dif]
                    print d['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif]

                    print d['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1]
                    print d['_pdbx_poly_seq_scheme.seq_id'][seq_num-1]

                    if d['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1] != '?':
                        d['_pdbx_poly_seq_scheme.pdb_mon_id'][seq_num-1] = d['_struct_ref_seq_dif.db_mon_id'][i_struct_ref_seq_dif]

        ## "remove" missing terminal residues from sequence
        l_wt_seq_apo = list(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_apo:index2_seq_apo])
        l_wt_seq_holo = list(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'][index1_seq_holo:index2_seq_holo])

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
            x = commonOverlapNaive('-'.join(l_wt_seq_apo),'-'.join(l_wt_seq_holo))
            if x == 0:
                x = commonOverlapNaive('-'.join(l_wt_seq_holo),'-'.join(l_wt_seq_apo),)
                index1_coord_holo = len(l_wt_seq_holo)-(((x-3)/4)+1)
                index2_coord_holo = len(l_wt_seq_holo)
                index1_coord_apo = 0
                index2_coord_apo = (((x-3)/4)+1)
            else:
                index1_coord_holo = 0
                index2_coord_holo = (((x-3)/4)+1)
                index1_coord_apo = len(l_wt_seq_apo)-(((x-3)/4)+1)
                index2_coord_apo = len(l_wt_seq_apo)
            if x == 0:
                print l_wt_seq_apo
                print l_wt_seq_holo
                print pdb_apo
                print pdb_holo
                for y in range(len(l_wt_seq_apo)):
                    print y, l_wt_seq_apo[y] == l_wt_seq_holo[y], l_wt_seq_apo[y], l_wt_seq_holo[y]
                    if l_wt_seq_apo[y] != l_wt_seq_holo[y]:
                        print l_wt_seq_apo[:y+1]
                        print l_wt_seq_holo[:y+1]
                        break
                print pdb_holo
                stop

##                    if l_wt_seq_holo[index1_coord_holo:index2_coord_holo] != l_wt_seq_apo[index1_coord_apo:index2_coord_apo]:
##                        print '#####'
##                        print x, (x-3)/4
##                        print '-'.join(l_wt_seq_apo)[:x].split('-')
##                        print '-'.join(l_wt_seq_holo)[-x:].split('-')
##                        print index1_coord_holo
##                        print index2_coord_holo
##                        print index1_coord_apo
##                        print index2_coord_apo
##                        print l_wt_seq_holo[index1_coord_holo:index2_coord_holo]
##                        print l_wt_seq_apo[index1_coord_apo:index2_coord_apo]
##                        print l_wt_seq_holo[index1_coord_holo:index2_coord_holo] == l_wt_seq_apo[index1_coord_apo:index2_coord_apo]
##                        stop

        l_coords_alpha_apo = l_coords_alpha_apo[index1_coord_apo:index2_coord_apo]
        l_coords_alpha_holo = l_coords_alpha_holo[index1_coord_holo:index2_coord_holo]
        if len(l_coords_alpha_apo) != len(l_coords_alpha_holo):
            print len(l_coords_alpha_apo)
            print len(l_coords_alpha_holo)
            print len(d_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])
            print len(d_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])
            print pdb_apo, pdb_holo
            print index1_coord_apo, index2_coord_apo
            print index1_coord_holo, index2_coord_holo
            stop_diff_len

    return l_coords_alpha_apo, l_coords_alpha_holo


def commonOverlapNaive(text1,text2):

##    x = commonOverlapNaive('Fire at Will','William Riker is number one',)
##    print x

    x = min(len(text1), len(text2))
    while x > 0:
        if text1[-x:] == text2[:x]:
            break
        x -= 1

    return x


if __name__ == '__main__':
    main()

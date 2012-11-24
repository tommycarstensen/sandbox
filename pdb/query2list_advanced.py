import parse_mmCIF

path = '/media/Tommy/mmCIF'
path = '/media/WDMyBook1TB/2TB/mmCIF'

def main():

    loop('one_polypeptide','_one_polypeptide')
    stop_finished

    loop('one_polysaccharide','_one_polysaccharide')
    stop_finished

    cluster_by_sequence()

    loop('modres_not_MSE','_pdbx_struct_mod_residue_notMSE')

##    unobs_nonterminal_residues()
    
    ## this method is not entirely correct... e.g. 1kwr...
    unobs_nonterminal_atoms_alpha()

    return


def one_polysaccharide(pdb,):

    l_data_categories = [
        '_entity',
        '_chem_comp',
        '_entity_poly',
        ]
    d = parse_mmCIF.main(
        pdb,
        l_data_categories = l_data_categories,
        )

    bool_append = False

    bool_polysaccharide = False
    if '_chem_comp.type' in d.keys():
        for chem_comp_type in d['_chem_comp.type']:
            if chem_comp_type.lower() in [
                'd-saccharide 1,4 and 1,4 linking', # 3amm
                'l-saccharide','d-saccharide','saccharide'
                ]:
                bool_polysaccharide = True
                break
##            elif 'acchar' in chem_comp_type.lower():
##                print d
##                print chem_comp_type
##                print pdb
##                print set(['D-saccharide','saccharide'])&set(d['_chem_comp.type'])
##                stop
##    else:
##        print pdb
##        stop

    count_polymer_sugar = 0
    bool_monosaccharide = False ## included to exclude 1a14 which contains polymers and monomers
    for i in range(len(d['_entity.type'])):
        entity_type = d['_entity.type'][i]
        if entity_type in [
            'polymer',
            ]:
            if d['_entity.pdbx_description'][i][:7] == 'SUGAR (':
                count_polymer_sugar += int(d['_entity.pdbx_number_of_molecules'][i])
                continue
##            ## polypeptide or polynucleotide (just a check)
##            elif d['_entity.pdbx_description'][i][:5] == 'SUGAR': ## eg 2c49
##                if d['_entity.id'][i] not in d['_entity_poly.entity_id']:
##                    print pdb
##                    stop
        elif entity_type == 'non-polymer' and d['_entity.pdbx_description'][i][:5] == 'SUGAR':
            bool_monosaccharide = True
##            ## just a check
##            if d['_entity.pdbx_description'][i][:7] != 'SUGAR (' and pdb not in ['1iuc',]:
##                print pdb
##                print d['_entity.pdbx_description'][i]
##                stop
##        ## anything else named SUGAR? just a check
##        elif entity_type != 'non-polymer' and d['_entity.pdbx_description'][i][:5] == 'SUGAR':
##            print d
##            print pdb
##            print entity_type
##            print d['_entity.pdbx_description'][i]
##            stop

    if bool_monosaccharide == False and bool_polysaccharide == True and count_polymer_sugar == 1:
        bool_append = True
##    elif pdb in ['3gvj','3gvk','3gvl','3hmy','3msg','1v0f',]:
##        bool_append = False
##    ## error check
##    elif bool_polysaccharide == False and count_polymer_sugar > 0:
##        print d
##        print bool_polysaccharide
##        print d['_entity.pdbx_description']
##        print count_polymer_sugar
##        print pdb
##        stop_no_poly_but_poly

    if pdb == '1dl2':
        print count_polymer_sugar
        print bool_append
        stop

    return bool_append


def one_polypeptide(pdb,):

    l_data_categories = ['_entity_poly',]
    d = parse_mmCIF.main(
        pdb,
        l_data_categories = l_data_categories,
        )

    bool_append = False

    ## make sure polymer is present (not vacomycin 1aa5)
    if '_entity_poly.type' in d.keys():
        ## one polypeptide?
        if d['_entity_poly.type'].count('polypeptide(L)') == 1:
            bool_append = True
##            if not ',' in ''.join(d['_entity_poly.pdbx_strand_id']):
##                bool_append = True
##            list_entity.pdbx_number_of_molecules__1.txt
    
    return bool_append
    

def modres_not_MSE(pdb,):

    l_data_categories = ['_pdbx_struct_mod_residue']
    d = parse_mmCIF.main(
        pdb,
        l_data_categories = l_data_categories,
        )

    bool_append = False

    ## has MODRES
    if '_pdbx_struct_mod_residue.id' in d.keys():
        if d['_pdbx_struct_mod_residue.label_comp_id'] != d['_pdbx_struct_mod_residue.auth_comp_id']:
            print pdb
            stop
        ## at least one MODRES is different from MSE
        if d['_pdbx_struct_mod_residue.auth_comp_id'] != len(d['_pdbx_struct_mod_residue.auth_comp_id'])*['MSE']:
            bool_append = True
    
    return bool_append


def loop(method,suffix,):

    import os, sys

    l_pdbs = []

    l_dn = os.listdir(path)
    l_dn.sort()
    for dn in l_dn:
        if not os.path.isdir('%s/%s' %(path,dn,)):
            continue
        if dn < sys.argv[-1]:
            continue
        print dn
        l_fn = os.listdir('%s/%s' %(path,dn,))
        l_fn.sort()
        for fn in l_fn:
            if fn[-3:] == '.gz':
                continue
            pdb = fn[:4]
##            if pdb != '3aav':
##                continue
            if method == 'modres_not_MSE':
                bool_append = modres_not_MSE(pdb,)
            elif method == 'one_polypeptide':
                bool_append = one_polypeptide(pdb,)
            elif method == 'one_polysaccharide':
                bool_append = one_polysaccharide(pdb,)
            elif method == 'non-polymer':
                bool_append = non_polymer(pdb,)
            if bool_append == True:
                l_pdbs += [pdb]

##            print bool_append
##            stop

    fd = open('%s/list%s.txt' %(path,suffix,),'w')
    fd.writelines('\n'.join(l_pdbs))
    fd.close()

    return l_pdbs


def cluster_by_sequence():

    import os, time

    if os.path.isfile('pdb_seqres.txt'):
        if (time.time()-os.path.getmtime('pdb_seqres.txt'))/(60*60*24) > 14:
            stop
            os.system('wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt')
    else:
        stop

    ## filter out low complexity sequences (XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX)
    ## http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#KarlinAltschul
    ## http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#LCR
    fd = open('pdb_seqres.txt','r')
    lines = fd.readlines()
    fd.close()
    lines_filtered = []
    for i in range(0,len(lines),2):
        if i % 100000 == 0:
            print i
##        if i < 688:
##            continue
##        if i > 4150:
##            break
        ## exclude DNA/RNA
        mol = lines[i].split()[1][4:]
        if mol == 'na':
            continue
        ## exclude short sequences
        length = int(lines[i].split()[2][7:])
        if length < 50:
            continue

        ## exclude low complexity sequences
        ## http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#KarlinAltschul
        ## http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#LCR
        seq = lines[i+1][:-1]
        if seq == length*'X':
            continue

        ## http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/demo/blastclust.c
        ## if (search->prog_number == blast_type_blastp)
        ##    id1 = SeqId2OrdinalId(search->rdfp, search->query_id);
        ## else 
        ##    id1 = -1;        
        ## if (id1 < id2) {
        ## } else {
        ## [blastclust] FATAL ERROR: Blastclust cannot process input files with non-unique sequence identifiers
        ID = lines[i][:7]
        lines[i] = '%s%s%06i%s' %(lines[i][:5],lines[i][6],i,lines[i][7:])

        ## append
        lines_filtered += lines[i:i+2]

    fd = open('pdb_seqres_filtered.txt','w')
    fd.writelines(lines_filtered)
    fd.close()

    if (time.time()-os.path.getmtime('pdb_seqres_filtered.out'))/(60*60*24) > 28:
        stop
        ## http://www.ncbi.nlm.nih.gov/Web/Newsltr/Spring04/blastlab.html
        ## To create a stringent non-redundant protein sequence set, use the following command line:
        ## blastclust -i infile -o outfile -p T -L 1 -b T -S 100
        ## -p(rotein) T(rue)
        ## -L(ength) 1(00%)
        ## -b 
        s = '/home/tc/Downloads/blast-2.2.25/bin/blastclust -i pdb_seqres_filtered.txt -o pdb_seqres_filtered.out -p T -L 1 -b T -S 100'
        os.system(s)

    fd = open('pdb_seqres_filtered.out','r')
    lines = fd.readlines()
    fd.close()
    for i in range(len(lines)):
        lines[i] = ' '.join([ID[:5] for ID in lines[i].split()])+'\n'
    fd = open('bc-100.out','w')
    fd.writelines(lines)
    fd.close()

    return


def unobs_nonterminal_atoms_alpha():

    ## this method is not entirely correct... e.g. 1kwr...

    category = fn = '_pdbx_unobs_or_zero_occ_atoms'

    fd = open('%s/list%s.txt' %(path,fn))
    s = fd.read()
    fd.close()
    l_pdbs_include = s.split()

    ## if a whole residue is missing, then all of it's atoms are also missing
    fd = open('%s/list_pdbx_unobs_residues__NONTERMINAL.txt' %(path))
    s = fd.read()
    fd.close()
    l_pdbs_exclude = s.split()

    l_data_categories = [
        '_pdbx_poly_seq_scheme',
        '_pdbx_unobs_or_zero_occ_atoms',
        '_entity_poly',
        '_struct', ## .pdbx_model_type_details
        '_exptl',
        ]
    d_breaks = {'_exptl.method':['SOLUTION NMR','SOLID-STATE NMR']}

    fn_out = 'list_pdbx_unobs_atoms__CA.txt'

    l_pdbs_out = []
    for pdb in l_pdbs_include:

##        if pdb[1:3] < 'fe':
##            continue
##        if pdb == '2kzt': ## takes too long...
##            continue
        if pdb != '3e3d':
            continue

        if pdb in l_pdbs_exclude:
            continue

        print pdb

        d = parse_mmCIF.main(pdb,l_data_categories=l_data_categories,d_breaks=d_breaks,)

        ## something has to be missing in the first place for it to be terminal/nonterminal
        if not category in d.keys():
            continue
        ## it has to be a polymer in the first place for anything to be terminal/nonterminal
        if not '_pdbx_poly_seq_scheme' in d.keys():
            continue
        ## don't deal with NMR models for now... (too many unobs records when hydrogen...)
        if d['_exptl.method'] != ['X-RAY DIFFRACTION']:
            continue
        if '_struct.pdbx_model_type_details' in d.keys():
            if d['_struct.pdbx_model_type_details'] in [
                ['?'],
                ['minimized average'],
                ['MINIMIZED AVERAGE'],
                ]:
                pass
            ## if residues are not missing, and model is CA only, then no CA are missing!!!
            elif 'CA ATOMS ONLY' in d['_struct.pdbx_model_type_details'][0]:
                continue
            else:
                print d['_struct.pdbx_model_type_details']
                stop
##        if not 'CA' in d['_pdbx_unobs_or_zero_occ_atoms.auth_atom_id']:
##            continue

        bool_append = False
        for i_unobs in range(len(d['_pdbx_unobs_or_zero_occ_atoms'])):
            if (
                d['_pdbx_unobs_or_zero_occ_atoms.auth_atom_id'][i_unobs] == 'CA'
                and
                d['_pdbx_unobs_or_zero_occ_atoms.polymer_flag'][i_unobs] == 'Y'
                and
                ## unobs (1), zero_occ (0)
                d['_pdbx_unobs_or_zero_occ_atoms.occupancy_flag'][i_unobs] == '1'
                ):
                l_pdbs_out += [pdb]
                print '***', pdb
                break

        continue

    print l_pdbs_out
    stop
    fd = open('%s/%s' %(path,fn_out,),'w')
    fd.write('\n'.join(l_pdbs_out))
    fd.close()

    for x in []:

        l_indexes_unobs = []
        bool_append = False
        s = ''.join(d['_pdbx_poly_seq_scheme.pdb_strand_id'])
        for i_unobs in range(len(d['_pdbx_unobs_or_zero_occ_atoms'])):

            ## skip if not alpha carbon
            if d['_pdbx_unobs_or_zero_occ_atoms.auth_atom_id'][i_unobs] != 'CA':
                continue
            ## skip if zero occupancy
            if d['_pdbx_unobs_or_zero_occ_atoms.occupancy_flag'][i_unobs] == '0':
                continue

            if 'HA' in d['_pdbx_unobs_or_zero_occ_atoms.auth_atom_id']:
                print pdb
                print len(d['_pdbx_unobs_or_zero_occ_atoms.auth_seq_id'])
                stop2
            if d['_pdbx_unobs_or_zero_occ_atoms.polymer_flag'].count('Y') > 800:
                print pdb
                print len(d['_pdbx_unobs_or_zero_occ_atoms.auth_seq_id'])
                stop1

            asymID_unobs = d['_pdbx_unobs_or_zero_occ_atoms.auth_asym_id'][i_unobs]
            seqID_unobs = d['_pdbx_unobs_or_zero_occ_atoms.auth_seq_id'][i_unobs]

            index1 = s.index(asymID_unobs)
            index2 = s.rindex(asymID_unobs)+1
            for i_poly in range(index1,index2,):

                asymID_poly = d['_pdbx_poly_seq_scheme.pdb_strand_id'][i_poly]
                seqID_poly = d['_pdbx_poly_seq_scheme.auth_seq_num'][i_poly]

                if seqID_poly == seqID_unobs:

                    if d['_pdbx_poly_seq_scheme.pdb_ins_code'][i_poly] == '.' and d['_pdbx_unobs_or_zero_occ_atoms.PDB_ins_code'][i_unobs] == '?':
                        pass
                    elif d['_pdbx_poly_seq_scheme.pdb_ins_code'][i_poly] == d['_pdbx_unobs_or_zero_occ_atoms.PDB_ins_code'][i_unobs]:
                        pass
                    elif d['_pdbx_unobs_or_zero_occ_atoms.PDB_ins_code'][i_unobs] == '?' and d['_pdbx_poly_seq_scheme.pdb_ins_code'][i_poly] != '.':
                        continue
                    elif not d['_pdbx_unobs_or_zero_occ_atoms.PDB_ins_code'][i_unobs] in d['_pdbx_poly_seq_scheme.pdb_ins_code']:
                        print d['_pdbx_poly_seq_scheme.pdb_ins_code'][i_poly]
                        print insCode_unobs
                        print pdb
                        print seqID_unobs, asymID_unobs
                        stop
                    else:
                        continue

                    if asymID_unobs != asymID_poly:
                        stop_add_with_check_of_identiiical_seqID

                    ## tmp!!! check!!!
                    if d['_pdbx_unobs_or_zero_occ_atoms.auth_comp_id'][i_unobs] != d['_pdbx_poly_seq_scheme.pdb_mon_id'][i_poly]:
                        print pdb
                        stop

##                    ## last residue
##                    if index2-i_poly == 0:
##                        pass ## should append...
##                    ## first residue
##                    elif i_poly-index1 == 0:
##                        pass ## should append...
####                    elif i_poly-index1 > 1 and bool_unobs_prev == False:
####                        bool_append = True
                    ## previous residues are missing
                    elif d['_pdbx_poly_seq_scheme.auth_seq_num'][index1:i_poly] == (i_poly-index1)*['?']:
                        bool_append = True
                    ## next residues are missing
                    elif d['_pdbx_poly_seq_scheme.auth_seq_num'][i_poly+1:index2] == (index2-i_poly-1)*['?']:
                        bool_append = True
                    ## zero occupancy residue prior to residue with unobserved atom(s)
                    elif pdb in ['7adh']:
                        bool_append = False
                        pass
                    else:
                        if len( set(range(index1,i_poly)) - set(l_indexes_unobs) ) == 0:
                            l_indexes_unobs += [i_poly]
                            stop1
                            pass
                        elif len( set(range(i_poly+1,index2)) - set(l_indexes_unobs) ) == 0:
                            l_indexes_unobs += [i_poly]
                            print pdb
                            print l_indexes_unobs
                            print i_poly, index1, index2
                            stop2
                            pass
                        else:
                            ## this method is not entirely correct... e.g. 1kwr...
                            if i_poly-index1 < 10 or index2-i_poly < 10:
                                print pdb
                                print i_poly-index1
##                        print index2-i_poly
                                print seqID_unobs
                                print pdb
                                print d['_pdbx_poly_seq_scheme.auth_seq_num'][index1:i_poly]
                                print d['_pdbx_poly_seq_scheme.auth_seq_num'][i_poly:index2]
                                print pdb
##                                stop
                            bool_append = True
                            break

            if bool_append == True:
                break

        if bool_append == True:
            print pdb
            l_pdbs_out += [pdb]
            continue

        if l_indexes_unobs != []:
            print l_indexes_unobs
            stop

    fd = open('%s/%s' %(path,fn_out,),'w')
    fd.write('\n'.join(l_pdbs_out))
    fd.close()

    return


def unobs_nonterminal_residues():

    ##
    ## unobs or zero occup not at terminals!!! (combination...)
    ## eg dont exlude 200l w 163,164 missing
    ## dont exclude 201l w 163,164 missing, but internally in _pdbx_poly_seq_scheme because 2 chains
    ##
    category = fn = '_pdbx_unobs_or_zero_occ_residues'
    fd = open('%s/list%s.txt' %(path,fn))
    s = fd.read()
    fd.close()
    l_pdbs_in = s.split()
    l_data_categories = [
        '_pdbx_poly_seq_scheme',
        '_pdbx_unobs_or_zero_occ_residues',
        '_entity_poly',
        ]

    fn_out = 'list_pdbx_unobs_residues__NONTERMINAL'

    loop_residues(category,fn_out,)

    l_pdbs_out = []
    for pdb in l_pdbs_in:

##        if pdb[1:3] < 'oa':
##            continue
##        if pdb != '2hub':
##            continue

        ## no residues are present! (e.g. 1oax, 1oay)
        if pdb in ['1oax','1oay',]:
            continue

        d = parse_mmCIF.main(pdb,l_data_categories=l_data_categories,)

##        print pdb

        if not category in d.keys():
            continue

        bool_append = False
        s = ''.join(d['_pdbx_poly_seq_scheme.pdb_strand_id'])
        for chains in d['_entity_poly.pdbx_strand_id']:
            for chain in chains.split(','):
                index1 = s.index(chain)
                index2 = s.rindex(chain)
##                print chain
                l_auth_seq_num = d['_pdbx_poly_seq_scheme.auth_seq_num'][index1:index2+1]
                while l_auth_seq_num[0] == '?':
                    l_auth_seq_num = l_auth_seq_num[1:]
                while l_auth_seq_num[-1] == '?':
                    l_auth_seq_num = l_auth_seq_num[:-1]
                ## non-terminal residues missing?
                if '?' in l_auth_seq_num:
                    print '****', pdb
                    bool_append = True
                    break
            if bool_append == True:
                break
        if bool_append == True:
            print pdb
            l_pdbs_out += [pdb]
            ## continue

    fd = open('%s/%s' %(path,fn_out,),'w')
    fd.write('\n'.join(l_pdbs_out))
    fd.close()

    return


if __name__ == '__main__':
    main()

## change this script to look for inclusions and especially exclusions relative to a prior run - what caused the change?!!!

import os, parse_mmCIF

path = '/media/Elements/mmCIF'
path = '/media/Tommy/mmCIF'


l_data_categories_skip = ['_atom_site','_atom_site_anisotrop']

l_data_categories = []
d_breaks = None

d_include = None
d_exclude_subset = None

##
bool_bool = False

####
#### _pdbx_struct_mod_residue category present
####
##suffix = '_pdbx_struct_mod_residue'
##l_data_categories = ['_pdbx_struct_mod_residue']


######
###### entity = non-polymer
######
##suffix = '_entity.type__non-polymer'
##d_include = {
##    '_entity.type':'non-polymer',
##    }

######
######
######
##suffix = '_chem_comp.type__saccharide'
##d_include = {
##    'chem_comp.type':'saccharide',
##    'chem_comp.type':'D-saccharide',
##    'chem_comp.type':'L-saccharide',
##    }

####
#### oligomeric_details = monomeric
####
##suffix = '_pdbx_struct_assembly.oligomeric_details__monomeric'
##d_include = {'_pdbx_struct_assembly.oligomeric_details':'monomeric',#}
####    '_entity.type':'polymer',
##    }


####
#### _pdbx_unobs_or_zero_occ_residues.id category present
####
##suffix = '_pdbx_unobs_or_zero_occ_residues'
##l_data_categories = ['_pdbx_unobs_or_zero_occ_residues']
##bool_bool = True

####
#### _pdbx_unobs_or_zero_occ_atoms.id category present
####
##suffix = '_pdbx_unobs_or_zero_occ_atoms'
##l_data_categories = ['_pdbx_unobs_or_zero_occ_atoms']
##bool_bool = True

####
####
####
##suffix = '_entity_poly.type__polyribonucleotide'
##d_include = {'_entity_poly.type':'polyribonucleotide',}

####
####
####
##suffix = '_entity_poly.type__polydeoxyribonucleotide'
##d_include = {'_entity_poly.type':'polydeoxyribonucleotide',}

####
####
####
##suffix = '_entity_poly.type__polysaccharide(D)'
##d_include = {'_entity_poly.type':'polysaccharide(D)',}

####
#### '_entity_poly.type' = 'polypeptide(L)'
####
##suffix = '_entity_poly.type__polypeptide(L)'
##d_include = {'_entity_poly.type':'polypeptide(L)'}

####
####
####
##suffix = '_exptl.method__X-RAY_DIFFRACTION'
##d_include = {'_exptl.method':'X-RAY DIFFRACTION',}
##d_breaks = {'_exptl.method':['SOLUTION NMR','SOLID-STATE NMR']}

####
####
####
##suffix = '_exptl.method__SOLUTION_NMR'
##d_include = {'_exptl.method':['SOLUTION_NMR',]}
##d_breaks = {'_exptl.method':['SOLUTION NMR','SOLID-STATE NMR']}


if d_include:
    for item_include,value_include in d_include.items():
        l_data_categories += [item_include[:item_include.index('.')]]
if d_exclude_subset:
    for item_exclude,value_exclude in d_exclude_subset.items():
        l_data_categories += [item_exclude[:item_exclude.index('.')]]


l = []

l_dn = os.listdir(path)
l_dn.sort()
for dn in l_dn:
##    if dn != 'hu':
##        continue
    if not os.path.isdir('%s/%s' %(path,dn,)):
        continue
    print dn
    l_fn = os.listdir('%s/%s' %(path,dn,))
    l_fn.sort()
    for fn in l_fn:
        pdb = fn[:4]
        if fn[-3:] == '.gz':
            continue
########        if pdb in ['2fl9','3gau','3gav','3gaw',]: ## tmp!!!
########            continue
##        print pdb
        fd = open('%s/%s/%s' %(path,dn,fn), 'r')
        lines = fd.readlines()
        fd.close()
        d = parse_mmCIF.main(
            pdb,lines,
            l_data_categories = l_data_categories,
            d_breaks = d_breaks,
            )

        if d_exclude_subset:
            bool_continue = False
            for item_exclude,l_values_exclude in d_exclude_subset.items():
                if not item_exclude in d.keys():
                    bool_continue = True
                    fd = open('%s/remediation_%s.txt' %(path,item_exclude,),'a')
                    fd.write('%s\n' %(pdb))
                    fd.close()
                    continue
                if len( set(d[item_exclude]) & set(l_values_exclude) ) > 0:
                    bool_continue = True
                    break
            if bool_continue == True:
                continue

##        if d_exclude_equal:
##            bool_continue = False
##            for item_exclude,value_exclude in d_exclude_equal.items():
##                if not item_exclude in d.keys():
##                    bool_continue = True
##                    fd = open('%s/expected_but_missing_%s.txt' %(path,item_exclude,),'a')
##                    fd.write('%s\n' %(pdb))
##                    fd.close()
##                    continue
##                if d[item_exclude] == value_exclude:
##                    bool_continue = True
##                    stopstop
##                    break
##            if bool_continue == True:
##                print pdb
##                stop
##                continue

        if d_include:
            bool_continue = False
            for item_include,value_include in d_include.items():
                if not item_include in d.keys():
                    bool_continue = True
                    fd = open('%s/remediation_%s.txt' %(path,item_include,),'a')
                    fd.write('%s\n' %(pdb))
                    fd.close()
                    continue
                if not value_include in d[item_include]:
                    bool_continue = True
                    break
            if bool_continue == True:
                continue

        if bool_bool == True:
            for data_category in l_data_categories:
                if data_category in d.keys():
                    l += ['%s\n' %(pdb)]
                    break
        else:
            l += ['%s\n' %(pdb)]

fd = open('%s/list%s.txt' %(path,suffix,), 'w')
fd.writelines(l)
fd.close()

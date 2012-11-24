import os, numpy, math
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import parse_mmCIF
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import gnuplot

def main():

    d_phipsi_res, d_phipsi_ss = parse_dihedrals()

    d_count = count_dihedrals(d_phipsi_res,d_phipsi_ss,)

    plot_dihedrals(d_count)

    print 'prePro (not Gly)', len(d_phipsi_res['prePRO_notGLY'])
    print 'prePro', len(d_phipsi_res['prePRO'])

    return


def plot_dihedrals(d_count,):

    phipsi_step = 5
    d_titles = {
        'sheet':'{/Symbol b} sheet',
        'helix_alpha':'{/Symbol a} helix',
        'helix_pi':'{/Symbol p} helix',
        'helix_310':'3_1_0 helix',
        'prePRO_notGLY':'prePro (excl. Gly)',
        'prePRO_GLY':'Gly (prePro)',
        'all_notgly_notpro_notprepro':'all (excl. Gly, Pro and pre-Pro)',
        'turns_notgly_notpro_notprepro':'Turn (excl. Gly, Pro and pre-Pro)',
        }

    for k in d_count.keys():
        print 'plot', k
        l = []
        for phi in range(-180,180,phipsi_step,):
            for psi in range(-180,180,phipsi_step,):
                l += ['%s %s %s\n' %(phi,psi,d_count[k][phi][psi])]
            l += ['\n']

        if k in d_titles.keys():
            title = d_titles[k]
        else:
            if k in ['ALA','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
##            if k not in [
##                'Turn','prePRO_notGLY','prePRO_GLY',
##                'GLY', ## not prePro
##                'PRO', ## not prePro
##                'cisPro','transPro',
##                ]:
                continue
            title = k
        
        gnuplot.contour_plot(
            k, l,
            title=title, xlabel='{/Symbol f}', ylabel='{/Symbol y}',
            x1 = -180, x2 = 180,
            y1 = -180, y2 = 180,
            z1 = 0,
            bool_remove = False,
            )

    return


def round_angle(angle,phipsi_step,):

##    angle = phipsi_step*int(angle/phipsi_step)
##    if angle == 180.:
##        angle = -180.
    if angle == 180.:
        angle = -180.
    else:
        angle = phipsi_step*int(angle/phipsi_step)

    return angle


def count_dihedrals(d_phipsi_res,d_phipsi_ss,):

    phipsi_step = 5

    d_count = {}
    for d in [d_phipsi_res,d_phipsi_ss,]:
        for k in d.keys():
            d_count[k] = {}
            for phi in range(-180,180):
                d_count[k][phi] = {}
                for psi in range(-180,180):
                    d_count[k][phi][psi] = 0
            for phi,psi in d[k]:
                phi = round_angle(phi,phipsi_step,)
                psi = round_angle(psi,phipsi_step,)
                d_count[k][phi][psi] += 1

    return d_count


def parse_dihedrals():

    import sys

    path = '/data/mmCIF'

    d_phipsi_res = {
        'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],
        'GLY':[],'HIS':[],'ILE':[],'LYS':[],'LEU':[],
        'MET':[],'ASN':[],'PRO':[],'GLN':[],'ARG':[],
        'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[],
        'prePRO':[],'prePRO_notGLY':[],'prePRO_GLY':[],
        'cisPro':[],'transPro':[],
        'all_notgly_notpro_notprepro':[],
        }

    d_phipsi_ss = {
        'sheet':[], ## _struct_sheet_order.sense
        ##_struct_conf.pdbx_PDB_helix_class
        'helix_alpha':[], ## i+4 # 1
        'helix_pi':[], ## i+5 # 3
        'helix_310':[], ## i+3 # 5
        'Turn':[], ## i+?
        ##
        'turns_notgly_notpro_notprepro':[],
        }

    d_counts = {
        'cisProALA':0,
        'cisProCYS':0,
        'cisProASP':0,
        'cisProGLU':0,
        'cisProPHE':0,
        'cisProGLY':0,
        'cisProHIS':0,
        'cisProILE':0,
        'cisProLYS':0,
        'cisProLEU':0,
        'cisProMET':0,
        'cisProASN':0,
        'cisProPRO':0,
        'cisProGLN':0,
        'cisProARG':0,
        'cisProSER':0,
        'cisProTHR':0,
        'cisProVAL':0,
        'cisProTRP':0,
        'cisProTYR':0,
        'cisPro_helix':0,
        'cisPro_sheet':0,
        'cisPro_turn':0,
        'cisPro_random':0,
        }

    l_dn = os.listdir(path)
    l_dn.sort()
    l_dn.remove('mmCIF.py')
    for dn in l_dn:
        if dn < sys.argv[-2]:
            continue
        if dn > sys.argv[-1]:
            continue
        print '*',dn
        l_fn = os.listdir('%s/%s' %(path,dn,))
        l_fn.sort()
        for fn in l_fn:
            pdb = fn[:4]
            print pdb
            d_mmCIF = parse_mmCIF.main(
                pdb,
                d_breaks = {'_exptl.method':['SOLUTION NMR']},
                l_data_categories = [
                    '_exptl',
                    '_refine',

                    '_struct_conf', ## HELIX
                    '_struct_sheet_range', ## SHEET

                    '_entity',
                    '_entity_poly',
                    '_entity_poly_seq',

                    '_atom_site',
                    ],
                )

            ## skip NMR models
            if ''.join(d_mmCIF['_exptl.method']) in [
                'SOLUTION NMR',
                'POWDER DIFFRACTION',
                'ELECTRON MICROSCOPY',
                ]:
                continue

            if not '_refine.ls_d_res_high' in d_mmCIF.keys():
                print d_mmCIF['_exptl.method']
                continue

            ## skip if multiple resolutions
            if len(d_mmCIF['_refine.ls_d_res_high']) > 1:
                continue

            ## skip if no resolution
            if ''.join(d_mmCIF['_refine.ls_d_res_high']) == '?':
                continue

            ## skip low resolution structures
            if float(''.join(d_mmCIF['_refine.ls_d_res_high'])) > 2:
                continue

            if not 'polymer' in d_mmCIF['_entity.type']:
                continue
            if not '_entity_poly.type' in d_mmCIF.keys(): ## e.g. 1hhu
                continue
            if d_mmCIF['_entity_poly.type'] == ['polydeoxyribonucleotide/polyribonucleotide hybrid']:
                continue
            if d_mmCIF['_entity_poly.type'] == ['polydeoxyribonucleotide']:
                continue

            d_sequence = {}
            for i_entity_poly_seq in range(len(d_mmCIF['_entity_poly_seq.entity_id'])):
                entity_id = int(d_mmCIF['_entity_poly_seq.entity_id'][i_entity_poly_seq])
                if not entity_id in d_sequence.keys():
                    d_sequence[entity_id] = []
                res_no = int(d_mmCIF['_entity_poly_seq.num'][i_entity_poly_seq])
                res_name = d_mmCIF['_entity_poly_seq.mon_id'][i_entity_poly_seq]
                d_sequence[entity_id] += [{'res_no':res_no,'res_name':res_name,}]

            l_entities_poly = []
            for i_entity_poly in range(len(d_mmCIF['_entity_poly.entity_id'])):
                ## skip if not polypeptide
                entity_poly_type = d_mmCIF['_entity_poly.type'][i_entity_poly]
                if entity_poly_type != 'polypeptide(L)':
                    continue
                ## skip if nonstd linkages
                if d_mmCIF['_entity_poly.nstd_linkage'][i_entity_poly] == 'yes':
                    print pdb
                    stop
                    continue
                ## parse entity_id and chains
                entity_id = int(d_mmCIF['_entity_poly.entity_id'][i_entity_poly])
                l_entities_poly += [entity_id]
            ## skip if no polypeptide chains
            if l_entities_poly == []:
                continue

            d_coords = {}
            for i_atom_site in range(len(d_mmCIF['_atom_site.id'])):

                entity_id = int(d_mmCIF['_atom_site.label_entity_id'][i_atom_site])
                ## not a polymer
                if not entity_id in l_entities_poly:
                    continue
                ## polymer, append
                elif not entity_id in d_coords.keys():
                    d_coords[entity_id] = {}

                model = int(d_mmCIF['_atom_site.pdbx_PDB_model_num'][i_atom_site])
                if model > 1:
                    continue

                chain = d_mmCIF['_atom_site.label_asym_id'][i_atom_site]
                if not chain in d_coords[entity_id].keys():
                    d_coords[entity_id][chain] = {}
                res_no = int(d_mmCIF['_atom_site.label_seq_id'][i_atom_site])
                if not res_no in d_coords[entity_id][chain].keys():
                    d_coords[entity_id][chain][res_no] = {}
                atom_name = d_mmCIF['_atom_site.label_atom_id'][i_atom_site]

                altloc = d_mmCIF['_atom_site.label_alt_id'][i_atom_site]
                if altloc not in ['.','A','1',]:
                    continue

                ## skip if zero occupancy
                occupancy = float(d_mmCIF['_atom_site.occupancy'][i_atom_site])
                if altloc == '.' and occupancy == 0:
                    continue

                if atom_name in ['CA','C','O','N',] and atom_name in d_coords[entity_id][chain][res_no].keys():
                    print pdb, chain, res_no, atom_name
                    print d_mmCIF['_atom_site.Cartn_x'][i_atom_site], d_mmCIF['_atom_site.Cartn_y'][i_atom_site]
                    print d_coords[entity_id][chain][res_no][atom_name]
                    stop
                x = float(d_mmCIF['_atom_site.Cartn_x'][i_atom_site])
                y = float(d_mmCIF['_atom_site.Cartn_y'][i_atom_site])
                z = float(d_mmCIF['_atom_site.Cartn_z'][i_atom_site])
                coord = numpy.array([x,y,z,])
                d_coords[entity_id][chain][res_no][atom_name] = coord

            d_helices = {}
            ## helices or turns present?
            if '_struct_conf.id' in d_mmCIF.keys():
                for i_struct_conf in range(len(d_mmCIF['_struct_conf.id'])):
                    chain1 = d_mmCIF['_struct_conf.beg_label_asym_id'][i_struct_conf]
                    chain2 = d_mmCIF['_struct_conf.end_label_asym_id'][i_struct_conf]
                    res_no1 = int(d_mmCIF['_struct_conf.beg_label_seq_id'][i_struct_conf])
                    res_no2 = int(d_mmCIF['_struct_conf.end_label_seq_id'][i_struct_conf])
                    conf_type_id = d_mmCIF['_struct_conf.conf_type_id'][i_struct_conf]
                    if chain1 != chain2:
                        print chain1, chain2, pdb
                        stop
                    if conf_type_id == 'HELX_P':
                        helix_class = int(d_mmCIF['_struct_conf.pdbx_PDB_helix_class'][i_struct_conf])
                    elif conf_type_id == 'TURN_P':
                        helix_class = 99
                    else:
                        print conf_type_id
                        print pdb
                        stop
                    l_res_nos = range(res_no1,res_no2+1,)
                    if not chain1 in d_helices.keys():
                        d_helices[chain1] = {}
                    for res_no in l_res_nos:
                        d_helices[chain1][res_no] = helix_class

            d_sheets = {}
            ## sheet present?
            if '_struct_sheet_range.sheet_id' in d_mmCIF.keys():
                for i_struct_sheet_range in range(len(d_mmCIF['_struct_sheet_range.sheet_id'])):
                    chain1 = d_mmCIF['_struct_sheet_range.beg_label_asym_id'][i_struct_sheet_range]
                    chain2 = d_mmCIF['_struct_sheet_range.end_label_asym_id'][i_struct_sheet_range]
                    res_no1 = int(d_mmCIF['_struct_sheet_range.beg_label_seq_id'][i_struct_sheet_range])
                    res_no2 = int(d_mmCIF['_struct_sheet_range.end_label_seq_id'][i_struct_sheet_range])
                    l_res_nos = range(res_no1,res_no2+1,)
                    if chain1 != chain2:
                        print chain1, chain2, pdb
                        stop
                    if not chain1 in d_sheets.keys():
                        d_sheets[chain1] = []
                    for res_no in l_res_nos:
                        d_sheets[chain1] += l_res_nos

            for entity_id in l_entities_poly:
                for chain in d_coords[entity_id].keys():
                    ## skip if short peptide (e.g. 13gs)
                    if len(d_sequence[entity_id]) <= 3:
                        continue
                    for i_res_no in range(1,len(d_sequence[entity_id])-1):
                        res_no_prev = int(d_sequence[entity_id][i_res_no-1]['res_no'])
                        res_no = int(d_sequence[entity_id][i_res_no]['res_no'])
                        res_no_next = int(d_sequence[entity_id][i_res_no+1]['res_no'])
                        res_name = d_sequence[entity_id][i_res_no]['res_name']
                        if res_name == 'MSE':
                            res_name = 'MET'
                        res_name_next = d_sequence[entity_id][i_res_no+1]['res_name']

                        ## not a standard residue
                        if not res_name in d_phipsi_res.keys():
                            continue

                        ## residue not observed
                        if not res_no_prev in d_coords[entity_id][chain].keys():
                            continue
                        if not res_no in d_coords[entity_id][chain].keys():
                            continue
                        if not res_no_next in d_coords[entity_id][chain].keys():
                            continue

                        ## atom not observed
                        if not 'C' in d_coords[entity_id][chain][res_no_prev]:
                            continue
                        if not 'N' in d_coords[entity_id][chain][res_no]:
                            continue
                        if not 'CA' in d_coords[entity_id][chain][res_no]:
                            continue
                        if not 'C' in d_coords[entity_id][chain][res_no]:
                            continue
                        if not 'N' in d_coords[entity_id][chain][res_no_next]:
                            continue
                        
                        C_prev = d_coords[entity_id][chain][res_no_prev]['C']
                        N = d_coords[entity_id][chain][res_no]['N']
                        CA = d_coords[entity_id][chain][res_no]['CA']
                        C = d_coords[entity_id][chain][res_no]['C']
                        N_next = d_coords[entity_id][chain][res_no_next]['N']
                        phi = calc_dihedral(C_prev,N,CA,C,)
                        psi = calc_dihedral(N,CA,C,N_next,)

                        if 'CA' in d_coords[entity_id][chain][res_no_prev].keys():
                            CA_prev = d_coords[entity_id][chain][res_no_prev]['CA']
                            omega = calc_dihedral(CA_prev,C_prev,N,CA,)
                        else:
                            omega = None

                        
                        if omega:
                            if (
                                omega
                                and
                                omega < 150
                                and
                                omega > -150
                                ): ## 12e8, PRO44D
                                if abs(omega) > 30: ## 12e8 PRO196D, 1a44 GLU82A
                                    omega = None
                                ## cis
                                else:
                                    omega = 'cis'
                                    pass
                            ## trans
                            else:
                                omega = 'trans'
                                pass
                        else:
                            omega = None
                        
                        bool_helix = False
                        if chain in d_helices.keys():
                            if res_no in d_helices[chain].keys():
                                bool_helix = True
                                helix_class = d_helices[chain][res_no]

                        bool_sheet = False
                        if chain in d_sheets.keys():
                            if res_no in d_sheets[chain]:
                                bool_sheet = True

##                        if bool_helix == True and bool_sheet == True and helix_class != 99:
##                            print pdb, chain, res_no, 'sheet and helix'
####                            stop
                        
                        if res_name_next == 'PRO':
                            d_phipsi_res['prePRO'] += [[phi,psi,]]
                            if res_name != 'GLY':
                                d_phipsi_res['prePRO_notGLY'] += [[phi,psi,]]
                            else:
                                d_phipsi_res['prePRO_GLY'] += [[phi,psi,]]
                        else:
                            d_phipsi_res[res_name] += [[phi,psi,]]
                            if res_name not in ['GLY','PRO',]:
                                d_phipsi_res['all_notgly_notpro_notprepro'] += [[phi,psi,]]
                            elif res_name == 'PRO' and omega:
                                d_phipsi_res['%sPro' %(omega)] += [[phi,psi,]]
                                if omega == 'cis':
                                    d_counts['cisPro%s' %(res_name)] += 1
                                    if bool_helix == True:
                                        if helix_class == 1:
                                            d_counts['cisPro_helix'] += 1
                                        elif helix_class == 99:
                                            d_counts['cisPro_turn'] += 99
                                    elif bool_sheet == True:
                                        d_counts['cisPro_sheet'] += 1
                                    else:
                                        d_counts['cisPro_random'] += 1
                                        

                        if bool_helix == True:
##                            if helix_class not in [1,3,5,99,]:
##                                print pdb, chain, res_no, helix_class
##                                print 'unexpected helix class'
####                                stop_helix_class
                            if helix_class == 1:
                                d_phipsi_ss['helix_alpha'] += [[phi,psi,]]
                            elif helix_class == 3:
                                d_phipsi_ss['helix_pi'] += [[phi,psi,]]
                            elif helix_class == 5:
                                d_phipsi_ss['helix_310'] += [[phi,psi,]]
                            elif helix_class == 99:
                                d_phipsi_ss['Turn'] += [[phi,psi,]]
                                if (
                                    res_name_next != 'PRO'
                                    and
                                    res_name not in ['GLY','PRO',]
                                    ):
                                    d_phipsi_ss['turns_notgly_notpro_notprepro'] += [[phi,psi,]]
                        if bool_sheet == True:
                            d_phipsi_ss['sheet'] += [[phi,psi,]]

    l = []
    for k in d_counts.keys():
        count = d_counts[k]
        l += ['%s %s\n' %(k,count,)]
    fd = open('count.txt','w')
    fd.writelines(l)
    fd.close()

    return d_phipsi_res, d_phipsi_ss


def calc_dihedral(c1,c2,c3,c4,):

    import numpy,math

    v1 = c2-c1
    v2 = c3-c2
    v3 = c4-c3

    angle = math.atan2(
        numpy.dot(
            math.sqrt(sum(v2*v2))*v1,
            cross(v2,v3),
            ),
        numpy.dot(
            cross(v1,v2),
            cross(v2,v3),
            ),
        )
    angle *= 180./math.pi

    return angle


def cross(v1,v2):

    import numpy

    n = numpy.array([
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0],
        ])

    return n


if __name__ == '__main__':
    main()

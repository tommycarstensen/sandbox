## Tommy Carstensen, UCD, 2009-2010

## imports
import math, numpy
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import parse_mmCIF
sys.path.append('/home/people/tc/svn/Protool/')
import geometry

## instances
instance_geometry = geometry.geometry()

d_symmetry = {
    'GLY':[],
    'ALA':[],
    'ILE':[],
    'CYS':[],
    'LYS':[],
    'MET':[],
    'PRO':[],
    'SER':[],
    'THR':[],
    'TRP':[],
    'ASP':['OD1','OD2',],
    'GLU':['OE1','OE2',],
    'ASN':['OD1','ND2',],
    'GLN':['OE1','NE2',],
    'PHE':['CD1','CD2','CE1','CE2',],
    'TYR':['CD1','CD2','CE1','CE2',],
    'VAL':['CG1','CG2',],
    'LEU':['CD1','CD2',],
    'ARG':['NH1','NH2',],
    'HIS':['CD2','ND1','CE1','NE2',],
    'MSE':[],
    }

def main():

    return


def parse_pdbs_from_cluster(i_cluster,):

##    fd = open('pdbS95bF.out','r')
    fd = open('bc-95.out','r')
    lines = fd.readlines()
    fd.close()
    l_pdbs = lines[i_cluster].split()
##    for i in range(len(l_pdbs)):
##        l_pdbs[i] = l_pdbs[i].lower()
##    l_pdbs = l_pdbs[:50]

    return l_pdbs


def parse_cifs(
    l_pdbs,
    ref_seq, l_db_codes,
    n_mutations_max,
    resolution_min,
    bool_multiple_entities = False,
    ):

    print 'parse cifs'

    n_mutants = 0
    l_wts = []
    l_wts_cysfree = []
    d_mutants = {}

    d_mmCIF_main = {}
    for pdb in l_pdbs:

        if pdb[:4].lower() in d_mmCIF_main.keys():
            continue

        d_mmCIF = parse_mmCIF.main(pdb[:4].lower(),)

        ## not an x-ray structure
        if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
            print pdb, d_mmCIF['_exptl.method']
            continue

        ## more than one type of polymer present
        n_entities = len(d_mmCIF['_entity_poly.entity_id'])
        if bool_multiple_entities == False:
            if n_entities > 1:
                print pdb, 'entities', n_entities #, d_mmCIF['_struct.title']
                continue

        ## low resolution
        if d_mmCIF['_refine.ls_d_res_high'] != d_mmCIF['_refine_hist.d_res_high']:
            print d_mmCIF['_refine.ls_d_res_high']
            print d_mmCIF['_refine_hist.d_res_high']
            stop
        if resolution_min:
##            if float(d_mmCIF['_refine.ls_d_res_high'][0]) >= resolution_min:
            if float(d_mmCIF['_refine.ls_d_res_high'][0]) > resolution_min:
                print pdb, 'resolution', d_mmCIF['_refine.ls_d_res_high']
                continue

        ## get entity ID from chain ID
        for i_entity in range(len(d_mmCIF['_entity_poly.entity_id'])):
            entity_id = d_mmCIF['_entity_poly.entity_id'][i_entity]
            s_chain_ids = d_mmCIF['_entity_poly.pdbx_strand_id'][i_entity]
            if pdb[-1] in s_chain_ids:
                break
        if pdb[-1] not in s_chain_ids:
            print pdb
            print s_chain_ids
            stop
        ## get sequence from entity ID
        seq = []
        for i in range(len(d_mmCIF['_entity_poly_seq.entity_id'])):
            if d_mmCIF['_entity_poly_seq.entity_id'][i] == entity_id:
                mon_id = d_mmCIF['_entity_poly_seq.mon_id'][i]
                if pdb[:4] == '1RCM' and i == 126:
                    if mon_id != 'CYS':
                        stop
                    mon_id = 'CCS'
                seq += [mon_id]

        ## wrong chain length
        if ref_seq:
            if len(seq) != len(ref_seq):
                if ''.join(ref_seq) in ''.join(seq):
                    print ref_seq
                    print seq
                    stop
                ## unobserved atoms not in seqres
                elif ''.join(seq) in ''.join(ref_seq):
                    pass
                ## last two residues unobserved
                elif len(seq) == 162 and pdb in [
                    '1KS3_A','1KW5_A','1KW7_A','1KY0_A','1KY1_A','1L0J_A','1LOK_A','1LPY_A','1LW9_A','1LWG_A','1LWK_A',
                    ]:
                    pass
                ## last two residues unobserved
                elif len(seq) == 162 and seq[-1] == 'LYS':
                    pass
                else:
                    print pdb, 'seqlen', len(seq)
                    continue

        ## not from Gallus gallus
        ## check not necessary, because sequence checked against ref seq
        entity_id = d_mmCIF['_entity_poly.entity_id'][i_entity]
        db_code = d_mmCIF['_struct_ref.db_code'][d_mmCIF['_struct_ref.entity_id'].index(entity_id)]
        if db_code not in l_db_codes:
            print pdb, 'uniprot', db_code
            continue

        ## more than 1 mutation?
        if n_mutations_max != None:
            l_mutations = []
            for i_seq in range(len(seq)):
                res_id_mmCIF = seq[i_seq]
                res_id_uniprot = ref_seq[i_seq]
                if res_id_mmCIF != res_id_uniprot:
                    l_mutations += ['%3s%i%3s' %(res_id_uniprot,i_seq+1,res_id_mmCIF,)]
##            if len(l_mutations) == 1:
            if len(l_mutations) > n_mutations_max:
                print pdb, 'muts', len(l_mutations)
                continue
            elif len(l_mutations) > 0:
                n_mutants += 1
                startmodel = parse_mmCIF_item(d_mmCIF,'_refine.pdbx_starting_model',pdb,)
                    

        ## append to lists and dictionaries
        d_mmCIF_main[pdb[:4]] = d_mmCIF
        if len(l_mutations) > 0:
            if l_mutations == ['CYS54THR', 'CYS97ALA']:
                l_wts_cysfree += [pdb]
            d_mutants[pdb] = {'mutations':l_mutations,'startmodel':startmodel}
        else:
            l_wts += [pdb]

##    print 'd_mutants', d_mutants

    return d_mmCIF_main, l_wts, d_mutants, l_wts_cysfree


def parse_mmCIF_item(d_mmCIF,item,pdb,):

    d_correct_startmodel = {
        ## Japanese Quail
        '2IHL': '3LYM',
        ## BWQL
        '1BQL': '6LYZ', '1DKJ': '6LYZ', '1DKK': '6LYZ',
##        ## HEWL w modres
##        '132L': '6LYZ', ## MLY
        ## 2LZH
##        '1HEL': '2LZH', '1LZT': '2LZH', '2LYZ': '2LZH', '6LYZ': '2LZH', '1LYZ': '2LZH', '5LYZ': '2LZH', '3LYZ': '2LZH', '4LYZ': '2LZH',
        ## 1HEL
        '1HEQ': '1HEL', '1HEP': '1HEL', '1HEO': '1HEL', '1HEN': '1HEL', '1HEM': '1HEL', '1HER': '1HEL',
        '1DPW': '1HEL', '1DPX': '1HEL',
        ## 6LYZ
        '1LZD': '6LYZ', '1LZB': '6LYZ', '1LZE': '6LYZ', '1LZC': '6LYZ', '1LZA': '6LYZ', '1LZG': '6LYZ',
        ## HEWL
        '1LJ3': '1UCO', '1LJ4': '1UCO', '1LJE': '1UCO', '1LJF': '1UCO', '1LJG': '1UCO', '1LJH': '1UCO', '1LJI': '1UCO', '1LJJ': '1UCO', '1LJK': '1UCO',
        ## 1RFP
        '1FLY': '1RFP', '1IOQ': '1RFP', '1FN5': '1RFP', '1IOR': '1RFP', '1IOT': '1RFP', '1IOS': '1RFP', '1FLQ': '1RFP', '1FLW': '1RFP', '1FLU': '1RFP',
        ## other
        '1V7S': '2LZT', '1V7T': '2LZT', '1F10': '1F0W',
        '1WTN': '1BGI', '1WTM': '1BGI',
        '1XEI': '1LMA', '1XEK': '1LMA', '1XEJ': '1XEI',
        ## T4L
        '1KS3': '1L63', '1KW5': '1L63', '1KW7': '1L63', '1KY0': '1L63', '1KY1': '1L63', '1L0J': '1L63', '1L0K': '1L63',
        '1LPY': '1L63', ## cryo structure, but room temperature starting model (not 1lw9)
        '1L35': '2LZM', ## why 2lzm and not 3lzm???
        '1PQJ': '1L63',
        '1L55': '1L63', '1L59': '1L63', '1L61': '1L63', '1L62': '1L63',
        '1L63': '4LZM', '1L57': '4LZM',
        '2LZM': '1LZM',
        '4LZM': '2LZM', '5LZM': '4LZM', '6LZM': '4LZM',
        '257L': '1L63', '258L': '1L63', '259L': '1L63', '260L': '1L63',
        ## 1L63 according to pdb, 4lzm-7lzm according to ref, 1l63 acc to text
        '3F8V': '1L63', ## Asn68 doesn't fit with 3lzm, good fit with 1l63
        '3FA0': '1L63', '3FAD': '1L63', '3F9L': '1L63',
        ## 1L63 according to pdb, 4lzm-7lzm according to ref
        '3C7W': '1L63', '3C7Y': '1L63', '3C7Z': '1L63',
        '3C80': '1L63', '3C81': '1L63', '3C82': '1L63', '3C83': '1L63',
        '3C8Q': '1L63', '3C8R': '1L63', '3C8S': '1L63',
        '3CDQ': '1L63', '3CDR': '1L63', '3CDT': '1L63', '3CDV': '1L63',
        ## T4L guesses
        '1LGU': '2LZM', ## "wild type T4 lysozyme"

        '1D9W': '3LZM', ## no SM, 2lzm doesnt fit with Asp20,Ser38,Leu39
        '226L': '3LZM', ## no SM, 2lzm doesn't fit with Asp20,glu64
        '256L': '3LZM', ## no SM, 2lzm or 3lzm or 7lzm (256l from 1998) 2lzm does not fit with Asp20,Ser38,Leu39

        '1DYA': '3LZM', ## no SM, 2lzm doesn't fit
        '1DYB': '3LZM', ## no SM, 2lzm doesn't fit
        '1DYC': '3LZM', ## no SM, 2lzm doesn't fit
        '1DYD': '3LZM', ## no SM, 2lzm doesn't fit
        '1DYE': '3LZM', ## no SM, 2lzm doesn't fit
        '1DYF': '3LZM', ## no SM, 2lzm doesn't fit
        '1DYG': '3LZM', ## no SM, 2lzm doesn't fit
        '1L02': '3LZM', ## with 1l03-1l15
        '1L03': '3LZM', ## no SM, 2lzm doesn't fit with asp20,leu39
        '1L04': '3LZM', ## with 1l03-1l15
        '1L05': '3LZM', ## with 1l03-1l15
        '1L06': '3LZM', ## with 1l03-1l15
        '1L07': '3LZM', ## with 1l03-1l15
        '1L08': '3LZM', ## with 1l03-1l15
        '1L09': '3LZM', ## with 1l03-1l15
        '1L10': '3LZM', ## 
        '1L11': '3LZM', ## with 1l03-1l15
        '1L12': '3LZM', ## with 1l03-1l15
        '1L13': '3LZM', ## with 1l03-1l15
        '1L14': '3LZM', ## with 1l03-1l15
        '1L15': '3LZM', ## with 1l03-1l15
        '1L16': '3LZM', ## 
        '1L17': '3LZM', ## no SM, 2lzm doesn't fit with asp20,glu64
        '1L18': '3LZM', ## no SM, 2lzm doesn't fit with asp20,glu64
        '1L19': '3LZM', ## no SM, 2lzm doesn't fit with asp20,leu39
        '1L20': '3LZM', ## no SM, 2lzm or 3lzm or 7lzm (256l from 1998) 2lzm does not fit with Asp20,Ser38,Leu39
        '1L21': '3LZM', ## with 1l22,1l33
        '1L22': '3LZM', ## no SM, 2lzm doesn't fit with asp20, lys43
        '1L23': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39, lys43, glu64
        '1L24': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39
        '1L25': '3LZM', ## no SM, 2lzm doesn't fit with arg8,lys135.gln141,lys147
        '1L26': '3LZM', ## with 1l25
        '1L27': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39
        '1L28': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39
        '1L29': '3LZM', ## with 1l25
        '1L30': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39, glu64
        '1L31': '3LZM', ## with 1l25
        '1L32': '3LZM', ## no SM, 2lzm does not fit with Asp20,Leu39
        '1L33': '3LZM', ## no SM, 2lzm does not fit with Asp20,Ser38,Leu39
        '1L34': '3LZM', ## no SM, 2lzm does not fit with Asp20,Leu39,lys43,glu64
        '1L37': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L38': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L42': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L43': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L44': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L45': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39
        '1L46': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L47': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39, lys43
        '1L48': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39, glu64
        '1L52': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39
        '1L53': '3LZM', ## no SM, 2lzm doesn't fit with arg8,met106,glu108,lys135,gln141,lys147
        '1L56': '3LZM', ## no SM, 2lzm doesn't fit with asp20, lys43
        '1L57': '3LZM', ## 
        '1L58': '3LZM', ## no SM, 2lzm or 3lzm or 7lzm (256l from 1998) 2lzm does not fit with Asp20,Ser38,Leu39
        '1L60': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39, lys43, glu64
        '1L69': '3LZM', ## no SM, 2lzm doesn't fit with arg8,ser117,lys135,gln141,lys147
        '1L96': '3LZM', ## no SM, 2lzm doesn't fit with asp20, leu39

##        '1L70': '3LZM', ## with 1l69
##        '1L71': '3LZM', ## with 1l69
##        '1L72': '3LZM', ## with 1l69
##        '1L73': '3LZM', ## with 1l69
##        '1L74': '3LZM', ## with 1l69
##        '1L75': '3LZM', ## with 1l69
        '1L98': '3LZM', ## no SM, 2lzm doesn't fit with asp20
        '1L99': '3LZM', ## no SM, 2lzm doesn't fit with arg8,met106,glu108,met120,lys135,arg145,lys147
        '1L00': '3LZM', ## no SM, 2lzm doesn't fit with arg8,thr115,lys135.gln141
        ## T4L assume error
        '1LLH': '1L63', ## not 2LZM, wt*
        '206L': '1L63', ## not 3LZM, wt*
        }

    d_correct_T = {
        '1LSA':120,'1LSB':180,'1LSC':250,'1LSD':280,'1LSE':295,'1LSF':95,
        '3LYT':100,'5LYT':100,'4LYT':298,'6LYT':298,
        }

    d_correct_pH = {
        '2G4P':4.5,'2G4Q':8.0,
        '1LZB':4.7,'1LZC':4.7,'1LZE':4.7,'1LZG':4.7,
        }

    ## parse nonmandatory mmCIF items

    if item in d_mmCIF.keys():
        l_values = d_mmCIF[item]
        if len(l_values) > 1:
            value = None
        elif l_values == ['?']:
            value = None
        else:
            value = l_values[0]
    else:
        value = None

    if item == '_refine.pdbx_starting_model' and value != None:
        value = value.upper()
        value = value.replace('PDB ENTRY ','')
        value = value.replace('PBD ENTRY ','')
        value = value.replace('PDB CODE ','')
        value = value.replace('.PDB','')
        ## HEWL
        value = value.replace('MONOMER A FROM 5LYM PDB ENTRY (MONOCLINIC LYSOZYME)','5LYM') ## 1LCN
        value = value.replace('MONOMER A FROM 5LYM (MONOCLINIC LYSOZYME)','5LYM')
        value = value.replace('PDB ENTRIES 1LYS AND 1IGM','1LYS')
        value = value.replace('PDB ENTRIES 1LKS AND 1QLE (RESIDUES 1-85)','1LKS')
        value = value.replace('PDB ENTRIES 1MEL AND 1BZQ','1MEL')
        value = value.replace('HULYS STRUCTURE, AND LYSOZYME FROM 1VFB','1VFB')

        ## T4L
##        value = value.replace('WILD TYPE T4 LYSOZYME','2LZM')
        value = value.replace('T4 LYSOZYME 2LZM','2LZM')

        value = value.replace('CYS-FREE WILD TYPE LYSOZYME','1L63')
        value = value.replace('CYS-FREE VERSION OF T4 LYSOZYME','1L63')
        value = value.replace('WT* T4 LYSOZYME','1L63')
        value = value.replace('ISOMORPHOUS WT L163','1L63')
        value = value.replace('L163','1L63')
        value = value.replace('1L63 WITH RESIDUES 106 - 114 DELETED','1L63')

        ## guessing...
        if pdb[:4] == '235L':
            value = value.replace('WILD TYPE','1L63')
        if pdb[:4] == '223L':
            value = value.replace('WILD TYPE T4 LYSOZYME','3LZM')
        if pdb[:4] == '225L':
            value = value.replace('WILD TYPE T4 LYSOZYME','3LZM')

        if pdb[:4] not in d_correct_startmodel.keys() and (len(value) != 4 or value.upper() == 'NONE'):
            if value.upper() != 'NONE':
                print pdb, item, value
            value = None

    ##
    ## corrections
    ##
            
    if item == '_refine.pdbx_starting_model':
        if pdb[:4] in d_correct_startmodel:
            value = d_correct_startmodel[pdb[:4]]

    if item == '_symmetry.space_group_name_H-M' and pdb[:4] == '1LKS':
        value = 'P 1'

    if item == '_diffrn.ambient_temp':
        if pdb[:4] in d_correct_T.keys():
            value = d_correct_T[pdb[:4]]

    if item == '_exptl_crystal_grow.pH':
        if pdb[:4] in d_correct_pH.keys():
            value = d_correct_pH[pdb[:4]]
            

    return value


def parse_coordinates(
    l_pdbs,d_mmCIF_main,method,
    protein,
    bool_exclusion_zero_occupancy = True,
    bool_exclusion_symmetry = False,
    bool_exclusion_altlocs = False,
    bool_exclusion_high_temp_factors = False,
    ):

    print 'parse coordinates', method

    d_coordinates = {}

    for pdb in l_pdbs:

        if not pdb[:4] in d_mmCIF_main.keys():
            continue

##        print 'parse coordinates', pdb

        d_mmCIF = d_mmCIF_main[pdb[:4]]

        d_coords = {}
        for i in range(len(d_mmCIF['_atom_site.id'])):

            ## chain
            if d_mmCIF['_atom_site.auth_asym_id'][i] != pdb[-1].upper(): ## not label_asym_id...
                continue
            ## ligand
            if d_mmCIF['_atom_site.group_PDB'][i] == 'HETATM': ## not used for modres...
                continue

            atom_name = d_mmCIF['_atom_site.label_atom_id'][i]
            occupancy = float(d_mmCIF['_atom_site.occupancy'][i])
            res_no = int(d_mmCIF['_atom_site.label_seq_id'][i])
            res_name = d_mmCIF['_atom_site.label_comp_id'][i]
            x = float(d_mmCIF['_atom_site.Cartn_x'][i])
            y = float(d_mmCIF['_atom_site.Cartn_y'][i])
            z = float(d_mmCIF['_atom_site.Cartn_z'][i])
            tempFactor = float(d_mmCIF['_atom_site.B_iso_or_equiv'][i])
            coord = numpy.array([x,y,z,])

            if method == 'alpha':
                ## alpha carbon
                if atom_name != 'CA':
                    continue
            elif method == 'heavy':
                ## heavy atoms
                if atom_name[0] == 'H':
                    continue
            else:
                print method
                stop

            altloc = d_mmCIF['_atom_site.label_alt_id'][i]

##            ## choose only the first altloc
##            if altloc not in ['.','A','1',]:
##                continue

##            ## first altloc has less than 50% occupancy
##            if altloc != '.' and occupancy < 0.5:
##                print 'altloc', altloc, 'oocupancy', occupancy, pdb, res_no, res_name, atom_name

##            ## alpha carbon atom has two altlocs
##            if altloc != '.' and atom_name == 'CA':
##                print 'altloc', altloc, pdb, res_no, res_name, atom_name

            ## altloc has zero occupancy
            if pdb[:4] not in [
                '2HUB', ## HEWL?
                '3CDO','3F9L','3FAD','3FI5',
                '3HTD',
                '2RB1',
                '1XEP',
                '3HH3',
                ] and occupancy == 0 and altloc != '.':
                print pdb, res_no, res_name, atom_name, altloc
                stop

            ## exclude zero occupancy atoms
            if bool_exclusion_zero_occupancy == True and occupancy == 0 and altloc == '.':
                continue
            ## exclude altlocs
            if bool_exclusion_altlocs == True and altloc != '.':
                continue
##            ## exclude high temp factor atoms
##            if bool_exclusion_high_temp_factors == True and tempFactor > 50:
##                continue
            ## exclude symmetry atoms
            if bool_exclusion_symmetry == True and res_name in d_symmetry.keys():
                if atom_name in d_symmetry[res_name]:
                    continue

            ##
            ## append coordinate
            ##
            if not res_no in d_coords.keys():
                d_coords[res_no] = {'res_name':res_name,'atoms':{},}
            if not altloc in d_coords[res_no]['atoms'].keys():
                d_coords[res_no]['atoms'][altloc] = {}
            d_coords[res_no]['atoms'][altloc][atom_name] = {'coord':coord,'tempFactor':tempFactor,'occupancy':occupancy,}

        ## select altloc
        l_res_nos = []
        for res_no in d_coords.keys():

            if res_no in l_res_nos:
                continue

            l_altlocs = d_coords[res_no]['atoms'].keys()

            if l_altlocs == ['.']:

                d_coords[res_no]['atoms'] = d_coords[res_no]['atoms']['.']

            else:

                if protein == 'T4L' and d_mmCIF['_symmetry.space_group_name_H-M'] != ['P 32 2 1']:
                    altloc = 'A'
                    l_res_nos = [res_no]

                elif protein == 'HEWL':

                    l_res_nos = [res_no]

                    if pdb == '252L_A':
                        altloc = 'B'
                    elif pdb == '3HH6_A' and res_no in range(108,115):
                        altloc = 'B'
                    else:
                        altloc = 'A'

                else:

                    ## loop over neigboring residues (e.g. 1lw9 54-56)
                    l_res_nos = []
                    d_tempFactors = {}
                    d_occupancies = {}
                    min_tempFactor = ['N/A','N/A',]
                    max_occupancy = [0,'N/A',]
                    for res_no in range(res_no,max(d_coords.keys())+1):
                        ## break if next residue has no altloc
                        if d_coords[res_no]['atoms'].keys() == ['.']:
                            break
                        else:
                            res_no_altloc = res_no
                        l_altlocs = d_coords[res_no_altloc]['atoms'].keys()
                        l_res_nos += [res_no_altloc]
                        if len(l_altlocs) < 3:
                            if l_altlocs != ['A','B',]:
                                if pdb == '1PQM_A' and res_no_altloc == 55:
                                    pass
                                elif pdb == '3C7Y_A' and res_no_altloc in [72,76,]:
                                    pass
                                elif pdb == '3CDQ_A' and res_no_altloc in [16,76,]:
                                    pass
                                elif pdb == '3C7W_A' and res_no_altloc in [76,]:
                                    pass
                                elif pdb == '1SX2_A' and res_no_altloc in [123,]:
                                    pass
                                else:
                                    print l_altlocs
                                    print pdb, res_no_altloc
##                                    stop
                        for altloc in l_altlocs:
                            if altloc == '.':
                                continue
                            d_tempFactors[altloc] = []
                            d_occupancies[altloc] = []
                            for atom_name in d_coords[res_no_altloc]['atoms'][altloc].keys():
                                tempFactor = d_coords[res_no_altloc]['atoms'][altloc][atom_name]['tempFactor']
                                occupancy = d_coords[res_no_altloc]['atoms'][altloc][atom_name]['occupancy']
                                d_occupancies[altloc] += [occupancy]
                                d_tempFactors[altloc] += [tempFactor]

                            l_tempFactors = d_tempFactors[altloc]
                            tempFactor_average = sum(l_tempFactors)/len(l_tempFactors)
                            if tempFactor_average < min_tempFactor[0]:
                                min_tempFactor = [tempFactor_average,altloc,]

                            l_occupancies = d_occupancies[altloc]
                            occupancy_average = sum(l_occupancies)/len(l_occupancies)
                            if occupancy != 0.5 and occupancy_average > max_occupancy[0]:
                                max_occupancy = [occupancy_average,altloc,]

                            ## check that all occupancies are the same
                            if len(l_occupancies)*[l_occupancies[0]] != l_occupancies and pdb not in ['1P2L_A','1SWZ_A','1ZWN_A',]:
                                print pdb, res_no_altloc, altloc, l_occupancies
##                                stop_different_occupancies

                    ## 50/50 occupancies
                    if max_occupancy[1] == 'N/A':
                        if pdb == '1PQM_A' and res_no_altloc in [76,106,]:
                            altloc = 'B'
                        elif pdb == '3C7Y_A' and res_no_altloc == 76:
                            altloc = 'B' ## lower bfactor
                        elif pdb == '3DMV_A' and res_no_altloc in [56,108,]:
                            altloc = 'B'
                        elif pdb == '3CDQ_A' and res_no_altloc in [1,16,76,]:
                            altloc = 'B' ## should it be A if it existed? in the case of 16,76.
                        elif pdb == '3F8V_A' and res_no_altloc in [106,]:
                            altloc = 'B' ## lower bfactor
                        elif pdb == '3FA0_A' and res_no_altloc in [44,]:
                            altloc = 'B' ## lower bfactor
                        elif pdb == '3FAD_A' and res_no_altloc in [109,]:
                            altloc = 'B' ## lower bfactor
                        elif pdb == '3CDQ_A' and res_no_altloc in [106,]:
                            altloc = 'B'
                        elif pdb == '3CDR_A' and res_no_altloc in [16,76,]:
                            altloc = 'B' ## should it be A if it existed?
                        elif pdb == '3CDV_A' and res_no_altloc in [16,76,]:
                            altloc = 'B' ## should it be A if it existed?
                        elif pdb == '3C7W_A' and res_no_altloc == 76:
                            altloc = 'B' ## should it be A if it existed?
                        elif pdb == '3C7Z_A' and res_no_altloc == 76:
                            altloc = 'B' ## should it be A if it existed?
                        elif pdb == '3C81_A' and res_no_altloc == 76:
                            altloc = 'B' ## should it be A if it existed?
                        elif pdb == '3C8Q_A' and res_no_altloc == 76:
                            altloc = 'B' ## should it be A if it existed?
                        elif pdb == '3C8S_A' and res_no_altloc in [16,76,]:
                            altloc = 'B' ## should it be A if it existed?
                        else:
                            altloc = 'A'

                    ## not 50/50 occupancies
                    elif max_occupancy[1] != 'N/A':
                        altloc = max_occupancy[1]
                        if altloc == 'B':
                            if pdb == 'xxx1LW9_A' and res_no_altloc == 137:
                                altloc = max_occupancy[1]

                            ## error in pdb?
                            elif pdb == '3CDQ_A' and res_no_altloc in [6,124,]:
                                altloc = 'A' ## ???
                            elif pdb == '3DMV_A' and res_no_altloc in [16,44,68,106,125,]:
                                altloc = 'A' ## ???
                            elif pdb == '3F9L_A' and res_no_altloc in [44,108,]:
                                altloc = 'A' ## ???
                            elif pdb == '3FAD_A' and res_no_altloc in [68,]:
                                altloc = 'A' ## ???

                            elif pdb == '1PQM_A' and res_no_altloc == 44:
                                altloc = 'B'
                            elif pdb == '3CDO_A' and res_no_altloc in [1,80,83,119,135,]:
                                altloc = 'B'
                            elif pdb == '3F9L_A' and res_no_altloc == 106:
                                pass  ## ???
                            elif pdb == '3FAD_A' and res_no_altloc == 57:
                                pass ### ???
                            elif pdb == '1LPY_A' and res_no_altloc in [91,99,]:
                                pass ## mutated residue

                            else:
                                print d_tempFactors
                                print d_occupancies
                                print 'A', sum(d_tempFactors['A'])/len(d_tempFactors['A'])
                                print 'B', sum(d_tempFactors['B'])/len(d_tempFactors['B'])
                                print max_occupancy
                                print pdb, res_no_altloc
                                print 'diff occ bfac'
                                pass
##                                stop_different_occupancies

                        if pdb == '1PQM_A' and res_no_altloc == 55:
                            altloc = 'B' ## should be A if it existed

                if res_no == 109 and pdb == '3HH6_A':
                    print altloc
                    stop

                ##
                ## replace 2 or more altlocs with 1 altloc
                ##
                for res_no_altloc in l_res_nos:
                    d = {}
                    ## append non-altloc
                    if '.' in d_coords[res_no_altloc]['atoms'].keys():
                        for atom_name in d_coords[res_no_altloc]['atoms']['.'].keys():
                            d[atom_name] = d_coords[res_no_altloc]['atoms']['.'][atom_name]
                    ## append altloc
                    try:
                        for atom_name in d_coords[res_no_altloc]['atoms'][altloc].keys():
                            d[atom_name] = d_coords[res_no_altloc]['atoms'][altloc][atom_name]
                        d_coords[res_no_altloc]['atoms'] = d
                    except:
                        print pdb, res_no, res_no_altloc, altloc
                        print max_occupancy
                        print d_coords[res_no_altloc]['atoms'][altloc].keys()
                        stop

        d_coordinates[pdb] = d_coords

    return d_coordinates


def get_coordinates(
    method, d_coordinates, res_no,
    pdb1, pdb2,
    ):

    l_coords1 = []
    l_coords2 = []
    l_bfacs1 = []
    l_bfacs2 = []
    try:
        set_res_names1 = set([d_coordinates[pdb1][res_no]['res_name']])
        set_res_names2 = set([d_coordinates[pdb2][res_no]['res_name']])
    except:
        print d_coordinates.keys()
        print pdb1, pdb2
        print res_no
        set_res_names1 = set([d_coordinates[pdb1][res_no]['res_name']])
        set_res_names2 = set([d_coordinates[pdb2][res_no]['res_name']])
        stop
    set_res_names = set_res_names1 & set_res_names2

    ##
    ## no mutation
    ##
    if len(set_res_names) == 1:
        res_name = ''.join(set_res_names)
        if list(set_res_names)[0] != ''.join(set_res_names):
            print list(set_res_names)[0]
            print ''.join(set_res_names)
            stop
        for res_name in set_res_names:
            set_atom_names1 = set(d_coordinates[pdb1][res_no]['atoms'].keys())
            set_atom_names2 = set(d_coordinates[pdb2][res_no]['atoms'].keys())
            set_atom_names = set_atom_names1 & set_atom_names2
            if method == 'alpha':
                set_atom_names = set(['CA'])
            for atom_name in set_atom_names:

                ## exclude symmetry atoms
                if method == 'heavy' and res_name in d_symmetry.keys():
                    if atom_name in d_symmetry[res_name]:
                        continue
                        
                if atom_name[0] == 'H':
                    print res_name, atom_name
                    stop_hydrogen
                coord1 = d_coordinates[pdb1][res_no]['atoms'][atom_name]['coord']
                try:
                    coord2 = d_coordinates[pdb2][res_no]['atoms'][atom_name]['coord']
                except:
                    print pdb1, pdb2
                    print res_no, atom_name
                    print res_no, d_coordinates[pdb2].keys()
                    print atom_name, d_coordinates[pdb2][res_no]['atoms'].keys()
                    print d_coordinates[pdb2][res_no]['atoms'][atom_name].keys()
                    print d_coordinates[pdb2][res_no]['atoms'][atom_name]['coord']
                    stop
                l_bfacs1 += [d_coordinates[pdb1][res_no]['atoms'][atom_name]['tempFactor']]
                l_bfacs2 += [d_coordinates[pdb2][res_no]['atoms'][atom_name]['tempFactor']]
                l_coords1 += [coord1]
                l_coords2 += [coord2]
    ##
    ## mutation
    ##
    elif len(set_res_names) == 0:
        if len(set_res_names1) == 1 and len(set_res_names2) == 1:
            res_name1 = list(set_res_names1)[0]
            res_name2 = list(set_res_names2)[0]
            set_atom_names1 = set(d_coordinates[pdb1][res_no]['atoms'].keys())
            set_atom_names2 = set(d_coordinates[pdb2][res_no]['atoms'].keys())

            ## use backbone atoms only if mutation
            if method == 'heavy':
                set_atom_names = set(['N','CA','CB','C','O',])
            elif method == 'alpha':
                set_atom_names = set(['CA'])

            set_atom_names = set_atom_names & set_atom_names1 & set_atom_names2

            for atom_name in set_atom_names:
                if atom_name[0] == 'H':
                    stop
                coord1 = d_coordinates[pdb1][res_no]['atoms'][atom_name]['coord']
                coord2 = d_coordinates[pdb2][res_no]['atoms'][atom_name]['coord']
                l_bfacs1 += [d_coordinates[pdb1][res_no]['atoms'][atom_name]['tempFactor']]
                l_bfacs2 += [d_coordinates[pdb2][res_no]['atoms'][atom_name]['tempFactor']]
                l_coords1 += [coord1]
                l_coords2 += [coord2]
        else:
            print set_res_names1
            print set_res_names2
            stop
    ##
    ## altloc
    ##
    else:
        print set_res_names
        stop

    return l_coords1,l_coords2, l_bfacs1, l_bfacs2


def calculate_rmsds(
    l_pdbs,d_coordinates,method,protein,
    suffix_exclusion,
    d_properties=None,do_subset=False,
    bool_exclusion_zero_occupancy = True,
    bool_exclusion_symmetry = False,
    bool_exclusion_altlocs = False,
    bool_exclusion_high_temp_factors = False,
    do_rmsd_per_residue=False,
    do_bfactors=False,
    ):

    print 'calculate rmsds'

    ##
    ## initiate dictionaries and lists
    ##
    d_rmsds_overall = {}
    for pdb in l_pdbs:
        d_rmsds_overall[pdb] = {}
        
    d_rmsds_per_residue = {}
    for i in range(1,129+1):
        d_rmsds_per_residue[i] = []
        
    l_gnuplot_bfactors = []

    ##
    ## loop over pdbs
    ##
    for i1 in range(len(l_pdbs)-1):
        pdb1 = l_pdbs[i1]

        print 'calculate rmsds', pdb1

        for i2 in range(i1+1,len(l_pdbs)):
            pdb2 = l_pdbs[i2]

##            print 'calc rmsd', pdb1, pdb2

            ##
            ## overall alignment and calculation of overall rmsd and transformation matrix
            ##
            l_coordinates1 = []
            l_coordinates2 = []
            l_bfactors1 = []
            l_bfactors2 = []
            d_coords_per_residue1 = {}
            d_coords_per_residue2 = {}

            set_res_nos1 = set(d_coordinates[pdb1].keys())
            set_res_nos2 = set(d_coordinates[pdb2].keys())
            set_res_nos = set_res_nos1 & set_res_nos2

            for res_no in set_res_nos:

##                if res_no in [1]+range(45,51)+range(67,74)+range(100,104)+range(128,130):
##                    continue

                l_coords1,l_coords2,l_bfacs1,l_bfacs2 = get_coordinates(
                    method, d_coordinates, res_no,
                    pdb1, pdb2,
                    )

                l_coordinates1 += l_coords1
                l_coordinates2 += l_coords2
                l_bfactors1 += l_bfacs1
                l_bfactors2 += l_bfacs2

                if do_rmsd_per_residue == True:
                    d_coords_per_residue1[res_no] = list(l_coords1)
                    d_coords_per_residue2[res_no] = list(l_coords2)

            rmsd_overall = instance_geometry.superpose(l_coordinates1,l_coordinates2,)

            if rmsd_overall == 0:
                stop
            if pdb1 == pdb2:
                stop

            fd = open('rmsds_%s_%s_%s.txt' %(protein,method,suffix_exclusion,),'a')
            fd.write('%s %s %s\n' %(pdb1,pdb2,rmsd_overall,))
            fd.close()

            tv1 = instance_geometry.fitcenter
            rm = instance_geometry.rotation
            tv2 = instance_geometry.refcenter
            rm_inv = numpy.linalg.inv(rm)

            d_rmsds_overall[pdb1][pdb2] = {
                'rmsd':rmsd_overall,'tv1':tv1,'tv2':tv2,'rm':rm,
                }
            d_rmsds_overall[pdb2][pdb1] = {
                'rmsd':rmsd_overall,'tv1':tv2,'tv2':tv1,'rm':rm_inv,
                }


            ##
            ## diff v bfactor, calculate average bfactor and append to list
            ##
            if do_bfactors == True:

                for i in range(len(l_coordinates1)):

                    c1 = l_coordinates1[i]
                    c2 = l_coordinates2[i]
                    c2 = numpy.dot(c2-tv1,rm)+tv2
                    dist = math.sqrt(sum((c2-c1)**2))

                    b1 = l_bfactors1[i]
                    b2 = l_bfactors2[i]
                    b_average = .5*(b1+b2)

                    l_gnuplot_bfactors += ['%s %s\n' %(dist,b_average,)]


            if do_rmsd_per_residue == True:
                for res_no in range(1,129+1):
                    if not res_no in d_coords_per_residue1.keys():
                        continue
                    if not res_no in d_coords_per_residue2.keys():
                        continue
                    l_coords1 = d_coords_per_residue1[res_no]
                    l_coords2 = d_coords_per_residue2[res_no]

                    ## no atom names in common (e.g. SNN101 in 1AT5)
                    if len(l_coords1) == 0 and len(l_coords2) == 0:
                        continue

                    for i in range(len(l_coords2)):
                        coord = l_coords2[i]
                        coord = numpy.dot(coord-tv1,rm)+tv2
                        l_coords2[i] = coord
                    try:
                        rmsd = calc_rmsd_of_pre_aligned_coordinates(l_coords1,l_coords2,)
                    except:
                        print pdb1, pdb2, res_no, l_coords1, l_coords2
                        stop
                    d_rmsds_per_residue[res_no] += [rmsd]

        if len(l_coordinates1) < 100:
            print '*****', pdb1, len(l_coordinates1)


    ##                  
    ## plot diff v bfactor per atom
    ##
    if do_bfactors == True:
        fd = open('bfactor.gnuplot_data','w')
        fd.writelines(l_gnuplot_bfactors)
        fd.close()


    return d_rmsds_overall, d_rmsds_per_residue


def calc_rmsd_of_pre_aligned_coordinates(l_coords1,l_coords2,):

    l_diff_sq = []
    for i_coord in range(len(l_coords1)):
        diff = l_coords1[i_coord]-l_coords2[i_coord]
        diff_sq = sum(diff**2)
        l_diff_sq += [diff_sq]

    mean = sum(l_diff_sq)/len(l_diff_sq)
    rmsd = math.sqrt(mean)
    
    return rmsd


if __name__ == '__main__':
    main()

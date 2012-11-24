## HEWL contour plot of RMSDs

import numpy, os, math
import core
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import gnuplot, statistics
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import parse_mmCIF
sys.path.append('/home/people/tc/svn/Protool/')
import geometry
instance_geometry = geometry.geometry()
import parse_MD_NMR_coordinates

d_proteins = {
    5:{
        'protein':'HEWL',
        'db_code':['LYSC_CHICK','LYSC_CHIC','LYC_CHICK',],
        'ref_seq':
        ['LYS', 'VAL', 'PHE', 'GLY', 'ARG', 'CYS', 'GLU', 'LEU', 'ALA', 'ALA', 'ALA', 'MET', 'LYS', 'ARG', 'HIS', 'GLY', 'LEU', 'ASP', 'ASN', 'TYR', 'ARG', 'GLY', 'TYR', 'SER', 'LEU', 'GLY', 'ASN', 'TRP', 'VAL', 'CYS', 'ALA', 'ALA', 'LYS', 'PHE', 'GLU', 'SER', 'ASN', 'PHE', 'ASN', 'THR', 'GLN', 'ALA', 'THR', 'ASN', 'ARG', 'ASN', 'THR', 'ASP', 'GLY', 'SER', 'THR', 'ASP', 'TYR', 'GLY', 'ILE', 'LEU', 'GLN', 'ILE', 'ASN', 'SER', 'ARG', 'TRP', 'TRP', 'CYS', 'ASN', 'ASP', 'GLY', 'ARG', 'THR', 'PRO', 'GLY', 'SER', 'ARG', 'ASN', 'LEU', 'CYS', 'ASN', 'ILE', 'PRO', 'CYS', 'SER', 'ALA', 'LEU', 'LEU', 'SER', 'SER', 'ASP', 'ILE', 'THR', 'ALA', 'SER', 'VAL', 'ASN', 'CYS', 'ALA', 'LYS', 'LYS', 'ILE', 'VAL', 'SER', 'ASP', 'GLY', 'ASN', 'GLY', 'MET', 'ASN', 'ALA', 'TRP', 'VAL', 'ALA', 'TRP', 'ARG', 'ASN', 'ARG', 'CYS', 'LYS', 'GLY', 'THR', 'ASP', 'VAL', 'GLN', 'ALA', 'TRP', 'ILE', 'ARG', 'GLY', 'CYS', 'ARG', 'LEU'],
        },
    1:{
        'protein':'T4L',
        'db_code':['LYS_BPT4','LYCV_BPT4',],
        'ref_seq':
        ['MET', 'ASN', 'ILE', 'PHE', 'GLU', 'MET', 'LEU', 'ARG', 'ILE', 'ASP', 'GLU', 'GLY', 'LEU', 'ARG', 'LEU', 'LYS', 'ILE', 'TYR', 'LYS', 'ASP', 'THR', 'GLU', 'GLY', 'TYR', 'TYR', 'THR', 'ILE', 'GLY', 'ILE', 'GLY', 'HIS', 'LEU', 'LEU', 'THR', 'LYS', 'SER', 'PRO', 'SER', 'LEU', 'ASN', 'ALA', 'ALA', 'LYS', 'SER', 'GLU', 'LEU', 'ASP', 'LYS', 'ALA', 'ILE', 'GLY', 'ARG', 'ASN', 'CYS', 'ASN', 'GLY', 'VAL', 'ILE', 'THR', 'LYS', 'ASP', 'GLU', 'ALA', 'GLU', 'LYS', 'LEU', 'PHE', 'ASN', 'GLN', 'ASP', 'VAL', 'ASP', 'ALA', 'ALA', 'VAL', 'ARG', 'GLY', 'ILE', 'LEU', 'ARG', 'ASN', 'ALA', 'LYS', 'LEU', 'LYS', 'PRO', 'VAL', 'TYR', 'ASP', 'SER', 'LEU', 'ASP', 'ALA', 'VAL', 'ARG', 'ARG', 'CYS', 'ALA', 'LEU', 'ILE', 'ASN', 'MET', 'VAL', 'PHE', 'GLN', 'MET', 'GLY', 'GLU', 'THR', 'GLY', 'VAL', 'ALA', 'GLY', 'PHE', 'THR', 'ASN', 'SER', 'LEU', 'ARG', 'MET', 'LEU', 'GLN', 'GLN', 'LYS', 'ARG', 'TRP', 'ASP', 'GLU', 'ALA', 'ALA', 'VAL', 'ASN', 'LEU', 'ALA', 'LYS', 'SER', 'ARG', 'TRP', 'TYR', 'ASN', 'GLN', 'THR', 'PRO', 'ASN', 'ARG', 'ALA', 'LYS', 'ARG', 'VAL', 'ILE', 'THR', 'THR', 'PHE', 'ARG', 'THR', 'GLY', 'THR', 'TRP', 'ASP', 'ALA', 'TYR', 'LYS', 'ASN', 'LEU'],
        },
    }

d_321 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
    'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
    'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
    'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
    }

d_123 = {}
for k,v in d_321.items():
    d_123[v] = k

## selection/restriction
i_cluster = int(sys.argv[sys.argv.index('-cluster')+1])
#### accept ligands?
bool_ligands = True
## accept mutants?
n_mutations_max = 100
## accept multiple sequence different chains?
bool_multiple_entities = False
## spacegroups identical when plotting?
bool_spacegroups_identical = True
## minimum resolution
resolution_min = 3.9
resolution_min = 2.00
## exclusion for calculating RMSD?
bool_exclusion_zero_occupancy = True ## True
bool_exclusion_symmetry = True ## True (relevant if heavy atom RMSD or Val chi1 diff)
bool_exclusion_high_temp_factors = False ## True (otherwise unobserved residues used)
max_bfactor = 30
bool_exclusion_altlocs = False ## False! do both altlocs instead

##bool_exclusion_altlocs = False ## False
##bool_exclusion_high_temp_factors = False ## False

suffix_exclusion = 'occ%s_alt%s_temp%s_symm%s' %(
    bool_exclusion_zero_occupancy,
    bool_exclusion_altlocs,
    bool_exclusion_high_temp_factors,
    bool_exclusion_symmetry,
    )

## methods
method = sys.argv[sys.argv.index('-method')+1]

protein = d_proteins[i_cluster]['protein']
ref_seq = d_proteins[i_cluster]['ref_seq']
l_db_codes = d_proteins[i_cluster]['db_code']

if protein == 'HEWL':
    d_properties = {
        'flexibility':{
            'fixed':range(2,45)+range(51,67)+range(74,100)+range(104,128),
            'flexible':[1]+range(45,51)+range(67,74)+range(100,104)+range(125,130),
            }
        }
elif protein == 'T4L':
    d_properties = {
        'flexibility':{
            'fixed':range(2,106)+range(115,163),
            ## 106-114 contains Gly107, Gly110, Gly113
            'flexible':[1]+range(106,115)+range(162,165),
            }
        }

def main():

    l_pdbs = core.parse_pdbs_from_cluster(i_cluster,)
    if protein == 'HEWL':
        for pdb in [
            ## glycosylated
            '1H6M_A','1UC0_A','2B5Z_A',
            ## covalently crosslinked
            '2HTX_A','2HU1_A',
##            '2LYO_A','3LYO_A', ## CCN
            ## radiation decay
            '3LYT_A','3LYT_B','4LYT_A','4LYT_B', ## monoclinic
            '5LYT_A', ## 100K, tetragonal
            '6LYT_A', ## 298K, tetragonal
            '1QIO_A',
            ## derived from 6Angstrom data
            '1LZT_A','1AKI_A',
            ## no structure factors (real reason, 1974, 2.0Angstrom, VISIBLE problems with geometry and stereochemistry, Ramachandran46,49,50,71)
##            '1LYZ_A','2LYZ_A','3LYZ_A','4LYZ_A','5LYZ_A','6LYZ_A',##'1HEW_A',
##            '1LSA_A', ## different from 1iee (due to C-terminal?)
            ## MODRES
            '132L_A', # methylated lysine (DM0, MLY)
            '1AT5_A', # succinimide (SNN)
            '1AT6_A', # isoaspartate (ASP)
            '1RCM_A','1RCM_B', # carboxymethylated cysteine (CCS)
            ]:
            l_pdbs.remove(pdb)
    elif protein == 'T4L':
        for pdb in [
##            ## T4L intermolecular disulfide bond
##            '2HUK_A','2HUL_A','2HUM_A','2HUM_B',
##            '3HWL_A',
            '3GUK_A','3GUK_B', ## average b-factor high
            ## spin-labeled
            '1ZYT_A','2CUU_A','3G3V_A','3G3W_A','3G3X_A',
            '2Q9D_A','2Q9E_A',
            '2IGC_A','2NTG_A','2NTH_A','2OU8_A','2OU9_A',
            ]:
            l_pdbs.remove(pdb)

##    l_pdbs = [ ## tmp!!!
##        ## wt
##        '2LZM_A','3LZM_A','3FA0_A','3FAD_A',
##        '4LZM_A','5LZM_A','6LZM_A','7LZM_A','1LYD_A',
##        '3F8V_A', ## R96H
##        ## wt*
##        '1LW9_A','219L_A','1L63_A',
####        ## random structures
####        '107L_A','3DMV_A','1L90_A','2Q9D_A','2HUK_A','1C6T_A','1PQM_A','3CDO_A','3CDQ_A','3F9L_A','1P2L_A','1LPY_A','1SWZ_A','1XEP_A','2OTY_X','3C7W_A',
##        ## medium chi1 changes
##        '1G1V_A','3C7Y_A',
##        ## large chi1 changes
##        '200L_A','1QUG_A','1LLH_A','1TLA_A',
##        ] ## tmp!!!
##    l_pdbs = ['107L_A', '108L_A', '109L_A', '110L_A', '111L_A', '112L_A', '113L_A', '114L_A', '115L_A', '118L_A', '119L_A', '120L_A', '122L_A', '123L_A', '125L_A', '126L_A', '127L_A', '128L_A', '129L_A', '130L_A', '131L_A', '137L_A', '137L_B', '138L_A', '139L_A', '141L_A', '142L_A', '143L_A', '145L_A', '146L_A', '147L_A', '152L_A', '155L_A', '156L_A', '157L_A', '158L_A', '159L_A', '160L_A', '161L_A', '162L_A', '163L_A', '164L_A', '165L_A', '166L_A', '172L_A', '173L_A', '180L_A', '180L_B', '181L_A', '182L_A', '183L_A', '184L_A', '185L_A', '186L_A', '187L_A', '188L_A', '190L_A', '191L_A', '192L_A', '195L_A', '198L_A', '199L_A', '1B6I_A', '1C60_A', '1C61_A', '1C63_A', '1C64_A', '1C65_A', '1C69_A', '1C6C_A', '1C6D_A', '1C6E_A', '1C6F_A', '1C6G_A', '1C6H_A', '1C6I_A', '1C6J_A', '1C6K_A', '1C6L_A', '1C6P_A', '1C6Q_A', '1C6T_A', '1CU2_A', '1CUP_A', '1CV3_A', '1CV4_A', '1CV5_A', '1CV6_A', '1CVK_A', '1CX7_A', '1D2W_A', '1D3J_A', '1D3N_A', '1D9W_A', '1DYA_A', '1DYB_A', '1DYE_A', '1DYF_A', '1EPY_A', '1G06_A', '1G07_A', '1G0G_A', '1G0J_A', '1G0K_A', '1G0L_A', '1G0M_A', '1G0P_A', '1G0Q_A', '1G1V_A', '1G1W_A', '1I6S_A', '1KNI_A', '1KW5_A', '1KW7_A', '1KY0_A', '1L00_A', '1L01_A', '1L02_A', '1L03_A', '1L04_A', '1L05_A', '1L06_A', '1L07_A', '1L08_A', '1L09_A', '1L0J_A', '1L10_A', '1L11_A', '1L12_A', '1L13_A', '1L14_A', '1L15_A', '1L16_A', '1L17_A', '1L18_A', '1L19_A', '1L20_A', '1L21_A', '1L22_A', '1L23_A', '1L24_A', '1L25_A', '1L26_A', '1L27_A', '1L28_A', '1L29_A', '1L30_A', '1L31_A', '1L32_A', '1L33_A', '1L34_A', '1L35_A', '1L36_A', '1L37_A', '1L38_A', '1L39_A', '1L40_A', '1L41_A', '1L42_A', '1L43_A', '1L44_A', '1L45_A', '1L46_A', '1L47_A', '1L48_A', '1L49_A', '1L50_A', '1L51_A', '1L52_A', '1L53_A', '1L54_A', '1L55_A', '1L56_A', '1L57_A', '1L58_A', '1L59_A', '1L60_A', '1L61_A', '1L62_A', '1L63_A', '1L64_A', '1L65_A', '1L66_A', '1L67_A', '1L68_A', '1L69_A', '1L70_A', '1L71_A', '1L72_A', '1L73_A', '1L74_A', '1L75_A', '1L76_A', '1L79_A', '1L80_A', '1L81_A', '1L83_A', '1L84_A', '1L85_A', '1L86_A', '1L87_A', '1L88_A', '1L89_A', '1L90_A', '1L91_A', '1L92_A', '1L93_A', '1L94_A', '1L95_A', '1L96_A', '1L97_A', '1L97_B', '1L98_A', '1L99_A', '1LGU_A', '1LGW_A', '1LGX_A', '1LI2_A', '1LI3_A', '1LI6_A', '1LLH_A', '1LPY_A', '1LW9_A', '1LWG_A', '1LYD_A', '1LYE_A', '1LYF_A', '1LYG_A', '1LYH_A', '1LYI_A', '1LYJ_A', '1NHB_A', '1OV7_A', '1OVH_A', '1OVJ_A', '1OWY_A', '1OWZ_A', '1P2L_A', '1P2R_A', '1P36_A', '1P37_A', '1P3N_A', '1P46_A', '1P64_A', '1P6Y_A', '1P7S_A', '1PQD_A', '1PQI_A', '1PQJ_A', '1PQK_A', '1PQK_B', '1PQK_C', '1PQM_A', '1PQO_A', '1QS9_A', '1QSB_A', '1QSQ_A', '1QT3_A', '1QT5_A', '1QT6_A', '1QT7_A', '1QT8_A', '1QTB_A', '1QTH_A', '1QTH_B', '1QTZ_A', '1QUD_A', '1QUG_A', '1QUH_A', '1QUO_A', '1SWY_A', '1SWZ_A', '1SX2_A', '1SX7_A', '1T8G_A', '1TLA_A', '1XEP_A', '1ZUR_A', '1ZWN_A', '1ZYT_A', '200L_A', '206L_A', '217L_A', '219L_A', '220L_A', '221L_A', '222L_A', '223L_A', '224L_A', '225L_A', '226L_A', '227L_A', '228L_A', '229L_A', '230L_A', '232L_A', '233L_A', '234L_A', '235L_A', '236L_A', '237L_A', '238L_A', '239L_A', '240L_A', '241L_A', '242L_A', '243L_A', '244L_A', '245L_A', '246L_A', '247L_A', '248L_A', '249L_A', '250L_A', '253L_A', '254L_A', '255L_A', '256L_A', '257L_A', '258L_A', '259L_A', '260L_A', '2A4T_A', '2CUU_A', '2HUK_A', '2HUL_A', '2IGC_A', '2L78_A', '2LZM_A', '2NTG_A', '2NTH_A', '2OTY_X', '2OU0_X', '2OU8_A', '2OU9_A', '2Q9D_A', '2RAY_X', '2RAZ_X', '2RB0_X', '2RB1_X', '2RB2_X', '2RBN_A', '2RBO_A', '2RBP_A', '2RBQ_A', '2RBR_A', '2RBS_A', '3C7W_A', '3C7Y_A', '3C7Z_A', '3C80_A', '3C81_A', '3C82_A', '3C83_A', '3C8Q_A', '3C8R_A', '3C8S_A', '3CDO_A', '3CDO_B', '3CDO_C', '3CDO_D', '3CDQ_A', '3CDR_A', '3CDT_A', '3CDV_A', '3DKE_X', '3DMV_A', '3DMX_A', '3DMZ_A', '3DN0_A', '3DN1_A', '3DN2_A', '3DN3_A', '3DN4_A', '3DN6_A', '3DN8_A', '3DNA_A', '3F8V_A', '3F9L_A', '3FA0_A', '3FAD_A', '3FI5_A', '3FI5_B', '3FI5_C', '3FI5_D', '3G3X_A', '3GUI_A', '3GUJ_A', '3GUN_A', '3GUN_B', '3GUP_A', '3GUP_B', '3HH3_A', '3HH4_A', '3HH5_A', '3HH6_A', '3HT6_A', '3HT7_A', '3HT8_A', '3HTB_A', '3HTD_A', '3HTF_A', '3HTG_A', '3HU8_A', '3HU9_A', '3HUA_A', '3HUK_A', '3HUQ_A', '3HWL_A', '3L64_A', '3LZM_A', '4LZM_A', '5LZM_A', '6LZM_A', '7LZM_A',]
##    l_pdbs = list(set(['2VB1_A','2LZT_A','193L_A','1BWJ_A',]+l_pdbs[:10]))

    d_mmCIF, l_wts, d_mutants, l_wts_cysfree = core.parse_cifs(
        l_pdbs,
        ref_seq, l_db_codes,
        n_mutations_max,
        resolution_min,
        bool_multiple_entities = bool_multiple_entities,
        )

    l_mutants = d_mutants.keys()
    l_pdbs = list(l_wts+l_mutants)
    l_pdbs.sort()

    d_pdbs = {}
    for pdb in l_pdbs:
        if not pdb[:4] in d_pdbs.keys():
            d_pdbs[pdb[:4]] = []
        d_pdbs[pdb[:4]] += [pdb]

##    for pdb in l_pdbs:
##        if '_struct_conn.conn_type_id' not in d_mmCIF[pdb[:4]].keys():
##            continue
##        if not 'disulf' in d_mmCIF[pdb[:4]]['_struct_conn.conn_type_id']:
##            continue
##        intermolecular = False
##        for i in range(len(d_mmCIF[pdb[:4]]['_struct_conn.conn_type_id'])):
##            if not (
##                d_mmCIF[pdb[:4]]['_struct_conn.ptnr1_label_asym_id'][i] == d_mmCIF[pdb[:4]]['_struct_conn.ptnr2_label_asym_id'][i]
##                and
##                d_mmCIF[pdb[:4]]['_struct_conn.ptnr1_symmetry'][i] == '1_555'
##                and
##                d_mmCIF[pdb[:4]]['_struct_conn.ptnr2_symmetry'][i] == '1_555'
##                ):
##                intermolecular = True
##        if intermolecular == True:
##            print pdb
##    stop

##    for pdb in l_pdbs:
##        if '_struct_conn.conn_type_id' not in d_mmCIF[pdb[:4]].keys():
##            continue
##        if not 'covale' in d_mmCIF[pdb[:4]]['_struct_conn.conn_type_id']:
##            continue
##        proteinligand = False
##        for i in range(len(d_mmCIF[pdb[:4]]['_struct_conn.conn_type_id'])):
##            if d_mmCIF[pdb[:4]]['_struct_conn.conn_type_id'][i] != 'covale':
##                continue
##            comp_id1 = d_mmCIF[pdb[:4]]['_struct_conn.ptnr1_label_comp_id'][i]
##            comp_id2 = d_mmCIF[pdb[:4]]['_struct_conn.ptnr2_label_comp_id'][i]
##            if not (
##                comp_id1 == 'NAG'
##                and
##                comp_id2 == 'NAG'
##                ):
##                print pdb, comp_id1, comp_id2
##                proteinligand = True
##        if proteinligand == True:
##            print pdb
##    stop

    ##
    ## calculate matthews coefficient
    ##
    d_mmCIF, d_matthews, d_Z = matthews(d_mmCIF,l_pdbs,)

    ##
    ## tabulate properties
    ##
    tabulate(d_mmCIF,ref_seq,l_wts,)

    ## parse citation
    d_citation = parse_citation(l_pdbs,d_mmCIF,)

    ## parse other properties
    (
        d_ph, d_ph_reverse,
        d_temperature, d_temperature_reverse,
        d_startingmodel, d_startingmodel_reverse,
        d_spacegroups, d_spacegroups_reverse,
        d_resolutions,d_resolutions_reverse,
        d_author,d_author_reverse,
        d_countentities, d_countentities_reverse,
        ) = parse_properties(d_mmCIF,l_pdbs,)

##    ## worst resolution?
##    print max(d_resolutions.keys())
##    print d_resolutions[max(d_resolutions.keys())]
##    stop

##    ## trace back starting model
##    (
##        d_startingmodel, d_startingmodel_reverse,
##        ) = trace_starting_model(d_startingmodel,d_startingmodel_reverse,)

    ##
    ## write properties to file
    ##
    d = {
        'startmodel':d_startingmodel_reverse,
        'spacegroup':d_spacegroups_reverse,
        'author':d_author_reverse,
        'entities':d_countentities_reverse,
        'matthews':d_matthews,
        'citation':d_citation,
        'Z':d_Z,
        }

    fd = open('d_%s.txt' %(protein),'w')
    fd.write(str(d))
    fd.close()

    ##
    ## parse coordinates
    ##
    d_coordinates = core.parse_coordinates(
        l_pdbs,d_mmCIF,method,
        protein,
        bool_exclusion_zero_occupancy = bool_exclusion_zero_occupancy,
        bool_exclusion_symmetry = bool_exclusion_symmetry,
        bool_exclusion_altlocs = bool_exclusion_altlocs,
        bool_exclusion_high_temp_factors = bool_exclusion_high_temp_factors,
        )

    ##
    ## parse NMR (1e8l) and MDs (2vb1 Glu35NEU/Asp52CHA)
    ##
    if protein == 'HEWL':
        for t in range(999,10000,1000,):
            s_MD = 'MD%02i_A' %((t+1)/1000.)
            l_pdbs += [s_MD]
            d_mmCIF[s_MD[:4]] = {
                '_symmetry.space_group_name_H-M':d_mmCIF['2VB1']['_symmetry.space_group_name_H-M'],
                }
            fd = open('/local/tc/MD_2vb1/amber99sb_CYM_35C/NEUCHA/trjconv/2vb1_MD%i.pdb' %(t),'r')
            lines = fd.readlines()
            fd.close()
            d_coordinates[s_MD] = parse_MD_NMR_coordinates.main(lines)
        fd = open('/data/pdb-v3.2/e8/pdb1e8l.ent','r')
        lines = fd.readlines()
        fd.close()
        for i_line in range(len(lines)):
            line = lines[i_line]
            record = line[:6].strip()
            if record == 'MODEL':
                i_line1 = i_line+1
                model = int(line[10:14])
                s_NMR = '1E8L_A%02i' %(model)
                l_pdbs += [s_NMR]
                d_mmCIF[s_NMR[:4]] = {
                    '_symmetry.space_group_name_H-M':None,
                    }
            elif record == 'ENDMDL':
                i_line2 = i_line-1
                d_coordinates[s_NMR] = parse_MD_NMR_coordinates.main(lines[i_line1:i_line2])

    ##
    ## calculate chi1 dihedrals
    ##
    d_chi1 = calc_chi1(l_pdbs,d_coordinates)
##    stop

##    ##
##    ## plot selected phi/psi dihedrals
##    ##    
##    if method == 'heavy' and bool_exclusion_altlocs == False and bool_exclusion_high_temp_factors == False:
##        l_range = d_properties['flexibility']['flexible']
##        plot_dihedral(d_coordinates,d_mmCIF,ref_seq,l_pdbs,l_range,l_mutants,{},)
##    stop

    ##
    ## parse rmsds
    ##
    if os.path.isfile('rmsds_%s_%s_%s.txt' %(protein,method,suffix_exclusion,)):
        print 'read previous rmsd calculations'

        d_rmsds_per_residue = None

        d_rmsds_overall = {}
        for pdb in l_pdbs:
            d_rmsds_overall[pdb] = {}

        fd = open('rmsds_%s_%s_%s.txt' %(protein,method,suffix_exclusion,),'r')
        lines = fd.readlines()
        fd.close()

        print suffix_exclusion
        for line in lines:
            (
                pdb1,pdb2,rmsd,
                rm11,rm12,rm13,
                rm21,rm22,rm23,
                rm31,rm32,rm33,
                tv11,tv12,tv13,
                tv21,tv22,tv23,
                ) = line.strip().split()
            if not pdb1 in l_pdbs:
                continue
            if not pdb2 in l_pdbs:
                continue
            tv1 = numpy.array([float(tv11),float(tv12),float(tv13),])
            tv2 = numpy.array([float(tv21),float(tv22),float(tv23),])
            rm = numpy.array([
                [float(rm11),float(rm12),float(rm13),],
                [float(rm21),float(rm22),float(rm23),],
                [float(rm31),float(rm32),float(rm33),],
                ])
            rm_inv = numpy.linalg.inv(rm)
            d1 = {'rm':rm,'tv1':tv1,'tv2':tv2,'rmsd':float(rmsd),}
##            d2 = {'rm':rm,'tv1':tv1,'tv2':tv2,'rmsd':float(rmsd),}
            d2 = {'rm':rm_inv,'tv1':tv2,'tv2':tv1,'rmsd':float(rmsd),}
            d_rmsds_overall[pdb1][pdb2] = d1
            d_rmsds_overall[pdb2][pdb1] = d2

    ##
    ## calculate rmsds
    ##
    else:
        
        d_rmsds_overall, d_rmsds_per_residue = core.calculate_rmsds(
            l_pdbs,d_coordinates,method,protein,
            suffix_exclusion,
##            do_subset=True, d_properties=d_properties,
            bool_exclusion_zero_occupancy = bool_exclusion_zero_occupancy,
            bool_exclusion_symmetry = bool_exclusion_symmetry,
            bool_exclusion_altlocs = bool_exclusion_altlocs,
            bool_exclusion_high_temp_factors = bool_exclusion_high_temp_factors,
            do_rmsd_per_residue=True,
            )

        l = []
        for i1 in range(len(l_pdbs)-1):
            pdb1 = l_pdbs[i1]
            for i2 in range(i1+1,len(l_pdbs)):
                pdb2 = l_pdbs[i2]
                rm = d_rmsds_overall[pdb1][pdb2]['rm']
                rm11 = rm[0][0]; rm12 = rm[0][1]; rm13 = rm[0][2]
                rm21 = rm[1][0]; rm22 = rm[1][1]; rm23 = rm[1][2]
                rm31 = rm[2][0]; rm32 = rm[2][1]; rm33 = rm[2][2]
                tv1 = d_rmsds_overall[pdb1][pdb2]['tv1']
                tv11 = tv1[0]; tv12 = tv1[1]; tv13 = tv1[2]
                tv2 = d_rmsds_overall[pdb1][pdb2]['tv2']
                tv21 = tv2[0]; tv22 = tv2[1]; tv23 = tv2[2]
                rmsd = d_rmsds_overall[pdb1][pdb2]['rmsd']
                s = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n' %(
                    pdb1,pdb2,rmsd,
                    rm11,rm12,rm13,
                    rm21,rm22,rm23,
                    rm31,rm32,rm33,
                    tv11,tv12,tv13,
                    tv21,tv22,tv23,
                    )
                l += [s]
        
        fd = open('rmsds_%s_%s_%s.txt' %(protein,method,suffix_exclusion,),'w')
        fd.writelines(l)
        fd.close()

##    if protein == 'T4L' and method == 'heavy':
    if 'jens' == 'jens':
        plot_T4_mutants(
            d_chi1,d_rmsds_overall,d_pdbs,d_coordinates,l_wts,l_wts_cysfree,
            )
 
    ##
    ## plots
    ##
    print 'plotting'

    if method == 'heavy':
        plot_rmsd_vs_property(
            l_pdbs,d_properties,d_rmsds_overall,'_exptl_crystal_grow.pH',d_mmCIF,'pH',method,
            )
        plot_rmsd_vs_property(
            l_pdbs,d_properties,d_rmsds_overall,'matthews',d_mmCIF,'MV',method,
            )

    if (protein == 'T4L' or n_mutations_max == 0) and bool_multiple_entities == False:

##        l_pdbs_startingmodel_grouped = group_pdbs(d_startingmodel,l_pdbs,)
##        l_pdbs_spacegroup_grouped = group_pdbs(d_spacegroups,l_pdbs,)
##        l_pdbs_ph_grouped = group_pdbs(d_ph,l_pdbs,)
##        l_pdbs_temperature_grouped = group_pdbs(d_temperature,l_pdbs,)
##    ##    l_pdbs_author_grouped = group_pdbs(d_authors,l_pdbs,)
##        l_pdbs_countentities_grouped = group_pdbs(d_countentities,l_pdbs,)
##
##        ## plot spacegroup if representative space group unknown
##        if spacegroup == None:
##            do_contour_plot(d_rmsds_overall,l_pdbs_spacegroup_grouped,d_spacegroups_reverse,'a_spacegroup',l_wts,)
##        do_contour_plot(d_rmsds_overall,l_pdbs_ph_grouped,d_ph_reverse,'b_ph',l_wts,)
##        do_contour_plot(d_rmsds_overall,l_pdbs_temperature_grouped,d_temperature_reverse,'c_temperature',l_wts,)
##    ##    do_contour_plot(d_rmsds_overall,l_pdbs_author_grouped,d_authors_reverse,'e_authors',l_wts,)
##    ##    do_contour_plot(d_rmsds_overall,l_pdbs_resolution_grouped,d_resolutions_reverse,'h_resolutions',l_wts,)
##        do_contour_plot(d_rmsds_overall,l_pdbs_startingmodel_grouped,d_startingmodel_reverse,'g_start_model',l_wts,)
##        do_contour_plot(d_rmsds_overall,l_pdbs_countentities_grouped,d_countentities_reverse,'h_count_entities',l_wts,)
##        do_contour_plot(d_rmsds_overall,l_wts+l_mutants,d_seq,'0_wt_v_mutant',l_wts,) ## temporary

        ## this should only be done for heavy atoms, as side chain conformations are dependent on the starting model...
        if protein == 'HEWL':
            l = list(l_wts)
        elif protein == 'T4L':
            l = list(l_pdbs)
        d_rmsds_subset = plot_startmodel_vs_wts_and_derived_wts(
            l,
            d_properties,d_startingmodel,d_coordinates,
            d_rmsds_overall, d_mmCIF,
            )

        if d_rmsds_per_residue:
            
            plot_rmsds_per_residue_per_startingmodel(
                d_rmsds_per_residue,
                d_startingmodel, d_startingmodel_reverse,
                d_coordinates,
                d_rmsds_overall,
                l_pdbs,
                d_mmCIF,
                )

            plot_rmsds_per_residue(d_rmsds_per_residue)

    elif n_mutations_max > 0:

        l = [
            ## model
            '1RFP_A',
            ## wt
            '1UIH_A',
            ## mut
            '1FN5_A', '1FLQ_A', '1IOR_A', '1FLW_A', '1FLY_A', '1IOQ_A', '1FLU_A', '1IOT_A', '1IOS_A',
            ## model
            '6LYZ_A',
            ## wt
            '1AZF_A',
            ## mut
            '1LZD_A', '1LZE_A', '1LZG_A',
            ]

        d_rmsds_subset = calculate_rmsd_per_subset(l,d_properties,d_rmsds_overall,d_coordinates,)

        if protein == 'T4L':
            plot_start_model_vs_derived_wt_and_derived_mut(d_rmsds_subset)

    else:

        pass

    return


def plot_T4_mutants(
    d_chi1,d_rmsds_overall,d_pdbs,d_coordinates,
    l_wts,l_wts_cysfree,
    ):

    dist_sq_min = 25

    d_domains = {
        'large':range(73,165)+range(3,12),
        'small':range(12,73),
        }
    d_buried = {
        'buried':[
            ## Gly/Ala
            28,42,63,74,98,129,146,
            ## small domain
            6,7,10,27,29,33,46,58,66,
            ## large domain
            84,91,99,101,121,138,153,161, ## 117 only slightly buried, but sidechain fully buried!!!
            ## core
            78,87,102,111,133,149,152,
            ],
        'exposed':[
            ## Gly/Ala
            12,23,30,41,49,51,56,73,77,82,93,94,107,110,112,113,130,134,156,160,
            ## small domain
            1,2,3,4,5,8,9,11,13,14,15,16,17,18,19,20,21,22,24,25,26,31,32,34,35,36,37,38,39,40,43,44,45,47,48,50,52,53,55,57,59,60,61,62,64,65,67,68,69,70,71,72,75,76,
            54, ## slighlty buried (altlocs in 1lw9)
            ## large domain
            79,80,81,83,85,86,88,89,90,92,95,96,97,104,105,108,109,114,115,116,119,123,124,125,126,127,128,131,132,135,136,137,139,140,141,142,143,144,145,147,148,150,151,154,155,157,158,159,162,163,164,
            ## core
            100,103,106,118,120,
            122, ## slightly buried
            117, ## slightly buried, give rise to large chi1 changes
            ],
        }
    d_contacts = {
        2: [4, 68], 3: [4], 4: [2], 135: [45, 48], 136: [48], 140: [69], 141: [65], 14: [125], 143: [96], 144: [90], 19: [127], 21: [124, 126], 132: [45, 48], 158: [68], 162: [61], 41: [116], 44: [115, 116], 45: [116], 48: [114], 53: [55], 55: [53, 55], 61: [162], 65: [94], 68: [2, 93], 69: [93], 72: [88, 96], 76: [79], 79: [76], 142: [124], 88: [72], 89: [76], 90: [144], 92: [141], 93: [68, 69], 94: [65], 95: [141], 96: [72], 114: [48], 115: [44], 116: [41, 45], 117: [48], 122: [144], 124: [21], 125: [14], 126: [21], 127: [19]
        }

    d_mutations, d_comparisons = parse_T4_table(l_wts,l_wts_cysfree,)

##    for pdb in d_comparisons.keys():
##        if len(d_comparisons[pdb]) == 2:
##            l_comparisons = []
##            for i in range(2):
##                if d_comparisons[pdb][i] in d_mutations.keys():
##                    l_comparisons += d_mutations[mut]
##                    print pdb, mut, d_mutations[mut]
##            d_comparisons[pdb] = l_comparisons

    l_gnuplot_chi1 = []
    l_gnuplot_coord_alpha = []
    l_gnuplot_coord_heavy = []

    d_dist_v_rmsd = {}
    for k1 in ['coord_alpha','coord_heavy','chi1',]:
        d_dist_v_rmsd[k1] = {}
        for k2 in ['mut','wt',]:
            d_dist_v_rmsd[k1][k2] = {'mut_dist':[],'diff':[],}

    l_keys = [
##        ## wt
##        'large_buried',
##        'small_buried',
##        'large_exposed',
##        'small_exposed',
##
##        ## protprot_protsolv_mutprot
##        'large_buried_large',
##        'small_buried_large',
##        'large_exposed_large',
##        'small_exposed_large',
##        'large_buried_small',
##        'small_buried_small',
##        'large_exposed_small',
##        'small_exposed_small',

        ## v2 mut
        'large_large',
        'large_small',
        'small_large',
        'small_small',

        ## v2 wt
        'large','small',
        ]

    l_res_nos_mut = []

    d_count_rmsds = {}
    for key in l_keys:
        d_count_rmsds[key] = []

    ## residue rmsd distributions per pdb
    d_distributions = {'coord_alpha':{},'coord_heavy':{},'chi1':{},}

    ## coordinate difference distribution per residue
    d_distributions_residue = {}
    for k1 in ['wt','mut_distant','mut_vicinal',]:
        d_distributions_residue[k1] = {}
        for k2 in ['coord_alpha','coord_heavy','chi1',]:
            d_distributions_residue[k1][k2] = {}
            for res_no in range(1,len(ref_seq)+1):
                d_distributions_residue[k1][k2][res_no] = []

    d_jens = {}

##    print len( set(d_comparisons.keys()) - set(l_wts+l_wts_cysfree) )
##    print set(d_comparisons.keys()) & set(l_wts), l_wts
##    print set(d_comparisons.keys()) & set(l_wts_cysfree), l_wts_cysfree
##    stop

##    count_1l63 = 0
##    count_3lzm =0
##    for pdb in d_comparisons.keys():
##        if len(d_comparisons[pdb]) > 1:
##            continue
##        if d_comparisons[pdb] == ['1L63']:
##            count_1l63 += 1
##        elif d_comparisons[pdb] == ['3LZM']:
##            count_3lzm += 1
##        else:
##            print pdb, d_comparisons[pdb]
##    print count_1l63
##    print count_3lzm
##    stop

    for pdb in d_comparisons.keys():

        if d_mutations[pdb] == '':

            mutation = 'NA'
            res_no_mut = 'NA'
            res_name_mut = 'NA'
            domain_mut = 'NA'
            k_mut = 'wt'

        else:

            mutation = d_mutations[pdb]
            res_no_mut = int(mutation[1:-1])
            l_res_nos_mut += [res_no_mut]
            res_name_mut = d_123[mutation[-1]]
            k_mut = 'mut'

            if res_name_mut == 'GLY':
                atom_name_mut = 'CA'
            else:
                atom_name_mut = 'CB'

            if res_no_mut in d_domains['large']:
                domain_mut = 'large'
            else:
                domain_mut = 'small'

        ##
        ## loop over chains
        ##
        for pdb_chain in d_pdbs[pdb]:

            d_jens[pdb_chain] = {}
            for res_no in range(1,len(ref_seq)+1):
                d_jens[pdb_chain][res_no] = {}

            for key in ['coord_alpha','coord_heavy','chi1',]:
                d_distributions[key][pdb_chain] = []

            if d_mutations[pdb] != '':

                coord_mut = d_coordinates[pdb_chain][res_no_mut]['atoms'][atom_name_mut]['coord']

                if res_no_mut in d_contacts:
                    l_contacts = d_contacts[res_no_mut]
                else:
                    l_contacts = []

            else:

                coord_mut = 'NA'

            d_res_nos = {}
            for key in l_keys:
                d_res_nos[key] = []

            l_vicinal = []

            for res_no in d_coordinates[pdb_chain].keys():

                ## skip if residue is mutated
                if res_no == res_no_mut:
                    continue
                ## continue if ion or water
                elif res_no > 164:
                    continue

                if d_mutations[pdb] == '':
                    pass
                ## vicinal by sequence
                elif res_no in [res_no_mut-1,res_no_mut,res_no_mut+1,]:
                    bool_vicinal = True
                ## vicinal by space
                else:
                    bool_vicinal = False
                    for atom_name_mut in d_coordinates[pdb_chain][res_no_mut]['atoms'].keys():
                        if atom_name_mut in ['N','C','O',]:
                            continue
                        coord_mut = d_coordinates[pdb_chain][res_no_mut]['atoms'][atom_name_mut]['coord']
                        for atom_name in d_coordinates[pdb_chain][res_no]['atoms'].keys():
                            if atom_name in ['N','C','O',]:
                                continue
##                            if res_name != 'GLY' and atom_name == 'CA':
##                                continue
                            coord = d_coordinates[pdb_chain][res_no]['atoms'][atom_name]['coord']
                            dist_sq = sum((coord-coord_mut)**2)
                            if dist_sq < dist_sq_min:
                                bool_vicinal = True
                                break ## break loop over atom_names
                        if bool_vicinal == True:
                            break

                if protein == 'T4L':
                    if bool_vicinal == True:
                        l_vicinal += [res_no]
                
                if res_no in d_buried['buried']:
                    s_solvent = 'buried'
                else:
                    s_solvent = 'exposed'
                if res_no in d_domains['large']:
                    s_domain = 'large'
                else:
                    s_domain = 'small'

                ## mut
                if d_mutations[pdb] != '':
##                    d_res_nos['%s_%s_%s' %(s_domain,s_solvent,domain_mut,)] += [res_no]
                    d_res_nos['%s_%s' %(s_domain,domain_mut,)] += [res_no]

                ## wt
                elif d_mutations[pdb] == '' and res_no not in [54,97,]:
##                    d_res_nos['%s_%s' %(s_domain,s_solvent,)] += [res_no]
                    d_res_nos['%s' %(s_domain,)] += [res_no]

            ##
            ## only compare to selected structures
            ##
            for pdb_comparison in d_comparisons[pdb]:

                ##
                ## loop over chains of comparative models
                ##
                for pdb_comparison_chain in d_pdbs[pdb_comparison]:

                    if not pdb_comparison_chain in d_rmsds_overall[pdb_chain].keys():
                        print 'not compared???', pdb_chain, pdb_comparison_chain
                        continue

                    if len(d_pdbs[pdb_comparison]) > 1:
                        print pdb_comparison
                        print d_pdbs[pdb_comparison]
                        print pdb
                        stop

                    s_chi1 = '%s(%s) ' %(d_mutations[pdb],pdb_chain.replace('_',':'),)
                    s_coord_alpha = '%s(%s) ' %(d_mutations[pdb],pdb_chain.replace('_',':'),)
                    s_coord_heavy = '%s(%s) ' %(d_mutations[pdb],pdb_chain.replace('_',':'),)

                    tv1 = d_rmsds_overall[pdb_chain][pdb_comparison_chain]['tv1']
                    rm = d_rmsds_overall[pdb_chain][pdb_comparison_chain]['rm']
                    tv2 = d_rmsds_overall[pdb_chain][pdb_comparison_chain]['tv2']

                    if d_mutations[pdb] != '':
                        for res_no in range(1,1+len(ref_seq)):

                            if not res_no in d_coordinates[pdb_comparison_chain].keys():
                                continue
                            if not res_no in d_coordinates[pdb_chain].keys():
                                continue

                            for atom_name in d_coordinates[pdb_chain][res_no]['atoms'].keys():

                                if not atom_name in d_coordinates[pdb_comparison_chain][res_no]['atoms'].keys():
                                    continue

                                coord = d_coordinates[pdb_chain][res_no]['atoms'][atom_name]['coord']
                                l_dist = []
                                for atom_name_mut in d_coordinates[pdb_chain][res_no_mut]['atoms'].keys():
                                    coord_mut = d_coordinates[pdb_chain][res_no_mut]['atoms'][atom_name_mut]['coord']
                                    dist = math.sqrt(sum((coord-coord_mut)**2))
                                    l_dist += [dist]
                                ## distance to mutation
                                dist_average_mut = sum(l_dist)/len(l_dist)

                                coord_comparison = d_coordinates[pdb_comparison_chain][res_no]['atoms'][atom_name]['coord']
                                coord_comparison = numpy.dot(coord_comparison-tv1,rm)+tv2
                                ## mut v comparison pdb
                                coord_diff = math.sqrt(sum((coord-coord_comparison)**2))

                                d_jens[pdb_chain][res_no][atom_name] = [coord_diff,dist_average_mut,]

                    bool_outlier = False
                    bool_outlier_coord_alpha = False
                    bool_outlier_coord_heavy = False

                    ##
                    ## loop over buried/exposed residues in large/small domain
                    ##
                    for key in l_keys:
                        s_outliers_exposed = ''
                        s_outliers_buried = ''
                        l_res_nos = d_res_nos[key]
                        l_chi1_diff = []
                        l_coord_diff_alpha = []
                        l_coord_diff_heavy = []

    ##                    l_coords_heavy1 = []
    ##                    l_coords_heavy2 = []

                        for res_no in l_res_nos:

                            if k_mut == 'mut':
                                if res_no in l_vicinal:
                                    k2_mut = 'mut_vicinal'
                                else:
                                    k2_mut = 'mut_distant'
                            else:
                                k2_mut = 'wt'

                            ## unobserved residue
                            if not res_no in d_coordinates[pdb_chain].keys():
                                continue
                            if not res_no in d_coordinates[pdb_comparison_chain].keys():
                                continue

                            if (
                                res_no in d_chi1.keys()
                                and
                                pdb_chain in d_chi1[res_no].keys()
                                and
                                pdb_comparison_chain in d_chi1[res_no].keys()
                                ):
                                chi1_diff = abs(d_chi1[res_no][pdb_chain]-d_chi1[res_no][pdb_comparison_chain])
                                if chi1_diff > 180:
                                    chi1_diff = 360-chi1_diff
                            else:
                                chi1_diff = None

                            res_name1 = d_coordinates[pdb_chain][res_no]['res_name']
                            res_name2 = d_coordinates[pdb_comparison_chain][res_no]['res_name']

                            ## continue if modified residue
                            if res_name1 not in ['MSE','CME',]+d_321.keys():
                                continue
                            if res_name2 not in ['MSE','CME',]+d_321.keys():
                                continue

                            l_diff_residue = []
                            for atom_name in d_coordinates[pdb_chain][res_no]['atoms'].keys():
                                ## atom not present in comparative pdb
                                if not atom_name in d_coordinates[pdb_comparison_chain][res_no]['atoms'].keys():
                                    continue
                                ## only compare backbone atoms if mutation
                                if res_name1 != res_name2:
                                    if atom_name not in ['N','CA','C','O',]:
                                        continue
                                if bool_exclusion_high_temp_factors == True:
                                    if d_coordinates[pdb_chain][res_no]['atoms'][atom_name]['tempFactor'] > max_bfactor:
                                        continue
                                    if d_coordinates[pdb_comparison_chain][res_no]['atoms'][atom_name]['tempFactor'] > max_bfactor:
                                        continue
                                if atom_name[0] == 'H':
                                    stop
                                coord = d_coordinates[pdb_chain][res_no]['atoms'][atom_name]['coord']
                                coord_comparison = d_coordinates[pdb_comparison_chain][res_no]['atoms'][atom_name]['coord']
                                coord_comparison = numpy.dot(coord_comparison-tv1,rm)+tv2
                                coord_diff = math.sqrt(sum((coord-coord_comparison)**2))

                                l_diff_residue += [coord_diff]

                                if k_mut == 'mut':

                                    mut_dist = math.sqrt(sum((coord_mut-coord)**2))

                                    if atom_name == 'CA':
                                        d_dist_v_rmsd['coord_alpha'][k_mut]['mut_dist'] += [mut_dist]
                                        d_dist_v_rmsd['coord_alpha'][k_mut]['diff'] += [coord_diff]
                                        d_distributions['coord_alpha'][pdb_chain] += [coord_diff]
                                        d_distributions_residue[k2_mut]['coord_alpha'][res_no] += [coord_diff]
                                    elif atom_name == 'CB':
                                        d_dist_v_rmsd['chi1'][k_mut]['mut_dist'] += [mut_dist]
                                        if chi1_diff != None:
                                            d_dist_v_rmsd['chi1'][k_mut]['diff'] += [chi1_diff]
                                    d_dist_v_rmsd['coord_heavy'][k_mut]['mut_dist'] += [mut_dist]
                                    d_dist_v_rmsd['coord_heavy'][k_mut]['diff'] += [coord_diff]

                                l_coord_diff_heavy += [coord_diff]

    ##                            l_coords_heavy1 += [coord]
    ##                            l_coords_heavy2 += [coord_comparison]

                                if atom_name == 'CA':
                                    l_coord_diff_alpha += [coord_diff]

                            rmsd_residue = statistics.do_rmsd(l_diff_residue)
                            if chi1_diff != None:
                                d_distributions['chi1'][pdb_chain] += [chi1_diff]
                                d_distributions_residue[k2_mut]['chi1'][res_no] += [chi1_diff]
                            d_distributions['coord_heavy'][pdb_chain] += [rmsd_residue]
                            d_distributions_residue[k2_mut]['coord_heavy'][res_no] += [rmsd_residue]

                            if chi1_diff > 30 and key[6:13] == 'exposed' and 'vicinal' not in key:
                                s_outliers_exposed += 'single %s %s %s %s %s %s %s %s %s\n' %(ref_seq[res_no-1], res_no, mutation, pdb, pdb_comparison_chain, round(chi1_diff,1), round(d_chi1[res_no][pdb_chain],1), round(d_chi1[res_no][pdb_comparison_chain],1), key)
                            if chi1_diff > 15 and key[6:12] == 'buried' and 'vicinal' not in key:
                                s_outliers_buried += 'single %s %s %s %s %s %s %s %s %s' %(ref_seq[res_no-1], res_no, mutation, pdb, pdb_comparison_chain, round(chi1_diff,1), round(d_chi1[res_no][pdb_chain],1), round(d_chi1[res_no][pdb_comparison_chain],1), key)

                            if chi1_diff != None:
                                l_chi1_diff += [chi1_diff]

                        ## end of loop over res_nos

                        if len(l_chi1_diff) == 0:

                            rmsd_chi1 = 'NA'
                            rmsd_coord_alpha = 'NA'
                            rmsd_coord_heavy = 'NA'

                        else:
                                
                            chi1_diff_average = sum(l_chi1_diff)/len(l_chi1_diff)
                            rmsd_coord_overall = float(d_rmsds_overall[pdb_chain][pdb_comparison_chain]['rmsd'])
                            rmsd_coord_alpha = statistics.do_rmsd(l_coord_diff_alpha)
                            rmsd_coord_heavy = statistics.do_rmsd(l_coord_diff_heavy)
                            rmsd_chi1 = statistics.do_rmsd(l_chi1_diff)

        ##                    d_count_rmsds[key] += [len(l_chi1_diff)]
                            d_count_rmsds[key] += [len(l_coord_diff_heavy)]

        ##                        rmsd_coord_heavy = instance_geometry.superpose(l_coords_heavy1,l_coords_heavy2,)

                            if (
                                (rmsd_chi1 > 10 and key[6:13] == 'exposed' and len(key) > 13 and 'vicinal' not in key and 'different' not in key)
                                ):
                                if pdb_chain not in [
        ##                                '1LLH_A', ## Leu79 chi2 180 rot
                                    ]:
                                    print 'multiple exposed', pdb_chain, mutation, pdb_comparison_chain, key, len(l_res_nos), len(l_chi1_diff), round(rmsd_chi1,1)
                                    print s_outliers_exposed
                                    bool_outlier = True

                            if (
                                (rmsd_chi1 > 10 and key[6:12] == 'buried' and len(key) > 13 and 'vicinal' not in key)
                                ):
                                if pdb_chain not in [
        ##                                '200L_A','1TLA_A','1QUG_A','1G1V_A','1LLH_A','3C7Y_A',
        ##                                '224L_A', ## Thr21 180 rot
                                    ]:
                                    print 'multiple buried', pdb_chain, mutation, pdb_comparison_chain, key, len(l_res_nos), len(l_chi1_diff), round(rmsd_chi1,1)
                                    print s_outliers_buried
                                    bool_outlier = True

                            if (
                                (rmsd_coord_alpha > 0.30 and 'vicinal' not in key)
                                or
                                (rmsd_coord_alpha > 0.25 and key[6:12] == 'buried' and 'vicinal' not in key)
                                ):
                                bool_outlier_coord_alpha = True

                            if (
                                (rmsd_coord_heavy > 0.35)
                                ):
                                bool_outlier_coord_heavy = True

                        s_chi1 += '%s ' %(rmsd_chi1)
                        s_coord_alpha += '%s ' %(rmsd_coord_alpha)
                        s_coord_heavy += '%s ' %(rmsd_coord_heavy)

                        continue

                    ## end of loop over group of residues

                    if bool_outlier == True:
                        s_chi1 += '%s(%s-%s)' %(d_mutations[pdb],pdb_chain.replace('_',':'),pdb_comparison_chain.replace('_',':'),)
    ##                    s_chi1 += '%s' %(d_mutations[pdb],)
                    if bool_outlier_coord_alpha == True:
                        s_coord_alpha += '%s(%s-%s)' %(d_mutations[pdb],pdb_chain.replace('_',':'),pdb_comparison_chain.replace('_',':'),)
    ##                    s_coord_alpha += '%s' %(d_mutations[pdb],)
                    if bool_outlier_coord_heavy == True:
                        s_coord_heavy += '%s(%s-%s)' %(d_mutations[pdb],pdb_chain.replace('_',':'),pdb_comparison_chain.replace('_',':'),)
    ##                    s_coord_heavy += '%s' %(d_mutations[pdb],)
                    s_chi1 += '\n'
                    s_coord_alpha += '\n'
                    s_coord_heavy += '\n'
                    l_gnuplot_chi1 += [s_chi1]
                    l_gnuplot_coord_alpha += [s_coord_alpha]
                    l_gnuplot_coord_heavy += [s_coord_heavy]

                    continue

                ## end of loop over pdb comparison

    if protein == 'HEWL':
        for pdb in list(d_jens.keys()):
            if d_jens[pdb][1] == {}:
                print pdb
                del d_jens[pdb]
    fd = open('d_jens.txt','w')
    fd.write(str(d_jens))
    fd.close()
    import pickle
    fd = open('d_jens.pickle','w')
    pickle.dump(d_jens,fd)
    fd.close()
    stop_jens_done

    for key in l_keys:
        print key, sum(d_count_rmsds[key])/float(len(d_count_rmsds[key])), min(d_count_rmsds[key]), max(d_count_rmsds[key])
##    stop

    ##
    ##
    ##
    print 'histograms per residue'
##    for key in d_distributions_residue.keys():
##        print key, len(d_distributions_residue[key]['coord_heavy'][92])
##    for res_no in range(1,161):
##        print res_no, len(d_distributions_residue['mut_vicinal']['coord_heavy'][res_no])
##    stop
    for key,xtic_freq,xtic_max,ymax, in [
        ['coord_alpha',0.1,5.0,150,],
        ['chi1',3.0,180.0,170,],
        ['coord_heavy',0.1,5.0,180,],
##        ['coord_alpha',0.1,4.0,],
        ]:
##        d_xtics = {'wt':{},'mut':{},}
        d_xtics = {'wt':{},'mut_distant':{},'mut_vicinal':{},}
        l_res_nos = []
        lines_plot_probability = []
        for res_no in range(1,161):
            if res_no in [54,97,]:
                continue
            if ref_seq[res_no-1] in ['GLY','ALA',] and key == 'chi1':
                continue
            for col_header in d_xtics.keys():
                d_xtics[col_header][res_no] = {}
                for x in range(int(xtic_max/xtic_freq)):
                    xtic = x*xtic_freq
                    d_xtics[col_header][res_no][xtic] = 0

            l_rmsds_residue_distant = d_distributions_residue['mut_distant'][key][res_no]
            l_rmsds_residue_vicinal = d_distributions_residue['mut_vicinal'][key][res_no]
            if len(l_rmsds_residue_vicinal) > 1:
                p = statistics.Kolmogorov_Smirnov(l_rmsds_residue_distant,l_rmsds_residue_vicinal,verbose=False,)
##                mean_distant,mean_vicinal,stderr,p = statistics.twosamplettest(l_rmsds_residue_distant,l_rmsds_residue_vicinal,verbose=False,)
                lines_plot_probability += ['%s%s %s\n' %(ref_seq[res_no-1],res_no,p,)]
                mean_distant = sum(l_rmsds_residue_distant)/len(l_rmsds_residue_distant)
                mean_vicinal = sum(l_rmsds_residue_vicinal)/len(l_rmsds_residue_vicinal)
                if p < 0.05 and mean_vicinal > mean_distant:
##                if p < 0.10 and mean_vicinal > mean_distant:
                    print '***', key, ref_seq[res_no-1], res_no, round(p,3), round(mean_distant,3), round(mean_vicinal,3)
                    l_res_nos += [res_no]
                else:
                    print key, ref_seq[res_no-1], res_no, p
##            else:
##                print key, ref_seq[res_no-1], res_no, 'N/A'

            d_count = {'wt':0,'mut_vicinal':0,'mut_distant':0,}
            for key_mut in ['wt','mut_distant','mut_vicinal',]:
                l_rmsds_residue = d_distributions_residue[key_mut][key][res_no]
                for rmsd in l_rmsds_residue:
                    d_xtics[key_mut][res_no][rmsd-rmsd%xtic_freq] += 1
                    d_count[key_mut] += 1

##            if d_count['wt'] == 0:
##                print ref_seq[res_no-1], key
##                continue

            l = ['xtics muts_distant muts_vicinal wts\n']
            for x in range(int(xtic_max/xtic_freq)):
                xtic = x*xtic_freq
                count_mut_distant = d_xtics['mut_distant'][res_no][xtic]
                count_mut_vicinal = d_xtics['mut_vicinal'][res_no][xtic]
##                count_wt = d_xtics['wt'][res_no][xtic]*d_count['mut_distant']/float(d_count['wt'])
                if count_mut_distant > ymax:
                    print key,res_no,count_mut,count_wt
                    stop
##                if count_wt > ymax:
##                    print key,res_no,count_mut,count_wt
##                    stop
                l += [
##                    '%s-%s %s %s %s\n' %(
                    '%s-%s %s %s\n' %(
                        xtic,xtic_freq*(xtic/xtic_freq+1),
                        count_mut_distant,
                        count_mut_vicinal,
##                        count_wt,
                        )
                    ]
            prefix = 'distribution_%s_%s_%03i%3s' %(
                dist_sq_min,
                key,
                res_no,ref_seq[res_no-1],
                )
            fd = open('%s.gnuplotdata' %(prefix,),'w')
            fd.writelines(l)
            fd.close()
            gnuplot.histogram(
                prefix,
                title='residue rmsd distribution, %s%s,  %s' %(
                    ref_seq[res_no-1], res_no, key,
                    ),
                xlabel = 'residue rmsd',
                ylabel = 'count',
                ymax = ymax,
##                bool_stacked = True,
##                bool_remove = False,
                )
            os.system('mv %s.png rmsd_distributions/%s/.' %(prefix,dist_sq_min,))
            print prefix

        prefix = 'prob_hist_%s_%s' %(dist_sq_min,key,)
        fd = open('%s.gnuplotdata' %(prefix,),'w')
        fd.writelines(lines_plot_probability)
        fd.close()
        gnuplot.histogram(
            prefix,
            title='probability that mutation does not have an effect, %s, cutoff=%iA' %(
                key,int(math.sqrt(dist_sq_min)),
                ),
            xlabel = 'residue',
            ylabel = 'probability',
            ymax = 1,
            bool_remove = False,
            )

        print l_res_nos
    stop_end

    ##
    ## residue rmsd distribution plots
    ##
    print 'histograms per pdb'
    for key,xtic_freq,xtic_max, in [
##        ['coord_heavy',0.1,5.0,],
        ['chi1',3.0,180.0,],
        ['coord_alpha',0.1,4.0,],
        ]:
        d_xtics = {}
        for pdb_chain in l_wts+l_wts_cysfree+['1L90_A',]:

            ## skip
            if not pdb_chain[:4] in d_comparisons.keys():
                continue

            d_xtics[pdb_chain] = {}
            for x in range(xtic_max/xtic_freq):
                xtic = x*xtic_freq
                d_xtics[pdb_chain][xtic] = 0

            l_rmsds_residue = d_distributions[key][pdb_chain]
            for rmsd in l_rmsds_residue:
                print key, pdb_chain, rmsd,xtic_freq
                d_xtics[pdb_chain][rmsd-rmsd%xtic_freq] += 1

        for l in [l_wts,l_wts_cysfree,]:
            for pdb_chain in l:
                ## skip
                if not pdb_chain[:4] in d_comparisons.keys():
                    continue
                count = 0
                for x in range(xtic_max/xtic_freq):
                    xtic = x*xtic_freq
                    d_xtics[pdb_chain][xtic] /= (float(len(l))-1)
                    count += d_xtics[pdb_chain][xtic]
                if count < 159 or count > 165:
                    print pdb_chain
                    print count
                    print len(l)
                    stop
            
        for pdb_chain in d_distributions[key].keys():
            if pdb_chain in l_wts+l_wts_cysfree:
                continue
            d_xtics[pdb_chain] = {}
            for x in range(xtic_max/xtic_freq):
                xtic = x*xtic_freq
                d_xtics[pdb_chain][xtic] = 0

            l_rmsds_residue = d_distributions[key][pdb_chain]
            count = 0
            for rmsd in l_rmsds_residue:
                d_xtics[pdb_chain][rmsd-rmsd%xtic_freq] += 1
                count += 1
            if count not in range(160,165):
                print pdb_chain
                print count
                print len(l_rmsds_residue)
                stop
            wt = ''.join(d_comparisons[pdb_chain[:4]])+'_A'
            l = ['xtics %s %s\n' %(pdb_chain,wt,)]
            for x in range(xtic_max/xtic_freq):
                xtic = x*xtic_freq
                count_mut = d_xtics[pdb_chain][xtic]
                count_wt = d_xtics[wt][xtic]
                l += [
                    '%s-%s %s %s\n' %(
                        xtic,xtic_freq*(xtic/xtic_freq+1),
                        count_mut,
                        count_wt,
                        )
                    ]
            prefix = 'distribution_%s_%s' %(pdb_chain,key,)
            fd = open('%s.gnuplotdata' %(prefix,),'w')
            fd.writelines(l)
            fd.close()
            gnuplot.histogram(
                prefix,
                title='residue rmsd distribution, %s (%s), %s' %(
                    pdb_chain.replace('_',''), d_mutations[pdb_chain[:4]], key,
                    ),
                xlabel = 'residue rmsd',
                ylabel = 'count',
                ymax = 160,
                )
            os.system('mv %s.png rmsd_distributions/.' %(prefix))
            print prefix
    stop_done


    ##
    ## sample wts as muts
    ##
    for pdb in l_wts+l_wts_cysfree:
        ## only compare most similar pdbs
        if pdb not in [
            '4LZM','219L',
            ]:
            continue
        print 'wt/wt*', pdb
        ## pdb was skipped previously; e.g. because of ligand
        if not pdb in d_comparisons.keys():
            continue
        for pdb_chain in d_pdbs[pdb]:
            for pdb_comparison in d_comparisons[pdb]:
                ## only do 1v1 comparison
                if pdb_comparison not in ['3LZM','1L63',]:
                    continue
                for pdb_comparison_chain in d_pdbs[pdb_comparison]:
                    if '1LYD_A' in [pdb_chain,pdb_comparison_chain,]:
                        continue
##                    if '3FA0_A' in [pdb_chain,pdb_comparison_chain,]:
##                        continue
                    if '2LZM_A' in [pdb_chain,pdb_comparison_chain,]:
                        continue
##                    if '5LZM_A' in [pdb_chain,pdb_comparison_chain,]: ## 3fa0 v 5lzm Arg80 CB
##                        continue
                    tv1 = d_rmsds_overall[pdb_chain][pdb_comparison_chain]['tv1']
                    rm = d_rmsds_overall[pdb_chain][pdb_comparison_chain]['rm']
                    tv2 = d_rmsds_overall[pdb_chain][pdb_comparison_chain]['tv2']
                    for res_no_mut in l_res_nos_mut: ## when regression calculation
##                    for res_no_mut in set(l_res_nos_mut): ## when plotting data points
                        res_name_mut = ref_seq[res_no_mut-1]
                        if res_name_mut == 'GLY':
                            atom_name_mut = 'CA'
                        else:
                            atom_name_mut = 'CB'
                        coord_mut = d_coordinates[pdb_chain][res_no_mut]['atoms'][atom_name_mut]['coord']
                        for res_no in d_coordinates[pdb_chain].keys():
                            if res_no == res_no_mut:
                                continue
                            if not res_no in d_coordinates[pdb_comparison_chain].keys():
                                continue
                            for atom_name in d_coordinates[pdb_chain][res_no]['atoms'].keys():
                                if bool_exclusion_high_temp_factors == True:
                                    if d_coordinates[pdb_chain][res_no]['atoms'][atom_name]['tempFactor'] > max_bfactor:
                                        continue
                                    if d_coordinates[pdb_comparison_chain][res_no]['atoms'][atom_name]['tempFactor'] > max_bfactor:
                                        continue
                                coord = d_coordinates[pdb_chain][res_no]['atoms'][atom_name]['coord']

                                mut_dist = math.sqrt(sum((coord_mut-coord)**2))

                                coord_comparison = d_coordinates[pdb_comparison_chain][res_no]['atoms'][atom_name]['coord']
                                coord_comparison = numpy.dot(coord_comparison-tv1,rm)+tv2
                                coord_diff = math.sqrt(sum((coord-coord_comparison)**2))
##                                if coord_diff > 0.85:
##                                    print pdb_chain, pdb_comparison_chain, res_no, ref_seq[res_no-1], atom_name
##                                    print 'coord_diff', coord_diff
##                                    stop1
##                                if mut_dist < 2.8:
##                                    print pdb_chain, res_no, ref_seq[res_no-1], atom_name
##                                    print pdb_chain, res_no_mut, res_name_mut, atom_name_mut
##                                    print 'mut_dist', mut_dist
##                                    stop2
##                                if mut_dist > 50:
##                                    print pdb_chain, res_no, ref_seq[res_no-1], atom_name
##                                    print pdb_chain, res_no_mut, res_name_mut, atom_name_mut
##                                    print 'mut_dist', mut_dist
##                                    stop3
                                if coord_diff > 5:
                                    print pdb_chain, pdb_comparison_chain, res_no, atom_name
                                    print coord_diff
                                    stop
                                d_dist_v_rmsd['coord_heavy']['wt']['mut_dist'] += [mut_dist]
                                d_dist_v_rmsd['coord_heavy']['wt']['diff'] += [coord_diff]


    ##
    ## mut dist v coord diff
    ##
    print 'mut dist v coord diff'
    for k,d,k2, in [
        ['coord_heavy',d_dist_v_residue_rmsd,'rmsd',],
        ['coord_heavy',d_dist_v_rmsd,'rmsd',],
        ['chi1',d_dist_v_rmsd,'diff',],
        ]:
        
        prefix = 'dist_v_%s_%s' %(k,k2,)
        s_plot = 'set xlabel "distance from alpha carbon atom of mutated residue"\nset ylabel "coordinate difference"\nplot [0:][0:]'
        s_plot += ' "dist_v_%s_wt.gnuplotdata" lc 9 pt 6 ps 2 t "", ' %(k)
        s_plot += ' "dist_v_%s_mut.gnuplotdata" lc 0 pt 7 ps 2 t "", ' %(k)
        s_plot += ' g(x) lt 1 lc 0 lw 20 t "mutants", f(x) lt 1 lc 9 lw 20 t "wt/wt*"\n'
        l_functions = []
        ## wt / f(x)
        a,b = statistics.do_regression(d[k]['wt']['mut_dist'],d[k]['wt'][k2],)
##        a = 0.149216785185
##        b = 0.000891537742289
        l_functions += ['f(x) = %s+%s*x' %(a,b,)]
        ## mutants / g(x)
        a,b = statistics.do_regression(d[k]['mut']['mut_dist'],d[k]['mut'][k2],)
        l_functions += ['g(x) = %s+%s*x' %(a,b,)]
        for k_mut in [
            'mut',
            'wt',
            ]:
            lines = []
            for i in range(len(d_dist_v_rmsd[k][k_mut]['mut_dist'])):
                lines += ['%s %s\n' %(d_dist_v_rmsd[k][k_mut]['mut_dist'][i],d_dist_v_rmsd[k][k_mut]['diff'][i],)]
            fd = open('%s_%s.gnuplotdata' %(prefix,k_mut),'w')
            fd.writelines(lines)
            fd.close()
        del d_dist_v_rmsd[k]['wt']
        del d_dist_v_rmsd[k]['mut']
        gnuplot.scatter_plot_2d(
            prefix, s_plot = s_plot,
            l_functions = l_functions,
            )

    ##
    ## calculate regression lines and confidence bands
    ##
    d_functions = {}
    for prop in [
        'coord_alpha',
        'coord_heavy',
##        'chi1',
        ]:
        l1 = [] ## large
        l2 = [] ## small
        fd = open('/local/tc/HEWL_variability/%s_T4L_%s.gnuplotdata' %(prop,suffix_exclusion,),'r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            for i in range(0,12,2,):
                if line.split()[i] != 'NA' and line.split()[i+2] != 'NA':
                    l1 += [float(line.split()[i+1])]
                    l2 += [float(line.split()[i+1+1])]
        t_crit = 1.960 ## table B3
        f,g,h = statistics.do_confidence_bands(l1,l2,t_crit,)
        d_functions[prop] = {'f':f,'g':g,'h':h,}


    for prefix,l_gnuplot,xmax,prop in [
        ## chi1
        [
##            'chi1',l_gnuplot_chi1,30,'{/Symbol c_1}',
            'chi1',l_gnuplot_chi1,55,'{/Symbol c_1}',
            ],
        ## alpha rmsd
        [
##            'coord_alpha',l_gnuplot_coord_alpha,0.4,'C_{/Symbol a}',
            'coord_alpha',l_gnuplot_coord_alpha,0.45,'C_{/Symbol a}',
            ],
        ## heavy rmsd
        [
##            'coord_heavy',l_gnuplot_coord_heavy,0.6,'heavy atom',
            'coord_heavy',l_gnuplot_coord_heavy,0.95,'heavy atom',
            ],
        ]:

        fn_prefix = '%s_%s_%s_v2' %(prefix,protein,suffix_exclusion,)
        fd = open('%s.gnuplotdata' %(fn_prefix), 'w')
        fd.writelines(l_gnuplot)
        fd.close()

        if prefix == 'chi1':
            l_functions = []
        else:
            l_functions = [d_functions[prefix]['f'],d_functions[prefix]['g'],d_functions[prefix]['h'],]

##        l_columns = [
##            ## mutations in the large domain
##                ## buried residues
##            {
##                'x':l_keys.index('large_buried_large')+2,'y':l_keys.index('small_buried_large')+2,
##    ##            't':'buried, mutation buried in large domain',
##                'label':len(l_keys)+2,
##                },
##                ## exposed residues
##            {
##                'x':l_keys.index('large_exposed_large')+2,'y':l_keys.index('small_exposed_large')+2,
##    ##            't':'exposed, mutation buried in large domain',
##                'label':len(l_keys)+2,
##                },
##            ## mutations in the small domain
##                ## buried residues
##            {
##                'x':l_keys.index('large_buried_small')+2,'y':l_keys.index('small_buried_small')+2,
##    ##            't':'buried, mutation buried in small domain',
##                'label':len(l_keys)+2,
##                },
##                ## exposed residues
##            {
##                'x':l_keys.index('large_exposed_small')+2,'y':l_keys.index('small_exposed_small')+2,
##    ##            't':'exposed, mutation buried in small domain',
##                'label':len(l_keys)+2,
##                },
##            ## wts
##                ## buried
##            {
##                'x':l_keys.index('large_buried')+2,'y':l_keys.index('small_buried')+2,
##    ##            't':'wt, buried',
##                'label':len(l_keys)+2,
##                },
##                ## exposed
##            {
##                'x':l_keys.index('large_exposed')+2,'y':l_keys.index('small_exposed')+2,
##    ##            't':'wt, exposed',
##                'label':len(l_keys)+2,
##                },
##            ]
##        l_pointtypes = [7,6,7,6,7,6,]
##        l_colors = [
##            'ff0000','ff0000',
##            '00ff00','00ff00',
##            ## wts
##            '0000ff','0000ff',
##            ]

        l_columns = [
            ## mutations in the large domain
            {
                'x':l_keys.index('large_large')+2,'y':l_keys.index('small_large')+2,
##                'label':len(l_keys)+2,
                },
            ## mutations in the small domain
            {
                'x':l_keys.index('large_small')+2,'y':l_keys.index('small_small')+2,
##                'label':len(l_keys)+2,
                },
            ## wts
            {
                'x':l_keys.index('large')+2,'y':l_keys.index('small')+2,
##                'label':len(l_keys)+2,
                },
            ]
        l_pointtypes = [7,7,7,]
        l_colors = [
            'ff0000',
            '00ff00',
            ## wts
            '0000ff',
            ]

        gnuplot.scatter_plot_2d(
            fn_prefix,
            l_functions = l_functions,
            l_columns = l_columns,
            l_pointtypes = l_pointtypes,
            l_colors = l_colors,
    ##        col_label = len(l_keys)+2,
            xmin = 0, ymin = 0,
            xmax = xmax, ymax = xmax,
            key_vert_pos = 'left',
            size = 'square',
            title = '%s RMSD' %(prop),
            xlabel = '%s RMSD, large domain' %(prop),
            ylabel = '%s RMSD, small domain' %(prop),
            bool_remove = False,
            )
    stop_fin1




    for col1,col2,title_suffix in [
        ## same v different
        [3,5,'domain',],[4,6,'domain',],
        ## buried v exposed
        [1,2,'exposure',],[3,4,'exposure',],[5,6,'exposure',],
        ## vicinal v vicinal crystal
        [1,7,'',],[2,7,'',],
        ## vicinal exposed v domains
        [1,3,'vicinal exposed',],[1,5,'vicinal exposed',],
        ## vicinal buried v domains
        [2,4,'vicinal buried',],[2,6,'vicinal buried',],
        ]:
        gnuplot.scatter_plot_2d(
            prefix,
            suffix = '%s%s' %(col1,col2,),
            colx = col1,
            coly = col2,
            col_label = len(l_keys)+1,
            xmin = 0, ymin = 0,
            title = 'chi1 RMSD %s' %(title_suffix),
            xlabel = l_keys[col1-1].replace('_',' '),
            ylabel = l_keys[col2-1].replace('_',' '),
            )

    stop_fin

    return


def parse_T4_table(l_wts,l_wts_cysfree,):

    fd = open('table_%s.txt' %(protein),'r')
    lines = fd.readlines()
    fd.close()

    d_HEWL_representative = {
        ## 193l (1995) is 1.33A in P 43 21 2
        ## 2vb1 (2007) is 0.65A in P 1
        ## 2lzt (1989) is 1.97A in P 1
        'C 1 2 1':'193L', ## 1ps5 from 1iee from 1lsa = p 43 21 2 = 193l
        'P 1':'2LZT', ## 1v7s, 1v7t, 3lzt, 4lzt from 2lzt
        'P 1 21 1':'193L', ## 1b2k and 1jj3 from 193l = p43 21 2
        'P 21 21 21':'2LZT', ## 1bgi start model
        'P 43 21 2':'193L', ##
        'P 61 2 2':'193L', ## 2fbb from 1dpw = p 43 21 2 = 193l
        }

    l_headers = [s.strip() for s in lines[0].split('\t')]

    d_mutations = {}
    d_comparisons = {}
    for line in lines[1:]:
        l = line.split('\t')
        d_row = {}
        for i in range(len(l_headers)):
            header = l_headers[i]
            s = l[i].strip()
            d_row[header] = s
        pdb = d_row['ID']
        if protein == 'T4L':
            n_mutations = d_row['mutation(s)'].count('+')
        elif protein == 'HEWL':
            if d_row['mutation(s)'] == 'wt':
                n_mutations = 0
            else:
                n_mutations = d_row['mutation(s)'].count('+')+1
        if protein == 'T4L':
            ## skip if not P 32 2 1
            if d_row['SG'] not in [
                'P 32 2 1', ## T4L
                ]:
                continue
        ## skip if multiple mutations in core
        if 'CORE' in d_row['mutation(s)']:
            continue
        ## skip if more than 2 mutations
        if n_mutations > 2:
            continue
        ## skip if ligands present
        if d_row['lig'] != '':
            continue
        ## skip if mutation to cys
        if 'C+' in d_row['mutation(s)']:
            continue
        if 'disulf' in d_row['mutation(s)']:
            continue
        ## skip if mutation to non-std res
        if 'X' in d_row['mutation(s)']:
            continue
        ## skip if high pressure or unusual salt
##        if d_row['comment'] != '':
        if 'ATM' in d_row['comment']:
            continue

        if n_mutations == 1:
##            mutation = d_row['mutation(s)']
##            if not mutation in d_mutations.keys():
##                d_mutations[mutation] = []
##            d_mutations[mutation] += pdb
            if protein == 'T4L':
                d_mutations[pdb] = d_row['mutation(s)'][d_row['mutation(s)'].index('+')+1:]
            elif protein == 'HEWL':
                d_mutations[pdb] = d_row['mutation(s)']
        elif n_mutations == 0:
            d_mutations[pdb] = ''
##            if pdb[:4] == '2LZM':
##                continue
##            continue
        elif n_mutations == 2:
            continue

        if n_mutations == 0:
            if protein == 'T4L':
                if '*' in d_row['mutation(s)']:
                    l_pdbs_compare = list(l_wts_cysfree)
                else:
                    l_pdbs_compare = list(l_wts)
            elif protein == 'HEWL':
##                l_pdbs_compare = ['2VB1']
                l_pdbs_compare = [d_HEWL_representative[d_row['SG']]]
            for i_pdb in range(len(l_pdbs_compare)):
                l_pdbs_compare[i_pdb] = l_pdbs_compare[i_pdb][:4]
            l_pdbs_compare = list(set(l_pdbs_compare))
            ## don't compare protein to itself
            if protein == 'T4L':
                l_pdbs_compare.remove(pdb)
        elif d_row['SM'] != 'N/A':
            l_pdbs_compare = [d_row['SM']]
            if protein == 'T4L':
                if d_row['mutation(s)'][:d_row['mutation(s)'].index('+')] == 'wt*':
                    if l_pdbs_compare not in [['1L63',],['1L90',],]:
                        print d_row['mutation(s)'][:d_row['mutation(s)'].index('+')], pdb, l_pdbs_compare
                else:
                    if l_pdbs_compare not in [['2LZM'],['3LZM',],]:
                        print d_row['mutation(s)'][:d_row['mutation(s)'].index('+')], pdb, l_pdbs_compare
##            elif protein == 'HEWL':
##                if l_pdbs_compare == ['6LYZ']:
##                    l_pdbs_compare = '2LZT'
        else:
            if n_mutations == 1:
                if protein == 'T4L':
                    if d_row['mutation(s)'][:d_row['mutation(s)'].index('+')] == 'wt*':
                        l_pdbs_compare = ['1L63']
                    else:
                        l_pdbs_compare = ['2LZM']
                        print pdb,'2lzm'
                        stop
                elif protein == 'HEWL':
##                    ## compare to atomic resolution structure if no starting model
##                    l_pdbs_compare = ['2VB1']
                    l_pdbs_compare = [d_HEWL_representative[d_row['SG']]]
            elif n_mutations == 2:
                continue
##                index1 = d_row['mutation(s)'].index('+')+1
##                index2 = index1+d_row['mutation(s)'][index1:].index('+')+1
##                mut1 = d_row['mutation(s)'][0:index1]+d_row['mutation(s)'][index1:index2]
##                mut2 = d_row['mutation(s)'][0:index1]+d_row['mutation(s)'][index2:]
##                l_pdbs_compare = [mut1,mut2,]

        d_comparisons[pdb] = l_pdbs_compare

##    print d_mutations
##    print d_comparisons
##    stop
    return d_mutations, d_comparisons


def parse_citation(l_pdbs,d_mmCIF,):

    d_citation = {}
    for pdb in l_pdbs:
        if 'primary' not in d_mmCIF[pdb[:4]]['_citation.id']:
            print d_mmCIF[pdb[:4]]['_citation.id']
            stop
        for i_citation in range(len(d_mmCIF[pdb[:4]]['_citation.id'])):
            if d_mmCIF[pdb[:4]]['_citation.id'][i_citation] == 'primary':
                break
        id_PubMed = d_mmCIF[pdb[:4]]['_citation.pdbx_database_id_PubMed'][i_citation]
        title = d_mmCIF[pdb[:4]]['_citation.title'][i_citation]
        journal = d_mmCIF[pdb[:4]]['_citation.journal_abbrev'][i_citation]
        if id_PubMed.lower() == 'none':
            stop
        if id_PubMed in ['?','-1',]:
            if journal.lower() == 'to be published':
                d_citation[pdb[:4]] = None
            else:
                d_citation[pdb[:4]] = title
        else:
            d_citation[pdb[:4]] = id_PubMed

    return d_citation


def plot_rmsd_vs_property(
    l_pdbs,d_properties,d_rmsds_overall,mmCIF_item,d_mmCIF,
    prefix,method,
    ):

    l = []
    l_correl = [[],[],]
    l_3d = []

    for i1 in range(len(l_pdbs)-1):

        pdb1 = l_pdbs[i1]

        value1 = core.parse_mmCIF_item(d_mmCIF[pdb1[:4]],mmCIF_item,pdb1,)
        MV1 = core.parse_mmCIF_item(d_mmCIF[pdb1[:4]],'matthews',pdb1,)

        if value1 == None:
            continue
        else:
            value1 = float(value1)

        for i2 in range(i1+1,len(l_pdbs)):

            pdb2 = l_pdbs[i2]

            if bool_spacegroups_identical == True:
                if d_mmCIF[pdb1[:4]]['_symmetry.space_group_name_H-M'] != d_mmCIF[pdb2[:4]]['_symmetry.space_group_name_H-M']:
                    continue
                if protein == 'T4L':
                    if d_mmCIF[pdb1[:4]]['_cell.Z_PDB'] != d_mmCIF[pdb2[:4]]['_cell.Z_PDB']:
                        continue
                    if pdb1[:4] == pdb2[:4]:
                        continue
                    if len( set([pdb1[:4],pdb2[:4],]) & set(['150L','3CDO','1JQU',]) ) == 2:
                        continue

            ## exclude ions
            if len(set([pdb1,pdb2,])&set(['1B2K_A','1LKR_B','1LCN_B',])) == 1:
                continue
            if len(set([pdb1,pdb2,])&set(['1V7T_A','1V7T_B','1XEK_A','2Z18_A','2D4J_A','1XEI_A','1XEJ_A','2Z12_A','2Z19_A','1LMA_A',])) == 1:
                continue

            ## exclude disulfide bonds
            if len( set([pdb1,pdb2,]) & set([
                ## intramolecular
                '152L_A','167L_A','167L_B','172L_A','178L_A',
                '1KNI_A','1L35_A',
                '3GUI_A','3GUJ_A','3GUK_A','3GUK_B','3GUL_A','3GUL_B','3GUM_A','3GUM_B','3GUN_A','3GUN_B','3GUO_A','3GUO_B','3GUP_A','3GUP_B',
                ## intermolecular
                ]) ) > 0:
                continue

            value2 = core.parse_mmCIF_item(d_mmCIF[pdb2[:4]],mmCIF_item,pdb2,)
            MV2 = core.parse_mmCIF_item(d_mmCIF[pdb2[:4]],'matthews',pdb2,)

            if value2 == None:
                continue
            else:
                value2 = float(value2)

##            if prefix == 'pH' and abs(
##                float(core.parse_mmCIF_item(d_mmCIF[pdb1[:4]],'matthews',pdb1,))
##                -
##                float(core.parse_mmCIF_item(d_mmCIF[pdb2[:4]],'matthews',pdb2,))
##                ) > 0.025:
##                continue

            diff = abs(value2-value1)

            rmsd = float(d_rmsds_overall[pdb1][pdb2]['rmsd'])

##            if prefix == 'pH' and rmsd > 0.6:
##                print 'pH', value1,value2,pdb1,pdb2,abs(value2-value1),rmsd

            if prefix == 'MV' and bool_spacegroups_identical == True and diff > 0.3:
                print prefix, 'diff', pdb1,pdb2, diff,rmsd

            if protein == 'T4L' and prefix == 'pH' and rmsd > 1.5:
                print 'pH', pdb1,pdb2,round(diff,1),round(rmsd,1)

            if protein == 'T4L' and prefix == 'MV' and rmsd > 1.5:
                print 'MV', pdb1,pdb2,round(diff,1),round(rmsd,1)

            if protein == 'HEWL' and prefix == 'MV' and rmsd > 1.2:# and diff < 0.1 and diff > 0.005:# and '2FBB_A' not in [pdb1,pdb2,]:
                print 'MV1', pdb1,pdb2,round(diff,1),round(rmsd,1)
                
            if prefix == 'MV' and rmsd < 0.05:
                print 'MV2', round(value1,1),round(value2,1),pdb1,pdb2,abs(value2-value1),round(rmsd,1)

            l += ['%s %s\n' %(diff,rmsd,)]
            l_correl[0] += [diff]
            l_correl[1] += [rmsd]

            l_3d += ['%s %s %s\n' %(diff,abs(MV1-MV2),rmsd,)]

    r = statistics.correlation(l_correl[0],l_correl[1],)
    print '\ncorrelation', prefix, r

    fn_prefix = 'prop_%s_%s_%s_spacegroupsidentical%s_%s_minres%s' %(prefix,protein,method,bool_spacegroups_identical,suffix_exclusion,resolution_min,)
    gnuplot.scatter_plot_2d(
        fn_prefix,l_data=l,
        xlabel='%s diff (Angstrom^3 / Da)' %(prefix),
        title='%s, %s, r = %.2f' %(protein,prefix,r,),
        ylabel='%s RMSD' %(method),
        bool_regression_linear=True,
        bool_remove=True,
        )

    print 'correlation', prefix, r

    if prefix == 'pH':
        fd = open('splot.txt','w')
        fd.writelines(l_3d)
        fd.close()

    if prefix == 'MV':
        stop_complete

    return


def tabulate(d_mmCIF_main,ref_seq,l_wts,):

    print '\nprepare table'

    l_sort = []

    l_pdbs = list(set(d_mmCIF_main.keys()))
    l_pdbs.sort()
    for i in range(len(l_pdbs)):
        pdb = l_pdbs[i]
        d_mmCIF = d_mmCIF_main[pdb]

        ## citation
        for i_citation_author in range(len(d_mmCIF['_citation_author.name'])):
            if d_mmCIF['_citation_author.citation_id'][i_citation_author] == 'primary':
                citation_title = d_mmCIF['_citation.title'][i_citation_author]
                break
        if citation_title[-1] == '.':
            citation_title = citation_title[:-1]

        ## mutations
        l_mutations = []
        if pdb in [s[:4] for s in l_wts]: ## not d_mmCIF['_entity_poly_seq.mon_id'] == ref_seq if 2 CTERM residues missing (3fa0)
            s_mutations = 'wt'
        else:
            for i in range(len(d_mmCIF['_entity_poly_seq.mon_id'])):
                res_id_mmCIF = d_mmCIF['_entity_poly_seq.mon_id'][i]
                res_id_ref = ref_seq[i]
                if res_id_mmCIF != res_id_ref:
                    ## selenomethionine substitution                    
                    if res_id_mmCIF == 'MSE' and res_id_ref == 'MET':
                        continue ## 1D3N,1C6X
                    res_no = i+1
                    if res_id_mmCIF in d_321.keys():
                        res_id_mmCIF = d_321[res_id_mmCIF]
                    elif res_id_mmCIF == 'MSE':
                        res_id_mmCIF = 'M'
                    else:
                        res_id_mmCIF = 'X'
                    l_mutations += ['%1s%i%1s' %(d_321[res_id_ref],res_no,res_id_mmCIF,)]
            ## error/correction
            if pdb == '1PQK':
                l_mutations.remove('G77A')
            s_mutations = '+'.join(l_mutations)

        ## ligands
        l_ligands = []
        for i_chem_comp in range(len(d_mmCIF['_chem_comp.id'])):
            chem_comp_id = d_mmCIF['_chem_comp.id'][i_chem_comp]
            chem_comp_type = d_mmCIF['_chem_comp.type'][i_chem_comp]
            if chem_comp_id == 'HOH':
                continue
            if chem_comp_type == 'L-peptide linking':
                continue
            if chem_comp_id == 'GLY' and chem_comp_type == 'peptide linking':
                continue
            l_ligands += [chem_comp_id]
        s_ligands = ','.join(list(set(l_ligands)-set([
            ## ligands to be ignored
            'BME','HED',#'CME', (CME is a MODRES, not a ligand) (HED is oxidized BME, necessary for fast crystallization cf. 1l90)
            'NA','CL','K','PO4','SO4','ACT','CO3','CA','MG','NO3',
            'AZI', ## protease inhibitor?
            'MPD','MRD', ## ("An overview on 2-methyl-2,4-pentanediol in crystallization and in crystals of biological macromolecules")
##            'AR','KR','XE',
            'EPE', ## HEPES buffer
            'TRS', ## Tris buffer
            'GOL', ## Glycerol
            'HEZ', ## BME replacement (1quh)
            ])))
                    
        spacegroup = core.parse_mmCIF_item(d_mmCIF,'_symmetry.space_group_name_H-M',pdb,)
        T = core.parse_mmCIF_item(d_mmCIF,'_diffrn.ambient_temp',pdb,)
        pH = core.parse_mmCIF_item(d_mmCIF,'_exptl_crystal_grow.pH',pdb,)
        resolution = core.parse_mmCIF_item(d_mmCIF,'_refine.ls_d_res_high',pdb,)
        starting_model = core.parse_mmCIF_item(d_mmCIF,'_refine.pdbx_starting_model',pdb,)
        method_to_determine_struct = core.parse_mmCIF_item(d_mmCIF,'_refine.pdbx_method_to_determine_struct',pdb,)

        if method_to_determine_struct == None:
            if starting_model == None:
                method_to_determine_struct = 'N/A'
            else:
                method_to_determine_struct = 'MR'
        method_to_determine_struct = method_to_determine_struct.replace('MOLECULAR REPLACEMENT','MR',)
        method_to_determine_struct = method_to_determine_struct.replace('MOLECULAR SUBSTITUTION','MR',)
        method_to_determine_struct = method_to_determine_struct.replace('STRUCTURE KNOWN','MR',)
        method_to_determine_struct = method_to_determine_struct.replace('Isomorphous replacement','MR',) ## 1llh
        method_to_determine_struct = method_to_determine_struct.replace('Rigid body','MR',) ## 1kni
        if protein == 'HEWL':
            method_to_determine_struct = method_to_determine_struct.replace('DIRECT REFINEMENT','MR',)
            method_to_determine_struct = method_to_determine_struct.replace('SHELX-97','MR',)
        if protein == 'T4L':
            method_to_determine_struct = method_to_determine_struct.replace('REFMAC','MR',)
        method_to_determine_struct = method_to_determine_struct.replace('FOURIER SYNTHESIS','FS',)
        method_to_determine_struct = method_to_determine_struct.replace('AB INITIO PHASING','AI',)
        method_to_determine_struct = method_to_determine_struct.replace('AB INITIO','AI',)
        method_to_determine_struct = method_to_determine_struct.replace('DIFFERENCE FOURIER','DF',)
        if starting_model != None:
            method_to_determine_struct = 'MR'
##        if method_to_determine_struct == 'MR' and starting_model == None and pdb not in ['1B6I','1LGU',]:
##            print pdb
##            stop

        matthews = core.parse_mmCIF_item(d_mmCIF,'matthews',pdb,)
        a = float(d_mmCIF_main[pdb]['_cell.length_a'][0])
        b = float(d_mmCIF_main[pdb]['_cell.length_b'][0])
        c = float(d_mmCIF_main[pdb]['_cell.length_c'][0])
        alpha = float(d_mmCIF_main[pdb]['_cell.angle_alpha'][0])
        beta = float(d_mmCIF_main[pdb]['_cell.angle_beta'][0])
        gamma = float(d_mmCIF_main[pdb]['_cell.angle_gamma'][0])
        Z = int(d_mmCIF_main[pdb]['_cell.Z_PDB'][0])

        R = d_mmCIF_main[pdb]['_refine.ls_R_factor_obs'][0]
        R_free = d_mmCIF_main[pdb]['_refine.ls_R_factor_R_free'][0]

        pubmed = d_mmCIF_main[pdb]['_citation.pdbx_database_id_PubMed'][0]

        ##
        ## title and descriptor
        ##
        struct_title = ''.join(d_mmCIF_main[pdb]['_struct.title']).upper()
        descriptor = ''.join(d_mmCIF_main[pdb]['_struct.pdbx_descriptor']).upper()
        
        title,descriptor = title_and_descriptor(
            d_mmCIF_main,pdb,struct_title,descriptor,
            citation_title,l_ligands,l_mutations,
            spacegroup,T,pH,resolution,
            )

        ##
        ## mutation
        ##
        if protein == 'T4L':
            if len(l_mutations) == 0:
                s_mutations = 'wt'
                mutation_class = 'wt'
            else:
                title, s_mutations, mutation_class = organize_mutants(l_mutations,title,s_mutations,)
        elif protein == 'HEWL':
            mutation_class = ''

        ## formatting
        if T != None:
            T = '%3i' %(int(float(T)))
        else:
            T = 'N/A'
        if pH != None:
            pH = '%3.1f' %(float(pH))
        else:
            pH = 'N/A'
        if starting_model == None:
            starting_model = 'N/A'
        if R != '?':
            R = '%4.2f' %(float(R))
        if R_free != '?':
            R_free = ' %4.2f' %(float(R_free))

        ## manual sorting
        if spacegroup == 'P 32 2 1':
            spacegroup_sort = '1'+spacegroup
        else:
            spacegroup_sort = spacegroup

        ## sort wt and wt* independently
        s_mutations_sort = s_mutations.replace('*','')[3:]

        ## sort by absence/presence of ligands
        if s_ligands == '':
            bool_ligands = 0
        else:
            bool_ligands = 1

        ## concatenation
        comment = title+descriptor

        if len(comment) > 9:
            print 'long comment', pdb, comment

        if 'MSE' in d_mmCIF['_entity_poly_seq.mon_id']:
            comment += 'SELENO'

        if protein == 'T4L':
            s_app = '%-9s\t%s' %(comment, s_mutations,)
            l_header_app = ['comment','mutation(s)',]
        elif protein == 'HEWL':
            s_app = '%s\t%s' %(s_mutations, comment,)
            l_header_app = ['mutation(s)','comment',]
        else:
            print protein
            stop

        l_sort += [
            [
                spacegroup_sort,
                mutation_class,
                pubmed,
                pdb,
                s_mutations_sort, bool_ligands,
                '%-4s\t%3s\t%3s\t%3.1f\t%4s\t%6s\t%-10s\t%-5s\t%-4s\t%3.1f\t%2i\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%-3s\t%s\n' %(
                    pdb,T,pH,float(resolution),R,R_free,spacegroup,
                    method_to_determine_struct,starting_model,
                    matthews,Z,a,b,c,alpha,beta,gamma,
                    s_ligands,
                    s_app,
                    ),
                ]
            ]
    
    ## sort lines
    l_sort.sort()

    ## add column headers
    lines = [
        '%-4s\t%3s\t%3s\tres\tR\tRfree\t%-10s\t%-5s\tSM  \t MV\t Z\t    a\t    b\t    c\talpha\t beta\tgamma\t%-3s\t%s\t%s\n' %tuple([
            'ID', 'T', 'pH', 'SG',
            'method',
            'lig',
            ]+l_header_app),
        ]
    ## loop over sorted lines
    for l in l_sort:
        ## append line
        lines += [l[-1]]

    ## write table
    fd = open('table_%s.txt' %(protein),'w')
    fd.writelines(lines)
    fd.close()

##    stop_finished

    return


def organize_mutants(l_mutations,title,s_mutations,):

    l_set_mutations = []
    l_mutation_IDs = []
    l_mutation_classes = []
    for set_mutations,mutation_ID,mutation_class in [
##                [set([       'I3L',       'S38D','A41V','A82P',                       'N116D',        'V131A',        'N144D',]),'wt+stab2',], ## only 189l
##        [set(['C97X',      'T21C','S38D',              'L99A','M102E','E108V',        'S117V',        'T142C','N144D',]),'wt+L99A+M102E+stab1+disulf','M102E',], ## 3gui-3gup
        [set(['C54T','C97A','L32T','T34K','K35V','S36D','P37G','S38N','L39S',]),'wt*+X32-39X','X32-39X',],
        [set(['I3C','I9C','T21C','C54T','T142C','L164C',]),'wt+3xdisulf','disulf',],
        [set(['I3C','I9C',       'C54T',        'L164C',]),'wt+2xdisulf','disulf',],
##        [set(['K16E','R119E','K135E','K147E',]),'wt+KRKK>EEEE','charge?',], ## QUAD
        ## wt* CORE 7, CORE 10
        [set(['C54T','C97A',       'I78V','V87M',                                                      'L118I','M120Y',        'L133F','V149I','T152V',]),'wt*+CORE7','core',],
        [set(['C54T','C97A',              'V87I',              'I100V','M102L','V103I','M106I','V111A',        'M120Y',        'L133F','V149I','T152V',]),'wt*+CORE10','core',],

        ## wt* + L99A + M102Q
        [set(['C54T','C97A','L99A','M102Q',]),'wt*+L99A+M102Q','cavity',],
        ## wt + L99A + M102Q
        [set([              'L99A','M102Q',]),'wt+L99A+M102Q','cavity',],
        ## wt* + L99A
        [set(['C54T','C97A','L99A',]),'wt*+L99A','cavity',],
        ## wt*
        [set(['C54T','C97A',]),'wt*','',],
        ]:
        l_set_mutations += [set_mutations]
        l_mutation_IDs += [mutation_ID]
        l_mutation_classes += [mutation_class]

    ##
    ## add prefix
    ##
    bool_break = False
    for i_mutation in range(len(l_set_mutations)):
        set_mutations = l_set_mutations[i_mutation]                
        mutation_ID = l_mutation_IDs[i_mutation]
        mutation_class = l_mutation_classes[i_mutation]
        ## all mutations within set present in pdb (set_mutations subset of l_mutations)
        if len( set_mutations - set(l_mutations) ) == 0:
            
            s_mutations = mutation_ID

##                    if len( set(l_mutations)-set_mutations ) > 0:
##                        s_mutations += '+'
            for mutation in l_mutations:
                if mutation in set_mutations:
                    continue
                s_mutations += '+%s' %(mutation)
            bool_break = True
            break

        ## one extra mutation
        elif len( set_mutations - set(l_mutations) ) == 1 and len( set(l_mutations) & set_mutations ) > 4:
            s_mutations_reverse = '-'.join(list(set_mutations-set(l_mutations)))
            s_mutations_extra = '+'.join(list(set(l_mutations)-set_mutations))
            if len(s_mutations_reverse) > 5:
                print s_mutations_reverse
                stop
            s_mutations = mutation_ID
            if len(s_mutations_extra) > 0:
                s_mutations += '+'+s_mutations_extra
            if len(s_mutations_reverse) > 0:
                s_mutations += '-'+s_mutations_reverse
                title = title.replace('%1s%i%1s' %(s_mutations_reverse[-1],int(s_mutations_reverse[1:-1]),s_mutations_reverse[0],),'')
            bool_break = True
            break

    if bool_break == False:
        s_mutations = 'wt+'+s_mutations

    ##
    ## group mutation
    ##
    if len( set(['L84M','L91M','L121M','L133M','V87M','V103M','V111M',]) & set(l_mutations) ) > 0:
        mutation_class = 'core'
##        if len( set(['L84M','L91M','L121M','L133M','V87M','V103M','V111M',]) & set(l_mutations) ) > 0:
##            res2 = 'M'
##        elif len( set(['L84X','L91X','L121X','L133X','V87X','V103X','V111X',]) & set(l_mutations) ) > 0:
##            res2 = 'X'
##        else:
##            print ':::::::', pdb, l_mutations
##            stop
##        for res1,set_mutations in [
##            ['V',set(['V87M','V103M','V111M','V87X','V103X','V111X',]),],
##            ['L',set(['L84M','L84X','L91M','L91X','L121M','L121X','L133M','L133X',]),],
##            ]:
##            if len( set_mutations & set(l_mutations) ) >= 2:
##                mutation_class = 'core'
##                short = res1
##                for res_no in [
##                    ## LEU --> NET
##                    84,91,99,118,121,133,
##                    ## VAL --> MET
##                    87,103,111,
##                    ]:
##                    if '%s%i%s' %(res1,res_no,res2,) in l_mutations:
##                        short += '%i,' %(res_no)
##                        s_mutations = s_mutations.replace('+%s%i%s' %(res1,res_no,res2,),'')
##                        s_mutations = s_mutations.replace('%s%i%s' %(res1,res_no,res2,),'')
##                short = '+'+short[:-1]+res2
##                s_mutations = s_mutations[:3]+short+s_mutations[3:]

    if mutation_class == '':
        if len( set([
            'L99A','M102Q',
            'M106A',
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'cavity'
        elif len( set(['A41X','S44X','T115X','L118X','V131X','T151X']) & set(l_mutations) ) > 0:
            mutation_class = 'spin'
        elif len( set([
            'I78A',
            'L99F','L99I','L99M','L99V',
            'M102L',
            'V111L',
            'L121A','L133F',
            'A129L','A129M','A129F','A129W',
            'F153L','F153M','F153I','F153A','F153V','F153W',
            ]) & set(l_mutations) ) > 0: ## L121A,M129L is size switch...
            mutation_class = 'core'
        elif len( set([
            'L99G',
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'cavity'
        elif len( set(['M6A','I17A','I27A','I29A','I50A','I58A','F67A','L84A','V87A','I100A','V103A','V111A','L121A','L133A','V149A',]) & set(l_mutations) ) > 0:
            mutation_class = 'L2S'
        elif len( set([
            'R96A','R96C','R96D','R96E','R96F','R96G','R96H','R96I','R96K','R96L','R96M','R96N','R96P','R96Q','R96S','R96T','R96V','R96W','R96Y',
            'D72A','K85A','D89A', ## R96H salt bridge
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'R96' ## R96H destabilization
        elif len( set([
            'T157A','T157C','T157D','T157E','T157F','T157G','T157H','T157I','T157L','T157N','T157R','T157S','T157V',
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'T157' ## T157 hydrogen bond stabilizations
        elif len( set(['P86A','P86C','P86D','P86G','P86H','P86L','P86R','P86S',]) & set(l_mutations) ) > 0:
            mutation_class = 'P86' ## P86 helix extension
        elif len( set(['I3P','M6I',]) & set(l_mutations) ) > 0:
            mutation_class = 'hinge' ## hinge bending
        elif len( set(['N55G','K124G',]) & set(l_mutations) ) > 0:
            mutation_class = '"Left-handed helical" residue' ## ## P86 helix extension
        elif len( set([
            ## helix propensity
            'S44E','S44G','S44I','S44K','S44L','S44N','S44P','S44R','S44T','S44V',
            'V131D','V131E','V131G','V131I','V131L','V131M','V131S',
            ## helix capping
            'T59A','T59D','T59G','T59N','T59S','T59V',
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'helix' ## helix propensity/capping
        elif len( set(['K16E','R119E','K135E','K147E','R154E',]) & set(l_mutations) ) > 0:
            mutation_class = 'RevCha'
        elif len( set([
            'A41S','A42S','A49S','A73S','A82S','A93S','A98S','A130S','A134S',
            'V75T','V87T','V149T',
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'burOH'
        elif len( set(['Q105E','Q105G','Q105A',]) & set(l_mutations) ) > 0:
            mutation_class = 'Trp138'
        elif len( set(['K60P','A82P','G77A','G113A',]) & set(l_mutations) ) > 0:
            mutation_class = 'stab_entropic' ## entropic stabilization
        elif len( set(['G156D',]) & set(l_mutations) ) > 0:
            mutation_class = 'stab' ## random stabilization
        elif len( set(['S38D','T109D','T115E','N116D','N144D',]) & set(l_mutations) ) > 0:
            mutation_class = 'stab_helixdipole' ## helix dipole stabilization
        elif len( set(['T115E','Q123E',]) & set(l_mutations) ) > 0:
            mutation_class = 'saltbridge' ## salt bridge
        elif len( set(['P143A',]) & set(l_mutations) ) > 0:
            mutation_class = 'pro'
        elif len( set(['K124D',]) & set(l_mutations) ) > 0:
            mutation_class = 'unpub'
        elif len( set(['L133G',]) & set(l_mutations) ) > 0:
            mutation_class = 'cavity'
        elif len( set(['A98V','T152S','V149I','V149C',]) & set(l_mutations) ) > 0:
            mutation_class = 'helix' ## helix packing analysis
        elif len( set(['I3Y','I3V','I3L',]) & set(l_mutations) ) > 0:
            mutation_class = 'hydrophobic' ## hydrophobic replacement
        elif len( set(['I3C','V131C',]) & set(l_mutations) ) > 0:
            mutation_class = 'disulf'
        ## mutation class not automatically assigned and more than one mutation
        elif len( set([
            ## 14-22 (sheet excl Leu15, Tyr18, Asp20)
            'R14A','K16A','I17A','K19A','T21A','E22A',
            ## 24-27 (sheet between Gly23 and Gly28)
            'Y24A','Y25A','T26A','I27A',
            ## helix115-123
            'T115A','N116A','S117A','R119A','M120A','Q122A','Q123A',
            ## 34-37 (loop)
            'T34A','K35A','S36A','P37A',
            ## 40-48 (helix)
            'N40A',       'S44A','E45A',       'D47A','K48A',
            ## 53-57 (turn, excl Gly52/Gly56)
            'N53A','N55A','V57A',
            ## 127-141
            'D127A','E128A','V131A','N132A',        'K135A','S136A','R137A','Y139A','N140A','Q141A',
            ]) & set(l_mutations) ) > 0:
            mutation_class = 'Ala'
        elif bool_break == False and '+' in s_mutations:
            print '@@@@@@', mutation_class, s_mutations
            mutation_class = 'wt'

    return title, s_mutations, mutation_class


def title_and_descriptor(
d_mmCIF_main,pdb,title,descriptor,
citation_title,l_ligands,l_mutations,
spacegroup,T,pH,resolution,
):

    d_mmCIF = d_mmCIF_main[pdb]

    title = title.replace(citation_title.upper(),'')
    title = title.replace(citation_title.upper().replace('----',' (RIGHT ARROW) '),'')
    title = title.replace(citation_title.upper().replace('CROSS-LINK','CROSSLINK'),'')
    title = title.replace(citation_title.upper().replace('-->','-> '),'')
    title = title.replace(citation_title.upper().replace('A RESOLUTION','ANGSTROMS RESOLUTION'),'')

    ## shorter
    title = title.replace('COMPARISON OF RADIATION-INDUCED DECAY AND STRUCTURE REFINEMENT FROM X-RAY DATA COLLECTED FROM LYSOZYME CRYSTALS AT LOW AND AMBIENT TEMPERATURES','RADIATION INDUCED DECAY')
    title = title.replace('SPECIFIC CHEMICAL AND STRUCTURAL DAMAGE CAUSED BY INTENSE SYNCHROTON RADIATION TO HEWL','RADIATION INDUCED DECAY') ## 1QIO

    ## structure title = shared citation title
    title = title.replace('RAPID CRYSTALLIZATION OF T4 LYSOZYME BY INTERMOLECULAR DISULFIDE CROSSLINKING','')

    ## mutants
    title = title.replace('WITH THE UNNATURAL AMINO ACID P-ACETYL-L-PHENYLALANINE INCORPORATED AT POSITION 131','')
    title = title.replace('HEN EGG WHITE LYSOZYME MUTANT WITH ALANINE SUBSTITUTED FOR GLYCINE','')
    title = title.replace('STRUCTURAL AND THERMODYNAMIC ANALYSIS OF COMPENSATING MUTATIONS WITHIN THE CORE OF CHICKEN EGG WHITE LYSOZYME','')
    title = title.replace('DISSECTION OF PROTEIN-CARBOHYDRATE INTERACTIONS IN MUTANT HEN EGG-WHITE LYSOZYME COMPLEXES AND THEIR HYDROLYTIC ACTIVITY','')
    title = title.replace('THERMAL STABILITY DETERMINANTS OF CHICKEN EGG-WHITE LYSOZYME CORE MUTANTS: HYDROPHOBICITY, PACKING VOLUME AND CONSERVED BURIED WATER MOLECULES','')
    title = title.replace('MULTIPLE METHIONINE SUBSTITUTIONS IN T4 LYSOZYME','')
    title = title.replace('METHIONINE CORE MUTATION','')
    title = title.replace('METHIONINE CORE MUTANT OF T4 LYSOZYME','')
    title = title.replace('STABILIZING DISULFIDE BRIDGE MUTANT OF T4 LYSOZYME','')
    title = title.replace('POLAR AND NON-POLAR CAVITIES IN PHAGE T4 LYSOZYME','')
    title = title.replace('T4 LYSOZYME METHIONINE CORE MUTANT','')
    title = title.replace('BACTERIOPHAGE LYSOZYME T4 LYSOZYME MUTANT','') ## 3c82
    title = title.replace('IM MUTANT OF LYSOZYME','')
    title = title.replace('PROTEIN STRUCTURE PLASTICITY EXEMPLIFIED BY INSERTION AND DELETION MUTANTS IN T4 LYSOZYME','')
    title = title.replace('STABILIZATION OF HEN EGG WHITE LYSOZYME BY A CAVITY-FILLING MUTATION','')
    title = title.replace('ANALYSIS OF THE STABILIZATION OF HEN LYSOZYME WITH THE HELIX DIPOLE AND CHARGED SIDE CHAINS','')
    title = title.replace('CONTRIBUTIONS OF ALL 20 AMINO ACIDS AT SITE 96 TO STABILITY AND STRUCTURE OF T4 LYSOZYME','')
    title = title.replace('ENHANCED PROTEIN THERMOSTABILITY FROM SITE-DIRECTED MUTATIONS THAT DECREASE THE ENTROPY OF UNFOLDING','')
    title = title.replace('ARE CARBOXY TERMINII OF HELICES CODED BY THE LOCAL SEQUENCE OR BY TERTIARY STRUCTURE CONTACTS','')
    title = title.replace('THE ENERGETIC COST AND THE STRUCTURAL CONSEQUENCES OF BURYING A HYDROXYL GROUP WITHIN THE CORE OF A PROTEIN DETERMINED FROM ALA TO SER AND VAL TO THR SUBSTITUTIONS IN T4 LYSOZYME','')
    title = title.replace('ROLE OF BACKBONE FLEXIBILITY IN THE ACCOMMODATION OF VARIANTS THAT REPACK THE CORE OF T4 LYSOZYME','')
    title = title.replace('STRUCTURAL BASIS OF ALPHA-HELIX PROPENSITY AT TWO SITES IN T4 LYSOZYME','')
    title = title.replace('STRUCTURAL BASIS OF AMINO ACID ALPHA HELIX PROPENSITY','')
    title = title.replace('STRUCTURES OF RANDOMLY GENERATED MUTANTS OF T4 LYSOZYME SHOW THAT PROTEIN STABILITY CAN BE ENHANCED BY RELAXATION OF STRAIN AND BY IMPROVED HYDROGEN BONDING VIA BOUND SOLVENT','')
    title = title.replace('CONTRIBUTIONS OF HYDROGEN BONDS OF THR 157 TO THE THERMODYNAMIC STABILITY OF PHAGE T4 LYSOZYME','')
    title = title.replace('NEW AZABORINE COMPOUNDS BIND TO THE T4 LYSOZYME L99A CAVITY','')
    title = title.replace('REPLACEMENTS OF PRO86 IN PHAGE T4 LYSOZYME EXTEND AN ALPHA-HELIX BUT DO NOT ALTER PROTEIN STABILITY','')
    title = title.replace('CUMULATIVE SITE-DIRECTED CHARGE-CHANGE REPLACEMENTS IN BACTERIOPHAGE T4 LYSOZYME SUGGEST THAT LONG-RANGE ELECTROSTATIC INTERACTIONS CONTRIBUTE LITTLE TO PROTEIN STABILITY','')
    title = title.replace('ANALYSIS OF THE INTERACTION BETWEEN CHARGED SIDE CHAINS AND THE ALPHA-HELIX DIPOLE USING DESIGNED THERMOSTABLE MUTANTS OF PHAGE T4 LYSOZYME','')
    title = title.replace('SIMILAR HYDROPHOBIC REPLACEMENTS OF LEU 99 AND PHE 153 WITHIN THE CORE OF T4 LYSOZYME HAVE DIFFERENT STRUCTURAL AND THERMODYNAMIC CONSEQUENCES','')
    title = title.replace('DISSECTION OF HELIX CAPPING IN T4 LYSOZYME BY STRUCTURAL AND THERMODYNAMIC ANALYSIS OF SIX AMINO ACID SUBSTITUTIONS AT THR 59','')
    title = title.replace('GENERATING LIGAND BINDING SITES IN T4 LYSOZYME USING DEFICIENCY-CREATING SUBSTITUTIONS','')
    title = title.replace('MULTIPLE STABILIZING ALANINE REPLACEMENTS WITHIN ALPHA-HELIX 126-134 OF T4 LYSOZYME HAVE INDEPENDENT, ADDITIVE EFFECTS ON BOTH STRUCTURE AND STABILITY','')
    title = title.replace('ENHANCEMENT OF PROTEIN STABILITY BY THE COMBINATION OF POINT MUTATIONS IN T4 LYSOZYME IS ADDITIVE','')
    title = title.replace('THE RESPONSE OF T4 LYSOZYME TO LARGE-TO-SMALL SUBSTITUTIONS WITHIN THE CORE AND ITS RELATION TO THE HYDROPHOBIC EFFECT','')
    title = title.replace('CONTRIBUTION OF ALL 20 AMINO ACIDS AT SITE 96 TO THE STABILITY AND STRUCTURE OF T4 LYSOZYME','') ## 3C8Q
    title = title.replace('THERMODYNAMIC AND STRUCTURAL COMPENSATION IN "SIZE-SWITCH" CORE-REPACKING VARIANTS OF T4 LYSOZYME','')
    title = title.replace('N-TERMINAL DOMAIN CORE METHIONINE MUTATION','')
    title = title.replace('CRYSTAL STRUCTURE T4 LYSOZYME INCORPORATING AN UNNATURAL AMINO ACID P-IODO-L-PHENYLALANINE AT POSITION 153','')
    title = title.replace('CONTRIBUTIONS OF ENGINEERED SURFACE SALT BRIDGES TO THE STABILITY OF T4 LYSOZYME DETERMINED BY DIRECTED MUTAGENESIS','')
    title = title.replace('STRUCTURAL AND THERMODYNAMIC ANALYSIS OF THE PACKING OF TWO ALPHA-HELICES IN BACTERIOPHAGE T4 LYSOZYME','')
    title = title.replace('DESIGN AND STRUCTURAL ANALYSIS OF ALTERNATIVE HYDROPHOBIC CORE PACKING ARRANGEMENTS IN BACTERIOPHAGE T4 LYSOZYME','')
    title = title.replace('TOLERANCE OF T4 LYSOZYME TO MULTIPLE XAA (RIGHT ARROW) ALA SUBSTITUTIONS: A POLYALANINE ALPHA-HELIX CONTAINING TEN CONSECUTIVE ALANINES','')
    title = title.replace('HYDROPHOBIC CORE REPACKING AND AROMATIC-AROMATIC INTERACTION IN THE THERMOSTABLE MUTANT OF T4 LYSOZYME SER 117 (RIGHT ARROW) PHE','')
    title = title.replace('DETERMINATION OF ALPHA-HELIX PROPENSITY WITHIN THE CONTEXT OF A FOLDED PROTEIN: SITES 44 AND 131 IN BACTERIOPHAGE T4 LYSOZYME','')
    title = title.replace('EVALUATION AT ATOMIC RESOLUTION OF THE ROLE OF STRAIN IN DESTABILIZING THE TEMPERATURE-SENSITIVE T4 LYSOZYME MUTANT ARG 96 --> HIS','')
    title = title.replace('EVAULAUTION AT ATOMIC RESOLUTION OF THE ROLE OF STRAIN IN DESTABILIZING THE TEMPERATURE SENSITIVE T4 LYSOZYME MUTANT ARG96-->HIS','') ## 3f8v
    title = title.replace('THE INTRODUCTION OF STRAIN AND ITS EFFECTS ON THE STRUCTURE AND STABILITY OF T4 LYSOZYME','')
    title = title.replace('CONTRIBUTIONS OF LEFT-HANDED HELICAL RESIDUES TO THE STRUCTURE AND STABILITY OF BACTERIOPHAGE T4 LYSOZYME','')
    title = title.replace('HYDROPHOBIC STABILIZATION IN  DETERMINED DIRECTLY BY MULTIPLE SUBSTITUTIONS OF ILE 3','')
    title = title.replace('T4 LYSOZYME M102E/L99A MUTANT WITH BURIED CHARGE IN APOLAR CAVITY','')
    title = title.replace('REPACKING OF THE CORE OF T4 LYSOZYME BY AUTOMATED DESIGN','')
    title = title.replace('T4 LYSOZYME CORE REPACKING MUTANT','')
    title = title.replace('T4 LYOSZYME CORE REPACKING MUTANT','')
    title = title.replace('AN ADAPTABLE METAL-BINDING SITE ENGINEERED INTO T4 LYSOZYME','')
    title = title.replace('T4 LYSOZYME SUBSTITUTED WITH SELENOMETHIONINE','') ## 1CX6
    title = title.replace('T4 LYSOZYME CORE REPACKING BACK-REVERTANT','')
    title = title.replace('CORE REDESIGN BACK-REVERTANT','')
    if 'CORE' in title:
        title = title.replace('/TA','')
    title = title.replace('/CORE7','')
    title = title.replace('/CORE10','')
    title = title.replace('CORE7','')
    title = title.replace('CORE10','')
    descriptor = descriptor.replace('LYSOZYME (E.C.3.2.1.17) MUTANT WITH THH 34, LYS 35, SER 36 AND PRO 37 REPLACED BY ALANINE','') ## 151l
    descriptor = descriptor.replace('S-GAMMA97-BETA-MERCAPTOETHANOL LYSOZYME (E.C.3.2.1.17)','')
    descriptor = descriptor.replace('SGAMMA86-BETA-MERCAPTOETHANOL-LYSOZYME (E.C.3.2.1.17)','') ## 1l26
    descriptor = descriptor.replace('S-GAMMA97-BETA-MERCAPTOETHANOL-LYSOZYME (E.C.3.2.1.17)','')
    descriptor = descriptor.replace('S-GAMMA157-BETA-MERCAPTOETHANOL-LYSOZYME (E.C.3.2.1.17)','')
    descriptor = descriptor.replace('(THREE ALANINE SUBSTITUTIONS)','')

    ## irrelevant or shared title
    if pdb[:4] in [
        ## SULFUR SAD PHASING/EXPERIMENTS
        '2W1L','2W1M','2W1X','2W1Y','3EXD',
        ## UV LASER EXCITED FLUORESCENCE
        '2C8O','2C8P',
        ## HEAVY WATER SOLUTION
        '2D4I','2D4J',
        ## LOW HUMIDITY
        '1LMA','2Z12','2Z18','2Z19',
        ## LYSOZYME GROWN AT BASIC PH AND ITS LOW HUMIDITY VARIANT
        '1HSX','1HSW',
        ## MAGNETIC FIELD
        ## ZERO GRAVITY
        ]:
        title = ''
        

    title = title.replace('EFFECT OF ALCOHOLS ON PROTEIN HYDRATION','') ## 1ykx...
    title = title.replace('HEW LYSOZYME: TRP...NA CATION-PI INTERACTION','') ## 1lpi
    title = title.replace('STRUCTURE OF TETRAGONAL HEN EGG WHITE LYSOZYME AT 0.94 A FROM CRYSTALS GROWN BY THE COUNTER-DIFFUSION METHOD','') ## 1iee
    title = title.replace('NOVEL BROMATE SPECIES TRAPPED WITHIN A PROTEIN CRYSTAL','') ## 2d6b
    title = title.replace('LYSOZYME STRUCTURE DERIVED FROM THIN-FILM-BASED CRYSTALS','')
    title = title.replace('X-RAY DIFFRACTION STUDIES OF ADDUCTS BETWEEN ANTICANCER PLATINUM DRUGS AND HEN EGG WHITE LYSOZYME','')
    title = title.replace('STUDIES OF MONOCLINIC HEN EGG WHITE LYSOZYME. IV. X-RAY REFINEMENT AT 1.8 ANGSTROM RESOLUTION AND A COMPARISON OF THE VARIABLE REGIONS IN THE POLYMORPHIC FORMS','')
    title = title.replace('PARENT STRUCTURE OF HEN EGG WHITE LYSOZYME GROWN IN ACIDIC PH 4.8. REFINEMENT FOR COMPARISON WITH CROSSLINKED MOLECULES OF LYSOZYME','')
    title = title.replace('HIGH ENERGY TETRAGONAL LYSOZYME X-RAY STRUCTURE','')
    title = title.replace('STRUCTURE OF HEN EGG-WHITE LYSOZYME DETERMINED FROM CRYSTALS GROWN IN PH 7.5','') ## 2hub
    title = title.replace('ORTHORHOMBIC LYSOZYME CRYSTALLIZED AT HIGH TEMPERATURE (310K)','') ## 1bgi
    title = title.replace('ANOMALOUS SIGNAL OF SOLVENT BROMINES USED FOR PHASING OF LYSOZYME','') ## 1lz9
    title = title.replace('LYSOZYME PHASED ON ANOMALOUS SIGNAL OF SULFURS AND CHLORINES','') ## 1lz8
    title = title.replace('THE KEDGE HOLMIUM DERIVATIVE OF HEN EGG-WHITE LYSOZYME AT HIGH RESOLUTION FROM SINGLE WAVELENGTH ANOMALOUS DIFFRACTION','') ## 2bpu
    title = title.replace('SIRAS STRUCTURE OF TETRAGONAL LYSOSYME USING DERIVATIVE DATA COLLECTED AT THE HIGH ENERGY REMOTE HOLMIUM KEDGE','') ## 2cgi
    title = title.replace('X-RAY STRUCTURE OF A MONOCLINIC FORM OF HEN EGG-WHITE LYSOZYME CRYSTALLIZED AT 313K. COMPARISON OF TWO INDEPENDENT MOLECULES','') ## 1lys
    title = title.replace('TRICLINIC LYSOZYME WITH LOW SOLVENT CONTENT OBTAINED BY PHASE TRANSITION','')
    title = title.replace('TRICARBONYLMANGANESE(I)-LYSOZYME COMPLEX : A STRUCTURALLY CHARACTERIZED ORGANOMETALLIC PROTEIN','')
    title = title.replace('THE 1.40 A STRUCTURE OF SPACEHAB-01 HEN EGG WHITE LYSOZYME','') ## 194l
    title = title.replace('THE STRUCTURE OF HEW LYSOZYME ORTHORHOMBIC CRYSTAL GROWTH UNDER A HIGH MAGNETIC FIELD','')
    title = title.replace("X-RAY STRUCTURE OF HEW LYSOZYME ORTHORHOMBIC CRYSTAL FORMED IN THE EARTH'S MAGNETIC FIELD",'')
    title = title.replace('REFINEMENT OF TRICLINIC LYSOZYME. II. THE METHOD OF STEREOCHEMICALLY RESTRAINED LEAST-SQUARES','') ## 2lzt
    title = title.replace('HEN EGG LYSOZYME CROSS-LINKED BY GLUTARALDEHYDE','') ## 2f2n not crosslinked
    title = title.replace('ON THE ROUTINE USE OF SOFT X-RAYS IN MACROMOLECULAR CRYSTALLOGRAPHY, PART III- THE OPTIMAL DATA COLLECTION WAVELENGTH','')
    title = title.replace('THE INFLUENCE OF TEMPERATURE ON LYSOZYME CRYSTALS. STRUCTURE AND DYNAMICS OF PROTEIN AND WATER','')
    title = title.replace('STUDIES OF HEWL. IV. X-RAY REFINEMENT AT 1.8 ANGSTROM RESOLUTION AND A COMPARISON OF THE VARIABLE REGIONS IN THE POLYMORPHIC FORMS','')
    title = title.replace('USE OF AN ION-BINDING SITE TO BYPASS THE 1000-ATOM LIMIT TO AB INITIO STRUCTURE DETERMINATION BY DIRECT METHODS','')
    title = title.replace('USE OF A HALIDE BINDING SITE TO BYPASS THE 1000-ATOM LIMIT TO AB INITIO STRUCTURE DETERMINATION','') ## 1swy
    title = title.replace('USE OF AN ION-BINDING SITE TO BYPASS THE 1000-ATOM LIMIT TO STRUCTURE DETERMINATION BY DIRECT METHODS','')
    title = title.replace('USE OF A HALIDE BINDING SITE TO BYPASS THE 1000-ATOM LIMIT TO STRUCTURE DETERMINATION BY DIRECT METHODS','') ## 
    title = title.replace('PROTEIN FLEXIBILITY AND ADAPTABILITY SEEN IN 25 CRYSTAL FORMS OF T4 LYSOZYME','')
    title = title.replace('SPECIFICITY OF LIGAND BINDING IN A BURIED NON-POLAR CAVITY OF T4 LYSOZYME: LINKAGE OF DYNAMICS AND STRUCTURAL PLASTICITY','')
    title = title.replace('CONSERVATION OF SOLVENT-BINDING SITES IN 10 CRYSTAL FORMS OF T4 LYSOZYME','')
    if pdb in ['1%iL' %(i) for i in range(55,67,)]:
        title = title.replace('CONTROL OF ENZYME ACTIVITY BY AN ENGINEERED DISULFIDE BOND','')
    if 'LABELED T4 LYSOZYME' in title and 'SPIN' in title:
        res_id2 = ''.join(d_mmCIF['_pdbx_struct_mod_residue.label_comp_id'])
        res_no = int(''.join(list(set(d_mmCIF['_pdbx_struct_mod_residue.label_seq_id']))))
        res1 = d_321[ref_seq[res_no-1]]
        title = title.replace('SPIN-LABELED','SPIN LABELED')
        title = title.replace('T4 LYSOZYME MUTANT','T4 LYSOZYME')
    ##            if 'T115' in title and pdb == '2IGC':
    ##                print pdb
    ##                print title
    ##                print 'SPIN LABELED T4 LYSOZYME (%s%iR1)' %(res1,res_no,)
    ##                stop
        title = title.replace('SPIN LABELED T4 LYSOZYME (%s%i%s)' %(res1,res_no,res_id2,),'')
        title = title.replace('SPIN LABELED T4 LYSOZYME (%s%i%s)' %(res1,res_no,res_id2[:2]),'')
        title = title.replace('SPIN LABELED T4 LYSOZYME %s%i%s' %(res1,res_no,res_id2,),'')
        title = title.replace('SPIN LABELED T4 LYSOZYME %s%i%s' %(res1,res_no,res_id2[:2],),'')

    ## descriptor redundancy
    title = title.replace('REFINEMENT OF AN ENZYME COMPLEX WITH INHIBITOR BOUND AT PARTIAL OCCUPANCY. ','')
    title = title.replace('HEN EGG-WHITE LYSOZYME AT A HYDROSTATIC PRESSURE OF 1000 ATMOSPHERES','')
    title = title.replace('COMPARISON OF BACTERIOPHAGE T4 LYSOZYME AT LOW, MEDIUM, AND HIGH IONIC STRENGTHS','')

    ## name, HEWL        
    title = title.replace('HEN LYSOZYME','') ## 1V7S
    title = title.replace('CHICKEN LYSOZYME','') ## 2lyo
    title = title.replace('HEW LYSOZYME','') ## 1C10
    title = title.replace('CHICKEN EGG WHITE LYSOZYME','') ## 1AZF
    title = title.replace('HEN EGG-WHITE LYSOZYME','')
    title = title.replace('HEN EGG WHITE LYSOZYME','')
    descriptor = descriptor.replace('HEN EGG WHITE LYSOZYME','') ## 2G4Q
    descriptor = descriptor.replace('HEN EGG-WHITE LYSOZYME','')
    descriptor = descriptor.replace('LYSOZYME C','')

    ## name, T4L
    title = title.replace('BACTERIOPHAGE T4 LYSOZYME','')
    title = title.replace('WILDTYPE PHAGE T4 LYSOZYME','') ## 3cdr
    title = title.replace('PHAGE T4 LYSOZYME','')
    title = title.replace('T4 LYSOZYME','')
    title = title.replace(' T4 ','') ## 2rbq
    descriptor = descriptor.replace('PROTEIN (T4 LYSOZYME)','')
    descriptor = descriptor.replace('BACTERIOPHAGE T4 LYSOZYME','')
    descriptor = descriptor.replace('T4 LYSOZYME','')
    descriptor = descriptor.replace('T4-LYSOZYME','')
    descriptor = descriptor.replace('PROTEIN (164-MER)','')

    ## name, both
    if title == 'HYDROLASE':
        title = ''
    if descriptor == 'PROTEIN':
        descriptor = ''
    title = title.replace('LYSOZYME','')
    title = title.replace('(MUCOPEPTIDE N-ACETYLMURAMYL HYDROLASE)','') ## 1HSW
    descriptor = descriptor.replace('(MUCOPEPTIDE N-ACETYLMURAMYL HYDROLASE)','')
    descriptor = descriptor.replace('(E.C.3.2.1.17)','')
    descriptor = descriptor.replace('(E.C. 3.2.1.17)','') ## 1QTK
    descriptor = descriptor.replace('(3.2.1.17)','')
    descriptor = descriptor.replace('PROTEIN (LYSOZYME) (3.2.1.17)','')
    descriptor = descriptor.replace('PROTEIN (LYSOZYME)','')
    descriptor = descriptor.replace('LYSOZYME T4','') ## 1QTD
    descriptor = descriptor.replace('LYSOZYME','') ## 1IR8
    descriptor = descriptor.replace('WILD TYPE','') ## 1hel

    ## prefix
    title = title.replace('THE CRYSTAL STRUCTURES OF ','')
    title = title.replace('THE CRYSTAL STRUCTURE OF ','')
    title = title.replace('CRYSTAL STRUCTURE ANALYSIS OF ','')
    title = title.replace('CRYSTAL STRUCTURE OF ','')
    title = title.replace('X-RAY STRUCTURE OF ','')
    title = title.replace('THE STRUCTURE OF ','')
    title = title.replace('MUTANT STRUCTURE OF ','')
    title = title.replace('THE %.2f A STRUCTURE OF ' %(float(resolution)),'')
    title = title.replace('THE %.1f A STRUCTURE OF ' %(float(resolution)),'')
    title = title.replace('ATOMIC RESOLUTION REFINEMENT OF ','') ## 4lzt
    title = title.replace('REFINEMENT OF ','')
    title = title.replace('AN X-RAY STUDY OF THE STRUCTURE AND BINDING PROPERTIES OF ','')
    ##        if title[:10] == 'STRUCTURE ':
    title = title.replace('ANOMALOUS SUBSTRUCTURE OF ','')
    title = title.replace('STRUCTURE OF ','')

    ## resolution
    title = title.replace(' AT ATOMIC RESOLUTION','')
    title = title.replace(' ANGSTROM RESOLUTION','A')
    title = title.replace(' ANGSTROMS RESOLUTION','A')
    title = title.replace('-ANGSTROMS RESOLUTION','A')
    title = title.replace(' A RESOLUTION','A')
    title = title.replace(' REFINED AT %.1fA' %(float(resolution)),'') ## 2lzm
    title = title.replace(' AT %.1fA' %(float(resolution)),'')
    title = title.replace(' AT %.2fA' %(float(resolution)),'')
    title = title.replace(' AT %.2f A' %(float(resolution)),'')

    ## pressure
    title = title.replace('HIGH-PRESSURE','')
    title = title.replace('HIGH PRESSURE','')
    title = title.replace('AT AMBIENT PRESSURE','')
    title = title.replace('AT 100 MPA','100 MPA')
    title = title.replace('AT 150 MPA','150 MPA')
    title = title.replace('AT 200 MPA','200 MPA')
    ##        title = title.replace(' 2 ATM','')
    ##        title = title.replace(' 4 ATM','')
    ##        title = title.replace(' 8 ATM','')
    ##        title = title.replace(' 16 ATM','')
    ##        title = title.replace(' 32 ATM','')
##    title = title.replace('1 ATMOSPHERE,','   1 ATM')
##    title = title.replace('1000 ATMOSPHERES,','1000 ATM')
    descriptor = descriptor.replace('1 ATMOSPHERE,','   1 ATM')
    descriptor = descriptor.replace('1000 ATMOSPHERES,','1000 ATM')
    title = title.replace('(8 BAR)','')

    ## pH
    title = title.replace(',PH %s' %(pH),'')
    title = title.replace('GROWN AT PH %s' %(pH),'')
    title = title.replace('AT PH %s' %(pH),'')
    title = title.replace('AT PH%s' %(pH),'')
    title = title.replace('(PH %s)' %(pH),'')
    descriptor = descriptor.replace('(PH %s)' %(pH),'')
    if pH == '9.6':
        title = title.replace('BASIC PH','')

    ## Temperature
    descriptor = descriptor.replace('(%s KELVIN)' %(T),'')
    descriptor = descriptor.replace('(%s K)' %(T),'')
    descriptor = descriptor.replace('AT %sK' %(T),'')
    title = title.replace('AT %sK' %(T),'')
    title = title.replace('AT %s K' %(T),'')
    if T != None:
        title = title.replace('AT %i K' %(int(float(T))),'')
    if T in ['298','298.0',]:
        title = title.replace('AT ROOM TEMPERATURE','')
    if T == '100.0':
        title = title.replace('AT LOW TEMPERATURE','')

    ## ligands and experimental conditions
    title = title.replace(' GROWN IN PRESENCE OF ','')
    title = title.replace(' GROWN IN PRESENCE ','')
    title = title.replace(' IN PRESENCE OF ','')
    title = title.replace(' IN THE PRESENCE OF ','')
    title = title.replace(' GROWN AT ',',')
    title = title.replace(' GROWN IN ',',')
    title = title.replace(' UNDER PRESSURE OF ','')
    title = title.replace('IN COMPLEX WITH','')
    title = title.replace(' COMPLEX WITH ','')
    title = title.replace(' SOAKED WITH ','')
    title = title.replace(' SOAKED IN ','')
    title = title.replace(' COCRYSTALLIZED WITH ','')
    title = title.replace(' IN THE HYDROPHOBIC CAVITY OF ','')
    title = title.replace(' BOUND BY ','')
    title = title.replace(' BOUND WITH ','')
    title = title.replace(' BINDING','')
    title = title.replace('APO STRUCTURE','') ## 3gui
    title = title.replace(' APO ','') ## 1qtv
    title = title.replace('FREE OF LIGAND','') ## 3dmv
    title = title.replace('THE INHIBITOR ','')
    title = title.replace('1.4 M NACL','') ## 3lym
    descriptor = descriptor.replace('1.4 M NACL','') ## 3lym
    descriptor = descriptor.replace('CO-CRYSTALLIZED WITH','')
    descriptor = descriptor.replace('COMPLEXED WITH','')
    descriptor = descriptor.replace('IN COMPLEX WITH','')
    descriptor = descriptor.replace('COMPLEX WITH','')
    descriptor = descriptor.replace('BINDING SITE IS OCCUPIED BY ','')
    for i_chem_comp in range(len(d_mmCIF['_chem_comp.id'])):
        chem_comp_id = d_mmCIF['_chem_comp.id'][i_chem_comp]
        chem_comp_name = d_mmCIF['_chem_comp.name'][i_chem_comp]
        chem_comp_synonyms = d_mmCIF['_chem_comp.pdbx_synonyms'][i_chem_comp]
        if chem_comp_name != chem_comp_name.upper():
            stop
        if chem_comp_id in ['AR','KR','XE',]:
            title = title.replace(chem_comp_name,chem_comp_id)
            descriptor = descriptor.replace(chem_comp_name,chem_comp_id)
            continue
        elif chem_comp_id in l_ligands:
            title = title.replace(chem_comp_name,'')
            descriptor = descriptor.replace(chem_comp_name,'')
            title = title.replace(chem_comp_synonyms,'')
            descriptor = descriptor.replace(chem_comp_synonyms,'')
            title = title.replace(chem_comp_name.replace(' ',''),'')
            descriptor = descriptor.replace(chem_comp_name.replace(' ',''),'')
            title = title.replace(chem_comp_name.replace('(2S)-',''),'')
            descriptor = descriptor.replace(chem_comp_name.replace('(2S)-',''),'')
    for chem_comp_id,chem_comp_name in [
        ## differently spelled ligands
        ['JZ9','4,5,6,7-TETRAHYDROINDOLE',],
        ['JZ0','2-METHYLPHENOL',],
        ['BRJ','BROMOETHANOL',],
        ['OXE','O-XYLENE',],
        ['PXY','P-XYLENE',],
        ['PYL','ETHYLBENZENE',],
        ['PIH','IODOBENZENE',],
        ['F5B','PENTAFLUOROBENZENE',],
        ['BBF','BROMOPENTAFLUOROBENZENE',],
        ['BCF','CHLOROPENTAFLUOROBENZENE',],
        ['IBF','IODOPENTAFLUOROBENZENE',],
        ['264','N-PHENYLGLYCINONITRILE',],
        ['263','3-METHYLBENZYLAZIDE',],
        ['ASR','PARA-ARSANILATE',],
        ['266','2-(N-PROPYLTHIO)ETHANOL',],
        ['259','4-(METHYLTHIO)NITROBENZENE',],
        ['260','2,6-DIFLUOROBENZYLBROMIDE',],
        ['F3B','1,3,5-TRIFLUORO-2,4,6-TRICHLOROBENZENE',],
        ['258','BETA-CHLOROPHENETOLE',],
        ['269','(R)(+)-3-CHLORO-1-PHENYL-1-PROPANOL',],
        ['MR3','1-METHYLPYRROLE',],
        ['MM2','CU2-XYLYLBICYCLAM',],
        ['MM6','NI-CYCLAM',],
        ['MM5','NI2-XYLYLBICYCLAM',],
        ## 
        ['I3C','THE MAGIC TRIANGLE',],
        ['EU','TRIS-DIPICOLINATE EU COMPLEX',],
        ['PDC','TRIS-DIPICOLINATE EU COMPLEX',],
        ## error
        ['JZ2','5-CHLORO-2-METHYLPHENOL',],
        ['FLM','3-FLUORO-2-METHYL_ANILINE',],
        ['IOD','TRI-',],
        ['IOD','PERIODATE',], ## 1hc0
        ['IOD','IODINE',], ## 1vat
        ['CMO','',],
        ## stereochemistry changed during remediation
        ['MRD','MPD',],
        ## ions without "ion" at the end
        ['SCN','THIOCYANATE',],
        ['IOD','IODIDE',],
        ['NO3','NITRATE',],
        ['BR','BROMIDE',],
        ['GD','GADOLINIUM',],
        ['XE','XE',],
        ]:
        if chem_comp_id in l_ligands:
            descriptor = descriptor.replace('%s DERIVATIVE OF' %(chem_comp_name),'')
            title = title.replace('%s DERIVATIVE OF' %(chem_comp_name),'')
            descriptor = descriptor.replace('%s SOLUTION' %(chem_comp_name),'')
            title = title.replace('%s SOLUTION' %(chem_comp_name),'')
            descriptor = descriptor.replace('WITH %s' %(chem_comp_name),'')
            title = title.replace('WITH %s' %(chem_comp_name),'')
            descriptor = descriptor.replace('%s COMPLEX' %(chem_comp_name),'')
            title = title.replace('%s COMPLEX' %(chem_comp_name),'')
            descriptor = descriptor.replace(chem_comp_name,'')
            title = title.replace(chem_comp_name,'')
    title = title.replace('ARGON','AR')
    title = title.replace('KRYPTON','KR')
    title = title.replace('XENON','XE')

    ##
    ## space group
    ##
    title = title.replace('THE ORTHORHOMBIC FORM OF ','')
    title = title.replace('THE TETRAGONAL FORM OF ','')
    title = title.replace('THE MONOCLINIC FORM OF ','')
    title = title.replace('THE MONOCLINIC C2 FORM OF ','') ## 1PS5
    title = title.replace('MONOCLINIC ','')
    title = title.replace('TETRAGONAL ','')
    title = title.replace('ORTHORHOMBIC ','')
    title = title.replace('TRICLINIC ','')
    title = title.replace('HEXAGONAL ','')
    descriptor = descriptor.replace('TRICLINIC CRYSTAL FORM','')
    descriptor = descriptor.replace('THE MONOCLINIC C2 FORM OF ','')
    descriptor = descriptor.replace('(SPACE GROUP %s' %(spacegroup),'')

    ##
    ## mutation
    ##
    for mutation in l_mutations:
        title = title.replace(mutation,'')
        descriptor = descriptor.replace(mutation,'')
        if mutation[-1] != 'X':
            title = title.replace('%3s %i REPLACED BY %3s' %(d_123[mutation[0]],int(mutation[1:-1]),d_123[mutation[-1]],),'')
            descriptor = descriptor.replace('%3s %i REPLACED BY %3s' %(d_123[mutation[0]],int(mutation[1:-1]),d_123[mutation[-1]],),'')
            descriptor = descriptor.replace('%3s %i %3s' %(d_123[mutation[0]],int(mutation[1:-1]),d_123[mutation[-1]],),'')

    title = title.replace('DOUBLE MUTANT WITH','')
    title = title.replace('CAVITY MUTANT','')
    title = title.replace('CAVITY CREATING MUTANT','')
    title = title.replace('CAVITY CREATING MUTATION','')
    title = title.replace('MUTANT OF','')
    title = title.replace('MUTANT WITH','')
    title = title.replace('MUTANT','')
    title = title.replace('SELENO VERSION','')
    title = title.replace('IN WILDTYPE BACKGROUND','')
    title = title.replace('SYNTHETIC DIMER','')
    if 'WT*' in title or 'PSEUDO-WT' in title or '/TA' in title:
        if 'C54T' in l_mutations and 'C97A' in l_mutations:
            title = title.replace('WT*','')
            title = title.replace('PSEUDO-WT','')
            title = title.replace('/TA','')
    descriptor = descriptor.replace('DOUBLE MUTANT WITH','') ## 1l01
    descriptor = descriptor.replace('SUBSTITUTIONS','')
    descriptor = descriptor.replace('MUTANT WITH','')
    descriptor = descriptor.replace('MUTANT','')

    ## common words
    title = title.replace('CRYSTAL','')

    ## errors
    if pdb == '1L35':
        descriptor = descriptor.replace('ILE 3 REPLACED BY TYR','')
    if pdb in ['1C69','1C6A','1C6B',]:
        title = title.replace('C54T/C97A','')
    if pdb == '1T8G':
        title = title.replace('T34A','')
    if pdb == '2OE4':
        title = title.replace('PSUEDO WILD TYPE','')

    descriptor = descriptor.replace('HIGH SALT','HI SALT')
    descriptor = descriptor.replace('MEDIUM SALT','MED SALT')
    descriptor = descriptor.replace('LOW SALT','LO SALT')
    descriptor = descriptor.replace('DITHIOTHREITOL','DTT')

    ## parentheses
    descriptor = descriptor.replace('()','')

    ## strip 1 (before common words)
    descriptor = descriptor.strip()
    title = title.strip()

    ## common words
    d = {'descriptor':descriptor,'title':title,}
    for k in ['descriptor','title',]:
        bool_special_characters_only = True
        for s in d[k].replace('AND','').replace('THE','').replace('OF','').replace('AS CONTROL','').replace('CRYSTAL',''):
            if s not in [',',' ','(',')','/','-',]:
                bool_special_characters_only = False
                break
        if bool_special_characters_only == True:
            d[k] = ''
    descriptor = d['descriptor']
    title = d['title']

    ## strip 2 (before redundancy)
    descriptor = descriptor.strip()
    title = title.strip()

    ## redundancy
    title = title.replace(descriptor,'')
    descriptor = descriptor.replace(title,'')

    if len(descriptor) > 0:
        ## leading comma
        if descriptor[0] == ',':
            descriptor = descriptor[1:]
    if len(descriptor) > 0:
        ## parentheses
        if descriptor[0] == '(' and descriptor[-1] == ')':
            descriptor = descriptor[1:-1]
    ## leading slash
    if len(title) > 0:
        while title[0] == '/':
            title = title[1:]
##        ## leading slash
##        if len(descriptor) > 0:
##            while descriptor[0] == '//':
##                descriptor = descriptor[1:]
    
    if descriptor != '':
        print pdb, descriptor

    return title, descriptor


def calc_chi1(l_pdbs,d_coordinates):

    print 'calc chi1'

    import dihedral

    d_chi1 = {}
    d_fourth = {
        'CYS':'SG','ILE':'CG1','SER':'OG','THR':'OG1',
        'VAL':'CG1', ## could easily be flipped
        'CME':'SG',
        }

    for res_no in range(1,1+len(ref_seq)):

##        print 'calc chi1, residue', res_no

        if res_no in [163,164,]:
            continue

        d_chi1[res_no] = {}

        for pdb in l_pdbs:

            if not res_no in d_coordinates[pdb].keys():
                if not (ref_seq[res_no-1] == 'VAL' and bool_exclusion_symmetry == False):
                    print 'unobserved residue', pdb, res_no
                continue

            res_name = d_coordinates[pdb][res_no]['res_name']

            ## no chi1 angle
            if res_name in ['GLY','ALA',]:
                continue
##            ## cg1 and cg2 can't be distinguished
##            if res_name == 'VAL':
##                continue
            ## modified residue
            if res_name not in ['MSE','CME',]+d_321.keys():
                continue
            if res_name in d_fourth.keys():
                atom4 = d_fourth[res_name]
            else:
                atom4 = 'CG'

            ## side chain atoms unobserved or modified residue
            if len( set(['N','CA','CB',atom4,]) & set(d_coordinates[pdb][res_no]['atoms'].keys()) ) < 4:
                if res_name != 'VAL':
                    print 'unobserved atoms', pdb, res_name, res_no, set(['N','CA','CB',atom4,]) - set(d_coordinates[pdb][res_no]['atoms'].keys())
                continue

            ## high temp factors
            if bool_exclusion_high_temp_factors == True:
                if float(d_coordinates[pdb][res_no]['atoms']['N']['tempFactor']) > max_bfactor:
                    continue
                if float(d_coordinates[pdb][res_no]['atoms']['CA']['tempFactor']) > max_bfactor:
                    continue
                if float(d_coordinates[pdb][res_no]['atoms']['CB']['tempFactor']) > max_bfactor:
                    continue
                if float(d_coordinates[pdb][res_no]['atoms'][atom4]['tempFactor']) > max_bfactor:
                    continue
            
            c1 = d_coordinates[pdb][res_no]['atoms']['N']['coord']
            c2 = d_coordinates[pdb][res_no]['atoms']['CA']['coord']
            c3 = d_coordinates[pdb][res_no]['atoms']['CB']['coord']
            c4 = d_coordinates[pdb][res_no]['atoms'][atom4]['coord']
            chi1 = dihedral.main(c1,c2,c3,c4)

            d_chi1[res_no][pdb] = chi1

    fd = open('d_chi1_%s_%s.txt' %(protein,suffix_exclusion,),'w')
    fd.write(str(d_chi1))
    fd.close()

    return d_chi1


def plot_dihedral(d_coordinates,d_mmCIF,ref_seq,l_pdbs,l_range,l_mutants,d_mutants,):

    import dihedral

    for res_no in l_range:

        ## no phi angle
        if res_no == 1:
            continue
        
        res_index = res_no-1
        print 'ramachandran', res_no
        l = []
        for pdb in l_pdbs:

##            ## MODRES
##            if pdb in [
##                ## HEWL
##                '132L_A', # methylated lysine (DM0, ...)
##                '1AT5_A', # succinimide (SNN)
##                '1AT6_A', # isoaspartate (ASP)
##                '1RCM_A','1RCM_B', # carboxymethylated cysteine (CCS)
####                ## T4L
####                '2NTH_A',
##                ]:
##                continue

            ## missing residues
            if res_no-1 not in d_coordinates[pdb].keys():
                continue
            if res_no not in d_coordinates[pdb].keys():
                continue
            if res_no+1 not in d_coordinates[pdb].keys():
                continue

            l_res_name_prev = d_coordinates[pdb][res_no-1].keys()
            if len(l_res_name_prev) == 1:
                res_name_prev = ''.join(l_res_name_prev)
            else:
                res_name_prev = ref_seq[res_index-1]

            l_res_name_curr = d_coordinates[pdb][res_no].keys()
            if len(l_res_name_curr) == 1:
                res_name_curr = ''.join(l_res_name_curr)
            else:
                res_name_curr = ref_seq[res_index]

            l_res_name_next = d_coordinates[pdb][res_no+1].keys()
            if len(l_res_name_next) == 1:
                res_name_next = ''.join(l_res_name_next)
            else:
                res_name_next = ref_seq[res_index+1]

            c1 = d_coordinates[pdb][res_index][res_name_prev]['C']['coord']
            c2 = d_coordinates[pdb][res_no][res_name_curr]['N']['coord']
            c3 = d_coordinates[pdb][res_no][res_name_curr]['CA']['coord']
            c4 = d_coordinates[pdb][res_no][res_name_curr]['C']['coord']
            phi = dihedral.main(c1,c2,c3,c4)

            c1 = d_coordinates[pdb][res_no][res_name_curr]['N']['coord']
            c2 = d_coordinates[pdb][res_no][res_name_curr]['CA']['coord']
            c3 = d_coordinates[pdb][res_no][res_name_curr]['C']['coord']
            c4 = d_coordinates[pdb][res_no+1][res_name_next]['N']['coord']
            psi = dihedral.main(c1,c2,c3,c4)

            spacegroup = core.parse_mmCIF_item(d_mmCIF[pdb[:4]],'_symmetry.space_group_name_H-M',pdb,)

            ## iodide/thiocyanate
            if pdb in ['1B2K_A','1LCN_B','1LKR_B',]:
                l += ['%s NA NA %s %s\n' %(pdb,phi,psi,)]
            ## dehydration
            elif protein == 'HEWL' and pdb in [
                '1V7T_A','1V7T_B','1XEK_A','2Z18_A','2D4J_A','1XEI_A','1XEJ_A','2Z12_A','2Z19_A','1LMA_A',
                ]:
                l += ['%s NA NA NA NA %s %s\n' %(pdb,phi,psi,)]
            elif protein == 'T4L' and spacegroup != 'P 32 2 1':
                l += ['%s NA NA NA NA %s %s\n' %(pdb,phi,psi,)]
            elif pdb in l_mutants:
                l += ['%s NA NA NA NA NA NA %s %s\n' %(pdb,phi,psi,)]
            ## all other            
            else:
                l += ['%s %s %s\n' %(pdb,phi,psi,)]

        if l != []:

            prefix = 'dihedral/dihedral_%s_%i' %(protein,res_no)
            fd = open('%s.gnuplotdata' %(prefix),'w')
            fd.writelines(l)
            fd.close()

            s_plot = 'plot [-180:180][-180:180] '
            if protein == 'HEWL':
                s_plot += '"%s.gnuplotdata" u 4:5:1 pt 7 ps 2 lc 1 t "%s", "%s.gnuplotdata" u 4:5:1 w labels t "", ' %(prefix, 'P21ions', prefix,)
                s_plot += '"%s.gnuplotdata" u 6:7:1 pt 7 ps 2 lc 2 t "%s", "%s.gnuplotdata" u 6:7:1 w labels t "", ' %(prefix, 'dehydrated', prefix,)
            elif protein == 'T4L':
                s_plot += '"%s.gnuplotdata" u 6:7:1 pt 7 ps 2 lc 2 t "%s", "%s.gnuplotdata" u 6:7:1 w labels t "", ' %(prefix, 'not P 32 2 1', prefix,)
            s_wt = '"%s.gnuplotdata" u 2:3:1 w labels t "", "%s.gnuplotdata" u 2:3:1 pt 7 ps 2 lc 3 t "%s", ' %(
                prefix, prefix,'WTs',
                )
            if l_mutants != []:
                s_mut = '"%s.gnuplotdata" u 8:9:1 w labels t "", "%s.gnuplotdata" u 8:9:1 pt 7 ps 2 lc 4 t "%s", ' %(
                    prefix, prefix,
                    'mutants'
                    )
            if protein == 'T4L':
                s_plot += s_mut+s_wt
            elif protein == 'HEWL':
                s_plot += s_wt+s_mut
            s_plot = s_plot[:-2]+'\n'
                
            gnuplot.scatter_plot_2d(
                prefix,
                bool_multiple_columns = True,
                s_plot = s_plot,
                size = 'square',
                bool_labels = True,
                xmin = -180, xmax = 180,
                ymin = -180, ymax = 180,
                xlabel = 'phi',
                ylabel = 'psi',
                title = '%s Ramachandran %3s%i' %(protein,ref_seq[res_no-1], res_no,),
                bool_remove = True,
                )

            os.system('rm %s.ps' %(prefix))
##            os.system('mv %s.png dihedral/.' %(prefix))

    return


def matthews(d_mmCIF,l_pdbs,):

    import matthews_coefficient

    d_matthews = {}
    d_Z = {}

    if protein == 'HEWL':
        mw = 14331.238
    elif protein == 'T4L':
        mw = 18662.609
    l = []
    MVmax = 0
    for pdb in l_pdbs:

        pdb = pdb[:4]

        if 'matthews' in d_mmCIF[pdb].keys():
            continue

        a = float(d_mmCIF[pdb]['_cell.length_a'][0])
        b = float(d_mmCIF[pdb]['_cell.length_b'][0])
        c = float(d_mmCIF[pdb]['_cell.length_c'][0])
        alpha = float(d_mmCIF[pdb]['_cell.angle_alpha'][0])
        beta = float(d_mmCIF[pdb]['_cell.angle_beta'][0])
        gamma = float(d_mmCIF[pdb]['_cell.angle_gamma'][0])
        Z = int(d_mmCIF[pdb]['_cell.Z_PDB'][0]) ## number of polymers in unit cell
        MV = matthews_coefficient.main(a,b,c,alpha,beta,gamma,mw,Z,)

        ## append
        l += ['%f %f %s\n' %(float(d_mmCIF[pdb]['_refine.ls_d_res_high'][0]), MV, pdb,)]
        d_mmCIF[pdb]['matthews'] = [MV]

        d_matthews[pdb] = MV
        d_Z[pdb] = Z
        
        if MV > MVmax:
            MVmax = MV

    prefix = 'matthews_%s' %(protein)

    fd = open('%s.gnuplotdata' %(prefix),'w')
    fd.writelines(l)
    fd.close()

    gnuplot.scatter_plot_2d(
        prefix,
        xlabel = 'Resolution (Angstrom)',
        ylabel = 'Matthews Coefficient (Angstrom^3 / Da)',
        col_label = 3,
        function = '1.55+0.1*x',
        xmin = 0,
        xmax = resolution_min,
        ymax = MVmax+0.1,
        ymin = 0,
        )

    os.system('rm %s.gnuplot*' %(prefix))
    os.system('rm %s.ps' %(prefix))

    return d_mmCIF, d_matthews, d_Z


def plot_start_model_vs_derived_wt_and_derived_mut(d_rmsds_subset):

    print 'plot start model vs derived wt and derived mut'

    prefix = 'mut_and_wt_v_startmodel_%s_%s' %(protein,method,)

    l_groups = ['flexible','fixed',]

    d = {
        '1RFP_A': {'wt': ['1UIH_A'], 'mut': ['1FN5_A', '1FLQ_A', '1IOR_A', '1FLW_A', '1FLY_A', '1IOQ_A', '1FLU_A', '1IOT_A', '1IOS_A'],},
        '6LYZ_A': {'wt': ['1AZF_A'], 'mut': ['1LZD_A', '1LZE_A', '1LZG_A'],}
        ## poor resolution
##        '2LZH_A': {'wt': ['1AKI_A'], 'mut': ['1HEQ_A', '1HEO_A', '1HEP_A', '1HEN_A', '1HEM_A', '1HER_A'],},
        }

    l = []
    for startmodel in d.keys():
        s = '%s ' %(startmodel)
        pdb_wt = d[startmodel]['wt'][0]
        l_pdbs_mut = d[startmodel]['mut']
        for group in l_groups:

            ## no need to do average, as there is only one wt
            rmsd_wt = d_rmsds_subset[startmodel][pdb_wt][group]

            l_rmsds_mut = []
            for pdb_mut in l_pdbs_mut:
                rmsd_mut = d_rmsds_subset[startmodel][pdb_mut][group]
                l_rmsds_mut += [rmsd_mut]

            average, stddev = statistics.do_stddev(l_rmsds_mut)
            s += '%s %s %s ' %(rmsd_wt,average,stddev,)
        s += '\n'
        l += [s]

    fd = open('%s.gnuplotdata' %(prefix,),'w')
    fd.writelines(l)
    fd.close()


    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(prefix,),
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
##            'set encoding iso_8859_1\n', ## postscript encoding for special characters
        'set xlabel "%s"\n' %('<RMSD_s_t_a_r_t_m_o_d_e_l _v _d_e_r_i_v_e_d _w_t_s>'),
        'set ylabel "%s"\n' %('<RMSD_s_t_a_r_t_m_o_d_e_l _v _d_e_r_i_v_e_d _m_u_t_a_n_t_s>'),
        'set title "%s"\n' %(method),
        'f(x) = x\n',
    ]
    
    line_plot = 'plot f(x) lc 0 t "", '
    for i in range(len(l_groups)):
        line_plot += '"%s.gnuplotdata" u %i:%i:%i t "%s" w errorb ps 2 pt 7 lt 1 lc %i lw 2, ' %(
            ## input
            prefix,
            ## u
            2+3*i,3+3*i,4+3*i,
            ## t
            l_groups[i],
            ## lc
            i+1,
            )
        line_plot += '"%s.gnuplotdata" u %i:%i:%i w labels font "Helevetica,36" left offset 0.25,0.25 t "", ' %(
            ## input
            prefix,
            ## x,y
            2+3*i,3+3*i,
            ## label
            1,
            )
    line_plot = line_plot[:-2]+'\n'

    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix,),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(prefix,))

    os.system('convert %s.ps %s.png' %(prefix,prefix,))

    os.system('rm %s.ps' %(prefix,))
    os.system('rm %s.gnuplot*' %(prefix,))

    return


def calculate_rmsd_per_subset(l_pdbs,d_properties,d_rmsds_overall,d_coordinates,):

    print 'calc rmsd per subset'

    d_rmsds_subset = {}
    for i1 in range(len(l_pdbs)-1):
        pdb1 = l_pdbs[i1]
        print 'calc rmsd subset', pdb1
        if not pdb1 in d_rmsds_subset.keys():
            d_rmsds_subset[pdb1] = {}
        for i2 in range(i1+1,len(l_pdbs)):
            pdb2 = l_pdbs[i2]
            d_rmsds_subset[pdb1][pdb2] = {}
            if not pdb2 in d_rmsds_subset.keys():
                d_rmsds_subset[pdb2] = {}
            d_rmsds_subset[pdb2][pdb1] = {}
            for group in d_properties['flexibility'].keys():
                l_coordinates1 = []
                l_coordinates2 = []
                for res_no in d_properties['flexibility'][group]:
                    
                    ## missing residues
                    if not res_no in d_coordinates[pdb1].keys():
                        continue
                    if not res_no in d_coordinates[pdb2].keys():
                        continue
                    
                    l_coords1,l_coords2,l_bfacs1,l_bfacs2 = core.get_coordinates(
                        method, d_coordinates, res_no,
                        pdb1, pdb2,
                        )
                    l_coordinates1 += l_coords1
                    l_coordinates2 += l_coords2
                tv1 = d_rmsds_overall[pdb1][pdb2]['tv1']
                rm = d_rmsds_overall[pdb1][pdb2]['rm']
                tv2 = d_rmsds_overall[pdb1][pdb2]['tv2']
                for i_coord in range(len(l_coordinates2)):
                    coord = l_coordinates2[i_coord]
                    coord = numpy.dot(coord-tv1,rm)+tv2
                    l_coordinates2[i_coord] = coord
                rmsd = core.calc_rmsd_of_pre_aligned_coordinates(l_coordinates1,l_coordinates2,)
                d_rmsds_subset[pdb1][pdb2][group] = rmsd
                d_rmsds_subset[pdb2][pdb1][group] = rmsd

    return d_rmsds_subset


def plot_startmodel_vs_wts_and_derived_wts(
    l_pdbs,
    d_properties,d_startingmodel,d_coordinates,
    d_rmsds_overall,
    d_mmCIF,
    ):

    print 'plot per startmodel'

    print 'calculate rmsd per fixed/flexible group'

    d_rmsds_subset = calculate_rmsd_per_subset(l_pdbs,d_properties,d_rmsds_overall,d_coordinates,)

    d_chaincount = {}
    for pdb in l_pdbs:
        for i_entity in range(len(d_mmCIF[pdb[:4]]['_entity_poly.entity_id'])):
            entity_id = d_mmCIF[pdb[:4]]['_entity_poly.entity_id'][i_entity]
            s_chain_ids = d_mmCIF[pdb[:4]]['_entity_poly.pdbx_strand_id'][i_entity]
            if len(s_chain_ids) in [2,4,6,8,10,]:
                stop
            if pdb[-1] in s_chain_ids:
                l_chain_ids = s_chain_ids.split(',')
                break
        d_chaincount[pdb] = len(l_chain_ids)
        

##    d_statistics = {}
##    for group in d_properties['flexibility'].keys():
##        l_rmsds = []
##        for i1 in range(len(l_pdbs)-1):
##            pdb1 = l_pdbs[i1]
##            for i2 in range(i1+1,len(l_pdbs)):
##                pdb2 = l_pdbs[i2]
##                rmsd = d_rmsds_subset[pdb1][pdb2][group]
##                l_rmsds += [rmsd]
##                average, stddev = statistics.do_stddev(l_rmsds)
##                d_statistics[group] = {'average':average,'stddev':stddev,}

    l_startingmodels = parse_list_of_starting_models(l_pdbs,d_startingmodel,d_mmCIF,)
    l_groups = [
##        'all',
        'flexible','fixed',
        ]
    l = []
    ## loop over starting models
    for i in range(len(l_startingmodels)):
        pdb1 = l_startingmodels[i]

        spacegroup1 = core.parse_mmCIF_item(d_mmCIF[pdb1[:4]],'_symmetry.space_group_name_H-M',pdb1,)

        if d_chaincount[pdb1] > 1:
            continue

        s = ' '
        ## flexible or fixed
        for group in l_groups:
            l_rmsds_all = []
            l_rmsds_derivedmodels = []
            ## loop over derived models
            for pdb2 in l_pdbs:

                if pdb1 == pdb2:
                    continue

                spacegroup2 = core.parse_mmCIF_item(d_mmCIF[pdb2[:4]],'_symmetry.space_group_name_H-M',pdb2,)

                if spacegroup1 != spacegroup2:
                    continue

                if d_chaincount[pdb2] > 1:
                    continue

                if pdb2 in [
                    ## monovalent ions
                    '1LCN_B','1LKR_B','1B2K_A',
                    ## dehydrated
                    '1V7T_A','1V7T_B','1XEK_A','2Z18_A','2D4J_A','1XEI_A','1XEJ_A','2Z12_A','2Z19_A','1LMA_A',
                    ]:
                    continue
    
                if group == 'all':
                    rmsd = d_rmsds_overall[pdb1][pdb2]['rmsd']
                else:
                    rmsd = d_rmsds_subset[pdb1][pdb2][group]

                l_rmsds_all += [rmsd]

                bool_derived = False
                if pdb2 in d_startingmodel[pdb1[:4]]:
                    bool_derived = True
                if bool_derived == True:
                    l_rmsds_derivedmodels += [rmsd]

            ## too few structures based on starting model (pdb1)
            if len(l_rmsds_derivedmodels) <= 1:
                continue

            average_all, stddev_all = statistics.do_stddev(l_rmsds_all)

            try:
                average_derived, stddev_derived = statistics.do_stddev(l_rmsds_derivedmodels)
            except:
                print l_rmsds_derivedmodels
                print pdb1, pdb2
                print d_startingmodel[pdb1[:4]]
                stop
            s += ' %s %s %s %s %s ' %(
                ## x,y
                average_all,average_derived,
                ## xyerrorbar
                stddev_all,stddev_derived,
                ## label
                pdb1.replace('_',':'),
                )
            if average_all < average_derived:
                print group, pdb1,
                print round(average_all,2), round(average_derived,2),
                print d_startingmodel[pdb1[:4]], l_rmsds_derivedmodels
        s += '\n'
        l += [s]

    prefix = 'startmodel_%s_%s' %(protein,method,)

    fd = open('%s.gnuplotdata' %(prefix),'w')
    fd.writelines(l)
    fd.close()

    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(prefix),
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
##            'set encoding iso_8859_1\n', ## postscript encoding for special characters
        'set xlabel "%s"\n' %('<RMSD_s_t_a_r_t_m_o_d_e_l _v _w_t_s>'),
        'set ylabel "%s"\n' %('<RMSD_s_t_a_r_t_m_o_d_e_l _v _d_e_r_i_v_e_d _m_o_d_e_l_s>'),
        'set title "%s"\n' %(method),
        'f(x) = x\n',
    ]
    
    line_plot = 'plot f(x) lc 0 t "", '
    for i in range(len(l_groups)):
        line_plot += '"%s.gnuplotdata" u %i:%i:%i:%i t "%s" w xyerrorbars ps 2 pt 7 lt 1 lc %i lw 2, ' %(
            ## input
            prefix,
            ## u
            1+5*i,2+5*i,3+5*i,4+5*i,
            ## t
            l_groups[i],
            ## lc
            i+1,
            )
        line_plot += '"%s.gnuplotdata" u %i:%i:%i w labels font "Helevetica,36" left offset 0.25,0.25 t "", ' %(
            prefix, 1+5*i,2+5*i,5+5*i,
            )
    line_plot = line_plot[:-2]+'\n'

    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(prefix))

    os.system('convert %s.ps %s.png' %(prefix,prefix,))

    print '/software/bin/gnuplot %s.gnuplotsettings' %(prefix)
    print 'convert %s.ps %s.png' %(prefix,prefix,)
    print method

    os.system('rm %s.ps' %(prefix))
    os.system('rm %s.gnuplot*' %(prefix))

    return d_rmsds_subset


def parse_list_of_starting_models(l_pdbs,d_startingmodel,d_mmCIF,):

    '''these exclusions should not precede the selection of valid models'''
    '''instead if too few derived models (rmsds) are available, then skip average/stddev calc'''

    ## which chains in which pdbs?
    d_chains = {}
    for pdb in l_pdbs:
        if not pdb[:4] in d_chains.keys():
            d_chains[pdb[:4]] = []
        d_chains[pdb[:4]] += [pdb]

    l_startingmodels = []
    for pdb in d_startingmodel.keys():
        ##
        if pdb == None:
            continue
        ## starting model has been excluded previously
        if not pdb in d_chains.keys():
            print '~~~~~~~~~~~~~', pdb, d_startingmodel[pdb]
            continue

        spacegroup_start = core.parse_mmCIF_item(d_mmCIF[pdb[:4]],'_symmetry.space_group_name_H-M',pdb,)

##        ## chains in P21 space group not equivalent (different crystal contacts)
##        if spacegroup_start == 'P 1 21 1':
##            continue

        ##
        set_derived_models = set(d_startingmodel[pdb[:4]])

        ## exclude selected chains from P21 space group
        set_derived_models -= set([
            ## IOD,SCN,NaCl
            '1LCN_B','1LKR_B','1B2K_A',
            '2Z18_A',
            ])

        ## exclude derived models from a different space group
        l_diff_spacegroup = []
        for pdb_derived in set_derived_models:
            spacegroup_derived = core.parse_mmCIF_item(d_mmCIF[pdb_derived[:4]],'_symmetry.space_group_name_H-M',pdb_derived,)
            ## exclude derived models from different space groups
            if spacegroup_start != spacegroup_derived:
                l_diff_spacegroup += [pdb_derived]
        set_derived_models -= set(l_diff_spacegroup)
            
        ## no stddev if less than two pdbs based on starting model
        if len(set_derived_models) < 2:
            continue
        
        ## append all chains of starting model
        l_startingmodels += d_chains[pdb]

    l_startingmodels.sort()
    if None in l_startingmodels:
        l_startingmodels.remove(None)

    return l_startingmodels


def plot_rmsds_per_residue_per_startingmodel(
    d_rmsds_per_residue,
    d_startingmodel, d_startingmodel_reverse,
    d_coordinates,
    d_rmsds_overall,
    l_pdbs,
    d_mmCIF,
    ):

    print 'plot'

    l_colors = [
        '000000','FF0000','FF8000','FFFF00','80FF00','00FF00','00FF80','00FFFF','0080FF','0000FF','8000FF','FF00FF','FF0080','DFDFDF','BFBFBF','9F9F9F','7F7F7F',
        ]

    l_startingmodels = parse_list_of_starting_models(l_pdbs,d_startingmodel,d_mmCIF,)
               
    l_columns = [[2,3,'wt',]]
    for i in range(len(l_startingmodels)):
        startingmodel = l_startingmodels[i]
        l_columns += [[2*(i+2),2*(i+2)+1,startingmodel,]]
        
    l = []
    for res_no in range(1,129+1):
        xxx
        l_rmsds = d_rmsds_per_residue[res_no]
        average, stddev = statistics.do_stddev(l_rmsds)
        s = '%s %s %s' %(res_no,average,stddev,)
        for i in range(len(l_startingmodels)):
            startingmodel = l_startingmodels[i]

            l_rmsds = []
            for model in d_startingmodel[startingmodel[:4]]:
                if startingmodel == model:
                    continue

                ## missing residues
                if not res_no in d_coordinates[model].keys():
                    continue
                if not res_no in d_coordinates[startingmodel].keys():
                    continue
                
                l_coords1,l_coords2,l_bfacs1,l_bfacs2 = core.get_coordinates(
                    method, d_coordinates, res_no,
                    model, startingmodel,
                    )
                tv1 = d_rmsds_overall[model][startingmodel]['tv1']
                rm = d_rmsds_overall[model][startingmodel]['rm']
                tv2 = d_rmsds_overall[model][startingmodel]['tv2']
                for i_coord in range(len(l_coords2)):
                    coord = l_coords2[i_coord]
                    coord = numpy.dot(coord-tv1,rm)+tv2
                    l_coords2[i_coord] = coord
                rmsd_residue = core.calc_rmsd_of_pre_aligned_coordinates(l_coords1,l_coords2,)
                ## tmp check of correct transformations
                if rmsd_residue > 11: ## 1lcn,5lym
                    print model, startingmodel, res_no, rmsd_residue
                    stop
                l_rmsds += [rmsd_residue]
            if len(l_rmsds) == 0:
                stop1
            if len(l_rmsds) == 1:
                print startingmodel
                print d_startingmodel[startingmodel[:4]]
                stop2
            average, stddev = statistics.do_stddev(l_rmsds)
            s += ' %s %s' %(average,stddev,)
        s += '\n'
        l += [s]

    prefix = 'rmsd_per_residue_per_startingmodel_%s' %(protein)

    fd = open('%s.gnuplotdata' %(prefix),'w')
    fd.writelines(l)
    fd.close()

    gnuplot.scatter_plot_2d(
        prefix,
        bool_errorbars = True,
        bool_multiple_columns = True,
        l_columns=l_columns,
        l_colors=l_colors,
        )

    return


def plot_rmsds_per_residue(d_rmsds_per_residue):

    l = []
    for res_no in range(1,129+1):
        xxx
        l_rmsds = d_rmsds_per_residue[res_no]
        average, stddev = statistics.do_stddev(l_rmsds)
        l += ['%s %s %s\n' %(res_no,average,stddev,)]

    prefix = 'rmsd_per_residue_%s' %(protein)

    fd = open('%s.gnuplotdata' %(prefix),'w')
    fd.writelines(l)
    fd.close()

    gnuplot.scatter_plot_2d(
        prefix,
        bool_errorbars=True,
        )

    return


def parse_properties(d_mmCIF_main,l_pdbs,):

    d_ph = {}
    d_ph_reverse = {}
    d_temperature = {}
    d_temperature_reverse = {}
    d_spacegroups = {}
    d_spacegroups_reverse = {}
    d_startingmodel = {} ## k=model,v=pdbs
    d_startingmodel_reverse = {} ## k=pdb,v=model
    d_resolutions = {}
    d_resolutions_reverse = {}
    d_countentities = {}
    d_countentities_reverse = {}
    d_author = {}
    d_author_reverse = {}

    for pdb in l_pdbs:

        d_mmCIF = d_mmCIF_main[pdb[:4]]

        spacegroup = core.parse_mmCIF_item(d_mmCIF,'_symmetry.space_group_name_H-M',pdb,)
        T = core.parse_mmCIF_item(d_mmCIF,'_diffrn.ambient_temp',pdb,)
        pH = core.parse_mmCIF_item(d_mmCIF,'_exptl_crystal_grow.pH',pdb,)
        startingmodel = core.parse_mmCIF_item(d_mmCIF,'_refine.pdbx_starting_model',pdb,)
        resolution = core.parse_mmCIF_item(d_mmCIF,'_refine.ls_d_res_high',pdb,)
        n_entities = len(d_mmCIF['_entity_poly.entity_id'])

        l_authors = []
        for i_citation_author in range(len(d_mmCIF['_citation_author.name'])):
            if d_mmCIF['_citation_author.citation_id'][i_citation_author] == 'primary':
                author = d_mmCIF['_citation_author.name'][i_citation_author]
                l_authors += [author]
        if protein == 'T4L' and 'Matthews, B.W.' in d_mmCIF['_citation_author.name']:
            author = 'Matthews, B.W.'
        else:
            author = l_authors[-1]
        author = l_authors

        if T != None:
            T = int(10*round(float(T)/10.))
        if pH != None:
            pH = round(float(pH),0)

##        if startingmodel == None:
##            startingmodel = pdb[:4]

        for prop,d_prop,d_prop_rev in [
            [pH,d_ph,d_ph_reverse,],
            [T,d_temperature,d_temperature_reverse,],
            [spacegroup,d_spacegroups,d_spacegroups_reverse,],
            [startingmodel,d_startingmodel,d_startingmodel_reverse,],
            [resolution,d_resolutions,d_resolutions_reverse,],
            [n_entities,d_countentities,d_countentities_reverse,],
            [author,d_author,d_author_reverse,],
            ]:
            d_prop_rev[pdb] = prop

            if prop != author:
                if not prop in d_prop.keys():
                    d_prop[prop] = []
                d_prop[prop] += [pdb]

    return (
        d_ph, d_ph_reverse,
        d_temperature, d_temperature_reverse,
        d_startingmodel, d_startingmodel_reverse,
        d_spacegroups, d_spacegroups_reverse,
        d_resolutions,d_resolutions_reverse,
        d_author, d_author_reverse,
        d_countentities, d_countentities_reverse,
        )


def do_contour_plot(d_rmsds,l_pdbs_grouped,d_reverse,prefix,l_wts,):

    print 'plot', prefix

    l_data = []
    d_xtics = {}
    for i1 in range(len(l_pdbs_grouped)):
        pdb1 = l_pdbs_grouped[i1]
        if pdb1 in l_wts:
            suffix_xtic = 'wt'
        else:
            suffix_xtic = 'mutant'
        d_xtics['%10s %s %s' %(d_reverse[pdb1],pdb1.replace('_',':'),suffix_xtic,)] = i1+.5
        for i2 in range(len(l_pdbs_grouped)):
            pdb2 = l_pdbs_grouped[i2]
            if pdb1 == pdb2:
                rmsd = 0
##            elif pdb1 not in d_rmsds.keys(): ## tmp!!
##                rmsd = 0
##            elif pdb2 not in d_rmsds[pdb1].keys(): ## tmp!!
##                rmsd = 0
            else:
                rmsd = d_rmsds[pdb1][pdb2]['rmsd']
            l_data += ['%i %i %f\n' %(i1,i2,rmsd,)]
        l_data += ['%i %i %f\n' %(i1,i2+1,rmsd,)]
        l_data += ['\n',]
    for i2 in range(len(l_pdbs_grouped)):
        l_data += ['%i %i %f\n' %(i1+1,i2,rmsd,)]
    l_data += ['%i %i %f\n' %(i1+1,i2+1,rmsd,)]

    gnuplot_prefix = '%s_%s_%s_%s_%s_%s_%s' %(protein,method,'ligands%s' %(bool_ligands),'mutations_max%s' %(n_mutations_max),'multipleentities%s' %(bool_multiple_entities),'spacegroup%s' %(spacegroup), prefix,)
    gnuplot.contour_plot(
        gnuplot_prefix, l_data,
        d_xtics=d_xtics, d_ytics=d_xtics,
        title = 'HEWL'
##        z1=0, z2=1.5,
##        bool_remove = False,
        )

    os.system('mv %s.png sphere/.' %(gnuplot_prefix,))

    return


def trace_starting_model(d_startingmodel,d_startingmodel_reverse,):

    print 'trace starting model'


    for k,v in d_startingmodel_reverse.items():
        del d_startingmodel_reverse[k]
        d_startingmodel_reverse[k[:4]] = v


    for startmodel in d_startingmodel.keys():
        ## no startmodel
        if startmodel == None:
            d_startingmodel[startmodel] += [startmodel]
        else:
            startmodel = startmodel[:4]
        startmodel1 = startmodel
        while True:

            ## current starting model has no starting model itself
            if not startmodel1 in d_startingmodel_reverse.keys():
                startmodel2 = startmodel1
            else:
                startmodel2 = d_startingmodel_reverse[startmodel1]

            ## current starting model has no starting model itself
            if startmodel2 == None or startmodel1 == startmodel2:
                ## replace if startmodel is a model of another startmodel itself
                if startmodel1 != startmodel:
                    d_startingmodel[startmodel1] += d_startingmodel[startmodel]
                    del d_startingmodel[startmodel]
                break
            ## current starting model (startmodel1) is a starting model of startmodel2
            ## trace back startmodel2
            else:
                startmodel1 = startmodel2


    d_startingmodel_reverse = {}
    for startingmodel,pdbs in d_startingmodel.items():
        for pdb in pdbs:
            d_startingmodel_reverse[pdb] = startingmodel


    return d_startingmodel, d_startingmodel_reverse


def group_pdbs(d,l_pdbs,):

##    ## sort by count
##    l = [[len(d[s]),s,] for s in d.keys()]
##    l.sort()
    ## sort alphabetically
    l = [[s,s,] for s in d.keys()]
    l.sort()

    l_pdbs_grouped = []
    for s in l:
        l_pdbs_group = d[s[1]]
        l_pdbs_group.sort()
        for pdb in l_pdbs_group:
            if not pdb in l_pdbs:
                continue
            l_pdbs_grouped += [pdb]

    return l_pdbs_grouped


def average_physprop(s,delimiter,):

    l = s.split(delimiter)
    if len(l) != 2:
        stop
    Sum = 0
    for s_val in l:
        Sum += float(s_val)
    average = Sum/len(l)
    
    return average


def find_min_rmsd(d_coordinates,):

##d_min_rmsd = {}
##for pdb_mutant in l_mutants:
##    coords_mutant = d_coordinates[pdb_mutant]
##    min_rmsd = [9.,None,]
##    l_rmsds = []
##    for pdb_wt in l_wts:
##        coords_wt = d_coordinates[pdb_wt]
##        rmsd = instance_geometry.superpose(coords_mutant,coords_wt,)
##        if rmsd < min_rmsd[0]:
##            min_rmsd = [rmsd,pdb_wt,]
##        if rmsd < 0.15:
##            l_rmsds += [pdb_wt]
##    d_min_rmsd[pdb_mutant] = min_rmsd
####    print pdb_mutant, min_rmsd
##    print pdb_mutant, l_rmsds
##print d_min_rmsd

    return

if __name__ == '__main__':
    main()

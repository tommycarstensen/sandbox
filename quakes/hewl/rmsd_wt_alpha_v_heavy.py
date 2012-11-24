import os
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc')
import statistics, gnuplot

##
## protein selection
##
protein = sys.argv[sys.argv.index('-protein')+1]

##
## parse dictionaries
##
fd = open('d_%s.txt' %(protein),'r')
s = fd.read()
fd.close()

d = eval(s)

d_startingmodel = d['startmodel']
d_spacegroups = d['spacegroup']
d_author = d['author']
d_entities = d['entities']
d_matthews = d['matthews']
d_citation = d['citation']
d_Z = d['Z']

## correction
d_author['1AKI_A'] = ['Carter, D.','He, J.','Ruble, J.R.','Wright, B.',]

fd = open('table_%s.txt' %(protein),'r')
lines = fd.readlines()
fd.close()
d_mutations = {}
d_ligands = {}
d_class = {}

d_table = {}
l_headers = [s.strip() for s in lines[0].split('\t')]
for s in l_headers:
    d_table[s] = {}

for line in lines[1:]:
    l = line.split('\t')
    pdb = l[0]
    for i in range(len(l_headers)):
        header = l_headers[i]
        s = l[i].strip()
        d_table[header][pdb] = s
    if protein == 'T4L':
        s_mutations = l[-1].strip()
        mutation_class = l[-2].strip()
    elif protein == 'HEWL':
        s_mutations = l[-2].strip()
        mutation_class = l[-2].strip()
    d_class[pdb] = mutation_class
    s_ligands = l[14].strip()
    d_mutations[pdb] = s_mutations
    s_ligands = s_ligands
    if s_ligands in ['XE','KR','AR',]:
        s_ligands = ''
    d_ligands[pdb] = s_ligands

l_pdbs_w_ligands = []
for pdb in d_ligands.keys():
    if d_ligands[pdb] != '':
        l_pdbs_w_ligands += [pdb]

##
## group chains by crystal contacts if multiple chains in ASU
##
d_asu = {
    ## P 21
    '1':[
        ],
    '2a':[
        ## IOD/SCN +
        '1B2K_A','1LCN_B','1LKR_B',
        ],
    '2b':[
        ## one chain (dehydrated)
        '2D4J_A','1XEI_A','1XEJ_A','1XEK_A','2Z12_A','2Z18_A','2Z19_A','1LMA_A',
        ## two chains, no iodide
        '3LYT_A','4LYT_A',
        '1LJ3_A','1LJ4_A','1LJE_A','1LJF_A','1LJG_A','1LJH_A','1LJI_A','1LJJ_A','1LJK_A','1UCO_A','1VDP_A','5LYM_A','1LYS_A','1HF4_A', ## A chain
        '2D4I_A','2D4K_A', ## A chain
        '1JJ3_A', ## A chain
        ## IOD/SCN -
        '1B2K_B','1LCN_A','1LKR_A',
        ],
    '2c':[
        ## two chains, no iodide
        '3LYT_B','4LYT_B',
        '1LJ3_B','1LJ4_B','1LJE_B','1LJF_B','1LJG_B','1LJH_B','1LJI_B','1LJJ_B','1LJK_B','1UCO_B','1VDP_B','5LYM_B','1LYS_B','1HF4_B', ## B chain
        '2D4I_B','2D4K_N', ## B chain
        '1JJ3_B', ## B chain
        ],
    ## P 1
    '3':[
        '1V7T_B',
        '1V7T_A',
        ],
    }

## HEWL mutants...
l_mutants = [
    '1FN5_A', '1LZD_A', '1IR8_A', '1UIF_A', '1LZE_A', '1UIC_A', '1HEQ_A', '1FLQ_A', '1IOR_A', '1HEO_A', '1HEP_A', '1LSM_A', '1FLW_A', '1LZG_A', '1FLY_A', '1IR9_A', '1HEN_A', '1LSN_A', '1H6M_A', '1HEM_A', '1IR7_A', '1UID_A', '1IOQ_A', '1KXY_A', '1UIE_A', '1KXX_A',
    '1FLU_A', '1IOT_A', '1HER_A', '1KXW_A', '3A3Q_A', '3A3R_X', '1IOS_A',
##    ## MODRES
##    '1AT5_A','1AT6_A','132L_A','1RCM_A',
    ]

##l_pressure = ['1C6%s_A' %(s) for s in range(0,10)+['A','B','C','D','E','F','G','H','I','J','K','L','M','N','P','Q','T',]]
##l_pressure += ['2B6T_A','2B6W_A','2B6X_A','2B6Y_A','2B6Z_A','2B70_A','2B72_A','2B73_A','2B74_A','2B75_A','2OE4_X','2OE7_X','2OE9_X','2OEA_X',]

##
## columns
##

d_colors = {'C0C0C0':'grey'}

if protein == 'HEWL':
    l_columns = [
        [ 2,'C0C0C0',7,1.4,'absence vs presence of iodide/thiocyanate ions (1B2K:A, 1LCN:B, 1LKR:B)',], ## grey
        ]
else:
    l_columns = [
        [ 2,'C0C0C0',7,1.0,'interm. disulf. bonds or intram. disulf. bond (excl. Cys9-Cys164)',], ## grey
        ]

l_columns += [
    [ 3,'00FF00',7,1.0,'diff. SG',], ## green
    ]

if protein == 'HEWL':
    l_columns += [
        [ 6,'FF0000',7,1.0,'low hydration (1V7T, 1XEI, 1XEJ, 1XEK, 2Z12, 2Z18, 2Z19, 2D4J, 1LMA)',], ## red
        [ 8,'000000',7,1.0,'same SG, two chains in ASU with diff. cryst. con.',], ## black
        [13,'808080',7,1.0,'NMR',], ## dark grey
        [14,'FF0080',7,1.0,'MD',], ## purple (between pink and red)
        ]
elif protein == 'T4L':
    l_columns += [
        [ 8,'000000',7,2.0,'same SG, two or more chains in ASU with diff. crystal contacts',], ## black
        [ 6,'FF0000',7,3.0,'same SG, A129W or T21H,T142H w ion or 3DKE',], ## red
##        [13,'808080',7,1.0,'same SG, apolar cavity L99A and polar ligand',], ## dark grey
##        [ 14,'000000',7,1,'180L (T26E mutant, P 1 21 1 space group) vs P 32 2 1 space group',],
##        [ 14,'80FF00',7,1,'173l (K16E,R119E,K135E,K147 reverse charge mutant)',], ## bright green (between yellow and green)
##        [ 15,'00FF80',7,1,'1l97 (I3P mutant)',], ## dark green (between green and cyan)
##        [ 14,'FF0080',7,1,'3DKE',], ## purple (between pink and red)
        ]

l_columns += [    
    ## starting models different
    [ 5,'0000FF',7,1.0,'same SG, diff. auth., diff. SM',], ## blue
    [10,'00FFFF',7,1.0,'same SG, same auth., diff. SM',], ## cyan

    ## starting models probably different because authors different
    [ 7,'FFFF00',7,1.0,'same SG, diff. auth., unkn. SM',], ## yellow

    ## starting models perhaps identical because authors identical
    [ 9,'FF8000',7,1.2,'same SG, same auth., unkn. SM, diff. prim. ref.',], ## orange
    [12,'0080FF',7,1.4,'same SG, same auth., unkn. SM, same prim. ref.',], ## light blue

    ## starting models identical
    [11,'FF00FF',7,1.4,'same SG, diff. auth., same SM',], ## pink
    [ 4,'8000FF',7,1.4,'same SG, same auth., same SM',], ## violet
    ]

##l_columns.reverse()

n_columns = len(l_columns)

##
## parse rmsds
##

suffix_exclusion = 'occTrue_altFalse_tempFalse_symmTrue'

##suffix_exclusion = 'occTrue_altFalse_tempTrue_symmTrue'

def main():

    if protein == 'HEWL':
        if 'tempFalse' in suffix_exclusion:
            xmax = 1.8; ymax = 40
            xmax = 3.8; ymax = 65 ## incl. MDs and NMR
        elif 'tempTrue' in suffix_exclusion:
            xmax = 1.8; ymax = 40
        title = 'alpha carbon RMSD and average chi1 difference between %s wt structures' %(protein)
    elif protein == 'T4L':
        if 'tempFalse' in suffix_exclusion:
            xmax = 4.2; ymax = 45
        elif 'tempTrue' in suffix_exclusion:
            xmax = 4.2; ymax = 40
        title = 'alpha carbon RMSD and average chi1 difference between %s wt, wt* and single point mutant structures' %(protein)

    fd = open('rmsds_%s_alpha_%s.txt' %(protein,suffix_exclusion),'r')
    lines_alpha = fd.readlines()
    fd.close()

    fd = open('rmsds_%s_heavy_%s.txt' %(protein,suffix_exclusion),'r')
    lines_heavy = fd.readlines()
    fd.close()

    fd = open('d_chi1_%s_%s.txt' %(protein,suffix_exclusion,),'r')
    s = fd.read()
    fd.close()
    d_chi1 = eval(s)

    if os.path.isfile('chi1_diff_avg_%s_%s.txt' %(protein,suffix_exclusion,)):
        fd = open('chi1_diff_avg_%s_%s.txt' %(protein,suffix_exclusion,),'r')
        lines = fd.readlines()
        fd.close()
        d_chi1_diff_avg = {}
        for line in lines:
            l = line.split()
            pdb1 = l[0]
            pdb2 = l[1]
            chi1_diff_avg = float(l[2])
            if not pdb1 in d_chi1_diff_avg.keys():
                d_chi1_diff_avg[pdb1] = {}
            d_chi1_diff_avg[pdb1][pdb2] = chi1_diff_avg
    else:
        d_chi1_diff_avg = {}

    ##
    ## loop
    ##

    ## gnuplot
    lines = []

    ## statistics
    d_statistics = {}
    for col in range(2,n_columns+3):
        d_statistics[col] = {'alpha':[],'heavy':[],'chi1':[],}
    d_statistics['sameSG'] = {'alpha':[],'heavy':[],'chi1':[],}
    d_statistics['diffSG'] = {'alpha':[],'heavy':[],'chi1':[],}

    ## statistics
    l_rmsds_alpha = []
    l_rmsds_heavy = []

    ## initiate loop
    n_lines = min(len(lines_alpha),len(lines_heavy))
    for i in range(n_lines):

        d = {}

        l_alpha = lines_alpha[i].strip().split()
        l_heavy = lines_heavy[i].strip().split()

        pdb1 = l_alpha[0]
        pdb2 = l_alpha[1]
        rmsd_alpha = float(l_alpha[2])
        rmsd_heavy = float(l_heavy[2])

        bool_MDNMR1 = False
        if pdb1[:2] == 'MD' or pdb1[:4] == '1E8L':
            bool_MDNMR1 = True
        bool_MDNMR2 = False
        if pdb2[:2] == 'MD' or pdb2[:4] == '1E8L':
            bool_MDNMR2 = True

        ## skip if excluded by previous script...
        if bool_MDNMR1 == False and pdb1 not in d_entities.keys():
            continue
        if bool_MDNMR2 == False and pdb2[:4] != '1E8L' and pdb2 not in d_entities.keys():
            continue

    ##    if i % 10000 == 0:
    ##        print i, pdb1, pdb2

        ## modified residues (don't accept mutations)
        if len( set([pdb1,pdb2,]) & set([
            ## HEWL
                ## MLY, DM0
                '132L_A',
                ## CCS
                '1RCM_A','1RCM_B',
                ## SNN,isoaspartate (ASP)
                '1AT5_A','1AT6_A',
                ## glycosylated, mutant, radiated
                '1H6M_A','1QIO_A','1UC0_A','2B5Z_A','2HTX_A','2HU1_A','3LYT_A','3LYT_B','4LYT_A','4LYT_B',
            ]) ) > 0:
            continue

        bool_calc_chi1_diff_avg = False
        if not pdb1 in d_chi1_diff_avg.keys():
            bool_calc_chi1_diff_avg = True
        elif not pdb2 in d_chi1_diff_avg[pdb1].keys():
            bool_calc_chi1_diff_avg = True
    ##    if pdb1 == '1QS5_A' and pdb2 == '1QS9_A':
    ##        bool_calc_chi1_diff_avg = True
        if bool_calc_chi1_diff_avg == True:
            l_chi1_diff = []
            for res_no in d_chi1.keys():

                ## skip due to one of the following reasons:
                ## 1) Gly, Ala, Val
                ## 2) missing
                ## 3) modres
                ## 4) zero occupancy
                ## 5) high tempfactor
                ## 6) altlocs
                if not pdb1 in d_chi1[res_no].keys():
                    continue
                if not pdb2 in d_chi1[res_no].keys():
                    continue

                chi1_1 = d_chi1[res_no][pdb1]
                chi1_2 = d_chi1[res_no][pdb2]
                chi1_diff = abs(d_chi1[res_no][pdb2]-d_chi1[res_no][pdb1])
                if chi1_diff > 180:
                    chi1_diff = 360-chi1_diff
                l_chi1_diff += [chi1_diff]
    ##            if pdb1 == '1QS5_A' and pdb2 == '1QS9_A':
    ##                print res_no, chi1_diff
##            if len(l_chi1_diff) == 0: ## tmp!!!
##                continue

##            print pdb1,pdb2
##            print d_chi1.keys()
##            print d_chi1[100].keys()
##            print 'chi1_diff_avg_%s_%s.txt' %(protein,suffix_exclusion,)

            chi1_diff_average = sum(l_chi1_diff)/len(l_chi1_diff)
    ##        if pdb1 == '1QS5_A' and pdb2 == '1QS9_A':
    ##            stop
            fd = open('chi1_diff_avg_%s_%s.txt' %(protein,suffix_exclusion,),'a')
            fd.write('%s %s %s\n' %(pdb1,pdb2,chi1_diff_average,))
            fd.close()
        else:
            chi1_diff_average = d_chi1_diff_avg[pdb1][pdb2]

        x = rmsd_alpha
    ##    y = rmsd_heavy
##        y_property = 'heavy'
        y = chi1_diff_average
        y_property = 'chi1'

        if x > xmax:
            print x,xmax,pdb1,pdb2
            stopx
        if y > ymax:
            print y,ymax,pdb1,pdb2
            stopy

    ##    ## radioactive damage
    ##    if len(set(l_alpha[1:]) & set([
    ####        '3LYT_A','3LYT_B', ## 120K
    ####        '4LYT_A','4LYT_B', ## 120K
    ##        '5LYT_A','6LYT_A', ## 298K
    ##        ])) > 0:
    ##        continue

    ##    ## modified residues (don't accept mutations)
    ##    if len( set([pdb1,pdb2,]) & set([
    ##        ## HEWL
    ##            ## MLY, DM0
    ##            '132L_A',
    ##            ## CCS
    ##            '1RCM_A','1RCM_B',
    ##            ## SNN,isoaspartate (ASP)
    ##            '1AT5_A','1AT6_A',
    ##        ]) ) > 0:
    ##        continue

        bool_disulfide_inter = False
        ## T4L intermolecular disulfide bond
        if len( set([pdb1,pdb2,]) & set([
            '2HUK_A', ## CYS131, C 1 2 1
            '2HUL_A', ## CYS44, C 1 2 1
            '2HUM_A','2HUM_B', ## CYS72, P 1 21 1
    ##        '3HWL_A', ## CYS68, CYS93, P 32 2 1
            ]) ) > 0: ## one structure has an intermolecular disulfide bond
            bool_disulfide_inter = True

        bool_disulfide_intra = False
        ## T4L intramolecular disulfide bond
        ## CYS21-CYS142 (hinge gap)
        if len( set([pdb1,pdb2,]) & set([
            ## CYS21-CYS142 (hinge gap)
            '3GUI_A','3GUJ_A','3GUK_A','3GUK_B','3GUL_A','3GUL_B','3GUM_A','3GUM_B','3GUN_A','3GUN_B','3GUO_A','3GUO_B','3GUP_A','3GUP_B','1KNI_A',
            ## CYS21-CYS142 (hinge gap), CYS3-CYS97 (hinge), CYS9-CYS164 (CTERM)
            '152L_A',
            ]) ) == 1:
            bool_disulfide_intra = True

        ## T4L intramolecular disulfide bond
        ## not CYS21-CYS142 (hinge gap)
        elif len( set([pdb1,pdb2,]) & set([
            ## CYS3-CYS97 (hinge)
            '172L_A', ## P 21 21 21
            ## CYS127-CYS154 (helices)
            '178L_A', ## P 4 21 2
    ##        ## CYS3-CYS97 (hinge), CYS9-CYS164
    ##        '167L_A','167L_B', ## P 21 21 21
    ##        ## CYS9-CYS164 (CTERM)
    ##        '1L35_A', ## P 32 2 1
            ]) ) == 1: ## one structure (not both) has an intramolecular disulfide bond
            bool_disulfide_intra = True

        ## V87M in combination w L84M,L91M,L99M,L118M,L121M,(V111M,L133M)
        bool_V87M = False
        if len( set([pdb1,pdb2,]) & set([
            '1LWK_A','1LWG_A','1LPY_A',
            ]) ) == 1:
            bool_V87M = True

        l = []
        for pdb in [pdb1,pdb2,]:
            if pdb[:4] in l_pdbs_w_ligands:
                l += [1]
            else:
                l += [0]
        d['ligands'] = l

##        ## more than 2 mutations
        ## more than 1 mutation
        if protein == 'T4L':
            if mutations1.count('+') > 1 or d_mutations[pdb2[:4]].count('+') > 1:
                continue
            if 'CORE' in mutations1 or 'CORE' in d_mutations[pdb2[:4]]:
                continue
            if 'disulf' in mutations1 or 'disulf' in d_mutations[pdb2[:4]]:
                continue

        if protein == 'T4L':
            for mut in [
                'L99A', ## apolar cavity
                'A129L','A129M','A129V','A129W', ## core
                'L121M','L121A', ## core
                'F153L','F153M','F153I','F153A','F153V', ## core
                'M102Q', ## polar cavity
                'M102A', ## apolar cavity
                ]:
                l = []
                for pdb in [pdb1,pdb2,]:
                    if mut in d_mutations[pdb[:4]]:
                        l += [1]
                    else:
                        l += [0]
                d[mut] = l

        ## mutants
        if protein == 'HEWL':
            bool_mutation = False
            if len( set([pdb1,pdb2,]) & set(l_mutants) ) > 0:
                bool_mutation = True
                continue

        ## other entities than HEWL/T4L
        if bool_MDNMR1 == False and d_entities[pdb1] > 1:
            continue
        if bool_MDNMR2 == False and d_entities[pdb2] > 1:
            continue

        if bool_MDNMR1 == True:
            mutations1 = 'wt'
            startmodel1 = None
            spacegroup1 = None
            l_authors1 = []
            citation1 = None
        else:
            mutations1 = d_mutations[pdb1[:4]]
            startmodel1 = d_startingmodel[pdb1]
            spacegroup1 = d_spacegroups[pdb1]
            l_authors1 = d_author[pdb1]
            citation1 = d_citation[pdb1[:4]]

        if bool_MDNMR2 == True:
            mutations2 = 'wt'
            startmodel2 = None
            spacegroup2 = None
            l_authors2 = []
            citation2 = None
        else:
            mutations2 = d_mutations[pdb2[:4]]
            startmodel2 = d_startingmodel[pdb2]
            spacegroup2 = d_spacegroups[pdb2]
            l_authors2 = d_author[pdb2]
            citation2 = d_citation[pdb1[:4]]


    ##    if not 'P 43 21 2' in [spacegroup1,spacegroup2,]:
    ##        continue

        if spacegroup1 == spacegroup2:
            if protein == 'HEWL':
                if spacegroup1 == 'P 1 21 1' and spacegroup2 == 'P 1 21 1':
                    if pdb1 in d_asu['1'] and pdb2 in d_asu['1']:
                        d['spacegroup'] = 'identical_P21'
                    elif pdb1 in d_asu['2a'] and pdb2 in d_asu['2a']:
                        d['spacegroup'] = 'identical_P21'
                    elif pdb1 in d_asu['2b'] and pdb2 in d_asu['2b']:
                        d['spacegroup'] = 'identical_P21'
                    elif pdb1 in d_asu['2c'] and pdb2 in d_asu['2c']:
                        d['spacegroup'] = 'identical_P21'
                    else:
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                elif spacegroup1 == 'P 1' and spacegroup2 == 'P 1':
                    if pdb1 in d_asu['3'] or pdb2 in d_asu['3']:
                        if pdb1 in d_asu['3'] and pdb2 in d_asu['3']:
                            d['spacegroup'] = 'identical'
                        else:
                            d['spacegroup'] = 'identical_different_crystalcontacts'
                    else:
                        d['spacegroup'] = 'identical'
                else:
                    d['spacegroup'] = 'identical'
            elif protein == 'T4L':

    ##            if spacegroup1 == spacegroup2 and spacegroup1 != 'P 32 2 1':
    ##                continue

                ## different number of chains in unit cell
                if d_Z[pdb1[:4]] != d_Z[pdb2[:4]]:
                    d['spacegroup'] = 'identical_different_crystalcontacts'
                    continue ## tmp!!!

                else:

                    ## dimer
                    if pdb1[:4] == '137L' and pdb2[:4] == '137L':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '167L' and pdb2[:4] == '167L': ## intermolecular disulfide bond
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '1L97' and pdb2[:4] == '1L97':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '1SSY' and pdb2[:4] == '1SSY':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '216L' and pdb2[:4] == '216L':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '174L' and pdb2[:4] == '174L':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '175L' and pdb2[:4] == '175L': ## R96A
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '1QTH' and pdb2[:4] == '1QTH': ## A98M
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1 == '3GUM_A' and pdb2 == '3GUN_B':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1 == '3GUM_B' and pdb2 == '3GUN_A':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    ## trimer
                    elif pdb1[:4] == '1PQK' and pdb2[:4] == '1PQK':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '2Q9E' and pdb2[:4] == '2Q9E':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    ## tetramer
                    elif pdb1[:4] == '150L' and pdb2[:4] in ['150L','1JQU','3CDO',]: ## M6I
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '1JQU' and pdb2[:4] in ['1JQU','3CDO',]: ## W158L
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '3CDO' and pdb2[:4] == '3CDO'and (pdb1[-1] != 'A' or pdb2[-1] != 'D'): ## R96V
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '3FI5' and pdb2[:4] == '3FI5':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    ## pentamer
                    elif pdb1[:4] == '168L' and pdb2[:4] == '168L':
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '169L' and pdb2[:4] == '169L':
                        d['spacegroup'] = 'identical_different_crystalcontacts'

                    ##
                    elif pdb1[:4] == '137L' and pdb2[:4] in ['176L','180L',]: ## S44F
                        d['spacegroup'] = 'identical_different_crystalcontacts'
                    elif pdb1[:4] == '176L' and pdb2[:4] in ['180L',]: ## X32-39X
                        d['spacegroup'] = 'identical_different_crystalcontacts'

                    elif pdb1[:4] == pdb2[:4] and rmsd_alpha > 0.5:
                        print pdb1, pdb2, round(rmsd_alpha,2)
                        stop

                    else:
                        d['spacegroup'] = 'identical'

        else:
            d['spacegroup'] = 'different'

        if startmodel1 == pdb2[:4] or startmodel2 == pdb1[:4]:
            d['startingmodel'] = 'identical'
        elif startmodel1 == None or startmodel2 == None:
            d['startingmodel'] = 'unknown'
    ##                if pdb1[:4] == startmodel2 or pdb2[:4] == startmodel1:
    ##                    bool_identical = True
    ##                else:
    ##                    bool_identical = None
        elif startmodel1 == startmodel2:
            d['startingmodel'] = 'identical'
        else:
            d['startingmodel'] = 'different'

        if len( set(l_authors1) & set(l_authors2) ) > 0:
    ##        if 'Matthews, B.W.' not in l_authors1 and len( set(l_authors1) & set(l_authors2) ) == 1:
    ##            print pdb1, pdb2
    ##            print l_authors1
    ##            print l_authors2
    ##            print set(l_authors1) & set(l_authors2)
    ##            stop
            d['author'] = 'identical'
        else:
            d['author'] = 'different'

        ## low hydration
        if protein == 'HEWL':
            if (
                (
                    len( set([
                        ## P 1
                        '1V7T_A',
                        '1V7T_B',
                        ## P 1 21 1
                        '1XEK_A', ## 9% RH
                        '2Z18_A', ## transformed
                        '1LMA_A',
                        ]) & set([pdb1,pdb2,]) ) > 0 
                    )
                or
                (
                    len( set([
                        ## P 1 21 1
                        '2D4J_A', ## transformed
                        '1XEI_A','1XEJ_A', ## 38% RH
                        '2Z12_A','2Z19_A', ## transformed
                        ]) & set([pdb1,pdb2,]) ) > 0
                    )
                ):
                bool_dehydrated = True
            else:
                bool_dehydrated = False
        else:
            bool_dehydrated = False

        ## 21-142 bridge
        bool_42_121 = False
        if len( set([pdb1,pdb2,]) & set([
    ##        '257L_A', ## T21H,T142H
            '258L_A', ## T21H,T142H Zinc 21-142 bridge
            '259L_A', ## T21H,T142H Cobolt 21-142 bridge
            '260L_A', ## T21H,T142H Nickel 21-142 bridge
            '1EPY_A', ## T21H,Q141H,T142H Cobolt 21-142 bridge
            ]) ) == 1:
            bool_42_121 = True

        ## ions
        if len(set([
            '1LKR_B', ## similar to 1lcn
            '1B2K_A', ## similar to 1lcn
            '1LCN_B', ## similar to 1hf4 if not high tempfactor atoms used
            ]) & set([pdb1,pdb2,])) == 1:
            bool_ions = True
        else:
            bool_ions = False

        ## citation
        if citation1 == None or citation2 == None:
            bool_identical_citation = False
        elif citation1 == citation2:
            bool_identical_citation = True
        else:
            bool_identical_citation = False
            
    ################################################################################

    ##    if bool_mutation == False:
    ##        continue

        ##
        ## pre space groups different
        ##

##        ## 180l (T26E mutant, P 1 21 1 space group) vs P 32 2 1 space group
##        ## high similarity
##        if '180L' in [pdb1[:4],pdb2[:4],] and 'P 32 2 1' in [spacegroup1,spacegroup2,]:
##            col = 14
####            continue

##        if '173L' in [pdb1[:4],pdb2[:4],]:
##            col = 14
##
##        elif '1L97' in [pdb1[:4],pdb2[:4],]:
##            col = 15
##
##        elif '1P7S' in [pdb1[:4],pdb2[:4],]:
##            col = 16

        if bool_disulfide_intra == True and '172L_A' not in [pdb1,pdb2]:
            print pdb1,pdb2
            stopstopstop

        if protein == 'HEWL' and (bool_MDNMR1 == True or bool_MDNMR2 == True):
            if pdb1[:2] == 'MD' or pdb2[:2] == 'MD':
##                if pdb1[:4] == '1E8L' or pdb2[:4] == '1E8L':
##                    continue
                col = 14
            else:
                col = 13

        elif (
            (protein == 'HEWL' and bool_ions == True)
            or
            (protein == 'T4L' and bool_disulfide_intra == True)
            or
            (protein == 'T4L' and bool_disulfide_inter == True)
            ):
            col = 2
            if protein == 'HEWL' and rmsd_alpha < 0.15: print 'ions', pdb1, pdb2, round(rmsd_alpha,2)
            elif protein == 'T4L' and rmsd_alpha < 0.4: print 'disulf', pdb1, pdb2, round(rmsd_alpha,2)

        ##
        ## space groups
        ##

        elif d['spacegroup'] == 'different':
            col = 3
            if protein == 'HEWL' and rmsd_alpha < 0.15: print '33333333', pdb1, pdb2, round(rmsd_alpha,2)
            elif protein == 'T4L' and rmsd_alpha < 0.30: print '33333333', pdb1, pdb2, round(rmsd_alpha,2)

        ##
        ## space groups not different
        ##

        elif protein == 'T4L' and len( set([pdb1,pdb2,]) & set([
    ##        '1L82_A', ## L99F+M102L+V111I+F153L
##            '1SSY_A','1SSY_B', ## G28A,I29A,G30A
##            '252L_A', ## Gly107/Gly110 phi/psi, M102A/M106A (not if altloc B)
            '3DKE_X', ## Gly107/Gly110 phi/psi, L99A/M102L, substrate; only outlier if high bfactor atoms included...
##            '192L_A', ## XXX40-48,127-132ALA
##            '1L75_A', ## D127A,E128A,V131A,N132A,L133A (why not 1l73 without L133A? and 1l74 without D127A)
##            '2NTH_A', ## Gly107/Gly110 phi/psi, spin label (L118[R1A]), steric clash with 
            ]) ) == 1:
            col = 6
            stoptmp
            if rmsd_alpha < 0.40 and '2NTH_A' not in [pdb1,pdb2,]:
##            if rmsd_alpha < 0.15:
                print pdb1,pdb2,round(rmsd_alpha,2)
                stop

##        elif protein == 'T4L' and len( set([pdb1,pdb2,]) & set([
####            '3HUK_A', ## Gly107/Gly110 phi/psi, L99A/M102Q, J0Z
####            '3HTG_A', ## Gly107/Gly110 phi/psi, L99A/M102Q, JZ7
####            '3HTF_A', ## Gly107/Gly110 phi/psi, L99A/M102Q, JZ6
####            '1QUD_A', ## Gly107/Gly110 phi/psi, L99G (similar to 3htf)
##
##    ##        '3HH6_A', ## L99A Gly107/Gly110 (use altloc B)
##            ]) ) == 1:
##            col = 14
##            if rmsd_alpha < 0.30:
####            if rmsd_alpha < 0.08:
##                print 'special', pdb1,pdb2,round(rmsd_alpha,2)
##                print mutations1, d_mutations[pdb2[:4]]
##                stop_remove_pdb_or_decrease_treshold

        ## His42-His121 bridge
        elif bool_42_121 == True:
            col = 6
            stoptmp
            if rmsd_alpha < 0.25:
                print pdb1,pdb2,rmsd_alpha
                stop

        ## mutations at core
        elif protein == 'T4L' and d['A129W'] in [[0,1],[1,0],]:
    ##        or
    ##        d['A129F'] != [0,0] ## 1QTC
    ##        ):
    ##        ## at core
    ##        (d_class[pdb1[:4]] == 'core' or d_class[pdb2[:4]] == 'core')
    ##        ## A129W
    ##        and
    ##        d['A129W'] != [0,0]
    ##        ## not a size switch mutation
    ##        and
    ##        ([d['L121A'][0],d['A129M'][0]] != [1,1] and [d['L121A'][1],d['A129M'][1]] != [1,1])
    ##        and
    ##        ([d['L121A'][0],d['A129V'][0]] != [1,1] and [d['L121A'][1],d['A129V'][1]] != [1,1])
    ##        ):
            col = 6
            stoptmp
##            if '2LZM_A' in [pdb1,pdb2,] and rmsd_alpha < 0.33:
            if rmsd_alpha < 0.40:
                print pdb1,pdb2,rmsd_alpha
                print d_class[pdb1[:4]], mutations1
                print d_class[pdb2[:4]], d_mutations[pdb2[:4]]
                print d['L121A']
                print d['A129M']
                print d['A129V']
                stop

    ##    ## polar cavity (L99A/M102Q)
    ##    elif (
    ##        protein == 'T4L' and d['L99A'] in [[1,0],[0,1]] and d['ligands'] == d['L99A'] and d['M102Q'] == d['L99A']#and d['M102Q'] in [d['L99A'],[1,1]]
    ##        and len( set([pdb1,pdb2,]) & set([
    ####            ## exclude
    ####            '1LGW_A', ## 1AN
    ####            '1LGX_A', ## 5AN
    ####            '1LI3_A', ## 3CH
    ####            '1XEP_A', ## CAQ
    ##            '2RBO_A', ## 265
    ##            '2RBS_A', ## 269
    ##            ]) ) == 1
    ##        ):
    ##        col = 13
    ##        if '2LZM_A' in [pdb1,pdb2,] and rmsd_alpha < 0.33:
    ##            print pdb1,pdb2,round(rmsd_alpha,2)
    ##            print d['L99A']
    ##            print d['M102Q']
    ##            stopsmall
    ##        if '2LZM_A' in [pdb1,pdb2,] and rmsd_alpha > 0.66:
    ##            print pdb1,pdb2,round(rmsd_alpha,2)
    ##            print d['L99A']
    ##            print d['M102Q']
    ##            stoplarge

##        elif (
##            protein == 'T4L' and d['L99A'] in [[1,0],[0,1]] and d['ligands'] == d['L99A'] and d['M102Q'] == [0,0]
##            and len( set([pdb1,pdb2,]) & set([
##                ## exclude apolar
##                '181L_A', ## BNZ (apolar)
##                '182L_A', ## BZF (apolar)
##                '183L_A', ## DEN
##                '184L_A', ## I4B
##                '185L_A', ## IND
##                '186L_A', ## N4B
##                '187L_A', ## PXY
##                '188L_A', ## OXE
##                '1L83_A', ## BNZ
##                '1L84_A', ## BNZ
##                '1NHB_A', ## PYL
##                '2OU0_X', ## MR3
##                '2OTZ_X', ## 1MR
##                '3DMX_A', ## BNZ
##                '3HH4_A', ## BNZ
##                ## neither polar/apolar
##                '2RAY_X', ## 258
##                '3HH3_A', ## B20
##                '3HH5_A', ## B24
##                ## polar
##    ##            '3DMZ_A', ## HFB
##    ##            '3DN0_A', ## F5B
##    ##            '3DN1_A', ## BCF
##    ##            '3DN2_A', ## BBF
##    ##            '3DN3_A', ## IBF
##    ##            '3DN4_A', ## PIH (3dna)
##    ##            '3DN6_A', ## F3B
##    ##            '3DN8_A', ## IBF (3dn3)
##    ##            '3DNA_A', ## PIH (3dn4)
##    ##            '2RAZ_X', ## 259
##    ##            '2RB0_X', ## 260
##    ##            '2RB1_X', ## 261
##    ##            '2RB2_X', ## 263
##    ##            '2OTY_X', ## YAN
##                ]) ) == 0
##            ):
##            col = 6 ## 13
##            print pdb1,pdb2
##            stop
##            if '2LZM_A' in [pdb1,pdb2,] and rmsd_alpha < 0.33:
##                print col,pdb1,pdb2,round(rmsd_alpha,2)
##                stop
##            if '2LZM_A' in [pdb1,pdb2,] and rmsd_alpha > 0.66:
##                print col,pdb1,pdb2,round(rmsd_alpha,2)
##                stop

        ## V87M in combination with L84-133M,V111M (1LWG,1LWK,1LPY)
        elif protein == 'T4L' and bool_V87M == True:
            col = 6
            stop
            if protein == 'T4L' and rmsd_alpha < 0.3:
                print pdb1,pdb2,round(rmsd_alpha,2)
                stop

        ## dehydrated
        elif protein == 'HEWL' and bool_dehydrated == True:
            col = 6
            if rmsd_alpha > 1.2 and 'tempFalse' in suffix_exclusion:
                print pdb1,pdb2
                stop

        elif d['spacegroup'] == 'identical_different_crystalcontacts':
            col = 8
            if protein == 'HEWL' and rmsd_alpha > 1.:
                print pdb1,pdb2
                print bool_ions
                stop
            if protein == 'T4L' and rmsd_alpha < 0.3:
                print '88888888', pdb1,pdb2,round(rmsd_alpha,2)

        ##
        ## space groups identical
        ##
                
        elif d['spacegroup'] in [
            'identical','identical_P21',
            ]:
            ## starting model different
            if d['author'] == 'different' and d['startingmodel'] == 'different':
                col = 5
                if chi1_diff_average < 2.5:
                    print 'diffauth,diffSM', pdb1, pdb2, 'small chi1', chi1_diff_average, mutations1, d_mutations[pdb2[:4]]
            elif d['author'] == 'identical' and d['startingmodel'] == 'different':
                col = 10
            ## starting model unknown
            elif d['author'] == 'different' and d['startingmodel'] == 'unknown':
                col = 7
                if rmsd_alpha < 0.05: print 'diffauth_unknSM', pdb1, pdb2, 'rmsd', round(rmsd_alpha,2)
            elif d['author'] == 'identical' and d['startingmodel'] == 'unknown':
                if bool_identical_citation == True:
                    col = 12
                    if chi1_diff_average > 20:
                        print 'sameauth_unknSM_sameCit', pdb1, pdb2, chi1_diff_average
    ##                if protein == 'HEWL' and rmsd_alpha > 0.5: print 'sameauth_unknSM', pdb1, pdb2, round(rmsd_alpha,2)
                    if protein == 'T4L' and rmsd_alpha > 0.85: print 'sameauth_unknSM', pdb1, pdb2, round(rmsd_alpha,2)
                elif bool_identical_citation == False:
                    col = 9
    ##                if protein == 'HEWL' and chi1_diff_average < 6:
    ##                    print 'sameauth_unknSM_diffCit', pdb1, pdb2, chi1_diff_average
    ##                if protein == 'HEWL' and rmsd_alpha < 0.1: print 'sameauth_unknSM_diffCit', pdb1, pdb2, spacegroup1, round(rmsd_alpha,2)
                    if protein == 'T4L' and rmsd_alpha > 0.85: print 'sameauth_unknSM_diffCit', pdb1, pdb2, round(rmsd_alpha,2)
            ## starting model identical
            elif d['author'] == 'different' and d['startingmodel'] == 'identical':
                col = 11
                if chi1_diff_average > 20:
                    print 'diffauth_sameSM', pdb1, pdb2, chi1_diff_average
            elif d['author'] == 'identical' and d['startingmodel'] == 'identical':
                col = 4
                if protein == 'HEWL' and rmsd_alpha > 0.6:
                    print '44444444', pdb1, pdb2, round(rmsd_alpha,2)
            else:
                print d
                stop
    ##        if d['spacegroup'] == 'identical_P21':
    ##            if col not in [4,]:
    ##                if rmsd_alpha < 0.15:
    ##                    print d
    ##                    print pdb1, pdb2
    ##                    stop
    ##                col = 8
        else:
            print pdb1,pdb2
            print d
            stop

    ##    if rmsd_alpha > rmsd_heavy:
    ##        print pdb1, pdb2, round(rmsd_alpha,2), rmsd_heavy

        ##
        ## outliers
        ##
        rmsd_alpha_max = {'HEWL':0.4,'T4L':0.5,}[protein]
        ## same SM
        if col in [
            11, ## auth diff, same SM
            4, ## auth same, same SM
##            13, 
            ] and rmsd_alpha > rmsd_alpha_max:
            print 'large', col, pdb1, pdb2, round(rmsd_alpha,2), spacegroup1, mutations1, d_mutations[pdb2[:4]],
            print d_class[pdb1[:4]], d_class[pdb2[:4]]
    ##        print d
        ## diff/unknown SM
        elif col in [
            7, ## auth diff, SM unkn
            5, ## auth diff, SM diff
            10, ## auth same, SM diff
            9, ## auth same, SM unkn, cit same
            12, ## auth same, SM unkn, cit diff
            ] and rmsd_alpha > rmsd_alpha_max+0.2:
            print 'large', col, pdb1, pdb2, round(rmsd_alpha,2), spacegroup1, mutations1, mutations2,
            print d_class[pdb1[:4]], d_class[pdb2[:4]]
        elif col in [
            2, # T4Ldisulf/HEWLions
            3, # space
            6, # special cases
            8, # crystal contacts
            13, # 
            15, # 
            14, # 
            16, # 
            ] and rmsd_alpha < 0.19:
##            ] and rmsd_alpha < 0.08:
            print 'small', col, pdb1, pdb2, round(rmsd_alpha,2), spacegroup1, mutations1, mutations2
    ##        print d
##            stop_small

    ##    if col in [4,7,9,12]:
    ##        continue

##        if col == 3 and rmsd_alpha > 1.5:
##            print pdb1, pdb2

        line = '%s %s%s\n' %(x, (col-2)*'N/A ', y,)
        lines += [line]

        d_statistics[col]['alpha'] += [rmsd_alpha]
        d_statistics[col]['heavy'] += [rmsd_heavy]
        d_statistics[col]['chi1'] += [chi1_diff_average]

        l_rmsds_alpha += [rmsd_alpha]
        l_rmsds_heavy += [rmsd_heavy]

        if col not in [2,6,]:
            if d['spacegroup'] == 'different':
                d_statistics['diffSG']['alpha'] += [rmsd_alpha]
                d_statistics['diffSG']['heavy'] += [rmsd_heavy]
                d_statistics['diffSG']['chi1'] += [chi1_diff_average]
            else:
                d_statistics['sameSG']['alpha'] += [rmsd_alpha]
                d_statistics['sameSG']['heavy'] += [rmsd_heavy]
                d_statistics['sameSG']['chi1'] += [chi1_diff_average]

    prefix = 'CA_v_%s_%s_%s' %(y_property,protein,suffix_exclusion,)

    fd = open('%s.gnuplotdata' %(prefix),'w')
    fd.writelines(lines)
    fd.close()

    average_alpha, stddev_alpha = statistics.do_stddev(l_rmsds_alpha)
    average_heavy, stddev_heavy = statistics.do_stddev(l_rmsds_heavy)
    print 'alpha rmsd', len(l_rmsds_alpha), 'average', average_alpha, 'stddev', stddev_alpha
    print 'heavy rmsd', len(l_rmsds_heavy), 'average', average_heavy, 'stddev', stddev_heavy

    ################################################################################

    #### mutants
    ##for i in range(n_columns):
    ##    l_columns += [[l_columns[i][0]+n_columns,'',l_columns[i][2],]]
    #### put mutants in the background behind wts
    ##l_columns = l_columns[n_columns:]+l_columns[:n_columns]

    l_colors = []
    l_pointsizes = []
    l_pointtypes = []
    for i in range(len(l_columns)):
        l_colors += [l_columns[i][1]]
        l_pointtypes += [l_columns[i][2]]
        if protein == 'HEWL':
            l_pointsizes += [1.5*l_columns[i][3]]
        else:
            l_pointsizes += [1.5*l_columns[i][3]]

    for i in range(len(l_columns)):
##        l_columns[i] = [l_columns[i][0],l_columns[i][-1]]
        l_columns[i] = [l_columns[i][0],'']

    gnuplot.scatter_plot_2d(
        prefix,
        xlabel = 'RMSD_C_{/Symbol a} / @^{/Symbol \ \260}A', ## {\305} is Angstrom if iso encoding
        ylabel = '<{/Symbol Dc_1}> / {/Symbol \260}',
        xmin = 0, ymin = 0,
        xmax = xmax,
        ymax = ymax,

    ##    ylabel = 'heavy atom RMSD',
    ##    ymax = ymax,
    ##    function = 'x',

        key_vert_pos = 'left',
##        title = 'alpha carbon and heavy atom RMSD between %s wt structures' %(protein),
        title = title,

        bool_multiple_columns = True,
    ##    d_columns = d_columns,
        l_columns = l_columns,
        l_colors = l_colors,
        l_pointtypes = l_pointtypes,
        l_pointsizes = l_pointsizes,

        pointsize = 1,

##        bool_title = False,

        bool_remove = False,
        )


    ##
    ## statistics
    ##
    print '\n-----statistics-----\n'
    l_sort = []
    d = {'alpha':{},'heavy':{},'chi1':{},}
    ##l_xtics = [l_columns[i][-1] for i in range(len(l_columns))]
    l_xtics = []
    l = []
    for col in [2,3,6,8,5,10,7,9,12,11,4,]:
        for i in range(len(l_columns)):
            if l_columns[i][0] != col:
                continue
            s = '%2i\t' %(col)
            for method in ['alpha','heavy','chi1',]:
                l_rmsds = d_statistics[col][method]
                n = len(l_rmsds)
                if method == 'alpha':
                    s += '%5i\t' %(n)
                if n in [0,1,]:
                    s += 'N/A\t'
                else:
                    average, stddev = statistics.do_stddev(l_rmsds)
                    if method == 'chi1':
                        s += '%5.2f +/- %5.2f\t' %(round(average,2), round(stddev,2),)
                    else:
                        s += '%4.2f +/- %4.2f\t' %(round(average,2), round(stddev,2),)
                    l_xtics += [l_columns[i][-1]]
                    d[method][l_columns[i][-1]] = l_rmsds
            s += '%s' %(l_columns[i][-1])
            if n > 0:
                l_sort += [[average,l_columns[i][-1],s]]
            l += ['%s\n' %(s)]
    fd = open('stats_%s.txt' %(protein),'w')
    fd.writelines(l)
    fd.close()

    l_xtics = []
    l_sort.sort()
    for average,s_col,s in l_sort:
        print s
    ##    l_xtics += [s.split('\t')[-1].strip()]
        l_xtics += [s_col]
    l_xtics.reverse()

    for method in ['alpha','heavy','chi1',]:
        prefix = 'rmsd_statistics_%s_%s' %(protein,method,)
        gnuplot.histogram(
            prefix, d[method], l_xtics,
            ylabel = '%s RMSD' %(method),
    ##        ymax = 5,
##            bool_remove = False,
            )

    ##
    ## statistics, excluding cross effects
    ##
    for prop,l_col1,l_col2, in [
        ['author',[5,7,11,],[10,9,12,4,],],
        ['model',[5,10,],[11,4,],],
        ['spacegroup',['sameSG',],['diffSG',],],
        ['A_spacegroup',[4,5,7,8,9,10,11,12,],[3,],],
        ['B_all_excl_special',range(3,6)+range(7,13),[],],
        ['C_all',range(2,13),[],],
        ]:
        for method in ['alpha','heavy','chi1',]:
            l_rmsds1 = []
            for col1 in l_col1:
                l_rmsds1 += d_statistics[col1][method]
            if len(l_col2) > 0:
                l_rmsds2 = []
                for col2 in l_col2:
                    l_rmsds2 += d_statistics[col2][method]
            else:
                l_rmsds2 = [0]
            print prop, method, len(l_rmsds1), round(sum(l_rmsds1)/len(l_rmsds1),2), len(l_rmsds2), round(sum(l_rmsds2)/len(l_rmsds2),2),
            print round( sum(l_rmsds1)/len(l_rmsds1)-sum(l_rmsds2)/len(l_rmsds2) , 2)

if __name__ == '__main__':
    main()

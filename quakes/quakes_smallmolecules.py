def main():

## identical
## ACE = ACY (same way as sugars...)
## NI = 3NI ## different oxidation states
## PEG = PG4 = PGE = 1PE = AE3
## 2xmonosacc = disacc
## HEM = HEC
## CLA = CL1
## MRD = MPD
## CIT = FLC
## 2HP = PI = PO4
## PTL = PE9
    
## different
## FE != FE2
## CU != CU1

    d = {}

    ## hetID:[chemical formula,charge]
    d_ions = {

        ## unknown
        'UNX':['',''], ## e.g. 1aqn
##
##        ## nitrate
        'NO3':['N1 O3',-1],
        ## ammonium
##        'NH4':['H4 N1',+1],
##        ## hydroxide
##        'OH' :['H1 O1',-1],
##        ## phosphate
##        '2HP':['O4 P1',-1],'PI' :['H1 O4 P1',-2],
        'PO4':['O4 P1',-3], ## different oxidation states
##        ## sulfate
        'SO4':['O4 S1',-2], ## different oxidation states
##        'SOH':[3,-1],
        ## carbonate
        'CO3':['C O3', -1],
        ## thiocyanate
        'SCN':['C N S', -1],
##        ## group1a
##        'LI' :['LI1',+1],
        'NA' :['NA1',+1],
        'K'  :['K1' ,+1],
##        'RB' :['RB1' ,+1],
##        'CS' :['CS1' ,+1],
##        ## group2a
##        'BEF':['BE F3',-1], ## BEF is a phosphate analog (e.g. 1w0j, 3b9b)
        'MG' :['MG1',+2],
        'CA' :['CA1',+2],
        'SR' :['SR' ,+2],
##        'BA' :['BA' ,+2],
##        ## group3a
##        'AL' :['AL1',+3], 'ALF':['AL F4',-1],'AF3':['AL F3', 0], ## ALF and AF3 are phosphate analogs (e.g. 3b9r), pH influences fluoride coordination number
##        'GA' :['GA1',+3],'GA1':['GA1',+2],
##        'TL' :['TL1',+1],
##        ## group4a?
##        'ARS':['AS1', 0],'ART':['O4 AS1',-3],'AST':-3,'TAS':['H3 O3 AS1', 0],
        'CAC':['C2 H6 AS O2', -1], ## different compounds
##        'PB' :['PB' ,+2],
##        ## group6a
##        'SE' :['SE1', 0],'SE4':['O4 SE1',-2], ## different compounds
##        ## group7a
        'CL' :['CL1',-1],
        'BR' :['BR1',-1],
        'IOD':['I1' ,-1], ## has an effect on HEWL structure 1b2k
##        ## group8a
##        'KR' :['KR1', 0],
        'XE' :['XE1', 0],
##        ## group3b
##        'V'  :+3,'VO4':['V1' ,-3], ## different oxidation states
##        ## group4b
##        'CR' :['CR1',+3],
##        'MO' :['MO1', 0],'4MO':['MO1',+4],'6MO':['MO1',+6],'2MO':['MO O2',-2], ## different compounds and different oxidation states
##        'WO4':['O4 W1',-2],
        ## group5b
        'MN' :['MN1',+2], ## different oxidation states
##        'MN3':['MN1',+3],
##        ## group6b
        'FE' :['FE1',+3], ## different oxidation states
        'FE2':['FE1',+2], ## different oxidation states
##        'FEO':['FE2 O1', 0],
##        
##        ## group7b
        'CO' :['CO1',+2],
##        '3CO':['CO1',+3], ## different oxidation states
##        ## group8b
        'NI' :['NI1',+2], ## different oxidation states
##        '3NI':['NI1',+3], 
##        'PD' :['PD1',+2],
##        'PT' :['PT1',+2],
##        ## group9b
        'CU' :['CU1',+2], ## different oxidation states
        'CU1':['CU1',+1], ## different oxidation states
##        ## group10b
        'ZN' :['ZN1',+2],
        'CD' :['CD1',+2],
        'HG' :['HG1',+2],
##        ## Lanthanides
##        'SM' :['SM1',+3],
##        'GD' :['GD1', 0],
##        'TB' :['TB1',+3],
##        'YB' :['YB1',+3],
##        ## Actinides
##        'IUM' :['O2 U',+4],
        }

    l_prosthetic_groups = [
        ## porphyrins (cyclic tetrapyrroles)
        ## Ferrochelatase catalyzes protophorphyrin+Fe(2+) --> protoheme + 2H(+)
        ## iron
        'HEM', ## protoporphyrin IX + Fe(II) (C3 vinyl,C8 vinyl,C18 methyl; tetramethyl,divinyl,dipropionate; *charge 2*)
        'HEC', ## Heme C (protoporphyrin IX; *charge 0*)
##        'HEA', ## Heme A (C3 hydroxyfarnesyl, C8 vinyl, C18 formyl)
##        'HEO', ## Heme O (C3 hydroxyfarnesyl, C8 vinyl, C18 methyl)
##        '2FH', ## 2-phenylheme
##        '1FH', ## 12-phenylheme
##        'DDH', ## dedivinyl,diacetyl heme (C3,C8)
##        'HEV', ## dedimethyl,divinyl heme
##        'HDM', ## tetramethyl,divinyl,dipropionate *ester* heme
##        'HAS', ## Heme-As (C3 geranylgeranyl, C8 vinyl, C18 formyl)
##        'VER', ## octaethylated porphyrin
##        'HEB', ## Heme B/C hybrid (tetramethyl,divinyl,dipropionate)
##        'DHE', ## Heme D (heme B derivative)
##        'HDD', ## cis-heme D hydroxychlorin gamma-spirolactone
##        'SRM', ## siroheme (partially reduced iron-porphyrin in e.g. nitrate reductase)
##        ## other metal
##        'HNI', ## protoporphyrin IX + Ni(II)
##        'HES', ## Zn substituted Heme C
        ## chlorophyll A
        'CLA',
        'CL1',
        'BCL',
        ]

    ## list of metals - also contains metalloids (e.g. Arsen)
    ## their connection is not checked (parse bond type from mmCIF instead)
    l_atoms_metal = [
        'LI','NA','K','RB','CS', ##1a
        'BE','MG','CA', ##2a
        'AL','GA','TL', ##3a
        'PB', ##4a?
        'AS', ##5a?
        'V', ##3b
        'CR','MO', ##4b
        'MN',     'RE', ##5b
        'FE', ##6b
        'CO', ##7b
        'NI','PT', ##8b
        'CU', ##9b
        'ZN','CD','HG', ##10b
        'SM','GD','YB', ## Lanthanides
        'U1', ## Actinides
        ]

    l_coenzymes = [
##        'RET', ## vitamin A
##        'TPP', ## vitamin B1
        'FMN', ## vitamin B2
        'FAD', ## vitamin B2
        'NAD', ## vitamin B3
        'NAP', ## vitamin B3, NADP
        'NDP', ## vitamin B3, NADPH
##        'NAJ','NAI',
        'COA', ## vitamin B5
        'PLP', ## vitamin B6
        'LLP', ## vitamin B6
##        'C2F', ## vitamin B9 (5-methyl THF)
        ]

    ## info from the PDB Ligand Depot (searched for cluster in chemical name)
    l_clusters = [
        ## iron clusters
        'SF4','SF3','FES',
        
##        'CFM','CFN','CLF','CLP','CN1','CNB','CNF','F3S','FS1','FS2','FSF','FSO','FSX','HC0','HC1','HF3','HF5','NFS','WCC','XCC',
##        ## iron-nickel clusters
##        'NFC',
##        ## copper clusters
##        'CUB','CUM','CUN','CUO','CUZ',
##        ## molybdenum "clusters"
##        'OMO',
##        ## hafnium clusters
##        'PHF',
##        ## zirconium clusters
##        'ZRC',
##        ## palladium clusters
##        'PLL',

        ## cobalt cluster (not really a cluster...)
        'NCO',
        ]

    ## wikipedia buffer solution, Good's buffers
    l_solutes = [

        ## water
        'HOH','DOD',

        ## unknown atom or ion
        'UNX',
        ## unknown ligand
        'UNL',
##
        ## ethylene glycol, protein precipitation
        'EDO',
        ## acetic acid (different protonation states)
        'ACT',
        'ACY',
        ## beta-mercapto-ethanol (reducing agent)
        'BME',
        ## poly-ethylene-glyocol (precipitant)
        'PEG', ## glycol (monomer Poly-EG)
        'PG4', ## Di-EG
        'PGE', ## Tri-EG
        '1PE', ## Penta-EG
##        'AE3',
        ## (amphiphile, surfactant)
        'C8E',
        ## Tris (buffering agent)
        'TRS',
        ## MES (buffering agent)
        'MES',
        ## citric acid (buffering agent, citrate synthase product...)
        'CIT', ## not charged
        'FLC', ## charged
        ## glycerol
        'GOL', ## glycerol (glycerol kinase substrate...)
        ## DMSO
        'DMS',
        ## isopropyl alcohol
        'IPA',
        ## HEPES (buffering agent)
        'EPE',
        ## imidazole
        'IMD', ## solvent

##        ## EEE, (not categorized...)
##        ##
##        ## methanol
##        'MOH',
        ## ethanol
        'EOH',
         ## 2-hydroxyethyl disulfide (solvent?)
        'HED',
        
##        ## di-thio-threitol (reducing agent)
##        'DTT', 
##        ## bis-tris methane (buffering agent)
##        'BTB',

        'FMT', ## solvent (wikipedia "Uses")
        'MPD', ## ("An overview on 2-methyl-2,4-pentanediol in crystallization and in crystals of biological macromolecules)
        'MRD', ## ("An overview on 2-methyl-2,4-pentanediol in crystallization and in crystals of biological macromolecules)

##
##        ## TAPS (buffering agent)
##        'T3A',
##        ## Bicine (buffering agent)
##        'BCN',
##        ## TES (buffering agent)
##        'NES',
##        ## MOPS
##        'MPO',
##        ## PIPES
##        'PIN',
##        ## Cacodylate (buffering agent)
##        'CAC',
##
        'BOG', ## membrane protein detergent (apolipoprotein substrate, 1l6l)
##        'SDS', ## detergent
        'LDA', ## surfactant (amphiphile)
##        'LI1', ## surfactant (amphiphile)
        'MYR', ## surfactant (amphiphile)
##        'SPM', ## spermine (lipophilic)
##        'DMF', ## solvent
##        'DTD', ## solvent? dithiane diol
##        'M2M', ## solvent (guess)
##        'BU3', ## EMBL Structural Bioinformatics Group
##        'PYR', ## pyruvate (pyruvate kinase substrate...)
##        'URE', ## denaturing agent (product of arginase...)
##
####        'POP', ## ion (why ignore???)

        'HGB', ## 4-(hydroxymercury)benzoate promotes growth of highly ordered crystals. How?
        ]

    ## only monosaccharides (NGA, FUC, GLC, NAG) involved in posttranslational modifications are included
    ## or sugars that are anomers and cannot be automatically distuinguished as it has been done during the remediation of the pdb archive
    d_saccharides = {

        ##
        ## monosaccharides, aldehydes
        ##

        ## pyranoses, hexoses
        'GLC':{'stereo':'GLC','derivate':['GLC']}, ## (alpha)-D-Glucose
        'BGC':{'stereo':'GLC','derivate':['GLC']}, ## beta-D-Glc

        'GAL':{'stereo':'GAL','derivate':['GAL']}, ## (beta)-D-Galactose
        'GLA':{'stereo':'GAL','derivate':['GAL']}, ## alpha-D-Gal

        'FUC':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, alpha-L-Fucose
        'FUL':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, beta-L-Fucose

        'MAN':{'stereo':'MAN','derivate':['MAN']}, ## alpha-D-Mannose
        'BMA':{'stereo':'MAN','derivate':['MAN']}, ## beta-D-Mannose

##        'ARA':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabinose
##        'ARB':{'stereo':'ARA','derivate':['ARA']}, ## beta-L-Arabinose

        ## pyranoses, pentoses
        'XYS':{'stereo':'XYS','derivate':['XYS']}, ## (alpha)-D-Xylose
        'XYP':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-Xylose
##        'LXC':{'stereo':'XYS','derivate':['XYS']}, ## beta-L-Xylose
##        'HSY':{'stereo':'XYS','derivate':['XYS']}, ## alpha-L-Xylose
##
####        ## furanoses, hexoses
####        'AHR':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabino*furano*se
####        ## furanoses, pentoses
####        'XYZ':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-xylo*furano*se
##
##        ## phosphorylated aldohexopyranoses
##        'G6P':{'stereo':'G6P','derivate':['GLC']}, ## alpha-D-Glc-6P
##        'BG6':{'stereo':'G6P','derivate':['GLC']}, ## beta-D-Glc-6P
##
####        ## deoxygenated aldohexopyranoses
####        'G6D':{'stereo':'G6D','derivate':['GLC']}, ## 6-deoxy-alpha-D-Glucose
####        ## oxygenated aldohexopyranoses
####        'KBG':{'stereo':'G6D','derivate':['GLC']}, ## 2-keto-beta-D-Glucose
##
        ## acetylated aldohexopyranose amines
        'NAG':{'stereo':'NAG','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine
        'NDG':{'stereo':'NAG','derivate':['GLC']}, ## 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
####        'NBG':{'stereo':'NBG','derivate':['GLC']}, ## 1-N-Acetyl-beta-D-Glucosamine
####        '16G':{'stereo':'16G','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine-6P
##        'NGA':{'stereo':'NGA','derivate':['GLC']}, ## (2)-N-Acetyl-D-Galactosamine
##
####        ##
####        ## monosaccharides, ketones
####        ##
####
####        ## furanoses, hexoses
####        'FRU':{'stereo':'FRU','derivate':['FRU']}, ## Fructose
####        'F6P':{'stereo':'F6P','derivate':['FRU']}, ## Fru-6P
##
        ##
        ## dissacharides
        ##
####        'SUC':{'stereo':'SUC','derivate':['GLC','FRU']}, ## GLC-a12-FRC, Sucrose
####        'LAT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, alpha-Lactose
####        'LBT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, beta-Lactose
####        'MAL':{'stereo':'MAL','derivate':['GLC','GLC']}, ## GLC-a14-GLC, Maltose
####        'TRE':{'stereo':'TRE','derivate':['GLC','GLC']}, ## GLC-a11a-GLC, Trehalose
####        'CBI':{'stereo':'CBI','derivate':['GLC','GLC']}, ## GLC-b14-GLC, Cellobiose
####        ##
####        ## polysaccharides
####        ##
####        'MTT':{'stereo':'MTT','derivate':['MAL','MAL']}, ## MAL-b14-MAL, Maltotetraose
####        ##
####        ## linear saccharides (neither furanoses nor pyranoses...) and their derivatives...
####        ##
####        'SOR':{'stereo':'SOR','derivate':['GLC']}, ## Sorbitol/Glucitol (reduced glucose)
####        'GLO':{'stereo':'GLO','derivate':['GLC']}, ## linear glucose
####        'XLS':{'stereo':'XLS','derivate':['XYL']}, ## linear xylose
####        'A5P':{'stereo':'ABF','derivate':['ARA']}, ## Arabinose-5-phosphate
####        'G6Q':{'stereo':'G6P','derivate':['GLC']}, ## Glc-6P (linear)
##
##        ##
##        ## conduritol (1,2,3,4-cyclohexenetetrol) derivatives
##        ##
##        'HMC':{'stereo':'HMC'}, ## 5-hydroxymethyl-chonduritol
##        'ACI':{'stereo':'ACI'}, ## 1-amino-5-hydroxymethyl-chonduritol

        ##
        ## Sialic acid (N-Acetylneuraminic acid, Neu5Ac, NANA)
        ##
        'SIA':{'stereo':'SIA','derivate':['SIA']}, ## (alpha)-sialic acid
        'SLB':{'stereo':'SIA','derivate':['SIA']}, ## beta-sialic acid
        }

    ## modified residues without parent aa/nucleotide residue
    l_peptidelinkers = [
##        ## l-peptide linkers (links nitrogens)
##        'DCL','FRD','2AO',
##        'BNO', ## 6lpr
##        'ADD', ## 2nyl
##        'GLM', ## 1khp
##        'DMT', ## 2rmb
##        'LOL', ## 1fkn
##        'ALQ', ## 1fkn
        ]

    ## functional groups (non-polymers)
    l_terminalmodres = [
        'ACE', ## N-terminal
        'NH2', ## C-terminal
        'PYR', ## N-terminal
##        'BOC','FMT', ## N-terminal
##        'PO0','PR0','PS0', ## N-terminal (after BOC...)
####        'OH',
##        'PHQ', ## 2az9,2azc
##        'CGN', ## 3cao
##        'KI2', ## 1u8g
##        'DIP', ## 1a1a (C-terminal)
##        'CH2', ## N- and C-terminal
##        'HOA', ## C-terminal
##        ## FOR,MYR,OHE,OME (N-term)
        ]

    ## d-peptide linkers
    l_dpeptidelinkers = [
        'DAL', ## ALA
        'DPN', ## PHE
        'DGN', ## GLN
        'ACB', ## ASP
        'FGA', ## NVA
        ]
    
    d['ions'] = d_ions
    d['prosthetic groups'] = l_prosthetic_groups
    d['metals'] = l_atoms_metal
    d['coenzymes'] = l_coenzymes
    d['clusters'] = l_clusters
    d['solutes'] = l_solutes
    d['saccharides'] = d_saccharides
    d['terminalmodres'] = l_terminalmodres
    d['dpeptidelinkers'] = l_dpeptidelinkers
    d['lpeptidelinkers'] = l_peptidelinkers

    return d

if __name__ == '__main__':
    main()

## most common hetIDs
## count by hetID
## [[299, 'CO3'], [300, 'SCN'], [302, 'LLP'], [302, 'SIA'], [307, 'CU1'], [310, 'MRD'], [314, 'MLE'], [314, 'OXY'], [315, 'C8E'], [320, 'CAC'], [332, 'GLA'], [337, ' XE'], [345, 'EOH'], [349, 'TRP'], [367, 'PGE'], [376, 'CSO'], [386, 'XYP'], [392, 'GTP'], [393, 'CMO'], [421, 'PTR'], [427, 'TPO'], [429, 'CGU'], [444, 'IMD'], [466, 'KCX'], [469, 'EPE'], [471, 'BCL'], [471, 'NCO'], [487, 'PCA'], [488, 'CME'], [491, 'COA'], [517, 'BOG'], [519, 'ANP'], [550, 'SEP'], [552, 'AMP'], [555, 'NDP'], [557, 'SAH'], [578, 'IPA'], [599, 'FES'], [622, 'NO3'], [627, 'PG4'], [650, 'FE2'], [672, 'GDP'], [676, 'CL1'], [677, 'CIT'], [688, 'HEC'], [699, 'MES'], [704, 'UNL'], [715, 'TRS'], [728, 'LDA'], [759, 'FUC'], [795, 'NH2'], [843, 'CLA'], [847, 'PEG'], [870, ' CO'], [899, 'ATP'], [900, 'FMN'], [955, 'BGC'], [955, 'GAL'], [956, 'NDG'], [1007, ' NI'], [1029, 'SF4'], [1059, 'PLP'], [1073, 'NAP'], [1074, 'ACE'], [1143, ' HG'], [1155, 'BME'], [1182, ' BR'], [1272, 'ACY'], [1365, 'BMA'], [1546, 'HYP'], [1574, 'MPD'], [1638, 'FMT'], [1648, 'GLC'], [1717, 'NAD'], [1748, ' CU'], [1761, 'UNX'], [1763, 'ADP'], [1827, 'FAD'], [2019, 'IOD'], [2132, ' FE'], [2289, 'DMS'], [2429, ' CD'], [3087, '  K'], [3146, 'MAN'], [3183, ' SR'], [3405, 'ACT'], [3685, 'MLY'], [4435, ' MN'], [4724, 'HEM'], [5641, 'PO4'], [9700, 'EDO'], [10867, ' CL'], [10966, ' NA'], [12012, 'NAG'], [13385, 'GOL'], [13935, ' CA'], [14122, ' ZN'], [25483, 'SO4'], [41799, ' MG'], [60262, 'MSE']]
## count by PDB
## [[108, 'LDA'], [108, 'SCN'], [109, 'PSU'], [110, 'HYP'], [111, 'MYR'], [116, 'XYP'], [118, 'F3S'], [119, 'FLC'], [119, 'POP'], [121, 'SUC'], [122, 'CU1'], [125, 'NCO'], [127, 'HED'], [128, 'OXY'], [129, 'CME'], [138, '1PE'], [138, 'MRD'], [140, 'SAM'], [141, 'CAC'], [141, 'GNP'], [145, 'UDP'], [146, 'LLP'], [147, ' BR'], [147, 'BOG'], [151, 'CSD'], [154, 'UNL'], [155, 'TYS'], [173, 'PGE'], [174, 'CO3'], [176, 'GTP'], [180, 'HEC'], [191, 'NO3'], [192, 'KCX'], [192, 'UNX'], [199, 'IPA'], [200, 'CMO'], [203, 'COA'], [204, 'DMS'], [205, 'CSO'], [210, 'IMD'], [213, 'IOD'], [215, 'PTR'], [236, 'EPE'], [242, 'TPO'], [243, 'AMP'], [243, 'NDP'], [252, 'ANP'], [258, 'FE2'], [278, 'FES'], [282, 'SEP'], [283, 'PG4'], [292, ' CO'], [311, 'BGC'], [322, 'PEG'], [324, ' HG'], [329, 'GAL'], [338, 'PCA'], [339, 'CIT'], [343, 'SAH'], [345, 'FMT'], [353, 'SF4'], [354, 'GLC'], [363, 'FUC'], [366, 'MES'], [385, 'NDG'], [386, 'TRS'], [396, 'GDP'], [399, ' NI'], [427, 'ACY'], [450, 'ATP'], [450, 'FMN'], [481, ' CD'], [490, 'MPD'], [519, 'NAP'], [546, 'PLP'], [555, 'NH2'], [565, 'BME'], [573, ' CU'], [578, 'BMA'], [632, 'NAD'], [648, 'ACE'], [695, ' FE'], [772, 'MAN'], [787, 'ADP'], [848, '  K'], [861, 'FAD'], [1270, ' MN'], [1271, 'ACT'], [1461, 'EDO'], [1992, 'PO4'], [2071, 'HEM'], [2147, 'NAG'], [2390, ' NA'], [3514, ' CL'], [3803, 'GOL'], [4509, 'MSE'], [4583, ' CA'], [4780, ' MG'], [5156, ' ZN'], [6815, 'SO4']]

## OBSOLETE
## SUL
## 'IPS' deprecated
## 'NAO','NA2','NA6','NA5','NAW' deprecated
## 'KO4' deprecated
## 'MO3','MO1','MO2','MO4','MO5','MO6' deprecated
## 'OC1','OC3','OC5','OC6','OC7','OC8','OC2','OC4' deprecated
## 'MW1','MW2','MW3','O4M','MN5','MN6' deprecated
## 'OF1','OF3','2OF' deprecated
## 'CO5','OCL','OCO','OCN','OCM' deprecated
## 'NI1','NI2','NI3','NIK' deprecated
## '1CU' deprecated
## 'ZN2','ZO3','ZN3','ZNO' deprecated

## MODRES for which MODRES is used
## MLY,HYP,CME,PCA,KCX,CGU,CSO,MLE,CSD,TYS
## SEP,TPO,PTR
## SAH (should be with A/ADP...)
## PSU

## frequently occuring hetIDs that have not been categorized
## ADP,GDP,UDP
## ATP,GTP
## AMP
## ANP,GNP
## POP
## CMO,OXY
## TRP
## SAM (should be with coenzymes...)

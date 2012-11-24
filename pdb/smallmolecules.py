def main():

    d = {}

    ## hetID:[chemical formula,charge]
    d_ions = {

        ## unknown
        'UNX':['',''], ## e.g. 1aqn

        ## nitrate, ammonium
        'NO3':['N1 O3',-1],'NH4':['H4 N1',+1],
        ## hydroxide
        'OH' :['H1 O1',-1],
        ## phosphate
        '2HP':['O4 P1',-1],'PI' :['H1 O4 P1',-2],'PO4':['O4 P1',-3], ## different oxidation states; 'IPS' deprecated
        ## sulfate
        'SO4':['O4 S1',-2],'SOH':[3,-1],'SUL':[3,-2], ## different oxidation states
        ## carbonate
        'CO3':['C O3', -1],
        ## group1a
        'LI' :['LI1',+1],
        'NA' :['NA1',+1], ## 'NAO','NA2','NA6','NA5','NAW' deprecated
        'K'  :['K1' ,+1], ## 'KO4' deprecated
        'RB' :['RB1' ,+1],
        'CS' :['CS1' ,+1],
        ## group2a
        'BEF':['BE F3',-1],
        'MG' :['MG1',+2], ## 'MO3','MO1','MO2','MO4','MO5','MO6' deprecated
        'CA' :['CA1',+2],'OC1':['CA1',+2], ## 'OC3','OC5','OC6','OC7','OC8','OC2','OC4' deprecated
        'SR' :['SR' ,+2],
        'BA' :['BA' ,+2],
        ## group3a
        'AL' :['AL1',+3],'ALF':['AL F4',-1],'AF3':['AL F3', 0],
        'GA' :['GA1',+3],'GA1':['GA1',+2],
        'TL' :['TL1',+1],
        ## group4a
        'ARS':['AS1', 0],'ART':['O4 AS1',-3],'AST':-3,'TAS':['H3 O3 AS1', 0],'CAC':['C2 H6 AS O2', -1], ## different compounds
        'PB' :['PB' ,+2],
        ## group6a
        'SE' :['SE1', 0],'SE4':['O4 SE1',-2], ## different compounds
        ## group7a
        'CL' :['CL1',-1],
        'BR' :['BR1',-1],
        'IOD':['I1' ,-1],
        ## group8a
        'KR' :['KR1', 0],
        ## group3b
        'V'  :+3,'VO4':['V1' ,-3], ## different oxidation states
        ## group4b
        'CR' :['CR1',+3],
        'MO' :['MO1', 0],'4MO':['MO1',+4],'6MO':['MO1',+6],'2MO':['MO O2',-2], ## different compounds and different oxidation states
        'WO4':['O4 W1',-2],
        ## group5b
        'MN' :['MN1',+2],'MN3':['MN1',+3], ## different oxidation states; 'MW1','MW2','MW3','O4M','MN5','MN6' deprecated
        ## group6b
        'FE2':['FE1',+2],'FE' :['FE1',+3],'FEO':['FE2 O1', 0], ## different oxidation states; 'OF1','OF3','2OF' deprecated
        ## group7b
        'CO' :['CO1',+2],'3CO':['CO1',+3], ## different oxidation states; 'CO5','OCL','OCO','OCN','OCM' deprecated
        ## group8b
        'NI' :['NI1',+2],'3NI':['NI1',+3], ## different oxidation states; 'NI1','NI2','NI3','NIK' deprecated
        'PD' :['PD1',+2],
        'PT' :['PT1',+2],
        ## group9b
        'CU1':['CU1',+1],'CU' :['CU1',+2], ## different oxidation states; '1CU' deprecated
        ## group10b
        'ZN' :['ZN1',+2], ## 'ZN2','ZO3','ZN3','ZNO' deprecated
        'CD' :['CD1',+2],
        'HG' :['HG1',+2],
        ## Lanthanides
        'SM' :['SM1',+3],
        'GD' :['GD1', 0],
        'TB' :['TB1',+3],
        'YB' :['YB1',+3],
        ## Actinides
        'IUM' :['O2 U',+4],
        }

    l_prosthetic_groups = [
        ## porphyrins (cyclic tetrapyrroles)
        ## Ferrochelatase catalyzes protophorphyrin+Fe(2+) --> protoheme + 2H(+)
        ## iron
        'HEM', ## protoporphyrin IX + Fe(II) (C3 vinyl,C8 vinyl,C18 methyl; tetramethyl,divinyl,dipropionate; *charge 2*)
        'HEC', ## Heme C (protoporphyrin IX; *charge 0*)
        'HEA', ## Heme A (C3 hydroxyfarnesyl, C8 vinyl, C18 formyl)
        'HEO', ## Heme O (C3 hydroxyfarnesyl, C8 vinyl, C18 methyl)
        '2FH', ## 2-phenylheme
        '1FH', ## 12-phenylheme
        'DDH', ## dedivinyl,diacetyl heme (C3,C8)
        'HEV', ## dedimethyl,divinyl heme
        'HDM', ## tetramethyl,divinyl,dipropionate *ester* heme
        'HAS', ## Heme-As (C3 geranylgeranyl, C8 vinyl, C18 formyl)
        'VER', ## octaethylated porphyrin
        'HEB', ## Heme B/C hybrid (tetramethyl,divinyl,dipropionate)
        'DHE', ## Heme D (heme B derivative)
        'HDD', ## cis-heme D hydroxychlorin gamma-spirolactone
        'SRM', ## siroheme (partially reduced iron-porphyrin in e.g. nitrate reductase)
        ## other metal
        'HNI', ## protoporphyrin IX + Ni(II)
        'HES', ## Zn substituted Heme C
        ]

    ## list of metals - also contains metalloids (e.g. Arsen)
    l_atoms_metal = [
        'LI','NA','K','CS', ##1a
        'BE','MG','CA', ##2a
        'AL','GA','TL', ##3a
        'PB', ##4a
        'AS', ##5a
        'V', ##3b
        'CR','MO', ##4b
        'MN',     'RE', ##5b
        'FE', ##6b
        'CO', ##7b
        'NI', ##8b
        'CU', ##9b
        'ZN','CD','HG', ##10b
        'GD', ## Lanthanides
        'U1', ## Actinides
        ]

    l_coenzymes = [
        ## vitamins
        'TPP','TDP', ## vitamin B1            
        'FMN','FAD', ## vitamin B2
        'NAD','NAI','NAJ','NAP','NDP', ## vitamin B3
        'COA', ## vitamin B5 / Coenzyme A
        'PLP', ## vitamin B6
        'BTN', ## vitamin B7 / biotin
##        'RET', ## vitamin A (not a coenzyme?)
        'C2F','THF','THL','DHF', ## vitamin B9 (5-methyl THF)
        'COB', ## vitamin B12 / Methylcobalamin / MeCbl
        'ASC', ## vitamin C / Ascorbic acid
        'PQN', ## vitamin K1 / Phylloquinone
        'F42', ## Coenzyme F420 (FMN derivative...)
        ## non-vitamins
        'U10', ## Coenzyme Q10
        'HEM',
        'MGD', ## Molybdopterin
        'SAH', ## SAM
        'H4B', ## Tetrahydrobiopterin
        ]

    ## info from the PDB Ligand Depot (searched for cluster in chemical name)
    l_clusters = [
        ## iron clusters
        'CFM','CFN','CLF','CLP','CN1','CNB','CNF','F3S','FES','FS1','FS2','FS4','FSF','FSO','FSX','HC0','HC1','HF3','HF5','NFS','SF3','SF4','WCC','XCC', ## 'FS3' deprecated ('F3S' maintained)
        ## iron-nickel clusters
        'NFC',
        ## copper clusters
        'CUB','CUM','CUN','CUO','CUZ',
        ## molybdenum "clusters"
        'OMO',
        ## hafnium clusters
        'PHF',
        ## zirconium clusters
        'ZRC',
        ## palladium clusters
        'PLL',
        ]

    ## wikipedia buffer solution, Good's buffers
    l_solutes = [

        ## unknown atom or ion
        'UNX',

        ## FMT,EEE,
        ##
        ## water
        'HOH','DOD',
        ## methanol
        'MOH',
        ## ethylene glycol, protein precipitation
        'EDO',
        ## acetic acid
        'ACY','ACT',
        ## citric acid (buffering agent, citrate synthase product...)
        'CIT', ## not charged
        'FLC', ## charged
        ## di-thio-threitol (reducing agent)
        'DTT', 
        ## beta-mercapto-ethanol (reducing agent)
        'BME',
        ## bis-tris methane (buffering agent)
        'BTB',

        ## TAPS (buffering agent)
        'T3A',
        ## Bicine (buffering agent)
        'BCN',
        ## Tris (buffering agent)
        'TRS',
        ## HEPES (buffering agent)
        'EPE',
        ## TES (buffering agent)
        'NES',
        ## MOPS
        'MPO',
        ## PIPES
        'PIN',
        ## Cacodylate (buffering agent)
        'CAC',
        ## MES
        'MES',

        ## Glycerol
        'GOL', ## glycerol (glycerol kinase substrate...)
        ## DMSO
        'DMS',
        ## isopropyl alcohol
        'IPA',

        'BOG', ## membrane protein detergent
        'LDA', ## surfactant (amphiphile)
        'LI1', ## surfactant (amphiphile)
        'MYR', ## surfactant (amphiphile)
        'SPM', ## spermine (lipophilic)
        'C8E',
        'FMT', ## solvent (wikipedia)
        'DMF', ## solvent
        'DTD', ## solvent? dithiane diol
        'IMD', ## solvent
        'M2M', ## solvent (guess)
        'BU3', ## EMBL Structural Bioinformatics Group
        'MPD', ## alcohol
        'MRD', ## diol
        'PEG', ## glycol (monomer of DEG and PEG)
        'PG4', ## DEG
        '1PE', ## penta-ethylene glycol
        'HED', ## alcohol, 2-hydroxyethyl disulfide
        'PYR', ## pyruvate (pyruvate kinase substrate...)
        'URE', ## denaturing agent (product of arginase...)

        'POP', ## ion
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

        'ARA':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabinose
        'ARB':{'stereo':'ARA','derivate':['ARA']}, ## beta-L-Arabinose

        ## pyranoses, pentoses
        'XYS':{'stereo':'XYS','derivate':['XYS']}, ## (alpha)-D-Xylose
        'XYP':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-Xylose
        'LXC':{'stereo':'XYS','derivate':['XYS']}, ## beta-L-Xylose
        'HSY':{'stereo':'XYS','derivate':['XYS']}, ## alpha-L-Xylose

##        ## furanoses, hexoses
##        'AHR':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabino*furano*se
##        ## furanoses, pentoses
##        'XYZ':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-xylo*furano*se

        ## phosphorylated aldohexopyranoses
##        'G1P':{'stereo':'G1P','derivate':['GLC']}, ## alpha-D-Glc-1P
        'G6P':{'stereo':'G6P','derivate':['GLC']}, ## alpha-D-Glc-6P
        'BG6':{'stereo':'G6P','derivate':['GLC']}, ## beta-D-Glc-6P
##        'BGP':{'stereo':'BGP','derivate':['GAL']}, ## beta-Gal-6P
##        'M1P':{'stereo':'M1P','derivate':['MAN']}, ## alpha-D-Man-1P
##        'M6P':{'stereo':'M6P','derivate':['MAN']}, ## alpha-D-Man-6P

##        ## deoxygenated aldohexopyranoses
##        'G6D':{'stereo':'G6D','derivate':['GLC']}, ## 6-deoxy-alpha-D-Glucose
##        ## oxygenated aldohexopyranoses
##        'KBG':{'stereo':'G6D','derivate':['GLC']}, ## 2-keto-beta-D-Glucose

        ## acetylated aldohexopyranose amines
        'NAG':{'stereo':'NAG','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine
        'NDG':{'stereo':'NAG','derivate':['GLC']}, ## 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
##        'NBG':{'stereo':'NBG','derivate':['GLC']}, ## 1-N-Acetyl-beta-D-Glucosamine
##        '16G':{'stereo':'16G','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine-6P
        'NGA':{'stereo':'NGA','derivate':['GLC']}, ## (2)-N-Acetyl-D-Galactosamine

##        ##
##        ## monosaccharides, ketones
##        ##
##
##        ## furanoses, hexoses
##        'FRU':{'stereo':'FRU','derivate':['FRU']}, ## Fructose
##        'F6P':{'stereo':'F6P','derivate':['FRU']}, ## Fru-6P

##        ##
##        ## dissacharides
##        ##
##        'SUC':{'stereo':'SUC','derivate':['GLC','FRU']}, ## GLC-a12-FRC, Sucrose
##        'LAT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, alpha-Lactose
##        'LBT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, beta-Lactose
##        'MAL':{'stereo':'MAL','derivate':['GLC','GLC']}, ## GLC-a14-GLC, Maltose
##        'TRE':{'stereo':'TRE','derivate':['GLC','GLC']}, ## GLC-a11a-GLC, Trehalose
##        'CBI':{'stereo':'CBI','derivate':['GLC','GLC']}, ## GLC-b14-GLC, Cellobiose
##        ##
##        ## polysaccharides
##        ##
##        'MTT':{'stereo':'MTT','derivate':['MAL','MAL']}, ## MAL-b14-MAL, Maltotetraose
##        ##
##        ## linear saccharides (neither furanoses nor pyranoses...) and their derivatives...
##        ##
##        'SOR':{'stereo':'SOR','derivate':['GLC']}, ## Sorbitol/Glucitol (reduced glucose)
##        'GLO':{'stereo':'GLO','derivate':['GLC']}, ## linear glucose
##        'XLS':{'stereo':'XLS','derivate':['XYL']}, ## linear xylose
##        'A5P':{'stereo':'ABF','derivate':['ARA']}, ## Arabinose-5-phosphate
##        'G6Q':{'stereo':'G6P','derivate':['GLC']}, ## Glc-6P (linear)

        ##
        ## conduritol (1,2,3,4-cyclohexenetetrol) derivatives
        ##
        'HMC':{'stereo':'HMC'}, ## 5-hydroxymethyl-chonduritol
        'ACI':{'stereo':'ACI'}, ## 1-amino-5-hydroxymethyl-chonduritol

        ##
        ## Sialic acid (N-Acetylneuraminic acid, Neu5Ac, NANA)
        ##
        'SIA':{'stereo':'SIA','derivate':['SIA']}, ## (alpha)-sialic acid
        'SLB':{'stereo':'SIA','derivate':['SIA']}, ## beta-sialic acid

        }

    ## functional groups
    l_functional = ['NH2','ACE','BOC','OH',]
    
    d['ions'] = d_ions
    d['prosthetic groups'] = l_prosthetic_groups
    d['metals'] = l_atoms_metal
    d['coenzymes'] = l_coenzymes
    d['clusters'] = l_clusters
    d['solutes'] = l_solutes
    d['saccharides'] = d_saccharides
    d['functional'] = l_functional

    return d

if __name__ == '__main__':
    main()

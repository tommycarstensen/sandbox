## todo: add new stereoisomers to quakes!!!

def main():

    import os, sys, urllib2

    path2 = '/local/data/pdb/'
    path1 = '/data/pdb/'

    pdbs = os.listdir(path1)
    pdbs.sort()

    errorpdbs = [
##        ## change from dnp (no missing atoms) to set (missing atom) for what reason?
##        '1bei','1jaa','1jc8','1jcp','1jd8',

        ## investigate
##        '1e0s', ## dtt > bme
##        '1kse', ## 5at
##        '1ek4', ## hyd > hoh (water)
##        '1lcs', ## 15p > toe (peg carbon/oxygen atom missing?)
##        '1lex','1ley', ## ipl > ilt (monomer to dimer?)
##        '1ent', ## pst > pta (incorrect hetID in nonr file?)
##        '1atn',

        ## tautomers (why change?)
##        '1dc4', ## G3H > G3P
        '1c7z', ## G3P > G3H
##        '1k8y', ## G3P > 13P
##        '1jld', ## PHM > PR0
    ##    ## carbonyl group to hydroxy group (move to tautomer section)
    ##    'DDB':'MDA','BSE':'ASA',

##        ## acid base (why change?)
##        '1lh6', ## 'a3d' > 'nio'
##
##        ## incorrect change from gdp to g25 because of missing atoms (beta-phosphate)
##        '1ije','1ijf',
##        ## incorrect change from nad to amp because of missing atoms
##        '1iso',
##        ## incorrect change from PGE to AE3 because of missing atom
##        '1j06',
##
##        ## rename from ips to po4, but should be to pi (q: why is pi in use?)
##        '1aal','1b93','1dam',
##
##        ## no het record for altloc
##        '1blc',
##
##        ## HETATM deletion
##        '1e3a',
##
##        ## standard residue in HET records (error or not?)
##        '1d40','1d41','1d76',
##        '1eg6', ## DU
##        '1lan', ## LEU
##
##        ## stereoinversion
##        ## dgl/glu
##        '1axd',
##
##        ## actinomycin (dth to thr)
##        '1a7y','1a7z',
##
##        ## incorrect d-amino acid change
##        '1jzp',
##
        ## check again (sugar conversion)
        '1dog',
##        '1au1',
####        '1ax0','1ax1','1ax2','1axy','1axz','1iew',
####        '1cap','1ioo','1iqq','1j86','1j8r','1dbg','1dbo','1dog','1dtu','1dwa',
####        '1dwf','1dwg',
##        ## nga to nag and man to rip (is latter correct?)
##        '1a7c',
##
##        ## oph > eno, leu > 2ac
##        '1a2c',
##
        ## o4 missing (g6d > gld)
        '1dbg','1dbo','1e3z',

        ## one and not multiple HET records
        '1cap',

        ]

    correctpdbs = [
    ##    ## methyl addition (from 5nc to 5cm)
    ##    '10mh',
    ##    ## renaming of units (small units to large unit) count atoms..?!
    ##    '1a2c','1a46',
    ##    ## renaming of units (large unit to small units) count atoms..?!
    ##    '1a85',
        ## correct change from 101 to da (due to incorrect bonding in non-rem)
        '1ii7',

        ## inccorect hetID in non-rem file
    ##    '106d', ## mc5>5cm
    ##    '1dak', ## dpp > dpu
    ##    '1kes', ## gls > g6s
    ##    '1lgn', ## apm > da
    ##    '1a46', ## rng (+bcy) > bic
    ##    '1a85', ## bnn > dbp
    ##    '1apm', ## pho > po3
        ## incorrect hetID in non-rem file (identical names and formulas)
        '1ad4','1aed',
        ## incorrect hetID in non-rem file (identical names and formulas except for hydrogen)
        '1at1','1bji','1bqy','1car',
        '1c82', ## change from gc4 (single bond) to dgc (double bond) for what reason? (hyaluronan?!)
        ## incorrect hetID in non-rem file (identical formulas but different names)
        '1ai6','1bp1','1atp','1a42','1apm','1cpc',
        '1d8e','1d8t',##'1lh1', ## 'ace' > 'act'
        ## incorrect hetID in non-rem file (identical formulas except for hydrogen but different names)
        '1e42', ## dtt > dtd (cyclization/oxidization)
        '1dr1','1dr3','1dr4','1dr6', ## change from bio (biopterin) to hbi (dihydro-biopterin) for what reason?
        ## incorrect hetID in non-rem file (more complex than 1to1 substitution)
        '1agm','1a46',
        ## incorrect hetID in non-rem file (different names and formulas)
        '1b6a','1ent',

        ## correct change from csx to cys (due to hydrolysis of differentiating atom)
        '1aj1',

        ## correct change from gld (dideoxy) to g6d (deoxy)
        '1a47',

        ## PO4 to SEP,TPO (SER,THR)
        '1cmk',
        ## PTR+TYR to PO4 (TYR)
        '1crx',
        ## CBM+GLU to CGA (GLU) in seqres?
        '1det',
        ## TMN+SER to TRS (SER) in seqres?
        '1jl0',
        ## PHS+HIS to NEP - in seqres??
        '1e58',

        ## hetID not mentioned in HET records in nonr file
        '1d21','1d22', ## AS
        '1e73', ## ZN

        ]


    ##
    ## multiple compounds to one compound
    ##
    d_changes = {
        ##
        ## nucleotide
        ##
        ## fluor
    ##    'FLO':['UFP','C49',],
        'F':['UFP','C49',],
        ## chlor
    ##    'CLO':['UCL',],
        ## brom
        'BR':['CBR','BRU',], ## '5BU',
        ## iod
##        'I':['C38','5IU',],
        'IOD':['5IU',],
        ## alkylation (not methylation)
        'ETH':['G36',],
        'SCC':['G47',],
        '3CN':['C31','U31',],
        ## methylation
        'C3M':['G31',],
        'C4M':['C34',],
        'C5M':['5CM',],
        'C6M':['6OG','6MA',],
        'CH3':['A2M','5CM','OMU','OMG','OMC','C49',],
        'CH2':['U34','G49','A40',],
        ## amide group
        'NH2':['A40','1AP',],
        ## methoxy
        'OME':['C45','A47',],
        ## "ethenylation"
        'ETD':['EDA',],
        ## oxygen
        'HYD':['64T','8OG',],
        'O'  :['8OG','A38',],
        ## sulfonation
        'SFO':['TYS',],
        ## phosphorylation (check 1atp,1j3v!!!)
        'PHS':['PTR',],
        ## s for o
        'S'  :['AS','SC',],
        ## other
        'DM6':['BNR',],
        'MEB':['BNR','BDA',],
        'DM1':['BDA',],
        ##
        ## peptide
        ##
        'GB':['PM3',], ## tyrosine
        ##
        ## water
        ##
        'HYD':['HOH',],
        }
    ##l_changes = []
    ##for k in d_changes.keys():
    ##    l_changes += d_changes[k]

    d_obsolete = {
        ## water
        'MTO':'HOH',
        'MGO':'MG','MO3':'MG','MO1':'MG','MO2':'MG','MO4':'MG','MO5':'MG','MO6':'MG',
        'OC1':'CA','OC2':'CA','OC5':'CA','OC6':'CA','OC7':'CA','OC3':'CA','OC4':'CA','543':'CA',
        'NAW':'NA',
        'MH3':'MN','MW2':'MN','O4M':'MN','MN5':'MN','MW1':'MN',
        'OF1':'FE','OF3':'FE','OF2':'FE',
        'ZNO':'ZN','ZN3':'ZN',
        'CD1':'CD','CD5':'CD','CD3':'CD',
        'NI1':'NI','NI2':'NI','NIK':'NI',
        'OCL':'CO','CO5':'CO','OCN':'CO',
        'KO4':'K',
        '1CU':'CU',
        ## standard residues
        'SEG':'SER','GGL':'GLU',
        ## phosphate double bond position
        '5IT':'5IU','ITS':'4IP',
        ## stereo-nonspecific to stereo-specific
        'PGY':'PG9',
        ## stereo-specific to stereo-nonspecific
        'GDB':'GDN',

        ## duplicates
        'SEO':'BME','GTO':'GCP','SBU':'NBU','5HP':'PCA','CYL':'ACI','MOL':'BIC',
        'FS3':'F3S','FS4':'SF4','NAH':'NAD','CBZ':'PHQ','HV5':'TBG','MGA':'MBG',
        'HAA':'DHY','2OG':'AKG','AMA':'PLA','OET':'ETH','BP' :'BAP','PGC':'PGH',
        'RON':'NVA','MHM':'HEM','DSP':'DAS','EHP':'MTY','HPG':'PDO','1PY':'PPY',
        'HPB':'PR0','NEM':'HIC','6AB':'BE2','IOH':'IPA','LLA':'HFA','ANE':'ADE',
        'SUL':'SO4','STY':'TYS','PNL':'CBG','2AS':'ACB','ICI':'ICT','EOX':'EOH',
        '2PI':'NVA','FA' :'FOL','CRY':'GOL','CEA':'CSO','AFL':'FUL','ISB':'ALQ',
        'NMO':'NO' ,'XL1':'SCC','PPM':'GB' ,'MGY':'SAR','TMB':'TBM','HGC':'MMC',
        'O2' :'OXY','EGL':'EDO','RIC':'RBZ','BH4':'H4B','HAC':'ALC','PAS':'PHD',
        'PLY':'PLM','OX' :'O'  ,'CBX':'FMT','PVL':'PYR','BUT':'NBU','THB':'H4B',
        'HMP':'HMI','BRI':'1GL','CL2':'CL1','COE':'MOT','E4N':'NET','KGR':'GLR',
        'ASQ':'PHD','CBM':'ACY','CN' :'CYN','DNJ':'NOJ','PTP':'THP','TMN':'TRS',
        'OLI':'OLA','450':'DMQ','MDG':'M7G','BDN':'BEN','NGL':'ASG','148':'BTB',
        'NAC':'A3D','MP3':'3PG','OTB':'BOC','OTE':'C8E','AA3':'ABA','7MQ':'MQ7',
        'BOX':'BEZ','HBL':'HBI','LIN':'AAE','BOT':'THZ','CPG':'5GP','GSA':'G4S',
        'GSD':'SGC','PAG':'2PG','2PL':'PGA','TIB':'TB9','C32':'CBR','UIC':'GRL',
        'DCG':'DGP','SAA':'APG','LP2':'LP1','PIG':'PGE','S4U':'4SU','NNS':'TYL',
        'NED':'IDG','SPG':'PGS','NTY':'GHP','2SI':'IDS','C25':'C5P','AGC':'GLC',
        'PGQ':'PGO','THQ':'TZP','TPH':'HPH','BFP':'FBP','HAP':'PLH','XPP':'HP2',
        'SI2':'SIA','AHA':'ACA','PA4':'IDG','GS4':'SGC','IDO':'IOD','13H':'243',
        'Q72':'672','SFN':'SO3','IPS':'PO4','LAU':'DAO','NIV':'NVP','OIP':'DI' ,
        'CMN':'ACI','CRC':'DKA','G42':'8OG','MP3':'3PG','FLO':'F'  ,'T37':'NYM',
        'PCC':'PCA','119':'P4P','BZO':'PHQ','MFA':'MFU','ZN2':'ZN' ,'ROB':'N'  ,
        'NWB':'CHH','T31':'U37','SCP':'DCS','DLA':'LAC','STO':'STU','NFB':'NFO',
        'DSX':'BAT','TRB':'TB9','GCM':'GLM','GCL':'XAO','TRG':'MLY','U25':'U5P',
        'OXO':'O'  ,'1NA':'MAG','DDM':'DMJ','DDB':'MDA','1NI':'LP1','IU6':'CHC',
        '0AT':'16G','INY':'CRP','0AV':'A2M','PC2':'PC2','MAM':'MMA','BTA':'NVA',
        'MPS':'VXA','BTC':'CYS','TS3':'GCG','GLI':'PSI','DOH':'BHD','GLB':'GAL',
        'OR5':'R5P','577':'SB2','CAY':'CCS','HNP':'H5P','DOX':'DIO','PEI':'LEA',
        'ACU':'ACE','NEV':'NVP','NEW':'PCQ','GLW':'G6D','SEU':'ITU','OLE':'1LU',
        '345':'CBP','FCY':'CYS','A39':'A2M','DGM':'ASO','SEA':'DHL','LTR':'TRP',
        'CLO':'CL' ,'DUM':'UNX','638':'XV6','GTN':'GNP','ADG':'TOA','PDL':'PP3',
        'T36':'5PY','0AU':'IU' ,'B1F':'B2F','OMB':'MOH','D1P':'ORP','4RR':'ROL',
        '5SA':'ARE','GUR':'GLL','BRO':'BR' ,'MEP':'T23','AC4':'AMZ','BPC':'BAP',
        'BAR':'TSA','395':'961','NC2':'NC1','RDP':'R1P','DIS':'HOH','CYP':'GPR',
        'HP3':'PGR','LOF':'HFA','SOM':'MPS','MBA':'DCI','FAT':'PLM','N2B':'PYL',
        'POC':'PC' ,'CYM':'SMC','DHO':'DXC','DHN':'AA4','G32':'6OG','ABK':'FKI',
        '2LU':'DCL','HPC':'PHM','DHU':'H2U','CNM':'ACM','U18':'F89','NAN':'SIA',
        'HHP':'PH2','TFH':'HDZ',

        }

    d_errors = {
    ##    ## incorrect hetID used OR maybe old hetID reused by new compound (should have own dic)
    ##    'BZO':'BZU','PTP':'HH2','2MA':'DTI','AGL':'AAL','APH':'4HP','PHO':'PO3',
    ##    'PAM':'PCT','ME' :'TPD','FUM':'FUG','TNP':'TN4',
    ##    'G21':'DPC','AD2':'ADE','PC' :'PC1','PGO':'2PE','AGS':'DGS',
    ##    'NMA':'CH2',
        }

    if len(set(d_errors.keys()) & set(d_obsolete.keys())) > 0:
        stop
    if len(set(d_errors.keys()) & set(d_changes.keys())) > 0:
        print set(d_errors.keys()) & set(d_changes.keys())
        stop
    if len(set(d_obsolete.keys()) & set(d_changes.keys())) > 0:
        print set(d_obsolete.keys()) & set(d_changes.keys())
        stop

    set_peptides = set(['SER','GLU',])

    d_sugars = {
        'XYS':'XYP','NAG':'NDG','MAN':'BMA','MPD':'MRD','FUC':'FUL',
        'NGA':'A2G','GLC':'BGC','SIA':'SLB','GAL':'GLA',
        }


    for i in range(len(pdbs)):
        pdb = pdbs[i][:4]
        if pdb in errorpdbs:
            continue
        if pdb in correctpdbs:
            continue
        if i < int(sys.argv[-1]):
            continue

        ## deprecated
        if not os.path.isfile('%s%s/pdb%s.ent' %(path2,pdb[1:3],pdb)):
            continue

        print
        print i,pdb,len(pdbs)
        fd = open('%s%s.pdb' %(path1,pdb),'r')
        lines1 = fd.readlines()
        fd.close()
        fd = open('%s%s/pdb%s.ent' %(path2,pdb[1:3],pdb),'r')
        lines2 = fd.readlines()
        fd.close()

        d_het = {}
        d_names = {}
        d_formulas = {}
        for j in range(2):
            d_het[j] = {}
            lines = [lines1,lines2][j]
            for line in lines:
                record = line[:6].strip()
                if record == 'HET':
                    hetID = line[7:10].strip()
                    if j == 0:
                        text = line[30:70].strip()
                        d_names[hetID] = text

                    ## check that old id not in new pdb
                    if j == 1 and hetID in d_obsolete.keys():
                        print pdb, hetID
                        stop1

##                    ## change obsolete hetID
##                    if j == 0 and hetID in d_obsolete.keys():
##                        hetID = d_obsolete[hetID]

                    ## add id to dic
                    if not hetID in d_het[j].keys():
                        d_het[j][hetID] = 0
                    d_het[j][hetID] += 1

                    print j, line,
                elif record == 'HETNAM':
                    hetID = line[11:14].strip()
                    continuation = line[8:10]
                    name = line[15:70].strip()
                    if continuation == '  ':
                        d_names[hetID] = name
                    else:
                        d_names[hetID] += ' '+name
                    print j, line,
                elif record == 'FORMUL':
                    hetID = line[12:15].strip()
                    formula = line[19:70].strip()
                    if formula[-1] == ')':
                        index1 = formula.index('(')+1
                        index2 = formula.index(')')
                        formula = formula[index1:index2]
##                        d_het[j][hetID] *= int(line[19:19+index1-1])
                    d_formulas[hetID] = formula
                    print j, line,


        print
        print i,pdb
        if d_het[0] != d_het[1]:


            ## check obsolete replacement
            d_het[2] = {}
            for hetID in d_het[0].keys():
                if hetID in d_obsolete.keys():
                    d_het[2][d_obsolete[hetID]] = d_het[0][hetID]
                else:
                    d_het[2][hetID] = d_het[0][hetID]
            if d_het[2] == d_het[1]:
                continue
            

            for hetID in d_het[0].keys():

                if hetID in [
##                    'SER','GLU',
##                    'DI',
##                    'HOH',
                    'IMP', ## 1bsu
                    ]:
                    continue

                if hetID not in d_het[1].keys():

                    ## check obsolete replacement
                    if hetID in d_obsolete.keys():
                        newID = d_obsolete[hetID]
                        if newID in set(d_het[1].keys())|set_peptides|set(['HOH']):
                            continue
                    else:
                        newID = hetID


                    print hetID


                    diff1 = set(d_het[0].keys()) - set(d_het[1].keys())
                    diff2 = set(d_het[1].keys()) - set(d_het[0].keys())

                    ## remove obsoletes
                    l_hets1 = list(diff1)
                    for het1 in l_hets1:
                        if het1 in d_obsolete.keys():
                            het2 = d_obsolete[het1]
                            print het1, het2, diff2
                            if het2 in d_het[1].keys():
                                diff1 -= set([het1])
                                diff2 -= set([het2])

                    ## sugars
                    sugars1 = set()
                    sugars2 = set()
                    for sugar1 in d_sugars:
                        sugar2 = d_sugars[sugar1]
                        sugars1 |= set([sugar1])
                        sugars2 |= set([sugar1])

                    ## amino acids
                    peptides1 = set([])
                    peptides2 = set([
                        'ASP',
                        'DSG',
                        ])

                    ## bases
                    nucleotides1 = set([])
                    nucleotides2 = set(['DU',])


                    ##
                    ## uniform stereconversion
                    ## 
                    if (
                        ## 1 stereocenter
                        (diff1-peptides1 == set(['XYS']) and diff2-peptides2 == set(['XYP'])) or
                        (diff1-peptides1 == set(['NAG']) and diff2-peptides2 == set(['NDG'])) or
                        (diff1-peptides1 == set(['MAN']) and diff2-peptides2 == set(['BMA'])) or
                        (diff1-peptides1 == set(['MPD']) and diff2-peptides2 == set(['MRD'])) or
                        (diff1-peptides1 == set(['FUC']) and diff2-peptides2 == set(['FUL'])) or
                        (diff1-peptides1 == set(['NGA']) and diff2-peptides2 == set(['A2G'])) or
                        (diff1-peptides1 == set(['GLC']) and diff2-peptides2 == set(['BGC'])) or
                        (diff1-peptides1 == set(['SIA']) and diff2-peptides2 == set(['SLB'])) or
                        (diff1-peptides1 == set(['GAL']) and diff2-peptides2 == set(['GLA'])) or
                        (diff1-peptides1 == set(['1AR']) and diff2-peptides2 == set(['ERI'])) or
                        ## 2 stereocenters
                        (diff1-peptides1 == set(['TAR']) and diff2-peptides2 == set(['TLA'])) or
                        (diff1-peptides1 == set(['YG' ]) and diff2-peptides2 == set(['YYG'])) or
                        ## other
                        (diff1-peptides1 == set(['PI' ]) and diff2-peptides2 == set(['PO4']))
                        ):
                        s1 = ''.join(list(diff1))
                        s2 = ''.join(list(diff2))
                        fd = open('hetID_stereo.txt','a')
                        fd.write('%4s %3s %3s\n' %(pdb,s1,s2))
                        fd.close()
                        continue

                    ##
                    ## mixed stereconversion
                    ## 
                    if hetID in d_sugars.keys():
                        het1 = hetID
                        het2 = d_sugars[het1]
                        count = 0
                        if het1 in d_het[1].keys():
                            count += d_het[1][het1]
                        if het2 in d_het[1].keys():
                            count += d_het[1][het2]
                        if d_het[0][het1] != count:
                            print het1,het2
                            print d_het[0]
                            print d_het[1]
                            stop_sugar
                        if len(diff1-sugars1) == 0:
                            print d_het[0]
                            print d_het[1]
                            continue


                    ## 1) hetID not obsolete and not changed
                    if not hetID in d_obsolete.keys()+d_changes.keys()+d_errors.keys():
                        print i, pdb
                        if (
                            len(diff1-set([
                                'XYS','SIA','FUC','MAN',
                                ])) == 0 and
                            len(diff2-set([
                                'XYP','SLB','FUL','BMA','NDG',
                                ])) == 0
                            ):
                            fd = open('hetID_multiplesugars.txt','a')
                            fd.write('%4s %s %s\n' %(pdb,str(diff1),str(diff2)))
                            fd.close()
                            continue
                        url = 'http://remediation.wwpdb.org/pyapps/IdHandler.py?query=%s' %(hetID)
                        urllines = urllib2.urlopen(url)
                        lines = urllines.readlines()
                        for j in range(len(lines)):
                            if '<b>Name</b>' in lines[j]:
                                index2 = lines[j+1].index('</td>')
                                index1 = lines[j+1].index('>')+1
                                name = lines[j+1][index1:index2]
                                print hetID, 'name', name
                            if 'Release status' in lines[j]:
                                index2 = lines[j+1].index('</td>')
                                index1 = lines[j+1].index('>')+1
                                status = lines[j+1][index1:index2]
                                print hetID, 'status', status
                            if 'Replaced by' in lines[j]: ## j+5
                                index2 = lines[j+1].index('</td>')
                                index1 = lines[j+1].index('>')+1
                                replaceID = lines[j+1][index1:index2]
                                print hetID, 'replaced by', replaceID
                                stop_obsolete
                        print d_het[0]
                        print d_het[1]
                        print 'old', diff1
                        print 'new', diff2

                        identical = compare_name_formula(
                            pdb,
                            diff1-sugars1-peptides1-nucleotides1,
                            diff2-sugars2-peptides2-nucleotides2,
                            d_names,d_formulas,
                            )
                        if identical == True:
                            continue

                        stop2_change_or_stereoisomer_or_wrongname

                    ## 2) ID obsolete but not in d_changes?
                    elif hetID in d_obsolete.keys() and hetID not in d_changes.keys() and newID not in d_changes.keys():

##                        ##
##                        ## obsolete ID
##                        ##
##                        newID = d_obsolete[hetID]
##                        if newID in d_het[1].keys():
##                            continue

                        print 'old', diff1
                        print 'new', diff2
                        identical = compare_name_formula(
                            pdb,diff1,diff2,d_names,d_formulas,
                            )
                        stop_new_change_or_incorrect_nonremID
                        ## which hetIDs in d_changes correspond to the previous hetID?
                        if prevID in d_changes.keys():
                            hetIDs = d_changes[prevID]
                        elif currentID in d_changes.keys():
                            hetIDs = d_changes[currentID]
                        else:
                            stop_expected_prev_or_current_ID_in_d_changes
                        ## prev changed ID correspond to current changed ID
                        if not (
                            len(set(d_het[1].keys())-set(hetIDs)) == 0
                            ):
                            print hetIDs
                            print d_het[1]
                            print hetID
                            stop3

                    ## 3) ID obsolete and in d_changes
                    elif hetID in d_changes.keys() or newID in d_changes.keys():
                        if newID in d_changes.keys():
                            currentID = newID
                        else:
                            stop
                            currentID = hetID
                        print d_het[0]
                        print d_het[1]
                        diff = set(d_het[1].keys()) - set(d_het[0].keys()) - set(d_changes[currentID]) - set([currentID]) - nucleotides2
                        for het in d_het[0].keys():
                            if het in d_obsolete.keys():
                                print diff
                                diff -= set([d_obsolete[het]])
                        for het in d_het[0].keys():
                            if het in d_changes.keys():
                                diff -= set(d_changes[het])
                        if len(diff) > 0:
                            print 'old', set(d_het[0].keys())-set(d_het[1].keys())
                            print diff
                            stop4_stdres_or_newchange
                        ## extra check
                        if not len(set(d_changes[currentID]) & set(d_het[1].keys()+['HOH'])) > 0:
                            print d_changes[currentID], d_het[1]
                            stop5
                    else:
                        print hetID
                        stop_end


def compare_name_formula(
    pdb,diff1,diff2,
    d_names,d_formulas,
    ):

    identical = False

    if len(diff1) == 1 and len(diff2) == 1:
        het1 = list(diff1)[0]
        het2 = list(diff2)[0]
        formula1 = d_formulas[het1]
        formula2 = d_formulas[het2]
        d_formula1 = formulastring2dic(formula1)
        d_formula2 = formulastring2dic(formula2)
        name1 = d_names[het1].replace('-','').replace(' ','')
        name2 = d_names[het2].replace('-','').replace(' ','')
        print het1, name1
        print het2, name2
        print d_formula1
        print d_formula2
        if name1 == name2 and d_formula1 == d_formula2:
            identical = True
            fd = open('hetID_incorrect_nonremediated.txt','a')
            fd.write('%4s %3s %3s\n' %(pdb,het1,het2))
            fd.close()

    else:
        for het in diff1 | diff2:
            if het in d_names.keys():
                print het, d_names[het]
            if het in d_names.keys():
                print het, d_formulas[het]

    return identical


def formulastring2dic(formula):
    
    d_formula = {}
    for x in formula.split():
        if x[-1] in ['+','-']:
            continue
        if len(x) == 1 or x[-1] in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            d_formula[x] = 1
        else:
            for i in range(len(x)):
                if x[i] not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    break
            print x, x[:i], x[i:]
            d_formula[x[:i]] = int(x[i:])

    return d_formula


if __name__ == '__main__':
    main()

##d_obs = {}
##s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890'
##for i in range(len(s)):
##    for j in range(len(s)):
##        for k in range(len(s)):
##            hetID = s[i]+s[j]+s[k]
##            print hetID
##            if hetID in d_obsolete.keys():
##                continue
##            url = 'http://remediation.wwpdb.org/pyapps/IdHandler.py?query=%s' %(hetID)
##            urllines = urllib2.urlopen(url)
##            lines = urllines.readlines()
##            for l in range(len(lines)):
##                if 'Release status' in lines[l]:
##                    index2 = lines[l+1].index('</td>')
##                    index1 = lines[l+1].index('>')+1
##                    status = lines[l+1][index1:index2]
##                if 'Replaced by' in lines[l]: ## j+5
##                    if status != 'OBS':
##                        stop
##                    index2 = lines[l+1].index('</td>')
##                    index1 = lines[l+1].index('>')+1
##                    replaceID = lines[l+1][index1:index2]
##                    print hetID, 'replaced by', replaceID
##                    d_obs[hetID] = replaceID
##print d_obs
##stop

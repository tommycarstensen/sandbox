## active sites
## CAA, ACC
## COA,COB are coenzymes

import sys

sublim = sys.argv[1]

import os
pdbpath = '/oxygenase_local/data/pdb/'
l_site_descriptions = [
    ## CATALYTIC TRIAD
    'CATALYTIC TRIAD',
    'TYR-LYS-SER CATALYTIC TRIAD',
    'CATALYTIC TRIAD ACTIVE SITE',
    'ACTIVE SITE, CATALYTIC TRIAD',
    'ACTIVE SITE CATALYTIC TRIAD',
    ## CATALYTIC
    'CATALYTIC RESIDUES',
    'CATALYTIC HISTIDINE',
    'CATALYTIC CYSTEINE',
    'CATALYTIC SITE', ## ASA
    'REFERS TO CATALYTIC SITE',
    'PUTATIVE CATALYTIC SITE', ## CTA
    'CATALYTIC RESIDUES OF THE ACTIVE SITE',
    'CATALYTIC ACID/BASE',
    'CATALYTIC CENTER',
    'CATALYTIC NUCLEOPHILE',
    'CATALYTIC BASE', ## ACB
    'CATALYTIC ACID', ## ACD
    'RESIDUE IMPLICATED IN CATALYSIS',
    'ACTIVE SITE CATALYTIC RESIDUES',
    'ACTIVE SITE',
    'ACTIVE-SITE',
    'THIS IS THE PROPOSED ACTIVE SITE FOR THE ENZYME', ## CAA
    'RESIDUE LINE THE ACTIVE SITE OF THE ENZYME',
    'ENZYME ACTIVE SITE',
    'PLP-BINDING LYSINE',
    'NUCLEOPHILE',
    'RESIDUES INVOLVED IN CATALYSIS',
    'RESIDUES INVOLVED IN CATALYSIS SITE', ## ASA
    'ACT BINDING SITE', ## NUC
    'ACID-BASE RESIDUE', ## NUC
    ]

set_site_names = set([])
d_site_names = {}
d_site_descriptions = {
    'CAT':'CATALYTIC SITE',
    'ACT':'ACTIVE SITE',
    'AC1':'CATALYTIC SITE',
    'AC2':'ACTIVE SITE', ## 3r1r
    'AC3':'ACTIVE SITE', ## 3r1r
    'ACI':'CATALYTIC NUCLEOPHILE', ## 1e0w
    'CA1':'CATALYTIC SITE',
    'CA2':'CATALYTIC CYSTEINE',
    'CA3':'CATALYTIC CYSTEINE',
    'CA4':'CATALYTIC CYSTEINE',
    'CA5':'ACTIVE SITE', ## 1h8v
    'CA6':'ACTIVE SITE', ## 1h8v
    'ACA':'CATALYTIC RESIDUES',
    'ACB':'CATALYTIC BASE',
    'ACC':'ACTIVE SITE', ## 6r1r
    'ACD':'CATALYTIC ACID',
    'ACE':'CATALYTIC RESIDUES', ## 8ruc
    'ACG':'CATALYTIC RESIDUES', ## 2dub,8ruc
    'CTA':'CATALYTIC TRIAD', ## 2ovw,1xzk
    'CTB':'CATALYTIC TRIAD', ## 1xzk
    'CTC':'CATALYTIC SITE', ## 1m6x
    'CTD':'CATALYTIC SITE', ## 1m6x
    'ASB':'ACTIVE SITE', ## 1a6e,1a0g
    'CAA':'CATALYTIC SITE', ## 5bca,1qm5
    'CAB':'CATALYTIC SITE', ## 1ajo
    'CAC':'CATALYTIC SITE', ## 1am7
    'CAD':'CATALYTIC HISTIDINE', ## 5bca,1auw
##    'S1':'CATALYTIC SITE',
##    'S4':'CATALYTIC SITE', ## 1kev
    'AS1':'ACTIVE SITE CATALYTIC TRIAD', ## 1qoz
    'AS2':'ACTIVE SITE CATALYTIC TRIAD', ## 1qoz
    'AS3':'ACTIVE SITE', ## 1hg1, chain C
    'AS4':'ACTIVE SITE', ## 1hg1, chain D
    'AS5':'ACTIVE SITE',
    'ASA':'ACTIVE SITE', ## 1a0g, subunit A
    'FEA':'ACTIVE SITE', ## 1bou
    'FEB':'ACTIVE SITE', ## 1bou
    'SA':'ACTIVE SITE', ## 1ips, chain A
    'SB':'ACTIVE SITE', ## 1j5s
    'GL1':'ACTIVE SITE', ## 1kas
    'ZN':'CATALYTIC CENTER',
    'ZNB':'ACTIVE SITE',
    'MNA':'CATALYTIC SITE', ## 1qnm
    '1':'ACTIVE SITE',
    '2':'ACTIVE SITE',
    '3':'ACTIVE SITE',
    'CR1':'CATALYTIC RESIDUE', ## 1dub
    'CR2':'CATALYTIC RESIDUE', ## 1dub
    }

pdbs = []
fd = open('REMARK800.txt','r')
lines = fd.readlines()
fd.close()
for line in lines:
    pdbs += [line[:4]]
fd = open('SITE.txt','r')
lines = fd.readlines()
fd.close()
for line in lines:
    pdbs += [line[:4]]

subdirs = os.listdir(pdbpath)
if sublim == '1':
    subdirs.sort()
elif sublim == '2':
    subdirs.sort()
    subdirs.reverse()
elif sublim == '3':
    None
else:
    subdirs.sort()
for subdir in subdirs:
    if subdir < sublim and sublim not in ['1','2','3']:
        continue
    files = os.listdir(pdbpath+subdir)
    for file in files:
        if file[:-2] == 'gz':
            continue
        pdb = file[3:7]
        SITE_siteIDs = set()
        set_siteIDs = set()
        if pdb in pdbs:
            continue
        if pdb in [
            ## errors
            '1upw', '1vz8', '1ng1', '1ha2', '1a2a', '1wcb', '1ahp', '1qbn', '2dpg', '1at6', '2jeb', '1nbm', '1at5', '1e7b', '2c7f', '2jgf', '1cuw', '2cix', '1wcv', '2gar', '2cdn', '1a9l', '1h9z', '3gar', '2uxn', '2bir', '1lb3', '2cmy', '2j4h', '1up0', '1okw', '1o7a', '1gq2', '1w1a', '1w1y', '1h3v', '1e56', '1w7m', '1qb1', '1ea7', '1oih', '2a9l', '1ay0', '1ol2', '1l9w', '1dzy', '1idk', '1hcx', '2ca5', '1slo', '1slp', '1tez', '1ier', '2ng1', '1gzm', '1akd', '1jnq', '1e90', '2ldz', '1dzx', '1uvj', '1uu4', '1w73', '2j3z', '1gwm', '1h8s', '2pah', '2jcm', '1dms', '1uxj', '1iib', '1h10', '1ldz', '1oj1', '1e9o', '1h6h', '1qc5', '1et4',
            ## COB
            '1cqe', ## HYDROPHOBIC CHANNEL RESPONSIBLE FOR CYCLOOGENASE ACTIVITY
            ## S4
            '1pau','1ibc', ## INHIBITOR BINDING SUB-SITE S4
            ## B1
            '3rsk','3rsp','3rsd','4rsd', ## PYRIMIDINE BINDING SITE ON 5' SIDE OF SCISSILE PHOSPHATE
            '2cam','1rav', ## RESIDUES INVOLVED IN BIOTIN HYDROGEN BONDING
            ## AB1
            '1rkd', ## RESIDUES MAKING DIRECT HYDROGEN BONDS WITH ADP
            ## AC1
            '1wau', ## SO4 BINDING SITE
            '1h8d', ## TRIPEPTIDE PHOSPHONATE INHIBITOR BINDING
            '1ndr','1nds', ## MONOMER A TYPE 1 COPPER SITE
            ## AC2
            '1gju', ## PO4 BINDING SITE FOR CHAIN A
            ## CAT
            '3pgh','1cx2', ## TYR 385 IS BELIEVED TO BE THE AMINO ACID THAT ABSTRACTS A HYDROGEN ATOM FROM THE SUBSTRATE. IT IS LOCATED CLOSE TO THE HEME. A TYROSINE RADICAL IS FORMED DURING THE COURSE OF THE REACTION
            '1oeo', ## OXIDIZED TO SULFONIC ACID
            '1oem', ## SULFENYL-AMIDE BOND
            '1qjj', ## HIS92, HIS96, HIS 102 ARE ZINC LIGANDS TYR 149, WHICH IS AN ADDITIONAL ZINC LIGAND IN THE UNINHIBITED ENZYME IS ONLY SLIGHTLY SHIFTED AS COMPARED TO THE PHOSPHINIC PEPTIDE BOUND ENZYME SEE SEPARATE PDB ENTRY, 1QJI)
            '1b5e', ## SG ATOM OF CYS 148 HAS DUAL CONFORMATION
            '4cox','5cox','6cox', ## TYR 385 IS BELIEVED TO BE THE AMINO ACID THAT ABSTRACTS A HYDROGEN ATOM FROM THE SUBSTRATE. IT IS LOCATED CLOSE TO THE HEME. A TYROSINE RADICAL IS FORMED DURING THE COURSE OF THE REACTION
            ## ACT
            '1zap', ## ASPARTIC PROTEINASES ARE CHARACTERIZED BY TWO ASP RESIDUES, ONE FROM EACH DOMAIN, CORRESPONDING TO PEPSIN ASP 32, ASP 218
            '1ame', ## ICE BINDING RESIDUES
            '1cci','1ccj', ## REMOVAL OF PHE 202 FORMS AN INTERNAL CAVITY ADJACENT TO ASP 235
            '1ceg','1cef', ## THOSE FOR S. R61 ARE LISTED BELOW
            '1auk', ## RESIDUE 69 WAS TREATED AS A GLYCINE DURING REFINEMENT TO AVOID BIAS IN THE INTERPRETATION OF THE DIFFERENCE ELECTRON DENSITY MAP. THE SIDE CHAIN OF RESIDUE 69 WAS INTERPRETED AS AN ALDEHYDE GROUP WITH TWO-FOLD DISORDERED ALDEHYDE FUNCTION (I.E. THE CARBONYL OXYGEN ATOM OCCUPIES TWO DIFFERENT POSITIONS). THE AUTHORS ASSUME THAT THE ALDEHYDE IS IN EQUILIBRIUM WITH ITS HYDRATED FORM, THE GEMINAL DIOL, WHICH IS IN ACCORDANCE WITH THE SHAPE OF THE ELECTRON DENSITY
            '1avu', ## P1 SITE
            '1gmy', ## FORMS THIOIMIDATE BOND WITH INHIBITOR
            ## CA2
            '1hhs','1hht','1hi8','1hi1','1hi0', ## THE THREE CONSERVED ACTIVE SITE ASPARTATES
            ## A
            '1bsx', ## NR-BOX BINDING SITE
            ## D
            '1bhn', ## C-GMP AND GDP ARE FOUND (EACH WITH HALF OCCUPANCY) IN THE NUCLEOTIDE BINDING SITE
            ## 3
            '2rbi', ## MAIN NUCLEOTIDE BINDING RESIDUE
            ## ASI
            '1rcy', ## DISTORTED TETRAHEDRAL COORDINATION BY 2 HIS, 1 CYS, AND 1 MET RESIDUES OF THE CU(II) ATOM (CU)
            ## NUC
            '1bvv', ## THE NUCLEOPHILE TO WHICH THE INHIBITOR IS BOUND IS GLU 78
            '1hlu', ##  ATP BINDING SITE, TAKEN FROM PDB ENTRY 2BTF
            ## ZN
            '1f3z', ## INCLUDES GLU 148 FROM A NEIGHBORING MOLECULE AND ONE SOLVENT MOLECULE
            ## COA
            '4gsa', ## REDUCED SCHIFF BASE LINKAGE TO LYS 273 IN CHAIN A
            ## PLA
            '1cs1', ## COFACTOR SITE CHAIN A
            ## NA1
            '1b14', ## NAD BINDING MOTIF IN DADHS G(A)XGXX
            ## NB1
            '1b15', ## NAD BINDING MOTIF IN DADHS G(A)XGXX
            ## FA1
            '1qlt','1qlu', ## COFACTOR
            ## MO
            '1dmr','3dmr','2dmr', ## SER SITE
            ## CIC
            '2plc', ## INOSITOL BINDING SITE
            '1a4l','1a4m', ## THE SITE BINDS ADENOSINE AND CONVERT IT TO INOSINE AND AMMONIAGENE.
            ## CDE
            '1aqh','1aqm', ## CHLORIDE IS AN ACTIVATOR FOR THIS ENZYME
            ## SPL
            '1slp', ## THIS MOLECULE IN VIVO IS CLEAVED BETWEEN G22 AND G23
            ## ASA
            '1bkj', ## FMN COFACTOR BINDING SITE-CHAIN A
            ## ASB
            '1a6e', ## RESIDUE NUMBERING ACCORDING TO ALPHA-TYPE SUBUNIT
            '1cjx', ## IRON-BINDING SITE
            '2daa', ## ESSENTIALLY IDENTICAL TO ASA
            '4daa', ## ESSENTIALLY THE SAME AS ASA
            ## ACC
            '1qnt', ## ALKYL ACCEPTOR
            ## AS1
            '1dw9', ## CHLORIDE IONS K19 K23 RESIDUES SER122A, ARG96D, ARG96I, SER122J
            ## ZNA
            '1an8', ## THIS SITE ALSO INCLUDES HIS 81 FROM A SYMMETRY-RELATED MOLECULE
            ## HMA
            '1ar1', ## AXIAL HEME A LIGANDS
            ## HTH
            '1uxd','1uxc', ## HTH_LACI FAMILY
            ## SBS
            '1shk', ## THE ELECTRON DENSITY FOR SHIKIMATE WAS AMBIGUOUS PREVENTING ITS INCLUSION IN THE MODEL. THE LISTED RESIDUES ARE GROUPED AROUND THE DENSITY AND ARE MOST LIKELY TO BE PART OF THE SHIKIMATE BINDING
            
            ## one word difference
            ## A MOIETY MODELLED AS A 12-CARBON FATTY ACID IS OBSERVED IN FULL OCCUPANCY IN THE HYDROPHOBIC POCKET WHERE ANTI-VIRAL COMPOUNDS BIND. THE HEAD GROUP OXYGENS CONTACT N1212 AND L1100
            ## A MOIETY MODELLED AS A 12-CARBON FATTY ACID IS OBSERVED IN PARTIAL OCCUPANCY IN THE HYDROPHOBIC POCKET WHERE ANTI-VIRAL COMPOUNDS BIND. THE HEAD GROUP OXYGENS CONTACT N1212 AND L1100
            '1ayn',
            ## remediation error introduction
            '2ahj','1bjo','1bjn','2bl2','1boy','1brm','1e1h','1e71',
            '1e4h','1gq1','1h7t','1h8e','1ha5',
            '1nlr',
            ## error
            '2jer','1rbl','1rsc','3tat','12as','1ak9',
            ## remediation merger
            '1azy','1byo','1qo5','1hxp',
            ## ':' in them
            '1e8m','1e8n','1e3d','2fus','1fur','1elv','1h2w','1i4q','1i4r',
            '1i4p','1kmn','1kmm','1qjk','1qjl','1qji','1a0e','1a3c','1a4x',
            '1ba9','1ao0','1aye','2cxg','1htt',
            ]:
            continue
        fd = open(pdbpath+subdir+'/'+file,'r')
        lines = fd.readlines()
        fd.close()
##        j = 0
##        for i in range(j, len(lines)):
        for i in range(len(lines)):
            line = lines[i]
            if line[6:].strip() == 'ATOM':
                break
            if line[:10] == 'REMARK 800':
                print line[:-1]
                if line[11:27] == 'SITE_IDENTIFIER:':
                    siteIDs = line[28:80].strip().split(',')
                    if siteIDs[-1][-1] == '.':
                        siteIDs[-1] = siteIDs[-1][:-1]
                    if siteIDs[-1] != 'AND' and siteIDs[-1].split()[0] == 'AND' and len(siteIDs[-1].split()) == 2:
                        siteIDs[-1] = siteIDs[-1].split()[1]
                    if len(siteIDs) == 1 and len(siteIDs[0].split()) == 3 and siteIDs[0].split()[-2] == 'AND':
                        siteIDs = [siteIDs[0].split()[0],siteIDs[0].split()[2]]
                    for j in range(len(siteIDs)):
                        siteID = siteIDs[j]
                        siteID = siteID.strip()
                        siteIDs[j] = siteID
                        if len(siteID) > 3 or len(siteID) < 1:
                            print pdb, siteID, siteIDs
                            stop
                    set_siteIDs |= set(siteIDs)
                elif line[11:28] == 'SITE_DESCRIPTION:':
                    print pdb
                    site_description = line[29:80].strip()
                    for j in range(i+1,len(lines)):
                        line = lines[j]
                        if line[:10] != 'REMARK 800' or line.strip() in ['REMARK 800','REMARK 800 SITE'] or line[11:28] == 'SITE_DESCRIPTION:' or line[11:27] == 'SITE_IDENTIFIER:':
                            if line[11:28] == 'SITE_DESCRIPTION:' or line.strip == 'REMARK 800 SITE':
                                print pdb, line
                                notexpected
                            line = lines[i]
                            break
                        site_description += ' '+line[11:80].strip()
                        if line[10] != ' ':
                            print line
                            notexpected
                    if ':' in site_description:
                        print site_description
                        print ': in site_description'
                        stop1
                    if site_description in ['','NULL']:
                        continue

################################################################################
                    print '1', site_description
                    ##
                    ## prefixes
                    ##
                    if site_description[:4] == 'THE ':
                        site_description = site_description[4:]
                    if site_description[:3] == 'AN ':
                        site_description = site_description[3:]
                    if site_description[:2] == 'A ':
                        site_description = site_description[2:]
                    ##
                    ## suffixes
                    ##
                    if site_description[-1] == '.':
                        site_description = site_description[:-1]
                    if site_description[-1] == ' ':
                        site_description = site_description[:-1]
                    if site_description[-15:-1] == ' FOR MOLECULE ':
                        site_description = site_description[:-15]
                    if site_description[-12:-1] == ' FOR CHAIN ':
                        site_description = site_description[:-12]
                    if site_description[-14:-1] == ' FOR SUBUNIT ':
                        site_description = site_description[:-14]
                    if site_description[-14:-1] == ' IN MOLECULE ':
                        site_description = site_description[:-14]
                    if site_description[-14:-1] == ' OF MOLECULE ':
                        site_description = site_description[:-14]
                    if site_description[-11:-1] == ' IN CHAIN ':
                        site_description = site_description[:-11]
                    if site_description[-13:-1] == ' IN MONOMER ':
                        site_description = site_description[:-13]
                    if site_description[-11:-1] == ' OF CHAIN ':
                        site_description = site_description[:-11]
                    if site_description[-11:-7] == ' IN ' and site_description[-6:] == ' CHAIN':
                        site_description = site_description[:-11]
                    if site_description[-14:-10] == ' IN ' and site_description[-9:] == ' MOLECULE':
                        site_description = site_description[:-14]
                    if site_description[-9:-7] == ', ' and site_description[-6:] == ' CHAIN':
                        site_description = site_description[:-9]
                    if site_description[-8:-1] == ' CHAIN ':
                        site_description = site_description[:-8]
                    if site_description[-10:-2] == ' (CHAIN ' and site_description[-1] == ')':
                        site_description = site_description[:-10]
                    if site_description[-7:-1] == ' SITE ':
                        site_description = site_description[:-2]
                    if site_description[-13:-1] == ' OF SUBUNIT ':
                        site_description = site_description[:-13]
                    if site_description[-1] == ',':
                        site_description = site_description[:-1]
                    if site_description[-1] == '.':
                        site_description = site_description[:-1]
##                    if 'SUBSTRATE (' in site_description and ')' in site_description:
##                        index1 = site_description.index('(')
##                        index2 = site_description.index(')')
##                        if site_description[index1-1] == ' ':
##                            print site_description
##                            site_description = site_description[0:index1-1]+site_description[index2+1:]
##                            print site_description
##                        else:
##                            stop2
                    print '2', site_description
################################################################################
                    if site_description in ['DESCRIPTION NOT PROVIDED','NO DESCRIPTION PROVIDED']:
                        continue
                    if site_description == 'ELECTRON DENSITY FOR SHIKIMATE WAS AMBIGUOUS PREVENTING ITS INCLUSION IN THE MODEL. THE LISTED RESIDUES ARE GROUPED AROUND THE DENSITY AND ARE MOST LIKELY TO BE PART OF THE SHIKIMATE BINDING' and siteID == 'SBS':
                        stop3
                    for siteID in siteIDs:
                        if siteID == 'AND':
                            stop
                        if not siteID in d_site_descriptions.keys():
                            if site_description[:11] == 'ACTIVE SITE':
                                site_description = 'ACTIVE SITE'
                            d_site_descriptions[siteID] = site_description
                        elif d_site_descriptions[siteID] != site_description:
                            if len(set(l_site_descriptions) & set([site_description,d_site_descriptions[siteID]])) != 2:
                                if '.' in site_description and '.' in d_site_descriptions[siteID] and site_description[:site_description.index('.')] == d_site_descriptions[siteID][:d_site_descriptions[siteID].index('.')]:
                                    None
                                elif (
                                    (
                                        ('SCHIFF-BASE' in site_description and 'LYSINE' in site_description and 'PLP' in site_description) or
                                        'SCHIFF BASE WITH SUBSTRATE' in site_description or
                                        'SCHIFF-BASE WITH SUBSTRATE' in site_description or
                                        'CLEAVAGE SITE' in site_description or
                                        ' CAT SITE' in site_description or
                                        'POCKET' in site_description or
                                        'GENERAL ACID BASE' in site_description or
                                        'ACID/BASE CATALYS' in site_description or
                                        'REACTION INTERMEDIATE' in site_description or
                                        'TRANSITION STATE' in site_description or
                                        'CATALY' in site_description or
                                        'METALLOCENTER' in site_description or
                                        'BINDING SITE' in site_description or
                                        'ACTIVE SITE' in site_description or 
                                        'ACTIVE-SITE' in site_description or 
                                        'SUBSTRATE BINDING RESIDUES' in site_description or
                                        site_description[:6] == 'ACTIVE' or
                                        site_description[-12:] == 'BINDING SITE' 
                                        ) and
                                    (
                                        d_site_descriptions[siteID] in l_site_descriptions
                                        )
                                    ):
                                    None
                                elif site_description.replace(' SITE','') in d_site_descriptions[siteID].replace(' SITE',''):
                                    d_site_descriptions[siteID] = site_description
                                elif d_site_descriptions[siteID].replace(' SITE','') in site_description.replace(' SITE',''):
                                    None
                                elif 'BINDING SITE' in site_description and 'BINDING SITE' in d_site_descriptions[siteID]:
                                    None
                                elif (
                                    siteID in ['CA','CAA','CAB','CAC','CA1','CA2','CA3','CA4','CAL',] and
                                    (
                                        'CA2+' in site_description or
                                        'CALCIUM' in site_description or
                                        ' CA BINDING SITE' in site_description or
                                        'CA BINDING SITE' == site_description[:len(' CA BINDING SITE')] or
                                        'CA ' == site_description[:3] or 
                                        ' CA ' in site_description or
                                        '(CA ' in site_description or
                                        site_description == 'SITE'
                                        )
                                    ):
                                    None
                                elif siteID in ['ZNA','ZNB','ZNC','ZN1','ZN2','ZN3','ZN4','ZN','ZNS','CTA','CTB','CTC','CTD','ZIN',] and (
                                    'ZINC' in site_description or
                                    'ZN' in site_description or
                                    'METAL' in site_description
                                    ):
                                    None
                                elif siteID in ['MNA','MNB','MN1','MN2','MN3','MN4',] and ('MANGANESE' in site_description or ('MN' in site_description and 'ION' in site_description)):
                                    None
                                elif siteID in ['MG','MGA','MG1',] and (' METAL ' in site_description or 'MAGNESIUM' in site_description or 'MG' == site_description[:2]):
                                    None
                                elif siteID in ['COB','CO2',] and ('COBALT' in site_description):
                                    None
                                elif siteID in ['COP','CUA','CUB','CU','CU1','CU2','CU3','CU4','CU5','CU6',] and ('COPPER' in site_description or 'CU' in site_description):
                                    None
                                elif siteID in ['CU2',] and ('TYPE II' in site_description):
                                    None
                                elif siteID in ['CD','CD1','CD2','CD4','CD5',] and ('CADMIUM' in site_description or 'CD ' == site_description[:3]):
                                    None
                                elif siteID in ['FE','FE1','FE2','FE3','FE4','FEA','FEB','FEC','FES',] and (
                                    ' FE BINDING SITE' in site_description or
                                    'FE' == site_description[:2] or
                                    'IRON' in site_description or 
                                    'METAL BINDING SITE' in site_description
                                    ):
                                    None
                                elif siteID in ['MO4','MO',] and 'MOLYBDENUM' in site_description:
                                    None
                                elif siteID in ['FMN',] and ('FLAVIN' in site_description or 'FLAVIN' or 'PROSTHETIC' in site_description):
                                    None
                                elif siteID in ['CSH',] and ('CHROMOPHORE' in site_description):
                                    None
                                elif siteID in ['FS4','22','FES','FSB',] and ('SULFUR CLUSTER' in site_description or 'IRON SULFUR' in site_description or 'IRON-SULFUR' in site_description or 'IRON SULPHUR' in site_description or 'CLUSTER' in site_description):
                                    None
                                elif siteID in ['PLP',] and ('LYSINE' in site_description or 'PYRIDOXAL' in site_description or 'SCHIFF' in site_description):
                                    None
                                elif siteID in ['K2',] and ('POTASSIUM' in site_description):
                                    None
                                elif siteID in ['NA','NA1','NAA','NAB',] and ('SODIUM' in site_description or 'NA BINDING SITE' in site_description):
                                    None
                                elif siteID in ['PO1','PO4','PBS',] and ('PO4' in site_description or 'PHOSPHATE' in site_description or 'PHOSPHORYLATION' in site_description):
                                    None
                                elif siteID in ['C93',] and ('UBIQUITIN' in site_description):
                                    None
                                elif siteID in ['FAD','FAA','FAB','FAC','FAE','FAF',] and ('COFACTOR' in site_description or 'FAD BINDING SITE' in site_description):
                                    None
                                elif siteID in ['NAD',] and ('COFACTOR' in site_description or 'NAD BINDING SITE' in site_description):
                                    None
                                elif siteID in ['FMN',] and ('FMN CO-FACTOR' in site_description):
                                    None
                                elif siteID == 'CT1' and ('CATION-BINDING SITE' in site_description):
                                    None
                                elif siteID == 'ACE' and ('ACETATE' in site_description):
                                    None
                                elif siteID in ['SB1','SB2','SB3',] and ('ACARBOSE' in site_description):
                                    None
                                elif siteID in ['OHB',] and (site_description[:3] == 'CD '): ## 1hf3
                                    None
                                elif siteID in ['COA',] and (site_description == 'COA BINDING SITE'):
                                    None
                                elif siteID in ['MBA','MBB',] and site_description[:5] == 'METAL':
                                    None
                                elif siteID == 'POC' and ('POCKET' in site_description):
                                    None
                                elif siteID == 'REG' and ('REGULATORY' in site_description):
                                    None
                                elif siteID == 'GL1' and ('GLYCOSYLATION SITE' == site_description):
                                    None
                                elif siteID == 'REA' and ('REACTIVE' in site_description):
                                    None
                                elif siteID == 'TRI' and ('CATALYTIC TRIAD' in site_description):
                                    None
                                elif siteID == 'NAD' and ('DINUCLEOTIDE' in site_description):
                                    None
                                elif siteID == 'OXY' and ('OXYANION' in site_description):
                                    None
                                elif siteID in ['HE1','HE2','HE3','HE4','HMA','HEM',] and ('HEME' in site_description or 'HAEM' in site_description):
                                    None
                                elif siteID == 'NUC' and ('CATALYTIC NUCLEOPHILE' in site_description):
                                    None
                                elif siteID == 'RNA' and ('RNA' == site_description[:3]):
                                    None
                                elif siteID in ['NUL','AVE','1','S1','2','S2']:
                                    None
                                else:
                                    fd = open('siteIDconflicts.txt','a')
                                    fd.write('pdb %s\nsID %s\nstr %s\ndic %s\n\n' %(pdb,siteID,site_description,d_site_descriptions[siteID]))
                                    fd.close()
                elif line.strip() in ['REMARK 800','REMARK 800 SITE']:
                    continue
            if line[:6].strip() == 'SITE':
                site_name = line[11:14].strip()
                SITE_siteIDs |= set([site_name])
                set_site_names |= set([site_name])
                if not site_name in d_site_names:
                    d_site_names[site_name] = set()
                d_site_names[site_name] |= set([pdb])

################################################################################
for site_name in d_site_names.keys():
    set_pdbs = d_site_names[site_name]
    if len(set_pdbs) == 1:
        del d_site_names[site_name]
l_sites = d_site_names.keys()
print d_site_names

##d_sites = {'ACC': 'ACTIVE SITE', 'ACB': 'ACTIVE SITE', 'ACA': 'ACTIVE SITE', 'ACG': 'ACTIVE SITE', 'ACE': 'ACTIVE SITE', 'ACH': 'ACTIVE SITE', 'ACL': 'ACTIVE SITE', 'ACT': 'ACTIVE SITE OF THE PPIASE DOMAIN IS MARKED BY THE BOUND ALA-PRO DIPEPTIDE (CHAIN B)', 'CAT': 'ACTIVE SITE', 'CAB': 'ACTIVE SITE', 'CAA': 'ACTIVE SITE'}
d_sites = d_site_descriptions
for site in d_sites.keys():
    if d_sites[site] == 'BINDING SITE':
        del d_sites[site]
        continue
    if d_sites[site].split()[0] not in ['ACTIVE','CATALYTIC']:
        print d_sites[site]
        del d_sites[site]
        continue
    if site not in l_sites:
        print site
        del d_sites[site]
        continue
print d_sites
stop

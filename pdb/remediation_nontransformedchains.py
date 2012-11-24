'''this script checks if a polymer is present (ATOM/SEQRES) which is not part of any of the transformed biomolecules (checks if present in the REMARK350 records)'''

import os, Numeric, sys, math
pdbpath = '/media/WDMyBook1TB/2TB/pdb/'

def main():

    subdirs = os.listdir(pdbpath)
    subdirs.sort()
    for subdir in subdirs:
        if subdir < sys.argv[-1]:
            continue
        print subdir
        files = os.listdir(pdbpath+subdir)
        files.sort()
        for file in files:
            if file[-2:] == 'gz':
                continue
            pdb = file[3:7]
            if pdb in [
##                ## nontransformations
##                '429d','43ca','13pk','454d','258d','474d','375d','378d','198d',
##                '1a0o','2a01','2a1m','2a1n','2a1o','2a27','2a2a','2a2o','2a2z',
##                '2a30',
##                ## other
##                '363d',
##                ## dna/rna
##                '2a04',
##                ## multiple chains, one ligand
##                '2a0q',
##                ## no ligands for one chain
##                '2a6c',
##                ## fixed
##                '2a3a','2a3b','2a3c',
##                ## not fixed
##                '2a38',
##                ## ligand contacts more than one chain
##                '1cqp','2b50',
##                ## standard error not fixed
##                '2a3e','2a3r','2a3w','1a4k',
##                ## standard error not fixed (ligands not a multiplum of chains)
##                '1com',
##                ## standard but ligands linked
##                '1b37',
                ]:
                continue
            fd = open(pdbpath+subdir+'/'+file,'r')
            lines = fd.readlines()
            fd.close()

            remark350chains = set()
            atomchains = set()
            d_biomolecules = {}
            d_coordsnotbm = {} ## by resno (and iCode and chain)
            d_coordsbm = {} ## by chain (or by biomolecule)
            d_link = {}
            longchain = False

            for i in range(len(lines)):
                line = lines[i]
                if line[:10] == 'REMARK 350':
                    if line[11:23] == 'BIOMOLECULE:':
                        biomolecules = set(line[23:80].replace(' ','').split(','))
                        for j in range(i+1,len(lines)):
                            if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
                                break
                            elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                                chains = set()
                                line_chains = lines[j][41:80]
                                chains |= parse_REMARK350_chains(line_chains)
                                remark350chains |= chains
                            elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
                                line_chains = lines[j][53:80]
                                chains |= parse_REMARK350_chains(line_chains)
                                remark350chains |= chains
                            elif ',' in lines[j][11:80]:
                                if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
                                    line_chains = lines[j][11:80]
                                    chains |= parse_REMARK350_chains(line_chains)
                                    remark350chains |= chains
                        for biomolecule in biomolecules:
                            d_biomolecules[biomolecule] = chains

                elif line[:6].strip() == 'SEQRES':
                    if line[19:80].strip().split() > 12:
                        longchain = True

                elif line[:6].strip() == 'LINK':

                    hetID1 = line[17:20].strip()
                    chain1 = line[21]
                    resno1 = int(line[22:26])
                    hetID2 = line[47:50].strip()
                    chain2 = line[51]
                    resno2 = int(line[52:56])
                    d_link[resno1] = resno2
                    d_link[resno2] = resno1

                elif line[:6].strip() in ['ATOM','HETATM']:
                    resname = line[17:20].strip()
                    chain = line[21]
                    resno = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord = Numeric.array([x, y, z])
                    if len(d_biomolecules.keys()) > 0 and chain == ' ' and resname not in ['HOH','ZN','NA','CL','SO4','MG']:
                        print line,
                    if resname not in [
                        ## solvents and ions (transfer to quakes.py)
                        'HOH','BME','DTT','UNX','MOH','IPA','EDO',
                        'FMT','GOL','EEE',
                        ]+{
##                            ## acetate
                            'ACT':['C2 H3 O2',-1],
##                            ## nitrate, ammonium
##                            'NO3':['N1 O3',-1],'NH4':['H4 N1',+1],
##                            ## hydroxide
##                            'OH' :['H1 O1',-1],
##                            ## phosphate
##                            '2HP':['O4 P1',-1],'PI' :['H1 O4 P1',-2],
                            'PO4':['O4 P1',-3],
##                            ## sulfate
                            'SO4':['O4 S1',-2],
##                            'SOH':[3,-1],'SUL':[3,-2], ## different oxidation states
##                            ## carbonate
##                            'CO3':['C O3', -1],
##                            ## group1a
##                            'LI' :['LI1',+1],
                            'NA' :['NA1',+1], 
                            'K'  :['K1' ,+1], 
                            'CS' :['CS' ,+1],
##                            ## group2a
##                            'BEF' :['BE F3',-1],
                            'MG' :['MG1',+2],
                            'CA' :['CA1',+2],'OC1':['CA1',+2], ## 'OC3','OC5','OC6','OC7','OC8','OC2','OC4' deprecated
                            'SR':[],
##                            'BA' :['BA' ,+2],
##                            ## group3a
##                            'AL' :['AL1',+3],'ALF' :['AL F4',-1],
##                            'GA' :['GA1',+3],'GA1' :['GA1',+2],
##                            'TL' :['TL1',+1],
##                            ## group4a
##                            'ARS':['AS1', 0],'ART':['O4 AS1',-3],'AST':-3,'TAS':['H3 O3 AS1', 0],'CAC':['C2 H6 AS O2', -1], ## different compounds
##                            'PB' :['PB' ,+2],
##                            ## group6a
##                            'SE' :['SE1', 0],'SE4':['O4 SE1',-2], ## different compounds
##                            ## group7a
                            'CL' :['CL1',-1],
                            'BR' :['BR1',-1],
##                            'IOD':['I1' ,-1],
##                            ## group8a
##                            'KR' :['KR1', 0],
##                            ## group3b
##                            'V'  :+3,'VO4':['V1' ,-3], ## different oxidation states
##                            ## group4b
##                            'CR' :['CR1',+3],
##                            'MO' :['MO1', 0],'4MO':['MO1', 0],'2MO':['MO O2',-2], ## different compounds and different oxidation states
##                            ## group5b
##                            'MN' :['MN1',+2],'MN3':['MN1',+3], ## different oxidation states; 'MW1','MW2','MW3','O4M','MN5','MN6' deprecated
##                            ## group6b
##                            'FE2':['FE1',+2],'FE' :['FE1',+3], ## different oxidation states; 'OF1','OF3','2OF' deprecated
##                            ## group7b
##                            'CO' :['CO1',+2],'3CO':['CO1',+3], ## different oxidation states; 'CO5','OCL','OCO','OCN','OCM' deprecated
##                            ## group8b
                            'NI' :['NI1',+2],'3NI':['NI1',+3], ## different oxidation states; 'NI1','NI2','NI3','NIK' deprecated
##                            ## group9b
##                            'CU1':['CU1',+1],'CU' :['CU1',+2], ## different oxidation states; '1CU' deprecated
##                            ## group10b
##                            'ZN' :['ZN1',+2], ## 'ZN2','ZO3','ZN3','ZNO' deprecated
##                            'CD' :['CD1',+2],
##                            'HG' :['HG1',+2],
##                            ## Lanthanides
##                            'TB' :['TB1',+3],
##                            'YB' :['YB1',+3],
                            }.keys():
                        atomchains |= set([chain])
                        if chain in remark350chains:
                            if not chain in d_coordsbm.keys():
                                d_coordsbm[chain] = []
                            d_coordsbm[chain] += [coord]
                        else:
                            if not resno in d_coordsnotbm.keys():
                                d_coordsnotbm[resno] = []
                            d_coordsnotbm[resno] += [coord]

            if longchain == False:
                continue

            if atomchains != remark350chains and len(d_biomolecules.keys()) > 1 and not len(d_biomolecules.keys()) > len(d_coordsnotbm.keys()):
                if 'RIBONUCLEIC' not in lines[0]:
                    print pdb
                    print 'ATOM', atomchains, 'REMARK350', remark350chains, d_biomolecules.keys()
                    print 'ATOM', atomchains-remark350chains, 'REMARK350', remark350chains-atomchains
                    d_dist = {}
                    for resno in d_coordsnotbm.keys():
                        d_dist[resno] = {}
                        for biomolecule in d_biomolecules.keys():
                            d_dist[resno][biomolecule] = 'N/A'
                            for chain in d_biomolecules[biomolecule]:
                                for c1 in d_coordsbm[chain]:
                                    for c2 in d_coordsnotbm[resno]:
                                        dist = math.sqrt(sum((c1-c2)**2))
                                        if dist < d_dist[resno][biomolecule]:
                                            d_dist[resno][biomolecule] = dist
                    for resno in d_dist.keys():
                        print resno, d_dist[resno]
                        count = 0
                        for bm in d_dist[resno].keys():
                            if d_dist[resno][bm] < 3.5:
                                count += 1
##                                print resno,bm,d_dist[resno][bm]
                            if count > 1:
                                print pdb
                                print lines[0]
                                oneligandmultiplechains
                        if count == 0 and resno not in d_link.keys():
                            print pdb
                            print 'resno', resno
                            print d_link
                            print 'distances between chains...', d_dist
                            print pdb
                            stop2
                                
                    print pdb
                    print lines[0]
                    if len(d_coordsnotbm.keys()) % len(d_biomolecules.keys()) != 0:
                        print len(d_coordsnotbm.keys()) % len(d_biomolecules.keys())
                        print len(d_coordsnotbm.keys()), len(d_biomolecules.keys())
                        stop3

    return

def parse_REMARK350_chains(line_chains):

    ## if sentence necessary due to e.g. 1qgc
    ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: 1 2 3 4 5
    if ',' not in line_chains:
        chains = line_chains.split()
    else:
        ## remove 'AND' from the line of chains (e.g. problem with 1rhi)
        ## replace '.' in the line of chains (e.g. problem with 1rbo and 1qgc)
        chains = line_chains.replace('AND',',').replace('.',',').replace(' ',',').split(',')

    ## loop removal of blank chains necessary due to e.g. 2g8g
    ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, ,
    for x in range(100):
        if '' in chains:
            chains.remove('')
        else:
            break

    for j in range(len(chains)):
        chain = chains[j]
        if chain == 'NULL':
            chains[j] = ' '

    return set(chains)

if __name__ == '__main__':
    main()

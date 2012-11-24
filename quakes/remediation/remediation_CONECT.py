#!/software/bin/python
#$Id: remediation_CONECT.py,v 1.1 2007/09/07 16:05:54 tc Exp $
#Tommy Carstensen, University College Dublin, 2007

def main():

    import os, math

    pdbpath = '/oxygenase_local/data/pdb/'

    pdbs = []
    subdirs = os.listdir(pdbpath)
    subdirs.sort()
    for subdir in subdirs:
        files = os.listdir(pdbpath+subdir)
        for file in files:
            pdbs += ['%s' %(file[3:7])]
    pdbs.sort()

    d_hydrogen_equivalents = {
        'C' :2,'O' :0,'H': 0,
        'F' :0,'CL':0,'BR':0,'I' :0,
        'N' :1,'P' :1,
        }
    d_valence = {
        'C' :4,'O' :2,'N' :3,
        }

    lines_out = []

##    pdbs = ['148l', '148l', '1a14', '1b9s', '1b9s', '1b9t', '1b9t', '1b9v', '1b9v', '1bcs', '1bcs', '1bcs', '1bcs', '1bcs', '1byb', '1ce7', '1ce7', '1cel', '1cel', '1cpy', '1cpy', '1cpy', '1cpy', '1cpy', '1cpy', '1dl2', '1dl2', '1dl2', '1dl2', '1dl2', '1dl2', '1dpj', '1dpj', '1e04', '1e04', '1e4k', '1e4k', '1e4m', '1e4m', '1e4m', '1e4m', '1e4m', '1e4m', '1e6q', '1e6q', '1e6q', '1e6q', '1e6q', '1e6q', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e6x', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e70', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e71', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e72', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1e73', '1ea5', '1ea5', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqg', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1eqh', '1f8b', '1f8b', '1f8c', '1f8c', '1f8d', '1f8e', '1fe8', '1fe8', '1fe8', '1ffr', '1fne', '1fne', '1fne', '1fne', '1fne', '1fne', '1fne', '1fne', '1fne', '1fne', '1fng', '1fng', '1fng', '1fng', '1fng', '1fng', '1fng', '1fng', '1fng', '1fng', '1fzc', '1fzc', '1fzd', '1fzd', '1fzd', '1fzd', '1fzd', '1fzd', '1fzf', '1fzf', '1g0v', '1g0v', '1gah', '1gah', '1gah', '1gah', '1gai', '1gai', '1gai', '1gai', '1gqh', '1gqh', '1gqh', '1gqh', '1h3v', '1h3v', '1h3v', '1h3v', '1h81', '1h81', '1h81', '1h81', '1h81', '1h81', '1h81', '1h81', '1h82', '1h82', '1h82', '1h82', '1h82', '1h82', '1h82', '1h82', '1h83', '1h83', '1h83', '1h83', '1h83', '1h83', '1h83', '1h83', '1h84', '1h84', '1h84', '1h84', '1h84', '1h84', '1h84', '1h84', '1h86', '1h86', '1h86', '1h86', '1h86', '1h86', '1h86', '1h86', '1hcn', '1hcn', '1hcn', '1hcn', '1hcn', '1hcn', '1hcn', '1hcn', '1hfu', '1hfu', '1hfu', '1hfu', '1hfu', '1hfu', '1hfu', '1hfu', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht5', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1ht8', '1hzh', '1hzh', '1hzh', '1hzh', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1i8l', '1ing', '1k5c', '1k5c', '1kbk', '1kbk', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kwf', '1kx0', '1kx0', '1kx0', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1kya', '1lqs', '1lqs', '1lqs', '1lqs', '1lwu', '1lwu', '1lwu', '1lwu', '1lwu', '1lwu', '1lwu', '1lwu', '1lwu', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1m2t', '1mh0', '1mh0', '1mh0', '1mh0', '1mpm', '1mpm', '1mpm', '1mpm', '1mpm', '1mpm', '1mpm', '1mpm', '1mpm', '1mwa', '1mwa', '1mwe', '1mwe', '1mx1', '1mx5', '1nfd', '1nfd', '1nm9', '1nn2', '1nn2', '1nn2', '1nnc', '1nql', '1nql', '1nql', '1nu8', '1nu8', '1nu8', '1nu8', '1nu8', '1nu8', '1nu8', '1nu8', '1nu8', '1ocn', '1ocn', '1ocn', '1ofl', '1ogq', '1ogq', '1ogq', '1ogq', '1ogq', '1ogq', '1ogq', '1oh4', '1onq', '1onq', '1op5', '1op5', '1op5', '1op5', '1op5', '1op5', '1op5', '1op5', '1ow0', '1ow0', '1ow0', '1ow0', '1pja', '1qon', '1qon', '1r9m', '1r9m', '1rem', '1rey', '1rez', '1sr5', '1sr5', '1sr5', '1sz6', '1sz6', '1sz6', '1sz6', '1sz6', '1sz6', '1sz6', '1sz6', '1u65']
##    print len(pdbs), len(set(pdbs))
##    stop

    for i in range(len(pdbs)):

        if i not in range(4*9000,5*9000):
            continue

        pdb = pdbs[i]
        if pdb in [
            ## missing atoms but mentioned in REMARK3
            '1a65',
            ## missing atoms or connected
            '1abr','1c3m','1cap','1cgu','1dot','1dva','1e05','1en2','1enm',
            '1ex1','2fyd','1js8','2msb',
            ## O1 missing
            '1bgc','1l1y','1l2a','1twx','1umz','1lof','1ven','1nma','1nmb','1nmc',
            ## O1 missing, NAG
            '1w1x','1w20','1w21','1rvx','1rvz','1zu8','1sbd','1sbe','1sbf',
            '1sgi','2alu','2aos','2b31','2c10','2c11','2c9a','1mx9','2cwg',
            '2dp8','2f83','2fa7','2g8z','2j8g','2ibx','2nwj','2o92','2ocv',
            '2po6','2sba','1fza','1fzb','1fze','2hpa',
            ## O1 missing, GAL/GLA/GLB
            '1lti','1tfm','1rvt','2c5c','2chb',
            ## O1 missing, GLC/BGC
            '1lax','1vfk','1yf8','1ryd','2aer','2bis','1mfu','1p2d','2fir',
            '2gjp','2j0y','2i5p','2osy',
            ## O1 missing, MAN/BMA
            '1zlw','2cml','2gn3','2gn7','2j0g','2j2p','2j3u','2man','3man',
            '2qwi','1bji',
            ## O1 missing, FUC
            '3kmb','4kmb',
            ## atoms other than O1 missing
            '1lfg','1uyx','1k3i','1lu1','1lzg','1muy','2dtx','1ovs','1q8v',
            '3chb','4hya',
            ]:
            continue
        fd = open('%s%s/pdb%s.ent' %(pdbpath,pdb[1:3],pdb),'r')
        lines = fd.readlines()
        fd.close()

        print i, len(pdbs), pdb

        d_FORMUL = {}
        d_HETATM = {}
        d_CONECT = {}
        l_iCodes = []
        l_MODRES = []
        d_LINK = {}

        for line in lines:

            record = line[:6].strip()

            if record == 'ATOM':
                d_HETATM = HETATM(line, d_HETATM)

            elif record == 'REMARK':
                continue

            elif record == 'HETATM':
                d_HETATM = HETATM(line, d_HETATM)

            elif record == 'FORMUL':
                d_FORMUL = FORMUL(line, d_FORMUL, d_hydrogen_equivalents)

            elif record == 'CONECT':
                d_CONECT, d_LINK = CONECT(line, d_FORMUL, d_HETATM, d_CONECT, d_LINK)

            elif record == 'MODRES':
                res_name = line[12:15].strip()
                l_MODRES += [res_name]

            elif record == 'HET':
                iCode = line[17]
                if iCode != ' ':
                    l_iCodes += [iCode]

##            elif record == 'LINK':
##                atom_name1 = line[12:16].strip()
##                res_name1 = line[17:20].strip()
##                chain1 = line[21]
##                res_no1 = int(line[22:26])
##                iCode1 = line[26]
##                atom_name2 = line[42:46].strip()
##                res_name2 = line[47:50].strip()
##                chain2 = line[51]
##                res_no2 = int(line[52:56])
##                iCode2 = line[56]
##                d_two = {
##                    1:{'atom_name':atom_name1,'res_name':res_name1,'chain':chain1,'res_no':res_no1,'iCode':iCode1},
##                    2:{'atom_name':atom_name2,'res_name':res_name2,'chain':chain2,'res_no':res_no2,'iCode':iCode2},
##                    }
##                for key in d_two.keys():
##                    atom_name = d_two[key]['atom_name']
##                    res_name = d_two[key]['res_name']
##                    chain = d_two[key]['chain']
##                    res_no = d_two[key]['res_no']
##                    iCode = d_two[key]['iCode']
##                    if not res_name in d_LINK.keys():
##                        d_LINK[res_name] = {}
##                    if not chain in d_LINK[res_name].keys():
##                        d_LINK[res_name][chain] = {}
##                    if not res_no in d_LINK[res_name][chain].keys():
##                        d_LINK[res_name][chain][res_no] = {}
##                    if not iCode in d_LINK[res_name][chain][res_no].keys():
##                        d_LINK[res_name][chain][res_no][iCode] = 0
##                    d_LINK[res_name][chain][res_no][iCode] += 1

                        
        for hetID in d_FORMUL.keys():
            if hetID not in d_CONECT.keys():
                continue
            if hetID in l_MODRES:
                continue
            if hetID not in ['GLC','AGC','BGC','G1P','G6P','BG6','G6Q','GAL','GLA','GLB','FUC','FUL','BGP','MAN','BMA','M1P','M6P','NAG','NBG','16G','NGA','FRU','F6P','SUC','LAT','LBT','MAL','TRE','CBI']:
                continue
            for chain in d_CONECT[hetID].keys():
                for res_no in d_CONECT[hetID][chain].keys():
                    for iCode in d_CONECT[hetID][chain][res_no].keys():
                        d_CONECT_elements = {'C':0,'N':0,'O':0}
                        for atom_name in d_CONECT[hetID][chain][res_no][iCode].keys():
                            if atom_name[0] in d_CONECT_elements.keys():
                                d_CONECT_elements[atom_name[0]] += 1

                        for element in ['C','N','O']:
                            if not element in d_FORMUL[hetID].keys():
                                continue
                            if d_CONECT_elements[element] != d_FORMUL[hetID][element]:

                                ##
                                ## count carbon connections if bound to anything but itself
                                ##
                                n_carbon_connections = 0
                                if chain in d_LINK.keys() and res_no in d_LINK[chain].keys() and iCode in d_LINK[chain][res_no].keys():
                                    for atom_name in d_LINK[chain][res_no][iCode].keys():
                                        if atom_name[0] == 'C':
                                            n_carbon_connections += 1

                                ##
                                ## incorrect connection
                                ##
                                if d_CONECT_elements[element] + n_carbon_connections != d_FORMUL[hetID][element]:
                                    atoms = d_CONECT[hetID][chain][res_no][iCode].keys()
                                    atoms.sort()
                                else:
                                    continue

                                ##
                                ## missing connection or atom
                                ##
##                                if chain not in d_LINK.keys() or res_no not in d_LINK[chain].keys() or iCode not in d_LINK[chain][res_no].keys():
                                if 'O1' not in d_CONECT[hetID][chain][res_no][iCode].keys():
                                    print pdb, hetID, chain, res_no, iCode
                                    print atoms
                                    print d_FORMUL[hetID]
                                    atom_no1 = d_CONECT[hetID][chain][res_no][iCode]['C1']['atom_no']
                                    coord1 = d_HETATM[atom_no1]['coordinate']
                                    found1 = False ## temp!!!
                                    for atom_no2 in d_HETATM.keys():
                                        if d_HETATM[atom_no2] == 'water':
                                            continue
                                        if (
                                            ## NAG N-glycosylation
                                            ## use ND2 *and* OD1 for ASN in case crystallographer assigned OD1,ND2 incorrectly
                                            (hetID == 'NAG' and d_HETATM[atom_no2]['res_name'] == 'ASN' and (d_HETATM[atom_no2]['atom_name'] == 'ND2' or d_HETATM[atom_no2]['atom_name'] == 'OD1')) or
                                            ## MAN/FUC O-glycosylation
                                            (hetID in ['MAN','FUC'] and d_HETATM[atom_no2]['res_name'] == 'SER' and d_HETATM[atom_no2]['atom_name'] == 'OG') or
                                            (hetID in ['MAN','FUC'] and d_HETATM[atom_no2]['res_name'] == 'THR' and d_HETATM[atom_no2]['atom_name'] == 'OG1') or
                                            ## MAN/GAL C-glycosylation
                                            (hetID == 'MAN' and d_HETATM[atom_no2]['res_name'] == 'TRP' and d_HETATM[atom_no2]['atom_name'] == 'CD1') or
                                            ## hetero connection
                                            (d_HETATM[atom_no2]['res_name'] not in ['GLY','ALA','VAL','LEU','ILE','SER','THR','GLU','ASP','ASN','GLN','ARG','LYS','TYR','HIS','TRP','CYS','PHE','MET','PRO'] and
                                             d_HETATM[atom_no2]['atom_name'][0] in ['O','S'] and not
                                             (
                                                 d_HETATM[atom_no2]['chain'] == chain and
                                                 d_HETATM[atom_no2]['res_no'] == res_no and
                                                 d_HETATM[atom_no2]['iCode'] == iCode
                                                 )
                                             )
                                            ):
                                            coord2 = d_HETATM[atom_no2]['coordinate']
                                            distance = math.sqrt(sum((coord1-coord2)**2))

                                            if distance < 3.5:

                                                found2 = False

                                                ## not C-glycosylation
                                                if d_HETATM[atom_no2]['res_name'] != 'TRP':

                                                    for atom_no3 in d_HETATM.keys():

                                                        if d_HETATM[atom_no3] == 'water':
                                                            continue
                                                        if (
                                                            d_HETATM[atom_no3]['chain'] == d_HETATM[atom_no2]['chain'] and
                                                            d_HETATM[atom_no3]['res_no'] == d_HETATM[atom_no2]['res_no'] and
                                                            d_HETATM[atom_no3]['iCode'] == d_HETATM[atom_no2]['iCode'] and
                                                            (
                                                                ## N-glycosylation
                                                                (d_HETATM[atom_no2]['res_name'] == 'ASN' and d_HETATM[atom_no3]['atom_name'] == 'CG') or
                                                                ## O-glycosylation
                                                                (d_HETATM[atom_no2]['res_name'] == 'SER' and d_HETATM[atom_no3]['atom_name'] == 'CB') or
                                                                (d_HETATM[atom_no2]['res_name'] == 'THR' and d_HETATM[atom_no3]['atom_name'] == 'CB') or
                                                                ## hetero connection
                                                                (d_HETATM[atom_no2]['res_name'] not in ['GLY','ALA','VAL','LEU','ILE','SER','THR','GLU','ASP','ASN','GLN','ARG','LYS','TYR','HIS','TRP','CYS','PHE','MET','PRO'] and
                                                                 d_HETATM[atom_no3]['atom_name'] == 'C'+d_HETATM[atom_no2]['atom_name'][1:])
                                                                )
                                                            ):
                                                            found2 = True
                                                            break
                                                    if found2 == False:
                                                        print d_HETATM[atom_no2]
                                                        stop
                                                    ## intramolecule C1-O1 vector
                                                    v1 = d_HETATM[atom_no1]['coordinate']-d_HETATM[atom_no2]['coordinate']
                                                    ## intermolecular O1-C2 vector
                                                    v2 = d_HETATM[atom_no3]['coordinate']-d_HETATM[atom_no2]['coordinate']

                                                    dotproduct = sum(v1*v2)
                                                    cosangle = dotproduct/math.sqrt(sum(v1**2)*sum(v2**2))
                                                    angle = math.acos(cosangle)*180/math.pi

                                                else:

                                                    found2 = True
                                                    angle = 180.

                                                found1 = True ## temp!!!
                                                line = '%s %5.3f %5.1f %5i %5i %3s %3s %4s %4s %1s %1s %4i %4i %1s %1s' %(
                                                    pdb, distance, angle, atom_no1, atom_no2,
                                                    d_HETATM[atom_no1]['res_name'], d_HETATM[atom_no2]['res_name'],
                                                    d_HETATM[atom_no1]['atom_name'], d_HETATM[atom_no2]['atom_name'],
                                                    d_HETATM[atom_no1]['chain'], d_HETATM[atom_no2]['chain'],
                                                    d_HETATM[atom_no1]['res_no'], d_HETATM[atom_no2]['res_no'],
                                                    d_HETATM[atom_no1]['iCode'], d_HETATM[atom_no2]['iCode'],
                                                    )
                                                lines_out += [line+'\n']
                                    if found1 == True: ## temp!!!
                                        continue
                                    else:
                                        print pdb, hetID, chain, res_no, iCode
                                        print atoms
                                        print d_FORMUL[hetID]
                                        stop_perhaps_missing_connection_or_atom


                                ##
                                ##
                                ##
                                print hetID, chain, res_no, iCode
                                print d_CONECT[hetID][chain][res_no][iCode]
                                print atoms
                                l_intra_atom_nos = []
                                for atom_name in d_CONECT[hetID][chain][res_no][iCode].keys():
                                    l_intra_atom_nos += d_CONECT[hetID][chain][res_no][iCode][atom_name]['atom_nos']
                                l_intra_atom_nos = list(set(l_intra_atom_nos))
                                d_minsqdistance = {'sqdist':'N/A'}
                                for intra_atom_no in l_intra_atom_nos:
                                    ## bond from carbon (in molecule with missing oxygen)
                                    if d_HETATM[intra_atom_no]['atom_name'][0] != 'C':
                                        continue
                                    intra_coordinate = d_HETATM[intra_atom_no]['coordinate']
                                    for inter_atom_no in d_HETATM.keys():
                                        if inter_atom_no in l_intra_atom_nos:
                                            continue
                                        if d_HETATM[inter_atom_no] == 'water':
                                            continue
                                        ## bond to oxygen
                                        if d_HETATM[inter_atom_no]['atom_name'][0] != 'O':
                                            continue
                                        hetID_inter = d_HETATM[inter_atom_no]['res_name']
                                        if hetID_inter not in ['GLC','AGC','BGC','G1P','G6P','BG6','G6Q','GAL','GLA','GLB','FUC','FUL','BGP','MAN','BMA','M1P','M6P','NAG','NBG','16G','NGA','FRU','F6P','SUC','LAT','LBT','MAL','TRE','CBI']:
                                            continue
                                        inter_coordinate = d_HETATM[inter_atom_no]['coordinate']
                                        sqdist = sum((inter_coordinate-intra_coordinate)**2)
                                        if sqdist < d_minsqdistance['sqdist']:
                                            d_minsqdistance['intra_atomno'] = intra_atom_no
                                            d_minsqdistance['inter_atomno'] = inter_atom_no
                                            d_minsqdistance['sqdist'] = sqdist
                                print d_minsqdistance
                                intra_atomno = d_minsqdistance['intra_atomno']
                                inter_atomno1 = d_minsqdistance['inter_atomno']
                                distance = math.sqrt(d_minsqdistance['sqdist'])
                                found = False
                                for inter_atomno2 in d_HETATM.keys():
                                    if (
                                        d_HETATM[inter_atomno2]['chain'] == d_HETATM[inter_atomno1]['chain'] and
                                        d_HETATM[inter_atomno2]['res_no'] == d_HETATM[inter_atomno1]['res_no'] and
                                        d_HETATM[inter_atomno2]['iCode'] == d_HETATM[inter_atomno1]['iCode'] and
                                        d_HETATM[inter_atomno2]['atom_name'] == 'C'+d_HETATM[inter_atomno1]['atom_name'][1:]
                                        ):
                                        found = True
                                        break
                                if found == False:
                                    stop
                                ## intramolecule C1-O1 vector
                                v1 = d_HETATM[inter_atomno2]['coordinate']-d_HETATM[inter_atomno1]['coordinate']
                                ## intermolecular O1-C2 vector
                                v2 = d_HETATM[inter_atomno1]['coordinate']-d_HETATM[intra_atomno]['coordinate']

                                dotproduct = sum(v1*v2)
                                cosangle = dotproduct/math.sqrt(sum(v1**2)*sum(v2**2))
                                angle = 180-math.acos(cosangle)*180/math.pi
                                print distance, angle
                                print intra_atomno, d_HETATM[intra_atomno]['atom_name']
                                print inter_atomno1, d_HETATM[inter_atomno1]['atom_name']
                                print inter_atomno2, d_HETATM[inter_atomno2]['atom_name']
                                print pdb, hetID, chain, res_no, iCode
                                print d_HETATM[inter_atomno1]['res_name'], d_HETATM[inter_atomno1]['chain'], d_HETATM[inter_atomno1]['res_no'], d_HETATM[inter_atomno1]['iCode']
                                error

                                ##
                                ##
                                ##
                                for atom_name in d_LINK[chain][res_no][iCode].keys():
                                    if len(d_LINK[chain][res_no][iCode][atom_name]['connections']) > 1:
                                        print d_LINK[chain][res_no][iCode][atom_name]
                                        print hetID, chain, res_no, iCode, atom_name, d_LINK[chain][res_no][iCode][atom_name]
                                        atom_no1 = d_LINK[chain][res_no][iCode][atom_name]['atom_no']
                                        for atom_no2 in d_LINK[chain][res_no][iCode][atom_name]['connections']:
                                            coord1 = d_HETATM[atom_no1]['coordinate']
                                            coord2 = d_HETATM[atom_no2]['coordinate']
                                            distance = math.sqrt(sum((coord1-coord2)**2))
                                            print atom_no2, distance, d_HETATM[atom_no2]['atom_name'], d_HETATM[atom_no2]['altloc']
                                        raise 'atom %s incorrectly connected' %atom_no1

                            
##        if d_CONECT != {}:
##            for hetID in d_FORMUL.keys():
##                d_elements = d_FORMUL[hetID]
##                hydrogens_FORMUL = expected_hydrogens(d_hydrogen_equivalents, d_elements)
##                for chain in d_CONECT[hetID].keys():
##                    for res_no in d_CONECT[hetID][chain].keys():
##                        hydrogens_HETATM = 0
##                        for iCode in d_CONECT[hetID][chain][res_no].keys():
##                            if iCode != ' ':
##                                notexpected
##                            for atom_name in d_CONECT[hetID][chain][res_no][iCode]:
##                                connections = d_CONECT[hetID][chain][res_no][iCode][atom_name]
##                                print atom_name, d_valence[atom_name[0]]-len(connections)
##                                hydrogens_HETATM += d_valence[atom_name[0]]-len(connections)
##                        print hetID, chain, res_no, hydrogens_HETATM, hydrogens_FORMUL, d_elements['H']
##                        if hydrogens_HETATM != d_elements['H']:
##                            print hetID, chain, res_no
##                            print d_CONECT[hetID][chain][res_no]
##                            print d_elements
##                            stop

    fd = open('output.txt','a')
    fd.writelines(lines_out)
    fd.close()

    return


def CONECT(line, d_FORMUL, d_HETATM, d_CONECT, d_LINK):

    atom_no = int(line[6:11])
    if d_HETATM[atom_no] == False:
        print atom_no
        stop
        return d_CONECT, d_LINK
    if d_HETATM[atom_no] == 'water':
        return d_CONECT, d_LINK

    hetID1 = d_HETATM[atom_no]['res_name']
    FORMULhetIDs = d_FORMUL.keys()
    if hetID1 in FORMULhetIDs:

        if hetID1 not in d_CONECT.keys():
            d_CONECT[hetID1] = {}

        ##
        ## parse atom nos
        ##
        atom_nos = []
        for j in range(6,31,5):
            if line[j:j+5] != '     ':
                atom_nos += [int(line[j:j+5])]

        ##
        ## translate atom nos to hetIDs
        ##
        atom_no1 = atom_nos[0]
        if hetID1 != d_HETATM[atom_no]['res_name']:
            notexpected
        chain1 = d_HETATM[atom_no]['chain']
        res_no1 = d_HETATM[atom_no]['res_no']
        iCode1 = d_HETATM[atom_no]['iCode']
        altloc1 = d_HETATM[atom_no]['altloc']
        atom_name1 = d_HETATM[atom_no]['atom_name']

        if atom_name1[0] == 'H':
            return d_CONECT, d_LINK

        if chain1 not in d_CONECT[hetID1].keys():
            d_CONECT[hetID1][chain1] = {}
        if res_no1 not in d_CONECT[hetID1][chain1].keys():
            d_CONECT[hetID1][chain1][res_no1] = {}
        if iCode1 not in d_CONECT[hetID1][chain1][res_no1].keys():
            d_CONECT[hetID1][chain1][res_no1][iCode1] = {}
        if atom_name1 not in d_CONECT[hetID1][chain1][res_no1][iCode1].keys():
            d_CONECT[hetID1][chain1][res_no1][iCode1][atom_name1] = {}

        l_atom_nos = []
        for atom_no2 in atom_nos[1:]:
            chain2 = d_HETATM[atom_no2]['chain']
            res_no2 = d_HETATM[atom_no2]['res_no']
            iCode2 = d_HETATM[atom_no2]['iCode']
            altloc2 = d_HETATM[atom_no2]['altloc']
            atom_name2 = d_HETATM[atom_no2]['atom_name']
            hetID2 = d_HETATM[atom_no2]['res_name']
            bond = 'intra'
            if (
                chain1 !=  chain2 or
                res_no1 != res_no2 or
                iCode1 != iCode2
                ):
                bond = 'inter'
                if not chain1 in d_LINK.keys():
                    d_LINK[chain1] = {}
                if not chain2 in d_LINK.keys():
                    d_LINK[chain2] = {}
                if not res_no1 in d_LINK[chain1].keys():
                    d_LINK[chain1][res_no1] = {}
                if not res_no2 in d_LINK[chain2].keys():
                    d_LINK[chain2][res_no2] = {}
                if not iCode1 in d_LINK[chain1][res_no1].keys():
                    d_LINK[chain1][res_no1][iCode1] = {}
                if not iCode2 in d_LINK[chain2][res_no2].keys():
                    d_LINK[chain2][res_no2][iCode2] = {}
                if not atom_name1 in d_LINK[chain1][res_no1][iCode1].keys():
                    d_LINK[chain1][res_no1][iCode1][atom_name1] = {'atom_no':atom_no1,'connections':set()}
                if not atom_name2 in d_LINK[chain2][res_no2][iCode2].keys():
                    d_LINK[chain2][res_no2][iCode2][atom_name2] = {'atom_no':atom_no1,'connections':set()}
                d_LINK[chain1][res_no1][iCode1][atom_name1]['connections'] |= set([atom_no2])
                d_LINK[chain2][res_no2][iCode2][atom_name2]['connections'] |= set([atom_no1])
                
            if d_HETATM[atom_no]['atom_name'][0] == 'H':
                continue
            l_atom_nos += [atom_no2]

##        if len(atom_nos) <= 1:
##            return d_CONECT

        ## append set of atom_nos (use set to avoid alternate locations)
        d_CONECT[hetID1][chain1][res_no1][iCode1][atom_name1]['atom_nos'] = list(set(l_atom_nos))
        d_CONECT[hetID1][chain1][res_no1][iCode1][atom_name1]['atom_no'] = atom_no1

    return d_CONECT, d_LINK


def HETATM(line, d_HETATM):

    import Numeric

    res_name = line[17:20].strip()
    atom_no = int(line[6:11])

    if res_name in ['HOH','H20','DOD','D2O']:
        d_HETATM[atom_no] = 'water'
        return d_HETATM

    atom_name = line[12:16].strip()
    altloc = line[16]
    res_name = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]
    if iCode != ' ' and line[:6] == 'REMARK':
        notexpected
    atom_x = float(line[30:38])
    atom_y = float(line[38:46])
    atom_z = float(line[46:54])
    coordinate = Numeric.array([atom_x, atom_y, atom_z])

    d_line = {
        'atom_no':atom_no,
        'atom_name':atom_name,
        'res_name':res_name,
        'chain':chain,
        'res_no':res_no,
        'iCode':iCode,
        'altloc':altloc,
        'coordinate':coordinate
        }

    d_HETATM[atom_no] = d_line

    return d_HETATM


def FORMUL(line, d_FORMUL, d_hydrogen_equivalents):

    hetID = line[12:15].strip()
    if hetID == 'HOH':
        return d_FORMUL
    continuation = line[16:18]
    if continuation != '  ':
        print line
        stop
    formula = line[19:70]
    if '(' in formula:
        index1 = formula.index('(')+1
        index2 = formula.index(')')
        formula = formula[index1:index2]
    formula = formula.split()

    d_elements = {}
    for atoms in formula:
        element = atoms
        for i in range(10):
            element = element.replace(str(i),'')
        count = atoms.replace(element,'')
        if count == '':
            count = 1
        d_elements[element] = int(count)

    ## other elements than carbon,hydrogen,oxygen,halogen
    if len(set(d_elements.keys()) - set(d_hydrogen_equivalents.keys())) > 0:
        return d_FORMUL
    ## exclude hetero compounds without carbon,hydrogen
    for element in ['C','H']:
        if not element in d_elements.keys():
            return d_FORMUL
    ## more than 8 carbon atoms (to include NAG and other acetyl substituted hexoses)
    if d_elements['C'] > 8:
        return d_FORMUL

    hydrogens = expected_hydrogens(d_hydrogen_equivalents, d_elements)

    if hydrogens != d_elements['H']:
        d_FORMUL[hetID] = d_elements

    return d_FORMUL


def expected_hydrogens(d_hydrogen_equivalents, d_elements):

    hydrogens = 2

    for element in d_hydrogen_equivalents.keys():
        if element in d_elements.keys():
            hydrogens += d_elements[element]*d_hydrogen_equivalents[element]

    return hydrogens


if __name__ == '__main__':
    main()

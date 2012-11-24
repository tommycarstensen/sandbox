#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

## replace [1] with [biounit] in the case of PDB transformations

## os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ %s' %(self.path_pdb))

## is there a risk PISA coordinates are formed from transformation instead of existing coordinates?! 1a99???

import sys, Numeric, biounit, os

class comparison:

    def main(self):

        '''Script for comparison of biological units from the PDB and from MSD-PISA'''

##        self.rsync()
##        self.gunzip()
        
        fd = open('pdbbind_2007.txt','r')
        s = fd.read()
        fd.close()
        l_pdbs = s.split()

        d_pdbs = {}

        for pdb in l_pdbs:
            if not pdb[1:3] in d_pdbs.keys():
                d_pdbs[pdb[1:3]] = [pdb]
            else:
                d_pdbs[pdb[1:3]] += [pdb]

        l_subs = d_pdbs.keys()
        l_subs.sort()
        l_pdbs = []
        for sub in l_subs:
            for pdb in d_pdbs[sub]:
                l_pdbs += [pdb]
            
        for pdb in l_pdbs:

            if l_pdbs.index(pdb)+1 < int(sys.argv[1]):
                continue
            print pdb, l_pdbs.index(pdb)+1, len(l_pdbs)

            if pdb in [
##                '1fiv','1hef','1heg','1e5a', ## ligand overlap upon transformation, two alternative ligand binding conformations?!' e.g. 1bm7!!!
                ## biounits should be identical, but different number of ligands in each PDB transformed biounit (water atoms often transformed incorrectly!!)
##                '1osv','1oxn','1oxq','2pcp','1u9l','1igj','1jq9','1jq8','1lpk','1lpg','1lpz','2a4m',
##                '1c1r', ## grey region
                '2a4m', ## ligand not in PDB
                '2a5b','2a5c','2a8g', ## some of the ligands not transformed in PDB
                '1a69', ## altloc used for alternative temperature factors (but identical coordinates!)
                '1a8i', ## v2 atom names in biounit

                ##
                ## multiple biounits
                ##
##                '1fpu',
                ## identical multimer in PISA and PDB
##                '1a4k','1a94','1pkx','1wc1','1fj4','1h1p','1cgl','1dqx','1eix','1fkn','1fl3','1h1s','1hsh','1i7z','2cht','1uw6','1jqy','1los','1lrh','1m4h','1a08','1is0','1mh5','1mjj','1njj','1p1q','1q4k','3tmk','1umw','1uv6','1uz8',
                ## tommy error
##                '2dqt', ## dimer in PISA and PDB

                ##
                ## other PISA multimers and/or interfaces are stable in solution
                ##
##                '1qhc','1fq5','1it6','1jn4','2jxr','1yei','1tuf','1e2k','1e2n','1oe7','1e2l','1nms','1afk','1fch','1fkf','1kc7','6std','7std','1fzj','1fzk','1slg','1sle','1jyq','1kyv','1o9d','1oe8','1os0','1p19','1qca','5std','1tyr','1ugp','1v48','1vfn','1vpo','1vwl','1aqc','1ghy','2izl',

                ##
                ## same size multimers, but different interfaces
                ##
                '1adl','1af2','2ans',##'1loq','1lyx','1kll','1trd','1w72','1b55','1oko','1lyb','3pck','3pcj',

                ##
                ## different multimers
                ##

                ## dimer in PISA, monomer in PDB
                '2ayr','1b11',##'1oar','1s39','1udt','1p6e','1caq','1rd4','1p9p','1q63','1qi0','1gpk','1lf9','1k1y','1k4g','1oim','1b8y','1gpn','1ow4','1h6h','1lee','1qy2','1q4w','1w3j','1c5s','1ciz','1d7j','1dy4','1f3e','1ghz','1imx','1j17','1jt1','1k4h','1kpm','1l2s','1lnm','1m13','1m48','1n2v','1njs','1nw7','1nw5','1oif','1p28','1q65','1q66','1q91','1qft','1qy1','1r5y','1s38','1sqn','1sw1','1uho','1uj6','1uj5','1uz1','1wm1','1xzx','5yas','1d7i',
                ## multimer in PISA, monomer in PDB
##                '1b8o','1b8n','1j4r','1g7v','1lf2','830c','1jn2','1fv0','4tmk','5tmp','2usn', '1usn',

                ## correct multimer might be in "grey area" (deltaG_dissociation ~< 0)
                ## PQS interfaces might be different from PDB interaces
                '1atl',##'1b42','1bky','1bra','3mag','3mct', ## monomer in PISA, dimer in PDB
##                '2csn','1iup','1gz9','1igb','1ii5','1p1n', ## monomer in PISA, dimer in PDB, dimer in PQS
##                '1kdk','1lhw', ## monomer in PISA, dimer in PDB, hexamer in PQS
##                '1f8d', ## monomer in PISA, tetramer in PDB, tetramer in PQS, tetramer in PISA upon removal of ligands
##                '1bm7','1n51','2tmn','4tmn','5tmn','1wht', ## dimer in PISA, tetramer in PDB, tetramer in PQS
                '1a1b','1a1c','1a1e', ## dimer in PISA, tetramer in PDB
##                '1m5w', ## dimer in PISA, octamer in PDB, octamer in PQS
##                '1ftm', ## trimer in PISA, monomer in PDB, dimer in PQS
##                '1awi', ## monomer/dimer in PISA, trimer in PDB, trimer in PQS
                '1a99', ## tetramer in PISA, dimer in PDB
                
                
                ##
                ## ligand positions
                ##
                ## identical multimers in PISA and PDB, but different ligand/ion/sugar position(s)
                '2aj8','2aoc','2aod','2aoe','2aog','1apv','1apw',##'2cgr','3gss','1elr','1elb','1hyo','1gvx','1gyx','1jet','2gss','1gyy','1gvu','1qkb','1ofz','1b9j','1jao','1jeu','1jev','10gs','11gs','1kui','1kuk','1kug','1obx','1ogx','1px4','1qka','1ur9','1e6s','1e6q','1e70',
                ## ACY not in PISA
                '2aac','2avm',
                ## MG positions
                '1af6',
                ## SO4 positions
                '2avo','2avs',
                ## GOL positions
                '2avq',
                ## CL positions
                '2avv',
                ## U1 positions
                '1b05','1b0h','1b2h','1b4h',
                ## IUM positions
                '1b1h','1b32','1b3f','1b3g','1b3h','1b3l','1b40','1b46','1b4z','1b51','1b52',

##                '1h22','1h23', ## acetylcholine esterase
                ## acetylcholine esterase
                ## dimer
                ## 4 helix bundle interface
                ## 0.490nm between LYS530NZ and ASP365ODD2
                ## 0.263nm between LYS530NZ and ASP369ODD2
                ## hydrophobic core involving LEU366,LEU373,PHE527,LEU531
                ## other AChE structures (1j06,1j07,1n5r,1n5m) have similar dimer interfaces with ligand in between interfaces
                ## PDB molecule of the month states it is a dimer, but no details about the dimer interface
                ]:
                continue


            ##
            ##
            ##
            d_transformations_PISA = biounit.biounit().parse_pisa_multimers(pdb)
            if d_transformations_PISA == {}:
                continue ## temporary!!!

            d_coordinates_PDB = {}
            d_lines_PDB = {}
            for bm in range(1,10000):
                if bm == 1 or os.path.isfile('/oxygenase_local/data/biounit/%s/%s.pdb%i' %(pdb[1:3],pdb,bm)):
                    print bm
                    fd = open('/oxygenase_local/data/biounit/%s/%s.pdb%i' %(pdb[1:3],pdb,bm),'r')
                    lines = fd.readlines()
                    fd.close()
                    d_coordinates_PDB[bm], d_lines_PDB[bm] = self.parse_coordinates(lines)
                    continue
                else:
                    if bm > 3:
                        print bm
                        stop
                    break

            l_biomolecules = d_coordinates_PDB.keys()
            l_assemblies = d_transformations_PISA.keys()

            for bm in l_biomolecules:
                biounits_identical = False
                chains_identical = False
                interfaces_identical = False
                l_chains_PDB = d_coordinates_PDB[bm]['ATOM'].keys()
                for assembly in l_assemblies:
                    print bm, assembly
                    l_chains_PISA = []
                    for chain in d_transformations_PISA[assembly]['chains'].keys():
                        if len(chain) == 1 and chain != '-':
                            l_chains_PISA += [chain]
##                    print d_transformations_PISA[assembly]['chains'].keys()
##                    print d_chains_PDB['ATOM'][bm]+d_chains_PDB['HETATM'][bm]

                    ## different chains (different IDs)
                    if len(set(l_chains_PISA) ^ set(l_chains_PDB)) > 0:
                        print 'different chains', l_chains_PDB, l_chains_PISA 
                        continue
                    else:
                        d_lines_PISA, d_coordinates_PISA = biounit.biounit().parse_pdb_coordinates(pdb, d_transformations_PISA, assembly)
                        ## different chains (different number of transformations)
                        if not (
                            len(d_lines_PISA[assembly]['ATOM']) % len(d_lines_PDB[bm]['ATOM']) == 0 and
                            len(d_lines_PDB[bm]['ATOM']) % len(d_lines_PISA[assembly]['ATOM'])== 0
                            ):
                            if (
                                len(d_lines_PISA[assembly]['ATOM']) % len(d_lines_PDB[bm]['ATOM']) != 0 and
                                len(d_lines_PDB[bm]['ATOM']) % len(d_lines_PISA[assembly]['ATOM']) != 0
                                ):
                                stop
                            chains_identical = False
                            print 'different chains', len(d_lines_PDB[bm]['ATOM']), len(d_lines_PISA[assembly]['ATOM'])
                            continue
                        ## identical chains (identical IDs, identical number of transformations)
                        else:
                            chains_identical = True

                    ATOM_identical = self.identical_d_coordinates('ATOM',d_coordinates_PISA,d_coordinates_PDB,assembly,bm,)
                    HETATM_identical = self.identical_d_coordinates('HETATM',d_coordinates_PISA,d_coordinates_PDB,assembly,bm,)
                    print ATOM_identical, HETATM_identical

                    if chains_identical == True and ATOM_identical == True and HETATM_identical == True:
                        biounits_identical = True
                        print bm,assembly,'identical'
                    elif chains_identical == True and ATOM_identical == False and HETATM_identical == True:
                        print bm,assembly,'different interfaces'
                        continue
                    elif chains_identical == True and ATOM_identical == True and  HETATM_identical == False:
                        interfaces_identical = True
                        print bm,assembly,'different ligands'
                                   
##                    if (
##                        len(set(d_lines_PISA[assembly]['ATOM'])^set(d_lines_PDB[bm]['ATOM'])) == 0 and
##                        len(set(d_lines_PISA[assembly]['HETATM'])^set(d_lines_PDB['HETATM'][bm])) == 0
##                        ):
##                        biounits_identical = True
##                        print assembly,bm, 'identical'
##                        break
##                    elif len(set(d_lines_PISA[assembly]['ATOM'])^set(d_lines_PDB[bm]['ATOM'])) != 0:
##                        set_PISA_ATOM = set(d_lines_PISA[assembly]['ATOM'])-set(d_lines_PDB[bm]['ATOM'])
##                        set_PDB_ATOM = set(d_lines_PDB[bm]['ATOM'])-set(d_lines_PISA[assembly]['ATOM'])
##                        print len(set_PISA_ATOM), len(set_PDB_ATOM)
##                        print bm, assembly, 'different interfaces'
##                        if assembly == 4:
##                            a =  list(set_PISA_ATOM)
##                            b = list(set_PDB_ATOM)
##                            a.sort()
##                            b.sort()
##                            print a[:10]
##                            print b[:10]
##                            stop
##                        continue
##                    else:
##                        interfaces_identical = True
##                        set_PISA_ATOM = set(d_lines_PISA[assembly]['ATOM'])-set(d_lines_PDB[bm]['ATOM'])
##                        set_PDB_ATOM = set(d_lines_PDB[bm]['ATOM'])-set(d_lines_PISA[assembly]['ATOM'])
##                        set_PISA_HETATM = set(d_lines_PISA[assembly]['HETATM'])-set(d_lines_PDB['HETATM'][bm])
##                        set_PDB_HETATM = set(d_lines_PDB['HETATM'][bm])-set(d_lines_PISA[assembly]['HETATM'])
##                        a = list(set_PISA_HETATM)
##                        a.sort()
##                        b = list(set_PDB_HETATM)
##                        b.sort()
##                        print a
##                        print b
##                        print len(set_PISA_ATOM), len(set_PDB_ATOM)
##                        print len(set_PISA_HETATM), len(set_PDB_HETATM)
##                        print assembly, bm, 'different ligands'
##                        continue
                if biounits_identical == False:
                    print pdb, l_pdbs.index(pdb)+1
                    if chains_identical == True and interfaces_identical == True:
                        stop_ligand_differences
                    elif chains_identical == True and interfaces_identical == False:  
                        stop_interfaces_different
                    elif chains_identical == False and interfaces_identical == False:
                        stop_different_multimers
                    else:
                        stop_not_expected
                elif biounits_identical == True:
                    continue
            continue


            ## monomer in asu and biou
            if d_transformations_REMARK350 == {} and d_transformations_PISA == {}:
                continue

            ## asu == biou in PDB, biou == asu in PISA
            if d_transformations_REMARK350 == {}:
                for assembly in d_transformations_PISA.keys():
                    for chain in d_transformations_PISA[assembly]['chains'].keys():
                        for molecule in d_transformations_PISA[assembly]['chains'][chain].keys():
                            if d_transformations_PISA[assembly]['chains'][chain][molecule]['r'] != Numeric.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]):
                                stop1
                            if d_transformations_PISA[assembly]['chains'][chain][molecule]['t'] != Numeric.array([0.,0.,0.]):
                                stop2
                continue

            ## biou=asu in PISA, biou!=asu in PDB
            if d_transformations_PISA == {}:
                for biou in d_transformations_REMARK350.keys():
                    chains = d_transformations_REMARK350[biou]['chains'].keys()
                    for chain in chains:
                        matrixnos = d_transformations_REMARK350[biou]['chains'][chain]
                        if len(matrixnos) != 1:
                            stop2
                        matrix = d_transformations_REMARK350[biou]['matrices'][list(matrixnos)[0]]
                        if matrix != [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]:
                            stop3
                continue

            biounits = d_transformations_REMARK350.keys()
            ## multimer in PISA and PDB
            if len(biounits) != 1 and pdb not in [
                ]: ## loop over biounits and replace [1] with [biounit] if this doesnt hold true!!!
                print d_transformations_PISA
                print d_transformations_REMARK350
                print pdb, biounits
                stop_multimer_in_PISA_and_PDB

##            print d_transformations_PISA
##            print d_transformations_REMARK350
##            print pdb

##            for assembly in d_transformations_PISA.keys():
##                size1 = 0
##                for chain_PISA in d_transformations_PISA[assembly]['chains'].keys():
##                    if len(chain_PISA) == 1:
##                        chain = chain_PISA
##                        if chain == '-':
##                            continue
##                    else:
##                        chain = chain_PISA[chain_PISA.index(']')+1]
##                        if chain == '-':
##                            continue
##                    molecules = d_transformations_PISA[assembly]['chains'][chain_PISA].keys()
##                    size1 += len(molecules)
##                print size1
##                break
##
##            size2 = 0
##            for chain in d_transformations_REMARK350[1]['chains'].keys():
##                matrixnos = d_transformations_REMARK350[1]['chains'][chain]
##                print chain, matrixnos
##                size2 += len(matrixnos)
##            print size2
##
##            if size2 > size1:
##                print d_transformations_REMARK350
##                stop_maybe_water_has_chain_id

            d_transformations = {'chains':{},'matrices':{}}
            for assembly in d_transformations_PISA.keys():
                matrices_identical = True
                chains_PISA = d_transformations_PISA[assembly]['chains'].keys()
                for chain_PISA in chains_PISA:

                    ## parse PISA matrix
                    molecules = d_transformations_PISA[assembly]['chains'][chain_PISA].keys()
                    for molecule in molecules:
                        r = d_transformations_PISA[assembly]['chains'][chain_PISA][molecule]['r']
                        t = d_transformations_PISA[assembly]['chains'][chain_PISA][molecule]['t']

                        ## convert PISA chain ID to default chain ID
                        if len(chain_PISA) == 1:
                            chain = chain_PISA
                            if chain == '-': ## temporary!!!
                                continue
                        else:
                            chain = chain_PISA[chain_PISA.index(']')+1]
                            if chain == '-': ## temporary!!!
                                continue

                        ## compare PISA and REMARK350 matrices
                        set_matrixnos = d_transformations_REMARK350[1]['chains'][chain]
                        for matrixno in set_matrixnos:
                            matrix_identical = True
                            matrix_REMARK350 = d_transformations_REMARK350[1]['matrices'][matrixno]
                            d_transformations['matrices'][matrixno] = matrix_REMARK350
                            for i in range(3):
                                if (
                                    round(float(matrix_REMARK350[i][0]),5) == round(r[i][0],5) and
                                    round(float(matrix_REMARK350[i][1]),5) == round(r[i][1],5) and
                                    round(float(matrix_REMARK350[i][2]),5) == round(r[i][2],5) and
                                    round(float(matrix_REMARK350[i][3]),5) == round(t[i],5)
                                    ):
                                    continue
                                else:
                                    matrix_identical = False
                                    break
                            ## continue loop over REMARK350 matrices
                            if matrix_identical == False:
                                continue
                            else:
                                if chain not in d_transformations['chains'].keys():
                                    d_transformations['chains'][chain] = set([matrixno])
                                else:
                                    d_transformations['chains'][chain] |= set([matrixno])
                                if matrix_identical == False:
                                    stop_temporary
                                matrix_identical = True
                                break
                        ## break loop over molecules
                        if matrix_identical == False:
                            for matrixno in set_matrixnos:
                                print d_transformations_REMARK350[1]['matrices'][matrixno]
                            print 'assembly', assembly
                            print 'molecule', molecule
                            print 'chain', chain_PISA
                            print d_transformations_PISA[assembly]['chains'][chain_PISA]
                            print float(matrix_REMARK350[i][0]) == round(r[i][0],6)
                            print float(matrix_REMARK350[i][1]) == round(r[i][1],6)
                            print float(matrix_REMARK350[i][2]) == round(r[i][2],6)
                            print float(matrix_REMARK350[i][3]), round(t[i],6)
                            if len(chain_PISA) == 1:
                                stop_multimer_difference
                            else:
                                stop_different_ligand_locations
                            matrices_identical = False
                            break
                    ## break loop over PISA chains
                    if matrices_identical == False:
                        break
                    if chain != '-':
                        if matrix_identical == False:
                            print assembly, molecule, chain_PISA
                            print matrices_identical
                            stop1
                ## continue loop over assemblies
                if matrices_identical == False:
                    continue
                if matrix_identical == False:
                    stop2
##            if matrix_identical == False:
##                stop3
            if d_transformations_PISA != {}:
                if matrices_identical == False:
                    print d_transformations_REMARK350[1]['matrices']
                    print d_transformations_PISA[assembly]['chains'][chain_PISA]
                    print assembly, molecule, chain_PISA
                    print d_transformations_PISA.keys()
                    stop4
                   
                            
                
            if d_transformations_REMARK350[1] != d_transformations:
##                if (
##                    d_transformations_REMARK350[1]['matrices'].keys() != 1 and
##                    d_transformations_REMARK350[1]['matrices'][1] != [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]
##                    ):
##                    print d_transformations_PISA
                    print d_transformations_REMARK350[1]
                    print d_transformations
                    print pdb
                    stop_PDB_larger_than_PISA

        return


    def identical_d_coordinates(self,record,d_coordinates_PISA,d_coordinates_PDB,assembly,bm,):

        ATOM_identical = True

        for chain in d_coordinates_PISA[assembly][record].keys():
            for res_no in d_coordinates_PISA[assembly][record][chain].keys():
                for iCode in d_coordinates_PISA[assembly][record][chain][res_no].keys():
                    for atom_name in d_coordinates_PISA[assembly][record][chain][res_no][iCode].keys():
                        l_coordinates_PISA = d_coordinates_PISA[assembly][record][chain][res_no][iCode][atom_name]
                        l_coordinates_PDB = d_coordinates_PDB[bm][record][chain][res_no][iCode][atom_name]
                        identical_PISA = self.identical_l_coordinates(l_coordinates_PISA,l_coordinates_PDB)
                        identical_PDB = self.identical_l_coordinates(l_coordinates_PDB,l_coordinates_PISA)
                        if identical_PISA == False or identical_PDB == False:
                            ATOM_identical = False
                            break
                    if ATOM_identical == False:
                        break
                if ATOM_identical == False:
                    break
            if ATOM_identical == False:
                break

        return ATOM_identical


    def identical_l_coordinates(self, l_coordinates_PISA,l_coordinates_PDB):

        ATOM_identical = True
        for c_PISA in l_coordinates_PISA:
            coordinate_identical = False
            for c_PDB in l_coordinates_PDB:
                for i in range(3):
                    if abs(c_PISA[i]-c_PDB[i]) >= 0.001:
                        break
                    if i == 2:
                        coordinate_identical = True
                if coordinate_identical == True:
                    break
            if coordinate_identical == True:
                continue
            else:
                ATOM_identical = False
                break

        return ATOM_identical


    def parse_coordinates(self,lines_in,):

        d_coordinates = {'ATOM':{},'HETATM':{},}
        d_lines = {'ATOM':[],'HETATM':[],}
##        set_chains_HETATM = set()

        for i in range(len(lines_in)):
            line = lines_in[i]
            record = line[:6].strip()
            if record in ['ATOM','HETATM',]:
                line_out = line[:6]+line[11:30]+line[30:37]+line[38:45]+line[46:53]+line[60:]
                line_out = line[:6]+line[11:30]+line[30:54]+line[60:]
                atom_name = line[12:16]
                res_name = line[17:20].strip()
                chain = line[21]
                res_no = int(line[22:26])
                iCode = line[26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinate = Numeric.array([x, y, z])
                if res_name == 'DOD':
                    stop
                if res_name == 'HOH':
                    continue
                
                ## append coordinate to dictionary
                if not chain in d_coordinates[record].keys():
                    d_coordinates[record][chain] = {}
                if not res_no in d_coordinates[record][chain].keys():
                    d_coordinates[record][chain][res_no] = {}
                if not iCode in d_coordinates[record][chain][res_no].keys():
                    d_coordinates[record][chain][res_no][iCode] = {}
                if not atom_name in d_coordinates[record][chain][res_no][iCode].keys():
                    d_coordinates[record][chain][res_no][iCode][atom_name] = []
                d_coordinates[record][chain][res_no][iCode][atom_name] += [coordinate]

                ## append line to dictionary
                d_lines[record] += [line_out]

##                ## append chain ID to list
##                if record == 'HETATM':
##                    if chain == ' ':
##                        chain = '-'
##                        chain_PISA = '[%s]%s:%i' %(res_name,chain,res_no)
##                        set_chains_HETATM |= set([chain_PISA])
##                    else:
##                        set_chains_HETATM |= set([chain])

##        l_chains_HETATM = list(set_chains_HETATM)

        return d_coordinates, d_lines


    def rsync(self):

        path_pdb = '/oxygenase_local/data/biounit/'
        os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/biounit/coordinates/divided/ %s' %(path_pdb))

        return
        

    def gunzip(self):

        path_pdb = '/oxygenase_local/data/biounit/'
        
        subdirs = os.listdir(path_pdb)
        subdirs.sort()
        for subdir in subdirs:
            print subdir
            files = os.listdir(path_pdb+subdir)
            for file in files:
                if file[-2:] == 'gz':
                    ## gunzip
##                    if os.path.isfile('%s%s/%s' %(path_pdb,subdir,file[:-3])):
##                        os.remove('%s%s/%s' %(path_pdb,subdir,file[:-3]))
                    os.system('gunzip -c %s%s/%s > %s%s/%s' %(path_pdb,subdir,file,path_pdb,subdir,file[:-3]))

        return


    def parse_REMARK350_biomolecules(self,d_transform, lines, i):

        biomolecules = lines[i][23:80].replace(' ','').split(',')
        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
                break

            elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                chains = set() ## e.g. 1upp.pdb (vs 8ruc.pdb)
                line_chains = lines[j][41:80]
                chains |= self.parse_REMARK350_chains(line_chains)

            elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
                line_chains = lines[j][53:80]
                chains |= self.parse_REMARK350_chains(line_chains)

            elif ',' in lines[j][11:80]:
                if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
                    line_chains = lines[j][11:80]
                    chains |= self.parse_REMARK350_chains(line_chains)

            ## count and parse chain transformations
            ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj.pdb)
            elif lines[j][13:19] in ['BIOMT3','SMTRY3']:

                matrixno = int(lines[j][19:24])
                ## parse transformation matrix
                matrixrow1 = lines[j-2][24:].split()
                matrixrow2 = lines[j-1][24:].split()
                matrixrow3 = lines[j-0][24:].split()
                matrixrows = [matrixrow1,matrixrow2,matrixrow3,]

                ## append transformation matrix to dictionary
                for biomolecule in biomolecules:

                    biomolecule = int(biomolecule)

                    ## biomolecule
                    if biomolecule not in d_transform.keys():
                        d_transform[biomolecule] = {}

                    ## biomolecule > matrices
                    if 'matrices' not in d_transform[biomolecule].keys():
                        d_transform[biomolecule]['matrices'] = {}
                    ## matrices > matrixno > matrix
                    d_transform[biomolecule]['matrices'][matrixno] = matrixrows

                    ## biomolecule > chains
                    if 'chains' not in d_transform[biomolecule].keys():
                        d_transform[biomolecule]['chains'] = {}
                    for chain in chains:
                        ## chains > chain
                        if chain not in d_transform[biomolecule]['chains'].keys():
                            d_transform[biomolecule]['chains'][chain] = set()
                        d_transform[biomolecule]['chains'][chain] |= set([matrixno])

        return d_transform
        

    def parse_REMARK350_chains(self, line_chains):

        ## if sentence necessary due to e.g. 1qgc.pdb
        ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: 1 2 3 4 5
        if ',' not in line_chains:
            chains = line_chains.split()
        else:
            ## remove 'AND' from the line of chains (e.g. problem with 1rhi.pdb)
            ## replace '.' in the line of chains (e.g. problem with 1rbo.pdb and 1qgc.pdb)
            chains = line_chains.replace('AND',',').replace('.',',').replace(' ',',').split(',')

        ## loop removal of blank chains necessary due to e.g. 2g8g.pdb
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
    comparison().main()

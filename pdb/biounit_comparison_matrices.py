#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

## replace [1] with [biounit] in the case of PDB transformations

import sys, Numeric, biounit

class comparison:

    def main(self):

        '''Script for comparison of biological units from the PDB and from MSD-PISA'''

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
            if pdb == '2bsu':
                pdb = '2v2w'

##            if pdb in [
##                '1fiv','1hef','1heg','1e5a', ## ligand overlap upon transformation, two alternative ligand binding conformations?!' e.g. 1bm7!!!
##                '1osv','1oxn','1oxq','2pcp','1u9l','1igj','1jq9','1jq8','1lpk','1lpg','1lpz', ## biounits should be identical, but different number of ligands in each PDB transformed biounit (water atoms often transformed incorrectly!!)
##                '1c1r', ## grey region
##
##                ##
##                ## multiple biounits
##                ##
####                '1fpu',
##                ## identical multimer in PISA and PDB
##                '1a4k','1a94','1pkx','1wc1','1fj4','1h1p','1cgl','1dqx','1eix','1fkn','1fl3','1h1s','1hsh','1i7z','2cht','1uw6','1jqy','1los','1lrh','1m4h','1a08','1b3l','1is0','1mh5','1mjj','1njj','1p1q','1q4k','3tmk','1umw','1uv6','1uz8',
##                ## tommy error
####                '2dqt', ## dimer in PISA and PDB
##
##                ##
##                ## other PISA multimers and/or interfaces are stable in solution
##                ##
##                '1qhc','1fq5','1it6','1jn4','2jxr','1yei','1tuf','1e2k','1e2n','1oe7','1e2l','1nms','1afk','1fch','1fkf','1kc7','6std','7std','1fzj','1fzk','1slg','1sle','1jyq','1kyv','1o9d','1oe8','1os0','1p19','1qca','5std','1tyr','1ugp','1v48','1vfn','1vpo','1vwl','1aqc','1ghy','2izl',
##
##                ##
##                ## same size multimers, but different interfaces
##                ##
##                '1loq','1lyx','1adl','2ans','1kll','1trd','1w72','1b55','1oko','1lyb','3pck','3pcj',
##
##                ##
##                ## different multimers
##                ##
##
##                ## dimer in PISA, monomer in PDB
##                '1oar','1s39','1udt','1p6e','1caq','1rd4','1p9p','1q63','1qi0','1gpk','1lf9','1k1y','1k4g','1oim','1b8y','1gpn','1ow4','1h6h','1lee','1qy2','1q4w','1w3j','1c5s','1ciz','1d7j','1dy4','1f3e','1ghz','1imx','1j17','1jt1','1k4h','1kpm','1l2s','1lnm','1m13','1m48','1n2v','1njs','1nw7','1nw5','1oif','1p28','1q65','1q66','1q91','1qft','1qy1','1r5y','1s38','1sqn','1sw1','1uho','1uj6','1uj5','1uz1','1wm1','1xzx','5yas','1d7i',
##                ## multimer in PISA, monomer in PDB
##                '1b8o','1b8n','1j4r','1g7v','1lf2','830c','1jn2','1fv0','4tmk','5tmp','2usn', '1usn',
##
##                ## correct multimer might be in "grey area" (deltaG_dissociation ~< 0)
##                ## PQS interfaces might be different from PDB interaces
##                '1b42','1bky','1bra','3mag','3mct', ## monomer in PISA, dimer in PDB
##                '2csn','1iup','1gz9','1igb','1ii5','1p1n', ## monomer in PISA, dimer in PDB, dimer in PQS
##                '1kdk','1lhw', ## monomer in PISA, dimer in PDB, hexamer in PQS
##                '1f8d', ## monomer in PISA, tetramer in PDB, tetramer in PQS, tetramer in PISA upon removal of ligands
##                '1bm7','1n51','2tmn','4tmn','5tmn','1wht', ## dimer in PISA, tetramer in PDB, tetramer in PQS
##                '1m5w', ## dimer in PISA, octamer in PDB, octamer in PQS
##                '1ftm', ## trimer in PISA, monomer in PDB, dimer in PQS
##                '1awi', ## monomer/dimer in PISA, trimer in PDB, trimer in PQS
##                '1a99', ## tetramer in PISA, dimer in PDB
##                
##                
##                ##
##                ## ligand positions
##                ##
##                ## identical multimers in PISA and PDB, but different ligand/ion/sugar position(s)
##                '2cgr','3gss','1elr','1elb','1hyo','1gvx','1gyx','1jet','2gss','1gyy','1gvu','1qkb','1ofz','1b9j','1jao','1jeu','1jev','1af6','10gs','11gs','1kui','1kuk','1kug','1obx','1ogx','1px4','1qka','1ur9','1e6s','1e6q','1e70',
##
##                '1h22','1h23', ## acetylcholine esterase
##                ## acetylcholine esterase
##                ## dimer
##                ## 4 helix bundle interface
##                ## 0.490nm between LYS530NZ and ASP365ODD2
##                ## 0.263nm between LYS530NZ and ASP369ODD2
##                ## hydrophobic core involving LEU366,LEU373,PHE527,LEU531
##                ## other AChE structures (1j06,1j07,1n5r,1n5m) have similar dimer interfaces with ligand in between interfaces
##                ## PDB molecule of the month states it is a dimer, but no details about the dimer interface
##                ]:
##                continue


            ##
            ##
            ##
            d_transformations_PISA, d_chains_PISA = biounit.biounit().parse_pisa_multimers(pdb)


            ##
            ##
            ##
            d_transformations_REMARK350 = {}
            set_water = set()
            set_nonwater = set()
            fd = open('/oxygenase_local/data/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb),'r')
            lines = fd.readlines()
            fd.close()
            for i in range(len(lines)):
                line = lines[i]
                record = line[:6].strip()
                if record == 'REMARK':
                    remark = int(line[7:10])
                    if remark == 350:
                        if line[11:23] == 'BIOMOLECULE:':
                            d_transformations_REMARK350 = self.parse_REMARK350_biomolecules(d_transformations_REMARK350, lines, i)
                elif record in ['ATOM','HETATM',]:
                    res_name = line[17:20]
                    chain = line[21]
                    if res_name == 'HOH':
                        set_water |= set([chain])
                    else:
                        set_nonwater |= set([chain])
            set_water -= set_nonwater

            ## remove water transformations
            if d_transformations_REMARK350 != {}:
                for chain in d_transformations_REMARK350[1]['chains'].keys():
                    if chain in set_water:
                        del d_transformations_REMARK350[1]['chains'][chain]

            fd = open('m3.txt')
            s = fd.read()
            fd.close()
            l_pdbs = s.split(',')
            print len(l_pdbs)
            l_pdbs = list(set(l_pdbs))
            print len(l_pdbs)
            stop
            if d_transformations_REMARK350.keys() not in [[],[1],]:
                print d_transformations_REMARK350.keys()
                set_matrices = set()
                for chain in d_transformations_REMARK350[2]['chains'].keys():
                    set_matrices |= d_transformations_REMARK350[2]['chains'][chain]
                m1 = False
                m2 = False
                for matrix_no in set_matrices:
                    matrix = d_transformations_REMARK350[2]['matrices'][matrix_no]
                    if matrix == [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]:
                        m1 = True
                    else:
                        m2 = True
                if m1 == True and m2 == True:
                    s = 'm3.txt'
                elif m1 == True and m2 == False:
                    s = 'm1.txt'
                elif m1 == False and m2 == True:
                    stop
                else:
                    stop
                fd = open(s,'a')
                fd.write('%s,' %(pdb))
                fd.close()
                continue
            else:
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

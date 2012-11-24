#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

## replace [1] with [biounit] in the case of PDB transformations

## is there a risk PISA coordinates are formed from transformation instead of existing coordinates?! 1a99???

import sys, Numeric, biounit, os, math

class comparison:


    def loop(self,path_pdb_bio,path_pdb_asu):

##        self.rsync(path_pdb_bio)
##        self.gunzip(path_pdb_bio)

        ## txt to list
        fd = open('pdbbind_2007.txt','r')
        s = fd.read()
        fd.close()
        l_pdbs = s.split()

        ## list to dic
        d_pdbs = {}
        for pdb in l_pdbs:
            if not pdb[1:3] in d_pdbs.keys():
                d_pdbs[pdb[1:3]] = [pdb]
            else:
                d_pdbs[pdb[1:3]] += [pdb]

        ## dic to sorted list
        l_subs = d_pdbs.keys()
        l_subs.sort()
        l_pdbs = []
        for sub in l_subs:
            for pdb in d_pdbs[sub]:
                l_pdbs += [pdb]
            
        for pdb in l_pdbs:


##            ## v2 atom names in PDB asu (not a problem when using matrices)
##            if pdb not in [
##                '1a8i','1c1u','1c1v','1c5n','1c5o','1epv','1f4e','1f4g','1f4f',
##                '1ghy','1ghv','1ghw','1gi7','1gi9','1gja','1gj7','1gjb','1gj4',
##                '1gjd','1gjc','1gj9','1gj8','1gj5','2izl','1kzk','1niu','1o3p',
##                '1p57','1toj','1toi','1tog','1tok','1tuf','1b55','1w0y','1xge',
##                ]:
##                continue

##            if l_pdbs.index(pdb)+1 < int(sys.argv[1]):
##                continue
            print pdb, l_pdbs.index(pdb)+1, len(l_pdbs)

            ## superseded structure(s)
            if pdb == '2bsu':
                pdb = '2v2w'

            if pdb in [

                ## v2 atom names in PDB asu (not a problem when using matrices)
                '1a8i','1c1u','1c1v','1c5n','1c5o','1epv','1f4e','1f4g','1f4f',
                '1ghy','1ghv','1ghw','1gi7','1gi9','1gja','1gj7','1gjb','1gj4',
                '1gjd','1gjc','1gj9','1gj8','1gj5','2izl','1kzk','1niu','1o3p',
                '1p57','1toj','1toi','1tog','1tok','1tuf','1b55','1w0y','1xge',

                ## ligand overlap upon transformation. two alternative ligand binding conformations?!' e.g. 1bm7!!!
##                '1e5a',##'1fiv','1heg',

                ## altloc
                    ## altloc used for alternative temperature factors (but identical coordinates!)
                    '1a69','1qk3',
                    ## altloc used for identical coordinate and temperature factors
                    '1cru',

                ##
                ## same size multimers, but different interfaces
                ##
                '1adl','1af2','2ans','1ct8','1kll','1loq','1lyx','1lyb','1oko',
                '3pck','3pcj',

                ##
                ## multimer difference
                ##

                ## monomer in PDB
                    ## dimer in PISA
                    '2ayr','1b11','1b8y','1c1r','1c5s','1caq','1ciz','1d7i',
                    '1d7j','1dy4','1f3e','1ghz','1gpk','1gpn','1h22','1h23',
                    '1h6h','1imx','1j17','1jt1','1k1y','1k4g','1k4h','1kpm',
                    '1l2s','1lee','1lf9','1lnm','1m13','1m48','1n2v','1njs',
                    '1nw7','1nw5','1oar','1oim','1oif','1ow4','1p28','1p6e',
                    '1p9p','1q4w','1q63','1q65','1q66','1q91','1qft','1qhc',
                    '1qi0','1qy2','1qy1','1r5y','1rd4','1s39','1s38','1sqn',
                    '1sw1','1udt','1uho','1uj6','1uj5','1uz1','1w3j','1wm1',
                    '1xzx','5yas',
                    ## multimer in PISA
                    '1b8n','1b8o','830c','1ftm','1fv0','1g7v','1j4r','1jn2',
                    '1lf2','1oxn','1oxq','5tmp','4tmk','2usn','1usn',

                ## monomer in PISA
                    ## dimer in PDB
                    '1hsl','1iup',

                ## multimer in PISA, multimer in PDB
                    ## incorrect PISA peptide ligand transformation
                    '1nh0','1o2g','1o3p',
                    ## dimer in PISA, tetramer in PDB
                    '1bm7','1n51','2tmn','4tmn','5tmn','1wht','1a1b','1a1c','1a1e','1bm7','1bq4','1ecq','1ec9','1rdt','1swg','1fm9','1fo0','1atl',
                    '1nw4', ## dimer in PISA, hexamer in PDB
                    '1m5w', ## dimer in PISA, octamer in PDB, octamer in PQS
                    '1a99','2bmz', ## tetramer in PISA, dimer in PDB

               
                ##
                ## ligand difference
                ##
                
                ## nontransformation
                    '1trd','1uz8',
                    ## some of the ligands not transformed in PDB (biou larger than asu)
                    '2a5b','2a5c','2a8g','1d09','1jcx',
                    ## one or more ligands of asu transformed to a single biounit in PDB (asu larger than biou)
                    '1a4k','2cht','1fl3','1hsh','1hsg','1igj','1lrh','2pcp',
                    ## biounits should be identical, but different number of ligands in each PDB transformed biounit (water atoms often transformed incorrectly!!)
                    '1jq9','1jq8','1lpk','1lpg','1lpz',
                    ## PISA fails to transform peptide ligand
                    '1bhx',
                ## positions
                    '2aj8','2aoc','2aod','2aoe','2aog','1apv','1apw','2bpy',
                    '1c5c','1c5y','2cgr','1dl7','1gi8','1gyx','1gyy','1jqy',
                    '1lgt','1mfi','1ogx','1tyr','2bjm','1ofz','1e6s','2bt9',
                    ## MG positions
                    '1af6',
                    ## SO4 positions
                    '2avo','2avs','1b6j','1b6k','1b6l','1b6m','1b6p','1d4l','1d4k','1elb','4er2','3gst','1gvu','1gvx','1hii','1hn4','1mtr','1sdu',
                    ## NA positions
                    '1sb1',
                    ## CL positions
                    '2avv','1sdv','1sdt',
                    ## CD positions
                    '1d6v','1kuk','1kui','1kug',
                    ## NI positions
                    '1elr','1hyo',
                    ## CA and/or ZN positions
                    '1ghy','1jao','1jaq','1p1q','7std','6std','1yei','1yej',
                    '1b55',
                    ## U1 positions
                    '1b05','1b0h','1b2h','1b4h','1b5h','1b5i','1b5j','1b6h',
                    ## IUM positions
                    '1b1h','1b32','1b3f','1b3g','1b3h','1b3l','1b40','1b46',
                    '1b4z','1b51','1b52','1b58','1b9j','1jet','1jeu','1jev',
                    '2olb','1qkb','1qka','2rkm',
                    ## IOD positions
                    '2bmk',
                    ## MES positions
                    '10gs','11gs','3gss','2gss','1nu3',
                    ## IPA positions
                    '1dqn',
                    ## GOL positions
                    '2avq','1e6q','1e70','1f74','2ien','1u0g','1ur9',
                    ## MPD positions
                    '1fzj','1ro7',
                    ## MRD positions
                    '1lan','1lcp',
                    ## CIT positions
                    '1gcz',
                    ## ACT positions
                    '1k6p','1k6c','1k6v','1k6t',
                    ## CO positions
                    '1obx',
                    ## DMS positions
                    '1os0','1px4','1qf0','1qf2','1qf1',
                    ## DMF positions
                    '1ppk',
                ]:
                continue

            else:
                s = self.main(pdb,path_pdb_bio,path_pdb_asu)
                print pdb, s

        return


    def alison(self,pdb,):

        path_pdb_bio = '/data/pdb_biounits/'
        path_pdb_asu = '/data/remediated_pdb/'
        try:
            s = self.main(pdb,path_pdb_bio,path_pdb_asu,verbose=False,)
        except:
            s = 'Syntax error in PDB. PISA and PDB multimer state cannot be compared.'

        return s


    def main(self,pdb,path_pdb_bio,path_pdb_asu,verbose=True):

        '''Script for comparison of biological units from the PDB and from MSD-PISA'''

        ##
        ##
        ##
        d_transformations_PISA,status = biounit.biounit().parse_pisa_multimers(pdb,verbose=verbose)
        if status == 'Broken composition in PA Graph':
            s = 'PISA unable to determine quarternary structure'
            return s
        if status != 'Ok':
            s = status
            return status
        if d_transformations_PISA == {}:
            s = 'PISA unable to identify any stable quarternary structures'
            return s
        l_assemblies = d_transformations_PISA.keys()

        d_coordinates_PDB = {}
        d_lines_PDB = {}
        l_chains_ATOM = []
        for bm in range(1,99):
            if bm == 1 or os.path.isfile('%s%s/%s.pdb%i' %(path_pdb_bio,pdb[1:3],pdb,bm)):
                fd = open('%s%s/%s.pdb%i' %(path_pdb_bio,pdb[1:3],pdb,bm),'r')
                lines = fd.readlines()
                fd.close()
                d_coordinates_PDB[bm], d_lines_PDB[bm] = self.parse_coordinates(lines)
                l_chains_ATOM += d_coordinates_PDB[bm]['ATOM'].keys()
                continue
            else:
                break

        l_chains_ATOM = list(set(l_chains_ATOM))
        l_biomolecules = d_coordinates_PDB.keys()

        for bm in l_biomolecules:
            biounits_identical = False
            chains_identical = False
            interfaces_identical = False
            l_chains_PDB = d_coordinates_PDB[bm]['ATOM'].keys()
            for assembly in l_assemblies:
                if verbose == True:
                    print bm, assembly
                l_chains_PISA = []
                for chain in d_transformations_PISA[assembly]['chains'].keys():
                    if len(chain) == 1 and chain != '-':
                        l_chains_PISA += [chain]
                ## exclude non-ATOM chains from PISA assembly chains
                l_chains_PISA = list(set(l_chains_ATOM)&set(l_chains_PISA))
##                print d_transformations_PISA[assembly]['chains'].keys()
##                print d_coordinates_PDB[bm]['ATOM'].keys()
##                print d_coordinates_PDB[bm]['HETATM'].keys()

                ## different chains (different IDs)
                if len(set(l_chains_PISA) ^ set(l_chains_PDB)) > 0:
                    if verbose == True:
                        print 'different chains', l_chains_PDB, l_chains_PISA
                    continue
                else:
                    d_lines_PISA, d_coordinates_PISA = biounit.biounit().parse_pdb_coordinates(pdb, d_transformations_PISA, assembly,path_pdb_asu)
                    len_PISA = float(len(d_lines_PISA[assembly]['ATOM']))
                    len_PDB = float(len(d_lines_PDB[bm]['ATOM']))
                    ## different chains (different number of transformations)
                    if not (
                        len_PISA % len_PDB == 0 and
                        len_PDB % len_PISA == 0
                        ):
                        ## len of coordinates not a multiplum of each other
                        if (
                            len_PISA % len_PDB != 0 and
                            len_PDB % len_PISA != 0
                            ):
                            ## print duplicate lines, which are most likely caused by transformation of coordinates with v2 atom names
                            for line in d_lines_PDB[bm]['ATOM']:
                                count = d_lines_PDB[bm]['ATOM'].count(line)
                                if count > 1:
                                    if verbose == True:
                                        print line
                            ## check that the number of coordinates are a multiplum of each other
                            if (
                                round(min(len_PDB,len_PISA) % math.modf(len_PISA/len_PDB)[0]*len_PDB,8) != 0. or
                                round(min(len_PDB,len_PISA) % math.modf(len_PDB/len_PISA)[0]*len_PISA,8) != 0.
                                ):
                                a = list(set(d_lines_PISA[assembly]['ATOM'])-set(d_lines_PDB[bm]['ATOM']))
                                b = list(set(d_lines_PDB[bm]['ATOM'])-set(d_lines_PISA[assembly]['ATOM']))
                                a.sort()
                                b.sort()
                                b.reverse()
                                if verbose == True:
                                    print a
                                    print b
                                    print len_PDB
                                    print len_PISA
                                    print min(len_PDB,len_PISA) % math.modf(len_PISA/len_PDB)[0]*len_PDB
                                    print min(len_PDB,len_PISA) % math.modf(len_PDB/len_PISA)[0]*len_PISA
                                stop_not_expected_or_v2_atom_names
                        if verbose == True:
                            print 'different chains', len(d_lines_PDB[bm]['ATOM']), len(d_lines_PISA[assembly]['ATOM'])
                        continue
                    ## identical chains (identical IDs, identical number of transformations)
                    else:
                        chains_identical = True

                ATOM_identical = self.identical_d_coordinates('ATOM',d_coordinates_PISA,d_coordinates_PDB,assembly,bm,verbose=verbose,)
                if 'HETATM' in d_coordinates_PISA[assembly].keys():
                    HETATM_identical = self.identical_d_coordinates('HETATM',d_coordinates_PISA,d_coordinates_PDB,assembly,bm,)
                else:
                    HETATM_identical = True
                if verbose == True:
                    print ATOM_identical, HETATM_identical
                
                ATOM_identical2 = self.identical_d_coordinates('ATOM',d_coordinates_PDB,d_coordinates_PISA,bm,assembly,)
                if ATOM_identical != ATOM_identical2:
                    if verbose == True:
                        print ATOM_identical, ATOM_identical2
                    stop

                if chains_identical == True and ATOM_identical == True and HETATM_identical == True:
                    biounits_identical = True
                    if verbose == True:
                        print bm,assembly,'identical'
                    break
                elif ATOM_identical == False:
                    if verbose == True:
                        print bm,assembly,'different interfaces'
                    continue
                elif ATOM_identical == True and HETATM_identical == False:
                    interfaces_identical = True
                    ## differences between sets of lines are not representative of differences if coordinates differ by less than 0.0001nm
                    a = list(set(d_lines_PISA[assembly]['HETATM'])-set(d_lines_PDB[bm]['HETATM']))
                    b = list(set(d_lines_PDB[bm]['HETATM'])-set(d_lines_PISA[assembly]['HETATM']))
                    a.sort()
                    b.sort()
                    if verbose == True:
                        print a
                        print b
                        print bm,assembly,'different ligands'
                else:
                    if verbose == True:
                        print chains_identical, ATOM_identical, HETATM_identical
                    stop_notexpected

            if biounits_identical == True:
                s = 'biounits identical'
            elif biounits_identical == False:
                if interfaces_identical == True:
                    s = 'different ligand transformations'
                elif chains_identical == True and interfaces_identical == False:
                    s = 'different peptide interfaces'
                elif chains_identical == False and interfaces_identical == False:
                    s = 'different multimers'
                else:
                    if verbose == True:
                        print chains_identical, interfaces_identical,
                    stop_not_expected

        return s


    def identical_d_coordinates(self,record,d_coordinates_PISA,d_coordinates_PDB,assembly,bm,verbose=True):

        ATOM_identical = True

        for chain in d_coordinates_PISA[assembly][record].keys():
            if not chain in d_coordinates_PDB[bm][record].keys():
                if record == 'ATOM':
                    if verbose == True:
                        print chain
                    stop_transformation_PDB
                if record == 'HETATM':
                    if verbose == True:
                        print chain
                    ATOM_identical = False
                    continue
            for res_no in d_coordinates_PISA[assembly][record][chain].keys():
                if not res_no in d_coordinates_PDB[bm][record][chain].keys():
                    if record == 'ATOM':
                        if verbose == True:
                            print chain
                        stop_transformation_PDB
                    if record == 'HETATM':
                        if verbose == True:
                            print res_no
                        ATOM_identical = False
                        continue
                for iCode in d_coordinates_PISA[assembly][record][chain][res_no].keys():
                    for altloc in d_coordinates_PISA[assembly][record][chain][res_no][iCode].keys():
                        res_name = d_coordinates_PISA[assembly][record][chain][res_no][iCode][altloc]['res_name']
                        atom_identical = True
                        for atom_name in d_coordinates_PISA[assembly][record][chain][res_no][iCode][altloc]['atom_names'].keys():
                            l_coordinates_PISA = d_coordinates_PISA[assembly][record][chain][res_no][iCode][altloc]['atom_names'][atom_name]
                            try:
                                l_coordinates_PDB = d_coordinates_PDB[bm][record][chain][res_no][iCode][altloc]['atom_names'][atom_name]
                            except:
                                if verbose == True:
                                    print chain, res_no, iCode, altloc, res_name
                                    print atom_name, d_coordinates_PDB[bm][record][chain][res_no][iCode][altloc]['atom_names'].keys()
                                continue
                                stop_v2_atom_name
                            identical_PISA = self.identical_l_coordinates(l_coordinates_PISA,l_coordinates_PDB)
                            identical_PDB = self.identical_l_coordinates(l_coordinates_PDB,l_coordinates_PISA)
                            if identical_PISA == False or identical_PDB == False:
                                ATOM_identical = False
                                atom_identical = False
                                break
                        if ATOM_identical == False:
                            ## print first occurence and break if ATOM record
                            if record == 'ATOM':
                                if verbose == True:
                                    print record, chain, res_no, iCode, altloc, res_name, atom_name
                                break
                            ## print all occurences if HETATM record
                            if record == 'HETATM' and atom_identical == False:
                                if verbose == True:
                                    print record, chain, res_no, iCode, altloc, res_name
                    if ATOM_identical == False and record == 'ATOM':
                        break
                if ATOM_identical == False and record == 'ATOM':
                    break
            if ATOM_identical == False and record == 'ATOM':
                break

        return ATOM_identical


    def identical_l_coordinates(self, l_coordinates_PISA,l_coordinates_PDB):

        ATOM_identical = True
        for c_PISA in l_coordinates_PISA:
            coordinate_identical = False
            for c_PDB in l_coordinates_PDB:
                for i in range(3):
                    if abs(c_PISA[i]-c_PDB[i]) > 0.001:
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
                altloc = line[16]
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
                    chain1=lines_in[i-1][21]
                    chain2=lines_in[i+1][21]
                    continue
                
                ## append coordinate to dictionary
                if not chain in d_coordinates[record].keys():
                    d_coordinates[record][chain] = {}
                if not res_no in d_coordinates[record][chain].keys():
                    d_coordinates[record][chain][res_no] = {}
                if not iCode in d_coordinates[record][chain][res_no].keys():
                    d_coordinates[record][chain][res_no][iCode] = {}
                if not altloc in d_coordinates[record][chain][res_no][iCode].keys():
                    d_coordinates[record][chain][res_no][iCode][altloc] = {'res_name':res_name,'atom_names':{}}
                if not atom_name.strip() in d_coordinates[record][chain][res_no][iCode][altloc]['atom_names'].keys():
                    d_coordinates[record][chain][res_no][iCode][altloc]['atom_names'][atom_name.strip()] = []
                d_coordinates[record][chain][res_no][iCode][altloc]['atom_names'][atom_name.strip()] += [coordinate]

                ## append line to dictionary
                d_lines[record] += [line_out]
            if record == 'SPRSDE':
                print line
                stop

        return d_coordinates, d_lines


    def rsync(self,path_pdb_bio):

        os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/biounit/coordinates/divided/ %s' %(path_pdb_bio))

        return
        

    def gunzip(self,path_pdb_bio):

        subdirs = os.listdir(path_pdb_bio)
        subdirs.sort()
        for subdir in subdirs:
            print subdir
            files = os.listdir(path_pdb_bio+subdir)
            for file in files:
                if file[-2:] == 'gz':
                    ## gunzip
##                    if os.path.isfile('%s%s/%s' %(path_pdb_bio,subdir,file[:-3])):
##                        os.remove('%s%s/%s' %(path_pdb_bio,subdir,file[:-3]))
                    os.system('gunzip -c %s%s/%s > %s%s/%s' %(path_pdb_bio,subdir,file,path_pdb_bio,subdir,file[:-3]))

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

    path_pdb_bio = '/oxygenase_local/data/biounit/'
    path_pdb_bio = '/data/pdb_biounits/'
    path_pdb_asu = '/data/remediated_pdb/'

    comparison().loop(path_pdb_bio,path_pdb_asu)

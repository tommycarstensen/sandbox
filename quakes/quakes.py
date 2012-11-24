#!/software/bin/python
#
#$Id: quakes.py 266 2007-11-01 13:15:43Z tc $
#
#Tommy Carstensen, University College Dublin, 2007

## built-ins
import os, sys, numpy, math, copy, time, re, urllib2, shutil
## non-built-ins
import quakes_pdb_parser, quakes_gif, quakes_dihedral, quakes_readwritehtml
sys.path.append('/home/tc/svn/tc_sandbox/misc/')
import gnuplot
sys.path.append('/home/tc/svn/tc_sandbox/pdb/')
import biounit,parse_pdb
sys.path.append('/home/tc/svn/Protool/')
import geometry
sys.path.append('/home/tc/svn/PEAT_DB/')
import sequence_alignment

class quakes:

    def main(self):

##        self.check_dirs()
##        self.remove_obsolete()

        ##
        self.step_one()

        ##
        ## description of the dataset
        ##
        quakes_readwritehtml.dataset_description()

        ##
        ## analysis of phi/psi changes due to mutations
        ##
        self.plot_phipsi_mutants()

        ##
        ## analysis of rmsd
        ##
        quakes_readwritehtml.analyze_rmsd()

        return


    def step_one(self):

        fd = open('errorpdbs.txt','r')
        s = fd.read()
        fd.close()
        errorpdbs = eval(s)

        if '-analyze_only' in sys.argv:
            pass

############# manual selection
        elif '-manual' in sys.argv:

            if len(sys.argv) > 2:
                for i in range(len(sys.argv)-1,-1,-1):
                    if sys.argv[i] == '-manual':
                        self.l_pdbs = sys.argv[i+1:]
                        break
            else:
                self.l_pdbs = list(set([
##                    ## SEQRES, ATOM, REMARK465 cases
##                    '3ee0','2gp9',
##                    '1abj','2thf',
                    ]))

##            fd = open('transform_error.txt','r')
##            lines = fd.readlines()
##            fd.close()
##            for line in lines:
##                self.l_pdbs += line.split()[:2]
##            self.l_pdbs = list(set(self.l_pdbs))

            ## removal
            for pdb in errorpdbs:
                if pdb in self.l_pdbs:
                    self.l_pdbs.remove(pdb)
                    print 'error', pdb

            if '-sphere' in sys.argv:
                bool_sphere = True
            else:
                bool_sphere = False

            bool_do_single_mutant = True
            do_gif_and_pdb = True
            if '-justdohtm' in sys.argv:
                bool_do_single_mutant = False
                do_gif_and_pdb = False

            if len(self.l_pdbs) == 1:
                pdb_main = ''.join(self.l_pdbs)
                fd = open('htm/%s.htm' %(''.join(self.l_pdbs)),'r')
                lines = fd.readlines()
                fd.close()
                set_pdbs = set()
                for line in lines[33:-1:28]:
                    index2 = line.index('</a>')
                    index1 = line[:index2].rindex('>')+1
                    pdb = line[index1:index2]
                    set_pdbs |= set([pdb])
                for line in lines[34:-1:28]:
                    index2 = line.index('</a>')
                    index1 = line[:index2].rindex('>')+1
                    pdb = line[index1:index2]
                    set_pdbs |= set([pdb])
                l_l_pdbs = []
                for pdb in set_pdbs:
                    l_l_pdbs += [[pdb_main,pdb,]]
            else:
                l_l_pdbs = [self.l_pdbs]

            for self.l_pdbs in l_l_pdbs:
                d_rmsd,d_header = self.analyze_pdbs(
                    verbose=True,
                    bool_sphere = bool_sphere,
                    bool_do_single_mutant = bool_do_single_mutant,
                    do_gif_and_pdb = do_gif_and_pdb,
                    )
        ##        d_rmsd = read_rmsd_from_file()
                quakes_readwritehtml.write_html(d_rmsd, d_header, prefix='quickrmsd')

            print 'manual finished'
            stop_manual_finished
            return

        ##
        ## singlemutants
        ##
        elif '-singlemutants' in sys.argv:

            ## 1) find all single mutants in single_mutants txt files
            ## 2) read html files of set of single mutants and find all non-mutants
            ## 3) calculate sphere rmsd for non-mutants and single mutants

            ## remove redundant lines
            lines_main = []
            for s in '0123456789abcddefghijklmnopqrstuvwxyz':
                fd = open('single_point_mutations/%s.txt' %(s),'r')
                lines = fd.readlines()
                fd.close()
                pdbs_main = []
                lines_nonredundant = []
                for i in range(len(lines)):
                    pdbs = lines[i].split()[0:2]

                    pdbs.sort()
                    pdbs = ''.join(pdbs)
                    if pdbs not in pdbs_main:
                        pdbs_main += [pdbs]
                        lines_nonredundant += [lines[i]]
                        
                if len(lines) != len(lines_nonredundant):
                    fd = open('single_point_mutations/%s.txt' %(s),'w')
                    fd.writelines(lines_nonredundant)
                    fd.close()
                lines_main += lines_nonredundant

            lines = lines_main
           
            for i in range(len(lines)):
                if i < int(sys.argv[-1]):
                    continue
                print '\n\n********\n', i, len(lines)
                line = lines[i]
                self.l_pdbs = line.split()[:2]
                if len(set(self.l_pdbs)-set(errorpdbs)) < 2:
                    continue
                d_rmsd,d_header = self.analyze_pdbs(
                    do_gif_and_pdb = False, verbose=False,
                    bool_sphere = True,
                    )
                quakes_readwritehtml.write_html(d_rmsd, d_header, prefix='rmsd_mutants')

            print 'xlabel "distance from mutation (from CA of mutated residue to CA sphere center)"'
            print 'plot [-1:]"heavy_mutant_comparison.txt" u 9:10 t "r=[0:5[", "heavy_mutant.txt" u 9:10 t "r=[0:5[", "heavy_mutant_comparison.txt" u 9:11 t "r=[5:10[", "heavy_mutant.txt" u 9:11 t "r=[5:10[", "heavy_mutant_comparison.txt" u 9:12 t "r=[10:20[", "heavy_mutant.txt" u 9:12 t "r=[10:20[", "heavy_mutant_comparison.txt" u 9:13 t "r=[20:50[", "heavy_mutant.txt" u 9:13 t "r=[20:50[", "heavy_mutant_comparison.txt" u 9:14 t "r=[50:]", "heavy_mutant.txt" u 9:14 t "r=[50:]"'
            print 'plot [-1:]"heavy_mutant_sameSG_sameauth_comparison.txt" u 9:10 t "r=[0:5[", "heavy_mutant_sameSG_sameauth.txt" u 9:10 t "r=[0:5[", "heavy_mutant_sameSG_sameauth_comparison.txt" u 9:11 t "r=[5:10[", "heavy_mutant_sameSG_sameauth.txt" u 9:11 t "r=[5:10[", "heavy_mutant_sameSG_sameauth_comparison.txt" u 9:12 t "r=[10:20[", "heavy_mutant_sameSG_sameauth.txt" u 9:12 t "r=[10:20[", "heavy_mutant_sameSG_sameauth_comparison.txt" u 9:13 t "r=[20:50[", "heavy_mutant_sameSG_sameauth.txt" u 9:13 t "r=[20:50[", "heavy_mutant_sameSG_sameauth_comparison.txt" u 9:14 t "r=[50:]", "heavy_mutant_sameSG_sameauth.txt" u 9:14 t "r=[50:]"'
            print 'plot [-1:]"newradii_heavy_mutant_sameSG_sameauth_comparison.txt" u 9:10 t "r=[0:10[", "newradii_heavy_mutant_sameSG_sameauth.txt" u 9:10 t "r=[0:10[", "newradii_heavy_mutant_sameSG_sameauth_comparison.txt" u 9:11 t "r=[10:20[", "newradii_heavy_mutant_sameSG_sameauth.txt" u 9:11 t "r=[10:20[", "newradii_heavy_mutant_sameSG_sameauth_comparison.txt" u 9:12 t "r=[20:40[", "newradii_heavy_mutant_sameSG_sameauth.txt" u 9:12 t "r=[20:40[", "newradii_heavy_mutant_sameSG_sameauth_comparison.txt" u 9:13 t "r=[40:[", "newradii_heavy_mutant_sameSG_sameauth.txt" u 9:13 t "r=[40:["'

            print 'ylabel "CA RMSD within sphere"'
            print 'plot [-1:]"out_sameSG_sameauth_allresidues.txt" u 9:10 t "r=[0:10[", "out_sameSG_sameauth_allresidues.txt" u 9:11 t "r=[10:20[", "out_sameSG_sameauth_allresidues.txt" u 9:12 t "r=[20:40[", "out_sameSG_sameauth_allresidues.txt" u 9:13 t "r=[40:["'
            print 'ylabel "heavy RMSD within sphere"'
            print 'plot [-1:]"out_sameSG_sameauth_allresidues.txt" u 9:14 t "r=[0:10[", "out_sameSG_sameauth_allresidues.txt" u 9:15 t "r=[10:20[", "out_sameSG_sameauth_allresidues.txt" u 9:16 t "r=[20:40[", "out_sameSG_sameauth_allresidues.txt" u 9:17 t "r=[40:["'
            print 'ylabel "phipsi RMSD within sphere"'
            print 'plot [-1:]"out_sameSG_sameauth_allresidues.txt" u 9:18 t "r=[0:10[", "out_sameSG_sameauth_allresidues.txt" u 9:19 t "r=[10:20[", "out_sameSG_sameauth_allresidues.txt" u 9:20 t "r=[20:40[", "out_sameSG_sameauth_allresidues.txt" u 9:21 t "r=[40:["'
            print 'ylabel "chi1 RMSD within sphere"'
            print 'plot [-1:]"out_sameSG_sameauth_allresidues.txt" u 9:22 t "r=[0:10[", "out_sameSG_sameauth_allresidues.txt" u 9:23 t "r=[10:20[", "out_sameSG_sameauth_allresidues.txt" u 9:24 t "r=[20:40[", "out_sameSG_sameauth_allresidues.txt" u 9:25 t "r=[40:["'
            print 'fit f(x) "heavy_mutant_comparison.txt" u 9:10 via a,b'
                
            stop_single_mutants

        ##        
        ## clusters 95
        ##
        elif '-clusters' in sys.argv:

##            os.system('wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/NR/clusters95.txt')

            print 'reading clusters'
            fd = open('pdbS95bF.out','r')
            lines = fd.readlines()
            fd.close()

            print 'reading lines of clusters'
            for i in range(len(lines)):
                if i < int(sys.argv[-2]):
                    continue
                if i > int(sys.argv[-1]):
                    continue
                line = lines[i]
                l_pdbs = line.split()
                cluster = i
                self.l_pdbs = []
                for pdb in l_pdbs:
                    self.l_pdbs += [pdb[:4].lower()]
                self.l_pdbs = list(set(self.l_pdbs))
                print '-------cluster--------', cluster, len(self.l_pdbs)
                if '-skip' in sys.argv:
                    self.l_pdbs = self.l_pdbs[int(sys.argv[-4]):]

                self.l_pdbs.sort()

                ## skip if error pdb
                for pdb in errorpdbs:
                    try:
                        self.l_pdbs.remove(pdb)
                    except:
                        None

                ## skip if only one pdb in cluster
                if len(self.l_pdbs) <= 1:
                    continue

                try:
                    d_rmsd,d_header = self.analyze_pdbs(do_gif_and_pdb = True)
                except:
                    self.l_pdbs = [self.pdb1,self.pdb2,]
                    self.analyze_pdbs()
                    print 'cluster', cluster
                    print self.pdb1, self.pdb2
                    error_cluster

            print 'finished clusters', cluster

            return
##########
        return


    def check_dirs(self):

        if os.getcwd() != self.topdir:
            print os.getcwd()
            print self.topdir
            stop_wrong_dir

        for s_dir in ['htm','ps','txt','pdb','tmp','dssp','src','log',]:
            if not os.path.isdir(s_dir):
                os.mkdir(s_dir)
            if s_dir in ['txt','pdb',]:
                for s in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
                    if not os.path.isdir('%s/%s' %(s_dir,s.lower(),)):
                        os.mkdir('%s/%s' %(s_dir,s.lower(),))
            elif s_dir in ['dssp',]:
                for s1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
                    for s2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
                        if not os.path.isdir('%s/%s%s' %(s_dir,s1.lower(),s2.lower(),)):
                            os.mkdir('%s/%s%s' %(s_dir,s1.lower(),s2.lower(),))
                

        return


    def analyze_pdbs(self, verbose = False, do_gif_and_pdb = True, bool_sphere = False, bool_do_single_mutant = True):

        self.verbose = verbose

        ## do not exclude pdb2 from pdb1 loop if sequence similar to previous pdb1
        ## since pdb A == B, A != C, B == C
        ## but do not analyze sequence of pdb A,B and then pdb B,A since A == B and B == A are equivalent

        d_rmsd_identical = {}
        d_coordinates = {}
        d_hetero = {}
        d_ATOMseq = {}
        d_header = {}
        d_chains_intrapdb_sequence_identical = {}

        ## parse transform errors
        l_transform_errors = []

        fd = open('transform_error.txt','r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            l_transform_errors += [line.split()[:4]]

        ##
        ## loop 1 over pdbs
        ##
        t1 = time.clock()
        for i1 in range(len(self.l_pdbs)):

            self.pdb1 = pdb1 = self.l_pdbs[i1]

            if not pdb1 in d_header.keys():
                d_header[pdb1]  = quakes_pdb_parser.parse_header(
                    pdb1,self.path_pdb,
                    self.min_len_chain,
                    )
                d_chains_intrapdb_sequence_identical[pdb1] = self.identify_identical_chains_from_sequence_intra(
                    d_header,pdb1,
                    )

            skippdb = self.pdbskip(d_header, pdb1)
            if skippdb == True:
                continue

            ## identify biomolecule(s)
            d_biomolecules1 = self.identify_biomolecule(pdb1, d_header)
            if d_biomolecules1.keys() == []:
                print d_header[pdb1]['REMARK350']
                stop1_check_REMARK350
                
            ## reset dictionary of coordinates to save memory
            d_coordinates = {}

            ##
            ## loop 2 over pdbs
            ##
            for i2 in range(i1+1,len(self.l_pdbs)):

                self.pdb2 = pdb2 = self.l_pdbs[i2]

                ## do not compare to self
                if pdb1 == pdb2:
                    continue

                if not pdb2 in d_header.keys():
                    d_header[pdb2]  = quakes_pdb_parser.parse_header(
                        pdb2,self.path_pdb,
                        self.min_len_chain,
                        )
                    d_chains_intrapdb_sequence_identical[pdb2] = self.identify_identical_chains_from_sequence_intra(
                        d_header,pdb2,
                        )

                skippdb = self.pdbskip(d_header, pdb2)
                if skippdb == True:
                    continue

                ## identify biomolecule(s)
                d_biomolecules2 = self.identify_biomolecule(pdb2, d_header)
                if d_biomolecules2.keys() == []:
                    stop2

                ##
                ## print status
                ##
                t2 = time.clock()
                if t2-t1 > self.time_status_update or i2 == i1+1:
                    print 'analyzing sequence of %s (%5i/%5i) and %s (%5i/%5i)' %(pdb1, i1+1, len(self.l_pdbs), pdb2, i2+1, len(self.l_pdbs),)
                    t1 = t2

                ##
                ## loop over biomolecule(s) of pdb1
                ##
                for biomolecule1 in d_biomolecules1.keys():

                    bm1 = biomolecule1
                    bmchains1 = d_biomolecules1[biomolecule1]['chains']
                    long_peptide_chains1 = self.find_long_peptide_chains(pdb1,d_header,)
                    bmpolymercount1 = d_biomolecules1[biomolecule1]['polymercount']

                    ##
                    ## loop over biomolecule(s) of pdb2
                    ##
                    for biomolecule2 in d_biomolecules2.keys():

                        bm2 = biomolecule2
                        bmchains2 = d_biomolecules2[biomolecule2]['chains']
                        long_peptide_chains2 = self.find_long_peptide_chains(pdb2,d_header,)
                        bmpolymercount2 = d_biomolecules2[biomolecule2]['polymercount']

                        ## continue if transform error
                        if [pdb1,pdb2,str(biomolecule1),str(biomolecule2),] in l_transform_errors:
                            print pdb1, pdb2, biomolecule1, biomolecule2, 'transform error'
                            continue

##                        ## skip if different number of chains in the biomolecule (R350 might be wrong...)
##                        if bmpolymercount1 != bmpolymercount2:
##                            continue

                        ## skip if different hetero compounds (1st check of hetIDs)
                        bool_different, d_hetIDs = self.different_hetero_compounds(pdb1,pdb2,d_header)
                        if bool_different == True:
                            if self.verbose == True:
                                print 'different hetIDs', pdb1, pdb2, d_hetIDs[pdb1], d_hetIDs[pdb2], d_hetIDs[pdb1]^d_hetIDs[pdb2]
                            continue

                        ## skip if different disulfide bonding (e.g. 2huk)
                        if 'SSBOND' in d_header[pdb1].keys() or 'SSBOND' in d_header[pdb2].keys():
                            if not 'SSBOND' in d_header[pdb1].keys():
                                print 'SSBOND', pdb2
                                continue
                            if not 'SSBOND' in d_header[pdb2].keys():
                                print 'SSBOND', pdb1
                                continue
                            if d_header[pdb1]['SSBOND'] == d_header[pdb2]['SSBOND']:
                                pass
                            else:
                                print d_header[pdb1]['SSBOND']
                                print d_header[pdb2]['SSBOND']
                                print pdb1, pdb2
                                stop

                        ## skip if different polymers (other than long peptides)
                        bool_different_polymer_ligands = self.different_nucleotides_saccharides_shortpeptides(pdb1,pdb2,d_header)
                        if bool_different_polymer_ligands == True:
                            continue



                        ##
                        ## identify sequence similar chains between pdbs (long peptides only)
                        ##
                        d_chains_interpdb_sequence_similar = self.identify_similar_chains_from_sequence_inter(
                            d_header,
                            pdb1, pdb2,
                            d_chains_intrapdb_sequence_identical,
                            bmchains1, bmchains2,
                            d_biomolecules1, d_biomolecules2,
                            )

                        ##
                        ## continue if there are no sequence similar chains between the two pdbs
                        ##
                        if d_chains_interpdb_sequence_similar == {}:
                            continue

                        ##
                        ## find chains which are not similar in between pdbs (if any)
                        ## also: are long peptide chains that are not in the *biomolecule( sequence similar to peptide chains that *are* in the biomolecule (e.g. 1i4o)
                        ##
                        (
##                            bmSEQRESchains1_not_similar_to_SEQRESchains2,
##                            bmSEQRESchains2_not_similar_to_SEQRESchains1
                            SEQRESchains1_not_similar_to_SEQRESchains2,
                            SEQRESchains2_not_similar_to_SEQRESchains1
                            ) = self.identify_chains_interpdb_not_sequence_similar(
                                pdb1, pdb2,
##                                bmchains1, bmchains2,
                                long_peptide_chains1,long_peptide_chains2,
                                d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
                                d_header,
                                )

                        ## check that non sequence similar chains (if any)
                        ## 1) are not long peptides and (already done)
                        ## 2) are sequence identical (already done)
                        if len(SEQRESchains1_not_similar_to_SEQRESchains2) > 0 or len(SEQRESchains2_not_similar_to_SEQRESchains1) > 0:
                            if self.verbose == True:
                                print biomolecule1, SEQRESchains1_not_similar_to_SEQRESchains2
                                print biomolecule2, SEQRESchains2_not_similar_to_SEQRESchains1

                            ## 1) check that non sequence similar chains (if any) are not long peptides
                            d_bmSEQRESchains_not_similar_to_SEQRESchains = {
                                pdb1:{'chains1':SEQRESchains1_not_similar_to_SEQRESchains2,'chains2':SEQRESchains2_not_similar_to_SEQRESchains1,'pdb2':pdb2},
                                pdb2:{'chains1':SEQRESchains2_not_similar_to_SEQRESchains1,'chains2':SEQRESchains1_not_similar_to_SEQRESchains2,'pdb2':pdb1},
                                }
                            for pdb in d_bmSEQRESchains_not_similar_to_SEQRESchains.keys():
                                if len(d_bmSEQRESchains_not_similar_to_SEQRESchains[pdb]['chains1']) > 0:
                                    for chain in d_bmSEQRESchains_not_similar_to_SEQRESchains[pdb]['chains1']:
                                        peptide = False
                                        long = False
                                        ## check if peptide
                                        if d_header[pdb]['SEQRES']['chains'][chain]['type'] == 'peptide':
                                            peptide = True
                                        ## check if long chain
                                        if len(d_header[pdb]['SEQRES']['chains'][chain]['seq']) > self.min_len_chain:
                                            long = True
                                        if peptide == True and long == True:
                                            break
                                    if peptide == True and long == True:
                                        break
                            ## continue if long peptide
                            if peptide == True and long == True:
                                continue
                            else:
                                stop_not_expected
                        ## this should be the only check!
                        if len(SEQRESchains1_not_similar_to_SEQRESchains2) > 0 or len(SEQRESchains2_not_similar_to_SEQRESchains1) > 0:
                            stop_not_expected

                        ##
                        ## parse coordinates now that they are needed
                        ##
                        if pdb1 not in d_coordinates.keys():
                            d_coordinates[pdb1], d_hetero[pdb1], d_ATOMseq[pdb1] = quakes_pdb_parser.parse_coordinates(
                                pdb1, d_header[pdb1],
                                path_pdb = self.path_pdb,
                                verbose=verbose,
                                )
                            ## append secondary structure
                            d_ATOMseq[pdb1] = self.append_ss(d_header[pdb1],d_ATOMseq[pdb1],)
                        if pdb2 not in d_coordinates.keys():
                            d_coordinates[pdb2], d_hetero[pdb2], d_ATOMseq[pdb2] = quakes_pdb_parser.parse_coordinates(
                                pdb2, d_header[pdb2],
                                path_pdb = self.path_pdb,
                                verbose=verbose,
                                )
                            ## append secondary structure
                            d_ATOMseq[pdb2] = self.append_ss(d_header[pdb2],d_ATOMseq[pdb2],)

                        print 'check if %s %s and %s %s are identical' %(pdb1,bm1,pdb2,bm2)
                        ## skip if only alpha carbon atoms (e.g. 1thi)
                        alpha1 = self.identify_atom_types_in_long_peptide_chains(d_header[pdb1], d_coordinates[pdb1])
                        if alpha1 == True:
                            print 'alpha', pdb1
                            if not d_header[pdb1]['MDLTYP'] == True:
                                stop1
                            continue
                        alpha2 = self.identify_atom_types_in_long_peptide_chains(d_header[pdb2], d_coordinates[pdb2])
                        if alpha2 == True:
                            print 'alpha', pdb2
                            if not d_header[pdb2]['MDLTYP'] == True:
                                stop2
                            continue

                        ##
                        ## skip if different hetero compounds (2nd check of connectivity)
                        ##
                        different = self.compare_hetero_compounds(
                            d_hetero,d_header,d_coordinates,
                            pdb1,pdb2,bmchains1,bmchains2,
                            d_chains_intrapdb_sequence_identical,d_chains_interpdb_sequence_similar,
                            d_ATOMseq,
                            )
                        if different == True: ## e.g. 1vbo,1vbp
                            print 'different connectivity', pdb1, pdb2, biomolecule1, biomolecule2
                            if verbose == True:
                                for root in d_hetero[pdb1].keys():
                                    print pdb1, root, d_hetero[pdb1][root]
                                for root in d_hetero[pdb2].keys():
                                    print pdb2, root, d_hetero[pdb2][root]
                            continue


                        ##
                        ## identify equivalent chains (interpdb) from structure
                        ##

                        print 'identify equivalent chains between %s %s and %s %s' %(pdb1,bm1,pdb2,bm2)
                        ## identify equivalent chains and calculate rmsd
                        (
                            l_equivalent_chains, rmsd, rmsd_expected, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                            ) = self.identify_interpdb_equivalent_chains_from_structure(
                                pdb1, pdb2,
                                d_chains_intrapdb_sequence_identical,
                                d_chains_interpdb_sequence_similar,
                                d_coordinates, d_header,
                                biomolecule1, biomolecule2,
                                d_biomolecules1, d_biomolecules2,
                                d_ATOMseq,
                                )

                        if rm == None:
                            if rmsd_expected == None:
                                rmsd_expected = 99.9
                            print pdb1,pdb2,'transform error'
                            fd = open('transform_error.txt','a')
                            fd.write(
                                '%4s %4s %1i %1i %2i %2i %10s %10s %.1f\n' %(
                                    pdb1,pdb2,
                                    int(biomolecule1),int(biomolecule2),
                                    len(bmchains1), len(bmchains2),
                                    d_header[pdb1]['CRYST1'].rjust(10),d_header[pdb2]['CRYST1'].rjust(10),
                                    rmsd_expected,
                                    )
                                )
                            fd.close()
                            continue

##                        ##
##                        ## calculate backboneRMSD
##                        ##
##                        tchains1 = l_equivalent_chains[0]
##                        tchains2 = l_equivalent_chains[1]
##                        l_atoms = ['N','CA','C','O']
##                        (
##                            rmsd_backbone, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
##                            ) = self.calculate_rmsd_for_multiple_chains(
##                                tchains1,tchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,
##                                d_chains_interpdb_sequence_similar,
##                                d_chains_intrapdb_sequence_identical,
##                                l_atoms = l_atoms
##                                )

                        ##
                        ## identify number and location of mutations
                        ##
                        n_mutations,d_mutations, lendiff = self.identify_mutations(
                            pdb1, pdb2, l_equivalent_chains,
                            d_chains_interpdb_sequence_similar, d_chains_intrapdb_sequence_identical,
                            )

                        if bool_do_single_mutant == True:
                            if n_mutations == 1:
                                self.phipsi_mutant(
                                    l_equivalent_chains,d_coordinates,
                                    pdb1,pdb2,biomolecule1,biomolecule2,
                                    d_mutations,d_chains_intrapdb_sequence_identical,d_ATOMseq,d_header,
                                    n_chains,
                                    )

                        if bool_sphere == True:
                            if n_chains == 1:
                                if n_mutations == 1:
                                    prefix = 'mutant'
                                    self.sphere(
                                        d_mutations,prefix,
                                        d_header,
                                        d_chains_intrapdb_sequence_identical,
                                        d_chains_interpdb_sequence_similar,
                                        d_coordinates,d_ATOMseq,
                                        l_equivalent_chains,
                                        rmsd, tv1, rm, tv2,
                                        pdb1,pdb2,bm1,bm2,
                                        )

                                elif n_mutations == 0:
                                    prefix = 'wt'
                                    self.sphere(
                                        d_mutations,prefix,
                                        d_header,
                                        d_chains_intrapdb_sequence_identical,
                                        d_chains_interpdb_sequence_similar,
                                        d_coordinates,d_ATOMseq,
                                        l_equivalent_chains,
                                        rmsd, tv1, rm, tv2,
                                        pdb1,pdb2,bm1,bm2,
                                        )

                        ##
                        ## add data to dictionary
                        ##
                        if pdb1 not in d_rmsd_identical.keys():
                            d_rmsd_identical[pdb1] = {}
                        if biomolecule1 not in d_rmsd_identical[pdb1].keys():
                            d_rmsd_identical[pdb1][biomolecule1] = {}
                        if pdb2 not in d_rmsd_identical[pdb1][biomolecule1].keys():
                            d_rmsd_identical[pdb1][biomolecule1][pdb2] = {}
                        if biomolecule2 not in d_rmsd_identical[pdb1][biomolecule1][pdb2].keys():
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2] = {}
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd'] = rmsd
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['chains'] = n_chains
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['residues'] = n_residues+n_mutations
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['coordinates'] = n_coordinates
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['l_equivalent_chains'] = l_equivalent_chains
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['mutations'] = n_mutations
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['transformations'] = transformations
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['lendiff'] = lendiff
                        if len( set(d_header[pdb1]['AUTHOR']) & set(d_header[pdb2]['AUTHOR']) ) > 0:
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['bool_identical_authors'] = True
                        else:
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['bool_identical_authors'] = False
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd4'] = rmsd4
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd8'] = rmsd8
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd16'] = rmsd16
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd32'] = rmsd32

                        d_biomolecules = {
                            pdb1:{'biomolecule':biomolecule1},
                            pdb2:{'biomolecule':biomolecule2},
                            }

                        ## color code structure by rmsd
                        if do_gif_and_pdb == True:
                            quakes_gif.rmsd2bfactor(
                                pdb1, pdb2, biomolecule1, biomolecule2, rmsd,
                                d_coordinates, d_header, tv1, rm, tv2, l_equivalent_chains,
                                bmchains1, bmchains2,
                                d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
                                d_ATOMseq,
                                self.path_pdb,
                                )

                ## save some memory
                if pdb2 in d_coordinates.keys():
                    del d_coordinates[pdb2]

        return d_rmsd_identical, d_header


    def sphere(
        self,d_mutations,prefix,
        d_header,
        d_chains_intrapdb_sequence_identical,
        d_chains_interpdb_sequence_similar,
        d_coordinates,d_ATOMseq,
        l_equivalent_chains,
        rmsd, tv1, rm, tv2,
        pdb1,pdb2,bm1,bm2,
        ):

        chains1 = l_equivalent_chains[0]
        chains2 = l_equivalent_chains[1]
        if len(chains1) != 1 or len(chains2) != 1:
            print chains1
            print chains2
            stop
        if len(chains1) > 1 or len(chains2) > 1:
            print pdb1, chains1
            print pdb2, chains2
            stop
        chain1 = chains1[0]
        chain2 = chains2[0]

        ##
        ## determine residue range
        ##
        rep_chain1 = self.chain2repchain(pdb1,chain1[0],d_chains_intrapdb_sequence_identical,)
        rep_chain2 = self.chain2repchain(pdb2,chain2[0],d_chains_intrapdb_sequence_identical,)
        l1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1']
        l2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2']
        r1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r1']
        r2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r2']
        if l1 == 0 and l2 == 0:
            i_range = range(len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-abs(r2-r1))
        elif l2 == 0 and l1 > 0:
            i_range = range(len(d_ATOMseq[pdb2][chain2[0]]['res_nos'])-l1-r1) 
        elif l1 == 0 and l2 > 0:
            i_range = range(len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-l2-r2) 
        else:
            stop_not_expected
        l_res_indexes = i_range

        ##
        ## determine residue index of mutation
        ##
        if prefix == 'mutant':

            l_mutations1 = d_mutations[pdb1][rep_chain1]
            l_mutation1 = l_mutations1[0]
            res_index1_mutation = l_mutation1[0]
            l_mutations2 = d_mutations[pdb2][rep_chain2]
            l_mutation2 = l_mutations2[0]
            res_index2_mutation = l_mutation2[1]

            if res_index1_mutation-l2 != res_index2_mutation-l1:
                print res_index1_mutation, res_index2_mutation
                print l1,l2
                stop

            if d_ATOMseq[pdb1][chain1[0]]['seq'][res_index1_mutation] == d_ATOMseq[pdb2][chain2[0]]['seq'][res_index2_mutation]:
                print d_ATOMseq[pdb1][chain1[0]]['seq'][res_index1_mutation]
                print d_ATOMseq[pdb2][chain2[0]]['seq'][res_index2_mutation]
                stop_identical

            res_index_mutation = res_index1_mutation-l2
            res_index_mutation = res_index2_mutation-l1

##            res_index1_mutation += l2
##            res_index2_mutation += l1

        ## set radii
        l_radii = [10,20,40,1000]

        ##
        ## calculate dihedrals
        ##
        d_dihedrals = {pdb1:{},pdb2:{},}
        for res_index_dihedral in l_res_indexes:
            res_index1_dihedral = res_index_dihedral+l2
            res_index2_dihedral = res_index_dihedral+l1
            for [pdb,chain,res_index_dihedral_loop,] in [
                [pdb1,chain1[0],res_index1_dihedral,],
                [pdb2,chain2[0],res_index2_dihedral,],
                ]:
                d = self.calculate_dihedrals_of_one_residue(d_coordinates,d_header,d_ATOMseq,pdb,chain,res_index_dihedral_loop,)

                if d == 'multiple_res_names':
                    d_dihedrals[pdb][res_index_dihedral_loop] = None
                else:
                    d_dihedrals[pdb][res_index_dihedral_loop] = {
                        'dihedrals':d['dihedrals'],
                        'phi':d['phi'],
                        'psi':d['psi'],
                        'ss':d['psi'],
                        'ss_prev':d['ss_prev'],
                        'ss_next':d['ss_next'],
                        }

        ##
        ## loop over type of comparison (atoms)
        ##
        (
            coordinates1, coordinates2, residue_count, d_lines, l_RMSDs,
            d_coordinates1, d_coordinates2,
            ) = self.prepare_coords_for_alignment(
                pdb1,pdb2,chains1,chains2,
                d_chains_intrapdb_sequence_identical,
                d_chains_interpdb_sequence_similar,
                d_coordinates,d_ATOMseq,
                l_atoms=[], ## important, otherwise CA
                rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2,
                bool_None_if_unobs = True,
                bool_return_transformed_coordinates = True,
                )

        ## residue missing (REMARK465) (e.g. 3dkf, 3dkg)
        if not res_index1_mutation in d_coordinates1[chain1].keys():
            return

        coord1_mutation = d_coordinates1[chain1][res_index1_mutation]['CA']
        try:
            coord2_mutation = d_coordinates2[chain2][res_index2_mutation]['CA']
        except:
            print d_coordinates2.keys()
            print d_coordinates2[chain2].keys()
            print d_coordinates1[chain1].keys()
            print chain2
            print rep_chain2
            print chains2
            print res_index2_mutation, res_index1_mutation
            coord2_mutation = d_coordinates2[chain2][res_index2_mutation]['CA']
            stop
        coord_mutation = (coord1_mutation+coord2_mutation)/2.

        if len(coordinates1) != len(coordinates2): ## tmp!!!
            print pdb1, pdb2
            print len(coordinates1), len(coordinates2)
            stop

        d_rmsd = {}

        lines_out_all = []
        lines_out_mutsite = []

        ##
        ## loop over sphere centers (residues, wt and mut...)
        ##
        print '    ',
        for res_index_ref in l_res_indexes:
            res_index1_ref = res_index_ref+l2
            res_index2_ref = res_index_ref+l1
            print '\b\b\b\b\b%4i' %(res_index_ref),

            ## residue (center of sphere) is REMARK465
            if not res_index1_ref in d_coordinates1[chain1].keys():
                continue
            if not res_index2_ref in d_coordinates2[chain2].keys():
                continue

            ## CA atom (center of sphere) is REMARK470
            if not 'CA' in d_coordinates1[chain1][res_index1_ref].keys():
                continue
            if not 'CA' in d_coordinates2[chain2][res_index2_ref].keys():
                continue

            coord_ref1 = d_coordinates1[chain1][res_index1_ref]['CA']
            coord_ref2 = d_coordinates2[chain2][res_index2_ref]['CA']

            ## unobserved residue
            if coord_ref1 == None:
                continue
            if coord_ref2 == None:
                continue

            coord_ref = (coord_ref1+coord_ref2)/2.
            dist_mutation = math.sqrt(sum((coord_ref-coord_mutation)**2))
            d_sqdiff = {
                'CA':{},
                'heavy':{},
                'phipsi':{},
                'chi1':{},
                }
            for key in d_sqdiff.keys():
                for radius in l_radii:
                    d_sqdiff[key][radius] = []

            ##
            ## loop over sphere radii (residues)
            ##
##                for res_index in range(len(d_header[pdb1]['SEQRES']['chains'][chain1[0]]['seq'])):
            for res_index in l_res_indexes:

                res_index1 = res_index+l2
                res_index2 = res_index+l1

                ## residue missing (REMARK465)
                if not res_index1 in d_coordinates1[chain1].keys():
                    continue
                if not res_index2 in d_coordinates2[chain2].keys():
                    continue

                phi1 = d_dihedrals[pdb1][res_index1]['phi']
                psi1 = d_dihedrals[pdb1][res_index1]['psi']
                phi2 = d_dihedrals[pdb2][res_index2]['phi']
                psi2 = d_dihedrals[pdb2][res_index2]['psi']
                ss1 = d_ATOMseq[pdb1][chain1[0]]['ss'][res_index1]
                ss2 = d_ATOMseq[pdb2][chain2[0]]['ss'][res_index2]
                ss1_prev = d_dihedrals[pdb1][res_index1]['ss_prev']
                ss2_prev = d_dihedrals[pdb2][res_index2]['ss_prev']
                ss1_next = d_dihedrals[pdb1][res_index1]['ss_next']
                ss2_next = d_dihedrals[pdb2][res_index2]['ss_next']
                if ss1 not in ['HELIX','SHEET','']:
                    print ss1
                    stop
                if (
                    ## random coil with very variable phi psi angles
                    ss1 == '' or ss2 == ''
                    or
                    ss1_prev == '' or ss2_prev == ''
                    or
                    ss1_next == '' or ss2_next == ''
                    or
                    ## atoms missing
                    'N/A' in [phi1,psi1,phi2,psi2,]
                    ):
                    diff_squared_phipsi = 'N/A'
                else:
                    diff_phi = abs(phi1-phi2)
                    diff_psi = abs(psi1-psi2)
                    if diff_phi > 180:
                        diff_phi = 360-diff_phi
                    if diff_psi > 180:
                        diff_psi = 360-diff_psi
                    diff_squared_phipsi = (diff_phi)**2 + (diff_psi)**2

                ## proline or atoms missing
                if not 'chi1' in d_dihedrals[pdb1][res_index1]['dihedrals'].keys():
                    diff_squared_chi1 = 'N/A'
                elif not 'chi1' in d_dihedrals[pdb2][res_index2]['dihedrals'].keys():
                    diff_squared_chi1 = 'N/A'
                else:
                    chi1 = d_dihedrals[pdb1][res_index1]['dihedrals']['chi1']
                    chi2 = d_dihedrals[pdb2][res_index2]['dihedrals']['chi1']
                    diff = abs(chi1-chi2)
                    if diff > 180:
                        diff = 360-diff
                    diff_squared_chi1 = diff**2

                ## loop over atoms of residue
                for atom_name in d_coordinates1[chain1][res_index1].keys():

                    ## only include heavy atoms
                    if atom_name[0] == 'H':
                        continue

                    ## atom missing (REMARK470)
                    if not atom_name in d_coordinates2[chain2][res_index2].keys():
                        continue
                    
                    coord1 = d_coordinates1[chain1][res_index1][atom_name]
                    coord2 = d_coordinates2[chain2][res_index2][atom_name]

                    dist1 = math.sqrt(sum((coord_ref-coord1)**2))
                    dist2 = math.sqrt(sum((coord_ref-coord2)**2))
                    dist = (dist1+dist2)/2.
                    diff_squared_coord = sum((coord1-coord2)**2)
                    for r in l_radii:
                        if dist < r:
                            d_sqdiff['heavy'][r] += [diff_squared_coord]
                            if atom_name == 'CA':
                                d_sqdiff['CA'][r] += [diff_squared_coord]
                                if diff_squared_phipsi != 'N/A':
                                    d_sqdiff['phipsi'][r] += [diff_squared_phipsi]
                                if diff_squared_chi1 != 'N/A':
                                    d_sqdiff['chi1'][r] += [diff_squared_chi1]
                            break
                        ## anything more than 50Angstrom away?
                        if r == l_radii[-1]:
                            print dist
                            stop

                    ## end of loop over atoms of residue

                ## end of loop over sphere radii (residues)

            s_rmsds = ''
            for key in ['CA','heavy','phipsi','chi1',]:
                for r in l_radii:
                    l_sqdiffs = d_sqdiff[key][r]
                    if (
##                        (not r in d_sqdiff[key].keys())
##                        or
                        len(l_sqdiffs) == 0
                        ):
                        s_rmsds += '     None '
                        continue
                    else:
##                        l_sqdiffs = d_sqdiff[key][r]
                        RMSD = math.sqrt(sum(l_sqdiffs)/len(l_sqdiffs))
##                        if (
##                            (r == 5 and RMSD > 9)
##                            or
##                            (r == 10 and RMSD > 18)
##                            or
##                            (r == 20 and RMSD > 20)
##                            ):
##                            print
##                            print r, RMSD
##                            print pdb1, pdb2
##                            print 'res_index'
##                            print res_index1_ref
##                            print res_index2_ref
##                            print res_index_ref
##                            print 'ref coord (sphere center)'
##                            print coord_ref1
##                            print coord_ref2
##                            print 'last coord (sphere radius)'
##                            print coord1
##                            print coord2
##                            print 'mut'
##                            print coord1_mutation
##                            print coord2_mutation
##                            print l_sqdiff
##                            print d_ATOMseq[pdb1][chain1[0]]['seq'][res_index1_ref]
##                            print d_ATOMseq[pdb2][chain2[0]]['seq'][res_index2_ref]
##                            stop_which_pdbs_cause_this_large_an_RMSD
                        s_rmsds += '%9.6f ' %(RMSD)

            ## identical authors
            if (
##                d_header[pdb1]['AUTHOR'] == d_header[pdb2]['AUTHOR'] ## all authors identical
##                or
                len( set(d_header[pdb1]['AUTHOR']) & set(d_header[pdb2]['AUTHOR']) ) >= 1 ## 1 or more identical authors (by name)
                ):
                bool_authors_identical = True
            else:
                bool_authors_identical = False

            line = '%4s %4s %2i %2i %1s %1s %4i %4i %4.1f %s %s\n' %(
                pdb1,pdb2,bm1,bm2,chain1[0],chain2[0],
                res_index_ref,
                res_index_mutation,
                dist_mutation,
                s_rmsds, ## rmsd for each key and each radius
                bool_authors_identical,
                )
            if prefix == 'mutant':
                if res_index_ref == res_index_mutation:
                    lines_out_mutsite += [line]
                lines_out_all += [line]
            else:
                stop

            ## end loop over sphere centers (residues)

##        fd = open('HEWL_%s.txt' %(prefix,),'a') ## mutant or wt
##        fd = open('sphere/%i_%s_%s.txt' %(self.cluster,s_atoms,prefix,),'a') ## mutant or wt
        if prefix == 'mutant':

            fn = 'sphere/out'

            if d_header[pdb1]['CRYST1'] == d_header[pdb2]['CRYST1']:
                fn += '_sameSG'
            else:
                fn += '_diffSG'
##            if bool_authors_identical == True:
##                fn += '_sameauth'
##            else:
##                fn += '_diffauth'

            fd = open('%s_mutsite.txt' %(fn),'a') ## mutant or wt
            fd.writelines(lines_out_mutsite)
            fd.close()

            fd = open('%s_allresidues.txt' %(fn),'a') ## mutant or wt
            fd.writelines(lines_out_all)
            fd.close()

        else:
            stop

        return


    def plot_phipsi_mutants(self,):

        print 'do you want to do phi/psi plots for mutated residues?'
        s = raw_input()
        if s != 'y':
            return

        phipsi_step = 5

        print 'preparing phi/psi python dictionary'
        d_ramachandran = {}
        for res1 in self.d_res1.keys():
            for res2 in self.d_res1.keys():
                if res1 == res2:
                    continue
                d_ramachandran[res1+'_'+res1+res2] = {}
                d_ramachandran[res2+'_'+res1+res2] = {}
                for phi in range(-180,180,phipsi_step,):
                    d_ramachandran[res1+'_'+res1+res2][phi] = {}
                    d_ramachandran[res2+'_'+res1+res2][phi] = {}
                    for psi in range(-180,180,phipsi_step,):
                        d_ramachandran[res1+'_'+res1+res2][phi][psi] = 0
                        d_ramachandran[res2+'_'+res1+res2][phi][psi] = 0
        print 'prepared dictionary'

        lines = []
        for s in '0123456789abcddefghijklmnopqrstuvwxyz':
            fd = open('single_point_mutations/%s.txt' %(s),'r')
            lines += fd.readlines()
            fd.close()

        print 'looping over', len(lines), 'lines'

        for line in lines:

            ## continue if phi/psi angles could not be determined
            if (
                'N/A' in [
                    line.split()[10],
                    line.split()[11],
                    line.split()[12],
                    line.split()[13],
                    ]
                ):
                continue
            res1 = line.split()[8]
            res2 = line.split()[9]
            phi1 = float(line.split()[10])
            psi1 = float(line.split()[11])
            phi2 = float(line.split()[12])
            psi2 = float(line.split()[13])

            if res1 not in self.d_res1.keys() or res2 not in self.d_res1.keys():
                print 'skip', line
                continue

    ##        count = sum_ramachandran(d_ramachandran,res1,phi1,psi1,phipsi_range,)
    ##        if count < count_min:
    ##            print line
    ##            print count
    ##            if res1 != 'G':
    ##                stop1
    ##
    ##        count = sum_ramachandran(d_ramachandran,res2,phi2,psi2,phipsi_range,)
    ##        if count < count_min:
    ##            print line
    ##            print count
    ##            if res2 != 'G':
    ##                stop2

            phi1 = self.round_angle(phi1,phipsi_step)
            psi1 = self.round_angle(psi1,phipsi_step)
            phi2 = self.round_angle(phi2,phipsi_step)
            psi2 = self.round_angle(psi2,phipsi_step)
            d_ramachandran[res1+'_'+res1+res2][phi1][psi1] += 1
            d_ramachandran[res2+'_'+res1+res2][phi2][psi2] += 1

        ##
        ## ramachandran plots
        ##
        for key in d_ramachandran.keys():
            print 'plotting ramachandran', key

            if os.path.isfile('phipsi/plot_phipsi_%s.png' %(key)):
                print key, 'plot'
                os.remove('phipsi/plot_phipsi_%s.png' %(key))
            l_gnuplot = []
            max_count = 0
            sum_count = 0
            for phi in range(-180,180,phipsi_step,):
                for psi in range(-180,180,phipsi_step,):
                     l_gnuplot += ['%s %s %s\n' %(phi,psi,d_ramachandran[key][phi][psi])]
                     if d_ramachandran[key][phi][psi] > max_count:
                         max_count = d_ramachandran[key][phi][psi]
                     sum_count += d_ramachandran[key][phi][psi]
##                l_gnuplot += ['%s %s %s\n' %(phi,psi,d_ramachandran[key][phi][psi])]
                l_gnuplot += ['\n']
##            for psi in range(-180,180,phipsi_step,):
##                l_gnuplot += ['%s %s %s\n' %(phi,psi,d_ramachandran[key][phi][psi])]
            if max_count > 1 and sum_count > 3:
                gnuplot.contour_plot(
                    'plot_phipsi_%s' %(key), l_gnuplot,
                    title='%s' %(key.replace('_',' ')), xlabel='{/Symbol f}', ylabel='{/Symbol y}',
                    x1=-180, x2=180, y1=-180, y2=180,
                    z1=0,
                    )

            ## was a plot generated?
            if not os.path.isfile('plot_phipsi_%s.png' %(key)):
                print 'no plot', key
                continue

            shutil.move(
                'plot_phipsi_%s.png' %(key),
                'phipsi/plot_phipsi_%s.png' %(key),
                )

        return


    def round_angle(self,angle,phipsi_step,):

        if angle == 180.:
            angle = -180.
        else:
            angle = phipsi_step*int(angle/phipsi_step)

        return angle


    def remove_obsolete(self,):

        url = 'ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat'
        print 'removing obsolete files'
        lines = urllib2.urlopen(url).readlines()
        for line in lines[1:]:
            pdb = line.split()[2].lower()

            ## delete pdb
            if os.path.isfile('%s/%s/pdb%s.ent' %(self.path_pdb,pdb[1:3],pdb,)):
                os.remove('%s/%s/pdb%s.ent' %(self.path_pdb,pdb[1:3],pdb,))
                print 'obsolete', pdb

            ## delete htm
            if os.path.isfile('htm/%s.htm' %(pdb,)):
                os.remove('htm/%s.htm' %(pdb,))
                print 'htm', pdb
            else:
                continue

            ## delete gif, pdb
            for x in range(1,3):
                for suffix in ['gif','pdb',]:
                    fn = '%s/%s/%s%02i*' %(suffix,pdb[1],pdb,x,)
##                    os.system('rm -rf %s' %(fn))
                    os.system('rm %s' %(fn))

        return


    def append_ss(self,d_header,d_ATOMseq):

        for chain in d_header['SEQRES']['chains'].keys():
            if chain not in d_ATOMseq.keys():
                continue
            for i in range(1,len(d_ATOMseq[chain]['seq'])):
                res_no = d_ATOMseq[chain]['res_nos'][i]
                iCode = d_ATOMseq[chain]['iCodes'][i]
                seqID = '%4i%1s' %(res_no,iCode,)
                ## check whether helix *and* sheet (remediation)
                if (
                    'HELIX' in d_header.keys()
                    and chain in d_header['HELIX'].keys()
                    and seqID in d_header['HELIX'][chain]
                    and 'SHEET' in d_header.keys()
                    and chain in d_header['SHEET'].keys()
                    and seqID in d_header['SHEET'][chain]
                    ):
                    print chain, seqID
                    stop_helix_and_sheet_same_residue
                if 'HELIX' in d_header.keys() and chain in d_header['HELIX'].keys() and seqID in d_header['HELIX'][chain]:
                    for j in range(i,i+d_header['HELIX'][chain][seqID]):
                        res_no_ss = d_ATOMseq[chain]['res_nos'][j]
                        iCode_ss = d_ATOMseq[chain]['iCodes'][j]
                        d_ATOMseq[chain]['ss'][j] = 'HELIX'
                elif 'SHEET' in d_header.keys() and chain in d_header['SHEET'] and seqID in d_header['SHEET'][chain]:
                    for j in range(i,len(d_ATOMseq[chain]['seq'])):
                        res_no_ss = d_ATOMseq[chain]['res_nos'][j]
                        iCode_ss = d_ATOMseq[chain]['iCodes'][j]
                        d_ATOMseq[chain]['ss'][j] = 'SHEET'
                        if (
                            res_no_ss == d_header['SHEET'][chain][seqID]['res_no']
                            and iCode_ss == d_header['SHEET'][chain][seqID]['iCode']
                            ):
                            break
                else:
                    continue

        return d_ATOMseq


    def phipsi_mutant(
        self,
        l_equivalent_chains,d_coordinates,
        pdb1,pdb2,bm1,bm2,
        d_mutations,d_chains_intrapdb_sequence_identical,d_ATOMseq,d_header,
        n_chains,
        ):

        d_loop = {
            pdb1:{},pdb2:{},
            }
        d_loop[pdb1]['bm'] = bm1
        d_loop[pdb2]['bm'] = bm2
        rep_chain1 = d_mutations[pdb1].keys()[0][0]
        rep_chain2 = d_mutations[pdb2].keys()[0][0]
        res_index1 = d_mutations[pdb1].values()[0][0][0]
        res_index2 = d_mutations[pdb2].values()[0][0][1]
        d_loop[pdb1]['rep_chain'] = rep_chain1
        d_loop[pdb2]['rep_chain'] = rep_chain2
        d_loop[pdb1]['res_index'] = res_index1
        d_loop[pdb2]['res_index'] = res_index2
        for pdb in d_loop.keys():
            rep_chain = d_loop[pdb]['rep_chain']
            res_index = d_loop[pdb]['res_index']
            l_chains = [rep_chain]+d_chains_intrapdb_sequence_identical[pdb][rep_chain]
            for chain in l_chains:
                if chain in d_mutations[pdb].keys():
                    l_mutations = d_mutations[pdb][chain]
                    if len(l_mutations) != 1:
                        print l_mutations
                        stop
                    single_mutation = l_mutations[0]
                    res_no = d_ATOMseq[pdb][chain]['res_nos'][res_index]
                    iCode = d_ATOMseq[pdb][chain]['iCodes'][res_index]
                    d_loop[pdb]['chain'] = chain
                    d_loop[pdb]['res_no'] = res_no
                    d_loop[pdb]['iCode'] = iCode
                    break

        ## skip if modified residue
        if single_mutation[2] == 'X' or single_mutation[3] == 'X':
            return

        chain1 = d_loop[pdb1]['chain']
        chain2 = d_loop[pdb2]['chain']
        res_no1 = d_loop[pdb1]['res_no']
        res_no2 = d_loop[pdb2]['res_no']
        iCode1 = d_loop[pdb1]['iCode']
        iCode2 = d_loop[pdb2]['iCode']

        for pdb in d_loop.keys():
            comment = ''
            if 'SEQADV' in d_header[pdb].keys():
                chain = d_loop[pdb]['chain']
                res_no = d_loop[pdb]['res_no']
                iCode = d_loop[pdb]['iCode']
                if chain in d_header[pdb]['SEQADV'].keys():
                    seq_ID = '%4i%1s' %(res_no,iCode,)
                    if seq_ID in d_header[pdb]['SEQADV'][chain].keys():
                        comment = d_header[pdb]['SEQADV'][chain][seq_ID]['comment']
            d_loop[pdb]['comment'] = comment
        comment1 = d_loop[pdb1]['comment']
        comment2 = d_loop[pdb2]['comment']

        print 'xxxxxxxxx', d_loop[pdb1]['comment'], d_loop[pdb2]['comment']

##        for res_name in self.d_res_names[single_mutation[-2]]:
##            if res_name in d_header[pdb1]['TITLE']:
##                bool_res_name1_in_title1 = True
##            if res_name in d_header[pdb2]['TITLE']:
##                bool_res_name1_in_title2 = True
##        for res_name in self.d_res_names[single_mutation[-1]]:
##            if res_name in d_header[pdb1]['TITLE']:
##                bool_res_name2_in_title1 = True
##            if res_name in d_header[pdb2]['TITLE']:
##                bool_res_name2_in_title2 = True

        bool_single_mutation,pdb_wt,pdb_mutant = self.is_single_mutation(
            single_mutation,d_header,d_coordinates,
            pdb1,pdb2,bm1,bm2,chain1,chain2,res_no1,res_no2,iCode1,iCode2,comment1,comment2,
            )
        if bool_single_mutation == False:
            return

        bm_wt = d_loop[pdb_wt]['bm']
        bm_mutant = d_loop[pdb_mutant]['bm']

        comment = '"%s"' %(comment)

        print 'SEQADV comment', comment

        d_loop = {
            pdb1:{
                'chains_equivalent':list(l_equivalent_chains[0]),
                'phi':{},'psi':{},'chains':[],'res_no':{},'iCode':{},'res_name':{},
                'r':{},'ss':{},'res_name_next':{},'ss_prev':{},'ss_next':{},
                },
            pdb2:{
                'chains_equivalent':list(l_equivalent_chains[1]),
                'phi':{},'psi':{},'chains':[],'res_no':{},'iCode':{},'res_name':{},
                'r':{},'ss':{},'res_name_next':{},'ss_prev':{},'ss_next':{},
                },
            }
##        for pdb in d_loop:
####            d_loop[pdb]['chains'] = d_loop[pdb]['chains']
##            for i in range(len(d_loop[pdb]['chains'])):
##                d_loop[pdb]['chains'][i] = d_loop[pdb]['chains'][i][0]
##            ## get mutated chain
##            rep_chain,mutation = d_mutations[pdb].items()[0]
##            chains = [rep_chain]+d_chains_intrapdb_sequence_identical[pdb][rep_chain]
##            chains_mutated = list(set(d_loop[pdb]['chains']) & set(chains))[0]
####            if len(list(set(d_loop[pdb]['chains']) & set(chains))) != 1:
####                print set(d_loop[pdb]['chains']), set(chains)
####                print d_mutations
####                print pdb, rep_chain, d_chains_intrapdb_sequence_identical[pdb][rep_chain]
####                print d_chains_intrapdb_sequence_identical[pdb]
####                print d_mutations[pdb]
####                stop_not_expected
##            d_loop[pdb]['mutated_chains'] = chains_mutated
##            ## get residue index of mutation
##            res_index = mutation[0][0]
##            d_loop[pdb]['res_index'] = res_index
##            print pdb,mutation

        for pdb in d_loop:

            rep_chain,mutation = d_mutations[pdb].items()[0]
            chains = [rep_chain]+d_chains_intrapdb_sequence_identical[pdb][rep_chain]
            ## first mutation in list of mutations (only 1 mutation!)
            res_index = mutation[0][[pdb1,pdb2,].index(pdb)]

            ## loop over chains sorted equivalently
            for chain in d_loop[pdb]['chains_equivalent']:

                ##
                ## continue if not mutated chain
                ##
                if chain not in chains:
                    continue

                ##
                ## calculate dihedrals
                ##
                d = self.calculate_dihedrals_of_one_residue(d_coordinates,d_header,d_ATOMseq,pdb,chain,res_index,)

                if d == 'multiple_res_names':
                    return

                d_loop[pdb]['ss'][chain] = d['ss'][0]
                d_loop[pdb]['ss_prev'][chain] = d['ss_prev'][0]
                d_loop[pdb]['ss_next'][chain] = d['ss_next'][0]
                for angle in ['phi','psi',]:
                    if d[angle] == 'N/A':
                        d_loop[pdb][angle][chain] = '   N/A'
                    else:
                        d_loop[pdb][angle][chain] = '%6.1f' %(d[angle])
                d_loop[pdb]['r'][chain] = d['r']
                d_loop[pdb]['res_no'][chain] = d['res_no']
                d_loop[pdb]['iCode'][chain] = d['iCode']
                d_loop[pdb]['res_name'][chain] = d['res_name']
                d_loop[pdb]['res_name_next'][chain] = d['res_name_next']
                d_loop[pdb]['chains'] += [d['chain']]
                print pdb,chain,'phi,psi',d['phi'],d['psi']

        ##
        ## read old data and prepare new data
        ##
        lines_new = []

##        for i in range(len(d_loop[pdb1]['chains'])):
        for i in range(len(l_equivalent_chains[0])):
            if pdb1 == pdb_wt:
                chain_wt = l_equivalent_chains[0][i]
                chain_mutant = l_equivalent_chains[1][i]
            else:
                chain_wt = l_equivalent_chains[1][i]
                chain_mutant = l_equivalent_chains[0][i]
            if not (chain_wt in d_loop[pdb_wt]['chains'] and chain_mutant in d_loop[pdb_mutant]['chains']):
                continue

            d_acc = {}
            for pdb,bm,chain in [
                [pdb1,bm1,l_equivalent_chains[0][i],],
                [pdb2,bm2,l_equivalent_chains[1][i],],
                ]:

                ## set res_no
                res_no = d_loop[pdb]['res_no'][chain]
                iCode = d_loop[pdb]['iCode'][chain]

                ## calculate SASA
                if not os.path.isfile('dssp/%s/%s%i.dssp' %(pdb[1:3],pdb,bm,)):
                    s = '/home/tc/Downloads/dsspcmbi %s/%s/%s.pdb%i dssp/%s/%s%i.dssp' %(
                        self.path_biounits, pdb[1:3],pdb,bm, pdb[1:3],pdb,bm,
                        )
                    if self.verbose == True:
                        print s
                    os.system(s)

                ## parse SASA of single residue
                ## ACC = residue water exposed surface in Angstrom**2
                fd = open('dssp/%s/%s%i.dssp' %(pdb[1:3],pdb,bm,),'r')
                lines = fd.readlines()
                fd.close()
                index_l = 1+lines.index('  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n')
                index_s = '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA '.index('ACC')
                acc = None
                for line in lines[index_l:]:
                    ## DSSP: Excessive C to N distance. Chain break residue inserted!!!
                    if line[13] == '!':
                        continue
                    if int(line[5:10]) == res_no and line[11] == chain:
                        acc = float(line[index_s:index_s+3])
                        break
                d_acc[pdb] = acc

            ## residue missing
            if d_acc[pdb1] == None or d_acc[pdb2] == None:
                return

            wt_accessibility = d_acc[pdb_wt]
            mutant_accessibility = d_acc[pdb_mutant]

##            ## calculate sasa for biounit for which no vicinal residues are missing
##            if (
##                d_loop[pdb_mutant]['res_name'][chain_mutant] != 'N/A'
##                and
##                d_loop[pdb_wt]['res_name'][chain_wt] != 'N/A'
##                and
##                (
##                    ## introduction of cavity?!
##                    (wt_accessibility == 0 and mutant_accessibility > 0.475 and d_loop[pdb_wt]['res_name'][chain_wt] not in ['GLY','ALA',] and d_loop[pdb_mutant]['res_name'][chain_mutant] not in ['GLY','ALA',])
##                    or
##                    ## removal of cavity?!
##                    (wt_accessibility > 0.45 and mutant_accessibility == 0 and d_loop[pdb_wt]['res_name'][chain_wt] not in ['GLY','ALA',] and d_loop[pdb_mutant]['res_name'][chain_mutant] not in ['GLY','ALA',])
##                    )
##                ):
##                print pdb1,chain_wt,d_loop[pdb_wt]['res_no'][chain_wt],d_loop[pdb_wt]['res_name'][chain_wt],wt_accessibility
##                print pdb2,chain_mutant,d_loop[pdb_mutant]['res_no'][chain_mutant],d_loop[pdb_mutant]['res_name'][chain_mutant],mutant_accessibility
##                stop

            line_new = '%4s %4s %2i %2i %1s %1s %4i %4i %3s %3s %6s %6s %6s %6s %1i %1i %4.1f %4.1f %1s %1s %3s %3s %1s %1s %1s %1s %2i %s %s\n' %(
                ## identifiers
                pdb_wt,pdb_mutant,
                bm_wt,bm_mutant,
                chain_wt,chain_mutant,
                d_loop[pdb_wt]['res_no'][chain_wt],
                d_loop[pdb_mutant]['res_no'][chain_mutant],
                d_loop[pdb_wt]['res_name'][chain_wt],
                d_loop[pdb_mutant]['res_name'][chain_mutant],
                ## phi/psi angles
                d_loop[pdb_wt]['phi'][chain_wt],
                d_loop[pdb_wt]['psi'][chain_wt],
                d_loop[pdb_mutant]['phi'][chain_mutant],
                d_loop[pdb_mutant]['psi'][chain_mutant],
                ## gauche of chi1 angle
                d_loop[pdb_wt]['r'][chain_wt],
                d_loop[pdb_mutant]['r'][chain_mutant],
                ## SASA
                wt_accessibility, mutant_accessibility,
                ## secondary structure
                d_loop[pdb_wt]['ss'][chain_wt],
                d_loop[pdb_mutant]['ss'][chain_mutant],
                ## next residue
                d_loop[pdb_wt]['res_name_next'][chain_wt],
                d_loop[pdb_mutant]['res_name_next'][chain_mutant],
                ## surrounding secondary structure
                d_loop[pdb_wt]['ss_prev'][chain_wt],
                d_loop[pdb_mutant]['ss_prev'][chain_mutant],
                d_loop[pdb_wt]['ss_next'][chain_wt],
                d_loop[pdb_mutant]['ss_next'][chain_mutant],
                ## number of chains
                n_chains,
                ## space groups identical
                d_header[pdb1]['CRYST1'].rjust(10) == d_header[pdb2]['CRYST1'].rjust(10),
                ## SEQADV
                comment,
                )
            lines_new += [line_new]

        ## include old lines not identical to new lines
        fd = open('single_point_mutations/%s.txt' %(pdb_wt[1],),'r')
        lines_old_redundant = fd.readlines()
        fd.close()
        lines_old_nonredundant = []
        for line in lines_old_redundant:
            if len(line.split()) in [27,28,29,30,31,]:
                if line.split()[:4] != [pdb_wt,pdb_mutant,str(bm_wt),str(bm_mutant),]:
                    lines_old_nonredundant += [line]
            elif len(line.split()) in [25,26,]: ## temp!!!
                print line
                stop
                if line.split()[:2] != [pdb_wt,pdb_mutant,]:
                    line2 = line[:9]+'  1  1'+line[9:]
                    lines_old_nonredundant += [line2]
            else:
##                continue
                print len(line.split())
                print line
                print 'single_point_mutations/%s.txt' %(pdb_wt[1],)
                stop

        fd = open('single_point_mutations/%s.txt' %(pdb_wt[1],),'w')
        fd.writelines(lines_old_nonredundant+lines_new)
        fd.close()

        return


    def calculate_dihedrals_of_one_residue(
        self,
        d_coordinates,d_header,d_ATOMseq,
        pdb,chain,res_index,
        ):

        d_out = {
            }

        ##
        ## get res_no,iCode if mutated chain
        ##
        res_no = d_ATOMseq[pdb][chain]['res_nos'][res_index]
        iCode = d_ATOMseq[pdb][chain]['iCodes'][res_index]
        record = d_ATOMseq[pdb][chain]['records'][res_index]
        ss = d_ATOMseq[pdb][chain]['ss'][res_index]

        r = 0
        d_dihedrals = {}

        ## residue missing
        if record == 'REMARK465':
            phi = 'N/A'
            psi = 'N/A'
            res_name = 'N/A'

        ## atoms missing
        elif (
            'N' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys() or
            'CA' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys() or
            'C' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys()
            ):
            phi = 'N/A'
            psi = 'N/A'
            altloc = min(d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys())
            res_name = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']

        ## atoms of current residues present
        else:
            N = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['N']['coordinate']
            CA = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA']['coordinate']
            C = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['C']['coordinate']
            ## calculate dihedrals
            if res_index-1 >= 0:
                res_no_prev = d_ATOMseq[pdb][chain]['res_nos'][res_index-1]
                iCode_prev = d_ATOMseq[pdb][chain]['iCodes'][res_index-1]
                record_prev = d_ATOMseq[pdb][chain]['records'][res_index-1]
                ## previous residue not present
                if record_prev == 'REMARK465':
                    phi = 'N/A'
                ## atom of previous residue not present
                elif 'C' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms'].keys():
                    phi = 'N/A'
                else:
                    C_prev = d_coordinates[pdb]['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['C']['coordinate']
##                    phi = '%6.1f' %(dihedral.main(C_prev,N,CA,C,))
                    phi = quakes_dihedral.main(C_prev,N,CA,C,)
                    dist = math.sqrt(sum((C_prev-N)**2))
                    if dist > 2.:
                        fd = open('remediation_CprevNdist.txt','a')
                        fd.write('%s %s %s %s\n' %(pdb, chain, res_no, dist))
                        fd.close()
            else:
                phi = 'N/A'
            if res_index+1 < len(d_ATOMseq[pdb][chain]['res_nos']):
                res_no_next = d_ATOMseq[pdb][chain]['res_nos'][res_index+1]
                iCode_next = d_ATOMseq[pdb][chain]['iCodes'][res_index+1]
                record_next = d_ATOMseq[pdb][chain]['records'][res_index+1]
                if record_next == 'REMARK465':
                    psi = 'N/A'
                elif 'N' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms'].keys():
                    psi = 'N/A'
                else:
                    N_next = d_coordinates[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['N']['coordinate']
##                    psi = '%6.1f' %(dihedral.main(N,CA,C,N_next))
                    psi = quakes_dihedral.main(N,CA,C,N_next)
                    dist = math.sqrt(sum((C-N_next)**2))
                    if dist > 2.:
                        fd = open('remediation_CNnextdist.txt','a')
                        fd.write('%s %s %s %s\n' %(pdb, chain, res_no, dist))
                        fd.close()
            else:
                psi = 'N/A'
            l_res_name = []
            for altloc in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                res_name = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']
                l_res_name += [res_name]
            ## return if multiple hetID for resID
            if len(set(l_res_name)) > 1:
                return 'multiple_res_names'

            ## sidechain dihedrals
            if res_name in self.d_dihedrals_atoms.keys():
                l_dihedrals = self.d_dihedrals_atoms[res_name].keys()
                l_dihedrals.sort()
                if self.verbose == True:
                    print res_name
                for s_dihedral in l_dihedrals:
                    if s_dihedral[:4] != 'chi1':
                        continue
                    atom_name1 = self.d_dihedrals_atoms[res_name][s_dihedral][0]
                    atom_name2 = self.d_dihedrals_atoms[res_name][s_dihedral][1]
                    atom_name3 = self.d_dihedrals_atoms[res_name][s_dihedral][2]
                    atom_name4 = self.d_dihedrals_atoms[res_name][s_dihedral][3]
                    Continue = False
                    for atom_name in [atom_name1,atom_name2,atom_name3,atom_name4,]:
                        if (
                            not atom_name in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys()
                            and atom_name in d_header[pdb]['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys()
                            ):
                            Continue = True
                    if Continue == True:
                        continue
                    coord1 = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name1]['coordinate']
                    coord2 = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name2]['coordinate']
                    coord3 = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name3]['coordinate']
                    coord4 = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name4]['coordinate']
                    angle = quakes_dihedral.main(coord1,coord2,coord3,coord4,)
                    d_dihedrals[s_dihedral] = angle
                    if res_name != 'PRO':
                        if angle >= 0 and angle < 120:
                            r = 1
                        elif (angle >= 120 and angle < 180) or (angle < -120 and angle > -180):
                            r = 2
                        elif angle >= -120 and angle < 0:
                            r = 3
                        else:
                            stop
                    elif res_name == 'PRO':
                        if angle >= 0 and angle < 90:
                            r = 1
                        elif angle >= -90 and angle < 0:
                            r = 2
                        else:
                            stop
                    if self.verbose == True:
                        print pdb, chain, res_no, iCode, res_name, s_dihedral, phi, psi, r, angle

        ## is next residue a proline?
        if self.verbose == True:
            print pdb, chain, res_no
        if res_index+1 != len(d_ATOMseq[pdb][chain]['records']) and d_ATOMseq[pdb][chain]['records'][res_index+1] != 'REMARK465':
            res_no_next = d_ATOMseq[pdb][chain]['res_nos'][res_index+1]
            iCode_next = d_ATOMseq[pdb][chain]['iCodes'][res_index+1]
            l_res_name_next = []
            for altloc_next in d_coordinates[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['altlocs'].keys():
                res_name_next = d_coordinates[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['altlocs'][altloc_next]['res_name']
                l_res_name_next += [res_name_next]
            ## multiple hetIDs per resID
            if len(set(l_res_name_next)) > 1:
                res_name_next = 'N/A'
        else:
            res_name_next = 'N/A'

        if res_index+1 != len(d_ATOMseq[pdb][chain]['records']):
            ss_next = d_ATOMseq[pdb][chain]['ss'][res_index+1]
        else:
            ss_next = 'RANDOM'
        ss_prev = d_ATOMseq[pdb][chain]['ss'][res_index-1]
        if ss_prev == '':
            ss_prev = 'RANDOM'
        if ss_next == '':
            ss_next = 'RANDOM'

        ## append and finish loop over chain
        if ss == '':
            ss = 'RANDOM'

        d_out['ss'] = ss
        d_out['ss_prev'] = ss_prev
        d_out['ss_next'] = ss_next
        d_out['phi'] = phi
        d_out['psi'] = psi
        d_out['dihedrals'] = d_dihedrals
        d_out['res_name'] = res_name
        d_out['res_name_next'] = res_name_next
        d_out['r'] = r
        d_out['chain'] = chain
        d_out['res_no'] = res_no
        d_out['iCode'] = iCode

        return d_out


    def is_single_mutation(
        self,single_mutation,d_header,d_coordinates,
        pdb1,pdb2,bm1,bm2,chain1,chain2,res_no1,res_no2,iCode1,iCode2,
        comment1,comment2,
        ):

        ## single point mutation is just a cloning artifact
        if comment1 in self.l_comments_false or comment2 in self.l_comments_false:
            return False,None,None

        ## mutation or modified residue
        elif comment1 == '' and comment2 in self.l_comments_true:
            pdb_wt = pdb1
            pdb_mutant = pdb2
            bm_wt = bm1
            bm_mutant = bm2
            comment = comment2
        elif comment2 == '' and comment1 in self.l_comments_true:
            pdb_wt = pdb2
            pdb_mutant = pdb1
            bm_wt = bm2
            bm_mutant = bm1
            comment = comment1

        ## both proteins mutated at the same site (e.g. 1xgq,1xgp)
        elif (
            '%i%1s' %(res_no1, single_mutation[-2],) in d_header[pdb1]['TITLE']
            and
            '%i%1s' %(res_no2, single_mutation[-1],) in d_header[pdb2]['TITLE']
            ):
            pdb_wt = pdb1
            pdb_mutant = pdb2
            bm_wt = bm1
            bm_mutant = bm2
            fd = open('single_point_mutations/%s.txt' %(pdb_wt[1],),'r')
            lines1 = fd.readlines()
            fd.close()
            lines3 = []
            found = False
            for line in lines1:
                if line.split()[:4] != [pdb_wt,pdb_mutant,str(bm_wt),str(bm_mutant),]:
                    lines3 += [line]
                else:
                    found = True
            if found == True:
                fd = open('deleted_lines.txt','a')
                fd.write('%s %s %s %s\n' %(pdb_wt,pdb_mutant,bm_wt,bm_mutant,))
                fd.close()
                fd = open('single_point_mutations/%s.txt' %(pdb_wt[1],),'w')
                fd.writelines(lines3)
                fd.close()
            return False,None,None
        ## both proteins are mutated at the same site
        ## e.g. 5ptd,6ptd
        elif comment1 in self.l_comments_true+['CONFLICT'] and comment2 in self.l_comments_true+['CONFLICT']:
            return False,None,None

        ## mutation (derived from title)
        ## D52N
        elif '%1s%i%1s' %(single_mutation[-1], res_no1, single_mutation[-2],) in d_header[pdb1]['TITLE']:
            pdb_wt = pdb2
            pdb_mutant = pdb1
            bm_wt = bm2
            bm_mutant = bm1
            comment = comment1
            res_no_mutant = res_no1
            fd = open('remediation_SEQADV_CONFLICT.txt','a')
            fd.write('%4s %4i %s _struct.title %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
            fd.close()
        elif '%1s%i%1s' %(single_mutation[-2], res_no2, single_mutation[-1],) in d_header[pdb2]['TITLE']:
            pdb_wt = pdb1
            pdb_mutant = pdb2
            bm_wt = bm1
            bm_mutant = bm2
            comment = comment2
            res_no_mutant = res_no2
            fd = open('remediation_SEQADV_CONFLICT.txt','a')
            fd.write('%4s %4i %s _struct.title %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
            fd.close()
        ## ASN52 MUTANT
        elif '%3s%i MUTANT' %(self.d_res3[single_mutation[-1]],res_no1,) in d_header[pdb1]['TITLE']:
            pdb_wt = pdb2
            pdb_mutant = pdb1
            bm_wt = bm2
            bm_mutant = bm1
            comment = comment1
            stop_reverse3
            res_no_mutant = res_no1
            fd = open('remediation_SEQADV_CONFLICT.txt','a')
            fd.write('%4s %4i %s _struct.title %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
            fd.close()
        elif '%3s%i MUTANT' %(self.d_res3[single_mutation[-2]],res_no2,) in d_header[pdb2]['TITLE']:
            pdb_wt = pdb1
            pdb_mutant = pdb2
            bm_wt = bm1
            bm_mutant = bm2
            comment = comment2
            stop_reverse4
            res_no_mutant = res_no2
            fd = open('remediation_SEQADV_CONFLICT.txt','a')
            fd.write('%4s %4i %s _struct.title %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
            fd.close()
##        ## ASP TO ASN
##        elif '%3s TO %3s' %(self.d_res3[single_mutation[-2]],self.d_res3[single_mutation[-1]],) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            stop_reverse5
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif '%3s TO %3s' %(self.d_res3[single_mutation[-1]],self.d_res3[single_mutation[-1]],) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            stop_reverse6
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ASP52
##        elif comment2 == '' and '%3s%i' %(self.d_res3[single_mutation[-1]],res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and '%3s%i' %(self.d_res3[single_mutation[-2]],res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ASP 52
##        elif comment2 == '' and '%3s %i' %(self.d_res3[single_mutation[-2]],res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            stop_reverse7
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and '%3s %i' %(self.d_res3[single_mutation[-1]],res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            stop_reverse8
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ASPARTIC ACID 52
##        elif comment2 == '' and '%s %i' %(self.d_res_names[single_mutation[-1]],res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and '%s %i' %(self.d_res_names[single_mutation[-2]],res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ' 52 '
##        elif comment2 == '' and ' %i ' %(res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and ' %i ' %(res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ' MUTA*'
##        elif comment2 == '' and ' MUTA' in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and ' MUTA' in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('remediation_SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()

        elif comment1 in ['','SEE REMARK 999',] and comment2 in ['','SEE REMARK 999',]:

            print chain1, res_no1
            print chain2, res_no2
            print comment1
            print comment2

            ## different organisms, both could be wt (e.g. 3fwq,2r7i)
            mol_id1 = d_header[pdb1]['COMPND'][chain1]
            mol_id2 = d_header[pdb2]['COMPND'][chain2]
            if d_header[pdb1]['SOURCE'][mol_id1] != d_header[pdb2]['SOURCE'][mol_id2]:
                return False,None,None

            ## temporary...
            elif res_no1 != res_no2 and not (res_index1 == 0 or res_index2 == 0):
                print comment1, comment2
                print '%1s%i%1s' %(single_mutation[-1], res_no1, single_mutation[-2],)
                print '%1s%i%1s' %(single_mutation[-1], res_no2, single_mutation[-2],)
                print chain1, chain2
                print d_header[pdb1]['SOURCE'], d_header[pdb2]['SOURCE']
                print res_index1, res_index2
                stop_res_no_example_curious

            ## can't determine which pdb is the wildtype (if any!)
            else:
                print d_header[pdb1]['TITLE']
                print d_header[pdb2]['TITLE']
                print '%1s%i%1s' %(single_mutation[-1], res_no1, single_mutation[-2],)
                print '%1s%i%1s' %(single_mutation[-2], res_no2, single_mutation[-1],)
                fd = open('remediation_SEQADV_missing.txt','a')
                fd.write('%s %s %s %s %s %s %s --- %s --- %s\n' %(pdb1,pdb2,chain1,chain2,res_no1,res_no2,single_mutation,d_header[pdb1]['TITLE'],d_header[pdb2]['TITLE'],))
                fd.close()
                return False,None,None

        ## heavy and light chain of variant antibody
        elif chain1 == chain2 == 'H' or chain1 == chain2 == 'L':
            return False,None,None
        elif 'MHC' in d_header[pdb1]['TITLE'] and 'MHC' in d_header[pdb2]['TITLE']:
            return False,None,None

##        ## alanine mutation
##        elif (
##            comment1 == 'CONFLICT' and comment2 == ''
##            and
##            'ALA' == d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['altlocs'][' ']['res_name']
##            and
##            single_mutation[-1] == 'A'
##            ):
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##        elif (
##            comment2 == 'CONFLICT' and comment1 == ''
##            and
##            'ALA' == d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['altlocs'][' ']['res_name']
##            and
##            single_mutation[-1] == 'A'
##            ):
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2

        ## spontaneous (de)amidation?
        elif (
            (single_mutation[-1] == 'D' and single_mutation[-2] == 'N') ## deamidation
            or
            (single_mutation[-1] == 'N' and single_mutation[-2] == 'D') ## amidation
            ):
            return False,None,None
        

        ## not expected
        else:
            bool_mutation = False
            for pdb,chain,res_no,res_wt,res_mutant in [
                [pdb1,chain1,res_no1,single_mutation[-1],single_mutation[-2],],
                [pdb2,chain2,res_no2,single_mutation[-2],single_mutation[-1],],
                ]:
                bool_mutation, s = quakes_pdb_parser.parse_mmCIF_mutation(
                    pdb,chain,res_no, res_wt,res_mutant,
                    self.path_mmCIF,
                    )
                ## break pdb loop
                if bool_mutation == True:
                    break

            if bool_mutation == True:
                if pdb == pdb1:
                    pdb_wt = pdb2
                    pdb_mutant = pdb1
                    bm_wt = bm2
                    bm_mutant = bm1
                    comment = comment1
                    fd = open('remediation_SEQADV_CONFLICT.txt','a')
                    fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no2,comment,s,))
                    fd.close()
                elif pdb == pdb2:
                    pdb_wt = pdb1
                    pdb_mutant = pdb2
                    bm_wt = bm1
                    bm_mutant = bm2
                    comment = comment2
                    fd = open('remediation_SEQADV_CONFLICT.txt','a')
                    fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no1,comment,s,))
                    fd.close()
            else:
                print d_header[pdb1]['TITLE']
                print d_header[pdb2]['TITLE']
                try:
                    print pdb1, chain1, res_no1, d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['altlocs'][' ']['res_name'], comment1
                except:
                    print pdb1, chain1, res_no1, d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['altlocs']['A']['res_name'], comment1
                try:
                    print pdb2, chain2, res_no2, d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['altlocs'][' ']['res_name'], comment2
                except:
                    print pdb2, chain2, res_no2, d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['altlocs']['A']['res_name'], comment2
                fd = open('remediation_SEQADV_unknownwhichismutant.txt','a')
                fd.write(
                    '%4s %4s %1s %1s %4i %4i %s %s %s\n' %(
                        pdb1,pdb2,
                        chain1,chain2,
                        res_no1,res_no2,
                        comment1,comment2,
                        single_mutation,
                        )
                    )
                fd.close()
                return False,None,None

        return True,pdb_wt,pdb_mutant


    def compare_hetero_compounds(
        self,d_hetero,d_header,d_coordinates,
        pdb1,pdb2,bmchains1,bmchains2,
        d_chains_intrapdb_sequence_identical,d_chains_interpdb_sequence_similar,
        d_ATOMseq,
        ):

        different = False
        pdbpairs = [[pdb1,pdb2],[pdb2,pdb1]]
        for pdbpair in pdbpairs:
            pdba = pdbpair[0]
            pdbb = pdbpair[1]
            rootsa = d_hetero[pdba].keys()
            for roota in rootsa:

                chaina = roota[0]
                ## check peptide connections only
##                if chaina not in d_header[pdba]['SEQRES']['chains'].keys():
##                    continue
##                if d_header[pdba]['SEQRES']['chains'][chaina]['type'] != 'peptide':
##                    continue
                rep_chaina = self.chain2repchain(pdba,chaina,d_chains_intrapdb_sequence_identical,)
                if rep_chaina not in d_chains_interpdb_sequence_similar.keys():
                    continue

                ## continue if root is not in biounit
                if (
                    (pdba == pdb1 and roota[0] != ' ' and roota[0] not in bmchains1) or
                    (pdba == pdb2 and roota[0] != ' ' and roota[0] not in bmchains2)
                    ):
                    continue

                ## continue if root is a short nucleotide (to avoid reversly numbered nucleotides etc.)
                if roota[0] in d_header[pdba]['chains'].keys() and d_header[pdba]['chains'][roota[0]]['type'] == 'nucleotide' and len(d_header[pdba]['chains'][roota[0]]['seq']) <= 3:
                    continue

                compounda = d_hetero[pdba][roota]
                identical = False
                rootsb = d_hetero[pdbb].keys()
                for rootb in rootsb:

                    chainb = rootb[0]
                    ## check peptide connections only
##                    if chainb not in d_header[pdbb]['SEQRES']['chains'].keys():
##                        continue
##                    if d_header[pdba]['SEQRES']['chains'][chaina]['type'] != 'peptide':
##                        continue
                    rep_chainb = self.chain2repchain(pdbb,chainb,d_chains_intrapdb_sequence_identical,)
                    if rep_chainb not in d_chains_interpdb_sequence_similar[rep_chaina].keys():
                        continue

                    ## continue if root is not in biounit
                    if (
                        (pdbb == pdb1 and rootb[0] != ' ' and rootb[0] not in bmchains1) or
                        (pdbb == pdb2 and rootb[0] != ' ' and rootb[0] not in bmchains2)
                        ):
                        continue

                    ## continue if root is a short nucleotide (to avoid reversly numbered nucleotides etc.)
                    if rootb[0] in d_header[pdbb]['chains'].keys() and d_header[pdbb]['chains'][rootb[0]]['type'] == 'nucleotide' and len(d_header[pdbb]['chains'][rootb[0]]['seq']) <= 3:
                        continue

                    compoundb = d_hetero[pdbb][rootb]

                    ## identical compounds, but also identical roots?
                    if compounda == compoundb:

                        ## hetero compounds not connected to peptide
                        if (
                            d_coordinates[pdba]['chains'][roota[0]]['residues'][int(roota[1:5])]['d_iCodes'][roota[-1]]['record'] == 'HETATM' and
                            d_coordinates[pdbb]['chains'][rootb[0]]['residues'][int(rootb[1:5])]['d_iCodes'][rootb[-1]]['record'] == 'HETATM'
                            ):
                            identical = True
                            break

                        if not (
                            chaina in d_ATOMseq[pdba].keys()
                            and
                            chainb in d_ATOMseq[pdba].keys()
                            ):

                            pass

                        ## e.g. 1ebk 2aoe
                        elif (
                            len(d_ATOMseq[pdba][chaina]['seq']) < 50
                            and
                            len(d_ATOMseq[pdbb][chainb]['seq']) < 50
                            ):

                            la = 0
                            lb = 0
                            ra = 0
                            rb = 0

                        ## hetero compounds connected to peptide
                        else:

                            res_indexa = d_ATOMseq[pdba][chaina]['indexes'].index(roota)
                            res_indexb = d_ATOMseq[pdbb][chainb]['indexes'].index(rootb)
 
                            if pdb1 == pdba and pdb2 == pdbb:
                                la = d_chains_interpdb_sequence_similar[rep_chaina][rep_chainb]['l1']
                                lb = d_chains_interpdb_sequence_similar[rep_chaina][rep_chainb]['l2']
                                ra = d_chains_interpdb_sequence_similar[rep_chaina][rep_chainb]['r1']
                                rb = d_chains_interpdb_sequence_similar[rep_chaina][rep_chainb]['r2']
                            elif pdb1 == pdbb and pdb2 == pdba:
                                la = d_chains_interpdb_sequence_similar[rep_chainb][rep_chaina]['l2']
                                lb = d_chains_interpdb_sequence_similar[rep_chainb][rep_chaina]['l1']
                                ra = d_chains_interpdb_sequence_similar[rep_chainb][rep_chaina]['r2']
                                rb = d_chains_interpdb_sequence_similar[rep_chainb][rep_chaina]['r1']
                            
                            if res_indexa+la != res_indexb+lb: ## e.g. 2jbj,2pvw;1z4v,1z4w,1z4z;1d6q,1re2;1o8a,2iux
                                identical = False
                                continue ## loop over rootsb
                            if d_header[pdba]['SEQRES']['chains'][chaina]['seq'][res_indexa] != d_header[pdbb]['SEQRES']['chains'][chainb]['seq'][res_indexb]:
                                stop_expected3

                        ## if not different, then identical
                        identical = True
                        break

                if identical == False:
                    print 'different', roota
                    print compounda
                    different = True
                    break
            if different == True:
                print rootsa, rootsb
                different = True
                fd = open('CONECT.txt','a')
                fd.write('%s %s\n%s\n%s\n' %(pdb1, pdb2, d_header[pdb1]['TITLE'], d_header[pdb2]['TITLE']))
                fd.close()
                break

        return different


    def identify_atom_types_in_long_peptide_chains(self, d_header, d_coordinates):

        atoms = set()
        for chain in d_header['SEQRES']['chains'].keys():
            if d_header['SEQRES']['chains'][chain]['type'] != 'peptide' or len(d_header['SEQRES']['chains'][chain]['seq']) < self.min_len_chain:
                continue
            for res_no in d_coordinates['chains'][chain]['residues'].keys():
                for iCode in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'].keys():
                    if 'REMARK' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                        continue
                    if 'HETATM' == d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record']:
                        continue
                    for atom in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        if 'REMARK' not in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom].keys():
                            atoms |= set([atom])
                    if len(atoms-set(['CA'])) > 0:
                        return False
        if atoms != set(['CA']):
            print atoms
            print d_header['SEQRES']['chains'][chain]['type']
            notexected
        return True


    def identify_mutations(self, pdb1, pdb2, l_equivalent_chains, d_chains_interpdb_sequence_similar, d_chains_intrapdb_sequence_identical):

        lendiff = False

        d_mutations = {pdb1:{},pdb2:{}}
        n_mutations = 0
        chains1 = l_equivalent_chains[0]
        chains2 = l_equivalent_chains[1]
        n_chains = len(chains1)
        for i in range(n_chains):
            chain1 = chains1[i][0]
            chain2 = chains2[i][0]
            d_chains = {pdb1:chain1,pdb2:chain2}
            for pdb in d_chains.keys():
                chain = d_chains[pdb]
                rep_chain = self.chain2repchain(pdb,chain,d_chains_intrapdb_sequence_identical,)
                d_chains[pdb] = rep_chain
            rep_chain1 = d_chains[pdb1]
            rep_chain2 = d_chains[pdb2]
            n_mutations_per_chain = len(d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l_mutations'])
            n_mutations += n_mutations_per_chain
            
            if len(d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l_mutations']) > 0:
                d_mutations[pdb1][rep_chain1] = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l_mutations']
                d_mutations[pdb2][rep_chain2] = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l_mutations']

            if (
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1'] != 0 or
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2'] != 0 or
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r1'] != 0 or
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r2'] != 0
                ):
                lendiff = True

        return n_mutations, d_mutations, lendiff


    def chain2repchain(self,pdb,chain,d_chains_intrapdb_sequence_identical,):

        for rep_chain in d_chains_intrapdb_sequence_identical[pdb]:
            if chain == rep_chain or chain in d_chains_intrapdb_sequence_identical[pdb][rep_chain]:
                return rep_chain

        return


    def different_nucleotides_saccharides_shortpeptides(self,pdb1,pdb2,d_header):

##        print 'identifying different polymers (not long peptides)'

        d_chains = {
            pdb1:{'peptide':set(),'saccharide':set(),'nucleotide':set(),'N/A':set(),},
            pdb2:{'peptide':set(),'saccharide':set(),'nucleotide':set(),'N/A':set(),},
            }
        for pdb in d_chains.keys():
            for chain in d_header[pdb]['SEQRES']['chains'].keys():
                chain_type = d_header[pdb]['SEQRES']['chains'][chain]['type']
                chain_seq = d_header[pdb]['SEQRES']['chains'][chain]['seq']
                chain_len = len(chain_seq)
                if chain_type != 'peptide' or chain_len < self.min_len_chain:
                    d_chains[pdb][chain_type] |= set([chain_seq])
        for polymer in ['peptide','saccharide','nucleotide']:
            if len(d_chains[pdb1][polymer] ^ d_chains[pdb2][polymer]) > 0:
                return True
        
        return False


    def different_hetero_compounds(self,pdb1,pdb2,d_header):

        ##
        ## skip if different hetero compounds
        ##
        d_hetIDs = {
            pdb1:set(),
            pdb2:set(),
            }
        ## parse hetIDs from HET records
        for pdb in d_hetIDs:
            for chain in d_header[pdb]['HET'].keys():
                for res_no in d_header[pdb]['HET'][chain].keys():
                    for iCode in d_header[pdb]['HET'][chain][res_no].keys():
                        d_hetIDs[pdb] |= d_header[pdb]['HET'][chain][res_no][iCode]
            d_hetIDs[pdb] -= set(self.d_modres.keys()) ## modified residues
        ## compare hetero compounds assuming the following two statements to be correct
        ## 1) "A particular HET group is represented in the PDB archives with a *unique* hetID."
        ## 2) Depositors specify *all* hetero atoms observed in the electron density map.
        ## ignore ions (not important) and ignore saccharides (treated elsewhere)
        if (
            d_hetIDs[pdb1]-set(self.d_ions.keys())-set(self.d_saccharides.keys())-set(self.l_solutes)
            !=
            d_hetIDs[pdb2]-set(self.d_ions.keys())-set(self.d_saccharides.keys())-set(self.l_solutes)
            ):
            return True, d_hetIDs

        ## skip if different saccharides (accept different stereoisomers e.g. 1d7b.pdb,1d7d.pdb)
        saccharides1_stereo = set(self.d_saccharides.keys()) & d_hetIDs[pdb1]
        saccharides2_stereo = set(self.d_saccharides.keys()) & d_hetIDs[pdb2]
        saccharides1 = set()
        saccharides2 = set()
        for saccharide in saccharides1_stereo:
            saccharides1 |= set([self.d_saccharides[saccharide]['stereo']])
        for saccharide in saccharides2_stereo:
            saccharides2 |= set([self.d_saccharides[saccharide]['stereo']])
        if saccharides1 != saccharides2:
            return True, d_hetIDs

##        ## skip if different charges (e.g. 1ext.pdb, 1ncf.pdb but also 1ncy,1ncz)
##        ions1 = set(self.d_ions.keys()) & d_hetIDs[pdb1]
##        ions2 = set(self.d_ions.keys()) & d_hetIDs[pdb2]
##        plus1 = 0
##        minus1 = 0
##        for ion in ions1:
##            formula = self.d_ions[ion][0]
##            if len(formula.split()) == 1:
##                charge = self.d_ions[ion][1]
##                if charge > 0:
##                    plus1 += charge
##                if charge < 0:
##                    minus1 += charge
##        plus2 = 0
##        minus2 = 0
##        for ion in ions2:
##            formula = self.d_ions[ion][0]
##            if len(formula.split()) == 1:
##                charge = self.d_ions[ion][1]
##                if charge > 0:
##                    plus2 += charge
##                if charge < 0:
##                    minus2 += charge
##        if plus1 != plus2 or minus1 != minus2:
##            return True

        return False, d_hetIDs


    def pdbskip(self, d_header, pdb):

        pdbskip = False
        SEQRESchains = d_header[pdb]['SEQRES']['chains'].keys()

        ## continue if rerefined structure
        if 'REMARK0' in d_header[pdb].keys():
            pdbskip = True
            return pdbskip

        ## continue if structure split in multiple files
        if 'SPLIT' in d_header[pdb].keys():
            pdbskip = True
            return pdbskip

##        ## continue if alpha carbon atoms only (need to be specific chains...)
##        if 'MDLTYP' in d_header[pdb].keys():
##            pdbskip = True
##            return pdbskip

        ## continue if no polymer chains (e.g. 1qd8.pdb)
        if SEQRESchains == []:
            pdbskip = True
            return pdbskip

        ## continue if no (long) protein chains
        if d_header[pdb]['proteinchains'] == []:
            pdbskip = True
            return pdbskip

        ## continue if NMR or EM structure
        if d_header[pdb]['EXPDTA'] != 'X-RAY':
            pdbskip = True
            return pdbskip

        ## continue if multiple models (NMR or non-NCS-averaged x-ray structures - e.g. 1htq.pdb)
        if 'MODEL' in d_header[pdb].keys():
            pdbskip = True
            return pdbskip

        ## continue if low resolution (e.g. BAD EXAMPLE 1jgp.pdb, 2w0c)
        if d_header[pdb]['REMARK2'] != 'N/A':
            if d_header[pdb]['REMARK2'] > self.minres:
                ## e.g. 1lzh,2qij
                pdbskip = True
                return pdbskip

        return pdbskip


    def find_long_peptide_chains(self, pdb, d_header,):

        l_long_peptide_chains = []
        for chain in d_header[pdb]['SEQRES']['chains'].keys():
            chain_type = d_header[pdb]['SEQRES']['chains'][chain]['type']
            chain_seq = d_header[pdb]['SEQRES']['chains'][chain]['seq']
            chain_len = len(chain_seq)
            if chain_type == 'peptide' and chain_len >= self.min_len_chain:
                l_long_peptide_chains += [chain]

        return l_long_peptide_chains
     

    def identify_biomolecule(self, pdb, d_header):

        SEQRESchains = d_header[pdb]['SEQRES']['chains'].keys()
        if 'REMARK350' in d_header[pdb].keys():
            d_biomolecules = {}
            biomolecules = d_header[pdb]['REMARK350'].keys()
            for biomolecule in biomolecules:
                chains = d_header[pdb]['REMARK350'][biomolecule]['chains'].keys()
                d_biomolecules[biomolecule] = {}
                d_biomolecules[biomolecule]['chains'] = []
                d_biomolecules[biomolecule]['polymercount'] = 0
                for chain in chains:
                    d_biomolecules[biomolecule]['chains'] += [chain]
                    ## only count if not hetero (e.g. not water not mentioned in remark525 records)
                    if chain in SEQRESchains:
                        d_biomolecules[biomolecule]['polymercount'] += len(d_header[pdb]['REMARK350'][biomolecule]['chains'][chain])

        ## assume everything to be the biomolecule
        else:
##        if d_biomolecules == {}:
            d_biomolecules = {
                1:{
                    'chains':SEQRESchains,
                    'polymercount':len(SEQRESchains),
                    }
                }
            d_header[pdb]['REMARK350'] = {1: {'chains': {}, 'matrices': {1: [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]},}}
            for chain in SEQRESchains:
                d_header[pdb]['REMARK350'][1]['chains'][chain] = set([1])

        return d_biomolecules


    def identify_chains_interpdb_not_sequence_similar(
        self, pdb1, pdb2,
##        bmchains1, bmchains2,
        long_peptide_chains1,long_peptide_chains2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_header,
        ):

##        print 'identifying different long peptides'

        ## pdb1
        SEQRESchains1_similar_to_SEQRESchains2 = set()
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            ## e.g. 1nvw.pdb,1nvv.pdb
            SEQRESchains1_similar_to_SEQRESchains2 |= set([rep_chain1])
            SEQRESchains1_similar_to_SEQRESchains2 |= set(d_chains_intrapdb_sequence_identical[pdb1][rep_chain1])
##        bmSEQRESchains1_similar_to_SEQRESchains2 = SEQRESchains1_similar_to_SEQRESchains2 & set(bmchains1)

        ## pdb2
        SEQRESchains2_similar_to_SEQRESchains1 = set()
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            rep_chains2 = d_chains_interpdb_sequence_similar[rep_chain1].keys()
            ## e.g. 1nvw.pdb,1nvv.pdb
            SEQRESchains2_similar_to_SEQRESchains1 |= set(rep_chains2)
            for rep_chain2 in rep_chains2:
                SEQRESchains2_similar_to_SEQRESchains1 |= set(d_chains_intrapdb_sequence_identical[pdb2][rep_chain2])
##        bmSEQRESchains2_similar_to_SEQRESchains1 = SEQRESchains2_similar_to_SEQRESchains1 & set(bmchains2)

        if set(d_header[pdb1]['SEQRES']['chains'].keys()) & (set(long_peptide_chains1) - SEQRESchains1_similar_to_SEQRESchains2) != set(long_peptide_chains1) - SEQRESchains1_similar_to_SEQRESchains2:
            stop1
        if set(d_header[pdb2]['SEQRES']['chains'].keys()) & (set(long_peptide_chains2) - SEQRESchains2_similar_to_SEQRESchains1) != set(long_peptide_chains2) - SEQRESchains2_similar_to_SEQRESchains1:
            stop2
##        SEQRESchains1_not_similar_to_SEQRESchains2 = set(d_header[pdb1]['SEQRES']['chains'].keys()) & (set(long_peptide_chains1) - SEQRESchains1_similar_to_SEQRESchains2)
##        SEQRESchains2_not_similar_to_SEQRESchains1 = set(d_header[pdb2]['SEQRES']['chains'].keys()) & (set(long_peptide_chains2) - SEQRESchains2_similar_to_SEQRESchains1)
        SEQRESchains1_not_similar_to_SEQRESchains2 = set(long_peptide_chains1) - SEQRESchains1_similar_to_SEQRESchains2
        SEQRESchains2_not_similar_to_SEQRESchains1 = set(long_peptide_chains2) - SEQRESchains2_similar_to_SEQRESchains1

        return SEQRESchains1_not_similar_to_SEQRESchains2, SEQRESchains2_not_similar_to_SEQRESchains1


    def calculate_rmsd_for_multiple_chains(
        self,
        chains1,chains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,
        d_header,
        d_chains_interpdb_sequence_similar,
        d_chains_intrapdb_sequence_identical,
        d_ATOMseq,
        verbose=True,
        l_atoms = ['CA'],
        ):

        instance_geometry = geometry.geometry()

        if len(chains1) != len(chains2):
            print pdb1, pdb2, biomolecule1, biomolecule2
            if len(chains1) < 60 or len(chains2) < 60:
                print chains1, chains2
            print len(chains1), len(chains2)
            notexpected

        ##
        ## parse coordinates for the specific combinations of chain1 and chain2
        ##
        (
            coordinates1, coordinates2, residue_count, d_lines, l_RMSDs,
            d_coordinates1, d_coordinates2,
            ) = self.prepare_coords_for_alignment(
            pdb1,pdb2,chains1,chains2,
            d_chains_intrapdb_sequence_identical,
            d_chains_interpdb_sequence_similar,
            d_coordinates,d_ATOMseq,
            l_atoms=l_atoms,
            )

        if coordinates1 == None and coordinates2 == None:
            return None, None, None, None, None, None, None

        rmsd = instance_geometry.superpose(coordinates1,coordinates2)
            
        if rmsd == 0:
            print pdb1, pdb2
            notexpected_structuresidentical
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter

##        ##
##        ## redo rmsd calculation excluding deviating parts of the structure
##        ##
##        if rmsd < self.max_rmsd and rmsd > 2.25: ## 2.25 was empirically chosen
##            coordinates1, coordinates2, residue_count, d_lines, l_RMSDs = self.prepare_coords_for_alignment(
##                pdb1,pdb2,chains1,chains2,
##                d_chains_intrapdb_sequence_identical,
##                d_chains_interpdb_sequence_similar,
##                d_coordinates,d_ATOMseq,
##                l_atoms=l_atoms,
##                rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2,
##                )
##
##            rmsd = instance_geometry.superpose(coordinates1,coordinates2)
##            if rmsd == 0:
##                print pdb1, pdb2
##                notexpected_structuresidentical
##            tv1 = instance_geometry.fitcenter
##            rm = instance_geometry.rotation
##            tv2 = instance_geometry.refcenter
        
        if verbose == True:
            if len(chains1) < 60 and len(chains2) < 60:
                if self.verbose == True:
                    print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s, %s, %s' %(
                        rmsd, len(chains1), residue_count, len(coordinates1), pdb1, pdb2, biomolecule1, biomolecule2,
                        )
                    print chains1, chains2
            else:
                print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s' %(rmsd, len(chains1), residue_count, len(coordinates1), pdb1, pdb2)

        return rmsd, len(chains1), residue_count, len(coordinates1), tv1, rm, tv2


    def prepare_coords_for_alignment(
        self,pdb1,pdb2,chains1,chains2,
        d_chains_intrapdb_sequence_identical,
        d_chains_interpdb_sequence_similar,
        d_coordinates,d_ATOMseq,
        ## only use alpha carbon atoms for structural alignment for speed purposes
        l_atoms=['CA'],
        rmsd=None, tv1=None, rm=None, tv2=None,
        bool_None_if_unobs = False,
        bool_return_transformed_coordinates = False,
        ):

        d_lines = {pdb1:{},pdb2:{},}

        coordinates1 = []
        coordinates2 = []
        d_coordinates1 = {}
        d_coordinates2 = {}
        residue_count = 0

        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            rep_chain1 = self.chain2repchain(pdb1,chain1[0],d_chains_intrapdb_sequence_identical,)
            rep_chain2 = self.chain2repchain(pdb2,chain2[0],d_chains_intrapdb_sequence_identical,)

            if not rep_chain2  in d_chains_interpdb_sequence_similar[rep_chain1].keys():
                ## e.g. 1bgy 1be3
                rmsd = 'N/A'
                return None, None, None, d_lines, None

            l1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1']
            l2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2']
            r1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r1']
            r2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r2']

            (
                coords1, coords2, rescount, lines1, lines2, l_RMSDs,
                d_coordinates1[chain1], d_coordinates2[chain2],
                ) = self.ATOMrecords2coordinates(
                    d_coordinates, pdb1, pdb2, chain1, chain2,
                    l1, l2, r1, r2, d_ATOMseq,
                    rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2,
                    l_atoms=l_atoms,
                    bool_None_if_unobs = bool_None_if_unobs,
                    bool_return_transformed_coordinates = bool_return_transformed_coordinates,
                    )

            ## append coordinates
            coordinates1 += coords1
            coordinates2 += coords2
            residue_count += rescount

            ##
            ## append lines by model to dictionary
            ##
            if len(chain1) == 1:
                model = 'wt'
            else:
                model = chain1[2:]
            if not model in d_lines[pdb1].keys():
                d_lines[pdb1][model] = []
            d_lines[pdb1][model] += lines1

            if len(chain2) == 1:
                model = 'wt'
            else:
                model = chain2[2:]
            if not model in d_lines[pdb2].keys():
                d_lines[pdb2][model] = []
            d_lines[pdb2][model] += lines2

        return coordinates1, coordinates2, residue_count, d_lines, l_RMSDs, d_coordinates1, d_coordinates2


    def ATOMrecords2coordinates(self,
        d_coordinates, pdb1, pdb2, chain1, chain2,
        l1, l2, r1, r2, d_ATOMseq,
        ## transform coordinates
        rmsd=None, tv1=None, rm=None, tv2=None,
        ## parse coordinates of selected atom types
        l_atoms=[],
        bool_None_if_unobs = False,
        bool_return_transformed_coordinates = False,
        ):

        rescount = 0
        coordinates1 = []
        coordinates2 = []
        d_coordinates1 = {}
        d_coordinates2 = {}
        lines1 = []
        lines2 = []
        l_RMSDs = []

        ## determine i range
        if l1 != 0 and l2 != 0:
            stop_not_expected
        if r1 != 0 and r2 != 0:
            stop_not_expected
        if len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-l2-r2 != len(d_ATOMseq[pdb2][chain2[0]]['res_nos'])-l1-r1:
            print len(d_ATOMseq[pdb1][chain1[0]]['res_nos']), l2, r2
            print len(d_ATOMseq[pdb2][chain2[0]]['res_nos']), l1, r1
            stop
        if l1 == 0 and l2 == 0:
            i_range = range(len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-abs(r2-r1))
        elif l2 == 0 and l1 > 0:
##            i_range = range(len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-abs(r2-r1)) ## wrong!
##            i_range = range(len(d_ATOMseq[pdb2][chain2[0]]['res_nos'])-abs(r2-r1))
            i_range = range(len(d_ATOMseq[pdb2][chain2[0]]['res_nos'])-l1-r1) 
        elif l1 == 0 and l2 > 0:
##            i_range = range(len(d_ATOMseq[pdb2][chain2[0]]['res_nos'])-abs(r2-r1)) ## wrong!
##            i_range = range(len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-abs(r2-r1))
            i_range = range(len(d_ATOMseq[pdb1][chain1[0]]['res_nos'])-l2-r2) 
        else:
            stop_not_expected

        if min(i_range)-1+l2 >= 0 and min(i_range)-1+l1 >= 0:
            print i_range
            stop
        if max(i_range)+1+l2 > len(d_ATOMseq[pdb1][chain1[0]]['res_nos']) and min(i_range)+1+l1 > len(d_ATOMseq[pdb2][chain2[0]]['res_nos'][i2]):
            print i_range
            stop

        ## loop over i range
        for i in i_range:
            i1 = i+l2
            i2 = i+l1
            res_no1 = d_ATOMseq[pdb1][chain1[0]]['res_nos'][i1]
            res_no2 = d_ATOMseq[pdb2][chain2[0]]['res_nos'][i2]
            iCode1 = d_ATOMseq[pdb1][chain1[0]]['iCodes'][i1]
            iCode2 = d_ATOMseq[pdb2][chain2[0]]['iCodes'][i2]
            if (
                d_ATOMseq[pdb1][chain1[0]]['records'][i1] == 'REMARK465'
                or
                d_ATOMseq[pdb2][chain2[0]]['records'][i2] == 'REMARK465'
                ):
                if bool_None_if_unobs == True:
                    coordinates1 += [None]
                    coordinates2 += [None]
                continue

            ## parse res_names associated with resID
            l_res_name1 = []
            for altloc1 in d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['altlocs'].keys():
                res_name1 = d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['altlocs'][altloc1]['res_name']
                l_res_name1 += [res_name1]
            l_res_name2 = []
            for altloc2 in d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['altlocs'].keys():
                res_name2 = d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['altlocs'][altloc2]['res_name']
                l_res_name2 += [res_name2]

            ## unknown if multiple different hetIDs per resID
            if len(set(l_res_name1)) > 1:
                res_name1 = 'UNK'
            if len(set(l_res_name2)) > 1:
                res_name2 = 'UNK'
                
            if res_name1 in self.d_modres.keys():
                res_name1 = self.d_modres[res_name1]
            if res_name2 in self.d_modres.keys():
                res_name2 = self.d_modres[res_name2]

            mutation = False
            if res_name1 != res_name2:
                mutation = True ## also e.g. LYS v LYZ (1xf6 v 1qgw)
                if res_name1 not in self.d_res1.keys() or res_name2 not in self.d_res1.keys():
                    print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2
##                    stop_temporary

            d_atoms1 = copy.deepcopy(d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'])
            d_atoms2 = copy.deepcopy(d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'])

            rescount += 1

            ##
            ## transformation of coordinates
            ##
            if rmsd:
##                for atom_name in d_atoms1.keys():
##                    if 'REMARK' in d_atoms1[atom_name].keys():
##                        continue
##                    coordinate1 = d_atoms1[atom_name]['coordinate']
##                    coordinates1 += [coordinate1]
                for atom_name in d_atoms2.keys():
                    if 'REMARK' in d_atoms2[atom_name].keys():
                        continue
                    coordinate2 = d_atoms2[atom_name]['coordinate']
                    coordinate2 = numpy.dot(coordinate2-tv1,rm)+tv2
                    d_atoms2[atom_name]['coordinate'] = coordinate2
##                    coordinates2 += [coordinate2]

            coordinates1_residue = []
            coordinates2_residue = []
            coordinates2_residue_nottransformed = []
            if rmsd:
                SS = []
            d_coordinates1[i1] = {}
            d_coordinates2[i2] = {}
            for atom_name in d_atoms1.keys():
                ## use only selected atoms
                if l_atoms != [] and atom_name not in l_atoms:
                    continue
                ## use backbone atoms only if mutated residue
                if mutation == True and atom_name not in ['N','CA','C','O']:
                    continue
                if atom_name not in d_atoms2.keys():
                    continue
                ## append coordinates to list of coordinates
                coordinate1 = d_atoms1[atom_name]['coordinate']
                coordinate2 = d_atoms2[atom_name]['coordinate']
                sqdiff = sum((coordinate2-coordinate1)**2)
##                ## only use secondary structure elements for structural alignment because of random coil super flexibility (e.g. 1ald v 1zal)
##                if (
##                    d_ATOMseq[pdb1][chain1[0]]['ss'][i1] != '' ## RANDOM COIL
##                    and
##                    d_ATOMseq[pdb2][chain2[0]]['ss'][i2] != '' ## RANDOM COIL
##                    ):
                ## only use alpha carbon atoms for structural alignment for speed purposes
                if l_atoms == [] or atom_name in l_atoms:
##                if atom_name == 'CA':
                    if rmsd:
                        coordinates1_residue += [coordinate1]
                        coordinates2_residue += [coordinate2]
                        coordinates2_residue_nottransformed += [d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name]['coordinate']]
                    else:
                        coordinates1 += [coordinate1]
                        coordinates2 += [coordinate2]
                d_coordinates1[i1][atom_name] = coordinate1
                d_coordinates2[i2][atom_name] = coordinate2
                if rmsd:
                    SS += [sqdiff]

            ##
            ## append lines
            ##
            if rmsd and SS != []:
                if d_ATOMseq[pdb2][chain2[0]]['ss'][i2] == '' and d_ATOMseq[pdb1][chain1[0]]['ss'][i1] == '':
##                    RMSD = 0.0
##                    RMSD = 100.0*rmsd
                    RMSD = math.sqrt(sum(SS)/len(SS))
                else:
                    RMSD = math.sqrt(sum(SS)/len(SS))
                l_RMSDs += [RMSD]

                ## append residue coordinates if residue rmsd below treshold
                if RMSD/rmsd < 9999999999.0:
                    coordinates1 += list(coordinates1_residue)
                    if bool_return_transformed_coordinates == False:
                        coordinates2 += list(coordinates2_residue_nottransformed)
                    else:
                        coordinates2 += list(coordinates2_residue)
                elif rmsd > 0.001:
                    print res_no1, res_no2, RMSD/rmsd
                    print RMSD, rmsd
                    stop

                d_line = {
                    pdb1:{
                        'res_name':res_name1,'chain':chain1,'res_no':res_no1,'d_atoms':d_atoms1,'iCode':iCode1,
                        },
                    pdb2:{
                        'res_name':res_name2,'chain':chain2,'res_no':res_no2,'d_atoms':d_atoms2,'iCode':iCode2,
                        },
                    }
                for pdb in d_line:
                    lines = []
                    d_atoms = d_line[pdb]['d_atoms']
                    res_name = d_line[pdb]['res_name']
                    chain = d_line[pdb]['chain']
                    res_no = d_line[pdb]['res_no']
                    iCode = d_line[pdb]['iCode']
##                    element = d_line[pdb]['element']
                    for atom_name in d_atoms.keys():
                        if 'REMARK' in d_atoms[atom_name].keys():
                            continue
                        coordinate = d_atoms[atom_name]['coordinate']
                        occupancy = RMSD
                        tempfactor = RMSD/rmsd
                        if tempfactor > 999.99:
                            tmpfactor = 999.99
                        if (
                            tempfactor > 10
                            and
                            (
                                (d_ATOMseq[pdb1][chain1[0]]['ss'][i1] == '' and d_ATOMseq[pdb2][chain2[0]]['ss'][i2] == '' and occupancy > 55) ## 1ex5
                                or
                                (d_ATOMseq[pdb1][chain1[0]]['ss'][i1] == 'HELIX' and d_ATOMseq[pdb2][chain2[0]]['ss'][i2] == 'HELIX' and occupancy > 17)
                                or
                                (d_ATOMseq[pdb1][chain1[0]]['ss'][i1] == 'SHEET' and d_ATOMseq[pdb2][chain2[0]]['ss'][i2] == 'SHEET' and occupancy > 10)
                                or
                                (d_ATOMseq[pdb1][chain1[0]]['ss'][i1] != '' and d_ATOMseq[pdb2][chain2[0]]['ss'][i2] != '' and occupancy > 17)
                                )
                            ):
                            print 'occupancy', occupancy
                            print 'tempfactor', tempfactor
                            print pdb1, chain1, res_no1, iCode1
                            print pdb2, chain2, res_no2, iCode2
                            print coordinate1
                            print coordinate2
                            print d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name]
                            print 'atom_name', atom_name
                            print d_ATOMseq[pdb1][chain1[0]]['ss'][i1]
                            print d_ATOMseq[pdb2][chain2[0]]['ss'][i2]
                            stop_occupancy
                        if occupancy > 999.99:
                            print 'occupancy', occupancy
                            stop
                        line = self.coordinates2ATOMline(res_name, chain[0], res_no, coordinate, iCode, occupancy, tempfactor, atom_name,)
                        lines += line
                    d_line[pdb]['lines'] = lines
                lines1 += d_line[pdb1]['lines']
                lines2 += d_line[pdb2]['lines']

        if len(coordinates1) == 0:
            print chain1, chain2
            stop ## e.g. 1ft8:E v 1koh:B

        return coordinates1, coordinates2, rescount, lines1, lines2, l_RMSDs, d_coordinates1, d_coordinates2


    def identify_interpdb_equivalent_chains_from_structure(
        self, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_coordinates, d_header,
        biomolecule1, biomolecule2,
        d_biomolecules1, d_biomolecules2,
        d_ATOMseq,
        ):

        import biounit, copy

##        print 'identifying structurally equivalent long peptides by structure between pdbs'

        bmchains1 = d_biomolecules1[biomolecule1]['chains']
        bmchains2 = d_biomolecules2[biomolecule2]['chains']

        n_chains = n_residues = n_coordinates = 0
        tv1 = 'N/A'
        rm = 'N/A'
        tv2 = 'N/A'

        ##
        ## find corresponding residue numbers
        ## this step after confirming rmsd is to be compared
        ## this step before transformation of chains and rmsd combinatorics
        ##
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            rep_chains2 = d_chains_interpdb_sequence_similar[rep_chain1].keys()
            ## e.g. 1nvv.pdb,1nvw.pdb
            for rep_chain2 in rep_chains2:
                l1SEQRES = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1']
                l2SEQRES = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2']
            
        l_chains1 = []
        l_chains2 = []
        chains1 = []
        chains2 = []

        ## loop over representative chains1
        for rep_chain1a in d_chains_interpdb_sequence_similar.keys():

            chains1seqid = [rep_chain1a]+d_chains_intrapdb_sequence_identical[pdb1][rep_chain1a]
            if len(set(chains1seqid) & set(chains1)) > 0:
                continue

            d_mutations = {}

            ##
            ## identify rep chains2 seq sim to rep chain1
            ##
            rep_chains2a = d_chains_interpdb_sequence_similar[rep_chain1a].keys()
            for rep_chain2 in rep_chains2a:
                n_mutations = len(d_chains_interpdb_sequence_similar[rep_chain1a][rep_chain2]['l_mutations'])
                if n_mutations not in d_mutations.keys():
                    d_mutations[n_mutations] = [] ## e.g. 1cb4.pdb,1e9p.pdb
                d_mutations[n_mutations] += [rep_chain2]
            l_mutations = d_mutations.keys()
            l_mutations.sort()


            ##
            ## identify other rep chains1 with seq sim to rep chains2 ## e.g. 1nvv.pdb,1nvw.pdb (Q;R;S == Q,R;S, Y64A mutant in one chain only)
            ##
            bmchains1seqsim = list(set([rep_chain1a]+d_chains_intrapdb_sequence_identical[pdb1][rep_chain1a]) & set(bmchains1))
            bmchains1seqsim.sort()
            for rep_chain1b in d_chains_interpdb_sequence_similar.keys():
                if rep_chain1a == rep_chain1b:
                    continue
                rep_chains2b = d_chains_interpdb_sequence_similar[rep_chain1b].keys()
                if len(set(rep_chains2a) & set(rep_chains2b)) > 0:
                    if set(rep_chains2a) ^ set(rep_chains2b) != set():
                        print pdb1,pdb2
                        print rep_chains2a, rep_chains2b
                        for chain in d_chains_interpdb_sequence_similar.keys():
                            print chain, d_chains_interpdb_sequence_similar[chain].keys()
                        if pdb1 != '1kj1' and pdb2 != '1bwu':
                            notexpected
                    bmchains1seqid = list(set([rep_chain1b]+d_chains_intrapdb_sequence_identical[pdb1][rep_chain1b]) & set(bmchains1))
                    bmchains1seqid.sort()
                    bmchains1seqsim += bmchains1seqid
                    
            ##
            ## identify chains seq id to the seq sim rep_chains2
            ##
            bmchains2seqsim = []
            ## loop over sorted number of mutations
            for n_mutation in l_mutations:
                for rep_chain2 in d_mutations[n_mutation]:
                    bmchains2seqid = list(set([rep_chain2]+d_chains_intrapdb_sequence_identical[pdb2][rep_chain2]) & set(bmchains2))
                    bmchains2seqid.sort()
                    bmchains2seqsim += bmchains2seqid
            
            ## identify chains identical to the representative chains
            chains1seqid = set(d_chains_intrapdb_sequence_identical[pdb1][rep_chain1a]+[rep_chain1a])
            ## exclude chains which are not in the biomolecule
            bmchains1seqid = list(set(chains1seqid) & set(bmchains1))
            ## sort chains
            bmchains1seqid.sort()
##            print rep_chain1a, bmchains1seqsim
##            print rep_chain1a, bmchains2seqsim
##            stop

            ## append chains to list of chains
            chains1 += bmchains1seqsim
            chains2 += bmchains2seqsim
            if bmchains1seqid != []:
                if bmchains1seqid not in l_chains1: ## e.g. 1sxa.pdb (vs 1e9o.pdb)
                    l_chains1 += [bmchains1seqsim]
            if bmchains2seqid != []:
                if bmchains2seqid not in l_chains2: ## e.g. 1sxa.pdb (vs 1e9o.pdb)
                    l_chains2 += [bmchains2seqsim]

##        print pdb1, biomolecule1, l_chains1
##        print pdb2, biomolecule2, l_chains2

        d_biomolecules = {
            pdb1:{'biomolecule':biomolecule1,'l_chains':l_chains1},
            pdb2:{'biomolecule':biomolecule2,'l_chains':l_chains2},
            }

        ##
        ## apply REMARK350 transformations
        ##
        ## identical number of chains after REMARK350 transformation
        ## e.g. A1,B1,A2,B2,A3,B3 == A,B,C,D,E,F of 1xnv.pdb,1xo6.pdb
        ## e.g. B,B,B,B == B,D,B,D of 1vwr.pdb,1vwi.pdb
        ## e.g. A,C = B,B of 1my3.pdb,1mxu.pdb

        ##
        ## 
        ##
        transformations = False
        for pdb in d_biomolecules.keys():
##                if 'REMARK350' not in d_header[pdb].keys():
##                    d_biomolecules[pdb]['tchains'] = d_biomolecules[pdb]['chains']
##                    continue
            biomolecule = d_biomolecules[pdb]['biomolecule']
            l_chains = d_biomolecules[pdb]['l_chains']
            l_tchains = []
            for chains in l_chains:
                tchains = []
                for chain in chains:
                    if 'REMARK350' not in d_header[pdb]:
                        tchains += [chain]
                    else:
                        matrix_nos = d_header[pdb]['REMARK350'][biomolecule]['chains'][chain]
                        for matrix_no in matrix_nos:
                            matrix = d_header[pdb]['REMARK350'][biomolecule]['matrices'][matrix_no]
                            d_coordinates, tchain = self.matrixtransformation(d_coordinates,pdb,chain,matrix,matrix_no)
                            tchains += [tchain]
                            if len(tchain) > 1:
                                if int(tchain[2:]) > 1:
                                    ## "transformations" True if matrix_no > 1 (trivial/identity matrix)
                                    ## a transformation hasnt taken place if second (or third...) matrix is also an identity matrix
                                    transformations = True
                l_tchains += [tchains]
            d_biomolecules[pdb]['l_tchains'] = l_tchains


        ##
        ## order chains by groups of sequence similar chains
        ##
        l_tchains1 = d_biomolecules[pdb1]['l_tchains']
        l_tchains2 = d_biomolecules[pdb2]['l_tchains']
        tchains1 = []
        for tchains in l_tchains1:
            tchains1 += list(tchains)
        tchains2 = []
        for tchains in l_tchains2:
            tchains2 += list(tchains)
        

        ##
        ## check if the expected correct combination of chains gives a low rmsd
        ##
        if len(tchains1) == len(tchains2):
            if self.verbose == True:
                print 'expected chain combination of REMARK350 transformations', tchains1, tchains2, l_tchains1, l_tchains2
            (
                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2
                ) = self.calculate_rmsd_for_multiple_chains(
                     tchains1,tchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,
                     d_header,
                     d_chains_interpdb_sequence_similar,
                     d_chains_intrapdb_sequence_identical,
                     d_ATOMseq,
                     )
            if self.verbose == True:
                print 'rmsd', rmsd
            rmsd_expected = rmsd
        else:
            rmsd = None
            if self.verbose == True:
                print 'rmsd cannot be determined for', tchains1, tchains2
        rmsd_expected = rmsd


        if n_chains != 0 and n_chains % 60 == 0:
            rmsd_max = self.d_rmsd_max['virus']
        else:
##            rmsd_max = self.max_rmsd
            rmsd_max = self.max_rmsd_wrong
            if len(tchains1) in self.d_rmsd_max.keys():
                rmsd_max = self.d_rmsd_max[len(tchains1)]
            else:
                rmsd_max = self.d_rmsd_max['multiple']
        ## correct transformation
        if rmsd != None and rmsd < rmsd_max:
            l_equivalent_chains = [tchains1,tchains2]
            return (
                l_equivalent_chains, rmsd, rmsd_expected,
                n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                )
        ## incorrect transformation and 1 chain v 1 chain
        elif len(tchains1) == 1 and len(tchains2) == 1:
            ## e.g. 1k53,1jml
            if rmsd > self.max_rmsd_wrong:
                self.incorrecttransformation(d_header,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,rmsd,)
            l_equivalent_chains = [tchains1,tchains2]
            return (
                l_equivalent_chains, rmsd, rmsd_expected,
                n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                )

        ## incorrect transformation and not 1 chain v 1 chain
        else:

##            ##
##            ## checks
##            ##
##            if len(l_tchains1) != len(l_tchains2):
##                print pdb1, pdb2
##                print l_tchains1, l_tchains2
##                notexpected
##            if len(tchains1) != len(tchains2):
##                print pdb1, pdb2
##                print tchains1, tchains2
##                notexpected
##            for i in range(len(l_tchains1)):
##                if len(l_tchains1[i]) != len(l_tchains2[i]):
##                    print pdb1, pdb2
##                    print l_tchains1, l_tchains2
##                    notexpected


            ##
            ## do different combinations of chain pairing than the expected (using remark350 transformations!)
            ##
            print 'different chain combinations of remark350 transformations',pdb1,pdb2,biomolecule1,biomolecule2
            if self.verbose == True:
                print 'l_tchains1', l_tchains1
                print 'l_tchains2', l_tchains2
            (
                l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                ) = self.shuffle_chains_and_calculate_rmsd(
                    d_coordinates,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
                    pdb1,pdb2,biomolecule1,biomolecule2,l_tchains1,l_tchains2,
                    d_ATOMseq,
                    )
            print rmsd, l_equivalent_chains
            if rmsd < rmsd_max:
                ## e.g. 2dtz,2hq5 for which PISA transformations does not exist and rmsd > 2.5
                ## e.g. 1u94,1u98
                return (
                    l_equivalent_chains, rmsd, rmsd_expected,
                    n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                    )

            if n_chains != None and n_chains % 60 == 0:
                print n_chains
                stop_not_expected

            ## parse PISA transformations
            d_transformations_PISA1,status1 = biounit.biounit().parse_pisa_multimers(pdb1, d_header[pdb1],)
            d_transformations_PISA2,status2 = biounit.biounit().parse_pisa_multimers(pdb2, d_header[pdb2],)

            ## PISA if not virus
            if n_chains == None or n_chains % 60 != 0:
                ##
                ## do PISA if 350 symmetry operations do not yield identical biounits
                ##
                l_chains1 = d_biomolecules[pdb1]['l_chains']
                l_chains2 = d_biomolecules[pdb2]['l_chains']
                chains1 = []
                for chains in l_chains1:
                    chains1 += chains
                chains2 = []
                for chains in l_chains2:
                    chains2 += chains

                if self.verbose == True:
                    print 'expected PISA transformations'
                for assembly1 in d_transformations_PISA1.keys():
    ##                ## no shared chains between PISA and REMARK350
    ##                if len(set(chains1)&set(d_transformations_PISA1[assembly1]['chains'].keys())) == 0:
    ##                    continue
                    d_coordinates, tchains1_PISA = self.apply_PISA_transformations(d_coordinates,pdb1,assembly1,d_transformations_PISA1, chains1)

                    ## PISA biounit built by chains other than the ones in the current biounit
                    if len(tchains1_PISA) == 0:
                        continue

    ##                ## require PISA and R350 multimers to have same size
    ##                if not (len(tchains1_PISA) == len(tchains1) or len(tchains1_PISA) == len(tchains2)):
    ##                    continue

                    for assembly2 in d_transformations_PISA2.keys():
                        print 'assemblies', assembly1, assembly2
    ##                    ## no shared chains between PISA and REMARK350
    ##                    if len(set(chains2)&set(d_transformations_PISA2[assembly2]['chains'].keys())) == 0:
    ##                        continue
                        d_coordinates, tchains2_PISA = self.apply_PISA_transformations(d_coordinates,pdb2,assembly2,d_transformations_PISA2, chains2)

                        ## PISA biounit built by chains other than the ones in the current biounit
                        if len(tchains2_PISA) == 0:
                            continue

    ##                    ## require PISA and R350 multimers to have same size
    ##                    if not (len(tchains2_PISA) == len(tchains2) or len(tchains2_PISA) != len(tchains2)):
    ##                        continue

                        ## different number of chains in the PISA multimers
                        if len(tchains1_PISA) != len(tchains2_PISA):
                            continue


                        ##
                        ## organize transformed chains by sequence identical chains (e.g. 1bgy 1be3)
                        ##
                        tchains1_PISA_seqsimorg = []
                        for chain1 in chains1:
                            rep_chain1 = self.chain2repchain(pdb1,chain1[0],d_chains_intrapdb_sequence_identical,)
                            for tchain1 in list(tchains1_PISA):
                                if tchain1[0] in d_chains_intrapdb_sequence_identical[pdb1][rep_chain1]+[rep_chain1]:
                                    tchains1_PISA_seqsimorg += [tchain1]
                                    tchains1_PISA.remove(tchain1)
                        tchains2_PISA_seqsimorg = []
                        for chain2 in chains2:
                            rep_chain2 = self.chain2repchain(pdb2,chain2[0],d_chains_intrapdb_sequence_identical,)
                            for tchain2 in list(tchains2_PISA):
                                if tchain2[0] in d_chains_intrapdb_sequence_identical[pdb2][rep_chain2]+[rep_chain2]:
                                    tchains2_PISA_seqsimorg += [tchain2]
                                    tchains2_PISA.remove(tchain2)
                        if len(tchains1_PISA) > 0:
                            stop1
                        if len(tchains2_PISA) > 0:
                            stop2
                        tchains1_PISA = tchains1_PISA_seqsimorg
                        tchains2_PISA = tchains2_PISA_seqsimorg


                        print 'PISA assemblies', assembly1, assembly2
                        print 'PISA', pdb1, pdb2
                        print tchains1_PISA
                        print tchains2_PISA
                        print 'expected chain combination of PISA transformations'
                        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                            tchains1_PISA,tchains2_PISA,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                            )
                        print 'PISA assemblies', assembly1, assembly2, rmsd

                        ##
                        ## check if the expected correct combination of chains gives a low rmsd
                        ##
                        if rmsd < rmsd_max: ## e.g. 1hqy,1ht2
                            l_equivalent_chains = [tchains1_PISA,tchains2_PISA]

                            if len(tchains1_PISA) == 1 and len(tchains2_PISA) == 1:
                                fn = 'pisa_suggested_combination_one_chain.txt' ## prev pisa1.txt
                            if len(tchains1_PISA) == 2 and len(tchains2_PISA) == 2:
                                fn = 'pisa_suggested_combination_two_chains.txt' ## prev pisa1.txt
                            else:
                                fn = 'pisa_suggested_combination_several_chains.txt' ## prev pisa1.txt

                            if os.path.isfile(fn):
                                fd = open(fn,'r')
                                lines = fd.readlines()
                                fd.close()
                            else:
                                lines = []

                            lines2 = []
                            for line in lines:
                                if [pdb1, pdb2, str(assembly1), str(assembly2),] != line.split()[:4]:
                                    lines2 += [line]
                            lines2 += ['%4s %4s %i %i %i %i %i %i %s\n' %(pdb1, pdb2, biomolecule1, biomolecule2, assembly1, assembly2, len(tchains1_PISA), len(tchains2_PISA), rmsd,)]

                            fd = open(fn,'w')
                            fd.writelines(lines2)
                            fd.close()

                            return (
                                l_equivalent_chains, rmsd, rmsd_expected,
                                n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                                )


                ##
                ## PISA v REMARK350 (e.g. 1fvk,1a2l;1a8m,4tsv)
                ##
                for assembly2 in d_transformations_PISA2.keys():
                    print 'assembly2', assembly2
                    d_coordinates, tchains2_PISA = self.apply_PISA_transformations(
                        d_coordinates, pdb2, assembly2, d_transformations_PISA2, chains2,
                        )

                    ## PISA biounit built by chains other than the ones in the current biounit
                    if len(tchains2_PISA) == 0:
                        continue

                    if len(tchains2_PISA) != len(tchains1):
                        continue
                    rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                        tchains1,tchains2_PISA,d_coordinates,pdb1,pdb2,
                        biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                        )
                    if rmsd < rmsd_max:
                        ## 1a8m,4tsv
                        l_equivalent_chains = [tchains1,tchains2_PISA]

                        fn = 'pisa_comb_w_remark350.txt' ## prev pisa3.txt

                        if os.path.isfile(fn):
                            fd = open(fn,'r')
                            lines = fd.readlines()
                            fd.close()
                        else:
                            lines = []

                        lines2 = []
                        for line in lines:
                            if [pdb1, pdb2, str(assembly1), str(assembly2),] != line.split()[:4]:
                                lines2 += [line]
                        lines2 += ['%s %s %s %s %s %s %s %s %s\n' %(pdb1, pdb2, assembly1, assembly2, biomolecule1, biomolecule2, len(tchains1), len(tchains2_PISA), rmsd,)]

                        fd = open(fn,'w')
                        fd.writelines(lines2)
                        fd.close()

                        return (
                            l_equivalent_chains, rmsd, rmsd_expected,
                            n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                            )
                for assembly1 in d_transformations_PISA1.keys():
                    print 'assembly1', assembly1
                    d_coordinates, tchains1_PISA = self.apply_PISA_transformations(
                        d_coordinates, pdb1, assembly1, d_transformations_PISA1, chains1,
                        )

                    ## PISA biounit built by chains other than the ones in the current biounit
                    if len(tchains1_PISA) == 0:
                        continue

                    if len(tchains1_PISA) != len(tchains2):
                        print tchains1_PISA, tchains2
                        continue
                    rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                        tchains1_PISA,tchains2,d_coordinates,pdb1,pdb2,
                        biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                        )
                    if rmsd < rmsd_max:
                        ## 2oz0,1fcb
                        l_equivalent_chains = [tchains1_PISA,tchains2]

                        fn = 'pisa_comb_w_remark350.txt' ## prev pisa3.txt

                        if os.path.isfile(fn):
                            fd = open(fn,'r')
                            lines = fd.readlines()
                            fd.close()
                        else:
                            lines = []

                        lines2 = []
                        for line in lines:
                            if [pdb1, pdb2, str(assembly1), str(assembly2),] != line.split()[:4]:
                                lines2 += [line]
                        lines2 += ['%s %s %s %s %s %s %s %s %s\n' %(pdb1, pdb2, assembly1, assembly2, biomolecule1, biomolecule2, len(tchains1_PISA), len(tchains2), rmsd,)]

                        fd = open(fn,'w')
                        fd.writelines(lines2)
                        fd.close()

                        return (
                            l_equivalent_chains, rmsd, rmsd_expected,
                            n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                            )


            if n_chains == None or n_chains % 60 != 0:
                ##
                ## do different combinations of chain pairing than the expected (using pisa transformations!)
                ##
                print 'different chain combinations of PISA transformations',pdb1,pdb2,biomolecule1,biomolecule2
                for assembly1 in d_transformations_PISA1.keys():
                    d_coordinates, tchains1_PISA = self.apply_PISA_transformations(d_coordinates,pdb1,assembly1,d_transformations_PISA1, chains1)
                    if len(tchains1_PISA) == 0:
                        continue
                    for assembly2 in d_transformations_PISA2.keys():
                        d_coordinates, tchains2_PISA = self.apply_PISA_transformations(d_coordinates,pdb2,assembly2,d_transformations_PISA2, chains2)
                        if len(tchains2_PISA) == 0:
                            continue
                        if len(tchains1_PISA) != len(tchains2_PISA):
                            continue
                        print 'PISA', pdb1, pdb2, tchains1_PISA, tchains2_PISA, chains1, chains2

                        d_l_tchains_PISA = {
                          pdb1:{'tchains':tchains1_PISA},
                            pdb2:{'tchains':tchains2_PISA},
                          }
                        for pdb in d_biomolecules.keys():
                            l_tchains_PISA = []
                            for l_chains in d_biomolecules[pdb]['l_chains']:
                                tchains_PISA = []
                                for chains in l_chains:
                                    for chain_PISA in d_l_tchains_PISA[pdb]['tchains']:
                                        if chain_PISA[0] in chains:
                                            tchains_PISA += [chain_PISA]
                                l_tchains_PISA += [tchains_PISA]
                            d_l_tchains_PISA[pdb]['l_tchains'] = l_tchains_PISA
                        l_tchains1_PISA = d_l_tchains_PISA[pdb1]['l_tchains']
                        l_tchains2_PISA = d_l_tchains_PISA[pdb2]['l_tchains']

                        tmp_chains1 = []
                        for tchains in l_tchains1_PISA:
                            tmp_chains1 += list(tchains)
                        tmp_chains2 = []
                        for tchains in l_tchains2_PISA:
                            tmp_chains2 += list(tchains)
                        ## continue if PISA multimers of different size
                        if len(tmp_chains1) != len(tmp_chains2):
                            continue
                        print tchains1
                        print tchains2
                        print tmp_chains1
                        print tmp_chains2
                        print n_chains

                        (
                            l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                            ) = self.shuffle_chains_and_calculate_rmsd(
                                d_coordinates,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
                                pdb1,pdb2,biomolecule1,biomolecule2,l_tchains1_PISA,l_tchains2_PISA,
                                d_ATOMseq,
                                )
                        if rmsd < rmsd_max: ## e.g. 1ryz,1t0u

                            if len(tmp_chains1) == 1 and len(tmp_chains2) == 1:
                                fn = 'pisa_diff_combinations_single_chain.txt' ## prev pisa2
                            else:
                                fn = 'pisa_diff_combinations_multiple_chains.txt' ## prev pisa2

                            if os.path.isfile(fn):
                                fd = open(fn,'r')
                                lines = fd.readlines()
                                fd.close()
                            else:
                                lines = []

                            lines2 = []
                            for line in lines:
                                if [pdb1, pdb2, str(assembly1), str(assembly2),] != line.split()[:4]:
                                    lines2 += [line]
                            lines2 += ['%s %s %s %s %s\n' %(pdb1, pdb2, assembly1, assembly2, rmsd,)]

                            fd = open(fn,'w')
                            fd.writelines(lines2)
                            fd.close()

                            return (
                                l_equivalent_chains, rmsd, rmsd_expected,
                                n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                                )

                ##
                ## do 290 symmetry operations if PISA does not yield identical biounits
                ##
                if rmsd > rmsd_max:

                    ##
                    ## REMARK 290 transformations, all combinations
                    ##
                    print 'apply REMARK290 transformations'
                    (
                        l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                        ) = self.remark290combinations(
                            d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,d_header,rmsd,
                            d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
                            d_ATOMseq,
                            )
                    if rmsd < rmsd_max: ## e.g. 1b8e,1bsy; 2lve,5lve
                        return (
                            l_equivalent_chains, rmsd, rmsd_expected,
                            n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                            )


        if rmsd == None or rmsd > rmsd_max: ## e.g. 2hqe, 2hqx
            l_tchains1 = d_biomolecules[pdb1]['l_tchains']
            l_tchains2 = d_biomolecules[pdb2]['l_tchains']
            tchains1 = []
            for tchains in l_tchains1:
                tchains1 += list(tchains)
            tchains2 = []
            for tchains in l_tchains2:
                tchains2 += list(tchains)
            l_equivalent_chains = [tchains1,tchains2]
            ##
            ## REMARK350, PISA, REMARK290 does not yield a correct transformation
            ##
            if rmsd == None or rmsd > rmsd_max:
                self.incorrecttransformation(d_header,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,rmsd,)
                if 'SWAP' in d_header[pdb1]['TITLE'] or 'SWAP' in d_header[pdb2]['TITLE']:
                    'rmsd', rmsd
                    print pdb1, d_header[pdb1]['TITLE']
                    print pdb2, d_header[pdb2]['TITLE']
                    stop

        if rmsd_expected != None and rmsd > rmsd_max:
            if len(tchains1) != len(tchains2):
                print rmsd, rmsd_expected
                print tchains1, tchains2
                stop_example
            l_equivalent_chains = [tchains1,tchains2]
            rmsd = rmsd_expected

        ##
        ## return 
        ##
        return (
            l_equivalent_chains, rmsd, rmsd_expected,
            n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
            )


    def shuffle_chains_and_calculate_rmsd(
        self,d_coordinates,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
        pdb1,pdb2,biomolecule1,biomolecule2,l_tchains1,l_tchains2,
        d_ATOMseq,
        ):

        time1 = time.clock()

        rmsd = None
        prev_rmsd = 999.99

        fd = open('virus.txt','r')
        lines = fd.readlines()
        fd.close()
        bool_prev_shuffle = False
        for line in lines:
            if pdb1 == line.split()[0] and pdb2 == line.split()[1]:
                chains1 = eval(line[line.index('['):line.index(']')+1])
                chains2 = eval(line[line.rindex('['):line.rindex(']')+1])
                bool_prev_shuffle = True
                break


        if bool_prev_shuffle == False:

            chains2 = []
            chains1 = []
            ksort = []
            for i in range(len(l_tchains1)):
                tchains1 = list(l_tchains1[i])
                tchains2 = list(l_tchains2[i])
                ## different biounits (PISA)
                if len(tchains1) != len(tchains2):
                    print l_tchains1, l_tchains2
                    print tchains1, tchains2
                    rmsd = 'N/A'
                    return None, rmsd, None, None, None, None, None, None
                ## one combination
                if len(tchains2) == 1:
                    chains1 += tchains1
                    chains2 += tchains2
                    continue
                ## multiple combinations
                else:
                    ## multiple different seq sim chains
                    if len(l_tchains1) > 1 and len(tchains1) == 60:
                        rmsdchains1 = chains1+tchains1
                        rmsdchains2 = chains2+tchains2
                        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                            rmsdchains1,rmsdchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                            )
                        ## expected combination of seq sim chains
                        if rmsd < self.max_rmsd:
                            chains2 += tchains2
                    else:
                        jskip = []
                        ##
                        ## loop over first chain of tchains1 (necessary when e.g. 1u94 v 1u98)
                        ##
                        for i_chain1 in range(1,len(tchains1)+1):
                            tchains2 = list(l_tchains2[i]) ## reset tchains2
                            chains1_tmp = list(chains1)
                            chains2_tmp = list(chains2)
                            ##
                            ## loop over second chain of tchains1
                            ##
                            for j in range(1,len(tchains1)+1):
                                if j in jskip:
                                    continue

                                ##
                                ## virus
                                ##
                                if len(tchains1) % 60 == 0 and (j-1) % 60 == 0 and len(tchains1) >= 60:

                                    rmsdchains1 = chains1+tchains1[:j-1+60]
        ##                            print 'rmsdchains1', rmsdchains1

                                    ## e.g. chains A,B,C of 1w39.pdb, 2fz1.pdb
                                    if ksort != range(1,61) and ksort != []:
                                        rmsdchains2 = list(chains2_tmp)
                                        for k in ksort:
                                            rmsdchains2 += [tchains2[k-1]]

        ##                                print 'rmsdchains2', rmsdchains2
                                        (
                                            rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                                            ) = self.calculate_rmsd_for_multiple_chains(
                                                rmsdchains1,rmsdchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False,
                                                )
                                        print rmsd, 'previous sequence'
                                        print ksort
                                        if rmsd < self.d_rmsd_max['virus']:
                                            chains1_tmp = rmsdchains1
                                            chains2_tmp = rmsdchains2
                                            tchains2 = tchains2[60:]
                                            jskip = range(j,j+60)
                                            ksort = ksort
                                            continue
        ## maybe A-transformation == C-transformation != B-transformation... then list of ksorts instead of resetting...

                                    ## e.g. chain C of 1aq4.pdb, 2bq5.pdb
                                    if len(tchains1) > 60:
                                        rmsdchains2 = list(chains2_tmp)+tchains2[:60]

        ##                                print 'rmsdchains2', rmsdchains2
                                        (
                                            rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                                            ) = self.calculate_rmsd_for_multiple_chains(
                                                rmsdchains1,rmsdchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                                                )
                                        print rmsd, 'integer sequence'
                                        if rmsd < self.d_rmsd_max['virus']:
                                            chains1_tmp = rmsdchains1
                                            chains2_tmp = rmsdchains2
                                            tchains2 = tchains2[60:]
                                            jskip = range(j,j+60)
                                            ksort = [] ## reset ksort
                                            continue
##                                        if (
##                                            (pdb1 == '1w39' and pdb2 == '2fz1') or
##                                            (pdb1 == '1auy' and pdb2 == '1w39') or
##                                            (pdb1 in ['1aq3','1aq4','1zdi',] and pdb2 in ['2b2e','2b2g','2bny',])
##                                            ):
##                                            if pdb1 == '1w39' and pdb2 == '2fz1':
##                                                seq = [1, 2, 3, 4, 5, 35, 31, 32, 33, 34, 29, 30, 26, 27, 28, 18, 19, 20, 16, 17, 10, 6, 7, 8, 9, 51, 52, 53, 54, 55, 48, 49, 50, 46, 47, 14, 15, 11, 12, 13, 23, 24, 25, 21, 22, 44, 45, 41, 42, 43, 60, 56, 57, 58, 59, 36, 37, 38, 39, 40]
##                                            if pdb1 == '1auy' and pdb2 == '1w39':
##                                                seq = [1, 2, 3, 4, 5, 22, 23, 24, 25, 21, 38, 39, 40, 36, 37, 19, 20, 16, 17, 18, 44, 45, 41, 42, 43, 13, 14, 15, 11, 12, 7, 8, 9, 10, 6, 56, 57, 58, 59, 60, 48, 49, 50, 46, 47, 34, 35, 31, 32, 33, 26, 27, 28, 29, 30, 52, 53, 54, 55, 51]
##                                            if pdb1 in ['1aq3','1aq4','1zdi',] and pdb2 in ['2b2e','2b2g','2bny',]:
##                                                seq = [2, 3, 4, 5, 1, 7, 8, 9, 10, 6, 12, 13, 14, 15, 11, 17, 18, 19, 20, 16, 22, 23, 24, 25, 21, 27, 28, 29, 30, 26, 32, 33, 34, 35, 31, 37, 38, 39, 40, 36, 42, 43, 44, 45, 41, 47, 48, 49, 50, 46, 52, 53, 54, 55, 51, 57, 58, 59, 60, 56]
##                                            rmsdchains2 = list(chains2)
##                                            for k in seq:
##                                                rmsdchains2 += [tchains2[k-1]]
##                                            (
##                                                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
##                                                ) = self.calculate_rmsd_for_multiple_chains(
##                                                    rmsdchains1,rmsdchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
##                                                    )
##                                            print rmsd, 'manual sequence'
##                                            if rmsd < self.max_rmsd:
##                                                chains1_tmp = rmsdchains1
##                                                chains2_tmp = rmsdchains2
##                                                tchains2 = tchains2[60:]
##                                                jskip = range(j,j+60)
##        ##                                        ksort = [1, 2, 3, 4, 5, 35, 31, 32, 33, 34, 29, 30, 26, 27, 28, 18, 19, 20, 16, 17, 10, 6, 7, 8, 9, 51, 52, 53, 54, 55, 48, 49, 50, 46, 47, 14, 15, 11, 12, 13, 23, 24, 25, 21, 22, 44, 45, 41, 42, 43, 60, 56, 57, 58, 59, 36, 37, 38, 39, 40]
##                                                continue
                                                

                                    jskip = []
                                    ksort = []

                                    ## end of if virus
                                    
                                minrmsd = ['N/A','N/A']
                                ##
                                ## loop over chains2 (virus and nonvirus)
                                ##
                                bool_rmsd_below_treshold = False
                                for k in range(len(tchains1)-j+1):
                                    chain2 = tchains2[k]
                                    if i_chain1 == j:
                                        rmsdchains1 = chains1+[tchains1[i_chain1-1]]+tchains1[:min(j,i_chain1-1)]+tchains1[max(j-1,i_chain1-1):min(j-1,i_chain1)]
                                    elif i_chain1 < j:
                                        rmsdchains1 = chains1+[tchains1[i_chain1-1]]+tchains1[0:i_chain1-1]+tchains1[i_chain1:j]
                                    elif i_chain1 > j:
                                        rmsdchains1 = chains1+[tchains1[i_chain1-1]]+tchains1[0:j-1]+tchains1[i_chain1:max(j,i_chain1)]
                                    rmsdchains2 = list(chains2_tmp)+[chain2]
                                    if len(rmsdchains1) != len(rmsdchains2) or len(rmsdchains1) != len(set(rmsdchains1)):
                                        print i_chain1, j
                                        print chains2, chains2_tmp, [chain2], rmsdchains2
                                        print
                                        print rmsdchains1
                                        print rmsdchains2
                                        stop_not_expected
                                    ## e.g. 2frp.pdb, 2gp1.pdb (chain X1 == chain X1)
                                    if len(tchains1) > 60 and j == 1 and (k % 60) != 0:
                                        continue
                                    ## e.g. 2g33.pdb, 2g34.pdb
                                    if len(tchains1) > 60 and (j-1) % 60 != 0 and chains2_tmp[-1][0] != chain2[0]:
                                        continue

                                    ## VIRUS: only align the past 3 chains (for two reasons: 1) virus dilemma max rmsd, 2) faster)
                                    if len(rmsdchains1) > 3:
                                        rmsdchains1 = rmsdchains1[-3:]
                                        rmsdchains2 = rmsdchains2[-3:]

                                    (
                                        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                                        ) = self.calculate_rmsd_for_multiple_chains(
                                            rmsdchains1,rmsdchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                                            )
                                    if len(rmsdchains1) >= 1: ## e.g. 1j4z,2eu1
                                        print 'rmsd %s %s %1i %1i %4s %4s %5.1f %3i of %3i %.1f' %(
                                            pdb1, pdb2, biomolecule1, biomolecule2,
                                            tchains1[j-1], chain2, rmsd, j, len(tchains1), rmsd/prev_rmsd,
                                            )

                                    if len(rmsdchains1) in self.d_rmsd_max.keys():
                                        rmsd_max = self.d_rmsd_max[len(rmsdchains1)]
                                    else:
                                        rmsd_max = self.d_rmsd_max['multiple']

                                    ## break if rmsd *beneath* treshold
                                    if (
                                        rmsd <= rmsd_max
                                        and
                                        ## solves dilemma: 2zp9,2zp8 (small rmsd if incorrect chain) v 1hrd,1aup (large rmsd if correct chain)
                                        ## factor 9 (2bc5,3c62)
                                        rmsd < 9*prev_rmsd
                                        ):
                                        ## append last chain to sequence of chains
                                        chains1_tmp += [rmsdchains1[-1]]
                                        chains2_tmp += [chain2]
                                        tchains2.remove(chain2)
                                        if len(tchains1) >= 60  and len(ksort) < 60:
                                            if len(chain2) > 1:
                                                ksort += [int(chain2[2:])]
                                            else:
                                                ksort += [1]
                                        prev_rmsd = rmsd
                                        bool_rmsd_below_treshold = True
                                        ## break loop over tchains2
                                        break
                                    ## determine if lower rmsd
                                    if rmsd < minrmsd[1]:
                                        minrmsd = [chain2,rmsd]
                                ## append chain with the lowest rmsd (because no chain was appended)
                                ## if it's wrong, then rmsd will go above the allowed treshold during one of the next steps...
                                if bool_rmsd_below_treshold == False:
                                    if k != len(tchains1)-j:
                                        notexpected
                                    chain2 = minrmsd[0]
                                    chains2_tmp += [chain2]
                                    tchains2.remove(chain2)
                                    if len(tchains1) >= 60 and len(ksort) < 60:
                                        if len(chain2) > 1:
                                            ksort += [int(chain2[2:])]
                                        else:
                                            ksort += [1]
                                    if len(tchains1) >= 60:
                                        print 'minrmsd %s %s %4s %4s %5.1f %3s of %3s' %(pdb1, pdb2, tchains1[j-1], chain2, minrmsd[1], j, len(tchains1))

                                ## e.g. 1a8r,1a9c; 1lcu,2q31 (unfortunately also 1u94 v 1u98)
                                if rmsd > self.max_rmsd_wrong:
                                    break ## end of loop over chains1
                                    return None, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

##                                if rmsd < self.max_rmsd and j == len(tchains1) and len(tchains2) == 0:
##                                    break

                                ## end of loop over chains1

                                continue

                            ## end of loop over i_chain1

                            if rmsd > self.max_rmsd_wrong:
                                continue ## continue i_chain1 loop

                            if rmsd < self.d_rmsd_max['multiple'] and j == len(tchains1) and len(tchains2) == 0:
                                break ## break i_chain1 loop

                        ## end of else

                        if rmsd > self.d_rmsd_max['multiple']:
                            l_equivalent_chains = None
                            return None, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

                        chains2 = list(chains2_tmp)

                if rmsd: ## rmsd calculated?
                    ## break loop over l_tchains1
                    if rmsd > self.max_rmsd_wrong:
                        stop_example_2_pdbs
                        break
                                
##                chains1 += tchains1
                chains1 = list(chains1_tmp)

                ## end of loop over l_tchains

            if rmsd < self.d_rmsd_max['virus'] and len(chains1)%60 == 0: ## temporary!!!
                fd = open('virus.txt','a')
                fd.write('%s %s %s %s\n' %(pdb1,pdb2,chains1,chains2))
                fd.close()

            ## end of if bool_prev_shuffle == False

        ## calculate final rmsd of all spatially identical chains
        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
            chains1,chains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True,
            )

        l_equivalent_chains = [chains1, chains2,]

        time2 = time.clock()

        print time2, time1
        if time2-time1 > 180:
            print time2-time1
            print rmsd
            print chains1
            print chains2
            stop_temp_write_txt_file_to_avoid_doing_all_combinations_again
            stop_this_took_too_long_dont_repeat

        return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2


    def apply_PISA_transformations(self,d_coordinates,pdb,assembly,d_transformations_PISA, chains):

        t_chains = []

##        if len(set(chains)&set(d_transformations_PISA[assembly]['chains'].keys())) == 0:
##            print pdb
##            print chains
##            print d_transformations_PISA[assembly]['chains'].keys()
##            print set(chains)-set(d_transformations_PISA[assembly]['chains'].keys())
##            stop_3h45_3d7e_excludinghetero_REMARK350_should_be_a_subset_of_PISA___at_least_one_shared_chain

##        for chain in chains:
        for chain in d_transformations_PISA[assembly]['chains'].keys():
##        ## loop over all chains in case the chains constituting the PISA biomolecule are different from those constituting the REMARK350 biomolecule (there should be no difference!!!)
##        for chain in d_coordinates[pdb]['chains'].keys():
            if chain not in d_transformations_PISA[assembly]['chains'].keys(): ## e.g. 1j4z
                continue
            if chain not in chains:
                continue
##                tchains = []
##                return d_coordinates, tchains

            ## sort molecules before loop (e.g. 1bgy, 1be3)
            l_molecules = d_transformations_PISA[assembly]['chains'][chain].keys()
            l_molecules.sort()
            for molecule in l_molecules:
                r = d_transformations_PISA[assembly]['chains'][chain][molecule]['r']
                t = d_transformations_PISA[assembly]['chains'][chain][molecule]['t']
                matrix = [
                    [r[0][0],r[0][1],r[0][2],t[0],],
                    [r[1][0],r[1][1],r[1][2],t[1],],
                    [r[2][0],r[2][1],r[2][2],t[2],],
                    ]
                d_coordinates, tchain = self.matrixtransformation(d_coordinates,pdb,chain,matrix,molecule,prefix='pisa',)
                t_chains += [tchain]

        return d_coordinates, t_chains



    def remark290combinations(
        self,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,d_header,rmsd,
        d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,
        ):

        sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
        import combinatorics

        fd = open('remark290.txt','r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            if line.split()[0] == pdb1 and line.split()[1] == pdb2 and int(line.split()[2]) == biomolecule1 and int(line.split()[3]) == biomolecule2:
                operator_combination1 = eval(line[line.index('['):line.index(']')+1])
                operator_combination2 = eval(line[line.rindex('['):line.rindex(']')+1])
                tchains1 = []
                for chain1 in chains1:
                    for operator1 in operator_combination1:
                        matrix1 = d_header[pdb1]['REMARK290'][operator1]
                        d_coordinates, tchain1 = self.matrixtransformation(d_coordinates,pdb1,chain1,matrix1,operator1)
                        tchains1 += [tchain1,]
                tchains2 = []
                for chain2 in chains2:
                    for operator2 in operator_combination2:
                        matrix2 = d_header[pdb2]['REMARK290'][operator2]
                        d_coordinates, tchain2 = self.matrixtransformation(d_coordinates,pdb2,chain2,matrix2,operator2)
                        tchains2 += [tchain2,]
                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                    tchains1,tchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                    )
                if rmsd < self.d_rmsd_max['multiple']:
                    l_equivalent_chains = [tchains1,tchains2]
                    return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2
                else:
                    stop

##        chains1 = d_header[pdb1]['REMARK350'][biomolecule1]['chains'].keys()
##        chains2 = d_header[pdb2]['REMARK350'][biomolecule2]['chains'].keys()
        operators1 = d_header[pdb1]['REMARK290'].keys()
        operators2 = d_header[pdb2]['REMARK290'].keys()

        if 'REMARK350' not in d_header[pdb1].keys() and 'REMARK350' not in d_header[pdb2].keys():
            if len(operators1) == 1 or len(operators2) == 1:
                n_matrices1 = 1
                n_matrices2 = 1
            else:
                print 'operators', operators1, operators2
                print 'chains', chains1, chains2
        else:
            if 'REMARK350' not in d_header[pdb1].keys():
                if len(operators1) == 1:
                    n_matrices1 = 1
                else:
                    n_matrices1 = len(chains1)
            else:
                n_matrices1 = len(d_header[pdb1]['REMARK350'][biomolecule1]['matrices'].keys())
##                if len(operators1) == 1:
##                    n_matrices1 = 1
##                else:
##                    n_matrices1 = len(chains1)
            if 'REMARK350' not in d_header[pdb2].keys():
                if len(operators2) == 1:
                    n_matrices2 = 1
                else:
                    n_matrices2 = len(chains2)
            else:
                n_matrices2 = len(d_header[pdb2]['REMARK350'][biomolecule2]['matrices'].keys())
##                if len(operators2) == 1:
##                    n_matrices2 = 1
##                else:
##                    n_matrices2 = len(chains2)

        print 'matrices', n_matrices1, n_matrices2
        print 'operators', operators1, operators2
        if (
            n_matrices1*len(chains1) != n_matrices2*len(chains2)
            ):
            print chains1, n_matrices1, operators1
            print chains2, n_matrices2, operators2
            print 'different_sized_multimers'
            return None, rmsd, None, None, None, None, None, None

        if len(operators1) == 1 and n_matrices1 > 1:
            operator_combinations1 = n_matrices1*[operators1]
        else:
            operator_combinations1 = combinatorics.permutation_wo_rep(operators1, n_matrices1)
        print 'sort1'
        operator_combinations1.sort()

        print operators2, n_matrices2
        stop
        if len(operators2) == 1 and n_matrices2 > 1:
            operator_combinations2 = n_matrices2*[operators2]
        else:
            operator_combinations2 = combinatorics.permutation_wo_rep(operators2, n_matrices2)
        print 'sort2'
        operator_combinations2.sort()

        if len(chains1) != len(chains2):
            if (
                (len(chains1) > 2 or len(chains2) > 2) and
                len(chains1) % len(chains2) != 0 and len(chains2) % len(chains1) != 0 and
                len(chains1)*len(d_header[pdb1]['REMARK350'][biomolecule1]['matrices'].keys()) != len(chains2)*len(d_header[pdb2]['REMARK350'][biomolecule2]['matrices'].keys())
                ):
                print chains1, chains2
                print operator_combinations1
                print operator_combinations2
                stop1

##        l_translations = combinatorics.permutation_w_rep([0,1,-1,],3)
        for operator_combination1 in operator_combinations1:
            tchains1 = []
            ## chain loop before operator loop for comparsion of interpdb sequence identical chains
            for chain1 in chains1:
                if len(chain1) > 1:
                    print chain1, chains2
                    notexpected
                for operator1 in operator_combination1:
                    matrix1 = d_header[pdb1]['REMARK290'][operator1]
                    d_coordinates, tchain1 = self.matrixtransformation(d_coordinates,pdb1,chain1,matrix1,operator1)
                    tchains1 += [tchain1,]

            for operator_combination2 in operator_combinations2:
                ## chain loop before operator loop for comparsion of interpdb sequence identical chains
                tchains2 = []
                for chain2 in chains2:
                    if len(chain2) > 1:
                        print chain2, chains2
                        notexpected
                    for operator2 in operator_combination2:
                        matrix2 = d_header[pdb2]['REMARK290'][operator2]
                        d_coordinates, tchain2 = self.matrixtransformation(d_coordinates,pdb2,chain2,matrix2,operator2)
                        tchains2 += [tchain2,]

                if len(tchains1) != len(tchains2):
##                    print operator_combinations1
##                    print operator_combinations2
##                    if (
##                        set([pdb1,pdb2]) != set(['1qve','1rg0']) and
##                        set([pdb1,pdb2]) != set(['1m3s','1viv']) and
##                        set([pdb1,pdb2]) != set(['1r4c','1tij']) and
##                        set([pdb1,pdb2]) != set(['2qkt','2qku']) ## 350biounit = monomer AND dimer
##                        ):
##                        print chains1, tchains1
##                        print chains2, tchains2
##                        print pdb1, pdb2
##                        print d_header[pdb1]['REMARK350'][biomolecule1]['matrices'].keys()
##                        print d_header[pdb2]['REMARK350'][biomolecule2]['matrices'].keys()
##                        print self.cluster, 'cluster'
                        notexpected
##                    else:
##                        continue
                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                    tchains1,tchains2,d_coordinates,pdb1,pdb2,biomolecule1,biomolecule2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                    )
                print rmsd, pdb1, pdb2, biomolecule1, biomolecule2, operator_combination1, operator_combination2
                if n_chains not in self.d_rmsd_max.keys():
                    rmsd_max = self.d_rmsd_max['multiple']
                else:
                    rmsd_max = self.d_rmsd_max[n_chains]
                if rmsd < rmsd_max:
                    fd = open('remark290.txt','a')
                    fd.write('%s %s %s %s %s %s %s\n' %(pdb1,pdb2,biomolecule1,biomolecule2,operator_combination1,operator_combination2,rmsd))
                    fd.close()
                    l_equivalent_chains = [tchains1,tchains2]
                    return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

        l_equivalent_chains = None
        n_chains = 0
        n_residues = 0
        n_coordinates = 0
        tv1 = None
        rm = None
        tv2 = None
        
        return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

    def incorrecttransformation(self,d_header,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,rmsd,):

        if rmsd in [None,'N/A',]:
            rmsd = 99.9

        ## toomany or just a high rmsd (due to incorrect transformation or something else...)
        toomany = []
        if os.path.isfile('incorrecttransformation.txt'):
            fd = open('incorrecttransformation.txt','r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                toomany += [ [line.split()[0],line.split()[1],int(line.split()[2]),int(line.split()[3])] ]
        if not [pdb1,pdb2,biomolecule1,biomolecule2] in toomany and not [pdb2,pdb1,biomolecule2,biomolecule1] in toomany:
            fd = open('incorrecttransformation.txt','a')
            fd.write('%s %s %s %s %3i %3i  %5.1f %5s %10s %10s %5s %s %s\n' %(
                pdb1, pdb2, biomolecule1, biomolecule2,
                len(chains1), len(chains2),
                rmsd,
                d_header[pdb1]['CRYST1']==d_header[pdb2]['CRYST1'], d_header[pdb1]['CRYST1'], d_header[pdb2]['CRYST1'],
                len(chains1)==len(chains2), len(chains1), len(chains2),
                ))
            fd.close()

        return


    def matrixtransformation(self,d_coordinates,s_pdb,chain,matrix,matrix_no,prefix='',):

        '''apply transformation'''

        ## matrix does not cause transformation
        if matrix == self.nontransformationmatrix:
            return d_coordinates,chain

##        print 'applying REMARK350 transformation to chain %s of pdb %s' %(chain, s_pdb)

        rmatrix, tvector = self.transformationmatrix2rotationmatrix_and_translationvector(matrix)

        tchain = chain+'_'+prefix+str(matrix_no)

        d_coordinates[s_pdb]['chains'][tchain] = {}
        if 'residues' not in d_coordinates[s_pdb]['chains'][tchain].keys():
            d_coordinates[s_pdb]['chains'][tchain]['residues'] = {}
        for res_no in d_coordinates[s_pdb]['chains'][chain]['residues'].keys():
            if res_no not in d_coordinates[s_pdb]['chains'][tchain]['residues'].keys():
                d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no] = {}
            if 'd_iCodes' not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'] = {}
            if 'l_iCodes' not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['l_iCodes'] = []
            for iCode in d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['l_iCodes']:
                if 'REMARK' in d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no] = d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]
                    continue
                if iCode not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'].keys():
                    d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode] = {}
                if 'altlocs' not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'] = {}
                for altloc in d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                    if altloc not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                        d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc] = {}
                    d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name'] = d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']
                if 'atoms' not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                for atom_name in d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                    if 'REMARK' in d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name].keys():
                        d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]
                        continue
                    if atom_name not in d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                    coord = d_coordinates[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate']
                    tcoord = numpy.dot(rmatrix, coord) + tvector
                    d_coordinates[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate'] = tcoord

        return d_coordinates, tchain


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):

        translationvector = numpy.array(
            [
                float(transformationmatrix[0][3]),
                float(transformationmatrix[1][3]),
                float(transformationmatrix[2][3]),
                ]
            )

        rotationmatrix = numpy.array(
            [
                [
                    float(transformationmatrix[0][0]),
                    float(transformationmatrix[0][1]),
                    float(transformationmatrix[0][2]),
                    ],
                [
                    float(transformationmatrix[1][0]),
                    float(transformationmatrix[1][1]),
                    float(transformationmatrix[1][2]),
                    ],
                [
                    float(transformationmatrix[2][0]),
                    float(transformationmatrix[2][1]),
                    float(transformationmatrix[2][2]),
                    ],
                ]
            )

        return rotationmatrix, translationvector


    def identify_identical_chains_from_sequence_intra(
        self, d_header, pdb,
        ):

        d_chains_intrapdb_sequence_identical = {}

        chains = d_header[pdb]['SEQRES']['chains'].keys()
        waterchains = set(d_header[pdb]['REMARK525'])

        ## return the following dictionary structure
        ## intrapdb: {pdb:{repchain:[seqidchains]}}

##        print 'identify_identical_chains_from_sequence', pdb1, pdb2

        d_chains_intrapdb_sequence_identical = {}

        for i in range(len(chains)):
            chaini = chains[i]
            if chaini in waterchains:
                continue

            ## continue if chaini identical to previous chaink
            identical = False
            for chaink in d_chains_intrapdb_sequence_identical.keys():
                if chaini in d_chains_intrapdb_sequence_identical[chaink]:
                    identical = True
            if identical == True:
                continue

            ## initiate list for chaini
            if chaini not in d_chains_intrapdb_sequence_identical.keys():
                d_chains_intrapdb_sequence_identical[chaini] = []

            d_chains_intrapdb_sequence_identical[chaini] = []

            seqi = d_header[pdb]['SEQRES']['chains'][chaini]['seq']

            for j in range(i+1,len(chains)):
                chainj = chains[j]
                if chainj in waterchains:
                    continue

                seqj = d_header[pdb]['SEQRES']['chains'][chainj]['seq']

                if seqi == seqj:
                    d_chains_intrapdb_sequence_identical[chaini] += [chainj]
##                    d_chains_intrapdb_sequence_identical[chainj] = chaini
                ## if last chain in loop and not identical to anything then append to list of representative chains
                elif i == len(chains)-2 and j == len(chains)-1:
                    identical = False
                    for chaink in d_chains_intrapdb_sequence_identical.keys():
                        if chainj in d_chains_intrapdb_sequence_identical[chaink]:
                            identical = True
                    if identical == False:
                        d_chains_intrapdb_sequence_identical[chainj] = []

        return d_chains_intrapdb_sequence_identical


    def identify_similar_chains_from_sequence_inter(
        self, d_header, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical,
        bmchains1, bmchains2,
        d_biomolecules1, d_biomolecules2,
        ):

##        print 'identifying similar chains'

        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        import sequence_alignment

        ## return the following dictionary structure
        ## interpdb: {repchain:{seqsimchain:{l1,l2}}}

        ## identify repchains
        repchains1 = []
        for repchain1 in d_chains_intrapdb_sequence_identical[pdb1].keys():
            chains1seqid = [repchain1]+d_chains_intrapdb_sequence_identical[pdb1][repchain1]
            for bmchain1 in bmchains1:
                if bmchain1 in chains1seqid:
                    repchains1 += [repchain1]
                    break
        repchains2 = []
        for repchain2 in d_chains_intrapdb_sequence_identical[pdb2].keys():
            chains2seqid = [repchain2]+d_chains_intrapdb_sequence_identical[pdb2][repchain2]
            for bmchain2 in bmchains2:
                if bmchain2 in chains2seqid:
                    repchains2 += [repchain2]
                    break

        ## set d_chains_interpdb_sequence_similar
        d_chains_interpdb_sequence_similar = {}

        ## only do sequential alignment for representative chains to save time
        for chain1 in repchains1:

            if d_header[pdb1]['SEQRES']['chains'][chain1]['type'] != 'peptide':
                continue

            seq1 = d_header[pdb1]['SEQRES']['chains'][chain1]['seq']

            if len(seq1) < self.min_len_chain:
                continue

            ## only do sequential alignment for representative chains to save time
            for chain2 in repchains2:

                if d_header[pdb2]['SEQRES']['chains'][chain2]['type'] != 'peptide':
                    continue

                seq2 = d_header[pdb2]['SEQRES']['chains'][chain2]['seq']

                if len(seq2) < self.min_len_chain:
                    continue

                if abs(len(seq1)-len(seq2)) > self.max_len_chain_difference:
                    continue

                ## fast sequence comparison
                if len(seq1) == len(seq2):

                    l1 = 0
                    l2 = 0
                    r1 = 0
                    r2 = 0
                    s1 = seq1
                    s2 = seq2

                    ## sequence identical
                    if seq1 == seq2:

                        n_chainmutations = 0
                        l_chainmutations = []

                    ## point mutation(s)
                    else:

                        s1 = seq1
                        s2 = seq2
                        n_chainmutations, l_chainmutations = self.point_mutations(seq1, seq2, l1, l2,)

                else:

                    if seq1 in seq2:

                        l1 = seq2.index(seq1)
                        l2 = 0
                        r1 = (seq2.rindex(seq1)-l1)+(len(seq2)-len(seq1)-l1)
                        r2 = 0
                        s1 = seq1[l1:len(seq1)-r1]
                        s2 = seq2
                        n_chainmutations = 0
                        l_chainmutations = []


                    elif seq2 in seq1:

                        l1 = 0
                        l2 = seq1.index(seq2)
                        r1 = 0
                        r2 = (seq1.rindex(seq2)-l2)+(len(seq1)-len(seq2)-l2)
                        s1 = seq1
                        s2 = seq2[l2:len(seq2)-r2]
                        n_chainmutations = 0
                        l_chainmutations = []

                    else:

                        ## 1st slow sequence comparison (SEQRESseq)

                        if self.verbose == True:
                            print pdb1, pdb2, chain1, chain2, 'begin seq aln of chains of len %s and %s' %(len(seq1),len(seq2))
                        instance = sequence_alignment.NW(seq1,seq2)
                        s1,s2 = instance.Align(verbose=False)[:2]
                        if self.verbose == True:
                            print 'end seq aln 1'

                        l1 = len(s1)-len(s1.lstrip('-'))
                        l2 = len(s2)-len(s2.lstrip('-'))
                        r1 = len(s1)-len(s1.rstrip('-'))
                        r2 = len(s2)-len(s2.rstrip('-'))

        ## change .1 to variable...
                        if l1 > self.max_len_chain_difference or l1/float(len(s1)) > .1:
                            continue
                        if r1 > self.max_len_chain_difference or r1/float(len(s1)) > .1:
                            continue
                        if l2 > self.max_len_chain_difference or l2/float(len(s2)) > .1:
                            continue
                        if r2 > self.max_len_chain_difference or r2/float(len(s2)) > .1:
                            continue

                        l = max(l1,l2)
                        r = max(r1,r2)

                        s1 = s1[l:len(s1)-r]
                        s2 = s2[l:len(s2)-r]

                        ## continue if insertions/deletions
                        if '-' in s1 or '-' in s2:
                            continue

                        n_chainmutations, l_chainmutations = self.point_mutations(s1, s2, l1, l2,)

                if n_chainmutations <= self.max_mutations:
                    if chain1 not in d_chains_interpdb_sequence_similar.keys():
                        d_chains_interpdb_sequence_similar[chain1] = {}
                    if chain2 in d_chains_interpdb_sequence_similar[chain1].keys():
                        notexpected
                    d_chains_interpdb_sequence_similar[chain1][chain2] = {
                        'l1':l1,'l2':l2,
                        's1':s1,'s2':s2,
                        'l_mutations':l_chainmutations,
                        'r1':r1,'r2':r2,
                        }

##        print d_chains_interpdb_sequence_similar
    
        return d_chains_interpdb_sequence_similar


    def point_mutations(self, s1, s2, l1, l2,):

        n_chainmutations = 0
        l_chainmutations = []
        for res in range(len(s1)):
            res1 = s1[res]
            res2 = s2[res]
            if res1 != res2:
                n_chainmutations += 1
                l_chainmutations += [[res+l2,res+l1,res1,res2]]
            if n_chainmutations > self.max_mutations:
                break

        return n_chainmutations, l_chainmutations

        


    def res_name2res_symbol(self, res_name):

        ## this function should take into account nucleotides
        ## unfortunately there is a conflict between:
        ## G: glycine and guanosine
        ## A: alanine and adenosine
        ## T: threonine and thymidine
        ## C: cysteine and cytidine
        ## N: aspargine and unknown nucleotide residue

        if res_name in self.d_res1.keys():
            symbol = self.d_res[res_name]
        elif res_name in self.l_nucleotides:
            symbol = res_name[-1]
        else:
            symbol = 'X'

        return symbol
    

    def ATOM2seq(self, d_coordinates, chain, d_header, res_no_max = None, stop_error = True):

        stop_this_still_in_use

        d_res_nos = {}
        seq = ''
        ATOMrespos = 0
        res_nos = d_coordinates['chains'][chain]['residues'].keys()
        res_nos.sort()
        for i in range(len(res_nos)):
            res_no = res_nos[i]
            if res_no_max != None and res_no >= res_no_max:
                continue

            for i in range(len(d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'])):
                iCode = d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'][i]
                try:
                    res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']
                except:
                    print d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
                    stop

                if res_name == 'HOH':
                    continue
                
                ##
                ## do not append to sequence if hetID is not a MODRES
                ## e.g. TRP in 1utv.pdb
                ##
                if chain in d_header['HET'].keys():
                    if res_no in d_header['HET'][chain].keys():
                        if iCode in d_header['HET'][chain][res_no].keys():
                            if res_name not in d_header['HET'][chain][res_no][iCode]:
                                if stop_error == True:
                                    print chain, res_no, iCode, res_name, d_header['HET'][chain][res_no][iCode], d_header['MODRES'][chain][res_no][iCode]
                                    notexpected
                            else:
                                if not chain in d_header['MODRES'].keys():
                                    continue
                                elif res_no not in d_header['MODRES'][chain].keys():
                                    continue
                                elif iCode not in d_header['MODRES'][chain][res_no].keys():
                                    continue
                                else:
                                    pass
                                
                d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':iCode}
                seq += self.res_name2res_symbol(res_name)
                ATOMrespos += 1

        return seq, d_res_nos


    def coordinates2ATOMline(self, res_name, chain, res_no, coordinate, iCode, occupancy, tempfactor, atom_name):

        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        atom_no = 1
        altloc = ''
        charge = ''
        element = ''
        if len(atom_name) < 4:
            atom_name = ' %3s' %(atom_name.ljust(3)) ## not sure if this is entirely correct
        line = [
            '%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n'
            %('ATOM'.ljust(6), atom_no, atom_name, altloc, res_name.ljust(3), chain, res_no, iCode, x, y, z, occupancy, tempfactor, element.rjust(2), charge.rjust(2))
            ]

        return line


    def parse_atom_no_range(self, d_conect, record, atom_no):
        
        if not record in d_conect.keys():
            d_conect[record] = [[atom_no,atom_no]]
        elif d_conect[record][-1][1] == atom_no-1:
            d_conect[record][-1][1] = atom_no
        else:
            d_conect[record] += [[atom_no,atom_no]]

        return d_conect


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):

        translationvector = numpy.array(
            [
                float(transformationmatrix[0][3]),
                float(transformationmatrix[1][3]),
                float(transformationmatrix[2][3]),
                ]
            )

        rotationmatrix = numpy.array(
            [
                [
                    float(transformationmatrix[0][0]),
                    float(transformationmatrix[0][1]),
                    float(transformationmatrix[0][2]),
                    ],
                [
                    float(transformationmatrix[1][0]),
                    float(transformationmatrix[1][1]),
                    float(transformationmatrix[1][2]),
                    ],
                [
                    float(transformationmatrix[2][0]),
                    float(transformationmatrix[2][1]),
                    float(transformationmatrix[2][2]),
                    ],
                ]
            )

        return rotationmatrix, translationvector


    def rewind_dictionary(self, root, path, d_molecules, d_adjacency_forward):

        ## forward in dictionary to previous branching point
        d_m_fwd = d_molecules[root]['bonds']
        for i in range(1,len(path)):
            m1 = path[i-1]
            m2 = path[i]
            bond = d_adjacency_forward[m1][m2]
            d_m_fwd = d_m_fwd[bond]['bonds']
        
        return d_m_fwd


    def identify_iCode_sequence(self, d_coordinates, chain, res_no, iCode, res_name, d_header):
        
        l_iCodes = list(d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'])
        d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
        for iCode_prev in l_iCodes:
            index_alphabet = self.s_alphabet.index(iCode)
            if (
                ## REMARK465
                'REMARK' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev].keys() or
                ## REMARK470
                {'REMARK':True} in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev]['atoms'].values()
                ):
                ATOMseq,d_res_nos = self.ATOM2seq(d_coordinates, chain, d_header, res_no_max = res_no)
                res_symbol = self.res_name2res_symbol(res_name)
                ## 1) REMARK residues before ATOM residues (N-terminal)
                ## e.g. 1jqz.pdb
                if (
                    len(ATOMseq) == 0
                    ):
                    None
                ## 2) REMARK465 residues after ATOM residues
                ## e.g. 1nuo.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES']['chains'][chain]['seq'][len(ATOMseq)] and
                    ATOMseq == d_header['SEQRES']['chains'][chain]['seq'][:len(ATOMseq)]
                    ):
                    l_iCodes = [iCode]+d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'][:-1]
                    d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 3) REMARK470 residues after ATOM residues
                ## e.g. 2lve.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES']['chains'][chain]['seq'][len(ATOMseq)+index_alphabet] and
                    ATOMseq == d_header['SEQRES']['chains'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' in l_iCodes
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 4) REMARK470 residues after ATOM residues
                ## e.g. 2j5q (chain B, res_no 54, iCode C)
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES']['chains'][chain]['seq'][len(ATOMseq)+index_alphabet-1] and
                    ATOMseq == d_header['SEQRES']['chains'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' not in l_iCodes 
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 5) REMARK465 residues before ATOM residues
                ## e.g. 1uij.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES']['chains'][chain]['seq'][len(ATOMseq)+len(l_iCodes)] and
                    ATOMseq == d_header['SEQRES']['chains'][chain]['seq'][:len(ATOMseq)]
                    ):
                    None
                else: ## e.g. 1fne
                    print '*****'
                    print iCode, iCode_prev, res_symbol, index_alphabet, l_iCodes
                    print d_header['SEQRES']['chains'][chain]['seq'][len(ATOMseq)+index_alphabet]
                    print
                    print len(ATOMseq) > 0
                    print res_symbol == d_header['SEQRES']['chains'][chain]['seq'][len(ATOMseq)+index_alphabet-1]
                    print ATOMseq == d_header['SEQRES']['chains'][chain]['seq'][:len(ATOMseq)]
                    print self.s_alphabet[index_alphabet-1] in l_iCodes
                    print 
                    print 1, ATOMseq, self.res_name2res_symbol(res_name)
                    print 2, d_header['SEQRES']['chains'][chain]['seq']
                    print chain, res_no, iCode, iCode_prev, res_name
                    expected
                break

        return d_coordinates


    def __init__(self):

        import sys
        sys.path.append('/home/people/tc/svn/tc_sandbox/pdb')
        import quakes_init
        import smallmolecules

        d = smallmolecules.main()

        ## connectivity *and* identity of solutes *and* ions not checked
        self.l_solutes = d['solutes']
        self.d_ions = d['ions']

        self.d_saccharides = d['saccharides']

        self.d_spacegroups, self.d_crystalsystems = quakes_init.main()

        self.d_res1 = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            'MSE':'M','UNK':'X','ASX':'X','GLX':'X',
            }

        ## HETATM res_names for which coordinates are parsed
        self.d_modres = {
            'MSE':'MET', ## selenomethionine
##            ## phosphorylation
##            'TPO':'THR',
##            'SEP':'SER',
##            'PHD':'ASP',
##            'PTR':'TYR',
            }

        self.d_res3 = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

        self.d_res_names = {
            'A':'ALANINE',
            'C':'CYSTEINE',
            'D':'ASPARTIC ACID', ## 'ASPARTATE',
            'E':'GLUTAMIC ACID', ## 'GLUTAMATE',
            'F':'PHENYLALANINE',
            'G':'GLYCINE',
            'H':'HISTIDINE',
            'I':'ISOLEUCINE',
            'K':'LYSINE',
            'L':'LEUCINE',
            'M':'METHIONINE',
            'N':'ASPARAGINE',
            'P':'PROLINE',
            'Q':'GLUTAMINE',
            'R':'ARGININE',
            'S':'SERINE',
            'T':'THREONINE',
            'V':'VALINE',
            'W':'TRYPTOPHAN',
            'Y':'TYROSINE',
            }

        self.time_status_update = 1

        self.nontransformationmatrix = [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]

        self.s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        self.l_comments_true = [
            'ENGINEERED','ENGINEERED MUTATION','MUTATION','SUBSTITUTION','VARIANT',
##            'CONFLICT SEE REMARK 9','CONFLICT',
            'RANDOM MUTAGENESIS',
            ## not a mutation, but a posttranslational modification
            'MODIFIED RESIDUE','HYDROXYLATION',
            ## natural mutant
            'POLYMORPHISM',
            ]
        self.l_comments_false = [
            'CLONING ARTIFACT','EXPRESSION TAG',
            'INITIATING METHIONINE','INITIATING MET',#'INSERTION AT N-TERMINUS',
            ]

        fd = open('d_dihedrals.txt','r')
        s = fd.read()
        fd.close()
        self.d_dihedrals_atoms = eval(s)

## delete this (not used here)
##############        self.l_expdta = [
##############            'X-RAY',
##############            'ELECTRON DIFFRACTION','NEUTRON DIFFRACTION',
##############            'NMR',
##############            'INFRARED SPECTROSCOPY',
##############            'CRYO-ELECTRON MICROSCOPY',
##############            'ELECTRON TOMOGRAPHY', ## e.g. 1o1a.pdb
##############            'SOLUTION SCATTERING', ## e.g. 1e07.pdb
##############            'FLUORESCENCE TRANSFER', ## e.g. 1rmn.pdb
##############            ]

        self.d_rmsd_max = { ## 1r7r,3cf3; 1nlf,1olo; 2fsy,2ft1 ## try other combinations if above this rmsd (must be below 5.0 for 1hrd,1aup)
            ## key = count of chains
            ## value = max rmsd
            1:4.69, ## 4.25 3dfp,1ald(3.63);2f3d,1rdz(4.69)
            2:6.19, ## 1hrd,1aup;1qpp,1qpx;2c1v,2c1u(5.56);1qku,1qkt(5.82);2qvf,1q1p(6.19)
            3:6.21, ## 1hrd,1aup;3eiz,3d63(6.21)
            'multiple':6.25, ## 1hrd,1aup(6.25)
            ## initial R350 transformation
            'virus':5.10, ## 2fsy,1ohg (dilemma: must be above 5.1 when 2fsy v 1ohg, must be below 4.2 when 2b2g v 1zdi; solution: align chains three-wise)
            }
        self.max_rmsd_wrong = 9.5 ## 2eu1 ## assume error if above this rmsd

        self.minres = 5.0 ## minimum resolution

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25
        self.max_mutations = 10

        self.path_pdb = '/data/pdb-v3.2'
        self.topdir = '/local/tc/quakes'

##        self.path_pdb = '/data/pdb_redo'
##        self.topdir = '/local/tc/quakes/pdb_redo'

##        self.path_pdb = '/media/39b8dbfd-b37f-472f-a030-c86a6e75c0d9/1TB/pdb'
##        self.path_mmCIF = '/media/39b8dbfd-b37f-472f-a030-c86a6e75c0d9/1TB/mmCIF'
##        self.path_biounits = '/media/39b8dbfd-b37f-472f-a030-c86a6e75c0d9/1TB/biounit'

        self.path_pdb = '/media/WDMyBook1TB/2TB/pdb'
        self.path_mmCIF = '/media/WDMyBook1TB/2TB/mmCIF'
        self.path_biounits = '/media/WDMyBook1TB/2TB/biounit'

        self.topdir = '/home/tc/UCD/quakes'

        self.path_cwd = os.getcwd()
        
        return

if __name__ == '__main__':
    instance_quakes = quakes()
    instance_quakes.main()

## answers
##
## q6
## compare structures if coordinates of atoms/residues missing? yes
##
## q7
## only compare chains of identical length ?
## or accept terminal appendixes (e.g. 1aon:a vs 1pf9:a) of a certain length
## (e.g. max 25 res or 10% of total length)?
## yes
##
## pdbs with different hetero compounds are not compared
## pdbs with different nonhetero compunds (peptide, nucleotide, saccharide) are not compared (e.g. 1ok7 vs 1mmi)

## structural alignment is only performed between long peptide chains

## ions are ignored, but sodium cause different oligomerization states (e.g. 1gt2.pdb vs 1o6z.pdb)

## large rmsd when multiple chains can be caused by different transformation matrices

            ## different oxidation states
            ## 1oc3:C

            ## different potassium ("2" vs "6") concentrations
            ## 2hvj FORMUL 4 K 2(K1 1+)
            ## 2hvk FORMUL 4 K 6(k1 1+)

            ## different quarternary structures
            ## 1c77, 1c78, 1c79

            ## incorrect remark 350 translation vector
            ## 1bks,2wsy

## 2007jun20
## It came to my attention that some chains have less than 95% sequence similarity, only because they are of different length. (e.g. 1jtn.pdb, 1qth.pdb)

## 2007jul05
## it occured to me that the number of hetero atoms can differ between two pdbs if different asymmetric units are given for the same biological unit in the two pdbs

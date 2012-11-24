#!/software/bin/python
#
#$Id: quakes.py,v 1.30 2007/06/02 14:23:13 tc Exp $
#
#Tommy Carstensen, University College Dublin, 2007

##
## questions of interest

## find correlation (if any) between RMSD of backbone/all atoms for all/neighboring(exponentional sphere radii)/surface residues

## rmsd not mathematically dependent on chain length! but rmsd physically dep on chain length? plot!

##
## questions of concern

## q1
## do not accept residue modifications (e.g. ser/thr phosphorylation)
## or consider residue modifications to be equal to mutations ?
## look up MODRES records..? LINK records???

## what to do with proteins with different saccharides?
## e.g. 1q8p LINK O3 (alpha 1-3 bond) and 1q8o LINK O2 (alpha 1-2 bond)
## what do to if LINK of sugar to protein (e.g. 1cb2)?
## what to do if sugar chains with same monomers but different lengths... same hetIDs but different structures...

## write a new faster seq aln alg

## if no remark350 record and high rmsd then try monomers...

## proteins with 10 or more mutations relative to wt
## ['7adh','1xac','1xad','1cx6','174l','1d3n','1hhl','192l','1a6i','1fbi','1lz2','1jhl','1sbt','2sbt']

## what to do with 1ft8:E and 1koh:B for which coordinates of aligned residues are not given???

## check identity of nonpeptides and small peptides *before* checking similarity of long peptides (speed upgrade!)

## solve the problem of too many chain combinations by sequential pairing
## this will yield a maximum of n! combinations to check
## sequential addition of sequence similar chains to reduce number of combinations further

class quakes:

    def main(self):

##        self.rsync()
##        self.gunzip()

#### out (rewrite and group paired pdbs)
##        import os ## temporary
##        self.l_pdbs = set()
##
##        dirs = os.listdir('/oxygenase_local/tc/quakes/pdb/')
##        for dir in dirs:
##            files = os.listdir('/oxygenase_local/tc/quakes/pdb/%s' %(dir))
##            for file in files:
##
##                self.l_pdbs |= set([file[0:4]]) ##
##                self.l_pdbs |= set([file[6:10]]) ##
##
##                if file[0:4] == file[6:10]:
##                    os.remove('/oxygenase_local/tc/quakes/pdb/%s/%s' %(dir,file))
##
##                if os.path.isfile('/oxygenase_local/tc/quakes/out/%s' %(file)):
##                    os.remove('/oxygenase_local/tc/quakes/out/%s' %(file))
##
####        files = os.listdir('/oxygenase_local/tc/quakes/out/')
####        for file in files:
####            time = os.path.getmtime('/oxygenase_local/tc/quakes/out/%s' %(file))
####
####            self.l_pdbs |= set([file[0:4]])
####            self.l_pdbs |= set([file[6:10]])
####
##        self.l_pdbs = list(self.l_pdbs)

#### clusters95
##
##        import os
##        fd = open('clusters95.txt','r')
##        lines = fd.readlines()
##        fd.close()
##        clusters = {}
##        for line in lines:
##            cluster = line.split()[0]
##            pdb = line.split()[2][:4].lower()
##            if cluster not in clusters.keys():
##                clusters[cluster] = []
##            clusters[cluster] += [pdb]
##        l_clusters = clusters.keys()
##        l_clusters.sort()
##        for i in range(len(l_clusters)):
##            if i < 10778:
##                continue
##            if i > 12001:
##                continue
##            if i == 12001:
##                i == 10777
##            cluster = l_clusters[i]
##            print '---------%s/%s--------------' %(i, len(l_clusters))
##            print cluster
##            self.l_pdbs = list(set(clusters[cluster]))
##            pdbs = self.l_pdbs
##            for pdb in pdbs:
##                if not os.path.isfile('/oxygenase_local/data/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb)):
##                    self.l_pdbs.remove(pdb)
##            try:
##                self.l_pdbs.remove('1fng')
##            except:
##                None
##            try:
##                self.l_pdbs.remove('1fne')
##            except:
##                None
##            if len(self.l_pdbs) == 1:
##                continue
##            self.pdbcount = len(self.l_pdbs)
##            d_seq, d_chains_intrapdb_sequence_identical = self.parse_sequences()
##            d_rmsd = self.analyze_sequences(d_seq, d_chains_intrapdb_sequence_identical)
##            self.write_rmsd_to_file(d_rmsd, d_seq, prefix='rmsd')

#### windows
##        import os ## temporary
##        pdbs = os.listdir('/oxygenase_local/data/pdb/')
##        self.l_pdbs = []
##        for pdb in pdbs:
##            self.l_pdbs += [pdb[:4]]

######## pair
##        self.l_pdbs = ['3btg','3btw'] ## temporary

#### temporary/toomany
        self.l_pdbs = set(['1a8r', '1a9c', '1acm', '1aq3', '1aq4', '1auy', '1bfr', '1bkh', '1c9s', '1d09', '1ezz', '1f1b', '1f30', '1f33', '1f9c', '1fbx', '1g3k', '1gh0', '1grl', '1gtf', '1gtn', '1gtp', '1ha7', '1hq6', '1hqy', '1ht1', '1ht2', '1i4y', '1i4z', '1i5o', '1ibt', '1ibu', '1j4z', '1jjw', '1jts', '1kpo', '1l8h', '1l8i', '1m8q', '1muc', '1mvw', '1n3t', '1n8j', '1o19', '1o1a', '1o1b', '1o1d', '1o1e', '1o1f', '1o1g', '1oel', '1ohg', '1oyx', '1oz3', '1q95', '1qaw', '1ss8', '1tth', '1tu0', '1upm', '1upp', '1utd', '1utf', '1utv', '1uvh', '1veq', '1w39', '1xjw', '1yep', '1yex', '1yf0', '1yf1', '1zdi', '2a0f', '2b2e', '2b2g', '2bny', '2bq5', '2c8t', '2cby', '2cwx', '2cxe', '2d69', '2eu1', '2eul', '2fkp', '2frp', '2fs3', '2fz1', '2g33', '2g34', '2gdq', '2gge', '2ggg', '2ggi', '2ggj', '2h66', '2hms', '2hmt', '2hmu', '2hmw', '2htn', '2i81', '2jd6', '2jd7', '2jd8', '2muc', '2nwc', '3at1', '3muc', '6at1', '8atc', '8ruc', '9atc'])
        fd = open('toomanytransformations.txt','r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            self.l_pdbs |= set([line[0:4]])
            self.l_pdbs |= set([line[5:9]])
        self.l_pdbs = list(self.l_pdbs)
        self.l_pdbs.sort()
        print len(self.l_pdbs)
        stop

        ## removal
        try:
            ## remove pdbs spanning multiple pdb files/IDs
            self.l_pdbs.remove('1p0t'); self.l_pdbs.remove('1otz')
        except:
            None
        for pdb in [
            '1ny7', ## typing errors (1ny7, REMARK465 vs SEQRES)
            '1ft8', ## too many missing coordinates
            ]:
            try:
                self.l_pdbs.remove(pdb)
            except:
                None

        ## sort and count after removal
        self.l_pdbs.sort()
        self.pdbcount = len(self.l_pdbs)

        d_seq, d_chains_intrapdb_sequence_identical = self.parse_sequences()
        d_rmsd = self.analyze_sequences(d_seq, d_chains_intrapdb_sequence_identical)
        self.write_rmsd_to_file(d_rmsd, d_seq, prefix='rmsd')

        return


    def gunzip(self):

        import os

        self.l_pdbs = []

        subdirs = os.listdir(self.pdbpath)
        subdirs.sort()
        for subdir in subdirs:
            files = os.listdir(self.pdbpath+subdir)
            for file in files:
                if file[-2:] == 'gz':
                    ## gunzip
                    if os.path.isfile('%s%s/%s' %(self.pdbpath,subdir,file[:-3])):
                        os.remove('%s%s/%s' %(self.pdbpath,subdir,file[:-3]))
                    os.system('gunzip %s%s/%s' %(self.pdbpath,subdir,file))
                self.l_pdbs += ['%s' %(file[3:7])]

        self.l_pdbs = list(set(self.l_pdbs))

        return


    def rsync(self):

        import os

        os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ %s' %(self.pdbpath))

        return


    def analyze_sequences(self, d_seq, d_chains_intrapdb_sequence_identical):

        import os, time

        print 'analyzing sequences'

        ## do not exclude pdb2 from pdb1 loop if sequence similar to previous pdb1
        ## since pdb A == B, A != C, B == C
        ## but do not analyze sequence of pdb A,B and then pdb B,A since A == B and B == A are equivalent

        d_rmsd_identical = {}
        d_representative_chains = {}

        if not os.path.isfile('status.txt'):
            fd = open('status.txt','a')
            fd.write('0\n')
            fd.close()
        i1status = 0
            
        ##
        ## loop 1 over pdbs
        ##
        t1 = time.clock()
        for i1 in range(self.pdbcount-1):

            if i1 in range(0,self.pdbcount,1000) and self.pdbcount > 10000:

                fd = open('status.txt','r')
                lines = fd.readlines()
                fd.close()
                i1status = int(lines[-1][:-1])

                if i1 == i1status:

                    fd = open('status.txt','a')
                    fd.write('%s\n' %(i1status+1000))
                    fd.close()

            ## continue if status updated while skipping i1 values
            if i1 < i1status and self.pdbcount > 10000:
                continue

            self.pdb1 = pdb1 = self.l_pdbs[i1]

            skippdb = self.pdbskip(d_seq, pdb1)
            if skippdb == True:
                continue

            ## identify biomolecule(s)
            d_biomolecules1 = self.identify_biomolecule(pdb1, d_seq)

            ## parse coordinate section and related records
            d_pdb = {}
            d_pdb = self.parse_coordinates(pdb1, d_pdb)

            ##
            ## loop 2 over pdbs
            ##
            for i2 in range(i1+1,self.pdbcount):
                self.pdb2 = pdb2 = self.l_pdbs[i2]
                SEQRESchains2 = d_seq[pdb2]['chains'].keys()

                skippdb = self.pdbskip(d_seq, pdb2)
                if skippdb == True:
                    continue

                ## identify biomolecule(s)
                d_biomolecules2 = self.identify_biomolecule(pdb2, d_seq)


                ##
                ## print status
                ##
                t2 = time.clock()
                if t2-t1 > self.time_status_update or i2 == i1+1:
                    print 'analyzing sequence of %s (%5i/%5i) and %s (%5i/%5i)' %(pdb1, i1+1, self.pdbcount, pdb2, i2+1, self.pdbcount)
                    t1 = t2


                ##
                ## loop over biomolecule(s) of pdb1
                ##
                for biomolecule1 in d_biomolecules1.keys():
                    bmchains1 = d_biomolecules1[biomolecule1]['chains']
                    bmpolymercount1 = d_biomolecules1[biomolecule1]['polymercount']

                    ##
                    ## loop over biomolecule(s) of pdb2
                    ##
##                    try: ## temporary comment out...
                    if 'a' == 'a': ## temporary
                        for biomolecule2 in d_biomolecules2.keys():

                            bmchains2 = d_biomolecules2[biomolecule2]['chains']
                            bmpolymercount2 = d_biomolecules2[biomolecule2]['polymercount']


                            ##
                            ## skip if different number of chains in the biomolecule
                            ##
                            if bmpolymercount1 != bmpolymercount2:
##                                print biomolecule1, biomolecule2, bmpolymercount1, bmpolymercount2, bmchains1, bmchains2 ## temporary
                                continue


                            ##
                            ## skip if different hetero compounds
                            ##
                            d_hetIDs = {
                                pdb1:set(),
                                pdb2:set(),
                                }
                            ## parse hetIDs from HET records
                            for pdb in d_hetIDs:
                                for chain in d_seq[pdb]['HET']:
                                    d_hetIDs[pdb] |= d_seq[pdb]['HET'][chain]
                                    d_hetIDs[pdb] -= set(self.d_ions.keys())
                                    d_hetIDs[pdb] -= set(self.l_modres)
                            ## compare hetero compounds assuming the following two statements to be correct
                            ## 1) "A particular HET group is represented in the PDB archives with a *unique* hetID."
                            ## 2) Depositors specify *all* hetero atoms observed in the electron density map.
                            if d_hetIDs[pdb1] != d_hetIDs[pdb2]:
##                                print d_hetIDs[pdb1] - d_hetIDs[pdb2] ## temporary
##                                print d_hetIDs[pdb2] - d_hetIDs[pdb1] ## temporary
                                continue


## check identity of nonpeptide polymers and small peptides *before* checking similary of long peptides!!! to save time!!!
                            ##
                            ## identify sequence similar chains between pdbs (long peptides only)
                            ##
                            d_chains_interpdb_sequence_similar, mutations = self.identify_similar_chains_from_sequence_inter(
                                d_seq,
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
                            ##
                            (
                                bmSEQRESchains1_not_similar_to_SEQRESchains2,
                                bmSEQRESchains2_not_similar_to_SEQRESchains1
                                ) = self.identify_chains_interpdb_not_sequence_similar(
                                    pdb1, pdb2, bmchains1, bmchains2,
                                    d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
                                    d_seq
                                    )

                            ## check that non sequence similar chains (if any)
                            ## 1) are not long peptides and
                            ## 2) are sequence identical
                            if len(bmSEQRESchains1_not_similar_to_SEQRESchains2) > 0 or len(bmSEQRESchains2_not_similar_to_SEQRESchains1) > 0:

                                ## 1) check that non sequence similar chains (if any) are not long peptides
                                d_bmSEQRESchains_not_similar_to_SEQRESchains = {
                                    pdb1:{'chains1':bmSEQRESchains1_not_similar_to_SEQRESchains2,'chains2':bmSEQRESchains2_not_similar_to_SEQRESchains1,'pdb2':pdb2},
                                    pdb2:{'chains1':bmSEQRESchains2_not_similar_to_SEQRESchains1,'chains2':bmSEQRESchains1_not_similar_to_SEQRESchains2,'pdb2':pdb1},
                                    }
##                                print bmSEQRESchains1_not_similar_to_SEQRESchains2 ## temporary
##                                print bmSEQRESchains2_not_similar_to_SEQRESchains1 ## temporary
                                for pdb in d_bmSEQRESchains_not_similar_to_SEQRESchains.keys():
                                    if len(d_bmSEQRESchains_not_similar_to_SEQRESchains[pdb]['chains1']) > 0:
                                        for chain in d_bmSEQRESchains_not_similar_to_SEQRESchains[pdb]['chains1']:
                                            peptide = False
                                            long = False
                                            ## check if peptide
                                            if d_seq[pdb]['chains'][chain]['type'] == 'peptide':
                                                peptide = True
                                            ## check if long chain
                                            if len(d_seq[pdb]['chains'][chain]['seq']) > self.min_len_chain:
                                                long = True
                                            if peptide == True and long == True:
                                                break
                                        if peptide == True and long == True:
                                            break
                                ## continue if long peptide
                                if peptide == True and long == True:
##                                    print chain, d_seq[pdb1]['chains'][chain]['seq'] ## temporary
##                                    print chain, d_seq[pdb2]['chains'][chain]['seq'] ## temporary
                                    continue

                                ## 2) check that nucleotides/saccharides or short peptides are sequence *identical*
## not true that identical for saccharide of 1q8p.pdb,1q8o.pdb (check LINK records)
                                for p1 in d_bmSEQRESchains_not_similar_to_SEQRESchains.keys():
                                    p2 = d_bmSEQRESchains_not_similar_to_SEQRESchains[p1]['pdb2']
                                    sequence_identical = True
                                    for chain1 in d_bmSEQRESchains_not_similar_to_SEQRESchains[p1]['chains1']:
                                        peptide = False
                                        long = False
                                        ## check if peptide
                                        if d_seq[p1]['chains'][chain1]['type'] == 'peptide':
                                            peptide = True
                                        ## check if long chain
                                        if len(d_seq[p1]['chains'][chain1]['seq']) > self.min_len_chain:
                                            long = True
                                        if (peptide == True and long == False) or peptide == False:
                                            sequence_identical = False
                                            for chain2 in d_bmSEQRESchains_not_similar_to_SEQRESchains[p1]['chains2']:
                                                if d_seq[p1]['chains'][chain1]['seq'] == d_seq[p2]['chains'][chain2]['seq']:
                                                    sequence_identical = True
                                                if sequence_identical == True:
                                                    break
                                            if sequence_identical == False:
##                                                print p1, chain1 ## temporary
                                                break
                                    ## continue if not sequence identical
                                    if sequence_identical == False:
                                        break
                                ## continue if not sequence identical
                                if sequence_identical == False:
                                    continue
                                    
                                        
                            ##
                            ## identify equivalent chains (interpdb) from structure
                            ##

                            ## parse coordinates
                            d_pdb = self.parse_coordinates(pdb2, d_pdb)

                            ## identify equivalent chains
                            d_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.identify_interpdb_equivalent_chains_from_structure(
                                pdb1, pdb2,
                                d_chains_intrapdb_sequence_identical,
                                d_chains_interpdb_sequence_similar,
                                d_pdb, d_seq,
                                biomolecule1, biomolecule2,
                                d_biomolecules1, d_biomolecules2,
                                )

                            ##
                            ## continue if there are no equivalent chains between the two pdbs
                            ##
                            if d_equivalent_chains == {}:
                                continue
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
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['residues'] = n_residues
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['coordinates'] = n_coordinates
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['d_equivalent_chains'] = d_equivalent_chains
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['mutations'] = mutations
                            if self.pdbcount < 10000:
                                d_quickrmsd = {pdb1:{biomolecule1:{pdb2:{biomolecule2:{
                                    'rmsd':rmsd,'chains':n_chains,'residues':n_residues,'coordinates':n_coordinates,'d_equivalent_chains':d_equivalent_chains,'mutations':mutations,
                                    }}}}}
                                self.write_rmsd_to_file(d_quickrmsd, d_seq, prefix='quickrmsd')
                                d_biomolecules = {
                                    pdb1:{'biomolecule':biomolecule1},
                                    pdb2:{'biomolecule':biomolecule2},
                                    }
                            self.rmsd2bfactor(pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_pdb, d_seq, tv1, rm, tv2, d_equivalent_chains, bmchains1, bmchains2)

##                    except:
##                        import sys
##                        fd = open('errors.txt','a')
##                        fd.write('%s %s %s\n' %(pdb1, pdb2, sys.exc_info()))
##                        fd.close()
##                        continue

        return d_rmsd_identical


    def pdbskip(self, d_seq, pdb):

        pdbskip = False
        SEQRESchains = d_seq[pdb]['chains'].keys()

        ## continue if superseded structure
        if 'SPRSDE' in d_seq[pdb].keys():
            pdbskip = True
            return pdbskip

        ## continue if no polymer chains (e.g. 1qd8.pdb)
        if SEQRESchains == []:
            pdbskip = True
            return pdbskip

        ## continue if no (long) protein chains
        if d_seq[pdb]['proteinchains'] == []:
            pdbskip = True
            return pdbskip

        ## continue if NMR or EM structure
        if d_seq[pdb]['EXPDTA'] in ['NMR','CRYO-ELECTRON MICROSCOPY','ELECTRON TOMOGRAPHY']:
##            print pdb, d_seq[pdb]['EXPDTA'] ## temporary
            pdbskip = True
            return pdbskip

        ## continue if multiple models (NMR or non-NCS-averaged x-ray structures - e.g. 1hto.pdb)
        if 'MODEL' in d_seq[pdb].keys():
            pdbskip = True
            return pdbskip

        return pdbskip


    def rmsd2bfactor(self, pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_pdb, d_seq, tv1, rm, tv2, d_equivalent_chains, bmchains1, bmchains2):

        biomolecule1 = str(biomolecule1).zfill(2)
        biomolecule2 = str(biomolecule2).zfill(2)

        import os

        pdblines1 = []
        pdblines2 = []
        for rep_chain1 in d_equivalent_chains.keys():
            chains1 = d_equivalent_chains[rep_chain1][0]
            chains2 = d_equivalent_chains[rep_chain1][1]
            if len(chains1) != len(chains2):
                print pdb1, pdb2
                notexpected
            for i in range(len(chains1)):
                chain1 = chains1[i]
                chain2 = chains2[i]
                d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(d_pdb, d_seq, pdb1, pdb2, chain1, chain2)
                coordinates1, coordinates2, rescount, lines1, lines2 = self.ATOMrecords2coordinates(d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2, rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2)
                pdblines1 += lines1
                pdblines2 += lines2

##        print d_equivalent_chains
##        print bmchains1, bmchains2
##        stop

        fd = open('pdb/%s/%s%s%s%s.pdb' %(pdb1[1], pdb1, str(biomolecule1).zfill(2), pdb2, str(biomolecule2).zfill(2)), 'w')
        fd.writelines(pdblines1)
        fd.close()
        fd = open('pdb/%s/%s%s%s%s.pdb' %(pdb2[1], pdb2, str(biomolecule2).zfill(2), pdb1, str(biomolecule1).zfill(2)), 'w')
        fd.writelines(pdblines2)
        fd.close()

        ##
        ## gif thumbnails
        ##
        d_biomolecules = {pdb1:biomolecule1,pdb2:biomolecule2}
        for pdb in d_biomolecules.keys():
            biomolecule = d_biomolecules[pdb]
            if pdb == pdb1:
                prefix = pdb1+str(biomolecule1).zfill(2)+pdb2+str(biomolecule2).zfill(2)
            if pdb == pdb2:
                prefix = pdb2+str(biomolecule2).zfill(2)+pdb1+str(biomolecule1).zfill(2)
            ## write rasmol script
            lines = [
                'rasmol -nodisplay pdb/%s/%s.pdb << EOF\n' %(pdb[1], prefix),
                'color temperature\n',
                'spacefill\n',
                'write tmp/%s.ppm\n' %(prefix),
                'exit\n',
                ]
            ## write rasmol script to file
            fd = open('tmp/%srasmol.src' %(prefix),'w')
            fd.writelines(lines)
            fd.close()
            ## execute rasmol script
            os.system('source tmp/%srasmol.src > tmp/%srasmol.log' %(prefix,prefix))
            ## convert rasmol output
            os.system('convert tmp/%s.ppm -resize x80 gif/%s/%s.gif' %(prefix,pdb[1],prefix))
            ## clean up
            os.remove('tmp/%s.ppm' %(prefix))
            os.remove('tmp/%srasmol.log' %(prefix))
            os.remove('tmp/%srasmol.src' %(prefix))
        
        return

        
    def write_rmsd_to_file(self, d_rmsd, d_seq, prefix):

        import os

        ## sorted list of parameters
        l_columns_html = ['gif1','gif2','pdb1', 'pdb2', 'bm1', 'bm2', 'rmsd', 'mutations', 'chains', 'residues', 'coordinates', 'pH1', 'pH2', 'T1', 'T2', 'res1', 'res2','spacegroup1', 'spacegroup2', 'title1','title2','hetIDs1', 'hetIDs2']
        ## keys=htmlkeys, values=htmlcolumns (table headings)
        d_columns_headers = {
            'gif1':'gif1','gif2':'gif2',
            'pdb1':'pdb1', 'pdb2':'pdb2',
            'bm1':'bm1', 'bm2':'bm2',
            'rmsd':'<a href="http://en.wikipedia.org/wiki/Protein_structural_alignment">rmsd</a>',
            'chains':'chains', 'residues':'residues', 'coordinates':'coordinates',
            'pH1':'pH1', 'pH2':'pH2', 'T1':'T1', 'T2':'T2',
            'res1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">res1</a>',
            'res2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">res2</a>',
            'spacegroup1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">spacegroup1</a>',
            'spacegroup2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">spacegroup2</a>',
            'title1':'title1','title2':'title2',
            'hetIDs1':'hetIDs1', 'hetIDs2':'hetIDs2',
            'mutations':'mutations',
            }

        ## initiate html lines
        th = '<tr>\n'
        for column in l_columns_html:
            th += '<td>%s</td>\n' %(d_columns_headers[column])
        th += '</tr>\n'

        l_tr = ''

        d_html = {}

        path = 'file:///oxygenase_local/tc/quakes/'

        ##
        ## loop over pdbs and biomolecules
        ##
        for pdb1 in d_rmsd:

            ## parse physiochemical properties (pH and temperature)
            if 'REMARK200' in d_seq[pdb1].keys():
                T1 = d_seq[pdb1]['REMARK200']['TEMPERATURE']
                pH1 = d_seq[pdb1]['REMARK200']['PH']
            else:
                T1 = 'N/A'
                pH1 = 'N/A'
            try:
                T1 = '%5.1f' %(float(T1))
            except:
                T1 = '%s' %(T1.rjust(5))
            try:
                pH1 = '%4.1f' %(float(pH1))
            except:
                pH1 = '%s' %(pH1.rjust(4))
            ## parse x-ray space group
            spacegroup1 = d_seq[pdb1]['CRYST1'].rjust(10)
            ## parse hetIDs
            hetIDs1 = set()
            ## parse resolution
            res1 = d_seq[pdb1]['REMARK2']
            try:
                res1 = '%5.2f' %(float(res1))
            except:
                res1 = '%s' %(res1.rjust(5))
            ## parse hetIDs
            for chain in d_seq[pdb1]['HET'].keys():
                hetIDs1 |= d_seq[pdb1]['HET'][chain]
            hetIDs1 = list(hetIDs1)

            for bm1 in d_rmsd[pdb1].keys():

                for pdb2 in d_rmsd[pdb1][bm1]:

                    ## parse physiochemical properties (pH and temperature)
                    if 'REMARK200' in d_seq[pdb2].keys():
                        T2 = d_seq[pdb2]['REMARK200']['TEMPERATURE']
                        pH2 = d_seq[pdb2]['REMARK200']['PH']
                    else:
                        T2 = 'N/A'
                        pH2 = 'N/A'
                    try:
                        T2 = '%5.1f' %(float(T2))
                    except:
                        T2 = '%s' %(T2.rjust(5))
                    try:
                        pH2 = '%4.1f' %(float(pH2))
                    except:
                        pH2 = '%s' %(pH2.rjust(4))
                    ## parse x-ray space group
                    spacegroup2 = d_seq[pdb2]['CRYST1'].rjust(10)
                    ## parse hetIDs
                    hetIDs2 = set()
                    ## parse resolution
                    res2 = d_seq[pdb2]['REMARK2']
                    try:
                        res2 = '%5.2f' %(float(res2))
                    except:
                        res2 = '%s' %(res2.rjust(5))
                    ## parse hetIDs
                    for chain in d_seq[pdb2]['HET'].keys():
                        hetIDs2 |= d_seq[pdb2]['HET'][chain]
                    hetIDs2 = list(hetIDs2)

                    for bm2 in d_rmsd[pdb1][bm1][pdb2].keys():

                        ## parse rmsd and related data
                        rmsd = '%5.2f' %(d_rmsd[pdb1][bm1][pdb2][bm2]['rmsd'])
                        n_chains = '%3i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['chains'])
                        n_residues = '%4i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['residues'])
                        n_coordinates = '%5i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['coordinates'])
                        mutations = '%2i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['mutations'])
##                        d_equivalent_chains = d_rmsd[pdb1][bm1][pdb2][bm2]['d_equivalent_chains']
##
##                        ## convert d_equivalent_chains to s_chains
##                        s_chains1 = ''
##                        s_chains2 = ''
##                        for rep_chain1 in d_equivalent_chains.keys():
##                            l_chains1 = d_equivalent_chains[rep_chain1][0]
##                            l_chains2 = d_equivalent_chains[rep_chain1][1]
##                            for chain in l_chains1:
##                                s_chains1 += chain+','
##                            s_chains1 = s_chains1[:-1]
##                            for chain in l_chains2:
##                                s_chains2 += chain+','
##                            s_chains2 = s_chains2[:-1]

                        ## write data to dictionary
                        htmlhetIDs1 = ''
                        for i in range(len(hetIDs1)):
                            hetID = hetIDs1[i]
                            htmlhetIDs1 += '<a href="http://ligand-depot.rcsb.org/pub/%s/%s.html">%s</a>, ' %(hetID,hetID,hetID)
                            if (i+1) % 3 == 0:
                                htmlhetIDs1 += '<br>'
                        htmlhetIDs1 = htmlhetIDs1[:-2]
                        htmlhetIDs2 = ''
                        for i in range(len(hetIDs2)):
                            hetID = hetIDs2[i]
                            htmlhetIDs2 += '<a href="http://ligand-depot.rcsb.org/pub/%s/%s.html">%s</a>, ' %(hetID,hetID,hetID)
                            if (i+1) % 3 == 0:
                                htmlhetIDs2 += '<br>'
                        htmlhetIDs2 = htmlhetIDs2[:-2]

                        prefix1 = '%s%s%s%s' %(pdb1, str(bm1).zfill(2), pdb2, str(bm2).zfill(2))
                        prefix2 = '%s%s%s%s' %(pdb2, str(bm2).zfill(2), pdb1, str(bm1).zfill(2))
                        d_columns_data = {
                            'title1':'<td style="font-size:80%%">%s' %(d_seq[pdb1]['TITLE']),
                            'title2':'<td style="font-size:80%%">%s' %(d_seq[pdb2]['TITLE']),
                            'gif1':'<td><a href="%spdb/%s/%s.pdb"><img src="%sgif/%s/%s.gif"></a>' %(path,pdb1[1],prefix1,path,pdb1[1],prefix1),
                            'gif2':'<td><a href="%spdb/%s/%s.pdb"><img src="%sgif/%s/%s.gif"></a>' %(path,pdb2[1],prefix2,path,pdb2[1],prefix2),
                            'pdb1':'<td><a href="%shtm/%s.htm">%s</a>' %(path,pdb1,pdb1),
                            'pdb2':'<td><a href="%shtm/%s.htm">%s</a>' %(path,pdb2,pdb2),
                            'bm1':'<td style="text-align: right">%s' %(bm1),
                            'bm2':'<td style="text-align: right">%s' %(bm2),
                            'pH1':'<td style="text-align: right">%s' %(pH1),
                            'pH2':'<td style="text-align: right">%s' %(pH2),
                            'T1':'<td style="text-align: right">%s' %(T1),
                            'T2':'<td style="text-align: right">%s' %(T2),
                            'res1':'<td style="text-align: right">%s' %(res1),
                            'res2':'<td style="text-align: right">%s' %(res2),
                            'spacegroup1':'<td nowrap>%s' %(spacegroup1),
                            'spacegroup2':'<td nowrap>%s' %(spacegroup2),
                            'hetIDs1':'<td style="font-size:80%%" nowrap>%s' %(htmlhetIDs1),
                            'hetIDs2':'<td style="font-size:80%%" nowrap>%s' %(htmlhetIDs2),
                            'rmsd':'<td style="text-align: right">%s' %(rmsd),
                            'chains':'<td style="text-align: right">%s' %(n_chains),
                            'residues':'<td style="text-align: right">%s' %(n_residues),
                            'coordinates':'<td style="text-align: right">%s' %(n_coordinates),
                            'mutations':'<td style="text-align: right">%s' %(mutations)
                            }

                        ## write data to html lines
                        tr = '<tr>\n'
                        for column in l_columns_html:
                            tr += '%s</td>\n' %(d_columns_data[column])
                        tr += '</tr>\n'

                        l_tr += tr

                        if not pdb1 in d_html.keys():
                            d_html[pdb1] = ''
                        d_html[pdb1] += tr
                        if not pdb2 in d_html.keys():
                            d_html[pdb2] = ''
                        d_html[pdb2] += tr
                        
        ## write html to file
        path = ''
        file = '%s.htm' %(prefix)
        self.append_table_rows(path,file,l_tr,th)

        ## write html to file
        if prefix == 'quickrmsd':
            path = 'htm/'
            for pdb in d_html.keys():
                file = '%s.htm' %(pdb)
                l_tr = d_html[pdb]
                self.append_table_rows(path,file,l_tr,th)

        return


    def append_table_rows(self,path, file,l_tr,th):

        import os

        if os.path.isfile('%s%s' %(path,file)):
            fd = open('%s%s' %(path, file),'r')
            html = fd.readlines()[:-1]
            fd.close()
            html += [l_tr+'</table>\n']
        else:
            html = ['<table border="1">\n'+th+l_tr+'</table>\n']
        fd = open('%s%s'%(path,file),'w')
        fd.writelines(html)
        fd.close()
     

    def identify_biomolecule(self, pdb, d_seq):

        SEQRESchains = d_seq[pdb]['chains'].keys()
        if 'REMARK350' in d_seq[pdb].keys():
            d_biomolecules = {}
            biomolecules = d_seq[pdb]['REMARK350'].keys()
            for biomolecule in biomolecules:
                chains = d_seq[pdb]['REMARK350'][biomolecule]['chains'].keys()
                d_biomolecules[biomolecule] = {}
                d_biomolecules[biomolecule]['chains'] = []
                d_biomolecules[biomolecule]['polymercount'] = 0
                for chain in chains:
                    d_biomolecules[biomolecule]['chains'] += [chain]
                    ## only count if not hetero (e.g. not water not mentioned in remark525 records)
                    if chain in SEQRESchains:
                        d_biomolecules[biomolecule]['polymercount'] += len(d_seq[pdb]['REMARK350'][biomolecule]['chains'][chain])
        ## assume everything to be the biomolecule
        else:
            d_biomolecules = {
                '1':{
                    'chains':SEQRESchains,
                    'polymercount':len(SEQRESchains),
                    }
                }

        return d_biomolecules


    def identify_chains_interpdb_not_sequence_similar(
        self, pdb1, pdb2, bmchains1, bmchains2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_seq
        ):

        ## pdb1
        SEQRESchains1_similar_to_SEQRESchains2 = set()
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            SEQRESchains1_similar_to_SEQRESchains2 |= set([rep_chain1])
            SEQRESchains1_similar_to_SEQRESchains2 |= set(d_chains_intrapdb_sequence_identical[pdb1][rep_chain1])
        bmSEQRESchains1_similar_to_SEQRESchains2 = SEQRESchains1_similar_to_SEQRESchains2 & set(bmchains1)

        ## pdb2
        SEQRESchains2_similar_to_SEQRESchains1 = set()
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']
            SEQRESchains2_similar_to_SEQRESchains1 |= set(rep_chain2)
            SEQRESchains2_similar_to_SEQRESchains1 |= set(d_chains_intrapdb_sequence_identical[pdb2][rep_chain2])
        bmSEQRESchains2_similar_to_SEQRESchains1 = SEQRESchains2_similar_to_SEQRESchains1 & set(bmchains2)

        SEQRESchains1_not_similar_to_SEQRESchains2 = set(d_seq[pdb1]['chains'].keys()) & (set(bmchains1) - SEQRESchains1_similar_to_SEQRESchains2)
        SEQRESchains2_not_similar_to_SEQRESchains1 = set(d_seq[pdb2]['chains'].keys()) & (set(bmchains2) - SEQRESchains2_similar_to_SEQRESchains1)

        return SEQRESchains1_not_similar_to_SEQRESchains2, SEQRESchains2_not_similar_to_SEQRESchains1


    def calculate_rmsd_for_multiple_chains(
        self,chains1,chains2,d_pdb,pdb1,pdb2,d_seq,
        verbose=True
        ):

##        print 'calculating rmsd for chains %s, %s of %s, %s' %(chains1, chains2, pdb1, pdb2)

        import Numeric
        import sys
        sys.path.append('/home/people/tc/python/Protool/')
        import geometry
        instance_geometry = geometry.geometry()

        if len(chains1) != len(chains2):
            print pdb1, pdb2
            print chains1, chains2
            notexpected

        coordinates1 = []
        coordinates2 = []
        residue_count = 0

        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            ## parse coordinates
            coords1, coords2, rescount = self.dcoordinates2lcoordinates(
                d_pdb,d_seq,
                pdb1, pdb2,
                chain1, chain2,
                )

            ## append coordinates
            coordinates1 += coords1
            coordinates2 += coords2
            residue_count += rescount

        rmsd = instance_geometry.superpose(coordinates1,coordinates2)
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter
        
        if verbose == True:
            if len(chains1) < 60 and len(chains2) < 60:
                print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s, %s, %s' %(rmsd, len(chains1), residue_count, len(coordinates1), chains1, chains2, pdb1, pdb2)
            else:
                print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s' %(rmsd, len(chains1), residue_count, len(coordinates1), pdb1, pdb2)

        return rmsd, len(chains1), residue_count, len(coordinates1), tv1, rm, tv2


    def identify_interpdb_equivalent_chains_from_structure(
        self, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_pdb, d_seq,
        biomolecule1, biomolecule2,
        d_biomolecules1, d_biomolecules2,
        ):

        bmchains1 = d_biomolecules1[biomolecule1]['chains']
        bmchains2 = d_biomolecules2[biomolecule2]['chains']
        bmpolymercount1 = d_biomolecules1[biomolecule1]['polymercount']
        bmpolymercount2 = d_biomolecules2[biomolecule2]['polymercount']

        import os

        ## return the following dictionary structure
        ## equivalent: {repchain1:[chains1,chains2]}

    ##    print 'identify_interpdb_equivalent_chains_from_structure', pdb1, pdb2

        n_chains = n_residues = n_coordinates = 0
        d_equivalent_chains = {}
        tv1 = 'N/A'
        rm = 'N/A'
        tv2 = 'N/A'

        spacegroup1 = d_seq[pdb1]['CRYST1']
        spacegroup2 = d_seq[pdb2]['CRYST1']
        crystalsystem1 = self.d_crystalsystems[spacegroup1]
        crystalsystem2 = self.d_crystalsystems[spacegroup2]

        l_chains1 = []
        l_chains2 = []
        chains1 = []
        chains2 = []

        ## group sequence similar chains to avoid a large number of permutations causing comparison of nonsequenceidentical chains
        ## instead do permutation of sequence identical chains and combination of permutations

        ## loop over representative chains1
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            ## identify rep chain2 seq sim to rep chain1
            rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']
            ## identify chains identical to the representative chains
            chains2seqid = [rep_chain2]+d_chains_intrapdb_sequence_identical[pdb2][rep_chain2]
            chains1seqid = d_chains_intrapdb_sequence_identical[pdb1][rep_chain1]+[rep_chain1]
            ## exclude chains which are not in the biomolecule
            bmchains1seqid = list(set(chains1seqid) & set(bmchains1))
            bmchains2seqid = list(set(chains2seqid) & set(bmchains2))
            ## sort chains
            bmchains1seqid.sort()
            bmchains2seqid.sort()
            ## append chains
            if bmchains1seqid != []:
                if bmchains1seqid not in l_chains1: ## e.g. 1sxa.pdb (vs 1e9o.pdb)
                    l_chains1 += [bmchains1seqid]
                    chains1 += bmchains1seqid
            if bmchains2seqid != []:
                if bmchains2seqid not in l_chains2: ## e.g. 1sxa.pdb (vs 1e9o.pdb)
                    l_chains2 += [bmchains2seqid]
                    chains2 += bmchains2seqid

        d_biomolecules = {
            pdb1:{'biomolecule':biomolecule1,'l_chains':l_chains1},
            pdb2:{'biomolecule':biomolecule2,'l_chains':l_chains2},
            }

        ##
        ## identical number of chains after REMARK350 transformation
        ##
        ## e.g. A1,B1,A2,B2,A3,B3 == A,B,C,D,E,F of 1xnv.pdb,1xo6.pdb
        ## e.g. B,B,B,B == B,D,B,D of 1vwr.pdb,1vwi.pdb
        ## e.g. A,C = B,B of 1my3.pdb,1mxu.pdb
        if bmpolymercount1 == bmpolymercount2:

## use this to reduce number of 1xnv,1xo6 permutations
##                ## group intrapdb seq id chains to reduced number of possible permutations
##                rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']
##                l_equivalent_bmchains = [
##                    [
##                        set(d_chains_intrapdb_sequence_identical[pdb1][rep_chain1]+[rep_chain1]) & set(bmchains1),
##                        set(d_chains_intrapdb_sequence_identical[pdb2][rep_chain2]+[rep_chain2]) & set(bmchains2),
##                        ]
##                    ]

            for pdb in d_biomolecules.keys():
##                if 'REMARK350' not in d_seq[pdb].keys():
##                    d_biomolecules[pdb]['tchains'] = d_biomolecules[pdb]['chains']
##                    continue
                biomolecule = d_biomolecules[pdb]['biomolecule']
                l_chains = d_biomolecules[pdb]['l_chains']
                l_tchains = []
                for chains in l_chains:
                    tchains = []
                    for chain in chains:
                        if 'REMARK350' not in d_seq[pdb]:
                            tchains += [chain]
                        else:
                            matrix_nos = d_seq[pdb]['REMARK350'][biomolecule]['chains'][chain]
                            for matrix_no in matrix_nos:
                                matrix = d_seq[pdb]['REMARK350'][biomolecule]['matrices'][matrix_no]
                                d_pdb, tchain = self.matrixtransformation(d_pdb,pdb,chain,matrix,matrix_no)
                                tchains += [tchain]
                    l_tchains += [tchains]
                d_biomolecules[pdb]['l_tchains'] = l_tchains

            l_tchains1 = d_biomolecules[pdb1]['l_tchains']
            l_tchains2 = d_biomolecules[pdb2]['l_tchains']
            tchains1 = []
            for tchains in l_tchains1:
                tchains1 += tchains
            tchains2 = []
            for tchains in l_tchains2:
                tchains2 += tchains

            ## check if the expected correct combination of chains gives a low rmsd
            rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(tchains1,tchains2,d_pdb,pdb1,pdb2,d_seq)
            if rmsd < self.maxrmsd or (len(tchains1) == 1 and len(tchains2) == 1):
                d_equivalent_chains[rep_chain1] = [tchains1,tchains2]
                return d_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

            if len(l_tchains1) != len(l_tchains2):
                print pdb1, pdb2
                notexpected
            if len(tchains1) != len(tchains2):
                print pdb1, pdb2
                notexpected

            ##
            ## calculate the number of combinations
            ##
            ncombinations = 1
            
            for tchains in l_tchains2:
                nchains = len(tchains)
                npermutations = self.faculty(nchains)
                ncombinations *= npermutations

            if ncombinations > 1024:
                toomany = []
                if os.path.isfile('toomanytransformations.txt'):
                    fd = open('toomanytransformations.txt','r')
                    lines = fd.readlines()
                    fd.close()
                    for line in lines:
                        toomany += [ [line.split()[0],line.split()[1]] ]
                if not [pdb1,pdb2] in toomany:
                    fd = open('toomanytransformations.txt','a')
                    if len(tchains1) < 16:
                        fd.write('%s %s %s %s %s %s\n' %(pdb1, pdb2, rmsd, ncombinations, tchains1, tchains2))
                    else:
                        fd.write('%s %s %s\n' %(pdb1, pdb2, rmsd))
                    fd.close()
                return d_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

            ##
            ## permutation of chains
            ##
            l_chains2permutations = []
            for i in range(len(l_tchains2)):
                bmchains2seqid = l_tchains2[i]
                chains2permutations = self.permutation(bmchains2seqid)
                l_chains2permutations += [chains2permutations]

            ##
            ## combination of permutation of chains
            ##
            for i in range(len(l_chains2permutations)-1):
                j = i+1
                ipermutations = l_chains2permutations[i]
                jpermutations = l_chains2permutations[j]
                ## combination
                combinations = []
                for k in range(len(ipermutations)):
                    for l in range(len(jpermutations)):
                        combinations += [ipermutations[k]+jpermutations[l]]
                ## replace permutations with combinations
                l_chains2permutations[j] = combinations
                l_chains2permutations[i] = []
            chains2combinations = l_chains2permutations[-1]

            ## identify a correct combination of chains
            minrmsd = ['N/A','N/A']
##            print 'chains2combinations', chains2combinations
            for i in range(len(chains2combinations)):
                chains2combination = chains2combinations[i]
                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(tchains1,chains2combination,d_pdb,pdb1,pdb2,d_seq,verbose=False)
                if n_chains < 60:
                    print '%s/%s %s %s' %(i+1, len(chains2combinations), chains2combination, rmsd)
                else:
                    print '%s/%s %s' %(i+1, len(chains2combinations), rmsd)
                if rmsd < minrmsd[1]:
                    minrmsd = [chains2combination,rmsd]
                    ## break to save time if many permutations
## change 4 to a variable...
                    if rmsd < self.maxrmsd and (len(tchains1) > 3 or i == 0):
                        break
            tchains2 = minrmsd[0]
            rmsd = minrmsd[1]
            d_equivalent_chains[rep_chain1] = [tchains1,tchains2]

        else:
            fd = open('notexpected_differentsized_biounits.txt','a')
            fd.write('%s %s %s %s\n' %(pdb1, pdb2, chains1, chains2))
            fd.close()
            rmsd = 'N/A'
            return d_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2

        return d_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2


    def matrixtransformation(self,d_pdb,s_pdb,chain,matrix,matrix_no):

        ## matrix does not cause transformation
        if matrix == self.nontransformationmatrix:
            return d_pdb,chain

        import Numeric

##        print 'applying REMARK350 transformation to chain %s of pdb %s' %(chain, s_pdb)

        tchains = []


##            tchains += [[]]
##                tchains[-1] += [[]]
##                    tchains[-1][-1] += [tchain]

        rmatrix, tvector = self.transformationmatrix2rotationmatrix_and_translationvector(matrix)

        tchain = chain+'_'+str(matrix_no)

        d_pdb[s_pdb]['chains'][tchain] = {}
        if 'residues' not in d_pdb[s_pdb]['chains'][tchain].keys():
            d_pdb[s_pdb]['chains'][tchain]['residues'] = {}
        for res_no in d_pdb[s_pdb]['chains'][chain]['residues'].keys():
            if res_no not in d_pdb[s_pdb]['chains'][tchain]['residues'].keys():
                d_pdb[s_pdb]['chains'][tchain]['residues'][res_no] = {}
            if 'd_iCodes' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'] = {}
            if 'l_iCodes' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['l_iCodes'] = []
            for iCode in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['l_iCodes']:
                if 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb[s_pdb]['chains'][tchain]['residues'][res_no] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]
                    continue
                if iCode not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'].keys():
                    d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode] = {}
                d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                if 'atoms' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                for atom_name in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                    if 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]
                        continue
                    if atom_name not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                    coord = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate']
                    tcoord = Numeric.matrixmultiply(rmatrix, coord) + tvector
                    d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate'] = tcoord


        return d_pdb, tchain


    def apply_remark350_transformation(self,d_pdb,s_pdb,d_seq,biomolecule,bmchainsseqid):

        import Numeric

##        print 'applying REMARK350 transformation to chain %s of pdb %s' %(chain, s_pdb)

        tchains = []

        ## loop over ordered matrices
        for i in range(len(d_seq[s_pdb]['REMARK350'][biomolecule]['matrices'])):

##            tchains += [[]]
##                tchains[-1] += [[]]
##                    tchains[-1][-1] += [tchain]

            matrix = d_seq[s_pdb]['REMARK350'][biomolecule]['matrices'][i]
            rmatrix, tvector = self.transformationmatrix2rotationmatrix_and_translationvector(matrix)

##            for j in range(len(l_chains)):
##                bmchainsseqid = l_chains[j]
##                tchains += [[]]

            for chain in bmchainsseqid:

                tchain = chain+'_'+str(i+2)
                tchains += [tchain]

                d_pdb[s_pdb]['chains'][tchain] = {}
                if 'residues' not in d_pdb[s_pdb]['chains'][tchain].keys():
                    d_pdb[s_pdb]['chains'][tchain]['residues'] = {}
                for res_no in d_pdb[s_pdb]['chains'][chain]['residues'].keys():
                    if res_no not in d_pdb[s_pdb]['chains'][tchain]['residues'].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no] = {}
                    if 'd_iCodes' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'] = {}
                    if 'l_iCodes' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['l_iCodes'] = []
                    for iCode in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['l_iCodes']:
                        if 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                            d_pdb[s_pdb]['chains'][tchain]['residues'][res_no] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]
                            continue
                        if iCode not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'].keys():
                            d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode] = {}
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                        if 'atoms' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode].keys():
                            d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                        for atom_name in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                            if 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name].keys():
                                d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]
                                continue
                            if atom_name not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                                d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                            coord = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate']
                            tcoord = Numeric.matrixmultiply(rmatrix, coord) + tvector
                            d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate'] = tcoord

        return d_pdb, tchains


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):

        import Numeric

        translationvector = Numeric.array(
            [
                float(transformationmatrix[0][3]),
                float(transformationmatrix[1][3]),
                float(transformationmatrix[2][3]),
                ]
            )

        rotationmatrix = Numeric.array(
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
        self, d_seq, pdb, d_chains_intrapdb_sequence_identical,
        ):

        chains = d_seq[pdb]['chains'].keys()
        waterchains = set(d_seq[pdb]['REMARK525'])

        ## return the following dictionary structure
        ## intrapdb: {pdb:{repchain:[seqidchains]}}

##        print 'identify_identical_chains_from_sequence', pdb1, pdb2

        d_chains_intrapdb_sequence_identical[pdb] = {}

        for i in range(len(chains)-1):
            chaini = chains[i]
            if chaini in waterchains:
                continue

            ## continue if chaini identical to previous chaink
            identical = False
            for chaink in d_chains_intrapdb_sequence_identical[pdb].keys():
                if chaini in d_chains_intrapdb_sequence_identical[pdb][chaink]:
                    identical = True
            if identical == True:
                continue

            ## initiate list for chaini
            if chaini not in d_chains_intrapdb_sequence_identical[pdb].keys():
                d_chains_intrapdb_sequence_identical[pdb][chaini] = []

            d_chains_intrapdb_sequence_identical[pdb][chaini] = []

            seqi = d_seq[pdb]['chains'][chaini]['seq']

            for j in range(i+1,len(chains)):
                chainj = chains[j]
                if chainj in waterchains:
                    continue

                seqj = d_seq[pdb]['chains'][chainj]['seq']
                
                if seqi == seqj:
                    d_chains_intrapdb_sequence_identical[pdb][chaini] += [chainj]
##                    d_chains_intrapdb_sequence_identical[pdb][chainj] = chaini
                ## if last chain in loop and not identical to anything then append to list of representative chains
                elif i == len(chains)-2 and j == len(chains)-1:
                    identical = False
                    for chaink in d_chains_intrapdb_sequence_identical[pdb].keys():
                        if chainj in d_chains_intrapdb_sequence_identical[pdb][chaink]:
                            identical = True
                    if identical == False:
                        d_chains_intrapdb_sequence_identical[pdb][chainj] = []

##        print 'intra', d_chains_intrapdb_sequence_identical
##        stop

        return d_chains_intrapdb_sequence_identical


    def identify_similar_chains_from_sequence_inter(
        self, d_seq, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical,
        bmchains1, bmchains2,
        d_biomolecules1, d_biomolecules2,
        ):

        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        sys.path.append('../../EAT_DB/')
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

            if d_seq[pdb1]['chains'][chain1]['type'] != 'peptide':
                continue

            seq1 = d_seq[pdb1]['chains'][chain1]['seq']

            if len(seq1) < self.min_len_chain:
                continue

            ## only do sequential alignment for representative chains to save time
            for chain2 in repchains2:

                if d_seq[pdb2]['chains'][chain2]['type'] != 'peptide':
                    continue

                seq2 = d_seq[pdb2]['chains'][chain2]['seq']

                if len(seq2) < self.min_len_chain:
                    continue

                if abs(len(seq1)-len(seq2)) > self.max_len_chain_difference:
                    continue

                instance = sequence_alignment.NW(seq1,seq2)
##                print 'aligning chain %s,%s of %s,%s' %(chain1,chain2,pdb1,pdb2)
                s1,s2 = instance.Align(verbose=False)[:2]

                l1 = len(s1)-len(s1.lstrip('-'))
                l2 = len(s2)-len(s2.lstrip('-'))
                l = max(l1,l2)
                r1 = len(s1)-len(s1.rstrip('-'))
                r2 = len(s2)-len(s2.rstrip('-'))
                r = max(r1,r2)
## change .1 to variable...
                if l1 > self.max_len_chain_difference or l1/float(len(s1)) > .1:
                    continue
                if r1 > self.max_len_chain_difference or r1/float(len(s1)) > .1:
                    continue
                if l2 > self.max_len_chain_difference or l2/float(len(s2)) > .1:
                    continue
                if r2 > self.max_len_chain_difference or r2/float(len(s2)) > .1:
                    continue

                s1 = s1[l:len(s1)-r]
                s2 = s2[l:len(s2)-r]

##                if s1 == s2:
##                    d_chains_interpdb_sequence_similar[chain1] = {'rep_chain2':chain2,'l1':l1,'l2':l2}
                n_chainmutations = 0
                l_chainmutations = []
                if '-' in s1 or '-' in s2:
                    continue
                for res in range(len(s1)):
                    res1 = s1[res]
                    res2 = s2[res]
                    if res1 != res2:
                        n_chainmutations += 1
                        l_chainmutations += [[res,res1,res2]]
                    if n_chainmutations > self.max_mutations:
                        break
                if n_chainmutations <= self.max_mutations:
                    if chain1 not in d_chains_interpdb_sequence_similar.keys():
                        d_chains_interpdb_sequence_similar[chain1] = []
                    d_chains_interpdb_sequence_similar[chain1] += [{
                        'rep_chain2':chain2,'l1':l1,'l2':l2,
                        'n_mutations':n_chainmutations, 'l_mutations':l_chainmutations
                    }]
    
        ## assign sequence similar chains correctly and count mutations (e.g. 1u9i.pdb vs 2gbl.pdb)
        ## this will not work if equal amount of mutations in sequence similar chains... (e.g. GAVLI, GAVLG, GAVGI)
        n_mutations = 0
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            min_mutations = ['N/A','N/A']
            for i in range(len(d_chains_interpdb_sequence_similar[rep_chain1])):
                mutations = d_chains_interpdb_sequence_similar[rep_chain1][i]['n_mutations']
                if mutations == min_mutations[1]:
                    l2a = d_chains_interpdb_sequence_similar[rep_chain1][i]
                    l2b = d_chains_interpdb_sequence_similar[rep_chain1][min_mutations[0]]
                    if l2a == l2b:
                        rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1][i]['rep_chain2']
                        print pdb1, pdb2, rep_chain1, rep_chain2, mutations
                        print min_mutations
                        print rep_chain1
                        print d_chains_interpdb_sequence_similar[rep_chain1]
                        expected
                    if l2a < l2b:
                        min_mutations = [i,mutations]
                if mutations < min_mutations[1]:
                    min_mutations = [i,mutations]
            i = min_mutations[0]
            n_mutations += min_mutations[1]
            d_chains_interpdb_sequence_similar[rep_chain1] = d_chains_interpdb_sequence_similar[rep_chain1][i]
                
##        print 'inter', d_chains_interpdb_sequence_similar

        return d_chains_interpdb_sequence_similar, n_mutations


    def res_name2res_symbol(self, res_name):

        if res_name in self.d_res.keys():
            symbol = self.d_res[res_name]
        else:
            symbol = 'X'

        return symbol
    

    def ATOM2seq(self, d_pdb, pdb, chain, SEQRESseq):

        d_res_nos = {}
        seq = ''
        ATOMrespos = 0
        res_nos = d_pdb[pdb]['chains'][chain]['residues'].keys()
        res_nos.sort()
        for i in range(len(res_nos)):
            res_no = res_nos[i]

            for i in range(len(d_pdb[pdb]['chains'][chain]['residues'][res_no]['l_iCodes'])):
                iCode = d_pdb[pdb]['chains'][chain]['residues'][res_no]['l_iCodes'][i]
                d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':iCode}
                res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                seq += self.res_name2res_symbol(res_name)
                ATOMrespos += 1

        return seq, d_res_nos


    def append_missing_residues_to_sequence(self, ATOMseq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq):

        for SEQRESpos in range(SEQRESrange1,SEQRESrange2):
            ATOMseq += SEQRESseq[SEQRESpos]
            d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}

        return ATOMseq, d_res_nos_SEQRES

    
    def identify_missing_nonterminal_residues(self, d_pdb, pdb, chain, SEQRESseq, ATOMseqgaplen, d_res_nos_ATOM):

        d_res_nos_SEQRES = {}
        seq = ''

        ## append N-terminal gap between ATOMseq and SEQRESseq
        SEQRESrange1 = 0
        SEQRESrange2 = ATOMseqgaplen
        seq, d_res_nos_SEQRES = self.append_missing_residues_to_sequence(seq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq)

        ## reset counters and indexes
        ATOMseqrespos = ATOMseqgaplen
        prevATOMseqrespos = 0
        index2 = 0
        ## initiate loop
        SEQRESrange1 = ATOMseqgaplen
        SEQRESrange2 = len(d_res_nos_ATOM.keys())+ATOMseqgaplen
        for SEQRESpos in range(SEQRESrange1,SEQRESrange2):
            ATOMpos = SEQRESpos-ATOMseqgaplen
            res_no = d_res_nos_ATOM[ATOMpos]['res_no']
            iCode = d_res_nos_ATOM[ATOMpos]['iCode']
            res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
            res_symbol = self.res_name2res_symbol(res_name)

            ## append gap between ATOMseq and SEQRESseq
            if seq+res_symbol != SEQRESseq[:SEQRESpos+1]:
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}
                ATOMseqgaplen += 1
                try:
                    seq += SEQRESseq[SEQRESpos]
                except:
                    print pdb, chain, res_no, iCode, res_name, SEQRESpos, SEQRESrange1, SEQRESrange2
                    print SEQRESseq
                    print self.pdb1, self.pdb2
                    stop
            ## append if no gap
            else:
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':res_no,'iCode':iCode}
                seq += res_symbol

        ## append C-terminal gap between ATOMseq and SEQRESseq
        if len(seq) != len(SEQRESseq):
            SEQRESrange1 = SEQRESpos+1
            SEQRESrange2 = len(SEQRESseq)
            seq, d_res_nos_SEQRES = self.append_missing_residues_to_sequence(seq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq)

        if seq != SEQRESseq:
            print pdb
            print seq
            print SEQRESseq
            notexpected

        return seq, d_res_nos_SEQRES

    def identify_missing_terminal_residues(self, d_pdb, pdb, chain, ATOMseq, SEQRESseq):

        for i in range(1,len(ATOMseq)+1):

            ## index the first occurence of the ATOMseq in the SEQRESseq
            try:
                index = SEQRESseq.index(ATOMseq[:i])
                ## index the next occurence of the ATOMseq in the SEQRESseq
                try:
                    SEQRESseq[index+1:].index(ATOMseq[:i])
                ## break if only one occurence of the ATOMseq in the SEQRESseq
                except:
                    break
            ## break if ATOMseq not occuring in the SEQRESseq
            except:
                print pdb, chain
                print ATOMseq[:i]
                print ATOMseq
                print SEQRESseq
                notexpected
                break

        return index


    def dcoordinates2lcoordinates(self, d_pdb, d_seq, pdb1, pdb2, chain1, chain2):

        d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(d_pdb, d_seq, pdb1, pdb2, chain1, chain2)
        coordinates1, coordinates2, rescount, lines1, lines2 = self.ATOMrecords2coordinates(
            d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2)

        return coordinates1, coordinates2, rescount


    def alignATOMseq(self, d_pdb, d_seq, pdb1, pdb2, chain1, chain2):

        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        sys.path.append('../../EAT_DB/')
        import sequence_alignment

        ##
        ## identify missing residues not mentioned in the REMARK465 records
        ## by alignment of SEQRESseq and ATOMseq
        ## and add gaps to the ATOMseq
        ##
        ## use the nontransformed chain IDs for sequence alignment
        d_ATOM_seqs = {
            pdb1:{'chain':chain1[0]},
            pdb2:{'chain':chain2[0]},
            }
        for pdb in d_ATOM_seqs.keys():
            chain = d_ATOM_seqs[pdb]['chain']
            SEQRESseq = d_seq[pdb]['chains'][chain]['seq']
            ATOMseq,d_res_nos = self.ATOM2seq(d_pdb, pdb, chain, SEQRESseq)
            ## find missing residues, Nterminal; align Nterminal SEQRESseq and ATOMseq
            ATOMseqindentation = self.identify_missing_terminal_residues(d_pdb, pdb, chain, ATOMseq, SEQRESseq)
            ## find missing residues, Cterminal or nonterminal; align SEQRESseq and ATOMseq
            ATOMseq,d_res_nos = self.identify_missing_nonterminal_residues(d_pdb, pdb, chain, SEQRESseq, ATOMseqindentation, d_res_nos)
            d_ATOM_seqs[pdb]['ATOMseq'] = ATOMseq
            d_ATOM_seqs[pdb]['d_res_nos'] = d_res_nos
        ATOMseq1 = d_ATOM_seqs[pdb1]['ATOMseq']
        ATOMseq2 = d_ATOM_seqs[pdb2]['ATOMseq']
        d_res_nos1 = d_ATOM_seqs[pdb1]['d_res_nos']
        d_res_nos2 = d_ATOM_seqs[pdb2]['d_res_nos']

        ##
        ## remove terminal residues from the ATOMseq
        ##
        if ATOMseq1 != ATOMseq2:

            instance = sequence_alignment.NW(ATOMseq1,ATOMseq2)
            s1,s2 = ATOMs1,ATOMs2 = instance.Align(verbose=False)[:2]

            l1 = len(s1)-len(s1.lstrip('-'))
            l2 = len(s2)-len(s2.lstrip('-'))
            r1 = len(s1)-len(s1.rstrip('-'))
            r2 = len(s2)-len(s2.rstrip('-'))
            if r2 == 0:
                ATOMseq1 = ATOMseq1[l2:]
            else:
                ATOMseq1 = ATOMseq1[l2:-r2]
            if r1 == 0:
                ATOMseq2 = ATOMseq2[l1:]
            else:
                ATOMseq2 = ATOMseq2[l1:-r1]

        else:

            l1 = 0
            l2 = 0

        return d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2


    def ATOMrecords2coordinates(
        self, d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2,
        rmsd=None, tv1=None, rm=None, tv2=None):

        import Numeric, math

        rescount = 0
        if len(ATOMseq1) != len(ATOMseq2):
            print pdb1, pdb2, chain1, chain2
            print ATOMseq1
            print ATOMseq2
            print l1, l2
            notexpected

        coordinates1 = []
        coordinates2 = []
        lines1 = []
        lines2 = []

        for SEQRESpos1 in range(l2,len(ATOMseq1)+l2):

            SEQRESpos2 = SEQRESpos1+l1-l2
            res_no1 = d_res_nos1[SEQRESpos1]['res_no']
            res_no2 = d_res_nos2[SEQRESpos2]['res_no']
            if res_no1 == '-' or res_no2 == '-':
                continue
            iCode1 = d_res_nos1[SEQRESpos1]['iCode']
            iCode2 = d_res_nos2[SEQRESpos2]['iCode']

            if 'REMARK' in d_pdb[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1].keys():
                continue
            if 'REMARK' in d_pdb[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2].keys():
                continue

            d_resname = {
                pdb1:{'chain':chain1,'res_no':res_no1,'iCode':iCode1},
                pdb2:{'chain':chain2,'res_no':res_no2,'iCode':iCode2},
                }
            for pdb in d_resname.keys():
                chain = d_resname[pdb]['chain']
                res_no = d_resname[pdb]['res_no']
                iCode = d_resname[pdb]['iCode']
                res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                d_atoms = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']
                d_resname[pdb]['res_name'] = res_name
                d_resname[pdb]['d_atoms'] = d_atoms
            res_name1 = d_resname[pdb1]['res_name']
            res_name2 = d_resname[pdb2]['res_name']
            d_atoms1 = d_resname[pdb1]['d_atoms']
            d_atoms2 = d_resname[pdb2]['d_atoms']

            if res_name1 in self.d_modres.keys():
                res_name1 = self.d_modres[res_name1]
            if res_name2 in self.d_modres.keys():
                res_name2 = self.d_modres[res_name2]
##            if res_name1 != res_name2:
##                print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2
##                print l1, l2
##                print 'SEQRES'
##                print ATOMseq1
##                print ATOMseq2
##                stop
##                fd = open('different_resnames.txt','a')
##                fd.write('%s %s %s %s %s %s %s %s\n' %(pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2))
##                fd.close()

            mutation = False
            if res_name1 != res_name2:
                mutation = True

            line470 = 'REMARK 470     %s %s%4i    ' %(res_name2, chain, res_no2)
            rescount += 1
            if rmsd:
                SS = []
            for atom_name in d_atoms1.keys():
                if mutation == True and atom_name not in ['N','CA','C','O']:
                    continue
                if atom_name not in d_atoms2.keys():
                    if atom_name[0] != 'H' and atom_name[:2] not in ['1H','2H','3H'] and atom_name != 'OXT':
                        line470 += '%s' %(atom_name.ljust(5))
                    continue
                if 'REMARK' in d_atoms1[atom_name].keys():
                    continue
                if 'REMARK' in d_atoms2[atom_name].keys():
                    continue
                ## append coordinates to list of coordinates
                coordinate1 = d_atoms1[atom_name]['coordinate']
                coordinate2 = d_atoms2[atom_name]['coordinate']
                coordinates1 += [coordinate1]
                coordinates2 += [coordinate2]
                if rmsd:
                    coordinate2 = Numeric.matrixmultiply(rm, coordinate2-tv1)+tv2
                    SS += [sum((coordinate2-coordinate1)**2)]

            if rmsd:
                try:
                    RMSD = math.sqrt(sum(SS)/len(SS))
                except:
                    print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2
                    stop
                for atom_name in d_atoms1.keys():
                    if atom_name not in d_atoms2.keys():
                        continue
                    if 'REMARK' in d_atoms1[atom_name].keys():
                        continue
                    if 'REMARK' in d_atoms2[atom_name].keys():
                        continue
                    coordinate1 = d_atoms1[atom_name]['coordinate']
                    coordinate2 = d_atoms2[atom_name]['coordinate']
                    occupancy = bfactor = RMSD/rmsd
                    line1 = self.coordinates2ATOMline(res_name1, chain1[0], res_no1, coordinate1, iCode1, bfactor, atom_name)
                    line2 = self.coordinates2ATOMline(res_name2, chain2[0], res_no2, coordinate2, iCode2, bfactor, atom_name)
                    lines1 += line1
                    lines2 += line1
            if line470 != 'REMARK 470     %s %s%4i    ' %(res_name2, chain, res_no2):

                line470 = "                '"+line470+"', ##%s.pdb %s.pdb\n" %(pdb2, pdb1)
##                fd = open('missingatoms.txt','a')
##                fd.write(line470)
##                fd.close()

        return coordinates1, coordinates2, rescount, lines1, lines2


    def coordinates2ATOMline(self, res_name, chain, res_no, coordinate, iCode, bfactor, atom_name):

        occupancy = bfactor
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        atom_no = 1
        altloc = ''
        charge = ''
        if 'H' in atom_name:
            element = 'H'
        else:
            element = atom_name[0]
        line = [
            '%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n'
            %('ATOM'.ljust(6), atom_no, atom_name.ljust(4), altloc, res_name.ljust(3), chain, res_no, iCode, x, y, z, occupancy, bfactor, element.rjust(2), charge.rjust(2))
            ]

        return line


    def permutation(self, group):

        count_group = len(group)
        count_combinations = self.faculty(count_group)
        combinations = []

        poplist = []

        if count_group > 6:
            stopandfindanothersolution

        for i in range(count_combinations):
            combinations += [[]]
            poplist += [list(group)]

        i = 0
        ## loop over n
        for i in range(count_group,0,-1):
            i_fac = self.faculty(i-1)
            j = 0
            ## loop over n!
            while j < count_combinations:
                ## loop over i ... to loop over i!=i*(i-1)
                for k in range(i):
                    ## loop over (i-1)! ... to loop over i!=i*(i-1)!
                    for l in range(i_fac):
                        combinations[j].append(poplist[j].pop(k))
                        j += 1

        return combinations


    def parse_coordinates(self, s_pdb, d_pdb):

        ## read lines
##        fd = open('%s%s.pdb' %(self.pdbpath, s_pdb.upper()),'r')
        fd = open('%s%s/pdb%s.ent' %(self.pdbpath, s_pdb.lower()[1:3], s_pdb.lower()),'r')
        lines = fd.readlines()
        fd.close()
        if s_pdb not in d_pdb.keys():
            d_coordinates = self.parse_pdbcoordinatesection(lines, s_pdb)
            d_pdb[s_pdb] = d_coordinates

        return d_pdb


    def parse_sequences(self):

        import time

        print 'parsing data from noncoordinate sections'
        
        d_seq = {}
        d_chains_intrapdb_sequence_identical = {}

        t1 = time.clock()
        for i in range(self.pdbcount):

            self.pdb = s_pdb = self.l_pdbs[i]

            t2 = time.clock()
            if t2-t1 > self.time_status_update:
                print 'parsing %s (%s/%s)' %(self.l_pdbs[i-1], i, self.pdbcount)
                t1 = t2

            ## read lines
##            fd = open('%s%s.pdb' %(self.pdbpath, s_pdb.upper()),'r')
            fd = open('%s%s/pdb%s.ent' %(self.pdbpath, s_pdb.lower()[1:3], s_pdb.lower()),'r')
            lines = fd.readlines()
            fd.close()

            ## parse data prior to the coordinate section
            d_noncoordinates = self.parse_pdbnoncoordinatesections(lines, s_pdb)
            d_seq[s_pdb] = d_noncoordinates
            d_chains_intrapdb_sequence_identical = self.identify_identical_chains_from_sequence_intra(
                d_seq,s_pdb,d_chains_intrapdb_sequence_identical,
                )
            

        return d_seq, d_chains_intrapdb_sequence_identical


    def parse_pdbnoncoordinatesections(self, lines, s_pdb):

        s_pdb = s_pdb.lower()
        ## import sequence alignment
        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        sys.path.append('../../EAT_DB/')
        import sequence_alignment

        ## parser written on the assumption that SEQRES is mandatory if ATOM records exist

        d_seq = {}
        d_conect = {}
        l_hetatms = []
        parse_atoms = False

        ## insertion chain, res_no, res_name
        prev_chain = ''
        d_insertions = {}
        biounit = 'N/A'

        for i in range(len(lines)):
            line = lines[i]

            record = line[:6].strip()

            if record == 'ATOM': ## section 9
                continue
##                if parse_atoms == False:
##                    if 'REMARK350' in d_seq.keys() and 'HET' in d_seq.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = self.parse_atom_no_range(d_conect, 'atom', atom_no)

            elif record == 'HETATM': ## section 9
                continue
##                if parse_atoms == False:
##                    if 'REMARK350' in d_seq.keys() and 'HET' in d_seq.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = self.parse_atom_no_range(d_conect, 'atom', atom_no)

            elif record == 'REMARK': ## section 2
                d_seq = self.parse_recordREMARK(d_seq, line, i, lines)

            elif record == 'SEQRES': ## section 3
                d_seq = self.parse_recordSEQRES(line, d_seq)

            elif record == 'HET': ## section 4
                hetID = line[7:10].strip()
                ## continue if water
                if hetID in ['H20','HOH','D2O','DOD']: ## D2O in 2JAJ
                    continue
                chain = line[12]
                if 'HET' not in d_seq.keys():
                    d_seq['HET'] = {}
                if chain not in d_seq['HET'].keys():
                    d_seq['HET'][chain] = set()
                d_seq['HET'][chain] |= set([hetID])

            elif record == 'MODEL':
                d_seq['MODEL'] = True

            elif record == 'TITLE': ## section 2
                if not 'TITLE' in d_seq.keys():
                    d_seq['TITLE'] = line[10:].strip()
                else:
                    if d_seq['TITLE'][-1] == '-':
                        d_seq['TITLE'] += line[10:].strip()
                    else:
                        d_seq['TITLE'] += ' '+line[10:].strip()

## #-MER / @-MER
## MULTIMER
## CHAINS A-C / A,B,C
## HOMO- / HETERO-
## OLIGOMER OF # SUBUNITS
##            elif record == 'COMPND': ## section 2
##                if 'BIOLOGICAL_UNIT' in line or 'BIOLOGICAL UNIT' in line:
##                    print lines[i-1], lines[i], lines[i+1]
##                    for s_biounit in self.l_biounits:
##                        if s_biounit in line:
##                            biounit = s_biounit
##                    if biounit == 'N/A':
##                        for s_biounit in self.d_biounits.keys():
##                            if s_biounit in line:
##                                biounit = s_biounit
##                    print self.pdb, biounit

            elif record == 'HEADER':
                d_seq['HEADER'] = line[10:50].strip()

            elif record == 'SPRSDE': ## section 2
                sIDcodes = line[31:70].lower().split()
                for sIDcode in sIDcodes:
                    if sIDcode not in d_seq.keys():
                        d_seq[sIDcode] = {}
                    d_seq[sIDcode]['SPRSDE'] = True

            elif record == 'EXPDTA': ## section 2
                methods = line[10:]
                if methods[:3] == 'NMR':
                    d_seq['EXPDTA'] = 'NMR'
                elif methods.strip() == 'CRYO-ELECTRON MICROSCOPY':
                    d_seq['EXPDTA'] = 'CRYO-ELECTRON MICROSCOPY'
                else:
                    d_seq['EXPDTA'] = 'N/A'

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                d_seq['CRYST1'] = spacegroup

        ## strings
        for key in ['TITLE','EXPDTA','REMARK2']:
            if key not in d_seq.keys():
                d_seq[key] = 'N/A'
        ## dics
        for key in ['chains','HET']:
            if key not in d_seq.keys():
                d_seq[key] = {}
        ## lists
        for key in ['REMARK525']:
            if key not in d_seq.keys():
                d_seq[key] = []

        proteinchains = []
        peptidechains = []
        nucleotidechains = []
        saccharidechains = []
        for chain in d_seq['chains'].keys():
            if d_seq['chains'][chain]['type'] == 'peptide':
                peptidechains += chain
                if len(d_seq['chains'][chain]['seq']) > self.min_len_chain:
                    proteinchains += chain
            elif d_seq['chains'][chain]['type'] == 'nucleotide':
                nucleotidechains += chain
            elif d_seq['chains'][chain]['type'] == 'saccharide':
                saccharidechains += chain

        d_seq['proteinchains'] = proteinchains

##        if self.pdbcount > 10000 and s_pdb not in ['1ady','1bhj']:
##            chains = d_seq['chains'].keys()
##            if 'REMARK350' not in d_seq.keys() and len(d_seq['chains'].keys()) > 1 and d_seq['EXPDTA'] != 'NMR' and len(saccharidechains) != len(d_seq['chains'].keys()):
##                ## biounit not specified as text
##                if biounit == 'N/A':
##                    fd = open('unknownbiounit.txt','a')
##                    fd.write('%s %s %s %s %s\n' %(s_pdb, len(proteinchains), len(peptidechains), len(nucleotidechains), len(saccharidechains)))
##                    fd.close()
##                ## unequal number of proteinchains and size of biounit
##                elif self.d_biounits[biounit] != len(chains):
##                    if biounit == 'MONOMER':
##                        d_seq['REMARK350'] = {}
##                        for i in range(len(proteinchains)):
##                            chain = proteinchains[i]
##                            d_seq['REMARK350'][i+1] = {'chains': [chain]}
##                    elif len(chains) % self.d_biounits[biounit] == 0:
#### add all combinations of chains to remark350 transformations ?! redudant hits if done for twin pdb as well...
##                        print biounit, self.d_biounits[biounit]
##                        print s_pdb
##                        print d_seq.keys()
##                        print d_seq['chains'].keys()
##                        print proteinchains
##                        print d_seq['HET']
##                        stop3
##                    else:
##                        print s_pdb, biounit, chains, proteinchains
##                        print d_seq['HET']
##                        stop3b
##                elif self.d_biounits[biounit] == len(chains):
##                    d_seq['REMARK350'] = {1:{'chains': [chains]}}
##                else:
##                    print s_pdb, biounit, proteinchains, chains
##                    stop5
####                d_seq['biounit'] = 'small' ## small opposed to monomer if multiple different chains
#### count number of similar chains by comparing SEQRESseqs

        return d_seq


    def parse_atom_no_range(self, d_conect, record, atom_no):
        
        if not record in d_conect.keys():
            d_conect[record] = [[atom_no,atom_no]]
        elif d_conect[record][-1][1] == atom_no-1:
            d_conect[record][-1][1] = atom_no
        else:
            d_conect[record] += [[atom_no,atom_no]]

        return d_conect


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):

        import Numeric

        translationvector = Numeric.array(
            [
                float(transformationmatrix[0][3]),
                float(transformationmatrix[1][3]),
                float(transformationmatrix[2][3]),
                ]
            )

        rotationmatrix = Numeric.array(
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


    def parse_recordSEQRES(self, line, d_seq):

        chain = line[11]

        if 'chains' not in d_seq:
            d_seq['chains'] = {}
        if chain not in d_seq['chains'].keys():
            d_seq['chains'][chain] = {}
        if not 'type' in d_seq['chains'][chain].keys():
            d_seq['chains'][chain]['type'] = 'unknown'

        residues = line[19:70].split()

        for i in range(len(residues)):
            residue = residues[i]
            if residue in self.d_res.keys():
                if d_seq['chains'][chain]['type'] == 'unknown':
                    d_seq['chains'][chain]['type'] = 'peptide'
                residues[i] = self.d_res[residue]
            elif residue in ['C','A','T','G','DC','DA','DT','DG']:
                if d_seq['chains'][chain]['type'] == 'unknown':
                    d_seq['chains'][chain]['type'] = 'nucleotide'
                residues[i] = residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_seq['chains'][chain]['type'] == 'unknown':
                    d_seq['chains'][chain]['type'] = 'saccharide'
                residues[i] = residue
            else:
                residues[i] = 'X'

        if 'seq' not in d_seq['chains'][chain].keys():
            d_seq['chains'][chain]['seq'] = ''
        d_seq['chains'][chain]['seq'] += ''.join(residues)

        return d_seq


    def parse_pdbcoordinatesection(self, lines, s_pdb):

##        print 'parsing coordinates of %s' %(s_pdb)

        s_pdb = s_pdb.lower()

        d_coordinates = {
            'chains':{},
            'HET':set(),
            }

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()
    
            if record == 'ATOM':
                d_coordinates = self.parse_recordATOM(line, d_coordinates, lines, i)[0]

            elif record == 'HETATM':
                res_name = line[17:20].strip()
                if res_name in self.d_modres.keys():
                    d_coordinates = self.parse_recordATOM(line, d_coordinates, lines, i)[0]

            elif record == 'REMARK':
                remark = int(line[6:10])
                if remark == 465:
                    d_coordinates = self.parse_recordREMARK465(line, d_coordinates, lines, i)
                elif remark == 470:
                    d_coordinates = self.parse_recordREMARK470(line, d_coordinates, lines, i)

            elif record == 'HET':
                HETID = line[7:10].strip()
                d_coordinates['HET'] |= set([HETID])

            elif record == 'MODEL':
                model = int(line.split()[1])

        return d_coordinates


    def parse_recordREMARK(self, d_seq, line, i, lines):

        remark = int(line[6:10])

        if remark == 200:

            if line[12:23].strip().upper() in ['TEMPERATURE','PH']:
                experimentaldetail_key = line[12:23].strip()
                experimentaldetail_value = line[44:].strip()
                if 'REMARK200' not in d_seq.keys():
                    d_seq['REMARK200'] = {}
                if experimentaldetail_key not in d_seq['REMARK200'].keys():
                    d_seq['REMARK200'][experimentaldetail_key] = experimentaldetail_value

        elif remark == 350:

            ## biological units
            ## (e.g. 2bq0.pdb, 1thj.pdb, 1m4x.pdb, 1d3i.pdb, 1qgc.pdb, 1rhi.pdb, 1rbo.pdb, 2g8g.pdb, 1h84.pdb)

            d_seq = self.parse_recordREMARK350(d_seq, i, lines)

        elif remark == 525: ## water association

            if line[11:].strip() == 'PROTEIN CHAIN  SOLVENT CHAIN':
                for j in range(i+1,len(lines)):
                    if lines[j][11:].strip() == '':
                        break
                    if lines[j][:10] != 'REMARK 525':
                        break
                    else:
                        solventchain = lines[j][11:].split()[1]
                        proteinchain = lines[j][11:].split()[0]
                        if solventchain == proteinchain: ## e.g. 2bq0.pdb
                            continue
                        if 'REMARK525' not in d_seq.keys():
                            d_seq['REMARK525'] = []
                        d_seq['REMARK525'] += [solventchain]

        elif remark == 2: ## resolution
            try:
                resolution = float(line[22:27])
            except:
                resolution = 'N/A'
            d_seq['REMARK2'] = resolution

        return d_seq


    def parse_recordREMARK350(self, d_seq, i, lines):

        line = lines[i]

        if 'REMARK350' not in d_seq.keys():
            d_seq['REMARK350'] = {}

        if line[11:23] == 'BIOMOLECULE:':
            biomolecules = line[23:80].replace(' ','').split(',')
            d_seq = self.loop_and_identify_chains_and_matrices(i, lines, d_seq, biomolecules)

        return d_seq


    def parse_REMARK350_chains(self, line_chains):

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


    def loop_and_identify_biomolecules(self, i, lines):

        for j in range(i-1,-1,-1):

            if lines[j][:10] != 'REMARK 350':
                return False
            if lines[j][11:23] == 'BIOMOLECULE:':
                return True


    def loop_and_identify_chains_and_matrices(self, i, lines, d_seq, biomolecules):

        chains = set()

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

                if self.pdb == '1m4x': ## 1m4x.pdb
                    index = lines[j][24:].index('.')
                    matrixno = int(lines[j][24+index-7:24+index-2])
                else:
                    matrixno = int(lines[j][19:24])
                ## parse transformation matrix
                matrixrow1 = lines[j-2][24:].split()
                matrixrow2 = lines[j-1][24:].split()
                matrixrow3 = lines[j-0][24:].split()
                matrixrows = [matrixrow1,matrixrow2,matrixrow3,]
##                ## find out whether transformation matrix yields a transformation
##                transformation = False
##                for k in range(3):
##                    ## add a zero translation vector if a translation vector is not given
##                    if len(matrixrows[k]) == 3:
##                        matrixrows[k] += [0.]
##                    if float(matrixrows[k][k]) == 1. and float(matrixrows[k][3]) == 0.:
##                        continue
##                    else:
##                        transformation = True

                ## append transformation matrix to dictionary
                for biomolecule in biomolecules:

                    biomolecule = int(biomolecule)

                    ## biomolecule
                    if biomolecule not in d_seq['REMARK350'].keys():
                        d_seq['REMARK350'][biomolecule] = {}

                    ## biomolecule > matrices
                    if 'matrices' not in d_seq['REMARK350'][biomolecule].keys():
                        d_seq['REMARK350'][biomolecule]['matrices'] = {}
                    ## matrices > matrixno > matrix
                    d_seq['REMARK350'][biomolecule]['matrices'][matrixno] = matrixrows

                    ## biomolecule > chains
                    if 'chains' not in d_seq['REMARK350'][biomolecule].keys():
                        d_seq['REMARK350'][biomolecule]['chains'] = {}
                    for chain in chains:
                        ## chains > chain
                        if chain not in d_seq['REMARK350'][biomolecule]['chains'].keys():
                            d_seq['REMARK350'][biomolecule]['chains'][chain] = set()
                        d_seq['REMARK350'][biomolecule]['chains'][chain] |= set([matrixno])

        return d_seq


    def parse_recordREMARK465(self, line, d_pdb, lines, i):

        ## missing residues

        if line[10:].strip() == 'M RES C SSSEQI':

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 465':
                    break

                try:
                    model = int(lines[j][12:14])
                except:
                    model = 'N/A'
                res_name = lines[j][15:18]
                chain = lines[j][19]
                res_no = int(lines[j][22:26])
                res_name = lines[j][15:18]
                iCode = lines[j][26]

                if not chain in d_pdb['chains'].keys():
                    d_pdb['chains'][chain] = {}
                if not 'residues' in d_pdb['chains'][chain].keys():
                    d_pdb['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['chains'][chain]['residues'].keys():
                    d_pdb['chains'][chain]['residues'][res_no] = {}

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}


                ## iCode > res_name
                if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name
                ## check that res_name is correct (e.g. 2fes:L:1)
                elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name and atom_altloc == ' ': ## 1fh2:A:30 atom_altloc
                    ## change the iCode
                    iCode_max = max(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
                    iCode_max = self.s_alphabet[self.s_alphabet.index(iCode_max)+1]
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_max] = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'].index(iCode)] = iCode_max

                    ## d_iCodes > iCode
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}


                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if not iCode in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## iCode > REMARK
                if not 'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] = True

        elif line[10:].strip().split() == ['M','RES','C','SSSEQI']:
            notexpected

        return d_pdb


    def parse_recordREMARK470(self, line, d_pdb, lines, i):

        ## missing atoms

        ## the latter equation is only to acommodate for 1fvk.pdb
        if line[10:].strip() == 'M RES CSSEQI  ATOMS' or line[10:].strip() == 'M RES C SEQI  ATOMS':

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 470':
                    break

                ## model M
                try:
                    model = int(lines[j][11:13])
                except:
                    model = 'N/A'

                ## res_name RES
                res_name = lines[j][15:18]
                if res_name not in self.d_res.keys():
                    continue

                ## chain C
                chain = lines[j][19]

                ## res_no SSEQ
                try:
                    res_no = int(lines[j][20:24])
                except:
                    res_no = lines[j][20:24].strip()

                ## iCode I
                iCode = lines[j][24]

                ## atoms ATOMS
                atoms = lines[j][25:].split()

                ##
                ## write to dictionary
                ##
                if not chain in d_pdb['chains'].keys():
                    d_pdb['chains'][chain] = {}
                if not 'residues' in d_pdb['chains'][chain].keys():
                    d_pdb['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['chains'][chain]['residues'].keys():
                    d_pdb['chains'][chain]['residues'][res_no] = {}

                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if not iCode in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

                ## iCode > atoms
                if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                ## atoms > atom_name > coordinate
                for atom_name in atoms:
                    if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['REMARK'] = 470

                ## iCode > res_name
                if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name

        elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
            print self.pdb1
            print self.pdb2
            notexpected

        return d_pdb


    def parse_recordATOM(self, line, d_pdb, lines, i):

        import Numeric

        atom_no = int(line[6:11])
        atom_name = line[12:16].strip()
        atom_altloc = line[16]
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        atom_x = float(line[30:38])
        atom_y = float(line[38:46])
        atom_z = float(line[46:54])
        coordinate = Numeric.array([atom_x, atom_y, atom_z])
        if not 'chains' in d_pdb.keys():
            d_pdb['chains'] = {}
        if not chain in d_pdb['chains'].keys():
            d_pdb['chains'][chain] = {}
        if not 'residues' in d_pdb['chains'][chain].keys():
            d_pdb['chains'][chain]['residues'] = {}

        ## res_no
        if not res_no in d_pdb['chains'][chain]['residues'].keys():
            d_pdb['chains'][chain]['residues'][res_no] = {}

        ## res_no > d_iCodes
        if not 'd_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
        ## d_iCodes > iCode
        if not iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes']:
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

        ## iCode > res_name
        if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name
        ## check that res_name is correct (e.g. 2fes:L:1)
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name and atom_altloc == ' ': ## 1fh2:A:30 atom_altloc
            ## change the iCode
            iCode_max = max(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
            iCode_max = self.s_alphabet[self.s_alphabet.index(iCode_max)+1]
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_max] = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'].index(iCode)] = iCode_max

            ## d_iCodes > iCode
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

        ## res_no > l_iCodes
        if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
        ## l_iCodes > iCode
        if not iCode in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

        ## iCode > atoms
        if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
        ## atoms > atom_name > coordinate
        if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}

        return d_pdb, {'chain':chain,'atom_no':atom_no}


    def faculty(self, n):

        fac = 1
        for i in range(1,n+1):
            fac *= i

        return fac


    def __init__(self):

        import os

        self.d_res = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            'MSE':'M',
            }
        
        ## HETATM res_names for which coordinates are parsed
        self.d_modres = {
            'MSE':'MET', ## selenomethionine
            }

        ## HET hetIDs which are ignored
        self.l_modres = [
            'MSE',
            'ACE', ## e.g. 1spd.pdb (vs 2c9v.pdb)
            ]
            
        self.d_res3 = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

        self.time_status_update = 1

        ## list of biounit strings in which other biounit string can be found
        self.l_biounits = [
            'DIMER OF DIMERS', ## DIMER
            'HEXADECAMER', ## DECAMER
            'DODECAMER', ## DECAMER
            ]
        self.d_biounits = {
            'MONOMER':1,
            'DIMER':2,
            'TRIMER':3,
            'TETRAMER':4,'DIMER OF DIMERS':4, ## DIMER
            'PENTAMER':5,
            'HEXAMER':6,
            'HEPTAMER':7,
            'OCTAMER':8,
            'DECAMER':10,
            'DODECAMER':12 ,'12-MER':12, ## DECAMER
            'HEXADECAMER':16, ## DECAMER
            'ICOSAHEDRAL':60,
            }

        ## http://en.wikipedia.org/wiki/Crystal_structure

        ## Hermann-Mauguin symbols *and* disallowed abbrevations *and* errors (e.g. A 2 space group of 1mbs.pdb)
        ## sorted from low to high symmetry

        ## errornous space groups assigned to a crystal system based on information about the unit cell
        ## from only one representative structure with the space group in question

        ## P 21 21 21, P 1 21 1, C 1 2 1 most common for proteins...

        d_spacegroups = {
            ## alpha,beta,gamma != 90
            'TRICLINIC':[
                'P 1',
                'P 1-', 'P1', ## neither HM symbols nor abbrevations
                'A 1', ## 1lks.pdb
                ],
            ## alpha != 90, beta,gamma==90
            'MONOCLINIC':[
                'P 1 21 1','C 1 2 1','P 1 2 1',
                'P 21', ## P 1 21 1 abbreviations
                'C 2', 'C 21', 'C 1 21 1', ## C 1 2 1 abbreviations
                'P 2', ## P 1 2 1 abbreviations
                'B 2', 'I 1 2 1', 'P 1 1 21', 'I 21', 'I 1 21 1', ## neither HM symbols nor abbrevations
                ],
            ## a != b != c (alpha,beta,gamma==90)
            'ORTHORHOMBIC':[
                'C 2 2 21','P 21 21 21','P 21 21 2','I 2 2 2','C 2 2 2','I 21 21 21','P 2 2 21','P 2 2 2','F 2 2 2',
                'P 2 21 21', ## neither HM symbols nor abbrevations
                'P 21 2 21',
                'P 21 21 2 A', ## P 21 21 2 error in 1b86.pdb
                'B 2 21 2', ## 1zna.pdb
                'B 1 1 2', ## 1qr6.pdb
                ],
            ## a != c (a == b, alpha,beta,gamma==90)
            'TETRAGONAL':[
                'P 43 21 2','I 41 2 2','I 41','I 4','P 42 21 2','P 41 21 2','I 4 2 2','P 41','P 43','P 4 21 2','P 4','P 4 2 2','P 41 2 2','P 43 2 2','P 42 2 2','P 42',
                ],
            ## RHOMBOHEDRAL (a=b=c, alpha,beta,gamma!=90)
            ## alpha,beta,gamma != 90
            'TRIGONAL':[
                'H 3 2','R 3','P 32 2 1','P 3','P 31 2 1','P 31','P 32 1 2','R 3 2','P 3 2 1','P 32','P 3 1 2','P 31 1 2',
                'H 3', ## R 3 equivalent
                'H 3 2', ## P 3 2 1 equivalent
                ],
            'HEXAGONAL':[
                'P 61 2 2','P 65','P 63','P 65 2 2','P 61','P 62 2 2','P 62','P 64 2 2','P 63 2 2','P 6 2 2','P 6','P 64',
                ],
            'CUBIC':[
                'F 41 3 2','P 21 3','I 4 3 2','I 2 3','P 2 3','P 41 3 2','P 4 3 2','F 4 3 2','P 43 3 2','I 21 3','F 2 3','P 42 3 2','I 41 3 2',
                ],
            'UNKNOWN':[
                'A 2',
                '', ## NMR structure 2ait constains a CRYST1 record
                ]
            }

        self.d_crystalsystems = {}
        for crystalsystem in d_spacegroups.keys():
            for spacegroup in d_spacegroups[crystalsystem]:
                self.d_crystalsystems[spacegroup] = crystalsystem

        self.nontransformationmatrix = [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]

        self.s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        ## info from the PDB Ligand Depot
        ## keys are hetIDs, values are charges (not oxidation states)
        ## chemical formulas are not specified
        self.d_ions = {
            'LI' :+1,
            'NA' :+1,'NAO':+1,'NA2':+1,'NAW':+1,'NA5':+1,'NA6':+1, ## different number of waters coordinated
            'MG' :+2,'MO1':+2,'MO2':+2,'MO3':+2,'MO4':+2,'MO5':+2,'MO6':+2, ## different number of waters coordinated
            'AL' :+3,
            '2HP':-1,'IPS':-2,'PI' :-2,'PO4':-3, ## synonyms and different oxidation states
            'SO4':-2,'SOH':-1,'SUL':-2, ## synonyms and different oxidation states
            'CL' :-1,
            'K'  :+1,'KO4':+1, ## different number of waters coordinated
            'CA' :+2,'OC1':+2,'OC2':+2,'OC3':+2,'OC4':+2,'OC5':+2,'OC6':+2,'OC7':+2,'OC8':+2, ## different number of waters coordinated
            'V'  :+3,'VO4':-3, ## different oxidation states
            'CR' :+3,
            'MN' :+2,'MN3':+3,'MW1':+2,'MW2':+2,'MW3':+2,'O4M':+2,'MN5':+2,'MN6':+2, ## different oxidation states and different number of waters coordinated
            'FE2':+2,'OF1':+2,'2OF':+2,'FE' :+3,'OF3':+3, ## different oxidation states and different number of waters coordinated
            'CO' :+2,'OCL':+2,'OCN':+2,'OCM':+2,'3CO':+3,'CO5':+3,'OCO':+3, ## different oxidation states and different number of waters coordinated
            'NI' :+2,'3NI':+3,'NI1':+2,'NI2':+2,'NI3':+2,'NIK':+2, ## different oxidation states and different number of waters coordinated
            'CU1':+1,'CU' :+2,'1CU':+2, ## different oxidation states and different number of waters coordinated
            'ZN' :+2,'ZN2':+2,'ZN3':+2,'ZNO':+2,'ZO3':+2, ## different number of waters coordinated and different crystal fold axes
            'GA' :+3,
            'ARS': 0,'ART':-3,'AST':-3,'TAS': 0, ## different compounds
            'SE' : 0,'SE4':-2, ## different compounds
            'BR' :-1,
            'KR' : 0,
            }

        ## saccharides returned from a search of the ligand depot for disaccharides and monosaccharides
        self.d_saccharides = {
            ## monosaccharides, aldehydes
            'GLC':['GLC'], ## (alpha)-D-Glucose
            'AGC':['GLC'], ## alpha-D-Glc
            'BGC':['GLC'], ## beta-D-Glc
            'BG6':['GLC'], ## beta-D-Glc-6P
            'G1P':['GLC'], ## alpha-D-Glc-1P
            'G6P':['GLC'], ## alpha-D-Glc-6P
            'G6Q':['GLC'], ## Glc-6P
            'GAL':['GAL'], ## (beta)-D-Galactose
            'GLA':['GAL'], ## alpha-D-Gal
            'GLB':['GAL'], ## beta-D-Gal
            'FUC':['GAL'], ## 6-deoxy-GAL, alpha-L-Fucose
            'FUL':['GAL'], ## 6-deoxy-GAL, beta-L-Fucose
            'BGP':['GAL'], ## beta-Gal-6P
            'MAN':['MAN'], ## alpha-D-Mannose
            'BMA':['MAN'], ## beta-D-Mannose
            'M1P':['MAN'], ## alpha-D-Man-1P
            'M6P':['MAN'], ## alpha-D-Man-6P
            ## monosaccharides, ketones
            'FRU':['FRU'], ## Fructose
            'F6P':['FRU'], ## Fru-6P
            ## dissacharides
            'SUC':['GLC','FRU'], ## GLC-a12-FRC, Sucrose
            'LAT':['GAL','GLC'], ## GAL-b14-GLC, alpha-Lactose
            'LBT':['GAL','GLC'], ## GAL-b14-GLC, beta-Lactose
            'MAL':['GLC','GLC'], ## GLC-a14-GLC, Maltose
            'TRE':['GLC','GLC'], ## GLC-a11a-GLC, Trehalose
            'CBI':['GLC','GLC'], ## GLC-b14-GLC, Cellobiose
            }

        self.maxrmsd = 2.5

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25
        self.max_mutations = 10

        self.pdbpath = '/oxygenase_local/data/pdb/'
        
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

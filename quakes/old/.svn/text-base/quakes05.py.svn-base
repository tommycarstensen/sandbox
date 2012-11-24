#!/software/bin/python
#
#$Id: quakes.py,v 1.21 2007/05/08 10:14:07 tc Exp $
#
# Tommy Carstensen, University College Dublin, 2007

##
## questions of interest

## find correlation (if any) between RMSD of backbone/all atoms for all/neighboring(exponentional sphere radii)/surface residues

## rmsd not mathematically dependent on chain length?! but rmsd physically dep on chain length?

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

## structural alignment of backbone atoms only if comparing chains with point mutations?

## write a new faster seq aln alg

## require identical ions to be present or only multiatom hetero atoms identical?

## ignore MSE upon comparison of hetIDs?

## use HETSYN records to find het synonyms if any...

## find a pdb file with SEQADV records about mutations
## and see if the ATOM records correspond to SEQRES or SEQADV !!

## add dposition year to rmsd.txt...
## add number of coordinates/atoms, residues, chains compared to rmsd.txt!
## add EC number if applicable to rmsd.txt

## do not outcomment "return 2" in python/EAT_DB/sequence_alignment.py/NW_gap

class quakes:

    def main(self):

        pdbs = []
        ## pdbs I find being good examples
        pdbs += [
            '6tim.pdb','1a4f.pdb','4ins.pdb','1aiy.pdb','1td3.pdb','2ktq.pdb','3ktq.pdb','2g34.pdb','2g33.pdb','2nu4.pdb','2nu2.pdb','2hjl.pdb','1tmc.pdb','1zhk.pdb','1jlp.pdb','1g7l.pdb','1vwi.pdb','1vwl.pdb','1vwj.pdb',
            '1kj3.pdb','1s7u.pdb','1n59.pdb','2pol.pdb','1mmi.pdb','2grr.pdb','2grn.pdb','1jy7.pdb','1lfl.pdb','1q8o.pdb','1q8p.pdb',
            '1bks.pdb','2wsy.pdb',
            ]
        ## GroEL-GroES multimers
        pdbs += ['1sjp.pdb','1we3.pdb','1wf4.pdb','1pf9.pdb','1aon.pdb','1grl.pdb','1iok.pdb','1j4z.pdb','1kp8.pdb','1kpo.pdb','1mnf.pdb','1oel.pdb','1pcq.pdb','1pf9.pdb']
        ## T4 lysozyme
        pdbs += ['2lzm.pdb','150l.pdb']
        ## hemoglobin with differences between asymnmetric unit and biological unit
        pdbs += ['2hhb.pdb','1hho.pdb','1hv4.pdb']
        ## identical
        pdbs += ['1cse.pdb', '1csg.pdb', '1bnc.pdb', '1ar9.pdb', '1d7d.pdb', '1bnj.pdb', '1bni.pdb', '1ar5.pdb', '1fnw.pdb', '1fnv.pdb', '1fnu.pdb', '1yhe.pdb', '1ca0.pdb', '1hn4.pdb', '1ca8.pdb', '1h7g.pdb', '1yhr.pdb', '1ez0.pdb', '1h26.pdb', '1f9q.pdb', '2gys.pdb', '1g5p.pdb', '1czg.pdb', '1caw.pdb', '1cav.pdb', '1cau.pdb', '1d7c.pdb', '1i4v.pdb', '1fn8.pdb', '1yh9.pdb', '1c77.pdb', '1djn.pdb', '1c79.pdb', '1c78.pdb', '1gut.pdb', '1gus.pdb', '1dsb.pdb', '1rzj.pdb', '1cl5.pdb', '1a7k.pdb', '1e5y.pdb', '1c7e.pdb', '2g4g.pdb', '1c7f.pdb', '1c08.pdb', '1uiw.pdb', '1dxt.pdb', '3bjl.pdb', '1zgo.pdb', '2bgf.pdb', '1bqq.pdb', '2ad6.pdb', '1sxa.pdb', '1sxb.pdb', '1sxc.pdb', '1ay9.pdb', '1d8a.pdb', '2nip.pdb', '1d8m.pdb', '1f5z.pdb', '1qty.pdb', '1v6l.pdb', '1brs.pdb', '1aym.pdb', '1ayn.pdb', '1ces.pdb', '1o7q.pdb', '1ayy.pdb', '1sx4.pdb', '1alj.pdb', '1ali.pdb', '1fx9.pdb', '1b2k.pdb', '1bi2.pdb', '1b2x.pdb', '1fmc.pdb', '1gho.pdb', '1h25.pdb', '1h24.pdb', '1h27.pdb', '3gct.pdb', '1bij.pdb', '1ggt.pdb', '2f74.pdb', '2dlf.pdb', '1ggx.pdb', '1c48.pdb', '1ggc.pdb', '1ggb.pdb', '1fqj.pdb', '1fqk.pdb', '1tta.pdb', '1gh7.pdb', '1h2u.pdb', '1h2t.pdb', '1asj.pdb', '1g6y.pdb', '1qia.pdb', '1qic.pdb', '1g6w.pdb', '1ffy.pdb', '1azb.pdb', '1azc.pdb', '2gaw.pdb', '1os7.pdb', '1cbj.pdb', '1gq9.pdb', '1qpz.pdb', '1cbl.pdb', '1cbm.pdb', '2hco.pdb', '1g69.pdb', '1n8o.pdb', '1g67.pdb', '1gzx.pdb', '1kd2.pdb', '1as8.pdb', '1bz1.pdb', '1as7.pdb', '1as6.pdb', '1lft.pdb', '1lfv.pdb', '1lfq.pdb', '1d7b.pdb', '1lfy.pdb', '1lfz.pdb', '1az0.pdb', '1az3.pdb', '1gqw.pdb', '1lfl.pdb', '1f13.pdb', '1gdn.pdb', '2fmp.pdb', '1bdi.pdb', '1gc0.pdb', '1c8v.pdb', '1gc2.pdb', '1bm3.pdb', '1bm0.pdb', '1pf9.pdb', '1b6z.pdb', '1bdt.pdb', '1ahi.pdb', '1a3n.pdb', '1fm0.pdb', '1efc.pdb', '1par.pdb', '1nrr.pdb', '1nrq.pdb', '1nrp.pdb', '1bmz.pdb', '1n5a.pdb', '1x1u.pdb', '1dv1.pdb', '1gct.pdb', '1hcj.pdb', '1rva.pdb', '1hco.pdb', '1rve.pdb', '1b27.pdb', '1qg6.pdb', '1fma.pdb', '1cm1.pdb', '1kxz.pdb', '1bz0.pdb', '1c14.pdb', '1afu.pdb', '1opg.pdb', '1c0m.pdb', '1dlf.pdb', '1fb2.pdb', '1fb0.pdb', '1fb6.pdb', '1c1a.pdb', '1xz2.pdb', '1bzw.pdb', '1qu2.pdb', '1gvj.pdb', '2f4y.pdb', '1aar.pdb', '1iv5.pdb', '1bze.pdb', '1bzd.pdb', '1gmr.pdb', '1gmq.pdb', '1gmp.pdb', '1jh6.pdb', '1hut.pdb', '1cqe.pdb', '1a0u.pdb', '1cqn.pdb', '1cqm.pdb', '1a0z.pdb', '1a0f.pdb', '1cqr.pdb', '1fwd.pdb', '1fwb.pdb', '1bhj.pdb', '1fph.pdb', '1d2c.pdb', '1fpj.pdb', '1hgb.pdb', '1hgc.pdb', '1w7s.pdb', '1hga.pdb', '1fpb.pdb', '1fpe.pdb', '1fpd.pdb', '1fpg.pdb', '1fpf.pdb', '1dg1.pdb', '1f41.pdb', '1g28.pdb', '1hbb.pdb', '1f49.pdb', '1rq3.pdb', '1a00.pdb', '1fia.pdb', '3fis.pdb', '1xpt.pdb', '1atz.pdb', '1qjz.pdb', '1cq9.pdb', '1dgc.pdb', '1f4a.pdb', '1n2a.pdb', '1f4h.pdb', '1f4j.pdb', '1dgr.pdb', '1dgw.pdb', '1k3o.pdb', '1buv.pdb', '1ffn.pdb', '1cjr.pdb', '1cjq.pdb', '1grc.pdb', '2ad7.pdb', '1abi.pdb', '1abj.pdb', '2ad8.pdb', '2tmd.pdb', '1d9q.pdb', '1c5x.pdb', '1rhp.pdb', '1c5w.pdb', '2f19.pdb', '1gix.pdb', '1bgs.pdb', '1ajd.pdb', '1ap6.pdb', '1ap5.pdb', '1g8l.pdb', '1hbu.pdb', '1hbs.pdb', '1xkj.pdb', '1enq.pdb', '1hbo.pdb', '1hbn.pdb', '1g8q.pdb', '1g8r.pdb', '1g8w.pdb', '1dc6.pdb', '1dc5.pdb', '1xps.pdb', '1dc3.pdb', '1flt.pdb', '1h31.pdb', '1f7y.pdb', '1f7z.pdb', '1feb.pdb', '1fec.pdb', '1f7n.pdb', '2avi.pdb', '1eqh.pdb', '1f7k.pdb', '1qvi.pdb', '1cbw.pdb', '1i7o.pdb', '1dcu.pdb', '3fbp.pdb', '1qqw.pdb', '1g72.pdb', '2sec.pdb', '1jy7.pdb', '1l3b.pdb', '1l3c.pdb', '1fz3.pdb', '1fz2.pdb', '1fz1.pdb', '1fz0.pdb', '1a50.pdb', '2dn2.pdb', '1fz4.pdb', '2dgc.pdb', '1bbb.pdb', '1e57.pdb', '1fzk.pdb', '1fzj.pdb', '1fzo.pdb', '1fzm.pdb', '1fza.pdb', '1fzg.pdb', '1fzf.pdb', '1fze.pdb', '1b5d.pdb', '1y85.pdb', '4bjl.pdb', '1c29.pdb', '1a5s.pdb', '1bks.pdb', '1bb0.pdb', '1fsi.pdb', '1oun.pdb', '1fsx.pdb', '1e5z.pdb', '1ky4.pdb', '1fhn.pdb', '1gx4.pdb', '2gct.pdb', '1dfl.pdb', '1i31.pdb', '1fwc.pdb', '1y4f.pdb', '2hbs.pdb', '1gxj.pdb', '1gxk.pdb', '1w7t.pdb', '1g44.pdb', '1fat.pdb', '1g40.pdb', '1fh2.pdb', '1bpy.pdb', '1fai.pdb', '1adv.pdb', '1fpl.pdb', '1rbb.pdb', '1y46.pdb', '4aah.pdb', '1f38.pdb', '1gjb.pdb', '1gjc.pdb', '1dpr.pdb', '1b94.pdb', '1b95.pdb', '1ajc.pdb', '1a6p.pdb', '1md0.pdb', '1hap.pdb', '3tgi.pdb', '1a64.pdb', '1a1b.pdb', '1fvk.pdb', '1gad.pdb', '1h4o.pdb', '1hao.pdb', '1v6o.pdb', '1v6n.pdb', '1v6m.pdb', '1adi.pdb', '1v6k.pdb', '1v6j.pdb', '1v6i.pdb', '2hpp.pdb', '2hpq.pdb', '1oc3.pdb', '1a07.pdb', '1bsl.pdb', '1bsm.pdb', '1y0d.pdb', '1f6p.pdb', '2gmf.pdb', '1y0w.pdb', '1y0t.pdb', '1xxt.pdb', '1dk1.pdb', '1eyu.pdb', '3bto.pdb', '1a01.pdb', '1eyy.pdb', '1pvi.pdb', '1lkr.pdb', '1cdm.pdb', '1cdd.pdb', '1gtt.pdb', '1aq8.pdb', '1gtq.pdb', '2dqj.pdb', '1bto.pdb', '1jgp.pdb', '1rfb.pdb', '1gsc.pdb', '1gsb.pdb', '1gsd.pdb', '1mtn.pdb', '1bjm.pdb', '1a2l.pdb', '1a2m.pdb', '1d01.pdb', '1fyz.pdb', '1ppb.pdb', '1g1m.pdb', '1h86.pdb', '1h84.pdb', '1g05.pdb', '1g08.pdb', '1g09.pdb', '1sar.pdb', '1avd.pdb', '1fy4.pdb', '1fy5.pdb', '1g9m.pdb', '2wsy.pdb', '1d0o.pdb', '1rv5.pdb', '1d0c.pdb', '1d0a.pdb', '1g0a.pdb', '1h33.pdb', '1fz5.pdb', '1d9c.pdb', '1d9g.pdb', '1ao6.pdb', '1gyp.pdb', '1ao3.pdb', '1tbe.pdb', '1adu.pdb', '1ade.pdb', '1qwh.pdb', '1jnu.pdb', '1aon.pdb', '1aom.pdb', '1aoj.pdb', '1jgo.pdb', '1bw8.pdb', '1jgq.pdb', '1b3r.pdb', '1gy6.pdb', '1aoq.pdb', '1i07.pdb', '1b5e.pdb']
        ##  proteins with 10 or more mutations relative to wt
        pdbs += ['7adh.pdb','1xac.pdb','1xad.pdb','1cx6.pdb','174l.pdb','1d3n.pdb','1hhl.pdb','192l.pdb','1a6i.pdb','1fbi.pdb','1lz2.pdb','1jhl.pdb','1sbt.pdb','2sbt.pdb']
        ## mutants
        pdbs += ['1vgk.pdb','2bst.pdb']
        pdbs += ['1lfv.pdb','1ye1.pdb']
        ## pdbs with more than 6 sequence identical chains after remark 350 transformation
        pdbs += ['1ht2.pdb','1hqy.pdb']
        ## pdbs with different peptide ligands
        pdbs += ['1t5z.pdb','1t79.pdb']
        pdbs += ['1p2c.pdb','1j1o.pdb']
        ## pdbs which caused problems with residue numbering
        pdbs += ['1bb0.pdb','1ca8.pdb','1a64.pdb','1a6p.pdb','1a3n.pdb','1bz1.pdb','1az0.pdb','1b94.pdb','1a50.pdb','1a5s.pdb','1a2l.pdb','1a2m.pdb','1fz3.pdb','1fz1.pdb','1hn4.pdb','1fx9.pdb','1dlf.pdb','2dlf.pdb','1q1r.pdb','1q1w.pdb','1csg.pdb','2gmf.pdb','1fia.pdb','3fis.pdb','1fph.pdb','1hao.pdb','1hap.pdb','1nrr.pdb']
        ## pdbs which caused problems with insertions in SEQRES sequence
        pdbs += ['1fph.pdb','103l.pdb']
##        self.l_pdbs = list(set(pdbs))

##        for i in range(len(self.l_pdbs)):
##            self.l_pdbs[i] = self.l_pdbs[i][:-4].lower()+'.pdb'

        self.l_pdbs.sort()

        self.l_pdbs = ['1h24.pdb','1h25.pdb'] ## temporary!!!

        try:
            ## remove pdbs spanning multiple pdb files/IDs
            self.l_pdbs.remove('1p0t.pdb'); self.l_pdbs.remove('1otz.pdb')
            ## remove pdbs with remark350 transformations refering to nonexisting chains
            self.l_pdbs.remove('1noj.pdb'); self.l_pdbs.remove('1b37.pdb')
            ## remove pdbs which only contain alpha carbon atoms
            self.l_pdbs.remove('1h6j.pdb')
        except:
            None

        self.pdbcount = len(self.l_pdbs)

        d_seq = self.parse_sequences()
        d_rmsd = self.analyze_sequences(d_seq)
        self.write_rmsd_to_file(d_rmsd, d_seq)
        print d_rmsd

        return


    def analyze_sequences(self, d_seq):

        import os, time
        import sys
        sys.path.append('/home/people/tc/python/Protool/')
        sys.path.append('../../Protool/')
        import geometry
        instance_geometry = geometry.geometry()

        print 'analyzing sequences'

        ## do not exclude pdb2 from pdb1 loop if sequence similar to previous pdb1
        ## since pdb A == B, A != C, B == C
        ## but do not analyze sequence of pdb A,B and then pdb B,A since A == B and B == A are equivalent

        d_rmsd_identical = {}
        d_rmsd_mutation = {}
        d_pdb = {}
        d_chains_intrapdb_sequence_identical = {}
        d_representative_chains = {}

        if not os.path.isfile('status.txt') and self.pdbcount > 10000:
            fd = open('status.txt','a')
            fd.write('0\n')
            fd.close()
        i1status = 0

        ##
        ## loop 1 over pdbs
        ##
        t1 = time.clock()
        for i1 in range(self.pdbcount-1):

            if i1 < i1status:
                continue

            if i1 in range(0,self.pdbcount,1000) and self.pdbcount > 10000:

                fd = open('status.txt','r')
                lines = fd.readlines()
                fd.close()
                i1status = int(lines[-1][:-1])

                if i1 == i1status:

                    fd = open('status.txt','a')
                    fd.write('%s\n' %(i1status+1000))
                    fd.close()
              
            pdb1 = self.l_pdbs[i1][:-4]
            SEQRESchains1 = d_seq[pdb1]['chains'].keys()

            ## continue if superseded structure
            if 'SPRSDE' in d_seq[pdb1].keys():
                continue

            ## continue if no polymer chains
            if SEQRESchains1 == []:
                continue

            ## continue if no protein chains
            if d_seq[pdb1]['proteinchains'] == []:
                continue

            ## continue if NMR or EM structure
            if d_seq[pdb1]['EXPDTA'] in ['NMR','CRYO-ELECTRON MICROSCOPY']:
                continue

            ## identify biomolecule(s)
            d_biomolecules1 = self.identify_biomolecule(pdb1, d_seq)
            biomolecules1 = d_biomolecules1.keys()

            ## identify water chains
            if 'REMARK525' in d_seq[pdb1].keys():
                waterchains1 = set(d_seq[pdb1]['REMARK525'])
            else:
                waterchains1 = set()

            ## identify sequence identical chains
            if pdb1 not in d_chains_intrapdb_sequence_identical.keys():
                d_chains_intrapdb_sequence_identical = self.identify_identical_chains_from_sequence_intra(
                    d_seq,pdb1,SEQRESchains1,waterchains1,
                    d_chains_intrapdb_sequence_identical,
                    )


            ##
            ## loop over biomolecule(s) of pdb1
            ##
            for biomolecule1 in biomolecules1:
                bmchains1 = d_biomolecules1[biomolecule1]['chains']
                bmtransformationcount1 = d_biomolecules1[biomolecule1]['transformants']

                
                ##
                ## loop 2 over pdbs
                ##
                for i2 in range(i1+1,self.pdbcount):
                    pdb2 = self.l_pdbs[i2][:-4]
                    SEQRESchains2 = d_seq[pdb2]['chains'].keys()

                    ## continue if superseded structure
                    if 'SPRSDE' in d_seq[pdb2].keys():
                        continue

                    ## continue if no polymer chains
                    if SEQRESchains2 == []:
                        continue

                    ## continue if no protein chains
                    if d_seq[pdb2]['proteinchains'] == []:
                        continue

                    ## continue if NMR or EM structure
                    if d_seq[pdb2]['EXPDTA'] == ['NMR','CRYO-ELECTRON MICROSCOPY']:
                        continue

                    ## identify biomolecule(s)
                    d_biomolecules2 = self.identify_biomolecule(pdb2, d_seq)
                    biomolecules2 = d_biomolecules2.keys()

                    ## identify water chains
                    if 'REMARK525' in d_seq[pdb2].keys():
                        waterchains2 = set(d_seq[pdb2]['REMARK525'])
                    else:
                        waterchains2 = set()

                    ## identify identical chains
                    if 'intra' not in d_seq[pdb2].keys():
                        d_chains_intrapdb_sequence_identical = self.identify_identical_chains_from_sequence_intra(
                            d_seq,pdb2,SEQRESchains2,waterchains2,d_chains_intrapdb_sequence_identical,
                            )
                        d_seq[pdb2]['intra'] = True

                    ## reset identified_chains_interpdb_sequence_similar for pdb2
                    identified_chains_interpdb_sequence_similar = False


                    ##
                    ## print status
                    ##
                    t2 = time.clock()
                    if t2-t1 > self.time_status_update:
                        print 'analyzing sequence of %s (%5i/%5i) and %s (%5i/%5i)' %(pdb1, i1+1, self.pdbcount, self.l_pdbs[i2][:-4], i2+1, self.pdbcount)
                        t1 = t2


                    ##
                    ## loop over biomolecule(s) of pdb2
                    ##
                    for biomolecule2 in biomolecules2:
                        bmchains2 = d_biomolecules2[biomolecule2]['chains']
                        bmtransformationcount2 = d_biomolecules2[biomolecule2]['transformants']


                        ##
                        ## skip if different number of chains in the biomolecule
                        ##
                        if len(bmchains1)*d_biomolecules1[biomolecule1]['transformants'] != len(bmchains2)*d_biomolecules2[biomolecule2]['transformants']:
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
                        ## compare hetero compounds assuming the following two statements to be correct
                        ## 1) "A particular HET group is represented in the PDB archives with a *unique* hetID."
                        ## 2) Depositors specify *all* hetero atoms observed in the electron density map.
                        if d_hetIDs[pdb1] != d_hetIDs[pdb2]:
                            continue


                        ##
                        ## identify sequence similar chains between pdbs (long peptides only)
                        ##
                        if identified_chains_interpdb_sequence_similar == False:
                            d_chains_interpdb_sequence_similar = self.identify_similar_chains_from_sequence_inter(
                                d_seq,
                                pdb1, pdb2,
                                d_chains_intrapdb_sequence_identical,
                                )
                            identified_chains_interpdb_sequence_similar = True


                        ##
                        ## continue if there are no sequence similar chains between the two pdbs
                        ##
                        if d_chains_interpdb_sequence_similar == {}:
                            continue

                        ## find chains which are not similar in between pdbs (if any)
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
##                            print pdb1, bmSEQRESchains1_not_similar_to_SEQRESchains2
##                            print pdb2, bmSEQRESchains2_not_similar_to_SEQRESchains1

                            ## 1) check that non sequence similar chains (if any) are not long peptides
                            d_bmSEQRESchains_not_similar_to_SEQRESchains = {
                                pdb1:bmSEQRESchains1_not_similar_to_SEQRESchains2,
                                pdb2:bmSEQRESchains2_not_similar_to_SEQRESchains1
                                }
                            peptide = True
                            short = False
                            for pdb in d_bmSEQRESchains_not_similar_to_SEQRESchains.keys():
                                for chain in d_bmSEQRESchains_not_similar_to_SEQRESchains[pdb]:
                                    ## check type of chain
                                    if chain not in d_seq[pdb]['chains'] and chain in d_seq[pdb]['HET'].keys():
                                        peptide = False
                                        break
                                    if d_seq[pdb]['chains'][chain]['type'] != 'peptide':
                                        peptide = False
                                        break
                                    ## check length of chain
                                    if len(d_seq[pdb]['chains'][chain]['seq']) < self.min_len_chain:
                                        short = True
                                        break
                                if peptide == False or short == True:
                                    break
                            ## continue if long peptide
                            if peptide == True and short == False:
                                continue

                            ## 2) check that nucleotides/saccharides or short peptides are sequence *identical*
                            sequence_identical = False
                            for chain1 in bmSEQRESchains1_not_similar_to_SEQRESchains2:
                                if chain1 not in d_seq[pdb1]['chains'] and chain1 in d_seq[pdb1]['HET'].keys():
                                    sequence_identical = True ## not true that identical for saccharide of 1q8p.pdb,1q8o.pdb (check LINK records)
                                    continue
                                for chain2 in bmSEQRESchains2_not_similar_to_SEQRESchains1:
                                    if chain2 not in d_seq[pdb1]['chains'] and chain2 in d_seq[pdb2]['HET'].keys():
                                        sequence_identical = True ## not true that identical for saccharide of 1q8p.pdb,1q8o.pdb (check LINK records)
                                        continue
                                    if d_seq[pdb1]['chains'][chain1]['seq'] == d_seq[pdb2]['chains'][chain2]['seq']:
                                        sequence_identical = True
                                    if sequence_identical == False:
                                        break
                                if sequence_identical == False:
                                    break
                            ## continue if not sequence identical
                            if sequence_identical == False:
                                continue
                                
                                    
                        ##
                        ## identify equivalent chains (interpdb) from structure
                        ##

                        ## parse coordinates
                        for s_pdb in [pdb1,pdb2]:
                            if s_pdb not in d_pdb.keys():
                                fd = open('%s%s.pdb' %(self.pdbpath, s_pdb.lower()),'r')
                                lines = fd.readlines()
                                fd.close()
                                d_pdb = self.parse_pdbcoordinatesection(lines, d_pdb, s_pdb)

                        ## identify equivalent chains
                        d_equivalent_chains, rmsd = self.identify_interpdb_equivalent_chains_from_structure(
                            pdb1, pdb2,
                            d_chains_intrapdb_sequence_identical,
                            d_chains_interpdb_sequence_similar,
                            d_pdb,
                            bmtransformationcount1, bmtransformationcount2,
                            d_seq,
                            bmchains1, bmchains2,
                            biomolecule1, biomolecule2,
                            )

                        ##
                        ## continue if there are no equivalent chains between the two pdbs
                        ##
                        if d_equivalent_chains == {}:
                            continue

                        set_equivalent_chains1 = set()
                        set_equivalent_chains2 = set()
                        for rep_chain1 in d_equivalent_chains.keys():
                            set_equivalent_chains1 |= set(d_equivalent_chains[rep_chain1][0])
                            set_equivalent_chains2 |= set(d_equivalent_chains[rep_chain1][1])
                        if len((set(bmchains1) - set_equivalent_chains1) ^ bmSEQRESchains1_not_similar_to_SEQRESchains2) > 0:
                            print pdb1, pdb2
                            print bmchains1
                            print set_equivalent_chains1
                            print set(bmchains1) - set_equivalent_chains1
                            print bmSEQRESchains1_not_similar_to_SEQRESchains2
                            print (set(bmchains1) - set_equivalent_chains1) ^ bmSEQRESchains1_not_similar_to_SEQRESchains2
                            print
                            print d_equivalent_chains
                            print bmchains1, bmchains2
                            print bmSEQRESchains1_not_similar_to_SEQRESchains2
                            print bmSEQRESchains2_not_similar_to_SEQRESchains1
                            expected
                        if len((set(bmchains2) - set_equivalent_chains2) ^ bmSEQRESchains2_not_similar_to_SEQRESchains1) > 0:
                            print set(bmchains2) ^ set_equivalent_chains2
                            print bmSEQRESchains2_not_similar_to_SEQRESchains1
                            print (set(bmchains2) ^ set_equivalent_chains2) ^ bmSEQRESchains2_not_similar_to_SEQRESchains1
                            print pdb1, pdb2
                            print d_equivalent_chains
                            print bmchains1, bmchains2
                            print bmSEQRESchains1_not_similar_to_SEQRESchains2
                            print bmSEQRESchains2_not_similar_to_SEQRESchains1
                            expected


                        if pdb1 not in d_rmsd_identical.keys():
                            d_rmsd_identical[pdb1] = {}
                        if biomolecule1 not in d_rmsd_identical[pdb1].keys():
                            d_rmsd_identical[pdb1][biomolecule1] = {}
                        if pdb2 not in d_rmsd_identical[pdb1][biomolecule1].keys():
                            d_rmsd_identical[pdb1][biomolecule1][pdb2] = {}
                        if biomolecule2 not in d_rmsd_identical[pdb1][biomolecule1][pdb2].keys():
                            d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2] = rmsd
                        fd = open('rmsdout.txt','a')
                        fd.write('%s %s %s\n' %(pdb1, pdb2, rmsd))
                        fd.close()
                        print pdb1, pdb2, rmsd, '\n'

    ##                            mutations = 0
    ##                            xcount = 0
    ##                            for i in range(len(seq1)):
    ##                                if seq1[i] != seq2[i]:
    ##                                    mutations += 1
    ##                                if seq1[i] == 'X':
    ##                                    xcount += 1
    ##                        
    ####                        if seq1aln.count('-') == 1 or seq2aln.count('-') == 1:
    ##                        if mutations == 1:
    ##                            fd = open('pointmutations.txt','a')
    ##                            fd.write('%s %s %s %s %s %s\n' %(s_pdb1, s_pdb2, chain1, chain2, xcount, len(seq1)))
    ##                            fd.close()
    ##                        elif mutations > 1 and mutations < 20:
    ##                            fd = open('multimutations.txt','a')
    ##                            fd.write('%s %s %s %s %s %s\n' %(s_pdb1, s_pdb2, chain1, chain2, mutations, len(seq1)))
    ##                            fd.close()

        return d_rmsd_identical


    def write_rmsd_to_file(self, d_rmsd, d_seq):

        import os

        lines = ['pdb1 pdb2 1 2  rmsd  pH1  pH2    T1    T2  res1  res2        SG1        SG2 HET1 HET2\n']

        for pdb1 in d_rmsd:

            ## parse physiochemical properties
            if 'REMARK200' in d_seq[pdb1].keys():
                T1 = d_seq[pdb1]['REMARK200']['TEMPERATURE']
                pH1 = d_seq[pdb1]['REMARK200']['PH']
            else:
                T1 = 'N/A'
                pH1 = 'N/A'
            ## parse x-ray space group
            spacegroup1 = d_seq[pdb1]['CRYST1']
            ## parse hetIDs
            hetIDs1 = set()
            ## parse resolution
            res1 = d_seq[pdb1]['REMARK2']
            for chain in d_seq[pdb1]['HET'].keys():
                hetIDs1 |= d_seq[pdb1]['HET'][chain]
            hetIDs1 = list(hetIDs1)

            for biomolecule1 in d_rmsd[pdb1].keys():

                for pdb2 in d_rmsd[pdb1][biomolecule1]:

                    ## parse physiochemical properties
                    if 'REMARK200' in d_seq[pdb2].keys():
                        T2 = d_seq[pdb2]['REMARK200']['TEMPERATURE']
                        pH2 = d_seq[pdb2]['REMARK200']['PH']
                    else:
                        T2 = 'N/A'
                        pH2 = 'N/A'
                    ## parse x-ray space group
                    spacegroup2 = d_seq[pdb2]['CRYST1']
                    ## parse hetIDs
                    hetIDs2 = set()
                    ## parse resolution
                    res2 = d_seq[pdb2]['REMARK2']
                    for chain in d_seq[pdb2]['HET'].keys():
                        hetIDs2 |= d_seq[pdb2]['HET'][chain]
                    hetIDs2 = list(hetIDs2)

                    for biomolecule2 in d_rmsd[pdb1][biomolecule1][pdb2].keys():

                        rmsd = d_rmsd[pdb1][biomolecule1][pdb2][biomolecule2]
                        lines += [
                            '%s %s %s %s %5.2f %s %s %5s %5s %5s %5s %s %s %s %s\n' %(
                                pdb1, pdb2, biomolecule1, biomolecule2, rmsd,
                                pH1.rjust(4), pH2.rjust(4),
                                T1.rjust(5), T2.rjust(5),
                                str(res1)[:5].rjust(5), str(res1)[:5].rjust(5),
                                spacegroup1.rjust(10), spacegroup2.rjust(10),
                                hetIDs1, hetIDs2,
                                )
                            ]

        if os.path.isfile('rmsd.txt'):
            lines = lines[1:]
        fd = open('rmsd.txt','a')
        fd.writelines(lines)
        fd.close()

        return


    def identify_biomolecule(self, pdb, d_seq):

        if 'REMARK350' in d_seq[pdb].keys():
            d_biomolecules = {}
            biomolecules = d_seq[pdb]['REMARK350'].keys()
            for biomolecule in biomolecules:
                d_biomolecules[biomolecule] = {}
                ## apply transformation to specific chains
                if 'chains' in d_seq[pdb]['REMARK350'][biomolecule].keys():
                    d_biomolecules[biomolecule]['chains'] = d_seq[pdb]['REMARK350'][biomolecule]['chains']
                ## apply transformations to all chains
                else:
                    d_biomolecules[biomolecule]['chains'] = d_seq[pdb]['chains'].keys()
                d_biomolecules[biomolecule]['transformants'] = d_seq[pdb]['REMARK350'][biomolecule]['transformants']
        ## assume everything to be the biomolecule
        else:
            d_biomolecules = {
                '1':{
                    'chains':d_seq[pdb]['chains'].keys(),'transformants':1
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

        SEQRESchains1_not_similar_to_SEQRESchains2 = set(bmchains1) - SEQRESchains1_similar_to_SEQRESchains2
        SEQRESchains2_not_similar_to_SEQRESchains1 = set(bmchains2) - SEQRESchains2_similar_to_SEQRESchains1

        return SEQRESchains1_not_similar_to_SEQRESchains2, SEQRESchains2_not_similar_to_SEQRESchains1


    def calculate_rmsd_for_multiple_chains(
        self,chains1,chains2,d_pdb,pdb1,pdb2,d_seq,
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
        print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s, %s, %s' %(rmsd, len(chains1), residue_count, len(coordinates1), chains1, chains2, pdb1, pdb2)

        return rmsd


    def identify_interpdb_equivalent_chains_from_structure(
        self, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_pdb,
        transformationcount1, transformationcount2,
        d_seq,
        bmchains1, bmchains2,
        biomolecule1, biomolecule2,
        ):

        ## return the following dictionary structure
        ## equivalent: {repchain1:[chains1,chains2]}

    ##    print 'identify_interpdb_equivalent_chains_from_structure', pdb1, pdb2

        d_equivalent_chains = {}

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
            l_chains1 += [bmchains1seqid]
            l_chains2 += [bmchains2seqid]

##        print pdb1, pdb2, l_chains1, l_chains2

        ## convert list of sequence identical chains to list of chains ordered by sequence identical chains
        chains1 = []
        for bmchains1seqid in l_chains1:
            chains1 += bmchains1seqid
        chains2 = []
        for bmchains2seqid in l_chains2:
            chains2 += bmchains2seqid

        ##
        ## identical number of chains and identical chain IDs
        ##
        if l_chains1 == l_chains2 and transformationcount1 == 1 and transformationcount2 == 1:

            rmsd = self.calculate_rmsd_for_multiple_chains(chains1,chains2,d_pdb,pdb1,pdb2,d_seq)
            d_equivalent_chains[rep_chain1] = [chains1, chains2]

        ##
        ## monomers
        ##
        elif (
            len(chains1) == 1 and len(chains2) == 1
            and
            transformationcount1 == 1 and transformationcount2 == 1
            ):

            rmsd = self.calculate_rmsd_for_multiple_chains(chains1,chains2,d_pdb,pdb1,pdb2,d_seq)
            d_equivalent_chains[rep_chain1] = [chains1, chains2]
            
        ##
        ## identical number of chains after REMARK350 transformation
        ##
        ## e.g. A1,B1,A2,B2,A3,B3 == A,B,C,D,E,F of 1xnv.pdb,1xo6.pdb
        ## e.g. B,B,B,B == B,D,B,D of 1vwr.pdb,1vwi.pdb
        ## e.g. A,C = B,B of 1my3.pdb,1mxu.pdb
        elif len(chains1)*transformationcount1 == len(chains2)*transformationcount2:

## use this to reduce number of 1xnv,1xo6 permutations
##                ## group intrapdb seq id chains to reduced number of possible permutations
##                rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']
##                l_equivalent_bmchains = [
##                    [
##                        set(d_chains_intrapdb_sequence_identical[pdb1][rep_chain1]+[rep_chain1]) & set(bmchains1),
##                        set(d_chains_intrapdb_sequence_identical[pdb2][rep_chain2]+[rep_chain2]) & set(bmchains2),
##                        ]
##                    ]

            d_biomolecules = {
                pdb1:{'biomolecule':biomolecule1,'chains':l_chains1,'tchains':[]},
                pdb2:{'biomolecule':biomolecule2,'chains':l_chains2,'tchains':[]},
                }
            for pdb in d_biomolecules.keys():
                if 'REMARK350' not in d_seq[pdb].keys():
                    d_biomolecules[pdb]['tchains'] = d_biomolecules[pdb]['chains']
                    continue
                biomolecule = d_biomolecules[pdb]['biomolecule']
                if 'matrices' not in d_seq[pdb]['REMARK350'][biomolecule].keys():
                    d_biomolecules[pdb]['tchains'] = d_biomolecules[pdb]['chains']
                    continue
                for bmchainsseqid in d_biomolecules[pdb]['chains']:
                    d_pdb, tchains = self.apply_remark350_transformation(d_pdb,pdb,d_seq,biomolecule,bmchainsseqid)
                    tchains = bmchainsseqid+tchains
                    d_biomolecules[pdb]['tchains'] += [tchains]

            l_tchains1 = d_biomolecules[pdb1]['tchains']
            l_tchains2 = d_biomolecules[pdb2]['tchains']
##            print pdb1, pdb2
##            print l_tchains1, l_tchains2
##            print l_chains1, l_chains2
##            stop

            ## convert list of sequence identical chains to list of chains ordered by sequence identical chains
            tchains1 = []
            for tbmchains1seqid in l_tchains1:
                tchains1 += tbmchains1seqid
            tchains2 = []
            for tbmchains2seqid in l_tchains2:
                tchains2 += tbmchains2seqid
            ## check if the expected correct combination of chains gives a low rmsd
            rmsd = self.calculate_rmsd_for_multiple_chains(tchains1,tchains2,d_pdb,pdb1,pdb2,d_seq)
            if rmsd < self.maxrmsd:
                d_equivalent_chains[rep_chain1] = [tchains1,tchains2]
                return d_equivalent_chains, rmsd

            ## permutations
            for i in range(len(l_tchains2)):
                bmchains2seqid = l_tchains2[i]
                if len(bmchains2seqid) > 6:
                    fd = open('toomanytransformations.txt','a')
                    fd.write('%s %s %s\n' %(pdb2, len(bmchains2seqid), pdb1))
                    fd.close()
                    return d_equivalent_chains, rmsd
                chains2permutations = self.permutation(bmchains2seqid)
                l_tchains2[i] = chains2permutations
            ## combination of permutations
            for i in range(len(l_tchains2)-1):
                j = i+1
                ipermutations = l_tchains2[i]
                jpermutations = l_tchains2[j]
                ## combination
                combinations = []
                for k in range(len(ipermutations)):
                    for l in range(len(jpermutations)):
                        combinations += [ipermutations[k]+jpermutations[l]]
                ## replace permutations with combinations
                l_tchains2[j] = combinations
                l_tchains2[i] = []
            chains2combinations = l_tchains2[-1]

            ## identify a correct combination of chains
            minrmsd = ['N/A','N/A']
            print 'chains2combinations', chains2combinations
            for i in range(len(chains2combinations)):
                chains2combination = chains2combinations[i]
                rmsd = self.calculate_rmsd_for_multiple_chains(tchains1,chains2combination,d_pdb,pdb1,pdb2,d_seq)
                if rmsd < minrmsd[1]:
                    minrmsd = [chains2combination,rmsd]
                    ## break to save time if many permutations
## change 4 to a variable...
                    if rmsd < self.maxrmsd and (len(tchains1) > 3 or i == 0):
                        break
            tchains2 = minrmsd[0]
            rmsd = minrmsd[1]
            d_equivalent_chains[rep_chain1] = [tchains1,tchains2]


        ##
        ## exceptions
        ##
        else:
            fd = open('notexpected_differentsized_biounits.txt','a')
            fd.write('%s %s %s %s\n' %(pdb1, pdb2, chains1, chains2))
            fd.close()
            rmsd = 'N/A'
            return d_equivalent_chains, rmsd

        return d_equivalent_chains, rmsd


    def apply_remark350_transformation(self,d_pdb,s_pdb,d_seq,biomolecule,bmchainsseqid):

        import Numeric

##        print 'applying REMARK350 transformation to chain %s of pdb %s' %(chain, s_pdb)

        tchains = []

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

                tchain = chain+str(i+2)
                tchains += [tchain]

                d_pdb[s_pdb]['chains'][tchain] = {}
                if 'residues' not in d_pdb[s_pdb]['chains'][tchain].keys():
                    d_pdb[s_pdb]['chains'][tchain]['residues'] = {}
                for res_no in d_pdb[s_pdb]['chains'][chain]['residues'].keys():
                    if 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]
                        continue
                    if res_no not in d_pdb[s_pdb]['chains'][tchain]['residues'].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no] = {}
                    d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['res_name'] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['res_name']
                    if 'atoms' not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no].keys():
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['atoms'] = {}
                    for atom_name in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'].keys():
                        if 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'][atom_name].keys():
                            d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['atoms'][atom_name] = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'][atom_name]
                            continue
                        if atom_name not in d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['atoms'].keys():
                            d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['atoms'][atom_name] = {}
                        coord = d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'][atom_name]['coordinate']
                        tcoord = Numeric.matrixmultiply(rmatrix, coord) + tvector
                        d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['atoms'][atom_name]['coordinate'] = tcoord

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
        self, d_seq, pdb, chains, waterchains,d_chains_intrapdb_sequence_identical,
        ):

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
                    d_chains_intrapdb_sequence_identical[pdb][chainj] = []

    ##    print 'intra', d_chains_intrapdb_sequence_identical

        return d_chains_intrapdb_sequence_identical


    def identify_similar_chains_from_sequence_inter(
        self, d_seq, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical
        ):

        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        sys.path.append('../../EAT_DB/')
        import sequence_alignment

        ## return the following dictionary structure
        ## interpdb: {repchain:{seqsimchain:{l1,l2}}}

        ## identify repchains
        repchains1 = d_chains_intrapdb_sequence_identical[pdb1].keys()
        repchains2 = d_chains_intrapdb_sequence_identical[pdb2].keys()

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

                if s1 == s2:
                    d_chains_interpdb_sequence_similar[chain1] = {'rep_chain2':chain2,'l1':l1,'l2':l2}
##                else:
##                    countmutations
##                    print chain1, chain2
##                    print s1, s2
##                    stop
##
##                    for i in range(len(s1)):

##                        if s1[i] != s2[i]:
##                            mutcount += 1
##                    if mutcount > 25 or mutcount/len(seq) > .25:
##                        notidentical

##        print 'inter', d_chains_interpdb_sequence_similar

        return d_chains_interpdb_sequence_similar


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

            ## ATOM record
            if 'insertion' in d_pdb[pdb]['chains'][chain]['residues'][res_no].keys():
                for i in range(len(d_pdb[pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes'])):
                    iCode = d_pdb[pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes'][i]
                    d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':iCode}
                    if iCode == ' ':
                        res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['res_name']
                    else:
                        res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['res_name']

                    seq += self.res_name2res_symbol(res_name)
                    ATOMrespos += 1
            ## REMARK record
            else:
                d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':' '}
                res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['res_name']
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
                        
            if iCode != ' ':
                res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['res_name']
            else:
                res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['res_name']

            res_symbol = self.res_name2res_symbol(res_name)

            ## append gap between ATOMseq and SEQRESseq
            if seq+res_symbol != SEQRESseq[:SEQRESpos+1]:
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}
                ATOMseqgaplen += 1
                seq += SEQRESseq[SEQRESpos]
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
                print SEQRESseq
                notexpected
                break

        return index


    def dcoordinates2lcoordinates(self, d_pdb, d_seq, pdb1, pdb2, chain1, chain2):

        coordinates1 = []
        coordinates2 = []

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
##            print ATOMseq1
##            print ATOMseq2
            if r2 == 0:
                ATOMseq1 = ATOMseq1[l2:]
            else:
                ATOMseq1 = ATOMseq1[l2:-r2]
            if r1 == 0:
                ATOMseq2 = ATOMseq2[l1:]
            else:
                ATOMseq2 = ATOMseq2[l1:-r1]
##            print ATOMseq1
##            print ATOMseq2
##            stop

        else:

            l1 = 0
            l2 = 0

        rescount = 0
        if len(ATOMseq1) != len(ATOMseq2):
            print pdb1, pdb2, chain1, chain2
            print ATOMseq1
            print ATOMseq2
            print l1, l2, r1, r2
            notexpected

        for SEQRESpos1 in range(l2,len(ATOMseq1)+l2):

            SEQRESpos2 = SEQRESpos1+l1-l2
            res_no1 = d_res_nos1[SEQRESpos1]['res_no']
            res_no2 = d_res_nos2[SEQRESpos2]['res_no']
            if res_no1 == '-' or res_no2 == '-':
                continue
            iCode1 = d_res_nos1[SEQRESpos1]['iCode']
            iCode2 = d_res_nos2[SEQRESpos2]['iCode']

            if 'REMARK' in d_pdb[pdb1]['chains'][chain1]['residues'][res_no1].keys():
                continue
            if 'REMARK' in d_pdb[pdb2]['chains'][chain2]['residues'][res_no2].keys():
                continue

            d_resname = {
                pdb1:{'chain':chain1,'res_no':res_no1,'iCode':iCode1},
                pdb2:{'chain':chain2,'res_no':res_no2,'iCode':iCode2},
                }
            for pdb in d_resname.keys():
                chain = d_resname[pdb]['chain']
                res_no = d_resname[pdb]['res_no']
                iCode = d_resname[pdb]['iCode']
                if iCode == ' ':
                    res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['res_name']
                    d_atoms = d_pdb[pdb]['chains'][chain]['residues'][res_no]['atoms']
                else:
                    res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['res_name']
                    d_atoms = d_pdb[pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['atoms']
                d_resname[pdb]['res_name'] = res_name
                d_resname[pdb]['d_atoms'] = d_atoms
            res_name1 = d_resname[pdb1]['res_name']
            res_name2 = d_resname[pdb2]['res_name']
            d_atoms1 = d_resname[pdb1]['d_atoms']
            d_atoms2 = d_resname[pdb2]['d_atoms']
            if self.pdbcount < 10000:
                print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2

            if res_name1 != res_name2:
                print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2
                print l1, l2, r1, r2
                print ATOMrespos1, ATOMrespos2
##                print resnogap1, resnogap2
                print 'SEQRES'
                print d_seq[pdb1]['chains'][chain1]['seq']
                print d_seq[pdb2]['chains'][chain2]['seq']
                print 'ATOM+i+gaps-aln'
                print ATOMseq1
                print ATOMseq2
                print 'ATOM+i+gaps+aln'
                print ATOMs1
                print ATOMs2
                print d_res_nos1[ATOMrespos1-1],d_res_nos1[ATOMrespos1],d_res_nos1[ATOMrespos2+1]
                print d_res_nos2[ATOMrespos2-1],d_res_nos2[ATOMrespos1],d_res_nos2[ATOMrespos2+1]
                stop
                fd = open('different_resnames.txt','a')
                fd.write('%s %s %s %s %s %s %s %s\n' %(pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2))
                fd.close()

            line470 = 'REMARK 470     %s %s%4i    ' %(res_name2, chain, res_no2)
            rescount += 1
            for atom_name in d_atoms1.keys():
                if atom_name not in d_atoms2.keys():
                    if atom_name[0] != 'H' and atom_name[:2] not in ['1H','2H','3H'] and atom_name != 'OXT':
                        line470 += '%s' %(atom_name.ljust(5))
                    continue
                if 'REMARK' in d_atoms1[atom_name].keys():
                    continue
                if 'REMARK' in d_atoms2[atom_name].keys():
                    continue
                ## append coordinates to list of coordinates
                coordinates1 += [d_atoms1[atom_name]['coordinate']]
                coordinates2 += [d_atoms2[atom_name]['coordinate']]
            if line470 != 'REMARK 470     %s %s%4i    ' %(res_name2, chain, res_no2):
                line470 = "                '"+line470+"', ##%s.pdb %s.pdb\n" %(pdb2, pdb1)
##                fd = open('missingatoms.txt','a')
##                fd.write(line470)
##                fd.close()

        return coordinates1, coordinates2, rescount


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


    def parse_sequences(self):

        import time

        print 'parsing data from noncoordinate sections'
        
        d_seq = {}
        d_seq_identical = {}

        t1 = time.clock()
        for i in range(self.pdbcount):

            pdb = self.l_pdbs[i][:-4]

            ## print status
            t2 = time.clock()
            if t2-t1 > self.time_status_update:
                print 'parsing %s (%s/%s)' %(self.l_pdbs[i-1][:-4], i, self.pdbcount)
                t1 = t2

            ## read lines
            fd = open('%s%s.pdb' %(self.pdbpath, pdb.lower()),'r')
            lines = fd.readlines()
            fd.close()

            ## parse data prior to the coordinate section
            d_seq = self.parse_pdbnoncoordinatesections(lines, d_seq, pdb)

        return d_seq


    def parse_pdbnoncoordinatesections(self, lines, d_seq, s_pdb):

        ## import sequence alignment
        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        sys.path.append('../../EAT_DB/')
        import sequence_alignment

        ## parser written on the assumption that SEQRES is mandatory if ATOM records exist

        if s_pdb in self.d_missinglines.keys():
            lines = self.d_missinglines[s_pdb]+lines

        ## correct errornous lines
        if s_pdb in self.d_noncoordinateerrors.keys():
            for i in range(len(lines)):
                line = lines[i]
                for j in range(len(self.d_noncoordinateerrors[s_pdb]['incorrect'])):
                    incorrectline = self.d_noncoordinateerrors[s_pdb]['incorrect'][j]
                    if line.strip() == incorrectline:
                        correctline = self.d_noncoordinateerrors[s_pdb]['correct'][j]
                        line = lines[i] = correctline
                        break


        if s_pdb not in d_seq.keys():
            d_seq[s_pdb] = {
                }

        ## insertion chain, res_no, res_name
        prev_chain = ''
        d_insertions = {}
        biounit = 'N/A'

        for i in range(len(lines)):
            line = lines[i]

##            if self.pdbcount > 10000:
##                if 'BIOLOGICAL_UNIT' in line or 'BIOLOGICAL UNIT' in line:
##                    for s_biounit in self.l_biounits:
##                        if s_biounit in line:
##                            biounit = s_biounit
##                    if biounit == 'N/A':
##                        for s_biounit in self.d_biounits.keys():
##                            if s_biounit in line:
##                                biounit = s_biounit
##                    if biounit == 'N/A':
##                        print lines[i-5:i+5]
##                        print s_pdb
##                    biologicalunit = 'unknown'
##                    for biounit in self.d_biounits.keys()+[
##                        'KNOWN','OLIGOMER','CHAIN','SUBUNITS','MULTIMER','-MER',
##                        '1MER','2MER','2MER','3MER','4MER','5MER','6MER','7MER','8MER','9MER',
##                        ]:
##                        for j in range(i-2,i+3):
##                            if biounit in lines[j]:
##                                biologicalunit = 'known'
##                                break
##                        if biologicalunit == 'known':
##                            break
##                    if biologicalunit == 'unknown' and 'REMARK 300' not in line:
##                        interestinglines = lines[i-2]+lines[i-1]+lines[i]+lines[i+1]+lines[i+2]
####                        fd = open('interestinglines.txt','a')
####                        fd.write('%s\n%s' %(s_pdb, interestinglines))
####                        fd.close()

            record = line[:6].strip()

            if record == 'ATOM': ## section 9
                continue

            elif record == 'HETATM': ## section 9
                continue

            elif record == 'REMARK': ## section 2
                d_seq = self.parse_recordREMARK(d_seq, line, s_pdb, i, lines)

            elif record == 'SEQRES': ## section 3
                d_seq = self.parse_recordSEQRES(line, d_seq, s_pdb)

            elif record == 'HET': ## section 4
                hetID = line[7:10].strip()
                ## continue if water
                if hetID in ['H20','HOH','D2O','DOD']: ## D2O in 2JAJ
                    continue
##                ## continue if ion
##                if hetID in self.d_ions.keys():
##                    continue
                chain = line[12]
                if 'HET' not in d_seq[s_pdb].keys():
                    d_seq[s_pdb]['HET'] = {}
                if chain not in d_seq[s_pdb]['HET'].keys():
                    d_seq[s_pdb]['HET'][chain] = set()
                d_seq[s_pdb]['HET'][chain] |= set([hetID])

            elif record == 'TITLE': ## section 2
                d_seq[s_pdb]['TITLE'] = line[10:].strip()

            elif record == 'HEADER':
                d_seq[s_pdb]['HEADER'] = line[10:].strip()

            elif record == 'SPRSDE': ## section 2
                sIDcodes = line[31:70].lower().split()
                for sIDcode in sIDcodes:
                    if sIDcode not in d_seq.keys():
                        d_seq[sIDcode] = {}
                    d_seq[sIDcode]['SPRSDE'] = True

            elif record == 'EXPDTA': ## section 2
                methods = line[10:]
                if methods[:3] == 'NMR':
                    d_seq[s_pdb]['EXPDTA'] = 'NMR'
                elif methods.strip() == 'CRYO-ELECTRON MICROSCOPY':
                    d_seq[s_pdb]['EXPDTA'] = 'CRYO-ELECTRON MICROSCOPY'
                else:
                    d_seq[s_pdb]['EXPDTA'] = 'N/A'

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                d_seq[s_pdb]['CRYST1'] = spacegroup

        for key in ['TITLE','EXPDTA','REMARK2']:
            if key not in d_seq[s_pdb].keys():
                d_seq[s_pdb][key] = 'N/A'
        for key in ['chains','HET']:
            if key not in d_seq[s_pdb].keys():
                d_seq[s_pdb][key] = {}

        proteinchains = []
        for chain in d_seq[s_pdb]['chains'].keys():
            if d_seq[s_pdb]['chains'][chain]['type'] == 'peptide' and len(d_seq[s_pdb]['chains'][chain]['seq']) > self.min_len_chain:
                proteinchains += chain
        d_seq[s_pdb]['proteinchains'] = proteinchains

##        if self.pdbcount > 10000:
##            if 'REMARK350' not in d_seq[s_pdb].keys() and len(proteinchains) > 1:
##                if biounit == 'N/A':
##                    fd = open('unknownbiounit.txt','a')
##                    fd.write('%s\n' %(s_pdb))
##                    fd.close()
##                elif self.d_biounits[biounit] != len(proteinchains):
##                    if biounit == 'MONOMER':
##                        d_seq[s_pdb]['REMARK350'] = {}
##                        for i in range(len(proteinchains)):
##                            chain = proteinchains[i]
##                            d_seq[s_pdb]['REMARK350'][i+1] = {'transformants': 1, 'chains': [chain]}
##                    else:
##                        print len(proteinchains) % self.d_biounits[biounit]
##                        print biounit, self.d_biounits[biounit]
##                        print s_pdb
##                        print d_seq[s_pdb].keys()
##                        print d_seq[s_pdb]['chains'].keys()
##                        print proteinchains
##                        stop3
####                d_seq[s_pdb]['biounit'] = 'small' ## small opposed to monomer if multiple different chains
#### count number of similar chains by comparing SEQRESseqs

        return d_seq


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):
            matrix = d_seq[s_pdb]['REMARK350'][biomolecule]['matrices'][i]
            rmatrix, tvector = self.transformationmatrix2rotationmatrix_and_translationvector(matrix)


    def parse_recordSEQRES(self, line, d_seq, s_pdb):

        chain = line[11]

        if 'chains' not in d_seq[s_pdb]:
            d_seq[s_pdb]['chains'] = {}
        if chain not in d_seq[s_pdb]['chains'].keys():
            d_seq[s_pdb]['chains'][chain] = {}
        if not 'type' in d_seq[s_pdb]['chains'][chain].keys():
            d_seq[s_pdb]['chains'][chain]['type'] = 'unknown'

        residues = line[19:70].split()

        for i in range(len(residues)):
            residue = residues[i]
            if residue in self.d_res.keys():
                if d_seq[s_pdb]['chains'][chain]['type'] == 'unknown':
                    d_seq[s_pdb]['chains'][chain]['type'] = 'peptide'
                residues[i] = self.d_res[residue]
            elif residue in ['C','A','T','G']:
                if d_seq[s_pdb]['chains'][chain]['type'] in ['unknown','peptide']:
                    d_seq[s_pdb]['chains'][chain]['type'] = 'nucleotide'
                residues[i] = residue
            elif residue in ['GLC']:
                if d_seq[s_pdb]['chains'][chain]['type'] in ['unknown','peptide']:
                    d_seq[s_pdb]['chains'][chain]['type'] = 'saccharide'
                residues[i] = residue
            else:
                residues[i] = 'X'

        if 'seq' not in d_seq[s_pdb]['chains'][chain].keys():
            d_seq[s_pdb]['chains'][chain]['seq'] = ''
        d_seq[s_pdb]['chains'][chain]['seq'] += ''.join(residues)

        return d_seq


    def parse_pdbcoordinatesection(self, lines, d_pdb, s_pdb):

##        print 'parsing coordinates of %s' %(s_pdb)

        ## append missing lines
        if s_pdb in self.d_missinglines.keys():
            lines = self.d_missinglines[s_pdb]+lines

        ## correct errornous lines
        if s_pdb in self.d_coordinateerrors.keys():
            for i in range(len(lines)):
                line = lines[i]
                for j in range(len(self.d_coordinateerrors[s_pdb]['incorrect'])):
                    incorrectline = self.d_coordinateerrors[s_pdb]['incorrect'][j]
                    if line.strip() == incorrectline:
                        correctline = self.d_coordinateerrors[s_pdb]['correct'][j]
                        line = lines[i] = correctline
                        break

        d_pdb[s_pdb] = {
            'chains':{},
            'HET':set(),
            }

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()
    
            if record == 'ATOM':
                d_pdb = self.parse_recordATOM(line, d_pdb, s_pdb, lines, i)

            elif record == 'REMARK':
                remark = int(line[6:10])
                if remark == 465:
                    d_pdb = self.parse_recordREMARK465(line, d_pdb, s_pdb, lines, i)
                elif remark == 470:
                    d_pdb = self.parse_recordREMARK470(line, d_pdb, s_pdb, lines, i)

            elif record == 'HET':
                HETID = line[7:10].strip()
                d_pdb[s_pdb]['HET'] |= set([HETID])

            elif record == 'MODEL':
                model = int(line.split()[1])

        return d_pdb


    def write_REMARK350_chains(self, line_chains,d_seq,lines,i,s_pdb):

        chains = self.parse_REMARK350_chains(line_chains,s_pdb,lines,i)
        biomolecules = ['1']
        for j in range(i-1,-1,-1):
            if lines[j][:10] != 'REMARK 350':
                break
            elif lines[j][11:23] == 'BIOMOLECULE:':
                biomolecules = lines[j][23:80].replace(' ','').split(',')
                break
        for biomolecule in biomolecules:
            if biomolecule not in d_seq[s_pdb]['REMARK350'].keys():
                d_seq[s_pdb]['REMARK350'][biomolecule] = {}
            if 'transformants' not in d_seq[s_pdb]['REMARK350'][biomolecule].keys():
                d_seq[s_pdb]['REMARK350'][biomolecule]['transformants'] = 1
            if 'chains' not in d_seq[s_pdb]['REMARK350'][biomolecule].keys():
                d_seq[s_pdb]['REMARK350'][biomolecule]['chains'] = []
            d_seq[s_pdb]['REMARK350'][biomolecule]['chains'] += chains

        return d_seq


    def parse_REMARK350_chains(self, line_chains,s_pdb,lines,i):

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
            if chain != 'NULL' and len(chain) > 1:
                print chains, s_pdb
                print lines[i]
                stop

        return chains


    def parse_recordREMARK(self, d_seq, line, s_pdb, i, lines):

        remark = int(line[6:10])

        if remark == 200:

            if line[12:23].strip().upper() in ['TEMPERATURE','PH']:
                experimentaldetail_key = line[12:23].strip()
                experimentaldetail_value = line[44:].strip()
                if 'REMARK200' not in d_seq[s_pdb].keys():
                    d_seq[s_pdb]['REMARK200'] = {}
                if experimentaldetail_key not in d_seq[s_pdb]['REMARK200'].keys():
                    d_seq[s_pdb]['REMARK200'][experimentaldetail_key] = experimentaldetail_value

        elif remark == 350: ## biological units

            if 'REMARK350' not in d_seq[s_pdb].keys():
                d_seq[s_pdb]['REMARK350'] = {}

    ##        if line[11:23] == 'BIOMOLECULE:':
    ##
    ##            biomolecule = int(line[23:80])
    ##            d_seq[s_pdb]['REMARK350'][biomolecule] = []

            if line[11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                line_chains = line[41:80]
                d_seq = self.write_REMARK350_chains(line_chains,d_seq,lines,i,s_pdb)

            ## count and parse chain transformations
            ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj)
            elif line[13:19] in ['BIOMT3','SMTRY3']:
                biomolecules = ['1']
                
                ## do not count chain transformation if no transformation
                ## upon application of transformation matrix
                matrixrow1 = lines[i-2][24:].split()
                matrixrow2 = lines[i-1][24:].split()
                matrixrow3 = lines[i-0][24:].split()
                matrixrows = [matrixrow1,matrixrow2,matrixrow3,]
                transformation = False
                for j in range(3):
                    ## add a zero translation vector if a translation vector is not given
                    if len(matrixrows[j]) == 3:
                        matrixrows[j] += [0.]
                    if float(matrixrows[j][j]) == 1. and float(matrixrows[j][3]) == 0.:
                        continue
                    else:
                        transformation = True
                if transformation == False:
                    return d_seq

                ## identify biomolecule(s)
                for j in range(i-1,-1,-1):
                    if lines[j][:10] != 'REMARK 350':
                        break
                    elif lines[j][11:23] == 'BIOMOLECULE:':
                        biomolecules = lines[j][23:80].replace(' ','').split(',')
                        break

                ## loop over biomolecule(s)
                for biomolecule in biomolecules:
                    if biomolecule not in d_seq[s_pdb]['REMARK350'].keys():
                        d_seq[s_pdb]['REMARK350'][biomolecule] = {}
                    ## count transformation matrix
                    if 'transformants' not in d_seq[s_pdb]['REMARK350'][biomolecule].keys():
                        d_seq[s_pdb]['REMARK350'][biomolecule]['transformants'] = 1
                    d_seq[s_pdb]['REMARK350'][biomolecule]['transformants'] += 1
                    ## parse transformation matrix
                    if 'matrices' not in d_seq[s_pdb]['REMARK350'][biomolecule].keys():
                        d_seq[s_pdb]['REMARK350'][biomolecule]['matrices'] = []
                    d_seq[s_pdb]['REMARK350'][biomolecule]['matrices'] += [matrixrows]


            elif ',' in line[11:80]:

                ## check if the line is a line of chains by checking
                ## if the line is preceded by 'APPLY THE FOLLOWING TO CHAINS:'
                ## and not preceded by a blank line (e.g. a problem with 1m4x, 1d3i)
                chains = False
                for j in range(i-1,-1,-1):
                    if lines[j][:10] != 'REMARK 350':
                        break
                    if lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                        chains = True
                        break
                    if lines[j][10:].strip() == '':
                        break
                if chains == False:
                    return d_seq

                ## add chains to biomolecule
                line_chains = line[11:80]
                d_seq = self.write_REMARK350_chains(line_chains,d_seq,lines,i,s_pdb)


        elif remark == 525: ## water association

            if line[11:].strip() == 'PROTEIN CHAIN  SOLVENT CHAIN':
                for j in range(i+1,len(lines)):
                    if lines[j][11:].strip() == '':
                        break
                    if lines[j][:10] != 'REMARK 525':
                        break
                    else:
                        solventchain = lines[j][11:].split()[1]
                        if 'REMARK525' not in d_seq[s_pdb].keys():
                            d_seq[s_pdb]['REMARK525'] = []
                        d_seq[s_pdb]['REMARK525'] += [solventchain]

        elif remark == 2: ## resolution
            try:
                resolution = float(line[22:27])
            except:
                resolution = 'N/A'
            d_seq[s_pdb]['REMARK2'] = resolution

        return d_seq


    def parse_recordREMARK465(self, line, d_pdb, s_pdb, lines, i):

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
                iCode = lines[j][26:28]

                if not chain in d_pdb[s_pdb]['chains'].keys():
                    d_pdb[s_pdb]['chains'][chain] = {}
                if not 'residues' in d_pdb[s_pdb]['chains'][chain].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'] = {}
                if not res_no in d_pdb[s_pdb]['chains'][chain]['residues'].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no] = {}
                ## res_no > insertion
                if not 'insertion' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion'] = {}
                ## insertion > l_iCodes
                if not 'l_iCodes' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion'].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes'] = []
                ## l_iCodes > iCode
                if not iCode in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes']:
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes'] += [iCode]
                ## insertion > d_iCodes
                if not 'd_iCodes' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion'].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'] = {}
                ## d_iCodes > iCode > res_name
                if not iCode in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode] = {
                        'res_name':res_name,'REMARK':True,
                        }
##                ## res_no > res_name
##                if not 'res_name' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
##                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['res_name'] = res_name
##                ## res_no > REMARK
##                if not 'REMARK' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
##                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['REMARK'] = res_name

        elif line[10:].strip().split() == ['M','RES','C','SSSEQI']:
            notexpected

        return d_pdb


    def parse_recordREMARK470(self, line, d_pdb, s_pdb, lines, i):

        ## missing atoms

        ## the latter equation is only to acommodate for 1fvk.pdb
        if line[10:].strip() == 'M RES CSSEQI  ATOMS' or line[10:].strip() == 'M RES C SEQI  ATOMS':

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 470':
                    break

                try:
                    model = int(lines[j][11:13])
                except:
                    model = 'N/A'
                res_name = lines[j][15:18]
                chain = lines[j][19]
                try:
                    res_no = int(lines[j][20:24])
                except:
                    res_no = lines[j][20:24].strip()
                atoms = lines[j][24:].split()

                if not chain in d_pdb[s_pdb]['chains'].keys():
                    d_pdb[s_pdb]['chains'][chain] = {}
                if not 'residues' in d_pdb[s_pdb]['chains'][chain].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'] = {}
                if not res_no in d_pdb[s_pdb]['chains'][chain]['residues'].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no] = {}
                if not 'atoms' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'] = {}
                if not 'res_name' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['res_name'] = res_name
                for atom_name in atoms:
                    if not atom_name in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'].keys():
                        d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'][atom_name] = {}
                    d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['atoms'][atom_name]['REMARK'] = 470

        elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
            notexpected

        return d_pdb


    def parse_recordATOM(self, line, d_pdb, s_pdb, lines, i):

        import Numeric
        
        atom_name = line[12:16].strip()
    ##    atom_altloc = line[16].strip()
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        atom_x = float(line[30:38].strip())
        atom_y = float(line[38:46].strip())
        atom_z = float(line[46:54].strip())
        coordinate = Numeric.array([atom_x, atom_y, atom_z])
        if not chain in d_pdb[s_pdb]['chains'].keys():
            d_pdb[s_pdb]['chains'][chain] = {}
        if not 'residues' in d_pdb[s_pdb]['chains'][chain].keys():
            d_pdb[s_pdb]['chains'][chain]['residues'] = {}

        ## res_no
        if not res_no in d_pdb[s_pdb]['chains'][chain]['residues'].keys():
            d_pdb[s_pdb]['chains'][chain]['residues'][res_no] = {}
        ## res_no > insertion
        if not 'insertion' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no].keys():
            d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion'] = {}
        ## insertion > l_iCodes
        if not 'l_iCodes' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion'].keys():
            d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes'] = []
        ## l_iCodes > iCode
        if not iCode in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes']:
            d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['l_iCodes'] += [iCode]


        if s_pdb == '1h25' and chain == 'D' and res_no == 175:
            print res_name
            stop

        ##
        ## insertion
        ##
        else:

            ## insertion > d_iCodes
            if not 'd_iCodes' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion'].keys():
                d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'] = {}
            ## d_iCodes > iCode
            d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode] = {}
            ## iCode > atoms
            if not 'atoms' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode].keys():
                d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['atoms'] = {}
            ## atoms > atom_name > coordinate
            if not atom_name in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['atoms'].keys():
                d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}
            ## iCode > res_name
            if not 'res_name' in d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode].keys():
                d_pdb[s_pdb]['chains'][chain]['residues'][res_no]['insertion']['d_iCodes'][iCode]['res_name'] = res_name

        return d_pdb


    def faculty(self, n):

        fac = 1
        for i in range(1,n+1):
            fac *= i

        return fac


    def __init__(self):

        import os

        self.d_res = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y'}
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

        self.d_alphabet = {
            'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8,'I':9,'J':10,'K':11,'L':12,'M':13,
            'N':14,'O':15,'P':16,'Q':17,'R':18,'S':19,'T':20,'U':21,'V':22,'W':23,'X':24,'Y':25,'Z':26,
            }
        self.s_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

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

        self.maxrmsd = 2.5

        self.d_missinglines = {

            '2g34': ## 2g33.pdb
            [
                'HET    HAP                                                                      ',
                ]
            ,
            '1bgs': ## 1brs.pdb
            [
                'REMARK 350 BIOMOLECULE: 1',
                'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, E',
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000',
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000',
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000',
                'REMARK 350 BIOMOLECULE: 2',
                'REMARK 350 APPLY THE FOLLOWING TO CHAINS: B, F',
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000',
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000',
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000',
                'REMARK 350 BIOMOLECULE: 3',
                'REMARK 350 APPLY THE FOLLOWING TO CHAINS: C, G',
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000',
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000',
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000',
                ]
            ,
            '1brs':
            [
                'REMARK 350 BIOMOLECULE: 1',
                'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, D',
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000',
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000',
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000',
                'REMARK 350 BIOMOLECULE: 2',
                'REMARK 350 APPLY THE FOLLOWING TO CHAINS: B, E',
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000',
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000',
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000',
                'REMARK 350 BIOMOLECULE: 3',
                'REMARK 350 APPLY THE FOLLOWING TO CHAINS: C, F',
                'REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000',
                'REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000',
                'REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000',
                ]
            }

        self.d_coordinateerrors = {
            '1l3b':
            {
                'incorrect':
                [
                    'ATOM   3519  N   ASP C   0      28.351 -21.054 125.451  1.00 26.35           N',
                    'ATOM   3520  CA  ASP C   0      28.823 -20.517 124.181  1.00 25.16           C',
                    'ATOM   3521  C   ASP C   0      28.676 -18.991 124.227  1.00 23.83           C',
                    'ATOM   3522  O   ASP C   0      27.843 -18.464 124.972  1.00 22.27           O',
                    'ATOM   3523  CB  ASP C   0      28.001 -21.094 123.023  1.00 29.03           C',
                    'ATOM   3524  CG  ASP C   0      27.871 -22.613 123.088  1.00 32.09           C',
                    'ATOM   3525  OD1 ASP C   0      28.844 -23.282 123.499  1.00 36.68           O',
                    'ATOM   3526  OD2 ASP C   0      26.798 -23.143 122.712  1.00 34.24           O',
                    ]
                ,
                'correct':
                [
                    'ATOM   3519  N   ASP C 100      28.351 -21.054 125.451  1.00 26.35           N',
                    'ATOM   3520  CA  ASP C 100      28.823 -20.517 124.181  1.00 25.16           C',
                    'ATOM   3521  C   ASP C 100      28.676 -18.991 124.227  1.00 23.83           C',
                    'ATOM   3522  O   ASP C 100      27.843 -18.464 124.972  1.00 22.27           O',
                    'ATOM   3523  CB  ASP C 100      28.001 -21.094 123.023  1.00 29.03           C',
                    'ATOM   3524  CG  ASP C 100      27.871 -22.613 123.088  1.00 32.09           C',
                    'ATOM   3525  OD1 ASP C 100      28.844 -23.282 123.499  1.00 36.68           O',
                    'ATOM   3526  OD2 ASP C 100      26.798 -23.143 122.712  1.00 34.24           O',
                    ]
                }
            ,
            '1zpt':
            {
                'incorrect':
                [
                    'REMARK 465     MET C     6',
                    'REMARK 465     SER C     7',
                    'REMARK 465     PHE C     8',
                    'REMARK 465     PHE C     9',
                    'REMARK 465     HIS C    10',
                    'REMARK 465     ALA C    11',
                    'REMARK 465     SER C    12',
                    'REMARK 465     GLN C    13',
                    'REMARK 465     ARG C    14',
                    'REMARK 465     ASP C    15',
                    'REMARK 465     ALA C    16',
                    'REMARK 465     LEU C    17',
                    'REMARK 465     ASN C    18',
                    'REMARK 465     GLN C    19',
                    'REMARK 465     SER C    20',
                    'REMARK 465     LEU C    21',
                    ]
                ,
                'correct':
                [
                    'REMARK 465     MET C     1',
                    'REMARK 465     SER C     2',
                    'REMARK 465     PHE C     3',
                    'REMARK 465     PHE C     4',
                    'REMARK 465     HIS C     5',
                    'REMARK 465     ALA C     6',
                    'REMARK 465     SER C     7',
                    'REMARK 465     GLN C     8',
                    'REMARK 465     ARG C     9',
                    'REMARK 465     ASP C    10',
                    'REMARK 465     ALA C    11',
                    'REMARK 465     LEU C    12',
                    'REMARK 465     ASN C    13',
                    'REMARK 465     GLN C    14',
                    'REMARK 465     SER C    15',
                    'REMARK 465     LEU C    16',
                    ]
                }
            ,
            }

        self.d_noncoordinateerrors = {
            '1h25':
            {
                'incorrect':
                [
                    'SEQRES   5 I   71  VAL ARG VAL PHE TYR ASN PRO GLY THR ASN VAL VAL ASN  2SEC 116',
                    ]
                ,
                'correct':
                [
                    'SEQRES   5 I   71  VAL ARG VAL PHE TYR ASN PRO     THR ASN VAL VAL ASN  2SEC 116',
                    ]
                }
##            '2sec':
##            {
##                'incorrect':
##                [
##                    'SEQRES   5 I   71  VAL ARG VAL PHE TYR ASN PRO GLY THR ASN VAL VAL ASN  2SEC 116',
##                    ]
##                ,
##                'correct':
##                [
##                    'SEQRES   5 I   71  VAL ARG VAL PHE TYR ASN PRO     THR ASN VAL VAL ASN  2SEC 116',
##                    ]
##                }
            }
                

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25

        self.pdbpath = '/oxygenase_local/data/pdb/'
        self.l_pdbs = os.listdir(self.pdbpath)
        
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
## use biological molecules and not asymmetric units nor monomers for rmsd calculation (e.g. to accomodate for 2c1o with chains of different multimeric biological molecules in the asymmetric unit)
##
## asymmetric units in the same pdb are not compared
##
## pdbs with different hetero compounds are not compared
## pdbs with different nonhetero compunds (peptide, nucleotide, saccharide) are not compared (e.g. 1ok7 vs 1mmi)


            ## different oxidation states
            ## 1oc3:C

            ## different potassium ("2" vs "6") concentrations
            ## 2hvj FORMUL 4 K 2(K1 1+)
            ## 2hvk FORMUL 4 K 6(k1 1+)

            ## different quarternary structures
            ## 1c77, 1c78, 1c79

            ## incorrect remark 350 translation vector
            ## 1bks,2wsy

#!/software/bin/python
#
#$Id: quakes.py,v 1.2 2007/07/04 16:08:00 tc Exp $
#
#Tommy Carstensen, University College Dublin, 2007

##
## questions of interest

## find correlation (if any) between RMSD of backbone/all atoms for all/neighboring(exponentional sphere radii)/surface residues

## rmsd not mathematically dependent on chain length! but rmsd physically dep on chain length? plot!

## rmsd plot of structures with "bad" and "good" geometry!

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

## maltotetraose (GLC) in 2fhf; maltotriose (MLR=3xGLC) in 2fhc; 2fhb, 2fh8, 2fh6
## parse some sugar molecules parsed and leave the rest if hetID of monomers and not multimer specified
## look up sugar dictionary...

## what to do with 1k56, 1k57?

## what to do with 2huk (v131c mutant) with ssbond between monomers  (cys131) of T4 lysozyme dimer?

## incorrect d_res_nos will be returned from ATOM2seq for 2bfk (vs 2bfl), chain A becauseof ASN61A in REMARK465 records!!!
## iCode ' ' from ATOM records should be put before iCode 'A' from REMARK465 records by comparsion to the SEQRES records upon parsing the ATOM records

## write a new faster seq aln alg

## if no remark350 record and high rmsd then try monomers...

## don't do seq aln during rmsd2bfactor...
## don't do 1st SEQRESseq aln for seq id chains?!

## accept neutron diffraction structures? accept 3ins which is x-ray AND neutron?

## proteins with 10 or more mutations relative to wt
## ['7adh','1xac','1xad','1cx6','174l','1d3n','1hhl','192l','1a6i','1fbi','1lz2','1jhl','1sbt','2sbt']

## what to do with 1ft8:E and 1koh:B for which coordinates of aligned residues are not given???

## solve the problem of too many chain combinations by sequential pairing
## this will yield a maximum of n(n+1)/2 combinations to check and with fewer coordinates per rmsd calculation
## sequential addition of sequence similar chains to reduce number of combinations further

## plot of rmsd vs disulphide bonds (y/n)
## plot of rmsd vs EC class
## plot of rmsd vs CATH class

import os, sys, Numeric, LinearAlgebra, math

class quakes:

    def main(self):

##        self.rsync()
##        self.gunzip()

        errorpdbs = [
##            ## pdb errors
            '1fng','1fne', ## incorrect use of residue numbers and iCodes
            '1yra','1yrb', ## residue renumbering
            '1t41', ## residue renumbering
            '2ht8', ## residue renumbering
            '1t7b', '1t7e', ## residue renumbering
            '1tud', ## residue renumbering
            '2qwa','2qwb','2qwc','2qwd','2qwe','2qwf','2qwg','2qwh','2qwi','2qwj','2qwk', ## residue renumbering
            '1ow6','1ow7', ## should be dimers according to paper?
            ## tommy errors
##            '1k56','1k57', ## post translational modification of selected chains
            '1ft8', ## too many missing coordinates
            ]

#########
####        ## incorrect transformations
##        self.l_pdbs = ['137l', '180l', '1a4b', '1a4c', '1a8r', '1a9c', '1ab0', '1acm', '1aiz', '1anb', '1anc', '1and', '1aos', '1asu', '1azb', '1azc', '1azn', '1azr', '1b0c', '1b5z', '1b71', '1b7t', '1b8e', '1be7', '1bfs', '1bm9', '1boy', '1bp0', '1bp6', '1bpj', '1bra', '1bsv', '1bsy', '1bws', '1bz5', '1bzd', '1bze', '1c09', '1c3i', '1c76', '1c77', '1c78', '1c79', '1c8t', '1cdp', '1chh', '1chi', '1chj', '1cie', '1cif', '1cih', '1ckg', '1cqj', '1crj', '1csu', '1cz7', '1czb', '1d09', '1d1l', '1d1m', '1d2t', '1dc3', '1dc5', '1dik', '1dvb', '1dvq', '1e2o', '1e3f', '1e5y', '1e5z', '1e67', '1eeq', '1emc', '1eym', '1ezl', '1f13', '1f1b', '1f30', '1f33', '1f41', '1fh2', '1fhn', '1fie', '1fkk', '1fza', '1fze', '1g1o', '1g3k', '1g9f', '1ggo', '1ggt', '1gqw', '1grl', '1h2n', '1hcj', '1hfy', '1hmk', '1hqy', '1ht1', '1ht2', '1huf', '1hzb', '1hzc', '1i2n', '1i2o', '1i2p', '1i2q', '1i2r', '1i5f', '1i5o', '1ib5', '1ibb', '1ibd', '1ibf', '1ibh', '1ihg', '1iin', '1iip', '1ils', '1iro', '1is6', '1iw8', '1izy', '1j0r', '1j1j', '1jde', '1jfa', '1jjw', '1jll', '1jre', '1jts', '1jvo', '1jyb', '1jzo', '1k3z', '1k46', '1k62', '1kaw', '1kbl', '1key', '1kiy', '1kly', '1km0', '1km1', '1km2', '1km3', '1lfa', '1lib', '1lkm', '1lve', '1lyx', '1lzo', '1m3s', '1m8f', '1m8g', '1m8j', '1m8k', '1mm9', '1mou', '1mov', '1moy', '1mp4', '1mp5', '1msp', '1mzf', '1mzy', '1n3t', '1n6m', '1n70', '1ncg', '1nch', '1ni0', '1nni', '1nwc', '1nxo', '1nxs', '1nxv', '1nxx', '1nzr', '1onr', '1os7', '1oy3', '1oyx', '1oz3', '1p2p', '1pv2', '1pvu', '1q95', '1qac', '1qg5', '1qom', '1qvc', '1qve', '1qwh', '1qyb', '1r4c', '1r95', '1rap', '1raq', '1rfq', '1rg0', '1rtv', '1ryt', '1s5g', '1s7y', '1s98', '1sbd', '1sbe', '1sbf', '1scz', '1sok', '1soq', '1svp', '1t4a', '1t6k', '1ta4', '1tfh', '1tij', '1tjd', '1trm', '1tsh', '1tta', '1ttb', '1ttc', '1tth', '1twj', '1txq', '1txy', '1u36', '1u3j', '1u3y', '1u3z', '1u41', '1u42', '1u74', '1umo', '1umr', '1uos', '1uz2', '1v4e', '1v4h', '1v4i', '1v4j', '1v4k', '1v5h', '1vg2', '1vg3', '1vg4', '1vg6', '1vg7', '1viv', '1vlx', '1vse', '1vwi', '1vwj', '1vwk', '1vwl', '1vwq', '1vwr', '1w7s', '1w7t', '1wkz', '1wl0', '1wl1', '1wl2', '1wl3', '1wlc', '1woc', '1x7s', '1x7t', '1xjw', '1xub', '1yai', '1ycc', '1ygz', '1yj4', '1yph', '1zoo', '1zop', '1zv2', '2a3t', '2akq', '2azu', '2bcn', '2bfk', '2bfl', '2blg', '2bqx', '2d6k', '2d6l', '2dik', '2dqr', '2dyu', '2dyv', '2e2k', '2eu1', '2g3x', '2g3z', '2g4e', '2g4g', '2gch', '2gsw', '2hkq', '2huk', '2i0b', '2i0c', '2if9', '2ipr', '2itj', '2iwo', '2iwp', '2ixj', '2j0n', '2jaa', '2lve', '2msp', '2nod', '2noy', '2p4m', '2pab', '2phi', '2q31', '2sba', '2sic', '2snw', '2trh', '2try', '2tsa', '3azu', '3sic', '4azu', '4blc', '4cha', '4cpv', '4lve', '4p2p', '5azu', '5cha', '5cpv', '5cro', '5gch', '5sic', '5ttr', '7cat', '8atc', '8cat']
####        ## viruses
##        self.l_pdbs = [
##            '1aq3', '1aq4', '1bms', '1kuo', '1msc', '1mst', '1mva', '1mvb', '1u1y', '1zdh', '1zdi', '1zdj', '1zdk', '1zse', '2b2d', '2b2e', '2b2g', '2bny', '2bq5', '2bs0', '2bs1', '2bu1', '2c4q', '2c4y', '2c4z', '2c50', '2c51', '2iz8', '2iz9', '2izn', '2ms2', '5msf', '6msf', '7msf',
##            ]
        self.l_pdbs = [
##            '1uyt','1uyv', ## check their rotation in firefox!!!
            ]
        ## different hetero compounds (mostly different saccharides)
        self.l_pdbs = ['101m', '103m', '104m', '105m', '106m', '107l', '107m', '108l', '108m', '109l', '109m', '110l', '110m', '111m', '112l', '112m', '113l', '114l', '115l', '118l', '119l', '11gs', '120l', '122l', '123l', '125l', '126l', '127l', '128l', '129l', '130l', '131l', '138l', '139l', '140l', '141l', '142l', '143l', '144l', '147l', '155l', '156l', '158l', '159l', '160l', '161l', '162l', '163l', '164l', '165l', '166l', '190l', '195l', '196l', '197l', '198l', '199l', '1a4f', '1a6g', '1a6m', '1abs', '1aes', '1aj9', '1ajg', '1ajh', '1anf', '1axm', '1azx', '1b3v', '1b3w', '1b3y', '1b3z', '1bbb', '1bfb', '1bfc', '1bfn', '1bos', '1byb', '1byc', '1bzr', '1c4q', '1c60', '1c61', '1c63', '1c64', '1c67', '1c6d', '1c6g', '1c6j', '1c6l', '1c6m', '1c6p', '1c6q', '1ccg', '1cqe', '1cqf', '1crl', '1cvn', '1cxk', '1d1i', '1d1k', '1d3c', '1d6p', '1d6q', '1do1', '1do3', '1do4', '1do7', '1ds4', '1dse', '1dsg', '1dso', '1dsp', '1dtu', '1dwr', '1dws', '1dwt', '1dxc', '1dxd', '1e4o', '1e5v', '1e60', '1e61', '1ebc', '1ehh', '1ehn', '1eib', '1en2', '1enm', '1eo5', '1eo7', '1eqh', '1ez9', '1f2q', '1f65', '1f9d', '1fae', '1faw', '1fbw', '1ffr', '1fqa', '1fqb', '1fqc', '1fqd', '1fsx', '1fzf', '1fzg', '1fzk', '1fzo', '1g08', '1g09', '1g0a', '1g1v', '1g1w', '1g94', '1g9f', '1g9h', '1gim', '1gin', '1gmn', '1gmo', '1gwl', '1h12', '1h1x', '1h3t', '1h3u', '1hcn', '1hco', '1hew', '1hrp', '1hvq', '1hzv', '1i3h', '1i82', '1i8a', '1iow', '1irc', '1isv', '1isw', '1isx', '1itc', '1iuc', '1iw6', '1iw9', '1ix4', '1ixf', '1j0i', '1j0j', '1j0k', '1j0y', '1j0z', '1j18', '1j3y', '1j3z', '1j40', '1j41', '1j7y', '1j88', '1j89', '1jin', '1jip', '1js4', '1jw5', '1jw8', '1k1k', '1k9t', '1kck', '1kni', '1kqz', '1kr0', '1kr1', '1ksz', '1l1y', '1l2a', '1l55', '1l59', '1l5w', '1l61', '1l62', '1l63', '1l66', '1l67', '1l6i', '1l72', '1l76', '1l77', '1l79', '1l80', '1l81', '1l82', '1l83', '1l84', '1l85', '1l86', '1l87', '1l88', '1l89', '1l90', '1l91', '1l92', '1l93', '1l94', '1l95', '1la6', '1lax', '1lgu', '1ljw', '1ljy', '1llh', '1lmp', '1lmq', '1loh', '1ltw', '1lxk', '1lye', '1lyf', '1lyg', '1lyh', '1lyi', '1lyj', '1lzb', '1lzc', '1lze', '1lzg', '1m18', '1m19', '1m1a', '1m34', '1m9p', '1mbc', '1mbi', '1mbo', '1mcy', '1md9', '1mdb', '1mfb', '1mfc', '1mfv', '1mko', '1mlf', '1mlj', '1mlm', '1mlq', '1mlu', '1moc', '1mpm', '1mpn', '1mpo', '1mqm', '1mqn', '1mwe', '1mym', '1myz', '1n2c', '1n7r', '1nca', '1ncb', '1ncc', '1nej', '1nh6', '1nht', '1nkx', '1nm9', '1nno', '1np1', '1nq9', '1nv0', '1nv1', '1nv2', '1nv3', '1nv4', '1nv5', '1nv6', '1nws', '1nwt', '1nwu', '1ofz', '1ona', '1ot1', '1ot2', '1our', '1ovs', '1ozl', '1ozw', '1p2y', '1p4m', '1p7r', '1pbx', '1pez', '1q6c', '1q6d', '1q6e', '1q6f', '1q6g', '1q8v', '1q9q', '1q9r', '1q9s', '1q9t', '1qdo', '1qnw', '1qoo', '1qos', '1qt4', '1qt6', '1qt7', '1qtv', '1re2', '1rey', '1rez', '1rp9', '1rpk', '1rvw', '1sbd', '1sbe', '1sbf', '1sf4', '1sf6', '1sf7', '1sfb', '1sfg', '1sl4', '1sla', '1slb', '1spe', '1sqy', '1sr0', '1syt', '1t8g', '1tk3', '1trh', '1u18', '1u7r', '1u8e', '1ubz', '1ucq', '1uh4', '1uih', '1uks', '1ukt', '1ulx', '1ur1', '1ur2', '1usw', '1uu4', '1uu6', '1uy2', '1uy3', '1uy4', '1uyx', '1uyy', '1uza', '1v3h', '1v3i', '1v3l', '1v3m', '1v6w', '1v6x', '1vb9', '1vbo', '1vbp', '1vem', '1ven', '1veo', '1vep', '1vfk', '1vfm', '1vfo', '1vfu', '1vxc', '1vxf', '1w2u', '1w3k', '1w3l', '1w8u', '1wd4', '1wds', '1wu5', '1wu6', '1x0k', '1x0s', '1xk0', '1xk3', '1y1e', '1y1i', '1yff', '1yvq', '1z32', '1zbk', '1zbv', '1zbw', '1zlu', '1zlv', '1zlw', '1zu8', '200l', '206l', '217l', '220l', '221l', '222l', '224l', '227l', '228l', '230l', '232l', '233l', '234l', '253l', '254l', '255l', '2aft', '2alu', '2ar6', '2arb', '2av6', '2aw3', '2axm', '2b31', '2b4f', '2b6t', '2b6w', '2b6x', '2b6y', '2b6z', '2b70', '2b72', '2b73', '2b74', '2b75', '2bjj', '2blj', '2bwa', '2bwc', '2bwh', '2c3w', '2c3x', '2c4a', '2c4l', '2c6p', '2civ', '2ciw', '2cnc', '2d0f', '2d0g', '2d0h', '2d22', '2d24', '2d3l', '2d3n', '2d43', '2dij', '2dln', '2dn3', '2dp8', '2dpe', '2ds9', '2dsu', '2dsv', '2dsw', '2dsz', '2dt0', '2dt1', '2dt3', '2dwh', '2f9n', '2f9p', '2fa7', '2g0r', '2g0s', '2g41', '2gmm', '2gn7', '2gnb', '2gnm', '2hco', '2hyu', '2hyv', '2it6', '2j18', '2j19', '2j1t', '2j1v', '2j72', '2j73', '2jaf', '2jag', '2l78', '2man', '2mga', '2mgc', '2mgf', '2mgk', '2mgm', '2mya', '2myb', '2myc', '2myd', '2mye', '2nmo', '2nn8', '2nwj', '2o1l', '2o92', '2oe4', '2oe7', '2oe9', '2oea', '2oot', '2sba', '2spl', '2spn', '2uvh', '2uvj', '3cel', '3gss', '3man', '3mbp', '3tf4', '4mbp', '4tf4', '5cel', '6cel', '7api', '7cel', '8api', '9api']
        self.l_pdbs = ['2c4a','2c4l','1mwe']
                       
        ## removal
        try:
            ## remove pdbs spanning multiple pdb files/IDs
            self.l_pdbs.remove('1p0t'); self.l_pdbs.remove('1otz')
            self.l_pdbs.remove('2j00'); self.l_pdbs.remove('2jo1'); self.l_pdbs.remove('2j02'); self.l_pdbs.remove('2jo3')
        except:
            None
        for pdb in errorpdbs:
            try:
                self.l_pdbs.remove(pdb)
            except:
                None

        ## sort and count after removal
        self.l_pdbs.sort()
##        self.l_pdbs.reverse()
        self.pdbcount = len(self.l_pdbs)

        d_seq, d_chains_intrapdb_sequence_identical = self.parse_sequences()
        d_rmsd = self.analyze_sequences(d_seq, d_chains_intrapdb_sequence_identical)
##        d_rmsd = read_rmsd_from_file()
        self.write_rmsd_to_file(d_rmsd, d_seq, prefix='rmsd')
##########
        

#### out
##        sets_pdbs = []
##        self.l_pdbs = set()
##
##        dirs = os.listdir('/oxygenase_local/tc/quakes/pdb/')
##        for dir in dirs:
##            files = os.listdir('/oxygenase_local/tc/quakes/pdb/%s' %(dir))
##            for file in files:
##
##                pdb1 = file[0:4]
##                pdb2 = file[6:10]
##                added = False
##                for i in range(len(sets_pdbs)):
##                    set_pdb = sets_pdbs[i]
##                    if pdb1 in set_pdb or pdb2 in set_pdb:
##                        sets_pdbs[i] |= set([pdb1,pdb2])
##                        added = True
##                if added == False:
##                    sets_pdbs += [set([pdb1,pdb2])]
##                self.l_pdbs |= set([pdb1]) ##
##                self.l_pdbs |= set([pdb2]) ##
##
##        self.l_pdbs = list(self.l_pdbs)
##        for set_pdb in sets_pdbs:
##            self.l_pdbs = list(set_pdb)
##            for pdb in errorpdbs:
##                try:
##                    self.l_pdbs.remove(pdb)
##                except:
##                    None
##            self.pdbcount = len(self.l_pdbs)
##            d_seq, d_chains_intrapdb_sequence_identical = self.parse_sequences()
##            d_rmsd = self.analyze_sequences(d_seq, d_chains_intrapdb_sequence_identical)
##            self.write_rmsd_to_file(d_rmsd, d_seq, prefix='rmsd')

## clusters95

        import os
        fd = open('clusters95.txt','r')
        lines = fd.readlines()
        fd.close()
        clusters = {}
        for line in lines:
            cluster = int(line.split()[0])
            pdb = line.split()[2][:4].lower()
            if cluster not in clusters.keys():
                clusters[cluster] = []
            clusters[cluster] += [pdb]
        l_clusters = clusters.keys()
        l_clusters.sort()
        for cluster in l_clusters:
            if cluster not in range(9950,10000):
                continue
            print '-------cluster--------', cluster
            self.l_pdbs = list(set(clusters[cluster]))
            pdbs = self.l_pdbs
            pdbs.sort()
            for pdb in pdbs:
                if not os.path.isfile('/oxygenase_local/data/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb)):
                    self.l_pdbs.remove(pdb)
            for pdb in errorpdbs:
                try:
                    self.l_pdbs.remove(pdb)
                except:
                    None
            if len(self.l_pdbs) == 1:
                continue
            self.pdbcount = len(self.l_pdbs)
            d_seq, d_chains_intrapdb_sequence_identical = self.parse_sequences()
            d_rmsd = self.analyze_sequences(d_seq, d_chains_intrapdb_sequence_identical)
            self.write_rmsd_to_file(d_rmsd, d_seq, prefix='rmsd')
##########

##        self.analyze_rmsd()

        return


    def analyze_rmsd(self):

        print 'analyze rmsd'

        import os

        n_lines = 1+len(self.l_columns_html)+1

        htmls = os.listdir('htm/')
        lines = []
        for html in htmls:
            fd = open('htm/%s' %(html),'r')
            lines += fd.readlines()[1+n_lines:-1]
            fd.close()

        d_data_ratio_discrete = {
            'mutations':{},
            'chains':{},
            'residues':{},
            'coordinates':{},
            }
        d_data_ratio_continuous = {
            'pH':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
            'T':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
            'res':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
            }
        d_data_nominal = {
            'spacegroup':{'identical':{},'different':{}},
            'hetIDs':{'identical':[],'different':[]},
            'REMARK465':{'True':[],'False':[]},
            'REMARK470':{'True':[],'False':[]},
            'transformations':{'True':[],'False':[]},
            }

        ##
        ## parse data
        ##
        ## loop over table rows
        for i in range(len(lines)/(n_lines)):
            ## parse pdbs
            pdb1_line = lines[i*n_lines+3]
            pdb2_line = lines[i*n_lines+4]
            pdbs = [pdb1_line,pdb2_line]
            for j in range(2):
                pdb_line = pdbs[j]
                pdb_index2 = pdb_line.index('</a>')
                pdb_index1 = pdb_line[:pdb_index2].rindex('>')+1
                pdb = pdb_line[pdb_index1:pdb_index2]
                pdbs[j] = pdb
            pdb1 = pdbs[0]
            pdb2 = pdbs[1]

            ##
            ## parse non-hyperlinkreferenced data
            ##
            d_parameters = {
                'rmsd':7,
                'mutations':8,
                'chains':9,
                'residues':10,
                'coordinates':11,
                'pH1':12,
                'pH2':13,
                'T1':14,
                'T2':15,
                'res1':16,
                'res2':17,
                'spacegroup1':18,
                'spacegroup2':19,
                'REMARK465':20,
                'REMARK470':21,
                'transformations':22,
                }
            for parameter in d_parameters.keys():
                i_line = d_parameters[parameter]
                s_line = lines[i*n_lines+i_line]
                index2 = s_line.rindex('</td>')
                index1 = s_line[:index2].rindex('>')+1
                value = s_line[index1:index2].strip()
                d_parameters[parameter] = value
            if float(d_parameters['residues']) < 50:
                errorpdbs = set(['2bfk','2bfl'])
                if len(set([pdb1,pdb2])-errorpdbs) > 0:
                    print set([pdb1,pdb2])-errorpdbs
                    print pdb1, pdb2
                    stop1 ## temp!!!

            ##
            ## parse hyperlinkreferenced data
            ##
            d_hetIDs = {
                'hetIDs1':25,
                'hetIDs2':26,
                }
            
            for hetIDs in d_hetIDs.keys():
                i_line = d_hetIDs[hetIDs]
                s_line = lines[i*n_lines+i_line]
                d_hetIDs[hetIDs] = []
                index3 = s_line.index('</td>')
                if s_line[index3-4:index3] == '</a>':
                    index2 = 0
                    while index2 != index3-4:
                        index2 += 1
                        index2 += s_line[index2:].index('</a>')
                        index1 = s_line[:index2].rindex('>')+1
                        hetID = s_line[index1:index2].strip()
                        d_hetIDs[hetIDs] += [hetID]

            ## set rmsd
            rmsd = d_parameters['rmsd']
            if float(rmsd) > 5.:
                continue

            ##
            ## continuous ratio scale data
            ##
            for parameter in d_data_ratio_continuous.keys():
                values = []
                for no in ['1','2']:
                    value = d_parameters[parameter+no]
                    if value not in ['NULL','N/A']:
                        if ';' in value:
                            errorpdbs = [ ## temp!!!
                                '1znb','1zlx','1i3h','1j5o','1f33','1qvc',
                                '1ypu','1uc7',
                                ]
                            if (no == '1' and pdb1 in errorpdbs) or (no == '2' and pdb2 in errorpdbs):
                                value = value.split(';')[-1]
                            else:
                                l_values = value.split(';')
                                set_values = set()
                                for value in l_values:
                                    set_values |= set([value.strip()])
                                    if len(set_values) > 1:
                                        print pdb1, pdb2, no
                                        print d_parameters['T'+no]
                                        print d_parameters['pH'+no]
                                        print set_values
                                        stop
                                    else:
                                        value = list(set_values)[0]
                        value = float(value)
                        if (value < 50 or value > 300) and parameter == 'T':
                            errorpdbs = set([
                                '1f1g','1f18','1abs','1dxd','1dvb',
                                '1dmm','1dmn','1dmq',
                                '195l','196l','200l','197l','199l','198l',
                                ## unknown
                                '1za1','1qtv','1khh',
                                ## correct (cryogen)
                                '2fbb','1j3y','1j41','1ajg','1ajh','1a6k',
                                ])
                            if len(set([pdb1,pdb2])-errorpdbs) == 2:
                                print value, no
                                print pdb1, pdb2
                                stop2 ## temp!!!
                        if (value < 1 or value > 12) and parameter == 'pH':
                            errorpdbs = set(['1fgn'])
                            if len(set([pdb1,pdb2])-errorpdbs) == 2:
                                print pdb1, pdb2
                                print d_parameters['pH1'], d_parameters['pH2']
                                stop3 ## temp!!!
                        if value > 11.5 and parameter == 'res':
                            print pdb1, pdb2
                            print d_parameters['res1'], d_parameters['res2']
                            stop4 ## temp!!
                        values += [value]
                        d_data_ratio_continuous[parameter]['single'] += [[value,rmsd]]
                if len(values) == 2:
                    d_data_ratio_continuous[parameter]['max'] += [[max(values),rmsd]]
                    d_data_ratio_continuous[parameter]['min'] += [[min(values),rmsd]]
                    d_data_ratio_continuous[parameter]['average'] += [[sum(values)/2.,rmsd]]
                    d_data_ratio_continuous[parameter]['difference'] += [[abs(values[1]-values[0]),rmsd]]

            ##
            ## discrete ratio scale data
            ##
            for parameter in d_data_ratio_discrete.keys():
                value = d_parameters[parameter]
                if value not in d_data_ratio_discrete[parameter].keys():
                    d_data_ratio_discrete[parameter][value] = []
                d_data_ratio_discrete[parameter][value] += [rmsd]

            ######
            ## nominal scale data
            ######
            ##
            ## spacegroup
            ##
            spacegroup1 = d_parameters['spacegroup1']
            spacegroup2 = d_parameters['spacegroup2']
            if spacegroup1 == spacegroup2:
                if not spacegroup1 in d_data_nominal['spacegroup']['identical'].keys():
                    d_data_nominal['spacegroup']['identical'][spacegroup1] = []
                d_data_nominal['spacegroup']['identical'][spacegroup1] += [rmsd]
            else:
                if not spacegroup1 in d_data_nominal['spacegroup']['different'].keys():
                    d_data_nominal['spacegroup']['different'][spacegroup1] = []
                d_data_nominal['spacegroup']['different'][spacegroup1] += [rmsd]
                if not spacegroup2 in d_data_nominal['spacegroup']['different'].keys():
                    d_data_nominal['spacegroup']['different'][spacegroup2] = []
                d_data_nominal['spacegroup']['different'][spacegroup2] += [rmsd]
            ##
            ## hetIDs
            ##
            difference_hetID = abs(len(d_hetIDs['hetIDs2'])-len(d_hetIDs['hetIDs1']))
            if difference_hetID > 0:
                d_data_nominal['hetIDs']['different'] += [rmsd]
            else:
                d_data_nominal['hetIDs']['identical'] += [rmsd]
            ##
            ## remarks, transformations
            ##
            for parameter in ['REMARK465','REMARK470','transformations']:
                if d_parameters[parameter] == 'True':
                    d_data_nominal[parameter]['True'] += [rmsd]
                elif d_parameters[parameter] == 'False':
                    d_data_nominal[parameter]['False'] += [rmsd]
               

        ## continuous ratio scale data
        for parameter in d_data_ratio_continuous.keys():
            for type in d_data_ratio_continuous[parameter]:
                gnuplotdata = ''
                for values in d_data_ratio_continuous[parameter][type]:
                    value = values[0]
                    rmsd = values[1]
                    gnuplotdata += '%s %s\n' %(value, rmsd)
                prefix_gnuplot = parameter+type
                fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
                fd.write(gnuplotdata)
                fd.close()
                self.gnuplot_plot(prefix_gnuplot)

        ## discrete ratio scale data
        for parameter in d_data_ratio_discrete.keys():

            gnuplotdata = ''
            for value in d_data_ratio_discrete[parameter].keys():
                for rmsd in d_data_ratio_discrete[parameter][value]:
                    gnuplotdata += '%s %s\n' %(value, rmsd)

            prefix_gnuplot = parameter
            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata)
            fd.close()

            if parameter in ['chains','residues','coordinates']:
                logarithmic = True
            else:
                logarithmic = False
            self.gnuplot_plot(prefix_gnuplot, logarithmic=logarithmic)

        ####
        ## nominal scale data
        ####
        ##
        ## hetIDs
        ##
        l1 = d_data_nominal['hetIDs']['identical']
        l2 = d_data_nominal['hetIDs']['different']
        gnuplotdata = ''
        for rmsd in l1:
            gnuplotdata += '1 %s\n' %(rmsd)
        for rmsd in l2:
            gnuplotdata += '2 %s\n' %(rmsd)
        prefix_gnuplot = 'hetIDs'
        fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
        fd.write(gnuplotdata)
        fd.close()
        d_xtics = {'identical':1.,'different':2.}
        self.gnuplot_plot(prefix_gnuplot, d_xtics=d_xtics)
        print prefix_gnuplot
        self.twosamplettest(l1,l2)
        ##
        ## remarks, transformations
        ##
        for parameter in ['REMARK465','REMARK470','transformations']:
            l1 = d_data_nominal[parameter]['True']
            l2 = d_data_nominal[parameter]['False']
            gnuplotdata = ''
            for rmsd in l1:
                gnuplotdata += '1 %s\n' %(rmsd)
            for rmsd in l2:
                gnuplotdata += '2 %s\n' %(rmsd)
            prefix_gnuplot = parameter
            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata)
            fd.close()
            d_xtics = {'True':1.,'False':2.}
            self.gnuplot_plot(prefix_gnuplot, d_xtics=d_xtics)
            print prefix_gnuplot
            self.twosamplettest(l1,l2)
        ##
        ## space groups
        ##
        ## xtics
        set_spacegroups = set()
        set_spacegroups |= set(d_data_nominal['spacegroup']['identical'].keys())
        set_spacegroups |= set(d_data_nominal['spacegroup']['different'].keys())
        l_spacegroups = list(set_spacegroups)
        l_spacegroups.sort()
        d_spacegroups = {}
        for i in range(len(l_spacegroups)):
            spacegroup = l_spacegroups[i]
            d_spacegroups[spacegroup] = float(i+1)
        ## gnuplotdata
        d_gnuplotdata = {'identical':'','different':''}
        for type in d_gnuplotdata.keys():
            for spacegroup in d_data_nominal['spacegroup'][type].keys():
                xval = d_spacegroups[spacegroup]
                for rmsd in d_data_nominal['spacegroup'][type][spacegroup]:
                    d_gnuplotdata[type] += '%s %s\n' %(xval, rmsd)
            prefix_gnuplot = 'spacegroups%s' %(type)
            gnuplotdata = d_gnuplotdata[type]
            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata)
            fd.close()
            self.gnuplot_plot(prefix_gnuplot, d_xtics=d_spacegroups, xlabel='spacegroups')

        return


    def twosamplettest(self,l1,l2):

        import math

        d_statistics = {
            1:{'sample':l1},
            2:{'sample':l2},
            }

        for i_sample in d_statistics.keys():
            l_sample = d_statistics[i_sample]['sample']
            n = len(l_sample)
            sumx = 0.
            sumxx = 0.
            for i in range(len(l_sample)):
                x = float(l_sample[i])
                sumx += x
                sumxx += x**2
            SS = sumxx-sumx/n
            mean = sumx/n
            d_statistics[i_sample]['mean'] = mean
            d_statistics[i_sample]['SS'] = SS
            d_statistics[i_sample]['n'] = n

        n1 = d_statistics[1]['n']
        n2 = d_statistics[2]['n']
        ss1 = d_statistics[1]['SS']
        ss2 = d_statistics[2]['SS']
        mean1 = d_statistics[1]['mean']
        mean2 = d_statistics[2]['mean']

        var_pooled = (ss1+ss2)/(n1+n2-2)
        stderr = stddev = math.sqrt(var_pooled/n1+var_pooled/n2)
        t = ((mean1-mean2)-mean)/stddev
        p = self.tdist(t,n1+n2-2)

        print 'n', n1, n2
        print 'var', ss1/n1, ss2/n2, 'pooled', var_pooled
        print 'mean', mean1, mean2
        print 't', t
        print 'p', p

        return p


    def gnuplot_plot(
        self,prefix, d_xtics=None, logarithmic=False, xlabel=''
        ):

        ylabel = 'RMSD'

        gnuplotsettings = []
        gnuplotsettings += [
            'set terminal postscript eps enhanced color "Helvetica" 48\n',
            'set output "ps/%s.ps"\n' %(prefix),
            'set size 4,4\n', ## scale 400%
            'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
##            'set encoding iso_8859_1\n', ## postscript encoding for special characters
##            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
            'set xlabel "%s"\n' %(xlabel),
            'set ylabel "%s"\n' %(ylabel),
        ]
        if logarithmic == True:
            gnuplotsettings += ['set logscale x\n']
        if d_xtics:
            line_xtic = 'set xtics ('
            for xtic in d_xtics.keys():
                line_xtic += '"%s" %s, ' %(xtic, d_xtics[xtic])
            line_xtic = line_xtic[:-2]+')\n'
            gnuplotsettings += [
                line_xtic,
                'set xtics rotate\n',
            ]
        gnuplotsettings += [
            'plot "%s.gnuplotdata" lt 0 ps 2 pt 2\n' %(prefix),
        ]
        fd = open('%s.gnuplotsettings' %(prefix),'w')
        fd.writelines(gnuplotsettings)
        fd.close()

        os.system('/software/bin/gnuplot %s.gnuplotsettings' %(prefix))
        os.remove('%s.gnuplotdata' %(prefix))
        os.remove('%s.gnuplotsettings' %(prefix))

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
        d_pdb = {}
        d_hetero = {}

        ##
        ## loop 1 over pdbs
        ##
        t1 = time.clock()
        for i1 in range(self.pdbcount-1):

            self.pdb1 = pdb1 = self.l_pdbs[i1]

            skippdb = self.pdbskip(d_seq, pdb1)
            if skippdb == True:
                continue

            ## identify biomolecule(s)
            d_biomolecules1 = self.identify_biomolecule(pdb1, d_seq)

            ## reset dictionary of coordinates to save memory
            d_pdb = {}

            ##
            ## loop 2 over pdbs
            ##
            for i2 in range(i1+1,self.pdbcount):

                self.pdb2 = pdb2 = self.l_pdbs[i2]

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
                    for biomolecule2 in d_biomolecules2.keys():

                        bmchains2 = d_biomolecules2[biomolecule2]['chains']
                        bmpolymercount2 = d_biomolecules2[biomolecule2]['polymercount']

                        ## skip if different number of chains in the biomolecule
                        if bmpolymercount1 != bmpolymercount2:
                            continue

                        ## skip if different hetero compounds
                        different = self.different_hetero_compounds(pdb1,pdb2,d_seq)
                        if different == True:
                            continue

                        ## skip if different polymers (other than long peptides)
                        different = self.different_nucleotides_saccharides_shortpeptides(pdb1,pdb2,d_seq)
                        if different == True:
                            continue

                        ##
                        ## identify sequence similar chains between pdbs (long peptides only)
                        ##
                        d_chains_interpdb_sequence_similar = self.identify_similar_chains_from_sequence_inter(
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
                            print biomolecule1, bmSEQRESchains1_not_similar_to_SEQRESchains2
                            print biomolecule2, bmSEQRESchains2_not_similar_to_SEQRESchains1

                            ## 1) check that non sequence similar chains (if any) are not long peptides
                            d_bmSEQRESchains_not_similar_to_SEQRESchains = {
                                pdb1:{'chains1':bmSEQRESchains1_not_similar_to_SEQRESchains2,'chains2':bmSEQRESchains2_not_similar_to_SEQRESchains1,'pdb2':pdb2},
                                pdb2:{'chains1':bmSEQRESchains2_not_similar_to_SEQRESchains1,'chains2':bmSEQRESchains1_not_similar_to_SEQRESchains2,'pdb2':pdb1},
                                }
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
                                            print pdb, chain
                                            break
                                    if peptide == True and long == True:
                                        break
                            ## continue if long peptide
                            if peptide == True and long == True:
                                continue


                        ##
                        ## identify equivalent chains (interpdb) from structure
                        ##

                        ## parse coordinates
                        if pdb1 not in d_pdb.keys():
                            d_pdb,d_hetero = self.parse_coordinates(pdb1, d_pdb, d_seq[pdb1],d_hetero)
                        d_pdb,d_hetero = self.parse_coordinates(pdb2, d_pdb, d_seq[pdb2],d_hetero)

                        ## skip if different hetero compounds (2nd check)
                        if d_hetero[pdb1] != d_hetero[pdb2]:
                            print d_hetero[pdb1]
                            print d_hetero[pdb2]
                            different = False
                            pdbpairs = [[pdb1,pdb2],[pdb2,pdb1]]
                            for pdbpair in pdbpairs:
                                pdba = pdbpair[0]
                                pdbb = pdbpair[1]
                                for compound in d_hetero[pdba]:
                                    if compound not in d_hetero[pdbb]:
                                        print compound
                                        different = True
                                        fd = open('CONECT.txt','a')
                                        fd.write('%s %s %s %s %s\n%s\n%s\n' %(pdb1, pdb2, compound, d_hetero[pdb1], d_hetero[pdb2], d_seq[pdb1]['TITLE'], d_seq[pdb2]['TITLE']))
                                        fd.close()
                                        break
                                if different == True:
                                    break
                            if different == True:
                                continue

                        ## identify equivalent chains and calculate rmsd
                        l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations = self.identify_interpdb_equivalent_chains_from_structure(
                            pdb1, pdb2,
                            d_chains_intrapdb_sequence_identical,
                            d_chains_interpdb_sequence_similar,
                            d_pdb, d_seq,
                            biomolecule1, biomolecule2,
                            d_biomolecules1, d_biomolecules2,
                            )

                        n_mutations = self.count_mutations(pdb1, pdb2, l_equivalent_chains, d_chains_interpdb_sequence_similar, d_chains_intrapdb_sequence_identical)

##                        ## calculate rmsd in spheres around the mutation site
##                        if mutations == 1 and n_chains == 1:
##                            rmsd4, rmsd8, rmsd16, rmsd32 = self.spherermsd(
##                                pdb1, pdb2,
##                                d_seq, d_pdb,
##                                l_equivalent_chains,
##                                d_chains_interpdb_sequence_similar,
##                                tv1, rm, tv2,
##                                )
##                        else:
##                            rmsd4 = rmsd8 = rmsd16 = rmsd32 = 'N/A'


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
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['residues'] = n_residues
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['coordinates'] = n_coordinates
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['l_equivalent_chains'] = l_equivalent_chains
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['mutations'] = n_mutations
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['transformations'] = transformations
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd4'] = rmsd4
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd8'] = rmsd8
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd16'] = rmsd16
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd32'] = rmsd32

                        ## write data to file
                        d_quickrmsd = {pdb1:{biomolecule1:{pdb2:{biomolecule2:d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]}}}}
                        self.write_rmsd_to_file(d_quickrmsd, d_seq, prefix='quickrmsd')

                        d_biomolecules = {
                            pdb1:{'biomolecule':biomolecule1},
                            pdb2:{'biomolecule':biomolecule2},
                            }

                        ## color code structure by rmsd
                        self.rmsd2bfactor(pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_pdb, d_seq, tv1, rm, tv2, l_equivalent_chains, bmchains1, bmchains2)

                        del d_pdb[pdb2]

        return d_rmsd_identical


    def count_mutations(self, pdb1, pdb2, l_equivalent_chains, d_chains_interpdb_sequence_similar, d_chains_intrapdb_sequence_identical):

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
                for rep_chain in d_chains_intrapdb_sequence_identical[pdb]:
                    if chain == rep_chain or chain in d_chains_intrapdb_sequence_identical[pdb][rep_chain]:
                        d_chains[pdb] = rep_chain
                        break
            rep_chain1 = d_chains[pdb1]
            rep_chain2 = d_chains[pdb2]
            n_mutations_per_chain = len(d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l_mutations'])
            n_mutations += n_mutations_per_chain

        return n_mutations


    def spherermsd(
        self,
        pdb1, pdb2, ## pdbs
        d_seq, ## sequences
        d_pdb, ## coordinates
        l_equivalent_chains, ## equivalent chains
        d_chains_interpdb_sequence_similar, ## mutations
        tv1, rm, tv2, ## transformation
        ):

        import sys
        sys.path.append('/home/people/tc/python/Protool/')
        import geometry
        instance_geometry = geometry.geometry()

        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']

            ATOMseq1,d_res_nos1 = self.ATOM2seq(d_pdb[pdb1], rep_chain1)
            ATOMseq2,d_res_nos2 = self.ATOM2seq(d_pdb[pdb2], rep_chain2)
            l1 = d_chains_interpdb_sequence_similar[rep_chain1]['l1']
            l2 = d_chains_interpdb_sequence_similar[rep_chain1]['l2']

            if l1 > 0 or l2 > 0: ## e.g. 1ftg.pdb,1dx9.pdb
                print ATOMseq1
                print ATOMseq2
                s1 = d_chains_interpdb_sequence_similar[rep_chain1]['s1']
                s2 = d_chains_interpdb_sequence_similar[rep_chain1]['s2']
                print s1
                print s2
                print pdb1, pdb2
                expected

            (
                coordinates1, coordinates2, rescount,
                ) = self.ATOMrecords2coordinates(
                    d_pdb, pdb1, pdb2, rep_chain1, rep_chain2, d_res_nos1, d_res_nos2,
                    l1, l2, len(ATOMseq1), tv1=tv1, rm=rm, tv2=tv2
                    )

            l_mutations = d_chains_interpdb_sequence_similar[rep_chain1]['l_mutations']
            for mutation in l_mutations:
                res_no1 = d_res_nos1[mutation[0]]['res_no']
                res_no2 = d_res_nos2[mutation[0]]['res_no']
                iCode1 = d_res_nos1[mutation[0]]['iCode']
                iCode2 = d_res_nos2[mutation[0]]['iCode']
                hypocenter1 = d_pdb[pdb1]['chains'][rep_chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms']['CA']['coordinate']
                hypocenter2 = d_pdb[pdb2]['chains'][rep_chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms']['CA']['coordinate']

                d_rmsds = {4:0,8:0,16:0,32:0}
                for dist in d_rmsds.keys():

                    sqdist = dist**2
                    sphere_coordinates1 = []
                    sphere_coordinates2 = []

                    for i in range(len(coordinates1)):

                        coordinate1 = coordinates1[i]
                        coordinate2 = coordinates2[i]

                        sqdist1 = sum((hypocenter1-coordinate1)**2)
                        sqdist2 = sum((hypocenter2-coordinate2)**2)

                        if min(sqdist1,sqdist2) < sqdist:

                            sphere_coordinates1 += [coordinate1]
                            sphere_coordinates2 += [coordinate2]

                    rmsd = instance_geometry.superpose(
                        sphere_coordinates1,sphere_coordinates2
                        )

                    print rep_chain1, mutation, dist, rmsd, float(len(sphere_coordinates1))/len(coordinates1)

                    d_rmsds[dist] = rmsd

        rmsd4 = d_rmsds[4]
        rmsd8 = d_rmsds[8]
        rmsd16 = d_rmsds[16]
        rmsd32 = d_rmsds[32]

        return rmsd4, rmsd8, rmsd16, rmsd32


    def different_nucleotides_saccharides_shortpeptides(self,pdb1,pdb2,d_seq):

##        print 'identifying different polymers (not long peptides)'

        d_chains = {
            pdb1:{'peptide':set(),'saccharide':set(),'nucleotide':set()},
            pdb2:{'peptide':set(),'saccharide':set(),'nucleotide':set()},
            }
        for pdb in d_chains.keys():
            for chain in d_seq[pdb]['chains']:
                chain_type = d_seq[pdb]['chains'][chain]['type']
                chain_seq = d_seq[pdb]['chains'][chain]['seq']
                chain_len = len(chain_seq)
                if chain_type != 'peptide' or chain_len < self.min_len_chain:
                    d_chains[pdb][chain_type] |= set([chain_seq])
        for polymer in ['peptide','saccharide','nucleotide']:
            if len(d_chains[pdb1][polymer] ^ d_chains[pdb2][polymer]) > 0:
                return True
        
        return False


    def different_hetero_compounds(self,pdb1,pdb2,d_seq):

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
            d_hetIDs[pdb] -= set(self.d_modres.keys()) ## modified residues
        ## compare hetero compounds assuming the following two statements to be correct
        ## 1) "A particular HET group is represented in the PDB archives with a *unique* hetID."
        ## 2) Depositors specify *all* hetero atoms observed in the electron density map.
        ## ignore ions and ignore saccharides
        if (
            d_hetIDs[pdb1]-set(self.d_ions.keys())-set(self.d_saccharides.keys())-set(self.d_stereoisomers.keys())
            !=
            d_hetIDs[pdb2]-set(self.d_ions.keys())-set(self.d_saccharides.keys())-set(self.d_stereoisomers.keys())
            ):
            return True

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
            return True

        ## skip if different heteroatoms (for which stereoisomers exist)
## change stereoisomers earlier and skip this step!!!
        stereoisomers1_stereo = set(self.d_stereoisomers.keys()) & d_hetIDs[pdb1]
        stereoisomers2_stereo = set(self.d_stereoisomers.keys()) & d_hetIDs[pdb2]
        stereoisomers1 = set()
        stereoisomers2 = set()
        for stereoisomer in stereoisomers1_stereo:
            stereoisomers1 |= set([self.d_stereoisomers[stereoisomer]])
        for stereoisomer in stereoisomers2_stereo:
            stereoisomers2 |= set([self.d_stereoisomers[stereoisomer]])
        if stereoisomers1 != stereoisomers2:
            return True

        ## skip if different charges (e.g. 1ext.pdb, 1ncf.pdb but also 1ncy,1ncz)
        ions1 = set(self.d_ions.keys()) & d_hetIDs[pdb1]
        ions2 = set(self.d_ions.keys()) & d_hetIDs[pdb2]
        plus1 = 0
        minus1 = 0
        for ion in ions1:
            formula = self.d_ions[ion][0]
            if len(formula.split()) == 1:
                charge = self.d_ions[ion][1]
                if charge > 0:
                    plus1 += charge
                if charge < 0:
                    minus1 += charge
        plus2 = 0
        minus2 = 0
        for ion in ions2:
            formula = self.d_ions[ion][0]
            if len(formula.split()) == 1:
                charge = self.d_ions[ion][1]
                if charge > 0:
                    plus2 += charge
                if charge < 0:
                    minus2 += charge
        if plus1 != plus2 or minus1 != minus2:
            return True

        return False


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
        if d_seq[pdb]['EXPDTA'] != 'X-RAY':
            pdbskip = True
            return pdbskip

        ## continue if multiple models (NMR or non-NCS-averaged x-ray structures - e.g. 1hto.pdb)
        if 'MODEL' in d_seq[pdb].keys():
            pdbskip = True
            return pdbskip

##        ## continue if low resolution (e.g. 1jgp.pdb)
##        if d_seq[pdb]['REMARK2'] != 'N/A':
##            if d_seq[pdb]['REMARK2'] > self.minres:
##                pdbskip = True
##                return pdbskip

        return pdbskip


    def rmsd2bfactor(self, pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_pdb, d_seq, tv1, rm, tv2, l_equivalent_chains, bmchains1, bmchains2):

        biomolecule1 = str(biomolecule1).zfill(2)
        biomolecule2 = str(biomolecule2).zfill(2)

        import os

        pdblines1 = []
        pdblines2 = []
        chains1 = l_equivalent_chains[0]
        chains2 = l_equivalent_chains[1]
        if len(chains1) != len(chains2):
            print pdb1, pdb2
            notexpected
        ##
        ## align residues for the specific combinations of chain1 and chain2
        ##
        d_seqalnparam = self.align_residues(chains1, chains2, d_pdb, d_seq, pdb1, pdb2)

        ##
        ## parse coordinates for the specific combinations of chain1 and chain2
        ##
        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            d_res_nos1 = d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos1']
            d_res_nos2 = d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos2']
            l1 = d_seqalnparam[chain1[0]][chain2[0]]['l1']
            l2 = d_seqalnparam[chain1[0]][chain2[0]]['l2']
            ATOMseq1 = d_seqalnparam[chain1[0]][chain2[0]]['seqATOM1']
            ATOMseq2 = d_seqalnparam[chain1[0]][chain2[0]]['seqATOM2']

            coordinates1, coordinates2, rescount, lines1, lines2 = self.ATOMrecords2coordinates(d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, len(ATOMseq1), rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2)
            pdblines1 += lines1
            pdblines2 += lines2

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


    def tdist(self,t,df):

        ## http://mathworld.wolfram.com/Studentst-Distribution.html (eq. 7)
        import math
        z = float(df)/(float(df)+float(t)**2)
        a = .5*float(df)
        b = .5
        lim1 = 0.
        d_lim2 = {'numerator':z,'denominator':1.}
        inv_stepsize = 1000000.
        function_beta_regularized = {}
        for key in d_lim2:
            lim2 = d_lim2[key]
            deltax = (lim2-lim1)/inv_stepsize
            area = 0.
            x = lim1
            while x <= lim2:
                y = ((x)**(a-1)) * ((1-x)**(b-1))
                area += y*deltax
                x += deltax
            function_beta_regularized[key] = area
        p = (function_beta_regularized['numerator']/function_beta_regularized['denominator'])

        return p        

        
    def write_rmsd_to_file(self, d_rmsd, d_seq, prefix):

        import os

        ## keys=htmlkeys, values=htmlcolumns (table headings)
        d_columns_headers = {
            'REMARK465':'missing residues',
            'REMARK470':'missing atoms',
            'transformations':'transformations',
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
        for column in self.l_columns_html:
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
                        remark465 = max(d_seq[pdb1]['REMARK465'],d_seq[pdb2]['REMARK465'])
                        remark470 = max(d_seq[pdb1]['REMARK470'],d_seq[pdb2]['REMARK470'])
                        transformations = d_rmsd[pdb1][bm1][pdb2][bm2]['transformations']

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
                            'transformations':'<td>%s' %(transformations),
                            'REMARK465':'<td>%s' %(remark465),
                            'REMARK470':'<td>%s' %(remark470),
                            'title1':'<td style="font-size:70%%">%s' %(d_seq[pdb1]['TITLE']),
                            'title2':'<td style="font-size:70%%">%s' %(d_seq[pdb2]['TITLE']),
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

                        d_columns_data_txt = {
                            'transformations':transformations,
                            'REMARK465':remark465,
                            'REMARK470':remark470,
                            'pH1':pH1,
                            'pH2':pH2,
                            'T1':T1,
                            'T2':T2,
                            'res1':res1,
                            'res2':res2,
                            'spacegroup1':spacegroup1,
                            'spacegroup2':spacegroup2,
                            'hetIDs1':','.join(hetIDs1),
                            'hetIDs2':','.join(hetIDs2),
                            'rmsd':rmsd,
                            'chains':n_chains,
                            'residues':n_residues,
                            'coordinates':n_coordinates,
                            'mutations':mutations,
                            }

                        lines_txt = []
                        for key in d_columns_data_txt.keys():
                            value = d_columns_data_txt[key]
                            lines_txt += ['%15s %s\n' %(key, value)]
                        for fd in [
                            'txt/%s/%s.txt' %(pdb1[1],prefix1),
                            'txt/%s/%s.txt' %(pdb2[1],prefix2),
                            ]:
                            fd = open(fd,'w')
                            fd.writelines(lines_txt)
                            fd.close()
                        
                        ## write data to html lines
                        tr = '<tr>\n'
                        for column in self.l_columns_html:
                            tr += '%s</td>\n' %(d_columns_data[column])
                        tr += '</tr>\n'

                        l_tr += tr

                        if not pdb1 in d_html.keys():
                            d_html[pdb1] = ''
                        d_html[pdb1] += tr
                        if not pdb2 in d_html.keys():
                            d_html[pdb2] = ''
                        d_html[pdb2] += tr
                        
        ## write html to global file
        if prefix != 'quickrmsd':
            path = ''
            file = '%s.htm' %(prefix)
            self.append_table_rows(path,file,l_tr,th)

        ## write html to local file
## overwrite if already that pdb comb and bm comb!!!
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

##        print 'identifying different long peptides'

        ## pdb1
        SEQRESchains1_similar_to_SEQRESchains2 = set()
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            ## e.g. 1nvw.pdb,1nvv.pdb
            SEQRESchains1_similar_to_SEQRESchains2 |= set([rep_chain1])
            SEQRESchains1_similar_to_SEQRESchains2 |= set(d_chains_intrapdb_sequence_identical[pdb1][rep_chain1])
        bmSEQRESchains1_similar_to_SEQRESchains2 = SEQRESchains1_similar_to_SEQRESchains2 & set(bmchains1)

        ## pdb2
        SEQRESchains2_similar_to_SEQRESchains1 = set()
        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            rep_chains2 = d_chains_interpdb_sequence_similar[rep_chain1].keys()
            ## e.g. 1nvw.pdb,1nvv.pdb
            SEQRESchains2_similar_to_SEQRESchains1 |= set(rep_chains2)
            for rep_chain2 in rep_chains2:
                SEQRESchains2_similar_to_SEQRESchains1 |= set(d_chains_intrapdb_sequence_identical[pdb2][rep_chain2])
        bmSEQRESchains2_similar_to_SEQRESchains1 = SEQRESchains2_similar_to_SEQRESchains1 & set(bmchains2)

        SEQRESchains1_not_similar_to_SEQRESchains2 = set(d_seq[pdb1]['chains'].keys()) & (set(bmchains1) - SEQRESchains1_similar_to_SEQRESchains2)
        SEQRESchains2_not_similar_to_SEQRESchains1 = set(d_seq[pdb2]['chains'].keys()) & (set(bmchains2) - SEQRESchains2_similar_to_SEQRESchains1)

        return SEQRESchains1_not_similar_to_SEQRESchains2, SEQRESchains2_not_similar_to_SEQRESchains1


    def calculate_rmsd_for_multiple_chains(
        self,chains1,chains2,d_pdb,pdb1,pdb2,d_seq,
        d_chains_interpdb_sequence_similar,
        d_chains_intrapdb_sequence_identical,
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
            if len(chains1) < 60 or len(chains2) < 60:
                print chains1, chains2
            print len(chains1), len(chains2)
            notexpected

        coordinates1 = []
        coordinates2 = []
        residue_count = 0

        ##
        ## align residues for the specific combinations of chain1 and chain2
        ##
        d_seqalnparam = self.align_residues(chains1, chains2, d_pdb, d_seq, pdb1, pdb2)

        ##
        ## parse coordinates for the specific combinations of chain1 and chain2
        ##
        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            d_res_nos1 = d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos1']
            d_res_nos2 = d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos2']
            l1 = d_seqalnparam[chain1[0]][chain2[0]]['l1']
            l2 = d_seqalnparam[chain1[0]][chain2[0]]['l2']
            ATOMseq1 = d_seqalnparam[chain1[0]][chain2[0]]['seqATOM1']
            ATOMseq2 = d_seqalnparam[chain1[0]][chain2[0]]['seqATOM2']

            (coords1, coords2, rescount
             ) = self.ATOMrecords2coordinates(
                 d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2,
                 l1, l2, len(ATOMseq1),
                 )

            ## append coordinates
            coordinates1 += coords1
            coordinates2 += coords2
            residue_count += rescount

        if len(coordinates1) > 500000:
            print 'structural alignment start'
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


    def align_residues(self, chains1, chains2, d_pdb, d_seq, pdb1, pdb2):

        d_seqalnparam = {}

        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            if chain1[0] in d_seqalnparam.keys():
                if chain2[0] in d_seqalnparam[chain1[0]].keys():
                    continue

            d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(d_pdb, d_seq, pdb1, pdb2, chain1, chain2)
            if not chain1[0] in d_seqalnparam.keys():
                d_seqalnparam[chain1[0]] = {}
            if not chain2[0] in d_seqalnparam[chain1[0]].keys():
                d_seqalnparam[chain1[0]][chain2[0]] = {}
            d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos1'] = d_res_nos1
            d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos2'] = d_res_nos2
            d_seqalnparam[chain1[0]][chain2[0]]['l1'] = l1
            d_seqalnparam[chain1[0]][chain2[0]]['l2'] = l2
            d_seqalnparam[chain1[0]][chain2[0]]['seqATOM1'] = ATOMseq1
            d_seqalnparam[chain1[0]][chain2[0]]['seqATOM2'] = ATOMseq2

        return d_seqalnparam


    def identify_interpdb_equivalent_chains_from_structure(
        self, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_pdb, d_seq,
        biomolecule1, biomolecule2,
        d_biomolecules1, d_biomolecules2,
        ):

        import os

##        print 'identifying structurally equivalent long peptides by structure between pdbs'

        bmchains1 = d_biomolecules1[biomolecule1]['chains']
        bmchains2 = d_biomolecules2[biomolecule2]['chains']

        n_chains = n_residues = n_coordinates = 0
        tv1 = 'N/A'
        rm = 'N/A'
        tv2 = 'N/A'
        transformations = False

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
                d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(
                    d_pdb, d_seq, pdb1, pdb2, rep_chain1, rep_chain2,
                    )
                if l1SEQRES != l1 or l2SEQRES != l2:
                    print pdb1, pdb2
                    expected_not_critical
                if ATOMseq1 != d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['s1']:
                    print pdb1, pdb2
                    expected_not_critical
                if ATOMseq2 != d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['s2']:
                    print pdb1, pdb2
                    expected_not_critical
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['d_res_nos1'] = d_res_nos1
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['d_res_nos2'] = d_res_nos1
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1ATOM'] = l1
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2ATOM'] = l2
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['seqATOM1'] = l1
                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['seqATOM2'] = l2
            
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
            ## identify other rep chains1 with seq sim to rep chains2 ## e.g. 1nvv.pdb,1nvw.pdb (Q;R;S == Q,R;S)
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
                            if len(tchain) > 1:
                                if int(tchain[2:]) > 1:
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
            tchains1 += tchains
        tchains2 = []
        for tchains in l_tchains2:
            tchains2 += tchains

        ## check if the expected correct combination of chains gives a low rmsd
        (rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2
         ) = self.calculate_rmsd_for_multiple_chains(
             tchains1,tchains2,d_pdb,pdb1,pdb2,d_seq,
             d_chains_interpdb_sequence_similar,
             d_chains_intrapdb_sequence_identical,
             )
        if rmsd < self.maxrmsd or (len(tchains1) == 1 and len(tchains2) == 1):
            if rmsd > 10.:
                self.incorrecttransformation(pdb1,pdb2,biomolecule1,biomolecule2,d_seq,rmsd,tchains1,tchains2)
            l_equivalent_chains = [tchains1,tchains2]
            return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations

        if len(l_tchains1) != len(l_tchains2):
            print pdb1, pdb2
            print l_tchains1, l_tchains2
            notexpected
        if len(tchains1) != len(tchains2):
            print pdb1, pdb2
            print tchains1, tchains2
            notexpected
        for i in range(len(l_tchains1)):
            if len(l_tchains1[i]) != len(l_tchains2[i]):
                print pdb1, pdb2
                print l_tchains1, l_tchains2
                notexpected

        chains2 = []
        chains1 = []
        ksort = []
        for i in range(len(l_tchains1)):
            tchains1 = l_tchains1[i]
            tchains2 = l_tchains2[i]
            if len(tchains1) != len(tchains2):
                notexpected
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
                    rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_seq,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,verbose=False)
                    ## expected combination of seq sim chains
                    if rmsd < self.maxrmsd:
                        chains2 += tchains2
                else:
                    jskip = []
                    ##
                    ## loop over chains1
                    ##
                    for j in range(1,len(tchains1)+1):
                        if j in jskip:
                            continue
                        if len(tchains1) % 60 == 0 and (j-1) % 60 == 0 and len(tchains1) >= 60:

                            rmsdchains1 = chains1+tchains1[:j-1+60]
##                            print 'rmsdchains1', rmsdchains1

                            ## e.g. chains A,B,C of 1w39.pdb, 2fz1.pdb
                            if ksort != range(1,61) and ksort != []:
                                rmsdchains2 = list(chains2)
                                for k in ksort:
                                    rmsdchains2 += [tchains2[k-1]]

##                                print 'rmsdchains2', rmsdchains2
                                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_seq,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,verbose=False)
                                print rmsd, 'previous sequence'
                                print ksort
                                if rmsd < self.maxrmsd:
                                    chains2 = rmsdchains2
                                    tchains2 = tchains2[60:]
                                    jskip = range(j,j+60)
                                    continue
## maybe A-transformation == C-transformation != B-transformation... then list of ksorts instead of resetting...

                            ## e.g. chain C of 1aq4.pdb, 2bq5.pdb
                            if len(tchains1) > 60:
                                rmsdchains2 = list(chains2)+tchains2[:60]

##                                print 'rmsdchains2', rmsdchains2
                                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_seq,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,verbose=False)
                                print rmsd, 'integer sequence'
                                if rmsd < self.maxrmsd:
                                    chains2 = rmsdchains2
                                    tchains2 = tchains2[60:]
                                    jskip = range(j,j+60)
                                    ksort = []
                                    continue

                            jskip = []
                            ksort = []
                            
                        minrmsd = ['N/A','N/A']
                        ##
                        ## loop over chains2
                        ##
                        for k in range(len(tchains1)-j+1):
                            chain2 = tchains2[k]
                            rmsdchains1 = chains1+tchains1[:j]
                            rmsdchains2 = chains2+[chain2]
                            ## e.g. 2frp.pdb, 2gp1.pdb (chain X1 == chain X1)
                            if len(tchains1) > 60 and j == 1 and (k % 60) != 0:
                                continue
                            ## e.g. 2g33.pdb, 2g34.pdb
                            if len(tchains1) > 60 and (j-1) % 60 != 0 and chains2[-1][0] != chain2[0]:
                                continue
                            rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_seq,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,verbose=False)
                            if len(tchains1) >= 60:
                                print 'rmsd %s %s %4s %4s %5.1f %3s of %3s' %(pdb1, pdb2, tchains1[j-1], chain2, rmsd, j, len(tchains1))
                            ## break if rmsd beneath treshold
                            if rmsd <= self.maxrmsd:
                                chains2 += [chain2]
                                tchains2.remove(chain2)
                                if len(tchains1) >= 60  and len(ksort) < 60:
                                    if len(chain2) > 1:
                                        ksort += [int(chain2[2:])]
                                    else:
                                        ksort += [1]
                                break
                            ## determine if lower rmsd
                            if rmsd < minrmsd[1]:
                                minrmsd = [chain2,rmsd]
                        ## append chain with the lowest rmsd
                        if rmsd > self.maxrmsd:
                            if k != len(tchains1)-j:
                                notexpected
                            chain2 = minrmsd[0]
                            chains2 += [chain2]
                            tchains2.remove(chain2)
                            if len(tchains1) >= 60 and len(ksort) < 60:
                                if len(chain2) > 1:
                                    ksort += [int(chain2[2:])]
                                else:
                                    ksort += [1]
                            if len(tchains1) >= 60:
                                print 'minrmsd %s %s %4s %4s %5.1f %3s of %3s' %(pdb1, pdb2, tchains1[j-1], chain2, minrmsd[1], j, len(tchains1))
                            
            chains1 += tchains1

        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(chains1,chains2,d_pdb,pdb1,pdb2,d_seq,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,verbose=True)
        l_equivalent_chains = [chains1,chains2]

        if rmsd > 10.:
            self.incorrecttransformation(pdb1,pdb2,biomolecule1,biomolecule2,d_seq,rmsd,chains1,chains2)

        return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations


    def incorrecttransformation(self,pdb1,pdb2,biomolecule1,biomolecule2,d_seq,rmsd, chains1, chains2):

##        pdbs = ['1chh','1chj']
##        if pdb1 in pdbs or pdb2 in pdbs:
##            return
        
        ## toomany or just a high rmsd (due to incorrect transformation or something else...)
        toomany = []
        if os.path.isfile('incorrecttransformation.txt'):
            fd = open('incorrecttransformation.txt','r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                toomany += [ [line.split()[0],line.split()[1],line.split()[2],line.split()[3]] ]
        if not [pdb1,pdb2,biomolecule1,biomolecule2] in toomany and not [pdb2,pdb1,biomolecule2,biomolecule1] in toomany:
            fd = open('incorrecttransformation.txt','a')
            fd.write('%s %s %s %s %5.1f %5s %10s %10s %s %s %s\n' %(pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_seq[pdb1]['CRYST1']==d_seq[pdb2]['CRYST1'], d_seq[pdb1]['CRYST1'], d_seq[pdb2]['CRYST1'], len(chains1), chains1[:6], chains2[:6]))
            fd.close()

        return


    def matrixtransformation(self,d_pdb,s_pdb,chain,matrix,matrix_no):

        '''apply REMARK350 transformation'''

        ## matrix does not cause transformation
        if matrix == self.nontransformationmatrix:
            return d_pdb,chain

        import Numeric

##        print 'applying REMARK350 transformation to chain %s of pdb %s' %(chain, s_pdb)

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

        for i in range(len(chains)):
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

        return d_chains_intrapdb_sequence_identical


    def identify_similar_chains_from_sequence_inter(
        self, d_seq, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical,
        bmchains1, bmchains2,
        d_biomolecules1, d_biomolecules2,
        ):

##        print 'identifying similar chains'

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

                ## fast sequence comparison
                if len(seq1) == len(seq2):

                    ## sequence identical
                    if seq1 == seq2:

                        l1 = 0
                        l2 = 0
                        r1 = 0
                        r2 = 0
                        s1 = seq1
                        s2 = seq2
                        n_chainmutations = 0
                        l_chainmutations = []

                    ## point mutation(s)
                    else:

                        l1 = 0
                        l2 = 0
                        r1 = 0
                        r2 = 0
                        s1 = seq1
                        s2 = seq2
                        n_chainmutations, l_chainmutations = self.point_mutations(seq1, seq2)

                else:

                    ## 1st slow sequence comparison (SEQRESseq)

                    print pdb1, pdb2, chain1, chain2, 'begin seq aln of chains of len %s and %s' %(len(seq1),len(seq2))
                    instance = sequence_alignment.NW(seq1,seq2)
##                    print 'aligning chain %s,%s of %s,%s' %(chain1,chain2,pdb1,pdb2)
                    s1,s2 = instance.Align(verbose=False)[:2]
                    print 'end seq aln'

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

                    ## continue if insertions/deletions
                    if '-' in s1 or '-' in s2:
                        continue

                    n_chainmutations, l_chainmutations = self.point_mutations(s1, s2)

                if n_chainmutations <= self.max_mutations:
                    if chain1 not in d_chains_interpdb_sequence_similar.keys():
                        d_chains_interpdb_sequence_similar[chain1] = {}
                    if chain2 in d_chains_interpdb_sequence_similar[chain1].keys():
                        notexpected
                    d_chains_interpdb_sequence_similar[chain1][chain2] = {
                        'l1':l1,'l2':l2,
                        's1':s1,'s2':s2,
                        'l_mutations':l_chainmutations,
                        'r1':0,'r2':0,
                        }

##        print d_chains_interpdb_sequence_similar
    
        return d_chains_interpdb_sequence_similar


    def point_mutations(self, s1, s2):

        n_chainmutations = 0
        l_chainmutations = []
        for res in range(len(s1)):
            res1 = s1[res]
            res2 = s2[res]
            if res1 != res2:
                n_chainmutations += 1
                l_chainmutations += [[res,res1,res2]]
            if n_chainmutations > self.max_mutations:
                break

        return n_chainmutations, l_chainmutations

        


    def res_name2res_symbol(self, res_name):

        if res_name in self.d_res.keys():
            symbol = self.d_res[res_name]
        else:
            symbol = 'X'

        return symbol
    

    def ATOM2seq(self, d_pdb, chain):

## incorrect d_res_nos will be returned from ATOM2seq for 2bfk (vs 2bfl), chain A becauseof ASN61A in REMARK465 records!!!

        d_res_nos = {}
        seq = ''
        ATOMrespos = 0
        res_nos = d_pdb['chains'][chain]['residues'].keys()
        res_nos.sort()
        for i in range(len(res_nos)):
            res_no = res_nos[i]

            for i in range(len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])):
                iCode = d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][i]
                d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':iCode}
                res_name = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
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
                ATOMseqgaplen += 1 ## not a gap but a reversal in the case of 2bfk.pdb!!!
                try:
                    seq += SEQRESseq[SEQRESpos]
                except:
                    print pdb, chain, SEQRESrange1, SEQRESrange2
                    print SEQRESseq
                    print seq
                    print self.pdb1, self.pdb2
                    print SEQRESseq[SEQRESrange1:]
                    print seq[SEQRESrange1:]
                    for i in range(len(seq)):
                        if seq[:i] != SEQRESseq[:i]:
                            print seq[:i]
                            print SEQRESseq[:i]
                            stop1
                    stop2
            ## append if no gap
            else:
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':res_no,'iCode':iCode}
                seq += res_symbol

        ## append C-terminal gap between ATOMseq and SEQRESseq
        if len(seq) != len(SEQRESseq):
            SEQRESrange1 = SEQRESpos+1
            SEQRESrange2 = len(SEQRESseq)
            seq, d_res_nos_SEQRES = self.append_missing_residues_to_sequence(seq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq)

##        if seq != SEQRESseq:
##            print pdb
##            print seq
##            print SEQRESseq
##            notexpected

        return seq, d_res_nos_SEQRES

    def identify_missing_terminal_residues(self, d_pdb, pdb, chain, ATOMseq, SEQRESseq):

##        print pdb, chain

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
            ATOMseq,d_res_nos = self.ATOM2seq(d_pdb[pdb], chain)
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
        if len(ATOMseq1) == len(ATOMseq2):

            l1 = 0
            l2 = 0

        else:

            ## 2nd slow sequence comparison (ATOMseq)

            print pdb1, pdb2, chain1, chain2, 'begin seq aln'
            instance = sequence_alignment.NW(ATOMseq1,ATOMseq2)
            s1,s2 = ATOMs1,ATOMs2 = instance.Align(verbose=False)[:2]
            print 'end seq aln'

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

        return d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2


    def ATOMrecords2coordinates(
        self, d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, len_ATOMseq,
        rmsd=None, tv1=None, rm=None, tv2=None):

        import Numeric, math, copy

        rescount = 0

        coordinates1 = []
        coordinates2 = []
        lines1 = []
        lines2 = []

        for SEQRESpos1 in range(l2,len_ATOMseq+l2):

            SEQRESpos2 = SEQRESpos1+l1-l2
##            print pdb1, pdb2, chain1, chain2, SEQRESpos1, SEQRESpos2
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
                d_atoms = copy.deepcopy(d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'])
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
                if rm:
                    coordinate2 = Numeric.matrixmultiply(rm, coordinate2-tv1)+tv2
                    d_atoms2[atom_name]['coordinate'] = coordinate2
                coordinates1 += [coordinate1]
                coordinates2 += [coordinate2]
                if rmsd:
                    SS += [sum((coordinate2-coordinate1)**2)]

            if rmsd:
                RMSD = math.sqrt(sum(SS)/len(SS))
                d_line = {
                    pdb1:{'res_name':res_name1,'chain':chain1,'res_no':res_no1,'d_atoms':d_atoms1,'iCode':iCode1},
                    pdb2:{'res_name':res_name2,'chain':chain2,'res_no':res_no2,'d_atoms':d_atoms2,'iCode':iCode2},
                    }
                for pdb in d_line:
                    lines = []
                    d_atoms = d_line[pdb]['d_atoms']
                    res_name = d_line[pdb]['res_name']
                    chain = d_line[pdb]['chain']
                    res_no = d_line[pdb]['res_no']
                    iCode = d_line[pdb]['iCode']
                    for atom_name in d_atoms.keys():
                        if 'REMARK' in d_atoms[atom_name].keys():
                            continue
                        coordinate = d_atoms[atom_name]['coordinate']
                        occupancy = bfactor = RMSD/rmsd
                        line = self.coordinates2ATOMline(res_name, chain[0], res_no, coordinate, iCode, bfactor, atom_name)
                        lines += line
                    d_line[pdb]['lines'] = lines
                lines1 += d_line[pdb1]['lines']
                lines2 += d_line[pdb2]['lines']
            if line470 != 'REMARK 470     %s %s%4i    ' %(res_name2, chain, res_no2):
                line470 = "                '"+line470+"', ##%s.pdb %s.pdb\n" %(pdb2, pdb1)
##                fd = open('missingatoms.txt','a')
##                fd.write(line470)
##                fd.close()

        if rmsd:
            return coordinates1, coordinates2, rescount, lines1, lines2
        else:
            return coordinates1, coordinates2, rescount


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


    def parse_coordinates(self, s_pdb, d_pdb, d_seq, d_hetero):

##        print 'parsing coordinates'

        ## read lines
##        fd = open('%s%s.pdb' %(self.pdbpath, s_pdb.upper()),'r')
        fd = open('%s%s/pdb%s.ent' %(self.pdbpath, s_pdb.lower()[1:3], s_pdb.lower()),'r')
        lines = fd.readlines()
        fd.close()
        if s_pdb not in d_pdb.keys():
            d_coordinates,l_connections = self.parse_pdbcoordinatesection(lines, s_pdb, d_seq)
            d_pdb[s_pdb] = d_coordinates
            d_hetero[s_pdb] = l_connections

        return d_pdb, d_hetero


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
                if hetID in ['HOH','H20','DOD','D2O']: ## D2O in 2JAJ
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
                methods = line[10:].strip().split(',')[0]
                if methods[:3] == 'NMR':
                    methods = 'NMR'
                elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                    methods = 'X-RAY'
                d_seq['EXPDTA'] = methods

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                d_seq['CRYST1'] = spacegroup

        ##
        ## missing records
        ##
        d_dics = {
            'chains':{},
            'HET':{},
            'REMARK200':{'TEMPERATURE':'N/A','PH':'N/A'},
            'REMARK525':[],
            'REMARK465':False,
            'REMARK470':False,
            'TITLE':'N/A',
            'EXPDTA':'N/A',
            'REMARK2':'N/A',
            }
        for key in d_dics.keys():
            if key not in d_seq.keys():
                d_seq[key] = d_dics[key]

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

##        if s_pdb not in ['1ady','1bhj']:
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
            d_seq['chains'][chain]['type'] = 'N/A'

        residues = line[19:70].split()

        for i in range(len(residues)):
            residue = residues[i]
            if residue in self.d_res.keys():
                if d_seq['chains'][chain]['type'] == 'N/A':
                    d_seq['chains'][chain]['type'] = 'peptide'
                residues[i] = self.d_res[residue]
            elif residue in ['C','A','U','G','I','DC','DA','DT','DG','DI','N']: ## N is any 5'-monophosphate nucleotide
                if d_seq['chains'][chain]['type'] == 'N/A':
                    d_seq['chains'][chain]['type'] = 'nucleotide'
                residues[i] = residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_seq['chains'][chain]['type'] == 'N/A':
                    d_seq['chains'][chain]['type'] = 'saccharide'
                residues[i] = residue
            else:
                if residue == 'UNK': ## e.g. 1pny.pdb
                    if d_seq['chains'][chain]['type'] == 'N/A':
                        d_seq['chains'][chain]['type'] = 'peptide'
                residues[i] = 'X'

        if 'seq' not in d_seq['chains'][chain].keys():
            d_seq['chains'][chain]['seq'] = ''
        d_seq['chains'][chain]['seq'] += ''.join(residues)

        return d_seq


    def parse_pdbcoordinatesection(self, lines, s_pdb, d_seq):

##        print 'parsing coordinates of %s' %(s_pdb)

        s_pdb = s_pdb.lower()

        ##
        ## set dictionaries
        ##
        d_atomnos = {}
        d_CONECT = {}

        d_coordinates = {
            'chains':{},
            'HET':set(),
            }

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'ATOM':
                d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, s_pdb)
                d_atomnos[d_line['atom_no']] = d_line

            elif record == 'HETATM':
                res_name = line[17:20].strip()
                ## modified residue of polypeptide or polynucleotide
                if res_name in ['HOH','H20','DOD','D2O']:
                    continue
                elif res_name in self.d_modres.keys(): ## e.g. 1gcj.pdb
                    d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, s_pdb)
                ## water
                ## (poly)saccharide or other hetero compound
                else:
                    d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, s_pdb)
                atom_no = d_line['atom_no']
                d_atomnos[atom_no] = d_line

            elif record == 'CONECT':
                d_CONECT = self.parse_recordCONECT(line, d_atomnos, d_CONECT)
                
            elif record == 'REMARK':
                remark = int(line[6:10])
                if remark == 465:
                    d_coordinates = self.parse_recordREMARK465(line, d_coordinates, lines, i)
                elif remark == 470:
                    d_coordinates = self.parse_recordREMARK470(line, d_coordinates, lines, i)

            elif record == 'HET':
                hetID = line[7:10].strip()
                d_coordinates['HET'] |= set([hetID])
##                chain = line[12]
##                res_no = int(line[13:17])
##                iCode = line[17]
##                if iCode != ' ':
##                    print s_pdb
##                    stop
##                n_atoms = int(line[20:25])
##                description = line[30:70]
##                if chain not in d_hetero['chains'].keys():
##                    d_hetero['chains'][chain] = {}
##                if 'residues' not in d_hetero['chains'][chain].keys():
##                    d_hetero['chains'][chain]['residues'] = {}
##                if res_no not in d_hetero['chains'][chain]['residues'].keys():
##                    d_hetero['chains'][chain]['residues'][res_no] = {}
##                if 'iCodes' not in d_hetero['chains'][chain]['residues'][res_no].keys():
##                    d_hetero['chains'][chain]['residues'][res_no]['iCodes'] = {}
##                d_hetero['chains'][chain]['residues'][res_no]['iCodes'][iCode] = hetID
                
            elif record == 'MODEL':
                model = int(line.split()[1])


        d_adjacency_backward = {} ## C --> O
        d_adjacency_forward = {} ## O --> C
        matrix_adjacency = []
        l_connections = []
        set_monomers = set()
        set_ATOMatomnos = set()
        d_connections = {}
        l_monomers = []
        for hetID1 in d_CONECT.keys():
            for chain1 in d_CONECT[hetID1].keys():
                for res_no1 in d_CONECT[hetID1][chain1].keys():
                    for iCode1 in d_CONECT[hetID1][chain1][res_no1].keys():
                        for atom_no1 in d_CONECT[hetID1][chain1][res_no1][iCode1].keys():
                            atom_name1 = d_atomnos[atom_no1]['atom_name']
                            for atom_no2 in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]:
                                hetID2 = d_atomnos[atom_no2]['res_name']
                                chain2 = d_atomnos[atom_no2]['chain']
                                res_no2 = d_atomnos[atom_no2]['res_no']
                                iCode2 = d_atomnos[atom_no2]['iCode']
                                atom_name2 = d_atomnos[atom_no2]['atom_name']
                                ## glycoside bond 1
                                if atom_name1[0] == 'C':
                                    monomer1 = chain1+str(res_no1)+iCode1
                                    monomer2 = chain2+str(res_no2)+iCode2
                                    if  atom_name2[0] == 'O':
                                        bond = atom_name1[1:]+','+atom_name2[1:]
                                    elif atom_name2 == 'ND2':
                                        bond = atom_name1[1:]+','+'N'
                                        if hetID2 != 'ASN':
                                            stop
                                    else:
                                        stop
                                ## glycoside bond 2
                                elif atom_name2[0] == 'C':
                                    monomer1 = chain2+str(res_no2)+iCode2
                                    monomer2 = chain1+str(res_no1)+iCode1
                                    if  atom_name1[0] == 'O':
                                        bond = atom_name2[1:]+','+atom_name1[1:]
                                    elif atom_name1 == 'ND2':
                                        bond = atom_name2[1:]+','+'N'
                                        if hetID1 != 'ASN':
                                            stop
                                    else:
                                        stop
                                ## disulphide bond
                                elif atom_name1 == 'SG' and atom_name2 == 'SG':
                                    bond = 'SS'
                                    monomer1 = chain1+str(res_no1)+iCode1
                                    monomer2 = chain2+str(res_no2)+iCode2
                                else:
                                    print atom_name1, atom_name2
                                    stop

                                connection = '"'+monomer1+'"-'+bond+'-"'+monomer2+'"'

                                if atom_name1 != 'SG' and atom_name2 != 'SG':
                                    if hetID1 == 'ASN':
                                        set_ATOMatomnos |= set([atom_no1])
                                    set_monomers |= set([monomer1,monomer2])
                                    d_connections[monomer1] = monomer2
                                    
                                    if not monomer2 in d_adjacency_forward.keys():
                                        d_adjacency_forward[monomer2] = set([monomer1])
                                    else:
                                        d_adjacency_forward[monomer2] |= set([monomer1])
                                        
                                    if monomer1 in d_adjacency_backward.keys() and d_adjacency_backward[monomer1] != monomer2:
                                        print monomer1, d_adjacency_backward
                                        notexpected
                                    d_adjacency_backward[monomer1] = monomer2
                                    
                                print hetID1,chain1,res_no1,iCode1,atom_name1,hetID2,chain2,res_no2,iCode2,atom_name2,connection


        ## identify roots of saccharides (simplified back trace)
        set_roots = set_monomers-set(d_adjacency_backward.keys())

        ## trace saccharide tree forward and identify the saccharide molecule
        d_molecules = {}
        l_branches = []
        for root in set_roots:
            path = [root]
            d_molecules[root] = {}
            d_m_fwd = d_molecules[root]
            l_monomers = [root]
            monomer1 = root
            while True:
                print 'xxx', monomer1, d_molecules[root]

                ## end of branch and no branching points                    
                if monomer1 not in d_adjacency_forward.keys() and len(l_branches) == 0:
                    break

                ## end of branch but branching points
                elif monomer1 not in d_adjacency_forward.keys() and len(l_branches) > 0:

                    ## rewind path to previous branching point
                    path = path[:path.index(l_branches[-1])+1]
##                    print 'a0path', path, monomer1, l_branches
                    ## forward in dictionary to previous branching point
                    d_m_fwd = d_molecules
                    for monomer in path:
                        d_m_fwd = d_m_fwd[monomer]

                    l_monomers.append(monomer1)

                    monomers2 = list(set(d_m_fwd.keys())-set(l_monomers))

                    ## end of branch
                    if len(monomers2) == 0:
                        monomer1 = path[-2]
                        l_branches.remove(path[-1])
                        ## rewind path
                        path = path[:-1]
##                        print 'a1path', path, monomer1, l_branches
                        ## forward dictionary
                        d_m_fwd = d_molecules
                        for monomer in path:
                            d_m_fwd = d_m_fwd[monomer]

                    else:
                        monomer1 = monomers2[0]
                        ## forward path
                        path.append(monomer1)
##                        print 'a2path', path, monomer1, l_branches
                        ## forward dictionary
                        d_m_fwd = d_m_fwd[monomer1]


                else:

                    l_monomers.append(monomer1)

                    monomers2 = list(d_adjacency_forward[monomer1]-set(l_monomers))

                    if len(monomers2) == 0 and len(l_branches) == 0:
                        break
                    if len(monomers2) == 0 and len(l_branches) == 1 and monomer1 == l_branches[0]:
                        break
                    if len(monomers2) == 0:
                        print l_branches
                        print d_adjacency_forward[monomer1]
                        notexpected

                    for monomer2 in monomers2:
                        d_m_fwd[monomer2] = {}

                    if len(monomers2) > 1:
                        l_branches += [monomer1]

                    monomer1 = monomers2[0]
                    
                    ## forward path
                    path.append(monomer1)
##                    print 'b0path', path, monomer1, l_branches
                    ## forward dictionary
                    d_m_fwd = d_m_fwd[monomer1]


        print d_molecules

        print set_monomers
        print d_adjacency
        l_monomers_connected = []
        for monomer1 in d_adjacency.keys():
            monomers1a = d_adjacency[monomer2a]
            for monomer2b in d_adjacency.keys():
                if monomer2a == monomer2b:
                    continue
                monomers1b = d_adjacency[monomer2b]
                
                print monomer2a,monomers1a,monomer2b,monomers1b
                stop
            monomers3 = set([])
            for monomer1 in monomers1:
                if monomer1 in d_adjacency.keys():
                    monomers = d_adjacency[monomer1]
                    
                           

        set_connections = set([])
        for j in range(len(l_connections)):
            set_connections |= l_connections[j]
        print l_connections

        for j in range(len(l_connections)):
            connections = list(l_connections[j])
            l_hetIDs = []
            for conect in connections:
                chain = conect[0]
                res_no = int(conect[1:-1])
                iCode = conect[-1]
                print conect
                print l_connections[j]
                res_name = d_atomnos['res_name']
                if res_name in self.d_saccharides.keys():
                    res_name = self.d_saccharides[res_name]['stereo']
                l_hetIDs += [res_name]
            set_hetIDs = set(l_hetIDs)
            d_hetIDs = {}
            for hetID in set_hetIDs:
                count = l_hetIDs.count(hetID)
                d_hetIDs[hetID] = count
            l_connections[j] = d_hetIDs

        l_connections.sort()
        print l_connections
        stop

        return d_coordinates, l_connections


    def monomertranslation(self,monomer,d_coordinates,d_molecules,path):

        chain = monomer[0]
        res_no = int(monomer[1:-1])
        iCode = monomer[-1]
        monomer_translated = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']

        d_m_fwd = d_molecules
        for m in path[:-1]:
            d_m_fwd = d_m_fwd[m]
        d_m_fwd[monomer_translated] = d_m_fwd[monomer]
        path = path[:-1]+[monomer_translated]
        del d_m_fwd[monomer]
        d_m_fwd = d_m_fwd[monomer_translated]

        return d_m_fwd,path


    def parse_recordCONECT(self,line,d_atomnos,d_CONECT):

        print line[:-1]

        ## parse atom nos
        atom_nos = []
        for j in range(6,31,5):
            if line[j:j+5] == '     ':
                break
            atom_nos += [int(line[j:j+5])]

        ## parse res nos of *hetero* atoms
        atom_no1 = atom_nos[0]

        ## continue if water
        if atom_no1 not in d_atomnos.keys():
            return d_CONECT
        
        chain1 = d_atomnos[atom_no1]['chain']
        res_no1 = d_atomnos[atom_no1]['res_no']
        iCode1 = d_atomnos[atom_no1]['iCode']
        res_name1 = d_atomnos[atom_no1]['res_name']
        atom_name1 = d_atomnos[atom_no1]['atom_name']
        if atom_name1[0] == 'H':
            print line
            stop
        ## append res/seq no if HETATM is not a metal ion to which molecules are not *covalently* bound
        ## and if HETATM is not a modified peptide residue
        if res_name1 in self.d_ions.keys()+self.l_clusters:
            return d_CONECT

        l_atom_nos = []
        for atom_no2 in atom_nos[1:]:
            chain2 = d_atomnos[atom_no2]['chain']
            res_no2 = d_atomnos[atom_no2]['res_no']
            iCode2 = d_atomnos[atom_no2]['iCode']
            res_name2 = d_atomnos[atom_no2]['res_name']
            if res_name2 in self.d_ions.keys()+self.l_clusters:
                continue
            if (
                chain1 !=  chain2 or
                res_no1 != res_no2 or
                iCode1 != iCode2
                ):
                l_atom_nos += [atom_no2]

        if len(l_atom_nos) == 0:
            return d_CONECT

        if res_name1 not in d_CONECT.keys():
            d_CONECT[res_name1] = {}
        if chain1 not in d_CONECT[res_name1].keys():
            d_CONECT[res_name1][chain1] = {}
        if res_no1 not in d_CONECT[res_name1][chain1].keys():
            d_CONECT[res_name1][chain1][res_no1] = {}
        if iCode1 not in d_CONECT[res_name1][chain1][res_no1].keys():
            d_CONECT[res_name1][chain1][res_no1][iCode1] = {}
        
        if res_name1 not in ['NAG','ASN','CYS','BMA','MAN',]:
            print atom_no1, atom_no2, l_atom_nos
            print d_atomnos[atom_no1]
            print d_atomnos[atom_no2]
            stop

        d_CONECT[res_name1][chain1][res_no1][iCode1][atom_no1] = l_atom_nos

        return d_CONECT


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

        elif remark == 465: ## missing residues

            d_seq['REMARK465'] = True

        elif remark == 470: ## missing atoms

            d_seq['REMARK470'] = True

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


    def parse_recordATOM(self, line, d_pdb, lines, i, d_seq, s_pdb):

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
        ## temp!!!
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name: ## temp!!!
            print res_name, d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] ## temp!!!
            print chain, res_no, iCode ## temp!!!
            stop
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
##########
##        l_iCodes = d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']
##        l_iCodes_REMARK465 = []
##        for iCode_REMARK465 in l_iCodes:
##            if 'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_REMARK465].keys():
##                l_iCodes_REMARK465 += [iCode_REMARK465]
##                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'].remove(iCode_REMARK465)
##                ATOMseq,d_res_nos = self.ATOM2seq(d_pdb, chain)
##                print ATOMseq
##                print d_seq['chains'][chain]['seq']
##                print s_pdb, chain, res_no
##                if s_pdb == '2bfk':
##                    stop
####        if not iCode in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
####            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
####        d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += l_iCodes_REMARK465
##########
        if not iCode in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

        ## iCode > atoms
        if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
        ## atoms > atom_name > coordinate
        if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}

        return d_pdb, {'chain':chain,'atom_no':atom_no,'atom_name':atom_name,'res_name':res_name,'res_no':res_no,'iCode':iCode}


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

        ## info from the PDB Ligand Depot (searched for cluster in chemical name)
        self.l_clusters = [
            ## iron clusters
            'CFM','CFN','CLF','CLP','CN1','CNB','CNF','F3S','FES','FS1','FS2','FS3','FS4','FSF','FSO','FSX','HC0','HC1','HF3','HF5','NFS','SF3','SF4','WCC','XCC',
            ## copper clusters
            'CUB','CUM','CUN','CUO',
            ## molybdenum "clusters"
            'OMO',
            ## hafnium clusters
            'PHF',
            ## zirconium clusters
            'ZRC',
            ]

        ## info from the PDB Ligand Depot
        ## keys are hetIDs, values are charges (not oxidation states)
        ## hetID:[chemical formula,charge,number of associated water molecules]
        self.d_ions = {
            ## acetate
            'ACT':['C2 H3 O2',-1],
            ## nitrate, ammonium
            'NO3':['N1 O3',-1],'NH4':['H4 N1',+1],
            ## hydroxide
            'OH' :['H1 O1',-1],
            ## phosphate
            '2HP':['O4 P1',-1],'IPS':[3,-2],'PI' :['H1 O4 P1',-2],'PO4':['O4 P1',-3], ## synonyms and different oxidation states
            ## sulfate
            'SO4':['O4 S1',-2],'SOH':[3,-1],'SUL':[3,-2], ## synonyms and different oxidation states
            ## group1a
            'LI' :['LI1',+1],
            'NA' :['NA1',+1],'NAO':['NA1',+1],'NA2':['NA1',+1],'NAW':['NA1',+1],'NA5':['NA1',+1],'NA6':['NA1',+1], ## different number of waters coordinated
            'K'  :['K1' ,+1],'KO4':['K1' ,+1], ## different number of waters coordinated
            ## group2a
            'MG' :['MG1',+2],'MO1':[1,+2],'MO2':[1,+2],'MO3':[1,+2],'MO4':[1,+2],'MO5':[1,+2],'MO6':[1,+2], ## different number of waters coordinated
            'CA' :['CA1',+2],'OC1':['CA1',+2],'OC2':+2,'OC3':+2,'OC4':+2,'OC5':+2,'OC6':+2,'OC7':+2,'OC8':+2, ## different number of waters coordinated
            ## group3a
            'AL' :['AL1',+3],
            'GA' :['GA1',+3,0],'GA1' :['GA1',+2,0],
            ## group4a
            'ARS':['AS1', 0,0],'ART':['O4 AS1',-3,0],'AST':-3,'TAS':['H3 O3 AS1', 0], ## different compounds
            ## group6a
            'SE' :['SE1', 0,0],'SE4':['O4 SE1',-2,0], ## different compounds
            ## group7a
            'CL' :['CL1',-1],
            'BR' :['BR1',-1],
            'IOD':['I1' ,-1],
            ## group8a
            'KR' :['KR1', 0],
            ## group3b
            'V'  :+3,'VO4':['V1' ,-3], ## different oxidation states
            ## group4b
            'CR' :['CR1',+3],
            'MO' :['MO1', 0],
            ## group5b
            'MN' :['MN1',+2],'MN3':['MN1',+3],'MW1':+2,'MW2':+2,'MW3':+2,'O4M':+2,'MN5':+2,'MN6':+2, ## different oxidation states and different number of waters coordinated
            ## group6b
            'FE2':['FE1',+2],'OF1':+2,'2OF':+2,'FE' :['FE1',+3],'OF3':+3, ## different oxidation states and different number of waters coordinated
            ## group7b
            'CO' :['CO1',+2],'OCL':['CO1',+2,1],'OCN':['CO1',+2,2],'OCM':['CO1',+2,3],'3CO':['CO1',+3],'CO5':+3,'OCO':+3, ## different oxidation states and different number of waters coordinated
            ## group8b
            'NI' :['NI1',+2],'3NI':['NI1',+3],'NI1':+2,'NI2':+2,'NI3':+2,'NIK':+2, ## different oxidation states and different number of waters coordinated
            ## group9b
            'CU1':['CU1',+1],'CU' :['CU1',+2],'1CU':+2, ## different oxidation states and different number of waters coordinated
            ## group10b
            'ZN' :['ZN1',+2],'ZN2':+2,'ZN3':+2,'ZNO':+2,'ZO3':+2, ## different number of waters coordinated and different crystal fold axes
            'CD' :['CD1',+2],
            'HG' :['HG1',+2],
            ## Lanthanides
            'TB' :['TB1',+3],
            }

        ## saccharides returned from a search of the ligand depot for disaccharides and monosaccharides
        self.d_saccharides = {
            ##
            ## monosaccharides, aldehydes
            ##
            'GLC':{'stereo':'GLC','derivate':['GLC']}, ## (alpha)-D-Glucose
            'AGC':{'stereo':'GLC','derivate':['GLC']}, ## alpha-D-Glc
            'BGC':{'stereo':'GLC','derivate':['GLC']}, ## beta-D-Glc
            'G1P':{'stereo':'G1P','derivate':['GLC']}, ## alpha-D-Glc-1P
            'G6P':{'stereo':'G6P','derivate':['GLC']}, ## alpha-D-Glc-6P
            'BG6':{'stereo':'G6P','derivate':['GLC']}, ## beta-D-Glc-6P
            'G6Q':{'stereo':'G6P','derivate':['GLC']}, ## Glc-6P
            'GAL':{'stereo':'GAL','derivate':['GAL']}, ## (beta)-D-Galactose
            'GLA':{'stereo':'GAL','derivate':['GAL']}, ## alpha-D-Gal
            'GLB':{'stereo':'GAL','derivate':['GAL']}, ## beta-D-Gal
            'FUC':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, alpha-L-Fucose
            'FUL':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, beta-L-Fucose
            'BGP':{'stereo':'BGP','derivate':['GAL']}, ## beta-Gal-6P
            'MAN':{'stereo':'MAN','derivate':['MAN']}, ## alpha-D-Mannose
            'BMA':{'stereo':'MAN','derivate':['MAN']}, ## beta-D-Mannose
            'M1P':{'stereo':'M1P','derivate':['MAN']}, ## alpha-D-Man-1P
            'M6P':{'stereo':'M6P','derivate':['MAN']}, ## alpha-D-Man-6P
            ## acetylated hexose amines
            'NAG':{'stereo':'NAG','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine
            'NBG':{'stereo':'NAG','derivate':['GLC']}, ## 1-N-Acetyl-beta-D-Glucosamine
            'NDG':{'stereo':'NAG','derivate':['GLC']}, ## 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
            '16G':{'stereo':'16G','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine-6P
            'NGA':{'stereo':'NGA','derivate':['GLC']}, ## (2)-N-Acetyl-D-Galactosamine
            ##
            ## monosaccharides, ketones
            ##
            'FRU':{'stereo':'FRU','derivate':['FRU']}, ## Fructose
            'F6P':{'stereo':'F6P','derivate':['FRU']}, ## Fru-6P
            ##
            ## dissacharides
            ##
            'SUC':{'stereo':'SUC','derivate':['GLC','FRU']}, ## GLC-a12-FRC, Sucrose
            'LAT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, alpha-Lactose
            'LBT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, beta-Lactose
            'MAL':{'stereo':'MAL','derivate':['GLC','GLC']}, ## GLC-a14-GLC, Maltose
            'TRE':{'stereo':'TRE','derivate':['GLC','GLC']}, ## GLC-a11a-GLC, Trehalose
            'CBI':{'stereo':'CBI','derivate':['GLC','GLC']}, ## GLC-b14-GLC, Cellobiose
            }

        self.d_stereoisomers = {
            }

        self.l_expdta = [
            'X-RAY',
            'ELECTRON DIFFRACTION','NEUTRON DIFFRACTION',
            'NMR',
            'INFRARED SPECTROSCOPY',
            'CRYO-ELECTRON MICROSCOPY',
            'ELECTRON TOMOGRAPHY', ## e.g. 1o1a.pdb
            'SOLUTION SCATTERING', ## e.g. 1e07.pdb
            'FLUORESCENCE TRANSFER', ## e.g. 1rmn.pdb
            ]

        self.l_columns_html = [
            'gif1','gif2','pdb1', 'pdb2', 'bm1', 'bm2',
            'rmsd', 'mutations', 'chains', 'residues', 'coordinates',
            'pH1', 'pH2', 'T1', 'T2', 'res1', 'res2','spacegroup1', 'spacegroup2',
            'REMARK465', 'REMARK470','transformations',
            'title1','title2','hetIDs1', 'hetIDs2'
            ]

        self.maxrmsd = 2.75 ## 2fsy.pdb,2ft1.pdb

        self.minres = 5.0

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25
        self.max_mutations = 10

        self.pdbpath = '/oxygenase_local/data/pdb/'
##        self.pdbpath = 'pdbin/' ## temp!!! isomerase cluster jobs...
        
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

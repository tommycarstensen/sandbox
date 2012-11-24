#!/software/bin/python
#
#$Id: quakes.py 266 2007-11-01 13:15:43Z tc $
#
#Tommy Carstensen, University College Dublin, 2007

##
## questions of interest

## find correlation (if any) between RMSD of backbone/all atoms for all/neighboring(exponentional sphere radii)/surface residues

## rmsd not mathematically dependent on chain length! but rmsd physically dep on chain length? plot!

## rmsd plot of structures with "bad" and "good" geometry! quality of structure (resolution, r, r_free-r, bond angle/length deviations)

##
## questions of concern

## maltotetraose (GLC) in 2fhf; maltotriose (MLR=3xGLC) in 2fhc; 2fhb, 2fh8, 2fh6
## parse some sugar molecules parsed and leave the rest if hetID of monomers and not multimer specified
## look up sugar dictionary...

## what to do with 1k56, 1k57?

## 1eia, 2eia does not align but not added to transform_error.txt ????

## ignore 1c2b because PISA can't handle it! put it in the code if other similar structures.

## accept neutron diffraction structures? accept 3ins which is x-ray AND neutron?

## solve the problem of too many chain combinations by sequential pairing
## this will yield a maximum of n(n+1)/2 combinations to check and with fewer coordinates per rmsd calculation
## sequential addition of sequence similar chains to reduce number of combinations further

## plot of rmsd vs disulphide bonds (y/n)
## plot of rmsd vs EC class
## plot of rmsd vs CATH class

## whatif accessibility should be calculated for the biou and not the asyu..!!!

## Sep26 - number of residues should be split in residues1 and residues2 and should include missing residues and mutated residues

## Oct12 - should NAD (NAD+) and NAI (NADH) be considered the same compound, since hydrogen atoms cant be seen in x-ray structures?
## Oct12 - should NAP (NADP+) and NDP (NADPH) be considered the same compound, since hydrogen atoms cant be seen in x-ray structures?
## Oct12 - should NAD (NAD+) and NAJ (NAD+ acidic form) be considered the same compound, since hydrogen atoms cant be seen in x-ray structures?

## incorrect d_res_nos will be returned from ATOM2seq for 2bfk (vs 2bfl), chain A becauseof ASN61A in REMARK465 records!!!

## search and replace self.parse_record with parse_pdb.parse_record

## todo - avoid repeating a cluster - merge all clusters into unique sets - avoid comparing structures with small chain identical chains

## go over incorrecttransformation.txt!!!

## 2009Sep16 - go over errorpdbs in function parse_html

## 2009Oct05 - exclude files if nonterminal residues are unobserved!!! solves 

## 2009Oct06 - find molecules with less than 10 atoms in them and categorize them in smallmolecules

## 2009Oct06 - if *only* software *only* determined biounits, then all biomolecules should be deleted and a flag be set to send directly to pisa upon structural comparison

## 2009Oct07 - BOG causes large conformational change between 1l6l and 2ou1


## built-ins
import os, sys, numpy, math, copy, time, re, urllib2
## non-built-ins
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import statistics,gnuplot
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import biounit,parse_pdb
sys.path.append('/home/people/tc/svn/Protool/')
import geometry
sys.path.append('/home/people/tc/svn/PEAT_DB/')
import sequence_alignment

class quakes:

    def main(self):

        import os

        if os.getcwd() != self.topdir:
            print os.getcwd()
            print self.topdir
            stop_wrong_dir

        for s_dir in ['htm','ps','txt','pdb','tmp',]:
            if not os.path.isdir(s_dir):
                os.mkdir(s_dir)
            if s_dir in ['txt','pdb',]:
                for s in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
                    if not os.path.isdir('%s/%s' %(s_dir,s.lower(),)):
                        os.mkdir('%s/%s' %(s_dir,s.lower(),))

##        self.obsolete()
##        self.rsync()
##        self.gunzip()

################################################################################

        errorpdbs = [
            ## residue renumbering (altloc, not diff res)
            '9abp','1bap','7abp','1abp','1apb','6abp','8abp',
            '2nn8','4xis','1q6f','1e2f',
            '1axa','1gnn','1gno','1gnm','2gsp',
##            ## residue renumbering (N-terminal residues)
##            '1hto','1htq','1fkn','2hk2','2hk3','2qri','2qrs','1t41','1yr9',
##            '1yrb','1pbh','2pbh','1fng','1fne','2z5t',
            '1ef0',
##            ## residue renumbering (internal insertion)
##            '2ht8',##'1jrt','1jrs','1g2w','1mbq','2tbs',
##            ## residue renumbering (reverse)
##            '1zzn','1u6b',

            ## REMARK465 header missing/duplicate
            '1unp','2cjc','1uwg',
            ## REMARK465 records missing
            '2vs6','1h5w','2wfx','2wg4','1w6r','1w1b','1gun','2uv4',
            '2uv7','1uc9','1w9a','1h7n','1w49','2bjg','2wb1','1gmv',
##            '1deq','1uf2','2v0z','2v7n','2vs4','1zbb','1zlv',
            ## REMARK465 iCode change
            '1p63',
            ## REMARK465 residue renumbering (REMARK465)
            '1unr','1e7r','2p8m','3d0l','2gk0','2c24',
            ## REMARK465 residue renaming (REMARK465)
            '1uzm','1dzy','1e4b','1usm',
            ## REMARK465 residue deletion (REMARK465)
            '1suc','1ut0','1v00','1o6v',

##            ## MODRES missing
##            '1cwb', ## DMT
            '2bhg','1h6g', ## MSE (being remediated, MET)
            '1ha7', ## MEN (being remediated, ASN)
            '1qng', ## ABA,SAR (being remediated)
            '2bws', ## CSO (being remediated)

            ## SEQADV missing
            '2bqk','2p9k','1abs','2i3v','2i3w','2mgk','2mgl','2mgm','2mbw',
            '1tye',

            ## SEQADV mislabel (CONFLICT --> MUTATION)
            '1uxw','1oc9','1wt3','1wsz','1wt0','1wt1','1wt2','1pll','2q21',
##            '1jaw','3ink','1ha5','1p09','7lpr','1p10','2lhm','3lhm','1laa',
##            '1tay','1tby','1tcy','1tdy','1kip','1kiq','1kir','1xz6',
##            '1lz4','1lzd','1lze','1her','1hem','1heo','1mpb','2trm',
##            '1lhl','1lhj','1lhh','1lhi','1lhk',
##            '1cgx','1cgw','1cgv','1cgy','1gil',
##            '1lye','1lyf','1lyg','1lyh','1lyi','1lyj',
##            '1fdd','1frk','1frj','1fri','1frh','1frm','1frl','1frx','1fd2',
##            

##            '1l12',
##            '1l19',
##            '1l21','1l22',
##            '1l31',
##            '1l39',
##            '1l40',
##            '1l46',
##            '1l47','1l48',
##            '1l55','1l57','1l59','1l61','1l62','1l63','1l54',
##            '1l82',
##            '1l85','1l86','1l87','1l88','1l89','1l90','1l91','1l92','1l93',
##            '1l94','1l95',
##            '1lzg',

##            '107l','108l','109l','110l','111l','112l','130l','113l','114l','115l',
##            '129l','131l','137l',
##            '149l','150l',
##            '172l',
##            '195l','196l','197l','198l','199l','200l',
##            '216l','217l',

            ## SEQADV VARIANT
            '1thu',

            ## SEQADV CONFLICT (UniProtKB difference, not mutant, unexplained)
            '1b3e','1f1f','1xta','1gam','1aij','1aig','1mkd','1h04','1h2p',
            '4dfr','1bi4','1gk1','1gk0','3fjp','1ha5',

            ## duplicate SEQADV records
            '3ck7','1gag','2om9','2rca','3bfv','3b2z','1yl5','2b9b','3euk',
            '1hjv','1hjx','2ofw','1q6y','1vr1','1xm3','2ig7','1rxj','1v58',
            '1ktk','1b6p','2c9t','1xo5','1sg2','1uc9','1txu','3b94','3gts',
            '2ff2','2aqz','1nzk','2hwa','1w5q','2hwm','1z4s','2qm8','1z2v',
            '2h1h','1jt5','2ntd','2hw9','2hz9','2o2n','1yto','1aw8','1aw8',
            '1r1z','3etm','3evo','3evm','3ena',

            ## SEQADV CONFLICT (ASN/ASP)
            '2foa',

            ## CAVEAT no continuation
            '3dyu',
##            ## REMARK350 errors
##            '1a8r','1a9c','1qom','2scu','1onr','2pab','2v5l','1jfa',
            '1ow6','2c7d','2c7c','1ht2','1ut1','3d7e','1jaw','1r3m',
##            ## pisa entry not found
##            '3bcm',
            ## residue present in REMARK465 *and* ATOM records
            '1ols',

##            ## incorrect atom names
##            '2ric','2exj','2exk',
####             ## vs 1d7e
####            '1x83', ## vs '1x84' ## BR atom missing in hetID SBH
##            ## incorrect residue name?
##            '1agy',
##
####            ## N-terminal posttranslational modification not in SEQRES records
####            '1vba','1al2','1ar6','1ar9','1asj',
####            ## altloc flag missing
##            '1xpk','3xis','1ofz','1abe','1abf','5abp','1axz',
##
##            ## SEQADV incorrect seqID (should be remediated quickly...)
##            '1jt4','1jt7','2yfp','1jtc','1q04','1pzz','1m16','1z93',
##            '1z97','2hw9','2hwa','2hwm','2ntd','2hz9','1nzk',
##
            ##
            ## tommy errors
            ##
            '2gp1', ## virus yields too big an rmsd...
            '1b6q','1gmg', ## yield too big an rmsd because of pro mutation in hinge
            '1ft8','1koh', ## too many missing coordinates (SEQRES overlap but no ATOM overlap)
##            ## alternate locations
##            '1uwn','1uwp', (rmsd = 0)
            ## HETATM moved from beginning of file to end of file
            '1dym', ## PCA
            '2bg7', ## CSO
            '2bwv','1w2m', ## CSO
            ## unexpected SEQADV comment (cant determine wt and mutant)
            '1e5v','1e61','1h5n', ## DMSO - SEQUENCED FROM MAP
            ## should this be split in chain A and B??? compare 1k34.
            '1k33',
            ## excluded because GroEl is rotated by a 22' angle due to a E461K mutation
            '2eu1',
            ## RMSD problem
            '1zbl','1zbi',

            ## HELIX and SHEET errors
            '1qai','6yas','1uoj','1oko','2vxj','1e7n','1e9g','1e0z','1qj4',

            ##
            ## connections
            ##

            ## incorrect connection - oxygen valence exceeded
            '1sf7','1xrl','2a74','2nt0','2po6','1fy1','2fs8','2znh','2dpe',
            '1ljy','1o8a','1h86','2jef','1fv3','2c4a','2c4l','1jne','1d7c',
            '1rvz','1fzj','1gwm','1vbo','1ib4','1etb','1eta','1pty','1y2x',
            '1sk0','1sfq','1h81','1e56',
            ## incorrect connection - sulfur valence exceeded
            '2uul','1aym',
            ## incorrect connection - nitrogen valence exceeded
            '1i8b',
            ## incorrect connection - carbon valence exceeded
            '1zlw','1e4i','1gh7','2wht',
            ## incorrect connection - altlocs
            '1rr8','1one','1zcs','1e0p',
            ## incorrect connection - distance
            '1e3a',
            ## incorrect connection - angle
            '2whu',
            ## incorrect connection - nterminal
            '1spd','2avq','1t1n',
            ## incorrect connection - cterminal
            '1lek','1k2c',
            ## incorrect connection - l-peptide linker
            '1daz',
            ## incorrect connection - other
            '1i89','2gbl','1ms8','2by3','2gaa','1h8v','1k2b','2wgm','1ueu',
            '2nsx','2jl0','1t9u',
            '1uiw', ## LYS:NZ-2FU:O8
##            '2c03','1ljy','1fzg','1n2c','2bod','1p8h','2fa7','2cfg',
##            '2bhy','1sbd','1sbe','1ppv','1x83','1ea5','2c3w',
##            '1ayn','2hrr',
            '1r9m',
            '1fmu', ## NDG:O4-NAG:O1
            '2bwc', ## GLC:O4-BGC:O1 (remediation in progress...)
            '1vfk', ## GLC:O1-GLC:O6 (remediation in progress...)

            ## missing connections
            '1eb4',

            ## unexpected ARG connection
            '4atj', ## ARG:NH2-BHO:O1
            '2sns', ## ARG:NH2-THP:O6P
            '1vha', ## ARG:NH2-POP:O5
            '1ake', ## ARG:NH1-AP5:O1D
            '2hg5', ## ARG:NH1-B3H:O11 (possibly incorrect, introduced during remediation)
            '2pyp', ## ARG:NH2-HC4:C3'

            ## unexpected oxygen-oxygen connection
            '1k7d', ## SER:OG-GRO:O2
            '1epv', ## SER:OG-DCS:O3P
            '2jax', ## SER:OG-ATP:O3G
            '1o3a', ## SER:OG-780:O6'
            '1o3m', ## SER:OG-785:O6'
            '1o2m', ## SER:OG-762:O6'
            '1ghx', ## SER:OG-BMZ:O6'
            '1fbp', ## THR:OG1-AMP:O3P
            '2gsu', ## THR:OG1-AMP:O1P
            '1w3l', ## TYR:OH-OXZ:O5
            '1ykp', ## TYR:OH-DHB:O1
            '1ykl', ## TYR:OH-DHB:O4

            '2uyv', ## GLY:O-TLA:O4
            '4gsp', ## 3GP:O3P-SGP:O3'
            ## unexpected nitrogen-nitrogen connection
            '1r1c', ## HIS:NE2-REP:N2
            ## unexpected nitrogen-oxygen connection
            '1lfg', ## GLN:NE2-NAG:O7
            '1rdy', ## LYS:NZ-AMP:O3P

            '1g18', ## THR:N-ADP:O1B
            '5at1', ## ILE:O-CTP:N4
            ## unexpected carbon-carbon connection
            '2fxh', ## TRP:CH2-TYR:CE1
            '1odw', ## DMN-HPH
            '1x83', ## CYS:CB-SBH:C5
            ## unexpected carbon-sulfur/oxygen connection
            '3gss', ## EAA:C11-GTT:SG2
            '1jl0', ## GLU:C-SER:O
            '2bog', ## SGC:S4-BGC:C1
            '2vzm', ## MET:SD-NRB:C13
            '3c7f','3c7h', ## XYP:C5B-XYP:O2B
            '1md9', ## GLY:CA-AMP:O2P
            '1mxt', ## MET:CE-OXY:O1 ## oxy... rename residue???
            '2exk', ## XYS:O4-XYS:C5 (O4-C1 in 2d22)
            '2fxj', ## MET:SD-TYR:CE2
            '1vyr', ## TRP:CH2-TNF:O22
            '2jjj', ## TSM:S-DPH:CM
            ## unexpected connection
            '1dmr','3dmr', ## PGD-O
            '2v5s', ## multiple connections from 1 residue and 2 altlocs of another residue
            '2hh5', ## GNQ:F3-GNQ:F1
            '2vvh', ## CYS:SG-SO3:O2
            '1qov', ## BCL-BPH
            '4rcr', ## ALA-BPH
            '1av8', ## HIS-FEO
            '1j3y','1j3z', ## LYS-2FU
            '1q82', ## C-PPU
            '2al2', ## LYS:NZ-ALA:N (possibly incorrect, introduced during remediation)
            '1y21', ## CYS:SG-NO:N
            '2ric', ## GMH:O3-GMH:C3
####            '1k7d', ## SER-GRO
####            '2hlp', ## CONECT record for a standard peptide bond
##            '1run',
##            '2nzy','1fmu','2ea0','1k7d',##
            '1mpm', ## GLC:O2-GLC:O3
            '1xwq', ## XYP:O4B-XYP:C5B bond
            '1l6m', ## DBH:O6-DBH:O3
            '1k3w','2opf', ## PRO:N-PED:C1'
            '2bd8','2bd5','2bd7', ## ILE:C-SER:O,ILE:C-ARG:N
            '1njt', ## SER:O-CFT:C (add CFT to smallmolecule.py)
##            '1h2z', ## SIN
##            '1sbh','1sbi','1sbn','1yja','1yjb', ## CYS-HYD
##            '1nhs', ## CYS-CYO
##            '2al2', ## LYS:NZ-ALA:N
##            '3by4', ## 3CN
##            '2bd2','2bd4','2bd7', ## casomorphin...
##            '1run', ## CMP:O5'-SER:OG
            '2bcn', ## ILE:CD1-LEU:O (not mentioned in primary ref, only present in chain B, cytochrome C peroxidase...)
##            '1uyx', ## GLN:OE1,BGC:O1
            '2vav',
            '3by4',
            '3foo', ## PRO:O-PXX:NAL
            '1m18', ## PYB:C-ABU:N
            ## hydrogen bonds
            '1bwn', ## SER:OG-4IP:OPH hydrogen bond

##            ## unexpected but correct connections
##            '1b8f','1gkj','2nxy', ## GLY:N/ALA:C
##            '1ggk','1gge','1ggj', ## TYR:CB-HIS:ND1
            '1mae', ## HDZ-N6A
##            '1mwv','2fxg','2b2o','2b2r','2dv1', ## TYR:CE1-TRP:CH2
            '2h5o', ## chromophore
##            '1odw', ## DMN-HPH
            '1mtb','2fgu', ## DIQ:CM-HPH:C
##            '1w3p', ## PYR:O1-HIS:NE2 hydrogen bond
##            '1y21', ## NO (nitrogen oxide)
##            '2azc', ## PHE:c-PHE:C
            '1ca8', ## GLY-NVA
##            ## small molecule connection (e.g. pyruvate, glycerol, hydrogen peroxide PEO)
####            '1qs7','1zbg','1fff','1fg8','1k2b','1vwb', ## vs 1qtx,1zj7,1ffi (NH2)
####            '222l','228l','220l','227l','1iii', ## BME
##
####            ## REMARK 465 record and/or SEQRES record, no MODRES record
####            '1gct','2iff','1x9s','1x9w','1xkf','1wu1','1tk5','1s10',
####            '1skw','1sal','1ry1','1pts',
####            ## residue renumbering (not involving icodes? and therefore no errors)
####            '1t7b','1t7e','1tud','1yra','2qwa','2qwb','2qwc','2qwd','2qwe','2qwf','2qwg','2qwh','2qwi','2qwj','2qwk',

            ## spelling/syntax/integer mistake
            '2w8k','1gaz','2w63','2w5a','1ilz','1im0','2zfo','2bvw','2z2z',
            '2i6k','1fah','1c7p','3dsw','2jh2','1gnz',

            ## sequence shifted?
            '2gn5',

            ## incorrect chain ID
            '1ogg','1w48',

            ## incorrect atom name
            '2b5t','1gmp','1gmq','1gmr','1rsn','1sar','2sar',

            ## buffering agent is ligand that induces conformational change
            '1l6l','2ou1',

            ## domain swapping
            '1jwj',
            '3bcm','1r5d','3bcp',
            '1k51',
            '1tij','1g96',
            ]

############# manual selection
        if '-manual' in sys.argv:

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

            d_rmsd,d_header = self.analyze_pdbs(
                verbose=True,do_gif_and_pdb = True,
                bool_sphere = bool_sphere,
                )
    ##        d_rmsd = read_rmsd_from_file()
##            self.write_html(d_rmsd, d_header, prefix='rmsd')
            print 'manual finished'
            return
    ##########

################### singlemutants
        if '-rerun' in sys.argv:
            fd = open('pisa3.old','r')
            lines = fd.readlines()
            fd.close()
            for i in range(len(lines)):
                if i < int(sys.argv[-1]):
                    continue
                print '@@@@@@@@@@@@@@@@@@@', i, len(lines)
                line = lines[i]
                self.l_pdbs = line.split()[:2]
                if '2vfn' in self.l_pdbs:
                    continue
                if len(set(self.l_pdbs)-set(errorpdbs)) < 2:
                    continue
                d_rmsd,d_header = self.analyze_pdbs(verbose=False,)
                self.write_html(d_rmsd, d_header, prefix='rmsd')
            stop_rerun
#########################

################### singlemutants
        if '-singlemutants' in sys.argv:
            ## 1) find all single mutants in single_mutants txt files
            ## 2) read html files of set of single mutants and find all non-mutants
            ## 3) calculate sphere rmsd for non-mutants and single mutants
            lines = []
            for s in '0123456789abcddefghijklmnopqrstuvwxyz':
                fd = open('single_point_mutations/%s.txt' %(s),'r')
                lines += fd.readlines()
                fd.close()
            for i in range(len(lines)):
                if i < int(sys.argv[-1]):
                    continue
                print '@@@@@@@@@@@@@@@@@@@', i, len(lines)
                line = lines[i]
                self.l_pdbs = line.split()[:2]
                if len(set(self.l_pdbs)-set(errorpdbs)) < 2:
                    continue
                d_rmsd,d_header = self.analyze_pdbs(verbose=False)
                self.write_html(d_rmsd, d_header, prefix='rmsd')
            stop_single_mutants
#########################

################### html
        if '-html' in sys.argv:

            path = 'htm/'
            htmls = os.listdir(path)
            htmls.sort()
            l_pdbs = []
            for i in range(len(htmls)):
                if i < int(sys.argv[-1]):
                    continue
                print i
                html = htmls[i]
                print html
                pdb = html[:4]
                if i % 100 == 0:
                    print '%4.1f%%' %(100*float(i)/len(htmls))
                fd = open('%s%s' %(path,html),'r')
                lines = fd.readlines()[2:-1]
                fd.close()

                ## this looks like a mess...
                for j in range(8+len(self.l_columns_html)+2,len(lines),len(self.l_columns_html)+2,): ## 2 for tr and /tr, 2 for table and jscript
                    line = lines[j]
                    n_mutations = int(line[line.index('>')+1:line.rindex('<')])
                    res1 = lines[j+8][lines[j+8].index('>')+1:lines[j+8].rindex('<')].strip()
                    res2 = lines[j+9][lines[j+9].index('>')+1:lines[j+9].rindex('<')].strip()
                    bm1 = int(lines[j-3][lines[j-3].index('>')+1:lines[j-3].rindex('<')])
                    bm2 = int(lines[j-2][lines[j-2].index('>')+1:lines[j-2].rindex('<')])
                    if res1 == 'N/A':
                        continue
                    if res2 == 'N/A':
                        continue
                    if n_mutations == 1 and bm1 == 1 and bm2 == 2 and float(res1) < 2. and float(res2) < 2.: ## bm2 = 2???
                        pdb1 = lines[j-5][20:24]
                        pdb2 = lines[j-4][20:24]
                        l_pdbs += [pdb1]
                        l_pdbs += [pdb2]
                        l_pdbs = [pdb1,pdb2,]
                        if pdb1 in errorpdbs or pdb2 in errorpdbs:
                            continue
                        self.l_pdbs = l_pdbs
                        d_rmsd,d_header = self.analyze_pdbs(verbose=False)
                        self.write_html(d_rmsd, d_header, prefix='rmsd')


#######################
        
############ out (use previous pairs as input)
        if '-out' in sys.argv:

            sets_pdbs = []
            self.l_pdbs = set()

            dirs = os.listdir('%s/pdb/' %(topdir))
            for dir in dirs:
                files = os.listdir('%s/pdb/%s' %(topdir,dir))
                for file in files:

##                    if sys.argv[-1] not in file:
##                        continue
##                    print '***', file

                    pdb1 = file[0:4]
                    pdb2 = file[6:10]
                    added = False
                    for i in range(len(sets_pdbs)):
                        set_pdb = sets_pdbs[i]
                        if pdb1 in set_pdb or pdb2 in set_pdb:
                            sets_pdbs[i] |= set([pdb1,pdb2])
                            added = True
                    if added == False:
                        sets_pdbs += [set([pdb1,pdb2])]
                    self.l_pdbs |= set([pdb1]) ##
                    self.l_pdbs |= set([pdb2]) ##

    ##        for pdb in self.l_pdbs:
    ##            os.remove('/oxygenase_local/tc/quakes/htm/%s.htm' %(pdb))

            self.l_pdbs = list(self.l_pdbs)
            for set_pdb in sets_pdbs:
                self.l_pdbs = list(set_pdb)

                for pdb in errorpdbs:
                    try:
                        self.l_pdbs.remove(pdb)
                    except:
                        None
                d_rmsd,d_header = self.analyze_pdbs(verbose=False)
                self.write_html(d_rmsd, d_header, prefix='rmsd')
##############

        if '-sphere' in sys.argv:

            if 1 == 1:

                print 'reading clusters'
                fd = open('pdbS95bF.out','r')
                lines = fd.readlines()
                fd.close()

                self.cluster = cluster = int(sys.argv[sys.argv.index('-cluster')+1])

                line = lines[cluster]
                l_pdbs = line.split()
                self.l_pdbs = []
                for pdb in l_pdbs:
                    self.l_pdbs += [pdb[:4].lower()]
                self.l_pdbs = list(set(self.l_pdbs))

                for pdb in errorpdbs:
                    if pdb in self.l_pdbs:
                        self.l_pdbs.remove(pdb)

                self.l_pdbs.sort()
                if '-skip' in sys.argv:
                    self.l_pdbs = self.l_pdbs[int(sys.argv[-1]):]

                if '1sf7' in self.l_pdbs:
                    print self.l_pdbs
                    stop

                d_rmsd,d_header = self.analyze_pdbs(do_gif_and_pdb = False, bool_sphere = True, bool_do_single_mutant = False)

            return
##########

                
############## clusters 95
        if '-clusters' in sys.argv:

##            os.system('wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/NR/clusters95.txt')

            print 'reading clusters'
            fd = open('pdbS95bF.out','r')
            lines = fd.readlines()
            fd.close()

##            l_pdbsets = []
##            for i in range(len(lines)):
##                if i % 1000 == 0:
##                    print i
##                xpdbs = []
##                line = lines[i]
##                bool_prev = False
##                for s in line.split():
##                    pdb = s[:4]
##                    for j in range(len(l_pdbsets)):
##                        if pdb in l_pdbsets[j]:
##                            bool_prev = True
##                            break
##                    xpdbs += [pdb]
##                if bool_prev == True:
##                    l_pdbsets[j] |= set(pdbs)
##                else:
##                    l_pdbsets += [set(pdbs)]
##            print len(l_pdbsets), len(lines)
##
##            s = ''
##            for i in range(len(l_pdbsets)):
##                for pdb in l_pdbsets:
##                    s += '%s ' %(pdb)
##                s += '\n'
##            fd = open('clusters.txt','w')
##            fd.write(s)
##            fd.close()
##            stop
##
##            fd = open('clusters.txt','r')
##            lines = fd.readlines()
##            fd.close()
                
            print 'finished reading clusters'
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

##                ## skip if *not* errorpdb
##                if len(set(pdbs)&set(errorpdbs)) == 0:
##                    continue

    ##            print cluster, pdbs
    ##            stop

                self.l_pdbs.sort()

                ## skip if pdb not rsynced
                for pdb in list(self.l_pdbs):
                    if not os.path.isfile('%s/%s/pdb%s.ent' %(self.path_pdb,pdb[1:3],pdb)):
                        self.l_pdbs.remove(pdb)

                ## skip if error pdb
                for pdb in errorpdbs:
                    try:
                        self.l_pdbs.remove(pdb)
                    except:
                        None

                ## skip if only one pdb in cluster
                if len(self.l_pdbs) <= 1:
                    continue

##                ## how far have we come?
##                if '1gby' not in self.l_pdbs:
##                    continue
##                else:
##                    print 'cluster', cluster
##                    stop
                    
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

        ##
        ## description of the dataset
        ##
        self.dataset_description()

        ##
        ## analysis of phi/psi changes due to mutations
        ##
        self.plot_phipsi()

        ##
        ## analysis of rmsd
        ##
        self.analyze_rmsd(conversion = False)

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
                d_header[pdb1]  = self.parse_header(pdb1)
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
                    d_header[pdb2]  = self.parse_header(pdb2)
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
                                continue
                            if not 'SSBOND' in d_header[pdb2].keys():
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
                            d_coordinates[pdb1], d_hetero[pdb1], d_ATOMseq[pdb1] = self.parse_coordinates(
                                pdb1, d_header[pdb1], verbose=verbose
                                )
                            ## append secondary structure
                            d_ATOMseq[pdb1] = self.append_ss(d_header[pdb1],d_ATOMseq[pdb1],)
                        if pdb2 not in d_coordinates.keys():
                            d_coordinates[pdb2], d_hetero[pdb2], d_ATOMseq[pdb2] = self.parse_coordinates(
                                pdb2, d_header[pdb2], verbose=verbose
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
                                self.phipsi(
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

##                        ## calculate solvent accessible surface area
##                        d_coordinates[pdb1] = self.whatif(pdb1,biomolecule1,d_coordinates[pdb1])

                        d_biomolecules = {
                            pdb1:{'biomolecule':biomolecule1},
                            pdb2:{'biomolecule':biomolecule2},
                            }

                        ## color code structure by rmsd
                        if do_gif_and_pdb == True:
                            self.rmsd2bfactor(
                                pdb1, pdb2, biomolecule1, biomolecule2, rmsd,
                                d_coordinates, d_header, tv1, rm, tv2, l_equivalent_chains,
                                bmchains1, bmchains2,
                                d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
                                d_ATOMseq,
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

        l_atoms = ['CA']
        l_atoms = []
        if l_atoms == []:
            s_atoms = 'heavy'
        elif l_atoms == ['CA']:
            s_atoms = 'CA'

        chains1 = l_equivalent_chains[0]
        chains2 = l_equivalent_chains[1]
        if len(chains1) != 1 or len(chains2) != 1:
            print chains1
            print chains2
            stop
        chain1 = chains1[0]
        chain2 = chains2[0]
        rep_chain1 = self.chain2repchain(pdb1,chain1[0],d_chains_intrapdb_sequence_identical,)
        rep_chain2 = self.chain2repchain(pdb2,chain2[0],d_chains_intrapdb_sequence_identical,)
        l1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1']
        l2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2']

        (
            coordinates1, coordinates2, residue_count, d_lines, l_RMSDs,
            d_coordinates1, d_coordinates2,
            ) = self.prepare_coords_for_alignment(
                pdb1,pdb2,chains1,chains2,
                d_chains_intrapdb_sequence_identical,
                d_chains_interpdb_sequence_similar,
                d_coordinates,d_ATOMseq,
                l_atoms=l_atoms,
                rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2,
                bool_None_if_unobs = True,
                bool_return_transformed_coordinates = True,
                )

        l_radii = [5,10,20,50,100,]
        
        ##
        ## HEWL
        ##
##        if len(d_ATOMseq[pdb1][chain1[0]]['seq']) != len(d_header['2lzm']['SEQRES']['chains']['A']['seq']) or len(d_ATOMseq[pdb2][chain2[0]]['seq']) != len(d_header['2lzm']['SEQRES']['chains']['A']['seq']):
##            stop
##            return
##
##        if len(coordinates1) != len(d_header['2lzm']['SEQRES']['chains']['A']['seq']) or len(coordinates2) != len(d_header['2lzm']['SEQRES']['chains']['A']['seq']):
##            print len(coordinates1)
##            print len(coordinates2)
##            stop

        ## sequences of different length
        if len(d_header[pdb1]['SEQRES']['chains'][chain1[0]]['seq']) != len(d_header[pdb2]['SEQRES']['chains'][chain2[0]]['seq']):
            return

        if prefix == 'mutant':
            res_index1 = d_mutations[pdb1].values()[0][0][0]
            res_index2 = d_mutations[pdb2].values()[0][0][1]
            if res_index1 != res_index2:
                print res_index1, res_index2
                print l1,l2
                stop
            l_res_indexes = [res_index1,]
        else:
            l_res_indexes = range(len(d_header[pdb1]['SEQRES']['chains'][chain1[0]]['seq']))

        d_rmsd = {}

        l = []
        print len(coordinates1), len(coordinates2)
        print d_header[pdb1]['SEQRES']['chains'][chain1[0]]['seq']
        print len(d_header[pdb2]['SEQRES']['chains'][chain2[0]]['seq'])
        for res_index_ref in l_res_indexes:
            d_rmsd[res_index_ref] = {}
            coord_ref1 = coordinates1[res_index_ref]
            coord_ref2 = coordinates2[res_index_ref]

            ## unobserved residue
            if coord_ref1 == None:
                continue
            if coord_ref2 == None:
                continue

            coord_ref = (coord_ref1+coord_ref2)/2.
            d_sqdiff = {}
            for res_index in range(len(d_header[pdb1]['SEQRES']['chains'][chain1[0]]['seq'])):

##                for x in [1]:
##                    coord1 = coordinates1[res_index]
##                    coord2 = coordinates2[res_index]

                ## residue missing
                if not res_index in d_coordinates1[chain1].keys():
                    continue
                if not res_index in d_coordinates2[chain2].keys():
                    continue
                
                for atom_name in d_coordinates1[chain1][res_index].keys():
                    ## only heavy atoms
                    if atom_name[0] == 'H':
                        continue
                    coord1 = d_coordinates1[chain1][res_index][atom_name]
                    coord2 = d_coordinates2[chain2][res_index][atom_name]

                    dist1 = math.sqrt(sum((coord_ref-coord1)**2))
                    dist2 = math.sqrt(sum((coord_ref-coord2)**2))
                    dist = (dist1+dist2)/2.
                    diff_squared = sum((coord1-coord2)**2)
                    for r in l_radii:
                        if dist < r:
                            if not r in d_sqdiff.keys():
                                d_sqdiff[r] = []
                            d_sqdiff[r] += [diff_squared]
                            break
                        ## anything more than 50Angstrom away?
                        if r == l_radii[-1]:
                            print dist
                            stop
            s_rmsd = ''
            for r in l_radii:
                if not r in d_sqdiff.keys():
                    s_rmsd += 'None '
                    continue
                l_sqdiff = d_sqdiff[r]
                RMSD = math.sqrt(sum(l_sqdiff)/len(l_sqdiff))
                s_rmsd += '%f ' %(RMSD)
                if bm1 > 9 or bm2 > 9:
                    print bm1,bm2
                    stop_increase_two_digits
            if (
##                d_header[pdb1]['AUTHOR'] == d_header[pdb2]['AUTHOR'] ## all authors identical
##                or
                len( set(d_header[pdb1]['AUTHOR']) & set(d_header[pdb2]['AUTHOR']) ) >= 1 ## 1 or more identical authors (by name)
                ):
                bool_authors_identical = True
            else:
                bool_authors_identical = False
            l += ['%4s %4s %1i %1i %1s %1s %4i %s %s\n' %(pdb1,pdb2,bm1,bm2,chain1,chain2,res_index_ref,s_rmsd,bool_authors_identical,)]

##        fd = open('HEWL_%s.txt' %(prefix,),'a') ## mutant or wt
        fd = open('sphere/%i_%s_%s.txt' %(self.cluster,s_atoms,prefix,),'a') ## mutant or wt
        fd.writelines(l)
        fd.close()

        return

    def dataset_description(self):

        print 'make sure nothing is written to html files. continue?'
        s = raw_input()
        if s != 'y':
            return


        ##
        ## table of single point mutations
        ##
        print 'generating table with statistics for single point mutations'
        lines = []
        for s in '0123456789abcddefghijklmnopqrstuvwxyz':
            fd = open('single_point_mutations/%s.txt' %(s),'r')
            lines += fd.readlines()
            fd.close()
        d_counts = {}
        d_counts_wt = {}
        d_counts_mut = {}
        count = 0
        ## temp!!! tmp!!!
        for res in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]: ## temp!!! tmp!!!
            d_counts_wt[res] = 0 ## temp!!! tmp!!!
            d_counts[res] = 0 ## temp!!! tmp!!!
            d_counts_mut[res] = 0 ## temp!!! tmp!!!
        for line in lines:
            res_name1 = line.split()[8]
            res_name2 = line.split()[9]
            if res_name1 == 'N/A':
                continue
            if res_name2 == 'N/A':
                continue
            res = res_name1+res_name2
            if not res in d_counts.keys():
                d_counts[res] = 0
            d_counts[res] += 1
            count += 1
            if not res_name1 in d_counts_wt.keys():
                d_counts_wt[res_name1] = 0
            d_counts_wt[res_name1] += 1
            if not res_name2 in d_counts_mut.keys():
                d_counts_mut[res_name2] = 0
            d_counts_mut[res_name2] += 1
        print 'total number of single point mutations:', count
        s = '   \t'
        for res in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
            s += '%5s\t' %(res)
        s += 'total\n'
        for res1 in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
            s += '%5s\t' %(res1)
            for res2 in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
                if res1 == res2:
                    s += '\t'
                elif not res1+res2 in d_counts.keys():
                    s += '\t'
                else:
                    s += '%5i\t' %(d_counts[res1+res2])
            s += '%5i' %(d_counts_wt[res1])
            s += '\n'
        s += 'total\t'
        for res2 in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
            s += '%5i\t' %(d_counts_mut[res2])
        s += '%5i' %(count)
        s += '\n'
        fd = open('mutation_table.txt','w')
        fd.write(s)
        fd.close()
##        stop


        ##
        ## statistics
        ##
        print 'generating statistics'

        htmls = os.listdir('htm/')
        n_lines = 1+len(self.l_columns_html)+1
        d_values = {
            'rmsd':[],
            'chains':{},
            'residues':[],
            'coordinates':[],
            'mutations':{},
            }
        d_parameters = {
            'rmsd':7,
            'mutations':8,
            'chains':9,
            'residues':10,
            'coordinates':11,
            }
        for j in range(len(htmls)):
            if j % 100 == 0:
                print '%4.1f%%' %(100*float(j)/len(htmls))
            html = htmls[j]
            pdb = html[:4]
            fd = open('htm/%s' %(html),'r')
            lines = fd.readlines()
            fd.close()
            if lines[0][:8] == '<script ':
                lines = lines[2+n_lines:-1]
            else:
                lines = lines[1+n_lines:-1]
            for i in range(len(lines)/(n_lines)):
                for parameter in ['mutations','residues','coordinates','chains','rmsd',]:
                    i_line = d_parameters[parameter]
                    s_line = lines[i*n_lines+i_line]
                    index2 = s_line.rindex('</td>')
                    index1 = s_line[:index2].rindex('>')+1
                    value = s_line[index1:index2].strip()
                    if parameter == 'rmsd':
                        if float(value) < 10.:
                            d_values['rmsd'] += [float(value)]
                        else:
                            break
                    if parameter == 'chains':
                        chains = int(value)
                        if not int(value) in d_values['chains'].keys():
                            d_values['chains'][chains] = 0
                        d_values['chains'][chains] += 1
                        if chains == 1:
                            if not mutations in d_values['mutations'].keys():
                                d_values['mutations'][mutations] = 0
                                if mutations >= 60:
                                    print html, htmls[j-1]
                                    stop
                            d_values['mutations'][mutations] += 1
                            d_values['residues'] += [int(residues)]
                            d_values['coordinates'] += [coordinates]
                    if parameter == 'mutations':
                        mutations = int(value)
                    if parameter == 'residues':
                        residues = int(value)
                    if parameter == 'coordinates':
                        coordinates = int(value)
        print 'avg rmsd', sum(d_values['rmsd'])/len(d_values['rmsd'])
        print 'avg residues', sum(d_values['residues'])/len(d_values['residues'])
        print 'avg coordinates', sum(d_values['coordinates'])/len(d_values['coordinates'])
        for k,v in d_values['chains'].items():
            print 'chains', k,v
        for k,v in d_values['mutations'].items():
            print 'mutations', k,v


        ## number of paired PDBs
        l_fn = os.listdir('%s/htm' %(self.topdir))
        print 'number of paired *pdb files*:', len(l_fn)

        ## protein with most partners
        l_fn = os.listdir('%s/htm' %(self.topdir))
        size_max = [0,'N/A',]
        for fn in l_fn:
            if os.path.getsize('%s/htm/%s' %(self.topdir,fn,)) > size_max[0]:
                size_max = [os.path.getsize('%s/htm/%s' %(self.topdir,fn,)),fn,]
        fd = open('%s/htm/%s' %(self.topdir,size_max[1],),'r')
        lines = fd.readlines()
        fd.close()
        i = -1
        for line in lines:
            if '<tr>' in line:
                i += 1
        print '%s has %i partners (incl. structures with multiple biounits)' %(size_max[1],i,)

        ## number of paired biomolecules
        n_pairs = 0
        n_biomolecules = 0
        for s in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
            l_fn = os.listdir('%s/pdb/%s' %(self.topdir,s.lower(),))
            n_pairs += len(l_fn)
            set_biomolecules = set()
            for s in l_fn:
                set_biomolecules |= set([s[:6]])
                set_biomolecules |= set([s[6:12]])
            n_biomolecules += len(set_biomolecules)
        print 'number of *biomolecules*', n_biomolecules
        print 'number of paired *biomolecules*', n_pairs

        return


    def plot_phipsi(self,):

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
                print line
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
            continue ## tmp!!!
            if os.path.isfile('plot_phipsi_%s.png' %(key)):
                os.remove('plot_phipsi_%s.png' %(key))
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

        return


    def round_angle(self,angle,phipsi_step,):

        if angle == 180.:
            angle = -180.
        else:
            angle = phipsi_step*int(angle/phipsi_step)

        return angle


    def analyze_rmsd(self, conversion = True):

        print 'analyze rmsd'

        d_data_ratio_discrete = {
            'mutations':{},
            'chains':{},
            'residues':{},
            'coordinates':{},
            }
        d_data_ratio_continuous = {
            'pH':{'single':[],'min':[],'max':[],'average':[],'difference':[],}, ## discrete?
            'T':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
            'res':{'single':[],'min':[],'max':[],'average':[],'difference':[],}, ## discrete?
            }
        d_data_nominal = {
            'spacegroup':{'identical':{},'different':{}},
            'hetIDs':{'identical':[],'different':[]},
            'REMARK465':{'True':[],'False':[]},
            'REMARK470':{'True':[],'False':[]},
            'transformations':{'True':[],'False':[]},
            'bool_identical_authors':{'True':[],'False':[]},
            }
        d_combined = {
            'chains_no_transformations':{},
            'residues_no_transformations':{},
            'mutations_no_transformations':{},
            'residues_1_chain':{},
            'mutations_1_chain':{},
            'pHdiff_nomutation':{},
            }

        terminal = 'postscript'

        n_lines = 1+len(self.l_columns_html)+1

        htmls = os.listdir('htm/')
        htmls.sort()
        lines = []
        d_pdbs_skip = {}
        l_gnuplot = []

        for i in range(len(htmls)):
            html = htmls[i]
            pdb = html[:4]
            if i % 100 == 0:
                print '%4.1f%%' %(100*float(i)/len(htmls))
            fd = open('htm/%s' %(html),'r')
            lines = fd.readlines()
            fd.close()
            if lines[0][:8] == '<script ':
                lines = lines[2+n_lines:-1]
            else:
                lines = lines[1+n_lines:-1]
            

            ##
            ## parse data
            ##
            (
                d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,set_pdbs,
                d_combined,
                l_gnuplot,
                ) = self.parse_html(
                    lines,n_lines,d_pdbs_skip,
                    d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,
                    d_combined,
                    html,
                    l_gnuplot,
                    )

        fd = open('data.gnu','w')
        fd.writelines(l_gnuplot)
        fd.close()

##        ##
##        ## plot continuous ratio scale data
##        ##
##        d_steps = {
##            'pH':{'xstep':0.1,'ystep':0.1,},
##            'T':{'xstep':5.0,'ystep':0.1,},
##            'res':{'xstep':0.1,'ystep':0.1,},
##            }
##        for parameter in d_data_ratio_continuous.keys():
##            step_x = d_steps[parameter]['xstep']
##            step_y = d_steps[parameter]['ystep']
##            for type in d_data_ratio_continuous[parameter]:
##                gnuplotdata = ''
##                d_discrete = {}
##                rmsd_discrete_max = 0
##                for values in d_data_ratio_continuous[parameter][type]:
##                    value = values[0]
##                    rmsd = values[1]
##                    gnuplotdata += '%s %s\n' %(value, rmsd)
##                    value_discrete = value-value%step_x
##                    rmsd_discrete = rmsd-rmsd%step_y
##                    if rmsd_discrete > rmsd_discrete_max:
##                        rmsd_discrete_max = rmsd_discrete
##                    if not value_discrete in d_discrete.keys():
##                        d_discrete[value_discrete] = {}
##                    if not rmsd_discrete in d_discrete[value_discrete].keys():
##                        d_discrete[value_discrete][rmsd_discrete] = 0
##                    d_discrete[value_discrete][rmsd_discrete] += 1
##                l_gnuplotdata_contour = []
##                for value_discrete in range(0,int(max(d_discrete.keys())/step_x)+1,):
##                    value_discrete *= step_x
##                    for rmsd_discrete in range(0,int(rmsd_discrete_max/step_y),):
##                        rmsd_discrete *= step_y
##                        count = 0
##                        if value_discrete in d_discrete.keys():
##                            if rmsd_discrete in d_discrete[value_discrete].keys():
##                                count = d_discrete[value_discrete][rmsd_discrete]
##                        l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete, rmsd_discrete, count,)]
####        gnuplot_splot_data.append('%4i %4i %16.13f\n\n' %(x+1, len(data[x])+1, data[0][0]))
####    for i in range(len(data)+1):
####        gnuplot_splot_data.append('%4i %4i %16.13f\n' %(len(data)+1, i+1, data[0][0]))
##                    l_gnuplotdata_contour += ['\n']
##
##                prefix_gnuplot = 'ps/%s%s' %(parameter,type,)
##                fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
##                fd.write(gnuplotdata)
##                fd.close()
##                ## linear regression
##                if type in ['difference','single',] or parameter == 'res':
##                    regression = True
##                else:
##                    regression = False
##                gnuplot.scatter_plot_2d(
##                    prefix_gnuplot, regression=regression, errorbars=True, terminal=terminal, xlabel = parameter,
##                    )
##                gnuplot.contour_plot(
##                    prefix_gnuplot, l_gnuplotdata_contour,
##                    xlabel = parameter+type, ylabel = 'RMSD', title = '%s%s v RMSD' %(parameter,type,),
##                    )
##                if conversion == True:
##                    print '***convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
##                    os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
##                    if os.path.isfile('%s.ps' %(prefix_gnuplot)):
##                        os.remove('%s.ps' %(prefix_gnuplot))
##                print prefix_gnuplot
##
##
##        ##
##        ## plot discrete ratio scale data
##        ##
##        for parameter in d_data_ratio_discrete.keys():
##
##            gnuplotdata = ''
##            for value in d_data_ratio_discrete[parameter].keys():
##                for rmsd in d_data_ratio_discrete[parameter][value]:
##                    gnuplotdata += '%s %s\n' %(value, rmsd)
##
##            prefix_gnuplot = 'ps/%s' %(parameter)
##            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
##            fd.write(gnuplotdata)
##            fd.close()
##
##            if parameter in ['residues','coordinates']:
##                logarithmic = True
##            else:
##                logarithmic = False
##            gnuplot.scatter_plot_2d(
##                prefix_gnuplot, regression=True, logarithmic=logarithmic, errorbars=True, terminal=terminal, xlabel = parameter
##                )
##            if conversion == True:
##                print '***convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
##                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
##                os.remove('%s.ps' %(prefix_gnuplot))
##
##
##        ####
##        ## nominal scale data
##        ####
##
##        ##
##        ## hetIDs, remarks, transformations
##        ##
##        d_nominal = {
##            'hetIDs':['identical','different'],
##            'REMARK465':['True','False'],
##            'REMARK470':['True','False'],
##            'transformations':['True','False'],
##            }
##        for parameter in d_nominal.keys():
##            l1 = d_data_nominal[parameter][d_nominal[parameter][0]]
##            l2 = d_data_nominal[parameter][d_nominal[parameter][1]]
####            gnuplotdata = ''
####            for rmsd in l1:
####                gnuplotdata += '1 %s\n' %(rmsd)
####            for rmsd in l2:
####                gnuplotdata += '2 %s\n' %(rmsd)
####            prefix_gnuplot = 'ps/%s' %(parameter)
####            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
####            fd.write(gnuplotdata)
####            fd.close()
####            d_xtics = {d_nominal[parameter][0]:1.,d_nominal[parameter][1]:2.}
####            gnuplot.scatter_plot_2d(prefix_gnuplot, d_xtics=d_xtics, terminal=terminal, xlabel = parameter)
####            if conversion == True:
####                print '***convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
####                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
####                os.remove('%s.ps' %(prefix_gnuplot))
##
##            print parameter
##            statistics.twosamplettest(l1,l2)
##
##
##        ##
##        ## space groups
##        ##
##        ## statistics
##        l1 = []
##        l2 = []
##        ## xtics
##        set_spacegroups = set()
##        set_spacegroups |= set(d_data_nominal['spacegroup']['identical'].keys())
##        set_spacegroups |= set(d_data_nominal['spacegroup']['different'].keys())
##        l_spacegroups = list(set_spacegroups)
##        l_spacegroups.sort()
##        d_spacegroups = {}
##        d_xtics = {}
##        for i in range(len(l_spacegroups)):
##            spacegroup = l_spacegroups[i]
##            d_spacegroups[spacegroup] = float(i)
##            d_xtics[spacegroup] = float(i+.5)
##        ## identical, different
##        for type in ['identical','different',]:
##            prefix_gnuplot = 'ps/spacegroups%s' %(type)
##            ## parse gnuplotdata
##            gnuplotdata = ''
##            ## discrete
##            d_discrete = {}
##            rmsd_discrete_max = 0
##            step_x = 1.0
##            step_y = 0.1
##            for spacegroup in d_data_nominal['spacegroup'][type].keys():
##                xval = d_spacegroups[spacegroup]
##                ## data for plot
##                for rmsd in d_data_nominal['spacegroup'][type][spacegroup]:
##                    ## non-discrete
##                    gnuplotdata += '%s %s\n' %(xval, rmsd)
##                    ## discrete
##                    value_discrete = xval
##                    rmsd_discrete = rmsd-rmsd%step_y
##                    if rmsd_discrete > rmsd_discrete_max:
##                        rmsd_discrete_max = rmsd_discrete
##                    if not value_discrete in d_discrete.keys():
##                        d_discrete[value_discrete] = {}
##                    if not rmsd_discrete in d_discrete[value_discrete].keys():
##                        d_discrete[value_discrete][rmsd_discrete] = 0
##                    d_discrete[value_discrete][rmsd_discrete] += 1
##                ## data for statistics
##                if type == 'identical':
##                    l1 += d_data_nominal['spacegroup'][type][spacegroup]
##                elif type == 'different':
##                    l2 += d_data_nominal['spacegroup'][type][spacegroup]
##
##            l_gnuplotdata_contour = []
####            for value_discrete in range(0,int(max(d_discrete.keys())/step_x)+1,):
##            for value_discrete in range(0,int((len(l_spacegroups)+1)/step_x)+1,):
##                value_discrete *= step_x
##                for rmsd_discrete in range(0,int(rmsd_discrete_max/step_y),):
##                    rmsd_discrete *= step_y
##                    count = 0
##                    if value_discrete in d_discrete.keys():
##                        if rmsd_discrete in d_discrete[value_discrete].keys():
##                            count = d_discrete[value_discrete][rmsd_discrete]
##                    l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete, rmsd_discrete, count,)]
####        gnuplot_splot_data.append('%4i %4i %16.13f\n\n' %(x+1, len(data[x])+1, data[0][0]))
####    for i in range(len(data)+1):
####        gnuplot_splot_data.append('%4i %4i %16.13f\n' %(len(data)+1, i+1, data[0][0]))
##                l_gnuplotdata_contour += ['\n']
##
##            ## write gnuplotdata
##            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
##            fd.write(gnuplotdata)
##            fd.close()
##            ## plot 2d scatter (non-discrete)
##            gnuplot.scatter_plot_2d(
##                prefix_gnuplot, d_xtics = d_spacegroups, xlabel='spacegroups', terminal=terminal,
##                )
##            ## plot 2d contour (discrete)
##            gnuplot.contour_plot(
##                prefix_gnuplot, l_gnuplotdata_contour,
##                xlabel = 'spacegroups', ylabel = 'RMSD', title = '%s%s v RMSD' %('spacegroups',type,),
##                d_xtics = d_xtics,
##                x2 = len(l_spacegroups),
##                )
##            if conversion == True:
##                print '***convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
##                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
##                os.remove('%s.ps' %(prefix_gnuplot))
##
##        ## statistics
##        print prefix_gnuplot
##        statistics.twosamplettest(l1,l2)


        ##
        ## plot combined discrete date
        ##
        for parameter in d_combined.keys():

            if parameter != 'mutations_1_chain': ## tmp!!!
                continue

            ## non-discrete data
            s_averages = ''
            gnuplotdata = ''
            for value in d_combined[parameter].keys():
                for rmsd in d_combined[parameter][value]:
                    gnuplotdata += '%s %s\n' %(value, rmsd)
                average = sum(d_combined[parameter][value])/len(d_combined[parameter][value])
                s_averages += '%s %s\n' %(value,average,)

            ## discrete data
            step_x = 1
            step_y = 0.1
            if parameter == 'mutations_1_chain':
                d_discrete = {}
                rmsd_discrete_max = 0
                for value in d_combined[parameter].keys():
                    value_discrete = value
                    for rmsd in d_combined[parameter][value]:
                        rmsd_discrete = rmsd-rmsd%step_y
                        if rmsd_discrete > rmsd_discrete_max:
                            rmsd_discrete_max = rmsd_discrete
                        if not value_discrete in d_discrete.keys():
                            d_discrete[value_discrete] = {}
                        if not rmsd_discrete in d_discrete[value_discrete].keys():
                            d_discrete[value_discrete][rmsd_discrete] = 0
                        d_discrete[value_discrete][rmsd_discrete] += 1
                l_gnuplotdata_contour = []
                for value_discrete in range(0,int(max(d_discrete.keys())/step_x)+1,):
                    value_discrete *= step_x
                    for rmsd_discrete in range(0,int(rmsd_discrete_max/step_y),):
                        rmsd_discrete *= step_y
                        count = 0
                        if value_discrete in d_discrete.keys():
                            if rmsd_discrete in d_discrete[value_discrete].keys():
                                count = d_discrete[value_discrete][rmsd_discrete]
                        l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete, rmsd_discrete, count,)]
                    l_gnuplotdata_contour += ['\n']

                fd = open('mutations_1_chain_discrete.gnuplotdata','w')
                fd.writelines(l_gnuplotdata_contour)
                fd.close()

            prefix_gnuplot = 'ps/%s' %(parameter)
            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata)
            fd.close()

            fd = open('ps/%s_averages.gnuplotdata' %(parameter),'w')
            fd.write(s_averages)
            fd.close()

            ## plot 2d contour (discrete)
            if parameter == 'mutations_1_chain':
                prefix_gnuplot = 'ps/mutations_1_chain_discrete'
                gnuplot.contour_plot(
                    prefix_gnuplot, l_gnuplotdata_contour,
                    xlabel = 'mutations 1 chain', ylabel = 'RMSD', title = '%s v RMSD' %('mutations 1 chain',),
                    x2 = 10+1,
                    )
            elif not parameter[parameter.index('_'):] in ['chains','residues',]:
                gnuplot.scatter_plot_2d(
                    prefix_gnuplot, regression=True, errorbars=True, terminal=terminal, xlabel = parameter.replace('_',' '),
                    averages = 'ps/%s_averages.gnuplotdata' %(parameter)
                    )
            else:
                gnuplot.scatter_plot_2d(
                    prefix_gnuplot, regression=True, errorbars=True, terminal=terminal, xlabel = parameter.replace('_',' '),
                    )
                
            if conversion == True:
                print '***convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
                print prefix_gnuplot
                os.remove('%s.ps' %(prefix_gnuplot))
        

        return


    def parse_html(
        self,lines,n_lines,d_pdbs_skip,
        d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,
        d_combined,
        html,
        l_gnuplot,
        ):

        ## loop over table rows
        for i in range(len(lines)/(n_lines)):

            s_gnuplot = ''

            ## parse pdbs
            pdb1_line = lines[i*n_lines+3]
            pdb2_line = lines[i*n_lines+4]
            pdbs = [pdb1_line,pdb2_line]
            for j in range(2):
                pdb_line = pdbs[j]
                try:
                    pdb_index2 = pdb_line.index('</a>')
                except:
                    print html
                    print pdbs
                    print pdb_line
                    print lines[i*n_lines:i*(n_lines+1)]
                    stop
                pdb_index1 = pdb_line[:pdb_index2].rindex('>')+1
                s_pdb = pdb_line[pdb_index1:pdb_index2]
                pdbs[j] = s_pdb
            pdb1 = pdbs[0]
            pdb2 = pdbs[1]

            if pdb1 < pdb2:
                if not pdb1 in d_pdbs_skip.keys():
                    d_pdbs_skip[pdb1] = []
                if pdb2 in d_pdbs_skip[pdb1]:
                    continue
                else:
                    d_pdbs_skip[pdb1] += [pdb2]
            elif pdb2 < pdb1:
                if not pdb2 in d_pdbs_skip.keys():
                    d_pdbs_skip[pdb2] = []
                if pdb1 in d_pdbs_skip[pdb2]:
                    continue
                else:
                    d_pdbs_skip[pdb2] += [pdb1]

            errorpdbs = set([
####                    ## error
####                    '2bfk','2bfl',
####                    ## less than 50 residues due to mutations
####                    '1p7i','1p7j',
##                ## less than 50 residues due to REMARK465 records
##                '2yxq','2yxr','1bx7','1bx8',
##                ## less than 50 residues due to non-overlapping residues
##                '1ft8','1koh',
##                ##
##                '1znb','1zlx','1i3h','1j5o','1f33','1qvc',
##                '1ypu','1uc7','1jr8','1g6l','1ixc','1xri',
##                '1xfi','2q40','2q47','2qd0','1m70',
##                ##
##                '1f1g','1f18','1abs','1dxd','1dvb',
##                '1dmm','1dmn','1dmq',
##                '195l','196l','200l','197l','199l','198l',
##                ## unknown whether error or not
##                '1za1','1qtv','1khh','1cvm','1ez9','1ulx','2der','1x98',
                ## correct (cryogen/cryo-cooled)
                '2fbb','1j3y','1j41','1ajg','1ajh','1a6k','1c1g','2i16','2qxw',
                '2j7n','3ecl','3e4n',
##                ## incorrect (celsius)
##                '1ade','1c0e',
##                ## remediation change
##                '1dos',
##                ## high pH
##                '3c90',
####                                '1fgn',
                ])

            if pdb1 in errorpdbs or pdb2 in errorpdbs:
                continue

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
                if parameter in [
                    'rmsd',
##                    'pH1','pH2','T1','T2','res1','res2',
                    ]:
                    try:
                        value = float(value)
                    except:
                        print parameter, value
                        print pdb1, pdb2
                        print html
                        stop
                if parameter in ['mutations','chains','residues','coordinates',]:
                    try:
                        value = int(value)
                    except:
                        print parameter, value
                        print pdb1,pdb2
                        print i
                        print lines[i*n_lines+i_line-1:i*n_lines+i_line+2]
                        print pdb1_line
                        print pdb2_line
                        stop
                d_parameters[parameter] = value

            ## skip if virus (might exclude other proteins...)
            if int(d_parameters['chains']) % 60 == 0:
                continue

            ## check resolution
            if len(set([pdb1,pdb2]) & set(['2c32','1ye1','2zqp',])) == 0:
                if d_parameters['res1'] != 'N/A':
                    if float(d_parameters['res1']) > 5.:
                        print pdb1,pdb2
                        stop_resolution
                if d_parameters['res2'] != 'N/A':
                    if float(d_parameters['res2']) > 5.:
                        print pdb2,pdb1
                        stop_resolution

            if float(d_parameters['residues']) < 50:
                errorpdbs = set([
####                    ## error
####                    '2bfk','2bfl',
####                    ## less than 50 residues due to mutations
####                    '1p7i','1p7j',
##                    ## less than 50 residues due to REMARK465 records
##                    '2yxq','2yxr','1bx7','1bx8',
##                    ## less than 50 residues due to non-overlapping residues
##                    '1ft8','1koh',
                    ])
                if len(set([pdb1,pdb2])-errorpdbs) > 0:
                    print set([pdb1,pdb2])-errorpdbs
                    print pdb1, pdb2
                    print d_parameters['residues']
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
                try:
                    index3 = s_line.index('</td>')
                except:
                    print pdb1, pdb2
                    print s_line
                    stop
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
            if float(rmsd) > 99.: ## temp!!! due to incorrect transformations...
                print rmsd, pdb1,pdb2
                continue
            if float(rmsd) > 5. and d_parameters['spacegroup1'] == d_parameters['spacegroup2']:
                print html, pdb1, pdb2, rmsd
                fd = open('tmp.txt','a')
                fd.write("'%s','%s,', " %(pdb1,pdb2,))
                fd.close()


            ######
            ## continuous ratio scale data
            ######
            d_diff = {}
            d_max = {}
            for parameter in d_data_ratio_continuous.keys():
                values = []
                for no in ['1','2']:
                    value = d_parameters[parameter+no]
                    if value not in ['NULL','N/A']:
                        if '-' in value: ## e.g. 4lzt pH
                            separator = '-'
                        if ';' in value:
                            separator = ';'
                        if '-' in value or ';' in value:
                            l_values = value.split(separator)
                            while 'NULL' in l_values:
                                    l_values.remove('NULL')
                            while ' NULL' in l_values:
                                    l_values.remove(' NULL')
                            value = 0
                            try:
                                for s_value in l_values:
                                    value += float(s_value)
                            except:
                                print pdb1, pdb2, s_value, l_values
                                stop
                            if len(l_values) == 0:
                                break
                            else:
                                value /= len(l_values)
##                            else:
##                                l_values = value.split(';')
##                                set_values = set()
##                                for value in l_values:
##                                    set_values |= set([value.strip()])
##                                    if len(set_values) > 1:
##                                        print pdb1, pdb2, no
##                                        print d_parameters['T'+no]
##                                        print d_parameters['pH'+no]
##                                        print set_values
##                                        stop
##                                    else:
##                                        value = list(set_values)[0]

                        try:
                            value = float(value)
                        except:
                            print value
                            print no, pdb1, pdb2
                            stop

                        if (value < 50 or value > 303) and parameter == 'T':
                            errorpdbs = set([
##                                '1f1g','1f18','1abs','1dxd','1dvb',
##                                '1dmm','1dmn','1dmq',
##                                '195l','196l','200l','197l','199l','198l',
##                                ## high temperature
##                                '1fah','1q0o',
##                                ## unknown whether error or not
##                                '1za1','1qtv','1khh','1cvm','1ez9','1ulx',
##                                '2der','1x98','1a3d','2v19','2zvz',
##                                ## correct (cryogen/cryo-cooled)
##                                '2fbb','1j3y','1j41','1ajg','1ajh','1a6k','1c1g','2i16','2qxw','2j7n','2j8t',
##                                ## incorrect (celsius)
##                                '1zin',
                                ])
                            if len(set([pdb1,pdb2])-errorpdbs) == 2:
                                print value, no
                                print pdb1, pdb2
                                stop2temperature ## temp!!!
                        if (value < 1 or value > 12) and value != 'NULL' and parameter == 'pH':
                            errorpdbs = set([
##                                ## high pH
##                                '3c90',
####                                '1fgn',
                                ])
                            if len(set([pdb1,pdb2])-errorpdbs) == 2:
                                print pdb1, pdb2
                                print 'pH', value
                                print d_parameters['pH1'], d_parameters['pH2']
                                stop3pH ## temp!!!
                        if value > 11.5 and parameter == 'res':
                            if len(set([pdb1,pdb2])-set(['1c1g','2tma',])) == 2:
                                print pdb1, pdb2
                                print d_parameters['res1'], d_parameters['res2']
                                stop4resolution ## temp!!

                        values += [value]
                        d_data_ratio_continuous[parameter]['single'] += [[value,rmsd]]


                if len(values) == 2:

                    ## check resolution
                    if parameter == 'res' and abs(values[0] - values[1]) > 4.61:
                        print abs(values[0] - values[1])
                        print pdb1,pdb2
                        stop

                    if parameter == 'pH' and abs(values[0] - values[1]) > 9.8:
                        print abs(values[0] - values[1])
                        print pdb1,pdb2
                        stop

                    d_data_ratio_continuous[parameter]['max'] += [[max(values),rmsd]]
                    d_data_ratio_continuous[parameter]['min'] += [[min(values),rmsd]]
                    d_data_ratio_continuous[parameter]['average'] += [[sum(values)/2.,rmsd]]
                    d_data_ratio_continuous[parameter]['difference'] += [[abs(values[1]-values[0]),rmsd]]
                    d_diff[parameter] = abs(values[1]-values[0])
                    d_max[parameter] = max(values)


            ######
            ## discrete ratio scale data
            ######
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
            binary_spacegroups_identical = None
            if spacegroup1 == spacegroup2:
                if not spacegroup1 in d_data_nominal['spacegroup']['identical'].keys():
                    d_data_nominal['spacegroup']['identical'][spacegroup1] = []
                d_data_nominal['spacegroup']['identical'][spacegroup1] += [rmsd]
                binary_spacegroups_identical = 1
            else:
                if not spacegroup1 in d_data_nominal['spacegroup']['different'].keys():
                    d_data_nominal['spacegroup']['different'][spacegroup1] = []
                d_data_nominal['spacegroup']['different'][spacegroup1] += [rmsd]
                if not spacegroup2 in d_data_nominal['spacegroup']['different'].keys():
                    d_data_nominal['spacegroup']['different'][spacegroup2] = []
                d_data_nominal['spacegroup']['different'][spacegroup2] += [rmsd]
                binary_spacegroups_identical = 0

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
                    if parameter == 'transformations':
                        binary_transformations = 1
                elif d_parameters[parameter] == 'False':
                    d_data_nominal[parameter]['False'] += [rmsd]
                    if parameter == 'transformations':
                        binary_transformations = 0


            ######
            ## combined data
            ######
            if d_parameters['transformations'] == 'False':
                chains = int(d_parameters['chains'])
                mutations = int(d_parameters['mutations'])
                for s in [
                    'chains_no_transformations',
                    'residues_no_transformations',
                    'mutations_no_transformations',
                    ]:
                    value = d_parameters[s[:s.index('_')]]
                    if value not in d_combined[s].keys():
                        d_combined[s][value] = []
                    d_combined[s][value] += [rmsd]
                if chains == 1:
                    for s in [
                        'residues_1_chain',
                        'mutations_1_chain',
                        ]:
                        value = d_parameters[s[:s.index('_')]]
                        if value not in d_combined[s].keys():
                            d_combined[s][value] = []
                        d_combined[s][value] += [rmsd]
                if mutations == 0:
                    if 'pH' in d_diff.keys():
                        value = d_diff['pH']
                        if value not in d_combined['pHdiff_nomutation'].keys():
                            d_combined['pHdiff_nomutation'][value] = []
                        d_combined['pHdiff_nomutation'][value] += [rmsd]

            if 'pH' in d_diff.keys() and 'T' in d_diff.keys():
                s_gnuplot = '%f %i %i %i %f %f %f %i %i %f\n' %(
                    rmsd,
                    int(d_parameters['mutations']),
                    int(d_parameters['chains']),
                    int(d_parameters['residues']),
                    d_diff['pH'],
                    d_diff['T'],
                    d_max['res'],
                    binary_spacegroups_identical,
                    binary_transformations,
                    (
                        0.59483772597790929*int(d_parameters['mutations'])
                        -0.53773619667629347*int(d_parameters['chains'])
                        -0.00024641340425029173*int(d_parameters['residues'])
                        -0.084648260889715471*d_diff['pH']
                        -0.0011915915999097186*d_diff['T']
                        +0.24537238074043699*d_max['res']
                        -0.62982428479768826*binary_spacegroups_identical
                        -0.046597220190475144*binary_transformations
                        )
                    )
                l_gnuplot += [s_gnuplot]


        return (
            d_data_ratio_discrete,
            d_data_ratio_continuous,
            d_data_nominal,
            d_pdbs_skip,
            d_combined,
            l_gnuplot,
            )


    def gunzip(self):

        for path in [
            self.path_pdb,
            '/data/mmCIF/',
            ]:

            subdirs = os.listdir(path)
            subdirs.sort()
            for subdir in subdirs:
                print 'gunzip', path, subdir
                files = os.listdir('%s/%s' %(path,subdir,))
                for file in files:
                    if file[-2:] == 'gz':
                        ## remove zipped file before gunzip
                        if os.path.isfile('%s/%s/%s' %(path,subdir,file[:-3])):
                            os.remove('%s/%s/%s' %(path,subdir,file[:-3]))
                        ## gunzip
                        os.system('gunzip %s/%s/%s' %(path,subdir,file))

        return


    def obsolete(self,):

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
##                    os.system('sudo rm -rf %s' %(fn))
                    os.system('sudo rm %s' %(fn))

        return

    def rsync(self):

        cmd = 'rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ %s' %(self.path_pdb)
        print cmd

        os.system(cmd)
        ## rsync -avz --port=8730 guest@rcsb-rsync-4.rutgers.edu::pdb-v3.2 .
##        os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ %s' %(self.path_pdb))

        return


    def append_ss(self,d_header,d_ATOMseq):

        for chain in d_header['SEQRES']['chains'].keys():
            if chain not in d_ATOMseq.keys():
                continue
            for i in range(1,len(d_ATOMseq[chain]['seq'])):
                res_no = d_ATOMseq[chain]['res_nos'][i]
                iCode = d_ATOMseq[chain]['iCodes'][i]
                seqID = '%4i%1s' %(res_no,iCode,)
                ## check whether helix *and* sheet (temp!!!) tmp
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


    def phipsi(
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
                'phi':{},'psi':{},'chains':[],'res_no':{},'res_name':{},
                'r':{},'ss':{},'res_name_next':{},'ss_prev':{},'ss_next':{},
                },
            pdb2:{
                'chains_equivalent':list(l_equivalent_chains[1]),
                'phi':{},'psi':{},'chains':[],'res_no':{},'res_name':{},
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
                ## get res_no,iCode if mutated chain
                ##
                res_no = d_ATOMseq[pdb][chain]['res_nos'][res_index]
                iCode = d_ATOMseq[pdb][chain]['iCodes'][res_index]
                record = d_ATOMseq[pdb][chain]['records'][res_index]
                ss = d_ATOMseq[pdb][chain]['ss'][res_index]

                ##
                ## get coordinates and calculate dihedrals
                ##
                r = 0

                ## residue missing
                if record == 'REMARK465':
                    phi = '   N/A'
                    psi = '   N/A'
                    res_name = 'N/A'

                ## atoms missing
                elif (
                    'N' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys() or
                    'CA' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys() or
                    'C' not in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys()
                    ):
                    phi = '   N/A'
                    psi = '   N/A'
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
                        if record_prev == 'REMARK465':
                            phi = '   N/A'
                        else:
                            C_prev = d_coordinates[pdb]['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['C']['coordinate']
                            phi = '%6.1f' %(self.dihedral(C_prev,N,CA,C,))
                            dist = math.sqrt(sum((C_prev-N)**2))
                            if dist > 2.:
                                print dist
                                stop_dist
                    else:
                        phi = '   N/A'
                    if res_index+1 < len(d_ATOMseq[pdb][chain]['res_nos']):
                        res_no_next = d_ATOMseq[pdb][chain]['res_nos'][res_index+1]
                        iCode_next = d_ATOMseq[pdb][chain]['iCodes'][res_index+1]
                        record_next = d_ATOMseq[pdb][chain]['records'][res_index+1]
                        if record_next == 'REMARK465':
                            psi = '   N/A'
                        else:
                            N_next = d_coordinates[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['N']['coordinate']
                            psi = '%6.1f' %(self.dihedral(N,CA,C,N_next))
                            dist = math.sqrt(sum((C-N_next)**2))
                            if dist > 2.:
                                print dist
                                stop_dist
                    else:
                        psi = '   N/A'
                    l_res_name = []
                    for altloc in d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                        res_name = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']
                        l_res_name += [res_name]
                    ## return if multiple hetID for resID
                    if len(set(l_res_name)) > 1:
                        return

                    ## sidechain dihedrals
                    if res_name in self.d_dihedrals.keys():
                        l_dihedrals = self.d_dihedrals[res_name].keys()
                        l_dihedrals.sort()
                        print res_name
                        for s_dihedral in l_dihedrals:
                            if s_dihedral[:4] != 'chi1':
                                continue
                            atom_name1 = self.d_dihedrals[res_name][s_dihedral][0]
                            atom_name2 = self.d_dihedrals[res_name][s_dihedral][1]
                            atom_name3 = self.d_dihedrals[res_name][s_dihedral][2]
                            atom_name4 = self.d_dihedrals[res_name][s_dihedral][3]
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
                            angle = self.dihedral(coord1,coord2,coord3,coord4,)
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
                            print pdb, chain, res_no, iCode, res_name, s_dihedral, phi, psi, r, angle

                ## is next residue a proline?
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
                d_loop[pdb]['ss'][chain] = ss[0]
                d_loop[pdb]['ss_prev'][chain] = ss_prev[0]
                d_loop[pdb]['ss_next'][chain] = ss_next[0]
                d_loop[pdb]['phi'][chain] = phi
                d_loop[pdb]['psi'][chain] = psi
                d_loop[pdb]['r'][chain] = r
                d_loop[pdb]['res_no'][chain] = res_no
                d_loop[pdb]['res_name'][chain] = res_name
                d_loop[pdb]['res_name_next'][chain] = res_name_next
                d_loop[pdb]['chains'] += [chain]
                print pdb,chain,'phi,psi',phi,psi

        ##
        ## read old data and prepare new data
        ##
        lines2 = []

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
            wt_accessibility = self.whatif_accessibility(pdb_wt,chain_wt,d_loop[pdb_wt]['res_no'][chain_wt],)
            mutant_accessibility = self.whatif_accessibility(pdb_mutant,chain_mutant,d_loop[pdb_mutant]['res_no'][chain_mutant],)
            ## calculate sasa for biounit for which no vicinal residues are missing
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
            line = '%4s %4s %2i %2i %1s %1s %4i %4i %3s %3s %6s %6s %6s %6s %1i %1i %4.1f %4.1f %1s %1s %3s %3s %1s %1s %1s %1s %2i %s %s\n' %(
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
                wt_accessibility,mutant_accessibility,
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
            lines2 += [line]

        ## include old lines not identical to new lines
        fd = open('single_point_mutations/%s.txt' %(pdb_wt[1],),'r')
        lines1 = fd.readlines()
        fd.close()
        lines3 = []
        for line in lines1:
            if len(line.split()) in [27,28,29,30,31,]:
                if line.split()[:4] != [pdb_wt,pdb_mutant,str(bm_wt),str(bm_mutant),]:
                    lines3 += [line]
            elif len(line.split()) in [25,26,]: ## temp!!!
                print line
                stop
                if line.split()[:2] != [pdb_wt,pdb_mutant,]:
                    line2 = line[:9]+'  1  1'+line[9:]
                    lines3 += [line2]
            else:
##                continue
                print len(line.split())
                print line
                print 'single_point_mutations/%s.txt' %(pdb_wt[1],)
                stop

        fd = open('single_point_mutations/%s.txt' %(pdb_wt[1],),'w')
        fd.writelines(lines3+lines2)
        fd.close()

        return


    def is_single_mutation(
        self,single_mutation,d_header,d_coordinates,
        pdb1,pdb2,bm1,bm2,chain1,chain2,res_no1,res_no2,iCode1,iCode2,comment1,comment2,
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
            fd = open('SEQADV_CONFLICT.txt','a')
            fd.write('%4s %4i %s _struct.title %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
            fd.close()
        elif '%1s%i%1s' %(single_mutation[-2], res_no2, single_mutation[-1],) in d_header[pdb2]['TITLE']:
            pdb_wt = pdb1
            pdb_mutant = pdb2
            bm_wt = bm1
            bm_mutant = bm2
            comment = comment2
            res_no_mutant = res_no2
            fd = open('SEQADV_CONFLICT.txt','a')
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
            fd = open('SEQADV_CONFLICT.txt','a')
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
            fd = open('SEQADV_CONFLICT.txt','a')
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
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif '%3s TO %3s' %(self.d_res3[single_mutation[-1]],self.d_res3[single_mutation[-1]],) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            stop_reverse6
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ASP52
##        elif comment2 == '' and '%3s%i' %(self.d_res3[single_mutation[-1]],res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and '%3s%i' %(self.d_res3[single_mutation[-2]],res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('SEQADV_CONFLICT.txt','a')
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
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and '%3s %i' %(self.d_res3[single_mutation[-1]],res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            stop_reverse8
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ASPARTIC ACID 52
##        elif comment2 == '' and '%s %i' %(self.d_res_names[single_mutation[-1]],res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and '%s %i' %(self.d_res_names[single_mutation[-2]],res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ' 52 '
##        elif comment2 == '' and ' %i ' %(res_no1,) in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and ' %i ' %(res_no2,) in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        ## ' MUTA*'
##        elif comment2 == '' and ' MUTA' in d_header[pdb1]['TITLE']:
##            pdb_wt = pdb2
##            pdb_mutant = pdb1
##            bm_wt = bm2
##            bm_mutant = bm1
##            comment = comment1
##            fd = open('SEQADV_CONFLICT.txt','a')
##            fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no_mutant,comment,d_header[pdb_mutant]['TITLE'],))
##            fd.close()
##        elif comment1 == '' and ' MUTA' in d_header[pdb2]['TITLE']:
##            pdb_wt = pdb1
##            pdb_mutant = pdb2
##            bm_wt = bm1
##            bm_mutant = bm2
##            comment = comment2
##            fd = open('SEQADV_CONFLICT.txt','a')
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
                fd = open('SEQADV_missing.txt','a')
                fd.write('%s %s %s %s %s %s %s --- %s --- %s\n' %(pdb1,pdb2,chain1,chain2,res_no1,res_no2,single_mutation,d_header[pdb1]['TITLE'],d_header[pdb2]['TITLE'],))
                fd.close()
                return False,None,None

        ## heavy and light chain of variant antibody
        elif chain1 == chain2 == 'H' or chain1 == chain2 == 'L':
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
                bool_mutation, s = self.parse_mmCIF(pdb,chain,res_no,res_wt,res_mutant,)
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
                    fd = open('SEQADV_CONFLICT.txt','a')
                    fd.write('%4s %4i %s %s\n' %(pdb_mutant,res_no2,comment,s,))
                    fd.close()
                elif pdb == pdb2:
                    pdb_wt = pdb1
                    pdb_mutant = pdb2
                    bm_wt = bm1
                    bm_mutant = bm2
                    comment = comment2
                    fd = open('SEQADV_CONFLICT.txt','a')
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
                stop_unknown_which_is_mutant_or_both_mutants

        return True,pdb_wt,pdb_mutant


    def parse_mmCIF(self,pdb,chain,res_no,res_wt,res_mutant,):

        bool_mutation = False
        s = ''

        d_characters = {
            "'":"'",
            '"':'"',
            '(':')',
            ';':';',
            }

        fd = open('/data/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
    ##    print '/data/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,)
        lines = fd.readlines()
        fd.close()
        j_line = 0
        for i_line in range(len(lines)):
            if i_line < j_line:
                continue
            ## new set of variables
            if lines[i_line].strip() == '#':
                k_line = j_line
                for j_line in range(i_line+1,len(lines)):
                    if j_line < k_line:
                        continue
                    if lines[j_line].strip() == 'loop_':
                        loop = True
                        ## loop over variables
                        l_variables = []
                        for k_line in range(j_line+1,len(lines)):
                            if lines[k_line][0] == '_':
                                l_variables += [lines[k_line].strip()]
                            else:
                                break
                        continue
                    ## not a loop
                    if lines[j_line][0] == '_':
                        variable = lines[j_line].split()[0]
                        if variable == '_struct.pdbx_descriptor':
                            if lines[j_line].strip() == variable:
                                value = lines[j_line+1].strip()
                            else:
                                value = lines[j_line][len(variable):].strip()
                            s = '_struct.pdbx_descriptor "%s"' %(value)
                            if (
                                '%1s%i%1s' %(res_wt,res_no,res_mutant,) in value
                                or
                                '%1s(%1s %i)%1s' %(res_wt,chain,res_no,res_mutant,) in value
                                ):
                                bool_mutation = True
                                break
                    ## line break
                    elif lines[j_line] == '':
                        stop
                        continue
                    ## #
                    elif lines[j_line].strip() == '#':
                        continue
                    ## loop
                    else:

                        variable = '_entity.pdbx_mutation'
                        if not variable in l_variables:
                            continue
                        index = l_variables.index(variable)

                        ## loop over values
                        l_values = []
                        l_characters = []
                        value = ''
                        for k_line in range(j_line,len(lines)):
                            if len(l_values) >= len(l_variables):
                                break
                            line = lines[k_line].rstrip()
                            l = lines[k_line].rstrip().split("'")
                            if (
                                '"' in lines[k_line]
                                or
                                "'" in lines[k_line]
                                or
                                "(" in lines[k_line]
                                or
                                ';' == lines[k_line][0]
                                or
                                l_characters != [] ## 1abw
                                ): ## 1wt1, 1a0a
                                i_character2 = 0
                                for i_character1 in range(len(line)):
                                    if i_character1 < i_character2:
                                        continue
                                    if line[i_character1] == ' ':
                                        continue
    ##                                if i_character1 == 0:
    ##                                    if line[i_character1] in [' ','"',"'",'(',]:
    ##                                        print line
    ##                                        stop
    ##                                    for i_character2 in range(i_character1+1,len(line)):
    ##                                        if line[i_character2] == ' ':
    ##                                            value = line[i_character1:i_character2]
    ##                                            l_values += [value]
    ##                                            break

                                    ## multi-line value
                                    if line[i_character1] == ';':

                                        if i_character1 != 0:
                                            print line
                                            stop

                                        ## semicolon end
                                        if l_characters == [';']:
                                            if line.strip() != ';':
                                                print line
                                                stop
                                            l_characters = []
                                            value += line[1:].strip()
                                            l_values += [value]
                                            value = ''
                                            break ## next line

                                        ## semicolon start
                                        elif l_characters == []:
                                            l_characters += [line[i_character1]]
                                            value += line[1:].strip()
                                            break ## next line

                                        else:
                                            print l_characters
                                            print line
                                            stop

                                    ## multi-line continuation
                                    elif l_characters == [';'] and ';' not in line:
                                        value += line.strip()
                                        break ## next line
        
                                    ## quoted value
                                    elif line[i_character1] in ['"',"'",'(',]:

    ##                                    if len(l_characters) > 0 and line[i_character1] == d_characters[line[i_character1]]:
    ##                                        l_characters = []
    ##                                        print line
    ##                                        stop
    ##                                    elif len(l_characters) > 0:
    ##                                        stop
    ##                                    else:
                                        l_characters += [line[i_character1]]

                                        if lines[k_line].strip() == ';':
                                            print lines[k_line-1]
                                            print l_characters
                                            print value
                                            stop
                                        if len(l_characters) > 1:
                                            print l_characters
                                            stop

                                        ## quoted value
                                        for i_character2 in range(i_character1+1,len(line)):
                                            if line[i_character2] == d_characters[line[i_character1]]:
                                                if len(l_characters) == 1:
                                                    i_character2 += 1
                                                    value = line[i_character1+1:i_character2-1]
                                                    if value.strip() == '' or value.strip() == ';':
                                                        print value
                                                        print line
                                                        stop
                                                    l_values += [value]
                                                    l_characters = []
                                                    value = ''
                                                    break
                                                else:
                                                    l_characters = l_characters[:-1]
                                            elif line[i_character2] == line[i_character1]:
                                                l_characters += [line[i_character2]]

                                    ## space
                                    elif line[i_character1] == ' ':
                                        continue

                                    ## non-quoted value
                                    else:
                                        if i_character1+1 == len(line):
                                            if line[i_character1] == ' ':
                                                stop
                                            value = line[i_character1]
                                            l_values += [value]
                                            value = ''
                                            if len(l_characters) > 0:
                                                print l_characters
                                                print value
                                                print line
                                                stop
                                        for i_character2 in range(i_character1+1,len(line)):
                                            if line[i_character2] == ' ':
                                                value = line[i_character1:i_character2]
                                                l_values += [value]
                                                if len(l_characters) > 0:
                                                    print l_characters
                                                    print value
                                                    print line
                                                    stop
                                                value = ''
                                                break

                            elif len(l) == 1:
                                l_values += lines[k_line].split()

                            else:
                                for s in l:
                                    if s == '': ## first character is '
                                        continue
                                    if s[0] == ' ' or s[-1] == ' ':
                                        l_values += s.split()
                                    else:
                                        l_values += [s]

                        if len(l_variables) != len(l_values):
                            print lines[j_line]
                            print lines[j_line+1]
                            print lines[j_line+2]
                            print lines[k_line]
                            print
                            print 'variables', l_variables
                            print 'values', l_values
    ##                                    print lines[j_line].split("'")
    ##                                    print lines[j_line].split()
                            print
                            for i in range(min(len(l_variables),len(l_values))):
                                print l_variables[i], l_values[i]
                            print pdb
                            stop

                        value = l_values[index]
                        s = '%s "%s"' %(variable,value,)
                        if value != '?':
                            if (
                                '%1s%i%1s' %(res_wt,res_no,res_mutant,) in value
                                or
                                '%1s(%1s %i)%1s' %(res_wt,chain,res_no,res_mutant,) in value
                                ):
                                bool_mutation = True
                                break

                        ## break j_line loop
                        if bool_mutation == True:
                            stop
                            break

                    ## break j_line loop
                    if bool_mutation == True:
                        break
            ## break i_line loop
            if bool_mutation == True:
                break

        return bool_mutation, s


    def whatif_accessibility(self,pdb,chain,res_no):

        l_acc = []
        source = '%s/whatif/%s.txt' %(self.path_cwd,pdb)
        if 1 == 1:
##        if not os.path.isfile(source):
            fd = open(source,'w')
            fd.writelines([
##                '/software/whatif/DO_WHATIF.COM <<EOF\n',
                '/local/software/whatif/DO_WHATIF.COM <<EOF\n',
                'GETMOL %s/%s/pdb%s.ent\n' %(self.path_pdb, pdb[1:3], pdb),
                '%s\n' %(pdb),
                '%DELWAT\n',
                ##'WATRAD=1.4\n',
                ##'ACCTYP=1\n',
                '%SETACC\n',
                'NOWAT 0\n',
                'NOWAT 0\n',
                '%LISTA\n',
                'TOT 0\n',
                'STOP\n',
                'Y\n',
                ])
            fd.close()
            os.system('source %s > whatif/%s.out' %(source, pdb))

            ## clean up the mess that whatif left
            for fn in ['ALTERR.LOG','PDBFILE','WHATIF.FIG','DAVADRUG.PDB',]:
                time.sleep(0.05)
                if os.path.isfile(fn):
                    try:
                        os.remove(fn)
                    except:
                        None
            l = os.listdir(os.getcwd())
            for s in l:
                if s[:3] == 'DRG' and s[6] == '.':
                    time.sleep(0.01)
                    if os.path.isfile(s):
                        try:
                            os.remove(s)
                        except:
                            None

            fd = open('whatif/%s.out' %(pdb),'r')
            lines = fd.readlines()
            fd.close()
            for i in range(len(lines)):
                line = lines[i]
                if line[:37] == 'New accessibilities calculated ... : ':
                    k = i
                    bool_break = False
                    ## loop until residue
                    for j in range(i+2,len(lines)):
                        line = lines[j]
                        if 'Residue:' in line:
                            index = line.index('Residue:')
                            res_name = lines[j][index+15:index+18]
                            if res_name == 'HOH':
                                continue
                            if lines[j][index+21:index+25].strip() == 'OXT':
                                continue
                            if not (
                                chain == lines[j][index+29]
                                and
                                res_no == int(lines[j][index+21:index+25])
                                ):
                                continue
##                            res_index = int(lines[j][index+10:index+14])
##                            res_no = int(lines[j][index+21:index+25])
                            ## loop until column headers of residue
                            for k in range(j+1,len(lines)):
                                line = lines[k]
                                if 'Atom     X     Y     Z   Acc   B   WT   VdW  Colr   AtOK  Val' in line:
                                    ## loop over data lines
                                    for l in range(k+1,len(lines)):
                                        line = lines[l]
                                        if line == ' \n':
                                            continue
                                        ## next residue
                                        if 'Residue:' in line:
                                            break
                                        ## next residue (WHATIF failed..?)
                                        if 'Option not found' in line:
                                            break
                                        ## next residue
                                        if 'Which in turn is attached to' in line:
                                            break
                                        atom_name = line[:4].strip()
                                        ## WHATIF2PDB_nomenclature
                                        if atom_name == "O'":
                                            atom_name = 'O'
                                        if (
                                            chain == lines[j][index+29]
                                            and res_no == int(lines[j][index+21:index+25])
##                                            and atom_name not in ['N','CA','C','O',]
##                                    and (
##                                        (res_name not in ['VAL','ILE','THR','SER','CYS',] and atom_name in ['CB','CG',])
##                                        or
##                                        (res_name in ['VAL','ILE',] and atom_name in ['CB','CG1',])
##                                        or
##                                        (res_name in ['CYS',] and atom_name in ['CB','SG',])
##                                        or
##                                        (res_name in ['SER',] and atom_name in ['CB','OG',])
##                                        or
##                                        (res_name in ['THR',] and atom_name in ['CB','CG2',])
##                                        )
                                            ):
                                            if line[24:28].strip() == '---':
                                                continue
##                                            print pdb
##                                            print 'i', lines[i]
##                                            print 'j', lines[j]
##                                            print 'k', lines[k]
##                                            print 'l', lines[l]
                                            acc = float(line[24:28])
                                            l_acc += [acc]
                                            if len(l_acc) > 20:
                                                stop
                                            bool_break = True
##                                            print 'acc', pdb, chain, int(lines[j][index+21:index+25]), res_name, atom_name, acc
                                if bool_break == True:
                                    break
                        if bool_break == True:
                            break
                    if bool_break == True:
                        break

        if len(l_acc) == 0: ## residue not present
            accessibility = 0
        elif len(l_acc) > 14: ## TRP
            print l_acc
            print len(l_acc)
            stop
        else:
            accessibility = sum(l_acc)/len(l_acc)

        return accessibility

    def dihedral(self,c1,c2,c3,c4):

        v1 = c2-c1
        v2 = c3-c2
        v3 = c4-c3

        angle = math.atan2(
            numpy.dot(
                math.sqrt(sum(v2*v2))*v1,
                self.cross(v2,v3),
                ),
            numpy.dot(
                self.cross(v1,v2),
                self.cross(v2,v3),
                ),
            )
        angle *= 180./math.pi

        return angle


    def cross(self,v1,v2):

        n = numpy.array([
            v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0],
            ])

        return n


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


    def whatif(self,pdb,biomolecule,d_coordinates):

        source = '%s/whatif/%s.txt' %(self.path_cwd,pdb1)
        if not os.path.isfile(source):
            fd = open(source,'w')
            fd.writelines([
                '/software/whatif/DO_WHATIF.COM <<EOF\n',
                'GETMOL %s/%s/pdb%s.ent\n' %(self.path_pdb, pdb1[1:3], pdb1),
                '%s\n' %(pdb1),
                '%DELWAT\n',
                ##'WATRAD=1.4\n',
                ##'ACCTYP=1\n',
                '%SETACC\n',
                'NOWAT 0\n',
                'NOWAT 0\n',
                '%LISTA\n',
                'TOT 0\n',
                'STOP\n',
                'Y\n',
                ])
            fd.close()
            os.system('source %s > whatif/%s_%s.out' %(source, pdb1, biomolecule1))
            fd = open('whatif/%s_%s.out' %(pdb1, biomolecule1),'r')
            lines = fd.readlines()
            fd.close()
            for i in range(len(lines)):
                line = lines[i]
                if line[:37] == 'New accessibilities calculated ... : ':
                    k = i
                    for j in range(i+2,len(lines)):
                        print j
                        if j < k+1:
                            continue
                        line = lines[j]
                        if 'Residue:' in line:
                            index = line.index('Residue:')
                            chain = line[index+29]
                            res_no = int(line[index+21:index+25])
                            for k in range(j+2,len(lines)):
                                line = lines[k]
                                if line == ' \n':
                                    break
                                atom_name = line[:4].strip()
                                acc = float(line[24:28])
                                if atom_name == "O'" and line[67] == 'O':
                                    atom_name = 'O'
                                d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][' ']['atoms'][atom_name]['acc'] = acc
                        else:
                            break
                    break

        os.remove('ALTERR.LOG')
        os.remove('PDBFILE')
        os.remove('WHATIF.FIG')

        return d_coordinates


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
        ## ignore ions and ignore saccharides
        if (
            d_hetIDs[pdb1]-set(self.d_ions.keys())-set(self.d_saccharides.keys())-set(self.d_stereoisomers.keys())-set(self.l_solutes)
            !=
            d_hetIDs[pdb2]-set(self.d_ions.keys())-set(self.d_saccharides.keys())-set(self.d_stereoisomers.keys())-set(self.l_solutes)
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
            return True, d_hetIDS

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

        ## continue if superseded structure
        if 'SPRSDE' in d_header[pdb].keys():
            pdbskip = True
            stop_delete_this_if_statement_because_no_cases
            return pdbskip

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


    def rmsd2bfactor(
        self, pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_coordinates, d_header,
        tv1, rm, tv2,
        l_equivalent_chains, bmchains1, bmchains2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_ATOMseq,
        ):

        print 'rmsd2bfactor'

        biomolecule1 = str(biomolecule1).zfill(2)
        biomolecule2 = str(biomolecule2).zfill(2)

        ##
        ## append secondary structure information
        ##
        d_lines_ss = {
            pdb1:[],
            pdb2:[],
            }
        for pdb in d_lines_ss:
            ##
            ## read lines
            ##
            fd = open('%s/%s/pdb%s.ent' %(self.path_pdb, pdb.lower()[1:3], pdb.lower()),'r')
##            fd = open('C:\Users\Tommy Carstensen\Downloads\pdb\%s.pdb' %(pdb.lower(),),'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                record = line[:6].strip()
                if record in ['HELIX','SHEET','TURN',]:
                    d_lines_ss[pdb] += line


        chains1 = l_equivalent_chains[0]
        chains2 = l_equivalent_chains[1]
        if len(chains1) != len(chains2):
            print pdb1, pdb2
            print chains1, chains2
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
            rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2,
            )

        ##
        ## append lines by model to final lines
        ##
        pdblines1 = d_lines_ss[pdb1]
        pdblines2 = d_lines_ss[pdb2]
        for model in d_lines[pdb1].keys():
            pdblines1 += ['MODEL     %4s                                                                  \n' %(model)]
            pdblines1 += d_lines[pdb1][model]
            pdblines1 += ['ENDMDL                                                                          \n']
        for model in d_lines[pdb2].keys():
            pdblines2 += ['MODEL     %4s                                                                  \n' %(model)]
            pdblines2 += d_lines[pdb2][model]
            pdblines2 += ['ENDMDL                                                                          \n']

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
            print 'gif thumbnail', pdb
            biomolecule = d_biomolecules[pdb]
            if pdb == pdb1:
                prefix = pdb1+str(biomolecule1).zfill(2)+pdb2+str(biomolecule2).zfill(2)
            if pdb == pdb2:
                prefix = pdb2+str(biomolecule2).zfill(2)+pdb1+str(biomolecule1).zfill(2)
            ## write rasmol script
            lines = [
                'rasmol -nodisplay pdb/%s/%s.pdb <<EOF\n' %(pdb[1], prefix),
                'color temperature\n',
                'wireframe 0\n',
                'cartoon\n',
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

        
    def write_html(self, d_rmsd, d_header, prefix):

        ## keys=htmlkeys, values=htmlcolumns (table headings)
        d_columns_headers = {
            'REMARK465':'missing<br>residues',
            'REMARK470':'missing<br>atoms',
            'transformations':'transformations',
            'identical authors':'identical<br>authors',
            'gif1':'gif1','gif2':'gif2',
            'pdb1':'pdb1', 'pdb2':'pdb2',
            'bm1':'bm1', 'bm2':'bm2',
            'rmsd':'<a href="http://en.wikipedia.org/wiki/Protein_structural_alignment">rmsd</a>',
            'chains':'chains (biounit)', 'residues':'residues (biounit)', 'coordinates':'coords (biounit)',
            'pH1':'pH1', 'pH2':'pH2', 'T1':'T1', 'T2':'T2',
            'res1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">resolution1</a>',
            'res2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">resolution2</a>',
            'spacegroup1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">space1</a>',
            'spacegroup2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">space2</a>',
            'title1':'title1','title2':'title2',
            'hetIDs1':'hetIDs1', 'hetIDs2':'hetIDs2',
            'mutations':'mutations (biounit)',
            }

        ## initiate html lines
        th = '<tr>\n'
        for column in self.l_columns_html:
            th += '<td>%s</td>\n' %(d_columns_headers[column])
        th += '</tr>\n'

        l_tr = ''

        d_html = {}

        ##
        ## loop over pdbs and biomolecules
        ##
        for pdb1 in d_rmsd:

            ## parse physiochemical properties (pH and temperature)
            if 'REMARK200' in d_header[pdb1].keys():
                T1 = d_header[pdb1]['REMARK200']['TEMPERATURE']
                pH1 = d_header[pdb1]['REMARK200']['PH']
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
            spacegroup1 = d_header[pdb1]['CRYST1'].rjust(10)
            ## parse hetIDs
            hetIDs1 = set()
            ## parse resolution
            res1 = d_header[pdb1]['REMARK2']
            try:
                res1 = '%5.2f' %(float(res1))
            except:
                res1 = '%s' %(res1.rjust(5))
            ## parse hetIDs
            for chain in d_header[pdb1]['HET'].keys():
                for res_no in d_header[pdb1]['HET'][chain].keys():
                    for iCode in d_header[pdb1]['HET'][chain][res_no].keys():
                        hetIDs1 |= d_header[pdb1]['HET'][chain][res_no][iCode]
            hetIDs1 = list(hetIDs1)

            for bm1 in d_rmsd[pdb1].keys():

                for pdb2 in d_rmsd[pdb1][bm1]:

                    ## parse physiochemical properties (pH and temperature)
                    if 'REMARK200' in d_header[pdb2].keys():
                        T2 = d_header[pdb2]['REMARK200']['TEMPERATURE']
                        pH2 = d_header[pdb2]['REMARK200']['PH']
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
                    spacegroup2 = d_header[pdb2]['CRYST1'].rjust(10)
                    ## parse hetIDs
                    hetIDs2 = set()
                    ## parse resolution
                    res2 = d_header[pdb2]['REMARK2']
                    try:
                        res2 = '%5.2f' %(float(res2))
                    except:
                        res2 = '%s' %(res2.rjust(5))
                    ## parse hetIDs
                    for chain in d_header[pdb2]['HET'].keys():
                        for res_no in d_header[pdb2]['HET'][chain].keys():
                            for iCode in d_header[pdb2]['HET'][chain][res_no].keys():
                                hetIDs2 |= d_header[pdb2]['HET'][chain][res_no][iCode]
                    hetIDs2 = list(hetIDs2)

                    for bm2 in d_rmsd[pdb1][bm1][pdb2].keys():

                        ## parse rmsd and related data
                        rmsd = '%5.2f' %(d_rmsd[pdb1][bm1][pdb2][bm2]['rmsd'])
                        n_chains = '%3i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['chains'])
                        n_residues = '%4i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['residues'])
                        n_coordinates = '%5i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['coordinates'])
                        mutations = '%2i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['mutations'])
                        if 'REMARK465' in d_header[pdb1].keys() or 'REMARK465' in d_header[pdb2].keys():
                            remark465 = True
                        else:
                            remark465 = False
                        if 'REMARK470' in d_header[pdb1].keys() or 'REMARK470' in d_header[pdb2].keys():
                            remark470 = True
                        else:
                            remark470 = False
                        transformations = d_rmsd[pdb1][bm1][pdb2][bm2]['transformations']
                        bool_identical_authors = d_rmsd[pdb1][bm1][pdb2][bm2]['bool_identical_authors']

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
                            'identical authors':'<td>%s' %(bool_identical_authors),
                            'REMARK465':'<td>%s' %(remark465),
                            'REMARK470':'<td>%s' %(remark470),
                            'title1':'<td style="font-size:70%%">%s' %(d_header[pdb1]['TITLE']),
                            'title2':'<td style="font-size:70%%">%s' %(d_header[pdb2]['TITLE']),
                            'gif1':'<td><a href="../pdb/%s/%s.pdb"><img src="../gif/%s/%s.gif"></a>' %(pdb1[1],prefix1,pdb1[1],prefix1),
                            'gif2':'<td><a href="../pdb/%s/%s.pdb"><img src="../gif/%s/%s.gif"></a>' %(pdb2[1],prefix2,pdb2[1],prefix2),
                            'pdb1':'<td><a href="../htm/%s.htm">%s</a>' %(pdb1,pdb1),
                            'pdb2':'<td><a href="../htm/%s.htm">%s</a>' %(pdb2,pdb2),
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
                            'bool_identical_authors':bool_identical_authors,
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
            path = 'tmp/'
            file = '%s.htm' %(prefix)
            self.append_table_rows(path,file,l_tr,th,d_rmsd)

        return


    def append_table_rows(self,path, file,l_tr,th,d_rmsd):

        if os.path.isfile('%s%s' %(path,file)):
            fd = open('%s%s' %(path, file),'r')
            lines_prev1 = fd.readlines()
            fd.close()
            for i in range(len(lines_prev1)-1,-1,-1,):
                if lines_prev1[i] in ['<table border="1" class="sortable" id="sortable">\n','<table border="1">\n',]:
                    lines_prev1 = lines_prev1[i+1:-1]
                    break
                
            n_lines = 1+len(self.l_columns_html)+1
            lines_prev2 = []

            for i_row in range(1,len(lines_prev1)/(n_lines),1,):
                ## parse pdbs
                pdb_lines = lines_prev1[i_row*n_lines+3:i_row*n_lines+5]
                bm_lines = lines_prev1[i_row*n_lines+5:i_row*n_lines+7]
                l_pdbs = []
                l_bms = []
                for j in range(2):
                    pdb_line = pdb_lines[j]
                    bm_line = bm_lines[j]
                    try:
                        pdb_index2 = pdb_line.index('</a>')
                    except:
                        print lines_prev1[i_row*n_lines:(i_row+1)*n_lines]
                        print lines_prev1[0*n_lines:1*n_lines]
                        print lines_prev1[0]
                        print lines_prev1[-1]
                        print path,file, pdb_line
                        print len(lines_prev1), n_lines
                        stop
                    pdb_index1 = pdb_line[:pdb_index2].rindex('>')+1
                    pdb = pdb_line[pdb_index1:pdb_index2]
                    bm_index2 = bm_line.index('</td>')
                    bm_index1 = bm_line[:bm_index2].rindex('>')+1
                    bm = bm_line[bm_index1:bm_index2]
                    l_pdbs += [pdb]
                    l_bms += [bm]
                pdb1 = l_pdbs[0]
                pdb2 = l_pdbs[1]
                try:
                    bm1 = int(l_bms[0])
                except:
                    print lines_prev1[i_row*n_lines:(i_row+1)*n_lines]
                    print l_bms
                    print bm_lines
                    print file
                    stop
                bm2 = int(l_bms[1])
                current = False
                if pdb1 in d_rmsd.keys():
                    if bm1 in d_rmsd[pdb1].keys():
                        if pdb2 in d_rmsd[pdb1][bm1].keys():
                            if bm2 in d_rmsd[pdb1][bm1][pdb2].keys():
                                current = True
                if pdb2 in d_rmsd.keys():
                    if bm2 in d_rmsd[pdb2].keys():
                        if pdb1 in d_rmsd[pdb2][bm2].keys():
                            if bm1 in d_rmsd[pdb2][bm2][pdb1].keys():
                                current = True
                if current == False:
                    lines_prev2 += lines_prev1[i_row*n_lines:(i_row+1)*n_lines]

            html = []
            ## javascript
            html += ['<script type="text/javascript" src="../sortable.js"></script>\n']
            ## table init
            html += ['<table border="1" class="sortable" id="sortable">\n']
            ## table header
            html += lines_prev1[:n_lines]
            ## previous lines
            html += lines_prev2
            ## new/updated lines
            html += [l_tr]
            ## table term
            html += ['</table>\n']
        else:
            html = ['<table border="1">\n'+th+l_tr+'</table>\n']
        fd = open('%s%s'%(path,file),'w')
        fd.writelines(html)
        fd.close()

        return


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
        l_atoms = [],
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
        l_atoms=[],
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
                if atom_name == 'CA':
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


    def determine_if_modres(self, d_header, d_coordinates, chain, res_no, iCode, res_name):

        modres = False
        if chain in d_header['MODRES'].keys():
            if res_no in d_header['MODRES'][chain].keys():
                if iCode in d_header['MODRES'][chain][res_no].keys():
                    l_hetIDs = []
                    for altloc in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                        l_hetIDs += [d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']]
                    if len(set(l_hetIDs) & d_header['MODRES'][chain][res_no][iCode]) > 0:
                        modres = True
                    else:
                        print chain, res_no, iCode, res_name
                        print d_header['MODRES'][chain][res_no][iCode]
                        print l_hetIDs
                        print 
                        print self.cluster, 'cluster'
                        notexpected

        return modres


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


    def parse_header(self, s_pdb, stop_error = True):

        print 'parsing header of', s_pdb

        ##
        ## read lines
        ##
        fd = open('%s/%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
##        fd = open('C:\Users\Tommy Carstensen\Downloads\pdb\%s.pdb' %(s_pdb.lower(),),'r')
        lines = fd.readlines()
        fd.close()

        ##
        ## set dictionaries
        ##
        d_header = {}
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
##                    if 'REMARK350' in d_header.keys() and 'HET' in d_header.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = self.parse_atom_no_range(d_conect, 'atom', atom_no)

            elif record == 'HETATM': ## section 9
                continue
##                if parse_atoms == False:
##                    if 'REMARK350' in d_header.keys() and 'HET' in d_header.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = self.parse_atom_no_range(d_conect, 'atom', atom_no)

            elif record == 'COMPND': ## section xxx
                d_header = self.parse_recordCOMPND(d_header, line, i, lines)

            elif record == 'SOURCE':

                if not 'SOURCE' in d_header.keys():
                    d_header['SOURCE'] = {}

                if lines[i][10:18].strip() == 'MOL_ID:':
                    mol_id = int(lines[i][18:].replace(';',''))
                else:
                    continue
                
                for j in range(i+1,len(lines)):
                    
                    if lines[j][:6] != 'SOURCE':
                        break
                    if lines[j][10:18].strip() == 'MOL_ID:':
                        break

                    if 'ORGANISM SCIENTIFIC' in lines[j]:
                        print lines[j]
                        stop
                    if 'ORGANISM_SCIENTIFIC' in lines[j]:
                        if lines[j][10:32] != ' ORGANISM_SCIENTIFIC: ':
                            print lines[j]
                            stop
                        organism = lines[j][32:80].strip()[:-1]
                        if ';' in organism:
                            print lines[j]
                            stop
                        if mol_id in d_header['SOURCE'].keys():
                            print d_header['SOURCE']
                            print lines[j]
                            stop
                        d_header['SOURCE'][mol_id] = organism

            elif record == 'MDLTYP':
                if line[10:].strip() == 'MINIMIZED AVERAGE':
                    continue
                if (
                    line[9] not in ['2','3',]
                    ## protein
                    and line[:31] != 'MDLTYP    CA ATOMS ONLY, CHAIN '
                    ## rna/dna
                    and line[:30] != 'MDLTYP    P ATOMS ONLY, CHAIN '
                    ## NMR
                    and line[:50] != 'MDLTYP    MINIMIZED AVERAGE; CA ATOMS ONLY, CHAIN '
                    ):
                    print line
                    print line[:31]
                    stop
##                chain = line[31]
##                d_header['MDLTYP'] = [chain]
                d_header['MDLTYP'] = True

            elif record == 'SSBOND':
                chain1 = line[15]
                res_no1 = int(line[17:21])
                iCode1 = line[21]
                chain2 = line[29]
                res_no2 = int(line[31:35])
                iCode2 = line[35]
                if line[59:65].strip() != '1555':
                    print line
                    stop
                if line[66:72].strip() != '1555':
                    if not 'SSBOND' in d_header.keys():
                        d_header['SSBOND'] = []
                    d_header['SSBOND'] += [[chain1,res_no1,iCode1,chain2,res_no2,iCode2,]]

            elif record == 'REMARK': ## section 2
                d_header = self.parse_recordREMARK(d_header, line, i, lines)

            elif record == 'SEQADV': ## section 3
                d_header = self.parse_recordSEQADV(line, d_header)

            elif record == 'CAVEAT': ## Title Section
                if line[6:10].strip() == '' and 'CHIRAL' not in line and 'GEOMETR' not in line and 'STEREOCHEMISTRY' not in line:
                    if s_pdb not in [
                        '1r9m',
                        '3ckz', ## PEPTIDE BACKBONE LINKAGE PROBLEM
                        '3cl0', ## PEPTIDE BACKBONE LINKAGE PROBLEM
                        '1g4b', ## NUMEROUS CLOSE CONTACTS. SEE REMARK 3
                        '3eyo', ## CLOSE CONTACTS BETWEEN SIDE CHAINS
                        '1heg', ## ABNORMALLY CLOSE CONTACTS PRESENT
                        '3dwg', ## ABNORMALLY SHORT LINK BETWEEN LYS A51 AND PLP A401
                        '3dwi', ## ABNORMALLY SHORT LINK BETWEEN LYS A51 AND PLP A401
                        '3fg4', ## C-N PEPTIDE BOND LENGTH ERROR D414-D418
                        '1fx7', ## continuation across two lines
                        '1euq', ## continuation across two lines
                        '3dwf', ## continuation across three lines
                        ]:
                        print line
                        stop

            elif record == 'SPLIT': ## Title Section
                d_header['SPLIT'] = True

            elif record == 'SEQRES': ## section 3
                d_header = self.parse_recordSEQRES(line, d_header)

            elif record == 'HET': ## section 4
                hetID = line[7:10].strip()
                ## continue if water
                if hetID in ['HOH','H2O','DOD','D2O']: ## D2O in 2JAJ
                    continue
                chain = line[12]
                res_no = int(line[13:17])
                iCode = line[17]
                if 'HET' not in d_header.keys():
                    d_header['HET'] = {}
                if chain not in d_header['HET'].keys():
                    d_header['HET'][chain] = {}
                if res_no not in d_header['HET'][chain].keys():
                    d_header['HET'][chain][res_no] = {}
                if iCode not in d_header['HET'][chain][res_no].keys():
                    d_header['HET'][chain][res_no][iCode] = set() ## multiple hetIDs possible if altlocs
                d_header['HET'][chain][res_no][iCode] |= set([hetID])

            elif record == 'MODRES': ## section 3
                hetID = line[12:15].strip()
                chain = line[16]
                res_no = int(line[18:22])
                iCode = line[22]
                res_name = line[24:27].strip()
                txt = line[29:80].strip()
                if hetID in set(self.d_res1.keys())-set(['MSE']) and res_name in set(self.d_res1.keys())-set(['MSE']):
                    continue
##                if txt not in [] and hetID in set(self.d_res1.keys()+self.l_nucleotides)-set(['MSE']):
##                    print line
##                    print s_pdb
##                    stop_duplicate_modres_records_or_std_res_name
                if 'MODRES' not in d_header.keys():
                    d_header['MODRES'] = {}
                if chain not in d_header['MODRES'].keys():
                    d_header['MODRES'][chain] = {}
                if res_no not in d_header['MODRES'][chain].keys():
                    d_header['MODRES'][chain][res_no] = {}
                if iCode not in d_header['MODRES'][chain][res_no].keys():
                    d_header['MODRES'][chain][res_no][iCode] = set() ## multiple hetIDs possible if altlocs
                d_header['MODRES'][chain][res_no][iCode] |= set([hetID])

            elif record == 'TITLE': ## section 2
                if not 'TITLE' in d_header.keys():
                    d_header['TITLE'] = line[10:].strip()
                else:
                    if d_header['TITLE'][-1] == '-':
                        d_header['TITLE'] += line[10:].strip()
                    else:
                        d_header['TITLE'] += ' '+line[10:].strip()

            elif record == 'HELIX': ## section 5
                d_header = parse_pdb.parse_recordHELIX(line,d_header,)

            elif record == 'SHEET': ## section 5
                d_header = parse_pdb.parse_recordSHEET(line,d_header,)

            elif record == 'HEADER':
                d_header['HEADER'] = line[10:50].strip()

            elif record == 'SPRSDE': ## section 2
                sIDcodes = line[31:70].lower().split()
                for sIDcode in sIDcodes:
                    if sIDcode not in d_header.keys():
                        d_header[sIDcode] = {}
                    if os.path.isfile('htm/%s' %(sIDcode)):
                        os.remove('htm/%s' %(sIDcode))
                        print sIDcode
                        stop_test_if_works
                    d_header[sIDcode]['SPRSDE'] = True

            elif record == 'EXPDTA': ## section 2
                methods = line[10:].strip().split(',')[0]
                if methods[:3] == 'NMR':
                    methods = 'NMR'
                elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                    methods = 'X-RAY'
                d_header['EXPDTA'] = methods

            elif record == 'AUTHOR':
                l_authors = line[10:].strip().split(',')
                if not 'AUTHOR' in d_header.keys():
                    d_header['AUTHOR'] = []
                d_header['AUTHOR'] += l_authors

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                if spacegroup in self.d_crystalsynonyms.keys():
                    spacegroup = self.d_crystalsynonyms[spacegroup]
                d_header['CRYST1'] = spacegroup

            else:
                continue

        ##
        ## missing records
        ##
        d_dics = {
            'chains':{},
            'HET':{},
            'MODRES':{},
            'REMARK200':{'TEMPERATURE':'N/A','PH':'N/A'},
            'REMARK525':[],
            'TITLE':'N/A',
            'EXPDTA':'N/A',
            'REMARK2':'N/A',
            }
        for key in d_dics.keys():
            if key not in d_header.keys():
                d_header[key] = d_dics[key]

        proteinchains = []
        peptidechains = []
        nucleotidechains = []
        saccharidechains = []
        for chain in d_header['SEQRES']['chains'].keys():
            if d_header['SEQRES']['chains'][chain]['type'] == 'peptide':
                peptidechains += chain
                if len(d_header['SEQRES']['chains'][chain]['seq']) > self.min_len_chain:
                    proteinchains += chain
            elif d_header['SEQRES']['chains'][chain]['type'] == 'nucleotide':
                nucleotidechains += chain
            elif d_header['SEQRES']['chains'][chain]['type'] == 'saccharide':
                saccharidechains += chain

        d_header['proteinchains'] = proteinchains

##        if s_pdb not in ['1ady','1bhj']:
##            chains = d_header['SEQRES']['chains'].keys()
##            if 'REMARK350' not in d_header.keys() and len(d_header['SEQRES']['chains'].keys()) > 1 and d_header['EXPDTA'] != 'NMR' and len(saccharidechains) != len(d_header['SEQRES']['chains'].keys()):
##                ## biounit not specified as text
##                if biounit == 'N/A':
##                    fd = open('unknownbiounit.txt','a')
##                    fd.write('%s %s %s %s %s\n' %(s_pdb, len(proteinchains), len(peptidechains), len(nucleotidechains), len(saccharidechains)))
##                    fd.close()
##                ## unequal number of proteinchains and size of biounit
##                elif self.d_biounits[biounit] != len(chains):
##                    if biounit == 'MONOMER':
##                        d_header['REMARK350'] = {}
##                        for i in range(len(proteinchains)):
##                            chain = proteinchains[i]
##                            d_header['REMARK350'][i+1] = {'chains': [chain]}
##                    elif len(chains) % self.d_biounits[biounit] == 0:
#### add all combinations of chains to remark350 transformations ?! redudant hits if done for twin pdb as well...
##                        print biounit, self.d_biounits[biounit]
##                        print s_pdb
##                        print d_header.keys()
##                        print d_header['SEQRES']['chains'].keys()
##                        print proteinchains
##                        print d_header['HET']
##                        stop3
##                    else:
##                        print s_pdb, biounit, chains, proteinchains
##                        print d_header['HET']
##                        stop3b
##                elif self.d_biounits[biounit] == len(chains):
##                    d_header['REMARK350'] = {1:{'chains': [chains]}}
##                else:
##                    print s_pdb, biounit, proteinchains, chains
##                    stop5
####                d_header['biounit'] = 'small' ## small opposed to monomer if multiple different chains
#### count number of similar chains by comparing SEQRESseqs

        return d_header


    def parse_recordCOMPND(self, d_header, line, i, lines,):

        if not 'COMPND' in d_header.keys():
            d_header['COMPND'] = {}
        if lines[i][10:18].strip() == 'MOL_ID:':
            mol_id = int(lines[i][18:].strip()[:-1])
        else:
            return d_header
        
        for j in range(i+1,len(lines)):
            if lines[j][:6] != 'COMPND':
                break
            if lines[j][10:18].strip() == 'MOL_ID:':
                break
            if lines[j][10:17].strip() == 'CHAIN:':
                chains = lines[j][17:80].strip().replace(':','').replace(';','').replace(' ','').split(',')
                for chain in chains:
                    if chain in d_header['COMPND'].keys():
                        print chain,d_header['COMPND']
                        stop
                    d_header['COMPND'][chain] = mol_id
##            if 'MUTATION: YES' in lines[j]:
##                mutation = True
##                d_header['COMPND'][chain] = True
##            elif 'MUTANT;' in lines[j]: ## e.g. 2jc2
##                mutation = True
##                d_header['COMPND'][chain] = True
##            elif 'MUTANT;' in lines[j]: ## e.g. 2jc2
##                mutation = True
##                d_header['COMPND'][chain] = True
##            elif '1 MUTATION' in lines[j]: ## e.g. 2vnd
##                mutation = True
##                d_header['COMPND'][chain] = True
##            elif 'ENGINEERED: YES' in lines[j]:
##                continue
####                mutation = True
####                d_header['COMPND'][chain] = True
##            elif 'MUTA' in lines[j] and 'MUTASE' not in lines[j]:
##                print lines[j]
##                stop_unexpected_line
##            elif 'CHAIN: ' in lines[j]:
##                chain = lines[j].split()[-1][0]
##            else:
##                continue

        return d_header


    def parse_recordSEQADV(self, line, d_header,):

        comment = line[49:70].strip()

        ## deletion
        if comment in ['DELETION','ENGINEERED DELETION',] and line[12:15].strip() == '' and line[18:23].strip() == '':
            return d_header
        ## incorrect seq in seq db?
        if comment == 'SEE REMARK 999' and line[12:15].strip() == '' and line[18:23].strip() == '':
            return d_header
        ## residue not in ATOM records (no residue number)
        if line[18:22].strip() == '':
            return d_header

        ## chromophore
        if line[12:15].strip() in [
            ## chromophore in name
            'CRO','CR2','CSY','MDO','XYG','DYG','NYG','CR8','RC7','CH6','NRQ',
            'CFY','CWR','CQR','CRQ','AYG','IIC','CRX','C99','CLV','C12','XXY',
            'CR0','GYS','X9Q','CRW','GYC','CR7','QLG','CH7',
            ## chromophore not in name
            'IEY','CRK','5ZA','CSH','CRG','4F3','MFC','PIA','CR5','CRF','NYC',
            'CZO','CCY','CJO',
            ]: ## or 'CHROMOPHORE' in comment
            return d_header

        chain = line[16]
        seq_no = int(line[18:22])
        iCode = line[22]
        res_name_ref = line[39:42].strip()
        if not 'SEQADV' in d_header.keys():
            d_header['SEQADV'] = {}
        if not chain in d_header['SEQADV'].keys():
            d_header['SEQADV'][chain] = {}
        if line[12:15].strip() not in [
            'SUI', ## ASP,GLY
            '175', ## ALA, SER, GLY
            ]:
            if line[7:11] not in [
                ## altlocs / microheterogeneity
                '2CI1','2JGE','2JGJ','3DL7',
                ## two residues
                '2E3V',
                ]:
                if '%4i%1s' %(seq_no,iCode,) in d_header['SEQADV'][chain].keys():
                    print seq_no, iCode
                    print d_header['SEQADV'][chain].keys()
                    print d_header['SEQADV'][chain]
                    print line
                    stop_duplicate_entry_maybe_chromophore
        d_header['SEQADV'][chain]['%4i%1s' %(seq_no,iCode,)] = {'res_name_seqdb':res_name_ref,'comment':comment,}

        return d_header


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


    def parse_recordSEQRES(self, line, d_header):

        chain = line[11]

        if 'SEQRES' not in d_header:
            d_header['SEQRES'] = {'chains':{},}
        if chain not in d_header['SEQRES']['chains'].keys():
            d_header['SEQRES']['chains'][chain] = {}
        if not 'type' in d_header['SEQRES']['chains'][chain].keys():
            d_header['SEQRES']['chains'][chain]['type'] = 'N/A'

        l_residues = line[19:70].split()

        s_residues = ''
        for i in range(len(l_residues)):
            residue = l_residues[i]
            if residue in self.d_res1.keys():
                if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                    d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
                elif d_header['SEQRES']['chains'][chain]['type'] != 'peptide': ## e.g. 1vq6
                    d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
                s_residues += self.d_res1[residue]
            elif residue in self.l_nucleotides:
                if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                    d_header['SEQRES']['chains'][chain]['type'] = 'nucleotide'
                elif d_header['SEQRES']['chains'][chain]['type'] != 'nucleotide':
                    stop
                s_residues += residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                    d_header['SEQRES']['chains'][chain]['type'] = 'saccharide'
                elif d_header['SEQRES']['chains'][chain]['type'] != 'saccharide':
                    stop
                s_residues += residue
            else:
                if residue == 'UNK': ## e.g. 1pny.pdb
                    if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                        d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
                s_residues += 'X'

        if 'seq' not in d_header['SEQRES']['chains'][chain].keys():
            d_header['SEQRES']['chains'][chain]['seq'] = ''
        d_header['SEQRES']['chains'][chain]['seq'] += s_residues

        if 'seq3' not in d_header['SEQRES']['chains'][chain].keys():
            d_header['SEQRES']['chains'][chain]['seq3'] = []
        d_header['SEQRES']['chains'][chain]['seq3'] += l_residues

        return d_header


    def parse_coordinates(
        self, s_pdb, d_header, verbose = False, parse_molecules = True,
        ):

        print 'parsing coordinates of', s_pdb

        ## deep copy header, because REMARK465 records are deleted during loop over biomolecules
        d_header = copy.deepcopy(d_header)

        ##
        ## read lines
        ##
        fd = open('%s/%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
##        fd = open('C:\Users\Tommy Carstensen\Downloads\pdb\%s.pdb' %(s_pdb.lower(),),'r')
        lines = fd.readlines()
        fd.close()

        ##
        ## set dictionaries
        ##
        d_atomnos = {}
        d_CONECT = {}
        d_coordinates = {}
        d_ATOMseq = {}
        model = 1

        ##
        ## loop over lines
        ##
        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'ATOM':
                if model != 1:
                    continue
                d_coordinates, d_line = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
                d_ATOMseq = self.build_ATOMseq(line, lines, i, d_ATOMseq, d_header, d_coordinates,)
                d_atomnos[d_line['atom_no']] = d_line

            elif record == 'HETATM':
                if model != 1:
                    continue
                res_name = line[17:20].strip()
                if res_name == 'UNL':
                    stop_write_code
                ## water
                if res_name in ['D2O','H2O',]:
                    print s_pdb, res_name
                    stop
                ## water must be parsed in case there is a connection to it (e.g. 2bpb)
                if res_name in ['HOH','DOD',]: ## DOD in 2d4j
                    d_coordinates, d_line, = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
                ## MSE
                elif res_name in self.d_modres.keys():
                    d_coordinates, d_line, = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
                ## (poly)saccharide or other hetero compound
                else:
                    d_coordinates, d_line, = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
                atom_no = d_line['atom_no']
                d_atomnos[atom_no] = d_line

                d_ATOMseq = self.build_ATOMseq(line, lines, i, d_ATOMseq, d_header, d_coordinates,) ## 1omw

            elif record == 'CONECT':
                d_CONECT = self.parse_recordCONECT(line, d_atomnos, d_CONECT)
                
            elif record == 'HET':
                hetID = line[7:10].strip()
                if 'HET' not in d_coordinates.keys():
                    d_coordinates['HET'] = set()
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

            elif record == 'ENDMDL': ## e.g. 2q4t
                break


        if parse_molecules == True:
            d_molecules = self.build_dictionary_of_molecules(d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=verbose,)
        else:
            d_molecules = {}

        for chain in d_header['SEQRES']['chains'].keys():

            ## skip if not peptide chain
            type = d_header['SEQRES']['chains'][chain]['type']
            if type != 'peptide':
                continue

            ## skip if short chain
            if len(d_header['SEQRES']['chains'][chain]['seq3']) < 5:
                continue

            ## skip if all residues are unknown
            if len(d_header['SEQRES']['chains'][chain]['seq3'])*['UNK'] == d_header['SEQRES']['chains'][chain]['seq3']:
                continue

            ## append remaining REMARK465 residues
            if 'REMARK465' in d_header.keys():
                if chain in d_header['REMARK465']['chains'].keys():
                    if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                        l_REMARK465_seq = list(d_header['REMARK465']['chains'][chain]['seq'])
                        d_ATOMseq,d_header = self.append_ATOMseq(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_seq,True,False,)

            if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3']:
                print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3']
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print chain
                print 'SEQRES', len(d_header['SEQRES']['chains'][chain]['seq3'])
                print 'ATOM  ', len(d_ATOMseq[chain]['seq'])
                for i in range(len(d_ATOMseq[chain]['seq'])):
                    if d_header['SEQRES']['chains'][chain]['seq3'][:i] != d_ATOMseq[chain]['seq'][:i]:
                        print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][:i]
                        print 'ATOM', d_ATOMseq[chain]['seq'][:i]
                        break
                print s_pdb,chain
                stop_different_SEQRESseq_ATOMseq

        return d_coordinates, d_molecules, d_ATOMseq


    def append_ATOMseq(
        self,record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_seq,append_REMARK465,append_ATOM
        ):

        if altloc == None:
            altloc = ' '

        if append_REMARK465 == True:
            for i_seq_REMARK465 in range(len(l_REMARK465_seq)):
                seq_REMARK465 = l_REMARK465_seq[i_seq_REMARK465]
                res_no_REMARK465 = int(seq_REMARK465[:4])
                iCode_REMARK465 = seq_REMARK465[4]
                altloc_REMARK465 = ' '
                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][altloc_REMARK465]['res_name']
                d_ATOMseq[chain]['seq'] += [res_name_REMARK465]
                d_ATOMseq[chain]['res_nos'] += [res_no_REMARK465]
                d_ATOMseq[chain]['iCodes'] += [iCode_REMARK465]
                d_ATOMseq[chain]['altlocs'] += [' ']
                d_ATOMseq[chain]['records'] += ['REMARK465']
                d_ATOMseq[chain]['indexes'] += ['%1s%4s%1s' %(chain,str(res_no).zfill(4),iCode)]
                d_ATOMseq[chain]['ss'] += ['']
                ## remove res_no for each iCode
                d_header['REMARK465']['chains'][chain]['seq'].remove(seq_REMARK465)
            for seq_REMARK465 in l_REMARK465_seq:
                res_no_REMARK465 = int(seq_REMARK465[:4])
                if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                    del d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]

        if append_ATOM == True:
            d_ATOMseq[chain]['seq'] += [res_name_ATOM]
            d_ATOMseq[chain]['res_nos'] += [res_no]
            d_ATOMseq[chain]['iCodes'] += [iCode]
            d_ATOMseq[chain]['altlocs'] += [altloc]
            d_ATOMseq[chain]['records'] += [record]
            d_ATOMseq[chain]['indexes'] += ['%1s%4s%1s' %(chain,str(res_no).zfill(4),iCode)]
            d_ATOMseq[chain]['ss'] += ['']

        return d_ATOMseq, d_header


    def build_ATOMseq(self,line,lines,i,d_ATOMseq,d_header, d_coordinates,):

        record = line[:6].strip()
        altloc = line[16]
        res_name_ATOM = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]

        skip = False

        ## water
        if res_name_ATOM == 'HOH':
            skip = True

        ## not water
        else:
            ## modified residue?
            SEQRESres = self.check_if_SEQRESres(res_name_ATOM,record,d_header,chain,res_no,iCode,)
            if SEQRESres == False:
                skip = True

            ## peptide chain?
            if chain in d_header['SEQRES']['chains'].keys():
                type = d_header['SEQRES']['chains'][chain]['type']
                if type != 'peptide':
                    skip = True
            else:
                skip = True

        if skip == True:
            return d_ATOMseq

        ## initiate new chain in dic
        if not chain in d_ATOMseq.keys():
            d_ATOMseq[chain] = {
                'seq':[],'iCodes':[],'res_nos':[],'altlocs':[],'records':[],
                'indexes':[],'ss':[],
                }

        if lines[i-1][:6].strip() in ['ATOM','HETATM','ANISOU','SIGUIJ','SIGATM',]:
            chain_prev = lines[i-1][21]
            if chain == chain_prev:
                res_no_prev = int(lines[i-1][22:26])
                iCode_prev = lines[i-1][26]
                altloc_prev = lines[i-1][16]
            else: ## e.g. 3bve
                res_no_prev = None
                iCode_prev = None
                altloc_prev = None
        else:
            res_no_prev = None
            iCode_prev = None
            altloc_prev = None

        ## skip if altloc or not the first atom in residue
        ## 2zfo (atlloc B in SEQRES)
        skip = False
        if chain in d_ATOMseq.keys():
            if len(d_ATOMseq[chain]['seq']) > 0:
                if (
                    d_ATOMseq[chain]['res_nos'][-1] == res_no and
                    d_ATOMseq[chain]['iCodes'][-1] == iCode
                    ):
                    ## skip if altloc A not added
                    skip = True

        if (
            (
                len(d_ATOMseq[chain]['seq']) == 0 or
                (
                    not (
                        res_no == res_no_prev and
                        iCode == iCode_prev and
                        altloc == altloc_prev ## 2zfo (atlloc B in SEQRES)
                        )
                    )
                ) and
            skip == False
            ):

            res_name_REMARK465 = None
            REMARK465_before_ATOM = False
            if 'REMARK465' in d_header.keys():
                if chain in d_header['REMARK465']['chains'].keys():

                    pass_if = False

                    if res_no_prev != None:
                        ## REMARK465 before ATOM (gap between REMARK465 and next ATOM but not prev ATOM)
                        ## 2nv7
                        if res_no_prev+1 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                            pass_if = True
                        ## REMARK465 before ATOM (gap between REMARK465 and next ATOM and prev ATOM)
                        ## 2pqj (not 2qqh)
                        if set(range(res_no_prev+1,res_no)) & set(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                            pass_if = True

                    if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                        ## REMARK465 before or after ATOM
                        if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                            pass_if = True
                        ## REMARK465 before ATOM
                        if res_no-1 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                            pass_if = True
                        ## REMARK465 before ATOM (no zero residue)
                        ## 1b9n,2fxm
                        if res_no_prev == None and res_no >= 1 and -1 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                            pass_if = True
                        ## REMARK465 before ATOM (first residue = 0 and second residue > 1) ## e.g 2asd,1ca5
                        if res_no_prev == None and res_no > 1 and 0 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                            pass_if = True
##                        ## REMARK465 before ATOM (first residue > 1 and second residue > 1) ## e.g 2h27
##                        if res_no_prev == None and res_no > 1 and res_no > min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
##                            print res_no
##                            print d_header['REMARK465']['chains'][chain]['residues'].keys()
##                            stop
##                            pass_if = True
                        ## REMARK465 before ATOM ## e.g 3bve
                        if res_no_prev == None and res_no > 1:
                            pass_if = True
                        ## REMARK465 before ATOM
                        ## 1sgf,1nu0
                        if res_no_prev == min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                            pass_if = True

                    if pass_if == True:

                        ## REMARK465 before ATOM?
                        ## e.g. 3bef
##                        if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
##
##                            l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']
##                            index1 = self.s_alphabet.index(min(l_iCodes_REMARK465))
##                            index2 = self.s_alphabet.index(max(l_iCodes_REMARK465))+1
##                            l_iCodes_ascending = ','.join(self.s_alphabet[index1:index2]).split(',')
##                            l_iCodes_descending = list(l_iCodes_ascending)
##                            l_iCodes_descending.reverse()
##                            if l_iCodes_REMARK465 == l_iCodes_ascending:
##                                ascending = True
##                                descending = False
##                            elif l_iCodes_REMARK465 == l_iCodes_descending:
##                                ascending = False
##                                descending = True
##
##                            if len(l_iCodes_REMARK465) > 1 and iCode == ' ' and ascending == True:
##                                if min(d_header['REMARK465']['chains'][chain]['residues'].keys()) == res_no: ## e.g. 3e5v
##                                    l_REMARK465_res_nos = [res_no]
##                                else:
##                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no) ## e.g. 1b8m
##                                stop1
##                            elif len(l_iCodes_REMARK465) > 1 and iCode == ' ' and descending == True:
##                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 2ass
##                                stop2
##                            elif iCode != ' ':
##                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 3bef
##                                stop3
##                            ## e.g. 2bvs
##                            elif len(l_iCodes_REMARK465) == 1:
##                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1)
##                                l_REMARK465_res_names = []
##                                for res_no_REMARK465 in l_REMARK465_res_nos:
##                                    if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
##                                        l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
##                                        for iCode_REMARK465 in l_iCodes_REMARK465:
##                                            res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][' ']['res_name']
##                                            l_REMARK465_res_names += [res_name_REMARK465]
##                                iCode_REMARK465 = l_iCodes_REMARK465[0]
##                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_REMARK465]['altlocs'][' ']['res_name']
##                                SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(l_REMARK465_res_names)]
##                                if SEQRES_seq == l_REMARK465_res_names:
##                                    pass
##                                else:
##                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)
##                                stop4
##                            else:
##                                stop
##                        else:
##                            if min(d_header['REMARK465']['chains'][chain]['residues'].keys()) == res_no:
##                                stop_example
####                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)
##                            l_REMARK465_seq = list(d_header['REMARK465']['chains'][chain]['seq'])
                        l_REMARK465_seq = list(d_header['REMARK465']['chains'][chain]['seq'])

                        ##
                        ## how many residues from REMARK465 records are to be used?
                        ##
                        ## e.g. 1bd7
                        bool_break = False
                        l_REMARK465_res_names = []
                        for i_seq_REMARK465 in range(len(l_REMARK465_seq)):
                            seq_REMARK465 = l_REMARK465_seq[i_seq_REMARK465]
                            res_no_REMARK465 = int(seq_REMARK465[:4])
                            iCode_REMARK465 = seq_REMARK465[4]
                            altloc_REMARK465 = ' ' ## altlocs not used in REMARK465 records
                            res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][altloc_REMARK465]['res_name']
                            bool_res_no_REMARK465_descending = False
                            bool_res_no_REMARK465_ascending = False
                            if i_seq_REMARK465 < len(l_REMARK465_seq)-1:
                                if l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465+1,' ',):
                                    bool_res_no_REMARK465_ascending = True
                                if l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465-1,' ',):
                                    bool_res_no_REMARK465_descending = True
                            if i_seq_REMARK465 > 0:
                                if l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465-1,' ',):
                                    bool_res_no_REMARK465_ascending = True
                                if l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465+1,' ',):
                                    bool_res_no_REMARK465_descending = True
                                
                            ## e.g. 1ef0, 3bve, 2iez, many others
                            if res_name_REMARK465 == res_name_ATOM:

##                                if chain == 'A' and res_no == 1 and iCode == ' ' and res_no_REMARK465 != 9999:
##                                    print chain, res_no, iCode
##                                    print 'res_name', res_name_REMARK465, res_name_ATOM
##                                    print 'res_no  ', res_no_REMARK465, res_no
##                                    print 'iCode   ', iCode_REMARK465, iCode
##                                    print
##                                    print '***REM465', l_REMARK465_seq
##                                    print '***ATOM  ', d_coordinates['chains'][chain]['seq']
##                                    print
####                                    stop

                                ##
                                ## not the first ATOM residue
                                ##
                                if len(d_coordinates['chains'][chain]['seq']) > 1:
                                    print 'b', chain, res_no, res_no_REMARK465

                                    ##
                                    ## 1g2w,1l4d
                                    if i_seq_REMARK465+1 == len(l_REMARK465_seq):
                                        bool_REMARK465_continuation_forward = 'Unknown'
                                    else:
                                        if (
                                            i_seq_REMARK465 != 0 ## 1g2w
                                            and
                                            res_no_REMARK465 == int(l_REMARK465_seq[i_seq_REMARK465+1][:4])-1
                                            ):
                                            bool_REMARK465_continuation_forward = True
                                        elif ( ## 2ven
                                            i_seq_REMARK465 == 0 ## 2ven
                                            and
                                            res_no_REMARK465 > res_no_prev
                                            and
                                            iCode == ' '
                                            and
                                            '%4i%1s' %(res_no-1,' ',) in l_REMARK465_seq
                                            ):
                                            bool_REMARK465_continuation_forward = True
                                        elif ( ## 1o6z
                                            i_seq_REMARK465 != 0
                                            and
                                            bool_res_no_REMARK465_ascending == True
                                            and
                                            '%4i%1s' %(res_no-1,' ',) in l_REMARK465_seq
                                            ):
                                            bool_REMARK465_continuation_forward = True
                                        else:
                                            bool_REMARK465_continuation_forward = False
                                    if i_seq_REMARK465 == 0:
                                        bool_REMARK465_continuation_backward = 'Unknown'
                                    else:
                                        if res_no_REMARK465 != int(l_REMARK465_seq[i_seq_REMARK465-1][:4])+1:
                                            bool_REMARK465_continuation_backward = True
                                        else:
                                            bool_REMARK465_continuation_backward = False
                                    if res_no_REMARK465 in [
                                        res_no-1,res_no_prev+1,
                                        ]:
                                        bool_REMARK465_continuation_ATOM = True
                                    else:
                                        bool_REMARK465_continuation_ATOM = False

                                    ##
                                    ## increasing ATOM res_no
                                    if iCode == ' ' and d_coordinates['chains'][chain]['seq'][-2] == '%4i%1s' %(res_no-1,' ',):
                                        ## 1abj, 1sg8
                                        bool_break = True
##                                        s1
                                        break
                                    ## increasing ATOM iCode
                                    elif res_no == int(d_coordinates['chains'][chain]['seq'][-2][:4]) and d_coordinates['chains'][chain]['seq'][-2] == '%4i%1s' %(res_no,self.s_alphabet[self.s_alphabet.index(iCode)-1],):
                                        ## 3ee0
                                        bool_break = True
##                                        s2
                                        break
                                    elif not (
                                        ## increasing res_no from ATOM to REMARK465
                                        bool_REMARK465_continuation_ATOM == True
                                        or
                                        bool_REMARK465_continuation_backward == True
                                        or
                                        bool_REMARK465_continuation_forward == True
                                        ):
                                        bool_break = True
                                        print i_seq_REMARK465
                                        print res_no_REMARK465
                                        print res_no_REMARK465, l_REMARK465_seq
                                        print bool_REMARK465_continuation_ATOM
                                        print bool_REMARK465_continuation_backward
                                        print bool_REMARK465_continuation_forward
##                                        stop
##                                        s3
                                        break
                                    else:
                                        bool_break = False
                                        pass

                                ##
                                ## the first ATOM residue
                                ##
                                elif len(d_coordinates['chains'][chain]['seq']) == 1:
                                    print 'a', chain, res_no, res_no_REMARK465
##                                    if chain == 'B' and res_no == 161:
##                                        print line
##                                        print res_no_REMARK465, res_no
##                                        stop

                                    ##
                                    ## break (current residue is ATOM)
                                    ##

                                    ## res_no ascending
                                    if res_no_REMARK465 > res_no and bool_res_no_REMARK465_descending == False:
                                        ## 2cmj, 2ge8, 2z2q, 3bve
                                        bool_break = True
                                        break
                                    ## iCode descending
                                    elif (
                                        res_no == res_no_REMARK465
                                        and
                                        self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)-1] == iCode
                                        and
                                        (iCode != ' ' or i_seq_REMARK465 != 0) ## e.g. 3e5v
                                        ):
##                                        print self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)+1]
##                                        print l_REMARK465_seq[i_seq_REMARK465+1]
##                                        print '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)-1],)
##                                        stop
##                                        print '%4i%1s' %(res_no,self.s_alphabet[self.s_alphabet.index(iCode)-1],), 'in', l_REMARK465_seq[i_seq_REMARK465:]
##                                        print l_REMARK465_seq[i_seq_REMARK465-1], '==', '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)-1],)
##                                        print i_seq_REMARK465, '==', 0
##                                        stop
                                        ## 1sgi, 3e5v
                                        bool_break = True
                                        break
                                    ## iCode ascending
                                    elif (
                                        l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no,self.s_alphabet[self.s_alphabet.index(iCode)-1],)
                                        ):
                                        ## e.g. 1jy0
                                        bool_break = True
                                        break
##                                    ## prev res_no in REMARK465, next res_no not in REMARK465; WHAT THE HELL????
##                                    elif iCode == ' ' and '%4i%1s' %(res_no-1,' ',) in l_REMARK465_seq and '%4i%1s' %(res_no-1,' ',) not in l_REMARK465_seq:
##                                        bool_break = False
##                                        pass

                                    ##
                                    ## don't break (current residue is REMARK465)
                                    ##
                                    
                                    ## res_no ascending (skip zero) (don't break)
                                    elif (
                                        res_no == 1 and res_no_REMARK465 == -1
                                        and
                                        l_REMARK465_seq[:l_REMARK465_seq.index('%4i%1s' %(-1,' ',))+1] == ['%4i%1s' %(i,' ',) for i in range(int(l_REMARK465_seq[0][:4]),0)]
                                        ):
                                        ## e.g. 1gou, 3csp
                                        bool_break = False
                                        pass
                                    ## res_no ascending (no skip zero) (don't break)
                                    elif (
                                        ## res_no_prev == res_no-1
                                        (
                                            l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465-1,' ',)
##                                            or
##                                            l_REMARK465_seq[i_seq_REMARK465-1] == -1 and res_no_REMARK465 == 1 ## 1oag
                                            )
                                        and
                                        (
                                            ## terminal
                                            res_no-1 == res_no_REMARK465
                                            or
                                            ## not terminal, res_no_next == res_no+1
                                            (
                                                ## last REMARK465 record
                                                i_seq_REMARK465+1 == len(l_REMARK465_seq) ## e.g. 1e1c
                                                or
                                                l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465+1,' ',)
                                                or
                                                (l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(1,' ',) and res_no_REMARK465 == -1) ## 1oag
                                                )
                                            )
                                        ):
                                        bool_break = False
                                        pass
                                    ## res_no descending (don't break)
                                    elif bool_res_no_REMARK465_descending == True:
                                        ## e.g. 1fze,3bve
                                        bool_break = False
                                        pass
                                    ## iCode descending (don't break)
                                    elif (
                                        res_no == res_no_REMARK465
                                        and
                                        ## next REMARK465 = current REMARK465 and lower iCode
                                        l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)-1],)
##                                        and
##                                        ## prev REMARK465 = current REMARK465 and higher iCode
##                                        l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)+1],)
                                        ):
                                        ## e.g. 1fze
                                        bool_break = False
                                        pass
                                    ## iCode ascending (don't break)
                                    elif (
                                        res_no == res_no_REMARK465
                                        and
                                        (
                                            ## prev expected ATOM in REMARK465
                                            ## e.g. 1sgi
                                            (iCode != ' ' and '%4i%1s' %(res_no,self.s_alphabet[self.s_alphabet.index(iCode)-1],) in l_REMARK465_seq[i_seq_REMARK465:])
                                            or
                                            ## next REMARK465 = current REMARK465 and higher iCode
                                            ## e.g. 3e5v
                                            (iCode == ' ' and l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)+1],))
                                            )
                                        and
                                        (
                                            ## prev REMARK465 = current REMARK465 and lower iCode (or current REMARK465 = prev REMARK465 and higher iCode)
                                            ## e.g. 1jy0
                                            l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)-1],)
                                            or
                                            ## first REMARK465
                                            ## e.g. 1jy0
                                            i_seq_REMARK465 == 0
                                            )
                                        ):
                                        bool_break = False
                                        pass
                                    elif i_seq_REMARK465 == 0 and res_no-1 == res_no_REMARK465:
                                        ## e.g. 1oqo
                                        bool_break = False
                                        pass
                                    elif res_no_REMARK465 < res_no and '%4i%1s' %(res_no_REMARK465+1,' ',) == l_REMARK465_seq[i_seq_REMARK465+1]:
                                        ## e.g. 3glm
                                        bool_break = False
                                        pass
                                    else:
                                        print chain, res_no, iCode
                                        print 'res_name', res_name_REMARK465, res_name_ATOM
                                        print 'res_no  ', res_no_REMARK465, res_no
                                        print 'iCode   ', iCode_REMARK465, iCode
                                        print
                                        print '***REM465', l_REMARK465_seq
                                        print '***ATOM  ', d_coordinates['chains'][chain]['seq']
                                        print
                                        print l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)-1],)
                                        print l_REMARK465_seq[i_seq_REMARK465-1],  '%4i%1s' %(res_no_REMARK465,self.s_alphabet[self.s_alphabet.index(iCode_REMARK465)+1],)
                                        print 
                                        stop_new_case
                                if res_no_REMARK465 > res_no and int(l_REMARK465_seq[i_seq_REMARK465-1][:4]) <= res_no: ## e.g. 2iez,2gp9
                                    bool_break = True
                                    break
                            
                                ## end if res_name_REMARK465 == res_name_ATOM
                            if l_REMARK465_res_names+[res_name_REMARK465] != d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)+1]:
                                bool_break = True
                                break
                            l_REMARK465_res_names += [res_name_REMARK465]
                            ## end for i_seq_REMARK465 in range(len(l_REMARK465_seq))
                        if bool_break == True:
                            l_REMARK465_seq = l_REMARK465_seq[:i_seq_REMARK465]

##                        if chain == 'A' and res_no == 14 and iCode == ' ':
##                            print 'REM465', l_REMARK465_res_names+[res_name_REMARK465]
##                            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)+1]
##                            print chain, res_no_REMARK465,iCode_REMARK465,altloc_REMARK465
##                            print l_REMARK465_seq
##                            stop2

                        ##
                        ## REMARK465 before ATOM ?
                        ##
                        if len(l_REMARK465_res_names) > 0:
                            l_SEQRES_res_names = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)]
                            ## multiple residue insertion
                            ## 4htc
                            if len(l_REMARK465_res_names) > 1 and l_REMARK465_res_names == l_SEQRES_res_names:
                                REMARK465_before_ATOM = True
                            ## single residue insertion
                            elif (
                                len(l_REMARK465_res_names) == 1 and
                                l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                                res_name_ATOM != l_SEQRES_res_names[0]
                                ):
                                REMARK465_before_ATOM = True
                            ## single residue insertion, res_no_465 < res_no_ATOM
                            ## 2hu9, 1oqo
                            elif (
                                len(l_REMARK465_res_names) == 1 and
                                l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                                res_name_ATOM == l_SEQRES_res_names[0]# and
##                                res_no > min(l_REMARK465_res_nos)
                                ):
                                REMARK465_before_ATOM = True
                            ## REMARK465 not before ATOM
                            ## e.g. 3bef,1bd7,4htc
                            else:
                                REMARK465_before_ATOM = False
                                s4_example_shouoldnt_be_any_if_none_then_delete_this_if_statement
                            REMARK465_before_ATOM = True
                        ## REMARK465 not before ATOM
                        ## e.g. 2a0q, 3ee0
                        else:
                            REMARK465_before_ATOM = False

                        if not REMARK465_before_ATOM == True: ## e.g. 103l
                            if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                res_no_REMARK465 = res_no
                                l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                iCode_REMARK465 = l_iCodes_REMARK465[0]
                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][' ']['res_name']

            try:
                res_name_SEQRES = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
            except:
                res_name_SEQRES = 'N/A'

            ##
            ## REMARK465 after ATOM
            ##
            if REMARK465_before_ATOM == False and res_name_ATOM == res_name_SEQRES:
                d_ATOMseq,d_header = self.append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,None,False,True,)
            ##
            ## REMARK465 before ATOM (certain)
            ##
            elif REMARK465_before_ATOM == True:# and res_name_ATOM != res_name_SEQRES:# and res_name_REMARK465 == res_name_SEQRES:
                d_ATOMseq,d_header = self.append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_seq,True,True,)
            ##
            ##
            ##
            else:
                
                atom_name = line[12:16].strip()
                print '---'
                print chain, res_no, iCode
                SEQRES_res_name = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
                SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print chain, res_no, iCode
                print line
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3']
                print 'maybe', d_header['SEQRES']['chains'][chain]['seq3'][0], chain, res_no, 'is a MODRES and no MODRES record is present?'
                if 'REMARK465' not in d_header.keys() and atom_name == 'CA' and d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] == {}:
                    stop_remark470_records_missing
                    pass
                elif 'REMARK465' not in d_header.keys() and atom_name == 'N' and res_no == 2 and res_name_SEQRES == 'MET':
                    stop_remark465_records_missing_met1
                    pass
                ## 2zfo
                elif 'REMARK465' not in d_header.keys() and altloc != ' ':
                    if res_name_SEQRES == res_name_ATOM:
                        stop
                    pass
##                elif res_name_REMARK465 == None and min(d_header['REMARK465']['chains'][chain]['residues'].keys()) < res_no: ## e.g. 3c5w
##                    print chain,res_no
##                    print d_header['REMARK465']['chains'][chain]['residues'].keys()
##                    stop_temp_broken
##                    pass
                ## REMARK465 after ATOM
                elif (
                    'REMARK465' in d_header.keys() and
                    (iCode == ' ' or len(d_ATOMseq[chain]['seq']) == 0) and
                    res_no <= min(d_header['REMARK465']['chains'][chain]['residues'].keys()) and
                    res_name_ATOM == SEQRES_res_name
                    ):
                    stop_temp_example_commentout
                    d_ATOMseq[chain]['seq'] += [res_name_ATOM]
                    d_ATOMseq[chain]['res_nos'] += [res_no]
                    d_ATOMseq[chain]['iCodes'] += [iCode]
                    d_ATOMseq[chain]['records'] += [record]
                else:
                    print '*******'
                    print 'ATOM  ', d_ATOMseq[chain]['seq']
                    print 'SEQRES', SEQRES_seq
                    print line
                    print chain,res_no
                    print 'SEQRES', SEQRES_res_name
                    print 'ATOM  ', res_name_ATOM
                    print 'REMARK', res_name_REMARK465
                    print
##                    print l_REMARK465_seq
                    print d_ATOMseq[chain]['res_nos']
                    stop_N_terminal

            if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]:
                print '*******'
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print line
                print res_name_ATOM
                print res_name_REMARK465
                print chain, res_no, iCode
                stop_sequence_difference

        return d_ATOMseq        


    def check_if_SEQRESres(self,res_name,record,d_header,chain,res_no,iCode):

        if res_name not in self.d_res1.keys()+self.l_nucleotides and record == 'ATOM':
            print res_name,record
            stop

        if res_name in self.d_res1.keys()+self.l_nucleotides and record in ['ATOM','REMARK465',]:
            return True

        MODRES = False
        if 'MODRES' in d_header.keys():
            if chain in d_header['MODRES'].keys():
                if res_no in d_header['MODRES'][chain].keys():
                    if iCode in d_header['MODRES'][chain][res_no].keys():
                        if res_name in d_header['MODRES'][chain][res_no][iCode]:
                            MODRES = True
                        else:
                            print chain, res_no,iCode
                            print res_name,d_header['MODRES'][chain][res_no][iCode]
                            stop

        ## cases for which MODRES is not used
        if MODRES == False:
            ## functional N-terminal groups
            if res_name in self.l_terminalmodres and res_name in [d_header['SEQRES']['chains'][chain]['seq3'][0],d_header['SEQRES']['chains'][chain]['seq3'][1],d_header['SEQRES']['chains'][chain]['seq3'][-1],]:
                MODRES = True
            ## l-peptide linkers
            if res_name in self.l_lpeptidelinkers and res_name in d_header['SEQRES']['chains'][chain]['seq3']:
                MODRES = True
            ## d-amino acids
            if res_name in self.l_dpeptidelinkers and res_name in d_header['SEQRES']['chains'][chain]['seq3']:
                MODRES = True

        return MODRES


    def build_dictionary_of_molecules(self,d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=True,):

        d_adjacency_backward = {} ## C --> O
        d_adjacency_forward = {} ## O --> C
        matrix_adjacency = []
        l_connections = []
        set_monomers = set()
        d_connections = {}
        l_monomers = []
        for hetID1 in d_CONECT.keys():
            for chain1 in d_CONECT[hetID1].keys():
                for res_no1 in d_CONECT[hetID1][chain1].keys():
                    for iCode1 in d_CONECT[hetID1][chain1][res_no1].keys():
                        for atom_no1 in d_CONECT[hetID1][chain1][res_no1][iCode1].keys():
                            
                            atom_name1 = d_atomnos[atom_no1]['atom_name']
                            element1 = d_atomnos[atom_no1]['element']
                            if element1 in self.l_atoms_metal:
                                continue
##                            if atom_name1 in self.l_atoms_metal:
##                                continue
##                            if atom_name1[:2] in self.l_atoms_metal:
##                                continue
##                            if atom_name1[:1] in self.l_atoms_metal:
##                                print atom_name1
##                                stop
##                                continue

                            ## check ...
                            if s_pdb not in ['1nkm','1nju',]:
                                if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']) > 1:
                                    atom_names = set([d_atomnos[atom_no1]['atom_name']])
                                    for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                        atom_names |= set([d_atomnos[atom_no]['atom_name']])
                                    if atom_names != set(['SG']):
                                        ## one atom of one residue connected to more than one atom of other residue(s)
                                        if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']) == 2:
                                            tmpatom_no1 = d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][0]
                                            tmpatom_no2 = d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][1]
                                            ## if not altloc of same residue
                                            if not (
                                                d_atomnos[tmpatom_no1]['chain'] == d_atomnos[tmpatom_no2]['chain'] and
                                                d_atomnos[tmpatom_no1]['res_no'] == d_atomnos[tmpatom_no2]['res_no'] and
                                                d_atomnos[tmpatom_no1]['iCode'] == d_atomnos[tmpatom_no2]['iCode'] and
                                                d_atomnos[tmpatom_no1]['atom_name'] == d_atomnos[tmpatom_no2]['atom_name'] and
                                                (d_atomnos[tmpatom_no1]['altloc'] == 'A' and d_atomnos[tmpatom_no2]['altloc'] == 'B')
                                                ):
                                                print 'DOUBLE CONNECTION'
                                                print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                                print 'res_names', d_atomnos[tmpatom_no1]['res_name'], d_atomnos[tmpatom_no2]['res_name']
                                                print 'chains', d_atomnos[tmpatom_no1]['chain'], d_atomnos[tmpatom_no2]['chain']
                                                print 'res_nos', d_atomnos[tmpatom_no1]['res_no'], d_atomnos[tmpatom_no2]['res_no']
                                                print 'iCodes', d_atomnos[tmpatom_no1]['iCode'], d_atomnos[tmpatom_no2]['iCode']
                                                print 'atom_names', d_atomnos[tmpatom_no1]['atom_name'], d_atomnos[tmpatom_no2]['atom_name']
                                                print 'altlocs', d_atomnos[tmpatom_no1]['altloc'], d_atomnos[tmpatom_no2]['altloc']
                                                print '*****'
                                                print (d_atomnos[tmpatom_no1]['altloc'] == 'A' and d_atomnos[tmpatom_no2]['altloc'] == 'B')
                                                print hetID1, chain1, res_no1, iCode1, atom_no1
                                                print hetID1, chain1, res_no1, iCode1, atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']
                                                print s_pdb
                                                for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                    print atom_no, d_atomnos[atom_no]
                                                try:
                                                    print self.cluster, 'cluster'
                                                except:
                                                    None
                                                notexpected_different_atllocs_connected_or_valence_exceeded
                                        else:
                                            if (
                                                (not d_atomnos[d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][-1]]['altloc'] == self.s_alphabet[len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'])])
                                                ):
                                                print atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                                print hetID1, chain1, res_no1, iCode1, atom_no1
                                                atom_names = set([d_atomnos[atom_no1]['atom_name']])
                                                for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                    atom_names |= set([d_atomnos[atom_no]['atom_name']])
                                                atom_names = atom_names - set(self.l_atoms_metal)
                                                print hetID2, atom_names
                                                error = True
                                                for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                    print d_atomnos[atom_no]
                                                    if d_atomnos[atom_no]['res_name'] in ['TML','CYO',]:
                                                        error = False
                                                subtract = 0

                                                ## metal connections
                                                for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                    print atom_no, d_atomnos[atom_no]['element']
                                                    if d_atomnos[atom_no]['element'] in self.l_atoms_metal:
                                                        subtract += 1
                                                if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'])-subtract == 1:
                                                    error = False

                                                print s_pdb
                                                if error == True:
                                                    print self.cluster, 'cluster'
                                                    notexpected

                            l_hetIDs_long_atom_names = [
                                ## pyranoses
                                'XYP','G6D','SIA',
                                'FCT','DAN',
                                ## furanoses
                                'AHR','HPD','3DR','FUB','AAB',
                                ## dissacharides
                                'DAF','DCB','FXP',
                                ## benzoxazinoids (hydroxamic acid)
                                'HBO', ## DIMBOA from Maize
                                ## (tetra)pyrrole
                                'DBV','PEB','OPP','ZNH','BCL',
                                ## pyrrole ring
                                'YRR',
                                ## phosphate group
                                'FHP','4IP',
                                ## other
                                'DPM','OAS','CSC','780','PED',
                                ## benzene ring
                                'TMM','PLR','785','762','BMZ','AEN',
                                ## p-Coumaric acid (phenyl propanoid)
                                'HC4',
                                ## nucleobases/nucleosides/nucleotides
                                'DA','UMP','PGD','DT','BZG','8OG','ADP','A','TSP','CMP','AMP','FOX','PPU',
                                ## adenine
                                '1MA','MA7',
                                ## guanine
                                'DG','GDP','OMG','GTP','G','SGP',
                                ## cytidine
                                'C','DC','OMC','DOC','CSF','5CM','ME6','CTP',
                                ## uridine
                                '5IU','PSU','OMU','U','UR3','5BU',
##                                        ## dinucleotides
##                                        'NAP',
                                ]
                            ## check 2a
                            if atom_name1[-1] in ['A','B','H',"'"] and len(atom_name1) > 2:
                                if hetID1 not in l_hetIDs_long_atom_names:
                                    print 'hetID', hetID1
                                    print 'atom_name', atom_name1
                                    print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                    print 'atom_no', atom_no1
                                    print 'res_no', res_no1
                                    stop_new_long_atom_name
                                atom_name1_no = atom_name1[1:-1]
                            else:
                                atom_name1_no = atom_name1[1:]

                            for atom_no2 in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                hetID2 = d_atomnos[atom_no2]['res_name']
                                chain2 = d_atomnos[atom_no2]['chain']
                                res_no2 = d_atomnos[atom_no2]['res_no']
                                iCode2 = d_atomnos[atom_no2]['iCode']

                                atom_name2 = d_atomnos[atom_no2]['atom_name']
                                element2 = d_atomnos[atom_no2]['element']
                                if element2 in self.l_atoms_metal:
                                    continue
##                                if atom_name2 in self.l_atoms_metal:
##                                    continue
##                                if atom_name2[:2] in self.l_atoms_metal:
##                                    continue
##                                if atom_name2[:1] in self.l_atoms_metal:
##                                    print atom_name2
##                                    stop
##                                    continue

                                ## check 2b
                                if atom_name2[-1] in ['A','B','H',"'"] and len(atom_name2) > 2:
                                    if hetID2 not in l_hetIDs_long_atom_names:
                                        print chain2, res_no2, iCode2, hetID2, atom_name2, d_CONECT[hetID2][chain2][res_no2][iCode2][atom_no2]
                                        print hetID2
                                        stop_new_long_atom_name
                                    atom_name2_no = atom_name2[1:-1]
                                else:
                                    atom_name2_no = atom_name2[1:]

                                modres1 = self.determine_if_modres(d_header, d_coordinates, chain1, res_no1, iCode1, hetID1)
                                if modres1 == True:
                                    continue
                                modres2 = self.determine_if_modres(d_header, d_coordinates, chain2, res_no2, iCode2, hetID2)
                                if modres2 == True:
                                    continue

                                if element1 == 'C' and element2 == 'C':
                                    ## unique carbon-carbon bonds
                                    if (
                                        (hetID1 == 'TYR' and atom_name1 == 'CE1' and hetID2 == 'TRP' and atom_name2 == 'CH2')
                                        or
                                        (hetID2 == 'TYR' and atom_name2 == 'CE1' and hetID1 == 'TRP' and atom_name1 == 'CH2')

                                        or

                                        (hetID2 == 'PHE' and atom_name2 == 'C' and hetID1 == 'PHE' and atom_name1 == 'C')

                                        or

                                        (hetID2 == 'ARG' and atom_name2 == 'C' and hetID1 == 'CH2' and atom_name1 == 'C')
                                        or
                                        (hetID1 == 'ARG' and atom_name1 == 'C' and hetID2 == 'CH2' and atom_name2 == 'C')

                                        or

                                        (hetID2 == 'ALA' and atom_name2 == 'C' and hetID1 == 'BAS' and atom_name1 == 'C1')
                                        or
                                        (hetID1 == 'ALA' and atom_name1 == 'C' and hetID2 == 'BAS' and atom_name2 == 'C1')

                                        or

                                        (hetID1 == 'LOL' and atom_name1 == 'C' and hetID2 == 'ALQ' and atom_name2 == 'CM')
                                        or
                                        (hetID2 == 'LOL' and atom_name2 == 'C' and hetID1 == 'ALQ' and atom_name1 == 'CM')
                                        ):
                                        pass
                                    else:
                                        print hetID1, chain1, res_no1, atom_name1
                                        print hetID2, chain2, res_no2, atom_name2
                                        print atom_no1, atom_no2
                                        print s_pdb
                                        try:
                                            print self.cluster, 'cluster'
                                        except:
                                            None
                                        ## unexpected carbon-carbon
                                        notexpected_carboncarbon_connection

                                if verbose == True and hetID1 != 'CYS' and hetID2 != 'CYS' and atom_name1 != 'SG' and atom_name2 != 'SG':
                                    print hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_name1_no, atom_name2_no
                                #####################################
                                ## posttranslational modifications ##
                                #####################################
                                ##
                                ## peptide bond
                                ##
                                if (
                                    hetID1 in self.d_res1.keys()+self.l_dpeptidelinkers+self.l_lpeptidelinkers and hetID2 in self.d_res1.keys()+self.l_dpeptidelinkers+self.l_lpeptidelinkers and
                                    atom_name1 == 'C' and atom_name2[0] == 'N'
                                    ):
                                    bond = 'C,N'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    hetID1 in self.d_res1.keys()+self.l_dpeptidelinkers+self.l_lpeptidelinkers and hetID2 in self.d_res1.keys()+self.l_dpeptidelinkers+self.l_lpeptidelinkers and
                                    atom_name2 == 'C' and atom_name1[0] == 'N'
                                    ):
                                    bond = 'C,N'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ##
                                ## cysteine/glutathione disulphide bond
                                ##
                                elif hetID1 in ['CYS','GSH',] and atom_name1[:2] == 'SG' and atom_name2[0] == 'S':
                                    bond = 'S,S'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 in ['CYS','GSH',] and atom_name2[:2] == 'SG' and atom_name1[0] == 'S':
                                    bond = 'S,S'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ##
                                ## lysine, cysteine (e.g. palmitoylation), histidine, tryptophan, tyrosine (e.g. adenyaltion)
                                ##
                                elif (
                                    (hetID1 == 'LYS' and atom_name1 == 'NZ'  and atom_name2[0] in ['C','P']) or
                                    (hetID1 == 'CYS' and atom_name1 == 'SG'  and atom_name2[0] == 'C') or
                                    (hetID1 == 'TRP' and atom_name1 == 'CZ3' and atom_name2[0] == 'O') or
##                                    (hetID1 == 'PRO' and atom_name1 == 'C' and atom_name2[0] == 'N') or
                                    (hetID1 == 'HIS' and atom_name1 == 'NE2' and atom_name2[0] in ['C','P']) or
                                    (hetID1 == 'HIS' and atom_name1 in ['ND1',] and atom_name2[0] in ['O',]) or
                                    (hetID1 == 'GLU' and atom_name1 in ['OE1','OE2',]) or
                                    (hetID1 == 'ASP' and atom_name1 in ['OD','OD1','OD2',]) or
                                    (hetID1 == 'TYR' and atom_name1 in ['OH',] and atom_name2[0] in ['C','P'])
                                    ):
                                    bond = atom_name1[0]+','+atom_name2[0]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    (hetID2 == 'LYS' and atom_name2 == 'NZ'  and atom_name1[0] in ['C','P']) or
                                    (hetID2 == 'CYS' and atom_name2 == 'SG'  and atom_name1[0] == 'C') or
                                    (hetID2 == 'TRP' and atom_name2 == 'CZ3' and atom_name1[0] == 'O') or
##                                    (hetID2 == 'PRO' and atom_name2 == 'C' and atom_name1[0] == 'N') or
                                    (hetID2 == 'HIS' and atom_name2 == 'NE2' and atom_name1[0] in ['C','P']) or
                                    (hetID2 == 'HIS' and atom_name2 in ['ND1',] and atom_name1[0] in ['O',]) or
                                    (hetID2 == 'GLU' and atom_name2 in ['OE1','OE2',]) or
                                    (hetID2 == 'ASP' and atom_name2 in ['OD','OD1','OD2',]) or
                                    (hetID2 == 'TYR' and atom_name2 in ['OH',] and atom_name1[0] in ['C','P'])
                                    ):
                                    bond = atom_name2[0]+','+atom_name1[0]
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2

                                ##
                                ## unique bonds
                                ##
                                    
                                ## unique tryptophan peroxidase pi electron porphyrin interaction (make more general solution for small molecules)
                                elif hetID1 == 'TRP' and atom_name1 == 'NE1' and hetID2 == 'PEO' and atom_name2[0] == 'O':
                                    bond = 'N'+','+atom_name2[1:]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'TRP' and atom_name2 == 'NE1' and hetID1 == 'PEO' and atom_name1[0] == 'O':
                                    bond = 'N'+','+atom_name1[1:]
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (1mk8, cytochrome c peroxidase)
                                elif hetID1 == 'TYR' and atom_name1 == 'CE2' and hetID2 == 'MET' and atom_name2 == 'SD':
                                    bond = 'CE1,SD'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'TYR' and atom_name2 == 'CE2' and hetID1 == 'MET' and atom_name1 == 'SD':
                                    bond = 'CE1,SD'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                elif (hetID1 == 'TYR' and atom_name1 == 'CE1' and hetID2 == 'TRP' and atom_name2 == 'CH2'):
                                    bond = 'CE,CH2'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'TYR' and atom_name2 == 'CE1' and hetID1 == 'TRP' and atom_name1 == 'CH2'):
                                    bond = 'CE,CH2'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (1gge, catalase/peroxidase)
                                elif hetID1 == 'TYR' and atom_name1 == 'CB' and hetID2 == 'HIS' and atom_name2 == 'ND1':
                                    bond = 'CB,ND1'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'TYR' and atom_name2 == 'CB' and hetID1 == 'HIS' and atom_name1 == 'ND1':
                                    bond = 'CB,ND1'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (2azc, inhibitor)
                                elif (hetID1 == 'PHE' and atom_name1 == 'C' and hetID2 == 'PHE' and atom_name2 == 'C'):
                                    bond = 'C,C'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'PHE' and atom_name2 == 'C' and hetID1 == 'PHE' and atom_name1 == 'C'):
                                    bond = 'C,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (1nkm, inhibitor)
                                elif (hetID1 == 'ALA' and atom_name1 == 'C' and hetID2 == 'BAS' and atom_name2 == 'C1'):
                                    bond = 'C,C'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'ALA' and atom_name2 == 'C' and hetID1 == 'BAS' and atom_name1 == 'C1'):
                                    bond = 'C,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

##                                ## unique (1nkm, inhibitor)
##                                elif (hetID1 == 'SER' and atom_name1 == 'OG' and hetID2 == 'ALA' and atom_name2 == 'C'):
##                                    bond = 'O,C'
##                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
##                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
##                                elif (hetID2 == 'SER' and atom_name2 == 'OG' and hetID1 == 'ALA' and atom_name1 == 'C'):
##                                    bond = 'O,C'
##                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
##                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (1abj, C-terminal alkylation)
                                elif (hetID1 == 'ARG' and atom_name1 == 'C' and hetID2 == 'CH2' and atom_name2 == 'C'):
                                    bond = 'C,C'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'ARG' and atom_name2 == 'C' and hetID1 == 'CH2' and atom_name1 == 'C'):
                                    bond = 'C,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (1fkn, L-peptide linkers without parent aa IDs)
                                elif (hetID1 == 'LOL' and atom_name1 == 'C' and hetID2 == 'ALQ' and atom_name2 == 'CM'):
                                    bond = 'C,C'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'LOL' and atom_name2 == 'C' and hetID1 == 'ALQ' and atom_name1 == 'CM'):
                                    bond = 'C,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ## unique (2hg5, GOA linker)
                                elif (hetID1 == 'ASP' and atom_name1 == 'N' and hetID2 == 'GOA' and atom_name2 == 'C1'):
                                    bond = 'C,N'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'ASP' and atom_name2 == 'N' and hetID1 == 'GOA' and atom_name1 == 'C1'):
                                    bond = 'C,N'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                elif (hetID1 == 'TYR' and atom_name1 == 'C' and hetID2 == 'GOA' and atom_name2 == 'O2'):
                                    bond = 'C,O'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (hetID2 == 'TYR' and atom_name2 == 'C' and hetID1 == 'GOA' and atom_name1 == 'O2'):
                                    bond = 'C,O'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                                ##
                                ## errors
                                ##
                                    
                                ## unexpected carbon-carbon
                                elif atom_name1[0] == 'C' and atom_name2[0] == 'C':
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected_conn1
                                ## unexpected oxygen-oxygen
                                elif atom_name1[0] == 'O' and atom_name2[0] == 'O':
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    if hetID1 in ['GAL',] and hetID2 in ['NAG'] and atom_name1 == 'O5' and atom_name2 == 'O4':
                                        fd = open('remediation_oxygen_valence_exceeded.txt','a')
                                        fd.write('%s %s %s %s %s %s %s %s %s %s %s\n' %(
                                            s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2,
                                            ))
                                        fd.close()
                                        pass
                                    else:
                                        expected_conn2_oxygenoxygen
                                ## unexpected other
                                elif atom_name1[0] == atom_name2[0]:
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected_conn3
                                ##
                                ## O-glycosylation
                                ##
                                elif (
                                    hetID1 == 'SER' and atom_name1 == 'OG' or
                                    hetID1 == 'THR' and atom_name1 == 'OG1'
                                    ):
                                    bond = 'O'+','+atom_name2[1:]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    hetID2 == 'SER' and atom_name2 == 'OG' or
                                    hetID2 == 'THR' and atom_name2 == 'OG1'
                                    ):
                                    bond = 'O'+','+atom_name1[1:]
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                ##
                                ## N-glycosylation
                                ##
                                elif atom_name1 == 'ND2' and hetID1 == 'ASN':
                                    bond = 'N'+','+atom_name2[1:]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif atom_name2 == 'ND2' and hetID2 == 'ASN':
                                    bond = 'N'+','+atom_name1[1:]
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
##                                ##
##                                ## temporary... hopefully... dependent on the remediation team
##                                ## 1-1 connections?
##                                ##
##                                elif hetID1 in ['KDO','KDA','ABU','BAL','DIB','PYB','IMT'] and hetID2 in ['KDO','KDA','ABU','BAL','DIB','PYB','IMT']:
##                                    print atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
##                                    print atom_no2, d_CONECT[hetID2][chain2][res_no2][iCode2][atom_no2]
##                                    if (
##                                        (atom_name2[0] == 'C' and atom_name1[0] in ['O']) or
##                                        (atom_name1[0] == 'C' and atom_name2[0] in ['N'])
##                                        ):
##                                        bond = atom_name1_no+','+atom_name2_no
##                                        monomer1 = chain1+str(res_no1).zfill(4)+iCode1
##                                        monomer2 = chain2+str(res_no2).zfill(4)+iCode2
##                                        print bond, monomer1, monomer2, atom_name1, atom_name2
##                                    elif (
##                                        (atom_name1[0] == 'C' and atom_name2[0] in ['O']) or
##                                        (atom_name2[0] == 'C' and atom_name1[0] in ['N'])
##                                        ):
##                                        bond = atom_name2_no+','+atom_name1_no
##                                        monomer1 = chain2+str(res_no2).zfill(4)+iCode2
##                                        monomer2 = chain1+str(res_no1).zfill(4)+iCode1
##                                    else:
##                                        print hetID1, hetID2, atom_name1, atom_name2
##                                        notexpected
##                                    if monomer1 > monomer2:## and not (hetID1 in ['BAL','DIB'] or hetID2 in ['BAL','DIB']):
##                                        print hetID1, hetID2, monomer1, monomer2
##                                        notexpected
                                ##
                                ## nucleotide phospo(thio)diester bonds (PTR,DT;)
                                ##
                                elif (atom_name1[0] in ['O','S'] and atom_name2 in ['P','P2']) or (atom_name1 == 'O3P' and atom_name2[0] == 'C'):
                                    bond = 'O,P'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (atom_name2[0] in ['O','S'] and atom_name1 in ['P','P2']) or (atom_name2 == 'O3P' and atom_name1[0] == 'C'):
                                    bond = 'O,P'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ##
                                ## bond (move to bottom after error finding...)
                                ##
                                elif atom_name1 in ['S'] and atom_name2 == 'N':
                                    bond = atom_name1[0]+','+atom_name2[0]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    if monomer1 > monomer2:
                                        print monomer1, monomer2
                                        notexpected
                                elif atom_name2 in ['S'] and atom_name1 == 'N':
                                    bond = atom_name2[0]+','+atom_name1[0]
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    if monomer1 > monomer2:
                                        print monomer1, monomer2
                                        notexpected
                                ##
                                ## N-terminal modification (methylation, acetylation, myristoylation, formylation)
                                ##
                                elif (
                                    hetID1 in ['ACE','MYR','OHE','OME','CH2','FMT',] and
##                                    chain1 == chain2 and
                                    (
                                        res_no2 == 1 or
                                        res_no1 == min(d_coordinates['chains'][chain1]['residues'].keys())
                                        ) and
                                    atom_name2 == 'N' and
                                    atom_name1[0] == 'C'
                                    ):
                                    bond = 'N,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                elif (
                                    hetID2 in ['ACE','MYR','OHE','OME','CH2','FMT',] and
##                                    chain1 == chain2 and
                                    (
                                        res_no1 == 1 or
                                        res_no2 == min(d_coordinates['chains'][chain2]['residues'].keys())
                                        ) and
                                    atom_name1 == 'N' and
                                    atom_name2[0] == 'C'
                                    ):
                                    bond = 'N,C'
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                ##
                                ## C-terminal modification (amination/amidation)
                                ##
                                elif (
                                    hetID1 in ['NH2',] and
                                    chain1 == chain2 and
                                    res_no1 == max(d_coordinates['chains'][chain1]['residues'].keys()) and
                                    res_no2 == max(set(d_coordinates['chains'][chain2]['residues'].keys())-set([max(d_coordinates['chains'][chain2]['residues'].keys())])) and
                                    atom_name2 == 'C' and
                                    atom_name1[0] == 'N'
                                    ):
                                    bond = 'N,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                elif (
                                    hetID2 in ['NH2',] and
                                    chain1 == chain2 and
                                    res_no2 == max(d_coordinates['chains'][chain2]['residues'].keys()) and
                                    res_no1 == max(set(d_coordinates['chains'][chain1]['residues'].keys())-set([max(d_coordinates['chains'][chain1]['residues'].keys())])) and
                                    atom_name1 == 'C' and
                                    atom_name2[0] == 'N'
                                    ):
                                    bond = 'N,C'
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                ##
                                ## cyclization
                                ##
                                elif (
                                    hetID1 in self.d_res1.keys() and
                                    hetID2 in self.d_res1.keys() and
                                    chain1 == chain2 and
                                    res_no1 == max(d_coordinates['chains'][chain1]['residues'].keys()) and
                                    res_no2 == min(d_coordinates['chains'][chain2]['residues'].keys()) and
                                    atom_name1 == 'C' and
                                    atom_name2 == 'N'
                                    ):
                                    bond = 'N,C'
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                elif (
                                    hetID1 in self.d_res1.keys() and
                                    hetID2 in self.d_res1.keys() and
                                    chain1 == chain2 and
                                    res_no2 == max(d_coordinates['chains'][chain2]['residues'].keys()) and
                                    res_no1 == min(d_coordinates['chains'][chain1]['residues'].keys()) and
                                    atom_name2 == 'C' and
                                    atom_name1 == 'N'
                                    ):
                                    bond = 'N,C'
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                ## sulphenyl-amide intermediate (e.g. 1oem)
                                elif (
                                    hetID1 == 'CYS' and
                                    hetID2 == 'SER' and
                                    chain1 == chain2 and
                                    res_no2-res_no1 == 1 and
                                    atom_name1 == 'SG' and
                                    atom_name2 == 'N'
                                    ):
                                    bond = 'S,N'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    hetID2 == 'CYS' and
                                    hetID1 == 'SER' and
                                    chain1 == chain2 and
                                    res_no1-res_no2 == 1 and
                                    atom_name2 == 'SG' and
                                    atom_name1 == 'N'
                                    ):
                                    bond = 'S,N'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ## trp-tyr crosslink (e.g. 1mk8)
                                elif (
                                    hetID1 == 'TRP' and
                                    hetID2 == 'TYR' and
                                    chain1 == chain2 and
                                    res_no2-res_no1 == 1 and
                                    atom_name1 == 'NE1' and
                                    atom_name2 == 'CE1'
                                    ):
                                    bond = 'N,C'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    hetID2 == 'TRP' and
                                    hetID1 == 'TYR' and
                                    chain1 == chain2 and
                                    res_no1-res_no2 == 1 and
                                    atom_name2 == 'NE1' and
                                    atom_name1 == 'CE1'
                                    ):
                                    bond = 'N,C'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ## N,C bond
                                elif (
                                    chain1 == chain2 and
                                    atom_name1 == 'N' and atom_name2[0] == 'C'
                                    ):
                                    bond = 'N,C'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    chain1 == chain2 and
                                    atom_name2 == 'N' and atom_name1[0] == 'C'
                                    ):
                                    bond = 'N,C'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ## error
                                elif (
                                    atom_name1_no != '1' and atom_name2_no != '1' and
                                    hetID1 not in ['SIA','SLB','XYP','KDO','KDA','KDB',] and hetID2 not in ['SIA','SLB','XYP','KDO','KDA','KDB',]
                                    ): ##  and atom_name1[0] == atom_name2[0]
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, iCode1, iCode2
                                    print atom_name1, atom_name2, atom_no1, atom_no2
                                    print min(d_coordinates['chains'][chain1]['residues'].keys())
                                    print min(d_coordinates['chains'][chain2]['residues'].keys())
                                    print max(d_coordinates['chains'][chain1]['residues'].keys())
                                    print max(d_coordinates['chains'][chain2]['residues'].keys())
                                    try:
                                        print 'cluster', self.cluster
                                    except:
                                        None
                                    notexpected_bond_not11
                                ##
                                ## glycosyl and peptide bonds
                                ##
                                ## 1,1-glycoside bond (e.g. trehalose in 1do1)
                                elif int(atom_name1_no) == 1 and int(atom_name2_no) == 1:
                                    if atom_name1[0] == 'O' and atom_name2[0] == 'C' and hetID2 in ['GLC','BGC',]:
                                        monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                        monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    elif atom_name1[0] == 'C' and atom_name2[0] == 'O' and hetID1 in ['GLC','BGC',]:
                                        monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                        monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    else:
                                        print atom_name1, atom_name2
                                        stopstop
                                    print monomer1, monomer2
                                    if hetID2 in self.d_saccharides.keys() and hetID1 not in self.d_saccharides.keys():
                                        if atom_name2[:2] == 'C1' and atom_name1[0] == 'O':
                                            if monomer1 != chain1+str(res_no1).zfill(4)+iCode1: ## temp!!!
                                                notexpected
                                            monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                            monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                            bond = atom_name2_no+','+atom_name1_no
                                            if self.verbose == True:
                                                print 'a', bond
                                        else:
                                            print hetID1, hetID2, atom_no2
                                            notexpected
                                    elif hetID1 in self.d_saccharides.keys() and hetID2 not in self.d_saccharides.keys():
                                        if atom_name1[:2] == 'C1' and atom_name2[0] == 'O':
                                            if monomer1 != chain2+str(res_no2).zfill(4)+iCode2: ## temp!!!
                                                notexpected
                                            monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                            monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                            bond = atom_name1_no+','+atom_name2_no
                                            if self.verbose == True:
                                                print 'b', bond
                                        else:
                                            print hetID1, hetID2, atom_name1, atom_name2
                                            notexpected
                                    elif hetID1 == 'GLC' and hetID2 == 'GLC':
                                        bond = '1,1'
                                        if self.verbose == True:
                                            print 'c', bond
                                        pass ## trehalose (e.g. 1do1)
                                    else:
                                        print s_pdb
                                        print hetID1,hetID2
                                        stop
                                    print monomer1, monomer2, bond
                                ## glycosyl 1
                                elif (
                                    int(atom_name1_no) != 1 and
                                    (
                                        (atom_name1[0] in ['O','N','S',] and atom_name2[:2] == 'C1')
                                        or
                                        (atom_name1[0] == 'C' and atom_name2[:2] in ['O1','N1',])
                                        or
                                        ## sialic acid
                                        (atom_name1[0] == 'O' and atom_name2 == 'C2' and hetID2 in ['SIA','SLB'])
                                        or
                                        (atom_name1[0] == 'C' and atom_name2 == 'O2' and hetID2 in ['SIA','SLB'])
                                        or
                                        (atom_name1 == 'C4B' and atom_name2 == 'O4A' and hetID1 in ['XYP',] and hetID2 in ['XYP',])
                                        or
                                        (atom_name1 in ['C2','O2',] and atom_name2 in ['O4','O8','C4','C8',] and hetID1 in ['KDO','KDA','KDB',] and hetID2 in ['KDA','KDO','KDB',])
                                        )
                                    ):
                                    bond = atom_name2_no+','+atom_name1_no
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    if self.verbose == True:
                                        print 'd', bond
                                ## glycosyl 2
                                elif (
                                    int(atom_name2_no) != 1 and
                                    (
                                        (atom_name2[0] in ['O','N','S',] and atom_name1[:2] == 'C1')
                                        or
                                        (atom_name2[0] == 'C' and atom_name1[:2] in ['O1','N1'])
                                        or
                                        ## sialic acid
                                        (atom_name2[0] == 'O' and atom_name1 == 'C2' and hetID1 in ['SIA','SLB'])
                                        or
                                        (atom_name2[0] == 'C' and atom_name1 == 'O2' and hetID1 in ['SIA','SLB'])
                                        or
                                        (atom_name2 == 'C4B' and atom_name1 == 'O4A' and hetID2 in ['XYP',] and hetID1 in ['XYP',])
                                        or
                                        (atom_name2 in ['C2','O2',] and atom_name1 in ['O4','O8','C4','C8',] and hetID2 in ['KDO','KDA','KDB',] and hetID1 in ['KDA','KDO','KDB',])
                                        )
                                    ):
                                    bond = atom_name1_no+','+atom_name2_no
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    if self.verbose == True:
                                        print 'e', bond
                                ## error
                                else:
                                    print s_pdb, hetID1, hetID2
                                    print chain1, chain2, res_no1, res_no2, iCode1, iCode2
                                    print atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    notexpected


                                ##
                                ## write adjacency dictionaries if not disulphide bonds
                                ##
                                if atom_name1 != 'SG' and atom_name2 != 'SG':
                                    set_monomers |= set([monomer1,monomer2])
                                    d_connections[monomer1] = monomer2

                                    if not monomer1 in d_adjacency_forward.keys():
                                        d_adjacency_forward[monomer1] = {}
                                    if monomer2 in d_adjacency_forward[monomer1]: ## temp!!!
                                        if bond != d_adjacency_forward[monomer1][monomer2]: ## temp!!!
                                            print s_pdb, bond, d_adjacency_forward[monomer1][monomer2]
                                            print 'monomers', monomer1, monomer2
                                            print 'atoms', atom_no1, atom_no2
                                            print d_adjacency_forward
                                            print d_adjacency_backward
                                            print s_pdb
                                            notexpected ## temp!!!
                                    d_adjacency_forward[monomer1][monomer2] = bond
                                        
                                    if not monomer2 in d_adjacency_backward.keys():
                                        d_adjacency_backward[monomer2] = []
                                    if not monomer1 in d_adjacency_backward[monomer2]:
                                        d_adjacency_backward[monomer2] += [monomer1]

                                ## end of loop over atom_no2
                            ## end of loop over atom_no1
            ## end of loop over hetID1


        ## identify roots of saccharides (simplified back trace)
        set_roots = set_monomers-set(d_adjacency_backward.keys())

        ##
        ## trace saccharide tree forward and identify the saccharide molecule
        ## e.g. 1h3u.pdb, 2c4a.pdb
        ##
        ## monomers (combination of chain,resno,iCode) are unique
        ## monomer hetIDs and bonds are not unique
        ##
        ## path of bonds are unique, but path of monomers are used
        ## to avoid traveling down branches multiple times
        ##
        ## that deeply nested dictionaries can be compared was checked
        ## by modifying nested values of the dictionaries a,b and then comparing them
        ##a = {'A86 ': {'bonds': {'1,N': {'monomer': ' 86A', 'bonds': {}}}}, 'A146 ': {'bonds': {'1,N': {'monomer': ' 146A', 'bonds': {}}}}, 'A200 ': {'bonds': {'1,N': {'monomer': ' 200A', 'bonds': {'1,4': {'monomer': ' 200B', 'bonds': {'1,4': {'monomer': ' 200C', 'bonds': {'1,6': {'monomer': ' 200G', 'bonds': {'1,6': {'monomer': ' 200H', 'bonds': {}}, '1,3': {'monomer': ' 200I', 'bonds': {}}}}, '1,3': {'monomer': ' 200D', 'bonds': {'1,2': {'monomer': ' 200E', 'bonds': {'1,2': {'monomer': ' 200F', 'bonds': {}}}}}}}}}}}}}}}
        ##b = {
        ##    'A146 ': {'bonds': {'1,N': {'monomer': ' 146A', 'bonds': {}}}},
        ##    'A86 ': {'bonds': {'1,N': {'monomer': ' 86A', 'bonds': {}}}},
        ##    'A200 ': {'bonds': {'1,N': {'monomer': ' 200A', 'bonds': {'1,4': {'monomer': ' 200B', 'bonds': {'1,4': {'monomer': ' 200C', 'bonds': {
        ##        '1,3': {'monomer': ' 200D', 'bonds': {'1,2': {'monomer': ' 200E', 'bonds': {'1,2': {'monomer': ' 200F', 'bonds': {}}}}}},
        ##        '1,6': {'monomer': ' 200G', 'bonds': {
        ##            '1,3': {'monomer': ' 200I', 'bonds': {}},
        ##            '1,6': {'monomer': ' 200H', 'bonds': {}},
        ##            }},
        ##        }}}}}}}}}
        ##print a == b

        d_molecules = {}
        for root in set_roots:
            path = [root]
            d_molecules[root] = {'bonds':{}}
            d_m_fwd = d_molecules[root]['bonds']
            l_monomers = [root]
            l_branches = []
            monomer1 = root
            while True:
##                print 'xxx', monomer1, d_molecules[root]

                ##
                ## end of branch and no branching points
                ##
                if monomer1 not in d_adjacency_forward.keys() and len(l_branches) == 0:
                    break

                ##
                ## end of branch but branching points
                ##
                elif monomer1 not in d_adjacency_forward.keys() and len(l_branches) > 0:

                    ## rewind path to previous branching point
                    path = path[:path.index(l_branches[-1])+1]
##                    print 'a0path', path, monomer1, l_branches
                    ## forward in dictionary to previous branching point
                    d_m_fwd = self.rewind_dictionary(root, path, d_molecules, d_adjacency_forward)

                    l_monomers.append(monomer1)

                    try:
                        monomers2 = list(set(d_adjacency_forward[path[-1]].keys())-set(l_monomers))
                    except:
                        print monomer1, d_adjacency_forward
                        print l_branches
                        print s_pdb
                        print path
                        stop

                    ## end of branch
                    if len(monomers2) == 0:
                        ## branch == root
                        if len(path) == 1 and len(l_branches) == 1 and path[0] == l_branches[0]:
                            break
                        elif len(path) == 1:
                            print path, l_branches
                            notexpected
                        ## branch != root
                        else:
                            ## remove branch point from list of branch points
                            l_branches.remove(path[-1])
                            ## rewind path to previous branching point
                            path = path[:-1]
                            monomer1 = path[-1]
##                            print 'a1path', path, monomer1, l_branches
                            ## forward in dictionary to previous branching point
                            d_m_fwd = self.rewind_dictionary(root, path, d_molecules, d_adjacency_forward)

                    else:
                        bond = d_adjacency_forward[path[-1]][monomers2[0]]
                        monomer1 = monomers2[0]
                        ## forward path
                        path.append(monomer1)
##                        print 'a2path', path, monomer1, l_branches
                        ## forward dictionary
                        d_m_fwd = d_m_fwd[bond]['bonds']

                ##
                ## monomer(s) forward of the current monomer (monomer1)
                ##
                else:

                    l_monomers.append(monomer1)

                    ## determine forward monomers that are in branches which have not been traveled
                    monomers2 = list(set(d_adjacency_forward[monomer1].keys())-set(l_monomers))

                    ## forward monomers traveled and no branches
                    if len(monomers2) == 0 and len(l_branches) == 0: ## e.g. 1h3t
                        break
                    ## forward monomers traveled and last branch point
                    elif len(monomers2) == 0 and len(l_branches) == 1 and [monomer1] == l_branches: ## e.g. 1h3u
                        break
                    ## jump backwards towards previous branch point
                    elif len(monomers2) == 0 and len(l_branches) > 0: ## e.g. 1h3u.pdb

                        ## rewind path to previous branching point
                        path = path[:path.index(l_branches[-1])+1]
                        monomer1 = path[-1]
                        ## remove branch point from list of branch points
                        l_branches.remove(path[-1])
##                        print 'a3path', path, monomer1, l_branches
                        ## forward in dictionary to previous branching point
                        d_m_fwd = self.rewind_dictionary(root, path, d_molecules, d_adjacency_forward)
                        continue

                    elif len(monomers2) == 0:

                        print 'l_branches', l_branches
                        print 'monomer1', monomer1
                        print 'monomers2', monomers2
                        print d_adjacency_forward
                        print d_molecules
                        print 'path', path
                        print s_pdb
                        notexpected2


                    ## translate monomers
                    for monomer2 in monomers2:
                        bond = d_adjacency_forward[monomer1][monomer2]
                        monomer = self.monomertranslation(monomer2,d_coordinates)
                        d_m_fwd[bond] = {'monomer':monomer,'bonds':{}}

                    if len(monomers2) > 1:
                        l_branches += [monomer1]

                    bond = d_adjacency_forward[monomer1][monomers2[0]]
                    monomer1 = monomers2[0]
                    
                    ## forward path
                    path.append(monomer1)
##                    print 'b0path', path, monomer1, l_branches
                    ## forward dictionary
                    d_m_fwd = d_m_fwd[bond]['bonds']

            root_hetID = self.monomertranslation(root,d_coordinates)
            d_molecules[root] = {root_hetID:d_molecules[root]}

##        print d_molecules

        return d_molecules


    def rewind_dictionary(self, root, path, d_molecules, d_adjacency_forward):

        ## forward in dictionary to previous branching point
        d_m_fwd = d_molecules[root]['bonds']
        for i in range(1,len(path)):
            m1 = path[i-1]
            m2 = path[i]
            bond = d_adjacency_forward[m1][m2]
            d_m_fwd = d_m_fwd[bond]['bonds']
        
        return d_m_fwd


    def monomertranslation(self,monomer,d_coordinates):

        chain = monomer[0]
        res_no = int(monomer[1:-1])
        iCode = monomer[-1]
        l_res_names = []
        for altloc in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
            res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']
            if res_name in self.d_saccharides.keys():
                res_name = self.d_saccharides[res_name]['stereo']
            l_res_names += [res_name]
        if len(set(l_res_names)) > 1:
            print l_res_names
            stop

        return res_name


    def parse_recordCONECT(self,line,d_atomnos,d_CONECT):

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
        element1 = d_atomnos[atom_no1]['element']

        ## skip if hydrogen atom
        if atom_name1[0] == 'H':
            return d_CONECT
        ## skip if a cofactor to which molecules are not *covalently* bound (use mmCIF to distinguish bond type instead..!)
        if res_name1 in self.l_cofactors+['HOH']+self.l_solutes:
            return d_CONECT
        ## skip if metal atom, because not covalent bond (use mmCIF to distinguish bond type instead..!)
        if element1 in self.l_atoms_metal: ## e.g. 1a6l
            return d_CONECT

        d_atom_nos = {'intra':[],'inter':[]}

        for atom_no2 in atom_nos[1:]:
            chain2 = d_atomnos[atom_no2]['chain']
            res_no2 = d_atomnos[atom_no2]['res_no']
            iCode2 = d_atomnos[atom_no2]['iCode']
            res_name2 = d_atomnos[atom_no2]['res_name']
            atom_name2 = d_atomnos[atom_no2]['atom_name']
            element2 = d_atomnos[atom_no2]['element']
            
            ## skip if hydrogen atom
            if atom_name2[0] == 'H':
                continue
            ## skip if a cofactor to which molecules are not *covalently* bound (use mmCIF to distinguish bond type instead..!)
            if res_name2 in self.l_cofactors+['HOH']+self.l_solutes:
                continue
            ## skip if metal atom, because not covalent bond (use mmCIF to distinguish bond type instead..!)
            if element2 in self.l_atoms_metal:
                continue
            if (
                chain1 !=  chain2 or
                res_no1 != res_no2 or
                iCode1 != iCode2
                ):
                d_atom_nos['inter'] += [atom_no2]
            else:
                d_atom_nos['intra'] += [atom_no2]


        ## skip if no *covalent* connections to non-hydrogen atoms
        ## do not ignore intra bonds (possibly disulfide bonds to other space groups)
        if len(d_atom_nos['inter']) == 0:
            return d_CONECT

        if res_name1 not in d_CONECT.keys():
            d_CONECT[res_name1] = {}
        if chain1 not in d_CONECT[res_name1].keys():
            d_CONECT[res_name1][chain1] = {}
        if res_no1 not in d_CONECT[res_name1][chain1].keys():
            d_CONECT[res_name1][chain1][res_no1] = {}
        if iCode1 not in d_CONECT[res_name1][chain1][res_no1].keys():
            d_CONECT[res_name1][chain1][res_no1][iCode1] = {}
        
        d_CONECT[res_name1][chain1][res_no1][iCode1][atom_no1] = d_atom_nos

        return d_CONECT


    def parse_recordREMARK(self, d_header, line, i, lines):

        remark = int(line[6:10])

        if remark == 200:

            experimentaldetail_key = line[12:23].strip().upper()
            if experimentaldetail_key in ['TEMPERATURE','PH']:
                experimentaldetail_value = line[44:].strip()
                if '(' in experimentaldetail_value and ')' in experimentaldetail_value: ## e.g. 1ar4
                    experimentaldetail_value = experimentaldetail_value[:experimentaldetail_value.index('(')]+experimentaldetail_value[experimentaldetail_value.index(')')+1:]
                    experimentaldetail_value = experimentaldetail_value.strip()
                    if experimentaldetail_value == '':
                        print experimentaldetail_value
                        stop
                if 'REMARK200' not in d_header.keys():
                    d_header['REMARK200'] = {}
                if experimentaldetail_key not in d_header['REMARK200'].keys():
                    d_header['REMARK200'][experimentaldetail_key] = experimentaldetail_value

        elif remark == 290:

            d_header = self.parse_recordREMARK290(d_header, i, lines)

        elif remark == 350:

            ## biological units
            ## (e.g. 2bq0.pdb, 1thj.pdb, 1m4x.pdb, 1d3i.pdb, 1qgc.pdb, 1rhi.pdb, 1rbo.pdb, 2g8g.pdb, 1h84.pdb)

            d_header = self.parse_recordREMARK350(d_header, i, lines)

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
                        if 'REMARK525' not in d_header.keys():
                            d_header['REMARK525'] = []
                        d_header['REMARK525'] += [solventchain]

        elif remark == 2: ## resolution
            if line[10:].strip() == '':
                return d_header
            if line.strip() == 'REMARK   2 RESOLUTION. NOT APPLICABLE.':
                return d_header
            if line[23:41] == '   NULL ANGSTROMS.':
                return d_header
            resolution = float(line[23:].replace('ANGSTROMS.',''))
            d_header['REMARK2'] = resolution

        elif remark == 0: ## rerefinement
            d_header['REMARK0'] = True

        elif remark == 465: ## missing residues (ATOM)
            d_header = parse_pdb.parse_recordREMARK465(line, lines, i, d_header, 465)

        elif remark == 470: ## missing atoms (ATOM)
            d_header = parse_pdb.parse_recordREMARK470(line, lines, i, d_header, 470)

        elif remark == 475: ## residues (ATOM) with zero occupancy
            d_header = parse_pdb.parse_recordREMARK465(line, lines, i, d_header, 475)

        elif remark == 480: ## polymer atoms (ATOM) with zero occupancy
            d_header = parse_pdb.parse_recordREMARK470(line, lines, i, d_header, 480)

        elif remark == 610: ## non-polymer residues (HETATM) with missing atoms
            pass

        elif remark == 615: ## non-polymer atoms (HETATM) with zero occupancy
            pass

        return d_header

    def parse_recordREMARK290(self, d_header, i, lines):

        line = lines[i]

        if 'REMARK290' not in d_header.keys():
            d_header['REMARK290'] = {}

        if line[13:18] == 'SMTRY':
            operator = int(line[21:23])
            if operator == 0:
                print line
                stop
            dimension = int(line[18])
            if dimension == 1:
                matrixrow1 = lines[i+0][24:].split()
                matrixrow2 = lines[i+1][24:].split()
                matrixrow3 = lines[i+2][24:].split()
                d_header['REMARK290'][operator] = [matrixrow1,matrixrow2,matrixrow3,]

        return d_header


    def parse_recordREMARK350(self, d_header, i, lines):

        line = lines[i]

        if 'REMARK350' not in d_header.keys():
            d_header['REMARK350'] = {}

        if line[11:23] == 'BIOMOLECULE:':
            biomolecules = line[23:80].replace(' ','').split(',') ## multiple biomolecules in e.g. 1wa3
            d_header = self.loop_and_identify_chains_and_matrices(i, lines, d_header, biomolecules)

        return d_header


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


    def loop_and_identify_chains_and_matrices(self, i, lines, d_header, biomolecules):

        chains = set()

        method = 'unknown'

        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
                break

            elif 'AUTHOR DETERMINED' in lines[j]:
                method = 'author'

            ## ignore software determined quarternary structures (parse directly from PISA later if needed) - problem with 5r1r,1r1r otherwise
            elif 'SOFTWARE DETERMINED' in lines[j]:
                if len(chains) > 0:
                    print lines[j]
                    stop
                if not (
                    'SOFTWARE USED: ' in lines[j+1]
                    or
                    'TOTAL BURIED SURFACE AREA:' in lines[j+1]
                    or
                    'APPLY THE FOLLOWING TO CHAINS:' in lines[j+1]
                    ):
                    print lines[j]
                    print lines[j+1]
                    stop
                if 'AUTHOR DETERMINED BIOLOGICAL UNIT:' in lines[j-1]:
                    method = 'combined'
                    continue
                elif 'AUTHOR PROVIDED BIOLOGICAL UNIT:' in lines[j-1]: ## 2vpl
                    method = 'combined'
                    continue
                else:
                    method = 'software'
                    continue
##                    if len(d_header['REMARK350'].keys()) == 0:
##                        ## do not break if only software determined (e.g. 2qkt,2qku; 1dze,1qm8)
##                        continue
##                    else:
##                        ## break if software determined and already author determined (e.g. 1y10,1y11?)
##                        if not 'BIOMOLECULE: ' in lines[j-1]:
##                            print lines[j-1]
##                            stop
##                        break

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

                if (
                    method == 'software'
                    and
                    int(lines[j][19:23]) == 1 ## only perform check for first matrix, otherwise break after first matrix if software determined (e.g. 1qkt)
                    ):
                    bool_break = False
                    for prev_bm in d_header['REMARK350'].keys():
                        ## one or more chains already part of another biomolecule?
                        if len(set(d_header['REMARK350'][prev_bm]['chains'].keys())&chains) > 0:
                            bool_break = True
                            break
                    if bool_break == True:
                        ## break if software determined and already author determined (e.g. 1y10,1y11?)
                        break ## break loop over lines
                    else:
                        ## do not break if only software determined (e.g. 2qkt,2qku; 1dze,1qm8)
                        pass

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
                    if biomolecule not in d_header['REMARK350'].keys():
                        d_header['REMARK350'][biomolecule] = {}

                    ## biomolecule > method
                    d_header['REMARK350'][biomolecule]['method'] = method

                    ## biomolecule > matrices
                    if 'matrices' not in d_header['REMARK350'][biomolecule].keys():
                        d_header['REMARK350'][biomolecule]['matrices'] = {}
                    ## matrices > matrixno > matrix
                    d_header['REMARK350'][biomolecule]['matrices'][matrixno] = matrixrows

                    ## biomolecule > chains
                    if 'chains' not in d_header['REMARK350'][biomolecule].keys():
                        d_header['REMARK350'][biomolecule]['chains'] = {}
                    for chain in chains:
                        ## chains > chain
                        if chain not in d_header['REMARK350'][biomolecule]['chains'].keys():
                            d_header['REMARK350'][biomolecule]['chains'][chain] = set()
                        d_header['REMARK350'][biomolecule]['chains'][chain] |= set([matrixno])

        return d_header


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
        import smallmolecules
        import quakes_init

        self.d_spacegroups, self.d_crystalsystems, self.d_crystalsynonyms = quakes_init.main()

        self.d_res1 = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            'MSE':'M','UNK':'X','ASX':'X','GLX':'X',
            }

        self.l_nucleotides = [
             'A', 'C', 'G', 'U',      'I', ## ribonucleotides
            'DA','DC','DG',     'DT','DI', ## deoxyribonucleotides
            'N', ## N is any 5'-monophosphate nucleotide
            ]
        
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

        self.d_dihedrals = {
            ## chi0
            'GLY':{},
            'ALA':{},
            ## chi1
            'VAL':{
##                'chi1a':['N','CA','CB','CG1',],'chi1b':['N','CA','CB','CG2',],
                'chi1':['N','CA','CB','CG1',],
                },
            'THR':{
##                'chi1a':['N','CA','CB','OG1',],'chi1b':['N','CA','CB','CG2',],
                'chi1':['N','CA','CB','CG2',],
                },
            'SER':{
                'chi1':['N','CA','CB','OG',],
                },
            'CYS':{
                'chi1':['N','CA','CB','SG',],
                },
            ## chi2
            'LEU':{
                'chi1':['N','CA','CB','CG',],
                'chi2a':['CA','CB','CG','CD1',],'chi2b':['CA','CB','CG','CD2',],
                },
            'ILE':{
##                'chi1a':['N','CA','CB','CG1',],'chi1b':['N','CA','CB','CG2',],
                'chi1':['N','CA','CB','CG1',],
                'chi2':['CA','CB','CG1','CD1',],
                },
            'ASP':{
                'chi1':['N','CA','CB','CG',],
                'chi2a':['CA','CB','CG','OD1',],'chi2b':['CA','CB','CG','OD2',],
                },
            'ASN':{
                'chi1':['N','CA','CB','CG',],
                'chi2a':['CA','CB','CG','OD1',],'chi2b':['CA','CB','CG','ND2',],
                },
            'HIS':{
                'chi1':['N','CA','CB','CG',],
                'chi2a':['CA','CB','CG','ND1',],'chi2b':['CA','CB','CG','CD2',],
                },
            'PRO':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','CD',],
                },
            'TYR':{
                'chi1':['N','CA','CB','CG',],
                'chi2a':['CA','CB','CG','CD1',],'chi2b':['CA','CB','CG','CD2',],
                },
            'PHE':{
                'chi1':['N','CA','CB','CG',],
                'chi2a':['CA','CB','CG','CD1',],'chi2b':['CA','CB','CG','CD2',],
                },
            'TRP':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','CD1',],
                },
            ## chi3
            'MET':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','SD',],
                'chi3':['CB','CG','SD','CE',],
                },
            'GLU':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','CD',],
                'chi3a':['CB','CG','CD','OE1',],'chi3b':['CB','CG','CD','OE2',],
                },
            'GLN':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','CD',],
                'chi3a':['CB','CG','CD','OE1',],'chi3b':['CB','CG','CD','NE2',],
                },
            ## chi4
            'LYS':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','CD',],
                'chi3':['CB','CG','CD','CE',],
                'chi4':['CG','CD','CE','NZ',],
                },
            ## chi4
            'ARG':{
                'chi1':['N','CA','CB','CG',],
                'chi2':['CA','CB','CG','CD',],
                'chi3':['CB','CG','CD','NE',],
                'chi4':['CG','CD','NE','CZ',],
                },
            }


        d = smallmolecules.main()

        ## connectivity of metal *atoms* are not checked (too many examples)
        self.l_atoms_metal = d['metals']

        self.l_clusters = d['clusters']
        self.l_prosthetic_groups = d['prosthetic groups']
        self.l_coenzymes = d['coenzymes']
        self.d_ions = d['ions']
        ## connectivity of cofactor *molecules* are not checked
        self.l_cofactors = self.l_clusters+self.l_prosthetic_groups+self.l_coenzymes+self.d_ions.keys()

        ## connectivity *and* identity of solutes *and* ions not checked
        self.l_solutes = d['solutes']

        self.d_saccharides = d['saccharides']
        self.d_stereoisomers = {
            }
        self.l_lpeptidelinkers = d['lpeptidelinkers']
        self.l_dpeptidelinkers = d['dpeptidelinkers']
        self.l_terminalmodres = d['terminalmodres']


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
            'title1','title2','hetIDs1', 'hetIDs2',
            ]

        self.d_rmsd_max = { ## 1r7r,3cf3; 1nlf,1olo; 2fsy,2ft1 ## try other combinations if above this rmsd (must be below 5.0 for 1hrd,1aup)
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

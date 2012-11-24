#!/software/bin/python
#
#$Id: quakes.py 266 2007-11-01 13:15:43Z tc $
#
#Tommy Carstensen, University College Dublin, 2007

## Log
## 2008May01 - Realized REMARK 290 transformations would have to be combined with the 27 types of translation (permutations with repetition of x,y,z and 0,+,-) for all combinations of chains in order to yield similar biounits. Some of these transformations would be identical to each other. Instead I chose to use the PISA transformations.

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

## iCode ' ' from ATOM records should be put before iCode 'A' from REMARK465 records by comparsion to the SEQRES records upon parsing the ATOM records

## write a new faster seq aln alg

## if no remark350 record and high rmsd then try monomers...

## don't do seq aln during rmsd2bfactor...
## don't do 1st SEQRESseq aln for seq id chains?!

## accept neutron diffraction structures? accept 3ins which is x-ray AND neutron?

## what to do with 1ft8:E and 1koh:B for which coordinates of aligned residues are not given???

## what to do if sugars are missing? they are not described in remark465 records. this can results in both false negatives and fakse positives...

## solve the problem of too many chain combinations by sequential pairing
## this will yield a maximum of n(n+1)/2 combinations to check and with fewer coordinates per rmsd calculation
## sequential addition of sequence similar chains to reduce number of combinations further

## make sure that glycosides are not only identical but also stemming from the same asn/ser/thr residue!

## plot of rmsd vs disulphide bonds (y/n)
## plot of rmsd vs EC class
## plot of rmsd vs CATH class

## do not use pdb if only CA atoms in ATOM records (make a set...)

## whatif accesibility should be calculated for the biou and not the asyu

## Sep26 - number of residues should be split in residues1 and residues2 and should include missing residues and mutated residues

## Oct12 - should NAD (NAD+) and NAI (NADH) be considered the same compound, since hydrogen atoms cant be seen in x-ray structures?
## Oct12 - should NAP (NADP+) and NDP (NADPH) be considered the same compound, since hydrogen atoms cant be seen in x-ray structures?
## Oct12 - should NAD (NAD+) and NAJ (NAD+ acidic form) be considered the same compound, since hydrogen atoms cant be seen in x-ray structures?

## 2008May02 - add info about sequence different peptide chains to html output

## 2008jul25 - what sequences to compare when two alternative res_names?

## incorrect d_res_nos will be returned from ATOM2seq for 2bfk (vs 2bfl), chain A becauseof ASN61A in REMARK465 records!!!

import os, sys, numpy, LinearAlgebra, math, copy
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import biounit,statistics,gnuplot
sys.path.append('/home/people/tc/svn/Protool/')
import geometry
sys.path.append('/home/people/tc/svn/EAT_DB/')
import sequence_alignment

class quakes:

    def main(self):

##        self.rsync()
##        self.gunzip()

##        import os
##        for s in 'abcdefghijklmnopqrstuvwxyz890':
##            print s
##            os.system('mkdir pdb/%s' %(s))
##        stop

        errorpdbs = [
            ## BAFF-BAFF-R complex
            '1p0t','1otz',
            ## TRAP
            '1utf','1utv',
            ## 30S subunits
            '1jgo','1jgp','1jgq',
            ## 30S,50S subunits of 70S ribosome
            '1gix','1giy',
            '1pns','1pnu','1pnx','1pny',
            '1s1h','1s1i',
            '1voq','1vor','1vos','1vou','1vov','1vow','1vox','1voy','1voz','1vp0',
            '1vs5','1vs6','1vs7','1vs8',
            '1vsa','2ow8',
            '2avy','2aw4','2aw7','2awb',
            '2b64','2b66','2b9m','2b9n','2b9o','2b9p',
            '2i2p','2i2t','2i2u','2i2v',
            '2j00','2j01','2j02','2j03',
            '2qou','2qov','2qow','2qox','2qoy','2qoz','2qp0','2qp1',
            '2v46','2v47','2v48','2v49',
            '2hgi','2hgj','2hgp','2hgq',
            ## fatty acid synthase
            '2uv9','2uva','2uvb','2uvc',

            ## identical ATOM/HETATM IDs (coordinate section)
            '1h9h','2olb','11gs','121p','12gs','13gs','16gs','16pk','17gs','185d','18gs','193d','19gs','1a05','1a0i','1a25','1a3l','1a44','1a48','1a4a','1a4b','1a4c','1a4f','1a4g','1a4k','1a4q','1a52','1a5a','1a5b','1a5z','1a65','1a69','1a6g','1a6m','1a6q','1a71','1a78','1a79','1a8i', '1a8s', '1a8u', '1a9c', '1a9m', '1a9x', '1a9y', '1a9z', '1aal', '1aax', '1aaz', '1aba', '1abw', '1aby', '1ad8', '1ad9', '1adb', '1adc', '1adf', '1adg', '1adj', '1aec', '1af6', '1afa', '1afb', '1ag0', '1ah8', '1ahe', '1ahf', '1ahx', '1ahy', '1ai4', '1ai5', '1ai6', '1ai7', '1ai8', '1aix', '1aj6', '1aj9', '1ajn', '1ajp', '1ajq', '1aku', '1akv', '1all', '1ami', '1ao5', '1aok', '1aor', '1apm', '1aq1', '1aq6', '1aq7', '1ar1', '1ari', '1art', '1arx', '1asm', '1asn', '1aso', '1asp', '1asq', '1at1', '1atg', '1atl', '1atn', '1atp', '1au1', '1aua', '1auj', '1aus', '1av4', '1av6', '1avb', '1avf', '1avq', '1awb', '1axa', '1axd', '1axs', '1aya', '1ayo', '1ayr', '1azr', '1azs', '1b08', '1b0m', '1b25', '1b4k', '1b4n', '1b55', '1b6h', '1b8t', '1b92', '1b9f', '1bb1', '1bch', '1bcj', '1bcr', '1bcs', '1bd0', '1be3', '1beh', '1ben', '1bfn', '1bgy', '1bh3', '1bij', '1biq', '1biz', '1bj3', '1bja', '1bk5', '1bks', '1bnl', '1bow', '1bps', '1bqa', '1bqd', '1bqh', '1brh', '1brk', '1brt', '1bsk', '1btc', '1bvc', '1bvd', '1bwn', '1bwo', '1cdk', '1cdm', '1cel', '1ces', '1cf5', '1cgk', '1ch0', '1ch4', '1cls', '1cml', '1cnt', '1coh', '1con', '1cow', '1cpc', '1cpq', '1cpr', '1crx', '1ctf', '1ctp', '1cxe', '1cxf', '1cxh', '1czf', '1d33', '1d35', '1d39', '1d40', '1d61', '1d8h', '1daj', '1dan', '1dbj', '1dbk', '1dbp', '1dbr', '1dff', '1dhy', '1dif', '1djc', '1dkk', '1dl4', '1dlr', '1dls', '1dmr', '1dms', '1dpe', '1dpm', '1drf', '1drj', '1drk', '1dut', '1eas', '1ecc', '1ecg', '1efg', '1eg6', '1elp', '1ent', '1epm', '1epn', '1epp', '1epq', '1eta', '1etb', '1eth', '1etj', '1fax', '1fbt', '1fdh', '1fgh', '1fgj', '1fig', '1fiv', '1fkb', '1fki', '1fn1', '1fpc', '1fpi', '1fpj', '1fpl', '1fsa', '1fui', '1fuo', '1fup', '1fuq', '1fur', '1fv7', '1fyl', '1g5l', '1gac', '1gaf', '1gan', '1gbu', '1gdi', '1gic', '1gj2', '1glc', '1gld', '1gle', '1gli', '1gmk', '1gmp', '1gmr', '1gnm', '1gnn', '1gno', '1gnw', '1gpa', '1gpd', '1gqh', '1gsu', '1gsy', '1gti', '1gtm', '1gtp', '1hab', '1hac', '1hah', '1hbg', '1hbi', '1hbs', '1hbt', '1hcn', '1hcs', '1hct',  '1hds', '1hdt', '1hfp', '1hfq', '1hfr', '1hga', '1hgb', '1hgc', '1hgt', '1hho', '1hih', '1hii', '1hiv', '1hmd', '1hmo', '1hpc', '1hr3', '1hro', '1hrp', '1htp', '1huj', '1huk', '1hur', '1hvq', '1hvr', '1hxp', '1hya', '1i3t', '1i47', '1i7v', '1iak', '1ibg', '1icj', '1idc', '1ide', '1iea', '1ieb', '1igt', '1igy', '1ils', '1ilu', '1imc', '1imd', '1ipw', '1irn', '1iro', '1ith', '1ixx', '1izb', '1j8l', '1jaw', '1jdb', '1jfd', '1jka', '1jkb', '1jkc', '1jkd', '1jlx', '1jpc', '1jsa', '1jst', '1jsu', '1kao', '1kpe', '1kr7', '1krb', '1krc', '1krn', '1ksa', '1kvq', '1kvr', '1kvs', '1kvt', '1kvu', '1lam', '1lcj', '1lck', '1lcp', '1ldp', '1len', '1lia', '1lkk', '1lkl', '1llo', '1lml', '1loa', '1lob', '1lpa', '1lpb', '1lth', '1ltr', '1lu2', '1lya', '1lyb', '1lyw', '1m5k', '1mae', '1maf', '1mbd', '1mfr', '1mhe', '1mhk', '1mmb', '1mmo', '1mmp', '1mmq', '1mmr', '1mpa', '1mpg', '1mpq', '1mpr', '1mwe', '1myf', '1myh', '1myj', '1nah', '1nai', '1nas', '1nbb', '1nbc', '1nbe', '1nci', '1nco', '1nec', '1nfd', '1nfp', '1nlk', '1nn2', '1nnc', '1np1', '1nrs', '1nsa', '1nsc', '1nsd', '1nzr', '1oaa', '1oat', '1ohj', '1olc', '1ord', '1otg', '1ouu', '1ova', '1ovw', '1p35', '1pag', '1pam', '1pbo', '1pea', '1pgt', '1php', '1ply', '1pnk', '1pnl', '1pnm', '1poi', '1ppm', '1ppr', '1psc', '1psh', '1ptv', '1pty', '1pyt', '1qbi', '1qd6', '1qdc', '1qf8', '1qge', '1qha', '1qhf', '1qpr', '1qrd', '1qs8', '1qsg', '1raa', '1rab', '1rac', '1rad', '1rae', '1raf', '1rag', '1rah', '1rai', '1rbl', '1rcm', '1rcp', '1rdg', '1rdi', '1rdj', '1rdk', '1rdl', '1rdm', '1rdn', '1rdv', '1rdx', '1rdy', '1rdz', '1rem', '1rge', '1rgf', '1rgg', '1rgh', '1rn1', '1rnc', '1rsy', '1rtf', '1rth', '1rti', '1rtj', '1rtm', '1rvw', '1sac', '1scu', '1sda', '1sdk', '1sdl', '1sep', '1sfi', '1sft', '1sgc', '1shb', '1shd', '1skj', '1sky', '1slb', '1slc', '1sli', '1slt', '1slu', '1smp', '1sps', '1stc', '1sty', '1swd', '1swe', '1swn', '1swp', '1swr', '1taq', '1tar', '1tat', '1taw', '1tcb', '1tcf', '1tco', '1tei', '1tet', '1thb', '1ths', '1tlc', '1tlg', '1tli', '1tlp', '1tmb', '1tmn', '1tmu', '1tn4', '1tpf', '1tpk', '1try', '1trz', '1tsd', '1ttp', '1ttq', '1tub', '1tyl', '1tym', '1tyu', '1tyw', '1tyx', '1udg', '1udh', '1uma', '1uz8', '1vdr', '1vrt', '1vru', '1vwt', '1wav', '1wgi', '1wgj', '1whs', '1wht', '1wyk', '1xel', '1xgm', '1xgn', '1xgs', '1xik', '1xso', '1yag', '1yec', '1yef', '1yeg', '1yeh', '1yrq', '1zeg', '1zeh', '1zni', '1znj', '20gs', '211d', '219d', '21gs', '221p', '22gs', '253d', '256b', '258d', '25c8', '277d', '2a3h', '2aaa', '2aac', '2aae', '2ahj', '2arc', '2at1', '2atc', '2ay1', '2ay2', '2ay3', '2ay4', '2ay5', '2ay6', '2ay7', '2ay8', '2ay9', '2azu', '2btf', '2bvw', '2c7e', '2cah', '2cht', '2cmm', '2cpk', '2cst', '2d34', '2d95', '2da8', '2dhb', '2dhn', '2dmr', '2dri', '2ecp', '2fbj', '2fus', '2gli', '2glr', '2gls', '2gsa', '2gss', '2hck', '2hhb', '2hhd', '2hk6', '2hmq', '2hmz', '2hr7', '2iep', '2kau', '2lal', '2lhb', '2mas', '2mcp', '2mhb', '2mpa', '2mpr', '2msb', '2nac', '2np1','2otc', '2ovw', '2oxi', '2pax', '2pcp', '2pfl', '2pgh', '2pgt', '2phk', '2pk4', '2pld', '2ple', '2prc', '2prg', '2q41', '2qwa', '2qwb', '2qwc', '2qwd', '2qwe', '2qwf', '2qwg', '2qwh', '2qwi', '2qwj', '2qwk', '2r2f', '2ran', '2rkm', '2rus', '2scp', '2sec', '2shk', '2sli', '2sn3', '2sni', '2sod', '2taa', '2tgd', '2tli', '2tmn', '2tn4', '2trc', '2tsc', '2wea', '2web', '2wec', '308d', '367d', '380d', '381d', '382d', '3a3h', '3at1', '3azu', '3chb', '3cyt', '3dmr', '3eng', '3fru', '3fx2', '3grs', '3gss', '3gst', '3hat', '3hhb', '3ljr', '3lkf', '3np1', '3p2p', '3pax', '3pca', '3pcb', '3pcc', '3pcd', '3pce', '3pcf', '3pcg', '3pch', '3pci', '3pcj', '3pck', '3pcn', '3pmg', '3prc', '3prn', '3pro', '3sc2', '3sli', '3sqc', '3tli', '3tmn', '421p', '455d', '456c', '4a3h', '4at1', '4cts', '4dfr', '4dmr', '4eng', '4fx2', '4gr1', '4gsa', '4gss', '4hhb', '4hmg', '4mdh', '4np1', '4ovw', '4pax', '4rub', '4sli', '4tgl', '4tmn', '5at1', '5bir', '5dnb', '5fx2', '5gss', '5hmg', '5mdh', '5prn', '5tmn', '621p', '6abp', '6adh', '6at1', '6gsp', '6gss', '6gsu', '6gsv', '6gsw', '6gsx', '6ins', '6prn', '6tmn', '721p', '7aat', '7abp', '7at1', '7gss', '7odc', '7prn', '7taa', '830c', '8aat', '8abp', '8at1', '8gss', '8prn', '966c', '9gss', '9pap', '9rub','2q44',
            ## identical ATOM/HETATM IDs (remark 465 section)
            '2nvq','2nvr','2nvt','2nvv','2nvx','2nvz','1ofr','1fah','1ksa','1qf8','2ay1','2ay2','2ay3','2ay4','2ay5','2ay6','2ay7','2ay8','2ay9','2bc2','2gsm','2hk6','2iep','2q41','3sqc','10gs','11gs','14gs','16gs','17gs','18gs','1a25','1a7k','1ady','1afs','1amu','1azs','1azt','1azx','1b2y','1bav','1bc2','1bd0','1beb','1bls','1bqh','1byf','1c09','1cdl','1cx2','1dao','1dbr','1dcm','1ddo','1dps','1ds6','1e0e','1e2y','1e6v','1e8a','1eh6','1eh7','1eh8','1f7n','1f7r','1f9c','1fag','1fat','1fpd','1fpe','1fpf','1fpg','1fpi','1fpj','1fpk','1fpl','1fpp','1frf','1fta','1fwq','1g20','1gcw','1ggy','1ggz','1glb','1gn3','1gnw','1gxf','1gyl','1hef','1heg','1hh8','1hkb','1hkc','1hro','1hrp','1hw7','1ids','1ieb','1ima','1imb','1imc','1imd','1ime','1j8t','1j8u','1jaf','1jk7','1jni','1jpz','1ju3','1jv4','1jw9','1jwb','1k9i','1ke1','1kf6','1kfy','1khv','1kif','1kw0','1kwp','1l0v','1lgc','1lld','1lm6','1lnq','1lq8','1lr0','1lrm','1lth','1mcx','1mlw','1mqf','1mu5','1muc','1mwo','1mz4','1n4p','1n4q','1n4r','1n4s','1n5k','1n5l','1n97','1nik','1nkx','1nlt','1nm0','1nqe','1nqg','1nqh','1nyk','1oas','1odb','1oqc','1oqj','1oqm','1ovx','1piy','1piz','1pj0','1pqv','1pxe','1pxx','1pyg','1q2l','1q55','1q5a','1q5b','1q5c','1q74','1qfk','1qj3','1qnl','1r33','1r5u','1r9s','1r9t','1rgq','1rly','1rqq','1sdd','1sfo','1si0','1si1','1sid','1sie','1sky','1smi','1smj','1sml','1t1g','1t1i','1tfp','1tfz','1tkw','1tlg','1tqn','1tqs','1tqt','1tqu','1tqv','1tqw','1tv2','1tv3','1tyl','1tym','1uaz','1uiv','1vyw','1w0e','1w0f','1w0g','1wiy','1x9r','1x9u','1xre','1xu3','1xuq','1xvd','1xvg','1y0j','1y67','1y6w','1ybv','1ymm','1yq4','1yro','1yu1','1ywh','1zo4','1zo9','1zoa','1zr6','1zt4','20gs','2a5k','2a8a','2a8h','2aa2','2aio','2alt','2alu','2ani','2ar3','2are','2atf','2axr','2b76','2b8k','2ccy','2cwg','2e2h','2e2i','2e2j','2eun','2euo','2eup','2euq','2eur','2eus','2eut','2euu','2exr','2fe6','2fer','2feu','2fta','2g73','2g74','2gmr','2gn7','2gnd','2gnm','2gss','2gsz','2gug','2h34','2h4g','2h4k','2h6a','2h7q','2h7r','2h7s','2h8h','2h98','2h9q','2hfz','2hg0','2hhm','2hkk','2hmj','2hml','2hmm','2hmn','2hmo','2hsi','2i8e','2i9z','2idb','2iga','2iid','2ij4','2iss','2min','2np5','2ns9','2nw8','2nw9','2ofp','2ogj','2ovw','2ox4','2oyc','2p0a','2p27','2p69','2p8e','2paj','2pms','2qqp','2r2f','2sqc','2yu9','3bls','3gss','3min','3pgh','3pgm','3sdh','3sdp','4blc','4cox','4gss','4sbv','4sdh','5cox','5gss','6cox','6gss','7cat','8cat','8ruc','9gss','2zb6','2rhb','2r8v','2rdz','1xmb','2zcg','2z2s','2rk3','2rk4','2qyf','2res','2rdb',
            ## identical ATOM/HETATM IDs (HET section)
            '1gah',
            ## residue renumbering (altloc, not diff res)
            '2nn8',
            ## residue renumbering (N-terminal residues)
            '1hto','1htq','1fkn','2hk2','2hk3','2qri','2qrs','1t41','1yr9','1yrb','1pbh','2pbh','1fng','1fne','2z5t',
            ## residue renumbering (internal insertion)
            '2ht8',##'1jrt','1jrs','1g2w','1mbq','2tbs',
            ## residue renumbering (reverse)
            '1zzn','1u6b',
##            ## duplicate MODRES records
##            '1ex1',
##            ## MODRES, multiple standard residue derivatives
##            '1fkp',
####            '1a14','1vjh','1zii','1zij','1xri','1o18','1tmq','2ja5','2ja6',
####            '2ja7','2ja8','1nc2','1u6p','1way','2bvr','2bxt','5gds','1bei',
####            '1c4e','1mqx','1mqy','1mqz','2h1o','1c30','1pu4','2c10',
##            ## misc
##            '1en2', ## identical ATOM/HETATM IDs (SER/GLY A 10 )
##            '1enm', ## identical ATOM/HETATM IDs (SER/GLY A 10 )
##            '1d8s','1i3q','1i50','1iw7', ## REMARK465 records missing
##            '1abj', ## REMARK465 record but residue not missing
##            ## TML vs 3ML
##            '1chh','1ccr','1chi','1chj','1cie','1cif','1cig','1chi','1crg',
##            '1crh','1cri','1crj','1csu','1csv','1csw','1csx','1cty','1ctz',
##            '1fhb','1rap','1raq','1ycc','1yea','1yeb','1ytc','2cln','2ycc',
            ## REMARK290 errors
            '1iph',
            ## REMARK350 errors
            '1pkh','1cr0','8atc',## '1a8r','1a9c','1qom','2scu','1onr','2pab','2v5l','1jfa',
            ## pisa entry not found
            '3bcm',

            ## should be dimers according to paper
            '1ow6','1ow7',
            ## large rmsd
            '1ygz','2bqx',
            

            ## incorrect atom names
            '2ric','2exj',
##             ## vs 1d7e
##            '1x83', ## vs '1x84' ## BR atom missing in hetID SBH

##            ## N-terminal posttranslational modification not in SEQRES records
##            '1vba','1al2','1ar6','1ar9','1asj',
##            ## altloc flag missing
            '1xpk','3xis','1ofz','1abe','1abf','5abp','1axz',

            ## duplicate atom position
            '2axm',

##            ##
##            ## tommy errors
##            ##
            '2gp1', ## virus yields too big an rmsd...
            '1ft8','1koh', ## too many missing coordinates (SEQRES overlap but no ATOM overlap)
##            ## alternate residues (alternate hetID/altloc)
##            '1fh2','2fi5','1eis',
##            ## alternate locations
##            '1uwn','1uwp',
##            '1k55','1k56', ## LYS/KCX
##            '2ftm', ## ASP/IAS
##            ## carbon atom (of small molecule) to backbone amide nitrogen
            '2h8p','2h5g', ## GOA linker
####            '1k56','1k57', ## post translational modification of selected chains
##            '2agx', ## TRQ er en MODRES, TSH er ikke. Derfor er raekkefoelgen givet. scan MODRES records for at faa info om TRQ...

            ## differently sized biounits
            '2qkt','2qku',

            ## cyclic peptide
            '1vwi','1vwj','1vwk','1vwl','1vwq','1vwr',

            ##
            ## connections
            ##

##            ## other unexpected connections
##            '1run',

            ## missing connections
            '1kx0',

            ## unexpected incorrect connections
            '2c03','1ljy','1fzg','1n2c','1r9m','2bod','1p8h','2fa7','2cfg','2bhy','1sbd','1sbe','1ppv','1x83','2v5s','3car','1ea5','2c3w','1ayn',

            ## unexpected connections (unknown if correct)
            '2nzy','1mpm','1fmu','4atj','1aym','1j3y','2ea0','1k7d',##'1e56',
            '1xwq', ## XYP O4B C5B bond
            '1h2z', ## SIN
            '1sbh','1sbi','1sbn','1yja','1yjb', ## CYS-HYD
            '1nhs', ## CYS-CYO
            '2uyv', ## GLY-TLA
            '2al2', ## LYS:NZ-ALA:N
            '3by4', ## 3CN
            '2bd2','2bd4','2bd7', ## casomorphin...
            '1run', ## CMP:O5'-SER:OG
            '2bcn', ## ILE:CD1-LEU:O
            '1uyx', ## GLN:OE1,BGC:O1
            '2pyp', ## ARG:NH2-HC4:C3'

##            '1k7d', ## SER-GRO
##            '2hlp', ## CONECT record for a standard peptide bond

            ## unexpected but correct connections
            '1b8f','1gkj','2nxy', ## GLY:N/ALA:C
            '1ggk','1gge','1ggj', ## TYR:CB-HIS:ND1
            '1mwv','2fxg','2b2o','2b2r','2dv1', ## TYR:CE1-TRP:CH2
            '1odw', ## DMN-HPH
            '1mtb','2fgu', ## DIQ:CM-HPH:C
            '1w3p', ## PYR:O1-HIS:NE2 hydrogen bond
            '1y21', ## NO (nitrogen oxide)
            '2azc', ## PHE:c-PHE:C
##            '1ca8','1bb0' ## vs  ## GLY-NVA
            ## small molecule connection (e.g. pyruvate, glycerol, hydrogen peroxide PEO)
##            '1qs7','1zbg','1fff','1fg8','1k2b','1vwb', ## vs 1qtx,1zj7,1ffi (NH2)
##            '222l','228l','220l','227l','1iii', ## BME

##            ## REMARK 465 record and/or SEQRES record, no MODRES record
##            '1gct','2iff','1x9s','1x9w','1xkf','1wu1','1tk5','1s10',
##            '1skw','1sal','1ry1','1pts',
##            ## REMARK 465 record *and* HETATM record
##            '1hao',
##            ## residue renumbering (not involving icodes? and therefore no errors)
##            '1t7b','1t7e','1tud','1yra','2qwa','2qwb','2qwc','2qwd','2qwe','2qwf','2qwg','2qwh','2qwi','2qwj','2qwk',

            ]

        if '-manual' in sys.argv:
    ######### manual selection

            self.l_pdbs = [sys.argv[-2],sys.argv[-1]]

            ## removal
            for pdb in errorpdbs:
                try:
                    self.l_pdbs.remove(pdb)
                except:
                    None

            self.pdbcount = len(self.l_pdbs)

            d_rmsd = self.analyze_sequences(verbose=True)
    ####        d_rmsd = read_rmsd_from_file()
    ##        self.write_rmsd_to_file(d_rmsd, d_header, prefix='rmsd')
            stop_manual
    ##########
        
############ out
        if '-out' in sys.argv:

            import os

            sets_pdbs = []
            self.l_pdbs = set()

            dirs = os.listdir('/oxygenase_local/tc/quakes/pdb/')
            for dir in dirs:
                files = os.listdir('/oxygenase_local/tc/quakes/pdb/%s' %(dir))
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
                self.pdbcount = len(self.l_pdbs)
                d_rmsd = self.analyze_sequences(verbose=False)
                self.write_rmsd_to_file(d_rmsd, d_header, prefix='rmsd')
##############
                
############## clusters 95
        if '-clusters' in sys.argv:

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
                self.cluster = cluster
                if cluster < int(sys.argv[-2]):
                    continue
                if cluster > int(sys.argv[-1]):
                    continue
                print '-------cluster--------', cluster, len(l_clusters)
                self.l_pdbs = list(set(clusters[cluster]))
                pdbs = self.l_pdbs

##                if '1bsr' not in pdbs:
##                    continue
    ##            print cluster, pdbs
    ##            stop

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
                d_rmsd = self.analyze_sequences()
                try:
                    d_rmsd = self.analyze_sequences()
                except:
                    print cluster
                    print self.pdb1, self.pdb2
                    stop
##########

        self.analyze_rmsd(conversion = True)

        return


    def analyze_rmsd(self, conversion = True):

        print 'analyze rmsd'

        import os

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

        terminal = 'postscript'

        n_lines = 1+len(self.l_columns_html)+1

        htmls = os.listdir('htm/')
        htmls.sort()
        lines = []
        d_pdbs_skip = {}
        for i in range(len(htmls)):
            html = htmls[i]
            pdb = html[:4]
            if i % 100 == 0:
                print '%4.1f%%' %(100*float(i)/len(htmls))
            fd = open('htm/%s' %(html),'r')
            lines = fd.readlines()[1+n_lines:-1]
            fd.close()

            ##
            ## parse data
            ##
            d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,set_pdbs = self.parse_html(
                lines,n_lines,d_pdbs_skip,
                d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal
                )
            
        print

        ## continuous ratio scale data
        for parameter in d_data_ratio_continuous.keys():
            for type in d_data_ratio_continuous[parameter]:
                gnuplotdata = ''
                for values in d_data_ratio_continuous[parameter][type]:
                    value = values[0]
                    rmsd = values[1]
                    gnuplotdata += '%s %s\n' %(value, rmsd)
                prefix_gnuplot = 'ps/%s%s' %(parameter,type)
                fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
                fd.write(gnuplotdata)
                fd.close()
                ## linear regression
                if type in ['difference','single',] or parameter == 'res':
                    regression = True
                else:
                    regression = False
                gnuplot.scatter_plot_2d(prefix_gnuplot, regression=regression, errorbars=True, terminal=terminal)
                if conversion == True:
                    print 'convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
                    os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
                    os.remove('%s.ps' %(prefix_gnuplot))

        ## discrete ratio scale data
        for parameter in d_data_ratio_discrete.keys():

            gnuplotdata = ''
            for value in d_data_ratio_discrete[parameter].keys():
                for rmsd in d_data_ratio_discrete[parameter][value]:
                    gnuplotdata += '%s %s\n' %(value, rmsd)

            prefix_gnuplot = 'ps/%s' %(parameter)
            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata)
            fd.close()

            if parameter in ['residues','coordinates']:
                logarithmic = True
            else:
                logarithmic = False
            gnuplot.scatter_plot_2d(prefix_gnuplot, regression=True, logarithmic=logarithmic, errorbars=True, terminal=terminal)
            if conversion == True:
                print 'convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
                os.remove('%s.ps' %(prefix_gnuplot))

        ####
        ## nominal scale data
        ####

        ##
        ## hetIDs, remarks, transformations
        ##
        d_nominal = {
            'hetIDs':['identical','different'],
            'REMARK465':['True','False'],
            'REMARK470':['True','False'],
            'transformations':['True','False'],
            }
        for parameter in d_nominal.keys():
            l1 = d_data_nominal[parameter][d_nominal[parameter][0]]
            l2 = d_data_nominal[parameter][d_nominal[parameter][1]]
##            gnuplotdata = ''
##            for rmsd in l1:
##                gnuplotdata += '1 %s\n' %(rmsd)
##            for rmsd in l2:
##                gnuplotdata += '2 %s\n' %(rmsd)
##            prefix_gnuplot = 'ps/%s' %(parameter)
##            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
##            fd.write(gnuplotdata)
##            fd.close()
##            d_xtics = {d_nominal[parameter][0]:1.,d_nominal[parameter][1]:2.}
##            gnuplot.scatter_plot_2d(prefix_gnuplot, d_xtics=d_xtics, terminal=terminal)
##            if conversion == True:
##                print 'convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
##                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
##                os.remove('%s.ps' %(prefix_gnuplot))

            print parameter
            statistics.twosamplettest(l1,l2)

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
                if type == 'identical':
                    l1 += d_data_nominal['spacegroup'][type][spacegroup]
                elif type == 'different':
                    l2 += d_data_nominal['spacegroup'][type][spacegroup]
            prefix_gnuplot = 'ps/spacegroups%s' %(type)
            gnuplotdata = d_gnuplotdata[type]
            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata)
            fd.close()
            gnuplot.scatter_plot_2d(prefix_gnuplot, d_xtics=d_spacegroups, xlabel='spacegroups', terminal=terminal)
            if conversion == True:
                print 'convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot)
                os.system('convert %s.ps %s.png' %(prefix_gnuplot,prefix_gnuplot))
                os.remove('%s.ps' %(prefix_gnuplot))

        ## statistics
        print prefix_gnuplot
        statistics.twosamplettest(l1,l2)

        return


    def parse_html(
        self,lines,n_lines,d_pdbs_skip,
        d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,
        ):

        ## loop over table rows
        for i in range(len(lines)/(n_lines)):
            ## parse pdbs
            pdb1_line = lines[i*n_lines+3]
            pdb2_line = lines[i*n_lines+4]
            pdbs = [pdb1_line,pdb2_line]
            for j in range(2):
                pdb_line = pdbs[j]
                try:
                    pdb_index2 = pdb_line.index('</a>')
                except:
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
##                    ## error
##                    '2bfk','2bfl',
##                    ## less than 50 residues due to mutations
##                    '1p7i','1p7j',
                ## less than 50 residues due to REMARK465 records
                '2yxq','2yxr','1bx7','1bx8',
                ## less than 50 residues due to non-overlapping residues
                '1ft8','1koh',
                ##
                '1znb','1zlx','1i3h','1j5o','1f33','1qvc',
                '1ypu','1uc7','1jr8','1g6l','1ixc','1xri',
                '1xfi','2q40','2q47','2qd0','1m70',
                ##
                '1f1g','1f18','1abs','1dxd','1dvb',
                '1dmm','1dmn','1dmq',
                '195l','196l','200l','197l','199l','198l',
                ## unknown whether error or not
                '1za1','1qtv','1khh','1cvm','1ez9','1ulx','2der','1x98',
                ## correct (cryogen/cryo-cooled)
                '2fbb','1j3y','1j41','1ajg','1ajh','1a6k','1c1g','2i16','2qxw','2j7n',
                ## incorrect (celsius)
                '1ade','1c0e',
                ## remediation change
                '1dos',
                ## high pH
                '3c90',
##                                '1fgn',
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
                d_parameters[parameter] = value

            ## skip if virus (might exclude other proteins...)
            if int(d_parameters['chains']) % 60 == 0:
                continue

            ## check resolution
            if len(set([pdb1,pdb2]) & set(['2c32','1ye1',])) == 0:
                if d_parameters['res1'] != 'N/A':
                    if float(d_parameters['res1']) > 5.:
                        print pdb1,pdb2
                        stop
                if d_parameters['res2'] != 'N/A':
                    if float(d_parameters['res2']) > 5.:
                        print pdb2,pdb1
                        stop

            if float(d_parameters['residues']) < 50:
                errorpdbs = set([
##                    ## error
##                    '2bfk','2bfl',
##                    ## less than 50 residues due to mutations
##                    '1p7i','1p7j',
                    ## less than 50 residues due to REMARK465 records
                    '2yxq','2yxr','1bx7','1bx8',
                    ## less than 50 residues due to non-overlapping residues
                    '1ft8','1koh',
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
            if float(rmsd) > 5.: ## temp!!! due to incorrect transformations...
                continue

            ######
            ## continuous ratio scale data
            ######
            for parameter in d_data_ratio_continuous.keys():
                values = []
                for no in ['1','2']:
                    value = d_parameters[parameter+no]
                    if value not in ['NULL','N/A']:
                        if ';' in value:
                            errorpdbs = [ ## temp!!!
                                '1znb','1zlx','1i3h','1j5o','1f33','1qvc',
                                '1ypu','1uc7','1jr8','1g6l','1ixc','1xri',
                                '1xfi','2q40','2q47','2qd0','1m70',
                                ]
                            if (no == '1' and pdb1 in errorpdbs) or (no == '2' and pdb2 in errorpdbs):
                                l_values = value.split(';')
                                value = 0
                                try:
                                    for s_value in l_values:
                                        value += float(s_value)
                                except:
                                    print pdb1, pdb2, s_value, l_values
                                    stop
                                value /= len(l_values)
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

                        try:
                            value = float(value)
                        except:
                            print value
                            print no, pdb1, pdb2

                        if (value < 50 or value > 303) and parameter == 'T':
                            errorpdbs = set([
                                '1f1g','1f18','1abs','1dxd','1dvb',
                                '1dmm','1dmn','1dmq',
                                '195l','196l','200l','197l','199l','198l',
                                ## unknown whether error or not
                                '1za1','1qtv','1khh','1cvm','1ez9','1ulx','2der','1x98',
                                ## correct (cryogen/cryo-cooled)
                                '2fbb','1j3y','1j41','1ajg','1ajh','1a6k','1c1g','2i16','2qxw','2j7n',
                                ## incorrect (celsius)
                                '1ade','1c0e',
                                ## remediation change
                                '1dos',
                                ])
                            if len(set([pdb1,pdb2])-errorpdbs) == 2:
                                print value, no
                                print pdb1, pdb2
                                stop2temperature ## temp!!!
                        if (value < 1 or value > 12) and value != 'NULL' and parameter == 'pH':
                            errorpdbs = set([
                                ## high pH
                                '3c90',
##                                '1fgn',
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

        return d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,d_pdbs_skip


    def gunzip(self):

        import os

        self.l_pdbs = []

        subdirs = os.listdir(self.path_pdb)
        subdirs.sort()
        for subdir in subdirs:
            files = os.listdir(self.path_pdb+subdir)
            for file in files:
                if file[-2:] == 'gz':
                    ## gunzip
                    if os.path.isfile('%s%s/%s' %(self.path_pdb,subdir,file[:-3])):
                        os.remove('%s%s/%s' %(self.path_pdb,subdir,file[:-3]))
                    os.system('gunzip %s%s/%s' %(self.path_pdb,subdir,file))
                self.l_pdbs += ['%s' %(file[3:7])]

        self.l_pdbs = list(set(self.l_pdbs))

        return


    def rsync(self):

        import os

        os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ %s' %(self.path_pdb))

        return


    def analyze_sequences(self, verbose=False):

        import os, time

        print 'analyzing sequences'

        ## do not exclude pdb2 from pdb1 loop if sequence similar to previous pdb1
        ## since pdb A == B, A != C, B == C
        ## but do not analyze sequence of pdb A,B and then pdb B,A since A == B and B == A are equivalent

        d_rmsd_identical = {}
        d_pdb = {}
        d_hetero = {}
        d_ATOMseq = {}
        d_header = {}
        d_chains_intrapdb_sequence_identical = {}

        ##
        ## loop 1 over pdbs
        ##
        t1 = time.clock()
        for i1 in range(self.pdbcount-1):

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

            ## reset dictionary of coordinates to save memory
            d_pdb = {}

            ##
            ## loop 2 over pdbs
            ##
            for i2 in range(i1+1,self.pdbcount):

                self.pdb2 = pdb2 = self.l_pdbs[i2]

                if not pdb2 in d_header.keys():
                    d_header[pdb2]  = self.parse_header(pdb2)
                    d_chains_intrapdb_sequence_identical[pdb2] = self.identify_identical_chains_from_sequence_intra(
                        d_header,pdb2,
                        )

                skippdb = self.pdbskip(d_header, pdb2)
                if skippdb == True:
                    continue

                ##
                ## parse coordinates
                ##
                if pdb1 not in d_pdb.keys():
                    d_pdb[pdb1], d_hetero[pdb1], d_ATOMseq[pdb1] = self.parse_coordinates(
                        pdb1, d_header[pdb1], verbose=verbose
                        )
                if pdb2 in d_pdb.keys():
                    stop
                d_pdb[pdb2], d_hetero[pdb2], d_ATOMseq[pdb2] = self.parse_coordinates(
                    pdb2, d_header[pdb2], verbose=verbose
                    )

                ## identify biomolecule(s)
                d_biomolecules2 = self.identify_biomolecule(pdb2, d_header)

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

                        ## skip if different hetero compounds (1st check of hetIDs)
                        different = self.different_hetero_compounds(pdb1,pdb2,d_header)
                        if different == True:
                            print 'different hetIDs', pdb1, pdb2
                            continue

                        ## skip if different polymers (other than long peptides)
                        different = self.different_nucleotides_saccharides_shortpeptides(pdb1,pdb2,d_header)
                        if different == True:
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
                        ##
                        (
                            bmSEQRESchains1_not_similar_to_SEQRESchains2,
                            bmSEQRESchains2_not_similar_to_SEQRESchains1
                            ) = self.identify_chains_interpdb_not_sequence_similar(
                                pdb1, pdb2, bmchains1, bmchains2,
                                d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
                                d_header,
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
                                        if d_header[pdb]['SEQRES'][chain]['type'] == 'peptide':
                                            peptide = True
                                        ## check if long chain
                                        if len(d_header[pdb]['SEQRES'][chain]['seq']) > self.min_len_chain:
                                            long = True
                                        if peptide == True and long == True:
                                            print pdb, chain
                                            break
                                    if peptide == True and long == True:
                                        break
                            ## continue if long peptide
                            if peptide == True and long == True:
                                continue

                        ## skip if only alpha carbon atoms (e.g. 1thi)
                        alpha1 = self.identify_atom_types_in_long_peptide_chains(d_header[pdb1], d_pdb[pdb1])
                        if alpha1 == True: 
                            print 'alpha', pdb1
                            continue
                        alpha2 = self.identify_atom_types_in_long_peptide_chains(d_header[pdb2], d_pdb[pdb2])
                        if alpha2 == True:
                            print 'alpha', pdb2
                            continue

                        ##
                        ## skip if different hetero compounds (2nd check of connectivity)
                        ##
                        different = self.compare_hetero_compounds(
                            d_hetero,d_header,d_pdb,
                            pdb1,pdb2,bmchains1,bmchains2,
                            d_chains_intrapdb_sequence_identical,
                            )
                        if different == True: ## e.g. 1vbo,1vbp
                            print 'different connectivity', pdb1, pdb2
                            if verbose == True:
                                for root in d_hetero[pdb1].keys():
                                    print pdb1, root, d_hetero[pdb1][root]
                                for root in d_hetero[pdb2].keys():
                                    print pdb2, root, d_hetero[pdb2][root]
                            continue


                        ##
                        ## identify equivalent chains (interpdb) from structure
                        ##

                        ## identify equivalent chains and calculate rmsd
                        (
                            l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations,
                            ) = self.identify_interpdb_equivalent_chains_from_structure(
                                pdb1, pdb2,
                                d_chains_intrapdb_sequence_identical,
                                d_chains_interpdb_sequence_similar,
                                d_pdb, d_header,
                                biomolecule1, biomolecule2,
                                d_biomolecules1, d_biomolecules2,
                                d_ATOMseq,
                                )

                        if rm == None:
                            fd = open('transform_error.txt','a')
                            fd.write('%s %s\n' %(pdb1,pdb2))
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
##                                tchains1,tchains2,d_pdb,pdb1,pdb2,d_header,
##                                d_chains_interpdb_sequence_similar,
##                                d_chains_intrapdb_sequence_identical,
##                                l_atoms = l_atoms
##                                )

                        ##
                        ## identify number and location of mutations
                        ##
                        n_mutations,d_mutations = self.identify_mutations(
                            pdb1, pdb2, l_equivalent_chains,
                            d_chains_interpdb_sequence_similar, d_chains_intrapdb_sequence_identical,
                            )

                        if n_mutations == 1:
                            self.phipsi(l_equivalent_chains,d_pdb,pdb1,pdb2,d_mutations,d_chains_intrapdb_sequence_identical,d_ATOMseq)

##                        ## calculate rmsd in spheres around the mutation site
##                        if mutations == 1 and n_chains == 1:
##                            rmsd4, rmsd8, rmsd16, rmsd32 = self.spherermsd(
##                                pdb1, pdb2,
##                                d_header, d_pdb,
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
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['residues'] = n_residues+n_mutations
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['coordinates'] = n_coordinates
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['l_equivalent_chains'] = l_equivalent_chains
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['mutations'] = n_mutations
                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['transformations'] = transformations
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd4'] = rmsd4
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd8'] = rmsd8
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd16'] = rmsd16
##                        d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]['rmsd32'] = rmsd32

##                        ## calculate solvent accesible surface area
##                        d_pdb[pdb1] = self.whatif(pdb1,biomolecule1,d_pdb[pdb1])

                        ## write data to file
                        d_quickrmsd = {pdb1:{biomolecule1:{pdb2:{biomolecule2:d_rmsd_identical[pdb1][biomolecule1][pdb2][biomolecule2]}}}}
                        self.write_rmsd_to_file(d_quickrmsd, d_header, prefix='quickrmsd')

                        d_biomolecules = {
                            pdb1:{'biomolecule':biomolecule1},
                            pdb2:{'biomolecule':biomolecule2},
                            }

                        ## color code structure by rmsd
                        self.rmsd2bfactor(
                            pdb1, pdb2, biomolecule1, biomolecule2, rmsd,
                            d_pdb, d_header, tv1, rm, tv2, l_equivalent_chains,
                            bmchains1, bmchains2,
                            d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
                            d_ATOMseq,
                            )

                d_pdb[pdb2]

        return d_rmsd_identical


    def phipsi(self,l_equivalent_chains,d_pdb,pdb1,pdb2,d_mutations,d_chains_intrapdb_sequence_identical,d_ATOMseq):

        d_loop = {
            pdb1:{'chains':l_equivalent_chains[0],},
            pdb2:{'chains':l_equivalent_chains[1],},
            }
        for pdb in d_loop:
            ## get mutated chain
            rep_chain,mutation = d_mutations[pdb].items()[0]
            chains = [rep_chain]+d_chains_intrapdb_sequence_identical[pdb][rep_chain]
            chain = list(set(d_loop[pdb]['chains']) & set(chains))[0]
            if len(list(set(d_loop[pdb]['chains']) & set(chains))) != 1:
                print set(d_loop[pdb]['chains']), set(chains)
                print d_mutations
                stop_not_expected
            d_loop[pdb]['chain'] = chain
            ## get residue index of mutation
            res_index = mutation[0][0]
            d_loop[pdb]['res_index'] = res_index
            print pdb,mutation
        chain1 = d_loop[pdb1]['chain']
        chain2 = d_loop[pdb2]['chain']

        for pdb in d_loop:
            chain = d_loop[pdb]['chain']
            res_no = d_ATOMseq[pdb][chain]['res_nos'][res_index]
            iCode = d_ATOMseq[pdb][chain]['iCodes'][res_index]
            record = d_ATOMseq[pdb][chain]['records'][res_index]

            ## get coordinates
            if record == 'REMARK465':
                phi = '   N/A'
                psi = '   N/A'
            elif (
                'REMARK' in d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['N'].keys() or
                'REMARK' in d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA'].keys() or
                'REMARK' in d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['C'].keys()
                ):
                phi = '   N/A'
                psi = '   N/A'
            else:
                print d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']
                N = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['N']['coordinate']
                CA = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA']['coordinate']
                C = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['C']['coordinate']
                ## calculate dihedrals
                if res_index-1 in d_ATOMseq[pdb][chain]['res_nos']:
                    res_no_prev = d_ATOMseq[pdb][chain]['res_nos'][res_index-1]
                    iCode_prev = d_ATOMseq[pdb][chain]['iCodes'][res_index-1]
                    if d_pdb[pdb]['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['record'] == 'REMARK465':
                        phi = '   N/A'
                    else:
                        C_prev = d_pdb[pdb]['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['C']['coordinate']
                        phi = '%6.1f' %(self.dihedral(C_prev,N,CA,C,))
                        dist = math.sqrt(sum((C_prev-N)**2))
                        if dist > 2.:
                            print dist
                            stop_dist
                else:
                    phi = '   N/A'
                if res_index+1 in d_ATOMseq[pdb][chain]['res_nos']:
                    res_no_next = d_ATOMseq[pdb][chain]['res_nos'][res_index+1]
                    iCode_next = d_ATOMseq[pdb][chain]['iCodes'][res_index+1]
                    if d_pdb[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['record'] == 'REMARK465':
                        psi = '   N/A'
                    else:
                        N_next = d_pdb[pdb]['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['N']['coordinate']
                        psi = '%6.1f' %(self.dihedral(N,CA,C,N_next))
                        dist = math.sqrt(sum((C-N_next)**2))
                        if dist > 2.:
                            print dist
                            stop_dist
                else:
                    psi = '   N/A'
            d_loop[pdb]['phi'] = phi
            d_loop[pdb]['psi'] = psi
            d_loop[pdb]['res_no'] = res_no
            print pdb,'phi,psi',phi,psi
        res_name1 = d_mutations[pdb1].values()[0][0][1]
        res_name2 = d_mutations[pdb2].values()[0][0][2]
        fd = open('phipsi.txt','a')
        fd.write('%4s %4s %1s %1s %4i %4i %1s %1s %6s %6s %6s %6s\n' %(
            pdb1,pdb2,
            chain1,chain2,
            d_loop[pdb1]['res_no'],d_loop[pdb2]['res_no'],
            res_name1,res_name2,
            d_loop[pdb1]['phi'],d_loop[pdb1]['psi'],
            d_loop[pdb2]['phi'],d_loop[pdb2]['psi'],
            ))
        fd.close()

        return


    def dihedral(self,c1,c2,c3,c4):

        import numpy,math

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

        import numpy

        n = numpy.array([
            v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0],
            ])

        return n


    def compare_hetero_compounds(
        self,d_hetero,d_header,d_pdb,
        pdb1,pdb2,bmchains1,bmchains2,
        d_chains_intrapdb_sequence_identical,
        ):

        different = False
        pdbpairs = [[pdb1,pdb2],[pdb2,pdb1]]
        for pdbpair in pdbpairs:
            pdba = pdbpair[0]
            pdbb = pdbpair[1]
            for roota in d_hetero[pdba].keys():

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
                for rootb in d_hetero[pdbb].keys():

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
                            d_pdb[pdba]['chains'][roota[0]]['residues'][int(roota[1:])]['d_iCodes'][roota[-1]]['record'] == 'HETATM' and
                            d_pdb[pdbb]['chains'][rootb[0]]['residues'][int(rootb[1:])]['d_iCodes'][rootb[-1]]['record'] == 'HETATM'
                            ):
                            break

                        print d_header[pdba]['chains'][roota[0]]
                        print d_header[pdbb]['chains'][rootb[0]]
                        rootsa = d_hetero[pdb1].keys()
                        rootsb = d_hetero[pdb2].keys()
                        rootsa.sort()
                        rootsb.sort()
                        print rootsa
                        print rootsb
                        print pdba, roota, compounda
                        print pdbb, rootb, compoundb
                        print bmchains1
                        print bmchains2
                        print '*************'
                        dontgetdresnostwice
                        d_res_nosa, d_res_nosb, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(d_pdb, d_header, pdba, pdbb, roota[0], rootb[0])
                        for k,v in d_res_nosa.items():
                            if v['res_no'] == int(roota[1:]):
                                ka = k
                                break
                        for k,v in d_res_nosb.items():
                            if v['res_no'] == int(rootb[1:]):
                                kb = k
                                break
                        if ka+l1 != kb+l2: ## e.g. 2jbj,2pvw;1z4v,1z4w,1z4z;1d6q,1re2;1o8a,2iux
                            print ka,l1
                            print kb,l2
                            print pdba,roota,rootsa,d_res_nosa[ka],d_header[pdba]['chains'][roota[0]]['seq'][ka]
                            print pdbb,rootb,rootsb,d_res_nosb[kb],d_header[pdbb]['chains'][rootb[0]]['seq'][kb]
                            print self.cluster, 'cluster'
                            identical = False
                            continue
                        if d_header[pdba]['chains'][roota[0]]['seq'][ka] != d_header[pdbb]['chains'][rootb[0]]['seq'][kb]:
                            stop_expected3

                        ## if not different, then identical
                        identical = True
                        break
                    
                if identical == False:
                    different = True
                    break
            if different == True:
                different = True
                fd = open('CONECT.txt','a')
                fd.write('%s %s\n%s\n%s\n' %(pdb1, pdb2, d_header[pdb1]['TITLE'], d_header[pdb2]['TITLE']))
                fd.close()
                break

        return different


    def whatif(self,pdb,biomolecule,d_pdb):

        source = '%s/whatif/%s.txt' %(self.path_cwd,pdb1)
        if not os.path.isfile(source):
            fd = open(source,'w')
            fd.writelines([
                '/software/whatif_debugged/DO_WHATIF.COM <<EOF\n',
                'GETMOL %s%s/pdb%s.ent\n' %(self.path_pdb, pdb1[1:3], pdb1),
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
                                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][' ']['atoms'][atom_name]['acc'] = acc
                        else:
                            break
                    break

        os.remove('ALTERR.LOG')
        os.remove('PDBFILE')
        os.remove('WHATIF.FIG')

        return d_pdb


    def identify_atom_types_in_long_peptide_chains(self, d_header, d_pdb):

        atoms = set()
        for chain in d_header['SEQRES'].keys():
            if d_header['SEQRES'][chain]['type'] != 'peptide' or len(d_header['SEQRES'][chain]['seq']) < self.min_len_chain:
                continue
            for res_no in d_pdb['chains'][chain]['residues'].keys():
                for iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'].keys():
                    if 'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                        continue
                    for atom in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        if 'REMARK' not in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom].keys():
                            atoms |= set([atom])
                    if len(atoms-set(['CA'])) > 0:
                        return False
        if atoms != set(['CA']):
            print atoms
            print d_header['SEQRES'][chain]['type']
            notexected
        return True


    def identify_mutations(self, pdb1, pdb2, l_equivalent_chains, d_chains_interpdb_sequence_similar, d_chains_intrapdb_sequence_identical):

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

        return n_mutations, d_mutations


    def chain2repchain(self,pdb,chain,d_chains_intrapdb_sequence_identical,):

        for rep_chain in d_chains_intrapdb_sequence_identical[pdb]:
            if chain == rep_chain or chain in d_chains_intrapdb_sequence_identical[pdb][rep_chain]:
                return rep_chain

        return


    def different_nucleotides_saccharides_shortpeptides(self,pdb1,pdb2,d_header):

##        print 'identifying different polymers (not long peptides)'

        d_chains = {
            pdb1:{'peptide':set(),'saccharide':set(),'nucleotide':set()},
            pdb2:{'peptide':set(),'saccharide':set(),'nucleotide':set()},
            }
        for pdb in d_chains.keys():
            for chain in d_header[pdb]['SEQRES']:
                chain_type = d_header[pdb]['SEQRES'][chain]['type']
                chain_seq = d_header[pdb]['SEQRES'][chain]['seq']
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
                        d_hetIDs[pdb] |= set([d_header[pdb]['HET'][chain][res_no][iCode]])
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


    def pdbskip(self, d_header, pdb):

        pdbskip = False
        SEQRESchains = d_header[pdb]['SEQRES'].keys()

        ## continue if superseded structure
        if 'SPRSDE' in d_header[pdb].keys():
            pdbskip = True
            return pdbskip

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

##        ## continue if low resolution (e.g. 1jgp.pdb)
##        if d_header[pdb]['REMARK2'] != 'N/A':
##            if d_header[pdb]['REMARK2'] > self.minres:
##                pdbskip = True
##                return pdbskip

        return pdbskip


    def rmsd2bfactor(
        self, pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_pdb, d_header,
        tv1, rm, tv2,
        l_equivalent_chains, bmchains1, bmchains2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_ATOMseq,
        ):

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
        ## parse coordinates for the specific combinations of chain1 and chain2
        ##
        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            rep_chain1 = self.chain2repchain(pdb1,chain1[0],d_chains_intrapdb_sequence_identical,)
            rep_chain2 = self.chain2repchain(pdb2,chain2[0],d_chains_intrapdb_sequence_identical,)

            l1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1']
            l2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2']
            r1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r1']
            r2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r2']

            ## symmetry operation
            ## alignment transformation
            (
                coordinates1, coordinates2, rescount, lines1, lines2
                ) = self.ATOMrecords2coordinates2(
                    d_pdb, pdb1, pdb2, chain1[0], chain2[0],
                    l1, l2, r1, r2, d_ATOMseq,
                    rmsd = rmsd, tv1=tv1, rm=rm, tv2=tv2,
                    )

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

        
    def write_rmsd_to_file(self, d_rmsd, d_header, prefix):

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
            'chains':'chains (biounit)', 'residues':'residues (biounit)', 'coordinates':'coords (biounit)',
            'pH1':'pH1', 'pH2':'pH2', 'T1':'T1', 'T2':'T2',
            'res1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">resolution1</a>',
            'res2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">resolution2</a>',
            'spacegroup1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">spacegroup1</a>',
            'spacegroup2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">spacegroup2</a>',
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
                        hetIDs1 |= set([d_header[pdb1]['HET'][chain][res_no][iCode]])
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
                                hetIDs2 |= set([d_header[pdb2]['HET'][chain][res_no][iCode]])
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

        ## write html to local file
## overwrite if already that pdb comb and bm comb!!!
        if prefix == 'quickrmsd':
            path = 'htm/'
            for pdb in d_html.keys():
                file = '%s.htm' %(pdb)
                l_tr = d_html[pdb]
                self.append_table_rows(path,file,l_tr,th,d_rmsd)

        return


    def append_table_rows(self,path, file,l_tr,th,d_rmsd):

        import os

        print file

        if os.path.isfile('%s%s' %(path,file)):
            fd = open('%s%s' %(path, file),'r')
            lines_prev1 = fd.readlines()[1:-1]
            fd.close()
            n_lines = 1+len(self.l_columns_html)+1
            lines_prev2 = []

            for i in range(len(lines_prev1)/(n_lines)-1,-1+1,-1):
                ## parse pdbs
                pdb_lines = lines_prev1[i*n_lines+3:i*n_lines+5]
                bm_lines = lines_prev1[i*n_lines+5:i*n_lines+7]
                l_pdbs = []
                l_bms = []
                for j in range(2):
                    pdb_line = pdb_lines[j]
                    bm_line = bm_lines[j]
                    pdb_index2 = pdb_line.index('</a>')
                    pdb_index1 = pdb_line[:pdb_index2].rindex('>')+1
                    pdb = pdb_line[pdb_index1:pdb_index2]
                    bm_index2 = bm_line.index('</td>')
                    bm_index1 = bm_line[:bm_index2].rindex('>')+1
                    bm = bm_line[bm_index1:bm_index2]
                    l_pdbs += [pdb]
                    l_bms += [bm]
                pdb1 = l_pdbs[0]
                pdb2 = l_pdbs[1]
                bm1 = int(l_bms[0])
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
                    lines_prev2 += lines_prev1[i*n_lines:(i+1)*n_lines]

            html = []
            ## table init
            html += ['<table border="1">\n']
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
     

    def identify_biomolecule(self, pdb, d_header):

        SEQRESchains = d_header[pdb]['SEQRES'].keys()
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
            d_biomolecules = {
                1:{
                    'chains':SEQRESchains,
                    'polymercount':len(SEQRESchains),
                    }
                }

        return d_biomolecules


    def identify_chains_interpdb_not_sequence_similar(
        self, pdb1, pdb2, bmchains1, bmchains2,
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

        SEQRESchains1_not_similar_to_SEQRESchains2 = set(d_header[pdb1]['SEQRES'].keys()) & (set(bmchains1) - SEQRESchains1_similar_to_SEQRESchains2)
        SEQRESchains2_not_similar_to_SEQRESchains1 = set(d_header[pdb2]['SEQRES'].keys()) & (set(bmchains2) - SEQRESchains2_similar_to_SEQRESchains1)

        return SEQRESchains1_not_similar_to_SEQRESchains2, SEQRESchains2_not_similar_to_SEQRESchains1


    def calculate_rmsd_for_multiple_chains(
        self,chains1,chains2,d_pdb,pdb1,pdb2,d_header,
        d_chains_interpdb_sequence_similar,
        d_chains_intrapdb_sequence_identical,
        d_ATOMseq,
        verbose=True,
        l_atoms = [],
        ):

##        print 'calculating rmsd for chains %s, %s of %s, %s' %(chains1, chains2, pdb1, pdb2)

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
        ## parse coordinates for the specific combinations of chain1 and chain2
        ##
        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            rep_chain1 = self.chain2repchain(pdb1,chain1[0],d_chains_intrapdb_sequence_identical,)
            rep_chain2 = self.chain2repchain(pdb2,chain2[0],d_chains_intrapdb_sequence_identical,)

            l1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1']
            l2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2']
            r1 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r1']
            r2 = d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['r2']

            (
                coords1, coords2, rescount
                ) = self.ATOMrecords2coordinates2(
                    d_pdb, pdb1, pdb2, chain1[0], chain2[0],
                    l1, l2, r1, r2, d_ATOMseq,
                    l_atoms=l_atoms,
                    )

            ## append coordinates
            coordinates1 += coords1
            coordinates2 += coords2
            residue_count += rescount

        if len(coordinates1) > 500000:
            print 'structural alignment start'
        rmsd = instance_geometry.superpose(coordinates1,coordinates2)
        if rmsd == 0:
            print pdb1, pdb2
            notexpected_structuresidentical
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter
        
        if verbose == True:
            if len(chains1) < 60 and len(chains2) < 60:
                print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s, %s, %s' %(rmsd, len(chains1), residue_count, len(coordinates1), chains1, chains2, pdb1, pdb2)
            else:
                print 'rmsd=%4.1f, %2i chains, %4i residues, %5i coordinates, %s, %s' %(rmsd, len(chains1), residue_count, len(coordinates1), pdb1, pdb2)

        return rmsd, len(chains1), residue_count, len(coordinates1), tv1, rm, tv2


    def ATOMrecords2coordinates2(self,
        d_pdb, pdb1, pdb2, chain1, chain2,
        l1, l2, r1, r2, d_ATOMseq,
        ## transform coordinates
        rmsd=None, tv1=None, rm=None, tv2=None,
        ## parse coordinates of selected atom types
        l_atoms=[],
        ):

        rescount = 0
        coordinates1 = []
        coordinates2 = []
        lines1 = []
        lines2 = []

        ## determine i range
        if l1 != 0 and l2 != 0:
            stop_not_expected
        if r1 != 0 and r2 != 0:
            stop_not_expected
        if l1 == 0 and l2 == 0:
            i_range = range(len(d_ATOMseq[pdb1][chain1]['res_nos'])-abs(r2-r1))
        elif l2 == 0 and l1 > 0:
            i_range = range(len(d_ATOMseq[pdb1][chain1]['res_nos'])-abs(r2-r1))
        elif l1 == 0 and l2 > 0:
            i_range = range(len(d_ATOMseq[pdb2][chain2]['res_nos'])-abs(r2-r1))
        else:
            stop_not_expected

        ## loop over i range
        for i in i_range:
            i1 = i+l2
            i2 = i+l1
            res_no1 = d_ATOMseq[pdb1][chain1]['res_nos'][i1]
            res_no2 = d_ATOMseq[pdb2][chain2]['res_nos'][i2]
            iCode1 = d_ATOMseq[pdb1][chain1]['iCodes'][i1]
            iCode2 = d_ATOMseq[pdb2][chain2]['iCodes'][i2]
            if d_ATOMseq[pdb1][chain1]['records'][i1] == 'REMARK465':
                continue
            if d_ATOMseq[pdb2][chain2]['records'][i2] == 'REMARK465':
                continue

            res_name1 = d_pdb[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['res_name']
            res_name2 = d_pdb[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['res_name']
                
            if res_name1 in self.d_modres.keys():
                res_name1 = self.d_modres[res_name1]
            if res_name2 in self.d_modres.keys():
                res_name2 = self.d_modres[res_name2]

            mutation = False
            if res_name1 != res_name2:
                mutation = True
                if res_name1 not in self.d_res.keys() or res_name2 not in self.d_res.keys():
                    stop_temporary

            d_atoms1 = copy.deepcopy(d_pdb[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'])
            d_atoms2 = copy.deepcopy(d_pdb[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'])

            rescount += 1

            ##
            ## transformation of coordinates
            ##
            if rmsd:
                SS = []
                for atom_name in d_atoms1.keys():
                    if 'REMARK' in d_atoms1[atom_name].keys():
                        continue
                    coordinate1 = d_atoms1[atom_name]['coordinate']
                    coordinates1 += [coordinate1]
                for atom_name in d_atoms2.keys():
                    if 'REMARK' in d_atoms2[atom_name].keys():
                        continue
                    coordinate2 = d_atoms2[atom_name]['coordinate']
                    coordinate2 = numpy.dot(coordinate2-tv1,rm)+tv2
                    d_atoms2[atom_name]['coordinate'] = coordinate2
                    coordinates2 += [coordinate2]

            for atom_name in d_atoms1.keys():
                ## use only selected atoms
                if l_atoms != [] and atom_name not in l_atoms:
                    continue
                ## use backbone atoms only if mutated residue
                if mutation == True and atom_name not in ['N','CA','C','O']:
                    continue
                if atom_name not in d_atoms2.keys():
                    continue
                if 'REMARK' in d_atoms1[atom_name].keys():
                    stop1
                    continue
                if 'REMARK' in d_atoms2[atom_name].keys():
                    stop2
                    continue
                ## append coordinates to list of coordinates
                coordinate1 = d_atoms1[atom_name]['coordinate']
                coordinate2 = d_atoms2[atom_name]['coordinate']
                if not rmsd:
                    coordinates1 += [coordinate1]
                    coordinates2 += [coordinate2]
                if rmsd:
                    SS += [sum((coordinate2-coordinate1)**2)]

            ##
            ## append lines
            ##
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

        if rmsd:
            return coordinates1, coordinates2, rescount, lines1, lines2
        else:
            return coordinates1, coordinates2, rescount

        return coordinates1, coordinates2


    def identify_interpdb_equivalent_chains_from_structure(
        self, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
        d_pdb, d_header,
        biomolecule1, biomolecule2,
        d_biomolecules1, d_biomolecules2,
        d_ATOMseq,
        ):

        import os, biounit, copy

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
##                dontalignseqtwice
##                d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(
##                    d_pdb, d_header, pdb1, pdb2, rep_chain1, rep_chain2,
##                    )
##                if l1SEQRES != l1 or l2SEQRES != l2:
##                    print pdb1, pdb2
##                    expected_not_critical
##                if ATOMseq1 != d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['s1']:
##                    print pdb1, pdb2
##                    expected_not_critical
##                if ATOMseq2 != d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['s2']:
##                    print pdb1, pdb2
##                    expected_not_critical
##                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['d_res_nos1'] = d_res_nos1
##                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['d_res_nos2'] = d_res_nos1
##                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l1ATOM'] = l1
##                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['l2ATOM'] = l2
##                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['seqATOM1'] = l1
##                d_chains_interpdb_sequence_similar[rep_chain1][rep_chain2]['seqATOM2'] = l2
            
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
                            d_pdb, tchain = self.matrixtransformation(d_pdb,pdb,chain,matrix,matrix_no)
                            tchains += [tchain]
                            if len(tchain) > 1:
                                if int(tchain[2:]) > 1:
                                    transformations = True
                l_tchains += [tchains]
            d_biomolecules[pdb]['l_tchains'] = l_tchains


        print d_biomolecules
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
        print 'expected chain combination of REMARK350 transformations'
        (
            rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2
            ) = self.calculate_rmsd_for_multiple_chains(
                 tchains1,tchains2,d_pdb,pdb1,pdb2,d_header,
                 d_chains_interpdb_sequence_similar,
                 d_chains_intrapdb_sequence_identical,
                 d_ATOMseq,
                 )

        ## correct transformation
        if rmsd < self.maxrmsd:
            l_equivalent_chains = [tchains1,tchains2]
            return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations
        ## incorrect transformation and 1 chain v 1 chain
        elif len(tchains1) == 1 and len(tchains2) == 1:
            if rmsd > self.maxrmsd_wrong:
                self.incorrecttransformation(d_header,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,rmsd,)
            l_equivalent_chains = [tchains1,tchains2]
            return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations

        ## incorrect transformation and not 1 chain v 1 chain
        else:

            ##
            ## checks
            ##
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

            ##
            ## do PISA if 350 symmetry operations do not yield identical biounits
            ##
            if rmsd > self.maxrmsd_wrong:
                l_chains1 = d_biomolecules[pdb1]['l_chains']
                l_chains2 = d_biomolecules[pdb2]['l_chains']
                chains1 = []
                for chains in l_chains1:
                    chains1 += chains
                chains2 = []
                for chains in l_chains2:
                    chains2 += chains
                ##
                ## PISA transformations
                d_transformations_PISA1,status1 = biounit.biounit().parse_pisa_multimers(pdb1)
                d_transformations_PISA2,status2 = biounit.biounit().parse_pisa_multimers(pdb2)
                for assembly1 in d_transformations_PISA1.keys():
                    d_pdb, tchains1_PISA = self.apply_PISA_transformations(d_pdb,pdb1,assembly1,d_transformations_PISA1, chains1)
                    print '*1*',assembly1,tchains1_PISA
                    if len(tchains1_PISA) == 0:
                        continue
                    for assembly2 in d_transformations_PISA2.keys():
                        d_pdb, tchains2_PISA = self.apply_PISA_transformations(d_pdb,pdb2,assembly2,d_transformations_PISA2, chains2)
                        print '*2*',assembly2,tchains2_PISA
                        if len(tchains2_PISA) == 0:
                            continue
                        if len(tchains1_PISA) != len(tchains2_PISA):
                            continue
                        print 'PISA', pdb1, pdb2, tchains1_PISA, tchains2_PISA, chains1, chains2
                        print 'expected chain combination of PISA transformations'
                        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                            tchains1_PISA,tchains2_PISA,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                            )

                        ##
                        ## check if the expected correct combination of chains gives a low rmsd
                        ##
                        if rmsd < self.maxrmsd_wrong: ## e.g. 1hqy,1ht2
                            l_equivalent_chains = [tchains1_PISA,tchains2_PISA]
                            fd = open('pisa1.txt','a')
                            fd.write('%s %s %s %s\n' %(pdb1, pdb2, assembly1, assembly2,))
                            fd.close()
                            return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations

            ##
            ## do different combinations of chain pairing than the expected (using remark350 transformations!)
            ##
            print 'different chain combinations of remark350 transformations', pdb1,pdb2,biomolecule1,biomolecule2
            (
                l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                ) = self.shuffle_chains_and_calculate_rmsd(
                    d_pdb,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
                    pdb1,pdb2,l_tchains1,l_tchains2,
                    d_ATOMseq,
                    )
            if rmsd < self.maxrmsd_wrong: ## e.g. 2dtz,2hq5 for which PISA transformations does not exist and rmsd > 2.5
                return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations


            ##
            ## do different combinations of chain pairing than the expected (using pisa transformations!)
            ##
            for assembly1 in d_transformations_PISA1.keys():
                d_pdb, tchains1_PISA = self.apply_PISA_transformations(d_pdb,pdb1,assembly1,d_transformations_PISA1, chains1)
                if len(tchains1_PISA) == 0:
                    continue
                for assembly2 in d_transformations_PISA2.keys():
                    d_pdb, tchains2_PISA = self.apply_PISA_transformations(d_pdb,pdb2,assembly2,d_transformations_PISA2, chains2)
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

                    print 'different chain combinations of PISA transformations'
                    (
                        l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                        ) = self.shuffle_chains_and_calculate_rmsd(
                            d_pdb,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
                            pdb1,pdb2,l_tchains1_PISA,l_tchains2_PISA,
                            d_ATOMseq,
                            )
                    if rmsd < self.maxrmsd_wrong: ## e.g. 1ryz,1t0u
                        fd = open('pisa2.txt','a')
                        fd.write('%s %s %s %s\n' %(pdb1, pdb2, assembly1, assembly2,))
                        fd.close()
                        return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations

            ##
            ## do 290 symmetry operations if PISA does not yield identical biounits
            ##
            if rmsd > self.maxrmsd_wrong:

                ##
                ## REMARK 290 transformations, all combinations
                ##
                print 'apply REMARK290 transformations'
                (
                    l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                    ) = self.remark290combinations(
                        d_pdb,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,d_header,rmsd,
                        d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
                        )
                if rmsd < self.maxrmsd_wrong: ## e.g. 1b8e, 1bsy
                    return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations


        if rmsd > self.maxrmsd_wrong: ## e.g. 2hqe, 2hqx
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
            if rmsd > self.maxrmsd_wrong:
                self.incorrecttransformation(d_header,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,rmsd,)

        ##
        ## return 
        ##
        return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2, transformations


    def shuffle_chains_and_calculate_rmsd(
        self,d_pdb,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
        pdb1,pdb2,l_tchains1,l_tchains2,
        d_ATOMseq,
        ):

        chains2 = []
        chains1 = []
        ksort = []
        for i in range(len(l_tchains1)):
            tchains1 = list(l_tchains1[i])
            tchains2 = list(l_tchains2[i])
            if len(tchains1) != len(tchains2):
                print tchains1, tchains2
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
                    rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                        rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                        )
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

                        ##
                        ## virus
                        ##
                        if len(tchains1) % 60 == 0 and (j-1) % 60 == 0 and len(tchains1) >= 60:

                            rmsdchains1 = chains1+tchains1[:j-1+60]
##                            print 'rmsdchains1', rmsdchains1

                            ## e.g. chains A,B,C of 1w39.pdb, 2fz1.pdb
                            if ksort != range(1,61) and ksort != []:
                                rmsdchains2 = list(chains2)
                                for k in ksort:
                                    rmsdchains2 += [tchains2[k-1]]

##                                print 'rmsdchains2', rmsdchains2
                                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                                    rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False,
                                    )
                                print rmsd, 'previous sequence'
                                print ksort
                                if rmsd < self.maxrmsd:
                                    chains2 = rmsdchains2
                                    tchains2 = tchains2[60:]
                                    jskip = range(j,j+60)
                                    ksort = ksort
                                    continue
## maybe A-transformation == C-transformation != B-transformation... then list of ksorts instead of resetting...

                            ## e.g. chain C of 1aq4.pdb, 2bq5.pdb
                            if len(tchains1) > 60:
                                rmsdchains2 = list(chains2)+tchains2[:60]

##                                print 'rmsdchains2', rmsdchains2
                                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                                    rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                                    )
                                print rmsd, 'integer sequence'
                                if rmsd < self.maxrmsd:
                                    chains2 = rmsdchains2
                                    tchains2 = tchains2[60:]
                                    jskip = range(j,j+60)
                                    ksort = [] ## reset ksort
                                    continue
                                if (
                                    (pdb1 == '1w39' and pdb2 == '2fz1') or
                                    (pdb1 == '1auy' and pdb2 == '1w39') or
                                    (pdb1 in ['1aq3','1aq4','1zdi',] and pdb2 in ['2b2e','2b2g','2bny',])
                                    ):
                                    if pdb1 == '1w39' and pdb2 == '2fz1':
                                        seq = [1, 2, 3, 4, 5, 35, 31, 32, 33, 34, 29, 30, 26, 27, 28, 18, 19, 20, 16, 17, 10, 6, 7, 8, 9, 51, 52, 53, 54, 55, 48, 49, 50, 46, 47, 14, 15, 11, 12, 13, 23, 24, 25, 21, 22, 44, 45, 41, 42, 43, 60, 56, 57, 58, 59, 36, 37, 38, 39, 40]
                                    if pdb1 == '1auy' and pdb2 == '1w39':
                                        seq = [1, 2, 3, 4, 5, 22, 23, 24, 25, 21, 38, 39, 40, 36, 37, 19, 20, 16, 17, 18, 44, 45, 41, 42, 43, 13, 14, 15, 11, 12, 7, 8, 9, 10, 6, 56, 57, 58, 59, 60, 48, 49, 50, 46, 47, 34, 35, 31, 32, 33, 26, 27, 28, 29, 30, 52, 53, 54, 55, 51]
                                    if pdb1 in ['1aq3','1aq4','1zdi',] and pdb2 in ['2b2e','2b2g','2bny',]:
                                        seq = [2, 3, 4, 5, 1, 7, 8, 9, 10, 6, 12, 13, 14, 15, 11, 17, 18, 19, 20, 16, 22, 23, 24, 25, 21, 27, 28, 29, 30, 26, 32, 33, 34, 35, 31, 37, 38, 39, 40, 36, 42, 43, 44, 45, 41, 47, 48, 49, 50, 46, 52, 53, 54, 55, 51, 57, 58, 59, 60, 56]
                                    rmsdchains2 = list(chains2)
                                    for k in seq:
                                        rmsdchains2 += [tchains2[k-1]]
                                    rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                                        rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                                        )
                                    print rmsd, 'manual sequence'
                                    if rmsd < self.maxrmsd:
                                        chains2 = rmsdchains2
                                        tchains2 = tchains2[60:]
                                        jskip = range(j,j+60)
                                        ksort = [1, 2, 3, 4, 5, 35, 31, 32, 33, 34, 29, 30, 26, 27, 28, 18, 19, 20, 16, 17, 10, 6, 7, 8, 9, 51, 52, 53, 54, 55, 48, 49, 50, 46, 47, 14, 15, 11, 12, 13, 23, 24, 25, 21, 22, 44, 45, 41, 42, 43, 60, 56, 57, 58, 59, 36, 37, 38, 39, 40]
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
                            (
                                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2,
                                ) = self.calculate_rmsd_for_multiple_chains(
                                    rmsdchains1,rmsdchains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=False
                                    )
                            if len(tchains1) >= 7: ## e.g. 1j4z,2eu1
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
                        if k == len(tchains1)-j and len(tchains1) >= 60: ## temporary!!!
                            fd = open('virus.txt','a')
                            fd.write('%s %s %s\n' %(pdb1,pdb2,chains2))
                            fd.close()
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

                        if rmsd > self.maxrmsd_wrong:
                            return None, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2
                            
            chains1 += tchains1

        rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
            chains1,chains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True,
            )

        l_equivalent_chains = [chains1, chains2,]

        return l_equivalent_chains, rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2


    def apply_PISA_transformations(self,d_pdb,pdb,assembly,d_transformations_PISA, chains):

        t_chains = []

        for chain in chains:
            if chain not in d_transformations_PISA[assembly]['chains'].keys(): ## e.g. 1j4z
                continue
##                tchains = []
##                return d_pdb, tchains
            for molecule in d_transformations_PISA[assembly]['chains'][chain].keys():
                r = d_transformations_PISA[assembly]['chains'][chain][molecule]['r']
                t = d_transformations_PISA[assembly]['chains'][chain][molecule]['t']
                matrix = [
                    [r[0][0],r[0][1],r[0][2],t[0],],
                    [r[1][0],r[1][1],r[1][2],t[1],],
                    [r[2][0],r[2][1],r[2][2],t[2],],
                    ]
                d_pdb, tchain = self.matrixtransformation(d_pdb,pdb,chain,matrix,molecule,prefix='pisa',)
                t_chains += [tchain]

        return d_pdb, t_chains



    def remark290combinations(
        self,d_pdb,pdb1,pdb2,biomolecule1,biomolecule2,chains1,chains2,d_header,rmsd,
        d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,
        ):

        sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
        import combinatorics

##        chains1 = d_header[pdb1]['REMARK350'][biomolecule1]['chains'].keys()
##        chains2 = d_header[pdb2]['REMARK350'][biomolecule2]['chains'].keys()
        operators1 = d_header[pdb1]['REMARK290'].keys()
        operators2 = d_header[pdb2]['REMARK290'].keys()

        if 'REMARK350' not in d_header[pdb1].keys() and 'REMARK350' not in d_header[pdb2].keys():
            if len(operators1) == 1 or len(operators2) == 1:
                n_matrices1 = 1
                n_matrices2 = 1
            else:
                print operators1, operators2
                print chains1, chains2
        else:
            if 'REMARK350' not in d_header[pdb1].keys():
                if len(operators1) == 1:
                    n_matrices1 = 1
                else:
                    n_matrices1 = len(chains1)
            else:
                n_matrices1 = len(d_header[pdb1]['REMARK350'][biomolecule1]['matrices'].keys())
            if 'REMARK350' not in d_header[pdb2].keys():
                if len(operators2) == 1:
                    n_matrices2 = 1
                else:
                    n_matrices2 = len(chains2)
            else:
                n_matrices2 = len(d_header[pdb2]['REMARK350'][biomolecule2]['matrices'].keys())

        print n_matrices1, n_matrices2
        print operators1, operators2
        if (
            n_matrices1*len(chains1) != n_matrices2*len(chains2)## and
##            set([pdb1,pdb2]) != set(['1qve','1rg0']) and
##            set([pdb1,pdb2]) != set(['1m3s','1viv']) and
##            set([pdb1,pdb2]) != set(['1r4c','1tij'])
            ):
            print d_header[pdb1]
            print d_header[pdb2]
            print chains1, n_matrices1, operators1
            print chains2, n_matrices2, operators2
            print d_header[pdb1].keys()
            print d_header[pdb2].keys()
##            stop_different_sized_multimers_check_remark350records

        if len(operators1) == 1 and n_matrices1 > 1:
            operator_combinations1 = n_matrices1*[operators1]
        else:
            operator_combinations1 = combinatorics.permutation_wo_rep(operators1, n_matrices1)
        operator_combinations1.sort()

        if len(operators2) == 1 and n_matrices2 > 1:
            operator_combinations2 = n_matrices2*[operators2]
        else:
            operator_combinations2 = combinatorics.permutation_wo_rep(operators2, n_matrices2)
        operator_combinations2.sort()

        ## previous knowledge about a correct REMARK290 transformation?
        if pdb1 in self.d_symop.keys():
            operator_combinations1 = [self.d_symop[pdb1]]+operator_combinations1
        if pdb2 in self.d_symop.keys():
            operator_combinations2 = [self.d_symop[pdb2]]+operator_combinations2

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
                    d_pdb, tchain1 = self.matrixtransformation(d_pdb,pdb1,chain1,matrix1,operator1)
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
                        d_pdb, tchain2 = self.matrixtransformation(d_pdb,pdb2,chain2,matrix2,operator2)
                        tchains2 += [tchain2,]

                if len(tchains1) != len(tchains2):
##                    print operator_combinations1
##                    print operator_combinations2
                    if (
                        set([pdb1,pdb2]) != set(['1qve','1rg0']) and
                        set([pdb1,pdb2]) != set(['1m3s','1viv']) and
                        set([pdb1,pdb2]) != set(['1r4c','1tij']) and
                        set([pdb1,pdb2]) != set(['2qkt','2qku']) ## 350biounit = monomer AND dimer
                        ):
                        print chains1, tchains1
                        print chains2, tchains2
                        print pdb1, pdb2
                        print d_header[pdb1]['REMARK350'][biomolecule1]['matrices'].keys()
                        print d_header[pdb2]['REMARK350'][biomolecule2]['matrices'].keys()
                        print self.cluster, 'cluster'
                        notexpected
                    else:
                        continue
                print tchains1, tchains2
                rmsd, n_chains, n_residues, n_coordinates, tv1, rm, tv2 = self.calculate_rmsd_for_multiple_chains(
                    tchains1,tchains2,d_pdb,pdb1,pdb2,d_header,d_chains_interpdb_sequence_similar,d_chains_intrapdb_sequence_identical,d_ATOMseq,verbose=True
                    )
                if rmsd < self.maxrmsd_wrong:
                    fd = open('symop.txt','a')
                    fd.write('%s %s %s %s %s\n' %(pdb1,pdb2,operator_combination1,operator_combination2,rmsd))
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
            fd.write('%s %s %s %s %3i   %5.1f %5s %10s %10s %s\n' %(pdb1, pdb2, biomolecule1, biomolecule2, len(chains1), rmsd, d_header[pdb1]['CRYST1']==d_header[pdb2]['CRYST1'], d_header[pdb1]['CRYST1'], d_header[pdb2]['CRYST1'], len(chains1)))
            fd.close()

        return


    def matrixtransformation(self,d_pdb,s_pdb,chain,matrix,matrix_no,prefix='',):

        '''apply transformation'''

        ## matrix does not cause transformation
        if matrix == self.nontransformationmatrix:
            return d_pdb,chain

##        print 'applying REMARK350 transformation to chain %s of pdb %s' %(chain, s_pdb)

        rmatrix, tvector = self.transformationmatrix2rotationmatrix_and_translationvector(matrix)

        tchain = chain+'_'+prefix+str(matrix_no)

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
                    tcoord = numpy.dot(rmatrix, coord) + tvector
                    d_pdb[s_pdb]['chains'][tchain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['coordinate'] = tcoord


        return d_pdb, tchain


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

        chains = d_header[pdb]['SEQRES'].keys()
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

            seqi = d_header[pdb]['SEQRES'][chaini]['seq']

            for j in range(i+1,len(chains)):
                chainj = chains[j]
                if chainj in waterchains:
                    continue

                seqj = d_header[pdb]['SEQRES'][chainj]['seq']

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

            if d_header[pdb1]['SEQRES'][chain1]['type'] != 'peptide':
                continue

            seq1 = d_header[pdb1]['SEQRES'][chain1]['seq']

            if len(seq1) < self.min_len_chain:
                continue

            ## only do sequential alignment for representative chains to save time
            for chain2 in repchains2:

                if d_header[pdb2]['SEQRES'][chain2]['type'] != 'peptide':
                    continue

                seq2 = d_header[pdb2]['SEQRES'][chain2]['seq']

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
                        n_chainmutations, l_chainmutations = self.point_mutations(seq1, seq2)

                else:

                    if seq1 in seq2:

                        l1 = seq2.index(seq1)
                        l2 = 0
                        r1 = seq2.rindex(seq1)+(len(seq2)-len(seq1))
                        r2 = 0
                        if l1 != 0:
                            stop_check_if_correctl1l2r1r2
                        s1 = seq1[l1:len(seq1)-r1]
                        s2 = seq2
                        n_chainmutations = 0
                        l_chainmutations = []

                        print l1,l2,r1,r2
                        instance = sequence_alignment.NW(seq1,seq2)
                        s1,s2 = instance.Align(verbose=False)[:2]
                        l1 = len(s1)-len(s1.lstrip('-'))
                        l2 = len(s2)-len(s2.lstrip('-'))
                        r1 = len(s1)-len(s1.rstrip('-'))
                        r2 = len(s2)-len(s2.rstrip('-'))
                        print l1,l2,r1,r2
                        stop

                    elif seq2 in seq1:

                        l1 = 0
                        l2 = seq1.index(seq2)
                        r1 = 0
                        r2 = seq1.rindex(seq2)+(len(seq1)-len(seq2))
                        if l2 != 0:
                            stop_check_if_correctl1l2r1r2
                        s1 = seq1
                        s2 = seq2[l2:len(seq2)-r2]
                        n_chainmutations = 0
                        l_chainmutations = []

                    else:

                        ## 1st slow sequence comparison (SEQRESseq)

                        print pdb1, pdb2, chain1, chain2, 'begin seq aln of chains of len %s and %s' %(len(seq1),len(seq2))
                        instance = sequence_alignment.NW(seq1,seq2)
                        s1,s2 = instance.Align(verbose=False)[:2]
                        print 'end seq aln 1'

                        print l1,l2,r1,r2
                        l1 = len(s1)-len(s1.lstrip('-'))
                        l2 = len(s2)-len(s2.lstrip('-'))
                        r1 = len(s1)-len(s1.rstrip('-'))
                        r2 = len(s2)-len(s2.rstrip('-'))
                        print l1,l2,r1,r2

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
                        'r1':r1,'r2':r2,
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

        ## this function should take into account nucleotides
        ## unfortunately there is a conflict between:
        ## G: glycine and guanosine
        ## A: alanine and adenosine
        ## T: threonine and thymidine
        ## C: cysteine and cytidine
        ## N: aspargine and unknown nucleotide residue

        if res_name in self.d_res.keys():
            symbol = self.d_res[res_name]
        elif res_name in self.l_nucleotides:
            symbol = res_name[-1]
        else:
            symbol = 'X'

        return symbol
    

    def ATOM2seq(self, d_pdb, chain, d_header, res_no_max = None, stop_error = True):

        d_res_nos = {}
        seq = ''
        ATOMrespos = 0
        res_nos = d_pdb['chains'][chain]['residues'].keys()
        res_nos.sort()
        for i in range(len(res_nos)):
            res_no = res_nos[i]
            if res_no_max != None and res_no >= res_no_max:
                continue

            for i in range(len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])):
                iCode = d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][i]
                try:
                    res_name = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                except:
                    print d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
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
                            if res_name != d_header['HET'][chain][res_no][iCode]:
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


    def determine_if_modres(self, d_header, d_pdb, chain, res_no, iCode, res_name):

        modres = False
        if chain in d_header['MODRES'].keys():
            if res_no in d_header['MODRES'][chain].keys():
                if iCode in d_header['MODRES'][chain][res_no].keys():
                    if d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] == d_header['MODRES'][chain][res_no][iCode]:
                        modres = True
                    else:
                        print chain, res_no, iCode, res_name, d_header['MODRES'][chain][res_no][iCode]
                        print 
                        print self.cluster
                        notexpected

        return modres


    def append_missing_residues_to_sequence(self, ATOMseq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq):

        for SEQRESpos in range(SEQRESrange1,SEQRESrange2):
            ATOMseq += SEQRESseq[SEQRESpos]
            d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}

        return ATOMseq, d_res_nos_SEQRES

    
    def identify_missing_nonterminal_residues(self, d_pdb, chain, SEQRESseq, ATOMseqgaplen, d_res_nos_ATOM):

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
            res_name = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
            if res_name == 'HOH':
                continue
            res_symbol = self.res_name2res_symbol(res_name)

            ## append gap between ATOMseq and SEQRESseq
            if seq+res_symbol != SEQRESseq[:SEQRESpos+1]:
                ## e.g. 2i0b.pdb
##                print d_pdb[pdb]['chains'][chain]['residues'][res_no]['l_iCodes']
##                print pdb, chain, res_no, iCode, res_name
##                print seq+res_symbol
##                print SEQRESseq[:SEQRESpos+1]
                if d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] == 'HETATM':
                    continue
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}
                ATOMseqgaplen += 1 ## not a gap but a reversal in the case of 2bfk.pdb!!!
                try:
                    seq += SEQRESseq[SEQRESpos]
                except:
                    print
                    print pdb, chain, res_no, iCode, res_name, SEQRESpos, SEQRESrange1, SEQRESrange2
                    print d_res_nos_ATOM[ATOMpos-1]['res_no'], d_res_nos_ATOM[ATOMpos-1]['iCode']
                    print len(SEQRESseq), SEQRESseq
                    print len(seq), seq
                    print self.pdb1, self.pdb2
                    print SEQRESseq[SEQRESrange1:]
                    print seq[SEQRESrange1:]
                    for i in range(len(seq)):
                        if seq[:i] != SEQRESseq[:i]:
                            print seq[:i]
                            print SEQRESseq[:i]
                            print self.cluster
                            stop1
                    print self.cluster
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
##        print ATOMseq
##        print SEQRESseq

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


    def alignATOMseq(self, d_pdb, d_header, pdb1, pdb2, chain1, chain2):

        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
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
            SEQRESseq = d_header[pdb]['SEQRES'][chain]['seq']
            ATOMseq,d_res_nos = self.ATOM2seq(d_pdb[pdb], chain, d_header[pdb])
            ## find missing residues, Nterminal; align Nterminal SEQRESseq and ATOMseq
            ATOMseqindentation = self.identify_missing_terminal_residues(d_pdb, pdb, chain, ATOMseq, SEQRESseq)
            ## find missing residues, Cterminal or nonterminal; align SEQRESseq and ATOMseq
            ATOMseq,d_res_nos = self.identify_missing_nonterminal_residues(d_pdb[pdb], chain, SEQRESseq, ATOMseqindentation, d_res_nos)
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

        else: ## pre-align SEQRESsequences instead!!!

            ## 2nd slow sequence comparison (ATOMseq)

            print pdb1, pdb2, chain1, chain2, 'init seq aln, chain len', len(ATOMseq1), len(ATOMseq2)
            instance = sequence_alignment.NW(ATOMseq1,ATOMseq2)
            s1,s2 = ATOMs1,ATOMs2 = instance.Align(verbose=False)[:2]
            print 'end seq aln 2'

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


    def parse_header(self, s_pdb, stop_error = True):

        ##
        ## read lines
        ##
        fd = open('%s%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
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

            elif record == 'REMARK': ## section 2
                d_header = self.parse_recordREMARK(d_header, line, i, lines)

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
                    d_header['HET'][chain][res_no][iCode] = hetID
                elif d_header['HET'][chain][res_no][iCode] != hetID:
                    if stop_error == True:
                        print d_header['HET'][chain][res_no][iCode], hetID
                        print chain, res_no, iCode
                        print s_pdb, line
                        notexpected

            elif record == 'MODRES':
                hetID = line[12:15].strip()
                chain = line[16]
                res_no = int(line[18:22])
                iCode = line[22]
                res_name = line[24:27].strip()
                txt = line[29:80].strip()
                if hetID in set(self.d_res.keys())-set(['MSE']) and res_name in set(self.d_res.keys())-set(['MSE']):
                    continue
##                if txt not in [] and hetID in set(self.d_res.keys()+self.l_nucleotides)-set(['MSE']):
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
                    d_header['MODRES'][chain][res_no][iCode] = hetID
                elif hetID != d_header['MODRES'][chain][res_no][iCode]:
                    print line, s_pdb
                    stop

            elif record == 'TITLE': ## section 2
                if not 'TITLE' in d_header.keys():
                    d_header['TITLE'] = line[10:].strip()
                else:
                    if d_header['TITLE'][-1] == '-':
                        d_header['TITLE'] += line[10:].strip()
                    else:
                        d_header['TITLE'] += ' '+line[10:].strip()

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
                d_header['HEADER'] = line[10:50].strip()

            elif record == 'SPRSDE': ## section 2
                sIDcodes = line[31:70].lower().split()
                for sIDcode in sIDcodes:
                    if sIDcode not in d_header.keys():
                        d_header[sIDcode] = {}
                    d_header[sIDcode]['SPRSDE'] = True

            elif record == 'EXPDTA': ## section 2
                methods = line[10:].strip().split(',')[0]
                if methods[:3] == 'NMR':
                    methods = 'NMR'
                elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                    methods = 'X-RAY'
                d_header['EXPDTA'] = methods

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                if spacegroup in self.d_crystalsynonyms.keys():
                    spacegroup = self.d_crystalsynonyms[spacegroup]
                d_header['CRYST1'] = spacegroup

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
        for chain in d_header['SEQRES'].keys():
            if d_header['SEQRES'][chain]['type'] == 'peptide':
                peptidechains += chain
                if len(d_header['SEQRES'][chain]['seq']) > self.min_len_chain:
                    proteinchains += chain
            elif d_header['SEQRES'][chain]['type'] == 'nucleotide':
                nucleotidechains += chain
            elif d_header['SEQRES'][chain]['type'] == 'saccharide':
                saccharidechains += chain

        d_header['proteinchains'] = proteinchains

##        if s_pdb not in ['1ady','1bhj']:
##            chains = d_header['SEQRES'].keys()
##            if 'REMARK350' not in d_header.keys() and len(d_header['SEQRES'].keys()) > 1 and d_header['EXPDTA'] != 'NMR' and len(saccharidechains) != len(d_header['SEQRES'].keys()):
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
##                        print d_header['SEQRES'].keys()
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


    def parse_atom_no_range(self, d_conect, record, atom_no):
        
        if not record in d_conect.keys():
            d_conect[record] = [[atom_no,atom_no]]
        elif d_conect[record][-1][1] == atom_no-1:
            d_conect[record][-1][1] = atom_no
        else:
            d_conect[record] += [[atom_no,atom_no]]

        return d_conect


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):

        import numpy

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
            d_header['SEQRES'] = {}
        if chain not in d_header['SEQRES'].keys():
            d_header['SEQRES'][chain] = {}
        if not 'type' in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['type'] = 'N/A'

        l_residues = line[19:70].split()

        s_residues = ''
        for i in range(len(l_residues)):
            residue = l_residues[i]
            if residue in self.d_res.keys():
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'peptide'
                elif d_header['SEQRES'][chain]['type'] != 'peptide':
                    stop
                s_residues += self.d_res[residue]
            elif residue in self.l_nucleotides:
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'nucleotide'
                elif d_header['SEQRES'][chain]['type'] != 'nucleotide':
                    stop
                s_residues += residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'saccharide'
                elif d_header['SEQRES'][chain]['type'] != 'saccharide':
                    stop
                s_residues += residue
            else:
                if residue == 'UNK': ## e.g. 1pny.pdb
                    if d_header['SEQRES'][chain]['type'] == 'N/A':
                        d_header['SEQRES'][chain]['type'] = 'peptide'
                s_residues = 'X'

        if 'seq' not in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['seq'] = ''
        d_header['SEQRES'][chain]['seq'] += s_residues

        if 'seq3' not in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['seq3'] = []
        d_header['SEQRES'][chain]['seq3'] += l_residues

        return d_header


    def parse_coordinates(
        self, s_pdb, d_header, verbose = False, parse_molecules = True,
        ):

        print s_pdb

        ## deep copy header, because REMARK465 records are deleted during loop over biomolecules
        d_header = copy.deepcopy(d_header)

        ##
        ## read lines
        ##
        fd = open('%s%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
        lines = fd.readlines()
        fd.close()

        ##
        ## set dictionaries
        ##
        d_atomnos = {}
        d_CONECT = {}
        d_coordinates = {}
        d_ATOMseq = {}

        ##
        ## loop over lines
        ##
        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'ATOM':
                d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
                d_ATOMseq = self.build_ATOMseq(line, lines, i, d_ATOMseq, d_header,)
                d_atomnos[d_line['atom_no']] = d_line

            elif record == 'HETATM':
                res_name = line[17:20].strip()
                ## water
                if res_name in ['D2O','H2O',]:
                    print s_pdb, res_name
                    stop
                ## water must be parsed in case there is a connection to it (e.g. 2bpb)
                if res_name in ['HOH','DOD',]: ## DOD in 2d4j
                    d_coordinates, d_line, = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
                ## MSE
                elif res_name in self.d_modres.keys():
                    d_coordinates, d_line, = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
                ## (poly)saccharide or other hetero compound
                else:
                    d_coordinates, d_line, = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
                atom_no = d_line['atom_no']
                d_atomnos[atom_no] = d_line

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


        if parse_molecules == True:
            d_molecules = self.build_dictionary_of_molecules(d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=verbose,)
        else:
            d_molecules = {}

        for chain in d_header['SEQRES']:

            ## skip if not peptide chain
            type = d_header['SEQRES'][chain]['type']
            if type != 'peptide':
                continue

            ## skip if all residues are unknown
            if len(d_header['SEQRES'][chain]['seq3'])*['UNK'] == d_header['SEQRES'][chain]['seq3']:
                continue

            ## append remaining REMARK465 residues
            if 'REMARK465' in d_header.keys():
                if chain in d_header['REMARK465']['chains'].keys():
                    if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                        l_REMARK465_res_nos = d_header['REMARK465']['chains'][chain]['residues'].keys() ## e.g. 1cd0
                        l_REMARK465_res_nos.sort()
                        d_ATOMseq,d_header = self.append_ATOMseq(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,False,)

            if d_ATOMseq[chain]['seq'] != d_header['SEQRES'][chain]['seq3']:
                print 'SEQRES', d_header['SEQRES'][chain]['seq3']
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print chain
                print 'SEQRES', len(d_header['SEQRES'][chain]['seq3'])
                print 'ATOM'  , len(d_ATOMseq[chain]['seq'])
                print pdb,chain
                stop_different_length_SEQRES_ATOM

        return d_coordinates, d_molecules, d_ATOMseq


    def append_ATOMseq(
        self,record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,append_REMARK465,append_ATOM
        ):

        if append_REMARK465 == True:
            for res_no_REMARK465 in l_REMARK465_res_nos:
                if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                    l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                    for iCode_REMARK465 in l_iCodes_REMARK465:
                        res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                        d_ATOMseq[chain]['seq'] += [res_name_REMARK465]
                        d_ATOMseq[chain]['res_nos'] += [res_no_REMARK465]
                        d_ATOMseq[chain]['iCodes'] += [iCode_REMARK465]
                        d_ATOMseq[chain]['altlocs'] += [' ']
                        d_ATOMseq[chain]['records'] += ['REMARK465']
                    del d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]
    ##            else: ## not 3b95
    ##                break

        if append_ATOM == True:
            d_ATOMseq[chain]['seq'] += [res_name_ATOM]
            d_ATOMseq[chain]['res_nos'] += [res_no]
            d_ATOMseq[chain]['iCodes'] += [iCode]
            d_ATOMseq[chain]['altlocs'] += [altloc]
            d_ATOMseq[chain]['records'] += [record]

        return d_ATOMseq, d_header


    def build_ATOMseq(self,line,lines,i,d_ATOMseq,d_header):

        record = line[:6].strip()
        altloc = line[16]
        res_name_ATOM = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]

        skip = False

        ## modified residue?
        SEQRESres = self.check_if_SEQRESres(res_name_ATOM,record,d_header,chain,res_no,iCode,)
        if SEQRESres == False:
            skip = True

        ## peptide chain?
        if chain in d_header['SEQRES']:
            type = d_header['SEQRES'][chain]['type']
            if type != 'peptide':
                skip = True
        else:
            skip = True

        if skip == True:
            return d_ATOMseq

        if not chain in d_ATOMseq.keys():
            d_ATOMseq[chain] = {
                'seq':[],'iCodes':[],'res_nos':[],'altlocs':[],'records':[],
                }

        if lines[i-1][:6].strip() in ['ATOM','HETATM','ANISOU','SIGUIJ','SIGATM',]:
            res_no_prev = int(lines[i-1][22:26])
            iCode_prev = lines[i-1][26]
            altloc_prev = lines[i-1][16]
        else:
            res_no_prev = None
            iCode_prev = None
            altloc_prev = None

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
                        ## REMARK465 before ATOM (first residue > 1 and second residue > 1) ## e.g 2h27
                        if res_no_prev == None and res_no > 1 and res_no > min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                            pass_if = True
                        ## REMARK465 before ATOM
                        ## 1sgf,1nu0
                        if res_no_prev == min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                            pass_if = True

                    if pass_if == True:

                        ## REMARK465 before ATOM?
                        ## e.g. 3bef
                        if min(d_header['REMARK465']['chains'][chain]['residues'].keys()) <= res_no:
                            if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():

                                l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']
                                index1 = s_alphabet.index(min(l_iCodes_REMARK465))
                                index2 = s_alphabet.index(max(l_iCodes_REMARK465))+1
                                l_iCodes_ascending = ','.join(s_alphabet[index1:index2]).split(',')
                                l_iCodes_descending = list(l_iCodes_ascending)
                                l_iCodes_descending.reverse()
                                if l_iCodes_REMARK465 == l_iCodes_ascending:
                                    ascending = True
                                    descending = False
                                elif l_iCodes_REMARK465 == l_iCodes_descending:
                                    ascending = False
                                    descending = True

                                if len(l_iCodes_REMARK465) > 1 and iCode == ' ' and ascending == True:
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no) ## e.g. 1b8m
                                elif len(l_iCodes_REMARK465) > 1 and iCode == ' ' and descending == True:
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 2ass
                                elif iCode != ' ':
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 3bef
                                ## e.g. 2bvs
                                elif len(l_iCodes_REMARK465) == 1:
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1)
                                    l_REMARK465_res_names = []
                                    for res_no_REMARK465 in l_REMARK465_res_nos:
                                        if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                            l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                            for iCode_REMARK465 in l_iCodes_REMARK465:
                                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                                                l_REMARK465_res_names += [res_name_REMARK465]
                                    iCode_REMARK465 = l_iCodes_REMARK465[0]
                                    res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_REMARK465]['res_name']
                                    SEQRES_seq = d_header['SEQRES'][chain]['seq3'][:len(l_REMARK465_res_names)]
                                    if SEQRES_seq == l_REMARK465_res_names:
                                        pass
                                    else:
                                        l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)
                                else:
                                    stop
                            else:
                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)

                            ## e.g. 1bd7
                            l_REMARK465_res_names = []
                            for res_no_REMARK465 in l_REMARK465_res_nos:
                                if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                    l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                    for iCode_REMARK465 in l_iCodes_REMARK465:
                                        res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                                        l_REMARK465_res_names += [res_name_REMARK465]
    ##                            else: ## not 3b95
    ##                                break

                            if len(l_REMARK465_res_names) > 0:
                                l_SEQRES_res_names = d_header['SEQRES'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)]
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
                                ## 2hu9
                                elif (
                                    len(l_REMARK465_res_names) == 1 and
                                    l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                                    res_name_ATOM == l_SEQRES_res_names[0] and
                                    res_no > min(l_REMARK465_res_nos)
                                    ):
                                    REMARK465_before_ATOM = True
                                ## REMARK465 not before ATOM
                                ## e.g. 3bef,1bd7,4htc
                                else:
                                    REMARK465_before_ATOM = False
                            ## REMARK465 not before ATOM
                            ## e.g. 2a0q
                            else:
                                REMARK465_before_ATOM = False

                        if not REMARK465_before_ATOM == True: ## e.g. 103l
                            if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                res_no_REMARK465 = res_no
                                l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                iCode_REMARK465 = l_iCodes_REMARK465[0]
                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']

            try:
                res_name_SEQRES = d_header['SEQRES'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
            except:
                res_name_SEQRES = 'N/A'

            if REMARK465_before_ATOM == False and res_name_ATOM != res_name_SEQRES and res_name_ATOM in ['FOR','HOA','ACE',] and (res_no == min(d_coordinates['chains'][chain]['residues'].keys()) or res_no == max(d_coordinates['chains'][chain]['residues'].keys())):
                fd = open('formyl.txt','a')
                fd.write('%s %s %s\n' %(pdb, chain, res_name_ATOM))
                fd.close()
                return d_coordinates, {
                    'chain':chain,
                    'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
                    'atom_no':atom_no,'atom_name':atom_name,'element':element,
                    }, d_ATOMseq

            ## REMARK465 after ATOM
            if REMARK465_before_ATOM == False and res_name_ATOM == res_name_SEQRES:# and res_name_REMARK465 != res_name_SEQRES:
                d_ATOMseq,d_header = self.append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,None,False,True,)
            ## REMARK465 before ATOM (certain)
            elif REMARK465_before_ATOM == True:# and res_name_ATOM != res_name_SEQRES:# and res_name_REMARK465 == res_name_SEQRES:
                d_ATOMseq,d_header = self.append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,True,)
            else:
                print '---'
                SEQRES_res_name = d_header['SEQRES'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
                SEQRES_seq = d_header['SEQRES'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print chain, res_no, iCode
                print line
                print d_ATOMseq[chain]['seq']
                print d_header['SEQRES'][chain]['seq3']
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
                elif res_name_REMARK465 == None and min(d_header['REMARK465']['chains'][chain]['residues'].keys()) < res_no: ## e.g. 3c5w
                    print chain,res_no
                    print d_header['REMARK465']['chains'][chain]['residues'].keys()
                    stop_temp_broken
                    pass
                ## REMARK465 after ATOM
                elif (
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
                    ## 1jly
                    if res_name_ATOM in ['FOR','HOA','ACE',] and (res_no == min(d_coordinates['chains'][chain]['residues'].keys()) or res_no == max(d_coordinates['chains'][chain]['residues'].keys())):
                        fd = open('formyl.txt','a')
                        fd.write('%s %s %s\n' %(pdb, chain, res_name_ATOM))
                        fd.close()
                        return d_coordinates, {
                            'chain':chain,
                            'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
                            'atom_no':atom_no,'atom_name':atom_name,'element':element,
                            }, d_ATOMseq
                    print '*******'
                    print 'ATOM  ', d_ATOMseq[chain]['seq']
                    print 'SEQRES', SEQRES_seq
                    print line
                    print chain,res_no
                    print 'SEQRES', SEQRES_res_name
                    print 'ATOM  ', res_name_ATOM, '***iCode***', iCode
                    print 'REMARK', res_name_REMARK465, iCode_REMARK465
                    stop_N_terminal

            if d_ATOMseq[chain]['seq'] != d_header['SEQRES'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]:
                print '*******'
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print 'SEQRES', d_header['SEQRES'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print res_name_ATOM
                print line
                print res_name_REMARK465
                print pdb, chain, res_no, iCode
                stop_sequence_difference

        return d_ATOMseq        


    def check_if_SEQRESres(self,res_name,record,d_header,chain,res_no,iCode):

        if res_name not in self.d_res.keys()+self.l_nucleotides and record == 'ATOM':
            print res_name,record
            stop

        if res_name in self.d_res.keys()+self.l_nucleotides and record in ['ATOM','REMARK465',]:
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
                                            print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                            print 'res_name', d_atomnos[tmpatom_no1]['res_name'], d_atomnos[tmpatom_no2]['res_name']
                                            print 'chain', d_atomnos[tmpatom_no1]['chain'], d_atomnos[tmpatom_no2]['chain']
                                            print 'res_no', d_atomnos[tmpatom_no1]['res_no'], d_atomnos[tmpatom_no2]['res_no']
                                            print 'iCode', d_atomnos[tmpatom_no1]['iCode'], d_atomnos[tmpatom_no2]['iCode']
                                            print 'atom_name', d_atomnos[tmpatom_no1]['atom_name'], d_atomnos[tmpatom_no2]['atom_name']
                                            print 'altloc', d_atomnos[tmpatom_no1]['altloc'], d_atomnos[tmpatom_no2]['altloc']
                                            print '*****'
                                            print (d_atomnos[tmpatom_no1]['altloc'] == 'A' and d_atomnos[tmpatom_no2]['altloc'] == 'B')
                                            print hetID1, chain1, res_no1, iCode1, atom_no1
                                            print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']
                                            print s_pdb
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                print atom_no, d_atomnos[atom_no]
                                            notexpected_different_atllocs_connected
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
                                                print self.cluster
                                                notexpected

                            l_hetIDs_long_atom_names = [
                                ## pyranoses
                                'XYP','G6D','SIA',
                                'FCT',
                                ## furanoses
                                'AHR','HPD',
                                ## dissacharides
                                'DAF',
                                ## benzoxazinoids (hydroxamic acid)
                                'HBO', ## DIMBOA from Maize
                                ## (tetra)pyrrole
                                'DBV','PEB','OPP','ZNH',
                                ## other
                                'DPM',
                                ## p-Coumaric acid (phenyl propanoid)
                                'HC4',
                                ## nucleobases/nucleosides/nucleotides
                                'DA','UMP','PGD','DT','BZG','8OG','ADP','A','TSP','CMP',
                                ## guanine
                                'DG','GDP','OMG',
                                ## cytidine
                                'C','DC','OMC',
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
                                    stop
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
                                        stop
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
                                    print hetID1, chain1, res_no1, atom_name1
                                    print hetID2, chain2, res_no2, atom_name2
                                    print atom_no1, atom_no2
                                    print s_pdb
                                    print self.cluster
                                    notexpected

                                if verbose == True and hetID1 != 'CYS' and hetID2 != 'CYS' and atom_name1 != 'SG' and atom_name2 != 'SG':
                                    print hetID1, hetID2, atom_name1, atom_name2, atom_name1_no, atom_name2_no
                                #####################################
                                ## posttranslational modifications ##
                                #####################################
                                ##
                                ## peptide bond
                                ##
                                if (
                                    hetID1 in self.d_res.keys() and hetID2 in self.d_res.keys() and
                                    atom_name1 == 'C' and atom_name2 == 'N'
                                    ):
                                    bond = 'C,N'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    hetID1 in self.d_res.keys() and hetID2 in self.d_res.keys() and
                                    atom_name2 == 'C' and atom_name1 == 'N'
                                    ):
                                    bond = 'C,N'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ##
                                ## cysteine disulphide bond
                                ##
                                if hetID1 == 'CYS' and atom_name1 == 'SG' and atom_name2[0] == 'S':
                                    bond = 'S,S'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'CYS' and atom_name2 == 'SG' and atom_name1[0] == 'S':
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
                                ## tryptophan peroxidase pi electron porphyrin interaction (make more general solution for small molecules)
                                ##
                                elif hetID1 == 'TRP' and atom_name1 == 'NE1' and hetID2 == 'PEO' and atom_name2[0] == 'O':
                                    bond = 'N'+','+atom_name2[1:]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'TRP' and atom_name2 == 'NE1' and hetID1 == 'PEO' and atom_name1[0] == 'O':
                                    bond = 'N'+','+atom_name1[1:]
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                ##
                                ## error
                                ##
                                elif atom_name1[0] == 'C' and atom_name2[0] == 'C':
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected1
                                elif atom_name1[0] == 'O' and atom_name2[0] == 'O':
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected2
                                elif atom_name1[0] == atom_name2[0]:
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected3
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
                                    chain1 == chain2 and
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
                                    chain1 == chain2 and
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
                                    hetID1 in self.d_res.keys() and
                                    hetID2 in self.d_res.keys() and
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
                                    hetID1 in self.d_res.keys() and
                                    hetID2 in self.d_res.keys() and
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
                                    hetID1 not in ['SIA','SLB','XYP',] and hetID2 not in ['SIA','SLB','XYP',]
                                    ): ##  and atom_name1[0] == atom_name2[0]
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, iCode1, iCode2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print min(d_coordinates['chains'][chain1]['residues'].keys())
                                    print min(d_coordinates['chains'][chain2]['residues'].keys())
                                    print max(d_coordinates['chains'][chain1]['residues'].keys())
                                    print max(d_coordinates['chains'][chain2]['residues'].keys())
                                    print self.cluster
                                    notexpected_bond_not11
                                ##
                                ## glycosyl and peptide bonds
                                ##
                                ## 1,1-glycoside bond (e.g. trehalose in 1do1)
                                elif int(atom_name1_no) == 1 and int(atom_name2_no) == 1:
                                    if atom_name1[0] == 'O' and atom_name2[0] == 'C' and hetID1 == 'GLC' and hetID2 == 'GLC':
                                        monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                        monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    elif atom_name1[0] == 'C' and atom_name2[0] == 'O' and hetID1 == 'GLC' and hetID2 == 'GLC':
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
                                            print 'b', bond
                                        else:
                                            print hetID1, hetID2, atom_name1, atom_name2
                                            notexpected
                                    elif hetID1 == 'GLC' and hetID2 == 'GLC':
                                        bond = '1,1'
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
                                        (atom_name1[0] in ['O','N'] and atom_name2[:2] == 'C1') or
                                        (atom_name1[0] == 'C' and atom_name2[:2] in ['O1','N1']) or
                                        ## sialic acid
                                        (atom_name1[0] == 'O' and atom_name2 == 'C2' and hetID2 in ['SIA','SLB']) or
                                        (atom_name1[0] == 'C' and atom_name2 == 'O2' and hetID2 in ['SIA','SLB']) or
                                        (atom_name1 == 'C4B' and atom_name2 == 'O4A' and hetID1 in ['XYP',] and hetID2 in ['XYP',])
                                        )
                                    ):
                                    bond = atom_name2_no+','+atom_name1_no
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    print 'd', bond
                                ## glycosyl 2
                                elif (
                                    int(atom_name2_no) != 1 and
                                    (
                                        (atom_name2[0] in ['O','N'] and atom_name1[:2] == 'C1') or
                                        (atom_name2[0] == 'C' and atom_name1[:2] in ['O1','N1']) or
                                        ## sialic acid
                                        (atom_name2[0] == 'O' and atom_name1 == 'C2' and hetID1 in ['SIA','SLB']) or
                                        (atom_name2[0] == 'C' and atom_name1 == 'O2' and hetID1 in ['SIA','SLB']) or
                                        (atom_name2 == 'C4B' and atom_name1 == 'O4A' and hetID2 in ['XYP',] and hetID1 in ['XYP',])
                                        )
                                    ):
                                    bond = atom_name1_no+','+atom_name2_no
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    print 'e', bond
                                ## error
                                else:
                                    print s_pdb, hetID1, hetID2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
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

##        print d_adjacency_forward
##        print d_adjacency_backward
                                    
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
        res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
        if res_name in self.d_saccharides.keys():
            res_name = self.d_saccharides[res_name]['stereo']

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

        ## skip if hydrogen atom
        if atom_name1[0] == 'H':
            return d_CONECT
        ## skip if a cofactor to which molecules are not *covalently* bound
        if res_name1 in self.l_cofactors+['HOH']+self.l_solutes:
            return d_CONECT

        d_atom_nos = {'intra':[],'inter':[]}

        for atom_no2 in atom_nos[1:]:
            chain2 = d_atomnos[atom_no2]['chain']
            res_no2 = d_atomnos[atom_no2]['res_no']
            iCode2 = d_atomnos[atom_no2]['iCode']
            res_name2 = d_atomnos[atom_no2]['res_name']
            atom_name2 = d_atomnos[atom_no2]['atom_name']
            ## skip if hydrogen atom
            if atom_name2[0] == 'H':
                continue
            ## skip if a cofactor to which molecules are not *covalently* bound
            if res_name2 in self.l_cofactors+['HOH']+self.l_solutes:
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

            if line[12:23].strip().upper() in ['TEMPERATURE','PH']:
                experimentaldetail_key = line[12:23].strip()
                experimentaldetail_value = line[44:].strip()
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

        elif remark == 465: ## missing residues

            d_header = self.parse_recordREMARK465(line, d_header, lines, i)

        elif remark == 470: ## missing atoms

            d_header = self.parse_recordREMARK470(line, d_header, lines, i)

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
            try:
                resolution = float(line[22:27])
            except:
                resolution = 'N/A'
            d_header['REMARK2'] = resolution

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
            biomolecules = line[23:80].replace(' ','').split(',')
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


    def parse_recordREMARK465(self, line, d_pdb, lines, i):

        ## missing residues

        if line[10:].strip() in ['M RES C SSSEQI','M RES C  SSEQI']:

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 465':
                    break

                try:
                    model = int(lines[j][12:14])
                except:
                    model = 'N/A'
                res_name = lines[j][15:18].strip()
                chain = lines[j][19]
                res_no = int(lines[j][22:26])
                iCode = lines[j][26]


                if not 'REMARK465' in d_pdb.keys():
                    d_pdb['REMARK465'] = {}
                if not 'chains' in d_pdb['REMARK465'].keys():
                    d_pdb['REMARK465']['chains'] = {}
                if not chain in d_pdb['REMARK465']['chains'].keys():
                    d_pdb['REMARK465']['chains'][chain] = {}
                if not 'residues' in d_pdb['REMARK465']['chains'][chain].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['REMARK465']['chains'][chain]['residues'].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no] = {}

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['REMARK465']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}


                ## iCode > res_name
                d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name

                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['REMARK465']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if not iCode in d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']:
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## iCode > REMARK
                if not 'REMARK' in d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] = True

                ## iCode > record
                d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = 'REMARK465'

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
                res_name = lines[j][15:18].strip()

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
                if not 'REMARK470' in d_pdb.keys():
                    d_pdb['REMARK470'] = {}
                if not 'chains' in d_pdb.keys():
                    d_pdb['REMARK470']['chains'] = {}
                if not chain in d_pdb['REMARK470']['chains'].keys():
                    d_pdb['REMARK470']['chains'][chain] = {}
                if not 'residues' in d_pdb['REMARK470']['chains'][chain].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['REMARK470']['chains'][chain]['residues'].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no] = {}

                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if len(d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes']) == 0:
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
                elif iCode not in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes']:
                    ## e.g. 2fs4 (chain A, res_no 162, iCode " ")
                    if iCode == ' ':
                        d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] = [iCode]+d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']
                    else:
                        d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

                ## iCode > atoms
                if not 'atoms' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                ## atoms > atom_name > coordinate
                for atom_name in atoms:
                    if not atom_name in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                    ## iCode > REMARK
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['REMARK'] = True

                ## iCode > res_name
                if not 'res_name' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name


        elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
            print self.pdb1
            print self.pdb2
            notexpected

        return d_pdb


    def parse_recordATOM(self, line, d_pdb, lines, i, d_header, record,):

        import numpy

        atom_no = int(line[6:11])
        atom_name = line[12:16].strip()
        altloc = line[16]
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        element = line[76:78].strip()
        coordinate = numpy.array([x, y, z])

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


        res_name_conflict = False
        ## iCode > res_name
        if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name
##        ## temp!!! temporary until pdb errors are fixed (otherwise split ATOM and HETATMs)
##        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name: ## temp!!!
##            print res_name, d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] ## temp!!!
##            print chain, res_no, iCode, altloc ## temp!!!
##            print line
##            stop_duplicate_remark465_and_coord_section
        ## check that res_name is correct (e.g. 2fes:L:1)
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name and altloc == ' ': ## 1fh2:A:30 altloc
            if record == 'HETATM':
                res_name_conflict = True
            else:
                ## change the iCode
                iCode_max = max(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
                iCode_max = self.s_alphabet[self.s_alphabet.index(iCode_max)+1]
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_max] = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'].index(iCode)] = iCode_max

                ## d_iCodes > iCode
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

        if res_name_conflict == False:
            ## res_no > l_iCodes
            if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
            ## l_iCodes > iCode
            if len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']) == 0:
                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
            elif iCode not in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
                d_pdb = self.identify_iCode_sequence(d_pdb, chain, res_no, iCode, res_name, d_header)

            ## iCode > atoms
            if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
            ## atoms > atom_name > coordinate
            if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}

            ## iCode > record
            if not 'record' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = record

        
        return d_pdb, {
            'chain':chain,
            'res_name':res_name,'res_no':res_no,'iCode':iCode,'altloc':altloc,
            'atom_no':atom_no,'atom_name':atom_name,'element':element,
            }


    def identify_iCode_sequence(self, d_pdb, chain, res_no, iCode, res_name, d_header):
        
        l_iCodes = list(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
        d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
        for iCode_prev in l_iCodes:
            index_alphabet = self.s_alphabet.index(iCode)
            if (
                ## REMARK465
                'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev].keys() or
                ## REMARK470
                {'REMARK':True} in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev]['atoms'].values()
                ):
                ATOMseq,d_res_nos = self.ATOM2seq(d_pdb, chain, d_header, res_no_max = res_no)
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
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)]
                    ):
                    l_iCodes = [iCode]+d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][:-1]
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 3) REMARK470 residues after ATOM residues
                ## e.g. 2lve.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' in l_iCodes
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 4) REMARK470 residues after ATOM residues
                ## e.g. 2j5q (chain B, res_no 54, iCode C)
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet-1] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' not in l_iCodes 
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 5) REMARK465 residues before ATOM residues
                ## e.g. 1uij.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+len(l_iCodes)] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)]
                    ):
                    None
                else: ## e.g. 1fne
                    print '*****'
                    print iCode, iCode_prev, res_symbol, index_alphabet, l_iCodes
                    print d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet]
                    print
                    print len(ATOMseq) > 0
                    print res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet-1]
                    print ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)]
                    print self.s_alphabet[index_alphabet-1] in l_iCodes
                    print 
                    print 1, ATOMseq, self.res_name2res_symbol(res_name)
                    print 2, d_header['SEQRES'][chain]['seq']
                    print chain, res_no, iCode, iCode_prev, res_name
                    expected
                break

        return d_pdb


    def spherermsd(
        self,
        pdb1, pdb2, ## pdbs
        d_header, ## sequences
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

            ATOMseq1,d_res_nos1 = self.ATOM2seq(d_pdb[pdb1], rep_chain1, d_header[pdb1])
            ATOMseq2,d_res_nos2 = self.ATOM2seq(d_pdb[pdb2], rep_chain2, d_header[pdb2])
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
                    l1, l2, len(ATOMseq1), d_ATOMseq, tv1=tv1, rm=rm, tv2=tv2
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


    def __init__(self):

        import os

        self.d_res = {
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

        ## chiral space groups

        d_spacegroups = {
            ## alpha,beta,gamma != 90
            'TRICLINIC':[
                'P 1', ## 1
##                'P 1-','P1', ## neither HM symbols nor abbrevations
##                'A 1', ## 1lks.pdb
                ],
            ## alpha != 90, beta,gamma==90
            'MONOCLINIC':[
                'P 1 2 1', ## 3
                'P 1 21 1', ## 4
                'C 1 2 1', ## 5
##                'C 2', 'C 21', 'C 1 21 1', ## C 1 2 1 abbreviations
##                'P 2', ## P 1 2 1 abbreviations
##                'P 21', ## P 1 21 1 abbreviations
##                'B 2', 'I 1 2 1', 'P 1 1 21', 'I 21', 'I 1 21 1', ## neither HM symbols nor abbrevations
                ],
            ## a != b != c (alpha,beta,gamma==90)
            'ORTHORHOMBIC':[
                'P 2 2 2', ## 16
                'P 2 2 21',
                'P 21 21 2',
                'P 21 21 21',
                'C 2 2 21',
                'C 2 2 2',
                'F 2 2 2',
                'I 2 2 2',
                'I 21 21 21', ## 24
##                'P 2 21 21', ## neither HM symbols nor abbrevations
##                'P 21 2 21', ## not a chiral space group?!
##                'P 21 21 2 A', ## P 21 21 2 error in 1b86.pdb
##                'B 2 21 2', ## 1zna.pdb
##                'B 1 1 2', ## 1qr6.pdb
                ],
            ## a != c (a == b, alpha,beta,gamma==90)
            'TETRAGONAL':[
                'P 4', ## 75
                'P 41', ## 76
                'P 42', ## 77
                'P 43', ## 78
                'I 4', ## 79
                'I 41', ## 80
                'P 4 2 2', ## 89
                'P 4 21 2', ## 90
                'P 41 2 2', ## 91
                'P 41 21 2',
                'P 42 2 2',
                'P 42 21 2',
                'P 43 2 2',
                'P 43 21 2',
                'I 4 2 2',
                'I 41 2 2',
                ],
            ## RHOMBOHEDRAL (a=b=c, alpha,beta,gamma!=90)
            ## alpha,beta,gamma != 90
            'TRIGONAL':[
                'P 3',
                'P 31',
                'P 32',
                'R 3',
                'P 3 1 2', ## 149
                'P 3 2 1',
                'P 31 1 2', ## 151
                'P 31 2 1', ## 152
                'P 32 1 2',
                'P 32 2 1',
                'R 3 2',
##                'H 3', ## R 3 equivalent
##                'H 3 2', ## P 3 2 1 equivalent
                ],
            'HEXAGONAL':[
                'P 61 2 2','P 65','P 63','P 65 2 2','P 61','P 62 2 2','P 62','P 64 2 2','P 63 2 2','P 6 2 2','P 6','P 64',
                ],
            'CUBIC':[
                'F 41 3 2','P 21 3','I 4 3 2','I 2 3','P 2 3','P 41 3 2','P 4 3 2','F 4 3 2','P 43 3 2','I 21 3','F 2 3','P 42 3 2','I 41 3 2',
                ],
##            'UNKNOWN':[
##                'A 2',
##                ]
            }

        self.d_crystalsystems = {}
        for crystalsystem in d_spacegroups.keys():
            for spacegroup in d_spacegroups[crystalsystem]:
                self.d_crystalsystems[spacegroup] = crystalsystem

        self.d_crystalsynonyms = {
            'P 2 21 21':'P 21 21 2',

            'P 21':'P 1 21 1',

            'C 2':'C 1 2 1',
            'C 21':'C 1 2 1',
            'C 1 21 1':'C 1 2 1',

            'P 2':'P 1 2 1',
            ## errors, temporary until fixed/remediated
            'P 21 21 2 A':'P 21 21 2',
            'P 2 21 21':'P 2 21 21', ## 2pnj
            }

        self.nontransformationmatrix = [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]

        self.s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        ## info from the PDB Ligand Depot (searched for cluster in chemical name)
        self.l_clusters = [
            ## iron clusters
            'CFM','CFN','CLF','CLP','CN1','CNB','CNF','F3S','FES','FS1','FS2','FS4','FSF','FSO','FSX','HC0','HC1','HF3','HF5','NFS','SF3','SF4','WCC','XCC', ## 'FS3' deprecated ('F3S' maintained)
            ## copper clusters
            'CUB','CUM','CUN','CUO',
            ## molybdenum "clusters"
            'OMO',
            ## hafnium clusters
            'PHF',
            ## zirconium clusters
            'ZRC',
            ]

        self.l_prosthetic_groups = [
            ## porphyrins (cyclic tetrapyrroles)
            ## Ferrochelatase catalyzes protophorphyrin+Fe(2+) --> protoheme + 2H(+)
            ## iron
            'HEM', ## protoporphyrin IX + Fe(II) (C3 vinyl,C8 vinyl,C18 methyl; tetramethyl,divinyl,dipropionate; *charge 2*)
            'HEC', ## Heme C (protoporphyrin IX; *charge 0*)
            'HEA', ## Heme A (C3 hydroxyfarnesyl, C8 vinyl, C18 formyl)
            'HEO', ## Heme O (C3 hydroxyfarnesyl, C8 vinyl, C18 methyl)
            '2FH', ## 2-phenylheme
            '1FH', ## 12-phenylheme
            'DDH', ## dedivinyl,diacetyl heme (C3,C8)
            'HEV', ## dedimethyl,divinyl heme
            'HDM', ## tetramethyl,divinyl,dipropionate *ester* heme
            'HAS', ## Heme-As (C3 geranylgeranyl, C8 vinyl, C18 formyl)
            'VER', ## octaethylated porphyrin
            'HEB', ## Heme B/C hybrid (tetramethyl,divinyl,dipropionate)
            'DHE', ## Heme D (heme B derivative)
            'HDD', ## cis-heme D hydroxychlorin gamma-spirolactone
            'SRM', ## siroheme (partially reduced iron-porphyrin in e.g. nitrate reductase)
            ## other metal
            'HNI', ## protoporphyrin IX + Ni(II)
            'HES', ## Zn substituted Heme C
            ]
        self.l_coenzymes = [
            'RET', ## vitamin A
##            'TPP', ## vitamin B1            
            'FMN','FAD', ## vitamin B2
            'NAD','NAP', ## vitamin B3
            'COA', ## vitamin B5
            'PLP', ## vitamin B6
            'C2F', ## vitamin B9 (5-methyl THF)
            ]

        ## list of metals also contains metalloids (e.g. Arsen)
        self.l_atoms_metal = [
            'LI','NA','K','CS', ##1a
            'BE','MG','CA', ##2a
            'AL','GA','TL', ##3a
            'PB', ##4a
            'AS', ##5a
            'V', ##3b
            'CR','MO', ##4b
            'MN', ##5b
            'FE', ##6b
            'CO', ##7b
            'NI', ##8b
            'CU', ##9b
            'ZN','CD','HG', ##10b
            ]

        ## info from the PDB Ligand Depot
        ## keys are hetIDs, values are charges (not oxidation states)
        ## hetID:[chemical formula,charge]
        self.d_ions = {

            ## unknown
            'UNX':['',''], ## e.g. 1aqn

            ## nitrate, ammonium
            'NO3':['N1 O3',-1],'NH4':['H4 N1',+1],
            ## hydroxide
            'OH' :['H1 O1',-1],
            ## phosphate
            '2HP':['O4 P1',-1],'PI' :['H1 O4 P1',-2],'PO4':['O4 P1',-3], ## different oxidation states; 'IPS' deprecated
            ## sulfate
            'SO4':['O4 S1',-2],'SOH':[3,-1],'SUL':[3,-2], ## different oxidation states
            ## carbonate
            'CO3':['C O3', -1],
            ## group1a
            'LI' :['LI1',+1],
            'NA' :['NA1',+1], ## 'NAO','NA2','NA6','NA5','NAW' deprecated
            'K'  :['K1' ,+1], ## 'KO4' deprecated
            'CS' :['CS' ,+1],
            ## group2a
            'BEF':['BE F3',-1],
            'MG' :['MG1',+2], ## 'MO3','MO1','MO2','MO4','MO5','MO6' deprecated
            'CA' :['CA1',+2],'OC1':['CA1',+2], ## 'OC3','OC5','OC6','OC7','OC8','OC2','OC4' deprecated
            'SR' :['SR' ,+2],
            'BA' :['BA' ,+2],
            ## group3a
            'AL' :['AL1',+3],'ALF' :['AL F4',-1],
            'GA' :['GA1',+3],'GA1' :['GA1',+2],
            'TL' :['TL1',+1],
            ## group4a
            'ARS':['AS1', 0],'ART':['O4 AS1',-3],'AST':-3,'TAS':['H3 O3 AS1', 0],'CAC':['C2 H6 AS O2', -1], ## different compounds
            'PB' :['PB' ,+2],
            ## group6a
            'SE' :['SE1', 0],'SE4':['O4 SE1',-2], ## different compounds
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
            'MO' :['MO1', 0],'4MO':['MO1', 0],'2MO':['MO O2',-2], ## different compounds and different oxidation states
            ## group5b
            'MN' :['MN1',+2],'MN3':['MN1',+3], ## different oxidation states; 'MW1','MW2','MW3','O4M','MN5','MN6' deprecated
            ## group6b
            'FE2':['FE1',+2],'FE' :['FE1',+3], ## different oxidation states; 'OF1','OF3','2OF' deprecated
            ## group7b
            'CO' :['CO1',+2],'3CO':['CO1',+3], ## different oxidation states; 'CO5','OCL','OCO','OCN','OCM' deprecated
            ## group8b
            'NI' :['NI1',+2],'3NI':['NI1',+3], ## different oxidation states; 'NI1','NI2','NI3','NIK' deprecated
            ## group9b
            'CU1':['CU1',+1],'CU' :['CU1',+2], ## different oxidation states; '1CU' deprecated
            ## group10b
            'ZN' :['ZN1',+2], ## 'ZN2','ZO3','ZN3','ZNO' deprecated
            'CD' :['CD1',+2],
            'HG' :['HG1',+2],
            ## Lanthanides
            'TB' :['TB1',+3],
            'YB' :['YB1',+3],
            }

        self.l_cofactors = self.l_clusters+self.l_prosthetic_groups+self.l_coenzymes+self.d_ions.keys()

        ## wikipedia buffer solution, Good's buffers
        self.l_solutes = [

            ## unknown atom or ion
            'UNX',

            ## IPA,FMT,GOL,EEE,EDO
            ##
            ## ethylene glycol, protein precipitation
            'EDO',
            ## acetic acid
            'ACY',
            ## water
            'HOH','DOD',
            ## methanol
            'MOH',
            ## di-thio-threitol (reducing agent)
            'DTT', 
            ## beta-mercapto-ethanol (reducing agent)
            'BME',
            ## bis-tris methane (buffering agent)
            'BTB',

            ## TAPS (buffering agent)
            'T3A',
            ## Bicine (buffering agent)
            'BCN',
            ## Tris (buffering agent)
            'TRS',
            ## HEPES (buffering agent)
            'EPE',
            ## TES (buffering agent)
            'NES',
            ## MOPS
            'MPO',
            ## PIPES
            'PIN',
            ## Cacodylate
            'CAC',
            ## MES
            'MES',
            ## Acetate
            'ACT',

            ## Glycerol
            'GOL',
            ]

        ## saccharides returned from a search of the ligand depot for disaccharides and monosaccharides
        self.d_saccharides = {
            ##
            ## monosaccharides, aldehydes
            ##
            ## pyranoses, hexoses
            'GLC':{'stereo':'GLC','derivate':['GLC']}, ## (alpha)-D-Glucose
##            'AGC':{'stereo':'GLC','derivate':['GLC']}, ## alpha-D-Glc (deprecated)
            'BGC':{'stereo':'GLC','derivate':['GLC']}, ## beta-D-Glc
            'GAL':{'stereo':'GAL','derivate':['GAL']}, ## (beta)-D-Galactose
            'GLA':{'stereo':'GAL','derivate':['GAL']}, ## alpha-D-Gal
##            'GLB':{'stereo':'GAL','derivate':['GAL']}, ## beta-D-Gal
            'FUC':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, alpha-L-Fucose
            'FUL':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, beta-L-Fucose
            'MAN':{'stereo':'MAN','derivate':['MAN']}, ## alpha-D-Mannose
            'BMA':{'stereo':'MAN','derivate':['MAN']}, ## beta-D-Mannose
            'ARA':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabinose
            'ARB':{'stereo':'ARA','derivate':['ARA']}, ## beta-L-Arabinose
            ## pyranoses, pentoses
            'XYS':{'stereo':'XYS','derivate':['XYS']}, ## (alpha)-D-Xylose
            'XYP':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-Xylose
            'LXC':{'stereo':'XYS','derivate':['XYS']}, ## beta-L-Xylose
            ## furanoses, hexoses
            'AHR':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabino*furano*se
            ## furanoses, pentoses
            'XYZ':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-xylo*furano*se
            ## phosphorylated aldohexopyranoses
            'G1P':{'stereo':'G1P','derivate':['GLC']}, ## alpha-D-Glc-1P
            'G6P':{'stereo':'G6P','derivate':['GLC']}, ## alpha-D-Glc-6P
            'BG6':{'stereo':'G6P','derivate':['GLC']}, ## beta-D-Glc-6P
            'G6Q':{'stereo':'G6P','derivate':['GLC']}, ## Glc-6P
            'BGP':{'stereo':'BGP','derivate':['GAL']}, ## beta-Gal-6P
            'M1P':{'stereo':'M1P','derivate':['MAN']}, ## alpha-D-Man-1P
            'M6P':{'stereo':'M6P','derivate':['MAN']}, ## alpha-D-Man-6P
            ## phosphorylated aldopentofuranoses
            'ABF':{'stereo':'ABF','derivate':['ARA']}, ## beta-D-Arabino*furano*se-5-phosphate
            ## methylated aldohexopyranoses
            'MMA':{'stereo':'MMA','derivate':['MAN']}, ## O1-methyl-mannose
            ## aminated aldohexopyranoses amine
            'AGL':{'stereo':'AGL','derivate':['GLC']}, ## 4,6-dideoxy-4-amino-alpha-D-Glucose
            ## deoxygenated aldohexopyranoses
            'G6D':{'stereo':'G6D','derivate':['GLC']}, ## 6-deoxy-alpha-D-Glucose
            ## oxygenated aldohexopyranoses
            'KBG':{'stereo':'G6D','derivate':['GLC']}, ## 2-keto-beta-D-Glucose
            ## acetylated aldohexopyranose amines
            'NAG':{'stereo':'NAG','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine
            'NBG':{'stereo':'NAG','derivate':['GLC']}, ## 1-N-Acetyl-beta-D-Glucosamine
            'NDG':{'stereo':'NAG','derivate':['GLC']}, ## 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
            '16G':{'stereo':'16G','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine-6P
            'NGA':{'stereo':'NGA','derivate':['GLC']}, ## (2)-N-Acetyl-D-Galactosamine
            ##
            ## monosaccharides, ketones
            ##
            ## furanoses, hexoses
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
            ##
            ## polysaccharides
            ##
            'MTT':{'stereo':'MTT','derivate':['MAL','MAL']}, ## MAL-b14-MAL, Maltotetraose
            ##
            ## linear saccharides (neither furanoses nor pyranoses...) and their derivatives...
            ##
            'SOR':{'stereo':'SOR','derivate':['GLC']}, ## Sorbitol/Glucitol (reduced glucose)
            'GLO':{'stereo':'GLO','derivate':['GLC']}, ## linear glucose
            'XLS':{'stereo':'XLS','derivate':['XYL']}, ## linear xylose
            'A5P':{'stereo':'ABF','derivate':['ARA']}, ## Arabinose-5-phosphate
            ##
            ## conduritol (1,2,3,4-cyclohexenetetrol) derivatives
            ##
            'HMC':{'stereo':'HMC'}, ## 5-hydroxymethyl-chonduritol
            'ACI':{'stereo':'ACI'}, ## 1-amino-5-hydroxymethyl-chonduritol
            ## Sialic acid (N-Acetylneuraminic acid, Neu5Ac, NANA)
            'SIA':{'stereo':'SIA','derivate':['SIA']}, ## (alpha)-sialic acid
            'SLB':{'stereo':'SIA','derivate':['SIA']}, ## beta-sialic acid
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

        self.d_symop = {
            '1mzy':[1,4,7,],
            '1n70':[9,6,3,],
            }

        self.maxrmsd = 2.75 ## 2fsy.pdb,2ft1.pdb ## try other combinations if above this rmsd
        self.maxrmsd_wrong = 9.5 ## 2eu1 ## assume error if above this rmsd

        self.minres = 5.0

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25
        self.max_mutations = 10

        self.path_pdb = '/oxygenase_local/data/pdb/'
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

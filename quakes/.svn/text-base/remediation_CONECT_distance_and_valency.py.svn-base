#!/software/bin/python
#$Id$
#Tommy Carstensen, University College Dublin, 2007

def main():

    import os, math, Numeric, sys

    pdbpath = '/oxygenase_local/data/pdb/'

    d_coordination = {
        'X' :1,
        'H' :1,
               'BE':1,                                                                      'B' :4,'C' :4,'N' :3,'O' :2,'F' :2,
        'NA':4,'MG':4,                                                                      'AL':4,       'P' :4,'S' :4,'CL':1,
        'K' :4,'CA':4,              'V' :4,'CR':1,'MN':4,'FE':4,'CO':4,'NI':4,'CU':4,'ZN':4,'GA':4,       'AS':1,'SE':1,'BR':1,
               'SR':4,'Y' :4,              'MO':4,       'RU':1,'RH':2,       'AG':2,'CD':4,                            'I' :1,
        'CS':4,'BA':4,              'TA':1,'W' :1,              'IR':1,'PT':4,'AU':2,'HG':4,'TL':2,'PB':4,
        'SM':4,
        'U' :2,
        }

    pdbs = []
    subdirs = os.listdir(pdbpath)
    subdirs.sort()
    for subdir in subdirs:
        if subdir < sys.argv[1]:
            continue
        print subdir
        files = os.listdir(pdbpath+subdir)
        files.sort()
        for file in files:
            pdb = file[3:7]
            if pdb in [
                ## CONECT errors
                '403d','406d','121d','127d','130d','264d','269d','1b16','2a4r',
                '1b2h','1b23','1b3j','1a61','2a73','1a72','1b46','2b5t','2c04',
                '3b6i','1c2w','1c2n','2b76','1b86','1c2t','2c3u','1c3q','2c3w',
                '2c3p','1b9j','1bei','3amv','2ca7','2bhd','1bhq','2bhc','2bha',
                '2bhj','2cdp','1d17','2cdz','2bh3','2bhb','1ce8','3cev','2d1g',
                '1ar1','1d21','1bkv','1d3s','3cgt','1aut','1ci7','2aw1','1cjk',
                '1d43','2bn7','1d44','1d45','1d46','2cl7','2clc','1ay6','1cml',
                '8bna','1cnr','1ayn','3ayk','1e08','1e0c','1e14','1cta','1bp8',
                '1e18','1cwj','1cwm','1bs1','1a08','821p','1a52','2c92','2c94',
                '2c97','2c9b','1ahv','1cwz','1cxe','1cxf','1cxh','1cxl','2bvr',
                '1e39','1d61','1e42','1e4k','1d7c','1d7d','5daa','2c2d','2c2h',
                '1dba','1e5f','1dbw','1c30','1c3o',
                ## CONECT errors involving carbon to metal ion
                '1cc1','2bj7','2bj8','2d5n','1e35','2e46',
                ## CONECT errors involving altlocs
                '2anq','2aot','2bhp','1d2u','1aym','2c9v','2ca4','2bw4','2bwi',
                '2by2','2by3',
                ## CONECT errors involving NMR models
                '2ccx',
                ## CONECT errors to self
                '2bhk','1bmp',
                ## different atoms with same coordinate (not altlocs; e.g. O1,O4 of sugar)
                '2axm',
                ## altlocs should have been used
                '2bky',
                ## bond shorter than minimum distance --> delete atoms
                '1c4k',
                ## delete atom 2790 and rename atom 2841 from HOA2 to O2A (or similar)
                '1c0i','2c03','1c0k','1c0l',
                ## rename atom 1140 in res_name EOH from C2 to oxygen
                '1cll',
                ## rename C4 to C2 (1125), C2 to C4 (1122); disconnect 1119,1122; 1125,1126
                '1cxu',
                ## rename 3359 H21 to O
                '1dba',
                ## CONECT error - rename atoms
                '2bvl',
                ## CONECT error?????
                '1brw',
                ## CONECT error??? carbon to metal
                '2al2','1ao0','2bs3','1au1',

##                ## ala-carbon-cu bond
##                '1ag0',
##                ## CONECT errors - water distance
##                '400d','200l',
##                ## SF4 iron sulfur cluster misconnection of diagonal elements (FE1,S4)
##                '1b0t','1b0j','1b0m','1b0v','1b0k','1b0p','1b0y',##'2a5h',
##                ## cyanide (cyn) reversal?
##                '1b0b','1aby',
##                ## cys connected to carbon of heme (correctly...)
##                '2b12',
                ## extreme but correct distance in NMR model
                '1aqs',
                ## extreme but correct distance in ligand
                '1d3f','1d3m','1d3n','1csc','2csc','3csc','4csc','4bp2','221l',
                '224l','336d','1cbn','1btu','1d6p','1d6q',
                ## extreme but correct distance to ligand
                '1a0r','1e55',
                ## extreme but correct distance involving altloc
                '445d','1chw',
                ## extre but correct distance involving beta-lactam ring
                '1b12',
                ## extreme but correct distance in low resolution structure
                '1c51',
                ## extreme but correct distance in epoxide
                '1b6a',

                ## nucleoside phosphate phsophor and not oxygen connected to metal ion
                '430d',
                ## carboxylate carbon and not oxygen connected to metal ion
                '2c6n','2bgg','2al1',
                
##                ## not an error of oxygen valence
##                '1a3w','1a3x',
##                ## ether o-o oxygen-oxygen bond
##                '2a7r','1aer','2b0c','1a80','2b8o',
##                ## valence number of boron exceeded
##                '1a3b',
##                ## misplaced atom causing large distance
##                '2akz',
##                ## correct trp-tyr bond in oxoferryl catalase-peroxidase
##                '2b2r',
##                ## correct weird wolfram connection
##                '1b25',
##                ## correct aluminium fluorid connection to magnesium
##                '1a6e',
##                ## asp/glu carbon connected to metal (mn,ca,mg) instead of oxygen
##                '1a47','1aa0','2afk','1afa','2afh','2afi','2c4f',
##                ## carbon-metal connection (not glu/asp)
##                '2b3x','2c2d','2c2h','1a9w','1aij','1aio','1c4c','1bbb','1aj9',
##                ## correct long lys-plp n-c bond
##                '1b4d',
##                ## incorrect hetero ligand id????
##                '2abb','3b51','3b52','3b53','2aba',

##                ## Fe-C connection
##                '1c4a', ## HC1
##                '1ajg','1do1','1do7','1cbm','1cg8', ## HEM-CMO
##                '1c51', ## SF4-UNK
##                '1biy','2bi4','7atj',
##                ## Zn-C
##                '2c6n','1au1',
##                ## Mg-C
##                '2al2','2al1', ## MG-PEP
##                '1ao0',
##                ## Co-C
##                '1cb7','1ccw','1bmt',
##                ## Se-Se connection
##                '2bc7','2bc8',
##                ## N-F
##                '2bef','1bjg',
                ]:
                continue
            fd = open('%s%s/%s' %(pdbpath,subdir,file), 'r')
            lines = fd.readlines()
            fd.close()
            d_atoms = {}
            for line in lines:
                record = line[:6].strip()
                if record in ['ATOM','HETATM',]:
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
                    coordinate = Numeric.array([x, y, z])
                    element = line[76:78].strip()
                    d_atoms[atom_no] = {'coordinate':coordinate,'atom_name':atom_name,'chain':chain,'res_no':res_no,'element':element,'altloc':altloc,'res_name':res_name,'iCode':iCode,}
                elif record == 'EXPDTA':
                    s_EXPDTA = line
                elif record == 'CONECT':
                    atom_no1 = int(line[6:11])
                    element1 = d_atoms[atom_no1]['element']
                    coordination = d_coordination[element1]
                    atom_nos = [int(line[i:i+5]) for i in range(11,len(line.strip()),5)]
                    res_no1 = d_atoms[atom_no1]['res_no']
                    res_name1 = d_atoms[atom_no1]['res_name']
                    atom_name1 = d_atoms[atom_no1]['atom_name']
                    altloc1 = d_atoms[atom_no1]['altloc']
                    chain1 = d_atoms[atom_no1]['chain']
                    iCode1 = d_atoms[atom_no1]['iCode']
                    if len(atom_nos) > coordination:
                        error = False
                        hydrogenbonds = 0
                        altlocs = 0
                        print d_atoms[atom_no1]
                        for i in range(len(atom_nos)):
                            atom_no2 = int(atom_nos[i])
                            res_name2 = d_atoms[atom_no2]['res_name']
                            element2 = d_atoms[atom_no2]['element']
                            res_no2 = d_atoms[atom_no2]['res_no']
                            altloc2 = d_atoms[atom_no2]['altloc']
                            chain2 = d_atoms[atom_no2]['chain']
                            if len(set([res_name1,res_name2])-set(['HOH'])) == 1:
                                hydrogenbonds += 1
                                print d_atoms[atom_no2]
                                pass
                            elif res_no1 != res_no2 and element1 == 'H' and element2 in ['O',]:
                                hydrogenbonds += 1
                                print d_atoms[atom_no2]
                                pass
                            elif len(set([altloc1,altloc2])) == 2:
                                altlocs += 1
                                pass
                            else:
                                print d_atoms[atom_no2]
                                error = True
                        if error == True and len(atom_nos)-hydrogenbonds-altlocs>coordination:
                            print hydrogenbonds
                            print d_atoms[atom_no1]
                            print len(atom_nos), coordination
                            print line
                            print file
                            stop_coordination
                    res_name1 = d_atoms[atom_no1]['res_name']
                    for atom_no2 in atom_nos:
                        atom_no2 = int(atom_no2)
                        element2 = d_atoms[atom_no2]['element']
                        res_no2 = d_atoms[atom_no2]['res_no']
                        res_name2 = d_atoms[atom_no2]['res_name']
                        atom_name2 = d_atoms[atom_no2]['atom_name']
                        altloc2 = d_atoms[atom_no2]['altloc']
                        chain2 = d_atoms[atom_no2]['chain']
                        iCode2 = d_atoms[atom_no2]['iCode']
                        if res_no1 == res_no2 and chain1 == chain2 and altloc1 != ' ' and altloc2 != ' ' and altloc1 != altloc2:
                            print atom_no1, d_atoms[atom_no1]
                            print atom_no2, d_atoms[atom_no2]
                            print pdb
                            stop_altlocs
                        if 'HOH' in set([res_name1,res_name2]):
                            maxlen = 6.096
                        ########################################################
                        elif set([element1,element2]) == set(['C','H']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.448
                            else:
                                maxlen = 1.499
                                print atom_no1, d_atoms[atom_no1]
                                print atom_no2, d_atoms[atom_no2]
                                print pdb
                                stop
                        elif set([element1,element2]) == set(['C','BE']):
                            maxlen = 1.93
                        elif set([element1,element2]) == set(['C','MG']):
                            maxlen = 2.07
                        elif set([element1,element2]) == set(['C','B']):
                            maxlen = 1.792 ## 1.56 ##1.741
                        elif set([element1,element2]) == set(['C','AL']):
                            maxlen = 2.24
                        elif set([element1,element2]) == set(['C','FE']):
                            if res_no1 == res_no2 and chain1 == chain2:
                                maxlen = 2.125
                            else:
                                maxlen = 2.368
                        elif set([element1,element2]) == set(['C']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 2.046
                            else:
                                maxlen = 1.779
                        elif set([element1,element2]) == set(['C','N']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.843 ##1.47-2.10
                            else:
                                maxlen = 2.045
                        elif set([element1,element2]) == set(['C','P']):
                            maxlen = 1.911
                        elif set([element1,element2]) == set(['C','AS']):
                            maxlen = 2.347 ##1.98 ##2.347
                        elif set([element1,element2]) == set(['C','O']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.843 ##-2.15 ## 2.546 ## epoxide/ethanol
                            else:
                                maxlen = 1.962
                        elif set([element1,element2]) == set(['C','S']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.988 ## -255
                            else:
                                maxlen = 2.991
                        elif set([element1,element2]) == set(['C','CR']):
                            maxlen = 1.92
                        elif set([element1,element2]) == set(['C','AU']):
                            maxlen = 1.742
                        elif set([element1,element2]) == set(['C','SE']):
                            maxlen = 2.464 ##1.98-2.71
                        elif set([element1,element2]) == set(['C','F']):
                            maxlen = 1.501 ## 1.34
                        elif set([element1,element2]) == set(['C','CL']):
                            maxlen = 2.034
                        elif set([element1,element2]) == set(['C','BR']):
                            maxlen = 2.807
                        elif set([element1,element2]) == set(['C','RU']):
                            maxlen = 2.257
                        elif set([element1,element2]) == set(['C','I']):
                            maxlen = 2.417  
                        elif set([element1,element2]) == set(['C','HG']):
                            maxlen = 3.134
                        elif set([element1,element2]) == set(['C','CO']):
                            maxlen = 3.372
                        elif set([element1,element2]) == set(['C','PT']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 3.080
                            else:
                                maxlen = 3.182
                        ########################################################
                        elif set([element1,element2]) == set(['O','H']):
                            if res_no1 == res_no2:
                                maxlen = 1.091
                            else:
                                maxlen = 1.918
                        elif set([element1,element2]) == set(['O','N']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.577
                            else:
                                maxlen = 2.863
                        elif set([element1,element2]) == set(['O','O']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.000
                            else:
                                maxlen = 1.000
                        elif set([element1,element2]) == set(['O','S']):
                            maxlen = 1.989
                        elif set([element1,element2]) == set(['O','SE']):
                            maxlen = 1.477
                        elif set([element1,element2]) == set(['O','B']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.673
                            else:
                                maxlen = 1.982
                        elif set([element1,element2]) == set(['O','P']):
                            maxlen = 2.133
                        elif set([element1,element2]) == set(['O','K']):
                            maxlen = 3.470
                        elif set([element1,element2]) == set(['O','NA']):
                            maxlen = 3.073
                        elif set([element1,element2]) == set(['O','F']):
                            maxlen = 2.389
                        elif set([element1,element2]) == set(['O','BR']):
                            maxlen = 2.060
                        elif set([element1,element2]) == set(['O','V']):
                            maxlen = 4.596
                        elif set([element1,element2]) == set(['O','ZN']):
                            maxlen = 4.453
                        elif set([element1,element2]) == set(['O','CA']):
                            maxlen = 4.839
                        elif set([element1,element2]) == set(['O','FE']):
                            maxlen = 4.167
                        elif set([element1,element2]) == set(['O','BE']):
                            maxlen = 1.934
                        elif set([element1,element2]) == set(['O','CO']):
                            maxlen = 3.501
                        elif set([element1,element2]) == set(['O','MN']):
                            maxlen = 3.104
                        elif set([element1,element2]) == set(['O','HG']):
                            maxlen = 3.024
                        elif set([element1,element2]) == set(['O','MO']):
                            maxlen = 2.225
                        elif set([element1,element2]) == set(['O','CD']):
                            maxlen = 2.728
                        elif set([element1,element2]) == set(['O','AS']):
                            maxlen = 2.628
                        elif set([element1,element2]) == set(['O','PB']):
                            maxlen = 3.032
                        elif set([element1,element2]) == set(['O','CU']):
                            maxlen = 2.653
                        elif set([element1,element2]) == set(['O','NI']):
                            maxlen = 2.746
                        elif set([element1,element2]) == set(['O','MG']):
                            maxlen = 4.030
                        elif set([element1,element2]) == set(['O','W']):
                            maxlen = 3.676
                        elif set([element1,element2]) == set(['O','AL']):
                            maxlen = 3.444
                        elif set([element1,element2]) == set(['O','CS']):
                            maxlen = 3.543
                        elif set([element1,element2]) == set(['O','BA']):
                            maxlen = 3.132
                        elif set([element1,element2]) == set(['O','Y']):
                            maxlen = 2.683
                        elif set([element1,element2]) == set(['O','U']):
                            maxlen = 2.397
                        elif set([element1,element2]) == set(['O','SM']):
                            maxlen = 2.967
                        elif set([element1,element2]) == set(['O','TL']):
                            maxlen = 3.013
                        elif set([element1,element2]) == set(['O','SR']):
                            maxlen = 2.955
                        elif set([element1,element2]) == set(['O','PT']):
                            maxlen = 1.651
                        elif set([element1,element2]) == set(['O','CR']):
                            maxlen = 1.978
##                        ########################################################
                        elif set([element1,element2]) == set(['N','H']):
                            maxlen = 1.104
                        elif set([element1,element2]) == set(['N','S']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.906
                            else:
                                maxlen = 2.051
                        elif set([element1,element2]) == set(['N','N']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.6+5
                            else:
                                maxlen = 1.000
                                print atom_no1, d_atoms[atom_no1]
                                print atom_no2, d_atoms[atom_no2]
                                print pdb
                                stop
                        elif set([element1,element2]) == set(['N','P']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.792
                            else:
                                maxlen = 1.638
                        elif set([element1,element2]) == set(['N','B']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 1.446
                            else:
                                maxlen = 1.470
                        elif set([element1,element2]) == set(['N','FE']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 2.294
                            else:
                                maxlen = 2.912
                        elif set([element1,element2]) == set(['N','RU']):
                            maxlen = 2.145
                        elif set([element1,element2]) == set(['N','MG']):
                            maxlen = 2.860
                        elif set([element1,element2]) == set(['N','PT']):
                            maxlen = 2.920
                        elif set([element1,element2]) == set(['N','NA']):
                            maxlen = 2.999
                        elif set([element1,element2]) == set(['N','F']):
                            if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                maxlen = 2.089
                            else:
                                maxlen = 2.399
                        elif set([element1,element2]) == set(['N','CL']):
                            maxlen = 1.509
                        elif set([element1,element2]) == set(['N','K']):
                            maxlen = 3.158
                        elif set([element1,element2]) == set(['N','ZN']):
                            maxlen = 2.946
                        elif set([element1,element2]) == set(['N','CA']):
                            maxlen = 3.071
                        elif set([element1,element2]) == set(['N','CD']):
                            maxlen = 2.527
                        elif set([element1,element2]) == set(['N','CU']):
                            maxlen = 2.570
                        elif set([element1,element2]) == set(['N','AU']):
                            maxlen = 2.065
                        elif set([element1,element2]) == set(['N','NI']):
                            maxlen = 2.616
                        elif set([element1,element2]) == set(['N','MN']):
                            maxlen = 2.768
                        elif set([element1,element2]) == set(['N','MO']):
                            maxlen = 2.370
                        elif set([element1,element2]) == set(['N','AL']):
                            maxlen = 2.136
                        elif set([element1,element2]) == set(['N','CO']):
                            maxlen = 2.681
                        elif set([element1,element2]) == set(['N','RH']):
                            maxlen = 2.135
                        elif set([element1,element2]) == set(['N','IR']):
                            maxlen = 2.318
                        elif set([element1,element2]) == set(['N','HG']):
                            maxlen = 2.739
                        ########################################################
                        elif set([element1,element2]) == set(['S']):
                            if res_no1 == res_no2 and chain1 == chain2:
                                maxlen = 2.337
                            else:
                                maxlen = 3.453
                        elif set([element1,element2]) == set(['S','X']):
                            maxlen = 2.015
                        elif set([element1,element2]) == set(['S','P']):
                            maxlen = 2.122
                        elif set([element1,element2]) == set(['S','ZN']):
                            maxlen = 3.440
                        elif set([element1,element2]) == set(['S','W']):
                            maxlen = 2.526
                        elif set([element1,element2]) == set(['S','MO']):
                            maxlen = 2.720
                        elif set([element1,element2]) == set(['S','CO']):
                            maxlen = 2.379
                        elif set([element1,element2]) == set(['S','NI']):
                            maxlen = 2.621
                        elif set([element1,element2]) == set(['S','H']):
                            maxlen = 1.331
                        elif set([element1,element2]) == set(['S','F']):
                            maxlen = 1.561
                        elif set([element1,element2]) == set(['S','AS']):
                            maxlen = 2.448
                        elif set([element1,element2]) == set(['S','PB']):
                            maxlen = 1.642
                        elif set([element1,element2]) == set(['S','FE']):
                            if res_name1 == 'SF4' and set([atom_name1,atom_name2,]) in [
                                set(['FE1','S4',]),set(['FE2','S3',]),set(['FE3','S1',]),set(['FE4','S2',]),
                                set(['FE1','S2',]),set(['FE3','S4',]),set(['FE3','S2',]),set(['FE4','S1',]),
                                ]:
                                maxlen = 4.332
                            else:
                                maxlen = 2.904
                        elif set([element1,element2]) == set(['S','CU']):
                            maxlen = 3.099
                        elif set([element1,element2]) == set(['S','CD']):
                            maxlen = 2.551
                        elif set([element1,element2]) == set(['S','AU']):
                            maxlen = 2.502
                        elif set([element1,element2]) == set(['S','AG']):
                            maxlen = 2.722
                        elif set([element1,element2]) == set(['S','HG']):
                            maxlen = 3.130
                        elif set([element1,element2]) == set(['S','PT']):
                            maxlen = 1.972
                        elif set([element1,element2]) == set(['S','GA']):
                            maxlen = 2.371
                        ########################################################
                        elif set([element1,element2]) == set(['P','AU']):
                            maxlen = 2.269
                        ########################################################
                        elif 'C' in [element1,element2] and res_no1 == res_no2:
                            print d_atoms[atom_no1]
                            print d_atoms[atom_no2]
                            print pdb
                            stopC
                        elif 'O' in [element1,element2] and res_no1 == res_no2:
                            print d_atoms[atom_no1]
                            print d_atoms[atom_no2]
                            print pdb
                            stopO
                        elif 'H' in [element1,element2] and res_no1 == res_no2:
                            print d_atoms[atom_no1]
                            print d_atoms[atom_no2]
                            stopH
                        elif 'N' in [element1,element2] and res_no1 == res_no2:
                            print d_atoms[atom_no1]
                            print d_atoms[atom_no2]
                            stopN
                        elif 'P' in [element1,element2] and res_no1 == res_no2:
                            print d_atoms[atom_no1]
                            print d_atoms[atom_no2]
                            stopP
                        elif 'S' in [element1,element2] and res_no1 == res_no2:
                            print d_atoms[atom_no1]
                            print d_atoms[atom_no2]
                            stopS
                        ########################################################
                        elif set([element1,element2]) == set(['F','MG']):
                            maxlen = 2.521
                        elif set([element1,element2]) == set(['F','AL']):
                            maxlen = 1.844
                        elif set([element1,element2]) == set(['F','BE']):
                            maxlen = 1.710
                        elif set([element1,element2]) == set(['F','FE']):
                            maxlen = 1.891
                        ########################################################
                        elif set([element1,element2]) == set(['BR','TA']):
                            maxlen = 2.673
                        ########################################################
                        elif set([element1,element2]) == set(['MG','AL']):
                            maxlen = 3.470
                        elif set([element1,element2]) == set(['NI','FE']):
                            maxlen = 2.968
                        ########################################################
                        elif set([element1,element2]) == set(['X','X']):
                            maxlen = 1.797
                        elif set([element1,element2]) == set(['CS','CS']):
                            maxlen = 1.931
                        elif set([element1,element2]) == set(['MG','MG']):
                            maxlen = 2.611
                        elif set([element1,element2]) == set(['HG','HG']):
                            maxlen = 1.841
                        elif set([element1,element2]) == set(['CU','CU']):
                            maxlen = 2.514
                        elif set([element1,element2]) == set(['B','B']):
                            maxlen = 1.790
                        elif set([element1,element2]) == set(['SE','SE']):
                            maxlen = 2.323
                        elif set([element1,element2]) == set(['TA','TA']):
                            maxlen = 4.087
                        elif set([element1,element2]) == set(['FE','FE']):
                            maxlen = 3.172
                        ########################################################
                        elif set([element1,element2]) == set(['CL','NA']):
                            maxlen = 2.765
                        elif set([element1,element2]) == set(['CL','FE']):
                            maxlen = 2.383
                        elif set([element1,element2]) == set(['CL','PT']):
                            maxlen = 2.547
                        elif set([element1,element2]) == set(['CL','NI']):
                            maxlen = 2.362
                        elif set([element1,element2]) == set(['CL','ZN']):
                            maxlen = 2.129
                        elif set([element1,element2]) == set(['CL','SR']):
                            maxlen = 2.594
                        elif set([element1,element2]) == set(['CL','CU']):
                            maxlen = 3.108
                        ########################################################
                        else:
                            maxlen = 1.6
                        dist = math.sqrt(sum((d_atoms[atom_no1]['coordinate']-d_atoms[atom_no2]['coordinate'])**2))
                        if dist < 0.140:
                            print dist
                            print atom_no1, d_atoms[atom_no1]
                            print atom_no2, d_atoms[atom_no2]
                            print line
                            print file
                            stop_minlen
                        if dist > maxlen:
                            if element1 == 'O' and element2 == 'O':
                                fd = open('ether.txt','r')
                                s = fd.read()
                                fd.close()
                                if pdb not in s:
                                    fd = open('ether.txt','a')
                                    s = fd.write('%s,' %(pdb))
                                    fd.close()
                                continue
                            if len(set([element1,element2,])&set(['C',])) == 1 and len(set([element1,element2,])&set(['CA','MG','MN','ZN',])) == 1 and len(set([res_name1,res_name2,])&set(['ASP','GLU',])) == 1 and len(set([res_name1,res_name2,])&set(['CA','MG','MN','ZN',])) == 1:
                                fd = open('aspglu_metal_bonds.txt','r')
                                s = fd.read()
                                fd.close()
                                if pdb not in s:
                                    fd = open('aspglu_metal_bonds.txt','a')
                                    s = fd.write('%s,' %(pdb))
                                    fd.close()
                                continue
                            print
                            print file, maxlen, dist
                            print atom_no1, d_atoms[atom_no1]
                            print atom_no2, d_atoms[atom_no2]
                            print s_EXPDTA
                            stop_distance

    d_hydrogen_equivalents = {
        'C' :2,'O' :0,'H': 0,
        'F' :0,'CL':0,'BR':0,'I' :0,
        'N' :1,'P' :1,
        }


if __name__ == '__main__':
    main()

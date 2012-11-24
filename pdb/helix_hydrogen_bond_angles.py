def main():

    import os, numpy, math
    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import parse_pdb

    path = '/oxygenase_local/data/pdb/'
    dirs = os.listdir(path)
    dirs.sort()
    for dir in dirs:
        if dir < sys.argv[1][1:3]:
            continue
        print dir
        ents = os.listdir('%s%s' %(path,dir))
        ents.sort()
        for ent in ents:
            pdb = ent[3:7]
            if dir == sys.argv[-1][1:3] and pdb < sys.argv[-1]:
                continue
            print pdb
            if pdb in [
                '1a0k', ## 310 helix?
                '2aeb', ## helix kink
                '2aoc','2aod','2aog', ## some error...
                '2c03','2c04','3bxd', ## bent helix
                '4fbp','5fbp', ## pro residue in i+6 position
                '1a2f','1a2g', ## gly residue in i+1,i+3 position (310 helix?)
                ]:
                continue
            if pdb in [ ## from phipsiparseall
                ## std_res_name in neither REMARK465 nor SEQRES (large peptide length)
                '1a3q','1c3w','1c3x','1ad5',
                ## hetID in SEQRES and HETATM    but not MODRES
                '2b2u','2a2x','2b7f','2c2k','2c2m','2c2o','2c2z','2aal','1bdu',
                '2ag3','2age','2agg','1an5','2ci1','2ank',
##                ## hetID in SEQRES and REMARK465 but not MODRES
##                '1aco','7acn','8acn','2aig','3aig','1g1f','1g1s','2f4i',
                ## hetID in MODRES               but not SEQRES
                '2a4o','3bbd','1orw',
                ## std_res in REMARK465 but not SEQRES
                '2uva','2uva','2uvc',
                ## hetID in SEQRES but not REMARK465,ATOM
                '2vhn',
                ## REMARK465 records missing
                '1c04','3b9v','1fka','1deq','2jcc','2jj4','1i3q','1i50','1i6h',
                '1iw7','1p0t','1smy','1uf2','2v0z','2v7n','2vs4','1zbb','1zlv',
                ## remark470 records missing (alpha carbon only in most cases)
                '1cc0','1f1o','1ffk','3b5d','3b5w','3b5x','3b5y','3b5z','3b61',
                '3b62','1bdx','1d3l','1gix','1e8s','1j5a','4cro','1i9w','1giy',
                '2v9l',

                ## SEQRES/ATOM conflict (incorrect res_name)
                '2c38','2j01','2j03','4icd','2aew','2plv','1d2q',

                ## remark465 initiation line missing
                '2iwq','2v07',

                ## incorrect residue numbers (N-terminal)
                '1ef0','2ji5','1eu3',
                ## incorrect residue numbers (reversed)
                '9lpr',
                ## incorrect residue numbers (insertion)
                '1a4k','2ged','1cj0',
                ## incorrect residue numbers (altloc)
                '2bb3',
                ## incorrect residue numbers (C-terminal)
                '1a7l','1ce0','2ci0','2cib','1f32',
                ## inccorect residue numbers (type error)
                '1f7o','1f7p','1ke8','1kj4','1n3f','1nlq','1nv8','1nv9','1pxx',
                '2by6','1ll0','2bfk','2gnk','1e5r','1lox','1ca8','1x11',
                ## incorrect residue numbers (REMARK465)
                '1jn6',
                ## incorrect residue numbers (ATOM)
                '2qqh',
                ## incorrect residue numbers ("reset"/"broken"... 28,29,30,1,2,3)
                '1jpl','1juq','1k4j','1q6j','1q6m','2oxg','2oxr','2qzl',
                '1swf','1swg','2z5s','1yr6','1yr7','1yr8','1yr9',
                ## incorrect residue numbers (iCodes)
                '1iao','2iad','1dki','1es0','7pck','1pfz','1qdm','1ygp',
                '2qri','2qrs','2qrt','3dgv',
                ## incorrect residue numbers (none)
                '2om7',

                ## incorrect chain ID (ATOM)
                '1bml','1h74','1jwt',
                ## incorrect iCode (ATOM)
                '2ius',
                ## incorrect SEQRES sequence
                '1clw',

                ## identical ATOM/HETATM IDs (coordinate section)
                '1h9h','427d','2olb','11gs','121p','12gs','13gs','16gs','16pk','17gs','185d','18gs','193d','19gs','1a05','1a0i','1a25','1a3l','1a44','1a48','1a4a','1a4b','1a4c','1a4f','1a4g','1a4k','1a4q','1a52','1a5a','1a5b','1a5z','1a65','1a69','1a6g','1a6m','1a6q','1a71','1a78','1a79','1a8i', '1a8s', '1a8u', '1a9c', '1a9m', '1a9x', '1a9y', '1a9z', '1aal', '1aax', '1aaz', '1aba', '1abw', '1aby', '1ad8', '1ad9', '1adb', '1adc', '1adf', '1adg', '1adj', '1aec', '1af6', '1afa', '1afb', '1ag0', '1ah8', '1ahe', '1ahf', '1ahx', '1ahy', '1ai4', '1ai5', '1ai6', '1ai7', '1ai8', '1aix', '1aj6', '1aj9', '1ajn', '1ajp', '1ajq', '1aku', '1akv', '1all', '1ami', '1ao5', '1aok', '1aor', '1apm', '1aq1', '1aq6', '1aq7', '1ar1', '1ari', '1art', '1arx', '1asm', '1asn', '1aso', '1asp', '1asq', '1at1', '1atg', '1atl', '1atn', '1atp', '1au1', '1aua', '1auj', '1aus', '1av4', '1av6', '1avb', '1avf', '1avq', '1awb', '1axa', '1axd', '1axs', '1aya', '1ayo', '1ayr', '1azr', '1azs', '1b08', '1b0m', '1b25', '1b4k', '1b4n', '1b55', '1b6h', '1b8t', '1b92', '1b9f', '1bb1', '1bch', '1bcj', '1bcr', '1bcs', '1bd0', '1be3', '1beh', '1ben', '1bfn', '1bgy', '1bh3', '1bij', '1biq', '1biz', '1bj3', '1bja', '1bk5', '1bks', '1bnl', '1bow', '1bps', '1bqa', '1bqd', '1bqh', '1brh', '1brk', '1brt', '1bsk', '1btc', '1bvc', '1bvd', '1bwn', '1bwo', '1cdk', '1cdm', '1cel', '1ces', '1cf5', '1cgk', '1ch0', '1ch4', '1cls', '1cml', '1cnt', '1coh', '1con', '1cow', '1cpc', '1cpq', '1cpr', '1crx', '1ctf', '1ctp', '1cxe', '1cxf', '1cxh', '1czf', '1d33', '1d35', '1d39', '1d40', '1d61', '1d8h', '1daj', '1dan', '1dbj', '1dbk', '1dbp', '1dbr', '1dff', '1dhy', '1dif', '1djc', '1dkk', '1dl4', '1dlr', '1dls', '1dmr', '1dms', '1dpe', '1dpm', '1drf', '1drj', '1drk', '1dut', '1eas', '1ecc', '1ecg', '1efg', '1eg6', '1elp', '1ent', '1epm', '1epn', '1epp', '1epq', '1eta', '1etb', '1eth', '1etj', '1fax', '1fbt', '1fdh', '1fgh', '1fgj', '1fig', '1fiv', '1fkb', '1fki', '1fn1', '1fpc', '1fpi', '1fpj', '1fpl', '1fsa', '1fui', '1fuo', '1fup', '1fuq', '1fur', '1fv7', '1fyl', '1g5l', '1gac', '1gaf', '1gan', '1gbu', '1gdi', '1gic', '1gj2', '1glc', '1gld', '1gle', '1gli', '1gmk', '1gmp', '1gmr', '1gnm', '1gnn', '1gno', '1gnw', '1gpa', '1gpd', '1gqh', '1gsu', '1gsy', '1gti', '1gtm', '1gtp', '1hab', '1hac', '1hah', '1hbg', '1hbi', '1hbs', '1hbt', '1hcn', '1hcs', '1hct',  '1hds', '1hdt', '1hfp', '1hfq', '1hfr', '1hga', '1hgb', '1hgc', '1hgt', '1hho', '1hih', '1hii', '1hiv', '1hmd', '1hmo', '1hpc', '1hr3', '1hro', '1hrp', '1htp', '1huj', '1huk', '1hur', '1hvq', '1hvr', '1hxp', '1hya', '1i3t', '1i47', '1i7v', '1iak', '1ibg', '1icj', '1idc', '1ide', '1iea', '1ieb', '1igt', '1igy', '1ils', '1ilu', '1imc', '1imd', '1ipw', '1irn', '1iro', '1ith', '1ixx', '1izb', '1j8l', '1jaw', '1jdb', '1jfd', '1jka', '1jkb', '1jkc', '1jkd', '1jlx', '1jpc', '1jsa', '1jst', '1jsu', '1kao', '1kpe', '1kr7', '1krb', '1krc', '1krn', '1ksa', '1kvq', '1kvr', '1kvs', '1kvt', '1kvu', '1lam', '1lcj', '1lck', '1lcp', '1ldp', '1len', '1lia', '1lkk', '1lkl', '1llo', '1lml', '1loa', '1lob', '1lpa', '1lpb', '1lth', '1ltr', '1lu2', '1lya', '1lyb', '1lyw', '1m5k', '1mae', '1maf', '1mbd', '1mfr', '1mhe', '1mhk', '1mmb', '1mmo', '1mmp', '1mmq', '1mmr', '1mpa', '1mpg', '1mpq', '1mpr', '1mwe', '1myf', '1myh', '1myj', '1nah', '1nai', '1nas', '1nbb', '1nbc', '1nbe', '1nci', '1nco', '1nec', '1nfd', '1nfp', '1nlk', '1nn2', '1nnc', '1np1', '1nrs', '1nsa', '1nsc', '1nsd', '1nzr', '1oaa', '1oat', '1ohj', '1olc', '1ord', '1otg', '1ouu', '1ova', '1ovw', '1p35', '1pag', '1pam', '1pbo', '1pea', '1pgt', '1php', '1ply', '1pnk', '1pnl', '1pnm', '1poi', '1ppm', '1ppr', '1psc', '1psh', '1ptv', '1pty', '1pyt', '1qbi', '1qd6', '1qdc', '1qf8', '1qge', '1qha', '1qhf', '1qpr', '1qrd', '1qs8', '1qsg', '1raa', '1rab', '1rac', '1rad', '1rae', '1raf', '1rag', '1rah', '1rai', '1rbl', '1rcm', '1rcp', '1rdg', '1rdi', '1rdj', '1rdk', '1rdl', '1rdm', '1rdn', '1rdv', '1rdx', '1rdy', '1rdz', '1rem', '1rge', '1rgf', '1rgg', '1rgh', '1rn1', '1rnc', '1rsy', '1rtf', '1rth', '1rti', '1rtj', '1rtm', '1rvw', '1sac', '1scu', '1sda', '1sdk', '1sdl', '1sep', '1sfi', '1sft', '1sgc', '1shb', '1shd', '1skj', '1sky', '1slb', '1slc', '1sli', '1slt', '1slu', '1smp', '1sps', '1stc', '1sty', '1swd', '1swe', '1swn', '1swp', '1swr', '1taq', '1tar', '1tat', '1taw', '1tcb', '1tcf', '1tco', '1tei', '1tet', '1thb', '1ths', '1tlc', '1tlg', '1tli', '1tlp', '1tmb', '1tmn', '1tmu', '1tn4', '1tpf', '1tpk', '1try', '1trz', '1tsd', '1ttp', '1ttq', '1tub', '1tyl', '1tym', '1tyu', '1tyw', '1tyx', '1udg', '1udh', '1uma', '1uz8', '1vdr', '1vrt', '1vru', '1vwt', '1wav', '1wgi', '1wgj', '1whs', '1wht', '1wyk', '1xel', '1xgm', '1xgn', '1xgs', '1xik', '1xso', '1yag', '1yec', '1yef', '1yeg', '1yeh', '1yrq', '1zeg', '1zeh', '1zni', '1znj', '20gs', '211d', '219d', '21gs', '221p', '22gs', '253d', '256b', '258d', '25c8', '277d', '2a3h', '2aaa', '2aac', '2aae', '2ahj', '2arc', '2at1', '2atc', '2ay1', '2ay2', '2ay3', '2ay4', '2ay5', '2ay6', '2ay7', '2ay8', '2ay9', '2azu', '2btf', '2bvw', '2c7e', '2cah', '2cht', '2cmm', '2cpk', '2cst', '2d34', '2d95', '2da8', '2dhb', '2dhn', '2dmr', '2dri', '2ecp', '2fbj', '2fus', '2gli', '2glr', '2gls', '2gsa', '2gss', '2hck', '2hhb', '2hhd', '2hk6', '2hmq', '2hmz', '2hr7', '2iep', '2kau', '2lal', '2lhb', '2mas', '2mcp', '2mhb', '2mpa', '2mpr', '2msb', '2nac', '2np1','2otc', '2ovw', '2oxi', '2pax', '2pcp', '2pfl', '2pgh', '2pgt', '2phk', '2pk4', '2pld', '2ple', '2prc', '2prg', '2q41', '2qwa', '2qwb', '2qwc', '2qwd', '2qwe', '2qwf', '2qwg', '2qwh', '2qwi', '2qwj', '2qwk', '2r2f', '2ran', '2rkm', '2rus', '2scp', '2sec', '2shk', '2sli', '2sn3', '2sni', '2sod', '2taa', '2tgd', '2tli', '2tmn', '2tn4', '2trc', '2tsc', '2wea', '2web', '2wec', '308d', '367d', '380d', '381d', '382d', '3a3h', '3at1', '3azu', '3chb', '3cyt', '3dmr', '3eng', '3fru', '3fx2', '3grs', '3gss', '3gst', '3hat', '3hhb', '3ljr', '3lkf', '3np1', '3p2p', '3pax', '3pca', '3pcb', '3pcc', '3pcd', '3pce', '3pcf', '3pcg', '3pch', '3pci', '3pcj', '3pck', '3pcn', '3pmg', '3prc', '3prn', '3pro', '3sc2', '3sli', '3sqc', '3tli', '3tmn', '421p', '455d', '456c', '4a3h', '4at1', '4cts', '4dfr', '4dmr', '4eng', '4fx2', '4gr1', '4gsa', '4gss', '4hhb', '4hmg', '4mdh', '4np1', '4ovw', '4pax', '4rub', '4sli', '4tgl', '4tmn', '5at1', '5bir', '5dnb', '5fx2', '5gss', '5hmg', '5mdh', '5prn', '5tmn', '621p', '6abp', '6adh', '6at1', '6gsp', '6gss', '6gsu', '6gsv', '6gsw', '6gsx', '6ins', '6prn', '6tmn', '721p', '7aat', '7abp', '7at1', '7gss', '7odc', '7prn', '7taa', '830c', '8aat', '8abp', '8at1', '8gss', '8prn', '966c', '9gss', '9pap', '9rub','2q44',
                ## identical IDs in REMARK465 and ATOM records
                '3c9u',

                ## altloc not in MODRES records
                '354d',
                ]:
                continue            
            file = '%s%s/%s' %(path,dir,ent)
            fd = open(file,'r')
            lines = fd.readlines()
            fd.close()

            ##
            ## parse header
            ##
            d_header = parse_pdb.parse_header(lines)
            if not 'HELIX' in d_header.keys():
                continue
            if d_header['EXPDTA'] != 'X-RAY':
                continue

            ##
            ## parse coordinates
            ##
            d_coordinates, d_ATOMseq = parse_pdb.parse_coordinates(lines, d_header,)

            for chain in d_header['SEQRES']['chains']:
                if chain not in d_header['HELIX'].keys():
                    continue
                for helix_no in d_header['HELIX'][chain].keys():
                    res_no = d_header['HELIX'][chain][helix_no]['res_no']
                    iCode = d_header['HELIX'][chain][helix_no]['iCode']
                    helix_len = d_header['HELIX'][chain][helix_no]['len']
                    if not res_no in d_coordinates['chains'][chain]['residues'].keys():
                        continue
                    seq_no = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['seq_no']
                    if seq_no == -1:
                        continue ## hoeker-loesning
                    if seq_no >= len(d_ATOMseq[chain]['res_nos']):
                        continue ## hoeker-loesning
                    if d_ATOMseq[chain]['res_nos'][seq_no] != res_no:
                        stop
                    if d_ATOMseq[chain]['iCodes'][seq_no] != iCode:
                        stop
##                    if d_ATOMseq[chain]['altlocs'][seq_no] != ' ':
##                        stop
                    if helix_len < 11:
                        continue
                    l_angles = []
                    for i in range(seq_no+1,min(seq_no+helix_len,len(d_ATOMseq[chain]['res_nos']))-4-1):
                        res_no1 = d_ATOMseq[chain]['res_nos'][i]
                        res_no2 = d_ATOMseq[chain]['res_nos'][i+4]
                        iCode1 = d_ATOMseq[chain]['iCodes'][i]
                        iCode2 = d_ATOMseq[chain]['iCodes'][i+4]
                        record1 = d_ATOMseq[chain]['records'][i]
                        record2 = d_ATOMseq[chain]['records'][i+4]

                        if not res_no1 in d_coordinates['chains'][chain]['residues'].keys() and record1 == 'REMARK465':
                            continue
                        if not res_no2 in d_coordinates['chains'][chain]['residues'].keys() and record2 == 'REMARK465':
                            continue

                        res_name1 = d_coordinates['chains'][chain]['residues'][res_no1]['d_iCodes'][iCode1]['res_name']

                        if not 'O' in d_coordinates['chains'][chain]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'].keys():
                            continue
                        if not 'H' in d_coordinates['chains'][chain]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'].keys():
                            continue

                        O = d_coordinates['chains'][chain]['residues'][res_no1]['d_iCodes'][iCode1]['atoms']['O']['coordinate']
                        N = d_coordinates['chains'][chain]['residues'][res_no2]['d_iCodes'][iCode2]['atoms']['N']['coordinate']
                        H = d_coordinates['chains'][chain]['residues'][res_no2]['d_iCodes'][iCode2]['atoms']['H']['coordinate']
                        vNH = N-H
                        vOH = O-H
                        angle = (180/math.pi)*math.acos(
                            numpy.dot(vNH,vOH)/math.sqrt(
                                sum(vNH**2)*sum(vOH**2)
                                )
                            )
                        print pdb, chain, res_no1, res_no2, angle
                        if res_name1 != 'GLY' and angle < 78: ## 113
                            print res_name1
                            stop1
                        if res_name1 == 'GLY' and angle < 101:
                            stop1
                        l_angles += [angle]
                    if len(l_angles) > 0:
                        print helix_len, sum(l_angles)/len(l_angles)
                        if sum(l_angles)/len(l_angles) > 171:
                            stop3
                        if sum(l_angles)/len(l_angles) < 149:
                            stop4

                        
##            d_helix = {}
##            l_modres = []
##            d_angles = {}
##            for line in lines:
##                record = line[:6].strip()
##                if record in ['ATOM','HETATM',]:
##                    if record == 'HETATM' and line[21:27] not in l_modres:
##                        continue
##                    chain = line[21]
##                    if chain not in d_helix.keys():
##                        continue
##                    res_no = int(line[22:26])
##                    if res_no not in d_helix[chain]:
##                        continue
##                    iCode = line[26]
##                    if iCode != ' ':
##                        print line
##                        stop
##                    atom_name = line[12:16].strip()
##                    x = float(line[30:38])
##                    y = float(line[38:46])
##                    z = float(line[46:54])
##                    coord = [numpy.array([x,y,z,])]
##                    if res_no+4 in d_helix[chain] and atom_name == 'O':
##                        d_helix[chain][res_no][atom_name] = coord
##                    if res_no-4 in d_helix[chain] and atom_name in ['N','H',]:
##                        d_helix[chain][res_no][atom_name] = coord
##                        if atom_name == 'H':
##                            N = d_helix[chain][res_no]['H']
##                            H = d_helix[chain][res_no]['N']
##                            O = d_helix[chain][res_no-4]['H']
##                            vNH = N-H
##                            vOH = O-H
##                            angle = math.acos(numpy.dot(vNH,vOH))
##                            d_angles[res_no] = angle
##                            print angle
##                elif record == 'HELIX':
##                    chain = line[19]
##                    res_no1 = int(line[21:25])
##                    res_no2 = int(line[33:37])
##                    if chain not in d_helix.keys():
##                        d_helix[chain] = {}
##                    for res_no in range(res_no1,res_no2+1):
##                        d_helix[chain][res_no] = {}
##                    if line[19] != line[31]:
##                        print line
##                        stop
##                elif record == 'MODRES':
##                    chain = line[16]
##                    res_no = int(line[18:22])
##                    iCode = line[22]
##                    l_modres += ['%s%4i%s' %(chain,res_no,iCode,)]
##
##            if d_angles != {}:
##                print d_helix
##                print d_angles
##                print file
##                stop


if __name__ == '__main__':
    main()

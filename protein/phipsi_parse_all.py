import os,sys,math

path_pdb = '/oxygenase_local/data/pdb/'

d_pep = {
    'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
    'UNK':'X','ASX':'X','GLX':'X',
    'MSE':'M', ## ligand in 2e1a
    }

l_nuc = [
     'A', 'C', 'G', 'U',      'I', ## ribonucleotides
    'DA','DC','DG',     'DT','DI', ## deoxyribonucleotides
    'N', ## N is any 5'-monophosphate nucleotide
    ]

s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def main():

##    pdbs = []
##    files = [
##        'phi_mit_gap.txt','psi_mit_gap.txt',
##        'phi_ohne_gap.txt','psi_ohne_gap.txt'
##        ]
##    for file in files:
##        fd = open(file,'r')
##        lines = fd.readlines()
##        fd.close()
##        for line in lines:
##            pdbs += [line.split()[1]]
##        os.remove(file)
##    pdbs = list(set(pdbs))
##    pdbs.sort()
##    print pdbs
##    pdbs = ['12ca', '1am5', '1aov', '1aro', '1b35', '1b76', '1bgw', '1bkc', '1bln', '1bxw', '1c12', '1c88', '1c8s', '1ca3', '1cm7', '1cz8', '1d2s', '1d5a', '1dca', '1deu', '1dfi', '1dfl', '1dts', '1dvf', '1dvp', '1dxl', '1e0k', '1e94', '1eb1', '1ecx', '1esb', '1etq', '1exz', '1fe8', '1fem', '1fen', '1fh5', '1fhw', '1fi8', '1flr', '1fn4', '1frt', '1fwu', '1g9x', '1ggm', '1gku', '1gl9', '1gqo', '1gt9', '1gxd', '1gyg', '1h9d', '1haw', '1hca', '1hdm', '1hdq', '1hdu', '1hea', '1hec', '1hee', '1hgu', '1hj9', '1hxj', '1i9y', '1i9z', '1im3', '1j4a', '1j4s', '1jio', '1jip', '1jl0', '1jl2', '1kek', '1kln', '1krl', '1ktk', '1ldt', '1lik', '1ll1', '1lnl', '1lqp', '1m7d', '1m7i', '1m7z', '1mdp', '1mpt', '1mvm', '1myp', '1o6o', '1o6v', '1o6y', '1oau', '1oay', '1ohf', '1om3', '1op3', '1op5', '1ot7', '1pfg', '1pj8', '1pnx', '1q15', '1q19', '1qbz', '1qe5', '1qgj', '1qmb', '1qub', '1qun', '1r4m', '1r4n', '1rba', '1rf5', '1rlr', '1ru7', '1rvx', '1rvz', '1rxt', '1s7c', '1scd', '1sid', '1sie', '1sxi', '1t71', '1tau', '1tmf', '1tr1', '1trr', '1tsr', '1tt5', '1tyv', '1ua0', '1un0', '1uwy', '1uxl', '1v00', '1v0u', '1v0v', '1v1p', '1voq', '1vor', '1vos', '1vou', '1vov', '1vow', '1vox', '1voy', '1voz', '1vp0', '1vz6', '1vz7', '1vz8', '1vzq', '1w1b', '1w48', '1w5o', '1w5q', '1w6s', '1w8x', '1w91', '1wao', '1wcb', '1wcm', '1wiq', '1x9f', '1xpp', '1ydd', '1yuh', '1zgl', '1zi1', '1zoq', '1zr0', '2apt', '2aq3', '2b66', '2b9n', '2b9p', '2bfe', '2bjg', '2bsk', '2bw1', '2bw4', '2c08', '2c4w', '2c57', '2c9v', '2cdb', '2cfa', '2chv', '2ckz', '2cmp', '2co7', '2cua', '2fhe', '2fpq', '2g3p', '2h6a', '2ixo', '2ixt', '2iyv', '2iz7', '2j1o', '2j1p', '2j6e', '2j7o', '2j8u', '2j98', '2jac', '2jcg', '2jd5', '2jie', '2lyn', '2p1y', '2pf8', '2pk9', '2r92', '2r93', '2rf2', '2sbt', '2uuv', '2uwa', '2uwb', '2uwe', '2v0n', '2v2d', '2v3t', '2v47', '2v49', '2v5k', '2v7z', '2vwi', '3b6w', '3ckz', '3cl0', '3ovo', '3pga', '3ypi', '4ca2', '4cac', '4crx', '4ovo', '4tsv', '5ca2', '5cac', '5csc', '6ca2', '7ca2', '8ca2', '9ca2', '9gpb', '9pai']

##    fd = open('SEQRES_not_MODRES.txt','r')
##    lines = fd.readlines()
##    fd.close()
##    pdbs = []
##    for line in lines:
##        pdbs += [line.split()[0]]
##    pdbs = list(set(pdbs))
##    pdbs.sort()
##    print pdbs
##    os.remove('SEQRES_not_MODRES.txt')

    d_phipsi = {
        'A':[],'C':[],'D':[],'E':[],'F':[],'G':[],'H':[],'I':[],'K':[],'L':[],
        'M':[],'N':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'V':[],'W':[],'Y':[],
        }

    d_stdres = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
        'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
        'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
        'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
        'MSE':'M',
        }

    subdirs = os.listdir(path_pdb)
    subdirs.sort()
    for subdir in subdirs:
        if len(sys.argv) > 1 and subdir < sys.argv[-1][1:3]:
            continue

        print subdir
        files = os.listdir(path_pdb+subdir)
        files.sort()
        for file in files:
            if file[-3:] == '.gz':
                continue
            pdb = file[3:7]
            if len(sys.argv) > 1 and subdir == sys.argv[-1][1:3] and pdb < sys.argv[-1]:
                continue

##            if pdb not in pdbs:
##                continue

            if '-skip' in sys.argv:
                if pdb != sys.argv[sys.argv.index('-skip')+1]:
                    continue

            if pdb in [
                ## std_res_name in neither REMARK465 nor SEQRES (large peptide length)
                '1a3q','1c3w','1c3x','1ad5',
                ## hetID in SEQRES and HETATM    but not MODRES
                '2b2u','2a2x','2b7f','2c2k','2c2m','2c2o','2c2z','2aal','1bdu',
                '2ag3','2age','2agg','1an5','2ci1',
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

##                ## glutamate gamma bond (exclude these cases from phi/psi txt files - gamma CONECT?)
##                '2bh8',

##                ## chromophore CRO
##                '2h5p','2h5o','2h5r',
##                ## chromophore CH6 (MET,TYR,GLY)
##                '2h5q',
##                ## chromophore CR7 (LYS,TYR,GLY)
##                '2a46','2a47','2a48',
                ]:
                continue

            fd = open('%s%s/%s' %(path_pdb,subdir,file),'r')
            lines = fd.readlines()
            fd.close()

            ##
            ## parse header
            ##
            d_header = parse_header(lines)

            ## skip if not x-ray structure
            if d_header['EXPDTA'] != 'X-RAY':
                continue

            ## skip if no chains
            if not 'SEQRES' in d_header.keys():
                continue

            ## check if SEQRES hetIDs in MODRES if not standard residues
            skip = SEQRES_vs_MODRES(pdb,d_header)
            if skip == True:
                continue

            if 'REMARK2' in d_header.keys():
                if d_header['REMARK2'] > 2.:
                    continue
            else:
                print pdb
                stop

            ##
            ## parse coordinates
            ##
            d_coordinates, d_ATOMseq = parse_coordinates(lines, d_header, pdb,)

            for chain in d_header['SEQRES']['chains']:

                ## skip if not peptide chain
                type = d_header['SEQRES']['chains'][chain]['type']
                if type != 'peptide':
                    continue

                ## skip if all residues are unknown
                if len(d_header['SEQRES']['chains'][chain]['seq3'])*['UNK'] == d_header['SEQRES']['chains'][chain]['seq3']:
                    continue

                ## append remaining REMARK465 residues
                if 'REMARK465' in d_header.keys():
                    if chain in d_header['REMARK465']['chains'].keys():
                        if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                            l_REMARK465_res_nos = d_header['REMARK465']['chains'][chain]['residues'].keys() ## e.g. 1cd0
                            l_REMARK465_res_nos.sort()
                            d_ATOMseq,d_header = append_sequence(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,False,)

                if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3']:
                    print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3']
                    print 'ATOM  ', d_ATOMseq[chain]['seq']
                    print chain
                    print 'SEQRES', len(d_header['SEQRES']['chains'][chain]['seq3'])
                    print 'ATOM'  , len(d_ATOMseq[chain]['seq'])
                    print pdb,chain
                    stop_different_length_SEQRES_ATOM

                for i in range(1,len(d_ATOMseq[chain]['seq'])-1):

                    record = d_ATOMseq[chain]['records'][i]
                    res_no = d_ATOMseq[chain]['res_nos'][i]
                    iCode = d_ATOMseq[chain]['iCodes'][i]
                    altloc = d_ATOMseq[chain]['altlocs'][i]
                    if record == 'REMARK465':
                        continue
                    res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']

                    if res_name not in d_stdres.keys():
                        continue

                    if 'REMARK465' == d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record']:
                        continue

                    if 'N' not in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        try:
                            if d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['N']['REMARK'] != True:
                                stop
                            continue
                        except:
                            if 'REMARK470' not in d_header.keys() and d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys() == ['CA']:
                                fd = open('alpha_only_no_REMARK470_records.txt','a')
                                fd.write('%4s %1s %4i\n' %(pdb,chain,res_no))
                                fd.close()
                                break
                            else:
                                stop_alpha_only
                            
                    if 'CA' not in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        if d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA']['REMARK'] != True:
                            stop
                        continue
                    if 'C' not in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        if d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['C']['REMARK'] != True:
                            stop
                        continue
                    
                    N = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['N']['coordinate']
                    CA = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA']['coordinate']
                    C = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['C']['coordinate']

                    ## phi
                    record_prev = d_ATOMseq[chain]['records'][i-1]
                    res_no_prev = d_ATOMseq[chain]['res_nos'][i-1]
                    iCode_prev = d_ATOMseq[chain]['iCodes'][i-1]
                    altloc_prev = d_ATOMseq[chain]['altlocs'][i-1]
                    res_name_prev = d_ATOMseq[chain]['seq'][i-1]
                    if record_prev == 'REMARK465':
                        phi = '   N/A'
                    elif res_name_prev not in ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','MET','ASP','GLU','ASN','GLN','HIS','LYS','ARG','PRO','PHE','TRP','TYR',]:
                        phi = '   N/A'
                    elif 'C' not in d_coordinates['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms'].keys():
                        if d_header['REMARK470']['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['C']['REMARK'] != True:
                            stop
                        phi = '   N/A'
                    else:
                        C_prev = d_coordinates['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['C']['coordinate']
                        phi = '%6.1f' %(dihedral(C_prev,N,CA,C,))

                    ## psi
                    record_next = d_ATOMseq[chain]['records'][i+1]
                    res_no_next = d_ATOMseq[chain]['res_nos'][i+1]
                    iCode_next = d_ATOMseq[chain]['iCodes'][i+1]
                    altloc_next = d_ATOMseq[chain]['altlocs'][i+1]
                    res_name_next = d_ATOMseq[chain]['seq'][i+1]
                    if record_next == 'REMARK465':
                        psi = '   N/A'
                    elif res_name_next not in ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','MET','ASP','GLU','ASN','GLN','HIS','LYS','ARG','PRO','PHE','TRP','TYR',]:
                        psi = '   N/A'
                    elif 'N' not in d_coordinates['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms'].keys():
                        if d_header['REMARK470']['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['N']['REMARK'] != True:
                            stop
                        psi = '   N/A'
                    else:
                        N_next = d_coordinates['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['N']['coordinate']
                        psi = '%6.1f' %(dihedral(N,CA,C,N_next))


                    if phi != '   N/A':
                        dist_phi = math.sqrt(sum((C_prev-N)**2))
                    if psi != '   N/A':
                        dist_psi = math.sqrt(sum((C-N_next)**2))

                    ## phi
                    if phi != '   N/A' and dist_phi > 1.99:
                        print 'phi'
                        if len(set(['N','CA'])&set(d_coordinates['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms'].keys())) == 2:
                            N_prev = d_coordinates['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['N']['coordinate']
                            CA_prev = d_coordinates['chains'][chain]['residues'][res_no_prev]['d_iCodes'][iCode_prev]['atoms']['CA']['coordinate']
                            psi_prev = '%6.1f' %(dihedral(N_prev,CA_prev,C_prev,N))
                        else:
                            psi_prev = '   N/A'
                        gap = [res_no_prev,res_no,res_no_next] != range(res_no_prev,res_no_next+1) and [iCode_prev,iCode,iCode_next] == [' ',' ',' ',]
                        if gap == True:
                            if [res_no_prev,res_no] == range(res_no_prev,res_no+1):
                                stop
                            fd = open('phi_mit_gap.txt','a')
                            fd.write(
                                '%4s %1s %4i %4i %1s %1s %1s %1s %6s %6s %4.1f\n' %(
                                    pdb,chain,res_no_prev,res_no,iCode_prev,iCode,altloc_prev,altloc,psi_prev,phi,dist_phi
                                    )
                                )
                            fd.close()
                        else:
                            fd = open('phi_ohne_gap.txt','a')
                            fd.write(
                                '%4s %1s %4i %4i %1s %1s %1s %1s %6s %6s %4.1f\n' %(
                                    pdb,chain,res_no_prev,res_no,iCode_prev,iCode,altloc_prev,altloc,psi_prev,phi,dist_phi
                                    )
                                )
                            fd.close()
                    ## psi
                    if psi != '   N/A' and dist_psi > 1.99:
                        print 'psi'
                        if len(set(['CA','C'])&set(d_coordinates['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms'].keys())) == 2:
                            CA_next = d_coordinates['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['CA']['coordinate']
                            C_next = d_coordinates['chains'][chain]['residues'][res_no_next]['d_iCodes'][iCode_next]['atoms']['C']['coordinate']
                            phi_next = '%6.1f' %(dihedral(C,N_next,CA_next,C_next,))
                        else:
                            phi_next = '   N/A'
                        gap = [res_no_prev,res_no,res_no_next] != range(res_no_prev,res_no_next+1) and [iCode_prev,iCode,iCode_next] == [' ',' ',' ',]
                        if gap == True:
                            if [res_no,res_no_next] == range(res_no,res_no_next+1):
                                stop
                            fd = open('psi_mit_gap.txt','a')
                            fd.write(
                                '%4s %1s %4i %4i %1s %1s %1s %1s %6s %6s %4.1f\n' %(
                                    pdb,chain,res_no,res_no_next,iCode,iCode_next,altloc,altloc_next,psi,phi_next,dist_psi
                                    )
                                )
                            fd.close()
                        else:
                            fd = open('psi_ohne_gap.txt','a')
                            fd.write(
                                '%4s %1s %4i %4i %1s %1s %1s %1s %6s %6s %4.1f\n' %(
                                    pdb,chain,res_no,res_no_next,iCode,iCode_next,altloc,altloc_next,psi,phi_next,dist_psi
                                    )
                                )
                            fd.close()

##                    print '%4s %s %4i %s %3s %6s %6s' %(pdb,chain,res_no,iCode,res_name,phi,psi)

                    if psi == 0. or phi == 0.:
                        fd = open('zerovalues.txt','a')
                        fd.write('%4s %4i %1s %1s %6s %6s' %(pdb, res_no, iCode, altloc, phi, psi))
                    d_phipsi[d_stdres[res_name]] += ['%6s %6s\n' %(phi,psi,)]
##                    if res_name_prev in d_stdres.keys():
##                        d_phipsi[d_stdres[res_name_prev]] += ['%6s %6s\n' %(phi,psi,)]
##                    if res_name_next in d_stdres.keys():
##                        d_phipsi[d_stdres[res_name_next]] += ['%6s %6s\n' %(phi,psi,)]

    for res in d_phipsi.keys():
        fd = open('phipsi_all_%s.txt' %(res),'w')
        fd.writelines(d_phipsi[res])
        fd.close()

    return

def SEQRES_vs_MODRES(pdb,d_header):

    skip = False

    for chain in d_header['SEQRES']['chains'].keys():
        l_hetIDs_MODRES = []
        if 'MODRES' in d_header.keys():
            if chain in d_header['MODRES'].keys():
                for res_no in d_header['MODRES'][chain].keys():
                    for iCode in d_header['MODRES'][chain][res_no].keys():
                        l_hetIDs_MODRES += [d_header['MODRES'][chain][res_no][iCode]]
        l_hetIDs = d_header['SEQRES']['chains'][chain]['seq3']
        ## e.g. 1c58
        if l_hetIDs == len(l_hetIDs)*['GLC']:
##            skip = True
            continue
        for hetID in set(l_hetIDs):
            if hetID not in d_pep.keys()+l_nuc:
                if hetID == 'GLC':
                    stop
                if not hetID in l_hetIDs_MODRES:
                    fd = open('SEQRES_not_MODRES.txt','a')
                    fd.write('%s %s %s\n' %(pdb,chain,hetID))
                    fd.close()
                    skip = True

    return skip


def check_if_SEQRESres(res_name,record,d_header,chain,res_no,iCode):

    if res_name not in d_pep.keys()+l_nuc and record == 'ATOM':
        print res_name,record
        stop

    if res_name in d_pep.keys()+l_nuc and record in ['ATOM','REMARK465',]:
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


def dihedral(c1,c2,c3,c4):

    import numpy,math

    v1 = c2-c1
    v2 = c3-c2
    v3 = c4-c3

    angle = math.atan2(
        numpy.dot(
            math.sqrt(sum(v2*v2))*v1,
            cross(v2,v3),
            ),
        numpy.dot(
            cross(v1,v2),
            cross(v2,v3),
            ),
        )
    angle *= 180./math.pi

    return angle


def cross(v1,v2):

    import numpy

    n = numpy.array([
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0],
        ])

    return n


def parse_coordinates(lines, d_header, pdb,):

    d_atomnos = {}
    d_coordinates = {}
    d_ATOMseq = {}

    for i in range(len(lines)):
        line = lines[i]
        record = line[:6].strip()

        if record == 'ATOM':
            d_coordinates, d_line, d_ATOMseq = parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM', d_ATOMseq, pdb)
            d_atomnos[d_line['atom_no']] = d_line

        elif record == 'HETATM':
            res_name = line[17:20].strip()
            chain = line[21]
            ## water
            if res_name in ['D2O','H2O',]:
                print pdb, res_name
                stop
            elif res_name in ['HOH','DOD',]: ## DOD e.g. 2d4i
                continue
            ## modified residue of polypeptide or polynucleotide
            elif res_name == 'MSE' and 'MSE' in d_header['SEQRES']['chains'][chain]['seq3']: ## e.g. 1gcj,2e1a
                d_coordinates, d_line, d_ATOMseq = parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM', d_ATOMseq, pdb)
            ## (poly)saccharide or other hetero compound
            else:
                d_coordinates, d_line, d_ATOMseq = parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM', d_ATOMseq, pdb)
            atom_no = d_line['atom_no']
            d_atomnos[atom_no] = d_line

        elif record == 'MODEL':
            model = int(line.split()[1])

        elif record == 'ENDMDL':
            break

    return d_coordinates, d_ATOMseq


def parse_header(lines):

    d_header = {}

    for i in range(len(lines)):
        line = lines[i]

        record = line[:6].strip()

        if record == 'EXPDTA': ## section 2
            methods = line[10:].strip().split(',')[0]
            if methods[:3] == 'NMR':
                methods = 'NMR'
            elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                methods = 'X-RAY'
            d_header['EXPDTA'] = methods

        elif record == 'REMARK': ## section 2
            d_header = parse_recordREMARK(d_header, line, i, lines)

        elif record == 'SEQRES': ## section 3
            d_header = parse_recordSEQRES(line, d_header)

        elif record == 'MODRES': ## section 3
            parse_recordMODRES(line, d_header,)

        elif record == 'HET': ## section 4
            parse_recordHET(line, d_header,)

    return d_header


def parse_recordATOM(line, d_coordinates, lines, i, d_header, record, d_ATOMseq, pdb):

    import Numeric

    atom_no = int(line[6:11])
    atom_name = line[12:16].strip()
    altloc = line[16]
    res_name_ATOM = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    element = line[76:78].strip()
    coordinate = Numeric.array([x, y, z])

    skip = False

    ## modified residue?
    SEQRESres = check_if_SEQRESres(res_name_ATOM,record,d_header,chain,res_no,iCode,)
    if SEQRESres == False:
        skip = True

    if chain in d_header['SEQRES']['chains']:
        type = d_header['SEQRES']['chains'][chain]['type']
        if type != 'peptide':
            skip = True
    else:
        skip = True

    if skip == True:
        return d_coordinates, {
            'chain':chain,
            'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
            'atom_no':atom_no,'atom_name':atom_name,'element':element,
            }, d_ATOMseq

    if not 'chains' in d_coordinates.keys():
        d_coordinates['chains'] = {}
    if not chain in d_coordinates['chains'].keys():
        d_coordinates['chains'][chain] = {}
    if not 'residues' in d_coordinates['chains'][chain].keys():
        d_coordinates['chains'][chain]['residues'] = {}

    ## res_no
    if not res_no in d_coordinates['chains'][chain]['residues'].keys():
        d_coordinates['chains'][chain]['residues'][res_no] = {}

    ## res_no > d_iCodes
    if not 'd_iCodes' in d_coordinates['chains'][chain]['residues'][res_no].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
    ## d_iCodes > iCode
    if not iCode in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes']:
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}


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
                                SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(l_REMARK465_res_names)]
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
            res_name_SEQRES = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
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
            d_ATOMseq,d_header = append_sequence(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,None,False,True,)
        ## REMARK465 before ATOM (certain)
        elif REMARK465_before_ATOM == True:# and res_name_ATOM != res_name_SEQRES:# and res_name_REMARK465 == res_name_SEQRES:
            d_ATOMseq,d_header = append_sequence(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,True,)
        else:
            print '---'
            SEQRES_res_name = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
            SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
            print chain, res_no, iCode
            print line
            print d_ATOMseq[chain]['seq']
            print d_header['SEQRES']['chains'][chain]['seq3']
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

        if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]:
            print '*******'
            print 'ATOM  ', d_ATOMseq[chain]['seq']
            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
            print res_name_ATOM
            print line
            print res_name_REMARK465
            print pdb, chain, res_no, iCode
            stop_sequence_difference

    ## iCode > record
    d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = record

    ## iCode > res_name
    if not 'res_name' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name_ATOM
    ## check that res_name is correct (e.g. 2fes:L:1)
    elif d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name_ATOM and altloc == ' ': ## 1fh2:A:30 altloc
        print d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
        print res_name_ATOM
        print altloc
        print chain, res_no, iCode
        print line
        stop_add_altloc

    ## res_no > l_iCodes
    if not 'l_iCodes' in d_coordinates['chains'][chain]['residues'][res_no].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] = []
    ## l_iCodes > iCode
    if not iCode in d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes']:
        d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

    ## iCode > atoms
    if not 'atoms' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
    ## atoms > atom_name > coordinate
    if not atom_name in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}

    ## iCode > record
    if not 'record' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = record

  
    return d_coordinates, {
        'chain':chain,
        'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
        'atom_no':atom_no,'atom_name':atom_name,'element':element,
        }, d_ATOMseq


def append_sequence(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,append_REMARK465,append_ATOM):

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


def parse_recordHET(line, d_header,):

    hetID = line[7:10].strip()
    ## continue if water
    if hetID in ['HOH','DOD',]:
        return d_header
    if hetID in ['D2O','DOD','H2O',]:
        print hetID
        stop
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
        d_header['HET'][chain][res_no][iCode] = [hetID]
    elif hetID not in d_header['HET'][chain][res_no][iCode]:
        d_header['HET'][chain][res_no][iCode] += [hetID]

    return d_header

def parse_recordMODRES(line, d_header,):

    hetID = line[12:15].strip()
    chain = line[16]
    res_no = int(line[18:22])
    iCode = line[22]
    res_name = line[24:27].strip()
    txt = line[29:80].strip()

    ## return if e.g. glycosylation site (e.g. 3c43)
    if hetID in set(d_pep.keys())-set(['MSE']) and res_name in set(d_pep.keys())-set(['MSE']):
        return d_header

    if 'MODRES' not in d_header.keys():
        d_header['MODRES'] = {}
    if chain not in d_header['MODRES'].keys():
        d_header['MODRES'][chain] = {}
    if res_no not in d_header['MODRES'][chain].keys():
        d_header['MODRES'][chain][res_no] = {}
    if iCode not in d_header['MODRES'][chain][res_no].keys():
        d_header['MODRES'][chain][res_no][iCode] = [hetID]
    elif hetID not in d_header['MODRES'][chain][res_no][iCode]:
        d_header['MODRES'][chain][res_no][iCode] += [hetID]

    return d_header


def parse_recordREMARK(d_header, line, i, lines):

    remark = int(line[6:10])

    if remark == 2:
        if (
            line[11:38] == 'RESOLUTION. NOT APPLICABLE.'
            or
            line[11:38] == 'RESOLUTION. NULL ANGSTROMS.'
            ):
            d_header['REMARK2'] = 'N/A'
            return d_header
        try:
            resolution = float(line[22:27])
            d_header['REMARK2'] = resolution
        except:
            pass

    elif remark == 465: ## missing residues

        d_header = parse_recordREMARK465(line, lines, i, d_header)

    elif remark == 470: ## missing atoms

        d_header = parse_recordREMARK470(line, lines, i, d_header)

    return d_header


def parse_recordREMARK465(line, lines, i, d_header):

    ## missing residues

    if line[10:].strip() in [
        'M RES C SSSEQI',
        'M RES C  SSEQI',
        ]:

        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK 465':
                break

            try:
                model = int(lines[j][12:14])
            except:
                model = 1
            res_name = lines[j][15:18].strip()
            chain = lines[j][19]
            res_no = int(lines[j][22:26])
            iCode = lines[j][26]

            if model != 1: ## e.g. 1ohh
                return d_header

            if not 'REMARK465' in d_header.keys():
                d_header['REMARK465'] = {}
            if not 'chains' in d_header['REMARK465'].keys():
                d_header['REMARK465']['chains'] = {}
            if not chain in d_header['REMARK465']['chains'].keys():
                d_header['REMARK465']['chains'][chain] = {}
            if not 'residues' in d_header['REMARK465']['chains'][chain].keys():
                d_header['REMARK465']['chains'][chain]['residues'] = {}
            if not res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                d_header['REMARK465']['chains'][chain]['residues'][res_no] = {}

            ## res_no > d_iCodes
            if not 'd_iCodes' in d_header['REMARK465']['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
            ## d_iCodes > iCode
            if not iCode in d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes']:
                d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

            ## iCode > record
            d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = 'REMARK465'

            ## iCode > res_name
            d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name

            ## res_no > l_iCodes
            if not 'l_iCodes' in d_header['REMARK465']['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes'] = []
            ## l_iCodes > iCode
            if not iCode in d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']:
                d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

            ## iCode > REMARK
            if not 'REMARK' in d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] = True

            ## iCode > record
            d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = 'REMARK465'


    return d_header


def parse_recordREMARK470(line, lines, i, d_header):

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
            if not 'REMARK470' in d_header.keys():
                d_header['REMARK470'] = {}
            if not 'chains' in d_header['REMARK470'].keys():
                d_header['REMARK470']['chains'] = {}
            if not chain in d_header['REMARK470']['chains'].keys():
                d_header['REMARK470']['chains'][chain] = {}
            if not 'residues' in d_header['REMARK470']['chains'][chain].keys():
                d_header['REMARK470']['chains'][chain]['residues'] = {}
            if not res_no in d_header['REMARK470']['chains'][chain]['residues'].keys():
                d_header['REMARK470']['chains'][chain]['residues'][res_no] = {}

            ## res_no > l_iCodes
            if not 'l_iCodes' in d_header['REMARK470']['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] = []
            ## l_iCodes > iCode
            if not iCode in d_header['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes']:
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

            ## res_no > d_iCodes
            if not 'd_iCodes' in d_header['REMARK470']['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
            ## d_iCodes > iCode
            if not iCode in d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes']:
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

            ## iCode > atoms
            if not 'atoms' in d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
            ## atoms > atom_name > coordinate
            for atom_name in atoms:
                if not atom_name in d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                    d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                ## iCode > REMARK
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['REMARK'] = True

            ## iCode > res_name
            if not 'res_name' in d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name


    elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
        print pdb1
        print pdb2
        notexpected

    return d_header


def parse_recordSEQRES(line, d_header):

    chain = line[11]

    if 'SEQRES' not in d_header.keys():
        d_header['SEQRES'] = {}
    if 'chains' not in d_header['SEQRES'].keys():
        d_header['SEQRES']['chains'] = {}
    if chain not in d_header['SEQRES']['chains'].keys():
        d_header['SEQRES']['chains'][chain] = {}
    if not 'type' in d_header['SEQRES']['chains'][chain].keys():
        d_header['SEQRES']['chains'][chain]['type'] = 'N/A'

    l_residues = line[19:70].split()

    s_residues = ''
    for i in range(len(l_residues)):
        residue = l_residues[i]
        if residue in d_pep.keys():
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
            s_residues += d_pep[residue]
        elif residue in l_nuc:
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'nucleotide'
            s_residues += residue
        elif residue in ['GLC','GAL','MAN','FRU']:
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'saccharide'
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


if __name__ == '__main__':
    main()

def main():

    import os
    d_bfactors1 = init_bfactor_dic()

    path = '/oxygenase_local/data/pdb/'
    subdirs = os.listdir(path)
    for i in range(len(subdirs)):
        subdir = subdirs[i]
##        if subdir != 'jj':
##            continue
        print '%s/%s %s' %(i+1,len(subdirs), subdir)
        files = os.listdir(path+subdir)
        for file in files:
            l_bfactors = []
            d_bfactors2 = init_bfactor_dic()
            d_ss = {}
            fd = open('%s%s/%s' %(path, subdir, file), 'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                record = line[:6].strip()
                if record == 'EXPDTA':
                    if line[10:13] == 'NMR':
                        break
                    if 'NMR' in line:
                        print file
                        stop
                elif record == 'HELIX':
                    if line[19] != line[31]:
                        stop
                    chain = line[19]
                    res_no1 = int(line[21:25])
                    iCode1 = line[25]
                    res_no2 = int(line[33:37])
                    iCode2 = line[37]
                    d_ss = SSrecord2dic(record,chain,res_no1,res_no2,iCode1,iCode2,d_ss)
                elif record == 'SHEET':
                    if line[21] != line[32]:
                        stop
                    chain = line[21]
                    res_no1 = int(line[22:26])
                    iCode1 = line[26]
                    res_no2 = int(line[33:37])
                    iCode2 = line[37]
                    d_ss = SSrecord2dic(record,chain,res_no1,res_no2,iCode1,iCode2,d_ss)
                elif record == 'TURN':
                    if line[19] != line[30]:
                        stop
                    chain = line[19]
                    res_no1 = int(line[20:24])
                    iCode1 = line[24]
                    res_no2 = int(line[31:35])
                    iCode2 = line[35]
                    d_ss = SSrecord2dic(record,chain,res_no1,res_no2,iCode1,iCode2,d_ss)
                elif record == 'ATOM':

                    res_name = line[17:20]
                    ## amino acid residues only
                    if res_name not in d_bfactors2.keys():
                        continue

                    atom_name = line[12:16].strip()
                    ## backbone atoms only
                    if atom_name not in ['N','CA','C','O',]:
                        continue

                    chain = line[21]
                    res_no = int(line[22:26])
                    iCode = line[26]    
                    bfactor = float(line[60:66])

                    l_bfactors += [bfactor]
                    ss = ATOMrecord2SS(chain,res_no,iCode,d_ss,)
                    d_bfactors2[res_name][ss][atom_name] += [bfactor]
            ## no amino acids
            if len(l_bfactors) == 0:
                continue
            average = sum(l_bfactors)/len(l_bfactors)
            ## NMR structure
            if average == 0:
                continue
            ## identical temperature factors for all atoms
            if len(l_bfactors)*[average] == l_bfactors:
                print file
                continue
            if average in [2,3,4,5,6,7,8,9,99,90,50,20,25,1,100,10,0]:
                print average
                print file
                stop

           
            for res_name in d_bfactors2.keys():
                for ss in d_bfactors2[res_name].keys():
                    for atom_name in d_bfactors2[res_name][ss].keys():
                        for bfactor in d_bfactors2[res_name][ss][atom_name]:
                            d_bfactors1[res_name][ss][atom_name] += [bfactor/average]

    for res_name in d_bfactors1.keys():
        for ss in d_bfactors1[res_name].keys():
            for atom_name in d_bfactors1[res_name][ss].keys():
                if len(d_bfactors1[res_name][ss][atom_name]) == 0:
                    print res_name, ss, atom_name, 0, 'N/A'
                else:
                    print res_name, ss, atom_name, len(d_bfactors1[res_name][ss][atom_name]), sum(d_bfactors1[res_name][ss][atom_name])/len(d_bfactors1[res_name][ss][atom_name])

    ## n, average, stderr
    ##GLU 1.10011273863
    ##LYS 1.0942662686
    ##ASP 1.07577365273
    ##PRO 1.04444999734
    ##ASN 1.04158659424
    ##GLN 1.04053967041
    ##SER 1.03339041002
    ##GLY 1.02060932502
    ##ARG 1.01252616164
    ##THR 0.992003037451
    ##HIS 0.978963717117
    ##ALA 0.970650226677
    ##MET 0.967888835552
    ##LEU 0.946739562158
    ##CYS 0.943281173243
    ##VAL 0.928932470807
    ##PHE 0.92506499047
    ##TYR 0.924553016351
    ##ILE 0.924389814215
    ##TRP 0.913661994108

def init_bfactor_dic():

    d_bfactors = {
        'ALA':{},'CYS':{},'ASP':{},'GLU':{},'PHE':{},
        'GLY':{},'HIS':{},'ILE':{},'LYS':{},'LEU':{},
        'MET':{},'ASN':{},'PRO':{},'GLN':{},'ARG':{},
        'SER':{},'THR':{},'VAL':{},'TRP':{},'TYR':{},
        }
    for res_name in d_bfactors.keys():
        for ss in ['HELIX','SHEET','TURN','other',]:
            d_bfactors[res_name][ss] = {}
            for atom_name in ['N','CA','C','O',]:
                d_bfactors[res_name][ss][atom_name] = []

    return d_bfactors

def ATOMrecord2SS(chain,res_no,iCode,d_ss,):
    if chain not in d_ss.keys():
        return 'other'
    if res_no not in d_ss[chain].keys():
        return 'other'
    if 'min' in d_ss[chain][res_no]:
        if d_ss[chain][res_no]['min'] > iCode:
            return 'other'
    if 'max' in d_ss[chain][res_no]:
        if d_ss[chain][res_no]['max'] < iCode:
            return 'other'
    return d_ss[chain][res_no]['ss']

def SSrecord2dic(ss,chain,res_no1,res_no2,iCode1,iCode2,d_ss,):
    if chain not in d_ss.keys():
        d_ss[chain] = {}
    for res_no in range(res_no1,res_no2+1):
        d_ss[chain][res_no] = {'ss':ss}
        if res_no == res_no2:
            d_ss[chain][res_no]['min'] = iCode1
        if res_no == res_no2:
            d_ss[chain][res_no]['max'] = iCode2
    return d_ss
    
if __name__=='__main__':
    main()

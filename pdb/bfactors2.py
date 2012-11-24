## calculate average normalized bfactors of all pdbs

import statistics

def main():

    import os
    d_bfactors1 = init_bfactor_dic()

    path = '/oxygenase_local/data/pdb'
    path = '/media/WDMyBook1TB/2TB/pdb'
    l_dns = os.listdir(path)
    l_dns.sort()
    for i in range(len(l_dns)):
        continue
        dn = l_dns[i]
##        if dn != 'e7':
##            continue
##        if dn < 'xk':
##            continue
        print '%s/%s %s' %(i+1,len(l_dns), dn)
        l_fns = os.listdir('%s/%s' %(path,dn))
        for fn in l_fns:
            continue
            pdb = fn[3:7]
            if pdb in ['2wto','2wtp',]: ## remediation
                continue
            l_bfactors = []
            d_bfactors2 = init_bfactor_dic()
            d_ss = {}
            fd = open('%s/%s/%s' %(path, dn, fn), 'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                record = line[:6].strip()
                if record == 'EXPDTA':
                    if line[10:13] == 'NMR':
                        break
                    elif line[10:22] == 'SOLUTION NMR':
                        break
                    elif line[10:25] == 'SOLID-STATE NMR':
                        break
                    elif line[10:41] == 'X-RAY DIFFRACTION; SOLUTION NMR': ## 1iob
                        break
                    elif line[10:43] == 'SOLUTION SCATTERING; SOLUTION NMR': ## 2klj
                        break
                    elif line[10:29] == 'ELECTRON MICROSCOPY':
                        break
                    elif line[10:27] != 'X-RAY DIFFRACTION':
                        print line
                        break
                    if 'NMR' in line:
                        print fn
                        print line
                        stop
                elif record == 'HEADER':
                    year = int(line[50:59][-2:])
                elif record == 'HELIX':
                    if line[19] != line[31]:
                        print line
                        print pdb
                        stop
                    chain = line[19]
                    res_no1 = int(line[21:25])
                    iCode1 = line[25]
                    res_no2 = int(line[33:37])
                    iCode2 = line[37]
                    d_ss = SSrecord2dic(record,chain,res_no1,res_no2,iCode1,iCode2,d_ss)
                elif record == 'SHEET':
                    if line[21] != line[32]:
                        print line
                        print pdb
                        stop
                    chain = line[21]
                    res_no1 = int(line[22:26])
                    iCode1 = line[26]
                    res_no2 = int(line[33:37])
                    iCode2 = line[37]
                    d_ss = SSrecord2dic(record,chain,res_no1,res_no2,iCode1,iCode2,d_ss)
                elif record == 'TURN':
                    stop
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

            ## no amino acids (or break)
            if len(l_bfactors) == 0:
                continue
            average = sum(l_bfactors)/len(l_bfactors)
            ## NMR structure
            if average == 0:
                continue
            ## identical temperature factors for all atoms
            if len(l_bfactors)*[average] == l_bfactors:
                print 'all bfactors same', pdb, 'year', year, 'bfac', average
                continue
            if average in [2,3,4,5,6,7,8,9,99,90,50,20,25,1,100,10,0]:
                print average
                print fn
                stop

            for res_name in d_bfactors2.keys():
                for ss in d_bfactors2[res_name].keys():
                    for atom_name in d_bfactors2[res_name][ss].keys():
                        l_bfactors = d_bfactors2[res_name][ss][atom_name]
                        if len(l_bfactors) == 0:
                            continue
                        lines_out = ['%s\n' %(bfac/average) for bfac in l_bfactors]
                        fd = open('bfac_%s_%s_%s.txt' %(res_name,ss,atom_name),'a')
                        fd.writelines(lines_out)
                        fd.close()
##                        for bfactor in d_bfactors2[res_name][ss][atom_name]:
##                            bfactor_normalized = bfactor/average
##                            d_bfactors1[res_name][ss][atom_name] += [bfactor_normalized]

    for res_name in d_bfactors1.keys():
        l_bfactors_res = []
        for ss in ['HELIX','SHEET','OTHER',]:
            for atom_name in d_bfactors1[res_name][ss].keys():
                fd = open('bfac_%s_%s_%s.txt' %(res_name,ss,atom_name),'r')
                lines = fd.readlines()
                fd.close()
                l_bfactors = [float(line) for line in lines]
                l_bfactors_res += l_bfactors
##                average, stderr = statistics.do_stderr(l_bfactors)
##                print '%s\t%s\t%s\t%s\t%s\t%s' %(
##                    res_name, ss, atom_name,
##                    len(l_bfactors),
##                    sum(l_bfactors)/len(l_bfactors),
##                    stderr,
##                    )
        average, stderr = statistics.do_stderr(l_bfactors_res)
        print '%s\t%s\t%s\t%s' %(
            res_name,
            len(l_bfactors_res),
            sum(l_bfactors_res)/len(l_bfactors_res),
            stderr,
            )
##                if len(d_bfactors1[res_name][ss][atom_name]) == 0:
##                    print res_name, ss, atom_name, 0, 'N/A'
##                else:
##                    print '%s\t%s\t%s\t%s\t%s' %(
##                        res_name, ss, atom_name,
##                        len(d_bfactors1[res_name][ss][atom_name]),
##                        sum(d_bfactors1[res_name][ss][atom_name])/len(d_bfactors1[res_name][ss][atom_name]),
##                        )

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
        d_bfactors[res_name] = {}
        for ss in ['HELIX','SHEET','OTHER',]:
            d_bfactors[res_name][ss] = {}
            for atom_name in ['N','CA','C','O',]:
                d_bfactors[res_name][ss][atom_name] = []

    return d_bfactors

def ATOMrecord2SS(chain,res_no,iCode,d_ss,):
    if chain not in d_ss.keys():
        return 'OTHER'
    if res_no not in d_ss[chain].keys():
        return 'OTHER'
    if 'min' in d_ss[chain][res_no]:
        if d_ss[chain][res_no]['min'] > iCode:
            return 'OTHER'
    if 'max' in d_ss[chain][res_no]:
        if d_ss[chain][res_no]['max'] < iCode:
            return 'OTHER'
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

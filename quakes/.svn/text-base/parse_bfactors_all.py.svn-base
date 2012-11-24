import os
d_bfactors1 = {
    'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],
    'GLY':[],'HIS':[],'ILE':[],'LYS':[],'LEU':[],
    'MET':[],'ASN':[],'PRO':[],'GLN':[],'ARG':[],
    'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[],
    }
path = '/oxygenase_local/data/pdb/'
subdirs = os.listdir(path)
for i in range(len(subdirs)):
    subdir = subdirs[i]
    print '%s/%s %s' %(i+1,len(subdirs), subdir)
    files = os.listdir(path+subdir)
    for file in files:
        l_bfactors = []
        d_bfactors2 = {
            'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],
            'GLY':[],'HIS':[],'ILE':[],'LYS':[],'LEU':[],
            'MET':[],'ASN':[],'PRO':[],'GLN':[],'ARG':[],
            'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[],
            }
        fd = open('%s%s/%s' %(path, subdir, file), 'r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            record = line[:6].strip()
            if record == 'ATOM':
                res_name = line[17:20]
                ## amino acid residues only
                if res_name not in d_bfactors2.keys():
                    continue
                atom_name = line[12:16].strip()
                ## backbone atoms only
                if atom_name not in ['N','CA','C','O',]:
                    continue
                bfactor = float(line[60:66])
                l_bfactors += [bfactor]
                d_bfactors2[res_name] += [bfactor]
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
            stop
            continue
        if average in [99,90,50,25,1,100,10,0]:
            print average
            stop
        for res_name in d_bfactors2.keys():
            for bfactor in d_bfactors2[res_name]:
                d_bfactors1[res_name] += [bfactor/average]
for res_name in d_bfactors1.keys():
    print res_name, sum(d_bfactors1[res_name])/len(d_bfactors1[res_name])
stop

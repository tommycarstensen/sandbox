import os

pdbpath = '/oxygenase_local/data/pdb/'

l_pdbs = []

l_modres = []

subdirs = os.listdir(pdbpath)
subdirs.sort()
for subdir in subdirs:
    if subdir < 'ag':
        continue
    print subdir,
    files = os.listdir(pdbpath+subdir)
    for file in files:
        if file[3:7] in [
            '2bl7','2bxt','2c63','1w9n','1wco',
            '1a0h','1af0',]:
            continue
        l_pdbs += ['%s' %(file[3:7])]
        fd = open('%s%s/%s' %(pdbpath,subdir,file),'r')
        lines = fd.readlines()
        for line in lines:
            if line[:6] == 'MODRES':
                res_name = line[24:27].strip()
                if res_name in ['C','A','G','U','T','I','DA','DG']:
                    continue
                if res_name in ['MSE','TYS','TPO','SEP',]:
                    continue
                if res_name in [
                    '','TRP','GLN','PRO','ARG','LYS','GLY','ALA','VAL','LEU','ILE','SER','THR','ASP','GLU','PHE','HIS','TYR','MET','ASN','CYS'
                    ]:
                    if line[12:15].strip() == res_name and line[29:].strip() not in [
                        'GLYCOSYLATION SITE','N-GLYCOSYLATION',
                        'COVALENTLY ATTACHED TO PLP',
                        ]:
                        print
                        print file, line[:-1]
                        stop
                    res_name = line[12:15].strip()
                    l_modres += [res_name]
                else:
##                    if line[12:15].strip() in [
##                        'AMU', ## substrate
##                        'QUI','HQU', ## nucleotide
##                        'CPC', ## non-AA-derivative (cyclic peptide)
##                        'DIP','DIY','DIX', ## non-AA-derivative (non-peptide)
##                        'PSE', ## non-AA-derivative (synthetic peptide)
##                        ]:
##                        continue
                    print
                    print line[:-1]

l_pdbs = list(set(l_pdbs))

print len(l_modres)
l_modres = list(set(l_modres))
print len(l_modres)
l_modres.sort()
print l_modres

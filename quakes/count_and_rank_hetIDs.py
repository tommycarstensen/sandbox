## which hetero compounds occur most often and should therefore be dealt with in smallmolecules.py?

import os

d1 = {}
d2 = {}

path_pdb = '/data/pdb-v3.2'
l_dn = os.listdir(path_pdb)
l_dn.sort()
d_resolutions = {
    'False':[],
    'True':[],
    }
for dn in l_dn:
##    if dn < 'by':
##        continue
    print dn
    l_fn = os.listdir('%s/%s' %(path_pdb,dn,))
    l_fn.sort()
    for fn in l_fn:
        pdb = fn[3:7]
        fd = open('%s/%s/%s' %(path_pdb,dn,fn,),'r')
        lines = fd.readlines()
        fd.close()
        l_hetIDs = []
        for line in lines:
            if line[:6] == 'ATOM  ':
                break
            elif line[:6] == 'HET   ':
                hetID = line[7:10]

                if not hetID in d1.keys():
                    d1[hetID] = 0
                d1[hetID] += 1

                l_hetIDs += [hetID]

        set_hetIDs = set(l_hetIDs)
        for hetID in set_hetIDs:
            if not hetID in d2.keys():
                d2[hetID] = 0
            d2[hetID] += 1

for d in [d1,d2,]:
    l = []
    for hetID,count in d.items():
        if count > 1:
            l += [[count,hetID,]]
    l.sort()
    print l[-100:]

import os

fd = open('../ligasite.csv','r')
lines = fd.readlines()
fd.close()
d_residues = {}
d_holo = {}
for line in lines[1:]:
    l = line.split(',')
    apoID = l[0]
    residue = l[1]
    holoIDs = l[2].strip().split('-')
    if not apoID in d_residues.keys():
        d_residues[apoID] = []
    d_residues[apoID] += [residue]
    d_holo[apoID] = holoIDs
l_pdbs = d_residues.keys()
print len(l_pdbs)
stop

fd = open('../bc-50.out','r')
lines = fd.readlines()
fd.close()
d = {}
for i_line in range(len(lines)):
    line = lines[i_line]
    l_cluster = line.split()
    for i in range(len(l_cluster)):
        l_cluster[i] = l_cluster[i][:4].lower()
    del l
    l_apoIDs = list( set(l_cluster) & set(l_pdbs) )
    if len(l_apoIDs) > 1:
        for pdb in l_apoIDs:
            l_pdbs.remove(pdb)
        for apoID in l_apoIDs:
            l_holoIDs = d_holo[apoID]
            if len( set(l_holoIDs) & set(l_cluster) ) == 0:
                print l_cluster
                print 'apo', apoID, l_apoIDs
                print 'holo', l_holoIDs
                print i_line
                print '2lzt' in l_cluster
                stop_ok1
            ## all holoIDs in cluster
            elif len( set(l_holoIDs) & set(l_cluster) ) == len(l_holoIDs):
                pass
            else:
                print l_holoIDs
                print l_cluster
                stop_wrong

            
        d[i_line] = l
##    l = line.split()
##    cluster = int(l[0])
##    rank = int(l[1])
##    chainID = l[2]
print d
print l_pdbs
stop
    

set_pdbs = set()
for line in lines[1:]:
    pdb = line[:4]
    set_pdbs |= set([pdb])

l_pdbs = list(set_pdbs)
l_pdbs.sort()

for i_pdb in range(3,len(l_pdbs),4,):
    pdb = l_pdbs[i_pdb]
    fd = open('/data/pdb-v3.2/%2s/pdb%s.ent' %(pdb[1:3],pdb,), 'r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            chain = line[21]
            break
    print pdb, chain
    os.system('nicer python /home/people/tc/svn/GoodVibes/goodvibes_ligand.py --pdb %s --chain %s' %(pdb,chain,))
    

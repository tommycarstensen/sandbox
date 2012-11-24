import os, numpy, math

cwd = os.getcwd()

fd = open('../LIGSITEcsc.csv','r')
lines = fd.readlines()
fd.close()
d_pdb = {}
for line in lines:
    ligand = line[:3]
    pdb = line[4:8].lower()
    chain_ligand = line[11]
    d_pdb[pdb] = {'ligand':ligand,'chain_ligand':chain_ligand,}

count_all = 0
count_hit = 0

l_fn = os.listdir(cwd)
for fn in l_fn:
    if fn[-4:] != '.pdb':
        continue
    if os.path.getsize(fn) == 1:
        continue
    pdb = fn[:4]
    print pdb, d_pdb.keys()
    if not pdb in d_pdb.keys():
        continue
    count_all += 1

    l_coords_ligand = []
##    fd = open('/data/pdb-v3.2/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
    fd = open('/media/WDMyBook1TB/2TB/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()
    l_coords_ligand = []
    for line in lines:
        record = line[:6].strip()
        if record != 'HETATM':
            continue
        res_name = line[17:20]
        if res_name != d_pdb[pdb]['ligand']:
            continue
        chain = line[21]
        if chain != d_pdb[pdb]['chain_ligand']:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        l_coords_ligand += [coord]
    coord_ligand_center = sum(l_coords_ligand)/len(l_coords_ligand)
    
    fd = open(fn,'r')
    lines = fd.readlines()
    fd.close()
    l_coords_probe = []
    l_overlaps = []
    for line in lines:
        record = line[:6].strip()
        if record != 'HETATM':
            continue
        res_name = line[17:20]
        if res_name != 'EXT':
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        l_coords_probe += [coord]
        tempfactor = float(line[60:66])
        overlap = (100-tempfactor)/100.
        l_overlaps += [overlap]

    l = []
    for i in range(len(l_overlaps)):
        l += [[
            l_overlaps[i],
            l_coords_probe[i][0],l_coords_probe[i][1],l_coords_probe[i][2],
            ]]
    l.sort()

    for i in range(len(l)):
        if i >= 3:
            continue
        x = l[i][1]
        y = l[i][2]
        z = l[i][3]
        for j in range(len(l_coords_ligand)):
            coord_ligand = l_coords_ligand[j]
            diff = coord_ligand-numpy.array([x,y,z,])
            dist = math.sqrt(sum(diff**2))
            if dist < 4:
                count_hit += 1
                break
        if dist > 4:
            diff = coord_ligand_center-numpy.array([x,y,z,])
            print pdb, i, math.sqrt(sum(diff**2))
        else:
            break

print count_hit, count_all, count_hit/float(count_all)

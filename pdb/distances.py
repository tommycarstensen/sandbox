import math, sys, numpy

fn = sys.argv[-3]
chainID1 = sys.argv[-1]
chainID2 = sys.argv[-2]

fd = open(fn,'r')
lines = fd.readlines()
fd.close()

d_coords = {
    chainID1:{},
    chainID2:{},
    }
for i in range(len(lines)):
    line = lines[i]
    record = line[:6].strip()
    if record in ['ATOM',]:
        atom_name = line[12:16].strip()
        if atom_name in ['N','C','O',]:
            continue
        res_name = line[16:19]
        if res_name != 'GLY' and atom_name == 'CA':
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        chain = line[21]
        res_no = int(line[22:26])
        if not res_no in d_coords[chain]:
            d_coords[chain][res_no] = []
        d_coords[chain][res_no] += [coord]

d_contacts = {}
for res_no1 in d_coords[chainID1]:
    bool_break = False
    for res_no2 in d_coords[chainID1]:
        for coord1 in d_coords[chainID1][res_no1]:
            for coord2 in d_coords[chainID2][res_no2]:
                sq_dist = sum((coord1-coord2)**2)
                if sq_dist < 25:
                    print res_no1, res_no2, sq_dist
                    bool_break = True
                    if not res_no1 in d_contacts.keys():
                        d_contacts[res_no1] = []
                    d_contacts[res_no1] += [res_no2]
                    break
            if bool_break == True:
                break
##        if bool_break == True:
##            break

print d_contacts

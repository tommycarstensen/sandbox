import numpy

d_atom_nos = {}
for i in range(0,100,1):
    fd = open('/local/tc/MD_2vb1/amber99sb_CYM_35C/CHANEU/trjconv/2vb1_MD%i.pdb' %(i),'r')
    lines = fd.readlines()
    fd.close()

    ## active site residue coordinates
    l_coords = []

    ## water oxygen atoms vicinal to active site residue
##    d_atom_nos = {}

    for line in lines:
        record = line[:6]
        if record != 'ATOM  ':
            continue
        res_name = line[17:20]
        if res_name == 'SOL':
##            print lines[5+prev_atom_no-1]
##            stop
            atom_name = line[12:16]
            if atom_name != ' OW ':
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord_water = numpy.array([x,y,z,])
            ## loop over protein atoms
            for coord in l_coords:
                dist_sq = sum((coord-coord_water)**2)
                if dist_sq < 9:
                    atom_no = int(line[6:11])
                    print i, atom_no, dist_sq
                    if not atom_no in d_atom_nos:
                        d_atom_nos[atom_no] = []
                    d_atom_nos[atom_no] += [i]
                    break
        elif res_name not in ['SOL',' Cl',]:
            res_no = int(line[22:26])
            if res_no in [35,52,]:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coord = numpy.array([x,y,z,])
                l_coords += [coord]
                if x < 5 or x > 65 or y < 5 or y > 65 or z < 5 or z > 65:
                    print coord
                    stop_edge

print d_atom_nos

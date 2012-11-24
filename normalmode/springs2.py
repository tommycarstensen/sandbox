def main():

    '''write script for drawing cylinders between atoms in vmd'''

    import Numeric, math

    fd = open('/oxygenase_local/data/pdb/lz/pdb2lzm.ent','r')
    lines = fd.readlines()
    fd.close()

    d_coordinates = {}
    for i in range(164):
        d_coordinates[i+1] = {}
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            atom_name = line[12:16].strip()
            res_no = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coordinate = Numeric.array([x, y, z])
            d_coordinates[res_no][atom_name] = coordinate

    lines0 = []
    lines0 += ['draw color green\n']
    lines1 = []
    lines1 += ['draw color violet\n']
    lines2 = []
    lines2 += ['draw color orange\n']
    lines3 = []
    lines3 += ['draw color black\n']
    l_backbone = ['N','CA','C','O','OXT','CB',]
    for i in range(164):
        for j in range(i+1,164):
            sq_dist_ala_both = []
            sq_dist_ala_one = []
            sq_dist_ala_none = []
            for atom_name1 in d_coordinates[i+1].keys():
                c1 = d_coordinates[i+1][atom_name1]
                for atom_name2 in d_coordinates[j+1].keys():
                    c2 = d_coordinates[j+1][atom_name2]
                    sq_dist = sum((c1-c2)**2)
                    sq_dist_ala_none += [sq_dist]
                    if atom_name1 in l_backbone and atom_name2 in l_backbone:
                        sq_dist_ala_both += [sq_dist]
                    if atom_name1 in l_backbone or atom_name2 in l_backbone:
                        sq_dist_ala_one += [sq_dist]
                    if atom_name1 == 'CA' and atom_name2 == 'CA':
                        CA_dist = sq_dist
                        x1 = c1[0]
                        y1 = c1[1]
                        z1 = c1[2]
                        x2 = c2[0]
                        y2 = c2[1]
                        z2 = c2[2]
            line = 'draw cylinder {%f %f %f} {%f %f %f} radius 0.1\n' %(
                x1,y1,z1,x2,y2,z2,
                )
            if min(sq_dist_ala_both) < min(sq_dist_ala_one):
                print min(sq_dist_ala_both)
                print min(sq_dist_ala_one)
                stop
            if CA_dist < 36:
                lines0 += [line]
            if min(sq_dist_ala_none) < 36:
                lines1 += [line]
            if min(sq_dist_ala_both) < 36:
                lines2 += [line]
            if min(sq_dist_ala_none) < 36 and min(sq_dist_ala_both) > 36:
                lines3 += [line]

    fd = open('vmd0.vmd','w')
    fd.writelines(lines0)
    fd.close()
    fd = open('vmd1.vmd','w')
    fd.writelines(lines1)
    fd.close()
    fd = open('vmd2.vmd','w')
    fd.writelines(lines2)
    fd.close()
    fd = open('vmd3.vmd','w')
    fd.writelines(lines3)
    fd.close()

    return


if __name__=='__main__':
    main()

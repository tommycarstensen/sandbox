def main():

    import Numeric, math

    fd = open('/oxygenase_local/data/pdb/lz/pdb2lzm.ent','r')
    lines = fd.readlines()
    fd.close()

    l_coordinates = []
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coordinate = Numeric.array([x, y, z])
            l_coordinates += [coordinate]

    lines = []
    lines += ['draw color green\n']
    for i in range(len(l_coordinates)):
        c1 = l_coordinates[i]
        for j in range(i+1,len(l_coordinates)):
            c2 = l_coordinates[j]
            d = sum(((c2-c1)**2))
            if d < 100:
                x1 = c1[0]
                y1 = c1[1]
                z1 = c1[2]
                x2 = c2[0]
                y2 = c2[1]
                z2 = c2[2]
                line = 'draw cylinder {%f %f %f} {%f %f %f} radius 0.1\n' %(
                    x1,y1,z1,x2,y2,z2,
                    )
                lines += [line]

    fd = open('vmd.vmd','w')
    fd.writelines(lines)
    fd.close()

    return


if __name__=='__main__':
    main()

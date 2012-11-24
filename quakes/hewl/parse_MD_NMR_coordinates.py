def main(lines):

    import numpy

    d_residues = {
        'NLY':'LYS',
        'CLE':'LEU',
        'GLH':'GLU',
##        'ASH':'ASP',
        'HID':'HIS',
##        'HIE':'HIS',
##        'HIP':'HIS',
        }

    d_coords = {}
    for line in lines:

        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue

        res_name = line[17:20]
        if res_name == 'SOL':
            break
        if res_name in d_residues.keys():
            res_name = d_residues[res_name]

        element = line[76:78].strip()
        if element == 'H':
            continue

        atom_name = line[12:16].strip()
        if atom_name[0] == 'H':
            continue
        
        res_no = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        occupancy = float(line[54:60])
        tempFactor = float(line[60:66])
        if not res_no in d_coords.keys():
            d_coords[res_no] = {'res_name':res_name,'atoms':{},}
        d_coords[res_no]['atoms'][atom_name] = {'tempFactor':tempFactor,'occupancy':occupancy,'coord':coord,}

    return d_coords


if __name__ == '__main__':
    main()

## copy ligand from pdb2 to pdb1

import numpy
import sys
sys.path.append('/home/people/tc/svn/Protool/')
import geometry
instance_geometry = geometry.geometry()

pdb1 = sys.argv[sys.argv.index('-pdb1')+1]
pdb2 = sys.argv[sys.argv.index('-pdb2')+1]
chain1 = sys.argv[sys.argv.index('-chain1')+1]
chain2 = sys.argv[sys.argv.index('-chain2')+1]
pdb2 = '/data/pdb-v3.2/%s/pdb%s.ent' %(pdb2[1:3],pdb2,)


def main():

    l_coordinates1 = parse_protein_coordinates(pdb1,chain1,)
    l_coordinates2 = parse_protein_coordinates(pdb2,chain2,)

    if len(l_coordinates1) != len(l_coordinates2):
        print len(l_coordinates1)
        print len(l_coordinates2)
        stop

    rmsd = instance_geometry.superpose(l_coordinates1,l_coordinates2,)
    print 'rmsd', rmsd

    tv1 = instance_geometry.fitcenter
    rm = instance_geometry.rotation
    tv2 = instance_geometry.refcenter

    l_coords,lines = parse_ligand(pdb2,chain2,)

    lines = transform_ligand(tv1,rm,tv2,l_coords,lines,)

    add_ligand(pdb1,lines,pdb2,)

    return


def add_ligand(pdb1,lines_ligand,pdb2,):

    fd = open(pdb1,'r')
    lines = fd.readlines()
    fd.close()

    lines += lines_ligand

    fd = open('%s_w_%s_ligand.pdb' %(pdb1[:-4],pdb2[-8:-4],),'w')
    fd.writelines(lines)
    fd.close()

    return


def transform_ligand(tv1,rm,tv2,l_coords,lines,):

    for i in range(len(l_coords)):

        coord = l_coords[i]
        coord = numpy.dot(coord-tv1,rm)+tv2

        lines[i] = '%30s%8.3f%8.3f%8.3f%6s%6s%14s\n' %(
            lines[i][:30],coord[0],coord[1],coord[2],lines[i][54:60],'100.00',lines[i][66:80],
            )

    return lines


def parse_ligand(pdb,chain,):

    fd = open(pdb,'r')
    lines = fd.readlines()
    fd.close()

    l_coords = []
    l_lines = []

    for line in lines:
        record = line[:6].strip()
        if record == 'HETATM':
            hetID = line[17:20]
            if hetID == 'HOH':
                continue
            if line[21] != chain:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            l_coords += [coord]
            l_lines += [line]

    return l_coords, l_lines


def parse_protein_coordinates(pdb,chain,):

    fd = open(pdb,'r')
    lines = fd.readlines()
    fd.close()

    l_coords = []

    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            if line[21] != chain:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            l_coords += [coord]

    return l_coords


if __name__ == '__main__':
    main()

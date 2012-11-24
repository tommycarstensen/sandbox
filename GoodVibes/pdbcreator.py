#!/bin/env python
#
# $Id: pdbcreator.py 28 2006-07-10 11:52:24Z tc $
#
## all lengths are measured in Angstrom

def main(atoms_initial = 200., cubicangstrom_per_residue = 100., interdistance = 1.75):

    import math, random, goodvibes_autohingeidentification

    instance_vibration = goodvibes_autohingeidentification.vibration()

    ## general stuff
    volume_initial = (atoms_initial * cubicangstrom_per_residue)
    print volume_initial
    R = math.pow(volume_initial / ((4./3)*math.pi), 1./3)
    print R

    coordinates_sphere = []
    atom_sphere = 0
    dist_sphere_min = R**2
    while atom_sphere < atoms_initial:
        x = 2*R*random.random() - R
        y = 2*R*random.random() - R
        z = 2*R*random.random() - R
        dist_sphere_center = x**2+y**2+z**2
        if dist_sphere_center > R**2:
            continue
        dist_sphere_interatom_min = R
        for prev_res in coordinates_sphere:
            deltax = prev_res[0]-x
            deltay = prev_res[1]-y
            deltaz = prev_res[2]-z
            dist_sphere_interatom = deltax**2+deltay**2+deltaz**2
            if dist_sphere_interatom < dist_sphere_interatom_min:
                dist_sphere_interatom_min = dist_sphere_interatom
        if dist_sphere_interatom_min < (2*interdistance)**2:
            continue
        coordinates_sphere.append([x,y,z,dist_sphere_center])
        atom_sphere += 1

##        ## write coordinates_sphere to file
##        lines_sphere = []
##        for atom in range(len(coordinates_sphere)):
##            line = 'ATOM  %5i CA   GLY A%4i    %8.3f%8.3f%8.3f      %6.2f              \n' %(atom, atom, coordinates_sphere[atom][0], coordinates_sphere[atom][1], coordinates_sphere[atom][2], coordinates_sphere[atom][3])
##            lines_sphere.append(line)
##        fd = open('triangle1.pdb', 'w')
##        fd.writelines(lines_sphere)
##        fd.close()

        ## read coordinates_sphere from file
        coordinates_sphere = []
        fd = open('triangle1.pdb', 'r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            coordinates_sphere.append([float(line[30:38]), float(line[38:46]), float(line[46:54]), 0])

    bandwidth = 2*interdistance
    ## loops applicable to triangular pocket
    for angle in range(0,180,180/60):
        if angle == 0:
            continue
        angle = float(angle)

        for depth in range(0,int(2*R),1):
            if depth == 0:
                continue
            depth = float(depth)

            print angle, depth

            ## remove coordinates from sphere
            coordinates_write = list(coordinates_sphere)
            for i in range(len(coordinates_write)-1, -1, -1):
                x = coordinates_write[i][0]
                y = coordinates_write[i][1]
                bandheight = bandwidth/math.sin(math.radians(angle/2))
                if y > R-depth and x < (y-(R-depth))*math.tan(math.pi*angle/360.) and x > -(y-(R-depth))*math.tan(math.pi*angle/360.):
                    del coordinates_write[i]

            ## prepare modified coordinates as input
            lines_job = []
            for atom in range(len(coordinates_write)):
                line = 'ATOM  %5i CA   GLY A%4i    %8.3f%8.3f%8.3f      %6.2f              \n' %(atom, atom, coordinates_write[atom][0], coordinates_write[atom][1], coordinates_write[atom][2], coordinates_write[atom][3])
                lines_job.append(line)

            job = 'triangle_%s%s_%s%s.pdb' %('angle', angle, 'depth', depth)
            results = instance_vibration.main(
                atoms_initial, lines_job, [], None,
                ['CA'], float(1), [6,7,8,9,10,11,12,13,14,15,16,18,20,24,28,36], 3, 'monomeric', job, 50,
                None, None, None, 
                )

    return

if __name__=='__main__':
    main()
##    sheet()

##    for angle in range(10., 180., 18.):
##        for depth in range(.1*2*r, 1.*2*r, .1*2*r):
##            chord = 2 * math.sqrt(radius**2 - (radius-depth)**2)
##            volume_cap = (1./6.) * math.pi * (3*chord**2 + depth**2) * depth
##            volume_cone = (1./3.) * math.pi * radius**2 * (radius-depth)
##            angle = .5 * math.asin(.5*chord / radius)
##            print volume_cap, volume_cone
##            h = radius-math.sqrt(radius**2-(radius*math.sin(.5*angle))**2)
##
##            cuboidz = math.sqrt(radius**2 - (.5*xlen)**2)
##            volume_cuboid = xlen*(2*radius-2*ylen)*cuboidz
##            sphericalcapzh = radius - cuboidz
##            volume_sphericalcapz = (1./6.) * math.pi * (3*(.5*xlen)**2 + sphericalcapzh**2) * sphericalcapzh ## 2 af dem

##    ## loops applicable to spheric pocket
##    for r in range(0.,R,3.):
##        if r == 0:
##            continue
##        for d in range(R-r,R+r,3.):
##            if d == R-r:
##                continue
##            print R,r,d
##            v = volume_initial - (1./(12.*d)) * math.pi * (R+r-d)**2 * (d**2 + 2*d*r - 3*r**2 + 2*d*R + 6*r*R - 3*R**2)
##           
##            atoms = round(v / cubicangstrom_per_residue)

##    for ylen_percentage in range(26., 100., 2.):
##        print ylen_percentage
##        for xlen_percentage in range(26., 30., 1.):
##            print xlen_percentage
##
##            xlen = (xlen_percentage / 100.) * 2*radius
##            ylen = (ylen_percentage / 100.) * 2*radius
##            print xlen
##            print ylen
##            
##            if radius - math.sqrt(radius**2 - xlen**2) > ylen:
##                continue
##            if radius + math.sqrt(radius**2 - xlen**2) < ylen:
##                continue
##    
##            chordy = 2 * math.sqrt(radius**2 - (.5*xlen)**2)
##            volume_spheresegmenty = (1./6.) * math.pi * (3*radius**2 + 3*chordy**2 + (.5*xlen)**2)*(.5*xlen) ## 2 af dem
##            print 2*volume_spheresegmenty
##
##            cuboidz = math.sqrt(radius**2 - (.5*xlen)**2)
##            volume_cuboid = xlen*(2*radius-2*ylen)*cuboidz
##            sphericalcapzh = radius - cuboidz
##            volume_sphericalcapz = (1./6.) * math.pi * (3*(.5*xlen)**2 + sphericalcapzh**2) * sphericalcapzh ## 2 af dem
##            print volume_cuboid
##            print 2*volume_sphericalcapz
##
##            volume_new = volume - (1./2.) * (2*volume_spheresegmenty - 2*volume_sphericalcapz - volume_cuboid)
##            print 'volume', volume, 'A^3'


##def sheet(strand_atomdist_intra = 1.75, strand_atomdist_inter = 1.75, strand_len = 5, strands = 8, strand_angle_inter = 0, circle_angles = 180):
##    for strand in range(strands):
##        y = 
##
##    coordinates_sheet = []
##    for strand in range(strands):
##        y = strand*strand_atomdist_inter*2
##        for residue in range(strand_len):
##            x = residue*strand_atomdist_intra*2
##            coordinates_sheet.append([x, y, 0])
##
##    lines_sheet = []
##    for atom in range(len(coordinates_sheet)):
##        line = 'ATOM  %5i CA   GLY A%4i    %8.3f%8.3f%8.3f                          \n' %(atom, atom, coordinates_sheet[atom][0], coordinates_sheet[atom][1], coordinates_sheet[atom][2])
##        lines_sheet.append(line)
##    fd = open('sheet.pdb', 'w')
##    fd.writelines(lines_sheet)
##    fd.close()
##    
##    return



#####################################cuboid#####################################

##            print 'residues', atoms
##
##            coordinates_sphere = []
##            coordinates_cube = []
##            atom_sphere = 0
##            atom_cube = 0
##            dist_sphere_min = R**2
##            while atom_sphere < atoms:
####                print atom_sphere, atoms
##                if atoms-atom_sphere < 5:
##                    print atom_sphere
##                x = 2*R*random.random() - R
##                y = 2*R*random.random() - R
##                z = 2*R*random.random() - R
##                dist_sphere_center = x**2+y**2+z**2
##                if dist_sphere_center > R**2:
##                    continue
##                if (x-d)**2+y**2+z**2 < r**2:
##                    continue
##                dist_sphere_interatom_min = R
##                for prev_res in coordinates_sphere:
##                    deltax = prev_res[0]-x
##                    deltay = prev_res[1]-y
##                    deltaz = prev_res[2]-z
##                    dist_sphere_interatom = deltax**2+deltay**2+deltaz**2
##                    if dist_sphere_interatom < dist_sphere_interatom_min:
##                        dist_sphere_interatom_min = dist_sphere_interatom
##                if atom_cube < atoms:
##                    coordinates_cube.append([x,y,z])
##                    atom_cube += 1
##                if dist_sphere_interatom_min < (2*interdistance)**2:
##                    continue
##                ## continue statement that applies to the cuboid only
####                if  (x > -xlen and x < xlen and y > radius-ylen):
####                    continue
##                coordinates_sphere.append([x,y,z,dist_sphere_center])
##                atom_sphere += 1

##            lines_sphere = []
##            for atom in range(len(coordinates_write)):
##                line = 'ATOM  %5i CA   GLY A%4i    %8.3f%8.3f%8.3f      %6.2f              \n' %(atom, atom, coordinates_write[atom][0], coordinates_write[atom][1], coordinates_write[atom][2], coordinates_write[atom][3])
##                lines_sphere.append(line)
##            fd = open('triangle_%s%s_%s%s.pdb' %('angle', angle, 'depth', depth), 'w')
##            fd.writelines(lines_sphere)
##            fd.close()

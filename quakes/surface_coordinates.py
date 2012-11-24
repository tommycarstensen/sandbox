import numpy

d_dimensions = {
    'x':{
        'min': 999.999,
        'max':-999.999,
        },
    'y':{
        'min': 999.999,
        'max':-999.999,
        },
    'z':{
        'min': 999.999,
        'max':-999.999,
        },
    }

fd = open('1E8L.pdb','r')
lines = fd.readlines()
fd.close()

l_coords_protein = []
for line in lines:
    record = line[:6].strip()
    if record == 'ATOM':
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        l_coords_protein += [coord]
        for coord,k_coord in [
            [x,'x',],
            [y,'y',],
            [z,'z',],
            ]:
            if coord < d_dimensions[k_coord]['min']:
                d_dimensions[k_coord]['min'] = coord
            if coord > d_dimensions[k_coord]['max']:
                d_dimensions[k_coord]['max'] = coord

dist_max = 4
dist_max_sq = dist_max**2
dist_min = .5
dist_min_sq = dist_min**2

i = 0
l_coords_solvent = []
for x in range(int(d_dimensions['x']['min'])-2,int(d_dimensions['x']['max'])+2+1,1):
    for y in range(int(d_dimensions['y']['min'])-2,int(d_dimensions['y']['max'])+2+1,1):
        for z in range(int(d_dimensions['z']['min'])-2,int(d_dimensions['z']['max'])+2+1,1):

            bool_vicinal = False
            for coord_protein in l_coords_protein:

                bool_distant = False

                for i_coord in range(3):
                    if abs(coord_protein[i_coord]-[x,y,z,][i_coord]) > dist_max:
                        bool_distant = True
                        break

                if bool_distant == True:
                    continue

                dist_sq = (
                    (x-coord_protein[0])**2
                    +
                    (y-coord_protein[1])**2
                    +
                    (z-coord_protein[2])**2
                    )
                if dist_sq < dist_min_sq:
                    bool_vicinal = True
                    break

            i += 1
            if i % 20 == 0:
                print i
            if bool_vicinal == False and bool_distant == False:
                l_coords_solvent += [numpy.array([x,y,z,])]

            
print len(l_coords_solvent)
print len(l_coords_protein)
print l_coords_solvent

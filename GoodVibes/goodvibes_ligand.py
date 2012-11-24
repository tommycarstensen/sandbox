##!/bin/env /usr/bin/python
##
##$Id$
##
##Tommy Carstensen, University College Dublin, December 2009


import numpy, math, sys, optparse, os
sys.path.append('/home/tc/svn/Protool/')
import geometry
import NMA

##
## settings
##

## cutoff distance (Angstrom)
dist_co = 6
##l_atoms = ['N','CA','C','O',]
l_atoms = ['CA',]
## shell size (Angstrom)
shell = 4
## grid size (Angstrom)
grid = 2

def main(
    pdb,chain,dist_max,dist_min,mode='single',v_apoholo=None,l_coords_probe=None,
    l_coords_protein_alpha=None,
    ):

    ##
    ## settings
    ##
    dist_min_sq = dist_min**2
    dist_max_sq = dist_max**2

    ## parse coordinates
    d_coords = parse_pdb_coordinates(pdb,chain,)

    if l_coords_protein_alpha == None:
        ## parse alpha carbon atoms
        l_coords_protein_alpha = parse_alpha_carbon_atoms(d_coords,)

    ## calulate hessian matrix
    matrix_hessian_protein = do_interactions(l_coords_protein_alpha,)
    ## diagonalize hessian matrix
    eigenvectors_protein, eigenvalues_protein = NMA.diagonalize_hessian(matrix_hessian_protein,)

    if v_apoholo != None:
        mode_max_apoholo, overlap_max_apoholo, l_factors = find_max_mode_apo_holo(
            pdb,eigenvectors_protein,v_apoholo,
            eigenvalues_protein,
            )
##    ## tmp!!!
##    mode_max_apoholo = 6
##    v1 = v_apoholo
##    v2 = eigenvectors_protein[mode_max_apoholo]
##    overlap_max_apoholo = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))

##    if 2+2 == 4: ## tmp!!!
##        print 'tmp!!!'
##        return l_factors

    ## determine dimensions of protein
    d_dimensions = determine_protein_dimensions(l_coords_protein_alpha,)
    ## add probe atoms
    fn = '/home/tc/UCD/GV_ligand_binding_site_identification/output/GoodVibes/distmax6_distmin3/%s_%s_probe.pdb' %(pdb,chain,)
    if l_coords_probe:
        print 'a'
        pass
    elif os.path.isfile(fn):
        print 'b'
        l_coords_probe = []
        fd = open(fn)
        lines = fd.readlines()
        fd.close()
        for line in lines:
            record = line[:6].strip()
            if record == 'HETATM' and line[17:20] == 'EXT':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coord = numpy.array([x,y,z,])
                l_coords_probe += [coord]
    else:
        l_coords_probe = add_probe_atoms(d_coords,d_dimensions,dist_min_sq,dist_max_sq,)

    ## calculate overlaps
    print 'looping over', len(l_coords_probe), 'probe coordinates'
    l_overlaps = []

    for i in range(len(l_coords_probe)):

        print i, len(l_coords_probe)

        coord_holo = l_coords_probe[i]

        l_coords = l_coords_protein_alpha+[coord_holo]

##        matrix_hessian_holo = do_interactions(l_coords,bool_extra=True,) ## tmp!!!
##        matrix_hessian_holo = do_interactions(l_coords,bool_strong=True) ## tmp!!!
        matrix_hessian_holo = do_interactions(l_coords)

        try:
            eigenvectors_holo, eigenvalues_holo = NMA.diagonalize_hessian(matrix_hessian_holo)
        except:
            print 'exception'
            l_overlaps += [1.]
            continue

        ## compare to x-ray motion
        if v_apoholo != None:
            v1 = v_apoholo
            v2 = eigenvectors_holo[mode_max_apoholo][:-3]
            overlap = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
            print 'overlap', overlap

            max_overlap = overlap
            max_mode = mode_max_apoholo

            ## check neigboring modes for max overlap...
            if 2+2 == 5:
                max_mode = mode_max_apoholo
                switch_max = 3
                for mode in range(max(6,mode_max_apoholo-switch_max),mode_max_apoholo+switch_max+1):
                    if mode == mode_max_apoholo:
                        continue
                    v2 = eigenvectors_holo[mode][:len(v_apoholo)]
                    overlap = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
                    if overlap > max_overlap:
                        print '********', mode, round(overlap,3), mode_max_apoholo, round(max_overlap,3), mode_max_apoholo, round(overlap_max_apoholo,3), pdb
                        max_overlap = overlap
                        max_mode = mode
                        if mode_max_apoholo < 12 and overlap > 1.2*overlap_max_apoholo:
                            print '******** induced fit?'
            l_overlaps += [max_overlap]

            ## perturb elastic netwrok and recalculate mode contribution
            if 2+2 == 5:
                eigenvectors_holo = numpy.transpose(eigenvectors_holo)
                vector = numpy.array([0.,0.,0.,])
                v_apoholo = numpy.array(list(v_apoholo)+[0.,0.,0.,])
                l_factors_holo = numpy.linalg.solve(eigenvectors_holo,v_apoholo,)
                l_factors_holo_abs = [abs(factor) for factor in l_factors_holo]

                if mode_max_apoholo != list(l_factors_holo_abs).index(max(l_factors_holo_abs)):
                    print mode_max_apoholo, list(l_factors_holo_abs).index(max(l_factors_holo_abs))
                    print mode_max_apoholo, overlap_max_apoholo, overlap
                    print l_factors_holo_abs[mode_max_apoholo], max(l_factors_holo)
                    s = '# mode factor absfactor eigenvalue\n'
                    for i in range(len(l_factors_holo)):
                        s += '%s %s %s\n' %(i+1, l_factors_holo[i],abs(l_factors_holo[i]),)
                    fd = open('facs_eigvals_%s_perturbed.txt' %(pdb),'w')
                    fd.write(s)
                    fd.close()
                    write_pdb(l_overlaps,l_coords_probe,pdb,chain,)
    ##                stop_mode

##            ## tmp!!!
##            v2 = eigenvectors_holo[6][:-3]
##            overlap6 = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
##            if overlap6 > 1.1*overlap:
##                print mode_max_apoholo, overlap
##                print 6, overlap6
##                v2 = eigenvectors_holo[7][:-3]
##                overlap7 = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
##                print 7, overlap7
##                v2 = eigenvectors_holo[8][:-3]
##                overlap8 = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
##                print 8, overlap8
##                stop
            
        elif mode == 'single':
            eigenvectors_holo = eigenvectors_holo[:-3]
            l = []
            ## check first 3 modes in case eigenvalues have swapped
            for mode_holo in range(6,10,):
                overlap = calc_overlap(
                    eigenvectors_protein, eigenvectors_holo,
##                    eigenvalues_protein, eigenvalues_holo,
                    mode_holo = mode_holo,
                    )
                l += [overlap]
                if overlap > 0.9:
                    break
            overlap_max = max(l)

            print pdb, i, len(l_coords_probe), overlap_max

##            ## go for mode 7
##            if overlap_max < 0.9:
##                overlap_max = l[0]
            l_overlaps += [overlap_max]
##            if overlap_max < 0.90:
##                print pdb, i+1, len(l_coords_probe), overlap
##                print calc_overlap(
##                    eigenvectors_protein,eigenvectors_holo,
##                    eigenvalues_protein, eigenvalues_holo,
##                    mode_holo = 6,
##                    )
##                stop
        elif mode == 'multiple':
            eigenvectors_holo = eigenvectors_holo[:-3]
            overlap = calc_overlap(
                eigenvectors_protein,eigenvectors_holo,
                eigenvalues_protein, eigenvalues_holo,
                l_factors = l_factors,
                )
            l_overlaps += [overlap]
            print overlap, i, len(l_coords_probe)
        else:
            print sys.argv
            stop
            
##    fd = open('l_overlaps.txt','r')
##    s = fd.read()
##    fd.close()
##    l_overlaps = s.split()
##    l_overlaps = l_overlaps[1::2]

    ##
    ## combine protein and probe coordinates and add bfactors
    ##

##    d_coords_holo = parse_pdb_coordinates(pdb_holo,chain_holo,)
##    if len(d_coords.keys()) != len(d_coords_holo.keys()):
##        print len(d_coords.keys())
##        print len(d_coords_holo.keys())
##        stop
##    l_coords_protein_alpha_holo = parse_alpha_carbon_atoms(d_coords_holo,)
##    instance_geometry = geometry.geometry()
##    rmsd = instance_geometry.superpose(l_coords_protein_alpha,l_coords_protein_alpha_holo,)
##    tv1 = instance_geometry.fitcenter
##    rm = instance_geometry.rotation
##    tv2 = instance_geometry.refcenter
##    parse_ligand_coordinates(pdb_holo,chain_holo,ligand_ID,)

    if v_apoholo != None and len(l_overlaps) > 1:
        print l_overlaps
        l_overlaps = fix_overlaps(l_overlaps)
        print max(l_overlaps), min(l_overlaps)

    if (
        v_apoholo == None
        or
        (v_apoholo != None and len(l_overlaps) > 1)
        ):
        write_pdb(l_overlaps,l_coords_probe,pdb,chain,)

    if v_apoholo != None:
        d = {
            'mode_max_apoholo':mode_max_apoholo,
            'overlap_max_apoholo':overlap_max_apoholo,
            'l_overlaps':l_overlaps,
            'l_factors':l_factors,
            'eigenvectors':eigenvectors_protein,
            }
        if 2+2 == 5:
            d['l_factors_probe'] = l_factors_holo
            d['max_mode'] = max_mode
            
        return d
    else:
        print 'how much to return to function that called me? just l_overlaps?'
        return l_overlaps


def find_max_mode_apo_holo(
    pdb,eigenvectors_protein,v_apoholo,
    eigenvalues_protein,
    ):

    eigenvectors_protein = numpy.transpose(eigenvectors_protein)
    l_factors = numpy.linalg.solve(eigenvectors_protein,v_apoholo,)
    eigenvectors_protein = numpy.transpose(eigenvectors_protein)

    if max(l_factors[:6]) > 0.001:
        print l_factors[:7]
        stop

#        x,residues,rank,s = numpy.linalg.lstsq(eigenvectors_protein,v_apoholo,)

    l_factors_abs = [abs(factor) for factor in l_factors]
    max_mode = list(l_factors_abs).index(max(l_factors_abs))
    print 'max contributing mode', max_mode
    if max_mode < 6:
        print max_mode
        stop

##    s = '# mode factor absfactor eigenvalue\n'
##    for i in range(len(l_factors)):
##        s += '%s %s %s %s\n' %(i+1, l_factors[i],abs(l_factors[i]),eigenvalues_protein[i],)
##    fd = open('facs_eigvals_%s.txt' %(pdb),'w')
##    fd.write(s)
##    fd.close()

    v1 = v_apoholo
    v2 = eigenvectors_protein[max_mode]
    overlap = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
    print 'overlap of max mode', overlap
    v2 = eigenvectors_protein[6]
    overlap6 = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
    print 'overlap of mode 6', overlap6

    return max_mode, overlap, l_factors


##def parse_ligand_coordinates(pdb,chain,ligand_ID,):
##
##    fd = open('/data/pdb-v3.2/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
##    lines = fd.readlines()
##    fd.close()
##
##    d_coordinates = {}
##    for line in lines:
##        record = line[:6].strip()
##
##
##    return


def fix_overlaps(l_overlaps):

    '''fix overlaps so b factors range from 0 to 100'''

    min_overlap = min(l_overlaps)
    max_overlap = max(l_overlaps)
    for i_overlap in range(len(l_overlaps)):
        overlap = l_overlaps[i_overlap]
        overlap = (overlap-min_overlap)/(max_overlap-min_overlap)
        l_overlaps[i_overlap] = overlap

    return l_overlaps


def calc_overlap(
    eigenvectors_protein,eigenvectors_extra,
    mode_holo=None,
    l_factors = None,
    ):

    if mode_holo:

        mode_apo = 6
        v1 = eigenvectors_protein[mode_apo]
        v2 = eigenvectors_extra[mode_holo][:-3]
        overlap = (
            abs(numpy.dot(v1,v2))
            /
            math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
            )

##    else:
##
##        v1 = eigenvectors_protein[6]
##        v2 = eigenvectors_extra[6]
##
##        for i_vector in range(2):
##
##            v,eigenvalues,eigenvectors = [
##                [v1,eigenvalues_protein,eigenvectors_protein,],
##                [v2,eigenvalues_extra,eigenvectors_extra,],
##            ][i_vector]
##
##            ## loop over modes
##            for i in range(7,len(eigenvectors)-1):
##
##                ## scale length of mode i relative to length of mode 7
##                len7 = 1
##                leni = 1
##                if l_factors:
##                    lenfactor = l_factors[i]
##                else:
##                    lenfactor = (len7/leni)/(eigenvalues[i]/eigenvalues[6])
####                v += eigenvectors[i]*lenfactor
##                eigenvectors[i] *= lenfactor
##                if lenfactor > 1:
##                    print lenfactor
##                    stop
        
    return overlap


def parse_alpha_carbon_atoms(d_coords,):

    l_coords_protein_alpha = []
    l_atom_nos = d_coords.keys()
    l_atom_nos.sort()
    
    for atom_no in l_atom_nos:
        d = d_coords[atom_no]
        atom_name = d['atom_name']
        coord = d['coord']
        if atom_name == 'CA':
            l_coords_protein_alpha += [coord]

    return l_coords_protein_alpha


def determine_protein_dimensions(l_coords_alpha,):

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

    for coord_set in l_coords_alpha:

        x = coord_set[0]
        y = coord_set[1]
        z = coord_set[2]

        for coord,k_coord in [
            [x,'x',],
            [y,'y',],
            [z,'z',],
            ]:
            if coord < d_dimensions[k_coord]['min']:
                d_dimensions[k_coord]['min'] = coord
            if coord > d_dimensions[k_coord]['max']:
                d_dimensions[k_coord]['max'] = coord

    return d_dimensions


def do_interactions(l_coords,bool_extra=False,bool_strong=False):

    N = 3*len(l_coords)
    matrix_hessian = numpy.zeros((N,N), dtype=float)

    row_sup = 0
    col_sup = 0
    for row_sup in range(len(l_coords)):
        xi = l_coords[row_sup][0]
        yi = l_coords[row_sup][1]
        zi = l_coords[row_sup][2]
        for col_sup in range(row_sup+1,len(l_coords)):
            xj = l_coords[col_sup][0]
            yj = l_coords[col_sup][1]
            zj = l_coords[col_sup][2]

            vx = xj-xi
            vy = yj-yi
            vz = zj-zi
            vector = [vx,vy,vz,]

            dist_sq = vx**2+vy**2+vz**2
            dist = math.sqrt(dist_sq)

##            ## make interaction stronger to "ligand"
##            if bool_strong == True and (row_sup >= len(l_coords)-1-3 or col_sup >= len(l_coords)-1-3):
##                dist /= 3
##                dist_sq /= 9

            factor = do_sigmoid_factor(dist)

##            ## make interaction stronger to "ligand"
##            if bool_strong == True and (row_sup >= len(l_coords)-1-3 or col_sup >= len(l_coords)-1-3):
##                factor *= 10
            
            for row_sub in range(3):
                for col_sub in range(3):

                    if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                        if dist_sq == 0:
                            print row_sub, col_sub, row_sup, col_sup
                            stop
                        value = factor*-vector[row_sub]*vector[col_sub]/dist_sq

                        matrix_hessian[3*row_sup+row_sub,3*col_sup+col_sub] = value ##upper super off-diagonal; xixj, xiyj, xizj, yiyj, yizj, zizj
                        matrix_hessian[3*col_sup+col_sub,3*row_sup+row_sub] = value ##lower super off-diagonal; xjxi, yjxi, zjxi, yjyi, zjyi, zjzi
                        matrix_hessian[3*row_sup+row_sub,3*row_sup+col_sub] -= value ##super diagonal (row); xixi, xiyi, xizi, yiyi, yizi, zizi
                        matrix_hessian[3*col_sup+col_sub,3*col_sup+row_sub] -= value ##super diagonal (col); xjxj, yjxj, zjxj, yjyj, zjyj, zjzj
                        if col_sub > row_sub: #fill lower subsymmetrical elements
                            matrix_hessian[3*row_sup+col_sub,3*col_sup+row_sub] = value #upper super off-diagonal; yixj, zixj, ziyj
                            matrix_hessian[3*col_sup+row_sub,3*row_sup+col_sub] = value #lower super off-diagonal; xjyi, xjzi, yjzi
                            matrix_hessian[3*row_sup+col_sub,3*row_sup+row_sub] -= value ##super diagonal; yixi, zixi, ziyi
                            matrix_hessian[3*col_sup+row_sub,3*col_sup+col_sub] -= value ##super diagonal; yjxj, zjxj, zjyj

    if bool_extra == True:
        matrix_hessian_new = numpy.zeros((N-3,N-3), dtype=float)
        for i in range(N-3):
            for j in range(N-3):
                matrix_hessian_new[i][j] = matrix_hessian[i][j]
        matrix_hessian = matrix_hessian_new


    return matrix_hessian


def do_sigmoid_factor(x):

    slope = 1
    y = 1. / ( 1. + math.exp( slope*(x-dist_co) ) )

    return y


def parse_pdb_coordinates(pdb,chain,):

    fd = open('/media/WDMyBook1TB/2TB/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()

    l_modres = []

    d_coordinates = {}
    for line in lines:
        record = line[:6].strip()

        if record in ['HETATM','ATOM',]:
            atom_no = int(line[6:11])
            atom_name = line[12:16].strip()
            altloc = line[16]
            chain_atom = line[21]
            res_no = int(line[22:26])
            iCode = line[26]
            if chain_atom != chain and chain_atom not in chain:
                continue
            if record == 'HETATM' and '%1s%4i%1s' %(chain_atom,res_no,iCode,) not in l_modres:
                continue
            if altloc not in [' ','A','1',]:
                print line
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            d_coordinates[atom_no] = {
                'atom_name':atom_name,
                'coord':coord,
                }

        elif record == 'MODRES':
            modres_chain = line[16]
            res_no = int(line[18:22])
            iCode = line[22]
            if modres_chain != chain:
                continue
            l_modres += ['%1s%4i%1s' %(modres_chain,res_no,iCode,)]

    return d_coordinates
            

def add_probe_atoms(d_coords,d_dimensions,dist_min_sq,dist_max_sq,):

    dist_max = math.sqrt(dist_max_sq)
    dist_min = math.sqrt(dist_min_sq)

    bool_test = True
    bool_test = False

    l_atom_nos = d_coords.keys()
    l_atom_nos.sort()

##    ## avoid placing probe next to terminal mobile residues (should be by residue number)
##    l_atom_nos = l_atom_nos[20:-20]

    l_coords_probe = []
    l_pdb_coords_solvent = []

    x_min = int(d_dimensions['x']['min'])-shell
    x_max = int(d_dimensions['x']['max'])+shell+1
    y_min = int(d_dimensions['y']['min'])-shell
    y_max = int(d_dimensions['y']['max'])+shell+1
    z_min = int(d_dimensions['z']['min'])-shell
    z_max = int(d_dimensions['z']['max'])+shell+1


##    ##
##    ## faster method but doesnt calculate angle
##    ##
##    d_class = {}
##    for x in range(x_min,x_max,grid):
##        d_class[x] = {}
##        for y in range(y_min,y_max,grid):
##            d_class[x][y] = {}
##            for z in range(z_min,z_max,grid):
##                d_class[x][y][z] = None
##    for k,v in d_coords.items():
##
##        atom_name_protein = v['atom_name']
##        if atom_name_protein[0] == 'H':
##            continue
##
##        coord = v['coord']
##
##        x_min_local = int(coord[0]-coord[0]%grid+x_min%grid-dist_max-grid)
##        x_max_local = int(coord[0]-coord[0]%grid+x_min%grid+dist_max+grid)
##        y_min_local = int(coord[1]-coord[1]%grid+y_min%grid-dist_max-grid)
##        y_max_local = int(coord[1]-coord[1]%grid+y_min%grid+dist_max+grid)
##        z_min_local = int(coord[2]-coord[2]%grid+z_min%grid-dist_max-grid)
##        z_max_local = int(coord[2]-coord[2]%grid+z_min%grid+dist_max+grid)
##        x_max_local = min(x_max,x_max_local)
##        y_max_local = min(y_max,y_max_local)
##        z_max_local = min(z_max,z_max_local)
##        x_min_local = max(x_min,x_min_local)
##        y_min_local = max(y_min,y_min_local)
##        z_min_local = max(z_min,z_min_local)
##
##        for x in range(x_min_local,x_max_local,grid):
##            for y in range(y_min_local,y_max_local,grid):
##                for z in range(z_min_local,z_max_local,grid):
##                    ## if too close once, then always too close
##                    if d_class[x][y][z] == 'tooclose':
##                        continue
##                    coord_grid = numpy.array([x,y,z])
##                    dist_sq = sum((coord-coord_grid)**2)
##                    if dist_sq > dist_max_sq:
##                        d_class[x][y][z] = 'toofar'
####                        print k, 'a'
##                    elif dist_sq < dist_min_sq:
##                        d_class[x][y][z] = 'tooclose'
####                        print k, 'b'
##                    else:
##                        d_class[x][y][z] = 'extra'
##    l_coords = []
##    for x in range(x_min,x_max,grid):
##        for y in range(y_min,y_max,grid):
##            for z in range(z_min,z_max,grid):
##                if d_class[x][y][z] == 'extra':
##                    coord = numpy.array([x,y,z])
##                    l_coords += [coord]
##    print len(l_coords)
##    l_overlaps = [0]*len(l_coords)
##    write_pdb(l_overlaps,l_coords,'2exo','A',suffix='_test_distmin%s_distmax%s' %(dist_min,dist_max,))
##    stop


    d_class = {}
    for x in range(x_min,x_max,grid):
        print x, x_min, x_max
        d_class[x] = {}
        for y in range(y_min,y_max,grid):
            d_class[x][y] = {}
            for z in range(z_min,z_max,grid):
                d_class[x][y][z] = {}

                bool_vicinal = False
                bool_not_distant = False
##                bool_not_distant = 0

                l_coords_not_distant = []
                for atom_no_protein in l_atom_nos:

                    d = d_coords[atom_no_protein]
                    atom_name_protein = d['atom_name']
                    coord_protein = d['coord']

                    ## skip if hydrogen atoms
                    if atom_name_protein[0] == 'H':
                        continue

##                    if not atom_name_protein in l_atoms:
##                        continue

                    dist_sq = (
                        (x-coord_protein[0])**2
                        +
                        (y-coord_protein[1])**2
                        +
                        (z-coord_protein[2])**2
                        )
                    if dist_sq < dist_min_sq:
                        bool_vicinal = True
##                    if dist_sq < dist_max_sq:
##                    if atom_name_protein == 'CA' and dist_sq < dist_max_sq:
##                    if atom_name_protein in ['N','CA','CB','C','O',] and dist_sq < dist_max_sq:
                    if atom_name_protein in l_atoms and dist_sq < dist_max_sq:
                        bool_not_distant = True
                        ## only do angle between alpha carbon atoms
                        if atom_name_protein == 'CA':
                            l_coords_not_distant += [coord_protein]
##                        bool_not_distant += 1
                    if bool_vicinal == True:
                        break

                if bool_vicinal == False:
                    if bool_not_distant == True:

##                        if bool_test != True:
##
##                            l_coords_probe += [numpy.array([x,y,z,])]
##                            d_class[x][y][z]['class'] = 'possible_pocket'
##
##                        else:
                            
                            d_class[x][y][z]['class'] = 'not_pocket'
                            if len(l_coords_not_distant) > 1:
                                for i in range(len(l_coords_not_distant)-1):
                                    coord1 = l_coords_not_distant[i]
                                    v1 = coord1-numpy.array([x,y,z,])
                                    for j in range(i+1,len(l_coords_not_distant)):
                                        coord2 = l_coords_not_distant[j]
                                        v2 = coord2-numpy.array([x,y,z,])
                                        
                                        ## angle larger than xxx'
                                        if numpy.dot(v1,v2) < 0: ## angle > 90
                                            d_class[x][y][z]['class'] = 'possible_pocket'
                                            l_coords_probe += [[x,y,z,]]
                                            break
                                    if d_class[x][y][z]['class'] == 'possible_pocket':
                                        break
                                    
                    elif bool_not_distant == False:
                        d_class[x][y][z]['class'] = 'distant'
                elif bool_vicinal == True:
                    d_class[x][y][z]['class'] = 'protein'

    if bool_test == True:

        sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
        import combinatorics
        l_translations = combinatorics.permutation_w_rep([-1,0,1,],3)
    ##    l_translations = combinatorics.permutation_w_rep([-2,-1,0,1,2,],5)
        l_translations.remove([0,0,0,])

        l_overlaps = []

        d_class_grid = {}
        for x in range(int(d_dimensions['x']['min'])-shell,int(d_dimensions['x']['max'])+shell+1,grid):
            print x
            d_class_grid[x] = {}
            for y in range(int(d_dimensions['y']['min'])-shell,int(d_dimensions['y']['max'])+shell+1,grid):
                d_class_grid[x][y] = {}
                for z in range(int(d_dimensions['z']['min'])-shell,int(d_dimensions['z']['max'])+shell+1,grid):
                    count_pocket = 1 ## point itself is pocket
                    count_protein = 0
                    count_distant = 0

                    if d_class[x][y][z]['class'] != 'possible_pocket':
                        continue

                    for translation in l_translations:
                        try:
                            Class = d_class[x+translation[0]*grid][y+translation[1]*grid][z+translation[2]*grid]['class']
                            if Class == 'possible_pocket':
                                count_pocket += 1
                            elif Class == 'protein':
                                count_protein += 1
                            elif Class == 'distant':
                                count_distant += 1
                        except:
                            None

                    count_pocket_min = 16
                    count_protein_min = 18
                    count_distant_max = 27
                    if (count_pocket >= count_pocket_min or count_protein >= count_protein_min) and count_distant <= count_distant_max:
                        print count_distant
                        ## near pocket and semi buried in protein (not surface)
                        if   count_pocket >= count_pocket_min and count_protein >  count_protein_min/2.:
                            count_pocket = 13 ## green, include
                        ## not vicinal to many pocket grid points, but very buried in protein
                        elif count_protein >= count_protein_min:
                            count_pocket = 27 ## blue, include
                        ## not near other pocket points
                        elif count_pocket <= count_pocket_min/2.:
                            count_pocket = 24 ## light blue, exclude
    ##                    ## surface of protein
    ##                    elif count_protein <= count_protein_min/2.:
    ##                        count_pocket = 6 ## orange, exclude
                        ## not buried = in a large pocket = possibly ligand binding site
                        else:
                            count_pocket = 21 ## cyan
                    else:
    ##                    print count_pocket+count_protein+count_distant, count_pocket, count_protein, count_distant
                        count_pocket = 0 ## red, exclude
                        print x,y,z

    ##                count_pocket_min = 25
    ##                count_protein_min = 70
    ##                print count_protein, count_pocket
    ##                if count_pocket >= count_pocket_min or count_protein >= count_protein_min:
    ##                    ## not vicinal to many pocket grid points, but very buried in protein
    ##                    if count_pocket < count_pocket_min and count_protein >= count_protein_min:
    ##                        count_pocket = 27 ## blue, include
    ##                    ## near pocket grid points, but surface of protein
    ##                    elif count_pocket >= count_pocket_min and count_protein <= count_protein_min/2.:
    ##                        count_pocket = 6 ## orange, exclude
    ##                    ## near pocket and semi buried in protein (not surface)
    ##                    elif count_pocket >= count_pocket_min and count_protein > count_protein_min/2.:
    ##                        count_pocket = 13 ## green, include
    ##                    else:
    ##                        stopstop
    ##                else:
    ##                    count_pocket = 0 ## red, exclude
                    l_overlaps += [count_pocket/27.]

    ##            d_class_grid[x][y][z] = count_pocket

    ##                if count > 0:
    ##                    print x,y,z,d_class_grid[x][y][z]
    ##                    stop
    ##    print d_class_grid[-10][-1]
    ##    stop

        print 'first coord', l_coords_probe[0]
        print 'last coord', l_coords_probe[-1]
        print 'coords', len(l_coords_probe)
        print 'bfactors', len(l_overlaps)
        if len(l_coords_probe) != len(l_overlaps):
            print len(l_coords_probe)
            print len(l_overlaps)
            stop

    return l_coords_probe


def write_pdb(l_overlaps,l_coords_probe,pdb,chain,suffix='',):

    fd = open('/media/WDMyBook1TB/2TB/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
    lines_in = fd.readlines()
    fd.close()

    lines_out = []
    l_modres = []

    for line in lines_in:
        record = line[:6].strip()
        if record in ['HETATM','ATOM',]:
            res_name = line[17:20].strip()
            altloc = line[16]
            chain_atom = line[21]
            res_no = int(line[22:26])
            iCode = line[26]
            if res_name == 'HOH':
                continue
            if chain_atom != chain:
                continue
            if record == 'HETATM' and '%1s%4i%1s' %(chain_atom,res_no,iCode,) not in l_modres:
                print line
                continue
            if altloc not in [' ','A','1',]:
                print line
                continue
            atom_no = int(line[6:11])
            bfactor = 100.
            bfactor = 0.
            line = '%60s%6.2f%s' %(line[:60],bfactor,line[66:],)

        elif record == 'MODRES':
            modres_chain = line[16]
            res_no = int(line[18:22])
            iCode = line[22]
            if modres_chain != chain:
                continue
            l_modres += ['%1s%4i%1s' %(modres_chain,res_no,iCode,)]

        lines_out += [line]

    lines_out_probe = []
    atom_no = 0
    for i in range(len(l_overlaps)):
        atom_no_probe = atom_no+i+1
        res_no_probe = atom_no+i+1
        bfactor = 100-100*float(l_overlaps[i])
        coord = l_coords_probe[i]
        x = coord[0]
        y = coord[1]
        z = coord[2]
        lines_out_probe += [
            'HETATM%5i   H  EXT O%4i    %8.3f%8.3f%8.3f  1.00%6.2f           H  \n' %(
                atom_no_probe,res_no_probe,x,y,z,bfactor,
                )
            ]

    lines_out += lines_out_probe

    fd = open('%s_%s_probe%s.pdb' %(pdb,chain,suffix,),'w')
    fd.writelines(lines_out)
    fd.close()

    return


def calc_dotproduct(v1,v2,):

    ## if numerator negative then angle > 90' !!! simplify!!!

    if len(v1) != len(v2):
        print len(v1), len(v2)
        stop

    numerator = 0
    for i in range(len(v1)):
        numerator += v1[i]*v2[i]
    dotproduct = numerator

##    denominator1 = 0
##    denominator2 = 0
##    for i in range(len(v1)):
##        denominator1 += v1[i]*v1[i]
##        denominator2 += v2[i]*v2[i]
##    denominator = math.sqrt(denominator1*denominator2)
##
##    if denominator == 0: ## temp!!!
##        cosang = 0.
##        stop
##    else:
##        cosang = numerator / denominator

    return dotproduct


if __name__ == '__main__':

    ## maximum distance from selected atoms (Angstrom)
    ##dist_max = 2
    dist_max = 6 ## maximum distance from selected atoms (e.g. alpha, heavy backbone) in protein (should be equal to or less than spring cutoff if not sigmoid...)
    dist_max = float(sys.argv[sys.argv.index('--dist_max')+1])
    dist_max_sq = dist_max**2
    ## minimum distance from any atom (Angstrom)
    ##dist_min = 1.5
    dist_min = 3 ## minimum distance from any atom in protein (3 fewer probes than 4; 4 more likely in centroid of large pockets than 3)
    dist_min = float(sys.argv[sys.argv.index('--dist_min')+1])

    ##
    ## input
    ##
    bool_error = False
    if not '--pdb' in sys.argv:
        print 'specify input pdb file with --pdb'
        bool_error = True
    if not '--chain' in sys.argv:
        print 'specify input chain with --chain'
        bool_error = True
    if bool_error == True:
        stop_input_missing
    pdb = sys.argv[sys.argv.index('--pdb')+1]
    chain = sys.argv[sys.argv.index('--chain')+1]
    ##pdb_holo = sys.argv[sys.argv.index('--pdb_holo')+1]
    ##chain_holo = sys.argv[sys.argv.index('--chain_holo')+1]
    ##ligand_ID = sys.argv[sys.argv.index('--ligand_ID')+1]

    main(pdb, chain, dist_max, dist_min,)

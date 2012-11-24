## Fpocket (/usr/local/bin/fpocket)
## Concavity (/usr/local/bin/concavity)
## POCASA (/usr/local/bin/POCASA)
## "null test" - center of protein as measured by average alpha carbon position

## Ligsitecsc (can't compile BALL)
## PASS (doesn't work)
## SURFNET (can't compile)
## Q-SiteFinder (only web)
## Pocket-Finder (only web)
## LIGSITE (abandoned)

## GoodVibes

import numpy, math
import sys
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('/home/tc/svn/GoodVibes/')
import NMA
sys.path.append('/home/tc/svn/tc_sandbox/math')
import statistics

path_pdb = '/media/WDMyBook1TB/2TB'

def main():

    l_pdbs = [
        ## hinge
        '1aj0A',
        '135lA','1og1A','1a2tA','1bqcA','1ra2A',
        ## potential hinge
        '1esoA', '1ak0A','1tmlA',
        '1l9xA','1pjaA','1n29A','1fobA','1cdeA','1ga8A',
        '1bolA','1smlA','1akoA','1lbaA', ## highest scoring grid point near flex residue
        ## multimer, but monomer hinge
        '1jh6A','1vzxA','2rnfA',

        ## not hinge
        '2lipA','1xqwA','1bu7A','1bzcA','1i78A','1qe3A','1d3gA','1rddA',
        '2aceA','1gojA','2cpoA',
        ## probably not hinge
        '1p1xA','1ptdA','1j00A','1u5uA',
        ## multimer and no hinges
        '1a8qA','1a4lA',

        ## unknown
        '1ru4A',
        '1ljlA',
        '1mrqA',
        '1w2nA',
        '1cz1A',
        '1mugA', '1o8aA', '1rtuA', '1uchA', '1b6gA', '1chdA', 
        '1oxaA', '1r44A', '1w0hA', '1expA', '1bt1A', '1i9aA', '1im5A', '1pp4A',
         '2ebnA', '1d1qA', '1ehyA', '1geqA', '1ca2A', '1gnsA', '1eh5A',
        '1l6pA', '1r4zA', '1a2tA', '1di1A', '2ayhA', '1astA',
##        '1cm0A',
        '1cv2A',
         '1dveA',   '1un1A', '1btlA', '2pecA', '1fcqA',
        '1czfA', '1thgA', '1booA', '1iu4A',  '206lA',  '1snzA',
        '1gq8A', '1aqlA', '1ps1A', '1s95A', '1pylA',  '1b6bA', '1pntA',
        '1e1aA', '2f9rA', '1v04A', '2nlrA', '1pbgA', '5cpaA', '1agmA',
        '1byaA',
##        '1vidA', '1h4gA', '1r76A',
##        '1akdA', '1fy2A', '1xqdA',
##        '1d6oA', '1qv0A', '1qjeA', '1fvaA', '1bp2A', '1ah7A', '2pthA', '2engA',
        '2acyA', '1qazA', '2a0nA', '1dl2A', '1gp5A', '1onrA', '1cwyA', '1pudA',
        '1bs9A', '1dinA', '1xyzA', '1bwlA', '1eugA', '1idjA', '1g24A', '1oygA',
        '1hzfA', '9papA', '1eb6A', '1ghsA', '1rbnA', '1bixA', '1bs4A', '1celA',
        '1hkaA', '1b02A', '1qibA', '1u3fA', '1agyA', '1zioA', '1pa9A', '2tpsA',
        '2plcA', '1qk2A', '1j53A', '1m21A',
        ]

    fd = open('datasets/CathDomainList','r')
    lines = fd.readlines()
    fd.close()
    d_CATH = {}
    for line in lines:
        if line[0] == '#':
            continue
        chainID = line[:5]
        if chainID not in l_pdbs:
            continue
        d_CATH[chainID] = '.'.join(line.split()[1:5])

    fd = open('datasets/CSA_2_2_12.dat','r')
    lines = fd.readlines()
    fd.close()

    d_catres = {}
    for i_line in range(1,len(lines)):
        line = lines[i_line]
##        if line[:4] != line[-5:].strip():
##            continue

        if not line[:4]+line[11] in l_pdbs:
            continue

        l = line.strip().split(',')

        pdbID = l[0]
    ##    site_number = l[1]
        chainID = l[3]
        res_name = l[2]
        res_no = int(l[4])
        ref = l[7]
        ## ion
        if chainID == '':
            continue
        if res_name == '':
            continue

        if not pdbID+chainID in d_catres.keys():
            d_catres[pdbID+chainID] = []
        d_catres[pdbID+chainID] += ['%3s%4i' %(res_name,res_no,)]


    l_best = []
    l_worst = []
    l_cutoff = []
    l_methods = [
        'POCASA',
        'Q-SiteFinder',
        'Pocket-Finder',
        'LIGSITE',
        'GoodVibes',
        'ConCavity',
##        'POCASA_w_GV_mode7',
##        'POCASA_w_GV_modes_combined',
        ]
    d_methods = {}
    for method in l_methods:
        d_methods[method] = 0
    s_table = '## minimum distance from all atoms of catalytics site residues to centre of calculated pocket\n'
    for method in l_methods:
        s_table += '\t%s' %(method)
    s_table += '\tbest\tworst\tCATH\tHETATMs\tcatalytic residues'
    s_table += '\n'
    l_table_header = [s_table]
    l_table_body = []
    set_correct_pdbs = set()
    for pdb in l_pdbs:

        s_table = '%s' %(pdb)
        ## PocketPicker - http://gecco.org.chemie.uni-frankfurt.de/pocketpicker/index.html
        ## CASTp - http://sts.bioengr.uic.edu/castp/calculation.php
        ## SURFNET
        ## fpocket
        coords_catalytic, l_ligands = parse_coords_CSA(pdb,d_catres,)

        ##
        ## parse predictions
        ##
        l_coords = []
        try:
            stop
        except:
            if 'POCASA' in l_methods:
                coord_POCASA = parse_POCASA(pdb)
                l_coords += [coord_POCASA]
            if 'Q-SiteFinder' in l_methods:
                coord_QSF = parse_QSiteFinder(pdb)
                l_coords += [coord_QSF]
            if 'Pocket-Finder' in l_methods:
                coord_PF = parse_PocketFinder(pdb)
                l_coords += [coord_PF]
            if 'LIGSITE' in l_methods:
                coord_LS = parse_LIGSITE(pdb)
                l_coords += [coord_LS]
            if 'GoodVibes' in l_methods:
##                path = 'GoodVibes_CSA_dataset'
                path = 'GoodVibes/distmax6_distmin3'
                coord_GV = parse_GoodVibes(pdb,path,)
##                coord_GV = parse_GoodVibes_exclude_flexible(pdb,path) ## tmp!!!
                l_coords += [coord_GV]
            if 'POCASA_w_GV_mode7' in l_methods:
                path = 'POCASA_w_GV_mode7'
                coord_POCASA_w_GV = parse_GoodVibes(pdb,path,)
                l_coords += [coord_POCASA_w_GV]
            if 'POCASA_w_GV_modes_combined' in l_methods:
                path = 'POCASA_w_GV_modes_combined'
                coord_POCASA_w_GV = parse_GoodVibes(pdb,path,)
                l_coords += [coord_POCASA_w_GV]
            if 'ConCavity' in l_methods:
                coord_CC = parse_ConCavity(pdb[:4])
                l_coords += [coord_CC]
        try: None
        except:
            print pdb
            continue

        print pdb

        l_dists_sq_min = []
        for i_coord in range(len(l_coords)):
            coord = l_coords[i_coord]

            dist_sq_min = 1000.
            for coord_cat in coords_catalytic:
                dist_sq = sum((coord_cat-coord)**2)
                if dist_sq < dist_sq_min:
                    dist_sq_min = dist_sq
            l_dists_sq_min += [dist_sq_min]

            ## stats
            if dist_sq_min < 36:
                l_cutoff += [i_coord]
                if l_methods[i_coord] == 'GoodVibes':
                    set_correct_pdbs |= set([pdb])
            if dist_sq_min < 64:
                d_methods[l_methods[i_coord]] += 1

            print i_coord, l_methods[i_coord], round(math.sqrt(dist_sq_min),1)
            s_table += '\t%4.1f' %(round(math.sqrt(dist_sq_min),1))
        i_min = l_dists_sq_min.index(min(l_dists_sq_min))
        i_max = l_dists_sq_min.index(max(l_dists_sq_min))
        print pdb, round(math.sqrt(l_dists_sq_min[i_min]),0), i_min,
        print l_methods[i_min]
        if i_min == l_methods.index('GoodVibes'):
            second_best = round(math.sqrt(min(l_dists_sq_min[:i_min]+l_dists_sq_min[i_min+1:])),0)
            if second_best > 6:
                print round(math.sqrt(l_dists_sq_min[i_min]),0) < 4
                print '***', pdb
                if pdb not in [
                    ## GoodVibes much better
                    '1i78A',
                    ## GoodVibes slightly better
                    '1btlA',
                    ## all wrong
                    '1rddA','1ru4A','1iu4A',
                    ]:
                    stop
        if pdb in d_CATH.keys():
            CATH = d_CATH[pdb]
        else:
            CATH = 'not in CATH'
        s_table += '\t%-13s\t%-13s\t%s\t%s\t%s\n' %(
            l_methods[i_min],
            l_methods[i_max],
            CATH,
            ','.join(list(set(l_ligands))),
            str(','.join(list(set(d_catres[pdb])))),
            )
        l_table_body += [[CATH,s_table,],]

        ## stats
        l_best += [i_min]
        l_worst += [i_max]

    print 'set_correct_pdbs', set_correct_pdbs
    set_correct_pdbs3 = set(['1n29A', '1aj0A', '1og1A', '1a2tA', '1cdeA', '1esoA', '1ak0A', '1fobA', '1ga8A', '1pjaA', '135lA'])
    set_correct_pdbs2 = set(['1esoA', '1ak0A', '1pjaA', '1n29A', '1ga8A', '1ra2A', '1og1A', '1a2tA', '135lA', '1cdeA'])
    set_correct_pdbs1 = set(['1n29A', '1og1A', '1bqcA', '1a2tA', '1cdeA', '1esoA', '1ak0A', '1ga8A', '1ra2A', '1pjaA', '135lA'])
    set_correct_pdbs0 = set(['1n29A', '1tmlA', '1og1A', '1bqcA', '1a2tA', '1fobA', '1cdeA', '1esoA', '1ak0A', '1aj0A', '1ga8A', '1ra2A', '1pjaA', '135lA'])
    print set_correct_pdbs3-set_correct_pdbs0
    print set_correct_pdbs2-set_correct_pdbs0
    print set_correct_pdbs1-set_correct_pdbs0

    for k,v in d_methods.items():
        print 'cutoff (8A)', k,v

    print
    print 'cutoff (6A)'
    for i in range(len(l_methods)):
        print 'cutoff', l_methods[i], l_cutoff.count(i), len(l_best)
    print
    print 'best'
    for i in range(len(l_methods)):
        print 'best', l_methods[i], l_best.count(i), len(l_best)
    print
    print 'worst'
    for i in range(len(l_methods)):
        print 'worst', l_methods[i], l_worst.count(i), len(l_best)

    l_table_body.sort()
    l_table = l_table_header
    for l in l_table_body:
        l_table += [l[1]]

    fd = open('table.tsv','w')
    fd.writelines(l_table)
    fd.close()

    return


def parse_ConCavity(pdb):

##    fd = open('/home/tc/Downloads/1ak1_cc-ligsite_search_blur_pocket.pdb','r')
##    lines = fd.readlines()
##    fd.close()
##    l_coords = []
##    for line in lines:
##        x = float(line[30:38])
##        y = float(line[38:46])
##        z = float(line[46:54])
##        coord = numpy.array([x,y,z,])
##        l_coords += [coord]
##    position = sum(l_coords)/len(l_coords)
##    print position

    fd = open('output/Concavity/%s_concavity_pocket.pdb' %(pdb),'r')
    lines = fd.readlines()
    fd.close()
    l_coords = []
    for line in lines:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        l_coords += [coord]
    position = sum(l_coords)/len(l_coords)

    return position


def parse_GoodVibes_exclude_flexible(pdb,path,):

    ##
    ## calculate amplitudes
    ##
    d_mmCIF = parse_mmCIF.main(pdb[:4],)
    d_coords, l_coords_alpha = mmCIF2coords.main(pdb[:4],d_mmCIF,query_chain=pdb[-1])
    print len(l_coords_alpha)
    ##
    ## eigenvector
    ##
    cutoff = 10
    matrix_hessian = NMA.hessian_calculation(l_coords_alpha,cutoff,)
    eigenvectors, eigenvalues = NMA.diagonalize_hessian(matrix_hessian)
    l_amplitudes = [
        math.sqrt(
            eigenvectors[6][i]**2+eigenvectors[6][i+1]**2+eigenvectors[6][i+2]**2
            )
        for i in range(0,len(eigenvectors[6]),3)
        ]

##    ## write pdb (color by bfactor)
##    l_bfactors = [100*(l_amplitudes[i]-min(l_amplitudes))/(max(l_amplitudes)-min(l_amplitudes)) for i in range(len(l_amplitudes))]
##    fd = open('output/%s/%s_%s_probe.pdb' %(path,pdb[:4],pdb[-1],),'r')
##    lines = fd.readlines()
##    fd.close()
##    index = [-1,None,]
##    lines_out = []
##    for line in lines:
##        record = line[:6].strip()
##        if record != 'ATOM':
##            lines_out += [line]
##        else:
##            res_no = int(line[22:26])
##            if res_no != index[1]:
##                index = [index[0]+1,res_no,]
##                bfactor = l_bfactors[index[0]]
##            line_out = '%s%6.2f%s' %(line[:60],bfactor,line[66:],)
##            lines_out += [line_out]
##    fd = open('output/%s/%s_%s_probe_color_by_amplitude.pdb' %(path,pdb[:4],pdb[-1],),'w')
##    fd.writelines(lines_out)
##    fd.close()

    ## average amplitude
    average = sum(l_amplitudes)/len(l_amplitudes)
    average,stddev = statistics.do_stddev(l_amplitudes)
    ##
    l_coords_rigid = []
    for i in range(len(l_coords_alpha)):
        if l_amplitudes[i] < average:
            l_coords_rigid += [l_coords_alpha[i]]
    l_coords_flexible = []
    for i in range(len(l_coords_alpha)):
        if l_amplitudes[i] > average+0.5*stddev:
            l_coords_flexible += [l_coords_alpha[i]]

    ## parse output
    fd = open('output/%s/%s_%s_probe.pdb' %(path,pdb[:4],pdb[-1],),'r')
    lines = fd.readlines()
    fd.close()

    max_bfactor = None
    coord = None
    for line in lines:
        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue
        res_name = line[17:20]
        if res_name != 'EXT':
            continue

        bfactor = float(line[60:66])

        if bfactor > max_bfactor:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

##            coord_tmp = numpy.array([x,y,z,])

##            bool_vicinal_to_rigid = False
##            for coord_rigid in l_coords_rigid:
##                dist_from_rigid = math.sqrt(sum((coord_rigid-coord_tmp)**2))
##                if dist_from_rigid < 6:
##                    bool_vicinal_to_rigid = True
##                    break
##            if bool_vicinal_to_rigid == False:
##                continue

##            bool_vicinal_to_flexible = False
##            for coord_flexible in l_coords_flexible:
##                dist_from_flexible = math.sqrt(sum((coord_flexible-coord_tmp)**2))
##                if dist_from_flexible < 6:
##                    bool_vicinal_to_flexible = True
##                    break
##            if bool_vicinal_to_flexible == True:
##                continue

##            min_dist = [1000.,None,]
##            for i_coord_alpha in range(len(l_coords_alpha)):
##                coord_alpha = l_coords_alpha[i_coord_alpha]
##                dist_from_alpha = math.sqrt(sum((coord_alpha-coord_tmp)**2))
##                if dist_from_alpha < min_dist[0]:
##                    min_dist = [dist_from_alpha,i_coord_alpha,]
##            if l_amplitudes[min_dist[1]] > average+stddev:
##                continue

            coord = numpy.array([x,y,z,])
            max_bfactor = bfactor

    return coord


def parse_GoodVibes(pdb,path,):

    fd = open('output/%s/%s_%s_probe.pdb' %(path,pdb[:4],pdb[-1],),'r')
    lines = fd.readlines()
    fd.close()

    max_bfactor = None
    coord = None
    for line in lines:
        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue
        res_name = line[17:20]
        if res_name != 'EXT':
            continue
        bfactor = float(line[60:66])
        if bfactor > max_bfactor:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            max_bfactor = bfactor
            coord = numpy.array([x,y,z,])

    return coord


def parse_coords_CSA(pdb,d_catres,):

    l_ligands = []
    l_coords = []
##    ## tmp!!!
##    import os, urllib2
##    if not os.path.isfile('%s/pdb/%s/pdb%s.ent' %(path_pdb,pdb[1:3],pdb[:4],)):
##        if not os.path.isdir('%s/pdb/%s' %(path_pdb,pdb[1:3])):
##            os.mkdir('%s/pdb/%s' %(path_pdb,pdb[1:3]))
##        url = 'http://www.pdb.org/pdb/files/%s.pdb' %(pdb[:4])
##        lines = urllib2.urlopen(url)
##        lines = lines.readlines()
##        fd = open('%s/pdb/%s/pdb%s.ent' %(path_pdb,pdb[1:3],pdb[:4],),'w')
##        fd.writelines(lines)
##        fd.close()
    fd = open('%s/pdb/%s/pdb%s.ent' %(path_pdb,pdb[1:3],pdb[:4],),'r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        record = line[:6].strip()
        if not record in ['ATOM','HETATM',]:
            continue
        res_name = line[17:20].strip()
        if res_name == 'HOH':
            continue
        if record == 'HETATM':
            if res_name not in [
                'BME','GOL','ACT','MSE','EDO','PCA',
                'SO4','PO4',
                'MG','ZN','CA','CL','AU','MN','CU','K','FE','IOD','FE2',
                ]:
                l_ligands += [res_name]
        chain = line[21]
        if chain != pdb[-1]:
            continue
        res_name = line[17:20]
        res_no = line[22:26]
        residue = res_name+res_no
        if residue not in d_catres[pdb]:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        l_coords += [coord]

    return l_coords, l_ligands


def parse_LIGSITE(pdb):

    fd = open('output/LIGSITE/%s_pocket.pdb' %(pdb[:4]),'r')
    lines = fd.readlines()
    fd.close()
    s = lines[0]
    x = float(s[30:38])
    y = float(s[38:46])
    z = float(s[46:54])
    coord = numpy.array([x,y,z,])

    return coord


def parse_PocketFinder(pdb):

    fd = open('output/Pocket-Finder/%s.pdb' %(pdb[:4]),'r')
    lines = fd.readlines()
    fd.close()
    x = []
    y = []
    z = []
    for line in lines:
        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue
##        if line[21] != pdb[-1] and line[:4] == 'ATOM' and line[17:20] in [
##            'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',
##            'MSE',
##            ]:
##            print pdb
##            print line
##            stop
        if line[17:20] != 'YAA':
            continue
        x += [float(line[30:38])]
        y += [float(line[38:46])]
        z += [float(line[46:54])]

    x = sum(x)/len(x)
    y = sum(y)/len(y)
    z = sum(z)/len(z)
    coord = numpy.array([x,y,z,])

    return coord


def parse_QSiteFinder(pdb):

    fd = open('output/Q-SiteFinder/%s.pdb' %(pdb[:4]),'r')
    lines = fd.readlines()
    fd.close()
    x = []
    y = []
    z = []
    for line in lines:
        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue
##        if line[21] != pdb[-1] and line[:4] == 'ATOM' and line[17:20] in [
##            'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',
##            'MSE',
##            ]:
##            print pdb
##            print line
##            stop
        if line[17:20] != 'YAA':
            continue
        x += [float(line[30:38])]
        y += [float(line[38:46])]
        z += [float(line[46:54])]

    x = sum(x)/len(x)
    y = sum(y)/len(y)
    z = sum(z)/len(z)
    coord = numpy.array([x,y,z,])

    return coord


def parse_POCASA(pdb):

    fd = open('output/POCASA/%s_Pocket_DepthCenters.pdb' %(pdb[:4]),'r')
    lines = fd.readlines()
    fd.close()
    s = lines[0]
    l = s.split()
    x = float(l[6])
    y = float(l[7])
    z = float(l[8])
    coord = numpy.array([x,y,z,])

    return coord


if __name__ == '__main__':
    main()

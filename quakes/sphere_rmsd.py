## sphere rmsds for all single point mutations in quakes dataset
## creates 2d contour plot
## all quakes single point mutations

## i cant remember what this script does...

import numpy, math

lines = []
for s in '01234567890abcdefghijklmnopqrstuvwxyz':
    print s
    fd = open('single_point_mutations/%s.txt' %(s),'r')
    lines += fd.readlines()
    fd.close()

##l_HEWL_wt = ['1dpx', '1qio', '1lpi', '2w1y', '1dpw', '1yl0', '7lyz', '1z55', '2w1l', '2w1m', '1lcn', '1wtn', '1azf', '1hf4', '2w1x', '2d4j', '2d4k', '2d4i', '2d91', '1bvx', '1lyo', '2fbb', '1bgi', '1lyz', '2vb1', '194l', '1ps5', '1f10', '1hc0', '3exd', '2cgi', '3lyz', '193l', '3lyo', '3lym', '1sf6', '8lyz', '2z19', '1lzt', '2z12', '1lza', '1lzc', '1lzb', '1lkr', '1lks', '1ved', '2g4q', '2g4p', '1dkj', '1dkk', '1wtm', '1lsf', '2f2n', '1lsd', '1lse', '1lsb', '1lsc', '1lsa', '1aki', '3lzt', '1lys', '2zq3', '1lz9', '1lz8', '2blx', '2bly', '2lyo', '1uig', '1jj3', '1hsw', '1jj1', '1hsx', '1rfp', '1f0w', '1xei', '1xej', '1xek', '5lym', '3lyt', '4lzt', '5lyt', '1c10', '2epe', '5lyz', '1ja2', '1rcm', '1ja6', '1ja4', '1lj3', '1lj4', '6lyt', '4lyz', '4lyt', '1hew', '4lyo', '1vdt', '4lym', '1hel', '1lma', '1vds', '2yvb', '1lje', '1ljg', '1ljf', '1lji', '1ljh', '1ljk', '1ljj', '1jis', '1jit', '1jiy', '1bwh', '1bwi', '1bwj', '1ykx', '6lyz', '2a6u', '2lym', '1b2k', '1bhz', '2aub', '1uco', '1vat', '1vau', '2hso', '1v7t', '1v7s', '1iee', '2c8p', '2cds', '2lzt', '2a7f', '2a7d', '2c8o', '2lyz', '2hs9', '2hs7', '2zq4', '1vdq', '1jpo', '1vdp']
##l_HEWL_mutant = ['1flq', '1flu', '1flw', '1fly', '1fn5', '1hem', '1heo', '1her', '1ios', '1iot', '1ir7', '1ir8', '1ir9', '1kxw', '1kxx', '1kxy', '1lzd', '1uic', '1uid', '1uie', '1uif']

l_gnuplot = []

for i in range(len(lines)):
##    if i < 15700:
##        continue
    if i % 100 == 0:
        print 'mutation', i
    l_gnuplot_pdb = []

    l = lines[i].split()

    try:
        n_chains = int(l[26])
    except:
        print lines[i]
        continue
    bool_spacegroups_identical = l[27]
    if n_chains > 1:
        continue
    pdb1 = l[0]
    pdb2 = l[1]

##    ## HEWL structures only
##    if len(set([pdb1,pdb2,]) & set(l_HEWL_mutant)) == 0:
##        continue

    bm1 = int(l[2])
    bm2 = int(l[3])
    chain1 = l[4]
    chain2 = l[5]
    res_no1 = int(l[6])
    res_no2 = int(l[7])

    pdb = 'pdb/%s/%s%02i%s%02i.pdb' %(pdb1[1],pdb1,bm1,pdb2,bm2,)

    fd = open(pdb,'r')
    lines_pdb = fd.readlines()
    fd.close()

    coord_ref = None
    line_ref = None

    for line in lines_pdb:
        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue
        if line[21] == chain1 and int(line[22:26]) == res_no1 and line[12:16].strip() == 'CA':
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord_ref = numpy.array([x,y,z,])
            if line[26] != ' ':
                print pdb
                print line
                print lines[i]
                stop_run_all_mutants_and_add_biomolecule_and_icode
            line_ref = line
            break

    ## coordinate missing in pdb
    if coord_ref == None:
        continue
            
    for line in lines_pdb:
        record = line[:6].strip()
        if record not in ['ATOM','HETATM',]:
            continue
        if line[21] != chain1:
            print pdb1,pdb2,bm1,bm2,'multiple chains *OR* biomolecule should be added to txt file!'
            stop
            l_gnuplot_pdb = []
            break
            
        if line[12:16].strip() == 'CA':
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            dist = math.sqrt(sum((coord_ref-coord)**2))
            rmsd = occupancy = float(line[54:60])
            relative_rmsd = occupancy = float(line[60:66])
            if rmsd > 60:
                print pdb1,pdb2,bm1,bm2,res_no1,'rmsd',rmsd,pdb
                print line
                print line_ref
                print 'rmsd', rmsd
                stop_rmsd
            if dist > 90: ## 1qwn,1qwu
                print pdb1,pdb2,bm1,bm2,'mut res', res_no1,'rmsd',rmsd,pdb
                print line
                print line_ref
                print 'dist', dist
                stop_dist
            l_gnuplot_pdb += ['%f %f\n' %(dist,rmsd,)]

    l_gnuplot += l_gnuplot_pdb

fd = open('sphere.gnu','w')
fd.writelines(l_gnuplot)
fd.close()


print 'reading lines'
fd = open('sphere.gnu','r')
lines = fd.readlines()
fd.close()
d = {}
rmsd_max = 0
bin_size = 1.
for line in lines:
    l = line.split()
    rmsd_rel = float(l[1])
    dist = float(l[0])
    rmsd_disc = rmsd_rel-rmsd_rel%bin_size
    if rmsd_disc > rmsd_max:
        rmsd_max = rmsd_disc
    dist_disc = dist-dist%bin_size
    if not dist_disc in d.keys():
        d[dist_disc] = {}
    if not rmsd_disc in d[dist_disc].keys():
        d[dist_disc][rmsd_disc] = 0
    d[dist_disc][rmsd_disc] += 1
print 'lines read'
l_out = []
for dist_discrete in range(0,int(max(d.keys())/bin_size)+1,):
    dist_discrete *= bin_size
    for rmsd_discrete in range(0,int(rmsd_max/bin_size)+1,):
        rmsd_discrete *= bin_size
        if not dist_discrete in d.keys():
            count = 0
        elif not rmsd_discrete in d[dist_discrete].keys():
            count = 0
        else:
            count = d[dist_discrete][rmsd_discrete]
        l_out += ['%s %s %s\n' %(dist_discrete,rmsd_discrete,count,)]
    l_out += ['\n']
fd = open('discrete.gnu','w')
fd.writelines(l_out)
fd.close()


##    def spherermsd(
##        self,
##        pdb1, pdb2, ## pdbs
##        d_header, ## sequences
##        d_coordinates, ## coordinates
##        l_equivalent_chains, ## equivalent chains
##        d_chains_interpdb_sequence_similar, ## mutations
##        tv1, rm, tv2, ## transformation
##        ):
##
##        rewritethisfunction
##
##        instance_geometry = geometry.geometry()
##
##        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
##            rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']
##
##            ATOMseq1,d_res_nos1 = self.ATOM2seq(d_coordinates[pdb1], rep_chain1, d_header[pdb1])
##            ATOMseq2,d_res_nos2 = self.ATOM2seq(d_coordinates[pdb2], rep_chain2, d_header[pdb2])
##            l1 = d_chains_interpdb_sequence_similar[rep_chain1]['l1']
##            l2 = d_chains_interpdb_sequence_similar[rep_chain1]['l2']
##
##            if l1 > 0 or l2 > 0: ## e.g. 1ftg.pdb,1dx9.pdb
##                print ATOMseq1
##                print ATOMseq2
##                s1 = d_chains_interpdb_sequence_similar[rep_chain1]['s1']
##                s2 = d_chains_interpdb_sequence_similar[rep_chain1]['s2']
##                print s1
##                print s2
##                print pdb1, pdb2
##                expected
##
##            (
##                coordinates1, coordinates2, rescount,
##                ) = self.xxxATOMxxxrecordsxxx2xxxcoordinates(
##                    d_coordinates, pdb1, pdb2, rep_chain1, rep_chain2, d_res_nos1, d_res_nos2,
##                    l1, l2, len(ATOMseq1), d_ATOMseq, tv1=tv1, rm=rm, tv2=tv2
##                    )
##
##            l_mutations = d_chains_interpdb_sequence_similar[rep_chain1]['l_mutations']
##            for mutation in l_mutations:
##                res_no1 = d_res_nos1[mutation[0]]['res_no']
##                res_no2 = d_res_nos2[mutation[0]]['res_no']
##                iCode1 = d_res_nos1[mutation[0]]['iCode']
##                iCode2 = d_res_nos2[mutation[0]]['iCode']
##                hypocenter1 = d_coordinates[pdb1]['chains'][rep_chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms']['CA']['coordinate']
##                hypocenter2 = d_coordinates[pdb2]['chains'][rep_chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms']['CA']['coordinate']
##
##                d_rmsds = {4:0,8:0,16:0,32:0}
##                for dist in d_rmsds.keys():
##
##                    sqdist = dist**2
##                    sphere_coordinates1 = []
##                    sphere_coordinates2 = []
##
##                    for i in range(len(coordinates1)):
##
##                        coordinate1 = coordinates1[i]
##                        coordinate2 = coordinates2[i]
##
##                        sqdist1 = sum((hypocenter1-coordinate1)**2)
##                        sqdist2 = sum((hypocenter2-coordinate2)**2)
##
##                        if min(sqdist1,sqdist2) < sqdist:
##
##                            sphere_coordinates1 += [coordinate1]
##                            sphere_coordinates2 += [coordinate2]
##
##                    rmsd = instance_geometry.superpose(
##                        sphere_coordinates1,sphere_coordinates2
##                        )
##
##                    print rep_chain1, mutation, dist, rmsd, float(len(sphere_coordinates1))/len(coordinates1)
##
##                    d_rmsds[dist] = rmsd
##
##        rmsd4 = d_rmsds[4]
##        rmsd8 = d_rmsds[8]
##        rmsd16 = d_rmsds[16]
##        rmsd32 = d_rmsds[32]
##
##        return rmsd4, rmsd8, rmsd16, rmsd32

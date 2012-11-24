## 2e3m ingen noder i ligand binding site i midten af proteinet! tilfoej node i "gennemsnitsdistance" fra bindingsrester.
## 1c48 ligasite numbering refers to B chain
## 1g95/1tyv biounit is not a monomer. exclude multimers?
## 1q7m biounit is a monomer, but 2 chains in ASU. exclude if multiple chains in ASU?
## exclude probe atoms if not in a cavity? how determine if in cavity?

import os, numpy, math

import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import gnuplot

d_catalytic = {

    ## in ligasite
    '4ake':[13,123,156,158,159,167,],
    '2lzm':[11,20,],
##    '1rbx':[12,41,119,120,],
##    '1stn':[35,87,],
    '2nwd':[35,53,],
##    '1cex':[42,120,121,175,188,],
    '1ex6':[34,38,41,50,69,78,79,80,99,100,103,],

    ## not in ligasite
    '2vb1':[35,52,], ## HEWL
##    '1eea':[200,327,440,], ## AChE
    '5dfr':[5,27,28,31,54,94,], ## DHFR
    '1lkz':[28,30,31,81,84,94,95,96,97,103,121,],
    '1dv1':[116,159,163,164,165,166,169,201,202,203,204,236,278,287,288,347,],
##    '1j3h':[166,170,],
    '2p9q':[38,215,396,],
    '2ili':[106,199,],
    '1ejd':[23,115,305,397,],

    }

d_ligasite = {}
fd = open('../ligasite.csv','r')
lines = fd.readlines()
fd.close()
for line in lines[1:]:
    resID = line.split(',')[1]
    res_no = int(resID[3:])
    pdb = line.split(',')[0]
    if not pdb in d_ligasite.keys():
        d_ligasite[pdb] = []
    d_ligasite[pdb] += [res_no]
d_catalytic = d_ligasite
print len(d_ligasite.keys())
stop

##d_csa = {}
##fd = open('CSA_2_2_12.dat','r')
##lines = fd.readlines()
##fd.close()
##for line in lines[1:]:
##    l = line.split(',')
##    if l[0] != l[-1][:4]:
##        continue
##    res_no = int(l[4])
##    pdb = l[0]
##    if not pdb in d_csa.keys():
##        d_csa[pdb] = []
##    d_csa[pdb] += [res_no]
##d_catalytic = d_csa


l_gnuplot = []
l_fn = os.listdir(os.getcwd())
l_fn.sort()
l_vicinal = []
##minmindist = [999.,'',]
for i_fn in range(len(l_fn)):
    fn = l_fn[i_fn]
    pdb = fn[:4]
    if i_fn % 100 == 0:
        print i_fn, len(l_fn), fn
##    ## biounit not monomer
##    if pdb in ['1tyv','1q7m','1g95',]:
##        continue

    if pdb in [
        ## binding property of ligand determined by cys position
        ## covalent binding!!!???
        ## cant use to test our hypothesis that proteins evolve dynamics/binding sites to accomodate ligand
        '2lzm','3lzm','1lyd',
        ## "binding residues" close to binding site of other entity/protein in holo structure
        '1rgp','1aaj',
        ]:
        continue
    if not fn[-4:] == '.pdb':
        continue
    if not pdb in d_catalytic.keys():
        continue

    ##
    ## skip if not monomeric
    ##
    fd = open('/data/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()
    monomeric = False
    for line in lines:
        if '_pdbx_struct_assembly.oligomeric_details' in line:
            if line.split()[-1] == 'monomeric':
                monomeric = True
    if monomeric == False:
        continue
    ##
    ## skip if multiple domains
    ##
    fd = open('../CathDomall','r')
    lines = fd.readlines()
    fd.close()
    n_domains = None
    for line in lines:
        if line[:4] == pdb:
            n_domains = int(line[7:9])
            break
    if n_domains == None:
        print pdb, 'n_domains', n_domains
        continue
    if n_domains > 1:
        continue
    
    l_overlaps = []
    l_dist_min = []
    l_coords_probe = []
    l_tempfactors = []
    l_coords_protein = []
    fd = open(fn,'r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        record = line[:6].strip()
        if record in ['ATOM','HETATM',]:
            res_name = line[17:20]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            res_no = int(line[22:26])
            if res_name == 'EXT':
                l_coords_probe += [coord]
                tempfactor = float(line[60:66])
                l_tempfactors += [tempfactor]
            elif res_no in d_catalytic[fn[:4]]:
                if len(l_coords_protein) == 92:
                    print line
                l_coords_protein += [coord]

##    min_overlap = (100-max(l_tempfactors))/100.

    ## no A chain...
    if len(l_coords_probe) == 0:
        continue

    for i in range(len(l_coords_probe)):
        min_dist = 999.
        tempfactor = float(l_tempfactors[i])
        for j in range(len(l_coords_protein)):
            dist = math.sqrt(sum((l_coords_probe[i]-l_coords_protein[j])**2))
##            if pdb == '4ake' and i == 311:
##                print dist
            if dist < min_dist:
                min_dist = dist
##                if pdb == '1akz' and i in [60,67] and dist < 13.5:
                if pdb == '1akz' and dist == 13.486878252583137:
                    print i, j, dist, l_coords_protein[j]
##        if tempfactor > 5 and min_dist > 10:
##            print tempfactor, min_dist, fn, l_coords_probe[i]
        overlap = (100-tempfactor)/100.
##        if overlap == min_overlap and min_dist < minmindist[0]:
##            minmindist = [min_dist,fn,]

##        l_gnuplot += ['%s %s\n' %(min_dist,overlap,)]

        l_overlaps += [overlap]
        l_dist_min += [min_dist]
    if min_dist == 999:
        print pdb, min_dist
        stop
        continue

    l = []
    print len(l_coords_probe), len(l_overlaps), len(l_dist_min)
    for i in range(len(l_overlaps)):
        l += [[l_overlaps[i],l_dist_min[i],l_coords_probe[i][0],l_coords_probe[i][1],l_coords_probe[i][2],]]
    l.sort()

    ## cluster probes
    d_clusters = {0:{'coords':[numpy.array([l[0][2],l[0][3],l[0][4],])],'indexes':[0],}}
    l_clusters_index = [0]
    for i1 in range(20):
        cluster_growth = 0
        for i2 in range(20):
            if i2 in l_clusters_index:
                continue
            coord2 = numpy.array([l[i2][2],l[i2][3],l[i2][4],])
            for cluster in d_clusters.keys():
                for coord in d_clusters[cluster]['coords']:
                    if (
                        coord[0] in [coord2[0]-2,coord2[0],coord2[0]+2,]
                        and
                        coord[1] in [coord2[1]-2,coord2[1],coord2[1]+2,]
                        and
                        coord[2] in [coord2[2]-2,coord2[2],coord2[2]+2,]
                        ):
                        l_clusters_index += [i2]
                        d_clusters[cluster]['coords'] += [coord2]
                        d_clusters[cluster]['indexes'] += [i2]
                        cluster_growth += 1
                        break
                if cluster_growth == 1:
                    break
            if cluster_growth == 1:
                break
            elif cluster_growth == 0 and i2 not in l_clusters_index:
                d_clusters[max(d_clusters.keys())+1] = {'coords':[coord2],'indexes':[i2],}
                l_clusters_index += [i2]
            else:
                print cluster_growth, i2, l_clusters_index
                stop
            
    bool_vicinal = False
    ## ligasite, 151 of 185 (82%) if top 10
    ## ligasite, 143 of 185 (77%) if top 5
    ## ligasite, 136 of 186 (73%) if top 3
    ## ligasite,  99 of 122 (81%) if top 3 and not multi domain (CATH) proteins
    ## ligasite,  99 of 131 (76%) if top 1 and not multi domain (CATH) proteins
    ## ligasite,  115 of 150 (77%) if top 1 and not multi domain (CATH) proteins and within 5Angstrom
    ## ligasite,  110 of 150 (73%) if top 1 and not multi domain (CATH) proteins and within 4Angstrom
##    for i in range(3):
##        if l[i][1] < 5:
##            bool_vicinal = True
##            break

##    for cluster in range(3):
    for cluster in range(1):
        for i in d_clusters[cluster]['indexes']:
            if l[i][1] < 4: ## 4 Angstrom in Huang2006,Weisel2007,Yu2009
                bool_vicinal = True
                break
        if bool_vicinal == True:
            break
    if bool_vicinal == True:
        l_vicinal += [True]
    else:
        print pdb, 'not vicinal', l[0], l[1], l[2]
        l_vicinal += [False]
        
    for i in range(len(l_overlaps)):
        overlap = l_overlaps[i]
        overlap_normalized = (overlap-min(l_overlaps))/(1.-min(l_overlaps))
        min_dist = l_dist_min[i]
##        if overlap_normalized < 0.5 and min_dist > 10:
##            print pdb, overlap, overlap_normalized, min_dist
        if overlap_normalized < 0.5 and min_dist > 30:
            print pdb, overlap, overlap_normalized, min_dist
        l_gnuplot += ['%s %s %s\n' %(min_dist,overlap_normalized,pdb,)]

print l_vicinal.count(True), len(l_vicinal), l_vicinal.count(True)/float(len(l_vicinal))

##print minmindist

prefix = 'gnuplot'

fd = open('%s.gnuplotdata' %(prefix),'w')
fd.writelines(l_gnuplot)
fd.close()

gnuplot.scatter_plot_2d(
    prefix,
##    bool_regression_linear = True,
    xlabel='minimum distance to catalytic site residue(s)',
    ylabel='overlap between apo and holo eigenvectors',
    bool_remove = False,
    )

print l_vicinal.count(True), len(l_vicinal), l_vicinal.count(True)/float(len(l_vicinal))

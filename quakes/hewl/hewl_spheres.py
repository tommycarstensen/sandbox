## HEWL sphere RMSDs
## 2D scatter plot of RMSD and distance from sphere center (mutation site)

import os, random
import sys

l_radii = [5,10,20,50,]

atoms = 'CA'
atoms = 'heavy'

wt_comparison = 'all'
wt_comparison = 'closest'
wt_comparison = 'starting_model'
wt_comparison = 'cutoff'

if wt_comparison in ['all','closest','starting_model',]:
    l_cutoffs = [9999.9,]
else:
    l_cutoffs = [.2,.4,.8,1.6,3.2,6.4,]

l_keys = ['mutant','wt',]

cluster = 9

fd = open('pdbS95bF.out','r')
lines = fd.readlines()
fd.close()
line = lines[cluster]
l_pdbs = line.split()
for i in range(len(l_pdbs)):
    l_pdbs[i] = l_pdbs[i][:4].lower()

l_wts = [
##    ## hexagonal
##    '2fbb',
##    ## not P 43 21 2
##    '1xei', '1xej', '1xek', '1v7s', '5lym', '3lyt', '4lzt', '1uco', '7lyz', '2zq4', '1wtm', '1wtn', '1lzt', '2z12', '1ps5', '1lcn', '1lj3', '1lkr', '1hf4', '1v7t', '2lzt', '1ved', '1lks', '4lyt', '1lj4', '2d4j', '2d4k', '2d4i', '2f2n', '2z19', '1vdq', '1lma', '1aki', '1jpo', '3lzt', '1lje', '2zq3', '1ljg', '1ljf', '1lji', '1ljh', '1ljk', '1ljj', '1f0w', '1bgi', '2vb1', '1jj3', '1hsw', '1jj1', '1rcm', '1hsx', '1lys', '1f10', '1vdp',

    ## P 43 21 2
    '1bwh_a', '1bwi_a', '1dpx_a', '3exd_a', '2cgi_a', '1qio_a', '6lyz_a', '6lyt_a', '1lpi_a', '2w1y_a', '3lyz_a', '1yl0_x', '193l_a', '5lyt_a', '1ykx_x', '2aub_a', '1c10_a', '2epe_a', '3lym_a', '1vau_a', '2w1l_a', '2w1m_a', '8lyz_a', '1z55_a', '1jit_a', '1azf_a', '2a7d_a', '1lza_a', '1iee_a', '5lyz_a', '2w1x_a', '1bwj_a', '2g4q_a', '2g4p_a', '2blx_a', '4lyz_a', '2cds_a', '1lsf_a', '1dpw_a', '1lsd_a', '1lse_a', '1lsb_a', '1lsc_a', '1lsa_a', '4lyo_a', '1vdt_a', '4lym_a', '1bvx_a', '1vds_a', '2yvb_a', '2c8p_a', '1bhz_a', '2c8o_a', '1jiy_a', '1lyo_a', '1lz9_a', '1lz8_a', '2lym_a', '2bly_a', '1jis_a', '1lyz_a', '194l_a', '2lyz_a', '2a7f_a',
    '1hel_a',
    '1rfp_a',
    '1uig_a',
##    '1b2k_a','1hc0_a','1vat_a','2d91_a', ## iodide ions

##    '2hso', '2hs9', '2hs7', ## powder diffraction
##    '2a6u', ## powder diffraction
##    '1ja2', '1ja6', '1ja4', ## powder diffraction
    ]
d_mutants = {'1uif_a': 14, '1fn5_a': 48, '1kxw_a': 26, '1fly_a': 101, '1ir7_a': 77, '1flu_a': 66, '1ir9_a': 97, '1flw_a': 70, '1uid_a': 14, '1flq_a': 116, '1kxy_a': 17, '1uie_a': 14, '1ir8_a': 57, '1iot_a': 11, '1ios_a': 11, '1hem_a': 90, '1uic_a': 14, '1heo_a': 54, '1her_a': 39, '1lzd_a': 61}
l_mutants = d_mutants.keys()

d_pdbs = {'wt':l_wts,'mutant':l_mutants,}

print 'reading HEWL overall rmsds'
fd = open('%s_rmsds.txt' %(cluster),'r')
s = fd.read()
fd.close()
d_rmsds_overall = eval(s)

def main():

    d_rmsds = {} ## wt/mutant, pdb, r, res_index
    for k in l_keys:
        d_rmsds[k] = {}
        for pdb in d_pdbs[k]:
            d_rmsds[k][pdb] = {}
            for r in l_radii:
                d_rmsds[k][pdb][r] = {}
                for res_index in range(129):
                    d_rmsds[k][pdb][r][res_index] = {}

    d_rmsds = parse_rmsds(d_rmsds,cluster,)

    for s_mutant in l_mutants:
        sample_and_plot_rmsds(d_rmsds, s_mutant,)

    return


def sample_and_plot_rmsds(d_rmsds, s_mutant,):

##    if s_mutant != '1flq':
##        return

    print 'sample and plot'

    prefix = '%i_wtvmut_%s_%s_%s' %(cluster,atoms,wt_comparison,s_mutant,)

    if wt_comparison == 'starting_model':
        if s_mutant not in ['1fly_a','1kxw_a','1kxy_a','1uic_a','1uid_a','1uie_a','1uif_a',]:
            return

    for k in l_keys:
        print k, s_mutant
        l_data = []
        if k == 'wt':
            l_pdbs = l_wts
            l_res_indexes = [d_mutants[s_mutant]]
        elif k == 'mutant':
##            l_pdbs = len(l_wts)*[s_mutant]
            l_pdbs = 1*[s_mutant]
            l_res_indexes = [d_mutants[s_mutant]]
        for pdb in l_pdbs:
##            if k == 'wt' and pdb != '1rfp': ## tmp!!!
##                continue
            for i_rmsd_cutoff in range(len(l_cutoffs)):
                rmsd_cutoff = l_cutoffs[i_rmsd_cutoff]

                if wt_comparison == 'cutoff':
                    ## use all pdbs below rmsd treshold
                    l_pdbs_cutoff = []
                    for pdb_wt in l_wts:
                        if k == 'wt' and pdb == pdb_wt:
                            continue
                        if d_rmsds_overall[pdb][pdb_wt] < rmsd_cutoff:
                            l_pdbs_cutoff += [pdb_wt]
                    n_loop = int(len(l_pdbs_cutoff)/2.)
                    if len(l_pdbs_cutoff) <= 1:
                        continue
                elif wt_comparison == 'closest':
                    ## use only most similar pdb
                    l1 = list(l_wts)
                    if k == 'wt':
                        l1.remove(pdb)
                    l2 = [[d_rmsds_overall[pdb][pdb_wt],pdb_wt,] for pdb_wt in l1]
                    l_pdbs_cutoff = [min(l2)[1]]
                    n_loop = 1
                elif wt_comparison == 'all':
                    l_pdbs_cutoff = list(l_wts)
                    if k == 'wt':
                        l_pdbs_cutoff.remove(pdb)
                    n_loop = 50
                elif wt_comparison == 'starting_model':
                    if k == 'wt' and pdb == '1rfp_a':
                        continue
                    l_pdbs_cutoff = ['1rfp_a']
                    n_loop = 1
                else:
                    stop
              
                for i_radius in range(len(l_radii)):
                    r = l_radii[i_radius]
                    x = i_radius
                    if k == 'mutant':
                        x += .4
                    x += .05*i_rmsd_cutoff
                    for res_index in l_res_indexes:
                        if len(d_rmsds[k][pdb][r][res_index]) < n_loop:
                            print k, pdb, r, res_index, d_rmsds[k][pdb][r][res_index]
                            print n_loop
                            stop

                        ## number of points on plot
##                        for i in range(n_loop):
                        for i in range(len(l_pdbs_cutoff)):

##                            l_rmsds = list(d_rmsds[k][pdb][r][res_index].values())

                            l_rmsds = []
                            for pdb_wt in l_pdbs_cutoff:
                                try:
                                    if d_rmsds_overall[pdb][pdb_wt] < rmsd_cutoff:
                                        l_rmsds += [d_rmsds[k][pdb][r][res_index][pdb_wt]]
                                except:
                                    print k,pdb,r,res_index,pdb_wt
                                    print s_mutant
                                    print pdb_wt
                                    print l_pdbs_cutoff
                                    if d_rmsds_overall[pdb][pdb_wt] < rmsd_cutoff:
                                        l_rmsds += [d_rmsds[k][pdb][r][res_index][pdb_wt]]
                                    stop

                            ##
                            ## calculate average RMSD
                            ## number of RMSDs used for each average/point on the plot
                            ##
                            l_rmsds_n = []
                            ## n_loop wildtype comparisons
                            for j in range(n_loop):
                                index = int(random.random()*len(l_rmsds))
                                l_rmsds_n += [l_rmsds[index]]
                                ## make sure no pdb is sampled more than once
                                del l_rmsds[index]

                            ## plot average
                            y = average = sum(l_rmsds_n)/len(l_rmsds_n)

                            error = 0
    ##                        if average > 2.0:
    ##                            print average, pdb, s_mutant
    ##                            print 'res_no', d_mutants[s_mutant]
    ##                            print l_indexes
    ##                            print 'radius', r
    ##                            print l_rmsds
    ##                            stop
                            l_data += ['%s %s %s %s\n' %(x,y,error,pdb,)]

        fd = open('%s_%s.gnuplot_data' %(prefix,k,),'w')
        fd.writelines(l_data)
        fd.close()

    ylabel = 'average_rmsd'
    title = s_mutant,d_mutants[s_mutant]+1
    gnuplot(l_radii,prefix,ylabel,title,)

    os.remove('%s_wt.gnuplot_data' %(prefix))
    os.remove('%s_mutant.gnuplot_data' %(prefix))

    return


def parse_rmsds(d_rmsds,cluster,):

    print 'parsing rmsds'

    for k in l_keys:

        fd = open('sphere/%i_%s_%s.txt' %(cluster,atoms,k,),'r')
        lines = fd.readlines()
        fd.close()
        print k, 'looping over lines'
        for i_line in range(len(lines)):
            line = lines[i_line]
            l = line.split()
            pdb1 = l[0]
            pdb2 = l[1]
##            ## sugar
##            ## 1sf6,1hew,1lzc,1lzb
##            ## 2lyo,3lyo (acetonitrile cross linker)
##            ## other species
##            ## 1dkj,1dkk
##            ## double mutant
##            ## 1kxx
            if pdb1 in ['1sf6','1hew','1lzc','1lzb','2lyo','3lyo','1dkj','1dkk','1kxx',]:
                continue
            if pdb2 in ['1sf6','1hew','1lzc','1lzb','2lyo','3lyo','1dkj','1dkk','1kxx',]:
                continue
            bm1 = int(l[2])
            bm2 = int(l[3])
            chain1 = l[4]
            chain2 = l[5]
            pdb1 = '%s_%s' %(pdb1,chain1.lower())
            pdb2 = '%s_%s' %(pdb2,chain2.lower())
            if pdb1 not in l_wts+l_mutants:
                continue
            if pdb2 not in l_wts+l_mutants:
                continue
            ## only use biomolecule 1
            if bm1 != 1 or bm2 != 1:
                continue
            res_index = int(l[6])
            bool_author_identical = l[11]
            for i in range(len(l_radii)):
                r = l_radii[i]
                rmsd = float(l[7+i])
                for pdb_a,pdb_b in [[pdb1,pdb2,],[pdb2,pdb1,],]:
##                    if k == 'wt' and pdb_a != '1rfp': ## tmp!!!
##                        continue
##                    if k == 'mutant' and '1uig' not in [pdb1,pdb2,]: ## tmp!!!
##                        continue
                    ## continue if wildtype and wt v mutant
                    if k == 'mutant' and pdb_a not in d_rmsds[k].keys():
                        continue
                    d_rmsds[k][pdb_a][r][res_index][pdb_b] = rmsd

    return d_rmsds


def gnuplot(l_radii,prefix,ylabel,title,):

    ## plot w errorbars u 1:2:3
    l_settings = [
        'set terminal postscript eps enhanced color "Helvetica" 36\n',
        'set output "%s.ps"\n' %(prefix),
##        'set terminal png\n',
##        'set output "gnuplot.png"\n',
        'set size 3,3\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set title "%s"\n' %(title,),
        'set ylabel "%s"\n' %(ylabel),
        'set xlabel "radii=%s, RMSDcutoffs=%s"\n' %(l_radii,[str(round(cutoff,1))[:3] for cutoff in l_cutoffs],),
        'plot [:][0:1.5]"%s_wt.gnuplot_data" u 1:2:3 w errorb lw 2 ps 4 lc 1 t "wt", "%s_mutant.gnuplot_data" u 1:2:3 w errorb lw 2 ps 4 lc 2 t "mutant"\n' %(prefix,prefix,),
        ]

    fd = open('%s.gnuplot_settings' %(prefix),'w')
    fd.writelines(l_settings)
    fd.close()

    os.system('/software/bin/gnuplot %s.gnuplot_settings' %(prefix))

    print 'imagemagick', prefix
##    if not os.path.isfile('sphere_%s_%s_%s_%s.png' %(prefix,s_mutant,d_mutants[s_mutant],)):
    os.system('convert %s.ps sphere_%s.png' %(prefix,prefix,))
    os.system('rm %s.ps' %(prefix))
    os.remove('%s.gnuplot_settings' %(prefix))

##    os.system('mv gnuplot.png sphere_%s_%s_%s_%s.png' %(prefix,s_mutant,d_mutants[s_mutant],))

    return    


if __name__ == '__main__':
    main()

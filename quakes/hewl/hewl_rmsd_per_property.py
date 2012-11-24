## HEWL scatter plot of sampled RMSD for wt and mutant for diff props
## conclusion: more diff between wts than wts and muts

## imports
import os, random, numpy, math
import core
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb')
import whatif_surface

## appendix
d_mutants = {
    '1flq_a': {'mutation':'G117A','color':'FF0000',},
    '1flu_a': {'mutation':'G67A' ,'color':'FF8000',},
    '1flw_a': {'mutation':'G71A' ,'color':'FFFF00',},
    '1fly_a': {'mutation':'G102A','color':'80FF00',},
    '1fn5_a': {'mutation':'G49A' ,'color':'00FF00',},
    '1hem_a': {'mutation':'S91T' ,'color':'00FF80',},
    '1heo_a': {'mutation':'I55V' ,'color':'00FFFF',},
    '1her_a': {'mutation':'T40S' ,'color':'0080FF',},
    '1ios_a': {'mutation':'M12F' ,'color':'0000FF',},
    '1iot_a': {'mutation':'M12L' ,'color':'0000FF',},
    '1ir7_a': {'mutation':'I78M' ,'color':'8000FF',},
    '1ir8_a': {'mutation':'I58M' ,'color':'FF00FF',},
    '1ir9_a': {'mutation':'I98M' ,'color':'FF0080',},
    '1kxw_a': {'mutation':'N27D' ,'color':'DFDFDF',},
    '1kxy_a': {'mutation':'D18N' ,'color':'BFBFBF',},
    '1lzd_a': {'mutation':'W62Y' ,'color':'9F9F9F',},
    '1uic_a': {'mutation':'H15A' ,'color':'7F7F7F',},
    '1uid_a': {'mutation':'H15F' ,'color':'7F7F7F',},
    '1uie_a': {'mutation':'H15G' ,'color':'7F7F7F',},
    '1uif_a': {'mutation':'H15V' ,'color':'7F7F7F',},
    }
l_mutants = d_mutants.keys()
l_mutants.sort()

d_solvacc = {1: 'exposed', 2: 'exposed', 3: 'exposed', 4: 'exposed', 5: 'exposed', 6: 'exposed', 7: 'exposed', 8: 'buried', 9: 'buried', 10: 'exposed', 11: 'exposed', 12: 'exposed', 13: 'exposed', 14: 'exposed', 15: 'exposed', 16: 'exposed', 17: 'buried', 18: 'exposed', 19: 'exposed', 20: 'exposed', 21: 'exposed', 22: 'exposed', 23: 'exposed', 24: 'exposed', 25: 'exposed', 26: 'exposed', 27: 'exposed', 28: 'buried', 29: 'buried', 30: 'buried', 31: 'buried', 32: 'buried', 33: 'exposed', 34: 'exposed', 35: 'exposed', 36: 'exposed', 37: 'exposed', 38: 'exposed', 39: 'exposed', 40: 'exposed', 41: 'exposed', 42: 'exposed', 43: 'exposed', 44: 'exposed', 45: 'exposed', 46: 'exposed', 47: 'exposed', 48: 'exposed', 49: 'exposed', 50: 'buried', 51: 'exposed', 52: 'exposed', 53: 'exposed', 54: 'buried', 55: 'exposed', 56: 'exposed', 57: 'exposed', 58: 'exposed', 59: 'exposed', 60: 'buried', 61: 'exposed', 62: 'exposed', 63: 'exposed', 64: 'buried', 65: 'exposed', 66: 'exposed', 67: 'exposed', 68: 'exposed', 69: 'exposed', 70: 'exposed', 71: 'exposed', 72: 'exposed', 73: 'exposed', 74: 'exposed', 75: 'exposed', 76: 'exposed', 77: 'exposed', 78: 'exposed', 79: 'exposed', 80: 'exposed', 81: 'exposed', 82: 'exposed', 83: 'exposed', 84: 'exposed', 85: 'exposed', 86: 'exposed', 87: 'exposed', 88: 'exposed', 89: 'exposed', 90: 'exposed', 91: 'exposed', 92: 'exposed', 93: 'exposed', 94: 'exposed', 95: 'buried', 96: 'exposed', 97: 'exposed', 98: 'exposed', 99: 'exposed', 100: 'exposed', 101: 'exposed', 102: 'exposed', 103: 'exposed', 104: 'exposed', 105: 'buried', 106: 'exposed', 107: 'exposed', 108: 'exposed', 109: 'exposed', 110: 'exposed', 111: 'exposed', 112: 'exposed', 113: 'exposed', 114: 'exposed', 115: 'buried', 116: 'exposed', 117: 'exposed', 118: 'exposed', 119: 'exposed', 120: 'exposed', 121: 'exposed', 122: 'exposed', 123: 'exposed', 124: 'exposed', 125: 'exposed', 126: 'exposed', 127: 'exposed', 128: 'exposed', 129: 'exposed'}

## alpha helix
helix = range(5,15)+range(25,37)+range(89,100)+range(109,115)
## 3/10-helix
helix += range(80,85)+range(104,108)+range(120,124)
## sheet
sheet = range(43,46)+range(51,54)+range(58,60)

d_secondary = {}
for res_no in range(1,130):
    d_secondary[res_no] = 'other'
for secstruc,l_residues in [
    ['helix',helix,],
    ['sheet',sheet,],
    ]:
    for res_no in l_residues:
        d_secondary[res_no] = secstruc

d_domains = {}
for res_no in range(1,130):
    d_domains[res_no] = 'beta'
for res_no in range(40,87):
    d_domains[res_no] = 'alpha'


d_properties = {
    'secstruc':['helix','sheet','other',],
    'solvacc':['exposed','buried',],
##    'normBfactor':['0-20','20-80','80-100',],
    'domain':['alpha','beta',],
##    'mutdist':[],
    }

d_prop2prop1 = {}
for prop1 in d_properties.keys():
    for prop2 in d_properties[prop1]:
        d_prop2prop1[prop2] = prop1

method = 'alpha'
method = 'heavy'

i_cluster = 5
ref_seq = 'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL'
ref_seq = ['LYS', 'VAL', 'PHE', 'GLY', 'ARG', 'CYS', 'GLU', 'LEU', 'ALA', 'ALA', 'ALA', 'MET', 'LYS', 'ARG', 'HIS', 'GLY', 'LEU', 'ASP', 'ASN', 'TYR', 'ARG', 'GLY', 'TYR', 'SER', 'LEU', 'GLY', 'ASN', 'TRP', 'VAL', 'CYS', 'ALA', 'ALA', 'LYS', 'PHE', 'GLU', 'SER', 'ASN', 'PHE', 'ASN', 'THR', 'GLN', 'ALA', 'THR', 'ASN', 'ARG', 'ASN', 'THR', 'ASP', 'GLY', 'SER', 'THR', 'ASP', 'TYR', 'GLY', 'ILE', 'LEU', 'GLN', 'ILE', 'ASN', 'SER', 'ARG', 'TRP', 'TRP', 'CYS', 'ASN', 'ASP', 'GLY', 'ARG', 'THR', 'PRO', 'GLY', 'SER', 'ARG', 'ASN', 'LEU', 'CYS', 'ASN', 'ILE', 'PRO', 'CYS', 'SER', 'ALA', 'LEU', 'LEU', 'SER', 'SER', 'ASP', 'ILE', 'THR', 'ALA', 'SER', 'VAL', 'ASN', 'CYS', 'ALA', 'LYS', 'LYS', 'ILE', 'VAL', 'SER', 'ASP', 'GLY', 'ASN', 'GLY', 'MET', 'ASN', 'ALA', 'TRP', 'VAL', 'ALA', 'TRP', 'ARG', 'ASN', 'ARG', 'CYS', 'LYS', 'GLY', 'THR', 'ASP', 'VAL', 'GLN', 'ALA', 'TRP', 'ILE', 'ARG', 'GLY', 'CYS', 'ARG', 'LEU']

## exclusion/restriction
spacegroup = 'P 43 21 2'
n_mutations_max = 1

l_ligands_exclude = ['IOD','NAG','NDG','SCN',]

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


def main():

    l_pdbs = core.parse_pdbs_from_cluster(i_cluster,)
    l_pdbs = l_pdbs[320:]

    d_mmCIF, d_seq, l_pdbs = core.parse_cifs(
        l_pdbs,spacegroup,l_ligands_exclude,ref_seq,n_mutations_max,
        )

    d_coordinates = core.parse_coordinates(
        l_pdbs,d_mmCIF,method,
        )
    
##    d_bfactors = normalize_bfactors(d_coordinates,)

    d_rmsds_overall, d_rmsds_subset = core.calculate_rmsds(l_pdbs,d_coordinates,d_seq,method,do_subset=True)

    fd = open('d_rmsds_subset','w')
    fd.write('%s' %(d_rmsds_subset))
    fd.close()

    sample_and_plot_rmsds(d_rmsds_subset,d_seq,)

    return


def normalize_bfactors(d_coordinates,):

    d_bfactors = {}

    for pdb in d_coordinates.keys():
        d_bfactors[pdb] = {}
        for res_no in d_coordinates[pdb].keys():
            stop
        stop

    return d_bfactors


def sample_and_plot_rmsds(d_rmsds,d_wt,):

    print 'sample and plot'

    for prop1 in d_properties.keys():
        
        l_properties2 = d_properties[prop1]

        prefix = '%i_%s_%s' %(i_cluster,method,prop1,)
        
        for k in ['wt','mutant',]:

            l_data = []
            print prop1,k

            for pdb in d_rmsds[k].keys():
                ## new data file for each mutant
                if k == 'mutant':
                    l_data = []
                    print pdb
                    i_pdb = l_mutants.index(pdb)
                    if d_wt[pdb] == 'wt':
                        print pdb
                        stop
                
                for i_prop2 in range(len(l_properties2)):
                    prop2 = l_properties2[i_prop2]
                    l_rmsds = []

                    ## sample all wts
                    for pdb_wt in d_rmsds[k][pdb][prop1][prop2].keys():
                        if d_wt[pdb_wt] != 'wt':
                            print pdb_wt
                            stop
                        rmsd = d_rmsds[k][pdb][prop1][prop2][pdb_wt]
                        l_rmsds += [rmsd]

                    for i in range(len(l_rmsds)):
                        l_rmsds_sampling_pop = list(l_rmsds)
                        l_rmsds_sampling_append = []
                        for j in range(len(l_rmsds)/2):
                            index = int(random.random()*len(l_rmsds_sampling_pop))
                            l_rmsds_sampling_append += [l_rmsds_sampling_pop[index]]
                            del l_rmsds_sampling_pop[index]
                        rmsd_average = sum(l_rmsds_sampling_append)/len(l_rmsds_sampling_append)

                        if k == 'wt':
                            x = 2*i_prop2+0
                        else:
                            x = 2*i_prop2+float(i_pdb+1)/len(l_mutants)
                        y = rmsd_average

                        l_data += ['%s %s\n' %(x,y,)]

                ## write data file for each mutant
                if k == 'mutant':
                    fd = open('%s_%s.gnuplot_data' %(prefix,pdb,),'w')
                    fd.writelines(l_data)
                    fd.close()

            ## write one data for all the wts
            if k == 'wt':
                fd = open('%s_wt.gnuplot_data' %(prefix,),'w')
                fd.writelines(l_data)
                fd.close()

        ## plot wt and mutant
        ylabel = 'average rmsd'
        title = prefix.replace('_',' ')
        gnuplot(prefix,ylabel,title,l_properties2,)

    os.remove('%s_wt.gnuplot_data' %(prefix))
    for mutant in l_mutants:
        if os.path.isfile('%s_%s.gnuplot_data' %(prefix,mutant,)):
            os.remove('%s_%s.gnuplot_data' %(prefix,mutant,))

    return


def gnuplot(prefix,ylabel,title,l_properties,):

    l_settings = [
        'set terminal postscript eps enhanced color "Helvetica" 36\n',
        'set output "%s.ps"\n' %(prefix),
##        'set terminal png\n',
##        'set output "gnuplot.png"\n',
        'set size 3,3\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set title "%s"\n' %(title,),
        'set ylabel "%s"\n' %(ylabel),
        'set xlabel "%s"\n' %(l_properties,),
        ]
    ## scatter plot
    s = 'plot [-1:%i][0:]' %(2*len(l_properties))
    s += '"%s_wt.gnuplot_data" u 1:2 lw 2 ps 4 lc 0 pt 5 t "wt", ' %(prefix)
    for mutant in l_mutants:
        s += '"%s_%s.gnuplot_data" u 1:2 lw 2 ps 4 pt 7 lc rgb "#%6s" t "%s", ' %(
            prefix, mutant, d_mutants[mutant]['color'], d_mutants[mutant]['mutation'],
            )
    s = s[:-2]+'\n'
    l_settings += [s]

    fd = open('%s.gnuplot_settings' %(prefix),'w')
    fd.writelines(l_settings)
    fd.close()

    os.system('/software/bin/gnuplot %s.gnuplot_settings' %(prefix))

    print 'imagemagick', prefix
    os.system('convert %s.ps %s.png' %(prefix,prefix,))
    os.system('rm %s.ps' %(prefix))
    os.remove('%s.gnuplot_settings' %(prefix))

    return    


if __name__ == '__main__':
    main()

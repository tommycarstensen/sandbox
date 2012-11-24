## HEWL scatter plots of RMSD vs misc properties (pH diff, T diff, max resol)

## imports
import os,numpy
import core
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb')
import parse_mmCIF
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import gnuplot, statistics
sys.path.append('/home/people/tc/svn/Protool/')
import geometry

instance_geometry = geometry.geometry()

## appendix
d_321 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
    'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
    'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
    'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
    }

## exclude/restrict
l_ligands_exclude = [
    'IOD', ## 1b2k, 1lkr
    'NAG','NDG', ## substrate
    'SCN', ## 1lcn, thiocyanate
    'CCN', ## acetonitrile cross linker (2lyo,3lyo)
    ]

n_mutations_max = 0
n_mutations_max = 1

##spacegroup = None
spacegroup = 'P 43 21 2'

##
i_cluster = 5
##s_uniprot = 'P00698'
##fd = open('uniprot_sprot.fasta','r')
##lines = fd.readlines()
##fd.close()
##for i in range(len(lines)):
##    line = lines[i]
##    if line[0] == '>':
##        if line[4:10] == s_uniprot:
##            s_seq = ''
##            for j in range(i+1,len(lines)):
##                if lines[j][0] == '>':
##                    break
##                s_seq += lines[j].strip()

ref_seq = 'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL'
ref_seq = ['LYS', 'VAL', 'PHE', 'GLY', 'ARG', 'CYS', 'GLU', 'LEU', 'ALA', 'ALA', 'ALA', 'MET', 'LYS', 'ARG', 'HIS', 'GLY', 'LEU', 'ASP', 'ASN', 'TYR', 'ARG', 'GLY', 'TYR', 'SER', 'LEU', 'GLY', 'ASN', 'TRP', 'VAL', 'CYS', 'ALA', 'ALA', 'LYS', 'PHE', 'GLU', 'SER', 'ASN', 'PHE', 'ASN', 'THR', 'GLN', 'ALA', 'THR', 'ASN', 'ARG', 'ASN', 'THR', 'ASP', 'GLY', 'SER', 'THR', 'ASP', 'TYR', 'GLY', 'ILE', 'LEU', 'GLN', 'ILE', 'ASN', 'SER', 'ARG', 'TRP', 'TRP', 'CYS', 'ASN', 'ASP', 'GLY', 'ARG', 'THR', 'PRO', 'GLY', 'SER', 'ARG', 'ASN', 'LEU', 'CYS', 'ASN', 'ILE', 'PRO', 'CYS', 'SER', 'ALA', 'LEU', 'LEU', 'SER', 'SER', 'ASP', 'ILE', 'THR', 'ALA', 'SER', 'VAL', 'ASN', 'CYS', 'ALA', 'LYS', 'LYS', 'ILE', 'VAL', 'SER', 'ASP', 'GLY', 'ASN', 'GLY', 'MET', 'ASN', 'ALA', 'TRP', 'VAL', 'ALA', 'TRP', 'ARG', 'ASN', 'ARG', 'CYS', 'LYS', 'GLY', 'THR', 'ASP', 'VAL', 'GLN', 'ALA', 'TRP', 'ILE', 'ARG', 'GLY', 'CYS', 'ARG', 'LEU']

##
method = 'heavy'
##method = 'alpha'

def main():

    l_pdbs = core.parse_pdbs_from_cluster(i_cluster,)
    l_pdbs = l_pdbs[300:]

    d_mmCIF, d_wt = core.parse_cifs(
        l_pdbs,spacegroup,l_ligands_exclude,ref_seq,n_mutations_max,
        )

    tabulate(d_mmCIF,)

    d_coordinates = core.parse_coordinates(l_pdbs,d_mmCIF,method,)
    d_rmsds = calculate_rmsds(d_coordinates,)

    fd = open('rmsds.txt','w')
    fd.write('%s' %(d_rmsds))
    fd.close()

    plot(d_mmCIF,d_rmsds,)

    return


def tabulate(d_mmCIF_main,):

    l_pdbs = list(set(d_mmCIF_main.keys()))
    l_pdbs.sort()

    l = ['ID   mut       T    pH  res spacegroup startingmodel\n']
    for i in range(len(l_pdbs)):
        pdb = l_pdbs[i]
        d_mmCIF = d_mmCIF_main[pdb]
        if d_mmCIF['_entity_poly_seq.mon_id'] == ref_seq:
            s_mutation = 'wt'
        else:
            for i in range(len(d_mmCIF['_entity_poly_seq.mon_id'])):
                res_id_mmCIF = d_mmCIF['_entity_poly_seq.mon_id'][i]
                res_id_ref = ref_seq[i]
                if res_id_mmCIF != res_id_ref:
                    res_no = i+1
                    s_mutation = '%1s%3i%1s' %(d_321[res_id_ref],res_no,d_321[res_id_mmCIF],)
        spacegroup = core.parse_mmCIF_item(d_mmCIF,'_symmetry.space_group_name_H-M',)
        T = core.parse_mmCIF_item(d_mmCIF,'_diffrn.ambient_temp',)
        pH = core.parse_mmCIF_item(d_mmCIF,'_exptl_crystal_grow.pH',)
        starting_model = core.parse_mmCIF_item(d_mmCIF,'_refine.pdbx_starting_model',)
        resolution = core.parse_mmCIF_item(d_mmCIF,'_refine.ls_d_res_high',)
        l += ['%4s %-5s %5s %5s %3.1f %-10s %4s\n' %(pdb,s_mutation,T,pH,float(resolution),spacegroup,starting_model,)]

    fd = open('HEWL_table.txt','w')
    fd.writelines(l)
    fd.close()

    return


def plot(d_mmCIF_main,d_rmsds,):

    l_pdbs = d_rmsds.keys()
    l_pdbs.sort()

    l_temperature = []
    l_ph = []
    l_resolution = []
    d_spacegroup = {}
    d_starting_model = {}

    l_correl_T = [[],[],]
    l_correl_pH = [[],[],]
    l_correl_resol_max = [[],[],]

    d_histo_pH = {}
    d_histo_T = {}
    d_histo_resol = {}

    for i1 in range(len(l_pdbs)-1):
        pdb1 = l_pdbs[i1]
        spacegroup1 = core.parse_mmCIF_item(d_mmCIF_main[pdb1[:4]],'_symmetry.space_group_name_H-M',)
        T1 = core.parse_mmCIF_item(d_mmCIF_main[pdb1[:4]],'_diffrn.ambient_temp',)
        pH1 = core.parse_mmCIF_item(d_mmCIF_main[pdb1[:4]],'_exptl_crystal_grow.pH',)
        starting_model1 = core.parse_mmCIF_item(d_mmCIF_main[pdb1[:4]],'_refine.pdbx_starting_model',)
        resolution1 = core.parse_mmCIF_item(d_mmCIF_main[pdb1[:4]],'_refine.ls_d_res_high',)

        for i2 in range(i1+1,len(l_pdbs)):
            pdb2 = l_pdbs[i2]
            spacegroup2 = core.parse_mmCIF_item(d_mmCIF_main[pdb2[:4]],'_symmetry.space_group_name_H-M',)
            T2 = core.parse_mmCIF_item(d_mmCIF_main[pdb2[:4]],'_diffrn.ambient_temp',)
            pH2 = core.parse_mmCIF_item(d_mmCIF_main[pdb2[:4]],'_exptl_crystal_grow.pH',)
            starting_model2 = core.parse_mmCIF_item(d_mmCIF_main[pdb2[:4]],'_refine.pdbx_starting_model',)
            resolution2 = core.parse_mmCIF_item(d_mmCIF_main[pdb2[:4]],'_refine.ls_d_res_high',)

            rmsd = d_rmsds[pdb1][pdb2]
            if rmsd > 1:
                print pdb1, pdb2, rmsd

            if T1 and T2:
                T_diff = abs(float(T2)-float(T1))
                l_temperature += ['%s %s\n' %(T_diff,rmsd),]
                l_correl_T[0] += [T_diff]
                l_correl_T[1] += [rmsd]

                print T_diff, 10*round(T_diff/10.,0)
                if not 10*round(T_diff/10.,0) in d_histo_T.keys():
                    d_histo_T[10*round(T_diff/10.,0)] = 0
                d_histo_T[10*round(T_diff/10.,0)] += 1

            if pH1 and pH2:
                pH_diff = abs(float(pH2)-float(pH1))
                l_ph += ['%s %s\n' %(pH_diff,rmsd),]
                l_correl_pH[0] += [pH_diff]
                l_correl_pH[1] += [rmsd]

                if not pH_diff in d_histo_pH.keys():
                    d_histo_pH[pH_diff] = 0
                d_histo_pH[pH_diff] += 1

            resolution_max = max(resolution1,resolution2,)
            l_resolution += ['%s %s\n' %(resolution_max,rmsd),]
            if resolution_max != 'N/A':
                l_correl_resol_max[0] += [float(resolution_max)]
                l_correl_resol_max[1] += [rmsd]

                if not round(float(resolution_max),0) in d_histo_resol.keys():
                    d_histo_resol[round(float(resolution_max),0)] = 0
                d_histo_resol[round(float(resolution_max),0)] += 1

            d_spacegroup = append_to_dictionary(d_spacegroup,spacegroup1,spacegroup2,rmsd,)
            d_starting_model = append_to_dictionary(d_starting_model,starting_model1,starting_model2,rmsd,)

    r1 = statistics.correlation(l_correl_T[0],l_correl_T[1],)
    r2 = statistics.correlation(l_correl_pH[0],l_correl_pH[1],)
    r3 = statistics.correlation(l_correl_resol_max[0],l_correl_resol_max[1],)

    ##
    ## plot histograms
    ##
    for prefix,d in [
        ['deltapH',d_histo_pH,],
        ['deltaT',d_histo_T,],
        ['maxresolution',d_histo_resol,],
        ]:
        
        l = []
        l_diffs = d.keys()
        l_diffs.sort()
        for diff in l_diffs:
            l += ['%s %s\n' %(diff,d[diff],)]
        fd = open('histo_%s.txt' %(prefix),'w')
        fd.writelines(l)
        fd.close()

        l = [
            'set terminal postscript eps enhanced color "Helvetica"\n',
            'set output "gnuplot.ps"\n',
            'set size 3,3\n',
            'set style data histogram\n',
            'set xtics rotate\n',
            'set xlabel "%s\n' %(prefix),
            'set ylabel "count\n',
            'plot "histo_%s.txt" u 2:xtic(1) t ""\n' %(prefix)
            ]
        fd = open('tmp.txt','w')
        fd.writelines(l)
        fd.close()

        os.system('gnuplot tmp.txt')
        os.system('convert gnuplot.ps histo_%s.png' %(prefix))

    ##
    ## plot rmsd as a function of each property (2d)
    ##
    for prefix,data,xlabel in [
        ['pH',l_ph,'pH diff',],
        ['Temperature',l_temperature,'T diff',],
        ['resolution',l_resolution,'maximum resolution',],
        ]:
        prefix += method
        fd = open('%s.gnuplotdata' %(prefix),'w')
        fd.writelines(data)
        fd.close()
        gnuplot.scatter_plot_2d(
            prefix,xlabel=xlabel,ylabel='RMSD %s' %(method,),
##            averages=True,
            regression=True,
            )

    ##
    ## plot rmsd as a function of each property (contour)
    ##
    for d,prefix in [
        [d_spacegroup,'spacegroup',],
        [d_starting_model,'startingmodel',],
        ]:

        d_tics = {}
        l_tics = d.keys()
        l_tics.sort()
        for i in range(len(l_tics)):
            d_tics[l_tics[i]] = i+.5
        z1 = 9
        z2 = 0

        l_data = []
        for x in range(len(l_tics)):
            k1 = l_tics[x]
            for y in range(len(l_tics)):
                k2 = l_tics[y]
                if not k2 in d[k1].keys():
                    average = 9
                else:
                    l_rmsds = d[k1][k2]
                    average = sum(l_rmsds)/len(l_rmsds)
                    if average < z1:
                        z1 = average
                    if average > z2:
                        z2 = average
                l_data += ['%s %s %s\n' %(x,y,average,)]
            l_data += ['%s %s %s\n' %(x,y+1,1,)]
            l_data += ['\n']
        for y in range(len(l_tics)):
            l_data += ['%s %s %s\n' %(x+1,y,1,)]
        l_data += ['%s %s %s\n' %(x+1,y+1,1,)]
        l_data += ['\n']
        gnuplot.contour_plot(
            prefix,l_data,
            title='%s %s' %(prefix,method,),zlabel='RMSD %s' %(method),
            d_xtics = d_tics, d_ytics = d_tics,
            palette = '0 1 0 0, 0.9999 0 0 1, 0.9999 1 1 1, 1 1 1 1',
            z1 = z1, z2 = z2+0.1,
            bool_remove = False,
            )
        os.system('convert %s.ps %s_spacegroup%s_mutations%s_atoms%s.png' %(prefix,prefix,spacegroup.replace(' ',''),n_mutations_max,method,))
##        os.remove('%s.ps' %(prefix,))

    print d_spacegroup
    print d_starting_model

    print r1
    print r2
    print r3

    return


def append_to_dictionary(d,k1,k2,rmsd,):
    
    if not k1 in d.keys():
        d[k1] = {}
    if not k2 in d[k1].keys():
        d[k1][k2] = []
    if not k2 in d.keys():
        d[k2] = {}
    if not k1 in d[k2].keys():
        d[k2][k1] = []
    d[k1][k2] += [rmsd]
    d[k2][k1] += [rmsd]

    return d


def calculate_rmsds(d_coordinates,):

    print 'calculate rmsds'

    d_rmsds = {}
    l_pdbs = d_coordinates.keys()
    l_pdbs.sort()

    for pdb in l_pdbs:
        d_rmsds[pdb] = {}

    for i1 in range(len(l_pdbs)-1):
        pdb1 = l_pdbs[i1]

        print 'calculate rmsds', pdb1

        for i2 in range(i1+1,len(l_pdbs)):
            pdb2 = l_pdbs[i2]

            l_coordinates1 = []
            l_coordinates2 = []
            set_res_nos1 = set(d_coordinates[pdb1].keys())
            set_res_nos2 = set(d_coordinates[pdb2].keys())
            set_res_nos = set_res_nos1 & set_res_nos2
            for res_no in set_res_nos:
                set_res_names1 = set(d_coordinates[pdb1][res_no].keys())
                set_res_names2 = set(d_coordinates[pdb2][res_no].keys())
                set_res_names = set_res_names1 & set_res_names2
                ## no mutation
                if len(set_res_names) == 1:
                    for res_name in set_res_names:
                        set_atom_names1 = set(d_coordinates[pdb1][res_no][res_name].keys())
                        set_atom_names2 = set(d_coordinates[pdb2][res_no][res_name].keys())
                        set_atom_names = set_atom_names1 & set_atom_names2
                        for atom_name in set_atom_names:
                            coord1 = d_coordinates[pdb1][res_no][res_name][atom_name]['coord']
                            coord2 = d_coordinates[pdb2][res_no][res_name][atom_name]['coord']
                            l_coordinates1 += [coord1]
                            l_coordinates2 += [coord2]
                ## mutation
                elif len(set_res_names) == 0:
                    if len(set_res_names1) == 1 and len(set_res_names2) == 1:
                        res_name1 = list(set_res_names1)[0]
                        res_name2 = list(set_res_names2)[0]
                        set_atom_names1 = set(d_coordinates[pdb1][res_no][res_name1].keys())
                        set_atom_names2 = set(d_coordinates[pdb2][res_no][res_name2].keys())
                        set_atom_names = set(['N','CA','CB','C','O',]) & set_atom_names1 & set_atom_names2
                        for atom_name in set_atom_names:
                            coord1 = d_coordinates[pdb1][res_no][res_name1][atom_name]['coord']
                            coord2 = d_coordinates[pdb2][res_no][res_name2][atom_name]['coord']
                            l_coordinates1 += [coord1]
                            l_coordinates2 += [coord2]
                    else:
                        print set_res_names1
                        print set_res_names2
                        stop
                ## altloc
                else:
                    print set_res_names
                    stop
            rmsd = instance_geometry.superpose(l_coordinates1,l_coordinates2,)

            d_rmsds[pdb1][pdb2] = rmsd
            d_rmsds[pdb2][pdb1] = rmsd

    return d_rmsds


if __name__ == '__main__':
    main()

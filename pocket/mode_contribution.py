import numpy, math, os
import sys
sys.path.append('/home/tc/svn/GoodVibes')
import goodvibes_ligand
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('/home/tc/svn/Protool/')
import geometry
sys.path.append('/home/tc/svn/tc_sandbox/math')
import statistics

def main():

    ## dictionary of apo and holo structures (from what script???)
    d_apo2holo = {
        ## conformational selection
        ## RNase A, 1kf3 high resolution
##        '1kf3': {'holo': '1rpg','ligand':'CPA',},
        '1kf5': {
            'holo': '1eow','ligand':'U2G',
            'site':[
                11,43,44,
##                119,120,121,122, ## terminal flexible residues...
                ],
            'title':'Ribonuclease (RNase) A',
            },
        ## CypA, highest resolution room
##        '3k0n': {'holo': '1cwa','ligand':['DAL','MLE','MVA','BMT','ABA','SAR',],'site':[18-1,54-1,59-1,62-1,71-1,100-1,101-1,102-1,110-1,112-1,120-1,121-1,125-1,163-1,],'title':'Peptidyl-prolyl isomerase A (CypA)',},
        '1w8v': {'holo': '1w8m','ligand':['E1P',],'site':[54,62,100,101,112,125,],'title':'Peptidyl-prolyl isomerase A (CypA)',},
        ## DHFR
        '1ra9': {'holo': '1ra2','ligand':'FOL','site':[4,5,6,26,27,30,31,56,93,112,],'title':'Dihydrofolate reductase (DHFR)',},
        ## AdK
        '2rh5': {'holo': '2rgx','ligand':'AP5','coords_apo':[0,202],'coords_holo':[0,202],'site':[8,9,10,11,12,13,14,30,31,34,57,58,63,81,84,88,119,120,123,134,137,149,160,188,189,190,],'title':'Adenylate kinase (AdK)',},
        ## PKA
        '3iia': {'holo': '3pna', 'ligand':'CMP', 'coords_apo':[4,133],'coords_holo':[0,129],'site':[
            182-112, 198-112, 199-112, 200-112, 201-112, 202-112, 209-112, 211-112,
            ],
                 'title':'Protein Kinase A (PKA)',
                 },
        ## induced fit
        ## PEPCK
        '2qew': {'holo': '3dt4', 'ligand':'OXL', 'coords_apo':[1,620],'coords_holo':[0,619],'site':[240,260,307,401],'title':'PEPCK',},
        ## beta-lactoglobulin
        '3npo': {'holo': '3nq3', 'ligand':'DKA','site':[53,104,106,],'title':'beta-lactoglobulin',},
        }

    for pdb_apo in d_apo2holo.keys():

        pdb_holo = d_apo2holo[pdb_apo]['holo']

##        continue ## tmp!!!
        print pdb_apo, pdb_holo

        ##
        ## parse coordinates
        ##
        d_mmCIF_apo, l_coords_alpha_apo = parse_coords(pdb_apo)
        d_mmCIF_holo, l_coords_alpha_holo = parse_coords(pdb_holo)
        if 'coords_apo' in d_apo2holo[pdb_apo].keys():
            l_coords_alpha_apo = l_coords_alpha_apo[
                d_apo2holo[pdb_apo]['coords_apo'][0]:d_apo2holo[pdb_apo]['coords_apo'][1]
                ]
            l_coords_alpha_holo = l_coords_alpha_holo[
                d_apo2holo[pdb_apo]['coords_holo'][0]:d_apo2holo[pdb_apo]['coords_holo'][1]
                ]
        else:
            ## sequential alignment of coordinates
            index1_seq_apo = next((i for i,v in enumerate(d_mmCIF_apo['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
            index1_seq_holo = next((i for i,v in enumerate(d_mmCIF_holo['_pdbx_poly_seq_scheme.pdb_mon_id']) if v != '?'))
            ## last non-?
            index2_seq_apo = len(d_mmCIF_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d_mmCIF_apo['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))
            index2_seq_holo = len(d_mmCIF_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])-next((i for i,v in enumerate(reversed(d_mmCIF_holo['_pdbx_poly_seq_scheme.pdb_mon_id'])) if v != '?'))
            ## first common non-?
            index1_coord_apo = max(0,index1_seq_holo-index1_seq_apo)
            index1_coord_holo = max(0,index1_seq_apo-index1_seq_holo)
            ## last common non-?
            index2_coord_apo = len(l_coords_alpha_apo)+min(0,index2_seq_holo-index2_seq_apo)
            index2_coord_holo = len(l_coords_alpha_holo)+min(0,index2_seq_apo-index2_seq_holo)
            l_coords_alpha_apo = l_coords_alpha_apo[index1_coord_apo:index2_coord_apo]
            l_coords_alpha_holo = l_coords_alpha_holo[index1_coord_holo:index2_coord_holo]

        if pdb_apo == '2qew' and pdb_holo == '3dt4':
##            l_coords_alpha_holo = l_coords_alpha_holo[:459]+l_coords_alpha_holo[466:]
            l_coords_alpha_holo = l_coords_alpha_holo[:461]+l_coords_alpha_holo[461+7:]

        if len(l_coords_alpha_apo) != len(l_coords_alpha_holo):
            print pdb_apo, pdb_holo
            print len(l_coords_alpha_apo), len(l_coords_alpha_holo)
            stop

##        if pdb_holo == '1eow':
##            print l_coords_alpha_holo[d_apo2holo[pdb_apo]['site'][0]]
##            print pdb_holo
##            stop

        tv1, rm, tv2, l_coords_alpha_apo, l_coords_alpha_holo = get_transformation_matrix(
            l_coords_alpha_apo,
            l_coords_alpha_holo,
            )

        vector_apo2holo = get_apo_holo_vector(
            d_mmCIF_apo, l_coords_alpha_apo,
            d_mmCIF_holo, l_coords_alpha_holo,
            tv1, rm, tv2,
            )

        chain_apo = ''.join(d_mmCIF_apo['_entity_poly.pdbx_strand_id'])
        chain_holo = ''.join(d_mmCIF_holo['_entity_poly.pdbx_strand_id'])

        if pdb_holo == '1cwa':
            ligand_pos_holo = numpy.array([3.307729,36.55456,17.45886])
            ligand_pos_apo = numpy.dot(ligand_pos_holo-tv1,rm)+tv2
        else:
            ligand_pos_apo, ligand_pos_holo, lines_ligand_apo = get_ligand_pos(
                d_mmCIF_holo,
                tv1, rm, tv2,
                d_apo2holo[pdb_apo]['ligand'],
                pdb_holo,
                )

        dist_max = 6
        dist_min = 3
##        print len(vector_apo2holo), len(l_coords_alpha_apo)
##        stop
        for pdb, chain, l_coords_alpha, ligand_pos in [
##            [pdb_holo,chain_holo,l_coords_alpha_holo,],
            [pdb_apo,chain_apo,l_coords_alpha_apo,ligand_pos_apo],
            ]:
##            l_coords_protein_alpha = []
##            for i in range(len(l_coords_alpha)):
##                l_coords_protein_alpha += [l_coords_alpha[i][0]]
##                l_coords_protein_alpha += [l_coords_alpha[i][1]]
##                l_coords_protein_alpha += [l_coords_alpha[i][2]]
            fn = '/home/tc/UCD/GV_ligand_binding_site_identification/%s_%s_probe.pdb' %(pdb,chain,)
            if os.path.isfile(fn):
                continue
            d = goodvibes_ligand.main(
                pdb,chain,
                dist_max,dist_min,
                v_apoholo=vector_apo2holo,
                l_coords_protein_alpha = l_coords_alpha,
##                l_coords_probe = [ligand_pos],
                )
##            d = goodvibes_ligand.main(
##                pdb,chain,
##                dist_max,dist_min,
##                v_apoholo=vector_apo2holo,
##                l_coords_protein_alpha = l_coords_alpha,
##                l_coords_probe = [ligand_pos],
##                )
            l_factors = d['l_factors']
##            l_factors_perturbed = d['l_factors_probe']

        if os.path.isfile(fn):
            fd = open(fn,'r')
            lines = fd.readlines()
            fd.close()
            lines += lines_ligand_apo
            fd = open(fn+'2','w')
            fd.writelines(lines)
            fd.close()
            continue

        eigenvectors = d['eigenvectors']
        l_factors_abs = [abs(factor) for factor in l_factors]
        mode_max_contribution = l_factors_abs.index(max(l_factors_abs))
        print mode_max_contribution

        print d['l_overlaps']

        v1 = vector_apo2holo
        eigenvector = v2 = eigenvectors[mode_max_contribution]
        overlap_max = abs(numpy.dot(v1,v2))/math.sqrt(numpy.dot(v1,v1)*numpy.dot(v2,v2))
        print 'mode_max_contribution', mode_max_contribution
        print 'overlap_max', overlap_max

        ## write amplitudes
        lines = []
        l1 = []
        l2 = []
        for i in range(0,len(eigenvector),3):
            amplitude = math.sqrt(eigenvector[i+0]**2+eigenvector[i+1]**2+eigenvector[i+2]**2)
            amplitude7 = math.sqrt(eigenvectors[6][i+0]**2+eigenvectors[6][i+1]**2+eigenvectors[6][i+2]**2)
            l1 += [amplitude]
            l2 += [amplitude7]
            if i/3 in d_apo2holo[pdb_apo]['site']:
                bool_site = amplitude
            else:
                bool_site = -1
            lines += ['%s %s %s\n' %(amplitude,amplitude7,bool_site,)]
        fd = open('amplitudes_%s%s.txt' %(pdb_apo,pdb_holo,),'w')
        fd.writelines(lines)
        fd.close()

        r = statistics.correlation(l1,l2)

        xmin = 0
        if pdb_holo == '1w8m':
            xmin = 2
        if pdb_holo == '3dt4':
            xmin = 5
        if pdb_holo == '1eow':
            xmin = 1
        lines = [
            'set terminal png\n',
            'set output "%s%s_amplitudes.png"\n' %(pdb_apo,pdb_holo,),
            'set size 1,1\n',
            'set xlabel "residue index"\n',
            'set ylabel "amplitude (a.u.)\n',
            'set title "%s (r = %.2f)"\n' %(d_apo2holo[pdb_apo]['title'],r,),
            'set key out\n',
            'f(x) = %s\n' %(sum(l1)/len(l1)),
            'plot [%s:][0:]"amplitudes_%s%s.txt" u 1 t "mode %i" w l, "amplitudes_%s%s.txt" u 2 t "mode 7" w l, "amplitudes_%s%s.txt" u 3 t "binding site" ps 1 pt 7, f(x) t "average amplitude"\n' %(
##            'plot [%s:][0:]"amplitudes_%s%s.txt" u 1 t "mode %i" w l, "amplitudes_%s%s.txt" u 2 t "mode 7" w l, "amplitudes_%s%s.txt" u 3 t "binding site" ps 1 pt 7\n' %(
                xmin,pdb_apo,pdb_holo,mode_max_contribution+1, pdb_apo,pdb_holo, pdb_apo,pdb_holo,
                ),
            ]
        fd = open('gnuplot.settings','w')
        fd.writelines(lines)
        fd.close()
        os.system('/usr/bin/gnuplot gnuplot.settings')

        s = ''
        for i in range(len(l_factors)):
            s += '%s %s %s\n' %(i+1, l_factors[i],abs(l_factors[i]),)
        fd = open('facs_eigvals_%s%s.txt' %(pdb_apo,pdb_holo,),'w')
        fd.write(s)
        fd.close()

##        s = ''
##        for i in range(len(l_factors_perturbed)):
##            s += '%s %s %s\n' %(i+1, l_factors_perturbed[i],abs(l_factors_perturbed[i]),)
##        fd = open('facs_eigvals_%s%s_perturbed.txt' %(pdb_apo,pdb_holo,),'w')
##        fd.write(s)
##        fd.close()
    
    return


def get_ligand_pos(
    d_mmCIF_holo,
    tv1, rm, tv2,
    ligand,
    pdb_holo,
    ):

    l_coords_ligand = []
    l_coords_ligand_apo = []
    lines_ligand_apo = []
    for i in range(len(d_mmCIF_holo['_atom_site.id'])):
        if (
            d_mmCIF_holo['_atom_site.group_PDB'][i] == 'HETATM'
            and
            (
                d_mmCIF_holo['_atom_site.label_comp_id'][i] == ligand
                or
                d_mmCIF_holo['_atom_site.label_comp_id'][i] in ligand
                )
            ):
            x = float(d_mmCIF_holo['_atom_site.Cartn_x'][i])
            y = float(d_mmCIF_holo['_atom_site.Cartn_y'][i])
            z = float(d_mmCIF_holo['_atom_site.Cartn_z'][i])
            coord = numpy.array([x,y,z,])
            l_coords_ligand += [coord]
            ## holo2apo
            coord_apo = numpy.dot(coord-tv1,rm)+tv2
##            rm_transpose = numpy.transpose(rm)
##            coord_apo = numpy.dot(coord-tv2,rm_transpose)+tv1
            l_coords_ligand_apo += [coord_apo]
            line = 'HETATM%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
                int(d_mmCIF_holo['_atom_site.id'][i]),
                d_mmCIF_holo['_atom_site.label_atom_id'][i],
                d_mmCIF_holo['_atom_site.label_comp_id'][i],
                d_mmCIF_holo['_atom_site.label_asym_id'][i],
                int(d_mmCIF_holo['_atom_site.label_seq_id'][i].replace('.','0')),
                coord_apo[0],coord_apo[1],coord_apo[2],
                1.,0.,
                d_mmCIF_holo['_atom_site.type_symbol'][i],
                )
            lines_ligand_apo += [line]

    position_ligand_apo = sum(l_coords_ligand_apo)/len(l_coords_ligand_apo)
    position_ligand = sum(l_coords_ligand)/len(l_coords_ligand)

    return position_ligand_apo, position_ligand, lines_ligand_apo


def get_transformation_matrix(
    l_coords_alpha_apo,
    l_coords_alpha_holo,
    ):

    instance_geometry = geometry.geometry()
    rmsd = instance_geometry.superpose(l_coords_alpha_apo,l_coords_alpha_holo,)
    print 'rmsd', rmsd
    tv1 = instance_geometry.fitcenter
    rm = instance_geometry.rotation
    tv2 = instance_geometry.refcenter

    return tv1, rm, tv2, l_coords_alpha_apo, l_coords_alpha_holo


def get_apo_holo_vector(
    d_mmCIF_apo, l_coords_alpha_apo,
    d_mmCIF_holo, l_coords_alpha_holo,
    tv1,rm,tv2,
    ):

    ##
    ## apply transformation matrix (holo2apo)
    ##
    for i_coord in range(len(l_coords_alpha_holo)):
        coord = l_coords_alpha_holo[i_coord]
        coord = numpy.dot(coord-tv1,rm)+tv2
        l_coords_alpha_holo[i_coord] = coord

    ##
    ## apo/holo eigenvector
    ##
    vector_apo2holo = []
    for i in range(len(l_coords_alpha_holo)):
        vector_apo2holo += [
            l_coords_alpha_holo[i][0]-l_coords_alpha_apo[i][0],
            l_coords_alpha_holo[i][1]-l_coords_alpha_apo[i][1],
            l_coords_alpha_holo[i][2]-l_coords_alpha_apo[i][2],
            ]
    vector_apo2holo = numpy.array(vector_apo2holo)

    return vector_apo2holo


def parse_coords(pdb):

    d_mmCIF = parse_mmCIF.main(pdb,)
    d_coords, l_coords_alpha = mmCIF2coords.main(pdb,d_mmCIF)

    return d_mmCIF, l_coords_alpha


if __name__ == '__main__':
    main()
l_pdbs = [
    [''],
    ]

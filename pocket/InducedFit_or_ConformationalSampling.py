import numpy
import sys
sys.path.append('/home/tc/svn/GoodVibes')
import goodvibes_ligand
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('/home/tc/svn/Protool/')
import geometry

def main():

    d_cath = parse_CATH()

    ## dictionary of apo and holo structures (from what script???)
    d_apo2holo = {
        '1ak1': {'holo': '1c1h', 'rmsd': 0.64, 'ligand': 'MMP', 'overlap': 0.36},
        '1lp8': {'holo': '1lpc', 'rmsd': 0.14003320157149318, 'ligand': 'CMP', 'overlap': 0.1663339954531419}, '2exo': {'holo': '1j01', 'rmsd': 0.31317202378931264, 'ligand': 'XIL', 'overlap': 0.5464584014950491}, '1syc': {'holo': '1syd', 'rmsd': 0.42370903255794784, 'ligand': 'THP', 'overlap': 0.312528957121447}, '1cua': {'holo': '3esd', 'rmsd': 0.4040397349113579, 'ligand': 'SXC', 'overlap': 0.08139930257927935}, '1znw': {'holo': '1znx', 'rmsd': 0.24704558797447082, 'ligand': '5GP', 'overlap': 0.12436531581580901}, '1gbs': {'holo': '1lsp', 'rmsd': 0.264857015278693, 'ligand': 'BUL', 'overlap': 0.2934470556757598}, '1sqg': {'holo': '1sqf', 'rmsd': 0.6990763011194481, 'ligand': 'SAM', 'overlap': 0.4723483760294988}, '3npo': {'holo': '3nq3', 'rmsd': 0.7183351557421368, 'ligand': 'DKA', 'overlap': 0.04712778545695382}, '2e3m': {'holo': '2e3q', 'rmsd': 0.8849429260246913, 'ligand': '18C', 'overlap': 0.09031491180306017}, '3k0n': {'holo': '1w8m', 'rmsd': 0.39833839324838677, 'ligand': 'E1P', 'overlap': 0.06688885295094558}, '3aap': {'holo': '3aar', 'rmsd': 0.3572808619937969, 'ligand': 'ANP', 'overlap': 0.35015656527393507}, '1bk7': {'holo': '1ucd', 'rmsd': 0.5608983876422892, 'ligand': 'U5P', 'overlap': 0.6239173853318093}, '1tgn': {'holo': '1tni', 'rmsd': 1.937557415391448, 'ligand': 'PBN', 'overlap': 0.17747463719381318}, '1ahc': {'holo': '1ahb', 'rmsd': 0.19271320420696525, 'ligand': 'FMP', 'overlap': 0.2906115852531481}, '3mft': {'holo': '3mfu', 'rmsd': 0.5294319317376655, 'ligand': 'ANP', 'overlap': 0.5901312084425671}, '1erk': {'holo': '4erk', 'rmsd': 1.8403643776383745, 'ligand': 'OLO', 'overlap': 0.2606792215029042}, '1sll': {'holo': '4sli', 'rmsd': 0.27110444955664886, 'ligand': 'CNP', 'overlap': 0.09295343747516871}, '1rtc': {'holo': '1br6', 'rmsd': 0.7736294829905356, 'ligand': 'PT1', 'overlap': 0.02244838892564617}, '1gy0': {'holo': '1og1', 'rmsd': 0.26938443289620695, 'ligand': 'TAD', 'overlap': 0.11337131846497911}, '2hbj': {'holo': '2hbl', 'rmsd': 0.2791711771177489, 'ligand': 'AMP', 'overlap': 0.2654241484172416}, '1rd6': {'holo': '1x6n', 'rmsd': 0.29783872255245375, 'ligand': 'AO3', 'overlap': 0.020924868377188755}, '1mzl': {'holo': '1mzm', 'rmsd': 0.4977604676331388, 'ligand': 'PLM', 'overlap': 0.15358351922009683}, '1hka': {'holo': '1eqm', 'rmsd': 3.101848541109514, 'ligand': '', 'overlap': 0.29938268367048587}, '2ppn': {'holo': '1fkh', 'rmsd': 0.5498751329287473, 'ligand': 'SBX', 'overlap': 0.008728602838910722}, '2vfy': {'holo': '2vfk', 'rmsd': 0.5234309202023915, 'ligand': 'AMP', 'overlap': 0.4186721907995304}, '1yhv': {'holo': '2hy8', 'rmsd': 0.7580522933908197, 'ligand': '1ST', 'overlap': 0.6167885584722257}, '1eur': {'holo': '1eus', 'rmsd': 0.3471466880748409, 'ligand': 'DAN', 'overlap': 0.656212298345592}, '3a0x': {'holo': '3a0t', 'rmsd': 1.4778758006848136, 'ligand': '', 'overlap': 0.08855830101729956}, '1mri': {'holo': '1mrh', 'rmsd': 0.20580455271425652, 'ligand': 'FMC', 'overlap': 0.2616567944535416}, '2zj8': {'holo': '2zja', 'rmsd': 0.39958904846327736, 'ligand': 'ACP', 'overlap': 0.15477640189079234}, '2ggo': {'holo': '2ggq', 'rmsd': 0.41217362536688407, 'ligand': 'TTP', 'overlap': 0.21415433688009522}, '1bbc': {'holo': '1db4', 'rmsd': 1.3942309487855171, 'ligand': '8IN', 'overlap': 0.12052760594053587}, '3c0e': {'holo': '3c11', 'rmsd': 0.40358994049970154, 'ligand': 'GMY', 'overlap': 0.18401963629557413}, '1tje': {'holo': '1tkg', 'rmsd': 0.5953813587886109, 'ligand': 'SSA', 'overlap': 0.3496015526554969}, '1wvw': {'holo': '1wvy', 'rmsd': 0.9674023140460931, 'ligand': 'STU', 'overlap': 0.13588790082354602}, '1ifb': {'holo': '2ifb', 'rmsd': 0.36597921461938127, 'ligand': 'PLM', 'overlap': 0.2728603476085213}, '1qtr': {'holo': '1x2e', 'rmsd': 0.4771243544672125, 'ligand': 'ATX', 'overlap': 0.003613865518161083}, '1pdb': {'holo': '3nxv', 'rmsd': 0.7597621325331363, 'ligand': 'D2F', 'overlap': 0.07405896998672913}, '2qev': {'holo': '2qeh', 'rmsd': 0.26344417969606054, 'ligand': 'SRO', 'overlap': 0.14227481558897756}, '3g6l': {'holo': '3g6m', 'rmsd': 0.1606666839962277, 'ligand': 'CFF', 'overlap': 0.2849120042129244}, '1arl': {'holo': '2rfh', 'rmsd': 0.349832513662932, 'ligand': '23N', 'overlap': 0.314294597059561}, '2vfb': {'holo': '3ltw', 'rmsd': 0.2811833881753945, 'ligand': 'HLZ', 'overlap': 0.03642491522930895}, '1akz': {'holo': '3fck', 'rmsd': 0.5257761484670386, 'ligand': 'FCK', 'overlap': 0.43841384942941936}, '2paw': {'holo': '1a26', 'rmsd': 0.3302247982261714, 'ligand': 'CNA', 'overlap': 0.07651713423342446}, '1kqx': {'holo': '1kqw', 'rmsd': 0.6648748206255346, 'ligand': 'RTL', 'overlap': 0.2243499171280635}, '3pte': {'holo': '1yqs', 'rmsd': 0.19520439729404254, 'ligand': 'BSA', 'overlap': 0.23347037537804616}, '1iad': {'holo': '1qji', 'rmsd': 0.33860211897807524, 'ligand': 'PKF', 'overlap': 0.19714803737138323}, '3blm': {'holo': '1blh', 'rmsd': 0.21196668540983332, 'ligand': 'FOS', 'overlap': 0.35086757294487353}, '2d59': {'holo': '2d5a', 'rmsd': 1.191748760948045, 'ligand': '', 'overlap': 0.19626132868398616}, '1lmn': {'holo': '1bb7', 'rmsd': 0.12717943435514337, 'ligand': 'GUM', 'overlap': 0.1248915962046876}, '1kf5': {'holo': '1eow', 'rmsd': 0.18795641863024692, 'ligand': 'U2G', 'overlap': 0.24946252748542863}, '2sil': {'holo': '2sim', 'rmsd': 0.1425616125054784, 'ligand': 'DAN', 'overlap': 0.021600224935827056}, '2zco': {'holo': '2zcs', 'rmsd': 0.2685379543685112, 'ligand': 'B70', 'overlap': 0.6230371569077154}, '2ac4': {'holo': '2q2o', 'rmsd': 0.4309651624839198, 'ligand': 'H01', 'overlap': 0.02673223321151024}, '1ojq': {'holo': '1ojz', 'rmsd': 0.8206121849480196, 'ligand': '', 'overlap': 0.496036330611861}, '1ey0': {'holo': '1stg', 'rmsd': 0.701380747069386, 'ligand': 'THP', 'overlap': 0.16741304780628186}, '1yes': {'holo': '1byq', 'rmsd': 0.8384813140092745, 'ligand': '', 'overlap': 0.23598536851862675}, '1sye': {'holo': '1syf', 'rmsd': 0.699830584716955, 'ligand': 'THP', 'overlap': 0.03728842503269433}, '1tqo': {'holo': '1tr5', 'rmsd': 0.6441666663712677, 'ligand': 'THP', 'overlap': 0.09943304760719905}, '2gg4': {'holo': '2gg6', 'rmsd': 4.053610489733732, 'ligand': 'S3P', 'overlap': 0.32226681990743533}, '1xqz': {'holo': '1xr1', 'rmsd': 0.8735291349884569, 'ligand': 'ANP', 'overlap': 0.03322565994229599}, '1ri5': {'holo': '1ri1', 'rmsd': 0.4756848438429417, 'ligand': 'GTG', 'overlap': 0.40634410780764774}, '1p38': {'holo': '1bmk', 'rmsd': 0.49475592255264006, 'ligand': 'SB5', 'overlap': 0.16284355281518356}, '1z1i': {'holo': '2gx4', 'rmsd': 0.943572422028664, 'ligand': 'NOL', 'overlap': 0.4448793087833217}, '4ape': {'holo': '2v00', 'rmsd': 0.37870807087859537, 'ligand': 'V15', 'overlap': 0.18888234876200738}, '1wos': {'holo': '1wop', 'rmsd': 0.2366436759942632, 'ligand': 'FFO', 'overlap': 0.051811499365777516}, '1jam': {'holo': '1lp4', 'rmsd': 0.5294940012944506, 'ligand': 'ANP', 'overlap': 0.14077385733755635}, '1mtz': {'holo': '1mu0', 'rmsd': 0.5811436641707657, 'ligand': 'PHK', 'overlap': 0.07460944026829339}, '1jcf': {'holo': '1jcg', 'rmsd': 0.3839860003692609, 'ligand': 'ANP', 'overlap': 0.15740954956225073}, '1xqo': {'holo': '1xqp', 'rmsd': 0.7306725684611958, 'ligand': '8HG', 'overlap': 0.006124059539787701}, '3ewq': {'holo': '3ewr', 'rmsd': 0.4118729656576964, 'ligand': 'APR', 'overlap': 0.2245799191683995}, '2sga': {'holo': '1sgc', 'rmsd': 0.10841517059892908, 'ligand': 'CST', 'overlap': 0.25558985372325815}, '1f10': {'holo': '1n4f', 'rmsd': 0.5191688204021571, 'ligand': 'ASR', 'overlap': 0.6209751871169322}
        }

    for pdb_apo in d_apo2holo.keys():

        bool_continue = skip_pdb(pdb_apo,d_apo2holo,d_cath,)
        if bool_continue == True:
            continue

        pdb_holo = d_apo2holo[pdb_apo]['holo']

##        continue ## tmp!!!
        print pdb_apo, pdb_holo

        ##
        ## parse coordinates
        ##
        d_mmCIF_apo, l_coords_alpha_apo = parse_coords(pdb_apo)
        d_mmCIF_holo, l_coords_alpha_holo = parse_coords(pdb_holo)

        tv1, rm, tv2, l_coords_alpha_apo, l_coords_alpha_holo = get_transformation_matrix(
            d_mmCIF_apo, l_coords_alpha_apo,
            d_mmCIF_holo, l_coords_alpha_holo,
            )

        vector_apo2holo = get_apo_holo_vector(
            d_mmCIF_apo, l_coords_alpha_apo,
            d_mmCIF_holo, l_coords_alpha_holo,
            tv1, rm, tv2,
            )

        chain_apo = ''.join(d_mmCIF_apo['_entity_poly.pdbx_strand_id'])
        chain_holo = ''.join(d_mmCIF_holo['_entity_poly.pdbx_strand_id'])

        ligand_pos_apo, ligand_pos_holo = get_ligand_pos(
            d_mmCIF_holo,
            tv1, rm, tv2,
            d_apo2holo[pdb_apo]['ligand'],
            )

        dist_max = 6
        dist_min = 3
        for pdb, chain, l_coords_alpha, ligand_pos in [
            [pdb_holo,chain_holo,l_coords_alpha_holo,ligand_pos_holo,],
            [pdb_apo,chain_apo,l_coords_alpha_apo,ligand_pos_apo,],
            ]:
            mode_max_apoholo, overlap_max_apoholo, l_overlaps, max_mode = goodvibes_ligand.main(
                pdb,chain,
                dist_max,dist_min,
                v_apoholo=vector_apo2holo,
                l_coords_probe = [ligand_pos],
                l_coords_protein_alpha = l_coords_alpha,
                )
            print pdb, mode_max_apoholo, overlap_max_apoholo, l_overlaps[0]
            s = '%s %s %s %s %s %s\n' %(
                pdb,pdb_apo,pdb_holo,mode_max_apoholo,
                round(overlap_max_apoholo,3), round(l_overlaps[0],3),
                )
            switch_max = 3
            if mode_max_apoholo < 12: suffix_mode = 'low'
            else: suffix_mode = 'high'
            if mode_max_apoholo == max_mode: suffix_switch = 'same'
            else: suffix_switch = 'switched'
            if pdb == pdb_apo: suffix_pdb = 'apo'
            else: suffix_pdb = 'holo'
##            fd = open('induced_or_sampling_%s_%s_%s.txt' %(
##                switch_max,suffix_pdb,suffix_mode,
            fd = open('induced_or_sampling_%s_%s_%s_%s.txt' %(
                switch_max,suffix_pdb,suffix_mode,suffix_switch,
                ),'a')
            fd.write(s)
            fd.close()

    print 'f(x) = x'
    print 'g(x) = x-0.05'
    print 'h(x) = x+0.05'
    print 'set xlabel "max overlap between apo/holo motion"'
    print 'set ylabel "max overlap between apo/holo motion - perturbed"'
    s = 'plot [0.1:1][0.1:1]'
    s += '"induced_or_sampling_3_apo_low_same.txt" u 5:6 t "apo, mode of max contribution = [6:12], mode not switched", '
    s += '"induced_or_sampling_3_holo_low_same.txt" u 5:6 t "holo, mode of max contribution = [6:12], mode not switched", '
    s += '"induced_or_sampling_3_apo_low_switched.txt" u 5:6 t "apo, mode of max contribution = [6:12], mode switched", '
    s += '"induced_or_sampling_3_holo_low_switched.txt" u 5:6 t "holo, mode of max contribution = [6:12], mode switched", '
    s += '"induced_or_sampling_3_apo_high_same.txt" u 5:6 t "apo, mode of max contribution > 12, mode not switched", '
    s += '"induced_or_sampling_3_holo_high_same.txt" u 5:6 t "holo, mode of max contribution > 12, mode not switched", '
    s += '"induced_or_sampling_3_apo_high_switched.txt" u 5:6 t "apo, mode of max contribution > 12, mode switched", '
    s += '"induced_or_sampling_3_holo_high_switched.txt" u 5:6 t "holo, mode of max contribution > 12, mode switched", '
    s += 'f(x) t "", g(x) t "" lc 0, h(x) t "" lc 0'
    print s
    print 'set key out vert top right'
    print 'plot [0.1:1][0.1:1]"induced_or_sampling_3_apo_low.txt" u 5:6 t "apo, mode of max contribution = [6:12]", "induced_or_sampling_3_holo_low.txt" u 5:6 t "holo, mode of max contribution = [6:12]", "induced_or_sampling_3_apo_high.txt" u 5:6 t "apo, mode of max contribution > 12", "induced_or_sampling_3_holo_high.txt" u 5:6 t "holo, mode of max overlap > 12", f(x) t "", g(x) t "" lc 0, h(x) t "" lc 0'
    print 'end of main'
    
    return


def get_ligand_pos(
    d_mmCIF_holo,
    tv1, rm, tv2,
    ligand,
    ):

    l_coords_ligand = []
    l_coords_ligand_apo = []
    for i in range(len(d_mmCIF_holo['_atom_site.id'])):
        if (
            d_mmCIF_holo['_atom_site.group_PDB'][i] == 'HETATM'
            and
            d_mmCIF_holo['_atom_site.label_comp_id'][i] == ligand
            ):
            x = float(d_mmCIF_holo['_atom_site.Cartn_x'][i])
            y = float(d_mmCIF_holo['_atom_site.Cartn_y'][i])
            z = float(d_mmCIF_holo['_atom_site.Cartn_z'][i])
            coord = numpy.array([x,y,z,])
            l_coords_ligand += [coord]
            ## holo2apo
            coord_apo = numpy.dot(coord-tv1,rm)+tv2
##            coord_apo = numpy.dot(coord-tv2,rm_transpose)+tv1
            l_coords_ligand_apo += [coord_apo]

    position_ligand_apo = sum(l_coords_ligand_apo)/len(l_coords_ligand_apo)
    position_ligand = sum(l_coords_ligand)/len(l_coords_ligand)

    return position_ligand_apo, position_ligand


def parse_CATH():

    fd = open('datasets/CathDomall.v3.4.0','r')
    lines = fd.readlines()
    fd.close()
    d_cath = {}
    for line in lines:
        if line[0] == '#':
            continue
        pdb = line[:4]
        count_domains = int(line[7:9])
        d_cath[pdb] = count_domains

    return d_cath


def skip_pdb(pdb_apo,d_apo2holo,d_cath,):

    bool_continue = False
    if d_apo2holo[pdb_apo]['ligand'] == '':
        return True
    pdb_holo = d_apo2holo[pdb_apo]['holo']
    if pdb_holo in ['1ucd','1ri1','3nxv',]: ## 2 ligands
        return True
    ## exclude multi domain proteins
    if not (pdb_apo in d_cath.keys() and pdb_holo in d_cath.keys()):
        if pdb_holo in ['2zja','3aar','2hbl','2ggq',]:
            return True
        elif pdb_holo in [
            '3nq3','2e3q','1w8m','1fkh','2vfk','2qeh','3ewr',
            '2zcs', ## 1 or 2 domains...
            ]:
            pass
        else:
            print pdb_apo, pdb_holo
            stop
    else:
        if d_cath[pdb_apo] > 1 and d_cath[pdb_holo] > 1:
            return True
        if not(d_cath[pdb_apo] == 1 and d_cath[pdb_holo] == 1):
            print d_cath[pdb_apo], d_cath[pdb_holo]
            stop

    return bool_continue

def get_transformation_matrix(
    d_mmCIF_apo, l_coords_alpha_apo,
    d_mmCIF_holo, l_coords_alpha_holo,
    ):

    ## structural alignment
    ## solution that works in all cases
    ## also for 2d59 and 2d5a, which have residues missing at the Nterm and Cterm, respectively
    ## first non-?
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

    instance_geometry = geometry.geometry()
    rmsd = instance_geometry.superpose(l_coords_alpha_apo,l_coords_alpha_holo,)
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

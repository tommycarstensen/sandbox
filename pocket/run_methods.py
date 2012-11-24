import os

def main():

    set_pdbs_apo = ['1ak1', '1lp8', '2exo', '1syc', '1cua', '1znw', '1gbs', '1w8v', '3a0x', '1sqg', '3npo', '2e3m', '3k0n', '3aap', '1bk7', '1tgn', '1ahc', '3mft', '1erk', '1sll', '1rtc', '1gy0', '2hbj', '1rd6', '1mzl', '1hka', '2ppn', '2vfy', '1yhv', '1eyd', '1eur', '1bbc', '1mri', '2zj8', '2ggo', '2vfb', '1f10', '3c0e', '1tje', '1wvw', '1wos', '1qtr', '1pdb', '2qev', '3g6l', '1arl', '1akz', '2paw', '1kqx', '3pte', '1iad', '2d59', '1vds', '1lmn', '1kf5', '1jcf', '2zco', '2ac4', '1ojq', '1ey0', '1yes', '1sye', '1tqo', '2gg4', '1xqz', '1ri5', '1p38', '1z1i', '4ape', '1ifb', '1jam', '1mtz', '2sil', '1xqo', '3ewq', '2sga', '3blm']
    set_pdbs_holo = set([])
    d_apo2holo = {'1ak1': {'holo': '1c1h', 'rmsd': 0.6390067852762965, 'ligand': 'MMP', 'overlap': 0.36157936212343905}, '1lp8': {'holo': '1lpc', 'rmsd': 0.14003320157149318, 'ligand': 'CMP', 'overlap': 0.1663339954531419}, '2exo': {'holo': '1j01', 'rmsd': 0.31317202378931264, 'ligand': 'XIL', 'overlap': 0.5464584014950491}, '1syc': {'holo': '1syd', 'rmsd': 0.42370903255794784, 'ligand': 'THP', 'overlap': 0.312528957121447}, '1cua': {'holo': '3esd', 'rmsd': 0.4040397349113579, 'ligand': 'SXC', 'overlap': 0.08139930257927935}, '1znw': {'holo': '1znx', 'rmsd': 0.24704558797447082, 'ligand': '5GP', 'overlap': 0.12436531581580901}, '1gbs': {'holo': '1lsp', 'rmsd': 0.264857015278693, 'ligand': 'BUL', 'overlap': 0.2934470556757598}, '1sqg': {'holo': '1sqf', 'rmsd': 0.6990763011194481, 'ligand': 'SAM', 'overlap': 0.4723483760294988}, '3npo': {'holo': '3nq3', 'rmsd': 0.7183351557421368, 'ligand': 'DKA', 'overlap': 0.04712778545695382}, '2e3m': {'holo': '2e3q', 'rmsd': 0.8849429260246913, 'ligand': '18C', 'overlap': 0.09031491180306017}, '3k0n': {'holo': '1w8m', 'rmsd': 0.39833839324838677, 'ligand': 'E1P', 'overlap': 0.06688885295094558}, '3aap': {'holo': '3aar', 'rmsd': 0.3572808619937969, 'ligand': 'ANP', 'overlap': 0.35015656527393507}, '1bk7': {'holo': '1ucd', 'rmsd': 0.5608983876422892, 'ligand': 'U5P', 'overlap': 0.6239173853318093}, '1tgn': {'holo': '1tni', 'rmsd': 1.937557415391448, 'ligand': 'PBN', 'overlap': 0.17747463719381318}, '1ahc': {'holo': '1ahb', 'rmsd': 0.19271320420696525, 'ligand': 'FMP', 'overlap': 0.2906115852531481}, '3mft': {'holo': '3mfu', 'rmsd': 0.5294319317376655, 'ligand': 'ANP', 'overlap': 0.5901312084425671}, '1erk': {'holo': '4erk', 'rmsd': 1.8403643776383745, 'ligand': 'OLO', 'overlap': 0.2606792215029042}, '1sll': {'holo': '4sli', 'rmsd': 0.27110444955664886, 'ligand': 'CNP', 'overlap': 0.09295343747516871}, '1rtc': {'holo': '1br6', 'rmsd': 0.7736294829905356, 'ligand': 'PT1', 'overlap': 0.02244838892564617}, '1gy0': {'holo': '1og1', 'rmsd': 0.26938443289620695, 'ligand': 'TAD', 'overlap': 0.11337131846497911}, '2hbj': {'holo': '2hbl', 'rmsd': 0.2791711771177489, 'ligand': 'AMP', 'overlap': 0.2654241484172416}, '1rd6': {'holo': '1x6n', 'rmsd': 0.29783872255245375, 'ligand': 'AO3', 'overlap': 0.020924868377188755}, '1mzl': {'holo': '1mzm', 'rmsd': 0.4977604676331388, 'ligand': 'PLM', 'overlap': 0.15358351922009683}, '1hka': {'holo': '1eqm', 'rmsd': 3.101848541109514, 'ligand': '', 'overlap': 0.29938268367048587}, '2ppn': {'holo': '1fkh', 'rmsd': 0.5498751329287473, 'ligand': 'SBX', 'overlap': 0.008728602838910722}, '2vfy': {'holo': '2vfk', 'rmsd': 0.5234309202023915, 'ligand': 'AMP', 'overlap': 0.4186721907995304}, '1yhv': {'holo': '2hy8', 'rmsd': 0.7580522933908197, 'ligand': '1ST', 'overlap': 0.6167885584722257}, '1eur': {'holo': '1eus', 'rmsd': 0.3471466880748409, 'ligand': 'DAN', 'overlap': 0.656212298345592}, '3a0x': {'holo': '3a0t', 'rmsd': 1.4778758006848136, 'ligand': '', 'overlap': 0.08855830101729956}, '1mri': {'holo': '1mrh', 'rmsd': 0.20580455271425652, 'ligand': 'FMC', 'overlap': 0.2616567944535416}, '2zj8': {'holo': '2zja', 'rmsd': 0.39958904846327736, 'ligand': 'ACP', 'overlap': 0.15477640189079234}, '2ggo': {'holo': '2ggq', 'rmsd': 0.41217362536688407, 'ligand': 'TTP', 'overlap': 0.21415433688009522}, '1bbc': {'holo': '1db4', 'rmsd': 1.3942309487855171, 'ligand': '8IN', 'overlap': 0.12052760594053587}, '3c0e': {'holo': '3c11', 'rmsd': 0.40358994049970154, 'ligand': 'GMY', 'overlap': 0.18401963629557413}, '1tje': {'holo': '1tkg', 'rmsd': 0.5953813587886109, 'ligand': 'SSA', 'overlap': 0.3496015526554969}, '1wvw': {'holo': '1wvy', 'rmsd': 0.9674023140460931, 'ligand': 'STU', 'overlap': 0.13588790082354602}, '1ifb': {'holo': '2ifb', 'rmsd': 0.36597921461938127, 'ligand': 'PLM', 'overlap': 0.2728603476085213}, '1qtr': {'holo': '1x2e', 'rmsd': 0.4771243544672125, 'ligand': 'ATX', 'overlap': 0.003613865518161083}, '1pdb': {'holo': '3nxv', 'rmsd': 0.7597621325331363, 'ligand': 'D2F', 'overlap': 0.07405896998672913}, '2qev': {'holo': '2qeh', 'rmsd': 0.26344417969606054, 'ligand': 'SRO', 'overlap': 0.14227481558897756}, '3g6l': {'holo': '3g6m', 'rmsd': 0.1606666839962277, 'ligand': 'CFF', 'overlap': 0.2849120042129244}, '1arl': {'holo': '2rfh', 'rmsd': 0.349832513662932, 'ligand': '23N', 'overlap': 0.314294597059561}, '2vfb': {'holo': '3ltw', 'rmsd': 0.2811833881753945, 'ligand': 'HLZ', 'overlap': 0.03642491522930895}, '1akz': {'holo': '3fck', 'rmsd': 0.5257761484670386, 'ligand': 'FCK', 'overlap': 0.43841384942941936}, '2paw': {'holo': '1a26', 'rmsd': 0.3302247982261714, 'ligand': 'CNA', 'overlap': 0.07651713423342446}, '1kqx': {'holo': '1kqw', 'rmsd': 0.6648748206255346, 'ligand': 'RTL', 'overlap': 0.2243499171280635}, '3pte': {'holo': '1yqs', 'rmsd': 0.19520439729404254, 'ligand': 'BSA', 'overlap': 0.23347037537804616}, '1iad': {'holo': '1qji', 'rmsd': 0.33860211897807524, 'ligand': 'PKF', 'overlap': 0.19714803737138323}, '3blm': {'holo': '1blh', 'rmsd': 0.21196668540983332, 'ligand': 'FOS', 'overlap': 0.35086757294487353}, '2d59': {'holo': '2d5a', 'rmsd': 1.191748760948045, 'ligand': '', 'overlap': 0.19626132868398616}, '1lmn': {'holo': '1bb7', 'rmsd': 0.12717943435514337, 'ligand': 'GUM', 'overlap': 0.1248915962046876}, '1kf5': {'holo': '1eow', 'rmsd': 0.18795641863024692, 'ligand': 'U2G', 'overlap': 0.24946252748542863}, '2sil': {'holo': '2sim', 'rmsd': 0.1425616125054784, 'ligand': 'DAN', 'overlap': 0.021600224935827056}, '2zco': {'holo': '2zcs', 'rmsd': 0.2685379543685112, 'ligand': 'B70', 'overlap': 0.6230371569077154}, '2ac4': {'holo': '2q2o', 'rmsd': 0.4309651624839198, 'ligand': 'H01', 'overlap': 0.02673223321151024}, '1ojq': {'holo': '1ojz', 'rmsd': 0.8206121849480196, 'ligand': '', 'overlap': 0.496036330611861}, '1ey0': {'holo': '1stg', 'rmsd': 0.701380747069386, 'ligand': 'THP', 'overlap': 0.16741304780628186}, '1yes': {'holo': '1byq', 'rmsd': 0.8384813140092745, 'ligand': '', 'overlap': 0.23598536851862675}, '1sye': {'holo': '1syf', 'rmsd': 0.699830584716955, 'ligand': 'THP', 'overlap': 0.03728842503269433}, '1tqo': {'holo': '1tr5', 'rmsd': 0.6441666663712677, 'ligand': 'THP', 'overlap': 0.09943304760719905}, '2gg4': {'holo': '2gg6', 'rmsd': 4.053610489733732, 'ligand': 'S3P', 'overlap': 0.32226681990743533}, '1xqz': {'holo': '1xr1', 'rmsd': 0.8735291349884569, 'ligand': 'ANP', 'overlap': 0.03322565994229599}, '1ri5': {'holo': '1ri1', 'rmsd': 0.4756848438429417, 'ligand': 'GTG', 'overlap': 0.40634410780764774}, '1p38': {'holo': '1bmk', 'rmsd': 0.49475592255264006, 'ligand': 'SB5', 'overlap': 0.16284355281518356}, '1z1i': {'holo': '2gx4', 'rmsd': 0.943572422028664, 'ligand': 'NOL', 'overlap': 0.4448793087833217}, '4ape': {'holo': '2v00', 'rmsd': 0.37870807087859537, 'ligand': 'V15', 'overlap': 0.18888234876200738}, '1wos': {'holo': '1wop', 'rmsd': 0.2366436759942632, 'ligand': 'FFO', 'overlap': 0.051811499365777516}, '1jam': {'holo': '1lp4', 'rmsd': 0.5294940012944506, 'ligand': 'ANP', 'overlap': 0.14077385733755635}, '1mtz': {'holo': '1mu0', 'rmsd': 0.5811436641707657, 'ligand': 'PHK', 'overlap': 0.07460944026829339}, '1jcf': {'holo': '1jcg', 'rmsd': 0.3839860003692609, 'ligand': 'ANP', 'overlap': 0.15740954956225073}, '1xqo': {'holo': '1xqp', 'rmsd': 0.7306725684611958, 'ligand': '8HG', 'overlap': 0.006124059539787701}, '3ewq': {'holo': '3ewr', 'rmsd': 0.4118729656576964, 'ligand': 'APR', 'overlap': 0.2245799191683995}, '2sga': {'holo': '1sgc', 'rmsd': 0.10841517059892908, 'ligand': 'CST', 'overlap': 0.25558985372325815}, '1f10': {'holo': '1n4f', 'rmsd': 0.5191688204021571, 'ligand': 'ASR', 'overlap': 0.6209751871169322}}

    for pdb_apo in set_pdbs_apo:

        if not pdb_apo in d_apo2holo.keys():
            continue

        if d_apo2holo[pdb_apo]['ligand'] == '':
            continue

        pdb_holo = d_apo2holo[pdb_apo]['holo']

        pdb = pdb_apo
        pdb = pdb_holo

    for pdb in [
        ## hinge
        '135lA','1og1A','1a2t','1aj0A','1bqcA','1ra2A',
        ## potential hinge
        '1esoA', '1ak0A','1tmlA','1smlA','1akoA','1lbaA','1bolA','1l9xA',
        '1pjaA','1n29A','1fobA','1cdeA','1ga8A',
        ## multimer, but monomer hinge
        '1jh6A','1vzxA','2rnfA',

        '2lipA','1xqwA','1bu7A','1bzcA','1i78A','1qe3A','1d3gA','1rddA',
        '2aceA','1gojA',
        ## probably not hinge
        '1p1xA','1ptdA','1j00A','1u5uA',
        ## multimer and no hinges
        '1a8qA','1a4lA',

        ## unknown
          '2cpoA', '1ru4A', '1ljlA', '1mrqA', '1w2nA', '1cz1A',
        '1mugA', '1o8aA', '1rtuA', '1uchA', '1b6gA', '1chdA', 
        '1oxaA', '1r44A', '1w0hA', '1expA', '1bt1A', '1i9aA', '1im5A', '1pp4A',
         '2ebnA', '1d1qA', '1ehyA', '1geqA', '1ca2A', '1gnsA', '1eh5A',
        '1l6pA', '1r4zA', '1a2tA', '1di1A', '2ayhA', '1astA', '1cm0A', '1cv2A',
         '1dveA',   '1un1A', '1btlA', '2pecA', '1fcqA',
        '1czfA', '1thgA', '1booA', '1iu4A',  '206lA',  '1snzA',
        '1gq8A', '1aqlA', '1ps1A', '1s95A', '1pylA',  '1b6bA', '1pntA',
        '1e1aA', '2f9rA', '1v04A', '2nlrA', '1pbgA', '5cpaA', '1agmA',
        '1byaA', '1r76A', '1vidA', '1h4gA', '1akdA', '1fy2A', '1xqdA',
        '1d6oA', '1qv0A', '1qjeA', '1fvaA', '1bp2A', '1ah7A', '2pthA', '2engA',
        '2acyA', '1qazA', '2a0nA', '1dl2A', '1gp5A', '1onrA', '1cwyA', '1pudA',
        '1bs9A', '1dinA', '1xyzA', '1bwlA', '1eugA', '1idjA', '1g24A', '1oygA',
        '1hzfA', '9papA', '1eb6A', '1ghsA', '1rbnA', '1bixA', '1bs4A', '1celA',
        '1hkaA', '1b02A', '1qibA', '1u3fA', '1agyA', '1zioA', '1pa9A', '2tpsA',
        '2plcA', '1qk2A', '1j53A', '1m21A',
        ]:
        pdb = pdb[:4]

##        ## tmp!!! ################################
##        if not os.path.isfile('/media/Tommy/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb,)):
##            print pdb
##            import urllib2
##            url = urllib2.urlopen('http://www.pdb.org/pdb/files/%s.pdb' %(pdb))
##            lines = url.readlines()
##            if not os.path.isdir('/media/Tommy/pdb/%s' %(pdb[1:3])):
##                os.mkdir('/media/Tommy/pdb/%s' %(pdb[1:3]))
##            fd = open('/media/Tommy/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb,), 'w')
##            fd.writelines(lines)
##            fd.close()
##
##        continue

        os.system('cp /media/WDMyBook1TB/2TB/pdb/%s/pdb%s.ent %s.pdb' %(pdb[1:3],pdb,pdb,))
    
        chain = get_first_chain(pdb)

        print pdb, chain

        if not os.path.isdir('output/Fpocket/%s_out' %(pdb)):
            run_Fpocket(pdb)

        if not os.path.isfile('output/Concavity/%s_concavity_pocket.pdb' %(pdb)):
            run_Concavity(pdb,chain,)

        if not os.path.isfile('output/POCASA/%s_Pocket_DepthCenters.pdb' %(pdb)):
            run_POCASA(pdb,chain,)

        if not os.path.isfile('output/GoodVibes/distmax6_distmin3/%s_%s_probe.pdb' %(pdb,chain,)):
            run_GoodVibes(pdb,chain,)

        if os.path.isfile('%s.pdb' %(pdb)):
            os.remove('%s.pdb' %(pdb))

    return


def get_first_chain(pdb):

    fd = open('%s.pdb' %(pdb,),'r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            chain = line[21]
            break

    return chain


def run_GoodVibes(pdb,chain,):

    dist_max = 6
    dist_min = 3
    s = 'python /home/tc/svn/GoodVibes/goodvibes_ligand.py --pdb %s --chain %s --dist_max %s --dist_min %s' %(
        pdb,chain,dist_max,dist_min,
        )
    print s
    os.system(s)
    os.system('mv %s_%s_probe.pdb output/GoodVibes/distmax6_distmin3/.' %(pdb,chain,))
    
    return


def run_Fpocket(pdb):

    s = 'fpocket -f %s.pdb' %(pdb,)
    print s
    os.system(s)
    os.system('mv %s_out output/Fpocket/.' %(pdb))

    return


def run_Concavity(pdb,chain,):

    s = 'concavity -print_grid_pdb 1 %s.pdb concavity' %(pdb,)
    print s
    os.system(s)
    os.system('mv %s_concavity_pocket.pdb output/Concavity/.' %(pdb,))
    os.remove('%s_concavity.dx' %(pdb,))
    os.remove('%s_%s_concavity.scores' %(pdb,chain,))
    os.remove('%s_concavity.pml' %(pdb,))
    os.remove('%s_concavity_residue.pdb' %(pdb,))

    return


def run_POCASA(pdb,chain,):

    probe_radius = 2
    if pdb == '1jcf':
        probe_radius = 3
    s = 'POCASA %s.pdb %s 1 16 18 5 %s' %(pdb,probe_radius,chain,)
    print s
    os.system(s)

    os.system('mv %s_Pocket_DepthCenters.pdb output/POCASA/%s_Pocket_DepthCenters.pdb' %(pdb,pdb,))
    os.remove('%s_simple.pdb' %(pdb,))
    if os.path.isfile('%s_TopN_pockets.pdb' %(pdb,)):
        os.remove('%s_TopN_pockets.pdb' %(pdb,))
    os.remove('%s_Parameters.txt' %(pdb,))

    return


if __name__ == '__main__':
    main()

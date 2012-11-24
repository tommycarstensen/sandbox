import numpy, math
import sys
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('/home/tc/svn/Protool/')
import geometry
sys.path.append('/home/tc/svn/GoodVibes/')
import NMA

def main():

    l_pdbs_apo = ['1ak1', '1lp8', '2exo', '1syc', '1cua', '1znw', '1gbs', '1w8v', '3a0x', '1sqg', '3npo', '2e3m', '3k0n', '3aap', '1bk7', '1tgn', '1ahc', '3mft', '1erk', '1sll', '1rtc', '1gy0', '2hbj', '1rd6', '1mzl', '1hka', '2ppn', '2vfy', '1yhv', '1eyd', '1eur', '1bbc', '1mri', '2zj8', '2ggo', '2vfb', '1f10', '3c0e', '1tje', '1wvw', '1wos', '1qtr', '1pdb', '2qev', '3g6l', '1arl', '1akz', '2paw', '1kqx', '3pte', '1iad', '2d59', '1vds', '1lmn', '1kf5', '1jcf', '2zco', '2ac4', '1ojq', '1ey0', '1yes', '1sye', '1tqo', '2gg4', '1xqz', '1ri5', '1p38', '1z1i', '4ape', '1ifb', '1jam', '1mtz', '2sil', '1xqo', '3ewq', '2sga', '3blm']
    d_apo2holo = {'1ak1': {'holo': '1c1h', 'rmsd': 0.6390067852762965, 'ligand': 'MMP', 'overlap': 0.36157936212343905}, '1lp8': {'holo': '1lpc', 'rmsd': 0.14003320157149318, 'ligand': 'CMP', 'overlap': 0.1663339954531419}, '2exo': {'holo': '1j01', 'rmsd': 0.31317202378931264, 'ligand': 'XIL', 'overlap': 0.5464584014950491}, '1syc': {'holo': '1syd', 'rmsd': 0.42370903255794784, 'ligand': 'THP', 'overlap': 0.312528957121447}, '1cua': {'holo': '3esd', 'rmsd': 0.4040397349113579, 'ligand': 'SXC', 'overlap': 0.08139930257927935}, '1znw': {'holo': '1znx', 'rmsd': 0.24704558797447082, 'ligand': '5GP', 'overlap': 0.12436531581580901}, '1gbs': {'holo': '1lsp', 'rmsd': 0.264857015278693, 'ligand': 'BUL', 'overlap': 0.2934470556757598}, '1sqg': {'holo': '1sqf', 'rmsd': 0.6990763011194481, 'ligand': 'SAM', 'overlap': 0.4723483760294988}, '3npo': {'holo': '3nq3', 'rmsd': 0.7183351557421368, 'ligand': 'DKA', 'overlap': 0.04712778545695382}, '2e3m': {'holo': '2e3q', 'rmsd': 0.8849429260246913, 'ligand': '18C', 'overlap': 0.09031491180306017}, '3k0n': {'holo': '1w8m', 'rmsd': 0.39833839324838677, 'ligand': 'E1P', 'overlap': 0.06688885295094558}, '3aap': {'holo': '3aar', 'rmsd': 0.3572808619937969, 'ligand': 'ANP', 'overlap': 0.35015656527393507}, '1bk7': {'holo': '1ucd', 'rmsd': 0.5608983876422892, 'ligand': 'U5P', 'overlap': 0.6239173853318093}, '1tgn': {'holo': '1tni', 'rmsd': 1.937557415391448, 'ligand': 'PBN', 'overlap': 0.17747463719381318}, '1ahc': {'holo': '1ahb', 'rmsd': 0.19271320420696525, 'ligand': 'FMP', 'overlap': 0.2906115852531481}, '3mft': {'holo': '3mfu', 'rmsd': 0.5294319317376655, 'ligand': 'ANP', 'overlap': 0.5901312084425671}, '1erk': {'holo': '4erk', 'rmsd': 1.8403643776383745, 'ligand': 'OLO', 'overlap': 0.2606792215029042}, '1sll': {'holo': '4sli', 'rmsd': 0.27110444955664886, 'ligand': 'CNP', 'overlap': 0.09295343747516871}, '1rtc': {'holo': '1br6', 'rmsd': 0.7736294829905356, 'ligand': 'PT1', 'overlap': 0.02244838892564617}, '1gy0': {'holo': '1og1', 'rmsd': 0.26938443289620695, 'ligand': 'TAD', 'overlap': 0.11337131846497911}, '2hbj': {'holo': '2hbl', 'rmsd': 0.2791711771177489, 'ligand': 'AMP', 'overlap': 0.2654241484172416}, '1rd6': {'holo': '1x6n', 'rmsd': 0.29783872255245375, 'ligand': 'AO3', 'overlap': 0.020924868377188755}, '1mzl': {'holo': '1mzm', 'rmsd': 0.4977604676331388, 'ligand': 'PLM', 'overlap': 0.15358351922009683}, '1hka': {'holo': '1eqm', 'rmsd': 3.101848541109514, 'ligand': '', 'overlap': 0.29938268367048587}, '2ppn': {'holo': '1fkh', 'rmsd': 0.5498751329287473, 'ligand': 'SBX', 'overlap': 0.008728602838910722}, '2vfy': {'holo': '2vfk', 'rmsd': 0.5234309202023915, 'ligand': 'AMP', 'overlap': 0.4186721907995304}, '1yhv': {'holo': '2hy8', 'rmsd': 0.7580522933908197, 'ligand': '1ST', 'overlap': 0.6167885584722257}, '1eur': {'holo': '1eus', 'rmsd': 0.3471466880748409, 'ligand': 'DAN', 'overlap': 0.656212298345592}, '3a0x': {'holo': '3a0t', 'rmsd': 1.4778758006848136, 'ligand': '', 'overlap': 0.08855830101729956}, '1mri': {'holo': '1mrh', 'rmsd': 0.20580455271425652, 'ligand': 'FMC', 'overlap': 0.2616567944535416}, '2zj8': {'holo': '2zja', 'rmsd': 0.39958904846327736, 'ligand': 'ACP', 'overlap': 0.15477640189079234}, '2ggo': {'holo': '2ggq', 'rmsd': 0.41217362536688407, 'ligand': 'TTP', 'overlap': 0.21415433688009522}, '1bbc': {'holo': '1db4', 'rmsd': 1.3942309487855171, 'ligand': '8IN', 'overlap': 0.12052760594053587}, '3c0e': {'holo': '3c11', 'rmsd': 0.40358994049970154, 'ligand': 'GMY', 'overlap': 0.18401963629557413}, '1tje': {'holo': '1tkg', 'rmsd': 0.5953813587886109, 'ligand': 'SSA', 'overlap': 0.3496015526554969}, '1wvw': {'holo': '1wvy', 'rmsd': 0.9674023140460931, 'ligand': 'STU', 'overlap': 0.13588790082354602}, '1ifb': {'holo': '2ifb', 'rmsd': 0.36597921461938127, 'ligand': 'PLM', 'overlap': 0.2728603476085213}, '1qtr': {'holo': '1x2e', 'rmsd': 0.4771243544672125, 'ligand': 'ATX', 'overlap': 0.003613865518161083}, '1pdb': {'holo': '3nxv', 'rmsd': 0.7597621325331363, 'ligand': 'D2F', 'overlap': 0.07405896998672913}, '2qev': {'holo': '2qeh', 'rmsd': 0.26344417969606054, 'ligand': 'SRO', 'overlap': 0.14227481558897756}, '3g6l': {'holo': '3g6m', 'rmsd': 0.1606666839962277, 'ligand': 'CFF', 'overlap': 0.2849120042129244}, '1arl': {'holo': '2rfh', 'rmsd': 0.349832513662932, 'ligand': '23N', 'overlap': 0.314294597059561}, '2vfb': {'holo': '3ltw', 'rmsd': 0.2811833881753945, 'ligand': 'HLZ', 'overlap': 0.03642491522930895}, '1akz': {'holo': '3fck', 'rmsd': 0.5257761484670386, 'ligand': 'FCK', 'overlap': 0.43841384942941936}, '2paw': {'holo': '1a26', 'rmsd': 0.3302247982261714, 'ligand': 'CNA', 'overlap': 0.07651713423342446}, '1kqx': {'holo': '1kqw', 'rmsd': 0.6648748206255346, 'ligand': 'RTL', 'overlap': 0.2243499171280635}, '3pte': {'holo': '1yqs', 'rmsd': 0.19520439729404254, 'ligand': 'BSA', 'overlap': 0.23347037537804616}, '1iad': {'holo': '1qji', 'rmsd': 0.33860211897807524, 'ligand': 'PKF', 'overlap': 0.19714803737138323}, '3blm': {'holo': '1blh', 'rmsd': 0.21196668540983332, 'ligand': 'FOS', 'overlap': 0.35086757294487353}, '2d59': {'holo': '2d5a', 'rmsd': 1.191748760948045, 'ligand': '', 'overlap': 0.19626132868398616}, '1lmn': {'holo': '1bb7', 'rmsd': 0.12717943435514337, 'ligand': 'GUM', 'overlap': 0.1248915962046876}, '1kf5': {'holo': '1eow', 'rmsd': 0.18795641863024692, 'ligand': 'U2G', 'overlap': 0.24946252748542863}, '2sil': {'holo': '2sim', 'rmsd': 0.1425616125054784, 'ligand': 'DAN', 'overlap': 0.021600224935827056}, '2zco': {'holo': '2zcs', 'rmsd': 0.2685379543685112, 'ligand': 'B70', 'overlap': 0.6230371569077154}, '2ac4': {'holo': '2q2o', 'rmsd': 0.4309651624839198, 'ligand': 'H01', 'overlap': 0.02673223321151024}, '1ojq': {'holo': '1ojz', 'rmsd': 0.8206121849480196, 'ligand': '', 'overlap': 0.496036330611861}, '1ey0': {'holo': '1stg', 'rmsd': 0.701380747069386, 'ligand': 'THP', 'overlap': 0.16741304780628186}, '1yes': {'holo': '1byq', 'rmsd': 0.8384813140092745, 'ligand': '', 'overlap': 0.23598536851862675}, '1sye': {'holo': '1syf', 'rmsd': 0.699830584716955, 'ligand': 'THP', 'overlap': 0.03728842503269433}, '1tqo': {'holo': '1tr5', 'rmsd': 0.6441666663712677, 'ligand': 'THP', 'overlap': 0.09943304760719905}, '2gg4': {'holo': '2gg6', 'rmsd': 4.053610489733732, 'ligand': 'S3P', 'overlap': 0.32226681990743533}, '1xqz': {'holo': '1xr1', 'rmsd': 0.8735291349884569, 'ligand': 'ANP', 'overlap': 0.03322565994229599}, '1ri5': {'holo': '1ri1', 'rmsd': 0.4756848438429417, 'ligand': 'GTG', 'overlap': 0.40634410780764774}, '1p38': {'holo': '1bmk', 'rmsd': 0.49475592255264006, 'ligand': 'SB5', 'overlap': 0.16284355281518356}, '1z1i': {'holo': '2gx4', 'rmsd': 0.943572422028664, 'ligand': 'NOL', 'overlap': 0.4448793087833217}, '4ape': {'holo': '2v00', 'rmsd': 0.37870807087859537, 'ligand': 'V15', 'overlap': 0.18888234876200738}, '1wos': {'holo': '1wop', 'rmsd': 0.2366436759942632, 'ligand': 'FFO', 'overlap': 0.051811499365777516}, '1jam': {'holo': '1lp4', 'rmsd': 0.5294940012944506, 'ligand': 'ANP', 'overlap': 0.14077385733755635}, '1mtz': {'holo': '1mu0', 'rmsd': 0.5811436641707657, 'ligand': 'PHK', 'overlap': 0.07460944026829339}, '1jcf': {'holo': '1jcg', 'rmsd': 0.3839860003692609, 'ligand': 'ANP', 'overlap': 0.15740954956225073}, '1xqo': {'holo': '1xqp', 'rmsd': 0.7306725684611958, 'ligand': '8HG', 'overlap': 0.006124059539787701}, '3ewq': {'holo': '3ewr', 'rmsd': 0.4118729656576964, 'ligand': 'APR', 'overlap': 0.2245799191683995}, '2sga': {'holo': '1sgc', 'rmsd': 0.10841517059892908, 'ligand': 'CST', 'overlap': 0.25558985372325815}, '1f10': {'holo': '1n4f', 'rmsd': 0.5191688204021571, 'ligand': 'ASR', 'overlap': 0.6209751871169322}}

    d_distances = {'Fpocket':{},'POCASA':{},'Concavity':{},'null':{},'GoodVibes':{},}
    d_differences = {'Fpocket':{},'POCASA':{},'Concavity':{},'null':{},'GoodVibes':{},}

    ## two different ligands
##    l_pdbs_apo.remove('2exo')
##    l_pdbs_apo.remove('1sqg')
##    l_pdbs_apo.remove('1eur')
##    l_pdbs_apo.remove('1akz')
##    l_pdbs_apo.remove('1z1i')
    ## two or more identical ligands
    l_pdbs_apo.remove('2vfb') ## 3ltw, HLZ
    l_pdbs_apo.remove('1f10') ## 1n4f, ASR
    l_pdbs_apo.remove('1wos') ## 1wop, FFO
    l_pdbs_apo.remove('3g6l') ## 3g6m, CFF
    l_pdbs_apo.remove('3pte') ## 1yqs, BSA
    l_pdbs_apo.remove('1xqo') ## 1xqp, 8HG

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

    l_pdbs_holo = []
    for pdb_apo in list(l_pdbs_apo):
        if not pdb_apo in d_apo2holo.keys():
            l_pdbs_apo.remove(pdb_apo)
            continue
        if d_apo2holo[pdb_apo]['ligand'] == '':
            l_pdbs_apo.remove(pdb_apo)
            continue
        pdb_holo = d_apo2holo[pdb_apo]['holo']
        if pdb_holo in ['1ucd','1ri1','3nxv',]: ## 2 ligands
            l_pdbs_apo.remove(pdb_apo)
            continue
        ## exclude multi domain proteins
        if not (pdb_apo in d_cath.keys() and pdb_holo in d_cath.keys()):
            if pdb_holo in ['2zja','3aar','2hbl','2ggq',]:
                l_pdbs_apo.remove(pdb_apo)
                continue
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
                l_pdbs_apo.remove(pdb_apo)
                continue
            if not(d_cath[pdb_apo] == 1 and d_cath[pdb_holo] == 1):
                print d_cath[pdb_apo], d_cath[pdb_holo]
                stop
        l_pdbs_holo += [pdb_holo]

    ##
    ## loop
    ##
    for s_pdb in ['apo','holo',]:

        for prefix in ['analysis','comparison',]:
            fd = open('%s_%s.dat' %(prefix,s_pdb,),'w')
            fd.write('')
            fd.close()

    ##    for pdb_apo in d_apo2holo.keys():

##        d_distances_apo = {'POCASA': {'1ak1': 6.199449529617582, '2exo': 3.092818538880355, '2gg4': 3.263963045507762, '1cua': 3.947276140326851, '1gbs': 3.986649993630572, '1kf5': 4.834496488775966, '1sqg': 6.857717138845294, '3npo': 16.31846403263011, '2e3m': 1.4874884580729815, '3k0n': 9.7207517077744, '1akz': 7.516176825035834, '1bk7': 4.589443137004204, '1tgn': 1.5847483096146586, '1ahc': 18.87865614481024, '3mft': 2.31648083010426, '1erk': 7.83539471520612, '1rtc': 4.534131537790011, '1lp8': 2.6560551049003007, '2hbj': 4.7238129388940875, '1rd6': 1.626520822866456, '1mzl': 4.3062142808395265, '2ppn': 14.394389278525585, '2vfy': 4.840583322231743, '1yhv': 7.438932039657079, '1eur': 23.17458482502239, '1bbc': 1.7142531851812597, '1mri': 2.8523519407193234, '2zj8': 22.50346008036817, '2ggo': 2.8028848681416854, '2vfb': 7.999549096872514, '3c0e': 1.8774353149784135, '1tje': 16.35552772368349, '1wvw': 1.9662576815672046, '1wos': 10.702376595583777, '1qtr': 5.996833465773859, '1pdb': 2.437121149676831, '2qev': 17.519139445184422, '3g6l': 2.659602052491104, '1arl': 1.7830489896638724, '3pte': 4.467448729746946, '1znw': 9.507934900814027, '2paw': 12.375665347750903, '1kqx': 2.8331624773890223, '1sll': 1.9009744139347544, '1iad': 9.270834736558474, '3blm': 29.237391846523828, '1ri5': 3.1648895912264337, '1lmn': 5.54834844834817, '1gy0': 1.6520578306111504, '2sil': 24.965983246657736, '2zco': 5.026808936172409, '2ac4': 10.96821632804668, '1ey0': 5.695540215125622, '1sye': 6.640274864472619, '1tqo': 4.197534148536337, '1syc': 7.26348645543582, '1xqz': 12.00416165845169, '3aap': 4.365423594655016, '1p38': 6.770237004206174, '1z1i': 27.56167820610955, '4ape': 1.8458991181813589, '1ifb': 1.5290039475957449, '1jam': 3.087029085186702, '1mtz': 11.995868211890802, '1jcf': 2.517019560416303, '1xqo': 14.938912107603478, '3ewq': 1.1461503008298723, '2sga': 5.206116921610066, '1f10': 17.94048515839648}, 'GoodVibes': {'1ak1': 7.164883242180358, '2exo': 1.8720528864135584, '2gg4': 10.482140643503413, '1cua': 13.542178466051826, '1gbs': 13.405565221137103, '1kf5': 19.86002394717646, '1sqg': 13.576178722783311, '3npo': 17.04899814025863, '2e3m': 38.81344310671712, '3k0n': 26.482422928221343, '1akz': 8.257667277066014, '1bk7': 16.346059874367935, '1tgn': 28.368437815536275, '1ahc': 16.231684376276476, '3mft': 8.62441421836171, '1erk': 29.60000378210547, '1rtc': 27.75987218423402, '1lp8': 4.717177584951407, '2hbj': 19.312912790996588, '1rd6': 29.943921129299333, '1mzl': 8.265256800344819, '2ppn': 5.819044821621574, '2vfy': 11.919923847538584, '1yhv': 10.047365376847814, '1eur': 37.53387074159419, '1bbc': 0.933188831052968, '1mri': 4.867064642289185, '2zj8': 36.51738461413, '2ggo': 26.944030655162557, '2vfb': 7.093044385587308, '3c0e': 29.912660608930157, '1tje': 14.503864531830128, '1wvw': 3.330755764853194, '1wos': 22.974609457076895, '1qtr': 18.712735167038208, '1pdb': 4.349689291399841, '2qev': 19.893642498465482, '3g6l': 16.267560918038843, '1arl': 23.43129567427903, '3pte': 22.948914951707568, '1znw': 16.1783312053967, '2paw': 11.888029990991212, '1kqx': 7.997057645585126, '1sll': 20.916028386658343, '1iad': 3.6972420435369684, '3blm': 4.320210423673805, '1ri5': 28.961357254945263, '1lmn': 4.5498442313977385, '1gy0': 3.252891152097238, '2sil': 25.55868777964332, '2zco': 19.121078713284767, '2ac4': 2.5019541216544003, '1ey0': 22.33844800073198, '1sye': 8.996412942526412, '1tqo': 3.5131669398923466, '1syc': 23.266371787971032, '1xqz': 11.842087295855945, '3aap': 3.9070229276320307, '1p38': 21.3211260997327, '1z1i': 34.75352186677474, '4ape': 8.55111223760082, '1ifb': 10.213598110757566, '1jam': 4.1557628052155735, '1mtz': 15.453681509679294, '1jcf': 16.826679482804074, '1xqo': 14.850768054407201, '3ewq': 24.518704285922084, '2sga': 7.250574349303759, '1f10': 5.976540504285892}, 'Concavity': {'1ak1': 4.872716952676343, '2exo': 5.455657198208549, '2gg4': 4.922386850719265, '1cua': 7.329890449516455, '1gbs': 6.4480202116423, '1kf5': 6.596326440516866, '1sqg': 11.807240280859205, '3npo': 4.147190924182665, '2e3m': 1.395631829562909, '3k0n': 6.737874539006127, '1akz': 10.211088963869344, '1bk7': 7.9454247574087296, '1tgn': 5.814112303436863, '1ahc': 7.24744774141453, '3mft': 3.9931355025187547, '1erk': 10.057968705134023, '1rtc': 6.419082932937224, '1lp8': 7.996470004090961, '2hbj': 7.090491894416727, '1rd6': 2.0868150798394476, '1mzl': 1.5252063916865464, '2ppn': 2.03450939210674, '2vfy': 4.095569084204522, '1yhv': 3.815637361132967, '1eur': 18.474258451207714, '1bbc': 1.1480407556738568, '1mri': 3.9993892609178063, '2zj8': 22.255233561058123, '2ggo': 5.382851153974144, '2vfb': 8.459373838692729, '3c0e': 1.4058698263527571, '1tje': 7.24860788146362, '1wvw': 7.9749540711430225, '1wos': 14.228126824706685, '1qtr': 6.6905182310762426, '1pdb': 2.5827275478539886, '2qev': 6.4357631701741935, '3g6l': 3.3988354601346966, '1arl': 5.066364286887843, '3pte': 1.1465768941200067, '1znw': 1.8543691574656544, '2paw': 11.74554806063864, '1kqx': 2.118447527308808, '1sll': 10.416223951357185, '1iad': 4.687865472164579, '3blm': 10.294021654540769, '1ri5': 5.94816212372463, '1lmn': 3.1531508990922705, '1gy0': 1.3963318122290185, '2sil': 16.218424627097292, '2zco': 6.614121432155386, '2ac4': 8.493531582302607, '1ey0': 4.346962375105388, '1sye': 3.4170190828158113, '1tqo': 5.844676918503678, '1syc': 2.977701634657883, '1xqz': 5.958939819141642, '3aap': 6.545476592206798, '1p38': 6.769515731510223, '1z1i': 17.99377072611475, '4ape': 3.2298426510915563, '1ifb': 1.5915086501165134, '1jam': 6.299550671483679, '1mtz': 7.8769276683737175, '1jcf': 2.15383612641595, '1xqo': 14.166224150134429, '3ewq': 0.7639208082474526, '2sga': 6.010030782460209, '1f10': 5.29622209819112}, 'null': {'1ak1': 13.10687509735474, '2exo': 8.880512246308427, '2gg4': 4.885145093663333, '1cua': 18.878059730787548, '1gbs': 11.691482876996274, '1kf5': 7.722728007787313, '1sqg': 13.177791501564954, '3npo': 5.200052079280746, '2e3m': 10.392553599058024, '3k0n': 11.641298945365582, '1akz': 13.721727145737688, '1bk7': 10.408434742045994, '1tgn': 12.58105976895593, '1ahc': 7.38144243716267, '3mft': 10.846572746227642, '1erk': 14.701536499292287, '1rtc': 6.673803516621123, '1lp8': 8.121224690764045, '2hbj': 8.633161057925475, '1rd6': 12.684835388880646, '1mzl': 2.2090932749397023, '2ppn': 9.031132967426009, '2vfy': 8.506051211044836, '1yhv': 10.168367680793327, '1eur': 12.27300162717592, '1bbc': 5.690365330562548, '1mri': 6.762743805512736, '2zj8': 20.562633055178246, '2ggo': 13.143420879143193, '2vfb': 11.64167205203945, '3c0e': 8.980419018676974, '1tje': 8.707377428018683, '1wvw': 11.125815402133538, '1wos': 12.644454898879868, '1qtr': 4.828187456915729, '1pdb': 7.336989990679621, '2qev': 6.174963826501516, '3g6l': 9.23687695596161, '1arl': 11.718062804565788, '3pte': 9.794589812296495, '1znw': 8.485467355603367, '2paw': 14.509293226788591, '1kqx': 4.831229736277777, '1sll': 10.194790586381915, '1iad': 8.590573082905909, '3blm': 12.654837724388917, '1ri5': 10.81112410797997, '1lmn': 12.204957948503525, '1gy0': 9.674539714744633, '2sil': 11.606632940608062, '2zco': 2.0465975214576275, '2ac4': 17.312964128262788, '1ey0': 11.90742045072594, '1sye': 11.687873964715914, '1tqo': 12.386565877253055, '1syc': 11.58940120349537, '1xqz': 10.922482957168333, '3aap': 8.675459727841238, '1p38': 12.570698116701697, '1z1i': 17.892355457726953, '4ape': 7.4465983665295346, '1ifb': 3.691746558203225, '1jam': 12.365202836881624, '1mtz': 7.050287329369915, '1jcf': 5.074851104643149, '1xqo': 14.568188927886801, '3ewq': 8.549604196963855, '2sga': 14.506878613535019, '1f10': 4.689998910424163}, 'Fpocket': {'1ak1': 10.475316927251527, '2exo': 5.225674089493345, '2gg4': 16.789368416762034, '1cua': 37.62178477632966, '1gbs': 2.322629074636926, '1kf5': 1.7316395904126278, '1sqg': 20.373592755383815, '3npo': 5.8760151868127695, '2e3m': 1.1098599862526928, '3k0n': 8.962245645599053, '1akz': 7.931931178498968, '1bk7': 18.862204130233238, '1tgn': 7.9592525555253415, '1ahc': 2.0194931174598323, '3mft': 3.9678822940046015, '1erk': 5.70533559323113, '1rtc': 5.122221272070662, '1lp8': 2.518103300446247, '2hbj': 9.105504315211743, '1rd6': 5.520107478626559, '1mzl': 2.907786675193471, '2ppn': 2.5758141406316626, '2vfy': 5.572992658187641, '1yhv': 9.268192980288458, '1eur': 23.517033902333814, '1bbc': 1.9422384433840385, '1mri': 2.247259040586352, '2zj8': 12.039278261029455, '2ggo': 4.961720125204553, '2vfb': 10.220399445780242, '3c0e': 3.486179417057517, '1tje': 0.65335853505364, '1wvw': 5.97290274110169, '1wos': 10.491918734300427, '1qtr': 6.4979443640555345, '1pdb': 1.0536177373988538, '2qev': 17.310094171109597, '3g6l': 1.7618733236837505, '1arl': 2.1393154764497537, '3pte': 28.564251183000774, '1znw': 16.567506804865484, '2paw': 11.593090875346663, '1kqx': 2.4905520685612554, '1sll': 27.80491681215723, '1iad': 9.02241256495526, '3blm': 2.85868424034716, '1ri5': 7.22687035850455, '1lmn': 4.1207975235312695, '1gy0': 7.151038385633744, '2sil': 23.35625489350367, '2zco': 9.796493674463154, '2ac4': 14.555548996566433, '1ey0': 5.350625473331994, '1sye': 5.977545698767248, '1tqo': 17.267761093357947, '1syc': 7.3102374193824975, '1xqz': 11.19660512068995, '3aap': 29.356551254783536, '1p38': 14.86129974890589, '1z1i': 28.380989496309546, '4ape': 4.722032958377486, '1ifb': 0.6656420648066298, '1jam': 0.5993120247869059, '1mtz': 12.618754520464263, '1jcf': 3.064449288345906, '1xqo': 15.82102392480593, '3ewq': 4.829757993991083, '2sga': 3.0977022560404204, '1f10': 19.069597528465692}}
##        d_distances_holo = {'POCASA': {'1ak1': 4.313553054066148, '2exo': 3.045647979869617, '2gg4': 2.668172438234905, '1cua': 6.111396936756368, '1gbs': 2.2279943747410003, '1kf5': 3.356800731351207, '1sqg': 24.18299182193001, '3npo': 4.861689998944126, '2e3m': 0.9648513946458285, '3k0n': 6.206716906258723, '1akz': 5.797195397224419, '1bk7': 5.165467417545419, '1tgn': 3.4325638525802864, '1ahc': 18.478346629826568, '3mft': 2.311234537966683, '1erk': 5.63555917077936, '1rtc': 23.98903215310535, '1lp8': 2.0460727076019136, '2hbj': 6.34328349583293, '1rd6': 1.732605841234986, '1mzl': 2.6284280414604817, '2ppn': 2.9912409189294804, '2vfy': 3.818522197794229, '1yhv': 1.9694431787725077, '1eur': 22.998523197914256, '1bbc': 0.8583710134147318, '1mri': 19.38873344546108, '2zj8': 23.554771532493604, '2ggo': 2.961604372460544, '2vfb': 8.34665500263577, '3c0e': 2.315355569345666, '1tje': 7.553831705118332, '1wvw': 5.654896282778279, '1wos': 11.38739716490967, '1qtr': 5.500815913970192, '1pdb': 3.865934256657761, '2qev': 7.4370784948433535, '3g6l': 2.6609625588865247, '1arl': 0.4528936004552661, '3pte': 3.530072751235176, '1znw': 15.997900396636256, '2paw': 12.137680462041676, '1kqx': 1.3268309032226757, '1sll': 2.3802640184928143, '1iad': 7.836574832092735, '3blm': 3.299600016837058, '1ri5': 11.889655308850696, '1lmn': 5.646686218732685, '1gy0': 2.0961915526542807, '2sil': 25.258774673121025, '2zco': 3.6875714701171862, '2ac4': 6.081941118577914, '1ey0': 3.601771645398971, '1sye': 3.724078260402164, '1tqo': 3.671693140337306, '1syc': 3.5141680749787687, '1xqz': 12.34961767442518, '3aap': 7.097264715379188, '1p38': 6.701758063083585, '1z1i': 4.460665894188685, '4ape': 2.2408108543909386, '1ifb': 1.4602747891669379, '1jam': 4.538421850928176, '1mtz': 1.293122666559223, '1jcf': 2.3208805995880906, '1xqo': 13.59121197237851, '3ewq': 0.4109109497713292, '2sga': 4.9206767112654415, '1f10': 4.811703955843258}, 'GoodVibes': {'1ak1': 6.93216805218481, '2exo': 15.990792296942349, '2gg4': 11.41657322317677, '1cua': 12.894569320963793, '1gbs': 12.826757036629898, '1kf5': 1.2710826684366348, '1sqg': 26.449901157420445, '3npo': 14.313871930607501, '2e3m': 30.73167247667705, '3k0n': 26.304349087333147, '1akz': 18.18032525764047, '1bk7': 15.457302060031582, '1tgn': 31.36363467961465, '1ahc': 16.25722497037056, '3mft': 5.872736235856486, '1erk': 2.4846547014700437, '1rtc': 27.875990498122537, '1lp8': 4.682490098739485, '2hbj': 18.56197286962929, '1rd6': 31.155250962493707, '1mzl': 11.62827734314657, '2ppn': 9.432588631742828, '2vfy': 6.603608271985085, '1yhv': 6.845481136484798, '1eur': 39.00399363251153, '1bbc': 4.3605505153215, '1mri': 5.2913151360640605, '2zj8': 31.565117701291058, '2ggo': 10.388195743429417, '2vfb': 17.108206762308946, '3c0e': 30.559723025781828, '1tje': 12.319511898985729, '1wvw': 12.651150844222244, '1wos': 24.81782731075964, '1qtr': 8.649707922696532, '1pdb': 4.996142479633661, '2qev': 10.957251813647995, '3g6l': 15.58325773834842, '1arl': 19.509231129903608, '3pte': 23.50817935076287, '1znw': 16.30674427245859, '2paw': 28.080822867653907, '1kqx': 10.760480524368408, '1sll': 20.952021285025616, '1iad': 9.35442998110537, '3blm': 3.932497120716611, '1ri5': 12.13925516097596, '1lmn': 5.010784650511457, '1gy0': 7.855598534023033, '2sil': 25.841221100936774, '2zco': 10.5393619358047, '2ac4': 1.7809057137179152, '1ey0': 9.917453250991407, '1sye': 7.283511714111538, '1tqo': 22.964747647574974, '1syc': 6.777377756861425, '1xqz': 8.575136146457433, '3aap': 4.385143392090926, '1p38': 22.41781813272108, '1z1i': 30.440218358414416, '4ape': 10.213998765672349, '1ifb': 8.960108643555667, '1jam': 4.681614209507197, '1mtz': 16.596851325570068, '1jcf': 16.96245102637176, '1xqo': 14.453454876884106, '3ewq': 4.468810807235434, '2sga': 2.5876169575523233, '1f10': 6.249685406958527}, 'Concavity': {'1ak1': 4.113559740334862, '2exo': 4.433574015505679, '2gg4': 5.160122588660032, '1cua': 11.508522228068731, '1gbs': 7.521928834540698, '1kf5': 8.54199116216703, '1sqg': 12.73802566936112, '3npo': 4.143577302878453, '2e3m': 1.474429169118829, '3k0n': 6.50189935789069, '1akz': 9.454300294727439, '1bk7': 7.7108496402037, '1tgn': 8.928446055907123, '1ahc': 6.172612343690967, '3mft': 4.062553691778567, '1erk': 12.190528276990817, '1rtc': 6.448330467849562, '1lp8': 9.424607062121755, '2hbj': 6.626482447663723, '1rd6': 2.4035998239699614, '1mzl': 2.6222189727503418, '2ppn': 2.529262107935774, '2vfy': 3.992037731565631, '1yhv': 3.652736681372341, '1eur': 17.792166860696902, '1bbc': 0.7636542501169263, '1mri': 4.816408920951069, '2zj8': 21.832816006507077, '2ggo': 3.79834676443031, '2vfb': 9.602235361563228, '3c0e': 1.4647374773332946, '1tje': 4.219803156482689, '1wvw': 7.4431090067436605, '1wos': 11.809716929403747, '1qtr': 5.846418421584703, '1pdb': 2.470294738709091, '2qev': 5.972950116843509, '3g6l': 3.366660623217371, '1arl': 5.241775772123677, '3pte': 1.629190814375921, '1znw': 2.3462964570110154, '2paw': 11.603970855165617, '1kqx': 1.1574811911879224, '1sll': 9.6602507892037, '1iad': 4.580776548542596, '3blm': 7.825147410810537, '1ri5': 5.193654240596455, '1lmn': 3.258409043719367, '1gy0': 2.069857255543681, '2sil': 16.33153942348782, '2zco': 6.487508580727465, '2ac4': 8.896136721550338, '1ey0': 2.9004244154068286, '1sye': 2.781921103010747, '1tqo': 4.4757873285760255, '1syc': 2.0592799553138668, '1xqz': 3.7088617268714494, '3aap': 7.566400666655647, '1p38': 6.062997400616175, '1z1i': 12.1436131230906, '4ape': 2.9003666881410934, '1ifb': 1.054883490213916, '1jam': 4.863267728653775, '1mtz': 6.790221937479254, '1jcf': 2.0034619758080434, '1xqo': 14.266102036559046, '3ewq': 0.6438966749250639, '2sga': 6.125502297208853, '1f10': 6.283973579269688}, 'null': {'1ak1': 12.982574666947409, '2exo': 8.871024796737933, '2gg4': 4.978396209511019, '1cua': 18.820171631540656, '1gbs': 11.7031454990019, '1kf5': 7.708597803633162, '1sqg': 13.116205046554894, '3npo': 5.192182444268401, '2e3m': 10.250026221393659, '3k0n': 11.851681267261876, '1akz': 13.66454549648406, '1bk7': 10.413713397339436, '1tgn': 12.53620972278645, '1ahc': 7.3784145499575775, '3mft': 10.862180758187097, '1erk': 15.09542303553847, '1rtc': 6.668822346736799, '1lp8': 8.12037528865089, '2hbj': 8.632854241038011, '1rd6': 12.661501423080107, '1mzl': 2.2128053563983396, '2ppn': 9.144136418670223, '2vfy': 8.838054173480199, '1yhv': 10.186946407935368, '1eur': 12.15059096119368, '1bbc': 5.7849975849458675, '1mri': 6.762238957938146, '2zj8': 20.569014632847942, '2ggo': 13.144320104999746, '2vfb': 11.741131128640895, '3c0e': 8.990033797584305, '1tje': 8.701434768032248, '1wvw': 11.155667994969376, '1wos': 12.61777023608208, '1qtr': 4.7966628763831585, '1pdb': 7.370684093029069, '2qev': 6.162968801189478, '3g6l': 9.24764699742328, '1arl': 11.716345468949246, '3pte': 9.736128769837874, '1znw': 8.504437042124092, '2paw': 14.43016404943335, '1kqx': 4.84778702538274, '1sll': 10.193766891580239, '1iad': 8.610962425721265, '3blm': 12.647463153582938, '1ri5': 10.82258667960808, '1lmn': 12.20840018188835, '1gy0': 9.687757205251232, '2sil': 11.611222645295062, '2zco': 2.056493043247595, '2ac4': 17.166245542234932, '1ey0': 11.85439436420821, '1sye': 11.71351525781101, '1tqo': 12.435970054665606, '1syc': 11.617156501711097, '1xqz': 10.938979153943178, '3aap': 8.665168934328593, '1p38': 12.564626856544919, '1z1i': 18.337975680811514, '4ape': 7.454546412647683, '1ifb': 3.7036319850931507, '1jam': 12.350391528608139, '1mtz': 7.278268091743817, '1jcf': 5.0425941818481075, '1xqo': 14.504850837658312, '3ewq': 8.560977752456745, '2sga': 14.509580356084465, '1f10': 4.754536569118486}, 'Fpocket': {'1ak1': 11.039606024827611, '2exo': 5.70008746391297, '2gg4': 1.0864671820429832, '1cua': 16.51095991486476, '1gbs': 5.04845276124027, '1kf5': 3.7282689773534217, '1sqg': 15.61969387442707, '3npo': 2.48300909724542, '2e3m': 2.48910876641288, '3k0n': 8.292874099872707, '1akz': 8.071041444514988, '1bk7': 18.15303948068212, '1tgn': 3.0189489615028644, '1ahc': 4.140758531204782, '3mft': 2.2173650290934184, '1erk': 31.435192586362234, '1rtc': 0.37011817280364334, '1lp8': 1.3502200215937354, '2hbj': 17.85515138027955, '1rd6': 11.701165555106597, '1mzl': 1.4368722757440395, '2ppn': 2.55032909152173, '2vfy': 5.047206777184247, '1yhv': 1.0185549267482792, '1eur': 24.32248325074462, '1bbc': 2.6640048837778165, '1mri': 2.3925427648376028, '2zj8': 26.893649186998466, '2ggo': 3.044262216120429, '2vfb': 19.111662201771246, '3c0e': 19.476327582797886, '1tje': 2.706333353441735, '1wvw': 7.132146180377451, '1wos': 15.357487549934586, '1qtr': 5.296128454785188, '1pdb': 4.528968785428489, '2qev': 1.1559537073146031, '3g6l': 2.7254428535613164, '1arl': 2.8933078373830003, '3pte': 4.650040745812233, '1znw': 5.000491480929575, '2paw': 11.557005924591047, '1kqx': 3.5544438009395742, '1sll': 21.42142926357061, '1iad': 8.921306693041458, '3blm': 2.143135523869571, '1ri5': 5.829309968670052, '1lmn': 15.229778689753523, '1gy0': 7.2817916290986435, '2sil': 22.259225433499505, '2zco': 3.7773192329364838, '2ac4': 1.0546520964459056, '1ey0': 4.307672126039723, '1sye': 4.265896509691928, '1tqo': 5.235864948423655, '1syc': 3.5495585958514746, '1xqz': 1.8191633541773506, '3aap': 5.443946278798127, '1p38': 10.406574686817086, '1z1i': 3.3140263109301173, '4ape': 2.306413810124328, '1ifb': 0.25975332031276044, '1jam': 3.8768598959211875, '1mtz': 12.373694325834355, '1jcf': 5.812015920389089, '1xqo': 14.903846117735075, '3ewq': 2.676660023702034, '2sga': 2.9280977018730256, '1f10': 19.968746273522406}}

        for pdb_apo in l_pdbs_apo:

            print pdb_apo

            pdb_holo = d_apo2holo[pdb_apo]['holo']

            if s_pdb == 'apo':
                pdb = pdb_apo
            elif s_pdb == 'holo':
                pdb = pdb_holo

            (
                position_ligand, chain, n_residues_protein, n_atoms_ligand,
                ligand, overlap_site,
                ) = get_position_ligand(pdb, pdb_apo,d_apo2holo,)

            position_Fpocket = get_position_Fpocket(pdb)

            position_POCASA = get_position_POCASA(pdb)

            position_Concavity, pocket_size_Concavity = get_position_Concavity(pdb)

            position_GoodVibes, overlap_min = get_position_GoodVibes(pdb,chain,)

            position_null = get_position_null(pdb, ligand,)

            d_positions = {
                'Fpocket':position_Fpocket,
                'POCASA':position_POCASA,
                'Concavity':position_Concavity,
                'GoodVibes':position_GoodVibes,
                'null':position_null,
                }

            l_methods = d_positions.keys()
            l_methods.sort()
##            print
##            print 'overlap', d_apo2holo[pdb_apo]['overlap'], 'rmsd', d_apo2holo[pdb_apo]['rmsd']
            dist_GV = math.sqrt(sum((d_positions['GoodVibes']-position_ligand)**2))
##            if dist_GV > 15 and d_apo2holo[pdb_apo]['overlap'] > 0.6:
##            if dist_GV > 30:
##            if pdb == pdb_holo and dist_GV > 15 and d_apo2holo[pdb_apo]['overlap'] > 0.4:
            if pdb == pdb_holo and dist_GV > 15 and overlap_site > 0.4:
                print '**********', pdb, dist_GV, pdb_apo, pdb_holo, overlap_min, overlap_site, d_apo2holo[pdb_apo]['overlap']
##                stop
            l = [str(dist_GV)]

            d_dist = {}
            for method in l_methods:
                dist = math.sqrt(sum((d_positions[method]-position_ligand)**2))
                d_dist[method] = dist
##                print pdb, method, dist
                l += [str(dist)]
                if method == 'null':
                    dist_null = dist
                d_distances[method][pdb] = dist
                if s_pdb == 'holo' and pdb_apo in d_distances[method].keys():
                    d_differences[method][pdb] = abs(d_distances[method][pdb_apo]-d_distances[method][pdb_holo])

            fd = open('comparison_%s.dat' %(s_pdb),'a')
            fd.write(' '.join(l)+'\n')
            fd.close()

            fd = open('analysis_%s.dat' %(s_pdb),'a')
            fd.write('%f %f %f %f %f %f %f %f %f %f\n' %(
                dist_GV,
                d_apo2holo[pdb_apo]['overlap'], ## correlation (higher overlap better)
                d_apo2holo[pdb_apo]['rmsd'], ## no correlation
                n_residues_protein, ## weak correlation (larger proteins better)
                n_atoms_ligand, ## no correlation
                dist_null, ## no correlation
                d_dist['Concavity'], ## 7, no correlation
                pocket_size_Concavity, ## no correlation
                overlap_min, ## no correlation
                overlap_site,
                ))
            fd.close()

        if s_pdb == 'apo': l_pdbs = l_pdbs_apo
        else: l_pdbs = l_pdbs_holo

        ##
        ##
        ##
        for method in d_distances.keys():
            l = []
            for pdb in l_pdbs:
                try:
                    l += [d_distances[method][pdb]]
                except:
                    print pdb, method
                    print d_distances[method].keys()
                    stop
            print 'average dist', s_pdb, len(l), method, sum(l)/len(l)
        if s_pdb == 'holo':
            for method in d_differences.keys():
                l = d_differences[method].values()
                print 'average diff', len(l), method, sum(l)/len(l)
                print 'less than 6A', [l[i] < 6 for i in range(len(l))].count(True)
                

        print 'f(x) = x'
        print 'set xlabel "distance (null model)"'
        print 'set ylabel "distance (individual methods)"'
        print 'plot "comparison_%s.dat" u 6:2 t "Concavity", "comparison_%s.dat" u 6:3 t "Fpocket", "comparison_%s.dat" u 6:4 t "GoodVibes", "comparison_%s.dat" u 6:5 t "POCASA", "comparison_%s.dat" u 6:6 t "Null model", f(x) t""' %(s_pdb,s_pdb,s_pdb,s_pdb,s_pdb,)
        print 'horizontal line = average, g(x)=4'

        print 'g(x) = a*x+b'
        print 'fit g(x) "analysis_%s.dat" u 2:1 via a,b' %(s_pdb)
        print 'set xlabel "overlap between apo/holo motion and normal mode 7"'
        print 'set ylabel "distance between center of ligand and center of predicted pocket"'
        print 'plot [0:1][0:]g(x) t "regression", "analysis_%s.dat" u 2:1 t "GoodVibes"' %(s_pdb)

    return


def get_position_ligand(pdb,pdb_apo,d_apo2holo,):

    pdb_holo = d_apo2holo[pdb_apo]['holo']
    d_mmCIF_holo = parse_mmCIF.main(pdb_holo,)
    d_coords, l_coords_alpha_holo = mmCIF2coords.main(pdb_holo,d_mmCIF_holo)

    ##
    ##
    ##
    ligand = d_apo2holo[pdb_apo]['ligand']

    l_residues = []
    for i in range(len(d_mmCIF_holo['_struct_site.id'])):
        if not 'BINDING SITE FOR RESIDUE %s' %(ligand) in d_mmCIF_holo['_struct_site.details'][i]:
            continue
        if len(l_residues) > 0:
            print pdb, pdb_apo, pdb_holo
            print l_residues
            print d_mmCIF_holo['_struct_site.details'][i]
            stop
        struct_site_ID = d_mmCIF_holo['_struct_site.id'][i]
        for j in range(len(d_mmCIF_holo['_struct_site_gen.site_id'])):
            struct_site_gen_ID = d_mmCIF_holo['_struct_site_gen.site_id'][j]
            if struct_site_ID == struct_site_gen_ID:
                residue = int(d_mmCIF_holo['_struct_site_gen.auth_seq_id'][j])
##                l_residues += [residue]
                ## include neighboring residues
                l_residues += [residue-1,residue,residue+1]
    l_residues = list(set(l_residues))
    if len(l_residues) == 0:
        print pdb
        stop

    ## 
    l_coords_ligand = []
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


    d_mmCIF_apo = parse_mmCIF.main(pdb_apo,)
    d_coords, l_coords_alpha_apo = mmCIF2coords.main(pdb_apo,d_mmCIF_apo)   

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


    if pdb == pdb_apo:
        l_seq_num = d_mmCIF_apo['_pdbx_poly_seq_scheme.pdb_seq_num'][index1_coord_apo:index2_coord_apo]
        chain = ''.join(d_mmCIF_apo['_entity_poly.pdbx_strand_id'])
        n_residues = len(l_coords_alpha_apo)
        l_coords_alpha = l_coords_alpha_apo
    else:
        l_seq_num = d_mmCIF_holo['_pdbx_poly_seq_scheme.pdb_seq_num'][index1_coord_holo:index2_coord_holo]
        chain = ''.join(d_mmCIF_holo['_entity_poly.pdbx_strand_id'])
        n_residues = len(l_coords_alpha_holo)
        l_coords_alpha = l_coords_alpha_holo

    overlap_site = 1.
##    ##
##    ## eigenvector
##    ##
##    cutoff = 10
##    matrix_hessian = NMA.hessian_calculation(l_coords_alpha,cutoff,)
##    eigenvectors, eigenvalues = NMA.diagonalize_hessian(matrix_hessian)
##
##    ## apply transformation matrix
##    if pdb == pdb_apo:
##        instance_geometry = geometry.geometry()
##        rmsd = instance_geometry.superpose(l_coords_alpha_apo,l_coords_alpha_holo,)
##        tv1 = instance_geometry.fitcenter
##        rm = instance_geometry.rotation
##        tv2 = instance_geometry.refcenter
##        for i_coord in range(len(l_coords_ligand)):
##            l_coords_ligand[i_coord] = numpy.dot(l_coords_ligand[i_coord]-tv1,rm)+tv2
##
##    ##
##    ## apo/holo eigenvector
##    ##
##    vector_apo2holo = []
##    for i in range(len(l_coords_alpha_holo)):
##        vector_apo2holo += [
##            l_coords_alpha_holo[i][0]-l_coords_alpha_apo[i][0],
##            l_coords_alpha_holo[i][1]-l_coords_alpha_apo[i][1],
##            l_coords_alpha_holo[i][2]-l_coords_alpha_apo[i][2],
##            ]
##    vector_apo2holo = numpy.array(vector_apo2holo)
##
##    ##
##    ## calculate overlap between normal modes and difference vector
##    ## in the ligand binding site!!!
##    ##
##    vector_apo2holo_site = []
##    eigenvector_site = []
##    ## exclude coordinate not at the ligand binding site
##    for i_seq_num in range(len(l_seq_num)):
##        seq_num = int(l_seq_num[i_seq_num])
##        if seq_num in l_residues:
##            eigenvector_site += list(eigenvectors[6][3*i_seq_num:3*i_seq_num+3])
##            vector_apo2holo_site += list(vector_apo2holo[3*i_seq_num:3*i_seq_num+3])
##    ## calculate overlap
##    vector_apo2holo_site = numpy.array(vector_apo2holo_site)
##    eigenvector_site = numpy.array(eigenvector_site)
##    overlap_site = abs(
##        numpy.dot(eigenvector_site,vector_apo2holo_site)
##        /
##        math.sqrt(
##            numpy.dot(eigenvector_site,eigenvector_site)
##            *
##            numpy.dot(vector_apo2holo_site,vector_apo2holo_site)
##            )
##        )
##    if overlap_site > 0.8:
##        print vector_apo2holo_site
##        print eigenvector_site
##        print pdb
##        print l_residues

    position_ligand = sum(l_coords_ligand)/len(l_coords_ligand)

    n_atoms = len(l_coords_ligand)

    return position_ligand, chain, n_residues, n_atoms, ligand, overlap_site


def get_position_GoodVibes(pdb,chain,):

    fd = open('output/GoodVibes/distmax6_distmin3/%s_%s_probe.pdb' %(pdb,chain,),'r')
    lines = fd.readlines()
    fd.close()
    min_overlap = [1,'position',]
    for line in lines:
        record = line[:6].strip()
        if record != 'HETATM':
            continue
        res_name = line[17:20]
        if res_name != 'EXT':
            continue
        overlap = (100-float(line[60:66]))/100.
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = numpy.array([x,y,z,])
        if overlap < min_overlap[0]:
            min_overlap = [overlap,coord,]

    position = min_overlap[1]
    overlap = min_overlap[0]

    return position, overlap


def get_position_POCASA(pdb):

    fd = open('output/POCASA/%s_Pocket_DepthCenters.pdb' %(pdb),'r')
    lines = fd.readlines()
    fd.close()
    l_coords = []
    for line in lines:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        position = [x,y,z,]
        break

    return position


def get_position_null(pdb,ligand,):

    fd = open('/media/WDMyBook1TB/2TB/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()
    l_coords = []
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            iCode = line[16]
            if iCode not in [' ','A','1',]:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            l_coords += [coord]
        elif (
            record == 'HETATM'
            and
            line[17:20] != 'HOH'
            and
            line[17:20] != ligand
            and
            line[17:20].strip() not in [
                'MG','CA','CL','CO','MN','SO4','ZN','IOD','NA',
                'GOL','FMT','EDO','DMS','ACT',
##                'URA',
##                'NDP',
##                'SAH',
                ]
            ):
            print line
            print line[17:20]
            print pdb
            stop
            continue
    position = sum(l_coords)/len(l_coords)

    return position


def get_position_Concavity(pdb):

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

    return position, len(l_coords)


def get_position_Fpocket(pdb):

    fd = open('output/Fpocket/%s_out/pockets/pocket0_vert.pqr' %(pdb), 'r')
    lines = fd.readlines()
    fd.close()
    l_coords = []
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coord = numpy.array([x,y,z,])
            l_coords += [coord]
    position = sum(l_coords)/len(l_coords)
    
    return position


if __name__ == '__main__':
    main()

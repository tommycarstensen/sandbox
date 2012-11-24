import os,time

import sys

## HEWL wt
l_wts = [
    ## all spacegroups
    ## one chain
    '193l_a', '194l_a', '1aki_a', '1at6_a', '1azf_a', '1b0d_a', '1bgi_a', '1bhz_a', '1bvx_a', '1bwh_a', '1bwi_a', '1bwj_a', '1c10_a', '1dpw_a', '1dpx_a', '1f0w_a', '1f10_a', '1h87_a', '1hel_a', '1hew_a', '1hsw_a', '1hsx_a', '1iee_a', '1jis_a', '1jit_a', '1jiy_a', '1jj0_a', '1jj1_a', '1jpo_a',
    '2d4j_a', '1lks_a', '1lma_a', '1lpi_a', '1lsa_a', '1lsb_a', '1lsc_a', '1lsd_a', '1lse_a', '1lsf_a', '1lyo_a', '1lys_a', '1lyz_a', '1lz8_a', '1lz9_a', '1lza_a', '1lzb_a', '1lzc_a', '1lzh_a', '1lzt_a', '1n4f_a', '1ps5_a', '1qio_a', '1qtk_a', '1rcm_a', '1rfp_a', '1t3p_a', '1uc0_a', '1uco_a', '1uig_a', '1uih_a', '1v7s_a', '1v7t_a', '1vau_a', '1vdp_a', '1vdq_a', '1vds_a', '1vdt_a', '1ved_a', '1w6z_a', '1wtm_a', '1wtn_a', '1xei_a', '1xej_a', '1xek_a', '1yik_a', '1yil_a',
    '1z55_a', '2a7d_a', '2a7f_a', '2aub_a', '2b5z_a', '2blx_a', '2bly_a', '2bpu_a', '2c8o_a', '2c8p_a', '2cds_a', '2cgi_a', '2d6b_a', '2epe_a', '2f2n_a', '2f30_a', '2f4a_a', '2f4g_a', '2fbb_a', '2g4p_a', '2g4q_a', '2h9j_a', '2h9k_a', '2htx_a', '2hu1_a', '2hu3_a', '2hub_a', '2i6z_a', '2lym_a', '2lyo_a', '2lyz_a', '2lzh_a', '2lzt_a', '2pc2_a', '2vb1_a', '2w1l_a', '2w1m_a', '2w1x_a', '2w1y_a', '2yvb_a', '2z12_a', '2z18_a', '2z19_a', '2zq3_a', '2zq4_a', '2zxs_a', '2zyp_a', '3e3d_a', '3ems_a', '3exd_a', '3lym_a', '3lyo_a', '3lyt_a', '3lyz_a', '3lzt_a', '4lym_a', '4lyo_a', '4lyt_a', '4lyz_a', '4lzt_a', '5lym_a', '5lyt_a', '5lyz_a', '6lyt_a', '6lyz_a', '7lyz_a', '8lyz_a',
    ## two chains
    '1jj3_a', '1hf4_a', '1lj3_a', '1lj4_a', '1lje_a', '1ljf_a', '1ljg_a', '1ljh_a', '1lji_a', '1ljj_a', '1ljk_a', '2d4i_a', '2d4k_a', 
    '1jj3_b', '1hf4_b', '1lj3_b', '1lj4_b', '1lje_b', '1ljf_b', '1ljg_b', '1ljh_b', '1lji_b', '1ljj_b', '1ljk_b', '2d4i_b', '2d4k_n',
    ## one chain but not "a"
    '1ykx_x', '1yky_x', '1ykz_x', '1yl0_x', '1yl1_x', '2q0m_x',
    ]

d_mutants = {
    '1flq_a':'G117A',
    '1flu_a':'G67A' ,
    '1flw_a':'G71A' ,
    '1fly_a':'G102A',
    '1fn5_a':'G49A' ,
    '1hem_a':'S91T' ,
    '1heo_a':'I55V' ,
    '1her_a':'T40S' ,
    '1ios_a':'M12F' ,
    '1iot_a':'M12L' ,
    '1ir7_a':'I78M' ,
    '1ir8_a':'I58M' ,
    '1ir9_a':'I98M' ,
    '1kxw_a':'N27D' ,
    '1kxy_a':'D18N' ,
    '1lzd_a':'W62Y' ,
    '1uic_a':'H15A' ,
    '1uid_a':'H15F' ,
    '1uie_a':'H15G' ,
    '1uif_a':'H15V' ,
    }

l_pdbs = l_wts+d_mutants.keys()

if '-manual' in sys.argv:
    index = sys.argv.index('-manual')
    l_pdbs = sys.argv[index+1:]

## go to maindir
dir_main = '/home/tc/UFFBAPS_out'
os.chdir(dir_main)

for i in range(len(l_pdbs)):

    pdb = l_pdbs[i]

    if not os.path.isdir(pdb):
        os.mkdir(pdb)

    ## go to subdir
    os.chdir('%s/%s/' %(dir_main,pdb))

    ## copy remediated pdb
    os.system('cp /data/pdb-v3.2/%s/pdb%s.ent %s.pdb' %(pdb[1:3],pdb[:4],pdb,))

    ##
    ## 1) delete all chains except one
    ## 2) delete all altlocs except A or 1
    ##
    fd = open('%s.pdb' %(pdb),'r')
    lines_in = fd.readlines()
    fd.close()
    res_no_initial = None
    lines_out = []
    for line in lines_in:
        record = line[:6].strip()
        if record in ['ATOM','HETATM',]:
            chain = line[21]
            if chain != pdb[-1].upper():
                continue
            altloc = line[16]
            if altloc not in [' ','A','1',]:
                continue
            res_no = int(line[22:26])
            if res_no_initial == None:
                res_no_initial = res_no
            lines_out += [line]
        else:
            lines_out += [line]
    fd = open('%s.pdb' %(pdb),'w')
    fd.writelines(lines_out)
    fd.close()

    ##
    ## create mutation list
    ##
    l = []
    chain = pdb[-1].upper() ## UFFBAPS is case sensitive
    if pdb in d_mutants.keys():
        mutation = d_mutants[pdb]
        mutation = mutation[-1]+mutation[1:-1]+mutation[0]
        l_mutations = [mutation]
    else:
        l_mutations = d_mutants.values()
    for mutation in l_mutations:
        res_no = int(mutation[1:-1])+res_no_initial-1
        l += ['%1s%i%1s\n' %(chain,res_no,mutation[-1],)]
    fd = open('mutation.list','w')
    fd.writelines(l)
    fd.close()
        

    ##
    ## write pbs script
    ##
    s = '#PBS -l nodes=1:ppn=1\n'
##    s = '#PBS -l nodes=node%i:ppn=1\n' %(i%11+1)
##    s = '#PBS -l nodes=node1:ppn=1\n'
    s += '#PBS -d %s/%s\n' %(dir_main,pdb)
    s += '#PBS -o %s/%s\n' %(dir_main,pdb)
    s += '''\

#PBS -q short

#PBS -j oe

'''
    s += 'python /home/tc/svn/PEAT/PEAT_SA/Core/ProteinDesignTool.py -p %s.pdb -w . --stability --mutationList mutation.list\n' %(pdb)

    fd = open('%s.pbs' %(pdb,),'w')
    fd.write(s)
    fd.close()


    ##
    ## execute pbs script
    ##
    os.system('qsub %s.pbs' %(pdb,))


    ##
    ## go to main dir
    ##
    os.chdir(dir_main)

import os
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import gnuplot

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

d_ddG = {}
for mutation in d_mutants.values():
    d_ddG[mutation] = {'forward':[],'backward':[],}

l_fn = os.listdir('/home/people/tc/UFFBAPSout/)
for fn in l_fn:
    if fn[-10:] != '.stability':
        continue
    pdb = fn[:6]
    fd = open(fn,'r')
    lines = fd.readlines()
    fd.close()
    for line in lines[1:]:
        l = line.split(',')
        mutation = l[0]
        ddG = float(l[-1])
        res_no = int(mutation[1:-2])
        if pdb in ['2d4i_b','2d4k_n',]:
            res_no -= 200
        res_symbol1 = mutation[-2]
        res_symbol2 = mutation[-1]
        if pdb in d_mutants.keys():
            mutation = '%1s%i%1s' %(res_symbol2,res_no,res_symbol1,)
            d_ddG[mutation]['backward'] += [ddG]
        else:
            mutation = '%1s%i%1s' %(res_symbol1,res_no,res_symbol2,)
            d_ddG[mutation]['forward'] += [ddG]

l = []
for mutation in d_ddG.keys():
    for ddG_forward in d_ddG[mutation]['forward']:
        for ddG_backward in d_ddG[mutation]['backward']:
            l += ['%s %s\n' %(-ddG_backward,ddG_forward,)]

fd = open('UFFBAPS.gnuplotdata','w')
fd.writelines(l)
fd.close()

prefix = 'UFFBAPS'
xlabel = 'ddG backward'
ylabel = 'ddG forward'
gnuplot.scatter_plot_2d(
    prefix, xlabel=xlabel, ylabel=ylabel,
    bool_multiple_columns = False,
    function = 'x',
    bool_remove = False,
    )

l_pdbs = [
    '107L_A', '108L_A', '109L_A', '110L_A', '111L_A', '112L_A', '113L_A', '114L_A', '115L_A', '118L_A', '119L_A', '120L_A', '122L_A', '123L_A', '125L_A', '126L_A', '127L_A', '128L_A', '129L_A', '130L_A', '131L_A', '137L_A', '137L_B', '138L_A', '139L_A', '141L_A', '142L_A', '143L_A', '145L_A', '146L_A', '147L_A', '152L_A', '155L_A', '156L_A', '157L_A', '158L_A', '159L_A', '160L_A', '161L_A', '162L_A', '163L_A', '164L_A', '165L_A', '166L_A', '172L_A', '173L_A', '180L_A', '180L_B', '181L_A', '182L_A', '183L_A', '184L_A', '185L_A', '186L_A', '187L_A', '188L_A', '190L_A', '191L_A', '192L_A', '195L_A', '198L_A', '199L_A', '1B6I_A', '1C60_A', '1C61_A', '1C63_A', '1C64_A', '1C65_A', '1C69_A', '1C6C_A', '1C6D_A', '1C6E_A', '1C6F_A', '1C6G_A', '1C6H_A', '1C6I_A', '1C6J_A', '1C6K_A', '1C6L_A', '1C6P_A', '1C6Q_A', '1C6T_A', '1CU2_A', '1CUP_A', '1CV3_A', '1CV4_A', '1CV5_A', '1CV6_A', '1CVK_A', '1CX7_A', '1D2W_A', '1D3J_A', '1D3N_A', '1D9W_A', '1DYA_A', '1DYB_A', '1DYE_A', '1DYF_A', '1EPY_A', '1G06_A', '1G07_A', '1G0G_A', '1G0J_A', '1G0K_A', '1G0L_A', '1G0M_A', '1G0P_A', '1G0Q_A', '1G1V_A', '1G1W_A', '1I6S_A', '1KNI_A', '1KW5_A', '1KW7_A', '1KY0_A', '1L00_A', '1L01_A', '1L02_A', '1L03_A', '1L04_A', '1L05_A', '1L06_A', '1L07_A', '1L08_A', '1L09_A', '1L0J_A', '1L10_A', '1L11_A', '1L12_A', '1L13_A', '1L14_A', '1L15_A', '1L16_A', '1L17_A', '1L18_A', '1L19_A', '1L20_A', '1L21_A', '1L22_A', '1L23_A', '1L24_A', '1L25_A', '1L26_A', '1L27_A', '1L28_A', '1L29_A', '1L30_A', '1L31_A', '1L32_A', '1L33_A', '1L34_A', '1L35_A', '1L36_A', '1L37_A', '1L38_A', '1L39_A', '1L40_A', '1L41_A', '1L42_A', '1L43_A', '1L44_A', '1L45_A', '1L46_A', '1L47_A', '1L48_A', '1L49_A', '1L50_A', '1L51_A', '1L52_A', '1L53_A', '1L54_A', '1L55_A', '1L56_A', '1L57_A', '1L58_A', '1L59_A', '1L60_A', '1L61_A', '1L62_A', '1L63_A', '1L64_A', '1L65_A', '1L66_A', '1L67_A', '1L68_A', '1L69_A', '1L70_A', '1L71_A', '1L72_A', '1L73_A', '1L74_A', '1L75_A', '1L76_A', '1L79_A', '1L80_A', '1L81_A', '1L83_A', '1L84_A', '1L85_A', '1L86_A', '1L87_A', '1L88_A', '1L89_A', '1L90_A', '1L91_A', '1L92_A', '1L93_A', '1L94_A', '1L95_A', '1L96_A', '1L97_A', '1L97_B', '1L98_A', '1L99_A', '1LGU_A', '1LGW_A', '1LGX_A', '1LI2_A', '1LI3_A', '1LI6_A', '1LLH_A', '1LPY_A', '1LW9_A', '1LWG_A', '1LYD_A', '1LYE_A', '1LYF_A', '1LYG_A', '1LYH_A', '1LYI_A', '1LYJ_A', '1NHB_A', '1OV7_A', '1OVH_A', '1OVJ_A', '1OWY_A', '1OWZ_A', '1P2L_A', '1P2R_A', '1P36_A', '1P37_A', '1P3N_A', '1P46_A', '1P64_A', '1P6Y_A', '1P7S_A', '1PQD_A', '1PQI_A', '1PQJ_A', '1PQK_A', '1PQK_B', '1PQK_C', '1PQM_A', '1PQO_A', '1QS9_A', '1QSB_A', '1QSQ_A', '1QT3_A', '1QT5_A', '1QT6_A', '1QT7_A', '1QT8_A', '1QTB_A', '1QTH_A', '1QTH_B', '1QTZ_A', '1QUD_A', '1QUG_A', '1QUH_A', '1QUO_A', '1SWY_A', '1SWZ_A', '1SX2_A', '1SX7_A', '1T8G_A', '1TLA_A', '1XEP_A', '1ZUR_A', '1ZWN_A', '1ZYT_A', '200L_A', '206L_A', '217L_A', '219L_A', '220L_A', '221L_A', '222L_A', '223L_A', '224L_A', '225L_A', '226L_A', '227L_A', '228L_A', '229L_A', '230L_A', '232L_A', '233L_A', '234L_A', '235L_A', '236L_A', '237L_A', '238L_A', '239L_A', '240L_A', '241L_A', '242L_A', '243L_A', '244L_A', '245L_A', '246L_A', '247L_A', '248L_A', '249L_A', '250L_A', '253L_A', '254L_A', '255L_A', '256L_A', '257L_A', '258L_A', '259L_A', '260L_A', '2A4T_A', '2CUU_A', '2HUK_A', '2HUL_A', '2IGC_A', '2L78_A', '2LZM_A', '2NTG_A', '2NTH_A', '2OTY_X', '2OU0_X', '2OU8_A', '2OU9_A', '2Q9D_A', '2RAY_X', '2RAZ_X', '2RB0_X', '2RB1_X', '2RB2_X', '2RBN_A', '2RBO_A', '2RBP_A', '2RBQ_A', '2RBR_A', '2RBS_A', '3C7W_A', '3C7Y_A', '3C7Z_A', '3C80_A', '3C81_A', '3C82_A', '3C83_A', '3C8Q_A', '3C8R_A', '3C8S_A', '3CDO_A', '3CDO_B', '3CDO_C', '3CDO_D', '3CDQ_A', '3CDR_A', '3CDT_A', '3CDV_A', '3DKE_X', '3DMV_A', '3DMX_A', '3DMZ_A', '3DN0_A', '3DN1_A', '3DN2_A', '3DN3_A', '3DN4_A', '3DN6_A', '3DN8_A', '3DNA_A', '3F8V_A', '3F9L_A', '3FA0_A', '3FAD_A', '3FI5_A', '3FI5_B', '3FI5_C', '3FI5_D', '3G3X_A', '3GUI_A', '3GUJ_A', '3GUN_A', '3GUN_B', '3GUP_A', '3GUP_B', '3HH3_A', '3HH4_A', '3HH5_A', '3HH6_A', '3HT6_A', '3HT7_A', '3HT8_A', '3HTB_A', '3HTD_A', '3HTF_A', '3HTG_A', '3HU8_A', '3HU9_A', '3HUA_A', '3HUK_A', '3HUQ_A', '3HWL_A', '3L64_A', '3LZM_A', '4LZM_A', '5LZM_A', '6LZM_A', '7LZM_A',
    ]

ref_seq = ['MET', 'ASN', 'ILE', 'PHE', 'GLU', 'MET', 'LEU', 'ARG', 'ILE', 'ASP', 'GLU', 'GLY', 'LEU', 'ARG', 'LEU', 'LYS', 'ILE', 'TYR', 'LYS', 'ASP', 'THR', 'GLU', 'GLY', 'TYR', 'TYR', 'THR', 'ILE', 'GLY', 'ILE', 'GLY', 'HIS', 'LEU', 'LEU', 'THR', 'LYS', 'SER', 'PRO', 'SER', 'LEU', 'ASN', 'ALA', 'ALA', 'LYS', 'SER', 'GLU', 'LEU', 'ASP', 'LYS', 'ALA', 'ILE', 'GLY', 'ARG', 'ASN', 'CYS', 'ASN', 'GLY', 'VAL', 'ILE', 'THR', 'LYS', 'ASP', 'GLU', 'ALA', 'GLU', 'LYS', 'LEU', 'PHE', 'ASN', 'GLN', 'ASP', 'VAL', 'ASP', 'ALA', 'ALA', 'VAL', 'ARG', 'GLY', 'ILE', 'LEU', 'ARG', 'ASN', 'ALA', 'LYS', 'LEU', 'LYS', 'PRO', 'VAL', 'TYR', 'ASP', 'SER', 'LEU', 'ASP', 'ALA', 'VAL', 'ARG', 'ARG', 'CYS', 'ALA', 'LEU', 'ILE', 'ASN', 'MET', 'VAL', 'PHE', 'GLN', 'MET', 'GLY', 'GLU', 'THR', 'GLY', 'VAL', 'ALA', 'GLY', 'PHE', 'THR', 'ASN', 'SER', 'LEU', 'ARG', 'MET', 'LEU', 'GLN', 'GLN', 'LYS', 'ARG', 'TRP', 'ASP', 'GLU', 'ALA', 'ALA', 'VAL', 'ASN', 'LEU', 'ALA', 'LYS', 'SER', 'ARG', 'TRP', 'TYR', 'ASN', 'GLN', 'THR', 'PRO', 'ASN', 'ARG', 'ALA', 'LYS', 'ARG', 'VAL', 'ILE', 'THR', 'THR', 'PHE', 'ARG', 'THR', 'GLY', 'THR', 'TRP', 'ASP', 'ALA', 'TYR', 'LYS', 'ASN', 'LEU']

d_res_nos = {}
for pdb in l_pdbs:
    pdb = pdb[:4].lower()
    fd = open('/data/pdb-v3.2/%s/pdb%s.ent' %(pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()
    l_res_nos = []
    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            altloc = line[16]
            res_no = int(line[22:26])
            res_name = line[17:20]
            if res_name != ref_seq[res_no-1]:
                continue
            atom_name = line[12:16].strip()
            if altloc != ' ':
                if res_no in l_res_nos:
                    continue
                if not res_no in d_res_nos.keys():
                    d_res_nos[res_no] = 0
                d_res_nos[res_no] += 1
                l_res_nos += [res_no]
                if res_no in [109,110,111,112,113,]:
                    print pdb, res_no, atom_name, altloc

for res_no in range(165):
    if not res_no in d_res_nos.keys():
        continue
    if d_res_nos[res_no] > 30:
        print res_no, d_res_nos[res_no]

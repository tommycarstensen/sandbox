import os

pdbpath = '/oxygenase_local/data/pdb/'

set_465 = set()

subdirs = os.listdir(pdbpath)
subdirs.sort()
for subdir in subdirs:
    files = os.listdir(pdbpath+subdir)
    print subdirs.index(subdir), len(subdirs), subdir, len(files)
    for file in files:
        pdb = file[3:7]
##        if pdb in [
##            '148l','363d','1abj','2aig','3aig','1ao9','1apu','5apr','2asj',
##            '1at4','1ay6','1b5g','1ba8','1bda','1btw','1cmx','4cpa','1cpi',
##            '1dan','1dit','2dlc','1dva','1dwe','1efm','1epq','1eqq','4est',
##            '1fph','2gbh','3gch','1gct','2gct','1ghb','2gmt','2grb','1gwb',
##            '3hat','1hai','1hap','1hao','1hef','1heg','1hlt','2hpq','2hpp',
##            '1htt','7hvp','2iff','2ja5','2jnd','2neo','1nrr','1nu9','1nu7',
##            '4pad','1ppk','1ppb','1ppc','3prk','1pts','1ry1','1sa1','1skw',
##            '1sl0','1tk5','1tlp','1tmu','1wu1','1x9s','1x9w','1xkf',
##            ]:
##            continue
        fd = open('%s%s/%s' %(pdbpath,subdir,file),'r')
        lines = fd.readlines()
        fd.close()
        parse = False
        for line in lines:
            if line[:10] == 'REMARK 465':
                if line.strip() == 'REMARK 465   M RES C SSSEQI':
                    parse = True
                    continue
                if parse == True:
                    res_name = line[15:18].strip()
                    set_465 |= set([res_name])
                    if res_name not in [
                        ## MODRES
                        'IU','TPN','GPN','U33',
                        'MSE','0CS', '143', '175', '1AB', '1LU', '1PA', '1TQ', '1TY', '23F', '23S', '2DF', '2MR', '2MT', '2PP', '2TY', '32S', '32T', '3AH', '3DR', '3MD', '3PA', '3TY', '4BA', '4DP', '4F3', '4FB', '4FW', '4HT', '4IN', '5CS', '5ZA', '6CL', '6CW', '9AC', 'AA4', 'AA6', 'AAB', 'AAR', 'AB7', 'ABA', 'ABU', 'ACA', 'ACB', 'ACE', 'ACL', 'ACY', 'ADD', 'AEA', 'AEI', 'AFA', 'AGM', 'AGT', 'AHB', 'AHH', 'AHO', 'AHP', 'AHS', 'AHT', 'AIB', 'AKL', 'ALC', 'ALG', 'ALM', 'ALN', 'ALO', 'ALQ', 'ALS', 'ALT', 'ALY', 'AMU', 'ANI', 'APE', 'APH', 'API', 'APK', 'APM', 'APP', 'APY', 'AR2', 'AR4', 'ARM', 'ARO', 'ASA', 'ASB', 'ASI', 'ASK', 'ASM', 'ASU', 'AYA', 'AYG', 'AZK', 'B2A', 'B2F', 'B2I', 'B2V', 'B3A', 'B3D', 'B3E', 'B3K', 'B3S', 'B3X', 'B3Y', 'BAL', 'BBC', 'BCC', 'BCS', 'BCX', 'BE2', 'BFD', 'BHD', 'BIC', 'BIF', 'BLE', 'BLY', 'BMT', 'BNA', 'BNN', 'BNO', 'BOC', 'BOR', 'BP4', 'BPE', 'BTR', 'BUC', 'BUG', 'BZP', 'C12', 'C1X', 'C3Y', 'C5C', 'C6C', 'C99', 'CAB', 'CAF', 'CAL', 'CAS', 'CAV', 'CBG', 'CCS', 'CCY', 'CDE', 'CGN', 'CGU', 'CH3', 'CH6', 'CH7', 'CHF', 'CHG', 'CHP', 'CHS', 'CLB', 'CLD', 'CLE', 'CLG', 'CLH', 'CLT', 'CLV', 'CME', 'CMH', 'CMT', 'CPC', 'CPI', 'CPV', 'CPY', 'CQR', 'CR0', 'CR2', 'CR5', 'CR7', 'CR8', 'CRG', 'CRK', 'CRO', 'CRQ', 'CRU', 'CRW', 'CRX', 'CS3', 'CS4', 'CSA', 'CSB', 'CSD', 'CSE', 'CSH', 'CSI', 'CSO', 'CSP', 'CSR', 'CSS', 'CSU', 'CSW', 'CSX', 'CSY', 'CSZ', 'CTH', 'CUC', 'CWR', 'CXM', 'CY0', 'CY1', 'CY3', 'CY4', 'CYD', 'CYF', 'CYG', 'CYJ', 'CYQ', 'CYR', 'CZ2', 'CZZ', 'D3', 'D4P', 'DAB', 'DAH', 'DAL', 'DAM', 'DAR', 'DAS', 'DBU', 'DBY', 'DBZ', 'DCI', 'DCL', 'DCY', 'DCZ', 'DDE', 'DDX', 'DFI', 'DFO', 'DFT', 'DG', 'DGL', 'DGN', 'DHA', 'DHI', 'DHL', 'DIL', 'DIP', 'DIV', 'DIX', 'DIY', 'DLE', 'DLS', 'DLY', 'DMH', 'DMK', 'DMT', 'DNE', 'DNG', 'DNL', 'DNM', 'DOA', 'DPH', 'DPL', 'DPN', 'DPP', 'DPR', 'DPY', 'DRP', 'DRZ', 'DSE', 'DSG', 'DSN', 'DSY', 'DTG', 'DTR', 'DTY', 'DVA', 'DYG', 'EFC', 'EMR', 'EPO', 'ESC', 'ESD', 'ETA', 'ETO', 'FAG', 'FBA', 'FBE', 'FCL', 'FGL', 'FGP', 'FHL', 'FKI', 'FLE', 'FLT', 'FME', 'FOE', 'FOG', 'FOR', 'FRD', 'FRF', 'FTR', 'FTY', 'GHG', 'GHP', 'GL3', 'GLH', 'GLM', 'GLZ', 'GMA', 'GPL', 'GT9', 'GTH', 'GVL', 'GYC', 'GYS', 'H5M', 'HAQ', 'HFA', 'HHK', 'HIA', 'HIC', 'HIP', 'HIQ', 'HLU', 'HMA', 'HMB', 'HMF', 'HMI', 'HMR', 'HOA', 'HPD', 'HPE', 'HPH', 'HPQ', 'HQU', 'HRG', 'HSE', 'HSL', 'HSO', 'HTI', 'HTR', 'HV7', 'HV8', 'HY3', 'HYP', 'IAM', 'IAS', 'IGL', 'IIC', 'IIL', 'ILG', 'ILX', 'IML', 'IPG', 'IPN', 'ISO', 'IVA', 'IYR', 'KCX', 'KI2', 'KOR', 'KPH', 'KYN', 'LAC', 'LAL', 'LCX', 'LED', 'LEF', 'LET', 'LLP', 'LLY', 'LME', 'LML', 'LOL', 'LOV', 'LPD', 'LPL', 'LSO', 'LTA', 'LYM', 'LYN', 'LYP', 'LYX', 'LYZ', 'M3L', 'MAA', 'MAI', 'MAZ', 'MBQ', 'MBZ', 'MC1', 'MCL', 'MCS', 'MDH', 'MDO', 'MDP', 'MEA', 'MEG', 'MEN', 'MEU', 'MFC', 'MGG', 'MGN', 'MHL', 'MHO', 'MHS', 'MIS', 'MLE', 'MLL', 'MLY', 'MLZ', 'MME', 'MN1', 'MN2', 'MN7', 'MN8', 'MNL', 'MNV', 'MOP', 'MOR', 'MPQ', 'MPR', 'MPT', 'MRM', 'MSA', 'MSO', 'MSU', 'MTY', 'MVA', 'MYR', 'NAL', 'NAM', 'NBQ', 'NC1', 'NCB', 'NEP', 'NFA', 'NFP', 'NIT', 'NIY', 'NLE', 'NLN', 'NLO', 'NMC', 'NME', 'NOA', 'NOR', 'NP3', 'NPH', 'NPT', 'NRI', 'NRQ', 'NSP', 'NVA', 'NYC', 'NZH', 'OAS', 'OBS', 'OCS', 'OCY', 'ODA', 'ODS', 'OHS', 'OLT', 'OME', 'OMT', 'OPH', 'OPN', 'OPR', 'ORN', 'ORP', 'ORQ', 'OSE', 'OXX', 'P1L', 'P2Y', 'PAA', 'PAQ', 'PAT', 'PBB', 'PBF', 'PCA', 'PCH', 'PCS', 'PDI', 'PEA', 'PEC', 'PED', 'PF5', 'PFF', 'PG1', 'PG9', 'PGL', 'PHA', 'PHD', 'PHI', 'PHL', 'PHM', 'PHQ', 'PIA', 'PIP', 'PIV', 'PLE', 'PLP', 'PM3', 'PN2', 'PO0', 'POM', 'PPH', 'PPN', 'PR0', 'PR3', 'PRR', 'PRS', 'PS0', 'PSA', 'PSE', 'PTA', 'PTH', 'PTL', 'PTM', 'PTR', 'PVA', 'PXZ', 'PYA', 'PYC', 'PYP', 'PYR', 'PYT', 'PYX', 'PYY', 'QLG', 'QNC', 'QUI', 'R1A', 'R1B', 'R1F', 'R7A', 'RC7', 'RHS', 'RNG', 'S02', 'S1H', 'SAC', 'SAR', 'SBD', 'SBG', 'SBL', 'SC2', 'SCH', 'SCS', 'SCY', 'SDP', 'SEB', 'SEC', 'SEL', 'SEP', 'SET', 'SGB', 'SGR', 'SHC', 'SHP', 'SIN', 'SLE', 'SLZ', 'SMC', 'SME', 'SMF', 'SNC', 'SNN', 'SOC', 'SOY', 'STA', 'SUI', 'SUN', 'SVA', 'SVV', 'SVX', 'SVY', 'SVZ', 'SXE', 'TBG', 'TBM', 'TCQ', 'TCR', 'TEE', 'TFA', 'THC', 'THO', 'THX', 'TIH', 'TLX', 'TMD', 'TNB', 'TOX', 'TPL', 'TPO', 'TPQ', 'TQQ', 'TRF', 'TRJ', 'TRN', 'TRO', 'TRQ', 'TRW', 'TRX', 'TTQ', 'TTS', 'TYB', 'TYC', 'TYI', 'TYN', 'TYO', 'TYQ', 'TYS', 'TYT', 'TYX', 'TYY', 'TYZ', 'UMA', 'VAD', 'VAF', 'VDL', 'VLL', 'VOL', 'X9Q', 'XAO', 'XXY', 'XYG', 'YCM', 'YOF', 'YRR',
                        ## nucleotides
                        'DU','DC','DG','DA','C','G','A','DT','U','T',
                        'UNK',
                        ## standard residues
                        'GLX',
                        'CYS', 'GLN', 'ILE', 'SER', 'VAL', 'MET', 'ASN', 'PRO', 'LYS', 'THR', 'PHE', 'ALA', 'HIS', 'GLY', 'ASP', 'LEU', 'ARG', 'TRP', 'GLU', 'TYR',
                        ]:
                        print file, res_name
                        stop
    print set_465
set_465 = set(['A', 'C', 'ILE', 'DT', 'GLN', 'DG', 'IU', 'THR', 'DC', 'G', 'GLY', 'ASP', 'U', 'TRP', 'MSE', 'DA', 'GLU', 'CYS', 'PTR', 'HIS', 'SER', 'LYS', 'PRO', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'VAL', 'ASN', 'TYR'])

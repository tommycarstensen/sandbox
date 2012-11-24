import os, sys
sys.path.append('svn/tc_sandbox/pdb')
import parse_mmCIF
import matthews_coefficient
sys.path.append('svn/tc_sandbox/misc')
import statistics

s_pdbs = '''
1A87
1AL1
1B87
1EME
1EMF
1FY7
1GRQ
1GRR
1J6U
1JKY
1KIB
1MJ9
1MJA
1MJB
1NHY
1OIA
1ORY
1QAX
1QAY
1QHN
1QHS
1QHX
1QHY
1R31
1R7I
1RKC
1SJV
1T02
1U61
1VHG
1VHZ
1VIP
1XX4
1YLI
2AO7
2EMD
2FD7
2FP3
2FZL
2G2D
2GWW
2IAK
2JA9
2PMP
2XBZ
2ZDB
3BSD
3CDG
3CII
3E4H
3EYL
3FOW
3HON
3IO2
3MMS
'''

def main():

    d_MV = {}

    path = '/data/mmCIF'
    l_dn = os.listdir(path)
    l_dn.sort()
    for dn in l_dn:
        if dn == 'mmCIF.py':
            continue
        if dn < sys.argv[-2]:
            continue
        if dn > sys.argv[-1]:
            continue
        l_fn = os.listdir('%s/%s' %(path,dn))
        for fn in l_fn:
            pdb = fn[:4]
##            if pdb.upper() not in s_pdbs:
##                continue
            d_mmCIF = parse_mmCIF.main(
                pdb,
                d_breaks = {'_exptl.method':'SOLUTION NMR'},
                l_data_categories = [
                    '_cell','_entity','_exptl','_exptl_crystal',
                    '_entity_poly',
                    '_symmetry',
                    ## virus
                    '_pdbx_struct_assembly',
                    ## split structure
                    '_pdbx_database_related',
                    ],
                )

            ## x-ray structure
            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                continue

            ## polymer present
            if not '_entity_poly.type' in d_mmCIF.keys():
                continue

            ## only polymer present is protein
            if d_mmCIF['_entity_poly.type'] != len(d_mmCIF['_entity_poly.type'])*['polypeptide(L)']:
                continue

            if not '_pdbx_struct_assembly.oligomeric_count' in d_mmCIF.keys():
                continue

            if d_mmCIF['_pdbx_struct_assembly.oligomeric_count'] == len(d_mmCIF['_pdbx_struct_assembly.oligomeric_count'])*['?']:
                continue
            
            ## virus
            if int(d_mmCIF['_pdbx_struct_assembly.oligomeric_count'][0]) % 60 == 0:
                continue

            ## not monomer
            if d_mmCIF['_pdbx_struct_assembly.oligomeric_count'] != len(d_mmCIF['_pdbx_struct_assembly.oligomeric_count'])*['1']:
                continue

            ## split structure
            if '_pdbx_database_related' in d_mmCIF.keys():
                if 'split' in d_mmCIF['_pdbx_database_related']:
                    continue
                if 'SPLIT' in d_mmCIF['_pdbx_database_related']:
                    print pdb
                    stop

            if not '_cell.Z_PDB' in d_mmCIF.keys():
                continue

            if pdb in [
                ## treshold
                '1e54','1e9i',
                ## difference between calculated MV and MV in mmCIF
                '3eiq',
                ## The crystals diffracted to 1.7Angstrom and appeared to be I centered tetragonal with
                ## unit cell dimension a=198.42Angstrom and c=396.6Angstrom, however the data only merged successfully in P1
                ## unit cell a=196.61 b=196.48 c=240.63 alpha=65.91 beta=65.91 gamma=90.01.
                ## Toscana has published with Hellinga...
                '2cjf','2bt4',
                ]:
                continue

##            if not ''.join(d_mmCIF['_symmetry.space_group_name_H-M']) in [
##                'P 1','P 43 21 2','P 21 3','P 42 3 2','C 1 2 1','F 2 3','P 64 2 2','H 3',
##                ]:
##                continue ## tmp!!!

            a = float(d_mmCIF['_cell.length_a'][0])
            b = float(d_mmCIF['_cell.length_b'][0])
            c = float(d_mmCIF['_cell.length_c'][0])
            alpha = float(d_mmCIF['_cell.angle_alpha'][0])
            beta = float(d_mmCIF['_cell.angle_beta'][0])
            gamma = float(d_mmCIF['_cell.angle_gamma'][0])
            Z = int(d_mmCIF['_cell.Z_PDB'][0])
            mw = 0
            for i in range(len(d_mmCIF['_entity.id'])):
##                if d_mmCIF['_entity.type'][i] == 'polymer':
                    s = d_mmCIF['_entity.formula_weight'][i]
                    ## unknown ligand
                    if s == '?':
                        continue
                    mw += float(s)

            MV = matthews_coefficient.main(a,b,c,alpha,beta,gamma,mw,Z)

            spacegroup = ''.join(d_mmCIF['_symmetry.space_group_name_H-M'])

            if spacegroup not in [
                'F 4 3 2',
                'F 41 3 2',
                'I 41 3 2',
                ]:
                continue ## tmp!!!

            if MV > 10:
                print pdb
                print 'mw', mw
                print 'MV', MV, d_mmCIF['_exptl_crystal.density_Matthews']
                print 'Z', Z
                import math
                alpha *= math.pi/180.
                beta *= math.pi/180.
                gamma *= math.pi/180.
                V = a*b*c*math.sqrt(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
                print 'V', V
                continue
                stop_treshold
                stop
            if '_exptl_crystal.density_Matthews' in d_mmCIF.keys():
                if d_mmCIF['_exptl_crystal.density_Matthews'] not in [['?'],len(d_mmCIF['_exptl_crystal.density_Matthews'])*['?'],]:
                    if abs(MV-float(d_mmCIF['_exptl_crystal.density_Matthews'][0])) > 1:
                        print 'MV', MV
                        print 'MV', d_mmCIF['_exptl_crystal.density_Matthews']
                        print 'mw', mw
                        print 'Z', Z
                        continue
                        stop_difference


            if not spacegroup in d_MV.keys():
                d_MV[spacegroup] = []
            d_MV[spacegroup] += [MV]

            print pdb, round(MV,2), spacegroup

##    fd = open('MV_v_spacegroup.txt','w')
##    fd.write(str(d_MV))
##    fd.close()

    l = ['# MV_average MV_stddev n spacegroup\n']
    for spacegroup in d_MV.keys():
        l_MV = d_MV[spacegroup]
        if len(l_MV) <= 1:
            continue
        average, stddev = statistics.do_stddev(l_MV)
        average, stderr = statistics.do_stderr(l_MV)
##        l += ['%s %s %s %s\n' %(average,stddev,len(l_MV),spacegroup,)]
        l += ['%s %s %s %s\n' %(average,stderr,len(l_MV),spacegroup,)]

    fd = open('MV_v_spacegroup.txt','w')
    fd.writelines(l)
    fd.close()

    return

if __name__ == '__main__':
    main()

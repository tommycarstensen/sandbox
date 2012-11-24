import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb')
import matthews_coefficient, parse_mmCIF


for pdb in [
    '2hhb','1hho','1hv4',
    '3hl9','3hlb','3hlc','3hld','3hle','3hlf','3hlg',
    ]:

    d_mmCIF = parse_mmCIF.main(pdb)

    a = float(d_mmCIF['_cell.length_a'][0])
    b = float(d_mmCIF['_cell.length_b'][0])
    c = float(d_mmCIF['_cell.length_c'][0])
    alpha = float(d_mmCIF['_cell.angle_alpha'][0])
    beta = float(d_mmCIF['_cell.angle_beta'][0])
    gamma = float(d_mmCIF['_cell.angle_gamma'][0])
    Z = int(d_mmCIF['_cell.Z_PDB'][0]) ## number of polymers in unit cell
    mw = 0
    for i in range(len(d_mmCIF['_entity.id'])):
        if d_mmCIF['_entity.type'][i] == 'polymer':
            mw += float(d_mmCIF['_entity.formula_weight'][i])
    MV = matthews_coefficient.main(a,b,c,alpha,beta,gamma,mw,Z,)

    print pdb, MV

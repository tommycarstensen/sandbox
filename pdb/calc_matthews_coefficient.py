## Tommy Carstensen, UCD, February 2010
## Matthews 1968,1976

import math, sys
import parse_mmCIF

def main(pdb):

    ## speed up by not reading atom section...
    d_mmCIF = parse_mmCIF.main(pdb)

    a = float(d_mmCIF['_cell.length_a'][0])
    b = float(d_mmCIF['_cell.length_b'][0])
    c = float(d_mmCIF['_cell.length_c'][0])
    alpha = float(d_mmCIF['_cell.angle_alpha'][0])
    beta = float(d_mmCIF['_cell.angle_beta'][0])
    gamma = float(d_mmCIF['_cell.angle_gamma'][0])
    Z = int(d_mmCIF['_cell.Z_PDB'][0])
    mw = 0
    for i in range(len(d_mmCIF['_entity.id'])):
        if d_mmCIF['_entity.type'][i] == 'polymer':
            mw += float(d_mmCIF['_entity.formula_weight'][i])
    VM = calc(a,b,c,alpha,beta,gamma,mw,Z,)
    print pdb, VM

    return VM

def calc(a,b,c,alpha,beta,gamma,mw,Z,):

    alpha *= math.pi/180.
    beta *= math.pi/180.
    gamma *= math.pi/180.
    ## unit cell volume
    V = volume = a*b*c*math.sqrt(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))

    VM = (V/mw)/Z

    return VM

if __name__ == '__main__':

    if '-pdb' in sys.argv:
        pdb = sys.argv[-1]
        main(pdb)

    else:

        for pdb in [
    '1woy','2d5b',
            ]:

            main(pdb)

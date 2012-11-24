import sys
sys.path.append('/home/tc/svn/tc_sandbox/pdb')
import parse_mmCIF, mmCIF2coords
sys.path.append('/home/tc/svn/GoodVibes')
import NMA, visualization

d_mmCIF = parse_mmCIF.main('2lzm',)
d_coords, l_coords_alpha = mmCIF2coords.main('2lzm',d_mmCIF)

cutoff = 10
matrix_hessian = NMA.hessian_calculation(l_coords_alpha, cutoff)
eigenvectors, eigenvalues = NMA.diagonalize_hessian(matrix_hessian,)
visualization.vmd_trajectory('2lzm',l_coords_alpha,eigenvectors)

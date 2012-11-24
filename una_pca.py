def main():

    import sets, LinearAlgebra, Numeric
    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/rmsf_nmr.py')
    import cumulative_overlap,rmsd_two_conformations,pca

    pdb1 = '1atp'
    chain1 = 'E'
    pdb2 = '1j3h'
    chain2 = 'A'

    d_coordinates1,l_coordinates1 = rmsd_two_conformations.parse_coordinates(pdb1,chain1)
    d_coordinates2,l_coordinates2 = rmsd_two_conformations.parse_coordinates(pdb2,chain2)

    vector_difference = []
    d_differences = {}
    for res_no in d_coordinates1[chain1].keys():
        if res_no not in d_coordinates2[chain2].keys():
            continue
        iCode = ' '
        c1 = d_coordinates1[chain1][res_no][iCode]['CA']
        c2 = d_coordinates2[chain2][res_no][iCode]['CA']
        d_differences[res_no] = c2-c1
        if res_no in range(15,318+1):
            vector_difference += [c2[0]-c1[0],c2[1]-c1[1],c2[2]-c1[2],]
    vector_difference = Numeric.array(vector_difference)
    l_resnos = list(sets.Set(range(15,350))-sets.Set(d_differences.keys()))
    l_resnos.sort()
    l_resnos.reverse()

    matrix = parse_covariance_matrix(l_resnos)

    eigenvalues, eigenvectors = pca.diagonalization(matrix)

    eigenvectors_transposed = cumulative_overlap.transpose_rows_and_columns(eigenvectors)

####    matrix = parse_principal_components()

    print len(eigenvectors_transposed), len(vector_difference)
    l_contributions = LinearAlgebra.solve_linear_equations(eigenvectors_transposed,vector_difference)
    print l_contributions

    fd = open('contrib.tmp','w')
    fd.close()
    for i in range(len(l_contributions)):
        fd = open('contrib.tmp','a')
        fd.write('%s %s\n' %(i+1,l_contributions[i]))
        fd.close()

    fo = 'cumoverlap.tmp'
    mode_min = 0
    cumulative_overlap.calculate_cumulated_overlap(fo,l_contributions,vector_difference,eigenvectors,mode_min,)

    return


def parse_covariance_matrix(l_resnos):

    import Numeric

    fd = open('covar.dat','r')
    l = fd.read().split()
    fd.close()
    n = 3*(350-14)

    matrix = Numeric.zeros((3*(318-14),3*(318-14)),typecode='d')
    for i in range(len(l)):
        col = (i % n)
        row = (i-col)/n
        if row > 3*(318-14)-1 or col > 3*(318-14)-1:
            continue
        matrix[row][col] = float(l[i])

    return matrix


def parse_principal_components():

    import Numeric

    fd = open('eigcomp.xvg','r')
    lines = fd.readlines()
    fd.close()

    n = 336
    matrix = Numeric.zeros((3*n,3*n),typecode='d')
##    matrix = Numeric.zeros((1001,1001),typecode='d')

    d_vectors = {}
    for line in lines:
        if line[0] == '@':
            axis = -1
            if line.split()[1:3] == ['yaxis','label',]:
                index1 = line.index('"vec')+4
                index2 = len(line)-2
                i = int(line[index1:index2])-1
                d_vectors[i] = {0:[],1:[],2:[],}
            continue
        elif line[0] == '&':
            axis += 1
        elif axis == -1:
            continue
        else:
            res_no = int(float(line.split()[0]))
            matrix[i][3*(res_no-1)+axis] = float(line.split()[1])
##            d_vectors[i][axis] += [float(line.split()[1])]

    return matrix


if __name__ == '__main__':
    main()

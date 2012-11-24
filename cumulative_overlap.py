def main():

    import Numeric,LinearAlgebra
    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/rmsf_nmr.py')
    import rmsf_nmr

    pdb = '1e8l'
    chain = 'A'
    model1 = 1
    model2 = 2

    d_coordinates = rmsf_nmr.parse_coordinates(pdb,chain,)

    vector_difference = calculate_difference_vector(d_coordinates,model1,model2,)

    l_coordinates = d_coordinates[1]
    matrix_hessian = rmsf_nmr.calculate_hessian_matrix(l_coordinates)

    l_eigenvectors = rmsf_nmr.calculate_eigenvectors(matrix_hessian)

    l_eigenvectors_transposed = transpose_rows_and_columns(l_eigenvectors)

    vector_difference = Numeric.array(vector_difference)
    l_contributions = LinearAlgebra.solve_linear_equations(l_eigenvectors_transposed,vector_difference)

    fd = open('contrib_%s_%s.tmp' %(model1,model2),'w')
    fd.close()
    for i in range(len(l_contributions)):
##        print i, vector[i]
        fd = open('contrib_%s_%s.tmp' %(model1,model2),'a')
        fd.write('%s %s\n' %(i+1,l_contributions[i]))
        fd.close()

    fo = 'cumoverlap_%s_%s.tmp' %(model1,model2)
    mode_min = 6
    calculate_cumulated_overlap(fo,l_contributions,vector_difference,l_eigenvectors,mode_min)

    return


def calculate_cumulated_overlap(fo,l_contributions,vector_difference,l_eigenvectors,mode_min,):

    fd = open(fo,'w')
    fd.close()
    vector_cumulated = []
    for i in range(len(l_contributions)):
        vector_cumulated += [0]
    ## loop over mode
    for i in range(len(l_contributions)):
        if i < mode_min:
            overlap = 0
        else:
            ## loop over coordinate
            for j in range(len(vector_cumulated)):
                vector_cumulated[j] += l_contributions[i]*l_eigenvectors[i][j]
            overlap = cosangle(vector_cumulated,vector_difference)
        print i,overlap
        fd = open(fo,'a')
        fd.write('%s %s\n' %(i+1,overlap))
        fd.close()

    return


def cosangle(v1,v2):

    import math

    if len(v1) != len(v2):
        print len(v1), len(v2)
        stop

    numerator = 0
    for i in range(len(v1)):
        numerator += v1[i]*v2[i]

    denominator1 = 0
    denominator2 = 0
    for i in range(len(v1)):
        denominator1 += v1[i]*v1[i]
        denominator2 += v2[i]*v2[i]
    denominator = math.sqrt(denominator1*denominator2)

    cosang = numerator / denominator

    return cosang
    

def transpose_rows_and_columns(l_eigenvectors):

    import Numeric

    n = len(l_eigenvectors)

    l_eigenvectors2 = Numeric.zeros((n,n),typecode='d')

    for i in range(n):
        for j in range(n):
            l_eigenvectors2[i][j] = l_eigenvectors[j][i]

    return l_eigenvectors2


def calculate_difference_vector(d_coordinates,model1,model2):

    l_vectors = []

    n_coordinates = len(d_coordinates[1])

    for i in range(n_coordinates):
        c2 = d_coordinates[model2][i]
        c1 = d_coordinates[model1][i]
        l_vectors += [c2[0]-c1[0],c2[1]-c1[1],c2[2]-c1[2],]

    return l_vectors
    

##def gauss_elimination(matrix,vector):
##
##    n = len(matrix)
##
##    for col in range(n):
##        if col % 10 == 0:
##            print 'col', col
##        for row in range(col+1,n):
##            ## addition
##            multiple = float(matrix[row][col])/float(matrix[col][col])
##            for i in range(col,n):
##                matrix[row][i] -= multiple*matrix[col][i]
##            vector[row] -= multiple*vector[col]
##            ## multiplication
##            multiple = float(matrix[row][col+1])
##            if multiple != 0:
##                for i in range(col,n):
##                    matrix[row][i] /= multiple
##                vector[row] /= multiple
##
##    for row in range(n-1-1,-1,-1):
####        print 'row', row
##        for col in range(row+1,n):
##            ## addition
##            multiple = float(matrix[row][col])
####            for i in range(col,col-1,-1):
####                matrix[row][i] -= multiple*
##            vector[row] -= multiple*vector[col]
##
##    return vector

if __name__ == '__main__':
    main()

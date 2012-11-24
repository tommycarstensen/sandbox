##!/software/bin/python2.3
##
##$Id$
##
##Tommy Carstensen, University College Dublin, 2008

import os

def main(
    m = [
        [.69,-1.31,.39,.09,1.29,.49,.19,-.81,-.31,-.71,],
        [.49,-1.21,.99,.29,1.09,.79,-.31,-.81,-.31,-1.01,],
        ]
    ):

    import numpy

    matrix = numpy.zeros((len(m),len(m)))
    for i in range(len(m)):
        l1 = m[i]
        n = len(l1)
        average1 = sum(l1)/n
        for j in range(i,len(m)):
            l2 = m[j]
            if n != len(l2):
                raise 'data lists are of different length'
            average2 = sum(l2)/n

            SS = 0
            for k in range(n):
                SS += (l1[k]-average1)*(l2[k]-average2)
            covar = SS/(n-1)
            matrix[i][j] = covar
            matrix[j][i] = covar
        
    eigenvalues, eigenvectors = diagonalization(matrix)
    print 'eigenvalues', eigenvalues
    print 'eigenvectors', eigenvectors

    if len(m) == 2:
        gnuplot(m,eigenvectors)

    return eigenvalues, eigenvectors


def gnuplot(m,eigenvectors):

    lines = []
    n = len(m[0])
    for i in range(n):
        f1 = m[0][i]
        f2 = m[1][i]
        lines += ['%s %s\n' %(f1,f2)]
    fd = open('test.gpd','w')
    fd.writelines(lines)
    fd.close()

    average1 = sum(m[0])/n
    average2 = sum(m[1])/n
    b1 = eigenvectors[0][0]/eigenvectors[0][1]
    b2 = eigenvectors[1][0]/eigenvectors[1][1]
    a1 = average1-b1*average2
    a2 = average1-b2*average2

    lines = [
        'set terminal postscript eps enhanced color "Helvetica" 12\n',
        'set output "gnuplot.ps"\n',
        'set size square\n',
        'f(x) = %f+%f*x\n' %(a1,b1),
        'g(x) = %f+%f*x\n' %(a2,b2),
        'plot "test.gpd" ps 1 pt 7, f(x), g(x)',
        ]
    fd = open('test.gps','w')
    fd.writelines(lines)
    fd.close()

    os.system('/software/bin/gnuplot test.gps')

    return


def diagonalization(matrix):

    print 'diagonalizing %sx%s matrix' %(len(matrix),len(matrix))

    import LinearAlgebra

    eigen_tuple = LinearAlgebra.eigenvectors(matrix)
    ## parse eigenvalues and eigenvectors
    eigenvalues = list(eigen_tuple[0])
    eigenvectors = list(eigen_tuple[1])
    ## organize eigenvalues and eigenvectors in list
    eigen_list = zip(eigenvalues, eigenvectors)
    ## sort list
    eigen_list.sort()
    ## reverse list
    eigen_list.reverse()
    ## parse sorted eigenvalues and eigenvectors
    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]

    return eigenvalues, eigenvectors
    

if __name__ == '__main__':
    main()

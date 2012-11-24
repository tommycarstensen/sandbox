def hessian_calculation(l_coordinates, cutoff, verbose = True):

    '''This function calculates the Hessian matrix. At least its supposed to...'''

    N = len(l_coordinates)

    if verbose == True:
        print 'calculating the hessian/second-order derivate Kirchoff/Laplacian matrix'
    
    import math, Numeric, numpy

    cutoff_sq = cutoff**2
    sigmoid_slope = 1
    sigmoid_factor = 1 ## default

    matrix_hessian = Numeric.zeros((3*N,3*N), typecode='d')
    matrix_hessian = numpy.zeros((3*N,3*N))

    for row_sup in range(N):
        for col_sup in range(N):

            if col_sup > row_sup:
                xi = l_coordinates[row_sup][0]
                xj = l_coordinates[col_sup][0]
                yi = l_coordinates[row_sup][1]
                yj = l_coordinates[col_sup][1]
                zi = l_coordinates[row_sup][2]
                zj = l_coordinates[col_sup][2]
                x = xj-xi
                y = yj-yi
                z = zj-zi
                dist_sq = x**2+y**2+z**2

##                ## 1) sharp cutoff
##                if dist_sq > cutoff_sq:
##                    continue
                ## 2) sigmoid cutoff
                sigmoid_factor = 1. / ( 1. + math.exp( sigmoid_slope*(math.sqrt(dist_sq)-cutoff) ) )

                vector = [x,y,z]
                for row_sub in range(3):
                    for col_sub in range(3):

                        if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                            value = sigmoid_factor*-vector[row_sub]*vector[col_sub]/dist_sq
                            matrix_hessian[3*row_sup+row_sub,3*col_sup+col_sub] = value ##upper super off-diagonal; xixj, xiyj, xizj, yiyj, yizj, zizj
                            matrix_hessian[3*col_sup+col_sub,3*row_sup+row_sub] = value ##lower super off-diagonal; xjxi, yjxi, zjxi, yjyi, zjyi, zjzi
                            matrix_hessian[3*row_sup+row_sub,3*row_sup+col_sub] -= value ##super diagonal (row); xixi, xiyi, xizi, yiyi, yizi, zizi
                            matrix_hessian[3*col_sup+col_sub,3*col_sup+row_sub] -= value ##super diagonal (col); xjxj, yjxj, zjxj, yjyj, zjyj, zjzj
                            if col_sub > row_sub: #fill lower subsymmetrical elements
                                matrix_hessian[3*row_sup+col_sub,3*col_sup+row_sub] = value #upper super off-diagonal; yixj, zixj, ziyj
                                matrix_hessian[3*col_sup+row_sub,3*row_sup+col_sub] = value #lower super off-diagonal; xjyi, xjzi, yjzi
                                matrix_hessian[3*row_sup+col_sub,3*row_sup+row_sub] -= value ##super diagonal; yixi, zixi, ziyi
                                matrix_hessian[3*col_sup+row_sub,3*col_sup+col_sub] -= value ##super diagonal; yjxj, zjxj, zjyj

    return matrix_hessian


def diagonalize_hessian(matrix_hessian, verbose = True):

    import math, numpy

    '''Calculates eigenvectors and eigenvalues of a matrix.'''
    if verbose == True:
        print 'calculating eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'

##    t1 = time.clock()
##    eigenvalues, eigenvectors = numpy.linalg.eig(matrix_hessian)
##    print eigenvalues, eigenvectors
##    t2 = time.clock()
##    print 'numpy.eig', t2-t1
##    eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
##    eigen_list = zip(eigenvalues, eigenvectors)
##    eigen_list.sort()
##    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
##    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]
##
##    t1 = time.clock()
##    eigenvaluesh, eigenvectorsh = numpy.linalg.eigh(matrix_hessian)
##    print eigenvaluesh, eigenvectorsh
##    t2 = time.clock()
##    print 'numpy.eigh', t2-t1
##    eigenvectorsh = numpy.transpose(eigenvectorsh) ## transpose when diagonalization with numpy instead of Numeric
##    eigen_listh = zip(eigenvaluesh, eigenvectorsh)
##    eigen_listh.sort()
##    eigenvaluesh = [eigen_listh[eigen][0] for eigen in range(len(eigen_listh))]
##    eigenvectorsh = [eigen_listh[eigen][1] for eigen in range(len(eigen_listh))]
##
##    print 'matrix'
##    print matrix_hessian
##    
##    for x in range(len(eigenvalues)):
##        if round(eigenvalues[x],4) != round(eigenvaluesh[x],4):
##            print '***', x, eigenvalues[x], eigenvaluesh[x]
##
##    for x in range(6,len(eigenvectors)):
##        for y in range(0,len(eigenvectors)):
##            if round(eigenvectors[x][y],6) != round(eigenvectorsh[x][y],6)*(round(eigenvectorsh[x][0],6)/round(eigenvectors[x][0],6)):
##                print '***', x,y, round(eigenvectors[x][y],6), round(eigenvectorsh[x][y],6)*(round(eigenvectorsh[x][0],6)/round(eigenvectors[x][0],6))
##
##    stop

##    import time
##    t1 = time.clock()
    ## diagonalize hessian (symmetric/hermittian) matrix
    eigenvalues, eigenvectors = numpy.linalg.eigh(matrix_hessian)
    ## transpose when diagonalization with numpy instead of Numeric
    eigenvectors = numpy.transpose(eigenvectors)
    ## organize eigenvalues and eigenvectors in list for sorting pairwise
    eigen_list = zip(eigenvalues, eigenvectors)
    ## sort list (already sorted if numpy.linalg.eigh)
    try:
        eigen_list.sort()
    except:
        for i in range(1,len(eigen_list)+1):
            try: eigen_list[:i].sort()
            except: break
        print eigen_list[i]
        print i
        eigen_list[:i].sort()
        stop
    ## parse sorted eigenvalues and eigenvectors
    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]
##    t2 = time.clock()
##    print t2-t1
##    stop

##    ##
##    ## scale modes
##    ##
##
##    ## calculate length of mode 7 (equals 1 when using module linearalgebra)
##    len7 = math.sqrt(sum(numpy.array(eigenvectors[6])*numpy.array(eigenvectors[6])))
##
##    ## loop over modes
##    for i in range(7,len(eigenvalues)-1):
##
##        ## calculate length of mode i
##        leni = math.sqrt(sum(numpy.array(eigenvectors[i])*numpy.array(eigenvectors[i])))
##
##        ## scale length of mode i relative to length of mode 7
##        try:
##            lenfactor = (len7/leni)/(eigenvalues[i]/eigenvalues[6+3])
##        except:
##            print eigenvalues
##            print i, len7, leni, eigenvalues[i], eigenvalues[6+3]
##            stop_float_divison
##        for j in range(len(eigenvectors[i])):
##            eigenvectors[i][j] *= lenfactor

##    ## copy lists of eigenvectors to eigenvectors_combined
##    eigenvectors_combined = []
##    for mode in range(len(eigenvalues)):
##        eigenvectors_combined.append(list(eigenvectors[mode]))
##    ## change mode i to be the sum of modes 7 to i
##    for mode in range(7,len(eigenvalues)):
##        for coordinate in range(len(eigenvalues)):
##            eigenvectors_combined[mode][coordinate] += eigenvectors_combined[mode-1][coordinate]

    if verbose == True:
        print 'calculated eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'

    return eigenvectors, eigenvalues

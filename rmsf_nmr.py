def main():

    pdb = '2bru'
    chain = 'A'

    d_coordinates = parse_coordinates(pdb,chain,)

    l_rmsfs_NMR = calculate_nmr_rmsf(d_coordinates)

    l_coordinates = d_coordinates[1]
    matrix_hessian = calculate_hessian_matrix(l_coordinates)

    l_eigenvectors = calculate_eigenvectors(matrix_hessian)

    for mode in range(6,len(l_eigenvectors)):

        l_rmsfs_NMA = calculate_nma_rmsf(l_eigenvectors[mode])

        r = calculate_correlation(l_rmsfs_NMR,l_rmsfs_NMA)

##        print mode,r

        fd = open('%s.gnuplotdata' %(pdb),'a')
        fd.write('%s %s\n' %(mode+1,r))
        fd.close()

    plot(pdb,)

    return


def plot(pdb,):

    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import gnuplot

    prefix = pdb
    gnuplot.scatter_plot_2d(
        prefix,regression=True,xlabel='mode',ylabel='correlation',ymin=-1,ymax=1,
        )

    return


def calculate_nma_rmsf(eigenvector):

    import math

    l_rmsfs = []
    for i in range(len(eigenvector)/3):
        x = eigenvector[3*i+0]
        y = eigenvector[3*i+1]
        z = eigenvector[3*i+2]
        rmsf = math.sqrt(x**2+y**2+z**2)
        l_rmsfs += [rmsf]

    return l_rmsfs


def calculate_correlation(l1,l2):

    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import statistics
    r = statistics.correlation(l1,l2)

    return r


def calculate_hessian_matrix(l_coordinates):

    import Numeric

    matrix_hessian = Numeric.zeros( (3*len(l_coordinates),3*len(l_coordinates)), typecode='d')

    for row_sup in range(len(l_coordinates)):
        for col_sup in range(row_sup+1,len(l_coordinates)):

            xi = l_coordinates[row_sup][0]
            xj = l_coordinates[col_sup][0]
            yi = l_coordinates[row_sup][1]
            yj = l_coordinates[col_sup][1]
            zi = l_coordinates[row_sup][2]
            zj = l_coordinates[col_sup][2]
            x = xj-xi
            y = yj-yi
            z = zj-zi
            vector = [x,y,z]
            dist_sq = x**2+y**2+z**2
            if dist_sq > 100:
                continue

            for row_sub in range(3):
                for col_sub in range(3):

                    if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                        value = -vector[row_sub]*vector[col_sub]/dist_sq
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


def calculate_eigenvectors(matrix_hessian):

    import LinearAlgebra
    ## diagonalize hessian matrix
    eigen_tuple = LinearAlgebra.Heigenvectors(matrix_hessian)
    ## parse eigenvalues and eigenvectors
    eigenvalues = list(eigen_tuple[0])
    eigenvectors = list(eigen_tuple[1])
    ## organize eigenvalues and eigenvectors in list
    eigen_list = zip(eigenvalues, eigenvectors)
    ## sort list
    eigen_list.sort()
    ## parse sorted eigenvalues and eigenvectors
    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]

    return eigenvectors


def calculate_nmr_rmsf(d_coordinates):

    import math

    n_coordinates = len(d_coordinates[1])
    n_models = len(d_coordinates.keys())

    l_dist = []
    for i in range(n_coordinates):
        l_dist += [[]]
    
    for i in range(n_models):

        for j in range(i+1,n_models):

            if i == j:
                continue

            for k in range(n_coordinates):
                xi = d_coordinates[i+1][k][0]
                yi = d_coordinates[i+1][k][1]
                zi = d_coordinates[i+1][k][2]
                xj = d_coordinates[j+1][k][0]
                yj = d_coordinates[j+1][k][1]
                zj = d_coordinates[j+1][k][2]
                dist = math.sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
                l_dist[k] += [dist]

    l_rmsfs = []
    for i in range(len(l_dist)):
        sum_x = sum(l_dist[i])
        n = len(l_dist[i])
        average = sum_x/n
        sum_xx = 0
        for j in range(n):
            sum_xx += l_dist[i][j]**2
        SS = sum_xx-sum_x**2/n
        stddev = math.sqrt(SS/n)
##        print average, stddev/average
        l_rmsfs += [average]
    print l_rmsfs

    return l_rmsfs


def parse_coordinates(pdb,chain,):

    path_pdb = '/data/remediated_pdb/'

    fd = open('%s%s/pdb%s.ent' %(path_pdb,pdb[1:3],pdb),'r')
    lines = fd.readlines()
    fd.close()

    d_coordinates = {}

    for line in lines:

        record  = line[:6].strip()

        if record == 'ATOM':

            atom_name = line[12:16].strip()

            if atom_name == 'CA' and line[21] == chain:

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                d_coordinates[model] += [[x,y,z]]

        elif record == 'MODEL':

            model = int(line[10:14])

            d_coordinates[model] = []

    return d_coordinates


if __name__ == '__main__':
    main()

#!/software/bin/python2.3
#
#$Id: pca.py 253 2007-09-07 16:17:28Z tc $
#
#Tommy Carstensen, University College Dublin, 2007

import os

def main():

    path = '/oxygenase_local/tc/PCA_PKA/'

    pdbs = os.listdir('%saverage' %(path))

    ## calculate eigenvectors

    for pdb in pdbs:
        f_eigenvectors(path, pdb)

    file_hessian = '1ATP_eigenvectors.txt'
    fd = open(path+file_hessian)
    raw = fd.read().split()
    fd.close()

    eigenvectors_hessian = []
    for val in raw[2:]:
        if val == '[':
            eigenvectors_hessian += [[]]
        elif val[-1] != ']':
            eigenvectors_hessian[-1] += [float(val)]
        else:
            eigenvectors_hessian[-1] += [float(val[:-1])]

    files = os.listdir(path)
    files.remove(file_hessian)
    files.remove('table.html')
    d_overlap = {}
    for file in files:
        print file
        matrix = parse_g_covar_matrix(path, file)
        eigenvalues_covariance, eigenvectors_covariance = diagonalization(matrix)
        d_overlap[file] = {}
        for mode_hessian in range(6,10):
            i = mode_hessian-6
            eigenvector_hessian = eigenvectors_hessian[i]
            d_overlap[file][mode_hessian] = {}
            for mode_covariance in range(0+1,4+1):
                j = -mode_covariance
                eigenvector_covariance = eigenvectors_covariance[j]
                overlap = cosangle(eigenvector_hessian, eigenvector_covariance)
                d_overlap[file][mode_hessian][mode_covariance] = overlap
                print mode_hessian, mode_covariance, overlap

    writehtml(d_overlap)


def f_eigenvectors(path, pdb):

    import sys
    sys.path.append('/home/people/tc/python/goodvibes/')
    import goodvibes
    instance_goodvibes = goodvibes.vibration()

    fd = open('%saverage/%s' %(path,pdb))
    pdb1lines = fd.readlines()
    fd.close()
    chains1 = ['A']
    pdb1model = [1]
    atoms_hessian = ['CA']
    cutoff_distance = 10.
    verbose = True
    jobid = pdb[:-4]

    pdb1ATOM_all, HEADERdepDate1, REMARK350transformations1, REMARK2resolution1, COMPNDchains1, helices1, strands1 = instance_goodvibes.parse_pdb(pdb1lines, chains1, pdb1model)
    coordinates1 = []
    vectors_difference = []
    for chain1 in chains1:
        for residue in pdb1ATOM_all[chain1]['residues'].keys():
            for atom in pdb1ATOM_all[chain1]['residues'][residue]['atoms'].keys():
                if atom in atoms_hessian:
                    coordinates1.append(pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['coordinates'])
                    pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['hessian'] = 1
    matrix_hessian = instance_goodvibes.hessian_calculation(coordinates1, float(cutoff_distance), verbose) ## calculate with coordinates2 as well... and compare results to switching pdb1 and pdb2
    eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_combined_nonperturbed = instance_goodvibes.eigenv_calccomb(matrix_hessian, jobid, verbose)

    stop

    return


def writehtml(d_overlap):

    html = '<html>\n<head>\n<title>hessian vs covariance</title>\n</head>\n'
    html += '<body>\n\n'
    html += 'rows are hessian modes 7-10\n'
    html += 'columns are covariance modes 1-4\n'
    html += '<table border="1">\n'
    html += '<tr>\n<td>covariance1</td>\n<td>covariance2</td>\n<td>covariance3</td>\n<td>covariance4</td>\n</tr>\n'

    for file in d_overlap.keys():
        html += '<tr>\n<td rowspan="4">%s</td>\n' %(file)
        modes_hessian = d_overlap[file].keys()
        modes_hessian.sort()
        for mode_hessian in modes_hessian:
            if mode_hessian != modes_hessian[0]:
                html += '<tr>\n'
            modes_covariance = d_overlap[file][mode_hessian].keys()
            modes_covariance.sort()
            for mode_covariance in modes_covariance:
                html += '<td>%s</td>\n' %(d_overlap[file][mode_hessian][mode_covariance])
            html += '</tr>\n'
    html += '</table>\n\n</body>\n</html>'

    fd = open('table.html', 'w')
    fd.write(html)
    fd.close()

    return


def cosangle(v1, v2):

    import math

    numerator = 0
    denominator1 = 0
    denominator2 = 0
    for i in range(len(v1)):
        numerator += v1[i]*v2[i]
        denominator1 += v1[i]*v1[i]
        denominator2 += v2[i]*v2[i]
    denominator = math.sqrt(denominator1*denominator2)
    cosang = numerator / denominator

    return cosang


def parse_g_covar_matrix(path, file):

    import math, Numeric

    fd = open(path+file,'r')
    cells_matrix = fd.read().split()
    fd.close()
    N = int(math.sqrt(len(cells_matrix)))
    matrix_covariance = Numeric.zeros((N,N), typecode='d')

    for i in range(N**2):
        val = float(cells_matrix[i])
        row = i/N
        col = i%N
        matrix_covariance[row][col] = val

    return matrix_covariance


def diagonalization(matrix):

    import LinearAlgebra

    eigen_tuple = LinearAlgebra.Heigenvectors(matrix)
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

    return eigenvalues, eigenvectors
    

if __name__ == '__main__':
    main()

##covar_ATP-ALL.dat
##6 -1 0.268914476717
##6 -2 0.437093787906
##6 -3 0.64850354151
##6 -4 -0.105146692206
##7 -1 -0.27108123851
##7 -2 -0.580048057223
##7 -3 0.476473948695
##7 -4 0.0780052379431
##8 -1 -0.25849576569
##8 -2 0.0439245298214
##8 -3 0.271381685615
##8 -4 0.212737770757
##9 -1 0.221297385631
##9 -2 -0.158709209643
##9 -3 -0.0601499336994
##9 -4 0.263257784446
##covar_ATP-ALL.dat (-ref)
##6 -1 -0.565240546224
##6 -2 0.0723194829344
##6 -3 -0.315166459197
##6 -4 -0.348905540667
##7 -1 0.561657140027
##7 -2 -0.106421704624
##7 -3 -0.514629747959
##7 -4 -0.294746384239
##8 -1 0.0273200017243
##8 -2 0.265735239381
##8 -3 -0.26635523018
##8 -4 -0.189986044003
##9 -1 -0.0703243142751
##9 -2 -0.253495728543
##9 -3 -0.03178923461
##9 -4 -0.0520661810065

##covar_ATP_ALL.xpm.dat
##6 -1 -0.240549949569
##6 -2 0.0804732244635
##6 -3 -0.482002328939
##6 -4 -0.0618177921432
##7 -1 -0.0405020407658
##7 -2 0.126652020438
##7 -3 0.126505530167
##7 -4 -0.33237130928
##8 -1 0.0302035311296
##8 -2 0.00322205719543
##8 -3 0.247868515536
##8 -4 -0.147874775463
##9 -1 -0.0679369220745
##9 -2 0.240698138282
##9 -3 0.0678863465103
##9 -4 -0.18725762494
##covar_ATP_ALL.dat (-ref)
##6 -1 0.0047179883912
##6 -2 0.288393918631
##6 -3 0.0530337759244
##6 -4 -0.514915490669
##7 -1 -0.0106361154784
##7 -2 0.0800844726588
##7 -3 0.0808868940439
##7 -4 0.103229300356
##8 -1 -0.0599513470977
##8 -2 0.0165230964581
##8 -3 -0.0458569266272
##8 -4 0.223872117614
##9 -1 -0.133323725718
##9 -2 0.212421634029
##9 -3 0.146180992115
##9 -4 0.0321502275094


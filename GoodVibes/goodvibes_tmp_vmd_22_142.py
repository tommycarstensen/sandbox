##!/bin/env /software/bin/python2.3
##
##$Id: goodvibes2.py 118 2007-10-04 10:22:23Z tc $
##
##Tommy Carstensen, University College Dublin, 2007

## this version does not include helix and proline rigidity

class vibration:


    def main(
        self, jobid, lines, atoms_hessian = ['CA'], frames = 50,
        cutoff_distance = 10.,
        path_python = None, verbose = False, paralleldir = '',
        biomolecule = None, chains = [], model = None,
        ):

        '''
        Use first model if no model specified by user.
        chain(s): Y, biomolecule: Y; parse chains specified by user and apply transformation
        chain(s): Y, biomolecule: N; parse chains specified by user but don't apply transformation
        chain(s): N, biomolecule: Y; parse chains of biomolecule and apply transformation
        chain(s): N, biomolecule: N; parse chains of first biomolecule and apply transformation
        '''

        import os, Numeric, goodvibes_core

        results = []

        ## parse pdb
        (
            d_REMARK350,
            d_primary, ## i.e. SEQRES, MODRES
            d_secondary, ## i.e. HELIX, SHEET
            d_coordinates, ## i.e. ATOM, HETATM, TER, MODEL, ENDMDL
            d_ligands,
            ) = goodvibes_core.parse_pdb(lines, chains)

        ## assume multimeric biological unit if chains not specified by user
        if chains == []:
            chains = d_coordinates['chains'].keys()
            chains.sort()

        ##
        ## calculate N and convert coordinates from dic to list
        ##
        N, d_hessian, l_coordinates = goodvibes_core.parse_dictionary_of_coordinates(d_coordinates, chains, atoms_hessian)

        ##
        ## calculate distance matrix
        ##
        matrix_distances = goodvibes_core.calculate_distance_matrix(l_coordinates)

        ##
        ## calculate hessian matrix
        ##
        matrix_hessian = self.hessian_calculation(N, d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, matrix_distances, l_coordinates, verbose = verbose)

        ##
        ## diagonalize hessian matrix
        ##
        eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_comb_nonperturbed, = goodvibes_core.eigenv_calccomb(
                matrix_hessian, jobid, verbose,
                )

        ##
        ## visualize eigenvectors
        ##
        goodvibes_core.morph(
            eigenvectors_nonperturbed, frames, chains, d_coordinates,
            jobid, d_primary,
            )

        ##
        ## do plots prior to perturbation
        ##
##        goodvibes_core.pre_perturbation_plot(
##            jobid,cutoff_distance,chains,d_secondary,
##            eigenvectors_nonperturbed,eigenvectors_comb_nonperturbed,
##            )

        ##
        ## set data lists and append matrices to be plotted for each combination of modes 6-12 before initiating loops over the two axes of the plot
        ##
        datadic = goodvibes_core.datadic_return(N)

        prefix = 'wt'
        suffix = 'vmd'
        goodvibes_core.overlap2structure(
            eigenvectors_nonperturbed, eigenvectors_nonperturbed,
            d_coordinates,[],d_hessian,l_coordinates,prefix,suffix,
            )
        stop

################################################################################
##                        plotdata calculation section                        ##
################################################################################

        ## loop over remres1
        for remres1 in range(N):

            print remres1

            ## continue the remres1 loop if another processor is running calculations for that value of remres1
            filename = '%s/%s_residue%s.txt' %(paralleldir, jobid, remres1+1)
            datadic, Continue = goodvibes_core.multiple_cpus_read(datadic,filename,remres1)
            if Continue == True and remres1 not in [20,21,141,19,22,140,142,]:
                continue

            ## loop over remres2
            for remres2 in range(remres1+1,N):

                if remres2 not in [20,21,141,19,22,140,142,]:
                    continue

                print remres1, remres2

                l_rem = [remres1,remres2]

                l_coordinates_perturbed = l_coordinates[:remres1]+l_coordinates[remres1+1:remres2]+l_coordinates[remres2+1:]
            
                matrix_hessian_perturbed = self.hessian_calculation(
                    N-2, d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, matrix_distances, l_coordinates_perturbed, verbose = verbose,
                    )

                (
                    eigenvectors_perturbed, eigenvalues_perturbed, eigenvectors_perturbed_combined,
                    ) = goodvibes_core.eigenv_calccomb(
                        matrix_hessian_perturbed, jobid, verbose
                        )

                eigenvectors_nonperturbed_aligned, eigenvectors_perturbed_aligned = goodvibes_core.align_vectors_of_different_length(eigenvectors_nonperturbed, eigenvectors_perturbed, l_rem = l_rem)
                eigenvectors_nonperturbed_aligned, eigenvectors_perturbed_combined_aligned = goodvibes_core.align_vectors_of_different_length(eigenvectors_nonperturbed, eigenvectors_perturbed_combined, l_rem = l_rem)

                (
                    overlaps_single, max_overlaps_single, perturbed_modes_of_max_overlap_single, delta_perturbed_eigenvalues_of_max_overlap,
                    ) = goodvibes_core.overlap_calculation(
                        eigenvectors_perturbed_aligned, eigenvectors_nonperturbed_aligned, eigenvalues_perturbed, eigenvalues_nonperturbed,
                        )

                (
                    overlaps_combined, max_overlaps_combined, perturbed_modes_of_max_overlap_combined,
                    ) = goodvibes_core.overlap_calculation(
                        eigenvectors_perturbed_combined_aligned, eigenvectors_nonperturbed_aligned, eigenvalues_perturbed, eigenvalues_nonperturbed,
                        )[:-1]

                datadic_loop = {
                    'eigenvalues_perturbed': eigenvalues_perturbed,
                    'emo': delta_perturbed_eigenvalues_of_max_overlap,
                    'overlaps_single': overlaps_single,
                    'overlaps_max': max_overlaps_single,
                    'mmo': perturbed_modes_of_max_overlap_single,
                    'overlaps_combined': overlaps_combined
                    }
                for key in datadic:
                    for mode in range(6,12):
                        datadic[key]['data'][mode][remres1][remres2] = datadic_loop[key][mode]
                        datadic[key]['data'][mode][remres2][remres1] = datadic_loop[key][mode]
                datadic['overlaps_combined']['data'][-1][remres1][remres2] = overlaps_combined[-1]
                datadic['overlaps_combined']['data'][-1][remres2][remres1] = overlaps_combined[-1]

                print overlaps_single[6]

                prefix = '%s' %(jobid)
                suffix = '%s_%s' %(remres1+1, remres2+1)
                goodvibes_core.overlap2structure(
                    eigenvectors_nonperturbed, eigenvectors_perturbed,
                    d_coordinates,l_rem,d_hessian,l_coordinates,prefix,suffix,
                    )

                prefix = '%s' %(jobid)
                suffix = '%s_%s' %(remres1+1, remres2+1)
                goodvibes_core.vmdarrow(
                    eigenvectors_perturbed,d_hessian,
                    prefix,suffix,l_rem=l_rem,
                    )
                if overlaps_single[6] > 0.9:
                    print overlaps_single[6]
                    print remres1+1, remres2+1
##                    stop

##            ##
##            ## write data to txt files in case of crash during loop over residues/clusters
##            ##
##            filename = '%s/%s_residue%s.txt' %(paralleldir, jobid, remres1+1)
##            goodvibes_core.multiple_cpus_write(datadic,filename,remres1)

##        ##
##        ## do plots for perturbation results
##        ##
##        if verbose == True:
##            print 'generating plots'
##        goodvibes_core.post_perturbation_plot(
##            datadic, jobid, chains, cutoff_distance, d_hessian, N, d_secondary, d_coordinates,
##            )

        return


    def hessian_calculation(self, N, d_coordinates, chains, atoms_hessian, cutoff, d_secondary, matrix_distances, l_coordinates, l_rem = [], verbose = False):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

        if verbose == True:
            print 'calculating the hessian/second-order derivate Kirchoff/Laplacian matrix'
        
        import math, Numeric

        cutoff_sq = cutoff**2

        matrix_hessian = Numeric.zeros((3*N,3*N), typecode='d')
        
        for row_sup in range(N):
            for col_sup in range(N):

                if col_sup > row_sup:
                    #does the Numeric module feature some smart built-in function to calculate length of vectors? use math.sqrt(math.pow(vector, 2))
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
##                    if dist_sq <= cutoff_sq:
                    sigmoidfactor = self.sigmoid(math.sqrt(dist_sq), cutoff)
                    vector = [x,y,z]
                    for row_sub in range(3):
                        for col_sub in range(3):

                            if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                                value = sigmoidfactor*-vector[row_sub]*vector[col_sub]/dist_sq
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


    def sigmoid(self, x, cutoff, slope=1):
        import math
        y = 1. / ( 1. + math.exp( slope*(x-cutoff) ) )
        return y

    def __init__(self):
        self.residue_terminal_n = {'A':-9999, 'B':-9999, 'C':-9999, 'D':-9999, 'E':-9999, 'F':-9999, 'G':-9999, 'H':-9999, 'I':-9999, 'J':-9999, 'K':-9999, 'L':-9999, 'M':-9999, 'N':-9999, 'O':-9999, 'P':-9999, 'Q':-9999, 'R':-9999, 'S':-9999, 'T':-9999, 'U':-9999, 'V':-9999, 'W':-9999, 'X':-9999, 'Y':-9999, 'Z':-9999, ' ':-9999,}
        self.residue_terminal_c = {'A':9999, 'B':9999, 'C':9999, 'D':9999, 'E':9999, 'F':9999, 'G':9999, 'H':9999, 'I':9999, 'J':9999, 'K':9999, 'L':9999, 'M':9999, 'N':9999, 'O':9999, 'P':9999, 'Q':9999, 'R':9999, 'S':9999, 'T':9999, 'U':9999, 'V':9999, 'W':9999, 'X':9999, 'Y':9999, 'Z':9999, ' ':9999,}
        ##
        self.chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'] ## used for remark350build
        self.d_res = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            'PCA':'X','ACE':'X','SEP':'X','TPO':'X'
            }
        self.d_res20 = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y'}
        self.weekdaydic = {0:'Monday',1:'Tuesday',2:'Wednesday',3:'Thursday',4:'Friday',5:'Saturday',6:'Sunday'}
        self.monthdic = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
        self.path_pdb = '/oxygenase_local/data/pdb/'

if __name__=='__main__':

    import os, sys, goodvibes_core

    if '-pdb' not in sys.argv:
        raise 'use -pdb to select a pdb'
    pdb = sys.argv[sys.argv.index('-pdb')+1][:4]
    if '-chains' in sys.argv:
        chains = sys.argv[sys.argv.index('-chains')+1].split(',')
    else:
        chains = []
    if '-biomolecule' in sys.argv:
        biomolecule = int(sys.argv[sys.argv.index('-biomolecule')+1])
    else:
        biomolecule = None

    instance_vibration = vibration()
    lines = goodvibes_core.pdb_import(pdb, instance_vibration.path_pdb)

    instance_vibration.main(
        pdb,
        lines,
        chains = chains,
        biomolecule = biomolecule,
        verbose = True,
        paralleldir = os.getcwd(),
        )

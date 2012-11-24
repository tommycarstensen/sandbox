##!/bin/env /software/bin/python2.3
##
##$Id$
##
##Tommy Carstensen, University College Dublin, 2007

''' size of residues are increased by dividing the cutoff distance by a factor 16'''


class vibration:


    def main(
        self, jobid, lines, atoms_hessian = ['CA'], frames = 50,
        cutoff_distance = 10.,
        path_html = None, path_python = None, verbose = False, paralleldir = '',
        chains = None,
        ):

        import goodvibes_core

        ##
        ## parse pdb
        ##
        (
            d_REMARK350,
            d_primary, ## i.e. SEQRES, MODRES
            d_secondary, ## i.e. HELIX, SHEET
            d_coordinates, ## i.e. ATOM, HETATM, TER, MODEL, ENDMDL
            d_ligands,
            ) = goodvibes_core.parse_pdb(lines, chains)


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
        eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_combined_nonperturbed = goodvibes_core.eigenv_calccomb(matrix_hessian, jobid, verbose)

        ##
        ## visualize eigenvectors
        ##
        goodvibes_core.morph(eigenvectors_nonperturbed, frames, chains, d_coordinates, jobid, d_primary)

        ##
        ## do plots prior to perturbation
        ##
        goodvibes_core.pre_perturbation_plot(
            jobid,cutoff_distance,chains,d_secondary,
            eigenvectors_nonperturbed,eigenvectors_combined_nonperturbed,
            )
        
        ##
        ## set data lists and append matrices to be plotted for each combination of modes 6-12 before initiating loops over the two axes of the plot
        ##
        datadic = goodvibes_core.datadic_return(N)

        ##
        ## loop over remres1
        ##
        for remres1 in range(N-1):

            ## continue the remres1 loop if another processor is running calculations for that value of remres1
            filename = '%s/%s_residue%s.txt' %(paralleldir, jobid, remres1+1)
            datadic, Continue = goodvibes_core.multiple_cpus_read(datadic,filename,remres1)
            if Continue == True:
                continue

            ##
            ## loop over remres2
            ##
            for remres2 in range(remres1+1,N):

                print remres1, remres2
                l_mod = [remres1,remres2]

                ##
                ## calculate hessian matrix
                ##
                matrix_hessian_perturbed = self.hessian_calculation(
                    N, d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, matrix_distances, l_coordinates, verbose = verbose, l_mod = l_mod,
                    )

                ##
                ## diagonalize hessian matrix and calculate eigenvectors
                ##
                (
                    eigenvectors_perturbed, eigenvalues_perturbed, eigenvectors_perturbed_combined
                    ) = goodvibes_core.eigenv_calccomb(
                        matrix_hessian_perturbed, jobid, verbose
                        )

                ##
                ## calculate overlap between eigenvectors
                ##
                (
                    overlaps_single, max_overlaps_single, perturbed_modes_of_max_overlap_single, delta_perturbed_eigenvalues_of_max_overlap
                    ) = goodvibes_core.overlap_calculation(
                        eigenvectors_perturbed, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed,
                        )
                
                (
                    overlaps_combined, max_overlaps_combined, perturbed_modes_of_max_overlap_combined
                    ) = goodvibes_core.overlap_calculation(
                        eigenvectors_perturbed_combined, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed,
                        )[:-1]
                print overlaps_single[6]

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

            ## write data to txt files in case of crash during loop over residues/clusters
            filename = '%s/%s_residue%s.txt' %(paralleldir, jobid, remres1+1)
            goodvibes_core.multiple_cpus_write(datadic,filename,remres1)

        ##
        ## do plots for perturbation results
        ##
        if verbose == True:
            print 'generating plots'
        goodvibes_core.post_perturbation_plot(
            datadic, jobid, chains, cutoff_distance, d_hessian, N, d_secondary, d_coordinates,
            )


        return


    def hessian_calculation(self, N, d_coordinates, chains, atoms_hessian, cutoff, d_secondary, matrix_distances, l_coordinates, l_mod = [], verbose = False):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

        if verbose == True:
            print 'calculating the hessian/second-order derivate Kirchoff/Laplacian matrix'
        
        import math, Numeric, time

        matrix_hessian = Numeric.zeros( (3*N,3*N), typecode='d')

        ##
        ## write matrix
        ##
        row_sup = 0
        col_sup = 0
        for row_sup in range(N):
            for col_sup in range(row_sup+1,N):
                sqdist = matrix_distances[row_sup][col_sup]

                if row_sup in l_mod and col_sup in l_mod:
##                    sqdist = 1. ## large residue!!!
                    sqdist /= 16. ## large residue!!!
                elif row_sup in l_mod:
##                    sqdist = 3. ## large residue!!!
                    sqdist /= 8. ## large residue!!!
                elif col_sup in l_mod:
##                    sqdist = 3. ## large residue!!!
                    sqdist /= 8. ## large residue!!!

                vector = l_coordinates[row_sup]-l_coordinates[col_sup]

                matrix_hessian = self.hessian_sub(
                    matrix_hessian,row_sup,col_sup,cutoff,d_coordinates,d_secondary,sqdist,vector,
                    N,
                    )

        return matrix_hessian


    def hessian_sub(
        self,matrix_hessian,row_sup,col_sup,dist_cutoff,d_coordinates,d_secondary,sqdist,vector,
        N,
        ):

        import math, goodvibes_core

        sqdist_cutoff = dist_cutoff**2
##        if sqdist_cutoff < sqdist: ## temp!!! avoid rigidity etc...
##            return matrix_hessian

        factor = 1

        matrix_hessian = goodvibes_core.fill_matrix_hessian(matrix_hessian, factor, row_sup, col_sup, sqdist, vector)

        return matrix_hessian


    def __init__(self):
        self.residue_terminal_n = {'A':-9999, 'B':-9999, 'C':-9999, 'D':-9999, 'E':-9999, 'F':-9999, 'G':-9999, 'H':-9999, 'I':-9999, 'J':-9999, 'K':-9999, 'L':-9999, 'M':-9999, 'N':-9999, 'O':-9999, 'P':-9999, 'Q':-9999, 'R':-9999, 'S':-9999, 'T':-9999, 'U':-9999, 'V':-9999, 'W':-9999, 'X':-9999, 'Y':-9999, 'Z':-9999, ' ':-9999,}
        self.residue_terminal_c = {'A':9999, 'B':9999, 'C':9999, 'D':9999, 'E':9999, 'F':9999, 'G':9999, 'H':9999, 'I':9999, 'J':9999, 'K':9999, 'L':9999, 'M':9999, 'N':9999, 'O':9999, 'P':9999, 'Q':9999, 'R':9999, 'S':9999, 'T':9999, 'U':9999, 'V':9999, 'W':9999, 'X':9999, 'Y':9999, 'Z':9999, ' ':9999,}
        ##
        self.chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'] ## used for remark350build
        self.d_res = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            }
        self.weekdaydic = {0:'Monday',1:'Tuesday',2:'Wednesday',3:'Thursday',4:'Friday',5:'Saturday',6:'Sunday'}
        self.monthdic = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
        self.path_pdb = '/oxygenase_local/data/pdb/'

if __name__=='__main__':

    import sys, os, goodvibes_core

    if '-pdb' not in sys.argv:
        raise 'use -pdb to select a pdb'
    pdb = sys.argv[sys.argv.index('-pdb')+1][:4]
    if '-chains' in sys.argv:
        chains = sys.argv[sys.argv.index('-chains')+1].split(',')
    else:
        raise 'use -chains to select chain(s)'
    if '-atoms' in sys.argv:
        atoms = sys.argv[sys.argv.index('-atoms')+1].split(',')
    else:
        atoms = ['CA']
    if '-cutoff' in sys.argv:
        cutoff_distance = float(sys.argv[sys.argv.index('-cutoff')+1])
    else:
        cutoff_distance = 6.
    for arg in sys.argv:
        if '-' in arg and arg not in ['-atoms','-cutoff','-chains','-pdb']:
            raise '%s is an unknown flag' %(arg)

    instance_vibration = vibration()
    lines = goodvibes_core.pdb_import(pdb, instance_vibration.path_pdb)

    instance_vibration.main(
        pdb,
        lines,
        chains = chains,
        atoms_hessian = atoms,
        cutoff_distance = cutoff_distance,
        verbose = True,
        paralleldir = os.getcwd(),
        )

##!/bin/env /software/bin/python2.3
##
##$Id$
##
##Tommy Carstensen, University College Dublin, 2007


class vibration:


    def main(
        self, jobid, lines, atoms_hessian = ['CA'], frames = 50,
        cutoff_distance = 10.,
        path_html = None, path_python = None, verbose = False, paralleldir = '',
        chains = None,
        winsize = 1,
        pre_perturbation_plot = True,
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
        none, d_coordinates, none, none = self.pre_hessian(d_coordinates, chains, atoms_hessian)

        ##
        ## calculate distance matrix
        ##
        matrix_distances = goodvibes_core.calculate_distance_matrix(l_coordinates)

        ##
        ## calculate hessian matrix
        ##
        matrix_hessian = self.hessian_calculation(
            N, d_coordinates, chains, atoms_hessian, cutoff_distance, d_secondary, l_rem = [], verbose = verbose,
            )

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
        if pre_perturbation_plot == True:
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
        for remres1 in range((winsize-1)/2,N-winsize-(winsize-1)/2):

            ##
            ## continue the remres1 loop if another processor is running calculations for that value of remres1
            ##
            filename = '%s/%s_residue%s.txt' %(paralleldir, jobid, remres1+1)
            datadic, Continue = goodvibes_core.multiple_cpus_read(datadic,filename,remres1)
            if Continue == True:
                continue

            ##
            ## loop over remres2
            ##
            for remres2 in range(remres1+2*(winsize-1)/2+1,N-(winsize-1)/2):

                print remres1, remres2

                l_rem = []
                for i in range(-(winsize-1)/2,(winsize+1)/2):
                    l_rem += [remres1+i,remres2+i]
                l_rem.sort()

                ##
                ## remove selected alpha carbon atoms!!!
                ##
                l_coordinates_perturbed = l_coordinates[:remres1-(winsize-1)/2]+l_coordinates[remres1+1+(winsize-1)/2:remres2-(winsize-1)/2]+l_coordinates[remres2+1+(winsize-1)/2:]

                matrix_hessian_perturbed = self.hessian_calculation(
                    N-len(l_rem), d_coordinates, chains, atoms_hessian, cutoff_distance, d_secondary, l_rem = l_rem, verbose = verbose,
                    )

                (
                    eigenvectors_perturbed, eigenvalues_perturbed, eigenvectors_perturbed_combined
                    ) = goodvibes_core.eigenv_calccomb(
                        matrix_hessian_perturbed, jobid, verbose
                        )
                
                (
                    overlaps_single, max_overlaps_single, perturbed_modes_of_max_overlap_single, delta_perturbed_eigenvalues_of_max_overlap
                    ) = goodvibes_core.overlap_calculation(
                        eigenvectors_perturbed, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem = l_rem,
                        )
                
                (
                    overlaps_combined, max_overlaps_combined, perturbed_modes_of_max_overlap_combined
                    ) = goodvibes_core.overlap_calculation(
                        eigenvectors_perturbed_combined, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem = l_rem,
                        )[:-1]
                print overlaps_single[6]

        ##        ## do vmd of perturbed structure
        ##        self.morph(eigenvectors, frames, biomolecule, d_coordinates, atoms_hessian, cluster, matrix_hessian, jobid+'-'+str(xvalue)+'-'+str(yvalue), cutoff_distance)

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
        

    def hessian_calculation(self, N, d_coordinates, chains, atoms_hessian, cutoff, d_secondary, l_rem = [], verbose = False):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

        if verbose == True:
            print 'calculating the hessian/second-order derivate Kirchoff/Laplacian matrix'
        
        import math, Numeric, time

##        matrix_hessian = Numeric.zeros( (3*(N-len(l_rem)),3*(N-len(l_rem)) ), typecode='d')
        matrix_hessian = Numeric.zeros( (3*N,3*N), typecode='d')

        ##
        ## write matrix
        ##
        row_sup = 0
        for chain1 in chains:
            ## assume sequential numbering of residues
            res_nos1 = d_coordinates['chains'][chain1]['residues'].keys()
            res_nos1.sort()
            for res_no1 in res_nos1:
                if res_no1-1 in l_rem:
                    continue
                ## assume sequential iCodes
                iCodes1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'].keys()
                iCodes1.sort()
                for iCode1 in iCodes1:
                    altloc1_residue = min(d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['altlocs'].keys())
                    res_name1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['altlocs'][altloc1_residue]['res_name']
                    atom_names1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'].keys()
                    atom_names1.sort()
                    for atom_name1 in atom_names1:
                        if atom_name1 in atoms_hessian:
                            altloc1 = min(d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'].keys())
                            coordinate1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['coordinate']

                            col_sup = 0
                            for chain2 in chains:
                                ## assume sequential numbering of residues
                                res_nos2 = d_coordinates['chains'][chain2]['residues'].keys()
                                res_nos2.sort()
                                for res_no2 in res_nos2:
                                    if res_no2-1 in l_rem:
                                        continue
                                    ## assume sequentail iCodes
                                    iCodes2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'].keys()
                                    iCodes2.sort()
                                    for iCode2 in iCodes2:
                                        altloc2_residue = min(d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['altlocs'].keys())
                                        res_name2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['altlocs'][altloc2_residue]['res_name']
                                        atom_names2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'].keys()
                                        atom_names2.sort()
                                        for atom_name2 in atom_names2:
                                            if atom_name2 in atoms_hessian:
                                                altloc2 = min(d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'].keys())
                                                coordinate2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'][altloc2]['coordinate']

                                                if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                                    sqdist_main = 0
                                                    sqdist_side = 0
                                                    sqdist_ala = 0
                                                    sqdist_ala1 = 0
                                                    sqdist_ala2 = 0
                                                else:
                                                    d_sqdist1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['sqdist']['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]
                                                    sqdist_main1 = d_sqdist1['sqdist_main']
                                                    sqdist_side1 = d_sqdist1['sqdist_side']
                                                    sqdist_ala1 = d_sqdist1['sqdist_ala']
                                                    sqdist_ala11 = d_sqdist1['sqdist_ala1']
                                                    sqdist_ala21 = d_sqdist1['sqdist_ala2']
                                                    d_sqdist2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['sqdist']['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]
                                                    sqdist_main2 = d_sqdist2['sqdist_main']
                                                    sqdist_side2 = d_sqdist2['sqdist_side']
                                                    sqdist_ala2 = d_sqdist2['sqdist_ala']
                                                    sqdist_ala12 = d_sqdist2['sqdist_ala1']
                                                    sqdist_ala22 = d_sqdist2['sqdist_ala2']

                                                    if sqdist_main1 != sqdist_main2 or sqdist_side1 != sqdist_side2 or sqdist_ala1 != sqdist_ala2 or sqdist_ala11 != sqdist_ala22 or sqdist_ala21 != sqdist_ala12:
                                                        print chain1, res_no1, iCode1, res_name1, altloc1, atom_name1
                                                        print chain2, res_no2, iCode2, res_name2, altloc2, atom_name2
                                                        print sqdist_main1, sqdist_side1, sqdist_ala1, sqdist_ala11, sqdist_ala21
                                                        print sqdist_main2, sqdist_side2, sqdist_ala2, sqdist_ala12, sqdist_ala22
                                                        import math
                                                        print math.sqrt(sqdist_main1), math.sqrt(sqdist_side1)
                                                        print math.sqrt(sqdist_main2), math.sqrt(sqdist_side2)
                                                        notexpected
                                                    else:
                                                        sqdist_main = sqdist_main1
                                                        sqdist_side = sqdist_side1
                                                        sqdist_ala = sqdist_ala1
                                                        sqdist_ala1 = sqdist_ala11
                                                        sqdist_ala2 = sqdist_ala21

                                                    if sqdist_ala1 < min(sqdist_main,sqdist_side):
                                                        stop1
                                                    if sqdist_ala2 < min(sqdist_main,sqdist_side):
                                                        stop2
                                                    if sqdist_ala1 > sqdist_ala:
                                                        import math
                                                        print math.sqrt(sqdist_ala1), math.sqrt(sqdist_ala)
                                                        print chain1, chain2, res_no1, res_no2, atom_name1, atom_name2
                                                        stop3
                                                    if sqdist_ala2 > sqdist_ala:
                                                        import math
                                                        print math.sqrt(sqdist_ala2), math.sqrt(sqdist_ala)
                                                        stop4

                                                if l_rem == []:
                                                    sqdist = min(sqdist_main,sqdist_side)
                                                elif row_sup in l_rem and col_sup in l_rem:
                                                    sqdist = sqdist_ala
                                                elif row_sup in l_rem:
                                                    sqdist = sqdist_ala1
                                                elif col_sup in l_rem:
                                                    sqdist = sqdist_ala2
                                                else:
                                                    sqdist = min(sqdist_main,sqdist_side)

                                                if row_sup >= col_sup:
                                                    pass
                                                elif (
                                                    chain1 == chain2 and
                                                    res_no1 == res_no2 and
                                                    iCode1 == iCode2 and
                                                    atom_name1 == atom_name2
                                                    ):
                                                    print chain1, res_no1, iCode1, res_name1, atom_name1
                                                    print row_sup, col_sup
                                                    notexpected
                                                elif chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2:
                                                    stop
                                                    if atom_name1 == 'CA' and atom_name2 == 'CB':
                                                        matrix_hessian = self.hessian_sub(
                                                            matrix_hessian,row_sup,col_sup,cutoff,d_coordinates,d_secondary,sqdist,coordinate1,coordinate2,
                                                            chain1,res_no1,iCode1,res_name1,atom_name1,altloc1,
                                                            chain2,res_no2,iCode2,res_name2,atom_name2,altloc2,
                                                            )
                                                    elif atom_name1 == 'CB' and atom_name2 == 'CA' and chain1 == chain2:
                                                        print chain1, chain2, res_no1, res_no2, iCode1, iCode2
                                                        notexpected1
                                                    elif atom_name1 == 'CA' and atom_name2 == 'CA' and chain1 == chain2:
                                                        print chain1, chain2, res_no1, res_no2, iCode1, iCode2
                                                        notexpected2
                                                    else:
                                                        print atom_name1, atom_name2
                                                        notexpected3
                                                else:
                                                    if atom_name1 in ['CA','CB'] and atom_name2 in ['CA','CB']:
                                                        matrix_hessian = self.hessian_sub(
                                                            matrix_hessian,row_sup,col_sup,cutoff,d_coordinates,d_secondary,sqdist,coordinate1,coordinate2,
                                                            chain1,res_no1,iCode1,res_name1,atom_name1,altloc1,
                                                            chain2,res_no2,iCode2,res_name2,atom_name2,altloc2,
                                                            )
                                                    else:
                                                        print atom_name1, atom_name2
                                                        notexpected5

                                                col_sup += 1

                            row_sup += 1

        return matrix_hessian


    def pre_hessian(self, d_coordinates, chains, atoms_hessian):

        d_hessian = {}
        l_coordinates = []
        l_CA = []
        i = 0

        for chain1 in chains:
            ## assume sequential numbering of residues
            res_nos1 = d_coordinates['chains'][chain1]['residues'].keys()
            res_nos1.sort()
            for res_no1 in res_nos1:
                ## assume sequential iCodes
                iCodes1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'].keys()
                iCodes1.sort()
                for iCode1 in iCodes1:
                    altloc_residue1 = min(d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['altlocs'].keys())
                    res_name1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['altlocs'][altloc_residue1]['res_name']
                    atom_names1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'].keys()
                    atom_names1.sort()

                    ##
                    ## append to l_coordinates
                    ##
                    for atom_name1 in atom_names1:
                        altloc1 = min(d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'].keys())
                        coordinate1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['coordinate']
                        if atom_name1 == 'CA':
                            l_CA += [i]
                        if atom_name1 in atoms_hessian:

                            d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['i'] = i
                            d_hessian[i] = {'chain':chain1,'res_no':res_no1,'iCode':iCode1,'atom_name':atom_name1,'altloc':altloc1, 'res_name': res_name1}
                            l_coordinates += [coordinate1]

                            i += 1

                    for chain2 in chains:
                        ## assume sequential numbering of residues
                        res_nos2 = d_coordinates['chains'][chain2]['residues'].keys()
                        res_nos2.sort()
                        for res_no2 in res_nos2:
                            ## assume sequential iCodes
                            iCodes2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'].keys()
                            iCodes2.sort()
                            for iCode2 in iCodes2:
                                altloc_residue2 = min(d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['altlocs'].keys())
                                res_name2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['altlocs'][altloc_residue2]['res_name']
                                atom_names2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'].keys()
                                atom_names2.sort()

                                if (
                                    chain1 == chain2 and
                                    res_no1 == res_no2 and
                                    iCode1 == iCode2
                                    ):
                                    continue

                                ##
                                ## calculate distances
                                ##
                                l_sqdist_main = []
                                l_sqdist_side = []
                                l_sqdist_ala = []
                                l_sqdist_ala1 = []
                                l_sqdist_ala2 = []
                                for atom_name1 in atom_names1:
                                    ## exclude light atoms
                                    if atom_name1[0] == 'H':
                                        continue
                                    altloc1 = min(d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'].keys())
                                    coordinate1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['coordinate']

                                    for atom_name2 in atom_names2:
                                        ## exclude light atoms
                                        if atom_name2[0] == 'H':
                                            continue
                                        altloc2 = min(d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'].keys())
                                        coordinate2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'][altloc2]['coordinate']

                                        sqdist = sum((coordinate2-coordinate1)**2)
                                        if atom_name1 in ['N','CA','C','O','OXT'] and atom_name2 in ['N','CA','C','O','OXT']:
                                            l_sqdist_main += [sqdist]
                                        else:
                                            l_sqdist_side += [sqdist]
                                        if atom_name1 in ['N','CA','C','O','OXT','CB',]:
                                            l_sqdist_ala1 += [sqdist]
                                        if atom_name2 in ['N','CA','C','O','OXT','CB',]:
                                            l_sqdist_ala2 += [sqdist]
                                        if atom_name1 in ['N','CA','C','O','OXT','CB',] and atom_name2 in ['N','CA','C','O','OXT','CB',]:
                                            l_sqdist_ala += [sqdist]

                                sqdist_main = min(l_sqdist_main)
                                if res_name1 == 'GLY' and res_name2 == 'GLY':
                                    sqdist_side = sqdist_main
                                else:
                                    sqdist_side = min(l_sqdist_side)
                                sqdist_ala = min(l_sqdist_ala)
                                sqdist_ala1 = min(l_sqdist_ala1)
                                sqdist_ala2 = min(l_sqdist_ala2)

                                d_sqdist = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]
                                if 'sqdist' not in d_sqdist.keys():
                                    d_sqdist['sqdist'] = {'chains':{}}
                                if chain2 not in d_sqdist['sqdist']['chains'].keys():
                                    d_sqdist['sqdist']['chains'][chain2] = {'residues':{}}
                                if res_no2 not in d_sqdist['sqdist']['chains'][chain2]['residues'].keys():
                                    d_sqdist['sqdist']['chains'][chain2]['residues'][res_no2] = {'iCodes':{}}
                                if iCode2 not in d_sqdist['sqdist']['chains'][chain2]['residues'][res_no2]['iCodes'].keys():
                                    d_sqdist['sqdist']['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2] = {'res_names':{}}
                                if res_name2 not in d_sqdist['sqdist']['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'].keys():
                                    d_sqdist['sqdist']['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2] = {
                                        'sqdist_main':sqdist_main,'sqdist_side':sqdist_side,'sqdist_ala':sqdist_ala,'sqdist_ala1':sqdist_ala1,'sqdist_ala2':sqdist_ala2,
                                        }

        return l_coordinates, d_coordinates, d_hessian, l_CA


    def hessian_sub(
        self,matrix_hessian,row_sup,col_sup,dist_cutoff,d_coordinates,d_secondary,sqdist,c1,c2,
        chain1,res_no1,iCode1,res_name1,atom_name1,altloc1,
        chain2,res_no2,iCode2,res_name2,atom_name2,altloc2,
        ):

        import math

        sqdist_cutoff = dist_cutoff**2
##        if sqdist_cutoff < sqdist: ## temp!!! avoid rigidity etc...
##            return matrix_hessian

        vector = c2-c1
        dist_sq = sum(vector**2)

        factor = 0
        if (res_name1 == 'PRO' or res_name2 == 'PRO') and (res_no2 == res_no1-1 or res_no2 == res_no1 or res_no2 == res_no1+1) and (iCode1 != ' ' or iCode2 != ' ') :
            expected
        if chain1 == chain2 and (res_no1 == res_no2+4 or res_no2 == res_no1+4) and (iCode1 != ' ' or iCode2 != ' '):
            expected
##        if (atom_name1 == 'CA' or atom_name2 == 'CA') and sqdist_main < sqdist_cutoff:
##            factor = 1
##        if (atom_name1 == 'CB' or atom_name2 == 'CB') and min(sqdist_main,sqdist_side) < sqdist_cutoff:
##            factor = 1
        if (atom_name1 == 'CA' or atom_name2 == 'CA') and sqdist < sqdist_cutoff:
            factor = 1
        if (atom_name1 == 'CB' or atom_name2 == 'CB') and min(sqdist_main,sqdist_side) < sqdist_cutoff:
            factor = 1
            stop

        ## alpha to own beta
        if chain1 == chain2 and res_no1 == res_no2 and iCode1 == iCode2 and ((atom_name1 == 'CA' and atom_name2 == 'CB') or (atom_name2 == 'CA' and atom_name1 == 'CB')):
            if factor == 0:
                stop
            factor *= 2
            stop
        ## proline
        if (
            atom_name1 == 'CA' and atom_name2 == 'CA' and
            (res_name1 == 'PRO' and (res_no2 == res_no1-1 or res_no2 == res_no1+1)) or
            (res_name2 == 'PRO' and (res_no1 == res_no2-1 or res_no1 == res_no2+1))
            ):
            if res_name1 == 'PRO' and res_no2 == res_no1-1:
                notexpected
            if res_name2 == 'PRO' and res_no1 == res_no2+1:
                notexpected
            if factor == 0:
                stop
            factor *= 2
        ## helix
        if atom_name1 == 'CA' and atom_name2 == 'CA' and chain1 == chain2 and (res_no1 == res_no2+4 or res_no2 == res_no1+4) and res_no1 in d_secondary['HELIX'][chain1]['list'] and res_no2 in d_secondary['HELIX'][chain2]['list']:
            if factor == 0:
                factor = 2
            else:
                factor *= 2
        if factor > 2:
            print factor
            print chain1, res_no1, iCode1, res_name1, altloc1, atom_name1
            print chain2, res_no2, iCode2, res_name2, altloc2, atom_name2
            print c2, c1
            print sqdist_main, sqdist_side
            print sum((c2-c1)**2)
            stop

        for row_sub in range(3):
            for col_sub in range(3):

                if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                    value = factor*-vector[row_sub]*vector[col_sub]/dist_sq
                    if (
                        3*row_sup+row_sub == len(matrix_hessian) or
                        3*col_sup+col_sub == len(matrix_hessian[3*row_sup+row_sub])
                        ):
                        print 3*row_sup+row_sub,3*col_sup+col_sub,len(matrix_hessian),len(matrix_hessian[3*row_sup+row_sub])
                        print chain1, chain2, res_no1, res_no2, iCode1, iCode2, atom_name1, atom_name2
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
    if '-winsize' in sys.argv:
        winsize = int(sys.argv[sys.argv.index('-winsize')+1])
    else:
        winsize = 1
    if '-pre_perturbation_plot' in sys.argv:
        pre_perturbation_plot = False
    else:
        pre_perturbation_plot = True

    for arg in sys.argv:
        if '-' in arg and arg not in [
            '-atoms','-cutoff','-chains','-pdb','-winsize','-pre_perturbation_plot',
            ]:
            raise '%s is an unknown flag' %(arg)

    instance_vibration = vibration()
    lines = goodvibes_core.pdb_import(pdb, instance_vibration.path_pdb)

    instance_vibration.main(
        pdb,
        lines,
        chains = chains,
        atoms_hessian = atoms,
        cutoff_distance = cutoff_distance,
        winsize = winsize,
        pre_perturbation_plot = pre_perturbation_plot,
        verbose = True,
        paralleldir = os.getcwd(),
        )

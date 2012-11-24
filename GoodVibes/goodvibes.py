#!/bin/env /software/bin/python2.3
#
# $Id: goodvibes.py 128 2008-07-10 12:55:18Z tc $
#
## make sure that missng atoms and residues don't go ilent by
## hvad skal der ske med atomer med altloc?
## add info E.C. # from pdb.org or pdb (REMARK2?! eller COMPND...)
## return error if a user specified chain is not found in pdb input file
## use info from SEQADV lines if modified residue and MODRES lines

##                resreseigenvalue ## 6 plots for 6 perturbed modes
##                resresmaxoverlap_single ## 1 plot
##                resresmaxoverlap_combined ## 1 plot
##                resreseigenvalueofmaxoverlap ## 1 plot
##                resresmodeofmaxoverlap_single ## 1 plot
##                resresmodeofmaxoverlap_combined ## 1 plot

'''This documentation refers to version 1.55. This module contains only one class. The module is used to do some NMA stuff. No flags can be set.
The output of this script is pdb trajectories of normal modes 1-12 and a txt file with results (overlap etc.) for each mutant.'''

'''this version doesnt accept user input (pdb,chain)'''

class vibration:

    '''This class can be instantiated/initialized from other modules. blah blah some more documentation'''


    def overlap_calculation(self, eigenvectors_perturbed, eigenvectors_nonperturbed, cluster, eigenvalues_perturbed, eigenvalues_nonperturbed):

        ## calculate fewer overlaps to speed up things. all overlaps only calculated to determine mode of max overlap
        ## reduce moderanges to speed up things significantly

        import math, Numeric
        cluster.sort()
        cluster.reverse()
        overlaps = []
        max_overlaps = []
        perturbed_modes_of_max_overlap = []
        delta_perturbed_eigenvalues_of_max_overlap = []
        moderanges_nonperturbed = range(6,18)+[len(eigenvectors_nonperturbed)-3*len(cluster)-1]
        moderanges_perturbed = range(6,12)+[len(eigenvectors_nonperturbed)-3*len(cluster)-1]
        for mode_nonperturbed in range(len(eigenvectors_nonperturbed)):
            max_overlap = 0
            overlaps_per_mode_perturbed = []
            ## convert eigenvector of nonperturbed structure and delete appropiate coordinates
            vector_nonperturbed = list(eigenvectors_nonperturbed[mode_nonperturbed])
            for i in cluster:
                for j in range(3-1,-1,-1):
                    del vector_nonperturbed[3*i+j]
            overlaps_per_mode_perturbed = [0 for mode_perturbed in range(len(eigenvectors_perturbed))]
            for mode_perturbed in range(len(eigenvectors_perturbed)):
                if mode_nonperturbed in moderanges_nonperturbed and mode_perturbed in moderanges_perturbed:
                    ## calculate overlap between eigenvector of nonperturbed and perturbed structure
                    overlap = math.fabs(self.cosangle(vector_nonperturbed, eigenvectors_perturbed[mode_perturbed]))
                else:
                    overlap = 0
                ## append overlap to list of overlaps
                overlaps_per_mode_perturbed[mode_perturbed] = overlap
                if mode_nonperturbed == mode_perturbed:
                    overlaps.append(overlap)
            ## identify max overlap in list of overlaps
            max_overlap = max(overlaps_per_mode_perturbed)
            max_overlaps.append(max_overlap)
            ## identify mode of max overlap
            perturbed_mode_of_max_overlap = overlaps_per_mode_perturbed.index(max_overlap)
            perturbed_modes_of_max_overlap.append(perturbed_mode_of_max_overlap+1)
            ## identify eigenvalue of mode of max overlap
            delta_perturbed_eigenvalue_of_max_overlap = eigenvalues_perturbed[perturbed_mode_of_max_overlap]-eigenvalues_nonperturbed[mode_nonperturbed]
            delta_perturbed_eigenvalues_of_max_overlap.append(delta_perturbed_eigenvalue_of_max_overlap)

        return overlaps, max_overlaps, perturbed_modes_of_max_overlap, delta_perturbed_eigenvalues_of_max_overlap


    def faculty(self, n):
        faculty = 1
        for i in range(1, n+1):
            faculty = faculty*i
        return faculty


    def alignment_sequential_scoring(self, pdb1ATOM_all, pdb2ATOM_all, chain1, chain2):

        pdb1ATOM_chain = pdb1ATOM_all[chain1]['residues']
        pdb2ATOM_chain = pdb2ATOM_all[chain2]['residues']
        seq1,seq2,seq1seq2equiv,seq2seq1equiv = self.alignment_sequential_core(pdb1ATOM_chain, pdb2ATOM_chain)

        if (seq1.count('-')/float(len(seq1)) > 0.80 or seq2.count('-')/float(len(seq2)) > 0.80):
            
            count_of_nonaligned_residues = (len(seq1)-seq1.count('-'))+(len(seq2)-seq2.count('-'))
            count_of_nonidentical_aligned_residues = (len(seq1)-seq1.count('-'))+(len(seq2)-seq2.count('-'))

        else:
            
            for i in range(len(seq1seq2equiv)-1,-1,-1): ## remove '-' twice
                if seq1seq2equiv[i] != '-':
                    break
                del seq1seq2equiv[i]
            seq1seq2equiv.reverse()
            for i in range(len(seq1seq2equiv)):
                if seq1seq2equiv[i] != '-':
                    break
                del seq1seq2equiv[i]
            for i in range(len(seq2seq1equiv)-1,-1,-1):
                if seq2seq1equiv[i] != '-':
                    break
                del seq2seq1equiv[i]
            seq2seq1equiv.reverse()
            for i in range(len(seq2seq1equiv)):
                if seq2seq1equiv[i] != '-':
                    break
                del seq2seq1equiv[i]
            seq1seq2equiv.reverse()
            seq2seq1equiv.reverse()

            count_of_nonaligned_residues = seq1seq2equiv.count('-')+seq2seq1equiv.count('-')

            count_of_nonidentical_aligned_residues = 0
            count_of_missing_residues = 0
            if len(seq1) != len(seq2):
                notexpected
            alignment_len = len(seq1)
            for residue in range(alignment_len): ## foer loop skal startpunktet for identiske residues findes... fundet fordi terminale '-' slettes?!
                if (
                    seq1[residue] != seq2[residue] and 
                    (pdb1ATOM_chain.has_key(residue+1+min(seq2seq1equiv)) and pdb2ATOM_chain.has_key(residue+1+min(seq1seq2equiv)))
                    ):
                    count_of_nonidentical_aligned_residues += 1

            residues1 = pdb1ATOM_chain.keys()
            residues1.sort()
            residues2 = pdb2ATOM_chain.keys()
            residues2.sort()

            addition = residues1[min(seq2seq1equiv)]-residues2[min(seq1seq2equiv)]

            pdb2ATOM_all[chain2]['alignment'] = addition

        return count_of_nonaligned_residues, count_of_nonidentical_aligned_residues, pdb2ATOM_all


    def chain_combinator(self, chain_combinations, monomers):
        combination=0
        while len(chain_combinations) > combination:
            for i in range(monomers):
                for j in range(self.faculty(i-1)):
                    chain_combinations[k][1].append (chain_combinations[k][2][i])
                    chain_combinations[k][2] = chain_combinations[k][2][:i]+chain_combinations[k][2][i+1:]
                    combination += 1
        return chain_combinations, i-1


    def parse_pdb(self, pdb_lines, pdb_chain, pdb_model):

        '''
        pdb parser based on "PDB Format Description Version 2.2" that
        1a) from pdb only parses information from the coordinate section and parsing thus
            initiates when first MODEL or ATOM is met and
            terminates when first ENDMDL or last ATOM respectively is met
        1b) from pdb coordinate section only parses information from the records; model, atom, ter, endmdl (not hetatm)
        2 ) from pdb coordinate section only parses atoms of any specified chain and model
        3 ) from atom records only parses; atom name (13-16), alternate location indicator (17),
            residue name (18-20), residue sequence number (23-26),
            x- (31-38), y- (39-46), z-coordinates (47-54)
        4 ) renumbers residues as specified by user
        '''

##        print 'parsing info from coordinate section about residues, atoms and coordinates'

        import Numeric, sets

        dictionary_of_coordinates = {}
        parse = True
        REMARK350transformations = {}
        REMARK350record = 'biomt'
        REMARK350biomolecule = 1
        REMARK465 = {}
        resolution = 0.
        compndchains = []
        HEADERdepDate = 'N/A'
        helices = {}
        strands = {}
        modres = {}

        for pdb_line in pdb_lines:
            record = pdb_line[:6].strip()

            if record == 'HEADER':
                HEADERdepDate = pdb_line[50:59]

            elif record == 'COMPND':
                if pdb_line[11:16] == 'CHAIN':
                    chains = []
                    for i in range(18,80,3):
                        if pdb_line[i:i+4] == 'NULL':
                            chains.append(' ')
                        elif pdb_line[i-3:i+1] == 'NULL':
                            continue
                        elif pdb_line[i] == ' ':
                            compndchains.append(chains)
                            break
                        else:
                            chains.append(pdb_line[i])

            elif record == 'REMARK':

                remark = int(pdb_line[6:10].strip())

                if remark == 2:
                    if pdb_line[11:22] == 'RESOLUTION.' and pdb_line[28:38] == 'ANGSTROMS.':
                        resolution = float(pdb_line[22:27])

                if remark == 350:
                    if pdb_line[10:22] == ' BIOMOLECULE':
                        REMARK350biomolecule = int(pdb_line[23:80].strip()) ## 24:25 i 1H4V
                    if pdb_line[10:40] == ' APPLY THE FOLLOWING TO CHAINS':
                        if REMARK350record == 'biomt':
                            REMARK350chains = []
                        REMARK350record = 'chains'
                        for i in range(42,80,3):
                            if pdb_line[i:i+4] == 'NULL':
                                REMARK350chains.append(' ')
                            elif pdb_line[i-3:i+1] == 'NULL':
                                continue
                            elif pdb_line[i] == ' ':
                                break
                            else:
                                REMARK350chains.append(pdb_line[i])
                                
                    elif pdb_line[10:18] == '   BIOMT':
                        REMARK350record = 'biomt'
                        if pdb_line[18:19] == '1':
                            REMARK350m = Numeric.zeros((3,3),typecode='d')
                            REMARK350v = []
                        REMARK350m[int(pdb_line[18:19])-1][0] = float(pdb_line[24:33])
                        REMARK350m[int(pdb_line[18:19])-1][1] = float(pdb_line[34:43])
                        REMARK350m[int(pdb_line[18:19])-1][2] = float(pdb_line[44:53])
                        if pdb_line[54:68] == '              ':
                            REMARK350v.append(0)
                        else:
                            REMARK350v.append(float(pdb_line[54:68]))
                        if pdb_line[18:19] == '3':
                            if not REMARK350transformations.has_key(REMARK350biomolecule):
                                REMARK350transformations[REMARK350biomolecule] = []
                            REMARK350transformations[REMARK350biomolecule].append([REMARK350chains, REMARK350m, REMARK350v])

                    elif REMARK350record == 'chains':
                        for i in range(12,80,3): ## does chains always continue from column 12? use str.split instead...
                            if pdb_line[i] == ' ':
                                break
                            else:
                                REMARK350chains.append(pdb_line[i])
                                    
                if remark == 465 and pdb_line[11:12] == ' ' and pdb_line[15:18] in self.resdic.keys():
                    print 'warning, warning, some residues are missing'
                    chain = pdb_line[19:20]
                    resno = int(pdb_line[21:26])
                    resname = pdb_line[15:18]
                    if not REMARK465.has_key(chain):
                        REMARK465[chain] = [[resno, resname]]
                    else:
                        REMARK465[chain].append([resno, resname])

##                if pdb_line[6:10].strip() == '470' and pdb_line[11:12] == ' ' and pdb_line[15:18] in self.resdic.keys():
##                    print 'warning, warning, some residues are missing'
##                for atom in self.atoms_hessian:
##                    if ' '+atom+' ' in pdb_line[27:80]:
##                        if not missing_atoms.has_key(pdb_line[19:20]):
##                            missing_atoms[pdb_line[19:20]] = {}
##                        if not missing_atoms[pdb_line[19:20]].has_key(int(pdb_line[21:24])):
##                            missing_atoms[pdb_line[19:20]][int(pdb_line[21:24])] = {}
##                        if not missing_atoms[pdb_line[19:20]][int(pdb_line[21:24])].has_key('atoms'):
##                            missing_atoms[pdb_line[19:20]][int(pdb_line[21:24])]['resname'] = pdb_line[15:18]
##                            missing_atoms[pdb_line[19:20]][int(pdb_line[21:24])]['atoms'] = [atom]
##                        else:
##                            missing_atoms[pdb_line[19:20]][int(pdb_line[21:24])]['atoms'].append(atom)

            elif record == 'MODRES':
                chain = pdb_line[16]
                resno = int(pdb_line[18:22])
                resname = pdb_line[24:27]
                if not modres.has_key(chain):
                    modres[chain] = {}
                modres[chain][resno] = resname
                continue

            elif record == 'HELIX':
                chain = pdb_line[19]
                if not helices.has_key(chain):
                    helices[chain] = [[int(pdb_line[21:25]), int(pdb_line[33:37])]]
                else:
                    helices[chain].append([int(pdb_line[21:25]), int(pdb_line[33:37])])
                continue

            elif record == 'SHEET':
                chain = pdb_line[21]
                if not strands.has_key(chain):
                    strands[chain] = [[int(pdb_line[22:26]), int(pdb_line[33:37])]]
                else:
                    strands[chain].append([int(pdb_line[22:26]), int(pdb_line[33:37])])
                continue

            elif record == 'MODEL':
                current_model = pdb_line[10:14].strip()
                if not pdb_model:
                    pdb_model = pdb_line[10:14].strip()
                    parse = True
                elif int(pdb_line[10:14].strip()) not in pdb_model:
                    parse = False
                elif int(pdb_line[10:14].strip()) in pdb_model:
                    parse = True
                continue

            elif record == 'ATOM':
                if not parse == True:
                    continue
                resname = pdb_line[17:20].strip()
                if resname not in self.resdic.keys():
                    continue
                chain = pdb_line[21]
                if not (pdb_chain == [] or chain in pdb_chain):
                    continue
                dictionary_of_coordinates = self.parse_atom(pdb_line, dictionary_of_coordinates)
                continue

            elif record == 'HETATM':
                if not parse == True:
                    continue
                chain = pdb_line[21]
                if not ((pdb_chain == [] or chain in pdb_chain)):
                    continue
                if not modres.has_key(chain):
                    continue
                resno = int(pdb_line[22:26].strip())
                if not modres[chain].has_key(resno):
                    continue
                dictionary_of_coordinates = self.parse_atom(pdb_line, dictionary_of_coordinates)
                continue

        for biomolecule in REMARK350transformations.keys():
            for transformation in range(len(REMARK350transformations[biomolecule])-1,-1,-1):
                chains = list(sets.Set(REMARK350transformations[biomolecule][transformation][0]) & sets.Set(dictionary_of_coordinates.keys()))
                if chains != []:
                    REMARK350transformations[biomolecule][transformation][0] = chains
                else:
                    notexpected
                    REMARK350transformations[biomolecule].remove(transformation)

        for chain in REMARK465.keys():
            for residue in range(len(REMARK465[chain])):
                dictionary_of_coordinates[chain]['residues'][REMARK465[chain][residue][0]] = {'res_name': REMARK465[chain][residue][1], 'atoms': {}}

        return dictionary_of_coordinates, HEADERdepDate, REMARK350transformations, resolution, compndchains, helices, strands


    def parse_atom(self, pdb_line, dictionary_of_coordinates, res_name = None):

        '''This func parses coordinates from 'pdb_line' and returns "dictionary_of_coordinates"'''

        import Numeric
        chain = pdb_line[21]
        res_no = int(pdb_line[22:26].strip())
        if not res_name:
            res_name = pdb_line[17:20].strip()
        atom_name = pdb_line[12:16].strip()
##        atom_altloc = pdb_line[16].strip()
        atom_x_coord = float(pdb_line[30:38].strip())
        atom_y_coord = float(pdb_line[38:46].strip())
        atom_z_coord = float(pdb_line[46:54].strip())
        coordinate = Numeric.array([atom_x_coord, atom_y_coord, atom_z_coord])
## der skal kunne skelnes mellem pdb1 og pdb2...
        if (self.residue_terminal_c[chain] != None and res_no > self.residue_terminal_c[chain]) or (self.residue_terminal_n[chain] != None and res_no < self.residue_terminal_n[chain]):
            print self.residue_terminal_c[chain], res_no, self.residue_terminal_n[chain]
            notexpected
            return dictionary_of_coordinates
        if not dictionary_of_coordinates.has_key(chain):
            dictionary_of_coordinates[chain] = {'alignment': 0, 'transformant': pdb_line[21], 'residues': {}}
        if not dictionary_of_coordinates[chain]['residues'].has_key(res_no):
            dictionary_of_coordinates[chain]['residues'][res_no] = {'res_name':res_name, 'atoms':{atom_name:{'coordinates':coordinate}}}
        elif not dictionary_of_coordinates[chain]['residues'][res_no]['atoms'].has_key(atom_name):
            dictionary_of_coordinates[chain]['residues'][res_no]['atoms'][atom_name] = {'coordinates':coordinate}

        return dictionary_of_coordinates


    def alignment_sequential_core(self, pdb1ATOM_chain, pdb2ATOM_chain):

## first part prepares sequences
        ## get residue range (twice)
        residues1 = pdb1ATOM_chain.keys()
        residues2 = pdb2ATOM_chain.keys()
        residues1.sort()
        residues2.sort()
        
        import sys
        sys.path.append('../') ## local implicit path for superpos (windows)
        import superpos

        ## write sequence (twice)
        seq1 = ''
        for residue in residues1:
            res_name = pdb1ATOM_chain[residue]['res_name']
            res_abbr = self.resdic[res_name]
            seq1 += res_abbr
        seq2 = ''
        for residue in residues2:
            res_name = pdb2ATOM_chain[residue]['res_name']
            res_abbr = self.resdic[res_name]
            seq2 += res_abbr

        instance_NW = superpos.NW(seq1,seq2)
        seq1,seq2,seq1seq2equiv,seq2seq1equiv = instance_NW.Align()

        return seq1,seq2,seq1seq2equiv,seq2seq1equiv


    def morph(self, eigenvectors_all_modes, nframes, biomolecule, pdb1ATOM_all, atoms_hessian, cluster, matrix_hessian, job, cutoff, fileprefix):

        '''This function puts in frames between the maximum amplitudes of a
movement given by eigenvectors. The values of the diagonal elements of the
hessian matrix are used for B factors to make coloring during simulation
possible.'''

##        print 'visualize the two extreme projections along a trajectory and interpolate n frames between them'

        import time, math

        bfactors = self.parse_connectivity(matrix_hessian)

        for mode in range(6,12):

            eigenvectors = eigenvectors_all_modes[mode]

##            eigenvectors = length_adjustment(mode, eigenvectors)

            output_vmd = ['REMARK color by connectivity (b-factor) or squared displacement (temperature factor)\n']
            for frame in range(nframes):
                output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
                i = 0
                for chain in range(len(biomolecule[0])): ## flyt denne loekke til loop der bruges inden hessian_calc og angiv med flag, om der skal morphes eller ej ## denne loekke maa kun optraede en gang ## skriv for h... en dickey, der angiver om atomet skal morphes eller ej!!!!!!!!!!
                    chain1 = biomolecule[0][chain]
                    for residue in pdb1ATOM_all[chain1]['residues'].keys():
                        if residue in cluster:
                            continue
                        for atom in pdb1ATOM_all[chain1]['residues'][residue]['atoms'].keys():
                            if atom in atoms_hessian: ## skriv den latterlige key morph til pdb1coords dic'en!!!!!!!!!!!!
                                coordinate = pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['coordinates']
                                res_name = pdb1ATOM_all[chain1]['residues'][residue]['res_name']
                                x1 = 9*eigenvectors[i+0]
                                y1 = 9*eigenvectors[i+1]
                                z1 = 9*eigenvectors[i+2]
                                sqlength = x1**2+y1**2+z1**2
                                x2 = coordinate[0]+(1-2*float(frame)/float(nframes))*x1
                                y2 = coordinate[1]+(1-2*float(frame)/float(nframes))*y1
                                z2 = coordinate[2]+(1-2*float(frame)/float(nframes))*z1
                                    
                                output_vmd.append('ATOM        %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(atom, res_name, chain1, int(residue), x2,y2,z2, bfactors[i/3], sqlength))
                                i += 3

                output_vmd.append('TER\nENDMDL\n')
            fd = open(fileprefix+'mode'+str(mode+1)+'.pdb', 'w') ## implicit path
            fd.writelines(output_vmd)
            fd.close()

        return


    def hessian_calculation(self, coordinates, cutoff, verbose = False):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

        if verbose == True:
            print 'calculating the hessian/second-order derivate Kirchoff/Laplacian matrix'
        
        import math, Numeric, time

        cutoff_sq = cutoff**2

        N = (len(coordinates))

        matrix_hessian = Numeric.zeros((3*N,3*N), typecode='d')
        
        for row_sup in range(N):
            for col_sup in range(N):
                if col_sup > row_sup:
                    #does the Numeric module feature some smart built-in function to calculate length of vectors? use math.sqrt(math.pow(vector, 2))
                    xi = coordinates[row_sup][0]
                    xj = coordinates[col_sup][0]
                    yi = coordinates[row_sup][1]
                    yj = coordinates[col_sup][1]
                    zi = coordinates[row_sup][2]
                    zj = coordinates[col_sup][2]
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


    def eigenv_calccomb(self, matrix_hessian, jobid, verbose):

        '''Calculates eigenvectors and eigenvalues of a matrix.'''

        if verbose == True:
            print 'calculating eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'
        
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
        if verbose == True:
            lines = ['rows=modes, cols=coordinates\n']
            for mode in range(6,len(eigenvectors)):
                lines += [str(eigenvectors[mode])+'\n']
            fd = open('%s_eigenvectors.txt' %(jobid),'w')
            fd.writelines(lines)
            fd.close()

        import math, Numeric

        ## calculate length of mode 7 (equals 1 when using module linearalgebra)
        len7 = math.sqrt(sum(Numeric.array(eigenvectors[6])*Numeric.array(eigenvectors[6])))

        ## loop over modes
        for i in range(7,len(eigenvalues)-1):

            ## calculate length of mode i
            leni = math.sqrt(sum(Numeric.array(eigenvectors[i])*Numeric.array(eigenvectors[i])))

            ## scale length of mode i relative to length of mode 7
            lenfactor = (len7/leni)/(eigenvalues[i]/eigenvalues[6])
            for j in range(len(eigenvectors[i])):
                eigenvectors[i][j] *= lenfactor

        ## copy lists of eigenvectors to eigenvectors_combined
        eigenvectors_combined = []
        for mode in range(len(eigenvalues)):
            eigenvectors_combined.append(list(eigenvectors[mode]))
        ## change mode i to be the sum of modes 7 to i
        for mode in range(7,len(eigenvalues)):
            for coordinate in range(len(eigenvalues)):
                eigenvectors_combined[mode][coordinate] += eigenvectors_combined[mode-1][coordinate]
        
##        print 'calculated eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'

        return eigenvectors, eigenvalues, eigenvectors_combined


    def pdb_import(self, pdb_structure):

        '''This function imports a pdb file and searches in the following order: local, Conway PDB, PDB'''

        print 'importing pdb %s' %pdb_structure
        
        import os, urllib2
        if not os.path.isfile(pdb_structure+'.pdb'):
            if os.path.isfile('/data/remediated_pdb/%s/pdb%s.ent' %(pdb_structure.lower()[1:3], pdb_structure.lower())):
                fd = open('/data/remediated_pdb/%s/pdb%s.ent' %(pdb_structure.lower()[1:3], pdb_structure.lower()), 'r')
                pdb_lines = fd.readlines()
                fd.close()
            else:
##                raise('pdb file %s not found' %pdb_structure)
                url = urllib2.urlopen('http://www.pdb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId='+pdb_structure)
                pdb_lines = url.readlines()
        else:
            fd = open(pdb_structure+'.pdb', 'r')
            pdb_lines = fd.readlines()
            fd.close()
        return pdb_lines


    def align_structural(self, chains1, chains2, pdb1ATOM_all, pdb2ATOM_all):

        print 'aligning chains %s and %s' %(chains1, chains2)
        
        import sys
        sys.path.append('../EAT_DB/')
        sys.path.append('../')
        import Protool
        instance_geometry=Protool.geometry()

        if len(chains1) != len(chains2):
            raise('hvad skal alignes?')

        coordinates1 = []
        coordinates2 = []
        for chain in range(len(chains1)):
            chain1 = chains1[chain]
            chain2 = chains2[chain]
            alnseqadd = pdb1ATOM_all[chain1]['alignment']-pdb2ATOM_all[chain2]['alignment']
            for residue in pdb1ATOM_all[chain1]['residues'].keys():
                if residue+alnseqadd in pdb2ATOM_all[chain2]['residues'].keys():
                    for atom in pdb1ATOM_all[chain1]['residues'][residue]['atoms'].keys():
                        if (
                            atom in pdb2ATOM_all[chain2]['residues'][residue+alnseqadd]['atoms'].keys() and
                            (
                                (pdb1ATOM_all[chain1]['residues'][residue]['res_name'] == pdb2ATOM_all[chain2]['residues'][residue+alnseqadd]['res_name']) or
                                (pdb1ATOM_all[chain1]['residues'][residue]['res_name'] != 'GLY' and pdb2ATOM_all[chain2]['residues'][residue+alnseqadd]['res_name'] != 'GLY' and atom in ['C','CA','N','O','CB']) or
                                (atom in ['C','CA','N','O'])
                                )
                            ):
                            coordinates1.append(pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['coordinates'])
                            coordinates2.append(pdb2ATOM_all[chain2]['residues'][residue+alnseqadd]['atoms'][atom]['coordinates'])

        r,t,rmsd,transv1,transv2=instance_geometry.superpose(coordinates1, coordinates2)

        return rmsd,r,t,transv1,transv2


    def REMARK350transformation_application(self, REMARK350transformations, compndchains, inputchains, pdbcoords):

        print 'applying transformation matrix from remark350 to coordinates of chains %s %s' %(compndchains, inputchains)

        import Numeric, sets

        chains_biomolecules = []
        for biomolecule in REMARK350transformations.keys():
            chains_biomolecule = []
            chains = sets.Set(self.chains)
            for transformation in REMARK350transformations[biomolecule]:
                chains = chains-sets.Set(transformation[0])
            for transformation in range(len(REMARK350transformations[biomolecule])):
                matrix = REMARK350transformations[biomolecule][transformation][1]
                vector = REMARK350transformations[biomolecule][transformation][2]
                if matrix != Numeric.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]]) or vector != [0,0,0]:
                    for transformchain_input in REMARK350transformations[biomolecule][transformation][0]:
                        if transformchain_input not in inputchains:
                            continue
                        
                        transformchain_output = chains.pop()

##                        REMARK350transformations[biomolecule].append(
##                            [
##                                [transformchain_output],
##                                Numeric.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]]),
##                                [0,0,0],
##                                ]
##                            )

                        for compndchain in range(len(compndchains)):
                            if transformchain_input in compndchains[compndchain]:
                                compndchains[compndchain].append(transformchain_output)

                        inputchains.append(transformchain_output)
                        if transformchain_input not in chains_biomolecule:
                            chains_biomolecule.append(transformchain_input)
                        chains_biomolecule.append(transformchain_output)

                        pdbcoords[transformchain_output] = {'alignment': pdbcoords[transformchain_input]['alignment'], 'residues': {}}
                        for residue in pdbcoords[transformchain_input]['residues'].keys(): ## split til buildREMARK350transform
                            if not pdbcoords[transformchain_output]['residues'].has_key(residue):
                                pdbcoords[transformchain_output]['residues'][residue] = {'res_name': pdbcoords[transformchain_input]['residues'][residue]['res_name'], 'atoms': {}}
                            for atom in pdbcoords[transformchain_input]['residues'][residue]['atoms'].keys():
                                coordinate = Numeric.matrixmultiply(matrix, pdbcoords[transformchain_input]['residues'][residue]['atoms'][atom]['coordinates']) + vector
                                pdbcoords[transformchain_output]['residues'][residue]['atoms'][atom] = {'coordinates': coordinate}

                else:
                    for transformchain_input in REMARK350transformations[biomolecule][transformation][0]:
                        if transformchain_input in inputchains:
                            chains_biomolecule.append(transformchain_input)

            chains_biomolecules.append(chains_biomolecule)

        return compndchains, inputchains, pdbcoords, chains_biomolecules


    def identify_identical_chains_intrapdb(self, pdbATOM_all, chains):

        if len(chains) == 1:
            chains_identical = [chains]
        else:
            chains_pool = list(chains)
            chains_identical = []
            i = -1
            for row in range(len(chains)):
                for col in range(len(chains)):
                    if row < col and chains[col] in chains_pool:
                        if len(pdbATOM_all[chains[row]]['residues'].keys()) >= 3 and len(pdbATOM_all[chains[col]]['residues'].keys()) >= 3:
                            count_of_nonaligned_residues, count_of_nonidentical_residues = self.alignment_sequential_scoring(pdbATOM_all, pdbATOM_all, chains[row], chains[col])[0:2]
                            if count_of_nonidentical_residues == 0:
                                if chains[row] in chains_pool:
                                    i += 1
                                    chains_identical.append([chains[row]])
                                    chains_pool.remove(chains[row])
                                chains_identical[i].append(chains[col])
                                chains_pool.remove(chains[col])

            if chains_identical == []:
                for chain in chains:
                    chains_identical.append([chain])

        return chains_identical


    def identify_identical_chains_interpdb(self, pdb1ATOM_all, pdb2ATOM_all, chains_identical_1, chains_identical_2, mutations):

        print 'determining identical chains between %s and %s' %(chains_identical_1, chains_identical_2)
        chains_identical = []
        i = -1
        for row in range(len(chains_identical_1)):
            for col in range(len(chains_identical_2)):
                if row <= col:
                    if not(len(pdb1ATOM_all[chains_identical_1[row][0]]['residues'].keys()) >= 3 and len(pdb2ATOM_all[chains_identical_2[col][0]]['residues'].keys()) >= 3):
                        notexpected
                    count_of_nonaligned_residues, count_of_nonidentical_residues = self.alignment_sequential_scoring(pdb1ATOM_all, pdb2ATOM_all, chains_identical_1[row][0], chains_identical_2[col][0])[0:2]
                    if count_of_nonidentical_residues <= mutations: ## allow fewer number of mutations?
                        chains_identical.append([chains_identical_1[row],chains_identical_2[col]])

        ## skal der printes besked om, hvad der ikke er identisk? i hvert fald hvis intet er identisk! return error to server!!!

        return chains_identical


    def identify_equivalent_chains_among_identical_chains(self, pdb1ATOM_all, pdb2ATOM_all, monomers1, monomers2):

        import math
        
        equivalent_chains = [[],[]]

##########################
##        biomolecules = []
##        for i in range(len(chains_identical)):
##            if len(chains_identical[i][0]) == len(chains_identical[i][1]):
##                biomolecules.append(chains_identical[i])
##            elif (
##                math.fmod(
##                    len(chains_identical[i][0]),
##                    len(chains_identical[i][1])
##                    ) == 0
##                or
##                math.fmod(
##                    len(chains_identical[i][1]),
##                    len(chains_identical[i][0])
##                    ) == 0
##                ):
##                a=1
##            else:
##                a=1
##
##        chains_equivalent = []
##        for i in range(len(chains_identical)):
##            if len(chains_identical[i][0]) > 5:
##                if sets.Set(chains_identical[i][0]) ^ sets.Set(identical_chains[1]) != sets.Set():
##                    raise('please specify equivalent chains')
##                else:
##                    identical_chains = chains_identical.pop(i)
##                    identical_chains[0].sort()
##                    identical_chains[1].sort()
##                    chains_equivalent.append([identical_chains[0],identical_chains[1]])
##            else:
##                chains_equivalent.append(self.identify_equivalent_chains_among_identical_chains(pdb1ATOM_all, pdb2ATOM_all))
###########################

        chains = len(monomers1)
        chain_combinations = []

        if chains > 5:
            monomers1.sort()
            monomers2.sort()
            if monomers1 == monomers2:
                equivalent_chains[0] += monomers1
                equivalent_chains[1] += monomers2
            else:
                raise('equivalent chains could not be determined. please specify chain input.')

        else:

            for combination in range(self.faculty(chains)):
                chain_combinations.append([list(monomers1),[],list(monomers2)])

            while chains > 0:
                combination = 0
                while combination < len(chain_combinations):
                    for chain in range(chains):
                        for chain_minus_1_combination in range(self.faculty(chains-1)):
                            chain_combinations[combination][1].append(chain_combinations[combination][2].pop(chain))
                            combination += 1
                chains -= 1

            rmsds = []
            for combination in range(len(chain_combinations)):
                rmsd = self.align_structural(chain_combinations[combination][0], chain_combinations[combination][1], pdb1ATOM_all, pdb2ATOM_all)[0]
                rmsds.append(rmsd)
    
            combination = rmsds.index(min(rmsds))
            equivalent_chains[0] += (chain_combinations[combination][0])
            equivalent_chains[1] += (chain_combinations[combination][1])

        return equivalent_chains


    def cosangle(self, v1, v2):
        ## Numeric arrays are not used, because they are slow!
        import math, Numeric
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


    def distance_residue_intra(self, coordinates1):
        import Numeric
        matrix_distance_residue_intra = Numeric.zeros((len(coordinates1),len(coordinates1)),typecode='d') ## change coordinates1 to onyl contain ca-atoms or some solution  like that
        for rowres in range(len(matrix_distance_residue_intra)):
            for colres in range(len(matrix_distance_residue_intra)):
                if colres > rowres:
                    v = coordinates1[rowres]-coordinates1[colres]
                    len_sq = v[0]**2+v[1]**2+v[2]**2
                    matrix_distance_residue_intra[rowres][colres] = len_sq
                    matrix_distance_residue_intra[colres][rowres] = len_sq
        return matrix_distance_residue_intra


    def parse_connectivity(self, matrix_hessian):

        ## instead of connectivity then the contact number could be used...

        connectivity = []
        for i in range(0,len(matrix_hessian),3):
            connectivity.append(
                matrix_hessian[i+0][i+0]+
                matrix_hessian[i+1][i+1]+
                matrix_hessian[i+2][i+2]
                )

        return connectivity


    def topology(self, filename, chain, biomolecule, helices, strands, axlen, calctype):

        import os

        if calctype == 2:
            bl = [390,833] ## bottom left coordinate
            tr = [1050,173] ## top right coordinate
            s = 14 ## space between plot and topology
        if calctype == 1:
            bl = [208,833] ## bottom left coordinate
            tr = [1232,173] ## top right coordinate
            s = 22 ## space between plot and topology
        w = (tr[0]-bl[0])/axlen ## width of squares
        secelms = [helices, strands]

        lines = 'convert '+filename
        for chain in biomolecule:
            for i in range(len(secelms)):
                if i == 0 and secelms[i].has_key(chain):
                    lines += ' -fill "rgb(0,0,255)" '
                if i == 1 and secelms[i].has_key(chain):
                    lines += ' -fill "rgb(255,0,0)" '
                if secelms[i].has_key(chain):
                    for secelmrange in secelms[i][chain]:
                        ## bottom
                        lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), bl[1]+s, int(bl[0]+w*secelmrange[1]), bl[1]+s)
                        ## top
                        lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), tr[1]-s, int(bl[0]+w*secelmrange[1]), tr[1]-s)
                        if calctype == 2:
                            ## left
                            lines += ' -draw "line %s,%s %s,%s" ' %(bl[0]-s, int(bl[1]-w*(secelmrange[0]-1)), bl[0]-s, int(bl[1]-w*secelmrange[1]))
                            ## right
                            lines += ' -draw "line %s,%s %s,%s" ' %(tr[0]+s, int(bl[1]-w*(secelmrange[0]-1)), tr[0]+s, int(bl[1]-w*secelmrange[1]))
            if helices.has_key(chain) or strands.has_key(chain):
                lines += filename
                os.system(lines)
                
        return


    def plot(self, jobid, biolunitchains, cutoff_distance, eigenvectors, chain, biomolecule, helices1, strands1, fileprefix, calctype, modes, ymax):

        import math, Numeric

        for mode in range(6,12)+[len(eigenvectors)-1]:

            if modes == 'single' and mode == len(eigenvectors)-1:
                continue

            ## calculate residue fluctuation (atom,RMSF) for each mode

            lengths = []
            for i in range(0,len(eigenvectors[mode]),3):
                length = math.sqrt(
                    eigenvectors[mode][i+0]**2+
                    eigenvectors[mode][i+1]**2+
                    eigenvectors[mode][i+2]**2
                    )
                lengths.append(length)

            ## write lengths to file

            lines = []
            for i in range(len(lengths)):
                lines.append(str(lengths[i])+'\n')
            fd = open(jobid+'_'+str(mode)+'_lengths.txt', 'w')
            fd.writelines(lines)
            fd.close()

            ## plot residue fluctuation (atom,RMSF) for each mode

            if modes == 'single':
                title1 = 'mode %s' %(mode+1)
            if modes == 'combo':
                title1 = 'modes 7-%s' %(mode+1)

            self.gnuplot_plot(
                jobid, cutoff_distance, biolunitchains, title1,
                xtitle = 'residue', ytitle = 'RMSF', plotname = 'residue displacements', filename = '%sdisplacement_mode%s' %(fileprefix, mode+1),
                data = lengths,
                )

            ## prepare data for cross correlation maps
            ccdata = Numeric.zeros((len(eigenvectors[mode])/3,len(eigenvectors[mode])/3), typecode='d')
            for ccrow in range(len(ccdata)):
                for cccol in range(len(ccdata)):
                    ccdata[ccrow][cccol] = self.cosangle(eigenvectors[mode][3*ccrow:3*ccrow+3], eigenvectors[mode][3*cccol:3*cccol+3])

            ## plot cross-correlation map (residue,residue,cross-correlation) for each mode
            for i in range(len(ccdata)):
                ccdata[i][i] = 1
            self.gnuplot_splot(
                2, jobid, cutoff_distance, biolunitchains, title1,
                len(ccdata),
                xtitle = 'residue', ytitle = 'residue', ztitle = 'cross correlation',
                plotname = 'cross correlation map',
                filename = fileprefix+'crosscorrelation_mode%s' %(mode+1),
                data = ccdata,
                z1 = -1, z2 = 1,
                )

            ## add topology to margin of cross-correlations maps
            self.topology(
                fileprefix+'crosscorrelation_mode'+str(mode+1)+'.png',
                chain, biomolecule, helices1, strands1, len(ccdata), 2,
                )

        return


    def loop_axis_x(self, res1, jobid, winsize, datadic, fileprefix, statusdir):

        filename = '%s%ssize%s_residue%s.txt' %(statusdir, fileprefix, winsize, res1+1)
        lines = []
        for key in datadic:
            lines.append('data '+key+'\n')
            data = datadic[key]['data']
            for mode in range(6,12):
                lines.append('mode '+str(mode+1)+'\n')
                for res2 in range(len(data[mode][res1])):
                    lines.append(str(data[mode][res1][res2])+'\n')
            if key == 'overlaps_combined':
                mode = len(data)-1
                lines.append('mode '+str(mode+1)+'\n')
                for res2 in range(len(data[mode][res1])):
                    lines.append(str(data[mode][res1][res2])+'\n')
        fd = open(filename, 'w')
        fd.writelines(lines)
        fd.close()

        return


    def loop_axis_y(
        self, pdb1ATOM_all, coordinates1, cutoff_distance,
        jobid, biolunitchains, chain, biomolecule, helices1, strands1, frames,
        atoms_hessian, HEADERdepDate1, REMARK2resolution1,
        xvalue, yvalue, calctype, ## dependent on uni or duo cluster
        coordinates_cluster, cluster,
        results, eigenvectors_nonperturbed,
        winsize, datadic_main, eigenvalues_nonperturbed,
        verbose,
        path_html = None, path_python = None
        ):

        if path_html:
            import sys
            sys.path.append('/var/www/cgi-bin/goodvibes/python/')
##            sys.path.append('/var/www/cgi-bin/goodvibes/python/goodvibes/')
            import time, server
            t1 = time.clock()
    
        matrix_hessian_perturbed = self.hessian_calculation(coordinates_cluster, float(cutoff_distance), verbose) ## calculate with coordinates2 as well... and compare results to switching pdb1 and pdb2; actually a slower step than eigenvec calc... optimize!
        eigenvectors_perturbed, eigenvalues_perturbed, eigenvectors_perturbed_combined = self.eigenv_calccomb(matrix_hessian_perturbed, jobid, verbose)
        overlaps_single, max_overlaps_single, perturbed_modes_of_max_overlap_single, delta_perturbed_eigenvalues_of_max_overlap = self.overlap_calculation(eigenvectors_perturbed, eigenvectors_nonperturbed, cluster, eigenvalues_perturbed, eigenvalues_nonperturbed)
        overlaps_combined, max_overlaps_combined, perturbed_modes_of_max_overlap_combined = self.overlap_calculation(eigenvectors_perturbed_combined, eigenvectors_nonperturbed, cluster, eigenvalues_perturbed, eigenvalues_nonperturbed)[:-1]

##        ## do vmd of perturbed structure
##        self.morph(eigenvectors, frames, biomolecule, pdb1ATOM_all, atoms_hessian, cluster, matrix_hessian, jobid+'-'+str(xvalue)+'-'+str(yvalue), biolunitchains, cutoff_distance)
##
##        ## write eigenvalues to file
##        fd = open('eigenvalues.txt', 'a')
##        fd.write('%13s ' %str(cluster))
##        for mode in range(6,18):
##            fd.write('%8.5f ' %eigenvalues[mode])
##        fd.write('\n')
##        fd.close()

        datadic_loop = {
            'eigenvalues_perturbed': eigenvalues_perturbed,
            'emo': delta_perturbed_eigenvalues_of_max_overlap,
            'overlaps_single': overlaps_single,
            'overlaps_max': max_overlaps_single,
            'mmo': perturbed_modes_of_max_overlap_single,
            'overlaps_combined': overlaps_combined
            }
        for key in datadic_main:
            for mode in range(6,12):
                datadic_main[key]['data'][mode][xvalue][yvalue] = datadic_loop[key][mode]
                if calctype == 2:
                    datadic_main[key]['data'][mode][yvalue][xvalue] = datadic_loop[key][mode]
        datadic_main['overlaps_combined']['data'][-1][xvalue][yvalue] = overlaps_combined[-1]
        if calctype == 2:
            datadic_main['overlaps_combined']['data'][-1][yvalue][xvalue] = overlaps_combined[-1]

        ## calculate remaining time based only on time spent on diagonalization (excluding finalizing steps with graphics generation)
        if path_html:
            residuecount = len(coordinates_cluster)+len(cluster)
            if calctype == 2:
                rem_col = (residuecount-winsize-2.*(winsize-1)/2.) - (xvalue+1) + (winsize-1)/2.
                rem_matrices_col = (rem_col/2.)*(1+rem_col)
                rem_matrices_row = residuecount - (yvalue+1) - (winsize-1)/2.
            if calctype == 1:
                rem_col = residuecount/2.-(xvalue+1)
                rem_matrices_col = rem_col*residuecount
                rem_matrices_row = residuecount-(yvalue+1)
            remaining_matrices = rem_matrices_col + rem_matrices_row
            t2 = time.clock()
            rem_time = (t2-t1)*remaining_matrices
            rem_days = int( rem_time/float(60*60*24) )
            rem_hours = int( rem_time/float(60*60) - rem_days*24 )
            rem_minutes = int( rem_time/float(60) - rem_days*24*60 - rem_hours*60 )
            rem_seconds = int( rem_time - rem_days*24*60*60 - rem_hours*60*60 - rem_minutes*60 )
            if rem_time == 0:
                htmlbody = 'Dear user. Your input has been processed. Graphics and html is being created. Your results will be here in less than a minute. Thank you for using GoodVibes.'
            else:
                rem_time = '%s day(s), %s hour(s), %s minute(s), %s second(s)' %(rem_days, rem_hours, rem_minutes, rem_seconds)
                htmlbody = 'Dear user. Your input is being processed. Your results will be here in approximately %s.' %rem_time
            instance_server = server.server()
            instance_server.slowhtml(htmlbody, jobid, path_results = path_html)

        return datadic_main, results

    def cluster1(self, remres1, coordinates1, jobid, cutoff_distance, datadic, clustersize, matrix_distance_residue_intra):

        distances = list(matrix_distance_residue_intra[remres1])
        distances_sorted = list(distances)
        distances_sorted.sort()

        cluster = []
        for res in range(0,clustersize+1):
            cluster.append(distances.index(distances_sorted[res]))
        cluster.sort()

        coordinates_cluster = list(coordinates1)
        for res in range(len(cluster)-1,-1,-1):
            del coordinates_cluster[cluster[res]]

        return coordinates_cluster, cluster

    def cluster2(self, remres1, coordinates1, jobid, cutoff_distance, winsize, datadic, remres2, residue_max):

        import math

        ## continue if windows overlap with each other
        if math.fabs(remres1 - remres2) < winsize:
            return None, None
        ## continue if windows overlap with borders of the matrix
        if remres1 < (winsize-1)/2 or remres1 > (residue_max-1)-(winsize-1)/2 or remres2 < (winsize-1)/2 or remres2 > (residue_max-1)-(winsize-1)/2:
            return None, None
        ## continue if overlap values have already been calculated for 
        if remres2 < remres1:
            return None, None
        cluster = []
        for rest in range((winsize-1)/2,-(winsize-1)/2-1,-1):
            cluster.append(remres1+rest)
            cluster.append(remres2+rest)
        coordinates_cluster = list(coordinates1)
        for rest in range((winsize-1)/2,-(winsize-1)/2-1,-1):
            del coordinates_cluster[remres2+rest]
        for rest in range((winsize-1)/2,-(winsize-1)/2-1,-1):
            del coordinates_cluster[remres1+rest]

        return coordinates_cluster, cluster


    def main(
        self, pdb1lines, chains1, pdb1model, atoms_hessian, quarternary, jobid, frames,
        cutoff_distances, calctype, winsizes,
        pdb2lines = None, chains2 = None, pdb2model = None,
        path_html = None, path_python = None, verbose = False, statusdir = ''):

        '''This is the main function that calls all the other functions.'''

        import os, urllib2, Numeric, sys, sets, math, time
        import goodvibes_core

        results = []

        pdb1ATOM_all, HEADERdepDate1, REMARK350transformations1, REMARK2resolution1, COMPNDchains1, helices1, strands1 = self.parse_pdb(pdb1lines, chains1, pdb1model)

        ## assume multimeric biological unit if chains not specified by user
        if chains1 == []:
            chains1 = pdb1ATOM_all.keys()

#### apply REMARK350transformations to get new chains
##        chains_identical_1 = []
##        if REMARK350transformations1 != {}:
##            print 'building biomolecules of pdb1\n'
##            COMPNDchains1, chains1, pdb1ATOM_all, chains_biomolecules_1 = self.REMARK350transformation_application(REMARK350transformations1, COMPNDchains1, chains1, pdb1ATOM_all)
##            for biomolecule in range(len(chains_biomolecules_1)):
##                chains_identical_1.append(self.identify_identical_chains_intrapdb(pdb1ATOM_all, chains_biomolecules_1[biomolecule]))
##        else:
##            chains_identical_1.append(self.identify_identical_chains_intrapdb(pdb1ATOM_all, chains1))

        ## append individual chains to biomolecules (monomeric biological units)
        biomolecules = []
        for chain in chains1:
            biomolecules.append([chain])

        ## loop over biomolecules
        for biomolecule in biomolecules: ##biomolecules = chains eller first chain if chains not given by user

            ## join values from the list biomolecule to the string biolunitchains
            biolunitchains = ''.join(biomolecule).replace(' ','-')

            ## parse coordinates from the dictionary pdb1ATOM_all to the list coordinates1
            ## and change the flag 'hessian' to 'true' for each atom in the dictionary
            coordinates1 = []
            vectors_difference = []
            for chain1 in biomolecule:
                for residue in pdb1ATOM_all[chain1]['residues'].keys():
                    for atom in pdb1ATOM_all[chain1]['residues'][residue]['atoms'].keys():
                        if atom in atoms_hessian:
                            coordinates1.append(pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['coordinates'])
                            pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['hessian'] = 1

            ## abort if too few residues ## move this error message to server.py before sending mail to user
            if len(coordinates1) < 4:
                raise('too few coordinates in pdb file. invalid input.')

            ## calculate intra-residue distances
            matrix_distance_residue_intra = self.distance_residue_intra(coordinates1)

            ## set ymax
            residue_max = len(coordinates1)
            if calctype == 2:
                ymax = residue_max
            if calctype == 1:
                ymax = 10

            ## loop over cutoffs
            for cutoff_distance in cutoff_distances:

                fileprefix = '%s_chains%s_cutoff%s_' %(jobid, biolunitchains, int(cutoff_distance))

                ## for the non-disrupted structure: calculate and visualize eigenvectors
                matrix_hessian = self.hessian_calculation(coordinates1, float(cutoff_distance), verbose) ## calculate with coordinates2 as well... and compare results to switching pdb1 and pdb2
                eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_combined_nonperturbed = self.eigenv_calccomb(matrix_hessian, jobid, verbose)
                self.morph(eigenvectors_nonperturbed, frames, biomolecule, pdb1ATOM_all, atoms_hessian, [], matrix_hessian, jobid, cutoff_distance, fileprefix)

                ## do plots for individual modes
                self.plot(jobid, biolunitchains, cutoff_distance, eigenvectors_nonperturbed, chain, biomolecule, helices1, strands1, fileprefix, calctype, 'single', ymax)
                ## do plots for combined modes
                self.plot(jobid, biolunitchains, cutoff_distance, eigenvectors_combined_nonperturbed, chain, biomolecule, helices1, strands1, fileprefix+'combo_', calctype, 'combo', ymax)
                xxx

                results.append(
                    {
                        'number of non-zero eigenvalues': len(eigenvalues_nonperturbed)-6,
                        'depdate1': HEADERdepDate1,
                        'cutoff': cutoff_distance,
                        'res1': REMARK2resolution1,
                        'biomolecules': biomolecule, 'chains1': pdb1ATOM_all.keys(),
                        }
                    )

################################################################################
##                     plotdata initialization section                        ##
################################################################################

                N = len(matrix_hessian)/3
                datadic = goodvibes_core.datadic_return(N)

################################################################################
##                        plotdata calculation section                        ##
################################################################################

                for winsize in winsizes:

                    ## loop over remres1
                    for remres1 in range(residue_max):

                        ## continue the remres1 loop if another processor is running calculations for that value of remres1
                        filename = '%s%ssize%s_residue%s.txt' %(statusdir, fileprefix, winsize, remres1+1)
                        if os.path.isfile(filename):
                            ## read datafile to data variables before you continue
                            fd = open(filename, 'r')
                            lines_main = fd.readlines()
                            fd.close()
                            line_number_main = 0
                            for line_main in lines_main:
                                if lines_main[line_number_main][:4] == 'data':
                                    line_number_data = line_number_main+2
                                    key = line_main[5:-1]
                                elif line_main[:4] == 'mode':
                                    line_number_data = line_number_main+1
                                    mode = int(line_main[5:-1])-1
                                elif line_number_main >= line_number_data:
                                    remres2 = line_number_main - line_number_data
                                    datadic[key]['data'][mode][remres1][remres2] = float(line_main)
                                    if calctype == 2:
                                        datadic[key]['data'][mode][remres2][remres1] = float(line_main)
                                line_number_main += 1

                            continue
                        
                        ## create data files (necessary to create data files if multithreading/parallel processing dependent on existence of files to determine if calculation on certain residues has been performed)
                        else:
                            fd = open(filename, 'w')
                            fd.close()

                        if calctype == 2:
                            
                            ## loop over remres2
                            for remres2 in range(residue_max):

                                if not path_html:
                                    print 'residue combination %s,%s of %s,%s' %(remres1+1,remres2+1,residue_max,residue_max)

                                coordinates_cluster, cluster = self.cluster2(remres1, coordinates1, jobid, cutoff_distance, winsize, datadic, remres2, residue_max)
                                if not cluster:
                                    continue
                                else:
                                    ## calculate data
                                    datadic, results = self.loop_axis_y(
                                        pdb1ATOM_all, coordinates1, cutoff_distance,
                                        jobid, biolunitchains, chain,
                                        biomolecule, helices1, strands1, frames,
                                        atoms_hessian, HEADERdepDate1, REMARK2resolution1,
                                        remres1, remres2, calctype, ## unique for calctype 2
                                        coordinates_cluster, cluster,
                                        results, eigenvectors_nonperturbed,
                                        winsize, datadic, eigenvalues_nonperturbed,
                                        verbose,
                                        path_html, path_python,
                                        )

                        if calctype == 1:

                            ## loop over clustersizes
                            for clustersize in range(10):

                                coordinates_cluster, cluster = self.cluster1(remres1, coordinates1, jobid, cutoff_distance, datadic, clustersize, matrix_distance_residue_intra)
                                if not path_html:
                                    print 'residue %s of %s and %s surrounding residues %s' %(remres1+1, residue_max, clustersize, cluster)
                                ## append cluster to a list of clusters and check if calculations already performed for this cluster
                                ## if calculations already performed then copy appropriate data

                                datadic, results = self.loop_axis_y(
                                    pdb1ATOM_all, coordinates1, cutoff_distance,
                                    jobid, biolunitchains, chain,
                                    biomolecule, helices1, strands1, frames,
                                    atoms_hessian, HEADERdepDate1, REMARK2resolution1,
                                    remres1, clustersize, calctype, ## unique for calctype 1
                                    coordinates_cluster, cluster,
                                    results, eigenvectors_nonperturbed,
                                    winsize, datadic, eigenvalues_nonperturbed,
                                    path_html, path_python,
                                    )

                        ## write data to txt files in case of crash during loop over residues/clusters
                        self.loop_axis_x(
                            remres1, jobid, winsize, datadic, fileprefix, statusdir)

################################################################################
##                                plot section                                ##
################################################################################
                                
                axistitle = 'central residue in removed cluster'
                if calctype == 2:
                    ## set the diagonal elements to a standard value to get full contour
                    for mode in range(6,12):
                        for key in datadic.keys():
                            diagonal = datadic[key]
                            if datadic[key]['diagonal'] == 'max':
                                diagonal = max(datadic[key]['data'][mode][0])
                            elif datadic[key]['diagonal'] == 'mode':
                                diagonal = mode+1
                            else:
                                diagonal = datadic[key]['diagonal']
                            for i in range((winsize-1)/2,residue_max-(winsize-1)/2):
                                for j in range(winsize):
                                    for k in range(winsize):
                                        datadic[key]['data'][mode][i-(winsize-1)/2+j][i-(winsize-1)/2+k] = diagonal
                    ytitle = axistitle
                if calctype == 1:
                    ytitle = 'number of residues in removed cluster'

                ## plot the data with the adjusted values of the diagonal elements

##                ## overlaps_combined between nonperturbed and perturbed modes 7-3N
##                self.gnuplot_splot(
##                    calctype, jobid, cutoff_distance, biolunitchains, 'modes 7-3N',
##                    ymax,
##                    xtitle = axistitle, ytitle = axistitle, ztitle = '{/Symbol D}overlap',
##                    plotname = '{/Symbol D}overlap',
##                    filename = fileprefix+'overlapscombined_mode3N',
##                    data = resres_overlaps_combined[-1],
##                    z1 = '', z2 = '',
##                    )
##
##                ## add topology to margin of cross-correlations maps
##                self.topology(
##                    fileprefix+'overlapscombined_mode3N.png',
##                    chain, biomolecule, helices1, strands1, residue_max, calctype,
##                    )
                
                for mode in range(6,12):

                    for key in datadic:

                        self.gnuplot_splot(
                            calctype, jobid, cutoff_distance, biolunitchains, mode+1,
                            ymax,
                            xtitle = axistitle, ytitle = axistitle, ztitle = datadic[key]['name'],
                            plotname = datadic[key]['name'],
                            filename = '%s%s_mode%s' %(fileprefix, key, mode+1),
                            data = datadic[key]['data'][mode],
                            z1 = datadic[key]['zrange'][0], z2 = datadic[key]['zrange'][1],
                            )

                        ## add topology to margin of cross-correlations maps
                        self.topology(
                            '%s%s_mode%s.png' %(fileprefix, key, mode+1),
                            chain, biomolecule, helices1, strands1, residue_max, calctype,
                            )

##                    ## perturbed eigenvalue of 6 perturbed modes
##                    data = resres_eigenvalues_perturbed

##                print 'end of loop over cutoff distance %s' %cutoff_distance

            ## end of loop over biomolecules

##        print 'done'
        return results


    def gnuplot_plot(self, jobid, cutoff_distance, chains, title1, xtitle, ytitle, plotname, filename, data):

        import os

        ## write gnuplot data to txt file
        lines = []
        lines.append('%4i %16.13f\n' %(0, 0))
        for i in range(len(data)):
            lines.append('%4i %16.13f\n' %(i+1, data[i]))
        lines.append('%4i %16.13f\n' %(len(data)+1, 0))
        fd = open('%s.gnuplotdata' %(filename), 'w')
        fd.writelines(lines)
        fd.close()
        ## write gnuplot settings to txt file
        lines = [
            'set terminal postscript eps enhanced color "Helvetica" 48\n',
            'set output "%s.ps"\n' %(filename),
            'set size 4,4\n', ## scale 400%
            'set autoscale fix\n', ## scale axes
            'set encoding iso_8859_1\n', ## postscript encoding for special characters
            'set style data histeps\n', ## change plot style to histogram
            'set style function histeps\n', ## change plot style to histogram
            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" offset 0,1\n' %(plotname, title1, jobid, chains, cutoff_distance),
            'set xlabel "%s"\n' %(xtitle),
            'set ylabel "%s"\n' %(ytitle),
            'plot "%s.gnuplotdata" title "" lt 3 lw 16\n' %(filename), ## plot gnuplot data file
            ## set ticslevel
            ]

        fd = open('%s.gnuplotsettings' %(filename), 'w')
        fd.writelines(lines)
        fd.close()
        ## plot data with gnuplot plot
        os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(filename))
        ## convert postscript to portable network graphics
        os.system('convert %s.ps %s.png' %(filename, filename))
        os.remove('%s.ps' %(filename))
        os.remove('%s.gnuplotdata' %(filename))
        os.remove('%s.gnuplotsettings' %(filename))

        return


    def gnuplot_splot(self, calctype, jobid, cutoff_distance, chains, title1, ymax, xtitle, ytitle, ztitle, plotname, filename, data, z1, z2):

        import os

        ## write gnuplot data to txt file
        gnuplot_splot_data = []
        for x in range(len(data)):
            for y in range(ymax):
                if calctype == 2:
                    if x >= y:
                        z = data[x][y]
                    if y > x:
                        z = data[y][x]
                if calctype == 1:
                    z = data[x][y]
                gnuplot_splot_data.append('%4i %4i %16.13f\n' %(x+1, y+1, z))
            gnuplot_splot_data.append('%4i %4i %16.13f\n\n' %(x+1, ymax+1, data[0][0]))
        for i in range(ymax+1):
            gnuplot_splot_data.append('%4i %4i %16.13f\n' %(len(data)+1, i+1, data[0][0]))
        fd = open('%s.gnuplotdata' %(filename), 'w')
        fd.writelines(gnuplot_splot_data)
        fd.close()
        ## write gnuplot settings to txt file
        lines = []
        if calctype == 2:
            lines += ['set size square\n'] ## scale square
        lines += [
            'set terminal postscript eps enhanced color "Helvetica" 48\n',
            'set output "%s.ps"\n' %(filename),
            'set size 4,4\n', ## scale 400%
            'set view map\n', ## change orientation of plot
            'set autoscale fix\n', ## scale axes
            'set style data pm3d\n', ## set by default?
            'set style function pm3d\n', ## set by default?
            'set encoding iso_8859_1\n', ## postscript encoding for special characters
            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" offset 0,1\n' %(plotname, title1, jobid, chains, cutoff_distance),
            'set xlabel "%s"\n' %(xtitle),
            'set ylabel "%s"\n' %(ytitle),
            'set cblabel "%s"\n' %(ztitle),
            'set cbrange [%s:%s]\n' %(z1,z2), ## set colorbox range
            'set palette model CMY rgbformulae 7,5,15\n',
            'set pm3d map corners2color c1\n', ## generate a 2D surface rather than 3D points
            'splot "%s.gnuplotdata" title ""\n' %(filename), ## splot gnuplot data file
            ## set ticslevel
            ]

        fd = open('%s.gnuplotsettings' %(filename), 'w')
        fd.writelines(lines)
        fd.close()
        ## plot data with gnuplot splot
        os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(filename))
        ## convert postscript to portable network graphics
        os.system('convert %s.ps %s.png' %(filename, filename))
        os.remove('%s.ps' %(filename))
        os.remove('%s.gnuplotdata' %(filename))
        os.remove('%s.gnuplotsettings' %(filename))

        return

    def sigmoid(self, x, cutoff, slope=1):
        import math
        y = 1. / ( 1. + math.exp( slope*(x-cutoff) ) )
        return y

    def __init__(self):
        self.residue_terminal_n = {'A':-9999, 'B':-9999, 'C':-9999, 'D':-9999, 'E':-9999, 'F':-9999, 'G':-9999, 'H':-9999, 'I':-9999, 'J':-9999, 'K':-9999, 'L':-9999, 'M':-9999, 'N':-9999, 'O':-9999, 'P':-9999, 'Q':-9999, 'R':-9999, 'S':-9999, 'T':-9999, 'U':-9999, 'V':-9999, 'W':-9999, 'X':-9999, 'Y':-9999, 'Z':-9999, ' ':-9999,}
        self.residue_terminal_c = {'A':9999, 'B':9999, 'C':9999, 'D':9999, 'E':9999, 'F':9999, 'G':9999, 'H':9999, 'I':9999, 'J':9999, 'K':9999, 'L':9999, 'M':9999, 'N':9999, 'O':9999, 'P':9999, 'Q':9999, 'R':9999, 'S':9999, 'T':9999, 'U':9999, 'V':9999, 'W':9999, 'X':9999, 'Y':9999, 'Z':9999, ' ':9999,}
        ##
        self.chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'] ## used for remark350build
        self.resdic = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y','PCA':'X','ACE':'X'}
        self.resdic20 = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y'}
        self.weekdaydic = {0:'Monday',1:'Tuesday',2:'Wednesday',3:'Thursday',4:'Friday',5:'Saturday',6:'Sunday'}
        self.monthdic = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
        '''Wow, this script implements a phantom spring/elastic network model.
        The paramater "atoms" is used to select, which atom vectors are used when calculating an average vector of a residue. i.e. used in trjconv_parsing and diff_vec_calc'''

if __name__=='__main__':

    import os

    results = []

    instance_vibration = vibration()

    pdb = ['1zvq', ['A']]
    pdblines = instance_vibration.pdb_import(pdb[0])
    chain = pdb[1]
    job = pdb[0]
    statusdir = os.getcwd()+'/'

    results = instance_vibration.main(
        pdblines, chain, [1], ['CA'],
        'monomeric', job, 50,
        [10], 2, [1], verbose = True, statusdir = statusdir
        )

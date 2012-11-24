#!/bin/env /software/bin/python2.3
#
# $Id: goodvibes_windows.py 26 2006-06-27 12:24:01Z tc $
#
## make sure that missing atoms and residues don't go silent by
## hvad skal der ske med atomer med altloc?
## add info E.C. # from pdb.org or pdb (REMARK2?! eller COMPND...)
## return error if a user specified chain is not found in pdb input file
## use info from SEQADV lines if modified residue and MODRES lines

##        import time
##        t1 = time.clock()
##        t2 = time.clock()
##        print t2-t1

'''This documentation refers to version 1.55. This module contains only one class. The module is used to do some NMA stuff. No flags can be set.
The output of this script is pdb trajectories of normal modes 1-12 and a txt file with results (overlap etc.) for each mutant.'''

class vibration:

    '''This class can be instantiated/initialized from other modules. blah blah some more documentation'''

##    def overlap_calculation_old(self, eigenvectors, eigenvalues, vectors_difference):
##
##        import math, Numeric
##        overlaps_per_mode = []
##        max_overlap = 0
##        for mode in range(len(eigenvectors)):
##            eigenvectors_per_mode = eigenvectors[mode]
##            overlaps_per_atom = []
##            for atom in range(len(vectors_difference)): ## paabegynd kun loop hvis len(eigenvectors_per_mode) == 3*len(vectors_difference)
##                sqv1 = eigenvectors_per_mode[3*atom+0]**2+eigenvectors_per_mode[3*atom+1]**2+eigenvectors_per_mode[3*atom+2]**2
##                sqv2 = Numeric.dot(vectors_difference[atom], vectors_difference[atom])
##                dot_product = (
##                    eigenvectors_per_mode[3*atom+0]*vectors_difference[atom][0]+
##                    eigenvectors_per_mode[3*atom+1]*vectors_difference[atom][1]+
##                    eigenvectors_per_mode[3*atom+2]*vectors_difference[atom][2]
##                    )
##                if math.sqrt((sqv1*sqv2)) == 0:
##                    continue
##                overlap_per_atom = dot_product/math.sqrt((sqv1*sqv2))
##                overlaps_per_atom.append(overlap_per_atom)
##            overlap_per_mode = math.fabs(sum(overlaps_per_atom)/len(overlaps_per_atom))
##            overlaps_per_mode.append(overlap_per_mode)
##            if overlap_per_mode > max_overlap:
##                max_overlap = overlap_per_mode
##                mode_of_max_overlap = mode
##
##        eigenvalue_of_max_overlap = eigenvalues[mode_of_max_overlap]
##        
##        return overlaps_per_mode, max_overlap, mode_of_max_overlap, eigenvalue_of_max_overlap
##
##    def dislin_graf3_hessian(self, cutoff, job, chains, matrix_hessian):
##
##        import os, math, sys
##        sys.path.append('/software/lib/dislin/python')
##        import dislin
##
##        N = len(matrix_hessian)/3
##
##        chains = ''.join(chains).replace(' ','-')
##
##        for row in range(len(matrix_hessian)):
##            if 'max_val' not in vars().keys() or 'min_val' not in vars().keys():
##                max_val = max(matrix_hessian[row])
##                min_val = min(matrix_hessian[row])
##            elif max(matrix_hessian[row]) > max_val:
##                max_val = max(matrix_hessian[row])
##            elif min(matrix_hessian[row]) > min_val:
##                min_val = min(matrix_hessian[row])
##
##        zmat = []
##        for row in range(len(matrix_hessian)):
##            for col in range(len(matrix_hessian[row])):
##                zmat.append(matrix_hessian[row][col])
##
##        dislin.metafl ('png')
##        dislin.scrmod ('revers')
##        dislin.disini ()
##
##        dislin.pagera () ## plots a border around the page
##        dislin.bmpfnt ('SIMPLEX')
##        dislin.titlin (job, 1)
##        dislin.titlin ('chains '+chains, 3)
##        dislin.titlin ('cut-off distance '+str(round(float(cutoff)/10,2))+' nm', 4)
##
##        dislin.name   ('second derivative', 'Z')
##
##        dislin.autres (3*N, 3*N)
##        dislin.axspos (300, 1850)
##        dislin.ax3len (9*N, 9*N, 9*N)
##        
##        dislin.graf3  (1, 3*N, 1, 3*N-1,
##                       1, 3*N, 1, 3*N-1,
##                       min_val, max_val, min_val, max_val-min_val)
##        dislin.crvmat (zmat, 3*N, 3*N, 1, 1)
##        dislin.height (36)
##        dislin.title  ()
##        dislin.disfin ()
##
##        os.rename('dislin.png', job+'_'+chains+'_'+cutoff+'_hessian.png')
##
##        return

    def overlap_calculation(self, eigenvectors, eigenvalues, eigenvectors_clustersize0, cluster):

        ## calculate overlap between eigenvectors of harmonic/Hooke/vibrational elastic block/Calpha phantom/Gaussian network model and difference vectors of transition between two states
##        print 'calculating overlap between ENM eigenvectors and transition/difference vectors'

        import math, Numeric
        overlaps_per_mode_clustersize0 = []
        max_overlaps = []
        modes_of_max_overlap = []
        eigenvalues_of_max_overlap = []
        for mode_clustersize0 in range(len(eigenvectors_clustersize0)):
            max_overlap = 0
            overlaps_per_mode = []
            for mode in range(6,12):
                residue = 0
                overlaps_per_residue = []
                for residue_clustersize0 in range(len(eigenvectors_clustersize0[mode_clustersize0])/3):
                    if residue_clustersize0 not in cluster:
                        v1 = eigenvectors_clustersize0[mode_clustersize0][3*residue_clustersize0+0:3*residue_clustersize0+3]
                        v2 = eigenvectors[mode][3*residue+0:3*residue+3]
                        if v1 != [0,0,0] and v2 != [0,0,0]:
                            overlap_per_residue = self.cosangle(v1, v2)
                            overlaps_per_residue.append(overlap_per_residue)
                        residue += 1
                overlap_per_mode = math.fabs(sum(overlaps_per_residue)/len(overlaps_per_residue))
                overlaps_per_mode.append(overlap_per_mode)
                max_overlap = max(overlaps_per_mode)
                mode_of_max_overlap = overlaps_per_mode.index(max_overlap)
                eigenvalue_of_max_overlap = eigenvalues[mode_of_max_overlap]
            overlaps_per_mode_clustersize0.append(overlaps_per_mode)
            max_overlaps.append(max_overlap)
            modes_of_max_overlap.append(mode_of_max_overlap)
            eigenvalues_of_max_overlap.append(eigenvalue_of_max_overlap)

        return overlaps_per_mode_clustersize0, max_overlaps, modes_of_max_overlap, eigenvalues_of_max_overlap


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
                elif pdb_line[10:14].strip() not in pdb_model:
                    parse = False
                elif pdb_line[10:14].strip() in pdb_model:
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


    def morph(self, eigenvectors_all_modes, nframes, biomolecule, pdb1ATOM_all, atoms_hessian, amplitude_average, cluster, matrix_hessian, job, chains, cutoff):

        '''This function puts in frames between the maximum amplitudes of a
movement given by eigenvectors. The values of the diagonal elements of the
hessian matrix are used for B factors to make coloring during simulation
possible.'''

        print 'visualize the two extreme projections along a trajectory and interpolate n frames between them'

        import time, math

        bfactors = self.parse_connectivity(matrix_hessian)

##        lines_morph = []
        lengths = []
        for mode in range(6,12):
##            lines_morph.append([])
            lengths.append([])
            eigenvectors = eigenvectors_all_modes[mode]

            ## calculate average length
            lengthsum = 0
            for i in range(0,len(eigenvectors),3):
                length = math.sqrt(
                    eigenvectors[i+0]**2+
                    eigenvectors[i+1]**2+
                    eigenvectors[i+2]**2
                    )
                lengthsum += length
                lengths[mode-6].append(length)
            lengthaverage = lengthsum/(len(eigenvectors)/3)

            ## calculate standard deviation of length
            lengthdiffsq = 0
            for i in range(0,len(eigenvectors),3):
                lengthdiffsq += (
                    math.sqrt(
                        eigenvectors[i+0]**2+
                        eigenvectors[i+1]**2+
                        eigenvectors[i+2]**2
                        )-lengthaverage
                    )**2
            lengthstdev = math.sqrt(lengthdiffsq/(len(eigenvectors)/3))

            ## decrease or increase length if not in the range of average +/- 2 standard deviations
            lengthchange = False
            maxlen = lengthaverage+5*lengthstdev
            minlen = lengthaverage-5*lengthstdev
            for i in range(0,len(eigenvectors),3):
                length = math.sqrt(
                    eigenvectors[i+0]**2+
                    eigenvectors[i+1]**2+
                    eigenvectors[i+2]**2
                    )
                if length > maxlen:
                    eigenvectors[i+0] *= (maxlen)/length
                    eigenvectors[i+1] *= (maxlen)/length
                    eigenvectors[i+2] *= (maxlen)/length
                    lengthchange = True
                    print 'mode %s, residue %s, length %s > %s' %(mode+1, (i+3)/3, length, maxlen)
                if length < minlen:
                    eigenvectors[i+0] *= (minlen)/length
                    eigenvectors[i+1] *= (minlen)/length
                    eigenvectors[i+2] *= (minlen)/length
                    lengthchange = True
                    print 'mode %s, residue %s, length %s < %s' %(mode+1, (i+3)/3, length, minlen)

            if lengthchange:
                ## calculate new average length
                lengthsum = 0
                for i in range(0,len(eigenvectors),3):
                    lengthsum += math.sqrt(
                        eigenvectors[i+0]**2+
                        eigenvectors[i+1]**2+
                        eigenvectors[i+2]**2
                        )
                lengthaverage = lengthsum/(len(eigenvectors)/3)

                ## adjust new average length to average amplitude
                factor = (amplitude_average/lengthaverage)
                for i in range(len(eigenvectors)):
                    eigenvectors[i] *= factor

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
                                x1 = eigenvectors[i+0]
                                y1 = eigenvectors[i+1]
                                z1 = eigenvectors[i+2]
                                sqlength = x1**2+y1**2+z1**2
                                x2 = coordinate[0]+(1-2*float(frame)/float(nframes))*x1
                                y2 = coordinate[1]+(1-2*float(frame)/float(nframes))*y1
                                z2 = coordinate[2]+(1-2*float(frame)/float(nframes))*z1
                                    
##                                lines_morph[mode-6][frame].append('ATOM        %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(atom, res_name, chain1, int(residue), x2,y2,z2, bfactors[i/3], sqlength))

                                output_vmd.append('ATOM        %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(atom, res_name, chain1, int(residue), x2,y2,z2, bfactors[i/3], sqlength))
                                i += 3

                output_vmd.append('TER\nENDMDL\n')
            fd = open(job+'_'+chains+'_'+str(cutoff)+'_'+str(mode)+'.pdb', 'w') ## implicit path
            fd.writelines(output_vmd)
            fd.close()

        return lengths


    def hessian_calculation(self, coordinates, cutoff):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

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
                    if dist_sq <= cutoff_sq:
                        vector = [x,y,z]
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

    def function_eigenvector_calculation(self, matrix_hessian, amplitude_average, jobid):

        '''Calculates eigenvectors and eigenvalues of a matrix.'''

        print 'calculating eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'
        
        import LinearAlgebra
        eigen_tuple = LinearAlgebra.Heigenvectors(matrix_hessian)
        eigenvalues = list(eigen_tuple[0])
        eigenvectors = list(eigen_tuple[1])
        eigen_list = zip(eigenvalues, eigenvectors)
        eigen_list.sort()
        eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
        eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]

        import math

        for mode in range(len(eigenvectors)):
            if mode >= 6:

                lengthsum = 0
                for atom in range(0,len(eigenvectors[6]),3):
                    lengthsum += math.sqrt(
                        eigenvectors[mode][atom+0]**2+
                        eigenvectors[mode][atom+1]**2+
                        eigenvectors[mode][atom+2]**2)
                lengthaverage = lengthsum/(len(eigenvectors[mode])/3)

                factor = (amplitude_average/lengthaverage)
                for coordinate in range(len(eigenvectors[mode])):
                    eigenvectors[mode][coordinate] *= factor

##        ## write eigenvectors to file
##        for mode in range(6,12):
##            lines_eigenvectors = []
##            for coordinate in range(len(eigenvectors[mode])):
##                lines_eigenvectors.append(str(eigenvectors[mode][coordinate])+'\n')
##            fd = open('%s_eigenvectors_mode%s' %(jobid, mode), 'w')
##            fd.writelines(lines_eigenvectors)
##            fd.close()

        return eigenvectors, eigenvalues

    def pdb_import(self, pdb_structure):

        '''This function imports a pdb file and searches in the following order: local, Conway PDB, PDB'''

        print 'importing pdb %s' %pdb_structure
        
        import os, urllib2
        if not os.path.isfile(pdb_structure.lower()+'.pdb'):
            if os.path.isfile('/data/pdb/'+pdb_structure.lower()+'.pdb'):
                fd = open('/data/pdb/'+pdb_structure.lower()+'.pdb', 'r')
                pdb_lines = fd.readlines()
                fd.close()
            else:
                raise('pdb file %s not found' %pdb_structure)
                url = urllib2.urlopen('http://www.pdb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId='+pdb_structure)
                pdb_lines = url.readlines()
        else:
            fd = open(pdb_structure.lower()+'.pdb', 'r')
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
        import math, Numeric
        dot_product = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
        sqv1 = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]
        sqv2 = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]
        cosang = dot_product / math.sqrt(sqv1*sqv2)
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

    def identify_cluster(self, clustersize, residue_cluster_central, clusters, coordinates_clustersize0, matrix_distance_residue_intra):

        ## deletion of coordinates should not happen.instead coordinate sshould be skipped like eigenvectors if necessary.

        coordinates_cluster = list(coordinates_clustersize0)
        cluster = []
        ## delete coordinates from list of coordinates if within cluster
        if clustersize != 0:
            distances = list(matrix_distance_residue_intra[residue_cluster_central])
            distances_sorted = list(distances)
            distances_sorted.sort()
            for i in range(clustersize):
                cluster.append(distances.index(distances_sorted[i]))
            cluster.sort()
##            print 'clustersize', len(cluster), 'central residue', residue_cluster_central, 'cluster', cluster

            for i in range(len(cluster)-1, -1, -1):
                del coordinates_cluster[cluster[i]]

        return cluster, clusters, coordinates_cluster

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

    def dislin_graf3(self, xtitle, ytitle, ztitle, filename, data, title1 = '', title2 = '', title3 = '', title4 = '', ZOR = -1.):

        '''
        G R A F 3
        
        The routine GRAF3 plots a 3-D axis system where the Z-axis is plotted as a colour bar.

        The call is: CALL GRAF3 (XA, XE, XOR, XSTP, YA, YE, YOR, YSTP, ZA, ZE, ZOR, ZSTP) level 1 

        XA, XE are the lower and upper limits of the X-axis. 
        XOR, XSTP are the first X-axis label and the step between labels. 
        YA, YE are the lower and upper limits of the Y-axis. 
        YOR, YSTP are the first Y-axis label and the step between labels. 
        ZA, ZE are the lower and upper limits of the Z-axis. 
        ZOR, ZSTP are the first Z-axis label and the step between labels. 
        '''

        import os, math, sys
        sys.path.append('/software/lib/dislin/python')
        import dislin

        zmat = []
        for row in range(len(data)):
            for col in range(len(data[row])):
                zmat.append(data[row][col])

        ZE = max(zmat)
        ZA = min(zmat)
        ZSTP = float(ZE)/float(10**math.floor(math.log(ZE,10)))
        for step in [1,2,5,10]:
            if ZSTP <= step:
                ZSTP = float(step*10**(math.floor(math.log(ZE,10))-1))
                break

        XA = 1
        XE = len(data) ## use as input and get rid of clusterdata sentences in main...
        XSTP = float(XE)/float(10**math.floor(math.log(XE,10)))
        for step in [1,2,5,10]:
            if XSTP <= step:
                XSTP = float(step*10**(math.floor(math.log(XE,10))-1))
                break
        XOR = XSTP*math.ceil(float(XA)/float(XSTP))
        
        YA = 1
        YE = len(data[0])
        YSTP = float(YE)/float(10**math.floor(math.log(YE,10)))
        for step in [1,2,5,10]:
            if YSTP <= step:
                YSTP = int(float(step*10**(math.floor(math.log(YE,10))-1)))
                break
        YOR = YSTP*math.ceil(float(YA)/float(YSTP))

        max_val = max(zmat)
        min_val = min(zmat)

        dislin.metafl ('png')
        dislin.scrmod ('revers') ## black to white
##        dislin.winsiz (900*math.ceil(float(XE)/900.), 900*math.ceil(float(YE)/900.)) ##lvl0-4, window size in pixels (default 853x603)
##        dislin.page   (3000, 3000) ## lvl0, page size in plot coordinates (default 2970x2100)
        dislin.disini ()

        dislin.pagera () ## plots a border around the page
##        dislin.bmpfnt ('SIMPLEX')
        dislin.hwfont ()
        
        dislin.titlin (title1, 1)
        dislin.titlin (title2, 2)
        dislin.titlin (title3, 3)
        dislin.titlin (title4, 3)

        dislin.name   (xtitle, 'X') ## x-axis title
        dislin.name   (ytitle, 'Y') ## y-axis title
        dislin.name   (ztitle, 'Z')

        dislin.autres (XE, YE) ##  With a call to AUTRES, the size of coloured rectangles will be automatically calculated by GRAF3 or CRVMAT.
        dislin.axspos (300, 2100-300) ## lvl1, the position of the lower left corner, default is: horizontal (center), vertical (page height-300)
        dislin.ax3len (int((2./3.)*2100), int((2./3.)*2100), int((2./3.)*2100)) ##  The routine AX3LEN defines the axis lengths of a coloured axis system.
##        dislin.intax  () ## integer on all axes
        
        dislin.graf3  (XA, XE, XOR, XSTP,
                       YA, YE, YOR, YSTP,
                       ZA, ZE, ZOR, 0.5) ## lvl1, sets lvl3
        dislin.crvmat (zmat, XE, YE, 1, 1) ## CRVMAT plots a coloured surface according to a matrix
        dislin.height (36)
        dislin.title  ()
        dislin.disfin ()

        os.rename('dislin.png', filename)


        return

    def dislin_graf(self, xvalues, yvalues, filename, title1, title2, title3, title4, xtitle, ytitle, YE, YSTP, XE, quadratic = None):

        '''
        G R A F 

        GRAF plots a two-dimensional axis system.

        The call is: CALL GRAF (XA, XE, XOR, XSTP, YA, YE, YOR, YSTP) level 1 

        XA, XE are the lower and upper limits of the X-axis. 
        XOR, XSTP are the first X-axis label and the step between labels. 
        YA, YE are the lower and upper limits of the Y-axis. 
        YOR, YSTP are the first Y-axis label and the step between labels.
        '''

        import os, sys, math
        sys.path.append('/software/lib/dislin/python/')
        import dislin

        XSTP = float(max(xvalues))/float(10**math.floor(math.log(max(xvalues),10)))
        for step in [1,2,5,10]:
            if XSTP <= step:
                XSTP = float(step*10**(math.floor(math.log(max(xvalues),10))-1))
                break

        XOR = XSTP*math.ceil(float(min(xvalues))/XSTP)

        dislin.metafl ('png') ## lvl0, file output type
##        dislin.winsiz (900*math.ceil(float(XE)/float(900)), 900*math.ceil(float(YE)/float(900))) ##lvl0-4, window size in pixels (default is 853x603)
##        dislin.page   (2970, 2100) ## lvl0, page size in plot coordinates, default is 2970x2100
        dislin.scrmod ('revers') ## black to white

        dislin.disini () ## lvl0, change to lvl1
        
        dislin.pagera () ## lvl1-3, plots a border around the page
        dislin.bmpfnt ('SIMPLEX') ## font
        dislin.titlin (title1, 1) ## title line 1
        dislin.titlin (title2, 2) ## title line 2
        dislin.titlin (title3, 3) ## title line 3
        dislin.titlin (title4, 4) ## title line 4
        dislin.name   (xtitle, 'X') ## x-axis title
        dislin.name   (ytitle, 'Y') ## y-axis title
        dislin.incmrk (0) ## selects line and/or symbol mode for CURVE
        dislin.setgrf ('NAME', 'NAME', 'TICKS', 'TICKS') ## removes a part of an axis or a complete axis from an axis system
##        dislin.axspos (100, 900*math.ceil(float(YE)/float(900))+100) ## lvl1, the position of the lower left corner, default is: horizontal (center), vertical (page height-300)
##        dislin.axslen (900*math.ceil(float(XE)/float(900)), 900*math.ceil(float(YE)/float(900))) ## lvl1, the size of an axis system

        dislin.graf   (min(xvalues), XE, XOR, XSTP, 0, YE, 0, YSTP) ## lvl1, change to lvl2

        dislin.polcrv ('BARS') ## interpolation method
        dislin.curve  (xvalues, yvalues, len(xvalues))
        dislin.title  ()

        dislin.endgrf () ## change to lvl1 (only applicable if multiple axis systems on 1 page)

        dislin.disfin () ## change to lvl0

        os.rename('dislin.png', filename)

    def topology(self, filename, chain, biomolecule, helices1, strands1, ccdata):
        f1 = 2./(2970./853.+2100./603.)
        f2 = f1*1400./float(len(ccdata))
        lines = '/usr/bin/convert '+filename
        for chain in biomolecule:
            if helices1.has_key(chain):
                lines += ' -fill "rgb(0,0,255)" '
                for secelm in helices1[chain]:
                    ## bottom
                    lines += " -draw 'line %s,%s,%s,%s' " %(int(300*f1+f2*secelm[0]), 603-int(300*f1)+10, int(300*f1+f2*secelm[1]), 603-int(300*f1)+10)
                    ## top
                    lines += " -draw 'line %s,%s,%s,%s' " %(int(300*f1+f2*secelm[0]), 603-int((300+1400)*f1)-10, int(300*f1+f2*secelm[1]), 603-int((300+1400)*f1)-10)
                    ## left
                    lines += " -draw 'line %s,%s,%s,%s' " %(int(300*f1)-10, 603-int(300*f1+f2*secelm[0]), int(300*f1)-10, 603-int(300*f1+f2*secelm[1]))
                    ## right
                    lines += " -draw 'line %s,%s,%s,%s' " %(int((300+1400)*f1)+10, 603-int(300*f1+f2*secelm[0]), int((300+1400)*f1)+10, 603-int(300*f1+f2*secelm[1]))
            if strands1.has_key(chain):
                lines += ' -fill "rgb(255,0,0)" '
                for secelm in strands1[chain]:
                    ## bottom
                    lines += " -draw 'line %s,%s,%s,%s' " %(int(300*f1+f2*secelm[0]), 603-int(300*f1)+10, int(300*f1+f2*secelm[1]), 603-int(300*f1)+10)
                    ## top
                    lines += " -draw 'line %s,%s,%s,%s' " %(int(300*f1+f2*secelm[0]), 603-int((300+1400)*f1)-10, int(300*f1+f2*secelm[1]), 603-int((300+1400)*f1)-10)
                    ## left
                    lines += " -draw 'line %s,%s,%s,%s' " %(int(300*f1)-10, 603-int(300*f1+f2*secelm[0]), int(300*f1)-10, 603-int(300*f1+f2*secelm[1]))
                    ## right
                    lines += " -draw 'line %s,%s,%s,%s' " %(int((300+1400)*f1)+10, 603-int(300*f1+f2*secelm[0]), int((300+1400)*f1)+10, 603-int(300*f1+f2*secelm[1]))
            if helices1.has_key(chain) or strands1.has_key(chain):
                lines += filename
                os.system(lines)


    def main(self, pdb1lines, chains1, pdb1model, atoms_hessian, amplitude_average, cutoff_distances, mutations, quarternary, jobid, frames, pdb2lines = None, chains2 = None, pdb2model = None):

        '''This is the main function that calls all the other functions.'''

##        if atoms_hessian == []:
##            atoms_hessian = ['OXT']
##            for chemical_symbol in ['C','O','H','N','S']:
##                for remoteness_indicator in [' ','A','B','G','D','E','Z','H']:
##                    for branch_designator in range(10):
##                        atoms_hessian.append(chemical_symbol+remoteness_indicator+str(branch_designator))
##        if self.residue_terminal_1_n == None:
##            self.residue_terminal_n = {'A':-9999, 'B':-9999, 'C':-9999, 'D':-9999, 'E':-9999, 'F':-9999, 'G':-9999, 'H':-9999, 'I':-9999, 'J':-9999, 'K':-9999, 'L':-9999, 'M':-9999, 'N':-9999, 'O':-9999, 'P':-9999, 'Q':-9999, 'R':-9999, 'S':-9999, 'T':-9999, 'U':-9999, 'V':-9999, 'W':-9999, 'X':-9999, 'Y':-9999, 'Z':-9999, ' ':-9999,}
##        if self.residue_terminal_1_c == None:
##            self.residue_terminal_c = {'A':9999, 'B':9999, 'C':9999, 'D':9999, 'E':9999, 'F':9999, 'G':9999, 'H':9999, 'I':9999, 'J':9999, 'K':9999, 'L':9999, 'M':9999, 'N':9999, 'O':9999, 'P':9999, 'Q':9999, 'R':9999, 'S':9999, 'T':9999, 'U':9999, 'V':9999, 'W':9999, 'X':9999, 'Y':9999, 'Z':9999, ' ':9999,}

        import os, urllib2, Numeric, sys, sets, math

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

            biolunitchains = ''.join(biomolecule).replace(' ','-')

            coordinates1 = []
            vectors_difference = []
            for chain1 in biomolecule:
                for residue in pdb1ATOM_all[chain1]['residues'].keys():
                    for atom in pdb1ATOM_all[chain1]['residues'][residue]['atoms'].keys():
                        if atom in atoms_hessian:
                            coordinates1.append(pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['coordinates'])
                            pdb1ATOM_all[chain1]['residues'][residue]['atoms'][atom]['hessian'] = 1

            ## calculate intra-residue distances
            matrix_distance_residue_intra = self.distance_residue_intra(coordinates1)

            ## loop over cutoffs
            for cutoff_distance in cutoff_distances:

## loop over window residues (remember to fix 4 lines 200 lines below...)
##############################################################################
                clustersize = 1
                clustersize_max = 1
                residue_cluster_central = 1
                clusterplotdata = []
                winsize = 1

                for mode_clustersize0 in range(6,12):
                    clusterplotdata.append([])
                    for mode in range(6,12):
                        clusterplotdata[mode_clustersize0-6].append(
                            Numeric.zeros((len(coordinates1),len(coordinates1)),typecode='d')
                            )

                ## loop over remres1 and write to clusterplotdata
                for remres1 in range(len(coordinates1)):
                    print '----------------', remres1, '-----------------'
                    for mode_clustersize0 in range(6,12):
                        for mode in range(6,12):
                            fd = open('cluster_%s_%s_%s.txt' %(remres1, mode_clustersize0, mode), 'w')
                            fd.close()
                    for remres2 in range(len(coordinates1)):
                        print '----', remres2, '-----'

                        ## continue if windows overlap with each other
                        if math.fabs(remres1 - remres2) < winsize:
                            continue
                        ## continue if windows overlap with borders of the matrix
                        if remres1 < (winsize-1)/2 or remres1 > 163-(winsize-1)/2 or remres2 < (winsize-1)/2 or remres2 > 163-(winsize-1)/2:
                            continue
                        ## continue if overlap values have already been calculated for 
                        if remres2 < remres1:
                            for mode_clustersize0 in range(6,12):
                                for mode in range(6,12):
                                    if mode_clustersize0 == 6 and mode == 6:
                                        print remres1, remres2
                                        print clusterplotdata[mode_clustersize0-6][mode-6][remres2][remres1], clusterplotdata[mode_clustersize0-6][mode-6][remres1][remres2]
##                                    clusterplotdata[mode_clustersize0-6][mode-6][remres2][remres1] = clusterplotdata[mode_clustersize0-6][mode-6][remres1][remres2]
                            continue
                        cluster = []
                        for rest in range((winsize-1)/2-1,-(winsize-1)/2-1,-1):
                            cluster.append(remres1+rest)
                            cluster.append(remres2+rest)
                        coordinates_cluster = list(coordinates1)
                        for rest in range((winsize-1)/2-1,-(winsize-1)/2-1,-1):
                            del coordinates_cluster[remres2+rest]
                        for rest in range((winsize-1)/2-1,-(winsize-1)/2-1,-1):
                            del coordinates_cluster[remres1+rest]
##############################################################################

## loop over clustersizes
##############################################################################
##                clustersize_max = 0 ## make it a user defined variable and use approx 25 as the default
##                residue_max = len(matrix_distance_residue_intra)
##                if clustersize_max > 0.5*residue_max:
##                    clustersize_max = int(0.5*residue_max) ## (watch out cluster is not bigger than protein somehow...)
##                clusterplotdata = []
##                for mode_clustersize0 in range(6,12):
##                    clusterplotdata.append([])
##                    for mode in range(6,12):
##                        clusterplotdata[mode_clustersize0-6].append(
##                            Numeric.zeros((clustersize_max,residue_max),typecode='d')
##                            )
##                ## loop over number of residues in cluster
##                for clustersize in range(clustersize_max+1):
##
##                    if clustersize != 0:
##                        if os.path.isfile('cluster_%s_6_6.txt' %clustersize):
##                            continue
##                        for mode_clustersize0 in range(6,12):
##                            for mode in range(6,12):
##                                if not os.path.isfile('cluster_%s_%s_%s.txt' %(clustersize, mode_clustersize0, mode)):
##                                    fd = open('cluster_%s_%s_%s.txt' %(clustersize, mode_clustersize0, mode), 'w')
##                                    fd.close()
##                                else:
##                                    fd = open('cluster_%s_%s_%s.txt' %(clustersize, mode_clustersize0, mode), 'r')
##                                    lines = fd.readlines()
##                                    fd.close()
##                                    for iline in range(len(lines)):
##                                        clusterplotdata[mode_clustersize0-6][mode-6][clustersize-1][iline] = float(lines[iline])
#### cummulated scores for alle modes (testing) remove...
####                                        if mode_clustersize0 == mode and mode != 6:
####                                            clusterplotdata[0][0][clustersize-1][iline] += float(lines[iline])
##
##                    ## loop over residues (works only if only CA used...)
##                    clusters = []
##
##                    for residue_cluster_central in range(residue_max):
##
##                        ## 1) continue and do not loop over residues if cluster size is 0
##                        if clustersize == 0 and residue_cluster_central > 0:
##                            continue
##                        else:
##                            cluster, clusters, coordinates_cluster = self.identify_cluster(clustersize, residue_cluster_central, clusters, coordinates1, matrix_distance_residue_intra)
##
##                        ## 2) continue and do not calculate clusterplotdata if already calculated for an identical cluster
##                        if clustersize != 0 and cluster in clusters:
##                            for mode_clustersize0 in range(6,12):
##                                for mode in range(6,12):
##                                    clusterplotdata[mode_clustersize0-6][mode-6][clustersize-1][residue_cluster_central] = clusterplotdata[mode_clustersize0-6][mode-6][clustersize-1][clusters.index(cluster)]
##                            continue
##                        else:
##                            clusters.append(cluster)
##
##############################################################################

                        ## 3) clustersize == 0
                        if not vars().has_key('eigenvectors_clustersize0'):

                            matrix_hessian = self.hessian_calculation(coordinates1, float(cutoff_distance)) ## calculate with coordinates2 as well... and compare results to switching pdb1 and pdb2
                            eigenvectors, eigenvalues = self.function_eigenvector_calculation(matrix_hessian, amplitude_average, jobid)

                            eigenvectors_clustersize0 = eigenvectors[6:12]

                            lengths = self.morph(eigenvectors, frames, biomolecule, pdb1ATOM_all, atoms_hessian, amplitude_average, [], matrix_hessian, jobid, biolunitchains, cutoff_distance)

                            for mode in range(6,12):

##                                ## plot residue fluctuation (atom,RMSF) for each mode
##                                self.dislin_graf(
##                                    range(1,len(lengths[mode-6])+1), lengths[mode-6],
##                                    jobid+'_'+biolunitchains+'_'+str(cutoff_distance)+'_2D_atom_displacement'+str(mode+1)+'.png',
##                                    title1 = 'job: '+jobid,
##                                    title2 = 'chains: '+biolunitchains,
##                                    title3 = 'cutoff: '+str(cutoff_distance),
##                                    title4 = 'mode: '+str(mode+1),
##                                    xtitle = 'atom/residue', ytitle='RMSF'+str(mode+1),
##                                    YE = 1.1*max(lengths[mode-6]), YSTP = max(lengths[mode-6])/10, XE = len(lengths[mode-6])+1,
##                                    )

                                ccdata = Numeric.zeros((len(eigenvectors[mode])/3,len(eigenvectors[mode])/3), typecode='d')
                                for ccrow in range(len(ccdata)):
                                    for cccol in range(len(ccdata)):
                                        ccdata[ccrow][cccol] = self.cosangle(eigenvectors[mode][3*ccrow:3*ccrow+3], eigenvectors[mode][3*cccol:3*cccol+3])

                                filename = jobid+'_'+biolunitchains+'_'+str(cutoff_distance)+'_crosscorrelation'+str(mode+1)+'.png'

                                ## plot cross-correlation map (residue,residue,cross-correlation) for each mode
##                                self.dislin_graf3(
##                                    'residue', 'residue', 'cross-correlation',
##                                    filename,
##                                    ccdata,
##                                    title1 = 'job: '+jobid,
##                                    title2 = 'chains: '+biolunitchains,
##                                    title3 = 'cutoff: '+str(cutoff_distance),
##                                    title4 = 'mode: '+str(mode+1),
##                                    )

##                                ## add topology to margin of cross-correlations maps
##                                self.topology(filename, chain, biomolecule, helices1, strands1, ccdata)

                                ## end of loop over modes
                            
                            results.append(
                                {
                                    'number of non-zero eigenvalues': len(eigenvectors)-6,
                                    'depdate1': HEADERdepDate1,
                                    'cutoff': cutoff_distance,
                                    'res1': REMARK2resolution1,
                                    'biomolecules': biomolecule, 'chains1': pdb1ATOM_all.keys(),
                                    }
                                )

                            ## end of if clustersize == 0

                        ## 4) clustersize != 0
                        if clustersize != 0:

                            matrix_hessian = self.hessian_calculation(coordinates_cluster, float(cutoff_distance)) ## calculate with coordinates2 as well... and compare results to switching pdb1 and pdb2
                            eigenvectors, eigenvalues = self.function_eigenvector_calculation(matrix_hessian, amplitude_average, jobid)

                            overlaps, max_overlaps, modes_of_max_overlap, eigenvalues_of_max_overlap = self.overlap_calculation(eigenvectors, eigenvalues, eigenvectors_clustersize0, cluster)

                            for mode_clustersize0 in range(6,12):
                                for mode in range(6,12):
                                    ## use this line when variable clustersize and 1 cluster
##                                    clusterplotdata[mode_clustersize0-6][mode-6][clustersize-1][residue_cluster_central] = overlaps[mode_clustersize0-6][mode-6]
                                    ## use this line when fixed clustersize/windowsize and 2 clusters/windows
                                    clusterplotdata[mode_clustersize0-6][mode-6][remres1][remres2] = overlaps[mode_clustersize0-6][mode-6]

                        print 'end of loop over central cluster residue %s' %residue_cluster_central

                    ## write clusterplotdata to txt file in case of crash during loop over residues/clusters
                    if clustersize != 0:
                        for mode_clustersize0 in range(6,12):
                            for mode in range(6,12):
                                ## use this line when variable clustersize and 1 cluster
##                                fd = open('cluster_%s_%s_%s.txt' %(clustersize, mode_clustersize0, mode), 'a')
##                                for residue_cluster_central in range(len(clusterplotdata[mode_clustersize0-6][mode-6][clustersize-1])):
##                                    fd.write(str(clusterplotdata[mode_clustersize0-6][mode-6][clustersize-1][residue_cluster_central])+'\n')
                                ## use this line when fixed clustersize/windowsize and 2 clusters/windows
                                fd = open('cluster_%s_%s_%s.txt' %(remres1, mode_clustersize0, mode), 'a')
                                for remres2data in range(len(clusterplotdata[mode_clustersize0-6][mode-6][remres1])):
                                    fd.write(str(clusterplotdata[mode_clustersize0-6][mode-6][remres1][remres2data])+'\n')
                                fd.close()

                    print 'end of loop over cluster size %s' %clustersize

                ## plot clusterplotdata
                if clustersize_max > 0:
                    for mode_clustersize0 in range(6,12):
                        for mode_perturbed in range(6,12):
##remove these lines after testing
##                            print clusterplotdata[mode_clustersize0-6][mode_perturbed-6][40]
##                            stop
                            self.dislin_graf3(xtitle='number of residues in removed cluster',
                                              ytitle='central residue in removed cluster',
                                              ztitle='overlapx100',
                                              filename=jobid+'_'+biolunitchains+'_'+str(cutoff_distance)+'_'+str(mode_clustersize0)+'_'+str(mode_perturbed)+'_cluster.png',
                                              data=clusterplotdata[mode_clustersize0-6][mode_perturbed-6],
                                              title1 = 'job: '+jobid,
                                              title2 = 'chains: '+biolunitchains,
                                              title3 = 'cutoff: '+str(cutoff_distance),
                                              title4 = 'modes: '+str(mode_clustersize0+1)+', '+str(mode_perturbed+1),
                                              ZOR = 0.,
                                              )

                print 'end of loop over cutoff distance %s' %cutoff_distance


            ## end of loop over biomolecules

        print 'done'
        return results

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

    results = []

    instance_vibration = vibration()
    for pdbs in [
        ['8ACN', [], '1ATP', ['E']]
##pka 1atp:e, 1bkx:a
##chymotrypsin inhibitor 2 (1ypa, 1ypb, 1ypc; 3ci2; 2ci2:i)
#### no problems
##            ['1JEJ', [], '1JG6', [], 'D-h-2'], ['1DKX', [], '1DKY', [], 'D-h-2'], ['1CBU', [], '1C9K', [], 'F-s-2'], ['1CNP', [], '1A03', [], 'F-h-2'], ['1CEW', [], '1A67', [], 'F-h-2'], ['3ENL', [], '7ENL', [], 'F-h-2'], ['2HMY', [], '3MHT', [], 'F-h-2'], ['6LDH', [], '1LDM', [], 'F-h-2'], ['1F3Y', [], '1JKN', [], 'F-h-2'], ['1K9P', [], '1K9K', [], 'F-s-2'], ['1THV', [], '1THI', [], 'F-?-2'], ['1PRV', [], '1PRU', [], 'F-?-2'], ['2PHY', [], '3PYP', [], 'F-?-2'], ['1NBV', [], '1CBV', [], 'S-n-2'], ['1Q12', [], '1Q1B', [], 'S-n-2'], ['1SU4', [], '1T5S', [], 'F-s-2'], ['1G7S', [], '1G7T', [], 'D-?-2'], ['1B7T', [], '1DFK', [], 'D-?-2'], ['1VKX', [], '1RAM', [], 'D-?-2'], ['1D6M', [], '1I7D', [], 'D-?-2'], ['1MAR', [], '2ACQ', [], 'D-n-2'], ['1NCX', [], '1NCY', [], 'D-h-2'], ['1TTP', [], '1TTQ', [], 'D-h-2'], ['1SSP', [], '1AKZ', [], 'D-h-2'], ['1K92', [], '1KP2', [], 'D-?-2'], ['1QF5', [], '1HOO', [], 'D-?-2'], ['1GRN', [], '1CF4', [], 'D-?-2'], ['1IKU', [], '1JSA', [], 'D-h-2'], ['1KS2', [], '1KS2', [], 'D-h-2'], ['1QLN', [], '1MSW', [], 'D-h-2'], ['1OMP', [], '3MBP', [], 'D-h-2'], ['1LUA', [], '1LU9', [], 'D-h-2'], ['1RKM', [], '2RKM', [], 'D-h-2'], ['1JLU', [], '1CMK', [], 'D-h-2'], ['2LAO', [], '1LAF', [], 'D-h-2'], ['1G0X', [], '1P7Q', [], 'D-h-2'], ['1IBN', [], '1IBO', [], 'D-h-2'], ['1HKA', [], '1Q0N', [], 'D-h-2'], ['1OXS', [], '1OXU', [], 'D-h-2'], ['1FTO', [], '1FTM', [], 'D-h-2'], ['1N0V', [], '1N0U', [], 'D-h-2'], ['2NAC', [], '2NAD', [], 'D-h-2'], ['1D9V', [], '1MRP', [], 'D-h-2'], ['1E8B', [], '1E88', [], 'D-h-2'], ['1JBV', [], '1JBW', [], 'D-h-2'], ['2EIA', [], '1EIA', [], 'D-s-2'], ['4APE', [], '5ER2', [], 'D-s-2'], ['3EZA', [], '2EZA', [], 'D-s-2'], ['1G59', [], '1GLN', [], 'D-s-2'], ['1GD1', [], '2GD1', [], 'D-s-2'], ['1I6I', [], '1I5S', [], 'D-s-2'], ['1H9K', [], '1H9M', [], 'D-s-2'], ['1JYS', [], '1NC3', [], 'D-s-2'], ['1PVU', [], '1PVI', [], 'D-s-2'], ['1UJ1', [], '1UK2', [], 'D-s-2'], ['1FFH', [], '1NG1', [], 'D-s-2'], ['1WRP', [], '3WRP', [], 'D-s-2'], ['1EVL', [], '1EVK', [], 'D-s-2'], ['1EA5', [], '1AMN', [], 'D-h-2'], ['1AKE', [], '1ANK', [], 'D-h-2'], ['1I2D', [], '1M8P', [], 'D-h-2'], ['1L5B', [], '1L5E', [], 'D-h-2'], ['3DAP', [], '1DAP', [], 'D-h-2'], ['1COY', [], '3COX', [], 'D-s-2'], ['1NJG', [], '1NJF', [], 'D-s-2'], ['1LB4', [], '1LB5', [], 'F-?-2'], ['1FOX', [], '2FOW', [], 'F-n-2'], ['1LCC', [], '1LQC', [], 'F-n-2'], ['1J74', [], '1J7D', [], 'F-n-2'], ['2KTQ', [], '3KTQ', [], 'F-h-2'], ['1TGL', [], '4TGL', [], 'F-h-2'], ['1JMW', [], '1WAS', [], 'F-?-2'], ['3CHY', [], '1CHN', [], 'F-?-2'], ['3ERD', [], '1PCG', [], 'F-?-2'], ['1GU0', [], '1GU1', [], 'D-h-2'], 
#### sequence alignment returns 1 or more nonidentical residues (check seq identity and count nonidentical res)
##            ['1BYU', [], '1RRP', [], 'F-h-2'], ['1STO', [], '1ORO', [], 'F-h-2'], ['4MDH', [], '1BMD', [], 'F-h-2'], ['1CRL', [], '1THG', [], 'F-h-2'], ['1AXN', [], '2RAN', [], 'F-h-2'], ['1CQR', [], '1CAQ', [], 'F-?-2'], ['6TIM', [], '1TRE', [], 'F-h-2'], ['1PJR', [], '3PJR', [], 'C----'], ['1VDE', [], '1LWT', [], 'C----'], ['1IH7', [], '1IG9', [], 'C----'], ['1BGW', [], '1BJT', [], 'C----'], ['1JMJ', [], '1JMO', [], 'S-a-2'], ['1BUY', [], '1EER', [], 'S-a-2'], ['1FBT', [], '1TIP', [], 'S-a-2'], ['1LKF', [], '7AHL', [], 'D-f-2'], ['1DUJ', [], '1KLQ', [], 'D-f-2'], ['7API', [], '1PSI', [], 'D-f-2'], ['1GPW', [], '1THF', [], 'D-?-2'], ['1CJW', [], '1RAM', [], 'D-?-2'], ['1TBA', [], '1YTB', [], 'D-?-2'], ['1MCP', [], '1NCA', [], 'D-n-2'], ['1QAW', [], '1WAP', [], 'D-h-1'], ['1ETU', [], '1EFT', [], 'D-?-2'], ['1FGU', [], '1JMC', [], 'D-h-2'], ['1K20', [], '1K23', [], 'D-h-2'], ['1GTR', [], '1NYL', [], 'D-h-2'], ['1EJD', [], '1A2N', [], 'D-h-2'], ['1QUK', [], '1OIB', [], 'D-h-2'], ['13PK', [], '1PHP', [], 'D-h-2'], ['1IPD', [], '1OSJ', [], 'D-h-2'], ['1LFG', [], '1LFH', [], 'D-h-2'], ['8OHM', [], '1CU1', [], 'D-h-2'], ['1K89', [], '1HRD', [], 'D-h-2'], ['1DPP', [], '1DPE', [], 'D-h-2'], ['2EFG', [], '1FNM', [], 'D-h-2'], ['1EPS', [], '1G6S', [], 'D-h-2'], ['1ERK', [], '2ERK', [], 'D-h-2'], ['1GGG', [], '1WDN', [], 'D-h-2'], ['1GDH', [], '1PSD', [], 'D-h-2'], ['1B1A', [], '1BE1', [], 'D-s-2'], ['1N8Y', [], '1N8Z', [], 'D-s-2'], ['2YHX', [], '1HKG', [], 'D-s-2'], ['1BA3', [], '1LCI', [], 'D-s-2'], ['1DV7', [], '1DVJ', [], 'D-s-2'], ['1PAH', [], '1DVJ', [], 'D-s-2'], ['2POL', [], '1JQJ', [], 'D-s-2'], ['1BNC', [], '1DV2', [], 'D-h-2'], ['1HNF', [], '1HNG', [], 'D-h-2'], ['1GTM', [], '1HRD', [], 'D-h-2'], ['1CTR', [], '1CLL', [], 'D-h-2'], ['1AVK', [], '1A2V', [], 'D-s-2'], ['4PEP', [], '1PSN', [], 'D-s-2'], ['1BU7', [], '1JPZ', [], 'D-s-2'], ['1PIN', [], '1F8A', [], 'F-n-2'], ['1SER', [], '1SES', [], 'F-h-2'], ['1D5W', [], '1DCM', [], 'F-?-2'], ['2GLS', [], '1LGR', [], 'F-?-2'], ['1AB3', [], '1A32', [], 'F-?-2'], ['1BJY', [], '1BJZ', [], 'F-s-2'],
#### errors in pdbs
##            ['1CIP', [], '1GP2', [], 'D-f-2'], ['1G6O', [], '1NLZ', [], 'S-a-2'], ['1II0', [], '1IHU', [], 'S-n-2'], ['8ATC', [], '5AT1', [], 'S-a-2'],
#### DNA only
##            ['1HMH', [], '1HMH', [], 'N-R-2'], ['1GID', [], '1GID', [], 'N-R-2'], ['2TRA', [], '3TRA', [], 'N-R-2'],
#### missing information in pdbs about biological unit
##            ['1DQZ', [], '1DQY', [], 'F-h-2'], ['1BEB', [], '1B0O', [], 'F-h-2'], ['3TMS', [], '2TSC', [], 'F-s-2'], ['1I69', [], '1I6A', [], 'F-n-2'],
##            ['4CRX', [], '1CRX', [], 'D-h-2'], ['2CBL', [], '1B47', [], 'D-h-2'], ['1CKM', [], '1CKM', [], 'D-h-2'], ['1TDE', [], '1F6M', [], 'D-h-2'], ['2HMI', [], '3HVT', [], 'D-f-2'], ['4HHB', [], '2HCO', [], 'S-a-2'], ['1PFK', [], '2PFK', [], 'S-a-2'], ['1DKR', [], '1DKU', [], 'S-a-2'], ['9GPB', [], '1GPB', [], 'S-a-2'], ['1E0S', [], '1HFV', [], 'F-s-2'], ['1BAM', [], '1BHM', [], 'S-n-2'], ['1BRD', [], '2BRD', [], 'F-s-2'], ['1CLB', [], '4ICB', [], 'F-s-2'], ['4DFR', [], '5DFR', [], 'F-s-2'], ['1HDN', [], '1PFH', [], 'F-s-2'], ['1A8V', [], '2A8V', [], 'F-s-2'], ['1RCL', [], '1RCK', [], 'F-s-2'], ['5CRO', [], '6CRO', [], 'F-h-2'], ['1ECB', [], '1ECC', [], 'F-h-2'], ['4HVP', [], '3HVP', [], 'F-h-2'], ['3ICD', [], '1AI2', [], 'F-h-2'], ['1IK9', [], '1FU1', [], 'S-n-2'], ['4Q21', [], '6Q21', [], 'F-h-2'], ['2C4Q', [], '1DZS', [], 'F-n-2'], ['9AAT', [], '1AMA', [], 'D-s-2'], ['8ADH', [], '6ADH', [], 'D-s-2'], ['1J7N', [], '1JKY', [], 'D-s-2'], ['4CTS', [], '1CTS', [], 'D-s-2'], ['1G51', [], '1EFW', [], 'D-s-2'], ['1DDT', [], '1MDT', [], 'D-h-2'], ['1EX6', [], '1EX7', [], 'D-h-2'], ['1EX7', [], '1EX6', [], 'D-h-2'], ['1AA7', [], '1EA3', [], 'D-h-2'], ['1L96', [], '1L97', [], 'D-h-2'], ['1QAI', [], '1MML', [], 'D-h-2'], ['1BPD', [], '2BPG', [], 'D-h-2'], ['1URP', [], '2DRI', [], 'D-h-2'], ['1BP5', [], '1A8E', [], 'D-h-2'], ['1K93', [], '1K8T', [], 'D-?-2'], ['1GDT', [], '2RSL', [], 'D-?-2'], ['1IDG', [], '1IDI', [], 'D-f-2'], ['1AON', [], '1OEL', [], 'D-h-2'], ['1AON', [], '1OEL', [], 'C----'], ['1AON', [], '1EGS', [], 'F-?-2'],
##            ['2LZM',[' '],'150L',['A'],'--h-2'],
##            ['2LZM',[' '],'150L',['B'],'--h-2'],
##            ['2LZM',[' '],'150L',['C'],'--h-2'],
##            ['2LZM',[' '],'150L',['D'],'--h-2'],
##            ['150L',['A'],'2LZM',[' '],'--h-2'],
##            ['150L',['B'],'2LZM',[' '],'--h-2'],
##            ['150L',['C'],'2LZM',[' '],'--h-2'],
##            ['150L',['D'],'2LZM',[' '],'--h-2'],
        ]:

        pdblines1 = instance_vibration.pdb_import(pdbs[0])
        pdblines2 = instance_vibration.pdb_import(pdbs[2])

        import time, os
        timetuple = time.gmtime(time.time())
        job = str(timetuple[0])+str(timetuple[1]).zfill(2)+str(timetuple[2]).zfill(2)+str(timetuple[3]).zfill(2)+str(timetuple[4]).zfill(2)+str(timetuple[5]).zfill(2)

        results = instance_vibration.main(
            pdblines1, pdbs[1], None,
            ['CA'], float(1), [6,8,10,12,14,16,18,20,22], 3, 'monomeric', job, 50,
            pdblines2, pdbs[3], None, 
            )

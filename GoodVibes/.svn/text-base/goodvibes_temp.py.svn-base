##!/bin/env /software/bin/python2.3
##
##$Id: goodvibes2.py 118 2007-10-04 10:22:23Z tc $
##
##Tommy Carstensen, University College Dublin, 2007
## make sure that missng atoms and residues don't go ilent by
## hvad skal der ske med atomer med altloc?
## use info from MODRES lines if modified residue

## this version does not include helix and proline rigidity

class vibration:


    def morph(self, eigenvectors_all_modes, nframes, chains, d_coordinates, atoms_hessian, matrix_hessian, job,):

        '''This function puts in frames between the maximum amplitudes of a
movement given by eigenvectors. The values of the diagonal elements of the
hessian matrix are used for B factors to make coloring during simulation
possible.'''

##        print 'visualize the two extreme projections along a trajectory and interpolate n frames between them'

        for mode in range(6,12):

            eigenvectors = eigenvectors_all_modes[mode]

##            eigenvectors = length_adjustment(mode, eigenvectors)

            output_vmd = ['REMARK color by connectivity (b-factor) or squared displacement (temperature factor)\n']
            for frame in range(nframes):
                output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
                i = 0
                for chain in chains:
                    res_nos = d_coordinates['chains'][chain]['residues'].keys()
                    res_nos.sort()
                    for res_no in res_nos:
                        for iCode in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys():
                            for altloc in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].keys():
                                for atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'].keys():
                                    if 'hessian' in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'][atom_name].keys():
                                        coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'][atom_name]['coordinate']
                                        res_name = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['res_name']
                                        x1 = 9*eigenvectors[3*(i/winsize)+0]
                                        y1 = 9*eigenvectors[3*(i/winsize)+1]
                                        z1 = 9*eigenvectors[3*(i/winsize)+2]
                                        sqlength = x1**2+y1**2+z1**2
                                        x2 = coordinate[0]+(1-2*float(frame)/float(nframes))*x1
                                        y2 = coordinate[1]+(1-2*float(frame)/float(nframes))*y1
                                        z2 = coordinate[2]+(1-2*float(frame)/float(nframes))*z1
                                    
                                        output_vmd += [
                                            'ATOM        %4s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(
                                                atom_name, altloc, res_name, chain, int(res_no), iCode, x2,y2,z2, sqlength, sqlength
                                                )
                                            ]
                                        i += 1

                output_vmd.append('TER\nENDMDL\n')
            fd = open(job+'_mode'+str(mode+1).zfill(2)+'.pdb', 'w') ## implicit path
            fd.writelines(output_vmd)
            fd.close()

        return


    def main(
        self, jobid, lines, atoms_hessian = ['CA'], frames = 50,
        cutoff_distance = 10.,
        path_html = None, path_python = None, verbose = False, paralleldir = '',
        biomolecule = None, chains = [], model = None,
        winsize = 1, stepsize = 1,
        ):

        '''
        Use first model if no model specified by user.
        chain(s): Y, biomolecule: Y; parse chains specified by user and apply transformation
        chain(s): Y, biomolecule: N; parse chains specified by user but don't apply transformation
        chain(s): N, biomolecule: Y; parse chains of biomolecule and apply transformation
        chain(s): N, biomolecule: N; parse chains of first biomolecule and apply transformation
        '''

        import os, Numeric

        if stepsize not in [1,winsize]:
            raise 'stepsize must be 1 or equal to winsize'

        results = []

        ## parse pdb
        (
            d_REMARK350,
            d_primary, ## i.e. SEQRES, MODRES
            d_secondary, ## i.e. HELIX, SHEET
            d_coordinates, ## i.e. ATOM, HETATM, TER, MODEL, ENDMDL
            ) = self.parse_pdb(lines, chains, model)

        ## assume multimeric biological unit if chains not specified by user
        if chains == []:
            chains = d_coordinates['chains'].keys()
            chains.sort()

        ##
        ## d_coordinates to l_coordinates
        ##
        l_coordinates = []
        for chain in chains:
            ## assuming sequential numbering of residues
            res_nos = d_coordinates['chains'][chain]['residues'].keys()
            res_nos.sort()
            for res_no in res_nos:
                for iCode in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys():
                    altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].keys())
                    for atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'].keys():
                        if atom_name in atoms_hessian:
                            coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'][atom_name]['coordinate']
                            d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'][atom_name]['hessian'] = True
                            l_coordinates += [coordinate]

        ## cluster coordinates
        l_coordinates = [sum(l_coordinates[i:i+winsize])/winsize for i in range(0,len(l_coordinates),winsize)]
        N = len(l_coordinates)

        ## calculate intra-residue distances
        matrix_distance_residue_intra = self.distance_residue_intra(l_coordinates)

        ## for the non-disrupted structure: calculate and visualize eigenvectors
        matrix_hessian = self.hessian_calculation(l_coordinates, float(cutoff_distance), verbose)
        eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_combined_nonperturbed = self.eigenv_calccomb(matrix_hessian, jobid, verbose)
        import math
        for mode in range(6,12):
            for i in range(20):
                E = 0
                f = 0.00001*(2**i)
                if len(matrix_hessian) != len(eigenvectors_nonperturbed[mode]):
                    stop
                if len(matrix_hessian) != len(eigenvectors_nonperturbed):
                    stop
                for j in range(len(matrix_hessian)):
                    ## E = .5kx**2
                    E += .5*matrix_hessian[j][j]*(f*eigenvectors_nonperturbed[mode][j])**2
##                    print matrix_hessian[j][j], eigenvectors_nonperturbed[mode][j], eigenvectors_nonperturbed[mode][j]**2, .5*matrix_hessian[j][j]*(eigenvectors_nonperturbed[mode][j])**2
##                print f,E

            V = 0
            E = 0
            for j in range(0,len(matrix_hessian),3):
                vx = eigenvectors_nonperturbed[mode][j+0]
                vy = eigenvectors_nonperturbed[mode][j+1]
                vz = eigenvectors_nonperturbed[mode][j+2]
                v = math.sqrt(vx**2+vy**2+vz**2)
                kx = matrix_hessian[j+0][j+0]
                ky = matrix_hessian[j+1][j+1]
                kz = matrix_hessian[j+2][j+2]
                k = math.sqrt(kx**2+ky**2+kz**2)
                E = .5*kx*(vx**2)+.5*ky*(vy**2)+.5*kz*(vz**2)
                print k,v,E
                V += vx+vy+vz
##            print V
##                stop
        stop
        self.morph(eigenvectors_nonperturbed, frames, chains, d_coordinates, atoms_hessian, matrix_hessian, jobid,)

        ## do plots for individual modes
        self.plot(
            jobid,
            cutoff_distance,
            eigenvectors_nonperturbed,
            chains,d_secondary,
            fileprefix = 'individual',
            )
        ## do plots for combined modes
        self.plot(
            jobid,
            cutoff_distance,
            eigenvectors_combined_nonperturbed,
            chains,
            d_secondary,
            fileprefix = 'combined',
            )

        results.append(
            {
                'number of non-zero eigenvalues': len(eigenvalues_nonperturbed)-6,
                'cutoff': cutoff_distance,
                'biomolecules': biomolecule, 'chains1': d_coordinates['chains'].keys(),
                }
            )

################################################################################
##                     plotdata initialization section                        ##
################################################################################

        ## set data lists and append matrices to be plotted for each combination of modes 6-12 before initiating loops over the two axes of the plot
        resres_overlaps_single = []
        resres_overlaps_combined = []
        resres_eigenvalue_of_max_overlap = []
        resres_mode_of_max_overlap = []
        resres_eigenvalues_perturbed = []
        for mode in range(3*N):
            if mode not in range(6,18)+[3*N-1]:
                resres_eigenvalue_of_max_overlap.append('')
                resres_mode_of_max_overlap.append('')
                resres_eigenvalues_perturbed.append('')
                resres_overlaps_single.append('')
                resres_overlaps_combined.append('')
            else:
                resres_eigenvalue_of_max_overlap.append(Numeric.zeros((N,N),typecode='d'))
                resres_mode_of_max_overlap.append(Numeric.zeros((N,N),typecode='d'))
                resres_eigenvalues_perturbed.append(Numeric.zeros((N,N),typecode='d'))
                resres_overlaps_single.append(Numeric.zeros((N,N),typecode='d'))
                resres_overlaps_combined.append(Numeric.zeros((N,N),typecode='d'))
        datadic = {
            'emo': {'data': resres_eigenvalue_of_max_overlap, 'name': 'change in eigenvalue of max overlap', 'diagonal': 0, 'zrange': ['','']},
            'mmo': {'data': resres_mode_of_max_overlap, 'name': 'perturbed modes of max overlap', 'diagonal': 'mode', 'zrange': [7,12]},
            'overlaps_single': {'data': resres_overlaps_single, 'name': 'overlaps of single modes', 'diagonal': 'max', 'zrange': ['','']},
            'overlaps_combined': {'data': resres_overlaps_combined, 'name': 'overlaps of combined modes', 'diagonal': 'max', 'zrange': ['','']},
            'eigenvalues_perturbed': {'data': resres_eigenvalues_perturbed, 'name': 'eigenvalues of perturbed modes', 'diagonal': 'max', 'zrange': ['','']},
            }

################################################################################
##                        plotdata calculation section                        ##
################################################################################

        ## loop over remres1
        for remres1 in range(N):

            ## continue the remres1 loop if another processor is running calculations for that value of remres1
            filename = '%s/%ssize%s_residue%s.txt' %(paralleldir, jobid, winsize, remres1+1)
            if os.path.isfile(filename):
## REWRITE THIS SECTION
                ## read datafile to data variables before you continue
                fd = open(filename, 'r')
                lines = fd.readlines()
                fd.close()
                for i in range(len(lines)):
                    line = lines[i]
                    if line[:4] == 'data':
                        line_number_data = i+2
                        key = line[5:-1]
                    elif line[:4] == 'mode':
                        line_number_data = i+1
                        mode = int(line[5:-1])-1
                    elif i >= line_number_data:
                        remres2 = i - line_number_data
                        if remres2 < remres1:
                            continue
                        datadic[key]['data'][mode][remres1][remres2] = float(line)
                        datadic[key]['data'][mode][remres2][remres1] = float(line)
                    else:
                        print line, i, line_number_data, remres1, remres2
                        stop

                continue
            
            ## create data files (necessary to create data files if multithreading/parallel processing dependent on existence of files to determine if calculation on certain residues has been performed)
            else:
                fd = open(filename, 'w')
                fd.close()

            ## loop over remres2
            for remres2 in range(N):

                if not path_html:
                    print 'residue combination %s,%s of %s,%s' %(remres1+1,remres2+1,N,N)

                ## continue if windows overlap with each other
                if winsize*abs(remres1 - remres2) < stepsize:
                    continue
                ## continue if windows overlap with borders of the matrix
                if remres1 < (winsize-1)/2 or remres1 > (N-1)-(winsize-1)/2 or remres2 < (winsize-1)/2 or remres2 > (N-1)-(winsize-1)/2:
                    continue
                ## continue if overlap values have already been calculated for 
                if remres2 < remres1:
                    continue

                coordinates_cluster, cluster = self.residueclustering(
                    remres1, l_coordinates, jobid, cutoff_distance, winsize,
                    datadic, remres2, N, stepsize,
                    )

                ## calculate data
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
        ##        self.morph(eigenvectors, frames, biomolecule, d_coordinates, atoms_hessian, cluster, matrix_hessian, jobid+'-'+str(xvalue)+'-'+str(yvalue), cutoff_distance)
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
                    'mmo': perturbed_modes_of_max_overlap_single,
                    'overlaps_combined': overlaps_combined
                    }
                for key in datadic:
                    for mode in range(6,12):
                        datadic[key]['data'][mode][remres1][remres2] = datadic_loop[key][mode]
                        datadic[key]['data'][mode][remres2][remres1] = datadic_loop[key][mode]
                datadic['overlaps_combined']['data'][-1][remres1][remres2] = overlaps_combined[-1]
                datadic['overlaps_combined']['data'][-1][remres2][remres1] = overlaps_combined[-1]

                ## calculate remaining time based only on time spent on diagonalization (excluding finalizing steps with graphics generation)
                if path_html:
                    residuecount = len(coordinates_cluster)+len(cluster)
                    rem_col = (residuecount-winsize-2.*(winsize-1)/2.) - (xvalue+1) + (winsize-1)/2.
                    rem_matrices_col = (rem_col/2.)*(1+rem_col)
                    rem_matrices_row = residuecount - (yvalue+1) - (winsize-1)/2.
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

            ## write data to txt files in case of crash during loop over residues/clusters
            lines = []
            for key in datadic:
                lines.append('data '+key+'\n')
                data = datadic[key]['data']
                for mode in range(6,12):
                    lines.append('mode '+str(mode+1)+'\n')
                    for remres2 in range(len(data[mode][remres1])):
                        lines.append(str(data[mode][remres1][remres2])+'\n')
                if key == 'overlaps_combined':
                    mode = len(data)-1
                    lines.append('mode '+str(mode+1)+'\n')
                    for remres2 in range(len(data[mode][remres1])):
                        lines.append(str(data[mode][remres1][remres2])+'\n')
            fd = open(filename, 'w')
            fd.writelines(lines)
            fd.close()

################################################################################
##                                plot section                                ##
################################################################################
                        
        axistitle = 'central residue in removed cluster'
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
                for i in range((winsize-1)/2,N-(winsize-1)/2):
                    for j in range(winsize):
                        for k in range(winsize):
                            datadic[key]['data'][mode][i-(winsize-1)/2+j][i-(winsize-1)/2+k] = diagonal

        xtitle = axistitle
        ytitle = axistitle
        ztitle = '{/Symbol D}overlap'
        filename = jobid+'overlapscombined_mode3N'
        title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %('{/Symbol D}overlap', 'modes 7-3N', jobid, chains, cutoff_distance),

        ## overlaps_combined between nonperturbed and perturbed modes 7-3N
        self.gnuplot_splot(
            xtitle, ytitle, ztitle,
            filename,
            resres_overlaps_combined[-1],
            title,
            z1 = '', z2 = '',
            )

        ## add topology to margin of cross-correlations maps
        filename = jobid+'overlapscombined_mode3N.png'
        self.topology(
            filename,chains,d_secondary,N,
            )
        
        for mode in range(6,12):

            for key in datadic:

                ztitle = datadic[key]['name']
                filename = '%s%s_mode%s' %(jobid, key, mode+1)
                title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %(datadic[key]['name'], mode+1, jobid, chains, cutoff_distance),
                self.gnuplot_splot(
                    xtitle, ytitle, ztitle,
                    filename,datadic[key]['data'][mode],title,
                    z1 = datadic[key]['zrange'][0], z2 = datadic[key]['zrange'][1],
                    )

                ## add topology to margin of cross-correlations maps
                filename = '%s%s_mode%s.png' %(jobid, key, mode+1)
                self.topology(
                    filename,chains,d_secondary, N,
                    )

            ## perturbed eigenvalue of 6 perturbed modes
            data = resres_eigenvalues_perturbed

        return results


    def parse_pdb(self, lines, parse_chains, parse_model):

##        print 'parsing info from coordinate section about residues, atoms and coordinates'

        import Numeric, sets

        d_coordinates = {'chains':{}} ## ATOM, HETATM, MODEL
        d_REMARK350 = {} 
        d_secondary = {'HELIX':{},'SHEET':{},} ## HELIX, SHEET
        d_primary = {} ## MODRES
        model = None
        l_MODRES = []

        for i in range(len(lines)):

            line = lines[i]

            record = line[:6].strip()

            if record == 'REMARK':

                remark = int(line[6:10].strip())

                if remark == 350:

                    d_REMARK350 = self.parse_REMARK350(lines,i, d_REMARK350)
                                    
            elif record == 'MODRES':
                chain = line[16]
                res_no = int(line[18:22])
                res_name = line[24:27]
                l_MODRES += [res_name]
                continue

            elif record == 'HELIX':
                chain = line[19]
                if not d_secondary['HELIX'].has_key(chain):
                    d_secondary['HELIX'][chain] = []
                d_secondary['HELIX'][chain] += [[int(line[21:25]), int(line[33:37])]]
                continue

            elif record == 'SHEET':
                chain = line[21]
                if not d_secondary['SHEET'].has_key(chain):
                    d_secondary['SHEET'][chain] = []
                d_secondary['SHEET'][chain] += [[int(line[22:26]), int(line[33:37])]]
                continue

            elif record == 'MODEL':
                model = int(line[10:14])
                continue

            elif record == 'ATOM':
                if model != parse_model:
                    continue
                res_name = line[17:20].strip()
                d_coordinates = self.parse_ATOM(line, d_coordinates)
                continue

            elif record == 'HETATM':
                if model != parse_model:
                    continue
                chain = line[21]
                res_no = int(line[22:26].strip())
                res_name = line[17:20].strip()
                if not res_name in l_MODRES:
                    continue
                d_coordinates = self.parse_atom(line, d_coordinates)
                continue

        return d_REMARK350, d_primary, d_secondary, d_coordinates


    def parse_REMARK350(self, lines, i, d_REMARK350):

        import sets

        line = lines[i]

        if line[11:23] == 'BIOMOLECULE:':
            biomolecules = line[23:80].replace(' ','').split(',')

            chains = sets.Set()

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
                    break

                elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                    chains = sets.Set() ## e.g. 1upp.pdb (vs 8ruc.pdb)
                    line_chains = lines[j][41:80]
                    chains |= self.parse_REMARK350_chains(line_chains)

                elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
                    line_chains = lines[j][53:80]
                    chains |= self.parse_REMARK350_chains(line_chains)

                elif ',' in lines[j][11:80]:
                    if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
                        line_chains = lines[j][11:80]
                        chains |= self.parse_REMARK350_chains(line_chains)

                ## count and parse chain transformations
                ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj.pdb)
                elif lines[j][13:19] in ['BIOMT3','SMTRY3']:

                    matrixno = int(lines[j][19:24])
                    ## parse transformation matrix
                    matrixrow1 = lines[j-2][24:].split()
                    matrixrow2 = lines[j-1][24:].split()
                    matrixrow3 = lines[j-0][24:].split()
                    matrixrows = [matrixrow1,matrixrow2,matrixrow3,]
    ##                ## find out whether transformation matrix yields a transformation
    ##                transformation = False
    ##                for k in range(3):
    ##                    ## add a zero translation vector if a translation vector is not given
    ##                    if len(matrixrows[k]) == 3:
    ##                        matrixrows[k] += [0.]
    ##                    if float(matrixrows[k][k]) == 1. and float(matrixrows[k][3]) == 0.:
    ##                        continue
    ##                    else:
    ##                        transformation = True

                    ## append transformation matrix to dictionary
                    for biomolecule in biomolecules:

                        biomolecule = int(biomolecule)

                        ## biomolecule
                        if biomolecule not in d_REMARK350.keys():
                            d_REMARK350[biomolecule] = {}

                        ## biomolecule > matrices
                        if 'matrices' not in d_REMARK350[biomolecule].keys():
                            d_REMARK350[biomolecule]['matrices'] = {}
                        ## matrices > matrixno > matrix
                        d_REMARK350[biomolecule]['matrices'][matrixno] = matrixrows

                        ## biomolecule > chains
                        if 'chains' not in d_REMARK350[biomolecule].keys():
                            d_REMARK350[biomolecule]['chains'] = {}
                        for chain in chains:
                            ## chains > chain
                            if chain not in d_REMARK350[biomolecule]['chains'].keys():
                                d_REMARK350[biomolecule]['chains'][chain] = sets.Set()
                            d_REMARK350[biomolecule]['chains'][chain] |= sets.Set([matrixno])


        return d_REMARK350


    def parse_REMARK350_chains(self, line_chains):

        import sets

        ## if sentence necessary due to e.g. 1qgc
        ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: 1 2 3 4 5
        if ',' not in line_chains:
            chains = line_chains.split()
        else:
            ## remove 'AND' from the line of chains (e.g. problem with 1rhi)
            ## replace '.' in the line of chains (e.g. problem with 1rbo and 1qgc)
            chains = line_chains.replace('AND',',').replace('.',',').replace(' ',',').split(',')

        ## loop removal of blank chains necessary due to e.g. 2g8g
        ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, ,
        for x in range(100):
            if '' in chains:
                chains.remove('')
            else:
                break

        for j in range(len(chains)):
            chain = chains[j]
            if chain == 'NULL':
                chains[j] = ' '

        return sets.Set(chains)


    def loop_and_identify_biomolecules(self, i, lines):

        for j in range(i-1,-1,-1):

            if lines[j][:10] != 'REMARK 350':
                return False
            if lines[j][11:23] == 'BIOMOLECULE:':
                return True


    def parse_ATOM(self, line, d_coordinates):

        import Numeric

        atom_no = int(line[6:11])
        atom_name = line[12:16].strip()
        altloc = line[16]
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coordinate = Numeric.array([x, y, z])
        occupancy = float(line[56:60])

        if not d_coordinates['chains'].has_key(chain):
            d_coordinates['chains'][chain] = {'residues':{}}
        if not d_coordinates['chains'][chain]['residues'].has_key(res_no):
            d_coordinates['chains'][chain]['residues'][res_no] = {'iCodes':{}}
        if not d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].has_key(iCode):
            d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode] = {'altlocs':{}}
        if not d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].has_key(altloc):
            d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc] = {'atoms':{}}
        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'][atom_name] = {'coordinate':coordinate}
        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['res_name'] = res_name

        return d_coordinates


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

##        ## loop over modes
##        for i in range(7,len(eigenvalues)-1):
##
##            ## calculate length of mode i
##            leni = math.sqrt(sum(Numeric.array(eigenvectors[i])*Numeric.array(eigenvectors[i])))
##
##            ## scale length of mode i relative to length of mode 7
##            lenfactor = (len7/leni)/(eigenvalues[i]/eigenvalues[6])
##            for j in range(len(eigenvectors[i])):
##                eigenvectors[i][j] *= lenfactor

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


    def pdb_import(self, pdb):

        
        import os, urllib2

        ## local dir
        if os.path.isfile('%s.pdb' %(pdb)):
            print 'importing pdb %s from local dir' %pdb
            fd = open(pdb+'.pdb', 'r')
            lines = fd.readlines()
            fd.close()
        elif os.path.isfile('%s.pdb' %(pdb.lower())):
            print 'importing pdb %s from local dir' %pdb
            fd = open(pdb.lower()+'.pdb', 'r')
            lines = fd.readlines()
            fd.close()
        elif os.path.isfile('%s.pdb' %(pdb.upper())):
            print 'importing pdb %s from local dir' %pdb
            fd = open(pdb.upper()+'.pdb', 'r')
            lines = fd.readlines()
            fd.close()
        elif os.path.isfile('%s%s/pdb%s.ent' %(self.path_pdb, pdb.lower()[1:3], pdb.lower())):
            print 'importing pdb %s from local pdb repository' %pdb
            fd = open('%s%s/pdb%s.ent' %(self.path_pdb, pdb.lower()[1:3], pdb.lower()), 'r')
            lines = fd.readlines()
            fd.close()
        else:
            print 'importing pdb %s from www.pdb.org' %pdb
            url = urllib2.urlopen('http://www.pdb.org/pdb/files/%s.pdb' %(pdb))
            lines = url.readlines()

        return lines


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


    def topology(
        self, filename, chains, d_secondary, axlen
        ):

        import os

        bl = [390,833] ## bottom left coordinate
        tr = [1050,173] ## top right coordinate
        s = 14 ## space between plot and topology
        w = (tr[0]-bl[0])/axlen ## width of squares

        lines = 'convert '+filename
        for chain in chains:
            for secelm in d_secondary.keys():
                if secelm == 'SHEET':
                    lines += ' -fill "rgb(0,0,255)" '
                elif secelm == 'HELIX':
                    lines += ' -fill "rgb(255,0,0)" '
                else:
                    continue
                if chain in d_secondary[secelm]:
                    for secelmrange in d_secondary[secelm][chain]:
                        ## bottom
                        lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), bl[1]+s, int(bl[0]+w*secelmrange[1]), bl[1]+s)
                        ## top
                        lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), tr[1]-s, int(bl[0]+w*secelmrange[1]), tr[1]-s)
                        ## left
                        lines += ' -draw "line %s,%s %s,%s" ' %(bl[0]-s, int(bl[1]-w*(secelmrange[0]-1)), bl[0]-s, int(bl[1]-w*secelmrange[1]))
                        ## right
                        lines += ' -draw "line %s,%s %s,%s" ' %(tr[0]+s, int(bl[1]-w*(secelmrange[0]-1)), tr[0]+s, int(bl[1]-w*secelmrange[1]))
        if lines != 'convert %s' %(filename):
            lines += filename
            os.system(lines)
                
        return


    def plot(
        self, jobid, cutoff_distance, eigenvectors, chains, d_secondary,
        fileprefix = '',
        ):

        import math, Numeric

        for mode in range(6,12):

            ## calculate residue fluctuation (atom,RMSF) for each mode

            lengths = []
            for i in range(0,len(eigenvectors[mode]),3):
                length = math.sqrt(
                    eigenvectors[mode][i+0]**2+
                    eigenvectors[mode][i+1]**2+
                    eigenvectors[mode][i+2]**2
                    )
                lengths.append(length)

##            ## write lengths to file
##
##            lines = []
##            for i in range(len(lengths)):
##                lines.append(str(lengths[i])+'\n')
## the line below is for when working with nmr ensembles
##            fd = open(jobid+'_model'+str(model[0])+'_'+str(mode)+'_lengths.txt', 'w')
## the lien below is when not working with nmr ensembles
##            fd = open(jobid+'_'+str(mode)+'_lengths.txt', 'w')
##            fd.writelines(lines)
##            fd.close()

            ## plot residue fluctuation (atom,RMSF) for each mode

            if fileprefix == 'individual':
                title1 = 'mode %s' %(mode+1)
            elif fileprefix == 'combined':
                title1 = 'modes 7-%s' %(mode+1)

            self.gnuplot_plot(
                jobid, cutoff_distance, chains, title1,
                xtitle = 'residue', ytitle = 'RMSF', plotname = 'residue displacements', filename = '%sdisplacement_mode%s' %(fileprefix, mode+1),
                data = lengths,
                )

            ## prepare data for cross correlation maps
            ccdata = Numeric.zeros((len(eigenvectors[mode])/3,len(eigenvectors[mode])/3), typecode='d')
            for ccrow in range(len(ccdata)):
                for cccol in range(len(ccdata)):
                    ccdata[ccrow][cccol] = self.cosangle(eigenvectors[mode][3*ccrow:3*ccrow+3], eigenvectors[mode][3*cccol:3*cccol+3])

            ## plot cross-correlation map (residue,residue,cross-correlation) for each mode
            filename = fileprefix+'crosscorrelation_mode%s' %(mode+1)

            for i in range(len(ccdata)):
                ccdata[i][i] = 1
            xtitle = 'residue'
            ytitle = 'residue'
            ztitle = 'cross correlation'
            title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %('cross correlation map', mode+1, jobid, chains, cutoff_distance),

            self.gnuplot_splot(
                xtitle, ytitle, ztitle,
                filename,ccdata,title,
                z1 = -1, z2 = 1,
                )

            filename += '.png'
            ## add topology to margin of cross-correlations maps
            self.topology(
                filename, chains, d_secondary, len(ccdata),
                )

        return


    def gnuplot_splot(
        self, xtitle, ytitle, ztitle, filename, data,
        title,
        z1 = None, z2 = None,
        ):

        import os

        ## write gnuplot data to txt file
        gnuplot_splot_data = []
        for x in range(len(data)):
            for y in range(len(data[x])):
                if x >= y:
                    z = data[x][y]
                if y > x:
                    z = data[y][x]
                gnuplot_splot_data.append('%4i %4i %16.13f\n' %(x+1, y+1, z))
            gnuplot_splot_data.append('%4i %4i %16.13f\n\n' %(x+1, len(data[x])+1, data[0][0]))
        for i in range(len(data)+1):
            gnuplot_splot_data.append('%4i %4i %16.13f\n' %(len(data)+1, i+1, data[0][0]))
        fd = open('%s.gnuplotdata' %(filename), 'w')
        fd.writelines(gnuplot_splot_data)
        fd.close()
        ## write gnuplot settings to txt file
        lines = ['set size square\n'] ## scale square
        lines += [
            'set terminal postscript eps enhanced color "Helvetica" 48\n',
            'set output "%s.ps"\n' %(filename),
            'set size 4,4\n', ## scale 400%
            'set view map\n', ## change orientation of plot
            'set autoscale fix\n', ## scale axes
            'set style data pm3d\n', ## set by default?
            'set style function pm3d\n', ## set by default?
            'set encoding iso_8859_1\n', ## postscript encoding for special characters
            'set title "%s" ,4\n' %(title),
            'set xlabel "%s"\n' %(xtitle),
            'set ylabel "%s"\n' %(ytitle),
            'set cblabel "%s"\n' %(ztitle),
            ]
        if z1 != None and z2 != None:
            lines += [
                'set cbrange [%s:%s]\n' %(z1,z2), ## set colorbox range
                ]
        lines += [
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


    def residueclustering(
        self, remres1, coordinates1, jobid, cutoff_distance, winsize, datadic,
        remres2, N, stepsize,
        ):

        import math

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
            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
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

    import os, sys

    if '-pdb' not in sys.argv:
        raise 'use -pdb to select a pdb'
    pdb = sys.argv[sys.argv.index('-pdb')+1][:4]
    if '-chains' in sys.argv:
        chains = sys.argv[sys.argv.index('-chains')+1].split(',')
    else:
        chains = []
    if '-winsize' in sys.argv:
        winsize = int(sys.argv[sys.argv.index('-winsize')+1])
    else:
        winsize = 1
    if '-biomolecule' in sys.argv:
        biomolecule = int(sys.argv[sys.argv.index('-biomolecule')+1])
    else:
        biomolecule = None
    stepsize = winsize

    instance_vibration = vibration()
    lines = instance_vibration.pdb_import(pdb)

    results = instance_vibration.main(
        pdb,
        lines,
        chains = chains,
        winsize = winsize,
        stepsize = stepsize,
        biomolecule = biomolecule,
        verbose = True,
        paralleldir = os.getcwd(),
        )

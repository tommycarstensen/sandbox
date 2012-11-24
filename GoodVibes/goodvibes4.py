##!/bin/env /software/bin/python2.3
##
##$Id$
##
##Tommy Carstensen, University College Dublin, 2007

## version 1 systematically removes alpha carbon atoms
## version 2 build on version1
## version 2 allows choice of pdb and chain from terminal
## version 2 allows different winsizes
## version 3 allows for the inclusion of atoms other than alpha carbon atoms
## version 3 performs an alanine scan (removal of beta carbon atoms upon calc of distance)
## version 4 builds on version 3
## version 4 allows for ligand sidechain extension
## version "rigidity" builds on version 2
## version "large residue" is a clean up of version 3 (back to version 2)

class vibration:


    def morph(self, eigenvectors_all_modes, nframes, chains, d_coordinates, atoms_hessian, matrix_hessian, job, d_MODRES):

        '''This function puts in frames between the maximum amplitudes of a
movement given by eigenvectors. The values of the diagonal elements of the
hessian matrix are used for B factors to make coloring during simulation
possible.'''

##        print 'visualize the two extreme projections along a trajectory and interpolate n frames between them'

        for mode in range(6,12):

            eigenvectors = eigenvectors_all_modes[mode]

##            eigenvectors = length_adjustment(mode, eigenvectors)

            ## pdb header
            output_vmd = []

            ## loop over frames
            for frame in range(nframes):
                ## frame header
##                output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
                output_vmd.append('HEADER    frame t= %4.3f\nMODEL     %4i\n' %(frame,frame+1))
                ## loop over coordinates
                chains = d_coordinates['chains'].keys()
                chains.sort()
                for chain in chains:
                    ## assume sequential numbering of residues
                    res_nos = d_coordinates['chains'][chain]['residues'].keys()
                    res_nos.sort()
                    for res_no in res_nos:
                        ## assume sequential iCodes
                        iCodes = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys()
                        iCodes.sort()
                        for iCode in iCodes:
                            altloc_residue = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].keys())
                            res_name = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc_residue]['res_name']
                            atom_names = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys()
                            atom_names.sort()
                            for atom_name in atom_names:
                                altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
                                coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['coordinate']
                                atom_no = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['atom_no']

                                x = coordinate[0]
                                y = coordinate[1]
                                z = coordinate[2]

                                ## eigenvector?
                                if 'i' in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc].keys():
                                    i = int(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['i'])
                                    vx = 20*eigenvectors[3*i+0]
                                    vy = 20*eigenvectors[3*i+1]
                                    vz = 20*eigenvectors[3*i+2]
                                    sqlength = vx**2+vy**2+vz**2
                                    x += (1-2*float(frame)/float(nframes))*vx
                                    y += (1-2*float(frame)/float(nframes))*vy
                                    z += (1-2*float(frame)/float(nframes))*vz
##                                    x += (float(frame)/float(nframes))*vx
##                                    y += (float(frame)/float(nframes))*vy
##                                    z += (float(frame)/float(nframes))*vz
                                else:
                                    sqlength = 1
                                    continue

                                ## MODRES?
                                MODRES = False
                                if chain in d_MODRES.keys():
                                    if res_no in d_MODRES[chain].keys():
                                        if iCode in d_MODRES[chain][res_no].keys():
                                            if res_name != d_MODRES[chain][res_no][iCode]:
                                                print res_name, d_MODRES[chain][res_no][iCode]
                                                notexpected
                                            MODRES = True

                                ## ATOM or HETATM ?
                                if MODRES == False and res_name not in self.d_res.keys():
                                    record = 'HETATM'
                                elif MODRES == True:
                                    record = 'ATOM'
## add a flag to know if std res is atom or hetatm
                                elif MODRES == False and res_name in self.d_res.keys():
                                    record = 'ATOM'

                                ## append atom line
                                line = '%6s%5i  %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
                                    record.ljust(6),atom_no,
                                    atom_name.ljust(3), altloc, res_name, chain, res_no, iCode,
                                    x,y,z, 1.0, sqlength,atom_name[0].rjust(2)
                                    )
                                output_vmd += [line]

                output_vmd.append('TER\nENDMDL\n')
            fd = open(job+'_mode'+str(mode+1).zfill(2)+'.pdb', 'w') ## implicit path
            fd.writelines(output_vmd)
            fd.close()

        return


    def main(
        self, jobid, lines, atoms_hessian = ['CA'], frames = 50,
        cutoff_distance = 10.,
        path_html = None, path_python = None, verbose = False, paralleldir = '',
        chains = [],
        ):

        '''
        Use first model if no model specified by user.
        chain(s): Y, biomolecule: Y; parse chains specified by user and apply transformation
        chain(s): Y, biomolecule: N; parse chains specified by user but don't apply transformation
        chain(s): N, biomolecule: Y; parse chains of biomolecule and apply transformation
        chain(s): N, biomolecule: N; parse chains of first biomolecule and apply transformation
        '''

        import os, Numeric, sets

        ## parse pdb
        (
            d_REMARK350,
            d_primary, ## i.e. SEQRES, MODRES
            d_secondary, ## i.e. HELIX, SHEET
            d_coordinates, ## i.e. ATOM, HETATM, TER, MODEL, ENDMDL
            d_ligands,
            d_atomnos,
            ) = self.parse_pdb(lines, chains)

        ## assume multimeric biological unit if chains not specified by user
        if chains == []:
            chains = list(d_primary['SEQRES'])
            chains.sort()

        ## for the non-disrupted structure: calculate and visualize eigenvectors
        d_coordinates, d_hessian, l_CA = self.pre_hessian(d_coordinates, chains, atoms_hessian, d_ligands, d_atomnos)
        N = len(d_hessian.keys())
        matrix_hessian = self.hessian_calculation(N, d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, verbose = verbose)
        N = len(matrix_hessian)/3
        if len(d_hessian.keys()) != len(matrix_hessian)/3:
            stop
        eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_combined_nonperturbed = self.eigenv_calccomb(matrix_hessian, jobid, verbose)
        self.morph(eigenvectors_nonperturbed, frames, chains, d_coordinates, atoms_hessian, matrix_hessian, jobid, d_primary)

##        ##
##        ## do plots for individual modes
##        ##
##        self.plot(
##            jobid,
##            cutoff_distance,
##            eigenvectors_nonperturbed,
##            chains,d_secondary,
##            l_CA,
##            plottype = 'individual',
##            )
##        ##
##        ## do plots for combined modes
##        ##
##        self.plot(
##            jobid,
##            cutoff_distance,
##            eigenvectors_combined_nonperturbed,
##            chains,
##            d_secondary,
##            l_CA,
##            plottype = 'combined',
##            )
        
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

        ##
        ## loop over remres1
        ##
        for remres1 in range(len(l_CA)-1):
            remCA1 = l_CA[remres1]

            l_rem1 = self.CA2CB(remCA1,d_coordinates,d_hessian)

            ## continue the remres1 loop if another processor is running calculations for that value of remres1
            filename = '%s/%s_residue%s.txt' %(paralleldir, jobid, remres1+1)
            if os.path.isfile(filename):
## REWRITE THIS SECTION
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
                        datadic[key]['data'][mode][remres2][remres1] = float(line_main)
                    line_number_main += 1

                continue
            
            ## create data files (necessary to create data files if multithreading/parallel processing dependent on existence of files to determine if calculation on certain residues has been performed)
            else:
                fd = open(filename, 'w')
                fd.close()

            ##
            ## loop over remres2
            ##
            for remres2 in range(remres1+1,len(l_CA)):
                remCA2 = l_CA[remres2]

                l_rem2 = self.CA2CB(remCA2,d_coordinates,d_hessian)

                if not path_html:
                    print 'residue combination %s,%s of %s,%s' %(remres1+1,remres2+1,len(l_CA),len(l_CA))

                ## continue if windows overlap with borders of the matrix
                if remres1 < 0 or remres1 > (N-1) or remres2 < 0 or remres2 > (N-1):
                    stop1
                    continue
                ## continue if overlap values have already been calculated for 
                if remres2 < remres1:
                    stop2
                    continue

                print remCA1, l_rem1, d_hessian[remCA1]
                print remCA2, l_rem2, d_hessian[remCA2]
                l_rem = l_rem1 + l_rem2

                matrix_hessian_perturbed = self.hessian_calculation(
                    N, d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, verbose = verbose, l_rem = l_rem,
                    )
                eigenvectors_perturbed, eigenvalues_perturbed, eigenvectors_perturbed_combined = self.eigenv_calccomb(
                    matrix_hessian_perturbed, jobid, verbose
                    )
                overlaps_single, max_overlaps_single, perturbed_modes_of_max_overlap_single, delta_perturbed_eigenvalues_of_max_overlap = self.overlap_calculation(
                    eigenvectors_perturbed, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem
                    )
                overlaps_combined, max_overlaps_combined, perturbed_modes_of_max_overlap_combined = self.overlap_calculation(
                    eigenvectors_perturbed_combined, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem,
                    )[:-1]
                print overlaps_single[6]

        ##        ## do vmd of perturbed structure
        ##        self.morph(eigenvectors, frames, biomolecule, d_coordinates, atoms_hessian, cluster, matrix_hessian, jobid+'-'+str(xvalue)+'-'+str(yvalue), cutoff_distance)

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

            ## write data to txt files in case of crash during loop over residues/clusters
            lines = []
            for key in datadic:
                lines.append('data '+key+'\n')
                for mode in range(6,12):
                    lines.append('mode '+str(mode+1)+'\n')
                    for remres2 in range(len(datadic[key]['data'][mode][remres1])):
                        lines.append(str(datadic[key]['data'][mode][remres1][remres2])+'\n')
                if key == 'overlaps_combined':
                    mode = len(datadic[key]['data'])-1
                    lines.append('mode '+str(mode+1)+'\n')
                    for remres2 in range(len(datadic[key]['data'][mode][remres1])):
                        lines.append(str(datadic[key]['data'][mode][remres1][remres2])+'\n')
            fd = open(filename, 'w')
            fd.writelines(lines)
            fd.close()

################################################################################
##                                plot section                                ##
################################################################################

        if verbose == True:
            print 'generating plots'

        ##
        ## delta overlap vs residue number
        ##
        d_overlaps = {}
        for mode in range(6,12):
            l_overlaps = []
            for i in range(len(datadic['overlaps_single']['data'][mode])):
                if datadic['overlaps_single']['data'][mode][i][i] != 0:
                    stop
                overlap = sum(datadic['overlaps_single']['data'][mode][i])/(len(datadic['overlaps_single']['data'][mode][i])-1)
                l_overlaps += [overlap]
            plotname = 'average change in overlap'
            title1 = 'mode %s' %(mode+1)
            title = '"%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}"' %(plotname, title1, jobid, chains, cutoff_distance)
            self.gnuplot_plot(
                jobid, cutoff_distance, chains, title1,
                xtitle = 'residue', ytitle = 'average change in overlap', plotname = plotname,
                filename = '%s_deltaoverlap_mode%s' %(jobid, str(mode+1).zfill(2)),
                data = l_overlaps,
                title = title,
                )
            ## write overlaps to txt file
            lines = []
            for i in range(len(l_overlaps)):
                overlap = l_overlaps[i]
                lines += ['%3i %f\n' %(i, overlap)]
            fd = open('delta_overlaps_mode%s.txt' %(str(mode+1).zfill(2)),'w')
            fd.writelines(lines)
            fd.close()

            ## group data
            d_overlaps = {
                'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],'GLY':[],'HIS':[],
                'ILE':[],'LYS':[],'LEU':[],'MET':[],'ASN':[],'PRO':[],'GLN':[],
                'ARG':[],'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[],
                }
            if len(d_hessian.keys()) != len(l_overlaps):
                stop
            for i in range(len(l_overlaps)):
                res_name = d_hessian[i]['res_name']
                overlap = l_overlaps[i]
                d_overlaps[res_name] += [overlap]

            fileprefix = '%s_deltaoverlap_residue_mode%s' %(jobid, str(mode+1).zfill(2)),
            self.gnuplot_histogram(d_overlaps, fileprefix)
                        
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
                for i in range(N):
                    datadic[key]['data'][mode][i][i] = diagonal

        ##
        ## delta overlap
        ##
        xtitle = axistitle
        ytitle = axistitle
        ztitle = '{/Symbol D}overlap'
        fileprefix = jobid+'overlapscombined_mode3N'
        title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %('{/Symbol D}overlap', 'modes 7-3N', jobid, chains, cutoff_distance),

        ## overlaps_combined between nonperturbed and perturbed modes 7-3N
        self.gnuplot_splot(
            xtitle, ytitle, ztitle,
            fileprefix,
            resres_overlaps_combined[-1],
            title,
            z1 = '', z2 = '',
            )

        ## add topology to margin of cross-correlations maps
        fileprefix = jobid+'overlapscombined_mode3N'
        self.topology(
            fileprefix,chains,d_secondary,N,
            )
        
        for mode in range(6,12):

            for key in datadic:

                ztitle = datadic[key]['name']
                fileprefix = '%s%s_mode%s' %(jobid, key, str(mode+1).zfill(2))
                title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %(datadic[key]['name'], mode+1, jobid, chains, cutoff_distance),
                self.gnuplot_splot(
                    xtitle, ytitle, ztitle,
                    fileprefix,datadic[key]['data'][mode],title,
                    z1 = datadic[key]['zrange'][0], z2 = datadic[key]['zrange'][1],
                    )

                ## add topology to margin of cross-correlations maps
                self.topology(
                    fileprefix,chains,d_secondary, N,
                    )

        return


    def gnuplot_histogram(self, d_data, fileprefix):

        import math

        l_xtics = d_data.keys()
        l_xtics.sort()

        ##
        ## write data
        ##
        yrange = []
        gnuplotdata = []
        for i in range(len(l_xtics)):
            res_name = l_xtics[i]
            for y in d_data[res_name]:
                gnuplotdata += ['%f %f\n' %(float(i),y)]
                yrange += [y]
        ymin = min(yrange)
        ymax = max(yrange)
        fd = open('gnuplot.data','w')
        fd.writelines(gnuplotdata)
        fd.close()

        ##
        ## calculate statistics
        ##
        gnuplot_statistics = []
        for i in range(len(l_xtics)):
            res_name = l_xtics[i]
            n = len(d_data[res_name])
            if n <= 1:
                continue
            sumx = 0
            sumxx = 0
            for x in d_data[res_name]:
                sumx += x
                sumxx += x**2
            average = sumx/n
            SS = sumxx-(sumx**2)/n
            MSE = SS / (n-1)
            if MSE < 0:
                SE = 0 ## temp!!! check the equation!!!
            else:
                SE = math.sqrt(MSE/n)
            gnuplot_statistics += ['%f %f %f\n' %(float(i),average,SE)]
        ## write statistics
        fd = open('gnuplot.statistics','w')
        fd.writelines(gnuplot_statistics)
        fd.close()

        ##
        ## write gnuplot settings
        ##
        gnuplotsettings = []
        gnuplotsettings += [
            'set terminal postscript eps enhanced color "Helvetica" 48\n',
            'set output "gnuplot.ps"\n',
            'set size 4,4\n', ## scale 400%
            'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
            'set xlabel "classification"\n',
            'set ylabel "max overlap"\n',
            ]
        line_xtic = 'set xtics ('
        for xtic in l_xtics:
            line_xtic += '"%s" %s, ' %(xtic, l_xtics.index(xtic))
        line_xtic = line_xtic[:-2]+')\n'
        gnuplotsettings += [
            line_xtic,
            'set xtics rotate\n',
        ]
        gnuplotsettings += [
            'plot [-1:%i][%f:%f] "gnuplot.data" lt 0 ps 2 pt 2' %(len(l_xtics)+1, ymin, ymax),
            ', "gnuplot.statistics" lt 1 lc 0 ps 0 pt 0 w errorb\n',
            ]
        ## write gnuplot settings
        fd = open('gnuplot.settings','w')
        fd.writelines(gnuplotsettings)
        fd.close()

        ##
        ## execute gnuplot settings
        ##
        os.system('/software/bin/gnuplot gnuplot.settings')
        ## convert postscript to portable network graphics
        os.system('convert gnuplot.ps %s.png' %(fileprefix))
        os.remove('gnuplot.ps')
        os.remove('gnuplot.data')
        os.remove('gnuplot.settings')
        os.remove('gnuplot.statistics')

        return


    def CA2CB(self,remCA,d_coordinates,d_hessian):
        chain = d_hessian[remCA]['chain']
        res_no = d_hessian[remCA]['res_no']
        iCode = d_hessian[remCA]['iCode']
        res_name = d_hessian[remCA]['res_name']
        l_rem = []
        for atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys():
            altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
            if 'i' in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc].keys():
                i = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['i']
                l_rem += [i]
        return l_rem


    def parse_pdb(self, lines, parse_chains):

        import Numeric, sets

        d_coordinates = {'chains':{}} ## ATOM, HETATM, MODEL
        d_REMARK350 = {} 
        d_secondary = {'HELIX':{},'SHEET':{},} ## HELIX, SHEET
        d_primary = {'SEQRES':sets.Set(),'MODRES':{}}
        d_ligands = {'chains':{}}
        d_atomnos = {}
        accept_missing_residues = False
        accept_altloc = False

        for i in range(len(lines)):

            line = lines[i]

            record = line[:6].strip()

            if record == 'ATOM':
                altloc = line[16]
                if altloc != ' ' and accept_altloc == True:
                    print 'Alternative locations of atoms exist. Do you want to proceed (Yes/No)?'
                    proceed = raw_input()
                    if proceed.upper() in ['YES','Y']:
                        accept_altloc = True
                    elif proceed.upper() in ['NO','N']:
                        accept_altloc = False
                        raise
                    else:
                        raise 'invalid answer'

                d_coordinates, d_atomnos = self.parse_ATOM(line, d_coordinates, d_atomnos)
                continue

            elif record == 'REMARK':

                remark = int(line[6:10].strip())

                if remark == 350:

                    d_REMARK350 = self.parse_REMARK350(lines,i, d_REMARK350)

                elif remark == 465:

                    if accept_missing_residues == False:
                        print 'Residues are missing (i.e. REMARK465 records are present). Do you want to proceed (Yes/No)?'
                        proceed = raw_input()
                        if proceed.upper() in ['YES','Y']:
                            accept_missing_residues = True
                        elif proceed.upper() in ['NO','N']:
                            accept_missing_residues = False
                            raise
                        else:
                            raise 'invalid answer'

            elif record == 'SEQRES':
                chain = line[11]
                d_primary['SEQRES'] |= sets.Set(chain)
                                    
            elif record == 'MODRES':
                chain = line[16]
                res_no = int(line[18:22])
                iCode = line[22]
                res_name = line[12:15]
                if chain not in d_primary['MODRES'].keys():
                    d_primary['MODRES'][chain] = {}
                if res_no not in d_primary['MODRES'][chain].keys():
                    d_primary['MODRES'][chain][res_no] = {}
                if iCode not in d_primary['MODRES'][chain][res_no].keys():
                    d_primary['MODRES'][chain][res_no][iCode] = res_name

            elif record == 'HELIX':
                chain = line[19]
                if not d_secondary['HELIX'].has_key(chain):
                    d_secondary['HELIX'][chain] = {'range':[],'list':[]}
                d_secondary['HELIX'][chain]['range'] += [[int(line[21:25]), int(line[33:37])]]
                d_secondary['HELIX'][chain]['list'] += range(int(line[21:25]), int(line[33:37])+1)
                continue

            elif record == 'SHEET':
                chain = line[21]
                if not d_secondary['SHEET'].has_key(chain):
                    d_secondary['SHEET'][chain] = {'range':[]} 
                d_secondary['SHEET'][chain]['range'] += [[int(line[22:26]), int(line[33:37])]]
                continue

            elif record == 'MODEL':
                model = int(line[10:14])
                continue

            elif record == 'HETATM':
                MODRES = False
                chain = line[21]
                res_no = int(line[22:26].strip())
                iCode = line[26]
                res_name = line[17:20].strip()
                if res_name in ['HOH','DOD',]:
                    continue
                if chain in d_primary['MODRES'].keys():
                    if res_no in d_primary['MODRES'][chain].keys():
                        if iCode in d_primary['MODRES'][chain][res_no].keys():
                            if res_name != d_primary['MODRES'][chain][res_no][iCode]:
                                print res_name, d_primary['MODRES'][chain][res_no][iCode]
                                notexpected
                            MODRES = True
                if MODRES == True:
                    d_coordinates, d_atomnos = self.parse_ATOM(line, d_coordinates, d_atomnos)
                else:
                    d_ligands, d_atomnos = self.parse_ATOM(line,d_ligands, d_atomnos)

                continue

        return d_REMARK350, d_primary, d_secondary, d_coordinates, d_ligands, d_atomnos


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


    def parse_ATOM(self, line, d_coordinates, d_atomnos):

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
        occupancy = float(line[54:60])
        beta = float(line[60:66])
        
        element = line[76:78].strip()

        coordinate = Numeric.array([x, y, z])

        if not chain in d_coordinates['chains'].keys():
            d_coordinates['chains'][chain] = {'residues':{}}
        if not res_no in d_coordinates['chains'][chain]['residues'].keys():
            d_coordinates['chains'][chain]['residues'][res_no] = {'iCodes':{}}
        if not iCode in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys():
            d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode] = {'res_names':{},'altlocs':{}}
        if not res_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'].keys():
            d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name] = {'atoms':{}}
        if not atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys():
            d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name] = {'altlocs':{}}

        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc] = {
            'res_name':res_name
            }
        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc] = {
            'coordinate':coordinate,'element':element,'atom_no':atom_no,'occupancy':occupancy,'beta':beta,
            }

        d_atomnos[atom_no] = {'atom_name':atom_name,'coordinate':coordinate}

        return d_coordinates, d_atomnos


    def overlap_calculation(self, eigenvectors_perturbed, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem):

        ## calculate fewer overlaps to speed up things. all overlaps only calculated to determine mode of max overlap
        ## reduce moderanges to speed up things significantly

        import math, Numeric
##        l_rem.sort()
##        l_rem.reverse()
        overlaps = []
        max_overlaps = []
        perturbed_modes_of_max_overlap = []
        delta_perturbed_eigenvalues_of_max_overlap = []
##        moderanges_nonperturbed = range(6,18)+[len(eigenvectors_nonperturbed)-3*len(l_rem)-1]
##        moderanges_perturbed = range(6,12)+[len(eigenvectors_nonperturbed)-3*len(l_rem)-1]
        moderanges_nonperturbed = range(6,18)+[len(eigenvectors_nonperturbed)-1]
        moderanges_perturbed = range(6,12)+[len(eigenvectors_nonperturbed)-1]
        for mode_nonperturbed in range(len(eigenvectors_nonperturbed)):
            max_overlap = 0
            overlaps_per_mode_perturbed = []
            ## convert eigenvector of nonperturbed structure and delete appropiate coordinates
            vector_nonperturbed = list(eigenvectors_nonperturbed[mode_nonperturbed])
##            for i in l_rem:
##                for j in range(3-1,-1,-1):
##                    del vector_nonperturbed[3*i+j]
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


    def pre_hessian(self, d_coordinates, chains, atoms_hessian, d_ligands, d_atomnos):

        print 'calculating distances between atoms'

        d_hessian = {}
        l_CA = []
        i = 0

        ##
        ## vicinal ligands
        ##
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

                    ## set lists
                    d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_wt'] = []
                    d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_ala'] = []

                    ##
                    ## loop over atom_names and append atom_nos
                    ##
                    for atom_name1 in atom_names1:
                        altloc1 = min(d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'].keys())
                        coordinate1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['coordinate']
                        if atom_name1 == 'CA':
                            l_CA += [i]
                        if atom_name1 in atoms_hessian:
                            d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['i'] = i
                            d_hessian[i] = {'chain':chain1,'res_no':res_no1,'iCode':iCode1,'atom_name':atom_name1,'altloc':altloc1, 'res_name': res_name1}
                            i += 1
                        atom_no1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atoms'][atom_name1]['altlocs'][altloc1]['atom_no']
                        d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_wt'] += [atom_no1]
                        if atom_name1 in ['N','CA','C','O','OXT','CB',]:
                            d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_ala'] += [atom_no1]

                    ##
                    ## vicinal ligands
                    ##
                    for chain2 in d_ligands['chains'].keys():
                        ## assume sequential numbering of residues
                        res_nos2 = d_ligands['chains'][chain2]['residues'].keys()
                        res_nos2.sort()
                        for res_no2 in res_nos2:
                            ## assume sequential iCodes
                            iCodes2 = d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'].keys()
                            iCodes2.sort()
                            for iCode2 in iCodes2:
                                altloc_residue2 = min(d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['altlocs'].keys())
                                res_name2 = d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['altlocs'][altloc_residue2]['res_name']
                                atom_names2 = d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'].keys()
                                atom_names2.sort()

                                if res_name2 not in ['XYP']: ## temp!!
                                    continue

                                ##
                                ## calculate distances
                                ##
                                l_sqdist_wt = []
                                l_sqdist_ala = []
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
                                        altloc2 = min(d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'].keys())
                                        coordinate2 = d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'][altloc2]['coordinate']

                                        sqdist = sum((coordinate2-coordinate1)**2)
                                        l_sqdist_wt += [sqdist]
                                        if atom_name1 in ['N','CA','C','O','OXT','CB',]:
                                            l_sqdist_ala += [sqdist]

                                sqdist_wt = min(l_sqdist_wt)
                                sqdist_ala = min(l_sqdist_ala)
                                if sqdist_wt < 25: ## temp!!!
                                    for atom_name2 in atom_names2:
                                        altloc2 = min(d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'].keys())
                                        atom_no2 = d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'][altloc2]['atom_no']
                                        d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_wt'] += [atom_no2]

                                if sqdist_ala < 25: ## temp!!!
                                    for atom_name2 in atom_names2:
                                        altloc2 = min(d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'].keys())
                                        atom_no2 = d_ligands['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atoms'][atom_name2]['altlocs'][altloc2]['atom_no']
                                        d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_ala'] += [atom_no2]


        ##
        ## vicinal residues
        ##
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
                    atom_nos_wt1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_wt']
                    atom_nos_ala1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['atom_nos_ala']

                    ##
                    ## vicinal residues
                    ##
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
                                atom_nos_wt2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atom_nos_wt']
                                atom_nos_ala2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['atom_nos_ala']

                                if (
                                    chain1 == chain2 and
                                    res_no1 == res_no2 and
                                    iCode1 == iCode2
                                    ):
                                    continue

                                ##
                                ## calculate distances
                                ##
                                l_sqdist_wt = []
                                l_sqdist_ala = []
                                l_sqdist_ala1 = []
                                l_sqdist_ala2 = []
                                for atom_no1 in atom_nos_wt1:
                                    coordinate1 = d_atomnos[atom_no1]['coordinate']

                                    for atom_no2 in atom_nos_wt2:
                                        coordinate2 = d_atomnos[atom_no2]['coordinate']

                                        sqdist = sum((coordinate2-coordinate1)**2)
                                        if atom_no1 in atom_nos_ala1 and atom_no2 in atom_nos_ala2:
                                            l_sqdist_ala += [sqdist]
                                            l_sqdist_ala1 += [sqdist]
                                            l_sqdist_ala2 += [sqdist]
                                        elif atom_no1 in atom_nos_ala1:
                                            l_sqdist_ala1 += [sqdist]
                                        elif atom_no2 in atom_nos_ala2:
                                            l_sqdist_ala2 += [sqdist]
                                        l_sqdist_wt += [sqdist]

                                sqdist_wt = min(l_sqdist_wt)
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
                                        'sqdist_wt':sqdist_wt,'sqdist_ala':sqdist_ala,'sqdist_ala1':sqdist_ala1,'sqdist_ala2':sqdist_ala2,
                                        }

        return d_coordinates, d_hessian, l_CA
        

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
                                                    sqdist_wt = 0
                                                    sqdist_ala = 0
                                                    sqdist_ala1 = 0
                                                    sqdist_ala2 = 0
                                                else:
                                                    d_sqdist1 = d_coordinates['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]['sqdist']['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]
                                                    sqdist_wt1 = d_sqdist1['sqdist_wt']
                                                    sqdist_ala1 = d_sqdist1['sqdist_ala']
                                                    sqdist_ala11 = d_sqdist1['sqdist_ala1']
                                                    sqdist_ala21 = d_sqdist1['sqdist_ala2']
                                                    d_sqdist2 = d_coordinates['chains'][chain2]['residues'][res_no2]['iCodes'][iCode2]['res_names'][res_name2]['sqdist']['chains'][chain1]['residues'][res_no1]['iCodes'][iCode1]['res_names'][res_name1]
                                                    sqdist_wt2 = d_sqdist2['sqdist_wt']
                                                    sqdist_ala2 = d_sqdist2['sqdist_ala']
                                                    sqdist_ala12 = d_sqdist2['sqdist_ala1']
                                                    sqdist_ala22 = d_sqdist2['sqdist_ala2']

                                                    if sqdist_wt1 != sqdist_wt2 or sqdist_ala1 != sqdist_ala2 or sqdist_ala11 != sqdist_ala22 or sqdist_ala21 != sqdist_ala12:
                                                        print chain1, res_no1, iCode1, res_name1, altloc1, atom_name1
                                                        print chain2, res_no2, iCode2, res_name2, altloc2, atom_name2
                                                        print sqdist_main1, sqdist_side1, sqdist_ala1, sqdist_ala11, sqdist_ala21
                                                        print sqdist_main2, sqdist_side2, sqdist_ala2, sqdist_ala12, sqdist_ala22
                                                        import math
                                                        print math.sqrt(sqdist_main1), math.sqrt(sqdist_side1)
                                                        print math.sqrt(sqdist_main2), math.sqrt(sqdist_side2)
                                                        notexpected
                                                    else:
                                                        sqdist_wt = sqdist_wt1
                                                        sqdist_ala = sqdist_ala1
                                                        sqdist_ala1 = sqdist_ala11
                                                        sqdist_ala2 = sqdist_ala21

                                                    if sqdist_ala1 < sqdist_wt:
                                                        stop1
                                                    if sqdist_ala2 < sqdist_wt:
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
                                                    sqdist = sqdist_wt
                                                elif row_sup in l_rem and col_sup in l_rem:
                                                    sqdist = sqdist_ala
                                                elif row_sup in l_rem:
                                                    sqdist = sqdist_ala1
                                                elif col_sup in l_rem:
                                                    sqdist = sqdist_ala2
                                                else:
                                                    sqdist = sqdist_wt

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

        if verbose == True:
            print 'calculated eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'

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
        self, fileprefix, chains, d_secondary, axlen
        ):

        import os

        bl = [390,833] ## bottom left coordinate
        tr = [1050,173] ## top right coordinate
        s = 14 ## space between plot and topology
        w = (tr[0]-bl[0])/axlen ## width of squares

        lines = 'convert %s.png' %(fileprefix)
        for chain in chains:
            for secelm in d_secondary.keys():
                if secelm == 'SHEET':
                    lines += ' -fill "rgb(0,0,255)" '
                elif secelm == 'HELIX':
                    lines += ' -fill "rgb(255,0,0)" '
                else:
                    continue
                if chain in d_secondary[secelm]:
                    for secelmrange in d_secondary[secelm][chain]['range']:
                        ## bottom
                        lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), bl[1]+s, int(bl[0]+w*secelmrange[1]), bl[1]+s)
                        ## top
                        lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), tr[1]-s, int(bl[0]+w*secelmrange[1]), tr[1]-s)
                        ## left
                        lines += ' -draw "line %s,%s %s,%s" ' %(bl[0]-s, int(bl[1]-w*(secelmrange[0]-1)), bl[0]-s, int(bl[1]-w*secelmrange[1]))
                        ## right
                        lines += ' -draw "line %s,%s %s,%s" ' %(tr[0]+s, int(bl[1]-w*(secelmrange[0]-1)), tr[0]+s, int(bl[1]-w*secelmrange[1]))
        if lines != 'convert %s.png' %(fileprefix):
            lines += '%s.png' %(fileprefix)
            os.system(lines)
                
        return


    def plot(
        self, jobid, cutoff_distance, eigenvectors, chains, d_secondary,l_CA,
        plottype,
        ):

        import math, Numeric

        for mode in range(6,12):

            print 'plotting mode %s' %(mode+1)

            ## calculate residue fluctuation (atom,RMSF) for each mode

            lengths = []
            for i in range(0,len(eigenvectors[mode]),3):
                if i/3 in l_CA:
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
            fd = open('%s_%s_%s_lengths.txt' %(jobid,plottype,mode+1), 'w')
            fd.writelines(lines)
            fd.close()

            ## plot residue fluctuation (atom,RMSF) for each mode

            if plottype == 'individual':
                title1 = 'mode %s' %(mode+1)
            elif plottype == 'combined':
                title1 = 'modes 7-%s' %(mode+1)

            plotname = 'residue displacements'
            title = '"%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}"' %(plotname, title1, jobid, chains, cutoff_distance)
            self.gnuplot_plot(
                jobid, cutoff_distance, chains, title1,
                xtitle = 'residue', ytitle = 'RMSF', plotname = plotname,
                filename = '%s_%sdisplacement_mode%s' %(jobid, plottype, str(mode+1).zfill(2)),
                data = lengths, title = title,
                )

            ## prepare data for cross correlation maps
            ccdata = Numeric.zeros((len(eigenvectors[mode])/3,len(eigenvectors[mode])/3), typecode='d')
            for ccrow in range(len(ccdata)):
                for cccol in range(len(ccdata)):
                    ccdata[ccrow][cccol] = self.cosangle(eigenvectors[mode][3*ccrow:3*ccrow+3], eigenvectors[mode][3*cccol:3*cccol+3])

            ## plot cross-correlation map (residue,residue,cross-correlation) for each mode
            fileprefix = '%s_%s_crosscorrelation_mode%s' %(jobid, plottype, str(mode+1).zfill(2))

            for i in range(len(ccdata)):
                ccdata[i][i] = 1
            xtitle = 'residue'
            ytitle = 'residue'
            ztitle = 'cross correlation'
            title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %('cross correlation map', mode+1, jobid, chains, cutoff_distance),

            self.gnuplot_splot(
                xtitle, ytitle, ztitle,
                fileprefix,ccdata,title,
                z1 = -1, z2 = 1,
                )

            ## add topology to margin of cross-correlations maps
            self.topology(
                fileprefix, chains, d_secondary, len(ccdata),
                )

        return


    def gnuplot_splot(
        self, xtitle, ytitle, ztitle, fileprefix, data,
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
        fd = open('%s.gnuplotdata' %(fileprefix), 'w')
        fd.writelines(gnuplot_splot_data)
        fd.close()
        ## write gnuplot settings to txt file
        lines = ['set size square\n'] ## scale square
        lines += [
            'set terminal postscript eps enhanced color "Helvetica" 48\n',
            'set output "%s.ps"\n' %(fileprefix),
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
            'splot "%s.gnuplotdata" title ""\n' %(fileprefix), ## splot gnuplot data file
            ## set ticslevel
            ]

        fd = open('%s.gnuplotsettings' %(fileprefix), 'w')
        fd.writelines(lines)
        fd.close()
        ## plot data with gnuplot splot
        os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(fileprefix))
        ## convert postscript to portable network graphics
        os.system('convert %s.ps %s.png' %(fileprefix, fileprefix))
        os.remove('%s.ps' %(fileprefix))
        os.remove('%s.gnuplotdata' %(fileprefix))
        os.remove('%s.gnuplotsettings' %(fileprefix))

        return


    def gnuplot_plot(self, jobid, cutoff_distance, chains, title1, xtitle, ytitle, plotname, filename, data, title):

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
            'set title %s ,4\n' %(title),
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

    import os, sys

    if '-pdb' not in sys.argv:
        raise 'use -pdb to select a pdb'
    pdb = sys.argv[sys.argv.index('-pdb')+1][:4]
    if '-chains' in sys.argv:
        chains = sys.argv[sys.argv.index('-chains')+1].split(',')
    else:
        chains = []
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
    lines = instance_vibration.pdb_import(pdb)

    instance_vibration.main(
        pdb,
        lines,
        chains = chains,
        atoms_hessian = atoms,
        cutoff_distance = cutoff_distance,
        verbose = True,
        paralleldir = os.getcwd(),
        )

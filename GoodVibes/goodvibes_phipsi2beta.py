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


    def morph(
        self, eigenvectors_all_modes, nframes, chains, d_coordinates,
        atoms_hessian, matrix_hessian, job,
        ):

	'''This function puts in frames between the maximum amplitudes of a
movement given by eigenvectors. The values of the diagonal elements of the
hessian matrix are used for B factors to make coloring during simulation
possible.'''

        sys.path.append('/home/people/tc/svn/tc_sandbox/quakes/')
        import phipsi_comparison
        sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
        import phipsi

        d_residues = {
            'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
            'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
            'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
            'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
            }

        ##d_ramachandran = phipsi_comparison.parse_dihedrals(d_residues)
        d_ramachandran = {}
        for res in d_residues.values():
            d_ramachandran[res] = {}
            for phi in range(-180,180):
                d_ramachandran[res][phi] = {}
                for psi in range(-180,180):
                    d_ramachandran[res][phi][psi] = 0
        d_ramachandran = phipsi_comparison.parse_dihedrals(d_residues) ## slow

        d_max = {}
        for res in d_ramachandran.keys():
            max = 0
            for phi in d_ramachandran[res].keys():
                for psi in d_ramachandran[res][phi].keys():
                    count = d_ramachandran[res][phi][psi]
                    if count > max:
                        max = count
            d_max[res] = max

##        print 'visualize the two extreme projections along a trajectory and interpolate n frames between them'

        fd = open('/oxygenase_local/data/pdb/lz/pdb2lzm.ent','r')
        lines = fd.readlines()
        fd.close()
        ss_lines = []
        for line in lines:
            record = line[:6].strip()
            if record in ['HELIX','TURN','SHEET',]:
                ss_lines += [line]

        for mode in range(6,min(12,len(matrix_hessian))):

            print mode

            eigenvectors = eigenvectors_all_modes[mode]

##            eigenvectors = length_adjustment(mode, eigenvectors)

            output_vmd = ['REMARK color by connectivity (b-factor) or squared displacement (temperature factor)\n']
            for frame in range(nframes):

                lines_rasmol = [
                    'rasmol -nodisplay NMAframe%s.pdb << EOF\n' %(str(frame).zfill(2)),
                    'cartoon\n',
                    'wireframe 0\n',
                    ]

                prefix = 'NMA%s' %(frame)
                ## write rasmol script
                lines = [
                    'rasmol -nodisplay NMAframe%s.pdb << EOF\n' %(str(frame).zfill(2)),
                    'cartoon\n',
                    'wireframe 0\n',
##                    'ribbons\n',
##                    'color temperature\n',
                    ]

                d_colors = {}
                print mode,frame
                output_rasmol = []
                output_rasmol = list(ss_lines)
                output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
                for chain in chains:
                    res_nos = d_coordinates['chains'][chain]['residues'].keys()
                    res_nos.sort()
                    for i_res_no in range(1,len(res_nos)-1):
                        res_no = res_nos[i_res_no]
                        for iCode in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys():
                            for altloc in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].keys():
                                        res_name = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['res_name']
##                                for atom_name in ['C','CA','N','O',]:
##                                for atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'].keys():
##                                    if 'hessian' in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'][atom_name].keys():
                                        C = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms']['C']['coordinate']
                                        CA = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms']['CA']['coordinate']
                                        N = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms']['N']['coordinate']
                                        O = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms']['O']['coordinate']

                                        vx_C = 64*eigenvectors[4*3*(i_res_no)+0]
                                        vy_C = 64*eigenvectors[4*3*(i_res_no)+1]
                                        vz_C = 64*eigenvectors[4*3*(i_res_no)+2]
                                        vx_CA = 64*eigenvectors[4*3*(i_res_no)+3]
                                        vy_CA = 64*eigenvectors[4*3*(i_res_no)+4]
                                        vz_CA = 64*eigenvectors[4*3*(i_res_no)+5]
                                        vx_N = 64*eigenvectors[4*3*(i_res_no)+6]
                                        vy_N = 64*eigenvectors[4*3*(i_res_no)+7]
                                        vz_N = 64*eigenvectors[4*3*(i_res_no)+8]
                                        vx_O = 64*eigenvectors[4*3*(i_res_no)+9]
                                        vy_O = 64*eigenvectors[4*3*(i_res_no)+10]
                                        vz_O = 64*eigenvectors[4*3*(i_res_no)+11]

                                        x2_C = C[0]+(1-2*float(frame)/float(nframes))*vx_C
                                        y2_C = C[1]+(1-2*float(frame)/float(nframes))*vy_C
                                        z2_C = C[2]+(1-2*float(frame)/float(nframes))*vz_C
                                        x2_CA = CA[0]+(1-2*float(frame)/float(nframes))*vx_CA
                                        y2_CA = CA[1]+(1-2*float(frame)/float(nframes))*vy_CA
                                        z2_CA = CA[2]+(1-2*float(frame)/float(nframes))*vz_CA
                                        x2_N = N[0]+(1-2*float(frame)/float(nframes))*vx_N
                                        y2_N = N[1]+(1-2*float(frame)/float(nframes))*vy_N
                                        z2_N = N[2]+(1-2*float(frame)/float(nframes))*vz_N
                                        x2_O = O[0]+(1-2*float(frame)/float(nframes))*vx_O
                                        y2_O = O[1]+(1-2*float(frame)/float(nframes))*vy_O
                                        z2_O = O[2]+(1-2*float(frame)/float(nframes))*vz_O

                                        C2 = [x2_C,y2_C,z2_C,]
                                        CA2 = [x2_CA,y2_CA,z2_CA,]
                                        N2 = [x2_N,y2_N,z2_N,]
                                        O2 = [x2_O,y2_O,z2_O,]

                                        C_prev = d_coordinates['chains'][chain]['residues'][res_no-1]['iCodes'][iCode]['altlocs'][altloc]['atoms']['C']['coordinate']
                                        vx_C_prev = 64*eigenvectors[4*3*(i_res_no-1)+0]
                                        vy_C_prev = 64*eigenvectors[4*3*(i_res_no-1)+1]
                                        vz_C_prev = 64*eigenvectors[4*3*(i_res_no-1)+2]
                                        x2_C_prev = C_prev[0]+(1-2*float(frame)/float(nframes))*vx_C_prev
                                        y2_C_prev = C_prev[1]+(1-2*float(frame)/float(nframes))*vy_C_prev
                                        z2_C_prev = C_prev[2]+(1-2*float(frame)/float(nframes))*vz_C_prev
                                        C_prev2 = [x2_C_prev,y2_C_prev,z2_C_prev,]
                                        
                                        phi = int(self.dihedral(C_prev2,N2,CA2,C2,))

                                        N_next = d_coordinates['chains'][chain]['residues'][res_no+1]['iCodes'][iCode]['altlocs'][altloc]['atoms']['N']['coordinate']
                                        vx_N_next = 64*eigenvectors[4*3*(i_res_no+1)+0]
                                        vy_N_next = 64*eigenvectors[4*3*(i_res_no+1)+1]
                                        vz_N_next = 64*eigenvectors[4*3*(i_res_no+1)+2]
                                        x2_N_next = N_next[0]+(1-2*float(frame)/float(nframes))*vx_N_next
                                        y2_N_next = N_next[1]+(1-2*float(frame)/float(nframes))*vy_N_next
                                        z2_N_next = N_next[2]+(1-2*float(frame)/float(nframes))*vz_N_next
                                        N_next2 = [x2_N_next,y2_N_next,z2_N_next,]
                                        psi = int(self.dihedral(N2,CA2,C2,N_next2,))

                                        bfactor = (
                                            100
                                            *d_ramachandran[d_residues[res_name]][phi][psi]
                                            /float(d_max[d_residues[res_name]])
                                            )
##                                        index = int(100*bfactor/d_factors[mode][0])
##
##                                        if index > 100 or index < 0:
##                                            print index,bfactor
##                                            stop

                                        h = 159.+80*bfactor/100.
                                        s = 240.
                                        l = 120.+120.*(50-abs(bfactor-50))/50.
                                        r,g,b = self.hsl2rgb(h,s,l,)
                                        lines_rasmol += [
                                            'select %i\n' %(res_no),
                                            'color [%i,%i,%i]\n' %(
                                                r, g, b,
                                                )
                                            ]

                                        line_out = 'ATOM         %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(
                                            'C  ', altloc, res_name, chain, int(res_no), iCode,
                                            x2_C,y2_C,z2_C, bfactor, bfactor,
                                            )
                                        output_rasmol += [line_out]
                                        line_out = 'ATOM         %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(
                                            'CA ', altloc, res_name, chain, int(res_no), iCode,
                                            x2_CA,y2_CA,z2_CA, bfactor, bfactor,
                                            )
                                        output_rasmol += [line_out]
                                        line_out = 'ATOM         %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(
                                            'N  ', altloc, res_name, chain, int(res_no), iCode,
                                            x2_N,y2_N,z2_N, bfactor, bfactor,
                                            )
                                        output_rasmol += [line_out]
                                        line_out = 'ATOM         %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n' %(
                                            'O  ', altloc, res_name, chain, int(res_no), iCode,
                                            x2_O,y2_O,z2_O, bfactor, bfactor,
                                            )
                                        output_rasmol += [line_out]

                fd = open('NMAframe%s.pdb' %(str(frame).zfill(2)),'w')
                fd.writelines(output_rasmol)
                fd.close()

                output_vmd.append('TER\nENDMDL\n')


                lines_rasmol += [
        ##                    'ribbons\n',
                    'rotate x 100\n',
                    'rotate z 30\n',
                    'rotate x 100\n',
                    'rotate z 90\n',
                    'rotate x 40\n',
                    'rotate y -20\n',
                    'write %s.ppm\n' %(prefix),
                    'exit\n',
                    ]
                ## write rasmol script to file
                fd = open('%srasmol.src' %(prefix),'w')
                fd.writelines(lines_rasmol)
                fd.close()
                ## execute rasmol script
                os.system('source %srasmol.src > %srasmol.log' %(prefix,prefix))
                ## convert rasmol output
        ##                os.system('convert %s.ppm -resize x80 %s.gif' %(prefix,prefix))
                os.system('convert %s.ppm %s.gif' %(prefix,prefix))
                ## clean up
                os.remove('%s.ppm' %(prefix))
                os.remove('%srasmol.log' %(prefix))
                os.remove('%srasmol.src' %(prefix))
##                os.remove('NMAframe%s.pdb' %(str(frame).zfill(2)))


            line = 'convert '
            for frame in range(50):
                line += 'NMA%s.gif ' %(frame)
            for frame in range(50-1-1,-1+1,-1):
                line += 'NMA%s.gif ' %(frame)
            line += '-loop 0 mode%s.gif' %(mode)
            print line
            os.system(line)
            for frame in range(50):
                os.remove('NMA%s.gif' %(frame))

            fd = open(job+'_mode'+str(mode+1).zfill(2)+'.pdb', 'w') ## implicit path
            fd.writelines(output_vmd)
            fd.close()

        return


    def dihedral(self,c1,c2,c3,c4):

        import numpy,math

        v1 = numpy.array(c2)-numpy.array(c1)
        v2 = numpy.array(c3)-numpy.array(c2)
        v3 = numpy.array(c4)-numpy.array(c3)

        angle = math.atan2(
            numpy.dot(
                math.sqrt(sum(v2*v2))*v1,
                self.cross(v2,v3),
                ),
            numpy.dot(
                self.cross(v1,v2),
                self.cross(v2,v3),
                ),
            )
        angle *= 180./math.pi

        return angle


    def cross(self,v1,v2):

        import numpy

        n = numpy.array([
            v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0],
            ])

        return n


    def hsl2rgb(self,h,s,l):

        h = h/240.
        s = s/240.
        l = l/240.

        if l < .5:
            q = l+l*s
        else:
            q = l+s-l*s
        p = 2*l-q

        tr = h+1./3
        tg = h
        tb = h-1./3

        l_rgb = [0,0,0]
        l_tc = [tr,tg,tb]
        for i in range(len(l_tc)):
            tc = l_tc[i]
            if tc < 0:
                tc += 1
            if tc > 1:
                tc -= 1

            if tc < 1./6:
                c = p+(q-p)*6*tc
            elif 1./6 <= tc and tc < .5:
                c = q
            elif .5 <= tc and tc < 2./3.:
                c = p+(q-p)*(2./3-tc)*6
            else:
                c = p

            l_rgb[i] = c

        r = int(l_rgb[0]*255)
        g = int(l_rgb[1]*255)
        b = int(l_rgb[2]*255)

        return r,g,b


    def main(
        self, jobid, lines, atoms_hessian = [], frames = 51,
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
                    atom_names = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc]['atoms'].keys()
                    atom_names.sort()
                    for atom_name in atom_names:
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

        self.morph(
            eigenvectors_nonperturbed, frames, chains, d_coordinates, atoms_hessian, matrix_hessian, jobid,
            )

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
            print 'calculating the %sx%s hessian/second-order derivate Kirchoff/Laplacian matrix' %(len(coordinates),len(coordinates),)
        
        import math, Numeric, time

        cutoff_sq = cutoff**2

        N = len(coordinates)

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

        import numpy, math

        if verbose == True:
            print 'calculating eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'
        
        ## diagonalize hessian matrix
        import time

        eigen_tuple = numpy.linalg.eigh(matrix_hessian)
        ## parse eigenvalues and eigenvectors
        eigenvalues = list(eigen_tuple[0])
        eigenvectors = list(eigen_tuple[1])
        eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
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
        len7 = math.sqrt(sum(numpy.array(eigenvectors[6])*numpy.array(eigenvectors[6])))

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
        atoms_hessian = ['C','CA','N','O',],
        )

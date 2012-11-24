#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

import os, sys, Numeric, LinearAlgebra, math, urllib2

class quakes:

    def main(self):

        l_pdbs = []
        files = os.listdir(os.getcwd())
        for file in files:
            print file
            if file[:9] != 'bfaccorr_':
                continue
            else:
                fd = open(file,'r')
                lines = fd.readlines()
                fd.close()
                for line in lines:
                    pdb = line[:4]
                    l_pdbs += [pdb]

        fd = open('CathDomainList.v3.1.0.txt')
        lines = fd.readlines()
        fd.close()
        d_cath = {}
        for i in range(len(lines)):
            line = lines[i]
            if i % 1000 == 0:
                print i, len(lines)
            if line[0] == '#':
                continue
            pdb = line[:4]
            if pdb not in l_pdbs:
                continue
            C = int(line.split()[1])
            A = int(line.split()[2])
            CA = '%i.%i' %(C,A)
            if pdb not in d_cath.keys():
                d_cath[pdb] = []
            d_cath[pdb] += [CA]

        d_overlap = {}
        files = os.listdir(os.getcwd())
        for file in files:
            print file
            if file[:9] != 'bfaccorr_':
                continue
            else:
                fd = open(file,'r')
                lines = fd.readlines()
                fd.close()
                for line in lines:
                    pdb = line[:4]
                    overlap = line.split(';')[1]
                    if pdb not in d_cath.keys():
                        continue
                    l_CA = d_cath[pdb]
                    for CA in l_CA:
                        if not CA in d_overlap.keys():
                            d_overlap[CA] = []
                        d_overlap[CA] += [overlap]

        for CA in d_overlap.keys():
            f_sum = 0
            for overlap in d_overlap[CA]:
                f_sum += float(overlap)
            print '%s %4i %.2f' %(CA, len(d_overlap[CA]), f_sum/len(d_overlap[CA]))
        stop


        d_overlap = {}
        files = os.listdir(os.getcwd())
        for file in files:
            print file
            if file[:9] != 'bfaccorr_':
                continue
            else:
                fd = open(file,'r')
                lines = fd.readlines()
                fd.close()
                for line in lines:
                    pdb = line[:4]
                    overlap = line.split(';')[1]
                    fd = open('%s%s/pdb%s.ent' %(self.path_pdb, pdb[1:3], pdb), 'r')
                    lines = fd.readlines()
                    fd.close()
                    for line in lines:
                        if line[:6] == 'CRYST1':
                            sGroup = line[55:66].strip()
                    if not sGroup in d_overlap.keys():
                        d_overlap[sGroup] = []
                    d_overlap[sGroup] += [overlap]
        print d_overlap
        for sGroup in d_overlap.keys():
            f_sum = 0
            for overlap in d_overlap[sGroup]:
                f_sum += float(overlap)
            print '%10s %4i %.2f' %(sGroup.ljust(10), len(d_overlap[sGroup]), f_sum/len(d_overlap[sGroup]))
        stop


        omax = [0,0]
        files = os.listdir(os.getcwd())
        for file in files:
            print file
            if file[:9] != 'bfaccorr_':
                continue
            else:
                fd = open(file,'r')
                lines = fd.readlines()
                fd.close()
                for line in lines:
                    pdb = line[:4]
                    if pdb in ['5rxn']:
                        continue
                    overlap = float(line.split(';')[1])
                    if overlap > omax[0]:
                        omax = [overlap,pdb]
        print omax
        stop


        subdirs = os.listdir(self.path_pdb)
        subdirs.sort()
        for subdir in subdirs:
            if subdir < sys.argv[1]:
                continue
            if os.path.isfile('bfaccorr_%s.txt' %(subdir)):
                continue
            fd = open('bfaccorr_%s.txt' %(subdir),'w')
            fd.write('')
            fd.close()
            print subdir
            lines_out = []
            files = os.listdir(self.path_pdb+subdir)
            for file in files:
                fd = open('%s%s/%s' %(self.path_pdb, subdir, file), 'r')
                lines = fd.readlines()
                fd.close()
                pdb = file[3:7]
                if pdb in [
                    ## hetID not in MODRES records
                    '1aob','1det','1cv7',
                    ## hetID differes from MODRES record
                    '308d',
                    ## REMARK 2 records of X-ray structure used for other purposes than resolution
                    '190d','1dn6',
                    ]:
                    continue
##                print pdb
                if pdb == '273d':
                    sprsde
                d_header = {}
                d_coordinates = {}
                s_break = False
                for i in range(len(lines)):
                    line = lines[i]

                    record = line[:6].strip()

                    if record == 'HEADER':
                        d_header['HEADER'] = line[10:50].strip()
                    elif record in ['TITLE','CAVEAT','COMPND','SOURCE','KEYWDS',]:
                        continue
                    elif record == 'EXPDTA': ## section 2
                        d_header = self.parse_recordEXPDTA(line, d_header, pdb)
                        if d_header['EXPDTA'] != 'DIFFRACTION':
                            s_break = True
                            break
                    elif record in ['AUTHOR','REVDAT','SPRSDE','JRNL',]:
                        continue
                    elif record == 'REMARK': ## section 2
                        remark = int(line[6:10])
                        if remark in [0,1,]:
                            continue
                        elif remark == 2:
                            if line[11:21] == 'RESOLUTION' and line[23:37] not in ['NOT APPLICABLE','NULL ANGSTROMS',]:
                                try:
                                    resolution = float(line[22:27])
                                except:
                                    print line[11:21]
                                    print line
                                    if line[11:21] != '          ':
                                        stop
                            elif line[11:21] != '          ' and line[23:37] not in ['NOT APPLICABLE','NULL ANGSTROMS',]:
                                fd = open('remark2.txt','a')
                                fd.write('%s,' %(pdb))
                                fd.close()
                        elif remark in [3,4,5,6,7,9,10,11,12,13,14,15,16,42,99,]:
                            continue
                        elif remark in [
                            100,102,105,106,200,210,215,
                            230, ## experimental details
                            240, ## experimental details
                            245,247,250,
                            265, ## experimental details
                            280,285,290,295,300,350,351,375,376,400,
                            ]:
                            continue
                        elif remark == 465:
                            s_break = True
                            break
                        elif remark == 470:
                            d_header, alpha = self.parse_recordREMARK470(d_header, pdb, lines, i)
                            if alpha == True:
                                s_break = True
                                break
                        elif remark in [
                            475,480,500,525,600,610,
                            615, ## zero occupancy atom
                            620,650,700,800,900,950,999,
                            ]:
                            continue
                        else:
                            print lines[i-1], lines[i], lines[i+1]
                            stop
                    elif record in ['DBREF','SEQADV',]:
                        continue
                    elif record == 'SEQRES': ## section 3
                        d_header = self.parse_recordSEQRES(line, d_header, pdb)
                    elif record == 'FTNOTE':
                        continue
                    elif record in [
                        'HET','HETNAM','HETSYN','FORMUL', ## section 4
                        'HELIX','SHEET','TURN',
                        'SSBOND','LINK','SLTBRG','HYDBND','CISPEP','SITE','CRYST1',
                        'ORIGX1','ORIGX2','ORIGX3',
                        'SCALE1','SCALE2','SCALE3','TVECT',
                        'MTRIX1','MTRIX2','MTRIX3',
                        ]:
                        continue
                    elif record == 'MODRES':
                        s_break = True
                        break
                    elif record in ['MODEL','ENDMDL',]:
                        s_break = True
                        break
##                        d_header = self.parse_recordMODRES(d_header,line)
                    elif record == 'ATOM': ## section 9
                        d_coordinates = self.parse_recordATOM(line, d_coordinates, d_header, record)
                    elif record in ['SIGATM','TER',]:
                        continue
                    elif record == 'HETATM': ## section 9
                        d_coordinates = self.parse_recordATOM(line, d_coordinates, d_header, record)
                    elif record in ['ANISOU','SIGUIJ',]:
                        continue
                    elif record == 'CONECT':
                        continue
                    elif record == 'MASTER':
                        continue
                    elif record == 'END':
                        continue
                    else:
                        print lines[i-1], lines[i], lines[i+1], lines[i+2], lines[i+3], lines[i+4], lines[i+5]
                        print d_header
                        print pdb
                        stop

                if s_break == True:
                    continue

                if 'SEQRES' not in d_header.keys():
                    continue
                else:
                    chains = d_header['SEQRES'].keys()
                if len(chains) > 1:
                    continue
                else:
                    chain = chains[0]

                if d_header['SEQRES'][chain]['type'] != 'peptide':
                    continue

                l_coordinates = []
                l_bfactors = []
                i = 0
                l_res_nos = d_coordinates['chains'][chain]['residues'].keys()
                l_res_nos.sort()
                for res_no in l_res_nos:
                    l_iCodes = d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes']
                    for iCode in l_iCodes:
                        i += 1
                        if 'CA' not in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                            if len(l_coordinates) != 0: ## and len(l_coordinates)+1 != len(d_header['SEQRES'][chain]['seq']):
                                ## C-terminal missing alpha carbon atom
##                                if len(l_coordinates)+1 == len(d_header['SEQRES'][chain]['seq']) and 'CA' in d_header['REMARK470'][chain][res_no][iCode]:
##                                    continue
                                print len(l_coordinates), len(d_header['SEQRES'][chain]['seq'])
                                print pdb, chain, res_no, iCode
                                stop
                            ## N-terminal modification
                            elif len(l_coordinates) == 0 and d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] not in self.d_res.keys():
                                continue
                        else:
                            coordinate = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA']['coordinate']
##                            if d_header['SEQRES'][chain]['seq'][i-1] != d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']:
##                                print pdb, i, d_header['SEQRES'][chain]['seq'][i-1], d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
##                                print pdb, chain, res_no, iCode
##                                stop
                        bfactor = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']['CA']['bfactor']
                        l_coordinates += [coordinate]
                        l_bfactors += [bfactor]

                if l_bfactors == len(l_bfactors)*[l_bfactors[0]]:
                    continue

                if len(l_coordinates) != len(d_header['SEQRES'][chain]['seq']):
                    ## N-terminal modification
                    try:
                        if len(l_coordinates) == len(d_header['SEQRES'][chain]['seq'])-1 and d_header['SEQRES'][chain]['seq'][0] not in self.d_res.keys():
                            pass
                        elif 'CA' in d_header['REMARK470'][chain][res_no][iCode]:
                            pass
                        else:
                            print d_header['SEQRES'][chain]['seq']
                            print len(l_coordinates), len(d_header['SEQRES'][chain]['seq'])
                            print pdb
                            stop
                    except:
                        print d_header['SEQRES'][chain]['seq']
                        print len(l_coordinates), len(d_header['SEQRES'][chain]['seq'])
                        print pdb
                        stop
                        

                matrix_hessian = self.hessian_calculation(l_coordinates, 10., verbose = True)
                eigenvectors, eigenvalues, eigenvectors_combined = self.eigenv_calccomb(matrix_hessian, verbose = True)
                line_out = '%s;' %(pdb)
                for mode in range(6,12):
                    l_rmsfs = []
                    for i in range(0,len(eigenvectors[mode]),3):
                        x = eigenvectors[mode][i+0]
                        y = eigenvectors[mode][i+1]
                        z = eigenvectors[mode][i+2]
                        rmsf = math.sqrt(x**2+y**2+z**2)
                        l_rmsfs += [rmsf]

                    try:
                        r = self.correlation(l_rmsfs,l_bfactors)
                    except:
                        print pdb, mode
                        print l_rmsfs,l_bfactors
                        stop
                    print pdb, mode, r
                    line_out += '%s;' %(r)
                line_out += '\n'
                lines_out += [line_out]

            if lines_out == []:
                os.remove('bfaccorr_%s.txt' %(subdir))
            else:
                fd = open('bfaccorr_%s.txt' %(subdir),'w')
                fd.writelines(lines_out)
                fd.close()

        return


    def parse_recordREMARK470(self, d_header, pdb, lines, i):

        alpha = False

        if lines[i][10:].strip() == 'M RES C SEQI  ATOMS':
            print pdb
            stop

        if lines[i][10:].strip() == 'M RES CSSEQI  ATOMS':

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 470':
                    break

                ## model M
                if lines[j][11:13] == '  ':
                    model = None
                else:
                    model = int(lines[j][11:13])

                ## res_name RES
                res_name = lines[j][15:18].strip()

                ## chain C
                chain = lines[j][19]

                ## res_no SSEQ
##                try:
                res_no = int(lines[j][20:24])
##                except:
##                    res_no = lines[j][20:24].strip()

                ## iCode I
                iCode = lines[j][24]

                ## atoms ATOMS
                atoms = lines[j][25:].split()
                if 'CA' in atoms:
                    if res_name not in self.d_res.keys():
                        stop
                    alpha = True

                if 'REMARK470' not in d_header.keys():
                    d_header['REMARK470'] = {}
                if chain not in d_header['REMARK470'].keys():
                    d_header['REMARK470'][chain] = {}
                if res_no not in d_header['REMARK470'][chain].keys():
                    d_header['REMARK470'][chain][res_no] = {}
                if iCode not in d_header['REMARK470'][chain][res_no].keys():
                    d_header['REMARK470'][chain][res_no][iCode] = atoms
                else:
                    d_header['REMARK470'][chain][res_no][iCode] += [atoms]

        return d_header, alpha


    def parse_recordMODRES(self, d_header, line):

        hetID = line[12:15].strip()
        chain = line[16]
        res_no = int(line[18:22])
        iCode = line[22]
        res_name = line[24:27].strip()
        txt = line[29:80].strip()
        if 'MODRES' not in d_header.keys():
            d_header['MODRES'] = {}
        if chain not in d_header['MODRES'].keys():
            d_header['MODRES'][chain] = {}
        if res_no not in d_header['MODRES'][chain].keys():
            d_header['MODRES'][chain][res_no] = {}
        if iCode not in d_header['MODRES'][chain][res_no].keys():
            d_header['MODRES'][chain][res_no][iCode] = hetID
        elif hetID != d_header['MODRES'][chain][res_no][iCode]:
            print line, s_pdb
            stop

        return d_header


    def correlation(self,l1,l2):

        import math

        sum_x = 0
        sum_y = 0
        sum_xy = 0
        sum_xx = 0
        sum_yy = 0
        if len(l1) != len(l2):
            stop
        n = len(l1)
        for i in range(n):
            sum_x += l1[i]
            sum_y += l2[i]
            sum_xy += l1[i]*l2[i]
            sum_xx += l1[i]*l1[i]
            sum_yy += l2[i]*l2[i]

        r = (n*sum_xy-sum_x*sum_y)/math.sqrt((n*sum_xx-sum_x**2)*(n*sum_yy-sum_y**2))

        return r


    def parse_recordEXPDTA(self,line,d_header,pdb):
        l_methods = line[10:].strip().split(',')
        if (
            len(l_methods) > 1 and
            l_methods[0] not in ['NMR','INFRARED SPECTROSCOPY','SOLUTION SCATTERING',] and
            l_methods[-1] != 'STRUCTURES'
            ):
            print l_methods
            print line
            stop
        if l_methods[0][:3] == 'NMR':
            d_header['EXPDTA'] = 'NMR'
        elif 'X-RAY' in l_methods[0] or l_methods[0] == 'FIBER DIFFRACTION' or l_methods[0] == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
            d_header['EXPDTA'] = 'DIFFRACTION'
        elif l_methods in [['NEUTRON DIFFRACTION'],]:
            d_header['EXPDTA'] = 'DIFFRACTION'
        elif l_methods in [['CRYO-ELECTRON MICROSCOPY'],['ELECTRON DIFFRACTION'],['ELECTRON TOMOGRAPHY'],]:
            d_header['EXPDTA'] = 'EM'
        elif l_methods[0] == 'SOLUTION SCATTERING':
            d_header['EXPDTA'] = 'SOLUTION SCATTERING'
        elif l_methods[0] == 'INFRARED SPECTROSCOPY':
            d_header['EXPDTA'] = 'INFRARED SPECTROSCOPY'
        elif l_methods[0] == 'FLUORESCENCE TRANSFER':
            d_header['EXPDTA'] = 'FLUORESCENCE TRANSFER'
        else:
            print line
            print l_methods
            print pdb
            print d_header['REMARK2']
            stop
        return d_header


    def eigenv_calccomb(self,matrix_hessian, verbose = True):

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

        import math, Numeric

        ## calculate length of mode 7 (equals 1 when using module linearalgebra)
        len7 = math.sqrt(sum(Numeric.array(eigenvectors[6])*Numeric.array(eigenvectors[6])))

        ## loop over modes
        for i in range(7,len(eigenvalues)-1):

            ## calculate length of mode i
            leni = math.sqrt(sum(Numeric.array(eigenvectors[i])*Numeric.array(eigenvectors[i])))

            ## scale length of mode i relative to length of mode 7
            if eigenvalues[6] == 0:
                lenfactor = 1
            else:
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


    def sigmoid(self, x, cutoff, slope=1):
        import math
        y = 1. / ( 1. + math.exp( slope*(x-cutoff) ) )
        return y


    def parse_recordATOM(self, line, d_pdb, d_header, record):

        import Numeric

        atom_no = int(line[6:11])
        atom_name = line[12:16].strip()
        atom_altloc = line[16]
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        bfactor = float(line[60:66])
        coordinate = Numeric.array([x, y, z])

        if record == 'HETATM':
            if not 'MODRES' in d_header.keys():
                return d_pdb
            if not chain in d_header['MODRES'].keys():
                return d_pdb
            if not res_no in d_header['MODRES'][chain].keys():
                return d_pdb
            if not iCode in d_header['MODRES'][chain][res_no].keys():
                return d_pdb
##            if d_header['MODRES'][chain][res_no][iCode] != res_name:
##                print chain, res_no, iCode, res_name
##                print d_header['MODRES'][chain][res_no][iCode]
##                stop

        if not 'chains' in d_pdb.keys():
            d_pdb['chains'] = {}
        if not chain in d_pdb['chains'].keys():
            d_pdb['chains'][chain] = {}
        if not 'residues' in d_pdb['chains'][chain].keys():
            d_pdb['chains'][chain]['residues'] = {}

        ## res_no
        if not res_no in d_pdb['chains'][chain]['residues'].keys():
            d_pdb['chains'][chain]['residues'][res_no] = {}

        ## res_no > d_iCodes
        if not 'd_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
        ## d_iCodes > iCode
        if not iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes']:
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

        ## iCode > res_name
        if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name
        ## temp!!!
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name: ## temp!!!
            return d_pdb
##            print res_name, d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] ## temp!!!
##            print chain, res_no, iCode, atom_altloc ## temp!!!
##            print line
##            stop
        ## check that res_name is correct (e.g. 2fes:L:1)
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name and atom_altloc == ' ': ## 1fh2:A:30 atom_altloc
            ## change the iCode
            iCode_max = max(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
            iCode_max = self.s_alphabet[self.s_alphabet.index(iCode_max)+1]
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_max] = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'].index(iCode)] = iCode_max

            ## d_iCodes > iCode
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

        ## res_no > l_iCodes
        if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
        ## l_iCodes > iCode
        if iCode not in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

        ## iCode > atoms
        if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
        ## atoms > atom_name > coordinate
        if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate,'bfactor':bfactor}

        return d_pdb


    def parse_recordSEQRES(self, line, d_header, pdb):

        chain = line[11]

        if 'SEQRES' not in d_header.keys():
            d_header['SEQRES'] = {}
        if chain not in d_header['SEQRES'].keys():
            d_header['SEQRES'][chain] = {}
        if not 'type' in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['type'] = 'N/A'

        residues = line[19:70].split()

        for i in range(len(residues)):
            residue = residues[i]
            if residue in self.d_res.keys():
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'peptide'
                elif d_header['SEQRES'][chain]['type'] != 'peptide': ## and len(residues) > 0:
                    print d_header['SEQRES'][chain]
                    print residue, residues
                    print pdb
                    if pdb != '1ttt':
                        stop_peptide
                residues[i] = residue
            elif residue in ['C','A','U','G','I','DC','DA','DT','DG','DI','N']: ## N is any 5'-monophosphate nucleotide
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'nucleotide'
                elif d_header['SEQRES'][chain]['type'] != 'nucleotide': ## and len(residues) > 0:
                    print d_header['SEQRES'][chain]
                    print residue, residues
                    stop_nucleotide
                residues[i] = residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'saccharide'
                elif d_header['SEQRES'][chain]['type'] != 'saccharide': ## and len(residues) > 0:
                    print d_header['SEQRES'][chain]
                    print residue, residues
                    stop_saccharide
                    stop
                residues[i] = residue
            else:
                if residue == 'UNK': ## e.g. 1pny.pdb
                    if d_seq['chains'][chain]['type'] == 'N/A':
                        d_seq['chains'][chain]['type'] = 'peptide'
                residues[i] = residue

        if 'seq' not in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['seq'] = []
        d_header['SEQRES'][chain]['seq'] += residues

        return d_header


    def determine_if_modres(self, d_seq, chain, res_no, iCode, res_name):

        modres = False
        if chain in d_seq['MODRES'].keys():
            if res_no in d_seq['MODRES'][chain].keys():
                if iCode in d_seq['MODRES'][chain][res_no].keys():
                    if res_name == d_seq['MODRES'][chain][res_no][iCode]:
                        modres = True
                    else:
                        print chain, res_no, iCode, res_name, d_seq['MODRES'][chain][res_no][iCode]
                        notexpected

        return modres


    def __init__(self):

        import os

        self.d_res = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            'MSE':'M','UNK':'X','ASX':'X','GLX':'X',
            }

        self.l_nucleotides = [
            'A','C','G','U','I', ## ribonucleotides
            'DA','DC','DG','DT','DI', ## deoxyribonucleotides
            'N', ## wild card
            ]
        
        self.d_res3 = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

        self.s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        ## info from the PDB Ligand Depot (searched for cluster in chemical name)
        self.l_clusters = [
            ## iron clusters
            'CFM','CFN','CLF','CLP','CN1','CNB','CNF','F3S','FES','FS1','FS2','FS4','FSF','FSO','FSX','HC0','HC1','HF3','HF5','NFS','SF3','SF4','WCC','XCC', ## 'FS3' deprecated ('F3S' maintained)
            ## copper clusters
            'CUB','CUM','CUN','CUO',
            ## molybdenum "clusters"
            'OMO',
            ## hafnium clusters
            'PHF',
            ## zirconium clusters
            'ZRC',
            ]

        self.l_prosthetic_groups = [
            ## porphyrins (cyclic tetrapyrroles)
            ## Ferrochelatase catalyzes protophorphyrin+Fe(2+) --> protoheme + 2H(+)
            ## iron
            'HEM', ## protoporphyrin IX + Fe(II) (C3 vinyl,C8 vinyl,C18 methyl; tetramethyl,divinyl,dipropionate; *charge 2*)
            'HEC', ## Heme C (protoporphyrin IX; *charge 0*)
            'HEA', ## Heme A (C3 hydroxyfarnesyl, C8 vinyl, C18 formyl)
            'HEO', ## Heme O (C3 hydroxyfarnesyl, C8 vinyl, C18 methyl)
            '2FH', ## 2-phenylheme
            '1FH', ## 12-phenylheme
            'DDH', ## dedivinyl,diacetyl heme (C3,C8)
            'HEV', ## dedimethyl,divinyl heme
            'HDM', ## tetramethyl,divinyl,dipropionate *ester* heme
            'HAS', ## Heme-As (C3 geranylgeranyl, C8 vinyl, C18 formyl)
            'VER', ## octaethylated porphyrin
            'HEB', ## Heme B/C hybrid (tetramethyl,divinyl,dipropionate)
            'DHE', ## Heme D (heme B derivative)
            'HDD', ## cis-heme D hydroxychlorin gamma-spirolactone
            'SRM', ## siroheme (partially reduced iron-porphyrin in e.g. nitrate reductase)
            ## other metal
            'HNI', ## protoporphyrin IX + Ni(II)
            'HES', ## Zn substituted Heme C
            ]
        self.l_coenzymes = [
            'RET', ## vitamin A
##            'TPP', ## vitamin B1            
            'FMN','FAD', ## vitamin B2
            'NAD','NAP', ## vitamin B3
            'COA', ## vitamin B5
            'PLP', ## vitamin B6
            'C2F', ## vitamin B9 (5-methyl THF)
            ]

        ## list of metals also contains metalloids (e.g. Arsen)
        self.l_atoms_metal = [
            'LI','NA','K','CS', ##1a
            'BE','MG','CA', ##2a
            'AL','GA','TL', ##3a
            'PB', ##4a
            'AS', ##5a
            'V', ##3b
            'CR','MO', ##4b
            'MN', ##5b
            'FE', ##6b
            'CO', ##7b
            'NI', ##8b
            'CU', ##9b
            'ZN','CD','HG', ##10b
            ]

        ## info from the PDB Ligand Depot
        ## keys are hetIDs, values are charges (not oxidation states)
        ## hetID:[chemical formula,charge]
        self.d_ions = {
            ## nitrate, ammonium
            'NO3':['N1 O3',-1],'NH4':['H4 N1',+1],
            ## hydroxide
            'OH' :['H1 O1',-1],
            ## phosphate
            '2HP':['O4 P1',-1],'PI' :['H1 O4 P1',-2],'PO4':['O4 P1',-3], ## different oxidation states; 'IPS' deprecated
            ## sulfate
            'SO4':['O4 S1',-2],'SOH':[3,-1],'SUL':[3,-2], ## different oxidation states
            ## carbonate
            'CO3':['C O3', -1],
            ## group1a
            'LI' :['LI1',+1],
            'NA' :['NA1',+1], ## 'NAO','NA2','NA6','NA5','NAW' deprecated
            'K'  :['K1' ,+1], ## 'KO4' deprecated
            'CS' :['CS' ,+1],
            ## group2a
            'BEF':['BE F3',-1],
            'MG' :['MG1',+2], ## 'MO3','MO1','MO2','MO4','MO5','MO6' deprecated
            'CA' :['CA1',+2],'OC1':['CA1',+2], ## 'OC3','OC5','OC6','OC7','OC8','OC2','OC4' deprecated
            'SR' :['SR' ,+2],
            'BA' :['BA' ,+2],
            ## group3a
            'AL' :['AL1',+3],'ALF' :['AL F4',-1],
            'GA' :['GA1',+3],'GA1' :['GA1',+2],
            'TL' :['TL1',+1],
            ## group4a
            'ARS':['AS1', 0],'ART':['O4 AS1',-3],'AST':-3,'TAS':['H3 O3 AS1', 0],'CAC':['C2 H6 AS O2', -1], ## different compounds
            'PB' :['PB' ,+2],
            ## group6a
            'SE' :['SE1', 0],'SE4':['O4 SE1',-2], ## different compounds
            ## group7a
            'CL' :['CL1',-1],
            'BR' :['BR1',-1],
            'IOD':['I1' ,-1],
            ## group8a
            'KR' :['KR1', 0],
            ## group3b
            'V'  :+3,'VO4':['V1' ,-3], ## different oxidation states
            ## group4b
            'CR' :['CR1',+3],
            'MO' :['MO1', 0],'4MO':['MO1', 0],'2MO':['MO O2',-2], ## different compounds and different oxidation states
            ## group5b
            'MN' :['MN1',+2],'MN3':['MN1',+3], ## different oxidation states; 'MW1','MW2','MW3','O4M','MN5','MN6' deprecated
            ## group6b
            'FE2':['FE1',+2],'FE' :['FE1',+3], ## different oxidation states; 'OF1','OF3','2OF' deprecated
            ## group7b
            'CO' :['CO1',+2],'3CO':['CO1',+3], ## different oxidation states; 'CO5','OCL','OCO','OCN','OCM' deprecated
            ## group8b
            'NI' :['NI1',+2],'3NI':['NI1',+3], ## different oxidation states; 'NI1','NI2','NI3','NIK' deprecated
            ## group9b
            'CU1':['CU1',+1],'CU' :['CU1',+2], ## different oxidation states; '1CU' deprecated
            ## group10b
            'ZN' :['ZN1',+2], ## 'ZN2','ZO3','ZN3','ZNO' deprecated
            'CD' :['CD1',+2],
            'HG' :['HG1',+2],
            ## Lanthanides
            'TB' :['TB1',+3],
            'YB' :['YB1',+3],
            }

        self.l_cofactors = self.l_clusters+self.l_prosthetic_groups+self.l_coenzymes+self.d_ions.keys()

        ## wikipedia buffer solution, Good's buffers
        self.l_solutes = [
            ## IPA,FMT,GOL,EEE,EDO
            ##
            ## ethylene glycol, protein precipitation
            'EDO',
            ## acetic acid
            'ACY',
            ## water
            'HOH','H2O','DOD','D2O',
            ## methanol
            'MOH',
            ## di-thio-threitol (reducing agent)
            'DTT', 
            ## beta-mercapto-ethanol (reducing agent)
            'BME',
            ## bis-tris methane (buffering agent)
            'BTB',

            ## TAPS (buffering agent)
            'T3A',
            ## Bicine (buffering agent)
            'BCN',
            ## Tris (buffering agent)
            'TRS',
            ## HEPES (buffering agent)
            'EPE',
            ## TES (buffering agent)
            'NES',
            ## MOPS
            'MPO',
            ## PIPES
            'PIN',
            ## Cacodylate
            'CAC',
            ## MES
            'MES',
            ## Acetate
            'ACT',
            ## unknown atom or ion
            'UNX',
            ]

        ## saccharides returned from a search of the ligand depot for disaccharides and monosaccharides
        self.d_saccharides = {
            ##
            ## monosaccharides, aldehydes
            ##
            ## pyranoses, hexoses
            'GLC':{'stereo':'GLC','derivate':['GLC']}, ## (alpha)-D-Glucose
##            'AGC':{'stereo':'GLC','derivate':['GLC']}, ## alpha-D-Glc (deprecated)
            'BGC':{'stereo':'GLC','derivate':['GLC']}, ## beta-D-Glc
            'GAL':{'stereo':'GAL','derivate':['GAL']}, ## (beta)-D-Galactose
            'GLA':{'stereo':'GAL','derivate':['GAL']}, ## alpha-D-Gal
            'GLB':{'stereo':'GAL','derivate':['GAL']}, ## beta-D-Gal
            'FUC':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, alpha-L-Fucose
            'FUL':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, beta-L-Fucose
            'MAN':{'stereo':'MAN','derivate':['MAN']}, ## alpha-D-Mannose
            'BMA':{'stereo':'MAN','derivate':['MAN']}, ## beta-D-Mannose
            'ARA':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabinose
            'ARB':{'stereo':'ARA','derivate':['ARA']}, ## beta-L-Arabinose
            ## pyranoses, pentoses
            'XYS':{'stereo':'XYS','derivate':['XYS']}, ## (alpha)-D-Xylose
            'XYP':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-Xylose
            'XYP':{'stereo':'XYS','derivate':['XYS']}, ## beta-L-Xylose
            ## furanoses, hexoses
            'AHR':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabino*furano*se
            ## furanoses, pentoses
            'XYZ':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-xylo*furano*se
            ## phosphorylated aldohexopyranoses
            'G1P':{'stereo':'G1P','derivate':['GLC']}, ## alpha-D-Glc-1P
            'G6P':{'stereo':'G6P','derivate':['GLC']}, ## alpha-D-Glc-6P
            'BG6':{'stereo':'G6P','derivate':['GLC']}, ## beta-D-Glc-6P
            'G6Q':{'stereo':'G6P','derivate':['GLC']}, ## Glc-6P
            'BGP':{'stereo':'BGP','derivate':['GAL']}, ## beta-Gal-6P
            'M1P':{'stereo':'M1P','derivate':['MAN']}, ## alpha-D-Man-1P
            'M6P':{'stereo':'M6P','derivate':['MAN']}, ## alpha-D-Man-6P
            ## phosphorylated aldopentofuranoses
            'ABF':{'stereo':'ABF','derivate':['ARA']}, ## beta-D-Arabino*furano*se-5-phosphate
            ## methylated aldohexopyranoses
            'MMA':{'stereo':'MMA','derivate':['MAN']}, ## O1-methyl-mannose
            ## aminated aldohexopyranoses amine
            'AGL':{'stereo':'AGL','derivate':['GLC']}, ## 4,6-dideoxy-4-amino-alpha-D-Glucose
            ## deoxygenated aldohexopyranoses
            'G6D':{'stereo':'G6D','derivate':['GLC']}, ## 6-deoxy-alpha-D-Glucose
            ## oxygenated aldohexopyranoses
            'KBG':{'stereo':'G6D','derivate':['GLC']}, ## 2-keto-beta-D-Glucose
            ## acetylated aldohexopyranose amines
            'NAG':{'stereo':'NAG','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine
            'NBG':{'stereo':'NAG','derivate':['GLC']}, ## 1-N-Acetyl-beta-D-Glucosamine
            'NDG':{'stereo':'NAG','derivate':['GLC']}, ## 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
            '16G':{'stereo':'16G','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine-6P
            'NGA':{'stereo':'NGA','derivate':['GLC']}, ## (2)-N-Acetyl-D-Galactosamine
            ##
            ## monosaccharides, ketones
            ##
            ## furanoses, hexoses
            'FRU':{'stereo':'FRU','derivate':['FRU']}, ## Fructose
            'F6P':{'stereo':'F6P','derivate':['FRU']}, ## Fru-6P
            ##
            ## dissacharides
            ##
            'SUC':{'stereo':'SUC','derivate':['GLC','FRU']}, ## GLC-a12-FRC, Sucrose
            'LAT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, alpha-Lactose
            'LBT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, beta-Lactose
            'MAL':{'stereo':'MAL','derivate':['GLC','GLC']}, ## GLC-a14-GLC, Maltose
            'TRE':{'stereo':'TRE','derivate':['GLC','GLC']}, ## GLC-a11a-GLC, Trehalose
            'CBI':{'stereo':'CBI','derivate':['GLC','GLC']}, ## GLC-b14-GLC, Cellobiose
            ##
            ## polysaccharides
            ##
            'MTT':{'stereo':'MTT','derivate':['MAL','MAL']}, ## MAL-b14-MAL, Maltotetraose
            ##
            ## linear saccharides (neither furanoses nor pyranoses...) and their derivatives...
            ##
            'SOR':{'stereo':'SOR','derivate':['GLC']}, ## Sorbitol/Glucitol (reduced glucose)
            'GLO':{'stereo':'GLO','derivate':['GLC']}, ## linear glucose
            'XLS':{'stereo':'XLS','derivate':['XYL']}, ## linear xylose
            'A5P':{'stereo':'ABF','derivate':['ARA']}, ## Arabinose-5-phosphate
            ##
            ## conduritol (1,2,3,4-cyclohexenetetrol) derivatives
            ##
            'HMC':{'stereo':'HMC'}, ## 5-hydroxymethyl-chonduritol
            'ACI':{'stereo':'ACI'}, ## 1-amino-5-hydroxymethyl-chonduritol
            ## Sialic acid (N-Acetylneuraminic acid, Neu5Ac, NANA)
            'SIA':{'stereo':'SIA','derivate':['SIA']}, ## (alpha)-sialic acid
            'SLB':{'stereo':'SIA','derivate':['SIA']}, ## beta-sialic acid
            }

        self.d_stereoisomers = {
            }

        self.l_expdta = [
            'X-RAY',
            'ELECTRON DIFFRACTION','NEUTRON DIFFRACTION',
            'NMR',
            'INFRARED SPECTROSCOPY',
            'CRYO-ELECTRON MICROSCOPY',
            'ELECTRON TOMOGRAPHY', ## e.g. 1o1a.pdb
            'SOLUTION SCATTERING', ## e.g. 1e07.pdb
            'FLUORESCENCE TRANSFER', ## e.g. 1rmn.pdb
            ]

        self.l_columns_html = [
            'gif1','gif2','pdb1', 'pdb2', 'bm1', 'bm2',
            'rmsd', 'mutations', 'chains', 'residues', 'coordinates',
            'pH1', 'pH2', 'T1', 'T2', 'res1', 'res2','spacegroup1', 'spacegroup2',
            'REMARK465', 'REMARK470','transformations',
            'title1','title2','hetIDs1', 'hetIDs2'
            ]

        self.maxrmsd = 2.75 ## 2fsy.pdb,2ft1.pdb

        self.minres = 5.0

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25
        self.max_mutations = 10

        self.path_pdb = '/oxygenase_local/data/pdb/'
        self.path_cwd = os.getcwd()
        
        return

if __name__ == '__main__':
    instance_quakes = quakes()
    instance_quakes.main()


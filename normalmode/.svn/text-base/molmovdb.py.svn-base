#!/software/bin/python
#
#$Id: molmovdb.py 267 2007-11-01 13:56:25Z tc $
#
#Tommy Carstensen, University College Dublin, 2007

import os, sys, Numeric, LinearAlgebra, math, urllib2
sys.path.append('/home/people/tc/svn/Protool/')
import geometry
instance_geometry = geometry.geometry()
from sets import *
set = Set

class molmovdb:

    def main(self):

        l_inputs = [
##            ['1ra2','1ra1','A','A',],
##            ['6ldh','1ldm','A','A',],
##            ['1eqm','1hka','A','A',],
##            ['1ex6','1ex7','A','A',],
##            ['4ake','1ake','A','A',],
##            ['1ggg','1wdn','A','A',],
##            ['2lao','1lst','A','A',],
##            ['1ypi','2ypi','A','A',],
##            ['1urp','2dri','A','A',],
##            ['1d9v','1mrp','A','A',],
##            ['2lzm','189l','A','A',],
##            ['2lzt','1sfg','A','A',],
##            ['1atp','1j3h','E','A',],
            ['2lzm','150l','A','D',],
            ]
        l_max_overlaps = []
        for input in l_inputs:
            pdb1 = input[0]
            pdb2 = input[1]
            chain1 = input[2]
            chain2 = input[3]
               
            max_overlap = self.calculate_overlap(pdb1,pdb2,chain1,chain2,)
            print input, max_overlap

            l_max_overlaps += [max_overlap]

        print l_max_overlaps
        print [[8, 0.82049528114527259]]
        stop

        self.loop_molmovdb()

        self.gnuplot()

        return


    def loop_molmovdb(self):

        errorpdbs = [
            '1p0t','1otz','1utf','1utv','1jgo','1jgp','1jgq','1gix','1giy','1pns','1pnu','1pnx','1pny','1s1h','1s1i','1voq','1vor','1vos','1vou','1vov','1vow','1vox','1voy','1voz','1vp0','1vs5','1vs6','1vs7','1vs8','1vsa','2ow8','2avy','2aw4','2aw7','2awb','2b64','2b66','2b9m','2b9n','2b9o','2b9p','2i2p','2i2t','2i2u','2i2v','2j00','2j01','2j02','2j03','2qou','2qov','2qow','2qox','2qoy','2qoz','2qp0','2qp1','2v46','2v47','2v48','2v49','2hgi','2hgj','2hgp','2hgq','2uv9','2uva','2uvb','2uvc','2bod','1run','1n2c','1fzg','2fa7','1ljy','2c3w','1j3y','2axm','3xis','2bhy','2c03','3car',
            '1gct','1hao','2iff','1x9s','1x9w','1xkf','1wu1','1tk5','1s10','1skw','1sal','1ry1','1pts','1t7b','1t7e','1tud','2ht8','1yra','1yrb','1t41','2qwa','2qwb','2qwc','2qwd','2qwe','2qwf','2qwg','2qwh','2qwi','2qwj','2qwk','1pbh','2pbh','2hk3','1jrt','1jrs','1g2w','1mbq','2tbs','1xwq','1nhs','1a14','1cu1',
            '1vjh','1zii','1zij','1xri','1o18','1tmq','2ja5','2ja6','2ja7','2ja8','1nc2','1u6p','1way','2bvr','2bxt','5gds','1bei','1c4e','1mqx','1mqy','1mqz','2h1o','2v63','2v67','1fng','1fne','1abe','1abf','1ow6','1ow7','2alu','1nkx','1azx','2gnm','2gn7','1en2','1enm','2cfg','1d8s','1i3q','1i50','1iw7','1abj','1uyx','1mpm','1k7d','1fmu','4atj','1e56','1sbd','1sbe','1r9m','1ofz','2nn8','1chh','1ccr','1chi','1chj','1cie','1cif','1cig','1chi','1crg','1crh','1cri','1crj','1csu','1csv','1csw','1csx','1cty','1ctz','1fhb','1rap','1raq','1ycc','1yea','1yeb','1ytc','2cln','2ycc','2exj','1ibv','1mtb','2fgu','1odw','1qs7','1zbg','1fff','1fg8','1k2b','1vwb','1ogt','1aym','2bd2','2bd4','2bd7','2pyp','1ppv','1ca8','1sbh','1sbi','1sbn','1yja','1yjb','1xpk','5abp','1axz','1vba','1al2','1ar6','1ar9','1asj','2gp1','1do1','1do3','1do4','1do7','1fh2','2fi5','1eis','1uwn','1uwp','1k55','1k56','2ftm','1vwi','1vwj','1vwk','1vwl','1vwq','1vwr','2h8p','2h5g','1ft8','2agx',
            ]

        d_sprsde = {'1hfv':'2j5x','1e7x':'2c4q','1dzs':'2iz9','1eqo':'1q0n'}

        url = urllib2.urlopen('http://www.molmovdb.org/cgi-bin/browse.cgi')
        lines = url.readlines()
        urls = []
        for line in lines:
            if (
                line[:56] == '[<A HREF="http://www.molmovdb.org/cgi-bin/motion.cgi?ID=' and
                line[-14:] == '">motion</A>]\n'
                ):
                index1 = line.index('"')+1
                index2 = line.rindex('"')
                url = line[index1:index2]
                urls += [url]

        l_data = []

        urls.sort()
        for j in range(len(urls)):
            url = urls[j]
            suffix = url[url.rindex('=')+1:]
 
            print
            print j, len(urls), url

            ##
            ## skip if exists
            ##
            if os.path.isfile('%s.txt' %(suffix)):
                continue

            ##
            ## parse url
            ##
            try:
                url = urllib2.urlopen(url)
            except:
                continue
            lines = url.readlines()
            pdb1 = pdb2 = ''
            for i in range(len(lines)):
                line = lines[i]
                if '<B>Classification</B>' in line:
                    index1 = lines[i+3].index('>')+1
                    index2 = lines[i+3].rindex('<')
                    classification = lines[i+3][index1+1:index2-1]
                if '>Best representative<' in line:
                    if '</A>' not in lines[i+11] and 'upload' in lines[i+11]:
                        break
                    if '</A>' not in lines[i+12] and 'upload' in lines[i+12]:
                        break
                    index2 = lines[i+11].index('</A>')
                    index1 = lines[i+11][:index2].rindex('>')+1
                    pdb1 = lines[i+11][index1:index2].lower()
                    index1 = lines[i+11].index('[')+1
                    index2 = lines[i+11].rindex(']')
                    chain1 = lines[i+11][index1+1:index2-1]
                    index2 = lines[i+12].index('</A>')
                    index1 = lines[i+12][:index2].rindex('>')+1
                    pdb2 = lines[i+12][index1:index2].lower()
                    index1 = lines[i+12].index('[')+1
                    index2 = lines[i+12].rindex(']')
                    chain2 = lines[i+12][index1+1:index2-1]
                    if ((pdb1 == '' and pdb2 == '') or (pdb1 == 'new1' or pdb2 == 'new2')) and chain1 == '' and chain2 == '':
                        break
                    if len(pdb1) != 4 or len(pdb2) != 4:
                        print pdb1, pdb2
                        stop
            if ((pdb1 == '' and pdb2 == '') or (pdb1 == 'new1' or pdb2 == 'new2')):
                continue
            if pdb1 == pdb2:
                continue
            if chain1 == '':
                chain1 = 'A'
            if chain2 == '':
                chain2 = 'A'
            if pdb1 in d_sprsde.keys():
                pdb1 = d_sprsde[pdb1]
            if pdb2 in d_sprsde.keys():
                pdb2 = d_sprsde[pdb2]
            l_data += [[pdb1,pdb2,chain1,chain2,classification]]

            if pdb1 in errorpdbs:
                print pdb1
                continue
            if pdb2 in errorpdbs:
                print pdb2
                continue
            if classification[0] == 'N':
                continue
            print pdb1, pdb2, chain1, chain2, classification

            ##
            ## skip if error
            ##
            if (
                ## fragments
                (pdb1 == '1aon' and pdb2 == '1egs') or
                ## MODRES
                (pdb1 == '1jmj' and pdb2 == '1jmo') or ## TYS/TYR
                ## mutations (SEQADV HIS TAG)
                (pdb1 == '1q12' and pdb2 == '1q1b') or
                (pdb1 == '1d6m' and pdb2 == '1i7d') or
                ## insertions/deletions (e.g. different organisms, synthetic constructs)
                (pdb1 == '4mdh' and pdb2 == '1bmd') or
                (pdb1 == '1axn' and pdb2 == '2ran') or
                (pdb1 == '1crl' and pdb2 == '1thg') or
                (pdb1 == '1hnf' and pdb2 == '1hng') or
                (pdb1 == '1gtm' and pdb2 == '1hrd') or
                (pdb1 == '1avk' and pdb2 == '1a2v') or
                (pdb1 == '1etu' and pdb2 == '1eft') or
                (pdb1 == '1gdh' and pdb2 == '1psd') or
                (pdb1 == '1b1a' and pdb2 == '1be1') or
                (pdb1 == '1lkf' and pdb2 == '7ahl') or
                (pdb1 == '1n8y' and pdb2 == '1n8z') or
                (pdb1 == '1mcp' and pdb2 == '1nca') or
                (pdb1 == '1ejd' and pdb2 == '1a2n') or
                (pdb1 == '13pk' and pdb2 == '1php') or
                (pdb1 == '1jlu' and pdb2 == '1cmk') or
                (pdb1 == '1k20' and pdb2 == '1k23') or
                (pdb1 == '1ab3' and pdb2 == '1a32') or
                (pdb1 == '6tim' and pdb2 == '1tre') or
                (pdb1 == '1qaw' and pdb2 == '1wap') or
                (pdb1 == '1sto' and pdb2 == '1oro') or ## few mutations
                ## terminal extension (multiple residues)
                (pdb1 == '1byu' and pdb2 == '1rrp') or
                (pdb1 == '1bnc' and pdb2 == '1dv2') or
                (pdb1 == '1grn' and pdb2 == '1cf4') or
                (pdb1 == '2hmi' and pdb2 == '3hvt') or
                (pdb1 == '1dv7' and pdb2 == '1dvj') or
                (pdb1 == '1tde' and pdb2 == '1f6m') or
                (pdb1 == '1fgu' and pdb2 == '1jmc') or
                (pdb1 == '7api' and pdb2 == '1psi') or
                (pdb1 == '1ffh' and pdb2 == '1ng1') or
                (pdb1 == '1ik9' and pdb2 == '1fu1') or ## also ENGINEERED
                ## terminal extension (1 N-terminal residue)
                (pdb1 == '1erk' and pdb2 == '2erk') or
                (pdb1 == '1fbt' and pdb2 == '1tip') or
                ## mutations (no explanation)
                (pdb1 == '1eps' and pdb2 == '1g6s') or
                (pdb1 == '1aon' and pdb2 == '1oel') or
                (pdb1 == '1tgl' and pdb2 == '4tgl') or
                ## mutations (SEQADV ENGINEERED, phosphorylation)
                (pdb1 == '1d5w' and pdb2 == '1dcm') or
                ## mutations (SEQADV ENGINEERED)
                (pdb1 == '1jmw' and pdb2 == '1was') or
                (pdb1 == '4crx' and pdb2 == '1crx') or
                (pdb1 == '1buy' and pdb2 == '1eer') or
                (pdb1 == '1i69' and pdb2 == '1i6a') or
                (pdb1 == '1quk' and pdb2 == '1oib') or
                (pdb1 == '2pol' and pdb2 == '1jqj') or
                (pdb1 == '1ih7' and pdb2 == '1ig9') or
                ## mutations (SEQADV CONFLICT)
                (pdb1 == '1dpp' and pdb2 == '1dpe') or ## ecpdpbp
                (pdb1 == '1clb' and pdb2 == '4icb') or
                (pdb1 == '4dfr' and pdb2 == '5dfr') or
                (pdb1 == '3tms' and pdb2 == '2tsc') or
                (pdb1 == '1beb' and pdb2 == '1b0o') or
                (pdb1 == '1dkx' and pdb2 == '1dky') or
                (pdb1 == '1i2d' and pdb2 == '1m8p') or
                (pdb1 == '1lfg' and pdb2 == '1lfh') or
                (pdb1 == '1thv' and pdb2 == '1thi') or
                ## mutations (SEQADV CONFLICT,ENGINEERED)
                (pdb1 == '2efg' and pdb2 == '1fnm') or ## efg
                (pdb1 == '1k89' and pdb2 == '1hrd') or
                (pdb1 == '1gpw' and pdb2 == '1thf') or
                (pdb1 == '1ipd' and pdb2 == '1osj')
                ):
                continue

            ##
            ## create file
            ##
            fd = open('%s.txt' %(suffix),'w')
            fd.write('')
            fd.close()

            max_overlap = self.calculate_overlap(pdb1,pdb2,chain1,chain2,)
            print pdb1,pdb2,max_overlap

            ##
            ## write overlaps to file
            ##
            fd = open('%s.txt' %(suffix),'w')
            fd.writelines(lines)
            fd.close()

        return


    def gnuplot(self):

        ##
        ## gnuplot
        ##
        l_xtics = [
            'D-h-2',
            'F-h-2',
            'D-s-2',
            'F-s-2',
            'D-?-2',
            'F-?-2',
            'D-n-2',
            'F-n-2',
            'D-f-2',
            'S-a-2',
            'S-n-2',
            'C----',
            ]

        ##
        ## gnuplot.data
        ##
        d_data = {}
        gnuplotdata = []
        files = os.listdir(os.getcwd())
        for file in files:
            if file[-4:] != '.txt':
                continue
            ## classification
            url = urllib2.urlopen('http://www.molmovdb.org/cgi-bin/motion.cgi?ID=%s' %(file[:-4]))
            lines = url.readlines()
            for i in range(len(lines)):
                line = lines[i]
                if '<B>Classification</B>' in line:
                    index1 = lines[i+3].index('>')+1
                    index2 = lines[i+3].rindex('<')
                    classification = lines[i+3][index1+1:index2-1]
            ## max overlap
            fd = open(file,'r')
            lines = fd.readlines()
            fd.close()
            overlaps = []
            for line in lines:
                overlap = line.split(';')[1]
                overlaps += [overlap]
            overlaps = overlaps[:3] ## temp!!!
            max_overlap = float(max(overlaps))
            mode = overlaps.index(max(overlaps))+7
            if mode > 12:
                print file, mode, classification, max(overlaps)
            ## append line
            index_classification = l_xtics.index(classification)
            if not index_classification in d_data.keys():
                d_data[index_classification] = []
            d_data[index_classification] += [max_overlap]
            ## gnuplotdata
            gnuplotdata += ['%f %f\n' %(float(index_classification),max_overlap)]
            if classification == 'D-h-2':
                if max_overlap < 0.3:
                    print 'less', file
                if max_overlap > 0.8:
                    print 'more', file

        ## write gnuplot.data
        fd = open('gnuplot.data','w')
        fd.writelines(gnuplotdata)
        fd.close()

        ##
        ## gnuplot.statistics
        ##
        ## average and standard error
        gnuplot_statistics = []
        for index_classification in d_data.keys():
            n = len(d_data[index_classification])
            if n == 1:
                continue
            sumx = 0
            sumxx = 0
            for x in d_data[index_classification]:
                sumx += x
                sumxx += x**2
            average = sum(d_data[index_classification])/len(d_data[index_classification])
            average = sumx/n
            SS = sumxx-(sumx**2)/n
            MSE = SS / (n-1)
            SE = math.sqrt(MSE/n)
            gnuplot_statistics += ['%f %f %f\n' %(float(index_classification),average,SE)]

        ## write gnuplot.statistics
        fd = open('gnuplot.statistics','w')
        fd.writelines(gnuplot_statistics)
        fd.close()

        ##
        ## gnuplot.settings
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
            'plot [-1:%i][0:1] "gnuplot.data" lt 0 ps 2 pt 2' %(len(l_xtics)+1),
            ', "gnuplot.statistics" lt 1 lc 0 ps 0 pt 0 w errorb\n',
            ]

        ## write gnuplot.settings
        fd = open('gnuplot.settings','w')
        fd.writelines(gnuplotsettings)
        fd.close()

        os.system('/software/bin/gnuplot gnuplot.settings')

        return


    def calculate_overlap(self,pdb1,pdb2,chain1,chain2,):

        d_seq = {}
        d_chains_intrapdb_sequence_identical = {}
        for s_pdb in [pdb1,pdb2]:
            ## read lines
            fd = open('%s%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
            lines = fd.readlines()
            fd.close()
            ## parse data prior to the coordinate section
            d_noncoordinates = self.parse_pdbheader(lines, s_pdb)
            d_seq[s_pdb] = d_noncoordinates

        ##
        ## incorrect chain given by molmovdb.org or chain refers to nonremediated file
        ##
        if chain1 not in d_seq[pdb1]['chains'].keys():
            chain1 = min(d_seq[pdb1]['chains'].keys())
            if chain1 == ' ':
                stop1
        if chain2 not in d_seq[pdb2]['chains'].keys():
            chain2 = min(d_seq[pdb2]['chains'].keys())
            if chain2 == ' ':
                stop2

##            ##
##            ## not a protein chain
##            ##
##            if d_seq[pdb1]['chains'][chain1]['type'] != 'peptide':
##                continue
##            if d_seq[pdb2]['chains'][chain2]['type'] != 'peptide':
##                continue

        ## parse coordinates
        d_pdb = {}
        d_pdb = self.parse_coordinates(pdb1, d_pdb, d_seq[pdb1], chain1)
        d_pdb = self.parse_coordinates(pdb2, d_pdb, d_seq[pdb2], chain2)
        l_equivalent_chains = [[chain1],[chain2]]

        seq1, d_resnos1 = self.ATOM2seq(d_pdb[pdb1], chain1, d_seq[pdb1],)
        seq2, d_resnos2 = self.ATOM2seq(d_pdb[pdb2], chain2, d_seq[pdb2],)

        ##
        ## align residues for the specific combinations of chain1 and chain2
        ##
        d_seqalnparam = self.align_residues(chain1, chain2, d_pdb, d_seq, pdb1, pdb2)

        ##
        ## parse coordinates for the specific combinations of chain1 and chain2
        ##
        d_res_nos1 = d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos1']
        d_res_nos2 = d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos2']
        l1 = d_seqalnparam[chain1[0]][chain2[0]]['l1']
        l2 = d_seqalnparam[chain1[0]][chain2[0]]['l2']
        ATOMseq1 = d_seqalnparam[chain1[0]][chain2[0]]['seqATOM1']
        ATOMseq2 = d_seqalnparam[chain1[0]][chain2[0]]['seqATOM2']

        ##
        ## check if sequences identical
        ##
        if ATOMseq1 != ATOMseq2:
            print seq1
            print ATOMseq1
            print seq2
            print ATOMseq2
            print len(seq1), len(ATOMseq1), len(seq2), len(ATOMseq2)
            if ATOMseq1 != seq1:
                terminaldiff1
            if ATOMseq2 != seq2:
                terminaldiff2
            mutations = 0
            if len(ATOMseq1) != len(ATOMseq2):
                sequencesofdifferentlength
            for i in range(len(ATOMseq1)):
                if ATOMseq1[i] != ATOMseq2[i]:
                    mutations += 1
                    print i+1, ATOMseq1[i]
                    print i+1, ATOMseq2[i]
            if mutations > 0:
                print 'mutations', mutations
##                sequencesdifferent

        coords1, coords2 = self.ATOMrecords2coordinates(
             d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2,
             l1, l2, len(ATOMseq1),
             )

        if len(coords1) != len(coords2):
            stop

        rmsd = instance_geometry.superpose(coords1,coords2)
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter

        cutoff_distance = 10.
        verbose = True
        jobid = pdb1
        matrix_hessian = self.hessian_calculation(coords1, float(cutoff_distance), verbose)
        eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_combined_nonperturbed = self.eigenv_calccomb(matrix_hessian, jobid, verbose)

        difference_vector = []
        for i in range(len(coords1)):
            coord1 = coords1[i]
            coord2 = Numeric.matrixmultiply(coords2[i]-tv1,rm)+tv2
            difference_vector += [
                coord1[0]-coord2[0],
                coord1[1]-coord2[1],
                coord1[2]-coord2[2],
                ]

        lines = []
        max_overlap = [0,0,]
        for mode in range(6,len(eigenvectors_nonperturbed)):
            overlap = self.cosangle(difference_vector, eigenvectors_nonperturbed[mode])
            print mode+1, overlap
            lines += ['%i;%f\n' %(mode+1,overlap)]
            if overlap > max_overlap[1]:
                max_overlap = [mode+1,overlap]
        print max_overlap

        return max_overlap


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
        cosang = abs(numerator / denominator)
        return cosang


    def hessian_calculation(self, coordinates, cutoff, verbose = False):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

        N = (len(coordinates))

        if verbose == True:
            print 'calculating the %ix%i hessian/second-order derivate Kirchoff/Laplacian matrix' %(3*N,3*N)
        
        import math, Numeric, time

        cutoff_sq = cutoff**2

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


    def eigenv_calccomb(self, matrix_hessian, jobid, verbose):

        '''Calculates eigenvectors and eigenvalues of a matrix.'''

        if verbose == True:
            print 'calculating eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'
        
        import LinearAlgebra
        import numpy
        ## diagonalize hessian matrix
        eigen_tuple = numpy.linalg.eig(matrix_hessian)
        ## parse eigenvalues and eigenvectors
        eigenvalues = list(eigen_tuple[0])
        eigenvectors = list(eigen_tuple[1])
        eigenvectors = numpy.transpose(eigenvectors)
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
##            fd = open('%s_eigenvectors.txt' %(jobid),'w')
##            fd.writelines(lines)
##            fd.close()

        import math, Numeric

        ## calculate length of mode 7 (equals 1 when using module linearalgebra)
        len7 = math.sqrt(sum(numpy.array(eigenvectors[6])*numpy.array(eigenvectors[6])))
        print len7
        stop

        ## loop over modes
        for i in range(7,len(eigenvalues)-1):

            ## calculate length of mode i
            leni = math.sqrt(sum(numpy.array(eigenvectors[i])*numpy.array(eigenvectors[i])))

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


    def align_residues(self, chains1, chains2, d_pdb, d_seq, pdb1, pdb2):

        d_seqalnparam = {}

        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            if chain1[0] in d_seqalnparam.keys():
                if chain2[0] in d_seqalnparam[chain1[0]].keys():
                    continue

            d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = self.alignATOMseq(d_pdb, d_seq, pdb1, pdb2, chain1, chain2)
            if not chain1[0] in d_seqalnparam.keys():
                d_seqalnparam[chain1[0]] = {}
            if not chain2[0] in d_seqalnparam[chain1[0]].keys():
                d_seqalnparam[chain1[0]][chain2[0]] = {}
            d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos1'] = d_res_nos1
            d_seqalnparam[chain1[0]][chain2[0]]['d_res_nos2'] = d_res_nos2
            d_seqalnparam[chain1[0]][chain2[0]]['l1'] = l1
            d_seqalnparam[chain1[0]][chain2[0]]['l2'] = l2
            d_seqalnparam[chain1[0]][chain2[0]]['seqATOM1'] = ATOMseq1
            d_seqalnparam[chain1[0]][chain2[0]]['seqATOM2'] = ATOMseq2

        return d_seqalnparam


    def identify_similar_chains_from_sequence_inter(
        self, d_seq, pdb1, pdb2,
        d_chains_intrapdb_sequence_identical,
        bmchains1, bmchains2,
        d_biomolecules1, d_biomolecules2,
        ):

##        print 'identifying similar chains'

        import sys
        sys.path.append('/home/people/tc/svn/trunk/EAT_DB/')
        sys.path.append('../../EAT_DB/')
        import sequence_alignment

        ## return the following dictionary structure
        ## interpdb: {repchain:{seqsimchain:{l1,l2}}}

        ## identify repchains
        repchains1 = []
        for repchain1 in d_chains_intrapdb_sequence_identical[pdb1].keys():
            chains1seqid = [repchain1]+d_chains_intrapdb_sequence_identical[pdb1][repchain1]
            for bmchain1 in bmchains1:
                if bmchain1 in chains1seqid:
                    repchains1 += [repchain1]
                    break
        repchains2 = []
        for repchain2 in d_chains_intrapdb_sequence_identical[pdb2].keys():
            chains2seqid = [repchain2]+d_chains_intrapdb_sequence_identical[pdb2][repchain2]
            for bmchain2 in bmchains2:
                if bmchain2 in chains2seqid:
                    repchains2 += [repchain2]
                    break

        ## set d_chains_interpdb_sequence_similar
        d_chains_interpdb_sequence_similar = {}

        ## only do sequential alignment for representative chains to save time
        for chain1 in repchains1:

            if d_seq[pdb1]['chains'][chain1]['type'] != 'peptide':
                continue

            seq1 = d_seq[pdb1]['chains'][chain1]['seq']

            if len(seq1) < self.min_len_chain:
                continue

            ## only do sequential alignment for representative chains to save time
            for chain2 in repchains2:

                if d_seq[pdb2]['chains'][chain2]['type'] != 'peptide':
                    continue

                seq2 = d_seq[pdb2]['chains'][chain2]['seq']

                if len(seq2) < self.min_len_chain:
                    continue

                if abs(len(seq1)-len(seq2)) > self.max_len_chain_difference:
                    continue

                ## fast sequence comparison
                if len(seq1) == len(seq2):

                    l1 = 0
                    l2 = 0
                    r1 = 0
                    r2 = 0
                    s1 = seq1
                    s2 = seq2

                    ## sequence identical
                    if seq1 == seq2:

                        n_chainmutations = 0
                        l_chainmutations = []

                    ## point mutation(s)
                    else:

                        s1 = seq1
                        s2 = seq2
                        n_chainmutations, l_chainmutations = self.point_mutations(seq1, seq2)

                else:

                    ## 1st slow sequence comparison (SEQRESseq)

                    print pdb1, pdb2, chain1, chain2, 'begin seq aln of chains of len %s and %s' %(len(seq1),len(seq2))
                    instance = sequence_alignment.NW(seq1,seq2)
##                    print 'aligning chain %s,%s of %s,%s' %(chain1,chain2,pdb1,pdb2)
                    s1,s2 = instance.Align(verbose=False)[:2]
                    print 'end seq aln'

                    l1 = len(s1)-len(s1.lstrip('-'))
                    l2 = len(s2)-len(s2.lstrip('-'))
                    l = max(l1,l2)
                    r1 = len(s1)-len(s1.rstrip('-'))
                    r2 = len(s2)-len(s2.rstrip('-'))
                    r = max(r1,r2)
    ## change .1 to variable...
                    if l1 > self.max_len_chain_difference or l1/float(len(s1)) > .1:
                        continue
                    if r1 > self.max_len_chain_difference or r1/float(len(s1)) > .1:
                        continue
                    if l2 > self.max_len_chain_difference or l2/float(len(s2)) > .1:
                        continue
                    if r2 > self.max_len_chain_difference or r2/float(len(s2)) > .1:
                        continue

                    s1 = s1[l:len(s1)-r]
                    s2 = s2[l:len(s2)-r]

                    ## continue if insertions/deletions
                    if '-' in s1 or '-' in s2:
                        continue

                    n_chainmutations, l_chainmutations = self.point_mutations(s1, s2)

                if n_chainmutations <= self.max_mutations:
                    if chain1 not in d_chains_interpdb_sequence_similar.keys():
                        d_chains_interpdb_sequence_similar[chain1] = {}
                    if chain2 in d_chains_interpdb_sequence_similar[chain1].keys():
                        notexpected
                    d_chains_interpdb_sequence_similar[chain1][chain2] = {
                        'l1':l1,'l2':l2,
                        's1':s1,'s2':s2,
                        'l_mutations':l_chainmutations,
                        'r1':0,'r2':0,
                        }

##        print d_chains_interpdb_sequence_similar
    
        return d_chains_interpdb_sequence_similar


    def res_name2res_symbol(self, res_name):

        ## this function should take into account nucleotides
        ## unfortunately there is a conflict between:
        ## G: glycine and guanosine
        ## A: alanine and adenosine
        ## T: threonine and thymidine
        ## C: cysteine and cytidine
        ## N: aspargine and unknown nucleotide residue

        if res_name in self.d_res.keys():
            symbol = self.d_res[res_name]
        else:
            symbol = 'X'

        return symbol
    

    def ATOM2seq(self, d_pdb, chain, d_seq,):

## incorrect d_res_nos will be returned from ATOM2seq for 2bfk (vs 2bfl), chain A becauseof ASN61A in REMARK465 records!!!

        d_res_nos = {}
        seq = ''
        ATOMrespos = 0
        res_nos = d_pdb['chains'][chain]['residues'].keys()
        res_nos.sort()
        for i in range(len(res_nos)):
            res_no = res_nos[i]

            for i in range(len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])):
                iCode = d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][i]
                res_name = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                
                ##
                ## do not append to sequence if hetID is not a MODRES
                ## e.g. TRP in 1utv.pdb
                ##
                if chain in d_seq['HET'].keys():
                    if res_no in d_seq['HET'][chain].keys():
                        if iCode in d_seq['HET'][chain][res_no].keys():
                            if res_name != d_seq['HET'][chain][res_no][iCode]:
                                None ## duplicate IDs for ATOM and HETATM
                            else:
                                if not 'MODRES' in d_seq.keys():
                                    continue
                                if not chain in d_seq['MODRES'].keys():
                                    continue
                                elif res_no not in d_seq['MODRES'][chain].keys():
                                    continue
                                elif iCode not in d_seq['MODRES'][chain][res_no].keys():
                                    continue
                                else:
                                    None
                                
                d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':iCode}
                seq += self.res_name2res_symbol(res_name)
                ATOMrespos += 1

        return seq, d_res_nos


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


    def append_missing_residues_to_sequence(self, ATOMseq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq):

        for SEQRESpos in range(SEQRESrange1,SEQRESrange2):
            ATOMseq += SEQRESseq[SEQRESpos]
            d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}

        return ATOMseq, d_res_nos_SEQRES

    
    def identify_missing_nonterminal_residues(self, d_pdb, pdb, chain, SEQRESseq, ATOMseqgaplen, d_res_nos_ATOM):

        d_res_nos_SEQRES = {}
        seq = ''

        ## append N-terminal gap between ATOMseq and SEQRESseq
        SEQRESrange1 = 0
        SEQRESrange2 = ATOMseqgaplen
        seq, d_res_nos_SEQRES = self.append_missing_residues_to_sequence(seq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq)

        ## reset counters and indexes
        ATOMseqrespos = ATOMseqgaplen
        prevATOMseqrespos = 0
        index2 = 0
        ## initiate loop
        SEQRESrange1 = ATOMseqgaplen
        SEQRESrange2 = len(d_res_nos_ATOM.keys())+ATOMseqgaplen
        for SEQRESpos in range(SEQRESrange1,SEQRESrange2):
            ATOMpos = SEQRESpos-ATOMseqgaplen
            res_no = d_res_nos_ATOM[ATOMpos]['res_no']
            iCode = d_res_nos_ATOM[ATOMpos]['iCode']
            res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
            if res_name == 'HOH':
                continue
            res_symbol = self.res_name2res_symbol(res_name)

            ## append gap between ATOMseq and SEQRESseq
            if seq+res_symbol != SEQRESseq[:SEQRESpos+1]:
                ## e.g. 2i0b.pdb
                print d_pdb[pdb]['chains'][chain]['residues'][res_no]['l_iCodes']
                print pdb, chain, res_no, iCode, res_name
                print seq+res_symbol
                print SEQRESseq[:SEQRESpos+1]
##                if d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] == True:
##                    continue
                if d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] == 'HETATM':
                    continue
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':'-','iCode':'-'}
                ATOMseqgaplen += 1 ## not a gap but a reversal in the case of 2bfk.pdb!!!
                try:
                    seq += SEQRESseq[SEQRESpos]
                except:
                    print
                    print pdb, chain, res_no, iCode, res_name, SEQRESpos, SEQRESrange1, SEQRESrange2
                    print pdb, chain, d_res_nos_ATOM[ATOMpos-1]['res_no'], d_res_nos_ATOM[ATOMpos-1]['iCode']
                    print len(SEQRESseq), SEQRESseq
                    print len(seq), seq
                    print self.pdb1, self.pdb2
                    print SEQRESseq[SEQRESrange1:]
                    print seq[SEQRESrange1:]
                    for i in range(len(seq)):
                        if seq[:i] != SEQRESseq[:i]:
                            print seq[:i]
                            print SEQRESseq[:i]
                            stop1
                    stop2
            ## append if no gap
            else:
                d_res_nos_SEQRES[SEQRESpos] = {'res_no':res_no,'iCode':iCode}
                seq += res_symbol

        ## append C-terminal gap between ATOMseq and SEQRESseq
        if len(seq) != len(SEQRESseq):
            SEQRESrange1 = SEQRESpos+1
            SEQRESrange2 = len(SEQRESseq)
            seq, d_res_nos_SEQRES = self.append_missing_residues_to_sequence(seq, d_res_nos_SEQRES, SEQRESrange1, SEQRESrange2, SEQRESseq)

##        if seq != SEQRESseq:
##            print pdb
##            print seq
##            print SEQRESseq
##            notexpected

        return seq, d_res_nos_SEQRES

    def identify_missing_terminal_residues(self, d_pdb, pdb, chain, ATOMseq, SEQRESseq):

##        print pdb, chain
##        print ATOMseq
##        print SEQRESseq

        for i in range(1,len(ATOMseq)+1):

            ## index the first occurence of the ATOMseq in the SEQRESseq
            try:
                index = SEQRESseq.index(ATOMseq[:i])
                ## index the next occurence of the ATOMseq in the SEQRESseq
                try:
                    SEQRESseq[index+1:].index(ATOMseq[:i])
                ## break if only one occurence of the ATOMseq in the SEQRESseq
                except:
                    break
            ## break if ATOMseq not occuring in the SEQRESseq
            except:
                print pdb, chain
                print ATOMseq[:i]
                print ATOMseq
                print SEQRESseq
                notexpected
                break

        return index


    def alignATOMseq(self, d_pdb, d_seq, pdb1, pdb2, chain1, chain2):

        import sys
        sys.path.append('/home/people/tc/python/EAT_DB/')
        sys.path.append('../../EAT_DB/')
        import sequence_alignment

        ##
        ## identify missing residues not mentioned in the REMARK465 records
        ## by alignment of SEQRESseq and ATOMseq
        ## and add gaps to the ATOMseq
        ##
        ## use the nontransformed chain IDs for sequence alignment
        d_ATOM_seqs = {
            'pdb1':{'pdb':pdb1,'chain':chain1[0]},
            'pdb2':{'pdb':pdb2,'chain':chain2[0]},
            }
        for key in d_ATOM_seqs.keys():
            pdb = d_ATOM_seqs[key]['pdb']
            chain = d_ATOM_seqs[key]['chain']
            SEQRESseq = d_seq[pdb]['chains'][chain]['seq']
            ATOMseq,d_res_nos = self.ATOM2seq(d_pdb[pdb], chain, d_seq[pdb])
            ## find missing residues, Nterminal; align Nterminal SEQRESseq and ATOMseq
            ATOMseqindentation = self.identify_missing_terminal_residues(d_pdb, pdb, chain, ATOMseq, SEQRESseq)
            ## find missing residues, Cterminal or nonterminal; align SEQRESseq and ATOMseq
            ATOMseq,d_res_nos = self.identify_missing_nonterminal_residues(d_pdb, pdb, chain, SEQRESseq, ATOMseqindentation, d_res_nos)
            d_ATOM_seqs[key]['ATOMseq'] = ATOMseq
            d_ATOM_seqs[key]['d_res_nos'] = d_res_nos
        ATOMseq1 = d_ATOM_seqs['pdb1']['ATOMseq']
        ATOMseq2 = d_ATOM_seqs['pdb2']['ATOMseq']
        d_res_nos1 = d_ATOM_seqs['pdb1']['d_res_nos']
        d_res_nos2 = d_ATOM_seqs['pdb2']['d_res_nos']

        ##
        ## remove terminal residues from the ATOMseq
        ##
        if len(ATOMseq1) == len(ATOMseq2):

            l1 = 0
            l2 = 0

        else:

            ## 2nd slow sequence comparison (ATOMseq)

            print pdb1, pdb2, chain1, chain2, 'begin seq aln'
            instance = sequence_alignment.NW(ATOMseq1,ATOMseq2)
            s1,s2 = ATOMs1,ATOMs2 = instance.Align(verbose=False)[:2]
            print 'end seq aln'

            l1 = len(s1)-len(s1.lstrip('-'))
            l2 = len(s2)-len(s2.lstrip('-'))
            r1 = len(s1)-len(s1.rstrip('-'))
            r2 = len(s2)-len(s2.rstrip('-'))
            if r2 == 0:
                ATOMseq1 = ATOMseq1[l2:]
            else:
                ATOMseq1 = ATOMseq1[l2:-r2]
            if r1 == 0:
                ATOMseq2 = ATOMseq2[l1:]
            else:
                ATOMseq2 = ATOMseq2[l1:-r1]

        return d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2


    def ATOMrecords2coordinates(
        self, d_pdb, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, len_ATOMseq,):

        import Numeric, math, copy

        rescount = 0

        coordinates1 = []
        coordinates2 = []
        lines1 = []
        lines2 = []

        for SEQRESpos1 in range(l2,len_ATOMseq+l2):

            SEQRESpos2 = SEQRESpos1+l1-l2
##            print pdb1, pdb2, chain1, chain2, SEQRESpos1, SEQRESpos2
            res_no1 = d_res_nos1[SEQRESpos1]['res_no']
            res_no2 = d_res_nos2[SEQRESpos2]['res_no']
            if res_no1 == '-' or res_no2 == '-':
                continue
            iCode1 = d_res_nos1[SEQRESpos1]['iCode']
            iCode2 = d_res_nos2[SEQRESpos2]['iCode']

            if 'REMARK' in d_pdb[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1].keys():
                continue
            if 'REMARK' in d_pdb[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2].keys():
                continue

            d_resname = {
                'pdb1':{'pdb':pdb1,'chain':chain1,'res_no':res_no1,'iCode':iCode1},
                'pdb2':{'pdb':pdb2,'chain':chain2,'res_no':res_no2,'iCode':iCode2},
                }
            for key in d_resname.keys():
                pdb = d_resname[key]['pdb']
                chain = d_resname[key]['chain']
                res_no = d_resname[key]['res_no']
                iCode = d_resname[key]['iCode']
                res_name = d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                d_atoms = copy.deepcopy(d_pdb[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'])
                d_resname[key]['res_name'] = res_name
                d_resname[key]['d_atoms'] = d_atoms
            res_name1 = d_resname['pdb1']['res_name']
            res_name2 = d_resname['pdb2']['res_name']
            d_atoms1 = d_resname['pdb1']['d_atoms']
            d_atoms2 = d_resname['pdb2']['d_atoms']

            if res_name1 in self.d_modres.keys():
                res_name1 = self.d_modres[res_name1]
            if res_name2 in self.d_modres.keys():
                res_name2 = self.d_modres[res_name2]
##            if res_name1 != res_name2:
##                print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2
##                stop
##                fd = open('different_resnames.txt','a')
##                fd.write('%s %s %s %s %s %s %s %s\n' %(pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2))
##                fd.close()

            mutation = False
            if res_name1 != res_name2:
                mutation = True

            rescount += 1
            for atom_name in d_atoms1.keys():
                if 'CA' not in d_atoms1.keys():
                    stop
                if mutation == True and atom_name not in ['N','CA','C','O']:
                    continue
                if atom_name not in d_atoms2.keys():
                    continue
                if 'REMARK' in d_atoms1[atom_name].keys():
                    continue
                if 'REMARK' in d_atoms2[atom_name].keys():
                    continue
                if atom_name not in ['CA']:
                    continue
                ## append coordinates to list of coordinates
                coordinate1 = d_atoms1[atom_name]['coordinate']
                coordinate2 = d_atoms2[atom_name]['coordinate']
                coordinates1 += [coordinate1]
                coordinates2 += [coordinate2]

        return coordinates1, coordinates2


    def parse_coordinates(self, s_pdb, d_pdb, d_seq, s_chain):

##        print 'parsing coordinates'

        ## read lines
##        fd = open('%s%s.pdb' %(self.path_pdb, s_pdb.upper()),'r')
        fd = open('%s%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
        lines = fd.readlines()
        fd.close()
        if s_pdb not in d_pdb.keys():
            d_coordinates = self.parse_pdbcoordinatesection(lines, s_pdb, d_seq, s_chain)
            d_pdb[s_pdb] = d_coordinates

        return d_pdb



    def parse_pdbheader(self, lines, s_pdb):

        s_pdb = s_pdb.lower()
        ## import sequence alignment
        import sys
        sys.path.append('/home/people/tc/svn/EAT_DB/')
        import sequence_alignment

        ## parser written on the assumption that SEQRES is mandatory if ATOM records exist

        d_seq = {'HET':{}}
        d_conect = {}
        l_hetatms = []
        parse_atoms = False

        ## insertion chain, res_no, res_name
        prev_chain = ''
        d_insertions = {}
        biounit = 'N/A'

        for i in range(len(lines)):
            line = lines[i]

            record = line[:6].strip()

            if record == 'ATOM': ## section 9
                continue

            elif record == 'SOURCE':
                print s_pdb, line,

            elif record == 'DBREF':
                print s_pdb, line,

            elif record == 'SEQADV':
                print s_pdb, line,

            elif record == 'HETATM': ## section 9
                continue

            elif record == 'REMARK': ## section 2
                d_seq = self.parse_recordREMARK(d_seq, line, i, lines)

            elif record == 'SEQRES': ## section 3
                d_seq = self.parse_recordSEQRES(line, d_seq)

            elif record == 'HET': ## section 4
                hetID = line[7:10].strip()
                ## continue if water
                if hetID in ['HOH','H2O','DOD','D2O']: ## D2O in 2JAJ
                    continue
                chain = line[12]
                res_no = int(line[13:17])
                iCode = line[17]
                if 'HET' not in d_seq.keys():
                    d_seq['HET'] = {}
                if chain not in d_seq['HET'].keys():
                    d_seq['HET'][chain] = {}
                if res_no not in d_seq['HET'][chain].keys():
                    d_seq['HET'][chain][res_no] = {}
                if iCode not in d_seq['HET'][chain][res_no].keys():
                    d_seq['HET'][chain][res_no][iCode] = hetID

            elif record == 'MODRES':
                hetID = line[12:15].strip()
                chain = line[16]
                res_no = int(line[18:22])
                iCode = line[22]
                res_name = line[24:27].strip()
                txt = line[29:80].strip()
                if hetID in set(self.d_res.keys())-set(['MSE']) and res_name in set(self.d_res.keys())-set(['MSE']):
                    continue
                if 'MODRES' not in d_seq.keys():
                    d_seq['MODRES'] = {}
                if chain not in d_seq['MODRES'].keys():
                    d_seq['MODRES'][chain] = {}
                if res_no not in d_seq['MODRES'][chain].keys():
                    d_seq['MODRES'][chain][res_no] = {}
                if iCode not in d_seq['MODRES'][chain][res_no].keys():
                    d_seq['MODRES'][chain][res_no][iCode] = hetID
                elif hetID != d_seq['MODRES'][chain][res_no][iCode]:
                    print line, s_pdb
                    notexpected

            elif record == 'TITLE': ## section 2
                if not 'TITLE' in d_seq.keys():
                    d_seq['TITLE'] = line[10:].strip()
                else:
                    if d_seq['TITLE'][-1] == '-':
                        d_seq['TITLE'] += line[10:].strip()
                    else:
                        d_seq['TITLE'] += ' '+line[10:].strip()

            elif record == 'HEADER':
                d_seq['HEADER'] = line[10:50].strip()

        return d_seq


    def parse_recordSEQRES(self, line, d_seq):

        chain = line[11]

        if 'chains' not in d_seq:
            d_seq['chains'] = {}
        if chain not in d_seq['chains'].keys():
            d_seq['chains'][chain] = {}
        if not 'type' in d_seq['chains'][chain].keys():
            d_seq['chains'][chain]['type'] = 'N/A'

        residues = line[19:70].split()

        for i in range(len(residues)):
            residue = residues[i]
            if residue in self.d_res.keys():
                if d_seq['chains'][chain]['type'] == 'N/A':
                    d_seq['chains'][chain]['type'] = 'peptide'
                residues[i] = self.d_res[residue]
            elif residue in ['C','A','U','G','I','DC','DA','DT','DG','DI','N']: ## N is any 5'-monophosphate nucleotide
                if d_seq['chains'][chain]['type'] == 'N/A':
                    d_seq['chains'][chain]['type'] = 'nucleotide'
                residues[i] = residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_seq['chains'][chain]['type'] == 'N/A':
                    d_seq['chains'][chain]['type'] = 'saccharide'
                residues[i] = residue
            else:
                if residue == 'UNK': ## e.g. 1pny.pdb
                    if d_seq['chains'][chain]['type'] == 'N/A':
                        d_seq['chains'][chain]['type'] = 'peptide'
                residues[i] = 'X'

        if 'seq' not in d_seq['chains'][chain].keys():
            d_seq['chains'][chain]['seq'] = ''
        d_seq['chains'][chain]['seq'] += ''.join(residues)

        return d_seq


    def parse_pdbcoordinatesection(self, lines, s_pdb, d_seq, s_chain):

        print 'parsing coordinates of %s' %(s_pdb)

        s_pdb = s_pdb.lower()

        ##
        ## set dictionaries
        ##
        d_atomnos = {}
        d_CONECT = {}

        d_coordinates = {}

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'ATOM':
                chain = line[21]
                if chain != s_chain:
                    continue
                d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, 'ATOM', s_pdb)
                d_atomnos[d_line['atom_no']] = d_line

            elif record == 'HETATM':
                MODRES = False
                res_name = line[17:20].strip()
                chain = line[21]
                res_no = int(line[22:26])
                iCode = line[26]
                if chain != s_chain:
                    continue
                if 'MODRES' in d_seq.keys():
                    if chain in d_seq['MODRES'].keys():
                        if res_no in d_seq['MODRES'][chain].keys():
                            if iCode in d_seq['MODRES'][chain][res_no].keys():
                                if res_name != d_seq['MODRES'][chain][res_no][iCode]:
                                    print res_name, d_seq['MODRES'][chain][res_no][iCode]
                                    notexpected
                                MODRES = True
                if MODRES == False:
                    continue
                ## water
                if res_name in ['HOH','H2O','DOD','D2O']:
                    d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, 'HETATM', s_pdb)
                ## modified residue of polypeptide or polynucleotide
                elif res_name in self.d_modres.keys(): ## e.g. 1gcj.pdb
                    d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, 'ATOM', s_pdb)
                ## (poly)saccharide or other hetero compound
                else:
                    d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_seq, 'HETATM', s_pdb)
                atom_no = d_line['atom_no']
                d_atomnos[atom_no] = d_line

            elif record == 'REMARK':
                remark = int(line[6:10])
                if remark == 465:
                    d_coordinates = self.parse_recordREMARK465(line, d_coordinates, lines, i)
                elif remark == 470:
                    d_coordinates = self.parse_recordREMARK470(line, d_coordinates, lines, i)

            elif record == 'HET':
                hetID = line[7:10].strip()
                if 'HET' not in d_coordinates.keys():
                    d_coordinates['HET'] = set()
                d_coordinates['HET'] |= set([hetID])
##                chain = line[12]
##                res_no = int(line[13:17])
##                iCode = line[17]
##                if iCode != ' ':
##                    print s_pdb
##                    stop
##                n_atoms = int(line[20:25])
##                description = line[30:70]
##                if chain not in d_hetero['chains'].keys():
##                    d_hetero['chains'][chain] = {}
##                if 'residues' not in d_hetero['chains'][chain].keys():
##                    d_hetero['chains'][chain]['residues'] = {}
##                if res_no not in d_hetero['chains'][chain]['residues'].keys():
##                    d_hetero['chains'][chain]['residues'][res_no] = {}
##                if 'iCodes' not in d_hetero['chains'][chain]['residues'][res_no].keys():
##                    d_hetero['chains'][chain]['residues'][res_no]['iCodes'] = {}
##                d_hetero['chains'][chain]['residues'][res_no]['iCodes'][iCode] = hetID
                
            elif record == 'MODEL':
                model = int(line.split()[1])


        return d_coordinates


    def monomertranslation(self,monomer,d_coordinates):

        chain = monomer[0]
        res_no = int(monomer[1:-1])
        iCode = monomer[-1]
        res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
        if res_name in self.d_saccharides.keys():
            res_name = self.d_saccharides[res_name]['stereo']

        return res_name


    def parse_recordREMARK(self, d_seq, line, i, lines):

        remark = int(line[6:10])

        if remark == 465: ## missing residues

            d_seq['REMARK465'] = True

        elif remark == 470: ## missing atoms

            d_seq['REMARK470'] = True

        return d_seq


    def parse_recordREMARK465(self, line, d_pdb, lines, i):

        ## missing residues

        if line[10:].strip() in ['M RES C SSSEQI','M RES C  SSEQI']:

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 465':
                    break

                try:
                    model = int(lines[j][12:14])
                except:
                    model = 'N/A'
                res_name = lines[j][15:18].strip()
                chain = lines[j][19]
                res_no = int(lines[j][22:26])
                iCode = lines[j][26]


                if not 'chains' in d_pdb.keys():
                    d_pdb['chains'] = {}
                if not chain in d_pdb['chains'].keys():
                    d_pdb['chains'][chain] = {}
                if not 'residues' in d_pdb['chains'][chain].keys():
                    d_pdb['chains'][chain]['residues'] = {}
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
                if not iCode in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## iCode > REMARK
                if not 'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] = True

        return d_pdb


    def parse_recordREMARK470(self, line, d_pdb, lines, i):

        ## missing atoms

        ## the latter equation is only to acommodate for 1fvk.pdb
        if line[10:].strip() == 'M RES CSSEQI  ATOMS' or line[10:].strip() == 'M RES C SEQI  ATOMS':

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 470':
                    break

                ## model M
                try:
                    model = int(lines[j][11:13])
                except:
                    model = 'N/A'

                ## res_name RES
                res_name = lines[j][15:18].strip()

                ## chain C
                chain = lines[j][19]

                ## res_no SSEQ
                try:
                    res_no = int(lines[j][20:24])
                except:
                    res_no = lines[j][20:24].strip()

                ## iCode I
                iCode = lines[j][24]

                ## atoms ATOMS
                atoms = lines[j][25:].split()

                ##
                ## write to dictionary
                ##
                if not 'chains' in d_pdb.keys():
                    d_pdb['chains'] = {}
                if not chain in d_pdb['chains'].keys():
                    d_pdb['chains'][chain] = {}
                if not 'residues' in d_pdb['chains'][chain].keys():
                    d_pdb['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['chains'][chain]['residues'].keys():
                    d_pdb['chains'][chain]['residues'][res_no] = {}

                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']) == 0:
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
                elif iCode not in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
                    ## e.g. 2fs4 (chain A, res_no 162, iCode " ")
                    if iCode == ' ':
                        d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = [iCode]+d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']
                    else:
                        d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

                ## iCode > atoms
                if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                ## atoms > atom_name > coordinate
                for atom_name in atoms:
                    if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                    ## iCode > REMARK
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['REMARK'] = True

                ## iCode > res_name
                if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name


        elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
            print self.pdb1
            print self.pdb2
            notexpected

        return d_pdb


    def parse_recordATOM(self, line, d_pdb, lines, i, d_seq, record, s_pdb):

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
        coordinate = Numeric.array([x, y, z])

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
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name:
            print res_name, d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
            print chain, res_no, iCode, atom_altloc
            print line
            stop
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
        if len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']) == 0:
            d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
        elif iCode not in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
            d_pdb = self.identify_iCode_sequence(d_pdb, chain, res_no, iCode, res_name, d_seq, s_pdb)

        ## iCode > atoms
        if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
        ## atoms > atom_name > coordinate
        if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}

        ## iCode > record
        if not 'record' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = record

        
        return d_pdb, {'chain':chain,'atom_no':atom_no,'atom_name':atom_name,'res_name':res_name,'res_no':res_no,'iCode':iCode,'altloc':atom_altloc}


    def identify_iCode_sequence(self, d_pdb, chain, res_no, iCode, res_name, d_seq, s_pdb):
        
        l_iCodes = list(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
        d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
        print '"%s"' %(iCode), l_iCodes
        for iCode_prev in l_iCodes:
            print chain, res_no, iCode, iCode_prev, l_iCodes
            index_alphabet = self.s_alphabet.index(iCode)
            if (
                ## REMARK465
                'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev].keys() or
                ## REMARK470
                {'REMARK':True} in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev]['atoms'].values()
                ):
                ATOMseq,d_res_nos = self.ATOM2seq(d_pdb, chain, d_seq,)
                res_symbol = self.res_name2res_symbol(res_name)
                ## REMARK residues before ATOM residues (N-terminal)
                ## e.g. 1jqz.pdb
                if (
                    len(ATOMseq) == 0
                    ):
                    None
                ## REMARK465 residues after ATOM residues
                ## e.g. 1nuo.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_seq['chains'][chain]['seq'][len(ATOMseq)] and
                    ATOMseq == d_seq['chains'][chain]['seq'][:len(ATOMseq)]
                    ):
                    l_iCodes = [iCode]+d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][:-1]
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## REMARK470 residues after ATOM residues
                ## e.g. 2lve.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_seq['chains'][chain]['seq'][len(ATOMseq)+index_alphabet] and
                    ATOMseq == d_seq['chains'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' in l_iCodes
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## REMARK470 residues after ATOM residues
                ## e.g. 2j5q (chain B, res_no 54, iCode C)
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_seq['chains'][chain]['seq'][len(ATOMseq)+index_alphabet-1] and
                    ATOMseq == d_seq['chains'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' not in l_iCodes 
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## REMARK465 residues before ATOM residues
                ## e.g. 1uij.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_seq['chains'][chain]['seq'][len(ATOMseq)+len(l_iCodes)] and
                    ATOMseq == d_seq['chains'][chain]['seq'][:len(ATOMseq)]
                    ):
                    None
                else:
                    print iCode, res_symbol, index_alphabet, l_iCodes
                    print d_seq['chains'][chain]['seq'][len(ATOMseq)+index_alphabet]
                    print
                    print len(ATOMseq) > 0
                    print res_symbol == d_seq['chains'][chain]['seq'][len(ATOMseq)+index_alphabet-1]
                    print ATOMseq == d_seq['chains'][chain]['seq'][:len(ATOMseq)]
                    print self.s_alphabet[index_alphabet-1] in l_iCodes
                    print 
                    print 1, ATOMseq
                    print 2, d_seq['chains'][chain]['seq']
                    print chain, res_no, iCode, iCode_prev, res_name
                    expected
                break

        return d_pdb


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
        
        ## HETATM res_names for which coordinates are parsed
        self.d_modres = {
            'MSE':'MET', ## selenomethionine
            }

        self.d_res3 = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

        self.nontransformationmatrix = [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]

        self.s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25

        self.path_pdb = '/oxygenase_local/data/pdb/'
        
        return

if __name__ == '__main__':
    instance_molmovdb = molmovdb()
    instance_molmovdb.main()

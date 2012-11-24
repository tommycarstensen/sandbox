#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

def main():

    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import gnuplot

    min_count = 2
    phipsi_range = 1
    phipsi_step = 5

    d_residues = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
        'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
        'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
        'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
        }

    d_ramachandran = parse_dihedrals(d_residues, plot = True)

    print 'preparing dictionary'
    for res1 in d_residues.values():
        for res2 in d_residues.values():
            if res1 == res2:
                continue
            d_ramachandran[res1+res2] = {}
            for phi in range(-180,180,phipsi_step,):
                d_ramachandran[res1+res2][phi] = {}
                for psi in range(-180,180,phipsi_step,):
                    d_ramachandran[res1+res2][phi][psi] = 0
    print 'prepared dictionary'

    fd = open('phipsi.txt','r')
    lines = fd.readlines()
    fd.close()

    print 'looping over', len(lines), 'lines'

    for line in lines:

        if len(line.split()) == 12:
            if (
                'N/A' in [
                    line.split()[8],
                    line.split()[9],
                    line.split()[10],
                    line.split()[11],
                    ]
                ):
                continue
            res1 = line.split()[6]
            res2 = line.split()[7]
            phi1 = float(line.split()[8])
            psi1 = float(line.split()[9])
            phi2 = float(line.split()[10])
            psi2 = float(line.split()[11])
        else:
            print line
            stop

        if res1 == 'X' or res2 == 'X':
            print line
            continue

##        count = sum_ramachandran(d_ramachandran,res1,phi1,psi1,phipsi_range,)
##        if count < count_min:
##            print line
##            print count
##            if res1 != 'G':
##                stop1
##
##        count = sum_ramachandran(d_ramachandran,res2,phi2,psi2,phipsi_range,)
##        if count < count_min:
##            print line
##            print count
##            if res2 != 'G':
##                stop2

        phi1 = round_angle(phi1,phipsi_step)
        psi1 = round_angle(psi1,phipsi_step)
        phi2 = round_angle(phi2,phipsi_step)
        psi2 = round_angle(psi2,phipsi_step)
        d_ramachandran[res2+res1][phi1][psi1] += 1
        d_ramachandran[res1+res2][phi2][psi2] += 1


    ##
    ## ramachandran plots
    ##
    for key in d_ramachandran.keys():
        if len(key) == 1:
            continue
        import os
        if os.path.isfile('plot_phipsi_%s.png' %(key)):
            os.remove('plot_phipsi_%s.png' %(key))
        print key
        l_gnuplot = []
        max_count = 0
        sum_count = 0
        for phi in range(-180,180,phipsi_step,):
            for psi in range(-180,180,phipsi_step,):
                 l_gnuplot += ['%s %s %s\n' %(phi,psi,d_ramachandran[key][phi][psi])]
                 if d_ramachandran[key][phi][psi] > max_count:
                     max_count = d_ramachandran[key][phi][psi]
                     sum_count += d_ramachandran[key][phi][psi]
            l_gnuplot += ['\n']
        if max_count > 1 and sum_count > 3:
            gnuplot.contour_plot(
                'plot_phipsi_%s' %(key), l_gnuplot,
                title='%s' %(key), xtitle='{/Symbol f}', ytitle='{/Symbol y}',
                x1=-180, x2=180, y1=-180, y2=180,
                z1=0,
                )


    return


def round_angle(angle,phipsi_step,):

    if angle == 180.:
        angle = -180.
    else:
        angle = phipsi_step*int(angle/phipsi_step)

    return angle


def sum_ramachandran(d_ramachandran,res,phi,psi,phipsi_range,):

    if phi == 180.:
        phi = -180.
    else:
        phi = int(phi)
    if psi == 180.:
        psi = 180.
    else:
        psi = int(psi)

    count = 0
    
    for i in range(-phipsi_range,phipsi_range+1):
        for j in range(-phipsi_range,phipsi_range+1):
            if not (
                int(phi)+i < -180. or
                int(phi)+i >= 180. or
                int(psi)+j < -180. or
                int(psi)+j >= 180.
                ):
                count += d_ramachandran[res][int(phi)+i][int(psi)+j]

    return count


def parse_dihedrals(d_residues, plot = False):

    import math,sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import gnuplot

    d_ramachandran = {}
    for res in d_residues.values():
        d_ramachandran[res] = {}
        for phi in range(-180,180):
            d_ramachandran[res][phi] = {}
            for psi in range(-180,180):
                d_ramachandran[res][phi][psi] = 0

    for res in d_residues.values():
        print 'read', res
        ##
        ## read lines
        ##
        fd = open('phipsi_all_%s.txt' %(res),'r')
        lines = fd.readlines()
        fd.close()
        ##
        ## parse lines
        ##
        for line in lines:
            if 'N/A' in [line.split()[0],line.split()[1]]:
                continue
            phi = float(line.split()[0])
            psi = float(line.split()[1])
            if phi == 180.:
                phi = -180.
            if psi == 180.:
                psi = -180.
            d_ramachandran[res][int(phi)][int(psi)] += 1

    ##
    ## ramachandran plots
    ##
    if plot == True:
        for res in d_ramachandran.keys():
            print 'plot', res
            l_gnuplot_res = []
            for phi in range(-180,180):
                for psi in range(-180,180):
                    l_gnuplot_res += ['%s %s %s\n' %(phi,psi,d_ramachandran[res][phi][psi])]
                l_gnuplot_res += ['\n']
            gnuplot.contour_plot('plot_phipsi_%s' %(res), l_gnuplot_res, title='phipsi_%s' %(res), xtitle='{/Symbol f}', ytitle='{/Symbol y}')
        l_gnuplot = []
        for phi in range(-180,180):
            for psi in range(-180,180):
                for res in d_ramachandran.keys():
                    l_gnuplot += ['%s %s %s\n' %(phi,psi,d_ramachandran[res][phi][psi])]
            l_gnuplot += ['\n']
        gnuplot.contour_plot('plot_phipsi', l_gnuplot, title='phipsi', xtitle='{/Symbol f}', ytitle='{/Symbol y}')



##    d_dihedral = {
##        'A':{'phi':[],'psi':[]},
##        'C':{'phi':[],'psi':[]},
##        'D':{'phi':[],'psi':[]},
##        'E':{'phi':[],'psi':[]},
##        'F':{'phi':[],'psi':[]},
##        'G':{'phi':[],'psi':[]},
##        'H':{'phi':[],'psi':[]},
##        'I':{'phi':[],'psi':[]},
##        'K':{'phi':[],'psi':[]},
##        'L':{'phi':[],'psi':[]},
##        'M':{'phi':[],'psi':[]},
##        'N':{'phi':[],'psi':[]},
##        'P':{'phi':[],'psi':[]},
##        'Q':{'phi':[],'psi':[]},
##        'R':{'phi':[],'psi':[]},
##        'S':{'phi':[],'psi':[]},
##        'T':{'phi':[],'psi':[]},
##        'V':{'phi':[],'psi':[]},
##        'W':{'phi':[],'psi':[]},
##        'Y':{'phi':[],'psi':[]},
##        }
##
##    for line in lines:
##        residue = d_residues[line.split()[0]]
##        if 'N/A' in [line.split()[1],line.split()[2]]:
##            continue
##        phi = float(line.split()[1])
##        psi = float(line.split()[2])
##        d_dihedral[residue]['phi'] += [phi]
##        d_dihedral[residue]['psi'] += [psi]
##
##    for residue in d_dihedral.keys():
##        for dihedral in ['phi','psi',]:
##            sumx = 0
##            sumxx = 0
##            n = len(d_dihedral[residue][dihedral])
##            for x in d_dihedral[residue][dihedral]:
##                sumx += x
##                sumxx += x**2
##            mean = sumx/n
##            stddev = math.sqrt(n*sumxx-(sumx**2))/n
##            d_dihedral[residue][dihedral] = {'mean':mean,'stddev':stddev}
##    print d_dihedral

    return d_ramachandran


if __name__ == '__main__':
    main()

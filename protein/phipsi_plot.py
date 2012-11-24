def main():

    l_pdbs = [
##        '2lzm','150l',
##        '1l96','1l97',
##        '2lzt','3lzt',
        
        ]

##    l_pdbs = range(50)

##    l_pdbs = []
##    for i in range(101):
##        l_pdbs += ['2LZM_wt_trjconv_%i_rotated' %(i)]

    res_range = range(1,1+161+3)
    l_connectors = [0,49,]

    l_pdbs = ['9aat','1ama',]
    res_range = range(313,344+1)
    l_connectors = l_pdbs

    chain = 'A'
    d_data = {}

    for pdb in l_pdbs:

        d_data[pdb], lines = calculate_phipsi(pdb, chain = chain, res_range = res_range,)

        fd = open('%s.gnuplotdata' %(pdb),'w')
        fd.writelines(lines)
        fd.close()

##    ## write gnuplot files per residue
##    for res_no in range(3,163):
##        lines = []
##        for pdb in l_pdbs:
##            phi = d_data[pdb][res_no][0]
##            psi = d_data[pdb][res_no][1]
##            s = '%s %s %s\n' %(phi,psi,res_no)
##            lines += [s]
##            fd = open('resno%s.gnuplotdata' %(res_no),'w')
##            fd.writelines(lines)
##            fd.close()
##    l_pdbs = []
##    for i in range(3,163):
##        l_pdbs += ['resno%i' %(i)]

    lines = []
    for i in range(len(l_pdbs)-1):
        pdb1 = l_pdbs[i]
        if pdb1 not in l_connectors:
            continue
        for j in range(i+1,len(l_pdbs)):
            pdb2 = l_pdbs[j]
            if pdb2 not in l_connectors:
                continue

            for res_no in d_data[pdb1].keys():
##                if res_no not in d_data[pdb2].keys():
##                    continue
##                if res_no not in res_range:
##                    continue
                phi1 = d_data[pdb1][res_no][0]
                psi1 = d_data[pdb1][res_no][1]
                phi2 = d_data[pdb2][res_no][0]
                psi2 = d_data[pdb2][res_no][1]
                if abs(phi1-phi2) > 180:
                    if phi1 < phi2:
                        phi3 = phi2-360
                        phi4 = phi1+360
                    else:
                        phi3 = phi2+360
                        phi4 = phi1-360
                else:
                    phi3 = phi2
                    phi4 = phi1
                if abs(psi1-psi2) > 180:
                    if psi1 < psi2:
                        psi3 = psi2-360
                        psi4 = psi1+360
                    else:
                        psi3 = psi2+360
                        psi4 = psi1-360
                else:
                    psi3 = psi2
                    psi4 = psi1
                lines += ['%s %s\n%s %s\n\n' %(phi1,psi1,phi3,psi3)]
                lines += ['%s %s\n%s %s\n\n' %(phi2,psi2,phi4,psi4)]
                    
    fd = open('connectors.gnuplotdata','w')
    fd.writelines(lines)
    fd.close()

    scatter_plot_2d_multiple_colors(
        l_pdbs,
        xlabel = '{/Symbol f}',
        ylabel = '{/Symbol y}',
        xmin = -180., xmax = 180.,
        ymin = -180., ymax = 180.,
        )

    return


def calculate_phipsi(pdb, chain = None, res_range = []):

    d_coordinates = parse_backbone_coordinates(pdb,)

    lines = []
    d_data = {}

    for res_no in d_coordinates[chain].keys():
        if res_no not in res_range:
            continue
        res_name = d_coordinates[chain][res_no]['res_name']
        N = d_coordinates[chain][res_no]['atoms']['N']
        CA = d_coordinates[chain][res_no]['atoms']['CA']
        C = d_coordinates[chain][res_no]['atoms']['C']
        if res_no-1 in d_coordinates[chain].keys() and res_no+1 in d_coordinates[chain].keys():
            C_prev = d_coordinates[chain][res_no-1]['atoms']['C']
            phi = dihedral(C_prev,N,CA,C,)
            N_next = d_coordinates[chain][res_no+1]['atoms']['N']
            psi = dihedral(N,CA,C,N_next)
            s = '%s %s %s%s\n' %(phi,psi,res_name,res_no)
            lines += [s]

            d_data[res_no] = [phi,psi]

    return d_data, lines


def scatter_plot_2d_multiple_colors(
    l_prefixes, d_xtics=None,
    logarithmic=False, regression=False, errorbars=False,
    xlabel='',ylabel='RMSD',
    xmin = '', xmax = '',
    ymin = '', ymax = '',
    prefix_out='gnuplot',
    ):

    import os

    gnuplotsettings = []
    if regression == True:
        gnuplotsettings += [
            'f(x) = a*x+b\n',
            'fit f(x) "%s.gnuplotdata" via a,b\n' %(prefix),
            ]
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 36\n',
        'set output "%s.ps"\n' %(prefix_out),
        'set size 4,4\n', ## scale 400%
        'set size square\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
##            'set encoding iso_8859_1\n', ## postscript encoding for special characters
##            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
        'set xlabel "%s"\n' %(xlabel),
        'set ylabel "%s"\n' %(ylabel),
    ]
    if logarithmic == True:
        gnuplotsettings += ['set logscale x\n']
    if d_xtics:
        line_xtic = 'set xtics ('
        for xtic in d_xtics.keys():
            line_xtic += '"%s" %s, ' %(xtic, d_xtics[xtic])
        line_xtic = line_xtic[:-2]+')\n'
        gnuplotsettings += [
            line_xtic,
            'set xtics rotate\n',
        ]

    line_plot = 'plot '
    line_plot += '[%s:%s]' %(xmin,xmax)
    line_plot += '[%s:%s]' %(ymin,ymax)
    l_hexa = range(10)+['A','B','C','D','E','F',]
    for i in range(len(l_prefixes)):
        prefix = l_prefixes[i]
        h = 239.*i/len(l_prefixes)
        s = 240.
        l = 120.
        r,g,b = hsl2rgb(h,s,l,)
        r *= 255
        g *= 255
        b *= 255
        s_hexa = '#%s%s%s%s%s%s' %(
            l_hexa[int(r%16)], l_hexa[int((r-r%16)/16.)],
            l_hexa[int(g%16)], l_hexa[int((g-g%16)/16.)],
            l_hexa[int(b%16)], l_hexa[int((b-b%16)/16.)],
            )
        ## points
##        line_plot += '"%s.gnuplotdata" u 1:2 pt 7 ps 1 lc rgb "%s" t "%s"' %(prefix,s_hexa,prefix,)
        line_plot += '"%s.gnuplotdata" u 1:2 pt 7 ps 1 lc rgb "%s" t ""' %(prefix,s_hexa,)
        if errorbars == True:
            line_plot += ' w errorb'
        if regression == True:
            line_plot += ', f(x) lt 1 lc 0 lw 10 t ""'
        line_plot += ','

        ## labels
        if i == 0 or i == len(l_prefixes)-1:
            line_plot += '"%s.gnuplotdata" u 1:2:3 w labels font "Helvetica,16" textcolor rgb "%s" offset character 0,0.3 t ""' %(prefix,s_hexa,)
            line_plot += ','
    line_plot = line_plot[:-1]

    line_plot += ','
    s_hexa = '#000000'
    line_plot += '"%s.gnuplotdata" u 1:2:3 w labels font "Helvetica,16" textcolor rgb "%s" offset character 0,0.3 t ""' %(prefix,s_hexa,)
##    line_plot += '\n'

    ## connectors
    line_plot += ','
    line_plot += '"connectors.gnuplotdata" w l lc 0 lt 1 t ""'
    line_plot += '\n'

    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix_out),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(prefix_out))
    for prefix in l_prefixes:
        os.remove('%s.gnuplotdata' %(prefix))
##    os.remove('%s.gnuplotsettings' %(prefix_out))

    return


def hsl2rgb(h,s,l):

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

##    r = int(l_rgb[0]*255)
##    g = int(l_rgb[1]*255)
##    b = int(l_rgb[2]*255)
    r = l_rgb[0]
    g = l_rgb[1]
    b = l_rgb[2]

    return r,g,b


def dihedral(c1,c2,c3,c4):

    import numpy,math

    v1 = c2-c1
    v2 = c3-c2
    v3 = c4-c3

    angle = math.atan2(
        numpy.dot(
            math.sqrt(sum(v2*v2))*v1,
            cross(v2,v3),
            ),
        numpy.dot(
            cross(v1,v2),
            cross(v2,v3),
            ),
        )
    angle *= 180./math.pi

    return angle


def parse_backbone_coordinates(pdb,):

    import numpy

##    fd = open('/local/data/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb),'r')
    fd = open('%s.pdb' %(pdb),'r')
    lines = fd.readlines()
    fd.close()

    d_coordinates = {}

    for line in lines:

        record = line[:6].strip()

        if record == 'ATOM':

            atom_name = line[12:16].strip()
            if not atom_name in ['C','CA','O','N',]:
                continue

            altloc = line[16]
            if altloc not in [' ','A',]:
                continue
            
            iCode = line[26]
            altloc = line[16]
            if iCode != ' ':
                stop

            chain = line[21]
            if not chain in d_coordinates.keys():
                d_coordinates[chain] = {}

            res_no = int(line[22:26])
            if not res_no in d_coordinates[chain].keys():
                res_name = line[17:20].strip()
                d_coordinates[chain][res_no] = {'res_name':res_name,'atoms':{}}

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            coordinate = numpy.array([x, y, z])

            d_coordinates[chain][res_no]['atoms'][atom_name] = coordinate

    return d_coordinates


def cross(v1,v2):

    import numpy

    n = numpy.array([
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0],
        ])

    return n
    

if __name__ == '__main__':
    main()

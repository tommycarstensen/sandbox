def main():

    import os
    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import gnuplot

    fd = open('hewl_pdbs.txt','r')
    s = fd.read()
    fd.close()

    l_pdbs = s.split()
    l_pdbs = ['2lzt','2vb1','3lzt',]

    for pdb in l_pdbs:
        if not os.path.isfile('pdb/%s.pdb' %(pdb)):
            os.system('cp /data/remediated_pdb/%s/pdb%s.ent pdb/%s.pdb' %(pdb[1:3],pdb,pdb))

    d_shifts = {}

    for pdb in l_pdbs:

        d_shifts[pdb] = {}

        if not os.path.isfile('pred/%s.tab' %(pdb)):
            os.system('/home/people/tc/Desktop/SPARTA/sparta -in pdb/%s.pdb' %(pdb))
            os.system('mv pred/pred.tab pred/%s.tab' %(pdb))

        fd = open('pred/%s.tab' %(pdb),'r')
        lines = fd.readlines()
        fd.close()

        for i in range(11,len(lines)):
            line = lines[i]
            if len(line) == 1:
                continue
            res_no = int(line.split()[0])
            res_name = line.split()[1]
            atom_name = line.split()[2]
            shift = float(line.split()[4])
            if res_name == 'P':
                continue
            if res_no not in d_shifts[pdb].keys():
                d_shifts[pdb][res_no] = {'res_name':res_name,'shifts':{}}
            d_shifts[pdb][res_no]['shifts'][atom_name] = shift

    for pdb in d_shifts.keys():
        lines = []
        for res_no in d_shifts[pdb].keys():
            shift_N = d_shifts[pdb][res_no]['shifts']['N']
            shift_H = d_shifts[pdb][res_no]['shifts']['H']
            res_name = d_shifts[pdb][res_no]['res_name']
            lines += ['%s %s %s%s\n' %(shift_N,shift_H,res_name,res_no,)]

        fd = open('%s.gnuplotdata' %(pdb),'w')
        fd.writelines(lines)
        fd.close()
    
    scatter_plot_2d_multiple_colors(
        l_pdbs,
        xlabel = '{/Symbol w}_N',
        ylabel = '{/Symbol w}_H',
        )
        

    return


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
        line_plot += '"%s.gnuplotdata" u 1:2 pt 7 ps 1 lc rgb "%s" t "%s"' %(prefix,s_hexa,prefix,)
        if errorbars == True:
            line_plot += ' w errorb'
        if regression == True:
            line_plot += ', f(x) lt 1 lc 0 lw 10 t ""'
        line_plot += ','
        ## labels
        line_plot += '"%s.gnuplotdata" u 1:2:3 w labels font "Helvetica,16" textcolor rgb "%s" offset character 0,0.3 t ""' %(prefix,s_hexa,)
        line_plot += ','
    line_plot = line_plot[:-1]
##    line_plot += '"%s.gnuplotdata" u 1:2:3 w labels font "Helvetica,16" textcolor rgb "%s" offset character 0,0.3 t ""' %(prefix,s_hexa,)
    line_plot += '\n'
    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix_out),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(prefix_out))
    for prefix in l_prefixes:
        os.remove('%s.gnuplotdata' %(prefix))
    os.remove('%s.gnuplotsettings' %(prefix_out))

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


if __name__=='__main__':
    main()

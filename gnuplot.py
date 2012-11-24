import os

## Angstrom {\305}
## degree

def scatter_plot_2d(
    prefix,
    l1=None, l2=None,
    d_xtics=None,
    logarithmic=False,
    errorbars=False,
    regression=False, regression_data = None, regression_title = None,
    xlabel='',ylabel='RMSD',
    xmin = '', xmax = '',
    ymin = '', ymax = '',
    terminal = 'postscript'
    ):

    d_output = {
        'postscript':'ps',
        'png':'png',
        }
    d_terminal = {
        'postscript':'%s eps enhanced color "Helvetica" 48' %(terminal),
        'png':'png',
        }
    
    gnuplotsettings = []

##    regression = False ## tmp!!!

    if regression == True:
        gnuplotsettings += [
##            'set fit logfile "%s.log"\n' %(prefix),
            'f(x) = a*x+b\n',
            ]
        ## separate set of data to fit to? (e.g. if errorbars is main plot...)
        if regression_data:
            gnuplotsettings += [
                'fit f(x) "%s.gnuplotdata" via a,b\n' %(regression_data),
                ]
        else:
            gnuplotsettings += [
                'fit f(x) "%s.gnuplotdata" via a,b\n' %(prefix),
                ]

    ## write file to be plotted if data is provided
    if l1 != None and l2 != None:
        lines = ['%s %s\n' %(pair[0],pair[1],) for pair in zip(l1,l2)]
        fd = open('%s.gnuplotdata' %(prefix),'w')
        fd.writelines(lines)
        fd.close()
            
    gnuplotsettings += [
        'set terminal %s\n' %(d_terminal[terminal]),
        'set output "%s.%s"\n' %(prefix, d_output[terminal]),
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set encoding iso_8859_1\n', ## postscript encoding *necessary* for special characters (e.g. Angstrom)
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

    if regression_data:
        ps = 3
    else:
        ps = 1
    line_plot += '"%s.gnuplotdata" lc 0 lt 1 ps %s lw 3 pt 7 t ""' %(prefix,ps,)

    if errorbars == True:
        line_plot += ' w errorb'
    if regression == True:
        line_plot += ', f(x) lt 1 lc 0 lw 10'
        if regression_title:
            line_plot += ' t "%s"' %(regression_title)
        else:
            line_plot += ' t ""'
            
    line_plot += '\n'
    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

##    if os.path.isfile('wgnuplot.exe'):
##        os.system('wgnuplot.exe %s.gnuplotsettings' %(prefix))
##        os.system('convert.exe %s.ps %s.png' %(prefix, prefix))
##    else:
##        os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(prefix))
##        os.system('convert %s.ps %s.png' %(prefix, prefix))

##    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(prefix))
    os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(prefix))
    os.system('convert %s.ps %s.png' %(prefix, prefix))

    ## remove postscript if it was generated
    if os.path.isfile('%s.ps' %(prefix)):
        os.remove('%s.ps' %(prefix))
##    os.remove('%s.gnuplotdata' %(prefix))
##    os.remove('%s.gnuplotsettings' %(prefix))

    return


def contour_plot(
    fileprefix,
    lines_gnuplotdata,
    title='',xlabel='',ylabel='',zlabel='',
    x1=None,x2=None,
    y1=None,y2=None,
    z1=None,z2=None,
    d_xtics = None,
    d_ytics = None,
    bool_remove = True,
    size = 4,
    ):

    l_gnuplotdata = lines_gnuplotdata

    fd = open('%s.gnuplotdata' %(fileprefix), 'w')
    fd.writelines(l_gnuplotdata)
    fd.close()
    ## write gnuplot settings to txt file
    lines = ['set size square\n'] ## scale square
    lines += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(fileprefix),
        'set size %i,%i\n' %(size,size), ## scale 400%
        'set view map\n', ## change orientation of plot
        'set autoscale fix\n', ## scale axes
        'set style data pm3d\n', ## set by default?
        'set style function pm3d\n', ## set by default?
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
        'set title "%s"\n' %(title),
        ]
    if xlabel:
        lines += ['set xlabel "%s"\n' %(xlabel),]
    if ylabel:
        lines += ['set ylabel "%s"\n' %(ylabel),]
    if zlabel:
        lines += ['set cblabel "%s"\n' %(zlabel),]
    if z1 != None or z2 != None:
        ## set colorbox range
        line = 'set cbrange ['
        if z1: line += '%s' %(z1,z2)
        line += ':'
        if z2: line += '%s' %(z2)
        line += ']\n'
        lines += [line]
    lines += [
        'set palette model CMY rgbformulae 7,5,15\n',
        'set pm3d map corners2color c1\n', ## generate a 2D surface rather than 3D points
        ]

    if d_xtics:
        line_xtic = 'set xtics ('
        for xtic in d_xtics.keys():
            line_xtic += '"%s" %s, ' %(xtic, d_xtics[xtic])
        line_xtic = line_xtic[:-2]+')\n'
        lines += [
            line_xtic,
            'set xtics rotate\n',
##            'set xtics rotate by -90 offset 0,-1.5\n', ## tmp!!!
        ]

    if d_ytics:
        line = 'set ytics ('
        for tic in d_ytics.keys():
            line += '"%s" %s, ' %(tic, d_ytics[tic])
        line = line[:-2]+')\n'
        lines += [
            line,
##            'set xtics rotate by -90 offset 0,-1.5\n', ## tmp!!!
        ]

    line = 'splot '
    line += '['
    if x1: line += str(x1)
    line += ':'
    if x2: line += str(x2)
    line += ']'
    line += '['
    if y1: line += str(y1)
    line += ':'
    if y2: line += str(y2)
    line += ']'
    line += '"%s.gnuplotdata" title ""\n' %(fileprefix) ## splot gnuplot data file
    lines += [line]

    fd = open('%s.gnuplotsettings' %(fileprefix), 'w')
    fd.writelines(lines)
    fd.close()
    ## plot data with gnuplot splot
    os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(fileprefix))
    ## convert postscript to portable network graphics
    os.system('convert %s.ps %s.png' %(fileprefix, fileprefix))
    os.remove('%s.ps' %(fileprefix))
    if bool_remove == True:
        os.remove('%s.gnuplotdata' %(fileprefix))
        os.remove('%s.gnuplotsettings' %(fileprefix))
    
    return


def histogram2(
    fileprefix,
    l_data = None,
    x_min = None, x_max = None, x_step = None,
    xlabel = None, title = None,
    bool_remove = True,
    tic_step = 6,
    lines_extra = None,
    ):

    lines = ['%f\n' %(y) for y in l_data]
    fd = open('%s.gnuplotdata' %(fileprefix),'w')
    fd.writelines(lines)
    fd.close()

    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 24\n',
        'set output "%s.ps"\n' %(fileprefix),
        'set size 2,2\n', ## scale 400%
        'set encoding iso_8859_1\n',
        ]
    gnuplotsettings += [
        'min=%f\n' %(float(x_min)),
        'max=%f\n' %(float(x_max)),
        'width=%s	#interval width\n' %(x_step),
        '#function used to map a value to the intervals\n',
        'hist(x,width)=width*floor(x/width)+width/2.0\n',
        'set xrange [min:max]\n',
        'set yrange [0:]\n',
        '#to put an empty boundary around the\n',
        '#data inside an autoscaled graph.\n',
##        'set offset graph 0.05,0.05,0.05,0.0\n',
##        'set xtics min,(max-min)/%i,max\n' %(n_tics), ## min and max must be floats...
        'set xtics min,%f,max\n' %(tic_step), ## min and max must be floats...
        'set boxwidth width*0.9\n',
        'set style fill solid 0.5	#fillstyle\n',
        'set tics out nomirror\n',
        'set xlabel "%s"\n' %(xlabel),
        'set title "%s"\n' %(title),
##        'set ylabel "Frequency"\n',
        '#count and plot\n',
        'unset ytics\n',
        ]
    if lines_extra:
        gnuplotsettings += lines_extra
    gnuplotsettings += [
        'plot "%s.gnuplotdata" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"blue" notitle\n' %(fileprefix),
        ]

    ## write gnuplot settings
    fd = open('%s.gnuplotsettings' %(fileprefix),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
##    os.system('/software/bin/gnuplot gnuplot.settings')
    os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(fileprefix))
    ## convert postscript to portable network graphics
    if os.path.isfile('%s.png' %(fileprefix)):
        os.remove('%s.png' %(fileprefix))
    os.system('convert %s.ps %s.png' %(fileprefix,fileprefix))
    os.remove('%s.ps' %(fileprefix))
    if bool_remove == True:
        os.remove('%s.gnuplotdata' %(fileprefix))
        os.remove('%s.gnuplotsettings' %(fileprefix))

    return


def histogram(
    fileprefix,
    d_data = None, l_xtics = None,
    l_data = None,
    ylabel=None,xlabel=None,l_plotdatafiles=[],title=None,
    x_min = None, x_max = None, x_step = None,
    y_min = '', y_max = '',
    ):

    import math, os

    ##
    ## write data
    ##
    if l_data == None:
        yrange = []
        gnuplotdata = []
        for i in range(len(l_xtics)):
            res_name = l_xtics[i]
            for y in d_data[res_name]:
                gnuplotdata += ['%f %f\n' %(float(i),y)]
                yrange += [y]
        y_min = min(yrange)
        y_max = max(yrange)
        fd = open('%s.gnuplotdata' %(fileprefix),'w')
        fd.writelines(gnuplotdata)
        fd.close()
    elif l_data != None:
        d_count = {}
        for x in l_data:
            xbin = x_step*int(x/x_step)
            if not xbin in d_count.keys():
                d_count[xbin] = 0
            d_count[xbin] += 1
        gnuplotdata = []
        for xbin in range(int(x_min/x_step),int(x_max/x_step)+1,):
            xbin *= x_step
            if xbin in d_count.keys():
                count = d_count[xbin]
            else:
                count = 0
            gnuplotdata += ['%f %i\n' %(xbin,count,)]
        fd = open('%s.gnuplotdata' %(fileprefix),'w')
        fd.writelines(gnuplotdata)
        fd.close()

        if l_xtics != None:
            l_xtics = [xbin for xbin in range(int(x_min/x_step),int(x_max/x_step)+1,)]
    else:
        stop
        
    ##
    ## calculate statistics
    ##
    if l_xtics != None:
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
        'set terminal postscript eps enhanced color "Helvetica" 24\n',
        'set output "gnuplot.ps"\n',
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set style fill\n',
##        'set style histogram\n',
##        'set style data histograms\n',
        ]
    if title:
        gnuplotsettings += [
            'set title "%s"\n' %(title),
            ]
    if xlabel:
        gnuplotsettings += [
            'set xlabel "%s"\n' %(xlabel),
            ]
    gnuplotsettings += [
        'set ylabel "{/=48 %s}"\n' %(ylabel),
        ]
    if l_xtics != None:
        line_xtic = 'set xtics ('
        for xtic in l_xtics:
            line_xtic += '"%s" %s, ' %(xtic, l_xtics.index(xtic))
        line_xtic = line_xtic[:-2]+')\n'
        gnuplotsettings += [
            line_xtic,
            'set xtics rotate\n',
        ]
    gnuplotsettings += [
        'plot ',
##        '[-1:%i][%f:%f] "gnuplot.data" lt 0 ps 2 pt 2 t ""' %(len(l_xtics)+1, ymin, ymax),
        '[%f:%f][%s:%s] "%s.gnuplotdata" u 1:2 lt 0 t ""' %(
            x_min, x_max, y_min, y_max, fileprefix,
            ),
        ]
    if l_xtics:
        gnuplotsettings += [
            ', ',
            '"gnuplot.statistics" u 2 lt 1 lc 0 ps 0 pt 0 w errorb t ""',
            ]
    for plotdatafile in l_plotdatafiles:
        gnuplotsettings += [
            ', ',
            '"%s" lt 0 lc 1 ps 3 pt 7 t ""\n' %(plotdatafile),
            ]

    ## unset ytics
    if ylabel == None:
        gnuplotsettings += [
            'unset ytics\n',
            ]
            
    ## write gnuplot settings
    fd = open('gnuplot.settings','w')
    fd.writelines(gnuplotsettings)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
##    os.system('/software/bin/gnuplot gnuplot.settings')
    os.system('/usr/bin/gnuplot gnuplot.settings')
    ## convert postscript to portable network graphics
    os.system('convert gnuplot.ps %s.png' %(fileprefix))
    os.remove('gnuplot.ps')
    os.remove('%s.gnuplotdata' %(fileprefix))
    os.remove('gnuplot.settings')
    if l_xtics != None:
        os.remove('gnuplot.statistics')

    return

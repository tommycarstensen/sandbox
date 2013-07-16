#!/usr/bin/python3

## Tommy Carstensen, 2007-2013

import os, math, sys
import time

## Angstrom {\305}
## degree {\260} ?
## circumflex ?

def scatter_plot_2d(
    prefix,
    l1=None, l2=None,
    d_xtics=None,
    logarithmic=False,
    errorbars=False,
    regression=False, regression_data = None, regression_title = None,
    xlabel='',ylabel='',
    title=None,
    x_min = '', x_max = '',
    y_min = '', y_max = '',
    terminal = 'postscript',
    column1 = '1', column2 = '2',
    line_plot = None, ## manual plot line
    s_plot = None, ## manual def of stuff to be plotted
    bool_remove = True,
    lines_extra = None,
    prefix_out=None,
    bool_execute=True,
    bool_title_enhanced=True,
    path_gnuplot=None,
    bool_timestamp = False,
    ):

    if not path_gnuplot:
        path_gnuplot = '/usr/bin/gnuplot'

    if not prefix_out:
        prefix_out = prefix

    d_output = {
        'postscript':'ps',
        'png':'png',
        }
    d_terminal = {
        'postscript':'%s eps enhanced color "Helvetica" 48' %(terminal),
        'png':'png',
        }
    
    sett = []

##    regression = False ## tmp!!!

    if regression == True:
        sett += [
##            'set fit logfile "%s.log"\n' %(prefix),
            'f(x) = a*x+b\n',
            ]
        ## separate set of data to fit to? (e.g. if errorbars is main plot...)
        if regression_data:
            sett += [
                'fit f(x) "%s" via a,b\n' %(regression_data),
                ]
        else:
            sett += [
                'fit f(x) "%s" via a,b\n' %(prefix),
                ]

    ## write file to be plotted if data is provided
    if l1 != None and l2 != None and not os.path.isfile(prefix):
        lines = ['%s %s\n' %(pair[0],pair[1],) for pair in zip(l1,l2)]
        fd = open('%s' %(prefix),'w')
        fd.writelines(lines)
        fd.close()
            
    sett += [
        'set terminal %s\n' %(d_terminal[terminal]),
        'set output "%s.%s"\n' %(prefix_out, d_output[terminal]),
        'set size 4,4\n', ## scale 400%
        'set encoding iso_8859_1\n', ## postscript encoding *necessary* for special characters (e.g. Angstrom)
##        'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
        'set xlabel "%s"\n' %(xlabel),
        'set ylabel "%s"\n' %(ylabel),
    ]

    if lines_extra:
        sett += lines_extra

    if title:
        if bool_timestamp == True:
            title += '\\n%s' %(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))
        if bool_title_enhanced == True:
            sett += ['set title "%s"\n' %(title)]
        else:
            sett += ['set title "%s" noenhanced\n' %(title)]

    if logarithmic == True:
        sett += ['set logscale x\n']
    if d_xtics:
        line_xtic = 'set xtics ('
        for xtic in d_xtics.keys():
            line_xtic += '"%s" %s, ' %(xtic, d_xtics[xtic])
        line_xtic = line_xtic[:-2]+')\n'
        sett += [
            line_xtic,
            'set xtics rotate\n',
        ]

    if not line_plot:

        line_plot = 'plot '
        if not s_plot:
            line_plot += '[%s:%s]' %(x_min,x_max)
            line_plot += '[%s:%s]' %(y_min,y_max)
            line_plot += '"%s" u %s:%s ' %(prefix,column1,column2,)
        else:
            line_plot += s_plot

        if regression_data:
            ps = 3
        else:
            ps = 2
        line_plot += ' lc 0 lt 1 ps %s lw 3 pt 7 t ""' %(
            ps,
            )

        if errorbars == True:
            line_plot += ' w errorb'
        if regression == True:
            line_plot += ', f(x) lt 1 lc 1 lw 10'
            if regression_title:
                line_plot += ' t "%s"' %(regression_title)
            else:
                line_plot += ' t ""'
                
        line_plot += '\n'

    sett += [line_plot]

    fd = open('%s.plt' %(prefix_out),'w')
    fd.writelines(sett)
    fd.close()


    if bool_execute == True:
        cmd = '%s %s.plt' %(path_gnuplot,prefix_out)
        os.system(cmd)
        cmd = 'convert %s.ps %s.png' %(prefix_out, prefix_out,)
        os.system(cmd)

    ## remove postscript if it was generated
    if os.path.isfile('%s.ps' %(prefix_out)):
        os.remove('%s.ps' %(prefix_out))
    if bool_remove == True:
        if bool_execute == True:
            os.remove('%s.plt' %(prefix_out))

    return


def plot_and_convert(
    lines, prefix,
    x_min='',x_max='',y_min='',y_max='',
    xlabel='', ylabel='', title='',
    bool_remove = False,
    ):

    if not (
        'plot' in ''.join(lines)
        or
        'sp ' in ''.join(lines)
        ):
        sys.exit(0)

    plt = []

    plt += [
        'set terminal postscript eps enhanced color "Helvetica" 24\n',
        'set output "%s.ps"\n' %(prefix),
        'set size 5,5\n', ## scale (not text)
        'set encoding iso_8859_1\n',
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        ]
    plt += [
        'set boxwidth 1\n',
##        'set tics out nomirror\n',
        'set xlabel "{/=48 %s}"\n' %(xlabel),
        'set ylabel "{/=48 %s}"\n' %(ylabel),
##        'set title "%s" noenhanced\n' %(title),
        'set title "%s"\n' %(title),
        ]

## The `frequency` option makes the data monotonic in x; points with the same
## x-value are replaced by a single point having the summed y-values.  The
## resulting points are then connected by straight line segments.
## See also
## smooth.dem
    plt += lines

    ## write gnuplot settings
    fd = open('%s.plt' %(prefix),'w')
    fd.writelines(plt)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
##    os.system('/usr/bin/gnuplot %s.plt' %(prefix))
    os.system('/nfs/team149/Software/bin/gnuplot %s.plt' %(prefix))
    ## convert postscript to portable network graphics
    if os.path.isfile('%s.png' %(prefix)):
        os.remove('%s.png' %(prefix))
    os.system('convert %s.ps %s.png' %(prefix,prefix,))
    os.remove('%s.ps' %(prefix))
    if bool_remove == True:
        os.remove('%s.plt' %(prefix))

    return


def dict2dat(fn,d,x1,x2,y1,y2,bool_replace_diagonal,l_sequence,bool_log,step=1,):

    if l_sequence:
        l_x = l_sequence+[len(l_sequence)]
        l_y = l_sequence+[len(l_sequence)]
    else:
        l_x = list(range(x1,x2+step,step,))
        l_y = list(range(y1,y2+step,step,))

    with open(fn,'w') as f:
        for i in range(len(l_x)):
            x = l_x[i]
            for j in range(len(l_y)):
                y = l_y[j]
                if bool_replace_diagonal and x == y:
                    try:
                        z = float(d[x+1][y])
                    except:
                        z = float(d[x-1][y])
                else:
                    try:
                        z = float(d[x][y])
                    except:
                        z = 0
                if bool_log == True and z != 0:
                    z = math.log(z,10)
                if not l_sequence:
                    f.write('%s %s %s\n' %(x,y,z))
                else:
                    f.write('%s %s %s\n' %(i,j,z))
            f.write('\n')

    return


def contour_plot(
    fileprefix=None,
    path_dat=None, ## input=file
    d_dat=None, ## input=dict
    l_dat=None, ## input=lines
    title='',xlabel='',ylabel='',zlabel='',
    x1=None,x2=None,
    y1=None,y2=None,
    z1=None,z2=None,
    d_xtics = None,
    d_ytics = None,
    bool_remove = True,
    size = 4,
    line_splot = None,
    col1=0, col2=1, col3=2,
    bool_summed = False,
    s_sep=' ',
    bool_log = False,
    bool_replace_diagonal = False,
    l_sequence = [],
    ):

    ## at some point I should make it possible
    ## to just feed this function with 2 lists and zip the data automatically...

    ## this function is currently a bit of a mess...

    if not line_splot:
        if not fileprefix:
            fileprefix = 'contour'
        if d_dat:
            if not x1:
                x1 = min(list(d_dat.keys()))
            if not x2:
                x2 = max(list(d_dat.keys()))
            if not y1:
                y1 = min(list(d_dat[x1].keys()))
            if not y2:
                y2 = max(list(d_dat[x1].keys()))
            dict2dat(
                '%s.dat' %(fileprefix),d_dat,x1,x2,y1,y2,
                bool_replace_diagonal,l_sequence,
                bool_log,
                )
        elif path_dat:
            fileprefix = path_dat
            if bool_summed == False:
                d = {}
                with open(path_dat,'r') as file_dat:
                    for line_dat in file_dat:
                        l_dat = line_dat.strip().split(s_sep)
                        if l_dat == ['']: continue ## blank line
                        x = int(l_dat[col1])
                        y = int(l_dat[col2])
                        if x1 == None and x2 == None and y1 == None and y2 == None:
                            x1 = x
                            x2 = x
                            y1 = y
                            y2 = y
                        else:
                            if x < x1: x1 = x
                            elif x > x2: x2 = x
                            if y < y1: y1 = y
                            elif y > y2: y2 = y
                        v3 = float(l_dat[col3])
                        if bool_log == True and v3 != 0:
                            v3 = math.log(v3,10)
                        try:
                            d[x][y] += v3
                        except:
                            try:
                                d[x][y] = v3
                            except:
                                d[x] = {y:v3}
                dict2dat(
                    '%s.dat' %(path_dat),d,x1,x2,y1,y2,bool_replace_diagonal,
                    l_sequence,bool_log,
                    )
        elif l_dat:
            fd = open('%s.dat' %(fileprefix), 'w')
            fd.writelines(l_dat)
            fd.close()
        else:
            stopnoinput

    ## write gnuplot settings to txt file
    lines = ['set size square\n'] ## scale square
    lines += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(fileprefix),
        'set size %i,%i\n' %(size,size), ## scale 400%
        ]
    if not line_splot:
        lines += [
            'set view map\n', ## change orientation of plot
            'set style data pm3d\n', ## set by default?
            'set style function pm3d\n', ## set by default?
            ]
    lines += [
        'set autoscale fix\n', ## scale axes
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
        'set title "%s"\n' %(title),
        ]
    if xlabel:
        lines += ['set xlabel "%s"\n' %(xlabel),]
    if ylabel:
        lines += ['set ylabel "%s"\n' %(ylabel),]
    if zlabel:
        lines += ['set cblabel "%s"\n' %(zlabel),]
        lines += ['set zlabel "%s"\n' %(zlabel),]
    if z1 != None or z2 != None:
        ## set colorbox range
        line = 'set cbrange ['
        if z1: line += '%s' %(z1,z2)
        line += ':'
        if z2: line += '%s' %(z2)
        line += ']\n'
        lines += [line]
    if not line_splot:
        lines += [
##            'set palette model CMY rgbformulae 7,5,15\n',
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

    if line_splot:
        lines += [line_splot]
    else:
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
        line += '"%s.dat" title ""\n' %(fileprefix) ## splot gnuplot data file
        lines += [line]

    fd = open('%s.plt' %(fileprefix), 'w')
    fd.writelines(lines)
    fd.close()
    ## plot data with gnuplot splot
    os.system('/usr/bin/gnuplot %s.plt' %(fileprefix))
    ## convert postscript to portable network graphics
    print('convert')
    os.system('convert %s.ps %s.png' %(fileprefix, fileprefix))
    os.remove('%s.ps' %(fileprefix))
    if bool_remove == True:
        if line_splot == None:
            os.remove('%s.dat' %(fileprefix))
        os.remove('%s.plt' %(fileprefix))
    
    return


def histogram2(
    fileprefix,
    l_data = None,
    x_min = None, x_max = None,
    y_max = '',
    x_step = None, ## width of boxes
    xlabel = None, ylabel = None, title = None,
    bool_remove = True,
    tic_step = None, ## tics on axis
    s_plot = None,
    column = '$1',
    lines_extra = None,
    color = 'blue',
    prefix_out = None,
    bool_timestamp = False,
    ):

    if not prefix_out:
        prefix_out = fileprefix

    if l_data:
        lines = ['%f\n' %(y) for y in l_data]
        fd = open('%s' %(fileprefix),'w')
        fd.writelines(lines)
        fd.close()
    elif s_plot == None and not os.path.isfile(fileprefix):
        print('nothing to plot', fileprefix)
        return

    if bool_timestamp == True:
        title += '\\n%s' %(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()))

    sett = []
    sett += [
        'set terminal postscript eps enhanced color "Helvetica" 24\n',
        'set output "%s.ps"\n' %(prefix_out),
        'set size 2,2\n', ## scale 400%
        'set encoding iso_8859_1\n',
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        ]
    if x_min != None and x_max != None:
        sett += [
            'min=%f\n' %(float(x_min)),
            'max=%f\n' %(float(x_max)),
            'set xrange [min:max]\n',
            ]
        if tic_step:
            sett += ['set xtics min,%f,max\n' %(tic_step),] ## min and max must be floats...
    elif x_max != None:
        sett += ['set xrange [:%f]\n' %(float(x_max))]

    sett += [
        'width=%s	#interval width\n' %(x_step),
        '#function used to map a value to the intervals\n',
        'hist(x,width)=width*floor(x/width)+width/2.0\n',
        'set yrange [0:]\n',
        '#to put an empty boundary around the\n',
        '#data inside an autoscaled graph.\n',
##        'set offset graph 0.05,0.05,0.05,0.0\n',
##        'set xtics min,(max-min)/%i,max\n' %(n_tics), ## min and max must be floats...
        'set boxwidth width*0.9\n',
        'set style fill solid 0.5	#fillstyle\n',
        'set tics out nomirror\n',
        'set xlabel "{/=36 %s}"\n' %(xlabel),
##        'set title "%s" noenhanced\n' %(title),
        'set title "%s"\n' %(title),
##        'set ylabel "Frequency"\n',
        '#count and plot\n',
##        'unset ytics\n',
        ]
    if ylabel:
        sett += ['set ylabel "{/=36 %s}"\n' %(ylabel),]
    if lines_extra:
        sett += lines_extra

## The `frequency` option makes the data monotonic in x; points with the same
## x-value are replaced by a single point having the summed y-values.  The
## resulting points are then connected by straight line segments.
## See also
## smooth.dem
    if s_plot == None:
        sett += [
            'plot [:][:%s]"%s" u (hist(%s,width)):(1.0) smooth freq w boxes lc rgb"%s" notitle\n' %(
                y_max, fileprefix, column, color,
                ),
            ]
    else:
        sett += [s_plot]

    ## write gnuplot settings
    fd = open('%s.plt' %(prefix_out),'w')
    fd.writelines(sett)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/usr/bin/gnuplot %s.plt' %(prefix_out))
    ## convert postscript to portable network graphics
    if os.path.isfile('%s.png' %(prefix_out)):
        os.remove('%s.png' %(prefix_out))
    os.system('convert %s.ps %s.png' %(prefix_out,prefix_out,))
    os.remove('%s.ps' %(prefix_out))
    if bool_remove == True:
        os.remove('%s.plt' %(prefix_out))

    return


def histogram(
    fileprefix,
    d_data = None, l_xtics = None,
    l_data = None,
    ylabel=None,xlabel=None,l_plotdatafiles=[],title=None,
    x_min = None, x_max = None, x_step = None,
    y_min = '', y_max = '',
    ):

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
        fd = open('%s.dat' %(fileprefix),'w')
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
        fd = open('%s.dat' %(fileprefix),'w')
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
    sett = []
    sett += [
        'set terminal postscript eps enhanced color "Helvetica" 24\n',
        'set output "gnuplot.ps"\n',
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set style fill\n',
##        'set style histogram\n',
##        'set style data histograms\n',
        ]
    if title:
        sett += [
            'set title "%s"\n' %(title),
            ]
    if xlabel:
        sett += [
            'set xlabel "%s"\n' %(xlabel),
            ]
    sett += [
        'set ylabel "{/=48 %s}"\n' %(ylabel),
        ]
    if l_xtics != None:
        line_xtic = 'set xtics ('
        for xtic in l_xtics:
            line_xtic += '"%s" %s, ' %(xtic, l_xtics.index(xtic))
        line_xtic = line_xtic[:-2]+')\n'
        sett += [
            line_xtic,
            'set xtics rotate\n',
        ]
    sett += [
        'plot ',
##        '[-1:%i][%f:%f] "gnuplot.data" lt 0 ps 2 pt 2 t ""' %(len(l_xtics)+1, ymin, ymax),
        '[%f:%f][%s:%s] "%s.dat" u 1:2 lt 0 t ""' %(
            x_min, x_max, y_min, y_max, fileprefix,
            ),
        ]
    if l_xtics:
        sett += [
            ', ',
            '"gnuplot.statistics" u 2 lt 1 lc 0 ps 0 pt 0 w errorb t ""',
            ]
    for plotdatafile in l_plotdatafiles:
        sett += [
            ', ',
            '"%s" lt 0 lc 1 ps 3 pt 7 t ""\n' %(plotdatafile),
            ]
    sett += ['\n']

    ## unset ytics
    if ylabel == None:
        sett += [
            'unset ytics\n',
            ]
            
    ## write gnuplot settings
    fd = open('%s.plt' %(fileprefix),'w')
    fd.writelines(sett)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/usr/bin/gnuplot gnuplot.plt')
    ## convert postscript to portable network graphics
    os.system('convert gnuplot.ps %s.png' %(fileprefix))
    os.remove('gnuplot.ps')
    os.remove('%s.dat' %(fileprefix))
    os.remove('%s.plt' %(fileprefix))
    if l_xtics != None:
        os.remove('gnuplot.statistics')

    return

def bisection(r1,r2,A):

    bisection = .5
    d = r1
    ## the area of the intersection
    ## as a function of the distance between the centers of the circles
    ## is monotonic
    ## and thus a numerical solution can be found for d
    ## using a method of bisection
    while True:
        ## http://en.wikipedia.org/wiki/Circular_segment
        ## a) d=d1+d2
        ## b) d1**2=r1**2-(c/2)**2
        ## c) d2**2=r2**2-(c/2)**2
        ## d1**2-r1**2=d2**2-r2**2
        ## d1**2=d2**2+r1**2-r2**2
        ## d1**2+d1**2+2d1d2=d1**2+d2**2+2d1d2+r1**2-r2**2
        ## 2d1(d1+d2)=(d1+d2)**2+r1**2-r2**2
        ## d1 = (d**2+r1**2-r2**2)/2d
        d1 = (d**2+r1**2-r2**2)/(2*d)
        d2 = (d**2+r2**2-r1**2)/(2*d)
        alpha1 = 2*math.acos(d1/r1)
        alpha2 = 2*math.acos(d2/r2)
        ## A = circular sector - triangular portion
        ## A = (pi*r**2 * angle/2pi) - ((r**2 sinangle)/2))
        ## A = .5r**2(angle-sinangle)
        A_intersection = (
            ## area of circular segment 1
            .5*(r1**2)*(alpha1-math.sin(alpha1))
            +
            ## area of circular segment 2
            .5*(r2**2)*(alpha2-math.sin(alpha2))
            )
##            ## http://mathworld.wolfram.com/Circle-CircleIntersection.html
##            for x in xrange(1000000):
##                A_intersection = (
##                    (r1**2)*math.acos((d**2+r1**2-r2**2)/(2*d*r1))
##                    +
##                    (r2**2)*math.acos((d**2+r2**2-r1**2)/(2*d*r2))
##                    -
##                    .5*math.sqrt((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
##                    )
        diff = A_intersection - A
        if abs(diff) < 0.001:
            break
        elif diff > 0:
            d += bisection*r2
            bisection *= .5
            pass
        else:
            d -= bisection*r2
            bisection *= .5
            pass
        continue

    return d,d1,d2


def venn2(
    f1=None,f2=None,
    l1=None,l2=None,
    i1=None,i2=None,i3=None,
    text1='a',text2='b',
    suffix = 'venn2',
    verbose = False,
    ):

    ## Tommy Carstensen, November 2012

    '''
draw area proportional Venn diagram for 2 sets
set object circle only works with gnuplot 4.3 and higher
alternatively plot a parametric function with gnuplot 4.2 and lower...
the postscript terminal does not support transparency
'''

    if f1 != None:
        fd = open(f1,'r')
        l1 = fd.readlines()
        fd.close()
        fd = open(f2,'r')
        l2 = fd.readlines()
        fd.close()

    ## working with sets is definitely not the fastest method!!!
    if l1 != None:
        ## circles
        set_circle1 = set(l1)
        set_circle2 = set(l2)
        ## circle-circle intersections
        set_intersection12 = set_circle1&set_circle2
        ## lengths
        i10 = i1 = len(set_circle1)
        i01 = i2 = len(set_circle2)
        i11 = i3 = len(set_intersection12)

    ## radii of circles
    R1 = math.sqrt((i1)/math.pi)
    R2 = math.sqrt((i2)/math.pi)
    ## area of intersections between circles
    A_intersect12 = i3

    if verbose == True:
        print(R1, R2, R3)

    ##
    ## find distances between circle centers yielding correct area
    ## http://mathworld.wolfram.com/Circle-CircleIntersection.html
    ##
    ## distance between circle centers
    d,d1,d2 = bisection(R1,R2,A_intersect12,)

    ## height = r(1-cos(alpha/2))

    ## lengths of sides of triangle between centers of circles
    a = s110 = s12 = d

    ## coordinate centers of circles
    c1 = [0,0,]
    c2 = [a,0,]

    sett = []

    sett += [
        'set terminal pngcairo transparent enhanced size 1440,1080\n',
        ## transparency does not work for the postscript terminal...
        'set output "venn2_%s.png"\n' %(suffix),
        'set size 1,1\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        ## remove borders and tics
        'set noborder\n',
        'set noxtics\n',
        'set noytics\n',
        ## avoid elongated circles...
##        'set size square\n',
        'set size ratio -1\n', ## http://stackoverflow.com/questions/11138012/drawing-a-circle-of-radius-r-around-a-point
        ## add labels
        'set label 1 "%s (%i)" at %s, %s front nopoint tc rgbcolor "red" left font "Verdana,24" noenhanced\n' %(
##            text1, i1+i4+i5+i7, -0.7*R1, 0.7*R1,
            text1, i1, 'graph(0.05)','graph(0.95)',
            ),
        'set label 2 "%s (%i)" at %s, %s front nopoint tc rgb "green" right font "Verdana,24" noenhanced\n' %(
##            text2, i2+i4+i6+i7, a+0.7*R2, 0.7*R2,
            text2, i2, 'graph(0.95)','graph(0.95)',
            ),
        ]
    
    ## dual intersection (pos needs to be fixed...)
    sett += [
        'set label 4 "%i" at %f, %f front nopoint tc rgb "black" center font "Verdana,24" noenhanced\n' %(
            i3, d1, 0,
            ),
        ]

    sett += [
        'set obj 1 circle center 0,0 size %f fc rgb "red" fs transparent solid 0.5 noborder\n' %(R1,),
        'set obj 2 circle center %f,0 size %f fc rgb "green" fs transparent solid 0.5 noborder\n' %(c2[0],R2,),
        'set yrange [%f:%f]\n' %(0-max(R1,R2),0+max(R1,R2),),
        'plot [%f:%f][%f:%f] NaN notitle\n' %(
            c1[0]-R1,c2[0]+R2,
            -max(R1,R2),max(R1,R2),
            ),
        ]

    ## write gnuplot settings
    fd = open('%s.plt' %(suffix),'w')
    fd.writelines(sett)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/nfs/team149/Software/bin/gnuplot %s.plt' %(suffix))

    os.remove('%s.plt' %(suffix))

    return


def venn3(
    f1=None,f2=None,f3=None,
    l1=None,l2=None,l3=None,
    i1=None,i2=None,i3=None,i4=None,i5=None,i6=None,i7=None,
    text1='a',text2='b',text3='c',
    suffix = '',
    verbose=False,
    bool_labels=True,
    bool_sorted=False,
    ):

    ## Tommy Carstensen, August 2012, July 2013

    '''
draw area proportional Venn diagram for 3 sets
set object circle only works with gnuplot 4.3 and higher
alternatively plot a parametric function with gnuplot 4.2 and lower...
the postscript terminal does not support transparency
'''

    if f1 != None:
        for fn in (f1,f2,f3):
            cmd = 'cat %s' %(fn)
            if bool_sorted == False:
                cmd += ' | sort'
            cmd += ' > %s.sorted' %(fn)
            execmd(cmd)
        ## intersection set
        cmd = 'comm -12 %s.sorted %s.sorted > 11x' %(f1,f2)
        execmd(cmd)
        cmd = 'comm -23 %s.sorted %s.sorted > 10x' %(f1,f2)
        execmd(cmd)
        cmd = 'comm -13 %s.sorted %s.sorted > 01x' %(f1,f2)
        execmd(cmd)
        ## union set
        cmd = 'cat %s.sorted %s.sorted > 00x' %(f1,f2) ## 00x is a misleading name for a union set...
        execmd(cmd)
        i111 = int(os.popen('comm -12 %s.sorted 11x | wc -l' %(f3)).read())
        i110 = int(os.popen('comm -13 %s.sorted 11x | wc -l' %(f3)).read())
        i001 = int(os.popen('comm -23 %s.sorted 00x | wc -l' %(f3)).read())
        i010 = int(os.popen('comm -13 %s.sorted 01x | wc -l' %(f3)).read())
        i100 = int(os.popen('comm -13 %s.sorted 10x | wc -l' %(f3)).read())
        i011 = int(os.popen('comm -12 %s.sorted 01x | wc -l' %(f3)).read())
        i101 = int(os.popen('comm -12 %s.sorted 10x | wc -l' %(f3)).read())
##        fd = open(f1,'r')
##        l1 = fd.readlines()
##        fd.close()
##        fd = open(f2,'r')
##        l2 = fd.readlines()
##        fd.close()
##        fd = open(f3,'r')
##        l3 = fd.readlines()
##        fd.close()

    ## working with sets is definitely not the fastest method!!!
    if l1 != None:
        ## circles
        set_circle1 = set(l1)
        set_circle2 = set(l2)
        set_circle3 = set(l3)
        ## circle-circle intersections
        set_intersection12 = set_circle1&set_circle2
        set_intersection13 = set_circle1&set_circle3
        set_intersection23 = set_circle2&set_circle3
        set7 = set_intersection12&set_intersection13&set_intersection23
        ## subtract intersections
        set4 = set_intersection12-set7
        set5 = set_intersection13-set7
        set6 = set_intersection23-set7
        set1 = set_circle1-set4-set5-set7
        set2 = set_circle2-set4-set6-set7
        set3 = set_circle3-set5-set6-set7
        ## lengths
        i100 = i1 = len(set1)
        i010 = i2 = len(set2)
        i001 = i3 = len(set3)
        i12 = i110 = i4 = len(set4)
        i13 = i101 = i5 = len(set5)
        i23 = i011 = i6 = len(set6)
        i123 = i111 = i7 = len(set7)

    sum1 = i100+i110+i101+i111
    sum2 = i010+i110+i011+i111
    sum3 = i001+i101+i011+i111
    sum4 = i110
    sum5 = i101
    sum6 = i011
    sum7 = i111

    ## radii of circles
    R1 = math.sqrt((i100+i110+i101+i111)/math.pi)
    R2 = math.sqrt((i010+i110+i011+i111)/math.pi)
    R3 = math.sqrt((i001+i101+i011+i111)/math.pi)
    ## area of intersections between circles
    A_intersect12 = i110+i111
    A_intersect13 = i101+i111
    A_intersect23 = i011+i111

    if verbose == True:
        print(R1, R2, R3)

##    r_max = max(R1, R2, R3)
##    R1 /= r_max
##    R2 /= r_max
##    R3 /= r_max
##    A_intersect12 /= math.pi*r_max**2
##    A_intersect13 /= math.pi*r_max**2
##    A_intersect23 /= math.pi*r_max**2

    ##
    ## find distances between circle centers yielding correct area
    ## http://mathworld.wolfram.com/Circle-CircleIntersection.html
    ##
    ## list of distances between circle centers
    l_d = []
    for r1,r2,A in [
        [R1,R2,A_intersect12,],
        [R1,R3,A_intersect13,],
        [R2,R3,A_intersect23,],
        ]:
        d,d1,d2 = bisection(r1,r2,A)
        ## append distance to list of distance between circle centers
        l_d += [[d,d1,d2,]]
    del r1,r2

    ## height = r(1-cos(alpha/2))

    ## lengths of sides of triangle between centers of circles
    a = s110 = s12 = l_d[0][0]
    b = s101 = s13 = l_d[1][0]
    c = s011 = s23 = l_d[2][0]
    ## angles of triangle with corners at centers of circles
    C = angle_213 = math.acos((a**2+b**2-c**2)/(2*a*b))
    A = angle_132 = math.asin(a*math.sin(C)/c)
    B = angle_123 = math.asin(b*math.sin(C)/c)

    ## coordinate centers of circles
    c1 = [0,0,]
    c2 = [a,0,]
    ## coordinate of center of third circle
    x = b*math.cos(C)
    y = -b*math.sin(C)
    c3 = [x,y,]

    sett = []

    sett += [
        'set terminal pngcairo transparent enhanced size 1440,1080\n',
        ## transparency does not work for the postscript terminal...
        'set output "venn3_%s.png"\n' %(suffix),
        'set size 1,1\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        ## remove borders and tics
        'set noborder\n',
        'set noxtics\n',
        'set noytics\n',
        ## avoid elongated circles...
##        'set size square\n',
        'set size ratio -1\n', ## http://stackoverflow.com/questions/11138012/drawing-a-circle-of-radius-r-around-a-point
        ## add labels
        'set label 1 "%s (%i)" at %s, %s front nopoint tc rgbcolor "red" left font "Verdana,24" noenhanced\n' %(
##            text1, i1+i4+i5+i7, -0.7*R1, 0.7*R1,
            text1, sum1, 'graph(0.05)','graph(0.95)',
            ),
        'set label 2 "%s (%i)" at %s, %s front nopoint tc rgb "green" right font "Verdana,24" noenhanced\n' %(
##            text2, i2+i4+i6+i7, a+0.7*R2, 0.7*R2,
            text2, sum2, 'graph(0.95)','graph(0.95)',
            ),
##        'set label 3 "%s (%i)" at %f, %f front nopoint tc rgb "blue" center font "Verdana,24" noenhanced\n' %(
##            text3, i3+i5+i6+i7, c3[0], c3[1]-0.95*R3,
        'set label 3 "%s (%i)" at graph(0.95),graph(0.05) front nopoint tc rgb "blue" center font "Verdana,24" noenhanced\n' %(
            text3, sum3,
            ),
        ]

    if bool_labels == True:
        ## dual intersection (pos needs to be fixed...)
        if i110 > 0:
            sett += [
            'set label 4 "%i" at %f, %f front nopoint tc rgb "black" center font "Verdana,24" noenhanced\n' %(
    ##            i110, (R2*c1[0]+R1*c2[0])/(R1+R2), (R2*c1[1]+R1*c2[1])/(R1+R2),
                sum4, l_d[0][1], 0,
                ),
            ]
        if i101 > 0:
            alpha=2*math.acos(l_d[1][1]/R1)
            sett += [
                'set label 5 "%i" at %f, %f front nopoint tc rgb "black" center offset -2,1 font "Verdana,24" noenhanced\n' %(
    ##            i101, (R3*c1[0]+5*R1*c3[0])/(5*R1+R3), (R3*c1[1]+5*R1*c3[1])/(5*R1+R3),
    ##            i101, l_d[1][1]*math.cos(C), -l_d[1][1]*math.sin(C)
                sum5, R1*math.cos(C+alpha/4), -R1*math.sin(C+alpha/4),
                ),
            ]
        if i011 > 0:
            alpha=2*math.acos(l_d[2][1]/R2)
            sett += [
                'set label 6 "%i" at %f, %f front nopoint tc rgb "black" center offset 0,1.5 font "Verdana,24" noenhanced\n' %(
    ##            i011, (R3*c2[0]+5*R2*c3[0])/(5*R2+R3), (R3*c2[1]+5*R2*c3[1])/(5*R2+R3),
                sum6,
                l_d[0][0]+R2*math.cos(angle_123-alpha/4),
                0-R2*math.sin(angle_123-alpha/4),
                ),
            ]
        ## triple intersection (pos needs to be fixed...)
        if i111 > 0:
            alpha1 = 2*math.acos(l_d[1][1]/R1)
            alpha2 = 2*math.acos(l_d[2][1]/R2)
            x1=R1*math.cos(C-alpha1/2)
            y1=-R1*math.sin(C-alpha1/2)
            x2=l_d[0][0]+R2*math.cos(angle_123+alpha2/2)
            y2=0-R2*math.sin(angle_123+alpha2/2)
            x=.5*(x1+x2)
            y=.5*(y1+y2)
            sett += [
                'set label 7 "%i" at %f, %f front nopoint tc rgb "black" center offset -1,2 font "Verdana,24" noenhanced\n' %(
    ##            i111, (c1[0]+c2[0]+5*c3[0])/7., (c1[1]+c2[1]+5*c3[1])/7.,
                sum7, x,y,
                ),
    ##        'set arrow from %s,%s to %s,%s nohead lc 0 lw 3\n' %(
    ##            c1[0],c1[1],
    ##            x,y,
    ##            ),
    ##        'set arrow from %s,%s to %s,%s nohead lc 0 lw 6\n' %(
    ##            c1[0],c1[1],
    ##            x1,y1,
    ##            ),
    ##        'set arrow from %s,%s to %s,%s nohead lc 0 lw 9\n' %(
    ##            c1[0],c1[1],
    ##            x2,y2,
    ##            ),
            ]

##        'set label 11 "%s" at %f, %f front nopoint tc rgb "blue" center font "Verdana,24"\n' %(
##            i1, c3[0], c3[1]-0.95*R3,
##            ),

    sett += [
        'set obj 1 circle center 0,0 size %f fc rgb "red" fs transparent solid 0.5 noborder\n' %(R1,),
        'set obj 2 circle center %f,0 size %f fc rgb "green" fs transparent solid 0.5 noborder\n' %(c2[0],R2,),
        'set obj 3 circle center %f,%f size %f fc rgb "blue" fs transparent solid 0.5 noborder\n' %(c3[0],c3[1],R3,),
        'set yrange [%f:%f]\n' %(c3[1]-R3,0+max(R1,R2),),
        'plot [%f:%f][%f:%f] NaN notitle\n' %(
            c1[0]-R1,c2[0]+R2,
            min(c1[1]-R1,c3[1]-R3),0+max(R1,R2),
            ),
        ]

    ## write gnuplot settings
    fd = open('%s.plt' %(suffix),'w')
    fd.writelines(sett)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/nfs/team149/Software/bin/gnuplot %s.plt' %(suffix))

    os.remove('%s.plt' %(suffix))

##    print i110, i101, i011, i111

    return


def execmd(cmd):

    print(cmd)
    os.system(cmd)

    return


def venn5(
    t00001='',
    t00010='',t00011='',
    t00100='',t00101='',t00110='',t00111='',
    t01000='',t01001='',t01010='',t01011='',
    t01100='',t01101='',t01110='',t01111='',
    t10000='',t10001='',t10010='',t10011='',
    t10100='',t10101='',t10110='',t10111='',
    t11000='',t11001='',t11010='',t11011='',
    t11100='',t11101='',t11110='',t11111='',
    suffix = '',
    ):

    ## Tommy Carstensen, April 2013

    sett = []

    sett += [
        'set terminal pngcairo transparent enhanced size 1440,1080\n',
        ## transparency does not work for the postscript terminal...
        'set output "venn5_%s.png"\n' %(suffix),
        'set size 1,1\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        ## remove borders and tics
        'set noborder\n',
        'set noxtics\n',
        'set noytics\n',
        ## avoid elongated circles...
##        'set size square\n',
        'set size ratio -1\n', ## http://stackoverflow.com/questions/11138012/drawing-a-circle-of-radius-r-around-a-point
        ]

    sett += [
        'set obj 1 ellipse center 0.0,0.0 size 2.0,5.0 angle 0. fc rgb "red" fs transparent solid 0.5 noborder\n',
        'set obj 2 ellipse center 0.5,-1.0 size 2.0,5.0 angle 72. fc rgb "green" fs transparent solid 0.5 noborder\n',
        'set obj 3 ellipse center -0.25,0.0 size 2.0,5.0 angle 144. fc rgb "blue" fs transparent solid 0.5 noborder\n',
        'set obj 4 ellipse center -0.5,-0.25 size 2.0,5.0 angle 216. fc rgb "orange" fs transparent solid 0.5 noborder\n',
        'set obj 5 ellipse center -0.25,0.25 size 2.0,5.0 angle 288. fc rgb "violet" fs transparent solid 0.5 noborder\n',
        ]

    sett += [
        'set label 1 "1" at -1.5,1.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n',
        'set label 2 "2" at -2.5,3.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n',
        'set label 2 "3" at -0.5,2.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n',
        ]

    sett += [
        'plot [%f:%f][%f:%f] NaN notitle\n' %(-2.5,3.5,-2.5,3.5,),
        ]

    ## write gnuplot settings
    fd = open('%s.plt' %(suffix),'w')
    fd.writelines(sett)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/nfs/team149/Software/bin/gnuplot %s.plt' %(suffix))
    os.remove('%s.plt' %(suffix))

    return


def venn4(
    i0001=None,i0010=None,i0011=None,i0100=None,i0101=None,
    i0110=None,i0111=None,i1000=None,i1001=None,i1010=None,
    i1011=None,i1100=None,i1101=None,i1110=None,i1111=None,
    fn1=None,fn2=None,fn3=None,fn4=None,
    fn_intersection=None,
    text1='a',text2='b',text3='c',text4='d',
    text0001='',text0010='',text0011='',text0100='',text0101='',
    text0110='',text0111='',text1000='',text1001='',text1010='',
    text1011='',text1100='',text1101='',text1110='',text1111='',
    suffix = '',
    verbose=False,
    bool_remove=True,
    bool_percentage=False,
    ## are files already sorted?
    bool_sorted=False,
    i_column = 0,
    rgb1 = 'red', rgb2 = 'green', rgb3 = 'blue', rgb4 = 'gray',
    ):

    ## Tommy Carstensen, October-December 2012

    if fn1:
        l_fn = [fn1,fn2,fn3,fn4,]

        for fn in l_fn:
            if not os.path.isfile(fn):
                print('file does not exist:', fn)
                sys.exit()
        
        if fn_intersection: l_fn += [fn_intersection]
        for fn in l_fn:
            cmd = 'cat %s' %(fn)
            if i_column != 0:
                cmd += " | awk '{print $%i}'" %(i_column)
            if bool_sorted == False:
                cmd += ' | sort'
            cmd += ' > %s.sorted' %(fn)
            execmd(cmd)
        ## some additional file with which to take the intersection of all sets
        if fn_intersection:
            for fn in [fn1,fn2,fn3,fn4,]:
                execmd('comm -12 %s %s > %s.comm'  %(fn,fn_intersection,fn,))
                execmd('mv %s.comm %s.sorted' %(fn,fn))
        execmd('cat %s.sorted %s.sorted | sort -u > 00xx' %(fn1,fn2))
        execmd('comm -12 %s.sorted %s.sorted > 11xx' %(fn1,fn2))
        execmd('comm -23 %s.sorted %s.sorted > 10xx' %(fn1,fn2))
        execmd('comm -13 %s.sorted %s.sorted > 01xx' %(fn1,fn2))
        execmd('cat %s.sorted %s.sorted | sort -u > xx00' %(fn3,fn4))
        execmd('comm -12 %s.sorted %s.sorted > xx11' %(fn3,fn4))
        execmd('comm -23 %s.sorted %s.sorted > xx10' %(fn3,fn4))
        execmd('comm -13 %s.sorted %s.sorted > xx01' %(fn3,fn4))

        i0001 = int(os.popen('comm -13 00xx xx01 | wc -l').read())
        i0010 = int(os.popen('comm -13 00xx xx10 | wc -l').read())
        i0011 = int(os.popen('comm -13 00xx xx11 | wc -l').read())

        i0100 = int(os.popen('comm -23 01xx xx00 | wc -l').read())
        i0101 = int(os.popen('comm -12 01xx xx01 | wc -l').read())
        i0110 = int(os.popen('comm -12 01xx xx10 | wc -l').read())
        i0111 = int(os.popen('comm -12 01xx xx11 | wc -l').read())

        i1000 = int(os.popen('comm -23 10xx xx00 | wc -l').read())
        i1001 = int(os.popen('comm -12 10xx xx01 | wc -l').read())
        i1010 = int(os.popen('comm -12 10xx xx10 | wc -l').read())
        i1011 = int(os.popen('comm -12 10xx xx11 | wc -l').read())

        i1100 = int(os.popen('comm -23 11xx xx00 | wc -l').read())
        i1101 = int(os.popen('comm -12 11xx xx01 | wc -l').read())
        i1110 = int(os.popen('comm -12 11xx xx10 | wc -l').read())
        i1111 = int(os.popen('comm -12 11xx xx11 | wc -l').read())

        if bool_remove == True:
            for fn in [
                'xx01','xx10','xx11','xx00','01xx','10xx','11xx','00xx',]:
                os.remove(fn)
            for fn in l_fn:
                os.remove('%s.sorted' %(fn))

    if text0001=='':
        text0001 = str(i0001)
        text0010 = str(i0010)
        text0011 = str(i0011)
        text0100 = str(i0100)
        text0101 = str(i0101)
        text0110 = str(i0110)
        text0111 = str(i0111)
        text1000 = str(i1000)
        text1001 = str(i1001)
        text1010 = str(i1010)
        text1011 = str(i1011)
        text1100 = str(i1100)
        text1101 = str(i1101)
        text1110 = str(i1110)
        text1111 = str(i1111)

##    print i0001,i0010,i0011,i0100,i0101,i0110,i0111,i1000,i1001,i1010,i1011,i1100,i1101,i1110,i1111
##    print sum([i1000,i1001,i1010,i1011,i1100,i1101,i1110,i1111])
##    print sum([i0100,i0101,i0110,i0111,i1100,i1101,i1110,i1111])
##    print sum([i0010,i0011,i0110,i0111,i1010,i1011,i1110,i1111])
##    print sum([i0001,i0011,i0101,i0111,i1001,i1011,i1101,i1111])

    sum_union = sum([i0001,i0010,i0011,i0100,i0101,i0110,i0111,i1000,i1001,i1010,i1011,i1100,i1101,i1110,i1111])

    sett = []

    sett += [
        'set terminal pngcairo transparent enhanced size 1440,1080\n',
        ## transparency does not work for the postscript terminal...
        'set output "venn4_%s.png"\n' %(suffix),
        'set size 1,1\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        ## remove borders and tics
        'set noborder\n',
        'set noxtics\n',
        'set noytics\n',
        ## avoid elongated circles...
##        'set size square\n',
        'set size ratio -1\n', ## http://stackoverflow.com/questions/11138012/drawing-a-circle-of-radius-r-around-a-point
        ]

    if text1 and text2 and text3 and text4:

        sum4 = sum([i0001,i0011,i0101,i0111,i1001,i1011,i1101,i1111])
        sum3 = sum([i0010+i0011+i0110+i0111+i1010+i1011+i1110+i1111])
        sum2 = sum([i0100+i0101+i0110+i0111+i1100+i1101+i1110+i1111])
        sum1 = sum([i1000,i1001,i1010,i1011,i1100,i1101,i1110,i1111])
        l = [[text4,sum4,],[text3,sum3,],[text2,sum2,],[text1,sum1,],]
        for i in range(4):
            text = '%s (%i)' %(l[i][0],l[i][1],)
            if bool_percentage == True:
                percent = int(round(100*l[i][1]/float(sum_union)))
                text = text[:-1]+', %s%%)' %(percent)
            l[i][0] = text
        text4 = l[0][0]
        text3 = l[1][0]
        text2 = l[2][0]
        text1 = l[3][0]

        sett += [

        ## add labels
        'set label 4 "%s" at %s,%s front nopoint tc rgbcolor "%s" center font "Verdana,24" noenhanced\n' %(
            text4, -1.00,3.50, rgb4,
            ),
        'set label 3 "%s" at %s,%s front nopoint tc rgb "%s" center font "Verdana,24" noenhanced\n' %(
##            text2, i2+i4+i6+i7, a+0.7*R2, 0.7*R2,
            text3, 0.00,3.25, rgb3,
            ),
##        'set label 3 "%s (%i)" at %f, %f front nopoint tc rgb "blue" center font "Verdana,24" noenhanced\n' %(
##            text3, i3+i5+i6+i7, c3[0], c3[1]-0.95*R3,
        'set label 2 "%s" at %s,%s front nopoint tc rgb "%s" center font "Verdana,24" noenhanced\n' %(
            text2, 1.75,3.75, rgb2,
            ),
        'set label 1 "%s" at %s,%s front nopoint tc rgb "%s" center font "Verdana,24" noenhanced\n' %(
            text1, 3.00,4.00, rgb1,
            ),
        ]

    l_labels = [
        text0001,text0011,text0010,text0101,text0111,
        text1111,text1100,text0100,text1000,text1110,
        text1011,text1001,text0110,text1010,text1101,
        ]
    if bool_percentage == True:
        for i in range(15):
            l_labels[i] = '%s (%s%%)' %(l_labels[i],int(round(100*float(l_labels[i])/sum_union)))
    sett += [
        'set label 5 "%s" at -1.5,1.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[0]),
        'set label 6 "%s" at -0.75,1.75 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[1]),
        'set label 7 "%s" at -0.25,2.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[2]),
        'set label 8 "%s" at -0.5,-0.25 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[3]),
        'set label 9 "%s" at 0,1 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[4]),
        'set label 10 "%s" at 1,0 front nopoint tc rgb "black" center font "Verdana Bold,13" noenhanced\n' %(l_labels[5]),
        'set label 11 "%s" at 2.75,1.75 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[6]),
        'set label 12 "%s" at 2.25,2.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[7]),
        'set label 13 "%s" at 3.5,1.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[8]),
        'set label 14 "%s" at 2,1 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[9]),
        'set label 15 "%s" at 1.75,-0.75 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[10]),
        'set label 16 "%s" at 1,-1.5 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[11]),
        'set label 17 "%s" at 1,2.0 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[12]),
        'set label 18 "%s" at 2.5,-0.25 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[13]),
        'set label 19 "%s" at 0.25,-0.75 front nopoint tc rgb "black" center font "Verdana,12" noenhanced\n' %(l_labels[14]),
        ]

    sett += [
        'set label 20 "union (%i)" at 0.00,4.25 front nopoint tc rgb "black" center font "Verdana,18" noenhanced\n' %(sum_union),
        ]

    sett += [
        'set obj 1 ellipse center 0.0,0.0 size 3.0,5.0 angle 45. fc rgb "%s" fs transparent solid 0.5 noborder\n' %(rgb4),
        'set obj 2 ellipse center 1.0,1.0 size 3.0,5.0 angle 45. fc rgb "%s" fs transparent solid 0.5 noborder\n' %(rgb3),
        'set obj 3 ellipse center 1.0,1.0 size 3.0,5.0 angle -45. fc rgb "%s" fs transparent solid 0.5 noborder\n' %(rgb2),
        'set obj 4 ellipse center 2.0,0.0 size 3.0,5.0 angle -45. fc rgb "%s" fs transparent solid 0.5 noborder\n' %(rgb1),
        ]

    sett += [
        'plot [%f:%f][%f:%f] NaN notitle\n' %(-2.5,4.5,-2.5,4.5,),
        ]

    ## write gnuplot settings
    fd = open('%s.plt' %(suffix),'w')
    fd.writelines(sett)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/nfs/team149/Software/bin/gnuplot %s.plt' %(suffix))
    os.remove('%s.plt' %(suffix))

    return


if __name__ == '__main__':

    venn2(
        l1 = [1,2,],
        l2 = [2,3,],
        )
    stop

    venn4(
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
        text1='quad positions',text2='quad rsIDs',
        text3='octo rsIDs',text4='octo positions',
        text0001='N/A',
        text0010='N/A',text0011='76993\\nGA008510',
        text0100='N/A',text0101='N/A',text0110='N/A',text0111='1657255\\nrs8187692',
        text1000='N/A',text1001='N/A',text1010='N/A',text1011='470\\nrs10026985',
        text1100='15129\\n200610-1',text1101='470\\nrs10026985',text1110='1657255\\nkgp1000003',text1111='695287\\nGA008524',
        )
    stop

    venn3(
        l1=[
            1,2,3,4,
            8,9,10,11,12,13,
            ],
        l2=[
            1,2,5,6,
            8,9,10,
            #11,12,13,
            ],
        l3=[1,3,5,7,],
        )

#!/software/bin/python
#Tommy Carstensen, University College Dublin, 2011

import os
import sys
sys.path.append('/home/tc/svn/tc_sandbox/misc/')
import statistics, gnuplot

## HEWL
s_exclude = '2VB1, 3LZT, 1IEE, 4LZT, 3AJN, 1V7S, 1V7T, 2D4J, 2D4K, 2D4I, 1WTN, 2Z18, 2Z12, 3AGI, 2HUB, 2Z19, 1QIO, 2HU3, 2ZYP, 3D9A, 2XBR, 2D6B, 2BLY, 2BLX, 2XBS, 1WTM, 193L, 3P64, 2CGI, 3N9E, 2BPU, 3N9A, 2IHL, 194L, 3A8Z, 3P66, 2HTX, 2PC2, 2X0A, 3EXD, 3N9C, 2FBB, 3A90, 1LZB, 3AGH, 1RJC, 1SQ2, 2F2N, 1LKR, 3A92, 1VDQ, 1VAU, 2W1L, 2C8P, 1HF4, 3E3D, 2C8O, 1XFP, 3A95, 1LZA, 1A2Y, 2HU1, 1ZVH, 3A91, 1AKI, 3AGG, 3A93, 2F30, 1IR8, 3KAM, 1YIL, 3A94, 1GPQ, 1T3P, 1VDS, 1H6M, 1HEL, 1LZ8, 2ZQ3, 1DPW, 3P4Z, 3QY4, 2XJW, 1LSM, 2A7D, 1B2K, 2F4G, 1VAT, 1VDT, 1DPX, 1HEM, 1LCN, 1VDP, 2B5Z, 1W6Z, 3A34, 1LZC, 3P68, 1KXX, 1ZVY, 1HEN, 1UIH, 1RFP, 2W1X, 2W1Y, 1LZE, 1HEP, 1BGI, 1IOS, 1UIB, 1H87, 1IOT, 3IJU, 3EMS, 2WAR, 1HEQ, 1N4F, 1LSE, 1LZG, 1LMA, 1HEO, 1IOR, 1GWD, 1HER, 1LYS, 1YQV, 1LSC, 1YIK, 1T6V, 3OK0, 1FLQ, 1IOQ, 1FN5, 1HC0, 1LSB, 1LSD, 1LZD, 1UIA, 2H9K, 3RT5, 1FLW, 2YVB, 1FLU, 1F10, 2H9J, 1BWI, 2LZT, 3OJP, 1KXY, 1BWJ, 132L, 2AUB, 1LSA, 2W1M, 2Q0M, 1LZ9, 2G4P, 1LSN, 3A96, 2G4Q, 1LSF, 1FLY, 1G7J, 1BVX, 1AT5, 2A7F, 1BWH, 1AT6, 1VFB, 3A6C, 2DQD, 1DKK, 2I25, 3SP3, 1G7I, 3A67, 5LYM, 1RI8, 1UIF, 1J1X, 2DQC, 1J1O, 1J1P, 1YKZ, 1AZF, 3RNX, 3IJV, 3A6B, 2DQJ, 1KIQ, 1LJF, 1LJH, 2XTH, 1IR9, 6LYT, 2F4A, 1VED, 1IR7, 2LYM, 3LYM, 1UC0, 1JJ3, 5LYT, 4LYT, 3RZ4, 3RW8, 1UIG, 1YKX, 1UID, 1UIC, 1UIE, 1UUZ, 1HEW, 1KXW, 1B0D, 1YKY, 1RCM, 1JJ1, 2LYO, 1NBY, 3LYO, 1F0W, 1G7H, 2I6Z, 1JIY, 1JIS, 1HSX, 1JIT, 1LYO, 1KIR, 2DQE, 1JJ0, 1G7M, 2ZQ4, 1UCO, 1LPI, 3M18, 1LJ4, 1LJE, 3LYT, 1DKJ, 1QTK, 1C10, 1LJ3, 1LJG, 1NDG, 2CDS, 3A3R, 1YL1, 3AW6, 1LJI, 1LJJ, 1KIP, 4LYM, 3AW7, 2EIZ, 1NBZ, 1PS5, 4LYO, 1Z55, 2YBI, 2YDG, 1LJK, 1HSW, 1YL0, 2YBH, 1UA6, 1IC7, 2YBJ, 1JTT, 1P2C, 2YBL, 1G7L, 1XEI, 1JPO, 3A3Q, 2YBM, 1ZV5, 1DQJ, 1XEJ, 2DQI, 2YBN, 2EKS, 1IC5, 1NDM, 3P65, 1XGP, 3T6U, 3M3U, 1LZT, 1XGU, 1XGQ, 2ZNX, 1XEK, 2EPE, 3G3A, 2DQG, 2DQH, 1XGR, 2D91, 1XGT, 1IC4, 1MLC, 1FDL, 1C08, 2YSS, 2I26, 1BQL, 1JTO, 1MEL, 2DQF, 2ZNW, 3G3B, 9LYZ, 1BVK, 3F6Z, 1ZMY, 3HFM, 1BHZ, 2HSO, 2HS7, 2HS9, 2A6U, 1SF4, 1SF6, 1SF7, 1SFB, 1SFG, 1GXX, 1GXV, 1JA2, 1JA4, 1JA6, 1JA7, 1IO5, 1E8L, 1LZN, 1LKS, 2IFF, 1LZH, 2LZH, 8LYZ, 7LYZ, 1LYZ, 2LYZ, 3LYZ, 4LYZ, 5LYZ, 6LYZ'
## T4L
s_exclude += '1SX7,1SX2,1SWY,1SWZ,3F8V,3FA0,3F9L,3FAD,3DKE,3HH6,3HH5,3HH3,3HH4,3HTG,3HUK,2RBO,2RBN,2Q9D,2IGC,2NTG,3SBB,3HTD,3HUA,3HUQ,2RBR,3HU9,3GUI,2RB2,2RBP,1P36,1LW9,1PQM,1P7S,3K2R,1P2L,1P6Y,2RBS,3FI5,1XEP,1P3N,1P37,3GUP,2OU9,3HT6,119L,3GUN,1P2R,161L,1L58,1L92,217L,1L83,1L62,1L19,3HT8,1ZUR,1L47,1PQI,244L,1L46,1L17,1L18,1L23,1L33,3LZM,131L,221L,1L66,1L60,1L24,1L52,1P64,241L,237L,1L74,219L,1G0M,128L,110L,1L65,1L30,1L22,1LYH,139L,130L,1L45,1L36,1L03,3C8Q,1L63,129L,1L68,1L44,4LZM,1PQD,1L48,1L32,3DMV,1G07,3GUJ,3CDT,165L,2RAZ,166L,1P46,1L90,211L,138L,1L10,1L59,1L41,1L06,1L13,1L15,1KNI,1L12,1L14,212L,3C7Z,247L,240L,243L,3C82,1L07,1L16,1L94,1KW5,1L29,163L,118L,1D9W,1L09,159L,1L93,1L91,2A4T,1L11,156L,122L,108L,1G0Q,164L,123L,2O79,173L,1L02,1L04,1L05,232L,1L26,1L01,158L,160L,114L,2RBQ,3DNA,162L,1G0L,1DYE,115L,206L,1L08,1ZYT,102L,238L,126L,111L,1LYE,1LYG,7LZM,3HT7,188L,1L38,3CDQ,1PQO,1CVK,1L42,1L43,1L35,1L27,1DYB,181L,1L61,127L,1LPY,1LYJ,5LZM,125L,1QSB,255L,239L,245L,246L,250L,229L,6LZM,157L,2RB1,2LZM,1G0J,1G0P,112L,1LWG,242L,184L,120L,113L,1L98,180L,256L,1LYF,1CV3,186L,220L,3C8S,1C69,1QT7,1NHB,107L,1L86,1L87,1L80,137L,3DN8,1L49,1L37,1L39,1L40,3CDR,1QT5,1QUD,182L,1L88,1L20,185L,187L,1G0K,1CU2,1QT3,1G06,1QS9,3C8R,1G1W,258L,1L56,1L31,226L,109L,1CV4,2CUU,183L,3CDV,3C7W,1CV5,199L,155L,2HUL,3HTF,224L,249L,1L51,1L25,1L50,146L,228L,3DN2,3HU8,1QSQ,1L67,1L28,3DN3,1L72,1L73,3DMX,1G0G,1L54,1D2W,1T8G,1L71,1KW7,235L,3DN6,3G3X,1ZWN,1L00,1CUP,3L2X,3DN4,254L,222L,2RAY,1L21,1CV6,1DYA,3HWL,248L,236L,225L,2OU8,260L,195L,1LGW,1QTB,1L84,1L69,3DN1,1L53,1CX7,1L55,1L57,1G1V,1C6Q,3DN0,1C6H,1C6I,1QT6,2NTH,3C83,1C6L,1QUH,1L79,1L76,223L,148L,1L75,1EPY,1QT8,233L,192L,1DYF,1C6J,210L,1LLH,1LI3,234L,1L89,1C6E,1C6K,1QUO,3HTB,1L70,1I6S,257L,1C6G,3C81,1OWY,1LGX,1L34,2O4W,1L64,2OTY,227L,3CDO,1KY0,1C6P,141L,1LYI,1L0J,2RB0,1L99,103L,214L,1D3J,1PQJ,1LGU,230L,198L,143L,191L,1P56,1B6I,2L78,1CX6,1OWZ,1C6F,1L0K,142L,1T6H,3L64,1QUG,1QTH,1L85,200L,2OEA,172L,1C6C,253L,1L95,1L81,1C60,1C6T,1QTZ,259L,2OE9,145L,147L,1CUQ,1D3N,1TLA,1L96,1L97,1LI2,1C64,152L,1C61,1LI6,1C65,1C6D,190L,1D3F,1CU5,1L77,1C63,1OVH,144L,3G3V,140L,1DYC,3HT9,1KY1,2OU0,1D2Y,1CV1,3DMZ,1PQK,2HUK,1DYD,201L,1CV0,2B6T,2OE7,1CU3,1C6M,1D3M,3GUK,215L,1QT4,252L,197L,2OE4,2B6X,1CU6,1OV7,1C6A,2B75,1LYD,1JTM,205L,1DYG,1KS3,2B72,1CTW,3C80,1OVJ,1L82,3C7Y,2B74,175L,2B73,218L,2B6W,1OVK,1OV5,1C6B,213L,1SSW,216L,1CU0,1C66,196L,167L,2Q9E,176L,1QS5,2B6Z,1C6N,151L,3G3W,1C62,2B70,2OTZ,1C67,1LWK,2B6Y,1QTC,1QTV,3GUO,150L,3GUL,1SSY,1T8F,174L,251L,1C68,1QTD,3GUM,231L,3SB5,3SB9,1JTN,189L,2QAR,149L,209L,171L,2HUM,104L,3SB7,178L,169L,177L,170L,3SB8,1P5C,168L,1JQU,3SB6,3SBA,2B7X,3JR6,2LC9,2LCB'

def main():

    d_units = {'CA':'{\305}','heavy':'\305','phipsi':'/Symbol \260','chi1':'/Symbol \260'} ## octals
    d_octals = {'CA':'C_{/Symbol a}','heavy':'heavy atom','phipsi':'{/Symbol f}{/Symbol y}','chi1':'{/Symbol c}_1'} ## octals
    l_radii = [10,20,40,]
    d_columns = {'CA':9,'heavy':13,'phipsi':17,'chi1':21,}
    l_distances = [0,10,20,40,]
##    suffix = 'sameSG_sameauth'
    suffix = 'sameSG'

    d_rmsds = {}
    for radius in l_radii:
        d_rmsds[radius] = {}
        for key in d_columns.keys():
            d_rmsds[radius][key] = {}
            for dist_from_mut in range(0,83+1):
                d_rmsds[radius][key][dist_from_mut] = []

    ## the function sphere in quakes.py wrote this file when using the -singlemutants flag
    fn = 'sphere/out_%s_allresidues.txt' %(suffix)
    fd = open(fn,'r')
    lines = fd.readlines()
    fd.close()
    print 'read', fn

    print 'looping over lines'
    for line in lines:
        l = line.split()
        ## exclude HEWL and T4L
        if l[0].upper() in s_exclude or l[1].upper() in s_exclude:
            continue
        dist_from_mut = int(float(l[8]))
        for key, col in d_columns.items():
            for i_radii in range(len(l_radii)):
                radius = l_radii[i_radii]
                rmsd = l[col+i_radii]
                if rmsd == 'None':
                    continue
                else:
                    d_rmsds[radius][key][dist_from_mut] += [float(rmsd)]
    print 'looped over lines'


    ##
    ## histogram
    ##
    lines = []
    for dist_from_mut in d_rmsds[10]['CA'].keys():
        lines += ['%s %s\n' %(dist_from_mut, len(d_rmsds[10]['CA'][dist_from_mut]))]
    fd = open('histogram.gnuplotdata','w')
    fd.writelines(lines)
    fd.close()
##    plot_histogram('histogram',suffix,'distance (\305) from mutation','count of RMSDs',)
    plot_histogram('histogram',suffix,'distance ({\305}) from mutation','count of RMSDs',)


    ##
    ## histograms
    ##
    ## divide 1Angstrom into smaller units/tics
##    d_divs = {'CA':1000.,'heavy':500.,'phipsi':50.,'chi1':20.}
    d_divs = {'CA':500.,'heavy':500.,'phipsi':50.,'chi1':20.}
    for key in d_columns.keys():

        print 'temporarily dont plot histograms'
        break

        xlabel = '%s RMSD (%s)' %(d_octals[key],d_units[key])
        div = d_divs[key]
        l_plot_files = []
        for radius in l_radii:

            d_l_rmsds_histogram = {'all':[],}
            for dist in l_distances:
                d_l_rmsds_histogram[dist] = []

##            l_rmsds = []
            ## organize data
            for dist_from_mut in d_rmsds[radius][key].keys():
                for dist in l_distances:
                    if dist_from_mut <= dist:
                        d_l_rmsds_histogram[dist] += d_rmsds[radius][key][dist_from_mut]
                        break
##                l_rmsds += d_rmsds[radius][key][dist_from_mut]
                d_l_rmsds_histogram['all'] += d_rmsds[radius][key][dist_from_mut]

            ## count (yvalues)
            d_count = {'all':{},}
            l_range = []
            for dist in l_distances:
                d_count[dist] = {}
            for dist in d_count.keys():
                ## ???
                for rmsd in range(int(div*min(d_l_rmsds_histogram[dist])),int(div*max(d_l_rmsds_histogram[dist]))+1):
                    d_count[dist][rmsd] = 0
                for rmsd in d_l_rmsds_histogram[dist]:
                    d_count[dist][int(div*rmsd)] += 1
                for rmsd,count in d_count[dist].items():
                    if count > 10:
                        l_range += [rmsd]

            ## xtics
            d_xtics = {}
            for rmsd in range(
                int(min(l_range)/div),
                int(max(l_range)/div)+1,
                ):
                d_xtics['%s' %(rmsd)] = rmsd

            ## convert dict to txt
            for dist in d_count.keys():
                lines = []
                for rmsd in range(
                    int(min(l_range)),
                    int(max(l_range))+1,
                    ):
                    if not rmsd in d_count[dist].keys():
                        count = 0
                    else:
                        count = d_count[dist][rmsd]
                    lines += ['%s %s\n' %(rmsd/div, count)]
                fn = 'histogram_%s_distfrommut%s_radius%s.gnuplotdata' %(key, dist, radius)
                fd = open(fn,'w')
                fd.writelines(lines)
                fd.close()


##            plot_histogram(prefix,suffix,'RMSD (Angstrom)','count of RMSDs',d_xtics=d_xtics)
        for dist in d_count.keys():
            l_plot_files = []
            for radius in l_radii:
                fn = 'histogram_%s_distfrommut%s_radius%s.gnuplotdata' %(key, dist, radius)
                l_plot_files += [fn]
            prefix = 'histogram_%s_distfrommut%s' %(key,dist,)
            if dist == 'all':
                title = 'all distances from site of mutation'
            else:
                title = '%s {\305} from site of mutation' %(dist)
            plot_histogram(prefix,suffix,xlabel,'count of RMSDs',l_plot_files=l_plot_files,title=title)

    for key in d_columns.keys():
        l_plot_files = []
        

    ##        
    ##
    ##
    for key in d_columns.keys():
##        xlabel = 'distance {Symbol \305} from mutation'
        xlabel = 'distance ({\305}) from mutation'
        ylabel = 'RMSD ({Symbol %s}) within sphere' %(d_units[key],)
        title = '%s' %(d_octals[key])
        for i_radii in range(len(l_radii)):
            radius = l_radii[i_radii]
            lines = []
            for dist_from_mut in d_rmsds[radius][key].keys():
                l_rmsds = d_rmsds[radius][key][dist_from_mut]
##                if len(l_rmsds) == 0:
##                if len(l_rmsds) < 100:
##                if len(l_rmsds) < 400:
                if len(l_rmsds) < 1000:
                    if dist_from_mut < 10:
                        print 'skipping', dist_from_mut, key, radius
                    continue
                average,stderr = statistics.do_stderr(l_rmsds)
                lines += ['%s %s %s\n' %(dist_from_mut+0.1*i_radii,average,stderr)]
            fd = open('%s_%s.gnuplotdata' %(key,radius),'w')
            fd.writelines(lines)
            fd.close()
        plot_scatter(key,xlabel,ylabel,l_radii,title=title,)

    plot_scatter_combined(suffix)

    return


def plot_scatter_combined(suffix):

    prefix = 'scatter_combined_%s' %(suffix)
    xlabel = 'distance ({\305}) from mutation'
    ylabel = 'C_{/Symbol a} and heavy atom RMSD ({\305}) within 10{\305} sphere'
    y2label = '{/Symbol f}{/Symbol y} and {/Symbol c}_1 RMSD ({/Symbol \260}) within 10{\305} sphere'
    l_keys = ['CA','heavy','phipsi','chi1',]
    l_axes = ['x1y1','x1y1','x1y2','x1y2',]
    radius = 10
    d_octals = {'CA':'C_{/Symbol a}','heavy':'heavy atom','phipsi':'{/Symbol f}{/Symbol y}','chi1':'{/Symbol c}_1'} ## octals

    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(prefix),
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
##            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
        'set xlabel "%s"\n' %(xlabel),
        'set ylabel "%s"\n' %(ylabel),
    ]

    gnuplotsettings += ['set y2label "%s"\n' %(y2label)]
    gnuplotsettings += ['set yrange [0.2:0.9]\n']
    gnuplotsettings += ['set y2range [8:34]\n']
    gnuplotsettings += ['set ytics nomirror\n']
    gnuplotsettings += ['set y2tics\n']
    gnuplotsettings += ['set tics out\n']
##    gnuplotsettings += ['set autoscale y\n']
##    gnuplotsettings += ['set autoscale y2\n']

    line_plot = 'plot '
    line_plot += '[-1:][:]'
    for i_key in range(len(l_keys)):
        key = l_keys[i_key]
        fn = '%s_%s.gnuplotdata' %(key,radius)
        title = d_octals[key]
        line_plot += '"%s_%s.gnuplotdata" lt 1 ps 2 pt 7 lc %i lw 5 t "%s RMSD" w errorb axis %s, ' %(
            key,radius, i_key, title, l_axes[i_key]
            )
    line_plot = line_plot[:-2]
    line_plot += '\n'
    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(prefix))

    os.system('convert %s.ps %s.png' %(prefix,prefix,))

    stop_tmp
    for key in l_keys:
        os.remove('%s_%s.gnuplotdata' %(key,radius))
    os.remove('%s.gnuplotsettings' %(prefix))
    os.remove('%s.ps' %(prefix))
  
    return


def plot_histogram(prefix,suffix,xlabel,ylabel,d_xtics=None,l_plot_files=None,title=None,):

    gnuplotsettings = []
    gnuplotsettings += [
##        'set terminal postscript eps enhanced color "Helvetica" 48\n',
##        'set output "%s.ps"\n' %(prefix),
##        'set size 4,4\n', ## scale 400%
        'set terminal wxt 1 enhanced\n',
        ]
    gnuplotsettings += ['set style fill transparent solid 0.75 noborder\n']
    gnuplotsettings += [
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
##            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
        'set xlabel "%s"\n' %(xlabel),
        'set ylabel "%s"\n' %(ylabel),
    ]
##    gnuplotsettings += ['set style data histogram\n']
    if d_xtics:
        line_xtic = 'set xtics ('
        for xtic in d_xtics.keys():
            line_xtic += '"%s" %s, ' %(xtic, d_xtics[xtic])
        line_xtic = line_xtic[:-2]+')\n'
        gnuplotsettings += [
            line_xtic,
            'set xtics rotate\n',
        ]

    if title:
        gnuplotsettings += ['set title "%s"\n' %(title)]

    if l_plot_files:
##        gnuplotsettings += ['set terminal png\n']
##        gnuplotsettings += ['set output "%s.png" size "1600,1200"\n' %(prefix)]
##        gnuplotsettings += ['set terminal png truecolor nocrop enhanced font "Helvetica" 9 size 512,384\n']
##        gnuplotsettings += ['set output "%s.png"\n' %(prefix)]
        line = 'plot [0:]'
        l_titles = ['r=[0:10[','r=[10:20[','r=[20:40[',]
        for i_fn in range(len(l_plot_files)):
            fn = l_plot_files[i_fn]
            print fn
            line += '"%s" u 1:2 lt 1 lc %i t "%s" w boxes, ' %(
                fn,i_fn+1,l_titles[i_fn],
                )
        line = line[:-2]
        line += '\n'
        gnuplotsettings += [line]
    else:
        gnuplotsettings += [
            'plot [-1:]"%s.gnuplotdata" u 1:2 lt 1 lc 0 notitle w boxes\n' %(
    ##        'plot [-1:]"%s.gnuplotdata" u 1:2 lt 1 lc 0 notitle w boxes\n' %(
                prefix,
                )
            ]

    fd = open('%s.gnuplotsettings' %(prefix),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(prefix))

    if not l_plot_files:
        os.system('convert %s.ps %s_%s.png' %(prefix,prefix,suffix,))
        os.remove('%s.gnuplotdata' %(prefix,))
##        os.remove('%s.ps' %(prefix))
        os.remove('%s.gnuplotsettings' %(prefix))
##    else:
##        for fn in l_plot_files:
##            os.remove(fn)

    return


def plot_scatter(prefix,xlabel,ylabel,l_radii,title=None,):

    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(prefix),
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
##            'set title "%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}" ,4\n' %(plotname, title1, jobid, chains, cutoff_distance),
        'set xlabel "%s"\n' %(xlabel),
        'set ylabel "%s"\n' %(ylabel),
    ]

    if title:
        gnuplotsettings += ['set title "%s"\n' %(title)]

    l_titles = ['r=[0:10[','r=[10:20[','r=[20:40[',]
    line_plot = 'plot '
    line_plot += '[-1:][0:]'
    for i_radii in range(len(l_radii)):
        radius = l_radii[i_radii]
        line_plot += '"%s_%s.gnuplotdata" lt 1 ps 1 pt 7 lc %i t "%s" w errorb, ' %(prefix,radius,i_radii+1,l_titles[i_radii])
    line_plot = line_plot[:-2]
    line_plot += '\n'
    gnuplotsettings += [line_plot]

    fd = open('%s.gnuplotsettings' %(prefix),'w')
    fd.writelines(gnuplotsettings)
    fd.close()

    os.system('/usr/bin/gnuplot %s.gnuplotsettings' %(prefix))

    os.system('convert %s.ps %s.png' %(prefix,prefix,))

    for radius in l_radii:
        ## used for combined scatter plot...
        if radius == 10:
            continue
        os.remove('%s_%s.gnuplotdata' %(prefix,radius))
    os.remove('%s.gnuplotsettings' %(prefix))
##    os.remove('%s.ps' %(prefix))

    return


if __name__ == '__main__':
    main()

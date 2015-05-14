import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import collections
import math
import sys
##from colour import Color
import colorsys
import operator
import os
import argparse
import matplotlib.patheffects as path_effects
import getpass
import platform
import time

projection = 'merc'
resolution = 'h'

def main():

    args = parse_arguments()

    draw(args)

    return


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('--out')
    parser.add_argument('--input_file', '--in', required=True)
    ## left, lower, right, upper
##    ('Africa', -21, -36, 55, 38),
    parser.add_argument(
        '--llcrnrlon', '--western', '--W', '--left', type=float, default=-21)
    parser.add_argument(
        '--llcrnrlat', '--southern', '--S', '--lower', '--bottom',
        type=float, default=-36)
    parser.add_argument(
        '--urcrnrlon', '--eastern', '--E', '--right',
        type=float, default=55)
    parser.add_argument(
        '--urcrnrlat', '--northern', '--N', '--upper', '--top',
        type=float, default=38)

    ## columns/fields
    parser.add_argument('--colors', type=int, default=None)
    parser.add_argument('--labels', type=int)
    parser.add_argument('--text', type=int, default=0)
    parser.add_argument('--lat', '--latitude', type=int, default=1)
    parser.add_argument('--lon', '--longitude', type=int, default=2)
    parser.add_argument('--markersize', '--markersizes', type=int)
    parser.add_argument('--markerstyle', type=int)
    ## The text offset is DPI dependent.
    parser.add_argument('--offsetx', type=int)
    parser.add_argument('--offsety', type=int)

    parser.add_argument(
        '--mapbackground', choices=('fillcontinents', 'bluemarble'),
        default='fillcontinents')
    parser.add_argument(
        '--drawparallels', '--parallels', action='store_true')
    parser.add_argument(
        '--drawmeridians', '--meridians', action='store_true')

    parser.add_argument('--log', action='store_true')
    parser.add_argument('--loc', type=int, default=3, help='legend location')

    ## values
    parser.add_argument('--offset', type=float, default=0.15)
    parser.add_argument('--ncol', '--ncols', type=int, default=2)
    parser.add_argument('--fillcontinents_color', '--color_fillcontinents', default='0.85')
    parser.add_argument('--markerfacecolor', default=None)
    parser.add_argument('--text_size', default='x-small')
    parser.add_argument('--colorlist', nargs='+')
    parser.add_argument('--dpi', default=300, type=int)

    ## booleans
    parser.add_argument('--no_text_labels', action='store_true', default=False)
    parser.add_argument('--no_legend', '--no-legend', action='store_true', default=False)

    ## files
    parser.add_argument('--offsets')

    args = parser.parse_args()

    assert args.urcrnrlat > args.llcrnrlat
    assert args.urcrnrlon > args.llcrnrlon

    return args


def parse_args():

    return


def draw(args):
    
    # draw Basemap
    lat_0 = (args.llcrnrlat+args.urcrnrlat)/2
    lon_0 = (args.llcrnrlon+args.urcrnrlon)/2
    m = Basemap(
        projection=projection, resolution = 'l',
        lat_0=lat_0, lon_0=lon_0,
        llcrnrlon=args.llcrnrlon, llcrnrlat=args.llcrnrlat,
        urcrnrlon=args.urcrnrlon, urcrnrlat=args.urcrnrlat,
        )
    m.drawcoastlines()
    m.drawcountries()
    if args.mapbackground == 'bluemarble':
        m.bluemarble()
    elif args.mapbackground == 'fillcontinents':
        m.fillcontinents(color = args.fillcontinents_color)
    elif args.mapbackground == 'etopo':
        m.etopo()
    elif args.mapbackground == 'shadedrelief':
        m.shadedrelief()
    m.drawmapboundary()
    if args.drawmeridians:
        m.drawmeridians(np.arange(0, 360, 1), labels=[1, 1, 1, 1], fontsize='xx-small')
    if args.drawparallels:
        m.drawparallels(np.arange(-90, 90, 1), labels=[1, 1, 1, 1], fontsize='xx-small')

    l_colors_VGA16 = (
        (1,0,0), # red
        (0,1,0), # green/lime
        (0,0,1), # blue
        (0,1,1), # cyan/aqua
        (1,0,1), # magenta/fuchsia
        (1,1,0), # yellow
        (1,.5,0), # orange
        (1,1,1), # white
    #    '#800000', # maroon
        '#800080', # purple
    #    '#008000', # dark green
    #    '#808000', # olive
        '#000080', # navy
        '#008080', # teal
        (0,0,0), # black
        '#808080', # dark grey / silver?
        '#C0C0C0', # light grey / silver?
    #    (.5,1,0), # green yellow
        (0,1,.5), # green cyan
        (0,.5,1), # blue cyan
        (.5,0,1), # blue magenta
        (1,0,.5), # red magenta
        (1,.5,.5), ## red?
        (.5,1,.5), ## green?
        (.5,.5,1), ## blue?
        )
    l_colors_colorbrewer = [
            (166,206,227),
            (31,120,180),
            (178,223,138),
            (51,160,44),
            (251,154,153),
            (227,26,28),
            (253,191,111),
            (255,127,0),
            (202,178,214),
            (106,61,154),
            (255,255,153),
            (177,89,40),
            ]
    for i, color in enumerate(l_colors_colorbrewer):
        l_colors_colorbrewer[i] = (l_colors_colorbrewer[i][0]/255,l_colors_colorbrewer[i][1]/255,l_colors_colorbrewer[i][2]/255)

##    ## https://docs.python.org/3/library/colorsys.html
##    l_colors = [(1,1,1),'#808080','#C0C0C0']+[Color(hue=h/360., saturation=1, luminance=.5).rgb for h in range(0,360,int(360/10))]
    l_colors = [(1,1,1),'#808080','#C0C0C0']+l_colors_colorbrewer
    l_colors = l_colors+l_colors
#    l_colors = l_colors_VGA16
#    l_colors = ['#f1a340', '#998ec3']  # colorbrewer2.org, 3 data classes, diverging, scheme 4
#    l_colors = ['#e41a1c', '#984ea3']  # colorbrewer2.org, 4 data classes, qualitative, scheme 6
#    l_colors = ['#d01c8b', '#4dac26']  # colorbrewer2.org, 4 data classes, diverging, scheme 2
#    l_colors = ['#e41a1c', '#4daf4a']
#    l_colors = ['rgb(27,158,119)','rgb(217,95,2)','rgb(117,112,179)','rgb(231,41,138)','rgb(102,166,30)']  # 5 classes, qual, scheme 2
#    l_colors = ['rgb(27,158,119)','rgb(217,95,2)','rgb(117,112,179)','rgb(231,41,138)','rgb(102,166,30)','rgb(230,171,2)']
#    l_colors = ['rgb(27,158,119)','rgb(217,95,2)','rgb(117,112,179)','rgb(231,41,138)','rgb(102,166,30)','rgb(230,171,2)','rgb(166,118,29)']
#    l_colors = ['rgb(228,26,28)','rgb(55,126,184)','rgb(77,175,74)','rgb(152,78,163)','rgb(255,127,0)','rgb(255,255,51)','rgb(166,86,40)']
#    l_colors_blue = ['rgb(222,235,247)','rgb(158,202,225)','rgb(49,130,189)']  # 3, sequential, single hue 1 (blue)
#    l_colors_green = ['rgb(229,245,224)','rgb(161,217,155)','rgb(49,163,84)']  # 3, seq, single hue 2 (green)
#    l_colors_orange = ['rgb(254,230,206)','rgb(253,174,107)','rgb(230,85,13)']  # 3, seq, single hue 4 (orange)
#    l_colors_red = ['rgb(254,224,210)','rgb(252,146,114)','rgb(222,45,38)']  # 3, seq, single hue 6 (red)
#    l_colors = [l_colors_red[0], l_colors_red[2], l_colors_green[0], l_colors_green[2], l_colors_blue[0], l_colors_blue[2],]
#    l_colors = [l_colors_red[0], l_colors_red[2], l_colors_green[0], l_colors_green[2], l_colors_blue[2],]
#    l_colors = [tuple(int(i)/255 for i in color[4:-1].split(',')) for color in l_colors]  # tmp!!! for colorbrewer2.org export
#    ## 4, qual, colorblind safe
#    l_colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c']
    ## 9, qual, printer friendly (4,5,6,7,8 are subsets)
    l_colors = [
        '#e41a1c', '#377eb8', '#4daf4a',
        '#984ea3', '#ff7f00', '#ffff33',
        '#a65628', '#f781bf', '#999999',
        ]

    l_markerstyles = ['o','s','^','v','*','D','p']

#    for i, color in enumerate(l_colors_colorbrewer):
#        l_colors_colorbrewer[i] = (l_colors_colorbrewer[i][0]/255,l_colors_colorbrewer[i][1]/255,l_colors_colorbrewer[i][2]/255)

    if args.colorlist:
        l_colors = ['#{}'.format(color) for color in args.colorlist]

    with open(args.input_file) as f:
        lines = f.readlines()

    d_colors = {}
    d_markerstyles = {}

    set_labels = set()
    for line in lines:

        ## skip blank lines
        if not line.strip():
            continue
        ## skip comment lines
        if line[0] == '#':
            continue
        ## split line into list of column fields
        l = line.rstrip().split('\t')

        ## parse values from line
        lat = float(l[args.lat])
        lon = float(l[args.lon])
        label = l[args.labels]
        text = l[args.text]

        ## tmp!!!
        try:
                label = {
#                    '1': 'REC approval, SEQUENCED',
#                    '2': 'REC approval, SEQUENCING',
                    '1': 'Sequenced/Sequencing',
                    '2': 'Sequenced/Sequencing',
#                    '3': 'REC approval, OTHERS',
#                    '5': 'EXISTING STUDIES NOT YET APPROVED',
                    '3': 'Other collections, Approved',
                    '5': 'Other collections, Not yet approved',
#                    '4': 'NEW COLLECTIONS',
                    '4': 'Planned collections',
                    'MENTOR': 'MENTOR Initiative',
#                    '6': '1000G phase 3',
                    '6': 'Sequenced/Sequencing',
                    }[label]
        except:
            pass

        ## markersize
        if args.markersize:
            if l[args.markersize] != 'NA':
                markersize = float(l[args.markersize])
            else:
                markersize = 10
        else:
            markersize = 10
        if args.log:
            markersize = cnt2markersize(markersize)

        ## markerstyle
        if args.markerstyle:
            try:
                markerstyle = d_markerstyles[l[args.markerstyle]]
            except KeyError:
                markerstyle = l_markerstyles[len(list(d_markerstyles.keys()))]
                d_markerstyles[l[args.markerstyle]] = markerstyle
        else:
            markerstyle = 'o'

        ## color
        if args.colors:
            color = l[args.colors]
        elif args.markerfacecolor:
            color = args.markerfacecolor
        else:
            try:
                color = d_colors[label]
            except KeyError:
                color = l_colors[len(list(d_colors.keys()))]
                d_colors[label] = color

        ## convert lon and lat to x and y
        x, y = m(lon, lat)
        if args.offsetx and args.offsety:
            offset_extra = (float(l[args.offsetx]), float(l[args.offsety]))
        else:
            offset_extra = (0, 0)
        ## get text position
        x_text, y_text = m(lon+offset_extra[0], lat+args.offset+offset_extra[1])
        if args.no_text_labels == False:
            text = plt.text(
                x_text, y_text, text.replace('\\n','\n'),
                horizontalalignment='center', size=args.text_size,
                color=color)
            text.set_path_effects([
                path_effects.Stroke(linewidth=2, foreground='black'),
                path_effects.Normal()])
        ## plot data point
        m.plot(
            x, y, markerstyle, markersize=markersize, markerfacecolor=color, label="")
        ## plot outside map to have label appear (not very elegant solution...)
        x, y = m(-60, -60)
        if not label in set_labels:
            set_labels.add(label)
            m.plot(x, y, 'o', markersize=10, markerfacecolor=color, label=label)

    if args.markersize:
#        for i in range(1,3):
        for i in range(0,3): # make this log dependent and input range dependent...
            m.plot(
                x, y, 'o', markersize=cnt2markersize(pow(10, i)),
                markerfacecolor='white', label=pow(10, i))

    if args.markerstyle:
        for label, markerstyle in sorted(d_markerstyles.items()):
            m.plot(
                x, y, markerstyle, markersize=cnt2markersize(10),
                markerfacecolor='white', label='Sequencing depth of coverage '+label)

    distance_from_plot = 0.05
    print('ncol', args.ncol)
    plt.legend(
        ncol=args.ncol, shadow=True,
##        prop={'size':8.5},
        prop={'size':6.5},
        numpoints=1,
        labelspacing=1.0,
        loc=args.loc,
        fancybox=True,
##        handleheight=2.5,
        borderpad=0.7,
        framealpha=0.5,
        )

    if args.out:
        out = '{}.jpg'.format(args.out)
    else:
        out = '{}_{}_{}.jpg'.format(
            os.path.splitext(os.path.basename(args.input_file))[0], projection, args.mapbackground)
    plt.savefig(out, dpi=args.dpi)

    plt.title('{} {} {}'.format(
            getpass.getuser(),
            '.'.join(platform.node().split('.')[1:]),
            time.strftime('%a %Y%b%d', time.gmtime())))

    plt.close()

    return


def cnt2markersize(cnt):

    log = 5  # smaller number bigger differences
    smallest_size = 3
    assert type(cnt) in (int, float)
    markersize = smallest_size*math.log(log*cnt, log)
#    markersize = math.sqrt(cnt)

    return markersize


if __name__ == '__main__':
    main()

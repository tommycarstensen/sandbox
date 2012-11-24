#!/bin/env /software/bin/python
# $Id: CPMG.py 211 2007-03-28 17:00:19Z tc $
## modify spectra.list and finish pk.tcl section
## Copyright (C) Tommy Carstensen, University College Dublin, 2006-2007

def main():

    import math

    ##
    ## variable definition
    ##

##    center_x = 4.630; center_y = 21.000 ## find them or derivatives (width and start point) in procpar
##    p0 = '-66.00'; p1 = '0.00'

    OS = 'windows'
    file_prefix = 'ACBP'; sparky_peaklist = 'ACBP_H1C13_methyl.list' ##None
    igor_prefix = '500'
    tiles = [3,3]; borders = [10,10]
    
    ##
    ## parsing
    ##

    procpar = procpar_parse()

    if OS == 'linux':
        overflows = log_parse()
    elif OS == 'windows':
        overflows = {}

    ##
    ## processing
    ##

    if OS == 'linux':

        sort_pseudo3D(len(procpar['ncyc_cp']), procpar['ni'], procpar['np'])

        var2pipe(len(procpar['ncyc_cp']), procpar['sw'], procpar['sw1'], procpar['sfrq'], procpar['dfrq'], center_x, center_y)

        nmrPipe(len(procpar['ncyc_cp']),p0,p1,overflows)

        sparkypeaklist2seriestabpeaklist(sparky_peaklist, file_prefix)

        seriestab_create(file_prefix)

    lines_seriestab = seriestab_parse(file_prefix, procpar['ncyc_cp'], procpar['time_T2'], igor_prefix)

    sslinear = fit_linear(lines_seriestab, file_prefix, sparky_peaklist, OS)

    ssnonlinear = fit(lines_seriestab, file_prefix, sparky_peaklist, OS)

    ftest(sslinear,ssnonlinear, float(len(lines_seriestab)-2), float(len(lines_seriestab)-5))

    montage(file_prefix, tiles, borders, len(lines_seriestab))
    
    return


def seriestab_parse(file_prefix, ncycs, T2, igor_prefix):

    '''
    parse intensity from file_prefix.seriestabout,
    calculate R2 using I and I0,
    and write R2 to file_prefix.gnudata
    '''

    import math

    ## read seriestab output
    fd = open('%s.seriestabout' %(file_prefix), 'r')
    lines_seriestab = fd.readlines()[13:]
    fd.close()

    ## determine I0s from seriestab output
    I0s = []
    for line in lines_seriestab:
        I0 = []
        for i in range(len(ncycs)):
            if ncycs[i] == 0:
                I0.append(float(line[69+8*i:69+8*(i+1)]))
        I0s.append(sum(I0)/len(I0))

    ## write gnuplot data file from seriestab output (use assignment from peaks_Sparky)
    ## write assignment row to gnuplot data file
    lines_gnuplot = ['%s_#freqCPMG ' %(igor_prefix)]
    for i in range(len(lines_seriestab)):
        line = lines_seriestab[i].split()
        ass = line[6]
        lines_gnuplot[0] += '%s_%s/%s ' %(igor_prefix, i, ass)
    lines_gnuplot[0] += '\n'
    ## write freqCPMG column and intensity columns to gnuplot data file
    for i in range(len(ncycs)):
        if ncycs[i] != 0:
            ## write ncyc column
            lines_gnuplot.append('%4i' %(ncycs[i]/T2))
            ## write intensity columns
            for j in range(len(lines_seriestab)):
	        line = lines_seriestab[j].split()
                I = float(line[7+i])
                R2 = -(1/T2)*math.log( I/I0s[j] )
                lines_gnuplot[-1] += '%8.3f' %(R2)
            lines_gnuplot[-1] += ('\n')
    ## write gnuplot data file
    fd = open('%s.gnudata' %(file_prefix), 'w')
    fd.writelines(lines_gnuplot)
    fd.close()
    
    return lines_seriestab


def ftest(ss1, ss2, df1, df2):

    ##
    ## linear or nonlinear model fits better?
    ##

    print len(ss1), len(ss2)
    for i in range(len(ss1)):
        F = ((ss1[i]-ss2[i])/ss2[i])/((df1-df2)/float(df2))
        P = fdist(F,df1-df2,df1)
        print i, (ss1[i]-ss2[i])/float(ss1[i]-ss2[i]), F, P

    return


def fdist(F,df_n,df_d):

## http://mathworld.wolfram.com/F-Distribution.html (eq. 4)
## http://mathworld.wolfram.com/RegularizedBetaFunction.html
## http://mathworld.wolfram.com/IncompleteBetaFunction.html (eq. 1)
## http://mathworld.wolfram.com/BetaFunction.html (eq. 18)

    import math
    z = df_n*F/float(df_d+df_n*F)
    a = df_n/2.
    b = df_d/2.
    lim1 = 0.
    d_lim2 = {'numerator':z,'denominator':1.}
    inv_stepsize = 1000000.
    function_beta_regularized = {}

    for key in d_lim2:

        lim2 = d_lim2[key]

        deltax = (lim2-lim1)/inv_stepsize
        area = 0.

        x = lim1+deltax
        while x <= lim2:
            y = ((x)**(a-1)) * ((1-x)**(b-1))
            area += y*deltax
            x += deltax

        function_beta_regularized[key] = area

    P = 1-function_beta_regularized['numerator']/function_beta_regularized['denominator']
    
    return P


def procpar_parse():
   
    procpar = {}

    fd = open('procpar' ,'r')
    lines = fd.readlines()
    fd.close()
    parameters = ['ni', 'np', 'time_T2', 'sfrq', 'dfrq', 'ncyc_cp', 'sw', 'sw1']
    for i in range(len(lines)):
        line = lines[i]
	for parameter in parameters:
	    if '%s ' %(parameter) == line[:len(parameter)+1]:
	        count = int(lines[i+1].split()[0])
		if count > 1:
	            procpar[parameter] = lines[i+1].split()[1:count+1]
                    for j in range(len(procpar[parameter])):
                        procpar[parameter][j] = int(procpar['ncyc_cp'][j])
		elif count == 1:
		    par = lines[i+1].split()[1]
		    if '.' in par:
	                procpar[parameter] = float(par)
		    else:
		        procpar[parameter] = int(par)

    return procpar

    
def fit(lines_seriestab, file_prefix, sparky_peaklist, OS):

    import os, math, cmath

    ##
    ## (tau,R2eff)-fit and kex,R2-determ.
    ##

    lines_fit_o = ['#atom Ra Rb ka kb dw\n']

    fd = open('%s.gnudata' %(file_prefix),'r')
    gnudata = fd.readlines()
    fd.close()

    ## set a list of sums of squares
    SS = []

    for i in range(len(lines_seriestab)):

        if not sparky_peaklist:
            ass = assignments_Sparky[peakindexes_NMRPipe[i]]
        else:
            ass = lines_seriestab[i][58:68].strip()

        print 'nonlinear fit for %s' %(ass)

        lines_gnuplot_script = ([
	    '#freqCP = x\n'
	    'tauCP(x) = 1/(4.*x)\n',
##            'psi(x) = (Ra-Rb-ka+kb)**2-deltaomega**2+4.*ka*kb\n', ## tc
            'psi(x) = kex**2-deltaomega**2\n', ## kt
##            'epsilon(x) = 2*deltaomega*(Ra-Rb-ka+kb)\n', ## tc
            'epsilon(x) = -2*deltaomega*kex*(2*pA-1)\n', ## kt
##            'Dpos(x) = (1./2.)*( 1+(psi(x)+2*deltaomega**2)/sqrt(psi(x)**2+epsilon(x)**2))\n',
##            'Dneg(x) = (1./2.)*(-1+(psi(x)+2*deltaomega**2)/sqrt(psi(x)**2+epsilon(x)**2))\n',
            'Dpos(x) = (1./2.)*( 1+(psi(x)+2*deltaomega**2)/sqrt(psi(x)**2+epsilon(x)**2))\n',
            'Dneg(x) = (1./2.)*(-1+(psi(x)+2*deltaomega**2)/sqrt(psi(x)**2+epsilon(x)**2))\n',
            'etapos(x) = (tauCP(x)/sqrt(2))*sqrt( psi(x)+sqrt(psi(x)**2+epsilon(x)**2))\n',
            'etaneg(x) = (tauCP(x)/sqrt(2))*sqrt(-psi(x)+sqrt(psi(x)**2+epsilon(x)**2))\n',
##            'f(x) = (1./2.)*(   Ra+Rb+ka+kb-(1/tauCP(x))*(  acosh( Dpos(x)*cosh(etapos(x)) - Dneg(x)*cos(etaneg(x)) )  )   )\n',
            'f(x) = R+(1./2.)*(   kex-(1/tauCP(x))*(  acosh( Dpos(x)*cosh(etapos(x)) - Dneg(x)*cos(etaneg(x)) )  )   )\n',
##            'Ra = 17\n', ## tc
##            'Rb = 23\n', ## tc
##            'ka = 29\n', ## tc
##            'kb = 19\n', ## tc
            'kex = 2000\n', ## kt
            'pA = 0.9\n', ## kt
            'R = 7\n', ## kt
            'deltaomega = 1000\n',
##            'fit f(x) "%s.gnudata" using 1:%s via Ra,Rb,ka,kb,deltaomega\n' %(file_prefix, i+2),
            'fit f(x) "%s.gnudata" using 1:%s via kex,pA,R,deltaomega\n' %(file_prefix, i+2),
            'set terminal png\n',
            'set output "%s%s.png"\n' %(file_prefix, i),
            'plot [][0:]"%s.gnudata" using 1:%s title "%s", f(x)\n' %(file_prefix, i+2, ass)
            ])

        ## write gnuplot script file
        fd = open('%s.gnuscript%s' %(file_prefix, i), 'w')
        fd.writelines(lines_gnuplot_script)
        fd.close()
        ## execute gnuplot script file and write gnuplot fit log
        if OS == 'linux':
            os.system('/usr/bin/gnuplot %s.gnuscript%s' %(file_prefix, i))
        elif OS == 'windows':
            os.system('gnuplot %s.gnuscript%s' %(file_prefix, i))
	## remove gnuplot script file
	os.remove('%s.gnuscript%s' %(file_prefix, i))

	## parse gnuplot fit log
	if OS == 'linux':
            lines_fit_i = os.popen('tail -n18 fit.log').readlines()
        elif OS == 'windows':
            fd = open('fit.log','r')
            lines_fit_i = fd.readlines()[-18:]
            fd.close()
##	Ra = float(lines_fit_i[3].split()[2])
##	Rb = float(lines_fit_i[4].split()[2])
##	ka = float(lines_fit_i[5].split()[2])
##	kb = float(lines_fit_i[6].split()[2])
##	dw = float(lines_fit_i[7].split()[2])
	kex = float(lines_fit_i[5].split()[2])
	pA = float(lines_fit_i[6].split()[2])
	R = float(lines_fit_i[7].split()[2])
	dw = float(lines_fit_i[8].split()[2])
        
	## write parsed fit data to other fit log
##        lines_fit_o += ['%10s %6.3f %6.3f %6.3f %6.3f %8.3f\n' %(ass, Ra, Rb, ka, kb, dw)]
        lines_fit_o += ['%10s %6.3f %6.3f %6.3f %8.3f\n' %(ass, kex, pA, R, dw)]

        ##
        ## SS calculation for use in F-test comparing nonlinear and linear fits
        ##

        n = len(gnudata[1:])
        deltaomega = dw
        sumdiffsq = 0
        sumdiff = 0
        for line in gnudata[1:]:
            x = float(line.split()[0])
            y = float(line.split()[i+1])
	    tauCP = 1/(4.*x)
##            psi = (Ra-Rb-ka+kb)**2-deltaomega**2+4.*ka*kb
	    psi = kex**2-deltaomega**2
##            epsilon = 2*deltaomega*(Ra-Rb-ka+kb)
	    epsilon = -2*deltaomega*kex*(2*pA-1)
            Dpos = (1./2.)*( 1+(psi+2*deltaomega**2)/math.sqrt(psi**2+epsilon**2))
            Dneg = (1./2.)*(-1+(psi+2*deltaomega**2)/math.sqrt(psi**2+epsilon**2))
            etapos = (tauCP/math.sqrt(2))*math.sqrt( psi+math.sqrt(psi**2+epsilon**2))
            etaneg = (tauCP/math.sqrt(2))*math.sqrt(-psi+math.sqrt(psi**2+epsilon**2))
##            yhat = (1./2.)*(   Ra+Rb+ka+kb-(1/tauCP)*(  cmath.acosh( Dpos*cmath.cosh(etapos) - Dneg*math.cos(etaneg) )  )   )
            yhat = R+(1./2.)*(   kex-(1/tauCP)*(  cmath.acosh( Dpos*cmath.cosh(etapos) - Dneg*math.cos(etaneg) )  )   )
            diff = yhat-y
            sumdiff += diff ## abs??
            sumdiffsq += diff**2
        ss = abs(sumdiffsq-(sumdiff**2)/n)
        stop
        ## append ss calculated for 1 chemical group to list of sums of squares
        SS += [ss]

    os.remove('fit.log')

    fd = open('%s.fit' %(file_prefix), 'w')
    fd.writelines(lines_fit_o)
    fd.close()

    return SS


def fit_linear(lines_seriestab, file_prefix, sparky_peaklist, OS):

    ## standard linear regression (without replication) of Zar, Biostatistical Analysis, ch. 17

    import os, math

    fd = open('%s.gnudata' %(file_prefix),'r')
    gnudata = fd.readlines()
    fd.close()

    ## set a list of sums of squares
    SS = []
    
    for i in range(len(lines_seriestab)):

        if not sparky_peaklist:
            ass = assignments_Sparky[peakindexes_NMRPipe[i]]
        else:
            ass = lines_seriestab[i][58:68].strip()

        ## reset sums
        n = len(gnudata[1:])
        sumxiyi = 0
        sumxi = 0
        sumyi = 0
        sumxisq = 0
        ## loop over x and y values and calculate sums
        for line in gnudata[1:]:
            x = float(line.split()[0])
            y = float(line.split()[i+1])
            sumxiyi += x*y
            sumxi += x
            sumyi += y
            sumxisq += x**2
        ## calculate regression statistics from sums
        sumxsq = sumxisq-(sumxi**2)/n
        sumxy = sumxiyi - (sumxi*sumyi)/n
        b = sumxy/sumxsq
        a = sumyi/n-b*sumxi/n

        sumdiffsq = 0
        sumdiff = 0
        for line in gnudata[1:]:
            x = float(line.split()[0])
            y = float(line.split()[i+1])
            yhat = a+b*x
            diff = yhat-y
            sumdiff += diff ## abs??
            sumdiffsq += diff**2
        ss = sumdiffsq-(sumdiff**2)/n

        ## append ss calculated for 1 chemical group to list of sums of squares
        SS += [ss]

        print ass
    stop

    return SS


def montage(file_prefix, tiles, borders, looplen):

    import os, math

    lines_montage = 'montage '
    ## define input
    for i in range(looplen):
        lines_montage += '%s%s.png ' %(file_prefix, i)
    ## set up tiles
    lines_montage += '-tile %sx%s ' %(int(tiles[0]), int(tiles[1]))
    ## add border
    geometry = [borders[0],borders[1]]
    lines_montage += '-geometry +%s+%s ' %(int(geometry[0]),int(geometry[1]))
    ## define output
    lines_montage += '%s.png' %(file_prefix)

    ## execute ImageMagick montage script
    os.system(lines_montage)
    
    for i in range(looplen):
        os.remove('%s%s.png' %(file_prefix, i))

    return
    

def seriestab_create(file_prefix):

    import os

    ##
    ## intensity calculation (seriesTab -in peaklist -out peaklist&intensities -list spectra)
    ##

    ##
    ## spectra list output
    ##
    lines = []
    for i in range(len(ncycs)):
        lines.append('%s/test.ft2\n' %(i))
    fd = open('spectra.list', 'w')
    fd.writelines(lines)
    fd.close()

    ## write seriestab output
    os.system('seriesTab -in %s.tab -out %s.seriestabout -list spectra.list -sum -dx 1 -dy 1' %(file_prefix, file_prefix))

    return


def sort_pseudo3D(looplen,ni,np):

    import os

    ##
    ## sort_pseudo3D
    ##
    if not os.path.isfile('fid'):
        stop
    os.system('~kte/bin/sort_pseudo3D -in fid -plane %s -mode 1 -ni %s -np %s' %(looplen, int(ni), int(np)))
    for i in range(looplen):
        if not os.path.isdir('%s' %(i)):
	    os.mkdir('%s' %(i))
        os.system('mv %s.fid %s/%s.fid' %(i, i, i))
	os.system('cp procpar %s/.' %(i))

    return


def nmrPipe(looplen,p0,p1,overflows={}):

    import os

    ##
    ## phasing and FT (nmrDraw; in fid, out ft2)
    ##
    for i in range(looplen):

        os.chdir('%s' %(i))

        if i+1 in overflows.keys():
            x1 = i+1
            xn = overflows[x1]/2
            line_LP = '| nmrPipe  -fn LP -after -x1 %i -xn %i -pred 1 -ord 16 \\\n' %(x1, xn)
        else:
            line_LP = ''

        lines = [
	    '#!/bin/csh\n\n',
	    'nmrPipe -in test.fid \\\n',
	    '| nmrPipe  -fn SOL                                    \\\n',
	    '| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \\\n',
	    '| nmrPipe  -fn ZF -auto                               \\\n',
	    '| nmrPipe  -fn FT -auto                               \\\n',
	    '| nmrPipe  -fn PS -p0 %s -p1 %s -di \\\n' %(p0, p1),
	    '| nmrPipe  -fn TP                                     \\\n',
            line_LP,
	    '| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \\\n',
	    '| nmrPipe  -fn ZF -auto                               \\\n',
	    '| nmrPipe  -fn FT -auto                               \\\n',
	    '| nmrPipe  -fn PS -p0 0.00 -p1 0.00 -di               \\\n',
	    '   -ov -out test.ft2\n',
	    ]

        fd = open('nmrproc.com', 'w')
        fd.writelines(lines)
        fd.close()
       
        os.system('chmod +x nmrproc.com')
#        os.system('nmrTerm -e csh nmrproc.com')
        os.system('nmrproc.com')

        os.chdir('../')
	
    return

def sparkypeaklist2seriestabpeaklist(sparky_peaklist, file_prefix):

    import os

    ##
    ## calculation between points and ppm (showhdr) only necessary if subsequent alignment with Sparky peaks
    ##

    header = os.popen('showhdr 0/test.ft2').readlines()
    xppmpp = (float(header[8][12:24])/float(header[9][12:24]))/float(header[6][12:24])
    yppmpp = (float(header[8][25:37])/float(header[9][25:37]))/float(header[6][25:37])
    xppm_orig = (float(header[10][12:24])+float(header[8][12:24]))/float(header[9][12:24])
    yppm_orig = (float(header[10][25:37])+float(header[8][25:37]))/float(header[9][25:37])

    ##
    ## parse Sparky peak list and convert to seriesTab format
    ##
    if sparky_peaklist:
        fd = open(sparky_peaklist, 'r')
        lines_Sparky = fd.readlines()[2:]
        fd.close()
	lines_tab = [
	    'VARS   INDEX X_AXIS Y_AXIS X_PPM Y_PPM VOL ASS\n',
            'FORMAT %5d %9.3f %9.3f %8.3f %8.3f %+e %10s\n\n'
	    ]
        for i in range(len(lines_Sparky)):
	    ass = lines_Sparky[i][0:17].strip()
	    xppm = float(lines_Sparky[i][18:28])
	    yppm = float(lines_Sparky[i][29:39])
            x = (xppm_orig-xppm)/xppmpp
            y = (yppm_orig-yppm)/yppmpp
	    lines_tab.append('%5i %9.3f %9.3f %8.3f %8.3f %+e %10s\n' %(i, x, y, xppm, yppm, 0, ass))
	fd = open('%s.tab' %(file_prefix), 'w')
	fd.writelines(lines_tab)
	fd.close()
	
    return


def log_parse():

    fd = open('log', 'r')
    lines = fd.readlines()
    fd.close()
    overflow_FIDs = []
    for line in lines:
        if line[26:30] == 'FID ' and line[35:48] == ' ADC overflow':
            overflow_FIDs += [int(line[30:34])]

    overflows = {}
    for FID in overflow_FIDs:
        ncyc = int(math.fmod(FID,30.))
        phase = int((FID-ncyc)/30)
        overflows[ncyc] = phase

    return overflows


def var2pipe(looplen, width_x, width_y, freq_x, freq_y, center_x, center_y):

    import os

    ##
    ## var2pipe (varian; in fid, out fid)
    ##
    for i in range(looplen):

        os.chdir('%s' %(i))
    
        fd = open('fid.com', 'w')
        fd.writelines([
            '#!/bin/csh\n\n',
            'var2pipe -in %s.fid -noaswap  \\\n' %(i),
            '  -xN              1024  -yN               256  \\\n',
            '  -xT               512  -yT               128  \\\n',
            '  -xMODE        Complex  -yMODE         States  \\\n',
            '  -xSW         %8.3f  -ySW         %8.3f  \\\n' %(width_x, width_y),
            '  -xOBS         %7.3f  -yOBS         %7.3f  \\\n' %(freq_x, freq_y),
            '  -xCAR           %5.3f  -yCAR          %6.3f  \\\n' %(center_x, center_y),
            '  -xLAB              H1  -yLAB             C13  \\\n',
            '  -ndim               2  -aq2D          States  \\\n',
            '  -out test.fid -ov\n',
            ])
        fd.close()

        os.system('chmod +x fid.com')
        os.system('fid.com')

        os.chdir('../')
	
    return

if __name__ == '__main__':
    main()

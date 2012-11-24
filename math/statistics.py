#!/usr/bin/python

# Copyright, Tommy Carstensen, University of Copenhagen, 2006, University College Dublin 2007-2011


import math, numpy


class tmp():


    def do_confidence_bands(self,l1,l2,t_crit,):

        ## tcrit = 1.960 if 5% two tail and n is infinite??? or t from do_regression???

        n = len(l1)
        print 'n', n

        sum_xy = 0
        sum_x = 0
        sum_y = 0
        sum_xx = 0
        sum_yy = 0
        for i in range(n):
            sum_xy += l1[i]*l2[i]
            sum_x += l1[i]
            sum_y += l2[i]
            sum_xx += l1[i]**2
            sum_yy += l2[i]**2
        SS_x = sum_xx-(sum_x**2)/float(n)
        SS_y = sum_yy-(sum_y**2)/float(n)
        ## sum of cross products (17.3)
        SS_xy = sum_of_cross_products = sum_xy-(sum_x*sum_y)/float(n)
        mean_x = sum_x/float(n)
        mean_y = sum_y/float(n)

        ## regression coefficient (17.4)
        b = SS_xy/SS_x
        print 'b', b
        ## intercept (17.7)
        a = mean_y-b*mean_x
        print 'a', a

        ## 17.10
        SS_total = SS_y
    ##    print 'SS_total', SS_total
        ## 17.11
        SS_regression = b*SS_xy
    ##    print 'SS_regression', SS_regression
        ## 17.13
        s_yx = SS_residual = SS_total-SS_regression
    ##    print 'SS_residual', SS_residual
        ## 17.14
        MS_residual = SS_residual/float(n-2)
    ##    print 'MS_residual', MS_residual

        #### variance of b (17.20)
        ##sb = math.sqrt(MS_residual/SS_x)
        #### (17.21)
        ##t = (b-0)/sb

        line = '%s+%s*x' %(a,b,)
        error = '%s*sqrt(%s*(1/%i.+((x-%s)**2)/%s))' %(t_crit,MS_residual,n,mean_x,SS_x,)
        f = 'f(x) = %s' %(line)
        g = 'g(x) = %s+%s' %(line,error)
        h = 'h(x) = %s-%s' %(line,error)
        print g
        print h

    ##    x = 10
    ##    stderr = math.sqrt(MS_residual*(1/float(n)+((x-mean_x)**2)/SS_x))
    ##    print stderr
    ##    print a+b*x+t_crit*stderr

        return f,g,h


    def Wilcoxon(self,l1,l2,):

        data = [l1,l2,]

        '''Paired sample test. Nonparametric analogue to the paired-sample t test.'''
        d_ranks = {}
        d_vals = {}
        val = None
        ## calculate differences
        data_diff = [abs(data[0][i]-data[1][i]) for i in range(len(data[0]))]
        ## sort differences
        data_diff.sort()
        ## rank data    
        for i in range(len(data_diff)):
            val = data_diff[i]
            if not d_ranks.has_key(val):
                rank = ((i+1)+(i+data_diff.count(val)))/2.
                d_ranks[val] = rank
                d_vals[rank] = val
        ## sum signed ranks
        Tpos = []
        Tneg = []
        data_diff = [data[0][i]-data[1][i] for i in range(len(data[0]))]
        for i in range(len(data_diff)):
            val = data_diff[i]
            rank = d_ranks[abs(val)]
            if val > 0:
                Tpos += [rank]
            else:
                Tneg += [rank]
        Tpos = sum(Tpos)
        Tneg = sum(Tneg)
        print 'Tpos', Tpos, 'Tneg', Tneg
        print min(Tpos, Tneg)

        return Tpos, Tneg


class continuous():

    def do_stddev(self,l):

        n = float(len(l))
        sumx = 0
        sumxx = 0
        for x in l:
            sumx += x
            sumxx += x*x
        average = sumx/n
        SS = sumxx-(sumx**2)/n
        var = SS/(n-1) ## division with n if population, n-1 if sample
        stddev = math.sqrt(var)
    ##    stderr = math.sqrt(var/n)

        return average, stddev


    def do_stderr(self,l):

        n = float(len(l))
        sumx = 0
        sumxx = 0
        for x in l:
            sumx += x
            sumxx += x*x
        average = sumx/n
        SS = sumxx-(sumx**2)/n
        var = SS/(n-1) ## division with n if population, n-1 if sample
    ##    stddev = math.sqrt(var)
        stderr = math.sqrt(var/n)

        return average, stderr


class misc():


    def normaldist(self, pop_stddev, pop_mean, Zlimit):

        ## old slow method calculating area
        
        import math
        area = .5
        Zmax = 4.
        inv_stepsize = 1000000.
        width = abs(pop_mean/float(inv_stepsize))
        txt = ''
        for x in range(int(inv_stepsize)):
            sample_mean = pop_mean+width*float(x)
            Z = (sample_mean-pop_mean)/pop_stddev
            numerator = pow(math.e,-(sample_mean-pop_mean)**2/(2*pop_stddev**2))
            denominator = pop_stddev*math.sqrt(2*math.pi)
            height = numerator/denominator
            area -= height*width
            if math.fmod(x,1000) == 0:
                txt += '%f %f\n' %(pop_mean+width*float(x), height)
                txt += '%f %f\n' %(pop_mean-width*float(x), height)
            if Z > Zlimit:
                print 'Z > Zlimit (%3.1f > %3.1f)' %(Z,Zlimit)
                break
            if Z > Zmax:
                print 'Z > Zmax (%3.1f > %3.1f)' %(Z,Zmax)
                break
        fd = open('gnuplot.dat','w')
        fd.write(txt)
        fd.close()

        return area


    def calc_sum_of_squares(self, l):

        n = len(l)
        sumx = 0.
        sumxx = 0.
        for x in l:
            x = float(x)
            sumx += x
            sumxx += x**2
    ##    SS = sumxx-(sumx**2)/n
    ##    mean = sumx/n

        return sumx, sumxx


    def faculty(self, n):

        n = int(n)
        fac = 1
        for i in range(1,n+1):
            fac *= i

        return fac


    def samplesize_to_estimate_mean_difference(self, variance_pooled,confidence_level,confidence_interval_halfwidth):

        import math
        alpha = float(1-confidence_level)
        d = confidence_interval_halfwidth
        n_min = n_min_0 = 1; n_max = n_max_0 = 1000
        print 'Make an initial of the sample size, n'
        n_guess = int(raw_input())
        while n_min != n_max:
            print 'Look up t(alpha(2)=%s,DF=%s) in Table B.3 on page App19' %(alpha,2*(n_guess-1))
            t = float(raw_input())
            n_calc = 2*variance_pooled*(t**2)/(d**2)
            print n_calc
            if n_max == n_max_0 or n_min == n_min_0:
                if n_calc > n_guess:
                    n_max = math.ceil(n_calc)
                    n_guess = n_max
                elif n_calc < n_guess:
                    n_min = math.ceil(n_calc)
                    n_guess = n_min
            else:
                if math.ceil(n_calc) == n_guess:
                    break
                if n_calc > n_guess:
                    n_min = max(math.ceil(n_calc),n_guess)
                    n_guess = n_min
                elif n_calc < n_guess:
                    n_max = min(math.ceil(n_calc),n_guess)
                    n_guess = n_max
            print '%s <= n <= %s' %(n_min, n_max)
        n = math.ceil(n_calc)
        print 'A sample of at least size n=%s would need to be taken to estimate the mean to within %s with %s confidence.' %(n,d,confidence_level)

        return


    def samplesize_to_estimate_mean(self, variance,confidence_level,confidence_interval_halfwidth):

        import math
        alpha = float(1-confidence_level)
        d = confidence_interval_halfwidth
        n_min = n_min_0 = 1; n_max = n_max_0 = 1000
        print 'Make an initial of the sample size, n'
        n_guess = int(raw_input())
        while n_min != n_max:
            print 'Look up t(alpha(2)=%s,DF=%s) in Table B.3 on page App19' %(alpha,n_guess-1)
            t = float(raw_input())
            n_calc = variance*t**2/d**2
            print n_calc
            if n_max == n_max_0 or n_min == n_min_0:
                if n_calc > n_guess:
                    n_max = math.ceil(n_calc)
                    n_guess = n_max
                elif n_calc < n_guess:
                    n_min = math.ceil(n_calc)
                    n_guess = n_min
            else:
                if math.ceil(n_calc) == n_guess:
                    break
                if n_calc > n_guess:
                    n_min = max(math.ceil(n_calc),n_guess)
                    n_guess = n_min
                elif n_calc < n_guess:
                    n_max = min(math.ceil(n_calc),n_guess)
                    n_guess = n_max
            print '%s <= n <= %s' %(n_min, n_max)
        n = math.ceil(n_calc)
        print 'A sample of at least size n=%s would need to be taken to estimate the mean to within %s with %s confidence.' %(n,d,confidence_level)

        return


    def samplesize_to_test_hypothesis(self, variance,significance_level,power,mean_difference):

        import math
        alpha = significance_level
        beta = float(1-power)
        delta = mean_difference
        n_min = n_min_0 = 1; n_max = n_max_0 = 1000
        print 'Make an initial of the sample size, n'
        n_guess = int(raw_input())
        while n_min != n_max:
            print 'Look up t(alpha(2)=%s,DF=%s) in Table B.3 on page App19' %(alpha,n_guess-1)
            ta = float(raw_input())
            print 'Look up t(beta(1)=%s,DF=%s) in Table B.3 on page App19' %(beta,n_guess-1)
            tb = float(raw_input())
            n_calc = (variance/delta**2)*(ta+tb)**2
            print 'n=%s' %(round(n_calc,2))
            if n_max == n_max_0 or n_min == n_min_0:
                if n_calc > n_guess:
                    n_max = math.ceil(n_calc)
                    n_guess = n_max
                elif n_calc < n_guess:
                    n_min = math.ceil(n_calc)
                    n_guess = n_min
            else:
                if math.ceil(n_calc) == n_guess:
                    break
                if n_calc > n_guess:
                    n_min = max(math.ceil(n_calc),n_guess)
                    n_guess = n_min
                elif n_calc < n_guess:
                    n_max = min(math.ceil(n_calc),n_guess)
                    n_guess = n_max
            print '%s <= n <= %s' %(n_min, n_max)
        n = math.ceil(n_calc)
        print 'A sample of at least size n=%s would need to be taken to test H0 at a %s significance level with a %s probability of rejecting H0 if the mean is %s larger than the mean of the null hypothesis.' %(n,significance_level,power,delta)

        return


    def akaike(self,):

        '''should be in class regression...'''

        ## unit test
        ## 40.84,27.23,18,1,3
        ## returns
        ## p1 = 0.375281991447
        ## p2 = 0.624718008553
        ## evidence_ratio = 1.66466290094

        k1 = float(pars1+1)
        k2 = float(pars2+1)
        aic1 = n*math.log(ss1/n)+2*k1
        aicc1 = aic1 + (2*k1*(k1+1)) / (n-k1-1)
        aic2 = n*math.log(ss2/n)+2*k2
        aicc2 = aic2 + (2*k2*(k2+1)) / (n-k2-1)
        print aicc1, aicc2
        p1 = math.exp(-.5*(aicc1-aicc2)) / (1+math.exp(-.5*(aicc1-aicc2)))
        p2 = 1-p1
        evidence_ratio = p2/p1
        print p1, p2, evidence_ratio

        return p1, p2, evidence_ratio


    def do_rmsd(self,l_diff,):

        import math

        sum_sqdiff = 0.
        for diff in l_diff:
            sq_diff = diff**2
            sum_sqdiff += sq_diff
        rmsd = math.sqrt(sum_sqdiff/len(l_diff))

        return rmsd


class NumericalRecipes():


    def gammp(self,a,x,):

        '''Returns the incomplete gamma function P(a,x).'''

        ASWITCH = 100. ## When to switch to quadrature method

        if (x < 0.0 or a <= 0.0):
            print "bad args in gammp"
        if (x == 0.0): return 0.0
        elif (int(a) >= ASWITCH): return gammpapprox(a,x,1) ## Quadrature.
        elif (x < a+1.0): return gser(a,x) ## Use the series representation.
        else: return 1.0-self.gcf(a,x) ## Use the continued fraction representation.
        if x < 0 or a <= 0:
            print 'invalid arguments'
            stop

        return


    def gser(self, a, x):
        '''Returns the incomplete gamma function P.a;x/ evaluated by its series representation.'''
        EPS = 0.00000000000000000000000000000001 ## desired relative error
        gln=gammln(a)
        ap=a
        Del=sum=1.0/a
        while True:
            ap+=1
            Del *= x/ap;
            sum += Del
            if (math.fabs(Del) < math.fabs(sum)*EPS):
                return sum*math.exp(-x+a*math.log(x)-gln)


##    def p(x2,nu,):
##
##        '''Given chi squared critical value and DF, return probability density function.'''
##
##        if (x2 <= 0.):
##            print "bad x2 in Chisqdist"
##            stop
##        fac = 0.693147180559945309*(0.5*nu)+gammln(0.5*nu)
##        return math.exp(-0.5*(x2-(nu-2.)*math.log((x2),10))-fac)


    def cdf(self,x2,nu,):

        '''Given chi squared = x2 and DF = nu, return cumulative distribution function.'''

        if (x2 < 0.):
            print "bad x2 in Chisqdist"
            stop

        return self.gammp(0.5*nu,0.5*x2)


    def betai(self,a,b,x):

        ## Numerical Recipes (www.nr.com)
        ## incomplete beta function

        if x < 0. or x > 1.:
            print x
            stop
        if x == 0. or x == 1.:
            bt = 0.0
        else:
            bt = math.exp(
                self.gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*math.log(1.0-x)
                )

        if x < (a+1.)/(a+b+2.):
            f = bt*self.betacf(a,b,x)/float(a)
        else:
            f = 1.-bt*self.betacf(b,a,1.-x)/float(b)

        return f


    def gammln(self,xx):
        
        ## Numerical Recipes (www.nr.com)
        ## gamma function

        coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5,]
        x = xx - 1.0
        tmp = x + 5.5
        tmp = tmp - (x+0.5)*math.log(tmp)
        ser = 1.0
        for j in range(len(coeff)):
            x += 1
            ser += coeff[j]/x

        return -tmp + math.log(2.50662827465*ser)


    def betacf(self,a,b,x):

        ## Numerical Recipes (www.nr.com)

        ITMAX = 200
        EPS = 3.0e-7

        bm = az = am = 1.
        qab = a+b
        qap = a+1.
        qam = a-1.
        bz = 1.-qab*x/qap
        for i in range(ITMAX+1):
            em = float(i+1)
            tem = em + em
            d = em*(b-em)*x/((qam+tem)*(a+tem))
            ap = az + d*am
            bp = bz+d*bm
            d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
            app = ap+d*az
            bpp = bp+d*bz
            aold = az
            am = ap/bpp
            bm = bp/bpp
            az = app/bpp
            bz = 1.
            if (abs(az-aold)<(EPS*abs(az))):
                break

        if i == ITMAX:
            print 'a or b too big, or ITMAX too small in Betacf.'
            stop

        return az


    def gcf(self, a, x):

        '''Returns the incomplete gamma function Q.a; x/ evaluated by its continued fraction representation.'''
        '''Also sets ln a as gln. User should not call directly.'''
    ##    Int i;
    ##    Doub an,b,c,d,del,h
        EPS = 0.00000000000000000000000000000001 ## desired relative error
        FPMIN = 0.00000000000000000000000000000001 ## number close to the smallest representable floating-point number
        gln=self.gammln(a)
        b=x+1.0-a ## Set up for evaluating continued fraction by modified Lentz's method (5.2) with b0 = 0.
        c=1.0/FPMIN
        d=1.0/b
        h=d
        for i in xrange(1,99999999,1): ## Iterate to convergence.
            an = -i*(i-a)
            b += 2.0
            d=an*d+b
            if (math.fabs(d) < FPMIN): d=FPMIN
            c=b+an/c
            if (math.fabs(c) < FPMIN): c=FPMIN
            d=1.0/d
            Del=d*c
            h *= Del
            if (math.fabs(Del-1.0) <= EPS):
                break
        return math.exp(-x+a*math.log(x)-gln)*h ## Put factors in front.


class correlation_and_regression():

    def twofactor_anova(
        ## 2x2 matrix
        lr1c1,lr1c2,lr2c1,lr2c2,
        verbose = True
        ):

    #### faculty.vassar.edu/lowry/anova2u
    ##l11 = [20.4,17.4,20,18.4,24.5,21,19.7,22.3,17.3,23.3]
    ##l12 = [20.5,26.3,26.6,19.8,25.4,28.2,22.6,23.7,22.5,22.6]
    ##l21 = [22.4,19.1,22.4,25.4,26.2,25.1,28.8,21.8,26.3,25.2]
    ##l22 = [17.5,13.6,16.9,12.4,16.4,18.3,13.6,19.1,16.1,20.5]


        l11 = lr1c1
        l12 = lr1c2
        l21 = lr2c1
        l22 = lr2c2

        if len(l11) != len(l12):
            stop1
        if len(l12) != len(l21):
            stop2
        if len(l21) != len(l22):
            stop3

        n = len(l11)
        rows = 2
        cols = 2
        l_n = [10,10,10,10]

    ##    b = len(groups[0][0]) ## replication
    ##    a = len(groups) ## number of rows
    ##    c = len(groups[0]) ## number of cols

        l_sumx = []
        l_sumxx = []
        for l in [l11,l12,l21,l22,]:
            sumx, sumxx = calc_sum_of_squares(l)
            l_sumx += [sumx]
            l_sumxx += [sumxx]

        ## machine formulas from www4.uwsp.edu/psych/stat/13/anova-2w.htm
        I = sum(l_sumx)**2/(4*n)
        II = sum(l_sumxx)
        ## rows
        III = (
            ## row1
            ((l_sumx[0]+l_sumx[1])**2)/(2*n)
            +
            ## row2
            ((l_sumx[2]+l_sumx[3])**2)/(2*n)
            )
        ## cols
        IV = (
            ## col1
            ((l_sumx[0]+l_sumx[2])**2)/(2*n)
            +
            ## col2
            ((l_sumx[1]+l_sumx[3])**2)/(2*n)
            )
        V = sum(numpy.array(l_sumx)**2)/n
    ##    if verbose == True:
    ##        print 'I)', I
    ##        print 'II), II
    ##        print 'III)', III
    ##        print 'IV)', IV
    ##        print 'V)', V

        ## rows
        SS_A = SS_rows = (
            (l_sumx[0]+l_sumx[1])**2/(2*n)
            +
            (l_sumx[2]+l_sumx[3])**2/(2*n)
            -
            I
            )
        SS_A = III-I
        if round(SS_A,8) == 0:
            SS_A = 0
        DF_A = DF_rows = 2-1
        ## cols
        SS_B = (
            (l_sumx[0]+l_sumx[2])**2/(2*n)
            +
            (l_sumx[1]+l_sumx[3])**2/(2*n)
            -
            I
            )
        SS_B = IV-I
        if round(SS_B,8) == 0:
            SS_B = 0
        DF_B = 2-1

        ## interaction
    ##    SS_interaction = SS_RxC = SS_AxB = SS_rxc = SS_cells-SS_A-SS_B
        SS_interaction = SS_AxB = V+I-III-IV
        if round(SS_AxB,8) == 0:
            SS_AxB = 0
        DF_AxB = (2-1)*(2-1)

        ## within groups (error)
    ##    SS_error = SS_wg = sum(l_SS)
        SS_error = II-V
        DF_error = rows*cols*n-rows*cols

        ## total
    ##    SS_total = sum(l_sumxx)-(sum(l_sumx)**2)/sum(l_n)
        SS_total = II-I
        DF_total = rows*cols*n-1

    ##    ## between groups (not used)
    ##    SS_bg = SS_total - SS_error
    ##    DF_bg = 2*2-1
    ##    print 'SS bg', SS_bg
    ##    print SS_total, SS_error
    ##    stop

    ##    var_pooled = (ss1+ss2)/(n1+n2-2)
    ##    stderr = stddev = math.sqrt(var_pooled/n1+var_pooled/n2)

        MS_A = SS_A/DF_A
        MS_B = SS_B/DF_B
        MS_AxB = SS_AxB/DF_AxB
        MS_error = SS_error/DF_error

        F_A = MS_A/MS_error
        F_B = MS_B/MS_error
        F_AxB = MS_AxB/MS_error

        l_means = numpy.array(l_sumx)/n

        l_SS = []
        l_variances = []
        for i in range(len(l_sumx)):
            sumx = l_sumx[i]
            sumxx = l_sumxx[i]
            SS = sumxx-(sumx**2)/n
            l_SS += [SS]
            l_variances += [SS/n]

        try:
            p_A = fdist(F_A,DF_A,DF_error)
            p_B = fdist(F_B,DF_B,DF_error)
            p_interaction = fdist(F_AxB,DF_AxB,DF_error)
        except:
    ##        for l in [l11,l12,l21,l22,]:
    ##            print
    ##            for v in l:
    ##                print v
            print 'SS_AxB = V+I-III-IV', V,I,III,IV
            print V+I
            print III+IV
            print round(SS_AxB,12)
            print 'FA', F_A
            print 'FB', F_B
            print 'FAxB', F_AxB
            p_A = fdist(F_A,DF_A,DF_error)
            p_B = fdist(F_B,DF_B,DF_error)
            p_interaction = fdist(F_AxB,DF_AxB,DF_error)
            stop
    ##        return None, None, None

        if verbose == True:
            print
            print 'SS_A', SS_A
            print 'SS_B', SS_B
            print 'SS_interaction', SS_AxB
            print 'SS error/within', SS_error
            print 'SS total', SS_total
            print
            print 'DF_A', DF_A
            print 'DF_B', DF_B
            print 'DF_AxB', DF_AxB
            print 'DF_error', DF_error
            print
            print 'MS_A', MS_A
            print 'MS_B', MS_B
            print 'MS_AxB', MS_AxB
            print 'MS_error', MS_error
            print
            print 'F_A', F_A
            print 'F_B', F_B
            print 'F_AxB', F_AxB
            print
            print p_A
            print p_B
            print p_interaction
            print
            print 'n', n

        if verbose == True:
            print 'means', l_means
            print 'variances', l_variances

        return p_A, p_B, p_interaction, l_means, l_variances


class tests():


    def Ztest(data,mean,stddev,limits):

        for limit in limits:
            Z = (limit-mean)/stddev
            P = normaldist(stddev, limit, abs(Z))
            if Z > 0:
                print 'P( X > %s ) = P( Z > %4.2f.. ) = %6.4f%s' %(limit, Z, 100*P, '%')
            if Z < 0:
                print 'P( X < %s ) = P( Z < %4.2f.. ) = %6.4f%s' %(limit, Z, 100*P, '%')
        print 'my', my
        print 'sigma', sigma

        return Z


    def chi_square(self,matrix_count_observed,):

        '''The test is appropriate when the following conditions are met:
    The sampling method is simple random sampling.
    Each population is at least 10 times as large as its respective sample.
    The variables under study are each categorical.
    If sample data are displayed in a contingency table, the expected frequency count for each cell of the table is at least 5.
    '''

        rows,cols = numpy.shape(matrix_count_observed)
        ## calculate degrees of freedom
        DF = (rows-1)*(cols-1)
        ## calculate row and column sums
        l_sum_col = matrix_count_observed.sum(axis=0)
        l_sum_row = matrix_count_observed.sum(axis=1)
        sum_matrix = sum(l_sum_row)
        ## calculate expected frequency count
    ##    matrix_count_expected = numpy.zeros((rows,cols))
        matrix_count_ratio = numpy.zeros((rows,cols))
        for row in range(rows):
            sum_row = l_sum_row[row]
            for col in range(cols):
                sum_col = l_sum_col[col]
                count_expected = sum_row*sum_col/sum_matrix
    ##            matrix_count_expected[row][col] = count_expected
                matrix_count_ratio[row][col] = (count_expected-matrix_count_observed[row][col])**2/count_expected
        test_statistic = chi_squared_critical_value = matrix_count_ratio.sum()
        ## Numerical Recipes 3rd, 6.14.8 Chi-Square Distribution, p. 330
        ## The chi-square distribution is actually just a special case of the gamma distribution, below,
        ## so its cdf is given by an incomplete gamma function P.a;x/,
        print 'chi2 squared critical value', chi_squared_critical_value
        instance_NumericalRecipes = NumericalRecipes()
        p = probability_cumulative = 1-instance_NumericalRecipes.cdf(chi_squared_critical_value,DF,)
        print 'p', p
        ## decline hypothesis that there is no relationship between row and column parameters if probability less than significance level

        return p


    def tdist(self, t,df):
    
        '''return probability of t'''

    ##    ## http://mathworld.wolfram.com/Studentst-Distribution.html (eq. 7)
    ##    import math
    ##    z = float(df)/(float(df)+float(t)**2)
    ##    a = .5*float(df)
    ##    b = .5
    ##    lim1 = 0.
    ##    d_lim2 = {'numerator':z,'denominator':1.}
    ##    inv_stepsize = 1000000.
    ##    function_beta_regularized = {}
    ##    for key in d_lim2:
    ##        lim2 = d_lim2[key]
    ##        deltax = (lim2-lim1)/inv_stepsize
    ##        area = 0.
    ##        x = lim1
    ##        while x <= lim2:
    ##            y = ((x)**(a-1)) * ((1-x)**(b-1))
    ##            area += y*deltax
    ##            x += deltax
    ##        function_beta_regularized[key] = area
    ##    p = (function_beta_regularized['numerator']/function_beta_regularized['denominator'])

        t = float(t)
        df = float(df)
        p = betai(df/2., .5, df/(df+t**2))

        return p


    def fdist(F,df_n,df_d):

        a = df_n/2.
        b = df_d/2.
        x = df_n*F/float(df_n*F+df_d)

        p = 1-betai(a,b,x,)

    #### http://mathworld.wolfram.com/F-Distribution.html (eq. 4)
    #### http://mathworld.wolfram.com/RegularizedBetaFunction.html
    #### http://mathworld.wolfram.com/IncompleteBetaFunction.html (eq. 1)
    #### http://mathworld.wolfram.com/BetaFunction.html (eq. 18)
    ##    import math
    ##    x = df_n*F/float(df_d+df_n*F)
    ##    a = df_n/2.
    ##    b = df_d/2.
    ##    lim1 = 0.
    ##    d_lim2 = {'numerator':x,'denominator':1.}
    ##    inv_stepsize = 1000000.
    ##    d_function_beta_regularized = {}
    ##    for key in d_lim2:
    ##        lim2 = d_lim2[key]
    ##        dt = (lim2-lim1)/inv_stepsize
    ##        area = 0.
    ##        t = lim1+dt
    ##        while t <= lim2:
    ##            area += ((t)**(a-1)) * ((1-t)**(b-1))*dt
    ##            t += dt
    ##        d_function_beta_regularized[key] = area
    ##    p = 1-d_function_beta_regularized['numerator']/d_function_beta_regularized['denominator']

        return p


    def onesamplettest(data,mean):

        '''return p - tdist wrapper'''

        import math
        n = len(data)
        print 'n', n
        sample_mean = float(sum(data))/len(data)
        print 'sample_mean', sample_mean
        sample_var = sum([(float(x)-sample_mean)**2 for x in data])/(n-1)
        print 'sample_var', sample_var
        sample_stderr = math.sqrt(sample_var/n)
        print 'sample_stderr', sample_stderr
        t = (sample_mean-mean)/sample_stderr ## rewrite and calculate mean/my from critical t values
        p = self.tdist(t,n-1)
        print 'my', mean
        print 't', t
        print 'p', p

        return t


    def twosamplettest(l1,l2,verbose=True,):

        '''assume 1) normal distributions, 2) equal sample sizes, 3) equal variances'''

        import math

        d_statistics = {
            1:{'sample':l1},
            2:{'sample':l2},
            }

        for i_sample in d_statistics.keys():
            l_sample = d_statistics[i_sample]['sample']
            n = len(l_sample)
            sumx = 0.
            sumxx = 0.
            for i in range(len(l_sample)):
                x = float(l_sample[i])
                sumx += x
                sumxx += x**2
            SS = sumxx-(sumx**2)/n
            mean = sumx/n
    ##        ss = sum([(float(x)-mean)**2 for x in l_sample])
            d_statistics[i_sample]['mean'] = mean
            d_statistics[i_sample]['SS'] = SS
            d_statistics[i_sample]['n'] = n

        n1 = d_statistics[1]['n']
        n2 = d_statistics[2]['n']
        ss1 = d_statistics[1]['SS']
        ss2 = d_statistics[2]['SS']
        mean1 = d_statistics[1]['mean']
        mean2 = d_statistics[2]['mean']

        mean = 0
        var_pooled = (ss1+ss2)/(n1+n2-2)
        stderr = stddev = math.sqrt(var_pooled/n1+var_pooled/n2)
        t = ((mean1-mean2)-mean)/stddev
        p = self.tdist(t,n1+n2-2)

        if verbose == True:
            print 'n', n1, n2
            print 'mean', mean1, mean2
            print 'var', ss1/n1, ss2/n2, 'pooled sp2', var_pooled
            print 'std err of the diff between the sample means, s(x1-x2)', stderr
            print 't', round(t,3)
            print 'p', p

        return mean1,mean2,stderr,t,p


    def Kolmogorov_Smirnov(l1,l2,verbose=True,):

        l = l1+l2
        l.sort()

        l_set = list(set(l))
        l_set.sort()

        n1 = float(len(l1))
        n2 = float(len(l2))

        Fi1 = 0
        Fi2 = 0
        D = 0
        for i in range(len(l_set)):

            xi = l_set[i]
            ## frequency
            if xi in l1:
                fi1 = l1.count(xi)
            else:
                fi1 = 0
            if xi in l2:
                fi2 = l2.count(xi)
            else:
                fi2 = 0
            ## cumulative frequency
            Fi1 += fi1
            Fi2 += fi2
            ## cumulative relative frequency
            rel_Fi1 = Fi1/n1
            rel_Fi2 = Fi2/n2
    ##        print rel_Fi1, rel_Fi2
            ## KS statistic
            Di = abs(rel_Fi1-rel_Fi2)
            if Di > D:
                D = Di

        D_crit = 1.36*math.sqrt((n1*n2)/(n1+n2))
        if verbose == True:
            print 'D', D
            print 'D_crit', D_crit

        ## THIS IS UTTERLY WRONG!!!
        if D < D_crit:
            p = 0.049
        else:
            p = 0.051

        return p


    def MannWhitney(l1,l2,):

        '''Two sample test. Nonparametric analogue to the two-sample t test. Need to be fixed...'''
        d_ranks = {}
        val = None
        ## sort data
        data_pooled = l1+l2
        data_pooled.sort()
        ## rank data
        for i in range(len(data_pooled)):
            val = data_pooled[i]
            if not d_ranks.has_key(val):
                d_ranks[val] = (i+1+i+data_pooled.count(val))/2.
        print d_ranks
        ## sum ranks
        R = [0.,0.]
        for l in [l1,l2,]:
            for val in l:
                R[i] += d_ranks[val]
        print 'summed ranks R', R
        U = []
        for i in range(2):
            l = [l1,l2,][i]
            U += [len(l1)*len(l2)+len(l)*(len(l)+1)/2.-R[i]]
        print 'MannWhitney statistic U = n1n2 + n1(n1+1)/2 - R1',U

        return U


    def do_regression(self,l1,l2,verbose=True,):

        ## fit to f(x) = b*x + a

        n = len(l1)

        sum_xy = 0
        sum_x = 0
        sum_y = 0
        sum_xx = 0
        sum_yy = 0
        for i in range(n):
            sum_xy += l1[i]*l2[i]
            sum_x += l1[i]
            sum_y += l2[i]
            sum_xx += l1[i]**2
            sum_yy += l2[i]**2
        SS_x = sum_xx-(sum_x**2)/float(n)
        SS_y = sum_yy-(sum_y**2)/float(n)
        ## sum of cross products (17.3)
        SS_xy = sum_of_cross_products = sum_xy-(sum_x*sum_y)/float(n)
        mean_x = sum_x/float(n)
        mean_y = sum_y/float(n)

        r = (sum_xy-sum_x*sum_y/n)/math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))

        ## regression coefficient (17.4) ... slope
        b = SS_xy/SS_x
        ## intercept (17.7)
        a = mean_y-b*mean_x

        ## is there a machine formula?
        SE = s_b = (
            math.sqrt( sum([(l2[i]-(a+b*l1[i]))**2 for i in range(len(l1))]) / (n-2) )
            /
            math.sqrt(sum([(l1[i]-mean_x)**2 for i in range(len(l1))]))
            )

        ## Student's t-test
        DF = n-2
        t_slope = b/SE
        p_slope = self.tdist(t_slope,DF)
        t_correlation = r*math.sqrt((n-2)/(1-r**2)) ## same as t_slope?
        p_correlation = self.tdist(t_correlation,DF)

        if verbose == True:
            print 'y = a+b*x'
            print 'n', n
            print 'r', r
            print 'r^2', r**2
            print 'b', b
            print 'a', a
            print 't', t_correlation
            print 'p', p_correlation

        return a,b,r,p_correlation


    def correlation(self,l1,l2):

        import math

        if len(l1) != len(l2):
            stop
        n = len(l1)

        sum_xy = 0
        sum_xx = 0
        sum_yy = 0
        sum_x = 0
        sum_y = 0
        for i in range(n):
            x = l1[i]
            y = l2[i]
            sum_xy += x*y
            sum_xx += x*x
            sum_yy += y*y
            sum_x += x
            sum_y += y
        r = (sum_xy-sum_x*sum_y/n)/math.sqrt((sum_xx-sum_x**2/n)*(sum_yy-sum_y**2/n))

        return r


if __name__ == '__main__':

    '''Two sample tests:
chi_square
do_regression
do_confidence_bands
correlation
MannWhitney
Wilcoxon
twosamplettest
Kolmogorov_Smirnov

One sample tests:
Ztest
onesamplettest

Calculate distributions:


Calculate statistics:
do_stddev
do_stderr
do_rmsd

Other:
samplesize_to_estimate_mean_difference
samplesize_to_estimate_mean
samplesize_to_estimate_hypothesis
'''


    ## chi squared test
    count_matrix = numpy.zeros((2,3))
    count_matrix[0][0] = 200
    count_matrix[0][1] = 150
    count_matrix[0][2] = 50
    count_matrix[1][0] = 250
    count_matrix[1][1] = 300
    count_matrix[1][2] = 50
    instance = tests()
    instance.chi_square(count_matrix)


    ## http://en.wikipedia.org/wiki/File:Method_example_for_GWA_study_designs.png
    import time
    t1 = time.time()
    count_matrix = numpy.zeros((2,2))
    count_matrix[0][0] = 2104.
    count_matrix[0][1] = 4000-2104.
    count_matrix[1][0] = 2676.
    count_matrix[1][1] = 6000-2675.
    instance = tests()
    instance.chi_square(count_matrix)
    t2 = time.time()
    print 1000000*(t2-t1)/3600


    ## regression
    l1 = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,]
    l2 = [1.1,1.9,3.1,3.9,5.1,5.9,7.1,7.9,9.1,9.9,]
##    l1 = [95,85,80,70,60,]
##    l2 = [85,95,70,65,70,]
    l1 = [207,180,220,205,190,]
    l2 = [6907,5991,6810,6553,6190,]
##    instance = 
    do_regression(l1,l2)

    pass

import math

#!/usr/bin/python
def main():

    return


def betai(a,b,x):

    ## Numerical Recipes (www.nr.com)
    ## incomplete beta function

    if x < 0. or x > 1.:
        stop
    if x == 0. or x == 1.:
        bt = 0.0
    else:
        bt = math.exp(
            gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*math.log(1.0-x)
            )

    if x < (a+1.)/(a+b+2.):
        f = bt*betacf(a,b,x)/float(a)
    else:
        f = 1.-bt*betacf(b,a,1.-x)/float(b)

    return f


def gammln(xx):
    
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


def betacf(a,b,x):

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


def stderr(l):

    n = float(len(l))
    sumx = 0
    sumxx = 0
    for x in l:
        sumx += x
        sumxx += x*x
    average = sumx/n
    SS = sumxx-(sumx**2)/n
    var = SS/(n-1)
    stddev = math.sqrt(var)
    stderr = math.sqrt(var/n)

    return average, stddev


def rmsd(l_diff,):

    import math

    sum_sqdiff = 0.
    for diff in l_diff:
        sq_diff = diff**2
        sum_sqdiff += sq_diff
    rmsd = math.sqrt(sum_sqdiff/len(l_diff))

    return rmsd

def correlation(l1,l2):

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


def Wilcoxon(data):

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


def MannWhitney(data):

    '''Two sample test. Nonparametric analogue to the two-sample t test.'''
    d_ranks = {}
    val = None
    ## sort data
    data_pooled = data[0]+data[1]
    data_pooled.sort()
    ## rank data
    for i in range(len(data_pooled)):
        val = data_pooled[i]
        if not d_ranks.has_key(val):
            d_ranks[val] = (i+1+i+data_pooled.count(val))/2.
    print d_ranks
    ## sum ranks
    R = [0.,0.]
    for i in range(2):
        for val in data[i]:
            R[i] += d_ranks[val]
    print 'summed ranks R', R
    U = []
    for i in range(2):
        U += [len(data[0])*len(data[1])+len(data[i])*(len(data[i])+1)/2.-R[i]]
    print 'MannWhitney statistic U = n1n2 + n1(n1+1)/2 - R1',U

    return U


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


def twosamplettest(l1,l2,verbose=True,):

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
    p = tdist(t,n1+n2-2)

    if verbose == True:
        print 'n', n1, n2
        print 'var', ss1/n1, ss2/n2, 'pooled sp2', var_pooled
        print 'mean', mean1, mean2
        print 't', t
        print 'p', p
        print 'standard error of he difference between the sample means s(x1-x2)', stderr

    return mean1,mean2,stderr,p


def onesamplettest(data,mean):

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
    p = tdist(t,n-1)
    print 'my', mean
    print 't', t
    print 'p', p

    return t


def samplesize_to_estimate_mean_difference(variance_pooled,confidence_level,confidence_interval_halfwidth):

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


def samplesize_to_estimate_mean(variance,confidence_level,confidence_interval_halfwidth):

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


def samplesize_to_test_hypothesis(variance,significance_level,power,mean_difference):

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


def normaldist(pop_stddev, pop_mean, Zlimit):
    
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


def faculty(n):

    n = int(n)
    fac = 1
    for i in range(1,n+1):
        fac *= i

    return fac


def tdist(t,df):

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
    p = betai(df/2.,.5,df/(df+t**2))

    return p


if __name__ == '__main__':
    main()

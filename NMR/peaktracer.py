#!/bin/env /software/bin/python
# $Id: peaktracer.py 212 2007-04-03 14:58:04Z tc $
# Copyright (C) Tommy Carstensen, University College Dublin, 2007

## identify peaks based on
## 1a) proton chemical shift !
## 1b) nitrogen chemical shift !
## 2a) peak width in proton dimension ## change if overlap!
## 2b) peak width in nitrogen dimension ## change if overlap!
## 3) peak height ## change if overlap!
## 4) peak movement direction ...

## when peaks split in two after being overlapped then get overlap from previous spectra!

## reassign peaks that move too far because assigned as part of a group (D48)

## dont use for titration while overlapping

## 7) if isolated unassigned ref peaks present then change contour scale
## 8) if grouped unassigned ref peaks present then allow overlap (multiple sample peaks assigned to one ref peak)
##     mention in dic how many assignments are made per peak (e.g. asscount)
##     move peaks slightly according to direction if new overlap or according to previous pos if previosly also overlap

## allow peak dissapearance if previosly low S/N (i.e. < 28) max value for D66 at pH697

def main(
    pHs_and_peaklists,
    pH_and_peaklist_reference,
    max_diff_w1_move = 0.63, ## max value observed for M12 between pH5.02 and pH5.55
    max_diff_w2_move = 0.079, ## max value observed for T40 between pH6.97 and pH7.45
    max_diff_w1_group = 0.50, ## 0.28
    max_diff_w2_group = 0.08,
    max_diff_w1_peakoverlap = 0.50,
    max_diff_w2_peakoverlap = 0.06,
    min_SN = 10, ## 25, 20
    ratio_concentration = 1.0, ## use for peak height normalization if different NMR samples...
    max_height_diff = 0.125, ## max percentage difference between peak heights not due to an overlap (R14N-HN in HEWL)
    ):

    import os, math

    peaklist_reference = pH_and_peaklist_reference.values()[0]

    pHs = pHs_and_peaklists.keys()
    pHs.sort()

    peaklists = []
    for pH in pHs:
        peaklists += [pHs_and_peaklists[pH]]
    index = peaklists.index(peaklist_reference)
    peaklists2 = peaklists[:index+1]
    peaklists2.reverse()
    peaklists1 = peaklists[index:]

    if peaklist_reference not in pHs_and_peaklists.values():
        raise 'peaklist_reference not in pHs_and_peaklists'

##    ###########################
##    ## parse chemical shifts ##
##    ###########################
##    peaks = {}
##    for peaklist in peaklists:
##        peaks = parsepeaklist(peaks, peaklist, peaklist_reference, min_SN)
##
##    #################
##    ## trace peaks ##
##    #################
##    peaks = tracepeaks(
##        peaks, peaklists1,peaklists2,
##        max_diff_w1_move, max_diff_w2_move,
##        max_diff_w1_group, max_diff_w2_group,
##        max_diff_w1_peakoverlap, max_diff_w2_peakoverlap,
##        max_height_diff,
##        )

    ################################################################
    ## fit chemical shift curve to HH eq asssuming rapid exchange ##
    ################################################################

    ##
    ## parse chemical shifts from the assigned peaklists (peaks dictionary)
    ##

    chemicalshifts = {}
##    lines = ['##pH']
    for pH in pHs:
##        assignments = lines[0].split()[1:]
##        lines += ['%4.2f' %(pH)]
        fd = open('ass%s.list' %(int(round(100*pH))))
        lines = fd.readlines()[2:]
        fd.close()
        for line in lines:
            ass = line.split()[0]
            if ass == '?-?':
                continue

            if ass not in chemicalshifts.keys():
                chemicalshifts[ass] = {}
            chemicalshifts[ass][pH] = [float(line.split()[1]),float(line.split()[2])]

    ##
    ## calculate correlation coefficient
    ##
##    linear = fit_linear(chemicalshifts, pHs)
    linear95 = ['A107N-HN', 'A95N-HN', 'L17N-HN', 'G71N-HN', 'T43N-HN', 'G49N-HN', 'T40N-HN', 'Q57N-HN', 'T51N-HN', 'K97N-HN', 'A90N-HN', 'G22N-HN', 'F34N-HN', 'G102N-HN', 'W62N-HN', 'R128N-HN', 'S50N-HN', 'G16N-HN', 'W28NE1-HE1', 'R14N-HN', 'H15N-HN', 'N27N-HN', 'D101N-HN', 'Y20N-HN', 'E35N-HN', 'R73N-HN', 'W108NE1-HE1', 'L75N-HN', 'L83N-HN', 'F3N-HN', 'C94N-HN', 'W63N-HN', 'A10N-HN', 'L129N-HN', 'R21N-HN', 'T118N-HN', 'A42N-HN', 'L56N-HN', 'V29N-HN', 'N93N-HN', 'D66N-HN', 'K96N-HN', 'K33N-HN', 'D119N-HN', 'A31N-HN', 'R112N-HN', 'N39N-HN', 'S86N-HN', 'D18N-HN', 'W62NE1-HE1', 'R125N-HN', 'C115N-HN', 'G54N-HN', 'L8N-HN', 'Q121N-HN', 'F38N-HN', 'D52N-HN', 'S100N-HN', 'M12N-HN', 'R114N-HN', 'V92N-HN', 'D48N-HN', 'W111NE1-HE1', 'N59N-HN', 'I98N-HN', 'L84N-HN', 'W123N-HN', 'N106N-HN', 'V109N-HN', 'S72N-HN', 'R45N-HN', 'C76N-HN', 'A11N-HN', 'N19N-HN', 'W28N-HN', 'G117N-HN', 'G67N-HN', 'W111N-HN', 'N77N-HN', 'A110N-HN', 'V99N-HN', 'N74N-HN', 'I78N-HN', 'E7N-HN', 'W63NE1-HE1', 'S85N-HN', 'S36N-HN', 'R5N-HN', 'G4N-HN', 'A32N-HN', 'G104N-HN', 'Q41N-HN', 'N113N-HN', 'A122N-HN', 'S91N-HN', 'T89N-HN']
    linear99 = ['A107N-HN', 'A95N-HN', 'L17N-HN', 'G71N-HN', 'T43N-HN', 'G49N-HN', 'T40N-HN', 'Q57N-HN', 'T51N-HN', 'K97N-HN', 'A90N-HN', 'G22N-HN', 'F34N-HN', 'G102N-HN', 'W62N-HN', 'R128N-HN', 'S50N-HN', 'G16N-HN', 'W28NE1-HE1', 'R14N-HN', 'H15N-HN', 'N27N-HN', 'D101N-HN', 'E35N-HN', 'W108NE1-HE1', 'L83N-HN', 'F3N-HN', 'C94N-HN', 'W63N-HN', 'A10N-HN', 'L129N-HN', 'R21N-HN', 'A42N-HN', 'L56N-HN', 'V29N-HN', 'N93N-HN', 'D66N-HN', 'K96N-HN', 'K33N-HN', 'D119N-HN', 'A31N-HN', 'R112N-HN', 'N39N-HN', 'S86N-HN', 'D18N-HN', 'W62NE1-HE1', 'R125N-HN', 'C115N-HN', 'L8N-HN', 'Q121N-HN', 'F38N-HN', 'D52N-HN', 'S100N-HN', 'M12N-HN', 'R114N-HN', 'V92N-HN', 'D48N-HN', 'W111NE1-HE1', 'N59N-HN', 'I98N-HN', 'L84N-HN', 'W123N-HN', 'N106N-HN', 'V109N-HN', 'S72N-HN', 'R45N-HN', 'C76N-HN', 'A11N-HN', 'N19N-HN', 'G67N-HN', 'W111N-HN', 'N77N-HN', 'V99N-HN', 'N74N-HN', 'E7N-HN', 'W63NE1-HE1', 'S85N-HN', 'S36N-HN', 'G4N-HN', 'A32N-HN', 'Q41N-HN', 'N113N-HN', 'A122N-HN', 'S91N-HN', 'T89N-HN']
    linear = linear99

    ##
    ## write the chemical shifts of different assignments at differnt pH values to file
    ##

    l_ass = chemicalshifts.keys()
    l_ass = linear
    l_ass.sort()

    lines = ['##1/pH']
    for pH in pHs:
        lines += ['%4.2f' %(pH)]
    i = 2
    for ass in l_ass:
        lines[0] += ' %s/%s' %(i, ass)
        for j in range(len(pHs)):
            pH = pHs[j]
            if pH not in chemicalshifts[ass].keys():
                lines[j+1] += '     N/A'
            else:
                w1 = chemicalshifts[ass][pH][0]
                w2 = chemicalshifts[ass][pH][1]
                cs = math.sqrt(w1**2+(10*w2)**2)
                lines[j+1] += ' %7.3f' %(cs)
        i += 1

    lines[0] += '\n'
    for j in range(len(pHs)):
        lines[j+1] += '\n'

    fd = open('pH.gnudata', 'w')
    fd.writelines(lines)
    fd.close()

    ##
    ## plot and fit
    ##

    for i in range(len(l_ass)):
        ass = l_ass[i]
        fileprefix = ass[:ass[1:].index('N')+1]
        fileprefix = ass
        lines = [
##            'f(x) = (  cs_prot+cs_deprot*10**( nH*(x-pKa) )  ) / (  1+cs_deprot*10**( nH*(x-pKa) )  )\n',
            'f(x) = (  cs_prot+cs_deprot*10**( (x-pKa) )  ) / (  1+cs_deprot*10**( (x-pKa) )  )\n',

            'cs_prot = %s\n' %(chemicalshifts[ass][min(chemicalshifts[ass].keys())][0]),
            'cs_deprot = %s\n' %(chemicalshifts[ass][max(chemicalshifts[ass].keys())][0]),
            'pKa = 7\n',
##            'nH = 1\n',

##            'fit f(x) "%s.gnudata" u 1:2 via pKa, cs_prot, cs_deprot, nH\n' %(fileprefix),
            'fit f(x) "%s.gnudata" u 1:2 via pKa, cs_prot, cs_deprot\n' %(fileprefix),

            'set terminal png\n',
            'set output "%s.png"\n' %(ass),
            'plot "pH.gnudata" u 1:%s t "%s", f(x)\n' %(i+2,ass),
            ]

        fd = open('gnuscript.txt','w')
        fd.writelines(lines)
        fd.close()
        
##        os.system('gnuplot gnuscript.txt')

    return


def fit_linear(chemicalshifts, pHs):

    import math

    linear99 = []
    linear95 = []
    
    for ass in chemicalshifts.keys():
        sumxiyi = 0.
        sumxixi = 0.
        sumyiyi = 0.
        sumxi = 0.
        sumyi = 0.
        n = 0.
        for pH in chemicalshifts[ass]:
            n += 1
            x = chemicalshifts[ass][pH][0]
            y = chemicalshifts[ass][pH][1]
            sumxiyi += x*y
            sumxixi += x*x
            sumyiyi += y*y
            sumxi += x
            sumyi += y
        if n > 1:
            sumxy = sumxiyi - sumxi*sumyi/n
            sumxx = sumxixi - sumxi*sumxi/n
            sumyy = sumyiyi - sumyi*sumyi/n
            r = abs(sumxy / math.sqrt(sumxx*sumyy))
            SE = math.sqrt((1-r**2)/(n-2))
            t = r/SE
            p = tdist(t,n-2)
            print ass, p##, t, r, n
            if p < 0.01:
                linear99 += [ass]
            if p < 0.05:
                linear95 += [ass]

    print linear95
    print linear99
                
    return linear95


def tdist(t,df):

## http://mathworld.wolfram.com/Studentst-Distribution.html (eq. 7)

    import math
    z = float(df)/(float(df)+float(t)**2)
    a = .5*float(df)
    b = .5
    lim1 = 0.
    d_lim2 = {'numerator':z,'denominator':1.}
    inv_stepsize = 1000000.
    function_beta_regularized = {}

    for key in d_lim2:

        lim2 = d_lim2[key]

        deltax = (lim2-lim1)/inv_stepsize
        area = 0.

        x = lim1
        while x <= lim2:
            y = ((x)**(a-1)) * ((1-x)**(b-1))
            area += y*deltax
            x += deltax

        function_beta_regularized[key] = area

    P = (function_beta_regularized['numerator']/function_beta_regularized['denominator'])
    
    return P


def tracepeaks(
    peaks, peaklists1, peaklists2,
    max_diff_w1_move, max_diff_w2_move,
    max_diff_w1_group, max_diff_w2_group,
    max_diff_w1_peakoverlap, max_diff_w2_peakoverlap,
    max_height_diff,
    ):

    ##
    ## main loop over spectra
    ##
    for peaklists in [peaklists1, peaklists2]:
        for i in range(len(peaklists)-1):
            peaklist1 = peaklists[i]
            for j in range(i+1,i+2):
                peaklist2 = peaklists[j]

                print peaklist1, peaklist2

                ## 1) group sample peaks
                peaks = grouppeaks(peaks, peaklist2, max_diff_w1_group, max_diff_w2_group)

                ## 2) find reference peaks in the vicinity of the grouped sample peaks
                peaks = find_assignments_vicinal_to_grouped_peaks(peaks, peaklist1, peaklist2, [max_diff_w1_move,max_diff_w2_move])

                ## 3) assign grouped sample peaks to reference peaks in the vicinity
                peaks = assignment(peaks, peaklist1, peaklist2)

                ## 4) check for reference peaks with multiple sample peaks assigned to them
                peaks = multiassignments(peaks, peaklist1, peaklist2)

                ## 5) assign nonassigned sample peaks to reference peaks
                peaks = assign_nonassigned_peaks(peaks, peaklist1, peaklist2, max_diff_w1_move, max_diff_w2_move)

                ## 6) assign sample peaks to nonassigned reference peaks
    ##            for peak in peaks[peaklist1]:
    ##                if peaks[peaklist1][peak]['w1'] == 118.819:
    ##                    print peak, peaks[peaklist1][peak]
    ##            stop
                peaks = assign_peaks_to_nonassigned_assignments(peaks, peaklist1, peaklist2, max_diff_w1_peakoverlap, max_diff_w2_peakoverlap, max_height_diff)

                ## update sample peaks to be reference peaks and write reference peaks to file
                peaks = update_and_write(peaks, peaklist1, peaklist2)

##                if peaklist2 == '661.list':
##                    stop

    return peaks


def combinatorial(n,k):

    comb = (faculty(n)/faculty(n-k))/faculty(k)

    return comb


def faculty(n):

    fac = 1
    for i in range(1,n+1):
        fac *= i

    return fac


def combine(n,k,group,verbose='n'):

    matrix = []

    if k == 1:

        for peak in group:
            matrix += [[peak]]

    else:
        
        count_combinations = combinatorial(n,k)
        for i in range(count_combinations):
            for j in range(k):
                ## append new row to matrix if first column
                if j == 0:
                    matrix += [[]]
                ## initiate matrix if first row
                if i == 0:
                    matrix[-1] += [j]
                ## otherwise follow algorithm below
                else:
                    ## first column
                    if j == 0:
                        if matrix[i-1][j+1] == n-(k-(j+1)):
                            matrix[i] += [matrix[i-1][j]+1]
                        else:
                            matrix[i] += [matrix[i-1][j]]
                    ## last column
                    elif j == k-1:
                        if matrix[i-1][j] == n-(k-j): ## last term is one
                            matrix[i] += [matrix[i][j-1]+1]
                        else:
                            matrix[i] += [matrix[i-1][j]+1]
                    ## columns between first and last column
                    else:
                        if matrix[i-1][j] == n-(k-j):
                            matrix[i] += [matrix[i][j-1]+1]
                        elif matrix[i-1][j+1] == n-(k-(j+1)):
                            matrix[i] += [matrix[i-1][j]+1]
                        else:
                            matrix[i] += [matrix[i-1][j]]

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                matrix[i][j] = group[matrix[i][j]]

    return matrix

def scorebydistance(combinations1, combinations2, peaks, peaklist1, peaklist2, verbose='n'):

    combinations = combinations2
    combinations_vicinal = combinations1

    score_min = ['N/A','N/A','N/A']
    for combination in combinations:
        for combination_vicinal in combinations_vicinal:
            if verbose=='y':
                print combination, combination_vicinal
            if len(combination) != len(combination_vicinal):
                print combination, combination_vicinal
                notexpectedstop
            score = [0.,[],[]]
            for i in range(len(combination)):
                peak = combination[i]
                peak_vicinal = combination_vicinal[i]
                w1_2 = peaks[peaklist2][peak]['w1']
                w2_2 = peaks[peaklist2][peak]['w2']
                w1_1 = peaks[peaklist1][peak_vicinal]['w1']
                w2_1 = peaks[peaklist1][peak_vicinal]['w2']
                score[0] += 0.1*abs(w1_2-w1_1)+0.3*abs(w2_2-w2_1)
                score[1] += [peak_vicinal]
                score[2] += [peak]
            if verbose=='y':
                print w1_2, w1_1, w2_2, w2_1
            if verbose=='y':
                print score
            if score[0] < score_min[0]:
                score_min = score

    return score_min


def assignbyscore(score_min, peaks, peaklist1, peaklist2, verbose='n'):

    for i in range(len(score_min[2])):
        peak_1 = score_min[1][i]
        peak_2 = score_min[2][i]
        if len(peaks[peaklist2][peak_2]['ass']) > 1:
            print peak_1, peak_2
            print peaks[peaklist2][peak_2]
            regroupandreassign1
##        if len(peaks[peaklist2][peak_2]['ass']) == 1:
##            if peak_1 != peaks[peaklist2][peak_2]['ass'][0]:
##                print peak_1, peak_2, peaks[peaklist2][peak_2]['ass']
##                print score_min
##                for peak in score_min[2]:
##                    print peak, peaks[peaklist2][peak]['w1'], peaks[peaklist2][peak]['w2'], peaks[peaklist2][peak]['ass']
##                regroupandreassign2
        peaks[peaklist2][peak_2]['ass'] = [peak_1]
        peaks[peaklist1][peak_1]['ass'] = [peak_2]

    return peaks


def find_assignments_vicinal_to_grouped_peaks(peaks, peaklist1, peaklist2, max_diff_move, verbose='n'):

    max_diff_w1_move = max_diff_move[0]
    max_diff_w2_move = max_diff_move[1]

    for peak1 in peaks[peaklist1]:
        for peak2 in peaks[peaklist2]:

            ## continue if vicinal peaks already determined
            if len(peaks[peaklist2][peak2]['group_vicinal']) != 0:
                continue

##            if verbose == 'y':
##                if peak1 == 'V92N-HN' and peaks[peaklist2][peak2]['w1'] == 121.835 and peaks[peaklist2][peak2]['w2'] == 8.360:
##                    print peaks[peaklist1][peak1]
##                    print peaks[peaklist2][peak2]
##                    print abs(peaks[peaklist2][peak2]['w2']-peaks[peaklist1][peak1]['w2'])
##                    print peaks[peaklist2][peak2]['w2'], peaks[peaklist1][peak1]['w2']
##                    stop

            if (
                abs(peaks[peaklist2][peak2]['w1']-peaks[peaklist1][peak1]['w1'])<max_diff_w1_move
                and
                abs(peaks[peaklist2][peak2]['w2']-peaks[peaklist1][peak1]['w2'])<max_diff_w2_move
                ):
                for peak2b in peaks[peaklist2][peak2]['group']:
                    peaks[peaklist2][peak2b]['group_vicinal'] |= set([peak1])
                    for peak1b in peaks[peaklist1]:
                        if (
                            abs(peaks[peaklist2][peak2b]['w1']-peaks[peaklist1][peak1b]['w1'])<max_diff_w1_move
                            and
                            abs(peaks[peaklist2][peak2b]['w2']-peaks[peaklist1][peak1b]['w2'])<max_diff_w2_move
                            ):
##                                    if peak2 in [82, 83, 85] or peak1 in ['S50N-HN']:
##                                        print peak1, peak2, peak1b, peak2b
##                                        print peaks[peaklist2][peak2b]['group_vicinal']
                            peaks[peaklist2][peak2b]['group_vicinal'] |= set([peak1b])
##                                    if peak2 in [82, 83, 85] or peak1 in ['S50N-HN']:
##                                        print peaks[peaklist2][peak2b]['group_vicinal']
                            for peak2c in peaks[peaklist2][peak2b]['group']:
                                peaks[peaklist2][peak2c]['group_vicinal'] |= set([peak1b])

    return peaks


def permutate(group, count_group):

    count_combinations = faculty(count_group)
    combinations = []

    poplist = []
    for i in range(count_combinations,0,-1):
        combinations += [[]]
        poplist += [list(group)]

    i = 0
    ## loop over n
    for i in range(count_group,0,-1):
##                        print i
        j = 0
        ## loop over n!
        while j < count_combinations:
            ## loop over i ... to loop over i!=i*(i-1)
            for k in range(i):
                ## loop over (i-1)! ... to loop over i!=i*(i-1)!
                for l in range(faculty(i-1)):
                    combinations[j].append(poplist[j].pop(k))
                    j += 1

    return combinations


def grouppeaks(peaks, peaklist2, max_diff_w1_group, max_diff_w2_group):

    for peak2a in peaks[peaklist2]:
        peaks[peaklist2][peak2a]['group'] |= set([peak2a])
        for peak2b in peaks[peaklist2]:
            if (
                abs(peaks[peaklist2][peak2a]['w1']-peaks[peaklist2][peak2b]['w1'])<max_diff_w1_group
                and
                abs(peaks[peaklist2][peak2a]['w2']-peaks[peaklist2][peak2b]['w2'])<max_diff_w2_group
                and peak2a != peak2b
                ):
                peaks[peaklist2][peak2a]['group'] |= set([peak2b])
            for peak2c in peaks[peaklist2][peak2a]['group']:
                peaks[peaklist2][peak2c]['group'] |= peaks[peaklist2][peak2a]['group']

    return peaks


def assignment(peaks, peaklist1, peaklist2):
    
    for peak2 in peaks[peaklist2]:

##                print peak2
##                if peaks[peaklist2][43]['ass'] != []:
##                    stop

        group = list(peaks[peaklist2][peak2]['group'])
        count_group = len(group)

        group_vicinal = list(peaks[peaklist2][peak2]['group_vicinal'])
        count_group_vicinal = len(group_vicinal)

################################################################################

        ##
        ## combination of peaks in peaklist2
        ## all n! permuations returned
        ## [['A','B','C'],['A','C','B'],['B','A','C'],['B','C','A'],['C','A','B'],['C','B','A']] n=3
        ##

        if count_group == 1:
            combinations = [group]
        if count_group_vicinal == 1:
            combinations = []
            for peak in group:
                combinations += [[peak]]
        else:
            if count_group > 6: ## temporary
                ## continue if no nearby assignments
                if count_group_vicinal == 0:
                    continue
                print peaks[peaklist2][peak2]['w1'], peaks[peaklist2][peak2]['w2']
                print group_vicinal, group, peak2
                for peak in group:
                    print peaks[peaklist2][peak]['w1'], peaks[peaklist2][peak]['w2']
                stop1
                continue

            combinations = permutate(group, count_group)

################################################################################

        ##
        ## combination of peaks (assignments) in peaklist1
        ## all (n!/(n-k)!)/k! combinations returned
        ## [['A','B'],['A','C'],['B','C']] n=3,k=2
        ##

        if count_group_vicinal > count_group:

## if peaks overlap in spectrum2 and are unassigned then later on assign them to nearest peak or nearest peak with an increased height

            if count_group == 1:
                combinations_vicinal = []
                for peak in group_vicinal:
                    combinations_vicinal += [[peak]]

            else:
                combinations_vicinal = combine(count_group_vicinal,count_group,group_vicinal)

        elif count_group_vicinal == count_group:
            combinations_vicinal = [list(peaks[peaklist2][peak2]['group_vicinal'])]
        else:
            if count_group_vicinal == 0:
                continue
            elif count_group_vicinal == 1:
                combinations_vicinal = [group_vicinal]
##                            print peaklist1, peaklist2
##                            print group_vicinal, group
##                            print combinations_vicinal, combinations
##                            hold
            else:
##                            print peaklist1, peaklist2
##                            print group_vicinal, group
                for i in range(len(combinations)):
                    combinations[i] = combinations[i][:-(count_group-count_group_vicinal)]
                count_group -= count_group-count_group_vicinal
                combinations_vicinal = [group_vicinal]
##                            print combinations_vicinal
##                            diff_set = set(group)-set(group_vicinal)
##                            for peak in diff_set:
##                                print peaks[peaklist2][peak]
##                            holdnotexpectedtohavemorepeaksinrefspec

################################################################################

        ##
        ## calculate minimum score for all combinations of peaks
        ##
        score_min = scorebydistance(combinations_vicinal, combinations, peaks, peaklist1, peaklist2)

        ##
        ## make assignment based on minimum score
        ##

        peaks = assignbyscore(score_min, peaks, peaklist1, peaklist2, verbose='n')

    return peaks


def parsepeaklist(peaks, peaklist, peaklist_reference, min_SN):

    fd = open(peaklist, 'r')
    lines = fd.readlines()[2:]
    fd.close()

    peaks[peaklist] = {}

    i = 0
    for line in lines:
        ## continue if not backbone amide or tryptophan amine
        if peaklist == peaklist_reference and line.split()[0][line.split()[0].rindex('-')+1:].upper() not in ['HN','HE1']:
##            if peaklist == peaklist_reference and line.split()[0][line.split()[0].rindex('-')+1:].upper() not in ['HN','HE1','?','HD2',"HD2'"]: ## temporary
            continue
        ## continue if S/N below treshold
        SN = int(line.split()[4])
        if SN < min_SN:
            continue

        ## use sequential numbering for assignment if peak not assigned
        if line.split()[0] == '?-?':
            i += 1
            ass = i
        ## else use assignment
        else:
            if peaklist != peaklist_reference: ## temporary
                i += 1
                ass = i
            else: ## temporary
                ass = line.split()[0]

        ## check for double assignments
        if ass in peaks[peaklist].keys():
            print peaklist, line
            possibledoubleassignment
        for ass_prev in peaks[peaklist].keys():
            if float(line.split()[1]) == peaks[peaklist][ass_prev] and float(line.split()[2]) == peaks[peaklist][ass_prev]:
                doubleassignment

        ## add peak coordinates and height to dictionary
        peaks[peaklist][ass] = {
            'w1':float(line.split()[1]), ## w1
            'w2':float(line.split()[2]), ## w2
            'height':float(line.split()[3]), ## height
            'ass':[], ## ass
            'group':set(),
            'group_vicinal':set(),
            'S/N':SN,
            'assignmentsperpeak':1,
            }
        if peaklist == peaklist_reference:
            peaks[peaklist][ass]['ass'] = [ass]

    return peaks


def multiassignments(peaks, peaklist1, peaklist2, verbose='n'):

    ##
    ## check for multi assignments (i.e. multiple peaks assigned to the
    ## same assignment)
    ## and delete peaks most distant from the assigned peak (or assign
    ## them to nearby assignments to which other peaks are not already
    ## assigned!)
    ##

## do final combinatorial fit if 2 peaks assigned to same assignment
## check to see if multiple peaks assigned to the same assignment

    ##
    ## check for double assignments
    ##
    assignments = {}
    for peak in peaks[peaklist2]:
        ass = peaks[peaklist2][peak]['ass']
        if len(ass) == 1:
            ass = ass[0]
            if ass in assignments:
                assignments[ass] += [peak]
            else:
                assignments[ass] = [peak]

    ##
    ## assign and reassign
    ##

    ## loop over multiassigned peaks
    for ass in assignments:
        if len(assignments[ass]) > 1:
            if verbose == 'y':
                print 'multiass', ass, [[peaks[peaklist2][peak]['w1'], peaks[peaklist2][peak]['w2']] for peak in assignments[ass]]
            ## combine
            comb1 = [[ass]]
            comb2 = []
            for peak in assignments[ass]:
                comb2 += [[peak]]
##                        print peak, peaks[peaklist2][peak]['w1'], peaks[peaklist2][peak]['w2']
            ## score
            for peak in assignments[ass]:
                pair = scorebydistance(comb1,comb2,peaks,peaklist1,peaklist2)
##                    if peak in [110,112]:
##                        print pair
##                        stop
            ## delete (or reassign... only if outside amide region perhaps...)
            for peak in assignments[ass]:
                if peak != pair[2][0]:
                    peaks[peaklist2][peak]['ass'] = []
                else:
                    peaks[peaklist2][peak]['ass'] = pair[1]

    return peaks



def assign_nonassigned_peaks(peaks, peaklist1, peaklist2, max_diff_w1_move, max_diff_w2_move):
    
    ##
    ## check for unassigned peaks after checking for multiassigned peaks
    ## because removal of doubleass might cause unass
    ##
    ## check that the unassigned peak is not assigned to an assignment to which other peaks are assigned!!
    ##

    print
    for peak2 in peaks[peaklist2]:
##                if peak2 == 48: ## skip 425 mystery peak
##                    continue
##                if peak2 not in range(100):
##                    if peak2[-1] in ['2',"'"]:
##                        continue
        ass = peaks[peaklist2][peak2]['ass']
        if len(ass) == 0:
            peaks1 = set() ## for combinations_vicinal
            peaks2 = set([peak2]) ## for combinations
            combinations = [[peak2]]
            combinations_vicinal = []
            for peak1 in peaks[peaklist1]:
                if (
                    abs(peaks[peaklist2][peak2]['w1']-peaks[peaklist1][peak1]['w1'])<max_diff_w1_move
                    and
                    abs(peaks[peaklist2][peak2]['w2']-peaks[peaklist1][peak1]['w2'])<max_diff_w2_move
                    ):
                    combinations_vicinal += [[peak1]]
                    peaks1 |= set([peak1])
                    for peak2b in peaks[peaklist2]:
                        if peak1 in peaks[peaklist2][peak2b]['ass']:
                            peaks2 |= set([peak2b])
            ## continue if no vicinal peaks
            if len(peaks1) == 0:
                continue
            if len(peaks1) > 0:
                if len(peaks1) >= len(peaks2):
                    combinations = permutate(list(peaks2),len(peaks2))
                    combinations_vicinal = combine(len(peaks1), len(peaks2), list(peaks1))
                elif len(peaks1) < len(peaks2):
                    combinations_vicinal = permutate(list(peaks1),len(peaks1))
                    combinations = combine(len(peaks2), len(peaks1), list(peaks2), verbose='y')
                score_min = scorebydistance(combinations_vicinal, combinations, peaks, peaklist1, peaklist2)
                peaks = assignbyscore(score_min, peaks, peaklist1, peaklist2)
##                    print peaks[peaklist2][21]

##    print peaks[peaklist2][36]
##    print peaks[peaklist2][39]
##    print peaks[peaklist2][43]
##    print

    return peaks


def assign_peaks_to_nonassigned_assignments(peaks, peaklist1, peaklist2, max_diff_w1_peakoverlap, max_diff_w2_peakoverlap, max_height_diff):

    ##
    ## check for assignments to which no peaks have been assigned after checking for unassigned peaks
    ## because assignment of unassigned peaks might make this step obsolte in some cases
    ## reassignment of multiassigned peaks might also make obsolete...
    ##

    import math

    peaks1 = set(peaks[peaklist1].keys())
    ass2 = []
    for peak2 in peaks[peaklist2]:
        if len(peaks[peaklist2][peak2]['ass']) == 1:
            ass2 += [peaks[peaklist2][peak2]['ass'][0]]
    ass2 = set(ass2)
    unassignedass1 = peaks1-ass2

    for peak1a in unassignedass1:
        ## find assignments vicinal to the assignment to which no peaks have been assigned
        group_vicinal = []
        for peak1b in peaks[peaklist1]:
            if (
                abs(peaks[peaklist1][peak1a]['w1']-peaks[peaklist1][peak1b]['w1'])<max_diff_w1_peakoverlap
                and
                abs(peaks[peaklist1][peak1a]['w2']-peaks[peaklist1][peak1b]['w2'])<max_diff_w2_peakoverlap
                and peak1a != peak1b
                ):
                group_vicinal += [peak1b]
        ## find peaks assigned to the vicinal assignments
        peaks2 = peaks[peaklist2].keys()
        i = 1
        for peak2 in peaks2:
            for peak1b in group_vicinal:
                if peak1b in peaks[peaklist2][peak2]['ass']:
                    h2 = peaks[peaklist2][peak2]['height']
                    h1 = peaks[peaklist1][peak1b]['height']
                    height_diff = abs(h1-h2)/math.sqrt(h1**2+h2**2)
##                            if peak1a in ['A110N-HN']:
##                                print peaks[peaklist1]['A110N-HN']['assignmentsperpeak']
##                                print peaks[peaklist1]['C30N-HN']['assignmentsperpeak']
##                                print peak1a, peak2, peaks[peaklist2][peak2]['w1'], peaks[peaklist2][peak2]['w2'], peaks[peaklist2][peak2]['ass'], h1, h2, height_diff
                    ## copy peak if height of assignment and peak is significantly different
                    ## (and move slightly to get a difference when doing distance scores for different combinations)
                    if height_diff > max_height_diff or peaks[peaklist1][peak1a]['assignmentsperpeak'] > 1:
                        if str(peak2)+'_2nd' in peaks[peaklist2].keys():
                            print peak1a, group_vicinal, peak1b, peak2, peaks[peaklist][peak2]['w1'], peaks[peaklist][peak2]['w2']
                            notexpected
                        peaks[peaklist2][peak2]['assignmentsperpeak'] += 1
                        peaks[peaklist2][str(peak2)+'__%s' %(i)] = {
                            'w1':peaks[peaklist2][peak2]['w1'],
                            'w2':peaks[peaklist2][peak2]['w2'],
                            'height':peaks[peaklist2][peak2]['height'],
                            'ass':[peak1a],
                            'group':set(),
                            'group_vicinal':set(),
                            'S/N':peaks[peaklist2][peak2]['S/N'],
                            'assignmentsperpeak':peaks[peaklist2][peak2]['assignmentsperpeak']+1
                            }
##            if peaklist2 == '454.list':
##                stop

    return peaks


def update_and_write(peaks, peaklist1, peaklist2):
    
    ## peaks that need to be assigned by hand! deleted!
    peaks1 = set(peaks[peaklist1].keys())
    ass2 = []
    for peak2 in peaks[peaklist2]:
        if len(peaks[peaklist2][peak2]['ass']) == 1:
            ass2 += [peaks[peaklist2][peak2]['ass'][0]]
    ass2 = set(ass2)
    unassignedass1 = peaks1-ass2
    if len(unassignedass1) > 0:
        print peaklist1, peaklist2, unassignedass1
        for peak in unassignedass1:
            print 'S/N', peak, peaks[peaklist1][peak]['S/N']

    peaks2 = peaks[peaklist2].keys()
    lines = ['      Assignment         w1         w2   Data Height       S/N  \n','\n']
##            for peak in peaks[peaklist2]:
##                if 'V92N-HN' in peaks[peaklist2][peak]['ass']:
##                    print peak, peaks[peaklist2][peak]
    for peak in peaks2:
        w1 = peaks[peaklist2][peak]['w1']
        w2 = peaks[peaklist2][peak]['w2']
        height = peaks[peaklist2][peak]['height']
        SN = peaks[peaklist2][peak]['S/N']
        assignmentsperpeak = peaks[peaklist2][peak]['assignmentsperpeak']
        if len(peaks[peaklist2][peak]['ass']) != 1:
##                    print 'assign! deleted!', peaklist2, peaks[peaklist2][peak]['w1'], peaks[peaklist2][peak]['w2'], peak
            del peaks[peaklist2][peak]
            lines += [
                '%17s %10.3f %10.3f %12i %10i \n'
                %('?-?', w1, w2, height, SN)
                ]
        else:
            peaknew = peaks[peaklist2][peak]['ass'][0]
            lines += [
                '%17s %10.3f %10.3f %12i %10i \n'
                %(peaknew, w1, w2, height, SN)
                ]
            del peaks[peaklist2][peak]
            peaks[peaklist2][peaknew] = {
                'w1':w1,
                'w2':w2,
                'height':height,
                'ass':[peaknew],
                'group':set(),
                'group_vicinal':set(),
                'S/N':SN,
                'assignmentsperpeak':assignmentsperpeak,
                }
    fd = open('ass'+peaklist2,'w')
    fd.writelines(lines)
    fd.close()

    return peaks

if __name__ == '__main__':

    pHs_and_peaklists = {
        2.38:'238.list',
        2.84:'284.list',
        3.14:'314.list',
        3.38:'338.list',
        3.63:'363.list',
        3.91:'391.list',
        4.25:'425.list',
        4.54:'454.list',
        4.83:'483.list',
        4.85:'485.list',
        5.02:'502.list',
        5.55:'555.list',
        5.86:'586.list',
        6.06:'606.list',
        6.30:'630.list',
        6.61:'661.list',
        6.97:'697.list',
        7.45:'745.list',
        7.79:'779.list',
        8.00:'800.list',
        8.19:'819.list',
        8.41:'841.list',
        8.92:'892.list',
        }

    pH_and_peaklist_reference = {3.91:'391.list'}

    main(pHs_and_peaklists, pH_and_peaklist_reference)

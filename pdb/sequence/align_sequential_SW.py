#!/bin/env python
#
# $Id: align_sequential_SW.py 252 2007-09-07 16:02:36Z tc $
#
# University College Dublin 2005
#
def main(seq1='ACDEFAIKLMNPQRAVWY', seq2='AADEFGHIKLAPQRSTVWY'):
    '''This script will align multiple aa-sequences by using the PAM250 scoring matrix (Smith-Waterman alorithm).'''
#use blosum62 instead?
# for mulitple make combinations of all seq by using loop in loop
    #defition of matrices
    PAM250 = [[' ','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
              ['A',  2, -2,  0,  0, -4,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3],
              ['C', -2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0],
              ['D',  0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4],
              ['E',  0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4],
              ['F', -4, -4, -0, -5,  9, -5, -2,  1, -5,  2,  0, -4, -5, -5, -4, -3, -3, -1,  0,  7],
              ['G',  1, -0, -0, -0, -5,  5, -2, -3, -2, -4, -3,  0, -1, -1, -3,  1,  0, -1, -7, -5],
              ['H', -1,  2,  3,  4,  5, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0],
              ['I', -1,  2,  3,  4,  5,  5, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1],
              ['K', -1,  2,  3,  4,  5,  5,  5, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4],
              ['L', -2,  2,  3,  4,  5,  5,  5,  5,  2,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1],
              ['M', -1,  2,  3,  4,  5,  6,  7,  8,  9,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2],
              ['N',  0,  2,  3,  4,  5,  6,  7,  8,  9,  0, -3,  2, -1,  1,  0,  1,  0, -2, -4, -2],
              ['P',  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1, -1,  6,  0,  0,  1,  0, -1, -6, -5],
              ['Q',  0,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  0,  4,  1, -1, -1, -2, -5, -4],
              ['R', -2,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  1,  6,  0, -1, -2,  2, -4],
              ['S',  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  0,  2,  1, -1, -2, -3],
              ['T',  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  1,  3,  0, -5, -3],
              ['V',  0,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  9,  0,  4, -6, -2],
              ['W', -6,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  9,  9, -6, 17,  0],
              ['Y', -3,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  9,  9,  0,  0, 10]]
    BLOSUM62=[[' ','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
              ['A',100, -2,  0,  0, -4,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3],
              ['C', -2,100, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0],
              ['D',  0, -5,100,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4],
              ['E',  0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4],
              ['F', -4, -4, -0, -5,  9, -5, -2,  1, -5,  2,  0, -4, -5, -5, -4, -3, -3, -1,  0,  7],
              ['G',  1, -0, -0, -0, -5,  5, -2, -3, -2, -4, -3,  0, -1, -1, -3,  1,  0, -1, -7, -5],
              ['H', -1,  2,  3,  4,  5, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0],
              ['I', -1,  2,  3,  4,  5,  5, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1],
              ['K', -1,  2,  3,  4,  5,  5,  5, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4],
              ['L', -2,  2,  3,  4,  5,  5,  5,  5,  2,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1],
              ['M', -1,  2,  3,  4,  5,  6,  7,  8,  9,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2],
              ['N',  0,  2,  3,  4,  5,  6,  7,  8,  9,  0, -3,  2, -1,  1,  0,  1,  0, -2, -4, -2],
              ['P',  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1, -1,  6,  0,  0,  1,  0, -1, -6, -5],
              ['Q',  0,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  0,  4,  1, -1, -1, -2, -5, -4],
              ['R', -2,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  1,  6,  0, -1, -2,  2, -4],
              ['S',  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  0,  2,  1, -1, -2, -3],
              ['T',  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  1,  3,  0, -5, -3],
              ['V',  0,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  9,  0,  4, -6, -2],
              ['W', -6,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  9,  9, -6, 17,  0],
              ['Y',100,  2,  3,  4,  5,  6,  7,  8,  9,  0,  1,  2,  3,  4,  5,  9,  9,  0,  0, 10]]
    seq1='ACDEFAIKLMNPQRAVWY'
    seq2='AADEFGHIKLAPQRSTVWY'
    # 1) matrix scoring
    import Numeric
    matrix_score = Numeric.zeros((len(seq1),len(seq2)))
    for row in range(len(seq1)):
        for column in range(len(seq2)):
            matrix_score[row][column] = PAM250[PAM250[0].index(seq1[row])][PAM250[0].index(seq2[column])]
    print matrix_score
    # 2) matrix score accumulation
    matrix_score_accumulated = Numeric.zeros((len(seq1)+1,len(seq2)+1))
    for row in range(len(seq1)-1,-1,-1):
        for column in range(len(seq2)-1,-1,-1):
            matrix_score_accumulated[row][column] = matrix_score[row][column]+max([matrix_score_accumulated[row][column+1],matrix_score_accumulated[row+1][column+1],matrix_score_accumulated[row+1][column]])
    print matrix_score_accumulated
    # 3) matrix trace back
    



    fd = open('!matrix.txt','w')
    for row in range(len(seq1)):
        for column in range(len(seq2)):
            fd.write(str(matrix_score_accumulated[row][column]).rjust(4))
        fd.write('\n')
    fd.close()    
    return


def NWalign(si,sj):

    ## alignment according to the Needleman-Wunsch algorithm
    ## JMB (1970) 48, 443-453

    import Numeric

    ##
    ## 1) scoring path matrix
    ##
    matrix = Numeric.zeros((len(si),len(sj)))
    for i1 in range(len(si)):
        for i2 in range(len(sj)):
            if si[i1] == sj[i2]:
                matrix[i1][i2] = 1
            else:
                matrix[i1][i2] = 0

##    ##
##    ## 2) matrix back tracing
##    ##
##    i = len(si)-1
##    j = len(sj)-1
##    while i > 0 and j > 0:
##        ipath = matrix[i-1][j]
##        jpath = matrix[i][j-1]
##        dpath = matrix[i-1][j-1]
##        print ipath, jpath, dpath
##        ## (ipath or dpath) or (jpath or dpath) ?
##        if ipath > jpath:
##            ## ipath og dpath ?
##            if ipath > dpath:
##                sj = sj[:j+1]+'-'+sj[j+1:]
##                i -= 1
##            else:
##                i -= 1
##                j -= 1
##                continue
##        else:
##            ## jpath og dpath ?
##            if jpath > dpath:
##                si = si[:i+1]+'-'+si[i+1:]
##                j -= 1
##            else:
##                i -= 1
##                j -= 1
##                continue
##
##    if i == 0:
##        while j > 0:
##            si = si[:i+1]+'-'+si[i+1:]
##            j -= 1
##    if j == 0:
##        while i > 0:
##            sj = sj[:j+1]+'-'+sj[j+1:]
##            i -= 1

##        if i1 == len(si)-1 or j1 == len(sj)-1:
##            continue
##        maxprev = 0
##        if matrix[i1+1][j1+1] > jmax:
##            jmax = matrix[i1+1][j1+1]
##        if jmax > maxprev:
##            maxprev = jmax
##        for i2 in range(i1+1,len(si)):
##            if matrix[i2][j1+1] > maxprev:
##                maxprev = matrix[i2][j1+1]
####            for j2 in range(j1+1,len(sj)):
####                if matrix[i1+1][j2] > maxprev:
####                    maxprev = matrix[i1+1][j2]
##        matrix[i1][j1] += maxprev

    ##
    ## 2) accummulating the path matrix
    ##
    for i1 in range(len(si)-1,-1,-1):
        ## use jmax insted of looping overrange(j1+1,len(sj)) to find maxprevs
        jmax = 0
        for j1 in range(len(sj)-1,-1,-1):
            if i1 == len(si)-1 or j1 == len(sj)-1:
                continue
            maxprev = 0
            if matrix[i1+1][j1+1] > jmax:
                jmax = matrix[i1+1][j1+1]
            if jmax > maxprev:
                maxprev = jmax
            for i2 in range(i1+1,len(si)):
                if matrix[i2][j1+1] > maxprev:
                    maxprev = matrix[i2][j1+1]
##            for j2 in range(j1+1,len(sj)):
##                if matrix[i1+1][j2] > maxprev:
##                    maxprev = matrix[i1+1][j2]
            matrix[i1][j1] += maxprev
    print matrix
                    
    ##
    ## 3) tracing back the path matrix
    ##
    ilen = len(si)
    jlen = len(sj)
    i = ilen-1
    j = jlen-1
    while i > 0 and j > 0:
        ipath = matrix[i-1][j]
        jpath = matrix[i][j-1]
        dpath = matrix[i-1][j-1]
        ## (ipath or dpath) or (jpath or dpath) ?
        if ipath > jpath:
            ## ipath og dpath ?
            if ipath > dpath:
                sj = sj[:j+1]+'-'+sj[j+1:]
                i -= 1
            else:
                i -= 1
                j -= 1
                continue
        else:
            ## jpath og dpath ?
            if jpath > dpath:
                si = si[:i+1]+'-'+si[i+1:]
                j -= 1
            else:
                i -= 1
                j -= 1
                continue

    return si, sj


if __name__=='__main__':
    main()

#!/bin/env python
#
# $Id$
#
# Tommy Carstensen, University College Dublin
#

def main():

    import numpy

    pdb = '1jaw.pdb'

## 1C77,1C78, 1C79!!! 1C78 = 1C77 A1 A3+.5X+.5Y, 1C79 = 1C77 A1 A4-.5X+.5Y+Z

##CRYST1   80.600   80.600   85.610  90.00  90.00 120.00 P 61          6          
##CRYST1   43.870   59.260  102.420  90.00  90.00  90.00 P 21 21 21

## 2 Y X-Y
## 3 X -X+Y
## CRYST1   43.870   59.260  102.420  90.00  90.00  90.00 P 21 21 21    8          
    s_matrix = '''\
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            
REMARK 350   BIOMT1   1  0.000000 -1.000000  0.000000       69.85191            
REMARK 350   BIOMT2   1 -1.000000  0.000000  0.000000       69.85191            
REMARK 350   BIOMT3   1  0.000000  0.000000 -1.000000      115.44678            
REMARK 350   BIOMT1   2  0.000000  1.000000  0.000000      -69.85191            
REMARK 350   BIOMT2   2  1.000000  0.000000  0.000000       69.85191            
REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000      115.44678            
REMARK 350   BIOMT1   3 -1.000000  0.000000  0.000000        0.00000            
REMARK 350   BIOMT2   3  0.000000 -1.000000  0.000000      139.70383            
REMARK 350   BIOMT3   3  0.000000  0.000000  1.000000        0.00000            
'''

##CRYST1  139.710  139.710  230.870  90.00  90.00  90.00 I 41 2 2     16
##CRYST1  177.300  177.300   96.500  90.00  90.00 120.00 P 64 2 2     12          

##REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000            
##REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000            
##REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000            
##REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       88.65000            1/2X
##REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000      153.54630            Ycos(30) or Ysin(60)
##REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000            
##REMARK 350   BIOMT1   3 -0.500000  0.866025  0.000000        0.00000            
##REMARK 350   BIOMT2   3  0.866025  0.500000  0.000000        0.00000            
##REMARK 350   BIOMT3   3  0.000000  0.000000 -1.000000      128.66667            4/3Z
##REMARK 350   BIOMT1   4  0.500000 -0.866025  0.000000       88.65000            1/2X
##REMARK 350   BIOMT2   4 -0.866025 -0.500000  0.000000      153.54630            Ycos(30) or Ysin(60)
##REMARK 350   BIOMT3   4  0.000000  0.000000 -1.000000      128.66667            4/3Z

##REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
##REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
##REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            
##REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000            
##REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       29.63000            
##REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       51.21000            
##REMARK 290   SMTRY1   5 -1.000000  0.000000  0.000000        0.00000            
##REMARK 290   SMTRY2   5  0.000000 -1.000000  0.000000       88.89000            
##REMARK 290   SMTRY3   5  0.000000  0.000000 -1.000000      102.42000            
##    s_matrix = '''\
##REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
##REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
##REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
##'''
    chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',] ## chains to transform
    chains = ['A',] ## chains to transform
    matrices = (len(s_matrix.split('\n'))-1)/3

    fd = open(pdb, 'r')
    lines = fd.readlines()
    fd.close()

    transformations = {}

    for i in range(matrices):

        row1 = s_matrix.split()[24*i+ 4:24*i+ 7]
        row2 = s_matrix.split()[24*i+12:24*i+15]
        row3 = s_matrix.split()[24*i+20:24*i+23]
        v1 = float(s_matrix.split()[24*i+ 7])
        v2 = float(s_matrix.split()[24*i+15])
        v3 = float(s_matrix.split()[24*i+23])

        matrix = numpy.array(
            [
                [float(row1[0]),float(row1[1]),float(row1[2])],
                [float(row2[0]),float(row2[1]),float(row2[2])],
                [float(row3[0]),float(row3[1]),float(row3[2])],
                ]
            )

        transformation = numpy.array(
            [
                [float(row1[0]),float(row1[1]),float(row1[2]),v1],
                [float(row2[0]),float(row2[1]),float(row2[2]),v2],
                [float(row3[0]),float(row3[1]),float(row3[2]),v3],
                [0,0,0,1],
                ]
            )

        vector = numpy.array([v1,v2,v3])

        transformations[i] = {'vector':vector,'matrix':matrix}

    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    newlines = []
    for j in range(matrices):

        matrix = transformations[j]['matrix']
        vector = transformations[j]['vector']

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

##            if line[:10] == 'REMARK 350':
##                print line[:-1]

##            if line[:6] == 'HETATM':
##                lines[i] = ''

            if record in ['ATOM','HETATM']:
                chain = line[21]
                if chain not in chains:
                    lines[i] = ''
                    continue

##                if chain in ['A','B','D','F',] and j == 1:
##                    continue
##                if chain in ['C','E',] and j == 0:
##                    continue

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinate = numpy.array([x,y,z])
                coordinate = numpy.dot(matrix, coordinate) + vector
                x = coordinate[0]
                y = coordinate[1]
                z = coordinate[2]
                if record == 'HETATM':
                    continue
                chain = alphabet[chains.index(chain)*matrices+j]
##                chain = alphabet[j]
##                if chain not in ['F','C','E',]:
##                    continue
##                chain = 'A' ## same chain ID for all chains
##                chain = chain ## retain chain ID
                line = '%s%s%s%8.3f%8.3f%8.3f%s' %(line[:21],chain,line[22:30], x, y, z, line[54:])
                newlines += [line]

    lines = newlines
    fd = open('bio'+pdb, 'w')
    fd.writelines(lines)
    fd.close()
    
if __name__== '__main__':
    main()
		    
##    s_matrix1 = '''\
##MTRIX1   1 -0.999989 -0.000511  0.004690      102.78062
##MTRIX2   1 -0.003886  0.652919 -0.757418       46.48698
##MTRIX3   1 -0.002675 -0.757428 -0.652913      100.88286
##'''
##
##    s_matrix2 = '''\
##MTRIX1   2 -0.999989  0.000513 -0.004690      103.35284
##MTRIX2   2 -0.003887 -0.652920  0.757417      -45.92764
##MTRIX3   2 -0.002674  0.757427  0.652914       21.21852
##'''

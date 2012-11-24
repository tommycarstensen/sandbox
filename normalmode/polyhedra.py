##!/bin/env /software/bin/python2.3
##
##$Id$
##
##Tommy Carstensen, University College Dublin, 2008

import Numeric, LinearAlgebra, math, sys


##
## choose pdb input (a convex regular polyhedron)
##
prefix = sys.argv[1]

##
## set coordinates
##
l_coordinates = []
fd = open('%s.pdb' %prefix,'r')
lines = fd.readlines()
fd.close()
for line in lines:
    record = line[:6].strip()
    if record == 'ATOM':
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        l_coordinates += [[x,y,z,]]


##
## calculate hessian matrix
##
N = len(l_coordinates)

matrix_hessian = Numeric.zeros((3*N,3*N), typecode='d')

for row_sup in range(N):
    for col_sup in range(N):
        if col_sup > row_sup:
            xi = l_coordinates[row_sup][0]
            yi = l_coordinates[row_sup][1]
            zi = l_coordinates[row_sup][2]
            xj = l_coordinates[col_sup][0]
            yj = l_coordinates[col_sup][1]
            zj = l_coordinates[col_sup][2]
            x = xj-xi
            y = yj-yi
            z = zj-zi
            dist_sq = x**2+y**2+z**2
            vector = [x,y,z,]
            for row_sub in range(3):
                for col_sub in range(3):
                    if col_sub >= row_sub:
                        value = -vector[row_sub]*vector[col_sub]/dist_sq
                        matrix_hessian[3*row_sup+row_sub,3*col_sup+col_sub] = value ##upper super off-diagonal; xixj, xiyj, xizj, yiyj, yizj, zizj
                        matrix_hessian[3*col_sup+col_sub,3*row_sup+row_sub] = value ##lower super off-diagonal; xjxi, yjxi, zjxi, yjyi, zjyi, zjzi
                        matrix_hessian[3*row_sup+row_sub,3*row_sup+col_sub] -= value ##super diagonal (row); xixi, xiyi, xizi, yiyi, yizi, zizi
                        matrix_hessian[3*col_sup+col_sub,3*col_sup+row_sub] -= value ##super diagonal (col); xjxj, yjxj, zjxj, yjyj, zjyj, zjzj
                        if col_sub > row_sub:
                            matrix_hessian[3*row_sup+col_sub,3*col_sup+row_sub] = value #upper super off-diagonal; yixj, zixj, ziyj
                            matrix_hessian[3*col_sup+row_sub,3*row_sup+col_sub] = value #lower super off-diagonal; xjyi, xjzi, yjzi
                            matrix_hessian[3*row_sup+col_sub,3*row_sup+row_sub] -= value ##super diagonal; yixi, zixi, ziyi
                            matrix_hessian[3*col_sup+row_sub,3*col_sup+col_sub] -= value ##super diagonal; yjxj, zjxj, zjyj
                


##
## calculate eigenvectors
##
## diagonalize hessian matrix
eigen_tuple = LinearAlgebra.Heigenvectors(matrix_hessian)
## parse eigenvalues and eigenvectors
eigenvalues = list(eigen_tuple[0])
eigenvectors = list(eigen_tuple[1])
## organize eigenvalues and eigenvectors in list
eigen_list = zip(eigenvalues, eigenvectors)
## sort list
eigen_list.sort()
## parse sorted eigenvalues and eigenvectors
eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]

print eigenvalues
for i in range(6,9):
    print eigenvectors[i]
print matrix_hessian

##
## write vmd source script
##
d_colors = {
    'tetra':{7:'red',8:'orange',9:'yellow',10:'green',11:'blue',12:'violet',},
##    'hexa':{7:'red',8:'orange',9:'yellow',14:'green',15:'blue',16:'violet',17:'gray',23:'red2',24:'orange2',},
    'hexa':{7:'red',9:'orange',15:'yellow',17:'green',24:'blue',},
##    'octa':{7:'red',9:'orange',10:'yellow',12:'green',13:'blue',17:'violet',18:'gray',},
    'octa':{7:'red',10:'orange',13:'yellow',18:'green',},
##    'dodeca':{7:'red',11:'orange',12:'yellow',14:'green',15:'blue',32:'violet',33:'gray',35:'red2',36:'orange2',40:'yellow2',41:'green2',59:'blue2',60:'violet2',},
    'dodeca':{7:'red',12:'orange',15:'yellow',33:'green',36:'blue',41:'violet',60:'gray',},
##    'isoca':{7:'red',11:'orange',12:'yellow',15:'green',16:'blue',19:'violet',20:'gray',24:'red2',25:'orange2',35:'yellow2',36:'green2',},
    'isoca':{7:'red',12:'orange',16:'yellow',20:'green',25:'blue',36:'violet',},
    'CO2':{1:'red',2:'orange',7:'yellow',8:'green',9:'blue',},
    'H2O':{7:'yellow',8:'green',9:'blue',},
    }
for mode in range(0,3*N):
    print mode+1
    if mode+1 not in d_colors[prefix].keys():
        continue
    if mode+1 == 3*N:
        lines = ['draw color blue\n']
    else:
        lines = ['draw color %s\n' %(d_colors[prefix][mode+1])]
    for i in range(N):
        x = l_coordinates[i][0]
        y = l_coordinates[i][1]
        z = l_coordinates[i][2]
        vx = eigenvectors[mode][3*i+0]
        vy = eigenvectors[mode][3*i+1]
        vz = eigenvectors[mode][3*i+2]
        v_len = math.sqrt(vx**2+vy**2+vz**2)
        v_cone = [
            vx/v_len,
            vy/v_len,
            vz/v_len,
            ]
        v_cylinder = [
            vx-v_cone[0],
            vy-v_cone[1],
            vz-v_cone[2],
            ]
        ## cylinder
        line = 'draw cylinder {%f %f %f} {%f %f %f} radius 0.05\n' %(
            x+v_cylinder[0],y+v_cylinder[1],z+v_cylinder[2],
            x+v_cylinder[0]+v_cone[0],y+v_cylinder[1]+v_cone[1],z+v_cylinder[2]+v_cone[2],
            )
##        lines += [line]
        ## cone
        line = 'draw cone {%f %f %f} {%f %f %f} radius 0.05\n' %(
            x,y,z,
            x+v_cone[0],y+v_cone[1],z+v_cone[2],
            )
        lines += [line]
    fd = open('%s%02i.vmd' %(prefix, mode+1), 'w')
    fd.writelines(lines)
    fd.close()

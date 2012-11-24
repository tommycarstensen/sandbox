import math

def main():

    vmd_arrows()

    vmd_trajectory()

    return


def vmd_arrows(pdb,l_coordinates,eigenvectors,):

    lines = ['draw color white\n']

    for i in range(len(l_coordinates)):

        coordinate = l_coordinates[i]
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        vx = 20*eigenvectors[6][3*i+0]
        vy = 20*eigenvectors[6][3*i+1]
        vz = 20*eigenvectors[6][3*i+2]
        v_len = math.sqrt(vx**2+vy**2+vz**2)

        ## cylinder *and* cone
        if v_len > 1:
            v_cone = [
                20*eigenvectors[6][3*i+0]/v_len,
                20*eigenvectors[6][3*i+1]/v_len,
                20*eigenvectors[6][3*i+2]/v_len,
                ]
            v_cylinder = [
                20*eigenvectors[6][3*i+0]-v_cone[0],
                20*eigenvectors[6][3*i+1]-v_cone[1],
                20*eigenvectors[6][3*i+2]-v_cone[2],
                ]
        ## cone only
        else:
            v_cone = [
                20*eigenvectors[6][3*i+0],
                20*eigenvectors[6][3*i+1],
                20*eigenvectors[6][3*i+2],
                ]
            v_cylinder = [0,0,0,]

        ## cylinder
        if v_len > 1:
            line = 'draw cylinder {%f %f %f} {%f %f %f} radius 0.1\n' %(
                x,y,z,
                x+v_cylinder[0],y+v_cylinder[1],z+v_cylinder[2],
                )
            lines += [line]

        ## cone
        line = 'draw cone {%f %f %f} {%f %f %f} radius 0.15\n' %(
            x+v_cylinder[0],y+v_cylinder[1],z+v_cylinder[2],
            x+v_cylinder[0]+v_cone[0],y+v_cylinder[1]+v_cone[1],z+v_cylinder[2]+v_cone[2],
            )
        lines += [line]

    fd = open('%s.src' %(pdb,), 'w')
    fd.writelines(lines)
    fd.close()

    return


def vmd_trajectory(pdb,l_coordinates,eigenvectors):

    '''visualize the two extreme projections along a trajectory and interpolate n frames between them'''

    d_colors = {6:'red',7:'green',8:'blue',9:'yellow',10:'violet',11:'orange',}
    n_frames = 10

    for mode in range(6,12):

        eigenvector = eigenvectors[mode]

        ## pdb header
        output_vmd = []

        ## loop over frames
        for frame in range(n_frames):
            
##            output_frame = []
            
            ## frame header
##                output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
            output_vmd.append('HEADER    frame t= %4.3f\nMODEL     %4i\n' %(frame,frame+1))
            ## loop over coordinates
            for i_coord in range(len(l_coordinates)):

                coordinate = l_coordinates[i_coord]

                ## starting coordinate
                x1 = coordinate[0]
                y1 = coordinate[1]
                z1 = coordinate[2]

                ## coordinate change
                vx = 32*eigenvector[3*i_coord+0]
                vy = 32*eigenvector[3*i_coord+1]
                vz = 32*eigenvector[3*i_coord+2]

                ## ending coordinte
                x2 = x1+(1-2*float(frame)/float(n_frames))*vx
                y2 = y1+(1-2*float(frame)/float(n_frames))*vy
                z2 = z1+(1-2*float(frame)/float(n_frames))*vz

                sqlength = vx**2+vy**2+vz**2
                len_v = math.sqrt(sqlength)

                record = 'ATOM'
                atom_no = i_coord+1
                atom_name = 'CA'
                altloc = ' '
                res_name = 'UNK'
                chain = 'A'
                res_no = i_coord+1
                iCode = ' '

                ## append atom line
                line = '%6s%5i  %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
                    record.ljust(6),atom_no,
                    atom_name.ljust(3), altloc, res_name, chain, res_no, iCode,
                    x2,y2,z2, 1.0, sqlength,atom_name[0].rjust(2)
                    )
                output_vmd += [line]

            ## terminate frame/MODEL
            output_vmd.append('TER\nENDMDL\n')

        fd = open(pdb+'_mode'+str(mode+1).zfill(2)+'.pdb', 'w')
        fd.writelines(output_vmd)
        fd.close()

    return


if __name__ == '__main__':
    main()

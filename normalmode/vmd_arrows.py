def main():

    l_eigenvectors = parse_eigenvectors()

    l_coordinates = parse_coordinates()

    vmd_arrows(l_eigenvectors,l_coordinates,)

    return
    

def parse_eigenvectors():

    import Numeric

    fd = open('eigcomp.xvg','r')
    lines = fd.readlines()
    fd.close()

    d_rmsfs = {1:[],2:[],3:[],}

    for line in lines:
        if line[0] == '@':
            axis = 0
            continue
        elif line[0] == '&':
            if axis == 3:
                break
            axis += 1
        elif axis == 0:
            continue
        else:
            d_rmsfs[axis] += [float(line.split()[1])]

    l_eigenvectors = []
    for i in range(len(d_rmsfs[1])):
        x = d_rmsfs[1][i]
        y = d_rmsfs[2][i]
        z = d_rmsfs[3][i]
        l_eigenvectors += [Numeric.array([x,y,z,])]

    return l_eigenvectors


def parse_coordinates():

    import Numeric

    fd = open('/data/remediated_pdb/at/pdb1atp.ent','r')
    lines = fd.readlines()
    fd.close()

    l_coordinates = []

    for line in lines:
        record = line[:6].strip()
        if record == 'ATOM':
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            chain = line[21]
            if chain != 'E':
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coordinate = Numeric.array([x,y,z])
            l_coordinates += [coordinate]
            
    return l_coordinates


def vmd_arrows(l_eigenvectors, l_coordinates,):

    import math

    lines = ['draw color white\n']

    for i in range(len(l_coordinates)):

        coordinate = l_coordinates[i]
        eigenvector = 20*l_eigenvectors[i]

        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        vx = eigenvector[0]
        vy = eigenvector[1]
        vz = eigenvector[2]
        v_len = math.sqrt(vx**2+vy**2+vz**2)

        ## cylinder *and* cone
        if v_len > 1:
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
        ## cone only
        else:
            v_cone = [
                vx,
                vy,
                vz,
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

    fd = open('una.vmd','w')
    fd.writelines(lines)
    fd.close()

    return


if __name__=='__main__':
    main()

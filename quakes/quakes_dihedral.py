import math,numpy

def main(c1,c2,c3,c4):

    v1 = c2-c1
    v2 = c3-c2
    v3 = c4-c3

    angle = math.atan2(
        numpy.dot(
            math.sqrt(sum(v2*v2))*v1,
            cross(v2,v3),
            ),
        numpy.dot(
            cross(v1,v2),
            cross(v2,v3),
            ),
        )
    angle *= 180./math.pi

    return angle


def cross(v1,v2):

    n = numpy.array([
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0],
        ])

    return n

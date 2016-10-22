#!/usr/bin/env python3

import numpy
import random
import timeit
import profile


s1 = '''
def multiply_numpy_dot():

    m1 = numpy.random.rand(3, 3)
    m2 = numpy.random.rand(3, 3)

    return numpy.dot(m1, m2)
'''

s2 = '''
def multiply_pep465():

    m1 = numpy.random.rand(3, 3)
    m2 = numpy.random.rand(3, 3)

    return m1 @ m2
'''

s3 = '''
def multiply_numpy_matmul():

    m1 = numpy.random.rand(3, 3)
    m2 = numpy.random.rand(3, 3)

    return numpy.matmul(m1, m2)
'''


def main():

    m1 = numpy.random.rand(3, 4)
    m2 = numpy.random.rand(4, 3)

    m3a = numpy.dot(m1, m2)
    m3b = m1 @ m2
    m3c = numpy.matmul(m1, m2)

    print(m3a)
    print(m3b)
    print(m3c)

    n = 10**7
    print(timeit.timeit(s1, number=n))
    print(timeit.timeit(s2, number=n))
    print(timeit.timeit(s3, number=n))

    return

if __name__ == '__main__':
    main()

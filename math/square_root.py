def main(x,precision,):

    ## http://www.homeschoolmath.net/teaching/square-root-algorithm.php

    index = str(x).index('.')
    x2 = str(x*10**(precision-index))

    sqrt = ''
    d_sqrt = {
        1:1,
        4:2,
        9:3,
        16:4,
        25:5,
        36:6,
        49:7,
        64:8,
        81:9,
        100:10,
        }
    l_squares = d_sqrt.keys()
    l_squares.sort()
    
    for i in range(str(x2)):
        x3 = str(x2)[i:i+2]
        x2 = d_sqrt[int(s)]
        sqrt += str(x2)
        remainder = int(s)-x2
        print sqrt
        print s, x2, x3
        print remainder
        stop
        if len(sqrt) == precision:
            break

    sqrt = x/2.
    iteration = 0
    while str(int(10000*sqrt))[-1] == '0':
        iteration += 1
        if sqrt**2 > x:
            index = str(int(10000*sqrt)).index('0')
            sqrt -= 1./(10.**(index-1))
        else:
            index = str(int(10000*sqrt)).index('0')-1
            sqrt += 1./(10.**(index))
            print sqrt
            print sqrt**2
            stop2

    print iteration
    print sqrt**2

    return sqrt


if __name__ == '__main__':
    x = 10.2
    precision = 5
    main(x,precision,)

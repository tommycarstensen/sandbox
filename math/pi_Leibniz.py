## a small script for calculating Pi using formula of Leibniz

## pi = 4 * (1 - 1/3 + 1/5 - 1/7 + ...)

import time

## method calculating when last digit has been reached
t1 = time.time()
add = 1
i = 0
pi = 0
while abs(add) > 0.00001:
    add = 4*(
        (-1.)**i ## sign
        /
        (2*i+1) ## denominator
        )
    i += 1
    pi += add
##print i, add
t2 = time.time()
print 't', t2-t1

## method doing given xrange
t1 = time.time()
print 4*sum(
    (
        (-1.)**i ## sign
        /
        (2*i+1) for i in xrange(5**8)
        )
    )
t2 = time.time()
print 't', t2-t1

print xrange(4)
print 5**8

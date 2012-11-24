'''test2'''
x = 3

def odd(x):

    return x % 2 != 0

def func(x,):
    for x in range(3):
        yield 3-x

def main():

    import re

    '''test'''

    iterator = func(range(3))

##    x = [1,2,3]
##    y = ['a','b','c',]
##
##    print zip(x,y,)
##    print dict(zip(x,y,))

##    text = 'ha ha ha'
##    keyword = re.compile(r'ha (ha )*')
##    result = keyword.search(text)
##    if result:
##        print '3.4', result.groups()
##
##    print 10%2
##
##    print 10+'a'

##    d = {}
##    d[('a','b')] = 'a'
##    d[main] = 'c'

##    import math
##    print 'a', filter(lambda(x) : math.sqrt(x), range(0,100))
####    print 'b', math.sqrt([0,1,99])
##    print 'c', map(math.sqrt,range(100))
##    print 'd', apply(math.sqrt(), list(100))

##    x = 1
##    y = 0
##    z = 0
##    try:
##        z += 1
##        a = x/y
##    except ZeroDivisionError:
##        z += 1
##    finally:
##        z += 1
##    print z

##    l = [-1,-4,2,5,-3]
####    print 'b', sort(l,key=abs)
##    print 'd', sorted(l,key=abs)
##    print 'e', map(lambda x : x[1], sorted(map(lambda x: (abs(x),x), l)))

##    print 'foo'[:-1]

##    import os
##    import stat
##    import fileinfo
##    print 'a', os.stat('alice.txt')[stat.ST_SIZE]
##    print 'd', os.stat.size('alice.txt')

##    import traceback, sys, pdb
##    print 'a'
##    traceback.excepthook = debug
##    print 'b'
##    sys.excepthook = debug
##    print 'c'
##    sys.exitfunc = debug
##    print 'd'
##    raise debug
##    print 'e'
##    pdb.runcall = debug

##    x = isinstance(u'34\xc2\xb0',str)
##    print x
##    y = 'abc'.decode('ascii')
##    print y
##    print type(y)
##    z = ''.join([u'34\xb0',u"56'",])
##    print type(z)
##    w = 'a'+u'b'
##    print type(w)

##    import os
##    message = 'd'
##    fh = open('messages.log','r')
####    fh.seek(0,os.SEEK_END)
##    fh.write(message+'\n')

##    data = range(7)
##    print '0', map(f1,filter(f2,data))
##    print 'b', map(lambda r: r*3, filter(lambda f: f>3 or 0,data))
##    print 'e', map(lambda r: r*3, filter(f2,data))

##    stringlengths = {1:'a',3:'ccc',4:'dddd',11:'abcdefghijk',}
##    print 'a', dict(filter(oklength,stringlengths.items()))
####    for i in stringlengths:
####        if stringlengths[i] > 10 or stringlengths[i] < 3:
####            del stringlengths[i]
####    print 'b', stringlengths
####    for i in stringlengths.keys():
####        if stringlengths[i] > 10 or stringlengths[i] < 3:
####            del stringlengths[i]
####    print 'c', stringlengths
####    print 'd', stringlengths.delete_if(lambda value : value < 3 or value > 10)

##    x = reduce(lambda x,y : [y] + x, range(5), [])
##    print x
   
    return

def oklength((s,l)):
    return l >= 3 and l <= 10

##def foo(x) {
##    return x
##    }    

def f1(x):
    return 3*x

def f2(x):
    try:
        return x > 3
    except:
        return 0

def debug(e,i,tb,):
    traceback.print_exception(e,i,tb)
    pdb.pm()

if __name__ == '__main__':
    main()

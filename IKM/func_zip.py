'''test2'''
x = 3

def odd(x):

    return x % 2 != 0

def func(x,):
    for x in range(3):
        yield 3-x

def main():

    '''test'''

    iterator = func(range(3))

    x = [1,2,3]
    y = ['a','b','c',]

    print zip(x,y,)
    print dict(zip(x,y,))
    
    return

if __name__ == '__main__':
    main()

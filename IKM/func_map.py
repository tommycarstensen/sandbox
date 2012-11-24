'''test2'''
x = 3

class aaa():

    def odd(self,x):

        return x % 2 != 0

    def square(self,x):
        return x*x

    def main(self,):

        '''test'''

        x = 3.1
        y = x
        print filter(self.odd, range(10))

        print globals()
        print hash(x)
        print id(x)
        print id(y)

        print issubclass(aaa,aaa)

        print map(self.square,range(10))
        
        return

if __name__ == '__main__':
    instance = aaa()
    instance.main()

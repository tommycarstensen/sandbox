'''test2'''
x = 3

class aaa():

    def odd(self,x):

        return x % 2 != 0

    def func(self,):
        for x in range(3):
            yield 3-x

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

        print max(1,2,3,key=self.func)

        iterator = self.func()
        print next(iterator)
        
        return

if __name__ == '__main__':
    instance = aaa()
    instance.main()

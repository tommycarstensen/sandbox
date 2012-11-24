'''test2'''
x = 3

class aaa():

    def odd(self,x):

        return x % 2 != 0

    def func(self,x,):
        for x in range(3):
            yield 3-x

    def main(self,):

        '''test'''

        x = 3.1
        y = x
        print filter(self.odd, range(10))

        iterator = self.func(range(3))
        print next(iterator)
        print next(iterator)
        
        return

if __name__ == '__main__':
    instance = aaa()
    instance.main()

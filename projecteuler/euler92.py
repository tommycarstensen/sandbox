import cProfile

def main():
    set89 = set([89])
    set1 = set([1])
    cnt89 = 0
    d_squares = {'0':0,'1':1,'2':4,'3':9,'4':16,'5':25,'6':36,'7':49,'8':64,'9':81}
    for i in range(2,10**7):
        setsquares = set([i])
        while True:
            sumsquares = sum(d_squares[j] for j in str(i))
            i = sumsquares
            setsquares |= set([i])
            if i in set1:
                set1 |= setsquares
                break
            if i in set89:
                set89 |= setsquares
                cnt89 +=1
                break

    print(cnt89)

cProfile.run('main()')

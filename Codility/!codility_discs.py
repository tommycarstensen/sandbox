def main():

    # write your code here

    A = [1,5,2,1,4,0,]

    count = 0
    for i in range(len(A)-1):
        r1 = A[i]
        for j in range(i+1,len(A)):
            r2 = A[j]
            dist_sq = (i-j)**2
            if dist_sq <= r1**2+r2**2:
                count += 1
            if count > 10000000:
                return -1

    return count

if __name__ == '__main__':
    print main()

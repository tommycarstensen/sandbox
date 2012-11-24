def main():

    # write your code here

    A = [2,2,1,0,4,4,1,]

    set_A = set(A)

    for i in range(len(A)):
        set_A -= set([A[i]])
        if len(set_A) == 0:
            break

    return i

if __name__ == '__main__':
    print main()

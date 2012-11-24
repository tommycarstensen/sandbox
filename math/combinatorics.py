def main():

    n = 8; r = 5
    n = 3; r = 2
    permutations = permutation(n,r,)
    print permutations

    n = 13; r = 3
    n = 3; r = 2
    combinations = combination(n,r,)
    print combinations
    
    return


def faculty(n):

    fac = 1
    for i in range(1,n+1):
        fac *= i

    return fac


def permutation(n,r,):

    ## counting rule for permutations
    ## the number of ways to arrange in order n different objects within r positions

    nPr = permutations = faculty(n)/faculty(n-r) ## = n*(n-1)*...*(n-r+1)*(n-r)

    return permutations


def combination(n,r,):

    ## the number of ways to pick in random order r different objects from n different objects

    nCr = combinations = faculty(n)/(faculty(r)*faculty(n-r))

    return combinations


def perm_wo_rep(group):

    '''Permutation without repetition'''

    count_group = len(group)
    count_combinations = faculty(count_group)
    combinations = []

    poplist = []

    if count_group > 6:
        stopandfindanothersolution

    for i in range(count_combinations):
        combinations += [[]]
        poplist += [list(group)]

    i = 0
    ## loop over n
    for i in range(count_group,0,-1):
        i_fac = faculty(i-1)
        j = 0
        ## loop over n!
        while j < count_combinations:
            ## loop over i ... to loop over i!=i*(i-1)
            for k in range(i):
                ## loop over (i-1)! ... to loop over i!=i*(i-1)!
                for l in range(i_fac):
                    combinations[j].append(poplist[j].pop(k))
                    j += 1

    return combinations


def permutation_sub(combinations1):

    '''return list of combinations'''

    combinations2 = []

    for group in combinations1:
        
        count_group = len(group)
        count_permutations = faculty(count_group)
        permutations = []

        poplist = []

        if count_group > 6:
            stopandfindanothersolution

        for i in range(count_permutations):
            permutations += [[]]
            poplist += [list(group)]

        i = 0
        ## loop over n
        for i in range(count_group,0,-1):
            i_fac = faculty(i-1)
            j = 0
            ## loop over n!
            while j < count_permutations:
                ## loop over i ... to loop over i!=i*(i-1)
                for k in range(i):
                    ## loop over (i-1)! ... to loop over i!=i*(i-1)!
                    for l in range(i_fac):
                        permutations[j].append(poplist[j].pop(k))
                        j += 1

        combinations2 += permutations

    return combinations2


def permutation_wo_rep(group,size):

    '''Permutation without repetition'''

##    group = ['A','B','C','D',]
##    size = 3

    binary2 = size*[1]+(len(group)-size)*[0]
    binary = (len(group)-size)*[0]+size*[1]
    binaries = []
    while binary != binary2:
        if sum(binary) == size:
            binaries += [list(binary)]
        if binary[-1] == 0:
            binary[-1] = 1
            continue
        else:
            for i in range(len(binary)-1,-1,-1):
                if binary[i] == 0:
                    binary[i] = 1
                    break
                else:
                    binary[i] = 0
    binaries += [binary2]

    combinations = []
    for binary in binaries:
        combination = []
        for i in range(len(binary)):
            if binary[i] == 1:
                combination += [group[i]]
        combinations += [combination]

    combinations = permutation_sub(combinations)

    return combinations            


def permutation_w_rep(group,size):

    '''Permutation with repetition'''

##    group = ['A','B','C','D',]
##    size = 3

    binary2 = size*[len(group)-1]
    binary = size*[0]
    binaries = []
    while binary != binary2:
        binaries += [list(binary)]
        for i in range(size-1,-1,-1):
            if binary[i] != len(group)-1:
                binary[i] += 1
                break
            else:
                binary[i] = 0
                continue
    binaries += [binary]

    combinations = []
    for binary in binaries:
        combination = []
        for i in range(len(binary)):
            combination += [group[binary[i]]]
        combinations += [combination]

    return combinations            


if __name__ == '__main__':
    permutation_w_rep([
        'A','B','C',
        'D',
        ])
    main()

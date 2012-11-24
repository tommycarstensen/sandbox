## 2006

m = [
    [0,2,0,0,0,4,9,0,7],
    [0,7,9,8,0,5,0,2,0],
    [0,1,0,2,7,9,5,0,6],
    [7,0,2,5,4,0,0,1,9],
    [9,0,4,0,2,0,7,0,5],
    [3,5,1,9,8,7,0,0,0],
    [2,0,0,7,0,0,0,9,0],
    [0,9,0,4,0,2,3,5,0],
    [0,0,0,0,0,0,0,7,0],
    ]
m = [
    [4,5,8,1,0,0,0,2,0],
    [2,6,1,0,0,8,0,0,0],
    [9,7,3,0,5,0,0,1,0],
    [0,4,0,9,0,0,2,0,1],
    [1,8,2,0,0,5,0,9,0],
    [3,9,7,0,0,1,0,0,0],
    [0,1,0,0,7,0,9,0,2],
    [7,0,0,5,0,0,1,4,0],
    [0,0,0,0,1,0,0,0,7],
    ]

for r in range(len(m)):
    for c in range(len(m[r])):
        if m[r][c] == 0:
            m[r][c] = set((1,2,3,4,5,6,7,8,9))
        else:
            m[r][c] = set((m[r][c],))

for i in range(100):
    ## loop over major rows and columns
    for R1 in range(3):
        for C1 in range(3):
            ## loop over minor rows and columns
            for r1 in range(3):
                for c1 in range(3):
                    ## eliminate possible numbers if not determined
                    if len(m[3*R1+r1][3*C1+c1]) != 1:
                        ## eliminate from rows
                        for R2 in range(3):
                            for r2 in range(3):
                                if len(m[3*R2+r2][3*C1+c1]) == 1 and not (R1 == R2 and r1 == r2):
                                    m[3*R1+r1][3*C1+c1] = m[3*R1+r1][3*C1+c1]-m[3*R2+r2][3*C1+c1]
                                if len(m[3*R1+r1][3*C1+c1]) == 1:
                                    print 'rows', 3*R1+r1+1, 3*C1+c1+1, m[3*R1+r1][3*C1+c1]
                                    break
                            if len(m[3*R1+r1][3*C1+c1]) == 1:
                                break
                        if len(m[3*R1+r1][3*C1+c1]) == 1:
                            break
                        ## eliminate from columns
                        for C2 in range(3):
                            for c2 in range(3):
                                if len(m[3*R1+r1][3*C2+c2]) == 1 and not (C1 == C2 and c1 == c2):
                                    m[3*R1+r1][3*C1+c1] = m[3*R1+r1][3*C1+c1]-m[3*R1+r1][3*C2+c2]
                                    if len(m[3*R1+r1][3*C1+c1]) == 1:
                                        print 'columns', 3*R1+r1+1, 3*C1+c1+1, m[3*R1+r1][3*C1+c1]
                                        break
                            if len(m[3*R1+r1][3*C1+c1]) == 1:
                                break
                        if len(m[3*R1+r1][3*C1+c1]) == 1:
                            break
                        ## eliminate from box
                        for r3 in range(3):
                            for c3 in range(3):
                                if len(m[3*R1+r3][3*C1+c3]) == 1 and not (r1 == r3 and c1 == c3):
                                    m[3*R1+r1][3*C1+c1] = m[3*R1+r1][3*C1+c1]-m[3*R1+r3][3*C1+c3]
                                    if len(m[3*R1+r1][3*C1+c1]) == 1:
                                        print 'box', 3*R1+r1+1, 3*C1+c1+1, m[3*R1+r1][3*C1+c1]
                                        break
                            if len(m[3*R1+r1][3*C1+c1]) == 1:
                                break
                        if len(m[3*R1+r1][3*C1+c1]) == 1:
                            break
                        ## add to box if only 1 position possible
                        box = set()
                        for r3 in range(3):
                            for c3 in range(3):
                                if not (r1 == r3 and c1 == c3):
    ##                                if 3*R1+r1 == 7 and 3*C1+c1 == 0:
    ##                                    print m[7][0], box, r3, c3
    ##                                    print m[6][0]
    ##                                    print m[6][1]
    ##                                    print m[6][2]
    ##                                    print m[7][0]
    ##                                    print m[7][1]
    ##                                    print m[7][2]
    ##                                    print m[8][0]
    ##                                    print m[8][1]
    ##                                    print m[8][2]
    ##                                    stop
                                    box = box | m[3*R1+r3][3*C1+c3]
                            if len(m[3*R1+r1][3*C1+c1]) == 1:
                                break
                        if len(m[3*R1+r1][3*C1+c1]-box) == 1:
                            m[3*R1+r1][3*C1+c1] = m[3*R1+r1][3*C1+c1]-box

                    ## eliminate possible numbers if determined
    ##                else:
    ##                    for r2 in range(9):
    ##                        if r1 != r2:
    ##                            m[r2][c1] = m[r2][c1]-m[r1][c1]
    ##                    for c2 in range(9):
    ##                        if c1 != c2:
    ##                            m[r1][c2] = m[r1][c2]-m[r1][c1]


for i in range(9):
    print m[i]
print ' _______ _______ _______'
print '|       |       |       |'
for i in range(9):
    print '|',
    for j in range(9):
        if j % 3 == 0 and j != 0:
            print '|',
        print list(m[i][j])[0],
        if len(m[i][j]) > 1:
            print
            print m[i][j]
            stop
    print '|'
    if i in [2,5]:
        print '|_______|_______|_______|'
        print '|       |       |       |'
print '|_______|_______|_______|'

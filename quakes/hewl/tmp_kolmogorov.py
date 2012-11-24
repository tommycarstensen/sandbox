l1 = [8.8,8.4,7.9,8.7,9.1,9.6,]
l2 = [9.9,9.0,11.1,9.6,8.7,10.4,9.5,]

l = l1+l2
l.sort()

l_set = list(set(l))
l_set.sort()

n1 = float(len(l1))
n2 = float(len(l2))

Fi1 = 0
Fi2 = 0
D = 0
for i in range(len(l_set)):

    xi = l_set[i]
    ## frequency
    if xi in l1:
        fi1 = l1.count(xi)
    else:
        fi1 = 0
    if xi in l2:
        fi2 = l2.count(xi)
    else:
        fi2 = 0
    ## cumulative frequency
    Fi1 += fi1
    Fi2 += fi2
    ## cumulative relative frequency
    rel_Fi1 = Fi1/n1
    rel_Fi2 = Fi2/n2
    print rel_Fi1, rel_Fi2
    ## KS statistic
    Di = abs(rel_Fi1-rel_Fi2)
    if Di > D:
        D = Di

print D
print n1,n2
print 30./42.
    

l = [1.4,2.6,3.3,4.2,4.7,5.6,5.6,6.4,7.7,9.3,10.6,11.5,12.4,18.6,22.3,]

l.sort()
l_set = list(set(l))
l_set.sort()

n = float(len(l))

l_fi = []
Fi = 0
relFi = Fi/n
l_relFi = [relFi]
Min = float(min(l_set))
Max = float(max(l_set))
Min = 0
Max = 25
D = 0
for i in range(len(l_set)):

    ## frequency
    fi = l.count(l_set[i])
    ## cumulative frequency
    Fi += fi
    ## cumulative relative frequency
    relFi = Fi/n
    l_relFi += [relFi]
    ## cumulative relative expected frequency
    relFi_circumflex = (l_set[i]-Min)/Max
    ## KS statistic
    Di = abs(relFi-relFi_circumflex)
    Di_mark = abs(l_relFi[i+1-1]-relFi_circumflex)
    if Di > D:
        D = Di
    if Di_mark > D:
        D = Di_mark

print D

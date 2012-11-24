import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import statistics

l_MDs = ['MD_2vb1_Glu35Asp52_00','MD_2vb1_Glu35Asp52_01',]

l = [[],[],[],[],[],]
for MD in l_MDs:
    for i in range(5):
        l[i] += [[]]
    fd = open('35C/%s/energies.txt' %(MD),'r')
    lines = fd.readlines()
    fd.close()

    for line in lines:
        for i in range(4):
            l[i][-1] += [float(line.split()[i+1])]
        Sum = float(line.split()[1])+float(line.split()[2])+float(line.split()[3])+float(line.split()[4])
        l[4][-1] += [Sum]
        if MD == 'MD_2vb1_Glu35Asp52_01' and len(l[0][-1]) == 500:
            break

Sum = 0
##l = [[8.8,8.4,7.9,8.7,9.1,9.6,],[9.9,9.0,11.1,9.6,8.7,10.4,9.5],]
mean1,mean2,stderr,p = statistics.twosamplettest(l[0][0],l[0][1],verbose=False,)
Sum += mean1-mean2
print 'Asp52', mean1-mean2,1.962*stderr
mean1,mean2,stderr,p = statistics.twosamplettest(l[1][0],l[1][1],verbose=False,)
Sum += mean1-mean2
print 'protein', mean1-mean2,1.962*stderr
mean1,mean2,stderr,p = statistics.twosamplettest(l[2][0],l[2][1],verbose=False,)
Sum += mean1-mean2
print 'chloride', mean1-mean2,1.962*stderr
mean1,mean2,stderr,p = statistics.twosamplettest(l[3][0],l[3][1],verbose=False,)
Sum += mean1-mean2
print 'water', mean1-mean2,1.962*stderr
mean1,mean2,stderr,p = statistics.twosamplettest(l[4][0],l[4][1],verbose=False,)
print 'sum', mean1-mean2,1.962*stderr

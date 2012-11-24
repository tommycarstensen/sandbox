import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import  statistics

l1 = [48.2,54.6,58.3,47.8,51.4,52.0,55.2,49.1,49.9,52.6,]
l2 = [52.3,57.4,55.6,53.2,61.3,58.0,59.8,54.8,]

mean1,mean2,stderr,p = statistics.twosamplettest(l1,l2,)

print p

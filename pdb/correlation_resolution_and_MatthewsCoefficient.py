import sys
sys.path.append('/home/tc/svn/tc_sandbox/misc')
import statistics

fd = open('/home/tc/svn/tc_sandbox/pdb/db_MatthewsCoefficient.txt','r')
lines = fd.readlines()
fd.close()
d_MV = {}
for line in lines:
    l = line.strip().split()
    pdb = l[0]
    v = l[1]
    if float(v) > 10:
        print pdb, 'MV', v
        continue
    if float(v) < 1.01:
        print pdb, 'MV', v
        continue
    d_MV[pdb] = v

fd = open('/home/tc/svn/tc_sandbox/pdb/db_resolution.txt','r')
lines = fd.readlines()
fd.close()
d_resolution = {}
for line in lines:
    l = line.strip().split()
    pdb = l[0]
    v = l[1]
    if v == "['.']":
        continue
    if float(v[2:-2]) > 5.0:
##        print pdb, 'resolution', v
        continue
    d_resolution[pdb] = v

l_MV = []
l_resolution = []
lines_gnuplot = []
set_pdbs = set(d_MV.keys())&set(d_resolution.keys())
for pdb in set_pdbs:
    MV = float(d_MV[pdb])
    resolution = float(d_resolution[pdb][2:-2])
    l_MV += [MV]
    l_resolution += [resolution]
    lines_gnuplot += ['%s %s\n' %(resolution,MV,)]

t = statistics.do_regression(l_resolution,l_MV)

fd = open('tmp.txt','w')
fd.writelines(lines_gnuplot)
fd.close()

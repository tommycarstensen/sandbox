s = '''
1l35
1l55
1l57
1l59
1l61
1l62
1l63
1lpy
2lzm
4lzm
5lzm
6lzm
1pqj
1xei
1xej
1xek
3f8v
3f9l
1d3s
1dbw
1dck
1eqd
1eyq
1q2c
2uwu
2v1r
'''

l_pdbs = s.split()

import urllib2


print 'PDB\tDEPOSITION\tREVISION\n'

for pdb in l_pdbs:

    print
##    print pdb
##    urllines = urllib2.urlopen('http://www.pdb.org/pdb/files/%s.pdb' %(pdb))
    urllines = urllib2.urlopen('http://www.pdb.org/pdb/files/%s.pdb?headerOnly=YES' %(pdb))
    lines = urllines.readlines()
    for i in range(len(lines)):
        line = lines[i]
        record = line[:6].strip()

        if record == 'HEADER':
            print pdb, '\t',
            print line[50:59], '\t',
        if record == 'REVDAT' and lines[i-1][:6] != 'REVDAT':
            print line[13:22],

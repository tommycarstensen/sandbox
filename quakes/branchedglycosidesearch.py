import os

pdbpath = '/oxygenase_local/data/pdb/'

l_hetids2 = []
l_hetids1 = []

subdirs = os.listdir(pdbpath)
subdirs.sort()
for subdir in subdirs:
    print subdir,
    files = os.listdir(pdbpath+subdir)
    for file in files:
        fd = open('%s%s/%s' %(pdbpath,subdir,file),'r')
        lines = fd.readlines()
        fd.close()
        hetIDs = set([])
        for line in lines:
            if line[:6].strip() == 'HET':
                hetID = line[7:10].strip()
                hetIDs |= set([hetID])
        if len( set(['MAN','GLC','NAG']) & hetIDs ) == 3:
            print file, hetIDs

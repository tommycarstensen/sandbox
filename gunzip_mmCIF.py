import os,sys

##os.system('rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/mmCIF/')

path = '/media/Elements/mmCIF'
path = '/media/Tommy/mmCIF'
path = '/media/Tommy/pdb'
path = '/media/39b8dbfd-b37f-472f-a030-c86a6e75c0d9/1TB/biounit'
path = '/media/WDMyBook1TB/2TB/biounit'
path = '/media/WDMyBook1TB/2TB/pdb'
path = '/media/WDMyBook1TB/2TB/mmCIF'

## for removal purposes
if 'pdb' in path:
    extension = 'ent'
elif 'mmCIF' in path:
    extension = 'cif'
else:
    print path
    stop_extension

l_dn = os.listdir(path)
l_dn.sort()
for dn in l_dn:
    if not os.path.isdir('%s/%s' %(path,dn,)):
        continue
    print dn
    l_fn = os.listdir('%s/%s' %(path,dn,))
    l_fn.sort()
    for fn in l_fn:
        ## remove cif/ent
        if fn[-4:] == '.%s' %(extension):
            ## if gz no longer part of repository
            if not os.path.isfile('%s/%s/%s.gz' %(path,dn,fn,)):
                print 'removing', fn
                os.remove('%s/%s/%s' %(path,dn,fn))
            ## if zero size
            elif os.path.getsize('%s/%s/%s' %(path,dn,fn)) == 0:
                print 'size = 0', fn
                print 'removing', fn
                os.remove('%s/%s/%s' %(path,dn,fn))
##            ## if only part of file
##            elif os.popen('tail -1 %s/%s/%s' %(path,dn,fn)).read() != 'END%s\n' %(77*' '):
##                print os.popen('tail -1 %s/%s/%s' %(path,dn,fn)).read()
##                print 'removing', fn
##                os.remove('%s/%s/%s' %(path,dn,fn))
        ## gunzip
        elif fn[-3:] == '.gz':

            ## keep or delete old file
            if os.path.isfile('%s/%s/%s' %(path,dn,fn[:-3],)):

                ## either
                ## 1) keep old files and continue
                if not '--delete' in sys.argv:
                    continue

                ## or
                ## 2) delete old files and gunzip
                ## if gz was changed, then pdb/cif needs to be changed as well
                os.remove('%s/%s/%s' %(path,dn,fn[:-3],))

            ## gunzip file without removing zipped file
            s = 'gunzip -c %s/%s/%s -> %s/%s/%s' %(
                path,dn,fn,
                path,dn,fn[:-3],
                )
            os.system(s)

##        elif fn[-4:] == '.pdb':
##            os.remove('%s/%s/%s' %(path,dn,fn))

##        else:
##            print 'delete', fn ,'? (y/n)'
##            s = raw_input()
##            if s == 'y':
##                os.remove('%s/%s/%s' %(path,dn,fn))

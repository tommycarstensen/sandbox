import os
pdbpath = '/oxygenase_local/data/pdb/'

d_siterecordsmissing = {}
d_remark800recordsmissing = {}
l_missinglinebreaks = []

##pdbs = []
##fd = open('REMARK800.txt','r')
##lines = fd.readlines()
##fd.close()
##for line in lines:
##    pdbs += [line[:4]]
##fd = open('REMARK800.txt','r')
##lines = fd.readlines()
##fd.close()
##for line in lines:
##    pdbs += [line[:4]]
##pdbs.sort()

subdirs = os.listdir(pdbpath)
subdirs.sort()
for subdir in subdirs:
    print subdir
    files = os.listdir(pdbpath+subdir)
    files.sort()
    for file in files:
        pdb = file[3:7]
        if pdb in [
            '2bl2','1ha5','2jer',
            ]:
            continue
        SITE_siteIDs = set()
        REMARK800_siteIDs = set()
        fd = open(pdbpath+subdir+'/'+file,'r')
        lines = fd.readlines()
        fd.close()
        for i in range(len(lines)):
            line = lines[i]
            if line[6:].strip() == 'ATOM':
                break
            if line[:10] == 'REMARK 800':
                if line[11:27] == 'SITE_IDENTIFIER:':
                    siteIDs = line[28:80].strip().split(',')
                    if siteIDs[-1][-1] == '.':
                        siteIDs[-1] = siteIDs[-1][:-1]
                    if siteIDs[-1] != 'AND' and siteIDs[-1].split()[0] == 'AND' and len(siteIDs[-1].split()) == 2:
                        siteIDs[-1] = siteIDs[-1].split()[1]
                    if len(siteIDs) == 1 and len(siteIDs[0].split()) == 3 and siteIDs[0].split()[-2] == 'AND':
                        siteIDs = [siteIDs[0].split()[0],siteIDs[0].split()[2]]
                    for j in range(len(siteIDs)):
                        siteID = siteIDs[j]
                        siteID = siteID.strip()
                        siteIDs[j] = siteID
                        if len(siteID) > 3 or len(siteID) < 1:
                            print pdb, siteID, siteIDs
                            stop
                    REMARK800_siteIDs |= set(siteIDs)
                elif line[11:28] == 'SITE_DESCRIPTION:':
                    site_description = line[29:80].strip()
                    for j in range(i+1,len(lines)):
                        line = lines[j]
                        if line[:10] != 'REMARK 800' or line.strip() in ['REMARK 800','REMARK 800 SITE'] or line[11:28] == 'SITE_DESCRIPTION:' or line[11:27] == 'SITE_IDENTIFIER:':
                            if line[11:28] == 'SITE_DESCRIPTION:' or line.strip == 'REMARK 800 SITE':
                                print pdb, line
                                notexpected
                            line = lines[i]
                            break
                        site_description += ' '+line[11:80].strip()
                        if line[10] != ' ':
                            print line
                            notexpected
                        if 'IDENTIFIER:' in site_description:
                            l_missinglinebreaks += [pdb]
##                    if site_description in ['','NULL']:
##                        continue

                elif line.strip() in ['REMARK 800','REMARK 800 SITE']:
                    continue
            if line[:6].strip() == 'SITE':
                site_name = line[11:14].strip()
                SITE_siteIDs |= set([site_name])

        if pdb not in [
##            ## REMARK800 records deleted during remediation
##            '43ca','43c9','257l','260l','8a3h',
##            ## REMARK800 records deleted during remediation (nonremediated site description = NULL)
##            '1a7x','1a7v','1a9p','1a9t',
##            ## remediation shortening of site identifier line
##            '1a9x',
##            ## REMARK800 missing
##            '117e','11ba',
##            ## REMARK800 missing (not converted from REMARK5)
##            '172l','8aat',
##            ## SITE record missing
##            '2a97',
##            ## SITE record missing (should not be REMARK800 record)
##            '1914',
            ]:
            if len(REMARK800_siteIDs) == 1 and len(SITE_siteIDs) == 0:
                d_siterecordsmissing[pdb] = REMARK800_siteIDs
            elif len(REMARK800_siteIDs) == 0 and len(SITE_siteIDs) == 1:
                d_remark800recordsmissing[pdb] = SITE_siteIDs
            elif REMARK800_siteIDs-SITE_siteIDs != set():
                d_siterecordsmissing[pdb] = REMARK800_siteIDs-SITE_siteIDs
##                print REMARK800_siteIDs,SITE_siteIDs
##                print REMARK800_siteIDs-SITE_siteIDs
##                print SITE_siteIDs-REMARK800_siteIDs
##                SITErecordmissing
            elif SITE_siteIDs-REMARK800_siteIDs != set():
                d_remark800recordsmissing[pdb] = SITE_siteIDs-REMARK800_siteIDs
##                print REMARK800_siteIDs,SITE_siteIDs
##                print REMARK800_siteIDs-SITE_siteIDs
##                print SITE_siteIDs-REMARK800_siteIDs
##                REMARK800recordmissing
            elif REMARK800_siteIDs == set() and SITE_siteIDs == set():
                continue
            elif REMARK800_siteIDs == SITE_siteIDs:
                continue
            else:
                print REMARK800_siteIDs,SITE_siteIDs
                stop

print d_siterecordsmissing
print d_remark800recordsmissing
l = list(set(l_missinglinebreaks))
l.sort()
print l

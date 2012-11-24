def main():

##    l_pdb_pairs = main_loop()

    identify_pdb_pairs_with_ligand_differences(l_pdb_pairs)

    return


def main_loop():

    l_pdbs = []
    if os.path.isfile('molmovdb.txt'):
        fd = open('molmovdb.txt','r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            l = line.split('\t')
            l[-1] = l[-1].strip()
            l_pdbs += [l]
        if len(lines) > 0: l_continue = [True,l[0],]
        else: l_continue = [False,None,]
    else: l_continue = [False,None,]

    ## came to and end already, don't read url
    if l[0] == 'trna':
        return l_pdbs

    url = 'http://www.molmovdb.org/cgi-bin/browse.cgi'
    urllines = urllib2.urlopen(url)

    ## lines to loop
    lines = urllines.readlines()

    for i1 in range(len(lines)):

        if not '<A HREF="' in lines[i1]:
            continue
        if not '>motion</A>' in lines[i1]:
            continue
        
        index1 = lines[i1].index('<A HREF="')+len('<A HREF="')
        index2 = index1+lines[i1][index1:].index('"')
        url = lines[i1][index1:index2]
        urlID = url[46:]

        if urlID == l_continue[1]:
            l_continue[0] = False
            continue
        elif l_continue[0] == True:
            continue

        pdb1, chain1, pdb2, chain2 = parse_pdbIDs_from_url(url)
        if pdb1 == '' and pdb2 == '':
            continue
        elif pdb1 == 'new1' and pdb2 == 'new2':
            continue
        elif pdb1 == None and pdb2 == None:
            continue
        if len(pdb1) != 4: ## tmp!!!
            print url
            print line
            print pdb1
            stop
        if len(pdb2) != 4: ## tmp!!!
            print url
            print line
            print pdb2
            stop
        l_pdbs += [[urlID,pdb1,chain1,pdb2,chain2,]]
        fd = open('molmovdb.txt','a')
        fd.write('%s\t%s\t%s\t%s\t%s\n' %(urlID,pdb1,chain1,pdb2,chain2,))
        fd.close()

    return l_pdbs


def parse_pdbIDs_from_url(url):

    print url
    try:
        urllines = urllib2.urlopen(url)
    except:
        pdb1 = chain1 = pdb2 = chain2 = None
        urllines = None

    if urllines:
        lines = urllines.readlines()
        for i1 in range(len(lines)):
            if 'Best representative' in lines[i1]:
                line1 = lines[i1+11]
                line2 = lines[i1+12]
                if 'upload' in line1 or 'upload' in line2:
                    pdb1 = chain1 = pdb2 = chain2 = None
                else:
                    pdb1, chain1 = parse_pdbID_from_line(line1)
                    pdb2, chain2 = parse_pdbID_from_line(line2)

    return pdb1, chain1, pdb2, chain2


def parse_pdbID_from_line(line):

    index2 = line.index('</A>')
    index1 = line[:index2].rindex('>')+1
    pdb = line[index1:index2].lower()
    index = index2+4+3
    chain = line[index]

    return pdb, chain


if __name__ == '__main__':
    main()

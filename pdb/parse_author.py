## built-in
import os, sys
import parse_mmCIF

def main():

    fd = open('db_authors.txt','r')
    lines = fd.readlines()
    fd.close()

    d_authors = {}
    for line in lines:
        l = line.strip().split()
        pdb = l[0]
        s_authors = l[1:]
        d_authors[pdb] = s_authors

    lines_out = []

    path = '/media/WDMyBook1TB/2TB/mmCIF'
    l_dns = os.listdir(path)
    l_dns.sort()
    for i in range(len(l_dns)):
        dn = l_dns[i]

        if dn < sys.argv[-1]:
            continue

        if not os.path.isdir('%s/%s' %(path,dn)):
            continue

        print '%s/%s %s' %(i+1,len(l_dns), dn)
        l_fns = os.listdir('%s/%s' %(path,dn))
        l_fns.sort()
        for fn in l_fns:
            if fn[-3:] == '.gz':
                continue

            pdb = fn[0:4]

            if pdb in d_authors.keys():
                continue

            print pdb

            fd = open('%s/%s/%s' %(path, dn, fn), 'r')
            lines = fd.readlines()
            fd.close()

            d_mmCIF = parse_mmCIF.main(
                pdb,lines,
                l_data_categories = [
                    '_audit_author',
                    '_citation_author',
                    ], ## parse selected data categories
                l_data_categories_break = [
                    '_citation_author',
                    ],
                )

            l_authors = d_mmCIF['_audit_author.name']
            s_authors = ';'.join(l_authors)

            if d_mmCIF['_audit_author.name'] == []:
                print d_mmCIF['_citation_author.name']
                print d_mmCIF['_audit_author.name']
                stop

            line = '%s %s\n' %(pdb,s_authors,)
            lines_out += [line]

            fd = open('db_authors.txt','a')
            fd.write(line)
            fd.close()

            d_authors[pdb] = s_authors

    ##
    ## write to file
    ##
    lines_out = []
    for pdb,s_authors in d_authors.items():
        line = '%s %s\n' %(pdb,s_authors,)
        lines_out += [line]
    fd = open('db_authors.txt','w')
    fd.writelines(lines_out)
    fd.close()

    return

if __name__ == '__main__':
    main()

'''this script checks for continuous numbering of residues. really not that important...'''

import os,sys

def main():

    pdbpath = '/media/WDMyBook1TB/2TB/pdb/'
    pdbs = []

    l_dn = os.listdir(pdbpath)
    l_dn.sort()
    for dn in l_dn:
        if dn < sys.argv[-1]:
            continue
        print dn
        l_fn = os.listdir(pdbpath+dn)
        for fn in l_fn:
            if fn[-3:] == '.gz':
                continue
            pdb = fn[3:7]

            fd = open('%s%s/pdb%s.ent' %(pdbpath,pdb[1:3],pdb),'r')
            lines = fd.readlines()
            fd.close()

##            print k, pdb, lines[0][:-1]

            prevchain = 'N/A'
            prevresno = 'N/A'
            prevline = ''

            error = False

            for i in range(len(lines)):

                line = lines[i]

                record = line[:6].strip()

                if record == 'ATOM':

                    chain = line[21]
                    res_no = int(line[22:26])
                    res_name = line[17:20].strip()
                    if res_name in ['C','A','U','G','I','DC','DA','DT','DG','DI','N']:
                        continue

                    prevline, prevchain, prevresno, error = check(error,prevchain,chain,prevresno,res_no,i,lines,prevline)

    ##            if record == 'REMARK':
    ##
    ##                remark = int(line[7:10])
    ##
    ##                if remark == 465 and line[10:80].strip() not in [
    ##                    '','MISSING RESIDUES',
    ##                    'THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE',
    ##                    'EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN',
    ##                    'IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)',
    ##                    'M RES C SSSEQI',
    ##                    ]:
    ##
    ##                    chain = line[19]
    ##                    print line
    ##                    res_no = int(line[22:26])
    ##
    ##                    prevline, prevchain, prevresno, error = check(error,prevchain,chain,prevresno,res_no,i,lines, prevline)

                elif record == 'MODEL':

                    if prevchain != 'N/A':
                        break

            if error == True:
                print pdb
                stop

    return


def check(error,prevchain,chain,prevresno,res_no,i,lines,prevline):

    if prevchain == 'N/A' or chain != prevchain:

        prevchain = chain
        prevresno = res_no

    elif chain == prevchain:

        if res_no < prevresno and prevline[26] != 'P':

            print prevline[:-1]
            print lines[i][:-1]
            prevresno = res_no
            error = True

        else:

            prevresno = res_no

    prevline = lines[i]

    return prevline, prevchain, prevresno, error


if __name__ == '__main__':
    main()

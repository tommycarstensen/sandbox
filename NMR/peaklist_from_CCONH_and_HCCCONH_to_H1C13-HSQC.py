def main():

    fd = open('ACBP_HCCCONH_NHH.list', 'r')
    linesH = fd.readlines()[2:]
    fd.close()

    fd = open('ACBP_CCONH_NCH.list', 'r')
    linesC = fd.readlines()[2:]
    fd.close()

    d_ass = {}
    ## check if HCCCONH peaks are in CCONH peak list
    d_ass = identify_non_assignments(linesH,linesC,d_ass)
    ## check if CCONH peaks are in HCCCONH peak list
    d_ass = identify_non_assignments(linesC,linesH,d_ass)

    ## add proton and carbon chemical shifts to HSQC peak list
    write_2Dpeaklist(d_ass,'ACBP_HSQC_C13.list')

    return
    

def write_2Dpeaklist(d_ass,file_output):

    lines = ['      Assignment         w1         w2  \n\n']
    for ass in d_ass.keys():
        cs = d_ass[ass]
        line = '%17s%11.3f%11.3f\n' %(ass, float(cs[0]), float(cs[1]))
        lines += [line]

    lines += [
        '        M24?CE-HE     16.727      1.698\n',
        '        M46?CE-HE     16.187      2.037\n',
        '        M70?CE-HE     16.963      2.097\n',
        '       I86CG2-HG2     17.249      0.839\n',
        ]

    fd = open(file_output,'w')
    fd.writelines(lines)
    fd.close()
    return


def identify_non_assignments(lines1,lines2,d_ass):

    double_ass = []
    # loop over lines of the 1st peak list
    for line1 in lines1:
        found = 0
        # split line in assignment and chemical shifts
        l1 = line1.split()
        # split assignment in atoms
        ass1 = l1[0].split('-')
        # continue if peak not assigned
        if ass1[0] == '?':
            continue
        # find corresponding assignment in 2nd peak list
        for line2 in lines2:
            l2 = line2.split()
            ass2 = l2[0].split('-')
            if ass1[0] == ass2[0] and ass1[2][:ass1[2].index('H')+1] == ass2[2][:ass2[2].index('H')+1]:
##            if ass1[0] == ass2[0] and ass1[2] == ass2[2]:
	        if 'H' in ass1[1][1:] and 'C' in ass2[1][1:]:
                    assH = ass1; assC = ass2
	            indexH = assH[1][1:].index('H')
                    atomH = assH[1][indexH+2:].replace('a','').replace('b','')
	            indexC = assC[1][1:].index('C')
                    atomC = assC[1][indexC+2:].replace('a','').replace('b','')
                    ass = '%s%s%s' %(assC[1],'-',assH[1][indexH+1:])
                    cs = [l2[2],l1[2]]
                if 'H' in ass2[1][1:] and 'C' in ass1[1][1:]:
                    assH = ass2; assC = ass1
	            indexH = assH[1][1:].index('H')
                    atomH = assH[1][indexH+2:].replace('a','').replace('b','')
	            indexC = assC[1][1:].index('C')
                    atomC = assC[1][indexC+2:].replace('a','').replace('b','')
                    ass = '%s%s%s' %(assC[1],'-',assH[1][indexH+1:])
                    cs = [l1[2],l2[2]]
                if len(atomH) == 0 or len(atomC) == 0 or atomH[0].upper() != atomH[0] or atomC[0].upper() != atomC[0]:
                    print 'H', assH, 'C', assC, 'ass', ass
                    stop
                if atomH == atomC:
                    if ass in d_ass.keys() and d_ass[ass] != cs:
                        double_ass += [ass]
                    else:
                        d_ass[ass] = cs
                        found = 1
                        break
        if found == 0:
            if line1[line1.index('-')+1] != 'K':
                print line1[:-1]

    print double_ass

    return d_ass


if __name__ == '__main__':
    main()

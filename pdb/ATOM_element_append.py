#!/bin/env python
#
# $Id$
#
# Tommy Carstensen, University College Dublin
#

def main():

    pdb = 'lysmgm.pdb'

    fd = open(pdb, 'r')
    lines = fd.readlines()
    fd.close()

    for i in range(len(lines)):
        line = lines[i]
        if line[:4] == 'ATOM':
            atom_name = line[12:16].strip()
            element = atom_name[0]
            line = '%s          %2s%s\n' %(line[:66], element, line[78:80],)
            lines[i] = line

    fd = open('element'+pdb, 'w')
    fd.writelines(lines)
    fd.close()
    
if __name__== '__main__':
    main()
		    

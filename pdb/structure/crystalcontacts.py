#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2009

## script for generating crystal contacts of biological units
## very slow compared to WHATIF

import os, numpy, math
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import biounit,parse_pdb
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import combinatorics

d_radii_vdw = {
    'H':1.20,
    'C':1.70,
    'N':1.55,
    'O':1.52,
    'S':1.80,
    }

s_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

l_translations = combinatorics.permutation_w_rep([-1,0,1,],3)
l_translations.remove([0,0,0,])

def main():

    for pdb in [
    ##    '2hhb',

    ##    '1hho',
    ##    '1hv4',
        
        '2lzm',
##        '2lzt',

##        '1c76',
##        '1c77',
##        '1c78',
        ]:

    ##    os.system('cp /local/data/pdb/%s/pdb%s.ent %s.pdb' %(pdb[1:3],pdb,pdb,))
    ##
    ####    biounit.biounit().main(pdb, '/data/remediated_pdb/', exclude_ligands = True)
    ####
    ####    os.system('cp %s_1.pdb %s.pdb' %(pdb,pdb,))
    ##
    ##    fd = open('%s.pdb' %(pdb),'r')
    ##    lines = fd.readlines()
    ##    fd.close()

        fd = open('C:\Users\Tommy Carstensen\pdb\%s.pdb' %(pdb),'r')
        lines = fd.readlines()
        fd.close()

        d_header = parse_pdb.parse_header(lines)
        d_coordinates, d_ATOMseq = parse_pdb.parse_coordinates(lines,d_header,)
        l_coordinates = []

        a = d_header['CRYST1']['edges'][0]
        b = d_header['CRYST1']['edges'][1]
        c = d_header['CRYST1']['edges'][2]
        alpha = math.pi*d_header['CRYST1']['angles'][0]/180.
        beta = math.pi*d_header['CRYST1']['angles'][1]/180.
        gamma = math.pi*d_header['CRYST1']['angles'][2]/180.
        ## unit cell voumne
        volume = a*b*c*math.sqrt(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
        matrix_fractional2cartesian = numpy.array([
            [a, b*math.cos(gamma), c*math.cos(beta),],
            [0, b*math.sin(gamma), c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma),],
            [0,0,volume/(a*b*math.sin(gamma)),],
            ])

    ##    lines = []
        for symop in d_header['REMARK290'].keys():
            print 'symop', symop
            matrix_symop = d_header['REMARK290'][symop]['matrix']
            vector_symop = d_header['REMARK290'][symop]['vector']
            for i_translation in range(len(l_translations)):
                print 'symop', symop, 'translation', i_translation
                vector_translation = l_translations[i_translation]
                
                for chain2 in d_coordinates['chains'].keys():
                    for res_no2 in d_coordinates['chains'][chain2]['residues'].keys():
    ##                    print symop, i, chain2, res_no2
                        for iCode2 in d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'].keys():
                            for chain1 in d_coordinates['chains'].keys():
                                for res_no1 in d_coordinates['chains'][chain1]['residues'].keys():
    ##                                if res_no2 != 75:
    ##                                    continue
                                    for iCode1 in d_coordinates['chains'][chain1]['residues'][res_no1]['d_iCodes'].keys():
                                        for atom_name1 in d_coordinates['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'].keys():
                                            coordinate1 = d_coordinates['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'][atom_name1]['coordinate']
                                            for atom_name2 in d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'].keys():
                                                coordinate2 = d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name2]['coordinate']
                                                coordinate2 = numpy.dot(matrix_symop,coordinate2)+vector_symop
                                                coordinate2 += numpy.dot(matrix_fractional2cartesian,vector_translation)
                                                coordinate2[0] = round(coordinate2[0],3)
                                                coordinate2[1] = round(coordinate2[1],3)
                                                coordinate2[2] = round(coordinate2[2],3)

                                                vicinity = False
                                                distant = False
                                                dist = math.sqrt(sum((coordinate2-coordinate1)**2))
                                                dist = 0
    ##                                            if dist < 5:
    ##                                                print '%2i %.2f %.2f %4i %4s %s %4i %4s %s ' %(
    ##                                                    i, round(dist,2), dist_treshold,
    ##                                                    res_no2, atom_name2, coordinate2, res_no1, atom_name1, coordinate1
    ##                                                    )
                                                dist_treshold = d_radii_vdw[atom_name1[0]]+d_radii_vdw[atom_name2[0]]+0.25
                                                ## break atom_name2 loop
                                                if dist < dist_treshold:
                                                    vicinity = True
                                                    break
                                                ## break atom_name2 loop
                                                if dist > 10.: ## length of lysine is 7-8Angstrom
                                                    distant = True
                                                    break
                                            ## break atom_name1 loop (append line and check next iCode1)
                                            if vicinity == True:
                                                for atom_name3 in d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'].keys():
                                                    coordinate3 = d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name3]['coordinate']
                                                    coordinate3 = numpy.dot(matrix_symop,coordinate3,)+vector_symop
                                                    coordinate3 += numpy.dot(matrix_fractional2cartesian,vector_translation,)
                                                    line = build_line(
                                                        atom_name3,d_coordinates,coordinate3,
                                                        chain2,res_no2,iCode2,
                                                        symop, i_translation,
                                                        )
                                                    lines += [line]
                                                break
                                            ## break atom_name1 loop (check next iCode1)
                                            if distant == True:
                                                break
                                        ## break iCode1 loop (check next iCode2)
                                        if vicinity == True:
                                            break
                                    ## break resno1 loop (check next iCode2)
                                    if vicinity == True:
                                        break
                                ## break chain1 loop (check next iCode2)
                                if vicinity == True:
                                    break

        fd = open('%s_crystalcontacts.pdb' %(pdb),'w')
        fd.writelines(lines)
        fd.close()

        source = 'whatif.src'
        fd = open(source,'w')
        fd.writelines([
            '/software/whatif/DO_WHATIF.COM <<EOF\n',
            'GETMOL %s.pdb\n' %(pdb,),
            '%s\n' %(pdb,),
            '%DELWAT\n',
            '%DELLIG\n',
            '%SOUSHL\n', ## crystal contacts
            '%MAKMOL\n',
            '\n', ## The file header will be copied from a PDB file. Hit return for the default header that has no information in it.
            'soushl_%s.pdb\n' %(pdb),
            'TOT 0\n',
            '\n', ## REMARKS
            'STOP\n',
            'Y\n',
            ])
        fd.close()
        os.system('source %s > whatif_surface/%s.out' %(source, pdb))

    os.system('rm DRG* DAVADRUG.PDB ALTERR.LOG PDBFILE whatif.src FOR*.DAT WHATIF.FIG')

    return


def build_line(
    atom_name3,d_coordinates,coordinate3,
    chain2,res_no2,iCode2,
    symop, i_translation,
    ):

    element = d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name3]['element']
    atom_name = '%2s%2s' %(element.rjust(2), atom_name3[atom_name3.index(element)+1:].ljust(2))

    record = 'ATOM'
    atom_no = 1
    altloc = ' '
    res_name = d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['res_name']
    chain = s_alphabet[i_translation+1]
    chain = s_alphabet[symop]
    res_no = res_no2
    iCode = iCode2
    x = coordinate3[0]
    y = coordinate3[1]
    z = coordinate3[2]
    occupancy = 100.
    tempfactor = 10*symop
    charge = 0.
    line = '%6s%5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' %(
        record.ljust(6),atom_no,atom_name,altloc,res_name,
        chain,res_no,iCode,
        x,y,z,occupancy,tempfactor,element.rjust(2),charge,
        )

    return line


if __name__ == '__main__':
    main()

#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2009

## script for generating crystal contacts of biological units

import os, numpy, math
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import biounit,parse_pdb
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import combinatorics

path_pdb = '/data/remediated_pdb'

pdb = sys.argv[sys.argv.index('-pdb')+1]

d_radii_vdw = {
    'H':1.20,
    'C':1.70,
    'N':1.55,
    'O':1.52,
    'S':1.80,
    }

l_translations = combinatorics.permutation_w_rep([-1,0,1,],3)
l_translations.remove([0,0,0,])

s_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' # 0123456789

## Give the cutoff for "NEAR" contacts (Angstrom)
if '-sympar' in sys.argv:
    sympar = float(sys.argv[sys.argv.index('-sympar')+1])
else:
    sympar = 5.

def main(pdb):

    os.system('cp %s/%s/pdb%s.ent %s.pdb' %(path_pdb,pdb[1:3],pdb,pdb,))

    ##
    ## Tommy crystal contacts
    ##

    ## create biounit
    biounit.biounit().main(pdb, '/data/remediated_pdb/', exclude_ligands = True)

    ##
    ## parse header (just use the asu instead!!!)
    ##
    fd = open('%s/%s/pdb%s.ent' %(path_pdb,pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()
    for i in range(len(lines)):
        line = lines[i]
        record = line[:6].strip()
        if record in ['MODEL','ATOM',]:
            break
    lines_header = lines[:i]

    ##
    ## parse coordinates
    ##
    fd = open('%s_1.pdb' %(pdb),'r')
    lines_biounit = fd.readlines()
    fd.close()

    fd = open('%s.pdb' %(pdb),'r')
    lines_asu = fd.readlines()
    fd.close()

##        fd = open('C:\Users\Tommy Carstensen\pdb\%s.pdb' %(pdb),'r')
##        lines = fd.readlines()
##        fd.close()

    d_header = parse_pdb.parse_header(lines_header)
    d_coordinates_biounit, d_ATOMseq = parse_pdb.parse_coordinates(
        lines_biounit,d_header,
        parse_atom_seq = False, parse_ligands = False,
        )
    d_coordinates_asu, d_ATOMseq = parse_pdb.parse_coordinates(
        lines_asu,d_header,
        parse_atom_seq = False, parse_ligands = False,
        )

    ## set new chain IDs
    l_old_chains = d_coordinates_asu['chains'].keys()
    l_new_chains = list(
        set(s_alphabet)-set(l_old_chains)
        )
    l_old_chains.sort()
    l_new_chains.sort()
    d_chains = {}
    for i in range(len(l_old_chains)):
        d_chains[l_old_chains[i]] = l_new_chains[i]

##    a = d_header['CRYST1']['edges'][0]
##    b = d_header['CRYST1']['edges'][1]
##    c = d_header['CRYST1']['edges'][2]
##    alpha = math.pi*d_header['CRYST1']['angles'][0]/180.
##    beta = math.pi*d_header['CRYST1']['angles'][1]/180.
##    gamma = math.pi*d_header['CRYST1']['angles'][2]/180.
##    ## unit cell voumne
##    volume = a*b*c*math.sqrt(1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2+2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
##    matrix_fractional2cartesian = numpy.array([
##        [a, b*math.cos(gamma), c*math.cos(beta),],
##        [0, b*math.sin(gamma), c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma),],
##        [0,0,volume/(a*b*math.sin(gamma)),],
##        ])
    matrix_scale = d_header['SCALE']
    matrix_scalei = numpy.linalg.inv(matrix_scale)

##    lines = []
    l_symop = d_header['REMARK290'].keys()
    l_symop.sort()
##        l_symop = [1]
    for symop in l_symop:
##        if symop != 1:
##            continue
        matrix_symop = d_header['REMARK290'][symop]['4x4matrix']
        for i in range(len(l_translations)):
            vector_translation = l_translations[i]
##            if vector_translation != [-1,0,0]:
##                continue
            vector_translation = numpy.array([vector_translation[0],vector_translation[1],vector_translation[2],0,])
##                if i != 0:
##                    continue
##                if vector_translation[0] != 1:
##                    continue
##                if vector_translation[1] != -1:
##                    continue
##                if vector_translation[2] != -1:
##                    continue
            for chain2 in d_coordinates_asu['chains'].keys():
                for res_no2 in d_coordinates_asu['chains'][chain2]['residues'].keys():
##                    if res_no2 != 1:
##                        continue
                    print '%4s %2i/%2i %2i/26 %1s %4i' %(pdb, symop, len(l_symop), i+1, chain2, res_no2)
                    for iCode2 in d_coordinates_asu['chains'][chain2]['residues'][res_no2]['d_iCodes'].keys():
                        for chain1 in d_coordinates_biounit['chains'].keys():
                            for res_no1 in d_coordinates_biounit['chains'][chain1]['residues'].keys():
                                for iCode1 in d_coordinates_biounit['chains'][chain1]['residues'][res_no1]['d_iCodes'].keys():
                                    for atom_name1 in d_coordinates_biounit['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'].keys():
                                        coordinate1 = d_coordinates_biounit['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms'][atom_name1]['coordinate']
                                        for atom_name2 in d_coordinates_asu['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'].keys():
                                            coordinate2 = d_coordinates_asu['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name2]['coordinate']
##                                            coordinate2 = numpy.dot(matrix_symop,coordinate2)+vector_symop
##                                            coordinate2 += numpy.dot(matrix_fractional2cartesian,vector_translation)
##                                            coordinate2[0] = round(coordinate2[0],3)
##                                            coordinate2[1] = round(coordinate2[1],3)
##                                            coordinate2[2] = round(coordinate2[2],3)

                                            ## conversion from 3x1 vector to 4x1 vector
                                            coordinate2 = numpy.array([coordinate2[0],coordinate2[1],coordinate2[2],1.,])
                                            ## conversion from cartesian to fractional coordinates
                                            coordinate2 = numpy.dot(matrix_scale,coordinate2)
                                            ## symmetry operator (before unit cell translation???)
                                            coordinate2 = numpy.dot(matrix_symop,coordinate2)
                                            ## unit cell translation
                                            coordinate2 += vector_translation
                                            ## conversion from fractional to cartesian coordinates
                                            coordinate2 = numpy.dot(matrix_scalei,coordinate2)
                                            ## conversion from 4x1 vector to 3x1 vector
                                            coordinate2 = numpy.array([coordinate2[0],coordinate2[1],coordinate2[2],])

                                            vicinity = False
                                            distant = False
                                            dist = math.sqrt(sum((coordinate2-coordinate1)**2))
                                            if dist < 5:
                                                print vector_translation
                                                print numpy.dot(matrix_fractional2cartesian,vector_translation)
                                                print matrix_fractional2cartesian
                                                print '%2i %2i %.2f %.2f %1s %4i %4s %s %1s %4i %4s %s %s' %(
                                                    symop, i, round(dist,2), dist_treshold,
                                                    chain2, res_no2, atom_name2, coordinate2,
                                                    chain1, res_no1, atom_name1, coordinate1,
                                                    d_coordinates_asu['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name2]['coordinate'],
                                                    )

                                            ## break atom_name2 loop
                                            if dist > 80.:
                                                distant = True
                                                break

                                            dist = 0 ## tmp!!! temp!!! get *all* translations

                                            dist_treshold = d_radii_vdw[atom_name1[0]]+d_radii_vdw[atom_name2[0]]+.25
                                            ## break atom_name2 loop
                                            if dist < dist_treshold:
                                                vicinity = True
                                                break

                                        ## break atom_name1 loop (append line and check next iCode1)
                                        if vicinity == True:
                                            for atom_name3 in d_coordinates_asu['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'].keys():
                                                coordinate3 = d_coordinates_asu['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name3]['coordinate']
                                                ## conversion from 3x1 vector to 4x1 vector
                                                coordinate3 = numpy.array([coordinate3[0],coordinate3[1],coordinate3[2],1.,])
                                                ## conversion from cartesian to fractional coordinates
                                                coordinate3 = numpy.dot(matrix_scale,coordinate3)
                                                ## symmetry operator (before unit cell translation???)
                                                coordinate3 = numpy.dot(matrix_symop,coordinate3)
                                                ## unit cell translation
                                                coordinate3 += vector_translation
                                                ## conversion from fractional to cartesian coordinates
                                                coordinate3 = numpy.dot(matrix_scalei,coordinate3)
                                                ## conversion from 4x1 vector to 3x1 vector
                                                coordinate3 = numpy.array([coordinate3[0],coordinate3[1],coordinate3[2],])
                                                line = build_line(
                                                    atom_name3,d_coordinates_asu,coordinate3,
                                                    chain2,res_no2,iCode2,
                                                    d_chains,
                                                    )
                                                lines += [line]
                                            break
                                        ## break atom_name1 loop (check next iCode1)
                                        if distant == True:
                                            break
                                    ## break iCode1 loop (check next iCode2)
                                    if vicinity == True:
                                        break
                                    if distant == True:
                                        break
                                ## break resno1 loop (check next iCode2)
                                if vicinity == True:
                                    break
                                if distant == True:
                                    break
                            ## break chain1 loop (check next iCode2)
                            if vicinity == True:
                                break
                            if distant == True:
                                break

    fd = open('%s_biou_cc.pdb' %(pdb),'w')
    fd.writelines(lines)
    fd.close()

    return


def build_line(
    atom_name3,d_coordinates,coordinate3,
    chain2,res_no2,iCode2,
    d_chains,
    ):

    element = d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms'][atom_name3]['element']
    atom_name = '%2s%2s' %(element.rjust(2), atom_name3[atom_name3.index(element)+1:].ljust(2))

    record = 'ATOM'
    atom_no = 1
    altloc = ' '
    res_name = d_coordinates['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2]['altlocs'][' ']['res_name']
    chain = d_chains[chain2]
    res_no = res_no2
    iCode = iCode2
    x = coordinate3[0]
    y = coordinate3[1]
    z = coordinate3[2]
    occupancy = 100.
    tempfactor = 100.
    charge = 0
    line = '%6s%5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n' %(
        record.ljust(6),atom_no,atom_name,altloc,res_name,
        chain,res_no,iCode,
        x,y,z,occupancy,tempfactor,element.rjust(2),charge,
        )

    return line


if __name__ == '__main__':

        main(pdb)

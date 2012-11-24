def main():

    pdb = '2lzt'

    import Numeric
    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import combinatorics

    l_translations = combinatorics.permutation_w_rep([-1,0,1,],3)
    l_translations.remove([0,0,0,])

    d_coordinates,d_290,l_coordinates,cryst1 = parse_pdb(pdb)

    l_chains = [
        'B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q',
        'R','S','T','U','V','W','X','Y','Z','b','c','d','e','f','g','h',
        'i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x',
        'y','z','0','1','2','3','4','5','6','7','8','9',
        ]
    chain_index = -1
    atom_no = 1001
    prev_residue = ['N/A','N/A','N/A',]
    lines = []
    for operator in d_290.keys():
        matrix = d_290[operator]
        for i in range(len(l_translations)):
            translation = Numeric.array(l_translations[i])*cryst1
            print operator,i,chain_index
            for chain in d_coordinates.keys():
                l_res_nos = d_coordinates[chain].keys()
                l_res_nos.sort()
                for res_no in l_res_nos:
                    for atom_name in d_coordinates[chain][res_no]['atoms'].keys():
                        coordinate2 = d_coordinates[chain][res_no]['atoms'][atom_name]['coordinate']
                        coordinate2 = Numeric.array(list(coordinate2))
                        coordinate2 += translation
                        vicinal = False
                        for coordinate1 in l_coordinates:
                            sqdist = sum((coordinate2-coordinate1)**2)
                            if sqdist < 225:
                                vicinal = True
                                break
                        if vicinal == True:
                            break
                    if vicinal == True or prev_residue == [i,res_no-1,]:

                        ## new chain ID
                        if prev_residue not in [
                            [i,res_no-3,],
                            [i,res_no-2,],
                            [i,res_no-1,],
                            ]:
                            print i, res_no, prev_residue
                            chain_index += 1

                        ## add previous residue
                        if prev_residue not in [
                            [i,res_no-2,],
                            [i,res_no-1,],
                            ] and res_no-1 in d_coordinates[chain].keys():
                            lines,atom_no = newline(atom_no,translation,lines,l_chains,chain_index,d_coordinates,pdb,chain,res_no-1,)

                        ## add current residue
                        if prev_residue not in [
                            [i,res_no-1,],
                            ]:
                            lines,atom_no = newline(atom_no,translation,lines,l_chains,chain_index,d_coordinates,pdb,chain,res_no,)

                        ## add next residue
                        if res_no+1 in d_coordinates[chain].keys():
                            lines,atom_no = newline(atom_no,translation,lines,l_chains,chain_index,d_coordinates,pdb,chain,res_no+1,)

                        if vicinal == True:
                            prev_residue = [i,res_no,]

    fd = open('%s_contacts.pdb' %(pdb),'w')
    fd.writelines(lines)
    fd.close()

    print len(lines), len(set(lines))
                                
    return


def newline(atom_no,translation,lines,l_chains,chain_index,d_coordinates,pdb,chain,res_no,):

    import Numeric

    res_name = d_coordinates[chain][res_no]['res_name']

    for atom_name in d_coordinates[chain][res_no]['atoms'].keys():

        coordinate = d_coordinates[chain][res_no]['atoms'][atom_name]['coordinate']
        coordinate = Numeric.array(list(coordinate))
        coordinate += translation
        element = d_coordinates[chain][res_no]['atoms'][atom_name]['element']

        atom_name = '%2s%2s' %(element.rjust(2), atom_name[atom_name.index(element)+1:].ljust(2))

        atom_no += 1

        record = 'ATOM'
        altloc = ' '
        iCode = ' '
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        occupancy = 1
        tempfactor = 1
        charge = ' '
        line = '%6s%5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' %(
            record.ljust(6),atom_no,atom_name,altloc,res_name,
            l_chains[chain_index],res_no,iCode,
            x,y,z,occupancy,tempfactor,element.rjust(2),charge,
            )
        lines += [line]

    return lines, atom_no


def parse_pdb(pdb):

    import Numeric

    fd = open('/data/remediated_pdb/%s/pdb%s.ent' %(pdb[1:3],pdb),'r')
    lines = fd.readlines()
    fd.close()

    l_modres = []
    d_coordinates = {}
    l_coordinates = []
    d_290 = {}

    for line in lines:

        record = line[:6].strip()

        if record == 'REMARK':

            remark = int(line[6:10])

            if remark == 290:
                if line[13:18] == 'SMTRY':

                    smtry = int(line[18])
                    operator_number = int(line[22])
                    if smtry == 1:
                        d_290[operator_number] = Numeric.zeros((3,4),typecode='d')

                    rx = float(line[24:33])
                    ry = float(line[34:43])
                    rz = float(line[44:53])
                    t = float(line[60:68])
                    d_290[operator_number][smtry-1][0] = rx
                    d_290[operator_number][smtry-1][1] = ry
                    d_290[operator_number][smtry-1][2] = rz
                    d_290[operator_number][smtry-1][3] = t

        elif record == 'CRYST1':
            a = float(line.split()[1])
            b = float(line.split()[2])
            c = float(line.split()[3])
            cryst1 = Numeric.array([a,b,c,])

        elif record == 'MODRES':

            hetID = line[12:15].strip()
            chain = line[16]
            res_no = int(line[18:22])
            iCode = line[22]

            modres = '%s:%s:%s' %(chain,res_no,iCode)

            l_modres += [modres]

        elif record == 'ATOM':

            d_coordinates,l_coordinates = parse_record_ATOM(d_coordinates,line,record,l_modres,l_coordinates)

        elif record == 'HETATM':

            d_coordinates,l_coordinates = parse_record_ATOM(d_coordinates,line,record,l_modres,l_coordinates)

    return d_coordinates, d_290, l_coordinates, cryst1


def parse_record_ATOM(d_coordinates,line,record,l_modres,l_coordinates,):

    import Numeric

    '''this function ignores altlocs'''

    atom_no = int(line[6:11])
    atom_name = line[12:16].strip()
    altloc = line[16]
    res_name = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    coordinate = Numeric.array([x, y, z])
    element = line[76:78].strip()

    if altloc != ' ':
        print line
        stop_altloc
    if iCode != ' ':
        print line
        stop_iCode

    modres = '%s:%s:%s' %(chain,res_no,iCode)

    if record == 'ATOM' or (record == 'HETATM' and modres in l_modres):

        if not chain in d_coordinates.keys():
            d_coordinates[chain] = {}
        if not res_no in d_coordinates[chain].keys():
            d_coordinates[chain][res_no] = {}
        if not 'atoms' in d_coordinates[chain][res_no].keys():
            d_coordinates[chain][res_no]['atoms'] = {}
        d_coordinates[chain][res_no]['atoms'][atom_name] = {'coordinate':coordinate,'atom_no':atom_no,'element':element}
        d_coordinates[chain][res_no]['res_name'] = res_name

        l_coordinates += [Numeric.array([x, y, z])]
        
    return d_coordinates,l_coordinates


if __name__=='__main__':
    main()

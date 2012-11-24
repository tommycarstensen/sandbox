d_pep = {
    'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
    'UNK':'X','ASX':'X','GLX':'X',
    'MSE':'M', ## ligand in 2e1a
    }

l_nuc = [
     'A', 'C', 'G', 'U',      'I', ## ribonucleotides
    'DA','DC','DG',     'DT','DI', ## deoxyribonucleotides
    'N', ## N is any 5'-monophosphate nucleotide
    ]

s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'


def check_if_SEQRESres(res_name,record,d_header,chain,res_no,iCode):

    if res_name in d_pep.keys()+l_nuc and record in ['ATOM','REMARK465',]:
        return True

    MODRES = False
    if 'MODRES' in d_header.keys():
        if chain in d_header['MODRES'].keys():
            if res_no in d_header['MODRES'][chain].keys():
                if iCode in d_header['MODRES'][chain][res_no].keys():
                    if res_name in d_header['MODRES'][chain][res_no][iCode]:
                        MODRES = True
                    else:
                        print chain, res_no,iCode
                        print res_name,d_header['MODRES'][chain][res_no][iCode]
                        stop

    return MODRES


def append_sequence(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,append_REMARK465,append_ATOM):

    if append_REMARK465 == True:
        for res_no_REMARK465 in l_REMARK465_res_nos:
            if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                for iCode_REMARK465 in l_iCodes_REMARK465:
                    res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                    d_ATOMseq[chain]['seq'] += [res_name_REMARK465]
                    d_ATOMseq[chain]['res_nos'] += [res_no_REMARK465]
                    d_ATOMseq[chain]['iCodes'] += [iCode_REMARK465]
                    d_ATOMseq[chain]['altlocs'] += [' ']
                    d_ATOMseq[chain]['records'] += ['REMARK465']
                del d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]
##            else: ## not 3b95
##                break

    if append_ATOM == True:
        d_ATOMseq[chain]['seq'] += [res_name_ATOM]
        d_ATOMseq[chain]['res_nos'] += [res_no]
        d_ATOMseq[chain]['iCodes'] += [iCode]
        d_ATOMseq[chain]['altlocs'] += [altloc]
        d_ATOMseq[chain]['records'] += [record]

    return d_ATOMseq, d_header


def parse_header(lines):

    d_header = {}

    for i in range(len(lines)):
        line = lines[i]

        record = line[:6].strip()

        if record == 'EXPDTA': ## section 2
            methods = line[10:].strip().split(',')[0]
            if methods[:3] == 'NMR':
                methods = 'NMR'
            elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                methods = 'X-RAY'
            d_header['EXPDTA'] = methods

        elif record == 'REMARK': ## section 2
            d_header = parse_recordREMARK(d_header, line, i, lines)

        elif record == 'SEQRES': ## section 3
            d_header = parse_recordSEQRES(line, d_header)

        elif record == 'MODRES': ## section 3
            parse_recordMODRES(line, d_header,)

        elif record == 'HET': ## section 4
            d_header = parse_recordHET(line, d_header,)

        elif record == 'HELIX': ## section 5
            d_header = parse_recordHELIX(line,d_header,)

        elif record == 'SHEET': ## section 5
            d_header = parse_recordSHEET(line,d_header,)

        elif record == 'TURN': ## section 5
            stop_notobsoletewithv32

        elif record == 'CRYST1':
            d_header = parse_recordCRYST1(line,d_header,)

        elif record in ['SCALE1','SCALE2','SCALE3',]:
            d_header = parse_recordSCALE(line,d_header,)

    return d_header


def parse_recordSCALE(line, d_header,):

    import numpy

    r = int(line[5])
    x = float(line[10:20])
    y = float(line[20:30])
    z = float(line[30:40])
    t = float(line[45:55])
    if r == 1:
        d_header['SCALE'] = numpy.zeros((4,4))
        d_header['SCALE'][3][3] = 1.

    d_header['SCALE'][r-1][0] = x
    d_header['SCALE'][r-1][1] = y
    d_header['SCALE'][r-1][2] = z
    d_header['SCALE'][r-1][3] = t

    return d_header

def parse_recordCRYST1(line, d_header,):

    import numpy

    a = float(line.split()[1])
    b = float(line.split()[2])
    c = float(line.split()[3])
    edges = numpy.array([a,b,c,])
    alpha = float(line.split()[4])
    beta = float(line.split()[5])
    gamma = float(line.split()[6])
    angles = numpy.array([alpha,beta,gamma,])
    space_group = line[55:66].strip()
    d_header['CRYST1'] = {
        'edges':edges,
        'angles':angles,
        'space group':space_group,
        }

    return d_header

def parse_recordSHEET(line, d_header,):

    sheet_ID = line[11:14].strip()
    sheet_no = int(line[7:10])

    chain1 = line[21]
    res_no1 = int(line[22:26])
    iCode1 = line[26]
    chain2 = line[32]
    res_no2 = int(line[33:37])
    iCode2 = line[37]

    if chain1 != chain2:
        print line
        print chain1, chain2
        stop_chainIDs_different

    if 'SHEET' not in d_header.keys():
        d_header['SHEET'] = {}
    if chain1 not in d_header['SHEET'].keys():
        d_header['SHEET'][chain1] = {}
    if not '%4i%1s' %(res_no1,iCode1) in d_header['SHEET'][chain1].keys():
        d_header['SHEET'][chain1]['%4i%1s' %(res_no1,iCode1)] = {'res_no':res_no2,'iCode':iCode2,}

    return d_header


def parse_recordHELIX(line, d_header,):

    chain1 = line[19]
    res_no1 = int(line[21:25])
    iCode1 = line[25]
    chain2 = line[31]
    res_no2 = int(line[33:37])
    iCode2 = line[37]
    try:
        helix_len = int(line[71:76])
    except:
        print line
        print line[71:76]
        stop_helix
    helix_no = int(line[7:10])

    if 'HELIX' not in d_header.keys():
        d_header['HELIX'] = {}
    if chain1 not in d_header['HELIX'].keys():
        d_header['HELIX'][chain1] = {}
    if not '%4i%1s' %(res_no1,iCode1) in d_header['HELIX'][chain1].keys():
        d_header['HELIX'][chain1]['%4i%1s' %(res_no1,iCode1)] = helix_len

    return d_header


def parse_coordinates(
    lines, d_header,
    parse_atom_seq = True, parse_multiple_models = False, parse_ligands = True,
    ):

    model = 999

    d_atomnos = {}
    d_coordinates = {model:{}}
    d_ATOMseq = {}

    for i in range(len(lines)):
        line = lines[i]
        record = line[:6].strip()

        if record == 'ATOM':
            d_coordinates[model], d_line = parse_recordATOM(
                line, d_coordinates[model], lines, i, d_header,
                'ATOM',
                parse_ligands=parse_ligands,
                )
            d_atomnos[d_line['atom_no']] = d_line

        elif record == 'HETATM':
            res_name = line[17:20].strip()
            chain = line[21]
            ## water
            if res_name in ['D2O','H2O',]:
                print pdb, res_name
                stop
            elif res_name in ['HOH','DOD',]: ## DOD e.g. 2d4i
                continue
            ## modified residue of polypeptide or polynucleotide
            elif res_name == 'MSE' and 'MSE' in d_header['SEQRES']['chains'][chain]['seq3']: ## e.g. 1gcj,2e1a
                d_coordinates[model], d_line = parse_recordATOM(
                    line, d_coordinates[model], lines, i, d_header,
                    'ATOM',
                    parse_ligands=parse_ligands,
                    )
            ## (poly)saccharide or other hetero compound
            else:
                d_coordinates[model], d_line = parse_recordATOM(
                    line, d_coordinates[model], lines, i, d_header,
                    'HETATM',
                    parse_ligands=parse_ligands,
                    )
            atom_no = d_line['atom_no']
            d_atomnos[atom_no] = d_line

        elif record == 'MODEL':
            ## WHATIF bug
            if lines[i+1][:6].strip() in ['MODEL','TER','ENDMDL','MASTER',]:
                continue
            model = int(line.split()[1])
            d_coordinates[model] = {}

        elif record == 'ENDMDL':
            continue

    ## 1) append remaining REMARK465 residues
    ## 2) check that sequences are identical
    for chain in d_ATOMseq.keys():

        ## append remaining REMARK465 residues
        if 'REMARK465' in d_header.keys():
            if chain in d_header['REMARK465']['chains'].keys():
                if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                    l_REMARK465_res_nos = d_header['REMARK465']['chains'][chain]['residues'].keys() ## e.g. 1cd0
                    l_REMARK465_res_nos.sort()
                    d_ATOMseq,d_header = append_sequence(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,False,)

        ## check that sequences are identical
        if d_header['SEQRES']['chains'][chain]['seq3'] != d_ATOMseq[chain]['seq']:
            for i in range(len(d_ATOMseq[chain]['seq'])):
                if d_ATOMseq[chain]['seq'][i] != d_header['SEQRES']['chains'][chain]['seq3'][i]:
                    print i
                    break
            l_REMARK465_res_nos = d_header['REMARK465']['chains'][chain]['residues'].keys()
            l_REMARK465_res_nos.sort()
            d_ATOMseq,d_header = append_sequence(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,False,)
            if d_header['SEQRES']['chains'][chain]['seq3'] != d_ATOMseq[chain]['seq']:
                print l_REMARK465_res_nos
                print len(d_header['SEQRES']['chains'][chain]['seq3']), len(d_ATOMseq[chain]['seq'])
                for i in range(len(d_ATOMseq[chain]['seq'])):
                    if d_ATOMseq[chain]['seq'][i] != d_header['SEQRES']['chains'][chain]['seq3'][i]:
                        print i
                        print d_ATOMseq[chain]['seq'][i-1:]
                        print d_header['SEQRES']['chains'][chain]['seq3'][i-1:]
                        break
                print lines[0]
                stop1
            stop2

    if parse_multiple_models == True:
        return d_coordinates, d_ATOMseq
    else:
        model = min(d_coordinates.keys())
        return d_coordinates[model], d_ATOMseq


def parse_recordHET(line, d_header,):

    hetID = line[7:10].strip()
    ## continue if water
    if hetID in ['HOH','DOD',]:
        return d_header
    if hetID in ['D2O','DOD','H2O',]:
        print hetID
        stop
    chain = line[12]
    res_no = int(line[13:17])
    iCode = line[17]
    if 'HET' not in d_header.keys():
        d_header['HET'] = {}
    if chain not in d_header['HET'].keys():
        d_header['HET'][chain] = {}
    if res_no not in d_header['HET'][chain].keys():
        d_header['HET'][chain][res_no] = {}
    if iCode not in d_header['HET'][chain][res_no].keys():
        d_header['HET'][chain][res_no][iCode] = {}
    if hetID not in d_header['HET'][chain][res_no][iCode]:
        d_header['HET'][chain][res_no][iCode][hetID] = {}
    else:
        stop

    MODRES = False
    if 'MODRES' in d_header.keys():
        if chain in d_header['MODRES'].keys():
            if res_no in d_header['MODRES'][chain].keys():
                if iCode in d_header['MODRES'][chain][res_no].keys():
                    if hetID == d_header['MODRES'][chain][res_no][iCode]:
                        MODRES = True
    
    d_header['HET'][chain][res_no][iCode][hetID]['MODRES'] = MODRES

    return d_header

def parse_recordMODRES(line, d_header,):

    hetID = line[12:15].strip()
    chain = line[16]
    res_no = int(line[18:22])
    iCode = line[22]
    res_name = line[24:27].strip()
    txt = line[29:80].strip()

    ## return if e.g. glycosylation site (e.g. 3c43)
    if hetID in set(d_pep.keys())-set(['MSE']) and res_name in set(d_pep.keys())-set(['MSE']):
        return d_header

    if 'MODRES' not in d_header.keys():
        d_header['MODRES'] = {}
    if chain not in d_header['MODRES'].keys():
        d_header['MODRES'][chain] = {}
    if res_no not in d_header['MODRES'][chain].keys():
        d_header['MODRES'][chain][res_no] = {}
    if iCode not in d_header['MODRES'][chain][res_no].keys():
        d_header['MODRES'][chain][res_no][iCode] = hetID
    elif hetID != d_header['MODRES'][chain][res_no][iCode]:
        print chain, res_no, iCode, hetID, d_header['MODRES'][chain][res_no][iCode]
        stop

    return d_header


def parse_recordREMARK(d_header, line, i, lines):

    remark = int(line[6:10])

    if remark == 2:
        if (
            line[11:38] == 'RESOLUTION. NOT APPLICABLE.'
            or
            line[11:38] == 'RESOLUTION. NULL ANGSTROMS.'
            ):
            d_header['REMARK2'] = 'N/A'
            return d_header
        try:
            resolution = float(line[22:27])
            d_header['REMARK2'] = resolution
        except:
            pass

    elif remark == 290: ## symmetry transformations

        d_header = parse_recordREMARK290(line, lines, i, d_header)

    elif remark == 350: ## missing residues

        d_header = parse_recordREMARK350(line, lines, i, d_header)

    elif remark == 465: ## missing residues

        d_header = parse_recordREMARK465(line, lines, i, d_header, remark,)

    elif remark == 470: ## missing atoms

        d_header = parse_recordREMARK470(line, lines, i, d_header, remark,)

    return d_header


def parse_recordREMARK290(line, lines, i, d_header):

##    stop_and_finish_this_function

    import numpy

    if lines[i][15:32] == 'NNNMMM   OPERATOR':
        for j in range(i+1,len(lines)):
            if lines[j][10:].strip() == '':
                break
            matrix = numpy.zeros((4,4))
            matrix[3][3] = 1.
            sym_op_no = int(lines[j][15:18])
            sym_op = lines[j][24:].strip()
            print sym_op
            sym_op = sym_op.split(',')
            l = ['X','Y','Z',]
            for k in range(3):
                ## 
                if sym_op[k] == l[k]:
                    matrix[k][k] = 1.
                ##
                else:
                    ## loop over xyz
                    for s in l:
                        if s in sym_op[k]:
                            ## index xyz
                            index = sym_op[k].index(s)
                            if sym_op[k] == s or index == 0:
                                matrix[k][l.index(s)] = 1.
                            elif sym_op[k][index-1] == '+':
                                matrix[k][l.index(s)] = 1.
                            elif sym_op[k][index-1] == '-':
                                matrix[k][l.index(s)] = -1.
                            else:
                                print
                                print sym_op[k], index
                                print sym_op[k][index-1]
                                print s
                                print lines[j]
                                stop
                            if sym_op[k][0:index-1]+sym_op[k][index+1:] == '':
                                pass
                            elif sym_op[k][index+1:] == '+1/4':
                                matrix[k][3] = .25
                            elif sym_op[k][index+1:] == '+1/2':
                                matrix[k][3] = .5
                            elif sym_op[k][index+1:] == '+3/4':
                                matrix[k][3] = .75
                            else:
                                print
                                print sym_op[k][0:index-1]
                                print sym_op[k][index+1:]
                                print sym_op[k]
                                stop
            if sym_op_no == 1:
                d_header['REMARK290'] = {}
            d_header['REMARK290'][sym_op_no] = {'4x4matrix':matrix}
               
##    if lines[i][13:19] != 'SMTRY1':
##        return d_header
##
##    sym_op = int(lines[i][22])
##    matrix = numpy.zeros((3,3))
##    vector = numpy.array([0.,0.,0.,])
##    for j in range(i,i+3):
##        vx = float(lines[j][24:33])
##        vy = float(lines[j][34:43])
##        vz = float(lines[j][44:53])
##        if lines[j][23] != ' ':
##            stop1
##        if lines[j][33] != ' ':
##            stop1
##        if lines[j][43] != ' ':
##            stop3
##        if lines[j][57] != ' ':
##            stop4
##        t = float(lines[j][58:68])
##        matrix[j-i][0] = vx
##        matrix[j-i][1] = vy
##        matrix[j-i][2] = vz
##        vector[j-i] = t
##
##    if not 'REMARK290' in d_header.keys():
##        d_header['REMARK290'] = {}
##    d_header['REMARK290'][sym_op] = {'matrix':matrix,'vector':vector,}

    return d_header


def parse_recordREMARK350(line, lines, i, d_header):

##    if 'REMARK350' not in d_header.keys():
##        d_header['REMARK350'] = {}
##
##    if lines[i:i+5] != [
##        'REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN           \n',
##        'REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE                \n',
##        'REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS          \n',
##        'REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND                          \n',
##        'REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.                               \n',
##        ]:
##        return d_header
##
##    for j in range(i,len(lines)):
##
##        print lines[j]
##
##        if lines[j][11:23] == 'BIOMOLECULE:':
##
##            biomolecules = biomolecules = line[23:80].replace(' ','').split(',')
##            if len(biomolecules) > 1:
##                stop
##
##            
##
##        if lines[j][:10] != 'REMARK 350':
##            break
##
##    stop
##
##    if d_header['REMARK350'] == {}:
##        stop

##      elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
##                chains = set() ## e.g. 1upp.pdb (vs 8ruc.pdb)
##                line_chains = lines[j][41:80]
##                chains |= self.parse_REMARK350_chains(line_chains)
##
##            elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
##                line_chains = lines[j][53:80]
##                chains |= self.parse_REMARK350_chains(line_chains)
##
##            elif ',' in lines[j][11:80]:
##                if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
##                    line_chains = lines[j][11:80]
##                    chains |= self.parse_REMARK350_chains(line_chains)
##
##            ## count and parse chain transformations
##            ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj.pdb)
##            elif lines[j][13:19] in ['BIOMT3','SMTRY3']:
##
##                matrixno = int(lines[j][19:24])
##                ## parse transformation matrix
##                matrixrow1 = lines[j-2][24:].split()
##                matrixrow2 = lines[j-1][24:].split()
##                matrixrow3 = lines[j-0][24:].split()
##                matrixrows = [matrixrow1,matrixrow2,matrixrow3,]
####                ## find out whether transformation matrix yields a transformation
####                transformation = False
####                for k in range(3):
####                    ## add a zero translation vector if a translation vector is not given
####                    if len(matrixrows[k]) == 3:
####                        matrixrows[k] += [0.]
####                    if float(matrixrows[k][k]) == 1. and float(matrixrows[k][3]) == 0.:
####                        continue
####                    else:
####                        transformation = True

    return d_header


def parse_recordREMARK465(line, lines, i, d_header, remark):

    ## missing residues

    if line[10:].strip() in [
        'M RES C SSSEQI',
        'M RES C  SSEQI',
        ]:

        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK %s' %(remark):
                break

            try:
                model = int(lines[j][12:14])
            except:
                model = 1
            res_name = lines[j][15:18].strip()
            chain = lines[j][19]
            res_no = int(lines[j][22:26])
            iCode = lines[j][26]

            if model != 1: ## e.g. 1ohh
                return d_header

            if not 'REMARK%s' %(remark) in d_header.keys():
                d_header['REMARK%s' %(remark)] = {}
            if not 'chains' in d_header['REMARK%s' %(remark)].keys():
                d_header['REMARK%s' %(remark)]['chains'] = {}
            if not chain in d_header['REMARK%s' %(remark)]['chains'].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain] = {}
            if not 'residues' in d_header['REMARK%s' %(remark)]['chains'][chain].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'] = {}
            if not res_no in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no] = {}

            ## chain > seq
            if not 'seq' in d_header['REMARK%s' %(remark)]['chains'][chain].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['seq'] = []
            d_header['REMARK%s' %(remark)]['chains'][chain]['seq'] += ['%4i%1s' %(res_no,iCode,)]

            ## res_no > d_iCodes
            if not 'd_iCodes' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
            ## d_iCodes > iCode
            if not iCode in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes']:
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

            ## iCode > record
            d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = 'REMARK465'

            ## iCode > altlocs
            if not 'altlocs' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'] = {}
            ## altlocs > altloc
            altloc = ' '
            if not altloc in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc] = {}
            ## altloc > res_name
            d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name'] = res_name

            ## res_no > l_iCodes
            if not 'l_iCodes' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['l_iCodes'] = []
            ## l_iCodes > iCode
            if not iCode in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['l_iCodes']:
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

            ## iCode > REMARK
            if not 'REMARK' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] = True

            ## iCode > record
            d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = 'REMARK%s' %(remark)


    return d_header


def parse_recordREMARK470(line, lines, i, d_header, remark,):

    ## missing atoms

    ## the latter equation is only to acommodate for 1fvk.pdb
    if line[10:].strip() == 'M RES CSSEQI  ATOMS' or line[10:].strip() == 'M RES C SEQI  ATOMS':

        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK 470':
                break

            ## model M
            try:
                model = int(lines[j][11:13])
            except:
                model = 'N/A'

            ## res_name RES
            res_name = lines[j][15:18].strip()

            ## chain C
            chain = lines[j][19]

            ## res_no SSEQ
            try:
                res_no = int(lines[j][20:24])
            except:
                res_no = lines[j][20:24].strip()

            ## iCode I
            iCode = lines[j][24]

            ## atoms ATOMS
            atoms = lines[j][25:].split()

            ##
            ## write to dictionary
            ##
            if not 'REMARK%s' %(remark) in d_header.keys():
                d_header['REMARK%s' %(remark)] = {}
            if not 'chains' in d_header['REMARK%s' %(remark)].keys():
                d_header['REMARK%s' %(remark)]['chains'] = {}
            if not chain in d_header['REMARK%s' %(remark)]['chains'].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain] = {}
            if not 'residues' in d_header['REMARK%s' %(remark)]['chains'][chain].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'] = {}
            if not res_no in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no] = {}

            ## res_no > l_iCodes
            if not 'l_iCodes' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['l_iCodes'] = []
            ## l_iCodes > iCode
            if not iCode in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['l_iCodes']:
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

            ## res_no > d_iCodes
            if not 'd_iCodes' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
            ## d_iCodes > iCode
            if not iCode in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes']:
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

            ## iCode > atoms
            if not 'atoms' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
            ## atoms > atom_name > coordinate
            for atom_name in atoms:
                if not atom_name in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                    d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                ## iCode > REMARK
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['REMARK'] = True

            ## iCode > altlocs
            if not 'altlocs' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'] = {}
            ## altlocs > altloc
            altloc = ' '
            if not altloc in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc] = {}
            ## altloc > res_name
            if not 'res_name' in d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc].keys():
                d_header['REMARK%s' %(remark)]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name'] = res_name


    elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
        print pdb1
        print pdb2
        notexpected

    return d_header


def parse_recordSEQRES(line, d_header):

    chain = line[11]

    if 'SEQRES' not in d_header.keys():
        d_header['SEQRES'] = {}
    if 'chains' not in d_header['SEQRES'].keys():
        d_header['SEQRES']['chains'] = {}
    if chain not in d_header['SEQRES']['chains'].keys():
        d_header['SEQRES']['chains'][chain] = {}
    if not 'type' in d_header['SEQRES']['chains'][chain].keys():
        d_header['SEQRES']['chains'][chain]['type'] = 'N/A'

    l_residues = line[19:70].split()

    s_residues = ''
    for i in range(len(l_residues)):
        residue = l_residues[i]
        if residue in d_pep.keys():
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
            s_residues += d_pep[residue]
        elif residue in l_nuc:
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'nucleotide'
            s_residues += residue
        elif residue in ['GLC','GAL','MAN','FRU']:
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'saccharide'
            s_residues += residue
        else:
            if residue == 'UNK': ## e.g. 1pny.pdb
                if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                    d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
            s_residues += 'X'

    if 'seq' not in d_header['SEQRES']['chains'][chain].keys():
        d_header['SEQRES']['chains'][chain]['seq'] = ''
    d_header['SEQRES']['chains'][chain]['seq'] += s_residues

    if 'seq3' not in d_header['SEQRES']['chains'][chain].keys():
        d_header['SEQRES']['chains'][chain]['seq3'] = []
    d_header['SEQRES']['chains'][chain]['seq3'] += l_residues

    return d_header


def parse_recordATOM(
    line, d_coordinates, lines, i, d_header, record,
    parse_ligands = True,
    ):

    import numpy

    atom_no = int(line[6:11])
    atom_name = line[12:16].strip()
    altloc = line[16]
    res_name_ATOM = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    element = line[76:78].strip()
    coordinate = numpy.array([x, y, z])
    tempFactor = float(line[60:66])

    skip = False

    ## standard or modified residue? parse seq?
    bool_SEQRESres = check_if_SEQRESres(res_name_ATOM,record,d_header,chain,res_no,iCode,)
    if bool_SEQRESres == False:
        parse_atom_seq = False

    if parse_ligands == False and bool_SEQRESres == False:
        return d_coordinates, {
            'chain':chain,
            'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
            'atom_no':atom_no,'atom_name':atom_name,'element':element,
            }

##    if 'SEQRES' in d_header.keys():
##        if chain in d_header['SEQRES']['chains']:
##            type = d_header['SEQRES']['chains'][chain]['type']
##            if type != 'peptide':
##                parse_atom_seq = False
##        else:
##            if parse_atom_seq == True:
##                stop1
##            parse_atom_seq = False

    if not 'chains' in d_coordinates.keys():
        d_coordinates['chains'] = {}
    if not chain in d_coordinates['chains'].keys():
        d_coordinates['chains'][chain] = {}
    if not 'residues' in d_coordinates['chains'][chain].keys():
        d_coordinates['chains'][chain]['residues'] = {}

    ## chain > seq
    if not 'seq' in d_coordinates['chains'][chain].keys():
        d_coordinates['chains'][chain]['seq'] = []
    if not '%4i%1s' %(res_no,iCode,) in d_coordinates['chains'][chain]['seq']:
        d_coordinates['chains'][chain]['seq'] += ['%4i%1s' %(res_no,iCode,)]

    ## residues > res_no
    if not res_no in d_coordinates['chains'][chain]['residues'].keys():
        d_coordinates['chains'][chain]['residues'][res_no] = {}

    ## res_no > d_iCodes
    if not 'd_iCodes' in d_coordinates['chains'][chain]['residues'][res_no].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
    ## d_iCodes > iCode
    if not iCode in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes']:
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

    ## iCode > altlocs
    if not 'altlocs' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'] = {}
    ## altlocs > altloc
    if not altloc in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc] = {}
    ## iCode > res_name
    if not 'res_name' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name'] = res_name_ATOM
    ## iCode > record
    d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['record'] = record
    
    ## check that res_name is correct (e.g. 2fes:L:1)
    if d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name'] != res_name_ATOM: ## 1fh2:A:30 altloc
        print d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
        print res_name_ATOM
        print altloc
        print chain, res_no, iCode
        print line
        stop_add_altloc

    ## res_no > l_iCodes
    if not 'l_iCodes' in d_coordinates['chains'][chain]['residues'][res_no].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] = []
    ## l_iCodes > iCode
    if not iCode in d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes']:
        d_coordinates['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

    ## iCode > atoms
    if not 'atoms' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
    ## atoms > atom_name > coordinate,element
    if not atom_name in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate,'element':element,'tempFactor':tempFactor,}

    ## iCode > record
    if not 'record' in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = record

    return d_coordinates, {
        'chain':chain,
        'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
        'atom_no':atom_no,'atom_name':atom_name,'element':element,
        }

def parse_atom_sequence(
    d_ATOMseq,d_header,
    lines, i,
    record, chain, res_no, iCode, altloc, res_name_ATOM,
    ):

    if not chain in d_ATOMseq.keys():
        d_ATOMseq[chain] = {
            'seq':[],'iCodes':[],'res_nos':[],'altlocs':[],'records':[],
            }

    if lines[i-1][:6].strip() in ['ATOM','HETATM','ANISOU','SIGUIJ','SIGATM',]:
        res_no_prev = int(lines[i-1][22:26])
        iCode_prev = lines[i-1][26]
        altloc_prev = lines[i-1][16]
    else:
        res_no_prev = None
        iCode_prev = None
        altloc_prev = None

    ## 2zfo (atlloc B in SEQRES)
    skip = False
    if chain in d_ATOMseq.keys():
        if len(d_ATOMseq[chain]['seq']) > 0:
            if (
                d_ATOMseq[chain]['res_nos'][-1] == res_no and
                d_ATOMseq[chain]['iCodes'][-1] == iCode
                ):
                ## skip if altloc A not added
                skip = True

    if (
        (
            len(d_ATOMseq[chain]['seq']) == 0 or
            (
                not (
                    res_no == res_no_prev and
                    iCode == iCode_prev and
                    altloc == altloc_prev ## 2zfo (atlloc B in SEQRES)
                    )
                )
            ) and
        skip == False
        ):


        res_name_REMARK465 = None
        REMARK465_before_ATOM = False
        if 'REMARK465' in d_header.keys():
            if chain in d_header['REMARK465']['chains'].keys():

                pass_if = False

                if res_no_prev != None:
                    ## REMARK465 before ATOM (gap between REMARK465 and next ATOM but not prev ATOM)
                    ## 2nv7
                    if res_no_prev+1 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                        pass_if = True
                    ## REMARK465 before ATOM (gap between REMARK465 and next ATOM and prev ATOM)
                    ## 2pqj (not 2qqh)
                    if set(range(res_no_prev+1,res_no)) & set(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                        pass_if = True

                if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                    ## REMARK465 before or after ATOM
                    if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                        pass_if = True
                    ## REMARK465 before ATOM
                    if res_no-1 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                        pass_if = True
                    ## REMARK465 before ATOM (no zero residue)
                    ## 1b9n,2fxm
                    if res_no_prev == None and res_no >= 1 and -1 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                        pass_if = True
                    ## REMARK465 before ATOM (first residue = 0 and second residue > 1) ## e.g 2asd,1ca5
                    if res_no_prev == None and res_no > 1 and 0 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                        pass_if = True
                    ## REMARK465 before ATOM (first residue > 1 and second residue > 1) ## e.g 2h27
                    if res_no_prev == None and res_no > 1 and res_no > min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                        pass_if = True
                    ## REMARK465 before ATOM
                    ## 1sgf,1nu0
                    if res_no_prev == min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                        pass_if = True

                if pass_if == True:

                    l_REMARK465_res_nos = d_header['REMARK465']['chains'][chain]['residues'].keys()
                    l_REMARK465_res_nos.sort()
                    l_REMARK465_res_nos_grouped = [[]]
                    for j in range(len(l_REMARK465_res_nos)):
                        res_no_REMARK465 = l_REMARK465_res_nos[j]
                        l_REMARK465_res_nos_grouped[-1] += [res_no_REMARK465]
                        if j == len(l_REMARK465_res_nos)-1:
                            break
                        if not res_no_REMARK465 == l_REMARK465_res_nos[j+1]-1:
                            l_REMARK465_res_nos_grouped += [[]]
                    l_REMARK465_res_names_grouped = []
                    for l_REMARK465_res_nos in l_REMARK465_res_nos_grouped:
                        l_REMARK465_res_names_grouped += [[]]
                        for res_no_REMARK465 in l_REMARK465_res_nos:
                            l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                            for iCode_REMARK465 in l_iCodes_REMARK465:
                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                                l_REMARK465_res_names_grouped[-1] += [res_name_REMARK465]

                    ## find REMARK465 seq before ATOM seq (if any)
                    REMARK465_before_ATOM_loop = True
                    l_REMARK465_res_nos = []
                    while REMARK465_before_ATOM_loop == True:
                        REMARK465_before_ATOM_loop = False
                        for j in range(len(l_REMARK465_res_names_grouped)):
                            l_REMARK465_res_names = l_REMARK465_res_names_grouped[j]
                            l_SEQRES_res_names = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])+len(l_REMARK465_res_nos):][:len(l_REMARK465_res_names)]
##                            print l_REMARK465_res_names
##                            print l_SEQRES_res_names
##                            stop
##                            SEQRESseq_C = d_header['SEQRES']['chains'][chain]['seq3'][-len(l_REMARK465_res_names):]
                            if len(l_REMARK465_res_names) > 1 and l_REMARK465_res_names == l_SEQRES_res_names:
                                REMARK465_before_ATOM_loop = True
                            ## single residue insertion
                            elif (
                                len(l_REMARK465_res_names) == 1 and
                                l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                                res_name_ATOM != l_SEQRES_res_names[0]
                                ):
                                REMARK465_before_ATOM_loop = True
                            ## single residue insertion, res_no_465 < res_no_ATOM
                            ## 2hu9
                            elif (
                                len(l_REMARK465_res_names) == 1 and
                                l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                                res_name_ATOM == l_SEQRES_res_names[0] and
                                res_no > min(l_REMARK465_res_nos_grouped[j])
                                ):
                                REMARK465_before_ATOM_loop = True
                            ## REMARK465 not before ATOM
                            ## e.g. 3bef,1bd7,4htc
                            else:
                                None
                            if REMARK465_before_ATOM_loop == True:
                                REMARK465_before_ATOM = True
                                l_REMARK465_res_nos += l_REMARK465_res_nos_grouped[j]
                                del l_REMARK465_res_names_grouped[j]
                                del l_REMARK465_res_nos_grouped[j]
                                break

        try:
            res_name_SEQRES = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
        except:
            res_name_SEQRES = 'N/A'

        if REMARK465_before_ATOM == False and res_name_ATOM != res_name_SEQRES and res_name_ATOM in ['FOR','HOA','ACE',] and (res_no == min(d_coordinates['chains'][chain]['residues'].keys()) or res_no == max(d_coordinates['chains'][chain]['residues'].keys())):
            fd = open('formyl.txt','a')
            fd.write('%s %s %s\n' %(pdb, chain, res_name_ATOM))
            fd.close()
            bool_return = True
            return d_ATOMseq, d_header, bool_return

        ## REMARK465 after ATOM (not appended)
        if REMARK465_before_ATOM == False and res_name_ATOM == res_name_SEQRES:# and res_name_REMARK465 != res_name_SEQRES:
            d_ATOMseq,d_header = append_sequence(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,None,False,True,)
        ## REMARK465 before ATOM (appended)
        elif REMARK465_before_ATOM == True:# and res_name_ATOM != res_name_SEQRES:# and res_name_REMARK465 == res_name_SEQRES:
            d_ATOMseq,d_header = append_sequence(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,True,)
        else:
            try:
                print '---'
                SEQRES_res_name = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
                SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print chain, res_no, iCode
                print lines[i]
                print d_ATOMseq[chain]['seq']
                print d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print res_name_ATOM, res_name_SEQRES
                print 'altloc', altloc
                if 'REMARK465' not in d_header.keys() and atom_name == 'CA' and d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] == {}:
                    stop_remark470_records_missing
                    pass
                elif 'REMARK465' not in d_header.keys() and atom_name == 'N' and res_no == 2 and res_name_SEQRES == 'MET':
                    stop_remark465_records_missing_met1
                    pass
                ## 2zfo
                elif 'REMARK465' not in d_header.keys() and altloc != ' ':
                    if res_name_SEQRES == res_name_ATOM:
                        stop
                    pass
                elif 'REMARK465' not in d_header.keys() and altloc == ' ':
                    pass ## new, temp
                elif res_name_REMARK465 == None and chain in d_header['REMARK465']['chains'].keys() and min(d_header['REMARK465']['chains'][chain]['residues'].keys()) < res_no: ## e.g. 3c5w
                    print chain,res_no
                    print d_header['REMARK465']['chains'][chain]['residues'].keys()
                    stop_temp_broken
                    pass
                ## REMARK465 after ATOM
                elif (
                    (iCode == ' ' or len(d_ATOMseq[chain]['seq']) == 0) and
                    res_no <= min(d_header['REMARK465']['chains'][chain]['residues'].keys()) and
                    res_name_ATOM == SEQRES_res_name
                    ):
                    stop_temp_example_commentout
                    d_ATOMseq[chain]['seq'] += [res_name_ATOM]
                    d_ATOMseq[chain]['res_nos'] += [res_no]
                    d_ATOMseq[chain]['iCodes'] += [iCode]
                    d_ATOMseq[chain]['records'] += [record]
                else:
                    ## 1jly
                    if res_name_ATOM in ['FOR','HOA','ACE',] and (res_no == min(d_coordinates['chains'][chain]['residues'].keys()) or res_no == max(d_coordinates['chains'][chain]['residues'].keys())):
                        fd = open('formyl.txt','a')
                        fd.write('%s %s %s\n' %(pdb, chain, res_name_ATOM))
                        fd.close()
                        bool_return = True
                        return d_ATOMseq, d_header, bool_return
                    print '*******'
                    print 'ATOM  ', d_ATOMseq[chain]['seq']
                    print 'SEQRES', SEQRES_seq
                    print line
                    print chain,res_no
                    print 'SEQRES', SEQRES_res_name
                    print 'ATOM  ', res_name_ATOM, '***iCode***', iCode
                    print 'REMARK', res_name_REMARK465, iCode_REMARK465
                    stop_N_terminal
            except:
                print lines[i]
                print lines[0]
                print REMARK465_before_ATOM
                print res_name_ATOM, res_name_SEQRES
                print d_ATOMseq[chain]['seq']
                stop_temp
                pass
                

        if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]:
            print '*******'
            print 'ATOM  ', d_ATOMseq[chain]['seq']
            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
            print res_name_ATOM, res_name_SEQRES
            print lines[i]
            print len(d_ATOMseq[chain]['seq'])
            print res_name_REMARK465
            print lines[0].strip()[63:80], chain, res_no, iCode
            stop_sequence_difference

    bool_return = False
    return d_ATOMseq, d_header, bool_return

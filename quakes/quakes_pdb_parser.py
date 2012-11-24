import copy, os
import sys
sys.path.append('/home/tc/svn/tc_sandbox/pdb/')
import biounit,parse_pdb
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb')
import smallmolecules
import quakes_smallmolecules

d_res1 = {
    'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
    'MSE':'M','UNK':'X','ASX':'X','GLX':'X',
    }

s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

l_nucleotides = [
     'A', 'C', 'G', 'U',      'I', ## ribonucleotides
    'DA','DC','DG',     'DT','DI', ## deoxyribonucleotides
    'N', ## N is any 5'-monophosphate nucleotide
    ]

## HETATM res_names for which coordinates are parsed
d_modres = {
    'MSE':'MET', ## selenomethionine
##            ## phosphorylation
##            'TPO':'THR',
##            'SEP':'SER',
##            'PHD':'ASP',
##            'PTR':'TYR',
    }


d = quakes_smallmolecules.main()

print d.keys()
l_terminalmodres = d['terminalmodres']
l_lpeptidelinkers = d['lpeptidelinkers']
l_dpeptidelinkers = d['dpeptidelinkers']

## connectivity *and* identity of solutes *and* ions not checked
l_solutes = d['solutes']
d_ions = d['ions']

l_clusters = d['clusters']
l_prosthetic_groups = d['prosthetic groups']
l_coenzymes = d['coenzymes']

## connectivity of cofactor *molecules* are not checked
l_cofactors = l_clusters+l_prosthetic_groups+l_coenzymes+d_ions.keys()

## connectivity of metal *atoms* are not checked (too many examples)
l_atoms_metal = d['metals']

d_saccharides = d['saccharides']

def parse_mmCIF_mutation(
    pdb,chain,res_no, res_wt,res_mutant,
    path_mmCIF,
    ):

    bool_mutation = False
    s = ''

    d_characters = {
        "'":"'",
        '"':'"',
        '(':')',
        ';':';',
        }

    fd = open('%s/%s/%s.cif' %(path_mmCIF,pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()
    j_line = 0
    for i_line in range(len(lines)):
        if i_line < j_line:
            continue
        ## new set of variables
        if lines[i_line].strip() == '#':
            k_line = j_line
            for j_line in range(i_line+1,len(lines)):
                if j_line < k_line:
                    continue
                if lines[j_line].strip() == 'loop_':
                    loop = True
                    ## loop over variables
                    l_variables = []
                    for k_line in range(j_line+1,len(lines)):
                        if lines[k_line][0] == '_':
                            l_variables += [lines[k_line].strip()]
                        else:
                            break
                    continue
                ## not a loop
                if lines[j_line][0] == '_':
                    variable = lines[j_line].split()[0]
                    if variable == '_struct.pdbx_descriptor':
                        if lines[j_line].strip() == variable:
                            value = lines[j_line+1].strip()
                        else:
                            value = lines[j_line][len(variable):].strip()
                        s = '_struct.pdbx_descriptor "%s"' %(value)
                        if (
                            '%1s%i%1s' %(res_wt,res_no,res_mutant,) in value
                            or
                            '%1s(%1s %i)%1s' %(res_wt,chain,res_no,res_mutant,) in value
                            ):
                            bool_mutation = True
                            break
                ## line break
                elif lines[j_line] == '':
                    stop
                    continue
                ## #
                elif lines[j_line].strip() == '#':
                    continue
                ## loop
                else:

                    variable = '_entity.pdbx_mutation'
                    if not variable in l_variables:
                        continue
                    index = l_variables.index(variable)

                    ## loop over values
                    l_values = []
                    l_characters = []
                    value = ''
                    for k_line in range(j_line,len(lines)):
                        if len(l_values) >= len(l_variables):
                            break
                        line = lines[k_line].rstrip()
                        l = lines[k_line].rstrip().split("'")
                        if (
                            '"' in lines[k_line]
                            or
                            "'" in lines[k_line]
                            or
                            "(" in lines[k_line]
                            or
                            ';' == lines[k_line][0]
                            or
                            l_characters != [] ## 1abw
                            ): ## 1wt1, 1a0a
                            i_character2 = 0
                            for i_character1 in range(len(line)):
                                if i_character1 < i_character2:
                                    continue
                                if line[i_character1] == ' ':
                                    continue
##                                if i_character1 == 0:
##                                    if line[i_character1] in [' ','"',"'",'(',]:
##                                        print line
##                                        stop
##                                    for i_character2 in range(i_character1+1,len(line)):
##                                        if line[i_character2] == ' ':
##                                            value = line[i_character1:i_character2]
##                                            l_values += [value]
##                                            break

                                ## multi-line value
                                if line[i_character1] == ';':

                                    if i_character1 != 0:
                                        print line
                                        stop

                                    ## semicolon end
                                    if l_characters == [';']:
                                        if line.strip() != ';':
                                            print line
                                            stop
                                        l_characters = []
                                        value += line[1:].strip()
                                        l_values += [value]
                                        value = ''
                                        break ## next line

                                    ## semicolon start
                                    elif l_characters == []:
                                        l_characters += [line[i_character1]]
                                        value += line[1:].strip()
                                        break ## next line

                                    else:
                                        print l_characters
                                        print line
                                        stop

                                ## multi-line continuation
                                elif l_characters == [';'] and ';' not in line:
                                    value += line.strip()
                                    break ## next line
    
                                ## quoted value
                                elif line[i_character1] in ['"',"'",'(',]:

##                                    if len(l_characters) > 0 and line[i_character1] == d_characters[line[i_character1]]:
##                                        l_characters = []
##                                        print line
##                                        stop
##                                    elif len(l_characters) > 0:
##                                        stop
##                                    else:
                                    l_characters += [line[i_character1]]

                                    if lines[k_line].strip() == ';':
                                        print lines[k_line-1]
                                        print l_characters
                                        print value
                                        stop
                                    if len(l_characters) > 1:
                                        print l_characters
                                        stop

                                    ## quoted value
                                    for i_character2 in range(i_character1+1,len(line)):
                                        if line[i_character2] == d_characters[line[i_character1]]:
                                            if len(l_characters) == 1:
                                                i_character2 += 1
                                                value = line[i_character1+1:i_character2-1]
                                                if value.strip() == '' or value.strip() == ';':
                                                    print value
                                                    print line
                                                    stop
                                                l_values += [value]
                                                l_characters = []
                                                value = ''
                                                break
                                            else:
                                                l_characters = l_characters[:-1]
                                        elif line[i_character2] == line[i_character1]:
                                            l_characters += [line[i_character2]]

                                ## space
                                elif line[i_character1] == ' ':
                                    continue

                                ## non-quoted value
                                else:
                                    if i_character1+1 == len(line):
                                        if line[i_character1] == ' ':
                                            stop
                                        value = line[i_character1]
                                        l_values += [value]
                                        value = ''
                                        if len(l_characters) > 0:
                                            print l_characters
                                            print value
                                            print line
                                            stop
                                    for i_character2 in range(i_character1+1,len(line)):
                                        if line[i_character2] == ' ':
                                            value = line[i_character1:i_character2]
                                            l_values += [value]
                                            if len(l_characters) > 0:
                                                print l_characters
                                                print value
                                                print line
                                                stop
                                            value = ''
                                            break

                        elif len(l) == 1:
                            l_values += lines[k_line].split()

                        else:
                            for s in l:
                                if s == '': ## first character is '
                                    continue
                                if s[0] == ' ' or s[-1] == ' ':
                                    l_values += s.split()
                                else:
                                    l_values += [s]

                    if len(l_variables) != len(l_values):
                        print lines[j_line]
                        print lines[j_line+1]
                        print lines[j_line+2]
                        print lines[k_line]
                        print
                        print 'variables', l_variables
                        print 'values', l_values
##                                    print lines[j_line].split("'")
##                                    print lines[j_line].split()
                        print
                        for i in range(min(len(l_variables),len(l_values))):
                            print l_variables[i], l_values[i]
                        print pdb
                        stop

                    value = l_values[index]
                    s = '%s "%s"' %(variable,value,)
                    if value != '?':
                        if (
                            '%1s%i%1s' %(res_wt,res_no,res_mutant,) in value
                            or
                            '%1s(%1s %i)%1s' %(res_wt,chain,res_no,res_mutant,) in value
                            ):
                            bool_mutation = True
                            break

                    ## break j_line loop
                    if bool_mutation == True:
                        stop
                        break

                ## break j_line loop
                if bool_mutation == True:
                    break
        ## break i_line loop
        if bool_mutation == True:
            break

    return bool_mutation, s


def parse_header(
    s_pdb, path_pdb,
    min_len_chain,
    stop_error = True,
    ):

    d_crystalsynonyms = {
        'P 2 21 21':'P 21 21 2',

        'P 21':'P 1 21 1',

        'C 2':'C 1 2 1',
        'C 21':'C 1 2 1',
        'C 1 21 1':'C 1 2 1',

        'P 2':'P 1 2 1',
        ## errors, temporary until fixed/remediated
        'P 21 21 2 A':'P 21 21 2',
        'P 2 21 21':'P 2 21 21', ## 2pnj
        }

    print 'parsing header of', s_pdb

    ##
    ## read lines
    ##
    fd = open('%s/%s/pdb%s.ent' %(path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
##        fd = open('C:\Users\Tommy Carstensen\Downloads\pdb\%s.pdb' %(s_pdb.lower(),),'r')
    lines = fd.readlines()
    fd.close()

    ##
    ## set dictionaries
    ##
    d_header = {}
    d_conect = {}
    l_hetatms = []
    parse_atoms = False

    ## insertion chain, res_no, res_name
    prev_chain = ''
    d_insertions = {}
    biounit = 'N/A'

    for i in range(len(lines)):
        line = lines[i]

        record = line[:6].strip()

        if record == 'ATOM': ## section 9
            continue
##                if parse_atoms == False:
##                    if 'REMARK350' in d_header.keys() and 'HET' in d_header.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = parse_atom_no_range(d_conect, 'atom', atom_no)

        elif record == 'HETATM': ## section 9
            continue
##                if parse_atoms == False:
##                    if 'REMARK350' in d_header.keys() and 'HET' in d_header.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = parse_atom_no_range(d_conect, 'atom', atom_no)

        elif record == 'COMPND': ## section xxx
            d_header = parse_recordCOMPND(d_header, line, i, lines)

        elif record == 'SOURCE':

            if not 'SOURCE' in d_header.keys():
                d_header['SOURCE'] = {}

            if lines[i][10:18].strip() == 'MOL_ID:':
                mol_id = int(lines[i][18:].replace(';',''))
            else:
                continue
            
            for j in range(i+1,len(lines)):
                
                if lines[j][:6] != 'SOURCE':
                    break
                if lines[j][10:18].strip() == 'MOL_ID:':
                    break

                if 'ORGANISM SCIENTIFIC' in lines[j]:
                    print lines[j]
                    stop
                if 'ORGANISM_SCIENTIFIC' in lines[j]:
                    if lines[j][10:32] != ' ORGANISM_SCIENTIFIC: ':
                        print lines[j]
                        stop
                    organism = lines[j][32:80].strip()[:-1]
                    if ';' in organism:
                        print lines[j]
                        stop
                    if mol_id in d_header['SOURCE'].keys():
                        print d_header['SOURCE']
                        print lines[j]
                        stop
                    d_header['SOURCE'][mol_id] = organism

        elif record == 'MDLTYP':
            if line[10:].strip() == 'MINIMIZED AVERAGE':
                continue
            if (
                line[9] not in ['2','3',]
                ## protein
                and line[:31] != 'MDLTYP    CA ATOMS ONLY, CHAIN '
                ## rna/dna
                and line[:30] != 'MDLTYP    P ATOMS ONLY, CHAIN '
                ## NMR
                and line[:50] != 'MDLTYP    MINIMIZED AVERAGE; CA ATOMS ONLY, CHAIN '
                ):
                print line
                print line[:31]
                stop
##                chain = line[31]
##                d_header['MDLTYP'] = [chain]
            d_header['MDLTYP'] = True

        elif record == 'SSBOND':
            chain1 = line[15]
            res_no1 = int(line[17:21])
            iCode1 = line[21]
            chain2 = line[29]
            res_no2 = int(line[31:35])
            iCode2 = line[35]
            if line[59:65].strip() != '1555':
                print line
                stop
            if line[66:72].strip() != '1555':
                if not 'SSBOND' in d_header.keys():
                    d_header['SSBOND'] = []
                d_header['SSBOND'] += [[chain1,res_no1,iCode1,chain2,res_no2,iCode2,]]

        elif record == 'REMARK': ## section 2
            d_header = parse_recordREMARK(d_header, line, i, lines)

        elif record == 'SEQADV': ## section 3
            d_header = parse_recordSEQADV(line, d_header)

        elif record == 'CAVEAT': ## Title Section
            if line[6:10].strip() == '' and 'CHIRAL' not in line and 'GEOMETR' not in line and 'STEREOCHEMISTRY' not in line:
                if s_pdb not in [
                    '1r9m',
                    '3ckz', ## PEPTIDE BACKBONE LINKAGE PROBLEM
                    '3cl0', ## PEPTIDE BACKBONE LINKAGE PROBLEM
                    '1g4b', ## NUMEROUS CLOSE CONTACTS. SEE REMARK 3
                    '3eyo', ## CLOSE CONTACTS BETWEEN SIDE CHAINS
                    '1heg', ## ABNORMALLY CLOSE CONTACTS PRESENT
                    '3dwg', ## ABNORMALLY SHORT LINK BETWEEN LYS A51 AND PLP A401
                    '3dwi', ## ABNORMALLY SHORT LINK BETWEEN LYS A51 AND PLP A401
                    '3fg4', ## C-N PEPTIDE BOND LENGTH ERROR D414-D418
                    '1fx7', ## continuation across two lines
                    '1euq', ## continuation across two lines
                    '3dwf', ## continuation across three lines
                    ]:
                    print line
                    stop

        elif record == 'SPLIT': ## Title Section
            d_header['SPLIT'] = True

        elif record == 'SEQRES': ## section 3
            d_header = parse_recordSEQRES(line, d_header)

        elif record == 'HET': ## section 4
            hetID = line[7:10].strip()
            ## continue if water
            if hetID in ['HOH','H2O','DOD','D2O']: ## D2O in 2JAJ
                continue
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
                d_header['HET'][chain][res_no][iCode] = set() ## multiple hetIDs possible if altlocs
            d_header['HET'][chain][res_no][iCode] |= set([hetID])

        elif record == 'MODRES': ## section 3
            hetID = line[12:15].strip()
            chain = line[16]
            res_no = int(line[18:22])
            iCode = line[22]
            res_name = line[24:27].strip()
            txt = line[29:80].strip()
            if hetID in set(d_res1.keys())-set(['MSE']) and res_name in set(d_res1.keys())-set(['MSE']):
                continue
##                if txt not in [] and hetID in set(d_res1.keys()+l_nucleotides)-set(['MSE']):
##                    print line
##                    print s_pdb
##                    stop_duplicate_modres_records_or_std_res_name
            if 'MODRES' not in d_header.keys():
                d_header['MODRES'] = {}
            if chain not in d_header['MODRES'].keys():
                d_header['MODRES'][chain] = {}
            if res_no not in d_header['MODRES'][chain].keys():
                d_header['MODRES'][chain][res_no] = {}
            if iCode not in d_header['MODRES'][chain][res_no].keys():
                d_header['MODRES'][chain][res_no][iCode] = set() ## multiple hetIDs possible if altlocs
            d_header['MODRES'][chain][res_no][iCode] |= set([hetID])

        elif record == 'TITLE': ## section 2
            if not 'TITLE' in d_header.keys():
                d_header['TITLE'] = line[10:].strip()
            else:
                if d_header['TITLE'][-1] == '-':
                    d_header['TITLE'] += line[10:].strip()
                else:
                    d_header['TITLE'] += ' '+line[10:].strip()

        elif record == 'HELIX': ## section 5
            d_header = parse_pdb.parse_recordHELIX(line,d_header,)

        elif record == 'SHEET': ## section 5
            d_header = parse_pdb.parse_recordSHEET(line,d_header,)

        elif record == 'HEADER':
            d_header['HEADER'] = line[10:50].strip()

        elif record == 'EXPDTA': ## section 2
            methods = line[10:].strip().split(',')[0]
            if methods[:3] == 'NMR':
                methods = 'NMR'
            elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                methods = 'X-RAY'
            d_header['EXPDTA'] = methods

        elif record == 'AUTHOR':
            l_authors = line[10:].strip().split(',')
            if not 'AUTHOR' in d_header.keys():
                d_header['AUTHOR'] = []
            d_header['AUTHOR'] += l_authors

        elif record == 'CRYST1': ## section 8
            spacegroup = line[55:66].strip()
            if spacegroup in d_crystalsynonyms.keys():
                print 'remediate?', s_pdb
                print 'spacegroup', spacegroup
                spacegroup = d_crystalsynonyms[spacegroup]
            d_header['CRYST1'] = spacegroup

        else:
            continue

    ##
    ## missing records
    ##
    d_dics = {
        'chains':{},
        'HET':{},
        'MODRES':{},
        'REMARK200':{'TEMPERATURE':'N/A','PH':'N/A'},
        'REMARK525':[],
        'TITLE':'N/A',
        'EXPDTA':'N/A',
        'REMARK2':'N/A',
        }
    for key in d_dics.keys():
        if key not in d_header.keys():
            d_header[key] = d_dics[key]

    proteinchains = []
    peptidechains = []
    nucleotidechains = []
    saccharidechains = []
    for chain in d_header['SEQRES']['chains'].keys():
        if d_header['SEQRES']['chains'][chain]['type'] == 'peptide':
            peptidechains += chain
            if len(d_header['SEQRES']['chains'][chain]['seq']) > min_len_chain:
                proteinchains += chain
        elif d_header['SEQRES']['chains'][chain]['type'] == 'nucleotide':
            nucleotidechains += chain
        elif d_header['SEQRES']['chains'][chain]['type'] == 'saccharide':
            saccharidechains += chain

    d_header['proteinchains'] = proteinchains

    return d_header


def parse_recordCOMPND(d_header, line, i, lines,):

    if not 'COMPND' in d_header.keys():
        d_header['COMPND'] = {}
    if lines[i][10:18].strip() == 'MOL_ID:':
        mol_id = int(lines[i][18:].strip()[:-1])
    else:
        return d_header
    
    for j in range(i+1,len(lines)):
        if lines[j][:6] != 'COMPND':
            break
        if lines[j][10:18].strip() == 'MOL_ID:':
            break
        if lines[j][10:17].strip() == 'CHAIN:':
            chains = lines[j][17:80].strip().replace(':','').replace(';','').replace(' ','').split(',')
            for chain in chains:
                if chain in d_header['COMPND'].keys():
                    print chain,d_header['COMPND']
                    stop
                d_header['COMPND'][chain] = mol_id

    return d_header


def parse_recordSEQADV(line, d_header,):

    comment = line[49:70].strip()

    ## deletion
    if comment in ['DELETION','ENGINEERED DELETION',] and line[12:15].strip() == '' and line[18:23].strip() == '':
        return d_header
    ## incorrect seq in seq db?
    if comment == 'SEE REMARK 999' and line[12:15].strip() == '' and line[18:23].strip() == '':
        return d_header
    ## residue not in ATOM records (no residue number)
    if line[18:22].strip() == '':
        return d_header

    ## chromophore
    if line[12:15].strip() in [
        ## chromophore in name
        'CRO','CR2','CSY','MDO','XYG','DYG','NYG','CR8','RC7','CH6','NRQ',
        'CFY','CWR','CQR','CRQ','AYG','IIC','CRX','C99','CLV','C12','XXY',
        'CR0','GYS','X9Q','CRW','GYC','CR7','QLG','CH7',
        ## chromophore not in name
        'IEY','CRK','5ZA','CSH','CRG','4F3','MFC','PIA','CR5','CRF','NYC',
        'CZO','CCY','CJO',
        ]: ## or 'CHROMOPHORE' in comment
        return d_header

    chain = line[16]
    seq_no = int(line[18:22])
    iCode = line[22]
    res_name_ref = line[39:42].strip()
    if not 'SEQADV' in d_header.keys():
        d_header['SEQADV'] = {}
    if not chain in d_header['SEQADV'].keys():
        d_header['SEQADV'][chain] = {}
    if line[12:15].strip() not in [
        'SUI', ## ASP,GLY
        '175', ## ALA, SER, GLY
        ]:
        if line[7:11] not in [
            ## altlocs / microheterogeneity
            '2CI1','2JGE','2JGJ','3DL7',
            ## two residues
            '2E3V',
            ]:
            if '%4i%1s' %(seq_no,iCode,) in d_header['SEQADV'][chain].keys():
                print seq_no, iCode
                print d_header['SEQADV'][chain].keys()
                print d_header['SEQADV'][chain]
                print line
                stop_duplicate_entry_maybe_chromophore
    d_header['SEQADV'][chain]['%4i%1s' %(seq_no,iCode,)] = {'res_name_seqdb':res_name_ref,'comment':comment,}

    return d_header


def parse_coordinates(
    s_pdb, d_header,
    path_pdb,
    verbose = False, parse_molecules = True,
    ):

    print 'parsing coordinates of', s_pdb

    ## deep copy header, because REMARK465 records are deleted during loop over biomolecules
    d_header = copy.deepcopy(d_header)

    ##
    ## read lines
    ##
    fd = open('%s/%s/pdb%s.ent' %(path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
##        fd = open('C:\Users\Tommy Carstensen\Downloads\pdb\%s.pdb' %(s_pdb.lower(),),'r')
    lines = fd.readlines()
    fd.close()

    ##
    ## set dictionaries
    ##
    d_atomnos = {}
    d_CONECT = {}
    d_coordinates = {}
    d_ATOMseq = {}
    model = 1

    ##
    ## loop over lines
    ##
    for i in range(len(lines)):
        line = lines[i]
        record = line[:6].strip()

        if record == 'ATOM':
            if model != 1:
                continue
            d_coordinates, d_line = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
            d_ATOMseq = build_ATOMseq(line, lines, i, d_ATOMseq, d_header, d_coordinates,)
            d_atomnos[d_line['atom_no']] = d_line

        elif record == 'HETATM':
            if model != 1:
                continue
            res_name = line[17:20].strip()
            if res_name == 'UNL':
                stop_write_code
            ## water
            if res_name in ['D2O','H2O',]:
                print s_pdb, res_name
                stop
            ## water must be parsed in case there is a connection to it (e.g. 2bpb)
            if res_name in ['HOH','DOD',]: ## DOD in 2d4j
                d_coordinates, d_line, = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
            ## MSE
            elif res_name in d_modres.keys():
                d_coordinates, d_line, = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
            ## (poly)saccharide or other hetero compound
            else:
                d_coordinates, d_line, = parse_pdb.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
            atom_no = d_line['atom_no']
            d_atomnos[atom_no] = d_line

            d_ATOMseq = build_ATOMseq(line, lines, i, d_ATOMseq, d_header, d_coordinates,) ## 1omw

        elif record == 'CONECT':
            d_CONECT = parse_recordCONECT(line, d_atomnos, d_CONECT)
            
        elif record == 'HET':
            hetID = line[7:10].strip()
            if 'HET' not in d_coordinates.keys():
                d_coordinates['HET'] = set()
            d_coordinates['HET'] |= set([hetID])
            
        elif record == 'MODEL':
            model = int(line.split()[1])

        elif record == 'ENDMDL': ## e.g. 2q4t
            break


    if parse_molecules == True:
        d_molecules = build_dictionary_of_molecules(d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=verbose,)
    else:
        d_molecules = {}

    for chain in d_header['SEQRES']['chains'].keys():

        ## skip if not peptide chain
        type = d_header['SEQRES']['chains'][chain]['type']
        if type != 'peptide':
            continue

        ## skip if short chain
        if len(d_header['SEQRES']['chains'][chain]['seq3']) < 5:
            continue

        ## skip if all residues are unknown
        if len(d_header['SEQRES']['chains'][chain]['seq3'])*['UNK'] == d_header['SEQRES']['chains'][chain]['seq3']:
            continue

        ## append remaining REMARK465 residues
        if 'REMARK465' in d_header.keys():
            if chain in d_header['REMARK465']['chains'].keys():
                if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                    l_REMARK465_seq = list(d_header['REMARK465']['chains'][chain]['seq'])
                    d_ATOMseq,d_header = append_ATOMseq(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_seq,True,False,)

        if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3']:
            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3']
            print 'ATOM  ', d_ATOMseq[chain]['seq']
            print chain
            print 'SEQRES', len(d_header['SEQRES']['chains'][chain]['seq3'])
            print 'ATOM  ', len(d_ATOMseq[chain]['seq'])
            for i in range(len(d_ATOMseq[chain]['seq'])):
                if d_header['SEQRES']['chains'][chain]['seq3'][:i] != d_ATOMseq[chain]['seq'][:i]:
                    print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][:i]
                    print 'ATOM', d_ATOMseq[chain]['seq'][:i]
                    break
            print s_pdb,chain
            stop_different_SEQRESseq_ATOMseq

    return d_coordinates, d_molecules, d_ATOMseq


def parse_recordSEQRES(line, d_header):

    chain = line[11]

    if 'SEQRES' not in d_header:
        d_header['SEQRES'] = {'chains':{},}
    if chain not in d_header['SEQRES']['chains'].keys():
        d_header['SEQRES']['chains'][chain] = {}
    if not 'type' in d_header['SEQRES']['chains'][chain].keys():
        d_header['SEQRES']['chains'][chain]['type'] = 'N/A'

    l_residues = line[19:70].split()

    s_residues = ''
    for i in range(len(l_residues)):
        residue = l_residues[i]
        if residue in d_res1.keys():
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
            elif d_header['SEQRES']['chains'][chain]['type'] != 'peptide': ## e.g. 1vq6
                d_header['SEQRES']['chains'][chain]['type'] = 'peptide'
            s_residues += d_res1[residue]
        elif residue in l_nucleotides:
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'nucleotide'
            elif d_header['SEQRES']['chains'][chain]['type'] != 'nucleotide':
                stop
            s_residues += residue
        elif residue in ['GLC','GAL','MAN','FRU']:
            if d_header['SEQRES']['chains'][chain]['type'] == 'N/A':
                d_header['SEQRES']['chains'][chain]['type'] = 'saccharide'
            elif d_header['SEQRES']['chains'][chain]['type'] != 'saccharide':
                stop
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


def parse_recordCONECT(line,d_atomnos,d_CONECT):

    ## parse atom nos
    atom_nos = []
    for j in range(6,31,5):
        if line[j:j+5] == '     ':
            break
        atom_nos += [int(line[j:j+5])]

    ## parse res nos of *hetero* atoms
    atom_no1 = atom_nos[0]

    ## continue if water
    if atom_no1 not in d_atomnos.keys():
        return d_CONECT
    
    chain1 = d_atomnos[atom_no1]['chain']
    res_no1 = d_atomnos[atom_no1]['res_no']
    iCode1 = d_atomnos[atom_no1]['iCode']
    res_name1 = d_atomnos[atom_no1]['res_name']
    atom_name1 = d_atomnos[atom_no1]['atom_name']
    element1 = d_atomnos[atom_no1]['element']

    ## skip if hydrogen atom
    if atom_name1[0] == 'H':
        return d_CONECT
    ## skip if a cofactor to which molecules are not *covalently* bound (use mmCIF to distinguish bond type instead..!)
    if res_name1 in l_cofactors+['HOH']+l_solutes:
        return d_CONECT
    ## skip if metal atom, because not covalent bond (use mmCIF to distinguish bond type instead..!)
    if element1 in l_atoms_metal: ## e.g. 1a6l
        return d_CONECT

    d_atom_nos = {'intra':[],'inter':[]}

    for atom_no2 in atom_nos[1:]:
        chain2 = d_atomnos[atom_no2]['chain']
        res_no2 = d_atomnos[atom_no2]['res_no']
        iCode2 = d_atomnos[atom_no2]['iCode']
        res_name2 = d_atomnos[atom_no2]['res_name']
        atom_name2 = d_atomnos[atom_no2]['atom_name']
        element2 = d_atomnos[atom_no2]['element']
        
        ## skip if hydrogen atom
        if atom_name2[0] == 'H':
            continue
        ## skip if a cofactor to which molecules are not *covalently* bound (use mmCIF to distinguish bond type instead..!)
        if res_name2 in l_cofactors+['HOH']+l_solutes:
            continue
        ## skip if metal atom, because not covalent bond (use mmCIF to distinguish bond type instead..!)
        if element2 in l_atoms_metal:
            continue
        if (
            chain1 !=  chain2 or
            res_no1 != res_no2 or
            iCode1 != iCode2
            ):
            d_atom_nos['inter'] += [atom_no2]
        else:
            d_atom_nos['intra'] += [atom_no2]


    ## skip if no *covalent* connections to non-hydrogen atoms
    ## do not ignore intra bonds (possibly disulfide bonds to other space groups)
    if len(d_atom_nos['inter']) == 0:
        return d_CONECT

    if res_name1 not in d_CONECT.keys():
        d_CONECT[res_name1] = {}
    if chain1 not in d_CONECT[res_name1].keys():
        d_CONECT[res_name1][chain1] = {}
    if res_no1 not in d_CONECT[res_name1][chain1].keys():
        d_CONECT[res_name1][chain1][res_no1] = {}
    if iCode1 not in d_CONECT[res_name1][chain1][res_no1].keys():
        d_CONECT[res_name1][chain1][res_no1][iCode1] = {}
    
    d_CONECT[res_name1][chain1][res_no1][iCode1][atom_no1] = d_atom_nos

    return d_CONECT


def parse_recordREMARK(d_header, line, i, lines):

    remark = int(line[6:10])

    if remark == 200:

        experimentaldetail_key = line[12:23].strip().upper()
        if experimentaldetail_key in ['TEMPERATURE','PH']:
            experimentaldetail_value = line[44:].strip()
            if '(' in experimentaldetail_value and ')' in experimentaldetail_value: ## e.g. 1ar4
                experimentaldetail_value = experimentaldetail_value[:experimentaldetail_value.index('(')]+experimentaldetail_value[experimentaldetail_value.index(')')+1:]
                experimentaldetail_value = experimentaldetail_value.strip()
                if experimentaldetail_value == '':
                    print experimentaldetail_value
                    stop
            if 'REMARK200' not in d_header.keys():
                d_header['REMARK200'] = {}
            if experimentaldetail_key not in d_header['REMARK200'].keys():
                d_header['REMARK200'][experimentaldetail_key] = experimentaldetail_value

    elif remark == 290:

        d_header = parse_recordREMARK290(d_header, i, lines)

    elif remark == 350:

        ## biological units
        ## (e.g. 2bq0.pdb, 1thj.pdb, 1m4x.pdb, 1d3i.pdb, 1qgc.pdb, 1rhi.pdb, 1rbo.pdb, 2g8g.pdb, 1h84.pdb)

        d_header = parse_recordREMARK350(d_header, i, lines)

    elif remark == 525: ## water association

        if line[11:].strip() == 'PROTEIN CHAIN  SOLVENT CHAIN':
            for j in range(i+1,len(lines)):
                if lines[j][11:].strip() == '':
                    break
                if lines[j][:10] != 'REMARK 525':
                    break
                else:
                    solventchain = lines[j][11:].split()[1]
                    proteinchain = lines[j][11:].split()[0]
                    if solventchain == proteinchain: ## e.g. 2bq0.pdb
                        continue
                    if 'REMARK525' not in d_header.keys():
                        d_header['REMARK525'] = []
                    d_header['REMARK525'] += [solventchain]

    elif remark == 2: ## resolution
        if line[10:].strip() == '':
            return d_header
        if line.strip() == 'REMARK   2 RESOLUTION. NOT APPLICABLE.':
            return d_header
        if line[23:41] == '   NULL ANGSTROMS.':
            return d_header
        resolution = float(line[23:].replace('ANGSTROMS.',''))
        if resolution == 0:
            print lines[0]
            print 'resolution', resolution
            stop_zero_resolution
        d_header['REMARK2'] = resolution

    elif remark == 0: ## rerefinement
        d_header['REMARK0'] = True

    elif remark == 465: ## missing residues (ATOM)
        d_header = parse_pdb.parse_recordREMARK465(line, lines, i, d_header, 465)

    elif remark == 470: ## missing atoms (ATOM)
        d_header = parse_pdb.parse_recordREMARK470(line, lines, i, d_header, 470)

    elif remark == 475: ## residues (ATOM) with zero occupancy
        d_header = parse_pdb.parse_recordREMARK465(line, lines, i, d_header, 475)

    elif remark == 480: ## polymer atoms (ATOM) with zero occupancy
        d_header = parse_pdb.parse_recordREMARK470(line, lines, i, d_header, 480)

    elif remark == 610: ## non-polymer residues (HETATM) with missing atoms
        pass

    elif remark == 615: ## non-polymer atoms (HETATM) with zero occupancy
        pass

    return d_header


def parse_recordREMARK290(d_header, i, lines):

    line = lines[i]

    if 'REMARK290' not in d_header.keys():
        d_header['REMARK290'] = {}

    if line[13:18] == 'SMTRY':
        operator = int(line[21:23])
        if operator == 0:
            print line
            stop
        dimension = int(line[18])
        if dimension == 1:
            matrixrow1 = lines[i+0][24:].split()
            matrixrow2 = lines[i+1][24:].split()
            matrixrow3 = lines[i+2][24:].split()
            d_header['REMARK290'][operator] = [matrixrow1,matrixrow2,matrixrow3,]

    return d_header


def parse_recordREMARK350(d_header, i, lines):

    line = lines[i]

    if 'REMARK350' not in d_header.keys():
        d_header['REMARK350'] = {}

    if line[11:23] == 'BIOMOLECULE:':
        biomolecules = line[23:80].replace(' ','').split(',') ## multiple biomolecules in e.g. 1wa3
        d_header = loop_and_identify_chains_and_matrices(i, lines, d_header, biomolecules)

    return d_header


def loop_and_identify_chains_and_matrices(i, lines, d_header, biomolecules):

    chains = set()

    method = 'unknown'

    for j in range(i+1,len(lines)):

        if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
            break

        elif 'AUTHOR DETERMINED' in lines[j]:
            method = 'author'

        ## ignore software determined quarternary structures (parse directly from PISA later if needed) - problem with 5r1r,1r1r otherwise
        elif 'SOFTWARE DETERMINED' in lines[j]:
            if len(chains) > 0:
                print lines[j]
                stop
            if not (
                'SOFTWARE USED: ' in lines[j+1]
                or
                'TOTAL BURIED SURFACE AREA:' in lines[j+1]
                or
                'APPLY THE FOLLOWING TO CHAINS:' in lines[j+1]
                ):
                print lines[j]
                print lines[j+1]
                stop
            if 'AUTHOR DETERMINED BIOLOGICAL UNIT:' in lines[j-1]:
                method = 'combined'
                continue
            elif 'AUTHOR PROVIDED BIOLOGICAL UNIT:' in lines[j-1]: ## 2vpl
                method = 'combined'
                continue
            else:
                method = 'software'
                continue
##                    if len(d_header['REMARK350'].keys()) == 0:
##                        ## do not break if only software determined (e.g. 2qkt,2qku; 1dze,1qm8)
##                        continue
##                    else:
##                        ## break if software determined and already author determined (e.g. 1y10,1y11?)
##                        if not 'BIOMOLECULE: ' in lines[j-1]:
##                            print lines[j-1]
##                            stop
##                        break

        elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
            chains = set() ## e.g. 1upp.pdb (vs 8ruc.pdb)
            line_chains = lines[j][41:80]
            chains |= parse_REMARK350_chains(line_chains)

        elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
            line_chains = lines[j][53:80]
            chains |= parse_REMARK350_chains(line_chains)

        elif ',' in lines[j][11:80]:
            if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
                line_chains = lines[j][11:80]
                chains |= parse_REMARK350_chains(line_chains)

        ## count and parse chain transformations
        ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj.pdb)
        elif lines[j][13:19] in ['BIOMT3','SMTRY3']:

            if (
                method == 'software'
                and
                int(lines[j][19:23]) == 1 ## only perform check for first matrix, otherwise break after first matrix if software determined (e.g. 1qkt)
                ):
                bool_break = False
                for prev_bm in d_header['REMARK350'].keys():
                    ## one or more chains already part of another biomolecule?
                    if len(set(d_header['REMARK350'][prev_bm]['chains'].keys())&chains) > 0:
                        bool_break = True
                        break
                if bool_break == True:
                    ## break if software determined and already author determined (e.g. 1y10,1y11?)
                    break ## break loop over lines
                else:
                    ## do not break if only software determined (e.g. 2qkt,2qku; 1dze,1qm8)
                    pass

            matrixno = int(lines[j][19:24])
            ## parse transformation matrix
            matrixrow1 = lines[j-2][24:].split()
            matrixrow2 = lines[j-1][24:].split()
            matrixrow3 = lines[j-0][24:].split()
            matrixrows = [matrixrow1,matrixrow2,matrixrow3,]
##                ## find out whether transformation matrix yields a transformation
##                transformation = False
##                for k in range(3):
##                    ## add a zero translation vector if a translation vector is not given
##                    if len(matrixrows[k]) == 3:
##                        matrixrows[k] += [0.]
##                    if float(matrixrows[k][k]) == 1. and float(matrixrows[k][3]) == 0.:
##                        continue
##                    else:
##                        transformation = True

            ## append transformation matrix to dictionary
            for biomolecule in biomolecules:

                biomolecule = int(biomolecule)

                ## biomolecule
                if biomolecule not in d_header['REMARK350'].keys():
                    d_header['REMARK350'][biomolecule] = {}

                ## biomolecule > method
                d_header['REMARK350'][biomolecule]['method'] = method

                ## biomolecule > matrices
                if 'matrices' not in d_header['REMARK350'][biomolecule].keys():
                    d_header['REMARK350'][biomolecule]['matrices'] = {}
                ## matrices > matrixno > matrix
                d_header['REMARK350'][biomolecule]['matrices'][matrixno] = matrixrows

                ## biomolecule > chains
                if 'chains' not in d_header['REMARK350'][biomolecule].keys():
                    d_header['REMARK350'][biomolecule]['chains'] = {}
                for chain in chains:
                    ## chains > chain
                    if chain not in d_header['REMARK350'][biomolecule]['chains'].keys():
                        d_header['REMARK350'][biomolecule]['chains'][chain] = set()
                    d_header['REMARK350'][biomolecule]['chains'][chain] |= set([matrixno])

    return d_header


def parse_REMARK350_chains(line_chains):

    ## if sentence necessary due to e.g. 1qgc.pdb
    ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: 1 2 3 4 5
    if ',' not in line_chains:
        chains = line_chains.split()
    else:
        ## remove 'AND' from the line of chains (e.g. problem with 1rhi.pdb)
        ## replace '.' in the line of chains (e.g. problem with 1rbo.pdb and 1qgc.pdb)
        chains = line_chains.replace('AND',',').replace('.',',').replace(' ',',').split(',')

    ## loop removal of blank chains necessary due to e.g. 2g8g.pdb
    ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, ,
    for x in range(100):
        if '' in chains:
            chains.remove('')
        else:
            break

    for j in range(len(chains)):
        chain = chains[j]
        if chain == 'NULL':
            chains[j] = ' '

    return set(chains)


def build_ATOMseq(line,lines,i,d_ATOMseq,d_header, d_coordinates,):

    ## I should get rid of this function and parse mmCIF instead!

    record = line[:6].strip()
    altloc = line[16]
    res_name_ATOM = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]

    skip = False

    ## water
    if res_name_ATOM == 'HOH':
        skip = True

    ## not water
    else:
        ## modified residue?
        SEQRESres = check_if_SEQRESres(res_name_ATOM,record,d_header,chain,res_no,iCode,)
        if SEQRESres == False:
            skip = True

        ## peptide chain?
        if chain in d_header['SEQRES']['chains'].keys():
            type = d_header['SEQRES']['chains'][chain]['type']
            if type != 'peptide':
                skip = True
        else:
            skip = True

    if skip == True:
        return d_ATOMseq

    ## initiate new chain in dic
    if not chain in d_ATOMseq.keys():
        d_ATOMseq[chain] = {
            'seq':[],'iCodes':[],'res_nos':[],'altlocs':[],'records':[],
            'indexes':[],'ss':[],
            }

    if lines[i-1][:6].strip() in ['ATOM','HETATM','ANISOU','SIGUIJ','SIGATM',]:
        chain_prev = lines[i-1][21]
        if chain == chain_prev:
            res_no_prev = int(lines[i-1][22:26])
            iCode_prev = lines[i-1][26]
            altloc_prev = lines[i-1][16]
        else: ## e.g. 3bve
            res_no_prev = None
            iCode_prev = None
            altloc_prev = None
    else:
        res_no_prev = None
        iCode_prev = None
        altloc_prev = None

    ## skip if altloc or not the first atom in residue
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
##                        ## REMARK465 before ATOM (first residue > 1 and second residue > 1) ## e.g 2h27
##                        if res_no_prev == None and res_no > 1 and res_no > min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
##                            print res_no
##                            print d_header['REMARK465']['chains'][chain]['residues'].keys()
##                            stop
##                            pass_if = True
                    ## REMARK465 before ATOM ## e.g 3bve
                    if res_no_prev == None and res_no > 1:
                        pass_if = True
                    ## REMARK465 before ATOM
                    ## 1sgf,1nu0
                    if res_no_prev == min(d_header['REMARK465']['chains'][chain]['residues'].keys()):
                        pass_if = True

                if pass_if == True:

                    ## REMARK465 before ATOM?
                    ## e.g. 3bef
##                        if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
##
##                            l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']
##                            index1 = s_alphabet.index(min(l_iCodes_REMARK465))
##                            index2 = s_alphabet.index(max(l_iCodes_REMARK465))+1
##                            l_iCodes_ascending = ','.join(s_alphabet[index1:index2]).split(',')
##                            l_iCodes_descending = list(l_iCodes_ascending)
##                            l_iCodes_descending.reverse()
##                            if l_iCodes_REMARK465 == l_iCodes_ascending:
##                                ascending = True
##                                descending = False
##                            elif l_iCodes_REMARK465 == l_iCodes_descending:
##                                ascending = False
##                                descending = True
##
##                            if len(l_iCodes_REMARK465) > 1 and iCode == ' ' and ascending == True:
##                                if min(d_header['REMARK465']['chains'][chain]['residues'].keys()) == res_no: ## e.g. 3e5v
##                                    l_REMARK465_res_nos = [res_no]
##                                else:
##                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no) ## e.g. 1b8m
##                                stop1
##                            elif len(l_iCodes_REMARK465) > 1 and iCode == ' ' and descending == True:
##                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 2ass
##                                stop2
##                            elif iCode != ' ':
##                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 3bef
##                                stop3
##                            ## e.g. 2bvs
##                            elif len(l_iCodes_REMARK465) == 1:
##                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1)
##                                l_REMARK465_res_names = []
##                                for res_no_REMARK465 in l_REMARK465_res_nos:
##                                    if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
##                                        l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
##                                        for iCode_REMARK465 in l_iCodes_REMARK465:
##                                            res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][' ']['res_name']
##                                            l_REMARK465_res_names += [res_name_REMARK465]
##                                iCode_REMARK465 = l_iCodes_REMARK465[0]
##                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_REMARK465]['altlocs'][' ']['res_name']
##                                SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(l_REMARK465_res_names)]
##                                if SEQRES_seq == l_REMARK465_res_names:
##                                    pass
##                                else:
##                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)
##                                stop4
##                            else:
##                                stop
##                        else:
##                            if min(d_header['REMARK465']['chains'][chain]['residues'].keys()) == res_no:
##                                stop_example
####                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)
##                            l_REMARK465_seq = list(d_header['REMARK465']['chains'][chain]['seq'])
                    l_REMARK465_seq = list(d_header['REMARK465']['chains'][chain]['seq'])

                    ##
                    ## how many residues from REMARK465 records are to be used?
                    ##
                    ## e.g. 1bd7
                    bool_break = False
                    l_REMARK465_res_names = []
                    for i_seq_REMARK465 in range(len(l_REMARK465_seq)):
                        seq_REMARK465 = l_REMARK465_seq[i_seq_REMARK465]
                        res_no_REMARK465 = int(seq_REMARK465[:4])
                        iCode_REMARK465 = seq_REMARK465[4]
                        altloc_REMARK465 = ' ' ## altlocs not used in REMARK465 records
                        res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][altloc_REMARK465]['res_name']
                        bool_res_no_REMARK465_descending = False
                        bool_res_no_REMARK465_ascending = False
                        if i_seq_REMARK465 < len(l_REMARK465_seq)-1:
                            if l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465+1,' ',):
                                bool_res_no_REMARK465_ascending = True
                            if l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465-1,' ',):
                                bool_res_no_REMARK465_descending = True
                        if i_seq_REMARK465 > 0:
                            if l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465-1,' ',):
                                bool_res_no_REMARK465_ascending = True
                            if l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465+1,' ',):
                                bool_res_no_REMARK465_descending = True
                            
                        ## e.g. 1ef0, 3bve, 2iez, many others
                        if res_name_REMARK465 == res_name_ATOM:

##                                if chain == 'A' and res_no == 1 and iCode == ' ' and res_no_REMARK465 != 9999:
##                                    print chain, res_no, iCode
##                                    print 'res_name', res_name_REMARK465, res_name_ATOM
##                                    print 'res_no  ', res_no_REMARK465, res_no
##                                    print 'iCode   ', iCode_REMARK465, iCode
##                                    print
##                                    print '***REM465', l_REMARK465_seq
##                                    print '***ATOM  ', d_coordinates['chains'][chain]['seq']
##                                    print
####                                    stop

                            ##
                            ## not the first ATOM residue
                            ##
                            if len(d_coordinates['chains'][chain]['seq']) > 1:
                                print 'b', 'not the first ATOM residue', chain, res_no, res_no_REMARK465

                                ##
                                ## 1g2w,1l4d
                                if i_seq_REMARK465+1 == len(l_REMARK465_seq):
                                    bool_REMARK465_continuation_forward = 'Unknown'
                                else:
                                    if (
                                        i_seq_REMARK465 != 0 ## 1g2w
                                        and
                                        res_no_REMARK465 == int(l_REMARK465_seq[i_seq_REMARK465+1][:4])-1
                                        ):
                                        bool_REMARK465_continuation_forward = True
                                    elif ( ## 2ven
                                        i_seq_REMARK465 == 0 ## 2ven
                                        and
                                        res_no_REMARK465 > res_no_prev
                                        and
                                        iCode == ' '
                                        and
                                        '%4i%1s' %(res_no-1,' ',) in l_REMARK465_seq
                                        ):
                                        bool_REMARK465_continuation_forward = True
                                    elif ( ## 1o6z
                                        i_seq_REMARK465 != 0
                                        and
                                        bool_res_no_REMARK465_ascending == True
                                        and
                                        '%4i%1s' %(res_no-1,' ',) in l_REMARK465_seq
                                        ):
                                        bool_REMARK465_continuation_forward = True
                                    else:
                                        bool_REMARK465_continuation_forward = False
                                if i_seq_REMARK465 == 0:
                                    bool_REMARK465_continuation_backward = 'Unknown'
                                else:
                                    if res_no_REMARK465 != int(l_REMARK465_seq[i_seq_REMARK465-1][:4])+1:
                                        bool_REMARK465_continuation_backward = True
                                    else:
                                        bool_REMARK465_continuation_backward = False
                                if res_no_REMARK465 in [
                                    res_no-1,res_no_prev+1,
                                    ]:
                                    bool_REMARK465_continuation_ATOM = True
                                else:
                                    bool_REMARK465_continuation_ATOM = False

                                ##
                                ## increasing ATOM res_no
                                if iCode == ' ' and d_coordinates['chains'][chain]['seq'][-2] == '%4i%1s' %(res_no-1,' ',):
                                    ## 1abj, 1sg8
                                    bool_break = True
##                                        s1
                                    break
                                ## increasing ATOM iCode
                                elif res_no == int(d_coordinates['chains'][chain]['seq'][-2][:4]) and d_coordinates['chains'][chain]['seq'][-2] == '%4i%1s' %(res_no,s_alphabet[s_alphabet.index(iCode)-1],):
                                    ## 3ee0
                                    bool_break = True
##                                        s2
                                    break
                                elif not (
                                    ## increasing res_no from ATOM to REMARK465
                                    bool_REMARK465_continuation_ATOM == True
                                    or
                                    bool_REMARK465_continuation_backward == True
                                    or
                                    bool_REMARK465_continuation_forward == True
                                    ):
                                    bool_break = True
                                    print i_seq_REMARK465
                                    print res_no_REMARK465
                                    print res_no_REMARK465, l_REMARK465_seq
                                    print bool_REMARK465_continuation_ATOM
                                    print bool_REMARK465_continuation_backward
                                    print bool_REMARK465_continuation_forward
##                                        stop
##                                        s3
                                    break
                                else:
                                    bool_break = False
                                    pass

                            ##
                            ## the first ATOM residue
                            ##
                            elif len(d_coordinates['chains'][chain]['seq']) == 1:
                                print 'a', 'first ATOM residue', chain, res_no, res_no_REMARK465
##                                    if chain == 'B' and res_no == 161:
##                                        print line
##                                        print res_no_REMARK465, res_no
##                                        stop

                                ##
                                ## break (current residue is ATOM)
                                ##

                                ## res_no ascending
                                if res_no_REMARK465 > res_no and bool_res_no_REMARK465_descending == False:
                                    ## 2cmj, 2ge8, 2z2q, 3bve
                                    bool_break = True
                                    break
                                ## iCode descending
                                elif (
                                    res_no == res_no_REMARK465
                                    and
                                    s_alphabet[s_alphabet.index(iCode_REMARK465)-1] == iCode
                                    and
                                    (iCode != ' ' or i_seq_REMARK465 != 0) ## e.g. 3e5v
                                    ):
##                                        print s_alphabet[s_alphabet.index(iCode_REMARK465)+1]
##                                        print l_REMARK465_seq[i_seq_REMARK465+1]
##                                        print '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)-1],)
##                                        stop
##                                        print '%4i%1s' %(res_no,s_alphabet[s_alphabet.index(iCode)-1],), 'in', l_REMARK465_seq[i_seq_REMARK465:]
##                                        print l_REMARK465_seq[i_seq_REMARK465-1], '==', '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)-1],)
##                                        print i_seq_REMARK465, '==', 0
##                                        stop
                                    ## 1sgi, 3e5v
                                    bool_break = True
                                    break
                                ## iCode ascending
                                elif (
                                    l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no,s_alphabet[s_alphabet.index(iCode)-1],)
                                    ):
                                    ## e.g. 1jy0
                                    bool_break = True
                                    break
##                                    ## prev res_no in REMARK465, next res_no not in REMARK465; WHAT THE HELL????
##                                    elif iCode == ' ' and '%4i%1s' %(res_no-1,' ',) in l_REMARK465_seq and '%4i%1s' %(res_no-1,' ',) not in l_REMARK465_seq:
##                                        bool_break = False
##                                        pass

                                ##
                                ## don't break (current residue is REMARK465)
                                ##
                                
                                ## res_no ascending (skip zero) (don't break)
                                elif (
                                    res_no == 1 and res_no_REMARK465 == -1
                                    and
                                    l_REMARK465_seq[:l_REMARK465_seq.index('%4i%1s' %(-1,' ',))+1] == ['%4i%1s' %(i,' ',) for i in range(int(l_REMARK465_seq[0][:4]),0)]
                                    ):
                                    ## e.g. 1gou, 3csp
                                    bool_break = False
                                    pass
                                ## res_no ascending (no skip zero) (don't break)
                                elif (
                                    ## res_no_prev == res_no-1
                                    (
                                        l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465-1,' ',)
##                                            or
##                                            l_REMARK465_seq[i_seq_REMARK465-1] == -1 and res_no_REMARK465 == 1 ## 1oag
                                        )
                                    and
                                    (
                                        ## terminal
                                        res_no-1 == res_no_REMARK465
                                        or
                                        ## not terminal, res_no_next == res_no+1
                                        (
                                            ## last REMARK465 record
                                            i_seq_REMARK465+1 == len(l_REMARK465_seq) ## e.g. 1e1c
                                            or
                                            l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465+1,' ',)
                                            or
                                            (l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(1,' ',) and res_no_REMARK465 == -1) ## 1oag
                                            )
                                        )
                                    ):
                                    bool_break = False
                                    pass
                                ## res_no descending (don't break)
                                elif bool_res_no_REMARK465_descending == True:
                                    ## e.g. 1fze,3bve
                                    bool_break = False
                                    pass
                                ## iCode descending (don't break)
                                elif (
                                    res_no == res_no_REMARK465
                                    and
                                    ## next REMARK465 = current REMARK465 and lower iCode
                                    l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)-1],)
##                                        and
##                                        ## prev REMARK465 = current REMARK465 and higher iCode
##                                        l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)+1],)
                                    ):
                                    ## e.g. 1fze
                                    bool_break = False
                                    pass
                                ## iCode ascending (don't break)
                                elif (
                                    res_no == res_no_REMARK465
                                    and
                                    (
                                        ## prev expected ATOM in REMARK465
                                        ## e.g. 1sgi
                                        (iCode != ' ' and '%4i%1s' %(res_no,s_alphabet[s_alphabet.index(iCode)-1],) in l_REMARK465_seq[i_seq_REMARK465:])
                                        or
                                        ## next REMARK465 = current REMARK465 and higher iCode
                                        ## e.g. 3e5v
                                        (iCode == ' ' and l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)+1],))
                                        )
                                    and
                                    (
                                        ## prev REMARK465 = current REMARK465 and lower iCode (or current REMARK465 = prev REMARK465 and higher iCode)
                                        ## e.g. 1jy0
                                        l_REMARK465_seq[i_seq_REMARK465-1] == '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)-1],)
                                        or
                                        ## first REMARK465
                                        ## e.g. 1jy0
                                        i_seq_REMARK465 == 0
                                        )
                                    ):
                                    bool_break = False
                                    pass
                                elif i_seq_REMARK465 == 0 and res_no-1 == res_no_REMARK465:
                                    ## e.g. 1oqo
                                    bool_break = False
                                    pass
                                elif res_no_REMARK465 < res_no and '%4i%1s' %(res_no_REMARK465+1,' ',) == l_REMARK465_seq[i_seq_REMARK465+1]:
                                    ## e.g. 3glm
                                    bool_break = False
                                    pass
                                else:
                                    print chain, res_no, iCode
                                    print 'res_name', res_name_REMARK465, res_name_ATOM
                                    print 'res_no  ', res_no_REMARK465, res_no
                                    print 'iCode   ', iCode_REMARK465, iCode
                                    print
                                    print '***REM465', l_REMARK465_seq
                                    print '***ATOM  ', d_coordinates['chains'][chain]['seq']
                                    print
                                    print l_REMARK465_seq[i_seq_REMARK465+1] == '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)-1],)
                                    print l_REMARK465_seq[i_seq_REMARK465-1],  '%4i%1s' %(res_no_REMARK465,s_alphabet[s_alphabet.index(iCode_REMARK465)+1],)
                                    print 
                                    stop_new_case
                            if res_no_REMARK465 > res_no and int(l_REMARK465_seq[i_seq_REMARK465-1][:4]) <= res_no: ## e.g. 2iez,2gp9
                                bool_break = True
                                break
                        
                            ## end if res_name_REMARK465 == res_name_ATOM
                        if l_REMARK465_res_names+[res_name_REMARK465] != d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)+1]:
                            bool_break = True
                            break
                        l_REMARK465_res_names += [res_name_REMARK465]
                        ## end for i_seq_REMARK465 in range(len(l_REMARK465_seq))
                    if bool_break == True:
                        l_REMARK465_seq = l_REMARK465_seq[:i_seq_REMARK465]

##                        if chain == 'A' and res_no == 14 and iCode == ' ':
##                            print 'REM465', l_REMARK465_res_names+[res_name_REMARK465]
##                            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)+1]
##                            print chain, res_no_REMARK465,iCode_REMARK465,altloc_REMARK465
##                            print l_REMARK465_seq
##                            stop2

                    ##
                    ## REMARK465 before ATOM ?
                    ##
                    if len(l_REMARK465_res_names) > 0:
                        l_SEQRES_res_names = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)]
                        ## multiple residue insertion
                        ## 4htc
                        if len(l_REMARK465_res_names) > 1 and l_REMARK465_res_names == l_SEQRES_res_names:
                            REMARK465_before_ATOM = True
                        ## single residue insertion
                        elif (
                            len(l_REMARK465_res_names) == 1 and
                            l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                            res_name_ATOM != l_SEQRES_res_names[0]
                            ):
                            REMARK465_before_ATOM = True
                        ## single residue insertion, res_no_465 < res_no_ATOM
                        ## 2hu9, 1oqo
                        elif (
                            len(l_REMARK465_res_names) == 1 and
                            l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                            res_name_ATOM == l_SEQRES_res_names[0]# and
##                                res_no > min(l_REMARK465_res_nos)
                            ):
                            REMARK465_before_ATOM = True
                        ## REMARK465 not before ATOM
                        ## e.g. 3bef,1bd7,4htc
                        else:
                            REMARK465_before_ATOM = False
                            s4_example_shouoldnt_be_any_if_none_then_delete_this_if_statement
                        REMARK465_before_ATOM = True
                    ## REMARK465 not before ATOM
                    ## e.g. 2a0q, 3ee0
                    else:
                        REMARK465_before_ATOM = False

                    if not REMARK465_before_ATOM == True: ## e.g. 103l
                        if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                            res_no_REMARK465 = res_no
                            l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                            iCode_REMARK465 = l_iCodes_REMARK465[0]
                            res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][' ']['res_name']

        try:
            res_name_SEQRES = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
        except:
            res_name_SEQRES = 'N/A'

        ##
        ## REMARK465 after ATOM
        ##
        if REMARK465_before_ATOM == False and res_name_ATOM == res_name_SEQRES:
            d_ATOMseq,d_header = append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,None,False,True,)
        ##
        ## REMARK465 before ATOM (certain)
        ##
        elif REMARK465_before_ATOM == True:# and res_name_ATOM != res_name_SEQRES:# and res_name_REMARK465 == res_name_SEQRES:
            d_ATOMseq,d_header = append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_seq,True,True,)
        ##
        ##
        ##
        else:
            
            atom_name = line[12:16].strip()
            print '---'
            print chain, res_no, iCode
            SEQRES_res_name = d_header['SEQRES']['chains'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
            SEQRES_seq = d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
            print chain, res_no, iCode
            print line
            print 'ATOM  ', d_ATOMseq[chain]['seq']
            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3']
            print 'maybe', d_header['SEQRES']['chains'][chain]['seq3'][0], chain, res_no, 'is a MODRES and no MODRES record is present?'
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
##                elif res_name_REMARK465 == None and min(d_header['REMARK465']['chains'][chain]['residues'].keys()) < res_no: ## e.g. 3c5w
##                    print chain,res_no
##                    print d_header['REMARK465']['chains'][chain]['residues'].keys()
##                    stop_temp_broken
##                    pass
            ## REMARK465 after ATOM
            elif (
                'REMARK465' in d_header.keys() and
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
                print '*******'
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print 'SEQRES', SEQRES_seq
                print line
                print chain,res_no
                print 'SEQRES', SEQRES_res_name
                print 'ATOM  ', res_name_ATOM
                print 'REMARK', res_name_REMARK465
                print
##                    print l_REMARK465_seq
                print d_ATOMseq[chain]['res_nos']
                stop_N_terminal

        if d_ATOMseq[chain]['seq'] != d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]:
            print '*******'
            print 'ATOM  ', d_ATOMseq[chain]['seq']
            print 'SEQRES', d_header['SEQRES']['chains'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
            print line
            print res_name_ATOM
            print res_name_REMARK465
            print chain, res_no, iCode
            stop_sequence_difference

    return d_ATOMseq        


def check_if_SEQRESres(res_name,record,d_header,chain,res_no,iCode):

    if res_name not in d_res1.keys()+l_nucleotides and record == 'ATOM':
        print res_name,record
        stop

    if res_name in d_res1.keys()+l_nucleotides and record in ['ATOM','REMARK465',]:
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

    ## cases for which MODRES is not used
    if MODRES == False:
        ## functional N-terminal groups
        if res_name in l_terminalmodres and res_name in [d_header['SEQRES']['chains'][chain]['seq3'][0],d_header['SEQRES']['chains'][chain]['seq3'][1],d_header['SEQRES']['chains'][chain]['seq3'][-1],]:
            MODRES = True
        ## l-peptide linkers
        if res_name in l_lpeptidelinkers and res_name in d_header['SEQRES']['chains'][chain]['seq3']:
            MODRES = True
        ## d-amino acids
        if res_name in l_dpeptidelinkers and res_name in d_header['SEQRES']['chains'][chain]['seq3']:
            MODRES = True

    return MODRES


def append_ATOMseq(
    record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_seq,append_REMARK465,append_ATOM
    ):

    if altloc == None:
        altloc = ' '

    if append_REMARK465 == True:
        for i_seq_REMARK465 in range(len(l_REMARK465_seq)):
            seq_REMARK465 = l_REMARK465_seq[i_seq_REMARK465]
            res_no_REMARK465 = int(seq_REMARK465[:4])
            iCode_REMARK465 = seq_REMARK465[4]
            altloc_REMARK465 = ' '
            res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['altlocs'][altloc_REMARK465]['res_name']
            d_ATOMseq[chain]['seq'] += [res_name_REMARK465]
            d_ATOMseq[chain]['res_nos'] += [res_no_REMARK465]
            d_ATOMseq[chain]['iCodes'] += [iCode_REMARK465]
            d_ATOMseq[chain]['altlocs'] += [' ']
            d_ATOMseq[chain]['records'] += ['REMARK465']
            d_ATOMseq[chain]['indexes'] += ['%1s%4s%1s' %(chain,str(res_no).zfill(4),iCode)]
            d_ATOMseq[chain]['ss'] += ['']
            ## remove res_no for each iCode
            d_header['REMARK465']['chains'][chain]['seq'].remove(seq_REMARK465)
        for seq_REMARK465 in l_REMARK465_seq:
            res_no_REMARK465 = int(seq_REMARK465[:4])
            if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                del d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]

    if append_ATOM == True:
        d_ATOMseq[chain]['seq'] += [res_name_ATOM]
        d_ATOMseq[chain]['res_nos'] += [res_no]
        d_ATOMseq[chain]['iCodes'] += [iCode]
        d_ATOMseq[chain]['altlocs'] += [altloc]
        d_ATOMseq[chain]['records'] += [record]
        d_ATOMseq[chain]['indexes'] += ['%1s%4s%1s' %(chain,str(res_no).zfill(4),iCode)]
        d_ATOMseq[chain]['ss'] += ['']

    return d_ATOMseq, d_header


def build_dictionary_of_molecules(d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=True,):

    d_adjacency_backward = {} ## C --> O
    d_adjacency_forward = {} ## O --> C
    matrix_adjacency = []
    l_connections = []
    set_monomers = set()
    d_connections = {}
    l_monomers = []
    for hetID1 in d_CONECT.keys():
        for chain1 in d_CONECT[hetID1].keys():
            for res_no1 in d_CONECT[hetID1][chain1].keys():
                for iCode1 in d_CONECT[hetID1][chain1][res_no1].keys():
                    for atom_no1 in d_CONECT[hetID1][chain1][res_no1][iCode1].keys():
                        
                        atom_name1 = d_atomnos[atom_no1]['atom_name']
                        element1 = d_atomnos[atom_no1]['element']
                        if element1 in l_atoms_metal:
                            continue
##                            if atom_name1 in l_atoms_metal:
##                                continue
##                            if atom_name1[:2] in l_atoms_metal:
##                                continue
##                            if atom_name1[:1] in l_atoms_metal:
##                                print atom_name1
##                                stop
##                                continue

                        ## check ...
                        if s_pdb not in ['1nkm','1nju',]:
                            if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']) > 1:
                                atom_names = set([d_atomnos[atom_no1]['atom_name']])
                                for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                    atom_names |= set([d_atomnos[atom_no]['atom_name']])
                                if atom_names != set(['SG']):
                                    ## one atom of one residue connected to more than one atom of other residue(s)
                                    if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']) == 2:
                                        tmpatom_no1 = d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][0]
                                        tmpatom_no2 = d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][1]
                                        ## if not altloc of same residue
                                        if not (
                                            d_atomnos[tmpatom_no1]['chain'] == d_atomnos[tmpatom_no2]['chain'] and
                                            d_atomnos[tmpatom_no1]['res_no'] == d_atomnos[tmpatom_no2]['res_no'] and
                                            d_atomnos[tmpatom_no1]['iCode'] == d_atomnos[tmpatom_no2]['iCode'] and
                                            d_atomnos[tmpatom_no1]['atom_name'] == d_atomnos[tmpatom_no2]['atom_name'] and
                                            (d_atomnos[tmpatom_no1]['altloc'] == 'A' and d_atomnos[tmpatom_no2]['altloc'] == 'B')
                                            ):
                                            print 'DOUBLE CONNECTION'
                                            print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                            print 'res_names', d_atomnos[tmpatom_no1]['res_name'], d_atomnos[tmpatom_no2]['res_name']
                                            print 'chains', d_atomnos[tmpatom_no1]['chain'], d_atomnos[tmpatom_no2]['chain']
                                            print 'res_nos', d_atomnos[tmpatom_no1]['res_no'], d_atomnos[tmpatom_no2]['res_no']
                                            print 'iCodes', d_atomnos[tmpatom_no1]['iCode'], d_atomnos[tmpatom_no2]['iCode']
                                            print 'atom_names', d_atomnos[tmpatom_no1]['atom_name'], d_atomnos[tmpatom_no2]['atom_name']
                                            print 'altlocs', d_atomnos[tmpatom_no1]['altloc'], d_atomnos[tmpatom_no2]['altloc']
                                            print '*****'
                                            print (d_atomnos[tmpatom_no1]['altloc'] == 'A' and d_atomnos[tmpatom_no2]['altloc'] == 'B')
                                            print hetID1, chain1, res_no1, iCode1, atom_no1
                                            print hetID1, chain1, res_no1, iCode1, atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']
                                            print s_pdb
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                print atom_no, d_atomnos[atom_no]
                                            try:
                                                print cluster, 'cluster'
                                            except:
                                                None
                                            notexpected_different_atllocs_connected_or_valence_exceeded
                                    else:
                                        if (
                                            (not d_atomnos[d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][-1]]['altloc'] == s_alphabet[len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'])])
                                            ):
                                            print atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                            print hetID1, chain1, res_no1, iCode1, atom_no1
                                            atom_names = set([d_atomnos[atom_no1]['atom_name']])
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                atom_names |= set([d_atomnos[atom_no]['atom_name']])
                                            atom_names = atom_names - set(l_atoms_metal)
                                            print hetID2, atom_names
                                            error = True
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                print d_atomnos[atom_no]
                                                if d_atomnos[atom_no]['res_name'] in ['TML','CYO',]:
                                                    error = False
                                            subtract = 0

                                            ## metal connections
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                print atom_no, d_atomnos[atom_no]['element']
                                                if d_atomnos[atom_no]['element'] in l_atoms_metal:
                                                    subtract += 1
                                            if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'])-subtract == 1:
                                                error = False

                                            print s_pdb
                                            if error == True:
                                                print cluster, 'cluster'
                                                notexpected

                        l_hetIDs_long_atom_names = [
                            ## pyranoses
                            'XYP','G6D','SIA',
                            'FCT','DAN',
                            ## furanoses
                            'AHR','HPD','3DR','FUB','AAB',
                            ## dissacharides
                            'DAF','DCB','FXP',
                            ## benzoxazinoids (hydroxamic acid)
                            'HBO', ## DIMBOA from Maize
                            ## (tetra)pyrrole
                            'DBV','PEB','OPP','ZNH','BCL',
                            ## pyrrole ring
                            'YRR',
                            ## phosphate group
                            'FHP','4IP',
                            ## other
                            'DPM','OAS','CSC','780','PED',
                            ## benzene ring
                            'TMM','PLR','785','762','BMZ','AEN',
                            ## p-Coumaric acid (phenyl propanoid)
                            'HC4',
                            ## nucleobases/nucleosides/nucleotides
                            'DA','UMP','PGD','DT','BZG','8OG','ADP','A','TSP','CMP','AMP','FOX','PPU',
                            ## adenine
                            '1MA','MA7',
                            ## guanine
                            'DG','GDP','OMG','GTP','G','SGP',
                            ## cytidine
                            'C','DC','OMC','DOC','CSF','5CM','ME6','CTP',
                            ## uridine
                            '5IU','PSU','OMU','U','UR3','5BU',
##                                        ## dinucleotides
##                                        'NAP',
                            ]
                        ## check 2a
                        if atom_name1[-1] in ['A','B','H',"'"] and len(atom_name1) > 2:
                            if hetID1 not in l_hetIDs_long_atom_names:
                                print 'hetID', hetID1
                                print 'atom_name', atom_name1
                                print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                print 'atom_no', atom_no1
                                print 'res_no', res_no1
                                stop_new_long_atom_name
                            atom_name1_no = atom_name1[1:-1]
                        else:
                            atom_name1_no = atom_name1[1:]

                        for atom_no2 in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                            hetID2 = d_atomnos[atom_no2]['res_name']
                            chain2 = d_atomnos[atom_no2]['chain']
                            res_no2 = d_atomnos[atom_no2]['res_no']
                            iCode2 = d_atomnos[atom_no2]['iCode']

                            atom_name2 = d_atomnos[atom_no2]['atom_name']
                            element2 = d_atomnos[atom_no2]['element']
                            if element2 in l_atoms_metal:
                                continue
##                                if atom_name2 in l_atoms_metal:
##                                    continue
##                                if atom_name2[:2] in l_atoms_metal:
##                                    continue
##                                if atom_name2[:1] in l_atoms_metal:
##                                    print atom_name2
##                                    stop
##                                    continue

                            ## check 2b
                            if atom_name2[-1] in ['A','B','H',"'"] and len(atom_name2) > 2:
                                if hetID2 not in l_hetIDs_long_atom_names:
                                    print chain2, res_no2, iCode2, hetID2, atom_name2, d_CONECT[hetID2][chain2][res_no2][iCode2][atom_no2]
                                    print hetID2
                                    stop_new_long_atom_name
                                atom_name2_no = atom_name2[1:-1]
                            else:
                                atom_name2_no = atom_name2[1:]

                            modres1 = determine_if_modres(d_header, d_coordinates, chain1, res_no1, iCode1, hetID1)
                            if modres1 == True:
                                continue
                            modres2 = determine_if_modres(d_header, d_coordinates, chain2, res_no2, iCode2, hetID2)
                            if modres2 == True:
                                continue

                            if element1 == 'C' and element2 == 'C':
                                ## unique carbon-carbon bonds
                                if (
                                    (hetID1 == 'TYR' and atom_name1 == 'CE1' and hetID2 == 'TRP' and atom_name2 == 'CH2')
                                    or
                                    (hetID2 == 'TYR' and atom_name2 == 'CE1' and hetID1 == 'TRP' and atom_name1 == 'CH2')

                                    or

                                    (hetID2 == 'PHE' and atom_name2 == 'C' and hetID1 == 'PHE' and atom_name1 == 'C')

                                    or

                                    (hetID2 == 'ARG' and atom_name2 == 'C' and hetID1 == 'CH2' and atom_name1 == 'C')
                                    or
                                    (hetID1 == 'ARG' and atom_name1 == 'C' and hetID2 == 'CH2' and atom_name2 == 'C')

                                    or

                                    (hetID2 == 'ALA' and atom_name2 == 'C' and hetID1 == 'BAS' and atom_name1 == 'C1')
                                    or
                                    (hetID1 == 'ALA' and atom_name1 == 'C' and hetID2 == 'BAS' and atom_name2 == 'C1')

                                    or

                                    (hetID1 == 'LOL' and atom_name1 == 'C' and hetID2 == 'ALQ' and atom_name2 == 'CM')
                                    or
                                    (hetID2 == 'LOL' and atom_name2 == 'C' and hetID1 == 'ALQ' and atom_name1 == 'CM')
                                    ):
                                    pass
                                else:
                                    print hetID1, chain1, res_no1, atom_name1
                                    print hetID2, chain2, res_no2, atom_name2
                                    print atom_no1, atom_no2
                                    print s_pdb
                                    try:
                                        print cluster, 'cluster'
                                    except:
                                        None
                                    ## unexpected carbon-carbon
                                    notexpected_carboncarbon_connection

                            if verbose == True and hetID1 != 'CYS' and hetID2 != 'CYS' and atom_name1 != 'SG' and atom_name2 != 'SG':
                                print hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_name1_no, atom_name2_no
                            #####################################
                            ## posttranslational modifications ##
                            #####################################
                            ##
                            ## peptide bond
                            ##
                            if (
                                hetID1 in d_res1.keys()+l_dpeptidelinkers+l_lpeptidelinkers and hetID2 in d_res1.keys()+l_dpeptidelinkers+l_lpeptidelinkers and
                                atom_name1 == 'C' and atom_name2[0] == 'N'
                                ):
                                bond = 'C,N'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (
                                hetID1 in d_res1.keys()+l_dpeptidelinkers+l_lpeptidelinkers and hetID2 in d_res1.keys()+l_dpeptidelinkers+l_lpeptidelinkers and
                                atom_name2 == 'C' and atom_name1[0] == 'N'
                                ):
                                bond = 'C,N'
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                            ##
                            ## cysteine/glutathione disulphide bond
                            ##
                            elif hetID1 in ['CYS','GSH',] and atom_name1[:2] == 'SG' and atom_name2[0] == 'S':
                                bond = 'S,S'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif hetID2 in ['CYS','GSH',] and atom_name2[:2] == 'SG' and atom_name1[0] == 'S':
                                bond = 'S,S'
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                            ##
                            ## lysine, cysteine (e.g. palmitoylation), histidine, tryptophan, tyrosine (e.g. adenyaltion)
                            ##
                            elif (
                                (hetID1 == 'LYS' and atom_name1 == 'NZ'  and atom_name2[0] in ['C','P']) or
                                (hetID1 == 'CYS' and atom_name1 == 'SG'  and atom_name2[0] == 'C') or
                                (hetID1 == 'TRP' and atom_name1 == 'CZ3' and atom_name2[0] == 'O') or
##                                    (hetID1 == 'PRO' and atom_name1 == 'C' and atom_name2[0] == 'N') or
                                (hetID1 == 'HIS' and atom_name1 == 'NE2' and atom_name2[0] in ['C','P']) or
                                (hetID1 == 'HIS' and atom_name1 in ['ND1',] and atom_name2[0] in ['O',]) or
                                (hetID1 == 'GLU' and atom_name1 in ['OE1','OE2',]) or
                                (hetID1 == 'ASP' and atom_name1 in ['OD','OD1','OD2',]) or
                                (hetID1 == 'TYR' and atom_name1 in ['OH',] and atom_name2[0] in ['C','P'])
                                ):
                                bond = atom_name1[0]+','+atom_name2[0]
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (
                                (hetID2 == 'LYS' and atom_name2 == 'NZ'  and atom_name1[0] in ['C','P']) or
                                (hetID2 == 'CYS' and atom_name2 == 'SG'  and atom_name1[0] == 'C') or
                                (hetID2 == 'TRP' and atom_name2 == 'CZ3' and atom_name1[0] == 'O') or
##                                    (hetID2 == 'PRO' and atom_name2 == 'C' and atom_name1[0] == 'N') or
                                (hetID2 == 'HIS' and atom_name2 == 'NE2' and atom_name1[0] in ['C','P']) or
                                (hetID2 == 'HIS' and atom_name2 in ['ND1',] and atom_name1[0] in ['O',]) or
                                (hetID2 == 'GLU' and atom_name2 in ['OE1','OE2',]) or
                                (hetID2 == 'ASP' and atom_name2 in ['OD','OD1','OD2',]) or
                                (hetID2 == 'TYR' and atom_name2 in ['OH',] and atom_name1[0] in ['C','P'])
                                ):
                                bond = atom_name2[0]+','+atom_name1[0]
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2

                            ##
                            ## unique bonds
                            ##
                                
                            ## unique tryptophan peroxidase pi electron porphyrin interaction (make more general solution for small molecules)
                            elif hetID1 == 'TRP' and atom_name1 == 'NE1' and hetID2 == 'PEO' and atom_name2[0] == 'O':
                                bond = 'N'+','+atom_name2[1:]
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif hetID2 == 'TRP' and atom_name2 == 'NE1' and hetID1 == 'PEO' and atom_name1[0] == 'O':
                                bond = 'N'+','+atom_name1[1:]
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (1mk8, cytochrome c peroxidase)
                            elif hetID1 == 'TYR' and atom_name1 == 'CE2' and hetID2 == 'MET' and atom_name2 == 'SD':
                                bond = 'CE1,SD'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif hetID2 == 'TYR' and atom_name2 == 'CE2' and hetID1 == 'MET' and atom_name1 == 'SD':
                                bond = 'CE1,SD'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                            elif (hetID1 == 'TYR' and atom_name1 == 'CE1' and hetID2 == 'TRP' and atom_name2 == 'CH2'):
                                bond = 'CE,CH2'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'TYR' and atom_name2 == 'CE1' and hetID1 == 'TRP' and atom_name1 == 'CH2'):
                                bond = 'CE,CH2'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (1gge, catalase/peroxidase)
                            elif hetID1 == 'TYR' and atom_name1 == 'CB' and hetID2 == 'HIS' and atom_name2 == 'ND1':
                                bond = 'CB,ND1'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif hetID2 == 'TYR' and atom_name2 == 'CB' and hetID1 == 'HIS' and atom_name1 == 'ND1':
                                bond = 'CB,ND1'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (2azc, inhibitor)
                            elif (hetID1 == 'PHE' and atom_name1 == 'C' and hetID2 == 'PHE' and atom_name2 == 'C'):
                                bond = 'C,C'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'PHE' and atom_name2 == 'C' and hetID1 == 'PHE' and atom_name1 == 'C'):
                                bond = 'C,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (1nkm, inhibitor)
                            elif (hetID1 == 'ALA' and atom_name1 == 'C' and hetID2 == 'BAS' and atom_name2 == 'C1'):
                                bond = 'C,C'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'ALA' and atom_name2 == 'C' and hetID1 == 'BAS' and atom_name1 == 'C1'):
                                bond = 'C,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

##                                ## unique (1nkm, inhibitor)
##                                elif (hetID1 == 'SER' and atom_name1 == 'OG' and hetID2 == 'ALA' and atom_name2 == 'C'):
##                                    bond = 'O,C'
##                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
##                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
##                                elif (hetID2 == 'SER' and atom_name2 == 'OG' and hetID1 == 'ALA' and atom_name1 == 'C'):
##                                    bond = 'O,C'
##                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
##                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (1abj, C-terminal alkylation)
                            elif (hetID1 == 'ARG' and atom_name1 == 'C' and hetID2 == 'CH2' and atom_name2 == 'C'):
                                bond = 'C,C'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'ARG' and atom_name2 == 'C' and hetID1 == 'CH2' and atom_name1 == 'C'):
                                bond = 'C,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (1fkn, L-peptide linkers without parent aa IDs)
                            elif (hetID1 == 'LOL' and atom_name1 == 'C' and hetID2 == 'ALQ' and atom_name2 == 'CM'):
                                bond = 'C,C'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'LOL' and atom_name2 == 'C' and hetID1 == 'ALQ' and atom_name1 == 'CM'):
                                bond = 'C,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ## unique (2hg5, GOA linker)
                            elif (hetID1 == 'ASP' and atom_name1 == 'N' and hetID2 == 'GOA' and atom_name2 == 'C1'):
                                bond = 'C,N'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'ASP' and atom_name2 == 'N' and hetID1 == 'GOA' and atom_name1 == 'C1'):
                                bond = 'C,N'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                            elif (hetID1 == 'TYR' and atom_name1 == 'C' and hetID2 == 'GOA' and atom_name2 == 'O2'):
                                bond = 'C,O'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (hetID2 == 'TYR' and atom_name2 == 'C' and hetID1 == 'GOA' and atom_name1 == 'O2'):
                                bond = 'C,O'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1

                            ##
                            ## errors
                            ##
                                
                            ## unexpected carbon-carbon
                            elif atom_name1[0] == 'C' and atom_name2[0] == 'C':
                                print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                expected_conn1
                            ## unexpected oxygen-oxygen
                            elif atom_name1[0] == 'O' and atom_name2[0] == 'O':
                                print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                if hetID1 in ['GAL',] and hetID2 in ['NAG'] and atom_name1 == 'O5' and atom_name2 == 'O4':
                                    fd = open('remediation_oxygen_valence_exceeded.txt','a')
                                    fd.write('%s %s %s %s %s %s %s %s %s %s %s\n' %(
                                        s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2,
                                        ))
                                    fd.close()
                                    pass
                                else:
                                    expected_conn2_oxygenoxygen
                            ## unexpected other
                            elif atom_name1[0] == atom_name2[0]:
                                print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                expected_conn3
                            ##
                            ## O-glycosylation
                            ##
                            elif (
                                hetID1 == 'SER' and atom_name1 == 'OG' or
                                hetID1 == 'THR' and atom_name1 == 'OG1'
                                ):
                                bond = 'O'+','+atom_name2[1:]
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (
                                hetID2 == 'SER' and atom_name2 == 'OG' or
                                hetID2 == 'THR' and atom_name2 == 'OG1'
                                ):
                                bond = 'O'+','+atom_name1[1:]
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                            ##
                            ## N-glycosylation
                            ##
                            elif atom_name1 == 'ND2' and hetID1 == 'ASN':
                                bond = 'N'+','+atom_name2[1:]
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif atom_name2 == 'ND2' and hetID2 == 'ASN':
                                bond = 'N'+','+atom_name1[1:]
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
##                                ##
##                                ## temporary... hopefully... dependent on the remediation team
##                                ## 1-1 connections?
##                                ##
##                                elif hetID1 in ['KDO','KDA','ABU','BAL','DIB','PYB','IMT'] and hetID2 in ['KDO','KDA','ABU','BAL','DIB','PYB','IMT']:
##                                    print atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
##                                    print atom_no2, d_CONECT[hetID2][chain2][res_no2][iCode2][atom_no2]
##                                    if (
##                                        (atom_name2[0] == 'C' and atom_name1[0] in ['O']) or
##                                        (atom_name1[0] == 'C' and atom_name2[0] in ['N'])
##                                        ):
##                                        bond = atom_name1_no+','+atom_name2_no
##                                        monomer1 = chain1+str(res_no1).zfill(4)+iCode1
##                                        monomer2 = chain2+str(res_no2).zfill(4)+iCode2
##                                        print bond, monomer1, monomer2, atom_name1, atom_name2
##                                    elif (
##                                        (atom_name1[0] == 'C' and atom_name2[0] in ['O']) or
##                                        (atom_name2[0] == 'C' and atom_name1[0] in ['N'])
##                                        ):
##                                        bond = atom_name2_no+','+atom_name1_no
##                                        monomer1 = chain2+str(res_no2).zfill(4)+iCode2
##                                        monomer2 = chain1+str(res_no1).zfill(4)+iCode1
##                                    else:
##                                        print hetID1, hetID2, atom_name1, atom_name2
##                                        notexpected
##                                    if monomer1 > monomer2:## and not (hetID1 in ['BAL','DIB'] or hetID2 in ['BAL','DIB']):
##                                        print hetID1, hetID2, monomer1, monomer2
##                                        notexpected
                            ##
                            ## nucleotide phospo(thio)diester bonds (PTR,DT;)
                            ##
                            elif (atom_name1[0] in ['O','S'] and atom_name2 in ['P','P2']) or (atom_name1 == 'O3P' and atom_name2[0] == 'C'):
                                bond = 'O,P'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (atom_name2[0] in ['O','S'] and atom_name1 in ['P','P2']) or (atom_name2 == 'O3P' and atom_name1[0] == 'C'):
                                bond = 'O,P'
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                            ##
                            ## bond (move to bottom after error finding...)
                            ##
                            elif atom_name1 in ['S'] and atom_name2 == 'N':
                                bond = atom_name1[0]+','+atom_name2[0]
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                if monomer1 > monomer2:
                                    print monomer1, monomer2
                                    notexpected
                            elif atom_name2 in ['S'] and atom_name1 == 'N':
                                bond = atom_name2[0]+','+atom_name1[0]
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                if monomer1 > monomer2:
                                    print monomer1, monomer2
                                    notexpected
                            ##
                            ## N-terminal modification (methylation, acetylation, myristoylation, formylation)
                            ##
                            elif (
                                hetID1 in ['ACE','MYR','OHE','OME','CH2','FMT',] and
##                                    chain1 == chain2 and
                                (
                                    res_no2 == 1 or
                                    res_no1 == min(d_coordinates['chains'][chain1]['residues'].keys())
                                    ) and
                                atom_name2 == 'N' and
                                atom_name1[0] == 'C'
                                ):
                                bond = 'N,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                            elif (
                                hetID2 in ['ACE','MYR','OHE','OME','CH2','FMT',] and
##                                    chain1 == chain2 and
                                (
                                    res_no1 == 1 or
                                    res_no2 == min(d_coordinates['chains'][chain2]['residues'].keys())
                                    ) and
                                atom_name1 == 'N' and
                                atom_name2[0] == 'C'
                                ):
                                bond = 'N,C'
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                            ##
                            ## C-terminal modification (amination/amidation)
                            ##
                            elif (
                                hetID1 in ['NH2',] and
                                chain1 == chain2 and
                                res_no1 == max(d_coordinates['chains'][chain1]['residues'].keys()) and
                                res_no2 == max(set(d_coordinates['chains'][chain2]['residues'].keys())-set([max(d_coordinates['chains'][chain2]['residues'].keys())])) and
                                atom_name2 == 'C' and
                                atom_name1[0] == 'N'
                                ):
                                bond = 'N,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                            elif (
                                hetID2 in ['NH2',] and
                                chain1 == chain2 and
                                res_no2 == max(d_coordinates['chains'][chain2]['residues'].keys()) and
                                res_no1 == max(set(d_coordinates['chains'][chain1]['residues'].keys())-set([max(d_coordinates['chains'][chain1]['residues'].keys())])) and
                                atom_name1 == 'C' and
                                atom_name2[0] == 'N'
                                ):
                                bond = 'N,C'
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                            ##
                            ## cyclization
                            ##
                            elif (
                                hetID1 in d_res1.keys() and
                                hetID2 in d_res1.keys() and
                                chain1 == chain2 and
                                res_no1 == max(d_coordinates['chains'][chain1]['residues'].keys()) and
                                res_no2 == min(d_coordinates['chains'][chain2]['residues'].keys()) and
                                atom_name1 == 'C' and
                                atom_name2 == 'N'
                                ):
                                bond = 'N,C'
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                            elif (
                                hetID1 in d_res1.keys() and
                                hetID2 in d_res1.keys() and
                                chain1 == chain2 and
                                res_no2 == max(d_coordinates['chains'][chain2]['residues'].keys()) and
                                res_no1 == min(d_coordinates['chains'][chain1]['residues'].keys()) and
                                atom_name2 == 'C' and
                                atom_name1 == 'N'
                                ):
                                bond = 'N,C'
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                            ## sulphenyl-amide intermediate (e.g. 1oem)
                            elif (
                                hetID1 == 'CYS' and
                                hetID2 == 'SER' and
                                chain1 == chain2 and
                                res_no2-res_no1 == 1 and
                                atom_name1 == 'SG' and
                                atom_name2 == 'N'
                                ):
                                bond = 'S,N'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (
                                hetID2 == 'CYS' and
                                hetID1 == 'SER' and
                                chain1 == chain2 and
                                res_no1-res_no2 == 1 and
                                atom_name2 == 'SG' and
                                atom_name1 == 'N'
                                ):
                                bond = 'S,N'
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                            ## trp-tyr crosslink (e.g. 1mk8)
                            elif (
                                hetID1 == 'TRP' and
                                hetID2 == 'TYR' and
                                chain1 == chain2 and
                                res_no2-res_no1 == 1 and
                                atom_name1 == 'NE1' and
                                atom_name2 == 'CE1'
                                ):
                                bond = 'N,C'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (
                                hetID2 == 'TRP' and
                                hetID1 == 'TYR' and
                                chain1 == chain2 and
                                res_no1-res_no2 == 1 and
                                atom_name2 == 'NE1' and
                                atom_name1 == 'CE1'
                                ):
                                bond = 'N,C'
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                            ## N,C bond
                            elif (
                                chain1 == chain2 and
                                atom_name1 == 'N' and atom_name2[0] == 'C'
                                ):
                                bond = 'N,C'
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                            elif (
                                chain1 == chain2 and
                                atom_name2 == 'N' and atom_name1[0] == 'C'
                                ):
                                bond = 'N,C'
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                            ## error
                            elif (
                                atom_name1_no != '1' and atom_name2_no != '1' and
                                hetID1 not in ['SIA','SLB','XYP','KDO','KDA','KDB',] and hetID2 not in ['SIA','SLB','XYP','KDO','KDA','KDB',]
                                ): ##  and atom_name1[0] == atom_name2[0]
                                print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, iCode1, iCode2
                                print atom_name1, atom_name2, atom_no1, atom_no2
                                print min(d_coordinates['chains'][chain1]['residues'].keys())
                                print min(d_coordinates['chains'][chain2]['residues'].keys())
                                print max(d_coordinates['chains'][chain1]['residues'].keys())
                                print max(d_coordinates['chains'][chain2]['residues'].keys())
                                try:
                                    print 'cluster', cluster
                                except:
                                    None
                                notexpected_bond_not11
                            ##
                            ## glycosyl and peptide bonds
                            ##
                            ## 1,1-glycoside bond (e.g. trehalose in 1do1)
                            elif int(atom_name1_no) == 1 and int(atom_name2_no) == 1:
                                if atom_name1[0] == 'O' and atom_name2[0] == 'C' and hetID2 in ['GLC','BGC',]:
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif atom_name1[0] == 'C' and atom_name2[0] == 'O' and hetID1 in ['GLC','BGC',]:
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                else:
                                    print atom_name1, atom_name2
                                    stopstop
                                print monomer1, monomer2
                                if hetID2 in d_saccharides.keys() and hetID1 not in d_saccharides.keys():
                                    if atom_name2[:2] == 'C1' and atom_name1[0] == 'O':
                                        if monomer1 != chain1+str(res_no1).zfill(4)+iCode1: ## temp!!!
                                            notexpected
                                        monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                        monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                        bond = atom_name2_no+','+atom_name1_no
                                        if verbose == True:
                                            print 'a', bond
                                    else:
                                        print hetID1, hetID2, atom_no2
                                        notexpected
                                elif hetID1 in d_saccharides.keys() and hetID2 not in d_saccharides.keys():
                                    if atom_name1[:2] == 'C1' and atom_name2[0] == 'O':
                                        if monomer1 != chain2+str(res_no2).zfill(4)+iCode2: ## temp!!!
                                            notexpected
                                        monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                        monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                        bond = atom_name1_no+','+atom_name2_no
                                        if verbose == True:
                                            print 'b', bond
                                    else:
                                        print hetID1, hetID2, atom_name1, atom_name2
                                        notexpected
                                elif hetID1 == 'GLC' and hetID2 == 'GLC':
                                    bond = '1,1'
                                    if verbose == True:
                                        print 'c', bond
                                    pass ## trehalose (e.g. 1do1)
                                else:
                                    print s_pdb
                                    print hetID1,hetID2
                                    stop
                                print monomer1, monomer2, bond
                            ## glycosyl 1
                            elif (
                                int(atom_name1_no) != 1 and
                                (
                                    (atom_name1[0] in ['O','N','S',] and atom_name2[:2] == 'C1')
                                    or
                                    (atom_name1[0] == 'C' and atom_name2[:2] in ['O1','N1',])
                                    or
                                    ## sialic acid
                                    (atom_name1[0] == 'O' and atom_name2 == 'C2' and hetID2 in ['SIA','SLB'])
                                    or
                                    (atom_name1[0] == 'C' and atom_name2 == 'O2' and hetID2 in ['SIA','SLB'])
                                    or
                                    (atom_name1 == 'C4B' and atom_name2 == 'O4A' and hetID1 in ['XYP',] and hetID2 in ['XYP',])
                                    or
                                    (atom_name1 in ['C2','O2',] and atom_name2 in ['O4','O8','C4','C8',] and hetID1 in ['KDO','KDA','KDB',] and hetID2 in ['KDA','KDO','KDB',])
                                    )
                                ):
                                bond = atom_name2_no+','+atom_name1_no
                                monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                if verbose == True:
                                    print 'd', bond
                            ## glycosyl 2
                            elif (
                                int(atom_name2_no) != 1 and
                                (
                                    (atom_name2[0] in ['O','N','S',] and atom_name1[:2] == 'C1')
                                    or
                                    (atom_name2[0] == 'C' and atom_name1[:2] in ['O1','N1'])
                                    or
                                    ## sialic acid
                                    (atom_name2[0] == 'O' and atom_name1 == 'C2' and hetID1 in ['SIA','SLB'])
                                    or
                                    (atom_name2[0] == 'C' and atom_name1 == 'O2' and hetID1 in ['SIA','SLB'])
                                    or
                                    (atom_name2 == 'C4B' and atom_name1 == 'O4A' and hetID2 in ['XYP',] and hetID1 in ['XYP',])
                                    or
                                    (atom_name2 in ['C2','O2',] and atom_name1 in ['O4','O8','C4','C8',] and hetID2 in ['KDO','KDA','KDB',] and hetID1 in ['KDA','KDO','KDB',])
                                    )
                                ):
                                bond = atom_name1_no+','+atom_name2_no
                                monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                if verbose == True:
                                    print 'e', bond
                            ## error
                            else:
                                print s_pdb, hetID1, hetID2
                                print chain1, chain2, res_no1, res_no2, iCode1, iCode2
                                print atom_name1, atom_name2, atom_no1, atom_no2
                                print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                notexpected


                            ##
                            ## write adjacency dictionaries if not disulphide bonds
                            ##
                            if atom_name1 != 'SG' and atom_name2 != 'SG':
                                set_monomers |= set([monomer1,monomer2])
                                d_connections[monomer1] = monomer2

                                if not monomer1 in d_adjacency_forward.keys():
                                    d_adjacency_forward[monomer1] = {}
                                if monomer2 in d_adjacency_forward[monomer1]: ## temp!!!
                                    if bond != d_adjacency_forward[monomer1][monomer2]: ## temp!!!
                                        print s_pdb, bond, d_adjacency_forward[monomer1][monomer2]
                                        print 'monomers', monomer1, monomer2
                                        print 'atoms', atom_no1, atom_no2
                                        print d_adjacency_forward
                                        print d_adjacency_backward
                                        print s_pdb
                                        notexpected ## temp!!!
                                d_adjacency_forward[monomer1][monomer2] = bond
                                    
                                if not monomer2 in d_adjacency_backward.keys():
                                    d_adjacency_backward[monomer2] = []
                                if not monomer1 in d_adjacency_backward[monomer2]:
                                    d_adjacency_backward[monomer2] += [monomer1]

                            ## end of loop over atom_no2
                        ## end of loop over atom_no1
        ## end of loop over hetID1


    ## identify roots of saccharides (simplified back trace)
    set_roots = set_monomers-set(d_adjacency_backward.keys())

    ##
    ## trace saccharide tree forward and identify the saccharide molecule
    ## e.g. 1h3u.pdb, 2c4a.pdb
    ##
    ## monomers (combination of chain,resno,iCode) are unique
    ## monomer hetIDs and bonds are not unique
    ##
    ## path of bonds are unique, but path of monomers are used
    ## to avoid traveling down branches multiple times
    ##
    ## that deeply nested dictionaries can be compared was checked
    ## by modifying nested values of the dictionaries a,b and then comparing them
    ##a = {'A86 ': {'bonds': {'1,N': {'monomer': ' 86A', 'bonds': {}}}}, 'A146 ': {'bonds': {'1,N': {'monomer': ' 146A', 'bonds': {}}}}, 'A200 ': {'bonds': {'1,N': {'monomer': ' 200A', 'bonds': {'1,4': {'monomer': ' 200B', 'bonds': {'1,4': {'monomer': ' 200C', 'bonds': {'1,6': {'monomer': ' 200G', 'bonds': {'1,6': {'monomer': ' 200H', 'bonds': {}}, '1,3': {'monomer': ' 200I', 'bonds': {}}}}, '1,3': {'monomer': ' 200D', 'bonds': {'1,2': {'monomer': ' 200E', 'bonds': {'1,2': {'monomer': ' 200F', 'bonds': {}}}}}}}}}}}}}}}
    ##b = {
    ##    'A146 ': {'bonds': {'1,N': {'monomer': ' 146A', 'bonds': {}}}},
    ##    'A86 ': {'bonds': {'1,N': {'monomer': ' 86A', 'bonds': {}}}},
    ##    'A200 ': {'bonds': {'1,N': {'monomer': ' 200A', 'bonds': {'1,4': {'monomer': ' 200B', 'bonds': {'1,4': {'monomer': ' 200C', 'bonds': {
    ##        '1,3': {'monomer': ' 200D', 'bonds': {'1,2': {'monomer': ' 200E', 'bonds': {'1,2': {'monomer': ' 200F', 'bonds': {}}}}}},
    ##        '1,6': {'monomer': ' 200G', 'bonds': {
    ##            '1,3': {'monomer': ' 200I', 'bonds': {}},
    ##            '1,6': {'monomer': ' 200H', 'bonds': {}},
    ##            }},
    ##        }}}}}}}}}
    ##print a == b

    d_molecules = {}
    for root in set_roots:
        path = [root]
        d_molecules[root] = {'bonds':{}}
        d_m_fwd = d_molecules[root]['bonds']
        l_monomers = [root]
        l_branches = []
        monomer1 = root
        while True:
##                print 'xxx', monomer1, d_molecules[root]

            ##
            ## end of branch and no branching points
            ##
            if monomer1 not in d_adjacency_forward.keys() and len(l_branches) == 0:
                break

            ##
            ## end of branch but branching points
            ##
            elif monomer1 not in d_adjacency_forward.keys() and len(l_branches) > 0:

                ## rewind path to previous branching point
                path = path[:path.index(l_branches[-1])+1]
##                    print 'a0path', path, monomer1, l_branches
                ## forward in dictionary to previous branching point
                d_m_fwd = rewind_dictionary(root, path, d_molecules, d_adjacency_forward)

                l_monomers.append(monomer1)

                try:
                    monomers2 = list(set(d_adjacency_forward[path[-1]].keys())-set(l_monomers))
                except:
                    print monomer1, d_adjacency_forward
                    print l_branches
                    print s_pdb
                    print path
                    stop

                ## end of branch
                if len(monomers2) == 0:
                    ## branch == root
                    if len(path) == 1 and len(l_branches) == 1 and path[0] == l_branches[0]:
                        break
                    elif len(path) == 1:
                        print path, l_branches
                        notexpected
                    ## branch != root
                    else:
                        ## remove branch point from list of branch points
                        l_branches.remove(path[-1])
                        ## rewind path to previous branching point
                        path = path[:-1]
                        monomer1 = path[-1]
##                            print 'a1path', path, monomer1, l_branches
                        ## forward in dictionary to previous branching point
                        d_m_fwd = rewind_dictionary(root, path, d_molecules, d_adjacency_forward)

                else:
                    bond = d_adjacency_forward[path[-1]][monomers2[0]]
                    monomer1 = monomers2[0]
                    ## forward path
                    path.append(monomer1)
##                        print 'a2path', path, monomer1, l_branches
                    ## forward dictionary
                    d_m_fwd = d_m_fwd[bond]['bonds']

            ##
            ## monomer(s) forward of the current monomer (monomer1)
            ##
            else:

                l_monomers.append(monomer1)

                ## determine forward monomers that are in branches which have not been traveled
                monomers2 = list(set(d_adjacency_forward[monomer1].keys())-set(l_monomers))

                ## forward monomers traveled and no branches
                if len(monomers2) == 0 and len(l_branches) == 0: ## e.g. 1h3t
                    break
                ## forward monomers traveled and last branch point
                elif len(monomers2) == 0 and len(l_branches) == 1 and [monomer1] == l_branches: ## e.g. 1h3u
                    break
                ## jump backwards towards previous branch point
                elif len(monomers2) == 0 and len(l_branches) > 0: ## e.g. 1h3u.pdb

                    ## rewind path to previous branching point
                    path = path[:path.index(l_branches[-1])+1]
                    monomer1 = path[-1]
                    ## remove branch point from list of branch points
                    l_branches.remove(path[-1])
##                        print 'a3path', path, monomer1, l_branches
                    ## forward in dictionary to previous branching point
                    d_m_fwd = rewind_dictionary(root, path, d_molecules, d_adjacency_forward)
                    continue

                elif len(monomers2) == 0:

                    print 'l_branches', l_branches
                    print 'monomer1', monomer1
                    print 'monomers2', monomers2
                    print d_adjacency_forward
                    print d_molecules
                    print 'path', path
                    print s_pdb
                    notexpected2


                ## translate monomers
                for monomer2 in monomers2:
                    bond = d_adjacency_forward[monomer1][monomer2]
                    monomer = monomertranslation(monomer2,d_coordinates)
                    d_m_fwd[bond] = {'monomer':monomer,'bonds':{}}

                if len(monomers2) > 1:
                    l_branches += [monomer1]

                bond = d_adjacency_forward[monomer1][monomers2[0]]
                monomer1 = monomers2[0]
                
                ## forward path
                path.append(monomer1)
##                    print 'b0path', path, monomer1, l_branches
                ## forward dictionary
                d_m_fwd = d_m_fwd[bond]['bonds']

        root_hetID = monomertranslation(root,d_coordinates)
        d_molecules[root] = {root_hetID:d_molecules[root]}

##        print d_molecules

    return d_molecules


def determine_if_modres(d_header, d_coordinates, chain, res_no, iCode, res_name):

    modres = False
    if chain in d_header['MODRES'].keys():
        if res_no in d_header['MODRES'][chain].keys():
            if iCode in d_header['MODRES'][chain][res_no].keys():
                l_hetIDs = []
                for altloc in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
                    l_hetIDs += [d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']]
                if len(set(l_hetIDs) & d_header['MODRES'][chain][res_no][iCode]) > 0:
                    modres = True
                else:
                    print chain, res_no, iCode, res_name
                    print d_header['MODRES'][chain][res_no][iCode]
                    print l_hetIDs
                    print 
                    print cluster, 'cluster'
                    notexpected

    return modres


def monomertranslation(monomer,d_coordinates):

    chain = monomer[0]
    res_no = int(monomer[1:-1])
    iCode = monomer[-1]
    l_res_names = []
    for altloc in d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'].keys():
        res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['altlocs'][altloc]['res_name']
        if res_name in d_saccharides.keys():
            res_name = d_saccharides[res_name]['stereo']
        l_res_names += [res_name]
    if len(set(l_res_names)) > 1:
        print l_res_names
        stop

    return res_name

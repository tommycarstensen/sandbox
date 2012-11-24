#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2007-2009

import os, sys, numpy, math, urllib2, parse_pdb
s_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz' # 0123456789

class biounit:

    def main(self, pdb, path_pdb, exclude_ligands = False, exclude_altlocs = True):

        ##
        ## parse lines
        ##
        fd = open('%s%s/pdb%s.ent' %(path_pdb,pdb[1:3],pdb,),'r')
        lines = fd.readlines()
        fd.close()

        ##
        ## parse header
        ##
        d_header = parse_pdb.parse_header(lines,)

        ##
        ## parse PISA transformations
        ##
        d_transformations, status = self.parse_pisa_multimers(pdb, d_header,)

        print d_header.keys()
        d_chains = {}
        for assembly in d_transformations.keys():
            d_chains[assembly] = self.write_transformed_coordinates(
                pdb,assembly,d_transformations,lines,d_header,
                exclude_ligands, exclude_altlocs,
                )

        return d_chains


    def write_transformed_coordinates(
        self,pdb,assembly,d_transformations,lines,d_header,
        exclude_ligands, exclude_altlocs,
        ):

        d_lines_PISA, d_coordinates_PISA, d_chains = self.parse_pdb_coordinates(
            pdb, d_transformations, assembly, lines, exclude_altlocs, 
            )

        lines = []
        l_seqres = []
        l_modres = []
        if 'REMARK465' in d_header.keys():
            l_remark465 = [
                'REMARK 465                                                                      \n',
                'REMARK 465 MISSING RESIDUES                                                     \n',
                'REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       \n',
                'REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               \n',
                'REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                \n',
                'REMARK 465                                                                      \n',
                'REMARK 465   M RES C SSSEQI                                                     \n',
                ]
        else:
            l_remark465 = []
        if 'REMARK470' in d_header.keys():
            l_remark470 = [
                'REMARK 470                                                                      \n',
                'REMARK 470 MISSING ATOM                                                         \n',
                'REMARK 470 THE FOLLOWING RESIDUES HAVE MISSING ATOMS(M=MODEL NUMBER;            \n',
                'REMARK 470 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;          \n',
                'REMARK 470 I=INSERTION CODE):                                                   \n',
                'REMARK 470   M RES CSSEQI  ATOMS                                                \n',
                ]
        else:
            l_remark470 = []
        for chain in d_chains.keys():
            ## REMARK465
            if 'REMARK465' in d_header.keys():
                if chain[0] in d_header['REMARK465']['chains'].keys():
                    for res_no in d_header['REMARK465']['chains'][chain[0]]['residues'].keys():
                        for iCode in d_header['REMARK465']['chains'][chain[0]]['residues'][res_no]['l_iCodes']:
                            res_name = d_header['REMARK465']['chains'][chain[0]]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                            s = 'REMARK 465     %3s %1s  %4i%1s                                                    \n' %(
                                res_name,d_chains[chain],res_no,iCode,
                                )
                            if len(s) > 80:
                                print len(s)
                                stop
                            l_remark465 += [s]
            ## REMARK470
            if 'REMARK470' in d_header.keys():
                if chain[0] in d_header['REMARK470']['chains'].keys():
                    for res_no in d_header['REMARK470']['chains'][chain[0]]['residues'].keys():
                        for iCode in d_header['REMARK470']['chains'][chain[0]]['residues'][res_no]['l_iCodes']:
                            res_name = d_header['REMARK470']['chains'][chain[0]]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                            atoms = d_header['REMARK470']['chains'][chain[0]]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys()
                            for i in range(0,len(atoms),7):
                                s = 'REMARK 470     %3s %1s%4i%1s ' %(
                                        res_name,d_chains[chain],res_no,iCode,
                                        )
                                for j in range(i,min(len(atoms),i+7)):
                                    atom = atoms[j]
                                    if len(atom) > 3:
                                        stop
                                    s += ' %3s  ' %(atom)
                                if len(s) > 80:
                                    print len(s)
                                    print s
                                    print len(atoms)
                                    stop
                                s += (80-len(s))*' '+'\n'
                                l_remark470 += [s]
            ## MODRES
            if 'MODRES' in d_header.keys():
                if chain[0] in d_header['MODRES'].keys():
                    for res_no in d_header['MODRES'][chain[0]].keys():
                        for iCode in d_header['MODRES'][chain[0]][res_no].keys():
                            hetID = d_header['MODRES'][chain[0]][res_no][iCode]
                            s = 'MODRES %4s %3s %1s %4i  %3s  STD RES' %(pdb.upper(),hetID,d_chains[chain],res_no,hetID,)
                            s += (80-len(s))*' '+'\n'
                            l_modres += [s]
            ## SEQRES
            if chain[0] in d_header['SEQRES']['chains'].keys():
                l_seq3 = d_header['SEQRES']['chains'][chain[0]]['seq3']
                for i in range(0,len(l_seq3),13):
                    s = 'SEQRES %3i %1s %4i ' %(i/13+1,d_chains[chain],len(l_seq3),)
                    for j in range(i,min(len(l_seq3),i+13)):
                        s += ' %3s' %(l_seq3[j])
                    s += (80-len(s))*' '+'\n'
                    l_seqres += [s]
        lines += l_remark465
        lines += l_remark470
        lines += l_seqres
        lines += l_modres
        l_molecules = d_lines_PISA['polymer'][assembly].keys()
        l_molecules.sort()
        for molecule in l_molecules:
            lines += d_lines_PISA['polymer'][assembly][molecule]
        if exclude_ligands == False and assembly in d_lines_PISA['ligand'].keys():
            l_molecules = d_lines_PISA['ligand'][assembly].keys()
            l_molecules.sort()
            for molecule in l_molecules:
                lines += d_lines_PISA['ligand'][assembly][molecule]
        fd = open('%s_%s.pdb' %(pdb,assembly,),'w')
        fd.writelines(lines)
        fd.close()

        return d_chains


    def parse_pdb_coordinates(
        self, pdb, d_transformations, assembly, lines, exclude_altlocs, 
        ):

        d_output = {
            'polymer':{},
            'ligand':{},
            }
        d_coordinates = {}
        d_chains = {}

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'ATOM':
                (
                    d_output, d_coordinates, d_chains,
                    ) = self.parse_recordATOM(
                        lines, i, d_transformations, d_output, assembly, d_coordinates,
                        d_chains, exclude_altlocs, 
                        )

            elif record == 'HETATM':
                (
                    d_output, d_coordinates, d_chains,
                    ) = self.parse_recordATOM(
                        lines, i, d_transformations, d_output, assembly, d_coordinates,
                        d_chains, exclude_altlocs, 
                        )

            elif record == 'TER':
                (
                    d_output, d_coordinates, d_chains,
                    ) = self.parse_recordATOM(
                        lines, i, d_transformations, d_output, assembly, d_coordinates,
                        d_chains, exclude_altlocs, 
                        )

            ## break if not first model
            elif record == 'MODEL':
                model = int(line[10:14])
                if model != 1:
                    break

##        MODRES = False
##        l_chains = []
##        for i in range(len(lines)):
##            line = lines[i]
##            record = line[:6].strip()
##            if record == 'MODRES':
##                print pdb
##                print line
##                MODRES = True
##                break
####                stop
##
##        if MODRES == True:
##            print d_output['ligand'].keys()
##            print d_coordinates.keys()
##            iCode = ' '
##            altloc = ' '
##            print 'res_name', d_coordinates[1]['HETATM']['A'][90][iCode][altloc]['res_name']
####            print d_coordinates[1]['ATOM']['A'][90]
##            model = 2
##            print d_output['ligand'][1][model]
##            print d_output['polymer'][1][1]
##            stop

        return d_output, d_coordinates, d_chains


    def parse_recordATOM(
        self, lines, i, d_transformations, d_output, assembly, d_coordinates,
        d_chains, exclude_altlocs, 
        ):

        import numpy, decimal

        line = lines[i]

        if line == 'TER\n': ## WHATIF
            return d_output, d_coordinates, d_chains

        res_name = line[17:20].strip()
        if res_name in ['HOH','DOD',]:
            return d_output, d_coordinates, d_chains

        record = line[:6].rstrip()
        atom_no = int(line[6:11])
        atom_no = 0
        atom_name = line[12:16]
        altloc = line[16]
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        if record in ['ATOM','HETATM']:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coordinate = numpy.array([x, y, z])
            occupancy = float(line[54:60])
            bfactor = float(line[60:66])
            element = line[76:78].strip()
            charge = line[78:80].strip()
            ## WHATIF solution
            if altloc not in [' ','A','1',]:
                return d_output, d_coordinates, d_chains
            ## Tommy solution
##            if altloc != ' ':
##                if altloc == 'A' and occupancy == .5:
##                    altloc = ' '
##                    pass
##                elif occupancy > .5:
##                    altloc = ' '
##                    pass
##                elif occupancy <= .5:
##                    return d_output, d_coordinates, d_chains
##                else:
##                    print altloc, occupancy
##                    print line
##                    stop

        pisa_chain_id = '[%s]%s:%i%s' %(res_name,chain.replace(' ','-'),res_no,iCode.replace(' ',''))
        l_chains = d_transformations[assembly]['chains'].keys()

        ## chain not part of biological unit
        if not (
            chain in l_chains or
            pisa_chain_id in l_chains
            ):
            return d_output, d_coordinates, d_chains

        ##
        ## change pisa_chain_id
        ##

        ## ligand
        if pisa_chain_id in l_chains:
            polymer = False
            atom = 'ligand'
            None
        ## polymer
        elif chain in l_chains:
            pisa_chain_id = chain
            polymer = True
            atom = 'polymer'
        ## water
        elif res_name == 'HOH':
            stop
            pisa_chain_id = pisa_chain_id
            return d_output, d_coordinates
        ## not expected
        else:
            print pisa_chain_id, chains, record
            print chain, chains
            print res_name
            stop

        if not assembly in d_output[atom].keys():
            d_output[atom][assembly] = {}

        l_molecules = d_transformations[assembly]['chains'][pisa_chain_id].keys()
        ## sort molecules to retain chain IDs (assumption that lowest molecule number is original position - transformation matrix equals identity matrix)
        l_molecules.sort()
        for molecule in l_molecules:
            if not molecule in d_output[atom][assembly].keys():
                d_output[atom][assembly][molecule] = []
            matrix = d_transformations[assembly]['chains'][pisa_chain_id][molecule]['r']
            vector = d_transformations[assembly]['chains'][pisa_chain_id][molecule]['t']
##            ## skip if asymmetric unit == biological unit (case 2 of 2)
##            if (
##                matrix[0][0] == 1 and matrix[0][1] == 0 and matrix[0][2] == 0 and
##                matrix[1][0] == 0 and matrix[1][1] == 1 and matrix[1][2] == 0 and
##                matrix[2][0] == 0 and matrix[2][1] == 0 and matrix[2][2] == 1 and
##                vector[0] == 0 and vector[1] == 0 and vector[2] == 0
##                ):
##                continue

            if polymer == True:
                if not chain+str(molecule) in d_chains.keys():
                    ## exclude original chains except self
                    l_new_chains = list(s_alphabet)
                    for i_chain in range(len(l_new_chains)-1,-1,-1):
                        if chain != l_new_chains[i_chain] and l_new_chains[i_chain] in d_transformations[assembly]['chains'].keys():
                            del l_new_chains[i_chain]
                    for i_chain in range(l_new_chains.index(chain),len(l_new_chains)+1):
                        new_chain = l_new_chains[i_chain]
                        if not new_chain in d_chains.values():
                            break
                    d_chains[chain+str(molecule)] = new_chain
                new_chain = d_chains[chain+str(molecule)]
            else:
                new_chain = chain

            if record in ['ATOM','HETATM']:
##                for j in range(3):
##                    vector[j] = round(vector[j],3)
                coordinate_transformed = numpy.dot(matrix,coordinate)+vector
                x = coordinate_transformed[0]
                y = coordinate_transformed[1]
                z = coordinate_transformed[2]
                ## use module decimal to achieve correct rounding of float
                x = float(round(decimal.Decimal(str(x)),3))
                y = float(round(decimal.Decimal(str(y)),3))
                z = float(round(decimal.Decimal(str(z)),3))

                line = '%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' %(
                    record.ljust(6), atom_no, atom_name.ljust(4),
                    altloc, res_name.rjust(3), new_chain, res_no, iCode,
                    x, y, z, occupancy, bfactor, element.rjust(2), charge.rjust(2),
                    )

                ## append to d_coordinates
                if not assembly in d_coordinates.keys():
                    d_coordinates[assembly] = {}
                if not record in d_coordinates[assembly].keys():
                    d_coordinates[assembly][record] = {}
                if not new_chain in d_coordinates[assembly][record].keys():
                    d_coordinates[assembly][record][new_chain] = {}
                if not res_no in d_coordinates[assembly][record][new_chain].keys():
                    d_coordinates[assembly][record][new_chain][res_no] = {}
                if not iCode in d_coordinates[assembly][record][new_chain][res_no].keys():
                    d_coordinates[assembly][record][new_chain][res_no][iCode] = {}
                if not altloc in d_coordinates[assembly][record][new_chain][res_no][iCode].keys():
                    d_coordinates[assembly][record][new_chain][res_no][iCode][altloc] = {'res_name':res_name,'atom_names':{}}
                if not atom_name.strip() in d_coordinates[assembly][record][new_chain][res_no][iCode][altloc]['atom_names'].keys():
                    d_coordinates[assembly][record][new_chain][res_no][iCode][altloc]['atom_names'][atom_name.strip()] = []
                d_coordinates[assembly][record][new_chain][res_no][iCode][altloc]['atom_names'][atom_name.strip()] += [coordinate_transformed]

            elif record == 'TER':
                line = '%6s%5i %4s%1s%3s %1s%4i%1s                                                     \n' %(
                    record.ljust(6), atom_no, atom_name.ljust(4),
                    altloc, res_name.rjust(3), new_chain, res_no, iCode,
                    )

            else:
                stop

            ## append line
            d_output[atom][assembly][molecule] += [line]

            if (
                (record == 'ATOM' and lines[i][:6]+lines[i][11:21]+lines[i][22:30]+lines[i][60:78] != line[:6]+lines[i][11:21]+line[22:30]+line[60:78])
                or
                (record == 'HETATM' and lines[i][:6]+lines[i][11:21]+lines[i][26:30]+lines[i][60:78] != line[:6]+lines[i][11:21]+line[26:30]+line[60:78])
                ):
                print lines[i]
                print line
                print lines[i][:6]+lines[i][11:21]+lines[i][22:30]+lines[i][60:]
                print line[:6]+lines[i][11:21]+line[22:30]+line[60:]
                print len(line)
                print len(lines[i])
                stop

        return d_output, d_coordinates, d_chains


    def parse_pisa_multimers(self, pdb, d_header, verbose=True):

        ## <asm_set> (set of stable assemblies) <ser_no> >> <assembly> <id>

        url = urllib2.urlopen('http://www.ebi.ac.uk/msd-srv/prot_int/cgi-bin/multimers.pisa?%s' %(pdb))
        lines = url.readlines()
        assembly = 0
        score = ''
        d_transformations = {}
        for i in range(len(lines)):

            line = lines[i]

            lvl = 2
            s_index1 = '<status>'
            s_index2 = '</status>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
                index2 = line.index(s_index2)
                index1 = line[:index2].rindex(s_index1)+len(s_index1)
                status = line[index1:index2]
                if status not in [
                    'Ok','Overlapping structures','Entry not found',
                    'No symmetry operations',
                    'Broken composition in PA graph',
                    'Improper symop', ## e.g. 8atc
                    'Too big system', ## e.g. 7rsa with D2O
                    ]:
                    print 'status:', status
                    stopstatus
                if verbose == True:
                    print status

##            lvl = 3
##            s_index1 = '<ser_no>'
##            s_index2 = '</ser_no>'
##            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
##                index2 = line.index(s_index2)
##                index1 = line[:index2].rindex(s_index1)+len(s_index1)
##                ser_no = int(line[index1:index2])
##                if ser_no != 1: ## instaed use <score> or diss_energy > 0 and int_energy << 0
##                    break

            lvl = 3
            s_index1 = '<assembly>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1:
                molecule = 0
                score = ''
                assembly += 1
                if assembly in d_transformations.keys():
                    stop
                d_transformations[assembly] = {'chains':{},}

##            lvl = 4
##            s_index1 = '<id>'
##            s_index2 = '</id>'
##            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
##                index2 = line.index(s_index2)
##                index1 = line[:index2].rindex(s_index1)+len(s_index1)
##                id = int(line[index1:index2])

            lvl = 4
            s_index1 = '<score>'
            s_index2 = '</score>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1:
                if '</score>' in line:
                    score = line[line.index(s_index1)+len(s_index1):line.index(s_index2)]
                else:
                    score = ''
                    for j in range(i+1,len(lines)):
                        if '</score>' in lines[j]:
                            break
                        if j > i+1:
                            score += ' '
                        score += lines[j].strip()
                if score == 'This assembly appears to be stable in solution.':
                    continue
                if score == 'Analysis of crystal interfaces has not revealed any strong indications that this assembly may form stable complexes in solution.':
                    del d_transformations[assembly]
                    continue
                if score == 'This assembly falls into a grey region of complexation criteria. It may or may not be stable in solution.':
                    continue
                if score == 'This assembly falls into a grey region of complex formation criteria. It may or may not be stable in solution.':
                    continue
                if score == 'This assembly may be formed from crystallographic considerations, however it does not appear to be stable in solution.':
                    del d_transformations[assembly]
                    break
                else:
                    print assembly, score
                    stop_unknown

            if score == 'Analysis of crystal interfaces has not revealed any strong indications that this assembly may form stable complexes in solution.':
                continue

            lvl = 4
            s_index1 = '<molecule>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1:
                molecule += 1

            lvl = 5
            s_index1 = '<chain_id>'
            s_index2 = '</chain_id>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
                index2 = line.index(s_index2)
                index1 = line[:index2].rindex(s_index1)+len(s_index1)
                chain_id = line[index1:index2]
                if chain_id not in d_transformations[assembly]['chains'].keys():
                    d_transformations[assembly]['chains'][chain_id] = {}
                if molecule in d_transformations[assembly]['chains'][chain_id].keys():
                    notexpected
                d_transformations[assembly]['chains'][chain_id][molecule] = {
                    'r':numpy.zeros((3,3)),
                    't':numpy.array([0.,0.,0.,]),
                    }
                l_tags = [
                    'rxx','rxy','rxz','tx',
                    'ryx','ryy','ryz','ty',
                    'rzx','rzy','rzz','tz',
                    ]
                for j in range(i+1,i+1+12):
                    k = j-i-1
                    s_index1 = '<%s>' %(l_tags[k])
                    s_index2 = '</%s>' %(l_tags[k])
                    index2 = lines[j].index(s_index2)
                    index1 = lines[j][:index2].rindex(s_index1)+len(s_index1)
                    value = float(lines[j][index1:index2])
##                    if value != 0:
##                        value = round(value,7+math.log(abs(value),10))
                    if k%4 in range(3):
                        d_transformations[assembly]['chains'][chain_id][molecule]['r'][k/4][k%4] = value
                    else:
                        d_transformations[assembly]['chains'][chain_id][molecule]['t'][k/4] = value

        ##
        ## no PISA transformations (assume monomeric)
        ##
        if d_transformations == {}:
            d_transformations = self.build_transformations(lines,d_header,)

        return d_transformations, status


    def build_transformations(self,lines,d_header,):

        ##
        ## monomers
        ##
        l_chains = d_header['SEQRES']['chains'].keys()
        l_chains.sort()
        d_transformations = {}
        for assembly in range(1,len(l_chains)+1):
            chain = l_chains[assembly-1]
            d_transformations[assembly] = {'chains':{}}
            molecule = 1
            d_transformations[assembly]['chains'][chain] = {molecule:{
                'r': numpy.array([
                    [ 1.,  0.,  0.],
                    [ 0.,  1.,  0.],
                    [ 0.,  0.,  1.]
                    ]),
                't': numpy.array([ 0.,  0.,  0.])
                }}

##        ##
##        ## multimers
##        ##
##        assembly = 1
##        d_transformations = {
##            assembly:{'chains':{}}}
##        molecule = 1
##        print d_header['SEQRES'].keys()
##        for chain in d_header['SEQRES']['chains'].keys():
##            d_transformations[assembly]['chains'][chain] = {molecule:{
##                'r': numpy.array([
##                    [ 1.,  0.,  0.],
##                    [ 0.,  1.,  0.],
##                    [ 0.,  0.,  1.]
##                    ]),
##                't': numpy.array([ 0.,  0.,  0.])
##                }}

        ##
        ## ligands
        ##
        if 'HET' in d_header.keys():
            for chain in d_header['HET'].keys():
                for res_no in d_header['HET'][chain]:
                    for iCode in d_header['HET'][chain][res_no].keys():
                        for hetID in d_header['HET'][chain][res_no][iCode]:
                            pisa_chain_id = '[%s]%s:%i%s' %(hetID,chain.replace(' ','-'),res_no,iCode.replace(' ',''))
                            d_transformations[assembly]['chains'][pisa_chain_id] = {molecule:{
                                'r': numpy.array([
                                    [ 1.,  0.,  0.],
                                    [ 0.,  1.,  0.],
                                    [ 0.,  0.,  1.]
                                    ]),
                                't': numpy.array([ 0.,  0.,  0.])
                                }}

        return d_transformations


    def __init__(self):

        return


if __name__ == '__main__':

    if '-pdb' not in sys.argv:
        raise 'use -pdb to select a pdb'
    pdb = sys.argv[sys.argv.index('-pdb')+1].replace('.pdb','')
    if len(pdb) != 4:
        raise 'incorrect pdb code'

    if pdb in [
        ## (deoxy)nucleotides
        '2db3','1utd','1f6e','1i8m','1o0k','1rsb','1xf2','1utf','1utv','1zzi',
        ## amino acid chain (with many missing residues)
        '2bki','2fhs','1mow','1ktr',
        ## amino acid ligand
        '1dlk','1t3e','2v5w','1pts','1a37',
        ## magnesium
        '2a68','2a69','1smy','1uik','1xbz',
        ## cadmium
        '2g4h','2gyd','1xz3',
        ## platinum
        '2dqa',
        ## sulfate
        '1nfv','1nf4',
        ## zinc
        '1zfn',
        ## calcium
        '1v3w',
        ## alpha carbon atoms only
        '1lbg','1y1v',
        ## NAG connected to protein
        '1l7f','1l7g','1l7h',
        ## same chain id in multiple biounits (not virus)
        '1b02','1ddk','2bt4','2ipj',
        ## virus (same chain id in multiple assemblies)
        '2b2d','2bbv','2e0z','1h8t','2bpa','1f8v','2g8g','2buk',
        '1m06','2gsy','1laj','1p5w','7msf','1piv','1z7s','1pvc',
        '1ej6',
        ## ribosomal subunit
        '2b64','2b9m','2b9o','2j00','2j01','1ffk','1hnw','2j02',
        '1k73','1k8a','2j03','1hnx','1ibk','1k9m','1m1k','1hnz',
        '1ibl','1kc8','1kd1','1ibm','1m90','1n32','1n33','2f4v',
        '1kqs','1n8r','2ogo','1s72','1nji','1njm','2otj','1q7y',
        '1njp','2otl','1q81','1q82','1q86','1nkw','2v46','1sm1',
        '2v47','2v48','1w2b','2v49','1pnu','1pny','1nwx','1xbp',
        '1nwy','1yhq','1yi2','2uu9','1yij','2uua','1yit','2uub',
        '1qvf','2uuc','1qvg','2uxb','2uxc','1xmo','1yj9','1xmq',
        '1yjn','1xnq','1yjw','1xnr','486d',
        ## status Improper symop
        '1c50','3al1','5csc','8atc','2g32','1h9r','1h9s','2gpm',
        '1l55','1l57','1l58','1l59','1l61','1l62','1l63','2gq4',
        '2gq5','2gq6','2gq7','2ig2','1krl','3tra','1pup',
        ## status Broken composition in PA graph
        '2asc','1md2','1neg','2tci',
        ## status Too big system
        '2d4i','2o01','1lzn','5rsa','6rsa','7rsa',
        ## overlapping chains
        '1t3n',
        ]:
        raise 'problems with %s. ask tommy.' %(pdb)

    path_pdb = '/data/remediated_pdb/'

    instance_biounit = biounit()
    instance_biounit.main(pdb, path_pdb)

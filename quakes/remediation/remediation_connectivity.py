#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2007

import os

class quakes:

    def main(self):

        self.l_pdbs = []
        subdirs = os.listdir(self.path_pdb)
        subdirs.sort()
        for subdir in subdirs:
            files = os.listdir(self.path_pdb+subdir)
            for file in files:
                if file[-2:] == 'gz':
                    ## gunzip
                    if os.path.isfile('%s%s/%s' %(self.path_pdb,subdir,file[:-3])):
                        os.remove('%s%s/%s' %(self.path_pdb,subdir,file[:-3]))
                    os.system('gunzip %s%s/%s' %(self.path_pdb,subdir,file))
                self.l_pdbs += ['%s' %(file[3:7])]

        self.l_pdbs = list(set(self.l_pdbs))

        d_coordinates = {}
        for pdb in self.l_pdbs:
            d_header = self.parse_header(pdb, stop_error = True)
            d_coordinates, d_molecules, d_ATOMseq = self.parse_coordinates(pdb, d_header, verbose = True, parse_molecules = True,)

        return


    def rsync(self):

        import os

        os.system('rsync -a rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/ %s' %(self.path_pdb))

        return


    def res_name2res_symbol(self, res_name):

        ## this function should take into account nucleotides
        ## unfortunately there is a conflict between:
        ## G: glycine and guanosine
        ## A: alanine and adenosine
        ## T: threonine and thymidine
        ## C: cysteine and cytidine
        ## N: aspargine and unknown nucleotide residue

        if res_name in self.d_res.keys():
            symbol = self.d_res[res_name]
        elif res_name in self.l_nucleotides:
            symbol = res_name[-1]
        else:
            symbol = 'X'

        return symbol
    

    def ATOM2seq(self, d_pdb, chain, d_header, res_no_max = None, stop_error = True):

        d_res_nos = {}
        seq = ''
        ATOMrespos = 0
        res_nos = d_pdb['chains'][chain]['residues'].keys()
        res_nos.sort()
        for i in range(len(res_nos)):
            res_no = res_nos[i]
            if res_no_max != None and res_no >= res_no_max:
                continue

            for i in range(len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])):
                iCode = d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][i]
                try:
                    res_name = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                except:
                    print d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
                    stop

                if res_name == 'HOH':
                    continue
                
                ##
                ## do not append to sequence if hetID is not a MODRES
                ## e.g. TRP in 1utv.pdb
                ##
                if chain in d_header['HET'].keys():
                    if res_no in d_header['HET'][chain].keys():
                        if iCode in d_header['HET'][chain][res_no].keys():
                            if res_name != d_header['HET'][chain][res_no][iCode]:
                                if stop_error == True:
                                    print chain, res_no, iCode, res_name, d_header['HET'][chain][res_no][iCode], d_header['MODRES'][chain][res_no][iCode]
                                    notexpected
                            else:
                                if not chain in d_header['MODRES'].keys():
                                    continue
                                elif res_no not in d_header['MODRES'][chain].keys():
                                    continue
                                elif iCode not in d_header['MODRES'][chain][res_no].keys():
                                    continue
                                else:
                                    pass
                                
                d_res_nos[ATOMrespos] = {'res_no':res_no,'iCode':iCode}
                seq += self.res_name2res_symbol(res_name)
                ATOMrespos += 1

        return seq, d_res_nos


    def determine_if_modres(self, d_header, d_pdb, chain, res_no, iCode, res_name):

        modres = False
        if chain in d_header['MODRES'].keys():
            if res_no in d_header['MODRES'][chain].keys():
                if iCode in d_header['MODRES'][chain][res_no].keys():
                    if d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] == d_header['MODRES'][chain][res_no][iCode]:
                        modres = True
                    else:
                        print chain, res_no, iCode, res_name, d_header['MODRES'][chain][res_no][iCode]
                        print 
                        print self.cluster
                        notexpected

        return modres


    def coordinates2ATOMline(self, res_name, chain, res_no, coordinate, iCode, bfactor, atom_name):

        occupancy = bfactor
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        atom_no = 1
        altloc = ''
        charge = ''
        if 'H' in atom_name:
            element = 'H'
        else:
            element = atom_name[0]
        line = [
            '%6s%5i %2s%2s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n'
            %('ATOM'.ljust(6), atom_no, element.rjust(2), atom_name[len(element):].ljust(2), altloc, res_name.ljust(3), chain, res_no, iCode, x, y, z, occupancy, bfactor, element.rjust(2), charge.rjust(2))
            ]

        return line


    def parse_header(self, s_pdb, stop_error = True):

        ##
        ## read lines
        ##
        fd = open('%s%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
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
##                d_conect = self.parse_atom_no_range(d_conect, 'atom', atom_no)

            elif record == 'HETATM': ## section 9
                continue
##                if parse_atoms == False:
##                    if 'REMARK350' in d_header.keys() and 'HET' in d_header.keys():
##                        break
##                    else:
##                        parse_atoms = True
##                atom_no = int(line[6:11])
##                d_conect = self.parse_atom_no_range(d_conect, 'atom', atom_no)

            elif record == 'REMARK': ## section 2
                d_header = self.parse_recordREMARK(d_header, line, i, lines)

            elif record == 'SEQRES': ## section 3
                d_header = self.parse_recordSEQRES(line, d_header)

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
                    d_header['HET'][chain][res_no][iCode] = hetID
                elif d_header['HET'][chain][res_no][iCode] != hetID:
                    if stop_error == True:
                        print d_header['HET'][chain][res_no][iCode], hetID
                        print chain, res_no, iCode
                        print s_pdb, line
                        notexpected

            elif record == 'MODRES':
                hetID = line[12:15].strip()
                chain = line[16]
                res_no = int(line[18:22])
                iCode = line[22]
                res_name = line[24:27].strip()
                txt = line[29:80].strip()
                if hetID in set(self.d_res.keys())-set(['MSE']) and res_name in set(self.d_res.keys())-set(['MSE']):
                    continue
##                if txt not in [] and hetID in set(self.d_res.keys()+self.l_nucleotides)-set(['MSE']):
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
                    d_header['MODRES'][chain][res_no][iCode] = hetID
                elif hetID != d_header['MODRES'][chain][res_no][iCode]:
                    print line, s_pdb
                    stop

            elif record == 'TITLE': ## section 2
                if not 'TITLE' in d_header.keys():
                    d_header['TITLE'] = line[10:].strip()
                else:
                    if d_header['TITLE'][-1] == '-':
                        d_header['TITLE'] += line[10:].strip()
                    else:
                        d_header['TITLE'] += ' '+line[10:].strip()

## #-MER / @-MER
## MULTIMER
## CHAINS A-C / A,B,C
## HOMO- / HETERO-
## OLIGOMER OF # SUBUNITS
##            elif record == 'COMPND': ## section 2
##                if 'BIOLOGICAL_UNIT' in line or 'BIOLOGICAL UNIT' in line:
##                    print lines[i-1], lines[i], lines[i+1]
##                    for s_biounit in self.l_biounits:
##                        if s_biounit in line:
##                            biounit = s_biounit
##                    if biounit == 'N/A':
##                        for s_biounit in self.d_biounits.keys():
##                            if s_biounit in line:
##                                biounit = s_biounit
##                    print self.pdb, biounit

            elif record == 'HEADER':
                d_header['HEADER'] = line[10:50].strip()

            elif record == 'SPRSDE': ## section 2
                sIDcodes = line[31:70].lower().split()
                for sIDcode in sIDcodes:
                    if sIDcode not in d_header.keys():
                        d_header[sIDcode] = {}
                    d_header[sIDcode]['SPRSDE'] = True

            elif record == 'EXPDTA': ## section 2
                methods = line[10:].strip().split(',')[0]
                if methods[:3] == 'NMR':
                    methods = 'NMR'
                elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                    methods = 'X-RAY'
                d_header['EXPDTA'] = methods

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                if spacegroup in self.d_crystalsynonyms.keys():
                    spacegroup = self.d_crystalsynonyms[spacegroup]
                d_header['CRYST1'] = spacegroup

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
        for chain in d_header['SEQRES'].keys():
            if d_header['SEQRES'][chain]['type'] == 'peptide':
                peptidechains += chain
                if len(d_header['SEQRES'][chain]['seq']) > self.min_len_chain:
                    proteinchains += chain
            elif d_header['SEQRES'][chain]['type'] == 'nucleotide':
                nucleotidechains += chain
            elif d_header['SEQRES'][chain]['type'] == 'saccharide':
                saccharidechains += chain

        d_header['proteinchains'] = proteinchains

##        if s_pdb not in ['1ady','1bhj']:
##            chains = d_header['SEQRES'].keys()
##            if 'REMARK350' not in d_header.keys() and len(d_header['SEQRES'].keys()) > 1 and d_header['EXPDTA'] != 'NMR' and len(saccharidechains) != len(d_header['SEQRES'].keys()):
##                ## biounit not specified as text
##                if biounit == 'N/A':
##                    fd = open('unknownbiounit.txt','a')
##                    fd.write('%s %s %s %s %s\n' %(s_pdb, len(proteinchains), len(peptidechains), len(nucleotidechains), len(saccharidechains)))
##                    fd.close()
##                ## unequal number of proteinchains and size of biounit
##                elif self.d_biounits[biounit] != len(chains):
##                    if biounit == 'MONOMER':
##                        d_header['REMARK350'] = {}
##                        for i in range(len(proteinchains)):
##                            chain = proteinchains[i]
##                            d_header['REMARK350'][i+1] = {'chains': [chain]}
##                    elif len(chains) % self.d_biounits[biounit] == 0:
#### add all combinations of chains to remark350 transformations ?! redudant hits if done for twin pdb as well...
##                        print biounit, self.d_biounits[biounit]
##                        print s_pdb
##                        print d_header.keys()
##                        print d_header['SEQRES'].keys()
##                        print proteinchains
##                        print d_header['HET']
##                        stop3
##                    else:
##                        print s_pdb, biounit, chains, proteinchains
##                        print d_header['HET']
##                        stop3b
##                elif self.d_biounits[biounit] == len(chains):
##                    d_header['REMARK350'] = {1:{'chains': [chains]}}
##                else:
##                    print s_pdb, biounit, proteinchains, chains
##                    stop5
####                d_header['biounit'] = 'small' ## small opposed to monomer if multiple different chains
#### count number of similar chains by comparing SEQRESseqs

        return d_header


    def parse_atom_no_range(self, d_conect, record, atom_no):
        
        if not record in d_conect.keys():
            d_conect[record] = [[atom_no,atom_no]]
        elif d_conect[record][-1][1] == atom_no-1:
            d_conect[record][-1][1] = atom_no
        else:
            d_conect[record] += [[atom_no,atom_no]]

        return d_conect


    def transformationmatrix2rotationmatrix_and_translationvector(self,transformationmatrix):

        import numpy

        translationvector = numpy.array(
            [
                float(transformationmatrix[0][3]),
                float(transformationmatrix[1][3]),
                float(transformationmatrix[2][3]),
                ]
            )

        rotationmatrix = numpy.array(
            [
                [
                    float(transformationmatrix[0][0]),
                    float(transformationmatrix[0][1]),
                    float(transformationmatrix[0][2]),
                    ],
                [
                    float(transformationmatrix[1][0]),
                    float(transformationmatrix[1][1]),
                    float(transformationmatrix[1][2]),
                    ],
                [
                    float(transformationmatrix[2][0]),
                    float(transformationmatrix[2][1]),
                    float(transformationmatrix[2][2]),
                    ],
                ]
            )

        return rotationmatrix, translationvector


    def parse_recordSEQRES(self, line, d_header):

        chain = line[11]

        if 'SEQRES' not in d_header:
            d_header['SEQRES'] = {}
        if chain not in d_header['SEQRES'].keys():
            d_header['SEQRES'][chain] = {}
        if not 'type' in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['type'] = 'N/A'

        l_residues = line[19:70].split()

        s_residues = ''
        for i in range(len(l_residues)):
            residue = l_residues[i]
            if residue in self.d_res.keys():
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'peptide'
                elif d_header['SEQRES'][chain]['type'] != 'peptide': ## e.g. 1vq6
                    d_header['SEQRES'][chain]['type'] = 'peptide'
                s_residues += self.d_res[residue]
            elif residue in self.l_nucleotides:
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'nucleotide'
                elif d_header['SEQRES'][chain]['type'] != 'nucleotide':
                    stop
                s_residues += residue
            elif residue in ['GLC','GAL','MAN','FRU']:
                if d_header['SEQRES'][chain]['type'] == 'N/A':
                    d_header['SEQRES'][chain]['type'] = 'saccharide'
                elif d_header['SEQRES'][chain]['type'] != 'saccharide':
                    stop
                s_residues += residue
            else:
                if residue == 'UNK': ## e.g. 1pny.pdb
                    if d_header['SEQRES'][chain]['type'] == 'N/A':
                        d_header['SEQRES'][chain]['type'] = 'peptide'
                s_residues += 'X'

        if 'seq' not in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['seq'] = ''
        d_header['SEQRES'][chain]['seq'] += s_residues

        if 'seq3' not in d_header['SEQRES'][chain].keys():
            d_header['SEQRES'][chain]['seq3'] = []
        d_header['SEQRES'][chain]['seq3'] += l_residues

        return d_header


    def parse_coordinates(
        self, s_pdb, d_header, verbose = False, parse_molecules = True,
        ):

        print 'parsing coordinates of', s_pdb

        ##
        ## read lines
        ##
        fd = open('%s%s/pdb%s.ent' %(self.path_pdb, s_pdb.lower()[1:3], s_pdb.lower()),'r')
        lines = fd.readlines()
        fd.close()

        ##
        ## set dictionaries
        ##
        d_atomnos = {}
        d_CONECT = {}
        d_coordinates = {}
        d_ATOMseq = {}

        ##
        ## loop over lines
        ##
        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'ATOM':
                d_coordinates, d_line = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
                d_ATOMseq = self.build_ATOMseq(line, lines, i, d_ATOMseq, d_header,)
                d_atomnos[d_line['atom_no']] = d_line

            elif record == 'HETATM':
                res_name = line[17:20].strip()
                ## water
                if res_name in ['D2O','H2O',]:
                    print s_pdb, res_name
                    stop
                ## water must be parsed in case there is a connection to it (e.g. 2bpb)
                if res_name in ['HOH','DOD',]: ## DOD in 2d4j
                    d_coordinates, d_line, = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
                ## MSE
                elif res_name in self.d_modres.keys():
                    d_coordinates, d_line, = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'ATOM',)
                ## (poly)saccharide or other hetero compound
                else:
                    d_coordinates, d_line, = self.parse_recordATOM(line, d_coordinates, lines, i, d_header, 'HETATM',)
                atom_no = d_line['atom_no']
                d_atomnos[atom_no] = d_line

                d_ATOMseq = self.build_ATOMseq(line, lines, i, d_ATOMseq, d_header,) ## 1omw

            elif record == 'CONECT':
                d_CONECT = self.parse_recordCONECT(line, d_atomnos, d_CONECT)
                
            elif record == 'HET':
                hetID = line[7:10].strip()
                if 'HET' not in d_coordinates.keys():
                    d_coordinates['HET'] = set()
                d_coordinates['HET'] |= set([hetID])
##                chain = line[12]
##                res_no = int(line[13:17])
##                iCode = line[17]
##                if iCode != ' ':
##                    print s_pdb
##                    stop
##                n_atoms = int(line[20:25])
##                description = line[30:70]
##                if chain not in d_hetero['chains'].keys():
##                    d_hetero['chains'][chain] = {}
##                if 'residues' not in d_hetero['chains'][chain].keys():
##                    d_hetero['chains'][chain]['residues'] = {}
##                if res_no not in d_hetero['chains'][chain]['residues'].keys():
##                    d_hetero['chains'][chain]['residues'][res_no] = {}
##                if 'iCodes' not in d_hetero['chains'][chain]['residues'][res_no].keys():
##                    d_hetero['chains'][chain]['residues'][res_no]['iCodes'] = {}
##                d_hetero['chains'][chain]['residues'][res_no]['iCodes'][iCode] = hetID
                
            elif record == 'MODEL':
                model = int(line.split()[1])


        if parse_molecules == True:
            d_molecules = self.build_dictionary_of_molecules(d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=verbose,)
        else:
            d_molecules = {}

        for chain in d_header['SEQRES']:

            ## skip if not peptide chain
            type = d_header['SEQRES'][chain]['type']
            if type != 'peptide':
                continue

            ## skip if all residues are unknown
            if len(d_header['SEQRES'][chain]['seq3'])*['UNK'] == d_header['SEQRES'][chain]['seq3']:
                continue

            ## append remaining REMARK465 residues
            if 'REMARK465' in d_header.keys():
                if chain in d_header['REMARK465']['chains'].keys():
                    if len(d_header['REMARK465']['chains'][chain]['residues'].keys()) > 0:
                        l_REMARK465_res_nos = d_header['REMARK465']['chains'][chain]['residues'].keys() ## e.g. 1cd0
                        l_REMARK465_res_nos.sort()
                        d_ATOMseq,d_header = self.append_ATOMseq(None,None,None,None,None,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,False,)

            if d_ATOMseq[chain]['seq'] != d_header['SEQRES'][chain]['seq3']:
                print 'SEQRES', d_header['SEQRES'][chain]['seq3']
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print chain
                print 'SEQRES', len(d_header['SEQRES'][chain]['seq3'])
                print 'ATOM  ', len(d_ATOMseq[chain]['seq'])
                print pdb,chain
                stop_different_length_SEQRES_ATOM

        return d_coordinates, d_molecules, d_ATOMseq


    def append_ATOMseq(
        self,record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,append_REMARK465,append_ATOM
        ):

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
                        d_ATOMseq[chain]['indexes'] += ['%1s%4s%1s' %(chain,str(res_no).zfill(4),iCode)]
                    del d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]
    ##            else: ## not 3b95
    ##                break

        if append_ATOM == True:
            d_ATOMseq[chain]['seq'] += [res_name_ATOM]
            d_ATOMseq[chain]['res_nos'] += [res_no]
            d_ATOMseq[chain]['iCodes'] += [iCode]
            d_ATOMseq[chain]['altlocs'] += [altloc]
            d_ATOMseq[chain]['records'] += [record]
            d_ATOMseq[chain]['indexes'] += ['%1s%4s%1s' %(chain,str(res_no).zfill(4),iCode)]

        return d_ATOMseq, d_header


    def build_ATOMseq(self,line,lines,i,d_ATOMseq,d_header):

        record = line[:6].strip()
        altloc = line[16]
        res_name_ATOM = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]

        skip = False

        ## modified residue?
        SEQRESres = self.check_if_SEQRESres(res_name_ATOM,record,d_header,chain,res_no,iCode,)
        if SEQRESres == False:
            skip = True

        ## peptide chain?
        if chain in d_header['SEQRES']:
            type = d_header['SEQRES'][chain]['type']
            if type != 'peptide':
                skip = True
        else:
            skip = True

        if skip == True:
            return d_ATOMseq

        if not chain in d_ATOMseq.keys():
            d_ATOMseq[chain] = {
                'seq':[],'iCodes':[],'res_nos':[],'altlocs':[],'records':[],'indexes':[],
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

                        ## REMARK465 before ATOM?
                        ## e.g. 3bef
                        if min(d_header['REMARK465']['chains'][chain]['residues'].keys()) <= res_no:
                            if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():

                                l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']
                                index1 = self.s_alphabet.index(min(l_iCodes_REMARK465))
                                index2 = self.s_alphabet.index(max(l_iCodes_REMARK465))+1
                                l_iCodes_ascending = ','.join(self.s_alphabet[index1:index2]).split(',')
                                l_iCodes_descending = list(l_iCodes_ascending)
                                l_iCodes_descending.reverse()
                                if l_iCodes_REMARK465 == l_iCodes_ascending:
                                    ascending = True
                                    descending = False
                                elif l_iCodes_REMARK465 == l_iCodes_descending:
                                    ascending = False
                                    descending = True

                                if len(l_iCodes_REMARK465) > 1 and iCode == ' ' and ascending == True:
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no) ## e.g. 1b8m
                                elif len(l_iCodes_REMARK465) > 1 and iCode == ' ' and descending == True:
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 2ass
                                elif iCode != ' ':
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1) ## e.g. 3bef
                                ## e.g. 2bvs
                                elif len(l_iCodes_REMARK465) == 1:
                                    l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no+1)
                                    l_REMARK465_res_names = []
                                    for res_no_REMARK465 in l_REMARK465_res_nos:
                                        if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                            l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                            for iCode_REMARK465 in l_iCodes_REMARK465:
                                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                                                l_REMARK465_res_names += [res_name_REMARK465]
                                    iCode_REMARK465 = l_iCodes_REMARK465[0]
                                    res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_REMARK465]['res_name']
                                    SEQRES_seq = d_header['SEQRES'][chain]['seq3'][:len(l_REMARK465_res_names)]
                                    if SEQRES_seq == l_REMARK465_res_names:
                                        pass
                                    else:
                                        l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)
                                else:
                                    stop
                            else:
                                l_REMARK465_res_nos = range(min(d_header['REMARK465']['chains'][chain]['residues'].keys()),res_no)

                            ## e.g. 1bd7
                            l_REMARK465_res_names = []
                            for res_no_REMARK465 in l_REMARK465_res_nos:
                                if res_no_REMARK465 in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                    l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                    for iCode_REMARK465 in l_iCodes_REMARK465:
                                        res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']
                                        l_REMARK465_res_names += [res_name_REMARK465]
    ##                            else: ## not 3b95
    ##                                break

                            if len(l_REMARK465_res_names) > 0:
                                l_SEQRES_res_names = d_header['SEQRES'][chain]['seq3'][len(d_ATOMseq[chain]['seq']):][:len(l_REMARK465_res_names)]
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
                                ## 2hu9
                                elif (
                                    len(l_REMARK465_res_names) == 1 and
                                    l_REMARK465_res_names[0] == l_SEQRES_res_names[0] and
                                    res_name_ATOM == l_SEQRES_res_names[0] and
                                    res_no > min(l_REMARK465_res_nos)
                                    ):
                                    REMARK465_before_ATOM = True
                                ## REMARK465 not before ATOM
                                ## e.g. 3bef,1bd7,4htc
                                else:
                                    REMARK465_before_ATOM = False
                            ## REMARK465 not before ATOM
                            ## e.g. 2a0q
                            else:
                                REMARK465_before_ATOM = False

                        if not REMARK465_before_ATOM == True: ## e.g. 103l
                            if res_no in d_header['REMARK465']['chains'][chain]['residues'].keys():
                                res_no_REMARK465 = res_no
                                l_iCodes_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['l_iCodes']
                                iCode_REMARK465 = l_iCodes_REMARK465[0]
                                res_name_REMARK465 = d_header['REMARK465']['chains'][chain]['residues'][res_no_REMARK465]['d_iCodes'][iCode_REMARK465]['res_name']

            try:
                res_name_SEQRES = d_header['SEQRES'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
            except:
                res_name_SEQRES = 'N/A'

            if REMARK465_before_ATOM == False and res_name_ATOM != res_name_SEQRES and res_name_ATOM in ['FOR','HOA','ACE',] and (res_no == min(d_coordinates['chains'][chain]['residues'].keys()) or res_no == max(d_coordinates['chains'][chain]['residues'].keys())):
                fd = open('formyl.txt','a')
                fd.write('%s %s %s\n' %(pdb, chain, res_name_ATOM))
                fd.close()
                return d_coordinates, {
                    'chain':chain,
                    'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
                    'atom_no':atom_no,'atom_name':atom_name,'element':element,
                    }, d_ATOMseq

            ## REMARK465 after ATOM
            if REMARK465_before_ATOM == False and res_name_ATOM == res_name_SEQRES:# and res_name_REMARK465 != res_name_SEQRES:
                d_ATOMseq,d_header = self.append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,None,False,True,)
            ## REMARK465 before ATOM (certain)
            elif REMARK465_before_ATOM == True:# and res_name_ATOM != res_name_SEQRES:# and res_name_REMARK465 == res_name_SEQRES:
                d_ATOMseq,d_header = self.append_ATOMseq(record,res_no,iCode,altloc,res_name_ATOM,d_ATOMseq,chain,d_header,l_REMARK465_res_nos,True,True,)
            else:
                atom_name = line[12:16].strip()
                print '---'
                SEQRES_res_name = d_header['SEQRES'][chain]['seq3'][len(d_ATOMseq[chain]['seq'])]
                SEQRES_seq = d_header['SEQRES'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print chain, res_no, iCode
                print line
                print d_ATOMseq[chain]['seq']
                print d_header['SEQRES'][chain]['seq3']
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
                elif res_name_REMARK465 == None and min(d_header['REMARK465']['chains'][chain]['residues'].keys()) < res_no: ## e.g. 3c5w
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
                        return d_coordinates, {
                            'chain':chain,
                            'res_name':res_name_ATOM,'res_no':res_no,'iCode':iCode,'altloc':altloc,
                            'atom_no':atom_no,'atom_name':atom_name,'element':element,
                            }, d_ATOMseq
                    print '*******'
                    print 'ATOM  ', d_ATOMseq[chain]['seq']
                    print 'SEQRES', SEQRES_seq
                    print line
                    print chain,res_no
                    print 'SEQRES', SEQRES_res_name
                    print 'ATOM  ', res_name_ATOM, '***iCode***', iCode
                    print 'REMARK', res_name_REMARK465, iCode_REMARK465
                    stop_N_terminal

            if d_ATOMseq[chain]['seq'] != d_header['SEQRES'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]:
                print '*******'
                print 'ATOM  ', d_ATOMseq[chain]['seq']
                print 'SEQRES', d_header['SEQRES'][chain]['seq3'][:len(d_ATOMseq[chain]['seq'])]
                print res_name_ATOM
                print line
                print res_name_REMARK465
                print pdb, chain, res_no, iCode
                stop_sequence_difference

        return d_ATOMseq        


    def check_if_SEQRESres(self,res_name,record,d_header,chain,res_no,iCode):

        if res_name not in self.d_res.keys()+self.l_nucleotides and record == 'ATOM':
            print res_name,record
            stop

        if res_name in self.d_res.keys()+self.l_nucleotides and record in ['ATOM','REMARK465',]:
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


    def build_dictionary_of_molecules(self,d_CONECT,d_atomnos,d_header,d_coordinates,s_pdb,verbose=True,):

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
                            if element1 in self.l_atoms_metal:
                                continue
##                            if atom_name1 in self.l_atoms_metal:
##                                continue
##                            if atom_name1[:2] in self.l_atoms_metal:
##                                continue
##                            if atom_name1[:1] in self.l_atoms_metal:
##                                print atom_name1
##                                stop
##                                continue

                            ## check ...                            
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
                                            print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                            print 'res_name', d_atomnos[tmpatom_no1]['res_name'], d_atomnos[tmpatom_no2]['res_name']
                                            print 'chain', d_atomnos[tmpatom_no1]['chain'], d_atomnos[tmpatom_no2]['chain']
                                            print 'res_no', d_atomnos[tmpatom_no1]['res_no'], d_atomnos[tmpatom_no2]['res_no']
                                            print 'iCode', d_atomnos[tmpatom_no1]['iCode'], d_atomnos[tmpatom_no2]['iCode']
                                            print 'atom_name', d_atomnos[tmpatom_no1]['atom_name'], d_atomnos[tmpatom_no2]['atom_name']
                                            print 'altloc', d_atomnos[tmpatom_no1]['altloc'], d_atomnos[tmpatom_no2]['altloc']
                                            print '*****'
                                            print (d_atomnos[tmpatom_no1]['altloc'] == 'A' and d_atomnos[tmpatom_no2]['altloc'] == 'B')
                                            print hetID1, chain1, res_no1, iCode1, atom_no1
                                            print d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']
                                            print s_pdb
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                print atom_no, d_atomnos[atom_no]
                                            notexpected_different_atllocs_connected
                                    else:
                                        if (
                                            (not d_atomnos[d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'][-1]]['altloc'] == self.s_alphabet[len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'])])
                                            ):
                                            print atom_no1, d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]
                                            print hetID1, chain1, res_no1, iCode1, atom_no1
                                            atom_names = set([d_atomnos[atom_no1]['atom_name']])
                                            for atom_no in d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter']:
                                                atom_names |= set([d_atomnos[atom_no]['atom_name']])
                                            atom_names = atom_names - set(self.l_atoms_metal)
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
                                                if d_atomnos[atom_no]['element'] in self.l_atoms_metal:
                                                    subtract += 1
                                            if len(d_CONECT[hetID1][chain1][res_no1][iCode1][atom_no1]['inter'])-subtract == 1:
                                                error = False

                                            print s_pdb
                                            if error == True:
                                                print self.cluster
                                                notexpected

                            l_hetIDs_long_atom_names = [
                                ## pyranoses
                                'XYP','G6D','SIA',
                                'FCT',
                                ## furanoses
                                'AHR','HPD',
                                ## dissacharides
                                'DAF',
                                ## benzoxazinoids (hydroxamic acid)
                                'HBO', ## DIMBOA from Maize
                                ## (tetra)pyrrole
                                'DBV','PEB','OPP','ZNH',
                                ## other
                                'DPM',
                                ## p-Coumaric acid (phenyl propanoid)
                                'HC4',
                                ## nucleobases/nucleosides/nucleotides
                                'DA','UMP','PGD','DT','BZG','8OG','ADP','A','TSP','CMP',
                                ## guanine
                                'DG','GDP','OMG','GTP',
                                ## cytidine
                                'C','DC','OMC',
                                ## uridine
                                '5IU',
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
                                    stop
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
                                if element2 in self.l_atoms_metal:
                                    continue
##                                if atom_name2 in self.l_atoms_metal:
##                                    continue
##                                if atom_name2[:2] in self.l_atoms_metal:
##                                    continue
##                                if atom_name2[:1] in self.l_atoms_metal:
##                                    print atom_name2
##                                    stop
##                                    continue

                                ## check 2b
                                if atom_name2[-1] in ['A','B','H',"'"] and len(atom_name2) > 2:
                                    if hetID2 not in l_hetIDs_long_atom_names:
                                        print chain2, res_no2, iCode2, hetID2, atom_name2, d_CONECT[hetID2][chain2][res_no2][iCode2][atom_no2]
                                        print hetID2
                                        stop
                                    atom_name2_no = atom_name2[1:-1]
                                else:
                                    atom_name2_no = atom_name2[1:]

                                modres1 = self.determine_if_modres(d_header, d_coordinates, chain1, res_no1, iCode1, hetID1)
                                if modres1 == True:
                                    continue
                                modres2 = self.determine_if_modres(d_header, d_coordinates, chain2, res_no2, iCode2, hetID2)
                                if modres2 == True:
                                    continue

                                if element1 == 'C' and element2 == 'C':
                                    print hetID1, chain1, res_no1, atom_name1
                                    print hetID2, chain2, res_no2, atom_name2
                                    print atom_no1, atom_no2
                                    print s_pdb
                                    print self.cluster
                                    notexpected

                                if verbose == True and hetID1 != 'CYS' and hetID2 != 'CYS' and atom_name1 != 'SG' and atom_name2 != 'SG':
                                    print hetID1, hetID2, atom_name1, atom_name2, atom_name1_no, atom_name2_no
                                #####################################
                                ## posttranslational modifications ##
                                #####################################
                                ##
                                ## peptide bond
                                ##
                                if (
                                    hetID1 in self.d_res.keys() and hetID2 in self.d_res.keys() and
                                    atom_name1 == 'C' and atom_name2 == 'N'
                                    ):
                                    bond = 'C,N'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif (
                                    hetID1 in self.d_res.keys() and hetID2 in self.d_res.keys() and
                                    atom_name2 == 'C' and atom_name1 == 'N'
                                    ):
                                    bond = 'C,N'
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                ##
                                ## cysteine disulphide bond
                                ##
                                if hetID1 == 'CYS' and atom_name1 == 'SG' and atom_name2[0] == 'S':
                                    bond = 'S,S'
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'CYS' and atom_name2 == 'SG' and atom_name1[0] == 'S':
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
                                ## tryptophan peroxidase pi electron porphyrin interaction (make more general solution for small molecules)
                                ##
                                elif hetID1 == 'TRP' and atom_name1 == 'NE1' and hetID2 == 'PEO' and atom_name2[0] == 'O':
                                    bond = 'N'+','+atom_name2[1:]
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                elif hetID2 == 'TRP' and atom_name2 == 'NE1' and hetID1 == 'PEO' and atom_name1[0] == 'O':
                                    bond = 'N'+','+atom_name1[1:]
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                ##
                                ## error
                                ##
                                elif atom_name1[0] == 'C' and atom_name2[0] == 'C':
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected1
                                elif atom_name1[0] == 'O' and atom_name2[0] == 'O':
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected2
                                elif atom_name1[0] == atom_name2[0]:
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print atom_name1_no, atom_name2_no, len(set([hetID1,hetID2]) & set(['GOL']))
                                    expected3
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
                                    hetID1 in self.d_res.keys() and
                                    hetID2 in self.d_res.keys() and
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
                                    hetID1 in self.d_res.keys() and
                                    hetID2 in self.d_res.keys() and
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
                                    hetID1 not in ['SIA','SLB','XYP',] and hetID2 not in ['SIA','SLB','XYP',]
                                    ): ##  and atom_name1[0] == atom_name2[0]
                                    print s_pdb, hetID1, hetID2, chain1, chain2, res_no1, res_no2, iCode1, iCode2, atom_name1, atom_name2, atom_no1, atom_no2
                                    print min(d_coordinates['chains'][chain1]['residues'].keys())
                                    print min(d_coordinates['chains'][chain2]['residues'].keys())
                                    print max(d_coordinates['chains'][chain1]['residues'].keys())
                                    print max(d_coordinates['chains'][chain2]['residues'].keys())
                                    print self.cluster
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
                                    if hetID2 in self.d_saccharides.keys() and hetID1 not in self.d_saccharides.keys():
                                        if atom_name2[:2] == 'C1' and atom_name1[0] == 'O':
                                            if monomer1 != chain1+str(res_no1).zfill(4)+iCode1: ## temp!!!
                                                notexpected
                                            monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                            monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                            bond = atom_name2_no+','+atom_name1_no
                                            print 'a', bond
                                        else:
                                            print hetID1, hetID2, atom_no2
                                            notexpected
                                    elif hetID1 in self.d_saccharides.keys() and hetID2 not in self.d_saccharides.keys():
                                        if atom_name1[:2] == 'C1' and atom_name2[0] == 'O':
                                            if monomer1 != chain2+str(res_no2).zfill(4)+iCode2: ## temp!!!
                                                notexpected
                                            monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                            monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                            bond = atom_name1_no+','+atom_name2_no
                                            print 'b', bond
                                        else:
                                            print hetID1, hetID2, atom_name1, atom_name2
                                            notexpected
                                    elif hetID1 == 'GLC' and hetID2 == 'GLC':
                                        bond = '1,1'
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
                                        (atom_name1[0] in ['O','N'] and atom_name2[:2] == 'C1') or
                                        (atom_name1[0] == 'C' and atom_name2[:2] in ['O1','N1']) or
                                        ## sialic acid
                                        (atom_name1[0] == 'O' and atom_name2 == 'C2' and hetID2 in ['SIA','SLB']) or
                                        (atom_name1[0] == 'C' and atom_name2 == 'O2' and hetID2 in ['SIA','SLB']) or
                                        (atom_name1 == 'C4B' and atom_name2 == 'O4A' and hetID1 in ['XYP',] and hetID2 in ['XYP',])
                                        )
                                    ):
                                    bond = atom_name2_no+','+atom_name1_no
                                    monomer1 = chain1+str(res_no1).zfill(4)+iCode1
                                    monomer2 = chain2+str(res_no2).zfill(4)+iCode2
                                    print 'd', bond
                                ## glycosyl 2
                                elif (
                                    int(atom_name2_no) != 1 and
                                    (
                                        (atom_name2[0] in ['O','N'] and atom_name1[:2] == 'C1') or
                                        (atom_name2[0] == 'C' and atom_name1[:2] in ['O1','N1']) or
                                        ## sialic acid
                                        (atom_name2[0] == 'O' and atom_name1 == 'C2' and hetID1 in ['SIA','SLB']) or
                                        (atom_name2[0] == 'C' and atom_name1 == 'O2' and hetID1 in ['SIA','SLB']) or
                                        (atom_name2 == 'C4B' and atom_name1 == 'O4A' and hetID2 in ['XYP',] and hetID1 in ['XYP',])
                                        )
                                    ):
                                    bond = atom_name1_no+','+atom_name2_no
                                    monomer1 = chain2+str(res_no2).zfill(4)+iCode2
                                    monomer2 = chain1+str(res_no1).zfill(4)+iCode1
                                    print 'e', bond
                                ## error
                                else:
                                    print s_pdb, hetID1, hetID2, res_no1, res_no2, atom_name1, atom_name2, atom_no1, atom_no2
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

##        print d_adjacency_forward
##        print d_adjacency_backward
                                    
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
                    d_m_fwd = self.rewind_dictionary(root, path, d_molecules, d_adjacency_forward)

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
                            d_m_fwd = self.rewind_dictionary(root, path, d_molecules, d_adjacency_forward)

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
                        d_m_fwd = self.rewind_dictionary(root, path, d_molecules, d_adjacency_forward)
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
                        monomer = self.monomertranslation(monomer2,d_coordinates)
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

            root_hetID = self.monomertranslation(root,d_coordinates)
            d_molecules[root] = {root_hetID:d_molecules[root]}

##        print d_molecules

        return d_molecules


    def rewind_dictionary(self, root, path, d_molecules, d_adjacency_forward):

        ## forward in dictionary to previous branching point
        d_m_fwd = d_molecules[root]['bonds']
        for i in range(1,len(path)):
            m1 = path[i-1]
            m2 = path[i]
            bond = d_adjacency_forward[m1][m2]
            d_m_fwd = d_m_fwd[bond]['bonds']
        
        return d_m_fwd


    def monomertranslation(self,monomer,d_coordinates):

        chain = monomer[0]
        res_no = int(monomer[1:-1])
        iCode = monomer[-1]
        res_name = d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
        if res_name in self.d_saccharides.keys():
            res_name = self.d_saccharides[res_name]['stereo']

        return res_name


    def parse_recordCONECT(self,line,d_atomnos,d_CONECT):

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

        ## skip if hydrogen atom
        if atom_name1[0] == 'H':
            return d_CONECT
        ## skip if a cofactor to which molecules are not *covalently* bound
        if res_name1 in self.l_cofactors+['HOH']+self.l_solutes:
            return d_CONECT

        d_atom_nos = {'intra':[],'inter':[]}

        for atom_no2 in atom_nos[1:]:
            chain2 = d_atomnos[atom_no2]['chain']
            res_no2 = d_atomnos[atom_no2]['res_no']
            iCode2 = d_atomnos[atom_no2]['iCode']
            res_name2 = d_atomnos[atom_no2]['res_name']
            atom_name2 = d_atomnos[atom_no2]['atom_name']
            ## skip if hydrogen atom
            if atom_name2[0] == 'H':
                continue
            ## skip if a cofactor to which molecules are not *covalently* bound
            if res_name2 in self.l_cofactors+['HOH']+self.l_solutes:
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


    def parse_recordREMARK(self, d_header, line, i, lines):

        remark = int(line[6:10])

        if remark == 200:

            if line[12:23].strip().upper() in ['TEMPERATURE','PH']:
                experimentaldetail_key = line[12:23].strip()
                experimentaldetail_value = line[44:].strip()
                if 'REMARK200' not in d_header.keys():
                    d_header['REMARK200'] = {}
                if experimentaldetail_key not in d_header['REMARK200'].keys():
                    d_header['REMARK200'][experimentaldetail_key] = experimentaldetail_value

        elif remark == 290:

            d_header = self.parse_recordREMARK290(d_header, i, lines)

        elif remark == 350:

            ## biological units
            ## (e.g. 2bq0.pdb, 1thj.pdb, 1m4x.pdb, 1d3i.pdb, 1qgc.pdb, 1rhi.pdb, 1rbo.pdb, 2g8g.pdb, 1h84.pdb)

            d_header = self.parse_recordREMARK350(d_header, i, lines)

        elif remark == 465: ## missing residues

            d_header = self.parse_recordREMARK465(line, d_header, lines, i)

        elif remark == 470: ## missing atoms

            d_header = self.parse_recordREMARK470(line, d_header, lines, i)

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
            try:
                resolution = float(line[22:27])
            except:
                resolution = 'N/A'
            d_header['REMARK2'] = resolution

        return d_header

    def parse_recordREMARK290(self, d_header, i, lines):

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


    def parse_recordREMARK350(self, d_header, i, lines):

        line = lines[i]

        if 'REMARK350' not in d_header.keys():
            d_header['REMARK350'] = {}

        if line[11:23] == 'BIOMOLECULE:':
            biomolecules = line[23:80].replace(' ','').split(',')
            d_header = self.loop_and_identify_chains_and_matrices(i, lines, d_header, biomolecules)

        return d_header


    def parse_REMARK350_chains(self, line_chains):

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


    def loop_and_identify_chains_and_matrices(self, i, lines, d_header, biomolecules):

        chains = set()

        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
                break

            elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                chains = set() ## e.g. 1upp.pdb (vs 8ruc.pdb)
                line_chains = lines[j][41:80]
                chains |= self.parse_REMARK350_chains(line_chains)

            elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
                line_chains = lines[j][53:80]
                chains |= self.parse_REMARK350_chains(line_chains)

            elif ',' in lines[j][11:80]:
                if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
                    line_chains = lines[j][11:80]
                    chains |= self.parse_REMARK350_chains(line_chains)

            ## count and parse chain transformations
            ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj.pdb)
            elif lines[j][13:19] in ['BIOMT3','SMTRY3']:

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


    def parse_recordREMARK465(self, line, d_pdb, lines, i):

        ## missing residues

        if line[10:].strip() in ['M RES C SSSEQI','M RES C  SSEQI']:

            for j in range(i+1,len(lines)):

                if lines[j][:10] != 'REMARK 465':
                    break

                try:
                    model = int(lines[j][12:14])
                except:
                    model = 'N/A'
                res_name = lines[j][15:18].strip()
                chain = lines[j][19]
                res_no = int(lines[j][22:26])
                iCode = lines[j][26]


                if not 'REMARK465' in d_pdb.keys():
                    d_pdb['REMARK465'] = {}
                if not 'chains' in d_pdb['REMARK465'].keys():
                    d_pdb['REMARK465']['chains'] = {}
                if not chain in d_pdb['REMARK465']['chains'].keys():
                    d_pdb['REMARK465']['chains'][chain] = {}
                if not 'residues' in d_pdb['REMARK465']['chains'][chain].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['REMARK465']['chains'][chain]['residues'].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no] = {}

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['REMARK465']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}


                ## iCode > res_name
                d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name

                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['REMARK465']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if not iCode in d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes']:
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## iCode > REMARK
                if not 'REMARK' in d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['REMARK'] = True

                ## iCode > record
                d_pdb['REMARK465']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = 'REMARK465'

        return d_pdb


    def parse_recordREMARK470(self, line, d_pdb, lines, i):

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
                if not 'REMARK470' in d_pdb.keys():
                    d_pdb['REMARK470'] = {}
                if not 'chains' in d_pdb.keys():
                    d_pdb['REMARK470']['chains'] = {}
                if not chain in d_pdb['REMARK470']['chains'].keys():
                    d_pdb['REMARK470']['chains'][chain] = {}
                if not 'residues' in d_pdb['REMARK470']['chains'][chain].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'] = {}
                if not res_no in d_pdb['REMARK470']['chains'][chain]['residues'].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no] = {}

                ## res_no > l_iCodes
                if not 'l_iCodes' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] = []
                ## l_iCodes > iCode
                if len(d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes']) == 0:
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
                elif iCode not in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes']:
                    ## e.g. 2fs4 (chain A, res_no 162, iCode " ")
                    if iCode == ' ':
                        d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] = [iCode]+d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']
                    else:
                        d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]

                ## res_no > d_iCodes
                if not 'd_iCodes' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
                ## d_iCodes > iCode
                if not iCode in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes']:
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

                ## iCode > atoms
                if not 'atoms' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
                ## atoms > atom_name > coordinate
                for atom_name in atoms:
                    if not atom_name in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                        d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {}
                    ## iCode > REMARK
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name]['REMARK'] = True

                ## iCode > res_name
                if not 'res_name' in d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                    d_pdb['REMARK470']['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name


        elif line[10:].strip().split() == ['M','RES','CSSEQI','ATOMS']:
            print self.pdb1
            print self.pdb2
            notexpected

        return d_pdb


    def parse_recordATOM(self, line, d_pdb, lines, i, d_header, record,):

        import numpy

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
        element = line[76:78].strip()
        coordinate = numpy.array([x, y, z])

        if not 'chains' in d_pdb.keys():
            d_pdb['chains'] = {}
        if not chain in d_pdb['chains'].keys():
            d_pdb['chains'][chain] = {}
        if not 'residues' in d_pdb['chains'][chain].keys():
            d_pdb['chains'][chain]['residues'] = {}

        ## res_no
        if not res_no in d_pdb['chains'][chain]['residues'].keys():
            d_pdb['chains'][chain]['residues'][res_no] = {}

        ## res_no > d_iCodes
        if not 'd_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'] = {}
        ## d_iCodes > iCode
        if not iCode in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes']:
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}


        res_name_conflict = False
        ## iCode > res_name
        if not 'res_name' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
            d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] = res_name
##        ## temp!!! temporary until pdb errors are fixed (otherwise split ATOM and HETATMs)
##        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name: ## temp!!!
##            print res_name, d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] ## temp!!!
##            print chain, res_no, iCode, altloc ## temp!!!
##            print line
##            stop_duplicate_remark465_and_coord_section
        ## check that res_name is correct (e.g. 2fes:L:1)
        elif d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name'] != res_name and altloc == ' ': ## 1fh2:A:30 altloc
            if record == 'HETATM':
                res_name_conflict = True
            else:
                ## change the iCode
                iCode_max = max(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
                iCode_max = self.s_alphabet[self.s_alphabet.index(iCode_max)+1]
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_max] = d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]
                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'].index(iCode)] = iCode_max

                ## d_iCodes > iCode
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode] = {}

        if res_name_conflict == False:
            ## res_no > l_iCodes
            if not 'l_iCodes' in d_pdb['chains'][chain]['residues'][res_no].keys():
                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = []
            ## l_iCodes > iCode
            if len(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']) == 0:
                d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
            elif iCode not in d_pdb['chains'][chain]['residues'][res_no]['l_iCodes']:
                d_pdb = self.identify_iCode_sequence(d_pdb, chain, res_no, iCode, res_name, d_header)

            ## iCode > atoms
            if not 'atoms' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'] = {}
            ## atoms > atom_name > coordinate
            if not atom_name in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'].keys():
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms'][atom_name] = {'coordinate':coordinate}

            ## iCode > record
            if not 'record' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode].keys():
                d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['record'] = record

        
        return d_pdb, {
            'chain':chain,
            'res_name':res_name,'res_no':res_no,'iCode':iCode,'altloc':altloc,
            'atom_no':atom_no,'atom_name':atom_name,'element':element,
            }


    def identify_iCode_sequence(self, d_pdb, chain, res_no, iCode, res_name, d_header):
        
        l_iCodes = list(d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'])
        d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] += [iCode]
        for iCode_prev in l_iCodes:
            index_alphabet = self.s_alphabet.index(iCode)
            if (
                ## REMARK465
                'REMARK' in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev].keys() or
                ## REMARK470
                {'REMARK':True} in d_pdb['chains'][chain]['residues'][res_no]['d_iCodes'][iCode_prev]['atoms'].values()
                ):
                ATOMseq,d_res_nos = self.ATOM2seq(d_pdb, chain, d_header, res_no_max = res_no)
                res_symbol = self.res_name2res_symbol(res_name)
                ## 1) REMARK residues before ATOM residues (N-terminal)
                ## e.g. 1jqz.pdb
                if (
                    len(ATOMseq) == 0
                    ):
                    None
                ## 2) REMARK465 residues after ATOM residues
                ## e.g. 1nuo.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)]
                    ):
                    l_iCodes = [iCode]+d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'][:-1]
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 3) REMARK470 residues after ATOM residues
                ## e.g. 2lve.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' in l_iCodes
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 4) REMARK470 residues after ATOM residues
                ## e.g. 2j5q (chain B, res_no 54, iCode C)
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet-1] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)] and
                    self.s_alphabet[index_alphabet-1] in l_iCodes and
                    ' ' not in l_iCodes 
                    ):
                    index = l_iCodes.index(self.s_alphabet[index_alphabet-1])
                    l_iCodes = l_iCodes[:index+1]+[iCode]+l_iCodes[index+1:]
                    list(set(self.s_alphabet)-set(l_iCodes))
                    d_pdb['chains'][chain]['residues'][res_no]['l_iCodes'] = l_iCodes
                ## 5) REMARK465 residues before ATOM residues
                ## e.g. 1uij.pdb
                elif (
                    len(ATOMseq) > 0 and
                    res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+len(l_iCodes)] and
                    ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)]
                    ):
                    None
                else: ## e.g. 1fne
                    print '*****'
                    print iCode, iCode_prev, res_symbol, index_alphabet, l_iCodes
                    print d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet]
                    print
                    print len(ATOMseq) > 0
                    print res_symbol == d_header['SEQRES'][chain]['seq'][len(ATOMseq)+index_alphabet-1]
                    print ATOMseq == d_header['SEQRES'][chain]['seq'][:len(ATOMseq)]
                    print self.s_alphabet[index_alphabet-1] in l_iCodes
                    print 
                    print 1, ATOMseq, self.res_name2res_symbol(res_name)
                    print 2, d_header['SEQRES'][chain]['seq']
                    print chain, res_no, iCode, iCode_prev, res_name
                    expected
                break

        return d_pdb


    def spherermsd(
        self,
        pdb1, pdb2, ## pdbs
        d_header, ## sequences
        d_pdb, ## coordinates
        l_equivalent_chains, ## equivalent chains
        d_chains_interpdb_sequence_similar, ## mutations
        tv1, rm, tv2, ## transformation
        ):

        import sys
        sys.path.append('/home/people/tc/python/Protool/')
        import geometry
        instance_geometry = geometry.geometry()

        for rep_chain1 in d_chains_interpdb_sequence_similar.keys():
            rep_chain2 = d_chains_interpdb_sequence_similar[rep_chain1]['rep_chain2']

            ATOMseq1,d_res_nos1 = self.ATOM2seq(d_pdb[pdb1], rep_chain1, d_header[pdb1])
            ATOMseq2,d_res_nos2 = self.ATOM2seq(d_pdb[pdb2], rep_chain2, d_header[pdb2])
            l1 = d_chains_interpdb_sequence_similar[rep_chain1]['l1']
            l2 = d_chains_interpdb_sequence_similar[rep_chain1]['l2']

            if l1 > 0 or l2 > 0: ## e.g. 1ftg.pdb,1dx9.pdb
                print ATOMseq1
                print ATOMseq2
                s1 = d_chains_interpdb_sequence_similar[rep_chain1]['s1']
                s2 = d_chains_interpdb_sequence_similar[rep_chain1]['s2']
                print s1
                print s2
                print pdb1, pdb2
                expected

            (
                coordinates1, coordinates2, rescount,
                ) = self.ATOMrecords2coordinates(
                    d_pdb, pdb1, pdb2, rep_chain1, rep_chain2, d_res_nos1, d_res_nos2,
                    l1, l2, len(ATOMseq1), d_ATOMseq, tv1=tv1, rm=rm, tv2=tv2
                    )

            l_mutations = d_chains_interpdb_sequence_similar[rep_chain1]['l_mutations']
            for mutation in l_mutations:
                res_no1 = d_res_nos1[mutation[0]]['res_no']
                res_no2 = d_res_nos2[mutation[0]]['res_no']
                iCode1 = d_res_nos1[mutation[0]]['iCode']
                iCode2 = d_res_nos2[mutation[0]]['iCode']
                hypocenter1 = d_pdb[pdb1]['chains'][rep_chain1]['residues'][res_no1]['d_iCodes'][iCode1]['atoms']['CA']['coordinate']
                hypocenter2 = d_pdb[pdb2]['chains'][rep_chain2]['residues'][res_no2]['d_iCodes'][iCode2]['atoms']['CA']['coordinate']

                d_rmsds = {4:0,8:0,16:0,32:0}
                for dist in d_rmsds.keys():

                    sqdist = dist**2
                    sphere_coordinates1 = []
                    sphere_coordinates2 = []

                    for i in range(len(coordinates1)):

                        coordinate1 = coordinates1[i]
                        coordinate2 = coordinates2[i]

                        sqdist1 = sum((hypocenter1-coordinate1)**2)
                        sqdist2 = sum((hypocenter2-coordinate2)**2)

                        if min(sqdist1,sqdist2) < sqdist:

                            sphere_coordinates1 += [coordinate1]
                            sphere_coordinates2 += [coordinate2]

                    rmsd = instance_geometry.superpose(
                        sphere_coordinates1,sphere_coordinates2
                        )

                    print rep_chain1, mutation, dist, rmsd, float(len(sphere_coordinates1))/len(coordinates1)

                    d_rmsds[dist] = rmsd

        rmsd4 = d_rmsds[4]
        rmsd8 = d_rmsds[8]
        rmsd16 = d_rmsds[16]
        rmsd32 = d_rmsds[32]

        return rmsd4, rmsd8, rmsd16, rmsd32


    def __init__(self):

        import os

        self.d_res = {
            'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','SER':'S','THR':'T','CYS':'C','MET':'M','ASP':'D','GLU':'E','ASN':'N','GLN':'Q','HIS':'H','LYS':'K','ARG':'R','PRO':'P','PHE':'F','TRP':'W','TYR':'Y',
            'MSE':'M','UNK':'X','ASX':'X','GLX':'X',
            }

        self.l_nucleotides = [
             'A', 'C', 'G', 'U',      'I', ## ribonucleotides
            'DA','DC','DG',     'DT','DI', ## deoxyribonucleotides
            'N', ## N is any 5'-monophosphate nucleotide
            ]
        
        ## HETATM res_names for which coordinates are parsed
        self.d_modres = {
            'MSE':'MET', ## selenomethionine
##            ## phosphorylation
##            'TPO':'THR',
##            'SEP':'SER',
##            'PHD':'ASP',
            }

        self.d_res3 = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR'}

        self.time_status_update = 1

        ## list of biounit strings in which other biounit string can be found
        self.l_biounits = [
            'DIMER OF DIMERS', ## DIMER
            'HEXADECAMER', ## DECAMER
            'DODECAMER', ## DECAMER
            ]
        self.d_biounits = {
            'MONOMER':1,
            'DIMER':2,
            'TRIMER':3,
            'TETRAMER':4,'DIMER OF DIMERS':4, ## DIMER
            'PENTAMER':5,
            'HEXAMER':6,
            'HEPTAMER':7,
            'OCTAMER':8,
            'DECAMER':10,
            'DODECAMER':12 ,'12-MER':12, ## DECAMER
            'HEXADECAMER':16, ## DECAMER
            'ICOSAHEDRAL':60,
            }

        ## http://en.wikipedia.org/wiki/Crystal_structure

        ## Hermann-Mauguin symbols *and* disallowed abbrevations *and* errors (e.g. A 2 space group of 1mbs.pdb)
        ## sorted from low to high symmetry

        ## errornous space groups assigned to a crystal system based on information about the unit cell
        ## from only one representative structure with the space group in question

        ## P 21 21 21, P 1 21 1, C 1 2 1 most common for proteins...

        ## chiral space groups

        d_spacegroups = {
            ## alpha,beta,gamma != 90
            'TRICLINIC':[
                'P 1', ## 1
##                'P 1-','P1', ## neither HM symbols nor abbrevations
##                'A 1', ## 1lks.pdb
                ],
            ## alpha != 90, beta,gamma==90
            'MONOCLINIC':[
                'P 1 2 1', ## 3
                'P 1 21 1', ## 4
                'C 1 2 1', ## 5
##                'C 2', 'C 21', 'C 1 21 1', ## C 1 2 1 abbreviations
##                'P 2', ## P 1 2 1 abbreviations
##                'P 21', ## P 1 21 1 abbreviations
##                'B 2', 'I 1 2 1', 'P 1 1 21', 'I 21', 'I 1 21 1', ## neither HM symbols nor abbrevations
                ],
            ## a != b != c (alpha,beta,gamma==90)
            'ORTHORHOMBIC':[
                'P 2 2 2', ## 16
                'P 2 2 21',
                'P 21 21 2',
                'P 21 21 21',
                'C 2 2 21',
                'C 2 2 2',
                'F 2 2 2',
                'I 2 2 2',
                'I 21 21 21', ## 24
##                'P 2 21 21', ## neither HM symbols nor abbrevations
##                'P 21 2 21', ## not a chiral space group?!
##                'P 21 21 2 A', ## P 21 21 2 error in 1b86.pdb
##                'B 2 21 2', ## 1zna.pdb
##                'B 1 1 2', ## 1qr6.pdb
                ],
            ## a != c (a == b, alpha,beta,gamma==90)
            'TETRAGONAL':[
                'P 4', ## 75
                'P 41', ## 76
                'P 42', ## 77
                'P 43', ## 78
                'I 4', ## 79
                'I 41', ## 80
                'P 4 2 2', ## 89
                'P 4 21 2', ## 90
                'P 41 2 2', ## 91
                'P 41 21 2',
                'P 42 2 2',
                'P 42 21 2',
                'P 43 2 2',
                'P 43 21 2',
                'I 4 2 2',
                'I 41 2 2',
                ],
            ## RHOMBOHEDRAL (a=b=c, alpha,beta,gamma!=90)
            ## alpha,beta,gamma != 90
            'TRIGONAL':[
                'P 3',
                'P 31',
                'P 32',
                'R 3',
                'P 3 1 2', ## 149
                'P 3 2 1',
                'P 31 1 2', ## 151
                'P 31 2 1', ## 152
                'P 32 1 2',
                'P 32 2 1',
                'R 3 2',
##                'H 3', ## R 3 equivalent
##                'H 3 2', ## P 3 2 1 equivalent
                ],
            'HEXAGONAL':[
                'P 61 2 2','P 65','P 63','P 65 2 2','P 61','P 62 2 2','P 62','P 64 2 2','P 63 2 2','P 6 2 2','P 6','P 64',
                ],
            'CUBIC':[
                'F 41 3 2','P 21 3','I 4 3 2','I 2 3','P 2 3','P 41 3 2','P 4 3 2','F 4 3 2','P 43 3 2','I 21 3','F 2 3','P 42 3 2','I 41 3 2',
                ],
##            'UNKNOWN':[
##                'A 2',
##                ]
            }

        self.d_crystalsystems = {}
        for crystalsystem in d_spacegroups.keys():
            for spacegroup in d_spacegroups[crystalsystem]:
                self.d_crystalsystems[spacegroup] = crystalsystem

        self.d_crystalsynonyms = {
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

        self.nontransformationmatrix = [['1.000000', '0.000000', '0.000000', '0.00000'], ['0.000000', '1.000000', '0.000000', '0.00000'], ['0.000000', '0.000000', '1.000000', '0.00000']]

        self.s_alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        ## info from the PDB Ligand Depot (searched for cluster in chemical name)
        self.l_clusters = [
            ## iron clusters
            'CFM','CFN','CLF','CLP','CN1','CNB','CNF','F3S','FES','FS1','FS2','FS4','FSF','FSO','FSX','HC0','HC1','HF3','HF5','NFS','SF3','SF4','WCC','XCC', ## 'FS3' deprecated ('F3S' maintained)
            ## copper clusters
            'CUB','CUM','CUN','CUO',
            ## molybdenum "clusters"
            'OMO',
            ## hafnium clusters
            'PHF',
            ## zirconium clusters
            'ZRC',
            ]

        self.l_prosthetic_groups = [
            ## porphyrins (cyclic tetrapyrroles)
            ## Ferrochelatase catalyzes protophorphyrin+Fe(2+) --> protoheme + 2H(+)
            ## iron
            'HEM', ## protoporphyrin IX + Fe(II) (C3 vinyl,C8 vinyl,C18 methyl; tetramethyl,divinyl,dipropionate; *charge 2*)
            'HEC', ## Heme C (protoporphyrin IX; *charge 0*)
            'HEA', ## Heme A (C3 hydroxyfarnesyl, C8 vinyl, C18 formyl)
            'HEO', ## Heme O (C3 hydroxyfarnesyl, C8 vinyl, C18 methyl)
            '2FH', ## 2-phenylheme
            '1FH', ## 12-phenylheme
            'DDH', ## dedivinyl,diacetyl heme (C3,C8)
            'HEV', ## dedimethyl,divinyl heme
            'HDM', ## tetramethyl,divinyl,dipropionate *ester* heme
            'HAS', ## Heme-As (C3 geranylgeranyl, C8 vinyl, C18 formyl)
            'VER', ## octaethylated porphyrin
            'HEB', ## Heme B/C hybrid (tetramethyl,divinyl,dipropionate)
            'DHE', ## Heme D (heme B derivative)
            'HDD', ## cis-heme D hydroxychlorin gamma-spirolactone
            'SRM', ## siroheme (partially reduced iron-porphyrin in e.g. nitrate reductase)
            ## other metal
            'HNI', ## protoporphyrin IX + Ni(II)
            'HES', ## Zn substituted Heme C
            ]
        self.l_coenzymes = [
            'RET', ## vitamin A
##            'TPP', ## vitamin B1            
            'FMN','FAD', ## vitamin B2
            'NAD','NAP', ## vitamin B3
            'COA', ## vitamin B5
            'PLP', ## vitamin B6
            'C2F', ## vitamin B9 (5-methyl THF)
            ]

        ## list of metals also contains metalloids (e.g. Arsen)
        self.l_atoms_metal = [
            'LI','NA','K','CS', ##1a
            'BE','MG','CA', ##2a
            'AL','GA','TL', ##3a
            'PB', ##4a
            'AS', ##5a
            'V', ##3b
            'CR','MO', ##4b
            'MN', ##5b
            'FE', ##6b
            'CO', ##7b
            'NI', ##8b
            'CU', ##9b
            'ZN','CD','HG', ##10b
            ]

        ## info from the PDB Ligand Depot
        ## keys are hetIDs, values are charges (not oxidation states)
        ## hetID:[chemical formula,charge]
        self.d_ions = {

            ## unknown
            'UNX':['',''], ## e.g. 1aqn

            ## nitrate, ammonium
            'NO3':['N1 O3',-1],'NH4':['H4 N1',+1],
            ## hydroxide
            'OH' :['H1 O1',-1],
            ## phosphate
            '2HP':['O4 P1',-1],'PI' :['H1 O4 P1',-2],'PO4':['O4 P1',-3], ## different oxidation states; 'IPS' deprecated
            ## sulfate
            'SO4':['O4 S1',-2],'SOH':[3,-1],'SUL':[3,-2], ## different oxidation states
            ## carbonate
            'CO3':['C O3', -1],
            ## group1a
            'LI' :['LI1',+1],
            'NA' :['NA1',+1], ## 'NAO','NA2','NA6','NA5','NAW' deprecated
            'K'  :['K1' ,+1], ## 'KO4' deprecated
            'CS' :['CS' ,+1],
            ## group2a
            'BEF':['BE F3',-1],
            'MG' :['MG1',+2], ## 'MO3','MO1','MO2','MO4','MO5','MO6' deprecated
            'CA' :['CA1',+2],'OC1':['CA1',+2], ## 'OC3','OC5','OC6','OC7','OC8','OC2','OC4' deprecated
            'SR' :['SR' ,+2],
            'BA' :['BA' ,+2],
            ## group3a
            'AL' :['AL1',+3],'ALF' :['AL F4',-1],
            'GA' :['GA1',+3],'GA1' :['GA1',+2],
            'TL' :['TL1',+1],
            ## group4a
            'ARS':['AS1', 0],'ART':['O4 AS1',-3],'AST':-3,'TAS':['H3 O3 AS1', 0],'CAC':['C2 H6 AS O2', -1], ## different compounds
            'PB' :['PB' ,+2],
            ## group6a
            'SE' :['SE1', 0],'SE4':['O4 SE1',-2], ## different compounds
            ## group7a
            'CL' :['CL1',-1],
            'BR' :['BR1',-1],
            'IOD':['I1' ,-1],
            ## group8a
            'KR' :['KR1', 0],
            ## group3b
            'V'  :+3,'VO4':['V1' ,-3], ## different oxidation states
            ## group4b
            'CR' :['CR1',+3],
            'MO' :['MO1', 0],'4MO':['MO1', 0],'2MO':['MO O2',-2], ## different compounds and different oxidation states
            ## group5b
            'MN' :['MN1',+2],'MN3':['MN1',+3], ## different oxidation states; 'MW1','MW2','MW3','O4M','MN5','MN6' deprecated
            ## group6b
            'FE2':['FE1',+2],'FE' :['FE1',+3], ## different oxidation states; 'OF1','OF3','2OF' deprecated
            ## group7b
            'CO' :['CO1',+2],'3CO':['CO1',+3], ## different oxidation states; 'CO5','OCL','OCO','OCN','OCM' deprecated
            ## group8b
            'NI' :['NI1',+2],'3NI':['NI1',+3], ## different oxidation states; 'NI1','NI2','NI3','NIK' deprecated
            ## group9b
            'CU1':['CU1',+1],'CU' :['CU1',+2], ## different oxidation states; '1CU' deprecated
            ## group10b
            'ZN' :['ZN1',+2], ## 'ZN2','ZO3','ZN3','ZNO' deprecated
            'CD' :['CD1',+2],
            'HG' :['HG1',+2],
            ## Lanthanides
            'TB' :['TB1',+3],
            'YB' :['YB1',+3],
            }

        self.l_cofactors = self.l_clusters+self.l_prosthetic_groups+self.l_coenzymes+self.d_ions.keys()

        ## wikipedia buffer solution, Good's buffers
        self.l_solutes = [

            ## unknown atom or ion
            'UNX',

            ## IPA,FMT,GOL,EEE,EDO
            ##
            ## ethylene glycol, protein precipitation
            'EDO',
            ## acetic acid
            'ACY',
            ## water
            'HOH','DOD',
            ## methanol
            'MOH',
            ## di-thio-threitol (reducing agent)
            'DTT', 
            ## beta-mercapto-ethanol (reducing agent)
            'BME',
            ## bis-tris methane (buffering agent)
            'BTB',

            ## TAPS (buffering agent)
            'T3A',
            ## Bicine (buffering agent)
            'BCN',
            ## Tris (buffering agent)
            'TRS',
            ## HEPES (buffering agent)
            'EPE',
            ## TES (buffering agent)
            'NES',
            ## MOPS
            'MPO',
            ## PIPES
            'PIN',
            ## Cacodylate
            'CAC',
            ## MES
            'MES',
            ## Acetate
            'ACT',

            ## Glycerol
            'GOL',
            ]

        ## saccharides returned from a search of the ligand depot for disaccharides and monosaccharides
        self.d_saccharides = {
            ##
            ## monosaccharides, aldehydes
            ##
            ## pyranoses, hexoses
            'GLC':{'stereo':'GLC','derivate':['GLC']}, ## (alpha)-D-Glucose
##            'AGC':{'stereo':'GLC','derivate':['GLC']}, ## alpha-D-Glc (deprecated)
            'BGC':{'stereo':'GLC','derivate':['GLC']}, ## beta-D-Glc
            'GAL':{'stereo':'GAL','derivate':['GAL']}, ## (beta)-D-Galactose
            'GLA':{'stereo':'GAL','derivate':['GAL']}, ## alpha-D-Gal
##            'GLB':{'stereo':'GAL','derivate':['GAL']}, ## beta-D-Gal
            'FUC':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, alpha-L-Fucose
            'FUL':{'stereo':'FUC','derivate':['GAL']}, ## 6-deoxy-GAL, beta-L-Fucose
            'MAN':{'stereo':'MAN','derivate':['MAN']}, ## alpha-D-Mannose
            'BMA':{'stereo':'MAN','derivate':['MAN']}, ## beta-D-Mannose
            'ARA':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabinose
            'ARB':{'stereo':'ARA','derivate':['ARA']}, ## beta-L-Arabinose
            ## pyranoses, pentoses
            'XYS':{'stereo':'XYS','derivate':['XYS']}, ## (alpha)-D-Xylose
            'XYP':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-Xylose
            'LXC':{'stereo':'XYS','derivate':['XYS']}, ## beta-L-Xylose
            ## furanoses, hexoses
            'AHR':{'stereo':'ARA','derivate':['ARA']}, ## alpha-L-Arabino*furano*se
            ## furanoses, pentoses
            'XYZ':{'stereo':'XYS','derivate':['XYS']}, ## beta-D-xylo*furano*se
            ## phosphorylated aldohexopyranoses
            'G1P':{'stereo':'G1P','derivate':['GLC']}, ## alpha-D-Glc-1P
            'G6P':{'stereo':'G6P','derivate':['GLC']}, ## alpha-D-Glc-6P
            'BG6':{'stereo':'G6P','derivate':['GLC']}, ## beta-D-Glc-6P
            'G6Q':{'stereo':'G6P','derivate':['GLC']}, ## Glc-6P
            'BGP':{'stereo':'BGP','derivate':['GAL']}, ## beta-Gal-6P
            'M1P':{'stereo':'M1P','derivate':['MAN']}, ## alpha-D-Man-1P
            'M6P':{'stereo':'M6P','derivate':['MAN']}, ## alpha-D-Man-6P
            ## phosphorylated aldopentofuranoses
            'ABF':{'stereo':'ABF','derivate':['ARA']}, ## beta-D-Arabino*furano*se-5-phosphate
            ## methylated aldohexopyranoses
            'MMA':{'stereo':'MMA','derivate':['MAN']}, ## O1-methyl-mannose
            ## aminated aldohexopyranoses amine
            'AGL':{'stereo':'AGL','derivate':['GLC']}, ## 4,6-dideoxy-4-amino-alpha-D-Glucose
            ## deoxygenated aldohexopyranoses
            'G6D':{'stereo':'G6D','derivate':['GLC']}, ## 6-deoxy-alpha-D-Glucose
            ## oxygenated aldohexopyranoses
            'KBG':{'stereo':'G6D','derivate':['GLC']}, ## 2-keto-beta-D-Glucose
            ## acetylated aldohexopyranose amines
            'NAG':{'stereo':'NAG','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine
            'NBG':{'stereo':'NAG','derivate':['GLC']}, ## 1-N-Acetyl-beta-D-Glucosamine
            'NDG':{'stereo':'NAG','derivate':['GLC']}, ## 2-(acetylamino)-2-deoxy-alpha-D-glucopyranose
            '16G':{'stereo':'16G','derivate':['GLC']}, ## (2)-N-Acetyl-D-Glucosamine-6P
            'NGA':{'stereo':'NGA','derivate':['GLC']}, ## (2)-N-Acetyl-D-Galactosamine
            ##
            ## monosaccharides, ketones
            ##
            ## furanoses, hexoses
            'FRU':{'stereo':'FRU','derivate':['FRU']}, ## Fructose
            'F6P':{'stereo':'F6P','derivate':['FRU']}, ## Fru-6P
            ##
            ## dissacharides
            ##
            'SUC':{'stereo':'SUC','derivate':['GLC','FRU']}, ## GLC-a12-FRC, Sucrose
            'LAT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, alpha-Lactose
            'LBT':{'stereo':'LAT','derivate':['GAL','GLC']}, ## GAL-b14-GLC, beta-Lactose
            'MAL':{'stereo':'MAL','derivate':['GLC','GLC']}, ## GLC-a14-GLC, Maltose
            'TRE':{'stereo':'TRE','derivate':['GLC','GLC']}, ## GLC-a11a-GLC, Trehalose
            'CBI':{'stereo':'CBI','derivate':['GLC','GLC']}, ## GLC-b14-GLC, Cellobiose
            ##
            ## polysaccharides
            ##
            'MTT':{'stereo':'MTT','derivate':['MAL','MAL']}, ## MAL-b14-MAL, Maltotetraose
            ##
            ## linear saccharides (neither furanoses nor pyranoses...) and their derivatives...
            ##
            'SOR':{'stereo':'SOR','derivate':['GLC']}, ## Sorbitol/Glucitol (reduced glucose)
            'GLO':{'stereo':'GLO','derivate':['GLC']}, ## linear glucose
            'XLS':{'stereo':'XLS','derivate':['XYL']}, ## linear xylose
            'A5P':{'stereo':'ABF','derivate':['ARA']}, ## Arabinose-5-phosphate
            ##
            ## conduritol (1,2,3,4-cyclohexenetetrol) derivatives
            ##
            'HMC':{'stereo':'HMC'}, ## 5-hydroxymethyl-chonduritol
            'ACI':{'stereo':'ACI'}, ## 1-amino-5-hydroxymethyl-chonduritol
            ## Sialic acid (N-Acetylneuraminic acid, Neu5Ac, NANA)
            'SIA':{'stereo':'SIA','derivate':['SIA']}, ## (alpha)-sialic acid
            'SLB':{'stereo':'SIA','derivate':['SIA']}, ## beta-sialic acid
            }

        self.d_stereoisomers = {
            }

        self.l_expdta = [
            'X-RAY',
            'ELECTRON DIFFRACTION','NEUTRON DIFFRACTION',
            'NMR',
            'INFRARED SPECTROSCOPY',
            'CRYO-ELECTRON MICROSCOPY',
            'ELECTRON TOMOGRAPHY', ## e.g. 1o1a.pdb
            'SOLUTION SCATTERING', ## e.g. 1e07.pdb
            'FLUORESCENCE TRANSFER', ## e.g. 1rmn.pdb
            ]

        self.l_columns_html = [
            'gif1','gif2','pdb1', 'pdb2', 'bm1', 'bm2',
            'rmsd', 'mutations', 'chains', 'residues', 'coordinates',
            'pH1', 'pH2', 'T1', 'T2', 'res1', 'res2','spacegroup1', 'spacegroup2',
            'REMARK465', 'REMARK470','transformations',
            'title1','title2','hetIDs1', 'hetIDs2'
            ]

        self.d_symop = {
            '1mzy':[1,4,7,],
            '1n70':[9,6,3,],
            }

        self.maxrmsd = 2.75 ## 2fsy.pdb,2ft1.pdb ## try other combinations if above this rmsd
        self.maxrmsd_wrong = 9.5 ## 2eu1 ## assume error if above this rmsd

        self.minres = 5.0

        self.min_len_chain = 50 ## at least one chain in the pdb must be 50 residues or longer if pdb is to be processed
        self.max_len_chain_difference = 25
        self.max_mutations = 10

        self.path_pdb = '/oxygenase_local/data/pdb/'
        self.path_cwd = os.getcwd()
        
        return

if __name__ == '__main__':
    instance_quakes = quakes()
    instance_quakes.main()

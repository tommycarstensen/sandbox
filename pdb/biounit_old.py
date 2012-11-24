#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2007

import os, sys, Numeric, LinearAlgebra, math, urllib2

class biounit:


    def main(self, pdb):

        d_transformations, d_chain_ids = self.parse_pisa_multimers(pdb)
        print d_transformations
        print pdb

        if d_transformations == {}:

            stop

            ## skip if asymmetric unit == biological unit (case 1 of 2)
            fd = open('%s_1.pdb' %(pdb),'w')
            fd.close()
##            os.system('cp %s%s/pdb%s.ent %s_1.pdb' %(self.path_pdb,pdb[1:3],pdb,pdb))

        else:

            d_output = self.parse_pdb_coordinates(pdb, d_transformations, d_chain_ids)

            for assembly in d_transformations.keys():
                self.write_transformed_coordinates(pdb,assembly,d_output)
            
        return


    def write_transformed_coordinates(self,pdb,assembly,d_lines):

        l_lines = []
        for chain in d_lines.keys():
            l_lines += d_lines[chain]
        fd = open('%s_%s.pdb' %(pdb,assembly), 'w')
        fd.writelines(l_lines)
        fd.close()

        return


    def parse_pdb_coordinates(self, pdb, d_transformations, d_chain_ids):

        fd = open('%s%s/pdb%s.ent' %(self.path_pdb,pdb[1:3],pdb), 'r')
        lines = fd.readlines()
        fd.close()

        set_chains = set([
            'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
            'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
            '0','1','2','3','4','5','6','7','8','9',
            ])
        d_chains = {}
        d_output = {}

        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()

            if record == 'SEQRES':
                chain = line[11]
                set_chains -= set([chain])

            elif record == 'HET':
                chain = line[12]
                set_chains -= set([chain])

            if record == 'ATOM':
                d_chains, set_chains, d_output = self.parse_recordATOM(lines, i, d_transformations, d_chain_ids, set_chains, d_chains, d_output)

            elif record == 'HETATM':
                d_chains, set_chains, d_output = self.parse_recordATOM(lines, i, d_transformations, d_chain_ids, set_chains, d_chains, d_output)

        return d_output


    def parse_pisa_multimers(self, pdb):

        url = urllib2.urlopen('http://www.ebi.ac.uk/msd-srv/prot_int/cgi-bin/multimers.pisa?%s' %(pdb))
        lines = url.readlines()
        assembly = 0
        d_transformations = {}
        d_chain_ids = {}
        for i in range(len(lines)):

            line = lines[i]
##            if '<r' not in line and '<tx' not in line and '<ty' not in line and '<tz' not in line:
##                print line[:-1]
##            print line[:-1]

            lvl = 2
            s_index1 = '<status>'
            s_index2 = '</status>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
                index2 = line.index(s_index2)
                index1 = line[:index2].rindex(s_index1)+len(s_index1)
                status = line[index1:index2]
                if status not in [
                    'Ok','No symmetry operations','Entry not found',
                    'Overlapping structures',
                    ]:
                    print 'status:', status
                    stopstatus

            lvl = 3
            s_index1 = '<ser_no>'
            s_index2 = '</ser_no>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
                index2 = line.index(s_index2)
                index1 = line[:index2].rindex(s_index1)+len(s_index1)
                ser_no = int(line[index1:index2])
                if ser_no != 1:
                    break
                if ser_no > 1:
                    print pdb
                    stop

            lvl = 3
            s_index1 = '<assembly>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1:
                assembly += 1
                if assembly in d_transformations.keys():
                    stop
                d_transformations[assembly] = {'chains':{},}
                molecule = 0

            lvl = 4
            s_index1 = '<id>'
            s_index2 = '</id>'
            if line[2*lvl:2*lvl+len(s_index1)] == s_index1 and line[-1-len(s_index2):-1] == s_index2:
                index2 = line.index(s_index2)
                index1 = line[:index2].rindex(s_index1)+len(s_index1)
                id = int(line[index1:index2])

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
                if chain_id in d_chain_ids.keys():
                    if d_chain_ids[chain_id] != assembly:
                        print pdb
                        print chain_id, assembly
                        print d_chain_ids[chain_id]
                        notexpectedmaybevirus
                d_chain_ids[chain_id] = assembly
                if chain_id not in d_transformations[assembly]['chains'].keys():
                    d_transformations[assembly]['chains'][chain_id] = {}
                if molecule in d_transformations[assembly]['chains'][chain_id].keys():
                    notexpected
                d_transformations[assembly]['chains'][chain_id][molecule] = {
                    'r':Numeric.zeros((3,3), typecode='d'),
                    't':Numeric.array([0.,0.,0.,]),
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

        return d_transformations, d_chain_ids


    def parse_recordATOM(self, lines, i, d_transformations, d_chain_ids, set_chains, d_chains, d_output):

        import Numeric

        line = lines[i]

        record = line[:6].rstrip()
        atom_no = int(line[6:11])
        atom_name = line[12:16]
        altloc = line[16]
        res_name = line[17:20].strip()
        chain = line[21]
        res_no = int(line[22:26])
        iCode = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coordinate = Numeric.array([x, y, z])
        occupancy = float(line[54:60])
        bfactor = float(line[60:66])
        element = line[76:78].strip()
        charge = line[78:80].strip()

        pisa_chain_id = '[%s]%s:%i%s' %(res_name,chain.replace(' ','-'),res_no,iCode.replace(' ',''))
        if pisa_chain_id in d_chain_ids and record == 'HETATM':
            assembly = d_chain_ids[pisa_chain_id]
        elif chain in d_chain_ids.keys():# and record == 'ATOM':
            if chain == ' ':
                stop
            pisa_chain_id = chain
            assembly = d_chain_ids[pisa_chain_id]
        elif res_name in ['HOH','DOD',]:
            return d_chains, set_chains,d_output
        elif chain == ' ':
            print line
            print d_chain_ids.keys()
            if res_name == 'MSE':
                stop
            pisa_chain_id = '-'
            if pisa_chain_id not in d_chain_ids.keys():
                if res_name in [
                    'UNK','UNX', ## unknown residues and ions
                    'ACY','EDO', ## solutes
##                    'MG', ## ions
                    ]:
                    return d_chains, set_chains, d_output
                elif res_name in [
                    'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET',
                    'ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',
                    ] and record == 'HETATM':
                    stop_temporary
                    return d_chains, set_chains, d_output
            else:
                assembly = d_chain_ids[pisa_chain_id]
        elif res_name in ['UNK','UNX',]:
            stop
            return d_chains, set_chains, d_output
        else:
            print chain, d_chain_ids.keys()
            print line
            notexpected

        for molecule in d_transformations[assembly]['chains'][pisa_chain_id].keys():
            matrix = d_transformations[assembly]['chains'][pisa_chain_id][molecule]['r']
            vector = d_transformations[assembly]['chains'][pisa_chain_id][molecule]['t']
            ## skip if asymmetric unit == biological unit (case 2 of 2)
            if (
                matrix[0][0] == 1 and matrix[0][1] == 0 and matrix[0][2] == 0 and
                matrix[1][0] == 0 and matrix[1][1] == 1 and matrix[1][2] == 0 and
                matrix[2][0] == 0 and matrix[2][1] == 0 and matrix[2][2] == 1 and
                vector[0] == 0 and vector[1] == 0 and vector[2] == 0
                ):
                continue
            coordinate_transformed = Numeric.matrixmultiply(matrix,coordinate)+vector
            x = coordinate_transformed[0]
            y = coordinate_transformed[1]
            z = coordinate_transformed[2]

            if chain+'_'+str(molecule) in d_chains.keys():
                chain2 = d_chains[chain+'_'+str(molecule)]
            else:
                chain2 = list(set_chains)[0]
                set_chains -= set([chain2])
                d_chains[chain+'_'+str(molecule)] = chain2

            line = '%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' %(
                record.ljust(6), atom_no, atom_name.ljust(4),
                altloc, res_name.rjust(3), chain2, res_no, iCode,
                x, y, z, occupancy, bfactor, element.rjust(2), charge.rjust(2),
                )
            if lines[i][:12] != line[:12] or lines[i][16:21] != line[16:21] or lines[i][22:30] != line[22:30] or lines[i][54:] != line[54:]:
                print atom_name
                print lines[i]
                print line
                print lines[i][:30]
                print line[:30]
                print lines[i][54:]
                print line[54:]
                stop

            print assembly, chain, set_chains, chain2
            if chain2 not in d_output[assembly].keys():
                d_output[assembly][chain2] = []
            d_output[assembly][chain2] += [line]

        return d_chains, set_chains, d_output


    def parse_pdbheader(self, lines, s_pdb):

        s_pdb = s_pdb.lower()
        ## import sequence alignment
        import sys
        sys.path.append('/home/people/tc/svn/EAT_DB/trunk/')
        import sequence_alignment

        ## parser written on the assumption that SEQRES is mandatory if ATOM records exist

        d_seq = {}
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

            if record == 'REMARK': ## section 2
                d_seq = self.parse_recordREMARK(d_seq, line, i, lines)

            elif record == 'SEQRES': ## section 3
                d_seq = self.parse_recordSEQRES(line, d_seq)

            elif record == 'HET': ## section 4
                hetID = line[7:10].strip()
                ## continue if water
                if hetID in ['HOH','H2O','DOD','D2O']: ## D2O in 2JAJ
                    continue
                chain = line[12]
                res_no = int(line[13:17])
                iCode = line[17]
                if 'HET' not in d_seq.keys():
                    d_seq['HET'] = {}
                if chain not in d_seq['HET'].keys():
                    d_seq['HET'][chain] = {}
                if res_no not in d_seq['HET'][chain].keys():
                    d_seq['HET'][chain][res_no] = {}
                if iCode not in d_seq['HET'][chain][res_no].keys():
                    d_seq['HET'][chain][res_no][iCode] = hetID
                else:
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
                if txt not in [] and hetID in set(self.d_res.keys()+self.l_nucleotides)-set(['MSE']):
                    print line
                    stop
                if 'MODRES' not in d_seq.keys():
                    d_seq['MODRES'] = {}
                if chain not in d_seq['MODRES'].keys():
                    d_seq['MODRES'][chain] = {}
                if res_no not in d_seq['MODRES'][chain].keys():
                    d_seq['MODRES'][chain][res_no] = {}
                if iCode not in d_seq['MODRES'][chain][res_no].keys():
                    d_seq['MODRES'][chain][res_no][iCode] = hetID
                elif hetID not in ['CSY']: ## SER-TYR-GLY chromophore
                    print line, s_pdb
                    notexpected

            elif record == 'TITLE': ## section 2
                if not 'TITLE' in d_seq.keys():
                    d_seq['TITLE'] = line[10:].strip()
                else:
                    if d_seq['TITLE'][-1] == '-':
                        d_seq['TITLE'] += line[10:].strip()
                    else:
                        d_seq['TITLE'] += ' '+line[10:].strip()

            elif record == 'HEADER':
                d_seq['HEADER'] = line[10:50].strip()

            elif record == 'SPRSDE': ## section 2
                sIDcodes = line[31:70].lower().split()
                for sIDcode in sIDcodes:
                    if sIDcode not in d_seq.keys():
                        d_seq[sIDcode] = {}
                    d_seq[sIDcode]['SPRSDE'] = True

            elif record == 'EXPDTA': ## section 2
                methods = line[10:].strip().split(',')[0]
                if methods[:3] == 'NMR':
                    methods = 'NMR'
                elif 'X-RAY' in methods or methods == 'FIBER DIFFRACTION' or methods == 'POWDER DIFFRACTION': ## e.g. (synchrotron) x-ray (fiber, powder) diffraction
                    methods = 'X-RAY'
                d_seq['EXPDTA'] = methods

            elif record == 'CRYST1': ## section 8
                spacegroup = line[55:66].strip()
                d_seq['CRYST1'] = spacegroup


        return d_seq


    def parse_recordREMARK(self, d_seq, line, i, lines):

        remark = int(line[6:10])

        if remark == 525: ## water association

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
                        if 'REMARK525' not in d_seq.keys():
                            d_seq['REMARK525'] = []
                        d_seq['REMARK525'] += [solventchain]

        return d_seq


    def __init__(self):

##        self.path_pdb = '/oxygenase_local/data/pdb/'
        self.path_pdb = '/local/data/pdb/'

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

    instance_biounit = biounit()
    instance_biounit.main(pdb)

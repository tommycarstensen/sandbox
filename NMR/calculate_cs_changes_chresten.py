#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

class delta_cs:

    def main(
        self,l_wts,d_pred,l_xtics,
        ):

        import os, sys
        sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
        import gnuplot, statistics

        ## parse experimental data
        d_exp = self.dic2csv(l_xtics)

        ## get cwd
        dir_main = os.getcwd()

        l_r = []

        for pdb in l_wts:

            print pdb, l_wts.index(pdb)

            if not os.path.isdir('%s/%s' %(dir_main,pdb)):
                os.mkdir('%s/%s' %(dir_main,pdb))

            os.chdir('%s/%s' %(dir_main,pdb))

            self.pre_whatif(pdb)

            if pdb in ['2vb1','1vdp',]:
                os.system('cp %s_monomer.pdb %s_protonated.pdb' %(pdb,pdb))
##            else:
##                self.whatif(pdb)

##            self.calculate_chemical_shifts(pdb)

            ## parse computational predictions
            d_pred = self.parse_chemical_shifts(pdb,d_pred)

            ## calculate correlation coefficients
            l_exp = []
            l_pred = []
            for titgrp in d_exp.keys():
                res_number = int(titgrp[1:])
                res_symbol = titgrp[0]
                res_name = self.d_ressymbol2resname[res_symbol]
                for nucleus in d_exp[titgrp].keys():
                    cs_exp = d_exp[titgrp][nucleus]
                    l_exp += [cs_exp]
                    index = nucleus.index('N-HN')
                    cs_pred = d_pred['%s%i' %(res_name,res_number)][nucleus[:index]][-1]
                    l_pred += [cs_pred]
                r = statistics.correlation(l_exp,l_pred)
                l_r += [r]
##                print titgrp,r

##            print sum(l_r)/len(l_r), min(l_r), max(l_r)

        ## change from local dir to main dir
        os.chdir(dir_main)

        ## plots
        for titgrp1 in d_exp.keys()+['E35']:
            res_number = int(titgrp1[1:])
            res_symbol = titgrp1[0]
            res_name = self.d_ressymbol2resname[res_symbol]
            titgrp3 = '%s%i' %(res_name,res_number)
            prefix = 'delta_cs_%s' %(titgrp3)
            ylabel = '{/Symbol D}{/Symbol w}_H'
            title = titgrp3
            gnuplot.histogram(
                d_pred[titgrp3],prefix,l_xtics,
                ylabel=ylabel,title=title,
##                l_plotdatafiles=['E34.txt'],
                )

        return


    def dic2csv(self,l_xtics):

        d_in = {}
        d_out = {}

##        ##
##        ## convert experimental data
##        ##
##        fd = open('0034GLU.csv','r')
##        lines1 = fd.readlines()
##        fd.close()
##
##        lines2 = []
##        for line1 in lines1:
##            nucleus = line1.split(',')[0]
##            if nucleus[-4:] != 'N-HN':
##                continue
##            residue = nucleus[:-4]
##            i = l_xtics.index(residue)
##            cs_helen = float(line1.split(',')[1])
##            lines2 += ['%f %f\n' %(i, cs_helen)]
##
##            d_out[residue] = cs_helen
##
##        fd = open('E34.txt','w')
##        fd.writelines(lines2)
##        fd.close()

        l_files = [
            'wt__col__pKas',
            'E7Q__col__pKas',
            'D18N__col__pKas',
            'D48N__col__pKas',
            'D52N__col__pKas',
            'D66N__col__pKas',
            'D87N__col__pKas',
            'D101N__col__pKas',
            'D119N__col__pKas',
            ]

        for file in l_files:
            index = file.index('__col__pKas')
            if file == 'wt__col__pKas':
                residue = 'wt'
            else:
                residue = file[:index-1]
            fd = open(file,'r')
            s = fd.read()
            fd.close()
            d_in[residue] = eval(s)
            del d_in[residue]['__fit_matches__']
            del d_in[residue]['__datatabs_fits__']
            del d_in[residue]['__Exp_Meta_Dat__']
            if residue == 'wt':
                del d_in[residue]['__datatab_structmapping__']


        ## parse wt chemical shifts
        d_pHs_wt = {}
        for nucleus in d_in['wt'].keys():
            d_pHs_wt[nucleus] = {}
            for key in d_in['wt'][nucleus][0].keys():
                if key == 'label':
                    continue
                pH = float(d_in['wt'][nucleus][0][key]['var'])
                cs = float(d_in['wt'][nucleus][1][key]['var'])
                d_pHs_wt[nucleus][pH] = cs
                

        ## parse mutant chemical shifts
        for tit_grp in d_in.keys():
            if tit_grp == 'wt':
                continue
            d_out[tit_grp] = {}
            for nucleus in d_in[tit_grp].keys():
                if nucleus in ['data','?-?']:
                    continue
                if nucleus[-4:] != 'N-HN':
                    continue
##                d_out[tit_grp][nucleus] = {}
                d_pHs = {}
                for key in d_in[tit_grp][nucleus][0].keys():
                    if key == 'label':
                        continue
                    pH = float(d_in[tit_grp][nucleus][0][key]['var'])
                    d_pHs[pH] = key

                    cs = float(d_in[tit_grp][nucleus][1][key]['var'])
##                    d_out[tit_grp][nucleus][pH] = cs

                pH_min = min(d_pHs.keys())
                pH_max = max(d_pHs.keys())
                key_min = d_pHs[pH_min]
                key_max = d_pHs[pH_max]
                cs_min = float(d_in[tit_grp][nucleus][1][key_min]['var'])
                cs_max = float(d_in[tit_grp][nucleus][1][key_max]['var'])

                if not nucleus in d_pHs_wt.keys():
                    continue

                l_pHs_wt = d_pHs_wt[nucleus].keys()
                l_pHs_wt.sort()
                for i in range(len(l_pHs_wt)):
                    pH = l_pHs_wt[i]
                    if pH > pH_max:
                        break

                pH_wt1 = l_pHs_wt[i-1]
                pH_wt2 = l_pHs_wt[i]
                cs_wt1 = d_pHs_wt[nucleus][pH_wt1]
                cs_wt2 = d_pHs_wt[nucleus][pH_wt2]

##                print pH_max, pH_wt1, pH_wt2, cs_wt1, cs_wt2
                percentage1 = (pH_max-pH_wt1)/(pH_wt2-pH_wt1)
##                percentage2 = (pH_wt2-pH_max)/(pH_wt2-pH_wt1)
                cs_extrapolation1 = cs_wt1+percentage1*(cs_wt2-cs_wt1)
##                cs_extrapolation2 = cs_wt2-percentage2*(cs_wt2-cs_wt1)

                d_out[tit_grp][nucleus] = cs_extrapolation1-cs_max


        ##
        ## write experimental data to file
        ##
        lines = []

        for i in range(len(l_xtics)):
            residue = l_xtics[i]
            nucleus = residue+'N-HN'
            if nucleus not in d_in.keys():
                continue
            
            d_pHs = {}
            for key in d_in[nucleus][0].keys():
                if key == 'label':
                    continue
                pH = float(d_in[nucleus][0][key]['var'])
                d_pHs[pH] = key

            pH_min = min(d_pHs.keys())
            pH_max = max(d_pHs.keys())
            key_min = d_pHs[pH_min]
            key_max = d_pHs[pH_max]
            cs_min = float(d_in[nucleus][1][key_min]['var'])
            cs_max = float(d_in[nucleus][1][key_max]['var'])

            cs_diff = cs_max-cs_min

            lines += ['%f %f\n' %(i, cs_diff)]

        fd = open('E34A.csv','w')
        fd.writelines(lines)
        fd.close()

        return d_out


    def parse_chemical_shifts(self,pdb,d_data):

        for titgrp in d_data.keys():
            fd = open('%s_%s_%i.pdb.m' %(pdb,titgrp[:3],int(titgrp[3:])),'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                residue = line.split()[0]
                atom = line.split()[1]
                cs = line.split()[2]
                if atom != 'H' and cs != '0.00':
                    print line
                    print titgrp
                    stop
                if atom == 'H':
                    d_data[titgrp][residue] += [float(cs)]

        return d_data


    def calculate_chemical_shifts(self,pdb):

        import os

        ## copy sybyl
        if not os.path.isdir('parameters'):
            os.mkdir('parameters')
        os.system('cp ../parameters/sybyl_types.txt parameters/.')

        fd = open('cs.src','w')
        fd.writelines([
            'load %s_protonated.pdb %s\n' %(pdb,pdb),
            'task chemical_shift_changes {%s}\n' %(pdb),
            ])
        fd.close()
        os.system('../genialtNavn cs.src')

        return


    def pre_whatif(self,pdb):

        import os

        ## copy pdb
        os.system('cp /oxygenase_local/data/pdb/%s/pdb%s.ent %s.pdb' %(pdb[1:3],pdb,pdb))
        ## remove chains
        fd = open('%s.pdb' %(pdb),'r')
        lines1 = fd.readlines()
        fd.close()
        lines2 = []
        for line in lines1:
            record = line[:6].strip()
            if record == 'ATOM':
                chain = line[21]
                if chain == 'A':
                    lines2 += [line]
        fd = open('%s_monomer.pdb' %(pdb),'w')
        fd.writelines(lines2)
        fd.close()
        if len(lines2) < 1000:
            print len(lines2)
            stop

        return


    def whatif(self,pdb):

        import os

        ## remove protonated pdb
        if os.path.isfile('%s_protonated.pdb' %(pdb)):
            os.remove('%s_protonated.pdb' %(pdb))
        ## copy topology
        os.system('cp /software/whatif_debugged/dbdata/TOPOLOGY.H .')

        ##
        ## add hydrogens
        ##
        fd = open('whatif.src','w')
        fd.writelines([
            ## start WHATIF
            '/software/whatif_debugged/DO_WHATIF.COM <<EOF\n',
            ## get molecule
            'GETMOL %s_monomer.pdb\n' %(pdb),
                ## set name
                '%s\n' %(pdb),
            ## delete water
            '%DELWAT\n',
            ## set method of hydrogen addition
            'SETWIF 339 1\n',
            ## add hydrogen atoms
            '%ADDHYD\n',
                ## range
                'TOT 0\n',
            ## make molecule
            'MAKMOL\n',
                ## header
                '%s.pdb\n' %(pdb),
                ## output
                '%s_protonated.pdb\n' %(pdb),
                ## range
                'TOT 0\n',
                ## remark
                'hydrogen atoms added\n',
                ## eol
                '\n',
            ## stop WHATIF    
            'STOP\n',
            'Y\n',
            ])
        fd.close()
        os.system('source whatif.src > tmp.txt')

        return


    def __init__(self):

        ## with and without nitrate
        l_wts = [
            ## without hydrogen
            '2lzt', '1aki', '1bhz', '1bvx', '1bwh', '1bwi', '1bwj', '1f0w', '1f10', '1hel', '1hsw', '1hsx', '1ja2', '1ja4', '1ja6', '1jis', '1jit', '1jiy', '1jj1', '1jj3', '1jpo', '1lj3', '1lj4', '1lje', '1ljf', '1ljg', '1ljh', '1lji', '1ljj', '1ljk', '1lks', '1lma', '1lsa', '1lsb', '1lsc', '1lsd', '1lse', '1lsf', '1lyo', '1lys', '1lyz', '1lza', '1lzt', '1ps5', '1rfp', '1uco', '1uig', '1v7s', '1vdq', '1vds', '1vdt', '1ved', '1xei', '1xej', '1xek', '2a6u', '2aub', '2c8o', '2c8p', '2cds', '2d4j', '2epe', '2f2n', '2hs7', '2hs9', '2hso', '2lym', '2lyz', '2vb1', '2yvb', '3lym', '3lyt', '3lyz', '3lzt', '4lym', '4lyo', '4lyt', '4lyz', '4lzt', '5lym', '5lyt', '5lyz', '6lyt', '6lyz', '7lyz', '8lyz',
            ## with hydrogen
            '2vb1',
##            ## atoms missing
##            '1vdp',
            ]

        ##for pdb in l_wts:
        ##    os.system('cp /oxygenase_local/data/pdb/%s/pdb%s.ent %s.pdb' %(pdb[1:3],pdb,pdb))
        ##os.system('cp /software/whatif_debugged/dbdata/TOPOLOGY.H .')

        ## set 3/1 aa dic
        d_resname2ressymbol = {
            'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
            'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
            'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
            'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'
            }
        d_ressymbol2resname = {
            }
        for res_name in d_resname2ressymbol.keys():
            res_symbol = d_resname2ressymbol[res_name]
            d_ressymbol2resname[res_symbol] = res_name
        self.d_ressymbol2resname = d_ressymbol2resname

        ## parse sequence
        fd = open('/oxygenase_local/data/pdb/lz/pdb2lzt.ent','r')
        lines = fd.readlines()
        fd.close()
        l_seq = []
        for line in lines:
            record = line[:6].strip()
            if record == 'SEQRES':
                l_seq += line[19:].split()

        ## prepare data dictionary
        d_data = {}
        for i in range(len(l_seq)):
            res_name = l_seq[i]
            if res_name in ['ASP','GLU','ARG','LYS','HIS',]:
                titgrp = '%s%i' %(res_name,i+1)
                d_data[titgrp] = {}
                for j in range(len(l_seq)):
                    residue = d_resname2ressymbol[l_seq[j]]+str(j+1)
                    d_data[titgrp][residue] = []

        l_xtics = []
        for j in range(len(l_seq)):
            residue = d_resname2ressymbol[l_seq[j]]+str(j+1)
            l_xtics += [residue]

        self.d_data = d_data
        self.l_wts = l_wts
        self.l_xtics = l_xtics

        return


if __name__ == '__main__':
    delta_cs().main(
        delta_cs().l_wts,
        delta_cs().d_data,
        delta_cs().l_xtics,
        )

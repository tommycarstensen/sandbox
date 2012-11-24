## add new column for BGM lysine mutant effects - global unfold (green), more dynamic (red), no change (cyan), structural reorganization (blue)

import os, urllib2
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
import statistics

url = urllib2.urlopen('http://peat.ucd.ie/pkacoop/targets.csv')
s = url.read().lower()
l_proteins = s.split()

s_path = '/local/tc/pkacoop/predictions'
##s_path = '/home/people/farrell/python/damien_sandbox/pkacoop/predictions'
##s_path = '/var/www/html/pkacoop/predictions'
s_path = '/home/people/tc/svn//pkacoop/predictions'

s_alphabet = 'fghijkl'

s_sequence = 'ALA THR SER THR LYS LYS LEU HIS LYS GLU PRO ALA THR LEU ILE LYS ALA ILE ASP GLY ASP THR VAL LYS LEU MET TYR LYS GLY GLN PRO MET THR PHE ARG LEU LEU LEU VAL ASP THR PRO GLU THR LYS HIS PRO LYS LYS GLY VAL GLU LYS TYR GLY PRO GLU ALA SER ALA PHE THR LYS LYS MET VAL GLU ASN ALA LYS LYS ILE GLU VAL GLU PHE ASP LYS GLY GLN ARG THR ASP LYS TYR GLY ARG GLY LEU ALA TYR ILE TYR ALA ASP GLY LYS MET VAL ASN GLU ALA LEU VAL ARG GLN GLY LEU ALA LYS VAL ALA TYR VAL TYR LYS PRO ASN ASN THR HIS GLU GLN HIS LEU ARG LYS SER GLU ALA GLN ALA LYS LYS GLU LYS LEU ASN ILE TRP SER GLU ASP ASN ALA ASP SER GLY GLN'
l_sequence = s_sequence.split()

## 25 mutation sites
## mutations to lys,glu,asp,arg
## ((4-1)x25-8)+18 = 85 experimental values
d_experimental = {
    'mutants':{ ## wt, phs, delta+phs
         20:{'LYS':'>10.4','GLU':'<4.5','ASP':4.0,'ARG':'>10.4',}, ## 
         23:{'LYS':7.3,'GLU':7.1,'ASP':6.8,'ARG':'>10.4',}, ## v23k, v23e
         25:{'LYS':6.3,'GLU':7.5,'ASP':6.8,'ARG':'>10.4',}, ## 3erq, 3evq
         34:{'LYS':7.1,'GLU':7.3,'ASP':7.8,'ARG':'>10.4',}, ## f34k, f34e
         36:{'LYS':7.2,'GLU':8.7,'ASP':7.9,'ARG':'>10.4',}, ## 3eji, l36e
         37:{'LYS':'>10.4','GLU':5.2,'ASP':4.0,'ARG':'>10.4',}, ## l37e
         38:{'LYS':'>10.4','GLU':6.8,'ASP':6.8,'ARG':'>10.4',}, ## 3d6c
         39:{'LYS':9.0,'GLU':8.2,'ASP':8.1,'ARG':'>10.4',}, ## v39k, v39e
         41:{'LYS':9.3,'GLU':6.5,'ASP':4.0,'ARG':'>10.4',}, ## t41k, t41e
         58:{'LYS':'>10.4','GLU':7.7,'ASP':6.8,'ARG':'>10.4',}, ## a58e
         62:{'LYS':8.1,'GLU':7.7,'ASP':8.7,'ARG':'>10.4',}, ## 3dmu, t62e
         66:{'LYS':5.6,'GLU':8.5,'ASP':8.1,'ARG':'>10.4',}, ## 2snm, 1u9r
         72:{'LYS':8.6,'GLU':7.3,'ASP':7.6,'ARG':'>10.4',}, ## 2rbm, 3ero
         74:{'LYS':7.4,'GLU':7.8,'ASP':8.3,'ARG':'>10.4',}, ## v74k, v74e
         90:{'LYS':8.6,'GLU':6.4,'ASP':7.5,'ARG':'>10.4',}, ## a90k, a90e
         91:{'LYS':9.0,'GLU':7.1,'ASP':7.2,'ARG':'>10.4',}, ## y91k, 3d4d
         92:{'LYS':5.3,'GLU':9.0,'ASP':8.1,'ARG':'>10.4',}, ## 1tt2, 1tr5/1tqo
         99:{'LYS':6.5,'GLU':8.4,'ASP':8.5,'ARG':'>10.4',}, ## v99k, v99e
        100:{'LYS':8.6,'GLU':7.6,'ASP':6.9,'ARG':'>10.4',}, ## n100k, n100e
        103:{'LYS':8.2,'GLU':8.9,'ASP':8.7,'ARG':'>10.4',}, ## 3e5s, l103e
        104:{'LYS':7.7,'GLU':9.4,'ASP':9.7,'ARG':'>10.4',}, ## 3c1f, v104e
        109:{'LYS':9.2,'GLU':7.9,'ASP':7.5,'ARG':'>10.4',}, ## a109k, a109e
        118:{'LYS':'>10.4','GLU':'<4.5','ASP':7.0,'ARG':'>10.4',}, ## 
        125:{'LYS':6.2,'GLU':9.1,'ASP':7.6,'ARG':'>10.4',}, ## 3c1e, l125e
        132:{'LYS':'>10.4','GLU':7.0,'ASP':7.0,'ARG':'>10.4',}, ## a132e
         },
    'phs':{},
    'delta+phs':{
        ## Molecular determinants of the pKa values of
        ## The pKa Values of Acidic and Basic Residues Buried at
        ## 0.1M KCl
         19:{'ASP':2.21,},
         21:{'ASP':6.54,},
         40:{'ASP':3.87,},
         77:{'ASP':'<2.2',},
         83:{'ASP':'<2.2',},
         95:{'ASP':2.16,},
        143:{'ASP':3.80,},
        146:{'ASP':3.86,},
         10:{'GLU':2.82},
         43:{'GLU':4.32},
         52:{'GLU':3.93},
         57:{'GLU':3.49},
         67:{'GLU':3.76},
         73:{'GLU':3.31},
         75:{'GLU':3.26},
        101:{'GLU':3.81},
        122:{'GLU':3.89},
        129:{'GLU':3.75},
        135:{'GLU':3.76},
        142:{'GLU':4.49}, ## missing in pdb
        },
    'wt':{
        ## Distance dependence and salt sensitivity of pairwise,
        ## Electrostatic Effects in a Network of Polar and Ionizable
        ## Electrostatic Effects in Highly Charged Proteins: Salt Sensitivity of pKa Values of Histidines in Staphylococcal Nuclease
        ## 0.1M KCl
          8:{'HIS':6.52},
         46:{'HIS':5.86},
        121:{'HIS':5.30},
        124:{'HIS':5.73},
          },
    }

l_allowed = [
    ## mutants
    'LYS66','ASP66','GLU66','LYS92','ASP92','GLU92','LYS38','ASP38','GLU38',
    'LYS25','GLU25','LYS36','LYS62','LYS72','GLU72','ARG72','ARG90','GLU91','LYS103','LYS104','ARG109','LYS125',
    ## wt
    'HIS8','HIS46','HIS121','HIS124',
    ## delta+phs
    'ASP19','ASP21','ASP40','ASP77','ASP83','ASP95','ASP143','ASP146',
    'GLU10','GLU43','GLU52','GLU57','GLU67','GLU73','GLU75','GLU101','GLU122','GLU129','GLU135','GLU142',
    ]

## model pka values from wikipedia
d_null = {
    'ASP':3.8,
    'GLU':4.3,
    'LYS':10.5,
    'ARG':12.0,
    'HIS':6.08,
    }

d_nonmutants = {
    '1stn':'wt',
    '1snc':'wt',
    '3bdc':'delta+phs',
    '1ey8':'phs',
    }

d_parents = { ## non delta+phs
    '1stn':'wt',
    '1snc':'wt',
##    '3bdc':'delta+phs',
    '1ey8':'phs',
    '2rks':'phs',
    '3d6c':'phs',
    '2snm':'wt',
    'v66k_phs':'phs',
    '1u9r':'phs',
    '2oxp':'phs',
    '3dmu':'phs',
    '3c1e':'phs',
    '2pw5':'phs',
    '2pw7':'phs',
    '2pzw':'phs',
    '2pzu':'phs',
    '2pyk':'phs',
    '2pzt':'phs',
    '2qdb':'NVIAGA',
    '2rdf':'NVIAGA',
    }

min_x = 2 ## 5 if excl. Asp, 2 if incl. Asp
max_x = 11
min_y = -10
max_y = 20

min_x_shift = -7
max_x_shift = 7
min_y_shift = -15 ## excludes Gernot
max_y_shift = 15

## manually sorted list of PIs
l_PIs = [
    ## pka values
    'Jim Warwicker_A',    ## +wt, ?
    'Jim Warwicker_B',    ## +wt, ?
    'Francesca Milletti', ## -wt, all mutants
    'Wei Yang',           ## -wt, some mutants
    'Qiang Cui',          ## -wt, some mutants
    'Ernest Mehler',      ## -wt, some mutants, continuum (dielectric=4)
    ## tit curves
    'Gernot Kieseritzky', ## +wt, all mutants, PBE based (Carlsberg+, APBS)
    'JanJensen_Jdec19',   ## +wt, all mutants, empirical method
    'JanJensen_Wdec19',   ## +wt, all mutants, empirical method
    'JanJensen_Jjan15',   ## +wt, all mutants, empirical method
    'JanJensen_Wjan15',   ## +wt, all mutants, empirical method
    'WHATIF',             ## +wt, all mutants, PBE based (hydrogen bond optimization, delphi, dielectric=8)
    'PDB2PKA',            ## +wt, all mutants, PBE based
    'Emil Alexov_1',      ## -wt, all mutants
    'Emil Alexov_2',      ## -wt, all mutants
    'Yifan Song',         ## -wt, all mutants
    'Jana K. Shen',       ## -wt, all mutants, MD based
    'Mike Word',          ## +wt, ?
    'Cat Chenal',         ## +wt, some mutants
    'Antonio Baptista',   ## -wt, some mutants, MD based
    'Sarah Williams',     ## +wt, no mutants
    ]
l_PIs.sort()

d_PIs = {
    'Antonio Baptista':{'submission':'titcurv','short':'Antonio','long':'Antonio Baptista',},
    'Cat Chenal':{'submission':'titcurv','short':'Cat','long':'Cat Chenal',},
    'Emil Alexov_1':{'submission':'titcurv','short':'Emil1','long':'Emil Alexov',},
    'Emil Alexov_2':{'submission':'titcurv','short':'Emil2','long':'Emil Alexov',},
    'Ernest Mehler':{'submission':'pka','short':'Ernest','long':'Ernest Mehler',},
    'Francesca Milletti':{'submission':'pka','short':'Francesca','long':'Francesca Milletti',},
    'Gernot Kieseritzky':{'submission':'titcurv','short':'Gernot','long':'Gernot Kieseritzky',},
    'Jana K. Shen':{'submission':'titcurv','short':'Jana','long':'Jana K. Shen',},
    'JanJensen_Jdec19':{'submission':'titcurv','short':'Jan1','long':'Jan Jensen',},
    'JanJensen_Wdec19':{'submission':'titcurv','short':'Jan2','long':'Jan Jensen',},
    'JanJensen_Jjan15':{'submission':'titcurv','short':'Jan3','long':'Jan Jensen',},
    'JanJensen_Wjan15':{'submission':'titcurv','short':'Jan4','long':'Jan Jensen',},
    'Jim Warwicker_A':{'submission':'pka','short':'JimA','long':'Jim Warwicker',},
    'Jim Warwicker_B':{'submission':'pka','short':'JimB','long':'Jim Warwicker',},
    'Mike Word':{'submission':'titcurv','short':'Mike','long':'Mike Word',},
    'PDB2PKA':{'submission':'titcurv','short':'PDB2PKA','long':'PDB2PKA',},
    'Qiang Cui':{'submission':'pka','short':'Qiang','long':'Qiang Cui',},
    'Sarah Williams':{'submission':'titcurv','short':'Sarah','long':'Sarah Williams',},
    'Wei Yang':{'submission':'pka','short':'Wei','long':'Wei Yang',},
    'WHATIF':{'submission':'titcurv','short':'WHATIF','long':'WHATIF',},
    'Yifan Song':{'submission':'titcurv','short':'Yifan','long':'Yifan Song',},
    }

d_PIs['PDB2PKA']['color'] = 'FF0000' ## red
d_PIs['PDB2PKA']['pt'] = 5
d_PIs['WHATIF']['color'] = 'FF8000' ## orange
d_PIs['WHATIF']['pt'] = 7
d_PIs['Gernot Kieseritzky']['color'] = 'FFFF00' ## yellow
d_PIs['Gernot Kieseritzky']['pt'] = 5
d_PIs['JanJensen_Jdec19']['color'] = '80FF00' ## green1
d_PIs['JanJensen_Jdec19']['pt'] = 7
d_PIs['JanJensen_Wdec19']['color'] = '80FF00' ## green1
d_PIs['JanJensen_Wdec19']['pt'] = 7
d_PIs['JanJensen_Jjan15']['color'] = '80FF00' ## green1
d_PIs['JanJensen_Jjan15']['pt'] = 7
d_PIs['JanJensen_Wjan15']['color'] = '80FF00' ## green1
d_PIs['JanJensen_Wjan15']['pt'] = 7
d_PIs['Jana K. Shen']['color'] = '00FF00' ## green2
d_PIs['Jana K. Shen']['pt'] = 5
d_PIs['Emil Alexov_1']['color'] = '00FF80' ## green3
d_PIs['Emil Alexov_1']['pt'] = 7
d_PIs['Emil Alexov_2']['color'] = '00FF80' ## green3
d_PIs['Emil Alexov_2']['pt'] = 7
d_PIs['Sarah Williams']['color'] = '00FFFF' ## cyan
d_PIs['Sarah Williams']['pt'] = 5
d_PIs['Cat Chenal']['color'] = '0080FF' ## light blue
d_PIs['Cat Chenal']['pt'] = 7
d_PIs['Mike Word']['color'] = '9F9F9F' ## grey 3
d_PIs['Mike Word']['pt'] = 7
d_PIs['Qiang Cui']['color'] = 'FF0080' ## purple
d_PIs['Qiang Cui']['pt'] = 7
d_PIs['Antonio Baptista']['color'] = '7F7F7F' ## grey 4
d_PIs['Antonio Baptista']['pt'] = 5
d_PIs['Yifan Song']['color'] = 'DFDFDF' ## grey 1
d_PIs['Yifan Song']['pt'] = 7
d_PIs['Jim Warwicker_A']['color'] = '0000FF' ## blue
d_PIs['Jim Warwicker_A']['pt'] = 5
d_PIs['Jim Warwicker_B']['color'] = '0000FF' ## blue
d_PIs['Jim Warwicker_B']['pt'] = 5
d_PIs['Francesca Milletti']['color'] = '9000FF' ## violet
d_PIs['Francesca Milletti']['pt'] = 7
d_PIs['Wei Yang']['color'] = 'BFBFBF' ## grey 2
d_PIs['Wei Yang']['pt'] = 5
d_PIs['Ernest Mehler']['color'] = 'FF00FF' ## pink
d_PIs['Ernest Mehler']['pt'] = 5

for PI in d_PIs.keys():
    if d_PIs[PI]['pt'] == 5:
        d_PIs[PI]['ul'] = 'square'
    elif d_PIs[PI]['pt'] == 7:
        d_PIs[PI]['ul'] = 'disc'
    else:
        print d_PIs[PI]['pt']
        stop
del PI

d_bgm_residues = {
    '1stn':['HIS',],
    '1snc':[],
    '3bdc':['HIS','ASP','GLU',],
    '1ey8':[],
    '2rks':['LYS:A:0038',],
    '3d6c':['GLU:A:0038',],
    '2snm':['LYS:A:0066',],
    'v66k_phs':['LYS:A:0066',],
    'v66k_delta+phs':['LYS:A:0066',],
    '1u9r':['GLU:A:0066',],
    '2oxp':['ASP:A:0066',],
    '1tt2':['LYS:A:0092',],
    '1tr5':['GLU:A:0092',],
    '1tqo':['GLU:A:0092',],
    '2oeo':['ASP:A:0092',],
    '3erq':['LYS:A:0025',],
    '3evq':['GLU:A:0025',],
    '3eji':['LYS:A:0036',],
    '3dmu':['LYS:A:0062',],
    '2rbm':['LYS:A:0072',],
    '3ero':['GLU:A:0072',],
    '3d8g':['ARG:A:0072',],
    '3dhq':['ARG:A:0090',],
    '3d4d':['GLU:A:0091',],
    '3e5s':['LYS:A:0103',],
    '3c1f':['LYS:A:0104',],
    '3d4w':['ARG:A:0109',],
    '3c1e':['LYS:A:0125',],
    '2pw5':['TYR:A:0066',],
    '2pw7':['TYR:A:0066',],
    '2pzw':['ASN:A:0066',],
    '2pzu':['ASN:A:0066',],
    '2pyk':['GLN:A:0066',],
    '2pzt':['GLN:A:0066',],
    '2rdf':['HIS',],
    '2qdb':['HIS',],
    }

d_mutations = {
    '2rks':'LYS:A:0038',
    '3d6c':'GLU:A:0038',
    '2snm':'LYS:A:0066',
    'v66k_phs':'LYS:A:0066',
    'v66k_delta+phs':'LYS:A:0066',
    '1u9r':'GLU:A:0066',
    '2oxp':'ASP:A:0066',
    '1tt2':'LYS:A:0092',
    '1tr5':'GLU:A:0092',
    '1tqo':'GLU:A:0092',
    '2oeo':'ASP:A:0092',
    '3erq':'LYS:A:0025',
    '3evq':'GLU:A:0025',
    '3eji':'LYS:A:0036',
    '3dmu':'LYS:A:0062',
    '2rbm':'LYS:A:0072',
    '3ero':'GLU:A:0072',
    '3d8g':'ARG:A:0072',
    '3dhq':'ARG:A:0090',
    '3d4d':'GLU:A:0091',
    '3e5s':'LYS:A:0103',
    '3c1f':'LYS:A:0104',
    '3d4w':'ARG:A:0109',
    '3c1e':'LYS:A:0125',
    '2pw5':'TYR:A:0066',
    '2pw7':'TYR:A:0066',
    '2pzw':'ASN:A:0066',
    '2pzu':'ASN:A:0066',
    '2pyk':'GLN:A:0066',
    '2pzt':'GLN:A:0066',
    '2rdf':'GLN:A:0075',
    '2qdb':'ALA:A:0075',
    }

## point one mutant to multiple pdbs
d_yifan = {}
for k,v in d_mutations.items():
    if not v in d_yifan.keys():
        d_yifan[v] = []
    d_yifan[v] += [k]

d_res1res3 = {
    'd':'ASP',
    'e':'GLU',
    'k':'LYS',
    'r':'ARG',
    }

d_res3res1 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
    'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
    'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
    'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
    }

d_pdbs = {}
for protein in l_proteins:
    if len(protein) == 4:
        pdb = '%s' %(protein,)
    else:
        pdb = '3bdc_%s' %(protein[:protein.index('_')],)
    d_pdbs[protein] = pdb

f_fulltit = 0.05


class Telluride:

    def main(self):

##        os.system('sshfs tc@peat:/var/www/html/pkacoop/predictions predictions/')

        ##
        ## parse data
        ##
        d_pkas = self.parse_previous_fits()

        d_titcurves, d_pkas, d_residues = self.parse_data(d_pkas)

        d_slopes = self.calculate_slopes(d_titcurves)

        ##
        ## after parsing data
        ##
        d_pkas = self.fit_and_plot_titcurves_individual(d_titcurves,d_pkas,d_residues,)

        ##
        ## after calculation of pka
        ##
        d_rmsd, d_differences, d_counts = self.calculate_rmsds(d_pkas,d_residues,)

        self.plot_titcurves_combined(d_residues, d_titcurves,d_pkas,d_differences,d_slopes,)

        self.plot_histograms(d_pkas,d_rmsd,d_differences,)

        lines_htm_map, lines_htm_map_restricted, lines_htm_map_shift, lines_htm_map_shift_restricted = self.plot_exp_pka_vs_calc_pka(d_residues, d_pkas,)

        ##
        ## after calculation of statistics (differences, rmsds)
        ##
        self.write_htm(d_pkas,d_residues,lines_htm_map,d_rmsd,d_differences,d_counts,lines_htm_map_restricted,lines_htm_map_shift,lines_htm_map_shift_restricted,)

        self.copy_to_server()


        os.system('rm *.dat')

##        os.system('fusermount -u predictions/')

##        ##                    
##        ## check proteins
##        ##
##        for PI in d_titcurves.keys():
##            ## anything that's not supposed to be there?
##            if len(set(d_titcurves[PI].keys())-set(l_proteins)) > 0:
##                print PI, set(d_titcurves[PI].keys())-set(l_proteins)
##                stop
##            print 'not submitted', PI, set(l_proteins)-set(d_titcurves[PI].keys())

        ##
        ## check titgrps
        ##
        for PI in d_titcurves.keys():
            for protein in l_proteins:
                if not protein in d_titcurves[PI].keys():
                    continue
                if protein in ['v66k_phs','v66k_delta+phs',]:
                    PI_ref = 'Francesca Milletti'
                    PI_ref = 'JanJensen_Wjan15'
                elif protein in ['3bdc',]:
                    PI_ref = 'PDB2PKA'
                else:
                    PI_ref = 'PDB2PKA'
                set1 = set(d_titcurves[PI_ref][protein].keys())-set(d_titcurves[PI][protein].keys())
                set2 = set(d_titcurves[PI][protein].keys())-set(d_titcurves[PI_ref][protein].keys())
    ##            if len(set1) > 0
    ##                print PI, protein, set1
                if len(set2) > 0:
                    l1 = list(set1)
                    l1.sort()
                    l2 = list(set2)
                    l2.sort()
                    print
                    print PI_ref, d_titcurves[PI_ref][protein].keys()
                    print PI, d_titcurves[PI][protein].keys()
                    print PI
                    print protein
                    print PI_ref, l1
                    print PI, l2 ## PI not PI_ref
                    stop

##        ##
##        ## update damiens dictionary
##        ##
##        import pickle
##        d_damien = pickle.load(open('pkadict.obj','r'))
##        for protein in l_proteins:
##            if len(protein) > 4:
##                protein_damien = protein.upper().replace('DELTA','delta',)
##            else:
##                protein_damien = protein
##            for PI in l_PIs:
##                if not PI in d_damien[protein_damien].keys():
##                    d_damien[protein_damien][PI] = {}
##                if not protein in d_pkas[PI].keys():
##                    continue
##                for res_ID in d_residues[protein]:
##                    if not res_ID in d_pkas[PI][protein].keys():
##                        continue
##                    pka_gnuplot = d_pkas[PI][protein][res_ID]['pka1']
##                    ## append
##                    if not res_ID in d_damien[protein_damien][PI].keys():
##                        if pka_gnuplot != 'N/A':
##                            pka_gnuplot = float(pka_gnuplot)
##                        d_damien[protein_damien][PI][res_ID] = pka_gnuplot
##                    ## overwrite
##                    else:
##                        pka_damien = d_damien[protein_damien][PI][res_ID]
##                        if pka_damien in [7.0,'N/A',] and pka_gnuplot == 'N/A':
##                            continue
##                        else:
##                            print PI,protein,res_ID,pka_damien,pka_gnuplot
##                            pka_gnuplot = float(d_pkas[PI][protein][res_ID]['pka1'])
##                        ## don't overwirte if large difference
##                        if abs(pka_gnuplot-pka_damien) > 1:
##                            print PI, protein, res_ID, pka_gnuplot, pka_damien
##                            continue
##                        d_damien[protein_damien][PI][res_ID] = pka_gnuplot
##        pickle.dump(d_damien,open('pkadict.obj','w'))
##        os.system('scp pkadict.obj tc@peat:/var/www/html/pkacoop/.')

        return


    def calculate_slopes(self,d_titcurves,):

        import pickle

        print
        print 'calculating slopes'

        if os.path.isfile('slopes.pickle'):

            d_slopes = pickle.load(open('slopes.pickle','r'))

        else:

            import pKarun.pKa_general
            instance_pKanalyse = pKarun.pKa_general.pKanalyse()

            d_slopes = {}

            for PI in d_titcurves.keys():
                d_slopes[PI] = {}
                for protein in d_titcurves[PI].keys():
                    d_slopes[PI][protein] = {}
                    for res_ID in d_titcurves[PI][protein].keys():
                        d = d_titcurves[PI][protein][res_ID]
##                        for pH in d.keys():
##                            if pH < 0:
##                                del d[pH]
                        solution,sq = instance_pKanalyse.fit_to_henderson(d)
                        slope = float(solution[0])
                        d_slopes[PI][protein][res_ID] = slope

            pickle.dump(d_slopes,open('slopes.pickle','w'))

        return d_slopes


    def copy_to_server(self,):

        print
        print 'scp files to server'

        for folder in ['analysis','analysis_restricted',]:

            print 'scp', folder
        
##            os.system('scp -r pdb tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,)) ## temp!!!

            os.system('scp -r jmol_scripts tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,))

            if folder == 'analysis':
                for PI in l_PIs:
                    os.system('scp %s.htm tc@peat:/var/www/html/pkacoop/%s/%s/. > out.txt' %(PI.replace(' ',''),folder,d_PIs[PI]['long'].replace(' ',''),))
                    os.remove('%s.htm' %(PI.replace(' ',''),))

            ## htm
            if folder == 'analysis':
                os.system('scp *.htm gnuplot.png gnuplot_shift.png tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,))
            else:
                os.system('scp analysis_PI_restricted.htm analysis_residue_restricted.htm analysis_residue_restricted_experimental.htm gnuplot.png gnuplot_shift.png tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,))

            ## histograms
            if folder == 'analysis':
                os.system('scp -r png_histograms htm_histograms htm_PI tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,))
            elif folder == 'analysis_restricted':
                os.system('scp -r png_histograms htm_histograms_restricted htm_PI tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,))

            ## png_titcurves/combined
            if folder == 'analysis':
                print 'scp', folder, 'png_titcurves/combined'
                os.system('scp -r png_titcurves/combined tc@peat:/var/www/html/pkacoop/%s/png_titcurves/. > out.txt' %(folder,))
            else:
                print 'scp', folder, 'png_titcurves/combined_restricted'
                os.system('scp -r png_titcurves/combined_restricted tc@peat:/var/www/html/pkacoop/%s/png_titcurves/. > out.txt' %(folder,))

            ##
            ## htm_titcurves
            ##
            if folder == 'analysis':
                os.system('scp -r htm_titcurves tc@peat:/var/www/html/pkacoop/%s/. > out.txt' %(folder,))
            elif folder == 'analysis_restricted':
                ## histograms
                os.system('scp -r htm_titcurves/histograms_restricted tc@peat:/var/www/html/pkacoop/%s/htm_titcurves/. > out.txt' %(folder,))
                ## combined
                os.system('scp -r htm_titcurves/combined_restricted tc@peat:/var/www/html/pkacoop/%s/htm_titcurves/. > out.txt' %(folder,))
                ## per PI
                if '-scpall' in sys.argv:
                    for PI in l_PIs:
                        print 'scp', folder, 'htm_titcurves', PI
                        l = os.listdir(
                            'htm_titcurves/%s' %(PI)
                            )
                        for fn in l:
                            res_no = int(fn[-8:-4])
                            res_name = fn[-14:-11]
                            experimental = False
                            if res_no in d_experimental['mutants'].keys():
                                if res_name in d_experimental['mutants'][res_no].keys():
                                    experimental = True
                            if experimental == True and not '%s%i' %(res_name,res_no,) in l_allowed:
                                continue
                            s = 'scp ./htm_titcurves/%s/%s tc@peat:"/var/www/html/pkacoop/%s/htm_titcurves/%s/." > out.txt' %(
                                PI.replace(' ','\ ',), fn, folder, PI.replace(' ','\ ',),
                                )
                            os.system(s)

            ##
            ## png_titcurves
            ##
            if folder == 'analysis':
                if '-scpall' in sys.argv:
                    for PI in l_PIs:
                        print 'scp png_titcurves', folder, PI
                        os.system(
                            'scp -r png_titcurves/%s tc@peat:/var/www/html/pkacoop/%s/png_titcurves/. > out.txt' %(
                                PI.replace(' ','\ ',), folder,
                                )
                            )
            elif folder == 'analysis_restricted':
                if '-scpall' in sys.argv:
                    for PI in l_PIs:
                        l = os.listdir(
                            'png_titcurves/%s' %(PI)
                            )
                        print 'scp png_titcurves', folder, PI
                        for fn in l:
                            res_no = int(fn[-8:-4])
                            res_name = fn[-14:-11]
                            experimental = False
                            if res_no in d_experimental['mutants'].keys():
                                if res_name in d_experimental['mutants'][res_no].keys():
                                    experimental = True
                            if experimental == True and not '%s%i' %(res_name,res_no,) in l_allowed:
                                continue
                            s = 'scp ./png_titcurves/%s/%s tc@peat:"/var/www/html/pkacoop/%s/png_titcurves/%s/." > out.txt' %(
                                PI.replace(' ','\ ',), fn, folder, PI.replace(' ','\ ',),
                                )
                            os.system(s)

        return


    def calculate_rmsds(self,d_pkas,d_residues,):
        
        set_th = set([
            'HH','pka1','1or2pka',
            'wt','phs','delta+phs','mutants',
            'ASP','GLU','LYS','HIS','all',
            'all_pdb','mutant_pdb','mutant_model',
            ])

        d_differences = {}
        d_differences2 = {}
        d_counts = {}
        for PI in l_PIs+['null']:
            n = 0
            d_differences[PI] = {}
            d_counts[PI] = {}
            for k in ['all','>1','>1.5','>2','>3','>4',]:
                d_counts[PI][k] = 0
            for k in set_th:
                d_differences[PI][k] = []

        ##
        ## 1a) calculate differences per PI
        ##
        for PI in l_PIs:

            if not PI in d_differences2.keys():
                d_differences2[PI] = {}

            for protein in l_proteins:

                if not protein in d_pkas[PI].keys():
                    continue

                if protein in d_nonmutants.keys():
                    protein_type = d_nonmutants[protein]
                else:
                    protein_type = 'mutants'

                for res_ID in d_pkas[PI][protein].keys():

                    ## skip residue?
                    bool_skip = self.skip_residue(protein,res_ID,)
                    if bool_skip == True:
                        continue
    
                    res_no  = int(res_ID[-4:])
                    res_name = res_ID[:3]

                    ## skip if non-mutated residue in mutant protein
                    if protein_type == 'mutants' and res_name == l_sequence[res_no-1]:
                        continue

                    ## skip if no experimental data
                    if not res_no in d_experimental[protein_type].keys():
                        continue
                    ## skip if no experimental data
                    if not res_name in d_experimental[protein_type][res_no].keys():
                        continue
                    pka_exp = d_experimental[protein_type][res_no][res_name]
                    ## skip if limit value
                    if type(pka_exp).__name__ == 'str':
                        continue

                    d_counts[PI]['all'] += 1

##                    if PI == 'JanJensen_Jdec19' and res_no in [143,142,146,21,]:
##                        print protein, res_ID, protein_type, res_no, res_name, s_pka
##                        stop

                    s_pka = d_pkas[PI][protein][res_ID]['s_pka']
                    pka1 = d_pkas[PI][protein][res_ID]['pka1']

                    ## incomplete titration (not titrating within pH range)
                    if s_pka == 'N/A':
                        continue

                    ## HH titration
                    if not ';' in s_pka:
                        d_differences[PI]['HH'] += [float(s_pka)-pka_exp]

                    ## 1 pka (HH and non-HH titration)
                    diff = float(pka1)-pka_exp
                    d_differences[PI]['pka1'] += [diff]
                    d_differences[PI]['all'] += [diff]
                    for diff_limit in [1,1.5,2,3,4,]:
                        if abs(diff) > diff_limit:
                            d_counts[PI]['>%s' %(diff_limit)] += 1

                    ## 1 (HH) or 2 (non-HH) pka
                    diff_min = 99
                    for pka_calc in s_pka.split(';'):
                        diff = float(pka_calc)-pka_exp
                        if abs(diff) < abs(diff_min):
                            diff_min = diff
                    d_differences[PI]['1or2pka'] += [diff_min]

                    ## res_name
                    d_differences[PI][res_name] += [diff_min]

                    ## parent
                    d_differences[PI][protein_type] += [diff_min]

                    if protein_type == 'mutants':
                        if len(protein) > 4:
                            d_differences[PI]['mutant_model'] += [diff]
                        else:
                            d_differences[PI]['mutant_pdb'] += [diff]
                            d_differences[PI]['all_pdb'] += [diff]
                    else:
                        if len(protein) != 4:
                            stop
                        d_differences[PI]['all_pdb'] += [diff]

                    if not protein in d_differences2[PI].keys():
                        d_differences2[PI][protein] = {}
                    d_differences2[PI][protein][res_ID] = diff_min

        ##
        ## 1b) calculate differences for null model
        ##
        d_counts['null'] = {}
        for k in ['all','>1','>1.5','>2','>3','>4',]:
            d_counts['null'][k] = 0
        for protein_type in d_experimental.keys():
            for res_no in d_experimental[protein_type].keys():
                for res_name in d_experimental[protein_type][res_no].keys():
                    pka_exp = d_experimental[protein_type][res_no][res_name]
                    ## skip if limit
                    if type(pka_exp).__name__ == 'str':
                        continue
                    ## calc diff
                    diff = pka_exp-d_null[res_name]
                    d_counts['null']['all'] += 1

                    ## all
                    d_differences['null']['all'] += [diff]
                    ## diff limits
                    for diff_limit in [1,1.5,2,3,4,]:
                        if abs(diff) > diff_limit:
                            d_counts['null']['>%s' %(diff_limit)] += 1
                    ## all residues
                    d_differences['null']['1or2pka'] += [diff]

                    ## res_name
                    d_differences['null'][res_name] += [diff_min]

                    ## parent
                    d_differences['null'][protein_type] += [diff_min]


        ##
        ## 2) calculate rmsd
        ##
        d_rmsd = {}
        for PI in l_PIs+['null']:
            d_rmsd[PI] = {}
            for k in set_th:
                if len(d_differences[PI][k]) > 1:
                    rmsd = statistics.rmsd(d_differences[PI][k])
                    rmsd = '%.1f' %(rmsd)
                else:
                    rmsd = 'N/A'
                d_rmsd[PI][k] = {'rmsd':rmsd,'n':len(d_differences[PI][k])}

        return d_rmsd, d_differences2, d_counts


    def htm_per_PI(
        self, d_counts, d_rmsd, prefix, folder,
        lines_htm_map_shift, lines_htm_map,
        ):

        ## loop over PIs (rows)
        l_th = ['all','wt','phs','delta+phs','mutants','ASP','GLU','HIS','LYS','all_pdb','mutant_pdb','mutant_model',]

        ##
        ## htm header
        ##
        l_header = []
        l_header += ['<html>\n<head>\n']
        ## javascript
        l_header += ['<script type="text/javascript" src="../sortable.js"></script>\n']
        ## style sheet
        l_header += ['<link href="../table.css" rel="stylesheet" type="text/css">\n']
        l_header += ['<title>Telluride 2009</title>\n</head>\n']

        ##
        ## index htm
        ##
        l_htm = ['<body>\n']
        l_htm += ['<a href="analysis_PI.htm">Analysis per research group</a><br>\n']
        l_htm += ['<a href="analysis_residue_simple.htm">Analysis per residue (statistics only)</a><br>\n']
        l_htm += ['<a href="analysis_residue_advanced.htm">Analysis per residue (statistics and all submitted pKa values)</a><br>\n']
        l_htm = l_header+l_htm
        fd = open('index.htm','w')
        fd.writelines(l_htm)
        fd.close()


        ##
        ## initiate htm per PI
        ##
        l_htm = ['<body>\n']
        l_htm = ['All titration curves have been fitted to 1 pKa value.<br>\n']
        l_htm = ['A few residues do not titrate within the reported pH range. No pKa value has been calculated for these few residues.<br>\n']
        l_htm = ['Residues for which only the lower or upper limit of the pKa value has been determined are not used for calculating RMSDs.<br>\n']
        s = 'For the null model the following values have been used: '
        for k,v in d_null.items():
            s += '%1s%2s %s, ' %(k[:1].upper(),k[1:3].lower(),v,)
        s = s[:-2]+'\n'
        l_htm += [s]
        l_htm = ['RMSD of "mutant pdb" and "mutant model" refers to the RMSD of mutants with and without available PDB structures respectively.<br>\n']

        ##
        ## table of rmsds for experimental data
        ##
        l_htm += ['<table border="0" class="sortable" id="sortable">\n']
        l_htm += ['<tr>\n']
        l_htm += ['<td>name</td>\n']
        l_htm += ['<td>submitted</td>\n']
        l_diff_limits = [1,1.5,2,3,4,]
        for diff_limit in l_diff_limits:
            l_htm += ['<td>n (&Delta;pK<lower>a</lower> > %s)</td>\n' %(diff_limit,)]
        l_htm += ['<td>rmsd (all)</td>\n']
        l_htm += ['<td><a href="%s/wt.htm">rmsd (wt)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/phs.htm">rmsd (PHS)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/delta+phs.htm">rmsd (&Delta;+PHS)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/mutants.htm">rmsd (mutants)</a></td>\n' %(folder)]
        for res_name in ['ASP','GLU','HIS','LYS',]:
            l_htm += ['<td><a href="%s/%s.htm">rmsd (%s)</a></td>\n' %(folder,res_name,res_name,)]
        l_htm += ['<td><a href="%s/allpdb.htm">rmsd (all_pdb)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/pdb.htm">rmsd (mutant_pdb)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/model.htm">rmsd (mutant_model)</a></td>\n' %(folder)]
        l_htm += ['</tr>\n']

        for PI in l_PIs+['null']:
            ## initiate table row
            l_htm += ['<tr>\n']
            ## name
            if PI == 'null':
                l_htm += ['<td>null model</td>\n']
            else:
                l_htm += ['<td><a href="htm_PI/%s.htm">%s</a></td>\n' %(PI.replace(' ',''),PI,)]
            ## n submissions
            l_htm += ['<td>%s</td>\n' %(d_counts[PI]['all'],)]
            for diff_limit in l_diff_limits:
                l_htm += ['<td>%02i</td>\n' %(d_counts[PI]['>%s' %(diff_limit)]),]
            ## loop over columns
            for k in l_th:
                rmsd = d_rmsd[PI][k]['rmsd']
                n = d_rmsd[PI][k]['n']
                if PI in l_PIs:
                    if n > 0:
                        l_htm += ['<td><a href="%s/%s_%s.htm">%s (n = %s)</a></td>\n' %(folder,k,PI,rmsd,n,)]
                    else:
                        l_htm += ['<td>&nbsp;</td>\n']
                else:
                    l_htm += ['<td>%s (n = %s)</td>\n' %(rmsd,n,)]
            l_htm += ['</tr>\n']

        ## links to histograms
        l_htm += ['<tr>\n']
        l_htm += ['<td>&nbsp;</td>\n']
        l_htm += ['<td>&nbsp;</td>\n']
        l_htm += ['<td>&nbsp;</td>\n']
        l_htm += ['<td>&nbsp;</td>\n']
        l_htm += ['<td>&nbsp;</td>\n']
        l_htm += ['<td>&nbsp;</td>\n']
        l_htm += ['<td>&nbsp;</td>\n']
##        l_htm += ['<td colspan="%i"></td>\n' %(2+len(l_diff_limits))]
        l_htm += ['<td><a href="%s/all.htm">rmsd (all)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/wt.htm">rmsd (wt)</a></td>\n' %(folder)]
        l_htm += ['<td>&nbsp;</td>\n'] ## PHS
        l_htm += ['<td><a href="%s/delta+phs.htm">rmsd (&Delta;+PHS)</a></td>\n' %(folder)]
        l_htm += ['<td><a href="%s/mutants.htm">rmsd (mutants)</a></td>\n' %(folder)]
        for res_name in ['ASP','GLU','HIS','LYS',]:
            l_htm += ['<td><a href="%s/%s.htm">rmsd (%s)</a></td>\n' %(folder,res_name,res_name,)]
        l_htm += ['<td>&nbsp;</td>\n'] ## all pdb
        l_htm += ['<td>&nbsp;</td>\n'] ## mutant pdb
        l_htm += ['<td>&nbsp;</td>\n'] ## mutant model
        l_htm += ['</tr>\n']

        l_htm += ['</table>\n']

        ##
        ## plot
        ##
        l_htm += ['\n\n<br><br>\n\n']
        l_htm += ['experimental pKa v predicted pKa<br>\n']
        l_htm += ['<img src="gnuplot.png" usemap="#Telluride" border="0"><br>\n']
        l_htm += ['experimental pKa shift v predicted pKa shift<br>\n']
        l_htm += ['<img src="gnuplot_shift.png" usemap="#Telluride_shift" border="0"><br>\n']
        l_htm += lines_htm_map
        l_htm += lines_htm_map_shift
        l_htm += ['\n\n<br><br>\n\n']

        ##
        ## terminate htm per PI
        ##
        l_htm += ['</body>\n</html>\n']

        l_htm = l_header+l_htm

        fd = open('%s.htm' %(prefix),'w')
        fd.writelines(l_htm)
        fd.close()

        return


    def write_htm(self,d_pkas,d_residues,lines_htm_map,d_rmsd,d_differences,d_counts,lines_htm_map_restricted,lines_htm_map_shift,lines_htm_map_shift_restricted,):

        prefix = 'analysis_PI'
        folder = 'htm_histograms'
        self.htm_per_PI(
            d_counts, d_rmsd, prefix, folder,
            lines_htm_map_shift, lines_htm_map, 
            )

        prefix = 'analysis_PI_restricted'
        folder = 'htm_histograms_restricted'
        self.htm_per_PI(
            d_counts, d_rmsd, prefix, folder,
            lines_htm_map_shift_restricted, lines_htm_map_restricted,
            )

        ## 1) statistics and all PIs
        prefix = 'analysis_residue_advanced'
        title = 'Telluride 2009'
        l_PIs_htm = l_PIs
        restricted = False
        limited = True
        self.htm_per_residue(
            l_PIs_htm, prefix, d_residues, d_differences, title,
            d_pkas, restricted, limited,
            )

        ## 2) statistics and one PI
        for PI in l_PIs:
            prefix = '%s' %(PI.replace(' ',''))
            title = PI
            l_PIs_htm = [PI]
            restricted = True
            limited = True
            self.htm_per_residue(
                l_PIs_htm, prefix, d_residues, d_differences, title,
                d_pkas, restricted, limited,
                )

        ## 3) statistics only
        prefix = 'analysis_residue_simple'
        title = 'Telluride 2009'
        l_PIs_htm = []
        restricted = False
        limited = True
        self.htm_per_residue(
            l_PIs_htm, prefix, d_residues, d_differences, title,
            d_pkas, restricted, limited,
            )

        ## 4) restricted dataset (experimental only)
        prefix = 'analysis_residue_restricted_experimental'
        title = 'Telluride 2009 - Restricted data set - Residues with experimental values only'
        l_PIs_htm = l_PIs
        restricted = True
        limited = True
        self.htm_per_residue(
            l_PIs_htm, prefix, d_residues, d_differences, title,
            d_pkas, restricted, limited,
            )

        ## 5) restricted dataset (all residues)
        prefix = 'analysis_residue_restricted'
        title = 'Telluride 2009 - Restricted data set - All residues'
        l_PIs_htm = l_PIs
        restricted = True
        limited = False
        self.htm_per_residue(
            l_PIs_htm, prefix, d_residues, d_differences, title,
            d_pkas, restricted, limited,
            )

        return


    def htm_per_residue(
        self, l_PIs_htm, prefix, d_residues, d_differences, title,
        d_pkas,
        restricted, limited,
        ):

        if restricted == True:
            folder = 'combined_restricted'
        else:
            folder = 'combined'

        if len(l_PIs_htm) == 1:
            path_htm_titcurve = '../../analysis_restricted/htm_titcurves'
        else:
            path_htm_titcurve = 'htm_titcurves'

        ##
        ## htm header
        ##
        l_header = []
        l_header += ['<html>\n<head>\n']
        if title in l_PIs:
            ## style sheet
            l_header += ['<link href="../../table.css" rel="stylesheet" type="text/css">\n']
            ## javascript
            l_header += ['<script type="text/javascript" src="../../sortable.js"></script>\n']
        else:
            ## style sheet
            l_header += ['<link href="../table.css" rel="stylesheet" type="text/css">\n']
            ## javascript
            l_header += ['<script type="text/javascript" src="../sortable.js"></script>\n']
        l_header += ['<title>%s</title>\n</head>\n' %(title)]

        ## statistics header titles
        l_columns = [
            'protein','mutant pdb','res_name','res_no','n','average','stddev','min','max',
            'pka exp','pka shift','min diff','max diff','rmsd',
            ]
##        d_columns = {
##            'n':'number of submitted/calculated pKa values for a residue',
##            'average':'average of calculated pKa values',
##            'stderr':'stderr of calculated pKa values (excl. Gernot, Emil_2, Jim_B, 3xJan)',
##            'min''min of calculated pKa values':,
##            'max':'max of calculated pKa values',
##            'pka exp':'experimental pKa value',
##            'pka shift':'shift of experimental pKa value from model pKa value',
##            'min diff':,
##            'max diff':,
##            'rmsd':,
##            }

        ##
        ## initiate htm
        ##
        l_htm = ['<body>\n']
        l_htm += ['n = number of submitted/calculated pKa values for a residue<br>\n']
        l_htm += ['average = average of calculated pKa values<br>\n']
        l_htm += ['stddev = stddev of calculated pKa values (excl. Gernot, Emil_2, Jim_B, 3xJan)<br>\n']
        l_htm += ['min/max = min/max of calculated pKa values<br>\n']
        l_htm += ['pka exp = experimental pKa value<br>\n']
        s = '('
        for res_name,pka_model in d_null.items():
            s += '%s %s, ' %(res_name, pka_model)
        s = s[:-2]+')'
        l_htm += ['pka shift = shift of experimental pKa value from model pKa value %s<br>\n' %(s)]
        l_htm += ['min/max diff = smallest/largest difference between experimental pKa value and calculated pKa values<br>\n']
        l_htm += ['rmsd = RMSD of differences between experimental pKa value and calculated pKa values (excl. Gernot, Emil_2, Jim_B, 3xJan)<br>\n']
        l_htm += ['For each group/method/PI the submitted pKa value is in parenthesis and the difference from the experimental pKa value is in front outside the parenthesis - "pka_diff (pka_calc)"<br>\n']

##        ## fixed table headers
##        l_fixed = []
##        l_fixed += ['<table class="fixed">\n']
##        for col in l_columns:
##            l_fixed += ['<col style= "width: 5em">\n']
##        s = ''
##        for col in l_columns:
##            s += '<td>%s</td>' %(col)
##        s += '\n'
##        l_fixed += [s]
##        l_fixed += ['</table>\n\n']
##
##        ## append fixed table headers
##        l_htm += l_fixed

        ## initiate table
##        l_htm += ['<div id="scroll">\n']
        l_htm += ['<table border="0" class="sortable" id="sortable"><br>\n']

##        s = ''
##        for col in l_columns:
##            s += '<td>%s</td>' %(col)
##        s += '\n'
##        l_htm += [s]

        ## table headers
        l_thead = ['\n<thead>\n']
        l_thead += ['<tr>\n']
        for s in l_columns:
            l_thead += ['<td><b>%s</b></td>\n' %(s,)]
        for PI in l_PIs_htm:
            l_thead += ['<td><font color="%s">%s</font></td>\n' %(d_PIs[PI]['color'],d_PIs[PI]['short'],)]
        l_thead += ['<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>\n']
        l_thead += ['</tr>\n']
        l_thead += ['</thead>\n']


        l_tbody = ['\n<tbody>\n']

        for protein in l_proteins:

            ## loop over residues
            l_res_IDs = list(d_residues[protein])
            l_res_IDs.sort()
            for res_ID in l_res_IDs:

                ## skip residue?
                bool_skip = self.skip_residue(protein,res_ID,)
                if bool_skip == True:
                    continue

                res_no  = int(res_ID[-4:])
                res_name = res_ID[:3]

                ## experimental?
                experimental = False
                if protein in d_nonmutants.keys():
                    protein_type = d_nonmutants[protein]
                else:
                    protein_type = 'mutants'
                if res_no in d_experimental[protein_type].keys():
                    if res_name in d_experimental[protein_type][res_no].keys():
                        experimental = True
                        pka_exp = d_experimental[protein_type][res_no][res_name]
                ## skip if no experimental
                if limited == True and experimental == False:
                    continue

                ## restricted?
                residue_restricted = False
                if restricted == True and protein_type == 'mutants':
                    if '%s%i' %(res_name,res_no,) not in l_allowed:
                        residue_restricted = True

                ## protein
                if protein in d_mutations.keys(): ## non-delta-phs mutant
                    res_mut = d_mutations[protein]
                    if protein in d_parents.keys(): ## not delta+phs
                        protein_htm = '%1s%i%s %s' %(d_res3res1[l_sequence[res_no-1]],res_no,d_res3res1[res_mut[:3]],d_parents[protein].upper(),)
                    else: ## delta+phs
                        protein_htm = '%1s%i%s &Delta;+PHS' %(d_res3res1[l_sequence[res_no-1]],res_no,d_res3res1[res_mut[:3]],)
                else:
                    if protein in d_parents.keys(): ## non-delta-phs
                        protein_htm = d_parents[protein].upper()
                    else: ## delta-phs
                        if protein == '3bdc': ## non-mutant
                            protein_htm = '&Delta;+PHS'
                        else: ## mutant
                            protein_htm = '%s &Delta;+PHS' %(protein[:protein.index('_')].upper())
                ## pdb
                if len(protein) == 4:
                    pdb_htm = protein.lower()
                else:
                    pdb_htm = ''
##                protein_htm = protein_htm.replace(' &Delta;+PHS','')

                ## initiate htm row
                l_tbody += ['<tr>\n']

                ##
                ## 1) columns protein, res_name, res_no
                ##
                if experimental == True:
                    l_tbody += ['<td><b>%s</b></td>\n' %(protein_htm,)]
                    l_tbody += ['<td><b>%s</b></td>\n' %(pdb_htm,)]
                    l_tbody += ['<td><b>%3s</b></td>\n' %(res_name,)]
                    l_tbody += ['<td><b>%03i</b></td>\n' %(res_no,)]
                else:
                    l_tbody += ['<td>%s</td>\n' %(protein_htm,)]
                    l_tbody += ['<td>%s</td>\n' %(pdb_htm,)]
                    l_tbody += ['<td>%3s</td>\n' %(res_name,)]
                    l_tbody += ['<td>%03i</td>\n' %(res_no,)]

                ##
                ## 2a) calculate diff per residue (should be in a diff func!!!)
                ##
                l = []
                l_diff = []
                d_diff = {}
                for PI in l_PIs:
                    if not protein in d_pkas[PI].keys():
                        d_diff[PI] = 'N/A'
                        continue
                    if not res_ID in d_pkas[PI][protein].keys():
                        d_diff[PI] = 'N/A'
                        continue
                    s_pka = d_pkas[PI][protein][res_ID]['s_pka']
                    if s_pka == 'N/A' or ';' in s_pka:
                        d_diff[PI] = 'N/A'
                        continue
                    ##
                    ## stderr
                    ##
                    if PI not in [## temp!!!
                        'Gernot Kieseritzky', ## way off
                        'Emil Alexov_2', ## duplicates
                        'Jim_Warwicker_B', ## duplicates
                        'JanJensen_Jdec19', 'JanJensen_Wdec19', 'JanJensen_Jjan15', ## duplicates
                        ]: 
                        l += [float(s_pka)]
                    if (
                        experimental == True
                        and
                        ## not a limit value
                        type(pka_exp).__name__ != 'str'
                        ):
                        pka_calc = float(d_pkas[PI][protein][res_ID]['pka1'])
##                        pka_diff = '%+04.1f' %(pka_exp-pka_calc) ## diff w sign
                        pka_diff = '%04.1f' %(abs(pka_exp-pka_calc)) ## abs diff
                        if d_differences[PI][protein][res_ID] != pka_calc-pka_exp: ## temp!!! if no errors then no need to recalculate pkadiff in this function!
                            print d_differences[PI][protein][res_ID]
                            print pka_calc-pka_exp
                            print PI, protein, res_ID
                            print pka_exp, d_pkas[PI][protein][res_ID]
                            stop
                        if PI not in [## temp!!!
                            'Gernot Kieseritzky', ## way off
                            'Emil Alexov_2', ## duplicates
                            'Jim_Warwicker_B', ## duplicates
                            'JanJensen_Jdec19', 'JanJensen_Wdec19', 'JanJensen_Jjan15', #'Wjan15', ## duplicates
                            ]: 
                            l_diff += [float(s_pka)-pka_exp]
                    else:
                        pka_diff = 'N/A'
                    d_diff[PI] = pka_diff

                ##
                ## 2b) calculate rmsd per residue (should be in a diff func!!!)
                ##
                if len(l) > 1:
                    average, stderr = statistics.stderr(l) ## stddev and not stderr
                    n = len(l)
                    average = '%04.1f' %(average)
                    stderr = '%.1f' %(stderr)
                    minimum = '%.1f' %(min(l))
                    maximum = '%.1f' %(max(l))
                else:
                    n = average = stderr = minimum = maximum = 'N/A'
                if len(l_diff) > 1:
                    rmsd_diff = statistics.rmsd(l_diff)
                    max_diff = 0
                    min_diff = 99
                    for diff in l_diff:
                        if abs(diff) > max_diff:
                            max_diff = abs(diff)
                        if abs(diff) < min_diff:
                            min_diff = abs(diff)
                else:
                    rmsd_diff = '&nbsp;'
                    max_diff = '&nbsp;'
                    min_diff = '&nbsp;'

                ##
                ## 2c) columns statistics
                ##
                ## n
                l_tbody += ['<td>%02i</td>\n' %(n)]
                ## average
                l_tbody += ['<td>%s</td>\n' %(average,)]
                ## stderr
                if stderr == 'N/A':
                    l_tbody += ['<td>%s</td>\n' %(stderr)]
                elif float(stderr) > 1:
                    l_tbody += ['<td><font color="red">%s</font></td>\n' %(stderr)]
                else:
                    l_tbody += ['<td><font color="green">%s</font></td>\n' %(stderr)]
                ## max, min
                l_tbody += ['<td>%s</td>\n' %(minimum)]
                l_tbody += ['<td>%s</td>\n' %(maximum)]
                ## rmsd (clean up all the if statements!!!)
                if experimental == True and residue_restricted == False:
                    ## limit value
                    if type(pka_exp).__name__ == 'str':
                        ## pka_exp
                        l_tbody += ['<td><b><a href="%s/%s/%s_%s.htm">%s</a></b></td>\n' %(path_htm_titcurve,folder,protein,res_ID,pka_exp,)]
                        ## pka_shift
                        l_tbody += ['<td>%1s%0.1f</td>\n' %(pka_exp[0], float(pka_exp[1:])-d_null[res_name])]
                    ## absolute value
                    else:
                        ## pka_exp
                        l_tbody += ['<td><b><a href="%s/%s/%s_%s.htm">%04.1f</a></b></td>\n' %(path_htm_titcurve,folder,protein,res_ID,pka_exp,)]
                        ## pka_shift
                        l_tbody += ['<td>%0.1f</td>\n' %(pka_exp-d_null[res_name])]
                    ## max, min diff
                    l_tbody += ['<td>%s</td>\n' %(min_diff)]
                    l_tbody += ['<td>%s</td>\n' %(max_diff)]
                elif experimental == False and residue_restricted == False:
                    ## pka_exp
                    l_tbody += ['<td><a href="%s/%s/%s_%s.htm">%s</a></td>\n' %(path_htm_titcurve,folder,protein,res_ID,'N/A',)]
                    ## pka_shift
                    l_tbody += ['<td>%s</td>\n' %('N/A')]
                    ## max, min diff
                    l_tbody += ['<td>%s</td>\n' %(min_diff)]
                    l_tbody += ['<td>%s</td>\n' %(max_diff)]
                elif experimental == True and residue_restricted == True:
                    ## pka_exp
                    l_tbody += ['<td><a href="%s/combined_restricted/%s_%s.htm">%s</a></td>\n' %(path_htm_titcurve,protein,res_ID,'restricted',)]
                    ## pka_shift
                    l_tbody += ['<td>%s</td>\n' %('N/A')]
                    ## max, min diff
                    l_tbody += ['<td>%s</td>\n' %(min_diff)]
                    l_tbody += ['<td>%s</td>\n' %(max_diff)]

                ## rmsd
                if rmsd_diff in ['N/A','&nbsp;',]:
                    l_tbody += ['<td>%s</td>\n' %(rmsd_diff,)]
                elif float(rmsd_diff) > 1:
                    l_tbody += ['<td><font color="red">%.1f</font></td>\n' %(rmsd_diff,)]
                else:
                    l_tbody += ['<td><font color="green">%.1f</font></td>\n' %(rmsd_diff,)]

                ##
                ## 3) columns pka, pkadiff per PI
                ##
                for PI in l_PIs_htm:
                    ## no calculated values for protein
                    if not protein in d_pkas[PI].keys():
                        l_tbody += ['<td>&nbsp;</td>\n']
                    ## no calculated value for residue
                    elif res_ID not in d_pkas[PI][protein].keys():
                        l_tbody += ['<td>&nbsp;</td>\n']
                    ## no experimental value (not restricted)
                    elif residue_restricted == False and d_diff[PI] == 'N/A':
                        l_tbody += [
                            '<td><a href="%s/%s/%s_%s.htm">(%s)</a></td>\n' %(
                                path_htm_titcurve,PI,protein,res_ID,
                                d_pkas[PI][protein][res_ID]['s_pka'],
                                )
                            ]
                    ## experimental value (not restricted)
                    elif residue_restricted == False:
                        l_tbody += [
                            '<td><a href="%s/%s/%s_%s.htm">%s (%s)</a></td>\n' %(
                                path_htm_titcurve,PI,protein,res_ID,
                                d_diff[PI],
                                d_pkas[PI][protein][res_ID]['s_pka'],
                                )
                            ]
                    ## experimental value (restricted)
                    elif residue_restricted == True:
                        l_tbody += [
                            '<td>(%s)</td>\n' %(
                                d_pkas[PI][protein][res_ID]['s_pka'],
                                )
                            ]
                    else:
                        stop_not_expected

                l_tbody += ['<td>&nbsp;</td>\n']
                l_tbody += ['</tr>\n']
        l_tbody += ['</tbody>\n']

        l_htm += l_thead+l_tbody
        l_htm += ['</table>\n']
        l_htm += ['</div>\n']

##        ## append fixed table headers
##        l_htm += l_fixed

        ##
        ## terminate htm
        ##
        l_htm += ['</body>\n</html>\n']
        l_htm = l_header+l_htm

        ##
        ## write htm
        ##
        fd = open('%s.htm' %(prefix,),'w')
        fd.writelines(l_htm)
        fd.close()

        return


    def skip_residue(self,protein,res_ID,):

        ## 1) terminal residue
        if res_ID[:5] in ['NTERM','CTERM',]:
            return True
        ## 2) mutant and non-mutant proteins
        if protein in d_bgm_residues.keys():
            ## residue of interest?
            if not (
                res_ID[:3] in d_bgm_residues[protein] ## res_name
                or
                res_ID in d_bgm_residues[protein]
                ):
                return True
        ## 3) mutant protein but not mutated residue
        elif len(protein) > 4:
            res1 = protein[-11]
            res_name = d_res1res3[res1]
            res_no = int(protein[1:-11])
            chain = 'A'
            if res_ID != '%s:%1s:%04i' %(res_name,chain,res_no,):
                return True

        return False


    def plot_histograms(self,d_pkas,d_rmsd,d_differences,):

        bw = 1.
        set_keys = set([
            'all','ASP','GLU','HIS','LYS','ARG','wt','delta+phs','mutants',
            'mutant_pdb','mutant_model','all_pdb',
            ])
        print

        d_count = {}
        d_maps = {}
        for PI in l_PIs:

            print 'plotting histograms', PI
            d_count[PI] = {}
            d_maps[PI] = {}
            for k in set_keys:
                d_count[PI][k] = {}
                for box in range(-15,15+1):
                    d_count[PI][k][box] = []


            ## 1) count occurences
            for protein in l_proteins:

                if not protein in d_pkas[PI].keys():
                    continue

                for res_ID in d_pkas[PI][protein].keys():

                    ## skip residue?
                    bool_skip = self.skip_residue(protein,res_ID,)
                    if bool_skip == True:
                        continue
                    
                    res_no  = int(res_ID[-4:])
                    res_name = res_ID[:3]

                    if protein in d_nonmutants.keys():
                        protein_type = d_nonmutants[protein]
                    else:
                        protein_type = 'mutants'

                    ## skip if no experimental data
                    if not res_no in d_experimental[protein_type].keys():
                        continue
                    ## skip if no experimental data
                    if not res_name in d_experimental[protein_type][res_no].keys():
                        continue
                    pka_exp = d_experimental[protein_type][res_no][res_name]
                    ## skip if limit value
                    if type(pka_exp).__name__ == 'str':
                        continue
                    ## skip if non-mutated residue in mutant protein
                    if protein_type == 'mutants' and res_name == l_sequence[res_no-1]:
                        continue
                    ## no titration in pH range
                    if d_pkas[PI][protein][res_ID]['s_pka'] == 'N/A':
                        continue

                    diff = d_differences[PI][protein][res_ID]
                    diff += .5*bw
                    box = (diff-diff%bw)/bw
                    d_count[PI]['all'][box] += [[PI,protein,res_ID,]]
                    d_count[PI][protein_type][box] += [[PI,protein,res_ID,]]
                    d_count[PI][res_name][box] += [[PI,protein,res_ID,]]
                    if protein_type == 'mutants':
                        if len(protein) > 4:
                            d_count[PI]['mutant_model'][box] += [[PI,protein,res_ID,]]
                        else:
                            d_count[PI]['mutant_pdb'][box] += [[PI,protein,res_ID,]]
                            d_count[PI]['all_pdb'][box] += [[PI,protein,res_ID,]]
                    else:
                        if len(protein) != 4:
                            stop
                        d_count[PI]['all_pdb'][box] += [[PI,protein,res_ID,]]


            ## 2) plot histogram
            for k in d_count[PI].keys():

                if os.path.isfile('png_histograms/%s_%s.png' %(k,PI.replace(' ',''),)):
                    continue

                l = []
                i = 0
                for box in range(-15,15+1):
                    l += ['%s %s\n' %(box,len(d_count[PI][k][box]),)]
                    i += len(d_count[PI][k][box])

                if i < 2:
                    continue

                fd = open('%s%s.dat' %(k,PI.replace(' ',''),),'w')
                fd.writelines(l)
                fd.close()

                lines_gnuplot = []
                lines_gnuplot += [
                    'set terminal postscript eps enhanced color "Helvetica" 24\n',
                    'set output "gnuplot.ps"\n',
                    'set size square\n',
                    'set size 2,2\n', ## scale size
                    'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
                    'set xlabel "calculated pKa - experimental pKa"\n',
                    'set ylabel "count"\n',
                    'set title "%s (%s, RMSD = %s, n = %i)"\n' %(PI.replace('_','\_'), k, d_rmsd[PI]['pka1']['rmsd'], i,),
                    'set xtics 1\n',
##                    'set style data histograms\n',
                    'set style fill solid border -1\n',
                    'unset key\n',
##                    'bw = 1\n', ## histogram box width
##                    'bin(x,width) = width*floor(x/width)\n',
                    ]

##                line = 'plot [-10:10][0:20] "%s%s.dat" u (bin($1,bw)):(1.0) smooth freq w boxes\n' %(s,PI.replace(' ',''),)
                line = 'plot [-15:15][0:25] "%s%s.dat" w boxes lc rgb "#%6s"\n' %(k,PI.replace(' ',''),d_PIs[PI]['color'],)
                lines_gnuplot += [line]

                fd = open('gnuplot.settings','w')
                fd.writelines(lines_gnuplot)
                fd.close()

                os.system('gnuplot gnuplot.settings')
                os.system('convert gnuplot.ps png_histograms/%s_%s.png' %(k,PI.replace(' ',''),))
                os.remove('gnuplot.ps')
                os.remove('%s%s.dat' %(k,PI.replace(' ',''),))
                os.remove('gnuplot.settings')

            ##
            ## 3) htm maps per histogram box (histogram htm)
            ##
            for subfolder_htm_titcurves,folder_htm_histograms in [['histograms','htm_histograms',],['histograms_restricted','htm_histograms_restricted',],]:

                for k in d_count[PI].keys():

                    ##
                    ## map
                    ##
                    lines_htm_map = ['<map id="%s_%s" name="%s_%s">\n' %(PI,k,PI,k,)]
                    for box in range(-15,15+1):
                        count = len(d_count[PI][k][box])
                        if count == 0:
                            continue
                        x_min = 161
                        x_max = 592
                        y_min = 36
                        y_max = 469
                        lines_htm_map += [
                            '<area shape="rect" coords="%i,%i,%i,%i" href="../htm_titcurves/%s/%s_%s_%i.htm">\n' %(
                                x_min+(box+15-.5)*(x_max-x_min)/(30.),
                                y_max-count*(y_max-y_min)/25.,
                                x_min+(box+15-.5+1)*(x_max-x_min)/(30.),
                                y_max,
                                subfolder_htm_titcurves,
                                PI,k,box,
                                ),
                            ]
                    lines_htm_map += ['</map>\n',]
                    d_maps[PI][k] = lines_htm_map

                    ##
                    ## img
                    ##
                    lines_htm = []
                    lines_htm += d_maps[PI][k]
                    lines_htm += [
                        '<img src="../png_histograms/%s_%s.png" usemap="#%s_%s" border="0">\n<br>\n' %(
                            k,PI.replace(' ',''),
                            PI,k,
                            )
                        ]
                    fd = open('%s/%s_%s.htm' %(folder_htm_histograms,k,PI,),'w')
                    fd.writelines(lines_htm)
                    fd.close()

            for folder in ['histograms','histograms_restricted',]:

                ##
                ## 4) htm map targets (titcurves)
                ##
                for k in d_count[PI].keys():
                    for box in range(-15,15+1):
                        lines_htm = []
                        if folder == 'histograms_restricted':
                            lines_htm += ['Some titration curves might not be shown, if they are not part of the restricted data set<br>\n']
                        if len(d_count[PI][k][box]) == 0:
                            continue
                        for l in d_count[PI][k][box]:
                            PI = l[0]
                            protein = l[1]
                            res_ID = l[2]
                            res_no  = int(res_ID[-4:])
                            res_name = res_ID[:3]
                            if folder == 'histograms_restricted':
                                if not '%s%i' %(res_name,res_no,) in l_allowed:
                                    continue
                            lines_htm += [
                                '<img src="../../png_titcurves/%s/%s_%s.png"></img>\n<br>\n' %(
                                    PI,protein,res_ID,
                                    )
                                ]
                        fd = open('htm_titcurves/%s/%s_%s_%s.htm' %(folder,PI,k,box,),'w')
                        fd.writelines(lines_htm)
                        fd.close()

            ## 5) htm per PI
            lines_htm = []
            for k in set_keys:

                if not os.path.isfile('png_histograms/%s_%s.png' %(k,PI.replace(' ',''),)):
                    continue

                lines_htm += d_maps[PI][k]
                lines_htm += [
                    '<img src="../png_histograms/%s_%s.png" usemap="#%s_%s" border="0">\n<br>\n' %(
                        k,PI.replace(' ',''),
                        PI,k,
                        )
                    ]

                fd = open('htm_PI/%s.htm' %(PI.replace(' ',''),),'w')
                fd.writelines(lines_htm)
                fd.close()

        ## 6) htm per type
        for folder in ['htm_histograms','htm_histograms_restricted',]:
            for k in set_keys:
                lines_htm = []
                for PI in l_PIs:
                    if not os.path.isfile('png_histograms/%s_%s.png' %(k,PI.replace(' ',''),)):
                        continue
                    lines_htm += d_maps[PI][k]
                    lines_htm += [
                        '<img src="../png_histograms/%s_%s.png" usemap="#%s_%s" border="0">\n<br>\n' %(
                            k,PI.replace(' ',''),
                            PI,k,
                            )
                        ]
                fd = open('%s/%s.htm' %(folder,k,),'w')
                fd.writelines(lines_htm)
                fd.close()
                
        return                


    def parse_previous_fits(self):

        d_pkas = {}

        ##
        ## parse previously fitted pka values
        ##
        if os.path.isfile('fits.csv'):
            fd = open('fits.csv','r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                l = line.strip().split(',')
                PI = l[0]
                protein = l[1]
                res_ID = l[2]
                pka1 = l[3]
                s_pka = l[4].strip()
                if not PI in d_pkas.keys():
                    d_pkas[PI] = {}
                if not protein in d_pkas[PI].keys():
                    d_pkas[PI][protein] = {}
                d_pkas[PI][protein][res_ID] = {'s_pka':s_pka,'pka1':pka1,}

        return d_pkas


    def plot_titcurves_combined(self,d_residues,d_titcurves,d_pkas,d_differences,d_slopes,):

        for folder in ['combined','combined_restricted',]:

            if not os.path.isdir('png_titcurves/%s' %(folder)):
                os.mkdir('png_titcurves/%s' %(folder))

            ##
            ## plots per residue
            ##
            for protein in l_proteins:
                for res_ID in d_residues[protein]:

                    res_no  = int(res_ID[-4:])
                    res_name = res_ID[:3]

                    ## protein type
                    if protein in d_nonmutants.keys():
                        protein_type = d_nonmutants[protein]
                    else:
                        protein_type = 'mutants'

                    ## restricted residue?
                    residue_restricted = True
                    if (
                        protein_type != 'mutants'
                        or
                        (folder == 'combined')
                        or
                        (folder == 'combined_restricted' and '%s%i' %(res_name,res_no,) in l_allowed)
                        ):
                        residue_restricted = False
                    
                    ## experimental?
                    pka_exp = None
                    experimental = False
                    limit = False
                    if res_no in d_experimental[protein_type].keys():
                        if res_name in d_experimental[protein_type][res_no].keys():
                            experimental = True
                            pka_exp = d_experimental[protein_type][res_no][res_name]
                            ## limit value?
                            if type(pka_exp).__name__ == 'str':
                                pka_exp = float(d_experimental[protein_type][res_no][res_name][1:])
                                limit = True

                    ## skip residue?
                    bool_skip = self.skip_residue(protein,res_ID,)
                    if bool_skip == True:
                        continue

                    ## previously plotted?
                    if os.path.isfile('png_titcurves/%s/%s_%s.png' %(folder,protein,res_ID,)):
                        continue

                    print 'plotting %s' %(folder), protein, res_ID

                    l_plots = []
                    s = 'plot [-10:20][-1:1]'

                    ##
                    ## titcurv submitted
                    ##
                    for PI in l_PIs:
                        
                        if d_PIs[PI]['submission'] != 'titcurv':
                            continue
                        if not protein in d_titcurves[PI].keys():
                            continue
                        if not res_ID in d_titcurves[PI][protein].keys():
                            continue
                        l_plots += [PI]
                        l = []
                        for pH,charge in d_titcurves[PI][protein][res_ID].items():
                            l += ['%s %s\n' %(pH, charge,)]
                        fd = open('%s.tit' %(PI), 'w')
                        fd.writelines(l)
                        fd.close()

                    for PI in l_plots:
                        PI_postscript = PI.replace('_','\_')
                        s += '"%s.tit" u 1:2 pt %i ps 2 lc rgb "#%6s" t "%s", ' %(PI,d_PIs[PI]['pt'],d_PIs[PI]['color'],PI_postscript,)

                    ##
                    ## no titcurv submitted
                    ##
                    s_function = 'f'
                    l_functions = []
                    for PI in l_PIs:
                        if not protein in d_pkas[PI].keys():
                            continue
                        if not res_ID in d_pkas[PI][protein].keys():
                            continue
                        if d_PIs[PI]['submission'] == 'titcurv':
                            continue
                        if res_ID[:3] in ['LYS','ARG','HIS','NTE',]:
                            a = 1.
                            b = -1.
                        elif res_ID[:3] in ['ASP','GLU','TYR','CTE',]:
                            a = 0.
                            b = -1.
                        else:
                            print res_ID
                            stop
                        c = float(d_pkas[PI][protein][res_ID]['pka1'])
                        l_functions += ['%s(x) = %s+%s/(1+10**(%s-x))\n' %(s_function,a,b,c,)]
                        PI_postscript = PI.replace('_','\_')
                        s += '%s(x) lt 1 lw 8 lc rgb "#%6s" t "%s", ' %(s_function,d_PIs[PI]['color'],PI_postscript,)
                        s_function = s_alphabet[s_alphabet.index(s_function)+1]
                        l_plots += [PI]

                    if len(l_plots) == 0:
                        print 'nothing to plot for', protein, res_ID
                        continue

                    ##
                    ## gnuplot source
                    ##
                    protein_postscript = protein.replace('_','\_').upper().replace('DELTA','{/Symbol D}')
                    title = "'%s %s'" %(protein_postscript, res_ID) ## N.B. postscript has to be in single brackets!!!
                    lines = [
                        'set terminal postscript eps enhanced color "Helvetica" 30\n',
                        'set output "%s_%s.ps"\n' %(protein,res_ID,),
                        'set size 3,3\n',
                        'set grid\n',
                        'set mxtics 1\n',
                        'set xtics 1\n',
                        'set xlabel "pH"\n',
                        'set ylabel "charge"\n',
                        'set title %s\n' %(title),

                        ## black background
                        'set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back\n',
                        'set object 1 rect fc rgb "black" fillstyle solid 1.0\n',
                        
    ####                    'set key box out vert center top\n',
                        'unset key\n',
                        'set grid front lc rgb "#FFFFFF"\n', ## white grid
                        ]
                    if (
                        experimental == True
                        and
                        residue_restricted == False
                        ):
                        if limit == False:
                            lines += [
                                'set arrow from %f, graph 0 to %f, graph 1 nohead lt 1 lw 10 lc rgb "#FFFFFF"\n' %(pka_exp,pka_exp,),
                                ]
                        else:
                            lines += [
                                'set arrow from %f, graph 0 to %f, graph 1 nohead lt 3 lw 10 lc rgb "#808080"\n' %(pka_exp,pka_exp,),
                                ]
                    s = s[:-2]+'\n' ## finish plot line
                    lines += l_functions ## add line of functions
                    lines += [s] ## add plot line

                    fd = open('gnuplot_%s_%s.src' %(protein,res_ID,),'w')
                    fd.writelines(lines)
                    fd.close()

                    os.system('gnuplot gnuplot_%s_%s.src' %(protein,res_ID,))
                    os.system('convert %s_%s.ps %s_%s.png' %(protein,res_ID,protein,res_ID,))
                    os.system('mv %s_%s.png png_titcurves/%s/.' %(protein,res_ID,folder,))

    ##                print PI, protein, res_ID
    ##                if experimental == True:
    ##                    stop
                    os.remove('gnuplot_%s_%s.src' %(protein,res_ID,))
                    os.remove('%s_%s.ps' %(protein,res_ID,))

                    for PI in l_plots:
                        if os.path.isfile('%s.tit' %(PI)):
                            os.remove('%s.tit' %(PI))

                    ##
                    ## htm
                    ##
                    l_htm = ['<html>\n<head>\n']
                    l_htm += [' <title>%s %s</title>\n' %(protein, res_ID,)]
                    l_htm += [' <script src="../../jmol/Jmol.js" type="text/javascript"></script>\n']
                    l_htm += ['</head>\n']
                    l_htm += ['<body bgcolor="white">\n']
                    if experimental == True and residue_restricted == False:
                        if limit == False:
                            l_htm += ['<font color="black">The experimental pKa value (%s) is shown with a white vertical line on the plot.</font>' %(pka_exp)]
                        else:
                            l_htm += ['<font color="black">The limit of the experimental pKa value (%s) is shown with a white stippled vertical line on the plot</font>' %(pka_exp)]
                    elif experimental == False:
                        l_htm += ['<font color="black">No experimental results.</font>']
                    else:
                        l_htm += ['<font color="black">Experimental results not released yet.</font>']
                    l_htm += ['<br><br><a href="../../pdb/%s.pdb"><font color="black">Download PDB file</font></a><br>\n' %(d_pdbs[protein],),]
                    ## initiate table
                    l_htm += ['<table>\n<tr>\n']
                    l_htm += ['<td><img src="../../png_titcurves/%s/%s_%s.png"></td>\n' %(folder,protein,res_ID,)]
                    l_htm += ['<td valign="top">']
                    
                    for PI in l_plots:
                        s = '('
                        ## slope
                        if d_PIs[PI]['submission'] == 'titcurv':
                            if protein in d_slopes[PI].keys():
                                if res_ID in d_slopes[PI][protein].keys():
                                    s += 'slope = %.1f, ' %(d_slopes[PI][protein][res_ID],)
                        ## pka_diff
                        if (
                            experimental == True and limit == False and d_pkas[PI][protein][res_ID]['s_pka'] != 'N/A'
                            and
                            residue_restricted == False
                            ):
                            s += 'pka_diff = %.1f, ' %(d_differences[PI][protein][res_ID],)
                        ## pka_calc
                        s += 'pka_calc = %s)' %(d_pkas[PI][protein][res_ID]['s_pka'],)
                        l_htm += ['<font color="%s"><ul type="%s"><li><b>%s %s</b></li></ul></font>\n' %(
                            d_PIs[PI]['color'], d_PIs[PI]['ul'], d_PIs[PI]['short'], s,
                            )]

                    l_htm += ['</td>\n']
                    ## terminate table
                    l_htm += ['</tr>\n</table>\n']

                    ## jmol
                    l_htm += ['<form>\n']
                    l_htm += [' <script type="text/javascript">\n']
                    ## script initialization
                    l_htm += ['  jmolInitialize("../../jmol");\n']
                    ## jmol script location
                    l_htm += ['  jmolApplet(768, "script ../../jmol_scripts/%s_%s.jmol");\n' %(protein,res_ID,)]
                    l_htm += [' </script>\n']
                    l_htm += ['</form>\n']

                    ## terminate htm
                    l_htm += ['</body>\n</html>\n']
                    if not os.path.isdir('htm_titcurves/%s' %(folder)):
                        os.mkdir('htm_titcurves/%s' %(folder))
                    ## write htm
                    fd = open('htm_titcurves/%s/%s_%s.htm' %(folder,protein,res_ID,),'w')
                    fd.writelines(l_htm)
                    fd.close()

        return


    def fit_and_plot_titcurves_individual(self,d_titcurves,d_pkas,d_residues,):

        print

        for PI in l_PIs:

            print 'plotting titcurve', PI

            if not os.path.isdir('/local/tc/pkacoop/htm_titcurves/%s' %(PI)):
                os.mkdir('/local/tc/pkacoop/htm_titcurves/%s' %(PI))
            
            for protein in l_proteins:

                if d_PIs[PI]['submission'] == 'titcurv':
                    if not protein in d_titcurves[PI].keys():
                        continue
                elif d_PIs[PI]['submission'] == 'pka':
                    if not protein in d_pkas[PI].keys():
                        continue

                for res_ID in d_residues[protein]:

                    ## skip residue?
                    bool_skip = self.skip_residue(protein,res_ID,)
                    if bool_skip == True:
                        continue

                    if d_PIs[PI]['submission'] == 'titcurv':
                        if not res_ID in d_titcurves[PI][protein].keys():
                            continue
                    elif d_PIs[PI]['submission'] == 'pka':
                        if not res_ID in d_pkas[PI][protein].keys():
                            continue

                    res_no  = int(res_ID[-4:])
                    res_name = res_ID[:3]

                    if protein in d_nonmutants.keys():
                        protein_type = d_nonmutants[protein]
                    else:
                        protein_type = 'mutants'

                    pka_exp = None
                    experimental = False
                    if res_no in d_experimental[protein_type].keys():
                        if res_name in d_experimental[protein_type][res_no].keys():
                            pka_exp = d_experimental[protein_type][res_no][res_name]
                            experimental = True
                            ## limit value
                            if type(pka_exp).__name__ == 'str':
                                pka_exp = float(d_experimental[protein_type][res_no][res_name][1:])


                    ##
                    ## plot titration curves
                    ##
                    if d_PIs[PI]['submission'] == 'titcurv':

                        if protein in d_pkas[PI].keys():
                            if (
                                ## already plotted?
                                os.path.isfile('png_titcurves/%s/%s_%s.png' %(PI,protein,res_ID,))
                                and
                                ## already calculated?
                                res_ID in d_pkas[PI][protein].keys()
                                ):
                                continue

                        ## write gnuplot data file
                        l = []
                        for pH,charge in d_titcurves[PI][protein][res_ID].items():
                            l += ['%s %s\n' %(pH, charge,)]
                        fd = open('fit.tit','w')
                        fd.writelines(l)
                        fd.close()

                        print 'plotting', PI, protein, res_ID

                        ## data suitable for fitting?
                        bool_fit, full_titration = self.checkdata(d_titcurves[PI][protein][res_ID], PI, protein, res_ID, pka_exp, bool_plot=True,)
                        if bool_fit == False:

                            pka1 = 'N/A'
                            s_pka = 'N/A'
                            
                            ## add to dictionary
                            if not protein in d_pkas[PI].keys():
                                d_pkas[PI][protein] = {}
                            d_pkas[PI][protein][res_ID] = {'s_pka':s_pka,'pka1':pka1,}

                            ## add to csv file
                            fd = open('fits.csv','a')
                            fd.write('%s,%s,%s,%s,%s\n' %(PI,protein,res_ID,pka1,s_pka,))
                            fd.close()

                            continue

                        if full_titration == False or full_titration == True:

                            if res_ID[:3] in ['LYS','ARG','HIS',]:
                                a = 1
                            elif res_ID[:3] in ['ASP','GLU','TYR',]:
                                a = 0
                            else:
                                print res_ID
                                stop
                            b = -1
                            d_fit = {
                                0:{
                                    'function':'%s+%s/(1+10**(c-x))\n' %(a,b,),
                                    'variables':'c'
                                    }
                                }
                            d_fit[0] = self.gnuplot_fit(d_fit[0],d_titcurves[PI][protein][res_ID],)
                            d_fit[0]['a'] = a
                            d_fit[0]['b'] = b
                            if a not in [0,1,] or b != -1:
                                stop

                            s_pka = '%.1f' %(d_fit[0]['c'])
                            pka1 = '%.1f' %(d_fit[0]['c'])

                            s_ss = 'N/A'
                            P = 'N/A'
                            MODEL = 0
                            title = "'%s   %s %s p=%s pK_as=%i (%s) SS %s'" %(
                                PI,
                                protein.replace('_','\_').upper().replace('DELTA','{/Symbol D}'),
                                res_ID, P, MODEL+1, s_pka, s_ss,
                                ),
                            self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)

##                                ## temp!!! removal of previos calc value
##                                fd = open('fits.csv','r')
##                                lines = fd.readlines()
##                                fd.close()
##                                lines2 = []
##                                for line in lines:
##                                    l = line.strip().split(',')
##                                    if PI == l[0] and protein == l[1] and res_ID == l[2]:
##                                        continue
##                                    lines2 += line
##                                fd = open('fits.csv','w')
##                                fd.writelines(lines2)
##                                fd.close()

                            fd = open('fits.csv','a')
                            fd.write('%s,%s,%s,%s,%s\n' %(PI,protein,res_ID,pka1,s_pka,))
                            fd.close()

                            if not protein in d_pkas[PI].keys():
                                d_pkas[PI][protein] = {}
                            d_pkas[PI][protein][res_ID] = {'s_pka':s_pka,'pka1':pka1,}

                            os.remove('fit.tit')

                        else:

##                            continue ## temp!!!

                            ##
                            ## F-test
                            ##
                            MODEL, P, d_fit = self.get_model(d_titcurves[PI][protein][res_ID],)
                            if P != 'N/A':
                                P = '%.4f' %(P)

                            print PI, protein, res_ID

                            ##
                            ## plot non-HH curves
                            ##
                            if MODEL > 0:

                                s_pka = ''
                                s_ss = ''
                                for i in range(MODEL+1):
                                    s_pka += '%.1f; ' %(d_fit[MODEL][['c','e','g',][i]],)
                                    s_ss += '%.4f; ' %(d_fit[i]['SS'],)
                                s_pka = s_pka[:-2]
                                s_ss = s_ss[:-2]

                                title = "'%s   %s %s p=%s pK_as=%i (%s) SS %s'" %(
                                    PI, protein.replace('_','\_').upper().replace('DELTA','{/Symbol D}'),res_ID,P,MODEL+1,
                                    s_pka,
                                    s_ss,
                                    ),
                                print title

                                self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)

                            else:

                                s_pka = '%.1f' %(d_fit[MODEL]['c'])
                                s_ss = '%.4f' %(d_fit[0]['SS'])
                                title = "'%s   %s %s p=%s pK_as=%i (%s) SS %s'" %(
                                    PI, protein.replace('_','\_').upper().replace('DELTA','{/Symbol D}'),res_ID,P,MODEL+1,
                                    s_pka,
                                    s_ss,
                                    ),
                                self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)


                            pka1 = '%.1f' %(d_fit[0]['c'])

                            fd = open('fits.csv','a')
                            fd.write('%s,%s,%s,%s,%s\n' %(PI,protein,res_ID,pka1,s_pka,))
                            fd.close()

                            if not protein in d_pkas[PI].keys():
                                d_pkas[PI][protein] = {}
                            d_pkas[PI][protein][res_ID] = {'s_pka':s_pka,'pka1':pka1,}

                            os.remove('fit.tit')

                    ##
                    ## plot pka values
                    ##
                    else:

                        if os.path.isfile('png_titcurves/%s/%s_%s.png' %(PI,protein,res_ID,)):
                            continue

                        s_pka = '%.1f' %(float(d_pkas[PI][protein][res_ID]['pka1']))
                        s_ss = 'N/A'
                        P = 'N/A'
                        MODEL = 0
                        if res_ID[:3] in ['LYS','ARG','HIS',]:
                            a = 1.
                            b = -1.
                        elif res_ID[:3] in ['ASP','GLU','TYR',]:
                            a = 0.
                            b = -1.
                        else:
                            print res_ID
                            stop
                        c = float(d_pkas[PI][protein][res_ID]['pka1'])
                        d_fit = {0:{'a':a,'b':b,'c':c,},}
                        title = "'%s   %s %s p=%s pK_as=%i (%s) SS %s'" %(
                            PI,
                            protein.replace('_','\_').upper().replace('DELTA','{/Symbol D}'),
                            res_ID, P, MODEL+1, s_pka, s_ss,
                            ),
                        self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)

                    ##
                    ## jmol per PI and residue
                    ##
                    lines_jmol = ['load ../../pdb/%s.pdb\n' %(d_pdbs[protein],),]

                    if not protein in d_parents.keys(): ## delta+phs
                        ## missing terminal residues
                        if res_no <= 37:
                            res_no_pdb = res_no-6
                        ## gap between residues 38-43
                        elif res_no >= 38:
                            res_no_pdb = res_no-12
                    else:
                        res_no_pdb = res_no
                        
                    lines_jmol += [
                        'color [255,255,255]\n',
                        'select (asp,glu,lys,arg,his,tyr)\n',
                        'color cpk\n',
                        'select (%i:A)\n' %(res_no_pdb),
                        'color [0,255,0]\n',
                        ]
                    fd = open('/local/tc/pkacoop/jmol_scripts/%s_%s.jmol' %(protein,res_ID,), 'w')
                    fd.writelines(lines_jmol)
                    fd.close()

                    ##
                    ## htm per PI and residue
                    ##
                    lines_htm = [
                        '<html>\n',
                        '<head>\n',
                        '  <title>%s %s %s</title>\n' %(PI,protein,res_ID,),
                        ## script source
                        '  <script src="../../jmol/Jmol.js" type="text/javascript"></script>\n',
                        '</head>\n',
                        '<body bgcolor="black">\n',
                        ]
                    lines_htm += [
                        '<font color="#00FF00">%s %s %s</font><br>\n' %(PI, protein, res_ID,),
                        ]
                    if experimental == True:
                        lines_htm += [
                            '<font color="#FFFFFF">measured pKa = %.1f, predicted pKa = %s</font><br>\n' %(pka_exp, d_pkas[PI][protein][res_ID]['s_pka'],),
                            ]
                    else:
                        lines_htm += [
                            '<font color="#FFFFFF">measured pKa = N/A, predicted pKa = %s</font><br>\n' %(d_pkas[PI][protein][res_ID]['s_pka'],),
                            ]
                    lines_htm += [
                        '<font color="#FFFFFF"><a href="../../pdb/%s.pdb">Download PDB file</a></font><br>\n',
                        '<font color="#FFFFFF">The verical black line is the measured pKa value. The vertical gray line is the predicted pKa value.</font><br>\n',
                        '<img src="../../png_titcurves/%s/%s_%s.png"><br>\n' %(PI,protein,res_ID,),
                        '<font color="#00FF00">Residue of interest is marked with green.</font><br>\n',
                        '<font color="#FFFFFF">Other titratable residues are in CPK colors.</font><br>\n',
                        '<form>\n',
                        ' <script type="text/javascript">\n',
                        ## script initialization
                        '  jmolInitialize("../../jmol");\n',
                        ## jmol script location
                        '  jmolApplet(768, "script ../../jmol_scripts/%s_%s.jmol");\n' %(protein,res_ID,),
                        ' </script>\n',
                        '</form>\n',
                        '</body>\n',
                        '</html>\n',
                        ]
                        
                    fd = open('/local/tc/pkacoop/htm_titcurves/%s/%s_%s.htm' %(PI,protein,res_ID,), 'w')
                    fd.writelines(lines_htm)
                    fd.close()

        return d_pkas


    def parse_data(self,d_pkas):

        d_titcurves = {}
        d_residues = {}
        print

        for PI in l_PIs:

            print 'parsing', PI

            l_fn = os.listdir('%s/%s' %(s_path,PI,))
            l_fn.sort()
            d_titcurves[PI] = {}
            if not PI in d_pkas.keys():
                d_pkas[PI] = {}

            for fn in l_fn:

                ## directory (Mike)
                if os.path.isdir('%s/%s/%s' %(s_path,PI,fn,)):
                    continue
                ## gedit ~
                if fn[-1] == '~':
                    continue
                ## ekinprj (Damien)
                if fn[-8:] == '.ekinprj':
                    continue
                ## zip (Jens)
                if fn[-4:] == '.zip':
                    continue
                ## gz (Jana)
                if fn[-3:] == '.gz':
                    continue
                ## xls (Cat)
                if fn[-4:] == '.xls':
                    continue
                ## doc (Ernest)
                if fn[-4:] == '.doc':
                    continue
                ## zip (Emil)
                if fn[-4:] == '.txt':
                    continue
                ## tar (Jan)
                if fn[-4:] == '.tar':
                    continue

                ##
                ## exceptions
                ##
                if PI == 'Mike Word' and fn == '2pwk.csv':
                    continue

                ## csv
                if fn[-4:] == '.csv':
                    pass
                ## None (Jan, Jana)
                if '.' not in fn:
                    if PI in [
                        'Jana K. Shen',
                        'JanJensen_Jdec19','JanJensen_Wdec19','JanJensen_Jjan15','JanJensen_Wjan15',
                        ]:
                        pass
                    else:
                        print PI, fn
                        stop
                        continue

                ## read file
                fd = open('%s/%s/%s' %(s_path,PI,fn,),'r')
                lines = fd.readlines()
                fd.close()

                ##
                ## exceptions
                ##
                if PI == 'Yifan Song':
                    fn = '1stn' ## wt and not delta+phs because of His46 and His124
                    lines = lines[0].split('\r')
                    lines[-1] = lines[-1][:-1] ## equal sign at the end of the last line

                ## protein name from csv file
                protein = fn.lower().replace('.csv','')

                ##
                ## exceptions
                ##
                if PI in ['Jim Warwicker_A','Jim Warwicker_B',]:
                    if protein[:4] == '3bdc':
                        index1 = protein.index('_')+1
                        index2 = index1+protein[index1:].index('_')
                        mutation = protein[index1:index2]
                        if mutation == 'wt':
                            protein = '3bdc'
                        else:
                            protein = mutation+'_delta+phs'
                    else:
                        if protein[4] != 'a':
                            print protein
                            stop
                        protein = protein[:4]

                if protein not in l_proteins:
                    print l_proteins
                    print protein
                    print fn
                    stop

                ##
                ## pka values
                ##
                if d_PIs[PI]['submission'] == 'pka':
                    d_pkas, d_residues = self.parse_pka(
                        PI, protein, d_pkas, lines, d_residues,
                        )

                ##
                ## titration curves
                ##
                elif d_PIs[PI]['submission'] == 'titcurv':
                    d_titcurves, d_residues = self.parse_titcurv(
                        PI, protein, d_titcurves, lines, d_residues,
                        )

        return d_titcurves, d_pkas, d_residues


    def parse_pka(self, PI, protein, d_pkas, lines, d_residues,):

        for line in lines:
            if line == '\n':
                continue
            if line[0] == '#':
                continue
            l = line.strip().split(',')
            pKa = l[1].strip()
            res_ID = l[0].strip()

            l_proteins = [protein]

            ##
            ## exceptions
            ##
            if PI == 'Francesca Milletti':
                index = res_ID.index(':')
                rindex = res_ID.rindex(':')
                res_name = res_ID[:index]
                res_no = int(res_ID[rindex+1:])
                chain = res_ID[index+1:rindex]
                if chain == 'B':
                    continue
                ## res_no
                if len(protein) > 4:
                    if res_no == 49:
                        res_no -= 6
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)
##            if PI == 'Jim Warwicker':
##                if protein == '3d4d':
##                    l_proteins = ['3d4w']
            if PI == 'Qiang Cui':
                if protein == '3bdc':
                    l_proteins = d_yifan[res_ID]
                if protein == '1u9r' and res_ID == 'ASP:A:0066':
                    l_proteins = ['2oxp']
            if PI == 'Ernest Mehler' and protein not in ['1stn',]: ## 1stn correct
                pKa = l[-1]
                res_no = int(res_ID[1:-1])
                chain = 'A'
                res_name = d_res1res3[res_ID[-1].lower()]
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)

            ## check res_ID
            if len(res_ID) not in [10,12,]:
                print PI, protein, res_ID
                stop

            ## loop over proteins in case pka is associated with more than one pdb (Yifan)
            for protein in l_proteins:
                if not protein in d_pkas[PI].keys():
                    d_pkas[PI][protein] = {}
                d_pkas[PI][protein][res_ID] = {'s_pka':pKa,'pka1':pKa,}
                if not protein in d_residues.keys():
                    d_residues[protein] = set()
                d_residues[protein] |= set([res_ID])

        return d_pkas, d_residues


    def parse_titcurv(self, PI, protein, d_titcurves, lines, d_residues,):

        l_proteins = [protein]

        l_pHs = lines[0].split(',')
        ##
        ## exceptions
        ##
        if PI in ['Emil Alexov_1','Emil Alexov_2','Yifan Song',]:
            l_pHs = l_pHs[1:]

        for line in lines[1:]:
            if line == '\n':
                continue
            l = line.split(',')
            if l == ['']:
                continue
            l_charges = l[1:]
            if max(l_charges) > 0 and min(l_charges) < 0:
                stop
            res_ID = l[0].strip()

            ##
            ## exceptions
            ##
            if PI == 'PDB2PKA':
                res_name = res_ID[1:4]
                res_no = int(res_ID[4:])
                chain = res_ID[0]
                ## chain
                chain = 'A'
                ## res_name
                if res_name == 'TER':
                    if res_no < 10:
                        res_name = 'NTERM'
                    if res_no > 100:
                        res_name = 'CTERM'
                ## res_no
                if len(protein) > 4:
                    ## missing terminal residues
                    if res_no <= 37:
                        res_no += 6
                    ## gap between residues 38-43
                    elif res_no >= 38:
                        res_no += 12
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)
            if PI == 'WHATIF':
                res_name = res_ID[1:4]
                res_no = int(res_ID[4:])
                chain = res_ID[0]
                ## chain
                chain = 'A'
                ## res_name
                if res_name == 'TER':
                    if res_no < 10:
                        res_name = 'NTERM'
                    if res_no > 100:
                        res_name = 'CTERM'
                    ## not terminal, next to missing residue
                    if res_no in [42,43,44,45,51,53,54,]:
                        continue
                ## res_no
                if len(protein) > 4:
                    ## missing terminal residues
                    if res_no <= 37:
                        res_no += 6
                    ## gap between residues 38-43
                    elif res_no >= 38:
                        res_no += 12
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)
            if PI in ['Emil Alexov','Emil Alexov_1','Emil Alexov_2',]:
                res_name = res_ID[:-7]
                res_no = int(res_ID[-4:])
                chain = res_ID[-6]
                ## res_no
                if PI != 'Emil Alexov_1' and protein in ['l37k_delta+phs',]:
                    ## missing terminal residues
                    if res_no <= 37:
                        res_no += 6
                    ## gap between residues 38-43
                    elif res_no >= 38:
                        res_no += 12
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)
            if PI == 'Gernot Kieseritzky':
                chain = res_ID[res_ID.index(':')+1]
                if chain == 'B':
                    continue
            if PI == 'Sarah Williams':
                if protein == '1stn' and res_ID == 'HIS:A:0005':
                    res_ID = 'HIS:A:0008'
            if PI == 'Cat Chenal':
                if res_ID[-1] != '_':
                    print protein, res_ID
                    stop
                res_ID = '%s:%s:%s' %(res_ID[:-7],res_ID[-6],res_ID[-5:-1])
                if protein == '3bdc':
                    if res_ID == 'LYS:A:0090':
                        l_proteins = ['a90k_delta+phs']
                    if res_ID == 'GLU:A:0072':
                        l_proteins = ['3ero']
            if PI == 'Jana K. Shen':
                if len(res_ID)-res_ID.rindex(':') < 5:
                    if len(res_ID)-res_ID.rindex(':') != 4:
                        stop
                    res_name = res_ID[:-6]
                    res_no = int(res_ID[-3:])
                    chain = res_ID[-5]
                    res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)
            if PI == 'Yifan Song':
                res_name = res_ID[:3]
                res_no = int(res_ID[3:])
                chain = 'A'
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)
                if res_ID in d_yifan.keys():
                    l_proteins = d_yifan[res_ID]
                else:
                    res1 = d_res3res1[l_sequence[res_no-1]].lower()
                    res2 = d_res3res1[res_name].lower()
                    if res1 != res2: ## mutation
                        l_proteins = ['%s%i%s_delta+phs' %(res1,res_no,res2,)]
            if PI == 'Mike Word':
                chain = res_ID[-6]
                if chain == 'B':
                    continue
                res_no = int(res_ID[-4:])
                if len(protein) != 4:
                    ## missing terminal residues
                    if res_no <= 37:
                        res_no += 6
                    ## gap between residues 38-43
                    elif res_no >= 38:
                        res_no += 12
                res_name = res_ID[:-7]
                ## change res_name
                if len(protein) > 4:
                    if res_no == int(protein[1:-11]):
                        if res_name != d_res1res3[protein[-11]]:
                            print protein, res_ID
                            stop
                            res_name = d_res1res3[protein[-11]]
                res_ID = '%s:%1s:%04i' %(res_name,chain,res_no,)

            ## cehck res_ID
            if len(res_ID) not in [10,12,]:
                print PI, protein, res_ID
                stop

            ## check that there are the same number of pH values and fractional charges
            if len(l_pHs) != len(l_charges):
                print line
                print l_pHs
                print l_charges
                stop

            ## loop over proteins in case pka is associated with more than one pdb (Yifan)
            for protein in l_proteins:
                ## add protein to dic
                if not protein in d_titcurves[PI].keys():
                    d_titcurves[PI][protein] = {}
                ## add residue to dictionary
                d_titcurves[PI][protein][res_ID] = {}

                ## loop over pH values
                for i in range(len(l_pHs)):

                    pH = float(l_pHs[i])
                    if pH > 20:
                        stop1
                    if pH < -10:
                        stop2
                    charge = float(l_charges[i])

                    ##
                    ## exceptions
                    ##
                    if PI == 'Antonio Baptista':
                        if res_ID[:-7] in ['ASP','GLU','TYR','CTERM',]:
                            charge -= 1
                    if PI == 'Mike Word':
                        if res_ID[:-7] in ['ASP','GLU','TYR','CTERM',]:
                            charge *= -1
                    if PI == 'Gernot Kieseritzky':
                        if res_ID[:-7] in ['ASP','GLU','TYR','CTERM',]:
                            charge *= -1
                    if PI == 'Jana K. Shen':
                        if res_ID[:-7] in ['LYS','ARG','HIS','NTERM',]:
                            charge = 1-charge
                        elif res_ID[:-7] in ['GLU','ASP','TYR','CTERM',]:
                            charge *= -1
                        else:
                            print res_ID
                            stop
                    if PI in ['Emil Alexov','Emil Alexov_1','Emil Alexov_2',]:
                        if res_ID[:-7] in ['GLU','ASP','TYR','CTERM',]:
                            charge *= -1
                        ## typing errors?
                        if abs(charge) > 1000:
                            charge /= 10000.
                        if protein == '3c1f' and res_ID == 'LYS:A:0104' and i == 0:
                            charge *= 10

                    if res_ID[:-7] == 'TYR' and charge > 0:
                        print PI, charge
                        stop
                    if res_ID[:-7] in ['ASP','GLU','TYR',] and charge > 0:
                        print PI, protein, res_ID, charge
                        stop
                    if res_ID[:-7] in ['LYS','ARG',] and charge < 0:
                        print PI, protein, res_ID, charge
                        stop

                    d_titcurves[PI][protein][res_ID][pH] = charge

                    continue ## end of loop over pH values

                if not protein in d_residues.keys():
                    d_residues[protein] = set()
                d_residues[protein] |= set([res_ID])

        return d_titcurves, d_residues


    def plot_exp_pka_vs_calc_pka(self,d_residues,d_pkas,):

        print
        print 'plotting experimental pka v calculated pka'

        img_x_min = 265 ## 248, 256 wo title
        img_x_max = 1215 ## 1218, 1224 wo title
        img_y_min = 29 ## 11, 11 wo title, 29 w title
        img_y_max = 979 ## 982, 979 wo title

        lines_data = []
        lines_data_shift = []

        lines_htm_map = []
        lines_htm_map_restricted = []
        lines_htm_map += ['<map id ="Telluride" name="Telluride">\n']
        lines_htm_map_restricted += ['<map id ="Telluride" name="Telluride">\n']

        lines_htm_map_shift = ['<map id ="Telluride_shift" name="Telluride_shift">\n']
        lines_htm_map_shift_restricted = ['<map id ="Telluride_shift" name="Telluride_shift">\n']

        for protein in l_proteins:
            for res_ID in d_residues[protein]:

                res_no  = int(res_ID[-4:])
                res_name = res_ID[:3]

                if protein in d_nonmutants.keys():
                    protein_type = d_nonmutants[protein]
                else:
                    protein_type = 'mutants'

                ## skip if not experimental
                experimental = False
                if res_no in d_experimental[protein_type].keys():
                    if res_name in d_experimental[protein_type][res_no].keys():
                        experimental = True
                        pka_exp = d_experimental[protein_type][res_no][res_name]
                if experimental == False:
                    continue

                ## skip if limit value
                if type(pka_exp).__name__ == 'str':
                    continue

                ## skip if non-mutated residue in mutant protein
                if protein_type == 'mutants' and res_name == l_sequence[res_no-1]:
                    continue

                ## append measured pka value to line
                line_data = '%s %s ' %(res_ID,pka_exp,)
                line_data_shift = '%s %s ' %(res_ID,pka_exp-d_null[res_name],)

                for PI in l_PIs:

                    if not protein in d_pkas[PI].keys():
                        pka_calc = 'N/A'
                    elif not res_ID in d_pkas[PI][protein].keys():
                        pka_calc = 'N/A'
                    else:
                        pka_calc = d_pkas[PI][protein][res_ID]['pka1'] ## temp!!! use best of 2 pkas instead if applicable
                        if pka_calc != 'N/A':
                            pka_calc = '%.1f' %(float(pka_calc))
                    ## append theoretical pka value to line
                    line_data += '%s ' %(pka_calc)
                    if pka_calc == 'N/A':
                        line_data_shift += '%s ' %('N/A',)
                    else:
                        line_data_shift += '%s ' %(float(pka_calc)-d_null[res_name])

                    if pka_calc != 'N/A' and d_PIs[PI]['submission'] == 'titcurv':

                        href = 'htm_titcurves/%s/%s_%s.htm' %(PI, protein,res_ID,)
                        ## absolute
                        line = '<area shape="circle" coords="%i,%i,4" href="%s">\n' %(
                            
        ##                        int(239+(img_x_max-239+1)*(E1A-Min)/float(Max-Min)), ## exlcuding ylabel
        ##                        int(img_y_min+(990-img_y_min+1)*(Max-E1B)/float(Max-Min)), ## excluding xlabel
                            int(img_x_min+(img_x_max-img_x_min+1)*(pka_exp-min_x)/float(max_x-min_x)), ## including ylabel
                            int(img_y_min+(img_y_max-img_y_min+1)*(max_y-float(pka_calc))/float(max_y-min_y)), ## including xlabel
                            href,
                            )
                        lines_htm_map += [line]
                        if '%s%i' %(res_name,res_no,) in l_allowed:
                            lines_htm_map_restricted += [line]

                        ## shifted
                        line = '<area shape="circle" coords="%i,%i,4" href="htm_titcurves/%s/%s_%s.htm">\n' %(
                            
        ##                        int(239+(img_x_max-239+1)*(E1A-Min)/float(Max-Min)), ## exlcuding ylabel
        ##                        int(img_y_min+(990-img_y_min+1)*(Max-E1B)/float(Max-Min)), ## excluding xlabel
                            int(img_x_min+(img_x_max-img_x_min+1)*((pka_exp-d_null[res_name])-min_x_shift)/float(max_x_shift-min_x_shift)), ## including ylabel, excluding title
                            int(img_y_min+(img_y_max-img_y_min+1)*(max_y_shift-(float(pka_calc)-d_null[res_name]))/float(max_y_shift-min_y_shift)), ## including xlabel, excluding title
                            PI, protein,res_ID,
                            )
                        lines_htm_map_shift += [line]
                        if '%s%i' %(res_name,res_no,) in l_allowed:
                            lines_htm_map_shift_restricted += [line]

                ## link to combined plot on diagonal

                ## absolute
                pka_calc = pka_exp
                lines_htm_map += [
                    '<area shape="circle" coords="%i,%i,4" href="htm_titcurves/combined/%s_%s.htm">\n' %(
                        int(img_x_min+(img_x_max-img_x_min+1)*(pka_exp-min_x)/float(max_x-min_x)), ## including ylabel
                        int(img_y_min+(img_y_max-img_y_min+1)*(max_y-float(pka_calc))/float(max_y-min_y)), ## including xlabel
                        protein,res_ID,
                        )
                    ]

                ## absolute, restricted
                if '%s%i' %(res_name,res_no,) in l_allowed:
                    lines_htm_map_restricted += [
                        '<area shape="circle" coords="%i,%i,4" href="htm_titcurves/combined_restricted/%s_%s.htm">\n' %(
                            int(img_x_min+(img_x_max-img_x_min+1)*(pka_exp-min_x)/float(max_x-min_x)), ## including ylabel
                            int(img_y_min+(img_y_max-img_y_min+1)*(max_y-float(pka_calc))/float(max_y-min_y)), ## including xlabel
                            protein,res_ID,
                            )
                        ]
                ## shifted pka
                lines_htm_map_shift += [
                    '<area shape="circle" coords="%i,%i,4" href="htm_titcurves/combined/%s_%s.htm">\n' %(
                        int(img_x_min+(img_x_max-img_x_min+1)*((pka_exp-d_null[res_name])-min_x_shift)/float(max_x_shift-min_x_shift)), ## including ylabel
                        int(img_y_min+(img_y_max-img_y_min+1)*(max_y_shift-(float(pka_calc)-d_null[res_name]))/float(max_y_shift-min_y_shift)), ## including xlabel
                        protein,res_ID,
                        )
                    ]
                if '%s%i' %(res_name,res_no,) in l_allowed:
                    lines_htm_map_shift_restricted += [
                        '<area shape="circle" coords="%i,%i,4" href="htm_titcurves/combined_restricted/%s_%s.htm">\n' %(
                            int(img_x_min+(img_x_max-img_x_min+1)*((pka_exp-d_null[res_name])-min_x_shift)/float(max_x_shift-min_x_shift)), ## including ylabel
                            int(img_y_min+(img_y_max-img_y_min+1)*(max_y_shift-(float(pka_calc)-d_null[res_name]))/float(max_y_shift-min_y_shift)), ## including xlabel
                            protein,res_ID,
                            )
                        ]

                line_data += '\n'
                line_data_shift += '\n'

                lines_data += [line_data]
                lines_data_shift += [line_data_shift]

        lines_htm_map += ['</map>\n',]
        lines_htm_map_restricted += ['</map>\n',]

        fd = open('gnuplot.dat','w')
        fd.writelines(lines_data)
        fd.close()

        fd = open('gnuplot_shift.dat','w')
        fd.writelines(lines_data_shift)
        fd.close()

        for prefix in ['gnuplot','gnuplot_shift',]:

            lines_gnuplot = []
            lines_gnuplot += [
                'set terminal postscript eps enhanced color "Helvetica" 18\n',
                'set output "gnuplot.ps"\n',
                'set size square\n',
                'set size 4,4\n', ## scale size
                'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
                ]
            if prefix == 'gnuplot':
                lines_gnuplot += [
                    'set xlabel "experimental pKa"\n',
                    'set ylabel "predicted pKa"\n',
                    'set title "experimental pKa v predicted pKa"\n', ## recalculate plot dimensions if title added!!!
                    ]
            else:
                lines_gnuplot += [
                    'set xlabel "experimental pKa - model pKa"\n',
                    'set ylabel "predicted pKa - model pKa"\n',
                    'set title "experimental pKa shift v predicted pKa shift"\n', ## recalculate plot dimensions if title added!!!
                    ]
            lines_gnuplot += [
        ##        'unset key\n',
                'f(x) = x\n',
                'g(x) = x+1\n',
                'h(x) = x-1\n',
                'i(x) = x+2\n',
                'j(x) = x-2\n',
                'k(x) = x+3\n',
                'l(x) = x-3\n',
                'm(x) = x+4\n',
                'n(x) = x-4\n',
                ]

            if prefix == 'gnuplot':
                line = 'plot [%s:%s][%s:%s] ' %(min_x,max_x,min_y,max_y,)
            else:
                line = 'plot [%s:%s][%s:%s] ' %(min_x_shift,max_x_shift,min_y_shift,max_y_shift,)
            line += 'f(x) lw 5 lc 0 lt 1 t "", g(x) lc 0 lw 5 lt 0 t "", h(x) lc 0 lw 5 lt 0 t "", i(x) lc 0 lw 5 lt 0 t "", j(x) lc 0 lw 5 lt 0 t "", k(x) lc 0 lw 5 lt 0 t "", l(x) lc 0 lw 5 lt 0 t "", m(x) lc 0 lw 5 lt 0 t "", n(x) lc 0 lw 5 lt 0 t "", '
            line += '"%s.dat" u 2:2 pt 7 ps 2 lt rgb "#000000" t "experimental", ' %(prefix)
    ##        line += '"gnuplot.dat" u 2:2:1 lt rgb "#000000" w labels t "", '
            for i in range(len(l_PIs)):
                PI = l_PIs[i]
                line += '"%s.dat" u 2:%i pt %i ps 2 lt rgb "#%6s" t "%s", ' %(prefix,i+3,d_PIs[PI]['pt'],d_PIs[PI]['color'],l_PIs[i],)
    ##            line += '"gnuplot.dat" u 2:($%i+0.1):1 lt rgb "#%6s" w labels t "", ' %(i+3,d_PIs[PI]['color'],)
            line = line[:-2]+'\n'
            lines_gnuplot += [line]

            fd = open('gnuplot.settings','w')
            fd.writelines(lines_gnuplot)
            fd.close()

            os.system('gnuplot gnuplot.settings')
            os.system('convert gnuplot.ps %s.png' %(prefix))
            os.remove('gnuplot.ps')
            os.remove('gnuplot.settings')

        return lines_htm_map, lines_htm_map_restricted, lines_htm_map_shift, lines_htm_map_shift_restricted


    def checkdata(self, d_titcurv, PI, protein, res_ID, pka_exp, bool_plot=False,):

        d_fit = None
        MODEL = 3
        P = 'N/A'
        title = '"%s   %s %s not a full titration"' %(PI, protein, res_ID,)
        full_titration = True

        if max(d_titcurv.values())-min(d_titcurv.values()) < 1.-f_fulltit:
            full_titration = False
        
##        ## not a full titration range
##        if max(d_titcurv.values()) == 0 and min(d_titcurv.values()) > -1+f_fulltit:
####            if bool_plot == True:
####                plot_titcurve(pdb,d_fit,MODEL,i_res_ID,res_ID,title,dn,)
####            return False
##            full_titration = False
##        if max(d_titcurv.values()) == 1 and min(d_titcurv.values()) > 0+f_fulltit:
####            if bool_plot == True:
####                plot_titcurve(pdb,d_fit,MODEL,i_res_ID,res_ID,title,dn,)
####            return False
##            full_titration = False
##        if max(d_titcurv.values()) < 0-f_fulltit and min(d_titcurv.values()) == -1:
####            if bool_plot == True:
####                plot_titcurve(pdb,d_fit,MODEL,i_res_ID,res_ID,title,dn,)
####            return False
##            full_titration = False
##        if max(d_titcurv.values()) < 1-f_fulltit and min(d_titcurv.values()) == 0:
####            if bool_plot == True:
####                plot_titcurve(pdb,d_fit,MODEL,i_res_ID,res_ID,title,dn,)
####            return False
##            full_titration = False

        ## nothing to fit, full charge at all pH values
        if max(d_titcurv.values()) == 1 and min(d_titcurv.values()) == 1:
            self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)
            return False, False
        if max(d_titcurv.values()) == -1 and min(d_titcurv.values()) == -1:
            print d_titcurv
            print protein, res_ID
            stop2
            self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)
            return False, False
        if max(d_titcurv.values()) == 0 and min(d_titcurv.values()) == 0:
            self.plot_titcurve(protein,d_fit,MODEL,res_ID,title,PI,pka_exp,)
            return False, False

        return True, full_titration


    def get_model(self,d_titcurv,):

        ##
        ## set models
        ##
        d_fit = {
            0:{
                'function':'a+b/(1+10**(c-x))\n',
                'variables':'a,b,c'
                },
            1:{
                'function':'a+b/(1+10**(c-x))+d/(1+10**(e-x))\n',
                'variables':'a,b,c,d,e'
                },
            2:{
                'function':'a+b/(1+10**(c-x))+d/(1+10**(e-x))+f/(1+10**(g-x))\n',
                'variables':'a,b,c,d,e,f,g'
                },
            }

        for i in range(2):
            model1 = i
            model2 = i+1
            ##
            ## fit data for each model
            ##
            if not 'SS' in d_fit[model1].keys():
                d_fit[model1] = self.gnuplot_fit(d_fit[model1],d_titcurv,)
            if not 'SS' in d_fit[model2].keys():
                try:
                    d_fit[model2] = self.gnuplot_fit(d_fit[model2],d_titcurv,)
                ## Singular matrix in Invert_RtR
                except:
                    MODEL = model1
                    P = 'N/A'
                    break
            ss1 = d_fit[model1]['SS']
            ss2 = d_fit[model2]['SS']
            n1 = d_fit[model1]['n']
            n2 = d_fit[model2]['n']
            if n1 != n2:
                stop
            df1 = n1-len(d_fit[model1]['variables'].split(','))
            df2 = n2-len(d_fit[model2]['variables'].split(','))

            ## break 1
            if ss2 > ss1:
                print i, 'F is negative. Complicated model is worse than simple model'
                MODEL = model1
                SS2 = ss2
                SS1 = ss1
                P = 'N/A'
                break

            print 'h(x) = %s+%s/(1+10**(%s-x))' %(d_fit[0]['a'],d_fit[0]['b'],d_fit[0]['c'],)
            print 'i(x) = %s+%s/(1+10**(%s-x))+%s/(1+10**(%s-x))' %(d_fit[1]['a'],d_fit[1]['b'],d_fit[1]['c'],d_fit[1]['d'],d_fit[1]['e'],)
            if 'a' in d_fit[2].keys():
                print 'j(x) = %s+%s/(1+10**(%s-x))+%s/(1+10**(%s-x))+%s/(1+10**(%s-x))\n' %(d_fit[2]['a'],d_fit[2]['b'],d_fit[2]['c'],d_fit[2]['d'],d_fit[2]['e'],d_fit[2]['f'],d_fit[2]['g'],),
            F = ((ss1-ss2)/ss2)/((df1-df2)/float(df2))
            p = statistics.fdist(F,df1-df2,df1)
            print 'p', p, 'ss', ss1,ss2, 'df', df1,df2,df1-df2, 'F', F

            if p < 0.01:

                ## break2
                if ss1 < 0.001: ## (2cgaAHIS40 fit to 1 pKa if < 0.06) (2aas with 2 pka fit to 1 pka if > .005)
                    print i, 'SS of simple model is very small (%s %s) and complicated model is not justified.' %(ss1,ss2)
                    MODEL = model1
                    SS2 = ss2
                    SS1 = ss1
                    P = 'N/A'
                    break

                l_spans = []
                for k in ['b','d','f',]:
                    if d_fit[model2][k] != None:
                        l_spans += [abs(d_fit[model2][k])]
                l_spans.sort()
                l_spans = l_spans[:-1]
                ## break 3
                if max(l_spans) < .10:
                    print i, 'Span is too small to justify complicated model.'
                    MODEL = model1
                    SS2 = ss2
                    SS1 = ss1
                    P = 'N/A'
                    break

                l_pkas = []
                for k in ['c','e','g',]:
                    if d_fit[model2][k] != None:
                        l_pkas += [d_fit[model2][k]]
                l_pkas.sort()
                ## break 4
                if abs(l_pkas[0]-l_pkas[1]) < 1.7:
                    print 'pkas too close', l_pkas
                    MODEL = model1
                    SS2 = ss2
                    SS1 = ss1
                    P = 'N/A'
                    break
            
            ## better fit if more variables
            if p < 0.01:
                P = p
                MODEL = model1
                SS2 = ss2
                SS1 = ss1
                continue
            else:
                ## break 5
                MODEL = model1
                SS2 = ss2
                SS1 = ss1
                P = p
                break

        return MODEL, P, d_fit


    def gnuplot_fit(self,d_fit,d_titcurv,):

        if not 'a' in d_fit.keys():
            d_fit['a'] = 0.01
        d_fit['c'] = 5
        d_fit['e'] = 10
        d_fit['g'] = 15
        if not 'b' in d_fit.keys():
            d_fit['b'] = -1
        d_fit['d'] = -.5
        d_fit['f'] = -.33
        lines = [
            'y(x) = %s' %(d_fit['function']),
            'a = %f\n' %(d_fit['a']),
            'c = %f\n' %(d_fit['c']),
            'e = %f\n' %(d_fit['e']),
            'g = %f\n' %(d_fit['g']),
            'b = %f\n' %(d_fit['b']),
            'd = %f\n' %(d_fit['d']),
            'f = %f\n' %(d_fit['f']),
            'set fit logfile "fit.log"\n',
            'fit y(x) "fit.tit" u 1:2 via %s\n' %( ## not weighted
    ##        'fit y(x) "%s.tit" u 1:%s:%s via %s\n' %( ## weighted
                d_fit['variables'], ## not weighted
    ##            i_res_ID*2+2,i_res_ID*2+3,d_fit['variables'], ## weighted
                ),
            ]

        fd = open('gnufit.src','w')
        fd.writelines(lines)
        fd.close()

        print '**********************'
        os.system('gnuplot gnufit.src')
        print '**********************'

        l = os.popen('tail -n21 fit.log').readlines()
        os.remove('fit.log') ## remove log, otherwise appended to during next fit?
        index = l.index('=======================            ==========================\n')
        index1 = index+2
        index2 = index1+l[index1:].index('\n')
        d_fit['b'] = None
        d_fit['d'] = None
        d_fit['f'] = None
        d_fit['c'] = None
        d_fit['e'] = None
        d_fit['g'] = None
        for s in l[index1:index2]:
            k = s.split()[0]
            v = float(s.split()[2])
##            if k == 'a':
##                if abs(v-d_fit['a']) > .1:
##                    print l[index1:index2]
##                    print v, d_fit['a']
##                    stop
##                continue
            d_fit[k] = v

        os.remove('gnufit.src')

        ##
        ## calculate sum of squares
        ##
        SS = []
        for pH in d_titcurv.keys():
            x = pH
            a = d_fit['a']
            b = d_fit['b']
            c = d_fit['c']
            d = d_fit['d']
            e = d_fit['e']
            f = d_fit['f']
            g = d_fit['g']
            fraction_pdb2pka = d_titcurv[pH]
    ##        if fraction_pdb2pka in [-1,0,1,]:
    ##            continue
            fraction_fit = eval(d_fit['function'])
            SS += [(fraction_pdb2pka-fraction_fit)**2]
        d_fit['SS'] = sum(SS)/len(SS)
        d_fit['n'] = len(SS)

        return d_fit


    def plot_titcurve(self,pdb,d_fit,MODEL,res_ID,title,PI,pka_exp,):

        print 'plotting titcurve', PI, pdb

        lines = [
            'set terminal postscript eps enhanced color "Helvetica" 30\n',
            'set output "%s_%s.ps"\n' %(pdb,res_ID,),
            'set size 3,3\n',
            'set grid\n',
            'set mxtics 1\n',
            'set xtics 1\n',
            'set xlabel "pH"\n',
            'set ylabel "charge"\n',
            'set title %s\n' %(title),
            ]

        ## calculated pka vertical line
        if MODEL in range(3):
            for i in range(MODEL+1):
                pka = d_fit[MODEL][['c','e','g',][i]]
                lines += [
                    'set arrow from %f, graph 0 to %f, graph 1 nohead\n' %(pka,pka,),
                    ]

        ## experimental pka vertical line
        if pka_exp:
            lines += ['set arrow from %f, graph 0 to %f, graph 1 nohead lw 5 lc 0\n' %(pka_exp,pka_exp,),]

        s = 'plot [-10:20][-1:1] '
        ## plot titcurv datapoints
        if d_PIs[PI]['submission'] == 'titcurv':
            s += '"fit.tit" u 1:2 pt 7 ps 2 t "data points", '
        if MODEL in [0,1,2,]:
            lines += [
                'h(x) = %s+%s/(1+10**(%s-x))\n' %(d_fit[0]['a'],d_fit[0]['b'],d_fit[0]['c'],),
                ]
            ## 2 pKa valuesz
            if 1 in d_fit.keys():
                lines += [
                    'i(x) = %s+%s/(1+10**(%s-x))+%s/(1+10**(%s-x))\n' %(d_fit[1]['a'],d_fit[1]['b'],d_fit[1]['c'],d_fit[1]['d'],d_fit[1]['e'],),
                    ]
                s += 'h(x) lw 4 t "1 pKa", '
                s += 'i(x) lw 4 t "2 pKas", '
            ## pka submitted
            else:
                s += 'h(x) lw 4 t "theoretical HH curve", '
        if MODEL == 2:
            lines += [
                'j(x) = %s+%s/(1+10**(%s-x))+%s/(1+10**(%s-x))+%s/(1+10**(%s-x))\n' %(d_fit[2]['a'],d_fit[2]['b'],d_fit[2]['c'],d_fit[2]['d'],d_fit[2]['e'],d_fit[2]['f'],d_fit[2]['g'],),
                ]
            s += 'j(x) lw 4 t "3 pKas"'
        s = s[:-2]+'\n'
        lines += [s]

        fd = open('gnuplot_%s_%s.src' %(pdb,res_ID,),'w')
        fd.writelines(lines)
        fd.close()

        os.system('gnuplot gnuplot_%s_%s.src' %(pdb,res_ID,))
        os.system('convert %s_%s.ps %s_%s.png' %(pdb,res_ID,pdb,res_ID,))
        if not os.path.isdir('png_titcurves/%s' %(PI)):
            os.mkdir('png_titcurves/%s' %(PI))
        os.system('mv %s_%s.png png_titcurves/%s/.' %(pdb,res_ID,PI.replace(' ','\ '),))

        os.remove('%s_%s.ps' %(pdb,res_ID,))
        os.remove('gnuplot_%s_%s.src' %(pdb,res_ID,))

        return


    def __init__(self):
        
##        for PI,color in d_colors.items():
##            if color == 0:
##                s_hex = '000000'
##            else:
##                h = 0+((color-1)%5)*(240./5)
##                s = 240
##                l = 90+((color-1)/55555555)*30
##                r,g,b = self.hsl2rgb(h,s,l,)
##                r *= 255
##                g *= 255
##                b *= 255
##                r = hex(int(round(r,0)))
##                g = hex(int(round(g,0)))
##                b = hex(int(round(b,0)))
##                s_hex = r[r.index('x')+1:].zfill(2)+g[g.index('x')+1:].zfill(2)+b[b.index('x')+1:].zfill(2)
##                s_hex = s_hex.upper()
##            d_colors[PI] = s_hex
##
##        self.d_colors = d_colors

        return


if __name__ == '__main__':
    instance_telluride = Telluride()
    instance_telluride.main()

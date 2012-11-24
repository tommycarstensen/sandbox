import math, os, statistics, sys
##import numpy ## do not use numpy. it's slow! more than 3 times slower!

l_backbone_atoms = [
    'N','H','H1','H2','H3', ## H3 if N-terminal (H,HXT in PDB)
    'CA','HA','HA2', ## HA2 if Gly (HA3 in PDB)
    'C','O','OC1','OC2', ## OC1,OC2 if C-terminal (O,OXT in PDB)
    ]

## conformational state 35,52, protonation state 35,52
l_states = [
    ## 35 charged protonation state
        ## 35 neutral conformational state
            ## 52 charged protonation state
    'NEUNEUCHACHA',
    'NEUCHACHACHA',
            ## 52 neutral protonation state
    'NEUNEUCHANEU',
    'NEUCHACHANEU',
        ## 35 charged conformational state
            ## 52 charged protonation state
    'CHANEUCHACHA',
    'CHACHACHACHA',
            ## 52 neutral protonation state
    'CHANEUCHANEU',
    'CHACHACHANEU',
    ## 35 neutral protonation state
        ## 35 neutral conformational state
            ## 52 charged protonation state
    'NEUNEUNEUCHA',
    'NEUCHANEUCHA',
            ## 52 neutral protonation state
    'NEUNEUNEUNEU',
    'NEUCHANEUNEU',
        ## 35 charged conformational state
            ## 52 charged protonation state
    'CHANEUNEUCHA',
    'CHACHANEUCHA',
            ## 52 neutral protonation state
    'CHANEUNEUNEU',
    'CHACHANEUNEU',

    ]

l_colors = [
    '000000', ## black
    'FF0000', ## red
    'FF8000', ## orange
    'F0F000', ## yellow
    '00FF00', ## green 2*
    'FF00FF', ## purple
    '0000FF', ## blue 3*
    '808080', ## grey
    ]

l_colors_datapoints = [
    '808080', ## black
    'FF8080', ## red
    'FFC080', ## orange
    'FFFF80', ## yellow
    '80FF80', ## green 2*
    'FF80FF', ## purple
    '8080FF', ## blue 3*
    'C0C0C0', ## grey
    ]

d_colors = {
##    'CHACHACHACHA':'FF0000',
##    'NEUCHACHACHA':'FF4040',
##    'CHANEUCHACHA':'FF8080',
##    'NEUNEUCHACHA':'FFC0C0',
##    'CHACHANEUCHA':'80FF00',
##    'NEUCHANEUCHA':'80FF40',
##    'CHANEUNEUCHA':'80FF80',
##    'NEUNEUNEUCHA':'80FFC0',
##    'CHACHACHANEU':'00FFFF',
##    'NEUCHACHANEU':'40FFFF',
##    'CHANEUCHANEU':'80FFFF',
##    'NEUNEUCHANEU':'C0FFFF',
##    'CHACHANEUNEU':'8000FF',
##    'NEUCHANEUNEU':'8040FF',
##    'CHANEUNEUNEU':'8080FF',
##    'NEUNEUNEUNEU':'80C0FF',

##    'CHACHACHACHA':'FF0000',
##    'NEUCHACHACHA':'FF8000',
##    'CHANEUCHACHA':'FFFF00',
##    'NEUNEUCHACHA':'80FF00',
##    'CHACHANEUCHA':'00FF00',
##    'NEUCHANEUCHA':'00FF80',
##    'CHANEUNEUCHA':'00FFFF',
##    'NEUNEUNEUCHA':'0080FF',
##    'CHACHACHANEU':'0000FF',
##    'NEUCHACHANEU':'8000FF',
##    'CHANEUCHANEU':'FF00FF',
##    'NEUNEUCHANEU':'FF0080',
##    'CHACHANEUNEU':'000000',
##    'NEUCHANEUNEU':'808080',
##    'CHANEUNEUNEU':'80808F',
##    'NEUNEUNEUNEU':'C0C0C0',

##    'CHACHANEUCHA':'000000', ## black
##    'NEUCHANEUCHA':'FF0000', ## red
##    'CHANEUNEUCHA':'FF8000', ## orange
##    'NEUNEUNEUCHA':'FFFF00', ## yellow
##    'CHACHANEUNEU':'00FF00', ## green 2*
##    'NEUCHANEUNEU':'FF00FF', ## purple
##    'CHANEUNEUNEU':'0000FF', ## blue 3*
##    'NEUNEUNEUNEU':'808080', ## grey
##    'CHACHACHACHA':'000000',
##    'NEUCHACHACHA':'FF0000',
##    'CHANEUCHACHA':'FF8000',
##    'NEUNEUCHACHA':'FFFF00',
##
##    'CHACHACHANEU':'00FF00',
##    'NEUCHACHANEU':'FF00FF',
##    'CHANEUCHANEU':'0000FF',
##    'NEUNEUCHANEU':'808080',
    }

for i in range(len(l_states)):
    state = l_states[i]
    color = l_colors[i%8]
    d_colors[state] = {'lc':color,'pc':l_colors_datapoints[i%8],}

def main():

    cwd = os.getcwd()[-6:]

    ## parse charges
    d_charges = parse_charges_from_topology_file()

    ## calculate energies
    d_coords = parse_coords_from_trajectory_file(d_charges,cwd,)

    ## calculate statistics
    calculate_averages_and_plot(cwd,)

    return


def parse_charges_from_topology_file():

    d_charges = {}
    for topology in ['NEUNEU','NEUCHA','CHANEU','CHACHA',]:

        fd = open('/local/tc/MD_2vb1/amber99sb_CYM/%s/2vb1.top' %(topology),'r')
    ##    fd = open('../MD_2vb1_Glu35Asp52_%s/2vb1.top' %(topology),'r')
        lines = fd.readlines()
        fd.close()

        d_charges[topology] = {}
        for i in range(len(lines)):
            if lines[i].strip() == '[ atoms ]':
                for j in range(i+2,len(lines)):
                    if lines[j].strip() == '':
                        break
                    charge = float(lines[j].split()[6])
                    res_no = int(lines[j].split()[2])
                    atom_name = lines[j].split()[4]
                    if not res_no in d_charges[topology].keys():
                        d_charges[topology][res_no] = {}
                    if atom_name in d_charges[topology][res_no].keys():
                        stop
                    d_charges[topology][res_no][atom_name] = charge
    ##                atom_no = int(lines[j].split()[0])
    ##                d_charges[atom_no] = charge

    return d_charges


def calculate_averages_and_plot(cwd,):

    import statistics

    print 'calculate averages and plot'

    for topology in ['NEUNEU','NEUCHA','CHANEU','CHACHA',]:

        print 'protonation state', topology

        fd = open('energies_%s.txt' %(topology),'r')
        lines = fd.readlines()
        fd.close()

        lines2 = []
        l_asp52 = []
        l_protein = []
        l_chloride = []
        l_water = []
        l_sum = []
        l_sum_excl_ions = []
        for line in lines:
            i = int(line.split()[0])
            if i % 1000 == 0:
                print 'average', topology, i
    ##        if i < 100:
    ##            continue
            l_asp52 += [float(line.split()[1])]
            l_protein += [float(line.split()[2])]
            l_chloride += [float(line.split()[3])]
            l_water += [float(line.split()[4])]
            l_sum += [
                float(line.split()[1])+
                float(line.split()[2])+
                float(line.split()[3])+
                float(line.split()[4])
                ]
            l_sum_excl_ions += [
                float(line.split()[1])+
                float(line.split()[2])+
                float(line.split()[4])
                ]

            lines2 += [
                '%i %f %f %f %f %f %f\n' %(
                    i, ## 1
                    sum(l_asp52)/len(l_asp52),
                    sum(l_protein)/len(l_protein),
                    sum(l_chloride)/len(l_chloride),
                    sum(l_water)/len(l_water),
                    sum(l_sum)/len(l_sum), ## 13
                    sum(l_sum_excl_ions)/len(l_sum_excl_ions), ## 12
                    )
                ]

        fd = open('energies_averages_%s.txt' %(topology),'w')
        fd.writelines(lines2)
        fd.close()

        fd = open('energies_averages_%s.txt' %(topology),'r')
        lines2 = fd.readlines()
        fd.close()
        average = float(lines2[-1].split()[5])
        print '******** average', average
        ## calculate rmsd
        l_diff = []
        for i in range(len(lines)):
            Sum = float(line.split()[1])+float(line.split()[2])+float(line.split()[3])+float(line.split()[4])
            l_diff += [Sum-average]
        rmsd = statistics.do_rmsd(l_diff)
        print '******** rmsd', rmsd

        if len(l_sum) > 0:
            print 'INCLUDING IONS'
            print 'correl asp52', statistics.correlation(l_sum,l_asp52)
            print 'correl protein', statistics.correlation(l_sum,l_protein)
            print 'correl chloride', statistics.correlation(l_sum,l_chloride)
            print 'correl water', statistics.correlation(l_sum,l_water)
            print 'EXCLUDING IONS'
            print 'correl asp52', statistics.correlation(l_sum_excl_ions,l_asp52)
            print 'correl protein', statistics.correlation(l_sum_excl_ions,l_protein)
            print 'correl chloride', statistics.correlation(l_sum_excl_ions,l_chloride)
            print 'correl water', statistics.correlation(l_sum_excl_ions,l_water)

    ##
    ## combined plot 3 (4 conformational states x 4 protonation states and their averages)
    ##
    ## orange=black, blue=green, yellow=red, grey=purple
    ## *NEUNEUCHA*NEU = *NEUNEUCHA*CHA (ion,water,protein)
    ## *CHANEUCHA*NEU = *CHANEUCHA*CHA (ion,water,protein)
    ## *NEUCHACHA*CHA = *NEUCHACHA*NEU (ion,water)
    ## *CHACHACHA*NEU = *CHACHACHA*CHA (ion,water)
    ## overlaps - water, ions, (protein)
    ## CHACHANEUCHA not overlap when protein
    for s_col,title,suffix,y1,y2 in [
        ['$2+$3+$4+$5','energies of 4 conformational states at 4 different protonation states (all terms)','1all',-700,250,],
        ['$2+$3+$5','energies of 4 conformational states at 4 different protonation states (excluding ions)','1exclions',-700,250,],
        ['$3+$4+$5','energies of 4 conformational states at 4 different protonation states (excluding Asp52)','1exclasp52',-700,250,],
        ['$2','energies of 4 conformational states at 4 different protonation states (Asp52)','2asp52',-700,250,],
        ['$2','energies of 4 conformational states at 4 different protonation states (Asp52)','2asp52_zoom1',-80,23,],
        ['$2','energies of 4 conformational states at 4 different protonation states (Asp52)','2asp52_zoom2',23,160,],
        ['$3','energies of 4 conformational states at 4 different protonation states (protein)','3protein',-700,250,],
        ['$4','energies of 4 conformational states at 4 different protonation states (ions)','4ions',-700,250,],
        ['$5','energies of 4 conformational states at 4 different protonation states (water)','5water',-700,250,],
        ['$5','energies of 4 conformational states at 4 different protonation states (water)','5water_zoom',-80,40,],
        ]:
        print 'combined plot 16 states', suffix
        lines = [
            'set terminal postscript eps enhanced color "Helvetica" 32\n',
            'set size 3,3\n',
            'set output "combined_16states.ps"\n',
            'set xlabel "t / ps"\n',
            'set ylabel "E / kT"\n',
            'set title "%s"\n' %(title),
            ]
    ##    line = 'plot [0:][%s:%s]' %(Min,Max,)
        line = 'plot [0:30000][%s:%s]' %(y1,y2,)
        ## data points
        for i_state in range(16):
            state = l_states[i_state]
    ##        pt = [6,7,4,5,12,13][i_state % 4]
            if i_state < 8:
                pt = 7
                ps = 1
            else:
                pt = 5
                ps = 1
            ## data points
            line += '"../%s/energies_%s.txt" u 1:(%s) lc rgb "#%6s" ps %i pt %i t "%s", ' %(state[:6],state[-6:],s_col,d_colors[state]['pc'],ps,pt,state,)
        ## lines
        for i_state in range(16):
            if i_state < 8:
                if i_state in [0,1,4,5,]:
                    lw = 16
                else:
                    lw = 12
            else:
                lw = 4
            state = l_states[i_state]
            ## lines (averages)
            line += '"../%s/energies_averages_%s.txt" u 1:(%s) w l lt 1 lc rgb "#%6s" lw %i t "%s", ' %(state[:6],state[-6:],s_col,d_colors[state]['lc'],lw,state,)
        line = line[:-2]+'\n'
        lines += [line]
        
        fd = open('gnu.set','w')
        fd.writelines(lines)
        fd.close()

        os.system('gnuplot gnu.set')
        os.system('convert combined_16states.ps combined_16states_%s.png' %(suffix))

    ##
    ## combined plot 2 (2 states with individual terms and their averages)
    ##
    for combination in [['CHACHA','CHANEU',],['NEUCHA','NEUNEU',],]:

        print 'combined plot', combination
        lines = [
            'set terminal postscript eps enhanced color "Helvetica" 32\n',
            'set size 3,3\n',
            'set output "combined.ps"\n',
            'set xlabel "t / ps"\n',
            'set ylabel "E / kT"\n',
            'set title "%s"\n' %('%s v %s' %(combination[0],combination[1],)),
            ]
        line = 'plot [0:][-500:100]'
        ## data points
        line += '"../%s/energies_%s.txt" u 1:%s lc 3 t "%s protein", ' %(combination[0],combination[0],'($3)',combination[0],)
        line += '"../%s/energies_%s.txt" u 1:%s lc 4 t "%s protein", ' %(combination[1],combination[1],'($3)',combination[1],)
        line += '"../%s/energies_%s.txt" u 1:%s lc 5 t "%s water", ' %(combination[0],combination[0],'($5)',combination[0],)
        line += '"../%s/energies_%s.txt" u 1:%s lc 6 t "%s water", ' %(combination[1],combination[1],'($5)',combination[1],)
        line += '"../%s/energies_%s.txt" u 1:%s lc 1 t "%s Asp52", ' %(combination[0],combination[0],'($2)',combination[0],)
        line += '"../%s/energies_%s.txt" u 1:%s lc 2 t "%s Asp52", ' %(combination[1],combination[1],'($2)',combination[1],)
        line += '"../%s/energies_averages_%s.txt" u 1:%s w l lt 1 lc 3 lw 16 t "%s protein average", ' %(combination[0],combination[0],'($3)',combination[0],)
        line += '"../%s/energies_averages_%s.txt" u 1:%s w l lt 1 lc 4 lw 16 t "%s protein average", ' %(combination[1],combination[1],'($3)',combination[1],)
        line += '"../%s/energies_averages_%s.txt" u 1:%s w l lt 1 lc 5 lw 16 t "%s water average", ' %(combination[0],combination[0],'($5)',combination[0],)
        line += '"../%s/energies_averages_%s.txt" u 1:%s w l lt 1 lc 6 lw 16 t "%s water average", ' %(combination[1],combination[1],'($5)',combination[1],)
        line += '"../%s/energies_averages_%s.txt" u 1:%s w l lt 1 lc 1 lw 16 t "%s Asp52 average", ' %(combination[0],combination[0],'($2)',combination[0],)
        line += '"../%s/energies_averages_%s.txt" u 1:%s w l lt 1 lc 2 lw 16 t "%s Asp52 average", ' %(combination[1],combination[1],'($2)',combination[1],)
    ##    line += '"../NEUCHA/energies_NEUCHA.txt" u 1:%s lc 7 t "NEUCHA Asp52", ' %('($2)',)
    ##    line += '"../NEUNEU/energies_NEUNEU.txt" u 1:%s lc 8 t "NEUNEU Asp52", ' %('($2)',)
    ##    line += '"../NEUCHA/energies_NEUCHA.txt" u 1:%s lc 9 t "NEUCHA protein", ' %('($3)',)
    ##    line += '"../NEUNEU/energies_NEUNEU.txt" u 1:%s lc 10 t "NEUNEU protein", ' %('($3)',)
        line = line[:-2]+'\n'
        lines += [line]
        
        fd = open('gnu.set','w')
        fd.writelines(lines)
        fd.close()

        os.system('gnuplot gnu.set')
        os.system('convert combined.ps combined_%s_v_%s.png' %(combination[0],combination[1],))

    return


def parse_coords_from_trajectory_file(d_charges,cwd,):

    ## vacuum permittivity
    epsilon0 = 8.854187817*10**-12 ## C2 N-1 m-2
    epsilon0 = 5.52635*10**7 ## e V-1 m-1 (division of SI units by the elementary charge)
    ## coulomb constant
    kc = 1./(4*math.pi*epsilon0) ## e-1 V m
    kc *= (10**10) ## e-1 V Angstrom
    k = 8.617343*(10**-5) ## eV/K
    eV2kT = 1./(k*308.15) ## 35'C

    d_coords = {}

    d_v3 = {
        ## v2 keys
        ## v3 values
        '1HH1':'HH11', ## ARG
        '1HH2':'HH12', ## ARG
        '2HH1':'HH21', ## ARG
        '2HH2':'HH22', ## ARG
        '1HE2':'HE21', ## GLN
        '2HE2':'HE22', ## GLN
        '1HG1':'HG11', ## VAL
        '2HG1':'HG12', ## VAL
        '3HG1':'HG13', ## VAL
        '1HG2':'HG21', ## VAL
        '2HG2':'HG22', ## VAL
        '3HG2':'HG23', ## VAL
        '1HD1':'HD11', ## LEU
        '2HD1':'HD12', ## LEU
        '3HD1':'HD13', ## LEU
        '1HD2':'HD21', ## LEU,ASN
        '2HD2':'HD22', ## LEU,ASN
        '3HD2':'HD23', ## LEU
        }

    d_charges_solvent = {
        'OW':-.834, ## TIP3
        'HW1':.417, ## TIP3
        'HW2':.417, ## TIP3
        'CL':-1,
        'Cl':-1,
        }


    for topology in ['NEUNEU','NEUCHA','CHANEU','CHACHA',]:
        if os.path.isfile('energies_%s.txt' %(topology)):
            fd = open('energies_%s.txt' %(topology),'r')
            lines = fd.readlines()
            fd.close()
            i = int(lines[-1].split()[0])+1
        else:
            i = 0

        while os.path.isfile('trjconv/2vb1_MD%i.pdb' %(i)):

            if i > 20000:
                break
            print topology, i

            d_coords35 = {}
            d_coords52 = {}
            fd = open('trjconv/2vb1_MD%i.pdb' %(i),'r')
            lines_trjconv = fd.readlines()
            fd.close()
            for line in lines_trjconv:
                record = line[:6].strip()
                if record == 'ATOM':
                    res_name1 = line[17:20]
                    if res_name1 == 'SOL':
                        break
                    res_no1 = int(line[22:26])
                    if res_no1 in [35,52,] and res_name1 != 'SOL':
                        atom_name1 = line[12:16].strip()
                        x1 = float(line[30:38])
                        y1 = float(line[38:46])
                        z1 = float(line[46:54])
##                        coord1 = numpy.array([x1,y1,z1,])
                        coord1 = [x1,y1,z1,]
                        if res_no1 == 35:
                            d_coords35[atom_name1] = coord1
                        if res_no1 == 52:
                            d_coords52[atom_name1] = coord1

            E1 = 0 ## Asp52 (side chain)
            E2 = 0 ## protein (excl. 35, incl. Asp52 main chain)
            E3 = 0 ## chloride ions
            E4 = 0 ## water
##            for atom_name1 in list(set(d_coords35.keys()+['HE2'],)-set(l_backbone_atoms)):
            for atom_name1 in list(set(d_coords35.keys()+['HE2'],)):
                print topology, i, atom_name1

                ## HE2 not present if Glu35 charged
                if atom_name1 == 'HE2' and topology[:3] == 'CHA':
                    continue

                ## Glu35 is protonated, but MD does not contain proton
                if topology[:3] == 'NEU' and cwd[:3] == 'CHA' and atom_name1 == 'HE2':
                    ## place hydrogen atom HE2 between oxygen atoms OE1,OE2
##                    coord1 = (d_coords35['OE1']+d_coords35['OE2'])/2.
                    coord1 = [
                        .5*(d_coords35['OE1'][0]+d_coords35['OE2'][0]),
                        .5*(d_coords35['OE1'][1]+d_coords35['OE2'][1]),
                        .5*(d_coords35['OE1'][2]+d_coords35['OE2'][2]),
                        ]
                else:
                    coord1 = d_coords35[atom_name1]

                charge1 = d_charges[topology][35][atom_name1]

                ##
                ## loop over atoms2
                ##
                for line in lines_trjconv:
                    
                    record = line[:6].strip()
                    if record == 'ATOM':
                        atom_no2 = int(line[6:11])
                        atom_name2 = line[12:16].strip()
                        if atom_name2 in d_v3.keys():
                            atom_name2 = d_v3[atom_name2]
                        res_name2 = line[17:20].strip()
                        res_no2 = int(line[22:26])
                        x2 = float(line[30:38])
                        y2 = float(line[38:46])
                        z2 = float(line[46:54])
##                        coord2 = numpy.array([x,y,z,])
                        coord2 = [x2,y2,z2,]

                        ## do not calculate interaction energy with self
                        if res_no2 == 35 and res_name2 != 'SOL':
                            continue

                        ## HD2 not present if Asp52 charged
                        if res_no2 == 52 and atom_name2 == 'HD2' and topology[-3:] == 'CHA':
                            continue
                        
                        if res_name2 in ['SOL','CL-','Cl',]:
                            charge2 = d_charges_solvent[atom_name2]
                        else:
        ##                    ## translate from gromacs to amber
        ##                    if atom_name == 'O1':
        ##                        atom_name = 'OC1'
        ##                    if atom_name == 'O2':
        ##                        atom_name = 'OC2'
                            try:
                                charge2 = d_charges[topology][res_no2][atom_name2]
                            except:
                                print cwd, topology, res_no2, atom_name2
                                print d_charges[topology][res_no2]
                                stop
##                        dist = math.sqrt(sum((coord2-coord1)**2))
                        dist = calc_dist(coord1,coord2,)
                        ## solvent
                        if res_name2 in ['SOL','CL-','Cl',]:
                            E = eV2kT*kc*(charge1*charge2)/(dist)
                            if res_name2 in ['CL-','Cl',]:
                                E3 += E
                            else:
                                E4 += E
                        ## protein
                        else:
                            E = eV2kT*kc*(charge1*charge2)/(dist)
                            if res_no2 == 52:
##                                if atom_name2 not in [
##                                    'N','H','CA','C','O',
##                                    ]: ## alpha carbon hydrogen...
##                                    E1 += E
##                                else:
##                                    E2 += E
                                E1 += E
                            else:
                                E2 += E

                    ## end of loop over atoms2

                ## Asp52 is protonated, but MD does not contain proton
                if topology[-3:] == 'NEU' and cwd[-3:] == 'CHA':
                    atom_name2 = 'HD2'
                    ## place hydrogen atom HD2 between oxygen atoms OD1,OD2
                    coord2 = [
                        .5*(d_coords52['OD1'][0]+d_coords52['OD2'][0]),
                        .5*(d_coords52['OD1'][1]+d_coords52['OD2'][1]),
                        .5*(d_coords52['OD1'][2]+d_coords52['OD2'][2]),
                        ]
##                    dist = math.sqrt(sum((coord2-coord1)**2))
                    dist = calc_dist(coord1,coord2,)
                    charge2 = d_charges[topology][52]['HD2']
                    E = eV2kT*kc*(charge1*charge2)/(dist)
                    E1 += E
                        
                ## end of loop over atom_names1

            ## append energies to file
            if i == 0:
                if lines[0][0] != '0':
                    line = ['%i %f %f %f %f\n' %(i,E1,E2,E3,E4,)]
                    lines = line+lines
                    fd = open('energies_%s.txt' %(topology),'w')
                    fd.writelines(lines)
                    fd.close()
                else:
                    print 'skip first frame. already calculated.', topology
            else:
                fd = open('energies_%s.txt' %(topology),'a')
                fd.write('%i %f %f %f %f\n' %(i,E1,E2,E3,E4,))
                fd.close()

            ## next frame
            i += 1

            ## end of loop over frames

        ## end of loop over topologies

    return d_coords


def calc_dist(coord1,coord2,):

    x_diff = (coord2[0]-coord1[0])
    y_diff = (coord2[1]-coord1[1])
    z_diff = (coord2[2]-coord1[2])
    dist = math.sqrt(x_diff*x_diff+y_diff*y_diff+z_diff*z_diff)

    return dist


if __name__=='__main__':
    main()

##    ##
##    ## histogram
##    ##
##    step = 1
##    d_histogram = {}
##    l = []
##    count_max = 0
##    for i in range(Min,Max+step,step):
##        d_histogram[i] = 0
##    for line in lines:
####        if int(line.split()[0]) < 50:
####            continue
####        if int(line.split()[0]) > 500:
####            break
##        Sum = float(line.split()[1])+float(line.split()[2])+float(line.split()[3])+float(line.split()[4])
##        l += [Sum]
##        d_histogram[step*int(Sum/step)] += 1
##        if d_histogram[step*int(Sum/step)] > count_max:
##            count_max = d_histogram[step*int(Sum/step)]
##    mean,sigma = statistics.stderr(l)
##    factor = count_max/(0.398942280401433/sigma)
##
##    lines2 = []
##    for k in range(Min,Max+step,step):
##        v = d_histogram[k]
##        lines2 += ['%s %s\n' %(k,v)]
##
##    fd = open('histogram.txt','w')
##    fd.writelines(lines2)
##    fd.close()
##
##    fd = open('gnu.set','w')
##    fd.writelines([
##        'set terminal png\n',
##        'set output "histogram_%s.png"\n' %(topology),
##        'set nokey\n',
##        'set title "%s"\n' %(topology),
##        'invsqrt2pi = 0.398942280401433\n',
##        'normal(x,mu,sigma)=sigma<=0?1/0:invsqrt2pi/sigma*exp(-0.5*((x-mu)/sigma)**2)\n',
##        'plot [%s:%s]"histogram.txt" u 1:2 w boxes, %s*normal(x,%s,%s)\n' %(Min,Max,factor,mean,sigma,),
####        'plot [-150:0]normal(x,%s,%s)\n' %(mean,sigma,),
##        ])
##    fd.close()
##
##    os.system('gnuplot gnu.set')

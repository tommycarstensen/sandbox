#!/software/bin/python
#$Id$
#Tommy Carstensen, University College Dublin, 2008

def main():

    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import math, gnuplot,os

    l_aa3 = [
        'ALA','CYS','ASP','GLU','PHE',
        'GLY','HIS','ILE','LYS','LEU',
        'MET','ASN','PRO','GLN','ARG',
        'SER','THR','VAL','TRP','TYR',
        ]
    l_aa1 = [
        'A','C','D','E','F',
        'G','H','I','K','L',
        'M','N','P','Q','R',
        'S','T','V','W','Y',
        ]

    ##
    ## convert data
    ##
    fd = open('single_mutant_phipsi_changes.txt','r')
    lines = fd.readlines()
    fd.close()
    lines1 = []
    lines2 = []
    lines3 = []
    lines4 = []
    lines_helix = []
    lines_helix_gly = []
    lines_helix_pro = []
    lines_helix_prepro = []
    lines_sheet = []
    lines_sheet_gly = []
    lines_sheet_pro = []
    lines_sheet_prepro = []
    lines_random = []
    for line in lines:
        if len(line.split()) < 24:
            continue
        res_name1 = line.split()[6]
        res_name2 = line.split()[7]
        if res_name1 == 'MSE':
            res_name1 = 'MET'
        if res_name2 == 'MSE':
            res_name2 = 'MET'
        if line.split()[8] == 'N/A':
            continue
        if line.split()[9] == 'N/A':
            continue
        if line.split()[10] == 'N/A':
            continue
        if line.split()[11] == 'N/A':
            continue
        phi1 = float(line.split()[8])
        psi1 = float(line.split()[9])
        phi2 = float(line.split()[10])
        psi2 = float(line.split()[11])
        ss1 = line.split()[16]
        ss2 = line.split()[17]
        res_name_next1 = line.split()[18]
        res_name_next2 = line.split()[19]
        ss_prev1 = line.split()[20]
        ss_prev2 = line.split()[21]
        ss_next1 = line.split()[22]
        ss_next2 = line.split()[23]
        if res_name_next1 != res_name_next2:
            continue ## temporary!!!
        delta_phi = phi2-phi1
        delta_psi = psi2-psi1
        if delta_phi > 180:
            delta_phi = 360-delta_phi
        if delta_psi > 180:
            delta_psi = 360-delta_psi
        if delta_phi < -180:
            delta_phi = 360+delta_phi
        if delta_psi < -180:
            delta_psi = 360+delta_psi
        delta_phipsi = math.sqrt((delta_phi)**2+(delta_psi)**2)
        hold = False
##        if ss1 == 'H' and ss2 == 'H':
##            lines_helix += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
##        elif ss1 == 'S' and ss2 == 'S':
##            lines_sheet += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
##        elif ss1 == 'R' and ss2 == 'R':
##            lines_random += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
##        else:
##            print line[:-1]
        if (
            'GLY' in [res_name1,res_name2,]
            ):
            if ss1 == 'H' and ss2 == 'H':
                if delta_phipsi > 71:
                    hold = True
                lines_helix_gly += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            elif ss1 == 'S' and ss2 == 'S':
                if delta_phipsi > 53:
                    hold = True
                lines_sheet_gly += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            elif ss1 == 'R' or ss2 == 'R':
                if delta_phipsi > 207:
                    hold = True
                lines_random += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            else:
                print ss1,ss2
                notexpected
        elif (
            'PRO' in [res_name1,res_name2,]
            ):
            if ss1 == 'H' and ss2 == 'H':
                if delta_phipsi > 88:
                    hold = True
                lines_helix_pro += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            elif ss1 == 'R' or ss2 == 'R':
                if delta_phipsi > 170:
                    hold = True
                lines_random += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            else:
                notexpected
        elif (
            'PRO' in [res_name_next1,res_name_next2,]
            ):
            if ss1 == 'H' and ss2 == 'H':
                if delta_phipsi > 27:
                    hold = True
                lines_helix_prepro += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            elif ss1 == 'R' or ss2 == 'R' or (ss1 == 'S' and ss2 == 'H'):
                if delta_phipsi > 135:
                    hold = True
                lines_random += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            else:
                print line
                print ss1,ss2
                notexpected
        else:
            if 'R' in [ss1,ss2,]:
                if delta_phipsi > 238:
                    hold = True
            elif 'R' in [ss_next1,ss_next2,ss_prev1,ss_prev2,]:
                if ss1 == 'H' and ss2 == 'H':
                    if delta_phipsi > 111: ## R/H,H,R/H
                        hold = True
                    lines_random += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
                elif ss1 == 'S' and ss2 == 'S':
                    if delta_phipsi > 70: ## R/S,S,R/S
                        hold = True
                    lines_random += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
                else:
                    notexpected
            elif 'H' in [ss1,ss2,]:
                if delta_phipsi > 59:
                    hold = True
                lines_helix += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            elif 'S' in [ss1,ss2,]:
                if delta_phipsi > 48:
                    hold = True
                lines_sheet += ['%s %s %s\n' %(l_aa3.index(res_name1),l_aa3.index(res_name2),delta_phipsi,)]
            else:
                notexpected
                

        if hold == True: ## max sqrt(180**2+180**2) = 255
            print line
            print delta_phipsi
            print delta_phi,delta_psi
            print ss_prev1,ss_prev2
            print ss1, ss2
            print ss_next1,ss_next2
            print res_name1,res_name2
            print res_name_next1,res_name_next2
            stop
##        if res_name1 == 'GLY':
##            lines2 += ['%s %s %s\n' %(l_aa.index(res_name1),l_aa.index(res_name2),delta_phipsi,)]
##        elif res_name1 == 'ALA':
##            lines3 += ['%s %s %s\n' %(l_aa.index(res_name1),l_aa.index(res_name2),delta_phipsi,)]
##        elif res_name1 == 'PRO':
##            lines4 += ['%s %s %s\n' %(l_aa.index(res_name1),l_aa.index(res_name2),delta_phipsi,)]
##        else:
##            lines1 += ['%s %s %s\n' %(l_aa.index(res_name1),l_aa.index(res_name2),delta_phipsi,)]

    d_files = {
        'helix':lines_helix,
        'sheet':lines_sheet,
        'random':lines_random,
        'helix_gly':lines_helix_gly,
        'helix_prepro':lines_helix_prepro,
        'helix_pro':lines_helix_pro,
        'sheet_gly':lines_sheet_gly,
        'sheet_prepro':lines_sheet_prepro,
        'sheet_pro':lines_sheet_pro,
        }
    for key in d_files:
        fd = open('%s.txt' %(key),'w')
        fd.writelines(d_files[key])
        fd.close()

    l_data = d_files.keys()
    for i in range(len(l_data)):
        l_data[i] = l_data[i]+'.txt'
    l_data.sort()
    print l_data

    gnuplot.plot_4d(
        l_data=l_data,
        l_xtics = l_aa1,
        x1 = 0, x2 = 19,
        y1 = 0, y2 = 19,
        xlabel = 'wtres', ylabel = 'mutres', zlabel = 'sqrt(phi^2+psi^2)',
        )

##    for key in d_files:
##        os.remove('%s.txt' %(key))

    return

if __name__ == '__main__':
    main()

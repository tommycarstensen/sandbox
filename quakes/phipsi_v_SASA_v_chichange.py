#!/software/bin/python
#$Id$
#Tommy Carstensen, University College Dublin, 2008

def main():

    import os

    ##
    ## convert data
    ##
    fd = open('single_mutant_phipsi_changes.txt','r')
    lines = fd.readlines()
    fd.close()

    lines1 = []
    lines2 = []
    for line in lines:
        if len(line.split()) < 16:
            continue
        res_name1 = line.split()[6]
        res_name2 = line.split()[7]
        if res_name1 in ['PRO','GLY','ALA',]:
            continue
        if res_name2 in ['PRO','GLY','ALA',]:
            continue
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
        r1 = line.split()[12]
        r2 = line.split()[13]
        sasa1 = float(line.split()[14])
        sasa2 = float(line.split()[15])
        phi = (phi1+phi2)/2.
        psi = (psi1+psi2)/2.
        sasa = (sasa1+sasa)/2.
        ss1 = line.split()[16]
        ss2 = line.split()[17]
        if sasa_cbcg1 == 0. and sasa_cbcg2 == 0.:
            print line
        if ss1 != 'R' and ss2 != 'R':
            if r1 == r2:
##                lines1 += ['%s %s %s\n' %(phi,psi,sasa_cbcg,)]
                lines1 += ['%s %s %s\n' %(phi1,psi1,sasa,)]
                lines1 += ['%s %s %s\n' %(phi2,psi2,sasa,)]
            else:
##                lines2 += ['%s %s %s\n' %(phi,psi,sasa_cbcg,)]
                lines2 += ['%s %s %s\n' %(phi1,psi1,sasa,)]
                lines2 += ['%s %s %s\n' %(phi2,psi2,sasa,)]

    fd = open('data1.txt','w')
    fd.writelines(lines1)
    fd.close()
    fd = open('data2.txt','w')
    fd.writelines(lines2)
    fd.close()

    settings()

##    os.remove('data1.txt')
##    os.remove('data2.txt')

def settings():

    import os

    ##
    ## write settings to list
    ##
    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 36\n',
        'set output "gnuplot.ps"\n',
        ]
    gnuplotsettings += [
        'set size square\n',
        'set size 4,4\n', ## scale size
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set xlabel "phi"\n',
        'set ylabel "psi"\n',
        'set zlabel "<SASA>"\n',
##        'set nokey\n',
        ]

    gnuplotsettings += [
        'splot [-180:180][-180:180][0:]"data1.txt" lt 2 ps 1 pt 7 t "identical gauche", "data2.txt" lt 1 ps 1 pt 7 t "different gauche"', ## splot gnuplot data file
        ]

    ##
    ## write settings to file
    ##
    fd = open('gnuplot.settings','w')
    fd.writelines(gnuplotsettings)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    gnuplot_binary = '/software/bin/gnuplot'
    os.system('%s gnuplot.settings' %(gnuplot_binary))

    ##
    ## remove data and settings
    ##
    os.remove('gnuplot.settings')
    os.remove('data1.txt')
    os.remove('data2.txt')
    os.system('convert gnuplot.ps gnuplot.png')
    os.remove('gnuplot.ps')

    return


if __name__ == '__main__':
    main()

def main():

    import os

    lines2 = []

    d_residues = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
        'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
        'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
        'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
        }

    d_max = {}

    import sys
    sys.path.append('/home/people/tc/svn/tc_sandbox/quakes/')
    import phipsi_comparison
    sys.path.append('/home/people/tc/svn/tc_sandbox/misc/')
    import phipsi_plot

    d_ramachandran = {}
    for res in d_residues.values():
        d_ramachandran[res] = {}
        for phi in range(-180,180):
            d_ramachandran[res][phi] = {}
            for psi in range(-180,180):
    ##            d_ramachandran[res][phi][psi] = 0
    ##d_ramachandran = phipsi_comparison.parse_dihedrals(d_residues) ## slow
                d_ramachandran[res][phi][psi] = 1

    for res in d_ramachandran.keys():
        max = 0
        for phi in d_ramachandran[res].keys():
            for psi in d_ramachandran[res][phi].keys():
                count = d_ramachandran[res][phi][psi]
                if count > max:
                    max = count
        d_max[res] = max

    option = 'pdbs'

    if option == 'frames':
        l = range(50)
    elif option == 'pdbs':
        l = ['150l',]
        chain = 'D'
    elif option == 'gromacs':
        l = []
        for i in range(101):
            l += ['2LZM_wt_trjconv_%i' %(i)]

    for pdb in l:

        if option == 'gromacs':
            i = int(pdb[pdb.rindex('_')+1:])

        lines2 = []

        if option in ['frames','gromacs',]:
            lines2 += ['HEADER    frame t=%2i.000                                                       \n' %(i)]
            lines2 += ['MODEL       %2i                                                                 \n' %(i)]

            prefix = 'rasmol%s' %(i)
            ## write rasmol script
            lines_rasmol = [
                'rasmol -nodisplay frame%s.pdb << EOF\n' %(str(i).zfill(2)),
                'cartoon\n',
                'wireframe 0\n',
##                    'ribbons\n',
##                    'color temperature\n',
                ]
     
        d_phipsi, lines_gnuplot = phipsi_plot.calculate_phipsi(pdb, chain = chain)

        fd = open('%s.pdb' %(pdb))
        lines1 = fd.readlines()
        fd.close()

        for line1 in lines1:

            record = line1[:6].strip()

            if record != 'ATOM':
                continue

            atom_name = line1[12:16].strip()
    ##        if atom_name != 'CA':
    ##            continue
            if atom_name not in ['N','CA','C','O',]:
                continue

            bfactor = float(line1[60:66])
            res_no = int(line1[22:26])
            res_name = line1[17:20].strip()
            if res_no not in d_phipsi.keys():
                bfactor = 50.
                continue
            else:
                phi = int(d_phipsi[res_no][0])
                psi = int(d_phipsi[res_no][1])
                if phi == 180.:
                    phi = -180.
                if psi == 180.:
                    psi = -180.
                bfactor = (
                    100
                    *d_ramachandran[d_residues[res_name]][phi][psi]
                    /float(d_max[d_residues[res_name]])
                    )

            lines2 += ['%s%6.2f%s' %(line1[:60],bfactor,line1[66:])]

            if option in ['frames','gromacs',]:
                lines2 += ['TER                                                                             \n']
                lines2 += ['ENDMDL                                                                          \n']

            if option == 'pdbs':
                fd = open('%s_phipsi_colored.pdb' %(pdb),'w')
                fd.writelines(lines2)
                fd.close()

            if option in ['frames','gromacs',]:

                h = 159.+80*bfactor/100.
                s = 240.
                l = 120.+120.*(50-abs(bfactor-50))/50.
                r,g,b = hsl2rgb(h,s,l,)
                lines_rasmol += [
                    'select %i\n' %(res_no),
                    'color [%i,%i,%i]\n' %(
                        r, g, b,
                        )
                    ]

        if option in ['frames','gromacs',]:
            fd = open('%s.pdb' %(pdb),'w')
            fd.writelines(lines2)
            fd.close()

            lines_rasmol += [
##                    'ribbons\n',
                'rotate x 100\n',
                'rotate z 30\n',
                'rotate x 100\n',
                'rotate z 90\n',
                'rotate x 40\n',
                'rotate y -20\n',
                'write %s.ppm\n' %(prefix),
                'exit\n',
                ]
            ## write rasmol script to file
            fd = open('%srasmol.src' %(prefix),'w')
            fd.writelines(lines_rasmol)
            fd.close()
            ## execute rasmol script
            os.system('source %srasmol.src > %srasmol.log' %(prefix,prefix))
            ## convert rasmol output
##                os.system('convert %s.ppm -resize x80 %s.gif' %(prefix,prefix))
            os.system('convert %s.ppm %s.gif' %(prefix,prefix))
            ## clean up
            os.remove('%s.ppm' %(prefix))
            os.remove('%srasmol.log' %(prefix))
            os.remove('%srasmol.src' %(prefix))
##                os.remove('frame%s.pdb' %(str(i).zfill(2)))


    if option in ['frames','gromacs',]:

        line = 'convert '
        for frame in range(50):
            line += 'rasmol%s.gif ' %(i)
        for frame in range(50-1-1,-1+1,-1):
            line += 'rasmol%s.gif ' %(i)
        line += '-loop 0 mode%s.gif' %(mode)
        print line
        os.system(line)
        for frame in range(50):
            os.remove('rasmol%s.gif' %(i))

    return


def hsl2rgb(h,s,l):

    h = h/240.
    s = s/240.
    l = l/240.

    if l < .5:
        q = l+l*s
    else:
        q = l+s-l*s
    p = 2*l-q

    tr = h+1./3
    tg = h
    tb = h-1./3

    l_rgb = [0,0,0]
    l_tc = [tr,tg,tb]
    for i in range(len(l_tc)):
        tc = l_tc[i]
        if tc < 0:
            tc += 1
        if tc > 1:
            tc -= 1

        if tc < 1./6:
            c = p+(q-p)*6*tc
        elif 1./6 <= tc and tc < .5:
            c = q
        elif .5 <= tc and tc < 2./3.:
            c = p+(q-p)*(2./3-tc)*6
        else:
            c = p

        l_rgb[i] = c

    r = int(l_rgb[0]*255)
    g = int(l_rgb[1]*255)
    b = int(l_rgb[2]*255)

    return r,g,b

if __name__=='__main__':

    main()

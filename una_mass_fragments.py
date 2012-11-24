#!/usr/bin/python
def main():

    ## this script upon calculation of fragment masses assumes no posttranslational modifications, no isotope labeling (natural isotope abundance) and no disulfide bonds

    max_delta_mass = 0.5

    ## input - sequence as a string
    sequence = '''
GNAAAAKKG  SEQESVKEFL AKAKEDFLKK WETPSQNTAQ LDQFDRIKTL
GTGSFGRVML VKHKESGNHY AMKILDKQKV VKLKQIEHTL NEKRILQAVN FPFLVKLEFS 
FKDNSNLYMV MEYVAGGEMF SHLRRIGRFS EPHARFYAAQ IVLTFEYLHS LDLIYRDLKP 
ENLLIDQQGY IQVTDFGFAK RVKGRTWTLC GTPEYLAPEI ILSKGYNKAV DWWALGVLIY 
EMAAGYPPFF ADQPIQIYEK IVSGKVRFPS HFSSDLKDLL RNLLQVDLTK RFGNLKNGVN 
DIKNHKWFAT TDWIAIYQRK VEAPFIPKFK GPGDTSNFDD YEEEEIRVSI NEKCGKEFTE
F
'''

    ## input - space separated masses
    s_m_f = '''
2756.0186
2760.0669
'''

    ## convert string to list
    l_m_f = s_m_f.split()

    ## convert list to dictionary
    d_m_f = {}
    for mass in l_m_f:
        d_m_f[float(mass)] = []
    
    m_max = max(d_m_f.keys())+10; m_min = min(d_m_f.keys())-10

    sequence = sequence.replace(' ','').replace('\n','')

    ## calculate mass of natural abundance isotopes
    d_m_a = {
        'H':(99.985*1.00782504+0.015*2.01410178)/100.,
        'C':(98.89*12+1.11*13.00335483)/100.,
        'N':(99.63*14.003074+0.366*15.00010897)/100.,
        'O':(99.76*15.99491463+0.038*16.9991312+0.202*17.9991603)/100.,
        'S':(95.02*31.9720707+0.75*32.9714584+4.21*33.9678667+0.017*35.9670806)/100.,
        }
    H = d_m_a['H']
    C = d_m_a['C']
    N = d_m_a['N']
    O = d_m_a['O']
    S = d_m_a['S']

    ## calculate mass of residues
    d_m_r = {
        'A': 7*H+ 3*C+1*N+2*O+0*S-2*H-1*O,
        'C': 7*H+ 3*C+1*N+2*O+1*S-2*H-1*O,
        'D': 7*H+ 4*C+1*N+4*O+0*S-2*H-1*O,
        'E': 9*H+ 5*C+1*N+4*O+0*S-2*H-1*O,
        'F':11*H+ 9*C+1*N+2*O+0*S-2*H-1*O,
        'G': 5*H+ 2*C+1*N+2*O+0*S-2*H-1*O,
        'H': 9*H+ 6*C+3*N+2*O+0*S-2*H-1*O,
        'I':13*H+ 6*C+1*N+2*O+0*S-2*H-1*O,
        'K':14*H+ 6*C+2*N+2*O+0*S-2*H-1*O,
        'L':13*H+ 6*C+1*N+2*O+0*S-2*H-1*O,
        'M':11*H+ 5*C+1*N+2*O+1*S-2*H-1*O,
        'N': 8*H+ 4*C+2*N+3*O+0*S-2*H-1*O,
        'P': 9*H+ 5*C+1*N+2*O+0*S-2*H-1*O,
        'Q':10*H+ 5*C+2*N+3*O+0*S-2*H-1*O,
        'R':14*H+ 6*C+4*N+2*O+0*S-2*H-1*O,
        'S': 7*H+ 3*C+1*N+3*O+0*S-2*H-1*O,
        'T': 9*H+ 4*C+1*N+3*O+0*S-2*H-1*O,
        'V':11*H+ 5*C+1*N+2*O+0*S-2*H-1*O,
        'W':12*H+11*C+2*N+2*O+0*S-2*H-1*O,
        'Y':11*H+ 9*C+1*N+3*O+0*S-2*H-1*O,
        }

    ## calculate mass of fragments and add to dictionary
    d_m_f_calc = {}
    for res1 in range(len(sequence)):
        print res1
        mass = 2*H+1*O
        for res2 in range(res1+1,len(sequence)):
            res = sequence[res2-1]
            mass += d_m_r[res]
            if mass >= m_min and mass <= m_max:
                if mass not in d_m_f_calc.keys():
                    d_m_f_calc[mass] = [[res1,res2]]
                else:
                    d_m_f_calc[mass] += [[res1,res2]]

    ## identify equalities among input and output masses
    l_m_f_calc = d_m_f_calc.keys()
    l_m_f_calc.sort()
    l_m_f = d_m_f.keys()
    l_m_f.sort()
    for m in l_m_f:
        print
        for m_calc in l_m_f_calc:
            if m_calc > m-max_delta_mass and m_calc < m+max_delta_mass:
                for res in d_m_f_calc[m_calc]:
                    print round(m_calc,1), '%3i' %(res[0]+1), '%3i' %(res[1]), '%s/%s/%s' %(sequence[res[0]-1], sequence[res[0]:res[1]], sequence[res[1]])
                d_m_f[m] += d_m_f_calc[m_calc]

    return d_m_f

if __name__ == '__main__':
    main()

def main():

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
##                'A 2',
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

    d_crystalsystems = {}
    for crystalsystem in d_spacegroups.keys():
        for spacegroup in d_spacegroups[crystalsystem]:
            d_crystalsystems[spacegroup] = crystalsystem

    ## list of biounit strings in which other biounit string can be found
    l_biounits = [
        'DIMER OF DIMERS', ## DIMER
        'HEXADECAMER', ## DECAMER
        'DODECAMER', ## DECAMER
        ]
    d_biounits = {
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

    return d_spacegroups, d_crystalsystems

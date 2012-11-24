import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import parse_mmCIF

lines = []
for s in '01234567890abcdefghijklmnopqrstuvwxyz':
    print s
    fd = open('single_point_mutations/%s.txt' %(s),'r')
    lines += fd.readlines()
    fd.close()

for i_line in range(len(lines)):
    if i_line % 100 == 0:
        d_coordinates = {}
    line = lines[i_line]
    l = line.split()
    pdb1 = l[0]
    pdb2 = l[1]
    chain1 = l[4]
    chain2 = l[5]

    for pdb,chain in [[pdb1,chain1,],[pdb2,chain2,],]:

        if pdb in d_coordinates.keys():
            continue

        d_mmCIF = parse_mmCIF.main(pdb)


        if d_mmCIF['_pdbx_poly_seq_scheme.pdb_seq_num'] != d_mmCIF['_pdbx_poly_seq_scheme.author_seq_num']:
            print d_mmCIF['_pdbx_poly_seq_scheme.pdb_seq_num']
            print d_mmCIF['_pdbx_poly_seq_scheme.author_seq_num']
            stop

        d_coords = {}
        d_ndb_seq_num = {}
        for i_seq in range(len(d_mmCIF['_pdbx_poly_seq_scheme.ndb_seq_num'])):
            if d_mmCIF['_pdbx_poly_seq_scheme.pdb_strand_id'][i_seq] != chain:
                continue
            ndb_seq_num = d_mmCIF['_pdbx_poly_seq_scheme.ndb_seq_num'][i_seq]
            pdb_seq_num = d_mmCIF['_pdbx_poly_seq_scheme.pdb_seq_num'][i_seq]
            d_ndb_seq_num[pdb_seq_num] = ndb_seq_num
            d_coords[ndb_seq_num] = {}
            
        for i_atom in range(len(d_mmCIF['_atom_site.label_atom_id'])):

            if d_mmCIF['_atom_site.label_atom_id'][i] != 'CA':
                continue

            ## chain
            if d_mmCIF['_atom_site.auth_asym_id'][i] != pdb[-1].upper():
                continue
            ## skip non-polymer atoms
            if d_mmCIF['_atom_site.group_PDB'][i] == 'HETATM': ## not used for modres...
                continue
            ## altloc
            if d_mmCIF['_atom_site.label_alt_id'][i] not in ['.','A',]:
                continue

            res_no = pdb_seq_num = int(d_mmCIF['_atom_site.label_seq_id'][i])
            ndb_seq_num = d_ndb_seq_num[pdb_seq_num]
            x = d_mmCIF['_atom_site.Cartn_x'][i_atom]
            y = d_mmCIF['_atom_site.Cartn_y'][i_atom]
            z = d_mmCIF['_atom_site.Cartn_z'][i_atom]
            coord = numpy.array([x,y,z,])
            tempfact = d_mmCIF['_atom_site.B_iso_or_equiv'][i_atom]
            atom_id = d_mmCIF['label_atom_id'][i_atom]
            d_coords[ndb_seq_num][atom_id] = [coord,tempfact,]

        l_coordinates = []
        for res_no in range(1,1+len_seq):
            if res_no in d_coords.keys():
                l_coordinates += [d_coords[res_no],]
            else:
                l_coordinates += [[None,None,]]

        d_coordinates[pdb] = l_coordinates

    if len(d_coordinates[pdb1]) != len(d_coordinates[pdb2]):
        print
        print pdb1, len(d_coordinates[pdb1])
        print pdb2, len(d_coordinates[pdb2])
        stop


d_uniprot = {}
for i in range(len(lines)):
    line = lines[i]
    l = line.split()
    pdb1 = l[0]
    pdb2 = l[1]
    chain1 = l[4]
    chain2 = l[5]

    for pdb,chain in [[pdb1,chain1,],[pdb2,chain2,],]:

        if pdb in d_uniprot.keys():
            continue

        print pdb

        d = parse_mmCIF.main(pdb)

        s_uniprot = d['_struct_ref_seq.pdbx_db_accession'][d['_struct_ref_seq.pdbx_strand_id'].index(chain)]

        d_uniprot[pdb] = s_uniprot

##d_uniprot_reverse = {}
##for k,v in d_uniprot.items():
##    if not v in d_uniprot_reverse.keys():
##        d_uniprot_reverse[v] = []
##    d_uniprot_reverse[v] += [k]

##l_uniprot = []
##for k,v in d_uniprot_reverse.items():
##    l_uniprot += [[len(v),k,]]
##l_uniprot.sort()

l_uniprot = [[1, '11125386'], [1, '19005'], [1, '1DQJ'], [1, '1GWO'], [1, '1GWT'], [1, '1GWU'], [1, '1GX2'], [1, '1GYD'], [1, '1GYH'], [1, '1T2X'], [1, '1XGP'], [1, '1XGQ'], [1, '1XGT'], [1, '1XGU'], [1, '208092'], [1, '2BXW'], [1, '2DQC'], [1, '2DQD'], [1, '2DQE'], [1, '2DQF'], [1, '2DQG'], [1, '2DQH'], [1, '2DQI'], [1, '2DQJ'], [1, '2J45'], [1, '2J46'], [1, '2JI0'], [1, '2VLO'], [1, '2W2N'], [1, '2WKP'], [1, '2WKQ'], [1, '2WKR'], [1, '3561059'], [1, '55771634'], [1, 'O07347'], [1, 'P00274'], [1, 'P00647'], [1, 'P02675'], [1, 'P07062'], [1, 'P62974'], [1, 'P62990'], [1, 'P80380'], [1, 'Q29846'], [1, 'Q6ISD4'], [1, 'Q99QV3'], [2, '545862'], [2, '896294'], [2, '902938'], [2, 'A1L467'], [2, 'B3EYM8'], [2, 'O14965'], [2, 'O57883'], [2, 'P00651'], [2, 'P00711'], [2, 'P00974'], [2, 'P01112'], [2, 'P01899'], [2, 'P02679'], [2, 'P02754'], [2, 'P02788'], [2, 'P04637'], [2, 'P05082'], [2, 'P06132'], [2, 'P06875'], [2, 'P08037'], [2, 'P08799'], [2, 'P09012'], [2, 'P09960'], [2, 'P11350'], [2, 'P15207'], [2, 'P15917'], [2, 'P19491'], [2, 'P21397'], [2, 'P28012'], [2, 'P60010'], [2, 'Q42795'], [2, 'Q746K1'], [2, 'Q7SID2'], [2, 'Q7SIG4'], [2, 'Q90092'], [2, 'Q9A5I0'], [2, 'Q9RVU2'], [2, 'Q9S0S3'], [3, '15928840'], [3, 'O15527'], [3, 'P00915'], [3, 'P01053'], [3, 'P01130'], [3, 'P01635'], [3, 'P01642'], [3, 'P01730'], [3, 'P02618'], [3, 'P02928'], [3, 'P03989'], [3, 'P04547'], [3, 'P04745'], [3, 'P0A6M2'], [3, 'P0AE18'], [3, 'P11309'], [3, 'P11509'], [3, 'P14555'], [3, 'P29373'], [3, 'P30685'], [3, 'P33247'], [3, 'P35804'], [3, 'P37957'], [3, 'P46925'], [3, 'Q51912'], [3, 'Q70V27'], [3, 'Q93LD7'], [3, 'Q9Y233'], [4, 'A2NAI0'], [4, 'O68601'], [4, 'P01823'], [4, 'P08174'], [4, 'P08235'], [4, 'P10538'], [4, 'P11439'], [4, 'P12295'], [4, 'P12497'], [4, 'P17707'], [4, 'P24991'], [4, 'P26361'], [4, 'P26918'], [4, 'P28907'], [4, 'P29679'], [4, 'P49356'], [4, 'Q03347'], [4, 'Q24451'], [4, 'Q55080'], [4, 'Q5SLP1'], [4, 'Q8NBP7'], [4, 'Q9U9J6'], [5, '1333979'], [5, 'P02883'], [5, 'P05102'], [5, 'P07906'], [5, 'P0ABQ4'], [5, 'P14489'], [5, 'P23222'], [5, 'P36924'], [5, 'Q46822'], [5, 'Q9AJ26'], [5, 'Q9BY41'], [6, 'P04190'], [6, 'P11124'], [6, 'P11540'], [6, 'P62988'], [6, 'Q01745'], [6, 'Q02762'], [7, 'P00374'], [7, 'P00695'], [7, 'P07061'], [7, 'P09850'], [7, 'P0AA25'], [7, 'P10275'], [7, 'P11086'], [7, 'P15879'], [7, 'P40859'], [7, 'Q40059'], [8, 'P00782'], [8, 'P13479'], [8, 'P26281'], [8, 'Q72497'], [8, 'Q8U094'], [9, 'P00171'], [9, 'P00431'], [9, 'P00433'], [9, 'P00763'], [9, 'P62942'], [10, 'P00282'], [10, 'P06886'], [10, 'P28720'], [10, 'P58502'], [10, 'P62593'], [10, 'Q2M889'], [11, 'P00044'], [11, 'P00323'], [11, 'P05798'], [11, 'P14779'], [12, 'P68390'], [13, 'P02189'], [13, 'P09601'], [13, 'P19120'], [13, 'Q16539'], [13, 'Q16769'], [14, 'P00644'], [14, 'P10824'], [14, 'P14769'], [14, 'P22364'], [15, 'P00268'], [15, 'P00533'], [15, 'P02787'], [17, 'P19614'], [17, 'P68082'], [17, 'Q53WB3'], [18, 'P00593'], [18, 'P04746'], [18, 'P16442'], [19, 'P00778'], [20, 'P00772'], [21, 'P21890'], [22, 'P00590'], [22, 'P06873'], [24, 'P00214'], [27, 'P0A7Y4'], [28, 'Q59560'], [38, 'P02185'], [39, 'P16113'], [79, 'P61823'], [171, 'P00698'], [171, 'P61626'], [356, 'P00720']]

## 2nwd (human, but chemical synthesis)

## 1963 pdbs compared, 15563 comparisons
## 356 T4 lysozyme
## 171 HEWL
## 171 human lysozyme
## 79 bovine RNase
## 39 halophile photoactive yellow protein
## 38 sperm whale myoglobin
## 28 mycobacterium recombinase
## 27 coli RNase
## 24 azotobacter ferredoxin
## 22 tritirachium alkaline protease
## 22 fuasrium cutinase
## 986 other

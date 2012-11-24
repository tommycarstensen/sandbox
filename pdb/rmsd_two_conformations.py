#!/software/bin/python
#
#$Id$
#
#Tommy Carstensen, University College Dublin, 2008

def main():

    import os
    import sys
    sys.path.append('/home/people/tc/svn/Protool/trunk')
    import geometry
    instance_geometry = geometry.geometry()

##    fd = open('clusters95.txt','r')
##    lines = fd.readlines()
##    fd.close()
##    d_pdbs = {}
##    d_clusters = {}
##    for line in lines:
##        cluster = int(line.split()[0])
##        if cluster not in d_clusters.keys():
##            d_clusters[cluster] = []
##        pdb = line.split()[2]
##        d_clusters[cluster] += [pdb]
##    for cluster in d_clusters:
##        l_pdbs = d_clusters[cluster]
##        if '2LZM:A' in l_pdbs:
##            break

    pdb1 = '2lzt'

##    l_pdbs = ['189L:A', '1C69:A', '1C6A:A', '1C6B:A', '1D9W:A', '1DYA:A', '1DYB:A', '1DYC:A', '1DYD:A', '1DYE:A', '1DYF:A', '1DYG:A', '1L00:A', '1L02:A', '1L03:A', '1L04:A', '1L05:A', '1L06:A', '1L07:A', '1L08:A', '1L09:A', '1L10:A', '1L11:A', '1L12:A', '1L13:A', '1L14:A', '1L15:A', '1L16:A', '1L17:A', '1L18:A', '1L19:A', '1L20:A', '1L21:A', '1L22:A', '1L23:A', '1L24:A', '1L25:A', '1L26:A', '1L27:A', '1L28:A', '1L29:A', '1L30:A', '1L31:A', '1L32:A', '1L33:A', '1L34:A', '1L37:A', '1L38:A', '1L42:A', '1L43:A', '1L44:A', '1L45:A', '1L46:A', '1L47:A', '1L48:A', '1L52:A', '1L53:A', '1L56:A', '1L57:A', '1L58:A', '1L60:A', '1L69:A', '1L70:A', '1L71:A', '1L96:A', '1L97:A', '1L97:B', '1L98:A', '1L99:A', '1LYD:A', '1T6H:A', '223L:A', '225L:A', '226L:A', '256L:A', '2LZM:A', '3LZM:A', '4LZM:A', '5LZM:A', '6LZM:A', '7LZM:A']
##    l_pdbs.remove('%s:A' %(pdb1.upper()))
##    l_pdbs.sort()

    l_pdbs = ['1sfg:A','1ja7:A','1h6m:A',]

    l_pdbs2 = [
        ## antibody
##        '2eiz','2eks','2yss','1a2y','1bql','1bvk','1c08','1fdl','1g7h','1g7i',
##        '1g7j','1g7l','1g7m',
        '1mel','1nbz','1yqv','1zvh',
        ## MODRES
        '132l','1at5','1at6',
        ## inhibitor
        '1gpq','1uuz',
        ## deletion
        '1uia','1uib',
        ## glycosylated
        '2b5z',
        ## NMR
        '1e8l',
        ## multiple models
        '1hc0','2d6b',
        ## different length
        '1lsg',
        ]

    l_pdbs_2lzm = [
        ## modified residue
        '1t6h',
        ## beta-mercapto-ethanol
        '1l97',
        ]

    d_coordinates = {}
    d_coordinates[pdb1],l_coordinates1 = parse_coordinates(pdb1,'A')

    max_rmsd = [0,0]
    for i in range(len(l_pdbs)):
        pdb = l_pdbs[i][:4].lower()
        chain = l_pdbs[i][-1]
        if pdb in l_pdbs_2lzm:
            continue
##        if i+1 < int(sys.argv[1]):
##            continue
        d_coordinates[pdb],l_coordinates2 = parse_coordinates(pdb,chain)
        if len(l_coordinates1) == 164 and len(l_coordinates2) in [162,163,]:
            rmsd = instance_geometry.superpose(l_coordinates1[:len(l_coordinates2)],l_coordinates2)
        elif len(l_coordinates1) != len(l_coordinates2):
            print pdb, chain, len(l_coordinates1), len(l_coordinates2)
            stop
        else:
            rmsd = instance_geometry.superpose(l_coordinates1,l_coordinates2)
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter
        print pdb, chain, i+1, len(l_pdbs), rmsd
        if rmsd > max_rmsd[1]:
            max_rmsd = [pdb,rmsd]

    print max_rmsd
    
    return

def parse_coordinates(pdb,chain):

    fd = open('/local/data/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb),'r')
    lines = fd.readlines()
    fd.close()

    d_coordinates = {}
    l_coordinates = []
    l_modres = []
    for line in lines:
        record = line[:6].strip()
        if record == 'MODRES':
            hetID = line[12:15].strip()
            l_modres += [hetID]
        elif record == 'ATOM':
            d_coordinates,l_coordinates = parse_recordATOM(line,d_coordinates,l_coordinates,chain)
        elif record == 'HETATM':
            res_name = line[17:20].strip()
            if res_name in l_modres:
                d_coordinates,l_coordinates = parse_recordATOM(line,d_coordinates,l_coordinates,chain)
            
    return d_coordinates,l_coordinates

def parse_recordATOM(line,d_coordinates,l_coordinates,CHAIN):

    import Numeric

    atom_name = line[12:16].strip()
    altloc = line[16]
    res_name = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    coordinate = Numeric.array([x, y, z])
    if iCode != ' ':
        print line
        stop

    if not chain in d_coordinates.keys():
        d_coordinates[chain] = {}
    if not res_no in d_coordinates[chain].keys():
        d_coordinates[chain][res_no] = {}
    if not altloc in d_coordinates[chain][res_no].keys():
        d_coordinates[chain][res_no][altloc] = {}
    if not atom_name in d_coordinates[chain][res_no][altloc].keys():
        d_coordinates[chain][res_no][altloc][atom_name] = coordinate

    if chain == CHAIN and atom_name == 'CA' and altloc in [' ','A',]:
        l_coordinates += [coordinate]

    return d_coordinates,l_coordinates

if __name__ == '__main__':
    main()

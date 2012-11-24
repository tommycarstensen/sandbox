## alpha carbon atom alignment

def main(pdb1,pdb2,chains1_align,chains2_align,):

    chains1_apply = chains1_align
    chains2_apply = chains2_align

    import os, sys, math
    sys.path.append('/home/people/tc/svn/Protool/')
    import geometry
    instance_geometry = geometry.geometry()

    domain_range = range(0,9999)

    os.system('cp /data/pdb-v3.2/%s/pdb%s.ent %s.pdb' %(pdb1[1:3],pdb1,pdb1,))
    os.system('cp /data/pdb-v3.2/%s/pdb%s.ent %s.pdb' %(pdb2[1:3],pdb2,pdb2,))
    
    ss_range1, l_missing1, seqres1, l_modres = parse_header(pdb1, chains1_align,)
    ss_range2, l_missing2, seqres2, l_modres = parse_header(pdb2, chains2_align,)

    ss_range = list(set(ss_range1)&set(ss_range2))
    l_missing = list(set(l_missing1)|set(l_missing2))

    if len(seqres1) != len(seqres2):
        d_replace = {
            'TPO':'THR','PTR':'TYR',
##                'SER':'CYS', ## 1tde v 1f6m
            }
        for i in range(len(seqres1)):
            if seqres1[i] in d_replace.keys():
                seqres1[i] = d_replace[seqres1[i]]
        for i in range(len(seqres2)):
            if seqres2[i] in d_replace.keys():
                seqres2[i] = d_replace[seqres2[i]]
        if not (''.join(seqres1) in ''.join(seqres2) or ''.join(seqres2) in ''.join(seqres1)):
            import sys
            sys.path.append('/home/people/tc/svn/EAT_DB/')
            import sequence_alignment
            d_res = {
                'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
                'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
                'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
                'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
                }
            seq1 = ''
            for res in seqres1:
                seq1 += d_res[res]
            seq2 = ''
            for res in seqres2:
                seq2 += d_res[res]
            instance = sequence_alignment.NW(seq1,seq2)
            s1,s2 = instance.Align(verbose=False)[:2]
            l1 = len(s1)-len(s1.lstrip('-'))
            l2 = len(s2)-len(s2.lstrip('-'))
            r1 = len(s1)-len(s1.rstrip('-'))
            r2 = len(s2)-len(s2.rstrip('-'))
            print seqres1
            print seqres2
            print len(seqres1)
            print len(seqres2)
            print pdb1, pdb2
            print l1, l2, r1, r2
            print seqres2[l1:len(seqres2)-r1]
            print seqres1[l2:len(seqres1)-r2]
        else:
            s1 = ''.join(seqres1)
            s2 = ''.join(seqres2)
            if s1 in s2:
                seqres2 = seqres2[s2.index(s1)/3:]
            else:
                seqres1 = seqres1[s1.index(s2)/3:]
            if len(seqres1) != len(seqres2):
                print len(seqres1), len(seqres2)
                stop

    l_coordinates1 = parse_coordinates(pdb1,chains1_align,domain_range,ss_range,l_missing,)
    l_coordinates2 = parse_coordinates(pdb2,chains2_align,domain_range,ss_range,l_missing,)

##        l_coordinates1 = l_coordinates1[l2:len(l_coordinates1)-r2]
##        l_coordinates2 = l_coordinates2[l1:len(l_coordinates2)-r1]

    if len(l_coordinates1) == 0 or len(l_coordinates2) == 0:
        stop

    if len(l_coordinates1) != len(l_coordinates2):
        print len(l_coordinates1)
        print len(l_coordinates2)
        stop

    rmsd = instance_geometry.superpose(l_coordinates1,l_coordinates2)
    print pdb1, pdb2
    print 'rmsd', round(rmsd,1)
    print 'residues', len(seqres1), len(seqres2)
    print 'coordinates', len(l_coordinates1)
    tv1 = instance_geometry.fitcenter
    rm = instance_geometry.rotation
    tv2 = instance_geometry.refcenter

    lines1 = apply_transformation_matrix(
        pdb1,chains1_apply,l_modres,
        [0,0,0],[[1,0,0],[0,1,0],[0,0,1]],[0,0,0],
        )
    lines2 = apply_transformation_matrix(
        pdb2,chains2_apply,l_modres,
        tv1,rm,tv2,
        )

    fd = open('rotated_%s%s.pdb' %(pdb1,pdb2,),'w')
    fd.writelines(lines1+lines2)
    fd.close()

    l_coordinates1 = parse_coordinates('rotated_'+pdb1,chains1_apply,range(-9999,9999),ss_range,l_missing,)
    l_coordinates2 = parse_coordinates('rotated_'+pdb2,chains2_apply,range(-9999,9999),ss_range,l_missing,)

    SUM = 0.
    n = len(l_coordinates1)
    for i in range(n):
        SUM += sum((l_coordinates1[i]-l_coordinates2[i])**2)
    RMSD = math.sqrt(SUM/n)
    print 'RMSD all atoms', RMSD

    return RMSD,l_coordinates1,l_coordinates2


def apply_transformation_matrix(pdb,chains,l_modres,tv1,rm,tv2,):

    import numpy, os

    fd = open('%s.pdb' %(pdb),'r')
    lines = fd.readlines()
    fd.close()

    lines2 = []
    for line in lines:
        record = line[:6].strip()
        ## not coordinate section
        if record not in ['ATOM','HETATM',]:
            lines2 += [line]
            continue
        if line[21] not in chains:
            if record == 'ATOM':
                continue
            if record == 'HETATM':
                if line[17:20] == 'HOH':
                    continue
##                if line[21:26] not in l_modres:
##                    continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coordinate = numpy.array([x,y,z,])
        coordinate2 = numpy.dot(coordinate-tv1,rm)+tv2
        x2 = coordinate2[0]
        y2 = coordinate2[1]
        z2 = coordinate2[2]
        line2 = '%s%8.3f%8.3f%8.3f%s' %(line[:30],x2,y2,z2,line[54:],)
        lines2 += [line2]

    fd = open('rotated_%s.pdb' %(pdb),'w')
    fd.writelines(lines2)
    fd.close()

    return lines2


def parse_header(pdb, chains,):

    fd = open('%s.pdb' %(pdb),'r')
    lines = fd.readlines()
    fd.close()

    ss_range = []
    l_missing = []
    seqres = []
    l_modres = []
    for i in range(len(lines)):
        line = lines[i]
        record = line[:6].strip()
        if record == 'REMARK':
            remark = int(line[6:10])
            if remark == 465:
                if line.strip() == 'REMARK 465   M RES C SSSEQI':
                    for j in range(i+1,len(lines)):
                        if lines[j][:10] != 'REMARK 465':
                            break
                        res_no = int(lines[j][22:26])
                        l_missing += [res_no]
        elif record == 'SEQRES':
            if line[11] in chains:
                seqres += line[19:70].split()
        elif record == 'HELIX':
            res_no1 = int(line[21:25])
            res_no2 = int(line[33:37])
            ss_range += range(res_no1,res_no2+1)
        elif record == 'SHEET':
            res_no1 = int(line[22:26])
            res_no2 = int(line[33:37])
            ss_range += range(res_no1,res_no2+1)
        elif record == 'MODRES':
            chain = line[16]
            res_no = line[18:22]
            iCode = line[22]
            l_modres += [chain+res_no+iCode]
                        
    return ss_range, l_missing, seqres, l_modres


def parse_coordinates(pdb,chains,domain_range,ss_range,l_missing,):

    import numpy

    fd = open('%s.pdb' %(pdb),'r')
    lines = fd.readlines()
    fd.close()

    l_coordinates = []
    for line in lines:
        record = line[:6].strip()
        if record == 'MODEL':
            stop
        if record != 'ATOM':
            continue
        atom_name = line[12:16].strip()
        if atom_name not in ['CA',]: ## alpha carbon only at the moment!!!
            continue
        res_no = int(line[22:26])
        chain = line[21]
        ## exclude flexible domain
        if res_no not in domain_range:
            continue
##        ## exclude random coil
##        if res_no not in ss_range:
##            continue
        ## exclude missing residues
        if res_no in l_missing:
            continue
        ## exclude other chains than the one of interest
        if chain not in chains:
            continue
        ## exlude altlocs
        if line[16] not in [' ','A',]:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        l_coordinates += [numpy.array([x,y,z,])]

    return l_coordinates


if __name__ == '__main__':
    pdb1 = '3lyt'
    pdb2 = '3lyt'
    chains1_align = ['A',]
    chains2_align = ['B',]
    chains1_apply = ['A',]
    chains2_apply = ['B',]
    main()

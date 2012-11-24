def main():

    import os, sys, math
    sys.path.append('/home/people/tc/svn/Protool/')
    import geometry
    instance_geometry = geometry.geometry()

##    ## lactoferrin
##    pdb1 = '1lfg'
##    pdb2 = '1lfh'
##    domain_range = range(1,88+1)+range(253,333+1)
##    chain = 'A'

##    ## trp repressor
##    pdb1 = '1wrp'
##    pdb2 = '2oz9' ## 2wrp
##    pdb2 = '1zt9' 
##    domain_range = range(1,999+1)
##    chain1 = 'R'
##    chain2 = 'R'
##    chain2 = 'A'
##    exclude_chain = ''

##    ## luciferase
##    pdb1 = '1ba3'
##    pdb2 = '1lci'
##    domain_range = range(1,999+1)
##    chain1 = 'A'
##    chain2 = 'A'
##    exclude_chain = ''

##    ## G3P DH
##    pdb1 = '1gd1'
##    pdb2 = '2gd1'
##    domain_range = range(1,999+1)
##    chain1 = 'O'
##    chain2 = 'O'
##    exclude_chain = ''

##    ## hexokinase
##    pdb1 = '1hkg'
##    pdb2 = '2yhx'
##    domain_range = range(1,999+1)
##    chain1 = 'A'
##    chain2 = 'A'
##    exclude_chain = ''

##    ## adk
##    pdb1 = '1ake'
##    pdb2 = '4ake'
##    domain_range = range(1,999+1)
##    chain1 = 'A'
##    chain2 = 'A'
##    exclude_chain = 'B'

##    ## t4l
##    pdb1 = '2lzm'
##    pdb2 = '150l'
####    domain_range = range(15,59+1)
####    domain_range = range(60,80+1)
##    domain_range = range(81,162+1)
##    chain1 = 'A'
##    chain2 = 'D'
##    exclude_chain = 'B'

    l_input = [
##        ##
##        {'pdb1':'1ipd','pdb2':'1osj','chain1':'A','chain2':'A','range':range(1,98+1)+range(253,345+1)},
##        {'pdb1':'1ipd','pdb2':'1osj','chain1':'A','chain2':'A','range':range(99,108+1)+range(109,252+1)},
###### shears
##        ## aspartate amino transferase
##        {'pdb1':'9aat','pdb2':'1ama','chain1':'A','chain2':'A','range':range(15,36+1)+range(349,410+1)},
##        {'pdb1':'9aat','pdb2':'1ama','chain1':'A','chain2':'A','range':range(50,312+1)},
        ## alcohol dehydrogenase
        {'pdb1':'6adh','pdb2':'8adh1','chain1':'A','chain2':'A','range':range(1,174+1)+range(322,374+1)},
        {'pdb1':'6adh','pdb2':'8adh2','chain1':'A','chain2':'A','range':range(193,317+1)},
##        ## citrate synthase
##        {'pdb1':'1cts','pdb2':'4cts','chain1':'A','chain2':'A','range':range(1,276+1)+range(386,999+1)},
###### hinges
##        ## atpsulf
##        {'pdb1':'1i2d','pdb2':'1m8p','chain1':'A','chain2':'A','range':range(1,389+1)},
##        ## dnak (different spacegroups)
##        {'pdb1':'1dkx','pdb2':'1dky','chain1':'A','chain2':'A','range':range(389,509+1)},
##        ## dnak (different spacegroups)
##        {'pdb1':'1ddt','pdb2':'1mdt','chain1':'A','chain2':'A','range':range(1,376+1)},
##        ## ecpdpbp
##        {'pdb1':'1dpp','pdb2':'1dpe','chain1':'A','chain2':'A','range':range(1,260+1)+range(479,999+1)},
##        ## ef2
##        {'pdb1':'1n0v','pdb2':'1n0u','chain1':'C','chain2':'A','range':range(1,478+1)}, ## large
##        {'pdb1':'1n0v','pdb2':'1n0u','chain1':'C','chain2':'A','range':range(479,560+1)}, ## independent
##        {'pdb1':'1n0v','pdb2':'1n0u','chain1':'C','chain2':'A','range':range(561,9999+1)}, ## small
##        ## febp
##        {'pdb1':'1d9v','pdb2':'1mrp','chain1':'A','chain2':'A','range':range(109,227+1)+range(292,309+1)},
##        {'pdb1':'1d9v','pdb2':'1mrp','chain1':'A','chain2':'A','range':range(1,96+1)+range(228,262+1)},
##        ## folylpolyglutamate synthetase
##        {'pdb1':'1jbv','pdb2':'1jbw','chain1':'A','chain2':'A','range':range(1,295+1)},
##        {'pdb1':'1jbv','pdb2':'1jbw','chain1':'A','chain2':'A','range':range(296,386+1)},
##        ## glucose ABC transporter ATPase subunit (different spacegroups)
##        {'pdb1':'1oxs','pdb2':'1oxu','chain1':'C','chain2':'C','range':range(1,209+1)},
##        {'pdb1':'1oxs','pdb2':'1oxu','chain1':'C','chain2':'C','range':range(244,999+1)},
##        ## groel domain
##        {'pdb1':'1aon','pdb2':'1oel','chain1':'A','chain2':'A','range':range(1,137+1)+range(410,999+1)},
##        {'pdb1':'1aon','pdb2':'1oel','chain1':'A','chain2':'A','range':range(192,374+1)},
##        {'pdb1':'1aon','pdb2':'1oel','chain1':'A','chain2':'A','range':range(138,190+1)+range(375,409+1)},
##        ## lao bp
##        {'pdb1':'2lao','pdb2':'1laf','chain1':'A','chain2':'E','range':range(1,90+1)+range(192,238+1)},
##        {'pdb1':'2lao','pdb2':'1laf','chain1':'A','chain2':'E','range':range(91,191+1)},
##        ## t4l
##        {'pdb1':'1l96','pdb2':'1l97','chain1':'A','chain2':'A','range':range(13,59+1)},
##        {'pdb1':'1l96','pdb2':'1l97','chain1':'A','chain2':'A','range':range(81,164+1)},
##        ## maltodextrin bp
##        {'pdb1':'1omp','pdb2':'3mbp','chain1':'A','chain2':'A','range':range(1,104+1)+range(268,313+1)},
##        {'pdb1':'1omp','pdb2':'3mbp','chain1':'A','chain2':'A','range':range(113,258+1)+range(314,370+1)},
##        ## mRNA capping enzyme
##        {'pdb1':'1ckm','pdb2':'1ckm','chain1':'A','chain2':'B','range':range(1,237+1)+range(319,327+1)},
##        {'pdb1':'1ckm','pdb2':'1ckm','chain1':'A','chain2':'B','range':range(241,303+1)},
##        ## mura
##        {'pdb1':'1ejd','pdb2':'1a2n','chain1':'A','chain2':'A','range':range(1,20+1)+range(230,417+1)},
##        {'pdb1':'1ejd','pdb2':'1a2n','chain1':'A','chain2':'A','range':range(20,230+1)},
##        ## oligopeptide bp
##        {'pdb1':'1rkm','pdb2':'2rkm','chain1':'A','chain2':'A','range':range(1,263+1)+range(491,517+1)},
##        {'pdb1':'1rkm','pdb2':'2rkm','chain1':'A','chain2':'A','range':range(277,477+1)},
##        ## protein kinase A
##        {'pdb1':'1jlu','pdb2':'1cmk','chain1':'E','chain2':'E','range':range(1,33+1)+range(125,310+1),
##        {'pdb1':'1jlu','pdb2':'1cmk','chain1':'E','chain2':'E','range':range(34,124+1)},
##        ## dna polymerase beta
##        {'pdb1':'1bpd','pdb2':'2bpg','chain1':'A','chain2':'A','range':range(1,82+1)},
##        {'pdb1':'1bpd','pdb2':'2bpg','chain1':'A','chain2':'A','range':range(106,132+1)},
##        {'pdb1':'1bpd','pdb2':'2bpg','chain1':'A','chain2':'A','range':range(148,262+1)},
##        {'pdb1':'1bpd','pdb2':'2bpg','chain1':'A','chain2':'A','range':range(262,335+1)},
##        ## ribose bp
##        {'pdb1':'1urp','pdb2':'2dri','chain1':'A','chain2':'A','range':range(1,98+1)+range(235,259+1)},
##        {'pdb1':'1urp','pdb2':'2dri','chain1':'A','chain2':'A','range':range(104,234+1)+range(265,271+1)},
##        ## thioredoxin reductase
##        {'pdb1':'1tde','pdb2':'1f6m','chain1':'A','chain2':'E','range':range(1,112+1)+range(248,320+1)},
##        {'pdb1':'1tde','pdb2':'1f6m','chain1':'A','chain2':'E','range':range(118,242+1)},
##        ## dna bp
##        {'pdb1':'1fgu','pdb2':'1jmc','chain1':'A','chain2':'A','range':range(183,283+1)},
##        ## transferrin
##        {'pdb1':'1bp5','pdb2':'1a8e','chain1':'A','chain2':'A','range':range(1,75+1)+range(249,316+1)},
##        {'pdb1':'1bp5','pdb2':'1a8e','chain1':'A','chain2':'A','range':range(103,242+1)},
##        ## uracil dna glycosylase
##        {'pdb1':'1ssp','pdb2':'1akz','chain1':'E','chain2':'A','range':range(82,144+1)+range(191,240+1)},
##        {'pdb1':'1ssp','pdb2':'1akz','chain1':'E','chain2':'A','range':range(166,182+1)+range(270,304+1)},
        ]

    for i in range(len(l_input)):

        pdb1 = l_input[i]['pdb1']
        pdb2 = l_input[i]['pdb2']
        chain1 = l_input[i]['chain1']
        chain2 = l_input[i]['chain2']
        domain_range = l_input[i]['range']

        os.system('cp /oxygenase_local/data/pdb/%s/pdb%s.ent %s.pdb' %(pdb1[1:3],pdb1,pdb1,))
        os.system('cp /oxygenase_local/data/pdb/%s/pdb%s.ent %s.pdb' %(pdb2[1:3],pdb2[:4],pdb2[:4],))

        ss_range1, l_missing1, seqres1, l_modres = parse_header(pdb1, chain1,)
        ss_range2, l_missing2, seqres2, l_modres = parse_header(pdb2[:4], chain2,)

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

        l_coordinates1 = parse_coordinates(pdb1,chain1,domain_range,ss_range,l_missing,)
        l_coordinates2 = parse_coordinates(pdb2[:4],chain2,domain_range,ss_range,l_missing,)

##        l_coordinates1 = l_coordinates1[l2:len(l_coordinates1)-r2]
##        l_coordinates2 = l_coordinates2[l1:len(l_coordinates2)-r1]

        if len(l_coordinates1) != len(l_coordinates2):
            print len(l_coordinates1)
            print len(l_coordinates2)
            stop

        rmsd = instance_geometry.superpose(l_coordinates1,l_coordinates2)
        print pdb1, pdb2, round(rmsd,1), len(l_coordinates1)/3.
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter

        apply_transformation_matrix(
            pdb1,chain1,l_modres,
            [0,0,0],[[1,0,0],[0,1,0],[0,0,1]],[0,0,0],
            )
        apply_transformation_matrix(
            pdb2,chain2,l_modres,
            tv1,rm,tv2,
            )

        l_coordinates1 = parse_coordinates(pdb1+'_rotated',chain1,range(1,9999),ss_range,l_missing,)
        l_coordinates2 = parse_coordinates(pdb2+'_rotated',chain2,range(1,9999),ss_range,l_missing,)

        SUM = 0.
        n = len(l_coordinates1)
        for i in range(n):
            SUM += sum((l_coordinates1[i]-l_coordinates2[i])**2)
        RMSD = math.sqrt(SUM/n)
        print RMSD

    return


def apply_transformation_matrix(pdb,chain,l_modres,tv1,rm,tv2,):

    import numpy, os

    fd = open('%s.pdb' %(pdb[:4]),'r')
    lines = fd.readlines()
    fd.close()

    lines2 = []
    for line in lines:
        record = line[:6].strip()
        ## not coordinate section
        if record not in ['ATOM','HETATM',]:
            lines2 += [line]
            continue
        if line[21] != chain:
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

    fd = open('%s_rotated.pdb' %(pdb),'w')
    fd.writelines(lines2)
    fd.close()

    return


def parse_header(pdb, chain,):

    fd = open('/oxygenase_local/data/pdb/%s/pdb%s.ent' %(pdb[1:3],pdb),'r')
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
            if chain == line[11]:
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


def parse_coordinates(pdb,chain,domain_range,ss_range,l_missing,):

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
        if atom_name not in ['N','CA','C',]:
            continue
        res_no = int(line[22:26])
        ## exclude flexible domain
        if res_no not in domain_range:
            continue
        ## exclude random coil
        if res_no not in ss_range:
            continue
        ## exclude missing residues
        if res_no in l_missing:
            continue
        ## exclude other chains than the one of interest
        if line[21] != chain:
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
    main()

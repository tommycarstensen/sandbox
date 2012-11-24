import os
import quakes

def rmsd2bfactor(
    pdb1, pdb2, biomolecule1, biomolecule2, rmsd, d_coordinates, d_header,
    tv1, rm, tv2,
    l_equivalent_chains, bmchains1, bmchains2,
    d_chains_intrapdb_sequence_identical, d_chains_interpdb_sequence_similar,
    d_ATOMseq,
    path_pdb,
    ):

    print 'rmsd2bfactor'

    biomolecule1 = str(biomolecule1).zfill(2)
    biomolecule2 = str(biomolecule2).zfill(2)

    ##
    ## append secondary structure information
    ##
    d_lines_ss = {
        pdb1:[],
        pdb2:[],
        }
    for pdb in d_lines_ss:
        ##
        ## read lines
        ##
        fd = open('%s/%s/pdb%s.ent' %(path_pdb, pdb.lower()[1:3], pdb.lower()),'r')
##            fd = open('C:\Users\Tommy Carstensen\Downloads\pdb\%s.pdb' %(pdb.lower(),),'r')
        lines = fd.readlines()
        fd.close()
        for line in lines:
            record = line[:6].strip()
            if record in ['HELIX','SHEET','TURN',]:
                d_lines_ss[pdb] += line


    chains1 = l_equivalent_chains[0]
    chains2 = l_equivalent_chains[1]
    if len(chains1) != len(chains2):
        print pdb1, pdb2
        print chains1, chains2
        notexpected


    ##
    ## parse coordinates for the specific combinations of chain1 and chain2
    ##
    instance_quakes = quakes.quakes()
    (
        coordinates1, coordinates2, residue_count, d_lines, l_RMSDs,
        d_coordinates1, d_coordinates2,
        ) = instance_quakes.prepare_coords_for_alignment(
        pdb1,pdb2,chains1,chains2,
        d_chains_intrapdb_sequence_identical,
        d_chains_interpdb_sequence_similar,
        d_coordinates,d_ATOMseq,
        rmsd=rmsd, tv1=tv1, rm=rm, tv2=tv2,
        )

    ##
    ## append lines by model to final lines
    ##
    pdblines1 = d_lines_ss[pdb1]
    pdblines2 = d_lines_ss[pdb2]
    for model in d_lines[pdb1].keys():
        pdblines1 += ['MODEL     %4s                                                                  \n' %(model)]
        pdblines1 += d_lines[pdb1][model]
        pdblines1 += ['ENDMDL                                                                          \n']
    for model in d_lines[pdb2].keys():
        pdblines2 += ['MODEL     %4s                                                                  \n' %(model)]
        pdblines2 += d_lines[pdb2][model]
        pdblines2 += ['ENDMDL                                                                          \n']

    fd = open('pdb/%s/%s%s%s%s.pdb' %(pdb1[1], pdb1, str(biomolecule1).zfill(2), pdb2, str(biomolecule2).zfill(2)), 'w')
    fd.writelines(pdblines1)
    fd.close()
    fd = open('pdb/%s/%s%s%s%s.pdb' %(pdb2[1], pdb2, str(biomolecule2).zfill(2), pdb1, str(biomolecule1).zfill(2)), 'w')
    fd.writelines(pdblines2)
    fd.close()

    ##
    ## gif thumbnails
    ##
    d_biomolecules = {pdb1:biomolecule1,pdb2:biomolecule2}
    for pdb in d_biomolecules.keys():
        print 'gif thumbnail', pdb
        biomolecule = d_biomolecules[pdb]
        if pdb == pdb1:
            prefix = pdb1+str(biomolecule1).zfill(2)+pdb2+str(biomolecule2).zfill(2)
        if pdb == pdb2:
            prefix = pdb2+str(biomolecule2).zfill(2)+pdb1+str(biomolecule1).zfill(2)
        ## write rasmol script
        lines = [
            'rasmol -nodisplay pdb/%s/%s.pdb <<EOF\n' %(pdb[1], prefix),
            'color temperature\n',
            'wireframe 0\n',
            'cartoon\n',
            'write ppm/%s.ppm\n' %(prefix),
            'exit\n',
            ]
        ## write rasmol script to file
        fd = open('src/%srasmol.src' %(prefix),'w')
        fd.writelines(lines)
        fd.close()
        ## execute rasmol script
        os.system('source src/%srasmol.src > log/%srasmol.log' %(prefix,prefix))
        ## convert rasmol output
        os.system('convert ppm/%s.ppm -resize x80 gif/%s/%s.gif' %(prefix,pdb[1],prefix))
        ## clean up
        os.remove('ppm/%s.ppm' %(prefix))
        os.remove('log/%srasmol.log' %(prefix))
        os.remove('src/%srasmol.src' %(prefix))
    
    return

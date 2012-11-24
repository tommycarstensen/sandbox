#!/software/bin/python
#
#$Id: maprmsd2struc.py 236 2007-05-11 16:17:11Z tc $
#
#Tommy Carstensen, University College Dublin, 2007

# check what is wrong with 1atz!!!

class rmsdcolor:


    def main(self, pdb1, pdb2, chains1, chains2):

        import Numeric, math, os

        import quakes
        instance_quakes = quakes.quakes()

        import sys
        sys.path.append('/home/people/tc/python/Protool/')
        import geometry
        instance_geometry = geometry.geometry()

        pdbs = [pdb1,pdb2]

        ##
        ## parse coordinates to dictionary
        ##
        d_noncoordinates = {}
        d_coordinates = {}
        for pdb in pdbs:

            ## read pdb file
            fd = open('%s%s.pdb' %(self.pdbpath, pdb),'r')
            lines = fd.readlines()
            fd.close()

            ## parse pdb file
            d_noncoords = instance_quakes.parse_pdbnoncoordinatesections(lines, pdb)
            d_coords = instance_quakes.parse_pdbcoordinatesection(lines, pdb)
            d_noncoordinates[pdb] = d_noncoords
            d_coordinates[pdb] = d_coords

        ##
        ## convert dictionary of coordinates to list of coordinates
        ##
        coordinates1 = []
        coordinates2 = []
        d_alignATOMseq = {}
        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2 = instance_quakes.alignATOMseq(d_coordinates, d_noncoordinates, pdb1, pdb2, chain1, chain2)
            d_alignATOMseq[chain2] = {}
            d_alignATOMseq[chain2]['d_res_nos1'] = d_res_nos1
            d_alignATOMseq[chain2]['d_res_nos2'] = d_res_nos2
            d_alignATOMseq[chain2]['l1'] = l1
            d_alignATOMseq[chain2]['l2'] = l2
            coords1, coords2, rescount = instance_quakes.ATOMrecords2coordinates(d_coordinates, pdb1, pdb2, chain1, chain2, d_res_nos1, d_res_nos2, l1, l2, ATOMseq1, ATOMseq2)

            coordinates1 += coords1
            coordinates2 += coords2

        ##
        ## align coordinates
        ##
        rmsd = instance_geometry.superpose(coordinates1,coordinates2)
        tv1 = instance_geometry.fitcenter
        rm = instance_geometry.rotation
        tv2 = instance_geometry.refcenter

        ##
        ## calculate RMSD of transformed coordinates of individual residues
        ##
        lines1 = []
        lines2 = []
        atom_no = 0
        for i in range(len(chains1)):

            chain1 = chains1[i]
            chain2 = chains2[i]

            for SEQRESpos1 in range(l2,len(ATOMseq1)+l2):

                SEQRESpos2 = SEQRESpos1+l1-l2
                res_no1 = d_res_nos1[SEQRESpos1]['res_no']
                res_no2 = d_res_nos2[SEQRESpos2]['res_no']
                if res_no1 == '-' or res_no2 == '-':
                    continue
                iCode1 = d_res_nos1[SEQRESpos1]['iCode']
                iCode2 = d_res_nos2[SEQRESpos2]['iCode']

                if res_no1 not in d_coordinates[pdb1]['chains'][chain1]['residues'].keys():
                    continue
                if res_no2 not in d_coordinates[pdb2]['chains'][chain2]['residues'].keys():
                    continue
                if 'REMARK' in d_coordinates[pdb1]['chains'][chain1]['residues'][res_no1]['d_iCodes'][iCode1].keys():
                    continue
                if 'REMARK' in d_coordinates[pdb2]['chains'][chain2]['residues'][res_no2]['d_iCodes'][iCode2].keys():
                    continue

                d_resname = {
                    pdb1:{'chain':chain1,'res_no':res_no1,'iCode':iCode1},
                    pdb2:{'chain':chain2,'res_no':res_no2,'iCode':iCode2},
                    }
                for pdb in d_resname.keys():
                    chain = d_resname[pdb]['chain']
                    res_no = d_resname[pdb]['res_no']
                    iCode = d_resname[pdb]['iCode']
                    res_name = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['res_name']
                    d_atoms = d_coordinates[pdb]['chains'][chain]['residues'][res_no]['d_iCodes'][iCode]['atoms']
                    d_resname[pdb]['res_name'] = res_name
                    d_resname[pdb]['d_atoms'] = d_atoms
                res_name1 = d_resname[pdb1]['res_name']
                res_name2 = d_resname[pdb2]['res_name']
                d_atoms1 = d_resname[pdb1]['d_atoms']
                d_atoms2 = d_resname[pdb2]['d_atoms']
##                    print pdb1, pdb2, chain1, chain2, res_no1, res_no2, res_name1, res_name2

                ##
                ## calculate rmsd for the atoms of the residue
                ##
                SS = []
                for atom_name in d_atoms1.keys():
                    if atom_name not in d_atoms2.keys():
                        continue
                    coordinate1 = d_atoms1[atom_name]['coordinate']
                    coordinate2 = Numeric.matrixmultiply(rm, d_atoms2[atom_name]['coordinate']-tv1)+tv2
                    SS += [sum((coordinate2-coordinate1)**2)]
                RMSD = math.sqrt(sum(SS)/len(SS))
                occupancy = bfactor = RMSD/rmsd
                ##
                ## append atoms and coordinates to lines
                ##
                for atom_name in d_atoms1.keys():
                    if atom_name not in d_atoms2.keys():
                        continue
                    coordinate1 = d_atoms1[atom_name]['coordinate']
                    coordinate2 = d_atoms2[atom_name]['coordinate']
                    x1 = coordinate1[0]; y1 = coordinate1[1]; z1 = coordinate1[2]
                    x2 = coordinate2[0]; y2 = coordinate2[1]; z2 = coordinate1[2]
                    atom_no += 1
                    altloc = ''
                    charge = ''
                    if 'H' in atom_name:
                        element = 'H'
                    else:
                        element = atom_name[0]
                    lines1 += [
                        '%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n'
                        %('ATOM'.ljust(6), atom_no, atom_name.ljust(4), altloc, res_name1.ljust(3), chain1, res_no1, iCode1, x1, y1, z1, occupancy, bfactor, element.rjust(2), charge.rjust(2))
                        ]
                    lines2 += [
                        '%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n'
                        %('ATOM'.ljust(6), atom_no, atom_name.ljust(4), altloc, res_name2.ljust(3), chain2, res_no2, iCode2, x2, y2, z2, occupancy, bfactor, element.rjust(2), charge.rjust(2))
                        ]

        fd = open('%s%s.pdb' %(self.path_out,pdb1),'w')
        fd.writelines(lines1)
        fd.close()
        fd = open('%s%s.pdb' %(self.path_out,pdb2),'w')
        fd.writelines(lines2)
        fd.close()

        for pdb in pdbs:
            ## write rasmol script
            lines = [
                'rasmol -nodisplay %s%s.pdb << EOF\n' %(self.path_out,pdb),
                'color temperature\n',
                'spacefill\n',
                'write %s%s.ppm\n' %(self.path_tmp,pdb),
                'exit\n',
                ]
            ## write rasmol script to file
            fd = open('%s%srasmol.src' %(self.path_tmp,pdb),'w')
            fd.writelines(lines)
            fd.close()
            ## execute rasmol script
            os.system('source %s%srasmol.src' %(self.path_tmp,pdb))
            ## convert rasmol output
            os.system('convert %s%s.ppm -resize x80 %s%s.gif' %(self.path_tmp,pdb, self.path_out,pdb))
            ## clean up
            os.system('rm %s%s.ppm' %(self.path_tmp,pdb))
        ## clean up
        os.system('rm %s%srasmol.src' %(self.path_tmp,pdb))
        
        os.system

        return


    def __init__(self):

        self.pdbpath = '/oxygenase_local/data/pdb/'
        self.path_out = 'out/'
        self.path_tmp = 'tmp/'

        return


if __name__ == '__main__':
    instance = rmsdcolor()
    instance.main({'1aq8':['A','B','C'],'1as7':['A','B','C']})

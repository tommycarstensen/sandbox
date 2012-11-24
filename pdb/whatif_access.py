    def whatif_accessibility(self,pdb,chain,res_no):

        l_acc = []
        source = '%s/whatif/%s.txt' %(self.path_cwd,pdb)
        if 1 == 1:
##        if not os.path.isfile(source):
            fd = open(source,'w')
            fd.writelines([
##                '/software/whatif/DO_WHATIF.COM <<EOF\n',
                '/local/software/whatif/DO_WHATIF.COM <<EOF\n',
                'GETMOL %s/%s/pdb%s.ent\n' %(self.path_pdb, pdb[1:3], pdb),
                '%s\n' %(pdb),
                '%DELWAT\n',
                ##'WATRAD=1.4\n',
                ##'ACCTYP=1\n',
                '%SETACC\n',
                'NOWAT 0\n',
                'NOWAT 0\n',
                '%LISTA\n',
                'TOT 0\n',
                'STOP\n',
                'Y\n',
                ])
            fd.close()
            os.system('source %s > whatif/%s.out' %(source, pdb))

            ## clean up the mess that whatif left
            for fn in ['ALTERR.LOG','PDBFILE','WHATIF.FIG','DAVADRUG.PDB',]:
                time.sleep(0.05)
                if os.path.isfile(fn):
                    try:
                        os.remove(fn)
                    except:
                        None
            l = os.listdir(os.getcwd())
            for s in l:
                if s[:3] == 'DRG' and s[6] == '.':
                    time.sleep(0.01)
                    if os.path.isfile(s):
                        try:
                            os.remove(s)
                        except:
                            None

            fd = open('whatif/%s.out' %(pdb),'r')
            lines = fd.readlines()
            fd.close()
            for i in range(len(lines)):
                line = lines[i]
                if line[:37] == 'New accessibilities calculated ... : ':
                    k = i
                    bool_break = False
                    ## loop until residue
                    for j in range(i+2,len(lines)):
                        line = lines[j]
                        if 'Residue:' in line:
                            index = line.index('Residue:')
                            res_name = lines[j][index+15:index+18]
                            if res_name == 'HOH':
                                continue
                            if lines[j][index+21:index+25].strip() == 'OXT':
                                continue
                            if not (
                                chain == lines[j][index+29]
                                and
                                res_no == int(lines[j][index+21:index+25])
                                ):
                                continue
##                            res_index = int(lines[j][index+10:index+14])
##                            res_no = int(lines[j][index+21:index+25])
                            ## loop until column headers of residue
                            for k in range(j+1,len(lines)):
                                line = lines[k]
                                if 'Atom     X     Y     Z   Acc   B   WT   VdW  Colr   AtOK  Val' in line:
                                    ## loop over data lines
                                    for l in range(k+1,len(lines)):
                                        line = lines[l]
                                        if line == ' \n':
                                            continue
                                        ## next residue
                                        if 'Residue:' in line:
                                            break
                                        ## next residue (WHATIF failed..?)
                                        if 'Option not found' in line:
                                            break
                                        ## next residue
                                        if 'Which in turn is attached to' in line:
                                            break
                                        atom_name = line[:4].strip()
                                        ## WHATIF2PDB_nomenclature
                                        if atom_name == "O'":
                                            atom_name = 'O'
                                        if (
                                            chain == lines[j][index+29]
                                            and res_no == int(lines[j][index+21:index+25])
##                                            and atom_name not in ['N','CA','C','O',]
##                                    and (
##                                        (res_name not in ['VAL','ILE','THR','SER','CYS',] and atom_name in ['CB','CG',])
##                                        or
##                                        (res_name in ['VAL','ILE',] and atom_name in ['CB','CG1',])
##                                        or
##                                        (res_name in ['CYS',] and atom_name in ['CB','SG',])
##                                        or
##                                        (res_name in ['SER',] and atom_name in ['CB','OG',])
##                                        or
##                                        (res_name in ['THR',] and atom_name in ['CB','CG2',])
##                                        )
                                            ):
                                            if line[24:28].strip() == '---':
                                                continue
##                                            print pdb
##                                            print 'i', lines[i]
##                                            print 'j', lines[j]
##                                            print 'k', lines[k]
##                                            print 'l', lines[l]
                                            acc = float(line[24:28])
                                            l_acc += [acc]
                                            if len(l_acc) > 20:
                                                stop
                                            bool_break = True
##                                            print 'acc', pdb, chain, int(lines[j][index+21:index+25]), res_name, atom_name, acc
                                if bool_break == True:
                                    break
                        if bool_break == True:
                            break
                    if bool_break == True:
                        break

        if len(l_acc) == 0: ## residue not present
            accessibility = 0
        elif len(l_acc) > 14: ## TRP
            print l_acc
            print len(l_acc)
            stop
        else:
            accessibility = sum(l_acc)/len(l_acc)

        return accessibility



    def whatif(self,pdb,biomolecule,d_coordinates):

        source = '%s/whatif/%s.txt' %(self.path_cwd,pdb1)
        if not os.path.isfile(source):
            fd = open(source,'w')
            fd.writelines([
                '/software/whatif/DO_WHATIF.COM <<EOF\n',
                'GETMOL %s/%s/pdb%s.ent\n' %(self.path_pdb, pdb1[1:3], pdb1),
                '%s\n' %(pdb1),
                '%DELWAT\n',
                ##'WATRAD=1.4\n',
                ##'ACCTYP=1\n',
                '%SETACC\n',
                'NOWAT 0\n',
                'NOWAT 0\n',
                '%LISTA\n',
                'TOT 0\n',
                'STOP\n',
                'Y\n',
                ])
            fd.close()
            os.system('source %s > whatif/%s_%s.out' %(source, pdb1, biomolecule1))
            fd = open('whatif/%s_%s.out' %(pdb1, biomolecule1),'r')
            lines = fd.readlines()
            fd.close()
            for i in range(len(lines)):
                line = lines[i]
                if line[:37] == 'New accessibilities calculated ... : ':
                    k = i
                    for j in range(i+2,len(lines)):
                        print j
                        if j < k+1:
                            continue
                        line = lines[j]
                        if 'Residue:' in line:
                            index = line.index('Residue:')
                            chain = line[index+29]
                            res_no = int(line[index+21:index+25])
                            for k in range(j+2,len(lines)):
                                line = lines[k]
                                if line == ' \n':
                                    break
                                atom_name = line[:4].strip()
                                acc = float(line[24:28])
                                if atom_name == "O'" and line[67] == 'O':
                                    atom_name = 'O'
                                d_coordinates['chains'][chain]['residues'][res_no]['d_iCodes'][' ']['atoms'][atom_name]['acc'] = acc
                        else:
                            break
                    break

        os.remove('ALTERR.LOG')
        os.remove('PDBFILE')
        os.remove('WHATIF.FIG')

        return d_coordinates

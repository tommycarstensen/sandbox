#!/usr/bin/env python
#$Id: gromacs.py 239 2007-05-16 14:16:02Z tc $
import os, sys, time, socket

hostname = socket.gethostname()

class gromacs:

    def main(self):

        pdb = self.pdb.replace('.pdb','')

        T = self.T
        t = self.t

        ## Gromacs force field
        self.path_gromacs = ''
        ff = 4

        if '-cluster' in sys.argv:
            self.path_gromacs = '/usr/local/cluster/gromacs_4/bin'
            ff = 10 ## AMBER force field when 4.0.5
            s_ff = 'amber'
        else:
            ## AMBER force field currently only available with Gromacs 3.3.1 and lower (2009)
            self.path_gromacs = '/local/software/gromacs-4.0.5/bin'
            ff = 10
            s_ff = 'amber'
##            self.path_gromacs = '/software/Gromacs/v.3.3.1/bin'
##            ff = 5 ## GROMOS96
##            s_ff = 'gromacs'

        index_pre_hydrogen = 'crystalcontacts.ndx'
        index = 'crystalcontacts_pdb2gmx.ndx' ## post hydrogen (pdb2gmx) index file

        ## calculate charge of protein by counting charged residues
        charge, s_pdb2gmx = self.preprocess_pdb(pdb,s_ff,)
        time.sleep(1)

        ##
        ## gromacs initialization
        ##
        
        ## add hydrogen atoms and generate a gromacs topology file
        self.pdb2gmx(pdb, ff=ff, s_ff=s_ff, s_pdb2gmx=s_pdb2gmx,)
        time.sleep(1)
        ## enlarge the box to accommodate for water (and center protein)
        self.editconf(pdb,)
        time.sleep(1)
        ## generate the box and fill it with water
        ## i.e. solvate protein and change topology file to include water
        self.genbox(pdb)
        time.sleep(1)

        ##
        ## energy minimization (minimize repulsive contacts between protein and solvent)
        ##

        ## generate mdrun input file (binary topology file .tpr)
        ## i.e. preprocess topology (.top), structure (.gro), parameters (.mdp)
        self.grompp(pdb, runtype='EM', i_suffix='genbox.gro',)
        time.sleep(1)

        if charge != 0:
            ## generate ions (i.e. generate net-zero total charge for the system)
            self.genion(pdb, charge, ff, s_ff,)
            self.check_ion_positions(pdb,charge,ff=ff,)
        else:
            stop_and_write_pdb2gmx_output_to_genion_output_slash_grompp_input

        ## rename cysteines
        if s_ff == 'amber':
            self.rename_cysteines(pdb,)

        ## preprocess generic structure from genion
        self.grompp(pdb, runtype='EM', i_suffix='genion.gro',)
        time.sleep(1)
        ## energy minimize protein and solvent
        ## i.e. minimize repulsive contacts between protein and solvent
        self.mdrun(pdb, runtype='EM')
        time.sleep(1)

        ##
        ## position restrained dynamics ("fill" protein cavities with water)
        ##
        self.grompp(pdb, runtype='PR', i_suffix='EM.gro', T=T,)
        time.sleep(1)
        self.mdrun(pdb, runtype='PR')
        time.sleep(1)

        ##
        ## molecular dynamics
        ##
        self.grompp(pdb, runtype='MD', i_suffix='PR.gro', t=t, T=T,) 
        time.sleep(1)
        self.mdrun(pdb, runtype='MD')
        time.sleep(1)



        ## check convergence
        print 'g_energy -f 2vb1_MDener.edr -o 2vb1_energy.xvg'
        print '9 0' ## potential energy

        ## convert trajectory
        print '-skip skip*dt*nstxout ps (default 5*2ps yielding 1000 frames if 10ns simulation)'
        print 'trjconv -f 2vb1_MD.trr -s 2vb1_MD.tpr -fit rot+trans -sep -o trjconv/2vb1_MD.pdb -skip'
        ## fit, backbone = 4
        ## output, system = 0

        ## multiple processors
        print 'mpirun -np xx mdrun'

        ##
        ## resume run
        ##
        ## tpr file generation
        print 'tpbconv -s 2vb1_MD1.tpr -f 2vb1_MD1.trr -e 2vb1_MDener1.edr -o 2vb1_MD2.tpr -time xxxx' ## -time (time found in 2nd column of log)
        ## run based on new tpr file
        print 'mdrun -s 2vb1_MD2.tpr -o 2vb1_MD2.trr -c 2vb1_MD.gro -e 2vb1_MDener2.edr -g 2vb1_MDmdrun2.log -v'

        ##
        ## extend run
        ##
        print 'tpbconv -extend 10000' ## ps
        print 'tpbconv -until 20000' ## ps
        print 'mdrun -s 2vb1_MD2.tpr -o 2vb1_MD2.trr -c 2vb1_MD.gro -e 2vb1_MDener2.edr -g 2vb1_MDmdrun.log -v'

        ##
        ## combine trr files
        ##
        print 'trjcat -o whole.trr -f part1.trr part2.trr'

        return


    def genpr(self,pdb,):

        source = '%s_genpr.src' %(pdb)
        lines = [
            '%s/genpr ' %(self.path_gromacs,),
            '-f %s_PR.gro ' %(pdb,),
            '-o %s_genpr.itp ' %(pdb,),
            ]
        lines += ['<< EOF\n4\n',]

        fd = open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))

        fd = open('%s.top' %(pdb),'r')
        lines = fd.readlines()
        fd.close()
        for i in range(len(lines)):
            if lines[i] == '; Include Position restraint file\n':
                lines = lines[:i+2]+['#include "%s_genpr.itp"\n' %(pdb)]+lines[i+3:]
        fd = open('%s_genpr.top' %(pdb),'w')
        fd.writelines(lines)
        fd.close()

        return


    def rename_cysteines(self,pdb,):

        fd = open('%s_genion.pdb' %(pdb),'r')
        lines_in = fd.readlines()
        fd.close()
        lines_out = []
        for line in lines_in:
            record = line[:6].strip()
            if record == 'ATOM':
                res_name = line[17:21].strip()
                if len(res_name) > 3:
                    stop
                if res_name == 'CYS':
                    line = line[:17]+'CYS2'+line[21:]
                lines_out += [line]
            else:
                lines_out += [line]
        fd = open('%s_genion.pdb' %(pdb),'w')
        fd.writelines(lines_out)
        fd.close()

        return
            

    def check_ion_positions(self,pdb,charge,ff=4,):

        import os

        source = '%s_pdb2gmx_genion.src' %(pdb)

        lines = [
            '%s/pdb2gmx ' %(self.path_gromacs),
            '-f %s_genion.gro ' %(pdb), ## input, generic structure
            '-o %s_genion.pdb ' %(pdb), ## output, generic structure
##            '-p %s.top ' %(pdb), ## output, topology file
##            '-i %s.itp ' %(pdb), ## output, include file for topology
            '<<EOF\n%s\n' %(ff), ## force field selection
            ]
        fd = open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))


##        fd = open('%s_genion.pdb' %(pdb),'r')
##        lines = fd.readlines(lines)
##        fd.close()
##
##        for line in lines:
##            record == line[:6].split():
##                if record == 'ATOM':
##                    if int(line[22:26])
##                    if line[17:20].strip() in ['Na','Cl',]:

        if not '-skipchecks' in self.arguments:
            print 'Check %s_genion.pdb to see if you are satisfied with the position of the %s ions? Press enter to continue.' %(pdb,abs(charge),)
            raw = raw_input()
        
        return


    def genion(self, pdb, charge, ff, s_ff,):

        if not os.path.isfile('%s_EM.tpr' %(pdb)):
            print '%s_EM.tpr' %(pdb), 'does not exist. A previous step failed.'
            stop

        source = '%s_genion.src' %(pdb)

        if charge > 0:
            chargeflag = 'nn'
            if s_ff == 'amber':
                iontype = 'Cl'
            elif s_ff == 'gromacs':
                iontype = 'CL-' ## Gromacs ff
        elif charge < 0:
            chargeflag = 'np'
            if s_ff == 'amber':
                iontype = 'Na'
            elif s_ff == 'gromacs':
                iontype = 'NA+' ## Gromacs ff

        lines = [
            '%s/genion ' %(self.path_gromacs),
            '-s %s_EM.tpr ' %(pdb),
            '-o %s_genion.gro ' %(pdb),
            '-%s %s ' %(chargeflag, abs(charge)),
            '-g %s_genion.log ' %(pdb),
            ]
        lines += ['<< EOF\n12\n',]

        fd=open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))


        ##
        ## modify topology file
        ##
        fd = open('%s.top' %(pdb),'r')
        lines = fd.readlines()
        fd.close()
        
        line_SOL = lines[-1]
        SOL_mols = int(line_SOL.split()[1])-abs(charge)
        line_sol = 'SOL     %13i\n' %(SOL_mols)

        line_ion = '%3s     %13i\n' %(iontype, abs(charge))

        lines = lines[:-1]+[line_sol]+[line_ion] ## ions after solvent when gro output
##        lines = lines[:-1]+[line_ion]+[line_sol] ## ions before solvent when pdb output

        fd = open('%s.top' %(pdb),'w')
        fd.writelines(lines)
        fd.close()


        return


    def mdrun(self, pdb, runtype=''):

        source = '%s_%smdrun.src' %(pdb, runtype)

        lines = []
        if self.parallel == True:
            lines += ['mpirun -np 4 ']
        lines += [
            '%s/mdrun ' %(self.path_gromacs),
            '-s %s_%s.tpr ' %(pdb, runtype), ## input, generic run input
            '-o %s_%s.trr ' %(pdb, runtype), ## output, trajectory
            '-c %s_%s.gro ' %(pdb, runtype), ## output, generic structure
            '-e %s_%sener.edr ' %(pdb, runtype), ## output, generic energy
            '-g %s_%smdrun.log ' %(pdb, runtype), ## output, log
            '-v ',
            '<< EOF\n0',
            ]

        fd=open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))

        return


    def grompp(self, pdb, runtype, i_suffix, index=None, t=2, T=298.15,):

        if not os.path.isfile('%s_%s' %(pdb,i_suffix,)):
            print '%s_%s' %(pdb,i_suffix,), 'does not exist. A previous step failed.'
            stop

        o_suffix = runtype

        if runtype == 'EM':

            d_mdp = {
##                'cpp':'/lib/cpp', ## obsolete in 4.0
                ## 7.3.2 preprocessing
                'define':'-DFLEX_SPC',
                ## 7.3.3 run control
                'integrator':'steep',
                'nsteps':10000., ## maxium number of steps to integrate
                ## 7.3.5 energy minimization
                'emtol':100.,
##                'emtol':1000.,
                'emstep':.01, ## maximum step size (nm)
                ## 7.3.10 electrostatics
                'coulombtype':'PME',
                ## 7.3.14 temperature coupling
                'tcoupl':'no',
                ## 7.3.15 pressure coupling
                'Pcoupl':'no',
                ## 7.3.17 velocity generation
                'gen_vel':'no',
                ## 7.3.18 bonds
                'constraints':'none',
                }

        if runtype == 'PR':

            d_mdp = {
##                'cpp':'/lib/cpp', ## obsolete in 4.0
                ## 7.3.2 preprocessing
                'define':'-DPOSRES',
                ## 7.3.3 run control
                'integrator':'md',
                'dt':.001, ## ps
##                'nsteps':50000,
                'nsteps':1,
##                ## 7.3.8 output control
##                'energygrps':'PROTEIN OTHER',
                ## 7.3.10 electrostatics
                'coulombtype':'PME',
                ## 7.3.14 temperature coupling
                'tcoupl':'v-rescale',
                'tc-grps':'PROTEIN OTHER', ## coupled group(s)
                'tau_t':'.1 .1', ## time constant(s) for coupling
                'ref_t':'%s %s' %(T,T,), ## reference temperature(s)
                ## 7.3.15 pressure coupling (NPT ensemble, maybe volume not constant to allow water to "spread out")
                'pcoupl':'no', ## No pressure coupling. This means a fixed box size.
##                'pcoupl':'berendsen',
##                'pcoupltype':'isotropic',
##                'compressibility':'4.5e-5',
##                'ref_p':1.,
                ## 7.3.18 bonds
##                'constraints':'all-bonds', ## restrain protein and only allow water to "relax"?
                'constraints':'none', ## restrain protein and only allow water to "relax"?
                }

        if runtype == 'MD':

            ## t (ns), 1000*t (ps)
            dt = .001 ## ps, 0.002 if Gromacs ff
            nsteps = 1000*t/dt
            d_mdp = {
##                'cpp':'/lib/cpp', ## obsolete in 4.0
                ## 7.3.3 run control
                'integrator':'md',
                'dt':dt, ## ps
                'nsteps':nsteps,
                ## 7.3.8 output control (## collect data every dt*nstxout ps)
                'nstxout':10000, ## 10000 frames = (1000*t/dt)/nstxout
                'nstvout':10000,
                'nstfout':0,
                'nstlog':10000,
                'nstenergy':10000,
                'energygrps':'PROTEIN OTHER', ## groups to write to energy file
                ## 7.3.10 electrostatics
                'coulombtype':'PME',
                ## 7.3.14 temperature coupling
                'tcoupl':'v-rescale', ## The Berendsen thermostat does not generate the correct kinetic energy distribution. You might want to consider using the V-rescale thermostat.
                'tc-grps':'PROTEIN OTHER',
                'tau_t':'.1 .1',
                'ref_t':'%s %s' %(T,T,),
                ## 7.3.15 pressure coupling (none, NVT ensemble)
                'pcoupl':'no', ## no pressure coupling, this means a fixed box size
                ## 7.3.17 velocity generation
                'gen_vel':'yes', ## yes or no?
                'gen_temp':T,
                ## 7.3.18 bonds
                'constraints':'none',
                }
            if self.posre == True:
                ## 7.3.2 preprocessing
                d_mdp['define'] = '-DPOSRES'

                self.genpr(pdb,)

##        ## freeze groups
##        if index != None:
##            d_mdp['freezegrps'] = 'crystalcontacts'
##            d_mdp['freezedim'] = 'Y Y Y'

        mdp = '%s_%sgrompp.mdp' %(pdb, runtype)
        fd=open(mdp,'w')
        for key in d_mdp.keys():
            fd.write('%s = %s\n' %(key,d_mdp[key]))
        fd.close()

        source = '%s_%sgrompp.src' %(pdb, i_suffix[:-3])

        lines = [
            '%s/grompp ' %(self.path_gromacs),
            '-f %s ' %(mdp), ## input, MD parameters
            '-po %s ' %(mdp), ## output, MD parameters
            '-c %s_%s ' %(pdb, i_suffix), ## input, generic structure
            ]
        if self.posre == True and runtype == 'MD':
            lines += [
                '-p %s_genpr.top ' %(pdb), ## input, topology file
                ]
        else:
            lines += [
                '-p %s.top ' %(pdb), ## input, topology file
                ]
        lines += [
            '-o %s_%s.tpr ' %(pdb, o_suffix), ## output, generic run input
            '-maxwarn 0 ',
            ]
##        if index != None:
##            lines += [
##                ' -n %s' %(index),
##                ]

        fd = open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))

        return


    def genbox(self, pdb):

        import os

        source = '%s_genbox.src' %(pdb)

        lines = [
            '%s/genbox ' %(self.path_gromacs),
            '-cp %s_editconf.gro ' %(pdb), ## input, generic structure (protein)
            '-cs ', ## input, generic structure (solvent)
            '-o %s_genbox.gro ' %(pdb), ## output, generic structure
            '-p %s.top' %(pdb), ## topology file
            ]

        fd = open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))

        if not os.path.isfile('%s_genbox.gro' %(pdb)):
            print '%s_genbox.gro' %(pdb), 'does not exist. This step failed.'
            stop

        return        


    def editconf(self, pdb,):

        import os

        source = '%s_editconf.src' %(pdb)

        if '-crystalcontacts' in sys.argv:
            fd = open('%s.pdb' %(pdb),'r')
            lines = fd.readlines()
            fd.close()
            for line in lines:
                record = line[:6].strip()
                if record == 'CRYST1':
                    print line
                    a = .1*float(line[6:15])
                    b = .1*float(line[15:24])
                    c = .1*float(line[24:33])
##                    a = .13*float(line[6:15])
##                    b = .13*float(line[15:24])
##                    c = .13*float(line[24:33])
                    alpha = float(line[33:40])
                    beta = float(line[40:47])
                    gamma = float(line[47:54])
                    spacegroup = line[55:66].strip()
                    if spacegroup != 'P 1':
                        stop_not_triclinic
                    break
            s = '%s/editconf -f %s_pdb2gmx.gro -o %s_editconf.gro -bt triclinic -box %f %f %f -angles %f %f %f' %(
                self.path_gromacs, pdb, pdb, a,b,c, alpha,beta,gamma,
                ) ## -box implies -c
        else:
            s = '%s/editconf -f %s_pdb2gmx.gro -o %s_editconf.gro -d 1.0 -bt cubic' %(self.path_gromacs, pdb, pdb) ## -d implies -c

        fd = open(source,'w')
        fd.write(s)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))

        return

    def pdb2gmx(self, pdb, ff=4, s_ff='gromacs', s_pdb2gmx = ''):

        import os

        source = '%s_pdb2gmx.src' %(pdb)

        lines = [
            '%s/pdb2gmx ' %(self.path_gromacs),
            '-f %s.pdb ' %(pdb), ## input, generic structure
            '-o %s_pdb2gmx.gro ' %(pdb), ## output, generic structure
            '-p %s.top ' %(pdb), ## output, topology file
            '-i %s.itp ' %(pdb), ## output, include file for topology
            '-water spc ', ## default
            '-ignh ', ## ignore hydrogen because of pdb/gromacs nomenclature differences
            ]
        ## if not amber then manual selection of protonation states from within gromacs
        if s_ff != 'amber':
            lines += [
                '-lys -asp -glu -his ', ## manual selection of protonation states (arg not in v3.3.1)
                ]
        ## if amber then manual selection of disulfides
        if s_ff == 'amber':
            lines += [
                '-ss ', ## manual selection of protonation states (arg not in v3.3.1)
                ]
        lines += [
            '<<EOF\n%s\n' %(ff), ## force field selection
            ]
        lines += [s_pdb2gmx]

        fd = open(source,'w')
        fd.writelines(lines)
        fd.close()

        if '-cluster' in sys.argv:
            os.system('source %s' %(source))
        else:
            os.system('source %s/%s' %(os.getcwd(),source,))

        return


    def preprocess_pdb(self,pdb,s_ff,):

        d_res = {}

        fd = open('%s.pdb' %(pdb),'r')
        lines = fd.readlines()
        fd.close()

        l_seqres = []
        lines2 = []
        for i in range(len(lines)):
            line = lines[i]
            record = line[:6].strip()
            if record == 'HETATM':
                if line[17:20] != 'HOH':
                    print line[:-1]
                continue
            elif record == 'SEQRES':
                l_seqres += line[19:70].split()
            elif record == 'ANISOU':
                if lines[i-1][:6] == 'HETATM':
                    if line[17:20] != 'HOH':
                        print line[:-1]
                continue
            elif record == 'CONECT':
                continue
            elif record == 'ATOM':
                altloc = line[16]
                if not altloc in [' ','A','1',]:
                    continue
                line = line[:16]+' '+line[17:]
            lines2 += [line]

        fd = open('%s.pdb' %(pdb),'w')
        fd.writelines(lines2)
        fd.close()

        ##
        ## calculate default overall charge
        ##
        charge = 0
        s = ''
        s_cys = ''
        i_cysteines = 0
        d_titratable = {
            'ASP':[],'GLU':[],
            'LYS':[],'ARG':[],
            'HIS':[],
            }
        for i in range(len(l_seqres)):
            res_name = l_seqres[i]
            if res_name in ['GLU','ASP',]:
                print res_name, i+1
                charge -= 1
                s += '%s%i,' %(res_name,i+1,)
                d_titratable[res_name] += [i]
            elif res_name in ['ARG','LYS',]:
                print res_name, i+1
                charge += 1
                if res_name == 'LYS':
                    s += '%s%i,' %(res_name,i+1,)
                    d_titratable[res_name] += [i]
            elif res_name in ['HIS',]:
                charge += 0
                s += '%s%i,' %(res_name,i+1,)
                d_titratable[res_name] += [i]
            elif res_name in ['CYS',]:
                s_cys += '%s%i,' %(res_name,i+1,)
                i_cysteines += 1
            elif res_name not in [
                'ALA','CYS','PHE','GLY','ILE','LEU','MET',
                'ASN','PRO','GLN','SER','THR','VAL','TRP','TYR',
                ]:
                print res_name
                stop

        print self.arguments
        if not '-nondefaultcharge' in self.arguments:
            print 'GLU,ASP,ARG,LYS are by default charged. HIS is by default not charged (HD1 proton present).'
            print 'N.B. The input pdb file is automatically edited, if the amber force field is being used.'
            print 'Use flag -nondefaultcharge'
            print "Type in the sequence number (separated by commas) of the following residues you want *not* to have a side chain with default charge:"
            print s
            raw = raw_input()
        else:
            raw = self.arguments[self.arguments.index('-nondefaultcharge')+1]
            if raw == 'None':
                raw = ''

        if not '-nodisulfide' in self.arguments:
            print 'Which cysteins do you want *not* to interact in disulfide bonds?'
            print 'Use flag -nodisulfide'
            print "Type in the sequence number (separated by commas) of the following residues you want *not* to interact in disulfide bonds:"
            print s_cys
            raw_cys = raw_input()
        else:
            raw_cys = self.arguments[self.arguments.index('-nodisulfide')+1]
            if raw_cys == 'None':
                raw_cys == ''
##            raise ('add flag to pdb2gmx to select disulfide bonds')

        l_seqIDs = raw.split(',')
        if l_seqIDs != [''] and l_seqIDs != ['None']:
            for i in range(len(l_seqIDs)):
                seqID = int(l_seqIDs[i])-1
                l_seqIDs[i] = seqID
                if l_seqres[seqID] in ['ASP','GLU',]:
                    charge += 1
                elif l_seqres[seqID] in ['HIS',]:
                    charge += 1
                elif l_seqres[seqID] in ['LYS','ARG',]:
                    charge -= 1
                else:
                    print seqID, l_seqres[seqID]
                    stop

        s_pdb2gmx = ''

        ## keep all disulfide bonds and assume all cysteines make disulfide bonds
        for i in range(0,i_cysteines,2):
            s_pdb2gmx += 'y\n'
        
        for res_name in ['LYS','ASP','GLU','HIS',]:
            for seqID in d_titratable[res_name]:
                ## default protonation
                if seqID not in l_seqIDs:
                    if res_name == 'LYS':
                        s_pdb2gmx += '1\n' ## protonated
                    elif res_name in ['ASP','GLU',]:
                        s_pdb2gmx += '0\n' ## not protonated
                    elif res_name == 'HIS':
                        s_pdb2gmx += '0\n' ##
                ## non-default protonation
                else:
                    if res_name == 'LYS':
                        s_pdb2gmx += '0\n' ## not protonated
                    elif res_name in ['ASP','GLU',]:
                        s_pdb2gmx += '1\n' ## protonated
                    elif res_name == 'HIS':
                        s_pdb2gmx += '2\n' ## protonated (HIP)

        ##
        ## write pdb for amber
        ##
        if s_ff == 'amber':

            fd = open('%s.pdb' %(pdb),'r')
            lines = fd.readlines()
            fd.close()

            d_titratable = {
                'ASP':{'nondefault':'ASH','default':'ASP',},
                'GLU':{'nondefault':'GLH','default':'GLU',},
                'LYS':{'nondefault':'LYN','default':'LYP',},
                'CYS':{'default':'CYN','nondefault':'CYM',},
                'HIS':{'default':'HID','nondefault':'HIP',}, ## HID if HD1, HIE if HE2
                }
            d_terminal = {
                'O':'OC1',
                'OXT':'OC2',
                }
            lines2 = []
            for i in range(len(lines)):
                line = lines[i]
                record = line[:6].strip()
                if record == 'ATOM':
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_no = int(line[22:26])
                    element = line[76:78].strip()
                    if res_name in d_titratable.keys():
                        if res_name == 'CYS':
                            res_name = 'CYS2'
                        elif res_no-1 in l_seqIDs:
                            res_name = d_titratable[res_name]['nondefault']
                        else:
                            res_name = d_titratable[res_name]['default']
                    if res_no == 1:
                        res_name = 'N'+res_name
                    if res_no == len(l_seqres):
                        res_name = 'C'+res_name
                        if atom_name in d_terminal.keys():
                            atom_name = d_terminal[atom_name]
                    if len(atom_name) < 4:
                        atom_name = '%2s%2s' %(element.rjust(2),atom_name[len(element):].ljust(2),)
                    line2 = '%12s%4s%1s%4s%59s\n' %(line[:12],atom_name,line[16:17],res_name.ljust(4),line[21:80])
                    lines2 += [line2]
                elif record == 'ANISOU':
                    line2 = '%12s%4s%1s%4s%59s\n' %(line[:12],atom_name,line[16:17],res_name.ljust(4),line[21:80])
                    lines2 += [line2]
                else:
                    lines2 += [line]

            fd = open('%s.pdb' %(pdb),'w')
            fd.writelines(lines2)
            fd.close()

        return charge, s_pdb2gmx


    def __init__(self):

        import optparse

##        parser = optparse.OptionParser()
##
##        parser.add_option(
##            '--posre',
##            dest='posre',
##            action="store_true",
##            default=False,
##            )
##
##        (options,args,) = parser.parse_args()
##
##        self.posre = options.posre

        if not '-pdb' in sys.argv:
            print 'Which protein do you want to use as input? Use flag -pdb'
            self.pdb = raw_input()
        else:
            self.pdb = sys.argv[sys.argv.index('-pdb')+1]

##        if not os.path.isfile('%s.pdb' %(self.pdb)):
##            print 'pdb %s not present in the current directory' %(self.pdb)

        if '-cluster' in sys.argv:
            os.system('cp /home/tc/pdb/%s.pdb %s.pdb' %(self.pdb[:4],self.pdb[:4],))

##            fd = open('%s.pdb' %(self.pdb),'r')
##            s = fd.read()
##            fd.close()
##            ## manual hydrogen bond optimization
##            s = s.replace(
##                'ATOM    714  OE1 GLU A  35       5.327  15.197  26.471  1.00  4.88           O',
##                'ATOM    714  OE2 GLU A  35       5.327  15.197  26.471  1.00  4.88           O',
##                )
##            s = s.replace(
##                'ATOM    715  OE2 GLU A  35       4.361  17.082  27.194  1.00  5.02           O',
##                'ATOM    715  OE1 GLU A  35       4.361  17.082  27.194  1.00  5.02           O',
##                )
####            s = s.replace(
####                'ATOM   1023  OD1 ASP A  52       5.098  10.151  31.543  1.00  5.85           O',
####                'ATOM   1023  OD2 ASP A  52       5.098  10.151  31.543  1.00  5.85           O',
####                )
####            s = s.replace(
####                'ATOM   1024  OD2 ASP A  52       6.886  10.645  30.401  1.00  6.57           O',
####                'ATOM   1024  OD1 ASP A  52       6.886  10.645  30.401  1.00  6.57           O',
####                )
##            if (
##                '52' in sys.argv[sys.argv.index('-nondefaultcharge')+1]
##                ):
####                ## delete hydrogen atoms (gromacs will add them again)
####                s = s.replace('ATOM   1158 HD21 ASN A  59       4.340   9.538  33.210  1.00  4.35           H  \n','',)
####                s = s.replace('ATOM   1159 HD22 ASN A  59       4.597   8.965  34.561  1.00  4.35           H  \n','',)
####                s = s.replace(
####                    'ATOM   1153  ND2 ASN A  59       4.070   9.360  34.007  1.00  3.63           N',
####                    'ATOM   1153  OD1 ASN A  59       4.070   9.360  34.007  1.00  3.63           O',
####                    )
####                s = s.replace(
####                    'ATOM   1152  OD1 ASN A  59       2.427   9.449  35.530  1.00  3.47           O',
####                    'ATOM   1152  ND2 ASN A  59       2.427   9.449  35.530  1.00  3.47           N',
####                    )
##                s = s.replace(
##                    'ATOM    889  OD1 ASN A  44      11.295   9.276  28.295  1.00  7.22           O',
##                    'ATOM    889  ND2 ASN A  44      11.295   9.276  28.295  1.00  7.22           N',
##                    )
##                s = s.replace(
##                    'ATOM    890  ND2 ASN A  44       9.972   9.708  30.042  1.00 10.35           N',
##                    'ATOM    890  OD1 ASN A  44       9.972   9.708  30.042  1.00 10.35           O',
##                    )
##            fd = open('%s.pdb' %(self.pdb),'w')
##            fd.write(s)
##            fd.close()

        else:
            os.system('cp /data/pdb-v3.2/%s/pdb%s.ent %s.pdb' %(self.pdb[1:3],self.pdb[:4],self.pdb[:4],))

        if not '-t' in sys.argv:
            print 'For how many nanoseconds do you want to run your simulation? Use flag -t'
            self.t = float(raw_input())
        else:
            self.t = float(sys.argv[sys.argv.index('-t')+1])

        if not '-T' in sys.argv:
            print 'At what Kelvin temperature do you want to run your simulation? Use flag -T'
            self.T = float(raw_input())
        else:
            self.T = float(sys.argv[sys.argv.index('-T')+1])

        self.posre = False
        if '-posre' in sys.argv:
            self.posre = True

        self.parallel = False
        if '-parallel' in sys.argv:
            self.parallel = True

        self.arguments = sys.argv

        return


if __name__=='__main__':

    instance_gromacs = gromacs()
    instance_gromacs.main()

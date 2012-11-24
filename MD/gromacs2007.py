#!/usr/bin/env python
#$Id: gromacs.py 239 2007-05-16 14:16:02Z tc $

## double precision ?

import os, sys

class gromacs:

    def main(self):

        pdb = self.pdb[:-4]
        index_pre_hydrogen = 'crystalcontacts.ndx'
        index = 'crystalcontacts_pdb2gmx.ndx' ## post hydrogen (pdb2gmx) index file

        ## calculate charge of protein by counting charged residues
        charge = self.calculate_charge(pdb,)

        ##
        ## gromacs initialization
        ##
        
        ## add hydrogen atoms and generate a gromacs topology file
        self.pdb2gmx(pdb, ff=4)
        ## enlarge the box to accommodate for water (and center protein)
        self.editconf(pdb,)
        ## generate the box and fill it with water
        ## i.e. solvate protein and change topology file to include water
        self.genbox(pdb)

        ##
        ## energy minimization (minimize repulsive contacts between protein and solvent)
        ##

        ## generate mdrun input file (binary topology file .tpr)
        ## i.e. preprocess topology (.top), structure (.gro), parameters (.mdp)
        self.grompp(pdb, runtype='EM', i_suffix='genbox', index=index,)

        if charge != 0:
            ## generate ions (i.e. generate net-zero total charge for the system)
            self.genion(pdb, charge,)

        ## preprocess generic structure from genion
##        self.grompp(pdb, runtype='EM', i_suffix='genion', index=index,)
        self.grompp(pdb, runtype='EM', i_suffix='genion',)
        ## energy minimize protein and solvent
        ## i.e. minimize repulsive contacts between protein and solvent
        self.mdrun(pdb, runtype='EM')

        ##
        ## position restrained dynamics ("fill" protein cavities with water)
        ##
        self.grompp(pdb, runtype='PR', i_suffix='EM', index=index,)
        self.mdrun(pdb, runtype='PR')

        ##
        ## molecular dynamics
        ##
        self.grompp(pdb, runtype='MD', i_suffix='PR', index=index,) 
        self.mdrun(pdb, runtype='MD')


    def genion(self, pdb, charge,):

        source = '%s_genion.src' %(pdb)

        if charge > 0:
            chargeflag = 'nn'
            iontype = 'CL-'
        elif charge < 0:
            chargeflag = 'np'
            iontype = 'NA+'

        lines = [
            'genion ',
            '-s %s_EM.tpr ' %(pdb),
            '-o %s_genion.gro ' %(pdb),
            '-%s %s ' %(chargeflag, abs(charge)),
            '-g %s_genion.log ' %(pdb),
            ]
        lines += ['<< EOF\n12\n',]

        fd=open(source,'w')
        fd.writelines(lines)
        fd.close()

        os.system('source %s' %(source))


        fd = open('%s.top' %(pdb),'r')
        lines = fd.readlines()
        fd.close()
        
        line_SOL = lines[-1]
        print line_SOL
        SOL_mols = int(line_SOL.split()[1])-abs(charge)
        line_sol = 'SOL     %13i\n' %(SOL_mols)

        line_ion = '%3s     %13i\n' %(iontype, abs(charge))

        lines = lines[:-1]+[line_sol]+[line_ion]

        fd = open('%s.top' %(pdb),'w')
        fd.writelines(lines)
        fd.close()


    def mdrun(self, pdb, runtype=''):

        source = '%s_%smdrun.src' %(pdb, runtype)

        fd=open(source,'w')
        fd.writelines([
            'mdrun ',
            '-s %s_%s.tpr ' %(pdb, runtype), ## input, generic run input
            '-o %s_%s.trr ' %(pdb, runtype), ## output, trajectory
            '-c %s_%s.gro ' %(pdb, runtype), ## output, generic structure
            '-e %s_%sener.edr ' %(pdb, runtype), ## output, generic energy
            '-g %s_%smdrun.log ' %(pdb, runtype), ## output, log
            '-v ',
            '<< EOF\n0',
            ])
        fd.close()

        os.system('source %s' %(source))

    def grompp(self, pdb, runtype, i_suffix, index=None):

        o_suffix = runtype

        if runtype == 'EM':

            d_mdp = {
                'cpp':'/lib/cpp',
                'define':'-DFLEX_SPC',
                'constraints':'none',
                'integrator':'steep',
                'nsteps':10000.,
                ## Energy minimizing stuff
                'emtol':100.,
                'emstep':.01,

                'nstcomm':1,
                'ns_type':'grid',
                'rlist':1,
                'rcoulomb':1.,
                'rvdw':1.,
                'Tcoupl':'no',
                'Pcoupl':'no',
                'gen_vel':'no',
                }

        if runtype == 'PR':

            d_mdp = {
                'cpp':'/lib/cpp',
                'define':'-DPOSRES',
                'constraints':'all-bonds',
                'integrator':'md',
                'dt':.002, ## ps timestep
                'nsteps':50000,
                'nstcomm':1,
                'nstxout':50,
                'nstvout':1000,
                'nstfout':0,
                'nstlog':10,
                'nstenergy':10,
                'nstlist':10,
                'ns_type':'grid',
                'rlist':1.0,
                'rcoulomb':1.0,
                'rvdw':1.0,
                ## Berendsen temperature coupling is on in two groups
                'Tcoupl':'berendsen',
                'tc-grps':'Protein OTHER',
                'tau_t':'.1 .1',
                'ref_t':'300 300',
                ## Energy monitoring
                'energygrps':'Protein OTHER',
                ## Pressure coupling is not on
                'Pcoupl':'no',
                'tau_p':.5,
                'compressibility':'4.5e-5',
                'ref_p':1.,
                }

        if runtype == 'MD':

            d_mdp = {
                'cpp':'/lib/cpp',
                'constraints':'all-bonds',
                'integrator':'md',
                'dt':.002, ## ps
                'nsteps':1000000,
                'nstcomm':1,
                'nstxout':1000,
                'nstvout':1000,
                'nstfout':0,
                'nstlog':100,
                'nstenergy':10,
                'nstlist':10,
                'ns_type':'grid',
                'rlist':1.,
                'rcoulomb':1.,
                'rvdw':1.,
                ## Berendsen temperature coupling is on in two groups
                'Tcoupl':'berendsen',
                'tc-grps':'Protein OTHER',
                'tau_t':'.1 .1',
                'ref_t':'300 300',
                ## Energy monitoring
                'energygrps':'Protein OTHER',
                ## Isotropic pressure coupling is now on
                'Pcoupl':'berendsen',
                'Pcoupltype':'isotropic',
                'tau_p':.5,
                'compressibility':'4.5e-5',
                'ref_p':1.,
                ## Generate velocites is off at 300 K.
                'gen_vel':'no',
                'gen_temp':300.,
                'gen_seed':173529,
                }

        ## freeze groups
        if index != None:
            d_mdp['freezegrps'] = 'crystalcontacts'
            d_mdp['freezedim'] = 'Y Y Y'

        mdp = '%s_%sgrompp.mdp' %(pdb, runtype)
        fd=open(mdp,'w')
        for key in d_mdp.keys():
            fd.write('%s = %s\n' %(key,d_mdp[key]))
        fd.close()

        source = '%s_%sgrompp.src' %(pdb, i_suffix)

        lines = [
            'grompp ',
            '-f %s ' %(mdp), ## input, MD parameters
            '-po %s ' %(mdp), ## output, MD parameters
            '-c %s_%s.gro ' %(pdb, i_suffix), ## input, generic structure
            '-p %s.top ' %(pdb), ## input, topology file
            '-o %s_%s.tpr ' %(pdb, o_suffix), ## output, generic run input
            ]
        if index != None:
            lines += [
                ' -n %s' %(index),
                ]

        fd = open(source,'w')
        fd.writelines(lines)
        fd.close()

        os.system('source %s' %(source))

        return

    def genbox(self, pdb):

        import os

        source = '%s_genbox.src' %(pdb)

        fd = open(source,'w')
        fd.writelines([
            'genbox ',
            '-cp %s_editconf.gro ' %(pdb), ## input, generic structure (protein)
            '-cs ', ## input, generic structure (solvent)
            '-o %s_genbox.gro ' %(pdb), ## output, generic structure
            '-p %s.top' %(pdb), ## topology file
            ]),
        fd.close()

        os.system('source %s' %(source))

        return        

    def editconf(self, pdb,):

        import os

        source = '%s_editconf.src' %(pdb)

        s = 'editconf -f %s_pdb2gmx.gro -o %s_editconf.gro -d 1.0 -bt cubic' %(pdb, pdb) ## -d implies -c

        fd = open(source,'w')
        fd.write(s)
        fd.close()

        os.system('source %s' %(source))

        return

    def pdb2gmx(self, pdb, ff=4):

        import os

        source = '%s_pdb2gmx.src' %(pdb)

        fd = open(source,'w')
        fd.writelines([
            'pdb2gmx ',
            '-f %s.pdb ' %(pdb), ## input, generic structure
            '-o %s_pdb2gmx.gro ' %(pdb), ## output, generic structure
            '-p %s.top ' %(pdb), ## output, topology file
            '-i %s.itp ' %(pdb), ## output, include file for topology
            '<<EOF\n%s\n' %(ff), ## force field selection
            ])
        fd.close()

        os.system('source %s' %(source))

        return


    def calculate_charge(self,pdb,):

        d_res = {}

        fd = open('%s.pdb' %(pdb),'r')
        lines = fd.readlines()
        fd.close()

        for line in lines:

            if line[:4] == 'ATOM':
                chain = line[21]
                res_name = line[17:20]
                res_no =  int(line[23:26])
                d_res['%s:%s' %(chain,res_no)] = res_name

        charge = 0
        for res_no in d_res.keys():
            res_name = d_res[res_no]
            if res_name in ['GLU','ASP']:
                charge -= 1
            elif res_name in ['ARG','LYS']:
                charge += 1

        return charge


    def __init__(self):

        if not '-pdb' in sys.argv:
            raise 'use -pdb flag to select pdb'
        self.pdb = sys.argv[sys.argv.index('-pdb')+1]
        if not os.path.isfile(self.pdb):
            raise 'pdb %s not present in the current directory' %(self.pdb)

if __name__=='__main__':

    instance_gromacs = gromacs()
    instance_gromacs.main()

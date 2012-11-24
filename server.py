#!/bin/env /software/bin/python
#
# $Id: server.py,v 1.6 2006/02/14 12:50:21 tc Exp $

class server:

    def main(self):
        import cgi, cgitb, time, os, sys, shutil
        cgitb.enable() ## enable traceback, which reports non-syntax errors


        ## set jobid
        timetuple = time.gmtime(time.time())
        timestr = str(timetuple[0])+str(timetuple[1]).zfill(2)+str(timetuple[2]).zfill(2)+str(timetuple[3]).zfill(2)+str(timetuple[4]).zfill(2)+str(timetuple[5]).zfill(2)

        form = cgi.FieldStorage()

        if form.has_key("timestr"):
            timestr = form["timestr"].value


        ## set paths
        self.path_tmp = '/var/www/cgi-bin/goodvibes/tmp/'
        self.path_results = '/var/www/html/goodvibes/results/'
        self.path_pdbs = '/data/pdb/'
        self.path_python = '/var/www/cgi-bin/goodvibes/python/'


        ## check form for erros
        if form["inputfiles"].value == '1':
            inputfiles = 1
            self.errorcheck1file(form)
        elif form["inputfiles"].value == '2':
            inputfiles = 2
            self.errorcheck2files(form)


        ## if no errors then queue job
        if os.path.isfile(self.path_tmp+'queue.txt'):
            fd = open(self.path_tmp+'queue.txt', 'a')
            fd.write(timestr+'\n')
            fd.close()
            queuemanager = False
        else:
            fd = open(self.path_tmp+'queue.txt', 'w')
            fd.write(timestr+'\n')
            fd.close()
            queuemanager = True


        ## write fast output and slow 1st temporary output
##        pid = os.fork()
        pid = 11
        if pid:
            htmlbody = 'Dear user. Your input data has been queued for processing. Your job is still in the queu.'
            self.slowhtml(htmlbody, timestr)

            htmlbody = 'Dear user. Your input data is being processed. Your results will be available from <a href="http://polymerase.ucd.ie/goodvibes/results/'+timestr+'.html">http://polymerase.ucd.ie/goodvibes/results/'+timestr+'.html</a> in a minute...'
            self.fasthtml(htmlbody)

##            sys.stdout.flush()
##            os.close(0)
##        
##        os.close(1)
##        os.close(2)

################################################################################

        ## parse html form

        if form.has_key("cutoffs"):
            cutoffs = form.getlist("cutoffs")
            for i in range(len(cutoffs)):
                cutoffs[i] = float(cutoffs[i])
        else: cutoffs = [10]
        if form.has_key("amplitude_average"): amplitude_average = form["amplitude"].value
        else: amplitude_average = 2
        if form.has_key("mutations"): mutations = form["mutations"].value
        else: mutations = 0
        if form.has_key("quarternary"): quarternary = form["quarternary"].value
        else: quarternary = 'monomeric'
        if form.has_key("frames"): frames = form['frames'].value
        else: frames = 50
        if form.has_key("pdb1_id"):
            fd = open(self.path_pdbs+form["pdb1_id"].value.lower()+'.pdb', 'r')
            pdb_lines = fd.readlines()
            fd.close()
            fd = open(self.path_tmp+timestr+'_reference.pdb', 'w')
            fd.writelines(pdb_lines)
            fd.close()
        if form.has_key("pdb1_file"):
            fd = open(self.path_tmp+timestr+'_reference.pdb', 'w')
            fd.write(form["pdb1_file"].value)
            fd.close()
        if form.has_key("chains1"): chains1 = form["chains1"].value.upper()
        else: chains1 = ''
        ## resranges need to refer to chains. resranges dont need to be updated when building biomolecules from remark350 because residues not parsed.
##        self.residue_terminal_n = {'A':-9999, 'B':-9999, 'C':-9999, 'D':-9999, 'E':-9999, 'F':-9999, 'G':-9999, 'H':-9999, 'I':-9999, 'J':-9999, 'K':-9999, 'L':-9999, 'M':-9999, 'N':-9999, 'O':-9999, 'P':-9999, 'Q':-9999, 'R':-9999, 'S':-9999, 'T':-9999, 'U':-9999, 'V':-9999, 'W':-9999, 'X':-9999, 'Y':-9999, 'Z':-9999, ' ':-9999,}
##        self.residue_terminal_c = {'A':9999, 'B':9999, 'C':9999, 'D':9999, 'E':9999, 'F':9999, 'G':9999, 'H':9999, 'I':9999, 'J':9999, 'K':9999, 'L':9999, 'M':9999, 'N':9999, 'O':9999, 'P':9999, 'Q':9999, 'R':9999, 'S':9999, 'T':9999, 'U':9999, 'V':9999, 'W':9999, 'X':9999, 'Y':9999, 'Z':9999, ' ':9999,}
        if form.has_key("resrange1start"): resrange1start = form["resrange1start"].value
        else: resrange1start = -99999
        if form.has_key("resrange1end"): resrange1end = form["resrange1end"].value
        else: resrange1end = 99999
        if form.has_key("model1"): model1 = form["model1"].value
        else: model1 = ''

        if not form.has_key("atoms"):
            atoms_hessian = ['CA']
        else:
            l_atoms = form.getlist("atoms")
            import sets
            all = sets.Set(['OXT'])
            for chemical_symbol in ['C','O','H','N','S']:
                all.add(chemical_symbol)
                for remoteness_indicator in ['A','B','G','D','E','Z','H']:
                    all.add(chemical_symbol+remoteness_indicator)
                    for branch_designator in range(10):
                        all.add(chemical_symbol+remoteness_indicator+str(branch_designator))
            mainchain = sets.Set(['C', 'O', 'N', 'CA'])
            sidechain = all-mainchain
            alpha = sets.Set(['CA'])
            beta = sets.Set(['CB'])

            atoms_hessian = sets.Set()
            for s_atoms in l_atoms:
                if s_atoms == 'yall':
                    atoms_hessian |= all
                    break
                if s_atoms == 'ymainchain':
                    atoms_hessian |= mainchain
                if s_atoms == 'ysidechain':
                    atoms_hessian |= sidechain
                if s_atoms == 'yalpha':
                    atoms_hessian |= alpha
                if s_atoms == 'ybeta':
                    atoms_hessian |= beta
            for s_atoms in l_atoms:
                if s_atoms == 'nmainchain':
                    atoms_hessian -= mainchain
                if s_atoms == 'nsidechain':
                    atoms_hessian -= sidechain
                if s_atoms == 'nalpha':
                    atoms_hessian -= alpha
                if s_atoms == 'nbeta':
                    atoms_hessian -= beta
            atoms_hessian = list(atoms_hessian)

        if form.has_key("pdb2_id"):
            fd = open(self.path_pdbs+form["pdb2_id"].value.lower()+'.pdb', 'r')
            pdb_lines = fd.readlines()
            fd.close()
            fd = open(self.path_tmp+timestr+'_conformer.pdb', 'w')
            fd.writelines(pdb_lines)
            fd.close()

        if form.has_key("pdb2_file"):
            fd = open(self.path_tmp+timestr+'_conformer.pdb', 'w')
            fd.write(form["pdb2_file"].value)
            fd.close()

        if form.has_key("chains2"): chains2 = form["chains2"].value.upper()
        else: chains2 = ''

        if form.has_key("resrange2start"): resrange2start = form["resrange2start"].value
        else: resrange2start = -99999

        if form.has_key("resrange2end"): resrange2end = form["resrange2end"].value
        else: resrange2end = 99999

        if form.has_key("model2"): model2 = form["model2"].value
        else: model2 = ''

        ## write job input to job file
        fd = open(self.path_tmp+'job_'+timestr+'.txt', 'w')
        fd.write(
            str(chains1)+'\n'+str(model1)+'\n'+str(resrange1start)+'\n'+str(resrange1end)+'\n'+
            str(cutoffs)+'\n'+str(amplitude_average)+'\n'+str(atoms_hessian)+'\n'+str(mutations)+'\n'+str(quarternary)+'\n'+
            str(frames)+'\n'+str(inputfiles)+'\n'+
            str(chains2)+'\n'+str(model2)+'\n'+str(resrange2start)+'\n'+str(resrange2end)+'\n'
            )
        fd.close()

        del timestr

################################################################################

        ## end if not queue manager
        if not queuemanager:
            return

################################################################################

        ## process queue if queue manager

        ## change dir in case of any goodvibes outputs (dislin etc.)
        os.chdir(self.path_results)
        ## append the python library to sys.paths after change of dir
        sys.path.append(self.path_python)
        import goodvibes
        instance_NMA = goodvibes.vibration()

        ## process queue while it exists
        while os.path.isfile(self.path_tmp+'queue.txt'):

            ## parse job id
            fd = open(self.path_tmp+'queue.txt', 'r')
            job = fd.readline().strip()
            fd.close()

            ## write slow 2nd temporary output
            htmlbody = 'Dear user. Your input is being processed. Your results - or an error if you gave an invalid input - will be here in less than 30 minutes.'
            self.slowhtml(htmlbody, job)

            ## parse parameters from job txt file to goodvibes
            fd = open(self.path_tmp+'job_'+job+'.txt', 'r')
            lines_job = fd.readlines()
            fd.close()

            if lines_job[0].strip() != '': chains1 = lines_job[0].strip().split(',')
            else: chains1 = []
            model1 = lines_job[1].strip()
            resrange1start = int(lines_job[2].strip())
            resrange1end = int(lines_job[3].strip())
            cutoffs = lines_job[4].strip()[1:-1].split(', ')
            for i in range(len(cutoffs)):
                cutoffs[i] = float(cutoffs[i])
            amplitude_average = float(lines_job[5].strip())
            atoms_hessian = lines_job[6].strip()[2:-2].split(',')
            mutations = int(lines_job[7].strip())
            quarternary = lines_job[8].strip()
            frames = float(lines_job[9].strip())
            inputfiles = int(lines_job[10].strip())

            fd = open(self.path_tmp+job+'_reference.pdb', 'r')
            pdblines1 = fd.readlines()
            fd.close()

            if inputfiles == 2:
                fd = open(self.path_tmp+job+'_conformer.pdb', 'r')
                pdblines2 = fd.readlines()
                fd.close()
                if lines_job[11].strip() != '': chains2 = lines_job[11].strip().split(',')
                else: chains2 = []
                model2 = lines_job[12].strip()
                resrange2start = int(lines_job[13].strip())
                resrange2end = int(lines_job[14].strip())
            else:
                pdblines2 = None
                chains2 = None
                model2 = None

            ## run goodvibes
            try:
                results = instance_NMA.main(pdblines1,chains1,model1, atoms_hessian, amplitude_average, cutoffs, mutations, quarternary, job, frames, pdblines2, chains2, model2)
            except 'chainerror', error:
                ## if expected error then remove from queue and report to user
                self.slowhtml(error, job)
                self.queue()
                continue
            except:
                ## if not expected error (not syntax errors) then remove from queue and report to user
                self.slowhtml('there was a bug. please email your input to <a href="mailto:tommy.carstensen@ucd.ie?subject=GoodVibes Server&body=Hi Tommy,">administrator</a> :-)', job)
                self.queue()
                results = instance_NMA.main(pdblines1,chains1,model1, atoms_hessian, amplitude_average, cutoffs, mutations, quarternary, job, frames, pdblines2, chains2, model2)
                continue

            ## write html if no exceptions
            html = '''
            <html>
            <head><title>GoodVibes results for 2 structures</title></head>
            <body>
            <table border="1">
            <tr>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">NUMBER OF NON-ZERO FREQUENCY MODES</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">DEPOSIT1</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">CUT-OFF DISTANCE</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">RESOLUTION1</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">CHAINS IN PDB1</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">VMD PDB MOVIE MODE 7</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">VMD PDB MOVIE MODE 8</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">VMD PDB MOVIE MODE 9</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">VMD PDB MOVIE MODE 10</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">VMD PDB MOVIE MODE 11</div></td>
            <td valign="bottom" align="left"><div style="writing-mode:tb-rl">VMD PDB MOVIE MODE 12</div></td>
            '''
            if inputfiles == 2:
                html += (
                    '\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">RMSD</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">DEPOSIT2</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">RESOLUTION2</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">OVERLAP OF MODE 7x12</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">OVERLAP OF MODE 8x7</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">OVERLAP OF MODE 9x12</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">OVERLAP OF MODE 10x7</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">OVERLAP OF MODE 11x12</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">OVERLAP OF MODE 12x7</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">EQUIVALENT CHAINS IN PDB1 AND PDB2</div></td>\
                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">CHAINS IN PDB2</div></td>\
                    '
                    )
##                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">MODES OF MAX OVERLAP</div></td>\
##                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">MAX OVERLAPS</div></td>\
##                    <td valign="bottom" align="left"><div style="writing-mode:tb-rl">EIGENVALUES OF MODE OF MAX OVERLAP</div></td>\
            html += '</tr>'

            for result in results:
                chains = ''.join(result['biomolecules'][0]).replace(' ','-')
                file_prefix = job+'_'+chains+'_'+str(result['cutoff'])+'_'
                file_vmd = job+'_'+chains+'_'+str(result['cutoff'])+'_mode7.pdb'
##                print result['number of non-zero eigenvalues'], result['mode of max overlap'], result['max overlap'],
##                print result['eigenvalue of max overlap'], result['rmsd'], result['depdate1'], result['depdate2'],
##                print result['cutoff'], result['res1'], result['res2'],
##                print result['overlap7'], result['overlap8'], result['overlap9'], result['overlap10'], result['overlap11'], result['overlap12'],
##                print result['biomolecules'], result['chains1'], result['chains2'],
##                print file_vmd,
                html += (
                    '\n<tr>\
                    \n<td align="right">%4i</td>\
                    \n<td align="left">%9s</td>\
                    \n<td align="right">%4.1f</td>\
                    \n<td align="right">%4.2f</td>\
                    \n<td align="left">%s</td>\
                    '
                    %(
                        result['number of non-zero eigenvalues'],  
                        result['depdate1'],
                        result['cutoff'], result['res1'],
                        result['chains1'],
                        )
                    )

                for mode in range(6,12):
                    html += '\n<td align="left"><a href="%s">mode %s</a></td>' %(file_prefix+str(mode+1)+'.pdb', str(mode+1))

                if inputfiles == 2:
                    html += (
                        '\
                        \n<td align="right">%3i</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%4.1f</td>\
                        \n<td align="left">%9s</td>\
                        \n<td align="right">%4.2f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="right">%5.3f</td>\
                        \n<td align="left">%s</td>\
                        \n<td align="left">%s</td>\
                        '
                        %(
                            result['mode of max overlap'],
                            result['max overlap'],
                            result['eigenvalue of max overlap'],
                            result['rmsd'],
                            result['depdate2'],
                            result['res2'],
                            result['overlap7'], result['overlap8'], result['overlap9'], result['overlap10'], result['overlap11'], result['overlap12'],
                            result['biomolecules'],
                            result['chains2'],
                            )
                        )

                html += '''
                </tr>
                </table>
                <br>
                '''

                ## vmd output
                for mode in range(6,12):
                    output_vmd = ['REMARK color by connectivity (b-factor) or squared displacement (temperature factor)\n']
                    for frame in range(len(result['morph lines'][mode-6])):
                        output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
                        output_vmd += result['morph lines'][mode-6][frame]
                        output_vmd.append('TER\nENDMDL\n')
                    fd = open(self.path_results+file_prefix+str(mode+1)+'.pdb', 'w') ## implicit path
                    fd.writelines(output_vmd)
                    fd.close()

                    ## rasmol input/output
                    for frame in range(len(result['morph lines'][mode-6])):
                        rasmol_pdb = job+'_'+str(frame)+'.pdb'
                        rasmol_script = self.path_tmp+'rasmol.script'
                        rasmol_img = job+'_'+str(frame)+'.ppm'
                        
                        fd = open(rasmol_pdb, 'w')
                        fd.writelines(result['morph lines'][mode-6][frame])
                        fd.close()
                        
                        fd = open(rasmol_script,'w')
                        fd.write('/software/bin/rasmol -nodisplay '+rasmol_pdb+' << EOF\n')
                        fd.write('color temperature\n') ## color by displacement
                        ## somehow center at center of wildtype... by translation only! no rotation!
                        fd.write('write '+self.path_tmp+rasmol_img+'\n')
                        fd.write('exit\n')
                        fd.close()
                        
                        os.system('source '+rasmol_script+' > '+self.path_tmp+'rasmol.log') ## write image
                        os.system('/usr/bin/convert '+self.path_tmp+rasmol_img+' -resize 15% '+self.path_tmp+job+'_'+str(frame)+'.gif') ## resize image
        
                    ImageMagick_script = '/usr/bin/convert -loop 0 -delay 2x'+str(len(result['morph lines'][mode-6])-2) ## -delay fps (sxf)
                    for frame in range(len(result['morph lines'][mode-6])):
                        rasmol_gif = job+'_'+str(frame)+'.gif'
                        ImageMagick_script += ' '+self.path_tmp+rasmol_gif
                    for frame in range(len(result['morph lines'][mode-6])-1-1,-1+1,-1):
                        rasmol_gif = job+'_'+str(frame)+'.gif'
                        ImageMagick_script += ' '+self.path_tmp+rasmol_gif
                    ImageMagick_script += ' '+self.path_results+file_prefix+str(mode+1)+'.gif'
                    os.system(ImageMagick_script)
                    html += '\n<img src="'+file_prefix+'%s.gif" alt="mode%s">' %(mode+1, mode+1)

                html += '\n<br><br>'

##                    fd = open(self.path_tmp+'ImageMagick.script', 'w')
##                    fd.write('/usr/bin/convert -loop 0 -delay 2x'+str(len(result['morph lines'][mode-6])-1)) ## -delay fps (sxf)
##                    for frame in range(len(result['morph lines'][mode-6])):
##                        rasmol_gif = job+'_'+str(frame)+'.gif'
##                        fd.write(' '+self.path_tmp+rasmol_gif)
##                    for frame in range(len(result['morph lines'][mode-6])-1-1,-1+1,-1):
##                        rasmol_gif = job+'_'+str(frame)+'.gif'
##                        fd.write(' '+self.path_tmp+rasmol_gif)
##                    fd.write(' '+self.path_results+file_prefix+str(mode+1)+'.gif << EOF\n')
##                    fd.close()
##                    
##                    os.system('source '+self.path_tmp+'ImageMagick.script'+' > '+self.path_tmp+'ImageMagick.log')

                ## image html output
                img_files = []
                for mode in range(6,12):
                    img_files.append('2D_atom_displacement'+str(mode+1)+'.png')
                    img_files.append('crosscorrelation'+str(mode+1)+'.png')
                if inputfiles == 2:
                    img_files += [
                    '2D_mode_overlap.png',
                    ]

                for filesuffix in img_files:
                    html += '\n<img src="'+file_prefix+filesuffix+'"><br><br>'

            if inputfiles == 1:
                html += (
                    '\n<br><hr><br>\
                    \n<form name="input" action="http://polymerase.ucd.ie/cgi-bin/goodvibes/python/goodvibes/server.py" method="post">\
                    \n<fieldset>\
                    \n<legend><h2>Secondary <b>sequence-identical</b> structure</h2></legend>\
                    \nEnter <label for="pdb2_id"><b>PDB ID</b></label>\
                    \n<input type="text" name="pdb2_id" value="" size="4" id="pdb2_id">\
                    \nor upload <label for="pdb2_file"><b>PDB file</b></label>\
                    \n<input type="file" name="pdb2_file" id="pdb2_file">\
                    \n<br><br>Enter <label for="chains2"><b>Chain ID(s)</b></label>\
                    \n<input type="text" name="chains2" id="chains2">\
                    \nseperated by commas.\
                    \n</fieldset>\
                    \n<input type="hidden" name="inputfiles" value="2">\
                    \n<input type="hidden" name="timestr" value="'+str(job)+'">\
                    \n<input type="hidden" name="chains1" value="'+lines_job[0]+'">\
                    \n<input type="hidden" name="model1" value="'+lines_job[1]+'">\
                    \n<input type="hidden" name="resrange1start" value="'+lines_job[2]+'">\
                    \n<input type="hidden" name="resrange1end" value="'+lines_job[3]+'">\
                    \n<br><input type="submit" value="Submit data and choices"><input type="reset" value="Reset all form fields">\
                    \n</form>\
                    '
                    )

            html += '''
            </body>
            </html>
            '''

            ## write the html results
            fd = open(self.path_results+job+'.html', 'w') ## implicit path
            fd.write(html)
            fd.close()


            ## remove the job from the queue and continue
            self.queue()
            continue

        ## return when no more jobs to process in queue
        return

    def errorcheck2files(self, form):

        import os
        
##        ## check for essential settings
##        if not (form.has_key("cutoffs") and form.has_key("amplitude") and
##                form.has_key("mutations") and form.has_key("atoms") and form.has_key("frames")):
##                self.fasthtml('Dear User. You need to set cutoff distance, amplitude, number of allowed mutations, frames and atoms to be used.', error = True)

        ## check for essential input
        if not (form.has_key("pdb2_id") or form.has_key("pdb2_file")):
            self.fasthtml('Dear user. You need to specify a pdb to be used.', error = True)

        ## check for valid pdb ids
        if form.has_key("pdb2_id"):
            if not os.path.isfile(self.path_pdbs+form["pdb2_id"].value.lower()+'.pdb'):
                self.fasthtml('The PDB ID you have specified is not found in the protein data bank.', error = True)

        ## check for identical pdb ids
        if form.has_key("pdb1_id") and form.has_key("pdb2_id"):
            pdb1_id = form["pdb1_id"].value.upper()
            pdb2_id = form["pdb2_id"].value.upper()
            if pdb1_id == pdb2_id:
                if form.has_key("chains1"): chains1 = form["chains1"].value.upper()
                else: chains1 = ''
                if form.has_key("chains2"): chains2 = form["chains2"].value.upper()
                else: chains2 = ''
                if chains1 == chains2:
                    self.fasthtml('Dear user. You used the same pdb id as reference structure and secondary structure.', error = True)

        ## check for valid number of frames
        if form.has_key("frames"):
            if int(form["frames"].value) > 500 or int(form["frames"].value) < 2:
                self.fasthtml('Please change the number of frames to be within the range 10 to 500.', error = True)

        ## acummulate error messages and then return them in errorhtml if there is any...

        return

    def errorcheck1file(self, form):

        import os
        
##        ## check for essential settings
##        if not (form.has_key("cutoffs") and form.has_key("amplitude") and
##                form.has_key("mutations") and form.has_key("atoms") and form.has_key("frames")):
##                self.fasthtml('Dear User. You need to set cutoff distance, amplitude, number of allowed mutations, frames and atoms to be used.', error = True)

        ## check for essential input
        if not (form.has_key("pdb1_id") or form.has_key("pdb1_file")):
            self.fasthtml('Dear user. You need to specify a pdb to be used.', error = True)

        ## check for valid pdb ids
        if form.has_key("pdb1_id"):
            if not os.path.isfile(self.path_pdbs+form["pdb1_id"].value.lower()+'.pdb'):
                self.fasthtml('The PDB ID you have specified is not found in the protein data bank.', error = True)

        ## check for valid number of frames
        if form.has_key("frames"):
            if int(form["frames"].value) > 500 or int(form["frames"].value) < 2:
                self.fasthtml('Please change the number of frames to be within the range 10 to 500.', error = True)

        ## acummulate error messages and then return them in errorhtml if there is any...

        return
    
    def queue(self):

        import os
        
        ## remove job from queue or remove queue
        fd = open(self.path_tmp+'queue.txt', 'r')
        lines_queue = fd.readlines()
        fd.close()
        if len(lines_queue) > 1:
            fd = open(self.path_tmp+'queue.txt', 'w')
            fd.writelines(lines_queue[1:])
            fd.close()
        else:
            os.remove(self.path_tmp+'queue.txt')
        
    def fasthtml(self, htmlbody, error = False):

        import sys, os
        
        print 'Content-type: text/html\n\n'
        print '''
        <html>
        <head><title>UCD GoodVibes</title></head>
        <body>
        '''
        print htmlbody
        print '''
        </body>
        </html>
        '''
        
        if error:
            sys.stdout.flush()
            os.close(0)
            os.close(1)
            os.close(2)

    def slowhtml(self, htmlbody, job):
        
        html = '''
        <html>
        <head><title>GoodVibes</title></head>
        <body>
        '''
        html += htmlbody
        html += '''
        </body>
        </html>
        '''
        fd = open(self.path_results+job+'.html', 'w')
        fd.write(html)
        fd.close()

if __name__=='__main__':
    instance_server = server()
    instance_server.main()

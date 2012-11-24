#!/bin/env /software/bin/python2.3
#
# $Id: server.py 111 2007-04-30 11:31:12Z tc $
#
# changes import goodvibes_dev to import goodvibes in nondev ver
## give option to choose number of frames and average amplitude

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
        self.path_python = '/var/www/cgi-bin/goodvibes/'


        ## initiate fasthtml
        print 'Content-type: text/html\n\n'
        print '''
        <html>
        <head><title>UCD GoodVibes</title></head>
        <body>
        '''


        ## check form for erros and terminate fasthtml if any errors
        self.errorcheckfile(form)


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


        ## if no errors and job queued, then write fast output and slow 1st temporary output
        sys.stdout.flush()
        sys.stderr.flush()
        pid = os.fork()
        if pid:
            htmlbody = 'Dear user. Your input data has been queued for processing. Your job is still in the queue.'
            self.slowhtml(htmlbody, timestr, self.path_results)

            htmlbody = 'Dear user. Your input data is being processed. Your results or an estimated time for the completion of your calculation will be available from <a href="http://peat.ucd.ie/goodvibes/results/'+timestr+'.html" target="_blank">http://peat.ucd.ie/goodvibes/results/'+timestr+'.html</a>'
            self.fasthtml(htmlbody)

##          flush the output from fasthtml
            sys.stdout.flush()
##            os._exit(0)
##
##        sys.stdout.flush() # double flush used in pKD
##        os.close(0)
##        os.close(1)
##        os.close(2)
##        fd=open('junkout','w')
##        sys.stderr=fd
##        sys.stdout=fd

################################################################################

        ## parse html form

        if form.has_key("amplitude_average"): amplitude_average = form["amplitude"].value
        else: amplitude_average = 2
        if form.has_key("mutations"): mutations = form["mutations"].value
        else: mutations = 0
        if form.has_key("quarternary"): quarternary = form["quarternary"].value
        else: quarternary = 'monomeric'
        if form.has_key("frames"): frames = form['frames'].value
        else: frames = 50
        if form["pdb1_id"].value != "":
            fd = open(self.path_pdbs+form["pdb1_id"].value.lower()+'.pdb', 'r')
            pdb_lines = fd.readlines()
            fd.close()
            fd = open(self.path_tmp+timestr+'_reference.pdb', 'w')
            fd.writelines(pdb_lines)
            fd.close()
        if form["pdb1_file"].value:
            fd = open(self.path_tmp+timestr+'_reference.pdb', 'w')
            fd.writelines(form["pdb1_file"].file.readlines())
            fd.close()
        if form["pdb1_txt"].value:
            fd = open(self.path_tmp+timestr+'_reference.pdb', 'w')
            fd.writelines(form["pdb1_txt"].file.readlines())
            fd.close()
        if form["chains1"].value: chains1 = form["chains1"].value.upper()
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
        lines_input = []
        lines_input.append('%s;%s\n' %('chains', chains1))
        lines_input.append('%s;%s\n' %('model', model1))
        lines_input.append('%s;%s\n' %('resrangestart', resrange1start))
        lines_input.append('%s;%s\n' %('resrangeend', resrange1end))
        lines_input.append('%s;%s\n' %('amplitude', amplitude_average))
        lines_input.append('%s;%s\n' %('atoms_hessian', ','.join(atoms_hessian)))
        lines_input.append('%s;%s\n' %('mutations', mutations))
        lines_input.append('%s;%s\n' %('quarternary', quarternary))
        lines_input.append('%s;%s\n' %('frames', frames))

        fd = open(self.path_tmp+'job_'+timestr+'.txt', 'w')
        fd.writelines(lines_input)
        fd.close()
        
        ## send mail via SMTP server to user about job being queued
        import smtplib
        from email.MIMEText import MIMEText
        mailfrom = 'GoodVibes'
        mailto = form['mail'].value
        msg = MIMEText('Dear user. Your input data is being processed. Your results or an estimated time for the completion of your calculation will be available from http://peat.ucd.ie/goodvibes/results/'+timestr+'.html <http://peat.ucd.ie/goodvibes/results/'+timestr+'.html>.')
        msg['Subject'] = 'Your GoodVibes job with the ID %s has been queued for processing' %timestr
        msg['From'] = mailfrom
        msg['Reply-To'] = 'tommy.carstensen@ucd.ie'
        msg['To'] = mailto

        server=smtplib.SMTP('mail.ucd.ie')
        server.set_debuglevel(1)
        server.sendmail(mailfrom, mailto, msg.as_string())
        server.close()

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
        import goodvibes3
        instance_NMA = goodvibes3.vibration()

        ## process queue while it exists
        while os.path.isfile(self.path_tmp+'queue.txt'):

            ## parse job id
            fd = open(self.path_tmp+'queue.txt', 'r')
            job = fd.readline().strip()
            fd.close()

            ## write slow 2nd temporary output
            htmlbody = 'Dear user. Your input is being processed. Return to this page in a few seconds to get an estimate of the end time of the calculations.'
            self.slowhtml(htmlbody, job, self.path_results)

            ## parse parameters from job txt file to goodvibes
            fd = open(self.path_tmp+'job_'+job+'.txt', 'r')
            lines_job = fd.readlines()
            fd.close()

            input_job = {}
            for line in lines_job:
                input_job[line.split(';')[0]] = line.split(';')[1][:-1]
            if input_job['chains'] != '': chains1 = input_job['chains'].split(',')
            else: chains1 = []
            model1 = input_job['model'].strip()
            resrange1start = int(input_job['resrangestart'].strip())
            resrange1end = int(input_job['resrangeend'].strip())
            amplitude_average = float(input_job['amplitude'].strip())
            atoms_hessian = input_job['atoms_hessian'].split(',')
            quarternary = input_job['quarternary'].strip()
            frames = float(input_job['frames'])

            fd = open(self.path_tmp+job+'_reference.pdb', 'r')
            pdblines1 = fd.readlines()
            fd.close()


            ## run goodvibes
##            try:
            self.queue()
            N = instance_NMA.main(
                job,pdblines1,
                atoms_hessian = atoms_hessian,
                frames = frames,
                chains = chains1,
                path_html = self.path_results, path_python = self.path_python,
                paralleldir = self.path_tmp,
                )
##            except 'chainerror', error:
##                ## if expected error then remove from queue and report to user
##                self.slowhtml(error, job, self.path_results)
##                self.queue()
##                continue
##            except IOError, (errno,strerror):
##                self.slowhtml(str(errno)+str(strerror), job, self.path_results)
##                self.queue()
##                results = instance_NMA.main(pdblines1,chains1,model1, atoms_hessian, quarternary, job, frames, [1], path_html = self.path_results, path_python = self.path_python)
##                continue
##            except:
##                print chains1,model1, atoms_hessian, quarternary, job, cutoffs, 
##                ## if not expected error (not syntax errors) then remove from queue and report to user
##                self.slowhtml('There was a bug. Please email your input to <a href="mailto:tommy.carstensen@ucd.ie?subject=GoodVibes Server&body=Hi Tommy,">the administrator</a>.', job, self.path_results)
##                self.queue()
##                results = instance_NMA.main(pdblines1,chains1,model1, atoms_hessian, quarternary, job, cutoffs, [1], path_html = self.path_results, path_python = self.path_python)
##                continue


            ###################################
            ## write output if no exceptions ##
            ###################################

            ## define dictionary of img files
            img_files = ['combo_displacement', 'displacement',
                         'combo_crosscorrelation', 'crosscorrelation',
                         'overlaps_combined', 'overlaps_single',
                         'mmo', 'emo']
#                         {'2D_atom_displacement': ['combo', 'atom displacement per coordinate'], 'crosscorrelation': ['combo', 'cross correlation'], 'emo': ['single', 'perturbed eigenvalues of max overlap'], 'mmo': ['single', 'perturbed modes of max overlap'], 'overlaps': ['combo', 'overlaps'], 'eigenvalues' : ['single', 'perturbed eigenvalues']}

            ## define html header
            header = [
                '<html>\n<head>\n<title>GoodVibes results for job %s</title>\n</head>\n<body>\n' %(job),
                '<a href="http://peat.ucd.ie/goodvibes"><img src="http://peat.ucd.ie/goodvibes/logo.png" alt="GoodVibes" border="0"></a>\n<br>\n<hr>\n<br>\n'
                ]

            ##
            ## initiate main html
            ##
            html = header+['<table border="1">\n']

            ##
            ## write html output
            ##
            html += [
                '<tr><td>Job</td><td>%s</td></tr>\n' %(job),
                '<tr><td>Chains</td><td>%s</td></tr>\n' %(chains1),
                '<tr><td>Number (3N-6) of Non-Zero Eigenvalues</td><td>%s</td></tr>\n' %(3*N-6),
                '<tr><td>Distance Cutoff (Angstrom)</td><td>%s</td></tr>\n' %(10),
                '</table>\n<table border="1">\n'
                ]
                     
            ## 
            ## write html, pdb and gif output
            ## write results to subhtml per mode
            ##
            html += [
                '<tr>\n'
                '<td>mode</td>\n',
                '<td>pdb trajectory</td>\n',
                '<td>calculation results</td>\n',
                '</tr>\n'
                ]
            for mode in range(6,12):

                fd = open('%s_mode%s.pdb' %(job, str(mode+1).zfill(2)), 'r')
                lines_vmd = fd.readlines()
                fd.close()

                ## loop over frames and read/write rasmol input/output
                for frame in range(frames):

                    ## write rasmol pdb for frame (speed up by avoiding reading the lines from the beginning every time)
                    rasmol_pdb = job+'_'+str(frame)+'.pdb'
                    line_first = 0
                    for line_current in range(len(lines_vmd)):
                        if lines_vmd[line_current] == 'HEADER    frame t=%2i.000\n' %frame:
                            line_first = line_current+2
                        if lines_vmd[line_current] == 'TER\n' and line_first != 0:
                            fd = open(rasmol_pdb, 'w')
                            fd.writelines(lines_vmd[line_first:line_current])
                            fd.close()
                            break

                    ## write rasmol script
                    rasmol_script = self.path_tmp+'rasmol.script'
                    fd = open(rasmol_script,'w')
                    fd.write('/software/bin/rasmol -nodisplay '+rasmol_pdb+' << EOF\n')
                    fd.write('color temperature\n') ## color by displacement
                    ## execute rasmol script and write rasmol img
                    ## somehow center at center of wildtype... by translation only! no rotation!
                    rasmol_img = job+'_'+str(frame)+'.ppm'
                    fd.write('write '+self.path_tmp+rasmol_img+'\n')
                    fd.write('exit\n')
                    fd.close()
                    ## execute rasmol script to open pdb frame file and do screen shot
                    os.system('source '+rasmol_script+' > '+self.path_tmp+'rasmol.log') ## write image
                    os.system('/usr/bin/convert '+self.path_tmp+rasmol_img+' -resize 15% '+self.path_tmp+job+'_'+str(frame)+'.gif') ## resize image

                ## convert screen shots to gif animation
                ImageMagick_script = '/usr/bin/convert -loop 0 -delay 2x'+str(frames) ## -delay fps (sxf)
                for frame in range(frames):
                    rasmol_gif = job+'_'+str(frame)+'.gif'
                    ImageMagick_script += ' '+self.path_tmp+rasmol_gif
                for frame in range(frames-1,-1+1,-1):
                    rasmol_gif = job+'_'+str(frame)+'.gif'
                    ImageMagick_script += ' '+self.path_tmp+rasmol_gif
                ImageMagick_script += ' %s%s_mode%s.gif' %(self.path_results,job,str(mode+1).zfill(2))
                os.system(ImageMagick_script)

                ## write results to table
                html += [
                    '<tr>\n'
                    '<td>mode %s</td>\n' %(mode+1),
                    '<td><a href="%s_mode%s.pdb"><img src="%s_mode%s.gif" border="0" alt="pdb trajectory mode%s"></a></td>\n' %(job, str(mode+1).zfill(2), job, str(mode+1), str(mode+1).zfill(2)),
                    '<td><a href="%s_mode%s.html">results</a></td>\n' %(job, mode+1),
                    '</tr>\n',
                    ]

                ## write results to subhtml per mode
                subhtml = header+['<table border="1">\n']
                for img_file in img_files:
                    subhtml += [
                        '<tr><td><img src="%s_%s_mode%s.png" border="0" alt="%s"></td></tr>\n' %(job, img_file, mode+1, img_file),
                        ]
                subhtml += ['</table>\n</body>\n</html>']
                fd = open('%s_mode%s.html' %(job, mode+1), 'w')
                fd.writelines(subhtml)
                fd.close()

            ## terminate mode result table and initiate plot result table
            html += ['</table>\n<table border="1">\n']

            ##
            ## write img output to main html
            ## write results to subhtml per img file type
            ##
            for img_file in img_files:
                ## write subhtml per img file type
                subhtml = header+['<table border="1">\n']
                for mode in range(6,12):
                    subhtml += ['<tr>\n']
                    subhtml += ['<td><img src="%s_%s_mode%s.png" border="0" alt="mode%s"></td>\n' %(job, img_file, str(mode+1).zfill(2), str(mode+1).zfill(2))]
                    subhtml += ['</tr>\n']
                if img_file == 'overlapscombined':
                    subhtml += [
                        '<tr><td><img src="%s_overlapscombined_mode3N" border="0" alt="mode3N"></td></tr>\n' %(job),
                        ]
                subhtml += ['</table>\n</body>\n</html>']
                fd = open('%s_%s.html' %(job, img_file), 'w')
                fd.writelines(subhtml)
                fd.close()
                ## write main  html
                html += ['<tr><td><a href="%s_%s.html">%s</a></td></tr>\n' %(job, img_file, img_file)]

            ##
            ## terminate main html
            ##
            html += ['</table>\n</body>\n</html>']

            ## delete temporary files
            for frame in range(frames):
                os.remove('%s%s_%s.ppm' %(self.path_tmp, job, frame))
                os.remove('%s%s_%s.gif' %(self.path_tmp, job, frame))

            ##
            ## write html to file (this output cannot be written prior to the subpages outputs)
            ##
            fd = open(self.path_results+job+'.html', 'w') ## implicit path
            fd.writelines(html)
            fd.close()

            ###############################################################
            ## remove the job from the queue and continue the while loop ##
            ###############################################################
            self.queue()
            continue

        ## return/end when no more jobs to process in queue
        return


    def errorcheckfile(self, form):

        import os, sys
        sys.path.append(self.path_python)
        import goodvibes3


        ## check for essential input
        if not ( form["pdb1_id"].value or form["pdb1_file"].value or form["pdb1_txt"].value ):
            self.fasthtml('Dear user. You need to specify a pdb to be used.', error = True)
            self.queue()


        ## check for multiple input
        if ( form["pdb1_id"].value and form["pdb1_file"].value ) or ( form["pdb1_id"].value and form["pdb1_txt"].value ) or ( form["pdb1_file"].value and form["pdb1_txt"].value ):
            self.fasthtml('Dear user. You must use *either* a PDB ID, an attached PDB file *or* a pasted PDB file as input.', error = True)
            self.queue()


        ## check for valid pdb ids if specified
        if form["pdb1_id"].value:
            if not os.path.isfile(self.path_pdbs+form["pdb1_id"].value.lower()+'.pdb'):
                self.fasthtml('The PDB ID you have specified is not found in the protein data bank.', error = True)
                self.queue()


        ## check that a mail address has been given
        if not form["mail"].value:
            self.fasthtml('Please specify an address where you want status about your job sent to.', error = True)
            self.queue()


        ## check for valid or missing chain id(s)

        ## parse chains
        if form["pdb1_id"].value:
            fd = open(self.path_pdbs+form["pdb1_id"].value.lower()+'.pdb', 'r')
            pdb_lines = fd.readlines()
            fd.close()
        if form["pdb1_file"].value != "":
            pdb_lines = form["pdb1_file"].file.readlines()
        if form["pdb1_txt"].value != "":
            pdb_lines = form["pdb1_txt"].file.readlines()
        instance_NMA = goodvibes3.vibration()
        d_coordinates = instance_NMA.parse_pdb(pdb_lines)[3]

        ## check for valid chain ids
        if form["chains1"].value:
            if not form["chains1"].value.upper() in d_coordinates['chains'].keys():
                self.fasthtml('Specified chain ID not present in PDB. Only chains (%s) are present in pdb.' %d_coordinates['chains'].keys(), error = True)
                self.queue()

        ## check that only 1 chain is present if none specified by user
        if form["chains1"].value == "":
            ## error if no chains present
            if len(d_coordinates['chains'].keys()) == 0:
                self.fasthtml('No protein chains present in pdb. GoodVibes does not compute polynucleotides at present.', error = True)
                self.queue()
            ## error if multiple chains
            if len(d_coordinates['chains'].keys()) > 1:
                self.fasthtml('Multiple chains (%s) present in pdb. Please specify which one to use.' %d_coordinates['chains'].keys(), error = True)
                self.queue()


        ## check if input exceeds a certain allwed length


        ## check for valid number of frames
        if form.has_key("frames"):
            if int(form["frames"].value) > 500 or int(form["frames"].value) < 2:
                self.fasthtml('Please change the number of frames to be within the range 10 to 500.', error = True)
                self.queue()

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

    def slowhtml(self, htmlbody, job, path_results):
        
        html = [
            '<html>\n<head>\n<title>GoodVibes</title>\n</head>\n<body>\n',
            '<img src="http://peat.ucd.ie/goodvibes/logo.png" alt="GoodVibes" border="0">\n'
            ]
        html += [htmlbody]
        html += ['</body>\n</html>']
        fd = open(path_results+job+'.html', 'w')
        fd.writelines(html)
        fd.close()

if __name__=='__main__':
    instance_server = server()
    instance_server.main()

#!/bin/env /software/bin/python2.3
# $Id: download.py 111 2007-04-30 11:31:12Z tc $
# Copyright (C) Tommy Carstensen, University College Dublin, 2007

def main():

    ## enable traceback, which reports non-syntax errors
    import cgitb
    cgitb.enable()

    ## read form
    import cgi
    form = cgi.FieldStorage()

    htmlheader = '''
<html>
<head>
 <title>GoodVibes Registration</title>
 <link rel="icon" href="favicon.ico">
 <link rel="shortcut icon" href="favicon.ico">
</head>

<body>

 <a href="http://polymerase.ucd.ie/goodvibes/"><img src="http://polymerase.ucd.ie/goodvibes/logo.png" border="0"></a>

 <hr>
'''

    htmltail = '\n</body>\n</html>'

    path_download = '/var/www/html/goodvibes/download/'

    if form['answer'].value == 'I agree to the terms and conditions for using the GoodVibes software':

        import tempfile
        tmpfile = tempfile.mkstemp(prefix='goodvibes_', suffix='.py', dir=path_download)

        ## copy python script
        import os
        os.system('cp /var/www/cgi-bin/goodvibes/python/goodvibes/goodvibes.py %s' %(tmpfile[1]))
        ## parse tmpfilename
        tmpfilename = tmpfile[1][-19:]
        ## link to tmpfile in html
        htmlbody = 'Download your copy of GoodVibes by right-clicking the link below and choosing "Save as..."\n<br><br>\n'
        htmlbody += '<a href="http://polymerase.ucd.ie/goodvibes/download/%s">DOWNLOAD by clicking here</a>' %(tmpfilename)

        fd = open('/var/www/cgi-bin/goodvibes/download/users.txt','a')
        fd.write(
            '\n%s\n%s\n%s\n%s\n'
            %(
                form['name'].value,
                form['organisation'].value,
                form['country'].value,
                form['email'].value,
                )
            )
        fd.close()
        
    elif form['answer'].value == 'I do not agree':
        htmlbody = 'You cannot download GoodVibes without accepting the agreement. Thanks for showing interest in GoodVibes.'

    print htmlheader+htmlbody+htmltail


if __name__ == '__main__':
    main()

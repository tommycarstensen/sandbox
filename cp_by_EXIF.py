#!/bin/python

## Tommy Carstensen, November 2012

import os
import pyexiv2
import shutil

def main():

    path = '/media/WDPassPort1TB/WDMyBook1TB/Photos'
##    path = '/media/WDMyBook1TB/WDPassPort1TBRed/Photos'
    l_dirs = []
    l_dirs = loop(path,l_dirs)

    path_out = '/home/tommy/Dropbox/billeder/'
##    path_out = '/home/tommy/Dropbox/billeder_Susanne/'
    copy_dirs(path,path_out,l_dirs)

    return


def copy_dirs(path_in,path_out,l_dirs,):

    for path_dn in l_dirs:
        path = path_dn[len(path_in)+1:]
        print path
        l = path.split('/')
        for i in xrange(len(l)):
            path = os.path.join(path_out,'/'.join(l[:i+1]),)
            if not os.path.isdir(path):
                os.mkdir(path)

        path = path_dn[len(path_in)+1:]
        src = os.path.join(path_in,path)
        dst = os.path.join(path_out,path)
        l = os.listdir(src)
        for fn in l:
            fp_out = os.path.join(dst,fn)
            if os.path.isdir(fp_out):
                continue
            if os.path.isfile(fp_out):
                continue
            fp_in = os.path.join(src,fn)
            if os.path.isdir(fp_in):
                continue
            shutil.copy(fp_in,fp_out)

    return


def loop(dn,l_dirs,):

    l = os.listdir(dn)
    bool_found = False
    for s in l:
        fp = os.path.join(dn,s)
        if os.path.isdir(fp):
            print fp
            l_dirs = loop(fp,l_dirs,)
            pass
        elif os.path.isfile(fp):
            extension = s[s.rindex('.')+1:].lower()
            if extension in [
                'zip','db','pdf', ## 
                'info','ini','txt','py','html', ## text
                'db-journal', ## digikam
##                'wmv','3gp','mp4', ## moved to Videos directory
                'wlmp', ## video
                ]:
                continue
            metadata = pyexiv2.ImageMetadata(fp)
            metadata.read()
##            l_keys += [
####                'Xmp.MicrosoftPhoto.LastKeywordXMP', ## list, e.g. People/Tommy
####                'Xmp.digiKam.TagsList', ## list, e.g. People/Tommy
####                'Xmp.lr.hierarchicalSubject', ## list, e.g. People|Tommy
##                ]
            for key in metadata.xmp_keys:
                if not 'Xmp.MP.RegionInfo' in key:
                    continue
                if not 'PersonDisplayName' in key:
                    continue
                tag = metadata[key].raw_value
                if tag in [
                    
##                    ## Susanne
####                    'People/Vinter' in tag
##                    'Susanne Vinter' in tag ## e.g. Xmp.MP.RegionInfo/MPRI:Regions[3]/MPReg:PersonDisplayName
##                    or
##                    'Jan Vinter' in tag
##                    or
##                    'Casper Vinter' in tag

                    ## Elisa
                    'Carstensen Elisa Carmen',
                    'Carstensen Tom',
                    'Salamanca Camacho Antonia',
                    'Carstensen Henry Hermann',
                    '(Carstensen) Olsen Grete Karen Sofie',
                    'Camacho Fernandez, Carmen',
                    'Joergensen Helga Marie',
                    'Olsen Olaf Henrik',
                    'Olsen Emma',

##                    ## Finn
##                    '(Carstensen) Olsen Grete Karen Sofie',
##                    'Olsen Emma',
##                    'Olsen Olaf Henrik',
##                    'Joergensen Helga Marie',
##                    'Carstensen Henry Hermann',
##                    Erik
##                    Elisabeth
##                    Finn

##                    ## Benny og Margrethe
##                    'Benny',
##                    'Margrethe',
##                    'Dorthe'
##                    'Pernille',
                    ]:
                    print fp, tag
                    bool_found = True
                    break
                elif 'Antonia' in str(tag):
                    print key,tag
                    print fp
                    stop
                elif 'Elisa' in str(tag) and 'Elisabeth' not in str(tag) and 'Noguiera' not in str(tag):
                    print key,tag
                    stop
                elif 'Henry' in str(tag):
                    print key,tag
                    stop
                elif 'Tom' in str(tag) and 'Tommy' not in str(tag):
                    print key,tag
                    stop
                elif 'Camacho' in str(tag) and 'Carmen' in tag and tag != 'Carmen Salamanca Camacho':
                    print key,tag
                    stop
                elif 'Grete' in str(tag):
                    print key,tag
                    stop
                elif 'Helga' in str(tag):
                    print key,tag
                    stop
                elif 'Olaf' in str(tag):
                    print key,tag
                    stop
##                elif 'Olsen' in str(tag) and 'Huse' not in str(tag) and 'Lisbeth' not in str(tag) and 'Glen' not in str(tag):
##                    print key,tag
##                    stop
                elif 'Salamanca' in str(tag) and 'Antonia' in str(tag):
                    print fp,key,tag
                    stop
##                if bool_found == True:
##                    break ## tmp exclusion
                continue ## continue loop over xmp_keys
            if 'Xmp.MicrosoftPhoto.LastKeywordXMP' in metadata: ## tmp
                if 'Antonia' in str(metadata['Xmp.MicrosoftPhoto.LastKeywordXMP'].raw_value):
                    if bool_found == False: ## tmp
                        print fp ## tmp
                        stop ## tmp
            pass
        else:
            print fp
            print extension
            stop
            pass

##        if bool_found == True:
##            break

        continue

    if bool_found == True:
        l_dirs += [dn]

    return l_dirs

if __name__ == '__main__':
    main()

##            ## http://zoia.org/files/python/xmppic_py
##            import libxml2
##            img = Image.open(fp)
##            rdf = img.app['APP1']
##            print rdf
##            if 'Carstensen' in rdf:
##                rdf.replace("\x00", "")  
##                rdf0 = rdf[rdf.find('<?'):]
##                rdf = rdf0.strip()
##                xml = libxml2.parseDoc(rdf)
##                print xml

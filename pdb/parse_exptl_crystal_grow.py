import os, parse_mmCIF, numpy, math, sys

def main():

    l_fn_out = [
        '_exptl_crystal_grow',
        '_exptl_crystal_grow_comp',
        ]

    d = {}
    for fn_out in l_fn_out:
        fd = open('db%s.txt' %(fn_out),'r')
        lines = fd.readlines()
        fd.close()
        d[fn_out] = {}
        for line in lines:
            if line == '\n':
                continue
            pdb = line[:4]
            s = line[5:]
            d[fn_out][pdb] = s

    fd = open('remediation_exptl_crystal_grow.pH.txt','r')
    lines = fd.readlines()
    fd.close()
    l_pdbs = [line[:4] for line in lines]

    path = '/media/WDMyBook1TB/2TB/mmCIF'
    l_dns = os.listdir(path)
    l_dns.sort()
    for i in range(len(l_dns)):
        dn = l_dns[i]

        if dn < sys.argv[-1]:
            continue

        if not os.path.isdir('%s/%s' %(path,dn)):
            continue

        print '%s/%s %s' %(i+1,len(l_dns), dn)
        l_fns = os.listdir('%s/%s' %(path,dn))
        l_fns.sort()
        for fn in l_fns:
            if fn[-3:] == '.gz':
                continue

            pdb = fn[0:4]

            ## continue if already in txt file from previous attempt to run loop
##            if pdb in d['_exptl_crystal_grow_comp'].keys():
##                continue
##            if pdb in d['_exptl_crystal_grow'].keys():
##                continue

##            print pdb

            if not pdb in l_pdbs:
                continue

            fd = open('%s/%s/%s' %(path, dn, fn), 'r')
            lines = fd.readlines()
            fd.close()

            d_mmCIF = parse_mmCIF.main(
                pdb,lines,
                d_breaks_negation = {
                    ## break if not x-ray diffraction
                    '_exptl.method':'X-RAY DIFFRACTION',
                    },
                l_data_categories_break = [
##                    '_atom',
                    '_diffrn',
                    ],
                l_data_categories = [
                    ## parse selected data categories
                    '_database_PDB_rev',
                    '_pdbx_database_status',
                    '_exptl',
                    '_exptl_crystal_grow',
                    '_exptl_crystal_grow_comp',
                    ],
                )

##            ## no polymers in structure?
##            if not '_entity_poly.entity_id' in d_mmCIF.keys():
##                continue

            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                continue

##            print pdb

            ##
            ##
            ##
            year = int(d_mmCIF['_database_PDB_rev.date'][0][:4])
            process_site = ''.join(d_mmCIF['_pdbx_database_status.process_site'])
            if (
                not '_exptl_crystal_grow.pdbx_details' in d_mmCIF.keys()
                and
                not '_exptl_crystal_grow_comp.name' in d_mmCIF.keys()
##                ''.join(d_mmCIF['_exptl_crystal_grow.pdbx_details']).strip() == '?'
                ):
                if process_site != '?':
                    print pdb, year, process_site
                continue

            ##
            if '_exptl_crystal_grow.pdbx_details' in d_mmCIF.keys():

                s_grow = ' '.join(d_mmCIF['_exptl_crystal_grow.pdbx_details']).strip()

                if (
                    ## pH not given
                    d_mmCIF['_exptl_crystal_grow.pH'] in [['?'],[''],['.'],]
                    and
                    d_mmCIF['_exptl_crystal_grow.pdbx_pH_range'] in [['?'],[''],['.'],]
                    and
                    ## but pH in growth details
                    (
                        ' PH ' in s_grow.upper()
                        or
                        'PH=' in s_grow.upper()
                        or
                        ',PH ' in s_grow.upper()
                        )
                    ):
                    fd = open('remediation_exptl_crystal_grow.pH.txt','a')
                    fd.write('%s\t%s\t%s\t%4i\t%s\t%s\n' %(
                        pdb,
                        ''.join(d_mmCIF['_exptl_crystal_grow.pH']),
                        ''.join(d_mmCIF['_exptl_crystal_grow.pdbx_pH_range']),
                        year,
                        process_site,
                        s_grow,
                        )
                             )
                    fd.close()

                if (
                    not '_exptl_crystal_grow_comp.name' in d_mmCIF.keys()
                    or
                    ''.join(d_mmCIF['_exptl_crystal_grow_comp.name']) in ['.','','?',]
                    ):

                    if '_exptl_crystal_grow_comp.name' in d_mmCIF.keys():
                        name = ''.join(d_mmCIF['_exptl_crystal_grow_comp.name'])
                    else:
                        name = 'N/A'

##                    ## remove end punctuation
##                    s = s_grow[:-1]+s_grow[-1].replace('.','')

                    ## split
##                    l_grow_punctuation = s_grow.upper().split('. ')
##                    l_grow = l_grow_comma = [s_grow.upper().split(',') for s in l_grow_punctuation]
                    l_grow = s_grow.upper().split(',')

                    ## strip space
                    l_grow = [x.strip() for x in l_grow]
                    ## remove empty
                    if '' in l_grow:
                        l_grow.remove('')
                    ## remove end punctuation
                    l_grow = [x[:-1]+x[-1].replace('.','') for x in l_grow]

                    ## remove selected words from elements of list
                    for x in [
                        'CRYSTALS OBTAINED BY CO-CRYSTALLIZATION AT ',
                        'PROTEIN SOLUTION (',
                        ]:
                        for i_grow in range(len(l_grow)):
                            l_grow[i_grow] = l_grow[i_grow].replace(x,'')

                    ## replace abbreviations
                    for i_grow in range(len(l_grow)):
                        l_grow[i_grow] = l_grow[i_grow].replace('MILLIMOLAR','MM')
                        
                    
                    ## remove selected words from list
                    l_remove = []
                    for x in [
                        'VAPOR DIFFUSION',
                        'VAPOUR DIFFUSION',
                        'HANGING DROP',
                        'SITTING DROP',
                        ]:
                        if x in l_grow:
                            l_remove += [x]
                            
                    ## removed other selected words from list
                    for i_grow in range(len(list(l_grow))):

                        ## remove physical conditions
                        bool_continue = False
                        for x in [
                            'TEMPERATURE',
                            'PH=',
                            'PH ',
                            'AT PH ',
                            ]:
                            if l_grow[i_grow][:len(x)] == x:
                                l_remove += [l_grow[i_grow]]
                                bool_continue = True
                                break
                        if bool_continue == True:
                            continue

                        ## remove long words (sentences)
                        if len(l_grow[i_grow]) > 50:
                            l_remove += [l_grow[i_grow]]
                            break
                    for remove in l_remove:
                        l_grow.remove(remove)
                    if len(l_grow) > 0:
                        ## write to file
                        line = '%s\t%s\t%s\t%4i\t%s\t%s\n' %(
                            pdb,
                            name,
                            l_grow,
                            year,
                            process_site,
                            s_grow,
                            )
                        fd = open('remediation_exptl_crystal_grow_comp.name.txt','a')
                        fd.write(line)
                        fd.close()

            else:
                s_grow = ''

            ##
            if '_exptl_crystal_grow_comp.name' in d_mmCIF.keys():
                l_grow_comp = d_mmCIF['_exptl_crystal_grow_comp.name']
            else:
                l_grow_comp = []

##            lines_out += [line]

            ## append to txt file in case loop doesn't finish
            d_lines = {}
            line = '%s %s\n' %(pdb,s_grow,)
            d_lines['_exptl_crystal_grow'] = line
            line = '%s %s\n' %(pdb,l_grow_comp,)
            d_lines['_exptl_crystal_grow_comp'] = line
            for fn_out in l_fn_out:
                fd = open('db%s.txt' %(fn_out),'a')
                fd.write(d_lines[fn_out])
                fd.close()

            ## append to dic for when loop finishes
            d['_exptl_crystal_grow'][pdb] = s_grow
            d['_exptl_crystal_grow_comp'][pdb] = l_grow_comp

    lines_out = []
    for pdb,s in d.items():
        line = '%s %s\n' %(pdb,s,)
        lines_out += [line]
    fd = open(fn_out,'w')
    fd.writelines(lines_out)
    fd.close()

    return

if __name__ == '__main__':
    main()

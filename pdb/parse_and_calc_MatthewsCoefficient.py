import os, parse_mmCIF, numpy, math, sys
import calc_matthews_coefficient

def main():

    fn_out = 'db_MatthewsCoefficient.txt'

    fd = open(fn_out,'r')
    lines = fd.readlines()
    fd.close()

    d = {}
    for line in lines:
        l = line.strip().split()
        pdb = l[0]
        v = l[1]
        if pdb == '2p51':
            v = '1.72610466393'
        d[pdb] = v

    lines_out = []

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

            if pdb in d.keys():
                continue

            ## Matthews Coefficient not calculated...
            if pdb in [
                '1vh7','1vho','1vhu','1vi3','1vi4','1vis',
                ]:
                continue

            ## Matthews Coefficient *wrong*
            if pdb in [
                '2p51',
                ## too high
                '1c5v','1q9i','1ut6','1x6x','1x6y','1xdn','1y63','1zix',
                ## too low
                '1t95','1jih','1t95','1d5t','1c7k',
                '1dbo','1d9x','1qt9','1ia5','1dcq',
                ]:
                continue

            fd = open('%s/%s/%s' %(path, dn, fn), 'r')
            lines = fd.readlines()
            fd.close()

            d_mmCIF = parse_mmCIF.main(
                pdb,lines,
                d_breaks = {
                    ## break if multiple polymer types (not monomeric)
                    '_entity_poly.entity_id':'2',
##                    '_exptl.method':'SOLUTION NMR', ## break if e.g. _exptl.method = SOLUTION NMR
                    ## break if multiple chains
                    '_entity_poly.pdbx_strand_id':',',
                    }, 
                d_breaks_negation = {
                    ## break if not x-ray diffraction
                    '_exptl.method':'X-RAY DIFFRACTION',
                    ## break if not monomeric
                    '_pdbx_struct_assembly.oligomeric_details':'monomeric',
                    },
                l_data_categories = [
                    '_exptl_crystal',
                    ], ## parse selected data categories
                l_data_categories_break = ['_exptl_crystal']
                )

            ## some unknown temporary error... or break before reaching this part when parsing...
            if not '_pdbx_struct_assembly.oligomeric_details' in d_mmCIF.keys():
                continue

            ## NMR structure?
            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                stop2
                continue

            ## no polymers in structure?
            if not '_entity_poly.entity_id' in d_mmCIF.keys():
                continue

            ## polymer(s) is/are not polypeptide(s)
            if d_mmCIF['_entity_poly.type'] != len(d_mmCIF['_entity_poly.type'])*['polypeptide(L)']:
                continue

            ## biounit not monomeric
            if d_mmCIF['_pdbx_struct_assembly.oligomeric_details'] != len(d_mmCIF['_pdbx_struct_assembly.oligomeric_details'])*['monomeric']:
                continue

            ## one polymer in assymetric unit
            if len(d_mmCIF['_entity_poly.entity_id']) > 1:
                continue

            if d_mmCIF['_exptl_crystal.density_Matthews'] == ['?']:
                v = VM = calc_matthews_coefficient.main(pdb)
##                continue
            else:
                v = float(''.join(d_mmCIF['_exptl_crystal.density_Matthews']))

            line = '%s %s\n' %(pdb,v,)

            fd = open(fn_out,'a')
            fd.write(line)
            fd.close()

            d[pdb] = v

    ##
    ## write calculated radii of gyration to file
    ##
    lines_out = []
    for pdb,v in d.items():
        line = '%s %s\n' %(pdb,v,)
        lines_out += [line]
    fd = open(fn_out,'w')
    fd.writelines(lines_out)
    fd.close()

    return

if __name__ == '__main__':
    main()

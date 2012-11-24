import os, parse_mmCIF, numpy, math, sys

def main():

    fd = open('remediation_negativeBiso.txt','r')
    lines = fd.readlines()
    fd.close()
    l_pdbs = []
    for line in lines:
        if line.strip() == '':
            continue
        if line[0] == '#':
            continue
        l = line.strip().split()
        pdb = l[0]
        l_pdbs += [pdb]

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

            if not pdb in l_pdbs:
                continue

            print pdb

            fd = open('%s/%s/%s' %(path, dn, fn), 'r')
            lines = fd.readlines()
            fd.close()

            d_mmCIF = parse_mmCIF.main(
                pdb,lines,
                d_breaks_negation = {
                    ## break if not x-ray diffraction
                    '_exptl.method':'X-RAY DIFFRACTION',
                    },
                l_data_categories = [
                    ## parse selected data categories
                    '_database_PDB_rev',
                    '_computing',
                    '_atom_site',
                    '_refine'
                    ],
                )

##            ## no polymers in structure?
##            if not '_entity_poly.entity_id' in d_mmCIF.keys():
##                continue

            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                continue

            print pdb

            ##
            ## parse bfactors
            ##
            for i_atom_site in range(len(d_mmCIF['_atom_site.id'])):

                bfactor = float(d_mmCIF['_atom_site.B_iso_or_equiv'][i_atom_site])

##                if bfactor == '?':
##                    continue

                element = d_mmCIF['_atom_site.type_symbol'][i_atom_site]
                comp_id = d_mmCIF['_atom_site.label_comp_id'][i_atom_site]

                if float(bfactor) < -0.01:
                    if (
                        element != 'H'
                        and
                        comp_id in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]
                        ):

                        print
                        print 'negative'
                        print

                        year = int(d_mmCIF['_database_PDB_rev.date'][0][:4])
                        atom_id = int(d_mmCIF['_atom_site.id'][i_atom_site])
                        refinement = ''.join(d_mmCIF['_computing.structure_refinement'])
                        solution = ''.join(d_mmCIF['_computing.structure_solution'])
                        resolution = float(''.join(d_mmCIF['_refine.ls_d_res_high']))
                        
                        fd = open('remediation_negativeBiso.txt','a')
                        fd.write(
##                            '%4s %4i %4i %3s %2s %6.2f %30s %20s\n' %(
                            '%4s\t%4i\t%4i\t%3s\t%2s\t%6.2f\t%6.2f\t%30s\t%20s\n' %(
                                pdb,year,atom_id,
                                comp_id,element,bfactor,resolution,
                                solution.ljust(30),refinement.ljust(20),
                                )
                            )
                        fd.close()
                        break

    return

if __name__ == '__main__':
    main()

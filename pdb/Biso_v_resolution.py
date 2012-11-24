import os, sys
import parse_mmCIF
sys.path.append('/home/tc/svn/tc_sandbox/misc/')
import gnuplot
sys.path.append('/home/tc/svn/tc_sandbox/math/')
import statistics

path = '/media/WDMyBook1TB/2TB/mmCIF'

def main():

    l_pdbs = []
    fd = open('Biso_v_resolution.gnuplotdata','r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        l = line.split()
        resolution = float(l[1])
        Biso = float(l[0])
        if resolution > 3.5 and Biso < 10:
            print line
        if resolution > 2.5 and Biso < 10:
            print line
        if resolution > 2.0 and Biso < 5:
            print line
##        if resolution > 1.5 and Biso < 5:
##            print line
        pdb = l[2]
        l_pdbs += [pdb]

    Biso_average_prev = 0

    l_dn = os.listdir(path)
    l_dn.sort()
    for dn in l_dn:
        if dn < sys.argv[-2]:
            continue
        if dn > sys.argv[-1]:
            continue
        if not os.path.isdir('%s/%s' %(path,dn)):
            continue
        l_fns = os.listdir('%s/%s' %(path,dn))
        l_fns.sort()
        for fn in l_fns:
            if fn[-3:] == '.gz':
                continue
            pdb = fn[0:4]

            if pdb in l_pdbs:
                continue

            if pdb in [
                '3bfn', ## PISA left out chains from biological unit
                '2jjg','1qjb', ## _pdbx_struct_assembly missing
                ]:
                continue

            ##
            ## parse header
            ##
            d_mmCIF = parse_mmCIF.main(
                pdb,
                l_data_categories = [
                    '_pdbx_struct_assembly',
                    '_entity_poly',
                    '_citation',
                    '_pdbx_database_related',
                    ], ## parse selected data categories
                d_breaks_negation = {
                    ## break if not x-ray diffraction
                    '_exptl.method':'X-RAY DIFFRACTION',
                    }
                )
            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                continue

            if not 'polypeptide(L)' in d_mmCIF['_entity_poly.type']:
                continue

            if '_pdbx_database_related.content_type' in d_mmCIF.keys():
                if 'split' in d_mmCIF['_pdbx_database_related.content_type']:
                    continue

            try:
                if d_mmCIF['_pdbx_struct_assembly.oligomeric_details'] != ['monomeric']:
                    continue
            except:
                print pdb
                stop

            if not '_citation.id' in d_mmCIF.keys():
                continue

            ##
            ## parse coordinate section
            ##
            d_mmCIF = parse_mmCIF.main(
                pdb,
                l_data_categories = [
                    '_database_PDB_rev',
                    '_refine',
                    '_refine_hist',
                    '_atom_site',
                    '_software',
                    '_entity','_entity_poly',
                    '_pdbx_struct_assembly',
                    '_pdbx_database_status',
                    ], ## parse selected data categories
                d_breaks_negation = {
                    ## break if not x-ray diffraction
                    '_exptl.method':'X-RAY DIFFRACTION',
                    }
                )

            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                continue

            if not 'polypeptide(L)' in d_mmCIF['_entity_poly.type']:
                continue

            if d_mmCIF['_pdbx_struct_assembly.oligomeric_details'] != ['monomeric']:
                continue

            resolution = float(''.join(d_mmCIF['_refine.ls_d_res_high']))

            if (
                int(d_mmCIF['_entity.pdbx_number_of_molecules'][0]) != 1
                or
                len(d_mmCIF['_entity_poly.pdbx_strand_id']) > 1
                or
                len(''.join(d_mmCIF['_entity_poly.pdbx_strand_id']).split(',')) > 1
                or
                len(d_mmCIF['_entity_poly.entity_id']) > 1
                ):
                print pdb
                print d_mmCIF['_entity.pdbx_number_of_molecules']
                print d_mmCIF['_entity_poly.pdbx_strand_id']
                stop

            entity_poly_id = int(''.join(d_mmCIF['_entity_poly.entity_id']))
            for i_entity_poly in range(len(d_mmCIF['_entity_poly.entity_id'])):
                entity_poly_id = d_mmCIF['_entity_poly.entity_id'][i_entity_poly]
                entity_poly_type = d_mmCIF['_entity_poly.entity_id'][i_entity_poly]

            l_Biso = []
            for i_atom_site in range(len(d_mmCIF['_atom_site.id'])):
                occupancy = float(d_mmCIF['_atom_site.occupancy'][i_atom_site])
                if occupancy != 1:
                    continue
                alt_id = d_mmCIF['_atom_site.label_alt_id'][i_atom_site]
                if alt_id != '.':
                    continue
                entity_id = d_mmCIF['_atom_site.label_entity_id'][i_atom_site]
                if entity_id != entity_poly_id:
                    continue
                comp_id = d_mmCIF['_atom_site.label_comp_id'][i_atom_site]
                if not comp_id in ['MSE','ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
                    continue
                type_symbol = d_mmCIF['_atom_site.type_symbol'][i_atom_site]
                if type_symbol == 'H':
                    continue
                atom_id = d_mmCIF['_atom_site.label_atom_id'][i_atom_site]
                if not atom_id in ['N','CA','C',]:
                    continue

                Biso = float(d_mmCIF['_atom_site.B_iso_or_equiv'][i_atom_site])
                l_Biso += [Biso]

            year = int(d_mmCIF['_database_PDB_rev.date'][0][:4])
            site = ''.join(d_mmCIF['_pdbx_database_status.process_site'])

            if len(l_Biso) == 0:
                continue

##            if l_Biso == len(l_Biso)*[l_Biso[0]]:
##                print pdb, year, l_Biso[0:3]
##                if year >= 2010:
##                    stop
##                continue

            Biso_average = sum(l_Biso)/len(l_Biso)

            bool_continue = False
            for Biso in set(l_Biso):
                count = l_Biso.count(Biso)
                if count > 20:
                    if '_software.name' in d_mmCIF.keys():
                        print pdb, Biso_average, Biso, count, d_mmCIF['_software.name']
                        s = '%s %6.2f %4i %6.2f %4i %s %s\n' %(
                            pdb,Biso,count,Biso_average,year,site,
                            str(d_mmCIF['_software.name']),
                            )
                    else:
                        print pdb, Biso_average, Biso, count
                        s = '%s %6.2f %4i %6.2f %4i %s\n' %(
                            pdb,Biso,count,Biso_average, year, site,
                            )
                    bool_continue = True
                    fd = open('remediation_Biso_duplicates.txt','a')
                    fd.write(s)
                    fd.close()
                    break
            if bool_continue == True:
                continue

##            if Biso_average in [2,3,4,5,6,7,8,9,99,90,50,20,25,1,100,10,0]:
            if Biso_average in range(0,100+1):
                print l_Biso
                print Biso_average
                print pdb
                print year
                stop

            if '_refine.pdbx_TLS_residual_ADP_flag' in d_mmCIF.keys():
                if ''.join(d_mmCIF['_refine.pdbx_TLS_residual_ADP_flag']) in ['UNVERIFIED','LIKELY RESIDUAL',]:
                    continue
                elif ''.join(d_mmCIF['_refine.pdbx_TLS_residual_ADP_flag']) in ['?',]:
                    pass
                else:
                    print d_mmCIF['_refine.pdbx_TLS_residual_ADP_flag']
                    print pdb, Biso_average
                    stop

            if round(Biso_average,4) == round(Biso_average_prev,4):
                print pdb, Biso_average, Biso_average_prev
                stop

            print pdb, round(Biso_average,2), resolution
            fd = open('Biso_v_resolution.gnuplotdata','a')
            fd.write('%s %s %s %s\n' %(Biso_average,resolution,pdb,year,))
            fd.close()

    plot()

def plot():

    fd = open('Biso_v_resolution.gnuplotdata','r')
    lines = fd.readlines()
    fd.close()
    d_splot = {}
    l_resolution = []
    l_Biso = []
    for line in lines:
        l = line.split()
        resolution = float(l[1])
        Biso = float(l[0])
        resolution_rounded = round(resolution,1)
##        resolution_rounded = 0.1*round(resolution/0.1,0)
        Biso_rounded = round(Biso,0)
        if not resolution_rounded in d_splot.keys():
            d_splot[resolution_rounded] = {}
        if not Biso_rounded in d_splot[resolution_rounded].keys():
            d_splot[resolution_rounded][Biso_rounded] = 0
        d_splot[resolution_rounded][Biso_rounded] += 1
        l_resolution += [resolution]
        l_Biso += [Biso]
    a,b,r,p = statistics.do_regression(l_resolution,l_Biso,verbose=False,)
    print 'correlation = r =', r
    tcrit = 1.960 ## 95% confidence 2 tail (student table...)
##    tcrit = 76.
##    f,g,h = statistics.do_confidence_bands(l_resolution,l_Biso,tcrit,)

    lines = []
##    for resolution in range(5,35+1):
    for resolution in range(5,35+1):
        resolution /= 10.
        for Biso in range(0,100+1):
            Biso /= 1.
            if not resolution in d_splot.keys():
                count = 0
            elif not Biso in d_splot[resolution].keys():
                count = 0
            else:
                count = d_splot[resolution][Biso]
            lines += ['%s %s %s\n' %(resolution,Biso,count,)]
        lines += ['\n']

    gnuplot.contour_plot(
        'validation_Biso_v_resolution',lines,
        xlabel='resolution (Angstrom)',ylabel='<Biso>',
        bool_remove = False,
        )

    return

if __name__ == '__main__':
    main()


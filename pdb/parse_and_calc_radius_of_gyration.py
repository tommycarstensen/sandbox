import os, parse_mmCIF, numpy, math, sys

d_mass = {
    'H':1.00794,
    'C':12.0107,
    'N':14.0067,
    'O':15.9994,
    'P':30.973762,
    'S':32.065,
    'SE':78.96, ## MSE / SeMet
    }

def main():

    fd = open('radius_of_gyration.txt','r')
    lines = fd.readlines()
    fd.close()
    d_radii = {}
    for line in lines:
        l = line.strip().split()
        pdb = l[0]
        r = l[1]
        d_radii[pdb] = r

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

            if pdb in d_radii.keys():
                continue

            print pdb

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
                    '_atom_site',
                    '_entity_poly',
                    '_pdbx_struct_assembly',
                    ], ## parse selected data categories
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

            print pdb

            ##
            ## calculate center of mass
            ##
            center_of_mass = numpy.array([0.,0.,0.,])
            l_coords = []
            l_masses = []
            for i_atom_site in range(len(d_mmCIF['_atom_site.id'])):

                if d_mmCIF['_atom_site.label_entity_id'][i_atom_site] not in d_mmCIF['_entity_poly.entity_id']:
                    continue

                element = d_mmCIF['_atom_site.type_symbol'][i_atom_site]

                ## only do heavy atoms
                if element == 'H':
                    continue
                if element not in d_mass.keys():
                    print pdb, d_mmCIF['_atom_site.type_symbol'][i_atom_site]
                    continue

                mass = d_mass[element]
                l_masses += [mass]

                x = float(d_mmCIF['_atom_site.Cartn_x'][i_atom_site])
                y = float(d_mmCIF['_atom_site.Cartn_y'][i_atom_site])
                z = float(d_mmCIF['_atom_site.Cartn_z'][i_atom_site])
                coord = numpy.array([x,y,z,])
                l_coords += [coord]

                center_of_mass += mass*coord

            center_of_mass /= sum(l_masses)

            ##
            ## calculate radius of gyration
            ##
            sum_r = 0
            for i_coord in range(len(l_coords)):
                coord = l_coords[i_coord]
                mass = l_masses[i_coord]
                sq_dist_from_center_of_mass = sum((coord-center_of_mass)**2)
                sum_r += mass*sq_dist_from_center_of_mass
            radius_of_gyration = math.sqrt(sum_r/sum(l_masses))

            print pdb, center_of_mass, radius_of_gyration

            line = '%s %s\n' %(pdb,radius_of_gyration,)
            lines_out += [line]

            fd = open('radius_of_gyration.txt','a')
            fd.write(line)
            fd.close()

            d_radii[pdb] = radius_of_gyration

    ##
    ## write calculated radii of gyration to file
    ##
    lines_out = []
    for pdb,radius_of_gyration in d_radii.items():
        line = '%s %s\n' %(pdb,radius_of_gyration,)
        lines_out += [line]
    fd = open('radius_of_gyration.txt','w')
    fd.writelines(lines_out)
    fd.close()

    return

if __name__ == '__main__':
    main()

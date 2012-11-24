def main(pdb, d_mmCIF, parse_ATOM = True, parse_HETATM = False, query_chain = None):

    import numpy

    d_coords = {}
    l_coords = []
    for i in range(len(d_mmCIF['_atom_site.id'])):
        element = d_mmCIF['_atom_site.type_symbol'][i]
        chain = d_mmCIF['_atom_site.label_asym_id'][i]
        if query_chain:
            if chain != query_chain:
                continue
        res_no = d_mmCIF['_atom_site.label_seq_id'][i]
        alt_id = d_mmCIF['_atom_site.label_alt_id'][i]
        x = float(d_mmCIF['_atom_site.Cartn_x'][i])
        y = float(d_mmCIF['_atom_site.Cartn_y'][i])
        z = float(d_mmCIF['_atom_site.Cartn_z'][i])
        coord = numpy.array([x,y,z,])
        atom_id = d_mmCIF['_atom_site.label_atom_id'][i]

        if parse_ATOM == True and atom_id != 'CA':
            continue

        ## break if multiple models
        if '_atom_site.pdbx_PDB_model_num' in d_mmCIF.keys():
            if int(d_mmCIF['_atom_site.pdbx_PDB_model_num'][i]) > 1:
                break

##        if d_mmCIF['_atom_site.group_PDB'] == 'HETATM':

        if alt_id in [
            'B','C','2',
            'D', ## e.g. 1o2i
            ]:
            ## tmp!!! because of error in 3eft... also error in 3gy7? no error in 3b7i!
            if alt_id == 'B' and d_mmCIF['_atom_site.label_alt_id'][i-1] != 'A':
                if pdb in ['3eft','3gy7',]:
                    pass
                elif pdb in ['3b7i']:
                    continue
                else:
                    print pdb
                    print chain, res_no, 'prev label_alt_id', d_mmCIF['_atom_site.label_alt_id'][i-1]
                    print d_mmCIF['_atom_site.id'][i]
                    stop
            else:
                continue
        elif alt_id not in ['.','A','1',]:
            print pdb
            print alt_id
            print chain, res_no
            print d_mmCIF['_atom_site'][i]
            stop_alt_id

        ##
        ## append
        ##
        if not chain in d_coords.keys():
            d_coords[chain] = {}
        if not res_no in d_coords[chain].keys():
            d_coords[chain][res_no] = {}
        d_coords[chain][res_no][alt_id] = {'atom_id':atom_id}
        ## alpha carbon atoms only...
        if parse_ATOM == True:
            if (
                d_mmCIF['_atom_site.group_PDB'][i] == 'ATOM'
                and
                element == 'C'
                and
                atom_id == 'CA'
                ):
                l_coords += [coord]
        if parse_HETATM == True:
            if (
                d_mmCIF['_atom_site.group_PDB'][i] == 'HETATM'
                and
                d_mmCIF['_atom_site.label_comp_id'][i] != 'HOH'
                ):
                l_coords += [coord]

    return d_coords, l_coords

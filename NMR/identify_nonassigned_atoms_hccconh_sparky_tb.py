fd = open('HCCCONH.shifts', 'r')
lines = fd.readlines()
fd.close()

atoms_expected_all = {
    'A': ['H','HA'],
    'R': ['H','HA',     'HBa','HBb',      'HGa','HGb',            'HDa',      'HE'],
    'D': ['H','HA',     'HBa','HBb'],
    'N': ['H','HA',     'HBa','HBb'],
    'C': ['H','HA',     'HBa','HBb',      'HGa','HGb'],
    'E': ['H','HA',     'HBa','HBb',      'HGa','HGb'],
    'Q': ['H','HA',     'HBa','HBb',      'HGa','HGb'],
    'G': ['H',          'HAa','HAb'],
    'H': ['H','HA',     'HBa','HBb'],
    'I': ['H','HA','HB',                 'HG1a','HG1b','HG2','HD1'],
    'L': ['H','HA',     'HBa','HBb','HG',                        'HD1','HD2'],
    'K': ['H','HA',     'HBa','HBb',     'HGa','HGb',            'HDa','HDb',     'HEa','HEb'],
    'M': ['H','HA',     'HBa','HBb',     'HGa','HGb'],
    'F': ['H','HA',     'HBa','HBb'],
    'P': [    'HA',     'HBa','HBb',     'HGa','HGb',            'HDa','HDb'],
    'S': ['H','HA',     'HBa','HBb'], #HG = hydroxyl
    'T': ['H','HA','HB',                                   'HG2'],
    'W': ['H','HA',     'HBa','HBb'], #more signals???
    'Y': ['H','HA',     'HBa','HBb'],
    'V': ['H','HA','HB',                             'HG1','HG2'],
    }

atoms_observed = []
for col in range(7,len(lines[0]),7):
    atoms_observed.append(lines[0][col:col+7].strip())

prolines = []
residues = []
for line in lines[1:]:
    res = line[:7].strip()[0]
    resno = int(line[:7].strip()[1:])
    residues.append(resno)
    if res == 'P':
        prolines.append(resno)
resminmax = [min(residues), max(residues)]

for line in lines[1:]:
    resno = int(line[:7].strip()[1:])
    if resno+1 in prolines or resno in resminmax:
        continue
    res = line[:7].strip()[0]
    atoms_expected_res = list(atoms_expected_all[res])
    atoms_removed_res = []
    for col in range(len(line)/7-1):
        cs = line[7*(col+1):7*(col+2)].strip()
        if cs != '-':
            atom_observed = atoms_observed[col]
            if atom_observed == 'N':
                continue
            if len(atom_observed) == 1:
                atoms_expected_res.remove(atom_observed)
            if len(atom_observed) == 2:
                if atom_observed == 'HA' and res != 'G':
                    atoms_expected_res.remove(atom_observed)
                else:
                    for suffix in ['','a','b']:
                        if atom_observed+suffix in atoms_expected_res and atom_observed+suffix not in atoms_removed_res:
                            atoms_expected_res.remove(atom_observed+suffix)
                            atoms_removed_res.append(atoms_observed[col][:2])
            if len(atom_observed) == 3:
                if atom_observed[-1] in ['a','b']:
                    atoms_expected_res.remove(atom_observed)
                    if atoms_observed[col][:2] in atoms_expected_res and atom_observed+suffix not in atoms_removed_res:
                        atoms_expected_res.remove(atoms_observed[col][:2])
                        atoms_removed_res.append(atoms_observed[col][:2])
                else:
                    for suffix in ['','a','b']:
                        if atom_observed+suffix in atoms_expected_res and atom_observed+suffix not in atoms_removed_res:
                            atoms_expected_res.remove(atom_observed+suffix)
                            atoms_removed_res.append(atoms_observed[col][:2])
    if not len(atoms_expected_res) == 0:
        print line
        print resno, res, atoms_expected_res
        stop

print 'all peaks assigned'

fd = open('/home/people/tc/1hkc_isomerase/1hkc.pdb','r')
lines1 = fd.readlines()
fd.close()

lines2 = []

d_models2chains = {2:{'A':'B'}}

for line1 in lines1:

    record = line1[:6].strip()
    if record in ['ATOM','HETATM','TER',]:
        if model in d_models2chains.keys():
            chain1 = line1[21]
            chain2 = d_models2chains[model][chain1]
            line2 = line1[:21]+chain2+line1[22:]
        else:
            line2 = line1
        lines2 += [line2]
        continue
    elif record == 'MODEL':
        model = int(line1[6:])
    elif record == 'ENDMDL':
        continue
    else:
        lines2 += [line1]

fd = open('/home/people/tc/1hkc_isomerase/1hkc_modified.pdb','w')
fd.writelines(lines2)
fd.close()

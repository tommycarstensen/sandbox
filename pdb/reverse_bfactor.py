import os
l_fn = os.listdir(os.getcwd())
for fn in l_fn:
    lines2 = []
    if not fn[-4:] == '.pdb':
        continue
##    if fn != '5dfr_A_probe.pdb':
##        continue
    fd = open(fn,'r')
    lines = fd.readlines()
    fd.close()
    for line in lines:
        record = line[:6].strip()
        if record in ['ATOM','HETATM',]:
            line = '%s%6.2f%s' %(line[:60],100.-float(line[60:66]),line[66:],)
        lines2 += [line]
    fd = open(fn,'w')
    fd.writelines(lines2)
    fd.close()

import sys

fn = sys.argv[-2]
chainID = sys.argv[-1]

fd = open(fn,'r')
lines = fd.readlines()
fd.close()

for i in range(len(lines)):
    line = lines[i]
    record = line[:6].strip()
    if record in ['ATOM','HETATM',]:
        lines[i] = line[:21]+chainID+line[22:]

fd = open(fn,'w')
fd.writelines(lines)
fd.close()

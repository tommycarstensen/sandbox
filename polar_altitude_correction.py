hrm = '08060701'

fd = open('%s.hrm' %(hrm),'r')
lines = fd.readlines()
fd.close()

alt2 = int(lines[-1].split()[1])

Continue = True
for i in range(len(lines)):
    line = lines[i]
    if line == '[HRData]\n':
        i1 = i
        n = len(lines)-i-1-1
        Continue = False
        continue
    if Continue == True:
        continue
    else:
        hr = int(line.split()[0])
        alt = float(line.split()[1])
        alt = int(alt-(float(i-i1)/n)*alt2)
        lines[i] = '%3s\t%s\n' %(hr,alt)

fd = open('test.hrm','w')
fd.writelines(lines)
fd.close()

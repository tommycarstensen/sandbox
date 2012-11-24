pdb = '2rh5'

fd = open('blastp_%s.txt' %(pdb),'r')
lines = fd.readlines()
fd.close()

d_seq = {}
for i in range(len(lines)):
    line = lines[i]
    if line[0] == '>':
        s = line[:-1]
        d_seq[s] = ''
    elif line[:5] == 'Sbjct':
##        if not 'alcohol dehydrogenase' in s.lower():
##            continue
        if expect < 0.001:
            seq = line[12:][:line[12:73].index(' ')]
            d_seq[s] += seq
    elif 'Expect' in line and line.split()[5] == 'Expect':
        expect = float(line.split()[7].replace(',',''))

l = []
for s in d_seq.keys():
    if d_seq[s] == '':
        continue
    l += ['%s\n' %(s)]
    l += ['%s\n' %(d_seq[s].replace('-',''))]

fd = open('fasta_%s.txt' %(pdb),'w')
fd.writelines(l)
fd.close()

for s in l:
    print s[:-1]

fd = open('HCCCONH.shifts', 'r')
lines = fd.readlines()
fd.close()

lines_gnuplot = []
for i in range(2,len(lines)):
    if len(lines[i]) > 140:
        for col in range(1,19):
            if lines[i][7*col+6:7*(col+1)] != ' ':
                lines[i] = lines[i][:7*col+6] + lines[i][7*col+6+lines[i][7*col+6:].index(' '):]
    res = lines[i][:7].strip()[0]
    if res == 'P':
        continue
    csN = lines[i][7*19:7*20].strip()
    csH = lines[i][7*1:7*2].strip()
    for col in range(1,19):
        cs = lines[i-1][7*col:7*(col+1)].strip()
        if cs != '-':
            lines_gnuplot.append('%s %s %s\n' %(csN,csH,cs))

fd = open('gnuplot.txt', 'w')
fd.writelines(lines_gnuplot)
fd.close()

print 'sparky2gnuplot completed'
## name coordinates in gnuplot

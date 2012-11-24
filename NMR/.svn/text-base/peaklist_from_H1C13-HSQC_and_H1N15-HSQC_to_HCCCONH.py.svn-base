fd = open('ACBP_H1N15.list', 'r')
linesNH = fd.readlines()[2:]
fd.close()

fd = open('ACBP_H1C13.list', 'r')
linesCH = fd.readlines()[2:]
fd.close()

assHBs = {}
for line in linesCH:
    assCH = line[0:17].strip()
    if assCH[assCH.index('-')+1:assCH.index('-')+3] == 'HB':
        cs = float(line[28:39])
	resno = int(assCH[1:assCH.index('-')-2])
	name = assCH[:assCH.index('-')-2]+assCH[assCH.index('-')+1:]
	if not assHBs.has_key(resno):
	    assHBs[resno] = {name: cs}
	else:
	    assHBs[resno][name] = cs

linesHCCCONH = ['      Assignment         w1         w2         w3  \n\n']
for line in linesNH:
    assNH = line[0:17].strip()
    if assNH[-2:] == 'HN':
        resno = int(assNH[1:assNH.index('-')-1])
	csN = float(line[17:28])
	csH = float(line[28:39])
	if not assHBs.has_key(resno-1):
	    continue
	for assHB in assHBs[resno-1]:
	    csHB = assHBs[resno-1][assHB]
	    assH = assNH[:assNH.index('-')]
	    assHCCCONH = assH+'-'+assHB+'-'+assH[:-1]+'H'
	    lineHCCCONH = '%17s%11.3f%11.3f%11.3f\n' %(assHCCCONH, csN, csHB, csH)
            linesHCCCONH.append(lineHCCCONH)

fd = open('ACBP_HCCCONH.list', 'w')
fd.writelines(linesHCCCONH)
fd.close()

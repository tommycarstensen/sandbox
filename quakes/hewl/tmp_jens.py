## create dict of muts

d_jens = {
    'T4L':{},
    'HEWL':{},
    }

for protein in d_jens.keys():
    fd = open('table_%s.txt' %(protein),'r')
    lines = fd.readlines()
    fd.close()

    for line in lines[1:]:
        l = line.split('\t')
        pdb = l[0]
        if protein == 'T4L':
            mutation = l[-1].strip()
        elif protein == 'HEWL':
            mutation = l[-2].strip()
            if mutation != 'wt':
                mutation = 'wt+'+mutation
        mutation = mutation.replace('wt*','wt+C54T+C97A',)
        mutation = mutation.replace('CORE10','V87I+I100V+M102L+V103I+M106I+V111A+M120Y+L133F+V149I+T152V',)
        mutation = mutation.replace('CORE7','I78V+V87M+L118I+M120Y+L133F+V149I+T152V',)
        mutation = mutation.replace('3xdisulf','I3C+I9C+T21C+C54T+T142C+L164C',)
        if '-' in mutation:
            m = mutation[mutation.index('-')+1:]
            mutation = mutation.replace('+%s' %(m),'',)
            mutation = mutation.replace('-%s' %(m),'',)
        print mutation
        l_mutations = mutation.split('+')
        for mutation in l_mutations:
            if mutation == 'wt':
                continue
            if '-' in mutation:
                print mutation
                stop
            if mutation[1] not in ['1','2','3','4','5','6','7','8','9',]:
                print mutation
                stop
            if mutation[-2] not in ['1','2','3','4','5','6','7','8','9','0',]:
                print mutation
                print l_mutations
                stop
            if 'X' in mutation:
                print pdb
                print mutation
                stio

        d_jens[protein][pdb] = mutation

import pickle

fd = open('d_jens.pickle','w')
pickle.dump(d_jens,fd)
fd.close()

print d_jens

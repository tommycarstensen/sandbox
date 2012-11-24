fd = open('d_jens.txt','r')
s = fd.read()
fd.close()

d_jens = eval(s)

import pickle

fd = open('d_jens.pickle','w')
pickle.dump(d_jens,fd)
fd.close()

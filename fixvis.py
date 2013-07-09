# Fix defects
import numpy as np

# Look for any surface opportunities
for j in range(len(Ddefect)):
    i = Ddefect[j,0]; #print i
    l = Ddefect[j,1]; #print i
    isurf = np.argwhere(nni[i]<0)
    if len(isurf)>0:
        print 'Possible surface opportunities for donor defects ...'
        print i,l,len(isurf)
        print nni[i]
        print nnitol[i]
        print nnltoi[i]
        print nni[l]
        print nnitol[l]
        print nnltoi[l]
               
for j in range(len(Adefect)):
    i = Adefect[j,0]; #print i
    l = Adefect[j,1]; #print i
    isurf = np.argwhere(nni[i]<0)
    if len(isurf)>0:
        print 'Possible surface opportunities for acceptor defects ...'
        print i,l,len(isurf)
        print nni[i]
        print nnitol[i]
        print nnltoi[i]
        print nni[l]
        print nnitol[l]
        print nnltoi[l]
 
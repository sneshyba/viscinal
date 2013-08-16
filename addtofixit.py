import numpy as np
import copy

#nnitol_new = copy.deepcopy(nnitol)
nnitol_new = nnitol

# Fixing the missing H2 cases
for i in range (nR):
        
    # Find index of H2 for ith residue
    H2ofi=np.argwhere(nnitol_new[i]==2)
    
    # See if we didn't find H2
    if len(H2ofi) == 0:
        
        # OK, let's find a good place for it
        test = np.argwhere(nni[i]==-1)
        if len(test) > 0:
            
            # Find a spot that's not being used
            for j in range(len(test)): 
            
                if nnitol_new[i,test[j]] == 0:
                    nnitol_new[i,test[j]]=2
                    print nni[i], nnitol[i], nnitol_new[i]
                    break
        else:
            print "Internal inconsistency ..."
            f = np.sqrt(-1.0)

# Fixing the missing H1 cases
for i in range (nR):

    # Find index of H1 for ith residue
    H1ofi=np.argwhere(nnitol_new[i]==1)
    
    # See if we didn't find H2
    if len(H1ofi) == 0:
        
        # OK, let's find a good place for it
        test = np.argwhere(nni[i]==-1)
        if len(test) > 0:
            
            # Find a spot that's not being used
            for j in range(len(test)): 
            
                if nnitol_new[i,test[j]] == 0:
                    nnitol_new[i,test[j]]=1
                    print nni[i], nnitol[i], nnitol_new[i]
                    break
        else:
            print "Internal inconsistency ..."
            f = np.sqrt(-1.0)

    # Use these vectors to construct the rotation matrix


import numpy as np
import vstuff as vs; reload(vs)
import random
import pdb

# Get an initial list of donors and acceptors
nnitol, nnltoi = vs.donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print len(Ddefect), len(Adefect)
print Ddefect
print Adefect

# Check for surface defects
nnitol,nnltoi = vs.fixsurface(nni,nnitol,nnltoi)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print len(Ddefect), len(Adefect)

# Checking ...
for i in range (nR):
    test = np.argwhere(nnitol[i]==0)
    if len(test)>2:
        print "There's a problem at i, nni, nnitol =", i, nni[i], nnitol[i]

#Specifiy number of times through the defect propogation loop
nprop=500

# Propagate a defect and check
for iprop in range(nprop):
    if len(Ddefect) > 0:
        print "****"
        print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
        m = random.randint(0,len(Ddefect)-1)
        i = Ddefect[m,0]
        l = Ddefect[m,1]
        print "working on ", i, l
        klofi = np.squeeze(np.argwhere(nni[i]==l))
        kiofl = np.squeeze(np.argwhere(nni[l]==i))
    
        # Figure out which Hydrogen of i (1 or 2) is pointing to l
        Hofi = nnitol[i,klofi]; #print Hofi    

        # Decide on a new nearest neighbor to point this Hydrogen to
        first=np.squeeze(np.argwhere(nnitol[i]==0)[0]) 
        second=np.squeeze(np.argwhere(nnitol[i]==0)[1])
        zeroorone = random.randint(0,1) #sloppy way to get random index
        if zeroorone == 0:
            klofip = first
        else:
            klofip = second
        print "first, second, klofip = ", first, second, klofip
        lp = nni[i,klofip]; #print lp
        kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

        # Point a lone pair to l instead
        nnitol[i,klofi] = 0
        nnltoi[l,kiofl] = 0
    
        # Point the offending hydrogen of i to next victim (nearest neighbor)
        nnitol[i,klofip] = Hofi
	nnltoi[lp,kioflp] = Hofi
    
        # Check for surface defects
        nnitol,nnltoi = vs.fixsurface(nni,nnitol,nnltoi)
        Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
        print "Ending iprop =   ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)

#Propagate Adefects
for iprop in range(nprop):
    if len(Adefect) > 0:
        print "****"
        
        print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
        m = random.randint(0,len(Adefect)-1)
        i = Adefect[m,0]
        l = Adefect[m,1]
        print "working on ", i, l
        klofi = np.squeeze(np.argwhere(nni[i]==l))
        kiofl = np.squeeze(np.argwhere(nni[l]==i))
        

        # Decide on a new nearest neighbor to point this lone pair to
        knonzeros=np.squeeze(np.argwhere(nnitol[i]!=0))
        ikrandom = random.randint(0,len(knonzeros)-1) #sloppy way to get random index
        klofip = knonzeros[ikrandom]
        print "nnitol[i], klofip = ", nnitol[i], klofip
        lp = nni[i,klofip]; #print lp
        kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

        # Point a hydrogen to l instead
        Hofi=nnitol[i,klofip]
        nnitol[i,klofi] = Hofi
        nnltoi[l,kiofl] = Hofi
    
        # Point the offending lone pair of i to next victim (nearest neighbor)
        nnitol[i,klofip] = 0
	nnltoi[lp,kioflp] = 0
        
        #pdb.set_trace()
        # Check for surface defects
        nnitol,nnltoi = vs.fixsurface(nni,nnitol,nnltoi)
        Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
        print "Ending iprop =   ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
        
        
# Assign missing H1s and H2s to slots
nnitol_old = copy.deepcopy(nnitol)

# Fixing the missing H2 cases
for i in range (nR):
        
    # Find index of H2 for ith residue
    H2ofi=np.argwhere(nnitol[i]==2)
    
    # See if we didn't find H2
    if len(H2ofi) == 0:
        
        # OK, let's find a good place for it
        test = np.argwhere(nni[i]==-1)
        if len(test) > 0:
            
            # Find a spot that's not being used
            for j in range(len(test)): 
            
                if nnitol[i,test[j]] == 0:
                    nnitol[i,test[j]]=2
                    #print nni[i], nnitol[i], nnitol[i]
                    break
        else:
            print "Internal inconsistency ..."
            f = np.sqrt(-1.0)

# Fixing the missing H1 cases
for i in range (nR):

    # Find index of H1 for ith residue
    H1ofi=np.argwhere(nnitol[i]==1)
    
    # See if we didn't find H2
    if len(H1ofi) == 0:
        
        # OK, let's find a good place for it
        test = np.argwhere(nni[i]==-1)
        if len(test) > 0:
            
            # Find a spot that's not being used
            for j in range(len(test)): 
            
                if nnitol[i,test[j]] == 0:
                    nnitol[i,test[j]]=1
                    #print nni[i], nnitol[i], nnitol[i]
                    break
        else:
            print "Internal inconsistency ..."
            f = np.sqrt(-1.0)

# Checking ...
for i in range (nR):
    test = np.argwhere(nnitol[i]==0)
    if len(test)>2:
        print "There's a problem at i, nni, nnitol =", i, nni[i], nnitol[i]

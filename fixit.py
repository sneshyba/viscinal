import numpy as np
import vstuff as vs; reload(vs)
import random

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

# Propagate a defect and check
for iprop in range(500):
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

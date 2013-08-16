import numpy as np
import vstuff as vs; reload(vs)
import random
import pdb
import copy

# Get an initial list of donors and acceptors
nnitol, nnltoi = vs.donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift)
nnitol_old = copy.deepcopy(nnitol)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print "Starting with defects:", len(Ddefect), len(Adefect)

# Check for surface defects
nnitol,nnltoi = vs.fixsurface(nni,nnitol,nnltoi)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print "After initial surface fix:", len(Ddefect), len(Adefect)

#Specify number of times through the defect propogation loop
nprop=1000

# Propagate a defect and check
icount = 0
for iprop in range(nprop):
    if len(Ddefect) > 0:
        icount += 1
        #print "****"
        #print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
        m = random.randint(0,len(Ddefect)-1)
        i = Ddefect[m,0]
        l = Ddefect[m,1]
        #print "working on ", i, l
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
        #print "first, second, klofip = ", first, second, klofip
        lp = nni[i,klofip]; #print lp
        kioflp = np.squeeze(np.argwhere(nni[lp]==i)); #print kioflp

        # Point a lone pair to l instead
        nnitol[i,klofi] = 0
        nnltoi[l,kiofl] = 0
    
        # Point the offending hydrogen of i to next victim (nearest neighbor)
        nnitol[i,klofip] = Hofi
	nnltoi[lp,kioflp] = Hofi
    
        # Check for surface defects
        #problem1 = vs.findthreezeros(nni,nnitol)
        nnitol,nnltoi = vs.fixsurface(nni,nnitol,nnltoi)
        Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
        #problem2 = vs.findthreezeros(nni,nnitol)
        #if (problem1==False) & (problem2==True):
        #    print "Fixsurface D is at fault"
            
        #    print "... at i = ", i
        #print "Ending iprop =   ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
print "After fixing Ddefects:", len(Ddefect), len(Adefect), "which took", icount, "iterations"

#Propagate Adefects
icount = 0
for iprop in range(nprop):
    if len(Adefect) > 0:
        icount += 1
        #print "****"
        
        #print "Starting iprop = ", iprop, " w/NDdefect, NAdefect = ", len(Ddefect), len(Adefect)
        m = random.randint(0,len(Adefect)-1)
        i = Adefect[m,0]
        l = Adefect[m,1]
        #print "working on ", i, l
        klofi = np.squeeze(np.argwhere(nni[i]==l))
        kiofl = np.squeeze(np.argwhere(nni[l]==i))  

        # Decide on a new nearest neighbor to point this lone pair to
        knonzeros=np.squeeze(np.argwhere(nnitol[i]!=0))
        ikrandom = random.randint(0,len(knonzeros)-1) #sloppy way to get random index
        klofip = knonzeros[ikrandom]
        #print "nnitol[i], klofip = ", nnitol[i], klofip
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
        #problem1 = vs.findthreezeros(nni,nnitol)
        nnitol,nnltoi = vs.fixsurface(nni,nnitol,nnltoi)
        Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
        #problem2 = vs.findthreezeros(nni,nnitol)
        #if (problem1==False) & (problem2==True):
        #    print "Fixsurface A is at fault"
print "After fixing Adefects:", len(Ddefect), len(Adefect), "which took", icount, "iterations"

# Just checking
for i in range (nR):
    test = np.argwhere(nnitol[i]==0)
    if len(test)!=2:
        print "There's a problem at i, nni, nnitol =", i, nni[i], nnitol[i]



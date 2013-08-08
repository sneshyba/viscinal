import numpy as np


def findnbad(nni):
    nR, nk = nni.shape
    nbad = 0
    for i in range (nR):
        if np.size(np.where(nni[i]<0))>0:
            nbad = nbad+1
    #print "nbad = ", nbad
    return nbad

def finddefects(nni,nnltoi,nnitol):
    nR, nk = nni.shape; #print nR, nk
    Ddefect=np.zeros((nR*10,2)).astype(np.int32); nDdefect = 0
    Adefect=np.zeros((nR*10,2)).astype(np.int32); nAdefect = 0
    for i in range(nR):
        iDdefect = 0
        iAdefect = 0
        for k in range(4):
            if (nni[i,k]>=0):
                if (nnltoi[i,k]!=0) & (nnitol[i,k]!=0):
                    Ddefect[nDdefect,0]=i
                    Ddefect[nDdefect,1]=nni[i,k]
                    nDdefect += 1
                elif (nnltoi[i,k]==0) & (nnitol[i,k]==0):
                    Adefect[nAdefect,0]=i
                    Adefect[nAdefect,1]=nni[i,k]
                    nAdefect += 1
        #print i, nni[i], nnltoi[i], nnitol[i], Ddefect, Adefect
    Ddefect_ret = Ddefect[0:nDdefect]
    Adefect_ret = Adefect[0:nAdefect]
    return Ddefect_ret, Adefect_ret

def donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift):

    # Pre-allocate the nnitol and nnltoi arrays
    nR, nk = nni.shape; #print nR, nk
    nnitol=np.zeros((nR,4), dtype='int32')
    nnltoi=np.zeros((nR,4), dtype='int32')

    # Minimium projection required in order call this a donor
    HBproject = 2.7
    
    # Loop over all the residues
    for i in range(nR):
    
        #Is l(k) donating to i?
        for k in range(4):
            if nni[i,k]>=0:
                l= nni[i,k]
                viO=xyzO[i]
                vkO=xyzO[l]+xyzshift[i,k]
                vkH1=xyzH1[l]+xyzshift[i,k]
                vkH2=xyzH2[l]+xyzshift[i,k]
                vOO=viO-vkO
                vH1O=vkH1-vkO
                vH2O=vkH2-vkO
                #print i,l,np.dot(vOO, vH1O), np.dot(vOO, vH2O)
                if np.dot(vOO, vH1O)>HBproject:
                    nnltoi[i,k]=1
                if np.dot(vOO, vH2O)>HBproject:
                    nnltoi[i,k]=2
    
        #Is i donating to l(k)?
        for k in range(4):
	    if nni[i,k]>=0:
                l= nni[i,k]
                viO=xyzO[i]
                vkO=xyzO[l]+xyzshift[i,k]
                viH1=xyzH1[i]
                viH2=xyzH2[i]
                vOO=vkO-viO
            	vH1O=viH1-viO
                vH2O=viH2-viO
                #print i,l,np.dot(vOO, vH1O), np.dot(vOO, vH2O)
                if np.dot(vOO, vH1O)>HBproject:
                    nnitol[i,k]=1
                if np.dot(vOO, vH2O)>HBproject:
                    nnitol[i,k]=2

    return nnitol, nnltoi
    
    
def fixsurface(nni,nnitol,nnltoi):


    # Check for defects
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    #print len(Ddefect), "Ddefects: "
    #print Ddefect
    #print len(Adefect), "Adefects: "
    #print Adefect


    #fix Ddefect with a "-1" nearest neighbor:
    #print "Now we are fixing Ddefects ..."
    for m in range(len(Ddefect)):
    
        # Pull out the index to the next residue that has a donor defect
        i = Ddefect[m,0]
        l = Ddefect[m,1]

        # See if this residue has any "-1" which means doesn't have a nearest neighbor    
        test = np.size(np.where(nni[i]<0))
    
        # If this defective residue is missing a nearest neighbor, point its hydrogen toward the space
        if test>0:
        
            # Get the positions in nni with the -1, and mutual pointer positions
            kzeroofi = np.squeeze(np.argwhere(nni[i]<0))
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))

            # Fix it by changing nnitol
            temp = nnitol[i,klofi] # Save which of i's Hydrogens is the donor defective guy
            nnitol[i,klofi] = 0 # Point i's lone pair to l
            nnitol[i,kzeroofi] = temp # Point i's Hydrogen to what was -1 
            nnltoi[l,kiofl] = 0 # Confirms that i no longer points its H to l
        
            # Report out
            print 'fixed ', i, l, ' donor defect'
        
    # Check for defects
    Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    #print len(Ddefect), "Ddefects: "
    #print Ddefect
                
    #Fix Adefect with a "-1" as a nearest neighbor(as long as nnltoi[i] does not have four 0's)
    #Same logic as fixing Ddefect with "-1" as nearest neighbor

    #print "Now we are fixing Adefects ..."
    for n in range(len(Adefect)):
        i = Adefect[n,0]
        l = Adefect[n,1]
        position = np.where(nni[i]<0)
        test = np.size(position)
        if test>0:
            klofi = np.squeeze(np.argwhere(nni[i]==l))
            kiofl = np.squeeze(np.argwhere(nni[l]==i))
        
            # Figure out which hydrogen is available (if any)
            timesH1isused = np.size(np.argwhere(nnitol[i]==1))
            timesH2isused = np.size(np.argwhere(nnitol[i]==2))
            whichone = 0
            if timesH1isused == 0:
                whichone = 1
            elif timesH2isused == 0:
                whichone = 2
            if whichone != 0:
                #print 'before: ', i,l,nni[i],nni[l],nnitol[i],nnltoi[l]
                nnitol[i,klofi] = whichone
                nnltoi[l,kiofl] = whichone
                #print 'after:  ', i,l,nni[i],nni[l],nnitol[i],nnltoi[l]
                            # Report out
                print 'fixed ', i, l, ' acceptor defect'
 

    #Ddefect, Adefect = finddefects(nni,nnltoi,nnitol)
    #print len(Adefect), "Adefects: "
    #print Adefect
    
    return nnitol,nnltoi

def getrefcoords(xyzO,xyzH1,xyzH2):
    lOH = np.sqrt(np.sum((xyzO-xyzH1)**2))
    lHH = np.sqrt(np.sum((xyzH2-xyzH1)**2))
    ycoord = lHH/2
    theta = np.arccos(ycoord/lOH)
    zcoord = lOH*np.sin(theta)
    phi = 2*(np.pi/2-theta)
    #print theta*180/np.pi, phi*180/np.pi
    vH1 = np.array([0.,-ycoord,-zcoord])
    vH2 = np.array([0.,ycoord,-zcoord])
    vO = np.array([0.,0.,0.])
    return vO, vH1, vH2
    
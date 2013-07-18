# Check for defects
nnitol, nnltoi = vs.donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print len(Ddefect), "Ddefects: "
print Ddefect



for m in range(len(Ddefect)):
    
    # Pull out the index to the next residue that has a donor defect
    i = Ddefect[m,0]
    l = Ddefect[m,1]

    # See if this residue has any "-1" which means doesn't have a nearest neighbor    
    test = size(np.where(nni[i]<0))
    
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
        
# Check for defects
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print len(Ddefect), "Ddefects: "
print Ddefect


        
        # Change nnitol[0,2] to 0, nnitol[0,3] to 2
        # Change nnltoi[68,3] to 0
        
        # Change nnitol[48,2] to 0, nnitol[48,3] to 2
        # Change nnltoi[116,2] to 0

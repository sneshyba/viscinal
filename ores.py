import numpy as np
# Orient residues in a desired config

# Get water residue in the reference configuration
import vstuff as vs; reload(vs); i=10; xyzO_ref, xyzH1_ref, xyzH2_ref = vs.getrefcoords(xyzO[i],xyzH1[i],xyzH2[i])

# Use nni and nnitol to construct the bisector, 21, and normal vectors

for i in range (nR):
    
    #Reference residue vectors
    v21not=vH2-vH1
    vBnot=(2*vO)-(vH1+vH2) #vO, vH1, vH2 from getrefcoords in vstuff...not sure if I should use xyzO_ref, etc. instead?
    vNnot=np.cross(v21not,vBnot)
    
    #ith residue vectors
    #Find index of H1 and H2 for ith residue
    H2ofi=np.squeeze(np.argwhere(nnitol[i]==2))
    H1ofi=np.squeeze(np.argwhere(nnitol[i]==1))
    
    #Find corresponding residues that H1 and H2 point at
    k2ofi=nni[i,H2ofi]
    k1ofi=nni[i,H1ofi]
    
    #bisector,21, and normal vectors
    v21=xyzO[k2ofi]-xyzO[k1ofi]
    vB=(2*xyzO[i])-(xyzO[k1ofi]+xyzO[k2ofi])
    vN=np.cross(v21,vB)


    # Use these vectors to construct the rotation matrix
    
    
    
    


# Use the rotation matrix to rotate the reference configuration into the desired configuration (ending up with updated xyzH1, and xyzH2)
import numpy as np
# Orient residues in a desired config

# Get water residue in the reference configuration
import vstuff as vs; reload(vs)
i=15; vO, vH1, vH2 = vs.getrefcoords(xyzO[i],xyzH1[i],xyzH2[i])

# Use nni and nnitol to construct the bisector, 21, and normal vectors

#Reference residue vectors
v21not=vH2-vH1; #v21not = v21not/np.linalg.norm(v21not)
vBnot=(2*vO)-(vH1+vH2); #vBnot = vBnot/np.linalg.norm(vBnot)
vNnot=np.cross(v21not,vBnot); #vNnot = vNnot/np.linalg.norm(vNnot)


for i in range (nR):
    
    
    #ith residue vectors
    #Find index of H1 and H2 for ith residue
    H2ofi=np.squeeze(np.argwhere(nnitol[i]==2))
    H1ofi=np.squeeze(np.argwhere(nnitol[i]==1))
    
    #Find corresponding residues that H1 and H2 point at
    k2ofi=nni[i,H2ofi]
    k1ofi=nni[i,H1ofi]
    
    # if neither k1ofi or k2ofi comes up negative, proceed
    if ((k2ofi!=-1) & (k1ofi!=-1)):
    
        #bisector,21, and normal vectors
        v21=xyzO[k2ofi]-xyzO[k1ofi]
        vB=(2*xyzO[i])-(xyzO[k1ofi]+xyzO[k2ofi])
        vN=np.cross(v21,vB)

    else:
        print i,nni[i],nnitol[i]

    # Use these vectors to construct the rotation matrix
    
    
    
    


# Use the rotation matrix to rotate the reference configuration into the desired configuration (ending up with updated xyzH1, and xyzH2)
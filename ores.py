# Orient residues in a desired config

# Get water residue in the reference configuration
import vstuff as vs; reload(vs); i=10; xyzO_ref, xyzH1_ref, xyzH2_ref = vs.getrefcoords(xyzO[i],xyzH1[i],xyzH2[i])

# Use nni and nnitol to construct the bisector, 21, and normal vectors

# Use these vectors to construct the rotation matrix

# Use the rotation matrix to rotate the reference configuration into the desired configuration (ending up with updated xyzH1, and xyzH2)
import vstuff as vs; reload(vs)

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
filename = 'spc_4_4_6.pdb'; nx = 4; ny = 4; nz = 6
#filename = 'spc_4_4_2.pdb'; nx = 4; ny = 4; nz = 2
viscinaldir = 'y'; nycel=1
xyzO, xyzH1, xyzH2, viscinaldir, shift, vshift, xbox, ybox, zbox = vs.loadit(filename, nx, ny, nz, viscinaldir, nycel)

# Get the nearest neighbor index and nearest-neighbor shift array
nni,xyzshift = vs.getnni(xyzO,shift) 

# Check for any defects based on the present xyzO, xyzH1, and xyzH2 arrays
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO,xyzH1,xyzH2,xyzshift)

# Make new nnitol and,nnltoi arrays that correct for defects
nnitol_fixed = vs.fixit(nni, nnitol, nnltoi, 1000)

# Use new nnitol array to reconstruct the xyzO, xyzH1, and xyzH2 arrays
xyzO_new, xyzH1_new, xyzH2_new = vs.reconstructit(xyzO, xyzH1, xyzH2, nni, nnitol_fixed, xyzshift)

# Check to make sure these are OK
nnitol, nnltoi, Ddefect, Adefect = vs.checkfordefects(nni,xyzO_new,xyzH1_new,xyzH2_new,xyzshift)

# Rotate reconstructed xyzO, xyzH1, and xyzH2 arrays to create the viscinal slab
xyzO_rot, xyzH1_rot, xyzH2_rot, xboxp, yboxp, zboxp = vs.rotateit(xyzO_new, xyzH1_new, xyzH2_new, viscinaldir, shift, vshift, xbox, ybox, zbox)

# Save the new slab
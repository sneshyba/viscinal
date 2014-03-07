import vstuff as vs; reload(vs)
import numpy as np
import copy

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
#filename = 'spc_4_4_6_v_withdefects.pdb'; xbox=17.9629248; ybox=29.3333332974; zbox=23.3345234043
filename = 'spc_4_4_6_v_withdefects_fixed.pdb'; xbox=17.9629248; ybox=29.3333332974; zbox=23.3345234043
namestem = filename.find('.pdb')

# Specify which is the exposed surface, and load the slab
viscinaldir = 'y'; nycel=0
xyzO, xyzH1, xyzH2, shift, structure = vs.loaditnew(filename, xbox, ybox, zbox, viscinaldir)
slab = vs.slab(filename, structure, xyzO, xyzH1, xyzH2, xbox, ybox, zbox)

# Get  nearest neighbor and defect information
nni,xyzshift = vs.getnni(xyzO,shift)
nnitol, nnltoi = vs.getnnitoletc(nni,xyzO,xyzH1,xyzH2,xyzshift)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print "Number of defects = ", len(Ddefect), len(Adefect)

# Save the surface and bulk slab subsets
surfacelist, dum = np.where(nni==-1)
slab_surface = copy.deepcopy(slab)
slab_surface.filename = filename[0:namestem]+'_surface.pdb'
slab_surface.keep(surfacelist)
slab_surface.saveit()
slab_bulk = copy.deepcopy(slab)
slab_bulk.filename = filename[0:namestem]+'_bulk.pdb'
slab_bulk.cullout(surfacelist)
slab_bulk.saveit()

# Save Ddefects
Ddefectlist = np.unique(Ddefect[:,0])
slab_Ddefect = copy.deepcopy(slab)
slab_Ddefect.filename = filename[0:namestem]+'_Ddefect.pdb'
slab_Ddefect.keep(Ddefectlist)
slab_Ddefect.saveit()

# Save Adefects
Adefectlist = np.unique(Adefect[:,0])
slab_Adefect = copy.deepcopy(slab)
slab_Adefect.filename = filename[0:namestem]+'_Adefect.pdb'
slab_Adefect.keep(Adefectlist)
slab_Adefect.saveit()

# Save instances where there are three -1's in nni
a1,b1 = vs.count_unique(surfacelist)
a1[np.argwhere(b1==1)]
a2,b2 = vs.count_unique(surfacelist)
a2[np.argwhere(b2==2)]
a3,b3 = vs.count_unique(surfacelist)
a3[np.argwhere(b3==3)]
a4,b4 = vs.count_unique(surfacelist)
a4[np.argwhere(b4==4)]

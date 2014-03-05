import vstuff as vs; reload(vs)
import numpy as np
import copy

# Get the spc slab (makes new xyzO, xyzH1, xyzH2, and shift-related arrays)
filename = 'spc_4_4_6_v_withdefects.pdb'; xbox=17.9629248; ybox=29.3333332974; zbox=23.3345234043

# Naming the output file
dum = filename.find('.pdb')
surfacefilename = filename[0:dum]+'_surface.pdb'
bulkfilename = filename[0:dum]+'_bulk.pdb'

# Specify which viscinal surface to generate, and load the slab
viscinaldir = 'y'; nycel=0
xyzO, xyzH1, xyzH2, shift, structure = vs.loaditnew(filename, xbox, ybox, zbox, viscinaldir)
slab = vs.slab(filename, structure, xyzO, xyzH1, xyzH2, xbox, ybox, zbox)

# Get some nearest neighbor information
nni,xyzshift = vs.getnni(xyzO,shift)

# Get surface list and bulk lists
surfacelist, dum = np.where(nni==-1)
alllist = range(slab.NR)
bulklist = np.squeeze(np.argwhere(~np.in1d(alllist,surfacelist)))

# Get slabs of bulk and surface residues
slab_surface = copy.deepcopy(slab)
slab_surface.filename = surfacefilename
slab_surface.keep(surfacelist)
slab_surface.saveit()

slab_bulk = copy.deepcopy(slab)
slab_bulk.filename = bulkfilename
slab_bulk.keep(bulklist)
slab_bulk.saveit()



import numpy as np
import vstuff as vs; reload(vs)

# Get an initial list of donors and acceptors
nnitol, nnltoi = vs.donorsandacceptors(nni,xyzO_new,xyzH1_new,xyzH2_new,xyzshift)
#nnitol, nnltoi = vs.donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print "Found defects:", len(Ddefect), len(Adefect)

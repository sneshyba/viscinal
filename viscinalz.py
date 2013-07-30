import numpy as np
import vstuff as vs; reload(vs)
import copy
import Bio
import Bio.PDB

# This is the start of the program

# These are cell dimensions
xcel = 4.4907312
ycel = 7.7781746
zcel = 3.6666666

# Read in the pdb structure & specify the box size
parser = Bio.PDB.PDBParser()
pdb1hlw = parser.get_structure('pdb', 'spc_4_4_2.pdb'); xbox = xcel*4; ybox = ycel*4; zbox = zcel*2
#pdb1hlw = parser.get_structure('pdb', 'spc_4_4_6.pdb'); xbox = xcel*4; ybox = ycel*4; zbox = zcel*6

shift = np.array([\
        [ xbox,       0,        0      ], \
        [  0,        ybox+10,   0,     ], \
        [  0,        ycel*1,     zbox    ]])
                
# A clumsy way to count the number of residues
nR = 0
for model in pdb1hlw:
    for chain in model:
        j=0
        for residue in chain:
            for atom in residue:
                j = j+1
nR = j/3
print "nR = ", nR

# Allocate space for the various arrays
xyz=np.zeros((nR*3,3))
xyzO=np.zeros((nR,3))
xyzH1=np.zeros((nR,3))
xyzH2=np.zeros((nR,3))

# Get all the coordinates
for model in pdb1hlw:
    for chain in model:
        j=0
        for residue in chain:
            for atom in residue:
                temp1 = atom.get_coord()
                #print temp1, j
                xyz[j] = temp1
                j = j+1

# Sort the coordinates into O, H1, H2
i=-1
for j in range(nR*3):
    test=np.mod(j,3)
    if (test == 0):
        i = i+1
        xyzO[i]=xyz[j]
    elif (test ==1):
        xyzH1[i]=xyz[j]
    elif (test == 2):
        xyzH2[i]=xyz[j]
        
# Calculate O-O distances
dist=np.zeros((nR,nR))
for i in range(nR):
    xi=xyzO[i][0]; yi=xyzO[i][1]; zi=xyzO[i][2]
    for j in range(nR):
        xj=xyzO[j][0]; yj=xyzO[j][1]; zj=xyzO[j][2] 
        dist[i,j]=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
        #print i,j
        
# Construct an initial nearest neighbor list, and arrays for related information 
OOdist = 3.0
nni=np.zeros((nR,4),dtype='int32')
nnd=np.zeros((nR,4))
xyzshift=np.zeros((nR,4,3))
for i in range (nR):
    onecol = np.argsort(dist[i])
    x= dist[i,onecol[1:5]]
    nnd[i,]=x
    nni[i,]=onecol[1:5]
    for k in range(4):
        if nnd[i,k]>=OOdist:
            nni[i,k]=-1  # This is flag saying there is a missing nearest neighbor

# Check out this list
print "nbad before searching across periodic boundaries = ", vs.findnbad(nni)

# Find nearest neighbors across periodic boundaries & fix the nearest neighbor list accordingly
for i in range (nR):
    if np.size(np.where(nni[i]<0))>0:
        xi=xyzO[i][0]; yi=xyzO[i][1]; zi=xyzO[i][2]
        for k in range (4):
            if nni[i,k]<0:
                for j in range (nR):
                    
                    # Look across x
                    xj=xyzO[j][0]-shift[0][0]; yj=xyzO[j][1]-shift[0][1]; zj=xyzO[j][2]-shift[0][2]; disttest1=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    xj=xyzO[j][0]+shift[0][0]; yj=xyzO[j][1]+shift[0][1]; zj=xyzO[j][2]+shift[0][2]; disttest2=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5

                    # Look across y
                    xj=xyzO[j][0]-shift[1][0]; yj=xyzO[j][1]-shift[1][1]; zj=xyzO[j][2]-shift[1][2]; disttest3=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    xj=xyzO[j][0]+shift[1][0]; yj=xyzO[j][1]+shift[1][1]; zj=xyzO[j][2]+shift[1][2]; disttest4=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    
                    # Look across z
                    xj=xyzO[j][0]-shift[2][0]; yj=xyzO[j][1]-shift[2][1]; zj=xyzO[j][2]-shift[2][2]; disttest5=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    xj=xyzO[j][0]+shift[2][0]; yj=xyzO[j][1]+shift[2][1]; zj=xyzO[j][2]+shift[2][2]; disttest6=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)**.5
                    
                    # Identify which cross-boundary search resulted in a nearest neighbor catch
                    T1 = disttest1<OOdist
                    T2 = disttest2<OOdist
                    T3 = disttest3<OOdist
                    T4 = disttest4<OOdist
                    T5 = disttest5<OOdist
                    T6 = disttest6<OOdist                   
                    whichone = np.where([T1,T2,T3,T4,T5,T6])
                    test = np.size(whichone)
                    
                    # Record the shift information; this logic must be consistent with other shift logic in this loop
                    if (test==0):
                        xyzshift[i,k] = np.array([0.,0.,0.])
                    elif T1:
                        xyzshift[i,k] = -shift[0] 
                    elif T2:
                        xyzshift[i,k] = shift[0]
                    elif T3:
                        xyzshift[i,k] = -shift[1]
                    elif T4:
                        xyzshift[i,k] = shift[1]
                    elif T5:
                        xyzshift[i,k] = -shift[2]
                    elif T6:
                        xyzshift[i,k] = shift[2]
                        
                        
                    # Record nearest neighbor information in the nni matrix, taking care to eliminate redundancies
                    if test>0:
                        if test>1:
                            print "Caught a double boundary case; the algorithm may not be valid"
                        taken = 0
                        for kp in range(4):
                            jp = nni[i,kp]
                            if jp==j:
                                #print "already taken ..."
                                taken = 1
                        if taken==0:
                            temp = copy.deepcopy(nni[i])
                            nni[i,k]=j
                            #print i, temp, "was replaced by", nni[i], "because of ", xyzshift[i,k] 
                            break

# Check out this list
print "nbad after searching across periodic boundaries = ", vs.findnbad(nni)

# Find defects
nnitol, nnltoi = vs.donorsandacceptors(nni,xyzO,xyzH1,xyzH2,xyzshift)
Ddefect, Adefect = vs.finddefects(nni,nnltoi,nnitol)
print len(Ddefect), "Ddefects: "
print Ddefect
print len(Adefect), "Adefects: "
print Adefect 
print "nnltoi", nnltoi
print "nnitol", nnitol

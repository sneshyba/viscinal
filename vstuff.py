import numpy as np; reload(np)


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
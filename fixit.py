for m in range(len(Ddefect)):
    n = Ddefect[m,0]
    test = size(np.where(nni[n]<0))
    if test>0:
        # Fix it by changing nnitol
        print n, test
        
        # Change nnitol[48,2] to 0, nnitol[48,3] to 2
        # Change nnltoi[116,2] to 0
import numpy

load("helpfunctions.sage")

def isDiagonal(A):
    # assumes A is triangular
    for i in range(A.ncols()):
        for j in range(i+1, A.ncols()):
            if A[i,j] != 0:
                return False
    return True

def HNFs(Delta, r):
    div = divisors(Delta)
    diags = []
    # compute all sorted diagonals that multiply to Delta
    for I in itertools.product(div, repeat=r):
        if prod(I) == Delta and tuple(sorted(I)) == I:
            diags.append(I)
            
    H = []
            
    # successively augment with columns of the form
    # [something<d[n], d[n], 0...0]
    for d in diags:
        I = [d[0]]+(r-1)*[0]
        D = matrix(ZZ, r, 1, I)
        matrixStack = [D]
        while matrixStack:
            D = matrixStack.pop()
            n = D.ncols()
            isDiag = isDiagonal(D)
            # if all other columns are just multiples of unit vectors, then the entries of this column can be reduced
            if isDiag: entries = list(range(0, d[n]//2+1))
            else: entries = list(range(0, d[n]))
            for I in itertools.product(entries, repeat=n):
                # gcd condition from sorting diagonal
                if gcd(I[n-1], d[n]) < D[n-1, n-1]: continue
                # can sort last column if diagonal is (1,...,1,Delta)
                if d[n]==Delta and tuple(sorted(I)) != I: continue
                J = list(I) + [d[n]] + (r-n-1)*[0]
                C = matrix(ZZ, r, 1, J)
                E = D.augment(C)
                if n+1 == r:
                    H.append(E)
                else:
                    matrixStack.append(E)

    H = remove(H)

    return H

# set values for which you want to calculate all HNFs
r = 3
startDelta = 2
endDelta = 50
for Delta in range(startDelta, endDelta+1):
    L = HNFs(Delta, r)
    #print to file
    fileName = "..data/hnfs/"+str(r)+"_"+str(Delta)+".txt"
    with open(fileName, "ab") as f:
        for l in L:
            a = l.numpy()
            numpy.savetxt(f, a, fmt = '%i', delimiter=" ")
            f.write(b"\n")
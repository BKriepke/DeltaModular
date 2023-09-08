import itertools

# removes equivalent HNFs by bruteforce
def remove(L):
    r = L[0].ncols()
    for x in L: x.set_immutable()
    L = set(L)
    # need to loop through a copy of the set, as we are changing it during the loop
    Lc = set(L)
    p = Permutations(r)
    d = [1, -1]
    # compute all possible options C=D*P with D=diag(+-1) and P permutation matrix
    CList = []
    for diag in itertools.product(d, repeat=r):
        # first entry of D can be positive
        # because ADP = A(-I)(-D)P = (-I)A(-D)P is equivalent, with -I unimodular
        if diag[0] == -1: continue
        D = diagonal_matrix(diag)
        for pi in p:
            P = pi.to_matrix()
            CList.append(D*P)
    for A in Lc:
        if A in L:
            # calculate all HNFs that are equivalent to A and remove them
            for C in CList:
                B = A*C
                H = B.hermite_form()
                if H!= A:
                    try:
                        L.remove(H)
                    except:
                        pass
                        
    return list(L)

# Computes all matrices in hermite normal form of size 2x2
# with determinant Delta
def allHNF2(Delta):
    divs = divisors(Delta)
    diags = []
    for div in divs:
        if div*div <= Delta:  # so diagonal is sorted
            diags.append([div, Delta/div])
            
    H = []
            
    for d in diags:
        for i in range(0, d[1]//2+1):
            A = matrix(ZZ, [[d[0], i],
                            [0, d[1]]])
            H.append(A)
                    
    H = remove(H)
    return H

def findKernelModDelta(A, Delta):
    # finds the solutions to Ax=0 mod Delta in Z_Delta^r
    # https://ask.sagemath.org/question/33890/how-to-find-kernel-of-a-matrix-in-mathbbzn/
    r = A.nrows()
    Zr = ZZ^r
    M = (Zr)/(Delta*Zr)
    # Sage calculates left kernel, so we have to transpose
    phi = M.hom([M(a) for a in A.transpose()])
    V = [M(b) for b in phi.kernel()]
    cols = []
    # calculate all integer representants, first entry positive
    for v in V:
        I = []
        if v[0] == 0: I.append([Delta])
        else: I.append([v[0]])
        for i in range(1,r):
            if v[i] == 0: I.append([-Delta, Delta])
            else: I.append([v[i], v[i]-Delta])
        for element in itertools.product(*I):
            cols.append(matrix(ZZ, r, 1, element))
    cols = removeMultiples(cols)
    return cols

def removeMultiples(cols):
    colsC = cols[:]
    for v in colsC:
        for mu in range(2,Delta+1):
            try:
                cols.remove(mu*v)
            except:
                pass
    return cols

def buildC(cols, Delta, single):
    # Builds the graph G_2
    # if single==True, then finds a single maximum clique in that graph
    # otherwise find all of them
    # return the corresponding totally generic Delta-bound matrices C
    
    D = cols[0]
    for v in cols[1:]: D = D.augment(v)
    n = D.ncols()
    
    # initialize 2-dim array to be used for adjacency matrix
    # M[i][j] = 1 iff (v_i, v_j) Delta-bound and totally generic
    M = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        v = cols[i]
        for j in range(i+1,n):
            B = v.augment(cols[j])
            if abs(B.det()) > Delta^2 or B.det() == 0:
                continue
            M[i][j] = 1
            M[j][i] = 1

    # Adjacency matrix for G_2 from M
    B = matrix(ZZ, n, n, M)
    G = Graph(B)
    
    if single:
        C = sage.graphs.cliquer.max_clique(G)
        return D.matrix_from_columns(C)
    else:
        Cs = sage.graphs.cliquer.all_max_clique(G)
        return [D.matrix_from_columns(K) for K in Cs]

def matrixInterval(k, ak, bk):
    # Computes matrix
    # [k   ... k]
    # [a_k ... b_k]
    # where the dots represent ever integer a_k <= i <= b_k coprime to k
    I = [i for i in range(ak, bk+1) if gcd(k, i) == 1]
    A = matrix(ZZ, [len(I)*[k],
                    I])
    return A

def M(a,b):
    # takes lists a, b
    # computes matrix M(a,b) introduced before Prop 3.2
    D = matrix(ZZ, [[0],
                    [1]])
    for i in range(len(a)):
        D = D.augment(matrixInterval(i+1, a[i], b[i]))
    return D

def families1to3(Delta):
    # checks if Delta belongs to families (F1)-(F3) and if yes, return a corresponding matrix
    # otherwise returns None
    if Delta+1 == Delta.next_prime():
        return M([0], [Delta])
    if Delta+2 == Delta.next_prime():
        return M([0, Delta], [Delta, Delta])
    if Delta % 3 == 2 and Delta+3 == Delta.next_prime():
        s = Delta//12
        if Delta % 12 == 2:
            return M([0, 4*s+1, 9*s+1], [7*s+1, 10*s+1, 12*s+2])
        if Delta % 12 == 8:
            return M([0, 4*s+3, 9*s+7], [7*s+5, 10*s+7, 12*s+8])
           
def sortMatrix(A):
    # scales all columns so that first entry is positive and then sorts lexicographically by rows
    I = identity_matrix(A.ncols())
    for i in range(A.ncols()):
        if A[0][i] < 0: I[i, i] = -1
        if A[0][i] == 0:
            if A[1][i] < 0: I[i, i] = -1
    A = A*I
    return matrix(ZZ, A.ncols(), 2, sorted(A.columns())).transpose()

def invariantVector(A, Delta):
    # computes list v
    # with v_i = number of minors of A that have absolute value i, i=1, .., Delta
    # vector v is invariant under equivalence relation defined in Definition 4.5
    n = A.ncols()
    r = A.nrows()
    I = list(range(n))
    S = Subsets(I, r)
    v = [0 for i in range(Delta)]
    for s in S:
        J = list(s)
        B = A.matrix_from_columns(J)
        d = abs(B.det())
        v[d-1] = v[d-1]+1
    return v

def maxDigits(A):
    # returns max digits needed to write every entry in matrix A
    m = max(max(A))
    return len(str(m))

def fetchTargetSize(Delta):
    # returns g(Delta, 2)-2, as this is the size of the maximum cliques that are relevant
    # assumes that this value exists in the file
    with open("../data/Generic/r=2/values.txt", 'r') as f:
        for line in f:
            [a, b] = [int(x) for x in line.split()]
            if a == Delta: return b-2

def singleExample(Delta, upperBound):
    # finds g(Delta, 2) and returns a single example matrix D achieving this value
    # Cf algorithm 1 in paper
    r = 2
    single = True
    H = allHNF2(Delta)
    print("Number of HNFs: {:}".format(len(H)))
    maxSoFar = 0
    counter = 0
    for A in H:
        if maxSoFar == upperBound: break
        cols = findKernelModDelta(A, Delta)
        C = buildC(cols, Delta, single)
        if C.ncols()+r > maxSoFar:
            maxSoFar = C.ncols()+r
            D = A.augment(A*C/Delta)
            print("Current max length found: {:} ".format(maxSoFar))
            print("Current max length matrix:  ")
            print(D)
            
        counter +=1
        if counter % 10 == 0: print("HNFs checked: {:} / {:} ".format(counter, len(H)))
    
    return D

def allExamples(Delta):
    # reads g(Delta, 2) from file and then returns a list of matrices D which achieve this value
    # this list is guaranteed to contain at least one representative of every equivalence class
    # but will in general contain more than a single representative per class
    r = 2
    single = False
    H = allHNF2(Delta)
    print("Number of HNFs: {:}".format(len(H)))
    Ds = []
    targetSize = fetchTargetSize(Delta)
    counter = 0
    for A in H:
        cols = findKernelModDelta(A, Delta)
        Cs = buildC(cols, Delta, single)
        if Cs[0].ncols() < targetSize: continue
        for C in Cs:
            D = A.augment(A*C/Delta)
            Ds.append(D)
            
        counter +=1
        if counter % 10 == 0: print("HNFs checked: {:} / {:} ".format(counter, len(H)))
    
    return Ds
import numpy

load('helpfunctions.sage')

# assumes g(Delta, 2) is already calculated for each startDelta <= Delta <= endDelta
# returns a list of matrices D which attain this value
# this list is guaranteed to contain at least one representative of every equivalence class
# but will in general contain more than a single representative per class
# edit to change the values the script calculates
startDelta = 2
endDelta = 400
# numberFile = open('numberstest.txt', 'a+')
for Delta in range(startDelta, endDelta+1):
    Delta = Integer(Delta)
    print("Current Delta: {:} ".format(Delta))
    
    Ds = allExamples(Delta)
    
    print("Number of matrices found: {:}".format(len(Ds)))

    # For each matrix computes the invariant vector, 
    # to find a lower bound on number of equivalence classes

    # Vs = set()
    # for D in Ds:
    #     v = invariantVector(D, Delta)
    #     Vs.add(tuple(v))
    # print("Number of different invariant vectors found: {:}".format(len(Vs)))
    # print("Invariant vectors: ")
    # for v in Vs:
    #     print(v)
    
    specificFile = open("Delta="+str(Delta)+"all.txt", 'w+')
    print("Matrices: ")
    for D in Ds:
        D = sortMatrix(D)
        d = D.numpy()
        digits = maxDigits(D)
        numpy.savetxt(specificFile, d, fmt = '%'+str(digits)+'i', delimiter=", ")
        specificFile.write("\n")
        print(D)
        
    specificFile.close()
    # numberFile.write("{:}\t{:}\t{:}\n".format(Delta, len(Ds), len(Vs)))
    # numberFile.flush()
    
# numberFile.close()
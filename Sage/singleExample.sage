import numpy

load('helpfunctions.sage')

# calculates g(Delta, 2) and returns one matrix D achieving this value
# for each startDelta <= Delta <= endDelta
# edit to change the values the script calculates
startDelta = 401
endDelta = 450
valueFile = open('../data/Generic/r=2/values.txt', 'a+')
for Delta in range(startDelta, endDelta+1):
    Delta = Integer(Delta)
    print("Current Delta: {:} ".format(Delta))
    D = families1to3(Delta)
    if D is None:
        D = singleExample(Delta, Delta.next_prime()+1)
        D = sortMatrix(D)
        print("Max length found: {:}".format(D.ncols()))
    else: 
        print("Is one of families (F1)-(F3)")
        print("Therefore max length: {:}".format(D.ncols()))
    print("One matrix with that length: ")
    print(D)
    
    specificFile = open("../data/Generic/r=2/SingleExample/Delta="+str(Delta)+".txt", 'w+')
    d = D.numpy()
    digits = maxDigits(D)
    numpy.savetxt(specificFile, d, fmt = '%'+str(digits)+'i', delimiter=", ")
    specificFile.close()
    
    valueFile.write("{:}\t{:}\n".format(Delta, D.ncols()))
    valueFile.flush()
valueFile.close()
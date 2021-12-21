from math import sqrt, log
from  typing import List, Tuple, Dict
import time

def P3inputCheck(P2Out_: Dict[Tuple, Tuple]) -> bool:
    """
    Checks the output of P2 and returns a True if malformed input and false otherwise
    """
    malInp = False
    if not isinstance(P2Out_, Dict):
        malInp = True
    # check if each element fits the type criteria
    for (x, y),(seqx, seqy) in P2Out_.items():
        if not isinstance(x, int) or not isinstance(y, int) or not isinstance(seqx, str) or not isinstance(seqy, str):
            malInp = True
    
    return malInp

def ComputeDistMatrix(P2Out: Dict[Tuple[int, int], Tuple[str, str]]) -> List[List[float]]:
    """ 
    This function reads in the Output of P2 witht the type Dict(Tuple : Tuple) and then:
    1. Creates a empty matrix (List of Lists) with the dimensions: number of Sequences x number of Sequences
    2. Compares each Element of the aligned sequences and calculates the distance and the adjusted distance value 
    """
    start = time.time()

    #check input
    if P3inputCheck(P2Out):
        print("P3: Malformed Input")
        return

    #calculate the number of sequences by the number of pairs (== lenght of the dict) with 
    # the inverse equation of: 
    # nrPairs = nrElements * (nrElements -1) / 2    -->    nrElements = sqrt(2 * nrPairs + 1/4) + 1/4
    numOfSeq = int(sqrt(2 * len(P2Out) + 1/4) + 1/2)

    #Init distance matrix
    distMa = [[0.0 for c in range(numOfSeq)] for r in range(numOfSeq)]


    #Calculate distance matrix
    for (x, y),(seqx, seqy) in P2Out.items():
        seqx = list(seqx)
        seqy = list(seqy)
        nrNukDiff = 0
        nrNukComp = 0
        # compare each Nukleotide
        for count, (xElement, yElement) in enumerate(zip(seqx, seqy)):

            #pass iteration if one of the sequences == "-"
            if xElement == "-" or yElement == "-":
                pass
            #if same
            elif xElement == yElement:
                nrNukComp += 1
            #add 1 to counter if elements are not equal
            elif xElement != yElement:
                nrNukDiff += 1
                nrNukComp += 1

        #Calculate the corrected evulutionary distance and dont make a deep copy, because matrix is symmetrical
        try: #check if the value inside the log is not negative
            distMa[y][x] = distMa[x][y] = -3/4 * log( 1 - 4/3 * (nrNukDiff / nrNukComp) )
        except RuntimeError:
            print("P3: Malformed Input")

    end = time.time()
    print(f"P3 time: {end-start} sec")
    return distMa

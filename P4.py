from collections import defaultdict
from  typing import List, Tuple, Dict
import time 
import re

def P4inputCheck(P3Out_: List[List], label_: List) -> bool:
    """
    Check if the output of P3 matches the type criteria 
    """
    malInp = False

    if not isinstance(label_, list) or not isinstance(P3Out_, List):
        malInp = True

    #assuming if first element is float rest is too
    if not isinstance(P3Out_[0][1], float):
        malInp = True

    for label in label_:
        if not isinstance(label, str):
            malInp = True
    
    return malInp

def connectP1AndP3(label: List, P3Out: List[List]) -> Tuple[Dict[Tuple, float], int]:
    """ 
    read in the output of P1 and P3 and connect the index of the label 
    with the values in the distance Matrix to create a dictionair with "labelIndex" : "Distance" 
    """
    #nr of pairs calculated by nr of sequences
    nPair = int(len(label) * (len(label) -1) / 2)

    dici = defaultdict(dict)

    #create all possible label Indeces with the same order as the real data
    labelLength = len(label)
    for line, lab in enumerate(label):
        for c in range(1,labelLength):
            #create 2d dictionair with labels, nr of labels and distMatrix values
            dici[nPair][(int(line), int(line + c))] = P3Out[line][line+c]
        labelLength -= 1
    return dici, nPair
 
def index2Label(indexTree_: Tuple, labelLength_: int, label_: List) -> str:
    """
     Transform the final clustered index tree to label tree by filtering all the elements in the index tree, 
     create the corresponding label list, and replace each occurancy of a element in the indextree with the label from the labellist
     """

    #create dict to replace string index to label
    labelDict = dict(zip( [str(i) for i in range(labelLength_)], label_))

    pat = r'[() ]'
    indexList = re.sub(pat, "", str(indexTree_)).split(",")
    pat2 = r'\d+'
    treeStructure = re.sub(pat2, "x", str(indexTree_)).replace(" ", "")

    labelList = []
    for ele in indexList:
        labelList.append(labelDict[ele])
    
    for lab in labelList:
        treeStructure = treeStructure.replace("x", lab, 1)

    return treeStructure

def caseSwitcher(first_: int or Tuple, second_: int or Tuple , keyUpdate_: int or Tuple) -> Tuple:
    """ 
    This function is needed to give a logic to the building of the label tree. 
    choose the variables a, cA, b and cB for the equation ((a,b), c) = (D(a,c) + D(b,c)) / 2 to then
    get the corresponding values in the dictionair "dict". The smaller hash value of an element gets the left 
    posistion (a or b)
    """

    if hash(first_) < hash(keyUpdate_):
        a, cA   = first_, keyUpdate_
    else:
        a, cA   = keyUpdate_, first_

    if hash(second_) < hash(keyUpdate_):
        b, cB   = second_, keyUpdate_
    else:
        b, cB   = keyUpdate_, second_

    return a, cA, b, cB

def Cluster(label: List[str], P3Out: List[List[float]]) -> str:
    """ 
    The Goal of this function is to build a label Tree from the P1 labels and the P3 Output with the help
    of the dictionair keys of the P2 Output
    1. Find the smallest Distance of the distance Matrix Dictionair.
    2. Add the labelIndeces which are not used in the current iteration to a set (to remove repetition)
    3. Split the keys of the current minimum, order them with caseSwitcher(),  calculate the new distances
       and add these to a new subDict()
    4. if a 3x3 matrix is left, do the final join and transform the indexLabel tree with the index2Label() function
       into a label Tree. 
    """
    start = time.time()

    #check correct type
    if P4inputCheck(P3Out, label):
        print("P4: Malformed Input")
        return

    #create dict
    dici, nPair = connectP1AndP3(label, P3Out)
    labelLength = len(label)

    exi = 0
    while exi == 0:
        #a 'for loop' inside a 'while loop' is needed so we can iterate over the new dict as soon as the new Matrix is calculated
        for count, dictValue in enumerate(dici[nPair].values()):
        #-------------- if 3x3 matrix left, exit with result -------------------------#
            if len(dici[nPair].values()) <= 3 and list(dici[nPair].values())[count] == min(dici[nPair].values()):

                #instead of a set, we just use one of the other two remaining possibilities and get the keys
                first, second = list(dici[nPair])[count-1]
                #
                if first not in list(dici[nPair])[count]:
                    indexTree = (first, list(dici[nPair])[count]) 
                    labelTree = index2Label(indexTree, labelLength, label)
                    
                else:
                    indexTree = (second, list(dici[nPair])[count]) 
                    labelTree = index2Label(indexTree, labelLength, label)
                
                end = time.time()
                print(f"P4 time: {end-start} sec")
                return labelTree

        #-------------- if not, finalize matrix reduction -------------------------#
            else:
                smallest = min(dici[nPair].values())

                if dictValue != smallest:
                    pass
                else:
                    break

        # to evaluate which sequences are not combined after finding the minimum Distance we insert all these 
        # key elements to a set
        filterSet = set()
        #for python >3.6  dicts are ordered and indexes (list(dict)) can be used without losing the order
        first, second = list(dici[nPair])[count]
        for key in dici[nPair].keys():
            if first not in key and second not in key:
                #add the pairs into the new dict which are not used for the calculation
                dici[nPair-1][key] = dici[nPair][key]
                #shorten these pairs into a set and use the numbers to create the new combinations
                firstElement, secondElement = key
                filterSet.add(firstElement)
                filterSet.add(secondElement)
    #-------------- Calculate new distance matrix -------------------------#
        for keyUpdate in filterSet:
            #use caseSwitcher() to define how to extract the values from dict to calculate new distances
            a, cA, b, cB = caseSwitcher(first, second, keyUpdate)

            # Formula for the new Distances: D((a,b), c) = (D(a,c) + D(b,c)) / 2
            Distance = (dici[nPair][(a, cA)] + dici[nPair][(b, cB)]) / 2

            #As a rule, save the new key with corresponding distnance depending on how many elements the tuple includes, the smaller one goes left
            # if they have the same number of elements, compare the sum of these elements, the smaller tuple goes left
            
            if hash((first, second)) < hash(keyUpdate):
                dici[nPair-1][((first, second), keyUpdate)] = Distance
            else:
                dici[nPair-1][(keyUpdate, (first, second))] = Distance
      
                
        #Iterate over next dict
        nPair -= 1
        count = 0
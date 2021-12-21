from typing import List, Tuple

#typing hint of python >3.8
def lineSpliter(path_: str) -> List:
    """ 
    This function opens the file by the given path, iterates through each line.
    If the line has a ">" sign, strip the strings infront of the label and split the lines by empyt characters.
    The first element in the list will be label if the input type is correct.  
    """
    fasta_ = open(path_, "r")
    fasta = fasta_.read()
    
    li = []
    for i in fasta.splitlines():
        if ">" in i:
            try:
                i = i.strip("> ")
                i = i.replace("\t", "")
                i = i.split(sep=" ")
                li.append(i)
            except:
                print("Malformed Input: Line19")
        else:
            pass
        
    fasta_.close()
    return li

def ParseSeqFile(path: str) -> List[Tuple[str, str]]:
    """ 
    This function reads in the file with the given path as string and does the following:
    1. Check if the input is type string
    2. Filters the labels and the sequences with the lineSpliter() function. 
    3. Checks the sequence Nukleotides other than A, C, G or T
    4. Returns the desired list ot tuples
    """
    #input type controle
    if not isinstance(path, str):
        return print("Malformed input")

    li = lineSpliter(path)
    li2 = []
    for i in li:
        lab = i[0]
        seq = "".join(i[1:])
        #check if all Nucleotides are in the sequences
        for i in "ACGT":
            if i in seq:
                continue
            else:
                print(f"{i} is missing in the sequence of {lab} ")
        #check if there are other Characters
        if seq.strip("ACGT") != "":
            oh = seq.strip("ACGT")
            print("malformed input")
        li2.append(
            (lab, seq)
        )

    return li2
    

def getLabel(path: str) -> List[str]:
    """
    This function extracts the labels from the P1 input for the usage of P4.
    """
    li = lineSpliter(path)
    label = []
    for i in li:
        label.append(i[0])
    
    return label

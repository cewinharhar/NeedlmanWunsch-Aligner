
import os
from time import time
ts = time()
path_ = r"C:\Users\kevin\OneDrive - ZHAW\KEVIN STUFF\ZHAW\MASTER\V5_1_Programming Algorithms and Data-Structures\finalProject\kevin_yar"
os.chdir(path_)
import P1, P2, P3, P4
#P1 
SequenceList = P1.ParseSeqFile('enzymeAlignment.txt')
print("#############")
#P2 
LabelSequenceDict = P2.AlingByDP(SequenceList)     
print("#############")
#P3 
DistMatrix = P3.ComputeDistMatrix(LabelSequenceDict)  
#P4
label = P1.getLabel("enzymeAlignment.txt")
BinaryTreeString = P4.Cluster(label, DistMatrix)
print(BinaryTreeString)
print("#############")
print(f"{time() -ts}")






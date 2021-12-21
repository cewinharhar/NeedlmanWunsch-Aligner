from copy import deepcopy
from  typing import List, Tuple, Dict
import time

def P2inputCheck(P1Out_: List[Tuple[str, str]]) -> bool:
  """
  Checks if the output of P1 is a List of tuples 
  """
  #check input
  malInp = False
  if not isinstance(P1Out_, list):
      malInp = True
  for (fi, se) in P1Out_:
    if not isinstance(fi, str) or not isinstance(se, str):
      malInp = True
  
  return malInp

def Aligner(sx_: str, sy_: str) -> Tuple:
  """
  The core algorithm based on the needlman wunsch algorithm for the alignment. 
  This function reads in the two current sequences as strings and returns the alignment as tuple by:
  1. Creating the initializer matrix (list of lists) for the algorithm and a shadow matrix with the same dimensions for the backtracing
  2. Run the Needlman-wunsch algorithm
  3. Run the backtracing
  """
   #-------------  declare variables --------------------#

  sx  = " " + sx_
  m   = len(sx) 
  sy  = " " + sy_
  n   = len(sy)

  MA      = 5
  MIMA    = -2
  INDEL   = -6

  #-------------  Needlman-wunsch algorithm --------------------#

  #create a m x n list-matrix and the corresponding shadow matrix to store the decisions
  S   = [[0 for c in range(n)] for r in range(m)]
  spy = deepcopy(S)

  #Fill the first column and the first row 
  for indx in range(m): 
    S[indx][0] = indx * INDEL
  for indy in range(n): 
    S[0][indy] = indy * INDEL

  for i in range(1, m):
    for j in range(1, n):
  #----------------------#
      #is match or not
      if sx[i] == sy[j]:
        ad = MA
      else: 
        ad = MIMA

      #check the values
      left  = S[i][j-1] + INDEL
      corn  = S[i-1][j-1] + ad
      above = S[i-1][j] + INDEL

      #compare which value is bigger, start with corner
      if corn >= left and corn >= above:
        S[i][j] = corn
        spy[i][j]= 1

      elif left >= above and left >= corn:
        S[i][j] = left 
        spy[i][j]= 0
    
      elif above >= left and above >= corn:
        S[i][j] = above 
        spy[i][j]= 2
  #------------- Back-tracing -----------------#

  sxA = []
  syA = []
  m_, n_ = m -1 , n -1 
  esc = 0
  while esc == 0:
    if n_ == 0 and m_ == 0:
      esc = 1
    #Security measur in case of algorithm failing
    elif n_ == -1 or m_ == -1:
      esc = 1
      #start with corner
    else:
      if spy[m_][n_] == 1:
          sxA.insert(0, sx[m_])
          syA.insert(0, sy[n_])
          n_ -= 1
          m_ -= 1
      elif spy[m_][n_] == 0:
          sxA.insert(0, "-")
          syA.insert(0, sy[n_])
          n_ -= 1

      elif spy[m_][n_] == 2:
          sxA.insert(0, sx[m_])
          syA.insert(0, "-")
          m_ -= 1
  # the replace command is just for security purposes when n_ or m_ gets to -1
#  return ("".join(sxA).replace(" ", "-"), "".join(syA).replace(" ", "-"))
  return ("".join(sxA), "".join(syA))
  

def AlingByDP(P1Out: List[Tuple[str, str]]) -> Dict[Tuple[int, int], Tuple[str, str]]:
  """ 
  This function reads in the output of P1 and:
  1. sets the dict keys of the P2 Output
  2. Aligns the sequences with the Aligner() function
  3. returns the final dictionair
  """
  start = time.time()

  #check input
  if P2inputCheck(P1Out):
    print("P2: Malformed Input")
    return

  intPair = []
  seq     = []
  l = len(P1Out)
  for line, (lab, seq_) in enumerate(P1Out):
    # The numbers start with 0 due to no exact clarification
    for c in range(1,l):
      intPair.append(
        (int(line), int(line + c))
      )
    l -= 1
    seq.append(seq_)

  #Align sequences and return the dict
  seqPair = []
  for (first, second) in intPair:
    #progress display
    print(f"AlignByDP: {(first, second)} / {intPair[-1]}", end='\r')
    seqPair.append(
      Aligner(seq[int(first)], seq[int(second)])
    )
  end = time.time()
  print(f"P2 time: {end-start} sec")
  return dict(zip(intPair, seqPair))


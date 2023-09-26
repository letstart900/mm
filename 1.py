import math

D = 100
K = 100
h = 0.02
L = 12

y_star = (2*K*D/h)**0.5
t0 = y_star/D #cycle length

if(L>t0):
  n = math.floor(L/t0)
  le = L - n*t0
  print("Reorder point thus occurs when inventory model drops to : ",le*D)
else:
  print("Reorder point thus occurs when inventory model drops to : ",L*D)

TCU_y = K*D/y_star + h*y_star/2
print("Inventory cost associated with the proposed policy is : ",TCU_y)

"""##Price Break

"""

from sympy import *

D = 187.5
h = 0.02
K = 20
L = 2
C1 = 3
C2 = 2.5
q = 1000

ym = (2*K*D/h)**0.5
TCU1_y = C1*D + K*D/ym + h*ym/2

def fun_Q():
  var('Q')
  return max(solve(Q**2 + (2*(C2*D - TCU1_y))*Q/h + 2*K*D/h,Q))

Q = fun_Q()

y_star = (2*K*D/h)**0.5
t0 = y_star/D #cycle length

if(L>t0):
  n = math.floor(L/t0)
  le = L - n*t0
  print("Reorder point thus occurs when inventory model drops to : ",le*D)
else:
  print("Reorder point thus occurs when inventory model drops to : ",L*D)

if(q<ym or q>Q):
  print("y* = ",ym)
else:
  print('y* = ',q)

"""## Multi item Storage limitation"""

import numpy as np
from scipy.optimize import minimize

K = [10,5,15]
D = [2,4,4]
h = [0.3,0.1,0.2]
a = [1,1,1]
A = 25

y_star = []
for i in range(len(K)):
  y_star.append((2*K[i]*D[i]/h[i])**0.5)

if(sum(y_star)>A):
  print("Violates!")

def objective(Y):
  val = 0
  for i in range(len(Y)):
    val+=(K[i]*D[i]/Y[i] + h[i]*Y[i]/2)
  return val

def constraint(Y):
  val = 0
  for i in range(len(Y)):
    val -= a[i]*Y[i]
  return val + A

Y = [2 for i in range(len(K))]
constr = {
    'type' : 'ineq',
    'fun' : constraint
}
bounds = ((0,float('inf')) for i in range(len(K)))
cons = [constr]
sol = minimize(objective,Y,method='SLSQP',constraints=cons)
print(sol)

"""## Needleman

"""

import numpy as np

def needlemanWunsch(seq1: str, seq2: str, matchScore: int, misMatch: int, gapScore: int):

  l1, l2 = len(seq1), len(seq2)

  # Empty array => l1 * l2
  arr = [[0 for i in range(l2+1)] for j in range(l1+1)]
  # First row
  for j in range(l2+1):
    arr[0][j] = gapScore * (j)
  # First col
  for i in range(l1+1):
    arr[i][0] = gapScore * (i)

  # Matrix construction
  for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
      arr[i][j] = max([
          arr[i-1][j] + gapScore,
          arr[i][j-1] + gapScore,
          arr[i-1][j-1] + (matchScore if seq1[i-1] == seq2[j-1] else misMatch)
        ])

  print('Matrix: \n', np.matrix(arr))

  # Traceback
  path, trace = [], []
  modSeq1, modSeq2 = '', ''
  i, j = l1, l2

  while i>0 and j>0:
    path.append(arr[i][j])

    if seq1[i-1] == seq2[j-1]:
      # match
      i -= 1
      j -= 1
      modSeq1 += seq1[i]
      modSeq2 += seq2[j]
      trace.append('M')
    else:
      # choose max neighbour
      left = arr[i][j-1]
      top = arr[i-1][j]
      diag = arr[i-1][j-1]

      if left > top and left > diag:
        j -= 1
        modSeq1 += '-'
        modSeq2 += seq2[j]
        trace.append('<')

      elif top > left and top > diag:
        i -= 1
        modSeq1 += seq1[i]
        modSeq2 += '-'
        trace.append('^')

      else:
          i -= 1
          j -= 1
          modSeq1 += seq1[i]
          modSeq2 += seq2[j]
          trace.append('M')

  while i>0:
    i -= 1
    trace.append('^')
    modSeq1.append(seq1[i])
    path.append(arr[i][0])

  while j>0:
    j -= 1
    trace.append('<')
    modSeq2.append(seq2[j])
    path.append(arr[0][j])

  print(f'Retrace scores: {path[::-1]}')
  print(f'Retrace chages: {trace[::-1]}\n')

  return modSeq1[::-1], modSeq2[::-1]

if __name__ == '__main__':
  str1 = 'AGCT'
  str2 = 'ATGCT'

  matchScore = 1
  misMatchScore = -1
  gapPenalty = -2

  s1, s2 = needlemanWunsch(str1, str2, matchScore, misMatchScore, gapPenalty)
  print(f'Needleman Wunsch: {s1} & {s2}')

"""## Smith"""

def maxIndex(list2d):
  flatList = [el for row in list2d for el in row]
  print(flatList)
  maxIndex = flatList.index(max(flatList))

  print(maxIndex)

  rowNo = int(maxIndex/len(list2d[0]))
  colNo = maxIndex%len(list2d[0])

  return rowNo, colNo

def smithWaterman(seq1: str, seq2: str, matchScore: int, misMatch: int, gapScore: int):

  l1, l2 = len(seq1), len(seq2)

  arr = [[0 for i in range(l2+1)] for j in range(l1+1)]
  for j in range(l2+1):
    arr[0][j] = gapScore * (j)
  for i in range(l1+1):
    arr[i][0] = gapScore * (i)

  for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
      arr[i][j] = max([
          arr[i-1][j] + gapScore,
          arr[i][j-1] + gapScore,
          arr[i-1][j-1] + (matchScore if seq1[i-1] == seq2[j-1] else misMatch),
          0
        ])

  print('Matrix: \n', np.matrix(arr))

  path, trace = [], []
  modSeq1, modSeq2 = '', ''

  i, j = maxIndex(arr)

  while arr[i][j] != 0:
    path.append(arr[i][j])

    if seq1[i-1] == seq2[j-1]:
      # match
      i -= 1
      j -= 1
      modSeq1 += seq1[i]
      modSeq2 += seq2[j]
      trace.append('M')
    else:
      # choose max neighbour
      left = arr[i][j-1]
      top = arr[i-1][j]
      diag = arr[i-1][j-1]

      if left > top and left > diag:
        j -= 1
        modSeq1 += '-'
        modSeq2 += seq2[j]
        trace.append('<')

      elif top > left and top > diag:
        i -= 1
        modSeq1 += seq1[i]
        modSeq2 += '-'
        trace.append('^')

      else:
          i -= 1
          j -= 1
          modSeq1 += seq1[i]
          modSeq2 += seq2[j]
          trace.append('M')

  print(f'Retrace scores: {path[::-1]}')

  return modSeq1[::-1], modSeq2[::-1]

if __name__ == '__main__':
  str1 = 'AGCT'
  str2 = 'ATGCT'

  matchScore = 1
  misMatchScore = -1
  gapPenalty = -2

  s1, s2 = smithWaterman(str1, str2, matchScore, misMatchScore, gapPenalty)
  print(f'\nSmith Waterman: {s1} & {s2}')

"""##Dynamic - With setup"""

inputDict = {
    "periods" : [1, 2, 3],
    "demand" : [3, 2, 4],
    "setup" : [3, 7, 6],
    "holding": [1, 3, 2]
}

def findCost(z, s):
  if z == 0:
    return 0
  if z <= 3:
    return (10*z) + s
  else:
    return (30 + 20*(z-3)) + s

def baseInventory():
  demand = inputDict['demand'][0]
  xList = [1]

  maxRange = sum(inputDict['demand']) - demand
  diffOrder = demand - xList[0]
  holdingCost = inputDict['holding'][0]
  setupCost = inputDict['setup'][0]
  baseMinInventCost = []

  for x in range(0,maxRange + 1):
    hx = holdingCost * x
    z = x + diffOrder
    CZ = findCost(z, setupCost)
    fx = CZ + hx
    baseMinInventCost.append(fx)

  return baseMinInventCost

def helper(period,MinInvestCost):
  demand = inputDict['demand'][period - 1]

  maxRange = 0
  for i in range(period, len(inputDict['demand'])):
    maxRange += inputDict['demand'][i]

  diffZ = demand

  zMaxRange = maxRange + diffZ
  setupCost = inputDict['setup'][period - 1]
  holdingCost = inputDict['holding'][period -1]

  finalFxList = []

  for x in range(0, maxRange + 1):
    fxList = []
    for z in range(0, x + diffZ + 1):
      hx = holdingCost * x
      CZ = findCost(z, setupCost)
      fxList.append(CZ + hx + MinInvestCost[x - z + demand])
    finalFxList.append(min(fxList))

  return finalFxList

MinInvestCost = baseInventory()
print("The optimum costs at period 1 is f1(x2)",MinInvestCost)
for p in range(2, len(inputDict['periods']) + 1):
  MinInvestCost = helper(p,MinInvestCost)
  print("The optimum costs at period",p,"is f( x",p+1,")",MinInvestCost)

print("The optimal cost is ",MinInvestCost[0])

"""## Dynamic - without setup

"""

#input dictionary
inputDict = {"Month":[1,2,3,4], "Regular":[90,100,120,110], "Overtime":[50,60,80,70], "Demand":[100,190,210,160]}

#create 2d matrix
rows = len(inputDict["Regular"])*2
cols = len(inputDict["Demand"])

#initialise the variables
mat = [[0 for i in range(cols)] for j in range(rows)]
regular = inputDict["Regular"]
overtime = inputDict["Overtime"]
demand = inputDict["Demand"]

#fill regular demand
for inx in range(0,len(regular)):
  output = min(regular[inx], demand[inx])
  demand[inx] -= output
  regular[inx] -= output
  mat[inx*2][inx] = output

  if(demand[inx] > 0):
    for jnx in reversed(range(inx)):
      output = min(regular[jnx], demand[inx])
      if(output > 0):
        demand[inx]-= output
        regular[jnx] -= output
        mat[(jnx*2)][inx] = output

#fill overtime demand
for inx in range(0, len(overtime)):
  output = min(overtime[inx], demand[inx]);
  demand[inx] -= output
  overtime[inx] -= output
  mat[(inx*2) + 1][inx] = output

  if(demand[inx] > 0):
    for jnx in reversed(range(inx)):
      output = min(overtime[jnx], demand[inx])
      if(output > 0):
        demand[inx]-= output
        overtime[jnx] -= output
        mat[(jnx*2)+1][inx] = output

#find surplus
surplus = 0
for inx in range(len(regular)):
  surplus += regular[inx]
  surplus += overtime[inx]

print("The surplus:",surplus)
print("\nMatrix:\n", mat)

price = 0

#find the total cost
for i in range(rows):
  idx = 0
  for j in range(i//2, cols):
    #print(i, j, mat[i][j])
    if i % 2 == 0:
      price += mat[i][j] * (6 + (idx * 0.1))
    else:
      price += mat[i][j] * (9 + (idx * 0.1))
    idx += 1
    #print(i, j, mat[i][j],price)
print(price)

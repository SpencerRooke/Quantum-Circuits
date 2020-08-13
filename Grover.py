#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spencer
"""
import numpy as np
import qClass
from qGates import *

def GDiffusionGate(N): #Returns Grover's Diffusion Gate on N qubits
    HN = Hadamard.matrix
    for i in range(1,N):
        HN = np.kron(HN,Hadamard.matrix)
    
    Mat = -1*nIdentity(N)
    Mat[0][0] = 1
    
    outArr = np.matmul(HN,Mat)
    outArr = np.matmul(outArr,HN)
    outLetter = 'G'+str(N) #might break if N>99
    
    return qClass.genericOperator(outArr,outLetter)

def Oracle(N,n): #Returns an Oracle for N qubits, with indicator at n (eg N,n=4, -> flips the sign of |0100> )
    Mat = nIdentity(N)
    Mat[n][n] = -1
    return qClass.genericOperator(Mat,'U_ω')

#------------------------------------------------------------------------------
print('Input number of qubits N:')
N = int(input())

print('Oracle indicator Number n:')
oracleNum = int(input())
  
newOracle = Oracle(N,oracleNum)
gGate = GDiffusionGate(N)
zeroqBit = qClass.qubit(1,0) #|0> state
Register = qClass.createRegister(np.full(N,zeroqBit))

circuitMat =  [np.full(N,None),np.full(N,Hadamard)]

for i in range(int(np.ceil(np.sqrt(N)))):
    circuitMat = circuitMat + [[newOracle]] + [[gGate]] 
circuitMat = circuitMat + [np.full(N,None)]

Grover = qClass.qCircuit(circuitMat,N)
Grover.draw()
#------------------------------------------------------------------------------
oracleState = bin(oracleNum)[2:]
while(len(oracleState)<N):oracleState = '0'+oracleState
oracleState = '|'+oracleState+'>'

count = 0
for i in range(10000):
    A = Grover.Simulate(Register)
    meas = qClass.Measurement(A)
    if(meas == oracleState): count = count+1

print("In 10,000 circuit simulations, the state "+oracleState+" was measured "+str(count)+" times ("+str(count/100)+"%)")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spencer
"""
import numpy as np
import qClass
from qGates import *
from collections import Counter

N = 15 #Number this will Factor

#All the controlled Nots we'll need
tNOT = controlled(PauliX,3)
NOTt = controlled(PauliX,3,False)
NOTtt = controlled(PauliX,4,False)
midNOT = controlled(NOTt,3)
bigNOT = controlled(NOTt,4,False)

#Init the system
zeroqBit = qClass.qubit(1,0) #|0> state
totRegister = qClass.createRegister(np.full(7,zeroqBit)) 

circInit = np.full(7,None)
circInit[:3] = Hadamard
circuitMat =  [np.full(7,None),circInit]

#Add circuit steps corresponding to the modular operation
A = [None,None,NOTt,None,None]
B = [None,None,NOTtt,None]
C = [None,None,None,NOTt,None]
D = [None,midNOT,None]
E = C
F = [None,None,None,None,tNOT]
G = [None,bigNOT]
H = F
circuitMat = circuitMat + [A,B,C,D,E,F,G,H]

#Add circuit steps corresponding to QFT^-1
shiftOne = PhaseShift(phi = np.pi/2,angleRep = '90')
shiftOne = controlled(shiftOne,2,False)
shiftTwo = PhaseShift(phi = np.pi/4,angleRep = '45')
shiftTwo = controlled(shiftTwo,3,False)

invQFT = [[Hadamard,None,None],
          [shiftOne,None],
          [None,Hadamard,None],
          [shiftTwo],
          [None,shiftOne],
          [None,None,Hadamard]]
for i in range(6):invQFT[i] = invQFT[i] + [None,None,None,None]

circuitMat = circuitMat + invQFT + [np.full(7,None)]

Shor = qClass.qCircuit(circuitMat,7)
Shor.draw()

outComes = []
for i in range(1000):
    A = Shor.Simulate(totRegister)
    meas = '|' + qClass.Measurement(A)[-4:]
    outComes = outComes + [meas]

print('In 1000 simulations of the above circuit, measurements on the first 3 qubits were:')
print(Counter(outComes))



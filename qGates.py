#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spencer
"""

#Some frequently used Operators
#Useful Characters, · │ ─ ┌───┐ └───┘ ┤ ├  x ┼ ┬ ┴

import numpy as np
import copy
import qClass

Identity = np.asarray([[1,0],[0,1]])
def UnityRoot(N,k=1): return np.exp(2*k*np.pi*1j / N)
def nIdentity(n=2):
    outMat = Identity
    for i in range(1,n):outMat = np.kron(outMat,Identity)
    return outMat

HArray = np.asarray([[1,1],[1,-1]])

#Single Qubit
PauliX = qClass.genericOperator([[0,1],[1,0]],' X ') #equivalent to the not gate
PauliY = qClass.genericOperator([[0,-1j],[1j,0]],' Y ')
PauliZ = qClass.genericOperator([[1,0],[0,-1]],' Z ')
PauliI = qClass.genericOperator([[1,0],[0,1]],' I ')
Hadamard = qClass.genericOperator(HArray/np.sqrt(2),' H ')

def PhaseShift(phi = np.pi,angleRep = '_π'): return qClass.genericOperator([[1,0],[0,np.exp(phi*1j)]],'R' + angleRep)

#Two Qubit
Swap = qClass.operator([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
Swap.cRep = ['     ', '──x──', '  │  ', '  │  ',  '──x──', '     ']

CNot = qClass.operator([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
CNot.cRep = ['     ', '──·──', '  │  ', '┌─┴─┐',  '┤ X ├', '└───┘']

#>2 Qubit
def Nswap(dist=2):#swaps qubits a distance dist apart.
    if(dist<2 or type(dist)!=int):raise ValueError('Swaps occur between two qubits an integer distance >=2 apart')
        
    nLines = 2**dist
    outArray = []
    for i in range(nLines):
        newLine = np.zeros(nLines,int)
        if(i<nLines/2): newLine[2*i] = 1
        else: newLine[(2*i+1)%nLines] = 1
        outArray = outArray+[newLine]
       
    circRep = ['     ', '──x──', '  │  ']    
    for k in range(1,dist):
        if (k == dist-1): circRep = circRep + ['  │  ',  '──x──', '     ']
        else: circRep = circRep + ['  │  ',  '──┼──', '  │  ']
        
    outOp = qClass.operator(np.asarray(outArray))
    outOp.cRep = circRep
    return outOp
    
def controlled(Operator,dist=2,OpOnTop = True):#Takes an Operator and returns a controlled version of it
    Op = copy.deepcopy(Operator)
    numIn = Op.numInputs
    preMat = Op.matrix
    preRep = Op.cRep
    del(Op)
    
    Id = Identity #Identity part of the control gate
    for i in range(numIn-1):Id = np.kron(Identity,Id)
    outShape = tuple([2*x for x in Id.shape])
    outArr = np.zeros(outShape,complex)
    if(OpOnTop):
        outArr[:preMat.shape[0],:preMat.shape[1]] = preMat
        outArr[-Id.shape[0]:,-Id.shape[1]:] = Id
    else:
        outArr[:Id.shape[0],:Id.shape[1]] = Id
        outArr[-preMat.shape[0]:,-preMat.shape[1]:] = preMat
    
    if (dist!=2):
    
        Swap = Nswap(dist-1).matrix
        if(OpOnTop):Swap = np.kron(Id,Swap)
        else:Swap = np.kron(Swap,Id)
        m = dist-2 
        Idm = Identity 
        for i in range(m-1):Idm = np.kron(Identity,Idm)
        if(OpOnTop):outArr = np.kron(outArr,Idm)
        else:outArr = np.kron(Idm,outArr)
        
        outArr = np.matmul(Swap,outArr)
        outArr = np.matmul(outArr,Swap)
        
    #Now make the circuit representation
    newCircRep = []
    if OpOnTop:
        newCircRep = preRep
        if(preRep[-1] == '└───┘'):newCircRep[-1] = '└─┬─┘'
        else:newCircRep[-1] = '  │  '        
    else:newCircRep = ['     ', '──·──', '  │  ']
    
    for i in range(dist-2):
        newCircRep = newCircRep + ['  │  ',  '──┼──', '  │  ']

    if OpOnTop:newCircRep = newCircRep + ['  │  ', '──·──','     ' ]
    else:
        if(preRep[0] == '┌───┐'):preRep[0] = '┌─┴─┐'
        else:preRep[0] = '  │  '        
        newCircRep = newCircRep + preRep
    
    outOp = qClass.operator(np.asarray(outArr))
    outOp.cRep = newCircRep
    return outOp

Toffoli = controlled(CNot,OpOnTop=False)
Fredkin = controlled(Swap)
CiNOT = controlled(PauliX,3)
twoSwap = Nswap(3)

def QFT(n):#returns the Quantum Fourier Transform for n qubits
    #Note that this can also be constructed with controlled rotations and Hadamards
    N = 2**n
    outArr = np.empty([N,N],complex)
    
    Omega = UnityRoot(N)
    for i in range(N):
        for j in range(N):
            outArr[i][j] = Omega**(i*j)
    outArr = outArr / np.sqrt(N)
    
    outOp = qClass.genericOperator(outArr,'QFT')
    return outOp
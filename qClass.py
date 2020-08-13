#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spencer
"""

#Class Definitions
import numpy as np

def Norm(a):return np.sqrt(a.real**2+a.imag**2)
Identity = np.asarray([[1,0],[0,1]])

class qubit:
        #state = alpha|0> + beta|1>, alpha & beta normalized
        def __init__(self, alpha,beta):    
            alpha = alpha+0j
            beta = beta+0j
            normFactor = np.sqrt(Norm(alpha)**2+Norm(beta)**2)
            alpha = alpha/normFactor
            beta = beta/normFactor
            
            self.alpha = alpha
            self.beta = beta
            
        def getBloch(self): #Blochsphere rep, using only relative phase
            theta = 2*np.arccos(self.alpha)
            temp = self.beta / np.sin(self.theta/2)
            if(temp == -1 or temp == -1+0j):phi = np.pi
            else: phi=np.log(temp)/1j
            return [theta,phi]
    
        def setBloch(self,theta,relPhi):#set qubit with Bloch Sphere rep.
            self.alpha = np.cos(theta/2)
            self.beta = np.sin(theta/2)*np.exp(1j*relPhi)
            
        def addTotalPhase(self,psi): 
            self.alpha = self.alpha * np.exp(1j*psi)
            self.beta = self.beta * np.exp(1j*psi)
        
def createRegister(qubitArray): #qubit registers will be numpy vectors <- tensor products of qubits
    totalState = np.asarray([qubitArray[0].alpha,qubitArray[0].beta],complex)
    for i in range(1,len(qubitArray)):
        qState = np.asarray([qubitArray[i].alpha,qubitArray[i].beta],complex)
        totalState = np.kron(totalState,qState)
    return totalState
#Can 'combine' two registers with np.kron(_,_)

def currentProbs(Register):#For the whole register
    probabilityArr = []
    for state in Register:
        probabilityArr = probabilityArr + [Norm(state)**2]
    return probabilityArr

def Measurement(Register): #For the whole register
    numBits = round(np.log2(len(Register)))
    probArr = currentProbs(Register)
    probArr = probArr/sum(probArr)
    outComes = np.arange(len(Register),dtype=int)
    
    outCome = np.random.choice(outComes, p = probArr)
    outCome = bin(outCome)[2:]
    while(len(outCome)<numBits):outCome = '0'+outCome
    outCome = '|'+outCome+'>'
    return outCome
    
class operator:
    
        def __init__(self,matrix,cRep=['     ', '     ', '     ']):
            if (len(matrix) != len(matrix[0])):raise ValueError('Need a Square Matrix.')
            elif (len(matrix)%2==1):raise ValueError('Need an even lengthed matrix')
            else: self.matrix = np.asarray(matrix,complex)
            self.numInputs = round(np.log2(len(matrix)))
            
        def times(self,operator2):return np.matmul(self.matrix, operator2.matrix)
        
        def tensorProd(self,operator2):
            return np.kron(self.matrix,operator2.matrix)

            
class genericOperator(operator): #Construction Methods for a generic operator for easy drawing.
        
        def __init__(self,matrix,circuitLetter='   '):
            if(len(circuitLetter)>3):circuitLetter=circuitLetter[:3]
            if(len(circuitLetter)==2):circuitLetter=circuitLetter+' '
            if(len(circuitLetter)==1):circuitLetter=' '+circuitLetter+' '
            
            operator.__init__(self,matrix)    
            self.letter = circuitLetter
            self.cRep = self.createRep(circuitLetter,self.numInputs)
                
        def createRep(self,letter,numInputs):
            Rep=[]
            lenRep = 3*numInputs
            if(lenRep%2==0):letterPlacement = lenRep/2-1
            else:letterPlacement = (lenRep-1)/2
            
            for i in range(lenRep):
                lEdge,rEdge = "│","│"
                if(i%3==1): lEdge,rEdge = "┤","├"
                
                if(i == 0):line="┌───┐"
                elif(i == lenRep-1):line="└───┘"
                elif(i==letterPlacement):line=lEdge+letter+rEdge
                else:line=lEdge+"   "+rEdge
                Rep = Rep + [line]
            return Rep
            
class qCircuit:

        def __init__(self, elements,registerLen):#elements gives a list of steps, which are lists of operators or none type 
            #((i,j) here defines column THEN row)
            repArray=[]
            for line in range(len(elements)):
                currentSum = 0
                currentArray = []
                for Op in elements[line]:
                    if (Op == None):
                        currentSum +=1
                        currentArray = currentArray+['     ', '─────', '     '] 
                    else:
                        currentSum += Op.numInputs
                        currentArray = currentArray+Op.cRep
                if(currentSum != registerLen):raise ValueError('Circuit Element Size Mismatch')
                repArray = repArray+[currentArray]
            
            self.opArray = elements
            self.repArray = np.asarray(repArray).T
            self.numInputs = registerLen
        
        def draw(self,Register = None):
            if(len(self.repArray)==0):raise ValueError('Nothing to draw')
            outString=''
            for i in range(len(self.repArray)):
                if(Register == None):currentString=''
                else:currentString = Register[i]
                for x in self.repArray[i]:
                    currentString = currentString+x
                outString = outString+currentString+'\n'
            print(outString)
    
        def Simulate(self, Register): #Register is length 2**N for N qubits
            if(len(Register) != 2**self.numInputs):raise ValueError('Circuit does not work on a Register of this size')
            newRegister = Register
            for tStep in range(len(self.opArray)):
                
                if(self.opArray[tStep][0] == None):currentTotOp = Identity
                else:currentTotOp = self.opArray[tStep][0].matrix
                
                firstOp = True
                for op in self.opArray[tStep]:
                    if(firstOp):firstOp = False
                    elif(op == None):currentTotOp = np.kron(currentTotOp,Identity)
                    else:currentTotOp = np.kron(currentTotOp,op.matrix)
                    
                newRegister = np.matmul(currentTotOp,newRegister)
            
            return newRegister
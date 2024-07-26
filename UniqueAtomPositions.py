# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 09:30:07 2024

@author: jgwbe
"""

import numpy as np

def UniqueAtomPositions(RBasisOccupancy,RBasisIsoDW,RBasisAtomPosition,
        RSymMat,RSymVec,RaVecM,RbVecM,RcVecM,
        SBasisAtomLabel,SBasisAtomName,
        IBasisAtomicNumber,
        IBasisAnisoDW,IMaxPossibleNAtomsUnitCell,INAtomsUnitCell,
        IErr):
    
    IAtomCount=len(RBasisOccupancy)
    RAtomCoordinate=np.zeros((IAtomCount,3),dtype="float")
    RAllAtomPosition=np.zeros((IAtomCount,3),dtype="float")
    RAllOccupancy=np.zeros(IAtomCount,dtype="float")
    RAllIsoDW=np.zeros(IAtomCount,dtype="float")
    SAllAtomLabel=np.zeros(IAtomCount,dtype="str")
    SAllAtomName=np.zeros(IAtomCount,dtype="str")
    IAllAtomicNumber=np.zeros(IAtomCount,dtype="int")
    IAllAnIsoDW=np.zeros(IAtomCount,dtype="int")
    
    RAtomPosition=np.zeros((IAtomCount,3),dtype="float")
    ROccupancy=np.zeros(IAtomCount,dtype="float")
    RIsoDW=np.zeros(IAtomCount,dtype="float")
    SAtomLabel=np.zeros(IAtomCount,dtype="str")
    SAtomName=np.zeros(IAtomCount,dtype="str")
    IAtomicNumber=np.zeros(IAtomCount,dtype="int")
    IAnIsoDW=np.zeros(IAtomCount,dtype="int")
    
    for i in range(IAtomCount):
        RAllAtomPosition=np.matmul(RSymMat[i,:,:],RBasisAtomPosition[i,:])+RSymVec[i,:]
        RAllOccupancy[i]=RBasisOccupancy[i]
        RAllIsoDW[i]=RBasisIsoDW[i]
        SAllAtomLabel[i]=SBasisAtomLabel[i]
        SAllAtomName[i]=SBasisAtomName[i]
        IAllAtomicNumber[i]=IBasisAtomicNumber[i]
        IAllAnIsoDW[i]=RBasisIsoDW[i]
    
    RAllAtomPosition=np.mod(RAllAtomPosition,1)
    
    RAtomPosition[0,:]=RAllAtomPosition[0,:]
    ROccupancy[0]=RAllOccupancy[0]
    RIsoDW[0]=RAllIsoDW[0]
    SAtomLabel[0]=SAllAtomLabel[0]
    SAtomName[0]=SAllAtomName[0]
    IAtomicNumber[0]=IAllAtomicNumber[0]
    IAnIsoDW[0]=IAllAnIsoDW[0]
    
    j=1
    for i in range(1,IMaxPossibleNAtomsUnitCell):
        Lunique=True
        for k in range(0,j-1):
            if np.sum(np.abs(RAllAtomPosition[i,:]-RAtomPosition[k,:])):
                if SAllAtomLabel[i]==SAtomLabel[k]:
                    Lunique=False
        if Lunique==True:
            RAtomPosition[j,:]=RAllAtomPosition[i,:]
            ROccupancy[j]=RAllOccupancy[i]
            RIsoDW[j]=RAllIsoDW[i]
            SAtomLabel[j]=SAllAtomLabel[i]
            SAllAtomName[j]=SAllAtomName[i]
            IAtomicNumber[j]=IAtomicNumber[i]
            IAnIsoDW[j]=IAllAnIsoDW[i]
    INAtomsUnitCell=j-1
    
    for i in range(INAtomsUnitCell):
        print("Atom:",str(i+1))
        print("Atom Position:",RAtomPosition[i,:])
        print(SAtomName)
        print("DWF, Occupancy:",RIsoDW[i],ROccupancy[i])
        print("")
    
    for j in range(INAtomsUnitCell):
        for i in range(0,3):
            RAtomCoordinate[j][i]=RAtomPosition[j][0]*RaVecM[i]+RAtomPosition[j][1]*RbVecM[i]+RAtomPosition[j][2]*RcVecM[i]
    
    
    return RAtomCoordinate
    
    
    
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 09:30:07 2024

@author: Jacob Watkiss
"""

"""
-----------------------------------------------
Imports any necessary packages, libraries, modules etc.
-----------------------------------------------
"""
import numpy as np
import globalvariables as g
IErr=g.IErr

"""
-----------------------------------------------
Defines the subprogram
-----------------------------------------------
"""
def UniqueAtomPositions(RBasisOccupancy,RBasisIsoDW,RBasisAtomPosition,
        RSymMat,RSymVec,RaVecM,RbVecM,RcVecM,
        SBasisAtomLabel,SBasisAtomName,
        IBasisAtomicNumber,
        IBasisAnisoDW,IMaxPossibleNAtomsUnitCell,INAtomsUnitCell):
    
    """
    -----------------------------------------------
    Declaring variables
    -----------------------------------------------
    """
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
    
    """
    -----------------------------------------------
    Create dummy arrays that can be modified and 
    compared to the original arrays of values
    -----------------------------------------------
    """
    for i in range(IAtomCount):
        """Apply symmetry to atom position
        to generate all equivalent positions"""
        RAllAtomPosition=np.matmul(RSymMat[i,:,:],RBasisAtomPosition[i,:])+RSymVec[i,:]
        RAllOccupancy[i]=RBasisOccupancy[i]
        RAllIsoDW[i]=RBasisIsoDW[i]
        SAllAtomLabel[i]=SBasisAtomLabel[i]
        SAllAtomName[i]=SBasisAtomName[i]
        IAllAtomicNumber[i]=IBasisAtomicNumber[i]
        IAllAnIsoDW[i]=RBasisIsoDW[i]
    
    """Ensure that atom positions are in
    between 0 and 1 (to find distance from
    next atom)"""
    RAllAtomPosition=np.mod(RAllAtomPosition,1)
    
    """
    -----------------------------------------------
    Reduce to the set of unique fractional atomic positions
    -----------------------------------------------
    """
    """First atom has to be in this set"""
    RAtomPosition[0,:]=RAllAtomPosition[0,:]
    ROccupancy[0]=RAllOccupancy[0]
    RIsoDW[0]=RAllIsoDW[0]
    SAtomLabel[0]=SAllAtomLabel[0]
    SAtomName[0]=SAllAtomName[0]
    IAtomicNumber[0]=IAllAtomicNumber[0]
    IAnIsoDW[0]=IAllAnIsoDW[0]
    
    """Work through all possible atom coordinates
    and check for duplicates"""
    j=1
    for i in range(1,IMaxPossibleNAtomsUnitCell):
        Lunique=True
        """Check against unique ones found so far"""
        for k in range(0,j-1):
            """If the position is the same"""
            if np.sum(np.abs(RAllAtomPosition[i,:]-RAtomPosition[k,:])<=10**-9):
                """And the label is too, they are not unique"""
                if SAllAtomLabel[i]==SAtomLabel[k]:
                    Lunique=False
        """If they are unique, add to the array"""
        if Lunique==True:
            RAtomPosition[j,:]=RAllAtomPosition[i,:]
            ROccupancy[j]=RAllOccupancy[i]
            RIsoDW[j]=RAllIsoDW[i]
            SAtomLabel[j]=SAllAtomLabel[i]
            SAllAtomName[j]=SAllAtomName[i]
            IAtomicNumber[j]=IAtomicNumber[i]
            IAnIsoDW[j]=IAllAnIsoDW[i]
    """Number of unique atoms"""
    INAtomsUnitCell=j-1
    

    for i in range(INAtomsUnitCell):
        print("Atom:",str(i+1))
        print("Atom Position:",RAtomPosition[i,:])
        print(SAtomName)
        print("DWF, Occupancy:",RIsoDW[i],ROccupancy[i])
        print("")
    
    """
    -----------------------------------------------
    Calculate atomic position vectors RAtomCoordinate
    -----------------------------------------------
    """
    """Calculate from Fractional Coordinate and 
    Lattice Vectors"""
    """In microscope reference frame, in Angstrom units"""
    for j in range(INAtomsUnitCell):
        for i in range(0,3):
            RAtomCoordinate[j][i]=RAtomPosition[j][0]*RaVecM[i]+RAtomPosition[j][1]*RbVecM[i]+RAtomPosition[j][2]*RcVecM[i]
    
    """Return new array with unique atom coordinates"""
    return RAtomCoordinate
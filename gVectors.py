# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:50:57 2024

@author: jgwbe
"""

"""
-----------------------------------------------
Imports any necessary packages, libraries, modules etc.
-----------------------------------------------
"""
import numpy as np
import ReciprocalLattice as RL

"""
-----------------------------------------------
Define the function
-----------------------------------------------
"""
def gVectors(RarVecM,RbrVecM,RcrVecM,RNormDirM,Rhkl,RConvergenceAngle,
             INhkl,IPixelCount,
             IErr):
    
    """
    -----------------------------------------------
    Declare/call the variables
    -----------------------------------------------
    """
    RgPool=np.zeros((INhkl,3),dtype="float")
    RgPoolMag=np.zeros(INhkl,dtype="float")
    RgDotNorm=np.zeros(INhkl,dtype="float")
    RgMatrix=np.zeros((INhkl,INhkl,3),dtype="float")
    
    """
    -----------------------------------------------
    Calculate the g-vector pool, the magnitudes and component
    parallel to specimen surface
    -----------------------------------------------
    """
    for j in range(0,INhkl):
        for i in range(0,3):
            RgPool[j][i]=(Rhkl[j][0]*RarVecM[i])+(Rhkl[j][1]*RbrVecM[i])+(Rhkl[j][2]*RcrVecM[i])
        RgPoolMag[j]=np.sqrt(np.dot(RgPool[j][:],RgPool[j][:]))
        RgDotNorm[j]=np.dot(RgPool[j][:],RNormDirM)
    
    """
    -----------------------------------------------
    Calculate matrix of g-vectors that corresponds to Ug matrix
    -----------------------------------------------
    """
    for j in range(0,INhkl):
        for i in range(0,INhkl):
            RgMatrix[j][i][:]=RgPool[j][:]-RgPool[i][:]
    
    """
    -----------------------------------------------
    Outputs and returns all
    -----------------------------------------------
    """
    print("First 16 g-vectors",RgMatrix[0:17][0][:])
    print("g-vectors and magnitude (1/A), in the microscope reference frame")
    for j in range(0,INhkl):
        print("hkl:",Rhkl[j][:])
        print("g mag:",RgPoolMag[j])
    print("g.n list")
    for j in range(0,INhkl):
        print("hkl:",Rhkl[j][:])
        print("g.n:",RgDotNorm[j])
    
    return RgPool,RgPoolMag,RgDotNorm,RgMatrix
    


"""
-----------------------------------------------
Call the function
-----------------------------------------------
"""









# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 12:00:35 2024

@author: Jacob Watkiss
"""

"""
-----------------------------------------------
Imports any necessary packages, libraries, modules etc.
-----------------------------------------------
"""
import numpy as np
from read_cif_mod import read_cif
import re
import globalvariables as g
IErr=g.IErr

"""
-----------------------------------------------
Define the function
-----------------------------------------------
"""
def ReciprocalLattice(IDiffractionFLAG,
                      RLengthX,RLengthY,RLengthZ,RAlpha,RBeta,RGamma,RVolume,RNormDirC,RXDirC,RZDirC,
                      SSpaceGroupName):
    
    """
    -----------------------------------------------
    Define and declare variables
    -----------------------------------------------
    """
    tiny=np.finfo(np.float32).tiny
    RarVecO=np.ones(3,dtype="float")
    RbrVecO=np.ones(3,dtype="float")
    RcrVecO=np.ones(3,dtype="float")
    RXDirO=np.ones(3,dtype="float")
    RYDirO=np.ones(3,dtype="float")
    RZDirO=np.ones(3,dtype="float")
    RTMatC2O=np.zeros((3,3),dtype="float")
    RTMatO2M=np.zeros((3,3),dtype="float")
    
    
    """
    -----------------------------------------------
    Direct lattice vectors in an orthogonal reference
    frame, Angstrom units
    -----------------------------------------------
    """
    RaVecO=np.array((RLengthX,0,0),dtype="float")
    RbVecO=np.array((RLengthY*g.cos(RGamma),RLengthY*g.sin(RGamma),0),dtype="float")
    RcVecO=np.array((RLengthZ*g.cos(RBeta),RLengthZ*((g.cos(RAlpha)-g.cos(RBeta)*g.cos(RGamma))/g.sin(RGamma)),(RLengthZ*RVolume)/(g.sin(RGamma)*RLengthX*RLengthY*RLengthZ)),dtype="float")
    
    """
    -----------------------------------------------
    "Some checks for rhombohedral cells?"?
    -----------------------------------------------
    """
    if IDiffractionFLAG==0:
        aMod=np.sqrt(np.dot(RaVecO,RaVecO))
        bMod=np.sqrt(np.dot(RbVecO,RbVecO))
        cMod=np.sqrt(np.dot(RcVecO,RcVecO))
        RTTest=(aMod*bMod*cMod)**(-4)*np.dot(RaVecO,RbVecO)*np.dot(RbVecO,RcVecO)*np.dot(RcVecO,RaVecO)
        
        if re.search("rR",SSpaceGroupName):
            if abs(RTTest)<tiny:
                SSpaceGroupName="V"
                """Crystal is either Obverse or Reverse
                Selection Rules are not in place to determine
                the difference, assume the crystal is Obverse"""
            else:
                SSpaceGroupName="P"
                """Primitive setting (Rhombohedral axes)"""
    
    """
    -----------------------------------------------
    Set up Reciprocal Lattice Vectors: orthogonal reference
    frame in 1/Angstrom units
    -----------------------------------------------
    """
    """
    RarDirO,RbrDirO,RcrDirO vectors are reiprocal lattice vectors
    2pi/a, 2pi/b, 2pi/c in an orthogonal frame
    Note that reciprocal lattice vectors have two pi included,
    we are using the Physics convention exp(i*g.r)
    """
    RarVecO=2*np.pi*np.cross(RbVecO,RcVecO)/np.dot(RbVecO,np.cross(RcVecO,RaVecO))
    RbrVecO=2*np.pi*np.cross(RcVecO,RaVecO)/np.dot(RcVecO,np.cross(RaVecO,RbVecO))
    RcrVecO=2*np.pi*np.cross(RaVecO,RbVecO)/np.dot(RaVecO,np.cross(RbVecO,RcVecO))
    
    for i in range(0,3):
        if abs(RarVecO[i]<tiny):
            RarVecO[i]=0
        if abs(RbrVecO[i]<tiny):
            RbrVecO[i]=0
        if abs(RcrVecO[i]<tiny):
            RcrVecO[i]=0
       
    """
    RTmat transforms from crystal (implicit units)
    to orthogonal reference frame (Angstrom)
    """
    for i in range(len(RaVecO)):
        RTMatC2O[i][0]=RaVecO[i]
        RTMatC2O[i][1]=RbVecO[i]
        RTMatC2O[i][2]=RcVecO[i]
        
    """
    RXDirC is the reciprocal lattice vector that defines the
    x-axis of the diffraction pattern and RZDirC the beam
    direction, coming from felix.inp. No check has been made
    to ensure that they are perpendicular, it is assumed.
    RXDirO,RYDirO,RZDirO vectors are UNIT reciprocal lattice
    vectors parallel to the above in an orthogonal frame
    """
    RXDirO=(RXDirC[0]*RarVecO)+(RXDirC[1]*RbrVecO)+(RXDirC[2]*RcrVecO)
    RXDirOMod=np.sqrt(np.dot(RXDirO,RXDirO))
    RXDirO=RXDirO/RXDirOMod
    RZDirO=(RZDirC[0]*RaVecO)+(RZDirC[1]*RbVecO)+(RZDirC[2]*RcVecO)
    RZDirOMod=np.sqrt(np.dot(RZDirO,RZDirO))
    RZDirO=RZDirO/RZDirOMod
    RYDirO=np.cross(RZDirO,RXDirO)
    
    """RTMatO2M transforms from orthogonal to microscope reference frame"""
    for i in range(len(RXDirO)):
        RTMatO2M[0][i]=RXDirO[i]
        RTMatO2M[1][i]=RYDirO[i]
        RTMatO2M[2][i]=RZDirO[i]
    
    """
    Unit normal to the specimen in REAL space
    This is used in diffraction pattern calculation
    """
    RNormDirM=np.matmul(RTMatO2M,np.matmul(RTMatC2O,RNormDirC))
    RNormDirMMod=np.sqrt(np.dot(RNormDirM,RNormDirM))
    RNormDirM=RNormDirM/RNormDirMMod
    
    """
    -----------------------------------------------
    Now transform from crystal reference frame to orthogonal 
    and then to microscope frame
    -----------------------------------------------
    """
    """
    RaVecM,RbVecM,RcVecM unit cell vectors in Angstrom
    units in the microscope frame
    """
    RaVecM=np.matmul(RTMatO2M,RaVecO)
    RbVecM=np.matmul(RTMatO2M,RbVecO)
    RcVecM=np.matmul(RTMatO2M,RcVecO)
    
    """
    -----------------------------------------------
    Create new set of reciprocal lattice vectors in microscope
    reference frame.
    Note that reciprocal lattice vectors have two pi included
    We are using the optical convention exp(i*g.r)
    -----------------------------------------------
    """

    RarVecM=2*np.pi*np.cross(RbVecM,RcVecM)/np.dot(RbVecM,np.cross(RcVecM,RaVecM))
    RbrVecM=2*np.pi*np.cross(RcVecM,RaVecM)/np.dot(RcVecM,np.cross(RaVecM,RbVecM))
    RcrVecM=2*np.pi*np.cross(RaVecM,RbVecM)/np.dot(RaVecM,np.cross(RbVecM,RcVecM))
    
    return (RTTest,SSpaceGroupName,
            RaVecO,RbVecO,RcVecO,RarVecO,RbrVecO,RcrVecO,RXDirO,RYDirO,RZDirO,
            RaVecM,RbVecM,RcVecM,RarVecM,RbrVecM,RcrVecM,RNormDirM,
            RTMatC2O,RTMatO2M
        )
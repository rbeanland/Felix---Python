# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:30:09 2024

@author: Jacob Watkiss
"""

"""
-----------------------------------------------
Imports any necessary packages, libraries, modules etc.
-----------------------------------------------
"""
import numpy as np
import re
from read_cif_mod import find
import globalvariables as g
IErr=g.IErr

"""
-----------------------------------------------
Defines the functions
-----------------------------------------------
"""
def strToReal(var):
    """
    -----------------------------------------------
    Cleans the values that include non-numerical 
    characters and turns them into individual floats
    ([1,2,3] --> 1 2 3)
    -----------------------------------------------
    """
    array=np.zeros(3,dtype="float")
    var=var.replace("[","")
    var=var.replace("]","")
    comma=re.search(",",var).start()
    array[0]=var[0:comma]
    var=var[comma+1:len(var)]
    comma=re.search(",",var).start()
    array[1]=var[0:comma]
    var=var[comma+1:len(var)]
    array[2]=var
    return array

def ReadInpFile():
    """
    -----------------------------------------------
    Declares variables of values to be extracted
    -----------------------------------------------
    """
    IWriteFLAG=0
    IScatterFactorMethodFLAG=0
    IHolzFLAG=0
    IAbsorbFLAG=0
    IByteSize=0
    INhkl=0
    IMinStrongBeams=0
    IMinWeakBeams=0
    ISimFLAG=0
    IRefineModeFLAG=0
    IWeightingFLAG=0
    IRefineMethodFLAG=0
    ICorrelationFLAG=0
    IImageProcessingFLAG=0
    INoofUgs=0
    IPrint=0
    IPixelCount=0
    RDebyeWallerConstant=0.0
    RAbsorptionPercentage=0.0
    RConvergenceAngle=0.0
    RZDirC=np.zeros(3,dtype="float")
    RXDirC=np.zeros(3,dtype="float")
    RNormDirC=np.zeros(3,dtype="float")
    RAcceleratingVoltage=0.0
    RAcceptanceAngle=0.0
    RInitialThickness=0.0
    RFinalThickness=0.0
    RDeltaThickness=0.0
    RBlurRadius=0.0
    RSimplexLengthScale=0.0
    RExitCriteria=0.0
    RPrecision=0.0
    RgLimit=0.0
    
    """
    -----------------------------------------------
    The types of refinement modes and methods
    -----------------------------------------------
    """
    modes=["Refining Structure Factors, A",
           "Refining Atomic Coordinates, B",
           "Refining Occupancies, C",
           "Refining Isotropic Debye Waller Factors, D",
           "Refining Anisotropic Debye Waller Factors, E",
           "Refining Lattice Parameters, F",
           "Refining Lattice Angles, G",
           "Refining Convergence Angle, H",
           "Refining Accelerating Voltage, I"
           ]
    SRefineMode=np.array(modes,dtype="str")
    methods=["Refining by simplex",
             "Refining by downhill (2-point) gradient",
             "Refining by maximum (3-point) gradient",
             "Refining by pairwise (2x2-point) gradient"
             ]
    SRefineMethod=np.array(methods,dtype="str")
    
    """
    -----------------------------------------------
    Opens the file
    -----------------------------------------------
    """
    filename="felix.inp"
    name=open(filename,"r")
    content=name.readlines()
    
    """
    -----------------------------------------------
    Extracts all of the simple values by using the 
    find function from read_cif_mod
    -----------------------------------------------
    """
    IWriteFLAG=find(content,"IWriteFLAG","int")
    if IWriteFLAG=="Absent":
        print("No IWriteFLAG found")
    IScatterFactorMethodFLAG=find(content,"IScatterFactorMethodFLAG","int")#
    if IScatterFactorMethodFLAG=="Absent":
        print("No IScatterFactorMethodFLAG found")
    IHolzFLAG=find(content,"IHolzFLAG","int")
    if IHolzFLAG=="Absent":
        print("No IHolzFLAG found")
    IAbsorbFLAG=find(content,"IAbsorbFLAG","int")
    if IAbsorbFLAG=="Absent":
        print("No IAbsorbFLAG found")
    IByteSize=find(content,"IByteSize","int")
    if IByteSize=="Absent":
        print("No IByteSize found")
    INhkl=find(content,"IMinReflectionPool","int")
    if INhkl=="Absent":
        print("No INhkl found")
    IMinStrongBeams=find(content,"IMinStrongBeams","int")
    if IMinStrongBeams=="Absent":
        print("No IMinStrongBeams found")
    IMinWeakBeams=find(content,"IMinWeakBeams","int")
    if IMinWeakBeams=="Absent":
        print("No IMinWeakBeams found")
    IRefineModeFLAG=find(content,"IRefineMode","int")
    if IRefineModeFLAG=="Absent":
        print("No IRefineModeFLAG found")
    else:
        print(SRefineMode[IRefineModeFLAG])
    IWeightingFLAG=find(content,"IWeightingFLAG","int")
    if IWeightingFLAG=="Absent":
        print("No IWeightingFLAG found")
    IRefineMethodFLAG=find(content,"IRefineMethodFLAG","int")
    if IRefineMethodFLAG=="Absent":
        print("No IRefineMethodFLAG found")
    else:
        print(SRefineMethod[IRefineMethodFLAG])
    ICorrelationFLAG=find(content,"ICorrelationFLAG","int")
    if ICorrelationFLAG=="Absent":
        print("No ICorrelationFLAG found")
    IImageProcessingFLAG=find(content,"IImageProcessingFLAG","int")
    if IImageProcessingFLAG=="Absent":
        print("No IImageProcessingFLAG found")
    INoofUgs=find(content,"INoofUgs","int")
    if INoofUgs=="Absent":
        print("No INoofUgs found")
    IPrint=find(content,"IPrint","int")
    if IPrint=="Absent":
        print("No IPrint found")
    IPixelCount=find(content,"IPixelCount","int")
    if IPixelCount=="Absent":
        print("No IPixelCount found")
    RDebyeWallerConstant=find(content,"RDebyeWallerConstant","real")
    if RDebyeWallerConstant=="Absent":
        print("No RDebyeWallerConstant found")
    RAbsorptionPercentage=find(content,"RAbsorptionPer","real")
    if RAbsorptionPercentage=="Absent":
        print("No RAbsorptionPercentage found")
    RConvergenceAngle=find(content,"RConvergenceAngle","real")
    if RConvergenceAngle=="Absent":
        print("No RConvergenceAngle found")
    temp1=find(content,"IIncidentBeamDirection","str")
    if temp1=="Absent":
        print("No RZDirC found")
    RZDirC=strToReal(temp1)
    temp2=find(content,"IXDirection","str")
    if temp2=="Absent":
        print("No RXDirC found")
    RXDirC=strToReal(temp2)
    temp3=find(content,"INormalDirection","str")
    if temp3=="Absent":
        print("No RNormDirC found")
    RNormDirC=strToReal(temp3)
    RAcceleratingVoltage=find(content,"RAcceleratingVoltage","real")
    if RAcceleratingVoltage=="Absent":
        print("No RAcceleratingVoltage found")
    RAcceptanceAngle=find(content,"RAcceptanceAngle","real")
    if RAcceptanceAngle=="Absent":
        print("No RAcceptanceAngle found")
    RInitialThickness=find(content,"RInitialThickness","real")
    if RInitialThickness=="Absent":
        print("No RInitialThickness found")
    RFinalThickness=find(content,"RFinalThickness","real")
    if RFinalThickness=="Absent":
        print("No RFinalThickness found")
    RDeltaThickness=find(content,"RDeltaThickness","real")
    if RDeltaThickness=="Absent":
        print("No RDeltaThickness found")
    RBlurRadius=find(content,"RBlurRadius","real")
    if RBlurRadius=="Absent":
        print("No RBlurRadius found")
    RSimplexLengthScale=find(content,"RSimplexLengthScale","real")
    if RSimplexLengthScale=="Absent":
        print("No RSimplexLengthScale found")
    RExitCriteria=find(content,"RExitCriteria","real")
    if RExitCriteria=="Absent":
        print("No RExitCriteria found")
    
    """
    -----------------------------------------------
    I am not sure where to find these in the .inp file,
    so the keys have been presumed --> please update
    with correct keys
    -----------------------------------------------
    """
    ISimFLAG=find(content,"ISimFLAG","int")
    if ISimFLAG=="Absent":
        print("No ISimFLAG found")
    RPrecision=find(content,"RPrecision","real")
    if RPrecision=="Absent":
        print("No RPrecision found")
    RgLimit=find(content,"RgLimit","real")
    if RgLimit=="Absent":
        print("No RgLimit found: RgLimit=INhkl")
        
    """
    -----------------------------------------------
    Closing the file
    -----------------------------------------------
    """
    name.close()
    
    """
    -----------------------------------------------
    Return all values
    -----------------------------------------------
    """
    return(IWriteFLAG,IScatterFactorMethodFLAG,IHolzFLAG,IAbsorbFLAG,IByteSize,
           INhkl,IMinStrongBeams,IMinWeakBeams,
           IRefineMethodFLAG,IRefineModeFLAG,IWeightingFLAG,ICorrelationFLAG,IImageProcessingFLAG,INoofUgs,
           IPrint,IPixelCount,
           RDebyeWallerConstant,RAbsorptionPercentage,
           RConvergenceAngle,RZDirC,RXDirC,RNormDirC,RAcceleratingVoltage,RAcceptanceAngle,
           RInitialThickness,RFinalThickness,RDeltaThickness,
           RBlurRadius,RSimplexLengthScale,RExitCriteria,
           RPrecision,RgLimit,ISimFLAG
        )

def ReadHklFile():
    """
    -----------------------------------------------
    Declares variables
    -----------------------------------------------
    """
    INoOfLacbedPatterns=0
    IOutputReflections=0
    RInputHKLs=0.0
    
    """
    -----------------------------------------------
    Opens the file
    -----------------------------------------------
    """
    filename="felix.hkl"
    name=open(filename,"r")
    content=name.readlines()
    
    """
    -----------------------------------------------
    Counts the amount of lines in the file
    --> this is the number of reflections to output
    -----------------------------------------------
    """
    INoOfLacbedPatterns=len(content)
    IOutputReflections=np.array(INoOfLacbedPatterns,dtype="int")
    RInputHKLs=np.zeros((INoOfLacbedPatterns,3),dtype="float")
    
    """
    -----------------------------------------------
    Read in the HKLs
    -----------------------------------------------
    """
    for x in range(INoOfLacbedPatterns):
        hkl=np.zeros(3,dtype="float")
        line=content[x]
        end=re.search("]",line).end()
        line=line[0:end]
        hkl=strToReal(line)
        RInputHKLs[x][0]=hkl[0]
        RInputHKLs[x][1]=hkl[1]
        RInputHKLs[x][2]=hkl[2]
    
    """
    -----------------------------------------------
    Closing the file
    -----------------------------------------------
    """
    name.close()
    
    """
    -----------------------------------------------
    Return values
    -----------------------------------------------
    """
    return(INoOfLacbedPatterns,IOutputReflections,RInputHKLs)
    
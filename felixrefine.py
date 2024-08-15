# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:07:16 2024

@author: Jacob Watkiss
"""

import numpy as np
from read_cif_mod import read_cif
from read_files_mod import ReadInpFile,ReadHklFile
from ReciprocalLattice import ReciprocalLattice
from gVectors import gVectors
from setup_reflections_mod import HKLMake
import globalvariables as g

IErr=g.IErr



IWriteFLAG,IScatterFactorMethodFLAG,IHolzFLAG,IAbsorbFLAG,IByteSize,INhkl,IMinStrongBeams,IMinWeakBeams,IRefineMethodFLAG,IRefineModeFLAG,IWeightingFLAG,ICorrelationFLAG,IImageProcessingFLAG,INoofUgs,IPrint,IPixelCount,RDebyeWallerConstant,RAbsorptionPercentage,RConvergenceAngle,RZDirC,RXDirC,RNormDirC,RAcceleratingVoltage,RAcceptanceAngle,RInitialThickness,RFinalThickness,RDeltaThickness,RBlurRadius,RSimplexLengthScale,RExitCriteria,RPrecision,RgLimit,ISimFLAG=ReadInpFile()
if RgLimit=="Absent" or RgLimit<=0:
    RgLimit=INhkl
INoOfLacbedPatterns,IOutputReflections,RInputHKLs=ReadHklFile() 
ISpaceGrp,IAtomCount,ILN,IVolumeFLAG,SBasisAtomLabel,SBasisAtomName,SSpaceGroupName,SSpaceGrp,SChemicalFormula,SSymString,RLengthX,RLengthY,RLengthZ,RAlpha,RBeta,RGamma,RVolume,RBasisAtomPosition,RBasisIsoDW,RBasisOccupancy,IBasisAtomicNumber,RSymMat,RSymVec=read_cif(IErr)







"""Some Maths to calculate velocity, wavelength and wavenumber"""

"""
To make the calculations easier to understand
(after all, that is one of the main reasons
 felix is being brought to Python), I have 
defined the values with their symbol counterparts
and calculated them from there.
"""
e=g.RElectronCharge
V=RAcceleratingVoltage
m=g.RElectronMass
c=g.RSpeedOfLight
h=g.RPlanckConstant
ang=g.RAngstromConversion

fraction=m*c**2/((e*V*1000)+(m*c**2))
RElectronVelocity=c*np.sqrt(1-fraction**2)

RElectronWavelength=h/np.sqrt((2*m*e*V*1000)+(e*V*1000/c)**2)
RElectronWavelength=RElectronWavelength*ang

RElectronWaveVectorMagnitude=2*np.pi/RElectronWavelength



RTTest,SSpaceGroupName,RaVecO,RbVecO,RcVecO,RarVecO,RbrVecO,RcrVecO,RXDirO,RYDirO,RZDirO,RaVecM,RbVecM,RcVecM,RarVecM,RbrVecM,RcrVecM,RNormDirM,RTMatC2O,RTMatO2M=ReciprocalLattice(0,RLengthX,RLengthY,RLengthZ,RAlpha,RBeta,RGamma,RVolume,RNormDirC,RXDirC,RZDirC,SSpaceGroupName)

"""
print("\nLengths:\nx",RLengthX,"\ny",RLengthY,"\nz",RLengthZ)
print("\nAngles in radians:\nAlpha",RAlpha,"\nBeta",RBeta,"\nGamma",RGamma)
print("\nAngles in degrees:\nAlpha",(RAlpha*180/np.pi),"\nBeta",(RBeta*180/np.pi),"\nGamma",(RGamma*180/np.pi))
print("\nDirCs:\nX",RXDirC,"\nNorm",RNormDirC,"\nZ",RZDirC)
print("\nVecOs:\na",RaVecO,"\nb",RbVecO,"\nc",RcVecO)
print("\nrVecOs:\na",RarVecO,"\nb",RbrVecO,"\nc",RcrVecO)
print("\nRTMatC2O:\n",RTMatC2O)
print("\nDirOs:\na",RXDirO,"\nb",RYDirO,"\nc",RZDirO)
print("\nRTMatO2M\n",RTMatO2M)
print("\nVecMs:\na",RaVecM,"\nb",RbVecM,"\nc",RcVecM)
print("\nrVecMs:\na",RarVecM,"\nb",RbrVecM,"\nc",RcrVecM)
"""

#INhkl,Rhkl=HKLMake(RZDirC,RarVecM,RbrVecM,RcrVecM,RInputHKLs,RElectronWaveVectorMagnitude,RgLimit,SSpaceGroupName,INhkl,INoOfLacbedPatterns,IMinStrongBeams)
print("\nHKLMake: Successful!")
#RgPool,RgPoolMag,RgDotNorm,RgMatrix=gVectors(RarVecM,RbrVecM,RcrVecM,RNormDirM,Rhkl,RConvergenceAngle,INhkl,IPixelCount)

#print(Rhkl)

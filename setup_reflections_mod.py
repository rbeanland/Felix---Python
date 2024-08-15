# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 15:13:47 2024

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
Defines SelectionRules which checks if the 
hkl values agree with expected crystal space groups
-----------------------------------------------
"""
def SelectionRules(SSpaceGroupName,Ih,Ik,Il,ISel):
    """F - Face Centred: All odd or all even"""
    """I - Body Centred"""
    """A - A-Face Centred"""
    """B - B-Face Centred"""
    """C - C-Face Centred"""
    """R - Rhombohedral Reverse"""
    """V - Rhombohedral Obverse"""
    """P - Primitive"""
    if SSpaceGroupName=="F": 
        if (np.mod(Ih+Ik,2)==0) and (np.mod(Ik+Il,2)==0) and (np.mod(Il+Ih,2)==0):
            ISel=1
    elif SSpaceGroupName=="I":
        if np.mod(Ih+Ik+Il,2)==0:
            ISel=1
    elif SSpaceGroupName=="A":
        if np.mod(Ik+Il,2)==0:
            ISel=1
    elif SSpaceGroupName=="B":
        if np.mod(Ih+Il,2)==0:
            ISel=1
    elif SSpaceGroupName=="C":
        if np.mod(Ih+Ik,2)==0:
            ISel=1
    elif SSpaceGroupName=="R":
        if np.mod(Ih-Ik+Il,3)==0:
            ISel=1
    elif SSpaceGroupName=="V":
        if np.mod(-Ih+Ik+Il,3)==0:
            ISel=1
    elif SSpaceGroupName=="P":
        ISel=1
    else:
        ISel=0
        IErr=1
        print("Error: SSpaceGroupName unrecognised")
    return ISel

"""
-----------------------------------------------
Defines HKLMake which fills the list of reciprocal
space vectors Rhkl
-----------------------------------------------
"""
def HKLMake(RZDirC,RarVecM,RbrVecM,RcrVecM,
            RInputHKLs,RElectronWaveVectorMagnitude,RgLimit,
            SSpaceGroupName,
            INhkl,INoOfLacbedPatterns,IMinStrongBeams):
    
    """A temporary function to write to"""
    RhklTemp=np.zeros((66666,3),dtype="float")
    
    """
    The upper limit for g-vector magnitudes
    Perhaps the tolerance for proximity to 
    the Ewald sphere needs increasing
    """
    """The k-vector for the incident beam"""
    """We are working in the microscope
    frame reference so k is along z"""
    Rk=(0.0,0.0,RElectronWaveVectorMagnitude)
    """Get the size of the reciprocal lattice basis vectors"""
    RarMag=np.sqrt(np.dot(RarVecM,RarVecM))
    RbrMag=np.sqrt(np.dot(RbrVecM,RbrVecM))
    RcrMag=np.sqrt(np.dot(RcrVecM,RcrVecM))
    """We work our way out from the origin in shells"""
    """The shell increment is the smallest basis vector"""
    RShell=np.min((RarMag,RbrMag,RcrMag))
    
    """First g is always 0 0 0"""
    RhklTemp[0,:]=np.zeros(3,dtype="float")
    """Number of reflections in the pool"""
    reflNum=1
    """Number of the shell"""
    shellNum=0
    """
    Maximum a*,b*,c* limit is determined by 
    the G magnitude limit
    """
    if RgLimit<10**-9:
        RgLimit=20*np.pi
    else:
        """Make INhkl large and RgLimit will be the cutoff"""
        INhkl=66666
    
    a=int(np.rint(RgLimit/RarMag))
    b=int(np.rint(RgLimit/RbrMag))
    c=int(np.rint(RgLimit/RcrMag))
    
    """Fill the Rhkl with beams near the Bragg condition"""
    while (reflNum<INhkl and shellNum*RShell<RgLimit):
        shellNum=shellNum+1
        """Make a hkl"""
        """
        This whole set of loops needs to be
        changed into a far more efficient process.
        At the time of writing, the current method
        takes ~2.5 years to run on my laptop...
        """
        for Ih in range(-a,a-1):
            for Ik in range(-b,b-1):
                for Il in range(-c,c-1):
                    print(Ih,Ik,Il)
                    ISel=0
                    """Check that it's allowed by selection rules"""
                    ISel=SelectionRules(SSpaceGroupName,Ih,Ik,Il,ISel)
                    """
                    Still need to check that we have space in Rhkl
                    because the while doesn't kick in until we finish
                    the for loops
                    """
                    if (ISel==1) and (reflNum<INhkl):
                        """Make a g-vector"""
                        """Miller indices:"""
                        RGtest=np.array((Ih,Ik,Il),dtype="float")
                        """In microscope frame"""
                        RGtestM=(Ih*RarVecM)+(Ik*RbrVecM)+(Il*RcrVecM)
                        RGtestMag=np.sqrt(np.dot(RGtestM,RGtestM))
                        """Is it in the shell?"""
                        if (RGtestMag>(shellNum-1)*RShell) and (RGtestMag<=shellNum*RShell):
                            """Is it near a Laue condition |k+g|=|k|"""
                            RGplusk=RGtestM+Rk
                            """Divide by |g| to get a measure of
                            incident beam tilt"""
                            RDev=np.abs(RElectronWaveVectorMagnitude-np.sqrt(np.dot(RGplusk,RGplusk)))/RGtestMag
                            """
                            Tolerance of 0.08 here is rather
                            arbitrary --> might need revisiting
                            """
                            if ((RDev-0.08)<10**-9):
                                """
                                Add it to the pool and increment
                                the counter. Need to check that we
                                have space in Rhkl because the while
                                doesn't kick in until we finish the
                                do loops.
                                """
                                reflNum=reflNum+1
                                RhklTemp.append(RGtest)
    """Check to make sure we have enough beams"""
    if reflNum<=IMinStrongBeams:
        print("Beam pool is too small, please increase RgLimit. Quiting...")
        IErr=1
        """QUIT PROGRAM"""
    """Make the beam pool"""
    INhkl=reflNum
    print("Using a beam pool of",str(INhkl),"reflections")
    Rhkl=np.zeros((INhkl,3),dtype="float")
    Rhkl=RhklTemp
    
    """Now check that the required output hkls are in this hkl list"""
    for i in range(INoOfLacbedPatterns):
        flag=0
        for k in range(INhkl):
            if (RInputHKLs[i,0]-Rhkl[k,0]<10**-9) and (RInputHKLs[i,1]-Rhkl[k,1]<10**-9) and (RInputHKLs[i,2]-Rhkl[k,2]<10**-9):
                flag=1
        if flag==0:
            print("Input hkl not found:",np.rint(RInputHKLs[i]))
    
    
    return(INhkl,Rhkl)
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:01:17 2024

@author: Jacob Watkiss
"""

"""
-----------------------------------------------
Imports any necessary packages, libraries, modules etc.
-----------------------------------------------
"""

from CifFile import ReadCif
import numpy as np
import re

"""
A little function to streamline the process to read the file.
All it does is reads the value associated with the inputted
label.
"""
def find(file,target,filename,dtype):
    found=False
    for line in file:
        if target in line:
            found=True
            line=line.replace(target,"")
            line=clean(line,dtype)
            return line
    if found==False:
        IErr=1
        print(target,"does not exist in CIF")
        return "Absent"

"""
A function to clean up each value (set the correct data type
and remove any blank spaces).
"""
def clean(var,dtype):
    var=var.strip()
    if dtype=="str":
        var=var.replace("'","")
        var=var.replace(" ","")
        var=str(var)
    elif dtype=="int":
        var=int(var)
    elif dtype=="real":
        if re.search(r"\(.\)",var):
            var=re.sub(r"\(.\)","",var)
        elif re.search(r"\(..\)",var):
            var=re.sub(r"\(..\)","",var)
        var=float(var.strip())
    return var

"""
Defines the subprogram
"""
def read_cif(IErr):
    """
    -----------------------------------------------
    Declaring variables
    -----------------------------------------------
    """
    
    """Integers"""
    ISpaceGrp=0
    ILN=0
    IAtomCount=0
    IVolumeFLAG=0
    
    """Reals"""
    RLengthX=0.0
    RLengthY=0.0
    RLengthZ=0.0
    RAlpha=0.0
    RBeta=0.0
    RGamma=0.0 
    """These angles will be read in degrees,
    and converted to radians when needed for calculations"""
    RVolume=0.0
    
    """Strings"""
    SSpaceGroupName=""
    SSpaceGrp=""
    SChemicalFormula=""
    SSymString=""
    
    """
    The Basis arrays will be defined once IAtomCount
    is defined as their lengths are dependent on the
    amount of atoms.
    """
    
    """
    -----------------------------------------------
    Opening the file
    -----------------------------------------------
    """
    filename="felix.cif"
    name=open(filename,"r")
    content=name.readlines()
    
    """
    -----------------------------------------------
    Reading the data from the file
    -----------------------------------------------
    """
    
    """Space Group information"""
    ISpaceGrp=find(content,"_symmetry_Int_Tables_number",name,"int")
    SSpaceGrp=find(content,"_symmetry_space_group_name_H-M",name,"str")
    """To remove the quotes in the space group name:"""
    SSpaceGrp=SSpaceGrp.replace("'","")
    SSpaceGroupName=SSpaceGrp[0]
    
    """Cell dimensions"""
    RLengthX=find(content,"_cell_length_a",name,"real")
    RLengthY=find(content,"_cell_length_b",name,"real")
    RLengthZ=find(content,"_cell_length_c",name,"real")
    RAlpha=find(content,"_cell_angle_alpha",name,"real")
    RBeta=find(content,"_cell_angle_beta",name,"real")
    RGamma=find(content,"_cell_angle_gamma",name,"real")
    """Convert angles from degrees to radians"""
    RAlpha=RAlpha*np.pi/180
    RBeta=RBeta*np.pi/180
    RGamma=RGamma*np.pi/180
    RVolume=find(content,"_cell_volume",name,"real")

    if RVolume=="Absent":
        IVolumeFLAG=0
        RVolume1=RLengthX*RLengthY*RLengthZ
        RVolume2=1-np.cos(RAlpha)**2-np.cos(RBeta)**2-np.cos(RGamma)**2
        RVolume3=2*np.cos(RAlpha)*np.cos(RBeta)*np.cos(RGamma)
        RVolume=RVolume1*np.sqrt(RVolume2+RVolume3)
    else:
        IVolumeFLAG=1
    
    """Chemical Formula"""
    SChemicalFormula=find(content,"_chemical_formula_moiety",name,"str")
    if SChemicalFormula=="Absent":
        SChemicalFormula=find(content,"_chemical_formula_sum",name,"str")
    ILN=int(len(SChemicalFormula.strip()))
    
    """Atoms, their positions, Debye-Waller factors and occupancy"""
    
    """
    I have written a lot of this code myself, but could see
    no simple soution to overcoming the format of the atoms
    in the CIF. As a result, PyCifRW has been used where I
    can not find a solution.
    Hence the file is opened twice in two different ways.
    """
    cif_data=ReadCif(filename)
    block=list(cif_data.keys())[0]
    atoms=cif_data[block]
    
    SBasisAtomLabel=np.array(atoms.get("_atom_site_label",[]),dtype="str")
    IAtomCount=len(SBasisAtomLabel)
    SBasisAtomName=np.array(atoms.get("_atom_site_type_symbol",[]),dtype="str")
    
    """
    For the position, a dummy array must be created in order
    to make an array of the correct length with tuples as
    the elements. The CIF is then read, and tuples with each
    atom's position is created and then added to the correct array.
    """
    zeros=np.zeros(IAtomCount,dtype="float")
    RBasisAtomPosition=np.array(zeros,dtype="object")
    x_positions=atoms.get("_atom_site_fract_x",[])
    y_positions=atoms.get("_atom_site_fract_y",[])
    z_positions=atoms.get("_atom_site_fract_z",[])
    for i in range(IAtomCount):
        x=clean(x_positions[i],"real")
        y=clean(y_positions[i],"real")
        z=clean(z_positions[i],"real")
        position=(x,y,z)
        RBasisAtomPosition[i]=position
    
    RBasisIsoDW=np.zeros(IAtomCount,dtype="float")
    tempArray=np.array(atoms.get("_atom_site_B_iso_or_equiv",[]),dtype="str")
    DWconstant=True
    for i in range(IAtomCount):
        if tempArray.size==0:
            tempArray=np.array(atoms.get("_atom_site_U_iso_or_equiv",[]),dtype="str")
            if tempArray.size==0:
                IErr=1
                DWconstant=False
        if IErr==0:    
            RBasisIsoDW[i]=clean(tempArray[i],"real")*8*np.pi**2
    if DWconstant==False:
        print("Error: No Debye-Waller Factors for atoms --> RBasisIsoDw=0 for all atoms")
        RBasisIsoDW=np.zeros(IAtomCount,dtype="float")
        IErr=1
        
    RBasisOccupancy=np.array(atoms.get("_atom_site_occupancy",[]),dtype="float")
    if RBasisOccupancy.size==0:
        print("Error: No Occupancy for atoms --> RBasisOccupancy=0 for all atoms")
        RBasisOccupancy=np.zeros(IAtomCount,dtype="float")
        IErr=1
    
    """Symmetry Operation"""
    Stext=np.array(atoms.get("_symmetry_equiv_pos_as_xyz",[]),dtype="str")
    if Stext.size==0:
        Stext=np.array(atoms.get("_space_group_symop_operation_xyz",[]),dtype="str")
        if Stext.size==0:
            print("Error: No Symmetry Groups")    
    ISymCount=len(Stext)
    length=0
    for i in range(ISymCount):
        if len(Stext[i])>length:
            length=len(Stext[i])
    length=str(length)
    SSymString=np.empty((ISymCount,2),dtype=f"<U"+length)
    for i in range(ISymCount):
        SSymString[i][0]=str(i+1)
        SSymString[i][1]=(Stext[i])
    
    """
    -----------------------------------------------
    Closing the file
    -----------------------------------------------
    """
    name.close()
    
    """
    -----------------------------------------------
    Return everything -- can be changed,
    but I wanted to avoid global variables
    -----------------------------------------------
    """
    return (ISpaceGrp,IAtomCount,ILN,IVolumeFLAG,
    SBasisAtomLabel,SBasisAtomName,SSpaceGroupName,SSpaceGrp,SChemicalFormula,SSymString,
    RLengthX,RLengthY,RLengthZ,RAlpha,RBeta,RGamma,RVolume,RBasisAtomPosition,RBasisIsoDW,RBasisOccupancy)
    
read_cif(0)
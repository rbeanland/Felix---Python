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
import globalvariables as g
IErr=g.IErr

"""
A little function to streamline the process to read the file.
All it does is reads the value associated with the inputted
label. --> Used in other file-reading programs for efficiency
"""
def find(file,target,dtype):
    found=False
    for line in file:
        if target in line:
            found=True
            line=line.replace(target,"")
            line=clean(line,dtype)
            return line
    if found==False:
        return "Absent"

"""
A function to clean up each value (set the correct data type
and remove any blank spaces). Removes: spaces, ', = and values
in brackets - the uncertainties - so that the values can be 
stored as floats or integers.
"""
def clean(var,dtype):
    var=var.strip()
    var=var.replace("=","")
    if dtype=="str":
        var=var.replace("'","")
        var=var.replace(" ","")
        var=str(var)
    elif dtype=="int":
        var=re.sub(r"\(.\)","",var)
        var=re.sub(r"\(..\)","",var)
        var=int(var)
    elif dtype=="real":
        var=re.sub(r"\(.\)","",var)
        var=re.sub(r"\(..\)","",var)
        var=float(var.strip())
    return var

"""
A function that turns "x/y" into a float value of the result
"""
def divide(var):
    var_str=str(var)
    div=re.search(r"/",var_str)
    if div != None:
        numerator=var[0:div.start()]
        denominator=var[div.end():len(var)]
        var=float(int(numerator)/int(denominator))
    return var

"""
Turns "ax+b" into two variables: a and b
to be used when creating tensors for multiplication
from strings 
"""
def split(var):
    plus=re.search(r"\+",var)
    if plus==None:
        minus=re.search(r"\-",var)
        if minus==None:
            constant=0
            coefficient=var
        else:
            constant=var[0:minus.start()]
            coefficient=var[minus.start():len(var)]
            if minus.start()==0:
                constant=0    
    else:
        constant=var[0:plus.start()]
        coefficient=var[plus.end():len(var)]
    constant=divide(constant)
    coefficient=coefficient.replace("x","")
    coefficient=coefficient.replace("y","")
    coefficient=coefficient.replace("z","")
    if coefficient=="-":
        coefficient=-1
    elif coefficient=="":
        coefficient=1
    coefficient=divide(coefficient)
    return constant,coefficient
    
r"""
A function to turn "(ax+d,by+e,cz+f)" into a matrix
and a vector that can be applied to a position (x,y,z)
and produce the same result. i.e. (ax+d,by+e,cz+f)
becomes:  /a 0 0  \        /d  \
         | 0 b 0  |  and  | e  |
         \ 0 0 c /        \ f /
"""
def convertToMatrix(SSymString):
    ISymCount=len(SSymString)
    mat=np.zeros((ISymCount,3,3),dtype="float")
    vec=np.zeros((ISymCount,3),dtype="float")
    for i in range(ISymCount):
        symOp=SSymString[i][1]
        firstComma=re.search(",",symOp).start()
        xComponent=symOp[0:firstComma].strip()
        symOp=symOp[firstComma+1:len(symOp)]
        secondComma=re.search(",",symOp).start()
        yComponent=symOp[0:secondComma].strip()
        symOp=symOp[secondComma+1:len(symOp)]
        zComponent=symOp.strip()
        
        xConstant,xCoefficient=split(xComponent)
        yConstant,yCoefficient=split(yComponent)
        zConstant,zCoefficient=split(zComponent)
        
        row1=np.array((xCoefficient,0,0),dtype="float")
        row2=np.array((0,yCoefficient,0),dtype="float")
        row3=np.array((0,0,zCoefficient),dtype="float")
        
        mat[i][0]=row1
        mat[i][1]=row2
        mat[i][2]=row3
        vec[i][0]=xConstant
        vec[i][1]=yConstant
        vec[i][2]=zConstant
        
    return mat,vec



"""
-----------------------------------------------
Defines the subprogram
-----------------------------------------------
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
    ISpaceGrp=find(content,"_symmetry_Int_Tables_number","int")
    SSpaceGrp=find(content,"_symmetry_space_group_name_H-M","str")
    """To remove the quotes in the space group name:"""
    SSpaceGrp=SSpaceGrp.replace("'","")
    SSpaceGroupName=SSpaceGrp[0]
    
    """Cell dimensions"""
    RLengthX=find(content,"_cell_length_a","real")
    RLengthY=find(content,"_cell_length_b","real")
    RLengthZ=find(content,"_cell_length_c","real")
    RAlpha=find(content,"_cell_angle_alpha","real")
    RBeta=find(content,"_cell_angle_beta","real")
    RGamma=find(content,"_cell_angle_gamma","real")
    """Convert angles from degrees to radians"""
    RAlpha=RAlpha*np.pi/180
    RBeta=RBeta*np.pi/180
    RGamma=RGamma*np.pi/180
    RVolume=find(content,"_cell_volume","real")

    if RVolume=="Absent":
        print("No Cell Volume found: Calculating from cell dimensions")
        IVolumeFLAG=0
        RVolume1=RLengthX*RLengthY*RLengthZ
        RVolume2=1-g.cos(RAlpha)**2-g.cos(RBeta)**2-g.cos(RGamma)**2
        RVolume3=2*g.cos(RAlpha)*g.cos(RBeta)*g.cos(RGamma)
        RVolume=RVolume1*np.sqrt(RVolume2+RVolume3)
    else:
        IVolumeFLAG=1
    
    """Chemical Formula"""
    SChemicalFormula=find(content,"_chemical_formula_moiety","str")
    if SChemicalFormula=="Absent":
        SChemicalFormula=find(content,"_chemical_formula_sum","str")
        if SChemicalFormula=="Absent":
            SChemicalFormula=find(content,"_chemical_formula_iupac","str")
            if SChemicalFormula=="Absent":
                SChemicalFormula=find(content,"_chemical_formula_structural","str")
                if SChemicalFormula=="Absent":
                    print("No Chemical Formula found: Check CIF")
                    IErr=1
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
    Find atomic number for atoms
    """
    elementalSymbols=np.array(("H","He",
                       "Li","Be","B","C","N","O","F","Ne",
                       "Na","Mg","Al","Si","P","S","Cl","Ar",
                       "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                       "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                       "Cs","Ba",
                       "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                       "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
                       "Fr","Ra",
                       "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
                       "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"),dtype="str")
    
    IBasisAtomicNumber=np.zeros(IAtomCount,dtype="int")
    for x in range(IAtomCount):
        for i in range(len(elementalSymbols)):
            if SBasisAtomName[x]==elementalSymbols[i]:
                IBasisAtomicNumber[x]=i+1
    
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
        print("No Debye-Waller Factors found: RBasisIsoDw=0 for all atoms")
        RBasisIsoDW=np.zeros(IAtomCount,dtype="float")
        IErr=1
        
    RBasisOccupancy=np.array(atoms.get("_atom_site_occupancy",[]),dtype="float")
    if RBasisOccupancy.size==0:
        print("No Occupancy found: RBasisOccupancy=1 for all atoms")
        RBasisOccupancy=np.ones(IAtomCount,dtype="float")
        IErr=1
    
    """Symmetry Operation"""
    Stext=np.array(atoms.get("_symmetry_equiv_pos_as_xyz",[]),dtype="str")
    if Stext.size==0:
        Stext=np.array(atoms.get("_space_group_symop_operation_xyz",[]),dtype="str")
        if Stext.size==0:
            print("No Symmetry Groups found")    
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

    """Creating RSymMat"""
    """
    RSymMat turns the symmetry operations into
    usable matrices
    """
    RSymMat=np.zeros((ISymCount,3,3),dtype="float")
    RSymVec=np.zeros((ISymCount,3),dtype="float")
    RSymMat,RSymVec=convertToMatrix(SSymString)
    
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
    RLengthX,RLengthY,RLengthZ,RAlpha,RBeta,RGamma,RVolume,RBasisAtomPosition,RBasisIsoDW,RBasisOccupancy,
    IBasisAtomicNumber,RSymMat,RSymVec)
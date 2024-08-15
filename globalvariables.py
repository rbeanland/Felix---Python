# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 15:21:27 2024

@author: Jacob Watkiss
"""

import numpy as np

"""
I have made this against my best nature,
I wish there were 0 global variables in 
the code, but I think that certain things
(such as IErr, the error checker) and 
universal constants like c, e and h should
be referable so that they don't have to be
declared every time.
I am aware that numpy is good for some constants
and scipy can do any others, but I am merely
translating this code, so I have just copied
the values for these constants from the original
felix.
All values in here should be constants except for,
of course, IErr.
- Jacob
"""

IErr=0
TINY=10**-9
HUGE=10**9
RSpeedOfLight=float(2.99762458*10**8)
RElectronMass=float(9.10938291*10**-31)
RElectronMassMeV=float(0.510998928)
RPlanckConstant=float(6.62606957*10**-34)
RElectronCharge=float(1.602176565*10**-19)
RAngstromConversion=float(10**10)

def cos(angle):
    if angle==0:
        result=1
    elif angle==np.pi/6:
        result=np.sqrt(3)/2
    elif angle==np.pi/4:
        result=1/np.sqrt(2)
    elif angle==np.pi/3:
        result=1/2
    elif angle==np.pi/2:
        result=0
    else:
        result=np.cos(angle)
    return result

def sin(angle):
    if angle==0:
        result=0
    elif angle==np.pi/6:
        result=1/2
    elif angle==np.pi/4:
        result=1/np.sqrt(2)
    elif angle==np.pi/3:
        result=np.sqrt(3)/2
    elif angle==np.pi/2:
        result=1
    else:
        result=np.sin(angle)
    return result





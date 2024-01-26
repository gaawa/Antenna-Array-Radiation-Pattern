# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 12:27:54 2022

@author: ThinkWin
"""
from enum import Enum, auto

class ArrayFormation(Enum):
    Linear = auto()
    Circular = auto()
    LinearCircular = auto()
    BowCircular = auto()
    
class PlotMode(Enum):
    plot2d = auto()
    plot3d = auto()
    
class BeamFormer(Enum):
    Projection = auto()
### LoopSheafFigure.py
### MIT LICENSE 2024 Alex Dowling

from .SheafFigure import *
from .ParameterLoopFigure import *

class LoopSheafFigure(SheafFigure, ParameterLoopFigure):
    def __init__(self, shf, stg_dict, parameter_graph=None, start=None):
        self.shf = shf
        ParameterLoopFigure.__init__(self, shf.P.children_, stg_dict, 
                                     parameter_graph=parameter_graph, 
                                     start=start)
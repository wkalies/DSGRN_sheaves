### PathSheafFigure.py
### MIT LICENSE 2024 Alex Dowling

from .SheafFigure import *
from .ParameterPathFigure import *

class PathSheafFigure(SheafFigure, ParameterPathFigure):
    def __init__(self, shf, stg_dict, parameter_graph=None):
        self.shf = shf
        ParameterPathFigure.__init__(self, shf.P.children_, stg_dict, 
                                     parameter_graph=parameter_graph)

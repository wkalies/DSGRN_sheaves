import DSGRN_utils
from .BifurcationQuery import *

class HysteresisQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())

        # ** ONE BISTABLE PARAMETER NODE **
        vertices = ['a','b','c']
        edges = [('a','b'), ('b','c')]
        match_grading = {1 : ['a','c'], 2 : ['b']}

        # ** TWO BISTABLE PARAMETER NODES **
        # vertices = ['a','b','c', 'd']
        # edges = [('a','b'), ('b','c'), ('c', 'd')]
        # match_grading = {1 : ['a','d'], 2 : ['b', 'c']}
        
        coho_criteria = [
                         # {'selection' : ['a','b'],
                         #  'predicate' : lambda sc : len(sc[0]) == 2, # SN                       
                         #  'dim' : 1,
                         #  'clean_stalks' : True
                         # }, 
                         # {'selection' : ['b','c'],
                         #  'predicate' : lambda sc : len(sc[0]) == 2, # SN                       
                         #  'dim' : 1,
                         #  'clean_stalks' : True
                         # }, 
                         {
                          'predicate' : lambda sc : len(sc[0]) == 1, # hysteresis
                          #'predicate' : lambda sc : len(sc[0]) == 2, # isola
                          'dim' : 1,
                          'clean_stalks' : True
                         }
                        ]    
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)

import DSGRN_utils
from .BifurcationQuery import *

class CuspQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())
    
        vertices = ['a','b','c','d']
        edges = [('a','b'), ('b','c'), ('c','d'), ('a','d')]
        match_grading = {1 : ['a','b','c'], 2 : ['d']}
        coho_criteria = [{
                          'predicate' : lambda sc : len(sc[0]) == 1,
                          'dim' : 1,
                          'clean_stalks' : True
                         }]
        
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)
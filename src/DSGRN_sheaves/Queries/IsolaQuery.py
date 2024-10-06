import DSGRN_utils
from .BifurcationQuery import *

class IsolaQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None,SN_check=True,ppath=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())
        
        vertices = ['a','b','c']
        edges = [('a','b'), ('b','c')]
        match_grading = {1 : ['a','c'], 2 : ['b']}

        if SN_check:
            coho_criteria = [
                             {'selection' : ['a','b'],
                              'predicate' : lambda sc : len(sc[0]) == 2, # SN                       
                              'dim' : 1,
                              'clean_stalks' : True
                             }, 
                             {'selection' : ['b','c'],
                              'predicate' : lambda sc : len(sc[0]) == 2, # SN                       
                              'dim' : 1,
                              'clean_stalks' : True
                             }, 
                             {
                              'predicate' : lambda sc : len(sc[0]) == 2,
                              'dim' : 1,
                              'clean_stalks' : True
                             }
                            ]    
        else:
            coho_criteria = [{
                              'predicate' : lambda sc : len(sc[0]) == 2,
                              'dim' : 1,
                              'clean_stalks' : True
                             }]    
        
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria, ppath=ppath)

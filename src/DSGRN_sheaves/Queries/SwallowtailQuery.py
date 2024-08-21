import DSGRN_utils
from .BifurcationQuery import *

class SwallowtailQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())
    
        vertices = ['a','b','c','d','e','f','g','h']
        edges = [('a','b'),('a','c'),('c','d'),('b','d'),
                 ('e','f'),('e','g'),('g','h'),('f','h'),
                 ('a','e'),('b','f'),('c','g'),('d','h')]
        match_grading = {1 : ['b'], 2 : ['a','c','d','e','f','g'], 3 : ['h']}
        coho_criteria = [
                         {'selection' : ['b','d','f','h'],
                          'predicate' : lambda sc : len(sc[0]) == 3 
                                                    and len(sc[1]) == 3, #bisolas
                          'dim' : 1,
                          'clean_stalks' : True}, 
                         {'selection' : ['a','c','e','g'], 
                          'predicate' : lambda sc : len(sc[0]) == 3 
                                                    and len(sc[1]) == 3, #stable
                          'dim' : 1,
                          'clean_stalks' : True}, 
                         {'selection' : ['c','d','g','h'], 
                          'predicate' : lambda sc : len(sc[0]) == 3 
                                                    and len(sc[1]) == 3, #cusp
                          'dim' : 1,
                          'clean_stalks' : True}, 
                         {'selection' : ['e','f','g','h'],
                          'predicate' : lambda sc : len(sc[0]) == 3 
                                                    and len(sc[1]) == 3, #cusp
                          'dim' : 1,
                          'clean_stalks' : True}, 
                         {'selection' : ['a','b','c','d'],
                          'predicate' : lambda sc : len(sc[0]) == 2 
                                                    and len(sc[1]) == 2, #isola
                          'dim' : 1,
                          'clean_stalks' : True}, 
                         {'selection' : ['a','b','e','f'], 
                          'predicate' : lambda sc : len(sc[0]) == 2 
                                                    and len(sc[1]) == 2, #isola
                          'dim' : 1,
                          'clean_stalks' : True}                 
                        ]

        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)
        


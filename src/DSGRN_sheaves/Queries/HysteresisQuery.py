import DSGRN_utils
from .BifurcationQuery import *

class HysteresisQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None,SN_check=True,num_bistable_nodes=1):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())

        # ** ONE BISTABLE PARAMETER NODE **
        # vertices = ['a','b','c']
        # edges = [('a','b'), ('b','c')]
        # match_grading = {1 : ['a','c'], 2 : ['b']}

        # ** TWO BISTABLE PARAMETER NODES **
        # vertices = ['a','b','c', 'd']
        # edges = [('a','b'), ('b','c'), ('c', 'd')]
        # match_grading = {1 : ['a','d'], 2 : ['b', 'c']}

        vertices = ['a','b1']
        edges = [('a','b1')]
        for i in range(2,num_bistable_nodes+1):
            vertices.append('b'+str(i))
            edges.append(('b'+str(i-1),'b'+str(i)))
        vertices.append('c')
        edges.append(('b'+str(num_bistable_nodes),'c'))
        bistable_vertices=vertices[1:num_bistable_nodes+1]
        match_grading = {1 : ['a','c'], 2 : bistable_vertices}
        selection_SN_left=['a','b1']
        selection_SN_right=[bistable_vertices[-1],'c']

        if SN_check:
            coho_criteria = [
                             {'selection' : ['a','b1'], 
                              'predicate' : lambda sc : len(sc[0]) == 2, # SN                       
                              'dim' : 1,
                              'clean_stalks' : True
                             }, 
                             # {'selection' : ['b'+str(num_bistable_nodes),'c'], 
                             #  'predicate' : lambda sc : len(sc[0]) == 2, # SN                       
                             #  'dim' : 1,
                             #  'clean_stalks' : True
                             # }, 
                             {
                              'predicate' : lambda sc : len(sc[0]) == 1, # hysteresis
                              'dim' : 1,
                              'clean_stalks' : True
                             }
                            ]   
        else:
            coho_criteria = [{
                              'predicate' : lambda sc : len(sc[0]) == 1, # hysteresis
                              'dim' : 1,
                              'clean_stalks' : True
                             }]   
            
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)

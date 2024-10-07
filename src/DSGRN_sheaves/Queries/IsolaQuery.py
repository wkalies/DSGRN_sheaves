import DSGRN_utils
from .BifurcationQuery import *

class IsolaQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None, SN_check=True, num_bistable_nodes=1, SN_SC=[2,2], ppath=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())
        
        # vertices = ['a','b','c']
        # edges = [('a','b'), ('b','c')]
        # match_grading = {1 : ['a','c'], 2 : ['b']}

        vertices = ['a','b1']
        edges = [('a','b1')]
        for i in range(2,num_bistable_nodes+1):
            vertices.append('b'+str(i))
            edges.append(('b'+str(i-1),'b'+str(i)))
        vertices.append('c')
        edges.append(('b'+str(num_bistable_nodes),'c'))
        bistable_vertices=vertices[1:num_bistable_nodes+1]
        match_grading = {1 : ['a','c'], 2 : bistable_vertices}
        
        if SN_check:
            coho_criteria = [
                             {'selection' : ['a','b1'],
                              'predicate' : lambda sc : len(sc[0]) == SN_SC[0], # 2 SN                       
                              'dim' : 1,
                              'clean_stalks' : True
                             }, 
                             {'selection' : ['b'+str(num_bistable_nodes),'c'],
                              'predicate' : lambda sc : len(sc[0]) == SN_SC[1], # 2 SN                       
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

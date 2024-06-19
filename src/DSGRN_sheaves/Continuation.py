### Continuation.py
### MIT LICENSE 2024 Alex Dowling

import DSGRN
import DSGRN_utils
import pychomp
import itertools

from DSGRN_sheaves.Sheaf import *
from DSGRN_sheaves.Cohomology import *
from DSGRN_sheaves.CechCell import *

class Interrupt(Exception):
    pass
    
class EqLandmine:
    def __eq__(self, other):
        raise Interrupt

def blank_stg(parameter_graph):
    parameter = parameter_graph.parameter(0)
    try:
        stg = DSGRN_utils.CubicalBlowupGraph.__new__(DSGRN_utils.CubicalBlowupGraph)
        stg.__init__(parameter, level=EqLandmine())
    except Interrupt:
        pass
    
    stg.digraph.graph_vertices = set(stg.blowup_complex(stg.dim))
    stg.digraph.adjacency_lists = {cell:set() for cell in stg.blowup_complex(stg.dim)}
    return stg

def stg_union(stg_list, level = 8):
    """ Inputs a list of RookRulesCubicalComplex objects with the same 
        sets of vertices. Outputs a RookRulesCubicalComplex whose digraph is
        the of those from the complexes of given list."""

    total_stg = blank_stg(DSGRN.ParameterGraph(stg_list[0].network))
    for v in total_stg.digraph.vertices():
        total_stg.digraph.adjacency_lists[v].update(*[stg.digraph.adjacency_lists[v] 
                                                      for stg in stg_list])
    return total_stg

def build_parameter_complex(parameter_graph, parameter_indices = None, 
                           dim = None, length_cap = 2, level = 3):
    """ Inputs a DSGRN parameter graph, a collection of parameter indices, and
        a maximum dimension. Outputs a directed acyclic graph of keys, and
        a dictionary of RookRulesCubicalComplex objects."""
    
    param_dim = parameter_graph.dimension()
    if parameter_indices is None: 
        parameter_indices = set(range(parameter_graph.size()))
    if dim is None: 
        dim = param_dim
    top_cells = {}
    parameter_complex = pychomp.DirectedAcyclicGraph()
    stg_dict = {}

    # Build the top dimensional cells and their state transition graphs first.
    # Populate a dictionary top_cells so later we won't have to build the 
    # CechCell again.
    for index in parameter_indices:
        parameter = parameter_graph.parameter(index)
        cell = top_cech_cell(parameter_graph, index, dim)
        top_cells.update({index : cell})
        parameter_complex.add_vertex(cell)
        stg_dict.update({cell : DSGRN_utils.CubicalBlowupGraph(
                                parameter, level)})

    # Build the remaining cells at every index from the top down (BFS ordering). 
    # If there aren't enough adjacent indices to coarsen with, stick in some 
    # copies of the top index.
    for codim, index in itertools.product(range(1, dim+1), parameter_indices):
        top_cell = top_cells[index]
        adj_indices = list(set(parameter_indices).intersection(
                           parameter_graph.adjacencies(index, "codim1")))
        adj_indices = adj_indices + max(0, dim - len(adj_indices) + 1)*[index]

        for combo in itertools.combinations(adj_indices, codim):
            adj_cells = [top_cells[index] for index in combo]
            cells = adj_cells + [top_cell]
            if (len(combo) >= length_cap and 
                any(len({cell[d] for cell in cells}) > length_cap
                    for d in range(param_dim))):
                continue
                
            cell_join = top_cell.join(adj_cells, dim - codim)
            cell = cell_join.permute()
            parameter_complex.add_vertex(cell)
            
            cell_stgs = [stg_dict[cell] for cell in cells]
            if cell in stg_dict:
                cell_stgs.append(stg_dict[cell])
            stg_dict.update({cell : stg_union(cell_stgs)})

            child_cells = set(cell_join.get_child_cells())
            child_cells &= parameter_complex.vertices()
            for child_cell in child_cells:
                parameter_complex.add_edge(cell, child_cell.permute())
            expand_cell = CechCell(cell.inequality_sets, cell.dim + 1, 
                                   cell.labels)
            if expand_cell in stg_dict:
                parameter_complex.add_edge(cell, expand_cell)

    return parameter_complex, stg_dict

def morse_restriction(stalk_0, stalk_1):
    m = len(stalk_1)
    n = len(stalk_0)
    
    R = [[0 for j in range(n)] for i in range(m)]
    for j, morse_0 in enumerate(stalk_0):      
        for i, morse_1 in enumerate(stalk_1):
            if morse_1.issubset(morse_0):
                R[i][j] = 1
    return R

def attractor_sheaf(parameter_complex, stg_dict, prune_grad='none'):
    """ Inputs a directed acyclic graph of keys and a dictionary of 
        RookRulesCubicalComplex objects with the same set of vertices. Outputs 
        the attractor sheaf on the corresponding poset."""
 
    dim = max([cell.dim for cell in parameter_complex.vertices()])
    stalk_dict = {}
    restr_dict = {}
    grading = [set() for i in range(dim+1)]
    
    for cell, stg in stg_dict.items():
        (scc_dag, graded_complex) = pychomp.FlowGradedComplex(stg.complex(), 
                                                              stg.digraph.adjacencies)
        connection_matrix = pychomp.ConnectionMatrix(graded_complex)
        if prune_grad == 'all' or (prune_grad == 'some' and cell.dim == dim):
            morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_complex, 
                                                 connection_matrix)
        else:
            morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_complex, 
                                                 connection_matrix, False)
        stalk_dict.update({cell : [frozenset({c for c in stg.digraph.vertices() 
                                             if graded_complex.value(c) == v}) 
                                  for v in morse_graph.vertices()]})
        grading[cell.dim].add(cell)

    for (cell_0, cell_1) in parameter_complex.edges():
        R = morse_restriction(stalk_dict[cell_0], stalk_dict[cell_1])
        restr_dict.update({(cell_0, cell_1) : R})
         
    return Sheaf(pychomp.Poset(parameter_complex), grading, stalk_dict, 
                 restr_dict)

import DSGRN
import DSGRN_utils
import math
import numpy as np
from numpy.linalg import matrix_rank
import time
import sys

from DSGRN_sheaves.Sheaf import *
from DSGRN_sheaves.Cohomology import *
from DSGRN_sheaves.Continuation import *
from DSGRN_sheaves.Attractors import *
from DSGRN_sheaves.SearchBifurcations import *
from DSGRN_sheaves.PlotAttractorSheaf import *
from DSGRN_sheaves.CechCell import *

class SaddleNodeQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
        
        vertices = ['a','b']
        edges = [('a','b')]
        match_grading = {1 : ['a'], 2 : ['b']}
        coho_criteria = [{
                          'predicate' : lambda sc : len(sc[0]) == 2,
                          'dim' : 1,
                          'clean_stalks' : True
                         }]
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)

class PitchforkQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
        
        vertices = ['a','b']
        edges = [('a','b')]
        match_grading = {1 : ['a'], 2 : ['b']}
        coho_criteria = [{
                          'predicate' : lambda sc : len(sc[0]) == 1,
                          'dim' : 1,
                          'clean_stalks' : True
                         }]
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)


class HysteresisQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
        
        vertices = ['a','b','c']
        edges = [('a','b'), ('b','c')]
        match_grading = {1 : ['a','c'], 2 : ['b']}
        coho_criteria = [{
                          'predicate' : lambda sc : len(sc[0]) == 1,
                          'dim' : 1,
                          'clean_stalks' : True
                         }]    
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)

class IsolaQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
        
        vertices = ['a','b','c']
        edges = [('a','b'), ('b','c')]
        match_grading = {1 : ['a','c'], 2 : ['b']}
        coho_criteria = [{
                          'predicate' : lambda sc : len(sc[0]) == 2,
                          'dim' : 1,
                          'clean_stalks' : True
                         }]    
        
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)

class CuspQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
    
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

class IsolaLoopQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
    
        vertices = ['a','b','c','d']
        edges = [('a','b'), ('b','c'), ('c','d'), ('a','d')]
        match_grading = {1 : ['a','b','c'], 2 : ['d']}
        coho_criteria = [{
                          'predicate' : lambda sc : len(sc[0]) == 2,
                          'dim' : 1,
                          'clean_stalks' : True
                         }]
        
        super().__init__(parameter_graph, vertices, edges, 
                         param_stability, match_grading, coho_criteria)

class SwallowtailQuery(BifurcationQuery):

    def __init__(self, parameter_graph, param_stability=None):
        if param_stability is None:
            param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network)
    
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

def general_hysteresis_query(parameter_graph, length, param_stability=None):
    
    if length < 3:
        return 

    if param_stability is None:
        param_stability = DSGRN_utils.StabilityQuery(parameter_graph.network())
    param_grading = {-1 : [i for i in range(parameter_graph.size())]}
    param_grading.update({1 : []})
    for key in param_stability.keys():
        if key > 1:
            param_grading[1] = param_grading[1] + param_stability[key]

    vertices = list(range(length))
    edges = [(i, i+1) for i in vertices[:-1]]
    match_grading = {-1 : [0, length-1], 1 : vertices[1:-1]}
    
    def in_img(M, v):
        A = np.concatenate((M, v), axis=1)
        return matrix_rank(M) == matrix_rank(A)

    def hysteresis_criteria(pg, match, ordering):
        #function for checking 
        l_match = match[:-1]
        r_match = match[:-2] + [match[-1]]
        c_match  = match[:-2]

        l_edge_key = top_cech_cell(parameter_graph, match[-2], 1)
        r_edge_key = top_cech_cell(parameter_graph, match[-1], 1)
        dummy_cell = CechCell(tuple(frozenset({('dummy',)})), 0, l_edge_key.labels)

        def sheaf_data(indices):
            pc, stg_dict = build_parameter_complex(parameter_graph, indices, 1)
            if len(indices) < 2:
                top_cell = top_cech_cell(parameter_graph, indices[0], 1)
                pc.add_edge(dummy_cell, top_cell)
                stg_dict.update({dummy_cell : stg_dict[top_cell]})
            shf = attractor_sheaf(pc, stg_dict)
            shf_cohomology = sheaf_cohomology(shf)
            rank = sum([len(shf.stalk(key)) for key in shf.grading[0]])
            return pc, stg_dict, shf, shf_cohomology, rank

        pc, stg_dict, shf, shf_cohomology, rank = sheaf_data(match)
        morse_dict = morse_dictionary(pc, stg_dict)
        att_secs = attractor_sections(shf, morse_dict)
        
        c_pc, c_stg_dict, c_shf, c_shf_cohomology, c_rank = sheaf_data(c_match, True)
        c_morse_dict = morse_dictionary(c_pc, c_stg_dict)
        c_att_secs = attractor_sections(c_shf, c_morse_dict)

        l_pc, l_stg_dict, l_shf, l_shf_cohomology, l_rank = sheaf_data(l_match)
        r_pc, r_stg_dict, r_shf, r_shf_cohomology, r_rank = sheaf_data(r_match)
        
        R_tc = shf.GF([[0 for j in range(rank)] for i in range(c_rank)])
        c_ranges = {}
        row = 0
        for key in c_shf.grading[0]:
            c_ranges.update({key : (row, row+len(c_shf.stalk(key)))})
            row = row + len(c_shf.stalk(key))
        col = 0
        for key in shf.grading[0]:
            if len(shf.P.children(key)) < 2:
                pass
            elif l_edge_key in shf.P.children(key):
                c_edge_key = [k for k in shf.P.children(key) 
                              if k != l_edge_key][0]
                target_key = CechCell(c_edge_key.inequality_sets, 0, c_edge_key.labels)
                R = morse_restriction(shf.stalk(key), c_shf.stalk(target_key))
                R_tc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(shf.stalk(key))] = R
            elif r_edge_key in shf.P.children(key):
                c_edge_key = [k for k in shf.P.children(key) 
                              if k != r_edge_key][0]
                target_key = CechCell(c_edge_key.inequality_sets, 0, c_edge_key.labels)
                if length == 3:
                    target_key = dummy_cell
                R = morse_restriction(shf.stalk(key), c_shf.stalk(target_key))
                R_tc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(shf.stalk(key))] = R
            else:
                target_key = key
                R = shf.GF(np.eye(len(shf.stalk(key))).astype(int))
                R_tc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(shf.stalk(key))] = R
            col = col + len(shf.stalk(key))

        R_lc = l_shf.GF([[0 for j in range(l_rank)] for i in range(c_rank)])
        col = 0
        for key in l_shf.grading[0]:
            if (l_edge_key in l_shf.P.children(key) 
                and len(l_shf.P.children(key)) > 1):
                c_edge_key = [k for k in l_shf.P.children(key) 
                              if k != l_edge_key][0]
                target_key = CechCell(c_edge_key.inequality_sets, 0, c_edge_key.labels)
                R = morse_restriction(l_shf.stalk(key), c_shf.stalk(target_key))
                R_lc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(l_shf.stalk(key))] = R
            elif l_edge_key not in l_shf.P.children(key):
                target_key = key
                if length == 3:
                    target_key = dummy_cell
                R = morse_restriction(l_shf.stalk(key), c_shf.stalk(target_key))
                #R = l_shf.GF(np.eye(len(l_shf.stalk(key))).astype(int))
                R_lc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(l_shf.stalk(key))] = R
            col = col + len(l_shf.stalk(key))

        R_rc = r_shf.GF([[0 for j in range(r_rank)] for i in range(c_rank)])
        col = 0
        for key in r_shf.grading[0]:
            if (r_edge_key in r_shf.P.children(key) 
                and len(r_shf.P.children(key))) > 1:
                c_edge_key = [k for k in r_shf.P.children(key) 
                              if k != r_edge_key][0]
                target_key = CechCell(c_edge_key.inequality_sets, 0, c_edge_key.labels)
                if length == 3:
                    target_key = dummy_cell
                R = morse_restriction(r_shf.stalk(key), c_shf.stalk(target_key))
                R_rc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(r_shf.stalk(key))] = R
            elif r_edge_key not in r_shf.P.children(key):
                target_key = key
                R = r_shf.GF(np.eye(len(r_shf.stalk(key))).astype(int))
                R_rc[c_ranges[target_key][0]:c_ranges[target_key][1], 
                     col:col+len(r_shf.stalk(key))] = R
            col = col + len(r_shf.stalk(key))
      
        check = False
        for section in att_secs:
            if len(att_secs.children(section)) != 1:
                continue
            zero = list(att_secs.children(section))[0]
            if any(a!=0 for a in zero):
                continue
            c_section = tuple([int(s==1) 
                               for s in np.matmul(R_tc, shf.GF(section))])
            if len(c_att_secs.children(c_section)) != 1:
                continue
            pred = list(c_att_secs.children(c_section))[0]
            if len(c_att_secs.children(pred)) != 2:
                continue

            s0 = c_shf.GF([[a] for a in list(c_att_secs.children(pred))[0]])
            s1 = c_shf.GF([[a] for a in list(c_att_secs.children(pred))[1]])
            K_l = l_shf.GF(l_shf_cohomology[0]).transpose()
            K_r = r_shf.GF(r_shf_cohomology[0]).transpose()       
            M_lc = np.matmul(R_lc, K_l)
            M_rc = np.matmul(R_rc, K_r)
            
            hys01 = (in_img(M_lc, s0) and in_img(M_rc, s1) 
                     and not in_img(M_lc, s1) and not in_img(M_rc, s0))
            hys10 = (in_img(M_lc, s1) and in_img(M_rc, s0) 
                     and not in_img(M_lc, s0) and not in_img(M_rc, s1))
            if hys01 or hys10:
                check = True
                break

        return check

    coho_criteria = [{"custom" : hysteresis_criteria}]

    return BifurcationQuery(parameter_graph, vertices, edges, param_grading,
                            match_grading, coho_criteria)
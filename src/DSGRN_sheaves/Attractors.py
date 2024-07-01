### Attractors.py
### MIT LICENSE 2024 Alex Dowling
###############################################################################
import itertools
import numpy as np
import pychomp
import DSGRN_utils

from .Sheaf import *

def morse_dictionary(parameter_complex, stg_dict, prune_grad='none'):
    dim = max(cell.dim for cell in parameter_complex.vertices())
    morse_dict = {}
    for cell in parameter_complex.vertices():
        stg = stg_dict[cell]
        (scc_dag, graded_comp) = pychomp.FlowGradedComplex(stg.complex(), 
                                                           stg.adjacencies())
        connection_matrix = pychomp.ConnectionMatrix(graded_comp)
        if prune_grad == 'all' or (prune_grad == 'some' and cell.dim == dim):
            morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_comp, 
                                                 connection_matrix)
        else:
            morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_comp, 
                                                 connection_matrix, False)
        morse_dict.update({cell : morse_graph})
    return morse_dict

def morse_sets_in_section(shf, morse_dict, section):
    n = 0
    section_morse_sets = {}
    for key in shf.grading[0]:
        m = len(shf.stalk(key))
        stalk_list = [v for v in morse_dict[key].vertices()]
        stalk_sets = [stalk_list[i] for i in range(m) if section[n+i] != 0]
        section_morse_sets.update({key : stalk_sets})
        n = n+m
    return section_morse_sets

def section_from_morse_sets(shf, section_morse_sets):
    section = []
    for key in shf.grading[0]:
        section = section + [int(M in section_morse_sets[key]) 
                             for M in shf.stalk(key)]
    return shf.GF(section)

def attractor_sections(shf, morse_dict, shf_cohomology=None):
    if shf_cohomology is None:
        shf_cohomology = sheaf_cohomology(shf)
    sections = shf_cohomology[0]
    att_secs = pychomp.DirectedAcyclicGraph()
    zero_sec = 0*sections[0]

    for c in itertools.product(*[[0, 1] for s in sections]):
        add = True
        section = sum([cv*sec for cv, sec in zip(c, sections)], zero_sec)
        mss = morse_sets_in_section(shf, morse_dict, section)
        for key in shf.grading[0]:
            selected = mss[key]
            P = pychomp.Poset(morse_dict[key])
            if set(mss[key]) != set(mss[key]).union(*[set(P.descendants(m)) 
                                                      for m in mss[key]]):
                add = False
                break
        if add:
            att_secs.add_vertex(tuple([int(s==1) for s in section]))

    for s1, s2 in itertools.combinations(att_secs.vertices(), 2):
        if all([s1[i] <= s2[i] for i in range(len(s1))]):
            att_secs.add_edge(s2, s1)
            continue
        elif all([s2[i] <= s1[i] for i in range(len(s1))]):
            att_secs.add_edge(s1, s2)

    return pychomp.Poset(att_secs)
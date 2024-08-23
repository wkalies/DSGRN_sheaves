### SheafFigure.py
### MIT LICENSE 2024 Alex Dowling

import galois
from itertools import chain
from .ParameterComplexFigure import *

class SheafFigure(ParameterComplexFigure):

    def propagate_section(self, section):
        assignment = {}
        dim = len(self.shf.grading) - 1
        
        n = 0
        for cell in self.shf.grading[0]:
            m = len(self.shf.stalk(cell))
            assignment[cell] = section[n:n+m]
            n = n+m
           
        for cell_0 in chain(*self.shf.grading):
            for cell_1 in self.shf.P.children(cell_0):
                if cell_1 in assignment:
                    continue
                assignment[cell_1] = self.shf.restriction(cell_0, cell_1,
                                                          assignment[cell_0])
        return assignment

    def build_section_dict(self, section):
        assignment = self.propagate_section(section)
        self.section_dict = {}
        
        for cell, v in assignment.items():
            morse_sets = [M for M, vi in zip(self.shf.stalk(cell), v)
                          if vi != 0]
            self.section_dict.update({cell : morse_sets})

    def plot_cell(self, cell, ax):
        stg = self.stg_dict[cell]
        if self.parameter_graph:
            label = str(cell.get_descendant_parameters(self.parameter_graph))
        else:
            label = None

        if self.prune_grad == 'all' or (self.prune_grad == 'some' 
                                        and cell.dim == self.dim):
            plot_stg(stg, morse_sets=self.section_dict[cell], 
                     plot_bdry_cells=self.plot_bdry_cells, 
                     prune_grad=True, ax=ax, label=label, 
                     figsize=self.figsize)
        else:
            plot_stg(stg, morse_sets=self.section_dict[cell],
                     plot_bdry_cells=self.plot_bdry_cells, 
                     prune_grad=False, ax=ax, label=label, 
                     figsize = self.figsize)

    def plot(self, section=None, fname=None):
        if section is None:
            rank = sum(len(self.shf.stalk(cell)) 
                       for cell in self.shf.grading[0])
            section = galois.GF2([1 for i in range(rank)])
        self.build_section_dict(section)
        super().plot(fname)
    
    def __init__(self, shf, stg_dict, parameter_graph=None, columns=3):
        self.shf = shf
        super().__init__(shf.P.children_, stg_dict, 
                         parameter_graph=parameter_graph, columns=columns)
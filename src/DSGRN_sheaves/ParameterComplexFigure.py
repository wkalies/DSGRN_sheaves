### ParameterComplexFigure.py
### MIT LICENSE 2024 Alex Dowling

import DSGRN
import DSGRN_utils
import pychomp
import math
import galois
import itertools

from DSGRN_sheaves.Sheaf import *
from DSGRN_sheaves.Continuation import *

import matplotlib
import matplotlib.pyplot as plt

def plot_stg(stg, morse_nodes=None, plot_bdry_cells=True, prune_grad=False, 
             ax=None, label=None, figsize=(7,7), fname=None, visible=True):
    """ Inputs a RookRulesCubicalComplex object. Builds matplotlib figure which 
        shows the Morse sets and state transition graph. """
    
    (scc_dag, graded_complex) = pychomp.FlowGradedComplex(
                                        stg.complex(), 
                                        stg.digraph.adjacencies)
    connection_matrix = pychomp.ConnectionMatrix(graded_complex)
    morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_complex, 
                                         connection_matrix, prune_grad)
    saveable = False
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
        saveable = fname is not None
    if label:
        ax.title.set_text(label)
    DSGRN_utils.PlotMorseSets(morse_graph, stg, graded_complex, ax=ax, 
                              morse_nodes=morse_nodes, 
                              plot_bdry_cells=plot_bdry_cells)
    if saveable:
        fig.savefig(fname=fname)
    if not visible:
        plt.close(fig)

class ParameterComplexFigure:

    def set_defaults(self):
        self.prune_grad = 'none'
        self.plot_bdry_cells = True
        self.figsize = (8,8)
        self.labels = True
        self.top_only = False
        self.visible = True

    def build_plot_list(self):
        if self.top_only:
            self.plot_list = [cell for cell 
                              in self.parameter_complex.vertices() 
                              if cell.dim == self.dim]
        else:
            self.plot_list = list(self.parameter_complex.vertices())

    def clean_axs(self):
        for ax in self.axs.reshape(-1):
            if not ax.lines:
                ax.axis('off')

    def plot_cell(self, cell, ax):
        stg = self.stg_dict[cell]
        if self.parameter_graph:
            label = str(cell.get_descendant_parameters(self.parameter_graph))
        else:
            label = None

        if self.prune_grad == 'all' or (self.prune_grad == 'some' 
                                        and cell.dim == self.dim):
            plot_stg(stg, plot_bdry_cells=self.plot_bdry_cells, 
                     prune_grad=True, ax=ax, label=label, 
                     figsize=self.figsize)
        else:
            plot_stg(stg, plot_bdry_cells=self.plot_bdry_cells, 
                     prune_grad=False, ax=ax, label=label, 
                     figsize = self.figsize)
    
    def plot(self, fname=None):
        self.rows = math.ceil(len(self.plot_list)/self.columns)
        self.fig, self.axs = plt.subplots(self.rows, self.columns, 
                                          figsize=self.figsize)
        self.axs = np.array(self.axs)
        self.ax_list = self.axs.reshape(-1)
        
        for i, cell in enumerate(self.plot_list):
            ax = self.ax_list[i]
            self.plot_cell(cell, ax)
        self.clean_axs()
        
        if fname:
            self.fig.savefig(fname=self.fname)
        if not self.visible:
            plt.close(self.fig)

    def __init__(self, parameter_complex, stg_dict, parameter_graph=None, 
                 columns=3):
        self.parameter_complex = parameter_complex
        self.stg_dict = stg_dict
        self.parameter_graph = parameter_graph
        self.columns = columns
        self.set_defaults()
        self.build_plot_list()

class ParameterPathFigure(ParameterComplexFigure):

    def build_plot_list(self):
        self.plot_list = []
        pc = self.parameter_complex
        pct = self.parameter_complex.transpose()
        
        current_cell = next(cell for cell in pc.vertices() 
                            if len(pc.adjacencies(cell)) == 1)
        self.plot_list.append(current_cell)
        adj_cells = list(pc.adjacencies(current_cell))

        while adj_cells:
            top_cell = adj_cells[0]
            self.plot_list.append(top_cell)
            next_cell = next(cell for cell in pct.adjacencies(top_cell)
                             if cell != current_cell)
            self.plot_list.append(next_cell)

            current_cell = next_cell
            adj_cells = list(pc.adjacencies(next_cell))
            adj_cells.remove(top_cell)

    def __init__(self, parameter_complex, stg_dict, parameter_graph=None):
        super().__init__(parameter_complex, stg_dict, 
                         parameter_graph=parameter_graph, 
                         columns=len(parameter_complex.vertices()))

class ParameterLoopFigure(ParameterComplexFigure):

    def build_plot_list(self):
        self.plot_list = []
        pc = self.parameter_complex
        pct = self.parameter_complex.transpose()

        if not self.start:
            self.start = next(cell for cell in pc.vertices() 
                              if cell.dim == 1)
        current_cell = self.start
        self.plot_list.append(current_cell)
        adj_cells = list(pct.adjacencies(current_cell))

        while True:
            bottom_cell = adj_cells[0]
            self.plot_list.append(bottom_cell)
            next_cell = next(cell for cell in pc.adjacencies(bottom_cell)
                             if cell != current_cell)
            self.plot_list.append(next_cell)

            current_cell = next_cell
            adj_cells = list(pct.adjacencies(next_cell))
            adj_cells.remove(bottom_cell)
            if current_cell == self.start: break

        self.plot_list.pop()

    def plot(self, fname=None):
        self.fig, self.axs = plt.subplots(self.rows, self.columns, 
                                          figsize=self.figsize)
        coord_list =  [(0, i) for i in range(self.columns)]
        coord_list += [(i, self.columns-1) for i in range(1, self.rows)]
        coord_list += [(self.rows-1, self.columns-i) 
                       for i in range(2, self.columns+1)]
        coord_list += [(self.rows-i, 0) for i in range(2, self.rows)]
       
        for n, cell in enumerate(self.plot_list):
            i, j = coord_list[n]
            ax = self.axs[i, j]
            self.plot_cell(cell, ax)
        self.clean_axs()
        
        if fname:
            self.fig.savefig(fname=self.fname)
        if not self.visible:
            plt.close(self.fig)


    def __init__(self, parameter_complex, stg_dict, parameter_graph=None, 
                 start=None):
        self.start = start
        super().__init__(parameter_complex, stg_dict, 
                         parameter_graph=parameter_graph)
        num_cells = len(self.parameter_complex.vertices())
        self.columns = (num_cells - 2) // 2
        self.rows = 3

class SheafFigure(ParameterComplexFigure):

    def propagate_section(self, section):
        assignment = {}
        dim = len(self.shf.grading) - 1
        
        n = 0
        for cell in self.shf.grading[0]:
            m = len(self.shf.stalk(cell))
            assignment[cell] = section[n:n+m]
            n = n+m
           
        for cell_0 in itertools.chain(*self.shf.grading):
            for cell_1 in self.shf.P.children(cell_0):
                if cell_1 in assignment:
                    continue
                assignment[cell_1] = self.shf.restriction(cell_0, cell_1,
                                                          assignment[cell_0])
        return assignment

    def build_morse_nodes(self, section):
        assignment = self.propagate_section(section)
        self.morse_nodes = {}
        
        for cell, v in assignment.items():
            self.morse_nodes.update({cell : [i for i, vi in enumerate(v) 
                                             if vi != 0]}) 

    def plot_cell(self, cell, ax):
        stg = self.stg_dict[cell]
        if self.parameter_graph:
            label = str(cell.get_descendant_parameters(self.parameter_graph))
        else:
            label = None

        if self.prune_grad == 'all' or (self.prune_grad == 'some' 
                                        and cell.dim == self.dim):
            plot_stg(stg, morse_nodes=self.morse_nodes[cell], 
                     plot_bdry_cells=self.plot_bdry_cells, 
                     prune_grad=True, ax=ax, label=label, 
                     figsize=self.figsize)
        else:
            plot_stg(stg, morse_nodes=self.morse_nodes[cell],
                     plot_bdry_cells=self.plot_bdry_cells, 
                     prune_grad=False, ax=ax, label=label, 
                     figsize = self.figsize)

    def plot(self, section=None, fname=None):
        if section is None:
            rank = sum(len(self.shf.stalk(cell)) 
                       for cell in self.shf.grading[0])
            section = galois.GF2([1 for i in range(rank)])
        self.build_morse_nodes(section)
        super().plot(fname)
    
    def __init__(self, shf, stg_dict, parameter_graph=None, columns=3):
        self.shf = shf
        super().__init__(shf.P.children_, stg_dict, 
                         parameter_graph=parameter_graph, columns=columns)

class PathSheafFigure(SheafFigure, ParameterPathFigure):
    def __init__(self, shf, stg_dict, parameter_graph=None):
        self.shf = shf
        ParameterPathFigure.__init__(self, shf.P.children_, stg_dict, 
                                     parameter_graph=parameter_graph)

class LoopSheafFigure(SheafFigure, ParameterLoopFigure):
    def __init__(self, shf, stg_dict, parameter_graph=None, start=None):
        self.shf = shf
        ParameterLoopFigure.__init__(self, shf.P.children_, stg_dict, 
                                     parameter_graph=parameter_graph, 
                                     start=start)
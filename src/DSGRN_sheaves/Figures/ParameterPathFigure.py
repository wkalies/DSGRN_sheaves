### ParameterPathFigure.py
### MIT LICENSE 2024 Alex Dowling

from .ParameterComplexFigure import *

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
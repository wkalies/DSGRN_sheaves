### ParameterLoopFigure.py
### MIT LICENSE 2024 Alex Dowling

from .ParameterComplexFigure import *

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
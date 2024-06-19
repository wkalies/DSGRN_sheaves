### CechCell.py
### MIT LICENSE 2024 Alex Dowling

import DSGRN
import itertools

class CechCell:

    def __init__(self, inequality_sets, dim, labels=None):
        self.inequality_sets = inequality_sets
        self.dim = dim
        self.labels = labels

    def __iter__(self):
        return iter(self.inequality_sets)
        
    def __getitem__(self, item):
        return self.inequality_sets[item]

    def __str__(self):
        ineq_strings = [f'cell dimension: {self.dim}\n']
        if self.labels:
            labels = self.labels
        else:
            labels = ['' for ineqs in self]
        for ineqs, label in zip(self, labels):
            if self.labels: ineq_strings.extend([label, " : "])
            ineq_strings.append("{")
            for ineq in ineqs:
                ineq_strings.extend([str(ineq).replace("'",""), ", "])
            ineq_strings[-1] = "}\n"
        ineq_strings[-1] = "}"
        return ''.join(ineq_strings)

    def __eq__(self, other):
        if not isinstance(other, CechCell):
            return False
        return (self.inequality_sets == other.inequality_sets 
                and self.dim == other.dim)

    def __hash__(self):
        return hash(self.inequality_sets + (self.dim,))

    def permute(self):
        if all([len(ineqs) < 3 for ineqs in self]): return self
        new_inequality_sets = []
        for ineqs in self:
            new_inequality_set = set(ineqs)
            for ineq_1, ineq_2 in itertools.combinations(ineqs, 2):
                indexes = [ineq_1.index(a) for a in ineq_2]
                permutation = lambda ineq : tuple(ineq[i] for i in indexes)
                new_inequality_set.update(map(permutation, ineqs))     
            new_inequality_sets.append(frozenset(new_inequality_set))
        return CechCell(tuple(new_inequality_sets), self.dim, self.labels)

    def join(self, others, dim=None):
        if dim is None:
            dim = self.dim - len(others)
        if not others: return self
        other_ineq_sets = tuple(zip(*[other.inequality_sets 
                                     for other in others]))
        all_ineq_sets = tuple(zip(self.inequality_sets, other_ineq_sets))
        new_ineq_sets = tuple(ineqs.union(*other_ineqs) 
                                    for ineqs, other_ineqs in all_ineq_sets)
        return CechCell(new_ineq_sets, dim, self.labels)

    def get_child_cells(self):
        child_cells = []
        for d, ineqs in enumerate(self):
            if len(ineqs) == 1:
                continue
            for ineq in ineqs:
                new_ineq_sets = (self.inequality_sets[:d] 
                                  + (ineqs.difference({ineq}),) 
                                  + self.inequality_sets[d+1:])
                child = CechCell(new_ineq_sets, self.dim + 1, self.labels)
                child_cells.append(child)
        return child_cells

    def get_descendant_parameters(self, parameter_graph):
        index = lambda p : DSGRN.index_from_partial_orders(parameter_graph, p)
        partial_orders = [[str(ineq).replace("'","") for ineq in ineqs] 
                          for ineqs in itertools.product(*self.inequality_sets)]
        return list(map(index, partial_orders))

def top_cech_cell(parameter_graph, index, dim=None):
    lines = parameter_graph.parameter(index).partialorders().split("\n")
    if dim is None:
        dim = parameter_graph.dimension()
    colons = [line.index(":") for line in lines]
    
    inequality_sets = tuple(frozenset({tuple(line[i+3:-1].split(", "))})
                            for line, i in zip(lines, colons))
    labels = tuple(line[:i-1] for line, i in zip(lines, colons))
    
    return CechCell(inequality_sets, dim, labels)
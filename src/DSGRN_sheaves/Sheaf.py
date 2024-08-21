### Sheaf.py
### MIT LICENSE 2024 Alex Dowling

import galois
import numpy as np
import pychomp.Poset

from .Cohomology import *

class Sheaf:
    """ Class for storing the stalks and restriction maps of a sheaf.

        Initialized with:
        - P:          a poset
        - grading:    a list of lists of elements of p such that each element
                      of P belongs to exactly one list in grading
        - stalk_dict: a dictionary of lists whose keys are the elements of P
        - restr_dict: a dictionary of matrices whose keys are incident pairs 
                      of elements of P.
    """
    # Currently only supports sheaves of vector spaces over Z2
    
    def __init__(self, P, grading, stalk_dict, restr_dict):
        self.P = P
        self.grading = grading

        self.GF = galois.GF(2)
        self.stalk_dict = stalk_dict
        self.restr_dict = {}
        for (p1, p2) in restr_dict:
            self.restr_dict.update({(p1, p2) : 
                                    self.GF(restr_dict[(p1, p2)])})          
          
    def stalk(self, p):
        """ Inputs an element p of the sheaf's underlying poset. 
        
            Outputs the stored basis for the stalk at p. 
        """

        return self.stalk_dict[p]
    
    def restriction(self, p1, p2, f = None):
        """ Inputs a pair of incident cells p1 and p2, and an optional element 
            to evaluate. 
            
            Outputs the matrix for the stalk restriction map, or
            the stalk restriction map evaluated at the inputted element. 
        """
        
        if p1 not in self.P.parents(p2):
            raise Exception("The input pair of cells are not incident.")

        R = self.restr_dict[(p1, p2)]
        if f is None:
            return R
        elif self.GF(f).shape[0] != R.shape[1]:
            raise Exception("Invalid input vector.")
        else:
            return np.matmul(R, self.GF(f))

def build_coboundary_map(shf, k):
    """ Inputs a cellular sheaf and an index k

        Outputs the kth coboundary map for the cochain complex of the sheaf.
    """

    # Initialize D as a matrix of zeroes
    cochainRank_1 = sum([len(shf.stalk(cell)) 
                         for cell in shf.grading[k]])
    cochainRank_2 = sum([len(shf.stalk(cell)) 
                         for cell in shf.grading[k+1]])     
    D = shf.GF([[0 for i in range(cochainRank_1)] 
                for j in range(cochainRank_2)])   
    
    col = 0
    for cell_1 in shf.grading[k]:
        row = 0
        n = len(shf.stalk(cell_1))
        for cell_2 in shf.grading[k+1]:
            m = len(shf.stalk(cell_2))
            # If (cell_1, cell_2) is an incident pair, load the 
            # corresponding block of D with the sheaf restriction.
            # Otherwise, leave it to be zero.
            if cell_1 in shf.P.parents(cell_2):
                D[row:row+m, col:col+n] = shf.restriction(cell_1, cell_2)
            row = row + m
        col = col + n

    return D

def build_coboundary_maps(shf):
    """ Inputs a cellular sheaf.

        Outputs a dictionary delta which stores the kth coboundary map.
    """

    dim = len(shf.grading) - 1
    # Build edge case coboundaries
    cochainRank_1 = sum([len(shf.stalk(cell)) for cell in shf.grading[0]])
    cochainRank_2 = sum([len(shf.stalk(cell)) for cell in shf.grading[dim]])
    delta = {-1 : shf.GF([[0] for i in range(cochainRank_1)])} 
    delta.update({dim : shf.GF([[0 for i in range(cochainRank_2)]])}) 
    
    # Build the rest
    for k in range(dim):
        delta.update({k : build_coboundary_map(shf, k)})
        
    return delta

def sheaf_cohomology(shf):
    """ Inputs a cellular sheaf.
    
        Outputs a list of lists of generators for each sheaf cohomology group. 
    """

    dim = len(shf.grading) - 1
    delta = build_coboundary_maps(shf)
    return [cohomology(delta[j], delta[j-1]) for j in range(dim+1)]

def betti_numbers(shf):
    return [len(generators) for generators in sheaf_cohomology(shf)]

def clean_stalks(shf):
    """ Inputs a cellular sheaf. 
        
        Outputs a sheaf with the same stalks at top dimensional cells. The 
        stalks at lower dimensional cells are stripped of elements which map 
        to zero under every restriction map. """
    
    dim = len(shf.grading) - 1
    P = pychomp.Poset(shf.P.children_.clone())
    grading = shf.grading.copy()
    stalk_dict = {cell : shf.stalk(cell) for cell in grading[dim]}
    restr_dict = {}
    # Keep all the indices at the top cells
    keep = {cell : list(range(len(shf.stalk(cell)))) 
            for cell in grading[dim]}
    
    # Look at stalks descending in dimension
    # Skip the top cells
    for codim in range(1,dim+1):
        for cell_1 in grading[dim-codim]:
            # Keep indices of elements in the stalk at cell_1 
            # which don't map to zero for at least one restriction map
            keep[cell_1] = []
            n = len(shf.stalk(cell_1))
            # Iterate through incident cells
            for cell_2 in P.children(cell_1):
                # Get restriction map
                R = shf.restr_dict[(cell_1, cell_2)]
                for j in range(n):
                    # Don't check an index if it's already being kept
                    if j in keep[cell_1]:
                        continue
                    # Keep the index if jth column of restriction map 
                    # has a nonzero entry
                    if any([R[i][j] != 0 for i in keep[cell_2]]):
                        keep[cell_1].append(j)
            keep[cell_1].sort()
            
            # New stalk consists of elements of old stalk whose index is kept
            stalk_dict[cell_1] = [shf.stalk(cell_1)[i] for i in keep[cell_1]]
            # Shrink domain and range for new restriction maps
            for cell_2 in P.children(cell_1):
                R = shf.restr_dict[(cell_1, cell_2)]
                restr_dict[(cell_1, cell_2)] = [[R[i][j] for j in keep[cell_1]]
                                                 for i in keep[cell_2]]
            
    return Sheaf(P, grading, stalk_dict, restr_dict)

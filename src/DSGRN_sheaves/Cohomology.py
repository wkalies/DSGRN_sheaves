### Cohomology.py
### MIT LICENSE 2024 Alex Dowling

### Based on code by Jeremy Kun: 
### http://jeremykun.com/2013/04/10/computing-homology/

import numpy
import galois

def cohomology(D, E):
    """ Inputs two matrices D, E with D*E=0. Outputs a linearly independent set
        of vectors in ker D which form a basis for ker D / img E under the 
        quotient map"""
    
    A, B, Q = simultaneous_reduce(D, E)
    B_pivot_rows = [i for i in range(B.shape[0]) 
                    if not numpy.all(B[i, :] == 0*B[i, :])]
    B_reduced = B[B_pivot_rows, :].copy().transpose().row_reduce()
    B_reduced_rows = []
    
    j = 0
    for i in range(B_reduced.shape[0]):
        while j < B_reduced.shape[1] and B_reduced[i,j]!=1:
            j += 1
        if j == B_reduced.shape[1]:
            break
        else:
            B_reduced_rows.append(j)
            
    B_span_rows = [B_pivot_rows[i] for i in B_reduced_rows]
    A_zero_cols = [j for j in range(A.shape[1]) 
                   if numpy.all(A[:, j] == 0*A[:, j])]
    generators = [Q[:, i] for i in A_zero_cols if i not in B_span_rows]
    
    return generators

def row_swap(A, i, j):
    temp = A[i, :].copy()
    A[i, :] = A[j, :]
    A[j, :] = temp

def col_swap(A, i, j):
    temp = A[:, i].copy()
    A[:, i] = A[:, j]
    A[:, j] = temp

def row_scale(A, i, c):
    A[i, :] = A[i, :]*c

def col_scale(A, i, c):
    A[:, i] = A[:, i]*c

def row_combine(A, add_to, row_scale, scale_amt):
    A[add_to, :] = A[add_to, :] + scale_amt*A[row_scale, :]
    
def col_combine(A, add_to, col_scale, scale_amt):
    A[:, add_to] = A[:, add_to] + scale_amt*A[:, col_scale]
    
def simultaneous_reduce(C, D):
    A = C.copy()
    B = D.copy()
    
    if A.shape[1] != B.shape[0]:
        raise Exception("Matrices have the wrong shape.")
        
    Q = A.copy()
    Q.resize(A.shape[1],A.shape[1])
    Q = Q*0
    for i in range(A.shape[1]):
        Q[i,i] = 1

    num_rows, num_cols = A.shape

    i = 0
    j = 0
    while True:
        if i >= num_rows or j >= num_cols:
            break

        if A[i][j] == 0:
            nonzero_col = j
            while nonzero_col < num_cols and A[i,nonzero_col] == 0:
                nonzero_col += 1

            if nonzero_col == num_cols:
                i += 1
                continue
            col_swap(A, j, nonzero_col)
            row_swap(B, j, nonzero_col)
            col_swap(Q, j, nonzero_col)

        pivot = A[i,j]
        col_scale(A, j, pivot**-1)
        row_scale(B, j, pivot**-1)
        col_scale(Q, j, pivot**-1)

        for other_col in range(num_cols):
            if other_col == j:
                continue
            if A[i, other_col] != 0:
                scale_amt = -A[i, other_col]
                col_combine(A, other_col, j, scale_amt)
                row_combine(B, j, other_col, -scale_amt)
                col_combine(Q, other_col, j, scale_amt)

        i += 1
        j += 1
        
    return A,B,Q

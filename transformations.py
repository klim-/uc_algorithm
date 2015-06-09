# -*- coding: utf-8 -*-

import sympy as sp
import numpy as np
import symb_tools as st
import non_commutative_tools as nct

from print_candy import *

from IPython import embed as IPS

t = nct.t
s = nct.s

def calculate_Q_matrix(myStack):
    print_line()
    print("Algorithmus am Ende\n")

    Q_relevant_matrices, Q_tilde_relevant_matrices = myStack.get_Q_relevant_matrices()

    def multiply_matrices_in_list(matrix_list):
        result = matrix_list[0]
        for i in xrange(1, len(matrix_list)):
            result = matrix_list[i]*result
        return result

    Q = multiply_matrices_in_list( Q_relevant_matrices )

    #IPS()

    if not Q_tilde_relevant_matrices==[]:
        Q_tilde = multiply_matrices_in_list( Q_tilde_relevant_matrices )
        Q_matrix = st.concat_rows(Q, Q_tilde)

    m, n = Q.shape
    assert n==len(myStack.vec_x), "Dimensions of Q-Matrix do not fit."

    print "Q-matrix = "
    print_nicely(Q_matrix)

    return Q_matrix

def calculate_G_matrix(myStack):
    pass
    #a = 1.3
    
    #IPS()

    #for i in xrange(myStack.iteration_steps()+1):
        #this_iter = myStack.get_iteration(i)
        #Gi = this_iter.P1_rpinv - this_iter.P1_roc*( this_iter.B_lpinv*s + this_iter.B_lpinv*this_iter.A )
        
